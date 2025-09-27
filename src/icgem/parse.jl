## Description #############################################################################
#
# Functions to parse ICGEM files.
#
## References ##############################################################################
#
# [1] Barthelmes, F., Förste, C (2011). The ICGEM-format. GFZ Postdam, Department 1
#     "Geodesy and Remote Sensing".
#
############################################################################################

############################################################################################
#                                        Functions                                         #
############################################################################################

"""
    parse_icgem(filename::AbstractString, T::DataType = Float64) -> IcgemFile

Parse the ICGEM file `filename` using the data type `T`.

!!! note

    `T` is converted to float to obtain the output type.
"""
function parse_icgem(filename::AbstractString, T::DataType = Float64)
    Tf = float(T)

    # Open the file and find the header.
    file = open(filename, "r")

    # == Header ============================================================================

    header_start_line   = 1
    header_end_line     = 0
    current_line        = 0
    begin_of_head_found = false
    end_of_head_found   = false

    # We need to first check the position of the header due to the `begin_of_head` keyword
    # that can define where the header starts.
    while !eof(file)
        current_line += 1
        tokens = split(readline(file))

        length(tokens) < 1 && continue

        # Search for the `begin_of_head`, which is optional and makes all previous lines to
        # be ignored.
        if tokens[1] == "begin_of_head"
            # We should not have two `begin_of_head`.
            if begin_of_head_found
                error("[Invalid ICGEM file] Two `begin_of_head` keywords were found!")
            end

            header_start_line = current_line
            begin_of_head_found = true
            continue

        elseif tokens[1] == "end_of_head"
            end_of_head_found = true
            header_end_line = current_line
            break
        end
    end

    # `end_of_head` keyword is mandatory.
    if !end_of_head_found
        error("[Invalid ICGEM file] The mandatory keyword `end_of_head` was not found!")
    end

    # Rewind file to read again.
    seek(file, 0)
    current_line = 0

    # == Get Keywords ======================================================================

    keywords = Dict{Symbol, String}()

    while current_line < header_end_line - 1
        line = readline(file)
        current_line += 1

        # Skip all lines until the beginning of the header.
        current_line < header_start_line && continue
        tokens = split(line)

        # We must have two keywords, otherwise we do not have a keyword. Here, we will just
        # skip the line because old versions of ICGEM files does not define well where the
        # header starts and comments are allowed.
        length(tokens) != 2 && continue

        keywords[Symbol(tokens[1])] = tokens[2]
    end

    # Read one mode line to take into account the "end_of_head" line.
    readline(file)
    current_line += 1

    # == Parse Mandatory Header Fields =====================================================

    mandatory_fields = (
        :product_type,
        :modelname,
        :earth_gravity_constant,
        :radius,
        :max_degree,
        :errors,
    )

    has_mandatory_fields = map(f -> haskey(keywords, f), mandatory_fields)

    if !prod(has_mandatory_fields)
        missing_fields = mandatory_fields[findall(!, has_mandatory_fields)]
        error("[Invalid ICGEM file] The following mandatory fields are missing: $missing_fields.")
    end

    product_type = Symbol(keywords[:product_type])
    model_name   = keywords[:modelname]
    max_degree   = parse(Int, keywords[:max_degree])
    errors       = Symbol(keywords[:errors])

    gravity_constant = _parse_icgem_float(Tf, keywords[:earth_gravity_constant])
    radius           = _parse_icgem_float(Tf, keywords[:radius])

    isnothing(gravity_constant) && error("[Invalid ICGEM file] Could not parse the gravity constant to $Tf.")
    isnothing(radius) && error("[Invalid ICGEM file] Could not parse the radius to $Tf.")

    # Check if some keywords are valid.
    gravity_constant <= 0 && error("[Invalid ICGEM file] The gravity constant must be positive.")

    radius <= 0 && error("[Invalid ICGEM file] The radius must be positive.")

    max_degree < 0 && error("[Invalid ICGEM file] The maximum degree must not be negative.")

    errors ∉ (:no, :calibrated, :calibrated_and_formal, :formal) &&
        error("[Invalid ICGEM file] An invalid value was found for the keyword `errors`.")

    # == Parse Optional Keywords ===========================================================

    tide_system = haskey(keywords, :tide_system) ? Symbol(keywords[:tide_system]) : :unknown
    norm        = haskey(keywords, :norm) ? Symbol(keywords[:norm]) : :fully_normalized

    # == Data ==============================================================================

    # Since we now have the maximum degree, we can pre-allocate and initialize the data
    # matrix.
    data_static = LowerTriangularStorage{IcgemGfcCoefficient{Tf}}(max_degree + 1,
        max_degree + 1)

    use_dynamic = false
    data_dynamic = LowerTriangularStorage{IcgemGfctCoefficient{Tf}}(max_degree + 1,
        max_degree + 1
    )


    # State of the parsing algorithm.
    state = :new

    # Auxiliary variables to build the coefficients.
    deg  = 0
    ord  = 0
    clm  = Tf(0)
    slm  = Tf(0)
    time = Dates.value(now() - _DT_J2000) / 1000

    has_trend = false
    trend_clm = Tf(0)
    trend_slm = Tf(0)
    asin_coefficients = NTuple{3, Tf}[]
    acos_coefficients = NTuple{3, Tf}[]

    line          = nothing
    read_new_line = true
    tokens        = nothing

    # Read the entire file and build the coefficients.
    while !eof(file)
        # Check if we need to read a new line from the file.
        if read_new_line
            # Read and tokenize each line.
            line = readline(file)
            current_line += 1
            tokens = split(line)
        end

        # Process the line according to the state.
        if state === :new
            if length(tokens) < 3
                @warn "[Line $current_line] Invalid data line."
                continue
            end

            # == `gfc` Data Line ===========================================================

            if tokens[1] == "gfc"
                ret = _parse_gfc_data_line(Tf, tokens, current_line)
                isnothing(ret) && continue

                deg, ord, clm, slm = ret
                if use_dynamic
                    data_dynamic[deg+1, ord+1] = IcgemGfctCoefficient(IcgemGfcCoefficient(clm, slm))
                else
                    data_static[deg+1, ord+1] = IcgemGfcCoefficient(clm, slm)
                end

                read_new_line = true

            # == `gfct` Data Line ==========================================================

            elseif tokens[1] == "gfct"
                ret = _parse_gfct_data_line(Tf, tokens, current_line)
                isnothing(ret) && continue
                deg, ord, clm, slm, time = ret

                # Now, we need to change the state to wait for the next terms.
                state = :gfct
                read_new_line = true

                has_trend = false
                trend_clm = T(0)
                trend_slm = T(0)
                empty!(asin_coefficients)
                empty!(acos_coefficients)
            end

        elseif state === :gfct

            # == `trnd` Data Line of a `gfct` Section ======================================

            if tokens[1] == "trnd"
                ret = _parse_trnd_data_line(Tf, tokens, current_line)
                if isnothing(ret)
                    read_new_line = true
                    continue
                end

                adeg, aord, trend_clm, trend_slm = ret

                ((adeg != deg) || (aord != ord)) &&
                    error("[Invalid ICGEM file] The degree or order of a `trnd` line is different from the corresponding `gfct` line.")

                has_trend = true
                read_new_line = true

            # == `asin` Data Line of a `gfct` Section ======================================

            elseif tokens[1] == "asin"
                ret = _parse_asin_acos_data_line(Tf, tokens, current_line)
                if isnothing(ret)
                    read_new_line = true
                    continue
                end

                adeg, aord, asin_amplitude_clm, asin_amplitude_slm, asin_period = ret

                ((adeg != deg) || (aord != ord)) &&
                    error("[Invalid ICGEM file] The degree or order of a `asin` line is different from the corresponding `gfct` line.")

                push!(asin_coefficients, (asin_amplitude_clm, asin_amplitude_slm, asin_period))
                read_new_line = true

            # == `acos` Data Line of a `gfct` Section ======================================

            elseif tokens[1] == "acos"
                ret = _parse_asin_acos_data_line(Tf, tokens, current_line)
                if isnothing(ret)
                    read_new_line = true
                    continue
                end
                adeg, aord, acos_amplitude_clm, acos_amplitude_slm, acos_period = ret

                ((adeg != deg) || (aord != ord)) &&
                    error("[Invalid ICGEM file] The degree or order of a `acos` line is different from the corresponding `gfct` line.")

                push!(acos_coefficients, (acos_amplitude_clm, acos_amplitude_slm, acos_period))
                read_new_line = true

            else
                # If we reach this part, the `gfct` section is over. Thus, we should create
                # the element related to `gfct` and proceed with the new information.
                if !use_dynamic
                    use_dynamic = true
                    for i in 1:(max_degree+1), j in 1:(max_degree+1)
                        data_dynamic[i, j] = IcgemGfctCoefficient(data_static[i, j])
                    end
                end

                data_dynamic[deg+1, ord+1] = IcgemGfctCoefficient(
                    clm,
                    slm,
                    time,
                    has_trend,
                    trend_clm,
                    trend_slm,
                    copy(asin_coefficients),
                    copy(acos_coefficients),
                )

                state = :new
                read_new_line = false
            end
        end
    end

    # Create the ICGEM object.
    icgem_file = IcgemFile(
        product_type,
        model_name,
        gravity_constant,
        radius,
        max_degree,
        errors,
        tide_system,
        Val(norm),
        use_dynamic ? data_dynamic : data_static
    )
    return icgem_file
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

#   _parse_icgem_float(T, input) -> Uniont{Nothing, T}
#
# Parse the `input` to float type `T` substituting all `D`s and `d`s  to `e`, so that we can
# convert numbers in FORTRAN format. If we cannot parse `input` to `T`, it returns
# `nothing`.
function _parse_icgem_float(T::DataType, input::AbstractString)
    data_str = replace(input, r"[D, d]" => "e")
    return tryparse(T, data_str)
end

# == Functions to Parse Data Lines =========================================================

#   _parse_degree_and_order(tokens, current_line) -> Int, Int
#
# Parse the degree in `tokens[2]` and order in `tokens[3]`. The `current_line` number is
# used for debugging purposes.
#
# # Returns
#
# - `Int`: Degree.
# - `Int`: Order.
function _parse_degree_and_order(tokens, current_line)
    deg = tryparse(Int, tokens[2])

    if isnothing(deg)
        @warn "[Line $current_line] Invalid degree: $(tokens[2])."
        return nothing
    end

    ord = tryparse(Int, tokens[3])

    if isnothing(ord)
        @warn "[Line $current_line] Invalid order: $(tokens[3])."
        return nothing
    end

    return deg, ord
end

#   _parse_gfc_data_line(Tf, tokens, current_line) -> Int, Int, Tf, Tf
#
# Parse the `gfc` data line in `tokens` using the data type `Tf` for the floating point
# fields. The `current_line` number is used for debugging purposes.
#
# # Returns
#
# - `Int`: Degree.
# - `Int`: Order.
# - `T`: `Clm` coefficient.
# - `T`: `Slm` coefficient.
function _parse_gfc_data_line(Tf, tokens, current_line)
    if length(tokens) < 5
        @warn "[Line $current_line] Invalid `gfc` data line."
        return nothing
    end

    ret = _parse_degree_and_order(tokens, current_line)
    isnothing(ret) && return nothing
    deg, ord = ret

    clm = _parse_icgem_float(Tf, tokens[4])

    if isnothing(clm)
        @warn "[Line $current_line] Could not parse `Clm` to $Tf: $(tokens[4])."
        return nothing
    end

    slm = _parse_icgem_float(Tf, tokens[5])

    if isnothing(slm)
        @warn "[Line $current_line] Could not parse `Slm` to $Tf: $(tokens[5])."
        return nothing
    end

    return deg, ord, clm, slm
end

#   _parse_gfc_data_line(Tf, tokens, current_line) -> Int, Int, Tf, Tf, Number
#
# Parse the `gfct` data line in `tokens` using the data type `Tf` for the floating point
# fields. The `current_line` number is used for debugging purposes.
#
# # Returns
#
# - `Int`: Degree.
# - `Int`: Order.
# - `T`: `Clm` coefficient.
# - `T`: `Slm` coefficient.
# - `Number`: Epoch (`t₀`) of the coefficients, expressed as the number of elapsed seconds
#   since J2000.0.
function _parse_gfct_data_line(Tf, tokens, current_line)
    if length(tokens) < 6
        @warn "[Line $current_line] Invalid `gfct` data line."
        return nothing
    end

    # The first 4 tokens after the key are the same as in `gfc` line.
    ret = _parse_gfc_data_line(Tf, tokens, current_line)
    isnothing(ret) && return nothing
    deg, ord, clm, slm = ret

    # Parse the time.
    time = Dates.value(DateTime(tokens[end], dateformat"yyyymmdd") - _DT_J2000) / 1000

    return deg, ord, clm, slm, time
end

#   _parse_trnd_data_line(Tf, tokens, current_line) -> Int, Int, Tf, Tf
#
# Parse the `trnd` data line in `tokens` using the data type `Tf` for the floating point
# fields. The `current_line` number is used for debugging purposes.
#
# # Returns
#
# - `Int`: Degree.
# - `Int`: Order.
# - `T`: `Clm` trend coefficient.
# - `T`: `Slm` trend coefficient.
function _parse_trnd_data_line(Tf, tokens, current_line)
    if length(tokens) < 5
        @warn "[Line $current_line] Invalid `trnd` data line."
        return nothing
    end

    # Parse the degree and order.
    ret = _parse_degree_and_order(tokens, current_line)
    isnothing(ret) && return nothing
    deg, ord = ret

    # Parse the other coefficients.
    trend_clm = _parse_icgem_float(Tf, tokens[4])

    if isnothing(trend_clm)
        @warn "[Line $current_line] Could not parse `trend_C` to $Tf: $(tokens[4])."
        return nothing
    end

    trend_slm = _parse_icgem_float(Tf, tokens[5])

    if isnothing(trend_slm)
        @warn "[Line $current_line] Could not parse `trend_S` to $Tf: $(tokens[5])."
        return nothing
    end

    return deg, ord, trend_clm, trend_slm
end

#   _parse_asin_acos_data_line(Tf, tokens, current_line) -> Int, Int, Tf, Tf, Tf
#
# Parse the `asin` or `acos` data line in `tokens` using the data type `Tf` for the floating
# point fields. The `current_line` number is used for debugging purposes.
#
# # Returns
#
# - `Int`: Degree.
# - `Int`: Order.
# - `T`: `Clm` amplitude.
# - `T`: `Slm` amplitude.
# - `T`: Period.
function _parse_asin_acos_data_line(Tf, tokens, current_line)
    if length(tokens) < 6
        @warn "[Line $current_line] Invalid `asin` or `acos` data line."
        return nothing
    end

    # Parse the degree and order.
    ret = _parse_degree_and_order(tokens, current_line)
    isnothing(ret) && return nothing
    deg, ord = ret

    # Parse the other coefficients.
    amplitude_clm = _parse_icgem_float(Tf, tokens[4])

    if isnothing(amplitude_clm)
        @warn "[Line $current_line] Could not parse `Clm` amplitude to $Tf: $(tokens[4])."
        return nothing
    end

    amplitude_slm = _parse_icgem_float(Tf, tokens[5])

    if isnothing(amplitude_slm)
        @warn "[Line $current_line] Could not parse `Slm` amplitude to $Tf: $(tokens[4])."
        return nothing
    end

    period = _parse_icgem_float(Tf, tokens[end])

    if isnothing(period)
        @warn "[Line $current_line] Could not parse period to $Tf: $(tokens[4])."
        return nothing
    end

    return deg, ord, amplitude_clm, amplitude_slm, period
end
