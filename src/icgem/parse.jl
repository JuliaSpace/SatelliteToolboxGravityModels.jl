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
    parse_icgem(filename::AbstractString, T::Type = Float64) -> IcgemFile
    parse_icgem(io::IO, T::Type = Float64) -> IcgemFile

Parse the ICGEM file `filename`, or the ICGEM data read from the stream `io`, using the
data type `T` and return an [`IcgemFile`](@ref) object with the parsed data. The file is
closed after parsing. The stream `io` must be seekable.

This function supports ICGEM gravity model files for Earth and other celestial bodies
(Moon, planets, etc.). The parser automatically detects whether the file uses
`earth_gravity_constant` (for Earth models) or `gravity_constant` (for non-Earth models).

The function throws an [`IcgemParseError`](@ref) if the file does not conform to the ICGEM
format, and logs a warning for each invalid data line, which is skipped.

!!! note

    `T` is converted to float to obtain the output type.

# References

- **[1]** Barthelmes, F., Förste, C (2011). *The ICGEM-format*. GFZ Potsdam, Department 1
    "Geodesy and Remote Sensing".
"""
function parse_icgem(filename::AbstractString, ::Type{T} = Float64) where {T}
    return open(filename, "r") do file
        return parse_icgem(file, T)
    end
end

function parse_icgem(file::IO, ::Type{T} = Float64) where {T}
    Tf = float(T)

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
                throw(
                    IcgemParseError(
                        "Two `begin_of_head` keywords were found.", current_line
                    ),
                )
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
        throw(IcgemParseError("The mandatory keyword `end_of_head` was not found."))
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
        # skip the line because old versions of ICGEM files do not define well where the
        # header starts and comments are allowed.
        length(tokens) != 2 && continue

        keywords[Symbol(tokens[1])] = tokens[2]
    end

    # Read one more line to take into account the "end_of_head" line.
    readline(file)
    current_line += 1

    # == Parse Mandatory Header Fields =====================================================

    mandatory_fields = (:product_type, :modelname, :radius, :max_degree, :errors)

    has_mandatory_fields = map(f -> haskey(keywords, f), mandatory_fields)

    if !prod(has_mandatory_fields)
        missing_fields = mandatory_fields[findall(!, has_mandatory_fields)]
        throw(
            IcgemParseError("The following mandatory fields are missing: $missing_fields.")
        )
    end

    # Check for gravity constant - either earth_gravity_constant (Earth) or gravity_constant (other bodies)
    has_earth_gravity_constant = haskey(keywords, :earth_gravity_constant)
    has_gravity_constant = haskey(keywords, :gravity_constant)

    if !has_earth_gravity_constant && !has_gravity_constant
        throw(
            IcgemParseError(
                "The gravity constant field is missing. Expected either `earth_gravity_constant` or `gravity_constant`.",
            ),
        )
    end

    product_type = Symbol(keywords[:product_type])
    model_name   = keywords[:modelname]
    max_degree   = tryparse(Int, keywords[:max_degree])
    errors       = Symbol(keywords[:errors])

    # Parse the gravity constant field (whichever one exists)
    gravity_constant_key =
        has_earth_gravity_constant ? :earth_gravity_constant : :gravity_constant
    gravity_constant = _parse_icgem_float(Tf, keywords[gravity_constant_key])
    radius = _parse_icgem_float(Tf, keywords[:radius])

    isnothing(gravity_constant) &&
        throw(IcgemParseError("Could not parse the gravity constant to $Tf."))
    isnothing(radius) && throw(IcgemParseError("Could not parse the radius to $Tf."))
    isnothing(max_degree) &&
        throw(IcgemParseError("Could not parse the maximum degree to an integer."))

    # Check if some keywords are valid.
    gravity_constant <= 0 &&
        throw(IcgemParseError("The gravity constant must be positive."))

    radius <= 0 && throw(IcgemParseError("The radius must be positive."))

    max_degree < 0 && throw(IcgemParseError("The maximum degree must not be negative."))

    errors ∉ (:no, :calibrated, :calibrated_and_formal, :formal) &&
        throw(IcgemParseError("An invalid value was found for the keyword `errors`."))

    # == Parse Optional Keywords ===========================================================

    tide_system = haskey(keywords, :tide_system) ? Symbol(keywords[:tide_system]) : :unknown
    norm_str    = haskey(keywords, :norm) ? keywords[:norm] : "fully_normalized"

    # Convert the normalization to the value expected by the Legendre functions.
    norm = if norm_str == "fully_normalized"
        Val(:full)
    elseif norm_str == "unnormalized"
        Val(:unnormalized)
    else
        throw(IcgemParseError("An invalid value was found for the keyword `norm`."))
    end

    # == Data ==============================================================================

    # Since we now have the maximum degree, we can pre-allocate and initialize the data
    # matrix. The storage must be filled with zeros because the file can omit coefficients,
    # such as those of degree 1, which are zero by definition.
    data = zeros(LowerTriangularStorage{RowMajor, IcgemGfcCoefficient{Tf}}, max_degree + 1)

    # Time-variable coefficients found in the file.
    time_variable_coefficients = IcgemTimeVariableCoefficient{Tf}[]

    # State of the parsing algorithm.
    state = :new

    # Auxiliary variables to build the time-variable coefficients.
    deg = 0
    ord = 0
    clm = Tf(0)
    slm = Tf(0)
    t₀  = Tf(0)
    t₁  = Tf(Inf)

    trend_clm      = Tf(0)
    trend_slm      = Tf(0)
    periodic_terms = IcgemPeriodicTerm{Tf}[]

    line          = nothing
    read_new_line = true
    tokens        = nothing

    # Flush the coefficient built from a `gfct` section and its subsequent lines to the
    # vector of time-variable coefficients.
    function flush_gfct_coefficient!()
        push!(
            time_variable_coefficients,
            IcgemTimeVariableCoefficient(
                deg, ord, clm, slm, t₀, t₁, trend_clm, trend_slm, copy(periodic_terms)
            ),
        )

        return nothing
    end

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
            # Every line processed in this state is consumed. Hence, we must read a new
            # line in the next iteration regardless of the outcome. Otherwise, we can
            # reach an infinite loop if this line was carried over from the end of a
            # `gfct` section.
            read_new_line = true

            if length(tokens) < 3
                @warn "[Line $current_line] Invalid data line."
                continue
            end

            # == `gfc` Data Line ===========================================================

            if tokens[1] == "gfc"
                ret = _parse_gfc_data_line(Tf, tokens, current_line)
                isnothing(ret) && continue

                deg, ord, clm, slm = ret
                _is_degree_and_order_valid(deg, ord, max_degree, current_line) || continue

                data[deg + 1, ord + 1] = IcgemGfcCoefficient(clm, slm)

                # == `gfct` Data Line ==========================================================

            elseif tokens[1] == "gfct"
                ret = _parse_gfct_data_line(Tf, tokens, current_line)
                isnothing(ret) && continue

                deg, ord, clm, slm, t₀ = ret
                _is_degree_and_order_valid(deg, ord, max_degree, current_line) || continue

                # Now, we need to change the state to wait for the next terms.
                state = :gfct

                t₁        = Tf(Inf)
                trend_clm = Tf(0)
                trend_slm = Tf(0)
                empty!(periodic_terms)
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

                ((adeg != deg) || (aord != ord)) && throw(
                    IcgemParseError(
                        "The degree or order of a `trnd` line is different from the corresponding `gfct` line.",
                        current_line,
                    ),
                )

                read_new_line = true

                # == `asin` and `acos` Data Lines of a `gfct` Section ==========================

            elseif (tokens[1] == "asin") || (tokens[1] == "acos")
                ret = _parse_asin_acos_data_line(Tf, tokens, current_line)
                if isnothing(ret)
                    read_new_line = true
                    continue
                end

                adeg, aord, amplitude_clm, amplitude_slm, period = ret

                ((adeg != deg) || (aord != ord)) && throw(
                    IcgemParseError(
                        "The degree or order of a `$(tokens[1])` line is different from the corresponding `gfct` line.",
                        current_line,
                    ),
                )

                _add_periodic_term!(
                    periodic_terms,
                    tokens[1] == "asin",
                    amplitude_clm,
                    amplitude_slm,
                    period,
                )

                read_new_line = true

            else
                # If we reach this part, the `gfct` section is over. Thus, we should create
                # the element related to `gfct` and proceed with the new information.
                flush_gfct_coefficient!()

                state = :new
                read_new_line = false
            end
        end
    end

    # If the file ended while we were processing a `gfct` section, we must flush the
    # pending coefficient. Otherwise, the last time-variable coefficient would be lost.
    state === :gfct && flush_gfct_coefficient!()

    # == Time-Variable Coefficients Index ==================================================

    # The time-variable coefficients are sorted by degree, order, and epoch so that the
    # validity intervals of the same coefficient are consecutive. The index maps the degree
    # and order to the first interval, and the static storage receives the values at the
    # epoch of the first interval.
    sort!(time_variable_coefficients; by = c -> (c.degree, c.order, c.t₀))

    max_time_variable_degree =
        isempty(time_variable_coefficients) ? -1 :
        maximum(c -> c.degree, time_variable_coefficients)

    time_variable_index = zeros(
        LowerTriangularStorage{RowMajor, Int32}, max(max_time_variable_degree + 1, 1)
    )

    for (k, c) in enumerate(time_variable_coefficients)
        (time_variable_index[c.degree + 1, c.order + 1] != 0) && continue

        time_variable_index[c.degree + 1, c.order + 1] = k
        data[c.degree + 1, c.order + 1] = IcgemGfcCoefficient(c.clm, c.slm)
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
        norm,
        data,
        max_time_variable_degree,
        time_variable_index,
        time_variable_coefficients,
    )

    return icgem_file
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _parse_icgem_float(::Type{T}, input::AbstractString) -> Union{Nothing, T}

Parse the `input` to the float type `T`, substituting all `D`s and `d`s by `e` so that
numbers in FORTRAN format can be converted. If `input` cannot be parsed to `T`, return
`nothing`.
"""
function _parse_icgem_float(::Type{T}, input::AbstractString) where {T}
    data_str = replace(input, r"[Dd]" => "e")
    return tryparse(T, data_str)
end

"""
    _is_degree_and_order_valid(degree::Int, order::Int, max_degree::Int, current_line::Int) -> Bool

Check if `degree` and `order` of a data line are valid, i.e. if `degree` is not higher than
`max_degree` and `order` is not higher than `degree`. If they are invalid, log a warning
with the `current_line` number and return `false`.
"""
function _is_degree_and_order_valid(
    degree::Int, order::Int, max_degree::Int, current_line::Int
)
    if (degree > max_degree) || (order > degree) || (order < 0)
        @warn "[Line $current_line] Invalid degree or order: $degree, $order."
        return false
    end

    return true
end

"""
    _add_periodic_term!(periodic_terms::Vector{IcgemPeriodicTerm{T}}, is_sine::Bool, amplitude_clm::T, amplitude_slm::T, period::T) -> Nothing

Add to `periodic_terms` the sine (`is_sine = true`) or cosine (`is_sine = false`) term with
the amplitudes `amplitude_clm` [-] and `amplitude_slm` [-] and the `period` [year]. If
`periodic_terms` already has a term with the same period, its amplitudes are updated.
Otherwise, a new term is pushed to the vector, which is modified in place.
"""
function _add_periodic_term!(
    periodic_terms::Vector{IcgemPeriodicTerm{T}},
    is_sine::Bool,
    amplitude_clm::T,
    amplitude_slm::T,
    period::T,
) where {T <: Number}
    k = findfirst(p -> p.period == period, periodic_terms)

    if isnothing(k)
        push!(periodic_terms, IcgemPeriodicTerm(zero(T), zero(T), zero(T), zero(T), period))
        k = lastindex(periodic_terms)
    end

    p = periodic_terms[k]

    periodic_terms[k] = if is_sine
        IcgemPeriodicTerm(
            amplitude_clm,
            amplitude_slm,
            p.amplitude_cos_clm,
            p.amplitude_cos_slm,
            period,
        )
    else
        IcgemPeriodicTerm(
            p.amplitude_sin_clm,
            p.amplitude_sin_slm,
            amplitude_clm,
            amplitude_slm,
            period,
        )
    end

    return nothing
end

# == Functions to Parse Data Lines =========================================================

"""
    _parse_degree_and_order(tokens, current_line) -> Union{Nothing, Tuple{Int, Int}}

Parse the degree in `tokens[2]` and the order in `tokens[3]`. If any of them cannot be
parsed, log a warning with the `current_line` number and return `nothing`.

# Arguments

- `tokens::AbstractVector{<:AbstractString}`: Tokens of the data line.
- `current_line::Int`: Number of the line being parsed, used in the warning messages.

# Returns

- `Int`: Degree.
- `Int`: Order.
"""
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

"""
    _parse_gfc_data_line(Tf, tokens, current_line) -> Union{Nothing, Tuple}

Parse the `gfc` data line in `tokens` using the type `Tf` for the floating point fields.
If any field cannot be parsed, log a warning with the `current_line` number and return
`nothing`.

# Arguments

- `Tf::Type`: Type used to parse the floating point fields.
- `tokens::AbstractVector{<:AbstractString}`: Tokens of the data line.
- `current_line::Int`: Number of the line being parsed, used in the warning messages.

# Returns

- `Int`: Degree.
- `Int`: Order.
- `Tf`: Coefficient `Clm` [-].
- `Tf`: Coefficient `Slm` [-].
"""
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

"""
    _parse_gfct_data_line(Tf, tokens, current_line) -> Union{Nothing, Tuple}

Parse the `gfct` data line in `tokens` using the type `Tf` for the floating point fields.
If any field cannot be parsed, log a warning with the `current_line` number and return
`nothing`. The function throws an `ArgumentError` if the epoch in the last token is not a
valid date in the `yyyymmdd` format.

# Arguments

- `Tf::Type`: Type used to parse the floating point fields.
- `tokens::AbstractVector{<:AbstractString}`: Tokens of the data line.
- `current_line::Int`: Number of the line being parsed, used in the warning messages.

# Returns

- `Int`: Degree.
- `Int`: Order.
- `Tf`: Coefficient `Clm` [-] at the epoch.
- `Tf`: Coefficient `Slm` [-] at the epoch.
- `Float64`: Epoch (`t₀`) of the coefficients, expressed as the number of elapsed seconds
    [s] since the J2000.0 epoch (2000-01-01T12:00:00).
"""
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
    time = _to_j2000_seconds(DateTime(tokens[end], dateformat"yyyymmdd"))

    return deg, ord, clm, slm, time
end

"""
    _parse_trnd_data_line(Tf, tokens, current_line) -> Union{Nothing, Tuple}

Parse the `trnd` data line in `tokens` using the type `Tf` for the floating point fields.
If any field cannot be parsed, log a warning with the `current_line` number and return
`nothing`.

# Arguments

- `Tf::Type`: Type used to parse the floating point fields.
- `tokens::AbstractVector{<:AbstractString}`: Tokens of the data line.
- `current_line::Int`: Number of the line being parsed, used in the warning messages.

# Returns

- `Int`: Degree.
- `Int`: Order.
- `Tf`: Linear trend of `Clm` [year⁻¹].
- `Tf`: Linear trend of `Slm` [year⁻¹].
"""
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

"""
    _parse_asin_acos_data_line(Tf, tokens, current_line) -> Union{Nothing, Tuple}

Parse the `asin` or `acos` data line in `tokens` using the type `Tf` for the floating
point fields. If any field cannot be parsed, log a warning with the `current_line` number
and return `nothing`.

# Arguments

- `Tf::Type`: Type used to parse the floating point fields.
- `tokens::AbstractVector{<:AbstractString}`: Tokens of the data line.
- `current_line::Int`: Number of the line being parsed, used in the warning messages.

# Returns

- `Int`: Degree.
- `Int`: Order.
- `Tf`: Amplitude of the periodic term for `Clm` [-].
- `Tf`: Amplitude of the periodic term for `Slm` [-].
- `Tf`: Period of the term [year].
"""
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
