# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Functions to show types.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

function show(io::IO, c::IcgemGfcCoefficient{T}) where T
    print(io, typeof(c), "(Clm = ", c.clm, ", Slm = ", c.slm, ")")
    return nothing
end

function show(io::IO, c::IcgemGfctCoefficient{T}) where T
    print(io, typeof(c), "(Clm₀ = ", c.clm, ", Slm₀ = ", c.slm, ")")
    return nothing
end

function show(io::IO, mime::MIME"text/plain", c::IcgemGfcCoefficient{T}) where T
    # Check for color support in the `io`.
    color = get(io, :color, false)
    b = color ? string(_B) : ""
    d = color ? string(_D) : ""

    println(io, typeof(c), ":")
    println(io, "$(b)  Clm :$(d) ", c.clm)
    print(  io, "$(b)  Slm :$(d) ", c.slm)

    return nothing
end

function show(io::IO, mime::MIME"text/plain", c::IcgemGfctCoefficient{T}) where T
    # Check for color support in the `io`.
    color = get(io, :color, false)
    b = color ? string(_B) : ""
    d = color ? string(_D) : ""

    println(io, typeof(c), ":")
    println(io, "$(b)    Clm₀ :$(d) ", c.clm)
    println(io, "$(b)    Slm₀ :$(d) ", c.slm)
    println(io, "$(b)   Epoch :$(d) ", c.time)
    println(io, "$(b)   Trend :$(d) Clm = ", c.trend_clm, ", Slm = ", c.trend_slm)
    print(  io, "$(b)    Sine :$(d) ")
    _print_asin_acos_vectors(io, c.asin_coefficients)
    println(io)
    print(  io, "$(b)  Cosine : $(d)")
    _print_asin_acos_vectors(io, c.acos_coefficients)

    return nothing
end

function show(io::IO, m::IcgemFile{T}) where T
    # Check for color support in the `io`.
    color = get(io, :color, false)
    b = color ? string(_B) : ""
    d = color ? string(_D) : ""

    print(io, "$(b)ICGEM ", m.model_name, "$(d) (Degree = ", m.max_degree, ") {", string(T), "}")
    return nothing
end

function show(io::IO, mime::MIME"text/plain", m::IcgemFile{T}) where T
    # Check for color support in the `io`.
    color = get(io, :color, false)
    b = color ? string(_B) : ""
    d = color ? string(_D) : ""

    println(io, typeof(m), ":")
    println(io, "$(b)      Product type :$(d) ", m.product_type)
    println(io, "$(b)       Model name  :$(d) ", m.model_name)
    println(io, "$(b)  Gravity constant :$(d) ", m.gravity_constant)
    println(io, "$(b)            Radius :$(d) ", m.radius)
    println(io, "$(b)    Maximum degree :$(d) ", m.max_degree)
    println(io, "$(b)            Errors :$(d) ", m.errors)
    println(io, "$(b)       Tide system :$(d) ", m.tide_system)
    println(io, "$(b)              Norm :$(d) ", m.norm)
    print(  io, "$(b)         Data type :$(d) ", string(T))

    return nothing
end

############################################################################################
#                                    Private Functions
############################################################################################

# Print the vectors with the coefficients related to the sine and cosine terms.
function _print_asin_acos_vectors(io::IO, v::Vector{NTuple{3, T}}) where T
    num_coefficients = length(v)
    for k in 1:num_coefficients
        c = v[k]
        print(io, "Period ", c[3], " y => ", "Amp. Clm = ", c[1], ", Amp. Slm = ", c[2])

        if k != num_coefficients
            println(io)
            print(io, "           ")
        end
    end

    return nothing
end
