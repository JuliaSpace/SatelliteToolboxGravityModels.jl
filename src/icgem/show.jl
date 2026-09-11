## Description #############################################################################
#
# Functions to show types.
#
############################################################################################

function show(io::IO, c::IcgemGfcCoefficient{T}) where {T}
    print(io, typeof(c), "(Clm = ", c.clm, ", Slm = ", c.slm, ")")
    return nothing
end

function show(io::IO, mime::MIME"text/plain", c::IcgemGfcCoefficient{T}) where {T}
    # Check for color support in the `io`.
    color = get(io, :color, false)
    b = color ? _B : ""
    d = color ? _D : ""

    println(io, typeof(c), ":")
    println(io, "$(b)  Clm :$(d) ", c.clm)
    print(io, "$(b)  Slm :$(d) ", c.slm)

    return nothing
end

function show(io::IO, c::IcgemTimeVariableCoefficient{T}) where {T}
    print(
        io,
        typeof(c),
        "(",
        c.degree,
        ", ",
        c.order,
        ", Clm₀ = ",
        c.clm,
        ", Slm₀ = ",
        c.slm,
        ")",
    )
    return nothing
end

function show(io::IO, mime::MIME"text/plain", c::IcgemTimeVariableCoefficient{T}) where {T}
    # Check for color support in the `io`.
    color = get(io, :color, false)
    b = color ? _B : ""
    d = color ? _D : ""

    println(io, typeof(c), ":")
    println(io, "$(b)    Degree :$(d) ", c.degree)
    println(io, "$(b)     Order :$(d) ", c.order)
    println(io, "$(b)      Clm₀ :$(d) ", c.clm)
    println(io, "$(b)      Slm₀ :$(d) ", c.slm)
    println(io, "$(b)     Epoch :$(d) ", _from_j2000_seconds(c.t₀))
    isfinite(c.t₁) && println(io, "$(b)  Valid to :$(d) ", _from_j2000_seconds(c.t₁))
    println(io, "$(b)     Trend :$(d) Clm = ", c.trend_clm, ", Slm = ", c.trend_slm)
    print(io, "$(b)  Periodic :$(d) ")
    _print_periodic_terms(io, c.periodic_terms)

    return nothing
end

function show(io::IO, m::IcgemFile{T}) where {T}
    # Check for color support in the `io`.
    color = get(io, :color, false)
    b = color ? _B : ""
    d = color ? _D : ""

    print(
        io,
        "$(b)ICGEM ",
        m.model_name,
        "$(d) (Degree = ",
        m.max_degree,
        ") {",
        string(T),
        "}",
    )
    return nothing
end

function show(io::IO, mime::MIME"text/plain", m::IcgemFile{T, Val{N}}) where {T, N}
    # Check for color support in the `io`.
    color = get(io, :color, false)
    b = color ? _B : ""
    d = color ? _D : ""

    println(io, "IcgemFile{", T, ", :", N, "}:")
    println(io, "$(b)      Product type :$(d) ", m.product_type)
    println(io, "$(b)       Model name  :$(d) ", m.model_name)
    println(io, "$(b)  Gravity constant :$(d) ", m.gravity_constant)
    println(io, "$(b)            Radius :$(d) ", m.radius)
    println(io, "$(b)     Angular speed :$(d) ", m.angular_speed)
    println(io, "$(b)    Maximum degree :$(d) ", m.max_degree)
    println(io, "$(b)            Errors :$(d) ", m.errors)
    println(io, "$(b)       Tide system :$(d) ", m.tide_system)
    println(io, "$(b)              Norm :$(d) ", N)
    print(io, "$(b)         Data type :$(d) ", string(T))

    return nothing
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _print_periodic_terms(io::IO, v::Vector{IcgemPeriodicTerm{T}}) -> Nothing

Print to `io` the periodic terms in `v`, one per line, with the continuation lines
indented to align with the first one.
"""
function _print_periodic_terms(io::IO, v::Vector{IcgemPeriodicTerm{T}}) where {T}
    isempty(v) && (print(io, "none"); return nothing)

    for (k, p) in enumerate(v)
        (k != 1) && print(io, "\n             ")
        print(
            io,
            "Period ",
            p.period,
            " y => Sine: Clm = ",
            p.amplitude_sin_clm,
            ", Slm = ",
            p.amplitude_sin_slm,
            "; Cosine: Clm = ",
            p.amplitude_cos_clm,
            ", Slm = ",
            p.amplitude_cos_slm,
        )
    end

    return nothing
end
