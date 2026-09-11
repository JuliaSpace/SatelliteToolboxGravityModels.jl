## Description #############################################################################
#
# Functions to show the ICGEM types using the representation helpers of
# SatelliteToolboxBase.jl.
#
############################################################################################

# == IcgemGfcCoefficient ===================================================================

function show(io::IO, c::IcgemGfcCoefficient)
    print(io, type_name(c), "(Clm = ", c.clm, ", Slm = ", c.slm, ")")
    return nothing
end

function show(io::IO, ::MIME"text/plain", c::IcgemGfcCoefficient)
    fields = PrintedField[("Clm", string(c.clm), ""), ("Slm", string(c.slm), "")]
    print_tree(io, type_name(c), fields, PrintedSection[])
    return nothing
end

# == IcgemTimeVariableCoefficient ==========================================================

function show(io::IO, c::IcgemTimeVariableCoefficient)
    print(
        io,
        type_name(c),
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

function show(io::IO, ::MIME"text/plain", c::IcgemTimeVariableCoefficient)
    fields = PrintedField[
        ("Degree", string(c.degree), ""),
        ("Order", string(c.order), ""),
        ("Clm₀", string(c.clm), ""),
        ("Slm₀", string(c.slm), ""),
        ("Epoch", string(_from_j2000_seconds(c.t₀)), ""),
    ]

    # The end of the validity interval is printed only if the file defines it.
    isfinite(c.t₁) && push!(fields, ("Valid To", string(_from_j2000_seconds(c.t₁)), ""))

    push!(fields, ("Trend Clm", string(c.trend_clm), "1/year"))
    push!(fields, ("Trend Slm", string(c.trend_slm), "1/year"))

    sections = PrintedSection[]

    for p in c.periodic_terms
        push!(
            sections,
            PrintedSection(
                "Periodic Term (Period = " * string(p.period) * " year)",
                PrintedField[
                    ("Sine Clm", string(p.amplitude_sin_clm), ""),
                    ("Sine Slm", string(p.amplitude_sin_slm), ""),
                    ("Cosine Clm", string(p.amplitude_cos_clm), ""),
                    ("Cosine Slm", string(p.amplitude_cos_slm), ""),
                ],
            ),
        )
    end

    print_tree(io, type_name(c), fields, sections)

    return nothing
end

# == IcgemFile =============================================================================

function show(io::IO, m::IcgemFile)
    print(io, type_name(m), ": ", m.model_name, " (Degree = ", m.max_degree, ")")
    return nothing
end

function show(io::IO, ::MIME"text/plain", m::IcgemFile{T, Val{N}}) where {T, N}
    time_variable = if isempty(m.time_variable_coefficients)
        "none"
    else
        string(
            length(m.time_variable_coefficients),
            " up to degree ",
            m.max_time_variable_degree,
        )
    end

    fields = PrintedField[
        ("Product Type", string(m.product_type), ""),
        ("Model Name", m.model_name, ""),
        ("Gravity Constant", format_value(m.gravity_constant), "m³/s²"),
        ("Radius", format_value(m.radius), "m"),
        ("Angular Speed", format_value(m.angular_speed), "rad/s"),
        ("Maximum Degree", string(m.max_degree), ""),
        ("Errors", string(m.errors), ""),
        ("Tide System", string(m.tide_system), ""),
        ("Normalization", string(N), ""),
        ("Time-Variable Coefficients", time_variable, ""),
    ]

    print_tree(io, type_name(m), fields, PrintedSection[])

    return nothing
end
