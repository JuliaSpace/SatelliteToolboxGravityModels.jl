## Description #############################################################################
#
# Tests related to the ICGEM file support.
#
############################################################################################

# == File: ./src/icgem/api.jl ==============================================================

@testset "API Support" verbose = true begin
    @testset "Unnormalized Coefficients" begin
        model = GravityModels.load(
            IcgemFile, "./icgem_test_files/unnormalized_coefficients.gfc"
        )
        @test GravityModels.coefficient_norm(model) == Val(:unnormalized)
    end
end

# == File: ./src/icgem/parse.jl ============================================================

@testset "Parsing IcgemFile from a Stream" verbose = true begin
    filename = "./icgem_test_files/unnormalized_coefficients.gfc"

    model_from_file   = GravityModels.load(IcgemFile, filename)
    model_from_stream = open(io -> SatelliteToolboxGravityModels.parse_icgem(io), filename)

    @test model_from_stream.model_name == model_from_file.model_name
    @test model_from_stream.max_degree == model_from_file.max_degree

    for n in 0:model_from_file.max_degree, m in 0:n
        @test GravityModels.coefficients(model_from_stream, n, m) ==
            GravityModels.coefficients(model_from_file, n, m)
    end
end

# == File: ./src/icgem/fetch.jl ============================================================

@testset "Fetching ICGEM files" verbose = true begin
    egm96_file = (@test_logs (
        :info,
        "Downloading the ICGEM file 'EGM96.gfc' from 'https://icgem.gfz-potsdam.de/getmodel/gfc/971b0a3b49a497910aad23cd85e066d4cd9af0aeafe7ce6301a696bed8570be3/EGM96.gfc'...",
    ) fetch_icgem_file(:EGM96))

    @test basename(egm96_file) == "EGM96.gfc"

    egm96_file_rerun = (@test_logs fetch_icgem_file(:EGM96))

    @test egm96_file_rerun == egm96_file
end

# == Files: ./src/icgem/compute.jl =========================================================

############################################################################################
#                                       Test Result                                        #
############################################################################################
#
# The EIGEN-6C file has time dependent coefficients. According to the documentation, we must
# compute those coefficients using:
#
#   G(t)=gfct + trnd*(t-t0) + asin1*sin(2pi/p1*(t-t0))+acos1*cos(2pi/p1*(t-t0))
#                           + asin2*sin(2pi/p2*(t-t0))+acos2*cos(2pi/p2*(t-t0))
#
# Let's take the information about the degree 2 and order 2:
#
#   gfct   2    2  2.43935818007e-06 -1.40028528390e-06 1.8060e-13 1.7953e-13 20050101
#   trnd   2    2  2.64270248646e-13 -3.70169986147e-12 3.1405e-14 3.1373e-14
#   asin   2    2  1.02121464558e-11 -3.01068425667e-11 1.7054e-13 1.6988e-13 1.0
#   acos   2    2  1.77729432622e-11  4.65203150438e-11 1.6648e-13 1.6601e-13 1.0
#   asin   2    2 -4.59036913656e-12  3.73609063085e-12 1.6698e-13 1.6664e-13 0.5
#   acos   2    2 -1.14657993582e-11 -1.83016952664e-12 1.6564e-13 1.6506e-13 0.5
#
# Thus, the `C` coefficient must be:
#
#   C_2_2(Δt) = +2.43935818007e-06 +
#               +2.64270248646e-13 * Δt +
#               +1.02121464558e-11 * sin(2π * Δt) +
#               +1.77729432622e-11 * cos(2π * Δt) +
#               -4.59036913656e-12 * sin(4π * Δt) +
#               -1.14657993582e-11 * cos(4π * Δt)
#
# Hence, we obtain the following value for the day 2023-06-19:
#
#   C_2_2(18.46132785763176) = 2.4393402707210633e-6
#
# For `S` coefficient, we have:
#
#   S_2_2(Δt) = -1.40028528390e-06 +
#               -3.70169986147e-12 * Δt +
#               -3.01068425667e-11 * sin(2π * Δt) +
#               +4.65203150438e-11 * cos(2π * Δt) +
#               +3.73609063085e-12 * sin(4π * Δt) +
#               -1.83016952664e-12 * cos(4π * Δt)
#
# Hence, we obtain the following value for the day 2023-06-19:
#
#   S_2_2(18.46132785763176) = -1.4004093829317802e-6
#
# Finally, for the degree 100 and order 1, we have:
#
#   gfc  100    1 -1.09755466854e-09  6.91287419630e-10 1.5840e-12 1.5848e-12
#
# or `C_100_1 = -1.09755466854e-09` and `S_100_1 = 6.91287419630e-10`.
#
############################################################################################

@testset "Computing ICGEM coefficients" verbose = true begin
    dt_J2000 = DateTime("2000-01-01T12:00:00.000")

    # We will fetch the EIGEN-6C model that has time dependent coefficients.
    eigen6c_file = fetch_icgem_file(
        "https://icgem.gfz-potsdam.de/getmodel/gfc/0776caed6c65af24051697a65147b59e436cb464cb0930c1863fee6ecfbc31b0/EIGEN-6C.gfc",
    )

    eigen6c = GravityModels.load(IcgemFile, eigen6c_file)

    Clm, Slm = GravityModels.coefficients(eigen6c, 2, 2, DateTime("2023-06-19"))

    @test Clm ≈ +2.4393402707210633e-6 atol = 1e-20
    @test Slm ≈ -1.4004093829317802e-6 atol = 1e-20

    Clm, Slm = GravityModels.coefficients(eigen6c, 100, 1, DateTime("2023-06-19"))

    @test Clm ≈ -1.09755466854e-09 atol = 1e-20
    @test Slm ≈ +6.91287419630e-10 atol = 1e-20

    # Testing the version without the time parameter, which defaults to J2000.0 epoch.
    Clm_j2000, Slm_j2000 = GravityModels.coefficients(
        eigen6c, 2, 2, DateTime("2000-01-01T12:00:00")
    )
    @test Clm_j2000 ≈ 2.439363161505511e-6 atol = 1e-20
    @test Slm_j2000 ≈ -1.4002219857412463e-6 atol = 1e-20

    Clm, Slm = GravityModels.coefficients(eigen6c, 2, 2)
    @test Clm == Clm_j2000
    @test Slm == Slm_j2000

    time = Dates.value(DateTime("2023-06-19") - dt_J2000) / 1000

    Clm, Slm = GravityModels.coefficients(eigen6c, 2, 2, time)

    @test Clm ≈ +2.4393402707210633e-6 atol = 1e-20
    @test Slm ≈ -1.4004093829317802e-6 atol = 1e-20

    Clm, Slm = GravityModels.coefficients(eigen6c, 100, 1, time)

    @test Clm ≈ -1.09755466854e-09 atol = 1e-20
    @test Slm ≈ +6.91287419630e-10 atol = 1e-20

    # Testing the version without the time parameter, which defaults to J2000.0 epoch.
    Clm_j2000, Slm_j2000 = GravityModels.coefficients(eigen6c, 2, 2, 0)
    @test Clm_j2000 ≈ 2.439363161505511e-6 atol = 1e-20
    @test Slm_j2000 ≈ -1.4002219857412463e-6 atol = 1e-20

    Clm, Slm = GravityModels.coefficients(eigen6c, 2, 2)
    @test Clm == Clm_j2000
    @test Slm == Slm_j2000

    # The time-variable coefficients must be stored sparsely, merging the sine and cosine
    # terms with the same period.
    @test eigen6c.max_time_variable_degree == 50
    @test size(eigen6c.time_variable_index) == (51, 51)
    @test length(eigen6c.time_variable_coefficients) == 1323
    @test eigen6c.time_variable_index[3, 3] != 0
    @test eigen6c.time_variable_index[2, 1] == 0

    c = eigen6c.time_variable_coefficients[eigen6c.time_variable_index[3, 3]]
    @test c.degree == 2
    @test c.order == 2
    @test c.clm == +2.43935818007e-06
    @test c.slm == -1.40028528390e-06
    @test c.t₀ == Dates.value(DateTime("2005-01-01") - dt_J2000) / 1000
    @test c.t₁ == Inf
    @test c.trend_clm == +2.64270248646e-13
    @test c.trend_slm == -3.70169986147e-12
    @test length(c.periodic_terms) == 2
    @test c.periodic_terms[1].period == 1.0
    @test c.periodic_terms[1].amplitude_sin_clm == +1.02121464558e-11
    @test c.periodic_terms[1].amplitude_sin_slm == -3.01068425667e-11
    @test c.periodic_terms[1].amplitude_cos_clm == +1.77729432622e-11
    @test c.periodic_terms[1].amplitude_cos_slm == +4.65203150438e-11
    @test c.periodic_terms[2].period == 0.5
    @test c.periodic_terms[2].amplitude_sin_clm == -4.59036913656e-12
    @test c.periodic_terms[2].amplitude_sin_slm == +3.73609063085e-12
    @test c.periodic_terms[2].amplitude_cos_clm == -1.14657993582e-11
    @test c.periodic_terms[2].amplitude_cos_slm == -1.83016952664e-12

    # The static storage must contain the value at the epoch.
    @test eigen6c.data[3, 3].clm == c.clm
    @test eigen6c.data[3, 3].slm == c.slm

    # The coefficients omitted in the file must be zero.
    egm96 = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

    for m in 0:1
        Clm, Slm = GravityModels.coefficients(egm96, 1, m)
        @test Clm == 0
        @test Slm == 0
    end
end

# == File: ./src/icgem/parse.jl ============================================================

@testset "Parsing IcgemFile in the Format 2.0" verbose = true begin
    dt_J2000 = DateTime("2000-01-01T12:00:00.000")
    model    = GravityModels.load(IcgemFile, "./icgem_test_files/icgem2_time_variable.gfc")

    # The header value of `errors` is followed by a comment.
    @test model.errors == :formal
    @test model.max_degree == 2
    @test model.max_time_variable_degree == 2
    @test length(model.time_variable_coefficients) == 3

    # The validity intervals must be sorted by epoch.
    c₁ = model.time_variable_coefficients[model.time_variable_index[3, 1]]
    c₂ = model.time_variable_coefficients[model.time_variable_index[3, 1] + 1]
    c₃ = model.time_variable_coefficients[model.time_variable_index[3, 2]]

    @test (c₁.degree, c₁.order) == (2, 0)
    @test c₁.t₀ == Dates.value(DateTime("2000-01-01T12:00:00") - dt_J2000) / 1000
    @test c₁.t₁ == Dates.value(DateTime("2010-01-01T00:00:00") - dt_J2000) / 1000
    @test (c₂.degree, c₂.order) == (2, 0)
    @test c₂.t₀ == Dates.value(DateTime("2010-01-01T00:00:00") - dt_J2000) / 1000
    @test c₂.t₁ == Dates.value(DateTime("2020-01-01T00:00:00") - dt_J2000) / 1000
    @test (c₃.degree, c₃.order) == (2, 1)
    @test c₃.t₀ == Dates.value(DateTime("2020-01-01T00:00:00") - dt_J2000) / 1000
    @test c₃.t₁ == Dates.value(DateTime("2030-01-01T00:00:00") - dt_J2000) / 1000

    # The static storage must contain the values of the first interval.
    @test model.data[3, 1].clm == -1.0e-3
    @test model.data[3, 2].clm == +7.0e-4

    # == Coefficients Inside the Validity Intervals ========================================

    expected(c, t) = begin
        Δt  = (t - c.t₀) / (86400 * 365.25)
        clm = c.clm + c.trend_clm * Δt
        slm = c.slm + c.trend_slm * Δt

        for p in c.periodic_terms
            clm +=
                p.amplitude_sin_clm * sin(2π * Δt / p.period) +
                p.amplitude_cos_clm * cos(2π * Δt / p.period)
            slm +=
                p.amplitude_sin_slm * sin(2π * Δt / p.period) +
                p.amplitude_cos_slm * cos(2π * Δt / p.period)
        end

        return clm, slm
    end

    for (c, date) in (
        (c₁, DateTime("2005-06-19T03:00:00")),
        (c₂, DateTime("2015-06-19T03:00:00")),
        # Before the first interval, the first coefficient is used.
        (c₁, DateTime("1990-06-19T03:00:00")),
        # After the last interval, the last coefficient is used.
        (c₂, DateTime("2035-06-19T03:00:00")),
    )
        t = Dates.value(date - dt_J2000) / 1000
        Clm, Slm = GravityModels.coefficients(model, 2, 0, date)
        Clm_e, Slm_e = expected(c, t)
        @test Clm ≈ Clm_e atol = 1e-20
        @test Slm ≈ Slm_e atol = 1e-20
    end

    # The boundary of an interval belongs to the next one.
    t = c₂.t₀
    Clm, Slm = GravityModels.coefficients(model, 2, 0, t)
    Clm_e, Slm_e = expected(c₂, t)
    @test Clm ≈ Clm_e atol = 1e-20
    @test Slm ≈ Slm_e atol = 1e-20

    t = Dates.value(DateTime("2025-01-01") - dt_J2000) / 1000
    Clm, Slm = GravityModels.coefficients(model, 2, 1, t)
    Clm_e, Slm_e = expected(c₃, t)
    @test Clm ≈ Clm_e atol = 1e-20
    @test Slm ≈ Slm_e atol = 1e-20

    # The constant coefficients are not affected.
    @test GravityModels.coefficients(model, 2, 2, t) == (1.0e-6, 2.0e-6)

    # == Printing ==========================================================================

    expected_str = """
SatelliteToolboxGravityModels.IcgemTimeVariableCoefficient{Float64}:
    Degree : 2
     Order : 1
      Clm₀ : 0.0007
      Slm₀ : 0.0008
     Epoch : 2020-01-01T00:00:00
  Valid to : 2030-01-01T00:00:00
     Trend : Clm = 0.0, Slm = 0.0
  Periodic : Period 0.5 y => Sine: Clm = 0.0, Slm = 0.0; Cosine: Clm = 1.0e-6, Slm = 2.0e-6"""

    @test sprint(show, MIME("text/plain"), c₃) == expected_str
end

@testset "Parsing IcgemFile [ERRORS]" verbose = true begin
    @test sprint(showerror, IcgemParseError("Message.")) == "IcgemParseError: Message."
    @test sprint(showerror, IcgemParseError("Message.", 8)) ==
        "IcgemParseError: [Line 8] Message."

    @test_throws(
        IcgemParseError("Two `begin_of_head` keywords were found.", 8),
        GravityModels.load(IcgemFile, "./icgem_test_files/two_begin_of_head.gfc")
    )

    @test_throws(
        IcgemParseError("The mandatory keyword `end_of_head` was not found."),
        GravityModels.load(IcgemFile, "./icgem_test_files/no_end_of_head.gfc")
    )

    @test_throws(
        IcgemParseError(
            "The following mandatory fields are missing: (:radius, :max_degree)."
        ),
        GravityModels.load(IcgemFile, "./icgem_test_files/missing_mandatory_fields.gfc")
    )

    @test_logs(
        (:warn, "[Line 18] Invalid data line."),
        GravityModels.load(IcgemFile, "./icgem_test_files/invalid_data_line.gfc")
    )

    @test_logs(
        (:warn, "[Line 18] Invalid degree: 2a."),
        (:warn, "[Line 19] Invalid order: 1a."),
        GravityModels.load(IcgemFile, "./icgem_test_files/invalid_degree_and_order.gfc")
    )

    @test_logs(
        (:warn, "[Line 17] Invalid `gfc` data line."),
        (:warn, "[Line 18] Could not parse `Clm` to Float64: -0.18a987635955e-09."),
        (:warn, "[Line 19] Could not parse `Slm` to Float64: -0.1400a6683654e-05."),
        GravityModels.load(IcgemFile, "./icgem_test_files/invalid_gfc_data_lines.gfc")
    )

    @test_logs(
        (:warn, "[Line 22] Invalid `gfct` data line."),
        (:warn, "[Line 28] Could not parse `Clm` to Float64: 9.57a11211877e-07."),
        (:warn, "[Line 34] Could not parse `Slm` to Float64: 0.00b000000000e+00."),
        (:warn, "[Line 41] Could not parse `trend_C` to Float64: -5.03a51696812e-12."),
        (
            :warn,
            "[Line 42] Could not parse `Clm` amplitude to Float64: -1.05b85537206e-10.",
        ),
        (:warn, "[Line 43] Could not parse `Clm` amplitude to Float64: 5.08c62560512e-12."),
        (:warn, "[Line 48] Invalid `asin` or `acos` data line."),
        (:warn, "[Line 49] Invalid `asin` or `acos` data line."),
        (:warn, "[Line 53] Invalid `trnd` data line."),
        (:warn, "[Line 59] Could not parse `trend_S` to Float64: 0.0b0000000000e+00."),
        (
            :warn,
            "[Line 86] Could not parse `Slm` amplitude to Float64: -1.07328392828e-12.",
        ),
        (:warn, "[Line 90] Could not parse period to Float64: -3.72637514028e-12."),
        GravityModels.load(IcgemFile, "./icgem_test_files/invalid_gfct_data_lines.gfc")
    )
end

# == File: ./src/icgem/show.jl =============================================================

@testset "Showing IcgemFile" verbose = true begin
    egm96 = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

    expected = "ICGEM EGM96 (Degree = 360) {Float64}"
    result = sprint(show, egm96)
    @test result == expected

    expected = """
IcgemFile{Float64, :full}:
      Product type : gravity_field
       Model name  : EGM96
  Gravity constant : 3.986004415e14
            Radius : 6.3781363e6
    Maximum degree : 360
            Errors : formal
       Tide system : tide_free
              Norm : full
         Data type : Float64"""

    result = sprint(show, MIME("text/plain"), egm96)

    @test result == expected
end

@testset "Showing IcgemGfcCoefficient" verbose = true begin
    egm96 = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

    expected = "SatelliteToolboxGravityModels.IcgemGfcCoefficient{Float64}(Clm = -0.000484165371736, Slm = 0.0)"
    result = sprint(show, egm96.data[3, 1])
    @test result == expected

    expected = """
SatelliteToolboxGravityModels.IcgemGfcCoefficient{Float64}:
  Clm : -0.000484165371736
  Slm : 0.0"""

    result = sprint(show, MIME("text/plain"), egm96.data[3, 1])
    @test result == expected
end

@testset "Showing IcgemTimeVariableCoefficient" verbose = true begin
    eigen6c_file = fetch_icgem_file(
        "https://icgem.gfz-potsdam.de/getmodel/gfc/0776caed6c65af24051697a65147b59e436cb464cb0930c1863fee6ecfbc31b0/EIGEN-6C.gfc",
    )

    eigen6c = GravityModels.load(IcgemFile, eigen6c_file)
    c = eigen6c.time_variable_coefficients[eigen6c.time_variable_index[3, 1]]

    expected = "SatelliteToolboxGravityModels.IcgemTimeVariableCoefficient{Float64}(2, 0, Clm₀ = -0.000484165299806, Slm₀ = 0.0)"
    result = sprint(show, c)
    @test result == expected

    expected = """
SatelliteToolboxGravityModels.IcgemTimeVariableCoefficient{Float64}:
    Degree : 2
     Order : 0
      Clm₀ : -0.000484165299806
      Slm₀ : 0.0
     Epoch : 2005-01-01T00:00:00
     Trend : Clm = -1.26060242677e-11, Slm = 0.0
  Periodic : Period 1.0 y => Sine: Clm = 5.32328946063e-11, Slm = 0.0; Cosine: Clm = 4.10012162817e-11, Slm = 0.0
             Period 0.5 y => Sine: Clm = -2.44339926664e-11, Slm = 0.0; Cosine: Clm = 3.33917546745e-11, Slm = 0.0"""

    result = sprint(show, MIME("text/plain"), c)
    @test result == expected
end
