## Description #############################################################################
#
# Tests related to the ICGEM file support.
#
############################################################################################

# == File: ./src/icgem/api.jl ==============================================================

@testset "API Support" verbose = true begin
    @testset "Unnormalized Coefficients" begin
        model = GravityModels.load(IcgemFile, "./icgem_test_files/unnormalized_coefficients.gfc")
        @test GravityModels.coefficient_norm(model) == :unnormalized
    end
end

# == File: ./src/icgem/fetch.jl ============================================================

@testset "Fetching ICGEM files" verbose = true begin
    egm96_file = (@test_logs (
        :info,
        "Downloading the ICGEM file 'EGM96.gfc' from 'http://icgem.gfz-potsdam.de/getmodel/gfc/971b0a3b49a497910aad23cd85e066d4cd9af0aeafe7ce6301a696bed8570be3/EGM96.gfc'..."
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
#   C_2_2(18.473972602739725) = 2.4393378057597012e-6
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
#   S_2_2(18.473972602739725) = -1.400407403685511e-6
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
        "http://icgem.gfz-potsdam.de/getmodel/gfc/0776caed6c65af24051697a65147b59e436cb464cb0930c1863fee6ecfbc31b0/EIGEN-6C.gfc"
    )

    eigen6c = GravityModels.load(IcgemFile, eigen6c_file)

    Clm, Slm = GravityModels.coefficients(eigen6c, 2, 2, DateTime("2023-06-19"))

    @test Clm ≈ +2.4393378057597012e-6 atol = 1e-20
    @test Slm ≈ -1.400407403685511e-6  atol = 1e-20

    Clm, Slm = GravityModels.coefficients(eigen6c, 100, 1, DateTime("2023-06-19"))

    @test Clm ≈ -1.09755466854e-09 atol = 1e-20
    @test Slm ≈ +6.91287419630e-10 atol = 1e-20

    # Testing the version without the time parameter, which defaults to J2000.0 epoch.
    Clm_j2000, Slm_j2000 = GravityModels.coefficients(
        eigen6c,
        2,
        2,
        DateTime("2000-01-01T12:00:00")
    )
    @test Clm_j2000 ≈ 2.4393631474296326e-6 atol = 1e-20
    @test Slm_j2000 ≈ -1.4002214986544143e-6 atol = 1e-20

    Clm, Slm = GravityModels.coefficients(eigen6c, 2, 2)
    @test Clm == Clm_j2000
    @test Slm == Slm_j2000

    time = Dates.value(DateTime("2023-06-19") - dt_J2000) / 1000

    Clm, Slm = GravityModels.coefficients(eigen6c, 2, 2, time)

    @test Clm ≈ +2.4393378057597012e-6 atol = 1e-20
    @test Slm ≈ -1.400407403685511e-6  atol = 1e-20

    Clm, Slm = GravityModels.coefficients(eigen6c, 100, 1, time)

    @test Clm ≈ -1.09755466854e-09 atol = 1e-20
    @test Slm ≈ +6.91287419630e-10 atol = 1e-20

    # Testing the version without the time parameter, which defaults to J2000.0 epoch.
    Clm_j2000, Slm_j2000 = GravityModels.coefficients(eigen6c, 2, 2, 0)
    @test Clm_j2000 ≈  2.4393631474296326e-6 atol = 1e-20
    @test Slm_j2000 ≈ -1.4002214986544143e-6 atol = 1e-20

    Clm, Slm = GravityModels.coefficients(eigen6c, 2, 2)
    @test Clm == Clm_j2000
    @test Slm == Slm_j2000
end

# == File: ./src/icgem/parse.jl ============================================================

@testset "Parsing IcgemFile [ERRORS]" verbose = true begin
    @test_throws(
        ErrorException("[Invalid ICGEM file] Two `begin_of_head` keywords were found!"),
        GravityModels.load(IcgemFile, "./icgem_test_files/two_begin_of_head.gfc")
    )

    @test_throws(
        ErrorException("[Invalid ICGEM file] The mandatory keyword `end_of_head` was not found!"),
        GravityModels.load(IcgemFile, "./icgem_test_files/no_end_of_head.gfc")
    )

    @test_throws(
        ErrorException("[Invalid ICGEM file] The following mandatory fields are missing: (:radius, :max_degree)."),
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
        (:warn, "[Line 42] Could not parse `Clm` amplitude to Float64: -1.05b85537206e-10."),
        (:warn, "[Line 43] Could not parse `Clm` amplitude to Float64: 5.08c62560512e-12."),
        (:warn, "[Line 48] Invalid `asin` or `acos` data line."),
        (:warn, "[Line 49] Invalid `asin` or `acos` data line."),
        (:warn, "[Line 53] Invalid `trnd` data line."),
        (:warn, "[Line 59] Could not parse `trend_S` to Float64: 0.0b0000000000e+00."),
        (:warn, "[Line 86] Could not parse `Slm` amplitude to Float64: -1.07328392828e-12."),
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
IcgemFile{Float64, :fully_normalized}:
      Product type : gravity_field
       Model name  : EGM96
  Gravity constant : 3.986004415e14
            Radius : 6.3781363e6
    Maximum degree : 360
            Errors : formal
       Tide system : tide_free
              Norm : fully_normalized
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

@testset "Showing IcgemGfctCoefficients" verbose = true begin
    eigen6c_file = fetch_icgem_file(
        "http://icgem.gfz-potsdam.de/getmodel/gfc/0776caed6c65af24051697a65147b59e436cb464cb0930c1863fee6ecfbc31b0/EIGEN-6C.gfc"
    )

    eigen6c = GravityModels.load(IcgemFile, eigen6c_file)

    expected = "SatelliteToolboxGravityModels.IcgemGfctCoefficient{Float64}(Clm₀ = -0.000484165299806, Slm₀ = 0.0)"
    result = sprint(show, eigen6c.data[3, 1])
    @test result == expected

    expected = """
SatelliteToolboxGravityModels.IcgemGfctCoefficient{Float64}:
    Clm₀ : -0.000484165299806
    Slm₀ : 0.0
   Epoch : 2005-01-01T00:00:00
   Trend : Clm = -1.26060242677e-11, Slm = 0.0
    Sine : Period 1.0 y => Amp. Clm = 5.32328946063e-11, Amp. Slm = 0.0
           Period 0.5 y => Amp. Clm = -2.44339926664e-11, Amp. Slm = 0.0
  Cosine : Period 1.0 y => Amp. Clm = 4.10012162817e-11, Amp. Slm = 0.0
           Period 0.5 y => Amp. Clm = 3.33917546745e-11, Amp. Slm = 0.0"""

    result = sprint(show, MIME("text/plain"), eigen6c.data[3, 1])
    @test result == expected
end
