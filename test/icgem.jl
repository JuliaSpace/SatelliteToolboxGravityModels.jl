# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Tests related to the ICGEM file support.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# File: ./src/icgem/fetch.jl
# ==========================================================================================

@testset "Fetching ICGEM files" verbose = true begin
    egm96_file = (@test_logs (
        :info,
        "Downloading the ICGEM file 'EGM96.gfc' from 'http://icgem.gfz-potsdam.de/getmodel/gfc/971b0a3b49a497910aad23cd85e066d4cd9af0aeafe7ce6301a696bed8570be3/EGM96.gfc'..."
    ) fetch_icgem_file(:EGM96))

    @test basename(egm96_file) == "EGM96.gfc"

    egm96_file_rerun = (@test_logs fetch_icgem_file(:EGM96))

    @test egm96_file_rerun == egm96_file
end

# Files: ./src/icgem/compute.jl
# ==========================================================================================

############################################################################################
#                                       Test Result
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
    # We will fetch the EIGEN-6C model that has time dependent coefficients.
    eigen6c_file = fetch_icgem_file(
        "http://icgem.gfz-potsdam.de/getmodel/gfc/0776caed6c65af24051697a65147b59e436cb464cb0930c1863fee6ecfbc31b0/EIGEN-6C.gfc"
    )

    eigen6c = GravityModels.load(Val(:ICGEM), eigen6c_file)

    Clm, Slm = GravityModels.coefficients(eigen6c, 2, 2, DateTime("2023-06-19"))

    @test Clm ≈ +2.4393378057597012e-6 atol = 1e-20
    @test Slm ≈ -1.400407403685511e-6  atol = 1e-20

    Clm, Slm = GravityModels.coefficients(eigen6c, 100, 1, DateTime("2023-06-19"))

    @test Clm ≈ -1.09755466854e-09 atol = 1e-20
    @test Slm ≈ +6.91287419630e-10 atol = 1e-20
end

# File: ./src/icgem/show.jl
# ==========================================================================================

@testset "Showing IcgemFile objects" verbose = true begin
    egm96 = GravityModels.load(Val(:ICGEM), fetch_icgem_file(:EGM96))

    expected = "ICGEM EGM96 (Degree = 360) {Float64}"
    result = sprint(show, egm96)
    @test result == expected

    expected = """
SatelliteToolboxGravityModels.IcgemFile{Float64}:
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
