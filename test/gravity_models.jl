## Description #############################################################################
#
# Tests related to the GravityModel API.
#
############################################################################################

# == File: ./src/GravityModels/api.jl ======================================================

@testset "API" verbose = true begin
    dt_J2000 = DateTime("2000-01-01T12:00:00.000")

    # We will fetch the EIGEN-6C model that has time dependent coefficients.
    eigen6c_file = fetch_icgem_file(
        "http://icgem.gfz-potsdam.de/getmodel/gfc/0776caed6c65af24051697a65147b59e436cb464cb0930c1863fee6ecfbc31b0/EIGEN-6C.gfc"
    )

    eigen6c = GravityModels.load(IcgemFile, eigen6c_file)

    time = Dates.value(DateTime("2023-06-19") - dt_J2000) / 1000

    Clm_expected, Slm_expected = GravityModels.coefficients(eigen6c, 2, 2, time)
    Clm, Slm = GravityModels.coefficients(eigen6c, 2, 2, DateTime("2023-06-19"))

    @test Clm == Clm_expected
    @test Slm == Slm_expected

    Clm_expected, Slm_expected = GravityModels.coefficients(eigen6c, 2, 2, 0)
    Clm, Slm = GravityModels.coefficients(eigen6c, 2, 2)

    @test Clm == Clm_expected
    @test Slm == Slm_expected
end

# == File: ./src/GravityModels/accelerations.jl ============================================

# TODO: Add tests with a model that has time-varying coefficients.

@testset "Gravitational Acceleration" verbose = true begin
    tests = (
        (:EGM96, "./test_results/gravitation/EGM96.gdf"),
        (:JGM2,  "./test_results/gravitation/JGM2.gdf"),
        (:JGM3,  "./test_results/gravitation/JGM3.gdf")
    )

    for t in tests
        @eval @testset $(string(t[1])) begin
            test_results = readdlm($t[2]; skipstart = 34)
            model = GravityModels.load(IcgemFile, fetch_icgem_file($t[1]))

            # Compare the results.
            for k in axes(test_results, 1)
                # Get latitude and longitude.
                lat = test_results[k, 1 + begin] |> deg2rad
                lon = test_results[k, 0 + begin] |> deg2rad

                # Get the expected result in m / s² (the online result is in mGal).
                expected_g_norm = test_results[k, 2 + begin] / 100000

                # Use the model to compute the gravity using all the coefficients.
                r_itrf = geodetic_to_ecef(lat, lon, 0)
                g_norm = norm(GravityModels.gravitational_acceleration(model, r_itrf))

                # Compare the results.
                # TODO: Check why the precision is worse in the poles.
                if abs(lat) ≈ π/2
                    @test g_norm ≈ expected_g_norm atol = 5e-8
                else
                    @test g_norm ≈ expected_g_norm atol = 1e-13
                end
            end
        end
    end

    # == Test the Signature with DateTime ==================================================

    dt_J2000 = DateTime("2000-01-01T12:00:00.000")

    # We will fetch the EIGEN-6C model that has time dependent coefficients.
    eigen6c_file = fetch_icgem_file(
        "http://icgem.gfz-potsdam.de/getmodel/gfc/0776caed6c65af24051697a65147b59e436cb464cb0930c1863fee6ecfbc31b0/EIGEN-6C.gfc"
    )

    eigen6c = GravityModels.load(IcgemFile, eigen6c_file)
    r_itrf = [7000.0e3, 0, 0]

    time = Dates.value(DateTime("2023-06-19") - dt_J2000) / 1000
    g_itrf_expected = GravityModels.gravitational_acceleration(eigen6c, r_itrf, time)
    g_itrf = GravityModels.gravitational_acceleration(eigen6c, r_itrf, DateTime("2023-06-19"))

    @test g_itrf == g_itrf_expected
end

@testset "Gravity Acceleration" verbose = true begin
    tests = (
        (:EGM96, "./test_results/gravity/EGM96.gdf"),
        (:JGM2,  "./test_results/gravity/JGM2.gdf"),
        (:JGM3,  "./test_results/gravity/JGM3.gdf")
    )

    for t in tests
        @eval @testset $(string(t[1])) begin
            test_results = readdlm($t[2]; skipstart = 34)
            model = GravityModels.load(IcgemFile, fetch_icgem_file($t[1]))

            # Compare the results.
            for k in axes(test_results, 1)
                # Get latitude and longitude.
                lat = test_results[k, 1 + begin] |> deg2rad
                lon = test_results[k, 0 + begin] |> deg2rad

                # Get the expected result in m / s² (the online result is in mGal).
                expected_g_norm = test_results[k, 2 + begin] / 100000

                # Use the model to compute the gravity using all the coefficients.
                r_itrf = geodetic_to_ecef(lat, lon, 0)
                g_norm = norm(GravityModels.gravity_acceleration(model, r_itrf))

                # Compare the results.
                @test g_norm ≈ expected_g_norm atol = 1e-8
            end
        end
    end

    # == Test the Signature with DateTime ==================================================

    dt_J2000 = DateTime("2000-01-01T12:00:00.000")

    # We will fetch the EIGEN-6C model that has time dependent coefficients.
    eigen6c_file = fetch_icgem_file(
        "http://icgem.gfz-potsdam.de/getmodel/gfc/0776caed6c65af24051697a65147b59e436cb464cb0930c1863fee6ecfbc31b0/EIGEN-6C.gfc"
    )

    eigen6c = GravityModels.load(IcgemFile, eigen6c_file)
    r_itrf = [7000.0e3, 0, 0]

    time = Dates.value(DateTime("2023-06-19") - dt_J2000) / 1000
    g_itrf_expected = GravityModels.gravity_acceleration(eigen6c, r_itrf, time)
    g_itrf = GravityModels.gravity_acceleration(eigen6c, r_itrf, DateTime("2023-06-19"))

    @test g_itrf == g_itrf_expected
end

# == File: ./src/GravityModels/gravitational_field_derivative.jl ===========================

# NOTE: We can skip most of the gravitational field derivative tests because it is
# indirectly tested in the previous sections.

@testset "Gravity Field Derivative" verbose = true begin
    dt_J2000 = DateTime("2000-01-01T12:00:00.000")

    # We will fetch the EIGEN-6C model that has time dependent coefficients.
    eigen6c_file = fetch_icgem_file(
        "http://icgem.gfz-potsdam.de/getmodel/gfc/0776caed6c65af24051697a65147b59e436cb464cb0930c1863fee6ecfbc31b0/EIGEN-6C.gfc"
    )

    eigen6c = GravityModels.load(IcgemFile, eigen6c_file)
    r_itrf = [7000.0e3, 0, 0]

    time = Dates.value(DateTime("2023-06-19") - dt_J2000) / 1000
    ∂U_∂r_expected, ∂U_∂ϕ_expected, ∂U_∂λ_expected =
        GravityModels.gravitational_field_derivative(eigen6c, r_itrf, time)

    ∂U_∂r, ∂U_∂ϕ, ∂U_∂λ =
        GravityModels.gravitational_field_derivative(eigen6c, r_itrf, DateTime("2023-06-19"))

    @test ∂U_∂r == ∂U_∂r_expected
    @test ∂U_∂ϕ == ∂U_∂ϕ_expected
    @test ∂U_∂λ == ∂U_∂λ_expected

    ∂U_∂r_expected, ∂U_∂ϕ_expected, ∂U_∂λ_expected =
        GravityModels.gravitational_field_derivative(eigen6c, r_itrf, 0)

    ∂U_∂r, ∂U_∂ϕ, ∂U_∂λ = GravityModels.gravitational_field_derivative(eigen6c, r_itrf)

    @test ∂U_∂r == ∂U_∂r_expected
    @test ∂U_∂ϕ == ∂U_∂ϕ_expected
    @test ∂U_∂λ == ∂U_∂λ_expected

    # Test when `max_degrees != max_order`.
    ∂U_∂r, ∂U_∂ϕ, ∂U_∂λ = GravityModels.gravitational_field_derivative(
        eigen6c,
        [6378.137e3, 0, 0];
        max_degree = 1,
        max_order = 0
    )

    @test ∂U_∂r ≈ -9.798285471812783
    @test ∂U_∂ϕ ≈ 0.0
    @test ∂U_∂λ ≈ 0.0
end

@testset "Gravity Field Derivative [ERRORS]" verbose = true begin
    egm96 = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))
    r_itrf = [7000.0e3, 0, 0]

    P = zeros(10, 10)
    @test_throws ArgumentError GravityModels.gravity_acceleration(egm96, r_itrf; P = P)

    P  = zeros(361, 361)
    dP = zeros(10, 10)
    @test_throws ArgumentError GravityModels.gravity_acceleration(egm96, r_itrf; dP = dP)
    @test_throws ArgumentError GravityModels.gravity_acceleration(egm96, r_itrf; P = P, dP = dP)
end

# == File: ./src/GravityModels/potential.jl ====================================================

@testset "Gravitational Potential" verbose = true begin
    dt_J2000 = DateTime("2000-01-01T12:00:00.000")

    # We will fetch the EIGEN-6C model that has time dependent coefficients.
    eigen6c_file = fetch_icgem_file(
        "http://icgem.gfz-potsdam.de/getmodel/gfc/0776caed6c65af24051697a65147b59e436cb464cb0930c1863fee6ecfbc31b0/EIGEN-6C.gfc"
    )

    eigen6c = GravityModels.load(IcgemFile, eigen6c_file)
    r_itrf = [7000.0e3, 0, 0]

    # == Test the Signature with DateTime ==================================================

    time = Dates.value(DateTime("2023-06-19") - dt_J2000) / 1000
    U_expected = GravityModels.gravitational_potential(eigen6c, r_itrf, time)
    U = GravityModels.gravitational_potential(eigen6c, r_itrf, DateTime("2023-06-19"))

    @test U == U_expected

    # == Test without time argument ========================================================

    U_expected = GravityModels.gravitational_potential(eigen6c, r_itrf, 0)
    U = GravityModels.gravitational_potential(eigen6c, r_itrf)

    @test U == U_expected

    # == Test with max_degree and max_order ================================================

    U = GravityModels.gravitational_potential(
        eigen6c,
        [6378.137e3, 0, 0];
        max_degree = 1,
        max_order = 0
    )

    # Verify the potential value is reasonable (should be positive and close to GM/r)
    μ = GravityModels.gravity_constant(eigen6c)
    r = 6378.137e3
    @test U > 0
    @test abs(U - (μ / r)) < 1e6  # Within 1 MJ/kg of the point mass approximation
end

@testset "Gravitational Potential Verification" verbose = true begin
    tests = (
        (:EGM96, "./test_results/potential/EGM96.gdf"),
        (:JGM2,  "./test_results/potential/JGM2.gdf"),
        (:JGM3,  "./test_results/potential/JGM3.gdf")
    )

    for t in tests
        @eval @testset $(string(t[1])) begin
            test_results = readdlm($t[2]; skipstart = 35)
            model = GravityModels.load(IcgemFile, fetch_icgem_file($t[1]))

            # Compare the results.
            for k in axes(test_results, 1)
                # Get latitude and longitude.
                lon = test_results[k, 0 + begin] |> deg2rad
                lat = test_results[k, 1 + begin] |> deg2rad

                # Get the expected result in m²/s².
                expected_U = test_results[k, 2 + begin]

                # Use the model to compute the gravitational potential using all coefficients.
                r_itrf = geodetic_to_ecef(lat, lon, 0)
                U = GravityModels.gravitational_potential(model, r_itrf)

                # Compare the results.
                @test U ≈ expected_U rtol = 1e-8
            end
        end
    end
end

@testset "Potential [ERRORS]" verbose = true begin
    egm96 = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))
    r_itrf = [7000.0e3, 0, 0]

    # Test that P matrix is too small
    P = zeros(10, 10)
    @test_throws ArgumentError GravityModels.gravitational_potential(egm96, r_itrf; P = P)
end
