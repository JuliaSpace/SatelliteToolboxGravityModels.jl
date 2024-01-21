## Description #############################################################################
#
# Tests related to the GravityModel API.
#
############################################################################################

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
            for k in 1:size(test_results, 1)
                # Get latitude and longitude.
                lat = test_results[k, 2] |> deg2rad
                lon = test_results[k, 1] |> deg2rad

                # Get the expected result in m / s² (the online result is in mGal).
                expected_g_norm = test_results[k, 3] / 100000

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
            for k in 1:size(test_results, 1)
                # Get latitude and longitude.
                lat = test_results[k, 2] |> deg2rad
                lon = test_results[k, 1] |> deg2rad

                # Get the expected result in m / s² (the online result is in mGal).
                expected_g_norm = test_results[k, 3] / 100000

                # Use the model to compute the gravity using all the coefficients.
                r_itrf = geodetic_to_ecef(lat, lon, 0)
                g_norm = norm(GravityModels.gravity_acceleration(model, r_itrf))

                # Compare the results.
                @test g_norm ≈ expected_g_norm atol = 1e-8
            end
        end
    end
end

# == File: ./src/GravityModels/gravitational_field_derivative.jl ===========================

# NOTE: We can skip the gravitational field derivative tests because it is indirectly tested
# in the previous sections.

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
