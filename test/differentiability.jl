## Description #############################################################################
#
# Tests related to automatic differentiation.
#
############################################################################################

const _GRAV_MODEL = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

@testset "State Automatic Differentiation" begin
    for backend in _BACKENDS
        testset_name = "Gravity Models " * string(backend[1])
        @testset "$testset_name" begin
            lat   = 27.5
            lon   = 235.3
            h     = 0.0

            # Use the model to compute the gravity using all the coefficients.
            r_itrf = geodetic_to_ecef(lat, lon, 0)

            f_fd, df_fd = value_and_jacobian(
                (x) -> GravityModels.gravitational_acceleration(_GRAV_MODEL, x),
                AutoFiniteDiff(),
                Array(r_itrf)
            )
        
            f_ad, df_ad = value_and_jacobian(
                (x) -> Array(GravityModels.gravitational_acceleration(_GRAV_MODEL, x)),
                backend[2],
                Array(r_itrf)
            )

            @test f_fd ≈ f_ad rtol=1e-14
            @test df_fd ≈ df_ad rtol=1e-4

            f_fd2, df_fd2 = value_and_gradient(
                (x) -> GravityModels.gravitational_potential(_GRAV_MODEL, x),
                AutoFiniteDiff(),
                Array(r_itrf)
            )
        
            f_ad2, df_ad2 = value_and_gradient(
                (x) -> GravityModels.gravitational_potential(_GRAV_MODEL, x),
                backend[2],
                Array(r_itrf)
            )

            @test f_fd2 ≈ f_ad2 rtol=1e-14
            @test df_fd2 ≈ df_ad2 rtol=1e-4

            f_fd3, df_fd3 = value_and_jacobian(
                (x) -> collect(GravityModels.gravitational_field_derivative(_GRAV_MODEL, x)),
                AutoFiniteDiff(),
                Array(r_itrf)
            )
        
            f_ad3, df_ad3 = value_and_jacobian(
                (x) -> collect(GravityModels.gravitational_field_derivative(_GRAV_MODEL, x)),
                backend[2],
                Array(r_itrf)
            )

            @test f_fd3 ≈ f_ad3 rtol=1e-14
            @test df_fd3 ≈ df_ad3 rtol=1e-4
        end
    end
end

@testset "Time Automatic Differentiation" begin
    for backend in _BACKENDS
        testset_name = "Gravity Models " * string(backend[1])
        @testset "$testset_name" begin
            lat   = 27.5
            lon   = 235.3
            h     = 0.0

            # Use the model to compute the gravity using all the coefficients.
            r_itrf = geodetic_to_ecef(lat, lon, 0)
            time = 0.0

            f_fd, df_fd = value_and_derivative(
                (x) -> Array(GravityModels.gravitational_acceleration(_GRAV_MODEL, r_itrf, x)),
                AutoFiniteDiff(),
                time
            )

            f_ad, df_ad = value_and_derivative(
                (x) -> Array(GravityModels.gravitational_acceleration(_GRAV_MODEL, r_itrf, x)),
                backend[2],
                time
            )

            @test f_fd ≈ f_ad rtol=1e-14
            @test df_fd ≈ df_ad rtol=1e-4

            f_fd2, df_fd2 = value_and_derivative(
                (x) -> GravityModels.gravitational_potential(_GRAV_MODEL, r_itrf, x),
                AutoFiniteDiff(),
                time
            )

            f_ad2, df_ad2 = value_and_derivative(
                (x) -> GravityModels.gravitational_potential(_GRAV_MODEL, r_itrf, x),
                backend[2],
                time
            )

            @test f_fd2 ≈ f_ad2 rtol=1e-14
            @test df_fd2 ≈ df_ad2 rtol=1e-4

            f_fd3, df_fd3 = value_and_derivative(
                (x) -> collect(GravityModels.gravitational_field_derivative(_GRAV_MODEL, r_itrf, x)),
                AutoFiniteDiff(),
                time
            )

            f_ad3, df_ad3 = value_and_derivative(
                (x) -> collect(GravityModels.gravitational_field_derivative(_GRAV_MODEL, r_itrf, x)),
                backend[2],
                time
            )

            @test f_fd3 ≈ f_ad3 rtol=1e-14
            @test df_fd3 ≈ df_ad3 rtol=1e-4
        end
    end
end

@testset "Rotation Rate Automatic Differentiation" begin
    for backend in _BACKENDS
        testset_name = "Gravity Models " * string(backend[1])
        @testset "$testset_name" begin
            lat   = 27.5
            lon   = 235.3
            h     = 0.0

            # Use the model to compute the gravity using all the coefficients.
            r_itrf = geodetic_to_ecef(lat, lon, 0)
            time = 0.0

            f_fd, df_fd = value_and_derivative(
                (x) -> Array(GravityModels.gravity_acceleration(_GRAV_MODEL, r_itrf; ω = x)),
                AutoFiniteDiff(),
                EARTH_ANGULAR_SPEED
            )

            f_ad, df_ad = value_and_derivative(
                (x) -> Array(GravityModels.gravity_acceleration(_GRAV_MODEL, r_itrf; ω = x)),
                backend[2],
                EARTH_ANGULAR_SPEED
            )

            f_ad2, df_ad2 = value_and_derivative(
                (x) -> Array(GravityModels.gravity_acceleration(_GRAV_MODEL, r_itrf, time; ω = x)),
                backend[2],
                EARTH_ANGULAR_SPEED
            )

            @test f_fd ≈ f_ad rtol=1e-14
            @test df_fd ≈ df_ad rtol=1e-2
            @test f_fd ≈ f_ad2 rtol=1e-14
            @test df_fd ≈ df_ad2 rtol=1e-2
        end
    end
end
