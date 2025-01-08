## Description #############################################################################
#
# Tests related to automatic differentiation.
#
############################################################################################

const _GRAV_MODEL = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

const _BACKENDS = (
    ("ForwardDiff", AutoForwardDiff()),
    ("Mooncake", AutoMooncake(;config=nothing)),
    ("PolyesterForwardDiff", AutoPolyesterForwardDiff()),
    ("Zygote", AutoZygote()),
)

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
                r_itrf
            )
        
            f_ad, df_ad = value_and_jacobian(
                (x) -> Array(GravityModels.gravitational_acceleration(_GRAV_MODEL, x)),
                backend[2],
                Array(r_itrf)
            )

            @test f_fd ≈ f_ad rtol=1e-14
            @test df_fd ≈ df_ad rtol=1e-4
        end
    end

    testset_name = "Gravity Models Enzyme"
    @testset "$testset_name" begin
        lat   = 27.5
        lon   = 235.3
        h     = 0.0

        # Use the model to compute the gravity using all the coefficients.
        r_itrf = geodetic_to_ecef(lat, lon, 0)

        f_fd, df_fd = value_and_jacobian(
            (x) -> GravityModels.gravitational_acceleration(_GRAV_MODEL, x),
            AutoFiniteDiff(),
            r_itrf
        )
    
        f_ad, df_ad = value_and_jacobian(
            (x) -> Array(GravityModels.gravitational_acceleration(_GRAV_MODEL, x)),
            AutoEnzyme(),
            Array(r_itrf)
        )

        @test f_fd ≈ f_ad rtol=1e-14
        @test df_fd ≈ df_ad rtol=1e-4
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
                (x) -> GravityModels.gravitational_acceleration(_GRAV_MODEL, r_itrf, x),
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
        end
    end

    testset_name = "Gravity Models Enzyme"
    @testset "$testset_name" begin
        lat   = 27.5
        lon   = 235.3
        h     = 0.0

        # Use the model to compute the gravity using all the coefficients.
        r_itrf = geodetic_to_ecef(lat, lon, 0)
        time = 0.0

        f_fd, df_fd = value_and_derivative(
            (x) -> GravityModels.gravitational_acceleration(_GRAV_MODEL, r_itrf, x),
            AutoFiniteDiff(),
            time
        )

        f_ad, df_ad = value_and_derivative(
            Const((x) -> Array(GravityModels.gravitational_acceleration(_GRAV_MODEL, r_itrf, x))),
            AutoEnzyme(),
            time
        )

        @test f_fd ≈ f_ad rtol=1e-14
        @test df_fd ≈ df_ad rtol=1e-4
    end
end
