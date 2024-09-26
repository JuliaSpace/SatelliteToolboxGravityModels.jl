## Description #############################################################################
#
# Tests related to automatic differentiation.
#
############################################################################################

const _GRAVITY_MODELS = (
    (:EGM96,),
    (:JGM2, ),
    (:JGM3, )
)

const _BACKENDS = (
    ("ForwardDiff", AutoForwardDiff()),
)

@testset "State Automatic Differentiation" begin

    for t in _GRAVITY_MODELS
        model = GravityModels.load(IcgemFile, fetch_icgem_file(t[1]))
        lat   = 27.5
        lon   = 235.3
        h     = 0.0

        # Use the model to compute the gravity using all the coefficients.
        r_itrf = geodetic_to_ecef(lat, lon, 0)

        f_fd, df_fd = value_and_jacobian(
            (x) -> GravityModels.gravitational_acceleration(model, x),
            AutoFiniteDiff(),
            r_itrf
        )

        for backend in _BACKENDS
            @eval @testset $(string(t[1]) * " " * string(backend[1])) begin
                f_ad, df_ad = value_and_jacobian(
                    (x) -> GravityModels.gravitational_acceleration($model, x),
                    $backend[2],
                    $r_itrf
                )

                @test $f_fd == f_ad
                @test $df_fd ≈ df_ad rtol=1e-4
            end
        end
    end
end

@testset "Time Automatic Differentiation" begin
    for t in _GRAVITY_MODELS
        model = GravityModels.load(IcgemFile, fetch_icgem_file(t[1]))
        lat   = 27.5
        lon   = 235.3
        h     = 0.0

        # Use the model to compute the gravity using all the coefficients.
        r_itrf = geodetic_to_ecef(lat, lon, 0)
        time = 0.0

        f_fd, df_fd = value_and_derivative(
            (x) -> GravityModels.gravitational_acceleration(model, r_itrf, x[1]),
            AutoFiniteDiff(),
            time
        )

        for backend in _BACKENDS
            @eval @testset $(string(t[1]) * " " * string(backend[1])) begin
                f_ad, df_ad = value_and_derivative(
                    (x) -> GravityModels.gravitational_acceleration($model, $r_itrf, x[1]),
                    $backend[2],
                    $time
                )

                @test $f_fd == f_ad
                @test $df_fd ≈ df_ad rtol=1e-4
            end
        end
    end
end
