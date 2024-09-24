@testset "Automatic Differentiation Tests" begin

    tests = (
        (:EGM96, ),
        (:JGM2,  ),
        (:JGM3,  )
    )

    for t in tests
        @eval @testset $(string(t[1])) begin
            model = GravityModels.load(IcgemFile, fetch_icgem_file($t[1]))

            lat = 27.5
            lon = 235.3
            h = 0.0

            # Use the model to compute the gravity using all the coefficients.
            r_itrf = geodetic_to_ecef(lat, lon, 0)

            g = ForwardDiff.jacobian((x) -> GravityModels.gravitational_acceleration(model, x), r_itrf)
            #g2 = Zygote.jacobian((x) -> GravityModels.gravitational_acceleration(model, x), r_itrf)

            print(g)
        end
    end
end