## Description #############################################################################
#
# Tests for the SatelliteToolboxGravityModelsZygoteExt extension, verifying that Zygote
# (via custom rrules backed by ForwardDiff) produces the same results as ForwardDiff alone.
#
############################################################################################

using DifferentiationInterface
using ForwardDiff, Zygote

const _GRAV_MODEL = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

lat    = 27.5
lon    = 235.3
r_itrf = Array(geodetic_to_ecef(lat, lon, 0))
time   = 0.0

@testset "gravitational_acceleration(model, r, t)" begin
    fn = (x) -> Array(GravityModels.gravitational_acceleration(_GRAV_MODEL, x[1:3], x[4]))
    input = [r_itrf; time]

    _, jac_fd = value_and_jacobian(fn, AutoForwardDiff(), input)
    _, jac_zy = value_and_jacobian(fn, AutoZygote(), input)

    @test jac_fd ≈ jac_zy rtol = 1e-12
end

@testset "gravitational_potential(model, r, t)" begin
    fn = (x) -> GravityModels.gravitational_potential(_GRAV_MODEL, x[1:3], x[4])
    input = [r_itrf; time]

    _, grad_fd = value_and_gradient(fn, AutoForwardDiff(), input)
    _, grad_zy = value_and_gradient(fn, AutoZygote(), input)

    @test grad_fd ≈ grad_zy rtol = 1e-12
end

@testset "gravitational_field_derivative(model, r, t)" begin
    fn = (x) -> collect(GravityModels.gravitational_field_derivative(_GRAV_MODEL, x[1:3], x[4]))
    input = [r_itrf; time]

    _, jac_fd = value_and_jacobian(fn, AutoForwardDiff(), input)
    _, jac_zy = value_and_jacobian(fn, AutoZygote(), input)

    @test jac_fd ≈ jac_zy rtol = 1e-12
end
