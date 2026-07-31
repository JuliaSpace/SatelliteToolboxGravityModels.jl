SatelliteToolboxGravityModels.jl Changelog
==========================================

Version 1.4.0
-------------

- ![Enhancement][badge-enhancement] The pre-configured ICGEM model URLs now use https.
- ![Enhancement][badge-enhancement] The functions `gravitational_field_derivative`,
  `gravitational_potential`, `gravitational_acceleration`, and `gravity_acceleration` now
  compute the spherical harmonics with concrete types when the user does not provide the
  matrices `P` and `dP`, improving the performance in this case.
- ![Enhancement][badge-enhancement] The data type passed to `GravityModels.load` and
  `parse_icgem` is now inferable, and the ICGEM coefficient computation is type stable for
  all combinations of model and time types.
- ![Bugfix][badge-bugfix] The function `icgem_coefficients` no longer throws an
  `UndefVarError` when called with a `DateTime` object.
- ![Bugfix][badge-bugfix] The regex used to parse numbers in FORTRAN format no longer
  replaces commas and spaces in the input.
- ![Bugfix][badge-bugfix] The parser no longer throws a `MethodError` when reading files
  with time-variable coefficients using a data type other than `Float64`.
- ![Bugfix][badge-bugfix] The parser no longer enters an infinite loop when an invalid
  line follows a `gfct` section.
- ![Bugfix][badge-bugfix] The last coefficient of an ICGEM file is no longer lost when the
  file ends inside a `gfct` section.
- ![Bugfix][badge-bugfix] Printing an `IcgemGfctCoefficient` whose epoch has fractional
  seconds no longer throws an `InexactError`.
- ![Bugfix][badge-bugfix] The function `gravity_acceleration` no longer returns `NaN` at
  the poles due to the centrifugal acceleration term, which is now computed using a
  simplified and faster expression.
- ![Bugfix][badge-bugfix] The elapsed time used to compute time-variable coefficients is
  now converted to years using the Julian year (365.25 days) instead of 365-day years.
  This modification slightly changes the results of models with time-variable
  coefficients.
- ![Bugfix][badge-bugfix] The function `fetch_icgem_file` now downloads the file
  atomically. Hence, an interrupted download is no longer treated as a valid cached file.
- ![Info][badge-info] The documentation of all functions, types, and structures was
  reviewed and improved, including many typo fixes.

Version 1.3.0
-------------

- ![Feature][badge-feature] The package can now compute the gravitational potential. This
  feature is also differentiable. (PR [#8][gh-pr-8])
- ![Feature][badge-feature] The package now supports ICGEM files for any body, not only
  Earth. (PR [#9][gh-pr-9])

Version 1.2.0
-------------

- ![Enhancement][badge-enhancement] The ICGEM structure and parsing algorithm was updated to
  avoid type-instabilities when storing the coefficients. This required to change the
  structure signature. However, this modification is internal to the package. Hence, this is
  not a breaking change. (PR [#6][gh-pr-6])
- ![Enhancement][badge-enhancement] If the user does not provide the storage matrices `P`
  and `dP`, the algorithm now uses the `LowerTriangularStorage` structure from the package
  `SatelliteToolbox.jl`. This modification provided a huge gain when evaluating large models
  such as the EGM2008, which decreased the time to compute the gravity acceleration by 40%.

Version 1.1.0
-------------

- ![Feature][badge-feature] The package now supports automatic differentiation using
  different backends. (PR [#4][gh-pr-4])
- ![Enhancement][badge-enhancement] Some allocations were removed. (PR [#4][gh-pr-4])

Version 1.0.0
-------------

- ![Info][badge-info] We dropped support for Julia 1.6. This version only supports the
  current Julia version and v1.10 (LTS).
- ![Info][badge-info] This version does not have breaking changes. We bump the version to
  1.0.0 because we now consider the API stable.

Version 0.1.6
-------------

- ![Enhancement][badge-enhancement] The pacakge is now compatible with auto-differentiation
  tools. (PR [#3][gh-pr-3])
- ![Bugfix][badge-bugfix] In the previous version, we documented that calling functions
  without the time information will use the J2000.0 epoch. However, we were using
  `2000-01-01T00:00:00` instead of `2000-01-01T12:00:00`. This bug was fixed and now we use
  the correct J2000.0 epoch.
- ![Info][badge-info] Due to external packages, we cannot test
  SatelliteToolboxGravityModels.jl against Julia 1.6 anymore. The support for this version
  will be removed in a future release.

Version 0.1.5
-------------

- ![Enhancement][badge-enhancement] Minor source-code updates.
- ![Enhancement][badge-enhancement] Documentation updates.

Version 0.1.4
-------------

- ![Enhancement][badge-enhancement] We now **truly** export `AbstractGravityModel`.

Version 0.1.3
-------------

- ![Enhancement][badge-enhancement] We exported `AbstractGravityModel`.
- ![Enhancement][badge-enhancement] The function `GravityModels.coefficients` can be called
  without the parameter `time`. In this case, J2000.0 epoch will be used.

Version 0.1.2
-------------

- ![Enhancement][badge-enhancement] We updated the dependency compatibility bounds.

Version 0.1.1
-------------

- ![Bugfix][badge-bugfix] In the previous version, we were accessing undefined memory
  regions when compute the gravitational field derivative if the maximum order is lower than
  maximum degree. We always need to compute `P` with one order higher than `dP` in those
  cases.

Version 0.1.0
-------------

- Initial version.
  - This version was based on the code in **SatelliteToolbox.jl**.

[badge-breaking]: https://img.shields.io/badge/Breaking-DC2626?style=flat-square
[badge-deprecation]: https://img.shields.io/badge/Deprecation-D97706?style=flat-square
[badge-feature]: https://img.shields.io/badge/Feature-16A34A?style=flat-square
[badge-enhancement]: https://img.shields.io/badge/Enhancement-0284C7?style=flat-square
[badge-bugfix]: https://img.shields.io/badge/Bugfix-DB2777?style=flat-square
[badge-info]: https://img.shields.io/badge/Info-475569?style=flat-square

[gh-pr-3]: https://github.com/JuliaSpace/SatelliteToolboxGravityModels.jl/pull/3
[gh-pr-4]: https://github.com/JuliaSpace/SatelliteToolboxGravityModels.jl/pull/4
[gh-pr-6]: https://github.com/JuliaSpace/SatelliteToolboxGravityModels.jl/pull/6
[gh-pr-8]: https://github.com/JuliaSpace/SatelliteToolboxGravityModels.jl/pull/8
[gh-pr-9]: https://github.com/JuliaSpace/SatelliteToolboxGravityModels.jl/pull/9

