SatelliteToolboxGravityModels.jl Changelog
==========================================

Version 2.0.0
-------------

- ![BREAKING][badge-breaking] The keywords `P` and `dP` of `gravitational_potential`,
  `gravitational_field_derivative`, `gravitational_acceleration`, and
  `gravity_acceleration` were replaced by `workspace`, which receives a
  `GravityModels.Workspace`. The workspace, created with `GravityModels.Workspace(model)`,
  holds the buffers of the associated Legendre functions and the precomputed recursion
  coefficients introduced in SatelliteToolboxLegendre.jl v1.2. The functions throw an
  `ArgumentError` if the workspace element type or its maximum degree and order do not
  match the computation.
- ![BREAKING][badge-breaking] `AbstractGravityModel` lost the norm type parameter and is
  now `AbstractGravityModel{T}`. The API function `GravityModels.coefficient_norm` must
  return the normalization wrapped in a `Val` (`Val(:full)`, `Val(:schmidt)`, or
  `Val(:unnormalized)`), which must be inferable from the model type.
- ![BREAKING][badge-breaking] The API function `GravityModels.angular_speed(model)` was
  added and must be implemented by the models. It returns the angular speed [rad/s] of
  the central body and is the default of the keyword `ω` of `gravity_acceleration`, which
  previously defaulted to Earth's value for any model.
- ![BREAKING][badge-breaking] `IcgemFile` is now `IcgemFile{T, N}`, where `N` is the `Val`
  with the normalization of the coefficients. The time-variable coefficients are stored
  sparsely in the new types `IcgemTimeVariableCoefficient` and `IcgemPeriodicTerm`, while
  the constant coefficients are stored for every degree and order. The types
  `AbstractIcgemCoefficient` and `IcgemGfctCoefficient` were removed.
- ![BREAKING][badge-breaking] The parser throws the new exception `IcgemParseError`, which
  carries the line number when applicable, instead of an `ErrorException`.
- ![BREAKING][badge-breaking] The rich representations of the ICGEM types follow the tree
  layout of SatelliteToolboxBase.jl v2.1, whose printing helpers are now used. Crayons.jl
  and ReferenceFrameRotations.jl are no longer dependencies.
- ![Feature][badge-feature] The ICGEM format 2.0 is supported. The parser reads the
  validity interval of the time-variable coefficients, whose epochs can be written as
  `yyyymmdd.hhmm`, stores one object per interval, and the evaluation selects the interval
  that contains the requested time, clamping to the first or last one outside the covered
  period. Header values followed by comments, such as
  `errors formal (sigma calibration factor = 1.00)`, are accepted.
- ![Feature][badge-feature] `IcgemFile` stores the angular speed of the central body,
  which can be set with the keyword `angular_speed` of `GravityModels.load` and
  `parse_icgem`, defaulting to Earth's value.
- ![Feature][badge-feature] `parse_icgem` accepts an `IO` stream.
- ![Enhancement][badge-enhancement] The gravitational acceleration of EGM96 is computed
  1.35 times faster with a workspace and, in this case, the evaluation functions do not
  allocate.
- ![Enhancement][badge-enhancement] The evaluation of models with time-variable
  coefficients, such as GOCO06s, is 3 times faster and uses 35% less memory, and the sine
  and cosine periodic terms with the same period share the trigonometric evaluations.
- ![Enhancement][badge-enhancement] The conversion of the time argument is centralized,
  and the methods without the time argument were merged into the methods with a default
  value.
- ![Enhancement][badge-enhancement] The parser skips data lines whose degree or order are
  out of range, or whose epoch is invalid, logging a warning instead of throwing.
- ![Bugfix][badge-bugfix] The acceleration at positions whose latitude rounds to ±π / 2,
  such as those obtained from `geodetic_to_ecef(±π / 2, λ, h)`, had the north component
  with the wrong sign, an error of 1.2e-4 m/s², and the east component on the polar axis
  was missing. The Legendre functions are now evaluated at the angle from the polar axis
  computed directly from the coordinates and folded to the northern hemisphere, and the
  east component on the axis is obtained from its limit.
- ![Bugfix][badge-bugfix] The coefficients omitted in an ICGEM file, such as those of
  degree 1, were uninitialized memory instead of 0.
- ![Bugfix][badge-bugfix] The ICGEM file handle was never closed after parsing.
- ![Bugfix][badge-bugfix] Differentiating with Zygote.jl a call that received the user
  buffers failed because they were forwarded to the ForwardDiff.jl pullbacks.
- ![Bugfix][badge-bugfix] The value π / 2 used in the kernels is now evaluated in the
  result type, improving the precision for types wider than `Float64`.
- ![Info][badge-info] Julia 1.13 was added to the supported versions.
- ![Info][badge-info] The documentation and the docstrings were reviewed for the new API.

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

