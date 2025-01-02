SatelliteToolboxGravityModels.jl Changelog
==========================================

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

[badge-breaking]: https://img.shields.io/badge/BREAKING-red.svg
[badge-deprecation]: https://img.shields.io/badge/Deprecation-orange.svg
[badge-feature]: https://img.shields.io/badge/Feature-green.svg
[badge-enhancement]: https://img.shields.io/badge/Enhancement-blue.svg
[badge-bugfix]: https://img.shields.io/badge/Bugfix-purple.svg
[badge-info]: https://img.shields.io/badge/Info-gray.svg

[gh-pr-3]: https://github.com/JuliaSpace/SatelliteToolboxGravityModels.jl/pull/3

