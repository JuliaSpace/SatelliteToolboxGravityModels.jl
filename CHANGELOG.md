SatelliteToolboxGravityModels.jl Changelog
==========================================

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

