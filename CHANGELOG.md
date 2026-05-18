# Changelog

All notable changes to this project will be documented in this file.

## [0.14.0] - 2025-04-13

### Added

- New `ColorCoordinates` object. Default color representation, carries color
  space information.
- Colormap generation.
- Add HWB color space.
- Add CIE1960UCS color space.
- Add CIEUVW color space.

### Fixed

- Be more precise in color difference computation.

### Changes

- Drop support for Python 3.8

## [0.13.0] - 2023-03-08

### Added

- HCT color space.
- OKLCH color space.
- ZCAM color appearance model.
- I<sub>G</sub>P<sub>G</sub>T<sub>G</sub> color space.

### Fixed

- HSL/HSV are now treated as full-fledged color spaces.
- Division-by-zero fix in LCH
- Various small fixes throughout

### Changes

- Color appearance models are now grouped under (`coloria.cam`)
- Color appearance models' `from_xyz100()` now return named tuples. You can
  access their values by, e.g.,

  ```python
  import coloria

  ciecam02 = coloria.cam.CIECAM02("average", 20, 100)
  corr = ciecam02.from_xyz100([19.31, 23.93, 10.14])
  corr.lightness
  # corr.lightness, corr.brightness, corr.chroma, ...
  ```

- Color appearance models now initialize with a literal surround parameters, e.g.,
  ```python
  import coloria
  ciecam02 = coloria.cam.CIECAM02("average", 20, 100)
  ```
  Possible values are `"average"`, `"dim"`, and `"dark"`.
