# colorio

Tools for color models.

[![CircleCI](https://img.shields.io/circleci/project/github/nschloe/colorio/master.svg)](https://circleci.com/gh/nschloe/colorio/tree/master)
[![codecov](https://img.shields.io/codecov/c/github/nschloe/colorio.svg)](https://codecov.io/gh/nschloe/colorio)
[![Codacy grade](https://img.shields.io/codacy/grade/b23fbc2af9884315bd7d6275aa2629b6.svg)](https://app.codacy.com/app/nschloe/colorio/dashboard)
[![awesome](https://img.shields.io/badge/awesome-yes-brightgreen.svg)](https://github.com/nschloe/colorio)
[![PyPi Version](https://img.shields.io/pypi/v/colorio.svg)](https://pypi.python.org/pypi/colorio)
[![GitHub stars](https://img.shields.io/github/stars/nschloe/colorio.svg?style=social&label=Stars)](https://github.com/nschloe/colorio)

All methods in colorio are fully vectorized.

### Color spaces

All color spaces implement the two methods
```python
vals = colorspace.from_xyz100(xyz)
xyz = colorspace.to_xyz100(vals)
```
for conversion from and to XYZ100. Adding new color maps is as easy as that.

The images below all show the SRGB gamut in the respective color space.

#### SRGB

<img src="https://nschloe.github.io/colorio/srgb.png" width="40%">

Linear SRGB space.

This class has the two additional methods
```
from_srgb1()
to_srgb1()
```
for conversion from and to standard RGB.

#### CIE XYZ
<img src="https://nschloe.github.io/colorio/xyz.png" width="40%">

```python
import colorio
import numpy

colorspace = colorio.XYZ()
xyz100 = numpy.random.rand(3)
vals = colorspace.from_xyz100(xyz100)
xyz100 = colorspace.to_xyz100(vals)
```

#### CIE XYY
<img src="https://nschloe.github.io/colorio/xyy.png" width="40%">

```python
import colorio
import numpy

colorspace = colorio.XYY()
xyz100 = numpy.random.rand(3)
vals = colorspace.from_xyz100(xyz100)
xyz100 = colorspace.to_xyz100(vals)
```

#### CIELAB
<img src="https://nschloe.github.io/colorio/cielab.png" width="40%">

```python
import colorio
import numpy

colorspace = colorio.CIELAB()
xyz100 = numpy.random.rand(3)
vals = colorspace.from_xyz100(xyz100)
xyz100 = colorspace.to_xyz100(vals)
```

#### CIELUV
<img src="https://nschloe.github.io/colorio/cieluv.png" width="40%">

```python
import colorio
import numpy

colorspace = colorio.CIELUV()
xyz100 = numpy.random.rand(3)
vals = colorspace.from_xyz100(xyz100)
xyz100 = colorspace.to_xyz100(vals)
```

#### CIECAM / CAM02-UCS
<img src="https://nschloe.github.io/colorio/cam02ucs.png" width="40%">

```python
import colorio
import numpy

L_A = 64 / numpy.pi / 5
ciecam02 = colorio.CIECAM02(0.69, 20, L_A)

xyz100 = numpy.random.rand(3)
J, C, H, h, M, s, Q = ciecam02.from_xyz100(xyz100)
xyz100 = ciecam02.to_xyz100(numpy.array([J, C, H]), 'JCH')

cam02 = colorio.CAM02('UCS', 0.69, 20, L_A)
vals = cam02.from_xyz100(xyz100)
xyz100 = cam02.to_xyz100(vals)
```

The implementation contains a few improvements over the CIECAM02 specification.

#### CAM16 / CAM16-UCS
<img src="https://nschloe.github.io/colorio/cam16ucs.png" width="40%">

From the article [Comprehensive color solutions: CAM16, CAT16, and
CAM16-UCS](https://doi.org/10.1002/col.22131) by Li et al.

```python
import colorio
import numpy

L_A = 64 / numpy.pi / 5
cs = colorio.CAM16(0.69, 20, L_A)

xyz100 = numpy.random.rand(3)
J, C, H, h, M, s, Q = cs.from_xyz100(xyz100)
xyz100 = ciecam02.to_xyz100(numpy.array([J, C, H]), 'JCH')

cs = colorio.CAM16UCS(0.69, 20, L_A)
vals = cs.from_xyz100(xyz100)
xyz100 = cs.to_xyz100(vals)
```

The implementation contains a few improvements over the CAM16 specification.

#### J<sub>z</sub>a<sub>z</sub>b<sub>z</sub>
<img src="https://nschloe.github.io/colorio/jzazbz.png" width="40%">

From the article [Perceptually uniform color space for image signals including
high dynamic range and wide gamut](https://doi.org/10.1364/OE.25.015131) by
Safdar et al.

```python
import colorio
import numpy

colorspace = colorio.JzAzBz()
xyz100 = numpy.random.rand(3)
vals = colorspace.from_xyz100(xyz100)
xyz100 = colorspace.to_xyz100(vals)
```

### Other tools

* To create the above gamut plots, write the corresponding mesh out to a file with`
  ```python
  colorio.show_srgb_gamut(colorspace, 'srgb.vtu', n=20)
  ```
  and open it with any program that knows how to handle VTU (e.g.,
  [ParaView](https://www.paraview.org/)). The `srgb` data set is defined in all
  points of the mesh.

* Color differences in any space can be computed with `colorio.delta(a, b)`.

### Installation

colorio is [available from the Python Package Index](https://pypi.python.org/pypi/colorio/), so with
```
pip install -U colorio
```
you can install/upgrade.

### Testing

To run the tests, simply check out this repository and run
```
pytest
```

### Distribution

To create a new release

1. bump the `__version__` number,

2. publish to PyPi and GitHub:
    ```
    $ make publish
    ```

### License
colorio is published under the [MIT license](https://en.wikipedia.org/wiki/MIT_License).
