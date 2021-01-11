<p align="center">
  <a href="https://github.com/nschloe/colorio"><img alt="colorio" src="https://nschloe.github.io/colorio/colorio-logo.svg" width="60%"></a>
  <p align="center">Tools for color models.</p>
</p>

[![PyPi Version](https://img.shields.io/pypi/v/colorio.svg?style=flat-square)](https://pypi.org/project/colorio)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/colorio.svg?style=flat-square)](https://pypi.org/pypi/colorio/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1172995.svg?style=flat-square)](https://doi.org/10.5281/zenodo.1172995)
[![GitHub stars](https://img.shields.io/github/stars/nschloe/colorio.svg?style=flat-square&logo=github&label=Stars&logoColor=white)](https://github.com/nschloe/colorio)
[![PyPi downloads](https://img.shields.io/pypi/dm/colorio.svg?style=flat-square)](https://pypistats.org/packages/colorio)

[![gh-actions](https://img.shields.io/github/workflow/status/nschloe/colorio/ci?style=flat-square)](https://github.com/nschloe/colorio/actions?query=workflow%3Aci)
[![codecov](https://img.shields.io/codecov/c/github/nschloe/colorio.svg?style=flat-square)](https://codecov.io/gh/nschloe/colorio)
[![LGTM](https://img.shields.io/lgtm/grade/python/github/nschloe/colorio.svg?style=flat-square)](https://lgtm.com/projects/g/nschloe/colorio)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg?style=flat-square)](https://github.com/psf/black)

### Color spaces

All color spaces implement the two methods
<!--exdown-skip-->
```python
vals = colorspace.from_xyz100(xyz)
xyz = colorspace.to_xyz100(vals)
```
for conversion from and to XYZ100. Adding new color spaces is as easy as
writing a class that provides those two methods.

The following color spaces are implemented:

 * [XYZ](https://en.wikipedia.org/wiki/CIE_1931_color_space) (`colorio.cs.XYZ1()`)
 * [xyY](https://en.wikipedia.org/wiki/CIE_1931_color_space#CIE_xy_chromaticity_diagram_and_the_CIE_xyY_color_space) (`colorio.cs.XYY1()`)
 * [Linear SRGB](https://en.wikipedia.org/wiki/SRGB)  (`colorio.SrgbLinear()`)
   This class has the two additional methods
   ```
   from_rgb1()
   to_rgb1()
   ```
   for conversion from and to standard RGB.
 * [HSL and HSV](https://en.wikipedia.org/wiki/HSL_and_HSV)  (`colorio.cs.HSL()`,
   `colorio.cs.HSV()`)
   These classes have the two methods
   ```
   from_srgb1()
   to_srgb1()
   ```
   for conversion from and to standard RGB.
 * [OSA-UCS](https://en.wikipedia.org/wiki/OSA-UCS) (`colorio.cs.OsaUcs()`)
 * [CIELAB](https://en.wikipedia.org/wiki/Lab_color_space) (`colorio.cs.CIELAB()`)
 * [CIELUV](https://en.wikipedia.org/wiki/CIELUV) (`colorio.cs.CIELUV()`)
 * [RLAB](https://doi.org/10.1002/(SICI)1520-6378(199610)21:5<338::AID-COL3>3.0.CO;2-Z)
   (`colorio.cs.RLAB()`)
 * [ICtCp](https://en.wikipedia.org/wiki/ICtCp) (`colorio.cs.ICtCp()`)
 * [IPT](http://www.ingentaconnect.com/content/ist/cic/1998/00001998/00000001/art00003)
   (`colorio.cs.IPT()`)
 * [CIECAM02 / CAM02-UCS](https://en.wikipedia.org/wiki/CIECAM02)
   ```python
   import math
   import colorio

   L_A = 64 / math.pi / 5
   ciecam02 = colorio.cs.CIECAM02(0.69, 20, L_A)
   cam02 = colorio.cs.CAM02("UCS", 0.69, 20, L_A)
   ```
   The implementation contains a few improvements over the CIECAM02 specification (see
   [here](https://arxiv.org/abs/1802.06067)).
 * [CAM16 / CAM16-UCS](https://doi.org/10.1002/col.22131)
   ```python
   import math
   import colorio

   L_A = 64 / math.pi / 5
   cam16 = colorio.cs.CAM16(0.69, 20, L_A)
   cam16ucs = colorio.cs.CAM16UCS(0.69, 20, L_A)
   ```
   The implementation contains a few improvements over the CAM16
   specification (see [here](https://arxiv.org/abs/1802.06067)).
 * [J<sub>z</sub>a<sub>z</sub>b<sub>z</sub>](https://doi.org/10.1364/OE.25.015131)
   (`colorio.cs.JzAzBz()`)
 * [OKLAB](https://bottosson.github.io/posts/oklab/)(`colorio.cs.OKLAB()`)

All methods in colorio are fully vectorized, i.e., computation is _really_ fast.

### Tools

colorio provides a number of useful tools for analyzing and visualizing color spaces.

#### Visualizing the SRGB gamut

<img src="https://nschloe.github.io/colorio/srgb-gamut-xyz.png" width="100%"> | <img src="https://nschloe.github.io/colorio/srgb-gamut-cielab.png" width="100%"> | <img src="https://nschloe.github.io/colorio/srgb-gamut-cam16.png" width="100%">
:---:|:-------:|:----------:|
XYZ  |  CIELAB |  CAM16-UCS |

The SRGB gamut is a perfect cube in SRGB space, and takes curious shapes when translated
into other color spaces. The above image shows the SRGB gamut in XYZ space. The image
data was created with
```python
import colorio

colorspace = colorio.cs.CIELAB()
colorspace.save_rgb_gamut("srgb.vtk", "srgb", n=50)

# The HDR (Rec.2100, Rec.2020) gamut works the same way
colorspace.save_rgb_gamut("hdr.vtk", "hdr", n=50)
```
The [VTK](https://www.vtk.org/VTK/img/file-formats.pdf) file can then be opened
in, e.g., ParaView, where the following instructions apply:
1. Open the file in ParaView and execute the following steps in the **Properties** tab
   to the left.
2. Press the **Apply** button.
3. Under the **Coloring** section, change `Solid Color` to `srgb`.
4. If necessary, press the gear button to the right of the search field to activate
   advanced options.
5. Under the **Scalar Coloring** section, uncheck **Map Scalars**.

More images are [in the gh-pages
branch](https://github.com/nschloe/colorio/tree/gh-pages).

The data can be written in most formats supported by
[meshio](https://github.com/nschloe/meshio). (You might have to install additional
packages for some formats.)

#### Visualizing the visible gamut

<img src="https://nschloe.github.io/colorio/visible-gamut-xyz.png" width="100%"> | <img src="https://nschloe.github.io/colorio/visible-gamut-cielab.png" width="100%"> | <img src="https://nschloe.github.io/colorio/visible-gamut-cam16.png" width="100%">
:---:|:-------:|:----------:|
XYZ  |  CIELAB |  CAM16-UCS |

Same as above, but with the gamut visible under a given illuminant.
```python
import colorio

illuminant = colorio.illuminants.d65()
observer = colorio.observers.cie_1931_2()

colorspace = colorio.cs.XYZ1()
colorspace.save_visible_gamut(observer, illuminant, "visible.vtk")
```
The gamut is shown in grey since SRGB screens are not able to display the colors anyway.

#### Slices through the color spaces

Instead of fiddling around with the proper 3D objects, colorio can plot slices through
all color spaces. One simply provides the slice index (typically the one that
corresponds to "lightness" in the respective color space, e.g., 2 in xyY and 0 in
CIELAB) and the slice level.

The solid line corresponds to monochromatic light; for convenience, the slice through
the SRGB gamut is also displayed.

<img src="https://nschloe.github.io/colorio/xyy-visible-slice.png" width="100%"> | <img src="https://nschloe.github.io/colorio/cielab-visible-slice.png" width="100%"> | <img src="https://nschloe.github.io/colorio/cam16ucs-visible-slice.png" width="100%">
:--------------:|:-------------------:|:---------------------:|
xyY (at Y=0.4)  |  CIELAB (at L=50)  |  CAM16-UCS (at J'=50) |

```python
import colorio

# xyy = colorio.cs.XYY1()
# xyy.show_visible_slice("xyy-visible-slice.png", 2, 0.4)

# cielab = colorio.cs.CIELAB()
# cielab.show_visible_slice(0, 50)

cam16 = colorio.cs.CAM16UCS(0.69, 20, 4.07)
cam16.show_visible_slice(0, 50)
# cam16.save_visible_slice("cam16ucs-visible-slice.png", 0, 50)
```

For convenience, it is also possible to show the classical visible gamut in xy with
[Planckian locus](https://en.wikipedia.org/wiki/Planckian_locus) and the SRGB colors (at
maximum luminosity).

<img src="https://nschloe.github.io/colorio/xy-gamut.png" width="30%">

```python
import colorio

colorio.show_flat_gamut()
```


#### Show experimental data

colorio contains lots of experimental data sets some of which can be used to assess
certain properties of color spaces.


###### MacAdam ellipses

<img src="https://nschloe.github.io/colorio/macadam-xyy.png" width="100%"> | <img src="https://nschloe.github.io/colorio/macadam-cielab.png" width="100%"> | <img src="https://nschloe.github.io/colorio/macadam-cam16.png" width="100%">
:--------------:|:------------------:|:----------------:|
xyY (at Y=0.4)  |  CIELAB (at L=50)  |  CAM16 (at L=50) |

The famous MacAdam ellipses (from [this
article](https://doi.org/10.1364%2FJOSA.32.000247)) can be plotted with

```python
import colorio

# xyy = colorio.cs.XYY1()
# colorio.data.macadam_1942.show(xyy, 0.4)
# colorio.data.macadam_1942.savefig("macadam-xyy.png", xyy, 0.4)

cieluv = colorio.cs.CIELUV()
colorio.data.macadam_1942.show(cieluv, 50.0)
```

###### Luo-Rigg

<img src="https://nschloe.github.io/colorio/luo-rigg-xyy.png" width="100%"> | <img src="https://nschloe.github.io/colorio/luo-rigg-cielab.png" width="100%"> | <img src="https://nschloe.github.io/colorio/luo-rigg-cam16.png" width="100%">
:--------------:|:------------------:|:----------------:|
xyY (at Y=0.4)  |  CIELAB (at L=50)  |  CAM16 (at L=50) |

Likewise for [Luo-Rigg](https://doi.org/10.1002/col.5080110107).

```python
import colorio

# xyy = colorio.cs.XYY1()
# colorio.data.luo_rigg.show(xyy, 0.4)
# colorio.data.luo_rigg.savefig("luo-rigg-xyy.png", xyy, 0.4)

cieluv = colorio.cs.CIELUV()
colorio.data.luo_rigg.show(cieluv, 50)
```

###### Ebner-Fairchild

<img src="https://nschloe.github.io/colorio/ebner-fairchild-xyy.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/ebner-fairchild-cielab.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/ebner-fairchild-cam16.svg" width="100%">
:--------------:|:---------------:|:---------------------:|
xyY             | CIELAB          |  CAM16             |

For example
```python
import colorio

colorspace = colorio.cs.JzAzBz()
colorio.data.ebner_fairchild.show(colorspace)
```
shows constant-hue data from [the Ebner-Fairchild
experiments](https://doi.org/10.1117/12.298269) in the hue-plane of some color spaces.
(Ideally, all colors in one set sit on a line.)


###### Hung-Berns
Likewise for [Hung-Berns](https://doi.org/10.1002/col.5080200506):

<img src="https://nschloe.github.io/colorio/hung-berns-xyy.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/hung-berns-cielab.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/hung-berns-cam16.svg" width="100%">
:--------------:|:---------------:|:---------------------:|
xyY             | CIELAB          |  CAM16             |

Note the dark blue distortion in CIELAB and CAM16.

```python
import colorio

colorspace = colorio.cs.JzAzBz()
colorio.data.hung_berns.show(colorspace)
```

###### Xiao et al.
Likewise for [Xiao et al.](https://doi.org/10.1002/col.20637):

<img src="https://nschloe.github.io/colorio/xiao-xyy.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/xiao-cielab.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/xiao-cam16.svg" width="100%">
:--------------:|:---------------:|:---------------------:|
xyY             | CIELAB          |  CAM16             |

```python
import colorio

colorspace = colorio.cs.CIELAB()
colorio.data.xiao.show(colorspace)
```

##### Munsell

<img src="https://nschloe.github.io/colorio/munsell-xyy.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/munsell-cielab.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/munsell-cam16.svg" width="100%">
:--------------:|:---------------:|:---------------------:|
xyY             | CIELAB          |  CAM16             |

[Munsell color data](https://www.rit.edu/cos/colorscience/rc_munsell_renotation.php) is
visualized with
```python
import colorio

cs = colorio.cs.CIELUV()
colorio.data.munsell.show(cs, V=5)
```

To retrieve the Munsell data in xyY format, use
```python
import colorio

h, V, C, xyy = colorio.data.munsell.load()
```

##### MacAdam color distances

TODO

#### Color differences

TODO

### Articles

 * [Algorithmic improvements for the CIECAM02 and CAM16 color appearance models, Nico
   Schlömer, 2018](https://arxiv.org/abs/1802.06067)
 * [On the conversion from OSA-UCS to CIEXYZ, Nico Schlömer,
   2019](https://arxiv.org/abs/1911.08323)

### Installation

colorio is [available from the Python Package Index](https://pypi.org/project/colorio/),
so just use
```
pip install colorio
```
to install.

### Testing

To run the tests, simply check out this repository and run
```
tox
```

### License
This software is published under the [GPLv3 license](https://www.gnu.org/licenses/gpl-3.0.en.html).
