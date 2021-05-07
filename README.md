<p align="center">
  <a href="https://github.com/nschloe/colorio"><img alt="colorio" src="https://nschloe.github.io/colorio/colorio-logo.svg" width="60%"></a>
  <p align="center">Tools for color models.</p>
</p>

[![PyPi Version](https://img.shields.io/pypi/v/colorio.svg?style=flat-square)](https://pypi.org/project/colorio)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/colorio.svg?style=flat-square)](https://pypi.org/pypi/colorio/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1172995.svg?style=flat-square)](https://doi.org/10.5281/zenodo.1172995)
[![GitHub stars](https://img.shields.io/github/stars/nschloe/colorio.svg?style=flat-square&logo=github&label=Stars&logoColor=white)](https://github.com/nschloe/colorio)
[![PyPi downloads](https://img.shields.io/pypi/dm/colorio.svg?style=flat-square)](https://pypistats.org/packages/colorio)

[![Discord](https://img.shields.io/static/v1?logo=discord&label=chat&message=on%20discord&color=7289da&style=flat-square)](https://discord.gg/hnTJ5MRX2Y)

[![gh-actions](https://img.shields.io/github/workflow/status/nschloe/colorio/ci?style=flat-square)](https://github.com/nschloe/colorio/actions?query=workflow%3Aci)
[![codecov](https://img.shields.io/codecov/c/github/nschloe/colorio.svg?style=flat-square)](https://codecov.io/gh/nschloe/colorio)
[![LGTM](https://img.shields.io/lgtm/grade/python/github/nschloe/colorio.svg?style=flat-square)](https://lgtm.com/projects/g/nschloe/colorio)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg?style=flat-square)](https://github.com/psf/black)

### Color spaces

All color spaces implement the two methods
<!--pytest-codeblocks:skip-->
```python
vals = colorspace.from_xyz100(xyz)
xyz = colorspace.to_xyz100(vals)
```
for conversion from and to XYZ100. Adding new color spaces is as easy as writing a class
that provides those two methods.
<!--pytest-codeblocks:skip-->
```python
colorspace.to_rgb_linear(vals)
colorspace.to_rgb1(vals)
colorspace.to_rgb255(vals)

# same for from_rgb*
```

The following color spaces are implemented:

 * [XYZ](https://en.wikipedia.org/wiki/CIE_1931_color_space) (`colorio.cs.XYZ(100)`, the
   parameter determining the scaling)
 * [xyY](https://en.wikipedia.org/wiki/CIE_1931_color_space#CIE_xy_chromaticity_diagram_and_the_CIE_xyY_color_space)
   (`colorio.cs.XYY(100)`, the paramter determining the scaling of `Y`)
 * [Linear sRGB](https://en.wikipedia.org/wiki/sRGB)  (`colorio.SrgbLinear()`)
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
 * [DIN99](https://de.wikipedia.org/wiki/DIN99-Farbraum) and its variants DIN99{b,c,d}
   (`colorio.cs.DIN99()`)
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
 * [Oklab](https://bottosson.github.io/posts/oklab/)(`colorio.cs.OKLAB()`)
 * [proLab](https://arxiv.org/abs/2012.07653)(`colorio.cs.PROLAB()`)

All methods in colorio are fully vectorized, i.e., computation is _really_ fast.


### Color difference formulas

 colorio implements the following color difference formulas:

 * [CIE76](https://en.wikipedia.org/wiki/Color_difference#CIE76)
   <!--pytest-codeblocks:skip-->
   ```python
   colorio.diff.cie76(lab1, lab2)
   ```
 * [CIE94](https://en.wikipedia.org/wiki/Color_difference#CIE94)
   <!--pytest-codeblocks:skip-->
   ```python
   colorio.diff.cie94(lab1, lab2)
   ```
 * [CIEDE2000](https://en.wikipedia.org/wiki/Color_difference#CIEDE2000)
   <!--pytest-codeblocks:skip-->
   ```python
   colorio.diff.ciede2000(lab1, lab2)
   ```
 * [CMC l:c](https://en.wikipedia.org/wiki/Color_difference#CMC_l:c_(1984))
   <!--pytest-codeblocks:skip-->
   ```python
   colorio.diff.cmc(lab1, lab2)
   ```

### Gamut visualization

colorio provides a number of useful tools for analyzing and visualizing color spaces.

#### sRGB gamut

CIELAB |  CAM16-UCS | Oklab  |
:---:|:-------:|:----------:|
<img src="https://nschloe.github.io/colorio/srgb-gamut-cielab.png" width="100%"> | <img src="https://nschloe.github.io/colorio/srgb-gamut-cam16.png" width="100%"> | <img src="https://nschloe.github.io/colorio/srgb-gamut-oklab.png" width="100%">
<img src="https://nschloe.github.io/colorio/srgb-gamut-slice-cielab.png" width="100%"> | <img src="https://nschloe.github.io/colorio/srgb-gamut-slice-cam16.png" width="100%"> | <img src="https://nschloe.github.io/colorio/srgb-gamut-slice-oklab.png" width="100%">

The sRGB gamut is a perfect cube in sRGB space, and takes curious shapes when translated
into other color spaces. The above images show the sRGB gamut in different color spaces.
<!--pytest-codeblocks:skip-->
```python
import colorio

colorspace = colorio.cs.CIELAB()
colorio.show_rgb_gamut(colorspace, n=51, show_grid=True, camera_position=None)
```
For more visualization options, you can store the sRGB data in a file 
```python
import colorio

colorspace = colorio.cs.CIELAB()
colorio.save_rgb_gamut("srgb.vtk", colorspace, n=51)
# all formats supported by https://github.com/nschloe/meshio
```
an open it with a tool of your choice. See
[here](https://github.com/nschloe/colorio/wiki/Visualizing-VTK-files) for how to open
the file in [ParaView](https://www.paraview.org/).

For lightness slices of the sRGB gamut, use
<!--pytest-codeblocks:skip-->
```python
import colorio

colorspace = colorio.cs.CIELAB()
colorio.show_rgb_slice(colorspace, lightness=50.0, n=51)
# or
# save_rgb_slice()
# plot_rgb_slice()
```
The `plot_rgb_slice()` method is especially useful for combining with other plots.

#### Surface color gamut

<img src="https://nschloe.github.io/colorio/surface-gamut-xyz.png" width="100%"> | <img src="https://nschloe.github.io/colorio/surface-gamut-cielab.png" width="100%"> | <img src="https://nschloe.github.io/colorio/surface-gamut-cam16.png" width="100%">
:---:|:-------:|:----------:|
XYZ  |  CIELAB |  CAM16-UCS |

Same as above, but with the surface color gamut visible under a given illuminant.
<!--pytest-codeblocks:skip-->
```python
import colorio

illuminant = colorio.illuminants.d65()
observer = colorio.observers.cie_1931_2()

colorspace = colorio.cs.XYZ(100)

colorio.save_surface_gamut("surface.vtk", colorspace, observer, illuminant)
colorio.show_surface_gamut(colorspace, observer, illuminant)
```
The gamut is shown in grey since sRGB screens are not able to display the colors anyway.

#### The visible gamut

xyY  |  JzAzBz |  Oklab |
:---:|:------:|:-------:|
<img src="https://nschloe.github.io/colorio/visible-gamut-xyy.png" width="100%"> | <img src="https://nschloe.github.io/colorio/visible-gamut-jzazbz.png" width="100%"> | <img src="https://nschloe.github.io/colorio/visible-gamut-oklab.png" width="100%">
<img src="https://nschloe.github.io/colorio/visible-gamut-slice-xyy.png" width="100%"> | <img src="https://nschloe.github.io/colorio/visible-gamut-slice-jzazbz.png" width="100%"> | <img src="https://nschloe.github.io/colorio/visible-gamut-slice-oklab.png" width="100%">

Same as above, but with the gamut of visible colors up to a given lightness `Y`.
<!--pytest-codeblocks:skip-->
```python
import colorio

observer = colorio.observers.cie_1931_2()

colorspace = colorio.cs.XYZ(100)

colorio.save_visible_gamut("visible.vtk", colorspace, observer, max_Y1=1)
colorio.show_visible_gamut(colorspace, observer, max_Y1=1)
```
The gamut is shown in grey since sRGB screens are not able to display the colors anyway.

For slices, use
```python
import colorio

colorspace = colorio.cs.CIELAB()
colorio.show_visible_slice(colorspace, lightness=0.5)
# or
# save_visible_slice()
# plot_visible_slice()
```


### Experimental data

colorio contains lots of experimental data sets some of which can be used to assess
certain properties of color spaces. Most data sets can also be visualized.


#### Color differences

<img src="https://nschloe.github.io/colorio/macadam1974-xyy.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/macadam1974-cielab.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/macadam1974-cam16.svg" width="100%">
:--------------:|:---------------:|:------------------:|
xyY             | CIELAB          |  CAM16             |

Color difference data from [MacAdam (1974)](https://doi.org/10.1364/JOSA.64.001691). The
above plots show the 43 color pairs that are of comparable lightness. The data is
matched perfectly if the arrow tips meet in one point.

```python
import colorio

cs = colorio.cs.CIELAB()
colorio.data.MacAdam1974().show(cs)
print(colorio.data.MacAdam1974().stress(cs))
```
```
24.531919167387624
```
The same is available for
```
colorio.data.BfdP()
colorio.data.Leeds()
colorio.data.RitDupont()
colorio.data.Witt()

colorio.data.COMBVD()  # a weighted combination of the above
```
which, combined and weighted, form the COMBVD color difference data set.


#### Munsell

<img src="https://nschloe.github.io/colorio/munsell-xyy.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/munsell-cielab.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/munsell-cam16.svg" width="100%">
:--------------:|:---------------:|:------------------:|
xyY             | CIELAB          |  CAM16             |

[Munsell color data](https://www.rit.edu/cos/colorscience/rc_munsell_renotation.php) is
visualized with
```python
import colorio

cs = colorio.cs.CIELUV()
colorio.data.Munsell().show(cs, V=5)
```

To retrieve the Munsell data in xyY format, use
```python
import colorio

munsell = colorio.data.Munsell()

# munsell.h
# munsell.V
# munsell.C
# munsell.xyy
```

#### Ellipses

##### MacAdam ellipses (1942)

<img src="https://nschloe.github.io/colorio/macadam1942-xyy.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/macadam1942-cielab.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/macadam1942-cam16.svg" width="100%">
:--------------:|:------------------:|:----------------:|
xyY (at Y=0.4)  |  CIELAB (at L=50)  |  CAM16 (at L=50) |

The famous MacAdam ellipses (from [this
article](https://doi.org/10.1364%2FJOSA.32.000247)) can be plotted with

```python
import colorio

# xyy = colorio.cs.XYY(100)
# colorio.data.MacAdam1942(50.0).show(xyy)
# colorio.data.MacAdam1942(50.0).savefig("macadam-xyy.png", xyy)

cieluv = colorio.cs.CIELUV()
colorio.data.MacAdam1942(50.0).show(cieluv)
```
The better the colorspace matches the data, the closer the ellipses are to circles of
the same size.

##### Luo-Rigg ellipses

<img src="https://nschloe.github.io/colorio/luo-rigg-xyy.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/luo-rigg-cielab.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/luo-rigg-cam16.svg" width="100%">
:---:|:------------------:|:----------------:|
xyY |  CIELAB |  CAM16 |

Likewise for [Luo-Rigg](https://doi.org/10.1002/col.5080110107).

```python
import colorio

# xyy = colorio.cs.XYY(100)
# colorio.data.LuoRigg(8).show(xyy, 0.4)
# colorio.data.LuoRigg(8).savefig("luo-rigg-xyy.png", xyy, 0.4)

cieluv = colorio.cs.CIELUV()
colorio.data.LuoRigg(8).show(cieluv, 50)
```

#### Hue linearity

##### Ebner-Fairchild

<img src="https://nschloe.github.io/colorio/ebner-fairchild-xyy.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/ebner-fairchild-cielab.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/ebner-fairchild-cam16.svg" width="100%">
:--------------:|:---------------:|:---------------------:|
xyY             | CIELAB          |  CAM16             |

For example
```python
import colorio

colorspace = colorio.cs.JzAzBz()
colorio.data.EbnerFairchild().show(colorspace)
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
colorio.data.HungBerns().show(colorspace)
```

###### Xiao et al.
Likewise for [Xiao et al.](https://doi.org/10.1002/col.20637):

<img src="https://nschloe.github.io/colorio/xiao-xyy.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/xiao-cielab.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/xiao-cam16.svg" width="100%">
:--------------:|:---------------:|:---------------------:|
xyY             | CIELAB          |  CAM16             |

```python
import colorio

colorspace = colorio.cs.CIELAB()
colorio.data.Xiao().show(colorspace)
```

#### Lightness

###### Fairchild-Chen

<img src="https://nschloe.github.io/colorio/fairchild-chen-xyy.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/fairchild-chen-cielab.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/fairchild-chen-cam16.svg" width="100%">
:--------------:|:---------------:|:---------------------:|
xyY             | CIELAB          |  CAM16             |

Lightness experiment by [Fairchild-Chen](https://doi.org/10.1117/12.872075).

```python
import colorio

cs = colorio.cs.CIELAB()
colorio.data.FairchildChen("SL2").show(cs)
```


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
