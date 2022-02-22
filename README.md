<p align="center">
  <a href="https://github.com/nschloe/colorio"><img alt="colorio" src="https://nschloe.github.io/colorio/colorio-logo.svg" width="60%"></a>
  <p align="center">Tools for color research.</p>
</p>

[![PyPi Version](https://img.shields.io/pypi/v/colorio.svg?style=flat-square)](https://pypi.org/project/colorio/)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/colorio.svg?style=flat-square)](https://pypi.org/project/colorio/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1172995.svg?style=flat-square)](https://doi.org/10.5281/zenodo.1172995)
[![GitHub stars](https://img.shields.io/github/stars/nschloe/colorio.svg?style=flat-square&logo=github&label=Stars&logoColor=white)](https://github.com/nschloe/colorio)
[![Downloads](https://pepy.tech/badge/colorio/month?style=flat-square)](https://pepy.tech/project/colorio)

<!--[![PyPi downloads](https://img.shields.io/pypi/dm/colorio.svg?style=flat-square)](https://pypistats.org/packages/colorio)-->

[![Discord](https://img.shields.io/static/v1?logo=discord&logoColor=white&label=chat&message=on%20discord&color=7289da&style=flat-square)](https://discord.gg/hnTJ5MRX2Y)

[![gh-actions](https://img.shields.io/github/workflow/status/nschloe/colorio/ci?style=flat-square)](https://github.com/nschloe/colorio/actions?query=workflow%3Aci)
[![codecov](https://img.shields.io/codecov/c/github/nschloe/colorio.svg?style=flat-square)](https://app.codecov.io/gh/nschloe/colorio)
[![LGTM](https://img.shields.io/lgtm/grade/python/github/nschloe/colorio.svg?style=flat-square)](https://lgtm.com/projects/g/nschloe/colorio)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg?style=flat-square)](https://github.com/psf/black)

### Illuminants, observers, white points

|                                Illuminants                                 |                                   CIE 1931 Observer                                    |
| :------------------------------------------------------------------------: | :------------------------------------------------------------------------------------: |
| <img src="https://nschloe.github.io/colorio/illuminants.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/cie-standard-observer-2.svg" width="100%"> |

```python
import colorio
import matplotlib.pyplot as plt

illu = colorio.illuminants.d65()
plt.plot(illu.lmbda_nm, illu.data)
plt.xlabel("wavelength [nm]")
plt.show()
```

The following illuminants are provided:

- Illuminant A ("indoor light", `colorio.illuminants.a(resolution_in_nm)`)
- Illuminant C (obsolete, "North sky daylight", `colorio.illuminants.c()`)
- Illuminants D ("natural daylight", `colorio.illuminants.d(nominal_temp)` or
  `colorio.illuminants.d65()`
  etc.)
- Illuminant E (equal energy, `colorio.illuminants.e()`)
- Illuminant series F ("fluorescent lighting", `colorio.illuminants.f2()` etc.)

Observers:

- CIE 1931 Standard 2-degree observer (`colorio.observers.colorio.observers.cie_1931_2()`)
- CIE 1964 Standard 10-degree observer (`colorio.observers.colorio.observers.cie_1964_10()`)

### Color coordinates and spaces

Color coordinates are handled as NumPy arrays or as `ColorCoordinates`, a thin
wrapper around the data that retains the color space information and has some
handy helper methods. Color spaces can be instantiated from the classes in
`colorio.cs`, e.g.,

```python
import colorio

colorio.cs.CIELAB()
```

Most methods that accept such a colorspace also accept a string, e.g.,
`cielab`.

As an example, to interpolate two sRGB colors in OKLAB, and return the sRGB:

```python
from colorio.cs import ColorCoordinates

# you can also plug in large numpy arrays instead of two lists here
c0 = ColorCoordinates([1.0, 1.0, 0.0], "srgb1")  # yellow
c1 = ColorCoordinates([0.0, 0.0, 1.0], "srgb1")  # blue

# naive interpolation gives [0.5, 0.5, 0.5], a mid gray

# convert to OKLAB
c0.convert("oklab")
c1.convert("oklab")

# interpolate
c2 = (c0 + c1) * 0.5

c2.convert("srgbhex", mode="clip")

print(c2.color_space)
print(c2.data)
```

<!--pytest-codeblocks:expected-output-->

```
<colorio color space sRGB-hex>
#6cabc7
```

All color spaces implement the two methods

<!--pytest-codeblocks:skip-->

```python
vals = colorspace.from_xyz100(xyz)
xyz = colorspace.to_xyz100(vals)
```

for conversion from and to XYZ100. Adding new color spaces is as easy as writing a class
that provides those two methods. The following color spaces are already implemented:

- [XYZ](src/colorio/cs/_xyz.py) (`colorio.cs.XYZ(100)`, the
  parameter determining the scaling)
- [xyY](src/colorio/cs/_xyy.py)
  (`colorio.cs.XYY(100)`, the parameter determining the scaling of `Y`)
- [sRGB](src/colorio/cs/_srgb.py) (`colorio.cs.SRGBlinear()`,
  `colorio.cs.SRGB1()`, `colorio.cs.SRGB255()`, `colorio.cs.SRGBhex()`)
- [HSL](src/colorio/cs/_hsl.py) and [HSV](src/colorio/cs/_hsv.py) (`colorio.cs.HSL()`,
  `colorio.cs.HSV()`)
  These classes have the two methods
  ```
  from_srgb1()
  to_srgb1()
  ```
  for conversion from and to standard RGB.
- [OSA-UCS](src/colorio/cs/_osa_ucs.py) (`colorio.cs.OsaUcs()`)
- [CIELAB](src/colorio/cs/_cielab.py) (`colorio.cs.CIELAB()`)
- [CIELUV](src/colorio/cs/_cieluv.py) (`colorio.cs.CIELUV()`)
- [RLAB](src/colorio/cs/_rlab.py) (`colorio.cs.RLAB()`)
- [DIN99](src/colorio/cs/_din99.py) and its variants DIN99{b,c,d} (`colorio.cs.DIN99()`)
- [ICtCp](src/colorio/cs/_ictcp.py) (`colorio.cs.ICtCp()`)
- [IPT](src/colorio/cs/_ipt.py) (`colorio.cs.IPT()`)
- [CIECAM02 / CAM02-UCS](src/colorio/cs/_ciecam02.py)

  ```python
  import math
  import colorio

  ciecam02 = colorio.cs.CIECAM02(0.69, 20, 100)
  cam02 = colorio.cs.CAM02("UCS", 0.69, 20, 100)
  ```

  The implementation contains a few improvements over the CIECAM02 specification (see
  [here](https://arxiv.org/abs/1802.06067)).

- [CAM16 / CAM16-UCS](src/colorio/cs/_cam16.py)

  ```python
  import math
  import colorio

  cam16 = colorio.cs.CAM16(0.69, 20, 100)
  cam16ucs = colorio.cs.CAM16UCS(0.69, 20, 100)
  ```

  The implementation contains a few improvements over the CAM16
  specification (see [here](https://arxiv.org/abs/1802.06067)).

- [J<sub>z</sub>a<sub>z</sub>b<sub>z</sub>](https://doi.org/10.1364/OE.25.015131)
  (`colorio.cs.JzAzBz()`)
- [Oklab](src/colorio/cs/_oklab.py) (`colorio.cs.OKLAB()`)
- [proLab](src/colorio/cs/_prolab.py) (`colorio.cs.PROLAB()`)
- [SRLAB2](src/colorio/cs/_srlab2.py) (`colorio.cs.SRLAB2()`)

All methods in colorio are fully vectorized, i.e., computation is _really_ fast.

### Color difference formulas

colorio implements the following color difference formulas:

- [CIE76](src/colorio/diff/_cie76.py)
  <!--pytest-codeblocks:skip-->
  ```python
  colorio.diff.cie76(lab1, lab2)
  ```
- [CIE94](src/colorio/diff/_cie94.py)
  <!--pytest-codeblocks:skip-->
  ```python
  colorio.diff.cie94(lab1, lab2)
  ```
- [CIEDE2000](src/colorio/diff/_ciede2000.py)
  <!--pytest-codeblocks:skip-->
  ```python
  colorio.diff.ciede2000(lab1, lab2)
  ```
- [CMC l:c](src/colorio/diff/_cmc.py)
  <!--pytest-codeblocks:skip-->
  ```python
  colorio.diff.cmc(lab1, lab2)
  ```

### Chromatic adaptation transforms

colorio implements the following CATs:

- [von Kries](src/colorio/cat/von_kries.py)
  <!--pytest-codeblocks:skip-->
  ```python
  cat, cat_inv = colorio.cat.von_kries(whitepoint_source, whitepoint_destination)
  xyz1 = cat @ xyz0
  ```
- [Bradford](src/colorio/cat/bradford.py) (`colorio.cat.bradford`)
- [sharp](src/colorio/cat/sharp.py) (`colorio.cat.sharp`)
- [CMCCAT2000](src/colorio/cat/cmccat2000.py) (`colorio.cat.cmccat2000`)
- [CAT02](src/colorio/cat/cat02.py) (`colorio.cat.cat02`)
- [CAT16](src/colorio/cat/cat16.py) (`colorio.cat.cat16`)
- [Bianco-Schettini](src/colorio/cat/bianco_schettini.py) (`colorio.cat.bianco_schettini`)

### Gamut visualization

colorio provides a number of useful tools for analyzing and visualizing color spaces.

#### sRGB gamut

|                                         CIELAB                                         |                                        CAM16-UCS                                         |                                         Oklab                                         |
| :------------------------------------------------------------------------------------: | :--------------------------------------------------------------------------------------: | :-----------------------------------------------------------------------------------: |
|    <img src="https://nschloe.github.io/colorio/srgb-gamut-cielab.png" width="100%">    |     <img src="https://nschloe.github.io/colorio/srgb-gamut-cam16.png" width="100%">      |    <img src="https://nschloe.github.io/colorio/srgb-gamut-oklab.png" width="100%">    |
| <img src="https://nschloe.github.io/colorio/srgb-gamut-slice-cielab.png" width="100%"> | <img src="https://nschloe.github.io/colorio/srgb-gamut-slice-cam16ucs.png" width="100%"> | <img src="https://nschloe.github.io/colorio/srgb-gamut-slice-oklab.png" width="100%"> |

The sRGB gamut is a perfect cube in sRGB space, and takes curious shapes when translated
into other color spaces. The above images show the sRGB gamut in different color spaces.

<!--pytest-codeblocks:skip-->

```python
import colorio

p = colorio.plot_rgb_gamut(
    "cielab",  # or colorio.cs.CIELAB()
    n=51,
    show_grid=True,
)
p.show()
```

For more visualization options, you can store the sRGB data in a file

```python
import colorio

colorio.save_rgb_gamut("srgb.vtk", "cielab", n=51)
# all formats supported by https://github.com/nschloe/meshio
```

and open it with a tool of your choice. See
[here](https://github.com/nschloe/colorio/wiki/Visualizing-VTK-files) for how to open
the file in [ParaView](https://www.paraview.org/).

For lightness slices of the sRGB gamut, use

<!--pytest-codeblocks:skip-->

```python
import colorio

p = colorio.plot_rgb_slice("cielab", lightness=50.0, n=51)
p.show()
# or
# p.screenshot("screenshot.png")
```

#### Surface color gamut

| <img src="https://nschloe.github.io/colorio/surface-gamut-xyz.png" width="100%"> | <img src="https://nschloe.github.io/colorio/surface-gamut-cielab.png" width="100%"> | <img src="https://nschloe.github.io/colorio/surface-gamut-cam16.png" width="100%"> |
| :------------------------------------------------------------------------------: | :---------------------------------------------------------------------------------: | :--------------------------------------------------------------------------------: |
|                                       XYZ                                        |                                       CIELAB                                        |                                     CAM16-UCS                                      |

Same as above, but with the surface color gamut visible under a given illuminant.

<!--pytest-codeblocks:skip-->

```python
import colorio

illuminant = colorio.illuminants.d65()
observer = colorio.observers.cie_1931_2()

p = colorio.plot_surface_gamut(
    "xyz100",  # or colorio.cs.XYZ(100)
    observer,
    illuminant,
)
p.show()
```

The gamut is shown in grey since sRGB screens are not able to display the colors anyway.

#### The visible gamut

|                                          xyY                                           |                                          JzAzBz                                           |                                          Oklab                                           |
| :------------------------------------------------------------------------------------: | :---------------------------------------------------------------------------------------: | :--------------------------------------------------------------------------------------: |
|    <img src="https://nschloe.github.io/colorio/visible-gamut-xyy.png" width="100%">    |    <img src="https://nschloe.github.io/colorio/visible-gamut-jzazbz.png" width="100%">    |    <img src="https://nschloe.github.io/colorio/visible-gamut-oklab.png" width="100%">    |
| <img src="https://nschloe.github.io/colorio/visible-gamut-slice-xyy.png" width="100%"> | <img src="https://nschloe.github.io/colorio/visible-gamut-slice-jzazbz.png" width="100%"> | <img src="https://nschloe.github.io/colorio/visible-gamut-slice-oklab.png" width="100%"> |

Same as above, but with the gamut of visible colors up to a given lightness `Y`.

<!--pytest-codeblocks:skip-->

```python
import colorio

observer = colorio.observers.cie_1931_2()

colorspace = colorio.cs.XYZ(100)

p = colorio.plot_visible_gamut(colorspace, observer, max_Y1=1)
p.show()
```

The gamut is shown in grey since sRGB screens are not able to display the colors anyway.

For slices, use

```python
import colorio

plt = colorio.plot_visible_slice("cielab", lightness=0.5)
plt.show()
```

### Color gradients

With colorio, you can easily visualize the basic color gradients of any color space.
This may make defects in color spaces obvious, e.g., the well-known blue-distortion of
CIELAB and related spaces. (Compare with [the hue linearity data
below](#hue-linearity).)

```python
import colorio

plt = colorio.plot_primary_srgb_gradients("cielab")
plt.show()
```

| <img src="https://nschloe.github.io/colorio/gradients-cielab.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/gradients-din99.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/gradients-oklab.svg" width="100%"> |
| :-----------------------------------------------------------------------------: | :----------------------------------------------------------------------------: | :----------------------------------------------------------------------------: |
|                                     CIELAB                                      |                                     DIN99                                      |                                     OKLAB                                      |

### Experimental data

colorio contains lots of experimental data sets some of which can be used to assess
certain properties of color spaces. Most data sets can also be visualized.

#### Color differences

| <img src="https://nschloe.github.io/colorio/macadam1974-xyy1.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/macadam1974-cielab.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/macadam1974-cam16ucs.svg" width="100%"> |
| :-----------------------------------------------------------------------------: | :-------------------------------------------------------------------------------: | :---------------------------------------------------------------------------------: |
|                                       xyY                                       |                                      CIELAB                                       |                                        CAM16                                        |

Color difference data from [MacAdam (1974)](https://doi.org/10.1364/JOSA.64.001691). The
above plots show the 43 color pairs that are of comparable lightness. The data is
matched perfectly if the facing line stubs meet in one point.

```python
import colorio

data = colorio.data.MacAdam1974()

cs = colorio.cs.CIELAB

plt = data.plot(cs)
plt.show()
print(colorio.data.MacAdam1974().stress(cs))
```

```
24.54774029343344
```

The same is available for

```
colorio.data.BfdP()
colorio.data.Leeds()
colorio.data.RitDupont()
colorio.data.Witt()

colorio.data.COMBVD()  # a weighted combination of the above
```

#### Munsell

| <img src="https://nschloe.github.io/colorio/munsell-xyy1.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/munsell-cielab.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/munsell-cam16ucs.svg" width="100%"> |
| :-------------------------------------------------------------------------: | :---------------------------------------------------------------------------: | :-----------------------------------------------------------------------------: |
|                                     xyY                                     |                                    CIELAB                                     |                                      CAM16                                      |

[Munsell color data](https://www.rit.edu/cos/colorscience/rc_munsell_renotation.php) is
visualized with

```python
import colorio

cs = colorio.cs.CIELUV
plt = colorio.data.Munsell().plot(cs, V=5)
plt.show()
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

| <img src="https://nschloe.github.io/colorio/macadam1942-xyy.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/macadam1942-cielab.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/macadam1942-cam16.svg" width="100%"> |
| :----------------------------------------------------------------------------: | :-------------------------------------------------------------------------------: | :------------------------------------------------------------------------------: |
|                                 xyY (at Y=0.4)                                 |                                 CIELAB (at L=50)                                  |                                 CAM16 (at L=50)                                  |

The famous MacAdam ellipses (from [this
article](https://doi.org/10.1364%2FJOSA.32.000247)) can be plotted with

```python
import colorio

cs = colorio.cs.CIELUV
plt = colorio.data.MacAdam1942(50.0).plot(cs)
plt.show()
```

The better the colorspace matches the data, the closer the ellipses are to circles of
the same size.

##### Luo-Rigg ellipses

| <img src="https://nschloe.github.io/colorio/luo-rigg-xyy.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/luo-rigg-cielab.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/luo-rigg-cam16ucs.svg" width="100%"> |
| :-------------------------------------------------------------------------: | :----------------------------------------------------------------------------: | :------------------------------------------------------------------------------: |
|                                     xyY                                     |                                     CIELAB                                     |                                      CAM16                                       |

Likewise for [Luo-Rigg](https://doi.org/10.1002/col.5080110107).

```python
import colorio

# xyy = colorio.cs.XYY(100)
# colorio.data.LuoRigg(8).show(xyy, 0.4)
# colorio.data.LuoRigg(8).savefig("luo-rigg-xyy.png", xyy, 0.4)

cieluv = colorio.cs.CIELUV()
plt = colorio.data.LuoRigg(8).plot(cieluv, 50)
plt.show()
```

#### Hue linearity

##### Ebner-Fairchild

| <img src="https://nschloe.github.io/colorio/ebner-fairchild-xyy1.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/ebner-fairchild-cielab.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/ebner-fairchild-cam16ucs.svg" width="100%"> |
| :---------------------------------------------------------------------------------: | :-----------------------------------------------------------------------------------: | :-------------------------------------------------------------------------------------: |
|                                         xyY                                         |                                        CIELAB                                         |                                          CAM16                                          |

For example

```python
import colorio

colorspace = colorio.cs.JzAzBz
plt = colorio.data.EbnerFairchild().plot(colorspace)
plt.show()
```

shows constant-hue data from [the Ebner-Fairchild
experiments](https://doi.org/10.1117/12.298269) in the hue-plane of some color spaces.
(Ideally, all colors in one set sit on a line.)

###### Hung-Berns

Likewise for [Hung-Berns](https://doi.org/10.1002/col.5080200506):

| <img src="https://nschloe.github.io/colorio/hung-berns-xyy.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/hung-berns-cielab.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/hung-berns-cam16.svg" width="100%"> |
| :---------------------------------------------------------------------------: | :------------------------------------------------------------------------------: | :-----------------------------------------------------------------------------: |
|                                      xyY                                      |                                      CIELAB                                      |                                      CAM16                                      |

Note the dark blue distortion in CIELAB and CAM16.

```python
import colorio

colorspace = colorio.cs.JzAzBz
plt = colorio.data.HungBerns().plot(colorspace)
plt.show()
```

###### Xiao et al.

Likewise for [Xiao et al.](https://doi.org/10.1002/col.20637):

| <img src="https://nschloe.github.io/colorio/xiao-xyy.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/xiao-cielab.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/xiao-cam16.svg" width="100%"> |
| :---------------------------------------------------------------------: | :------------------------------------------------------------------------: | :-----------------------------------------------------------------------: |
|                                   xyY                                   |                                   CIELAB                                   |                                   CAM16                                   |

```python
import colorio

colorspace = colorio.cs.CIELAB
plt = colorio.data.Xiao().plot(colorspace)
plt.show()
```

#### Lightness

###### Fairchild-Chen

| <img src="https://nschloe.github.io/colorio/fairchild-chen-xyy.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/fairchild-chen-cielab.svg" width="100%"> | <img src="https://nschloe.github.io/colorio/fairchild-chen-cam16.svg" width="100%"> |
| :-------------------------------------------------------------------------------: | :----------------------------------------------------------------------------------: | :---------------------------------------------------------------------------------: |
|                                        xyY                                        |                                        CIELAB                                        |                                        CAM16                                        |

Lightness experiment by [Fairchild-Chen](https://doi.org/10.1117/12.872075).

```python
import colorio

cs = colorio.cs.CIELAB
plt = colorio.data.FairchildChen("SL2").plot(cs)
plt.show()
```

### Articles

- [Algorithmic improvements for the CIECAM02 and CAM16 color appearance models, Nico
  Schlömer, 2018](https://arxiv.org/abs/1802.06067)
- [On the conversion from OSA-UCS to CIEXYZ, Nico Schlömer,
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
