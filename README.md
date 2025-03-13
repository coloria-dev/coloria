<p align="center">
  <a href="https://github.com/coloria-dev/coloria"><img alt="coloria" src="https://raw.githubusercontent.com/coloria-dev/coloria/main/logo/coloria-logo.png" width="60%"></a>
  <p align="center">Tools for color research.</p>
</p>

[![PyPi Version](https://img.shields.io/pypi/v/coloria.svg?style=flat-square)](https://pypi.org/project/coloria/)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/coloria.svg?style=flat-square)](https://pypi.org/project/coloria/)
[![GitHub stars](https://img.shields.io/github/stars/coloria-dev/coloria.svg?style=flat-square&logo=github&label=Stars&logoColor=white)](https://github.com/coloria-dev/coloria)
[![Downloads](https://pepy.tech/badge/coloria/month?style=flat-square)](https://pepy.tech/project/coloria)

<!--[![PyPi downloads](https://img.shields.io/pypi/dm/coloria.svg?style=flat-square)](https://pypistats.org/packages/coloria)-->

[![Discord](https://img.shields.io/static/v1?logo=discord&logoColor=white&label=chat&message=on%20discord&color=7289da&style=flat-square)](https://discord.gg/hnTJ5MRX2Y)

### Installation

Install Coloria [from PyPI](https://pypi.org/project/coloria/) with

```
pip install coloria
```

To run Coloria, you need a license. See [here](https://github.com/coloria-dev)
for more info.

### Illuminants, observers, white points

|                                                Illuminants                                                |                                                   CIE 1931 Observer                                                   |
| :-------------------------------------------------------------------------------------------------------: | :-------------------------------------------------------------------------------------------------------------------: |
| <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/illuminants.svg" width="100%"> | <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/cie-standard-observer-2.svg" width="100%"> |

```python
import coloria
import matplotlib.pyplot as plt

illu = coloria.illuminants.d65()
plt.plot(illu.lmbda_nm, illu.data)
plt.xlabel("wavelength [nm]")
plt.show()
```

The following illuminants are provided:

- Illuminant A ("indoor light", `coloria.illuminants.a(resolution_in_nm)`)
- Illuminant C (obsolete, "North sky daylight", `coloria.illuminants.c()`)
- Illuminants D ("natural daylight", `coloria.illuminants.d(nominal_temp)` or
  `coloria.illuminants.d65()`
  etc.)
- Illuminant E (equal energy, `coloria.illuminants.e()`)
- Illuminant series F ("fluorescent lighting", `coloria.illuminants.f2()` etc.)

Observers:

- CIE 1931 Standard 2-degree observer (`coloria.observers.coloria.observers.cie_1931_2()`)
- CIE 1964 Standard 10-degree observer (`coloria.observers.coloria.observers.cie_1964_10()`)

### Color appearance models

Color appearance models (CAMs) predicts all kinds of parameters in color perception,
e.g., lightness, brightness, chroma, colorfulness, saturation etc. Since these
values depend on various factors, such as the surrouning, the models are initialized
with various different parameters.

CAMs can be used to construct color _spaces_ (see below).

The color appearance models available in coloria are

- CIECAM02 / CAM02-UCS

  ```python
  import coloria

  ciecam02 = coloria.cam.CIECAM02(
      surround_type="average",  # or dim, dark
      background_luminance_percent=20,
      adapting_luminance_cd_m2=100,
  )

  xyz = [19.31, 23.93, 10.14]
  corr = ciecam02.from_xyz100(xyz)
  # then work with those values:
  corr.lightness
  corr.brightness
  corr.chroma
  corr.hue_composition
  corr.hue_angle_degrees
  corr.colorfulness
  corr.saturation
  ```

- CAM16 / CAM16-UCS

  ```python
  import coloria

  cam16 = coloria.cam.CAM16("average", 20, 100)
  ```

- ZCAM

  ```python
  import coloria

  cam16 = coloria.cam.ZCAM("average", 20, 100, 20)
  ```

### Color coordinates and spaces

Color coordinates are handled as NumPy arrays or as `ColorCoordinates`, a thin
wrapper around the data that retains the color space information and has some
handy helper methods. Color spaces can be instantiated from the classes in
`coloria.cs`, e.g.,

```python
import coloria

coloria.cs.CIELAB()
```

Most methods that accept such a colorspace also accept a string, e.g.,
`cielab`.

As an example, to interpolate two sRGB colors in OKLAB, and return the sRGB:

```python
from coloria.cs import ColorCoordinates

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
<coloria color space sRGB-hex>
#6cabc7
```

All color spaces implement the two methods

<!--pytest.mark.skip-->

```python
vals = colorspace.from_xyz100(xyz)
xyz = colorspace.to_xyz100(vals)
```

for conversion from and to XYZ100. Adding new color spaces is as easy as writing a class
that provides those two methods. The following color spaces are already implemented:

- XYZ (`coloria.cs.XYZ(100)`, the
  parameter determining the scaling)
- xyY
  (`coloria.cs.XYY(100)`, the parameter determining the scaling of `Y`)
- sRGB (`coloria.cs.SRGBlinear()`,
  `coloria.cs.SRGB1()`, `coloria.cs.SRGB255()`, `coloria.cs.SRGBhex()`)
- HSL and HSV (`coloria.cs.HSL()`,
  `coloria.cs.HSV()`)
  These classes also have the two methods
  ```
  from_srgb1()
  to_srgb1()
  ```
  for direct conversion from and to standard RGB.
- [OSA-UCS (`coloria.cs.OsaUcs()`)](https://en.wikipedia.org/wiki/OSA-UCS), 1947
- CIE 1960 UCS (`coloria.cs.CIE1960UCS()`), 1960
- CIEUVW (`coloria.cs.CIEUVW()`), 1964
- CIELAB (`coloria.cs.CIELAB()`), 1976
- CIELUV (`coloria.cs.CIELUV()`), 1976
- [RLAB (`coloria.cs.RLAB()`)](https://doi.org/10.1117/12.149061), 1993
- [IPT
  (`coloria.cs.IPT()`)](https://www.ingentaconnect.com/content/ist/cic/1998/00001998/00000001/art00003),
  1998
- DIN99 and its variants DIN99{b,c,d} (`coloria.cs.DIN99()`), 1999
- CAM02-UCS, 2002

  ```python
  import coloria

  cam02 = coloria.cs.CAM02("UCS", "average", 20, 100)
  ```

  The implementation contains a few improvements over the CIECAM02
  specification (see [here](https://arxiv.org/abs/1802.06067)).

- CAM16-UCS, 2016

  ```python
  import coloria

  cam16ucs = coloria.cs.CAM16UCS("average", 20, 100)
  ```

  The implementation contains a few improvements over the CAM16
  specification (see [here](https://arxiv.org/abs/1802.06067)).

- SRLAB2 (`coloria.cs.SRLAB2()`)
- [J<sub>z</sub>a<sub>z</sub>b<sub>z</sub>](https://doi.org/10.1364/OE.25.015131)
  (`coloria.cs.JzAzBz()`), 2017
- [ICtCp (`coloria.cs.ICtCp()`)](https://en.wikipedia.org/wiki/ICtCp), 2018
- [I<sub>G</sub>P<sub>G</sub>T<sub>G</sub>
  (`coloria.cs.IGPGTG()`)](https://doi.org/10.2352/J.Percept.Imaging.2020.3.2.020401),
  2020
- [proLab (`coloria.cs.PROLAB()`)](https://arxiv.org/abs/2012.07653), 2020
- [Oklab (`coloria.cs.OKLAB()`)](https://bottosson.github.io/posts/oklab/), 2020
- OkLCh (`coloria.cs.OKLCH()`), 2020
- [HCT (`coloria.cs.HCT()`/ HCTLAB
  (`coloria.cs.HCTLAB()`)](https://material.io/blog/science-of-color-design),
  2022

All methods in coloria are fully vectorized, i.e., computation is _really_
fast.

### Color difference formulas

coloria implements the following color difference formulas:

- CIE76
  <!--pytest.mark.skip-->
  ```python
  coloria.diff.cie76(lab1, lab2)
  ```
- CIE94
  <!--pytest.mark.skip-->
  ```python
  coloria.diff.cie94(lab1, lab2)
  ```
- CIEDE2000
  <!--pytest.mark.skip-->
  ```python
  coloria.diff.ciede2000(lab1, lab2)
  ```
- CMC l:c
  <!--pytest.mark.skip-->
  ```python
  coloria.diff.cmc(lab1, lab2)
  ```

### Chromatic adaptation transforms

coloria implements the following CATs:

- von Kries
  <!--pytest.mark.skip-->
  ```python
  cat, cat_inv = coloria.cat.von_kries(whitepoint_source, whitepoint_destination)
  xyz1 = cat @ xyz0
  ```
- Bradford (`coloria.cat.bradford`)
- sharp (`coloria.cat.sharp`)
- CMCCAT2000 (`coloria.cat.cmccat2000`)
- CAT02 (`coloria.cat.cat02`)
- CAT16 (`coloria.cat.cat16`)
- Bianco-Schettini (`coloria.cat.bianco_schettini`)

### Gamut visualization

coloria provides a number of useful tools for analyzing and visualizing color spaces.

#### sRGB gamut

|                                                        CIELAB                                                         |                                                        CAM16-UCS                                                        |                                                        Oklab                                                         |
| :-------------------------------------------------------------------------------------------------------------------: | :---------------------------------------------------------------------------------------------------------------------: | :------------------------------------------------------------------------------------------------------------------: |
|    <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/srgb-gamut-cielab.png" width="100%">    |     <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/srgb-gamut-cam16.png" width="100%">      |    <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/srgb-gamut-oklab.png" width="100%">    |
| <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/srgb-gamut-slice-cielab.png" width="100%"> | <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/srgb-gamut-slice-cam16ucs.png" width="100%"> | <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/srgb-gamut-slice-oklab.png" width="100%"> |

The sRGB gamut is a perfect cube in sRGB space, and takes curious shapes when translated
into other color spaces. The above images show the sRGB gamut in different color spaces.

<!--pytest.mark.skip-->

```python
import coloria

p = coloria.plot_rgb_gamut(
    "cielab",  # or coloria.cs.CIELAB()
    n=51,
    show_grid=True,
)
p.show()
```

For more visualization options, you can store the sRGB data in a file

```python
import coloria

coloria.save_rgb_gamut("srgb.vtk", "cielab", n=51)
# all formats supported by https://github.com/coloria-dev/meshio
```

and open it with a tool of your choice. See
[here](https://github.com/coloria-dev/coloria/wiki/Visualizing-VTK-files) for how to open
the file in [ParaView](https://www.paraview.org/).

For lightness slices of the sRGB gamut, use

<!--pytest.mark.skip-->

```python
import coloria

p = coloria.plot_rgb_slice("cielab", lightness=50.0, n=51)
p.show()
# or
# p.screenshot("screenshot.png")
```

#### Surface color gamut

| <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/surface-gamut-xyz.png" width="100%"> | <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/surface-gamut-cielab.png" width="100%"> | <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/surface-gamut-cam16.png" width="100%"> |
| :-------------------------------------------------------------------------------------------------------------: | :----------------------------------------------------------------------------------------------------------------: | :---------------------------------------------------------------------------------------------------------------: |
|                                                       XYZ                                                       |                                                       CIELAB                                                       |                                                     CAM16-UCS                                                     |

Same as above, but with the surface color gamut visible under a given illuminant.

<!--pytest.mark.skip-->

```python
import coloria

illuminant = coloria.illuminants.d65()
observer = coloria.observers.cie_1931_2()

p = coloria.plot_surface_gamut(
    "xyz100",  # or coloria.cs.XYZ(100)
    observer,
    illuminant,
)
p.show()
```

The gamut is shown in grey since sRGB screens are not able to display the colors anyway.

#### The visible gamut

|                                                          xyY                                                          |                                                          JzAzBz                                                          |                                                          Oklab                                                          |
| :-------------------------------------------------------------------------------------------------------------------: | :----------------------------------------------------------------------------------------------------------------------: | :---------------------------------------------------------------------------------------------------------------------: |
|    <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/visible-gamut-xyy.png" width="100%">    |    <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/visible-gamut-jzazbz.png" width="100%">    |    <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/visible-gamut-oklab.png" width="100%">    |
| <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/visible-gamut-slice-xyy.png" width="100%"> | <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/visible-gamut-slice-jzazbz.png" width="100%"> | <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/visible-gamut-slice-oklab.png" width="100%"> |

Same as above, but with the gamut of visible colors up to a given lightness `Y`.

<!--pytest.mark.skip-->

```python
import coloria

observer = coloria.observers.cie_1931_2()

colorspace = coloria.cs.XYZ(100)

p = coloria.plot_visible_gamut(colorspace, observer, max_Y1=1)
p.show()
```

The gamut is shown in grey since sRGB screens are not able to display the colors anyway.

For slices, use

```python
import coloria

plt = coloria.plot_visible_slice("cielab", lightness=0.5)
plt.show()
```

### Color gradients

With coloria, you can easily visualize the basic color gradients of any color space.
This may make defects in color spaces obvious, e.g., the well-known blue-distortion of
CIELAB and related spaces. (Compare with [the hue linearity data
below](#hue-linearity).)

```python
import coloria

plt = coloria.plot_primary_srgb_gradients("cielab")
plt.show()
```

| <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/gradients-cielab.svg" width="100%"> | <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/gradients-din99.svg" width="100%"> | <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/gradients-oklab.svg" width="100%"> |
| :------------------------------------------------------------------------------------------------------------: | :-----------------------------------------------------------------------------------------------------------: | :-----------------------------------------------------------------------------------------------------------: |
|                                                     CIELAB                                                     |                                                     DIN99                                                     |                                                     OKLAB                                                     |

### Experimental data

coloria contains lots of experimental data sets some of which can be used to assess
certain properties of color spaces. Most data sets can also be visualized.

#### Color differences

| <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/macadam1974-xyy1.svg" width="100%"> | <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/macadam1974-cielab.svg" width="100%"> | <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/macadam1974-cam16ucs.svg" width="100%"> |
| :------------------------------------------------------------------------------------------------------------: | :--------------------------------------------------------------------------------------------------------------: | :----------------------------------------------------------------------------------------------------------------: |
|                                                      xyY                                                       |                                                      CIELAB                                                      |                                                       CAM16                                                        |

Color difference data from [MacAdam (1974)](https://doi.org/10.1364/JOSA.64.0"average"1). The
above plots show the 43 color pairs that are of comparable lightness. The data is
matched perfectly if the facing line stubs meet in one point.

```python
import coloria

data = coloria.data.MacAdam1974()

cs = coloria.cs.CIELAB

plt = data.plot(cs)
plt.show()
print(coloria.data.MacAdam1974().stress(cs))
```

```
24.54774029343344
```

The same is available for

```
coloria.data.BfdP()
coloria.data.Leeds()
coloria.data.RitDupont()
coloria.data.Witt()

coloria.data.COMBVD()  # a weighted combination of the above
```

#### Munsell

| <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/munsell-xyy1.svg" width="100%"> | <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/munsell-cielab.svg" width="100%"> | <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/munsell-cam16ucs.svg" width="100%"> |
| :--------------------------------------------------------------------------------------------------------: | :----------------------------------------------------------------------------------------------------------: | :------------------------------------------------------------------------------------------------------------: |
|                                                    xyY                                                     |                                                    CIELAB                                                    |                                                     CAM16                                                      |

[Munsell color data](https://www.rit.edu/cos/colorscience/rc_munsell_renotation.php) is
visualized with

```python
import coloria

cs = coloria.cs.CIELUV
plt = coloria.data.Munsell().plot(cs, V=5)
plt.show()
```

To retrieve the Munsell data in xyY format, use

```python
import coloria

munsell = coloria.data.Munsell()

# munsell.h
# munsell.V
# munsell.C
# munsell.xyy
```

#### Ellipses

##### MacAdam ellipses (1942)

| <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/macadam1942-xyy.svg" width="100%"> | <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/macadam1942-cielab.svg" width="100%"> | <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/macadam1942-cam16.svg" width="100%"> |
| :-----------------------------------------------------------------------------------------------------------: | :--------------------------------------------------------------------------------------------------------------: | :-------------------------------------------------------------------------------------------------------------: |
|                                                xyY (at Y=0.4)                                                 |                                                 CIELAB (at L=50)                                                 |                                                 CAM16 (at L=50)                                                 |

The famous MacAdam ellipses (from [this
article](https://doi.org/10.1364%2FJOSA.32.000247)) can be plotted with

```python
import coloria

cs = coloria.cs.CIELUV
plt = coloria.data.MacAdam1942(50.0).plot(cs)
plt.show()
```

The better the colorspace matches the data, the closer the ellipses are to circles of
the same size.

##### Luo-Rigg ellipses

| <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/luo-rigg-xyy.svg" width="100%"> | <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/luo-rigg-cielab.svg" width="100%"> | <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/luo-rigg-cam16ucs.svg" width="100%"> |
| :--------------------------------------------------------------------------------------------------------: | :-----------------------------------------------------------------------------------------------------------: | :-------------------------------------------------------------------------------------------------------------: |
|                                                    xyY                                                     |                                                    CIELAB                                                     |                                                      CAM16                                                      |

Likewise for [Luo-Rigg](https://doi.org/10.1002/col.5080110107).

```python
import coloria

# xyy = coloria.cs.XYY(100)
# coloria.data.LuoRigg(8).show(xyy, 0.4)
# coloria.data.LuoRigg(8).savefig("luo-rigg-xyy.png", xyy, 0.4)

cieluv = coloria.cs.CIELUV()
plt = coloria.data.LuoRigg(8).plot(cieluv, 50)
plt.show()
```

#### Hue linearity

##### Ebner-Fairchild

| <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/ebner-fairchild-xyy1.svg" width="100%"> | <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/ebner-fairchild-cielab.svg" width="100%"> | <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/ebner-fairchild-cam16ucs.svg" width="100%"> |
| :----------------------------------------------------------------------------------------------------------------: | :------------------------------------------------------------------------------------------------------------------: | :--------------------------------------------------------------------------------------------------------------------: |
|                                                        xyY                                                         |                                                        CIELAB                                                        |                                                         CAM16                                                          |

For example

```python
import coloria

colorspace = coloria.cs.JzAzBz
plt = coloria.data.EbnerFairchild().plot(colorspace)
plt.show()
```

shows constant-hue data from [the Ebner-Fairchild
experiments](https://doi.org/10.1117/12.298269) in the hue-plane of some color spaces.
(Ideally, all colors in one set sit on a line.)

###### Hung-Berns

Likewise for [Hung-Berns](https://doi.org/10.1002/col.5080200506):

| <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/hung-berns-xyy.svg" width="100%"> | <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/hung-berns-cielab.svg" width="100%"> | <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/hung-berns-cam16.svg" width="100%"> |
| :----------------------------------------------------------------------------------------------------------: | :-------------------------------------------------------------------------------------------------------------: | :------------------------------------------------------------------------------------------------------------: |
|                                                     xyY                                                      |                                                     CIELAB                                                      |                                                     CAM16                                                      |

Note the dark blue distortion in CIELAB and CAM16.

```python
import coloria

colorspace = coloria.cs.JzAzBz
plt = coloria.data.HungBerns().plot(colorspace)
plt.show()
```

###### Xiao et al.

Likewise for [Xiao et al.](https://doi.org/10.1002/col.20637):

| <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/xiao-xyy.svg" width="100%"> | <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/xiao-cielab.svg" width="100%"> | <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/xiao-cam16.svg" width="100%"> |
| :----------------------------------------------------------------------------------------------------: | :-------------------------------------------------------------------------------------------------------: | :------------------------------------------------------------------------------------------------------: |
|                                                  xyY                                                   |                                                  CIELAB                                                   |                                                  CAM16                                                   |

```python
import coloria

colorspace = coloria.cs.CIELAB
plt = coloria.data.Xiao().plot(colorspace)
plt.show()
```

#### Lightness

###### Fairchild-Chen

| <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/fairchild-chen-xyy.svg" width="100%"> | <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/fairchild-chen-cielab.svg" width="100%"> | <img src="https://raw.githubusercontent.com/coloria-dev/coloria/main/plots/fairchild-chen-cam16.svg" width="100%"> |
| :--------------------------------------------------------------------------------------------------------------: | :-----------------------------------------------------------------------------------------------------------------: | :----------------------------------------------------------------------------------------------------------------: |
|                                                       xyY                                                        |                                                       CIELAB                                                        |                                                       CAM16                                                        |

Lightness experiment by [Fairchild-Chen](https://doi.org/10.1117/12.872075).

```python
import coloria

cs = coloria.cs.CIELAB
plt = coloria.data.FairchildChen("SL2").plot(cs)
plt.show()
```

### Articles

- [Algorithmic improvements for the CIECAM02 and CAM16 color appearance models,
  Nico Schlömer, 2018](https://arxiv.org/abs/1802.06067)
- [On the conversion from OSA-UCS to CIEXYZ, Nico Schlömer,
  2019](https://arxiv.org/abs/1911.08323)
