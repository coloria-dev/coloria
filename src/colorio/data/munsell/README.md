Text and data retrieved from
<https://www.rit.edu/cos/colorscience/rc_munsell_renotation.php>.

## Munsell Renotation Data

### Overview and History

In the 1940's the color science community recognized that the most
visually-uniform color space to date, the Munsell Color Order System, had
inconsistencies that required examination and remedy. Towards this goal, a
large-scale visual experiment was taken with many observers across several
continents. The results amounted to an adjustment of the target color
coordinates for the Munsell colors. The files here reflect that correction.

There are three files available for download. All are of the same format: six
columns of Munsell hue, Munsell value, Munsell chroma, CIE x, y, and Y. The
chromaticity coordinates were calculated using illuminant C and the CIE 1931 2
degree observer. In a sense, all three files represent the same set of data, in
that all depend on the scaling experiments of the late 1930's.

A report entitled "One Set of Munsell Re-renotations," by Deane B. Judd and
Dorothy Nickerson was issues by the National Bureau of Standards (now NIST) in
1967. As the title implies, they proposed an alternative to the original
renotation scheme. As far as we know, these did not receive much attention, and
their utility is uncertain. The report and associated data table have been
scanned. If you use this please let us know! We would be interested in any
useful application of the report of data.

**None of these data should be confused with actual measurements from a Munsell
Book of Color!**

* `all.dat`: real and unreal

  These are all the Munsell data, including the extrapolated colors. Note that
  extrapolated colors are in some cases unreal. That is, some lie outsize the
  Macadam limits.

  This file should be used for those performing multidimensional interpolation
  to/from Munsell data. You will need the unreal colors in order to completely
  encompass the real colors, which is required to do the interpolation when near
  the Macadam limits.

* `real.dat`: by the book

  These are real colors only, "real" being those lying inside the Macadam
  limits.  Specifically, these are those colors listed the original 1943
  renotation article (Newhall, Judd, and Nickerson, JOSA, 1943).

  This file should be used for a complete mapping between the Munsell system and
  its CIE equivalents. Note, however, that many of these colors were not used in
  the original scaling experiments, and are therefore extrapolated or at best
  interpolated from the test colors used.

  Flash! Here are sRGB values and CIELAB for most of the colors in the real.dat
  file. There are some important notes regarding these data in the spreadsheet.

* `1929.dat`: back to the source

  These are only those colors physically appearing in the 1929 Munsell Book of
  Color. These data might be of useful for those interested in the input colors
  used for the scaling experiments leading to the 1943 renotation. Remember
  though, these are renotation colors of those original patches, not necessarily
  the colors of the input data used in the visual experiment.


## Munsell lightness data

From

Final Report of the O.S.A. Subcommittee on the Spacing of the Munsell Colors,
SIDNEY M. NEWHALL, DOROTHY NICKERSON, DEANE B. JUDD,
JOURNAL OF THE OPTICAL SOCIETY OF AMERICA,
VOLUME 33, NUMBER 7, JULY, 1943.

TABLE II. I.C.I. (Y) equivalents (in percent relative to MgO) of the recommended Munsell
value scale (V) from 0/ to 10/.
