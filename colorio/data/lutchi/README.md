# LUTCHI data set

Retrieved from
https://web.archive.org/web/20031230164218/http://colour.derby.ac.uk:80/colour/info/lutchi/data/.

Also available at http://www.rit-mcsl.org/fairchild/CAM.html.


## Using the LUTCHI Colour Appearance Data

M. Ronnier Luo and Peter A. Rhodes

Design Research Centre
University of Derby
Derby UK

### Introduction
The LUTCHI data is a large set of psychophysical experimental data for
describing colour appearance. The main body of the results was obtained from
two consecutive research projects funded by the British Government: Alvey
(1987-1989) and IEAPT (1990-1992) programmes.  These results have previously
been published in Color Research Application1-4. The project consortium were:
LUTCHI, Crosfield Electronics (now Fuji Film Electronic Imaging (FFEI)), Coats
Viyella and Sigmax Displays Ltd (now MISYS). The data were produced at the
Loughborough University of Technology Computer and Human Interface Research
Centre and are hence named the LUTCHI colour appearance data.  Subsequently,
two new data sets were also accumulated: Kuo & Luo5, and BIT6. These data sets
were also included to form the full LUTCHI Colour Appearance Data Set.

The data were used to test various colour appearance models2, to refine the
Hunt colour appearance model7,8 and to derive the LLAB colour appearance
model9,10.  Most importantly, this data set  together with the others were used
for the development of the CIE colour appearance model, CIECAM97S11,12.

The experiments were conducted using a magnitude estimation method1 under
various conditions of illumination. The data are divided into 8 groups
according to the experimental viewing conditions summarised in Table I. They
are mainly categorised by the media used: reflective (R), self-luminous
(Monitor), 35mm projection (35mm), large cut-sheet transparency (LT) and BIT.
The latter group was complied by the Beijing Institute of Technology, China to
investigate colour appearance under mesopic region with an isolated viewing
field. It was treated as an independent group. For data in the reflective
group, they are further divided into 4 subgroups: high luminance (HL), low
luminance (LL), various luminance (VL) and textile samples (textile).  In
total, there are 59 different phases. Each phase represents a set of distinct
experimental conditions. A number of colours were assessed by a panel of 4-8
normal colour vision observers in each phase. Each colour was scaled in terms
of its lightness, colourfulness and hue attributes.  In addition, the
brightness attribute was used in the R-VL group. (There are 12 phases in the R-
VL group. The lightness, colourfulness and hue, and brightness, colourfulness
and hue results were obtained for the first and last 6 phases respectively.)

#### Table I Summary of the experimental groups

| Data Group    | No.of phases  | Description of each data group                                | Ref. No. |
| ------------- | ------------- | ------------------------------------------------------------- | -------- |
| R-HL          | 6             | Reflection media with luminances ranging 364-232 cd/m2        | 1        |
| R-LL          | 6             | Reflection media with luminances ranging 44-41 cd/m2          | 1        |
| R-VL          | 12            | Reflection media with luminances ranging 843-0.4 cd/m2        | 2        |
| R-textile     | 3             | Large textile samples with luminances ranging 730-340 cd/m2   | 4        |
| Self-luminous | 11            | Monitor colours with luminances ranging 45-20 cd/m2           | 1        |
| 35mm          | 6             | 35 mm transparency with luminances ranging 113-45 cd/m2       | 3        |
| LT            | 10            | cut-sheet transparency with luminances ranging 2100-320 cd/m2 | 3        |
| BIT           | 5             | Isolated viewing field with luminances ranging 90-3.6 cd/m2   | 6        |

### Data File Description
Table IIa summarises the data files and viewing conditions used in each
experimental phase.  Each phase includes 2 data files: visual and colorimetric.
The data in the visual file have been  arranged into 4 columns: the sample
number, mean visual results of the lightness or brightness, colourfulness and
hue.  The arithmetic mean was used for calculating lightness (0 for black, 100
for white) and hue (0 for red, 100 for yellow, 200 for green, 300 for blue and
400 for red) results, the geometric mean for brightness (0 for black) and
colourfulness (0 for neutral colours) results with open ended scales. The data
in the colorimetric file have been arranged into 3 columns: x, y and Y, which
were measured using a telespectroradiometer (TSR).  The sequence of
experimental samples is the same for both the visual and colorimetric files.
The other parameters are: the number of samples used in the phase, the start
and end numbers of the neutral samples, the scaling factor of the Y values used
for adjusting those in the colorimetric file, the Y value of the background and
the luminance of the reference white sample. The illuminant used in each phase
is given in Table IIb described using the Xo Yo and Zo tristimulus values.

The visual results for the neutral samples have no hue and a colourfulness of
zero. In the visual data file, they are not consistently expressed. Users need
to realise this before attempting to use these samples.  The hue results of the
neutral colours should not be used to test or derive colour appearance models.
The Y value of the background has already been scaled. Luminance, in cd/m2
units, was measured against the reference white by the TSR under a particular
phase. All colourfulness results are in the same visual scale except those of
the R-textile and BIT data. For testing colour models' performance, there is a
need to obtain a suitable scaling factor to adjust each model's chroma or
colourfulness predictions to make them have the same scale as the visual
results.

#### Table IIa: Summary of the data files and viewing conditions used in each phase.

| Group | Phase | Filename Visual | Filename Colorimetric | No. of Samples | Neutrals Start | Neutrals End | Scaling Y of Factor | Scaling Y of Background | Luminance |
| --------- | -- | ------------ | ------------ | --- | -- | -- | ---- | ----- | ----- |
| R-HL      | 1  | nlmean.wh    | cold50wnl    | 105 | 41 | 46 | 0.88 | 100.0 | 264.0 |
|           | 2  | nlmean.bh    | cold50gb     | 105 | 41 | 46 | 0.84 | 6.2   | 252.0 |
|           | 3  | nlmean.gh    | cold50gb     | 105 | 41 | 46 | 0.84 | 21.5  | 252.0 |
|           | 4  | nld65.gh     | cold65       | 105 | 41 | 46 | 0.81 | 21.5  | 243.0 |
|           | 5  | nlwf.gh      | colwf        | 105 | 41 | 46 | 0.84 | 21.5  | 252.0 |
|           | 6  | nla.gh       | colah        | 105 | 41 | 46 | 0.84 | 21.5  | 232.0 |
| R-LL      | 1  | nlmean.wl    | cold50wnl    | 105 | 41 | 46 | 0.88 | 100.0 | 44.0  |
|           | 2  | nlmean.bl    | cold50gb     | 105 | 41 | 46 | 0.84 | 6.2   | 42.0  |
|           | 3  | nlmean.gl    | cold50gb     | 105 | 41 | 46 | 0.84 | 21.5  | 42.0  |
|           | 4  | nld65.gl     | cold65       | 105 | 41 | 46 | 0.81 | 21.5  | 40.5  |
|           | 5  | nlwf.gl      | colwf        | 105 | 41 | 46 | 0.84 | 21.5  | 42.0  |
|           | 6  | nla.gl       | colal        | 105 | 41 | 46 | 0.84 | 21.4  | 42.0  |
| R-VL      | 1  | mean4.p1     | col.rf.p1    | 40  | 37 | 40 | 1.00 | 21.5  | 843.0 |
|           | 2  | mean4.p2     | col.rf.p2    | 40  | 37 | 40 | 1.00 | 21.5  | 200.0 |
|           | 3  | mean4.p3     | col.rf.p3    | 40  | 37 | 40 | 1.00 | 21.5  | 62.0  |
|           | 4  | mean4.p4     | col.rf.p4    | 40  | 37 | 40 | 1.00 | 21.5  | 17.0  |
|           | 5  | mean4.p5     | col.rf.p5    | 40  | 37 | 40 | 1.00 | 21.5  | 6.0   |
|           | 6  | mean4.p6     | col.rf.p6    | 40  | 37 | 40 | 1.00 | 21.5  | 0.4   |
|           | 7  | mean4.p7     | col.rf.p1    | 40  | 37 | 40 | 1.00 | 21.5  | 843.0 |
|           | 8  | mean4.p8     | col.rf.p2    | 40  | 37 | 40 | 1.00 | 21.5  | 200.0 |
|           | 9  | mean4.p9     | col.rf.p3    | 40  | 37 | 40 | 1.00 | 21.5  | 62.0  |
|           | 10 | mean4.p10    | col.rf.p4    | 40  | 37 | 40 | 1.00 | 21.5  | 17.0  |
|           | 11 | mean4.p11    | col.rf.p5    | 40  | 37 | 40 | 1.00 | 21.5  | 6.0   |
|           | 12 | mean4.p12    | col.rf.p6    | 40  | 37 | 40 | 1.00 | 21.5  | 0.4   |
| R-textile | 1  | kuo.d65.vis  | kuo.d65.col  | 240 | 1  | 12 | 0.74 | 16.0  | 250.0 |
|           | 2  | kuo.tl84.vis | kuo.tl84.col | 239 | 1  | 10 | 0.74 | 16.0  | 540.0 |
|           | 3  | kuo.a.vis    | kuo.a.col    | 239 | 1  | 11 | 0.74 | 16.0  | 250.0 |
| CRT       | 1  | lmean.ww     | cold50wl     | 94  | 39 | 44 | 0.89 | 100.0 | 40.0  |
|           | 2  | lmean.bb     | cold50gbl    | 100 | 39 | 44 | 0.89 | 5.0   | 44.5  |
|           | 3  | lmean.gg     | cold50gbl    | 100 | 39 | 44 | 0.89 | 20.0  | 44.5  |
|           | 4  | lmean.gw     | cold50gbl    | 100 | 39 | 44 | 0.89 | 20.0  | 44.5  |
|           | 5  | lmean.gb     | cold50gbl    | 100 | 39 | 44 | 0.89 | 20.0  | 44.5  |
|           | 6  | ld65.gg      | cold65.l     | 103 | 39 | 44 | 0.81 | 21.5  | 40.5  |
|           | 7  | ld65.gw      | cold65.l     | 103 | 39 | 44 | 0.81 | 21.5  | 40.5  |
|           | 8  | lwf.gg       | colwf.l      | 103 | 39 | 44 | 0.84 | 21.5  | 28.4  |
|           | 9  | lwf.gw       | colwf.l      | 103 | 39 | 44 | 0.84 | 21.5  | 28.4  |
|           | 10 | la.gg        | cola.l       | 61  | 29 | 34 | 0.84 | 21.5  | 20.3  |
|           | 11 | la.gw        | cola.l       | 61  | 29 | 34 | 0.84 | 21.5  | 20.3  |
| 35mm      | 1  | mean.35.p1   | col.35.p1    | 99  | 93 | 99 | 1.00 | 15.6  | 75.0  |
|           | 2  | mean.35.p2   | col.35.p2    | 99  | 93 | 98 | 1.00 | 14.7  | 75.0  |
|           | 3  | mean.35.p3   | col.35.p3    | 99  | 93 | 99 | 1.00 | 18.9  | 113.0 |
|           | 4  | mean.35.p4   | col.35.p1    | 99  | 93 | 99 | 1.00 | 18.9  | 45.0  |
|           | 5  | mean.35.p5   | col.35.p5    | 95  | 89 | 95 | 1.00 | 19.2  | 47.0  |
|           | 6  | mean.35.p6   | col.35.p6    | 36  | 26 | 30 | 1.00 | 18.9  | 113.0 |
| LT        | 1  | mean.p1      | col.p1       | 98  | 94 | 98 | 1.00 | 15.9  | 2259.0 |
|           | 2  | mean.p2      | col.p2       | 98  | 94 | 98 | 1.00 | 17.1  | 689.0 |
|           | 3  | mean.p3      | col.p3       | 98  | 94 | 98 | 1.00 | 16.9  | 325.0 |
|           | 4  | mean.p4      | col.p4       | 98  | 94 | 98 | 1.00 | 17.4  | 670.0 |
|           | 5  | mean.p5t     | col.p5t      | 97  | 93 | 97 | 1.00 | 9.6   | 1954.0|
|           | 6  | mean.p6t     | col.p6t      | 94  | 90 | 94 | 1.00 | 9.5   | 619.0  |
|           | 7  | mean.p7t     | col.p7t      | 93  | 89 | 93 | 1.00 | 9.8   | 319.0 |
|           | 8  | mean.p8      | col.p8       | 98  | 94 | 98 | 1.00 | 9.4   | 642.0 |
|           | 9  | mean.p10     | col.p10      | 98  | 94 | 98 | 1.00 | 9.6   | 658.0 |
|           | 10 | mean.p1t     | col.p1t      | 94  | 90 | 94 | 1.00 | 17.5  | 680.0 |
| BIT       | 1  | BIT.p1.vis   | BIT.p1.col   | 120 | 1  | 16 | 1.00 | 0.6   | 90.0  |
|           | 2  | BIT.p2.vis   | BIT.p2.col   | 120 | 1  | 18 | 1.00 | 0.6   | 3.6   |
|           | 3  | BIT.p3.vis   | BIT.p3.col   | 120 | 1  | 15 | 1.00 | 0.6   | 90.0  |
|           | 4  | BIT.p4.vis   | BIT.p4.col   | 90  | 1  | 18 | 1.00 | 0.6   | 90.0  |
|           | 5  | BIT.p5.vis   | BIT.p5.col   | 90  | 1  | 12 | 1.00 | 0.6   | 3.6   |

#### Table IIb Normalised tristimulus values (Xo Yo Zo) for each phase

| Group     | Phase | File Name | File Name    | No. of Visual Colorimetric  Samples | Xo | Yo | Zo |
| --------- | -- | ------------ | ------------ | --- | ------ | ------ | ----- |
| R-HL      | 1  | nlmean.wh    | cold50wnl    | 105 | 97.13  | 100.00 | 76.62 |
|           | 2  | nlmean.bh    | cold50gb     | 105 | 97.09  | 100.00 | 83.10 |
|           | 3  | nlmean.gh    | cold50gb     | 105 | 97.09  | 100.00 | 83.10 |
|           | 4  | nld65.gh     | cold65       | 105 | 94.52  | 100.00 | 114.98 |
|           | 5  | nlwf.gh      | colwf        | 105 | 102.50 | 100.00 | 47.93 |
|           | 6  | nla.gh       | colah        | 105 | 112.92 | 100.00 | 28.62 |
| R-LL      | 1  | nlmean.wl    | cold50wnl    | 105 | 97.13  | 100.00 | 76.62 |
|           | 2  | nlmean.bl    | cold50gb     | 105 | 97.09  | 100.00 | 83.10 |
|           | 3  | nlmean.gl    | cold50gb     | 105 | 97.09  | 100.00 | 83.10 |
|           | 4  | nld65.gl     | cold65       | 105 | 94.52  | 100.00 | 114.98 |
|           | 5  | nlwf.gl      | colwf        | 105 | 102.50 | 100.00 | 47.93 |
|           | 6  | nla.gl       | colal        | 105 | 117.26 | 100.00 | 22.44 |
| R-VL      | 1  | mean4.p1     | col.rf.p1    | 40  | 94.04  | 100.00 | 76.29 |
|           | 2  | mean4.p2     | col.rf.p2    | 40  | 93.04  | 100.00 | 72.24 |
|           | 3  | mean4.p3     | col.rf.p3    | 40  | 93.94  | 100.00 | 73.51 |
|           | 4  | mean4.p4     | col.rf.p4    | 40  | 93.44  | 100.00 | 72.49 |
|           | 5  | mean4.p5     | col.rf.p5    | 40  | 92.22  | 100.00 | 70.12 |
|           | 6  | mean4.p6     | col.rf.p6    | 40  | 90.56  | 100.00 | 58.59 |
|           | 7  | mean4.p7     | col.rf.p1    | 40  | 94.04  | 100.00 | 76.29 |
|           | 8  | mean4.p8     | col.rf.p2    | 40  | 93.04  | 100.00 | 72.24 |
|           | 9  | mean4.p9     | col.rf.p3    | 40  | 93.94  | 100.00 | 73.51 |
|           | 10 | mean4.p10    | col.rf.p4    | 40  | 93.44  | 100.00 | 72.49 |
|           | 11 | mean4.p11    | col.rf.p5    | 40  | 92.22  | 100.00 | 70.12 |
|           | 12 | mean4.p12    | col.rf.p6    | 40  | 90.56  | 100.00 | 58.59 |
| R-textile | 1  | kuo.d65.vis  | kuo.d65.col  | 240 | 96.46  | 100.00 | 108.62 |
|           | 2  | kuo.tl84.vis | kuo.tl84.col | 239 | 103.07 | 100.00 | 64.29 |
|           | 3  | kuo.a.vis    | kuo.a.col    | 239 | 115.19 | 100.00 | 23.75 |
| CRT       | 1  | lmean.ww     | cold50wl     | 94  | 97.13  | 100.00 | 76.62 |
|           | 2  | lmean.bb     | cold50gbl    | 100 | 97.13  | 100.00 | 76.62 |
|           | 3  | lmean.gg     | cold50gbl    | 100 | 97.13  | 100.00 | 76.62 |
|           | 4  | lmean.gw     | cold50gbl    | 100 | 97.13  | 100.00 | 76.62 |
|           | 5  | lmean.gb     | cold50gbl    | 100 | 97.13  | 100.00 | 76.62 |
|           | 6  | ld65.gg      | cold65.l     | 103 | 94.52  | 100.00 | 114.98 |
|           | 7  | ld65.gw      | cold65.l     | 103 | 94.52  | 100.00 | 114.98 |
|           | 8  | lwf.gg       | colwf.l      | 103 | 102.50 | 100.00 | 47.93 |
|           | 9  | lwf.gw       | colwf.l      | 103 | 102.50 | 100.00 | 47.93 |
|           | 10 | la.gg        | cola.l       | 61  | 117.26 | 100.00 | 22.44 |
|           | 11 | la.gw        | cola.l       | 61  | 117.26 | 100.00 | 22.44 |
| 35mm      | 1  | mean.35.p1   | col.35.p1    | 99  | 92.90  | 100.00 | 46.10 |
|           | 2  | mean.35.p2   | col.35.p2    | 99  | 86.71  | 100.00 | 75.49 |
|           | 3  | mean.35.p3   | col.35.p3    | 99  | 92.57  | 100.00 | 45.47 |
|           | 4  | mean.35.p4   | col.35.p1    | 99  | 92.90  | 100.00 | 46.10 |
|           | 5  | mean.35.p5   | col.35.p5    | 95  | 93.80  | 100.00 | 52.39 |
|           | 6  | mean.35.p6   | col.35.p6    | 36  | 95.32  | 100.00 | 53.37 |
| LT        | 1  | mean.p1      | col.p1       | 98  | 93.09  | 100.00 | 62.02 |
|           | 2  | mean.p2      | col.p2       | 98  | 93.23  | 100.00 | 61.15 |
|           | 3  | mean.p3      | col.p3       | 98  | 92.36  | 100.00 | 59.91 |
|           | 4  | mean.p4      | col.p4       | 98  | 93.05  | 100.00 | 58.58 |
|           | 5  | mean.p5t     | col.p5t      | 97  | 93.30  | 100.00 | 59.54 |
|           | 6  | mean.p6t     | col.p6t      | 94  | 93.20  | 100.00 | 60.48 |
|           | 7  | mean.p7t     | col.p7t      | 93  | 92.43  | 100.00 | 59.21 |
|           | 8  | mean.p8      | col.p8       | 98  | 93.34  | 100.00 | 57.98 |
|           | 9  | mean.p10     | col.p10      | 98  | 93.34  | 100.00 | 57.98 |
|           | 10 | mean.p1t     | col.p1t      | 94  | 93.34  | 100.00 | 57.98 |
| BIT       | 1  | BIT.p1.vis   | BIT.p1.col   | 120 | 100.60 | 100.00 | 113.20 |
|           | 2  | BIT.p2.vis   | BIT.p2.col   | 120 | 100.60 | 100.00 | 113.20 |
|           | 3  | BIT.p3.vis   | BIT.p3.col   | 120 | 100.60 | 100.00 | 113.20 |
|           | 4  | BIT.p4.vis   | BIT.p4.col   | 90  | 100.60 | 100.00 | 113.20 |
|           | 5  | BIT.p5.vis   | BIT.p5.col   | 90  | 100.60 | 100.00 | 113.20 |

An agreement was reached between members of the consortium to make the data
available in May 1997.  The data can be publicly accessed from the Internet web
site of the Design Research Centre (DRC), University of Derby
(http://ziggy.derby.ac.uk/colour/). Researchers or industrialists are welcome
to acquire this data set for further study.  We hope that more models will be
derived or refined with the availability of this data set.

### REFERENCES

1. M.R. Luo, A.A. Clarke, P.A.Rhodes, S.A.R. Scrivener, A. Schappo and C.J.
   Tait, Quantifying Colour Appearance. Part I. LUTCHI Colour Appearance Data',
   Color Res.  Appl., 16 166-180 (1991).

2. M.R. Luo, A.A. Clarke, P.A.Rhodes, S.A.R. Scrivener, A. Schappo and C.J.
   Tait, Quantifying Colour Appearance. Part II. Testing Colour Models
   Performance Using LUTCHI  Colour Appearance Data, Color Res. Appl., 16
   181-197 (1991).

3. M.R. Luo, X.W. Gao, P.A. Rhodes, J.H. Xin, A.A. Clarke and S.A.R. Scrivener,
   Quantifying Colour Appearance. Part III. Supplementary LUTCHI Colour
   Appearance Data, Color Res. Appl., 18 98-113 (1993).

4. M. R. Luo, X. W. Gao, P. A. Rhodes, H. J. Xin, A. A. Clarke, and S. A. R.
   Scrivener, Quantifying colour appearance: Part IV- Transmissive Media, Color
   Res. Appl., 18 191- 209 (1993).

5. W. G. Kuo, M.R. Luo and H. Bez, Various chromatic-adaptation transforms
   tested using new colour appearance data in textile, Color Res. Appl. 20
   313-327 (1995).

6. M. R. Luo, H. Xu, S-Q. Tang and F-K Zhou, Testing colour models performance
   using mesopic colour appearance data, Proceedings of the AIC '97 Kyoto,
   000-000 (1997).

7. R. W. G. Hunt, Revised Colour-Appearance Model for Related and Unrelated
   Colours, Color Res. Appl. 16,  146-165 (1991).

8. R.W.G. Hunt, An improved predictor of colourfulness in a model of colour
   vision, Color Res. Appl. 19 23-33 (1994).

9. M.R.Luo, M.C. Lo and W.G. Kuo, The LLAB(l:c) colour model, Color Res. Appl.,
   21 412-429 (1996).

10. J. Morovic and M.R. Luo, Two unsolved issues in colour management - colour
    appearance and gamut mapping, The 5th International Conference on High
    Technology on Imaging Science and Technology - Evolution and Promise,
    Chiba, Japan, 11-14 September, (1996).

11. M.R. Luo and R.W.G. Hunt, Testing colour appearance models using
    corresponding colour and magnitude estimation data sets, Proceedings of the
    CIE Expert Symposium '97 Colour Standards for Image Technology, CIE Publ.
    No. x000, Central Bureau of the CIE, Vienna, 000-000 (1997).

12. R.W.G. Hunt and M.R. Luo, The structure of the CIE 1997 Colour Appearance
    Model (CIECAM97S),   Proceedings of the CIE Expert Symposium '97 Colour
    Standards for Image Technology, CIE Publ. No. x000, Central Bureau of the
    CIE, Vienna, 000-000 (1997).
