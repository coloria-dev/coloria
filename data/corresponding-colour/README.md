# Corresponding-colour data sets

Retrieved from
<https://web.archive.org/web/20031123133629/http://colour.derby.ac.uk:80/colour/info/catweb/>
and
<https://web.archive.org/web/20031212112125/http://colour.derby.ac.uk:80/colour/info/catweb/table2.html>.


## CORRESPONDING-COLOUR DATA SETS
M. Ronnier Luo and Peter A. Rhodes
Colour & Imaging Institute
University of Derby
Derby
England

### INTRODUCTION
A chromatic adaptation transform is capable of predicting corresponding
colours. Corresponding colours are described by two sets of tristimulus values
that give rise to the same perceived colour when the two samples are viewed
under test and reference light sources or illuminants. The two light sources or
illuminants differ in terms of their colour temperatures (or chromaticity
coodinates). A chromatic adaptation transform can be effectively used for
numerous industrial applications such as the evaluation of colour inconstancy
for surface samples, the calculation of colour difference between pairs of
samples assessed under non-daylight sources or illuminants, the provision of a
colour rendering index for assessing the quality of light sources, or the
prediction of coloured images across different sources or illuminants.

In October 1998, the CIE formed a new technical committee, TC 1-52, on
Chromatic Adaptation Transforms during its interim meeting in Baltimore, USA
with Professor M. R. Luo as its chairman. The objective of this committee is to
review certain chromatic adaptation transforms with a view to making a CIE
recommendation. The performance of chromatic adaptation transforms is normally
evaluated using corresponding-colour experimental data sets in which each
colour is defined by two sets of tristimulus values under two illuminants. Many
experiments were carried out using a variety of psychophysical methods under
different viewing conditions. A comprehensive collection of these data sets has
been accumulated by Luo and Hunt [1] for the purposes of deriving and
evaluating the CIE colour appearance model, CIECAM97s [2], and the CMC
chromatic adaptation transform, CMCCAT97 [3]. The Committee has decided to make
these data sets available via the Internet for public assessment. This task has
been completed and the resulting database is now available via the world wide
web at http://colour.derby.ac.uk. Researchers or industrialists are welcome to
acquire this database for further study. This paper gives a brief description
of each data set and describes the format of the data.

### EXPERIMENTAL DATA SETS
Fourteen data sets have been accumulated from nine sources [4-11]: the Color
Science Association of Japan (CSAJ), Helson, Lam and Rigg, LUTCHI, Kuo and Luo,
Breneman, Braun and Fairchild, and McCann. Each data set includes a number of
corresponding-colour pairs in which both colours in a pair appear the same when
each is viewed under different viewing conditions. Table I summarises the
experimental conditions in each data set including the number of phases (as
defined by a set of viewing conditions), the number of corresponding-colour
pairs and the viewing parameters used. The parameters considered are the light
sources used for the test and reference conditions, illuminance (lux), the
luminance factor of the neutral background (Y%), sample size, media and
psychophysical method.

The CSAJ [4] data was divided into three sets: -C, -Hunt and -Stevens according
to studies on chromatic adaptation, Hunt and Stevens effects respectively. The
Helson [5], Lam and Rigg [6] data sets include corresponding colours between
the test source (A) and reference source (D65). The LUTCHI [7] data includes
three sets - A, D50 and WF - which are the test illuminants against a reference
D65 simulator. Similarly, there are two sets for Kuo and Luo [8] data: A and
TL84, which are the test light sources against a reference D65 simulator. The
only data set based upon transparent media in this category is the Breneman [9]
data which was divided into two sets: -C and -L according to investigations on
chromatic adaptation and illuminance effects respectively. The Braun and
Fairchild [10] data was accumulated by asking observers to adjust monitor
colours to match those presented on reflection prints. The McCann [11] data
were obtained by investigating the chromatic adaptation effect using a Mondrain
figure viewed under highly chromatic test illuminants with low illuminances.
Its original data was further analysed to obtain corresponding tristimulus
values by Nayatani et al [12].

In total, 746 corresponding-colour pairs were gathered from experiments
involving 38 phases of viewing conditions. The psychophysical methods used are
haploscopic matching, memory matching and magnitude estimation.

### DATA FILE DESCRIPTION
Table II summarises the data file names and number of samples in each
experimental data set. For each file, the data were arranged in a fixed format.
In the top row, there are six figures corresponding to the tristimulus values
of the reference and test illuminants, i.e. Xr, Yr, Zr and Xt, Yt, Zt. There is
only one figure in the second row denoting the number of samples in the file.
The other rows are also arranged in the same manner as the first row. These are
the corresponding tristimulus values under the reference and test illuminants
for each sample.

#### Table I: Summary of the corresponding-colour data sets
Experimental conditions

| Data Set | No. of Phases | No. of Samples | Illuminant Test | Illuminant Ref | Illuminance (lux) | Background (Y%) | Sample Size | Medium | Experimental Method |
| ----------------- | - | --- | -------- | ------------- | ------- | -- | - | ----------- | ----------- |
| CSAJ-C            | 1 | 87  | D65      | A             | 1000    | 20 | S | Refl.       | Haploscopic |
| CSAJ-Hunt         | 4 | 20  | D65      | D65           | 10-3000 | 20 | S | Refl.       | Haploscopic |
| CSAJ-Stevens      | 4 | 19  | D65      | D65           | 10-3000 | 20 | S | Refl.       | Haploscopic |
| Helson            | 1 | 59  | D65      | A             | 1000    | 20 | S | Refl.       | Memory      |
| Lam & Rigg        | 1 | 58  | D65      | A             | 1000    | 20 | L | Refl.       | Memory      |
| Lutchi (A)        | 1 | 43  | D65      | A             | 1000    | 20 | S | Refl.       | Magnitude   |
| Lutchi (D50)      | 1 | 44  | D65      | D50           | 1000    | 20 | S | Refl.       | Magnitude   |
| Lutchi (WF)       | 1 | 41  | D65      | WF            | 1000    | 20 | S | Refl.       | Magnitude   |
| Kuo & Luo (A)     | 1 | 0   | D65      | A             | 1000    | 20 | L | Refl.       | Magnitude   |
| Kuo & Luo (TL84)  | 1 | 1   | D65      | TL84          | 1000    | 20 | S | Refl.       | Magnitude   |
| Breneman-C        | 9 | 107 | D65, D55 | A, P, G       | 50-3870 | 30 | S | Trans.      | Magnitude   |
| Breneman-L        | 3 | 36  | D55      | D55           | 50-3870 | 30 | S | Trans.      | Haploscopic |
| Braun & Fairchild | 4 | 66  | D65      | D30, D65, D95 | 129     | 20 | S | Mon., Refl. | Matching    |
| McCann            | 5 | 85  | D65      | R, Y, G, B    | 14-40   | 30 | S | Refl.       | Haploscopic |

#### Table II: The data file names and number of samples in each experimental data set

| Data Set          | No. of Specimen | File Names     |
| ----------------- | --------------- | -------------- |
| CSAJ-C            | 87              | CSAJ.da.dat    |
| CSAJ-Stevens      | 5               | Steve.10.dat   |
|                   | 5               | Steve.50.dat   |
|                   | 5               | Steve.1000.dat |
|                   | 4               | Steve.3000.dat |
| CSAJ-Hunt         | 5               | CSAJ.10.dat    |
|                   | 5               | CSAJ.50.dat    |
|                   | 5               | CSAJ.1000.dat  |
|                   | 5               | CSAJ.3000.dat  |
| Helson            | 59              | helson.ca.dat  |
| Lam & Rigg        | 58              | lam.da.dat     |
| Lutchi (A)        | 43              | lutchi.da.dat  |
| Lutchi (D50)      | 44              | lutchi.dd.dat  |
| Lutchi (WF)       | 41              | lutchi.dw.dat  |
| Kuo & Luo (A)     | 40              | Kuo.da.dat     |
| Kuo & Luo (TL84)  | 41              | Kuo.dt.dat     |
| Breneman-C        | 12              | Brene.p1.dat   |
|                   | 12              | Brene.p2.dat   |
|                   | 12              | Brene.p3.dat   |
|                   | 12              | Brene.p4.dat   |
|                   | 11              | Brene.p6.dat   |
|                   | 12              | Brene.p8.dat   |
|                   | 12              | Brene.p9.dat   |
|                   | 12              | Brene.p11.dat  |
|                   | 12              | Brene.p12.dat  |
| Breneman-L        | 12              | Brene.p5.dat   |
|                   | 12              | Brene.p7.dat   |
|                   | 12              | Brene.p10.dat  |
| Braun & Fairchild | 17              | RIT.1.dat      |
|                   | 16              | RIT.2.dat      |
|                   | 17              | RIT.3.dat      |
|                   | 16              | RIT.4.dat      |
| McCann            | 17              | mcan.b.dat     |
|                   | 17              | mcan.g.dat     |
|                   | 17              | mcan.grey.dat  |
|                   | 17              | mcan.r.dat     |
|                   | 17              | mcan.y.dat     |


### SUMMARY
A data library has been established to include fourteen corresponding-colour
data sets and is available at from http://colour.derby.ac.uk. It is publicly
available for all interested parties who are doing research in the areas of
colour appearance and chromatic adaptation. Users are encouraged to report new
findings to the Chairman of CIE TC 1-52.

### REFERENCES
1. LUO M.R. and HUNT R.W.G, Testing colour appearance models using
   corresponding-colour and magnitude-estimation data sets, Color Res. Appl. 23
   147-153 (1998).
2. LUO M.R. and HUNT R.W.G, The structure of the CIE 1997 colour appearance
   model (CIECAM97), Color Res. Appl. 23 138-146 (1998).
3. LUO M.R. and HUNT R.W.G, A chromatic adaptation transform and a colour
   inconstancy index, Color Res. Appl. 23 154-158 (1998).
4. MORI L., SOBAGAKI H., KOMATSUBARA H. and IKEDA K., Field trials on the CIE
   chromatic adaptation formula, Proceedings of the CIE 22nd Session, 55-58
   (1991).
5. HELSON H., JUDD D. B., and WARREN M. H., Object-color changes from daylight
   to incandescent filament illumination, Illum. Eng. 47, 221-233 (1952).
6. LAM K. M., Metamerism and colour constancy, Ph.D. Thesis, University of
   Bradford, 1985.
7. LUO M. R., CLARKE A. A., RHODES P. A., SCRIVENER, S. A. R., SCHAPPO A. and
   TAIT C.J., Quantifying Colour Appearance. Part I. LUTCHI Colour Appearance
   Data, Color Res. Appl., 16 166-180 (1991).
8. KUO W.G., LUO M.R., BEZ H. E., Various chromatic-adaptation transforms
   tested using new colour appearance data in textiles, Color Res. Appl., 21
   313-327 (1995).
9. BRENEMAN E. J., Corresponding chromaticities for different states of
   adaptation to complex visual fields, J. Opt. Soc. Am. 4:6 1115-1129 (1987).
10. MCCANN J.J., MCKEE S. P. and TAYLOR T. H., Quantitative studies in retinex
    theory, Vision Res. 16 445-458 (1976).
11. BRAUN K.M. and FAIRCHILD M.D., Psychophysical generation of matching images
    for cross-media color reproduction, in IS&T and SID's 4th Color Imaging
    Conference: Color Science, Systems and Applications, 214-220, IS&T,
    Springfield, Va., (1996).
12. NAYATANI Y., TAKAHAMA K. and SOBAGAKI H., Prediction of color appearance of
    object colors in a complex visual field, J. Light & Vis. Env. 19 5-14
    (1995).
