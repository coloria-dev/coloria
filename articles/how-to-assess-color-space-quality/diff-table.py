import numpy

import colorio

color_spaces = [
    colorio.cs.CAM02("UCS", 0.69, 20, 64 / numpy.pi / 5),
    colorio.cs.CAM16UCS(0.69, 20, 64 / numpy.pi / 5),
    colorio.cs.CIELAB(),
    colorio.cs.CIELUV(),
    colorio.cs.IPT(),
    colorio.cs.JzAzBz(),
    colorio.cs.OKLAB(),
    colorio.cs.OsaUcs(),
    colorio.cs.RLAB(),
    colorio.cs.XYY1(),
]

for cs in color_spaces:
    vals = [
        colorio.data.macadam_1974.residual(cs),
    ]
    print(f"{cs.name} & {vals[0]:.3f}\\\\")
