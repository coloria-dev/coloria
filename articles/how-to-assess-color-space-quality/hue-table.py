import numpy

import colorio

color_spaces = [
    colorio.CAM02("UCS", 0.69, 20, 64 / numpy.pi / 5),
    colorio.CAM16UCS(0.69, 20, 64 / numpy.pi / 5),
    colorio.CIELAB(),
    colorio.CIELUV(),
    colorio.IPT(),
    colorio.JzAzBz(),
    colorio.OKLAB(),
    colorio.OsaUcs(),
    colorio.RLAB(),
    colorio.XYY(),
]

for cs in color_spaces:
    vals = [
        sum(colorio.data.hung_berns.residuals(cs)),
        sum(colorio.data.ebner_fairchild.residuals(cs)),
        sum(colorio.data.xiao.residuals(cs)),
    ]
    s = sum(vals)
    print(f"{cs.name} & {vals[0]:.3f} & {vals[1]:.3f} & {vals[2]:.3f} & {s:.3f}\\\\")
