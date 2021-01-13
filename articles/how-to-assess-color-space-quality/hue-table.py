import numpy

import colorio

color_spaces = [
    colorio.cs.CAM02("UCS", 0.69, 20, 64 / numpy.pi / 5),
    colorio.cs.CAM16UCS(0.69, 20, 64 / numpy.pi / 5),
    colorio.cs.CIELAB(),
    colorio.cs.CIELUV(),
    colorio.cs.IPT(),
    colorio.cs.ICtCp(),
    colorio.cs.JzAzBz(),
    colorio.cs.OKLAB(),
    colorio.cs.OsaUcs(),
    colorio.cs.RLAB(),
    colorio.cs.XYY1(),
]

for cs in color_spaces:
    vals = [
        colorio.data.hung_berns.stress(cs),
        colorio.data.ebner_fairchild.stress(cs),
        colorio.data.xiao.stress(cs),
    ]
    avg = numpy.average(numpy.concatenate(vals))
    vals = [numpy.average(val) for val in vals]
    # print(f"{cs.name}    {sum(vals)}")
    print(f"{cs.name} & {vals[0]:.1f} & {vals[1]:.1f} & {vals[2]:.1f} & {avg:.1f}\\\\")
