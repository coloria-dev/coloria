'''
Helper tool for converting XLSX data to YAML.
'''
import openpyxl
import numpy
import yaml


wb = openpyxl.load_workbook('xiao.xlsx')

index = 3
filename = 'unique_blue.yaml'

# first worksheet
ws = wb.worksheets[index]
session1_x = [str(ws['C{}'.format(k)].value) for k in range(4, 1669)]
session1_y = [str(ws['D{}'.format(k)].value) for k in range(4, 1669)]
session1_z = [str(ws['E{}'.format(k)].value) for k in range(4, 1669)]

session2_x = [str(ws['G{}'.format(k)].value) for k in range(4, 1669)]
session2_y = [str(ws['H{}'.format(k)].value) for k in range(4, 1669)]
session2_z = [str(ws['I{}'.format(k)].value) for k in range(4, 1669)]

session3_x = [str(ws['K{}'.format(k)].value) for k in range(4, 1669)]
session3_y = [str(ws['L{}'.format(k)].value) for k in range(4, 1669)]
session3_z = [str(ws['M{}'.format(k)].value) for k in range(4, 1669)]

data = numpy.array([
    [session1_x, session2_x, session3_x],
    [session1_y, session2_y, session3_y],
    [session1_z, session2_z, session3_z],
    ]).T

data = data.reshape(-1, 9, 3, 3)
data = numpy.moveaxis(data, 1, 2)

with open(filename, 'w') as outfile:
    yaml.dump(data.tolist(), outfile)
