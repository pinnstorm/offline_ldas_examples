#!/usr/bin/env python3
import os
import sys
import xarray as xr
import numpy as np
import metview as mv

args = sys.argv
ifile = args[1]
ofile = args[2]

grib = xr.open_dataset(ifile, engine='cfgrib')

lon = grib.variables['longitude'].data
lat = grib.variables['latitude'].data
scf = grib.variables['sd'].data
num = len(scf)

f = open('tmp.gpt', 'w')

f.write("#GEO\n")
f.write("#FORMAT XYV\n")
f.write("#DATA\n")

p_screen = np.where(
    (scf[:] >= 0.0) & (scf[:] <= 1.0)
)[0]

for x in p_screen:
    f.write(str("%%f  %%f  %%f\n" %% (lon[x], lat[x], scf[x])))
f.close()

grib = mv.read('tmp.grib')
geo  = mv.read('tmp.gpt')

num = mv.geo_to_grib(
    geopoints = geo,
    grid_definition_mode = "grib",
    template_grib = grib,
    parameter = 141,
    grib_table2_version = 128,
    interpolation_method = "NEAREST_GRIDPOINT_COUNT"
    )

sum = mv.geo_to_grib(
    geopoints = geo,
    grid_definition_mode = "grib",
    template_grib = grib,
    parameter = 141,
    grib_table2_version = 128,
    interpolation_method = "NEAREST_GRIDPOINT_SUM"
    )

num = mv.bitmap(num, 0)
scf = sum / num
scf = mv.nobitmap(scf, -1)

mv.write(ofile,scf)

os.remove('tmp.gpt')

