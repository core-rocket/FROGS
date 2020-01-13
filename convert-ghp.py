#!/bin/python

import csv
import pymap3d as pm
import json

launch_lat = 34.679730
launch_lon = 139.438373
launch_alt = 39.4422

e = 0.0
n = 0.0
u = 0.0

def enu2llh(e, n, u):
    x,y,z = pm.enu2ecef(e, n, u, launch_lat, launch_lon, launch_alt)
    return pm.ecef2geodetic(x, y, z)
    #print("LLH: (%f, %f, %f)" % (lat, lon, alt))

with open("ghp-3-70.csv") as f:
    ghp_js = open("ghp-output.js", "w")
    reader = csv.reader(f)
    data = [row for row in reader]
    for wind_vel in range(1,8):
        print("wind vel: %f" % wind_vel)
        ghp_js.write("var ghp%d = L.polygon([\n" % int(wind_vel))

        data_x = data[(wind_vel-1)*2]
        data_y = data[(wind_vel-1)*2+1]
        for i in range(16):
            x = float(data_x[i])
            y = float(data_y[i])
            lat,lon,alt = enu2llh(x, y, 0.0)
            print("\tpos: %f, %f -> %f, %f, %f" % (x, y, lat,lon,alt))
            ghp_js.write("\t[%f, %f],\n" % (lat, lon))
        ghp_js.write("],{\n")
        ghp_js.write("color: 'red',\n")
        ghp_js.write("fillOpacity: 0.0,\n")
        ghp_js.write("}).addTo(map);\n\n")
