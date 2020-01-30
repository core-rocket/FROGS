#!/bin/python

import sys
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
    ret = pm.ecef2geodetic(x, y, z)
    lat,lon,alt = ret
    print("ENU(%f, %f, %f) => LLH(%f, %f, %f)" % (e, n, u, lat, lon, alt))
    return ret

def csv2js(f):
    output = ""
    reader = csv.DictReader(f, skipinitialspace=True)
    data = {}
    for row in reader:
        wspeed = float(row["wspeed"])
        case = {"wdir": float(row["wdir"]), "ghp_x": float(row["ghp_x"]), "ghp_y": float(row["ghp_y"])}
        if not wspeed in data.keys():
            data[wspeed] = []
        data[wspeed].append(case)

    for wspeed in data.keys():
        cases = data[wspeed]
        output += "var ghp_%d = L.polygon([\n" % int(wspeed)
        for case in cases:
            #print(case)
            lat,lon,alt = enu2llh(case["ghp_x"], case["ghp_y"], 0.0)
            output += "\t[%f, %f],\n" % (lat, lon)
        output += "],{\n"
        output += "\tcolor: 'red',\n"
        output += "\tfillOpacity: 0.0\n"
        output += "}).addTo(map);\n\n"

    #print(output)
    return output

with open(sys.argv[1]) as f:
    js = csv2js(f)
    ghp_js = open("ghp-output.js", "w")
    ghp_js.write(js)
