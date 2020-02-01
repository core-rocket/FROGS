#!/bin/python

import sys
import csv
import pymap3d as pm
import json

RESTRICT_ALTITUDE = 914.0

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

def read_ghp(f):
    reader = csv.DictReader(f, skipinitialspace=True)
    data = {}
    for row in reader:
        wspeed = float(row["wspeed"])
        case = {"wdir": float(row["wdir"]), "ghp_e": float(row["ghp_e"]), "ghp_n": float(row["ghp_n"]), "max_altitude": float(row["max_altitude"])}
        if not wspeed in data.keys():
            data[wspeed] = []
        data[wspeed].append(case)
    return data

def ghp2js(data):
    output = ""
    for wspeed in data.keys():
        cases = data[wspeed]
        output += "var ghp_%d = L.polygon([\n" % int(wspeed)
        for case in cases:
            #print(case)
            lat,lon,alt = enu2llh(case["ghp_e"], case["ghp_n"], 0.0)
            output += "\t[%f, %f],\n" % (lat, lon)
        output += "],{\n"
        output += "\tcolor: 'blue',\n"
        output += "\tfillOpacity: 0.0\n"
        output += "}).addTo(map);\n\n"

    #print(output)
    return output

def check_restrict(ghp):
    for wspeed in ghp.keys():
        cases = ghp[wspeed]
        print("check %f m/s..." % wspeed)
        for case in cases:
            check_case(case)

def check_case(case):
    case["restrict_altitude"] = (case["max_altitude"] < RESTRICT_ALTITUDE)

def show_restrict_table(ghp):
    print("| 風向風速 |", end="")
    for wspeed in ghp.keys():
        print(" %s m/s |" % str(wspeed), end="")
    print("")
    print("|:-|", end="")
    for ws in ghp.keys():
        print("-|", end="")
    print("")

    wdir_list = set([])
    for wspeed in ghp.keys():
        for case in ghp[wspeed]:
            wdir_list.add(case["wdir"])
    wdir_list = sorted(wdir_list)
    #print(wdir_list)
    for wdir in wdir_list:
        print("| %s deg |" % str(wdir), end = "")
        for wspeed in ghp.keys():
            for case in ghp[wspeed]:
                if wdir == case["wdir"]:
                    output = "不可("
                    if case["restrict_altitude"] == False:
                        output += "高度)"
                    else:
                        output = ""
                    output += " |"
                    print(output, end="")
        print("")

with open(sys.argv[1]) as f:
    ghp = read_ghp(f)
    check_restrict(ghp)

    js = ghp2js(ghp)
    ghp_js = open("ghp-output.js", "w")
    ghp_js.write(js)

    show_restrict_table(ghp)
