#!/bin/python

import sys
import csv
import pymap3d as pm
import json
import simplekml

RM_RESTRICT_POINT = False

RESTRICT_ALTITUDE = 3000.0

launch_lat = 40.242865
launch_lon = 140.010450
launch_alt = 0.0

e = 0.0
n = 0.0
u = 0.0

max_altitude_case = {"max_altitude": 0.0}

def enu2llh(e, n, u):
    x,y,z = pm.enu2ecef(e, n, u, launch_lat, launch_lon, launch_alt)
    ret = pm.ecef2geodetic(x, y, z)
    lat,lon,alt = ret
    print("ENU(%f, %f, %f) => LLH(%f, %f, %f)" % (e, n, u, lat, lon, alt))
    return ret

def llh2enu(lat, lon, alt):
    x,y,z = pm.geodetic2ecef(lat, lon, alt)
    ret = pm.ecef2enu(x, y, z, launch_lat, launch_lon, launch_alt)
    e,n,u = ret
    print("LLH(%f, %f, %f) => ENU(%f, %f, %f)" % (lat,lon,alt, e,n,u))
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

def ghp2js(data, color):
    output = ""
    for wspeed in data.keys():
        cases = data[wspeed]
        output += "var ghp_%d = L.polygon([\n" % int(wspeed)
        for case in cases:
            if RM_RESTRICT_POINT and case["restrict_region"] == False:
                continue
            #print(case)
            lat,lon,alt = enu2llh(case["ghp_e"], case["ghp_n"], 0.0)
            output += "\t[%f, %f],\n" % (lat, lon)
        output += "],{\n"
        output += "\tcolor: '" + color + "',\n"
        output += "\tfillOpacity: 0.0\n"
        output += "}).addTo(map);\n\n"

    #print(output)
    return output

def check_restrict(ghp):
    for wspeed in ghp.keys():
        cases = ghp[wspeed]
        print("check %f m/s..." % wspeed)
        for case in cases:
            check_case(wspeed, case)

def check_region(e, n, u):
    ce,cn,cu = llh2enu(34.661857, 139.454987, 0.0)
    cr = 2500
    if (e-ce)*(e-ce) + (n-cn)*(n-cn) >= (cr*cr):
        return False
    p1e,p1n,p1u = llh2enu(34.684392, 139.454677, 0.0)
    p2e,p2n,p2u = llh2enu(34.667917, 139.428718, 0.0)

    if ((p2e-p1e)*(n-p1n) + (p2n-p1n)*(p1e-e)) * ((p2e-p1e)*(cn-p1n) + (p2n-p1n)*(p1e-ce)) > 0:
        return True
    return False

def check_case(wspeed, case):
    r_alt = (case["max_altitude"] < RESTRICT_ALTITUDE)
    r_reg = check_region(case["ghp_e"], case["ghp_n"], 0.0)
    case["restrict_altitude"] = r_alt
    #(case["max_altitude"] < RESTRICT_ALTITUDE)
    case["restrict_region"] = r_reg
    #check_region(case["ghp_e"], case["ghp_n"], 0.0)

    if r_alt and r_reg:
        global max_altitude_case
        if case["max_altitude"] > max_altitude_case["max_altitude"]:
            max_altitude_case = {"wspeed": wspeed}
            max_altitude_case.update(case)
        return True
    return False

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
                    elif case["restrict_region"] == False:
                        output += "領域)"
                    else:
                        output = ""
                    output += " |"
                    print(output, end="")
        print("")

def write_kml(ghp, color):
    kml = simplekml.Kml()
    for wspeed in ghp.keys():
        lname = str(wspeed) + "m/s"
        line = kml.newlinestring(name=lname)
        if color == "blue":
            c = simplekml.Color.blue
        elif color == "red":
            c = simplekml.Color.red
        line.style.linestyle.color = c
        line.style.linestyle.width = 8
        line.extrude = 1
        line.altitudemode = simplekml.AltitudeMode.absolute
        points = []
        for case in ghp[wspeed]:
            llh = enu2llh(case["ghp_e"], case["ghp_n"], 0.0)
            points.append((llh[1], llh[0], llh[2]))
        line.coords = points
    kml.save("output.kml")

with open(sys.argv[1]) as f:
    ghp = read_ghp(f)
    check_restrict(ghp)

    js = ghp2js(ghp, sys.argv[2])
    ghp_js = open("ghp-output.js", "w")
    ghp_js.write(js)

    show_restrict_table(ghp)

    write_kml(ghp, sys.argv[2])

    print("max altitude case: ", max_altitude_case)
