#!/usr/bin/python

import string
from sys import argv,stdout
from os import popen,system
from os.path import exists,basename
#from amino_acids import longer_names
from math import log, sqrt, pi
from numpy import cross, array, dot, vdot, arccos
import numpy as np
import copy


import scipy.spatial


rgroups = {}

#Sidechain pi-system planes are computed based on these atom lists
#The first atom is only used to orient the coordinate system, it is not considered part of the group
rgroups["ASN"] = ["CB",    "CG", "OD1", "ND2"]
rgroups["GLN"] = ["CG",    "CD", "OE1", "NE2"]
rgroups["ASP"] = ["CB",    "CG", "OD1", "OD2"]
rgroups["GLU"] = ["CG",    "CD", "OE1", "OE2"]
rgroups["ARG"] = ["NE",    "CZ",  "NH1", "NH2", "NE"] # This order is fairly important
rgroups["HIS"] = ["CB",    "CG", "ND1", "CD2", "CE1", "NE2"]
rgroups["PHE"] = ["CB",    "CG", "CD1", "CD2", "CE1", "CE2", "CZ"]
rgroups["TYR"] = ["CB",    "CG", "CD1", "CD2", "CE1", "CE2", "CZ"]
rgroups["TRP"] = ["CB",    "CG", "CD2", "CD1", "NE1", "CE2", "CZ2", "CH2", "CZ3", "CE3"]

rgroups["BB"]  = ["O",     "O",  "C",   "N"]

rgroups["U"] =   ["C1'",   "N1", "C6", "C2", "O2", "N3", "C4", "O4", "C5"]
rgroups["C"] =   ["C1'",   "N1", "C6", "C2", "O2", "N3", "C4", "N4", "C5"]
rgroups["A"] =   ["C1'",   "N9", "C8", "C4", "N3", "C2", "N1", "C6", "C5", "N7", "N6"]
rgroups["G"] =   ["C1'",   "N9", "C8", "C4", "N3", "C2", "N1", "C6", "C5", "N7", "N2", "O6"]

def unit_vector(vector):
    #Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)


def make_axis_vectors(A, B, C): # For np.array([x,y,z]) of connected atoms A->B->C
    
    Avec = A-B # vector along A->B bond
    
    r1cen = A - B#pdbdict[p1][2][rg1[0][0]] - pdbdict[p1][2][rg1[0][1]] 
    r2cen = C - B#pdbdict[p1][2][rg1[1][1]] - pdbdict[p1][2][rg1[0][1]]

    Nvec = cross(r1cen, r2cen) # vector parralel to plane

    r1cen = A - B
    r2cen = Nvec

    Cvec = cross(r1cen, r2cen) # 90 degree angle to Avec, in plane with A,B,C

    Alen = sqrt(Avec[0]**2+Avec[1]**2+Avec[2]**2)
    if Alen > 0.0:
        Avec /= Alen

    Nlen = sqrt(Nvec[0]**2+Nvec[1]**2+Nvec[2]**2)
    if Nlen > 0.0:
        Nvec /= Nlen

    Clen = sqrt(Cvec[0]**2+Cvec[1]**2+Cvec[2]**2)
    if Clen > 0.0:
        Cvec /= Clen

    return [Avec, Nvec, Cvec]

#   for both of two residues, gets the a list of points representing the 1.7 angstrom above the plane coordinates
#   closest to the center of mass for the other sidechain
def get_closepoints( udlist1, udlist2, center1, center2):

    l1 = {}
    for ik in udlist1.keys():
        i = udlist1[ik]
        dp = center2 - i[0]
        distanceU = sqrt(dp[0]**2 + dp[1]**2 + dp[2]**2)
        dp = center2 - i[1]
        distanceD = sqrt(dp[0]**2 + dp[1]**2 + dp[2]**2)

        if distanceU < distanceD:
            l1[ik] = i[0]
        else:
            l1[ik] = i[1]

    l2 = {}
    for ik in udlist2.keys():
        i = udlist2[ik]
        dp = center1 - i[0]
        distanceU = sqrt(dp[0]**2 + dp[1]**2 + dp[2]**2)
        dp = center1 - i[1]
        distanceD = sqrt(dp[0]**2 + dp[1]**2 + dp[2]**2)

        if distanceU < distanceD:
            l2[ik] = i[0]
        else:
            l2[ik] = i[1]

    return [l1,l2]


def compare_dict_to_dict( D1, nvec1, D2, nvec2 ):

#    vecs1 = make_axis_vectors( a1, D1[dk1[0]], D1[dk1[1]] )
#    vecs2 = make_axis_vectors( a2, D2[dk2[0]], D2[dk2[1]] )
#    nvec1 = vecs1[1]
#    nvec2 = vecs2[1]

    center1 = array([0.0,0.0,0.0])
    center2 = array([0.0,0.0,0.0])

    uddict1 = {}
    uddict2 = {}
    plen = 1.7
    for y in D1.keys():
        yi = D1[y]
        up = yi + (nvec1 * plen)
        down = yi - (nvec1 * plen)
        uddict1[y] = [up,down]
        center1 += yi
    center1 /= float(len(D1.keys()))

    for x in D2.keys():
        xi = D2[x]
        up = xi + (nvec2 * plen)
        down = xi - (nvec2 * plen)
        uddict2[x] = [up,down]
        center2 += xi
    center2 /= float(len(D2.keys()))

    dot_orthog = (vdot(nvec1, nvec2))
    pairmatch = calc_dict_match(uddict1,uddict2)
    closest_points = get_closepoints(uddict1, uddict2, center1, center2)

    return [dot_orthog, pairmatch, closest_points]

#   Finds the smallest number of atoms in two sidechains that satisfy the match criteria
#   where a match is when points 1.7 angstroms above the plane for each atom (ie: the carbon VDW surface) 
#   are within 1.5 angstroms of each other
def calc_dict_match( uddict1, uddict2 ):
    
    udlist1 = []
    for i in uddict1.keys():
        udlist1.append( uddict1[i] )

    udlist2 = []
    for i in uddict2.keys():
        udlist2.append( uddict2[i] )

    return calc_match( udlist1, udlist2 )


#   Finds the smallest number of atoms in two sidechains that satisfy the match criteria
#   where a match is when points 1.7 angstroms above the plane for each atom (ie: the carbon VDW surface) 
#   are within 1.5 angstroms of each other
def calc_match( udlist1, udlist2 ):
    
    return_points = False

    short = udlist1
    long = udlist2

    if len(udlist1) > len(udlist2):
        short = udlist2
        long = udlist1

    pointlist = []

    ntotal = 0.0
    nmatch = 0.0
    allmatch = []
    for n1 in range(len(long)):#ud1 in long:
        ud1 = long[n1]
        ntotal += 1.0

        vec1 = unit_vector(ud1[0] - ud1[1])

        c1 = (ud1[0]+ud1[1])
        c1[0] /= 2.0
        c1[1] /= 2.0
        c1[2] /= 2.0


        for n2 in range(len(short)):#short:
            min_dis = -1
            ud2 = short[n2]

            vec2 = unit_vector(ud2[0] - ud2[1])

            c2 = (ud2[0]+ud2[1])
            c2[0] /= 2.0
            c2[1] /= 2.0
            c2[2] /= 2.0

            #v1to2 = unit_vector(c2 - c1)
            #v2to1 = unit_vector(c1 - c2)

            #d1to2 = abs(vdot(v1to2, vec1))
            #d2to1 = abs(vdot(v2to1, vec2))

            #if d1to2 >= 0.5 and d2to1 >= 0.5:
            dp = ud1[0] - ud2[0]
            distance = sqrt(dp[0]**2 + dp[1]**2 + dp[2]**2)
            if distance < min_dis or min_dis == -1.0:
                min_dis = distance
                store_ud = [ud1[0],ud2[0]]

            dp = ud1[0] - ud2[1]
            distance = sqrt(dp[0]**2 + dp[1]**2 + dp[2]**2)
            if distance < min_dis:
                min_dis = distance
                store_ud = [ud1[0],ud2[1]]


            dp = ud1[1] - ud2[0]
            distance = sqrt(dp[0]**2 + dp[1]**2 + dp[2]**2)
            if distance < min_dis:
                min_dis = distance
                store_ud = [ud1[1],ud2[0]]


            dp = ud1[1] - ud2[1]
            distance = sqrt(dp[0]**2 + dp[1]**2 + dp[2]**2)
            if distance < min_dis:
                min_dis = distance
                store_ud = [ud1[1],ud2[1]]
            
            if min_dis < 1.5:#2.5:#1.6:#2.5:#1.7:#1.5:
                allmatch.append([min_dis, n1, n2, store_ud])
 
    allmatch.sort()
    allmatch.reverse()

    taken1 = {}
    taken2 = {}
    for a in allmatch:
        if (a[1] in taken1) == False and (a[2] in taken2) == False:
            pointlist.append(a[3])
        taken1[a[1]] = True
        taken2[a[2]] = True
            
    if len(taken1.keys()) < len(taken2.keys()):
        nmatch = len(taken1.keys())
    else:
        nmatch = len(taken2.keys())

    if return_points == False:
        return nmatch#100.0*(nmatch/ntotal)
    else:
        return pointlist # ha ha it's pointlist



class PlPDB(object):
    pdict = {}
    planes = {}
    atomnumbers = {}
    atomkeys = {}
    planeKD = ""
    npCenArray = np.array([])
    npKeyArray = np.array([])


#def make_pdb_map( pdb, frame1 ):


def write_pdb( pdb, MAP_XYZ2ATOM, xyzframe_AD, ofilename ):
    plpdb = PlPDB()
    
    plpdb.pdict = {}

    ofile = open(ofilename, 'w')

    #pdb = open( pfile ).readlines()
    for i in pdb:
        line = i.split()

        if len(line) > 5:
            if i[:4] == "ATOM":
                atomnumber = str(int(i[6:11]))
                atom = i[12:16].split()[0]
                res = i[17:20].split()[0]
                chain = i[21:22]
                #x = float(i[30:38])
                #y = float(i[38:46])
                #z = float(i[46:54])

                xyz = i[31:54]
                xyzATOM = MAP_XYZ2ATOM[xyz]
                x = xyzframe_AD[xyzATOM][0]
                y = xyzframe_AD[xyzATOM][1]
                z = xyzframe_AD[xyzATOM][2]

                ofile.write(i[:31]+"%7.3f %7.3f %7.3f\n" % (x, y, z))
                




WRONGRES = {"HID":"HIS", "HIE":"HIS", "CYX":"CYS"}
                
    
def read_pdb( pdb, MAP_XYZ2ATOM, xyzframe_AD ):
    plpdb = PlPDB()
    
    plpdb.pdict = {}
    

    #pdb = open( pfile ).readlines()
    for i in pdb:
        line = i.split()

        if len(line) > 5:
            if i[:4] == "ATOM":
                atomnumber = str(int(i[6:11]))
                atom = i[12:16].split()[0]
                res = i[17:20].split()[0]

                if res in WRONGRES:
                    res = WRONGRES[res]
                
                chain = i[21:22]
                #x = float(i[30:38])
                #y = float(i[38:46])
                #z = float(i[46:54])

                xyz = i[31:54]
                xyzATOM = MAP_XYZ2ATOM[xyz]
                x = xyzframe_AD[xyzATOM][0]
                y = xyzframe_AD[xyzATOM][1]
                z = xyzframe_AD[xyzATOM][2]

                resN = int(i[22:26])

                if chain == " ": chain = "A"
                
                key = chain+"."+str(resN)

                if (key in plpdb.pdict ) == False:
                    plpdb.pdict[key] = [res, {}]
                    plpdb.atomnumbers[key] = {}
                if (atom in plpdb.pdict[key][1] ) == False:
                    plpdb.pdict[key][1][atom] = array([x,y,z])
                    plpdb.atomnumbers[key][atom] = atomnumber
                    plpdb.atomkeys[atomnumber] = key+"."+res
                    
    plpdb.planes = {}
    CenArray = []
    KeyArray = []

    for key in plpdb.pdict.keys():
        if (plpdb.pdict[key][0] in rgroups):

            try:
                plane = {}

                rgroup = rgroups[plpdb.pdict[key][0]]
                
                vecs = make_axis_vectors( plpdb.pdict[key][1][rgroup[0]], \
                                          plpdb.pdict[key][1][rgroup[1]], \
                                          plpdb.pdict[key][1][rgroup[2]] )

                for atom in rgroups[plpdb.pdict[key][0]][1:]:
                    akey = atom+":"+plpdb.atomnumbers[key][atom]
                    #plane[atom] =  plpdb.pdict[key][1][atom]
                    plane[akey] =  plpdb.pdict[key][1][atom] 

                CEN =  np.array([0.0,0.0,0.0])
                for atom in plane.keys():
                    CEN += plane[atom]
                CEN /= float(len(plane.keys()))
                 
                planekey = key+":"+plpdb.pdict[key][0]

                CenArray.append(CEN)
                KeyArray.append(planekey)

                plpdb.planes[planekey] = [vecs, plane]

            except:
                pass

        resn = key.split('.')[0]+"."+str(int(key.split('.')[1])+1)
        if (resn in plpdb.pdict):

            try:

                vecs = make_axis_vectors( plpdb.pdict[key][1]["O"], \
                                          plpdb.pdict[key][1]["C"], \
                                          plpdb.pdict[resn][1]["N"] )
                plane = {}
                plane["O:"+plpdb.atomnumbers[key]['O']] = plpdb.pdict[key][1]["O"]
                plane["C:"+plpdb.atomnumbers[key]['C']] = plpdb.pdict[key][1]["C"]
                plane["N:"+plpdb.atomnumbers[resn]['N']] = plpdb.pdict[resn][1]["N"]
                
                CEN =  np.array([0.0,0.0,0.0])
                for atom in plane.keys():
                    CEN += plane[atom]
                CEN /= float(len(plane.keys()))

                planekey = key+":"+plpdb.pdict[key][0]+plpdb.pdict[resn][0]

                CenArray.append(CEN)
                KeyArray.append(planekey)

                plpdb.planes[planekey] = [vecs, plane]

            except:
                pass

    plpdb.npCenArray = np.asarray(CenArray)
    plpdb.npKeyArray = np.asarray(KeyArray)

    plpdb.planesKD = scipy.spatial.cKDTree( plpdb.npCenArray, leafsize=2 )

    return plpdb



def read_for_print( file ): # Fills dictionary with D{chain+"."res number} = [ aa, {atoms} = [x,y,z] ]
    pdb = open( file ).readlines()
    pdbdict = {}
    for i in pdb:
        line = i.split()

        if len(line) > 5:
            if i[:6] == "HETATM" or i[:4] == "ATOM":
                atom = i[12:16].split()[0]
                chain = i[21:22]
                resN = int(i[22:26])

                key = chain+"."+str(resN)

                if (key in pdbdict ) == False:
                    pdbdict[key] = {}

                if (atom in pdbdict[key]) == False:
                    pdbdict[key][atom] = i

    return pdbdict

""" 
plpdb = read_pdb("5f66.pdb")

for x in range(len(plpdb.npCenArray)):

    CENX = plpdb.npCenArray[x]

    yindexlist = plpdb.planesKD.query_ball_point( CENX, r=6.0 )

    for y in yindexlist:

        xkey = plpdb.npKeyArray[x]
        ykey = plpdb.npKeyArray[y]
        
        if xkey != ykey:
        
            D1 = plpdb.planes[xkey][1]
            nvec1 = plpdb.planes[xkey][0][1]
            D2 = plpdb.planes[ykey][1]
            nvec2 = plpdb.planes[ykey][0][1]
        
            print( compare_dict_to_dict( D1, nvec1, D2, nvec2 ) )
            exit()
"""
