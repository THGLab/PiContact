#!/usr/bin/python3

import string
from sys import argv,stdout
import os
from os import popen,system
from os.path import exists,basename
#from amino_acids import longer_names
from math import log, sqrt
from numpy import cross, array, dot, vdot, arccos
import numpy as np

import math


#import MakePiPlanes as pl
#from MakePiPlanes import agroups, rgroups

import AnnotatePiPlanesBASE as pl2

diCUT = float(4.9)
dCUT = "4.9"
CUTi = int(2)


def unit_vector(vector):
    #Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

  

tags_printed = {}

def print_out_atoms( ATOMS, points, match, T1, T2, TAG ):
    tags_printed[TAG] = 0
    Xp = points[0]
    Yp = points[1]

    Xm = match[0]
    Ym = match[1]
    
    closest = []

    hX = {}
    hY = {}
    for ix in Xp.keys():
        for iy in Yp.keys():
            xp = Xp[ix]
            yp = Yp[iy]

            dp = xp - yp
            dist = sqrt(dp[0]**2 + dp[1]**2 + dp[2]**2)##

            if dist < 1.5:
                closest.append( [dist, ix, iy] )
                if (ix in hX) == False:
                    hX[ix] = []
                if (iy in hY) == False:
                    hY[iy] = []
                hX[ix].append(iy)
                hY[iy].append(ix)

                
    for ix in hX.keys():
        dbX = Xm[ix]
        stringo = "%3s %24s %4s %8.3f %8.3f %8.3f %2s %8.3f %8.3f %8.3f\n" % (T1, TAG, ix, dbX[0], dbX[1], dbX[2], "to", Xp[ix][0], Xp[ix][1], Xp[ix][2])
        if (ix.split(':')[1] in ATOMS) == False:
        #    ATOMS[ix.split(':')[1]] = 0
        #ATOMS[ix.split(':')[1]] += 1
            ATOMS[ix] = 0
        ATOMS[ix] += 1
        
    for iy in hY.keys():
        dbY = Ym[iy]
        stringo = "%3s %24s %4s %8.3f %8.3f %8.3f %2s %8.3f %8.3f %8.3f\n" % (T2, TAG, iy, dbY[0], dbY[1], dbY[2], "to", Yp[iy][0], Yp[iy][1], Yp[iy][2])
        if (iy.split(':')[1] in ATOMS) == False:
        #    ATOMS[iy.split(':')[1]] = 0
        #ATOMS[iy.split(':')[1]] += 1
            ATOMS[iy] = 0
        ATOMS[iy] += 1



comp = {}
comp["ASN"] = ["CB", "CG", "OD1", "ND2"]
comp["GLN"] = ["CG", "CD", "OE1", "NE2"]
comp["ASP"] = ["CB", "CG", "OD1", "OD2"]
comp["GLU"] = ["CG", "CD", "OE1", "OE2"]
comp["ARG"] = ["NE", "CZ", "NH1", "NH2"]
comp["HIS"] = ["CB", "CG", "ND1", "CD2", "CE1", "NE2"]
comp["PHE"] = ["CB","CG","CD1","CD2","CE1","CE2","CZ"]
comp["TYR"] = ["CB","CG","CD1","CD2","CE1","CE2","CZ"]
comp["TRP"] = ["CB", "CG", "CD2", "CD1", "NE1", "CE2", "CZ2", "CH2", "CZ3", "CE3"]
comp["G"] = []
comp["A"] = []
comp["T"] = []
comp["U"] = []
comp["C"] = []

NACIDS = {"G":0,"U":0,"T":0,"C":0,"A":0}


#FNAME = "/home/rvernon/Dropbox/Projects2017/THG/MDRESULTS/DRIVE/5f66_amoebapro13water03_traj_nvt_10ns.pdb"
#OFILE = open("temp", 'w')

#FNAME = sys.argv[1]
#OFILE = open("nvm", 'w')

#FNAME = "/home/rvernon/Dropbox/Projects2017/THG/MDRESULTS/DRIVE/5f66_amoebapro13water03_traj_nvt_10ns.pdb"
#OFILE = open("Anno.1b", 'w')

#FNAME = "/home/rvernon/Dropbox/Projects2017/THG/MDRESULTS/DRIVE/5f66_ff99sbtip4pew_traj_nvt_10ns.pdb"
#OFILE = open("Anno.2b", 'w')

def run_on_readlines( modelnum, rlines, OLDCON ):

    plpdb = pl2.read_pdb( rlines )

    ATOMS = {}

    USEDPAIR = {}

    SC2SC = 0
    SC2BB = 0
    BB2BB = 0

    PICONTACTS = {}

    for x in range(len(plpdb.npCenArray)):

        CENX = plpdb.npCenArray[x]

        yindexlist = plpdb.planesKD.query_ball_point( CENX, r=6.0 )

        for y in yindexlist:

            xkey = plpdb.npKeyArray[x]
            ykey = plpdb.npKeyArray[y]

            if (xkey != ykey) and (xkey+"*"+ykey in USEDPAIR) == False:

                USEDPAIR[xkey+"*"+ykey] = 0
                USEDPAIR[ykey+"*"+xkey] = 0

                D1 = plpdb.planes[xkey][1]
                nvec1 = plpdb.planes[xkey][0][1]
                D2 = plpdb.planes[ykey][1]
                nvec2 = plpdb.planes[ykey][0][1]

                plcomp = pl2.compare_dict_to_dict( D1, nvec1, D2, nvec2 )
                if int(plcomp[1]) >= CUTi:
                    if abs(plcomp[0]) >= 0.8:

                        points = plcomp[2]
                        match_atoms = [ {}, {} ]
                        for i in points[0].keys():
                            match_atoms[0][i] = D1[i]

                        for i in points[1].keys():
                            match_atoms[1][i] = D2[i]

                        STYPE1 = "SC"
                        if len(xkey.split(':')[1]) == 6: STYPE1 = "BB"
                        if len(xkey.split(':')[1]) == 1: STYPE1 = "XNA"

                        STYPE2 = "SC"
                        if len(ykey.split(':')[1]) == 6: STYPE2 = "BB"
                        if len(ykey.split(':')[1]) == 1: STYPE2 = "XNA"

                        OTAG = xkey+"*"+ykey

                        if STYPE1 == "SC" and STYPE2 == "SC":
                            SC2SC += 1
                            
                        if (STYPE1 == "SC" and STYPE2 == "BB") or \
                           (STYPE1 == "BB" and STYPE2 == "SC"):
                            SC2BB += 1
                            
                        if STYPE1 == "BB" and STYPE2 == "BB":
                            BB2BB += 1

                        PICONTACTS[OTAG] = 0

                        print_out_atoms( ATOMS, points, match_atoms, STYPE1, STYPE2, OTAG )
    #OSTR = "Sums     ModelNum %8s   #Atoms %8i   #PiPi %8i   #SCSC %8i   #SCBB %8i   #BBBB %8i" % (modelnum, len(ATOMS.keys()), SC2SC+SC2BB+BB2BB, SC2SC, SC2BB, BB2BB)
    #print(OSTR)


    KIN = 0
    SCKIN = 0
    
    for K in OLDCON.keys():
        hasK = False
        if K in PICONTACTS.keys():
            hasK = True
        K2 = K.split('*')[1]+"*"+K.split('*')[0]
        if K2 in PICONTACTS.keys():
            hasK = True
        if hasK:
            KIN += 1
            if len(K.split('*')[0].split(':')[1]) == 3 and \
               len(K.split('*')[1].split(':')[1]) == 3:
                SCKIN += 1
            
    #print(KIN)

    if len(OLDCON.keys()) == 0:
        KP = 100.0
        KIN =  SC2SC+SC2BB+BB2BB
        SCKIN = SC2SC
    else:
        KP = 100.0*float(KIN)/float(len(OLDCON.keys()))

    OSTR = "Sums     ModelNum %8s   #Atoms %8i   %8i %8i %8.2f #PiPi %8i   #SCSC %8i   #SCBB %8i   #BBBB %8i" % (modelnum, len(ATOMS.keys()), SCKIN, KIN, KP, SC2SC+SC2BB+BB2BB, SC2SC, SC2BB, BB2BB)
    print(OSTR)
    
    OSTR = "Contacts ModelNum %8s    " % (modelnum)
    for P in PICONTACTS.keys():
        OSTR += " "+P
    print(OSTR)
    
    OSTR = "Atoms    ModelNum %8s    " % (modelnum)
    for A in ATOMS.keys():
        OSTR += " "+A
    print(OSTR)
    
    return PICONTACTS


uselines = open(argv[1]).readlines()

BASECON = run_on_readlines( "0", uselines, {} )

