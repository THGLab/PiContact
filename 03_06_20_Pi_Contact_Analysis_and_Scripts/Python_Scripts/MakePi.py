#!/usr/bin/python

import string
from sys import argv,stdout
import os
from os import popen,system
from os.path import exists,basename
#from amino_acids import longer_names
from math import log, sqrt
from numpy import cross, array, dot, vdot, arccos
import numpy as np
from pymol.cgo import *
from pymol import cmd
import math


import MakePiPlanes as pl
from MakePiPlanes import agroups, rgroups


diCUT = float(4.9)
dCUT = "4.9"
CUTi = int(2)


def unit_vector(vector):
    #Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)


def rough_distance( res1, res2 ):
    dp = res1[1][pl.rgroups[res1[0]][0]] - res2[1][pl.rgroups[res2[0]][0]]
    distance = sqrt(dp[0]**2 + dp[1]**2 + dp[2]**2)
            
    if distance < 30.0 and distance > 0.5:
        return True
    else:
        return False    

def mindist_count( l1, l2):
    
    xdone = {}
    ydone = {}
    for x in l1:

        for y in l2:
            dp = x - y
            distance = sqrt(dp[0]**2 + dp[1]**2 + dp[2]**2)
            
            if distance <= diCUT:
                xdone[str(x)] = 0
                ydone[str(y)] = 0

    mincount = len(xdone.keys())
    if len(ydone.keys()) < mincount:
        mincount = len(ydone.keys())

    return mincount

def get_bblist( bb, bbn):
    bblist = []
    bblist.append(bb[1]["C"])
    bblist.append(bb[1]["O"])
    bblist.append(bbn[1]["N"])
    return bblist

def get_sclist( sc ):
    sclist = []
    ag = pl.agroups[sc[0]]
    for a in ag:
        sclist.append(sc[1][a])
    return sclist

def Nallcompare_sc2sc( GC, sc1, sc2, TAG, sc2type ):
    hits = [0,0,0]

    distance = mindist_count( get_sclist(sc1), get_sclist(sc2) )
    if distance >= CUTi:
        plcomp = pl.NEWcompare_res_to_res( sc1, sc2 )

        if int(plcomp[1]) >= CUTi:

            if abs(plcomp[0]) >= 0.8:

                points = plcomp[2]

                match_atoms = [ {}, {} ]
                for i in points[0].keys():
                    match_atoms[0][i] = sc1[1][i]

                for i in points[1].keys():
                    match_atoms[1][i] = sc2[1][i]

                if len(sc1[0]) == 3:
                    sc1type = "SC"
                else:
                    sc1type = "XNA"

                if len(sc2[0]) == 3:
                    sc2type = "SC"
                else:
                    sc2type = "XNA"

                print_out_prepfile( GC, points, match_atoms, sc1type, sc2type, TAG )

def Nallcompare_bb2sc( GC, bb, bbn, sc, TAG, sctype ):
    hits = [0,0,0]

    distance = mindist_count( get_bblist(bb, bbn), get_sclist(sc) )
    if distance >= CUTi:
        plcomp = pl.NEWcompare_res_to_BB( bb, bbn, sc)

        if int(plcomp[1]) >= CUTi:

            if abs(plcomp[0]) >= 0.8:

                points = plcomp[2]

                match_atoms = [ {}, {} ]
                for i in points[0].keys():
                    match_atoms[0][i] = sc[1][i]

                for i in points[1].keys():
                    if i != "N":
                        match_atoms[1][i] = bb[1][i]
                    else:
                        match_atoms[1][i] = bbn[1][i]

                if len(sc[0]) == 3:
                    sctype = "SC"
                else:
                    sctype = "XNA"

                print_out_prepfile( GC, points, match_atoms, sctype, "BB", TAG )

def Nallcompare_bb2bb( GC, bb1, bbn1, bb2, bbn2, TAG ):
    hits = [0,0,0]
    
    distance = mindist_count( get_bblist(bb1, bbn1), get_bblist(bb2, bbn2) )

   
    if distance >= CUTi:
        plcomp = pl.NEWcompare_BB_to_BB( bb1, bbn1, bb2, bbn2)

        #print(TAG, distance, plcomp[0], plcomp[1])
        
        if int(plcomp[1]) >= CUTi:

            if abs(plcomp[0]) >= 0.8:

                points = plcomp[2]

                match_atoms = [ {}, {} ]
                for i in points[0].keys():
                    if i != "N":
                        match_atoms[0][i] = bb1[1][i]
                    else:
                        match_atoms[0][i] = bbn1[1][i]
                

                for i in points[1].keys():
                    if i != "N":
                        match_atoms[1][i] = bb2[1][i]
                    else:
                        match_atoms[1][i] = bbn2[1][i]

                print_out_prepfile( GC, points, match_atoms, "BB", "BB", TAG )

tags_printed = {}

def print_out_prepfile( GC, points, match, T1, T2, TAG ):
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

    closest.sort()
    usedX = {}
    usedY = {}
    for ic in closest:
        uX = ic[1]
        uY = ic[2]
        if ((uX+uY) in usedX) == False:
            usedX[uX+uY] = 0
            usedX[uY+uX] = 0
            stringo = "%3s %24s %4s %8.3f %8.3f %8.3f %2s %8.3f %8.3f %8.3f\n" % ("C", TAG, "-", Xp[uX][0], Xp[uX][1], Xp[uX][2], "to", Yp[uY][0], Yp[uY][1], Yp[uY][2])
            #ofile.write(stringo)
            GC.append(stringo)
            
    for ix in hX.keys():
        dbX = Xm[ix]
        stringo = "%3s %24s %4s %8.3f %8.3f %8.3f %2s %8.3f %8.3f %8.3f\n" % (T1, TAG, ix, dbX[0], dbX[1], dbX[2], "to", Xp[ix][0], Xp[ix][1], Xp[ix][2])
        #ofile.write(stringo)
        GC.append(stringo)
    for iy in hY.keys():
        dbY = Ym[iy]
        stringo = "%3s %24s %4s %8.3f %8.3f %8.3f %2s %8.3f %8.3f %8.3f\n" % (T2, TAG, iy, dbY[0], dbY[1], dbY[2], "to", Yp[iy][0], Yp[iy][1], Yp[iy][2])
        #ofile.write(stringo)
        GC.append(stringo)


def run( FNAME, tags=False ):

    GC = []
    
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

    if exists(FNAME) == False:

        if exists(FNAME+".pdb"):
            FNAME = FNAME+".pdb"
        else:
            print("ERROR: File does not exist - "+FNAME)
            exit()
    
    pdbdict, hasbase = pl.read_pdb( FNAME )
    pl.is_complete(pdbdict)
    forprint = pl.read_for_print( FNAME )

    tags_used = {}

    pdbkeys = list(pdbdict.keys())
    
    for xn in range(len(pdbkeys)):
       x = pdbkeys[xn]
       aa1 = pdbdict[x][0]
       res1 = int(x.split('.')[1])
       chain1 = x.split('.')[0]

       res1_next = int(x.split('.')[1])+1
       xnext = chain1+'.'+str(res1_next)

       if pdbdict[x][2]:
           for yn in range(len(pdbkeys)):
               y = pdbkeys[yn]
               aa2 = pdbdict[y][0]
               res2 = int(y.split('.')[1])
               chain2 = y.split('.')[0]

               res2_next = int(y.split('.')[1])+1
               ynext = chain2+'.'+str(res2_next)

               if pdbdict[y][2]:
                   if rough_distance(pdbdict[x], pdbdict[y]):

                       if 1==1:#len(aa2) == 3:

                           if (aa1 in comp) and (aa2 in comp):
                               TAG = x+"."+aa1+":"+y+"."+aa2
                               TAGalt = y+"."+aa2+":"+x+"."+aa1
                               if (TAG in tags_used) == False and (TAGalt in tags_used) == False:
                                   tags_used[TAG] = 0
                                   Nallcompare_sc2sc( GC, pdbdict[x], pdbdict[y], TAG, "SC" )

                           if (aa1 in comp) and (ynext in pdbdict) and len(aa2) == 3:
                               if pdbdict[ynext][2]:
                                   TAG =    x+"."+aa1+":"+y+".BB"
                                   TAGalt = y+".BB:"+x+"."+aa1
                                   if (TAG in tags_used) == False and (TAGalt in tags_used) == False:
                                       tags_used[TAG] = 0
                                       Nallcompare_bb2sc( GC, pdbdict[y], pdbdict[ynext], pdbdict[x], TAG, "SC" )


                           if (ynext in pdbdict) and (xnext in pdbdict) and len(aa1) == 3 and len(aa2) == 3:
                               if pdbdict[ynext][2] and pdbdict[xnext][2]:
                                   TAG =    x+".BB:"+y+".BB"
                                   TAGalt = y+".BB:"+x+".BB"
                                   if (TAG in tags_used) == False and (TAGalt in tags_used) == False:
                                       tags_used[TAG] = 0
                                       Nallcompare_bb2bb( GC, pdbdict[x], pdbdict[xnext], pdbdict[y], pdbdict[ynext], TAG )

                       else:
                           if (aa1 in comp) and (aa2 in NACIDS):
                               TAG = x+"."+aa1+":"+y+"."+aa2
                               TAGalt = y+"."+aa2+":"+x+"."+aa1
                               if (TAG in tags_used) == False and (TAGalt in tags_used) == False:
                                   tags_used[TAG] = 0
                                   Nallcompare_sc2sc( GC, pdbdict[x], pdbdict[y], TAG, "XNA" )

                           if (xnext in pdbdict):
                               if pdbdict[xnext][2]:
                                   TAG = x+".BB"+":"+y+"."+aa2
                                   if (TAG in tags_used) == False:
                                       tags_used[TAG] = 0
                                       Nallcompare_bb2sc( GC, pdbdict[x], pdbdict[xnext], pdbdict[y], TAG, "XNA")




    bb = {"C":0,"O":0,"N":0}

    spatout = {}


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


    pfile = GC#open('sp2_contacts.prep').readlines()


    #colordictionary
    color = {}
    color["purple"] = [0.29, 0.0, 0.51,     0.8, 0.0, 0.8]
    color["cyan"] = [0.1, 0.3, 0.3,     0.0, 0.7, 0.8]
    color["orange"] = [0.6, 0.2, 0.0,     1.0, 0.5, 0.0]
    color["gray"] = [0.3, 0.3, 0.3,     0.3, 0.3, 0.3]

    radius = 0.14

    obj = {}
    for p in pfile:
       line = p.split()
       if len(line) == 10:

          TAG = "planes_"+FNAME.replace('.pdb', '')
          if tags:
              TAG = line[1].replace(':','_')
              
          if (TAG in obj) == False:
              obj[TAG] = []

          x1 = float(line[3])
          y1 = float(line[4])
          z1 = float(line[5])

          x2 = float(line[7])
          y2 = float(line[8])
          z2 = float(line[9])

          ctype = line[0]
          c = color["purple"]
          if ctype == "BB":
              c = color["cyan"]
          if ctype == "HET":
              c = color["orange"]
          if ctype == "C":
             c = color["gray"]
          if ctype == "XNA":
             c = color["gray"]

          if ctype != "C":
             obj[TAG].extend( [ CYLINDER, x1, y1, z1, x2, y2, z2, radius, c[0], c[1], c[2], c[3], c[4], c[5] ] )
             obj[TAG].extend( [ COLOR, c[3], c[4], c[5] ] )
             obj[TAG].extend( [ SPHERE, x2, y2, z2, radius ] )
          else:
             a1 = array( [ x1, y1, z1 ] )
             a2 = array( [ x2, y2, z2 ] )

             d = a2 - a1
             mag = sqrt( d[0]*d[0] + d[1]*d[1] + d[2]*d[2] )

             u = unit_vector( d )

             obj[TAG].extend( [ COLOR, c[3], c[4], c[5] ] )
             for i in range(int(mag*10)):
                 obj[TAG].extend( [ SPHERE, a1[0]+0.1*i*u[0],  a1[1]+0.1*i*u[1],  a1[2]+0.1*i*u[2], 0.03 ] )

    print(len(obj.keys()))
    for TAG in obj.keys():
        print(TAG)
        cmd.delete(TAG)
        cmd.load_cgo(obj[TAG],TAG)

cmd.extend( "makepi", run );
