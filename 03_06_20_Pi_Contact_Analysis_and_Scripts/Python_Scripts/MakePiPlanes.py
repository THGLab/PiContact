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

def unit_vector(vector):
    #Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    #Returns the angle in radians between vectors 'v1' and 'v2'::
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    angle = np.arccos(np.dot(v1_u, v2_u))
    if np.isnan(angle):
        if (v1_u == v2_u).all():
            return 0.0
        else:
            return np.pi
    return angle
   


def read_pdb( file ): # Fills dictionary with D{chain+"."res number} = [ aa, {atoms} = [x,y,z] ]
    pdb = open( file ).readlines()
    pdbdict = {}
    hasbase = False
    base = {"U":0,"C":0,"A":0,"G":0}
    for i in pdb:
        line = i.split()

        if len(line) > 5:
            if i[:4] == "ATOM":
                atom = i[12:16].split()[0]
                res = i[17:20].split()[0]
                chain = i[21:22]
                x = float(i[30:38])
                y = float(i[38:46])
                z = float(i[46:54])
                resN = int(i[22:26])

                key = chain+"."+str(resN)

                if (res in rgroups ):

                    if (res in base ):
                        hasbase = True

                    if (key in pdbdict ) == False:
                        pdbdict[key] = [res, {}]

                    if (atom in pdbdict[key][1] ) == False:
                        pdbdict[key][1][atom] = array([x,y,z])
    return pdbdict, hasbase

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

def read_hetatm( file ):
    pdb = open( file ).readlines()
    pdbdict = {}
    for i in pdb:
        line = i.split()

        if len(line) > 5:
            if i[:6] == "HETATM":
                atom = i[12:16].split()[0]
                res = i[17:20].split()[0]
                chain = i[21:22]
                x = float(i[30:38])
                y = float(i[38:46])
                z = float(i[46:54])
                resN = int(i[22:26])

                key = chain+"."+str(resN)

                if (key in pdbdict) == False:
                    pdbdict[key] = [res, {}]

                if (atom in pdbdict[key][1]) == False:
                    pdbdict[key][1][atom] = array([x,y,z])

    return pdbdict



def is_complete( pdbdict ):
    for key in pdbdict.keys():
        complete = True
        aa = pdbdict[key][0]
        for a in rgroups[aa]:
            if (a in pdbdict[key][1]) == False:
                complete = False

        if aa != "A" and aa != "G" and aa != "C" and aa != "U":
            for a in rgroups["BB"]:
                if (a in pdbdict[key][1]) == False:
                    complete = False
            if ('N' in pdbdict[key][1]) == False:
                complete = False
        
        pdbdict[key].append(complete)

  
def peptide_details(res1,res2):#N, CA, C, O, Np, CAp, Cp):
    
    N = res1["N"]
    CA = res1["CA"]
    C = res1["C"]
    O = res1["O"]
    Np = res2["N"]
    CAp = res2["CA"]
    Cp = res2["C"]

    psi = torsion_angle(N,CA,C,Np)
    #while psi < 0:
    #    psi += pi

    phi_p1 = torsion_angle(C,Np,CAp,Cp)
    #while phi_p1 < 0:
    #    phi_p1 += pi

    omega = torsion_angle(CA,C,Np,CAp)
    #while omega < 0:
    #    omega += pi

    vecs = make_axis_vectors(O, C, CA)

    normvec = vecs[1]

    center = (CA+C+O+Np+CAp)/5.0

    alist = [C,O,Np]

    return [np.degrees(psi), np.degrees(phi_p1), np.degrees(omega), normvec, center, alist]

def compare_res_to_BB( bb, res ):
    return compare_res_to_BB_type( bb, res, "default")

def NEWcompare_res_to_res( res1, res2 ):

    vecs1 = make_axis_vectors( res1[1][rgroups[res1[0]][0]],  res1[1][rgroups[res1[0]][1]],  res1[1][rgroups[res1[0]][2]] )
    vecs2 = make_axis_vectors( res2[1][rgroups[res2[0]][0]],  res2[1][rgroups[res2[0]][1]],  res2[1][rgroups[res2[0]][2]] )

    nvec1 = vecs1[1]
    nvec2 = vecs2[1]

    resdict1 = {}
    for i in agroups[res1[0]]:
        resdict1[i] = res1[1][i]

    resdict2 = {}
    for i in agroups[res2[0]]:
        resdict2[i] = res2[1][i]

    return compare_dict_to_dict( resdict1, nvec1, resdict2, nvec2 )


def NEWcompare_res_to_BB( bb, bbn, res ):

    vecs1 = make_axis_vectors( res[1][rgroups[res[0]][0]],  res[1][rgroups[res[0]][1]],  res[1][rgroups[res[0]][2]] )

    vecs2 = make_axis_vectors( bb[1]["CA"],  bb[1]["C"],  bb[1]["O"] )
    nvec1 = vecs1[1]
    nvec2 = vecs2[1]

    resdict = {}
    for i in agroups[res[0]]:
        resdict[i] = res[1][i]


    bbdict = {}
    bbdict["C"] = bb[1]["C"]
    bbdict["O"] = bb[1]["O"]
    bbdict["N"] = bbn[1]["N"]

    return compare_dict_to_dict( resdict, nvec1, bbdict, nvec2 )


def NEWcompare_BB_to_BB( bb1, bbn1, bb2, bbn2 ):

    vecs1 = make_axis_vectors( bb1[1]["CA"],  bb1[1]["C"],  bb1[1]["O"] )
    vecs2 = make_axis_vectors( bb2[1]["CA"],  bb2[1]["C"],  bb2[1]["O"] )
    nvec1 = vecs1[1]
    nvec2 = vecs2[1]

    bbdict1 = {}
    bbdict1["C"] = bb1[1]["C"]
    bbdict1["O"] = bb1[1]["O"]
    bbdict1["N"] = bbn1[1]["N"]

    bbdict2 = {}
    bbdict2["C"] = bb2[1]["C"]
    bbdict2["O"] = bb2[1]["O"]
    bbdict2["N"] = bbn2[1]["N"]

    return compare_dict_to_dict( bbdict1, nvec1, bbdict2, nvec2 )



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
    

def compare_res_to_BB_type( bb, res, rtype ):

    rg = rgroups[res[0]]#rgroups[res2[0]]
    ag = agroups[res[0]]

    A = res[1][rg[0]]
    B = res[1][rg[1]]
    C = res[1][rg[2]]

    rvecs = make_axis_vectors(A, B, C) 

    #dot_toBB = vdot(r1vecs[0], r2vecs[0]) # A-B bond vectors
    dot_orthog = (vdot(bb[3], rvecs[1])) # vectors perpendicular to plane
    #NOTE: Absolute value because side chain plane vector is based on degenerate atoms

    #dot_90deg = vdot(r1vecs[2], r2vecs[2]) # vectors 90 degrees to A-B bond in plane

    center = array([0.0,0.0,0.0])
    for xi in rgroups[res[0]]:
        center += res[1][xi]
    center /= float(len(rgroups[res[0]]))

    vd = center - bb[4]
    vdlen = sqrt(vd[0]**2 + vd[1]**2 + vd[2]**2)
    vd[0] = vd[0] / vdlen
    vd[1] = vd[1] / vdlen
    vd[2] = vd[2] / vdlen

    dot_rcomp1 = vdot(vd, bb[3])# # res-bb vector orientation to bb_orthog
    #dot_rcomp2 = vdot(vd, r2vecs[1]) # res1B-res2B vector orientation to res2_orthog 

    min_distance = -1.0
    for xi in rg:
        #for yi in rg2:#res2[1].keys():
        dp = res[1][xi] - bb[4]#res2[1][yi]
        distance = sqrt(dp[0]**2 + dp[1]**2 + dp[2]**2)

        if distance < min_distance or min_distance == -1.0:
            min_distance = distance

    abmin_pi_distance = -1.0
    pairmatch = [0,0]
    rnvec = rvecs[1]
    bnvec = bb[3]

    udlist1 = []
    udlist2 = []
    plen = 1.7
    for y in range(len(bb[5])):#yi in bb[5]:
        yi = bb[5][y]
        bb_up = yi + (bnvec * plen)
        bb_down = yi - (bnvec * plen)
        udlist1.append([bb_up,bb_down])

    for x in range(len(ag)):#xi in rg:
        xi = ag[x]
        res_up = res[1][xi] + (rnvec * plen)
        res_down = res[1][xi] - (rnvec * plen)
        udlist2.append([res_up,res_down])

    pairmatch[0] = calc_match(udlist1,udlist2)

    if rtype == "points":
        udlistS = {}
        for xi in agroups[res[0]]:
            x = res[1][xi]
            plen = 1.7
            u = x + (rnvec*plen)
            d = x - (rnvec*plen)
            udlistS[xi] = [u,d]
             
        centerB = array([0.0,0.0,0.0])

        bnD = ["C","O","N"]
        
        udlistB = {}
        for yi in range(len(bb[5])):#yi in bb[5]:
            y = bb[5][yi]
            centerB += y
            plen = 1.7
            u = y + (bnvec*plen)
            d = y - (bnvec*plen)
            udlistB[bnD[yi]] = [u,d]
        centerB /= float(len(bb[5]))

        return get_closepoints(udlistS, udlistB, center, centerB)

    return [dot_orthog, dot_rcomp1, min_distance, abmin_pi_distance, pairmatch]


def torsion_angle(A, B, C, D):
        
    b1 = C - D
    b2 = B - C
    b3 = A - B

    if b2[0] != 0 and b2[1] != 0 and b2[2] != 0:
    #if abs(b2[0]) > 0 or abs(b2[1]) > 0 or abs(b2[2]) > 0:
        left = dot(cross( cross(b1,b2), cross(b2,b3) ), b2/abs(b2))
        right = dot(cross(b1,b2),cross(b2,b3))
        torangle = np.arctan2(left,right) 
        return torangle
    else:
        return 0.0

def make_axis_vectors(A, B, C): # For np.array([x,y,z]) of connected atoms A->B->C
    
    Avec = A-B # vector along A->B bond
    
    r1cen = A - B#pdbdict[p1][2][rg1[0][0]] - pdbdict[p1][2][rg1[0][1]] 
    r2cen = C - B#pdbdict[p1][2][rg1[1][1]] - pdbdict[p1][2][rg1[0][1]]

    Nvec = cross(r1cen, r2cen) # vector parralel to plane

    r1cen = A - B
    r2cen = Nvec

    Cvec = cross(r1cen, r2cen) # 90 degree angle to Avec, in plane with A,B,C

    Alen = sqrt(Avec[0]**2+Avec[1]**2+Avec[2]**2)
    Avec /= Alen

    Nlen = sqrt(Nvec[0]**2+Nvec[1]**2+Nvec[2]**2)
    Nvec /= Nlen

    Clen = sqrt(Cvec[0]**2+Cvec[1]**2+Cvec[2]**2)
    Cvec /= Clen

    return [Avec, Nvec, Cvec]

#   Tests whether or not a residue has all the atoms it needs... function is broken and irrelevant though
def all_compare( res, rgtype):
    rg = rgroups["BB"]

    rg.append("N")
    
    for i in rg:
        if (i in res[1]) == False:
            return False

    if can_compare( res, rgtype ) == False:
        return False
            
    return True

#   Tests the residue data to make sure there are coordinates for each of the
#   atom types in the rgtype list.
#
#   For "BB" (backbone) type this also checks if it has a nitrogen
#
#   This was, of course, meant to be the next nitrogen in the chain, but other completeness issues
#   caused this function to become irrelevant, so it was never fixed. Instead we complete the PDBs before
#   using them.
def can_compare( res, rgtype):
    rg = rgroups[rgtype]
    
    if rgtype == "BB":
        rg.append("N")
    
    for i in rg:
        if (i in res[1]) == False:
            return False
            
    return True
  
#   The distance between the second atoms in the rgroups lists
#   which for this function we can just assume is a random atom
#   specific to the residue type in question
def rough_distance( res1, res2 ):
    if (rgroups[res1[0]][1] in res1[1]) and (rgroups[res2[0]][1] in res2[1]):
        dp = res1[1][rgroups[res1[0]][1]] - res2[1][rgroups[res2[0]][1]]
        distance = sqrt(dp[0]**2 + dp[1]**2 + dp[2]**2)

        return distance
    else:
        return 1000000.0

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

def new_get_closepoints( udlist1, udlist2, center1, center2):

    l1 = []
    for i in udlist1:
        dp = center2 - i[0]
        distanceU = sqrt(dp[0]**2 + dp[1]**2 + dp[2]**2)
        dp = center2 - i[1]
        distanceD = sqrt(dp[0]**2 + dp[1]**2 + dp[2]**2)

        if distanceU < distanceD:
            l1.append( i[0] )
        else:
            l1.append( i[1] )

    l2 = []
    for i in udlist2:
        dp = center1 - i[0]
        distanceU = sqrt(dp[0]**2 + dp[1]**2 + dp[2]**2)
        dp = center1 - i[1]
        distanceD = sqrt(dp[0]**2 + dp[1]**2 + dp[2]**2)

        if distanceU < distanceD:
            l2.append( i[0] )
        else:
            l2.append( i[1] )

    return [l1,l2]



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
        return pointlist#ha ha it's pointlist

def compare_residues( res1, res2 ):
    return calc_residue_match(res1,res2,"")

def get_above( resC, resU ):
    rgC = rgroups[resC[0]]
    A = resC[1][rgC[0]]
    B = resC[1][rgC[1]]
    C = resC[1][rgC[2]]
    rCvecs = make_axis_vectors(A, B, C) 

    centerC = array([0.0,0.0,0.0])
    for xi in rgroups[resC[0]]:
        centerC += resC[1][xi]
    centerC /= float(len(rgroups[resC[0]]))

    centerU = array([0.0,0.0,0.0])
    for xi in rgroups[resU[0]]:
        centerU += resU[1][xi]
    centerU /= float(len(rgroups[resU[0]]))

    dot_vsorthogL = (vdot(rCvecs[1], unit_vector(centerU-centerC)))

    orthog = rCvecs[1]
    if dot_vsorthogL > 0:
        orthog = orthog*-1


    test_vsorthogL = (vdot(rCvecs[1], unit_vector((centerC+orthog) - centerC)))

    above_plane = centerC+orthog

    return centerC, above_plane

def get_orthog_points( res ):
    above = []
    center = []
    below = []

    acb_points = {}

    rg = rgroups[res[0]]
    A = res[1][rg[0]]
    B = res[1][rg[1]]
    C = res[1][rg[2]]
    rvecs = make_axis_vectors(A, B, C) 
    orthog = 1.7*rvecs[1]

    for ag in agroups[res[0]]:
        acb_points[ag] = []

        acb_points[ag].append(res[1][ag]+orthog)
        acb_points[ag].append(res[1][ag])
        acb_points[ag].append(res[1][ag]-orthog)

    return acb_points


def get_6_above( resC, resU ):
    rgC = rgroups[resC[0]]
    A = resC[1][rgC[0]]
    B = resC[1][rgC[1]]
    C = resC[1][rgC[2]]
    rCvecs = make_axis_vectors(A, B, C) 

    centerC = array([0.0,0.0,0.0])
    for xi in rgroups[resC[0]]:
        centerC += resC[1][xi]
    centerC /= float(len(rgroups[resC[0]]))

    centerU = array([0.0,0.0,0.0])
    for xi in rgroups[resU[0]]:
        centerU += resU[1][xi]
    centerU /= float(len(rgroups[resU[0]]))

    dot_vsorthogL = (vdot(rCvecs[1], unit_vector(centerU-centerC)))

    orthog = 1.7*rCvecs[1]
    if dot_vsorthogL > 0:
        orthog = orthog*-1

    test_vsorthogL = (vdot(rCvecs[1], unit_vector((centerC+orthog) - centerC)))

    above_plane = []#centerC+orthog
    center_match = []

    for xi in rgroups[resC[0]]:
        if xi != "CB":
            above_plane.append( resC[1][xi]+orthog )
            center_match.append( resC[1][xi] )

    return center_match, above_plane


def calc_sandwich( resL, res, resR ):
    rgC = rgroups[res[0]]
    A = res[1][rgC[0]]
    B = res[1][rgC[1]]
    C = res[1][rgC[2]]
    rCvecs = make_axis_vectors(A, B, C) 

    rg1 = rgroups[resL[0]]
    A1 = resL[1][rg1[0]]
    B1 = resL[1][rg1[1]]
    C1 = resL[1][rg1[2]]
    r1vecs = make_axis_vectors(A1, B1, C1) 

    rg2 = rgroups[resR[0]]
    A2 = resR[1][rg2[0]]
    B2 = resR[1][rg2[1]]
    C2 = resR[1][rg2[2]]
    r2vecs = make_axis_vectors(A2, B2, C2) 

    dot_vsorthogL = (vdot(rCvecs[1], unit_vector(B1-B)))
    dot_vsorthogR = (vdot(rCvecs[1], unit_vector(B2-B)))

    return abs(dot_vsorthogL - dot_vsorthogR)


#   Finds the smallest number of atoms in two sidechains that satisfy the match criteria
#   where a match is when points 1.7 angstroms above the plane for each atom (ie: the carbon VDW surface) 
#   are within 1.5 angstroms of each other
#
#   Default return style puts out a list of other useful variables
#   [dot_orthog, dot_rcomp1, dot_rcomp2, abmin_distance, num_match]
#
#       dot_orthog
#           dot product orientation between the two planes (1.0 = identical, 0.0 = 90 degree angle) 
#           (you can get the angle between the planes by putting this into acos(x))
#
#       dot_rcomp1
#           dot product between residue 1's above-the-plane unit vector and the unit vector between the 
#           two sidechain centers of mass. So this describes residue2's position relative to residue1's plane
#           with 1.0 = directly above the plane and 0.0 = alongside the plane
#
#       dot_rcomp2
#           same as above, but using residue 2's above-the-plane unit vector, so it describes residue1's position
#           relative to residue2 this time.
#
#       abmin_distance
#           smallest distance between two atoms
# 
#       num_match
#           the number of atoms satisfying the match criteria
#
#   the "style" option allows you to switch behavior.
#
#   match = just returns num_match
#
#   points = just spitting out a list of points representing the 1.7 angstrom above the plane coordinates closest 
#   to the other system for both sidechains... (no distance cutoff)
def calc_residue_match( res1, res2, style ):
    rg1 = rgroups[res1[0]]
    rg2 = rgroups[res2[0]]

    A1 = res1[1][rg1[0]]
    B1 = res1[1][rg1[1]]
    C1 = res1[1][rg1[2]]

    r1vecs = make_axis_vectors(A1, B1, C1) 

    A2 = res2[1][rg2[0]]
    B2 = res2[1][rg2[1]]
    C2 = res2[1][rg2[2]]

    r2vecs = make_axis_vectors(A2, B2, C2) 

    #dot_toBB = vdot(r1vecs[0], r2vecs[0]) # A-B bond vectors
    dot_orthog = (vdot(r1vecs[1], r2vecs[1])) # vectors perpendicular to plane
    #dot_90deg = vdot(r1vecs[2], r2vecs[2]) # vectors 90 degrees to A-B bond in plane

    center1 = array([0.0,0.0,0.0])
    for xi in rgroups[res1[0]]:
        center1 += res1[1][xi]
    center1 /= float(len(rgroups[res1[0]]))
        
    center2 = array([0.0,0.0,0.0])
    for xi in rgroups[res2[0]]:
        center2 += res2[1][xi]
    center2 /= float(len(rgroups[res2[0]]))

    #vd = res1[1][rg1[1]] - res2[1][rg2[1]]
    vd = center1 - center2

    vdlen = sqrt(vd[0]**2 + vd[1]**2 + vd[2]**2)
    vd[0] = vd[0] / vdlen
    vd[1] = vd[1] / vdlen
    vd[2] = vd[2] / vdlen

    dot_rcomp1 = vdot(vd, r1vecs[1]) # res1B-res2B vector orientation to res1_orthog 
    dot_rcomp2 = vdot(vd, r2vecs[1]) # res1B-res2B vector orientation to res2_orthog 

    abmin_distance = -1.0
    for xi in rg1:#res1[1].keys():
        for yi in rg2:#res2[1].keys():
            dp = res1[1][xi] - res2[1][yi]
            distance = sqrt(dp[0]**2 + dp[1]**2 + dp[2]**2)

            if distance < abmin_distance or abmin_distance == -1.0:
                abmin_distance = distance

    r1nvec = r1vecs[1]
    r2nvec = r2vecs[1]

    if (res1[0] in agroups) and (res2[0] in agroups):
        ag1 = agroups[res1[0]]
        ag2 = agroups[res2[0]]

        rg1_nlist = []
        for xi in ag1:
            x = res1[1][xi]
            plen = 1.7
            u = x + (r1nvec*plen)
            d = x - (r1nvec*plen)

            rg1_nlist.append([u,d])

        rg2_nlist = []
        for xi in ag2:
            x = res2[1][xi]
            plen = 1.7
            u = x + (r2nvec*plen)
            d = x - (r2nvec*plen)

            rg2_nlist.append([u,d])

    

    if style == "points":

        if (res1[0] in agroups) and (res2[0] in agroups):
            ag1 = agroups[res1[0]]
            ag2 = agroups[res2[0]]

            rg1_nlist = {}
            for xi in ag1:
                x = res1[1][xi]
                plen = 1.7
                u = x + (r1nvec*plen)
                d = x - (r1nvec*plen)

                rg1_nlist[xi] = [u,d]

            rg2_nlist = {}
            for xi in ag2:
                x = res2[1][xi]
                plen = 1.7
                u = x + (r2nvec*plen)
                d = x - (r2nvec*plen)

                rg2_nlist[xi] = [u,d]

        closepoints = get_closepoints(rg1_nlist, rg2_nlist, center1, center2)
        return closepoints
        
    num_match = calc_match( rg1_nlist, rg2_nlist )
    if style == "match":
        return num_match

    return [dot_orthog, dot_rcomp1, dot_rcomp2, abmin_distance, num_match]



rgroups = {}

#No SP3s, these atom sets are random planes I was using as a control in the hoary days of olde
#but I ended up using the keys in this dictionary as my standin list of amino acids so here they are
rgroups["ALA"] = ["C", "CB", "CA"]
rgroups["LEU"] = ["CD1", "CG", "CD2", "CB"]
rgroups["VAL"] = ["CG1", "CB", "CG2"]
rgroups["ILE"] = ["CD1", "CG1", "CB", "CG2"]
rgroups["MET"] = ["CE", "SD", "CG", "CB"]
rgroups["CYS"] = ["SG", "CB", "CA"]
rgroups["SER"] = ["OG", "CB", "CA"]
rgroups["THR"] = ["OG1", "CB", "CG2"]
rgroups["LYS"] = ["NZ", "CE", "CD", "CG", "CB"]

#This one exists to prevent index errors when testing against the GLY sidechain, which... uh, I'm not sure that came up...
rgroups["GLY"] = ["O", "C", "CA"]

#Sidechain pi-system planes are computed based on these atom lists
rgroups["ASN"] = ["CB", "CG", "OD1", "ND2"]
rgroups["GLN"] = ["CG", "CD", "OE1", "NE2"]
rgroups["ASP"] = ["CB", "CG", "OD1", "OD2"]
rgroups["GLU"] = ["CG", "CD", "OE1", "OE2"]
rgroups["ARG"] = ["NE", "CZ", "NH1", "NH2"]#NOTE THE MISSING CD
rgroups["PRO"] = ["CA", "N", "CD", "CG", "CB"]
rgroups["HIS"] = ["CB", "CG", "ND1", "CD2", "CE1", "NE2"]
rgroups["PHE"] = ["CB","CG","CD1","CD2","CE1","CE2","CZ"]
rgroups["TYR"] = ["CB","CG","CD1","CD2","CE1","CE2","CZ"]
rgroups["TRP"] = ["CB", "CG", "CD2", "CD1", "NE1", "CE2", "CZ2", "CH2", "CZ3", "CE3"]

#(one of the functions checks if the residue in question has all the atoms it needs before trying to compute
# pi-system planes and stuff, and when that function looks at backbones it uses this for the first residue and
# then a special N check against the second residue... I don't know man, this system kind of grew organically here...)
rgroups["BB"] = ["O", "C", "CA"]

rgroups["U"] = ["C1'", "N1", "C6", "C2", "O2", "N3", "C4", "O4", "C5"]
rgroups["C"] = ["C1'", "N1", "C6", "C2", "O2", "N3", "C4", "N4", "C5"]
rgroups["A"] = ["C1'", "N9", "C8", "C4", "N3", "C2", "N1", "C6", "C5", "N7", "N6"]
rgroups["G"] = ["C1'", "N9", "C8", "C4", "N3", "C2", "N1", "C6", "C5", "N7", "N2", "O6"]

agroups = {}

#SP3 Sidechain Atoms (basically rgroups minus the atom vector that points towards the backbone)
agroups["ASN"] = ["CG", "OD1", "ND2"]
agroups["GLN"] = ["CD", "OE1", "NE2"]
agroups["ASP"] = ["CG", "OD1", "OD2"]
agroups["GLU"] = ["CD", "OE1", "OE2"]
agroups["ARG"] = ["NE", "CZ", "NH1", "NH2"]
agroups["HIS"] = ["CG", "ND1", "CD2", "CE1", "NE2"]
agroups["PHE"] = ["CG","CD1","CD2","CE1","CE2","CZ"]
agroups["TYR"] = ["CG","CD1","CD2","CE1","CE2","CZ"]
agroups["TRP"] = ["CG", "CD2", "CD1", "NE1", "CE2", "CZ2", "CH2", "CZ3", "CE3"]

agroups["U"] = ["N1", "C6", "C2", "O2", "N3", "C4", "O4", "C5"]
agroups["C"] = ["N1", "C6", "C2", "O2", "N3", "C4", "N4", "C5"]
agroups["A"] = ["N9", "C8", "C4", "N3", "C2", "N1", "C6", "C5", "N7", "N6"]
agroups["G"] = ["N9", "C8", "C4", "N3", "C2", "N1", "C6", "C5", "N7", "N2", "O6"]
