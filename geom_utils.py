#!/usr/bin/env python

# geom_utils.py
# util functions for geometry in simple polygons

import random, math
from math import sqrt,cos,sin,atan2,pi
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from itertools import combinations
import visilibity as vis
from copy import copy, deepcopy

DEBUG = False
EPSILON = 0.00000001

def PolyToWindowScale(poly, ydim):
    newpoly = []
    for p in poly:
        newpoly.append((p[0]/2.0, ydim - p[1]/2.0))
    return newpoly

def PointToWindowScale(p, ydim):
     return (p[0]/2.0, ydim - p[1]/2.0)

# Do we have theta1 <= theta <= theta2 (mod 2pi) ?
# It is just a particular cyclic permutation
def AngleBetween(theta,theta1,theta2):
    return (((theta1 <= theta) and (theta <= theta2)) or \
            ((theta2 <= theta1) and (theta1 <= theta)) or \
            ((theta <= theta2) and (theta2 <= theta1)))

def AngleDistance(theta1,theta2):
    d = abs(theta1-theta2)
    return min(d,2.0*pi-d)

def AngleDifference(theta1,theta2):
    d = theta2 - theta1
    if (d < pi):
        d += 2.0*pi;
    if (d > pi):
        d -= 2.0*pi;
    return d

def FixAngle(theta):
    return theta % (2.0*pi)

def PointDistance(p,q):
    return sqrt((p[0]-q[0])*(p[0]-q[0])+(p[1]-q[1])*(p[1]-q[1]))

def Points2Vect(p1,p2):
    return (p2[0]-p1[0], p2[1]-p1[1])

def GetVectorLen(v):
    return math.sqrt(sum((a*b) for a, b in zip(v, v)))

# angle between two vectors in range [0,pi]
def GetVector2Angle(v1, v2):
    dot_prod = v1[0]*v2[0]+v1[1]*v2[1]
    ratio = max(min(dot_prod/(GetVectorLen(v1)*GetVectorLen(v2)), 1), -1)
    return np.arccos(ratio)

def GetAngleFromThreePoint(p1, p2, origin):
    v1 = (p1[0]-origin[0], p1[1]-origin[1])
    v2 = (p2[0]-origin[0], p2[1]-origin[1])
    return GetVector2Angle(v1, v2)

# r is left of the vector formed by p->q
def IsLeftTurn(p,q,r):
    return q[0]*r[1] + p[0]*q[1] + r[0]*p[1] - \
           (q[0]*p[1] + r[0]*q[1] + p[0]*r[1]) > 0

# r is right of the vector formed by p->q
def IsRightTurn(p,q,r):
    return q[0]*r[1] + p[0]*q[1] + r[0]*p[1] - \
           (q[0]*p[1] + r[0]*q[1] + p[0]*r[1]) < 0

def IsIntersectSegments(p1, p2, q1, q2):
    return (IsLeftTurn(p1, p2, q1) != IsLeftTurn(p1, p2, q2)) and (IsLeftTurn(q1, q2, p1) != IsLeftTurn(q1, q2, p2))

# checks if vector sp->bp points within edge p1p2
# start point, boundary point, edge p1, edge p2
def BouncePointInEdge(sp,bp,ep1,ep2):
    return ((IsLeftTurn(sp,bp,ep2) and
             IsRightTurn(sp,bp,ep1)) or
            (IsRightTurn(sp,bp,ep2) and
             IsLeftTurn(sp,bp,ep1)))

def IsThreePointsOnLine(p1, p2, p3):
    degeneracy_error = 0.00001
    cross_prod = (p1[1] - p2[1]) * (p3[0] - p2[0]) - (p1[0] - p2[0]) * (p3[1] - p2[1])
    return abs(cross_prod)<degeneracy_error

# true if p3 is in line with and in between p1 and p2
def IsThreePointsOnLineSeg(p1, p2, p3):
    if IsThreePointsOnLine(p1, p2, p3):
        v1 = (p1[0]-p3[0], p1[1]-p3[1])
        v2 = (p2[0]-p3[0], p2[1]-p3[1])
        return (v1[0]*v2[0]+v1[1]*v2[1])<0
    return False

# return true if the four given points are not in general position
def SegsInGeneralPos(p1, p2, q1, q2):
    tests = (IsThreePointsOnLine(p1, p2, q1) or
           IsThreePointsOnLine(p1, p2, q2) or
           IsThreePointsOnLine(q1, q2, p1) or
           IsThreePointsOnLine(q1, q2, p2))
    return (not tests)

def PolyInGeneralPos(poly):
    for (a,b,c) in combinations(poly,3):
        if IsThreePointsOnLine(a,b,c):
            return False
    return True

# find line intersection parameter of edge (v1,v2)
# state :: (x,y,theta) initial point and angle of ray
def ShootRay(state, v1, v2):
    x1 = state[0]
    y1 = state[1]
    theta = state[2]
    x2 = v1[0]
    y2 = v1[1]
    x3 = v2[0]
    y3 = v2[1]
    a = y3 - y2
    b = x2 - x3
    c = y2*(x3-x2) - x2*(y3-y2)
    den = (a*cos(theta)+b*sin(theta))
    # case when lines are parallel
    if abs(den) < 0.000001:
        if DEBUG:
            print("divide by zero in ray shoot")
            print("shot from ",state,"to",v1,v2)
        raise ValueError
    else:
        t = (-c - b*y1 - a*x1)/den
        pint = (x1 + cos(theta)*t, y1 + sin(theta)*t)
    return t, pint

# test if point p is in poly using crossing number
def IsInPoly(p, poly):
    intersects = 0
    theta = random.random()*2*pi
    state=(p[0],p[1],theta)
    psize = len(poly)
    for j in range(psize):
        v1, v2 = poly[j], poly[(j+1) % psize]
        try:
            t, pt = ShootRay(state, v1, v2)
            if t>0 and BouncePointInEdge(p, pt, v1, v2):
                intersects += 1
        except:
            pass
    return not (intersects%2 == 0)

def VertexExists(v, poly):
    for pt in poly:
        if PointDistance(v, pt) < EPSILON:
            return True
    return False


# shoot a ray starting at p1, along vector p1->p2
# find intersection with edge v1,v2
# will return first intersection *after* p2
def ShootRayFromVect(p1, p2, v1, v2):
    (x1,y1), (x2,y2) = p1, p2
    (x,y) = (x2-x1, y2-y1) # recenter points so p1 is at origin
    theta = FixAngle(math.atan2(y,x))
    state = (p2[0], p2[1], theta)
    return ShootRay(state, v1, v2)

# shoot ray from p1 toward p2 in poly, return closest intersect point
# will not return point on "last_bounce_edge"
def ClosestPtAlongRay(p1,p2,poly,last_bounce_edge=-1):
    psize = len(poly)
    closest_bounce = 100000000000
    bounce_point = (0.0, 0.0)
    found_coll = False
    bounce_edge = last_bounce_edge

    # check each edge for collision
    for j in range(psize):
        if (j != last_bounce_edge):
            v1, v2 = poly[j], poly[(j+1) % psize]
            x1, y1 = p2[0], p2[1]
            # The line parameter t; needs divide by zero check!
            try:
                t,pt = ShootRayFromVect(p1, p2, v1, v2)
                
                # Find closest bounce for which t > 0
                pdist = PointDistance(pt,p2)
                if ((t > 0) and (pdist < closest_bounce) and
                    BouncePointInEdge((x1,y1),pt,v1,v2)):
                    found_coll = True
                    bounce_point = pt
                    closest_bounce = pdist
                    bounce_edge = j
            # bounce was parallel to edge j
            # check if it's parallel and overlapping for degenerate polys
            except:
                pass
    if found_coll:
        return bounce_point, bounce_edge
    else:
        return False

# return indices of all reflex vertices in poly
def FindReflexVerts(poly):
    psize = len(poly)
    reflex_verts = []

    for j in range(psize):
        v1, v2, v3 = poly[(j-1) % psize], poly[j], poly[(j+1) % psize]
        if IsRightTurn(v1,v2,v3) and not IsThreePointsOnLine(v1,v2,v3):
            reflex_verts.append(j)

    return reflex_verts

# shoot two rays from each reflex vertex, along incident edges
# group by edge so we can insert them in correct order
def ShootRaysFromReflex(poly, j):
    psize = len(poly)
    p2 = poly[j]
    p1_ccw = poly[(j+1) % psize]
    p1_cw = poly[(j-1) % psize]

    int_pts = []

    # ClosestPtAlongRay returns False if we shoot ray into existing vertex
    if ClosestPtAlongRay(p1_ccw,p2,poly,j):
        pt, k = ClosestPtAlongRay(p1_ccw,p2,poly,j)
        if not VertexExists(pt, poly):
            int_pts.append((pt,k))

    if ClosestPtAlongRay(p1_cw,p2,poly,((j-1) % psize)):
        pt, k = ClosestPtAlongRay(p1_cw,p2,poly,((j-1) % psize))
        if not VertexExists(pt, poly):
            int_pts.append((pt,k))

    return int_pts

# shoot ray from visible vertices through reflex verts
# Poly -> Int -> [(Point, Int)]
def ShootRaysToReflexFromVerts(poly, j):
    psize = len(poly)
    r_v = poly[j]
    pts = []
    visible_verts = GetVisibleVertices(poly,j)

    # only ray shoot from non-adjacent vertices
    # previous and next neighbors always visible
    for v in visible_verts:
        if (v != (j-1)%psize) and (v != (j+1)%psize):
            #print("shooting ray from",v,"to",j)
            res = ClosestPtAlongRay(poly[v], r_v, poly)
            if res:
                pt, k = res
                if (IsInPoly(((pt[0]+r_v[0])/2, (pt[1]+r_v[1])/2), poly) and
                    not VertexExists(pt, poly)):
                    #print("successful insert")
                    pts.append((pt,k))
    return pts


# for polygons not in general position: returns all visible vertices along a ray
def GetVisibleVertices(poly, j):
    #print("At vertex",j)
    psize = len(poly)

    # Outer boundary polygon must be COUNTER-CLOCK-WISE(ccw)
    # Create the outer boundary polygon
    # Define an epsilon value (should be != 0.0)
    vpoly = [vis.Point(*pt) for pt in poly]
    walls = vis.Polygon(vpoly)
    walls.enforce_standard_form()

    # point from which to calculate visibility
    p1 = poly[j]
    vp1 = vis.Point(*p1)

    # Create environment, wall will be the outer boundary because
    # is the first polygon in the list. The other polygons will be holes
    env = vis.Environment([walls])
    # Necesary to generate the visibility polygon
    vp1.snap_to_boundary_of(env, EPSILON)
    vp1.snap_to_vertices_of(env, EPSILON)
    isovist = vis.Visibility_Polygon(vp1, env, EPSILON)
    vvs = [(isovist[i].x(), isovist[i].y()) for i in range(isovist.n())]
    #print(vvs)
    visibleVertexSet = []
    for i in range(psize):
        p = vis.Point(*poly[i])
        vp = copy(p).projection_onto_boundary_of(isovist)
        if (vis.distance(vp, p) < EPSILON) and i != j:
            visibleVertexSet.append(i)
    if DEBUG:
        print("vertices",visibleVertexSet,"visible from",j)
    return visibleVertexSet

# sort by distance and remove duplicates
def SortByDistance(p1, unsorted_vs):
    unsorted_vs.append(p1)
    verts = list(set(unsorted_vs))
    return sorted(verts, key = lambda v: PointDistance(v,p1))

# find all induced transition points, sort and insert into polygon
# probably the most naive way to do this
# return a list of edge indicators for each vertex
def InsertAllTransitionPts(poly):
    t_pts_grouped = {i:[] for i in range(len(poly))}
    rvs = FindReflexVerts(poly)
    # find all transition points, group by edge
    for p in rvs:
        t_pts = []
        t_pts.extend(ShootRaysFromReflex(poly, p))
        t_pts.extend(ShootRaysToReflexFromVerts(poly,p))
        for (pt,edge_i) in t_pts:
            t_pts_grouped[edge_i].append(pt)

    # sort transition points along edge
    new_poly_grouped = {}
    for i in range(len(poly)):
        new_poly_grouped[i] = SortByDistance(poly[i], t_pts_grouped[i])
        if DEBUG:
            print("Inserted verts", new_poly_grouped[i], "on edge",i)

    # insert into polygon
    new_poly = []
    for i in range(len(poly)):
        new_poly.extend(new_poly_grouped[i])
    return new_poly

# make all sets of vertices visible from everywhere along edge
def mkVizSets(poly):
    psize = len(poly)
    all_viz_vxs = [GetVisibleVertices(poly, i) for i in range(psize)]
    if DEBUG:
        print("All visible verts:")
        print(all_viz_vxs)
        print("\n")

    # each element in the map is indexed by its clockwise vertex (smaller index)
    vizSets = {i:[] for i in range(psize)}

    # get vertices that are visible to the current vertex and the next vertex
    for i in range(psize):
        viz_vxs = list(set(all_viz_vxs[i]) & set(all_viz_vxs[(i+1)%psize]))
        # if the following line included, allows transition to next edge by wall
        # following, even if reflex angle
        # TODO: figure out if we want to allow this behavior
        viz_vxs.append((i+1)%psize) 
        vizSets[i] = viz_vxs

    if DEBUG:
        print("viz sets")
        print("--------")
        print(vizSets)
        print("")
    return vizSets 

# Resolution is the number of sample points on each edge
def GetLinkDiagram(poly, resolution = 15):
    psize = len(poly)
    link_diagram = np.nan*np.ones((psize, resolution*psize))
    vizSets = mkVizSets(poly)
    for i in range(psize):
        # for each visible vertex, we need to calculate:
        #  (1) the view angle at the current vertex w.r.t the current edge and
        #  (2) the view angle at the sample points on the edge and the next vertex w.r.t the current edge
        #  (3) insert a np.nan for discontinuity otherwise matplotlib will try to connect them together
        for vx in vizSets[i]:
            curr_p = poly[i]
            next_p = poly[(i+1)%psize]
            interest_p = poly[vx]

            link_diagram[vx][resolution*i] = GetAngleFromThreePoint(next_p, interest_p, curr_p)
            # the if case is for when we are considering the "next vertex"
            if interest_p == next_p:
                for stride in range(resolution)[1:-1]:
                    link_diagram[vx][resolution*i+stride] = link_diagram[vx][resolution*i]
            else:
                for stride in range(resolution)[1:-1]:
                    # this is the point on the same line of the current edge but in front of the next vertex.
                    # Use this so that we don't get zero length vector from the next vertex perspective
                    extended_point = (curr_p[0]+2*(next_p[0]-curr_p[0]), curr_p[1]+2*(next_p[1]-curr_p[1]))
                    ratio = 1./(resolution-2)*stride
                    origin = (ratio*next_p[0]+(1-ratio)*curr_p[0], ratio*next_p[1]+(1-ratio)*curr_p[1])
                    link_diagram[vx][resolution*i+stride] = GetAngleFromThreePoint(interest_p, extended_point, origin)
            link_diagram[vx][resolution*i+resolution-1] = np.nan
    return link_diagram

# return minimum and maximum angles that allow transition from e1 to e2
# from *somewhere* on e1
# only need to check endpoints
def AnglesBetweenSegs(e1, e2):
    (p1,p2) = e1
    (p3,p4) = e2
    if DEBUG:
        print("finding angle between",e1,"and",e2)

    min_ang = GetVector2Angle(Points2Vect(*e1),Points2Vect(p1,p3))
    max_ang = GetVector2Angle(Points2Vect(*e1),Points2Vect(p2,p4))

    return min_ang, max_ang

# return angle range that allows transition from e1 to e2
# from *anywhere* on e1
# or return None if there is no such range
def SafeAngles(e1, e2):
    (p1,p2) = e1
    (p3,p4) = e2
    if p1 == p4:
        theta_l = pi
        theta_r = GetVector2Angle(Points2Vect(*e1),Points2Vect(p2,p3))
    elif p2 == p3:
        theta_r = 0.0
        theta_l = GetVector2Angle(Points2Vect(*e1),Points2Vect(p1,p4))
    else:
        theta_l = GetVector2Angle(Points2Vect(*e1),Points2Vect(p1,p4))
        theta_r = GetVector2Angle(Points2Vect(*e1),Points2Vect(p2,p3))
    not_parallel = abs(theta_l - theta_r) > 0.01
    if theta_l > theta_r and not_parallel: # declare lines parallel if within 1 deg
        return (theta_l, theta_r)
    else:
        return None

#def MaxPreimage(e1, e2):


