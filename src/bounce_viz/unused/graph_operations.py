''' Applications of the bounce visibility graph to robotic tasks
'''
from copy import copy

from geom_utils import *
from helper.bounce_graph_helper import *
from maps import *
import numpy.linalg as la

# a strategy should be an automata where inputs = sensor obs, and outputs = bounce angles

def interp(p1, p2, s):
    (x1,y1) = p1
    (x2,y2) = p2
    x_new = (1-s)*x1 + s*x2
    y_new = (1-s)*y1 + s*y2
    return (x_new, y_new)

def polyLens(poly):
    p = list(copy(poly))
    p.append(poly[0])
    edge_lens = [la.norm(p1-p2) for (p1, p2) in zip(p, p[1:])]
    perim_len = sum(edge_lens)
    return edge_lens, perim_len

# [0,1) -> dP
def s(s_param, poly):
    psize = len(poly)
    edge_lens, perim_len = polyLens(poly)
    s_len = s_param*perim_len

    for i in range(psize):
        if s_len - edge_lens[i] < 0.0:
            s_edge = s_len/edge_lens[i]
            if DEBUG:
                print('param', s_edge, 'on edge', poly[i], poly[(i+1)%psize])
            return interp(poly[i],poly[(i+1)%psize],s_edge)
        else:
            s_len = s_len - edge_lens[i]

# dp -> [0,1)
# assumes point is on boundary
def s_inv(pt, poly):
    psize = len(poly)
    edge_lens, perim_len = polyLens(poly)
    edge = 0
    s_edge = 0.0
    for i in range(psize):
        p1,p2 = poly[i], poly[(i+1)%psize]
        if IsThreePointsOnLineSeg(p1,p2,pt) or (p1==pt).all():
            edge = i
            s = sum(edge_lens[:i])/perim_len
            s += la.norm(p1-pt)/perim_len
            return s

def PropagatePath(poly, path, S):
    psize = len(poly)    
    (s1, s2) = S
    p1, p2 = s(s1, poly), s(s2, poly)
    ints = [(p1,p2)]
    print(path)
    for (i,j, ang_range) in path:
        theta_i = FixAngle( np.arctan2(poly[(i+1)%psize][1]-poly[i][1],\
                                 poly[(i+1)%psize][0]-poly[i][0]))
        print('orientation of edge',i,'is', theta_i)
        P1_FOUND = False
        P2_FOUND = False

        if abs(ang_range[0] - pi) < EPSILON:
            next_p2 = poly[i]
            P2_FOUND = True
        if abs(ang_range[1]) < EPSILON:
            next_p1 = poly[(i+1)%psize]
            P1_FOUND = True


        theta_l = FixAngle(theta_i + ang_range[0])
        print('theta_l', theta_l)
        theta_r = FixAngle(theta_i + ang_range[1])
        print('theta_r', theta_r)
        # heuristic to generate vector, hacky
        d = 0.1*la.norm(poly[i]-poly[j])
        (p1,p2) = ints[-1]
        p1theta = (p1[0]+d*cos(theta_l), p1[1]+d*sin(theta_l))
        p2theta = (p2[0]+d*cos(theta_r), p2[1]+d*sin(theta_r))
        if not P1_FOUND:
            print('shoot ray from',p2,'to',p2theta)
            next_p1, k = ClosestPtAlongRay(p2,p2theta,poly,last_bounce_edge=i)
        if not P2_FOUND:
            print('shoot ray from',p1,'to',p1theta)
            next_p2, k = ClosestPtAlongRay(p1,p1theta,poly,last_bounce_edge=i)
        ints.append((next_p1, next_p2))
    return ints

if __name__ == '__main__':
    navigate(square, (0.1,0.2), (0.8,0.9))
