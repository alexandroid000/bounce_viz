import numpy as np
from itertools import combinations

def IsThreePointsOnLine(p1, p2, p3):
    ''' Check whether the input three points are colinear
    '''
    degeneracy_error = 0.00001
    cross_prod = (p1[1] - p2[1]) * (p3[0] - p2[0]) - (p1[0] - p2[0]) * (p3[1] - p2[1])
    return abs(cross_prod)<degeneracy_error

def IsThreePointsOnLineSeg(p1, p2, p3):
    ''' Check whether a point is on the line segment defined two other points

    Parameters
    ----------
    p1:np.array(2)
        line segment point one
    p2:np.array(2)
        line segment point two
    p3:np.array(2)
        the point we want to know whether it is on the line segment

    Returns
    -------
    bool
        true if p3 is in line with and in between p1 and p2
    '''
    if IsThreePointsOnLine(p1, p2, p3):
        return np.dot(p1-p3, p2-p3)<0
    return False

def SegsInGeneralPos(p1, p2, q1, q2):
    ''' return true if the four given points are not in general position
    '''
    tests = (IsThreePointsOnLine(p1, p2, q1) or
           IsThreePointsOnLine(p1, p2, q2) or
           IsThreePointsOnLine(q1, q2, p1) or
           IsThreePointsOnLine(q1, q2, p2))
    return (not tests)

def PolyInGeneralPos(poly):
    ''' Check whether the input polygon is in general position
    '''
    for (a,b,c) in combinations(poly.complete_vertex_list,3):
        if IsThreePointsOnLine(a,b,c):
            return False
    return True
