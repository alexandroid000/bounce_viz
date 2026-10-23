''' Utility functions for geometry in the plane
'''
import numpy as np
import numpy.linalg as la

# Do we have theta1 <= theta <= theta2 (mod 2pi) ?
# It is just a particular cyclic permutation
def FixAngle(theta):
    return theta % (2.0*np.pi)

# angle between two vectors in range [0,pi]
def AngleBetween(v1, v2):
    ratio = max(min(np.dot(v1, v2)/(la.norm(v1)*la.norm(v2)), 1), -1)
    return np.arccos(ratio)

def Cross2d(p,q):
        return p[0]*q[1] - p[1]*q[0]

# r is left of the vector formed by p->q
def IsLeftTurn(p,q,r):
    return q[0]*r[1] + p[0]*q[1] + r[0]*p[1] - \
           (q[0]*p[1] + r[0]*q[1] + p[0]*r[1]) > 0

# r is right of the vector formed by p->q
def IsRightTurn(p,q,r):
    return q[0]*r[1] + p[0]*q[1] + r[0]*p[1] - \
           (q[0]*p[1] + r[0]*q[1] + p[0]*r[1]) < 0

def sort_by_distance(p1, unsorted_vs):
    ''' sort the input vertices by distance and remove duplicates

    '''
    if unsorted_vs.shape[0] == 0:
        unsorted_vs = p1
    else:
        unsorted_vs = np.vstack((unsorted_vs, p1))
    verts = list(np.unique(np.reshape(unsorted_vs, (-1, 2)), axis = 0))
    return sorted(verts, key = lambda v: la.norm(v-p1))

def get_angle_from_three_pts(p1, p2, origin):
    return AngleBetween(p1-origin, p2-origin)
