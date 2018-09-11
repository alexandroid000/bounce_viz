from geom_utils import *
from helper.general_position_helper import *
from settings import *
from helper.point_helper import *

def VertexExists(v, poly):
    for pt in poly:
        if la.norm(v-pt) < EPSILON:
            return True
    return False

def FindReflexVerts(poly):
    ''' return indices of all reflex vertices in poly
    '''
    psize = poly.shape[0]
    reflex_verts = []
    for j in range(psize):
        v1, v2, v3 = poly[(j-1) % psize], poly[j], poly[(j+1) % psize]
        if IsRightTurn(v1,v2,v3) and not IsThreePointsOnLine(v1,v2,v3):
            reflex_verts.append(j)

    return reflex_verts