from helper.geometry_helper import *
from helper.general_position_helper import *
from settings import *

def VertexExists(v, vertex_info):
    ''' Check whether a given vertex exists in the polygon
    '''

    vs = [v for i,v in vertex_info]
    for pt in vs:
        if la.norm(v-pt) < EPSILON:
            return True
    return False

