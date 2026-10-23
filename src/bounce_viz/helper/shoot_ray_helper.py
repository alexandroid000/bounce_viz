from helper.polygon_helper import *
from helper.visibility_helper import visibleVertices
from settings import *
import numpy as np

def ShootRay(state, v1, v2):
    '''
    # find line intersection parameter of edge (v1,v2)
    # state :: (x,y,theta) initial point and angle of ray
    # (x,y,theta) -> Point -> Point -> Maybe (Double, Point)
    # theta is defined with 0 along the y axis
    # https://stackoverflow.com/questions/14307158/how-do-you-check-for-intersection-between-a-line-segment-and-a-line-ray-emanatin/32146853
    # https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect/565282#565282
    '''
    pt = np.array([state[0], state[1]])
    theta = state[2]

    # points on ray are (x1,y1) + t*r
    r = (np.cos(theta), np.sin(theta))
    # points on segment are (x2,y2) + u*s
    s = v2-v1
    rXs = Cross2d(r,s)

    # if ray and target edge are parallel, will get divide by zero
    if abs(rXs) < EPSILON:
        if DEBUG:
            print('divide by zero in ray shoot')
            print('shot from ',state,'to',v1,v2)
        raise ValueError
    else:
        u = Cross2d(v1-pt, r)/rXs
        t = Cross2d(v1-pt, s)/rXs
        pint = np.array([pt[0] + np.cos(theta)*t, pt[1] + np.sin(theta)*t])
        return t, u, pint

def ShootRayFromVect(p1, p2, v1, v2):
    ''' 
    - shoot a ray starting at p1, along vector p1->p2
    - find intersection with edge v1,v2
    - will return first intersection *after* p2
    '''
    [x,y] = p2-p1 # recenter points so p1 is at origin
    theta = FixAngle(np.arctan2(y,x))
    state = (p2[0], p2[1], theta)
    return ShootRay(state, v1, v2)

def RayInEdge(sp,bp,ep1,ep2):
    '''
    checks if vector sp->bp points within edge p1p2
    Parameters
    ----------
    sp : numpy.array
        start point
    bp : numpy.array
        boundary point
    ep1 : numpy.array
        edge p1
    ep2 : numpy.array
        edge p2
    '''
    return ((IsLeftTurn(sp,bp,ep2) and
             IsRightTurn(sp,bp,ep1)) or
            (IsRightTurn(sp,bp,ep2) and
             IsLeftTurn(sp,bp,ep1)))

def ClosestPtFromPt(state, poly, last_bounce_edge=-1):
    psize = poly.size
    vs = poly.complete_vertex_list
    closest_bounce = 100000000000
    bounce_point = np.array([state[0], state[1]])
    found_coll = False
    bounce_edge = last_bounce_edge
    (x,y,theta) = state

    # check each edge for collision
    for j in range(psize):
        if (j != last_bounce_edge):
            v1, v2 = vs[j], vs[(j+1) % psize]
            try:
                t,u,pt = ShootRay(state, v1, v2)
                
                # Find closest bounce for which t > 0
                pdist = np.linalg.norm(pt-np.array([state[0], state[1]]))
                if (t > 0) and (0 < u) and (u < 1) and (pdist < closest_bounce):
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

# shoot ray from p1 toward p2 in poly, return closest intersect point
# will not return point on 'last_bounce_edge'
def ClosestPtAlongRay(p1,p2,poly,last_bounce_edge=-1):
    '''
    shoot ray from p1 toward p2 in poly, return closest intersect point will not return point on 'last_bounce_edge'
    '''
    vs = poly.complete_vertex_list
    closest_bounce = 100000000000
    bounce_point = np.array([0.0, 0.0])
    found_coll = False
    bounce_edge = last_bounce_edge
    for poly_vxs in poly.vertex_list_per_poly:
        psize = len(poly_vxs)
        # check each edge for collision
        for j in range(psize):
            v1, v2 = poly_vxs[j][1], poly_vxs[(j+1) % psize][1]
            try:
                t,u,pt = ShootRayFromVect(p1, p2, v1, v2)
                # Find closest bounce for which t > 0
                pdist = np.linalg.norm(pt-p2)
                if (t > 0) and (u > 0) and (u < 1) and (pdist < closest_bounce):
                    found_coll = True
                    bounce_point = pt
                    closest_bounce = pdist
                    bounce_edge = poly_vxs[j][0]
            # bounce was parallel to edge j
            # check if it's parallel and overlapping for degenerate polys
            except:
                pass
    if found_coll:
        return bounce_point, bounce_edge
    else:
        return False

def ShootRaysFromReflex(curr_poly_vxs, vertex_list_per_poly, j, poly):
    ''' shoot two rays from each reflex vertex, along incident edges
        group by edge so we can insert them in correct order
    '''
    curr_poly_size = len(curr_poly_vxs)
    p2 = curr_poly_vxs[j][1]
    p1_ccw = curr_poly_vxs[(j+1) % curr_poly_size][1]
    p1_cw = curr_poly_vxs[(j-1) % curr_poly_size][1]

    int_pts = []

    # ClosestPtAlongRay returns False if we shoot ray into existing vertex
    result = ClosestPtAlongRay(p1_ccw, p2, poly)
    if result:
        pt, k = result
        if not VertexExists(pt, curr_poly_vxs):
            int_pts.append((pt,k, (j+1) % curr_poly_size))

    result = ClosestPtAlongRay(p1_cw, p2, poly)
    if result:
        pt, k = result
        if not VertexExists(pt, curr_poly_vxs):
            int_pts.append((pt, k, (j-1) % curr_poly_size))

    return int_pts

def ShootRaysToReflexFromVerts(curr_poly_vxs, curr_poly_index, vertex_list_per_poly, poly, j):
    ''' shoot ray from visible vertices through a specific reflex vert
        j is the index of the reflex vertex
        Poly -> Int -> [(Point, Int)]
    '''
    r_v = curr_poly_vxs[j][1]
    prev_r_v = curr_poly_vxs[j-1][1]
    next_r_v = curr_poly_vxs[(j+1)%len(curr_poly_vxs)][1]
    bv1 = r_v - prev_r_v
    bv2 = r_v - next_r_v

    pts = []
    visible_verts = visibleVertices(curr_poly_vxs, vertex_list_per_poly, poly, j)

    # only ray shoot from non-adjacent vertices
    # previous and next neighbors always visible
    for index, curr_vis_vx_set in enumerate(visible_verts):
        vis_poly_vx = vertex_list_per_poly[index]
        vis_poly_size = len(vis_poly_vx)
        for v in curr_vis_vx_set:
            if (curr_poly_index == index and v != (j-1)%vis_poly_size) and (v != (j+1)%vis_poly_size) or (curr_poly_index != index):
                mid_vector = vis_poly_vx[v][1] - r_v
                if (vector_bwt_two_vector(mid_vector, bv1, bv2)):
                    continue
                res = ClosestPtAlongRay(vis_poly_vx[v][1], r_v, poly)
                if res:
                    pt, k = res
                    if not VertexExists(pt, curr_poly_vxs):
                        pts.append((pt, k, v))
    return pts

def vector_bwt_two_vector(input_vector, bv1, bv2):
    if np.cross(bv1, bv2) >= 0:
        return (np.cross(bv1, input_vector) >= 0 and np.cross(bv2, input_vector) <= 0)
    else:
        return (np.cross(bv1, input_vector) < 0 and np.cross(bv2, input_vector) > 0)

def IsInPolyNoHoles(p, vs):
    ''' test if point p is in poly using crossing number.
    Note: this does not work with holes
    '''
    intersects = 0
    theta = np.random.rand()*2*np.pi
    state=(p[0],p[1],theta)
    psize = len(vs)
    for j in range(psize):
        v1, v2 = vs[j], vs[(j+1) % psize]
        try:
            t, u, pt = ShootRay(state, v1, v2)
            if t>0 and (0 < u) and (u < 1):
                intersects += 1
        except:
            pass

    return not (intersects%2 == 0)

def IsInPoly(p, poly):
    ''' test if point p is in poly using crossing number.
    Checks first if p in in outer boundary, then checks that p is not in a hole.
    '''
    outer_vs = [v for i,v in poly.vertex_list_per_poly[0]]
    isInOuter = IsInPolyNoHoles(p, outer_vs)
    if not isInOuter:
        return False
    else:
        holes = poly.holes
        inHoles = []
        for hole in holes:
            vs = [v for i,v in hole]
            vs.reverse()
            inHoles.append(IsInPolyNoHoles(p, vs))
        return (not any(inHoles))


