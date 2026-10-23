from settings import *
from simple_polygon import Simple_Polygon
from helper.shoot_ray_helper import *

class Partial_Local_Sequence(object):
    ''' Partial local sequence is a sequence of points on the boundary of P associated with a vertex v0, and is constructed by shooting a ray through a reflex vertex v0 from every visible vertex and keeping the resulting sequence of intersections with the boundary of P. The partial local sequence for a convex vertex is empty.

    Attributes
    ----------
    polygon : numpy.array
        The polygon that the partial local sequence is built upon
    inserted_polygon : numpy.array
        The polygon with all new vertices inserted 
    sequence_info : list
        A 2D array stores the partial local sequence for each vertex in polygon in rows
    '''
    def compute_sequence(self, poly):
        ''' Compute the partial local sequence for all vertices of the polygon
        '''
        rvs = poly.reflex_vertices
        sequence_info = []
        def compute_seq_for_single_polygon(curr_polygon_vx, curr_polygon_index):
            curr_seq_info = []
            for i in range(len(curr_polygon_vx)):
                if not curr_polygon_vx[i][0] in rvs:
                    curr_seq_info.append([])
                    continue
                r_children = ShootRaysFromReflex(curr_polygon_vx, poly.vertex_list_per_poly, i, poly)
                transition_pts = ShootRaysToReflexFromVerts(curr_polygon_vx, curr_polygon_index, poly.vertex_list_per_poly, poly, i)
                # transition_pts = []
                transition_pts.extend(r_children)
                curr_seq_info.append(transition_pts)
            return curr_seq_info
        sequence_info.extend(compute_seq_for_single_polygon(poly.outer_boundary_vertices, 0))
        for index, hole in enumerate(poly.holes):
            sequence_info.extend(compute_seq_for_single_polygon(hole, index+1))
        return sequence_info


    def compute_inserted_polygon(self, polygon, sequence_info):
        ''' Compute the polygon with all new vertices from the partial local sequence inserted
        '''
        t_pts_grouped = {i:[] for i in range(polygon.size)}
        for i in range(len(sequence_info)):
            for (pt, edge, _) in sequence_info[i]:
                t_pts_grouped[edge].append(pt)
        
        def get_new_vertices_for_single_polygon(curr_polygon_vx):
            # sort transition points along edge
            new_poly_grouped = {}
            for i in range(len(curr_polygon_vx)):
                new_poly_grouped[i] = sort_by_distance(curr_polygon_vx[i][1], np.array(t_pts_grouped[curr_polygon_vx[i][0]]))
                if DEBUG:
                    print('Inserted verts', new_poly_grouped[i], 'on edge',i)

            inserted_polygon = []
            new_vertices = []
            for i in range(len(curr_polygon_vx)):
                new_vertices.extend(new_poly_grouped[i])
            return new_vertices
        new_vertices = get_new_vertices_for_single_polygon(polygon.outer_boundary_vertices)
        new_holes = []
        for hole in polygon.holes:
            new_holes.append(get_new_vertices_for_single_polygon(hole))
        inserted_polygon = Simple_Polygon('inserted_' + polygon.name, np.array(new_vertices), new_holes)
        return inserted_polygon

    # initialized with Simple_Polygon instance
    def __init__(self, polygon):
        self.polygon = polygon
        self.sequence_info = self.compute_sequence(self.polygon)
        self.inserted_polygon = self.compute_inserted_polygon(self.polygon, self.sequence_info) 
