from numpy.linalg import norm
import numpy as np
import matplotlib.pyplot as plt
import os.path as osp

import settings
from helper.geometry_helper import IsRightTurn
from helper.general_position_helper import IsThreePointsOnLine

class Simple_Polygon(object):
    ''' Partial local sequence is a sequence of points on the boundary of P associated with a vertex v0, and is constructed by shooting a ray through a reflex vertex v0 from every visible vertex and keeping the resulting sequence of intersections with the boundary of P. The partial local sequence for a convex vertex is empty.

    Attributes
    ----------
    name : str
        name of the environment; used for saving results
    outer_boundary_vertices : numpy.array
        A list of coordinates of vertices in outer most boundary of the polygon in counter clockwise order.
    holes : list of numpy.array 
        A list of list of coordinates, each list contains vertices in a hole in the polygon in clockwise order
    complete_vertex_list : list of (int, np.array)
        Assigned each vertex an index and store the the index and coordinates in a list. The tuples are stored according to the index increasing order
    vertex_list_per_poly : list of list of np.array
        Store the vertex list for each close boundary in the polygon into a list
    size : int
        the total number of vertices in the polygon
    reflex_vertices : list of numpy.array
        list of reflex vertices in the polygon using the index assigned in complete_vertex_list
    unit_interval_mapping : (list of floats, list of list of floats)
        see compute_unit_interval_mapping

        
    '''
    def compute_unit_interval_mapping(self):
        ''' Map the vertices of the polygon on to a given interval (length specified in settings) where both the start and the end of the interval identify with the first vertex in the polygon
        '''
        size = len(self.outer_boundary_vertices)
        edge_len = [norm(self.outer_boundary_vertices[i][1]-self.outer_boundary_vertices[(i+1)%size][1]) for i in range(size)]
        outer_acc_edge_len = np.array([sum(edge_len[:i]) for i in range(len(edge_len)+1)])
        holes_acc_edge_len = []
        total_edge_len = outer_acc_edge_len[-1]
        for hole in self.holes:
            size = len(hole)
            edge_len = [norm(hole[i][1]-hole[(i+1)%size][1]) for i in range(size)]
            acc_edge_len = np.array([sum(edge_len[:i]) for i in range(len(edge_len)+1)])
            holes_acc_edge_len.append(acc_edge_len)
            total_edge_len += acc_edge_len[-1]

        outer_acc_edge_len /= total_edge_len/settings.unit_interval_len
        holes_acc_edge_len = [acc_edge_len / (total_edge_len/settings.unit_interval_len) for acc_edge_len in holes_acc_edge_len]
        return outer_acc_edge_len, holes_acc_edge_len

    def reflex_verts_single_boundary(self, vs):
        ''' return indices of all reflex vertices in poly
        '''
        reflex_verts = []
        size = len(vs)
        for j in range(size):
            v1, v2, v3 = vs[(j-1) % size][1], vs[j][1], vs[(j+1) % size][1]
            if IsRightTurn(v1,v2,v3) and not IsThreePointsOnLine(v1,v2,v3):
                reflex_verts.append(vs[j][0])

        return reflex_verts

    def get_reflex_vertices(self):
        reflex_verts = self.reflex_verts_single_boundary(self.outer_boundary_vertices)
        for hole in self.holes:
            reflex_verts.extend(self.reflex_verts_single_boundary(hole))
        return reflex_verts

    def get_hole_size_list(self, start, holes):
        size_list = [start]
        curr_sum = start
        for hole in holes:
            curr_sum += len(hole)
            size_list.append(curr_sum)
        return size_list
    
    def visualize(self):
        ''' Draws polygon with numbered vertices
        '''
        jet = plt.cm.jet
        colors = jet(np.linspace(0, 1, self.size))

        # plot only polygon, no axis or frame
        fig = plt.figure(frameon=False)
        ax = plt.gca()
        ax.axis('off')
        ax.axis('equal')
        outer_vxs = list(list(zip(*self.outer_boundary_vertices))[1])
        wall = np.vstack((outer_vxs, outer_vxs[0]))
        wall = np.vstack((wall, np.array([np.nan, np.nan])))
        for hole in self.holes:
            hole_vxs = list(list(zip(*hole))[1])
            wall = np.vstack((wall, hole_vxs))
            wall = np.vstack((wall, hole_vxs[0]))
            wall = np.vstack((wall, np.array([np.nan, np.nan])))
        plt.plot(wall[:, 0], wall[:, 1], 'black')

        def plot_vertices(vertices):
            for vx in vertices:
                point = vx[1]
                plt.scatter(point[0], point[1], color=colors[vx[0]], s=60)
                plt.annotate(str(vx[0]), (point[0], point[1]), size = 10)

        plot_vertices(self.outer_boundary_vertices)
        for hole in self.holes:
            plot_vertices(hole)
        plt.savefig(osp.join(settings.image_save_folder,self.name+'.png'), dpi = 300, bbox_inches='tight')
        if settings.DEBUG:
            plt.show()

    def __init__(self, name, outer_boundary_vertices, holes = []):
        self.name = name
        N = outer_boundary_vertices.shape[0]
        self.complete_vertex_list = [v for v in outer_boundary_vertices]
        self.outer_boundary_vertices = [(i, np.array(outer_boundary_vertices[i])) for i in range(N)]
        size_list = self.get_hole_size_list(N, holes)
        self.holes = [[(i+size_list[j], np.array(holes[j][i])) for i in range(len(holes[j]))] for j in range(len(holes))]
        self.vertex_list_per_poly = [self.outer_boundary_vertices]
        for hole in self.holes:
            self.complete_vertex_list.extend([vx[1] for vx in hole])
            self.vertex_list_per_poly.append(hole)
        self.size = len(self.complete_vertex_list)
        self.reflex_vertices = self.get_reflex_vertices()
        self.unit_interval_mapping = self.compute_unit_interval_mapping()
