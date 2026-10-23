import numpy as np
import matplotlib.pyplot as plt
import os.path as osp

from helper.visibility_helper import *
from helper.geometry_helper import *
import settings


class Bounce_Visibility_Diagram(object):
    ''' Bounce visibility diagram is a representation of the vertex-edge visibility structure of a partitioned polygon
    
    Attributes
    ----------
    resolution : int
        the discretization constant for the angle function
    partial_local_sequence : :obj:`Partial_Local_Sequence`
        The partial local sequence for the polygon
    visible_vx_set_for_edges : list of list of list of list
        format: [poly1: [e1:[viz vx poly1:[], viz vx poly2:[]], e2:[viz vx poly1:[], viz vx poly2:[]]], poly2: [e1:[viz vx poly1:[], viz vx poly2:[]], e2:[viz vx poly1:[], viz vx poly2:[]]]]
    visible_angle_info : list of numpy.array
        The angle functions for each edge in each closed polygonal boundary
    '''
    def compute_visible_angle_info(self, pls, visible_vx_set_for_edges, resolution):
        ''' Compute the visible angle info for the input polygon
            For each visible vertex, we need to calculate:
             (1) the view angle at the current vertex w.r.t the current edge and
             (2) the view angle at the sample points on the edge and the next vertex w.r.t the current edge
             (3) insert a np.nan for discontinuity otherwise matplotlib will try to connect them together

        Parameters
        ----------
        pls: :obj:`Partial_Local_Sequence`
        resolution:int
            the number of sample points on each edge

        Returns
        -------
        np.array
            The visible angle info for the input polygon
        '''
        def compute_visible_angle_info_single_polygon(curr_poly, viz_set_for_curr_poly, complete_vertex_list, total_size):
            psize =  len(curr_poly)
            visible_angle_info = np.nan*np.ones((total_size, resolution*psize))
            for i, edge_visible_vxs in enumerate(viz_set_for_curr_poly):
                curr_p, next_p = curr_poly[i][1], curr_poly[(i+1)%psize][1]
                for poly_index, viz_poly_vx in enumerate(edge_visible_vxs):
                    for vx_index in viz_poly_vx:
                        interest_p = complete_vertex_list[vx_index]
                        visible_angle_info[vx_index][resolution*i] = get_angle_from_three_pts(next_p, interest_p, curr_p)
                        visible_angle_info[vx_index][resolution*i+resolution-1] = np.nan
                        # the if case is for when we are considering the 'next vertex'
                        if (interest_p == next_p).all():
                            for stride in range(resolution)[1:-1]:
                                visible_angle_info[vx_index][resolution*i+stride] = visible_angle_info[vx_index][resolution*i]
                            continue
                        for stride in range(resolution)[1:-1]:
                            # this is the point on the same line of the current edge but in front of the next vertex.
                            # Use this so that we don't get zero length vector from the next vertex perspective
                            extended_point = curr_p+2*(next_p-curr_p)
                            ratio = 1./(resolution-2)*stride
                            origin = ratio*next_p+(1-ratio)*curr_p
                            visible_angle_info[vx_index][resolution*i+stride] = get_angle_from_three_pts(interest_p, extended_point, origin)
            return visible_angle_info

        visible_angle_info = [compute_visible_angle_info_single_polygon(pls.inserted_polygon.outer_boundary_vertices, visible_vx_set_for_edges[0], pls.inserted_polygon.complete_vertex_list, pls.inserted_polygon.size)]

        for index, hole in enumerate(pls.inserted_polygon.holes):
            visible_angle_info.append(compute_visible_angle_info_single_polygon(hole, visible_vx_set_for_edges[index+1], pls.inserted_polygon.complete_vertex_list, pls.inserted_polygon.size))
        return visible_angle_info
    
    def visualize(self, hline = None):
        ''' Save the link diagram for a given polygon to image

        Parameters
        ----------
        bvd: :obj:`Bounce_Visibility_Diagram`
            The link diagram computed for the given polygon
        hline:float
            The angle of fix theta bouncing
        '''
        def set_y_labels(ax):
            ax.set_yticks([0., np.pi/6., np.pi/3., np.pi/2., 2*np.pi/3., 5*np.pi/6., np.pi])
            ax.set_yticklabels(["$0$", r"$\frac{1}{6}\pi$", r"$\frac{1}{3}\pi$", r"$\frac{1}{2}\pi$", r"$\frac{2}{3}\pi$", r"$\frac{5}{6}\pi$", r"$\pi$"])
            plt.xlabel(r"vertices on $\partial P'$", fontsize=15)
            plt.ylabel(r"bounce angle $\theta$", fontsize=15)

        def get_visualize_info_for_single_poly(curr_poly, acc_edge_len, curr_visible_angle_info, total_size):
            psize = len(curr_poly)
            x = []
            for i in range(psize):
                x.extend(list(np.linspace(acc_edge_len[i], acc_edge_len[i+1], self.resolution-1)))
                x.extend([acc_edge_len[i+1]])
                plt.axvline(x=acc_edge_len[i], linestyle='--', linewidth=0.8)
            plt.axvline(x = acc_edge_len[-1], linestyle='--', linewidth=0.8)

            for i in range(total_size):
                plt.plot(x, curr_visible_angle_info[i], label= '{}'.format(i), alpha=0.7, color = colors[i])
            if hline != None:
                plt.axhline(y = hline)
            leg = plt.legend(loc=9, bbox_to_anchor=(0.5, -0.15), ncol=5)
            x_labels = [curr_poly[i][0] for i in range(psize)]
            x_labels.append(curr_poly[0][0])
            return x_labels

        def visualize_single_close_boundary_bvd(curr_poly, acc_edge_len, index, curr_visible_angle_info, total_size):
            ''' Visualize the bvd for vertices on a single close boundary
            '''
            fig, ax = plt.subplots()
            x_labels = get_visualize_info_for_single_poly(curr_poly, acc_edge_len, curr_visible_angle_info, total_size)
            ax.set_xticks(acc_edge_len)
            ax.set_xticklabels(x_labels)
            set_y_labels(ax)
            plt.savefig(osp.join(settings.image_save_folder, 'bvd_'+self.partial_local_sequence.inserted_polygon.name+'_'+str(index)+'.png'), dpi = 300, bbox_inches='tight')
            if DEBUG:
                plt.show()

        def visualize_all_bvd():
            ''' Visualize all vertices in one bvd
            '''
            fig, ax = plt.subplots()
            total_acc_edge_len = []
            gap_bwt_polys = settings.unit_interval_len/15.
            total_x_labels = get_visualize_info_for_single_poly(outer_boundary_vertices, unit_interval_mapping[0], self.visible_angle_info[0], total_size)
            total_acc_edge_len.extend(unit_interval_mapping[0])
            start_x = unit_interval_mapping[0][-1] + gap_bwt_polys
            
            for index, hole in enumerate(holes):
                acc_edge_len = [x+start_x for x in unit_interval_mapping[1][index]]
                total_x_labels.extend(get_visualize_info_for_single_poly(hole, acc_edge_len, self.visible_angle_info[index+1], self.partial_local_sequence.inserted_polygon.size))
                total_acc_edge_len.extend(acc_edge_len)
                start_x = acc_edge_len[-1] + gap_bwt_polys

            ax.set_xticks(total_acc_edge_len)
            ax.set_xticklabels(total_x_labels, fontsize = 9)
            set_y_labels(ax)
            plt.savefig(osp.join(settings.image_save_folder, 'bvd_all_'+self.partial_local_sequence.inserted_polygon.name+'.png'), dpi = 300, bbox_inches='tight')
            if DEBUG:
                plt.show()

        total_size = self.partial_local_sequence.inserted_polygon.size
        jet = plt.cm.jet
        colors = jet(np.linspace(0, 1, total_size))
        unit_interval_mapping = self.partial_local_sequence.inserted_polygon.unit_interval_mapping
        outer_boundary_vertices = self.partial_local_sequence.inserted_polygon.outer_boundary_vertices
        holes = self.partial_local_sequence.inserted_polygon.holes

        visualize_single_close_boundary_bvd(outer_boundary_vertices, unit_interval_mapping[0], 0, self.visible_angle_info[0], total_size)
        for index, hole in enumerate(holes):
            visualize_single_close_boundary_bvd(hole, unit_interval_mapping[1][index], index+1, self.visible_angle_info[index+1], total_size)
        visualize_all_bvd()

    def __init__(self, pls):
        self.resolution = 15
        self.partial_local_sequence = pls
        self.visible_vx_set_for_edges = get_all_edge_visible_vertices(pls.inserted_polygon)
        self.visible_angle_info = self.compute_visible_angle_info(pls, self.visible_vx_set_for_edges, self.resolution)
