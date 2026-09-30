from environments.maps import *
from settings import *
from helper.geometry_helper import *
from helper.shoot_ray_helper import ShootRaysFromReflex, ShootRaysToReflexFromVerts
from partial_local_sequence import Partial_Local_Sequence

import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
import os.path as osp
import networkx as nx

def visualize_all_partial_order_sequence(poly_vx, inserted_poly_vx, sequence_info):
    ''' Draws all transition points in poly and saves to viz_test.pdf
    '''
    # Plot the outer boundary with black color
    plt.figure()
    wall = np.vstack((poly_vx, poly_vx[0]))
    ax = plt.gca()
    ax.axis('off')
    plt.plot(wall[:, 0], wall[:, 1], 'black')
    plt.plot(inserted_poly_vx[:, 0], inserted_poly_vx[:, 1], 'ro')
    for index, seq in enumerate(sequence_info):
        for (pt,k, _) in seq:
            plt.plot([poly_vx[index][0], pt[0]], [poly_vx[index][1], pt[1]], 'green')
    plt.savefig(osp.join(image_save_folder,'viz_test.pdf'))
    if DEBUG:
        plt.show()

def visualize_polygon(poly, fname):
    ''' Draws polygon with numbered vertices
    '''
    psize = poly.shape[0]
    jet = plt.cm.jet
    colors = jet(np.linspace(0, 1, psize))

    # plot only polygon, no axis or frame
    fig = plt.figure(frameon=False)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')
    wall = np.vstack((poly, poly[0]))
    plt.plot(wall[:, 0], wall[:, 1], 'black')
    for i in range(psize):
        point = poly[i]
        plt.scatter(point[0], point[1], color=colors[i], s=100)
        plt.annotate(str(i), (point[0]+5, point[1]+5), size = 20)
    plt.axis('equal')

    plt.savefig(osp.join(image_save_folder,fname+'.png'), dpi = 300, bbox_inches='tight')
    if DEBUG:
        plt.show()


def visualize_subset_polygon(poly, fname, subset):
    ''' Draws polygon with numbered vertices
        Colors selected subset of edges red
    '''
    psize = poly.shape[0]
    jet = plt.cm.jet
    colors = jet(np.linspace(0, 1, psize))

    # plot only polygon, no axis or frame
    fig = plt.figure(frameon=False)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')
    wall = np.vstack((poly, poly[0]))
    plt.plot(wall[:, 0], wall[:, 1], 'black')
    for i in subset:
        edge = np.array([wall[i], wall[i+1]])
        plt.plot(edge[:, 0], edge[:, 1], 'red', linewidth=2)
    for i in range(psize):
        point = poly[i]
        plt.scatter(point[0], point[1], color='black', s=20)
    plt.axis('equal')

    print("Saving", image_save_folder+'/'+fname+'.png')
    plt.savefig(osp.join(image_save_folder,fname+'.png'), dpi = 300, bbox_inches='tight')
    if DEBUG:
        plt.show()



def visualize_partial_local_sequence_for_one_vx(poly_vx, origin, sequence, fname = 'partial_local_sequence'):
    ''' Draw the partial local sequence for a given vertex
    '''
    plt.figure()
    ax = plt.gca()
    ax.axis('off')
    wall = np.vstack((poly_vx, poly_vx[0]))
    plt.plot(wall[:, 0], wall[:, 1], 'black')
    psize = poly_vx.shape[0]
    for i, pt in enumerate(poly_vx):
        plt.plot(pt[0], pt[1], 'ko', markersize = 5)
        plt.annotate("{}".format(i), xy = (pt[0], pt[1]+7), size =8, color = 'red')

    for (pt,k, v) in sequence:
        plt.plot([origin[0], pt[0]], [origin[1], pt[1]], 'b--',  markersize = 4)
        plt.plot(pt[0], pt[1], 'bo', markersize = 5)
        plt.annotate("{}'".format(v), xy = (pt[0]-8, pt[1]-20),  size =8, color = 'red')

    plt.plot(origin[0], origin[1], 'bo', markersize = 5)
    plt.savefig(osp.join(image_save_folder, fname+'.png'), dpi = 300)
    if DEBUG:
        plt.show()

def visualize_bounce_visibility_diagram(bvd, hline = None, fname = 'bvd.png'):
    ''' Save the link diagram for a given polygon to image

    Parameters
    ----------
    bvd: :obj:`Bounce_Visibility_Diagram`
        The link diagram computed for the given polygon
    hline:float
        The angle of fix theta bouncing
    fname:string
        The output file name for the link diagram
    '''
    i_poly_vx = np.array(bvd.partial_local_sequence.inserted_polygon.complete_vertex_list)
    psize = i_poly_vx.shape[0]
    acc_edge_len = bvd.partial_local_sequence.inserted_polygon.unit_interval_mapping
    jet = plt.cm.jet
    colors = jet(np.linspace(0, 1, psize))
    fig, ax = plt.subplots()
    x = []
    for i in range(psize):
        x.extend(list(np.linspace(acc_edge_len[i], acc_edge_len[i+1], bvd.resolution-1)))
        x.extend([acc_edge_len[i+1]])
    for i in range(psize):
        plt.plot(x, bvd.visible_angle_info[i], label= '{}'.format(i), alpha=0.7, color = colors[i])
        plt.axvline(x=acc_edge_len[i], linestyle='--')
    plt.axvline(x = acc_edge_len[-1], linestyle='--')
    if hline != None:
        plt.axhline(y = hline)
    leg = plt.legend(loc=9, bbox_to_anchor=(0.5, -0.15), ncol=5)
    ax.set_xticks(acc_edge_len)
    x_labels = list(range(psize))
    x_labels.append(0)
    ax.set_xticklabels(x_labels)
    ax.set_yticks([0., np.pi/6., np.pi/3., np.pi/2., 2*np.pi/3., 5*np.pi/6., np.pi])
    ax.set_yticklabels(["$0$", r"$\frac{1}{6}\pi$", r"$\frac{1}{3}\pi$", r"$\frac{1}{2}\pi$", r"$\frac{2}{3}\pi$", r"$\frac{5}{6}\pi$", r"$\pi$"])
    plt.xlabel(r"vertices on $\partial P'$", fontsize=15)
    plt.ylabel(r"bounce angle $\theta$", fontsize=15)
    plt.savefig(osp.join(image_save_folder,fname), bbox_inches='tight', dpi = 300)
    if DEBUG:
        plt.show()

def visualize_graph(G, fname = 'graph'):
    ''' Draw graph in circular shape
    '''
    plt.clf()
    pos = nx.circular_layout(G)
    labels = nx.get_edge_attributes(G,'weight')
    new_pos = dict()
    count = 0
    for i in pos:
        new_pos[count] = pos[i]
        count += 1
    nx.draw_networkx(G, with_labels=True,
        pos = new_pos,
        node_color='white',
        width=1.2,
        node_size = 600, 
        font_size=14)
    plt.axis('off')
    ax = plt.gca()
    ax.collections[0].set_edgecolor('black') 
    plt.savefig(osp.join(image_save_folder,fname+'.png'), bbox_inches='tight', dpi = 300)

# def VizPath(poly, intervals):
#     psize = len(poly)
#     jet = plt.cm.jet
#     colors = jet(np.linspace(0, 1, psize))

#     # plot only polygon, no axis or frame
#     fig = plt.figure(frameon=False)
#     ax = fig.add_axes([0, 0, 1, 1])
#     ax.axis('off')

#     wall_x = [x for (x,y) in poly]
#     wall_x.append(wall_x[0])
#     wall_y = [y for (x,y) in poly]
#     wall_y.append(wall_y[0])
#     plt.plot(wall_x, wall_y, 'black')
#     for i in range(psize):
#         point = poly[i]
#         plt.scatter(point[0], point[1], color=colors[i])
#         plt.annotate(str(i), (point[0]+10, point[1]+10), size = 'small')
#     plt.axis('equal')

#     print(intervals)

#     interval_ts = list(zip(intervals, intervals[1:]))
#     for i in interval_ts:
#         print(i)
#     p1, p2 = list(interval_ts)[0][0]
#     plt.plot([p1[0],p2[0]], [p1[1], p2[1]], 'red',linewidth=5)

#     for ((pt1, pt2), (new_pt1, new_pt2)) in interval_ts:
#         print('plotting', ((pt1, pt2), (new_pt2, new_pt1)))
#         plt.plot([new_pt1[0],new_pt2[0]], [new_pt1[1], new_pt2[1]], 'red',linewidth=5)
#         plt.plot([pt1[0],new_pt2[0]], [pt1[1], new_pt2[1]], 'green')
#         plt.plot([pt2[0],new_pt1[0]], [pt2[1], new_pt1[1]], 'green')

#     plt.savefig(osp.join(image_save_folder,'path.png'), dpi = 300)
#     if DEBUG:
#         plt.show()
