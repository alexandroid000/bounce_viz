#! /usr/bin/env python3
import sys
sys.path.append('../src')
from helper.visualization_helper import *
from environments.maps import simple_bit
from simple_polygon import Simple_Polygon
from bounce_visibility_diagram import Bounce_Visibility_Diagram
from bounce_graph import Bounce_Graph
from navigation import Navigation

if __name__ == '__main__':
    poly = Simple_Polygon("sp",simple_bit[0])
    poly_vx = poly.complete_vertex_list
    poly_name = "sp"
    pls = Partial_Local_Sequence(poly)
    bvd = Bounce_Visibility_Diagram(pls)
    bounce_graph = Bounce_Graph(bvd)

#    pls2 = Partial_Local_Sequence(poly2)
#    origin = poly2.vertices[10]
#    sequence = pls2.sequence_info[10]
    
    start = (0.1, 0.15)
    goal = (0.51, 0.66)
    nav_task = Navigation(start, goal, bounce_graph)
    path = nav_task.navigate()
    print(path)

    visualize_all_partial_order_sequence(pls.polygon.vertices, pls.inserted_polygon.vertices, pls.sequence_info)
    visualize_polygon(poly_vx, poly_name)
    visualize_bounce_visibility_diagram(bvd, hline = 1.4707)
    visualize_polygon(pls.inserted_polygon.vertices, 'inserted_'+poly_name)
    visualize_partial_local_sequence_for_one_vx(poly2.vertices, origin, sequence)
    visualize_graph(bounce_graph.visibility_graph, "bounce_visibility_graph")
    visualize_graph(bounce_graph.safe_action_graph, "bounce_safe_action_graph")
    # intervals = PropagatePath(pls.inserted_polygon.vertices, path, start)
    # VizPath(p1, intervals)
