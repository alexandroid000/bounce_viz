def get_edges_intersecting_interval(poly_unit_map, interval):
    ''' Find edges on the polygon that intersect with the interval
        Interval from (0,1) parameterization of polygon
        TODO: Only works for polygons with no obstacles/holes right now
        does not currently sort output, no need to for now
    '''
    (s1, s2) = interval
    vxs = []
    if len(poly_unit_map[1]) == 0:
        map = poly_unit_map[0]

        first_i = 0

        # first element of map always 0
        # last element of map always 1
        for i in range(len(map)-1):
            if map[i]<s1:
                first_i = i
                continue
            if map[i]>s2:
                break
            vxs.append(i)

        # add index of first edge intersecting with interval
        vxs.append(first_i)
    else:
        raise ValueError("Polygon contains holes")

    return vxs

def get_transition_over_path(path, safe_action_graph):
    ''' Get the transition with angle range information from a given path
    '''
    transitions = []
    path = list(path)[0]
    for (i,j) in zip(path, path[1:]):
        ang_range = safe_action_graph.get_edge_data(i,j)
        transitions.append((i,j,ang_range))
    return transitions

class Navigation(object):
    ''' Description of a navigation task for a given polygon
    Attributes
    ----------
    start_position : (float, float)
        The range of starting position in the unit interval mapping of the polygon
    end_position : (float, float)
        The range of ending position in the unit interval mapping of the polygon
    bvg : :obj:`Bounce_Visibility_Graph`
    bvd: :obj:`Bounce_Visibility_Diagram`
    pls: :obj:`Partial_Local_Sequence`
    polygon: :obj:`Simple_Polygon`
    inserted_polygon: :obj:`Simple_Polygon`
    '''
    def navigate(self):
        ''' Executing the navigation task with a given strategy
        '''
        paths = []
        start_nodes = get_edges_intersecting_interval(self.inserted_polygon.unit_interval_mapping, self.start_position)
        end_nodes = get_edges_intersecting_interval(self.inserted_polygon.unit_interval_mapping, self.end_position)
        if len(start_nodes) != 1 or len(end_nodes) != 1:
            print("Start nodes:", start_nodes)
            print("Goal nodes:", end_nodes)
            print("Vertex mapping:", self.inserted_polygon.unit_interval_mapping)
            raise ValueError("Please choose start/goal intervals contained within one edge")
        for g in end_nodes:
            for s in start_nodes:
                paths.append(self.bvg.get_shortest_path(s, g, 'safe'))
        return get_transition_over_path(paths[0], self.bvg.safe_action_graph)

    def __init__(self, start_position, end_position, bvg):
        self.start_position = start_position 
        self.end_position = end_position  
        self.bvg = bvg
        self.bvd = bvg.bounce_visibility_diagram
        self.pls = self.bvd.partial_local_sequence
        self.polygon = self.pls.polygon
        self.inserted_polygon = self.pls.inserted_polygon
