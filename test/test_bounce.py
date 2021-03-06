# #! /usr/bin/env python
import sys
sys.path.append("./src")
import unittest
from bounce_graph import *
from simple_polygon import *
from bounce_visibility_diagram import *
from partial_local_sequence import *
from helper.visibility_helper import *
from helper.bounce_graph_helper import *
from navigation import *
from maps import *
from classify import *
from math import pi
import numpy as np

def compare(a,b):
    if not isinstance(b, type(a)):
        raise TypeError("The two input object must be of same type")
    else:
        if isinstance(a, (list, tuple)):
            for aa, bb in zip(a,b):
                if not compare(aa, bb):
                    return False
            return True
        else:
            return np.allclose(a,b)

class TestGeomUtils(unittest.TestCase):

    def setUp(self):
         self.origin = np.array([0.0, 0.0])
         # square on xy axis
         self.p1 = np.array([ 10.0,   0.0])
         self.p2 = np.array([-10.0,   0.0])
         self.p3 = np.array([  0.0,  10.0])
         self.p4 = np.array([  0.0, -10.0])
         # right triangle
         self.t1 = np.array([6.0, 0.0])
         self.t2 = np.array([3.0, 4.0])
         self.t3 = np.array([3.0,0.0])
         self.x3 = [self.origin, self.t3]
         self.y4 = [self.t3,self.t2]
         # edges
         self.e1 = [self.origin, np.array([1.0, 0.0])]
         self.e2 = [np.array([2.0,1.0]), np.array([2.0, 2.0])]
         self.maxDiff = None

    def test_left(self):
         v1 = IsLeftTurn(self.p1, self.p2, self.p4)
         v2 = IsLeftTurn(self.p2, self.p1, self.p3)
         self.assertTrue(v1 and v2)

    def test_right(self):
         v1 = IsRightTurn(self.p1, self.p2, self.p3)
         v2 = IsRightTurn(self.p2, self.p1, self.p4)
         self.assertTrue(v1 and v2)

    def test_contains(self):
         self.assertTrue(IsInPoly((0.0,0.0), Simple_Polygon("sq",square[0])))

    def test_contains_w_hole(self):
         self.assertTrue(IsInPoly((5.0,3.0), Simple_Polygon("sqh",simple_holes[0], simple_holes[1])))

         self.assertFalse(IsInPoly((2.0,1.5), Simple_Polygon("sqh",simple_holes[0], simple_holes[1])))

         self.assertFalse(IsInPoly((8.0,8.5), Simple_Polygon("sqh",simple_holes[0], simple_holes[1])))

    def test_notIn(self):
         self.assertFalse(IsInPoly((500.0,0.0), Simple_Polygon("p1",poly1[0])))

    def test_general_pos(self):
        tr = Simple_Polygon("tr",tworooms2[0])
        self.assertFalse(PolyInGeneralPos(tr))
        self.assertTrue(Simple_Polygon("p1",poly1[0]))

    # def ShootRay(state, v1, v2):
    def test_shootRay(self):
         p1 = np.array([-10.0, 10.0])
         p2 = np.array([-10.0, -10.0])
         state = (0.0, 0.0, pi)
         t, u, (x,y) = ShootRay(state, p1, p2)
         self.assertAlmostEqual(u, 0.5)
         self.assertAlmostEqual(x, -10.0)
         self.assertAlmostEqual(y, 0.0)

    def test_shootRayFromVect(self):
         t, u, (x,y) = ShootRayFromVect(self.t1, self.t2, self.origin, self.p3)
         self.assertAlmostEqual(x, 0.0)
         self.assertAlmostEqual(y, 8.0)

    def test_intersect_interior(self):
         poly = Simple_Polygon("p1",poly1[0])
         (x,y), _ = ClosestPtAlongRay(self.origin, self.p1, poly)
         self.assertAlmostEqual(x, 139.1304347826087)
         self.assertAlmostEqual(y, 0.0)

    def test_intersect_thru_vertex(self):
         poly = Simple_Polygon("p1",poly1[0])
         (x,y), _ = ClosestPtAlongRay(self.origin, (150,50), poly)
         self.assertAlmostEqual(x, 198.21428571428572)
         self.assertAlmostEqual(y, 66.07142857142857)

    def test_reflex(self):
         p1 = Simple_Polygon("p1",poly1[0])
         self.assertEqual([3,7,10], p1.reflex_vertices)

#ShootRaysToReflexFromVerts(curr_poly_vxs, curr_poly_index, vertex_list_per_poly, j):
    def test_transition_pts(self):
         p1 = Simple_Polygon("p1",poly1[0])
         outer_vxs = p1.vertex_list_per_poly[0]

         t3 =  ShootRaysToReflexFromVerts(outer_vxs, 0, p1.vertex_list_per_poly,p1, 3)
         t7 =  ShootRaysToReflexFromVerts(outer_vxs, 0, p1.vertex_list_per_poly,p1, 7)
         t9 =  ShootRaysToReflexFromVerts(outer_vxs, 0, p1.vertex_list_per_poly,p1, 9)
         t10 = ShootRaysToReflexFromVerts(outer_vxs, 0, p1.vertex_list_per_poly,p1, 10)
         ts3 = [(np.array([-67.74193548, 182.25806452]), 4, 1), (np.array([109.82142857, 186.60714286]), 1, 5),
         (np.array([-60., 190.]), 4, 10)]
         ts7 = [(np.array([-206.18384401, -199.13649025]), 8, 4), (np.array([-127.36842105, -194.21052632]), 8,
         5)]
         ts10 = [(np.array([-177.77777778, 50.]), 5, 1), (np.array([207.44680851,  11.70212766]), 0, 3),
         (np.array([209.34065934,  40.10989011]), 0, 5), (np.array([205.3602812 ,  56.32688928]), 1, 6),
         (np.array([199.90439771,  63.7667304 ]), 1, 7), (np.array([190.12048193,  77.10843373]), 1,
         8)]

         compare(ts3, t3)
         compare(ts7, t7)
         compare([], t9)
         compare(ts10, t10)

    def test_pls_convex(self):
        pent = Simple_Polygon("pent",np.array([(550,450), (645,519), (609,631), (491,631), (455,519)], dtype=np.float))
        pls = Partial_Local_Sequence(pent)
        poly_prime = pls.inserted_polygon

        compare(pent.complete_vertex_list, poly_prime.complete_vertex_list)

    def test_viz_convex(self):

        pent = Simple_Polygon("pent",np.array([(550,450), (645,519), (609,631), (491,631), (455,519)], dtype=np.float))
        pls = Partial_Local_Sequence(pent)
        bvg = Bounce_Graph(Bounce_Visibility_Diagram(pls))
        viz_edges = bvg.visibility_graph.edges
        expected = [(0, 1), (0, 2), (0, 3), (0, 4),
                    (1, 0), (1, 2), (1, 3), (1, 4),
                    (2, 0), (2, 1), (2, 3), (2, 4),
                    (3, 0), (3, 1), (3, 2), (3, 4),
                    (4, 0), (4, 1), (4, 2), (4, 3)]

    def test_viz_verts(self):
        poly = Simple_Polygon("sb",simple_bit[0])
        poly_prime_vs = Partial_Local_Sequence(poly).inserted_polygon.complete_vertex_list
        poly_prime = Simple_Polygon("sb2", np.array(poly_prime_vs))
        vvs1 = visibleVertices(poly_prime.outer_boundary_vertices, poly_prime.vertex_list_per_poly, poly, 3)
        vvs2 = visibleVertices(poly_prime.outer_boundary_vertices, poly_prime.vertex_list_per_poly, poly, 2)
        self.assertEqual([[2, 4, 5, 6, 7]], vvs1)
        self.assertEqual([[0, 1, 3, 4, 5, 6, 7, 8, 10, 11]], vvs2)

    def test_viz_gp(self):
        poly = Simple_Polygon("tr",tworooms[0])
        poly_prime_vs = Partial_Local_Sequence(poly).inserted_polygon.complete_vertex_list
        poly_prime = Simple_Polygon("gptr", np.array(poly_prime_vs))
        vvs = visibleVertices(poly_prime.outer_boundary_vertices, poly_prime.vertex_list_per_poly, poly, 0)
        self.assertEqual([[1, 2, 11, 12, 13, 14, 15]], vvs)


    def test_viz_simple_gp(self):
        poly = Simple_Polygon("sqgp", square_gp[0])
        vvs = visibleVertices(poly.outer_boundary_vertices, poly.vertex_list_per_poly, poly, 2)
        self.assertEqual([[0, 1, 3, 4]], vvs)


#def visibleVertices(curr_poly_vx, vertex_list_per_poly, j):
    def test_viz_holes(self):
        poly = Simple_Polygon("simple_holes", simple_holes[0], simple_holes[1])
        vvs = visibleVertices(poly.outer_boundary_vertices, poly.vertex_list_per_poly, poly, 2)
        expected = [[1,3],[1,2]]
        self.assertEqual(expected, vvs)

    def test_angle_parallel(self):
        sq = Simple_Polygon("sq",square[0])
        a = angleBound(sq, 1, 3)
        self.assertAlmostEqual(a, 1.57079632679)

    def test_angle_collinear(self):
        p = Simple_Polygon("line", np.array([(0,0),(0,10),(0,20)]))
        a = angleBound(p, 0, 1)
        self.assertAlmostEqual(a, 0.0)

    # use triangle properties to check that angles are calculated properly
    def test_maxmin_angles(self):
        max_a = 63.4349488*(pi/180)
        min_a = 26.5650511*(pi/180)
        min_c, max_c = AnglesBetweenSegs(self.e1, self.e2)
        self.assertAlmostEqual(min_c, min_a)
        self.assertAlmostEqual(max_c, max_a)

    def test_maxmin_adjacent(self):
        min_a = 0.0
        max_a = pi/2
        min_c, max_c = AnglesBetweenSegs(self.x3, self.y4)
        self.assertAlmostEqual(min_c, min_a)
        self.assertAlmostEqual(max_c, max_a)

    def test_simple_domain(self):
        e1 = np.array([(0.0, 0.0), (2.0, 0.0)])
        e2 = np.array([(3.0, 1.0), (1.0, 1.0)])
        theta = np.pi/2
        subset = findDomain(e1, e2, theta)
        compare(subset, np.array([[1., 0.], [2., 0.]]))

    def test_weirder_domain(self):
        e1 = np.array([(1.0, 1.0), (0.0, 0.0)])
        e2 = np.array([(3.0, -1.0), (3.0, 3.0)])
        theta1 = 3.*np.pi/4.
        theta2 = np.pi/2.
        subset1 = findDomain(e1, e2, theta1)
        subset2 = findDomain(e2, e1, theta2)
        compare(subset1, np.array([[0., 0.], [1., 1.]]))
        compare(subset2, np.array([[3., 0.], [3., 1.]]))

    def test_miss(self):
        e1 = np.array([(0.0, 0.0), (1.0, 0.0)])
        e2 = np.array([(2.0, 1.0), (3.0, 1.0)])
        theta = np.pi/2.
        subset = findDomain(e1, e2, theta)
        self.assertTrue(subset.size == 0)

    def test_indexing(self):
        poly = Simple_Polygon("sh", simple_holes[0], simple_holes[1])
        self.assertEqual(next_v(0, poly), 1)
        self.assertEqual(next_v(3, poly), 0)
        self.assertEqual(next_v(4, poly), 5)
        self.assertEqual(next_v(6, poly), 4)
        poly2 = Simple_Polygon("pent", pent[0])
        self.assertEqual(next_v(4, poly2), 0)


    def test_contraction(self):
        e1 = np.array([[20.,7.],[10.,0.]])
        e2 = np.array([[10.,0.],[0.,7.]])
        e3 = np.array([[16.,18.],[20.,7.]])
        theta1 = 5.*np.pi/6.
        theta2 = np.pi/6.
        theta3 = np.pi/3.
        val1, c1 = isContraction(e1, e2, theta1)
        val2, c2 = isContraction(e1, e3, theta2)
        val3, c3 = isContraction(e1, e3, theta3)
        self.assertTrue(val1)
        self.assertTrue(val2)
        self.assertFalse(val3)

    def test_regpoly_classify(self):
        pent = Simple_Polygon("pent",np.array([(550,450), (645,519), (609,631), (491,631), (455,519)], dtype=np.float))

        output = { 0: [[True, np.array([[550., 450.], [645., 519.]]), 0.20004096545497796]],
                   1: [[True, np.array([[645., 519.], [609., 631.]]), 0.19989950974821585]],
                   2: [[True, np.array([[609., 631.], [491., 631.]]), 0.19989950974821585]],
                   3: [[True, np.array([[491., 631.], [455., 519.]]), 0.20004096545497796]],
                   4: [[True, np.array([[455., 519.], [550., 450.]]), 0.1999777758081791]]}
        pprime, segments = classifyBoundary(pent, 0.2)
        compare(list(segments), list(output))


    def test_edge_viz_graph_nonconv(self):
        init_poly = Simple_Polygon("sh", simple_nonconv[0])
        pls = Partial_Local_Sequence(init_poly)
        bvd = Bounce_Visibility_Diagram(pls)
        bvg = Bounce_Graph(bvd)
        result = list(bvg.visibility_graph.edges)
        expected = [(0, 3), (0, 4), (0, 5), (0, 6), (1, 2), (1, 3), (1, 4), (1,
        5), (2, 1), (2, 4), (2, 5), (2, 6), (3, 0), (3, 1), (3, 4), (3, 5), (3,
        6), (4, 0), (4, 1), (4, 2), (4, 3), (4, 5), (4, 6), (5, 0), (5, 1), (5,
        2), (5, 3), (5, 4), (6, 0), (6, 2), (6, 3), (6, 4)]

        self.assertEqual(sorted(result), sorted(expected))

    def test_edge_viz_tricky(self):
        init_poly = Simple_Polygon("sh", pinched_square[0])
        pls = Partial_Local_Sequence(init_poly)
        bvd = Bounce_Visibility_Diagram(pls)
        bvg = Bounce_Graph(bvd)
        result = list(bvg.visibility_graph.edges)
        expected = [(0,5), (0,6), (0,7), (0,11),
                    (1,5), (1,6), (1,7), (1,8), (1,11),
                    (2,5), (2,6), (2,7), (2,8), (2,9), (2,10), (2,11),
                    (3,5), (3,8), (3,9), (3,10), (3,11),
                    (4,5), (4,9), (4,10), (4,11),
                    (5,0), (5,1), (5,2), (5,3), (5,4), (5,10), (5,11),
                    (6,0), (6,1), (6,2), (6,8), (6,9), (6,10), (6,11),
                    (7,0), (7,1), (7,2), (7,8), (7,9), (7,10),
                    (8,1), (8,2), (8,3), (8,6), (8,7), (8,9), (8,10),
                    (9,2), (9,3), (9,4), (9,6), (9,7), (9,8),
                    (10,2), (10,3), (10,4), (10,5), (10,6), (10,7), (10,8),
                    (11,0), (11,1), (11,2), (11,3), (11,4), (11,5), (11,6)]

        missing = set(expected) - set(result)
        print("our algorithm should have seen but didn't:", missing)
        extra = set(result) - set(expected)
        print("our algorithm saw when it shouldn't:", extra)

        self.assertEqual(sorted(result), sorted(expected))

    def test_valid_transit(self):
        init_poly = Simple_Polygon("sh", pinched_square[0])
        pls = Partial_Local_Sequence(init_poly)
        test1 = check_valid_transit(10, 1, pls.inserted_polygon)
        self.assertFalse(test1)

    def test_mutual_vis(self):
        init_poly = Simple_Polygon("sh", pinched_square[0])
        pls = Partial_Local_Sequence(init_poly)
        vvs = get_all_edge_visible_vertices(pls.inserted_polygon)
        viz_from_10 = vvs[0][10][0]
        expected = [2,3,4,5,6,7,8,9]
        self.assertEqual(sorted(viz_from_10), sorted(expected))

    def test_edges_in_interval(self):
        poly = Simple_Polygon("sp",pinched_square[0])
        pls = Partial_Local_Sequence(poly)
        mapping = pls.inserted_polygon.unit_interval_mapping
        r1 = get_edges_intersecting_interval(mapping, (0.51, 0.66))
        r2 = get_edges_intersecting_interval(mapping, (0.05, 0.15))

        self.assertEqual(r1, [8])
        self.assertEqual(sorted(r2), [1,2,3])

    def test_SafeAngles(self):
        e1 = np.array([[0.0,0.0], [1.0, 0.0]])
        e2 = np.array([[1.0, 0.0], [1., 1.]])
        correct_thetas = (np.pi/4, 0.)
        compare(correct_thetas, SafeAngles(e1,e2))
        e3 = np.array([[1.,1.],[-1.,1.]])
        e4 = np.array([[-2.,-2.],[2.,-2.]])
        correct_thetas2 = (3.*np.pi/4., np.pi/4.)
        compare(correct_thetas2, SafeAngles(e3,e4))


#    def test_graph_reduce(self):
#        pls = Partial_Local_Sequence(square)
#        bvd = Bounce_Visibility_Diagram(pls)
#        G = Bounce_Graph(bvd)
#        G = mkGraph(square)
#        H = reduceGraphWrtAngle(G, 0.0, 0.2)
#        self.assertEqual(H.edge,
#        {0: {1: {'weight': [(0, 0.2)]}}, 1: {2: {'weight': [(0, 0.2)]}}, 2: {3: {'weight': [(0, 0.2)]}}, 3: {0: {'weight': [(0, 0.2)]}}}
#        )


#    def test_shootInterval(self):
#        # land inside interval
#        pt1, pt2 = ShootInterval(square, 3, 0, (-10,-250), (10,-250), pi/2, pi/2)
#        self.assertAlmostEqual(pt1[0], 10)
#        self.assertAlmostEqual(pt1[1], 250)
#        self.assertAlmostEqual(pt2[0], -10)
#        self.assertAlmostEqual(pt2[1], 250)
#        # land outside interval, clip to endpoints
#        pt3, pt4 = ShootInterval(square, 3, 0, (-10,-250), (10,-250), 0.1, pi-0.1)
#        self.assertAlmostEqual(pt3[0], 250)
#        self.assertAlmostEqual(pt3[1], 250)
#        self.assertAlmostEqual(pt4[0], -250)
#        self.assertAlmostEqual(pt4[1], 250)
#
#    def test_nodeParam(self):
#        self.assertEqual(nodesCovered(square,(0.1,0.6)),[0,1,2])
#        self.assertEqual(nodesCovered(square,(0.9,0.1)),[3,0])
#        self.assertEqual(nodesCovered(square,(0.1,0.2)),[0])

# if __name__ == '__main__':
#     unittest.main()
