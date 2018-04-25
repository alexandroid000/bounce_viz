#! /usr/bin/env python

from geom_utils import *
from maps import *
import numpy as np

# Used to plot the example
import matplotlib.pyplot as plt

def VizRay():
    poly = poly1
    # Set the title
    wall_x = [x for (x,y) in poly]
    wall_x.append(wall_x[0])
    wall_y = [y for (x,y) in poly]
    wall_y.append(wall_y[0])
    
    # Plot the outer boundary with black color
    plt.figure()
    plt.plot(wall_x, wall_y, 'black')

    pt1 = (0.0, 0.0)
    pt2 = (150.0, 50.0)

    plt.plot([pt1[0]], [pt1[1]], 'go')
    plt.plot([pt2[0]], [pt2[1]], 'go')

    # mark reflex vertices
    rvs = FindReflexVerts(poly)
    for p in rvs:
        plt.plot([poly[p][0]], [poly[p][1]], 'bo')
        r1, r2 = ShootRaysFromReflex(poly, p)
        transition_pts = ShootRaysToReflexFromVerts(poly,p)
        print("reflex ",p," transitions ",transition_pts)
        transition_pts.extend([r1,r2])
        for (pt,k) in transition_pts:
            plt.plot([poly[p][0], pt[0]], [poly[p][1], pt[1]], 'green')
    
    print(rvs)

    plt.savefig('test.pdf')

def bounce_viz_diagram():
    poly = poly1
    psize = len(poly)
    plt.figure()
    bounce_viz_heat_map = []
    for i in range(psize):
        p1, p2 = poly[i], poly[(i+1) % psize]
        v1 = (p2[0]-p1[0],p2[1]-p1[1])
        visibleVertices = GetVisibleVertices(poly, i)
        sweepTheta = [GetVector2Angle(v1, (poly[j][0]-p1[0],poly[j][1]-p1[1])) for j in visibleVertices[1:]]
        for j in range(len(sweepTheta)):
            bounce_viz_heat_map.append((i, sweepTheta[j], j))
    print(bounce_viz_heat_map)
    plt.plot(bounce_viz_heat_map)
    plt.savefig('visible.pdf')


if __name__ == '__main__':
    VizRay()
    # bounce_viz_diagram()
