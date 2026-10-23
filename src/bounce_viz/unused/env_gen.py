from math import cos, sin, pi
from random import random, uniform

def reg_poly(N, r, x, y):
    ''' generate regular polygon with n sides, radius r, centered at x,y

    Parameters
    ----------
    N : int
        The size of the polygon
    r : float
        The radius of the polygon
    (x,y) : (float, float)
        The coordinates of the center of the polygon
    
    Returns
    -------
    list of (float, float)
        The regular polygon represented by vertex coordinates
    '''
    return [(x + r * cos(2*pi*n/N), y + r * sin(2*pi*n/N)) for n in range(N)]


def room(r, x0, y0):
    top = 0.5*random()
    bottom = uniform(0.5, 1)

    base =  [(1.1*r, 1.1*r), (1.1*r, r*top), (r, r*top),
    (r,r),(-r,r),(-r,-r),(r,-r),
    (r,-r*bottom), (1.1*r, -r*bottom), (1.1*r, -1.1*r)]
    return [(x+x0, y+y0) for (x,y) in base]

def rooms(N, r, x, y):
    # room with r=a will have height h=2.2a
    env = []
    for i in range(N):
        env.extend(room(r, 0, -2.3*r*i))

    env.extend([(2*r,-2.1*r*N), (2*r,1.1*r)])
    return env

