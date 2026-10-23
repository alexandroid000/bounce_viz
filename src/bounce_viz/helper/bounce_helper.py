
def PolyToWindowScale(poly, ydim):
    newpoly = []
    for p in poly:
        newpoly.append((p[0]/2.0, ydim - p[1]/2.0))
    return newpoly

def PointToWindowScale(p, ydim):
     return (p[0]/2.0, ydim - p[1]/2.0)

def AngleBetween(theta,theta1,theta2):
    return (((theta1 <= theta) and (theta <= theta2)) or \
            ((theta2 <= theta1) and (theta1 <= theta)) or \
            ((theta <= theta2) and (theta2 <= theta1)))

def AngleDistance(theta1,theta2):
    d = abs(theta1-theta2)
    return min(d,2.0*pi-d)

def AngleDifference(theta1,theta2):
    d = theta2 - theta1
    if (d < pi):
        d += 2.0*pi;
    if (d > pi):
        d -= 2.0*pi;
    return d
