import numpy as np
import matplotlib.pyplot as plt

def FindLE_top(X):
    """Return index dividing upper and lower surface given MSES geometry.
    Search along upper surface until LE.
    MSES files start at rear of airfoil, and x diminishes until the leading
    edge, where it then increases back to the trailing edge.  This code finds
    the transition where x goes from decreasing to increasing.
    X --> MSES x coordinates
    """
    xold = X[0]
    for i, x in enumerate(X[1:]):
        if x >= xold:
            #If current x greater/equal to prev x, x is increasing (lower surf)
            return i #return index of Leading Edge (divides upper/lower surfs)
        else:
            #If current x less than prev x, x still diminishing (upper surf)
            xold = x

def FindLE_bot(X):
    """Return index dividing upper and lower surface given MSES geometry.
    Search along lower surface until LE.
    MSES files start at rear of airfoil, and x diminishes until the leading
    edge, where it then increases back to the trailing edge.  This code finds
    the transition where x goes from decreasing to increasing.
    X --> MSES x coordinates
    """
    Xreverse = X[::-1]
    xold = Xreverse[0]
    for i, x in enumerate(Xreverse[1:]):
        if x >= xold:
            #If current x greater/equal to prev x, x is increasing (on upper surf)
            return len(X) - 1 - i #return index of Leading Edge (divides upper/lower surfs)
        else:
            #If current x less than prev x, x still diminishing (still on lower surf)
            xold = x

def MsesSplit(x, y):
    """Split MSES format into upper and lower surfaces.
    Find LE from MSES x geometry coordinates,
    Split y at this index(s).
    If LE point is at y=0, include in both sets of data.
    Return y split into upper/lower surfaces, with LE overlapping
    x --> MSES x coordinates
    y --> Any other MSES parameter (e.g. x/c, z/c, Cp, etc)
    """
    #FIND LE FROM BOTH SIDES (DETECT SHARED LE POINT)
    #Get index of leading edge starting from upper surface TE
    iLE_top = FindLE_top(x)
    #Get index of leading edge starting from lower surface TE
    iLE_bot = FindLE_bot(x)
    #Split upper and lower surface, reverse order upper surface
    up = y[iLE_top::-1]
    lo = y[iLE_bot:]
    return up, lo

def MsesInterp(xout, xmses, ymses):
    """Split MSES format data into upper and lower surfaces.  Then
    interpolate data to match given xout vector.
    xout  --> desired x locations
    xmses --> original x MSES data
    ymses --> original x/c, z/c, Cp, etc MSES data
    """
    xup_mses, xlo_mses = MsesSplit(xmses, xmses)
    yup_mses, ylo_mses = MsesSplit(xmses, ymses)
    yup = np.interp(xout, xup_mses, yup_mses)
    ylo = np.interp(xout, xlo_mses, ylo_mses)
    return yup, ylo


x,y = np.loadtxt('clarky10.dat', unpack=True,skiprows=1)

xnew = np.linspace(0,1,21)
upper, lower = MsesInterp(xnew, x, y)

camber  = (upper+lower)/2


plt.figure()
plt.plot(xnew,upper)
plt.plot(xnew,lower)
plt.plot(xnew,camber)
plt.axis('equal')
plt.show()

file = open("camberclarky1.0txt", "w")
file.write("{}\n\n{}".format(xnew,camber))
file.close()