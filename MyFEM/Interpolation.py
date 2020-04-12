# Geomtric interpolator
def GeometricInterpolator(xi, eta, x):
    """ Allows to estimate the transformation of the point (xi, eta)
    in the reference cell to a global coordinate
    """
    # Shape functions interpolator
    N = ShapeTri(xi, eta)
    
    # Coordinates interpolator
    y = np.zeros(2)
    y[0] = N.dot(x[:,0])
    y[1] = N.dot(x[:,1])

    return y


# Baricentric coordinates
def BaricentricCoordinates(x, X):
    """ Returns the weigths lambda to intepolate a function
    at the point x, which lies into the triangulation of X
    """
    # Baricentric coordinates
    lmbda = np.zeros(3)
    lmbda[0] = ((X[1,1] - X[2,1])*(x[0] - X[2,0]) + (X[2,0] - X[1,0])*(x[1] - X[2,1])) \
             / ((X[1,1] - X[2,1])*(X[0,0] - X[2,0]) + (X[2,0] - X[1,0])*(X[0,1] - X[2,1]))
    lmbda[1] = ((X[2,1] - X[0,1])*(x[0] - X[2,0]) + (X[0,0] - X[2,0])*(x[1] - X[2,1])) \
             / ((X[1,1] - X[2,1])*(X[0,0] - X[2,0]) + (X[2,0] - X[1,0])*(X[0,1] - X[2,1]))
    lmbda[2] = 1 - lmbda[0] - lmbda[1]

    return lmbda

# Interpolate function
def interpolate(u, dof_indices, x, X):
    """ Interpolate the value of the Function u at the point x,
    which is inside of the triangulation given by the points 
    cell_coordinates
    """
    # Baricentric coordinates
    lmbda = BaricentricCoordinates(x, X)

    # Interpolation
    ui = lmbda.dot(u.vector()[dof_indices])

    return ui
