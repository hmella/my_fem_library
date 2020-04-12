import numpy as np
from MatrixOperations import *

#################################
#   Scalar space variables
#################################
def ShapeTri(p, shape):
    """ Returns a list of shape functions evaluated in the point
    (r, s) in the reference cell.
    """    
    # Point components
    r = p[0]
    s = p[1]

    # Vector dimension
    l = shape[0]

    # Basis functions
    phi = np.array([1-r-s, r, s])

    # Interpolation functions
    N = np.zeros((l, 3*l))
    for i in range(l):
        N[i, i::l] = phi
    
    return N


def ShapeQuad(p, shape):
    """ Returns a list of shape functions evaluated in the point
    (r, s) in the reference cell.
    """    
    # Point components
    r = p[0]
    s = p[1]

    # Vector dimension
    l = shape[0]

    # Basis functions
    phi = 0.25*np.array([(1 - r)*(1 - s),
                 (1 + r)*(1 - s),
                 (1 + r)*(1 + s),
                 (1 - r)*(1 + s)])

    # Interpolation functions
    N = np.zeros((l, 4*l))
    for i in range(l):
        N[i, i::l] = phi
    
    return N


def ShapeHex(p, shape):
    """ Returns a list of shape functions evaluated in the point
    (r, s, t) in the reference cell.
    """
    # Point components
    r = p[0]
    s = p[1]
    t = p[2]

    # Vector dimension
    l = shape[0]
    
    # Basis functions
    phi = 0.125*np.array([(1 - r)*(1 - s)*(1 - t),
                          (1 + r)*(1 - s)*(1 - t),
                          (1 + r)*(1 + s)*(1 - t),
                          (1 - r)*(1 + s)*(1 - t),
                          (1 - r)*(1 - s)*(1 + t),
                          (1 + r)*(1 - s)*(1 + t),
                          (1 + r)*(1 + s)*(1 + t),
                          (1 - r)*(1 + s)*(1 + t)])

    # Interpolation functions
    N = np.zeros((l, 8*l))
    for i in range(l):
        N[i, i::l] = phi
    
    return N

# Derivatives interpolator
def DerivativeInterpTri(p, x):
    """ TODO: Add description here!
    """
    # Point components
    r = p[0]
    s = p[1]

    # Derivative interpolator in quadrature points
    dhdx = np.array([
            [-1.0, 1.0, 0.0],
            [-1.0, 0.0, 1.0]])

    # determinant and jacobian matrix
    det, inv_J = JacobianTri(dhdx, x[:,0:2])        

    # Derivative interpolator on the global cell
    dhdx = inv_J.dot(dhdx)

    return det, dhdx


def DerivativeInterpQuad(p, x):
    """ TODO: Add description here!
    """
    # Point components
    r = p[0]
    s = p[1]

    # Derivative interpolator in quadrature points
    dhdx = 0.25*np.array([
            [s - 1, -s + 1, s + 1, -s - 1],
            [r - 1, -r - 1, r + 1, -r + 1]])

    # determinant and jacobian matrix
    det, inv_J = JacobianQuad(dhdx, x[:,0:2])        

    # Derivative interpolator on the global cell
    dhdx = inv_J.dot(dhdx)

    return det, dhdx


def DerivativeInterpHex(p, x):
    """ TODO: Add description here!
    """
    # Point components
    r = p[0]
    s = p[1]
    t = p[2]

    # Derivative interpolator in quadrature points
    dhdx = 0.125*np.array([[-(1 - s)*(1 - t),
                             (1 - s)*(1 - t), 
                             (1 + s)*(1 - t), 
                            -(1 + s)*(1 - t), 
                            -(1 - s)*(1 + t),
                             (1 - s)*(1 - t),
                             (1 + s)*(1 + t),
                            -(1 + s)*(1 + t)],
                           [-(1 - r)*(1 - t),
                            -(1 + r)*(1 - t),
                             (1 + r)*(1 - t),
                             (1 - r)*(1 - t),
                            -(1 - r)*(1 + t),
                            -(1 + r)*(1 + t),
                             (1 + r)*(1 + t),
                             (1 - r)*(1 + t)],
                           [-(1 - r)*(1 - s),
                            -(1 + r)*(1 - s),
                            -(1 + r)*(1 + s),
                            -(1 - r)*(1 + s),
                             (1 - r)*(1 - s),
                             (1 + r)*(1 - s),
                             (1 + r)*(1 + s),
                             (1 - r)*(1 + s)]])

    # determinant and jacobian matrix
    det, inv_J = JacobianQuad(dhdx, x)        

    # Derivative interpolator on the global cell
    dhdx = inv_J.dot(dhdx)

    return det, dhdx    


# Strain displacement interpolator
def stdHex(p, x):
    # Point components
    r = p[0]
    s = p[1]
    t = p[2]

    # Derivative interpolator in quadrature points
    dhdx = 0.125*np.array([[-(1 - s)*(1 - t),
                             (1 - s)*(1 - t), 
                             (1 + s)*(1 - t), 
                            -(1 + s)*(1 - t), 
                            -(1 - s)*(1 + t),
                             (1 - s)*(1 - t),
                             (1 + s)*(1 + t),
                            -(1 + s)*(1 + t)],
                           [-(1 - r)*(1 - t),
                            -(1 + r)*(1 - t),
                             (1 + r)*(1 - t),
                             (1 - r)*(1 - t),
                            -(1 - r)*(1 + t),
                            -(1 + r)*(1 + t),
                             (1 + r)*(1 + t),
                             (1 - r)*(1 + t)],
                           [-(1 - r)*(1 - s),
                            -(1 + r)*(1 - s),
                            -(1 + r)*(1 + s),
                            -(1 - r)*(1 + s),
                             (1 - r)*(1 - s),
                             (1 + r)*(1 - s),
                             (1 + r)*(1 + s),
                             (1 - r)*(1 + s)]])

    # Interpolator
    B = np.zeros((6, 24))
    B[0, 0::3] = dhdx[:,0]
    B[1, 1::3] = dhdx[:,1]
    B[2, 2::3] = dhdx[:,2]
    B[3, 0::3] = dhdx[:,1]
    B[3, 1::3] = dhdx[:,0]
    B[4, 1::3] = dhdx[:,2]
    B[4, 2::3] = dhdx[:,1]
    B[5, 0::3] = dhdx[:,2]
    B[5, 2::3] = dhdx[:,0]

    return B

#################################
#   Jacobians and derivatives
#################################    

# Jacobian matrix
def JacobianTri(dNdx, x):
    """ Jacobian of the shape functions on the reference cell
    """
    # Jacobian matrix
    J = np.matmul(dNdx, x)

    # Inverse Jacobian matrix
    det_J = det2(J)
    inv_J = inv2(J)

    return det_J, inv_J
    
def JacobianQuad(dNdx, x):
    """ Jacobian of the shape functions on the reference cell
    """
    # Jacobian matrix
    J = np.matmul(dNdx, x)

    # Inverse Jacobian matrix
    det_J = det2(J)
    inv_J = inv2(J)

    return det_J, inv_J
    
def JacobianHex(dNdx, x):
    """ Jacobian of the shape functions on the reference cell
    """
    # Jacobian matrix
    J = np.matmul(dNdx, x)

    # Inverse Jacobian matrix
    det_J = np.det3(J)
    inv_J = np.inv3(J)

    return det_J, inv_J
