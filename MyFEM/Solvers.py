import numpy as np
from scipy import sparse as sp
from scipy.sparse import linalg as spla


# Direct solver for linear system
def solve(A, u, b):
    print("Solving system using direct inversion")
    """ Solves the linear system A*u = b
        Input
            A: csr_matrix of NxN components (LHS)
            b: csr_matrix of Nx1 components (RHS)
            u: numpy array of Nx1 components (Solution)
    """
    
    # Change RHS representation
    A = sp.csr_matrix(A)
    b = sp.csr_matrix(b).T
    
    # Solve system
    u[:] = spla.spsolve(A, b)

    
# GMRES solver
def gmres_solve(A, u, b, tol=1e-08):
    print("Solving system using GMRES solver")
    """ Solves the linear system A*u = b
        Input
            A: numpy array of NxN components (LHS)
            b: numpy array of Nx1 components (RHS)
            u: numpy array of Nx1 components (Solution)
    """

    # Change RHS representation
    A = sp.csr_matrix(A)
        
    # Solve system
    u[:] = spla.gmres(A, b, tol=tol)[0]
    

# CG solver   
def cg_solve(A, u, b, tol=1e-08):
    print("Solving system using CG solver")
    """ Solves the linear system A*u = b
        Input
            A: numpy array of NxN components (LHS)
            b: numpy array of Nx1 components (RHS)
            u: numpy array of Nx1 components (Solution)
    """
    
    # Change RHS representation
    A = sp.csr_matrix(A)
        
    # Solve system
    u[:] = spla.cg(A, b, tol=tol)[0]


# CGS solver   
def cgs_solve(A, u, b, tol=1e-08):
    print("Solving system using CGS solver")
    """ Solves the linear system A*u = b
        Input
            A: numpy array of NxN components (LHS)
            b: numpy array of Nx1 components (RHS)
            u: numpy array of Nx1 components (Solution)
    """
    
    # Change RHS representation
    A = sp.csr_matrix(A)
        
    # Solve system
    u[:] = spla.cgs(A, b, tol=tol)[0]
    

# MINRES solver   
def minres_solve(A, u, b, tol=1e-08):
    print("Solving system using MINRES solver")
    """ Solves the linear system A*u = b
        Input
            A: numpy array of NxN components (LHS)
            b: numpy array of Nx1 components (RHS)
            u: numpy array of Nx1 components (Solution)
    """
    
    # Change RHS representation
    A = sp.csr_matrix(A)
        
    # Solve system
    u[:] = spla.minres(A, b, tol=tol)[0]
    
    
# BICG solver   
def bicg_solve(A, u, b, tol=1e-08):
    print("Solving system using BICG solver")
    """ Solves the linear system A*u = b
        Input
            A: numpy array of NxN components (LHS)
            b: numpy array of Nx1 components (RHS)
            u: numpy array of Nx1 components (Solution)
    """
    
    # Change RHS representation
    A = sp.csr_matrix(A)
        
    # Solve system
    u[:] = spla.bicg(A, b, tol=tol)[0]
    
    
# BICGSTAB solver   
def bicgstab_solve(A, u, b, tol=1e-08):
    print("Solving system using BICGSTAB solver")
    """ Solves the linear system A*u = b
        Input
            A: numpy array of NxN components (LHS)
            b: numpy array of Nx1 components (RHS)
            u: numpy array of Nx1 components (Solution)
    """
    
    # Change RHS representation
    A = sp.csr_matrix(A)
        
    # Solve system
    u[:] = spla.bicgstab(A, b, tol=tol)[0]
    
    
# QSR solver   
def qmr_solve(A, u, b, tol=1e-08):
    print("Solving system using QMR solver")
    """ Solves the linear system A*u = b
        Input
            A: numpy array of NxN components (LHS)
            b: numpy array of Nx1 components (RHS)
            u: numpy array of Nx1 components (Solution)
    """

    # Change RHS representation
    A = sp.csr_matrix(A)
        
    # Solve system
    u[:] = spla.qmr(A, b, tol=tol)[0]
