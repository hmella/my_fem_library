from MyFEM.Geometry import *
from MyFEM.Quadrature import *
from MyFEM.ShapeFunctions import *
import scipy.sparse as sp
from AssembleUtilities import *
import petsc4py.PETSc as PETSc

def assemble_(V):
    print("Assembling matrix for linear system with {:.0f} dofs".format(V.num_dofs()))

    # Finite element
    element = V.get_element()

    # Retrieve useful information
    CellDofs        = V.cells_dofmap()             # cell dofs
    DofsCoordinates = V.dof_coordinates()          # dof coordinates
    DofsInCell      = element.get_dofs_per_cell()  # number of dofs in each cell

    # Sizes
    ndofs  = V.num_dofs()          # number of global dofs
    ncells = V.mesh().num_cells()  # number of global cells

    # Quadrature rule
    w = element.quadrature_weights
    x = element.quadrature_points

    # Derivatives interpolator
    derivative_interpolator = V.DerivativeInterp

    # Shape functions
    shape_functions = V.shape_function

    # Shape of the function space
    function_space_shape = V.shape

    # Assemble
    M, K = Pybind11Assemble(V, CellDofs, DofsCoordinates,
                            DofsInCell, ndofs, shape_functions,
                            derivative_interpolator, w, x)

    return M, K



# Local cell assemble
def assemble_cell(qweights, 
                  qpoints,
                  dofs_per_cell,
                  derivative_interpolator,
                  shape_functions,
                  function_space_shape,
                  cell_x):
    """ Assemble local matrices of mass and stiffnes:
        Input
            x: global coordinates of the element
        Output
            me: local mass matrix
            ke: local stiffnes matrix
    """
    # Number of quadrature points
    n = qweights.size

    # Local stiffnes and mass matrices
    ke = np.zeros([dofs_per_cell, dofs_per_cell])
    me = np.zeros([dofs_per_cell, dofs_per_cell])

    # Assemble local
    for i in range(n):
        # Assemble
        det, B = derivative_interpolator(qpoints[i,:], cell_x)
        N  = shape_functions(qpoints[i,:], function_space_shape)
        ke += np.dot(B.T, B)*det*qweights[i]
        me += np.dot(N.T, N)*det*qweights[i]

    return ke, me

# TODO: Local facet assemble
def assemble_facet():
  return True

# Global assembler
# TODO: add list of cells for cell assembling
# TODO: add list of facets for facet assembling (Neumann conditions)
def assemble(V):
    print("Assembling matrix for linear system with {:.0f} dofs".format(V.num_dofs()))

    # Finite element
    element = V.get_element()

    # Retrieve useful information
    cell_dofs    = V.cells_dofmap()
    dof_coords   = V.dof_coordinates()
    dofs_per_cell = element.get_dofs_per_cell()

    # Global stiffnes matrices for shape functions and derivatives
    [m, _] = dof_coords.shape
    [n, N] = cell_dofs.shape
    N = N*N

    rows = np.zeros([dofs_per_cell**2*n,])   # row indices
    cols = np.zeros([dofs_per_cell**2*n,])   # col indices
    K = np.zeros([dofs_per_cell**2*n,])      # stiffnes matrix in vector format
    M = np.zeros([dofs_per_cell**2*n,])      # mass matrix in vector format

    # Quadrature rule
    quadrature_weights = element.quadrature_weights
    quadrature_points  = element.quadrature_points

    # Derivatives interpolator
    derivative_interpolator = V.DerivativeInterp

    # Shape functions
    shape_functions = V.shape_function

    # Shape of the function space
    function_space_shape = V.shape

    # Assembling loop
    for I in range(n):
        # Dofs indices in cell
        c = cell_dofs[I,:]

        # Extract global coordinates of cell nodes
        x = dof_coords[c]

        # Local matrices and row-col indices
        ke, me = assemble_cell(quadrature_weights, 
                  quadrature_points,
                  dofs_per_cell,
                  derivative_interpolator,
                  shape_functions,
                  function_space_shape,
                  x)
    
        RC = np.array([[i,j] for i in c for j in c])
           
        # Fill global matrices
        M[N*I:N*(I+1)] = me.flatten()
        K[N*I:N*(I+1)] = ke.flatten()
        rows[N*I:N*(I+1)] = RC[:,0]
        cols[N*I:N*(I+1)] = RC[:,1]

    # Sparse representation of matrices
    M = sp.coo_matrix((M, (rows, cols)), shape=(m, m)).tolil()
    K = sp.coo_matrix((K, (rows, cols)), shape=(m, m)).tolil()

    return M, K