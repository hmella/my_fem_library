import numpy as np

# Function
class Function(object):
    def __init__(self, V):
        self.V = V
        self.array = np.zeros([V.num_dofs(), ])

    # Function space
    def function_space(self):
        return self.V
        
    # Assign method
    def assign(self, u):
        self.array[:] = u
        
    # Vector
    def vector(self):
        return self.array
        
    # Derivative
    def grad(self):
        # Retrieve mesh and some parameters
        mesh     = self.V.mesh()
        vertices = mesh.vertex_coordinates()
        cells    = mesh.cells_connectivity()
        num_cells = mesh.num_cells()
        
        # Create gradient array
        geo_dim = mesh.geometric_dimension()
        grad_u  = np.zeros([geo_dim*num_cells,])

        # Quadrature points and weights
        W, P = QuadratureTri(self.function_space().quadrature_order)

        # Estimate derivatives
        for c in range(num_cells):
            # Cell connectivity and dofs coordinates
            cell = cells[c]
            x = vertices[cell]
            
            # Derivative interpolator in quadrature points
            for i in range(W.size):            
                # Derivative interpolator in global coordinates
                det, B = DerivativeInterp(P[i,0], P[i,1], x[:,0:2])

                grad_u[2*c]   += B[0,:].dot(self.array[cell])
                grad_u[2*c+1] += B[1,:].dot(self.array[cell])
            
        return grad_u
