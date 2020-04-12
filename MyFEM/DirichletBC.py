import numpy as np

# Dirichlet boundary condition class
class DirichletBC(object):
    def __init__(self, V, values, marked_facets, label, sub_spaces=[0]):
        self.V = V
        self.values = values
        self.label  = label
        self.marked_facets = marked_facets
        self.dofs  = self.find_dofs()
        self.sub_spaces = sub_spaces

    def find_dofs(self):
        """ Apply the boundary condition 'value' to the function 'u' on lines
        labeled as 'label'. The physical information about the boundaries is 
        contained in mesh 
        """
        # Facets connectivity
        facets_connectivity = self.V.mesh().facets_connectivity()
        
        # Find dofs with label == marker
        dofs = facets_connectivity[self.marked_facets[self.label]].flatten('F')

        return dofs.astype(np.int32)
        
    # Apply dirichlet bc and return sparse matrices
    def apply(self, A, b):
        """ Apply boundary conditions to the matrix A and vector b
        """
        # Number of dofs per node
        n = self.V.element_shape()[0]
        
        # Apply boundary conditions
        for i in self.sub_spaces:
            # Subspace dofs and value
            sub_dofs = self.dofs[i::n]

            # LHS
            if isinstance(self.values, list):
                if isinstance(self.values[i], float):
                    b.setValues(sub_dofs, self.values[i]*np.ones(sub_dofs.shape))
                elif self.values[i].size > 1: 
                    b.setValues(sub_dofs, self.values[i][sub_dofs])
            else:
                if isinstance(self.values, float):
                    b.setValues(sub_dofs, self.values*np.ones(sub_dofs.shape), addv=False)
                elif self.values.size > 1:
                    b.setValues(sub_dofs, self.values[sub_dofs], addv=False)

            # LHS
            A.setOption(18, True)
            A.zeroRows(sub_dofs, diag=1)

            # Finalize assembly
            A.assemblyBegin()
            A.assemblyEnd()

        return True
