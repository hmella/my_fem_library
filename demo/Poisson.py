import numpy as np
from MyFEM import *
import time

if __name__ == "__main__":

    # Mesh import
    mesh = Mesh("mesh/meshq.msh")

    # Cell and finite element
    cell = mesh.cell_type()
    FE = FiniteElement(cell, quadrature_degree=1)
#    FE = VectorElementElement(cell, quadrature_degree=3)   # TODO: This doesn't work!
#                                                          the problem is in the derivatives interpolator 
#                                                          (not defined for vector-valued finite elements)

    # Function space
    V = FunctionSpace(mesh, FE)

    # Assemble mass and stiffness matrices
    start = time.time()
    M, K = assemble_(V)
    end1 = time.time()-start

    start = time.time()
    M_, K_ = assemble(V)
    end2 = time.time()-start
    print("--- Speed-up: {:.2f}---".format(end2/end1))

    # Create RHS
    u, b = K.getVecs()
    u.set(0)
    b.set(0.001)
    
    # Dirichlet boundary condition
    dofmap = V.vertex_to_dof_map()
    x = V.dof_coordinates()

    # Boundary condition
    Omega = Domain(mesh, build_cells=False, build_facets=True)
    bc1 = DirichletBC(V, np.sin(10*x[dofmap,0]), Omega.facets(), 1)
    bc2 = DirichletBC(V, 0.0, Omega.facets(), 3)
    bcs = [bc1, bc2]

    # Apply Dirichlet boundary conditions
    [bc.apply(K, b) for bc in bcs]
        
    # Solve system
    # create linear solver
    ksp = PETSc.KSP()
    ksp.create(PETSc.COMM_WORLD)
    # # use conjugate gradients
    # ksp.setType('cg')
    # # and incomplete Cholesky
    # ksp.getPC().setType('icc')
    # and next solve
    ksp.setOperators(K)
    ksp.setFromOptions()
    ksp.solve(b, u)

    # Export solution
    U = Function(V)
    U.vector()[:] = u[...]
    write_vtk(U, path="output/u.vtk", name="u")