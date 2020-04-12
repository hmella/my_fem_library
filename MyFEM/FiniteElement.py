from MyFEM.Geometry import *
from MyFEM.ShapeFunctions import *
from MyFEM.Quadrature import *


# Base class for finite elements
class FiniteElementBase(object):
    def __init__(self, kind, degree=1, quadrature_degree=3):
        self.kind = kind
        self.quadrature_degree = quadrature_degree
        self.quadrature_weights, self.quadrature_points = self.quadrature()
        self.degree = degree

    def quadrature(self):
        if self.kind=="triangle":
            x, w = triangle_scheme(self.quadrature_degree)
        if self.kind=="quadrilateral":
            x, w = quadrilateral_scheme(self.quadrature_degree)
        if self.kind=="hexahedron":
            x, w = hexahedron_scheme(self.quadrature_degree)
        if self.kind=="tetrahedron":
            x, w = tetrahedron_scheme(self.quadrature_degree)
        return w, x

    # Get element degree
    def get_degree(self):
        return self.degree

    # Geometric dimension
    def dimension(self):
        if cells_dict[self.kind] in CELL_TYPES_2:
            dim = 2
        elif cells_dict[self.kind] in CELL_TYPES_3:
            dim = 3
        return dim

    # Element shape function interpolator
    def get_shape_functions(self):
        if self.degree==1 and self.kind=="triangle":
            shape_function = ShapeTri
        if self.degree==1 and self.kind=="quadrilateral":
            shape_function = ShapeQuad
        if self.degree==1 and self.kind=="hexahedron":
            shape_function = ShapeHex

        return shape_function

    # Element derivative function interpolator
    def get_derivative_interpolator(self):
        if self.degree==1 and self.kind=="triangle":
            derivative_interpolator = DerivativeInterpTri
        if self.degree==1 and self.kind=="quadrilateral":
            derivative_interpolator = DerivativeInterpQuad
        if self.degree==1 and self.kind=="hexahedron":
            derivative_interpolator = DerivativeInterpHex

        return derivative_interpolator

    # Number of dofs per cell
    def get_dofs_per_cell(self):
        return _get_dofs_per_cell(self.kind, self.shape())


# Finite element class
class FiniteElement(FiniteElementBase):
    def __init__(self, *args, **kwargs):
        super(FiniteElement, self).__init__(*args, **kwargs)

    # Shape of finite element
    def shape(self):
        return (1,)


# Vector element class
class VectorElement(FiniteElementBase):
    def __init__(self, *args, **kwargs):
        super(VectorElement, self).__init__(*args, **kwargs)

    def shape(self):
        return (self.dimension(),)


# Number of dofs per cell
def _get_dofs_per_cell(cell_type, shape):
  if cell_type == "triangle":
    n = TRIANGLE
  elif cell_type == "quadrilateral":
    n = QUADRILATERAL
  elif cell_type == "tetrahedron":
    n = TETRAHEDRON
  elif cell_type == "tetrahedron10":
    n = TETRAHEDRON10
  elif cell_type == "hexahedron":
    n = HEXAHEDRON
  return shape[0]*n
