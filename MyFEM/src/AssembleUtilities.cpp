#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <petsc/petscmat.h>
#include "petsc_caster.h"
#include <tuple>

namespace py = pybind11;
using RowMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

// PETSc Matrix
Mat PETScMatrix(const PetscInt rows, const PetscInt cols, const PetscInt nz){
  // Create matrix
  Mat M;
  MatCreateSeqAIJ(PETSC_COMM_WORLD, rows, cols, nz, NULL, &M);

  // Set options
  // MatSetOption(M, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);

  return M;
}

// TODO: PETSc Vector

void LocalAssemble(const py::object &Interpolator,
                   const py::object &ShapeFunctions,
                   const py::object &shape,
                   const Eigen::MatrixXd &x,
                   Eigen::MatrixXd &ke,
                   Eigen::MatrixXd &me,
                   const Eigen::VectorXd &W,
                   const Eigen::MatrixXd &P,
                   const int n) {

    // Output
    std::tuple<py::object, py::object> out; 
    Eigen::MatrixXd B, N;
    double det;
    
    // Set local matrices to zero
    ke.setZero();
    me.setZero();

    // Assemble
    for (int i=0; i < n; i++) {
        // Determinant and interpolator
        out = Interpolator(P.row(i), x).cast<std::tuple<py::object, py::object>>();
        det = std::get<0>(out).cast<double>();
        B   = std::get<1>(out).cast<Eigen::MatrixXd>();

        // Shape functions at quadrature points
        N   = ShapeFunctions(P.row(i), shape).cast<Eigen::MatrixXd>();

        // Sum local matrices
        ke += (B.transpose()*B)*det*W(i);
        me += (N.transpose()*N)*det*W(i);
    }
}


// Global assemble
std::tuple<Mat, Mat> Pybind11Assemble(const py::object &V,  // function space
                                      const Eigen::MatrixXi &CellDofs, // dofs in each cell
                                      const Eigen::MatrixXd &DofsCoordinates, // dofs coordinates
                                      const std::size_t &d, // number of dofs in each cell
                                      const std::size_t &m, // number of global dofs
                                      const py::object &ShapeFunctions, // Shape function
                                      const py::object &Interpolator, // derivative interpolator
                                      const Eigen::VectorXd &qw, // quadrature weights
                                      const Eigen::MatrixXd &qx){ // quadrature points)

  // Mass and stiffness matrices
  Mat M = PETScMatrix(m, m, 10);
  Mat K = PETScMatrix(m, m, 10);

  // Local coordinates matrix
  Eigen::MatrixXd x(d, 3);

  // Local mass and stiffness matrices
  Eigen::MatrixXd ke(d,d);
  Eigen::MatrixXd me(d,d);

  // Function space shape
  const py::object shape = V.attr("shape");
  const int n = qw.size();   // number of quadrature points

  // PETSc object for maps
  int c[d];
  double Ke[d*d];
  double Me[d*d];

  // Iterators
  std::size_t i, j;

  // Assembling loop
  for (i=0; i < CellDofs.rows(); i++){

    // Cell indices and coordinates
    for (j=0; j<d; j++) {
      c[j] = CellDofs(i,j);
      x.row(j) = DofsCoordinates.row(c[j]);
    }

    // Assemble local matrices
    LocalAssemble(Interpolator, ShapeFunctions, shape, x, ke, me, qw, qx, n);

    // Map eigen matrix to PetscScalar arrays
    Eigen::Map<RowMatrixXd>(&Ke[0], ke.rows(), ke.cols()) = ke;
    Eigen::Map<RowMatrixXd>(&Me[0], me.rows(), me.cols()) = me;

    // Fill global matrices
    MatSetValues(K, d, c, d, c, Ke, ADD_VALUES);
    MatSetValues(M, d, c, d, c, Me, ADD_VALUES);
  }
  
  // Finalize assembling
  MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);

  return std::make_tuple(M, K);
}


PYBIND11_MODULE(AssembleUtilities, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring
    m.def("Pybind11Assemble", &Pybind11Assemble, "Assembling test");
}
