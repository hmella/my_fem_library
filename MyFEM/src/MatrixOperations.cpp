#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <Eigen/LU>

namespace py = pybind11;

// Determinant of 2x2 matrix
double Determinant2x2(const Eigen::Matrix<double,2,2> &A) {

  // Determinant
  const double det = A(0,0)*A(1,1) - A(0,1)*A(1,0);

  return det;
}

// Inverse of a 2x2 matrix
Eigen::Matrix<double,2,2> Inverse2x2(const Eigen::Matrix<double,2,2> &A) {

  // Determinant
  const double det = Determinant2x2(A);

  // Inversion
  Eigen::Matrix<double,2,2> A_inv(2,2);
  A_inv(0,0) = A(1,1);
  A_inv(0,1) = -A(0,1);
  A_inv(1,0) = -A(1,0);
  A_inv(1,1) = A(0,0);
  A_inv /= det;

  return A_inv;
}


// Determinant of a nxn matrix (n>2)
double Determinant(const Eigen::MatrixXd &A) {

  // Determinant
  double det = A.determinant();

  return det;
}

// Inverse of a nxn matrix (n>2)
Eigen::MatrixXd Inverse(const Eigen::MatrixXd &A) {

  // Inversion
  Eigen::MatrixXd A_inv = A.inverse();

  return A_inv;
}



// PYBIND11 modules
PYBIND11_MODULE(MatrixOperations, m) {
    m.doc() = "pybind11 matrix operations plugin"; // optional module docstring
    m.def("det2", &Determinant2x2, "Determinant of a 2x2 matrix");
    m.def("det3", &Determinant, "Determinant of a NxN matrix (N>2)");
    m.def("inv2", &Inverse2x2, py::return_value_policy::reference, "Inverse of a 2x2 matrix");
    m.def("inv3", &Inverse, py::return_value_policy::reference, "Inverse of a NxN matrix (N>2)");
}
