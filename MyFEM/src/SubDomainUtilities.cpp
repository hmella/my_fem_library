#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

namespace py = pybind11;

// Midpoints estimation
Eigen::MatrixXd midpoints(const Eigen::MatrixXd &points,
                          const Eigen::MatrixXi &connectivity) {

  // Number of points per entitie
  const int n = connectivity.row(0).size();
  
  // Number of entitites
  const int m = connectivity.col(0).size();

  // Midpoints
  Eigen::MatrixXd M(m, 3);
  Eigen::MatrixXd c(n, 3);

  // Calculate midpoints
  for (int i=0; i<m; i++) {
    // Set midpoint to zero
    c.setZero();

    // Mean coordinate
    for (int j=0; j<n; j++) {
      c.row(j) = points.row(connectivity(i,j));
    }
    M.row(i) = c.colwise().sum()/n;

  }
   
  return M;

}


// Domain dictionary construction
std::map<int, std::vector<int>> build_domain(const std::vector<int> domains,
                                             const std::vector<int> physical_markers) {

  // Number of domains
  const int n = domains.size();

  // Number of entities
  const int m = physical_markers.size();

  // Dictionary with markers and entities
  std::map<int, std::vector<int>> domain_entities;

  // Create vector of vectors
  for (int i=0; i<n; i++) {
    // Entities of domain i
    std::vector<int> v;
    domain_entities[domains[i]] = v;
  }
  
  // Fill dictionary
  for (int j=0; j<m; j++) {
    domain_entities[physical_markers[j]].push_back(j);
  }
  
  
  return domain_entities;

}




PYBIND11_MODULE(SubDomainUtilities, m) {
    m.doc() = "Utilities for Domain and SubDomain classes"; // optional module docstring
    m.def("midpoints", &midpoints, py::return_value_policy::reference, "Midpoint of cells");
    m.def("build_domain", &build_domain, py::return_value_policy::reference, "Find entities in subdomains marked as physical group in gmsh mesh");
}
