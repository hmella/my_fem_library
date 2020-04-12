from SubDomainUtilities import *

class Domain:
  ''' Build the whole domain of the mesh from physical groups
  of the gmsh mesh
  '''
  def __init__(self, mesh, build_cells=True, build_facets=False):

    # Physical information on facets
    if build_facets:
      try:
        facet_markers = mesh.mesh[3][mesh.facet_type()]['gmsh:physical']
        self.facet_domains = list(set(facet_markers))
      except:
        print(mesh.mesh[3])
        print('WARNING: There is not physical information for boundary Facets in '+str(mesh))
        facet_markers = []
        self.facet_domains = []
        pass

      # Entities in subdomains
      self.marked_facets = build_domain(self.facet_domains, facet_markers)

    # Physical information on cells
    if build_cells:
      try:
        cell_markers = mesh.mesh[3][mesh.cell_type()]['gmsh:physical']
        self.cell_domains = list(set(cell_markers))
      except:
        print('WARNING: There is not physical information for Cells in '+str(mesh))
        cell_markers = []
        self.cell_domains = []
        pass

      # Entities in subdomains
      self.marked_cells = build_domain(self.cell_domains, cell_markers)

  # Marked facets
  def facets(self):
    return self.marked_facets
  
  # Marked cells
  def cells(self):
    return self.cell_domains



# Subdomain class
class SubDomain:
  ''' Build subdomains from Python functions
  '''
  def __init__(self, condition, mesh, on_boundary=False):

    # Cells connectivity
    if on_boundary:
      number_of_entities = mesh.num_boundary_facets()
      connectivity = mesh.facets_connectivity()
    else:
      number_of_entities = mesh.num_cells()
      connectivity = mesh.cells_connectivity()

    # Cells midpoint
    points  = mesh.vertex_coordinates()
    mpoints = midpoints(points, connectivity)

    # Find subdomain entities
    subdomain_entities = []
    for i in range(number_of_entities):
      if condition(mpoints[i]):
        subdomain_entities.append(i)

    # Cells in subdomain
    self.subdomain_cells = subdomain_entities

  # Get cells in subdomain
  def cells(self):
    return self.subdomain_cells
