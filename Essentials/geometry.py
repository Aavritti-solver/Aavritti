from dolfin import *

class Geometry:
    '''
    Geometry
    ========
    created by Vighnesh JR on 27/7/2024

    Features
    --------
    1. reads mesh from directory in  'xml' format
    2. assigns nomenclature and defines domain and boundary and their corresponding differential forms
    3. creates elements and functionspaces for the mesh (and also vectorfunctionspaces)
    4. Supports mesh of 1, 2 and 3 dimensional nature
    5. One can also specify higher order elements if required 
    '''

    def __init__(self, settings,degree = None):
        self.mesh        = Mesh(settings.mesh_folder + '/' + settings.mesh_name + '.xml')
        self.domains     = MeshFunction('size_t', self.mesh, self.mesh.topology().dim())
        self.dim         = self.mesh.geometry().dim()
        self.degree      = degree
        if(self.dim > 1):
            self.boundaries  = MeshFunction('size_t', self.mesh, settings.mesh_folder + '/' + settings.mesh_name + '_facet_region.xml')
            self.ds          = Measure('ds', domain=self.mesh, subdomain_data=self.boundaries)
            self.n           = FacetNormal(self.mesh)

        self.init_elements()

    def load_base_mesh(self, settings):
        self.mesh_base        = Mesh(settings.mesh_folder + '/' + settings.mesh_base_name + '.xml')
        self.domains_base     = MeshFunction('size_t', self.mesh_base, self.mesh_base.topology().dim())
        self.dim_base         = self.mesh_base.geometry().dim()
        if(self.dim > 1):
            self.boundaries_base  = MeshFunction('size_t', self.mesh_base, settings.mesh_folder + '/' + settings.mesh_base_name + '_facet_region.xml')
            self.ds_base          = Measure('ds', domain=self.mesh_base, subdomain_data=self.boundaries_base)
            self.n_base           = FacetNormal(self.mesh_base)

        self.init_elements_base()

    def init_elements_base(self):
        cell_list = {"1":interval,
                     "2":triangle,
                     "3":tetrahedron,}
        key = str(self.dim)
        cell= cell_list[key]

        self.P1_base = FiniteElement('P',cell,1)
        self.P2_base = FiniteElement('P',cell,2)
        self.R0_base = FiniteElement('R',cell,0) # I think we dont need this element for Helmholtz solver

        self.UR0_base = FunctionSpace(self.mesh_base, self.R0_base)
        self.UP1_base = FunctionSpace(self.mesh_base, self.P1_base)
        self.UP2_base = FunctionSpace(self.mesh_base, self.P2_base)
        self.VP2_base = VectorFunctionSpace(self.mesh_base, self.P2_base,self.dim)

        # Custom Higher order element
        if(self.degree != None):
            self.Pn_base  = FiniteElement('P',cell,self.degree)
            self.UPn_base = FunctionSpace(self.mesh_base,self.Pn_base)
            self.VPn_base = VectorFunctionSpace(self.mesh_base,self.Pn_base,self.dim)

    def init_elements(self):
        cell_list = {"1":interval,
                     "2":triangle,
                     "3":tetrahedron,}
        key = str(self.dim)
        cell= cell_list[key]

        self.P1 = FiniteElement('P',cell,1)
        self.P2 = FiniteElement('P',cell,2)
        self.R0 = FiniteElement('R',cell,0) # I think this element is not needed for Helmholtz solver
        
        self.UR0 = FunctionSpace(self.mesh, self.R0)
        self.UP1 = FunctionSpace(self.mesh, self.P1)
        self.UP2 = FunctionSpace(self.mesh, self.P2)
        self.VP2 = VectorFunctionSpace(self.mesh, self.P2, self.dim)

        # Custom Higher order element
        if(self.degree != None):
            self.Pn  = FiniteElement('P',cell,self.degree)
            self.UPn = FunctionSpace(self.mesh,self.Pn)
            self.VPn = VectorFunctionSpace(self.mesh,self.Pn,self.dim)