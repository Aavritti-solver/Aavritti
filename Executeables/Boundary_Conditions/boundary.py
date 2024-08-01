'''
Editable files

specify your boundary conditions here
'''

from Equations.unforced_acoustics import UnforcedAcoustics1
from dolfin import *

class ClassicalProblem(UnforcedAcoustics1):
    '''
    With inlet, oulet and walls

    No fancy reflectivity and stuff 
    '''
    def set_boundary_conditions(self,settings,geometry):
        BCmethod = "geometric"
        self.bcs_homogeneous = [

            # Inlet and Outlet
            DirichletBC(self.space, Constant(0.0),	geometry.boundaries, 	settings.boundary["inlet"],	    BCmethod),
            DirichletBC(self.space, Constant(0.0),	geometry.boundaries, 	settings.boundary["outlet"],	BCmethod)

            # Eat 5 star and do nothing for the walls
        ] 