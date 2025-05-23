'''
Editable files

specify your boundary conditions here
'''

from Equations.unforced_acoustics import UnforcedAcoustics1
from Equations.unforced_acoustics import UnforcedAcoustics1_1
from Equations.aero_acoustics import FWH_Acoustics
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

class ClassicalMeanFlowProblem(UnforcedAcoustics1_1):
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

class BasicRotorAcoustics(FWH_Acoustics):
    '''
    - acoustically transparent ends
    - walls
    - emission surface
    '''
    def set_boundary_conditions(self,settings,geometry):

        # Assigning names
        self.set_open_air_boundary(geometry.ds(settings.boundary["air"]))
        self.set_source_surface(geometry.ds(settings.boundary["source"]))
        self.set_wall_boundary(geometry.ds(settings.boundary["wall"]))

class ExperimentSetupOrifice(SpeakerAcoustics):
    '''
    Designed for BTP2 by vighnesh J.R. for impedance measurement of orifices
    - Orifice   : p = 0
    - Walls     : (∇p).n = 0
    - speaker   : (∇p).n = i(ω/c0)p - 2iωρF
    - computed surface : where orifice impedance is measured
    '''
    def set_boundary_conditions(self,settings,geometry):
        BCmethod = "geometric"
        self.bcs_homogeneous = [

            # Orifice
            DirichletBC(self.space, Constant(0.0), geometry.boundaries, settings.boundary["orifice"],    BCmethod)

            # Eat 5 star and do nothing for the walls
        ]
        self.set_speaker_boundary(geometry.ds(settings.boundary["speaker"]))
