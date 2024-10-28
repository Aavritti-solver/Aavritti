# Author  : Vighnesh JR
# contact : jrvigh@gmail.com
# -------------------------------------------------------------------------

'''
Aero Acoustics
==================
Module that sets up FEM problem for aerodynamically excited 
acoustic problems. 

created by Vighnesh J.R. on 15/10/2024

Features
---------
- `FWH_Acoustics()`
'''
#from equation import Equation
from dolfin import *
import numpy as np

class CompactOrifice():

    # TODO : Read Noirays paper on how flow across the  orifice generate sound and implement the model 
    # TODO : Bring the above model to navier stokes framework and do FEM formulation
    pass

class FWH_Acoustics():
    '''
    Ffwocks William and Hawking Model
    =================================
    simplified model for aeroacoustics of moving surfaces
    '''
    def __init__(self,settings,geometry):

        self.constants = settings.constants
        self.values    = settings.values

        # Initialising named boundaries
        self.open_air_boundary          = None
        self.wall_boundary              = None
        self.source_surface             = None

        # Boundaries where reflectivity is specified
        self.reflective_boundary_list   = None
        self.n_refl_boundary            = None

        # Initialise functionspaces
        self.init_function_spaces(geometry, settings)


    def init_function_spaces(self,geometry,settings):

        # Check if user requires special elements
        if(settings.is_high_order):
            self.space = geometry.UPn

            # Function objects for various solutions
            # Main solution modes
            self.p_real  = Function(geometry.UPn) 
            self.p_imag  = Function(geometry.UPn)

            # Adjoint Modes
            #self.pH_real = Function(geometry.UPn)
            #self.pH_imag = Function(geometry.UPn)
        else:
            self.space = geometry.UP2

            # Function objects for various solutions
            # Main solution modes
            self.p_real  = Function(geometry.UP2)
            self.p_imag  = Function(geometry.UP2)

            # Adjoint Modes
            #self.pH_real = Function(geometry.UP2)
            #self.pH_imag = Function(geometry.UP2)

        # trial and test solutions
        self.p_trial = TrialFunction(self.space)
        self.p_test  = TestFunction(self.space)


    # Assigning named boundaries on demand
    def set_open_air_boundary(self,boundaries):
        self.open_air_boundary = boundaries
    
    def set_wall_boundary(self,boundaries):
        self.wall_boundary = boundaries

    # Non permeable Ffwocks Hawking surface
    def set_source_surface(self,boundaries):
        self.source_surface = boundaries

    def build_Stiffness_matrix(self):

        # Jacobian matrix = -∫(∇v).(∇u)dV
        print('Building stiffness matrix :Ffwocks Hawking')
        L = -1*inner(grad(self.p_test),grad(self.p_trial))*dx
        self.L_weakform = L

    def build_Mass_matrix(self):
    
        # Mass matrix  = (1/c0²)∫(uv)dV
        print('Building Mass matrix :Ffwocks Hawking')
        A = ((self.values["omega_forcing"]/self.constants["c0"])**2)*inner(self.p_test,self.p_trial)*dx
        self.A_weakform = A
    
    def build_Boundary_matrix(self):
        c0      = self.constants['c0']
        omega_f = self.values["omega_forcing"]

        # Boundary matrix = i/(Z_n*c0) ∮ uv ds_n
        print('Building boundary matrices (Damping) :Ffwocks Hawking')

        B = 1*inner(self.p_test,self.p_trial)*dx - 1*inner(self.p_test,self.p_trial)*dx
        B = Constant(omega_f/c0)*inner(self.p_test,self.p_trial)*self.open_air_boundary
        self.B_weakform =B

    def bulid_Forcing_tuple(self):

        f  = self.values["Forcing_magnitude"]
        print('Build boundary forcing tuple :Ffwocks Hawking')
        F = self.p_test*Constant(f)*self.source_surface
        self.F_weakform = F