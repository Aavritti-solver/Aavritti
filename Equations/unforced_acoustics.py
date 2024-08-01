# Author  : Vighnesh JR
# contact : jrvigh@gmail.com
# -------------------------------------------------------------------------

import os, sys
file_dir = os.path.dirname(os.path.realpath(__file__))
'''
Unforced Acoustics
==================
Module that sets up FEM problem for unforced eigenvalue problems
with various unforced boundary conditions

created by Vighnesh J.R. on 29/7/2024

Contents
---------------------
- `UnforcedAcoustics1` : Uniform base field 
- `UnforcedAcoustics2` : Non Uniform base field 
- `UnforcedAcoustics3` : Non Uniform base field with freqency dependent Impedances
'''

#from equation import Equation
from dolfin import *
import numpy as np

class UnforcedAcoustics1():
    '''
    Unforced Acoustics ver 1.0
    ============================
    Created : Vighnesh J.R. on 29/7/2024
    inherited from Equation 
    -----------------------------------------
    Features
    --------
    - Can handle only uniform baseflow (speed of sound as a Field)
    - Classical boundary conditions applicable
    - Impedance boundary conditions applicable
    - Named Inlets and outlets
    '''
    def __init__(self,settings,geometry):

        self.constants = settings.constants
        self.values    = settings.values

        # Initialising named boundaries
        self.inlet_boundary             = None
        self.outlet_boundary            = None

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
            self.pH_real = Function(geometry.UPn)
            self.pH_imag = Function(geometry.UPn)
        else:
            self.space = geometry.UP1

            # Function objects for various solutions
            # Main solution modes
            self.p_real  = Function(geometry.UP1)
            self.p_imag  = Function(geometry.UP1)

            # Adjoint Modes
            self.pH_real = Function(geometry.UP1)
            self.pH_imag = Function(geometry.UP1) 

        pass

        self.p_trial = TrialFunction(self.space)
        self.p_test  = TestFunction(self.space)

    # Assigning named boundaries on demand
    def set_inlet_boundary(self,boundaries):
        self.inlet_boundary = boundaries
    
    def set_outlet_boundary(self,boundaries):
        self.outlet_boundary = boundaries

    def set_reflective_boundaries(self,
                                  boundaries =[],    # Pass a list of boundaries ( numbering must match with coefficients)
                                  reflectivities =[] # Pass a list of Reflectivity coefficients ( numbering must match)
                                  ):
        
        self.reflective_boundary_list = boundaries
        self.n_refl_boundary = len(boundaries)
        Rs = np.array(reflectivities,dtype=np.complex128)
        self.impedances  = (1 + Rs)/(1 - Rs) # If reflectivity is zero then impedance is 1 --> loss of energy


    def build_Stiffness_matrix(self):

        # Jacobian matrix = -∫(∇v).(∇u)dV
        print('Building stiffness matrix UnforcedAcoustics1')
        L = -1*inner(grad(self.p_test),grad(self.p_trial))*dx
        self.L_weakform = L

    def build_Mass_matrix(self):
        
        # Mass matrix  = (1/c0²)∫(uv)dV
        print('Building Mass matrix UnforcedAcoustics1')
        A = (1/self.constants["c0"]**2)*inner(self.p_test,self.p_trial)*dx
        self.A_weakform = A

    def build_Boundary_matrix(self):
        c0 = self.constants["c0"]

        # Boundary matrix = i/(Z_n*c0) ∮ uv ds_n
        print('Building boundary matrices (Damping)')

        # Initialisation (Need to figure out if there is a better way without collapse)
        B_real = 1*inner(self.P_test,self.p_trial) - 1*inner(self.P_test,self.p_trial)
        B_imag = 1*inner(self.P_test,self.p_trial) - 1*inner(self.P_test,self.p_trial)

        if(self.n_refl_boundary==0):
            pass
        else:
            for n in range(self.n_refl_boundary):
                B_real += np.real(1j/c0*self.impedances[n])*inner(self.p_test,self.p_trial)*self.reflective_boundary_list[n]
                B_imag += np.imag(1j/c0*self.impedances[n])*inner(self.p_test,self.p_trial)*self.reflective_boundary_list[n]

        self.B_real_weakform = B_real
        self.B_imag_weakform = B_imag

class UnforcedAcoustics2(UnforcedAcoustics1):
    '''
    Unforced Acoustics ver 2.0
    ============================
    details : To be updated
    inherited from version 1.0 

    New Features
    --------
    - Can handle non uniform baseflow
    '''
    pass
class UnforcedAcoustics3(UnforcedAcoustics2):
    '''
    Unforced Acoustics ver 3.0
    ============================
    details : To be updated
    inherited from version 2.0

    New Features
    --------
    - Can handle frequncy dependent impedances
    - Needs non linear eigenvalue solver
    '''
    pass















    

        