# Author  : Vighnesh JR
# contact : jrvigh@gmail.com
# -------------------------------------------------------------------------

'''
Unforced Acoustics
==================
Module that sets up FEM problem for unforced eigenvalue problems
with various unforced boundary conditions

created by Vighnesh J.R. on 29/7/2024

Contents
---------------------
- `UnforcedAcoustics1`   : Uniform base field 
- `UnforcedAcoustics1_1` : Non Uniform base field with mean flow
- `UnforcedAcoustics1_2` : Freqency dependent Impedances
- `SpeakerAcoustics`     : Speaker forced boundary conditions
'''
import os, sys
file_dir = os.path.dirname(os.path.realpath(__file__))


#from equation import Equation
from dolfin import *
import numpy as np

###########################################################################################################################
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
            self.space = geometry.UP2

            # Function objects for various solutions
            # Main solution modes
            self.p_real  = Function(geometry.UP2)
            self.p_imag  = Function(geometry.UP2)

            # Adjoint Modes
            self.pH_real = Function(geometry.UP2)
            self.pH_imag = Function(geometry.UP2) 

        

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
        B_real = 1*inner(self.p_test,self.p_trial)*dx - 1*inner(self.p_test,self.p_trial)*dx
        B_imag = 1*inner(self.p_test,self.p_trial)*dx - 1*inner(self.p_test,self.p_trial)*dx

        if(self.n_refl_boundary==0):
            pass
        else:
            for n in range(self.n_refl_boundary):
                B_real += np.real(1j/c0*self.impedances[n])*inner(self.p_test,self.p_trial)*self.reflective_boundary_list[n]
                B_imag += np.imag(1j/c0*self.impedances[n])*inner(self.p_test,self.p_trial)*self.reflective_boundary_list[n]

        self.B_real_weakform = B_real
        self.B_imag_weakform = B_imag

 ###########################################################################################################################

class UnforcedAcoustics1_1(UnforcedAcoustics1):
    '''
    Unforced Acoustics ver 1.1
    ============================
    Created : Vighnesh J.R. on 29/7/2024
    inherited from version 1.0
    -----------------------------------
    New Features
    --------
    - Can handle non uniform baseflow that is
      - Steady
      - Incompressible mean flow
      - Low mach number
    - Cannot handle impedance boundary conditions
    '''
    def init_function_spaces(self,geometry,settings):

        # Check if user requires special elements
        if(settings.is_high_order):
            self.space = geometry.UPn
            self.vspace = geometry.VPn

            # Function objects for various solutions
            # Main solution modes
            self.p_real  = Function(geometry.UPn) 
            self.p_imag  = Function(geometry.UPn)

            # Mean flow field
            self.mean_flow = Function(geometry.VPn)

            # Adjoint Modes
            self.pH_real = Function(geometry.UPn)
            self.pH_imag = Function(geometry.UPn)
        else:
            self.space = geometry.UP2
            self.vspace = geometry.VP2

            # Function objects for various solutions
            # Main solution modes
            self.p_real  = Function(geometry.UP2)
            self.p_imag  = Function(geometry.UP2)

            # Mean flow field
            self.mean_flow = Function(geometry.VP2)

            # Adjoint Modes
            self.pH_real = Function(geometry.UP2)
            self.pH_imag = Function(geometry.UP2)

        self.p_trial = TrialFunction(self.space)
        self.p_test  = TestFunction(self.space)

    def build_Stiffness_matrix(self):

        # Jacobian matrix = ∫(∇v).M ∇.(Mu)dV -∫(∇v).(∇u)dV
        print('Building stiffness matrix UnforcedAcoustics 1.1')
        L = inner(grad(self.p_test), self.mean_flow)*div(self.mean_flow*self.p_trial)*dx - inner(grad(self.p_test), grad(self.p_trial))*dx        
        self.L_weakform = L
    
    def build_Damping_matrix(self):
        
        # Damping Matrix = 2i∫ v∇.(Mu)dV
        # Imaginary part omitted
        print('Bulding damping matrix UnforcedAcoustics 1.1')
        B = 2*inner(self.p_test,div(self.mean_flow*self.p_trial))*dx
        self.B_weakform = B
    
    def build_Mass_matrix(self):

        # ∫(uv)dV
        print('Building mass matrix UnforcedAcoustics 1.1')
        A = inner(self.p_test,self.p_trial)*dx
        self.A_weakform = A

#########################################################################################################################################

class UnforcedAcoustics1_2(UnforcedAcoustics1):
    '''
    Unforced Acoustics ver 1.2
    ============================
    details : To be updated
    inherited from version 1.0 

    New Features
    --------
    - Can handle frequncy dependent impedances
    - Needs non linear eigenvalue solver
    '''
    def build_Boundary_matrix(self):
        c0 = Constant(self.constants["c0"])

        # Boundary matrix = 1/(c0) ∮ uv ds_n
        print('Building boundary matrices (Damping)')

        # Initialisation of list and form of size nxn
        Bs = []
        B0 = 1*inner(self.p_test,self.p_trial)*dx - 1*inner(self.p_test,self.p_trial)*dx

        # Appending boundary integrals for each boundary
        if(self.n_refl_boundary==0):
            pass
        else:
            for n in range(self.n_refl_boundary):
                Bs.append((1/c0)*inner(self.p_test,self.p_trial)*self.reflective_boundary_list[n]+ B0)
        
        # This shall be saved inn result folder as per the needs of the executeable file   
        self.B_weakform_list = Bs

    # Setting boundary endowed with Frequency dependent Impedance
    def set_reflective_boundaries(self,
                                  boundaries =[],    # Pass a list of boundaries ( numbering must match with coefficients)
                                  ):
        self.reflective_boundary_list = boundaries
        self.n_refl_boundary = len(boundaries)

#########################################################################################################
class SpeakerAcoustics(UnforcedAcoustics1):
    def __init__(self,settings,geometry):

        self.constants = settings.constants
        self.values    = settings.values

        # Initialising named boundaries
        self.open_boundary              = None
        self.computed_surface           = None
        self.speaker_boundary           = None

        # Boundaries where reflectivity is specified
        self.reflective_boundary_list   = None
        self.n_refl_boundary            = None

        # Initialise functionspaces
        self.init_function_spaces(geometry, settings)

    # Setting up speaker
    def set_speaker_boundary(self,boundaries):
        self.speaker_boundary = boundaries

    # Open boundary
    def set_open_boundary(self,boundaries):
        self.open_boundary = boundaries
    
    # Computed surface (for post processing)
    def set_computed_surface(self,area):
        self.computed_surface = area

    # Building matrices and tuples
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
        c0 = Constant(self.constants["c0"])

        # Boundary matrix = i/(Z_n*c0) ∮ uv ds_n
        print('Building boundary matrices (Damping)')

        # Initialisation (Need to figure out if there is a better way without collapse)
        B_real = 1*inner(self.p_test,self.p_trial)*dx - 1*inner(self.p_test,self.p_trial)*dx
        B_imag = 1*inner(self.p_test,self.p_trial)*dx - 1*inner(self.p_test,self.p_trial)*dx

        if(self.n_refl_boundary==0):
            pass
        else:
            for n in range(self.n_refl_boundary):
                B_real += np.real(1j/c0*self.impedances[n])*inner(self.p_test,self.p_trial)*self.reflective_boundary_list[n]
                B_imag += np.imag(1j/c0*self.impedances[n])*inner(self.p_test,self.p_trial)*self.reflective_boundary_list[n]

        self.B_real_weakform = B_real
        self.B_imag_weakform = B_imag
    
    def build_Speaker_matrix(self):
        print('Building Speaker B matrix and forcing tuple')
        f = self.values["Forcing_magnitude"]
        null =  1*inner(self.p_test,self.p_trial)*dx - 1*inner(self.p_test,self.p_trial)*dx
        Bs = null + Constant(1/self.constants["c0"])*inner(self.p_test,self.p_trial)*self.speaker_boundary
        F  = self.p_test*Constant(f)*self.speaker_boundary 

        self.B_speaker = Bs
        self.F_speaker = F 
        


















    

        
