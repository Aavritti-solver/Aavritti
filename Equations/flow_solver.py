# Author  : Vighnesh JR
# contact : jrvigh@gmail.com
# -------------------------------------------------------------------------


'''
Flow Solver
===========
Module that solves base flow using simplified 
linear fluid models like potential flow theory

created by Vighnesh J.R. on 12/9/2024

contents
--------
- `Potential` : solves simple potential flow field
'''


import os, sys
file_dir = os.path.dirname(os.path.realpath(__file__))

#from equation import Equation
from dolfin import *
import numpy as np

class Potential():

    def __init__(self,settings,geometry):

        self.constants = settings.constants
        self.values    = settings.values

        # Initialising named boundaries
        self.inlet_boundary             = None
        self.outlet_boundary            = None

        # Initialise functionspaces
        self.init_function_spaces(geometry, settings)
    
    def init_function_spaces(self,geometry,settings):

        # Check if user requires special elements
        if(settings.is_high_order):
            self.space = geometry.UPn
            self.vspace = geometry.VPn

            # Function objects for various solutions
            # Main solution modes
            self.phi  = Function(geometry.UPn) 

            # Velocity field
            self.vel = Function(geometry.VPn)

        else:
            self.space = geometry.UP2
            self.vspace = geometry.VP2

            # Function objects for various solutions
            # Main solution modes
            self.phi  = Function(geometry.UP2)

            # Velocity field
            self.vel = Function(geometry.VP2)

        pass

        self.phi_trial = TrialFunction(self.space)
        self.phi_test  = TestFunction(self.space)

        # Assigning named boundaries on demand
    def set_inlet_boundary(self,boundaries):
        self.inlet_boundary = boundaries
    
    def set_outlet_boundary(self,boundaries):
        self.outlet_boundary = boundaries

    def build_Stiffness_matrix(self):
        # Jacobian matrix = -∫(∇v).(∇u)dV
        print('Building stiffness matrix : Potential')
        L = -1*inner(grad(self.phi_test),grad(self.phi_trial))*dx
        self.L_weakform = L

    def bulid_Forcing_tuple(self):
        # Forcing tuple = ∫v(u_out*ds_out) - ∫v(u_in*ds_in)
        print('Build boundary forcing tuple : Potential')
        f = self.phi_test*self.constants['u_out']*self.outlet_boundary - self.phi_test*self.constants['u_in']*self.inlet_boundary
        self.f_weakform = f


    def set_boundary_conditions(self,settings,geometry):
        
        # eat 5 star and do nothing for the walls
        self.set_outlet_boundary(boundaries = geometry.ds(settings.boundary["outlet"]))
        self.set_inlet_boundary(boundaries = geometry.ds(settings.boundary["inlet"]))