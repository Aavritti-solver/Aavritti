'''
Equation
========
Basic template for all the FEM assembly modules
The solver modules are inherited from `Equation(settings,geometry)` class

-----------------------------------------------------------------------
created by Vighnesh JR
on 29/7/2024

Disclaimer : The structure is a borrowed from dynX solver
'''
class Equation:

    # Initialisation
    def __init__(self, settings, geometry):
        print("Constructor not implemented!")

    def init_function_spaces(fem):
        print("Initalisation of function spaces not implemented!")

    # Bulding FEM matrices
    # * Stiffness matrix = ∫(∇u).(∇v) dV
    # * Boundary matrix = ∮ u(∇v).ds
    # * Mass matrix = ∫uv dV
    # * Forcing Tuple = ∫uQdV
    def build_Boundary_matrix(self):
        print("Boundary matrix not implemented!")

    def build_Stiffness_matrix(self):
        print("Stiffness matrix not implemented!")

    def build_Mass_matrix(self):
        print("Mass matrix not implemented!")

    def build_Time_step(self, dt):
        print("Time stepping operator not implemented!")
    
    def bulid_Forcing_tuple(self):
        print("Forcing tuple not implemented")

    # Setting up boundary condition
    # This has to be defined not in solver module but an executeable module
    def set_boundary_conditions(self):
        print("Boundary conditions not implemented")
