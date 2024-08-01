import os, sys
file_dir = os.path.dirname(os.path.realpath(__file__))

sys.path.append(os.path.join(file_dir,'..','..')) # Path to Aavritti directory

from    Essentials.settings    import  Settings
from    petsc4py    import  PETSc   as  PETSc
from    slepc4py    import  SLEPc   as  SLEPc
import  numpy       as      np
import  scipy
from    scipy.io   import  savemat

number_of_modes = 10

if(sys.argv[1] == 'True'):      # Post processing flag.
    is_post_process = True
else:
    is_post_process = False

# Shift for eigenvalue solver
shift =  np.pi**2/100
settings = Settings(mesh_folder   ='mesh',
                    mesh_name     = 'rectangle',
                    result_folder = 'result',
                    result_name   = 'result',
                    boundary      = {"inlet":2,"outlet":3,"walls":4}
                    )

settings.constants ={"c0":1}

# Run with FEniCS anaconda Environment
if(str(PETSc.ScalarType) == "<class 'numpy.float64'>"):
    print("Real PETSc/FEniCS project.")

    from Executeables.Boundary_Conditions.boundary import ClassicalProblem
    from Essentials.geometry import Geometry
    from Essentials.io import Io
    from dolfin import *

    # Create discretisation of Mesh
    geometry = Geometry(settings)

    # Create equation here : imported from boundary file
    equation = ClassicalProblem(settings=settings,geometry=geometry)
    
    # To build and assemble matrices
    if(is_post_process == False):

        # Setting Boundary Condition
        equation.set_boundary_conditions(settings,geometry)

        # Build and Assemble the matrices
        equation.build_Stiffness_matrix()
        equation.build_Mass_matrix()

        L = PETScMatrix()
        A = PETScMatrix()

        print("Assemble Stiffness matrix. \n")
        assemble(equation.L_weakform, tensor = L)

        print("Assemble Mass matrix. \n")
        assemble(equation.A_weakform, tensor = A)

        [bc.apply(L) for bc in equation.bcs_homogeneous]
        [bc.apply(A) for bc in equation.bcs_homogeneous]

        L_pet  = as_backend_type(L).mat()
        A_pet  = as_backend_type(A).mat()

        # Save in scipy sparse format.
        ai, aj, av = A_pet.getValuesCSR()
        A = scipy.sparse.csr_matrix((av, aj, ai))
        scipy.sparse.save_npz(settings.result_folder + '/A', A)

        li, lj, lv = L_pet.getValuesCSR()
        L = scipy.sparse.csr_matrix((lv, lj, li))
        scipy.sparse.save_npz(settings.result_folder + '/L', L)

    if(is_post_process == True):

        print("Post processing ...")
        for iEv in range(number_of_modes):

            # Load direct and adjoint eigenmodes.
            evR_np = np.load(settings.result_folder + '/ev_' + str(iEv) + '.npz')["evR"]
            evL_np = np.load(settings.result_folder + '/ev_' + str(iEv) + '.npz')["evL"]

            equation.p_real.vector().set_local(np.real(evR_np[equation.space.dofmap().dofs()]))
            equation.p_imag.vector().set_local(np.imag(evR_np[equation.space.dofmap().dofs()]))
            equation.pH_real.vector().set_local(np.real(evL_np[equation.space.dofmap().dofs()]))
            equation.pH_imag.vector().set_local(np.imag(evL_np[equation.space.dofmap().dofs()]))


            fields_to_write = {}

            fields_to_write["p_real"]   = equation.p_real
            fields_to_write["p_imag"]   = equation.p_imag
            fields_to_write["pH_real"]  = equation.pH_real
            fields_to_write["pH_imag"]  = equation.pH_imag 

            # Write response.
            io = Io()
            io.write_paraview(geometry, settings, "mode=" + str(iEv), fields_to_write)

# Run with PETSc environment
if(str(PETSc.ScalarType) == "<class 'numpy.complex128'>"):
    print("Complex PETSc.")

    if(is_post_process == False):
        from LinearModule.ev_solver import Ev_solver
        ev_solver = Ev_solver()
        ev_solver.solve(settings = settings, shift = shift)
        ev_solver.write_result(settings = settings)

