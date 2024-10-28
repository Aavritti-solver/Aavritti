# Author  : Vighnesh JR
# Contact : jrvigh@gmail.com
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

if(sys.argv[2] == 'True'):      # Flow solver flag
    is_compute_flow = True
else:
    is_compute_flow = False

# Shift for eigenvalue solver
shift =  1*np.pi/1.0
settings = Settings(mesh_folder     ='mesh',
                    mesh_name       = 'duct',
                    result_folder   = 'result',
                    result_name     = 'result',
                    baseflow_folder = 'base_flow',
                    baseflow_name   = '/phi',
                    boundary        = {"inlet":1,"outlet":2,"walls":3}
                    )

settings.constants ={
                    "u_in"   :100.0, # inlet velocity
                    "u_out"  :100.0, # outlet velocity
                    "c0"     :352  # Speed of Sound
                     }

if(str(PETSc.ScalarType) == "<class 'numpy.float64'>"):
    print("Real PETSc/FEniCS project.")

    from Equations.flow_solver import Potential
    from Executeables.Boundary_Conditions.boundary import ClassicalMeanFlowProblem
    from Essentials.geometry import Geometry
    from Essentials.io import Io
    from dolfin import *

    # Create discretisation of Mesh
    geometry = Geometry(settings)

    # Object to complie matrices for flow solver
    Flow = Potential(settings   = settings,
                    geometry   = geometry)
    
    # Object to complie matrices for acoustic problems
    Sound = ClassicalMeanFlowProblem(settings    = settings,
                                    geometry    = geometry)
    # To build and assemble matrices
    if(is_post_process == False):
        if(is_compute_flow == True):

            # Create Flow here : imported from Sounds
            

            # Setting Boundary Condition
            Flow.set_boundary_conditions(settings,geometry)

            # Build and Assemble the matrices
            Flow.build_Stiffness_matrix()
            Flow.bulid_Forcing_tuple()

            L = PETScMatrix()

            print("Assemble Stiffness matrix.")
            assemble(Flow.L_weakform, tensor = L)

            print("Assemble Mass matrix.")
            f = assemble(Flow.f_weakform)

            L_pet  = as_backend_type(L).mat()

            # Save in scipy sparse format.
            li, lj, lv = L_pet.getValuesCSR()
            L = scipy.sparse.csr_matrix((lv, lj, li))
            scipy.sparse.save_npz(settings.baseflow_folder + '/L', L)

            # Save Forcing as npz file
            np.savez(settings.baseflow_folder + '/f', f = f.get_local())

            # Solving the problem ( we dont need to use complex petsc)
            print('\n solving flow ...')

            # Function for setting up KSP solver
            def setupKSP(P):
                P.setType('preonly')
                P.setTolerances(1e-10, 200)
                PC = P.getPC(); PC.setType('lu')
                PC.setFactorSolverType('mumps')
                P.setFromOptions()
            
            # Load matrix and forcing (bit redundant step, need to improve)
            L = scipy.sparse.load_npz(settings.baseflow_folder + '/L.npz')
            f = np.load(settings.baseflow_folder + '/f.npz')["f"]

            # Convert to PETSc format and make LHS,RHS and solution holder'
            L_pet = PETSc.Mat().createAIJ(size=L.shape, csr=(L.indptr, L.indices, L.data))
            f_pet,_ = L_pet.getVecs()
            f_pet.setValues(np.arange(len(f),dtype=np.int32),f)

            # Create solution vector
            phi_pet,_ = L_pet.getVecs()

            # setting up KSP object 
            LHS = PETSc.KSP().create()
            LHS.setOperators(L_pet)
            setupKSP(LHS)

            LHS.solve(f_pet,phi_pet) # L*phi = f

            # Mach no potential is saved here
            phi = phi_pet.getArray()/settings.constants['c0']

            # Base flow is saved in potential form as 'base_flow_name.npz)
            np.savez(settings.baseflow_folder + settings.baseflow_name +'.npz', phi = phi)
    
        
        
        # Loading baseflow as velocity field
        phi = np.load(settings.baseflow_folder + settings.baseflow_name +'.npz')["phi"]
        Flow.phi.vector().set_local(np.real(phi[Flow.space.dofmap().dofs()]))
        Sound.mean_flow =  project(grad(Flow.phi), Sound.vspace, solver_type="cg", preconditioner_type="ilu")

        fields_to_write ={"MeanFlow":Sound.mean_flow}

        # Write response.
        io = Io()
        io.write_paraview(geometry, settings,"base_flow", fields_to_write)

        # Assigning boundary conditions
        Sound.set_boundary_conditions(settings=settings,geometry=geometry)

        Sound.build_Stiffness_matrix()
        Sound.build_Mass_matrix()
        Sound.build_Damping_matrix()

        L = PETScMatrix()
        A = PETScMatrix()
        B = PETScMatrix() # Devoid of i

        print("Assemble Stiffness matrix.")
        assemble(Sound.L_weakform, tensor = L)

        print("Assemble Mass matrix.")
        assemble(Sound.A_weakform, tensor = A)

        print("Assemble Damping matrix.")
        assemble(Sound.B_weakform, tensor = B)
        
        [bc.apply(L) for bc in Sound.bcs_homogeneous]
        [bc.apply(A) for bc in Sound.bcs_homogeneous]
        [bc.apply(B) for bc in Sound.bcs_homogeneous]

        L_pet  = as_backend_type(L).mat()
        A_pet  = as_backend_type(A).mat()
        B_pet  = as_backend_type(B).mat()

        # Save in scipy sparse format.
        ai, aj, av = A_pet.getValuesCSR()
        A = scipy.sparse.csr_matrix((av, aj, ai))
        scipy.sparse.save_npz(settings.result_folder + '/A', A)

        li, lj, lv = L_pet.getValuesCSR()
        L = scipy.sparse.csr_matrix((lv, lj, li))
        scipy.sparse.save_npz(settings.result_folder + '/L', L)

        bi, bj, bv = B_pet.getValuesCSR()
        B = 1j*scipy.sparse.csr_matrix((bv, bj, bi))
        scipy.sparse.save_npz(settings.result_folder + '/B', B)

    if(is_post_process == True):
        print("Post processing ...")
        for iEv in range(number_of_modes):

            # Load eigenmodes.
            ev_np = np.load(settings.result_folder + '/ev_' + str(iEv) + '.npz')["ev"]
            Sound.p_real.vector().set_local(np.real(ev_np[Sound.space.dofmap().dofs()]))
            Sound.p_imag.vector().set_local(np.imag(ev_np[Sound.space.dofmap().dofs()]))

            fields_to_write = {}

            fields_to_write["p_real"]   = Sound.p_real
            fields_to_write["p_imag"]   = Sound.p_imag

            # Write response.
            io = Io()
            io.write_paraview(geometry, settings, "mode=" + str(iEv), fields_to_write)

# Run with PETSc environment
if(str(PETSc.ScalarType) == "<class 'numpy.complex128'>"):
    print("Complex PETSc.")

    if(is_post_process == False):
        from LinearModule.qev_solver import QEV_solver
        qev_solver = QEV_solver()
        qev_solver.solve(settings = settings, shift = shift)
        qev_solver.write_result(settings = settings)
