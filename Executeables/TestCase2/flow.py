import os, sys
file_dir = os.path.dirname(os.path.realpath(__file__))

sys.path.append(os.path.join(file_dir,'..','..')) # Path to Aavritti directory

from    Essentials.settings    import  Settings
from    petsc4py    import  PETSc   as  PETSc
from    slepc4py    import  SLEPc   as  SLEPc
import  numpy       as      np
import  scipy
from    scipy.io   import  savemat

if(sys.argv[1] == 'True'):      # Post processing flag.
    is_post_process = True
else:
    is_post_process = False

# Shift for eigenvalue solver
shift =  1*np.pi**2/(0.75**2)
settings = Settings(mesh_folder   ='mesh',
                    mesh_name     = 'duct',
                    result_folder = 'result',
                    result_name   = 'result',
                    boundary      = {"inlet":1,"outlet":2,"walls":3}
                    )

settings.constants ={
                    "u_in"   :1.0,# inlet velocity
                    "u_out"  :1.0 # outlet velocity
                     }

# Run with FEniCS anaconda Environment
if(str(PETSc.ScalarType) == "<class 'numpy.float64'>"):
    print("Real PETSc/FEniCS project.")

    from Equations.flow_solver import Potential
    from Essentials.geometry import Geometry
    from Essentials.io import Io
    from dolfin import *

    # Create discretisation of Mesh
    geometry = Geometry(settings)

    # Create equation here : imported from boundary file
    equation = Potential(settings=settings,geometry=geometry)
    
    # To build and assemble matrices
    if(is_post_process == False):

        # Setting Boundary Condition
        equation.set_boundary_conditions(settings,geometry)

        # Build and Assemble the matrices
        equation.build_Stiffness_matrix()
        equation.bulid_Forcing_tuple()

        L = PETScMatrix()

        print("Assemble Stiffness matrix.")
        assemble(equation.L_weakform, tensor = L)

        print("Assemble Mass matrix.")
        f = assemble(equation.f_weakform)

        L_pet  = as_backend_type(L).mat()

        # Save in scipy sparse format.
        li, lj, lv = L_pet.getValuesCSR()
        L = scipy.sparse.csr_matrix((lv, lj, li))
        scipy.sparse.save_npz(settings.result_folder + '/L', L)

        # Save Forcing as npz file
        np.savez(settings.result_folder + '/f', f = f.get_local())

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
        L = scipy.sparse.load_npz(settings.result_folder + '/L.npz')
        f = np.load(settings.result_folder + '/f.npz')["f"]

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
        phi = phi_pet.getArray()
        np.savez(settings.result_folder + '/phi.npz', phi = phi)

    if(is_post_process == True):
        print('postprocessing...')

        phi = np.load(settings.result_folder + '/phi.npz')["phi"]
        equation.phi.vector().set_local(np.real(phi[equation.space.dofmap().dofs()]))
        equation.vel=  project(grad(equation.phi), equation.vspace, solver_type="cg", preconditioner_type="ilu")

        fields_to_write = {}
        fields_to_write["phi"]        = equation.phi
        fields_to_write["velocity"]   = equation.vel

        # Write response.
        io = Io()
        io.write_paraview(geometry, settings, "solution", fields_to_write)
        
