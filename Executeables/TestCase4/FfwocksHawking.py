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

settings = Settings(mesh_folder   ='mesh',
                    mesh_name     = 'rotor',
                    result_folder = 'result',
                    result_name   = 'result',
                    boundary      = {"air":1,"source":2,"wall":3}
                    )

settings.constants ={"c0":352.0}

settings.values = {
                   "Forcing_magnitude" : 20.0,
                   "omega_forcing"     : 6000.0
                   }
if(str(PETSc.ScalarType) == "<class 'numpy.float64'>"):
    print("Real PETSc/FEniCS project.")

    from Executeables.Boundary_Conditions.boundary import BasicRotorAcoustics
    from Essentials.geometry import Geometry
    from Essentials.io import Io
    from dolfin import *

    # Create discretisation of Mesh
    geometry = Geometry(settings)

    # Create equation here : imported from boundary file
    equation = BasicRotorAcoustics(settings=settings,geometry=geometry)
    
    # To build and assemble matrices
    if(is_post_process == False):

        # Setting Boundary Condition
        equation.set_boundary_conditions(settings,geometry)

        # Build and Assemble the matrices
        equation.build_Stiffness_matrix()
        equation.build_Boundary_matrix()
        equation.build_Mass_matrix()
        equation.bulid_Forcing_tuple()

        L = PETScMatrix()
        A = PETScMatrix()
        B = PETScMatrix()

        print("Assemble Stiffness matrix.")
        assemble(equation.L_weakform, tensor = L)
        print("Assemble Mass matrix.")
        assemble(equation.A_weakform, tensor = A)
        print("Assemble Damping matrix.")
        assemble(equation.A_weakform, tensor = B)

        L_pet  = as_backend_type(L).mat()
        A_pet  = as_backend_type(A).mat()
        B_pet =  as_backend_type(B).mat()

        print("Assemble Forcing Tuple.")
        f = assemble(equation.F_weakform)

        L_pet  = as_backend_type(L).mat()

        # Save in scipy sparse format.
        li, lj, lv = L_pet.getValuesCSR()
        L = scipy.sparse.csr_matrix((lv, lj, li))
        scipy.sparse.save_npz(settings.result_folder + '/L', L)

        ai, aj, av = A_pet.getValuesCSR()
        A = scipy.sparse.csr_matrix((av, aj, ai))
        scipy.sparse.save_npz(settings.result_folder + '/A', A)

        bi, bj, bv = B_pet.getValuesCSR()
        B = scipy.sparse.csr_matrix((bv, bj, bi))
        scipy.sparse.save_npz(settings.result_folder + '/B', B)

        # Save Forcing as npz file
        np.savez(settings.result_folder + '/f', f = f.get_local())

    if(is_post_process == True):
        print('postprocessing...')

        p = np.load(settings.result_folder + '/p.npz')["p"]
        equation.p_real.vector().set_local(np.real(p[equation.space.dofmap().dofs()]))
        equation.p_imag.vector().set_local(np.imag(p[equation.space.dofmap().dofs()]))

        fields_to_write = {}
        fields_to_write["p_real"]   = equation.p_real
        fields_to_write["p_imag"]   = equation.p_imag

        # Write response.
        io = Io()
        io.write_paraview(geometry, settings, "solution", fields_to_write)

if(str(PETSc.ScalarType) == "<class 'numpy.complex128'>"):
    print("Complex PETSc.")

    # Function for setting up KSP solver
    def setupKSP(P):
            P.setType('preonly')
            P.setTolerances(1e-10, 200)
            PC = P.getPC(); PC.setType('lu')
            PC.setFactorSolverType('mumps')
            P.setFromOptions()

    # Load Matrix and Forcing
    L =     scipy.sparse.load_npz(settings.result_folder + '/L.npz')
    B = 1j *scipy.sparse.load_npz(settings.result_folder + '/B.npz')
    A =     scipy.sparse.load_npz(settings.result_folder + '/A.npz')
    f = np.load(settings.result_folder + '/f.npz')["f"]

    # Convert to PETSc format and make LHS,RHS and solution holder'
    L_pet = PETSc.Mat().createAIJ(size=L.shape, csr=(L.indptr, L.indices, L.data))
    B_pet = PETSc.Mat().createAIJ(size=B.shape, csr=(B.indptr, B.indices, B.data))
    A_pet = PETSc.Mat().createAIJ(size=A.shape, csr=(A.indptr, A.indices, A.data))
    f_pet,_ = L_pet.getVecs()
    f_pet.setValues(np.arange(len(f),dtype=np.int32),f)

    # Create solution vector
    p_pet,_ = L_pet.getVecs()

    # setting up KSP object 
    LHS = PETSc.KSP().create()
    LHS.setOperators(L_pet + B_pet + A_pet)
    setupKSP(LHS)

    # Solving  [L + B + A]*p = f
    LHS.solve(f_pet,p_pet) 
    pressure = p_pet.getArray()
    np.savez(settings.result_folder + '/p.npz', p = pressure)
    