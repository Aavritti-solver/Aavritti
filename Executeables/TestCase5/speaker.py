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
                    mesh_name     = 'tube',
                    result_folder = 'result',
                    result_name   = 'result',
                    boundary      = {"orifice":1,"walls":3,"speaker":2,"Domain":4}
                    )

settings.constants ={"c0":352.0}
settings.values = { 
                   "Forcing_magnitude" : 20.0,
                   "omega_forcing"     : 180,
                   "omega_sweep"       : np.linspace(100,300,21)*np.pi*2
                   }
if(str(PETSc.ScalarType) == "<class 'numpy.float64'>"):
    print("Real PETSc/FEniCS project.")

    from Executeables.Boundary_Conditions.boundary import ExperimentSetupOrifice
    from Essentials.geometry import Geometry
    from Essentials.io import Io
    from dolfin import *

    # Create discretisation of Mesh
    geometry = Geometry(settings)

    # Create equation here : imported from boundary file
    equation = ExperimentSetupOrifice(settings=settings,geometry=geometry)
    
    # To build and assemble matrices
    if(is_post_process == False):

        # Setting Boundary Condition
        equation.set_boundary_conditions(settings,geometry)

        # Build and Assemble the matrices
        equation.build_Stiffness_matrix()
        equation.build_Speaker_matrix()
        equation.build_Mass_matrix()

        L = PETScMatrix()
        A = PETScMatrix()
        B = PETScMatrix()

        print("Assemble Stiffness matrix.")
        assemble(equation.L_weakform, tensor = L)
        print("Assemble Mass matrix.")
        assemble(equation.A_weakform, tensor = A)
        print("Assemble Damping matrix.")
        assemble(equation.B_speaker, tensor = B)

        [bc.apply(L) for bc in equation.bcs_homogeneous]
        [bc.apply(A) for bc in equation.bcs_homogeneous]
        [bc.apply(B) for bc in equation.bcs_homogeneous]


        L_pet  = as_backend_type(L).mat()
        A_pet  = as_backend_type(A).mat()
        B_pet =  as_backend_type(B).mat()

        print("Assemble Forcing Tuple.")
        f = assemble(equation.F_speaker)
        [bc.apply(f) for bc in equation.bcs_homogeneous]


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

    if is_post_process:
        print("Postprocessing...")

        # Set Boundary Condition
        equation.set_boundary_conditions(settings, geometry)

        with open(settings.result_folder + "/reflectivity.txt", "w") as file:
        file.write("Frequency \t R \t phi \n")

        for w in settings.values["omega_sweep"]:
            f = 0.5 * w / np.pi
            print(f"Processing frequency: {f} Hz")

            # Load pressure solution
            p = np.load(settings.result_folder + f"/solutions/p{int(f)}.npz")["p"]
            equation.p_real.vector().set_local(np.real(p[equation.space.dofmap().dofs()]))
            equation.p_imag.vector().set_local(np.imag(p[equation.space.dofmap().dofs()]))

            # Interpolators
            p_real_interp = equation.p_real
            p_imag_interp = equation.p_imag
            p_real_interp.set_allow_extrapolation(True)
            p_imag_interp.set_allow_extrapolation(True)

            # Sample along z-axis: z from 0.75 to 0.0, x = 0, y = 0
            num_points = 100
            z_vals = np.linspace(0.0, 0.75, num_points)
            x_fixed = 0.0
            y_fixed = 0.0

            pressure_mags = []
            complex_pressures = []

            for z in z_vals:
                pt = Point(x_fixed, y_fixed, z)
                p_real = p_real_interp(pt)
                p_imag = p_imag_interp(pt)
                p_complex = p_real + 1j * p_imag
                complex_pressures.append(p_complex)
                pressure_mags.append(np.abs(p_complex))

            pressure_mags = np.array(pressure_mags)
            complex_pressures = np.array(complex_pressures)

            if len(pressure_mags) == 0:
                print(f"WARNING: No valid pressure points at frequency {f} Hz")
                file.write(f"{f} \t nan \t nan \n")
                continue

            p_max = np.max(pressure_mags)
            p_min = np.min(pressure_mags)
            0 = z_vals[np.argmin(pressure_mags)]

            # Computing reflectivity magnitude and phase
            a = (p_max + p_min) / 2
            b = (p_max - p_min) / 2

            R = b / a if a != 0 else 0
            k = w / settings.constants["c0"]
            phi = np.pi - 2 * k * (0.75 - z0)

            file.write(f"{f:.6f} \t {R:.6f} \t {phi:.6f}\n")


if str(PETSc.ScalarType) == "<class 'numpy.complex128'>":
    print("Complex PETSc.")

    # Function for setting up KSP solver
    def setupKSP(P):
        P.setType('preonly')
        P.setTolerances(1e-10, 200)
        PC = P.getPC()
        PC.setType('lu')
        PC.setFactorSolverType('mumps')
        P.setFromOptions()

    # Load Matrix and Forcing
    L = scipy.sparse.load_npz(settings.result_folder + '/L.npz')
    B0 = scipy.sparse.load_npz(settings.result_folder + '/B.npz')
    A0 = scipy.sparse.load_npz(settings.result_folder + '/A.npz')
    f0 = np.load(settings.result_folder + '/f.npz')["f"]

    # Convert to PETSc format
    L_pet = PETSc.Mat().createAIJ(size=L.shape, csr=(L.indptr, L.indices, L.data))
    B_pet = PETSc.Mat().createAIJ(size=B0.shape, csr=(B0.indptr, B0.indices, B0.data))
    A_pet = PETSc.Mat().createAIJ(size=A0.shape, csr=(A0.indptr, A0.indices, A0.data))
    
    f_pet, _ = L_pet.getVecs()
    p_pet, _ = L_pet.getVecs()

    print('check 3')
    # Create KSP solver once and reuse
    LHS = PETSc.KSP().create()
    setupKSP(LHS)

    for w in settings.values["omega_sweep"]:
        print('for omega =' + str(w))

        # Reset matrices before modifying
        B_temp = B_pet.duplicate(copy=True)
        A_temp = A_pet.duplicate(copy=True)

        B_temp.scale(1j * w)
        A_temp.scale(w**2)

        # Compute new LHS matrix
        LHS.setOperators(L_pet + B_temp + A_temp)

        # Scale the forcing vector
        f_pet.setValues(np.arange(len(f0), dtype=np.int32), (1j * w * f0))
        f_pet.assemble()

        # Solve the system
        LHS.solve(f_pet, p_pet)
        pressure = p_pet.getArray()

        # Save result
        np.savez(settings.result_folder + '/p' + str(int(0.5 * w / np.pi)) + '.npz', p=pressure)


