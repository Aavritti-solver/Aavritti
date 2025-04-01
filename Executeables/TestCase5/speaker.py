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

    # Setting Boundary Condition
    equation.set_boundary_conditions(settings, geometry)

    with open(settings.result_folder + "/impedances.txt", "w") as file:
        file.write("Frequency \t Re(Z) \t Im(Z) \n")

        for w in settings.values["omega_sweep"]:
            print(f"Processing frequency: {0.5 * w / np.pi} Hz")

            # Load pressure solution
            p = np.load(settings.result_folder + f"/p{int(0.5 * w / np.pi)}.npz")["p"]
            equation.p_real.vector().set_local(np.real(p[equation.space.dofmap().dofs()]))
            equation.p_imag.vector().set_local(np.imag(p[equation.space.dofmap().dofs()]))

            # Create interpolators
            p_real_interp = equation.p_real
            p_imag_interp = equation.p_imag
            p_real_interp.set_allow_extrapolation(True)
            p_imag_interp.set_allow_extrapolation(True)

            # Define the measurement plane (z = 0.75)
            measurement_plane_location = 0.75
            threshold = 1e-4  # Small tolerance

            # Get points inside the domain near z = 0.75
            probe_points = []
            for cell in cells(geometry.mesh):
                coords = np.array(cell.get_vertex_coordinates())  # Convert to NumPy array
                coords = coords.reshape(len(coords) // 3, 3)  # Manually reshape (each row = [x, y, z])

                if np.any(np.abs(coords[:, 2] - measurement_plane_location) < threshold):
                    for v in coords:
                        if np.abs(v[2] - measurement_plane_location) < threshold:
                            probe_points.append(v)

            probe_points = np.array(probe_points)

            if len(probe_points) == 0:
                print(f"WARNING: No valid probe points found near z = {measurement_plane_location}. Skipping this frequency.")
                continue

            print(f"Using {len(probe_points)} probe points for integration.")

            # Compute pressure at probe points
            p_measured_vals = np.array([
                p_real_interp(Point(pt)) + 1j * p_imag_interp(Point(pt))
                for pt in probe_points
            ])

            # Compute gradients using finite differences (only inside the domain)
            dp_measured_vals = []
            dx = 1e-4  # Small step size

            for pt in probe_points:
                x, y, z = pt
                dp_dx = (p_real_interp(Point(x + dx, y, z)) - p_real_interp(Point(x - dx, y, z))) / (2 * dx)
                dp_dy = (p_real_interp(Point(x, y + dx, z)) - p_real_interp(Point(x, y - dx, z))) / (2 * dx)
                dp_dz = (p_real_interp(Point(x, y, z + dx)) - p_real_interp(Point(x, y, z - dx))) / (2 * dx)
                dp_measured_vals.append(dp_dx + dp_dy + dp_dz)

            dp_measured_vals = np.array(dp_measured_vals)

            if np.all(np.abs(dp_measured_vals) == 0):
                print(f"WARNING, step skipped at w = {w}")
                file.write(f"{0.5 * w / np.pi} \t inf \t inf \n")
                continue

            # Compute impedance
            Z = 1j * w * np.sum(p_measured_vals) / (settings.constants['c0'] * np.sum(dp_measured_vals))
            file.write(f"{0.5 * w / np.pi} \t {np.real(Z)} \t {np.imag(Z)}\n")

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


