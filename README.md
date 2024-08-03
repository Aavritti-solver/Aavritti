Aavritti

Aavritti is a FEniCS-based Helmholtz solver developed at IIT Bombay by Vighnesh JR, who joined in 2021 as an Aerospace Engineering student. Derived from the Helmholtz package of the dynX module at TU Munich, Aavritti is designed to solve:

    Acoustic modes using the Helmholtz equation.
    Acoustic response towards flame (planned for future implementation).

Setup

    Download: Get the entire directory.
    Conda Environments:
        One for FEniCS dolfin and UFL packages.
        Another for PETSc and SLEPc with complex128 datatype support.
    TestCase Folder: Create with subfolders: mesh, baseflow, result.
    Boundary Conditions: Create/edit _boundary.py files in Executables/Boundary_Conditions.

Execution

    Run the executable in the FEniCS environment (post-processing disabled) to build weak form matrices and save them.
    Run the same file in the complex PETSc environment to solve the FEM problem and save the solution.
    Run the executable again in the FEniCS environment (post-processing enabled) to convert the solution to .xdmf format for ParaView visualization.

Acknowledgements

Developed under the guidance of Professor Wolfgang Polifke at TU Munich's Thermo-Fluid Dynamics group. The dynX module, from which Aavritti is derived, was developed by Philipp Brokof (MSc.) and Gregoire Varillon (PhD), who supervised the author's internship at TUM
