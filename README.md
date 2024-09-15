# Aavritti

Aavritti is a FEniCS-based Helmholtz solver developed at IIT Bombay by Vighnesh JR, who joined in 2021 as an Aerospace Engineering student. Derived from the Helmholtz package of the dynX module at TU Munich, Aavritti is designed to solve:

- Acoustic modes using the Helmholtz equation.
- Acoustic response towards flame (planned for future implementation).
# Latest Feature!!!
Added new potential flow solver module which will be integrated with acoustic solver for mean flow analysis.

## Additional Softwares required
1. **ParaView**: A Linux-based freeware for visualization, compatible with FEniCS. [ParaView](https://www.paraview.org/)
2. **Gmsh**: A tool for naming domains, boundaries, and meshing. [Gmsh](https://gmsh.info/)
3. **FreeCAD**: CAD software for creating geometries. [FreeCAD](https://www.freecad.org/)
## Setup
- Download: Get the entire directory.
- Conda Environments:
  1. One for FEniCS dolfin and UFL packages.
  2. Another for PETSc and SLEPc with `numpy.complex128` datatype support.
- TestCase Folder: Create with subfolders: `mesh`, `baseflow`, `result`.
- Boundary Conditions: `Create/edit _boundary.py` files in `Executables/Boundary_Conditions`.

## Execution

- Run the executable in the FEniCS environment (post-processing disabled) to build weak form matrices and save them.
- Run the same file in the complex PETSc environment to solve the FEM problem and save the solution.
- Run the executable again in the FEniCS environment (post-processing enabled) to convert the solution to .xdmf format for ParaView visualization.

## Acknowledgements

Developed under the guidance of Professor Wolfgang Polifke at TU Munich's Thermo-Fluid Dynamics group. The dynX module, from which Aavritti is derived, was developed by Philipp Brokof (MSc.) and Gregoire Varillon (PhD), who supervised the author's internship at TUM
