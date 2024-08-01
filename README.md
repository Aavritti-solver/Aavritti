# Aavritti

**Aavritti** is a FEniCS-based Helmholtz solver developed at IIT Bombay. This package was created by **Vighnesh JR**, a 2021 graduate in Aerospace Engineering, inspired by his **DAAD WISE** internship with the Thermo-Fluid Dynamics group at TU Munich. Aavritti is derived from the `Helmholtz` package (made by the authour during his internship) of the dynX module developed by the group. It is designed to solve:

- Acoustic modes using the Helmholtz equation.
- Acoustic response towards flame (Linearised, planned for future implementation).

## General Use Guidelines

1. Download the entire directory to your system.
2. Set up two conda environments:
   - One that supports FEniCS dolfin and UFL packages.
   - Another that supports the complex128 datatype for PETSc and SLEPc in Python.
3. All relevant files are located in the `Executables` directory.
4. Create a *TestCase* folder with the following subfolders:
   - `mesh`
   - `baseflow`
   - `result`
5. In the *TestCase* folder, create an executable file following the provided example format. To set the boundary conditions go to `Executeables/Boundary_Conditions` and create or edit any _boundary.py_ file
6. Run the executable Python file in the FEniCS environment with post-processing disabled. This builds the weak form matrices and saves them.
7. Rerun the same file in the complex PETSc environment to solve the FEM problem and save the solution.
8. Finally, run the executable file again in the FEniCS-supported environment with post-processing enabled. This converts the solution into a viewable `.xdmf` format for visualization in ParaView.

## Additional Software

1. **ParaView**: A Linux-based freeware for visualization, compatible with FEniCS. [ParaView](https://www.paraview.org/)
2. **Gmsh**: A tool for naming domains, boundaries, and meshing. [Gmsh](https://gmsh.info/)
3. **FreeCAD**: CAD software for creating geometries. [FreeCAD](https://www.freecad.org/)

## Extra Notes

### Mesh Generation

FEniCS can read `.xml` files. Follow these steps to process the geometry:

1. Create a geometry using FreeCAD (or Gmsh) and save it in `.brep` (boundary representation) format.
2. Import the `geometry.brep` file into Gmsh, name the important physical domains and boundaries, and save it as a `geometry.gmsh` file.
3. Using Gmsh's coding interface (C++ based), specify the meshing parameters and specifications. Mesh the geometry and export the mesh as a `geometry.msh` file.
   - **Note**: The version of FEniCS used supports only msh 2-ASCII format, not the default latest msh 4-ASCII format.
   - **Note**: Note down the numbering of physical regions in the .msh file
4. Convert your `.msh` file to `.xml` file using `dolfin-convert` command in your FEniCS environment. This can be done using your linux ubuntu terminal. Example is shown below:
```bash
conda activate fenicsproject
dolfin-convert geometry.msh geometry.xml
```
5. The above step should generate the following three `.xml` files. Make sure they are generated. Otherwise there might be errors in the above steps.
   - `geometry.xml` : the base file
   - `geometry_physical_region.xml` : contains your naming of physical regions like domain and boundary or surface or line etc-
   - `geometry_facet_region.xml` : contains topological boundary details,
### Conda Environment setup

## Acknowledgements

This project was developed under the guidance of Professor Wolfgang Polifke at TU Munich's Thermo-Fluid Dynamics group. The dynX module, from which Aavritti is derived, was developed by Philipp Brokof (MSc.) and Gregoire Varillon (PhD), who supervised the author's internship at TUM.

- [Thermo-Fluid Dynamics group, TU Munich](https://www.epc.ed.tum.de/en/tfd/home/)

