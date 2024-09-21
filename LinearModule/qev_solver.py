# Author  : Vighnesh JR
# Contact : jrvigh@gmail.com

from petsc4py import PETSc
from slepc4py import SLEPc

import numpy as np
import scipy.sparse
from scipy.io import savemat

class QEV_solver:
    '''
    Quadratic Eigenvalue Solver
    ===========================
    Class for solving special case of polynomial eigenvalue problem (quadratic) of the form 

    [L + ωB + ω²A]P = 0

    It is useful in cases involving 
    - Reflective Boundary Conditions
    - Mean Flow field

    Instructions
    ------------
    - Name the matrices as `L.npz` , `A.npz` and `B.npz`
    - Generalised solver, complex valued terms must be present within A,B and L, No further assembly will be done 
    '''

    def __init__(self):
        #check for complex build in PETSc
        print(PETSc.ScalarType)

    def __setup_General_Quadratic_PEP(self,dim=10):

        # Problem_type is a  general polyomial problem 
        self.QES.setProblemType(SLEPc.PEP.ProblemType.GENERAL)
        self.QES.setType(SLEPc.PEP.Type.TOAR)
        self.QES.setDimensions(dim, PETSc.DEFAULT, PETSc.DEFAULT) # Minimum no of eigenvalues to be computed

        # Setting the shift value for ω
        self.QES.setTarget(self.shift)

        # Essential for Non hermitian problems
        self.QES.setWhichEigenpairs(SLEPc.EPS.Which.TARGET_MAGNITUDE)
        self.QES.setTolerances(1e-9,50) # Before 200
        st = self.QES.getST()
        st.setType('sinvert')
        st.setShift(self.shift)
        ksp = st.getKSP()
        ksp.setType('preonly')
        pc = ksp.getPC()
        pc.setType('lu')
        pc.setFactorSolverType('mumps')

        # Modifications and review
        self.QES.setFromOptions

    def solve(self,settings,shift):

        # Set shift 
        self.shift = shift
        
        # Load matrices.
        L = scipy.sparse.load_npz(settings.result_folder + '/L.npz')
        A = scipy.sparse.load_npz(settings.result_folder + '/A.npz')
        B = scipy.sparse.load_npz(settings.result_folder + '/B.npz')

        # Convert to PETSc matrices 
        L_pet = PETSc.Mat().createAIJ(size=L.shape, csr=(L.indptr, L.indices, L.data))
        A_pet = PETSc.Mat().createAIJ(size=A.shape, csr=(A.indptr, A.indices, A.data))
        B_pet = PETSc.Mat().createAIJ(size=B.shape, csr=(B.indptr, B.indices, B.data))
        
        # Format of definition of Problem: [L + ωB + ω²A]P = 0, need to solve for ω
        Mat   = [L_pet,B_pet,A_pet]

        # quadratic eigen-problem solver object creation
        self.QES = SLEPc.PEP().create()

        # Inputing the matrices to be solved in the order [L,B,A]
        self.QES.setOperators(Mat)

        # tuning the settings for default dim = 10
        self.__setup_General_Quadratic_PEP()
        print("Solve for shift = " + str(shift))
        self.QES.solve()
        

        # Access solution.
        tol, maxit = self.QES.getTolerances()
        nconv      = self.QES.getConverged()
                    
        print("Stopping condition: tol=%.4g, maxit=%d" % (tol, maxit))
        print("Number of converged eigenpairs %d" % nconv)
        evR_np = []    
        self.spectrum = np.zeros(nconv, dtype=complex)
        self.eig_error = np.zeros(nconv)    
        if nconv > 0:
        
            # Create vectors of dimensions of matrix.
            evR_real, evR_imag = L_pet.getVecs()  
        

            print()
            print("        k          ||Ax-kx||/||kx|| ")
            print("----------------- ------------------")
                
            for iEv in range(nconv):
            
                k     = self.QES.getEigenpair(iEv,evR_real,evR_imag) #extracting eigenvalues and eigenvectors
                error = self.QES.getErrorEstimate(iEv) #extracting numerical errors
                        
                self.spectrum[iEv]  = k
                self.eig_error[iEv] = error
            
                if k.imag != 0.0:
                    print(" %9f%+9f j %12g" % (k.real, k.imag, error))
                else:
                    print(" %12f      %12g" % (k.real, error))
            
                evR_np.append(evR_real.getArray() + (1j)*evR_imag.getArray()) #making a list of eigenvectors 


            # Save spectrum of case.
            np.savez(settings.result_folder + '/spectrum', spectrum = self.spectrum)
            matlabdic = {'spectrum':self.spectrum}
            savemat(settings.result_folder + '/spectrum.mat',matlabdic)


#Writes down the eigenvalues in a tabulated human readable form
    def write_result(self,settings):
        ks = self.spectrum
        es = self.eig_error
        # Open the text file in write mode
        with open(settings.result_folder + '/eigenvalues.txt', 'w') as file:
        # Write some text to the file
            file.write("\n")
            file.write("        k          ||Ax-kx||/||kx|| \n")
            file.write("----------------- ------------------\n")
            for i in range(np.size(ks)):
                k       = ks[i]
                error   = es[i]
                if k.imag != 0.0:
                    file.write(" %9f%+9f j %12g\n" % (k.real, k.imag, error))
