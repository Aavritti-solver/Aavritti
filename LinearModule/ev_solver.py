from petsc4py import PETSc
from slepc4py import SLEPc

import numpy as np
import scipy.sparse
from scipy.io import savemat

class Ev_solver:
        
        def __init__(self):
             
            # Check for complex build of petsc.
            print(PETSc.ScalarType)     

        # Settings for linear eigenvalue problem.
        def setupGeneralizedNonHermitianEPS(self):
            self.LEP.setProblemType(SLEPc.EPS.ProblemType.GNHEP)
            self.LEP.setType(SLEPc.EPS.Type.KRYLOVSCHUR)             # Use KRYLOVSHUR to get an idea of the spectrum.
            self.LEP.setDimensions(10, PETSc.DEFAULT, PETSc.DEFAULT) # 2 only for testing --> production 10!
            self.LEP.setTarget(self.shift)
            self. LEP.setTwoSided(True)
            #self.LEP.setWhichEigenpairs(SLEPc.EPS.Which.TARGET_REAL)   # Target real for scanning??
            #self.LEP.setWhichEigenpairs(SLEPc.EPS.Which.TARGET_IMAGINARY)
            self.LEP.setWhichEigenpairs(SLEPc.EPS.Which.TARGET_MAGNITUDE)    # Real target seems to work better for Krylov Schur.
            self.LEP.setTolerances(1e-9,50) # Before 200
        
            # Setting up shift invert.
            st = self.LEP.getST()
            st.setType('sinvert')
            st.setShift(self.shift)
            ksp = st.getKSP()
            ksp.setType('preonly')
            pc = ksp.getPC()
            pc.setType('lu')
            pc.setFactorSolverType('mumps')

            self.LEP.setFromOptions
            self.LEP.view()

        def solve(self, settings, shift):

            # Set shift.
            self.shift = shift
            
            # Load matrices.
            L = scipy.sparse.load_npz(settings.result_folder + '/L.npz')
            A = scipy.sparse.load_npz(settings.result_folder + '/A.npz')

            # Important: Here equations are formulatet as (\sigma A + L)q = 0.
            # However, slepc expects Lq = \sigma A q.
            # => Multiply A matrix with -1, before stability analysis.
            A = A * (-1)

            # Conver to petsc matrices.
            L_pet = PETSc.Mat().createAIJ(size=L.shape, csr=(L.indptr, L.indices, L.data))
            A_pet = PETSc.Mat().createAIJ(size=A.shape, csr=(A.indptr, A.indices, A.data))

            self.LEP = SLEPc.EPS().create()
            self.LEP.setOperators(L_pet, A_pet)
            self.setupGeneralizedNonHermitianEPS()

            print("Solve for shift = " + str(shift))
            self.LEP.solve()

            # Access solution.
            tol, maxit = self.LEP.getTolerances()
            nconv      = self.LEP.getConverged()
                    
            print("Stopping condition: tol=%.4g, maxit=%d" % (tol, maxit))
            print("Number of converged eigenpairs %d" % nconv)
            
            self.spectrum = np.zeros(nconv, dtype=complex)
            self.eig_error = np.zeros(nconv)
            
            if nconv > 0:
            
                # Create vectors of dimensions of matrix.
                evR, evL = L_pet.getVecs()  

                print()
                print("        k          ||Ax-kx||/||kx|| ")
                print("----------------- ------------------")
                
                for iEv in range(nconv):
            
                    k     = self.LEP.getEigenvalue(iEv)
                    error = self.LEP.getErrorEstimate(iEv)
                        
                    self.LEP.getEigenvector(iEv, evR)
                    self.LEP.getLeftEigenvector(iEv, evL)
            
                    self.spectrum[iEv]  = k
                    self.eig_error[iEv] = error
            
                    if k.imag != 0.0:
                        print(" %9f%+9f j %12g" % (k.real, k.imag, error))
                    else:
                        print(" %12f      %12g" % (k.real, error))
            
                    evR_np = evR.getArray()
                    evL_np = evL.getArray()

                    # Save eigen vector data.
                    np.savez(settings.result_folder + '/ev_' + str(iEv), evR = evR_np, evL = evL_np)

                # Save spectrum of case.
                np.savez(settings.result_folder + '/spectrum', spectrum = self.spectrum)
                matlabdic = {'spectrum':self.spectrum}
                savemat(settings.result_folder + '/spectrum.mat',matlabdic)
                
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
                    else:
                        file.write(" %12f      %12g\n" % (k.real, error))
            

