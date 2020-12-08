c ====================================================================================
c Collection of Fortran subroutines that use FDVM2 based numerical method to solve generalised Serre - Green -Naghdi equations (gSGN)
c As well as some methods to calculate the total conserved quantities. 
c ====================================================================================


c ====
c Subroutine that runs a numerical experiment with all the parameters, 
c the subroutine generates discretisation, initial conditions and then solves the initial condition problem
c it then prints out the desired outputs h,u,G at cell centers, parameter values and convergence and conservation errors.
c It also takes a list of intermediate times to print results at.
c
c Inputs:   a0 - background water depth on which travelling wave propagates
c           a1 - amplitude of travelling wave in h
c           a2 - translation speed of wave (wave is function of x - a2*t)
c           a3 - width of travelling wave
c           a4 - amplitude of travelling wave in u
c           a5 - beta parameter : speed (betas are a function of x - a5*t)
c           b1a6 - beta parameter :
c           b1a7 - beta parameter : beta1 = b1a6*(x - a5*t) + b1a7
c           b2a6 - beta parameter :
c           b2a7 - beta parameter : beta2 = b2a6*(x - a5*t) + b2a7
c           tstart - start time
c           tend - time to produce numerical solution to
c           ga - acceleration due to gravity 
c           dx - cell width
c           dt - time step size
c           n_GhstCells - number of ghost cells at each boundary
c           xbc_len - total number of cells in discretisation
c           xbc - array of cell centers from discretisation
c           hbc_init - values of h at the cell centers (xbc) at t = tstart (equivalent to cell average, as second order)
c           Gbc_init - values of G at the cell centers (xbc) at t = tstart (equivalent to cell average, as second order)
c           ubc_init - values of u at the cell centers (xbc) at t = tstart (equivalent to cell average, as second order)
c           tlist - list of times to print solutions at
c           tlist_len - length of tlist
c           ExpWdir - string containing directory to print the individual experiments results to
c           expwdirlen - length of the string ExpWdir
c Outputs:
c           currenttime - current time of numerical solution
c           hbc_fin - values of h at the cell centers (xbc) at t = currenttime (equivalent to cell average, as second order)
c           Gbc_fin - values of h at the cell centers (xbc) at t = currenttime (equivalent to cell average, as second order)
c           ubc_fin - values of h at the cell centers (xbc) at t = currenttime (equivalent to cell average, as second order)
c           Energ_Init - Array of conserved quantities (h,G,uh,Energy) at t = tstart
c           Energ_Fin - Array of conserved quantities (h,G,uh,Energy) at t = currenttime

c ===== 
      SUBROUTINE NumericalSolveForcedTSPrint(a0,a1,a2,a3,a4,a5,
     . b1a6,b1a7,b2a6,b2a7,
     . tstart,tend,ga,dx,dt,n_GhstCells,xbc_len,xbc,
     . hbc_init,Gbc_init,ubc_init,
     . tlist,tlist_len,ExpWdir,ExpWdir_len,
     . currenttime,hbc_fin,Gbc_fin,ubc_fin)
     
     
      IMPLICIT NONE
      !Inputs  
      INTEGER n_GhstCells,xbc_len,ExpWdir_len,tlist_len
      CHARACTER(LEN=ExpWdir_len) ExpWdir
      DOUBLE PRECISION a0,a1,a2,a3,a4,a5,b1a6,b1a7,b2a6,b2a7,
     . tstart,tend,ga,dx,dt,currenttime
      DOUBLE PRECISION xbc(xbc_len),hbc_init(xbc_len),
     . Gbc_init(xbc_len), ubc_init(xbc_len),
     . tlist(tlist_len)
     
      !Outputs  
      DOUBLE PRECISION hbc_fin(xbc_len),Gbc_fin(xbc_len),
     . ubc_fin(xbc_len)
      
      !Local variables - Loop variable i, filecount and strct ensure that intermediate time prints write to different files     
      INTEGER i,filecount
      CHARACTER(LEN=2) strct
      
      !initial time
      currenttime  = tstart
      filecount = 1
      
      !Set current values of h,u,G to initial conditions (intialise)
      DO i = 1,xbc_len
         hbc_fin(i) = hbc_init(i) 
         Gbc_fin(i) = Gbc_init(i) 
         ubc_fin(i) = ubc_init(i)
      END DO
      
      !Evolve system through time
      DO WHILE (currenttime  .LT. tend )      
         !Determine if we should print results which happens - at initial time, and when currenttime is just before an intermediate time
         IF ((currenttime + dt .GE. tlist(filecount)) 
     .      .OR. (filecount .EQ. 1 ))  THEN
            
            ! Get unique name for print file, then open and print t,x,h,G,u  
            WRITE (strct,'(I2)') filecount
            OPEN(8, FILE = ExpWdir// strct //'.dat') 
            DO i = 1,xbc_len
               WRITE(8,*) currenttime,xbc(i),hbc_fin(i),
     .            Gbc_fin(i),ubc_fin(i)
            END DO
            CLOSE(8)
            
            filecount = filecount + 1
         END IF

         !Evolve system
         CALL EvolveStepWrap(a0,a1,a2,a3,a4,a5,
     .    b1a6,b1a7,b2a6,b2a7,currenttime,xbc,
     .    xbc_len,n_GhstCells,ga,
     .    dx,dt,hbc_fin,Gbc_fin,ubc_fin)
         !Update current time and print, so we can get some feedback about what its doing
         currenttime  = currenttime  + dt
         PRINT *, 'Current Time : ', currenttime 
      END DO
       
      !Calculate final energies  
      ! u solve
      CALL GetufromhG(xbc_len,xbc,hbc_fin,Gbc_fin,a5,
     . currenttime,b1a6,b1a7,
     . dx,n_GhstCells,ubc_fin)
                 
      END


c============
c  Finite Difference Solver Routines
c  Functions thet solve finite difference matrix Au = G
c  for u. A is tridiagonal (from finite different approximation to G equation), so we solve using Thomas algorithm
c============

c ====
c Function that builds the three diagonals of the matrix produced by the finite difference approximation to elliptic equation
c It then uses a matrix solve on the diagonals, and the value of G to solve for u (for gSGN equations).
c Inputs:  
c            xbc_len - total number of cells
c            xbc - cell centers of all cells
c            hbc - h at cell centers of all cells
c            Gbc - G at cell centers of all cells
c            a5 - beta1 parameter
c            t - current time
c            b1a6 - beta1 parameter
c            b1a7 - beta1 parameter : beta1 = b1a6*(x - a5*t) + b1a7
c            dx - cell width
c            n_GhstCells - number of ghost cells
c Input/Outputs:
c            ubc - In  : vector contains the known ghost cell values at the left and right boundaries
c                - Out : vector contains u at cell centers of all cells
c ====   
      SUBROUTINE GetufromhG(xbc_len,xbc,hbc,Gbc,a5,t,b1a6,b1a7,
     . dx,n_GhstCells,ubc)
      
      IMPLICIT NONE
      !Inputs  
      INTEGER n_GhstCells,xbc_len
      DOUBLE PRECISION xbc(xbc_len),hbc(xbc_len),Gbc(xbc_len)
      DOUBLE PRECISION dx,a5,t,b1a6,b1a7
      
      !Outputs  
      DOUBLE PRECISION ubc(xbc_len)
      
      !Local variables               
      DOUBLE PRECISION subdiag1(xbc_len),
     . diag(xbc_len),
     . supdiag1(xbc_len),
     . RHS(xbc_len)  
      DOUBLE PRECISION ht1,ht2,dhc,beta1
      INTEGER i

      ! Loop over all cell in interior and update diagonals and subdiagonals of A
      ! as well as right hand side B
      ! in Au = b matrix vector equation 
      !First cell in interior is n_GhstCells+1, as first n_GhstCells cells are ghost cells
      !Last cell in interior is xbc_len - n_GhstCells, as last n_GhstCells cells are ghost cells           
      DO i=n_GhstCells+1,xbc_len - n_GhstCells 

         beta1 = b1a6*(xbc(i) - a5*t) + b1a7
         ! calculate gradient dhc
         ! calculate common h terms in the diagonals        
         dhc = (hbc(i+1) - hbc(i-1)) /(2d0*dx)
         ht1 = (beta1/2d0)*(hbc(i)**3/(dx*dx))
         ht2 = (3d0/2d0)*(beta1)*
     .      (hbc(i)**2*dhc / (2*dx))
         
         ! write to diagonal arrays
         subdiag1(i)  = -ht1 + ht2
         diag(i) = hbc(i) + (2d0*ht1)
         supdiag1(i)  = -ht1 - ht2 
         
         !write to right hand side vector         
         RHS(i) = Gbc(i)
      END DO
      
      !Impose boundary conditions
      !first and last n_GhstCells  x n_GhstCells submatrices in tridiagonal matrix
      !Should be identity with the RHS being given by the known values of u in the ghost cells
      DO i = 1, n_GhstCells
         ! left boundary
         subdiag1(i) = 0d0
         diag(i) = 1d0
         supdiag1(i) = 0d0
         RHS(i) = ubc(i)
         
         !right boundary
         subdiag1(xbc_len - n_GhstCells + i) = 0d0
         diag(xbc_len - n_GhstCells + i) = 1d0
         supdiag1(xbc_len - n_GhstCells + i)= 0d0
         RHS(xbc_len - n_GhstCells + i) = 
     .      ubc(xbc_len - n_GhstCells + i)
      END DO

      ! Solve Au = b equation
      CALL ThomasTriSolve(xbc_len,subdiag1,diag,supdiag1,RHS,ubc)
      
      END


c ====
c Subroutine that solves a tridiagonal matrix equation Ax = d using Thomas Algorithm (direct solver)
c where A is tridiagonal with subdiagonal a, diagonal b and superdiagonal c
c Inputs:  
c            n - length of diagonal array (b)
c            a - subdiagonal of matrix A
c            b - diagonal of matrix A
c            c - superdiagonal of matrix A
c            d - right handside vector
c Outputs:
c            x - approximation to u at cell centers of all cells
c ====
      SUBROUTINE ThomasTriSolve(n,a,b,c,d,x)
      IMPLICIT NONE
      !Inputs  
      INTEGER n
      DOUBLE PRECISION a(n), c(n), b(n), d(n)
      
      !Outputs  
      DOUBLE PRECISION x(n)
      
      ! Local variables
      INTEGER i
      DOUBLE PRECISION q
      
      !Perform Elimination Step
      DO i = 2,n
         q = a(i)/b(i - 1)
         b(i) = b(i) - c(i - 1)*q
         d(i) = d(i) - d(i - 1)*q
      END DO
      
      !Perform Back Substitution Step
      q = d(n)/b(n)
      x(n) = q
      DO i = n - 1,1,-1
         q = (d(i) - c(i)*q)/b(i)
         x(i) = q
      END DO
      
      END


c============
c  Evolution Routines
c  Functions that evolve h^n, G^n  -> h^{n+1}, G^{n+1}
c  Using u solve, Kurganov method and SSP RK2 time stepping
c============


c ====
c Subroutine that combines other subroutines to evolve from h^n, G^n to h^n+1 and G^n+1
c
c Inputs:   a0 - background water depth on which travelling wave propagates
c           a1 - amplitude of travelling wave in h
c           a2 - translation speed of wave (wave is function of x - a2*t)
c           a3 - width of travelling wave
c           a4 - amplitude of travelling wave in u
c           a5 - beta parameter : speed (betas are a function of x - a5*t)
c           b1a6 - beta parameter :
c           b1a7 - beta parameter : beta1 = b1a6*(x - a5*t) + b1a7
c           b2a6 - beta parameter :
c           b2a7 - beta parameter : beta2 = b2a6*(x - a5*t) + b2a7 
c           t - current time
c           xbc - array of cell centers
c           xbc_len - total number of cells in discretisation
c           n_GhstCells - number of ghost cells at each boundary
c           ga  - acceleration due to gravity
c           dx - cell width
c           dt - time step size
c Input/Outputs:
c                  hbc - In  : array of h^n at all cells
c                        Out : array of h^n+1 at all cells
c                  Gbc - In  : array of G^n at all cells
c                        Out : array of G^n+1 at all cells
c                  ubc - In  : vector contains the known ghost cell values for u at the left and right boundaries
c                              middle is rubbish old u values
c                      - Out : vector contains the known ghost cell values for u at the left and right boundaries. 
c                              middle is rubbish old u values
c ====
      SUBROUTINE EvolveStepWrap(a0,a1,a2,a3,a4,a5,
     . b1a6,b1a7,b2a6,b2a7,t,xbc,
     . xbc_len,n_GhstCells,ga,
     . dx,dt,
     . hbc,Gbc,ubc)
      
      !Inputs    
      INTEGER n_GhstCells,xbc_len
      DOUBLE PRECISION xbc(xbc_len)
      DOUBLE PRECISION a0,a1,a2,a3,a4,a5,
     . b1a6,b1a7,b2a6,b2a7,t,ga,dx,dt
      
      !Outputs
      DOUBLE PRECISION hbc(xbc_len),Gbc(xbc_len),ubc(xbc_len)
     
      !Local variables
      DOUBLE PRECISION hpbc(xbc_len),Gpbc(xbc_len),
     . hppbc(xbc_len),Gppbc(xbc_len)
      INTEGER i
         
      !Solve for u
      CALL GetufromhG(xbc_len,xbc,hbc,Gbc,a5,t,b1a6,b1a7,
     . dx,n_GhstCells,ubc)

      !Update cell averages to h', G'
      CALL SingleEulerStepForced(xbc,t,a0,a1,a2,a3,a4,a5,
     . b1a6,b1a7,b2a6,b2a7,xbc_len,hbc,Gbc,ubc,ga,
     . n_GhstCells,dt,dx,hpbc,Gpbc)
     
      !update ghost cells, since EulerStep function, assumes they are constant
      CALL ForcedGaussian(xbc(1:n_GhstCells),n_GhstCells,
     .   t + dt,a0,a1,a2,a3,a4,a5,b1a6,b1a7,hpbc(1:n_GhstCells),
     .   ubc(1:n_GhstCells),Gpbc(1:n_GhstCells))
     
      CALL ForcedGaussian(xbc(xbc_len - n_GhstCells + 1:xbc_len),
     .   n_GhstCells, t + dt,a0,a1,a2,a3,a4,a5,b1a6,b1a7,
     .   hpbc(xbc_len - n_GhstCells + 1:xbc_len),
     .   ubc(xbc_len - n_GhstCells + 1:xbc_len),
     .   Gpbc(xbc_len - n_GhstCells + 1:xbc_len))
     
      !Solve for u - boundary conditions don't change for this analytic solution
      CALL GetufromhG(xbc_len,xbc,hpbc,Gpbc,a5,t+dt,b1a6,b1a7,
     . dx,n_GhstCells,ubc)
      
      !Update cell averages to h'', G''
      CALL SingleEulerStepForced(xbc,t+dt,a0,a1,a2,a3,a4,a5,
     . b1a6,b1a7,b2a6,b2a7,xbc_len,hpbc,Gpbc,ubc,ga,
     . n_GhstCells,dt,dx,hppbc,Gppbc)
      
      !Use RK timestepping to convexly combine (h^n,G^n), (h',G') and (h'',G'')
      !to producesecond order approximation to approximation to h^{n+1}, G^{n+1}
      !since boundary conditions are constant, the average will be the initial value
      DO i= 1,xbc_len
         hbc(i) = ( hbc(i) + hppbc(i))/2d0
         Gbc(i) = ( Gbc(i) + Gppbc(i))/2d0
      END DO


      !update ghost cells, since EulerStep function, assumes they are constant
      CALL ForcedGaussian(xbc(1:n_GhstCells),n_GhstCells,
     .   t + dt,a0,a1,a2,a3,a4,a5,b1a6,b1a7,hbc(1:n_GhstCells),
     .   ubc(1:n_GhstCells),Gbc(1:n_GhstCells))
     
      CALL ForcedGaussian(xbc(xbc_len - n_GhstCells + 1:xbc_len),
     .   n_GhstCells, t + dt,a0,a1,a2,a3,a4,a5,b1a6,b1a7,
     .   hbc(xbc_len - n_GhstCells + 1:xbc_len),
     .   ubc(xbc_len - n_GhstCells + 1:xbc_len),
     .   Gbc(xbc_len - n_GhstCells + 1:xbc_len))
     
      END
    

c ====
c Subroutine that performs an Euler time step for gSGN equations
c Inputs: 
c           xbc - array of cell centers
c           t - current time
c           a0 - background water depth on which travelling wave propagates
c           a1 - amplitude of travelling wave in h
c           a2 - translation speed of wave (wave is function of x - a2*t)
c           a3 - width of travelling wave
c           a4 - amplitude of travelling wave in u
c           a5 - beta parameter : speed (betas are a function of x - a5*t)
c           b1a6 - beta parameter :
c           b1a7 - beta parameter : beta1 = b1a6*(x - a5*t) + b1a7
c           b2a6 - beta parameter :
c           b2a7 - beta parameter : beta2 = b2a6*(x - a5*t) + b2a7 
c           xbc_len - total number of cells in discretisation
c           hbc - h values at cell centers (equivalent to cell averages, for second order)
c           Gbc - G values at cell centers (equivalent to cell averages, for second order)
c           ubc - u values at cell centers (equivalent to cell averages, for second order)
c           ga  - acceleration due to gravity
c           n_GhstCells - number of ghost cells at each boundary
c           dt - time step size
c           dx - cell width
c Outputs:
c            hpbc - h values at cell centers after Euler step (equivalent to cell averages, for second order)
c            Gpbc - G values at cell centers after Euler step (equivalent to cell averages, for second order)
c ====     
      SUBROUTINE SingleEulerStepForced(xbc,t,a0,a1,a2,a3,a4,a5,
     . b1a6,b1a7,b2a6,b2a7,xbc_len,hbc,Gbc,ubc,ga,
     . n_GhstCells,dt,dx,hpbc,Gpbc)

      !Inputs
      INTEGER n_GhstCells,xbc_len
      DOUBLE PRECISION xbc(xbc_len)
      DOUBLE PRECISION a0,a1,a2,a3,a4,a5,
     . b1a6,b1a7,b2a6,b2a7,t,ga,dx,dt
      DOUBLE PRECISION hbc(xbc_len),Gbc(xbc_len),ubc(xbc_len)
      
      !Outputs
      DOUBLE PRECISION hpbc(xbc_len),Gpbc(xbc_len)
     
      !Local variables
      DOUBLE PRECISION cdhi,cdGi,fih,fiG,foh,foG,beta1,beta2,
     . dfluxhdx,dfluxGdx,dhdt,dGdt
      INTEGER i
            
      !We fill in boundary conditions for hpbc and Gpbc assuming
      !boundary conditions are constant - which is true for analytic solution
      DO i = 1,n_GhstCells
      
         !left boundary
         hpbc(i) = hbc(i)
         Gpbc(i) = Gbc(i)
         
         !right boundary
         hpbc(xbc_len - n_GhstCells + i) = 
     .      hbc(xbc_len - n_GhstCells + i)
         Gpbc(xbc_len - n_GhstCells + i) = 
     .      Gbc(xbc_len - n_GhstCells + i)
      END DO
      
      !Update interior, first calculating flux across left boundary of interior
      !Then loop over cells to get flux in/out and thus updating cell average values
      
      !Left boundary (between last left ghost cell xbc(n_GhstCells) and first interior cell xbc(n_GhstCells) ) 
      i = n_GhstCells
      
      !Reconstructed Gradient of h,G,u across cell i (using slope limiting)
      CALL ReconLinLimGrad(hbc(i-1),hbc(i+1),cdhi)
      CALL ReconLinLimGrad(Gbc(i-1),Gbc(i+1),cdGi)
      CALL ReconLinLimGrad(ubc(i-1),ubc(i+1),cdui)
      
      !calculate beta values at flux edges and current time       
      beta1 = b1a6*(xbc(i) + 0.5*dx - a5*t) + b1a7
      beta2 = b2a6*(xbc(i) + 0.5*dx - a5*t) + b2a7
      
      !Calculate flux across boundary x_{i +1/2}
      !Routine also updates cdhi,cdGi to be gradient across cell i + 1 (next cell in loop)
      CALL Fluxxiph(xbc_len,hbc,Gbc,ubc,ga,beta1,beta2,dx,i,
     . cdhi,cdGi,cdui,foh,foG)
     
      !flux out becomes flux in on next cell
      fih = foh
      fiG = foG
      
      !loop over interior cells (do not update ghost cells)
      DO i = n_GhstCells + 1, xbc_len - n_GhstCells
      
         !calculate beta values at flux edges and current time       
         beta1 = b1a6*(xbc(i) + 0.5*dx - a5*t) + b1a7
         beta2 = b2a6*(xbc(i) + 0.5*dx - a5*t) + b2a7
         
         !Calculate flux across boundary x_{i +1/2}
         CALL Fluxxiph(xbc_len,hbc,Gbc,ubc,ga,beta1,beta2,dx,i,
     .                 cdhi,cdGi,cdui,foh,foG)
     
         !Calculate source term
         CALL Forcedht(xbc(i),t,a1,a2,a3,dhdt)
         CALL Forcedfluxhx(xbc(i),t,a0,a1,a2,a3,a4,dfluxhdx) 
         
         CALL ForcedGt(xbc(i),t,a0,a1,a2,a3,a4,a5,b1a6,b1a7,dGdt) 
         CALL ForcedfluxGx(xbc(i),t,ga,a0,a1,a2,a3,a4,a5,
     . b1a6,b1a7,b2a6,b2a7,dfluxGdx) 
         
         !Evolve system through time    
         hpbc(i) = hbc(i) + dt*dhdt 
     .    + dt*( dfluxhdx  - (foh - fih)/dx ) 
     
         Gpbc(i) =  Gbc(i)+ dt*dGdt 
     .    + dt*( dfluxGdx  - (foG - fiG)/dx ) 

         !flux out becomes flux in on next cell
         fih = foh
         fiG = foG
         
      END DO
           
      END
 
 
c ====
c Subroutine that given arrays, and i calculates flux across x_{i+1/2} for gSGN equations  
c Inputs:  
c            xbc_len - total number of cells in discretisation
c            hbc - h values at cell centers (equivalent to cell averages, for second order)
c            Gbc - G values at cell centers (equivalent to cell averages, for second order)
c            ubc - u values at cell centers (equivalent to cell averages, for second order)
c            ga  - acceleration due to gravity
c            beta1 - beta1 value of gSGN equations
c            beta2 - beta2 value of gSGN equations
c            dx - cell width           
c            i - cell index
c Input/Outputs:
c            cdhi - In  : reconstructed slope for h across cell i
c                   Out : reconstructed slope for h across cell i +1
c            cdGi - In  : reconstructed slope for G across cell i
c                   Out : reconstructed slope for G across cell i +1
c            cdui - In  : reconstructed slope for u across cell i
c                   Out : reconstructed slope for u across cell i +1
c Outputs:
c            foh - flux across x_{i+1/2} for h conservation equation
c            foG - flux across x_{i+1/2} for G conservation equation        
c ====
      SUBROUTINE Fluxxiph(xbc_len,hbc,Gbc,ubc,ga,beta1,beta2,dx,i,
     . cdhi,cdGi,cdui,foh,foG)
     
       !Inputs  
       INTEGER i,xbc_len
       DOUBLE PRECISION hbc(xbc_len),Gbc(xbc_len),ubc(xbc_len)
       DOUBLE PRECISION ga,beta1,beta2,dx

       !Inputs/Outputs
       DOUBLE PRECISION cdhi,cdGi,cdui
       
       !Outputs
       DOUBLE PRECISION foh,foG
     
       !Local variables
       DOUBLE PRECISION cdGip1,felG,felh,ferG,ferh,sr,sl,isrmsl,
     . hir,Gir,uir,duir,hip1l,Gip1l,uip1l,duip1l,
     . dhir,ddhir,dhip1l,ddhip1l, alpha
     
      !Use slopes * dx   in cell i to reconstruct at x_{i+1/2} from left    
      hir = hbc(i) + cdhi/2d0
      Gir = Gbc(i) + cdGi/2d0
      uir = ubc(i) + cdui/2d0
      
      !Approximate derivatives at x_{i+1/2}
      CALL ReconDq(ubc(i),ubc(i+1),dx,duir)
      CALL ReconDq(hbc(i),hbc(i+1),dx,dhir)
      CALL ReconDDq(hbc(i-1),hbc(i),hbc(i+1),hbc(i+2),dx,ddhir)     

      !Reconstruct slope * dx   on cell i+1
      CALL ReconLinLimGrad(hbc(i),hbc(i+2),cdhip1)
      CALL ReconLinLimGrad(Gbc(i),Gbc(i+2),cdGip1)
      CALL ReconLinLimGrad(ubc(i),ubc(i+2),cduip1)
      
      !Use slopes * dx   in cell i+1 to reconstruct at x_{i+1/2} from right   
      hip1l = hbc(i + 1) - cdhip1/2d0
      Gip1l = Gbc(i + 1) - cdGip1/2d0
      uip1l = ubc(i + 1) - cduip1/2d0
      
      !Approximate derivatives at x_{i+1/2}
      CALL ReconDq(ubc(i),ubc(i+1),dx,duip1l)
      CALL ReconDq(hbc(i),hbc(i+1),dx,dhip1l)
      CALL ReconDDq(hbc(i-1),hbc(i),hbc(i+1),hbc(i+2),dx,ddhip1l)    
     
      ! Wave speed bounds
      ! alpha = max(1,beta2/ beta) )
      ! only works if beta1 != 0
      ! if beta1 == 0, then must have beta2 =0 - enforced at user level
      IF  (beta1 < 10d0**(-10)) THEN
         alpha = 1
      ELSE
         alpha = MAX(1d0,beta2 / beta1)
      END IF
      
      sl  = MIN(0d0, uir - DSQRT(alpha*ga*hir) ,
     . uip1l - DSQRT(alpha*ga*hip1l)  )
      sr  = MAX(0d0, uir + DSQRT(alpha*ga*hir) ,
     . uip1l + DSQRT(alpha*ga*hip1l)  )
      
      !Flux evalutated at left
      felh = uir*hir
      felG = uir*Gir + ga*(hir**2)/2d0 
     .      - beta1*hir**3*duir**2
     .      - beta2/4d0*ga*(hir**2)*(2*hir*ddhir + (dhir**2))
     
      !Flux evalutated at right   
      ferh = uip1l*hip1l
      ferG = uip1l*Gip1l + ga*(hip1l**2)/2d0 
     .      - beta1*hip1l**3*duip1l**2
     .      - beta2/4d0*ga*(hip1l**2)*(2*hip1l*ddhip1l + (dhip1l**2))
     
      !Wave speed weighting 1/ (sr - sl) - return 0 if sr = sl
      IF (sr == sl) THEN
         isrmsl = 0.0
      ELSE
         isrmsl = 1.0 / (sr - sl)
      END IF
      
      !calculate flux from cell i to cell i + 1 (Kurganov)
      foh = isrmsl*(sr*felh - sl*ferh + sl*sr*(hip1l - hir))
      foG = isrmsl*(sr*felG - sl*ferG + sl*sr*(Gip1l - Gir))
      
      !return gradient of h,G across cell i + 1 (for next iteration)
      cdGi = cdGip1
      cdhi = cdhip1   
      cdui = cduip1 
      END

c ====
c Subroutine that calculates reconstruction slope without limiting central difference 
c Inputs:  
c         qim1 - q(i-1) - average of quantity in i-1
c         qi - q(i) - average of quantity in i
c         qip1 - q(i+1) - average of quantity in i+1
c Outputs:
c         cdq - calculated slope * dx    
c ====
      SUBROUTINE ReconLinLimGrad(qim1,qip1,cdq)
      !Inputs
      DOUBLE PRECISION qim1,qip1
      !Outputs
      DOUBLE PRECISION cdq

      cdq  = 0.5*(qip1 - qim1)

      END

c ====
c Subroutine that calculates derivative dq/dx at x_{i+1/2} 
c Inputs:  
c         qi - q(i) - average of quantity in i
c         qip1 - q(i+1) - average of quantity in i+1
c         dx - cell width
c Outputs:
c         cdq - calculated derivative dq/dx at x_{i+1/2}     
c ====
      SUBROUTINE ReconDq(qi,qip1,dx,cdq)
      !Input
      DOUBLE PRECISION qi,qip1,dx
      !Output
      DOUBLE PRECISION cdq
      cdq = (qip1 - qi) /dx
      END

c ====
c Subroutine that calculates derivative dq/dx at x_{i+1/2} 
c Inputs: 
c         qim1 - q(i-1) - average of quantity in i -1
c         qi - q(i) - average of quantity in i
c         qip1 - q(i+1) - average of quantity in i+1
c         dx - cell width
c Outputs:
c         cddq - calculated derivative d^2q/dx^2 at x_{i+1/2}     
c ====
      SUBROUTINE ReconDDq(qim1,qi,qip1,qip2,dx,cddq)
      !Input
      DOUBLE PRECISION qim1,qi,qip1,qip2,dx
      !Output
      DOUBLE PRECISION cddq
      cddq = (qip2  - qip1 - qi + qim1 )/(2*dx**2)   
      END

c ====
c Subroutine that calculates minmod(a,b,c)
c Inputs: 
c         a,b,c - double precision numbers
c Outputs:
c         d - minmod(a,b,c)  
c ====
      SUBROUTINE minmod(a,b,c,d)
      IMPLICIT NONE
      DOUBLE PRECISION a,b,c, d
      
      IF ((a .GT. 0d0) .AND. (b .GT. 0d0) .AND. (c .GT. 0d0)) THEN
         d = MIN(a,b,c)
      ELSE IF ((a .LT. 0d0) .AND. (b .LT. 0d0) .AND. (c .LT. 0d0)) THEN
         d = MAX(a,b,c)
      ELSE
         d = 0d0
      END IF
      
      END

c============
c  Source Term Routines - Calculates analytic values of dh/dt, dG/dt and df(u,G,h)/dx [for both conserved quantities h and G]
c============

c ====
c Subroutine that calculates analytic value of dh/dt for Gaussian forced solution
c Inputs: 
c         x - location to calculate flux dh/dt at
c         t - time to calculate flux dh/dt at
c         a1 - amplitude of wave in h
c         a2 - speed of translation (function of x - a2*t)
c         a3 - width of wave
c Outputs:
c         dhdt - analytic value of dh/dt 
c ====
      SUBROUTINE Forcedht(x,t,a1,a2,a3,dhdt) 
      !Input
      DOUBLE PRECISION x,t,a1,a2,a3
      !Output
      DOUBLE PRECISION dhdt
      !Local variables
      DOUBLE PRECISION :: EXPPHI1,PHI
      
      PHI  = x - a2*t
      EXPPHI1 = DEXP(-PHI**2 / (2*a3))
      dhdt = EXPPHI1*PHI*a1*a2/a3
      
      END 

c ====
c Subroutine that calculates analytic value of df(u,G,h)/dt for h in Gaussian forced solution
c Inputs: 
c         x - location to calculate flux dh/dt at
c         t - time to calculate flux dh/dt at
c         a1 - amplitude of wave in h
c         a2 - speed of translation (function of x - a2*t)
c         a3 - width of wave
c         a4 - amplitude of wave in u
c Outputs:
c         dfluxhdx - analytic value of df(u,G,h)/dt for h
c ====      
      SUBROUTINE  Forcedfluxhx(x,t,a0,a1,a2,a3,a4,dfluxhdx) 
      !Input
      DOUBLE PRECISION x,t,a1,a2,a4
      !Output
      DOUBLE PRECISION dfluxhdx
      !Local variables
      DOUBLE PRECISION EXPPHI1,PHI,hi,ui,dhdx,dudx
      
      PHI  = x - a2*t
      EXPPHI1 = DEXP(-PHI**2 / (2*a3))
      hi  = a0 + a1*EXPPHI1
      ui  = a4*EXPPHI1
      dhdx = -EXPPHI1*PHI*a1/a3
      dudx = -EXPPHI1*PHI*a4/a3
      
      dfluxhdx = ui*dhdx + hi*dudx

      END

c ====
c Subroutine that calculates analytic value of dG/dt for Gaussian forced solution
c Inputs: 
c         x - location to calculate flux dh/dt at
c         t - time to calculate flux dh/dt at
c         a0 - background water depth on which travelling wave propagates
c         a1 - amplitude of travelling wave in h
c         a2 - translation speed of wave (wave is function of x - a2*t)
c         a3 - width of travelling wave
c         a4 - amplitude of travelling wave in u
c         a5 - beta parameter : speed (betas are a function of x - a5*t)
c         b1a6 - beta parameter :
c         b1a7 - beta parameter : beta1 = b1a6*(x - a5*t) + b1a7
c Outputs:
c         dGdt - analytic value of dG/dt 
c ====
      SUBROUTINE ForcedGt(x,t,a0,a1,a2,a3,a4,a5,b1a6,b1a7,dGdt) 
      !Input
      DOUBLE PRECISION x,t,a0,a1,a2,a3,a4,a5,b1a6,b1a7
      !Output
      DOUBLE PRECISION dGdt
      !Local variables
      DOUBLE PRECISION :: EXPPHI1,PHI,hi,ui,beta1
      
      PHI  = x - a2*t
      EXPPHI1 = DEXP(-PHI**2 / (2*a3))
      hi  = a0 + a1*EXPPHI1
      ui  = a4*EXPPHI1
      beta1 = b1a6*(x - a5*t) + b1a7
      dGdt = -3*EXPPHI1**2*PHI**3*a1**2*a2*beta1*hi*ui/a3**3 
     . + 3*EXPPHI1*PHI**2*a1*a5*b1a6*hi**2*ui/(2*a3**2) 
     . + EXPPHI1*PHI*a1*a2*ui/a3 
     . - 9*PHI**3*a1*a2*beta1*hi**2*ui*DEXP(-PHI**2/(2*a3))/(2*a3**3) 
     . - PHI**3*a2*a4*beta1*hi**3*DEXP(-PHI**2/(2*a3))/(2*a3**3) 
     . + PHI**2*a4*a5*b1a6*hi**3*DEXP(-PHI**2/(2*a3))/(2*a3**2) 
     . + 9*PHI*a1*a2*beta1*hi**2*ui*DEXP(-PHI**2/(2*a3))/(2*a3**2) 
     . + PHI*a2*hi*ui/a3 
     . + 3*PHI*a2*a4*beta1*hi**3*DEXP(-PHI**2/(2*a3))/(2*a3**2) 
     . - a4*a5*b1a6*hi**3*DEXP(-PHI**2/(2*a3))/(2*a3)
      
      END 

c ====
c Subroutine that calculates analytic value of df(u,h,G)/dx for G for Gaussian forced solution
c Inputs: 
c         x - location to calculate flux dh/dt at
c         t - time to calculate flux dh/dt at
c         a0 - background water depth on which travelling wave propagates
c         a1 - amplitude of travelling wave in h
c         a2 - translation speed of wave (wave is function of x - a2*t)
c         a3 - width of travelling wave
c         a4 - amplitude of travelling wave in u
c         a5 - beta parameter : speed (betas are a function of x - a5*t)
c         b1a6 - beta parameter :
c         b1a7 - beta parameter : beta1 = b1a6*(x - a5*t) + b1a7
c         b2a6 - beta parameter :
c         b2a7 - beta parameter : beta2 = b2a6*(x - a5*t) + b2a7
c Outputs:
c         dfluxGdx - analytic value of df(u,h,G)/dx
c ====
      SUBROUTINE ForcedfluxGx(x,t,ga,a0,a1,a2,a3,a4,a5,
     . b1a6,b1a7,b2a6,b2a7,dfluxGdx) 
      !Input
      DOUBLE PRECISION x,t,ga,a0,a1,a2,a3,a4,a5,b1a6,b1a7,b2a6,b2a7
      !Output
      DOUBLE PRECISION dfluxGdx
      !Local variables
      DOUBLE PRECISION :: EXPPHI1,PHI,hi,ui,beta1,beta2,dudx,d2hdx2,
     .d2udx2, d3hdx3, d3udx3 ,dhdx 
      
      PHI  = x - a2*t
      EXPPHI1 = DEXP(-PHI**2 / (2*a3))
      hi  = a0 + a1*EXPPHI1
      ui  = a4*EXPPHI1
      
      dhdx = -EXPPHI1*PHI*a1/a3
      dudx = -EXPPHI1*PHI*a4/a3
      d2udx2 = -a4*(-PHI**2 + a3)*DEXP(-PHI**2/(2*a3))/a3**2
      d2hdx2 = -a1*(-PHI**2 + a3)*DEXP(-PHI**2/(2*a3))/a3**2
      d3hdx3 = PHI*a1*(-PHI**2 + 3*a3)*DEXP(-PHI**2/(2*a3))/a3**3
      d3udx3 = PHI*a4*(-PHI**2 + 3*a3)*DEXP(-PHI**2/(2*a3))/a3**3
      
      beta1 = b1a6*(x - a5*t) + b1a7
      beta2 = b2a6*(x - a5*t) + b2a7
      
      dfluxGdx = -b1a6*d2udx2*hi**3*ui/2 
     . - 3*b1a6*dhdx*dudx*hi**2*ui/2 
     . - b1a6*dudx**2*hi**3 
     . - b2a6*d2hdx2*ga*hi**3/2 
     . - b2a6*dhdx**2*ga*hi**2/4 
     . - 3*beta1*d2hdx2*dudx*hi**2*ui/2 
     . - 3*beta1*d2udx2*dhdx*hi**2*ui 
     . - 5*beta1*d2udx2*dudx*hi**3/2 
     . - beta1*d3udx3*hi**3*ui/2 
     . - 3*beta1*dhdx**2*dudx*hi*ui 
     . - 9*beta1*dhdx*dudx**2*hi**2/2 
     . - 2*beta2*d2hdx2*dhdx*ga*hi**2 
     . - beta2*d3hdx3*ga*hi**3/2 
     . - beta2*dhdx**3*ga*hi/2 
     . + dhdx*ga*hi 
     . + dhdx*ui**2 
     . + 2*dudx*hi*ui

      
      END 
