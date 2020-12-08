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
c Inputs:  
c           tstart - start time
c           tend - time to produce numerical solution to
c           ga  - acceleration due to gravity
c           beta1 - beta1 value of gSGN equations
c           beta2 - beta2 value of gSGN equations
c           theta - reconstruction parameter theta for generalised minmod limiter
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
      SUBROUTINE NumericalSolveTSPrint(tstart,tend,
     . ga,beta1,beta2,theta,dx,dt,n_GhstCells,xbc_len,
     . xbc,hbc_init,Gbc_init,ubc_init,
     . tlist,tlist_len,ExpWdir,expwdirlen,
     . currenttime,hbc_fin,Gbc_fin,ubc_fin,
     . Energ_Init, Energ_Fin)
     
     
      IMPLICIT NONE
      !Inputs  
      INTEGER n_GhstCells,xbc_len,expwdirlen,tlist_len
      CHARACTER(LEN=expwdirlen) ExpWdir
      DOUBLE PRECISION tstart,tend,ga,beta1,beta2,
     . theta,dx,dt,currenttime
      DOUBLE PRECISION xbc(xbc_len),hbc_init(xbc_len),
     . Gbc_init(xbc_len), ubc_init(xbc_len),
     . tlist(tlist_len)
     
      !Outputs  
      DOUBLE PRECISION Energ_Init(4), Energ_Fin(4),hbc_fin(xbc_len),
     . Gbc_fin(xbc_len), ubc_fin(xbc_len)
      
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
      
      !calculate initial Energies  
      CALL TotalEnergy(xbc_len,n_GhstCells,dx,ga,beta1,beta2,
     . hbc_init,ubc_init,Gbc_init,
     . Energ_Init)

      !Evolve system through time
      DO WHILE (currenttime  .LT. tend )   
c      DO WHILE (currenttime  .LT. 10.0**(-14) )         
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
         CALL EvolveStepWrap(xbc_len,n_GhstCells,ga,beta1,beta2,
     .    dx,dt,theta,hbc_fin,Gbc_fin,ubc_fin)
         !Update current time and print, so we can get some feedback about what its doing
         currenttime  = currenttime  + dt
         PRINT *, 'Current Time : ', currenttime 
      END DO
       
      !Calculate final energies  
      ! u solve
      CALL GetufromhG(xbc_len,hbc_fin,Gbc_fin,beta1,dx,n_GhstCells,
     .  ubc_fin)

      CALL TotalEnergy(xbc_len,n_GhstCells,dx,ga,beta1,beta2,
     . hbc_fin,ubc_fin,Gbc_fin,
     . Energ_Fin)

                 
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
c            hbc - h at cell centers of all cells
c            Gbc - G at cell centers of all cells
c            beta1 - beta1 parameter in gSGN
c            dx - cell width
c            n_GhstCells - number of ghost cells
c Input/Outputs:
c            ubc - In  : vector contains the known ghost cell values at the left and right boundaries
c                - Out : vector contains u at cell centers of all cells
c ====    
      SUBROUTINE GetufromhG(xbc_len,hbc,Gbc,beta1,dx,n_GhstCells,ubc)
      
      IMPLICIT NONE
      !Inputs  
      INTEGER n_GhstCells,xbc_len
      DOUBLE PRECISION hbc(xbc_len),Gbc(xbc_len)
      DOUBLE PRECISION dx,beta1
      
      !Outputs  
      DOUBLE PRECISION ubc(xbc_len)
      
      !Local variables               
      DOUBLE PRECISION subdiag1(xbc_len),
     . diag(xbc_len),
     . supdiag1(xbc_len),
     . RHS(xbc_len)  
      DOUBLE PRECISION ht1,ht2,dhc
      INTEGER i
      
      !print *, 'u Solve', beta1,dx,n_GhstCells 
                  
      ! Loop over all cell in interior and update diagonals and subdiagonals of A
      ! as well as right hand side B
      ! in Au = b matrix vector equation 
      !First cell in interior is n_GhstCells+1, as first n_GhstCells cells are ghost cells
      !Last cell in interior is xbc_len - n_GhstCells, as last n_GhstCells cells are ghost cells           
      DO i=n_GhstCells+1,xbc_len - n_GhstCells 

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
c Inputs:  
c           xbc_len - total number of cells in discretisation
c           n_GhstCells - number of ghost cells at each boundary
c           ga  - acceleration due to gravity
c           beta1 - beta1 value of gSGN equations
c           beta2 - beta2 value of gSGN equations
c           dx - cell width
c           dt - time step size
c           theta - reconstruction parameter theta for generalised minmod limiter
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

      SUBROUTINE EvolveStepWrap(xbc_len,n_GhstCells,ga,beta1,beta2,
     . dx,dt,theta,
     . hbc,Gbc,ubc)
      
      !Inputs    
      INTEGER n_GhstCells,xbc_len
      DOUBLE PRECISION ga,beta1,beta2,theta,dx,dt
      
      !Outputs
      DOUBLE PRECISION hbc(xbc_len),Gbc(xbc_len),ubc(xbc_len)
     
      !Local variables
      DOUBLE PRECISION hpbc(xbc_len),Gpbc(xbc_len),
     . hppbc(xbc_len),Gppbc(xbc_len)
      INTEGER i
      
         
      !Solve for u
      CALL GetufromhG(xbc_len,hbc,Gbc,beta1,dx,n_GhstCells,ubc)

      !Update cell averages to h', G'
      CALL SingleEulerStep(xbc_len,hbc,Gbc,ubc,ga,beta1,beta2
     . ,theta,n_GhstCells,dt,dx,hpbc,Gpbc)
     
      !Solve for u - boundary conditions don't change for this analytic solution
      CALL GetufromhG(xbc_len,hpbc,Gpbc,beta1,dx,n_GhstCells,ubc)
      
      !Update cell averages to h'', G''
      CALL SingleEulerStep(xbc_len,hpbc,Gpbc,ubc,ga,beta1,beta2
     . ,theta,n_GhstCells,dt,dx,hppbc,Gppbc)
      
      !Use RK timestepping to convexly combine (h^n,G^n), (h',G') and (h'',G'')
      !to producesecond order approximation to approximation to h^{n+1}, G^{n+1}
      !since boundary conditions are constant, the average will be the initial value
      DO i= 1,xbc_len
         hbc(i) = ( hbc(i) + hppbc(i))/2d0
         Gbc(i) = ( Gbc(i) + Gppbc(i))/2d0
      END DO

      END
    

c ====
c Subroutine that performs an Euler time step for gSGN equations
c Inputs:  
c            xbc_len - total number of cells in discretisation
c            hbc - h values at cell centers (equivalent to cell averages, for second order)
c            Gbc - G values at cell centers (equivalent to cell averages, for second order)
c            ubc - u values at cell centers (equivalent to cell averages, for second order)
c            ga  - acceleration due to gravity
c            beta1 - beta1 value of gSGN equations
c            beta2 - beta2 value of gSGN equations
c            theta - reconstruction parameter theta for generalised minmod limiter
c            n_GhstCells - number of ghost cells at each boundary
c            dt - time step size
c            dx - cell width
c Outputs:
c            hpbc - h values at cell centers after Euler step (equivalent to cell averages, for second order)
c            Gpbc - G values at cell centers after Euler step (equivalent to cell averages, for second order)
c ====
      SUBROUTINE SingleEulerStep(xbc_len,hbc,Gbc,ubc,ga,beta1,beta2,
     . theta,n_GhstCells,dt,dx,hpbc,Gpbc)

      !Inputs
      INTEGER n_GhstCells,xbc_len
      DOUBLE PRECISION dt,dx,ga,beta1,beta2,theta
      DOUBLE PRECISION hbc(xbc_len),Gbc(xbc_len),ubc(xbc_len)
      
      !Outputs
      DOUBLE PRECISION hpbc(xbc_len),Gpbc(xbc_len)
     
      !Local variables
      DOUBLE PRECISION cdhi,cdGi, fih,fiG,foh,foG
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
      CALL ReconLinLimGrad(hbc(i-1),hbc(i),hbc(i+1),theta,cdhi)
      CALL ReconLinLimGrad(Gbc(i-1),Gbc(i),Gbc(i+1),theta,cdGi)
      CALL ReconLinLimGrad(ubc(i-1),ubc(i),ubc(i+1),theta,cdui)
      
      !Calculate flux across boundary x_{i +1/2}
      !Routine also updates cdhi,cdGi to be gradient across cell i + 1 (next cell in loop)
      CALL Fluxxiph(xbc_len,hbc,Gbc,ubc,ga,beta1,beta2,theta,dx,i,
     . cdhi,cdGi,cdui,foh,foG)
     
      !flux out becomes flux in on next cell
      fih = foh
      fiG = foG
      
      !loop over interior cells (do not update ghost cells)
      DO i = n_GhstCells + 1, xbc_len - n_GhstCells
      
         !Calculate flux across boundary x_{i +1/2}
         CALL Fluxxiph(xbc_len,hbc,Gbc,ubc,ga,beta1,beta2,theta,dx,i,
     .                 cdhi,cdGi,cdui,foh,foG)
         
         !Evolve system through time    
         hpbc(i) = hbc(i) - dt*(foh - fih)/dx 
         Gpbc(i) = Gbc(i) - dt*(foG - fiG)/dx 

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
c            theta - reconstruction parameter theta for generalised minmod limiter
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
      SUBROUTINE Fluxxiph(xbc_len,hbc,Gbc,ubc,ga,beta1,beta2,theta,dx,i,
     . cdhi,cdGi,cdui,foh,foG)
     
       !Inputs  
       INTEGER i,xbc_len
       DOUBLE PRECISION hbc(xbc_len),Gbc(xbc_len),ubc(xbc_len)
       DOUBLE PRECISION ga,beta1,beta2,theta,dx

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
      CALL ReconLinLimGrad(hbc(i),hbc(i+1),hbc(i+2),theta,cdhip1)
      CALL ReconLinLimGrad(Gbc(i),Gbc(i+1),Gbc(i+2),theta,cdGip1)
      CALL ReconLinLimGrad(ubc(i),ubc(i+1),ubc(i+2),theta,cduip1)
      
      !Use slopes * dx   in cell i+1 to reconstruct at x_{i+1/2} from right   
      hip1l = hbc(i + 1) - cdhip1/2d0
      Gip1l = Gbc(i + 1) - cdGip1/2d0
      uip1l = ubc(i + 1) - cduip1/2d0
      
      !Approximate derivatives at x_{i+1/2}
      CALL ReconDq(ubc(i),ubc(i+1),dx,duip1l)
      CALL ReconDq(hbc(i),hbc(i+1),dx,dhip1l)
      CALL ReconDDq(hbc(i-1),hbc(i),hbc(i+1),hbc(i+2),dx,ddhip1l)    
     
      ! Wave speed bounds
      ! alpha = max(1,beta2/ (beta1) )
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
c Subroutine that calculates reconstruction slope with limiting (generalised minmmod limter)  
c Inputs:  
c         qim1 - q(i-1) - average of quantity in i-1
c         qi - q(i) - average of quantity in i
c         qip1 - q(i+1) - average of quantity in i+1
c         theta - reconstruction parameter for generalised minmod
c Outputs:
c         cdq - calculated slope * dx    
c ====
      SUBROUTINE ReconLinLimGrad(qim1,qi,qip1,theta,cdq)
      !Inputs
      DOUBLE PRECISION qim1,qi,qip1,theta
      !Outputs
      DOUBLE PRECISION cdq
      !Local Variable
      DOUBLE PRECISION fdq,mdq,bdq
      
      fdq = qip1 - qi
      mdq = 0.5*(qip1 - qim1)
      bdq = qi - qim1
      
      !Minmod limiter for slope
      CALL minmod(theta*fdq,mdq,theta*bdq,cdq)
      
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
c  Analysis Routines - Calculates Total Amount of Conserved Quantities
c  In particular h,G, uh and Energy function
c============
      
c =====
c Routine to generate interpolcation quartic polynomial coefficients using cell averages in cells
c j-2,j-1,j,j+1,j+2
c Inputs: 
c          qjm2 - j-2 cell average 
c          qjm1 - j-1 cell average 
c          qj - j cell average 
c          qjp1 - j+1 cell average 
c          qjp2  - j+2 cell average 
c          dx - cell width
c Outputs:
c         QuartCoeff - array of length 5 containing coefficients
c           that interpolate the function values.
c =====      
      SUBROUTINE QuarticInterp(qjm2,qjm1,qj,qjp1,qjp2,dx,QuartCoeff)
     
      !Input
      DOUBLE PRECISION qjm2,qjm1,qj,qjp1,qjp2,dx
      !Output
      DOUBLE PRECISION QuartCoeff(5)      
      QuartCoeff(1) = (6*qj - 4*qjm1 + qjm2 - 4*qjp1 + qjp2)/
     . (24*dx**4)
      QuartCoeff(2) = (2*qjm1 - qjm2 - 2*qjp1 + qjp2)/(12*dx**3)
      QuartCoeff(3) = (-22*qj + 12*qjm1 - qjm2 + 12*qjp1 - qjp2)/
     . (16*dx**2) 
      QuartCoeff(4) = (-34*qjm1 + 5*qjm2 + 34*qjp1 - 5*qjp2)/(48*dx)
      QuartCoeff(5) = 1067*qj/960 - 29*qjm1/480 + 3*qjm2/640 
     . - 29*qjp1/480 + 3*qjp2/640
      
      END


c =====
c Routine to use quartic polynomial coefficients to calculate value at x - x_i
c Inputs: 
c          QuarticCoeff - Array of coefficients
c          xmxj -  x- x_i location for evaluation
c Outputs:
c         qatxj - value of P(xmxj) where coeficients of polynomial P given by QuarticCoeff
c =====      
      SUBROUTINE QuarticCoeffEvalxj(QuarticCoeff,xmxj,qatxj)
      
      !Input
      DOUBLE PRECISION QuarticCoeff(5),xmxj
      
      !Output
      DOUBLE PRECISION qatxj
      
      qatxj = QuarticCoeff(1)*xmxj**4 + QuarticCoeff(2)*xmxj**3 +
     . QuarticCoeff(3)*xmxj**2 + QuarticCoeff(4)*xmxj +
     . QuarticCoeff(5)
     
      END

c =====
c Routine to use quartic polynomial coefficients to calculate value of derivative at x - x_i
c Inputs: 
c          QuarticCoeff - Array of coefficients
c          xmxj -  x- x_i location for evaluation
c Outputs:
c         dqatxj - value of dP/dx(xmxj) where coefficients of polynomial P given by QuarticCoeff
c =====         
      SUBROUTINE QuarticCoeffEvalGradxj(QuarticCoeff,xmxj,dqatxj)
      
      !Input
      DOUBLE PRECISION QuarticCoeff(5),xmxj
      
      !Output
      DOUBLE PRECISION dqatxj
      
      dqatxj = 4*QuarticCoeff(1)*xmxj**3 + 3*QuarticCoeff(2)*xmxj**2 +
     . 2*QuarticCoeff(3)*xmxj + QuarticCoeff(4)
     
      END

c =====
c Routine to calculate total amount of conserved quantities in cell
c Inputs: 
c            xbc_len - number of cells
c            h - h at cell centers (equivalent to cell average values because its second-order )
c            u - u at cell centers
c            G - G at cell centers
c            ga - acceleration due to gravity
c            beta1 - beta1 gSGN parameter
c            beta2 - beta2 gSGN parameter
c            j - cell index
c            dx - cell width
c Outputs:
c            CellEnergies - array of total amount of quantities in cell (h,G,uh,Energy)
c ===== 
      SUBROUTINE AllEnergiesIntegralCell(xbc_len,h,u,G,ga,beta1,beta2,j
     . ,dx,CellEnergies)
      
      !Input
      INTEGER j,xbc_len
      DOUBLE PRECISION h(xbc_len),u(xbc_len),G(xbc_len)
      DOUBLE PRECISION dx,ga,beta1,beta2
      
      !Output
      DOUBLE PRECISION CellEnergies(4)
      
      !Local variable
      INTEGER i
      DOUBLE PRECISION fGPe(4),sGPe(4),tGPe(4)
      DOUBLE PRECISION GPmxj,hGP,GGP,uGP,uxGP,hxGP
      DOUBLE PRECISION hCoeff(5), uCoeff(5), GCoeff(5)
      
      ! Coefficients of quatrics that interpolate h,u,G
      CALL QuarticInterp(h(j-2),h(j-1),h(j),h(j+1),h(j+2),dx,hCoeff)
      CALL QuarticInterp(u(j-2),u(j-1),u(j),u(j+1),u(j+2),dx,uCoeff)
      CALL QuarticInterp(G(j-2),G(j-1),G(j),G(j+1),G(j+2),dx,GCoeff)
      
      !Location of First Gauss Point
      GPmxj = -dx*DSQRT(3.0d0/5.0d0)/2
      
      !Evaluate required functions at gauss point (h,dh/dx,G,u,du/dx)
      CALL QuarticCoeffEvalxj(hCoeff,GPmxj,hGP)
      CALL QuarticCoeffEvalGradxj(hCoeff,GPmxj,hxGP)
      CALL QuarticCoeffEvalxj(GCoeff,GPmxj,GGP)
      CALL QuarticCoeffEvalxj(uCoeff,GPmxj,uGP)
      CALL QuarticCoeffEvalGradxj(uCoeff,GPmxj,uxGP)
      
      !Combine into array of values at first Gauss point
      fGPe(1) = hGP
      fGPe(2) = GGP
      fGPe(3) = hGP*uGP
      fGPe(4) = (hGP*uGP**2 + beta1/2d0*(uxGP**2)*(hGP**3)
     . + ga*hGP**2*(1d0 + beta2/2d0*hxGP**2 )   )/2d0

      !Location of Second Gauss Point
      GPmxj = 0.0 
      
      !Evaluate required functions at gauss point (h,dh/dx,G,u,du/dx)
      CALL QuarticCoeffEvalxj(hCoeff,GPmxj,hGP)
      CALL QuarticCoeffEvalGradxj(hCoeff,GPmxj,hxGP)
      CALL QuarticCoeffEvalxj(GCoeff,GPmxj,GGP)
      CALL QuarticCoeffEvalxj(uCoeff,GPmxj,uGP)
      CALL QuarticCoeffEvalGradxj(uCoeff,GPmxj,uxGP)
      
      !Combine into array of values at second Gauss point
      sGPe(1) = hGP
      sGPe(2) = GGP
      sGPe(3) = hGP*uGP
      sGPe(4) = (hGP*uGP**2 +  beta1/2d0*(uxGP**2)*(hGP**3)
     . + ga*hGP**2*(1d0 + beta2/2d0*hxGP**2 )   )/2d0
      
      !Location of Third Gauss Point
      GPmxj = dx*DSQRT(3.0d0/5.0d0)/2
      
      !Evaluate required functions at gauss point (h,dh/dx,G,u,du/dx)
      CALL QuarticCoeffEvalxj(hCoeff,GPmxj,hGP)
      CALL QuarticCoeffEvalGradxj(hCoeff,GPmxj,hxGP)
      CALL QuarticCoeffEvalxj(GCoeff,GPmxj,GGP)
      CALL QuarticCoeffEvalxj(uCoeff,GPmxj,uGP)
      CALL QuarticCoeffEvalGradxj(uCoeff,GPmxj,uxGP)
      
      !Combine into array of values at third Gauss point
      tGPe(1) = hGP
      tGPe(2) = GGP
      tGPe(3) = hGP*uGP
      tGPe(4) = (hGP*uGP**2 + beta1/2d0*(uxGP**2)*(hGP**3)
     . + ga*hGP**2*(1d0 + beta2/2d0*hxGP**2 )   )/2d0
      
      !Combine function evaluations to approximate integral using Gaussian quadrature.
      DO i = 1,4
         CellEnergies(i) = (dx /2d0)*( (5.0/9.0)*fgpe(i) 
     .                     + (8.0/9.0)*sgpe(i) + (5.0/9.0)*tgpe(i))
      END DO
      
      END
      

c =====
c Routine to calculate total conserved quantities in whole domain, 
c by summing the total amount in each cell
c
c Inputs: 
c            xbc_len - total number of cells
c            n_GhstCells - number of ghost cells 
c            dx - cell width
c            ga - acceleration due to gravity
c            beta1 - beta1 gSGN parameter
c            beta2 - beta2 gSGN parameter
c            hbc - h at cell centers (equivalent to cell average values because its second-order )
c            ubc - u at cell centers
c            Gbc - G at cell centers
c Outputs:
c            TotEnergVals - array of total amount of quantities in cell (h,G,uh,Energy)
c ===== 
      SUBROUTINE TotalEnergy(xbc_len,n_GhstCells,dx,ga,beta1,beta2,
     . hbc,ubc,Gbc,
     . TotEnergVals)
      
      !Inputs
      INTEGER xbc_len,n_GhstCells
      DOUBLE PRECISION hbc(xbc_len),ubc(xbc_len),Gbc(xbc_len)
      DOUBLE PRECISION dx,ga,beta1,beta2
      
      !Outputs
      DOUBLE PRECISION TotEnergVals(4)
      
      !Local variables
      DOUBLE PRECISION CellEnergVals(4)
      INTEGER i,j
            
      !running totals for energy values, start at 0
      DO i = 1,4
         TotEnergVals(i) = 0.0
      END DO
      
      !just loop over interior of hbc,Gbc, ubc which have interior values + ghost cell values
      DO j= n_GhstCells + 1, xbc_len - n_GhstCells
         CALL  AllEnergiesIntegralCell(xbc_len,hbc,ubc,Gbc,ga,
     .      beta1,beta2,j,dx,CellEnergVals)
     
         !add cell energy value to running total
         DO i = 1,4
            TotEnergVals(i) = TotEnergVals(i) + CellEnergVals(i)
         END DO
      END DO
      
      END

      


