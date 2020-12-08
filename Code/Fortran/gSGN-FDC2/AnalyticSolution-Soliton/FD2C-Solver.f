c ====================================================================================
c Collection of Fortran subroutines that use naive centered finite difference based numerical method to solve generalised Serre - Green -Naghdi equations (gSGN)
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
c           dx - cell width
c           dt - time step size
c           n_GhstCells - number of ghost cells at each boundary
c           xbc_len - total number of cells in discretisation
c           xbc - array of cell centers from discretisation
c           hbc_init - values of h at the cell centers (xbc) at t = tstart 
c           ubc_init - values of u at the cell centers (xbc) at t = tstart 
c           phbc_init - values of h at the cell centers (xbc) at t = tstart -dt 
c           pubc_init - values of u at the cell centers (xbc) at t = tstart -dt  
c           tlist - list of times to print solutions at
c           tlist_len - length of tlist
c           ExpWdir - string containing directory to print the individual experiments results to
c           expwdirlen - length of the string ExpWdir
c Outputs:
c           currenttime - current time of numerical solution
c           hbc_fin - values of h at the cell centers (xbc) at t = currenttime 
c           ubc_fin - values of h at the cell centers (xbc) at t = currenttime 
c           Energ_Init - Array of conserved quantities (h,uh,Energy) at t = tstart
c           Energ_Fin - Array of conserved quantities (h,uh,Energy) at t = currenttime

c ===== 
      SUBROUTINE NumericalSolveTSPrint(tstart,tend,
     . ga,beta1,beta2,dx,dt,n_GhstCells,xbc_len,
     . xbc,hbc_init,ubc_init,
     . phbc_init,pubc_init,
     . tlist,tlist_len,ExpWdir,expwdirlen,
     . currenttime,hbc_fin,ubc_fin,
     . Energ_Init, Energ_Fin)
     
     
      IMPLICIT NONE
      !Inputs  
      INTEGER n_GhstCells,xbc_len,expwdirlen,tlist_len
      CHARACTER(LEN=expwdirlen) ExpWdir
      DOUBLE PRECISION tstart,tend,ga,beta1,beta2,
     . dx,dt,currenttime
      DOUBLE PRECISION xbc(xbc_len),hbc_init(xbc_len),
     . ubc_init(xbc_len),phbc_init(xbc_len),
     . pubc_init(xbc_len),
     . tlist(tlist_len)
     
      !Outputs  
      DOUBLE PRECISION Energ_Init(3), Energ_Fin(3),hbc_fin(xbc_len),
     . ubc_fin(xbc_len)
      
      !Local variables - Loop variable i, filecount and strct ensure that intermediate time prints write to different files     
      INTEGER i,filecount
      CHARACTER(LEN=2) strct
      DOUBLE PRECISION phbc_fin(xbc_len),pubc_fin(xbc_len)
      
      !initial time
      currenttime  = tstart
      filecount = 1
      
      !Set current values of h,u,G to initial conditions (intialise)
      DO i = 1,xbc_len
         hbc_fin(i) = hbc_init(i) 
         ubc_fin(i) = ubc_init(i)
         phbc_fin(i) = phbc_init(i) 
         pubc_fin(i) = pubc_init(i) 
      END DO
      
      !calculate initial Energies  
      CALL TotalEnergy(xbc_len,n_GhstCells,dx,ga,beta1,beta2,
     . hbc_init,ubc_init,Energ_Init)

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
               WRITE(8,*) currenttime,xbc(i),hbc_fin(i),ubc_fin(i)
            END DO
            CLOSE(8)
            
            filecount = filecount + 1
         END IF

         !Evolve system
         CALL EvolveWrap(ga,beta1,beta2,dx,dt,
     . xbc_len,n_GhstCells,hbc_fin,ubc_fin,phbc_fin,pubc_fin) 
         !Update current time and print, so we can get some feedback about what its doing
         currenttime  = currenttime  + dt
         PRINT *, 'Current Time : ', currenttime 
      END DO
       
      !Calculate final energies  
      ! u solve

      CALL TotalEnergy(xbc_len,n_GhstCells,dx,ga,beta1,beta2,
     . hbc_fin,ubc_fin,
     . Energ_Fin)

                 
      END

c============
c  Evolution Routines
c  Functions that evolve h^n, u^n  -> h^{n+1}, u^{n+1}
c  Using finite difference approximations for all derivatives
c============

c ====
c Subroutine that evolves h^n to h^{n+1} and u^n to u^{n+1}
c Inputs:  
c         ga - acceleration due to gravity
c         beta1 - beta1 value
c         beta2 - beta2 value
c         dx  - cell width   
c         dt - time step
c         xbc_len - total number of cells in discretisation
c         n_GhstCells - number of ghost cells
c Inputs/Outputs:
c         hbc -  In : h values at cell centers at current time (h^{n})
c                Out: h values at cell centers at next time (h^{n+1})
c         phbc - In : h values at cell centers at previous time (h^{n-1})
c                Out: h values at cell centers at current time (h^{n})
c         ubc -  In : u values at cell centers at current time (u^{n})
c                Out: u values at cell centers at next time (u^{n+1})
c         pubc - In : u values at cell centers at previous time (u^{n-1})
c                Out: u values at cell centers at current time (u^{n})
c ====
      SUBROUTINE EvolveWrap(ga,beta1,beta2,dx,dt,
     . xbc_len,n_GhstCells,hbc,ubc,phbc,pubc)   
     
      !Inputs
      INTEGER n_GhstCells,xbc_len
      DOUBLE PRECISION beta1,beta2,dx,dt,ga
     
      !Input/Outputs
      DOUBLE PRECISION hbc(xbc_len),ubc(xbc_len),
     . phbc(xbc_len),pubc(xbc_len)
     
      !Local variables
      DOUBLE PRECISION nhbc(xbc_len),nubc(xbc_len) 
      INTEGER i
      
      CALL Evolveh(hbc,ubc,phbc,dx,dt,xbc_len,n_GhstCells,nhbc)
      CALL Evolveu(hbc,ubc,pubc,beta1,beta2,ga,dx,dt,xbc_len,
     . n_GhstCells,nubc)
     
     
      !Return updated h,u,ph,pu
      DO i = 1, xbc_len
         phbc(i) = hbc(i)
         pubc(i) = ubc(i)
         hbc(i) = nhbc(i)
         ubc(i) = nubc(i)
      END DO
      
      END   

c ====
c Subroutine that evolves h^n to h^{n+1}
c Inputs:  
c         hbc - h values at cell centers at current time
c         ubc - u values at cell centers at current time
c         phbc - h values at cell centers at previous time
c         dx  - cell width   
c         dt - time step
c         xbc_len - total number of cells in discretisation
c         n_GhstCells - number of ghost cells
c Outputs:
c         nhbc - u values at cell centers at next time
c ====
      SUBROUTINE Evolveh(hbc,ubc,phbc,dx,dt,xbc_len,n_GhstCells,nhbc)
   
      !Inputs
      INTEGER n_GhstCells,xbc_len
      DOUBLE PRECISION hbc(xbc_len),ubc(xbc_len),phbc(xbc_len)
      DOUBLE PRECISION dx,dt
     
      !Outputs
      DOUBLE PRECISION nhbc(xbc_len)
     
      !Local variables
      INTEGER i,j    
      DOUBLE PRECISION chx,cux
      
      !interior
      DO i = n_GhstCells + 1, xbc_len - n_GhstCells
        chx = (hbc(i+1) - hbc(i-1)) / (2*dx)
        cux = (ubc(i+1) - ubc(i-1))/ (2*dx)
        nhbc(i) = phbc(i) - 2*dt*(hbc(i)*cux + ubc(i)*chx)
      END DO
      
      ! boundaries
      DO i = 1,n_GhstCells
        nhbc(i) = hbc(i)
        j = xbc_len - n_GhstCells + i
        nhbc(j) = hbc(j)
      END DO
      
      
      END


c ====
c Subroutine that evolves u^n to u^{n+1}
c Inputs:  
c         hbc - h values at cell centers at current time
c         ubc - u values at cell centers at current time
c         pubc - u values at cell centers at previous time
c         beta1 - beta1 value of gSGN equations
c         beta2 - beta2 value of gSGN equations
c         ga - acceleration due to gravity
c         dx  - cell width   
c         dt - time step
c         xbc_len - total number of cells in discretisation
c         n_GhstCells - number of ghost cells
c Outputs:
c            nubc - u values at cell centers at next time
c ====
      
      SUBROUTINE EvolveU(hbc,ubc,pubc,beta1,beta2,ga,dx,dt,xbc_len,
     . n_GhstCells,nubc)
     
      !Inputs
      INTEGER n_GhstCells,xbc_len
      DOUBLE PRECISION hbc(xbc_len),ubc(xbc_len),pubc(xbc_len)
      DOUBLE PRECISION beta1,beta2,dx,dt,ga
     
      !Outputs
      DOUBLE PRECISION nubc(xbc_len)
     
      !Local variables
      INTEGER i
      DOUBLE PRECISION subdiag1(xbc_len), diag(xbc_len),
     . supdiag1(xbc_len), RHS(xbc_len)      
      DOUBLE PRECISION unterms,
     . chx,cux,chxx,cuxx,chxxx,cuxxx,pux,puxx     
     
      !Loop over interior of cells, building matrix A such that u^{n+1} = Au^{n}
      DO i = n_GhstCells + 1, xbc_len - n_GhstCells

         chx = 0.5*(hbc(i+1) - hbc(i-1))/dx
         cux = 0.5*(ubc(i+1) - ubc(i-1))/dx
         pux = 0.5*(pubc(i+1) - pubc(i-1))/dx
         chxx = (hbc(i+1) - 2*hbc(i) + hbc(i-1))/ (dx**2)
         cuxx = (ubc(i+1) - 2*ubc(i) + ubc(i-1))/ (dx**2)
         puxx = (pubc(i+1) - 2*pubc(i) + pubc(i-1))/ (dx**2)
         chxxx = (hbc(i+2) - 2*hbc(i+1)  + 2*hbc(i-1) - hbc(i-2))/
     .     (2*dx**3)
         cuxxx = (ubc(i+2) - 2*ubc(i+1)  + 2*ubc(i-1) - ubc(i-2))/
     .     (2*dx**3)
     
         unterms = 3*beta1*chx*cux**2*hbc(i)/2 
     .    - 3*beta1*chx*cuxx*hbc(i)*ubc(i)/2 
     .    + beta1*cux*cuxx*hbc(i)**2/2 
     .    - beta1*cuxxx*hbc(i)**2*ubc(i)/2 
     .    - beta2*chx**3*ga/2 
     .    - 2*beta2*chx*chxx*ga*hbc(i) 
     .    - beta2*chxxx*ga*hbc(i)**2/2
     .    + chx*ga 
     .    + cux*ubc(i)
         
         unm1terms = 3*beta1*chx*hbc(i)*pux/(4*dt) +
     .     beta1*hbc(i)**2*puxx/(4*dt) - pubc(i)/(2*dt)
         
         subdiag1(i) = beta1*hbc(i)*(3*chx*dx - 2*hbc(i))/(8*dt*dx**2)
         diag(i) = (beta1*hbc(i)**2 + dx**2)/(2*dt*dx**2)
         supdiag1(i) = -beta1*hbc(i)*(3*chx*dx + 2*hbc(i))/(8*dt*dx**2)
     
         RHS(i) = -unm1terms -unterms
     
      END DO
      
      !Boundary values
      DO i = 1, n_GhstCells
         ! left
         subdiag1(i) = 0d0
         diag(i) = 1d0
         supdiag1(i) = 0d0
         RHS(i) = ubc(i)
         
         !right
         subdiag1(xbc_len - n_GhstCells + i) = 0d0
         diag(xbc_len - n_GhstCells + i) = 1d0
         supdiag1(xbc_len - n_GhstCells + i)= 0d0

         RHS(xbc_len - n_GhstCells + i) = 
     .      ubc(xbc_len - n_GhstCells + i)
      END DO
      
      CALL ThomasTriSolve(xbc_len,subdiag1,diag,supdiag1,RHS,nubc)
      
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
c            ga - acceleration due to gravity
c            beta1 - beta1 gSGN parameter
c            beta2 - beta2 gSGN parameter
c            j - cell index
c            dx - cell width
c Outputs:
c            CellEnergies - array of total amount of quantities in cell (h,G,uh,Energy)
c ===== 
      SUBROUTINE AllEnergiesIntegralCell(xbc_len,h,u,ga,beta1,beta2,j
     . ,dx,CellEnergies)
      
      !Input
      INTEGER j,xbc_len
      DOUBLE PRECISION h(xbc_len),u(xbc_len)
      DOUBLE PRECISION dx,ga,beta1,beta2
      
      !Output
      DOUBLE PRECISION CellEnergies(3)
      
      !Local variable
      INTEGER i
      DOUBLE PRECISION fGPe(3),sGPe(3),tGPe(3)
      DOUBLE PRECISION GPmxj,hGP,uGP,uxGP,hxGP
      DOUBLE PRECISION hCoeff(5), uCoeff(5)
      
      ! Coefficients of quatrics that interpolate h,u,G
      CALL QuarticInterp(h(j-2),h(j-1),h(j),h(j+1),h(j+2),dx,hCoeff)
      CALL QuarticInterp(u(j-2),u(j-1),u(j),u(j+1),u(j+2),dx,uCoeff)
      
      !Location of First Gauss Point
      GPmxj = -dx*DSQRT(3.0d0/5.0d0)/2
      
      !Evaluate required functions at gauss point (h,dh/dx,G,u,du/dx)
      CALL QuarticCoeffEvalxj(hCoeff,GPmxj,hGP)
      CALL QuarticCoeffEvalGradxj(hCoeff,GPmxj,hxGP)
      CALL QuarticCoeffEvalxj(uCoeff,GPmxj,uGP)
      CALL QuarticCoeffEvalGradxj(uCoeff,GPmxj,uxGP)
      
      !Combine into array of values at first Gauss point
      fGPe(1) = hGP
      fGPe(2) = hGP*uGP
      fGPe(3) = (hGP*uGP**2 + beta1/2d0*(uxGP**2)*(hGP**3)
     . + ga*hGP**2*(1d0 + beta2/2d0*hxGP**2 )   )/2d0

      !Location of Second Gauss Point
      GPmxj = 0.0 
      
      !Evaluate required functions at gauss point (h,dh/dx,G,u,du/dx)
      CALL QuarticCoeffEvalxj(hCoeff,GPmxj,hGP)
      CALL QuarticCoeffEvalGradxj(hCoeff,GPmxj,hxGP)
      CALL QuarticCoeffEvalxj(uCoeff,GPmxj,uGP)
      CALL QuarticCoeffEvalGradxj(uCoeff,GPmxj,uxGP)
      
      !Combine into array of values at second Gauss point
      sGPe(1) = hGP
      sGPe(2) = hGP*uGP
      sGPe(3) = (hGP*uGP**2 +  beta1/2d0*(uxGP**2)*(hGP**3)
     . + ga*hGP**2*(1d0 + beta2/2d0*hxGP**2 )   )/2d0
      
      !Location of Third Gauss Point
      GPmxj = dx*DSQRT(3.0d0/5.0d0)/2
      
      !Evaluate required functions at gauss point (h,dh/dx,G,u,du/dx)
      CALL QuarticCoeffEvalxj(hCoeff,GPmxj,hGP)
      CALL QuarticCoeffEvalGradxj(hCoeff,GPmxj,hxGP)
      CALL QuarticCoeffEvalxj(uCoeff,GPmxj,uGP)
      CALL QuarticCoeffEvalGradxj(uCoeff,GPmxj,uxGP)
      
      !Combine into array of values at third Gauss point
      tGPe(1) = hGP
      tGPe(2) = hGP*uGP
      tGPe(3) = (hGP*uGP**2 + beta1/2d0*(uxGP**2)*(hGP**3)
     . + ga*hGP**2*(1d0 + beta2/2d0*hxGP**2 )   )/2d0
      
      !Combine function evaluations to approximate integral using Gaussian quadrature.
      DO i = 1,3
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
c Outputs:
c            TotEnergVals - array of total amount of quantities in cell (h,uh,Energy)
c ===== 
      SUBROUTINE TotalEnergy(xbc_len,n_GhstCells,dx,ga,beta1,beta2,
     . hbc,ubc,
     . TotEnergVals)
      
      !Inputs
      INTEGER xbc_len,n_GhstCells
      DOUBLE PRECISION hbc(xbc_len),ubc(xbc_len)
      DOUBLE PRECISION dx,ga,beta1,beta2
      
      !Outputs
      DOUBLE PRECISION TotEnergVals(3)
      
      !Local variables
      DOUBLE PRECISION CellEnergVals(3)
      INTEGER i,j
            
      !running totals for energy values, start at 0
      DO i = 1,3
         TotEnergVals(i) = 0.0
      END DO
      
      !just loop over interior of hbc, ubc which have interior values + ghost cell values
      DO j= n_GhstCells + 1, xbc_len - n_GhstCells
         CALL  AllEnergiesIntegralCell(xbc_len,hbc,ubc,ga,
     .      beta1,beta2,j,dx,CellEnergVals)
     
         !add cell energy value to running total
         DO i = 1,3
            TotEnergVals(i) = TotEnergVals(i) + CellEnergVals(i)
         END DO
      END DO
      
      END

      


