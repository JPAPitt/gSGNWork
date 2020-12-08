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
c         a0 - background water depth on which travelling wave propagates
c         a1 - amplitude of travelling wave in h
c         a2 - translation speed of wave (wave is function of x - a2*t)
c         a3 - width of travelling wave
c         a4 - amplitude of travelling wave in u
c         a5 - beta parameter : speed (betas are a function of x - a5*t)
c         t -  current time
c         xbc - array of cell centers
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
      SUBROUTINE Evolveh(a0,a1,a2,a3,a4,t,xbc,
     . hbc,ubc,phbc,dx,dt,xbc_len,n_GhstCells,nhbc)
   
      !Inputs
      INTEGER n_GhstCells,xbc_len
      DOUBLE PRECISION xbc(xbc_len),hbc(xbc_len),ubc(xbc_len),
     . phbc(xbc_len)
      DOUBLE PRECISION a0,a1,a2,a3,a4,t,dx,dt
     
      !Outputs
      DOUBLE PRECISION nhbc(xbc_len)
     
      !Local variables
      INTEGER i,j    
      DOUBLE PRECISION chx,cux,dhdt,dfluxhdx
      
      !interior
      DO i = n_GhstCells + 1, xbc_len - n_GhstCells
        chx = (hbc(i+1) - hbc(i-1)) / (2*dx)
        cux = (ubc(i+1) - ubc(i-1))/ (2*dx)
        
        CALL Forcedht(xbc(i),t,a1,a2,a3,dhdt) 
        CALL Forcedfluxhx(xbc(i),t,a0,a1,a2,a3,a4,dfluxhdx) 
        nhbc(i) = phbc(i) - 2*dt*(hbc(i)*cux + ubc(i)*chx)
     .   + 2*dt*(dhdt + dfluxhdx)
      END DO
      
      ! boundaries
      DO i = 1,n_GhstCells      
        nhbc(i) = a0 + a1*DEXP(-(xbc(i) - a2*t)**2 / (2*a3))
        j = xbc_len - n_GhstCells + i
        nhbc(j) = a0 + a1*DEXP(-(xbc(j) - a2*t)**2 / (2*a3))
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
     
         CALL Forcedut(x,t,a2,a3,a4,dudt) 
         CALL Forcedfluxux(x,t,ga,a0,a1,a2,a3,a4,a5,
     . b1a6,b1a7,b2a6,b2a7,dfluxudx) 
         
         RHS(i) = (dudt +dfluxudx)  -unm1terms -unterms
     
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
c Subroutine that calculates analytic value of du/dt for Gaussian forced solution
c Inputs: 
c         x - location to calculate flux dh/dt at
c         t - time to calculate flux dh/dt at
c         a2 - translation speed of wave (wave is function of x - a2*t)
c         a3 - width of travelling wave
c         a4 - amplitude of travelling wave in u
c Outputs:
c         dudt - analytic value of du/dt 
c ====
      SUBROUTINE Forcedut(x,t,a2,a3,a4,dudt) 
      !Input
      DOUBLE PRECISION x,t,a2,a3,a4
      !Output
      DOUBLE PRECISION dGdt
      !Local variables
      DOUBLE PRECISION :: EXPPHI1,PHI,ui
      
      PHI  = x - a2*t
      EXPPHI1 = DEXP(-PHI**2 / (2*a3))
      ui  = a4*EXPPHI1
      dudt = PHI*a2*ui/a3
      
      END 

c ====
c Subroutine that calculates analytic value of df(u,h)/dx for u for Gaussian forced solution
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
c         dfluxudx - analytic value of df(u,h)/dx
c ====
      SUBROUTINE Forcedfluxux(x,t,ga,a0,a1,a2,a3,a4,a5,
     . b1a6,b1a7,b2a6,b2a7,dfluxudx) 
      !Input
      DOUBLE PRECISION x,t,ga,a0,a1,a2,a3,a4,a5,b1a6,b1a7,b2a6,b2a7
      !Output
      DOUBLE PRECISION dfluxudx
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
      
      dfluxudx = -b1a6*d2udx2*hi**2*ui/2 
     . + b1a6*dudx**2*hi**2/2 
     . - b2a6*d2hdx2*ga*hi**2/2 
     . - b2a6*dhdx**2*ga*hi/4 
     . - 3*beta1*d2udx2*dhdx*hi*ui/2 
     . + beta1*d2udx2*dudx*hi**2/2 
     . - beta1*d3udx3*hi**2*ui/2 
     . + 3*beta1*dhdx*dudx**2*hi/2 
     . - 2*beta2*d2hdx2*dhdx*ga*hi 
     . - beta2*d3hdx3*ga*hi**2/2 
     . - beta2*dhdx**3*ga/2 
     . + dhdx*ga 
     . + dudx*ui
      
      END  
     


