c ====================================================================================
c Main program that runs a series of experiments to validate numerical solution using the Gaussian forced solution which occurs
c for all beta values, and even varying beta, as done in this example. 
c ====================================================================================        

c ====
c Subroutine that runs a numerical experiment with all the parameters, 
c the subroutine generates discretisation, initial conditions and then solves the initial condition problem
c it then prints out the desired outputs h,u,G at cell centers, parameter values and convergence and conservation errors
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
c           ga - acceleration due to gravity
c           xstart - first cell center in interior
c           xbc_len - total number of cells in discretisation
c           n_GhstCells - number of ghost cells at each boundary
c           dx - cell width
c           tstart - start time
c           tend - time to produce numerical solution to
c           dt - time step size
c           NormFile - number for Norm file (so individual experiments can write convergence error to the same file)
c           EnergFile- number for Energy file (so individual experiments can write conservation error to the same file)
c           tlist - list of times to print solutions at
c           tlist_len - length of tlist
c           ExpWdir - string containing directory to print the individual experiments results to
c           ExpWdir_len - length of the string ExpWdir
c ===== 
     
      SUBROUTINE SingleForcedGaussian(a0,a1,a2,a3,a4,a5,
     . b1a6,b1a7,b2a6,b2a7,ga,
     . xstart,xbc_len,n_GhstCells,dx,tstart,tend,
     . dt,NormFile,tlist,tlist_len,
     . ExpWdir,ExpWdir_len)

      IMPLICIT NONE
      !Inputs  
      DOUBLE PRECISION a0,a1,a2,a3,a4,a5,b1a6,b1a7,b2a6,b2a7,
     . ga,xstart,dx,
     . tstart,tend,dt
      INTEGER xbc_len,n_GhstCells,ExpWdir_len,NormFile,
     . tlist_len
      CHARACTER(LEN=ExpWdir_len) ExpWdir
      DOUBLE PRECISION tlist(tlist_len)
      
      !Local Variables
      DOUBLE PRECISION xbc(xbc_len),hbc_init(xbc_len),
     . Gbc_init(xbc_len),ubc_init(xbc_len),hbc_fin(xbc_len),
     . Gbc_fin(xbc_len),ubc_fin(xbc_len),hbc_fin_a(xbc_len),
     . Gbc_fin_a(xbc_len),ubc_fin_a(xbc_len),
     . Norms(3),sumerr(3),suma(3)
      INTEGER i
      DOUBLE PRECISION currenttime
      
      
      !generate cell center discretisation array - xbc
      CALL Generatexbc(xstart,dx,xbc_len,n_GhstCells,xbc)
      
      !Generate initial conditions of h,u and G at cell centers - xbc
      !hbc_init,ubc_init,Gbc_init will remain as the intial conditions values
      CALL ForcedGaussian(xbc,xbc_len,tstart,
     . a0,a1,a2,a3,a4,a5,b1a6,b1a7,hbc_init,ubc_init,Gbc_init)
     
      
      !Numerically Solve gSGN with beta values until currenttime > tend     
      CALL NumericalSolveForcedTSPrint(a0,a1,a2,a3,a4,a5,
     . b1a6,b1a7,b2a6,b2a7,
     . tstart,tend,ga,dx,dt,n_GhstCells,xbc_len,xbc,
     . hbc_init,Gbc_init,ubc_init,
     . tlist,tlist_len,ExpWdir,ExpWdir_len,
     . currenttime,hbc_fin,Gbc_fin,ubc_fin)

      !Get analytic values of h,u,G at currenttime which is the time the numerical solution is at
      CALL ForcedGaussian(xbc,xbc_len,currenttime,
     . a0,a1,a2,a3,a4,a5,b1a6,b1a7,hbc_fin_a,ubc_fin_a,Gbc_fin_a)
     
      ! Measure convergence error - L2      
      ! Initialise both numerator (h' - h)  and denominator h in relative L2 error = (h' - h) /h
      ! (where h' is numerical solution and h is analytic solution)
      DO i = 1,3
         sumerr(i) = 0
         suma(i) = 0
      END DO
      
      !Sum together numerator and denominator of relative L2 error
      DO i = 1, xbc_len
         sumerr(1) = sumerr(1) + (hbc_fin_a(i) - hbc_fin(i))**2
         suma(1) = suma(1) + (hbc_fin_a(i))**2
         
         sumerr(2) = sumerr(2) + (ubc_fin_a(i) - ubc_fin(i))**2
         suma(2) = suma(2) + (ubc_fin_a(i))**2
         
         sumerr(3) = sumerr(3) + (Gbc_fin_a(i) - Gbc_fin(i))**2
         suma(3) = suma(3) + (Gbc_fin_a(i))**2 
      END DO
      
      DO i = 1,3
         !if denominator is very small, just return absolute error
         IF (suma(i) .LT. 10d0**(-10)) THEN
            Norms(i) = DSQRT(sumerr(i))
         ELSE
            Norms(i) = DSQRT(sumerr(i)) / DSQRT(suma(i))
         END IF
      END DO
      
      ! Open files for result printing
      OPEN(1, FILE = ExpWdir//'Init.dat') 
      OPEN(2, FILE = ExpWdir//'End.dat') 
      OPEN(3, FILE = ExpWdir//'EndAna.dat') 
      OPEN(4, FILE = ExpWdir//'Params.dat') 
      
      ! Write out initial conditions,numerical and analytic solutions
      DO i = 1,xbc_len
         WRITE(1,*) tstart,xbc(i),hbc_init(i),Gbc_init(i),ubc_init(i)
         WRITE(2,*) currenttime,xbc(i),hbc_fin(i),Gbc_fin(i),ubc_fin(i)
         WRITE(3,*) currenttime,xbc(i),hbc_fin_a(i),Gbc_fin_a(i),
     .               ubc_fin_a(i)
      END DO
      
      ! Write out parameters used to generate experiment
      WRITE(4,*) 'xstart',xstart
      WRITE(4,*) 'xend',xbc(xbc_len)
      WRITE(4,*) 'x_len',xbc_len - 2*n_GhstCells
      WRITE(4,*) 'n_GhstCells',n_GhstCells
      WRITE(4,*) 'xbc_len',xbc_len
      WRITE(4,*) 'dx' , dx
      WRITE(4,*) 'tstart', tstart
      WRITE(4,*) 'tend',tend 
      WRITE(4,*) 'actual_end_time', currenttime
      WRITE(4,*) 'dt' , dt
      WRITE(4,*) 'gravity' , ga
      WRITE(4,*) 'a0' , a0
      WRITE(4,*) 'a1' , a1
      WRITE(4,*) 'a2' , a2
      WRITE(4,*) 'a3' , a3
      WRITE(4,*) 'a4' , a4
      WRITE(4,*) 'a5' , a5
      WRITE(4,*) 'b1a6' , b1a6
      WRITE(4,*) 'b1a7' , b1a7
      WRITE(4,*) 'b2a6' , b2a6
      WRITE(4,*) 'b2a7' , b2a7
      
      
      !Write out convergence error to combined file NormFile
      WRITE(NormFile,*) dx,Norms(1),Norms(2),Norms(3)
      
      
      CLOSE(1)
      CLOSE(2)
      CLOSE(3)
      CLOSE(4)
      END
      
c ====
c  Main program
c  That will loop over values of dx, to investigate convergence and conservation properties for analytic solution (Solitary travelling wave)
c  to validate the numerical method
c ===== 
 
      PROGRAM main
         
      IMPLICIT NONE
   
      ! wdirlen - maximum length for directory string
      ! NormFile,EnergFile - integers for files that combine the convergence and conservation errors
      ! tlist_len - list of what t values to print intermediate solutions at
      INTEGER wdirlen,NormFile,EnergFile,tlist_len
      PARAMETER(wdirlen= 300,NormFile = 98, EnergFile = 99,
     . tlist_len=11)
      
      ! wdir - directory address to write out to
      ! strdiri - used to produce directory names for the individual experiments
      CHARACTER(LEN =wdirlen) wdir
      CHARACTER(LEN =2) strdiri
     
      !expi - current experiment number 
      !x_len - number of cells in interior for current experiment
      !xbc_len - total number of cells for current experiment
      ! n_GhstCells- number of ghost cells at boundaries for all experiments (will depend on reconstruction, so doesnt need to change across experiments)
      ! lowestresx - number of cells when expi = 0 (the coarsest resolution) which is then doubled to see what happens when dx -> 0.
      INTEGER expi,x_len,xbc_len,n_GhstCells, lowestresx
      
      !a0 - background water depth on which the travelling wave propagated
      !a1 - amplitude of travelling wave
      !ga - acceleration due to gravity
      !xstart - center of first cell in interior
      !xend - center of last cell in interior
      !tstart - time of initial condition
      !tend - lower bound for final time of simulation
      !dx - cell width
      !dt - time step size
      !Cr - ensuring CFL satisfied
      !maxwavespeed - estimate for maximum wavespeed to satisfy CFL
      !beta1 - gSGN beta1 parameter
      !beta2 - gSGN beta2 parameter
      !alpha - used to bound wave speed when beta2 > beta1 (advancing dispersive waves)
      DOUBLE PRECISION a0,a1,a2,a3,a4,a5,b1a6,b1a7,b2a6,b2a7,
     . ga,xstart,xend,tstart,tend,
     . dx,dt,Cr,maxwavespeed,alpha,
     . beta1lb,beta1rb,beta1min,beta2lb,beta2rb,beta2max
     
      !list of intermediate times to print solution at
      DOUBLE PRECISION tlist(tlist_len)
     
      !effective length of directory address string
      INTEGER effeclenwdir
      
      wdir = "../Results/Validation/" // 
     . "ForcedSolution/GaussianISGN/" 
     
      !use LenTrim routine to calculate effeclenwdir
      CALL LenTrim(wdir,wdirlen,effeclenwdir)
      
      !Refresh directories dump data
      CALL SYSTEM('rm -rf '//wdir)
      CALL SYSTEM('mkdir -p '// wdir)
      
      !open Conservation and Convergence output files
      OPEN(NormFile, FILE = wdir(1:effeclenwdir)//'Norms.dat') 

      
      n_GhstCells = 6
      
      ga = 9.81d0
      
      !h and u parameters
      a0 = 1d0
      a1 = 0.5d0
      a2 = 5d0
      a3 = 20d0
      a4 = 0.3d0
      
      !beta parameters
      ! Must ensure that beta1 > 0 and beta2 >0 for all values. 
      ! So we have : beta1 = b1a6*(x(i) - a5*t) + b1a7
      ! and :        beta2 = b2a6*(x(i) - a5*t) + b2a7
      a5 = -0.1
      b1a6 = 0.005
      b1a7 = 1.0
      b2a6 = 0.006
      b2a7 = 1.1

      xstart = -100d0
      xend = 100d0
      
      
      tstart = 0d0
      tend = 10.0d0
      
      ! tlist is every 1s between 0s and 30s (inclusive)
      CALL EqualSpaced(tstart,tlist_len,1.0, tlist)
            
      lowestresx = 100
      
      !Perform the Solitary Travelling Wave experiment a number of times, decreasing \Delta x each time
      DO expi = 0,14
      
         !make directory for inidividual experiment
         WRITE (strdiri,'(I2.2)') expi
         CALL SYSTEM('mkdir -p '// wdir(1:effeclenwdir) //strdiri)
         
         !calculate number of cells, both interior and total
         x_len = lowestresx *2**(expi)
         xbc_len = x_len + 2 *n_GhstCells
         
         !alpha will now be a function of time
         !this alpha value is used to get CFL condition
         ! We can ensure that we pick the maximum alpha value by noticing that beta1 and beta2 are monotonic functions in (x - a5*t) where a5 > 0
         !So we  can pick the maximum possible beta2 value, and the smallest beta1 value, so that we bound alpha(t) for the CFL condition
         beta1lb = b1a6*MIN(xstart - a5*tend,xstart) + b1a7
         beta1rb = b1a6*MAX(xend - a5*tend,xend) + b1a7
         beta1min = MIN(beta1lb,beta1rb)
         
         beta2lb = b2a6*MIN(xstart - a5*tend,xstart) + b2a7
         beta2rb = b2a6*MAX(xend - a5*tend,xend) + b2a7
         beta2max = MAX(beta2lb,beta2rb)
         IF  (DABS(beta1min) < 10d0**(-10))  THEN
            alpha = 1d0
         ELSE
            alpha = MAX(1d0,beta2max / beta1min)
         END IF
         
         !calculate dx to obtain the right number of cells
         !then calculate dt using CFL condition with Cr
         dx = (xend - xstart) / (x_len -1)
         Cr = 0.5
         !really want to ensure that CFL condition isn't violated, so not only do I include max(u) and max(h), but also the 'speed' of the guassian a2 in (x - a2*t)
         maxwavespeed = a2 + a4 + DSQRT(alpha*ga*(a0 + a1))
         dt  = (Cr / maxwavespeed) *dx
         
         PRINT *,'++++++++++++ Experiment : ',expi ,' || ', '# Cells :',
     .    x_len , '++++++++++++'
         
         !Call routine that generates and then runs individual experiment
         CALL SingleForcedGaussian(a0,a1,a2,a3,a4,a5,
     .    b1a6,b1a7,b2a6,b2a7,ga,
     .    xstart,xbc_len,n_GhstCells,dx,tstart,tend,
     .    dt,NormFile,tlist,tlist_len,
     .    wdir(1:effeclenwdir)//strdiri//'/',effeclenwdir+3)

      
      END DO
      
      CLOSE(NormFile)
      
      END
