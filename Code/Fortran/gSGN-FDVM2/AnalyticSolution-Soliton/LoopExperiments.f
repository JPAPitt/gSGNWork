c ====================================================================================
c Main program that runs a series of experiments to validate numerical solution using the analytic solitary travelling
c wave solution which occurs when beta1= 2/3 and beta2 = 0.
c ====================================================================================        

c ====
c Subroutine that runs a numerical experiment with all the parameters, 
c the subroutine generates discretisation, initial conditions and then solves the initial condition problem
c it then prints out the desired outputs h,u,G at cell centers, parameter values and convergence and conservation errors
c
c Inputs:   a0 - background water depth on which travelling wave propagates
c           a1 - amplitude of travelling wave
c           ga - acceleration due to gravity
c           beta1 - beta1 value of gSGN equations
c           beta2 - beta2 value of gSGN equations
c           xstart - first cell center in interior
c           xbc_len - total number of cells in discretisation
c           n_GhstCells - number of ghost cells at each boundary
c           dx - cell width
c           tstart - start time
c           tend - time to produce numerical solution to
c           dt - time step size
c           theta - reconstruction parameter theta for generalised minmod limiter
c           NormFile - number for Norm file (so individual experiments can write convergence error to the same file)
c           EnergFile- number for Energy file (so individual experiments can write conservation error to the same file)
c           tlist - list of times to print solutions at
c           tlist_len - length of tlist
c           ExpWdir - string containing directory to print the individual experiments results to
c           ExpWdir_len - length of the string ExpWdir
c ===== 
      SUBROUTINE SingleSerreSoliton(a0,a1,ga,beta1,beta2,
     . xstart,xbc_len,n_GhstCells,dx,tstart,tend,
     . dt,theta,NormFile,EnergFile,tlist,tlist_len,
     . ExpWdir,ExpWdir_len)
     
     
     
      IMPLICIT NONE
      !Inputs  
      DOUBLE PRECISION a0,a1,ga,beta1,beta2,xstart,dx,
     . tstart,tend,dt,theta
      INTEGER xbc_len,n_GhstCells,ExpWdir_len,NormFile,EnergFile,
     . tlist_len
      CHARACTER(LEN=ExpWdir_len) ExpWdir
      DOUBLE PRECISION tlist(tlist_len)
      
      !Local Variables
      DOUBLE PRECISION xbc(xbc_len),hbc_init(xbc_len),
     . Gbc_init(xbc_len),ubc_init(xbc_len),hbc_fin(xbc_len),
     . Gbc_fin(xbc_len),ubc_fin(xbc_len),hbc_fin_a(xbc_len),
     . Gbc_fin_a(xbc_len),ubc_fin_a(xbc_len),
     . Energs_init(4), Energs_fin(4),Energ_Err(4),
     . Norms(3),sumerr(3),suma(3)
      INTEGER i
      DOUBLE PRECISION currenttime
      
      
      !generate cell center discretisation array - xbc
      CALL Generatexbc(xstart,dx,xbc_len,n_GhstCells,xbc)
      
      !Generate initial conditions of h,u and G at cell centers - xbc
      !hbc_init,ubc_init,Gbc_init will remain as the intial conditions values
      CALL SGNSolitaryTravellingWave(xbc,xbc_len,tstart,
     . ga,a0,a1,hbc_init,ubc_init,Gbc_init)
     
      
      !Numerically Solve gSGN with beta values until currenttime > tend     
      CALL NumericalSolveTSPrint(tstart,tend,
     . ga,beta1,beta2,theta,dx,dt,n_GhstCells,xbc_len,xbc,
     . hbc_init,Gbc_init,ubc_init,
     . tlist,tlist_len,ExpWdir,ExpWdir_len,
     . currenttime,hbc_fin,Gbc_fin,ubc_fin,
     . Energs_init, Energs_fin)

      !Get analytic values of h,u,G at currenttime which is the time the numerical solution is at
      call SGNSolitaryTravellingWave(xbc,xbc_len,currenttime,ga,a0,a1,
     . hbc_fin_a,ubc_fin_a,Gbc_fin_a)
     
     
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

      !Conservation Norm Tests  - C1 
      DO i = 1,4
         !if denominator small just return absolute error
         IF (DABS(Energs_init(i)) .LT. 10d0**(-10)) THEN
            Energ_Err(i) = DABS(Energs_fin(i) - Energs_init(i))
         ELSE
            Energ_Err(i) = DABS(Energs_fin(i) - Energs_init(i))/
     .       DABS(Energs_init(i))
         END IF
      END DO   

           
      ! Open files for result printing
      OPEN(1, FILE = ExpWdir//'Init.dat') 
      OPEN(2, FILE = ExpWdir//'End.dat') 
      OPEN(3, FILE = ExpWdir//'EndAna.dat') 
      OPEN(4, FILE = ExpWdir//'Params.dat') 
      OPEN(5, FILE = ExpWdir//'Energy.dat') 
      
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
      WRITE(4,*) 'theta' , theta
      WRITE(4,*) 'gravity' , ga
      WRITE(4,*) 'a0' , a0
      WRITE(4,*) 'a1' , a1
      WRITE(4,*) 'beta1' , beta1
      WRITE(4,*) 'beta2' , beta2
      
      ! Write out Initial and Final Energy
      WRITE(5,*) 'When , h , G ,  uh , Energy'
      WRITE(5,*) 'Initial ',Energs_init(1),Energs_init(2),
     .   Energs_init(3),Energs_init(4)
      WRITE(5,*) 'End ',Energs_fin(1),Energs_fin(2),
     .   Energs_fin(3),Energs_fin(4)
      
      !Write out convergence error to combined file NormFile
      WRITE(NormFile,*) dx,Norms(1),Norms(2),Norms(3)
      
      !Write out conservation error to combined file EnergFile
      WRITE(EnergFile,*) dx,Energ_Err(1),Energ_Err(2),
     .   Energ_Err(3),Energ_Err(4)
      
      
      CLOSE(1)
      CLOSE(2)
      CLOSE(3)
      CLOSE(4)
      CLOSE(5)
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
     . tlist_len=31)
      
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
      !theta - reconstruction parameter-  generalised minmod
      !Cr - ensuring CFL satisfied
      !maxwavespeed - estimate for maximum wavespeed to satisfy CFL
      !beta1 - gSGN beta1 parameter
      !beta2 - gSGN beta2 parameter
      !alpha - used to bound wave speed when beta2 > beta1 (advancing dispersive waves)
      DOUBLE PRECISION a0,a1,ga,xstart,xend,tstart,tend,
     . dx,dt,theta,Cr,maxwavespeed,beta1,beta2,alpha
     
      !list of intermediate times to print solution at
      DOUBLE PRECISION tlist(tlist_len)
     
      !effective length of directory address string
      INTEGER effeclenwdir
      
      wdir = "../Results/Validation/" // 
     . "AnalyticSolutions/SerreSoliton/" 
     
      !use LenTrim routine to calculate effeclenwdir
      CALL LenTrim(wdir,wdirlen,effeclenwdir)
      
      !Refresh directories dump data
      CALL SYSTEM('rm -rf '//wdir)
      CALL SYSTEM('mkdir -p '// wdir)
      
      !open Conservation and Convergence output files
      OPEN(EnergFile, FILE = wdir(1:effeclenwdir)//'Energy.dat') 
      OPEN(NormFile, FILE = wdir(1:effeclenwdir)//'Norms.dat') 

      
      n_GhstCells = 6
      
      ga = 9.81d0
      
      beta1 = 2.0/3.0
      beta2 = 0d0
      
      a0 = 1.0
      a1 = 0.7
      
      xstart = -200d0
      xend = 200d0
      
      theta = 1.2d0
      
      tstart = 0d0
      tend = 1.0!10.0**(-14)!30d0
      
      ! tlist is every 1s between 0s and 30s (inclusive)
      CALL EqualSpaced(tstart,tlist_len,1.0, tlist)
      
      IF  (DABS(beta1) < 10d0**(-10))  THEN
         alpha = 1d0
      ELSE
         alpha = MAX(1d0,beta2 / beta1)
      END IF
      
      lowestresx = 100
      
      !Perform the Solitary Travelling Wave experiment a number of times, decreasing \Delta x each time
      DO expi = 0,14
      
         !make directory for inidividual experiment
         WRITE (strdiri,'(I2.2)') expi
         CALL SYSTEM('mkdir -p '// wdir(1:effeclenwdir) //strdiri)
         
         !calculate number of cells, both interior and total
         x_len = lowestresx *2**(expi)
         xbc_len = x_len + 2 *n_GhstCells

         !calculate dx to obtain the right number of cells
         !then calculate dt using CFL condition with Cr
         dx = (xend - xstart) / (x_len -1)
         Cr = 0.5
         maxwavespeed = dsqrt(alpha*ga*(a0 + a1))
         dt  = (Cr / maxwavespeed) *dx
         
         PRINT *,'++++++++++++ Experiment : ',expi ,' || ', '# Cells :',
     .    x_len , '++++++++++++'
         
         !Call routine that generates and then runs individual experiment
         CALL SingleSerreSoliton(a0,a1,ga,beta1,beta2,
     . xstart,xbc_len,n_GhstCells,dx,tstart,tend,
     . dt,theta,NormFile,EnergFile,tlist,tlist_len,
     .      wdir(1:effeclenwdir)//strdiri//'/',effeclenwdir+3)

      
      END DO
      
      CLOSE(EnergFile)
      CLOSE(NormFile)
      
      END
