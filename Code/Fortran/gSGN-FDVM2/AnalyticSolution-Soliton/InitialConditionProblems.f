c ====================================================================================
c Collection of fortran subroutines to generate the initial condition problems
c Most of these routines also take a time input to generate the analytic solutions at a particular time
c ====================================================================================        

c ====
c Subroutine that generates the cell centres, which are our locations for the initial conditions h,G and u
c we generate both the interior and ghost cells. 
c
c Inputs :  x_start - first cell center in interior (excluding ghost cells)
c           dx - cell width
c           xbc_len - total number of cells (number of cells in interior [x_len] + 2 x number of ghost cells at boundary)
c           n_GhstCells - number of ghost cells at the boundary (will be 2 boundaries in 1D  (left and right))
c Outputs : xbc - cell centers from the discretisation of domain
c
c When xbc_len = x_len + 2*n_GhstCells (where x_len are the number of points in the inteior) and xbc_len is the number of interior
c points plus the ghost cells at each edge (hence the multiple of 2).
c When dx = (xend - xstart) / (x_len -1), the first cell centre in the interior at xbc(n_GhstCells + 1) will be xstart and the last 
c cell centre in the interior xbc(xbc_len - n_GhstCells) will be xend.
c ====      
      SUBROUTINE Generatexbc(x_start,dx,xbc_len,n_GhstCells,xbc)
          
      IMPLICIT NONE
      !Inputs
      DOUBLE PRECISION x_start,dx
      INTEGER xbc_len,n_GhstCells
      
      !Outputs
      DOUBLE PRECISION xbc(xbc_len)
      
      !Local variables
      INTEGER  i
      
      !Generate xbc, the cell centers in the interior and the ghost cells
      DO i=1,xbc_len  
         xbc(i)  = x_start + (i -1 - n_GhstCells )*dx
      END DO  
                    
      END

c ====
c Subroutine that generates the Serre Solitary Travelling Wave solution at the cell centers (equivalent to cell averages for second order methods)
c This is an analytic solution of the gSGN when beta1 = 2/3, beta2 = 0
c Inputs:xbc - array of cell centers of discretisation
c        xbc_len - length of xbc
c        t - time value to produce solution at
c        ga - acceleration due to gravity
c        a0 - background water depth on which solitary wave travels
c        a1 - amplitude of solitary wave
c
c Outputs:  h - array of h value at cell centers xbc (equivalent to cell average for second-order method)
c           u - array of u value at cell centers xbc
c           G - array of G value at cell centers xbc
c =====      
      SUBROUTINE SGNSolitaryTravellingWave(xbc,xbc_len,t,ga,a0,a1,
     . h,u,G)
      IMPLICIT NONE
      
      !Inputs
      INTEGER xbc_len
      DOUBLE PRECISION xbc(xbc_len)
      DOUBLE PRECISION t,a0,a1,ga
      
      !Outputs
      DOUBLE PRECISION h(xbc_len),G(xbc_len),u(xbc_len)
      
      !Local Variables
      DOUBLE PRECISION k,c,phi,sechkphi
      INTEGER i
      
      ! Derived parameters k (kappa in paper) and c
      k = DSQRT(3*a1) / (2*a0*DSQRT(a0 + a1))
      c = DSQRT(ga*(a0 + a1))
      
      DO i=1,xbc_len
         
         phi  = xbc(i) - c*t
         sechkphi = 1.0 / DCOSH(k*phi)
         h(i)  = a0 + a1*sechkphi**2
         u(i)  = c*(1 - a0 / h(i))
                        
         G(i) = u(i)*h(i) + 2.0/3*a0*a1*c*(k**2)*sechkphi**4*h(i)
     . - 4.0/3*a0*a1**2*c*k**2*sechkphi**4*DTANH(k*phi)**2 
     . - 4.0/3*a0*a1*c*k**2*sechkphi**2*h(i)*DTANH(k*phi)**2 
   
      END DO 
      
      END

c ====
c Subroutine that generates the SWWE Dam break solution at the cell centers (equivalent to cell averages for second order methods)
c This is an analytic solution of the gSGN when beta1 = 0, beta2 = 0
c Inputs:xbc - array of cell centers of discretisation
c        xbc_len - length of xbc
c        t - time value to produce solution at
c        ga - acceleration due to gravity
c        hl - the water depth left of x= 0
c        hr - the water depth right of x= 0
c
c Outputs:  h - array of h value at cell centers xbc (equivalent to cell average for second-order method)
c           u - array of u value at cell centers xbc
c           G - array of G value at cell centers xbc
c =====    
      SUBROUTINE Dambreak(xbc,xbc_len,t,ga,hl,hr,h,u,G)
      IMPLICIT NONE
      
      !Inputs
      INTEGER xbc_len
      DOUBLE PRECISION xbc(xbc_len)
      DOUBLE PRECISION t,hl,hr,ga
      
      !Outputs
      DOUBLE PRECISION h(xbc_len),G(xbc_len),u(xbc_len)
      
      !Local variables
      DOUBLE PRECISION h2,u2,S2,zmin,fzmin,z,fz,fzmax,zmax
      INTEGER i
      
      ! If t is small, we can just  use the simple initial condition problem to get h,u and G = uh (since beta1 = 0)
      IF (DABS(t) .LE. 10d0**(-10)) THEN
         DO i=1,xbc_len 

            IF (xbc(i) .LE. 0) THEN
               h(i) = hl
               u(i) = 0
            ELSE
               h(i) = hr
               u(i) = 0
            END IF
            
            G(i) = u(i)*h(i)
         END DO
      
      ! Otherwise we need to calculate analytic solution to Dambreak problem for SWWE
      ELSE
      
         !h2 - the height of middle state is solution of the following implicit equation.
         !h2 = hr/2*DSQRT(DSQRT(1 + 8*((2*h2/(h2 - hr))* ((DSQRT(ga*hl) - DSQRT(ga*h2))/DSQRT(ga*hr)) )**2 ) -1)
         !We use bisection method to calculate h2 (find zero of h2 - f(h2) = 0), since we only do this few times, we dont need to be more efficient.
         
         !Initialise
         !Left boundary value, and associated function value 
         zmin = hl
         fzmin = zmin - hr/2*(DSQRT(1 + 8*((2*zmin/(zmin - hr))*
     .      ((DSQRT(ga*hl) - DSQRT(ga*zmin))/ DSQRT(ga*hr)) )**2 ) -1)
         !right boundary value, and associated function value
         zmax = hr
         fzmax = zmax - hr/2*(DSQRT(1 + 8*((2*zmax/(zmax - hr))*
     .      ((DSQRT(ga*hl) - DSQRT(ga*zmax))/ DSQRT(ga*hr)) )**2 ) -1)
     
         !Loop bisection method, until left and right boundaries are smaller than the tolerance value 
         DO WHILE(DABS(zmax - zmin) .GT. 10d0**(-14))
           !bisection point and function value
           z = (zmin + zmax)/2.d0
           fz = z - hr/2*(DSQRT(1 + 8*((2*z/(z - hr))*
     .      ((DSQRT(ga*hl) - DSQRT(ga*z))/ sqrt(ga*hr)) )**2 ) -1)
        
           !Determine which half of interval to iterate on 
           IF(fz .GT. 0.d0)THEN
             zmin = z
             fzmin = fz
           ELSE
             zmax = z
             fzmax = fz
           END IF
           
         END DO
         
         ! After bisection, we have very good approximation to h2
         ! From which we can calculate shock speed  - S2
         ! and middle state velocity u2
         h2 = z
         S2 = 2d0*h2 / (h2 - hr) *(DSQRT(ga*hl) -  DSQRT(ga*h2))
         u2 = 2*(DSQRT(ga*hl) - DSQRT(ga*h2))
         
         ! Having h2,u2 and S2 we can give h,u,G for the dam break solution      
         DO i=1,xbc_len 
         
            IF (xbc(i) .LE. -t*DSQRT(ga*hl)) THEN
               h(i) = hl
               u(i) = 0
            ELSE IF ((xbc(i) .GT. -t*DSQRT(ga*hl)) 
     .         .AND. (xbc(i) .LE. t*(u2 - DSQRT(ga*h2)))) THEN
               h(i) = 4d0/(9d0*ga) *(DSQRT(ga*hl) - xbc(i)/(2d0*t))**2
               u(i) = 2d0/3d0*(DSQRT(ga*hl) + xbc(i)/t)
            ELSE IF ((xbc(i) .GT. t*(u2 - DSQRT(ga*h2)) ) 
     .         .AND. (xbc(i) .LT. t*S2)) THEN
               h(i) = h2    
               u(i) = u2  
            ELSE
               h(i) = hr
               u(i) = 0
            END IF
            
            G(i) = u(i)*h(i)
         
      
         END DO 
      
      END IF

      END
      



