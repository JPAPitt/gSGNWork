c ====================================================================================        
c Collection of Miscellaneous routines that are used
c ====================================================================================        

c LenTrim finds the length of string up to the first whitespace
c Inputs : String,n
c Output : resultval
c so a string 'abc  a' would produce  resultval = 3 
      SUBROUTINE LenTrim(String,n,resultval)
      !Inputs
      INTEGER n
      CHARACTER(len=n) String
      
      !Outputs
      INTEGER resultval

      resultval = 1

      ! Loop over String extracting charachter in resultval position
      ! If the charachter is not a space, keep iterating otherwise stop
      DO WHILE ((String(resultval:resultval) .NE. ' ' ) .AND.
     . (resultval .LE. n )  )
         resultval = resultval + 1
      END DO 
      
      resultval = resultval -1

      END 
      
      
c EqualSpaced produces an array of equally spaced locations/times
c Inputs : startv,num,stepsize
c Output : res
      SUBROUTINE EqualSpaced(startv,num,stepsize, res)
      !Inputs
      DOUBLE PRECISION startv,stepsize
      INTEGER num
      
      !Outputs
      DOUBLE PRECISION res(num)
      
      !Local variables
      INTEGER i
      
      ! iterate to produce res an array starting at startv and with spacing stepsize      
      DO  i = 1,num
         res(i) = startv + (i-1)*stepsize
      END DO 
      
      END
   
   
