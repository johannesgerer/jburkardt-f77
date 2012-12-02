      program root

c*********************************************************************72
c
c  This program estimates the square root of a number.  The user 
c  types in a value, and ten steps of Newton's method are used to 
c  estimate the square root.  Negative input values are not 
c  acceptable.

      integer i
      real x
      real xroot

c  Get the value of the number whose square root is desired.

10    continue

      write (*,*) 'Enter a number for which you want the square root.'
      read (*,*) x

c  Watch out for negative or zero values!

      if (x .lt. 0.0) then

        write (*,*) 'Your input was negative, and unacceptable!'
        go to 10

      else if (x.eq.0.0) then

        write (*,*) 'The answer is 0.0'
        stop

      endif

c  Use 10 steps of Newton's method

      xroot = x/2.0

      do 20 i=1,10

        write (*,*) 'Current estimate is ',xroot
        xroot = (xroot + (x/xroot) ) * 0.5

20    continue

c  Print the result.

      write (*,*) 'The square root of ', x,' is ', xroot
      write (*,*) 'We check by printing xroot**2=', xroot**2

      stop
      end
