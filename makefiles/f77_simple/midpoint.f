      subroutine midpoint ( a, b, f, int_num, quad )

c*******************************************************************************
c
cc MIDPOINT approximates an integral using the composite midpoint rule.
c
c  Modified:
c
c    03 December 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real A, B, the endpoints of the interval of integration.
c
c    Input, real external F, the name of the function to be integrated.
c
c    Input, integer INT_NUM, the number of intervals to be used.
c
c    Output, real QUAD, the approximate value of the integral.
c
      implicit none

      real a
      real b
      real f
      external f
      integer i
      integer int_num
      real int_width
      real quad
      real t

      quad = 0.0

      do i = 1, int_num
    
        t = ( real ( 2 * int_num - 2 * i + 1 ) * a   
     &      + real (               2 * i - 1 ) * b ) 
     &      / real ( 2 * int_num )
  
        quad = quad + f ( t )

      end do

      int_width = ( b - a ) / real ( int_num )

      quad = quad * int_width

      return
      end
