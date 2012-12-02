      function f ( t )

c*******************************************************************************
c
cc F evaluates the power consumption function F(T).
c
c  Modified:
c
c    03 December 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Rubin Landau,
c    A First Course in Scientific Computing,
c    Princeton, 2005,
c    ISBN: 0-691-12183-4
c    LC: Q183.9.L36
c
c  Parameters:
c
c    Input, real T, the argument of the function.
c
c    Output, real F, the value of the function.
c
      implicit none

      real f
      real pi
      parameter ( pi = 3.14159265 )
      real t

      f = ( 4.0 + t / 365.0 + 0.5 * sin ( pi * t / 91.0 ) ) * 
     &    ( 2.0 + exp ( - sin ( 2.0 * pi * t ) ) )

      return
      end
