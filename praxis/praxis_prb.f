      program main

c*********************************************************************72
c
cc MAIN is the main program for PRAXIS_PRB.
c
c  Discussion:
c
c    PRAXIS_PRB tests PRAXIS.
c
c  Modified:
c
c    02 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PRAXIS_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the PRAXIS library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
      call test05 ( )
      call test06 ( )
      call test07 ( )
      call test08 ( )
      call test09 ( )
      call test10 ( )
      call test11 ( )
      call test12 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PRAXIS_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 calls PRAXIS for the Beale function.
c
c  Modified:
c
c    02 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 2 )

      double precision f_01
      external f_01
      double precision fmin
      double precision h0
      integer i
      double precision machep
      double precision pr
      double precision praxis
      integer prin
      double precision r8_epsilon
      double precision t0
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  The Beale function.'

      t0 = 0.001D+00
      machep = r8_epsilon ( )
      h0 = 0.25D+00
      prin = 0

      x(1) = 0.1D+00
      x(2) = 0.1D+00

      call r8vec_print ( n, x, '  Initial point:' )

      write ( *, '(a,g14.6)' ) '  Function value = ', f_01 ( x, n )

      pr = praxis ( t0, machep, h0, n, prin, x, f_01, fmin )

      call r8vec_print ( n, x, '  Computed minimizer:' )
      
      write ( *, '(a,g14.6)' ) '  Function value = ', f_01 ( x, n )

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 calls PRAXIS for the Box function.
c
c  Modified:
c
c    02 February 2008
c
c  Author:
c
c    John Burkardt
c
c  Modified:
c
c    02 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      double precision f_02
      external f_02
      double precision fmin
      double precision h0
      integer i
      double precision machep
      double precision pr
      double precision praxis
      integer prin
      double precision r8_epsilon
      double precision t0
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  The Box function.'

      t0 = 0.001D+00
      machep = r8_epsilon ( )
      h0 = 20.0D+00
      prin = 0

      x(1) = 0.0D+00
      x(2) = 10.0D+00
      x(3) = 20.0D+00

      call r8vec_print ( n, x, '  Initial point:' )

      write ( *, '(a,g14.6)' ) '  Function value = ', f_02 ( x, n )

      pr = praxis ( t0, machep, h0, n, prin, x, f_02, fmin )

      call r8vec_print ( n, x, '  Computed minimizer:' )
      
      write ( *, '(a,g14.6)' ) '  Function value = ', f_02 ( x, n )

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 calls PRAXIS for the Chebyquad function.
c
c  Modified:
c
c    02 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 8 )

      double precision f_03
      external f_03
      double precision fmin
      double precision h0
      integer i
      double precision machep
      double precision pr
      double precision praxis
      integer prin
      double precision r8_epsilon
      double precision t0
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  The Chebyquad function.'

      t0 = 0.001D+00
      machep = r8_epsilon ( )
      h0 = 0.1D+00
      prin = 0

      do i = 1, n
        x(i) = dble ( i ) / dble ( n + 1 )
      end do

      call r8vec_print ( n, x, '  Initial point:' )

      write ( *, '(a,g14.6)' ) '  Function value = ', f_03 ( x, n )

      pr = praxis ( t0, machep, h0, n, prin, x, f_03, fmin )

      call r8vec_print ( n, x, '  Computed minimizer:' )
      
      write ( *, '(a,g14.6)' ) '  Function value = ', f_03 ( x, n )

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 calls PRAXIS for the Cube function.
c
c  Modified:
c
c    02 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 2 )

      double precision f_04
      external f_04
      double precision fmin
      double precision h0
      integer i
      double precision machep
      double precision pr
      double precision praxis
      integer prin
      double precision r8_epsilon
      double precision t0
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) '  The Cube function.'

      t0 = 0.001D+00
      machep = r8_epsilon ( )
      h0 = 1.0D+00
      prin = 0

      x(1) = -1.2D+00
      x(2) = -1.0D+00

      call r8vec_print ( n, x, '  Initial point:' )

      write ( *, '(a,g14.6)' ) '  Function value = ', f_04 ( x, n )

      pr = praxis ( t0, machep, h0, n, prin, x, f_04, fmin )

      call r8vec_print ( n, x, '  Computed minimizer:' )
      
      write ( *, '(a,g14.6)' ) '  Function value = ', f_04 ( x, n )

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 calls PRAXIS for the Fletcher-Powell Helix function.
c
c  Modified:
c
c    02 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      double precision f_05
      external f_05
      double precision fmin
      double precision h0
      integer i
      double precision machep
      double precision pr
      double precision praxis
      integer prin
      double precision r8_epsilon
      double precision t0
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) '  The Fletcher-Powell Helix function.'

      t0 = 0.001D+00
      machep = r8_epsilon ( )
      h0 = 1.0D+00
      prin = 0

      x(1) = -1.0D+00
      x(2) =  0.0D+00
      x(3) =  0.0D+00

      call r8vec_print ( n, x, '  Initial point:' )

      write ( *, '(a,g14.6)' ) '  Function value = ', f_05 ( x, n )

      pr = praxis ( t0, machep, h0, n, prin, x, f_05, fmin )

      call r8vec_print ( n, x, '  Computed minimizer:' )
      
      write ( *, '(a,g14.6)' ) '  Function value = ', f_05 ( x, n )

      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 calls PRAXIS for the Hilbert function.
c
c  Modified:
c
c    02 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      double precision f_06
      external f_06
      double precision fmin
      double precision h0
      integer i
      double precision machep
      double precision pr
      double precision praxis
      integer prin
      double precision r8_epsilon
      double precision t0
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' ) '  The Hilbert function.'

      t0 = 0.001D+00
      machep = r8_epsilon ( )
      h0 = 10.0D+00
      prin = 0

      do i = 1, n
        x(i) = 1.0D+00
      end do

      call r8vec_print ( n, x, '  Initial point:' )

      write ( *, '(a,g14.6)' ) '  Function value = ', f_06 ( x, n )

      pr = praxis ( t0, machep, h0, n, prin, x, f_06, fmin )

      call r8vec_print ( n, x, '  Computed minimizer:' )
      
      write ( *, '(a,g14.6)' ) '  Function value = ', f_06 ( x, n )

      return
      end
      subroutine test07 ( )

c*********************************************************************72
c
cc TEST07 calls PRAXIS for the Powell 3D function.
c
c  Modified:
c
c    02 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      double precision f_07
      external f_07
      double precision fmin
      double precision h0
      integer i
      double precision machep
      double precision pr
      double precision praxis
      integer prin
      double precision r8_epsilon
      double precision t0
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07'
      write ( *, '(a)' ) '  The Powell 3D function.'

      t0 = 0.001D+00
      machep = r8_epsilon ( )
      h0 = 1.0D+00
      prin = 0

      x(1) = 0.0D+00
      x(2) = 1.0D+00
      x(3) = 2.0D+00

      call r8vec_print ( n, x, '  Initial point:' )

      write ( *, '(a,g14.6)' ) '  Function value = ', f_07 ( x, n )

      pr = praxis ( t0, machep, h0, n, prin, x, f_07, fmin )

      call r8vec_print ( n, x, '  Computed minimizer:' )
      
      write ( *, '(a,g14.6)' ) '  Function value = ', f_07 ( x, n )


      return
      end
      subroutine test08 ( )

c*********************************************************************72
c
cc TEST08 calls PRAXIS for the Rosenbrock function.
c
c  Modified:
c
c    02 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 2 )

      double precision f_08
      external f_08
      double precision fmin
      double precision h0
      integer i
      double precision machep
      double precision pr
      double precision praxis
      integer prin
      double precision r8_epsilon
      double precision t0
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08'
      write ( *, '(a)' ) '  The Rosenbrock function.'

      t0 = 0.001D+00
      machep = r8_epsilon ( )
      h0 = 1.0D+00
      prin = 0

      x(1) = -1.2D+00
      x(2) =  1.0D+00

      call r8vec_print ( n, x, '  Initial point:' )

      write ( *, '(a,g14.6)' ) '  Function value = ', f_08 ( x, n )

      pr = praxis ( t0, machep, h0, n, prin, x, f_08, fmin )

      call r8vec_print ( n, x, '  Computed minimizer:' )
      
      write ( *, '(a,g14.6)' ) '  Function value = ', f_08 ( x, n )

      return
      end
      subroutine test09 ( )

c*********************************************************************72
c
cc TEST09 calls PRAXIS for the Powell Singular function.
c
c  Modified:
c
c    02 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      double precision f_09
      external f_09
      double precision fmin
      double precision h0
      integer i
      double precision machep
      double precision pr
      double precision praxis
      integer prin
      double precision r8_epsilon
      double precision t0
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST09'
      write ( *, '(a)' ) '  The Powell Singular function.'

      t0 = 0.001D+00
      machep = r8_epsilon ( )
      h0 = 1.0D+00
      prin = 0

      x(1) =  3.0D+00
      x(2) = -1.0D+00
      x(3) =  0.0D+00
      x(3) =  1.0D+00

      call r8vec_print ( n, x, '  Initial point:' )

      write ( *, '(a,g14.6)' ) '  Function value = ', f_09 ( x, n )

      pr = praxis ( t0, machep, h0, n, prin, x, f_09, fmin )

      call r8vec_print ( n, x, '  Computed minimizer:' )
      
      write ( *, '(a,g14.6)' ) '  Function value = ', f_09 ( x, n )

      return
      end
      subroutine test10 ( )

c*********************************************************************72
c
cc TEST10 calls PRAXIS for the Tridiagonal function.
c
c  Modified:
c
c    02 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      double precision f_10
      external f_10
      double precision fmin
      double precision h0
      integer i
      double precision machep
      double precision pr
      double precision praxis
      integer prin
      double precision r8_epsilon
      double precision t0
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST10'
      write ( *, '(a)' ) '  The Tridiagonal function.'

      t0 = 0.001D+00
      machep = r8_epsilon ( )
      h0 = 8.0D+00
      prin = 0

      do i = 1, n
        x(i) = 0.0D+00
      end do

      call r8vec_print ( n, x, '  Initial point:' )

      write ( *, '(a,g14.6)' ) '  Function value = ', f_10 ( x, n )

      pr = praxis ( t0, machep, h0, n, prin, x, f_10, fmin )

      call r8vec_print ( n, x, '  Computed minimizer:' )
      
      write ( *, '(a,g14.6)' ) '  Function value = ', f_10 ( x, n )

      return
      end
      subroutine test11 ( )

c*********************************************************************72
c
cc TEST11 calls PRAXIS for the Watson function.
c
c  Modified:
c
c    02 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 6 )

      double precision f_11
      external f_11
      double precision fmin
      double precision h0
      integer i
      double precision machep
      double precision pr
      double precision praxis
      integer prin
      double precision r8_epsilon
      double precision t0
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST11'
      write ( *, '(a)' ) '  The Watson function.'

      t0 = 0.001D+00
      machep = r8_epsilon ( )
      h0 = 1.0D+00
      prin = 0

      do i = 1, n
        x(i) = 0.0D+00
      end do

      call r8vec_print ( n, x, '  Initial point:' )

      write ( *, '(a,g14.6)' ) '  Function value = ', f_11 ( x, n )

      pr = praxis ( t0, machep, h0, n, prin, x, f_11, fmin )

      call r8vec_print ( n, x, '  Computed minimizer:' )
      
      write ( *, '(a,g14.6)' ) '  Function value = ', f_11 ( x, n )

      return
      end
      subroutine test12 ( )

c*********************************************************************72
c
cc TEST12 calls PRAXIS for the Wood function.
c
c  Modified:
c
c    02 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      double precision f_12
      external f_12
      double precision fmin
      double precision h0
      integer i
      double precision machep
      double precision pr
      double precision praxis
      integer prin
      double precision r8_epsilon
      double precision t0
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST12'
      write ( *, '(a)' ) '  The Wood function.'

      t0 = 0.001D+00
      machep = r8_epsilon ( )
      h0 = 10.0D+00
      prin = 0

      x(1) = -3.0D+00
      x(2) = -1.0D+00
      x(3) = -3.0D+00
      x(4) = -1.0D+00

      call r8vec_print ( n, x, '  Initial point:' )

      write ( *, '(a,g14.6)' ) '  Function value = ', f_12 ( x, n )

      pr = praxis ( t0, machep, h0, n, prin, x, f_12, fmin )

      call r8vec_print ( n, x, '  Computed minimizer:' )
      
      write ( *, '(a,g14.6)' ) '  Function value = ', f_12 ( x, n )

      return
      end
      function f_01 ( x, n )

c*********************************************************************72
c
cc F_01 evaluates the Beale function.
c
c  Discussion:
c
c    The function is the sum of the squares of three functions.
c
c    This function has a valley approaching the line X(2) = 1.
c
c  Reference:
c
c    E Beale,
c    On an Iterative Method for Finding a Local Minimum of a Function
c    of More than One Variable,
c    Technical Report 25, Statistical Techniques Research Group,
c    Princeton University, 1958.
c
c    Richard Brent,
c    Algorithms for Finding Zeros and Extrema of Functions Without
c    Calculating Derivatives,
c    Stanford University Technical Report STAN-CS-71-198.
c
c  Modified:
c
c    04 May 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X(N), the argument of the objection function.
c
c    Input, integer N, the number of variables.
c
c    Output, double precision F_01, the value of the objective function.
c
      implicit none

      integer n

      double precision c1
      parameter ( c1 = 1.5D+00 )
      double precision c2
      parameter ( c2 = 2.25D+00 )
      double precision c3
      parameter ( c3 = 2.625D+00 )
      double precision f_01
      double precision x(n)

      f_01 =    ( c1 - x(1) * ( 1.0D+00 - x(2)    ) )**2 
     &        + ( c2 - x(1) * ( 1.0D+00 - x(2)**2 ) )**2 
     &        + ( c3 - x(1) * ( 1.0D+00 - x(2)**3 ) )**2

      return
      end
      function f_02 ( x, n )

c*********************************************************************72
c
cc F_02 evaluates the Box function.
c
c  Discussion:
c
c    The function is formed by the sum of squares of 10 separate terms.
c
c  Modified:
c
c    04 May 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X(N), the argument of the objection function.
c
c    Input, integer N, the number of variables.
c
c    Output, double precision F_02, the value of the objective function.
c
      implicit none

      integer n

      double precision c
      double precision f_02
      double precision fi
      integer i
      double precision x(n)

      f_02 = 0.0D+00

      do i = 1, 10

        c = - dble ( i ) / 10.0D+00

        fi = exp ( c * x(1) ) - exp ( c * x(2) ) - x(3) * 
     &    ( exp ( c ) - exp ( 10.0D+00 * c ) )
       
        f_02 = f_02 + fi**2

      end do

      return
      end
      function f_03 ( x, n )

c*********************************************************************72
c
cc F_03 evaluates the Chebyquad function.
c
c  Discussion:
c
c    The function is formed by the sum of squares of N separate terms.
c
c  Modified:
c
c    26 February 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X(N), the argument of the objection function.
c
c    Input, integer N, the number of variables.
c
c    Output, double precision F_03, the value of the function.
c
      implicit none

      integer n

      double precision f_03
      double precision fvec(n)
      integer i
      integer j
      double precision t
      double precision t1
      double precision t2
      double precision th
      double precision x(n)

      fvec(1:n) = 0.0D+00

      do j = 1, n
        t1 = 1.0D+00
        t2 = 2.0D+00 * x(j) - 1.0D+00
        t = 2.0D+00 * t2
        do i = 1, n
          fvec(i) = fvec(i) + t2
          th = t * t2 - t1
          t1 = t2
          t2 = th
        end do
      end do

      do i = 1, n
        fvec(i) = fvec(i) / dble ( n )
        if ( mod ( i, 2 ) == 0 ) then
          fvec(i) = fvec(i) + 1.0D+00 / ( dble ( i )**2 - 1.0D+00 )
        end if
      end do
c
c  Compute F.
c
      f_03 = sum ( fvec(1:n)**2 )

      return
      end
      function f_04 ( x, n )

c*********************************************************************72
c
cc F_04 evaluates the Cube function.
c
c  Discussion:
c
c    The function is the sum of the squares of two functions.
c
c  Modified:
c
c    27 February 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X(N), the argument of the objection function.
c
c    Input, integer N, the number of variables.
c
c    Output, double precision F_04, the value of the objective function.
c
      implicit none

      integer n

      double precision f_04
      double precision x(n)

      f_04 = ( 10.0D+00 * ( x(2) - x(1)**3 ) )**2 
     &  + ( 1.0D+00 - x(1) )**2

      return
      end
      function f_05 ( x, n )

c*********************************************************************72
c
cc F_05 evaluates the Helix function.
c
c  Discussion:
c
c    The function is the sum of the squares of three functions.
c
c  Modified:
c
c    27 February 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X(N), the argument of the objection function.
c
c    Input, integer N, the number of variables.
c
c    Output, double precision F_05, the value of the objective function.
c
      implicit none

      integer n

      double precision f_05
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision theta
      double precision x(n)

      r = sqrt ( x(1)**2 + x(2)**2 )

      if ( x(1) >= 0.0D+00 ) then
        theta = 0.5D+00 * atan2 ( x(2), x(1) ) / pi
      else if ( x(1) < 0.0D+00 ) then
        theta = 0.5D+00 * ( atan2 ( x(2), x(1) ) + pi ) / pi
      end if

      f_05 = 
     &    ( 10.0D+00 * ( x(3) - 10.0D+00 * theta ) )**2 
     &  + ( 10.0D+00 * ( r - 1.0D+00 ) )**2 
     &  + x(3)**2

      return
      end
      function f_06 ( x, n )

c*********************************************************************72
c
cc F_06 evaluates the Hilbert function.
c
c  Discussion:
c
c    The function is a positive definite quadratic function of 
c    the form
c
c      f(x) = x' A x
c
c    where A is the Hilbert matrix.
c
c  Modified:
c
c    27 February 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X(N), the argument of the objection function.
c
c    Input, integer N, the number of variables.
c
c    Output, double precision F_06, the value of the objective function.
c
      implicit none

      integer n

      double precision f_06
      integer i
      integer j
      double precision x(n)

      f_06 = 0.0D+00

      do i = 1, n
        do j = 1, n
          f_06 = f_06 + x(i) * x(j) / dble ( i + j - 1 )
        end do
      end do

      return
      end
      function f_07 ( x, n )

c*********************************************************************72
c
cc F_07 evaluates the Powell 3D function.
c
c  Reference:
c
c    M J D Powell,
c    An Efficient Method for Finding the Minimum of a Function of
c      Several Variables Without Calculating Derivatives,
c    Computer Journal, 
c    Volume 7, Number 2, pages 155-162, 1964.
c    
c  Modified:
c
c    03 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X(N), the argument of the objection function.
c
c    Input, integer N, the number of variables.
c
c    Output, double precision F_07, the value of the objective function.
c
      implicit none

      integer n

      double precision f_07
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision x(n)

      f_07 = 3.0D+00 - 1.0D+00 / ( 1.0D+00 + ( x(1) - x(2) )**2 ) 
     &  - sin ( 0.5D+00 * pi * x(2) * x(3) ) 
     &  - exp ( - ( ( x(1) - 2.0D+00 * x(2) + x(3) ) / x(2) )**2 )

      return
      end
      function f_08 ( x, n )

c*********************************************************************72
c
cc F_08 evaluates the Rosenbrock function.
c
c  Modified:
c
c    01 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X(N), the argument of the objection function.
c
c    Input, integer N, the number of variables.
c
c    Output, double precision F_08, the value of the objective function.
c
      implicit none

      integer n

      double precision f_08
      integer j
      double precision x(n)

      f_08 = 0.0D+00
      do j = 1, n
        if ( mod ( j, 2 ) == 1 ) then
          f_08 = f_08 + ( 1.0D+00 - x(j) )**2
        else
          f_08 = f_08 + 100.0D+00 * ( x(j) - x(j-1)**2 )**2
        end if
      end do

      return
      end
      function f_09 ( x, n )

c*********************************************************************72
c
cc F_09 evaluates the Powell Singular function.
c
c  Modified:
c
c    01 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X(N), the argument of the objection function.
c
c    Input, integer N, the number of variables.
c
c    Output, double precision F_09, the value of the objective function.
c
      implicit none

      integer n

      double precision f_09
      double precision f1
      double precision f2
      double precision f3
      double precision f4
      integer j
      double precision x(n)
      double precision xjp1
      double precision xjp2
      double precision xjp3

      f_09 = 0.0D+00

      do j = 1, n, 4

        if ( j + 1 <= n ) then
          xjp1 = x(j+1)
        else
          xjp1 = 0.0D+00
        end if

        if ( j + 2 <= n ) then
          xjp2 = x(j+2)
        else
          xjp2 = 0.0D+00
        end if

        if ( j + 3 <= n ) then
          xjp3 = x(j+3)
        else
          xjp3 = 0.0D+00
        end if
     
        f1 = x(j) + 10.0D+00 * xjp1
        if ( j + 1 <= n ) then
          f2 = xjp2 - xjp3
        else
          f2 = 0.0D+00
        end if
        if ( j + 2 <= n ) then
          f3 = xjp1 - 2.0D+00 * xjp2
        else
          f3 = 0.0D+00
        end if
        if ( j + 3 <= n ) then
          f4 = x(j) - xjp3
        else
          f4 = 0.0D+00
        end if

        f_09 = f_09 + f1**2 + 5.0D+00 * f2**2 
     &  + f3**4 + 10.0D+00 * f4**4

      end do

      return
      end
      function f_10 ( x, n )

c*********************************************************************72
c
cc F_10 evaluates the tridiagonal function.
c
c  Modified:
c
c    01 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X(N), the argument of the objection function.
c
c    Input, integer N, the number of variables.
c
c    Output, double precision F_10, the value of the objective function.
c
      implicit none

      integer n

      double precision f_10
      integer i
      double precision x(n)

      f_10 = x(1)**2 + 2.0D+00 * sum ( x(2:n)**2 )

      do i = 1, n-1
        f_10 = f_10 - 2.0D+00 * x(i) * x(i+1)
      end do

      f_10 = f_10 - 2.0D+00 * x(1)

      return
      end
      function f_11 ( x, n )

c*********************************************************************72
c
cc F_11 evaluates the Watson function.
c
c  Modified:
c
c    01 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X(N), the argument of the objection function.
c
c    Input, integer N, the number of variables.
c
c    Output, double precision F_11, the value of the objective function.
c
      implicit none

      integer n

      double precision d
      double precision f_11
      integer i
      integer j
      double precision s1
      double precision s2
      double precision x(n)

      f_11 = 0.0D+00
      do i = 1, 29

        s1 = 0.0D+00
        d = 1.0D+00
        do j = 2, n
          s1 = s1 + dble ( j - 1 ) * d * x(j)
          d = d * dble ( i ) / 29.0D+00
        end do

        s2 = 0.0D+00
        d = 1.0D+00
        do j = 1, n
          s2 = s2 + d * x(j)
          d = d * dble ( i ) / 29.0D+00
        end do

        f_11 = f_11 + ( s1 - s2**2 - 1.0D+00 )**2

      end do

      f_11 = f_11 + x(1)**2 + ( x(2) - x(1)**2 - 1.0D+00 )**2

      return
      end
      function f_12 ( x, n )

c*********************************************************************72
c
cc F_12 evaluates the Wood function.
c
c  Modified:
c
c    01 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X(N), the argument of the objection function.
c
c    Input, integer N, the number of variables.
c
c    Output, double precision F_12, the value of the objective function.
c
      implicit none

      integer n

      double precision f_12
      double precision f1
      double precision f2
      double precision f3
      double precision f4
      double precision f5
      double precision f6
      double precision x(n)

      f1 = x(2) - x(1)**2
      f2 = 1.0D+00 - x(1)
      f3 = x(4) - x(3)**2
      f4 = 1.0D+00 - x(3)
      f5 = x(2) + x(4) - 2.0D+00
      f6 = x(2) - x(4)

      f_12 = 100.0D+00 * f1**2 + f2**2 + 90.0D+00 * f3**2 + f4**2 
     &  + 10.0D+00 * f5**2 + 0.1D+00 * f6**2

      return
      end
      function r8_epsilon ( )

c*********************************************************************72
c
cc R8_EPSILON returns the R8 roundoff unit.
c
c  Discussion:
c
c    The roundoff unit is a number R which is a power of 2 with the
c    property that, to the precision of the computer's arithmetic,
c      1 .lt. 1 + R
c    but
c      1 = ( 1 + R / 2 )
c
c    FORTRAN90 provides the superior library routine
c
c      EPSILON ( X )
c
c  Modified:
c
c    06 March 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision R8_EPSILON, the R8 roundoff unit.
c
      implicit none

      double precision r8
      double precision r8_epsilon
      double precision r8_test

      r8 = 1.0D+00
      r8_test = 1.0D+00 + ( r8 / 2.0D+00 )

10    continue

      if ( 1.0D+00 .lt. r8_test ) then
        r8 = r8 / 2.0D+00
        r8_test = 1.0D+00 + ( r8 / 2.0D+00 )
        go to 10
      end if

      r8_epsilon = r8

      return
      end
      subroutine r8vec_print ( n, a, title )

c*********************************************************************72
c
cc R8VEC_PRINT prints an R8VEC.
c
c  Discussion:
c
c    An R8VEC is an array of double precision real values.
c
c  Modified:
c
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of components of the vector.
c
c    Input, double precision A(N), the vector to be printed.
c
c    Input, character * ( * ) TITLE, an optional title.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      integer s_len_trim
      character ( len = * ) title
      integer title_length

      title_length = s_len_trim ( title )
      if ( 0 .lt. title_length ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) title(1:title_length)
      end if

      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,2x,g16.8)' ) i, a(i)
      end do

      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Discussion:
c
c    This FORTRAN77 version is made available for cases where the
c    FORTRAN90 version cannot be used.
c
c  Modified:
c
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    None
c
      implicit none

      character * ( 8 ) ampm
      integer d
      character * ( 8 ) date
      integer h
      integer m
      integer mm
      character * ( 9 ) month(12)
      integer n
      integer s
      character * ( 10 ) time
      integer y

      save month

      data month /
     &  'January  ', 'February ', 'March    ', 'April    ',
     &  'May      ', 'June     ', 'July     ', 'August   ',
     &  'September', 'October  ', 'November ', 'December ' /

      call date_and_time ( date, time )

      read ( date, '(i4,i2,i2)' ) y, m, d
      read ( time, '(i2,i2,i2,1x,i3)' ) h, n, s, mm

      if ( h .lt. 12 ) then
        ampm = 'AM'
      else if ( h .eq. 12 ) then
        if ( n .eq. 0 .and. s .eq. 0 ) then
          ampm = 'Noon'
        else
          ampm = 'PM'
        end if
      else
        h = h - 12
        if ( h .lt. 12 ) then
          ampm = 'PM'
        else if ( h .eq. 12 ) then
          if ( n .eq. 0 .and. s .eq. 0 ) then
            ampm = 'Midnight'
          else
            ampm = 'AM'
          end if
        end if
      end if

      write ( *,
     &  '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' )
     &  d, month(m), y, h, ':', n, ':', s, '.', mm, ampm

      return
      end
      function s_len_trim ( s )

c*********************************************************************72
c
cc S_LEN_TRIM returns the length of a string to the last nonblank.
c
c  Modified:
c
c    05 March 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*(*) S, a string.
c
c    Output, integer S_LEN_TRIM, the length of the string to the last nonblank.
c
      implicit none

      integer i
      character*(*) s
      integer s_len_trim

      do i = len ( s ), 1, -1

        if ( s(i:i) .ne. ' ' ) then
          s_len_trim = i
          return
        end if

      end do

      s_len_trim = 0

      return
      end
