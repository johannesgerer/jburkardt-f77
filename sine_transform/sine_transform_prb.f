      program main

c*********************************************************************72
c
cc SINE_TRANSFORM_TEST tests SINE_TRANSFORM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SINE_TRANSFORM_TEST'
      write ( *, '(a)' ) '  FORTRAN77 version.'
      write ( *, '(a)' ) '  Test the SINE_TRANSFORM library.'

      call sine_transform_test01 ( )
      call sine_transform_test02 ( )
      call sine_transform_test03 ( )
      call sine_transform_test04 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SINE_TRANSFORM_TEST'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine sine_transform_test01 ( )

c*********************************************************************72
c
cc SINE_TRANSFORM_TEST01 demonstrates that the transform is its own inverse.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer i
      integer seed
      double precision r(n)
      double precision s(n)
      double precision t(n)

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SINE_TRANSFORM_TEST01:'
      write ( *, '(a)' ) 
     &  '  SINE_TRANSFORM_DATA does a sine transform of data'
      write ( *, '(a)' ) '  defined by a vector.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Demonstrate that the transform is its own inverse.'
      write ( *, '(a)' ) '  Let R be a random N vector.'
      write ( *, '(a)' ) '  Let S be the transform of D.'
      write ( *, '(a)' ) '  Let T be the transform of E.'
      write ( *, '(a)' ) '  Then R and T will be equal.'

      call r8vec_uniform_01 ( n, seed, r )
      call sine_transform_data ( n, r, s )
      call sine_transform_data ( n, s, t )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I      R(I)        S(I)        T(I)'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,i4,2x,g10.4,2x,g10.4,2x,g10.4)' ) 
     &    i, r(i), s(i), t(i)
      end do

      return
      end
      subroutine sine_transform_test02 ( )

c*********************************************************************72
c
cc SINE_TRANSFORM_TEST02 uses the functional form of the routine.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 9 )

      double precision a
      double precision b
      double precision f1(n)
      double precision f2(n)
      double precision fa
      double precision fb
      integer i
      double precision poly5
      external poly5
      double precision s(n)
      double precision x(n)

      a = 1.0D+00
      b = 3.0D+00
c
c  Evenly spaced points between A and B, but omitting
c  A and B themselves.
c
      do i = 1, n
        x(i) = ( dble ( n - i + 1 ) * a   
     &         + dble (     i     ) * b ) 
     &         / dble ( n     + 1 )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SINE_TRANSFORM_TEST02:'
      write ( *, '(a)' ) 
     &  '  SINE_TRANSFORM_FUNCTION does a sine transform of data'
      write ( *, '(a)' ) 
     &  '  defined by a function F(X) evaluated at equally spaced'
      write ( *, '(a)' ) '  points in an interval [A,B].'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Demonstrate that the transform is its own inverse.'
      write ( *, '(a)' ) 
     &  '  Let X(0:N+1) be N+2 equally spaced points in [A,B].'
      write ( *, '(a)' ) '  Let S be the transform of F(X(1:N)).'
      write ( *, '(a)' ) 
     &  '  Let F1 be the linear interpolant of (A,F(A)), (B,F(B)).'
      write ( *, '(a)' ) '  Let F2 be the transform of S.'
      write ( *, '(a)' ) '  Then F(X(1:N)) = F1(X(1:N)) + F2(1:N).'

      call sine_transform_function ( n, a, b, poly5, s )

      fa = poly5 ( a )
      fb = poly5 ( b )
      do i = 1, n
        f1(i) = ( ( b - x(i)     ) * fa   
     &          + (     x(i) - a ) * fb ) 
     &          / ( b          - a )
      end do

      call sine_transform_data ( n, s, f2 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I      X(I)      F(X(I))       ' // 
     &  'S           F1          F2          F1+F2'
      write ( *, '(a)' ) ' '

      write ( *, 
     &  '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4,2x,' // 
     &  'f10.4,2x,f10.4,2x,f10.4)' ) 
     &  0, a, poly5 ( a ), 0.0D+00, fa, 0.0D+00, fa

      do i = 1, n
          write ( *, 
     &      '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4,2x,' // 
     &      'f10.4,2x,f10.4,2x,f10.4)' ) 
     &      i, x(i), poly5 ( x(i) ), s(i), f1(i), f2(i), 
     &      f1(i) + f2(i)
      end do

      write ( *, 
     &  '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4,2x,' // 
     &  'f10.4,2x,f10.4,2x,f10.4)' ) 
     &  n + 1, b, poly5 ( b ), 0.0D+00, fb, 0.0D+00, fb

      return
      end
      subroutine sine_transform_test03 ( )

c*********************************************************************72
c
cc SINE_TRANSFORM_TEST03 evaluates the sine transform interpolant.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 9 )
      integer n2
      parameter ( n2 = 1 + 2 * ( n + 1 ) )

      double precision a
      double precision b
      double precision f2(n2)
      double precision fa
      double precision fb
      integer i
      double precision poly5
      external poly5
      double precision s(n)
      double precision x(n)
      double precision x2(n2)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SINE_TRANSFORM_TEST03:'
      write ( *, '(a)' ) 
     &  '  SINE_TRANSFORM_FUNCTION does a sine transform of data'
      write ( *, '(a)' ) 
     &  '  defined by a function F(X) evaluated at N equally spaced'
      write ( *, '(a)' ) '  points in an interval [A,B].'
      write ( *, '(a)' ) 
     &  '  SINE_TRANSFORM_INTERPOLANT evaluates the interpolant.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  The interpolant will be 0 at the 0th and (N+1)-th points.'
      write ( *, '(a)' ) 
     &  '  It equals the function at points 1 through N.'
      write ( *, '(a)' ) 
     &  '  In between, it can approximate smooth functions,'
      write ( *, '(a)' ) '  and the approximation improves with N.'
c
c  N determines the number of data points, indexed by 1 to N.  
c  However, we essentially have N+2 data points, indexed 0 to N+1,
c  with the data value being 0 at the first and last auxilliary points.
c
      a = 1.0D+00
      b = 4.0D+00
c
c  Evenly spaced points between A and B, but omitting
c  A and B themselves.
c
      do i = 1, n
        x(i) = ( dble ( n - i + 1 ) * a   
     &         + dble (     i     ) * b ) 
     &         / dble ( n     + 1 )
      end do
c
c  Determine the interpolant coefficients.
c
      call sine_transform_function ( n, a, b, poly5, s )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I      X(I)      F(X(I))        S(I)'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4)' ) 
     &    i, x(i), poly5 ( x(i) ), s(i)
      end do
c
c  Evaluate the interpolant.
c
      fa = poly5 ( a )
      fb = poly5 ( b )
c
c  Evenly spaced points between A and B, including A and B,
c  and twice the density of the previous set of points.
c
      do i = 1, n2
        x2(i) = ( dble ( n2 - i     ) * a   
     &          + dble (      i - 1 ) * b ) 
     &          / dble ( n2     - 1 )
      end do

      call sine_transform_interpolant ( n, a, b, fa, fb, s, n2, x2, f2 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I      X            F(X)        FHAT(X)'
      write ( *, '(a)' ) ' '

      do i = 1, n2
        write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4)' ) 
     &    i, x2(i), poly5 ( x2(i) ), f2(i)
      end do

      return
      end
      subroutine sine_transform_test04 ( )

c*********************************************************************72
c
cc SINE_TRANSFORM_TEST04 evaluates the sine transform interpolant.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 15 )
      integer n2
      parameter ( n2 = 1 + 5 * ( n + 1 ) )

      double precision a
      double precision b
      double precision cosine_sum
      external cosine_sum
      double precision f2(n2)
      double precision fa
      double precision fb
      integer i
      double precision s(n)
      double precision x(n)
      double precision x2(n2)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SINE_TRANSFORM_TEST04:'
      write ( *, '(a)' ) 
     &  '  SINE_TRANSFORM_FUNCTION does a sine transform of data'
      write ( *, '(a)' ) 
     &  '  defined by a function F(X) evaluated at N equally spaced'
      write ( *, '(a)' ) '  points in an interval [A,B].'
      write ( *, '(a)' ) 
     &  '  SINE_TRANSFORM_INTERPOLANT evaluates the interpolant.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  The interpolant will be 0 at the 0th and (N+1)-th points.'
      write ( *, '(a)' ) 
     &  '  It equals the function at points 1 through N.'
      write ( *, '(a)' ) 
     &  '  In between, it can approximate smooth functions,'
      write ( *, '(a)' ) '  and the approximation improves with N.'
c
c  N determines the number of data points, indexed by 1 to N.  
c  However, we essentially have N+2 data points, indexed 0 to N+1,
c  with the data value being 0 at the first and last auxilliary points.
c
      a = 0.0D+00
      b = 7.0D+00
c
c  Evenly spaced points between A and B, but omitting
c  A and B themselves.
c
      do i = 1, n
        x(i) = ( dble ( n - i + 1 ) * a   
     &         + dble (     i     ) * b ) 
     &         / dble ( n     + 1 )
      end do
c
c  Determine the interpolant coefficients.
c
      call sine_transform_function ( n, a, b, cosine_sum, s )
c
c  Evaluate the interpolant.
c
      fa = cosine_sum ( a )
      fb = cosine_sum ( b )
c
c  Evenly spaced points between A and B, including A and B,
c  and twice the density of the previous set of points.
c
      do i = 1, n2
        x2(i) = ( dble ( n2 - i     ) * a   
     &          + dble (      i - 1 ) * b )
     &          / dble ( n2     - 1 )
      end do

      call sine_transform_interpolant ( n, a, b, fa, fb, s, n2, x2, f2 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Expect exact agreement every 5th sample.'
      write ( *, '(a)' ) ' '

      do i = 1, n2
        write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4)' ) 
     &    i, x2(i), cosine_sum ( x2(i) ), f2(i)
      end do

      return
      end
      function cosine_sum ( x )

c*********************************************************************72
c
cc COSINE_SUM evaluates a function which is a cosine sum.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision COSINE_SUM, the value.
c
      implicit none

      double precision cosine_sum
      double precision x

      cosine_sum =  cos (           x ) 
     &  + 5.0D+00 * cos ( 1.6D+00 * x ) 
     &  - 2.0D+00 * cos ( 2.0D+00 * x ) 
     &  + 5.0D+00 * cos ( 4.5D+00 * x ) 
     &  + 7.0D+00 * cos ( 9.0D+00 * x )

      return
      end
      function poly5 ( x )

c*********************************************************************72
c
cc POLY5 evaluates a particular fifth-degree polynomial.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision POLY5, the value of the polynomial at X.
c
      implicit none

      double precision poly5
      double precision x

      poly5 = ( x - 0.1D+00 ) * 
     &        ( x - 0.2D+00 ) * 
     &        ( x - 0.4D+00 ) * 
     &        ( x - 2.1D+00 ) * 
     &        ( x - 3.0D+00 )

      return
      end
