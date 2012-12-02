      program main

c*********************************************************************72
c
cc MINPACK_PRB runs the MINPACK tests.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MINPACK_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the MINPACK library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
      call test05 ( )
      call test06 ( )
      call test07 ( )
      call test08 ( )
      call test09 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MINPACK_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests CHKDER.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )
      integer m
      parameter ( m = n )
      integer ldfjac
      parameter ( ldfjac = n )

      double precision err(m)
      double precision fjac(ldfjac,n)
      double precision fvec(m)
      double precision fvecp(m)
      integer i
      integer ido
      integer iflag
      integer j
      integer mode
      double precision x(n)
      double precision xp(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  CHKDER compares a user supplied jacobian'
      write ( *, '(a)' ) '  and a finite difference approximation to it'
      write ( *, '(a)' ) '  and judges whether the jacobian is correct.'

      do ido = 1, 2

        if ( ido .eq. 1 ) then

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  On test 1, use a correct jacobian.'

        else if ( ido .eq. 2 ) then

           write ( *, '(a)' ) ' '
           write ( *, '(a)' ) '  On test 2, use a "bad" jacobian'
           write ( *, '(a)' ) '  and see if the routine notices!'

         end if
c
c  Set the point at which the test is to be made:
c
        do i = 1, n
          x(i) = 0.5D+00
        end do

        call r8vec_print ( n, x, '  Evaluation point X:' )

        mode = 1
        call chkder ( m, n, x, fvec, fjac, ldfjac, xp, fvecp, mode,
     &    err )

        iflag = 1

        call f01 ( n, x, fvec, fjac, ldfjac, iflag )
        call f01 ( n, xp, fvecp, fjac, ldfjac, iflag )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Sampled function values F(X) and F(XP)'
        write ( *, '(a)' ) ' '
        do i = 1, m
          write ( *, '(i3,2g14.6)' ) i, fvec(i), fvecp(i)
        end do

        iflag = 2
        call f01 ( n, x, fvec, fjac, ldfjac, iflag )
c
c  Here's where we put a mistake into the jacobian, on purpose.
c
        if ( ido == 2 ) then
          fjac(1,1) = 1.01D+00 * fjac(1,1)
          fjac(2,3) = - fjac(2,3)
        end if

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Computed jacobian'
        write ( *, '(a)' ) ' '
        do i = 1, m
          write ( *, '(5g14.6)' ) fjac(i,1:n)
        end do

        mode = 2
        call chkder ( m, n, x, fvec, fjac, ldfjac, xp, fvecp, mode,
     &  err )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  CHKDER gradient error estimates:'
        write ( *, '(a)' ) '     > 0.5, probably correct.'
        write ( *, '(a)' ) '     < 0.5, probably incorrect.'
        write ( *, '(a)' ) ' '
        do i = 1, m
          write ( *, '(i6,g14.6)' ) i, err(i)
        end do

      end do

      return
      end
      subroutine f01 ( n, x, fvec, fjac, ldfjac, iflag )

c*********************************************************************72
c
cc F01 is a function/jacobian routine.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 April 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of variables.
c
c    Input, double precision X(N), the variable values.
c
c    Output, double precision FVEC(N), the function values at X,
c    if IFLAG = 1.
c
c    Output, double precision FJAC(LDFJAC,N), the N by N jacobian
c    at X, if IFLAG = 2.
c
c    Input, integer LDFJAC, the leading dimension of FJAC, which must
c    be at least N.
c
c    Input, integer IFLAG:
c    1, please compute F(I) (X).
c    2, please compute FJAC(I,J) (X).
c
      implicit none

      integer ldfjac
      integer n

      double precision fjac(ldfjac,n)
      double precision fvec(n)
      integer i
      integer iflag
      integer j
      double precision x(n)
      double precision x_prod
      double precision x_sum

      x_prod = 1.0D+00
      x_sum = 0.0D+00
      do i = 1, n
        x_prod = x_prod * x(i)
        x_sum = x_sum + x(i)
      end do
c
c  If IFLAG is 1, we are supposed to evaluate F(X).
c
      if ( iflag .eq. 1 ) then

        do i = 1, n - 1
          fvec(i) = x(i) - dble ( n + 1 ) + x_sum
        end do

        fvec(n) = x_prod - 1.0D+00
c
c  If IFLAG is 2, we are supposed to evaluate FJAC(I,J) = d F(I)/d X(J)
c
      else if ( iflag .eq. 2 ) then

        do j = 1, n
          do i = 1, n - 1
            fjac(i,j) = 1.0D+00
          end do
        end do

        do i = 1, n - 1
          fjac(i,i) = 2.0D+00
        end do

        do j = 1, n
          fjac(n,j) = x_prod / x(j)
        end do

      end if

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests HYBRD1.
c
c  Discussion:
c
c    This is an example of what your main program would look
c    like if you wanted to use MINPACK to solve N nonlinear equations
c    in N unknowns.  In this version, we avoid computing the jacobian
c    matrix, and request that MINPACK approximate it for us.
c
c    The set of nonlinear equations is:
c
c      x1 * x1 - 10 * x1 + x2 * x2 + 8 = 0
c      x1 * x2 * x2 + x1 - 10 * x2 + 8 = 0
c
c    with solution x1 = x2 = 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 2 )
      integer lwa
      parameter ( lwa = ( n * ( 3 * n + 13 ) ) / 2 )

      external f02
      double precision fvec(n)
      integer i
      integer iflag
      integer info
      double precision tol
      double precision wa(lwa)
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' )
     &  '  HYBRD1 solves a nonlinear system of equations.'

      x(1) = 3.0D+00
      x(2) = 0.0D+00
      call r8vec_print ( n, x, '  Initial X:' )
      iflag = 1
      call f02 ( n, x, fvec, iflag )
      call r8vec_print ( n, fvec, '  F(X):' )

      tol = 0.00001D+00

      call hybrd1 ( f02, n, x, fvec, tol, info, wa, lwa )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Returned value of INFO = ', info
      call r8vec_print ( n, x, '  X:' )
      call r8vec_print ( n, fvec, '  F(X):' )

      return
      end
      subroutine f02 ( n, x, fvec, iflag )

c*********************************************************************72
c
cc F02 is a function routine.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n

      double precision fvec(n)
      integer iflag
      double precision x(n)

      fvec(1) = x(1) * x(1) - 10.0D+00 * x(1)
     &  + x(2) * x(2) + 8.0D+00

      fvec(2) = x(1) * x(2) * x(2) + x(1) - 10.0D+00 * x(2)
     &  + 8.0D+00

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests HYBRJ1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 2 )
      integer ldfjac
      parameter ( ldfjac = n )
      integer lwa
      parameter ( lwa = ( n * ( n + 13 ) ) / 2 )

      external f03
      double precision fjac(ldfjac,n)
      double precision fvec(n)
      integer i
      integer iflag
      integer info
      double precision tol
      double precision wa(lwa)
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' )
     &  '  HYBRJ1 solves a nonlinear system of equations.'

      x(1) = 3.0D+00
      x(2) = 0.0D+00
      call r8vec_print ( n, x, '  Initial X:' )
      iflag = 1
      call f02 ( n, x, fvec, iflag )
      call r8vec_print ( n, fvec, '  F(X):' )

      tol = 0.00001D+00

      call hybrj1 ( f03, n, x, fvec, fjac, ldfjac, tol, info, wa, lwa )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Returned value of INFO = ', info
      call r8vec_print ( n, x, '  X:' )
      call r8vec_print ( n, fvec, '  F(X):' )

      return
      end
      subroutine f03 ( n, x, fvec, fjac, ldfjac, iflag )

c*********************************************************************72
c
cc F03 is a function/jacobian routine.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer ldfjac
      integer n

      double precision fjac(ldfjac,n)
      double precision fvec(n)
      integer iflag
      double precision x(n)

      if ( iflag .eq. 1 ) then

        fvec(1) = x(1) * x(1) - 10.0D+00 * x(1) + x(2) * x(2) + 8.0D+00
        fvec(2) = x(1) * x(2) * x(2) + x(1) - 10.0D+00 * x(2) + 8.0D+00

      else if ( iflag .eq. 2 ) then

        fjac(1,1) = 2.0D+00 * x(1) - 10.0D+00
        fjac(1,2) = 2.0D+00 * x(2)
        fjac(2,1) = x(2) * x(2) + 1.0D+00
        fjac(2,2) = 2.0D+00 * x(1) * x(2) - 10.0D+00

      end if

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests LMDER1.
c
c  Discussion:
c
c    LMDER1 solves M nonlinear equations in N unknowns, with M >= N.
c
c    LMDER1 seeks a solution X minimizing the euclidean norm of the residual.
c
c    In this example, the set of equations is actually linear, but
c    normally they are nonlinear.
c
c    In this problem, we have a set of pairs of data points, and we
c    seek a functional relationship between them.  We assume the
c    relationship is of the form
c
c      y = a * x + b
c
c    and we want to know the values of a and b.  Therefore, we would like
c    to find numbers a and b which satisfy a set of equations.
c
c    The data points are (2,2), (4,11), (6,28) and (8,40).
c
c    Therefore, the equations we want to satisfy are:
c
c      a * 2 + b -  2 = 0
c      a * 4 + b - 11 = 0
c      a * 6 + b - 28 = 0
c      a * 8 + b - 40 = 0
c
c    The least squares solution of this system is a=6.55, b=-12.5,
c    In other words, the line y=6.55*x-12.5 is the line which "best"
c    models the data in the least squares sense.
c
c    Problems with more variables, or higher degree polynomials, would
c    be solved similarly.  For example, suppose we have (x,y,z) data,
c    and we wish to find a relationship of the form f(x,y,z).  We assume
c    that x and y occur linearly, and z quadratically.  Then the equation
c    we seek has the form:
c
c      a*x+b*y+c*z + d*z*z + e = 0
c
c    and, supposing that our first two points were (1,2,3), (1,3,8), our set of
c    equations would begin:
c
c      a*1+b*2+c*3 + d*9  + e = 0
c      a*1+b*3+c*8 + d*64 + e = 0
c
c    and so on.
c
c    M is the number of equations, which in this case is the number of
c    (x,y) data values.
c
c    N is the number of variables, which in this case is the number of
c    'free' coefficients in the relationship we are trying to determine.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 4 )
      integer n
      parameter ( n = 2 )
      integer ldfjac
      parameter ( ldfjac = m )
      integer lwa
      parameter ( lwa = 5 * n + m )

      external f04
      double precision fjac(ldfjac,n)
      double precision fvec(m)
      integer i
      integer iflag
      integer info
      integer ipvt(n)
      double precision tol
      double precision wa(lwa)
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' )
     &  '  LMDER1 minimizes M functions in N variables.'

      x(1) = 0.0D+00
      x(2) = 5.0D+00
      call r8vec_print ( n, x, '  Initial X:' )
      iflag = 1
      call f04 ( m, n, x, fvec, fjac, ldfjac, iflag )
      call r8vec_print ( m, fvec, '  F(X):' )

      tol = 0.00001D+00

      call lmder1 ( f04, m, n, x, fvec, fjac, ldfjac, tol, info,
     &  ipvt, wa, lwa )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Returned value of INFO = ', info
      call r8vec_print ( n, x, '  X:' )
      call r8vec_print ( m, fvec, '  F(X):' )

      return
      end
      subroutine f04 ( m, n, x, fvec, fjac, ldfjac, iflag )

c*********************************************************************72
c
cc F04 is a function/jacobian routine.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer ldfjac
      integer m
      integer n

      double precision fjac(ldfjac,n)
      double precision fvec(m)
      integer i
      integer iflag
      double precision x(n)
      double precision xdat(4)
      double precision ydat(4)

      save xdat
      save ydat

      data xdat / 2.0D+00,  4.0D+00,  6.0D+00,  8.0D+00 /
      data ydat / 2.0D+00, 11.0D+00, 28.0D+00, 40.0D+00 /

      if ( iflag .eq. 1 ) then

        do i = 1, m
          fvec(i) = x(1) * xdat(i) + x(2) - ydat(i)
        end do

      else if ( iflag .eq. 2 ) then

        do i = 1, m
          fjac(i,1) = xdat(i)
          fjac(i,2) = 1.0D+00
        end do

      end if

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 tests LMDER1.
c
c  Discussion:
c
c    LMDER1 solves M nonlinear equations in n unknowns, where M is greater
c    than N.  The functional fit is nonlinear this time, of the form
c
c      y=a+b*x**c,
c
c    with x and y data, and a, b and c unknown.
c
c    This problem is set up so that the data is exactly fit by by
c    a=1, b=3, c=2.  Normally, the data would only be approximately
c    fit by the best possible solution.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 10 )
      integer n
      parameter ( n = 3 )
      integer ldfjac
      parameter ( ldfjac = m )
      integer lwa
      parameter ( lwa = 5 * n + m )

      external f05
      double precision fjac(ldfjac,n)
      double precision fvec(m)
      integer i
      integer iflag
      integer info
      integer ipvt(n)
      double precision tol
      double precision wa(lwa)
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' )
     &  '  LMDER1 minimizes M functions in N variables.'

      x(1) = 0.0D+00
      x(2) = 5.0D+00
      x(3) = 1.3D+00
      call r8vec_print ( n, x, '  Initial X:' )
      iflag = 1
      call f05 ( m, n, x, fvec, fjac, ldfjac, iflag )
      call r8vec_print ( m, fvec, '  F(X):' )

      tol = 0.00001D+00

      call lmder1 ( f05, m, n, x, fvec, fjac, ldfjac, tol, info,
     &  ipvt, wa, lwa )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Returned value of INFO = ', info
      call r8vec_print ( n, x, '  X:' )
      call r8vec_print ( m, fvec, '  F(X):' )

      return
      end
      subroutine f05 ( m, n, x, fvec, fjac, ldfjac, iflag )

c*********************************************************************72
c
cc F05 is a function/jacobian routine.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer ldfjac
      integer m
      integer n

      double precision fjac(ldfjac,n)
      double precision fvec(m)
      integer i
      integer iflag
      double precision x(n)
      double precision xdat(10)
      double precision ydat(10)

      save xdat
      save ydat

      data xdat /
     &  1.0D+00, 2.0D+00, 3.0D+00, 4.0D+00, 5.0D+00,
     &  6.0D+00, 7.0D+00, 8.0D+00, 9.0D+00, 10.0D+00 /
      data ydat /
     &  4.0D+00, 13.0D+00, 28.0D+00, 49.0D+00, 76.0D+00,
     &  109.0D+00, 148.0D+00, 193.0D+00, 244.0D+00, 301.0D+00 /

      if ( iflag .eq. 1 ) then

        do i = 1, m
          fvec(i) = x(1) + x(2) * xdat(i)**x(3) - ydat(i)
        end do

      else if ( iflag .eq. 2 ) then

        do i = 1, m
          fjac(i,1) = 1.0D+00
          fjac(i,2) = xdat(i)**x(3)
          fjac(i,3) = x(2) * log ( xdat(i) ) * xdat(i)**x(3)
        end do

      end if

      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 tests LMDIF1.
c
c  Discussion:
c
c    LMDIF1 solves M nonlinear equations in N unknowns, where M is greater
c    than N.  Generally, you cannot get a solution vector x which will satisfy
c    all the equations.  That is, the vector equation f(x)=0 cannot
c    be solved exactly.  Instead, minpack seeks a solution x so that
c    the euclidean norm transpose(f(x))*f(x) is minimized.  The size
c    of the euclidean norm is a measure of how good the solution is.
c
c    In this example, the set of equations is actually linear, but
c    normally they are nonlinear.
c
c    In this problem, we have a set of pairs of data points, and we
c    seek a functional relationship between them.  We assume the
c    relationship is of the form
c
c      y=a*x+b
c
c    and we want to know the values of a and b.  Therefore, we would like
c    to find numbers a and b which satisfy a set of equations.
c
c    The data points are (2,2), (4,11), (6,28) and (8,40).
c
c    Therefore, the equations we want to satisfy are:
c
c      a * 2 + b -  2 = 0
c      a * 4 + b - 11 = 0
c      a * 6 + b - 28 = 0
c      a * 8 + b - 40 = 0
c
c    The least squares solution of this system is a=6.55, b=-12.5,
c    In other words, the line y=6.55*x-12.5 is the line which "best"
c    models the data in the least squares sense.
c
c    Problems with more variables, or higher degree polynomials, would
c    be solved similarly.  For example, suppose we have (x,y,z) data,
c    and we wish to find a relationship of the form f(x,y,z).  We assume
c    that x and y occur linearly, and z quadratically.  Then the equation
c    we seek has the form:
c
c      a*x+b*y+c*z + d*z*z + e = 0
c
c    and, supposing that our first two points were (1,2,3), (1,3,8), our set of
c    equations would begin:
c
c      a*1+b*2+c*3 + d*9  + e = 0
c      a*1+b*3+c*8 + d*64 + e = 0
c
c    and so on.
c
c    M is the number of equations, which in this case is the number of
c    (x,y) data values.
c
c    N is the number of variables, which in this case is the number of
c    'free' coefficients in the relationship we are trying to determine.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 4 )
      integer n
      parameter ( n = 2 )
      integer lwa
      parameter ( lwa = m * n + 5 * n + m )

      external f06
      double precision fvec(m)
      integer i
      integer iflag
      integer info
      integer iwa(n)
      double precision tol
      double precision wa(lwa)
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' )
     &  '  LMDIF1 minimizes M functions in N variables.'

      x(1) = 0.0D+00
      x(2) = 5.0D+00
      call r8vec_print ( n, x, '  Initial X:' )
      iflag = 1
      call f06 ( m, n, x, fvec, iflag )
      call r8vec_print ( m, fvec, '  F(X):' )

      tol = 0.00001D+00

      call lmdif1 ( f06, m, n, x, fvec, tol, info, iwa, wa, lwa )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Returned value of INFO = ', info
      call r8vec_print ( n, x, '  X:' )
      call r8vec_print ( m, fvec, '  F(X):' )

      return
      end
      subroutine f06 ( m, n, x, fvec, iflag )

c*********************************************************************72
c
cc F06 is a function routine.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      integer n

      double precision fvec(m)
      integer i
      integer iflag
      double precision x(n)
      double precision xdat(4)
      double precision ydat(4)

      save xdat
      save ydat

      data xdat / 2.0D+00,  4.0D+00,  6.0D+00,  8.0D+00 /
      data ydat / 2.0D+00, 11.0D+00, 28.0D+00, 40.0D+00 /

      do i = 1, m
        fvec(i) = x(1) * xdat(i) + x(2) - ydat(i)
      end do

      return
      end
      subroutine test07 ( )

c*********************************************************************72
c
cc TEST07 tests LMDIF1.
c
c  Discussion:
c
c    LMDIF1 solves M nonlinear equations in N unknowns, where M is greater
c    than N.  It is similar to test02, except that the functional fit is
c    nonlinear this time, of the form
c
c      y = a + b * x**c,
c
c    with x and y data, and a, b and c unknown.
c
c    This problem is set up so that the data is exactly fit by by
c    a=1, b=3, c=2.  Normally, the data would only be approximately
c    fit by the best possible solution.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 10 )
      integer n
      parameter ( n = 3 )
      integer lwa
      parameter ( lwa = m * n + 5 * n + m )

      external f07
      double precision fvec(m)
      integer i
      integer iflag
      integer info
      integer iwa(n)
      double precision tol
      double precision wa(lwa)
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07'
      write ( *, '(a)' )
     &  '  LMDIF1 minimizes M functions in N variables.'

      x(1) = 0.0D+00
      x(2) = 5.0D+00
      x(3) = 1.3D+00
      call r8vec_print ( n, x, '  Initial X:' )
      iflag = 1
      call f07 ( m, n, x, fvec, iflag )
      call r8vec_print ( m, fvec, '  F(X):' )

      tol = 0.00001D+00

      call lmdif1 ( f07, m, n, x, fvec, tol, info, iwa, wa, lwa )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Returned value of INFO = ', info
      call r8vec_print ( n, x, '  X:' )
      call r8vec_print ( m, fvec, '  F(X):' )

      return
      end
      subroutine f07 ( m, n, x, fvec, iflag )

c*********************************************************************72
c
cc F07 is a function routine.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      integer n

      double precision fvec(m)
      integer i
      integer iflag
      double precision x(n)
      double precision xdat(10)
      double precision ydat(10)

      data xdat /
     &  1.0D+00, 2.0D+00, 3.0D+00, 4.0D+00, 5.0D+00,
     &  6.0D+00, 7.0D+00, 8.0D+00, 9.0D+00, 10.0D+00 /
      data ydat /
     &  4.0D+00, 13.0D+00, 28.0D+00, 49.0D+00, 76.0D+00,
     &  109.0D+00, 148.0D+00, 193.0D+00, 244.0D+00, 301.0D+00 /

      do i = 1, m
        fvec(i) = x(1) + x(2) * xdat(i)**x(3) - ydat(i)
      end do

      return
      end
      subroutine test08 ( )

c*********************************************************************72
c
cc TEST08 tests LMSTR1.
c
c  Discussion:
c
c    LMSTR1 solves M nonlinear equations in N unknowns, where M is greater
c    than N.  Generally, you cannot get a solution vector x which will satisfy
c    all the equations.  That is, the vector equation f(x)=0 cannot
c    be solved exactly.  Instead, minpack seeks a solution x so that
c    the euclidean norm transpose(f(x))*f(x) is minimized.  The size
c    of the euclidean norm is a measure of how good the solution is.
c
c    In this example, the set of equations is actually linear, but
c    normally they are nonlinear.
c
c    In this problem, we have a set of pairs of data points, and we
c    seek a functional relationship between them.  We assume the
c    relationship is of the form
c
c      y=a*x+b
c
c    and we want to know the values of a and b.  Therefore, we would like
c    to find numbers a and b which satisfy a set of equations.
c
c    The data points are (2,2), (4,11), (6,28) and (8,40).
c
c    Therefore, the equations we want to satisfy are:
c
c      a * 2 + b -  2 = 0
c      a * 4 + b - 11 = 0
c      a * 6 + b - 28 = 0
c      a * 8 + b - 40 = 0
c
c    The least squares solution of this system is a=6.55, b=-12.5,
c    In other words, the line y=6.55*x-12.5 is the line which "best"
c    models the data in the least squares sense.
c
c    Problems with more variables, or higher degree polynomials, would
c    be solved similarly.  For example, suppose we have (x,y,z) data,
c    and we wish to find a relationship of the form f(x,y,z).  We assume
c    that x and y occur linearly, and z quadratically.  Then the equation
c    we seek has the form:
c
c      a*x + b*y + c*z + d*z*z + e = 0
c
c    and, supposing that our first two points were (1,2,3), (1,3,8), our set of
c    equations would begin:
c
c      a*1 + b*2 + c*3 + d*9  + e = 0
c      a*1 + b*3 + c*8 + d*64 + e = 0
c
c    and so on.
c
c    M is the number of equations, which in this case is the number of
c    (x,y) data values.
c
c    N is the number of variables, which in this case is the number of
c    'free' coefficients in the relationship we are trying to determine.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 4 )
      integer n
      parameter ( n = 2 )
      integer ldfjac
      parameter ( ldfjac = m )
      integer lwa
      parameter ( lwa = 5 * n + m )

      external f08
      double precision fjac(ldfjac,n)
      integer fjrow
      double precision fvec(m)
      integer i
      integer iflag
      integer info
      integer ipvt(n)
      double precision tol
      double precision wa(lwa)
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08'
      write ( *, '(a)' )
     &  '  LMSTR1 minimizes M functions in N variables.'

      x(1) = 0.0D+00
      x(2) = 5.0D+00
      call r8vec_print ( n, x, '  Initial X:' )
      iflag = 1
      call f08 ( m, n, x, fvec, fjrow, iflag )
      call r8vec_print ( m, fvec, '  F(X):' )

      tol = 0.00001D+00

      call lmstr1 ( f08, m, n, x, fvec, fjac, ldfjac, tol, info,
     &  ipvt, wa, lwa )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Returned value of INFO = ', info
      call r8vec_print ( n, x, '  X:' )
      call r8vec_print ( m, fvec, '  F(X):' )

      return
      end
      subroutine f08 ( m, n, x, fvec, fjrow, iflag )

c*********************************************************************72
c
cc F08 is a function/jacobian routine.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      integer n

      double precision fjrow(n)
      double precision fvec(m)
      integer i
      integer iflag
      double precision x(n)
      double precision xdat(4)
      double precision ydat(4)

      save xdat
      save ydat

      data xdat / 2.0D+00,  4.0D+00,  6.0D+00,  8.0D+00 /
      data ydat / 2.0D+00, 11.0D+00, 28.0D+00, 40.0D+00 /

      if ( iflag .eq. 1 ) then

        do i = 1, m
          fvec(i) = x(1) * xdat(i) + x(2) - ydat(i)
        end do

      else

        fjrow(1) = xdat(iflag-1)
        fjrow(2) = 1.0D+00

      end if

      return
      end
      subroutine test09 ( )

c*********************************************************************72
c
cc TEST09 tests LMSTR1.
c
c  Discussion:
c
c    LMSTR1 solves M nonlinear equations in N unknowns, where M is greater
c    than N.  This test is similar to test02, except that the functional fit
c    is nonlinear this time, of the form
c
c      y = a + b * x**c,
c
c    with x and y data, and a, b and c unknown.
c
c    This problem is set up so that the data is exactly fit by by
c    a=1, b=3, c=2.  Normally, the data would only be approximately
c    fit by the best possible solution.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 10 )
      integer n
      parameter ( n = 3 )
      integer ldfjac
      parameter ( ldfjac = m )
      integer lwa
      parameter ( lwa = 5 * n + m )

      external f09
      double precision fjac(ldfjac,n)
      double precision fjrow(n)
      double precision fvec(m)
      integer i
      integer iflag
      integer info
      integer ipvt(n)
      double precision tol
      double precision wa(lwa)
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST09'
      write ( *, '(a)' )
     &  '  LMSTR1 minimizes M functions in N variables.'

      x(1) = 0.0D+00
      x(2) = 5.0D+00
      x(3) = 1.3D+00
      call r8vec_print ( n, x, '  Initial X:' )
      iflag = 1
      call f09 ( m, n, x, fvec, fjrow, iflag )
      call r8vec_print ( m, fvec, '  F(X):' )

      tol = 0.00001D+00

      call lmstr1 ( f09, m, n, x, fvec, fjac, ldfjac, tol, info,
     &  ipvt, wa, lwa )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Returned value of INFO = ', info
      call r8vec_print ( n, x, '  X:' )
      call r8vec_print ( m, fvec, '  F(X):' )

      return
      end
      subroutine f09 ( m, n, x, fvec, fjrow, iflag )

c*********************************************************************72
c
cc F09 is a function/jacobian routine.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      integer n

      double precision fjrow(n)
      double precision fvec(m)
      integer i
      integer iflag
      double precision x(n)
      double precision xdat(10)
      double precision ydat(10)

      save xdat
      save ydat

      data xdat /
     &  1.0D+00, 2.0D+00, 3.0D+00, 4.0D+00, 5.0D+00,
     &  6.0D+00, 7.0D+00, 8.0D+00, 9.0D+00, 10.0D+00 /
      data ydat /
     &  4.0D+00, 13.0D+00, 28.0D+00, 49.0D+00, 76.0D+00,
     &  109.0D+00, 148.0D+00, 193.0D+00, 244.0D+00, 301.0D+00 /

      if ( iflag .eq. 1 ) then

        do i = 1, m
          fvec(i) = x(1) + x(2) * xdat(i)**x(3) - ydat(i)
        end do

      else

        fjrow(1) = 1.0D+00
        fjrow(2) = xdat(iflag-1)**x(3)
        fjrow(3) = x(2) * log ( xdat(iflag-1) ) * xdat(iflag-1)**x(3)

      end if

      return
      end
