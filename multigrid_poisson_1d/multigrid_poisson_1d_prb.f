      program main

c*********************************************************************72
c
cc MULTIGRID_POISSON_1D_PRB tests MULTIGRID_POISSON_1D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 December 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MULTIGRID_POISSON_1D_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) 
     &  '  Test the MULTIGRID_POISSON_1D multigrid library.'

      call test01_mono ( ) 
      call test01_multi ( ) 
      call test02_mono ( ) 
      call test02_multi ( ) 
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MULTIGRID_POISSON_1D_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01_mono ( ) 

c*********************************************************************72
c
cc TEST01_MONO tests MONOGRID_POISSON_1D on test case 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision difmax
      external exact1
      double precision exact1
      external force1
      double precision force1
      integer i
      integer it_num
      integer k
      integer n
      double precision u(33)
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01_MONO'
      write ( *, '(a)' ) 
     &  '  MONOGRID_POISSON_1D solves a 1D Poisson BVP'
      write ( *, '(a)' ) '  using the Gauss-Seidel method.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  -u"(x) = 1, for 0 < x < 1'
      write ( *, '(a)' ) '  u(0) = u(1) = 0.'
      write ( *, '(a)' ) '  Solution is u(x) = ( -x^2 + x ) / 2'

      do k = 5, 5

        n = 2**k

        write ( *, '(a)' ) ' '
        write ( *, '(a,i4)' ) '  Mesh index K = ', k
        write ( *, '(a,i6)' ) '  Number of intervals N=2^K = ', n
        write ( *, '(a,i6)' ) '  Number of nodes = 2^K+1 =   ', 2**k + 1

        call monogrid_poisson_1d ( n, force1, exact1, it_num, u )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &    '     I        X(I)      U(I)         U Exact(X(I))'
        write ( *, '(a)' ) ' '
        do i = 1, n + 1
          x = dble ( i - 1 ) / dble ( n )
          write(6,'(2x,i4,2x,f10.4,2x,g14.6,2x,g14.6)') 
     &      i, x, u(i), exact1 ( x )
        end do

        write ( *, '(a)' ) ' '

        difmax = 0.0D+00
        do i = 1, n + 1
          x = dble ( i - 1 ) / dble ( n )
          difmax = max ( difmax, abs ( u(i) - exact1 ( x ) ) )
        end do 
        write ( *, '(a,g14.6)' ) '  Maximum error = ', difmax
        write ( *, '(a,i6)' ) '  Number of iterations = ', it_num

      end do

      return
      end
      subroutine test01_multi ( ) 

c*********************************************************************72
c
cc TEST01_MULTI tests MULTIGRID_POISSON_1D on test case 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 128 )

      double precision difmax
      external exact1
      double precision exact1
      external force1
      double precision force1
      integer i
      integer it_num
      integer k
      integer n
      double precision u(n_max+1)
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01_MULTI'
      write ( *, '(a)' ) 
     &  '  MULTIGRID_POISSON_1D solves a 1D Poisson BVP'
      write ( *, '(a)' ) '  using the multigrid method.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  -u"(x) = 1, for 0 < x < 1'
      write ( *, '(a)' ) '  u(0) = u(1) = 0.'
      write ( *, '(a)' ) '  Solution is u(x) = ( -x^2 + x ) / 2'

      do k = 5, 5

        n = 2**k

        write ( *, '(a)' ) ' '
        write ( *, '(a,i4)' ) '  Mesh index K = ', k
        write ( *, '(a,i6)' ) '  Number of intervals N=2^K = ', n
        write ( *, '(a,i6)' ) '  Number of nodes = 2^K+1 =   ', n + 1

        call multigrid_poisson_1d ( n, force1, exact1, it_num, u )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &    '     I        X(I)      U(I)         U Exact(X(I))'
        write ( *, '(a)' ) ' '
        do i = 1, n + 1
          x = dble ( i - 1 ) / dble ( n )
          write(6,'(2x,i4,2x,f10.4,2x,g14.6,2x,g14.6)') 
     &      i, x, u(i), exact1 ( x )
        end do

        write ( *, '(a)' ) ' '

        difmax = 0.0D+00
        do i = 1, n + 1
          x = dble ( i - 1 ) / dble ( n )
          difmax = max ( difmax, abs ( u(i) - exact1 ( x ) ) )
        end do 
        write ( *, '(a,g14.6)' ) '  Maximum error = ', difmax
        write ( *, '(a,i6)' ) '  Number of iterations = ', it_num

      end do

      return
      end
      function exact1 ( x )

c*********************************************************************72
c                                                    
cc EXACT1 evaluates the exact solution.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 November 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input, double precision X, the evaluation point.
c
c    Output, double precision EXACT1, the value of the exact solution at X.
c
      double precision exact1
      double precision x

      exact1 = 0.5D+00 * ( - x * x + x )

      return
      end
      function force1 ( x )

c*********************************************************************72
c                                                    
cc FORCE1 evaluates the forcing function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 November 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input, double precision X, the evaluation point.
c
c    Output, double precision FORCE1, the value of the forcing function at X.
c
      double precision force1
      double precision x

      force1 = 1.0D+00

      return
      end
      subroutine test02_mono ( ) 

c*********************************************************************72
c
cc TEST02_MONO tests MONOGRID_POISSON_1D on test case 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision difmax
      external exact2
      double precision exact2
      external force2
      double precision force2
      integer i
      integer it_num
      integer k
      integer n
      double precision u(33)
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02_MONO'
      write ( *, '(a)' ) 
     &  '  MONOGRID_POISSON_1D solves a 1D Poisson BVP'
      write ( *, '(a)' ) '  using the Gauss-Seidel method.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  -u"(x) = - x * (x+3) * exp(x), for 0 < x < 1'
      write ( *, '(a)' ) '  u(0) = u(1) = 0.'
      write ( *, '(a)' ) '  Solution is u(x) = x * (x-1) * exp(x)'

      do k = 5, 5

        n = 2**k

        write ( *, '(a)' ) ' '
        write ( *, '(a,i4)' ) '  Mesh index K = ', k
        write ( *, '(a,i6)' ) '  Number of intervals N=2^K = ', n
        write ( *, '(a,i6)' ) '  Number of nodes = 2^K+1 =   ', n + 1

        call monogrid_poisson_1d ( n, force2, exact2, it_num, u )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &    '     I        X(I)      U(I)         U Exact(X(I))'
        write ( *, '(a)' ) ' '
        do i = 1, n + 1
          x = dble ( i - 1 ) / dble ( n )
          write(6,'(2x,i4,2x,f10.4,2x,g14.6,2x,g14.6)') 
     &      i, x, u(i), exact2 ( x )
        end do

        write ( *, '(a)' ) ' '

        difmax = 0.0D+00
        do i = 1, n + 1
          x = dble ( i - 1 ) / dble ( n )
          difmax = max ( difmax, abs ( u(i) - exact2 ( x ) ) )
        end do 
        write ( *, '(a,g14.6)' ) '  Maximum error = ', difmax
        write ( *, '(a,i6)' ) '  Number of iterations = ', it_num

      end do

      return
      end
      subroutine test02_multi ( ) 

c*********************************************************************72
c
cc TEST02_MULTI tests MULTIGRID_POISSON_1D on test case 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 128 )

      double precision difmax
      external exact2
      double precision exact2
      external force2
      double precision force2
      integer i
      integer it_num
      integer k
      integer n
      double precision u(n_max+1)
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02_MULTI'
      write ( *, '(a)' ) 
     &  '  MULTIGRID_POISSON_1D solves a 1D Poisson BVP'
      write ( *, '(a)' ) '  using the multigrid method.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  -u"(x) = - x * (x+3) * exp(x), for 0 < x < 1'
      write ( *, '(a)' ) '  u(0) = u(1) = 0.'
      write ( *, '(a)' ) '  Solution is u(x) = x * (x-1) * exp(x)'

      do k = 5, 5

        n = 2**k

        write ( *, '(a)' ) ' '
        write ( *, '(a,i4)' ) '  Mesh index K = ', k
        write ( *, '(a,i6)' ) '  Number of intervals N=2^K = ', n
        write ( *, '(a,i6)' ) '  Number of nodes = 2^K+1 =   ', n + 1

        call multigrid_poisson_1d ( n, force2, exact2, it_num, u )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &    '     I        X(I)      U(I)         U Exact(X(I))'
        write ( *, '(a)' ) ' '
        do i = 1, n + 1
          x = dble ( i - 1 ) / dble ( n )
          write(6,'(2x,i4,2x,f10.4,2x,g14.6,2x,g14.6)') 
     &      i, x, u(i), exact2 ( x )
        end do

        write ( *, '(a)' ) ' '

        difmax = 0.0D+00
        do i = 1, n + 1
          x = dble ( i - 1 ) / dble ( n )
          difmax = max ( difmax, abs ( u(i) - exact2 ( x ) ) )
        end do 
        write ( *, '(a,g14.6)' ) '  Maximum error = ', difmax
        write ( *, '(a,i6)' ) '  Number of iterations = ', it_num

      end do

      return
      end
      function exact2 ( x )

c*********************************************************************72
c                                                    
cc EXACT2 evaluates the exact solution.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 November 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input, double precision X, the evaluation point.
c
c    Output, double precision EXACT2, the value of the exact solution at X.
c
      double precision exact2
      double precision x

      exact2 = x * ( x - 1.0D+00 ) * exp ( x )

      return
      end
      function force2 ( x )

c*********************************************************************72
c                                                    
cc FORCE2 evaluates the forcing function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 November 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input, double precision X, the evaluation point.
c
c    Output, double precision FORCE2, the value of the forcing function at X.
c
      double precision force2
      double precision x

      force2 = - x * ( x + 3.0D+00 ) * exp ( x )

      return
      end
