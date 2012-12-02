      program main

c*********************************************************************72
c
cc TOMS454_PRB2 tests JCONSX, ACM TOMS algorithm 454.
c
c  Modified:
c
c    23 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer k
      integer l
      integer m
      integer n

      parameter ( k = 6 )

      parameter ( l = 4 )
      parameter ( m = 4 )
      parameter ( n = 3 )

      real alpha
      real beta
      real delta
      real f(k)
      real g(m)
      integer gamma
      real h(m)
      integer i
      integer iev2
      integer it
      integer itmax
      integer j
      integer k0
      real r(k,n)
      real rn
      integer seed
      real x(k,l)
      real xc(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS454_PRB2'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test TOMS algorithm 454, which'
      write ( *, '(a)' ) '  implements the complex method of'
      write ( *, '(a)' ) '  constrained optimization.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Here we consider the Post Office problem.'

      itmax = 100
      alpha = 1.3E+00
      beta = 0.01E+00
      gamma = 5
      delta = 0.0001E+00
      k0 = 6

      x(1,1) = 1.0E+00
      x(1,2) = 1.0E+00
      x(1,3) = 1.0E+00
      x(1,4) = x(1,1) + 2.0E+00 * x(1,2) + 2.0E+00 * x(1,3)

      do i = 2, k
        do j = 1, l
          x(i,j) = 0.0E+00
        end do
      end do

      seed = 123456789

      do i = 1, k
        do j = 1, n
          r(i,j) = rn ( seed )
        end do
      end do

      call jconsx ( n, m, k, itmax, alpha, beta, gamma, 
     &  delta, x, r, f, it, iev2, k0, g, h, xc, l )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The exact solution is:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    F(24,12,12) = 3456'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS454_PRB2'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
      subroutine jfunc ( n, m, k, x, f, i, l )

c*********************************************************************72
c
cc JFUNC evaluates the function F(I) at X(I,1:L).
c
c  Discussion:
c
c    This is the objective function for the Post Office problem.
c
      implicit none

      integer k
      integer l

      real f(k)
      integer i
      integer m
      integer n
      real x(k,l)

      f(i) = x(i,1) * x(i,2) * x(i,3)

      return
      end
      subroutine jcnst1 ( n, m, k, x, g, h, i, l )

c*********************************************************************72
c
cc JCNST1 applies the constraints.
c
      implicit none

      integer k
      integer l
      integer m
      integer n

      real g(m)
      real h(m)
      integer i
      real x(k,l)
c
c  Apparently, if you have "implicit variables", you need to set them
c  here, every time this routine is called.
c
      x(i,4) = x(i,1) + 2.0E+00 * x(i,2) + 2.0E+00 * x(i,3)

      g(1) = 0.0E+00
      g(2) = 0.0E+00
      g(3) = 0.0E+00
      g(4) = 0.0E+00

      h(1) = 42.0E+00
      h(2) = 42.0E+00
      h(3) = 42.0E+00
      h(4) = 72.0E+00

      return
      end
      function rn ( seed )

c*********************************************************************72
c
cc RN returns a unit single precision pseudorandom number.
c
c  Discussion:
c
c    This routine implements the recursion
c
c      seed = 16807 * seed mod ( 2**31 - 1 )
c      rn = seed / ( 2**31 - 1 )
c
c    The integer arithmetic never requires more than 32 bits,
c    including a sign bit.
c
c    If the initial seed is 12345, then the first three computations are
c
c      Input     Output      RN
c      SEED      SEED
c
c         12345   207482415  0.096616
c     207482415  1790989824  0.833995
c    1790989824  2035175616  0.947702
c
c  Modified:
c
c    11 August 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, L E Schrage,
c    A Guide to Simulation,
c    Springer Verlag, pages 201-202, 1983.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley Interscience, page 95, 1998.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, pages 362-376, 1986.
c
c    P A Lewis, A S Goodman, J M Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, pages 136-143, 1969.
c
c  Parameters:
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, real RN, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      integer k
      integer seed
      real rn

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + 2147483647
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      rn = real ( dble ( seed ) * 4.656612875D-10 )

      return
      end
