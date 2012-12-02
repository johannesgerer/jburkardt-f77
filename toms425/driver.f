      program main

c***********************************************************************
c
cc TOMS425_PRB tests RNVR.
c
      implicit none

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS425_PRB'
      write ( *, '(a)' ) '  Test TOMS algorithm 425, multivariate'
      write ( *, '(a)' ) '  random normal variable generation.'

      call test01
      call test02
      call test03

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS425_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
      subroutine test01

c***********************************************************************
c
cc TEST01 tests RNVR.
c
      integer ni
      integer nv
      parameter ( ni = 2 )
      parameter ( nv = 2 )

      real a(nv,nv)
      integer i, j
      integer iarg
      integer ient
      real x(ni)

      a(1,1) = 1.0
      a(1,2) = 0.0
      a(2,1) = 0.0
      a(2,2) = 1.0

      ient = -1
      iarg = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Test RNVR'
      write ( *, '(a)' ) '  Dimension is 2.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Covariance matrix:'
      write ( *, '(a)' ) '    1.0  0.0'
      write ( *, '(a)' ) '    0.0  1.0'
      write ( *, '(a)' ) ' '

      do i = 1, 10

        call rnvr ( x, a, nv, ni, ient, iarg )
        write ( *, '(2x,g14.6,2x,g14.6)' ) ( x(j), j = 1, ni )

      end do

      return
      end
      subroutine test02

c***********************************************************************
c
cc TEST02 tests RNVR.
c
      integer ni
      integer nv
      parameter ( ni = 2 )
      parameter ( nv = 2 )

      real a(nv,nv)
      integer i, j
      integer iarg
      integer ient
      real x(ni)

      a(1,1) = 1.0
      a(1,2) = 0.3
      a(2,1) = 0.3
      a(2,2) = 1.0

      ient = -1
      iarg = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  Test RNVR'
      write ( *, '(a)' ) '  Dimension is 2.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Covariance matrix:'
      write ( *, '(a)' ) '    1.0  0.3'
      write ( *, '(a)' ) '    0.3  1.0'
      write ( *, '(a)' ) ' '

      do i = 1, 10

        call rnvr ( x, a, nv, ni, ient, iarg )
        write ( *, '(2x,g14.6,2x,g14.6)' ) ( x(j), j = 1, ni )

      end do

      return
      end
      subroutine test03

c***********************************************************************
c
cc TEST03 tests RNVR.
c
      integer ni
      integer nv
      parameter ( ni = 2 )
      parameter ( nv = 2 )

      real a(nv,nv)
      integer i, j
      integer iarg
      integer ient
      real x(ni)

      a(1,1) = 1.00
      a(1,2) = 0.99
      a(2,1) = 0.99
      a(2,2) = 1.00

      ient = -1
      iarg = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  Test RNVR'
      write ( *, '(a)' ) '  Dimension is 2.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Covariance matrix:'
      write ( *, '(a)' ) '    1.00  0.99'
      write ( *, '(a)' ) '    0.99  1.00'
      write ( *, '(a)' ) ' '

      do i = 1, 10

        call rnvr ( x, a, nv, ni, ient, iarg )
        write ( *, '(2x,g14.6,2x,g14.6)' ) ( x(j), j = 1, ni )

      end do

      return
      end
      function rn ( seed )

c*******************************************************************************
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
