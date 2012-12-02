      program main

c*********************************************************************72
c
cc MAIN tests SORT.
c
c  Modified:
c
c    04 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer ii
      integer jj

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS347_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS347 library.'

      ii = 1
      jj = 20
      call test01 ( ii, jj )

      ii = 5
      jj = 18
      call test01 ( ii, jj )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS347_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( ii, jj )

c*********************************************************************72
c
cc TEST01 tests SORT on a particular range of indices.
c
c  Modified:
c
c    04 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 20 )

      integer a(n)
      integer i
      integer ii
      integer jj
      real rn
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  SORT ascending sorts an integer vector.'
      write ( *, '(a,i6)' ) '  Here we sort entries II = ', ii
      write ( *, '(a,i6)' ) '  through JJ = ', jj

      seed = 123456789
      do i = 1, n
        a(i) = int ( n * rn ( seed ) )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Unsorted array:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2x,i6)' ) i, a(i)
      end do

      call sort ( a, ii, jj )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Sorted array:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2x,i6)' ) i, a(i)
      end do

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
c    16 September 2005
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

      character ( len = 8 ) date
      character ( len = 10 ) time

      call date_and_time ( date, time )

      write ( *, '(a8,2x,a10)' ) date, time

      return
      end
