      program main

c*********************************************************************72
c
cc MAIN is the main program for COST_TEST.
c
c  Discussion:
c
c    COST_TEST tests COST1B, COST1F, COST1I.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 March 2005
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4096 )

      integer lenwrk
      parameter ( lenwrk = n - 1 )

      integer lensav
      parameter ( lensav = 8204 )

      integer ier
      integer inc
      integer lenr
      real r(n)
      integer seed
      real work(lenwrk)
      real wsave(lensav)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'COST_TEST'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  For real fast cosine transforms, 1D,'
      write ( *, '(a)' ) '  COST1I initializes the transforms,'
      write ( *, '(a)' ) '  COST1F does a forward transforms;'
      write ( *, '(a)' ) '  COST1B does a backward transforms.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  The number of data items is N = ', n
c
c  Set the data values.
c
      seed = 1973

      call r4vec_uniform_01 ( n, seed, r )

      call r4vec_print_some ( n, r, 10, '  The original data:' )
c
c  Allocate and initialize the WSAVE array.
c
      call cost1i ( n, wsave, lensav, ier )
c
c  Compute the FFT coefficients.
c
      inc = 1
      lenr = n

      call cost1f ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

      call r4vec_print_some ( n, r, 10, '  The FFT coefficients:' )
c
c  Compute inverse FFT of coefficients.  Should get back the
c  original data.
c
      call cost1b ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

      call r4vec_print_some ( n, r, 10, '  The retrieved data:' )

      stop
      end
      subroutine r4vec_uniform_01 ( n, seed, r )

c*********************************************************************72
c
cc R4VEC_UNIFORM_01 returns a unit pseudorandom R4VEC.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2005
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
c    Input, integer N, the number of entries in the vector.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, real R(N), the vector of pseudorandom values.
c
      implicit none

      integer n

      integer i
      integer k
      integer seed
      real r(n)

      do i = 1, n

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed < 0 ) then
          seed = seed + 2147483647
        end if

        r(i) = real ( seed ) * 4.656612875E-10

      end do

      return
      end
      subroutine r4vec_print_some ( n, a, max_print, title )

c*********************************************************************72
c
cc R4VEC_PRINT_SOME prints "some" of an R4VEC.
c
c  Discussion:
c
c    The user specifies MAX_PRINT, the maximum number of lines to print.
c
c    If N, the size of the vector, is no more than MAX_PRINT, then
c    the entire vector is printed, one entry per line.
c
c    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
c    followed by a line of periods suggesting an omission,
c    and the last entry.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 December 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries of the vector.
c
c    Input, real A(N), the vector to be printed.
c
c    Input, integer MAX_PRINT, the maximum number of lines to print.
c
c    Input, character ( len = * ) TITLE, an optional title.
c
      implicit none

      integer n

      real a(n)
      integer i
      integer max_print
      character ( len = * ) title

      if ( max_print <= 0 ) then
        return
      end if

      if ( n <= 0 ) then
        return
      end if

      if ( 0 < len_trim ( title ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) title
        write ( *, '(a)' ) ' '
      end if

      if ( n <= max_print ) then

        do i = 1, n
          write ( *, '(2x,i6,2x,g14.6)' ) i, a(i)
        end do

      else if ( 3 <= max_print ) then

        do i = 1, max_print-2
          write ( *, '(2x,i6,2x,g14.6)' ) i, a(i)
        end do
        write ( *, '(a)' ) '  ......  ..............'
        i = n
        write ( *, '(2x,i6,2x,g14.6)' ) i, a(i)

      else

        do i = 1, max_print - 1
          write ( *, '(2x,i6,2x,g14.6)' ) i, a(i)
        end do
        i = max_print
        write ( *, '(2x,i6,2x,g14.6,2x,a)' ) i, a(i), '...more...'

      end if

      return
      end
