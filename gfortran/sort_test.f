      program main

c*********************************************************************72
c
cc MAIN is the main program for SORT_TEST.
c
c  Discussion:
c
c    SORT_TEST tests a sorting routine.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      real a(n)
      real ahi
      real alo
      integer i
      integer seed

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SORT_TEST'
      write ( *, '(a)' ) '  FORTRAN77 version.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Demonstrate the use of R4VEC_SORT_BUBBLE'
      write ( *, '(a)' ) '  to sort a real array using bubble sort.'

      seed = 123456789

      do i = 1, 2

        alo = 10.0E+00
        ahi = 25.0E+00

        call r4vec_uniform ( n, alo, ahi, seed, a )

        call r4vec_print ( n, a, '  Unsorted array:' )

        call r4vec_sort_bubble_a ( n, a )

        call r4vec_print ( n, a, '  Sorted array:' )

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SORT_TEST:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine r4vec_print ( n, a, title )

c*********************************************************************72
c
cc R4VEC_PRINT prints an R4VEC.
c
c  Discussion:
c
c    An R4VEC is an array of real values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
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
c    Input, real A(N), the vector to be printed.
c
c    Input, character * ( * ) TITLE, an optional title.
c
      implicit none

      integer n

      real a(n)
      integer i
      integer s_len_trim
      character ( len = * ) title
      integer title_length

      title_length = len_trim ( title )
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
      subroutine r4vec_sort_bubble_a ( n, a )

c*********************************************************************72
c
cc R4VEC_SORT_BUBBLE_A ascending sorts an R4VEC using bubble sort.
c
c  Discussion:
c
c    An R4VEC is a vector of R4's.
c
c    Bubble sort is simple to program, but inefficient.  It should not
c    be used for large arrays.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 December 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input/output, real A(N).
c    On input, an unsorted array.
c    On output, the array has been sorted.
c
      implicit none

      integer n

      real a(n)
      integer i
      integer j
      real t

      do i = 1, n - 1 
        do j = i + 1, n
          if ( a(j) .lt. a(i) ) then
            t    = a(i)
            a(i) = a(j)
            a(j) = t
          end if
        end do
      end do

      return
      end
      subroutine r4vec_uniform ( n, a, b, seed, r )

c*********************************************************************72
c
cc R4VEC_UNIFORM returns a scaled pseudorandom R4VEC.
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
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Second Edition,
c    Springer, 1987,
c    ISBN: 0387964673,
c    LC: QA76.9.C65.B73.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, December 1986, pages 362-376.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley, 1998,
c    ISBN: 0471134031,
c    LC: T57.62.H37.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, Number 2, 1969, pages 136-143.
c
c  Parameters:
c
c    Input, integer M, the number of entries in the vector.
c
c    Input, real A, B, the lower and upper limits.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, real R(N), the vector of pseudorandom values.
c
      implicit none

      integer n

      real a
      real b
      integer i
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer k
      integer seed
      real r(n)

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R4VEC_UNIFORM - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      do i = 1, n

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed .lt. 0 ) then
          seed = seed + i4_huge
        end if
    
        r(i) = a + ( b - a ) * real ( seed ) * 4.656612875E-10

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
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
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
