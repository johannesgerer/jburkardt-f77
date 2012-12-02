      program main

c*********************************************************************72
c
cc MAIN is the main program for MULTITASK_OPENMP.
c
c  Discussion:
c
c    This program demonstrates how OpenMP can be used for multitasking, that 
c    is, a simple kind of parallel processing in which a certain number of 
c    perhaps quite unrelated tasks must be done.
c
c    The OpenMP SECTIONS directive identifies the portion of the program where
c    the code for these tasks is given.
c
c    The OpenMP SECTION directive is used repeatedly to divide this area of
c    the program into independent tasks.
c
c    The code will get the benefit of parallel processing up to the point where
c    there are as many threads as there are tasks.
c
c    The code will get a substantial speedup if the tasks take roughly the
c    same amount of time.  However, if one task takes substantially more time
c    than the others, this results in a limit to the parallel speedup that is
c    possible.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 October 2011
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

      include 'omp_lib.h'

      integer prime_num
      parameter ( prime_num = 20000 )
      integer sine_num
      parameter ( sine_num = 20000 )
 
      integer primes(prime_num)
      double precision sines(sine_num)
      double precision wtime
      double precision wtime1
      double precision wtime2

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MULTITASK_OPENMP:'
      write ( *, '(a)' ) '  FORTRAN77/OpenMP version'
      write ( *, '(a)' ) 
     &  '  Demonstrate how OpenMP can "multitask" by using the'
      write ( *, '(a)' ) 
     &  '  SECTIONS directive to carry out several tasks in parallel.'

      wtime = omp_get_wtime ( )

c$omp parallel shared ( primes, sines )

c$omp sections

c$omp section
        wtime1 = omp_get_wtime ( )
        call prime_table ( prime_num, primes )
        wtime1 = omp_get_wtime ( ) - wtime1

c$omp section
        wtime2 = omp_get_wtime ( )
        call sine_table ( sine_num, sines )
        wtime2 = omp_get_wtime ( ) - wtime2

c$omp end sections

c$omp end parallel

      wtime = omp_get_wtime ( ) - wtime

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) 
     &  '  Number of primes computed was ', prime_num
      write ( *, '(a,i12)' ) '  Last prime was ', primes(prime_num)
      write ( *, '(a,i6)' ) '  Number of sines computed was ', sine_num
      write ( *, '(a,g14.6)' ) 
     &  '  Last sine computed was ', sines(sine_num)
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Elapsed time = ', wtime
      write ( *, '(a,g14.6)' ) '  Task 1 time = ', wtime1
      write ( *, '(a,g14.6)' ) '  Task 2 time = ', wtime2
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MULTITASK_OPENMP:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine prime_table ( prime_num, primes )

c*********************************************************************72
c
cc PRIME_TABLE computes a table of the first PRIME_NUM prime numbers.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PRIME_NUM, the number of primes to compute.
c
c    Output, integer PRIMES(PRIME_NUM), the computed primes.
c
      implicit none

      integer prime_num

      integer i
      integer j
      integer p
      logical prime
      integer primes(prime_num)

      i = 2
      p = 0

10    continue

      if ( p .lt. prime_num ) then

        prime = .true.

        do j = 2, i - 1
          if ( mod ( i, j ) .eq. 0 ) then
            prime = .false.
            go to 20
          end if
        end do

20      continue

        if ( prime ) then
          p = p + 1
          primes(p) = i
        end if

        i = i + 1

        go to 10

      end if

      return
      end
      subroutine sine_table ( sine_num, sines )

c*********************************************************************72
c
cc SINE_TABLE computes a table of sines.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer SINE_NUM, the number of sines to compute.
c
c    Output, double precision SINES(SINE_NUM), the sines.
c
      implicit none

      integer sine_num

      double precision a
      integer i
      integer j
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision sines(sine_num)

      do i = 1, sine_num
        sines(i) =  0.0D+00
        do j = 1, i
          a = dble ( j - 1 ) * pi / dble ( sine_num - 1 )
          sines(i) = sines(i) + sin ( a )
        end do
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
