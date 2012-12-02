      program main

c*********************************************************************72
c
cc MAIN is the main program for RANDOM_OPENMP.
c
c  Discussion:
c
c    This program simply explores one issue in the generation of random
c    numbers in a parallel program.  If the random number generator uses
c    an integer seed to determine the next entry, then it is not easy for
c    a parallel program to reproduce the same exact sequence.
c
c    But what is worse is that it might not be clear how the separate
c    OpenMP threads should handle the SEED value - as a shared or private
c    variable?  It seems clear that each thread should have a private
c    seed that is initialized to a distinct value at the beginning of
c    the computation.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 September 2012
c
c  Author:
c
c    John Burkardt
c
      include 'omp_lib.h'

      integer n
      integer seed

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RANDOM_OPENMP'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  An OpenMP program using random numbers.'
      write ( *, '(a)' ) '  The random numbers depend on a seed.'
      write ( *, '(a)' ) '  We need to insure that each OpenMP thread'
      write ( *, '(a)' ) '  starts with a different seed.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 
     &  '  Number of processors available = ', omp_get_num_procs ( )
      write ( * ,'(a,i8)' ) 
     &  '  Number of threads =              ', omp_get_max_threads ( )

      n = 100
      seed = 123456789
      call monte_carlo ( n, seed )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RANDOM_OPENMP'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine monte_carlo ( n, seed )

c*********************************************************************72
c
cc MONTE_CARLO carries out a Monte Carlo calculation with random values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of values to generate.
c
c    Input/output, integer SEED, a seed for the random 
c    number generator.
c
      include 'omp_lib.h'

      integer n

      integer i
      integer my_id
      integer my_seed
      integer seed
      double precision x(n)

c$omp master
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Thread   Seed  I   X(I)'
      write ( *, '(a)' ) ' '
c$omp end master

c$omp parallel private ( i, my_id, my_seed ) shared ( n, x )
      my_id = omp_get_thread_num ( )
      my_seed = seed + my_id
      write ( *, '(2x,i6,2x,i12)' ) my_id, my_seed
c$omp do
      do i = 1, n
        call random_value ( my_seed, x(i) )
        write ( * , '(2x,i6,2x,i12,2x,i6,2x,g14.6)' ) 
     &    my_id, my_seed, i, x(i)
      end do
c$omp end do

c$omp end parallel

      return
      end
      subroutine random_value ( seed, r )

c*********************************************************************72
c
cc RANDOM_VALUE generates a random value R.
c
c  Discussion:
c
c    This is not a good random number generator.  It is a SIMPLE one.
c    It illustrates a model which works by accepting an integer seed value
c    as input, performing some simple operation on the seed, and then
c    producing a "random" real value using some simple transformation.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer SEED, a seed for the random 
c    number generator.
c
c    Output, double precision R, the random value.
c
      implicit none

      double precision r
      integer seed

      seed = mod ( seed, 65536 )
      seed = mod ( ( 3125 * seed ), 65536 )
      r = dble ( seed ) / 65536.0D+00

      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
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
