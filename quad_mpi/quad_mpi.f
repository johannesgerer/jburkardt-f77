      program main

c*********************************************************************72
c
cc MAIN is the main program for QUAD_MPI.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      include 'mpif.h'

      double precision a
      double precision b
      double precision error
      integer error_flag
      double precision exact
      external f
      double precision f
      integer i
      integer master
      double precision my_a
      double precision my_b
      integer my_id
      integer my_n
      double precision my_total
      integer p
      integer p_num
      integer n
      integer source
      integer status(MPI_STATUS_SIZE)
      integer tag
      integer target
      double precision total
      double precision wtime
      double precision x

      a =  0.0
      b = 10.0
      n = 10000000
      exact = 0.49936338107645674464D+00

      master = 0

      call MPI_Init ( error_flag )

      call MPI_Comm_size ( MPI_COMM_WORLD, p_num, error_flag )

      call MPI_Comm_rank ( MPI_COMM_WORLD, my_id, error_flag )
c
c  Process 0 reads in the quadrature rule, and parcels out the
c  evaluation points among the processes.
c
      if ( my_id .eq. 0 ) then
c
c  We want N to be the total number of evaluations.
c  If necessary, we adjust N to be divisible by the number of processors.
c
        my_n = n / ( p_num - 1 )
        n = ( p_num - 1 ) * my_n

        wtime = MPI_Wtime ( )

        call timestamp ( )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'QUAD_MPI'
        write ( *, '(a)' ) '  FORTRAN77/MPI version'
        write ( *, '(a)' ) '  Estimate an integral of f(x) from A to B.'
        write ( *, '(a)' ) '  f(x) = 50 / (pi * ( 2500 * x * x + 1 ) )'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) ' '
        write ( *, '(a,g14.6)' ) '  A        = ', a
        write ( *, '(a,g14.6)' ) '  B        = ', b
        write ( *, '(a,i12)' ) '  N        = ', n
        write ( *, '(a,g24.16)' ) '  Exact    = ', exact
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Use MPI to divide the computation among'
        write ( *, '(a)' ) '  multiple processes.'

      end if

      source = master

      call MPI_Bcast ( my_n, 1, MPI_INTEGER, source, MPI_COMM_WORLD, 
     &  error_flag )
c
c  Process 0 assigns each process a subinterval of [A,B].
c
      if ( my_id .eq. 0 ) then

        do p = 1, p_num - 1

          my_a = ( dble ( p_num - p     ) * a   
     &           + dble (         p - 1 ) * b ) 
     &           / dble ( p_num     - 1 )

          target = p
          tag = 1
          call MPI_Send ( my_a, 1, MPI_DOUBLE_PRECISION, target, tag, 
     &      MPI_COMM_WORLD, error_flag )

          my_b = ( dble ( p_num - p - 1 ) * a   
     &           + dble (         p     ) * b ) 
     &           / dble ( p_num     - 1 )

          target = p
          tag = 2
          call MPI_Send ( my_b, 1, MPI_DOUBLE_PRECISION, target, tag, 
     &      MPI_COMM_WORLD, error_flag )

        end do

        total = 0.0D+00
        my_total = 0.0D+00
c
c  Processes receive MY_A, MY_B, and compute their part of the integral.
c
      else

        source = master
        tag = 1

        call MPI_Recv ( my_a, 1, MPI_DOUBLE_PRECISION, source, tag, 
     &    MPI_COMM_WORLD, status, error_flag )

        source = master
        tag = 2

        call MPI_Recv ( my_b, 1, MPI_DOUBLE_PRECISION, source, tag, 
     &    MPI_COMM_WORLD, status, error_flag )

        my_total = 0.0D+00
        do i = 1, my_n
          x = ( ( my_n - i ) * my_a + ( i - 1 ) * my_b ) / ( my_n - 1 )
          my_total = my_total + f ( x )
        end do

        my_total = ( my_b - my_a ) * my_total / my_n

        write ( *, '(a,i8,a,g14.6)' ) 
     &    '  Process ', my_id, ' contributes MY_TOTAL = ', my_total

      end if
c
c  Each process sends its value of MY_TOTAL to the master process, to
c  be summed in TOTAL.
c
      call MPI_Reduce ( my_total, total, 1, MPI_DOUBLE_PRECISION, 
     &  MPI_SUM, master, MPI_COMM_WORLD, error_flag )
c
c  Compute the weighted estimate.
c
      if ( my_id .eq. master ) then
        error = abs ( total - exact )
        wtime = MPI_Wtime ( ) - wtime
        write ( *, '(a)' ) ' '
        write ( *, '(a,g24.16)' ) '  Estimate = ', total
        write ( *, '(a,g14.6)' ) '  Error    = ', error
        write ( *, '(a,g14.6)' ) '  Time     = ', wtime
      end if
c
c  Terminate MPI.
c
      call MPI_Finalize ( error_flag )
c
c  Terminate.
c
      if ( my_id .eq. master ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'QUAD_MPI:'
        write ( *, '(a)' ) '  Normal end of execution.'
        write ( *, '(a)' ) ' '
        call timestamp ( ) 
      end if

      stop
      end
      function f ( x )

c*********************************************************************72
c
cc F evaluates the function.
c
      double precision f
      double precision pi
      double precision x

      pi = 3.141592653589793D+00
      f = 50.0D+00 / ( pi * ( 2500.0D+00 * x * x + 1.0D+00 ) )

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
