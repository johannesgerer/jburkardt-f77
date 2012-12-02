      program main

c*********************************************************************72
c
cc MAIN is the main program for RING_MPI.
c
c  Discussion:
c
c    RING_MPI sends messages of various size from process 0 to 1 to 2 to
c    ...the last process and then back to 0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 October 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Peter Pacheco,
c    Parallel Programming with MPI,
c    Morgan Kaufman, 1996,
c    ISBN: 1558603395,
c    LC: QA76.642.P3.
c
      implicit none

      include 'mpif.h'

      integer error
      integer id
      integer p
c
c  Initialize MPI.
c
      call MPI_Init ( error )
c
c  Get the number of processes.
c
      call MPI_Comm_size ( MPI_COMM_WORLD, p, error )
c
c  Get the individual process ID.
c
      call MPI_Comm_rank ( MPI_COMM_WORLD, id, error )
c
c  Print a message.
c
      if ( id .eq. 0 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RING_MPI:'
        write ( *, '(a)' ) '  FORTRAN77/MPI version'
        write ( *, '(a)' ) '  Measure time required to transmit data'
        write ( *, '(a)' ) '  around a ring of processes'
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  The number of processes is ', p

      end if

      call ring_io ( p, id )
c
c  Shut down MPI.
c
      call MPI_Finalize ( error )
c
c  Terminate.
c
      if ( id .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RING_MPI:'
        write ( *, '(a)' ) '  Normal end of execution.'
      end if

      stop
      end
      subroutine ring_io ( p, id )

c*********************************************************************72
c
cc RING_IO carries out the tasks of process ID, of a total of P processes.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 October 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Peter Pacheco,
c    Parallel Programming with MPI,
c    Morgan Kaufman, 1996,
c    ISBN: 1558603395,
c    LC: QA76.642.P3.
c
      include 'mpif.h'

      integer n_test_num
      parameter ( n_test_num = 5 )

      integer n_max
      parameter ( n_max = 1000000 )

      integer dest
      integer error
      integer i
      integer id
      integer j
      integer n
      integer n_test(n_test_num)
      integer p
      integer source
      integer status(MPI_STATUS_SIZE)
      double precision tave
      integer test
      integer test_num
      parameter ( test_num = 10 )
      double precision tmax
      double precision tmin
      double precision wtime
      double precision x(n_max)

      save n_test

      data n_test / 100, 1000, 10000, 100000, 1000000 /

        if ( id .eq. 0 ) then

          write ( *, '(a)' ) ' '
          write ( *, '(a,i6,a)' ) 
     &      '  Timings based on ', test_num, ' experiments'
          write ( *, '(a,i6,a)' ) 
     &      '  N double precision values were sent'
          write ( *, '(a)' ) 
     &      '  in a ring transmission starting and ending at process 0'
          write ( *, '(a,i6,a)' ) 
     &      '  and using a total of ', p, ' processes.'  
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 
     &      '         N           T min           T ave           T max'
          write ( *, '(a)' ) ' '

      end if
c
c  Choose message size.
c
      do i = 1, n_test_num
        
        n = n_test(i)
c
c  Process 0 sends very first message, 
c  then waits to receive the "echo" that has gone around the world.
c
        if ( id .eq. 0 ) then

          dest = 1
          source = p - 1

          tave = 0.0D+00
          tmin = 1.0D+30
          tmax = 0.0D+00

          do test = 1, test_num
c
c  Just in case, set the entries of X in a way that identifies
c  which iteration of the test is being carried out.
c
            do j = 1, n
              x(j) = dble ( test + j - 1 )
            end do

            wtime = MPI_Wtime ( )
            call MPI_Send ( x, n, MPI_DOUBLE_PRECISION, dest,   0, 
     &        MPI_COMM_WORLD,         error )
            call MPI_Recv ( x, n, MPI_DOUBLE_PRECISION, source, 0, 
     &        MPI_COMM_WORLD, status, error )
            wtime = MPI_Wtime ( ) - wtime
c
c  Record the time it took.
c
            tave =       tave + wtime
            tmin = min ( tmin,  wtime )
            tmax = max ( tmax,  wtime )

          end do

          tave = tave / dble ( test_num )

          write ( *, '(2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &      n, tmin, tave, tmax
c
c  Worker ID must receive first from ID-1, then send to ID+1.
c
        else

          source = id - 1
          dest = mod ( id + 1, p )
     
          do test = 1, test_num
            call MPI_Recv ( x, n, MPI_DOUBLE_PRECISION, source, 0, 
     &        MPI_COMM_WORLD, status, error )
            call MPI_Send ( x, n, MPI_DOUBLE_PRECISION, dest,   0, 
     &        MPI_COMM_WORLD,         error )
          end do

        end if

      end do

      return
      end
