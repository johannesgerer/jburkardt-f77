      program main

c*********************************************************************72
c
cc MAIN is the main program for MATMAT.
c
c  Discussion:
c
c    MATMAT uses MPI to multiply two matrices.
c
c    The computation C = A * B is carried out by giving every process
c    a copy of the matrix A.  Then each process is given one column of
c    the matrix B, say "B(*,J)", and computes the product of A times
c    this column, returning the result, which is column J of C.
c
c    This code is "self scheduling", because as soon as a process finishes
c    working on one column, it goes back to the master process and requests
c    another one.  This is an appropriate way to schedule work when the
c    processes are likely to run at different speeds, or individual
c    tasks might take greatly varying times.  Neither problem is likely
c    to occur here, but this is meant as an illustration.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 October 2005
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    William Gropp, Ewing Lusk, Anthony Skjellum,
c    Using MPI: Portable Parallel Programming with the
c    Message-Passing Interface,
c    Second Edition,
c    MIT Press, 1999,
c    ISBN: 0262571323,
c    LC: QA76.642.G76.
c
c    Marc Snir, Steve Otto, Steven Huss-Lederman, David Walker, 
c    Jack Dongarra,
c    MPI: The Complete Reference,
c    Volume I: The MPI Core,
c    Second Edition,
c    MIT Press, 1998,
c    ISBN: 0-262-69216-3,
c     LC: QA76.642.M65.
c
      include 'mpif.h'
c
c implicit none
c
c  NA and MB must be equal.
c
      integer ma
      integer na
      integer mb
      integer nb

      parameter ( ma = 10 )
      parameter ( na = 20 )
      parameter ( mb = 20 )
      parameter ( nb = 25 )
c
c  All processes need to store a full copy of A.
c
      real a(ma,na)
c
c  Process 0 needs a full copy of B.
c  The other processes do not.
c  We should really do an allocatable array here.
c
      real b(mb,nb)
      real b_column(mb)
c
c  Process 0 needs a full copy of C.
c  The other processes do not.
c
      real c(ma,nb)
      real c_column(ma)
      integer dest
      integer i
      integer ierr
      integer j
      integer jhi
      integer master
      integer my_id
      integer num_procs
      integer num_received
      integer num_sent
      integer source
      integer status(mpi_status_size)
      integer tag
      real value

      master = 0
c
c  Initialize MPI.
c
      call MPI_Init ( ierr )
c
c  Get this process's ID.
c
      call MPI_Comm_rank ( mpi_comm_world, my_id, ierr )
c
c  Get the number of processes.
c
      call MPI_Comm_size (  MPI_COMM_WORLD, num_procs, ierr )

      if ( my_id .eq. master ) then
        call timestamp ( )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MATMAT'
        write ( *, '(a)' ) '  FORTRAN77 version'
        write ( *, '(a)' ) '  A program using MPI'
        write ( *, '(a)' ) '  to compute a matrix product C = A * B.'
        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' ) '  The number of processes is ', num_procs
        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' ) '  Number of rows    of matrix A = ', ma
        write ( *, '(a,i6)' ) '  Number of columns of matrix A = ', na
        write ( *, '(a,i6)' ) '  Number of rows    of matrix B = ', mb
        write ( *, '(a,i6)' ) '  Number of columns of matrix B = ', nb

        if ( num_procs .lt. 2 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'MATMAT - Fatal error!'
          write ( *, '(a)' ) '  Must have at least 2 processes!'
        end if

      end if

      if ( num_procs .lt. 2 ) then
        call MPI_Abort (  MPI_COMM_WORLD, 1, ierr )
        stop
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6,a)' ) '  Process ', my_id, ' is active.'
c
c  The master process initializes A and B.
c
      if ( my_id .eq. master ) then

        do i = 1, ma
          do j = 1, na
            if ( j .eq. i-1 ) then
              a(i,j) = -1.0E+00
            else if ( j .eq. i ) then
              a(i,j) = 2.0E+00
            else if ( j .eq. i+1 ) then
              a(i,j) = -1.0E+00
            else
              a(i,j) = 0.0E+00
            end if
          end do
        end do

        do i = 1, mb
          do j = 1, nb
            if ( i .le. j ) then
              b(i,j) = real ( i * ( nb + 1 - j ) ) / real ( nb + 1 )
            else
              b(i,j) = real ( j * ( mb + 1 - i ) ) / real ( mb + 1 )
            end if
          end do
        end do

      end if
c
c  We have rigged the game so that, if the matrices are square, then B is the 
c  inverse of A, and therefore the product C should be the identity matrix.
c
c  Just so we can do rectangular problems too, and still know we're on the
c  right track, compute the exact result here and print out a chunk of it.
c
      if ( my_id .eq. master ) then

        call mat_mul ( ma, na, nb, a, b, c )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MATMAT - Master process:'
        write ( *, '(a)' ) '  Initial 5 x 5 block '
        write ( *, '(a)' ) '  of exact product matrix C:'
        write ( *, '(a)' ) ' '

        do i = 1, 5
          write ( *, '(5g14.6)' ) ( c(i,j), j = 1, 5 )
        end do

      end if
c
c  The master process broadcasts a copy of the entire A matrix 
c  to all other processes.
c
      source = master

      call MPI_Bcast ( a, ma*nb, mpi_real, source, 
     &  MPI_COMM_WORLD, ierr )
c
c  Now the master process distributes columns of B to the other processes,
c  and collects the products C(1:MA,J) = A(1:MA,1:NA) * B(1:NA,J)
c
      if ( my_id .eq. master ) then

        num_sent = 0
c
c  Start by sending column J of B to process J.
c
        jhi = min ( num_procs-1, nb )

        do j = 1, jhi

          do i = 1, mb
            b_column(i) = b(i,j)
          end do

          dest = j
          tag = j

          call MPI_Send ( b_column, mb, MPI_REAL, dest, tag, 
     &      MPI_COMM_WORLD, ierr )

          num_sent = num_sent + 1

        end do
c
c  Process 0 waits to receive a result from any process, and sends that
c  process a new column, or a terminate signal if there are no more columns.
c
c  Once all columns have been received, process 0 exits from the loop.
c
        num_received = 0

95      continue

        if ( num_received .lt. nb ) then

          call MPI_Recv ( c_column, ma, MPI_REAL, 
     &      MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, 
     &      status, ierr )

          num_received = num_received + 1
          source = status(mpi_source)
          tag = status(mpi_tag)

          do i = 1, ma
            c(i,tag) = c_column(i)
          end do

          if ( num_sent .lt. nb ) then

            num_sent = num_sent + 1
            do i = 1, mb
              b_column(i) = b(i,num_sent)
            end do
            dest = source
            tag = num_sent

            call MPI_Send ( b_column, mb, MPI_REAL, dest, tag, 
     &        MPI_COMM_WORLD, ierr )

          else

            value = 1.0E+00
            dest = source
            tag = 0

            call MPI_Send ( value, 1,  MPI_REAL, dest, tag, 
     &        MPI_COMM_WORLD, ierr )

          end if

          go to 95

        end if
c
c  Each worker process receives a column of B.
c
      else if ( my_id .le. nb ) then

98      continue

          source = master

          call MPI_Recv ( b_column, mb, MPI_REAL, source, 
     &      MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr )
c
c  If the tag is 0, we're actually being told to terminate.
c
          tag = status(mpi_tag)

          if ( tag .eq. 0 ) then
            go to 99
          end if

          call mat_mul ( ma, na, 1, a, b_column, c_column )

          dest = master

          call MPI_Send ( c_column, ma, MPI_REAL, dest, tag, 
     &       MPI_COMM_WORLD, ierr )

        go to 98

99      continue

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Process ', my_id
        write ( *, '(a)' ) '  MPI has no work for mec'

      end if
c 
c  Process 0 prints out a bit of the answer.
c
      if ( my_id .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MATMAT - Master process:'
        write ( *, '(a)' ) '  Initial 5 x 5 block of '
        write ( *, '(a)' ) '  computed product matrix C:'
        write ( *, '(a)' ) ' '

        do i = 1, 5
          write ( *, '(5g14.6)' ) ( c(i,j), j = 1, 5 )
        end do

      end if
c
c  Terminate MPI.
c
      call MPI_Finalize ( ierr )
c
c  Terminate.
c
      if ( my_id .eq. master ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MATMAT - Master process:'
        write ( *, '(a)' ) '  Normal end of execution.'
        write ( *, '(a)' ) ' '
        call timestamp ( )
      end if

      stop
      end
      subroutine mat_mul ( n1, n2, n3, a, b, c )

c*********************************************************************72
c
cc MAT_MUL multiplies two matrices.
c
c  Discussion:
c
c    In FORTRAN90, a similar facility named MATMUL is available as a system library routine.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 October 2005
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, N3, dimensions of A, B and C.
c
c    Input, real A(N1,N2), B(N2,N3), the matrices to multiply.
c
c    Output, real C(N1,N3), the product matrix.
c
      integer n1
      integer n2
      integer n3

      real a(n1,n2)
      real b(n2,n3)
      real c(n1,n3)
      integer i
      integer j
      integer k

      do i = 1, n1
        do j = 1, n3
          c(i,j) = 0.0E+00
          do k = 1, n2
            c(i,j) = c(i,j) + a(i,k) * b(k,j)
          end do
        end do
      end do

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

