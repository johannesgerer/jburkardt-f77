      program main

c*********************************************************************72
c
cc MAIN is the main program for MATVEC.
c
c  Discussion:
c
c    MATVEC uses MPI to multiply a matrix times a vector.
c
c    The computation b = A * x is carried out by giving every process
c    a copy of the vector x.  Then each process is given one row of
c    the matrix A, say "A(I,*)", and computes the product of A times
c    x, returning the result, which is b(i).
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

      integer MAX_COLS
      integer MAX_ROWS

      parameter ( MAX_COLS = 1000 )
      parameter ( MAX_ROWS = 1000 )

      double precision a(MAX_ROWS,MAX_COLS)
      double precision ans
      integer anstype
      double precision b(MAX_ROWS)
      double precision buffer(MAX_COLS)
      integer i
      integer ierr
      integer j
      integer j_one
      integer m
      integer master
      integer my_id
      integer n
      integer num_procs
      integer numsent
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      integer row
      integer sender
      integer status(MPI_STATUS_SIZE)
      integer tag
      double precision x(MAX_COLS)
c
c  Initialize MPI.
c
      call MPI_INIT ( ierr )
c
c  Get this processor's ID.
c
      call MPI_COMM_RANK ( MPI_COMM_WORLD, my_id, ierr )
c
c  Get the number of processors.
c
      call MPI_COMM_SIZE ( MPI_COMM_WORLD, num_procs, ierr )

      master = 0
      m = 100
      n = 50
c
c  The master process initializes and then dispatches.
c
c  Initialize A and X.
c
      if ( my_id .eq. 0 ) then
        call timestamp ( )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MATVEC'
        write ( *, '(a)' ) '  FORTRAN77 version'
        write ( *, '(a)' ) '  A program using MPI to'
        write ( *, '(a)' ) '  compute a matrix vector product'
        write ( *, '(a)' ) '    b = A * x.'
        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' ) '  The number of processes is ', num_procs
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8,a)' ) '  Process ', my_id, ' is active.'

      call MPI_Bcast ( n, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
c
c  The master process initializes A and x.
c
      if ( my_id == 0 ) then

        do i = 1, m
          do j = 1, n
            a(i,j) = sqrt ( 2.0D+00 / dble ( n + 1 ) ) 
     &        * sin ( dble ( i * j ) * pi / dble ( n + 1 ) )
          end do
        end do
c
c  X is specially chosen so that b = A * x is known in advance.
c  The value of B will be zero, except that entry J_ONE will be 1.
c  Pick any value of J_ONE between 1 and M.
c
        j_one = 17

        do i = 1, n
          x(i) = sqrt ( 2.0D+00 / dble ( n + 1 ) ) 
     &      * sin ( dble ( i * j_one ) * pi / dble ( n + 1 ) )
        end do

      end if
c
c  Send X to each worker process.
c
      if ( my_id == 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Broadcasting vector X to all processes.'
      end if

      call MPI_BCAST ( x, n, MPI_DOUBLE_PRECISION, master, 
     &  MPI_COMM_WORLD, ierr )
c
c  Send a row to each worker process; tag with row number.
c
      if ( my_id == master ) then

        numsent = 0

        do i = 1, min ( num_procs - 1, m )
          do j = 1, n
            buffer(j) = a(i,j)
          end do
          tag = i
          call MPI_SEND ( buffer, n, MPI_DOUBLE_PRECISION, i, 
     &      tag, MPI_COMM_WORLD, ierr )
          numsent = numsent + 1
        end do
c
c  Wait to receive a result back from a processor, and immediately
c  send that processor another row to work with.
c
        do i = 1, m

          call MPI_RECV ( ans, 1, MPI_DOUBLE_PRECISION, 
     &       MPI_ANY_SOURCE, MPI_ANY_TAG, 
     &       MPI_COMM_WORLD, status, ierr )

          sender = status(MPI_SOURCE)    
          tag = status(MPI_TAG)
          b(tag) = ans

          if ( numsent .lt. m ) then

            numsent = numsent + 1

            do j = 1, n
              buffer(j) = a(numsent,j)
            end do

            tag = numsent
            call MPI_SEND ( buffer, n, MPI_DOUBLE_PRECISION, 
     &         sender, tag, MPI_COMM_WORLD, ierr )

          else

            tag =  m + 1
            call MPI_SEND ( MPI_BOTTOM, 0, MPI_DOUBLE_PRECISION, 
     &       sender, tag, MPI_COMM_WORLD, ierr )

          end if

        end do
c
c  Each worker process repeatedly receives rows of A (with TAG indicating 
c  which row it is), computes dot products A(I,1:N) * X(1:N) and returns 
c  the result (and TAG), until receiving the "DONE" message.
c
      else
c
c  Skip if more processes than work.
c
        if ( rank .gt. m ) then
          go to 200     
        end if

 90     continue

        call MPI_RECV ( buffer, n, MPI_DOUBLE_PRECISION, master, 
     &    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr )

        tag = status(MPI_TAG)

        if ( tag .eq. m + 1 ) then
          write ( *, '(a,i8,a)' ) '  Process ', my_id, ' shutting down.'
          go to 200
        end if

        ans = 0.0D+00
        do i = 1, n
          ans = ans + buffer(i) * x(i)
        end do

        call MPI_SEND ( ans, 1, MPI_DOUBLE_PRECISION, master, 
     &     tag, MPI_COMM_WORLD, ierr )

        go to 90

 200    continue

      end if
c
c  Print the answer.
c
      if ( my_id .eq. master ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MATVEC - Master process:'
        write ( *, '(a)' ) '  Product vector b = A * x:'
        write ( *, '(a,i8)' ) 
     &    '  (Should be zero, except for a 1 in entry ', j_one
        write ( *, '(a)' ) ' '
        do i = 1, m
          write ( *, '(i4,2x,f10.4)' ) i, b(i)
        end do
      end if
c
c  Terminate MPI.
c
      call MPI_FINALIZE ( ierr )
c
c  Terminate.
c
      if ( my_id .eq. master ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MATVEC - Master process:'
        write ( *, '(a)' ) '  Normal end of execution.'
        write ( *, '(a)' ) ' '
        call timestamp ( )
      end if

      stop
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
