      program main

c*********************************************************************72
c
cc TOMS467_PRB tests XPOSE.
c
c  Modified:
c
c    12 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer a_max
      integer moved_max

      parameter ( a_max = 3000 )
      parameter ( moved_max = 100 )

      real a(a_max)
      logical moved(moved_max)
      integer n1
      integer n12
      integer n2
      integer nwork

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS467_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS467 library.'

      n1 = 10
      n2 = 10
      n12 = n1 * n2
      nwork = ( n1 + n2 ) / 2

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Row dimension N1 =    ', n1
      write ( *, '(a,i6)' ) '  Column dimension N2 = ', n2
      write ( *, '(a,i6)' ) '  Total size N12 =      ', n12
      write ( *, '(a,i6)' ) '  Workspace NWORK =     ', nwork

      call set_a ( n1, n2, a )

      call print_a ( n1, n2, a, 1, 5, 1, 5 )

      call xpose ( a, n1, n2, n12, moved, nwork )

      call print_a ( n2, n1, a, 1, 5, 1, 5 )

      n1 = 7
      n2 = 30
      n12 = n1 * n2
      nwork = ( n1 + n2 ) / 2

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Row dimension N1 =    ', n1
      write ( *, '(a,i6)' ) '  Column dimension N2 = ', n2
      write ( *, '(a,i6)' ) '  Total size N12 =      ', n12
      write ( *, '(a,i6)' ) '  Workspace NWORK =     ', nwork

      call set_a ( n1, n2, a )

      call print_a ( n1, n2, a, 1, 5, 1, 5 )

      call xpose ( a, n1, n2, n12, moved, nwork )

      call print_a ( n2, n1, a, 1, 5, 1, 5 )

      n1 = 24
      n2 = 8
      n12 = n1 * n2
      nwork = ( n1 + n2 ) / 2

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Row dimension N1 =    ', n1
      write ( *, '(a,i6)' ) '  Column dimension N2 = ', n2
      write ( *, '(a,i6)' ) '  Total size N12 =      ', n12
      write ( *, '(a,i6)' ) '  Workspace NWORK =     ', nwork

      call set_a ( n1, n2, a )

      call print_a ( n1, n2, a, 1, 5, 1, 5 )

      call xpose ( a, n1, n2, n12, moved, nwork )

      call print_a ( n2, n1, a, 1, 5, 1, 5 )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS467_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine set_a ( n1, n2, a )

c*********************************************************************72
c
cc SET_A sets the matrix A.
c
      implicit none

      integer n1
      integer n2

      real a(n1,n2)
      integer i
      integer j

      do i = 1, n1
        do j = 1, n2
          a(i,j) = 1000 * i + j
        end do
      end do

      return
      end
      subroutine print_a ( m, n, a, i_lo, i_hi, j_lo, j_hi )

c*********************************************************************72
c
cc PRINT_A prints the matrix A.
c
      implicit none

      integer m
      integer n

      real a(m,n)
      integer i
      integer i_hi
      integer i_lo
      integer j
      integer j_hi
      integer j_lo

      write ( *, '(a)' ) ' '

      do i = i_lo, i_hi
        write ( *, '(2x,5f8.0)' ) ( a(i,j), j = j_lo, j_hi )
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
