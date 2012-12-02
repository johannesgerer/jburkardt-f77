      program main

c*********************************************************************72
c
cc MAIN is the main program for ASA144_PRB.
c
c  Discussion:
c
c    ASA144_PRB calls a set of problems for ASA144.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA144_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the ASA144 library.'

      call test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA144_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests RCONT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nrow
      parameter ( nrow = 5 )
      integer ncol
      parameter ( ncol = 5 )

      integer i
      integer ifault
      logical key
      integer matrix(nrow,ncol)
      integer ncolt(ncol)
      integer nrowt(nrow)
      integer nsubt(ncol)
      integer test
      integer test_num
      parameter ( test_num = 10 )

      data ncolt / 2, 2, 2, 2, 1 /
      data nrowt / 3, 2, 2, 1, 1 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  RCONT constructs a random matrix with'
      write ( *, '(a)' ) '  given row and column sums.'

      call i4vec_print ( nrow, nrowt, '  The rowsum vector:' )
      call i4vec_print ( ncol, ncolt, '  The columnsum vector: ' )

      key = .false.

      do test = 1, test_num

        call rcont ( nrow, ncol, nrowt, ncolt, nsubt, matrix,
     &    key, ifault )

        if ( ifault .ne. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a,i8)' ) '  RCONT returned IFAULT = ', ifault
          return
        end if

        call i4mat_print ( nrow, ncol, matrix,
     &    '  The rowcolsum matrix:' )

      end do

      return
      end
