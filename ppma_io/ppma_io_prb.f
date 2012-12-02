      program main

c*********************************************************************72
c
cc MAIN is the main program for PPMA_IO_PRB.
c
c  Discussion:
c
c    PPMA_IO_PRB calls the PPMA_IO test routines.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PPMA_IO_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the PPMA_IO library.'

      call test01 ( )
      call test02 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PPMA_IO_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests PPM_EXAMPLE and PPMA_WRITE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer ncol
      parameter ( ncol = 300 )
      integer nrow
      parameter ( nrow = 300 )

      integer b(nrow,ncol)
      character * ( 80 ) file_name
      parameter ( file_name = 'test01.ascii.ppm' )
      integer g(nrow,ncol)
      integer ierror
      integer r(nrow,ncol)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  PPMA_EXAMPLE sets up sample PPMA data.'
      write ( *, '(a)' ) '  PPMA_WRITE writes an ASCII PPMA file.'

      call ppma_example ( nrow, ncol, r, g, b )

      call ppma_write ( file_name, nrow, ncol, r, g, b, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST01 - Fatal errorc'
        write ( *, '(a,i6)' ) 'PPMA_WRITE returns IERROR = ', ierror
        return
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Wrote the header and data for "' 
     &  // trim ( file_name ) //'".'
      write ( *, '(a,i6)' ) '  Number of rows of data =    ', nrow
      write ( *, '(a,i6)' ) '  Number of columns of data = ', ncol

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests PPMA_READ_DATA and PPMA_READ_HEADER.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer b(300,300)
      character * ( 80 ) file_name
      parameter ( file_name = 'test02.ascii.ppm' )
      integer file_unit
      integer g(300,300)
      integer i
      integer ierror
      integer ios
      integer j
      integer k
      integer ncol
      integer nrow
      integer r(300,300)
      integer rgb_max

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  PPMA_READ_HEADER reads the header.'
      write ( *, '(a)' ) '  PPMA_READ_HEADER reads the data.'

      call ppma_write_test ( file_name )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  PPMA_WRITE_TEST created some data.'

      call get_unit ( file_unit )

      open ( unit = file_unit, file = file_name, status = 'old', 
     &  iostat = ios )

      if ( ios .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST02 - Fatal error!'
        write ( *, '(a)' ) '  Could not open the file.'
        return
      end if

      call ppma_read_header ( file_unit, nrow, ncol, rgb_max, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST02 - Fatal error!'
        write ( *, '(a)' ) '  Error while reading the header.'
        return
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  PPMA_READ_HEADER read the header.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Number of rows of data =    ', nrow
      write ( *, '(a,i6)' ) '  Number of columns of data = ', ncol
      write ( *, '(a,i6)' ) '  Maximum RGB value =         ', rgb_max

      call ppma_read_data ( file_unit, nrow, ncol, r, g, b, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST02 - Fatal error!'
        write ( *, '(a)' ) '  PPMA_READ_DATA failed.'
        return
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  PPMA_READ_DATA read the data.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Ten sample values:'
      write ( *, '(a)' ) ' '
      do k = 1, 10
        i = ( ( 10 - k ) * 1 + ( k - 1 ) * nrow ) / ( 10 - 1 )
        j = ( ( 10 - k ) * 1 + ( k - 1 ) * ncol ) / ( 10 - 1 )
        write ( *, '(2i4,4x,3i6)' ) i, j, r(i,j), g(i,j), b(i,j)
      end do

      call ppma_check_data ( nrow, ncol, rgb_max, r, g, b, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST02 - Error!'
        write ( *, '(a,i6)' ) 
     &    '  The data was not accepted by PPMA_CHECK_DATA.'
        return
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The data was accepted by PPMA_CHECK_DATA.'

      return
      end
