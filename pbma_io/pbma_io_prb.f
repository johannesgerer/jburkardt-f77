      program main

c*********************************************************************72
c
cc MAIN is the main program for PBMA_IO_PRB.
c
c  Discussion:
c
c    PBMA_IO_PRB calls the PBMA_IO test routines.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 June 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PBMA_IO_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the PBMA_IO library.'

      call test01 ( )
      call test02 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PBMA_IO_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests PBMA_EXAMPLE and PBMA_WRITE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 June 2010
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
      parameter ( file_name = 'pbma_io_prb_01.ascii.pbm' )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  PBMA_EXAMPLE sets up ASCII PBM data.'
      write ( *, '(a)' ) '  PBMA_WRITE writes an ASCII PBM file.'

      call pbma_example ( nrow, ncol, b )

      call pbma_write ( file_name, nrow, ncol, b )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Wrote the header and data for "'
     &  // trim ( file_name ) //'".'
      write ( *, '(a,i8)' ) '  Number of rows of data =    ', nrow
      write ( *, '(a,i8)' ) '  Number of columns of data = ', ncol

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests PBMA_READ.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 June 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m_max
      parameter ( m_max = 300 )
      integer n_max
      parameter ( n_max = 300 )

      integer b(m_max,n_max)
      character * ( 80 ) file_name
      parameter ( file_name = 'pbma_io_prb_02.ascii.pbm' )
      integer file_unit
      integer i
      integer ierror
      integer ios
      integer j
      integer k
      integer ncol
      integer nrow

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  PBMA_READ reads an ASCII PBM file.'

      call pbma_write_test ( file_name )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  PBMA_WRITE_TEST created some data.'

      call get_unit ( file_unit )

      open ( unit = file_unit, file = file_name, status = 'old',
     &  err = 10 )

      call pbma_read_header ( file_unit, nrow, ncol )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  PBMA_READ_HEADER has read the header.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of rows of data =    ', nrow
      write ( *, '(a,i8)' ) '  Number of columns of data = ', ncol

      if ( m_max * n_max .lt. nrow * ncol ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST02 - Error!'
        write ( *, '(a)' ) '  Internal array not large enough.'
        return
      end if

      call pbma_read_data ( file_unit, nrow, ncol, b )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  PBMA_READ_DATA has read the data.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Sample data:'
      write ( *, '(a)' ) ' '
      do k = 1, 30
        i = ( ( 30 - k ) * 1 + ( k - 1 ) * nrow ) / ( 30 - 1 )
        j = ( ( 30 - k ) * 1 + ( k - 1 ) * ncol ) / ( 30 - 1 )
        write ( *, '(2x,i8,2x,i8,2x,i8)' ) i, j, b(i,j)
      end do

      call pbma_check_data ( nrow, ncol, b )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The data was accepted by PBMA_CHECK_DATA.'

      return

10    continue

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02 - Error!'
      write ( *, '(a)' ) '  Could not open the file.'

      return
      end
