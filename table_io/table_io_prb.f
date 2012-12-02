      program main

c*********************************************************************72
c
cc MAIN is the main program for TABLE_IO_PRB.
c
c  Discussion:
c
c    TABLE_IO_PRB calls the TABLE_IO test routines.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 August 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TABLE_IO_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version:'
      write ( *, '(a)' ) '  Test the TABLE_IO library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
      call test05 ( )
      call test06 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TABLE_IO_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests R8MAT_WRITE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 August 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 20 )
      integer m
       parameter ( m = 5 )

      integer i
      integer j
      character * ( 80 ) output_filename
      double precision table(m,n)

      output_filename = 'r8mat_05_00020.txt'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) 
     &  '  R8MAT_WRITE writes an R8MAT file.'

      do i = 1, m
        do j = 1, n
          table(i,j) = dble ( 100 * j + i ) / 10.0D+00
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Spatial dimension M = ', m
      write ( *, '(a,i8)' ) '  Number of points N  = ', n

      call r8mat_print_some ( m, n, table, 1, 1, 5, 5, 
     &  '  5x5 portion of the data written to file:' )

      call r8mat_transpose_print_some ( m, n, table, 1, 1, 5, 5, 
     &  '  5x5 portion of the TRANSPOSED data:' )

      call r8mat_write ( output_filename, m, n, table )

      write ( *, '(a)' ) ' '
      write ( *, '(a,a,a)' ) '  Wrote the file "',
     &  output_filename, '".'

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests R8MAT_HEADER_READ, R8MAT_DATA_READ.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 August 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      character * ( 80 ) input_filename
      integer m
      integer n
      double precision table(5,20)

      input_filename = 'r8mat_05_00020.txt'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  For an R8MAT file,'
      write ( *, '(a)' ) '  R8MAT_HEADER_READ reads the header'
      write ( *, '(a)' ) '  (about the dimensions of the data);'
      write ( *, '(a)' ) '  R8MAT_DATA_READ reads the data.'

      call r8mat_header_read (  input_filename, m, n )

      write ( *, '(a)' ) ' '
      write ( *, '(a,a,a)' ) 
     &  '  Read the header of "', input_filename, '".'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Spatial dimension M = ', m
      write ( *, '(a,i8)' ) '  Number of points N  = ', n

      call r8mat_data_read ( input_filename, m, n, table )

      write ( *, '(a)' ) ' '
      write ( *, '(a,a,a)' )
     &  '  Read the data in "', input_filename, '".'

      call r8mat_print_some ( m, n, table, 1, 1, 5, 5, 
     &  '  5x5 portion of data read from file:' )

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests I4MAT_WRITE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 August 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 20 )
      integer m
      parameter ( m = 5 )

      integer i
      integer j
      character ( len = 80 ) output_filename
      integer table(m,n)

      output_filename = 'i4mat_05_00020.txt'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  I4MAT_WRITE writes an I4MAT file.'

      do i = 1, m
        do j = 1, n
          table(i,j) = 100 * j + i
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Spatial dimension M = ', m
      write ( *, '(a,i8)' ) '  Number of points N  = ', n

      call i4mat_print_some ( m, n, table, 1, 1, 5, 5, 
     &  '  5 x 5 portion of data written to file:' )

      call i4mat_write ( output_filename, m, n, table )

      write ( *, '(a)' ) ' '
      write ( *, '(a,a,a)' ) '  Wrote the header and data for "',
     &  output_filename, '".'

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests I4MAT_HEADER_READ, I4MAT_DATA_READ.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 August 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      character * ( 80 ) input_filename
      integer m
      integer n
      integer table(20,5)

      input_filename = 'i4mat_05_00020.txt'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) '  Foran I4MAT file,'
      write ( *, '(a)' ) '  I4MAT_HEADER_READ reads the header'
      write ( *, '(a)' ) '  (about the dimensions of the data);'
      write ( *, '(a)' ) '  I4MAT_DATA_READ reads the data.'

      call i4mat_header_read (  input_filename, m, n )

      write ( *, '(a)' ) ' '
      write ( *, '(a,a,a)' ) 
     &  '  Read the header of "', input_filename, '".'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Spatial dimension M = ', m
      write ( *, '(a,i8)' ) '  Number of points N  = ', n

      call i4mat_data_read (  input_filename, m, n, table )

      write ( *, '(a)' ) ' '
      write ( *, '(a,a,a)' ) 
     &  '  Read the data in "', input_filename, '".'

      call i4mat_print_some ( m, n, table, 1, 1, 5, 5, 
     &  '  5 x 5 portion of data read from file:' )

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 tests R8MAT_UNIFORM_01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 August 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 2 )
      integer n
      parameter ( n = 10 )

      integer seed
      double precision table(m,n)

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) '  R8MAT_UNIFORM_01 sets a random R8MAT.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Spatial dimension M = ', m
      write ( *, '(a,i8)' ) '  Number of points N  = ', n

      call r8mat_uniform_01 ( m, n, seed, table )

      call r8mat_print_some ( m, n, table, 1, 1, 5, 10, 
     &  '  5x10 portion of random real table dataset:' )

      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 tests I4MAT_BORDER_ADD, I4MAT_BORDER_CUT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 August 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 6 )
      integer n 
      parameter ( n = 4 )

      integer table(m,n)
      integer table2(m-2,n-2)
      integer table3(m,n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' ) '  I4MAT_BORDER_CUT cuts off the border;'
      write ( *, '(a)' ) '  I4MAT_BORDER_ADD adds a zero border.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Spatial dimension M = ', m
      write ( *, '(a,i8)' ) '  Number of points N  = ', n

      call i4mat_indicator ( m, n, table )

      call i4mat_print ( m, n, table, '  Initial dataset:' )

      call i4mat_border_cut ( m, n, table, table2 )

      call i4mat_print ( m-2, n-2, table2, '  "Cut" dataset:' )

      call i4mat_border_add ( m-2, n-2, table2, table3 )

      call i4mat_print ( m, n, table3, '  "Added" dataset:' )

      return
      end
