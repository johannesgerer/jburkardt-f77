      program main

c*********************************************************************72
c
cc LATINIZE_PRB tests the LATINIZE routines.
c
c  Discussion:
c
c    The dataset is presumed to be an M by N array of real numbers,
c    where M is the spatial dimension, and N is the number of sample points.
c
c    The dataset is presumed to be stored in a file, with N records,
c    one per each sample point.  (Comment records may be included,
c    which begin with '#'.)
c
c    The program reads the data file, "latinizes" the data, and writes
c    the latinized data to a new file.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LATINIZE_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the LATINIZE library.'

      call test01 ( 'cvt_02_00010.txt' )
      call test01 ( 'cvt_03_00007.txt' )
      call test01 ( 'cvt_03_00056.txt' )
      call test01 ( 'cvt_07_00100.txt' )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LATINIZE_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( input_filename )

c*********************************************************************72
c
cc TEST01 reads a datafile, latinizes it, and writes the new data out.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      character * ( * ) input_filename
      integer m
      integer n
      character * ( 255 ) output_filename
      double precision table(7,100)
c
c  Get the row and column dimensions of the dataset.
c
      call r8mat_header_read (  input_filename, m, n )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Read the header of "' // trim ( input_filename ) //'".'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Spatial dimension M = ', m
      write ( *, '(a,i6)' ) '  Number of points N  = ', n
c
c  Read the array.
c
      call r8mat_data_read ( input_filename, m, n, table )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Read the data in "' // trim ( input_filename ) //'".'

      call r8mat_transpose_print_some ( m, n, table, 1, 1, 5, 5, 
     &  '  Small portion of data read from file:' )
c
c  Latinize the array.
c
      call r8mat_latinize ( m, n, table )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Latinized the data.'
c
c  Print out a small sample of the array.
c
      call r8mat_transpose_print_some ( m, n, table, 1, 1, 5, 5, 
     &  '  Small portion of Latinized data:' )
c
c  Make up a name for the output file that is likely to be related
c  to the input file name, and unique.
c
      output_filename = input_filename
      call file_name_ext_swap ( output_filename, 'latin.txt' )
c
c  Write the latinized array to a file.
c
      call r8mat_write ( output_filename, m, n, table )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Wrote the latinized data to "' 
     &  // trim ( output_filename ) //'".'

      return
      end
