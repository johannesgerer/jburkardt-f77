      program main

c*********************************************************************72
c
cc MAIN is the main program for IMAGE_EDGE_PRB.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 December 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m_dim
      parameter ( m_dim = 246 )
      integer n_dim
      parameter ( n_dim = 300 )

      integer e(m_dim,n_dim)
      integer g(m_dim,n_dim)
      integer g_histo(0:255)
      integer g_max
      integer i
      character * ( 255 ) input_filename
      parameter ( input_filename = 'coins.ascii.pgm' )
      integer input_unit
      integer m
      integer n
      character * ( 255 ) output_filename
      parameter ( output_filename = 'coin_edges.ascii.pbm' )
      integer output_unit

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'IMAGE_EDGE_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the IMAGE_EDGE library.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Demonstrate the NEWS stencil for edge'
      write ( *, '(a)' ) '  detection in images.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The input file is "'
     &  // trim ( input_filename ) // '"'
c
c  Open the input file and read the data.
c
      call get_unit ( input_unit )

      open ( unit = input_unit, file = input_filename, status = 'old',
     &  err = 10 )

      go to 20

10    continue

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'IMAGE_EDGE_PRB - Fatal error.'
      write ( *, '(a)' ) '  Could not open the input file.'
      stop

20    continue

      call pgma_read_header ( input_unit, m, n, g_max )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of rows =          ', m
      write ( *, '(a,i8)' ) '  Number of columns =       ', n
      write ( *, '(a,i8)' ) '  Maximum pixel intensity = ', g_max

      call pgma_read_data ( input_unit, m, n, g )

      close ( unit = input_unit )

      call i4mat_histogram ( m, n, g, 255, g_histo )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) ' Gray     Count'
      write ( *, '(a)' ) ' '
      do i = 0, 255
        write ( *, '(2x,i3,2x,i8)' ) i, g_histo(i)
      end do

      call news ( m, n, g, e )
c
c  Write the edge information as a portable BIT map (0/1).
c
      call pbma_write ( output_filename, m, n, e )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Wrote edge information to "'
     &  // trim ( output_filename ) // '".'
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'IMAGE_EDGE_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
