      program main

c*********************************************************************72
c
cc MAIN is the main program for IMAGE_DENOISE_PRB.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'IMAGE_DENOISE_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the IMAGE_DENOISE library.'

      call test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'IMAGE_DENOISE_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      return
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests GRAY_MEDIAN_NEWS.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer g(428,320)
      integer g_max
      integer g2(428,320)
      character * ( 80 ) input_filename
      integer input_unit
      integer ios
      integer m
      integer n
      character * ( 80 ) output_filename

      input_filename = 'glassware_noisy.ascii.pgm'
      output_filename = 'glassware_median_news.ascii.pgm'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' ) '  GRAY_MEDIAN_NEWS uses a NEWS median filter'
      write ( *, '(a)' ) '  on a noisy grayscale image.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  The input file is "' // trim ( input_filename ) // '".'
c
c  Open the input file and read the data.
c
      call get_unit ( input_unit )

      open ( unit = input_unit, file = input_filename, status = 'old', 
     &  iostat = ios )

      if ( ios .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST01 - Fatal error!'
        write ( *, '(a)' ) 
     &    '  Could not open "' // trim ( input_filename ) // '".'
        stop
      end if

      call pgma_read_header ( input_unit, m, n, g_max )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  Number of rows =          ', m
      write ( *, '(a,i4)' ) '  Number of columns =       ', n
      write ( *, '(a,i4)' ) '  Maximum pixel intensity = ', g_max

      call pgma_read_data ( input_unit, m, n, g )

      close ( unit = input_unit )
 
      call gray_median_news ( m, n, g, g2 )
c
c  Write the denoised images.
c
      call pgma_write ( output_filename, m, n, g2 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Wrote denoised image to "' 
     &  // trim ( output_filename ) // '".'

      return
      end
