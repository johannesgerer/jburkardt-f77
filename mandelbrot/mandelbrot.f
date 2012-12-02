      program main

c*********************************************************************72
c
cc MAIN is the main program for MANDELBROT.
c
c  Discussion:
c
c    MANDELBROT computes an image of the Mandelbrot set.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Local Parameters:
c
c    Local, integer COUNT_MAX, the maximum number of iterations taken
c    for a particular pixel.
c
      implicit none

      integer n
      parameter ( n = 501 )

      integer b(n,n)
      integer c
      integer c_max
      integer count(n,n)
      integer count_max
      parameter ( count_max = 400 )
      character * ( 80 ) filename
      integer g(n,n)
      integer i
      integer ierror
      integer j
      integer k
      integer r(n,n)
      double precision x
      double precision x_max
      double precision x_min
      double precision x1
      double precision x2
      double precision y
      double precision y_max
      double precision y_min
      double precision y1
      double precision y2

      filename = 'mandelbrot.ppm'
      x_max =   1.25D+00
      x_min = - 2.25D+00
      y_max =   1.75D+00
      y_min = - 1.75D+00

      write ( *, '(a)' ) ' '
      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MANDELBROT'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) 
     &  '  Create an ASCII PPM image of the Mandelbrot set.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  For each point C = X + i*Y'
      write ( *, '(a,g14.6,a,g14.6,a)' ) 
     &  '  with X range [', x_min, ',', x_max, ']'
      write ( *, '(a,g14.6,a,g14.6,a)' ) 
     &  '  and  Y range [', y_min, ',', y_max, ']'
      write ( *, '(a,i8,a)' ) 
     &  '  carry out ', count_max, ' iterations of the map'
      write ( *, '(a)' ) '  Z(n+1) = Z(n)^2 + C.'
      write ( *, '(a)' ) 
     &  '  If the iterates stay bounded (norm less than 2)'
      write ( *, '(a)' ) '  then C is taken to be a member of the set.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  An ASCII PPM image of the set is created using'
      write ( *, '(a,i8,a)' ) 
     &  '    N = ', n, ' pixels in the X direction and'
      write ( *, '(a,i8,a)' ) 
     &  '    N = ', n, ' pixels in the Y direction.'
c
c  Carry out the iteration for each pixel, determining COUNT.
c
      do i = 1, n
        do j = 1, n

          x = ( dble (     j - 1 ) * x_max   
     &        + dble ( n - j     ) * x_min ) 
     &        / dble ( n     - 1 )

          y = ( dble (     i - 1 ) * y_max   
     &        + dble ( n - i     ) * y_min ) 
     &        / dble ( n     - 1 )

          count(i,j) = 0

          x1 = x
          y1 = y

          do k = 1, count_max

            x2 = x1 * x1 - y1 * y1 + x
            y2 = 2 * x1 * y1 + y

            if ( x2 .lt. -2.0D+00 .or. 2.0D+00 .lt. x2 .or. 
     &           y2 .lt. -2.0D+00 .or. 2.0D+00 .lt. y2 ) then
              count(i,j) = k
              go to 10
            end if

            x1 = x2
            y1 = y2

          end do

10        continue

        end do
      end do
c
c  Determine the coloring of each pixel.
c
      c_max = 0
      do i = 1, n
        do j = 1, n
          c_max = max ( c_max, count(i,j) )
        end do
      end do

      do i = 1, n
        do j = 1, n
          if ( mod ( count(i,j), 2 ) .eq. 1 ) then
            r(i,j) = 255
            g(i,j) = 255
            b(i,j) = 255
          else
            c = int ( 255.0D+00 * sqrt ( sqrt ( sqrt ( 
     &        ( dble ( count(i,j) ) / dble ( c_max ) ) ) ) ) )
            r(i,j) = 3 * c / 5
            g(i,j) = 3 * c / 5
            b(i,j) = c
          end if

        end do
      end do
c
c  Write an image file.
c
      filename = 'mandelbrot.ppm'

      call ppma_write ( filename, n, n, r, g, b, ierror )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  ASCII PPM image data stored in "' 
     &  // trim ( filename ) // '".'
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MANDELBROT'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine get_unit ( unit )

c*********************************************************************72
c
cc GET_UNIT returns a free FORTRAN unit number.
c
c  Discussion:
c
c    A "free" FORTRAN unit number is a value between 1 and 99 which
c    is not currently associated with an I/O device.  A free FORTRAN unit
c    number is needed in order to open a file with the OPEN command.
c
c    If UNIT = 0, then no free FORTRAN unit could be found, although
c    all 99 units were checked (except for units 5, 6 and 9, which
c    are commonly reserved for console I/O).
c
c    Otherwise, UNIT is a value between 1 and 99, representing a
c    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
c    are special, and will never return those values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 October 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer UNIT, the free unit number.
c
      implicit none

      integer i
      integer unit

      unit = 0

      do i = 1, 99

        if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

          open ( unit = i, err = 10, status = 'scratch' )
          close ( unit = i )

          unit = i

          return
        end if

10      continue

      end do

      return
      end
      subroutine ppma_write ( file_out_name, row_num, col_num, r, g, 
     &  b, ierror )

c*********************************************************************72
c
cc PPMA_WRITE writes an ASCII portable pixel map file.
c
c  Discussion:
c
c    PPM files can be viewed by XV.
c
c    Programs to convert files to this format include:
c
c      GIFTOPPM  - GIF file
c      PGMTOPPM  - Portable Gray Map file
c      PICTTOPPM - Macintosh PICT file
c      XPMTOPPM  - X11 pixmap file
c
c    Various programs can convert other formats to PPM format, including:
c
c      BMPTOPPM - Microsoft Windows BMP file.
c
c    A PPM file can also be converted to other formats, by programs:
c
c      PPMTOACAD - AutoCAD file
c      PPMTOGIF  - GIF file
c      PPMTOPGM  - Portable Gray Map file
c      PPMTOPICT - Macintosh PICT file
c      PPMTOPUZZ - X11 puzzle file
c      PPMTORGB3 - 3 Portable Gray Map files
c      PPMTOXPM  - X11 pixmap file
c      PPMTOYUV  - Abekas YUV file
c    
c  Example:
c
c    P3
c    # feep.ppma created by PBMLIB(PPMA_WRITE).
c    4 4
c    15
c     0  0  0    0  0  0    0  0  0   15  0 15
c     0  0  0    0 15  7    0  0  0    0  0  0
c     0  0  0    0  0  0    0 15  7    0  0  0
c    15  0 15    0  0  0    0  0  0    0  0  0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 February 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) FILE_OUT_NAME, the name of the file 
c    to which the data should be written.
c
c    Input, integer ROW_NUM, COL_NUM, the number of rows 
c    and columns of data.
c
c    Input, integer R(ROW_NUM,COL_NUM), G(ROW_NUM,COL_NUM),
c    B(ROW_NUM,COL_NUM), the red, green and blue values of each pixel.  These 
c    should be positive.
c
c    Output, integer IERROR, an error flag.
c    0, no error.
c    1, the data was illegal.
c    2, the file could not be opened.
c
      implicit none

      integer col_num
      integer row_num

      integer b(row_num,col_num)
      character ( len = * ) file_out_name
      integer file_out_unit
      integer g(row_num,col_num)
      integer i
      integer ierror
      integer ios
      integer j
      integer r(row_num,col_num)
      integer rgb_max

      ierror = 0
c
c  Compute the maximum color value.
c
      rgb_max = 0
      do j = 1, col_num
        do i = 1, row_num
          rgb_max = max ( rgb_max, r(i,j) )
          rgb_max = max ( rgb_max, g(i,j) )
          rgb_max = max ( rgb_max, b(i,j) )
        end do
      end do
c
c  Open the file.
c
      call get_unit ( file_out_unit )

      open ( unit = file_out_unit, file = file_out_name, 
     &  status = 'replace', form = 'formatted', access = 'sequential', 
     &  iostat = ios )

      if ( ios .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PPMA_WRITE - Fatal error!'
        write ( *, '(a)' ) '  Could not open the file.'
        ierror = 2
        return
      end if
c
c  Write the header.
c
      call ppma_write_header ( file_out_name, file_out_unit, row_num, 
     &  col_num, rgb_max, ierror )
c
c  Write the data.
c
      call ppma_write_data ( file_out_unit, row_num, col_num, r, 
     &  g, b, ierror )
c
c  Close the file.
c
      close ( unit = file_out_unit )

      return
      end
      subroutine ppma_write_data ( file_out_unit, row_num, col_num, 
     &  r, g, b, ierror )

c*********************************************************************72
c
cc PPMA_WRITE_DATA writes the data of a PPMA file.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer FILE_OUT_UNIT, the output file unit number.
c
c    Input, integer ROW_NUM, COL_NUM, the number of rows 
c    and columns of data.
c
c    Input, integer R(ROW_NUM,COL_NUM), G(ROW_NUM,COL_NUM), 
c    B(ROW_NUM,COL_NUM), the red, green and blue values of each pixel.  These
c    should be positive.
c
c    Output, integer IERROR, an error flag.
c    0, no error.
c    1, the data was illegal.
c    2, the file could not be opened.
c
      implicit none

      integer col_num
      integer row_num

      integer b(row_num,col_num)
      integer file_out_unit
      integer g(row_num,col_num)
      integer i
      integer ierror
      integer j
      integer jhi
      integer jlo
      integer r(row_num,col_num)

      ierror = 0
c
c  Write the header.
c
      do i = 1, row_num
        do jlo = 1, col_num, 4
          jhi = min ( jlo + 3, col_num )
          write ( file_out_unit, '(12i5)' ) 
     &      ( r(i,j), g(i,j), b(i,j), j = jlo, jhi )
        end do
      end do

      return
      end
      subroutine ppma_write_header ( file_out_name, file_out_unit, 
     &  row_num, col_num, rgb_max, ierror )

c*********************************************************************72
c
cc PPMA_WRITE_HEADER writes the header of a PPMA file.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) FILE_OUT_NAME, the name of the output file.
c
c    Input, integer FILE_OUT_UNIT, the output file unit number.
c
c    Input, integer ROW_NUM, COL_NUM, the number of rows and 
c    columns of data.
c
c    Input, integer RGB_MAX, the maximum value of any data component.
c
c    Output, integer IERROR, an error flag.
c    0, no error.
c    1, the data was illegal.
c    2, the file could not be opened.
c
      implicit none

      character * ( * )  file_out_name
      integer file_out_unit
      integer ierror
      integer col_num
      integer row_num
      integer rgb_max

      ierror = 0
c
c  Write the header.
c
      write ( file_out_unit, '(a2)' ) 'P3'
      write ( file_out_unit, '(a)' ) '# ' // trim ( file_out_name ) 
     &  // ' created by PPMA_WRITE.F.'
      write ( file_out_unit, '(i5,2x,i5)' ) col_num, row_num
      write ( file_out_unit, '(i5)' ) rgb_max

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
