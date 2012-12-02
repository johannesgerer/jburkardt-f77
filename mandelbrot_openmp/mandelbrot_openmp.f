      program main

c*****************************************************************************80
c
cc MAIN is the main program for MANDELBROT_OPENMP.
c
c  Discussion:
c
c    MANDELBROT_OPENMP computes an image of the Mandelbrot set.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 July 2010
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

      include 'omp_lib.h'

      integer m
      parameter ( m = 500 )
      integer n
      parameter ( n = 500 )

      integer b(m,n)
      integer c
      integer c_max
      integer count(m,n)
      integer count_max
      integer g(m,n)
      integer i
      integer ierror
      integer ios
      integer j
      integer jhi
      integer jlo
      integer k
      character * ( 80 ) output_filename
      integer output_unit
      integer r(m,n)
      double precision wtime
      double precision wtime_total
      double precision x_max
      double precision x_min
      double precision x
      double precision x1
      double precision x2
      double precision y_max
      double precision y_min
      double precision y
      double precision y1
      double precision y2

      count_max = 2000
      x_max =   1.25D+00
      x_min = - 2.25D+00
      y_max =   1.75D+00
      y_min = - 1.75D+00

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MANDELBROT_OPENMP'
      write ( *, '(a)' ) '  FORTRAN77/OpenMP version'
      write ( *, '(a)' ) ' '
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
     &  '    M = ', m, ' pixels in the X direction and'
      write ( *, '(a,i8,a)' ) 
     &  '    N = ', n, ' pixels in the Y direction.'

      wtime = omp_get_wtime ( )
c
c  Carry out the iteration for each pixel, determining COUNT.
c
c$omp parallel 
c$omp&  shared ( b, count, count_max, g, r, x_max, x_min, y_max, y_min ) 
c$omp&  private ( i, j, k, x, x1, x2, y, y1, y2 )

c$omp do

      do i = 1, m

        do j = 1, n

          x = ( dble (     j - 1 ) * x_max   
     &        + dble ( m - j     ) * x_min ) 
     &        / dble ( m     - 1 )

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

          if ( mod ( count(i,j), 2 ) .eq. 1 ) then
            r(i,j) = 255
            g(i,j) = 255
            b(i,j) = 255
          else
            c = int ( 255.0D+00 * sqrt ( sqrt ( sqrt ( 
     &        ( dble ( count(i,j) ) / dble ( count_max ) ) ) ) ) )
            r(i,j) = 3 * c / 5
            g(i,j) = 3 * c / 5
            b(i,j) = c
          end if

        end do

      end do
c$omp end do

c$omp end parallel

      wtime = omp_get_wtime ( ) - wtime
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Time = ', wtime
c
c  Write data to an ASCII PPM file.
c
      output_unit = 10
      output_filename = 'mandelbrot.ppm'

      open ( unit = output_unit, file = output_filename, 
     &  status = 'replace', form = 'formatted', access = 'sequential', 
     &  iostat = ios )

      write ( output_unit, '(a2)' ) 'P3'
      write ( output_unit, '(i5,2x,i5)' ) n, m
      write ( output_unit, '(i3)' ) 255
      do i = 1, m
        do jlo = 1, n, 4
          jhi = min ( jlo + 3, n )
          write ( output_unit, '(12i5)' ) 
     &      ( r(i,j), g(i,j), b(i,j), j = jlo, jhi )
        end do
      end do

      close ( unit = output_unit )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Graphics data written to "' // trim ( output_filename ) 
     &  // '".'
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MANDELBROT_OPENMP'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

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
