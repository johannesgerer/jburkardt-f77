      program main

c*********************************************************************72
c
cc MAIN is the main program for GNUFOR_PRB.
c
c  Discussion:
c
c    GNUFOR_PRB tests the GNUFOR routines.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GNUFOR_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the GNUFOR library.'

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
      write ( *, '(a)' ) 'GNUFOR_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 demonstrates the plotting of Y(X) data.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 101 )

      double precision angle
      double precision area
      character * ( 255 ) command_filename
      character * ( 255 ) data_filename
      integer i
      integer ierror
      double precision x(n)
      double precision xmin
      double precision xmax
      double precision y(n)

      command_filename = 'test01_commands.txt'
      data_filename = 'test01_data.txt'
      xmin = 0.0D+00
      xmax = 20.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  To plot a simple set of (X,Y) data,'
      write ( *, '(a)' ) '  WRITE_XY_DATA writes the data file,'
      write ( *, '(a)' ) '  WRITE_XY_PLOT writes the plot command file.'

      do i = 1, n

        x(i) = ( dble ( n - i     ) * xmin   
     &        + dble (     i - 1 ) * xmax ) 
     &          / dble ( n     - 1 )

        y(i) = sin ( x(i) ) * sin ( 4.0D+00 * x(i) )

      end do

      call write_xy_data ( data_filename, n, x, y, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST01'
        write ( *, '(a,i6)' ) 
     &    '  WRITE_XY_DATA returned IERROR = ', ierror
      end if

      call write_xy_plot ( command_filename, data_filename, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST01'
        write ( *, '(a,i6)' ) 
     &    '  WRITE_XY_PLOT returned IERROR = ', ierror
      end if

      call run_gnuplot ( command_filename )

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 demonstrates the plotting of a table of data.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nrow
      parameter ( nrow = 101 )
      integer ncol
      parameter ( ncol = 4 )

      integer lda
      parameter ( lda = nrow )

      double precision angle
      double precision area
      character * ( 255 ) command_filename
      character * ( 255 ) data_filename
      double precision height
      integer i
      integer ierror
      double precision r
      double precision width
      double precision x(lda,ncol)

      command_filename = 'test02_commands.txt'
      data_filename = 'test02_data.txt'
      r = 50.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  To plot X versus multiple sets of Y data,'
      write ( *, '(a)' ) '  WRITE_XYY_DATA writes the data file,'
      write ( *, '(a)' ) 
     &  '  WRITE_XYY_PLOT writes the plot command file.'

      do i = 1, nrow

        height = 2.0d+00 * r * dble ( i - 1 ) / dble ( nrow - 1 )
        width = 2.0D+00 * sqrt ( r**2 - ( r - height )**2 )
        angle = acos ( ( r - height ) / r )
        area = 0.5D+00 * r**2 * 2.0d+00 * acos ( ( r - height ) / r ) 
     &    - ( r - height ) * sqrt ( height * ( 2.0D+00 * r - height ) )

        x(i,1) = height
        x(i,2) = width
        x(i,3) = angle
        x(i,4) = area

      end do

      call write_xyy_data ( data_filename, lda, nrow, ncol, x, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST02'
        write ( *, '(a,i6)' ) 
     &    '  WRITE_XYY_DATA returned IERROR = ', ierror
      end if

      call write_xyy_plots ( command_filename, data_filename, ncol, 
     &  ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST02'
        write ( *, '(a,i6)' ) 
     &    '  WRITE_XYY_PLOTS returned IERROR = ', ierror
      end if

      call run_gnuplot ( command_filename )

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 plots parameter (X,Y,Z) data.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 101 )

      character * ( 255 ) command_filename
      character * ( 255 ) data_filename
      integer i
      integer ierror
      integer nturn 
      double precision pi
      double precision r
      double precision theta
      double precision x(n)
      double precision y(n)
      double precision z(n)
      double precision zmax

      command_filename = 'test03_commands.txt'
      data_filename = 'test03_data.txt'
      nturn = 5
      pi = 3.141592653589793D+00
      r = 5.0D+00
      zmax = 10.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  To plot a (parametric) set of (X,Y,Z) data,'
      write ( *, '(a)' ) '  WRITE_XYZ_DATA writes the data file,'
      write ( *, '(a)' ) 
     &  '  WRITE_XYZ_PLOT writes the plot command file.'

      do i = 1, n

        z(i) = zmax * dble ( i - 1 ) / dble ( n - 1 )

        theta = ( 2.0D+00 * pi ) * z(i) * dble ( nturn ) / zmax

        x(i) = r * cos ( theta )
        y(i) = r * sin ( theta )

      end do

      call write_xyz_data ( data_filename, n, x, y, z, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST03'
        write ( *, '(a,i6)' ) 
     &    '  WRITE_XYZ_DATA returned IERROR = ', ierror
      end if

      call write_xyz_plot ( command_filename, data_filename, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST03'
        write ( *, '(a,i6)' ) 
     &    '  WRITE_XYZ_PLOT returned IERROR = ', ierror
      end if

      call run_gnuplot ( command_filename )

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 plots vector data.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nx
      parameter ( nx = 21 )
      integer ny
      parameter ( ny = 21 )
      integer n
      parameter ( n = nx * ny )

      character * ( 255 ) command_filename
      character * ( 255 ) data_filename
      double precision dx(n)
      double precision dy(n)
      integer i
      integer ierror
      integer j
      integer k
      double precision x(n)
      double precision xmax
      double precision xmin
      double precision xx
      double precision y(n)
      double precision ymax
      double precision ymin
      double precision yy

      command_filename = 'test04_commands.txt'
      data_filename = 'test04_data.txt'
      xmax = 1.0D+00
      xmin = -1.0D+00
      ymax = 1.0D+00
      ymin = -1.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) '  To plot a vector field,'
      write ( *, '(a)' ) '  WRITE_VECTOR_DATA writes the data file,'
      write ( *, '(a)' ) 
     &  '  WRITE_VECTOR_PLOT writes the plot command file.'

      k = 0

      do i = 1, nx

        do j = 1, ny

          k = k + 1

          xx = ( dble ( nx - i     ) * xmin   
     &         + dble (      i - 1 ) * xmax ) 
     &         / dble ( nx     - 1 )

          yy = ( dble ( ny - j     ) * ymin   
     &         + dble (      j - 1 ) * ymax ) 
     &         / dble ( ny     - 1 )

          dx(k) = - 0.10D+00 * yy
          dy(k) =   0.10D+00 * xx

          x(k) = xx - 0.5D+00 * dx(k)
          y(k) = yy - 0.5D+00 * dy(k)

        end do

      end do

      call write_vector_data ( data_filename, n, x, y, dx, dy, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST04'
        write ( *, '(a,i6)' ) 
     &    '  WRITE_VECTOR_DATA returned IERROR = ', ierror
      end if

      call write_vector_plot ( command_filename, data_filename, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST04'
        write ( *, '(a,i6)' ) 
     &    '  WRITE_VECTOR_PLOT returned IERROR = ', ierror
      end if

      call run_gnuplot ( command_filename )

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 plots Z(X,Y) grid data as a surface.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nx
      parameter ( nx = 21 )
      integer ny
      parameter ( ny = 21 )

      integer nrow
      parameter ( nrow = nx * ny )

      character * ( 255 ) command_filename
      character * ( 255 ) data_filename
      integer i
      integer ierror
      integer j
      double precision x
      double precision xmax
      double precision xmin
      double precision xyz(3,nx,ny)
      double precision y
      double precision ymax
      double precision ymin
      double precision z

      command_filename = 'test05_commands.txt'
      data_filename = 'test05_data.txt'
      xmax = 1.0D+00
      xmin = 0.0D+00
      ymax = 1.0D+00
      ymin = 0.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) 
     &  '  To plot a gridded set of Z(X,Y) data as a surface,'
      write ( *, '(a)' ) '  WRITE_XYZGRID_DATA writes the data file,'
      write ( *, '(a)' ) 
     &  '  WRITE_XYZGRID_SURFACE writes the plot command file.'

      do i = 1, nx

        x = ( dble ( nx - i     ) * xmin   
     &      + dble (      i - 1 ) * xmax ) 
     &      / dble ( nx     - 1 )

        do j = 1, ny

          y = ( dble ( ny - j     ) * ymin   
     &        + dble (      j - 1 ) * ymax ) 
     &        / dble ( ny     - 1 )

          z = sin ( 64.0D+00 * ( x - 0.5D+00 )**2 * ( y - 0.5D+00 )**2 )

          xyz(1,i,j) = x
          xyz(2,i,j) = y
          xyz(3,i,j) = z

        end do

      end do

      call write_xyzgrid_data ( data_filename, nx, ny, xyz, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST05'
        write ( *, '(a,i6)' ) 
     &    '  WRITE_XYZGRID_DATA returned IERROR = ', ierror
      end if

      call write_xyzgrid_surface ( command_filename, data_filename, 
     &  ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST05'
        write ( *, '(a,i6)' ) 
     &    '  WRITE_XYZGRID_SURFACE returned IERROR = ', ierror
      end if

      call run_gnuplot ( command_filename )

      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 plots Z(X,Y) grid data as contours.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nx
      parameter ( nx = 41 )
      integer ny
      parameter ( ny = 41 )

      character * ( 255 ) command_filename
      character * ( 255 ) data_filename
      integer i
      integer ierror
      integer j
      double precision x
      double precision xmax
      double precision xmin
      double precision xyz(3,nx,ny)
      double precision y
      double precision ymax
      double precision ymin
      double precision z

      command_filename = 'test06_commands.txt'
      data_filename = 'test06_data.txt'
      xmax = 1.0D+00
      xmin = 0.0D+00
      ymax = 1.0D+00
      ymin = 0.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' ) '  To plot gridded Z(X,Y) data as contours,'
      write ( *, '(a)' ) '  WRITE_XYZGRID_DATA writes the data file,'
      write ( *, '(a)' ) 
     &  '  WRITE_XYZGRID_CONTOUR writes the plot command file.'

      do i = 1, nx

        x = ( dble ( nx - i     ) * xmin   
     &      + dble (      i - 1 ) * xmax ) 
     &      / dble ( nx     - 1 )

        do j = 1, ny

          y = ( dble ( ny - j     ) * ymin   
     &        + dble (      j - 1 ) * ymax ) 
     &        / dble ( ny     - 1 )

          z = sin ( 64.0D+00 * ( x - 0.5D+00 )**2 * ( y - 0.5D+00 )**2 )

          xyz(1,i,j) = x
          xyz(2,i,j) = y
          xyz(3,i,j) = z

        end do

      end do

      call write_xyzgrid_data ( data_filename, nx, ny, xyz, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST06'
        write ( *, '(a,i6)' ) 
     &    '  WRITE_XYZGRID_DATA returned IERROR = ', ierror
      end if

      call write_xyzgrid_contour ( command_filename, data_filename, 
     &  ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST06'
        write ( *, '(a,i6)' ) 
     &    '  WRITE_XYZGRID_CONTOUR returned IERROR = ', ierror
      end if

      call run_gnuplot ( command_filename )

      return
      end
