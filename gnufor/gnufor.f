      subroutine get_unit ( unit )

c*********************************************************************72
c
cc GET_UNIT returns a free FORTRAN unit number.
c
c  Discussion:
c
c    A "free" FORTRAN unit number is an integer between 1 and 99 which
c    is not currently associated with an I/O device.  A free FORTRAN unit
c    number is needed in order to open a file with the OPEN command.
c
c    If UNIT = 0, then no free FORTRAN unit could be found, although
c    all 99 units were checked (except for units 5, 6 and 9, which
c    are commonly reserved for console I/O).
c
c    Otherwise, UNIT is an integer between 1 and 99, representing a
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

        if ( i .ne. 5 .and. i .ne. 6 .and. i .ne. 9 ) then

          open ( unit = i, err = 10, status = 'scratch' )
          close ( unit = i )

          unit = i

          return
        end if

10      continue

      end do

      return
      end
      subroutine run_gnuplot ( command_filename )

c*********************************************************************72
c
cc RUN_GNUPLOT runs GNUPLOT with a given command file.
c
c  Discussion:
c
c    The GNUPLOT program must be available.  To check whether
c    this is so, try typing
c
c      which gnuplot
c
c    If the response is
c
c      gnuplot: command not found
c
c    then you're going to have to make GNUPLOT available.
c
c    You may need to set the environment variable GNUTERM:
c
c      setenv GNUTERM x11
c
c    so that GNUPLOT automatically displays to your X window terminal.
c
c
c    This routine expects that there is a text file containing the appropriate
c    commands to GNUPLOT to display your picture.  There are a number of
c    routines in this package that will do this for simple plotting tasks.
c    Most of them require that you also set up a file of data to be plotted.
c
c    Once this routine invokes GNUPLOT, a graphics window should open
c    up, and the FORTRAN program will pause.  Hitting RETURN should advance
c    to the next picture, or terminate the window at the end, allowing the
c    FORTRAN routine to proceed.
c
c
c    You can look at the data and command files created by the routines.
c    Moreover, you can easily modify the command file to change the options
c    used in GNUPLOT, and then run GNUPLOT interactively, as in:
c
c      gnuplot commands
c
c    In particular, if you want a PostScript version of your graphics files,
c    insert the command "set term postscript" at the beginning of the command
c    file and run gnuplot as follows:
c
c      gnuplot commands > mypicture.ps
c
c    You will also have to hit RETURN once for each plot that is made.
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
c  Parameters:
c
c    Input, character * ( * ) COMMAND_FILENAME, the name of the
c    command file.
c
      implicit none

      character * ( 255 ) command
      character * ( * ) command_filename
      integer status
      integer system

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GNUFOR:'
      write ( *, '(a)' ) '  GNUPLOT / FORTRAN77 command interface.'
c
c  Issue a command to the system that will start GNUPLOT, using
c  the file we just wrote as input.
c
c  The "&" will run GNUPLOT in the background, so the FORTRAN program
c  can continue execution independently, and the PERSIST switch tells
c  GNUPLOT that if there are multiple plots, they should each go in
c  a separate window.
c
c  Thanks to Morag Am-Shallem for suggesting these improvements.
c  17 October 2007
c
      write ( command, * ) 
     &  'gnuplot -persist ' // trim ( command_filename ) // ' &'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GNUFOR:'
      write ( *, '(a)' ) 
     &  '  Issuing the command:"' // trim ( command ) // '".'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Press RETURN to proceed.'

      status = system ( trim ( command ) )

      if ( status .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GNUFOR - Fatal error!'
        write ( *, '(a)' ) 
     &    '  An error code was returned when the GNUPLOT command'
        write ( *, '(a)' ) 
     &    '  was issued.  Perhaps GNUPLOT is not in your path.'
        write ( *, '(a)' ) '  Type "which gnuplot" to check this.'
        stop
      end if
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GNUFOR:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

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
      subroutine write_vector_data ( data_filename, n, x, y, dx, 
     &  dy, ierror )

c*********************************************************************72
c
cc WRITE_VECTOR_DATA writes vector data to a file, for plotting by GNUPLOT.
c
c  Discussion:
c
c    Each vector is described by 4 values, X, Y, dX, dY, indicating that
c    a vector is to be drawn from (X,Y) to (X+dX,Y+dY).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 February 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) DATA_FILENAME, the name of the data file.
c
c    Input, integer N, the number of vectors.
c
c    Input, double precision X(N), Y(N), DX(N), DY(N), the vector data.
c
c    Output, integer IERROR, nonzero if an error occurred.
c
      implicit none

      integer n

      character * ( * ) data_filename
      double precision dx(n)
      double precision dy(n)
      integer file_unit
      integer i
      integer ierror
      integer ios
      double precision x(n)
      double precision y(n)

      ierror = 0

      call get_unit ( file_unit )

      if ( file_unit .eq. 0 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_VECTOR_DATA - Fatal error!'
        write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
        return
      end if

      open ( unit = file_unit, file = data_filename, 
     &  status = 'replace', iostat = ios )

      if ( ios .ne. 0 ) then
        ierror = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_VECTOR_DATA - Fatal error!'
        write ( *, '(a)' ) '  Could not open the output file.'
        return
      end if

      do i = 1, n
        write ( file_unit, * ) x(i), y(i), dx(i), dy(i)
      end do

      close ( unit = file_unit )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'WRITE_VECTOR_DATA:'
      write ( *, '(a)' ) 
     &  '  Wrote the GNUPLOT vector data file "' // 
     &  trim ( data_filename ) // '"'

      return
      end
      subroutine write_vector_plot ( command_filename, data_filename, 
     &  ierror )

c*********************************************************************72
c
cc WRITE_VECTOR_PLOT writes GNUPLOT commands to plot vectors.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) COMMAND_FILENAME, the name of the
c    command file.
c
c    Input, character * ( * ) DATA_FILENAME, the name of the data file.
c
c    Output, integer IERROR, nonzero if an error occurred.
c
      implicit none

      character * ( * ) command_filename
      character * ( * ) data_filename
      integer file_unit
      integer ierror
      integer ios
c
c  Write the data file.
c
      ierror = 0

      call get_unit ( file_unit )

      if ( file_unit .eq. 0 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_VECTOR_PLOT - Fatal error!'
        write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
        return
      end if

      open ( unit = file_unit, file = command_filename, 
     &  status = 'replace', iostat = ios )

      if ( ios .ne. 0 ) then
        ierror = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_VECTOR_PLOT - Fatal error!'
        write ( *, '(a)' ) '  Could not open the output file.'
        return
      end if

      write ( file_unit, '(a)' ) 'set title "GNUFOR plot"'
      write ( file_unit, '(a)' ) 'set xlabel "x"'
      write ( file_unit, '(a)' ) 'set ylabel "y"'
      write ( file_unit, '(a)' ) 'set style data vector'
      write ( file_unit, '(a,i2,a)' ) 'plot "' // trim ( data_filename )
      write ( file_unit, '(a)' ) 'pause -1'
      write ( file_unit, '(a)' ) 'q'

      close ( unit = file_unit )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'WRITE_VECTOR_PLOT:'
      write ( *, '(a)' ) 
     &  '  Wrote the GNUPLOT table plots command file "' // 
     &  trim ( command_filename ) // '"'

      return
      end
      subroutine write_xy_data ( data_filename, n, x, y, ierror )

c*********************************************************************72
c
cc WRITE_XY_DATA writes X(1:N), Y(1:N) data to a file.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 February 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) DATA_FILENAME, the name of the data file.
c
c    Input, integer N, the number of data items.
c
c    Input, double precision X(N), Y(N), the X and Y data
c
c    Output, integer IERROR, nonzero if an error occurred.
c
      implicit none

      integer n

      character * ( * ) data_filename
      integer file_unit
      integer i
      integer ierror
      integer ios
      double precision x(n)
      double precision y(n)

      ierror = 0

      call get_unit ( file_unit )

      if ( file_unit .eq. 0 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_XY_DATA - Fatal error!'
        write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
        return
      end if

      open ( unit = file_unit, file = data_filename, 
     &  status = 'replace', iostat = ios )

      if ( ios .ne. 0 ) then
        ierror = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_XY_DATA - Fatal error!'
        write ( *, '(a)' ) '  Could not open the output file.'
        return
      end if

      do i = 1, n
        write ( file_unit, * ) x(i), y(i)
      end do

      close ( unit = file_unit )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'WRITE_XY_DATA:'
      write ( *, '(a)' ) '  Wrote the GNUPLOT XY data file "' // 
     &  trim ( data_filename ) // '"'

      return
      end
      subroutine write_xy_plot ( command_filename, data_filename, 
     &  ierror )

c*********************************************************************72
c
cc WRITE_XY_PLOT writes GNUPLOT commands to plot X(1:N), Y(1:N) data.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) COMMAND_FILENAME, the name of the
c    command file.
c
c    Input, character * ( * ) DATA_FILENAME, the name of the data file.
c
c    Output, integer IERROR, nonzero if an error occurred.
c
      implicit none

      character * ( * ) command_filename
      character * ( * ) data_filename
      integer file_unit
      integer i
      integer ierror
      integer ios
c
c  Write the data file.
c
      ierror = 0

      call get_unit ( file_unit )

      if ( file_unit .eq. 0 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_XY_PLOT - Fatal error!'
        write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
        return
      end if

      open ( unit = file_unit, file = command_filename, 
     &  status = 'replace', iostat = ios )

      if ( ios .ne. 0 ) then
        ierror = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_XY_PLOT - Fatal error!'
        write ( *, '(a)' ) '  Could not open the output file.'
        return
      end if

      write ( file_unit, '(a)' ) 'set title "GNUFOR plot"'
      write ( file_unit, '(a)' ) 'set xlabel "x"'
      write ( file_unit, '(a)' ) 'set ylabel "y"'
      write ( file_unit, '(a,i2,a)' ) 'plot "' 
     &  // trim ( data_filename ) // '" using 1:2 with lines'
      write ( file_unit, '(a)' ) 'pause -1'
      write ( file_unit, '(a)' ) 'q'

      close ( unit = file_unit )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'WRITE_XY_PLOT:'
      write ( *, '(a)' ) 
     &  '  Wrote the GNUPLOT XY plot command file "' // 
     &  trim ( command_filename ) // '"'

      return
      end
      subroutine write_xyy_data ( data_filename, lda, nrow, ncol, x, 
     &  ierror )

c*********************************************************************72
c
cc WRITE_XYY_DATA writes a table of data to a file, for plotting by GNUPLOT.
c
c  Discussion:
c
c    The first column of data is assumed to be the independent variable, X.
c    Separate plots are made of X versus all the other columns of data.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 February 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) DATA_FILENAME, the name of the data file.
c
c    Input, integer LDA, the leading dimension of X.
c
c    Input, integer NROW, NCOL, the dimensions of X.
c
c    Input, double precision X(LDA,NCOL), the NROW by NCOL data to be written.
c
c    Output, integer IERROR, nonzero if an error occurred.
c
      implicit none

      integer lda
      integer ncol

      character * ( * ) data_filename
      integer file_unit
      integer i
      integer ierror
      integer ios
      integer nrow
      double precision x(lda,ncol)

      ierror = 0

      call get_unit ( file_unit )

      if ( file_unit .eq. 0 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_XYY_DATA - Fatal error!'
        write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
        return
      end if

      open ( unit = file_unit, file = data_filename, 
     &  status = 'replace', iostat = ios )

      if ( ios .ne. 0 ) then
        ierror = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_XYY_DATA - Fatal error!'
        write ( *, '(a)' ) '  Could not open the output file.'
        return
      end if

      do i = 1, nrow
        write ( file_unit, * ) x(i,1:ncol)
      end do

      close ( unit = file_unit )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'WRITE_XYY_DATA:'
      write ( *, '(a)' ) '  Wrote the GNUPLOT XYY data file "' // 
     &  trim ( data_filename ) // '"'

      return
      end
      subroutine write_xyy_plots ( command_filename, data_filename, 
     &  ncol, ierror )

c*********************************************************************72
c
cc WRITE_XYY_PLOTS writes GNUPLOT commands to make multiple (X,Y) plots.
c
c  Discussion:
c
c    The first column of data is assumed to be the independent variable, X.
c    Separate plots are made of X versus all the other columns of data.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) COMMAND_FILENAME, the name of the
c    command file.
c
c    Input, character * ( * ) DATA_FILENAME, the name of the data file.
c
c    Input, integer NCOL, the number of columns of data.
c
c    Output, integer IERROR, nonzero if an error occurred.
c
      implicit none

      character * ( * ) command_filename
      character * ( * ) data_filename
      integer file_unit
      integer i
      integer ierror
      integer ios
      integer ncol
c
c  Write the data file.
c
      ierror = 0

      call get_unit ( file_unit )

      if ( file_unit .eq. 0 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_XYY_PLOTS - Fatal error!'
        write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
        return
      end if

      open ( unit = file_unit, file = command_filename, 
     &  status = 'replace', iostat = ios )

      if ( ios .ne. 0 ) then
        ierror = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_XYY_PLOTS - Fatal error!'
        write ( *, '(a)' ) '  Could not open the output file.'
        return
      end if

      write ( file_unit, '(a)' ) 'set title "GNUFOR plot"'
      write ( file_unit, '(a)' ) 'set xlabel "x"'
      write ( file_unit, '(a)' ) 'set ylabel "y"'
      do i = 2, ncol
        write ( file_unit, '(a,i2,a)' ) 'plot "' // 
     &    trim ( data_filename ) // '" using ', i, ' with lines'
        write ( file_unit, '(a)' ) 'pause -1'
      end do
      write ( file_unit, '(a)' ) 'q'

      close ( unit = file_unit )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'WRITE_XYY_PLOTS:'
      write ( *, '(a)' ) 
     &  '  Wrote the GNUPLOT XYY plots command file "' // 
     &  trim ( command_filename ) // '"'

      return
      end
      subroutine write_xyz_data ( data_filename, n, x, y, z, ierror )

c*********************************************************************72
c
cc WRITE_XYZ_DATA writes X(1:N), Y(1:N), Z(1:N) data to a file.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 February 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) DATA_FILENAME, the name of the data file.
c
c    Input, integer N, the number of data items.
c
c    Input, double precision X(N), Y(N), Z(N), the X, Y, Z data
c
c    Output, integer IERROR, nonzero if an error occurred.
c
      implicit none

      integer n

      character * ( * ) data_filename
      integer file_unit
      integer i
      integer ierror
      integer ios
      double precision x(n)
      double precision y(n)
      double precision z(n)

      ierror = 0

      call get_unit ( file_unit )

      if ( file_unit .eq. 0 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_XYZ_DATA - Fatal error!'
        write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
        return
      end if

      open ( unit = file_unit, file = data_filename, 
     &  status = 'replace', iostat = ios )

      if ( ios .ne. 0 ) then
        ierror = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_XYZ_DATA - Fatal error!'
        write ( *, '(a)' ) '  Could not open the output file.'
        return
      end if

      do i = 1, n
        write ( file_unit, * ) x(i), y(i), z(i)
      end do

      close ( unit = file_unit )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'WRITE_XYZ_DATA:'
      write ( *, '(a)' ) '  Wrote the GNUPLOT XYZ data file "' // 
     &  trim ( data_filename ) // '"'

      return
      end
      subroutine write_xyz_plot ( command_filename, data_filename, 
     &  ierror )

c*********************************************************************72
c
cc WRITE_XYZ_PLOT writes commands to plot parametric (X,Y,Z) data.
c
c  Discussion:
c
c    This routine tries to write a command file suitable for displaying
c    a 3D arc specified by points (X,Y,Z).  A grid data file, containing
c    values of X, Y and Z, will also be needed.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 February 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) COMMAND_FILENAME, the name of the
c    command file.
c
c    Input, character * ( * ) DATA_FILENAME, the name of the data file.
c
c    Output, integer IERROR, nonzero if an error occurred.
c
      implicit none

      character * ( * ) command_filename
      character * ( * ) data_filename
      integer file_unit
      integer i
      integer ierror
      integer ios
      integer ncol
c
c  Write the data file.
c
      ierror = 0

      call get_unit ( file_unit )

      if ( file_unit .eq. 0 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_XYZ_PLOT - Fatal error!'
        write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
        return
      end if

      open ( unit = file_unit, file = command_filename, 
     &  status = 'replace', iostat = ios )

      if ( ios .ne. 0 ) then
        ierror = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_XYZ_PLOT - Fatal error!'
        write ( *, '(a)' ) '  Could not open the output file.'
        return
      end if

      write ( file_unit, '(a)' ) 'set title "GNUFOR plot"'
      write ( file_unit, '(a)' ) 'set xlabel "x"'
      write ( file_unit, '(a)' ) 'set ylabel "y"'
      write ( file_unit, '(a)' ) 'set parametric'
      write ( file_unit, '(a)' ) 
     &  'splot "' // trim ( data_filename ) // 
     &    '" using 1:2:3 with lines'
      write ( file_unit, '(a)' ) 'pause -1'
      write ( file_unit, '(a)' ) 'q'

      close ( unit = file_unit )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'WRITE_XYZ_PLOT:'
      write ( *, '(a)' ) 
     &  '  Wrote the GNUPLOT SPLOT command file "' // 
     &  trim ( command_filename ) // '"'

      return
      end
      subroutine write_xyzgrid_contour ( command_filename, 
     &  data_filename, ierror )

c*********************************************************************72
c
cc WRITE_XYZGRID_CONTOUR writes commands to plot contours of Z(X,Y).
c
c  Discussion:
c
c    This routine tries to write a command file suitable for displaying
c    contours of Z(X,Y) gridded data.  A grid data file, containing values
c    of X, Y and Z, will also be needed.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) COMMAND_FILENAME, the name of the
c    command file.
c
c    Input, character * ( * ) DATA_FILENAME, the name of the data file.
c
c    Output, integer IERROR, nonzero if an error occurred.
c
      implicit none

      character * ( * ) command_filename
      character * ( * ) data_filename
      integer file_unit
      integer i
      integer ierror
      integer ios
      integer ncol
c
c  Write the data file.
c
      ierror = 0

      call get_unit ( file_unit )

      if ( file_unit .eq. 0 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_XYZGRID_CONTOUR - Fatal error!'
        write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
        return
      end if

      open ( unit = file_unit, file = command_filename, 
     &  status = 'replace', iostat = ios )

      if ( ios .ne. 0 ) then
        ierror = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_XYZGRID_CONTOUR - Fatal error!'
        write ( *, '(a)' ) '  Could not open the output file.'
        return
      end if

      write ( file_unit, '(a)' ) 'set title "GNUFOR plot"'
      write ( file_unit, '(a)' ) 'set xlabel "x"'
      write ( file_unit, '(a)' ) 'set ylabel "y"'
      write ( file_unit, '(a)' ) 'set parametric'
      write ( file_unit, '(a)' ) 'set nosurface'
      write ( file_unit, '(a)' ) 'set contour'
      write ( file_unit, '(a)' ) 'set cntrparam levels 10'
      write ( file_unit, '(a)' ) 
     &  'splot "' // trim ( data_filename ) // 
     &  '" using 1:2:3 with lines'

      write ( file_unit, '(a)' ) 'set table "table.txt"'
      write ( file_unit, '(a)' ) 'replot'
      write ( file_unit, '(a)' ) 'unset table'

      write ( file_unit, '(a)' ) 'set term x11'
      write ( file_unit, '(a)' ) 'plot "table.txt" with lines'
      write ( file_unit, '(a)' ) 'pause -1'
      write ( file_unit, '(a)' ) 'q'

      close ( unit = file_unit )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'WRITE_XYZGRID_CONTOUR:'
      write ( *, '(a)' ) 
     &  '  Wrote the GNUPLOT XYZGRID contour plot command file "' // 
     &  trim ( command_filename ) // '"'

      return
      end
      subroutine write_xyzgrid_data ( data_filename, nx, ny, xyz, 
     &  ierror )

c*********************************************************************72
c
cc WRITE_XYZGRID_DATA writes a file of XYZ grid data.
c
c  Discussion:
c
c    It is assumed that values of Z are available on a regular NX by NY grid
c    of (X,Y) points.
c
c    The form of the data file requires that all the data for a given value
c    of Y be listed, followed by a blank line, followed by the data for
c    another value of Y.
c
c  Example:
c
c    Here is a grid data file for a 3 by 3 grid, with Z = X + Y.
c
c    0.0 0.0 0.0
c    1.0 0.0 1.0
c    2.0 0.0 2.0
c
c    0.0 1.0 1.0
c    1.0 1.0 2.0
c    2.0 1.0 3.0
c
c    0.0 2.0 2.0
c    1.0 2.0 3.0
c    2.0 2.0 4.0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 February 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) DATA_FILENAME, the name of the data file.
c
c    Input, integer NX, NY, the dimensions of the grid.
c
c    Input, double precision XYZ(3,NX,NY), the XYZ grid data to be written.
c
c    Output, integer IERROR, nonzero if an error occurred.
c
      implicit none

      integer nx
      integer ny

      character * ( * ) data_filename
      integer file_unit
      integer i
      integer ierror
      integer ios
      integer j
      double precision xyz(3,nx,ny)

      ierror = 0

      call get_unit ( file_unit )

      if ( file_unit .eq. 0 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_XYZGRID_DATA - Fatal error!'
        write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
        return
      end if

      open ( unit = file_unit, file = data_filename, 
     &  status = 'replace', iostat = ios )

      if ( ios .ne. 0 ) then
        ierror = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_XYZGRID_DATA - Fatal error!'
        write ( *, '(a)' ) '  Could not open the output file.'
        return
      end if

      do j = 1, ny
        do i = 1, nx
          write ( file_unit, * ) xyz(1:3,i,j)
        end do
        write ( file_unit, '(a)' )
      end do

      close ( unit = file_unit )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'WRITE_XYZGRID_DATA:'
      write ( *, '(a)' ) '  Wrote the GNUPLOT XYZ grid data file "' // 
     &  trim ( data_filename ) // '"'

      return
      end
      subroutine write_xyzgrid_surface ( command_filename, 
     &  data_filename, ierror )

c*********************************************************************72
c
cc WRITE_XYZGRID_SURFACE writes a file of GNUPLOT commands to plot a 3D surface.
c
c  Discussion:
c
c    This routine tries to write a command file suitable for displaying
c    a surface Z(X,Y).  A grid data file, containing values of X, Y and Z,
c    will also be needed.  The routine WRITE_XYZGRID_DATA can write this file.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 February 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) COMMAND_FILENAME, the name of the
c    command file.
c
c    Input, character * ( * ) DATA_FILENAME, the name of the data file.
c
c    Output, integer IERROR, nonzero if an error occurred.
c
      implicit none

      character * ( * ) command_filename
      character * ( * ) data_filename
      integer file_unit
      integer i
      integer ierror
      integer ios
      integer ncol
c
c  Write the data file.
c
      ierror = 0

      call get_unit ( file_unit )

      if ( file_unit .eq. 0 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_XYZGRID_SURFACE - Fatal error!'
        write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
        return
      end if

      open ( unit = file_unit, file = command_filename, 
     &  status = 'replace', iostat = ios )

      if ( ios .ne. 0 ) then
        ierror = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'WRITE_XYZGRID_SURFACE - Fatal error!'
        write ( *, '(a)' ) '  Could not open the output file.'
        return
      end if

      write ( file_unit, '(a)' ) 'set title "GNUFOR plot"'
      write ( file_unit, '(a)' ) 'set xlabel "x"'
      write ( file_unit, '(a)' ) 'set ylabel "y"'
      write ( file_unit, '(a)' ) 'set parametric'
      write ( file_unit, '(a)' ) 'set hidden3d'
      write ( file_unit, '(a)' ) 'set contour'
      write ( file_unit, '(a)' ) 
     &  'splot "' // trim ( data_filename ) // 
     &  '" using 1:2:3 with lines'
      write ( file_unit, '(a)' ) 'pause -1'
      write ( file_unit, '(a)' ) 'q'

      close ( unit = file_unit )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'WRITE_SURFACE_COMMANDS:'
      write ( *, '(a)' ) 
     &  '  Wrote the GNUPLOT surface plot command file "' // 
     &  trim ( command_filename ) // '"'

      return
      end

