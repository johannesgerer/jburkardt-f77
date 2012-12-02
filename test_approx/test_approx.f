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
      subroutine p00_dat ( prob, data_num, x, y )

c*********************************************************************72
c
cc P00_DAT returns the data vector for any problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the number of the desired test problem.
c
c    Input, integer DATA_NUM, the number of data points,
c    as specified by P00_DATA_NUM.
c
c    Output, double precision X(DATA_NUM), the abscissa data.
c
c    Output, double precision Y(DATA_NUM), the ordinate data.
c
      implicit none

      integer data_num

      integer prob
      double precision x(data_num)
      double precision y(data_num)

      if ( prob .eq. 1 ) then
        call p01_dat ( data_num, x, y )
      else if ( prob .eq. 2 ) then
        call p02_dat ( data_num, x, y )
      else if ( prob .eq. 3 ) then
        call p03_dat ( data_num, x, y )
      else if ( prob .eq. 4 ) then
        call p04_dat ( data_num, x, y )
      else if ( prob .eq. 5 ) then
        call p05_dat ( data_num, x, y )
      else if ( prob .eq. 6 ) then
        call p06_dat ( data_num, x, y )
      else if ( prob .eq. 7 ) then
        call p07_dat ( data_num, x, y )
      else if ( prob .eq. 8 ) then
        call p08_dat ( data_num, x, y )
      else if ( prob .eq. 9 ) then
        call p09_dat ( data_num, x, y )
      else if ( prob .eq. 10 ) then
        call p10_dat ( data_num, x, y )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_DAT - Fatal error!'
        write ( *, '(a,i6)' ) '  Illegal problem number = ', prob
        stop
      end if

      return
      end
      subroutine p00_data_num ( prob, data_num )

c*********************************************************************72
c
cc P00_DATA_NUM returns the dimension of the data vector for any problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DATA_NUM, the dimension of the data vector.
c
c    Output, integer DATA_NUM, the dimension of the data vector.
c
      implicit none

      integer prob
      integer data_num

      if ( prob .eq. 1 ) then
        call p01_data_num ( data_num )
      else if ( prob .eq. 2 ) then
        call p02_data_num ( data_num )
      else if ( prob .eq. 3 ) then
        call p03_data_num ( data_num )
      else if ( prob .eq. 4 ) then
        call p04_data_num ( data_num )
      else if ( prob .eq. 5 ) then
        call p05_data_num ( data_num )
      else if ( prob .eq. 6 ) then
        call p06_data_num ( data_num )
      else if ( prob .eq. 7 ) then
        call p07_data_num ( data_num )
      else if ( prob .eq. 8 ) then
        call p08_data_num ( data_num )
      else if ( prob .eq. 9 ) then
        call p09_data_num ( data_num )
      else if ( prob .eq. 10 ) then
        call p10_data_num ( data_num )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_DATA_NUM - Fatal error!'
        write ( *, '(a,i6)' ) '  Illegal problem number = ', prob
        stop
      end if

      return
      end
      subroutine p00_prob_num ( prob_num )

c*********************************************************************72
c
cc P00_PROB_NUM returns the number of test problems.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer PROB_NUM, the number of test problems.
c
      implicit none

      integer prob_num

      prob_num = 10

      return
      end
      subroutine p00_story ( prob )

c*********************************************************************72
c
cc P00_STORY prints the "story" for any problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
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

      integer prob

      if ( prob .eq. 1 ) then
        call p01_story ( )
      else if ( prob .eq. 2 ) then
        call p02_story ( )
      else if ( prob .eq. 3 ) then
        call p03_story ( )
      else if ( prob .eq. 4 ) then
        call p04_story ( )
      else if ( prob .eq. 5 ) then
        call p05_story ( )
      else if ( prob .eq. 6 ) then
        call p06_story ( )
      else if ( prob .eq. 7 ) then
        call p07_story ( )
      else if ( prob .eq. 8 ) then
        call p08_story ( )
      else if ( prob .eq. 9 ) then
        call p09_story ( )
      else if ( prob .eq. 10 ) then
        call p10_story ( )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_STORY - Fatal error!'
        write ( *, '(a)' ) '  Unexpected input value of PROB.'
        stop
      end if

      return
      end
      subroutine p00_title ( prob, title )

c*********************************************************************72
c
cc P00_TITLE returns the title of any problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the number of the desired test problem.
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      integer prob
      character ( len = * ) title

      if ( prob .eq. 1 ) then
        call p01_title ( title )
      else if ( prob .eq. 2 ) then
        call p02_title ( title )
      else if ( prob .eq. 3 ) then
        call p03_title ( title )
      else if ( prob .eq. 4 ) then
        call p04_title ( title )
      else if ( prob .eq. 5 ) then
        call p05_title ( title )
      else if ( prob .eq. 6 ) then
        call p06_title ( title )
      else if ( prob .eq. 7 ) then
        call p07_title ( title )
      else if ( prob .eq. 8 ) then
        call p08_title ( title )
      else if ( prob .eq. 9 ) then
        call p09_title ( title )
      else if ( prob .eq. 10 ) then
        call p10_title ( title )
      else
        title = ' '
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_TITLE - Fatal error!'
        write ( *, '(a,i6)' ) '  Illegal problem number = ', prob
        stop
      end if

      return
      end
      subroutine p01_dat ( data_num, x, y )

c*********************************************************************72
c
cc P01_DAT returns the data vector for problem 1.
c
c  Discussion:
c
c    The X data is measured in days, and the Y data represents the
c    observed position of Mars in a heliocentric coordinate system.
c
c    The X data is equally spaced.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Samuel Conte, Carl deBoor,
c    Elementary Numerical Analysis,
c    McGraw Hill, 1972, page 217.
c
c  Parameters:
c
c    Input, integer DATA_NUM, the dimension of the data vector.
c
c    Output, double precision X(DATA_NUM), the abscissa data.
c
c    Output, double precision Y(DATA_NUM), the ordinate data.
c
      implicit none

      integer data_num

      double precision x(data_num)
      double precision x_save(10)
      double precision y(data_num)
      double precision y_save(10)

      data x_save /
     &  1250.5D+00,
     &  1260.5D+00,
     &  1270.5D+00,
     &  1280.5D+00,
     &  1290.5D+00, 
     &  1300.5D+00,
     &  1310.5D+00,
     &  1320.5D+00,
     &  1330.5D+00,
     &  1340.5D+00 /
      data y_save /
     &  1.39140D+00,
     &  1.37696D+00,
     &  1.34783D+00,
     &  1.30456D+00,
     &  1.24787D+00, 
     &  1.17862D+00,
     &  1.09776D+00,
     &  1.00636D+00,
     &  0.90553D+00,
     &  0.79642D+00 /

      call r8vec_copy ( data_num, x_save, x )
      call r8vec_copy ( data_num, y_save, y )

      return
      end
      subroutine p01_data_num ( data_num )

c*********************************************************************72
c
cc P01_DATA_NUM returns the dimension of the data vector for problem 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer DATA_NUM, the dimension of the data vector.
c
      implicit none

      integer data_num

      data_num = 10

      return
      end
      subroutine p01_story ( )

c*********************************************************************72
c
cc P01_STORY prints the "story" for problem 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
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

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  This example is due to deBoor.'
      write ( *, '(a)' ) 
     &  '  For this example, X is measured in days, and'
      write ( *, '(a)' ) 
     &  '  Y records the observed position of Mars in a heliocentric'
      write ( *, '(a)' ) 
     &  '  coordinate system.'

      return
      end
      subroutine p01_title ( title )

c*********************************************************************72
c
cc P01_TITLE returns the title of problem 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'DeBoor example, Mars position'

      return
      end
      subroutine p02_dat ( data_num, x, y )

c*********************************************************************72
c
cc P02_DAT returns the data vector for problem 2.
c
c  Discussion:
c
c    The data lies roughly along a straight line.  Polynomial
c    interpolation is inappropriate.  Instead, a least squares
c    approximation should be sought, of the form:
c
c      F(X) = A + B * X
c
c    The X data is equally spaced.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Samuel Conte, Carl deBoor,
c    Elementary Numerical Analysis,
c    McGraw Hill, 1972, page 217.
c
c  Parameters:
c
c    Input, integer DATA_NUM, the dimension of the data vector.
c
c    Output, double precision X(DATA_NUM), the abscissa data.
c
c    Output, double precision Y(DATA_NUM), the ordinate data.
c
      implicit none

      integer data_num

      double precision x(data_num)
      double precision x_save(11)
      double precision y(data_num)
      double precision y_save(11)

      data x_save /
     &  1.0D+00, 2.0D+00, 3.0D+00, 4.0D+00,  5.0D+00, 
     &  6.0D+00, 7.0D+00, 8.0D+00, 9.0D+00, 10.0D+00, 
     & 11.0D+00 /
      data y_save /
     &  0.00D+00, 0.60D+00, 1.77D+00, 1.92D+00, 3.31D+00, 
     &  3.52D+00, 4.59D+00, 5.31D+00, 5.79D+00, 7.06D+00, 
     &  7.17D+00 /

      call r8vec_copy ( data_num, x_save, x )
      call r8vec_copy ( data_num, y_save, y )

      return
      end
      subroutine p02_data_num ( data_num )

c*********************************************************************72
c
cc P02_DATA_NUM returns the dimension of the data vector for problem 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer DATA_NUM, the dimension of the data vector.
c
      implicit none

      integer data_num

      data_num = 11

      return
      end
      subroutine p02_story ( )

c*********************************************************************72
c
cc P02_STORY prints the "story" for problem 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
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

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  This example is due to deBoor.'
      write ( *, '(a)' ) 
     &  '  The data lies roughly along a straight line.  Polynomial'
      write ( *, '(a)' ) 
     &  '  interpolation is inappropriate.  Instead, a least squares'
      write ( *, '(a)' ) 
     &  '  approximation should be sought, of the form:'
      write ( *, '(a)' ) 
     &  '    F(X) = A + B * X'

      return
      end
      subroutine p02_title ( title )

c*********************************************************************72
c
cc P02_TITLE returns the title of problem 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'DeBoor example, roughly linear data'

      return
      end
      subroutine p03_dat ( data_num, x, y )

c*********************************************************************72
c
cc P03_DAT returns the data vector for problem 3.
c
c  Discussion:
c
c    The data is all zero except for a single value of 1 in the center.
c    This data set is interesting because an interpolation method that
c    is "local" will produce an interpolating curve that is exactly
c    zero over most of the outlying intervals, whereas a nonlocal
c    interpolation method may produce a curve that "wiggles" over the
c    entire interpolation interval.
c
c    The X data is equally spaced.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DATA_NUM, the dimension of the data vector.
c
c    Output, double precision X(DATA_NUM), the abscissa data.
c
c    Output, double precision Y(DATA_NUM), the ordinate data.
c
      implicit none

      integer data_num

      double precision x(data_num)
      double precision x_save(11)
      double precision y(data_num)
      double precision y_save(11)

      data x_save /
     &   0.0D+00, 1.0D+00, 2.0D+00, 3.0D+00, 4.0D+00, 
     &   5.0D+00, 6.0D+00, 7.0D+00, 8.0D+00, 9.0D+00, 
     &  10.0D+00 /
      data y_save /
     &  0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 
     &  1.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 
     &  0.0D+00 /

      call r8vec_copy ( data_num, x_save, x )
      call r8vec_copy ( data_num, y_save, y )

      return
      end
      subroutine p03_data_num ( data_num )

c*********************************************************************72
c
cc P03_DATA_NUM returns the dimension of the data vector for problem 3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer DATA_NUM, the dimension of the data vector.
c
      implicit none

      integer data_num

      data_num = 11

      return
      end
      subroutine p03_story ( )

c*********************************************************************72
c
cc P03_STORY prints the "story" for problem 3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
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

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  The data is all zero except for a single value of 1 in '
      write ( *, '(a)' ) 
     &  '  the center.  This data set is interesting because an '
      write ( *, '(a)' ) 
     &  '  interpolation method that is "local" will produce an '
      write ( *, '(a)' ) 
     &  '  interpolating curve that is exactly zero over most of '
      write ( *, '(a)' ) 
     &  '  the outlying intervals, whereas a nonlocal interpolation '
      write ( *, '(a)' ) 
     &  '  method may produce a curve that "wiggles" over the'
      write ( *, '(a)' ) 
     &  '  entire interpolation interval.'

      return
      end
      subroutine p03_title ( title )

c*********************************************************************72
c
cc P03_TITLE returns the title of problem 3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'The pulse data, 0 0 0 0 0 1 0 0 0 0 0'

      return
      end
      subroutine p04_dat ( data_num, x, y )

c*********************************************************************72
c
cc P04_DAT returns the data vector for problem 4.
c
c  Discussion:
c
c    Theoretically, the data is a step, 0 to the left of 5, and 1
c    to the right.  To keep things simple, the data is defined
c    to be 0 up to 5 - RADIUS, 1/2 at 5, 1 at 5 + RADIUS and beyond,
c    with RADIUS set to a "small" value, currently 0.01.
c    Some interpolation methods will violently overreact to this
c    jump.
c
c    The X data is NOT equally spaced because of the handling of the pulse.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DATA_NUM, the dimension of the data vector.
c
c    Output, double precision X(DATA_NUM), the abscissa data.
c
c    Output, double precision Y(DATA_NUM), the ordinate data.
c
      implicit none

      integer data_num

      double precision x(data_num)
      double precision x_save(13)
      double precision y(data_num)
      double precision y_save(13)

      data x_save /
     &  0.0D+00,  1.0D+00, 2.0D+00,  3.0D+00, 4.0D+00, 
     &  4.99D+00, 5.0D+00, 5.01D+00, 6.0D+00, 7.0D+00, 
     &  8.0D+00,  9.0D+00, 10.0D+00 /

      data y_save / 
     &  0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 
     &  0.0D+00, 0.5D+00, 1.0D+00, 1.0D+00, 1.0D+00, 
     &  1.0D+00, 1.0D+00, 1.0D+00 /

      call r8vec_copy ( data_num, x_save, x )
      call r8vec_copy ( data_num, y_save, y )

      return
      end
      subroutine p04_data_num ( data_num )

c*********************************************************************72
c
cc P04_DATA_NUM returns the dimension of the data vector for problem 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer DATA_NUM, the dimension of the data vector.
c
      implicit none

      integer data_num

      data_num = 13

      return
      end
      subroutine p04_story ( )

c*********************************************************************72
c
cc P04_STORY prints the "story" for problem 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
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

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Theoretically, the data is a step, 0 to the left of 5, '
      write ( *, '(a)' ) 
     &  '  and 1 to the right.  To keep things simple, the data is '
      write ( *, '(a)' ) 
     &  '  defined to be 0 up to 5 - RADIUS, 1/2 at 5, 1 at 5 + RADIUS'
      write ( *, '(a)' ) 
     &  '   and beyond, with RADIUS set to a "small" value, '
      write ( *, '(a)' ) 
     &  '  currently 0.01.  Some interpolation methods will violently '
      write ( *, '(a)' ) 
     &  '  overreact to this jump.'

      return
      end
      subroutine p04_title ( title )

c*********************************************************************72
c
cc P04_TITLE returns the title of problem 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'The jump data, 0 0 0 0 0 1/2 1 1 1 1 1'

      return
      end
      subroutine p05_dat ( data_num, x, y )

c*********************************************************************72
c
cc P05_DAT returns the data vector for problem 5.
c
c  Discussion:
c
c    This example is due to deBoor.
c    This data represents a property of titanium as a function of temperature.
c
c    The X data is equally spaced.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Carl deBoor,
c    A Practical Guide to Splines,
c    Springer Verlag.
c
c  Parameters:
c
c    Input, integer DATA_NUM, the dimension of the data vector.
c
c    Output, double precision X(DATA_NUM), the abscissa data.
c
c    Output, double precision Y(DATA_NUM), the ordinate data.
c
      implicit none

      integer data_num

      double precision x(data_num)
      double precision x_save(49)
      double precision y(data_num)
      double precision y_save(49)

      data x_save /
     &  595.0D+00, 
     &  605.0D+00,  615.0D+00,  625.0D+00,  635.0D+00,  645.0D+00, 
     &  655.0D+00,  665.0D+00,  675.0D+00,  685.0D+00,  695.0D+00, 
     &  705.0D+00,  715.0D+00,  725.0D+00,  735.0D+00,  745.0D+00, 
     &  755.0D+00,  765.0D+00,  775.0D+00,  785.0D+00,  795.0D+00, 
     &  805.0D+00,  815.0D+00,  825.0D+00,  835.0D+00,  845.0D+00, 
     &  855.0D+00,  865.0D+00,  875.0D+00,  885.0D+00,  895.0D+00, 
     &  905.0D+00,  915.0D+00,  925.0D+00,  935.0D+00,  945.0D+00, 
     &  955.0D+00,  965.0D+00,  975.0D+00,  985.0D+00,  995.0D+00, 
     & 1005.0D+00, 1015.0D+00, 1025.0D+00, 1035.0D+00, 1045.0D+00, 
     & 1055.0D+00, 1065.0D+00, 1075.0D+00 /
      data y_save /
     &  0.644D+00, 
     &  0.622D+00, 0.638D+00, 0.649D+00, 0.652D+00, 0.639D+00, 
     &  0.646D+00, 0.657D+00, 0.652D+00, 0.655D+00, 0.644D+00, 
     &  0.663D+00, 0.663D+00, 0.668D+00, 0.676D+00, 0.676D+00, 
     &  0.686D+00, 0.679D+00, 0.678D+00, 0.683D+00, 0.694D+00, 
     &  0.699D+00, 0.710D+00, 0.730D+00, 0.763D+00, 0.812D+00, 
     &  0.907D+00, 1.044D+00, 1.336D+00, 1.881D+00, 2.169D+00, 
     &  2.075D+00, 1.598D+00, 1.211D+00, 0.916D+00, 0.746D+00, 
     &  0.672D+00, 0.627D+00, 0.615D+00, 0.607D+00, 0.606D+00, 
     &  0.609D+00, 0.603D+00, 0.601D+00, 0.603D+00, 0.601D+00, 
     &  0.611D+00, 0.601D+00, 0.608D+00 /

      call r8vec_copy ( data_num, x_save, x )
      call r8vec_copy ( data_num, y_save, y )

      return
      end
      subroutine p05_data_num ( data_num )

c*********************************************************************72
c
cc P05_DATA_NUM returns the dimension of the data vector for problem 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer DATA_NUM, the dimension of the data vector.
c
      implicit none

      integer data_num

      data_num = 49

      return
      end
      subroutine p05_story ( )

c*********************************************************************72
c
cc P05_STORY prints the "story" for problem 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
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

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  This example is due to deBoor.'
      write ( *, '(a)' ) 
     &  '  This data represents a temperature dependent property '
      write ( *, '(a)' ) '  of titanium.'

      return
      end
      subroutine p05_title ( title )

c*********************************************************************72
c
cc P05_TITLE returns the title of problem 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'DeBoor''s Titanium property data'

      return
      end
      subroutine p06_dat ( data_num, x, y )

c*********************************************************************72
c
cc P06_DAT returns the data vector for problem 6.
c
c  Discussion:
c
c    The X data is equally spaced.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DATA_NUM, the dimension of the data vector.
c
c    Output, double precision X(DATA_NUM), the abscissa data.
c
c    Output, double precision Y(DATA_NUM), the ordinate data.
c
      implicit none

      integer data_num

      integer i
      integer j
      integer n
      integer num_int
      parameter ( num_int = 5 )
      double precision x(data_num)
      double precision y(data_num)

      n = 1
      x(n) = 0.0D+00
      y(n) = 0.0D+00

      do i = 1, num_int

        do j = 1, i
          n = n + 1
          x(n) = dble ( i - 1 ) + 0.5D+00 * dble ( j ) / dble ( i )
          y(n) = dble ( j ) / dble ( i )
        end do

        do j = 1, i
          n = n + 1
          x(n) = dble ( i - 1 ) + 0.5D+00 
     &      + 0.5D+00 * dble ( j ) / dble ( i )
          y(n) = 1.0D+00 - dble ( j ) / dble ( i )
        end do

      end do

      return
      end
      subroutine p06_data_num ( data_num )

c*********************************************************************72
c
cc P06_DATA_NUM returns the dimension of the data vector for problem 6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer DATA_NUM, the dimension of the data vector.
c
      implicit none

      integer data_num
      integer num_int
      parameter ( num_int = 5 )

      data_num = 1 + num_int * ( num_int + 1 )

      return
      end
      subroutine p06_story ( )

c*********************************************************************72
c
cc P06_STORY prints the "story" for problem 6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
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

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  This is a data vector.'

      return
      end
      subroutine p06_title ( title )

c*********************************************************************72
c
cc P06_TITLE returns the title of problem 6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'The Sawtooth data'

      return
      end
      subroutine p07_dat ( data_num, x, y )

c*********************************************************************72
c
cc P07_DAT returns the data vector for problem 7.
c
c  Discussion:
c
c    The X data is NOT equally spaced.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DATA_NUM, the dimension of the data vector.
c
c    Output, double precision X(DATA_NUM), the abscissa data.
c
c    Output, double precision Y(DATA_NUM), the ordinate data.
c
      implicit none

      integer data_num

      double precision x(data_num)
      double precision x_save(9)
      double precision y(data_num)
      double precision y_save(9)

      data x_save /
     &            0.0D+00,  0.1D+00,  0.2D+00,  0.3D+00, 0.4D+00, 
     &            0.5D+00,  0.6D+00,  0.8D+00,  1.0D+00 /
      data y_save /
     &            0.0D+00,  0.9D+00,  0.95D+00, 0.9D+00, 0.1D+00, 
     &            0.05D+00, 0.05D+00, 0.2D+00,  1.0D+00 /

      call r8vec_copy ( data_num, x_save, x )
      call r8vec_copy ( data_num, y_save, y )

      return
      end
      subroutine p07_data_num ( data_num )

c*********************************************************************72
c
cc P07_DATA_NUM returns the dimension of the data vector for problem 7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer DATA_NUM, the dimension of the data vector.
c
      implicit none

      integer data_num

      data_num = 9

      return
      end
      subroutine p07_story ( )

c*********************************************************************72
c
cc P07_STORY prints the "story" for problem 7
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
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

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  This is a data vector.'

      return
      end
      subroutine p07_title ( title )

c*********************************************************************72
c
cc P07_TITLE returns the title of problem 7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'Concavity test data'

      return
      end
      subroutine p08_dat ( data_num, x, y )

c*********************************************************************72
c
cc P08_DAT returns the data vector for problem 8.
c
c  Discussion:
c
c    This example is due to Pierre Blais.
c
c    Data is only available over the interval [0, 238], but extrapolation
c    is to be used to extend the approximate function to a maximum argument
c    of 1023.  The behavior of the extrapolated curve is of great interest.
c
c    The X data is NOT equally spaced.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DATA_NUM, the dimension of the data vector.
c
c    Output, double precision X(DATA_NUM), the abscissa data.
c
c    Output, double precision Y(DATA_NUM), the ordinate data.
c
      implicit none

      integer data_num

      double precision x(data_num)
      double precision x_save(12)
      double precision y(data_num)
      double precision y_save(12)

      data x_save /
     &    0.0D+00,  71.0D+00,  104.0D+00,  135.0D+00, 145.0D+00, 
     &  160.0D+00, 181.0D+00,  193.0D+00,  205.0D+00, 215.0D+00, 
     &  225.0D+00, 238.0D+00 /
      data y_save /
     &    0.0000D+00,
     &    7.7554D+00,
     &   19.7062D+00,
     &   35.53786D+00,
     &   42.91537D+00, 
     &   54.7752D+00,
     &   66.75865D+00,
     &   78.49286D+00,
     &   89.76833D+00,
     &  101.746D+00, 
     &  113.4824D+00,
     &  135.4566D+00 /

      call r8vec_copy ( data_num, x_save, x )
      call r8vec_copy ( data_num, y_save, y )

      return
      end
      subroutine p08_data_num ( data_num )

c*********************************************************************72
c
cc P08_DATA_NUM returns the dimension of the data vector for problem 8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer DATA_NUM, the dimension of the data vector.
c
      implicit none

      integer data_num

      data_num = 12

      return
      end
      subroutine p08_story ( )

c*********************************************************************72
c
cc P08_STORY prints the "story" for problem 8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
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

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  This example is due to Pierre Blais.'
      write ( *, '(a)' ) 
     &  '  Data is only available over the interval [0, 238], but '
      write ( *, '(a)' )
     &  '  extrapolation is to be used to extend the approximate '
      write ( *, '(a)' )
     &  '  function to a maximum argument of 1023.  The behavior of '
      write ( *, '(a)' )
     &  '  the extrapolated curve is of great interest.'

      return
      end
      subroutine p08_title ( title )

c*********************************************************************72
c
cc P08_TITLE returns the title of problem 8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'Extrapolation test data'

      return
      end
      subroutine p09_dat ( data_num, x, y )

c*********************************************************************72
c
cc P09_DAT returns the data vector for problem 9.
c
c  Discussion:
c
c    This example is due to Max Waldmeier.
c
c    This data represents a measure of sunspot activity over the
c    years 1700 to 1960.  The X value is the year, and the Y value
c    is a measure of the sunspot activity, which is usually, but
c    not always, an integer.
c
c    The X data is equally spaced.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Max Waldmeier,
c    The Sunspot-Activity in the Years 1610-1960,
c    Shulthess, Zurich, 1961.
c
c  Parameters:
c
c    Input, integer DATA_NUM, the dimension of the data vector.
c
c    Output, double precision X(DATA_NUM), the abscissa data.
c
c    Output, double precision Y(DATA_NUM), the ordinate data.
c
      implicit none

      integer data_num

      integer i
      double precision x(data_num)
      double precision y(data_num)
      double precision y_save(261)

      do i = 1, 261
        x(i) = dble ( 1699 + i )
      end do

      data y_save /
     &   5.0D+00,  11.0D+00,  16.0D+00,  23.0D+00,  36.0D+00, 
     &  58.0D+00,  29.0D+00,  20.0D+00,  10.0D+00,   8.0D+00, 
     &   3.0D+00,   0.0D+00,   0.0D+00,   2.0D+00,  11.0D+00, 
     &  27.0D+00,  47.0D+00,  63.0D+00,  60.0D+00,  39.0D+00, 
     &  28.0D+00,  26.0D+00,  22.0D+00,  11.0D+00,  21.0D+00, 
     &  40.0D+00,  78.0D+00, 122.0D+00, 103.0D+00,  73.0D+00, 
     &  47.0D+00,  35.0D+00,  11.0D+00,   5.0D+00,  16.0D+00, 
     &  34.0D+00,  70.0D+00,  81.0D+00, 111.0D+00, 101.0D+00, 
     &  73.0D+00,  40.0D+00,  20.0D+00,  16.0D+00,   5.0D+00, 
     &  11.0D+00,  22.0D+00,  40.0D+00,  60.0D+00,  80.9D+00, 
     &  83.4D+00,  47.7D+00,  47.8D+00,  30.7D+00,  12.2D+00, 
     &   9.6D+00,  10.2D+00,  32.4D+00,  47.6D+00,  54.0D+00, 
     &  62.9D+00,  85.9D+00,  61.2D+00,  45.1D+00,  36.4D+00, 
     &  20.9D+00,  11.4D+00,  37.8D+00,  69.8D+00, 106.1D+00, 
     & 100.8D+00,  81.6D+00,  66.5D+00,  34.8D+00,  30.6D+00, 
     &   7.0D+00,  19.8D+00,  92.5D+00, 154.4D+00, 125.9D+00, 
     &  84.8D+00,  68.1D+00,  38.5D+00,  22.8D+00,  10.2D+00, 
     &  24.1D+00,  82.9D+00, 132.0D+00, 130.9D+00, 118.1D+00, 
     &  89.9D+00,  66.6D+00,  60.0D+00,  46.9D+00,  41.0D+00, 
     &  21.3D+00,  16.0D+00,   6.4D+00,   4.1D+00,   6.8D+00, 
     &  14.5D+00,  34.0D+00,  45.0D+00,  43.1D+00,  47.5D+00, 
     &  42.2D+00,  28.1D+00,  10.1D+00,   8.1D+00,   2.5D+00, 
     &   0.0D+00,   1.4D+00,   5.0D+00,  12.2D+00,  13.9D+00, 
     &  35.4D+00,  45.8D+00,  41.1D+00,  30.1D+00,  23.9D+00, 
     &  15.6D+00,   6.6D+00,   4.0D+00,   1.8D+00,   8.5D+00, 
     &  16.6D+00,  36.3D+00,  49.6D+00,  64.2D+00,  67.0D+00, 
     &  70.9D+00,  47.8D+00,  27.5D+00,   8.5D+00,  13.2D+00, 
     &  56.9D+00, 121.5D+00, 138.3D+00, 103.2D+00,  85.7D+00, 
     &  64.6D+00,  36.7D+00,  24.2D+00,  10.7D+00,  15.0D+00, 
     &  40.1D+00,  61.5D+00,  98.5D+00, 124.7D+00,  96.3D+00, 
     &  66.6D+00,  64.5D+00,  54.1D+00,  39.0D+00,  20.6D+00, 
     &   6.7D+00,   4.3D+00,  22.7D+00,  54.8D+00,  93.8D+00, 
     &  95.8D+00,  77.2D+00,  59.1D+00,  44.0D+00,  47.0D+00, 
     &  30.5D+00,  16.3D+00,   7.3D+00,  37.6D+00,  74.0D+00, 
     & 139.0D+00, 111.2D+00, 101.6D+00,  66.2D+00,  44.7D+00, 
     &  17.0D+00,  11.3D+00,  12.4D+00,   3.4D+00,   6.0D+00, 
     &  32.3D+00,  54.3D+00,  59.7D+00,  63.7D+00,  63.5D+00, 
     &  52.2D+00,  25.4D+00,  13.1D+00,   6.8D+00,   6.3D+00, 
     &   7.1D+00,  35.6D+00,  73.0D+00,  85.1D+00,  78.0D+00, 
     &  64.0D+00,  41.8D+00,  26.2D+00,  26.7D+00,  12.1D+00, 
     &   9.5D+00,   2.7D+00,   5.0D+00,  24.4D+00,  42.0D+00, 
     &  63.5D+00,  53.8D+00,  62.0D+00,  48.5D+00,  43.9D+00, 
     &  18.6D+00,   5.7D+00,   3.6D+00,   1.4D+00,   9.6D+00, 
     &  47.4D+00,  57.1D+00, 103.9D+00,  80.6D+00,  63.6D+00, 
     &  37.6D+00,  26.1D+00,  14.2D+00,   5.8D+00,  16.7D+00, 
     &  44.3D+00,  63.9D+00,  69.0D+00,  77.8D+00,  64.9D+00, 
     &  35.7D+00,  21.2D+00,  11.1D+00,   5.7D+00,   8.7D+00, 
     &  36.1D+00,  79.7D+00, 114.4D+00, 109.6D+00,  88.8D+00, 
     &  67.8D+00,  47.5D+00,  30.6D+00,  16.3D+00,   9.6D+00, 
     &  33.2D+00,  92.6D+00, 151.6D+00, 136.3D+00, 134.7D+00, 
     &  83.9D+00,  69.4D+00,  31.5D+00,  13.9D+00,   4.4D+00, 
     &  38.0D+00, 141.7D+00, 190.2D+00, 184.8D+00, 159.0D+00, 
     & 112.3D+00 /

      call r8vec_copy ( data_num, y_save, y )

      return
      end
      subroutine p09_data_num ( data_num )

c*********************************************************************72
c
cc P09_DATA_NUM returns the dimension of the data vector for problem 9.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer DATA_NUM, the dimension of the data vector.
c
      implicit none

      integer data_num

      data_num = 261

      return
      end
      subroutine p09_story ( )

c*********************************************************************72
c
cc P09_STORY prints the "story" for problem 09
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
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

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  This example is due to Max Waldmeier.'
      write ( *, '(a)' ) 
     &  '  This data represents a measure of sunspot activity over the'
      write ( *, '(a)' ) 
     &  '  years 1700 to 1960.  The X value is the year, and the Y'
      write ( *, '(a)' ) 
     &  '  value is a measure of the sunspot activity, which is '
      write ( *, '(a)' ) 
     &  '  usually, but not always, an integer.'

      return
      end
      subroutine p09_title ( title )

c*********************************************************************72
c
cc P09_TITLE returns the title of problem 09.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'Sunspot data, 1700-1960.'

      return
      end
      subroutine p10_dat ( data_num, x, y )

c*********************************************************************72
c
cc P10_DAT returns the data vector for problem 10.
c
c  Discussion:
c
c    100 uniformly random X values between -2 and 5 were selected,
c    and the formula Y = 2 + 5 * X + 10 * N(0,1) was evaluated, where
c    N(0,1) represents random normal values with 0 mean and unit variance.
c
c    The original data was unsorted, but this caused problems for various
c    approximation codes, so the data has now been sorted.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DATA_NUM, the dimension of the data vector.
c
c    Output, double precision X(DATA_NUM), the abscissa data.
c
c    Output, double precision Y(DATA_NUM), the ordinate data.
c
      implicit none

      integer data_num

      double precision x(data_num)
      double precision x_save(100)
      double precision y(data_num)
      double precision y_save(100)

      save x_save
      save y_save

      data x_save /
     & -1.9935D+00, -1.8884D+00, -1.8813D+00, -1.8786D+00, -1.8714D+00, 
     & -1.8108D+00, -1.6614D+00, -1.5991D+00, -1.4423D+00, -1.3441D+00, 
     & -1.1730D+00, -1.1575D+00, -0.9886D+00, -0.8196D+00, -0.7771D+00, 
     & -0.7706D+00, -0.7270D+00, -0.7243D+00, -0.7228D+00, -0.6235D+00, 
     & -0.6023D+00, -0.1650D+00, -0.1117D+00, -0.0570D+00,  0.0923D+00, 
     &  0.1132D+00,  0.3157D+00,  0.3343D+00,  0.3356D+00,  0.3827D+00, 
     &  0.4019D+00,  0.4570D+00,  0.4866D+00,  0.6041D+00,  0.6803D+00, 
     &  0.8395D+00,  0.8734D+00,  0.9115D+00,  0.9190D+00,  1.0716D+00, 
     &  1.1649D+00,  1.2261D+00,  1.3118D+00,  1.3484D+00,  1.3698D+00, 
     &  1.5721D+00,  1.6258D+00,  1.7749D+00,  1.8234D+00,  1.9732D+00, 
     &  1.9784D+00,  2.0022D+00,  2.0231D+00,  2.0398D+00,  2.0401D+00, 
     &  2.0757D+00,  2.0824D+00,  2.2051D+00,  2.2634D+00,  2.3076D+00, 
     &  2.3660D+00,  2.3931D+00,  2.4175D+00,  2.5928D+00,  2.6235D+00, 
     &  2.8109D+00,  2.8898D+00,  2.9561D+00,  2.9669D+00,  3.1161D+00, 
     &  3.1430D+00,  3.1749D+00,  3.2373D+00,  3.4166D+00,  3.6131D+00, 
     &  3.6635D+00,  3.6669D+00,  3.7486D+00,  3.7530D+00,  3.7862D+00, 
     &  3.8021D+00,  3.8298D+00,  3.8696D+00,  3.9051D+00,  3.9509D+00, 
     &  4.0678D+00,  4.1038D+00,  4.1557D+00,  4.2456D+00,  4.2564D+00, 
     &  4.3228D+00,  4.3403D+00,  4.4157D+00,  4.4998D+00,  4.5767D+00, 
     &  4.6312D+00,  4.7029D+00,  4.7717D+00,  4.8424D+00,  4.8915D+00 /

      data y_save /
     & -1.8918D+00, -8.6202D+00, -0.4152D+00, -4.6967D+00, -2.4144D+00, 
     &-21.8853D+00,-16.5097D+00,-10.4657D+00, -4.1150D+00,  6.5665D+00, 
     & -6.7649D+00,  8.8276D+00,  1.8109D+00,  9.6428D+00, -0.6165D+00, 
     & -8.4213D+00,-16.4491D+00, -0.0667D+00,  6.5713D+00, -4.0436D+00, 
     & -6.4194D+00, -1.9114D+00, -9.5246D+00, -3.2154D+00,  0.6541D+00, 
     &  3.0247D+00,  2.9410D+00,  9.7848D+00,  4.7713D+00, 22.0539D+00, 
     &  7.1301D+00, 22.3302D+00, -2.7979D+00, 10.2864D+00,  2.7993D+00, 
     & 12.1989D+00, 12.3065D+00,-15.3026D+00, -6.6749D+00, -7.0521D+00, 
     & 11.8429D+00, 22.8325D+00,  5.2912D+00, 16.8654D+00, 14.3045D+00, 
     & -0.6553D+00, 14.1038D+00,  3.3556D+00, 26.2797D+00, 11.5406D+00, 
     & 28.2522D+00,  7.7605D+00, 18.0102D+00, 11.5710D+00, -8.0188D+00, 
     &  2.5576D+00, 18.5371D+00, 12.4771D+00,  2.1298D+00,  7.2745D+00, 
     & 16.3256D+00,  4.0354D+00, 23.8374D+00,  8.5572D+00, 33.2063D+00, 
     &  5.2561D+00, 18.4409D+00,  1.5702D+00,  9.5982D+00, 11.6483D+00, 
     & 21.7285D+00, 27.2958D+00, 21.1916D+00, 15.3527D+00, 28.2207D+00, 
     & 28.3065D+00, 21.5368D+00, 26.4556D+00, 24.8933D+00, 11.0617D+00, 
     & 28.6064D+00, 14.5770D+00, 15.3088D+00, 23.2952D+00, 18.6796D+00, 
     & 21.0208D+00, 28.4730D+00, 33.2469D+00, 21.2484D+00, 26.5592D+00, 
     & 21.2314D+00, 25.9979D+00, 28.4785D+00, 18.3307D+00, 27.6318D+00, 
     & 31.1673D+00, 26.4379D+00, 43.1573D+00, 20.1264D+00, 19.0873D+00 /

      call r8vec_copy ( data_num, x_save, x )
      call r8vec_copy ( data_num, y_save, y )

      return
      end
      subroutine p10_data_num ( data_num )

c*********************************************************************72
c
cc P10_DATA_NUM returns the dimension of the data vector for problem 10.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer DATA_NUM, the dimension of the data vector.
c
      implicit none

      integer data_num

      data_num = 100

      return
      end
      subroutine p10_story ( )

c*********************************************************************72
c
cc P10_STORY prints the "story" for problem 10.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
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

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  100 uniformly random X values between -2 and 5 were'
      write ( *, '(a)' )
     &  '  selected, and the formula Y = 2 + 5 * X + 10 * N(0,1) was'
      write ( *, '(a)' )
     &  '  evaluated, where N(0,1) represents random normal values ' 
      write ( *, '(a)' )
     &  '  with 0 mean and unit variance.'

      return
      end
      subroutine p10_title ( title )

c*********************************************************************72
c
cc P10_TITLE returns the title of problem 10.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'Y = 2 + 5*X + 10*N(0,1).'

      return
      end
      subroutine r8vec_copy ( n, a1, a2 )

c*********************************************************************72
c
cc R8VEC_COPY copies an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8 values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the length of the vectors.
c
c    Input, double precision A1(N), the vector to be copied.
c
c    Output, double precision A2(N), a copy of A1.
c
      implicit none

      integer n

      double precision a1(n)
      double precision a2(n)
      integer i

      do i = 1, n
        a2(i) = a1(i)
      end do

      return
      end
      subroutine r8vec2_print ( n, a1, a2, title )

c*********************************************************************72
c
cc R8VEC2_PRINT prints an R8VEC2.
c
c  Discussion:
c
c    An R8VEC2 is a dataset consisting of N pairs of R8s, stored
c    as two separate vectors A1 and A2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of components of the vector.
c
c    Input, double precision A1(N), A2(N), the vectors to be printed.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer n

      double precision a1(n)
      double precision a2(n)
      integer i
      character ( len = * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a1(i), a2(i)
      end do

      return
      end
      subroutine r8vec2_write ( output_filename, n, x, y )

c*********************************************************************72
c
cc R8VEC2_WRITE writes an R8VEC2 file.
c
c  Discussion:
c
c    An R8VEC2 is a pair of vectors of R8 values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 July 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character ( len = * ) OUTPUT_FILENAME, the output file name.
c
c    Input, integer N, the number of points.
c
c    Input, double precision X(N), Y(N), the data.
c
      implicit none

      integer n

      integer j
      character * ( * ) output_filename
      integer output_status
      integer output_unit
      double precision x(n)
      double precision y(n)
c
c  Open the file.
c
      call get_unit ( output_unit )

      open ( unit = output_unit, file = output_filename, 
     &  status = 'replace', iostat = output_status )

      if ( output_status .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8VEC_WRITE - Fatal error!'
        write ( *, '(a,i8)' ) '  Could not open the output file "' // 
     &    trim ( output_filename ) // '" on unit ', output_unit
        output_unit = -1
        stop
      end if

      if ( 0 .lt. n ) then
c
c  Write the data.
c
        do j = 1, n
          write ( output_unit, '(2x,g24.16,2x,g24.16)' ) x(j), y(j)
        end do

      end if
c
c  Close the file.
c
      close ( unit = output_unit )

      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Discussion:
c
c    This FORTRAN77 version is made available for cases where the
c    FORTRAN90 version cannot be used.
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
