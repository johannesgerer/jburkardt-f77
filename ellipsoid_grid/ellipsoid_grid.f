      subroutine ellipsoid_grid ( n, r, c, ng, xyz )

c*********************************************************************72
c
cc ELLIPSOID_GRID generates the grid points inside an ellipsoid.
c
c  Discussion:
c
c    The ellipsoid is specified as
c
c      ( ( X - C1 ) / R1 )^2 
c    + ( ( Y - C2 ) / R2 )^2 
c    + ( ( Z - C3 ) / R3 )^2 = 1
c
c    The user supplies a number N.  There will be N+1 grid points along
c    the shortest axis.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 November 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of subintervals.
c
c    Input, double precision R(3), the half axis lengths.
c
c    Input, double precision C(3), the center of the ellipsoid.
c
c    Input, integer NG, the number of grid points.
c
c    Output, double precision XYZ(3,NG), the grid point coordinates.
c
      implicit none

      integer ng

      double precision c(3)
      double precision h
      integer i
      integer i4_ceiling
      integer ii
      integer j
      integer jj
      integer k
      integer m
      integer n
      integer ni
      integer nj
      integer nk
      integer np
      double precision p(3,8)
      double precision r(3)
      double precision rmin
      double precision x
      double precision xyz(3,ng)
      double precision y
      double precision z

      call r8vec_min ( 3, r, rmin )

      if ( r(1) .eq. rmin ) then
        h = 2.0D+00 * r(1) / dble ( 2 * n + 1 )
        ni = n
        nj = i4_ceiling ( r(2) / r(1) ) * dble ( n )
        nk = i4_ceiling ( r(3) / r(1) ) * dble ( n )
      else if ( r(2) .eq. rmin ) then
        h = 2.0D+00 * r(2) / dble ( 2 * n + 1 )
        nj = n
        ni = i4_ceiling ( r(1) / r(2) ) * dble ( n )
        nk = i4_ceiling ( r(3) / r(2) ) * dble ( n )
      else
        h = 2.0D+00 * r(3) / dble ( 2 * n + 1 )
        nk = n
        ni = i4_ceiling ( r(1) / r(3) ) * dble ( n )
        nj = i4_ceiling ( r(2) / r(3) ) * dble ( n )
      end if

      ng2 = 0

      do k = 0, nk
        z = c(3) + dble ( k ) * h
        do j = 0, nj
          y = c(2) + dble ( j ) * h
          do i = 0, ni
            x = c(1) + dble ( i ) * h
c
c  If we have left the ellipsoid, the I loop is completed.
c
            if ( 1.0D+00 .lt. ( ( x - c(1) ) / r(1) )**2 
     &                      + ( ( y - c(2) ) / r(2) )**2 
     &                      + ( ( z - c(3) ) / r(3) )**2 ) then
              exit
            end if
c
c  At least one point is generated, but more possible by symmetry.
c
            np = 1
            p(1,np) = x
            p(2,np) = y
            p(3,np) = z

            if ( 0 .lt. i ) then
              do m = 1, np
                p(1,m+np) = 2.0D+00 * c(1) - p(1,m)
                p(2,m+np) = p(2,m)
                p(3,m+np) = p(3,m)
              end do
              np = 2 * np
            end if

            if ( 0 .lt. j ) then
              do m = 1, np
                p(1,m+np) = p(1,m)
                p(2,m+np) = 2.0D+00 * c(2) - p(2,m)
                p(3,m+np) = p(3,m)
              end do
              np = 2 * np
            end if

            if ( 0 .lt. k ) then
              do m = 1, np
                p(1,m+np) = p(1,m)
                p(2,m+np) = p(2,m)
                p(3,m+np) = 2.0D+00 * c(3) - p(3,m)
              end do
              np = 2 * np
            end if

            do jj = 1, np
              do ii = 1, 3
                xyz(ii,ng2+jj) = p(ii,jj)
              end do
            end do

            ng2 = ng2 + np

          end do
        end do
      end do

      return
      end
      subroutine ellipsoid_grid_count ( n, r, c, ng )

c*****************************************************************************80
c
cc ELLIPSOID_GRID_COUNT counts the grid points inside an ellipsoid.
c
c  Discussion:
c
c    The ellipsoid is specified as
c
c      ( ( X - C1 ) / R1 )^2 
c    + ( ( Y - C2 ) / R2 )^2 
c    + ( ( Z - C3 ) / R3 )^2 = 1
c
c    The user supplies a number N.  There will be N+1 grid points along
c    the shortest axis.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 November 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of subintervals.
c
c    Input, double precision R(3), the half axis lengths.
c
c    Input, double precision C(3), the center of the ellipsoid.
c
c    Output, integer NG, the number of grid points.
c
      implicit none

      double precision c(3)
      double precision h
      integer i
      integer i4_ceiling
      integer j
      integer k
      integer n
      integer ng
      integer ni
      integer nj
      integer nk
      integer np
      double precision r(3)
      double precision rmin
      double precision x
      double precision y
      double precision z

      call r8vec_min ( 3, r, rmin )

      if ( r(1) .eq. rmin ) then
        h = 2.0D+00 * r(1) / dble ( 2 * n + 1 )
        ni = n
        nj = i4_ceiling ( r(2) / r(1) ) * dble ( n )
        nk = i4_ceiling ( r(3) / r(1) ) * dble ( n )
      else if ( r(2) .eq. rmin ) then
        h = 2.0D+00 * r(2) / dble ( 2 * n + 1 )
        nj = n
        ni = i4_ceiling ( r(1) / r(2) ) * dble ( n )
        nk = i4_ceiling ( r(3) / r(2) ) * dble ( n )
      else
        h = 2.0D+00 * r(3) / dble ( 2 * n + 1 )
        nk = n
        ni = i4_ceiling ( r(1) / r(3) ) * dble ( n )
        nj = i4_ceiling ( r(2) / r(3) ) * dble ( n )
      end if

      ng = 0

      do k = 0, nk
        z = c(3) + dble ( k ) * h
        do j = 0, nj
          y = c(2) + dble ( j ) * h
          do i = 0, ni
            x = c(1) + dble ( i ) * h
c
c  If we have left the ellipsoid, the I loop is completed.
c
            if ( 1.0D+00 .lt. ( ( x - c(1) ) / r(1) )**2 
     &                      + ( ( y - c(2) ) / r(2) )**2 
     &                      + ( ( z - c(3) ) / r(3) )**2 ) then
              exit
            end if
c
c  At least one point is generated, but more possible by symmetry.
c
            np = 1

            if ( 0 .lt. i ) then
              np = 2 * np
            end if

            if ( 0 .lt. j ) then
              np = 2 * np
            end if

            if ( 0 .lt. k ) then
              np = 2 * np
            end if

            ng = ng + np

          end do
        end do
      end do

      return
      end
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
      function i4_ceiling ( r )

c*********************************************************************72
c
cc I4_CEILING rounds an R8 "up" to the nearest I4.
c
c  Example:
c
c     R     Value
c
c    -1.1  -1
c    -1.0  -1
c    -0.9   0
c     0.0   0
c     5.0   5
c     5.1   6
c     5.9   6
c     6.0   6
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 November 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the value to be rounded up.
c
c    Output, integer I4_CEILING, the rounded value.
c
      implicit none

      double precision r
      integer i4_ceiling
      integer value

      value = int ( r )
      if ( dble ( value ) .lt. r ) then
        value = value + 1
      end if

      i4_ceiling = value

      return
      end
      subroutine r83vec_print_part ( n, a, max_print, title )

c*********************************************************************72
c
cc R83VEC_PRINT_PART prints "part" of an R83VEC.
c
c  Discussion:
c
c    The user specifies MAX_PRINT, the maximum number of lines to print.
c
c    If N, the size of the vector, is no more than MAX_PRINT, then
c    the entire vector is printed, one entry per line.
c
c    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
c    followed by a line of periods suggesting an omission,
c    and the last entry.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 November 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries of the vector.
c
c    Input, double precision A(3,N), the vector to be printed.
c
c    Input, integer MAX_PRINT, the maximum number of lines
c    to print.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer n

      double precision a(3,n)
      integer i
      integer max_print
      character * ( * )  title

      if ( max_print .le. 0 ) then
        return
      end if

      if ( n .le. 0 ) then
        return
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '

      if ( n .le. max_print ) then

        do i = 1, n
          write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &      i, ':', a(1,i), a(2,i), a(3,i)
        end do

      else if ( 3 .le. max_print ) then

        do i = 1, max_print - 2
          write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &      i, ':', a(1,i), a(2,i), a(3,i)
        end do
        write ( *, '(a)' ) 
     &    '  ........  ..............  ..............  ..............'
        i = n
        write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &    i, ':', a(1,i), a(2,i), a(3,i)

      else

        do i = 1, max_print - 1
          write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &      i, ':', a(1,i), a(2,i), a(3,i)
        end do
        i = max_print
        write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6,2x,g14.6,2x,a)' ) 
     &    i, ':', a(1,i), a(2,i), a(3,i), '...more entries...'

      end if

      return
      end
      subroutine r8mat_write ( output_filename, m, n, table )

c*********************************************************************72
c
cc R8MAT_WRITE writes a R8MAT file.
c
c  Discussion:
c
c    An R8MAT is an array of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 October 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) OUTPUT_FILENAME, the output file name.
c
c    Input, integer M, the spatial dimension.
c
c    Input, integer N, the number of points.
c
c    Input, double precision TABLE(M,N), the data.
c
      implicit none

      integer m
      integer n

      integer j
      character * ( * ) output_filename
      integer output_unit
      character * ( 30 ) string
      double precision table(m,n)
c
c  Open the file.
c
      call get_unit ( output_unit )

      open ( unit = output_unit, file = output_filename,
     &  status = 'replace' )
c
c  Create the format string.
c
      if ( 0 .lt. m .and. 0 .lt. n ) then

        write ( string, '(a1,i8,a1,i8,a1,i8,a1)' )
     &    '(', m, 'g', 24, '.', 16, ')'
c
c  Write the data.
c
        do j = 1, n
          write ( output_unit, string ) table(1:m,j)
        end do

      end if
c
c  Close the file.
c
      close ( unit = output_unit )

      return
      end
      subroutine r8vec_min ( n, a, amin )

c*********************************************************************72
c
cc R8VEC_MIN returns the minimum value in an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input, double precision A(N), the array.
c
c    Output, double precision AMIN, the value of the smallest entry.
c
      implicit none

      integer n

      double precision a(n)
      double precision amin
      integer i

      amin = a(1)
      do i = 2, n
        amin = min ( amin, a(i) )
      end do

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
