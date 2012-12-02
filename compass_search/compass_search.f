      subroutine compass_search ( function_handle, m, x0, delta_tol, 
     &  delta_init, k_max, x, fx, k )

c*********************************************************************72
c
cc COMPASS_SEARCH carries out a direct search minimization algorithm.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Tamara Kolda, Robert Michael Lewis, Virginia Torczon,
c    Optimization by Direct Search: New Perspectives on Some Classical 
c    and Modern Methods,
c    SIAM Review,
c    Volume 45, Number 3, 2003, pages 385-482. 
c
c  Parameters:
c
c    Input, external double precision FUNCTION_HANDLE, the name of
c    a FORTRAN90 function which evaluates the function to be minimized, of the
c    form FUNCTION FUNCTION_HANDLE ( M, X ).
c
c    Input, integer M, the number of variables.
c
c    Input, double precision X0(M), a starting estimate for the minimizer.
c
c    Input, double precision DELTA_TOL, the smallest step size that is allowed.
c
c    Input, double precision DELTA_INIT, the starting stepsize.  
c
c    Input, integer K_MAX, the maximum number of steps allowed.
c
c    Output, double precision X(M), the estimated minimizer.
c
c    Output, double precision FX, the function value at X.
c
c    Output, integer K, the number of steps taken.
c
      implicit none

      integer m

      logical decrease
      double precision delta
      double precision delta_init
      double precision delta_tol
      double precision function_handle
      external function_handle
      double precision fx
      double precision fxd
      integer i
      integer ii
      integer k
      integer k_max
      double precision s
      double precision x(m)
      double precision x0(m)
      double precision xd(m)

      k = 0
      do i = 1, m
        x(i) = x0(i)
      end do
      fx = function_handle ( m, x )

      if ( delta_tol .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'COMPASS_SEARCH - Fatal error!'
        write ( *, '(a)' ) '  DELTA_TOL <= 0.0.'
        write ( *, '(a,g14.6)' ) '  DELTA_TOL = ', delta_tol
        stop
      end if

      if ( delta_init .le. delta_tol ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'COMPASS_SEARCH - Fatal error!'
        write ( *, '(a)' ) '  DELTA_INIT < DELTA_TOL.'
        write ( *, '(a,g14.6)' ) '  DELTA_INIT = ', delta_init
        write ( *, '(a,g14.6)' ) '  DELTA_TOL = ', delta_tol
        stop
      end if

      delta = delta_init

10    continue

      if ( k .lt. k_max ) then

        k = k + 1
c
c  For each coordinate direction I, seek a lower function value
c  by increasing or decreasing X(I) by DELTA.
c
        decrease = .false.
        s = + 1.0D+00
        i = 1

        do ii = 1, 2 * m

          xd = x
          xd(i) = xd(i) + s * delta
          fxd = function_handle ( m, xd )
c
c  As soon as a decrease is noticed, accept the new point.
c
          if ( fxd .lt. fx ) then
            x = xd
            fx = fxd
            decrease = .true.
            go to 20
          end if

          s = - s
          if ( s == + 1.0D+00 ) then
            i = i + 1
          end if

        end do

20      continue
c
c  If no decrease occurred, reduce DELTA.
c
        if ( .not. decrease ) then
          delta = delta / 2.0D+00
          if ( delta .lt. delta_tol ) then
            go to 30
          end if
        end if

        go to 10

      end if

30    continue

      return
      end
      subroutine r8vec_print ( n, a, title )

c*********************************************************************72
c
cc R8VEC_PRINT prints an R8VEC.
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
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of components of the vector.
c
c    Input, double precision A(N), the vector to be printed.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      character ( len = * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
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
