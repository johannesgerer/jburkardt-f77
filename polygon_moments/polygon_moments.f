      subroutine moment ( n, x, y, p, q, nu_pq )

c*********************************************************************72
c
cc MOMENT computes an unnormalized moment of a polygon.
c
c  Discussion:
c
c    Nu(P,Q) = Integral ( x, y in polygon ) x^p y^q dx dy
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 October 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Carsten Steger,
c    On the calculation of arbitrary moments of polygons,
c    Technical Report FGBV-96-05,
c    Forschungsgruppe Bildverstehen, Informatik IX,
c    Technische Universitaet Muenchen, October 1996.
c
c  Parameters:
c
c    Input, integer N, the number of vertices of the polygon.
c
c    Input, double precision X(N), Y(N), the vertex coordinates.
c
c    Input, integer P, Q, the indices of the moment.
c
c    Output, double precision NU_PQ, the unnormalized moment Nu(P,Q).
c
      implicit none

      integer n

      integer i
      integer k
      integer l
      double precision nu_pq
      integer p
      integer q
      double precision r8_choose
      double precision s_pq
      double precision x(n)
      double precision xi
      double precision xj
      double precision y(n)
      double precision yi
      double precision yj

      nu_pq = 0.0D+00

      xj = x(n)
      yj = y(n)

      do i = 1, n

        xi = x(i)
        yi = y(i)

        s_pq = 0.0D+00
        do k = 0, p
          do l = 0, q
            s_pq = s_pq 
     &        + r8_choose ( k + l, l ) 
     &        * r8_choose ( p + q - k - l, q - l ) 
     &        * xi ** k * xj ** ( p - k ) 
     &        * yi ** l * yj ** ( q - l )
          end do
        end do

        nu_pq = nu_pq + ( xj * yi - xi * yj ) * s_pq

        xj = xi
        yj = yi

      end do

      nu_pq = nu_pq / dble ( p + q + 2 ) 
     &  / dble ( p + q + 1 ) 
     &  / r8_choose ( p + q, p )

      return
      end
      subroutine moment_central ( n, x, y, p, q, mu_pq )

c*********************************************************************72
c
cc MOMENT_CENTRAL computes central moments of a polygon.
c
c  Discussion:
c
c    The central moment Mu(P,Q) is defined by
c
c      Mu(P,Q) = Integral ( polygon ) (x-Alpha(1,0))^p (y-Alpha(0,1))^q dx dy
c              / Area ( polygon )
c
c    where 
c
c      Alpha(1,0) = Integral ( polygon ) x dx dy / Area ( polygon )
c      Alpha(0,1) = Integral ( polygon ) y dx dy / Area ( polygon )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 October 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Carsten Steger,
c    On the calculation of arbitrary moments of polygons,
c    Technical Report FGBV-96-05,
c    Forschungsgruppe Bildverstehen, Informatik IX,
c    Technische Universitaet Muenchen, October 1996.
c
c  Parameters:
c
c    Input, integer N, the number of vertices of the polygon.
c
c    Input, double precision X(N), Y(N), the vertex coordinates.
c
c    Input, integer P, Q, the indices of the moment.
c
c    Output, double precision MU_PQ, the unnormalized moment Mu(P,Q).
c
      implicit none

      integer n

      double precision alpha_01
      double precision alpha_10
      double precision alpha_ij
      integer i
      integer j
      double precision mu_pq
      integer p
      integer q
      double precision r8_choose
      double precision r8_mop
      double precision x(n)
      double precision y(n)

      call moment_normalized ( n, x, y, 1, 0, alpha_10 )
      call moment_normalized ( n, x, y, 0, 1, alpha_01 )

      mu_pq = 0.0D+00

      do i = 0, p
        do j = 0, q

          call moment_normalized ( n, x, y, i, j, alpha_ij )

          mu_pq = mu_pq + r8_mop ( p + q - i - j ) 
     &      * r8_choose ( p, i ) * r8_choose ( q, j ) 
     &      * alpha_10 ** ( p - i ) * alpha_01 ** ( q - j ) * alpha_ij

        end do
      end do

      return
      end
      subroutine moment_normalized ( n, x, y, p, q, alpha_pq )

c*********************************************************************72
c
cc MOMENT_NORMALIZED computes a normalized moment of a polygon.
c
c  Discussion:
c
c    Alpha(P,Q) = Integral ( x, y in polygon ) x^p y^q dx dy / Area ( polygon )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 October 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Carsten Steger,
c    On the calculation of arbitrary moments of polygons,
c    Technical Report FGBV-96-05,
c    Forschungsgruppe Bildverstehen, Informatik IX,
c    Technische Universitaet Muenchen, October 1996.
c
c  Parameters:
c
c    Input, integer N, the number of vertices of the polygon.
c
c    Input, double precision X(N), Y(N), the vertex coordinates.
c
c    Input, integer P, Q, the indices of the moment.
c
c    Output, double precision ALPHA_PQ, the normalized moment Alpha(P,Q).
c
      implicit none

      integer n

      double precision alpha_pq
      double precision nu_00
      double precision nu_pq
      integer p
      integer q
      double precision x(n)
      double precision y(n)

      call moment ( n, x, y, p, q, nu_pq )
      call moment ( n, x, y, 0, 0, nu_00 )

      alpha_pq = nu_pq / nu_00

      return
      end
      function r8_choose ( n, k )

c*********************************************************************72
c
cc R8_CHOOSE computes the binomial coefficient C(N,K) as an R8.
c
c  Discussion:
c
c    The value is calculated in such a way as to avoid overflow and
c    roundoff.  The calculation is done in R8 arithmetic.
c
c    The formula used is:
c
c      C(N,K) = N! / ( K! * (N-K)! )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 June 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    ML Wolfson, HV Wright,
c    Algorithm 160:
c    Combinatorial of M Things Taken N at a Time,
c    Communications of the ACM,
c    Volume 6, Number 4, April 1963, page 161.
c
c  Parameters:
c
c    Input, integer N, K, are the values of N and K.
c
c    Output, double precision R8_CHOOSE, the number of combinations of N
c    things taken K at a time.
c
      implicit none

      integer i
      integer k
      integer mn
      integer mx
      integer n
      double precision r8_choose
      double precision value

      mn = min ( k, n - k )

      if ( mn .lt. 0 ) then

        value = 0.0D+00

      else if ( mn .eq. 0 ) then

        value = 1.0D+00

      else

        mx = max ( k, n - k )
        value = dble ( mx + 1 )

        do i = 2, mn
          value = ( value * dble ( mx + i ) ) / dble ( i )
        end do

      end if

      r8_choose = value

      return
      end
      function r8_mop ( i )

c*********************************************************************72
c
cc R8_MOP returns the I-th power of -1 as an R8.
c
c  Discussion:
c
c    An R8 is a double precision real value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the power of -1.
c
c    Output, double precision R8_MOP, the I-th power of -1.
c
      implicit none

      integer i
      double precision r8_mop

      if ( mod ( i, 2 ) .eq. 0 ) then
        r8_mop = + 1.0D+00
      else
        r8_mop = - 1.0D+00
      end if

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
