      subroutine pdeone ( t, u, udot, npde, npts, alpha, beta, dval,
     &  gamma, uavg, ux )

c*********************************************************************72
c
cc PDEONE converts a 1D PDE to a system of ODE's.
c
c  Discussion:
c
c    PDEONE is an interface subroutine which uses centered difference
c    approximations to convert one-dimensional systems of partial
c    differential equations into a system of ordinary differential equations:
c
c      UDOT = F(T,X,U).
c
c    This routine is intended to be used with a robust ODE integrator.
c
c  Modified:
c
c    04 February 2011
c
c  Author:
c
c    Richard Sincovec, Niel Madsen
c
c  Reference:
c
c    Richard Sincovec, Niel Madsen,
c    Software for Nonlinear Partial Differential Equations,
c    ACM Transactions on Mathematical Software,
c    Volume 1, Number 3, September 1975, pages 232-260.
c
c    Richard Sincovec, Niel Madsen,
c    Algorithm 494: PDEONE,
c    ACM Transactions on Mathematical Software,
c    Volume 1, Number 3, September 1975, pages 261-263.
c
c  Parameters:
c
c    Input, real T, the current value of time.
c
c    Input, real U(NPDE,NPTS), the computed solution at time T.
c
c    Output, real UDOT(NPDE,NPTS), the right hand side of the resulting
c    system of ODEs, F(T,X,U), obtained by discretizing the given PDE's.
c
c    Input, input NPDE, the number of partial differential equations.
c
c    Input, input NPTS, the number of spatial grid points.
c
      integer npde
      integer npts

      real alpha(npde)
      real beta(npde)
      real dval(npde,npde,2)
      real dxi
      real dxic
      real dxil
      real dxir
      real gamma(npde)
      integer i
      integer ibck
      integer icord
      integer icord1
      integer ilim
      integer itest
      integer k
      integer l
      real u(npde,npts)
      real uavg(npde)
      real udot(npde,npts)
      real ux(npde)
      real x(1)
      real xavgl
      real xavgr
c
c  common block mesh contains the user specified spatial
c  grid points.
c
c  common block coord contains 0,1, or 2 depending on whether
c  the problem is in cartesian, cylindrical, or spherical
c  coordinates, respectively.
c
      common /mesh/ x
      common /coord/ icord

      icord1 = icord + 1
c
c  Update the left boundary values.
c
      call bndry ( t, x(1), u, alpha, beta, gamma, npde )

      itest = 0
      dxi = 1.0E+00 / ( x(2) - x(1) )

      do k = 1, npde
        if ( beta(k) .eq. 0.0 ) then
          u(k,1) = gamma(k) / alpha(k)
          itest = itest + 1
        end if
      end do

      if ( itest .ne. npde ) then

        if ( itest .ne. 0 ) then
          call bndry ( t, x(1), u, alpha, beta, gamma, npde )
        end if
c
c  Evaluate D coefficients at the left boundary.
c
        call d ( t, x(1), u, dval, npde )
c
c  Form approximations to dU/dX at the left boundary.
c
        do k = 1, npde
          if (beta(k).eq.0.0) then
            ux(k) = dxi*(u(k,2)-u(k,1))
          else
            ux(k) = (gamma(k)-alpha(k)*u(k,1)) / beta(k)
          end if
        end do

      end if
c
c  Evaluate U-average in the first interval.
c
      do k=1,npde
        uavg(k) = 0.5*(u(k,2)+u(k,1))
      end do
c
c  Evaluate the D coefficients in the first interval.
c
      xavgr = 0.5 * ( x(2) + x(1) )
      call d ( t, xavgr, uavg, dval(1,1,2), npde )
      dxil = 1.0
      dxir = dxi

      if ( icord .ne. 0 ) then
        dxil = x(1)**icord
        dxir = dxir*xavgr**icord
      end if
c
c  Evaluate dUxx at the left boundary.
c
   70 if (itest.eq.npde) go to 100
      dxic = float(icord1) / (xavgr**icord1-x(1)**icord1)
      do l=1,npde
        do k=1,npde
          dval(k,l,1) = dxic*(dval(k,l,2)*(u(l,2)-u(l,1))*
     *     dxir-dval(k,l,1)*ux(l)*dxil)
        end do
      end do
c
c  Evaluate righthand side of PDE's at the left boundary
c
      call f ( t, x(1), u, ux, dval, udot, npde )
c
c  Set Udot = 0 for known left boundary values.
c
  100 do k=1,npde
        if (beta(k).eq.0.0) then
          udot(k,1) = 0.0
        end if
      end do
c
c  Update the right boundary values.
c
      call bndry ( t, x(npts), u(1,npts), alpha, beta, gamma, npde )
      itest = 0

      do k=1,npde
        if (beta(k).eq.0.0) then
          u(k,npts) = gamma(k) / alpha(k)
          itest = itest + 1
        end if
      end do

      ibck = 1
      ifwd = 2
      ilim = npts - 1
c
c  Main loop to form ordinary differential equations at the grid points.
c
      do i=2,ilim

        k = ibck
        ibck = ifwd
        ifwd = k
        xavgl = xavgr
        xavgr = 0.5*(x(i+1)+x(i))
        dxi = 1.0 / (x(i+1)-x(i-1))
        dxil = dxir
        if ( icord .eq. 0 ) then
          dxir = 1.0 / (x(i+1)-x(i))
        else
          dxir = xavgr**icord / ( x(i+1) - x(i) )
        end if
        dxic = float(icord1) / (xavgr**icord1-xavgl**icord1)
c
c  Evaluate dU/dX and U-average at the I-th grid point and
c  interval respectively.
c
        do k=1,npde
          ux(k) = dxi*(u(k,i+1)-u(k,i-1))
          uavg(k) = 0.5 * (u(k,i+1)+u(k,i))
        end do
c
c  Evaluate the D coefficients in the I-th interval.
c
        call d ( t, xavgr, uavg, dval(1,1,ifwd), npde )
c
c  Evaluate dUxx at the I-th grid point.
c
        do l=1,npde
          do k=1,npde
            dval(k,l,ibck) = (dval(k,l,ifwd)*(u(l,i+1)-u(l,i))*
     *       dxir-dval(k,l,ibck)*(u(l,i)-u(l,i-1))*dxil)*dxic
          end do
        end do
c
c  Evaluate righthand side of PDE's at the I-th grid point.
c
        call f ( t, x(i), u(1,i), ux, dval(1,1,ibck), udot(1,i), npde )

      end do
c
c  Finish updating the right boundary if necessary.
c
      if (itest.eq.npde) go to 220

      if ( itest .ne. 0 ) then
        call bndry ( t, x(npts), u(1,npts), alpha, beta, gamma, npde )
      end if
c
c  Evaluate the D coefficients at the right boundary.
c
      call d ( t, x(npts), u(1,npts), dval(1,1,ibck), npde )
c
c  Form approximations to dU/dX at the right boundary.
c
      dxi = 1.0 / (x(npts)-x(ilim))

      do k=1,npde
        if ( beta(k).eq.0.0) then
          ux(k) = dxi*(u(k,npts)-u(k,ilim))
        else
          ux(k) = (gamma(k)-alpha(k)*u(k,npts)) / beta(k)
        end if
      end do

      dxil = dxir
      if ( icord .eq. 0 ) then
        dxir = 1.0
      else
        dxir = x(npts)**icord
      end if
      dxic = float(icord1) / (x(npts)**icord1-xavgr**icord1)
c
c  Evaluate dUxx at the right boundary.
c
      do l = 1, npde
        do k = 1, npde
          dval(k,l,ibck) = dxic*(dval(k,l,ibck)*ux(l)*dxir-dval(k,l,
     *     ifwd)*(u(l,npts)-u(l,ilim))*dxil)
        end do
      end do
c
c  Evaluate righthand side of PDE's at the right boundary.
c
      call f ( t, x(npts), u(1,npts), ux, dval(1,1,ibck),
     *  udot(1,npts), npde )
c
c  Set Udot = 0 for known right boundary values.
c
  220 continue

      do k = 1, npde
        if ( beta(k) .eq. 0.0 ) then
          udot(k,npts) = 0.0
        end if
      end do

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
