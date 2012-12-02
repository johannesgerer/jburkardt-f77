      subroutine idcldp ( ndp, xd, yd, ncp, ipc )

c*********************************************************************72
c
cc IDCLDP finds several data points closest to a given point.
c
c  Discussion:
c
c    This subroutine selects several data points that are closest
c    to each of the data points.
c
c    This subroutine arbitrarily sets a restriction that NCP must
c    not exceed 25.
c
c  Modified:
c
c    31 December 2006
c
c  Author:
c
c    Hiroshi Akima,
c    modifications by Albrecht Preusser
c
c  Reference:
c
c    Hiroshi Akima,
c    Algorithm 526:
c    A Method of Bivariate Interpolation and Smooth Surface Fitting
c    for Values Given at Irregularly Distributed Points,
c    ACM Transactions on Mathematical Software,
c    Volume 4, Number 2, June 1978, pages 160-164.
c
c    Albrecht Preusser,
c    Algorithm 626:
c    TRICP: A Contour Plot Program for Triangular Meshes,
c    ACM Transactions on Mathematical Software,
c    Volume 10, Number 4, December 1984, pages 473-475.
c
c  Parameters:
c
c    Input, integer NDP, the number of data points.
c
c    Input, real XD(NDP), YD(NDP), the coordinates of the data points.
c
c    Input, integer NCP, the number of data points closest to each data
c    points.
c
c    Output, integer IPC(NCP*NDP), the indices of NCP data points 
c    closest to each of the NDP data points.
c
      implicit none

      integer ncp
      integer ncpmx
      parameter ( ncpmx = 25 )
      integer ndp

      real dsq0(ncpmx)
      real dsqf
      real dsqi
      real dsqmn
      real dsqmx
      real dx12
      real dx13
      real dy12
      real dy13
      integer ip1
      integer ip2
      integer ip2mn
      integer ip3
      integer ip3mn
      integer ipc(ncp*ndp)
      integer ipc0(ncpmx)
      integer j1
      integer j2
      integer j3
      integer j4
      integer jmx
      integer nclpt
      integer ncp0
      integer ndp0
      real u1
      real u2
      real v1
      real v2
      real x1
      real xd(ndp)
      real y1
      real yd(ndp)
c
c  Statement function
c
      dsqf ( u1, v1, u2, v2 ) = ( u2 - u1 )**2 + ( v2 - v1 )**2

      ndp0 = ndp
      ncp0 = ncp

      if ( ndp0 .lt. 2 ) then
        go to 90
      end if

      if ( ncp0 .lt. 1 .or. ncp0 .gt. ncpmx .or. ncp0 .ge. ndp0 ) then
        go to 90
      end if

      do ip1 = 1, ndp0
c
c  Select NCP points.
c
        x1 = xd(ip1)
        y1 = yd(ip1)
        j1 = 0
        dsqmx = 0.0E+00

        do ip2 = 1, ndp0

          if ( ip2 .ne. ip1 ) then

            dsqi = dsqf ( x1, y1, xd(ip2), yd(ip2) )
            j1 = j1+1
            dsq0(j1) = dsqi
            ipc0(j1) = ip2

            if ( dsqmx .lt. dsqi ) then
              dsqmx = dsqi
              jmx = j1
            end if

            if ( j1 .ge. ncp0 ) then
              go to 23
            end if

          end if

        end do

   23   continue

        ip2mn = ip2 + 1

        do ip2 = ip2mn, ndp0

          if ( ip2 .ne. ip1 ) then

            dsqi = dsqf ( x1, y1, xd(ip2), yd(ip2) )

            if ( dsqi .lt. dsqmx ) then

              dsq0(jmx) = dsqi
              ipc0(jmx) = ip2
              dsqmx = 0.0E+00

              do j1 = 1, ncp0
                if ( dsqmx .lt. dsq0(j1) ) then
                  dsqmx = dsq0(j1)
                  jmx = j1
                end if
              end do

            end if

          end if

        end do
c
c  Check if all the NCP+1 points are collinear.
c
        ip2 = ipc0(1)
        dx12 = xd(ip2) - x1
        dy12 = yd(ip2) - y1
c
c  0.06698 corresponds to an angle of about 15. degrees( = SIN**2)
c
        do j3 = 2, ncp0

          ip3 = ipc0(j3)
          dx13 = xd(ip3) - x1
          dy13 = yd(ip3) - y1

          if ( ( dy13 * dx12 - dx13 * dy12)**2 / ( dsq0(1) * dsq0(j3) )
     &       .gt. 0.06698E+00 ) then
            go to 50
          end if

        end do
c
c  Search for the closest noncollinear point.
c
        nclpt = 0

        do ip3 = 1, ndp0

          if ( ip3 .eq. ip1 ) then
            go to 43
          end if

          do j4 = 1, ncp0
            if ( ip3 .eq. ipc0(j4) ) then
              go to 43
            end if
          end do

          dx13 = xd(ip3) - x1
          dy13 = yd(ip3) - y1
          dsqi = dsqf ( x1, y1, xd(ip3), yd(ip3) )

          if ( ( dy13 * dx12 - dx13 * dy12 )**2 / ( dsq0(1) * dsqi )
     &        .le. 0.06698E+00 ) then
            go to 43
          end if

          if ( nclpt .eq. 0 ) then
            go to 42
          end if

          if ( dsqi .ge. dsqmn ) then
            go to 43
          end if

   42     continue

          nclpt = 1
          dsqmn = dsqi
          ip3mn = ip3

   43     continue

        end do

        if ( nclpt .eq. 0 ) then
          go to 91
        end if

        dsqmx = dsqmn
        ipc0(jmx) = ip3mn
c
c  Replace the local array for the output array.
c
   50   continue

        j1 = ( ip1 - 1 ) * ncp0

        do j2 = 1, ncp0
          j1 = j1 + 1
          ipc(j1) = ipc0(j2)
        end do

      end do

      return
c
c  Error exit
c
   90 write ( * ,2090 )
      go to 92
   91 write ( *, '(a)' ) '  The data points are collinear.'
   92 write ( *,2092)  ndp0,ncp0
      ipc(1) = 0
      return

 2090 format(1x/41h ***   improper input parameter value(s).)
 2092 format(8h   ndp = ,i5,5x,5hncp = ,i5/
     &   35h error detected in routine   idcldp/)
      end
      subroutine idpdrv ( ndp, xd, yd, zd, ncp, ipc, pd )

c*********************************************************************72
c
cc IDPDRV estimates partial derivatives at the data points.
c
c  Discussion:
c
c    This subroutine estimates partial derivatives of the first and
c    second order at the data points.
c
c  Modified:
c
c    29 December 2006
c
c  Author:
c
c    Hiroshi Akima,
c    modifications by Albrecht Preusser
c
c  Reference:
c
c    Hiroshi Akima,
c    Algorithm 526:
c    A Method of Bivariate Interpolation and Smooth Surface Fitting
c    for Values Given at Irregularly Distributed Points,
c    ACM Transactions on Mathematical Software,
c    Volume 4, Number 2, June 1978, pages 160-164.
c
c    Albrecht Preusser,
c    Algorithm 626:
c    TRICP: A Contour Plot Program for Triangular Meshes,
c    ACM Transactions on Mathematical Software,
c    Volume 10, Number 4, December 1984, pages 473-475.
c
c  Parameters:
c
c    Input, integer NDP, the number of data points.
c
c    Input, real XD(NDP), YD(NDP), ZD(NDP), the coordinates of the 
c    data points.
c
c    Input, integer NCPP, the number of additional data points used 
c    for estimating partial derivatives at each data point,
c
c    Input, integer IPC(NCP*NDP), the indices of NCP data points closest to
c    each of the NDP data points.
c
c    Output, real PD(5*NDP), the estimated ZX, ZY, ZXX, ZXY, and ZYY values.
c
      implicit none

      integer ncp
      integer ndp

      real dnmx
      real dnmxx
      real dnmxy
      real dnmy
      real dnmyx
      real dnmyy
      real dnmz
      real dx1
      real dx2
      real dy1
      real dy2
      real dz1
      real dz2
      real dzx1
      real dzx2
      real dzy1
      real dzy2
      integer ic1
      integer ic2
      integer ic2mn
      integer ip0
      integer ipc(ncp*ndp)
      integer ipi
      integer jipc
      integer jipc0
      integer jpd
      integer jpd0
      integer ncp0
      integer ncpm1
      integer ndp0
      real nmx
      real nmxx
      real nmxy
      real nmy
      real nmyx
      real nmyy
      real nmz
      real pd(5*ndp)
      real x0
      real xd(ndp)
      real y0
      real yd(ndp)
      real z0
      real zd(ndp)
      real zx0
      real zy0

      ndp0 = ndp
      ncp0 = ncp
      ncpm1 = ncp0 - 1
c
c  Estimation of ZX and ZY.
c
      do ip0 = 1, ndp0

        x0 = xd(ip0)
        y0 = yd(ip0)
        z0 = zd(ip0)

        nmx = 0.0E+00
        nmy = 0.0E+00
        nmz = 0.0E+00

        jipc0 = ncp0 * ( ip0 - 1 )

        do ic1 = 1, ncpm1

          jipc = jipc0 + ic1
          ipi = ipc(jipc)

          dx1 = xd(ipi) - x0
          dy1 = yd(ipi) - y0
          dz1 = zd(ipi) - z0

          ic2mn = ic1 + 1

          do ic2 = ic2mn, ncp0

            jipc = jipc0 + ic2
            ipi = ipc(jipc)
            dx2 = xd(ipi) - x0
            dy2 = yd(ipi) - y0
            dnmz = dx1 * dy2 - dy1 * dx2

            if ( dnmz .ne. 0.0E+00 ) then

              dz2 = zd(ipi) - z0
              dnmx = dy1 * dz2 - dz1 * dy2
              dnmy = dz1 * dx2 - dx1 * dz2

              if ( dnmz .lt. 0.0E+00 ) then
                dnmx = -dnmx
                dnmy = -dnmy
                dnmz = -dnmz
              end if

              nmx = nmx + dnmx
              nmy = nmy + dnmy
              nmz = nmz + dnmz

            end if

          end do

        end do

        jpd0 = 5 * ip0
        pd(jpd0-4) = -nmx / nmz
        pd(jpd0-3) = -nmy / nmz

      end do
c
c  Estimation of ZXX, ZXY, and ZYY.
c
      do ip0 = 1, ndp0

        jpd0 = jpd0 + 5
        x0 = xd(ip0)
        jpd0 = 5 * ip0
        y0 = yd(ip0)
        zx0 = pd(jpd0-4)
        zy0 = pd(jpd0-3)
        nmxx = 0.0E+00
        nmxy = 0.0E+00
        nmyx = 0.0E+00
        nmyy = 0.0E+00
        nmz = 0.0E+00
        jipc0 = ncp0 * ( ip0 - 1 )

        do ic1 = 1, ncpm1

          jipc = jipc0 + ic1
          ipi = ipc(jipc)
          dx1 = xd(ipi) - x0
          dy1 = yd(ipi) - y0
          jpd = 5 * ipi
          dzx1 = pd(jpd-4) - zx0
          dzy1 = pd(jpd-3) - zy0
          ic2mn = ic1+1
 
          do ic2 = ic2mn, ncp0

            jipc = jipc0 + ic2
            ipi = ipc(jipc)
            dx2 = xd(ipi) - x0
            dy2 = yd(ipi) - y0
            dnmz = dx1 * dy2 - dy1 * dx2

            if ( dnmz .ne. 0.0E+00 ) then

              jpd = 5 * ipi
              dzx2 = pd(jpd-4) - zx0
              dzy2 = pd(jpd-3) - zy0
              dnmxx = dy1 * dzx2 - dzx1 * dy2
              dnmxy = dzx1 * dx2 - dx1 * dzx2
              dnmyx = dy1 * dzy2 - dzy1 * dy2
              dnmyy = dzy1 * dx2 - dx1 * dzy2

              if ( dnmz .lt. 0.0E+00 ) then
                dnmxx = -dnmxx
                dnmxy = -dnmxy
                dnmyx = -dnmyx
                dnmyy = -dnmyy
                dnmz = -dnmz
              end if

              nmxx = nmxx + dnmxx
              nmxy = nmxy + dnmxy
              nmyx = nmyx + dnmyx
              nmyy = nmyy + dnmyy
              nmz = nmz + dnmz

            end if

          end do

        end do

        pd(jpd0-2) = - nmxx / nmz
        pd(jpd0-1) = - ( nmxy + nmyx ) / ( 2.0E+00 * nmz )
        pd(jpd0) = - nmyy / nmz

      end do

      return
      end
      subroutine idtang ( ndp, xd, yd, nt, ipt, nl, ipl, iwl, iwp, wk )

c*********************************************************************72
c
cc IDTANG determines a triangulation of points in the plane.
c
c  Discussion:
c
c    This subroutine performs triangulation.  It divides the X-Y
c    plane into a number of triangles according to given data
c    points in the plane, determines line segments that form the
c    border of data area, and determines the triangle numbers
c    corresponding to the border line segments.
c
c    At completion, point numbers of the vertexes of each triangle
c    are listed counter-clockwise.  Point numbers of the end points
c    of each border line segment are listed counter-clockwise,
c    listing order of the line segments being counter-clockwise.
c
c  Modified:
c
c    29 December 2006
c
c  Author:
c
c    Hiroshi Akima,
c    modifications by Albrecht Preusser
c
c  Reference:
c
c    Hiroshi Akima,
c    Algorithm 526:
c    A Method of Bivariate Interpolation and Smooth Surface Fitting
c    for Values Given at Irregularly Distributed Points,
c    ACM Transactions on Mathematical Software,
c    Volume 4, Number 2, June 1978, pages 160-164.
c
c    Albrecht Preusser,
c    Algorithm 626:
c    TRICP: A Contour Plot Program for Triangular Meshes,
c    ACM Transactions on Mathematical Software,
c    Volume 10, Number 4, December 1984, pages 473-475.
c
c  Parameters:
c
c    Input, integer NDP, the number of data points.
c
c    Input, real XD(NDP), YD(NDP), the coordinates of the data points.
c
c    Output, integer NT, the number of triangles.
c
c    Output, integer IPT(6*IDP-15), where the point numbers of the vertices
c    of triangle IT are stored in entries (3*IT-2), (3*IT-1) and 3*IT.
c
c    Output, integer NL, the number of border line segments,
c
c    Output, integer IPL(6*NDP), stores the point indices of the
c    IL-th border line segment and its respective triangle number 
c    as entries (3*IL-2), (3*IL-1), and (3*IL).
c
c    Workspace, integer IWL(18*NDP).
c
c    Workspace, integer IWP(NDP),
c
c    Workspace, real WK(NDP).
c
      implicit none

      integer ndp
      integer nt

      real ar
      real armx
      real armn
      real dsqf
      real dsqi
      real dsq12
      real dsqmn
      real dsqmx
      real dx
      real dx21
      real dx2l1
      real dxmn
      real dxmx
      real dy
      real dy21
      real dymn
      real dymx
      integer idxchg
      integer ilf
      integer ilft2
      integer ip
      integer ip1
      integer ip11
      integer ip12
      integer ip1p1
      integer ip2
      integer ip3
      integer ipl(600)
      integer ipl1
      integer ipl2
      integer iplj1
      integer iplj2
      integer ipll
      integer ipmn1
      integer ipmn2
      integer ipt(585)
      integer ipt1
      integer ipt2
      integer ipt3
      integer ipti
      integer ipti1
      integer ipti2
      integer irep
      integer it
      integer it1t3
      integer it2t3
      integer itf(2)
      integer its
      integer itt3
      integer itt3r
      integer iwl(1800)
      integer iwp(100)
      integer jlt3
      integer jp
      integer jp1
      integer jp2
      integer jp2t3
      integer jp3t3
      integer jpc
      integer jpmn
      integer jpmx
      integer jwl
      integer jwl1
      integer jwl1mn
      integer nl0
      integer ndp0
      integer ndpm1
      integer nl
      integer nlf
      integer nlfc
      integer nlft2
      integer nln
      integer nlnt3
      integer nlt3
      integer nrep
      integer nsh
      integer nsht3
      integer nt0
      integer ntf
      integer ntt3
      integer ntt3p3
      real ratio
      real side
      real u1
      real u2
      real u3
      real v1
      real v2
      real v3
      real wk(100)
      real x1
      real xd(ndp)
      real xdmp
      real y1
      real yd(ndp)
      real ydmp

      save nrep
      save ratio

      data nrep / 100 /
      data ratio / 1.0E-06 /
c
c  Statement functions
c
      dsqf(u1,v1,u2,v2) = ( u2 - u1 )**2 + ( v2 - v1 )**2

      side(u1,v1,u2,v2,u3,v3) = ( v3 - v1 ) * ( u2 - u1 ) 
     &                        - ( u3 - u1 ) * ( v2 - v1 )

      ndp0 = ndp
      ndpm1 = ndp0 - 1

      if ( ndp0 .lt. 4 ) then
        go to 90
      end if
c
c  Determine the closest pair of data points and their midpoint.
c
      dsqmn = dsqf ( xd(1), yd(1), xd(2), yd(2) )
      ipmn1 = 1
      ipmn2 = 2

      do ip1 = 1, ndpm1

        x1 = xd(ip1)
        y1 = yd(ip1)
        ip1p1 = ip1 + 1

        do ip2 = ip1p1, ndp0

          dsqi = dsqf ( x1, y1, xd(ip2), yd(ip2) )

          if ( dsqi .eq. 0.0E+00 ) then
            go to 91
          end if

          if ( dsqi .lt. dsqmn ) then
            dsqmn = dsqi
            ipmn1 = ip1
            ipmn2 = ip2
          end if

        end do

      end do

      dsq12 = dsqmn
      xdmp = ( xd(ipmn1) + xd(ipmn2) ) / 2.0E+00
      ydmp = ( yd(ipmn1) + yd(ipmn2) ) / 2.0E+00
c
c  Sort the other (NDP-2) data points in ascending order of
c  distance from the midpoint and store the sorted data point
c  numbers in the IWP array.
c
      jp1 = 2

      do ip1 = 1, ndp0
        if ( ip1 .ne. ipmn1 .and. ip1 .ne. ipmn2 ) then
          jp1 = jp1 + 1
          iwp(jp1) = ip1
          wk(jp1) = dsqf ( xdmp, ydmp, xd(ip1), yd(ip1) )
        end if
      end do

      do jp1 = 3, ndpm1

        dsqmn = wk(jp1)
        jpmn = jp1

        do jp2 = jp1, ndp0
          if ( wk(jp2) .lt. dsqmn ) then
            dsqmn = wk(jp2)
            jpmn = jp2
          end if
        end do

        its = iwp(jp1)
        iwp(jp1) = iwp(jpmn)
        iwp(jpmn) = its
        wk(jpmn) = wk(jp1)

      end do
c
c  If necessary, modify the ordering in such a way that the
c  first three data points are not collinear.
c
      ar = dsq12 * ratio
      x1 = xd(ipmn1)
      y1 = yd(ipmn1)
      dx21 = xd(ipmn2) - x1
      dy21 = yd(ipmn2) - y1

      do jp = 3, ndp0

        ip = iwp(jp)

        if ( abs ( ( yd(ip) - y1 ) * dx21 - ( xd(ip) - x1 ) * dy21 ) 
     &    .gt. ar ) then
          go to 37
        end if

      end do

      go to 92
 
37    continue

      if ( jp .ne. 3 ) then

        jpmx = jp
        jp = jpmx + 1

        do jpc = 4, jpmx
          jp = jp - 1
          iwp(jp) = iwp(jp-1)
        end do

        iwp(3) = ip

      end if
c
c  Form the first triangle.  Store point numbers of the vertexes of the 
c  triangle in the IPT array, and store point numbers of the border line 
c  segments and the triangle number in the IPL array.
c
      ip1 = ipmn1
      ip2 = ipmn2
      ip3 = iwp(3)

      if ( side ( xd(ip1), yd(ip1), xd(ip2), yd(ip2), xd(ip3), yd(ip3) )
     &     .lt. 0.0E+00 ) then
        ip1 = ipmn2
        ip2 = ipmn1
      end if

      nt0 = 1
      ntt3 = 3
      ipt(1) = ip1
      ipt(2) = ip2
      ipt(3) = ip3
      nl0 = 3
      nlt3 = 9
      ipl(1) = ip1
      ipl(2) = ip2
      ipl(3) = 1
      ipl(4) = ip2
      ipl(5) = ip3
      ipl(6) = 1
      ipl(7) = ip3
      ipl(8) = ip1
      ipl(9) = 1
c
c  Add the remaining (NDP-3) data points, one by one.
c
      do jp1 = 4, ndp0

        ip1 = iwp(jp1)
        x1 = xd(ip1)
        y1 = yd(ip1)
c
c  Determine the visible border line segments.
c
        ip2 = ipl(1)
        jpmn = 1
        dxmn = xd(ip2) - x1
        dymn = yd(ip2) - y1
        dsqmn = dxmn**2 + dymn**2
        armn = dsqmn * ratio
        jpmx = 1
        dxmx = dxmn
        dymx = dymn
        dsqmx = dsqmn
        armx = armn

        do jp2 = 2, nl0

          ip2 = ipl(3*jp2-2)
          dx = xd(ip2) - x1
          dy = yd(ip2) - y1
          ar = dy * dxmn - dx * dymn

          if ( ar .le. armn ) then

            dsqi = dx**2 + dy**2

            if ( ar .lt. (-armn) .or. dsqi .lt. dsqmn ) then
              jpmn = jp2
              dxmn = dx
              dymn = dy
              dsqmn = dsqi
              armn = dsqmn * ratio
            end if

          end if

          ar = dy * dxmx - dx * dymx

          if ( ar .ge. ( -armx ) ) then

            dsqi = dx**2 + dy**2

            if ( ar .gt. armx .or. dsqi .lt. dsqmx ) then

              jpmx = jp2
              dxmx = dx
              dymx = dy
              dsqmx = dsqi
              armx = dsqmx * ratio

            end if

          end if

        end do

        if ( jpmx .lt. jpmn ) then
          jpmx = jpmx + nl0
        end if

        nsh = jpmn - 1
c
c  Shift (rotate) the IPL array to have the invisible border
c  line segments contained in the first part of the IPL array.
c
        if ( 0 .lt. nsh ) then

          nsht3 = nsh * 3

          do jp2t3 = 3, nsht3, 3
            jp3t3 = jp2t3 + nlt3
            ipl(jp3t3-2) = ipl(jp2t3-2)
            ipl(jp3t3-1) = ipl(jp2t3-1)
            ipl(jp3t3) = ipl(jp2t3)
          end do

          do jp2t3 = 3, nlt3, 3
            jp3t3 = jp2t3 + nsht3
            ipl(jp2t3-2) = ipl(jp3t3-2)
            ipl(jp2t3-1) = ipl(jp3t3-1)
            ipl(jp2t3) = ipl(jp3t3)
          end do

          jpmx = jpmx - nsh

        end if
c
c  Add triangles to the IPT array, update border line
c  segments in the IPL array, and set flags for the border
c  line segments to be reexamined in the IWL array.
c
        jwl = 0

        do jp2 = jpmx, nl0

          jp2t3 = jp2 * 3
          ipl1 = ipl(jp2t3-2)
          ipl2 = ipl(jp2t3-1)
          it = ipl(jp2t3)
c
c  Add a triangle to the IPT array.
c
          nt0 = nt0 + 1
          ntt3 = ntt3 + 3
          ipt(ntt3-2) = ipl2
          ipt(ntt3-1) = ipl1
          ipt(ntt3) = ip1
c
c  Update border line segments in the IPL array.
c
          if ( jp2 .eq. jpmx ) then
            ipl(jp2t3-1) = ip1
            ipl(jp2t3) = nt0
          end if

          if ( jp2 .eq. nl0 ) then
            nln = jpmx + 1
            nlnt3 = nln * 3
            ipl(nlnt3-2) = ip1
            ipl(nlnt3-1) = ipl(1)
            ipl(nlnt3) = nt0
          end if
c
c  Determine the vertex that does not lie on the border
c  line segments.
c
          itt3 = it * 3
          ipti = ipt(itt3-2)

          if ( ipti .eq. ipl1 .or. ipti .eq. ipl2 ) then

            ipti = ipt(itt3-1)

            if ( ipti .eq. ipl1 .or. ipti .eq. ipl2 ) then
              ipti = ipt(itt3)
            end if

          end if
c
c  Check if the exchange is necessary.
c
c  63     continue

          if ( idxchg ( xd, yd, ip1, ipti, ipl1, ipl2 ) .ne. 0 ) 
     &    then
c
c  Modify the IPT array when necessary.
c
            ipt(itt3-2) = ipti
            ipt(itt3-1) = ipl1
            ipt(itt3) = ip1
            ipt(ntt3-1) = ipti

            if ( jp2 .eq. jpmx ) then
              ipl(jp2t3) = it
            end if

            if ( jp2 .eq. nl0 .and. ipl(3) .eq. it ) then
              ipl(3) = nt0
            end if
c
c  Set flags in the IWL array.
c
            jwl = jwl + 4
            iwl(jwl-3) = ipl1
            iwl(jwl-2) = ipti
            iwl(jwl-1) = ipti
            iwl(jwl) = ipl2

          end if

        end do

        nl0 = nln
        nlt3 = nlnt3
        nlf = jwl / 2

        if ( nlf .eq. 0 ) then
          go to 79
        end if
c
c  Improve triangulation.
c
        ntt3p3 = ntt3 + 3

        do irep = 1, nrep

          do ilf = 1, nlf

            ilft2 = ilf * 2
            ipl1 = iwl(ilft2-1)
            ipl2 = iwl(ilft2)
c
c  Locate in the IPT array two triangles on both sides of
c  the flagged line segment.
c
            ntf = 0

            do itt3r = 3, ntt3, 3

              itt3 = ntt3p3 - itt3r
              ipt1 = ipt(itt3-2)
              ipt2 = ipt(itt3-1)
              ipt3 = ipt(itt3)

              if ( ipl1.ne.ipt1 .and. ipl1 .ne. ipt2 .and.
     &           ipl1 .ne. ipt3 ) then
                go to 71
              end if

              if ( ipl2 .ne. ipt1 .and. ipl2 .ne. ipt2 .and.
     &           ipl2 .ne. ipt3 ) then
                go to 71
              end if

              ntf = ntf + 1
              itf(ntf) = itt3 / 3

              if ( ntf .eq. 2 ) then
                go to 72
              end if

   71         continue

            end do

            if ( ntf .lt. 2 ) then
              go to 76
            end if
c
c  Determine the vertexes of the triangles that do not lie
c  on the line segment.
c
   72       continue

            it1t3 = itf(1) * 3
            ipti1 = ipt(it1t3-2)

            if ( ipti1 .ne. ipl1 .and. ipti1 .ne. ipl2 ) then
              go to 73
            end if

            ipti1 = ipt(it1t3-1)

            if ( ipti1 .ne. ipl1 .and. ipti1 .ne. ipl2 ) then
              go to 73
            end if

            ipti1 = ipt(it1t3)

   73       continue

            it2t3 = itf(2) * 3
            ipti2 = ipt(it2t3-2)

            if ( ipti2 .ne. ipl1 .and. ipti2 .ne. ipl2 ) then
              go to 74
            end if

            ipti2 = ipt(it2t3-1)

            if ( ipti2 .ne. ipl1 .and. ipti2 .ne. ipl2 ) then
              go to 74
            end if

            ipti2 = ipt(it2t3)
c
c  Check if the exchange is necessary.
c
   74       continue

            if ( idxchg ( xd, yd, ipti1, ipti2, ipl1, ipl2 ) 
     &        .eq. 0 ) then
              go to 76
            end if
c
c  Modify the IPT array when necessary.
c
            ipt(it1t3-2) = ipti1
            ipt(it1t3-1) = ipti2
            ipt(it1t3) = ipl1
            ipt(it2t3-2) = ipti2
            ipt(it2t3-1) = ipti1
            ipt(it2t3) = ipl2
c
c  Set new flags.
c
            jwl = jwl + 8
            iwl(jwl-7) = ipl1
            iwl(jwl-6) = ipti1
            iwl(jwl-5) = ipti1
            iwl(jwl-4) = ipl2
            iwl(jwl-3) = ipl2
            iwl(jwl-2) = ipti2
            iwl(jwl-1) = ipti2
            iwl(jwl) = ipl1

            do jlt3 = 3, nlt3, 3

              iplj1 = ipl(jlt3-2)
              iplj2 = ipl(jlt3-1)

              if ( ( iplj1 .eq. ipl1 .and. iplj2 .eq. ipti2 ) .or.
     &            ( iplj2 .eq. ipl1 .and. iplj1 .eq. ipti2 ) ) then
                ipl(jlt3) = itf(1)
              end if

              if ( ( iplj1 .eq. ipl2 .and. iplj2 .eq. ipti1 ) .or.
     &             ( iplj2 .eq. ipl2 .and. iplj1 .eq. ipti1 ) ) then
                ipl(jlt3) = itf(2)
              end if

            end do

   76       continue

          end do

          nlfc = nlf
          nlf = jwl / 2

          if ( nlf .eq. nlfc ) then
            go to 79
          end if
c
c  Reset the IWL array for the next round.
c
          jwl = 0
          jwl1mn = ( nlfc + 1 ) * 2
          nlft2 = nlf * 2

          do jwl1 = jwl1mn, nlft2, 2
            jwl = jwl + 2
            iwl(jwl-1) = iwl(jwl1-1)
            iwl(jwl) = iwl(jwl1)
          end do

          nlf = jwl / 2

        end do

   79   continue

      end do
c
c  Rearrange the IPT array so that the vertexes of each triangle
c  are listed counter-clockwise.
c
      do itt3 = 3, ntt3, 3

        ip1 = ipt(itt3-2)
        ip2 = ipt(itt3-1)
        ip3 = ipt(itt3)

        if ( side ( xd(ip1), yd(ip1), xd(ip2), yd(ip2), xd(ip3), 
     &    yd(ip3) ) .lt. 0.0E+00 ) then
          ipt(itt3-2) = ip2
          ipt(itt3-1) = ip1
        end if

      end do

      nt = nt0
      nl = nl0
      return
c
c ERROR EXIT
c
   90 write ( *,2090)  ndp0
      go to 93
   91 write ( *,2091)  ndp0,ip1,ip2,x1,y1
      go to 93
   92 write ( *,2092)  ndp0
   93 write ( *,2093)
      nt = 0
      return

 2090 format(1x/23h ***   ndp less than 4./8h   ndp = ,i5)
 2091 format(1x/29h ***   identical data points./
     &   8h   ndp = ,i5,5x,5hip1 = ,i5,5x,5hip2 = ,i5,
     &   5x,4hxd = ,e12.4,5x,4hyd = ,e12.4)
 2092 format(1x/33h ***   all collinear data points./
     &   8h   ndp = ,i5)
 2093 format(35h error detected in routine   idtang/)
      end
      function idxchg ( x, y, i1, i2, i3, i4 )

c*********************************************************************72
c
cc IDXCHG determines whether to reform two triangles by swapping the diagonal.
c
c  Discussion:
c
c    This function determines whether or not the exchange of two
c    triangles is necessary on the basis of max-min-angle criterion
c    by C. L. Lawson.
c
c  Modified:
c
c    29 December 2006
c
c  Author:
c
c    Hiroshi Akima,
c    modifications by Albrecht Preusser
c
c  Reference:
c
c    Hiroshi Akima,
c    Algorithm 526:
c    A Method of Bivariate Interpolation and Smooth Surface Fitting
c    for Values Given at Irregularly Distributed Points,
c    ACM Transactions on Mathematical Software,
c    Volume 4, Number 2, June 1978, pages 160-164.
c
c    Albrecht Preusser,
c    Algorithm 626:
c    TRICP: A Contour Plot Program for Triangular Meshes,
c    ACM Transactions on Mathematical Software,
c    Volume 10, Number 4, December 1984, pages 473-475.
c
c  Parameters:
c
c    Input, real X(*), Y(*), the coordinates of the data points.
c
c    Input, integer I1, I2, I3, I4, the indices of four points that form a
c    quadrilateral, with P3 and P4 connected diagonally.
c
c    Output, integer IDXCHG, is 1 if an exchange is necessary,
c    and 0 otherwise.
c
      implicit none

      real a1sq
      real b1sq
      real c1sq
      real a2sq
      real b2sq
      real c3sq
      integer i1
      integer i2
      integer i3
      integer i4
      integer idx
      integer idxchg
      real s1sq
      real s2sq
      real s3sq
      real s4sq
      real u1
      real u2
      real u3
      real u4
      real x(*)
      real x1
      real x2
      real x3
      real x4
      real y(*)
      real y1
      real y2
      real y3
      real y4

      x1 = x(i1)
      y1 = y(i1)
      x2 = x(i2)
      y2 = y(i2)
      x3 = x(i3)
      y3 = y(i3)
      x4 = x(i4)
      y4 = y(i4)

      idx = 0

      u3 = ( y2 - y3 ) * ( x1 - x3 ) - ( x2 - x3 ) * ( y1 - y3 )
      u4 = ( y1 - y4 ) * ( x2 - x4 ) - ( x1 - x4 ) * ( y2 - y4 )

      if ( 0.0E+00 .lt. u3 * u4 ) then

        u1 = ( y3 - y1 ) * ( x4 - x1 ) - ( x3 - x1 ) * ( y4 - y1 )
        u2 = ( y4 - y2 ) * ( x3 - x2 ) - ( x4 - x2 ) * ( y3 - y2 )

        a1sq = ( x1 - x3 )**2 + ( y1 - y3 )**2
        b1sq = ( x4 - x1 )**2 + ( y4 - y1 )**2
        c1sq = ( x3 - x4 )**2 + ( y3 - y4 )**2

        a2sq = ( x2 - x4 )**2 + ( y2 - y4 )**2
        b2sq = ( x3 - x2 )**2 + ( y3 - y2 )**2
        c3sq = ( x2 - x1 )**2 + ( y2 - y1 )**2

        s1sq = u1 * u1 / ( c1sq * amax1 ( a1sq, b1sq ) )
        s2sq = u2 * u2 / ( c1sq * amax1 ( a2sq, b2sq ) )
        s3sq = u3 * u3 / ( c3sq * amax1 ( b2sq, a1sq ) )
        s4sq = u4 * u4 / ( c3sq * amax1 ( b1sq, a2sq ) )

        if ( amin1 ( s1sq, s2sq ) .lt. amin1 ( s3sq, s4sq ) ) then
          idx = 1
        end if

      end if

      idxchg = idx

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
c  Modified:
c
c    28 November 2006
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

      character*(8) date
      character*(10) time

      call date_and_time ( date, time )

      write ( *, '(a8,2x,a10)' ) date, time

      return
      end
      subroutine tricp ( xd, yd, zd, nd, c, nc, wk, iwk, mode )

c*********************************************************************72
c
cc TRICP is the user interface for the triangle contour plotting program.
c
c  Modified:
c
c    29 December 2006
c
c  Author:
c
c    Albrecht Preusser
c
c  Reference:
c
c    Albrecht Preusser,
c    Algorithm 626:
c    TRICP: A Contour Plot Program for Triangular Meshes,
c    ACM Transactions on Mathematical Software,
c    Volume 10, Number 4, December 1984, pages 473-475.
c
c  Parameters:
c
c    Input, real XD(ND), YD(ND), ZD(ND), the coordinates, and associated
c    function values, of the data points.
c
c    Input, integer ND, the number of data points.
c
c    Input, real C(NC), the contour levels.
c
c    Input, integer NC, the number of contour levels.
c
c    Workspace, real WK(5*ND), on output, contains the partial derivatives
c    (ZX(I),ZY(I),ZXX(I),ZXY(I),ZYY(I),I = 1,ND)
c
c     IWK        INTEGER WORK ARRAY. FOR MODE = 0, THE DIMENSION IS
c                                  31*ND  .
c                FOR OTHER MODES, SEE BELOW.
c
c                ON EXIT, IWK CONTAINS
c                A) IWK(1) = NT = number of triangles
c                B) IWK(2...3*NT+1) point numbers for the triangles.
c                   the first 3 numbers determine the vertices of the
c                   first triangle, the next 3 for the second and so on.
c                   The numbers coorespond to the indices of the
c                   XD,YD,ZD arrays. They are arranged counter-
c                   clockwise within a triangle.
c
c                C) IWK(3*NT+2  ...  3*NT+NCP*ND+1)
c                   A LIST OF NCP*ND POINT NUMBERS, REPRESENTING
c                   NCP POINTS FOR EACH DATA POINT THAT ARE
c                   THE CLOSEST TO THIS POINT. THE FIRST NCP
c                   NUMBERS ARE FOR THE FIRST DATA POINT, THE NEXT
c                   NCP FOR THE SECOND AND SO ON. THESE NUMBERS
c                   WERE USED FOR THE COMPUTATION OF THE PARTIAL
c                   DERIVATIVES.
c                   NCP IS AN INSTALLATION PARAMETER AND WILL BE SET
c                   TO 4.
c
c    Input, integer MODE, the usage mode.
c    * 0, normal mode, see above.
c    * 1, no triangulation requested.  IWK must contain the information about
c    the triangles as described under A) and B) on entry.
c    These locations of IWK will not be changed, and in addition, the information
c    under C) will be available on exit.
c    * 2, no triangulation, and no determination of the NCP closest points for
c    each data point.  On input, IWK must contain the information as described in
c    A), B) and C.  The contents of IWK will not be changed.
c    * 3, NO TRIANGULATION AND NO COMPUTATION OF THE PARTIAL DERIVATIVES.
c    IWK MUST CONTAIN THE INFORMATION A) AND B) AND WK MUST CONTAIN THE PARTIAL
c    DERIVATIVES (ZX(I),ZY(I),ZXX(I),ZXY(I),ZYY(I),I = 1,ND) ON ENTRY.
c    THE CONTENTS OF WK AND IWK WILL NOT BE CHANGED.  THIS MODE IS ESPECIALLY 
c    USEFUL WHEN TRICP IS CALLED AGAIN AFTER A PREVIOUS CALL. FOR INSTANCE,
c    ONLY THE CONTOUR LEVELS MAY HAVE CHANGED.  When designing a surface, it may 
c    be appropriate to change the XD,YD,ZD PARAMETERS AND THE PARTIAL
c    DERIVATIVES INTERACTIVELY AND TO CALL TRICP AGAIN USING THIS MODE.
c
c    FOR MODES 1, 2, AND 3 THE DIMENSION OF IWK CAN BE REDUCED TO
c    3*NT + NCP*ND + 1.
c
      implicit none

      integer nc
      integer nd

      real c(nc)
      integer iwk(*)
      integer mode
      integer n1
      integer n2
      integer n3
      integer n4
      integer ncp
      integer nl
      real wk(*)
      real xd(nd)
      real yd(nd)
      real zd(nd)
c
c  Set installation parameters.
c
      ncp = 4
c
c  Set starting addresses for IWK.
c
      n1 = 2
      n2 = n1 + 6 * nd - 15

      if ( 1 .le. mode ) then
        n2 = n1 + 3 * iwk(1)
      end if

      n3 = n2 + 6 * nd
      n4 = n3 + 18 * nd
c
c  Compute the triangulation if the user did not supply it.
c
      if ( mode .lt. 1 ) then

        call idtang ( nd, xd, yd, iwk, iwk(n1), nl, iwk(n2),
     &    iwk(n3), iwk(n4), wk )

        n2 = n1 + 3 * iwk(1)

      end if
c
c  Find NCP points closest to each data point.
c
      if ( mode .lt. 2 ) then
        call idcldp ( nd, xd, yd, ncp, iwk(n2) )
      end if
c
c  Compute the partial derivatives.
c
      if ( mode .lt. 3 ) then
        call idpdrv ( nd, xd, yd, zd, ncp, iwk(n2), wk )
      end if
c
c  Compute contours by successive solution of
c  quintic polynomial equations.
c
      if ( mode .lt. 4 ) then
        call trp00 ( xd, yd, zd, c, nc, iwk(n1), wk, iwk(1) )
      end if

      return
      end
      function tricpf ( t )

c*********************************************************************72
c
cc TRICPF evaluates the interpolating function or its derivatives.
c
c  Modified:
c
c    11 January 2007
c
c  Author:
c
c    Albrecht Preusser
c
c  Reference:
c
c    Albrecht Preusser,
c    Algorithm 626:
c    TRICP: A Contour Plot Program for Triangular Meshes,
c    ACM Transactions on Mathematical Software,
c    Volume 10, Number 4, December 1984, pages 473-475.
c
c  Parameters:
c
c    Input, real T, the evaluation point.
c
c    Output, real TRICPF, the value of the function at T.
c 
      implicit none

      real ap
      real bp
      real cl
      real cor
      real cp
      real dp
      real dx(3)
      real dx2(3)
      real dy(3)
      real dy2(3)
      real h0
      real h1
      real h2
      real h3
      real h4
      integer kk
      integer kse
      real p0(3)
      real p1(3)
      real p11
      real p12
      real p13
      real p14
      real p2(3)
      real p21
      real p22
      real p23
      real p3(3)
      real p31
      real p32
      real p4(3)
      real p41
      real p5(3)
      real q0(3)
      real q1(3)
      real q2(3)
      real q3(3)
      real q4(3)
      real r0(3)
      real r1(3)
      real r2(3)
      real r3(3)
      real s0(3)
      real s1(3)
      real s2(3)
      real sl(3)
      real sl2(3)
      real sir
      real t
      real t0(3)
      real t1(3)
      real tricpf
      real u
      real v
      real x(3)
      real xx4
      real xx4f
      real y(3)
      real yy4
      real yy4f
      real z(3)
      real zt(3,3)
      real ztt(3,3)
      real zuv(3)

      common /trpcof/ kk,kse,xx4f,yy4f,sir,cor,cl
c
c     KK      NUMBER OF FUNCTION TO BE EVALUATED
c     KK = 1    4TH DERIVATIVE ALONG SIDE KSE
c     KK = 2    3RD DERIVATIVE ALONG SIDE KSE
c     KK = 3    2ND DERIVATIVE ALONG SIDE KSE
c     KK = 4    1ST DERIVATIVE ALONG SIDE KSE
c     KK = 5    ORIGINAL POLYNOMIAL ALONG SIDE KSE
c     KK = 6    BIVARIATE POLYNOMIAL INSIDE TRIANGLE
c
      common /trpcop/  p0,p1,p2,p3,p4,p5,q0,q1,q2,q3,q4,
     &     r0,r1,r2,r3,s0,s1,s2,t0,t1
     &,     p11,p12,p13,p14,p21,p22,p23,p31,p32,p41

      common /trpcot/ x,y,z,dx,dy,dx2,dy2,sl,sl2,zt,ztt,zuv,
     &  ap,bp,cp,dp

      if ( kk .eq. 1 ) then

        tricpf = t0(kse) + t * t1(kse)

      else if ( kk .eq. 2 ) then

        tricpf = s0(kse) + t * ( s1(kse) + t * s2(kse) )

      else if ( kk .eq. 3 ) then

        tricpf = r0(kse) 
     &    + t * ( r1(kse) 
     &    + t * ( r2(kse) 
     &    + t *   r3(kse) ) )

      else if ( kk .eq. 4 ) then

        tricpf = q0(kse) 
     &    + t * ( q1(kse) 
     &    + t * ( q2(kse) 
     &    + t * ( q3(kse) 
     &    + t *   q4(kse))))

      else if ( kk .eq. 5 ) then

        tricpf = p0(kse) 
     &    + t * ( p1(kse)
     &    + t * ( p2(kse)
     &    + t * ( p3(kse)
     &    + t * ( p4(kse)
     &    + t *   p5(kse) )))) - cl

      else if ( kk .eq. 6 ) then

        xx4 = xx4f + cor * t
        yy4 = yy4f + sir * t

        u = ap * xx4 + bp * yy4
        v = cp * xx4 + dp * yy4

        h0 =       p0(1) 
     &     + v * ( p1(1) 
     &     + v * ( p2(1)
     &     + v * ( p3(1)
     &     + v * ( p4(1) 
     &     + v *   p5(1) ))))

        h1 =       p1(2) 
     &     + v * ( p11 
     &     + v * ( p12
     &     + v * ( p13
     &     + v *   p14 )))

        h2 =       p2(2)
     &     + v * ( p21 
     &     + v * ( p22
     &     + v *   p23 ))

        h3 =       p3(2) 
     &     + v * ( p31 
     &     + v *   p32 )

        h4 =       p4(2) 
     &     + v *   p41

        tricpf = h0 
     &   + u * ( h1 
     &   + u * ( h2
     &   + u * ( h3
     &   + u * ( h4
     &   + u *   p5(2) )))) - cl

      end if

      return
      end
      function tricpz ( ta, tb, f1, f2, er )

c*********************************************************************72
c
cc TRICPZ seeks a zero of a function between two points.
c
c  Discussion:
c
c    The method is a combination of the regula falsi
c    and the midpoint method.
c
c    IT IS A MODIFIED VERSION OF THE VIM- (CONTROL DATA
c    USER GROUP) ROUTINE WITH CATALOG IDENTIFICATION
c                C2BKYZERO
c    WRITTEN BY LOREN P. MEISSNER, 1965
c
c  Modified:
c
c    29 December 2006
c
c  Author:
c
c    Albrecht Preusser
c
c  Reference:
c
c    Albrecht Preusser,
c    Algorithm 626:
c    TRICP: A Contour Plot Program for Triangular Meshes,
c    ACM Transactions on Mathematical Software,
c    Volume 10, Number 4, December 1984, pages 473-475.
c
c  Parameters:
c
c    Input, real TA, TB, two points at which the function has
c    opposite sign.
c
c    Input, real F1, F2, the function value at TA and TB.
c
c    Input, real ER, the error tolerance.
c
c    Output, real TRICPZ, the location of the zero.
c
      implicit none

      real a
      real b
      real c
      real e
      real er
      real f1
      real f2
      real fa
      real fb
      real fc
      real fg
      real fs
      real fy
      real g
      real h
      real s
      real ta
      real tb
      real tricpf
      real tricpz
      real y

      a = ta
      b = tb
      fa = f1
      fb = f2
      c = a
      fc = fa
      s = c
      fs = fc

   10 continue

      h = 0.5E+00 * ( b + c )

      if ( abs ( h - b ) .le. er ) then
        tricpz = h
        return
      end if

      if ( abs ( fc ) .lt. abs ( fb )) then
        y = b
        fy = fb
        g = b
        fg = fb
        s = c
        fs = fc
      else
        y = s
        fy = fs
        g = c
        fg = fc
        s = b
        fs = fb
      end if

      if ( fy .eq. fs ) then

        b = h

      else

        e = ( s * fy - y * fs ) / ( fy - fs )

        if ( abs ( e - s ) .le. er ) then
          e = s + sign ( er, g - s )
        end if

        if ( ( e - h ) * ( s - e ) .lt. 0.0E+00 ) then
          b = h
        else
          b = e
        end if

      end if
c
c  Function call
c
      fb = tricpf ( b )

      if ( 0.0E+00 .le. fg * fb ) then
        c = s
        fc = fs
      else
        c = g
        fc = fg
      end if

      go to 10

      return
      end
      subroutine trp00 ( xd, yd, zd, cn, nc, ipt, pdd, nt )

c*********************************************************************72
c
cc TRP00 calls TRP001 and TRP002 for the computation of polynomial coefficients.
c
c  Discussion:
c
c    Plotting is done mainly by calling the Calcomp routine 
c
c      PLOT(X,Y,IPEN)
c
c    IPEN = 3 means move to (X,Y) with the pen up;
c    IPEN = 2 means move to (X,Y) with the pen down.
c
c  Modified:
c
c    29 December 2006
c
c  Author:
c
c    Albrecht Preusser
c
c  Reference:
c
c    Albrecht Preusser,
c    Algorithm 626:
c    TRICP: A Contour Plot Program for Triangular Meshes,
c    ACM Transactions on Mathematical Software,
c    Volume 10, Number 4, December 1984, pages 473-475.
c
c  Parameters:
c
c     XD,YD,ZD  COORDINATES OF DATA POINTS
c
c    Input, real CN(NC), the contour levels.
c
c    Input, integer NC, the number of contour levels.
c
c     IPT       POINT NUMBERS STORED AS 3*I-2, 3*I-1, 3*I TH
c               ELEMENT FOR TRIANGLE I
c
c     PDD       PARTIAL DERIVATIVES
c               (ZX,ZY,ZXX,ZXY,ZYY)
c
c    Input, integer NT, the number of triangles.
c
c  Local Parameters:
c
c    Local, real CMSCAL, variable for switching between CM and inches.
c
c    Local, real DSPMIN, the resolution of plotting device in use
c    (minimal distance of two points to be plotted)
c
c    Local, real DX, DY, the DIFFERENCES OF X and Y.
c
c    Local, real DX2,DY2    (DIFFERENCES OF X,Y)**2
c
c    Local, integer IN, NUMBER OF INTERVALS
c
c    Local, integer NPMAX, the maximum number of points per line for a triangle.
c
c    Local, real SL, the SIDE LENGTHS
c
c    Local, real SL2, the (SIDE LENGTHS)**2
c
c    Local, real TS1, U,V,W-VALUES AT ENDPOINTS OF INTERVALS
c
c    Local, real X,Y,Z, the COORDINATES AT THE THREE VERTICES OF A TRIANGLE.
c
c    Local, real Z1, Z-VALUES AT ENDPOINTS OF INTERVALS
c
c    ZT(IV,K)   FIRST PARTIAL DERIVATIVE WITH RESPECT TO VARIABLE IV
c                (U,V,W FOR IV = 1,2,3) AT POINT K
c
c    ZTT        SECOND PARTIAL DERIVATIVES
c
c    ZUV        MIXED PARTIAL DERIVATIVES
c
c    AP,BP,CP,DP  CONSTANTS DERIVED FROM DX,DY, NAMES DUE TO ALG. 526 CACM
c
c    ZMAX,ZMIN    MAXIMUM AND MINIMUM Z VALUE ALONG A SIDE
c    TS2          WORKING ARRAY FOR COMPUTING TS1
c    TZR          U,V,W VALUES FOR ZEROS FOUND ON SIDES
c    X0,Y0        X,Y COORDINATES FOR ZEROS FOUND
c    NZ           NUMBER OF ZEROS FOUND
c    TPER         PERMITTED POSITION ERROR FOR ZEROS ON SIDES
c    DT           DIFFERENCE IN U,V,W FOR ESTIMATING STARTING
c                  DIRECTION OF A CONTOUR LINE
c    THETAS       ANGLE OF DIRECTION FOR SIDES
c    SI,CO        COSINUS OF DIRECTION FOR SIDES
c
      integer nc
      integer nt

      real cmscal
      real cn(nc)
      real co(3)
      real dspmin
      real dt(3)
      integer in(3)
      integer ipt(*)
      integer it
      integer j
      integer ji
      integer npmax
      integer nz(3)
      real pdd(*)
      real si(3)
      real thetas(3)
      real tper(3)
      real ts1(6,3)
      real ts2(6)
      real tzr(5,3)
      real x0(5,3)
      real xd(*)
      real y0(5,3)
      real yd(*)
      real z1(6,3)
      real zd(*)
      real zmax(3)
      real zmin(3)

      common /trpcot/ x(3),y(3),z(3),dx(3),dy(3),dx2(3),dy2(3)
     &,               sl(3),sl2(3),zt(3,3),ztt(3,3),zuv(3)
     &,               ap,bp,cp,dp

      common /trpcof/ kk,kse,xx4f,yy4f,sir,cor,cl
c
c     /TRPCOF/ CONTAINS VARIABLES WHICH ARE PASSED TO FUNCTION
c              TRICPF AS PARAMETERS
c
c     KK          NUMBER OF FUNCTION TO BE EVALUATED BY TRICPF
c     KSE         ACTUAL SIDE NUMBER
c     XX4F,YY4F   COORDINATES FOR POINT P4F (PRELEMINARY POSITION
c                 OF POINT P4)
c     SIR,COR     COSINUS OF DIRECTION NORMAL TO CURVE DIRECTION
c     CL          ACTUAL SCALED CONTOUR LEVEL
c
      cmscal = 2.54E+00
      dspmin = 0.02E+00 / cmscal
      npmax = 500

      pi = 4.0E+00 * atan2 ( 1.0E+00, 1.0E+00 )

      do it = 1, nt
c
c  Load coordinates at vertices.
c
        ji = 3 * ( it - 1 )

        do j = 1, 3
          ji = ji + 1
          id = ipt(ji)
          x(j) = xd(id)
          y(j) = yd(id)
          z(j) = zd(id)
        end do

        zs = ( z(1) + z(2) + z(3) ) / 3.0E+00
c
c  Some basic geometry for the triangle.
c
        do j = 1, 3

          z(j) = z(j) - zs
          np1 = 3
          np2 = 2

          if ( j .eq. 3 ) then
            np1 = 1
          end if

          if ( j .eq. 2 ) then
            np2 = 1
          end if

          dx(j) = x(np2) - x(np1)
          dy(j) = y(np2) - y(np1)
          dx2(j) = dx(j) * dx(j)
          dy2(j) = dy(j) * dy(j)
          sl2(j) = dx2(j) + dy2(j)
          sl(j) = sqrt ( sl2(j) )
          co(j) = dx(j) / sl(j)
          si(j) = dy(j) / sl(j)

        end do

        co(1) = -co(1)
        si(1) = -si(1)
        slmax = amax1 ( sl(1), sl(2), sl(3) )
        di1 = abs ( dy(3) * co(1) - dx(3) * si(1) )
        di2 = abs ( dy(1) * co(2) - dx(1) * si(2) )
        di3 = abs ( dy(2) * co(3) - dx(2) * si(3) )
        hmin = amin1 ( di1, di2, di3 )
c
c  HMIN = shortest distance between a vertex and
c  its opposite side (minimum height).
c
c  Issue a warning message here, if HMIN is less than,
c  say 0.01/CMSCAL.
c
c  Define constants depending on HMIN.
c
        rmax = amin1 ( 0.02E+00 / cmscal, hmin * 0.01E+00 )
c
c  RMAX = distance normal to curve direction within which
c  a zero must be found.
c
        rmin = rmax * 0.2E+00
c
c  If zero has been found within a distance smaller
c  than RMIN, then stepsize DS is multiplied by 2.
c
        poserr = amin1 ( 1.0E-03 / cmscal, 
     &    1.0E-03 * hmin * hmin / slmax )
c
c  POSERR = permitted position error.
c
        dsmax = hmin * 0.2E+00
c
c  DSMAX = maximum step size.
c
        fstep = rmax * 10.0E+00
c
c  FSTEP = starting step size.
c
        eps = hmin * 0.01E+00
c
c  EPS = difference for estimating starting direction.
c
        rmax2 = rmax * 2.0E+00
c
c  RMAX2 = distance normal to curve direction. A zero must be
c  found within when crossing triangle border.
c
c  Compute coefficients for polynomials along sides.
c
        call trp001 ( it, ipt, pdd )

        cl = 0.0E+00
c
c  Loop over sides
c
        do jse = 1, 3

          kse = jse
          tper(kse) = poserr / sl(kse)
c
c  Compute endpoints of intervals
c
          ts1(1,kse) = 0.0E+00
          ts1(2,kse) = 1.0E+00
          i = 2

          do k = 1,4

            kk = k
            ts2(1) = 0.0E+00
            ii = 2
            tb = 0.0E+00
            f2 = tricpf ( tb )

            do j = 2, i

              ta = tb
              f1 = f2
              tb = ts1(j,kse)
              f2 = tricpf ( tb )

              if ( f1 * f2 .le. 0.0E+00 ) then

                if ( f1 .ne. 0.0E+00 .or. f2 .ne. 0.0E+00 ) then
                  ts2(ii) = tricpz ( ta, tb, f1, f2, tper(kse) )
                  ii = ii + 1
                end if

              end if

            end do

            ts2(ii) = 1.0E+00

            do j = 1, ii
              ts1(j,kse) = ts2(j)
            end do

            i = ii

          end do
c
c  IN(KSE) = number of intervals
c
          i = i - 1
          in(kse) = i
c
c  (E.G. IF IN(KSE) = 1, there is no point for which 1ST derivative. = 0)
c
c  Compute maxima and minima for each side.
c
          np1 = 3
          np2 = 2

          if ( kse .eq. 3 ) then
            np1 = 1
          end if

          if ( kse .eq. 2 ) then
            np2 = 1
          end if

          zmax(kse) = amax1 ( z(np1), z(np2) )
          zmin(kse) = amin1 ( z(np1), z(np2) )
          z1(1,kse) = z(np1)
          z1(i+1,kse) = z(np2)
          kk = 5

          do j = 2, i

            z1(j,kse) = tricpf ( ts1(j,kse) )

            if ( z1(j,kse) .gt. zmax(kse) ) then
              zmax(kse) = z1(j,kse)
            end if

            if ( z1(j,kse) .lt. zmin(kse) ) then
              zmin(kse) = z1(j,kse)
            end if

          end do

        end do

        nptri = 0
        ka = 1
        ke = 1
        kc = 0
c
c  Start loop over contour levels.
c
        do k = 1, nc

          cl = cn(k) - zs
c
c  Compute zeros (if any) on the three sides for level CL.
c
          kk = 5

          do jse = 1, 3

            nz(jse) = 0

            if ( cl .lt. zmin(jse) .or. cl .gt. zmax(jse) ) then
              go to 1200
            end if

            kse = jse
            jn = 0
            ni = in(kse)
c
c  Check intervals for zeros
c
            do j = 1, ni

              f1 = z1(j,kse) - cl
              f2 = z1(j+1,kse) - cl
              f1f2 = f1 * f2

              if ( f1f2 .gt. 0.0E+00 ) then
                go to 1100
              end if

              if ( f1f2 .lt. 0.0E+00 ) then
                go to 1090
              end if
c
c  Special situations.
c
              if ( ni .eq. 1 .and. 
     &          f1 .eq. 0.0E+00 .and. f2 .eq. 0.0E+00 ) then
                go to 1070
              end if

              if ( j .eq. 1  .and. f1 .eq. 0.0E+00 .or.
     &          j .eq. ni .and. f2 .eq. 0.0E+00 ) then
                go to 1080
              else
                go to 1090
              end if

 1070         continue
c
c  Contour line = side KSE.
c
              np1 = 3
              if ( kse .eq. 3 ) then
                np1 = 1
              end if

              np2 = 2
              if ( kse .eq. 2 ) then
                np2 = 1
              end if

              call plot ( x(np1), y(np1), 3 )
              call plot ( x(np2), y(np2), 2 )
              go to 1100

 1080         continue
c
c  Line passes through data point.
c  only one zero for this triangle, skip either side
c
              if ( kse .eq. 3 ) then
                go to 1100
              end if

              if ( f1 .eq. 0.0E+00 .and. kse .eq. 2 ) then
                go to 1100
              end if

              jn = jn + 1

              if ( f2 .eq. 0.0E+00 ) then
                tzr(jn,kse) = 1.0E+00 - tper(kse) * 0.5E+00
              end if

              if ( f1 .eq. 0.0E+00 ) then
                tzr(jn,kse) = tper(kse) * 0.5E+00
              end if

              go to 1100

 1090         continue

              jn = jn + 1
              tzr(jn,kse) = tricpz ( ts1(j,kse), ts1(j+1,kse), 
     &          f1, f2, tper(kse) )

 1100         continue

            end do

            nz(kse) = jn
c
c  NZ(KSE) = number of zeros for side KSE and contour level CL
c
 1200       continue

          end do

          if ( nz(1) + nz(2) + nz(3) .lt. 2 ) then
            go to 4000
          end if
c
c  If true then stop processing this contour level
c
c  Compute X0, Y0 for each zero (relative to X(3),Y(3))
c
          do jse = 1, 3

            jn = nz(jse)

            do jzr = 1, jn

              t = tzr(jzr,jse)

              if ( jse .ne. 3 ) then
                x0(jzr,jse) = dx(jse) * t
                y0(jzr,jse) = dy(jse) * t
              else
                x0(jzr,jse) = dx(2) + dx(3) * t
                y0(jzr,jse) = dy(2) + dy(3) * t
              end if

            end do

          end do

          kc = kc + 1

          if ( kc .eq. 1 ) then
c
c  Compute direction of sides.
c
            thetas(1) = atan2 (-dy(1),-dx(1) )
            thetas(2) = atan2 ( dy(2), dx(2) )
            thetas(3) = atan2 ( dy(3), dx(3) )
c
c  Compute differences for estimating start direction.
c
            dt(1) = -eps / sl(1)
            dt(2) =  eps / sl(2)
            dt(3) =  eps / sl(3)
c
c  Compute coefficients for polynomial inside triangle.
c
            call trp002

          end if
c
c  Follow contours
c
          rma = rmax
c
c  Loop over sides, KSE = side number
c
          do jse = 1, 3

            kse = ka
            jn = nz(kse)

            if ( jn .eq. 0 ) then
              go to 2010
            end if

            rma = sign ( rmax, -rma )
c
c  Loop over zeros
c
            do jzr = 1, jn

              if ( tzr(jzr,kse) .gt. 5.0E+00 ) then
                go to 2000
              end if
c
c  If true then this zero has already been processed
c
c  Estimate starting direction
c
c  Compute F1 on side using first derivative FD1
c
              kk = 4
              fd1 = tricpf ( tzr(jzr,kse) )
              f1 = fd1 * dt(kse)
c
c  Compute F2 normal to side.
c
              kk = 6
              xx4f = x0(jzr,kse)
              yy4f = y0(jzr,kse)
              sir = co(kse)
              cor = - si(kse)
              f2 = tricpf(eps)
c
c  Compute angle.
c
              if ( f2 .eq. 0.0E+00 ) then
                go to 1470
              end if

              thetsc = atan2 ( -f1, f2 )

              if ( thetsc .le. 0.0E+00 ) then 
                thetsc = thetsc + pi
              end if

              go to 1480

 1470         continue

              if ( fd1 .eq. 0.0 ) then
                go to 2000
              end if

              thetsc = pi * 0.5E+00

 1480         continue

              thetac = thetas(kse) + thetsc
c
c  Move pen to start.
c
              call plot ( x0(jzr,kse)+x(3), y0(jzr,kse)+y(3), 3 )
c
c  Initialize tracing
c
              ds = fstep
              dx12 = ds * cos ( thetac )
              dy12 = ds * sin ( thetac )
              xx3 = xx4f
              yy3 = yy4f
              xx2 = xx3 - dx12
              yy2 = yy3 - dy12
              xx1 = xx2 - dx12
              yy1 = yy2 - dy12
              ds23 = ds
              ds12 = ds
c
c  Points P1,P2,P3,P4 with coordinates XX1...XX4, YY1...YY4
c  are referred to as *queue*. DS12...DS34 are the distances
c  between points in the queue.
c
c  Every point normally must pass through the queue
c  before it is plotted.
c
c  First (oldest) point in the queue is P1, which is normally
c  (NPQ = 4) the first candidate to be plotted.
c  P4 is the next point to be computed.
c
              npq = 0
c
c  NPQ = number of points in the queue counted from the end
c  which can be plotted.
c  E.G. NPQ = 1, P4 is still to be plotted
c
              np = 0
c
c  NP = number of points computed for this contour line
c
              dsp = 0.0E+00
c
c  DSP = distance of point to be plotted to last point plotted
c
              nstop = 0
              nfound = 0
              nost = 0
c
c  Compute new point for contour line
c
c  Compute curve direction
c
 1500         continue

              ds13 = ds23 + ds12
              pl0 =  ds23 / ( ds12 * ds13 )
              pl1 = -ds13 / ( ds12 * ds23 )
              pl2 = ( ds13 + ds23 ) / ( ds13 * ds23 )
              dxds = pl0 * xx1 + pl1 * xx2 + pl2 * xx3
              dyds = pl0 * yy1 + pl1 * yy2 + pl2 * yy3
              sq = sqrt ( dxds * dxds + dyds * dyds )
              dxds = dxds / sq
              dyds = dyds / sq
              cor = -dyds
              sir =  dxds
c
c  Search for two points with opposite sign
c
              rma = sign ( rmax, rma )

 1550         continue

              xx4f = xx3 + dxds * ds
              yy4f = yy3 + dyds * ds
              f1 = tricpf ( 0.0E+00 )
              f2 = tricpf ( rma )

              if ( f1 * f2 .le. 0.0E+00 ) then
                go to 1600
              end if

              if ( abs(f2) .lt. abs(f1) ) then
                go to 1560
              end if

              rma = -rma
              f2 = tricpf ( rma )

              if ( f1 * f2 .le. 0.0E+00 ) then
                go to 1600
              end if
c
c  Divide stepsize in curve direction by 2.
c
 1560         continue

              ds = ds * 0.5E+00
              if ( ds .gt. poserr ) then
                go to 1550
              end if
c
c  Divide stepsize normal to curve by 2.
c
              nost = nost + 1
              ds = ds * 2.0E+00
              rma = rma * 0.5E+00

 1570         continue

              f2 = tricpf ( rma )

              if ( f1 * f2 .le. 0.0E+00 ) then
                go to 1600
              end if

              rma = -rma
              f2 = tricpf(  rma )

              if ( f1 * f2 .le. 0.0E+00 ) then
                go to 1600
              end if

              rma = rma * 0.5E+00

              if ( abs ( rma ) .gt. poserr ) then
                go to 1570
              end if

              nstop = 1
              go to 1900
c
c  Find zero for new point.
c
 1600         continue

              npq = npq + 1
              np = np + 1

              if ( np .gt. npmax ) then
                go to 1990
              end if

              r = tricpz ( 0.0E+00, rma, f1, f2, poserr )
              ds34 = sqrt ( ds * ds + r * r )
              xx4 = xx4f + cor * r
              yy4 = yy4f + sir * r
              u = ap * xx4  + bp * yy4
              v = cp * xx4  + dp * yy4
c
c  Check if point is outside the triangle.
c
              ke = 1
              if ( u .lt. 0.0E+00 ) then
                go to 1700
              end if

              ke = 2
              if ( v .lt. 0.0E+00 ) then
                go to 1700
              end if

              ke = 3
              if ( v .gt. 1.0E+00 - u ) then
                go to 1700
              end if
c
c  Point is inside the triangle.
c  Plot last point of queue and update queue.
c
              if ( npq .ne. 4 ) then
                go to 1610
              end if
c
c  Insert call for labeling routine here
c  Using XX1...XX4, YY1...YY4, DS12...DS34 as parameters
c  and reduce NPQ by the number of points which have been
c  used from the queue for labeling
c
              npq = npq - 1
              dsp = dsp + ds01

              if ( dspmin .le. dsp ) then
                call plot ( xx1+x(3), yy1+y(3), 2 )
                dsp = 0.0E+00
              end if

 1610         continue

              ds01 = ds12
              ds12 = ds23
              ds23 = ds34
              xx1 = xx2
              yy1 = yy2
              xx2 = xx3
              yy2 = yy3
              xx3 = xx4
              yy3 = yy4
c
c  Set new step size.
c
              if ( abs ( r ) .lt. rmin ) then
                ds = amin1 ( dsmax, ds * 2.0E+00 )
              end if

              go to 1500
c
c  Point is outside.
c
 1700         continue
c
c  No search if first point was computed and
c  start was at data point.
c
              if ( np .ne. 1 ) then
                go to 1710
              end if

              if ( tzr(jzr,kse) .gt. 1.0E+00 - tper(kse) .or.
     &          tzr(jzr,kse) .lt. tper(kse) ) then
                go to 2000
              end if

 1710         continue
c
c  Flag zero where line started.
c
              tzr(jzr,kse) = tzr(jzr,kse) + 10.0E+00
c
c  Search for corresponding zero on side KE.
c
              check = 5.0E+00
              dist = 99999.0E+00

              do n = 1, 2

                do jese = 1, 3

                  jj = nz(ke)

                  do j = 1, jj

                    if ( tzr(j,ke) .le. check ) then

                      di1 = abs ( ( y0(j,ke) - yy3 ) * dxds 
     &                          - ( x0(j,ke) - xx3 ) * dyds )

                      if ( di1 .lt. dist ) then
                        j1 = j
                        dist = di1
                      end if

                    end if

                  end do

                  if ( dist .lt. rmax2 ) then
                    go to 1800
                  end if

                  ke = mod ( ke, 3 ) + 1

                end do
c
c  For N = 2, disregard flags.
c
                check = 15.0E+00

              end do
c
c  No acceptable zero on all three sides.
c
              nfound = 1
              go to 1900
c
c  Flag zero found.
c
 1800         continue

              tzr(j1,ke) = tzr(j1,ke) + 10.0E+00
c
c  Clear queue and finish contour line.
c
 1900         continue

              if ( npq .gt. 3 ) then

                dsp = dsp + ds01

                if ( dspmin .le. dsp ) then
                  call plot ( xx1+x(3), yy1+y(3), 2 )
                  dsp = 0.0E+00
                end if

              end if

              if ( npq .gt. 2 ) then

                dsp = dsp + ds12

                if ( dspmin .le. dsp ) then
                  call plot ( xx2+x(3), yy2+y(3), 2 )
                  dsp = 0.0E+00
                end if

              end if

              if ( npq .gt. 1 ) then
  
                dsp = dsp + ds23

                if ( dspmin .le. dsp ) then
                  call plot ( xx3+x(3), yy3+y(3), 2 )
                end if

              end if

              if ( nfound + nstop  .eq. 0 ) then
                call plot ( x0(j1,ke)+x(3), y0(j1,ke)+y(3), 2 )
              end if

 1990         continue

              nptri = nptri + np

 2000         continue

            end do

 2010       continue

            ka = mod ( ka, 3 ) + 1

          end do

          ka = ke

 4000     continue

        end do

      end do

      return
      end
      subroutine trp001 ( it0, ipt, pdd )

c*********************************************************************72
c
cc TRP001 computes polynomial coefficients along sides.
c
c  Modified:
c
c    29 December 2006
c
c  Author:
c
c    Albrecht Preusser
c
c  Reference:
c
c    Albrecht Preusser,
c    Algorithm 626:
c    TRICP: A Contour Plot Program for Triangular Meshes,
c    ACM Transactions on Mathematical Software,
c    Volume 10, Number 4, December 1984, pages 473-475.
c
c  Parameters:
c
c    Input, integer IT0, the number of the triangle.
c
c     IPT     POINT NUMBERS STORED AS 3*I-2, 3*I-1, 3*I TH
c             ELEMENT FOR TRIANGLE I
c
c     PDD     PARTIAL DERIVATIVES AT DATA POINTS
c             (ZX,ZY,ZXX,ZXY,ZYY)
c
      real ad
      real adbc
      real ap
      real bc
      real bp
      real cp
      real dlt
      real dp
      real dxdy1
      real dxdy2
      real dxdy3
      integer idp
      integer ipt(*)
      integer it0
      integer jipt
      integer jpdd
      integer k
      real pd(5,3)
      real pdd(*)

      common /trpcop/ p0(3),p1(3),p2(3),p3(3),p4(3),p5(3)
     &,     q0(3),q1(3),q2(3),q3(3),q4(3)
     &,     r0(3),r1(3),r2(3),r3(3),s0(3),s1(3),s2(3),t0(3),t1(3)
     &,     p11,p12,p13,p14,p21,p22,p23,p31,p32,p41
c     P0...P5     COEFFICIENTS OF POLYNOMIALS ALONG THE SIDES
c     Q0...Q4     COEFFICIENTS OF 1ST DERIVATIVES ALONG THE SIDES
c     R0...R3     COEFFICIENTS OF 2ND DERIVATIVES ALONG THE SIDES
c     S0...S2     COEFFICIENTS OF 3RD DERIVATIVES ALONG THE SIDES
c     T0...T1     COEFFICIENTS OF 4TH DERIVATIVES ALONG THE SIDES
c     P11...P41   COEFFICIENTS FOR BIVARIATE POLYNOMIAL INSIDE TRIANGLE
c
      common /trpcot/ x(3),y(3),z(3),dx(3),dy(3),dx2(3),dy2(3)
     &,               sl(3),sl2(3),zt(3,3),ztt(3,3),zuv(3)
     &,               ap,bp,cp,dp
c
c  Load partial derivatives at the vertices.
c
      jipt = 3 * ( it0 - 1 )

      do k = 1, 3

        jipt = jipt + 1
        idp = ipt(jipt)
        jpdd = 5 * ( idp - 1 )

        do kpd = 1, 5
          jpdd = jpdd + 1
          pd(kpd,k) = pdd(jpdd)
        end do

      end do
c
c  Determine the coefficients for the coordinate system
c  transformation from the X-Y system to the U-V system.
c
      ad =  dx(2) * dy(1)
      bc =  dx(1) * dy(2)
      dlt = ad - bc
      ap =  dy(1) / dlt
      bp = -dx(1) / dlt
      cp = -dy(2) / dlt
      dp =  dx(2) / dlt
c
c  Convert the partial derivatives at the vertexes of the
c  triangle for the U-V coordinate system.
c
      ab =  dx(2) * dx(1)
      adbc = ad + bc
      cd =  dy(1) * dy(2)
      dxdy1 = 2.0E+00 * dx(1) * dy(1)
      dxdy2 = 2.0E+00 * dx(2) * dy(2)
      dxdy3 = 2.0E+00 * dx(3) * dy(3)

      do k = 1, 3
        zt(1,k) = dx(2) * pd(1,k) + dy(2) * pd(2,k)
        zt(2,k) = dx(1) * pd(1,k) + dy(1) * pd(2,k)
        ztt(1,k) = dx2(2) * pd(3,k) + dxdy2 * pd(4,k) + dy2(2) * pd(5,k)
        zuv(k) = ab * pd(3,k) + adbc * pd(4,k) + cd * pd(5,k)
        ztt(2,k) = dx2(1) * pd(3,k) + dxdy1 * pd(4,k) + dy2(1) * pd(5,k)
      end do

      do k = 1, 2
        zt(3,k) = dx(3) * pd(1,k)  + dy(3) * pd(2,k)
        ztt(3,k) = dx2(3) * pd(3,k) + dxdy3 * pd(4,k) + dy2(3) * pd(5,k)
      end do
c
c  Calculate the coefficients of the polynomials along
c  the three sides of the triangle.
c
      do jse = 1, 3

        np1 = 3
        if ( jse .eq. 3 ) then
          np1 = 1
        end if

        np2 = 2
        if ( jse .eq. 2 ) then
          np2 = 1
        end if

        iv = np2
        if ( jse .eq. 3 ) then
          iv = 3
        end if

        p0(jse) = z(np1)
        p1(jse) = zt(iv,np1)
        p2(jse) = 0.5E+00 * ztt(iv,np1)
        h1 = z(np2) - p0(jse) - p1(jse) - p2(jse)
        h2 = zt(iv,np2) - p1(jse) - ztt(iv,np1)
        h3 = ztt(iv,np2) - ztt(iv,np1)
        p3(jse) = 10.0E+00 * h1 - 4.0E+00 * h2 + 0.5E+00 * h3
        p4(jse) = -15.0E+00 * h1 + 7.0E+00 * h2 - h3
        p5(jse) =  6.0E+00 * h1 - 3.0E+00 * h2 + 0.5E+00 * h3
c
c  Calculate coefficients for derivatives along sides.
c
        q0(jse) =           p1(jse)
        q1(jse) = 2.0E+00 * p2(jse)
        q2(jse) = 3.0E+00 * p3(jse)
        q3(jse) = 4.0E+00 * p4(jse)
        q4(jse) = 5.0E+00 * p5(jse)
        r0(jse) =           q1(jse)
        r1(jse) = 2.0E+00 * q2(jse)
        r2(jse) = 3.0E+00 * q3(jse)
        r3(jse) = 4.0E+00 * q4(jse)
        s0(jse) =           r1(jse)
        s1(jse) = 2.0E+00 * r2(jse)
        s2(jse) = 3.0E+00 * r3(jse)
        t0(jse) =           s1(jse)
        t1(jse) = 2.0E+00 * s2(jse)

      end do

      return
      end
      subroutine trp002

c*********************************************************************72
c
cc TRP002 computes polynomial coefficients inside the triangle.
c
c  Discussion:
c
c    Computation of polynomial coefficients P11...P41 for
c    bivariate polynomial inside triangle.
c
c  Modified:
c
c    29 December 2006
c
c  Author:
c
c    Albrecht Preusser
c
c  Reference:
c
c    Albrecht Preusser,
c    Algorithm 626:
c    TRICP: A Contour Plot Program for Triangular Meshes,
c    ACM Transactions on Mathematical Software,
c    Volume 10, Number 4, December 1984, pages 473-475.
c
c  Parameters:
c
c    None
c
      implicit none

      real ap
      real bp
      real cp
      real dp
      real dx(3)
      real dx2(3)
      real dy(3)
      real dy2(3)
      real e1
      real g1
      real g2
      real g3
      real h1
      real h2
      real p0(3)
      real p1(3)
      real p11
      real p12
      real p13
      real p14
      real p2(3)
      real p21
      real p22
      real p23
      real p3(3)
      real p31
      real p32
      real p4(3)
      real p41
      real p5(3)
      real q0(3)
      real q1(3)
      real q2(3)
      real q3(3)
      real q4(3)
      real r0(3)
      real r1(3)
      real r2(3)
      real r3(3)
      real s0(3)
      real s1(3)
      real s2(3)
      real sl(3)
      real sl2(3)
      real t0(3)
      real t1(3)
      real x(3)
      real y(3)
      real z(3)
      real zt(3,3)
      real ztt(3,3)
      real zuv(3)

      common /trpcop/  p0,p1,p2,p3,p4,p5,q0,q1,q2,q3,q4,r0,r1,r2,
     &  r3,s0,s1,s2,t0,t1,p11,p12,p13,p14,p21,p22,p23,p31,p32,p41

      common /trpcot/ x,y,z,dx,dy,dx2,dy2,sl,sl2,zt,ztt,zuv,
     &  ap,bp,cp,dp

      p14 = p5(1) * ( 2.5E+00 * ( sl2(2) - sl2(3) ) / sl2(1) + 2.5E+00 )
      p41 = p5(2) * ( 2.5E+00 * ( sl2(1) -s l2(3) ) / sl2(2) + 2.5E+00 )
      p11 = zuv(3)

      h1 = zt(2,1) - p1(1) - p11 - p41
      h2 = zuv(1) - p11 - 4.0E+00 * p41

      p21 = 3.0E+00 * h1 - h2
      p31 = -2.0E+00 * h1 + h2

      h1 = zt(1,2) - p1(2) - p11 - p14
      h2 = zuv(2) - p11 - 4.0E+00 * p14

      p12 = 3.0E+00 * h1 - h2
      p13 = -2.0E+00 * h1 + h2

      h1 = 0.5E+00 * ztt(2,1) - p2(1) - p12
      h2 = 0.5E+00 * ztt(1,2) - p2(2) - p21
      e1 = 2.5E+00 * ( sl2(2) - sl2(1) ) / sl2(3) + 2.5E+00
      g1 = 3.0E+00 - e1
      g2 = -2.0E+00 + e1
      g3 = e1 * ( p5(1) - p5(2) + p41 - p14 )
     &  + p14 - 4.0E+00 * p41 + 5.0E+00 * p5(2)

      p22 = g1 * h1 + g2 * h2 + g3
      p32 = h1 - p22
      p23 = h2 - p22

      return
      end
