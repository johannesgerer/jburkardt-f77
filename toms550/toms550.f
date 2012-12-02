      subroutine r4_props ( kv, kfa, ax, ay, az, vx, vy, vz )

c*********************************************************************72
c
cc R4_PROPS calculates mass properties of a general 3D solid.
c
c  Discussion:
c
c    This routine calculates mass properties of a general three
c    dimensional solid using four point Gaussian quadrature
c    over triangles.
c
c    The routine calculates contribution of each surface to mass
c    property integrals.  Integration is performed in
c    subroutine SRFINT.  Mass property calculations are
c    then completed in PROPS.  Data is assumed to have
c    been previously set up in storage by a data entry routine.
c
c    Surfaces are processed one at a time and the mass property calculations 
c    are accumulated during the processing.  Surface processing involves 
c    starting with the initial face, and proceding through the secondary 
c    faces, if any.  Surface area vector components are summed for each face
c    of the surface.  The area vector, whose magnitude is the area of the 
c    surface, is then computed from the components and is accumulated in 
c    the total surface area of the solid.
c
c    In the face processing, faces are divided into all unique triangles 
c    having as vertices the the first one of the list and any other two 
c    adjacent vertices of the face polygon.  Four-point Gaussian quadrature 
c    is applied to each triangle, and the results are accumulated in the
c    total mass properties.  Thus each triangle of each face of each surface 
c    has a contribution in the final results.
c
c  Modified:
c
c    01 January 2007
c
c  Author:
c
c    Adrian Messner, G Taylor
c
c  Reference:
c
c    Adrian Messner, G Taylor,
c    Algorithm 550:
c    Solid Polyhedron Measures,
c    ACM Transactions on Mathematical Software,
c    Volume 6, Number 1, March 1980, pages 121-130.
c
c  Parameters:
c
c    Input, integer KV(*), the face list array.
c
c    Input, integer KFA(MAXS), contains pointers to the initial face
c    list in KV for each surface, ordered by surface number.
c    If L1, L2, and L3 are the locations in KV of face lists for 
c    surfaces 1, 2, 3, then the program would set
c      KFA(1) = L1
c      KFA(2) = L2
c      KFA(3) = L3
c    A pointer value of zero (KFA(I) = 0) indicates that surface I
c    was not utilized in defining the solid.  Hence gaps may
c    occur in the surface numbers.  Also, the order of data
c    entry of the face lists can be entirely arbitrary,
c    because of the afore-mentioned linkage mechanism.
c
c    Output, real AX(MAXS), AY(MAXS), AZ(MAXS), the components of
c    the surface area vector for each surface.
c
c    Input, real VX(MAXV), VY(MAXV), VZ(MAXV), the vertex coordinates.
c
c  Common Parameters:
c
c    Common, integer MAXS, the number of surfaces.
c
c    Common, integer MAXV, the number of vertices.
c
c    Common, real AREA, surface area of solid
c
c    Common, real WEIGHT, weight of solid
c
c    Common, real VOL, the volume of solid
c
c    Common, real DENS, the mass density of solid
c
c    Common, real CG(3), the center of gravity.
c
c    Common, real IO(3), the mass moments, with respect to the origin.
c
c    Common, real ICG(3), the mass moments, with respect to the center
c    of gravity.
c
c    Common, real PR(3), the product moments, with respect to the origin.
c
c    Common, real PRCG(3), the product moments, with respect to the
c    center of gravity.
c
c    Common, logical ERROR, error indicator
c
c    Common, integer NERR, error code
c
c  Local Parameters:
c
c    Local, integer LA, pointer to face list (initial or secondary) data.
c
c    Local, integer NS, the current surface being processed.
c
c solid polyhedron measures
c -------------------------
c
c solid
c   defined by a set of surfaces, where each surface is
c   composed of one or more coplanar closed polygonal
c   patches called faces.  disconnected faces and faces with
c   one or more holes in them may be used, if necessary,
c   to represent the solid.
c
c face
c   a face is represented by the set of vertex points of the
c   polygonal surface patch, ordered by starting with an
c   arbitrary vertex, and proceeding point by point around
c   the boundary until the first point is reached.
c
c face list
c   the list of the set of vertex numbers, ordered as above
c   for a face, is called a face list, which is the data
c   representation of a face.  a face list has the form
c        nsf,na,nb,nc,...,ni,...,na
c   where nsf is the surface number and na,nb,nc,ni are
c   vertex numbers.  note that the face list is a closed
c   list, starting and ending with the first point selected.
c   for any surface with one or more holes in
c   it, any interior boundary (hole) is represented by a
c   face list with vertices enumerated in a sequence
c   opposite to that given for an exterior boundary.
c
c initial, secondary face lists
c   when multiple face lists are needed to define a surface,
c   (ie, a surface with a hole or one with disconnected but
c   coplanar parts), the first list to appear in the data is
c   called the initial face list, and any other is called a
c   secondary face list.  these lists will have a common
c   surface number and will be linked together in storage by
c   the data entry routine.
c
c global data arrays
c
c linkage to secondary face lists
c   a face list in the array kv has the form
c        nsf,lsi,na,nb,nc,...,ni,nj,nk,...,na
c   where nsf ic the surface number, na,nb,... are vertex
c   numbers, and lsi is a link to a secondary face list.
c   lsi = 0 indicates this is the last face list for a surface.
c   note, it may be the only list.
c
      implicit none

      real area
      real ax(3)
      real ay(3)
      real az(3)
      real cg(3)
      real dens
      real dens2
      real dens3
      logical error
      real icg(3)
      real io(3)
      integer j
      integer kfa(*)
      integer kv(*)
      integer la
      integer maxs
      integer maxv
      integer nerr
      integer ns
      real pr(3)
      real prcg(3)
      real vol
      real vol2
      real vol3
      real vx(*)
      real vy(*)
      real vz(*)
      real weight

      common / r4_common / maxs, maxv, area, weight, vol, dens, 
     &  cg, io, icg ,pr, prcg, error, nerr
c
c  Clearing of accumulation parameters.
c
      area = 0.0E+00
      vol = 0.0E+00
      do j = 1, maxs
        ax(j) = 0.0E+00
        ay(j) = 0.0E+00
        az(j) = 0.0E+00
      end do

      do j = 1, 3
        cg(j) = 0.0E+00
        io(j) = 0.0E+00
        pr(j) = 0.0E+00
      end do
c
c  Each surface is processed.
c
      do ns = 1, maxs
c
c  LA is a pointer to the face list data
c
        la  =  kfa(ns)

        if ( la .ne. 0 ) then
c
c  Contribution of each face of current surface is calculated
c
3         continue

          if ( kv(la) .ne. ns ) then
            error = .true.
            nerr = 5
            return
          end if

          call r4_srfint ( kv(la+2), vx, vy, vz, ns, vol, cg, pr, io, 
     &      ax, ay, az )
c
c  Check is made for secondary face list.
c
          la = kv(la+1)

          if ( la .ne. 0 ) then
            go to 3
          end if
c
c  Area summation.
c
          area = area + sqrt ( ax(ns)**2 + ay(ns)**2 + az(ns)**2 )

        end if

      end do
c
c  Mass property calculations are completed.
c
      vol = vol / 3.0E+00
      weight = vol * dens

      vol2 = vol * 2.0E+00
      vol3 = vol / 3.0E+00

      dens2 = dens / 2.0E+00
      dens3 = dens / 3.0E+00

      cg(1) = cg(1) / vol2
      cg(2) = cg(2) / vol2
      cg(3) = cg(3) / vol2

      io(1) = io(1) * dens3
      io(2) = io(2) * dens3
      io(3) = io(3) * dens3

      icg(1) = io(1) - ( cg(2)**2 + cg(3)**2 ) * weight
      icg(2) = io(2) - ( cg(1)**2 + cg(3)**2 ) * weight
      icg(3) = io(3) - ( cg(1)**2 + cg(2)**2 ) * weight

      pr(1) = pr(1) * dens2
      pr(2) = pr(2) * dens2
      pr(3) = pr(3) * dens2

      prcg(1) = pr(1) - cg(1) * cg(2) * weight
      prcg(2) = pr(2) - cg(2) * cg(3) * weight
      prcg(3) = pr(3) - cg(3) * cg(1) * weight

      return
      end
      subroutine r4_srfint ( kv, vx, vy, vz, ns, vol, cg, pr, io, 
     &  ax, ay, az )

c*********************************************************************72
c
cc R4_SRFINT calculates surface integrals over a face of a solid.
c
c  Discussion:
c
c    The routine divides the face into triangles having as vertices the
c    first vertex of the polygon and two adjacent vertices
c    of the polygon.  
c
c    Numerical integration is performed using four-point Gaussian quadrature 
c    over the triangles.
c
c    The results are accumulated to get the total integrals over the surface.
c
c  Modified:
c
c    01 January 2007
c
c  Author:
c
c    Adrian Messner, G Taylor
c
c  Reference:
c
c    Adrian Messner, G Taylor,
c    Algorithm 550:
c    Solid Polyhedron Measures,
c    ACM Transactions on Mathematical Software,
c    Volume 6, Number 1, March 1980, pages 121-130.
c
c  Parameters:
c
c    Input, integer KV(*), the face list array.  Note that when we
c    call this routine, we are supplying the portion of KV associated
c    with surface NS.
c
c    Input, real VX(MAXV), VY(MAXV), VZ(MAXV), the vertex coordinates.
c
c    Input, integer NS, the current surface number.
c
c    Input/output, real VOL, the volume of solid
c
c    Input/output, real CG(3), center of gravity coordinate array
c
c    Input/output, real PR(3), coords of prod moments, with respect to origin
c
c    Input/output, real IO(3), coords of mass moments, with respect to origin
c
c    Output, real AX(MAXS), AY(MAXS), AZ(MAXS), the components of
c    the surface area vector for each surface.
c
c  Local Parameters:
c
c    Local, real ANX, ANY, ANZ, the area vector components for a 
c    single triangle
c
c    Local, real CX, CY, CZ, the 1st, 2nd or 3rd powers of quadrature 
c    coordinates, weighted and area normalized.
c
c    Local, real PX(4), PY(4), PZ(4), quadrature points.
c
c    Local, real W(4), quadrature weights.
c
c    Local, real X(3), Y(3), Z(3), the face triangle coordinates.
c
      implicit none

      integer quad_num
      parameter ( quad_num = 4 )

      real anx
      real any
      real anz
      real ax(3)
      real ay(3)
      real az(3)
      real cg(3)
      real cx
      real cy
      real cz
      real dx
      real dy
      real dz
      integer i
      real io(3)
      integer k
      integer kv(*)
      integer l
      integer ns
      real pr(3)
      real px(quad_num)
      real py(quad_num)
      real pz(quad_num)
      integer quad
      real w(quad_num)
      real vol
      real vx(*)
      real vy(*)
      real vz(*)
      real x(3)
      real y(3)
      real z(3)
c
c  w(i) are weighing factors for 4-point gaussian quadrature
c  over triangles.
c
      w(1) = -0.5625E+00
      w(2) = 0.5208333333333333E+00
      w(3) = 0.5208333333333333E+00
      w(4) = 0.5208333333333333E+00
c
c  The polygon is segmented into triangles, each having the
c  first vertex of the polygon as the first point of the triangle.
c
      k = kv(1)

      x(1) = vx(k)
      y(1) = vy(k)
      z(1) = vz(k)
c
c  Process the face list for this surface.
c
      l = 1

10    continue

        l = l + 1

        k = kv(l+1)
c
c  Testing for end of surface definition.
c
        if ( k .eq. kv(1) ) then
          return
        end if

        x(3) = vx(k)
        y(3) = vy(k)
        z(3) = vz(k)

        k = kv(l)

        x(2) = vx(k)
        y(2) = vy(k)
        z(2) = vz(k)
c
c  Data is now in terms of a single triangle.
c
        px(1) = ( x(1) + x(2) + x(3) ) / 3.0E+00
        py(1) = ( y(1) + y(2) + y(3) ) / 3.0E+00
        pz(1) = ( z(1) + z(2) + z(3) ) / 3.0E+00

        dx = 0.6E+00 * px(1)
        dy = 0.6E+00 * py(1)
        dz = 0.6E+00 * pz(1)

        do i = 1, 3
          k = i + 1
          px(k) = dx + 0.4E+00 * x(i)
          py(k) = dy + 0.4E+00 * y(i)
          pz(k) = dz + 0.4E+00 * z(i)
        end do
c
c  Area vectors.
c
        anx = 0.5E+00 * 
     &    ( z(1) * ( y(2) - y(3) ) 
     &    + z(2) * ( y(3) - y(1) )
     &    + z(3) * ( y(1) - y(2) ) )

        any = 0.5E+00 * 
     &    ( x(1) * ( z(2) - z(3) ) 
     &    + x(2) * ( z(3) - z(1) )
     &    + x(3) * ( z(1) - z(2) ) )

        anz = 0.5E+00 * 
     &    ( y(1) * ( x(2) - x(3) ) 
     &    + y(2) * ( x(3) - x(1) )
     &    + y(3) * ( x(1) - x(2) ) )

        ax(ns) = ax(ns) + anx
        ay(ns) = ay(ns) + any
        az(ns) = az(ns) + anz
c
c  Calculation of mass properties
c
        do quad = 1, quad_num
c
c  1st power of integration coordinates
c
          cx = anx * px(quad) * w(quad)
          cy = any * py(quad) * w(quad)
          cz = anz * pz(quad) * w(quad)

          vol = vol + ( cx + cy + cz )
c
c  2nd power of integration coordinates
c
          cx = cx * px(quad)
          cy = cy * py(quad)
          cz = cz * pz(quad)

          cg(1) = cg(1) + cx
          cg(2) = cg(2) + cy
          cg(3) = cg(3) + cz

          pr(1) = pr(1) + cx * py(quad)
          pr(2) = pr(2) + cy * pz(quad)
          pr(3) = pr(3) + cz * px(quad)
c
c  3rd power of integration coordinates
c
          cx = cx * px(quad)
          cy = cy * py(quad)
          cz = cz * pz(quad)

          io(1) = io(1) + cy + cz
          io(2) = io(2) + cz + cx
          io(3) = io(3) + cx + cy

        end do

      go to 10

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
