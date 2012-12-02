      program main

c*********************************************************************72
c
cc MAIN is the main program for SCPACK_PRB.
c
c  Modified:
c
c    07 December 2007
c
      write ( *, '(a)' ) ' '
      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCPACK_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the SCPACK library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCPACK_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 is the first test for SCPACK.
c
c  Discussion:
c
c    This test checks the computed forward map WSC and inverse map ZSC
c    against a known exact solution.  This is an example of a polygon
c    with a vertex at infinity.
c
c  Modified:
c
c    07 December 2007
c
c  Author:
c
c    Nick Trefethen
c
      implicit complex*16(c,w,z)
      implicit real*8(a-b,d-h,o-v,x-y)

      dimension z(20),w(20),betam(20),qwork(344)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Compare computed WSC and ZSC against'
      write ( *, '(a)' ) '  a known exact solution.'

      zero = (0.d0,0.d0)
      zi = (0.d0,1.d0)
c
c  set up problem:
c
      n = 4
      wc = dcmplx(0.d0,sqrt(2.d0))
      w(1) = zi
      w(2) = zero
      w(3) = (1.d20,1.d20)
      w(4) = zero
      betam(1) = 1.d0
      betam(2) = -.5d0
      betam(3) = -2.d0
      betam(4) = -.5d0
c
c  compute nodes and weights for parameter problem:
c
      nptsq =12
      call qinit(n,betam,nptsq,qwork)
c
c  solve parameter problem:
c  (initial guess must be given to avoid accidental exact solution)
c
      iprint = 0
      iguess = 1

      do k = 1,4
        z(k) = exp(dcmplx(0.d0,k-4.d0))
      end do

      tol = 1.d-14

      call scsolv(iprint,iguess,tol,errest,n,c,z,
     &  wc,w,betam,nptsq,qwork)
c
c  compare wsc(z) to exact values for various z:
c
      do i = 1,4
        zz = (.3d0,0.d0) * dcmplx(i-2.d0,.2d0*i+.5d0)
        ww = wsc(zz,0,zero,wc,0,n,c,z,betam,nptsq,qwork)
        ztmp = -zi * (zz-zi) / (zz+zi)
        wwex = zi * sqrt(-ztmp**2 + (1.d0,0.d0))
        err = abs(ww-wwex)
        write (6,201) zz,ww,wwex,err
      end do

      write (6,200)
c
c  compare zsc(w) to exact values for various w:
c
      do i = 1,6
        ww = dcmplx(i-2.d0,sqrt(i+1.d0))
        ier = 0
        zz = zsc(ww,0,zero,zero,wc,0,tol,ier,
     &    n,c,z,wc,w,betam,nptsq,qwork)
        wtmp = zi * sqrt((-1.d0,0.d0)-ww**2)
        zzex = -zi * (wtmp-zi) / (wtmp+zi)
        err = abs(zz-zzex)
        write (6,202) ww,zz,zzex,err
      end do

      return
  200 format (1x)
  201 format (' z,w,wex,err: ',3('(',f12.9,',',f12.9,') '),d11.4)
  202 format (' w,z,zex,err: ',3('(',f12.9,',',f12.9,') '),d11.4)
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 calls RESIST to check a resistance.
c
c  Discussion:
c
c    This routine computes the resistance of an L-shaped hexagon,
c    assuming that a voltage difference is applied between the
c    top and bottom edges.
c
c    The exact answer is sqrt(3) = 1.7320508.
c
c    The subroutine RESIST is used for the computation.  It
c    calls SCSOLV to compose two Schwarz-Christoffel maps to map
c    a polygon onto a rectangle and thereby determine the dimensions
c    of a rectangle that is electrically equivalent to the given
c    polygonal resistor.
c
c  Modified:
c
c    07 December 2007
c
c  Author:
c
c    Nick Trefethen
c
      implicit complex*16(c,w,z), real*8(a-b,d-h,o-v,x-y)
      dimension w(10),ibrk(4),qwork(300)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  Determine the resistance of an L-shaped'
      write ( *, '(a)' ) '  resistor.'

      n = 6
      w(1) = (0.d0,0.d0)
      w(2) = (2.d0,0.d0)
      w(3) = (2.d0,1.d0)
      w(4) = (1.d0,1.d0)
      w(5) = (1.d0,2.d0)
      w(6) = (0.d0,2.d0)
      wc = (.5d0,.5d0)
      ibrk(1) = 1
      ibrk(2) = 2
      ibrk(3) = 5
      ibrk(4) = 6
c
c  main loop: different accuracy specifications:
c
      do ndig = 1,11,5
        write (6,202) ndig
        r = resist(n,w,wc,ibrk,ndig,errest,qwork)
        write (6,201) r,errest
      end do

      return

  201 format ('   r,errest:',2d23.15)
  202 format (/' ndig =',i3,':')
      end
      function resist(n,w,wc,ibrk,ndig,errest,qwork)

c*********************************************************************72
c
cc RESIST returns the resistance of a polygonal resistor.
c
c  Discussion:
c
c    This function returns the resistance of a polygonally
c    shaped resistor.  computations are based on the schwarz-
c    christoffel transformation.  we normalize by assuming that
c    a square has resistance 1.
c
c  Modified:
c
c    07 December 2007
c
c  Author:
c
c    Nick Trefethen
c
c input parameters:
c
c   n      number of vertices of the polygon.  n must satisfy
c          4 .le. n .le. 20.
c
c   w      complex array of dimension at least n containing
c          the positions of the vertices, viewed as complex numbers.
c          the vertices must be listed in counterclockwise order
c          around the polygon.  it is a good idea to keep the
c          w(k) roughly on the scale of unity.
c
c   wc     a point in the interior of the polygon.  try to
c          pick wc as central as possible in the sense that
c          as few parts of the polygon as possible are shielded
c          from it by reentrant edges.
c
c   ibrk   array of dimension at least 4 containing indices of
c          the vertices which define the breaks
c          between constant-voltage and insulated portions
c          of the boundary.  the program will assume that
c          the voltage is applied between side w(ibrk(1))-w(ibrk(2))
c          and side w(ibrk(3))-w(ibrk(4)), with the other two
c          sides insulated.  the break vertices must be numbered
c          in counterclockwise order; thus the integers ibrk(i)
c          must increase with i, except that they may wrap
c          once around n.
c
c   ndig   input integer giving the desired number of digits
c          of accuracy in the result.  ndig must be at least 1.
c          it should be no more than a couple of digits less
c          than full-word precision.
c
c output parameter:
c
c   errest rough but conservative estimate of the size of
c          the error in the value returned.
c
c work space parameter:
c
c   qwork  real work array.  dimension at least (ndig+1)*(2n+3),
c          but no more than 460.
c
c
c some advice:
c
c   the program is not infallible.  if it yields strange
c   messages about not converging, it's probably best to
c   try a simpler geometry.
c
c   the amount of time required is roughly proportional to ndig
c   and also roughly proportional to n**3.
c   thus problems with many corners can
c   take quite a while -- several minutes of cpu time on
c   the ibm 370-168 for a problem with n = 20.
c
c lloyd n. trefethen
c department of mathematics
c massachusetts institute of technology
c november 1979, revised july 1983
c
      implicit complex*16(c,w,z), real*8(a-b,d-h,o-v,x-y)
c     map from disk to resistor:
      dimension z(20),w(1),ibrk(1),betam(20),qwork(1)
c     map from disk to rectangle with equal resistance:
      dimension z2(4),w2(4),betam2(4)

      zero = (0.d0,0.d0)
      pi = acos(-1.d0)
c
c compute angles and check input parameters:
c
      call angles(n,w,betam)

      if (ndig.lt.1) write (6,101)
      if (ndig.lt.1) stop
      if (n.lt.4) write (6,102)
      if (n.lt.4) stop
      do 1 k = 1,4
        if (ibrk(k).gt.n .or. ibrk(k).lt.1) goto 2
    1   continue
      goto 3
    2 write (6,103)
      stop
    3 continue
c
c set up first map:
c
      nptsq = ndig + 1
      tol = max(1.e-12, 10.d0**(-nptsq))
      call qinit(n,betam,nptsq,qwork)
      call scsolv(-2,0,tol,eest,n,c,z,wc,w,betam,nptsq,qwork)
c
c set up map to rectangle:
c (qwork is overwritten to save space)
c
      n2 = 4
      c2 = (1.d0,0.d0)
      do 9 k = 1,4
        betam2(k) = -.5d0
    9   z2(k) = z(ibrk(k))
      nptsq2 = ndig + 1
      call qinit(n2,betam2,nptsq2,qwork)
      do 12 k = 1,4
   12   w2(k)=wsc(z2(k),k,zero,zero,0,n2,c2,z2,betam2,nptsq2,qwork)
c
c compute length, width, resistance, and error estimate:
c
      xlen1 = abs(w2(2)-w2(3))
      xlen2 = abs(w2(4)-w2(1))
      errl = max(abs(xlen1-xlen2),2.d0*eest) / xlen1
      xwid1 = abs(w2(2)-w2(1))
      xwid2 = abs(w2(4)-w2(3))
      errw = max(abs(xwid1-xwid2),2.d0*eest) / xwid1
      resist = (xlen1+xlen2) / (xwid1+xwid2)
      errest = resist * (errl + errw)

      return

  101 format (' *** error in resist *** ndig should be at',
     &  ' least 1.'/)
  102 format (' *** error in resist *** n must be no',
     &  ' greater than 20'/' and no smaller than 4')
  103 format (' *** error in resist *** each ibrk(i) must',
     &  ' be in the',/,' range from 1 to n')
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 generates a rectilinear grid.
c
c  Discussion:
c
c    This routine computes the conformal map from a polygon (1) onto a
c    rectangle (2) and plots a corresponding grid of size nx by ny.
c    The corners of the rectangle are (0,i), (0,0), (h,0), (h,i).
c
c  Modified:
c
c    07 December 2007
c
c  Author:
c
c    Nick Trefethen
c
      implicit complex*16(c,w,z), real*8(a-b,d-h,o-v,x-y)

      dimension wgrid(0:41,0:41)
      dimension z1(12),w1(12),betam1(12),qwork1(270),ibrk(4)
      dimension z2(4),w2(4),betam2(4),qwork2(110)

      data zero /(0.0D+00, 0.0D+00)/

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  Generate a rectilinear grid.'
c
c  input data
c
c  Number of vertices:
c
      n1 = 6
c
c  Vertices:
c
      w1(1) = cmplx ( 0.0D+00,  1.0D+00 )
      w1(2) = cmplx ( 0.0D+00,  0.0D+00 )
      w1(3) = cmplx ( 2.0D+00,  0.0D+00 )
      w1(4) = cmplx ( 3.0D+00, -1.0D+00 )
      w1(5) = cmplx ( 3.5D+00, -0.5D+00 )
      w1(6) = cmplx ( 2.5D+00,  1.0D+00 )
c
c  Image of zero in the polygon.
c
      wc1 = cmplx ( 2.0D+00, 0.5D+00 )
c
c  Four distinguished vertices.
c
      ibrk(1) = 1
      ibrk(2) = 2
      ibrk(3) = 4
      ibrk(4) = 5
c
c  Number of streamlines, number of equipotential lines.
c
      nx = 7
      ny = 20

      nq1 = 3
      nq2 = nq1
      tol = 10.0D+00**(-(nq1+1))

      nxp = nx + 1
      nyp = ny + 1
c
c  compute parameters for map from disk to polygon
c
      call angles(n1,w1,betam1)
      call qinit(n1,betam1,nq1,qwork1)
      call scsolv(0,0,tol,eest,n1,c1,z1,wc1,w1,betam1,nq1,qwork1)
c
c  compute parameters for map from disk to rectangle
c
      n2 = 4
      c2 = ( 1.0D+00, 0.0D+00 )

      do k = 1,4
        betam2(k) = -0.5D+00
        z2(k) = z1(ibrk(k))
      end do

      call qinit(n2,betam2,nq2,qwork2)

      do k = 1,4
        w2(k) = wsc(z2(k),k,zero,zero,0,n2,c2,z2,betam2,nq2,qwork2)
      end do

      c2 = ( 0.0D+00, 1.0D+00 ) / (w2(1)-w2(2))
      wc2 = -c2*w2(2)

      do  k = 1,4
        w2(k) = c2*w2(k) + wc2
      end do

      write (6,102) (w2(k),k=1,4)
  102 format (' vertices of image rectangle: ',
     &  4(/'   (',e13.6,',',e13.6,')')/)
      h = abs(w2(3)-w2(2))
c
c  compute grid points
c
      do j = 0,nyp
        i1 = 0
        i2 = nxp
        if (j.eq.0.or.j.eq.nyp) i1 = 1
        if (j.eq.0.or.j.eq.nyp) i2 = nx
        do i = i1,i2
          ww2 = cmplx((i*h)/nxp,dble(j)/nyp)
          call nearw(ww2,zn2,wn2,kn2,n2,z2,wc2,w2,betam2)
          ier = 0
          iguess = 1
          if (i.eq.i1) iguess = 0
          zz = zsc(ww2,iguess,zz,zn2,wn2,kn2,1.0D-3,ier,n2,c2,
     &      z2,wc2,w2,betam2,nq2,qwork2)
          call nearz(zz,zn1,wn1,kn1,n1,z1,wc1,w1,betam1)
          wgrid(i,j) = wsc(zz,0,zn1,wn1,kn1,n1,c1,z1,betam1,nq1,qwork1)
          end do
          write (6,105) j,nyp
        end do

  105   format (' finished row',i3,'/',i2,' of grid points')
c
c  Draw plot
c
c   The " " is some kind of "end of line" marker.
c
      open ( unit = 10, file = 'scpack_prb_grid.txt' )

      do  k = 1,n1
        write (10,103) w1(k)
      end do

      write (10,103) w1(1)
      write (10, '(a)' ) ' '

  103 format (2f10.5)

      do j = 1,ny
        write (10,103) (wgrid(i,j),i=0,nx)
        write (10,103) wgrid(nxp,j)
        write (10,'(a)' ) ' '
      end do

      do i = 1,nx
        write (10,103) (wgrid(i,j),j=0,ny)
        write (10,103) wgrid(i,nyp)
        write (10,'(a)' ) ' '
      end do

      close ( unit = 10 )
c
c  Try again using XY and XYL format.
c
      open ( unit = 10, file = 'scpack_prb.xy' )

      do j = 1, ny
        do i = 1, nx
          write ( 10, '(2x,g14.6,2x,g14.6)' ) wgrid(i,j)
        end do
      end do

      close ( unit = 10 )

      open (  unit = 10, file = 'scpack_prb.xyl' )

      do j = 1, ny
        write ( 10, '(20i4)' ) ( i + ( j - 1 ) * nx, i = 1, nx )
      end do

      do i = 1, nx
        write ( 10, '(20i4)' ) ( i + ( j - 1 ) * nx, j = 1, ny )
      end do

      close ( unit = 10 )

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

