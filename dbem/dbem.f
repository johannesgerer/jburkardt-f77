        program main

c*********************************************************************72
c
cc MAIN is the main program for DBEM.
c
c  Discussion:
c
c    This program performs 2D elastic stress analysis in a homogeneous
c    isotropic medium using the direct boundary element method (BEM).
c
c    Capabilites include plane strain, plane stress and axisymmetry with
c    centrifugal and gravitational body forces.
c
c       this program will cater to :
c       maximum no. of nodes = 100
c       maximum no. of elements = 50  (quadratic elements)
c       maximum no. of surfaces = 10
c       maximum no. of interior points = 100
c       single region analysis
c       prescription of global boundary conditions only
c 
c        developed at :
c                     computational mechanics division
c                     department of civil engineering
c                     state university of new york at buffalo
c                     buffalo, new york 14260.
c
c        features added at :
c                     computer aided engineering center
c                     department of mechanical engineering
c                     worcester polytechnic institute
c                     worcester, massachusetts 01609.
c
c         ragevendra aithal :
c                     finite difference algorithm to calculate a,l b,l
c                     2d sensitivity analysis
c         jeffrey t borggaard :
c                     axisymmetric sensitivity analysis
c                     body force sensitivities
c                     axisymmetric potential analysis
c
c  Modified:
c
c    14 December 2007
c
c  Author:
c
c    Ragevendra Aithal,
c    Jeff Borggaard.
c
c  Reference:
c
c    Prasanta Banerjee, Roy Butterfield,
c    Boundary Element Methods in Engineering Science,
c    McGraw-Hill, 1994,
c    ISBN: 0070841209,
c    LC: TA347.B69.B36.
c
c  Parameters:
c
c    Common, integer IPLSTR, chooses the problem type.
c    0, plane strain,
c    1, plane stress,
c    2, axisymmetric elasticity,
c    3, axisymmetric potential.
c
      implicit real*8(a-h,o-z)

      parameter(mxsz = 405,mxsz2 = 600) 

      common /matrix/ bmat(mxsz,mxsz2),amat(mxsz,mxsz)

      common /matril/ bmatl(mxsz,mxsz2),amatl(mxsz,mxsz)  
      common /geomat/ cord(200,2),ntnode,ntelem,icon(100,3) 
      common /temp4/ t(600),u(600),tp(600),up(600) 

      common /geomal/ cordl(200,2),valuep
      common /inter/ nip,neli(200),xip(200),yip(200),nbp,nbn(200)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DBEM'
      write ( *, '(a)' ) '  Elasto-static stress analysis by'
      write ( *, '(a)' ) '  direct boundary element method.'
      write ( *, '(a)' ) ' '
 
      isens = 2
c
c  Perform an analysis without sensitivity.
c
      if ( isens .eq. 0 ) then

        call meshgn ( )
        call input ( nsize, ncolg )
        call gaussp
        call matcon
        call integra1 ( isens, bmat, amat, bmatl, amatl, nsize, ncolg )
        call partin ( 1 )
        call asmblr ( bmat, amat, nsize, ncolg, isens )
        call intstr
c
c  this algorithm is used to determine sensitivities using finite 
c  difference to determine amat,l and bmat,l
c
      else if ( isens .eq. 1 ) then
c
c  Read input data
c
        call meshgn ( )
        call input ( nsize, ncolg )
c
c  Set up gauss points and weights.
c
        call gaussp
c
c  set up material constants
c
        call matcon
c
c   direct integration  and assemble [f] and [g] matrices
c
        call integra1 ( isens, bmat, amat, bmatl, amatl, nsize, ncolg )

        call partin ( 1 )

        call store ( bmat, amat, bmatl, amatl, nsize, ncolg )
c
c   assemble system coefficient matrix and right hand vector;
c   finally, solve for the unknowns and print results
c
        call asmblr ( bmat, amat, nsize, ncolg, isens ) 
c
c  Calculate interior stresses and displacements if necessary
c
        call intstr 
c
c  Read the sensitivity analysis information.
c
        call inputl
c
c   calls the integration routine again
c
        call integra1 ( isens, bmat, amat, bmatl, amatl, nsize, ncolg )

        call partin ( 1 )
c
c   form the force vector
c
        call force ( bmat, amat, bmatl, amatl, nsize, ncolg )
c
c   Solve for the nodal sensitivity.
c
        call assm ( bmat, amat, nsize, ncolg )
c
c  Final nodal sensitivity recovery.
c
        call bunsnt

        close(15)
        close(12)
        close(13)
c
c   this algorithm is used for determining sensitivities using
c   implicit differentiation of the kernels
c
      else if ( isens .eq. 2 ) then

        call meshgn ( )
        call input ( nsize, ncolg )
        call gaussp
        call matcon
        call inputl1
        call integra1 ( isens, bmat, amat, bmatl, amatl, nsize, ncolg )
        call partin ( 2 )
        call asmblr ( bmat, amat, nsize, ncolg, isens )
        call intstr
        call force1 ( bmat, amat, bmatl, amatl, nsize, ncolg )
        call assm ( bmat, amat, nsize, ncolg )
        call bunsnt

      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DBEM:'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end       
      subroutine asmblr ( bmat, amat, nsize, ncolg, isens ) 

c*********************************************************************72
c
cc ASMBLR assembles and solves the system.
c
c  Modified:
c
c    14 December 2007
c
c  Author:
c
c    Ragevendra Aithal,
c    Jeff Borggaard.
c
      implicit real*8 (a-h,o-z)

      dimension bmat(nsize,ncolg),amat(nsize,nsize+1) 
      dimension res(3,4),n(3) 

      common / temp3 / bv(600),xv(600),bp(600)
      common / temp4 / t(600),u(600),tp(600),up(600) 
      common /temp6 / tpl(600),upl(600)
      common / geomat / cord(200,2),ntnode,ntelem,icon(100,3)
      common / boundc / bc(100,6),ibtype(100,6) 
      common / maprop / emod,pr,rho,shmod,iplstr,iprob
      common / trnode / in,rn(300,2) 
      common / intcon /cons1,cons2,cons3,cons4 
      common / fixdis / ifxd(200,2),fxdp(200,2) 
      common / segmen / nseg,isnod(10),isele(10) 
      common / parcns / c1,c2,c3,c4,c5,c6,c7,c8,c9,c10
      common / body / omega,accn  

      do i = 1,ntnode
        do k = 1,2
         ifxd(i,k) = 0
        end do
      end do

      uscal = shmod
      in = 0

      do i = 1,ncolg
        bp(i) = 0.0
      end do

      if (isens.ne.0) then
        write(13)((amat(i,j),j = 1,nsize),i = 1,nsize)
        write(13)((bmat(i,j),j = 1,ncolg),i = 1,nsize)
      end if

      if ( iprob .ne. 0 ) then

        do i = 1,nsize
          bp(i) = 0.0D+00
          do j = 1,ncolg
            bp(i) = bp(i)-bmat(i,j)*tp(j)
          end do
        end do

        do i = 1,nsize
          do j = 1,nsize
            bp(i) = bp(i)+amat(i,j)*up(j)
          end do
        end do
 
      end if

      do i = 1,ntelem 
        do l = 1,3 

          nf = icon(i,l)
          ng = (i-1)*3 + l
          nf1 = 2*nf - 1
          nf2 = nf1 + 1 
          ng1 = 2*ng - 1
          ng2 = ng1 + 1 
          l1 = 2*l - 1
          l2 = l1 + 1 
          id1 = ibtype(i,l1)
          id2 = ibtype(i,l2)
          xv(ng1) = bc(i,l1)
          xv(ng2) = bc(i,l2)

          if ( l .eq. 3 ) then

            ip = i + 1
            il = 1
            if ( i .eq. isele(1))  then
              ip = 1
            end if

            if ( 1 < nseg ) then

              do lip = 1,nseg
                if (i.eq.isele(lip)) then
                  ip = isele(lip-1) + 1
                  go to 16
                end if
              end do

            end if

16          continue

            idp1 = ibtype(ip,1)
            idp2 = ibtype(ip,2)

           else if ( l .ne. 2 ) then

             il = 3
             if (i.eq.1) then
               ip = isele(1)
             else
               ip = i - 1
             end if

             do lip = 2, nseg
               if (ip.eq.isele(lip-1)) then
                 ip = isele(lip)
                 go to 19
               end if
             end do

   19        continue

             idp1 = ibtype(ip,5)
             idp2 = ibtype(ip,6)

           end if

           ngp = 3*(ip-1) + il
           ngp1 = 2*ngp - 1
           ngp2 = 2*ngp

          if ( id1.ne.1 .and. ifxd(nf,1) .ne. 1) then

            do j = 1, nsize 
              bv(j) = amat(j,nf1) 
              amat(j,nf1) = -bmat(j,ng1)*uscal
              bmat(j,ng1) = -bv(j)
              if (l.ne.2.and.idp1.ne.1) then
                amat(j,nf1) = amat(j,nf1) - bmat(j,ngp1)*uscal
                bmat(j,ngp1) = -bv(j)
              end if
            end do
            ifxd(nf,1) = 1

          end if

          if (id2.ne.1) then
            if(ifxd(nf,2).ne.1) then

              do j = 1,nsize 
                bv(j) = amat(j,nf2) 
                amat(j,nf2) = -bmat(j,ng2)*uscal
                bmat(j,ng2) = -bv(j)
                if (l.ne.2.and.idp2.ne.1) then
                  amat(j,nf2) = amat(j,nf2) - bmat(j,ngp2)*uscal
                  bmat(j,ngp2) = -bv(j)
                end if
              end do

              ifxd(nf,2) = 1

            end if
          end if

        end do
      end do
c 
c  get the rhs vector
c
      do i = 1,nsize 
        bv(i) = bp(i)
        do j = 1,ncolg 
          bv(i) = bv(i) + bmat(i,j) * xv(j) 
        end do
      end do

      jn = nsize+1
      do i = 1,nsize
        amat(i,jn) = bv(i)
      end do

         do i = 1,nsize 
            xv(i) = 0.0D+00
         end do
c
c  call equation solver
c  save the amat before solving
c
        open(unit = 13, file = 'matrix.dat',form = 'unformatted')

        write(13)((amat(i,j),j = 1,nsize),i = 1,nsize)

        call solver ( amat, nsize ) 

        do i = 1,nsize
          xv(i) = amat(i,jn)
        end do
c
c  Print results.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Results of analysis:'
      write ( *, '(a)' ) '--------------------'

      write(*,900)
  900 format(//,
     1      ' ele. no  node no.   disp1 ' 
     2      ,9x,'disp2',12x,'trac1',9x,'trac2',/) 

      do i = 1, ntelem

        do j = 1,3
          j1 = 2*j - 1
          j2 = j1 + 1 
          id1 = ibtype(i,j1)
          id2 = ibtype(i,j2)
          val1 = bc(i,j1) 
          val2 = bc(i,j2) 
          nf = icon(i,j)
          n(j) = nf 
          nf1 = 2*nf - 1
          nf2 = nf1 + 1 

          if (id1.ne.2) then
            res(j,1) = xv(nf1)
            res(j,3) = val1 
          else
            res(j,1) = val1 
            res(j,3) = xv(nf1)*uscal
          end if

          if (id2.ne.2) then
            res(j,2) = xv(nf2)
            res(j,4) = val2 
          else
            res(j,2) = val2 
            res(j,4) = xv(nf2)*uscal
          end if

          if (ifxd(nf,1).eq.1) then
            res(j,1) = fxdp(nf,1) 
          end if

          if (ifxd(nf,2).eq.1) then
            res(j,2) = fxdp(nf,2) 
          end if

        end do
c
c  Store results in temp4 locations 
c
        do j = 1, 3
          uu2 = sqrt(res(j,1)*res(j,1)+res(j,2)*res(j,2))
          pp2 = sqrt(res(j,3)*res(j,3)+res(j,4)*res(j,4))
          write (6,910) i,n(j),(res(j,k),k = 1,4) 
          ng = 3*(i-1) + j
          nf1 = n(j)*2-1 
          nf2 = nf1+1
          ng1 = 2*ng-1 
          ng2 = ng1+1
          u(nf1) = res(j,1) 
          u(nf2) = res(j,2) 
          t(ng1) = res(j,3) 
          t(ng2) = res(j,4) 
        end do

        write ( *, '(a)' ) ' '

      end do

  910 format (1x,i6,2x,i7,3x,2(1p1e12.4,3x),2x,2(1p1e12.4,3x))
      return
      end
      subroutine assm ( gmat, fmat, nsize, ncolg ) 

c*********************************************************************72
c
cc ASSM assembles the sensitivity system.
c
c  Modified:
c
c    14 December 2007
c
c  Author:
c
c    Ragevendra Aithal,
c    Jeff Borggaard.
c
      implicit real*8 (a-h,o-z)

      dimension gmat(nsize,ncolg),fmat(nsize,nsize+1) 
      dimension res(3,4),n(3) 

      common / temp3 / bv(600),xv(600),bp(600)
      common / temp4 / t(600),u(600),tp(600),up(600) 
      common / temp5 / tl(600),ul(600)
      common / temp6 / tpl(600),upl(600)
      common / geomat / cord(200,2),ntnode,ntelem,icon(100,3)
      common / boundc / bc(100,6),ibtype(100,6) 
      common / maprop / emod,pr,rho,shmod,iplstr,iprob
      common / trnode / in,rn(300,2) 
      common / intcon /cons1,cons2,cons3,cons4 
      common / fixdis / ifxd(200,2),fxdp(200,2) 
      common / segmen / nseg,isnod(10),isele(10) 
      common / parcns / c1,c2,c3,c4,c5,c6,c7,c8,c9,c10
      common / body / omega,accn  

      uscal = shmod
      jn = nsize + 1 

      call solver ( fmat, nsize )

      do i = 1, nsize
        xv(i) = fmat(i,jn)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          Results of sensitivity analysis:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          ********************************'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &' ele.no nodeno.   disp1       disp2           trac1       trac2'
      write ( *, '(a)' ) ' '

      do i = 1, ntelem

        do j = 1, 3

          j1 = 2*j - 1
          j2 = j1 + 1 
          id1 = ibtype(i,j1)
          id2 = ibtype(i,j2)
          val1 = bc(i,j1) 
          val2 = bc(i,j2) 
          nf = icon(i,j)
          n(j) = nf 
          nf1 = 2*nf - 1
          nf2 = nf1 + 1 

          if (id1.ne.2) then
            res(j,1) = xv(nf1)
            res(j,3) = 0.0D+00
          else
            res(j,1) = 0.0D+00 
            res(j,3) = xv(nf1)*uscal
          end if

          if (id2.ne.2) then
            res(j,2) = xv(nf2)
            res(j,4) = 0.0D+00 
          else
            res(j,2) = 0.0D+00 
            res(j,4) = xv(nf2)*uscal
          end if

          if (ifxd(nf,1).eq.1) then
            res(j,1) = 0.0D+00 
          end if

          if (ifxd(nf,2).eq.1) then
            res(j,2)  = 0.0D+00
          end if

        end do

        do j = 1, 3
          write ( *, 
     &    '(1x,i3,2x,i4,3x,2(1p1e12.4,3x),2x,2(1p1e12.4,3x))' )
     &    i, n(j), ( res(j,k), k = 1, 4 )
          ng = 3 * (i-1) + j
          nf1 = n(j) * 2-1 
          nf2 = nf1 + 1
          ng1 = 2 * ng-1 
          ng2 = ng1 + 1
          ul(nf1) = res(j,1) 
          ul(nf2) = res(j,2) 
          tl(ng1) = res(j,3) 
          tl(ng2) = res(j,4) 
        end do

      end do

      return
      end
      subroutine autint ( ngp, xpt, ypt, xy, nsub, rmaxl )

c*********************************************************************72
c
cc AUTINT chooses the order of the integration rule.
c
c  Discussion:
c
c    See page 379 of the reference for a brief discussion and
c    additional references.
c
c  Modified:
c
c    29 December 2007
c
c  Author:
c
c    Ragevendra Aithal,
c    Jeff Borggaard.
c
c  Parameters:
c
c    Output, integer NGP, ?
c
c    Input, double precision XPT, ?
c
c    Input, double precision YPT, ?
c
c    Input, double precision XY(3,2), ?
c
c    Output, integer NSUB, ?
c
c    Input, double precision RMAXL, the maximum element length.
c
      implicit none

      double precision alg1
      double precision alg2
      double precision ertol
      integer ngp
      integer nord
      integer nsub
      double precision r1
      double precision r2
      double precision r3
      double precision rabc
      double precision rjc
      double precision rmaxl
      double precision rmin
      double precision rnod
      double precision xpt
      double precision xy(3,2)
      double precision ypt
c 
c  Error tolerance.
c
      nsub = 1
      rabc = sqrt ( ( xy(3,1) - xy(1,1) )**2 
     &            + ( xy(3,2) - xy(1,2) )**2 )
      ertol = 0.0000001D+00
c
c  Distance between field and nodes of segment.
c
      r1 = sqrt ( ( xy(1,1) - xpt )**2 + ( xy(1,2) - ypt )**2 )
      r2 = sqrt ( ( xy(2,1) - xpt )**2 + ( xy(2,2) - ypt )**2 )
      r3 = sqrt ( ( xy(3,1) - xpt )**2 + ( xy(3,2) - ypt )**2 )
c
c  Find the distance to the nearest node.
c
      rmin = min ( r1, min ( r2, r3 ) )

           if ( rmin .lt.           rmaxl ) then
        nsub = 4
      else if ( rmin .lt. 1.5D+00 * rmaxl ) then
        nsub = 3
      else if ( rmin .lt. 2.0D+00 * rmaxl ) then
        nsub = 2
      end if

      if ( rmin .lt. 1.0D-05 ) then
        ngp = 12
        return
      end if

      rjc = 0.25D+00 * rabc / rmin

      if ( abs ( rjc - 1.0D+00 ) .lt. 1.0D-05 ) then
        rjc = rjc + 0.00001D+00
      end if

      alg1 = dlog ( ertol / 8.0D+00 )
      alg2 = dlog ( rjc )
      rnod = 0.5D+00 * abs ( alg1 / alg2 - 1.0D+00 ) 
      nord = int ( rnod ) + 1

      if ( nord .gt. 8 ) then
        nord = 12 
      end if

      ngp = nord

      return
      end 
      subroutine axkrnl ( r, z, xy, bk, ak, bkl, akl, ngp, nsub, sgp,
     &  sl, ntype, ising, isens, ek, tk, xyl, rl, zl )

c*********************************************************************72
c
cc AXKRNL sets up integration kernels for the axisymmetric case.
c
c  Discussion:
c
c    This routine sets up kernels for integration for the 
c    axisymmetric case (axisymmetric loading assumed).
c
c    r and z are the coordinates of the load ring
c
c  Modified:
c
c    14 December 2007
c
c  Author:
c
c    Ragevendra Aithal,
c    Jeff Borggaard.
c
c  Reference:
c
c    DP Henry, DA Pape, PK Banerjee,
c    New axisymmetric BEM formulation for body forces using 
c    particular integrals, 
c    Journal of Engineering Mechanics,
c    Volume 113, Number 5, May 1987, pages 671-688
c
      implicit real*8 (a-h,o-z)

      dimension sl(4)

      dimension bk(2,6),ak(2,6),xy(3,2),gwpj(3),sgp(4,12) 
      dimension bkl(2,6),akl(2,6),xyl(3,2),gwpjl(3)
      common/gauswp/gp(12,12),gpw(12,12),gsp(8,8),gspw(8,8) 
      common/eliptc/a(5),b(5),c(4),ddd(4)
      common/maprop/ emod,pr,rho,shmod,iplstr,iprob
      common /axkern/ u(4,2),t(4,2),sk(4,4)

      shn1(s) =  2.0D+00*(s-0.5D+00)*(s-1.0D+00)
      shn2(s) =  -4.0D+00*s*(s-1.0D+00) 
      shn3(s) =  2.0D+00*s*(s-0.5D+00)
      dsh1(s) =  4.0D+00*s-3.0D+00
      dsh2(s) =  -8.0D+00*s+4.0D+00
      dsh3(s) =  4.0D+00*s-1.0D+00

      tol = 0.00001D+00
      pi = 3.141592653589793D+00
      xa = xy(1,1)
      xb = xy(2,1)
      xc = xy(3,1)
      ya = xy(1,2)
      yb = xy(2,2)
      yc = xy(3,2)

      pratio = pr
      g = emod/(2.0D+00*(1.0D+00+pratio))
c
c  Start Gaussian integration loop.
c
      do in = 1,nsub

        xl =  sl(in)

        do ig = 1,ngp 
c 
c  Find r,z coords of g-pt e,z respectivly
c 
          s = sgp(in,ig) 
          e = shn1(s)*xa +shn2(s)*xb +shn3(s)*xc
          f = shn1(s)*ya +shn2(s)*yb +shn3(s)*yc
          el = shn1(s)*xyl(1,1)+shn2(s)*xyl(2,1)+shn3(s)*xyl(3,1)
          fl = shn1(s)*xyl(1,2)+shn2(s)*xyl(2,2)+shn3(s)*xyl(3,2)
          zb = f - z
          zb2 = zb*zb
c
c  Find 1-d jacobian (scale factor ) 
c
          dxis =  dsh1(s)*xa +dsh2(s)*xb +dsh3(s)*xc
          dyis =  dsh1(s)*ya +dsh2(s)*yb +dsh3(s)*yc
          rjac = sqrt(dxis*dxis + dyis*dyis)
          eal = dsh1(s)*xyl(1,1)+dsh2(s)*xyl(2,1)+dsh3(s)*xyl(3,1)
          fal = dsh1(s)*xyl(1,2)+dsh2(s)*xyl(2,2)+dsh3(s)*xyl(3,2)
          rjacl = (dxis*eal+dyis*fal)/rjac
c
c  Outward normal at Gauss pt
c 
          xnur = dyis/rjac
          xnuz = -dxis/rjac
          xnurl = -xnur*rjacl/rjac+fal/rjac
          xnuzl = -xnuz*rjacl/rjac-eal/rjac
c
c  Integration kernels 
c
          gpwj = gpw(ngp,ig)*rjac*xl 
          gwpj(1) = shn1(s)*gpwj
          gwpj(2) = shn2(s)*gpwj
          gwpj(3) = shn3(s)*gpwj
          gpwjl = gpw(ngp,ig)*rjacl*xl
          gwpjl(1) = shn1(s)*gpwjl
          gwpjl(2) = shn2(s)*gpwjl
          gwpjl(3) = shn3(s)*gpwjl
    
          d2 = (e-r)**2+zb2
          d = dsqrt(d2)

          if(d.lt.tol)then
            write(*,*)'  this is an error diagnostic message :'
            write(*,1)r,z
 1          format(39hnumerical error likely from the point (,e11.4,1h,,e11.4,
     &      41h) being too close to an integration point)
          end if

          if ( r .gt. 0.0D+00 ) go to 2501
c
c  subroutine to calculate uij and tij kernels when load point
c  is on the axis of rotation by j t borggaard
c
          c1 = 1.0D+00/(8.0D+00*g*(1.0D+00-pratio))
          c2 = 3.0D+00 - 4.0D+00*pratio
          c3 = -0.25D+00/(1.0D+00-pratio)
          c4 = 1.0D+00 - 2.0D+00*pratio
          d3 = d*d2
          ezb = e*zb

          urr = 0.0D+00
          uzr = 0.0D+00
          urz = c1*e*ezb/d3
          uzz = c1*e*(c2/d + zb2/d3)

          x1 = (e*e*xnur+e*zb*xnuz)/d3
          x2 = (e*e*xnuz-e*zb*xnur)/d3
          trr = 0.0D+00
          tzr = 0.0D+00
          trz = c3*(3.0D+00*ezb/d2*x1+c4*x2)
          tzz = c3*(c4+3.0D+00*zb2/d2)*x1

          if (isens .eq. 2) then
            zbl = fl - zl
            dl = (e*el + zb*zbl)/d
            urrl = 0.0D+00
            uzrl = 0.0D+00
            urzl = c1*((2.0D+00*e*el*zb+e*e*zbl)/d3
     &        -3.0D+00*e*ezb*dl/d2/d2)
            uzzl = c1*(c2*(el/d-e*dl/d2)+(el*zb2+2.0D+00*ezb*zbl)/d3
     &        -3.0D+00*e*zb2*dl/d2/d2)

            x1l  = (2.0D+00*e*el*xnur+e*e*xnurl+(el*zb+e*zbl)*xnuz
     &        +ezb*xnuzl)/d3-3.0D+00*x1*dl/d
            x2l  = (2.0D+00*e*el*xnuz+e*e*xnuzl-(el*zb+e*zbl)*xnur
     &        -ezb*xnurl)/d3-3.0D+00*x2*dl/d
            trrl = 0.0D+00
            tzrl = 0.0D+00
            trzl = c3*(3.0D+00*(el*zb*x1+e*zbl*x1+e*zb*x1l)/d2
     &        -6.0D+00*ezb*x1*dl/d3
     &        + c4*x2l)
            tzzl = c3*(c4*x1l+3.0D+00*(zb2*x1l+2.0D+00*zb*zbl*x1)/d2
     &        -6.0D+00*zb2*x1*dl/d3)
          end if

          go to 450

2501  continue
c
c  Calculate kernel functions and their derivatives.
c
      q2 = (e+r)*(e+r) + zb2
      q = dsqrt(q2)
      q3 = q*q2
      d3 = d*d2
      rq = r*q
      c1 = 1.0D+00/(8.0D+00*pi*g*(1.0D+00-pratio))
      c2 = 3.0D+00 - 4.0D+00*pratio
      c3 = (1.0D+00 - pratio)/(1.0D+00 - 2.0D+00*pratio)
      c4 = pratio/(1.0D+00 - 2.0D+00*pratio)
      r2 = r*r
      e2 = e*e
      ezb = e*zb
      emr = e - r
      epr = e + r
      s2 = e2 + r2 + zb2
      s3 = r2 - e2 + zb2
      s4 = e2 - r2 + zb2
      x1 = (1.0D+00/e + (e+r)/q2)
      x2 = (e + r)/q2
      x4 = 8.0D+00*(1.0D+00-pratio)/rq
      x5 = (2.0D+00/d2 + 1.0D+00/q2)

      urr1 = (c2*s2 + zb2)/rq
      urr2 = (c2*q2 + zb2*s2/d2)/rq
      urz1 = zb/q
      urz2 = urz1*s3/d2
      uzr1 = ezb/rq
      uzr2 = ezb*s4/rq/d2
      uzz1 = c2*e/q
      uzz2 = e*zb2/d2/q
c
c  Calculation of elliptic integrals and their derivatives.
c
      vm = 4.0D+00*e*r/q2
      vm1 = 1.0D+00 - vm
      vlg = dlog(1.0D+00/vm1)
      ei1 = 1.38629436112D+00 + 0.5D+00*vlg
      ei1 = ei1 + (.096578619622D+00 + 0.12499929597D+00*vlg)*vm1
      ei1 = ei1 + (.031559431627D+00 + 0.070148757782D+00*vlg)*vm1**2
      ei1 = ei1 + (.023761224857D+00 + 0.044983875539D+00*vlg)*vm1**3
      ei1 = ei1 + (.02596288452D+00 + 0.018751660276D+00*vlg)*vm1**4
      ei1 = ei1 + (.0066398011146D+00 + 0.0018472341632D+00*vlg)*vm1**5
      ei2 = 1.0D+00
      ei2 = ei2 + (0.44315287472D+00 + 0.24999920273D+00*vlg)*vm1
      ei2 = ei2 + (0.057566998484D+00 + 0.093564907830D+00*vlg)*vm1**2
      ei2 = ei2 + (0.031761145524D+00 + 0.054260524448D+00*vlg)*vm1**3
      ei2 = ei2 + (0.030662347457D+00 + 0.021836021169D+00*vlg)*vm1**4
      ei2 = ei2 + (0.0076529606032D+00 + 0.0021247918284D+00*vlg)*vm1**5
      y1 = s3/2.0D+00/e
      y2 = -ei1/q2 + ei2/d2
      ei1r = y1*y2
      ei1z = -zb*y2
      ei2r = y1*(ei2-ei1)/q2
      ei2z = -zb*(ei2-ei1)/q2

      urr1r = 2.0D+00*c2*e/rq - urr1*x1
      urr2r = 2.0D+00*(c2*(epr)+zb2*(e/d2-(emr)*s2/d2/d2))/rq-urr2*x1
      urz1r = x1*urz1
      urz2r = 2.0D+00*ezb/q/d2 + (x1 + 2.0D+00*(emr)/d2)*urz2
      uzr1r = x2*uzr1
      uzr2r = 2.0D+00*e2*zb/rq/d2 - (x2 + 2.0D+00*(emr)/d2)*uzr2
      uzz1r = x2*uzz1
      uzz2r = (x2 + 2.0D+00*(emr)/d2)*uzz2
      urr1z = zb*(x4-urr1/q2)
      urr2z = zb*(2.0D+00*c2/rq + 2.0D+00*(s2-2.0D+00*zb2*e*r/d2)/d2/rq 
     &  - urr2/q2)
      urz1z = (q2 - zb2)/q3
      urz2z = (e2-r2-3.0D+00*zb2)/d2/q + zb*x5*urz2
      uzr1z = e/r*urz1z
      uzr2z = e*(e2-r2+3.0D+00*zb2)/rq/d2 - zb*x5*uzr2
      uzz1z = zb*uzz1/q2
      uzz2z = 2.0D+00*ezb/d2/q - zb*x5*uzz2

      urr = c1*(urr1*ei1 - urr2*ei2)
      urz = c1*(urz1*ei1 - urz2*ei2)
      uzr = -c1*(uzr1*ei1 - uzr2*ei2)
      uzz = 2.0D+00*c1*(uzz1*ei1 + uzz2*ei2)

      urrr = c1*(urr1*ei1r + urr1r*ei1 - urr2*ei2r - urr2r*ei2)
      urzr = c1*(urz1*ei1r - urz1r*ei1 - urz2*ei2r + urz2r*ei2)
      uzrr = -c1*(uzr1*ei1r - uzr1r*ei1 - uzr2*ei2r - uzr2r*ei2)
      uzzr = 2.*c1*(uzz1*ei1r - uzz1r*ei1 + uzz2*ei2r - uzz2r*ei2)
      urrz = c1*(urr1*ei1z + urr1z*ei1 - urr2*ei2z - urr2z*ei2)
      urzz = c1*(urz1*ei1z + urz1z*ei1 - urz2*ei2z + urz2z*ei2)
      uzrz = -c1*(uzr1*ei1z + uzr1z*ei1 - uzr2*ei2z - uzr2z*ei2)
      uzzz = 2.*c1*(uzz1*ei1z - uzz1z*ei1 + uzz2*ei2z + uzz2z*ei2)

      trr1 = c3*urrr + c4*(urr/e + uzrz)
      trr2 = 0.5D+00*(urrz + uzrr)
      tzr1 = c3*uzrz + c4*(urr/e + urrr)
      tzr2 = trr2
      trz1 = c3*urzr + c4*(urz/e + uzzz)
      trz2 = 0.5D+00*(urzz + uzzr)
      tzz1 = c3*uzzz + c4*(urz/e + urzr)
      tzz2 = trz2

      trr = 2.0D+00*g*(trr1*xnur + trr2*xnuz)
      tzr = 2.0D+00*g*(tzr1*xnuz + tzr2*xnur)
      trz = 2.0D+00*g*(trz1*xnur + trz2*xnuz)
      tzz = 2.0D+00*g*(tzz1*xnuz + tzz2*xnur)
c
c  skip to gaussian integration if isens ne 2
c
      if ( isens .eq. 2) then

      zbl = fl - zl
      zzl = zb*zbl
      ql = (e+r)*(el+rl)/q + zzl/q
      dl = (e-r)*(el-rl)/d + zzl/d
      s2l = e*el + r*rl + zzl
      rql = rl*q + r*ql
      rq2 = r2*q2
      s3l = r*rl - e*el + zzl
      s4l = e*el - r*rl + zzl
      x2l = (el+rl)/q2 - 2.0D+00*(epr)*ql/q3
      eprl = el + rl
      emrl = el - rl
      ezl = el*zb + e*zbl
      x6 = -el/e2 + x2l
      y3 = (el/e + rl/r)/2.0D+00 - ql/q

      ei1l = y2*q2*y3
      ei2l = y3*(ei2-ei1)

      y2l = -ei1l/q2 + 2.*ei1*ql/q3 + ei2l/d2 - 2.0D+00*ei2*dl/d3
      ei1rl = (s3l/e - s3*el/2.0D+00/e2)*y2 + y1*y2l
      ei1zl = -zbl*y2 - zb*y2l
      ei2rl = ((2.*s3l-(el*q+2.*e*ql)*s3/e/q)*(ei2-ei1)
     &        + s3*(ei2l-ei1l))/(2.*e*q2)
      ei2zl = -(zbl-2.*zb*ql/q)*(ei2-ei1)/q2 - zb*(ei2l-ei1l)/q2

      urr1l = 2.*(c2*s2l+zzl)/rq - (c2*s2+zb2)*rql/rq2
      urr2l = 2.*(c2*q*ql+(zzl*d-zb2*dl)*s2/d3+zb2*s2l/d2)/rq
     &        - (c2*q2 + zb2*s2/d2)*rql/rq2
      urz1l = zbl/q - zb*ql/q2
      urz2l = zbl*s3/d2/q + zb*(2.*s3l*d*q
     &        - s3*(2.*dl*q + d*ql))/d3/q2
      uzr1l = ezl/rq - ezb*rql/rq2
      uzr2l = ezl*s4/d2/rq + 2.0*ezb*s4l/rq/d2-ezb*s4*(rql*d+
     &       2.0*rq*dl)/rq2/d3
      uzz1l = c2*el/q - c2*e*ql/q2
      uzz2l = (el*zb2+2.0*e*zzl)/d2/q - e*zb2*(2.0*dl*q+d*ql)/d3/q2

      urr1rl = 2.0*c2*(el-e*rql/rq)/rq - x6*urr1 - x1*urr1l

      urr2rl = 2.0*(c2*eprl*d2+2.0*c2*epr*dl*d+2.*zzl*(e-emr*s2/d2)
     &   +zb2*(el-(emrl*s2+2.*emr*s2l)/d2+2.*emr*s2*dl/d3))/rq/d2
     &   -2.0*(c2*epr*d2+zb2*(e-emr*s2/d2))*(rql*d+2.0*rq*dl)/rq2/d3
     &   -x6*urr2 - x1*urr2l

      urz1rl = x6*urz1 + x1*urz1l

      urz2rl = 2.0*ezl/q/d2 - 2.0*ezb*(ql*d+2.0*q*dl)/q2/d3 + 
     &  (x1+2.0*emr/d2)*urz2l+(x6+2.*emrl/d2-4.0*emr*dl/d3)*urz2

      uzr1rl = (eprl/q2-2.0*epr*ql/q3)*uzr1+epr*uzr1l/q2

      uzr2rl = 2.0*((2.0*e*el*zb+e2*zbl)-e2*zb*(rql*d+2.0*rq*dl)/rq/d)
     &   /rq/d2 -(eprl/q2-2.0*epr*ql/q3+2.0*emrl/d2-4.0*emr
     &   *dl/d3)*uzr2 -(epr/q2+2.0*emr/d2)*uzr2l

      uzz1rl = x2l*uzz1 + x2*uzz1l

      uzz2rl = (x2l+2.0*emrl/d2-4.0*emr*dl/d3)*uzz2
     &  +(x2+2.0*emr/d2)*uzz2l
 
      urrl = c1*(urr1l*ei1 + urr1*ei1l - urr2l*ei2 - urr2*ei2l)
      urzl = c1*(urz1l*ei1 + urz1*ei1l - urz2l*ei2 - urz2*ei2l)
      uzrl = -c1*(uzr1l*ei1 + uzr1*ei1l - uzr2l*ei2 - uzr2*ei2l)
      uzzl = 2.0*c1*(uzz1l*ei1 + uzz1*ei1l + uzz2l*ei2 + uzz2*ei2l)

      urrrl = c1*(urr1l*ei1r + urr1*ei1rl + urr1rl*ei1 + urr1r*ei1l
     &        -urr2l*ei2r - urr2*ei2rl - urr2rl*ei2 - urr2r*ei2l)
      urzrl = c1*(urz1l*ei1r + urz1*ei1rl - urz1rl*ei1 - urz1r*ei1l
     &        -urz2l*ei2r - urz2*ei2rl + urz2rl*ei2 + urz2r*ei2l)
      uzrrl = -c1*(uzr1l*ei1r + uzr1*ei1rl - uzr1rl*ei1 - uzr1r*ei1l
     &        -uzr2l*ei2r - uzr2*ei2rl - uzr2rl*ei2 - uzr2r*ei2l)
      uzzrl = 2.*c1*(uzz1l*ei1r + uzz1*ei1rl - uzz1rl*ei1 - uzz1r*ei1l
     &        +uzz2l*ei2r + uzz2*ei2rl - uzz2rl*ei2 - uzz2r*ei2l)

      urr1zl = zbl*(x4-urr1/q2)+zb*(-x4*rql/rq-urr1l/q2+2.*urr1*ql/q3)
      urr2zl = zbl*(2.*c2/rq+2.*(s2-2.*zb2*e*r/d2)/d2/rq-urr2/q2)+
     &       zb*(-2.*c2*rql/rq2-urr2l/q2+2.*urr2*ql/q3)
     &          +2.*zb*((2.*s2l-(4.*zzl*e*r+2.*zb2*(el*r+e*rl))/d2 +
     &        4.0*zb2*e*r*dl/d3)/d2/rq-(s2-2.*zb2*e*r/d2)*
     &        (2.0*dl*rq+d*rql)/d3/rq2)
      urz1zl = 2.0*(q*ql - zzl)/q3 - 3.0*(q2-zb2)*ql/q2/q2
      urz2zl = 2.0*(e*el-r*rl-3.0*zzl)/d2/q - (e2-r2-3.0*zb2)*(2.0*dl*q
     & +d*ql)/d3/q2-(4.0*dl/d3+2.0*ql/q3)*urz2*zb+x5*(urz2l*zb+zbl*urz2)
      uzr1zl = el*(q2-zb2)/rq/q2 + e*(2.*(q*ql-zzl)/rq/q2 - (q2-zb2)*
     &        (rl*q + 3.*r*ql)/r2/q2/q2)
      uzr2zl = el*(e2-r2+3.*zb2)/rq/d2 + 2.0*e*(e*el-r*rl+3.0*zzl)/rq/d2
     &        -e*(e2-r2+3.*zb2)*(rql*d+2.0*r*dl*q)/rq2/d3
     &         - zbl*x5*uzr2 + zb*(4.0*dl/d3+2.*ql/q3)*uzr2 - 
     &        zb*x5*uzr2l
      uzz1zl = (zbl/q2 - 2.*zb*ql/q3)*uzz1 + zb/q2*uzz1l
      uzz2zl = 2.0*ezl/d2/q-2.0*ezb*(2.0*dl*q+d*ql)/d3/q2-
     &     zbl*x5*uzz2 + zb*(4.*dl/d3 + 2.0*ql/q3)*uzz2 - zb*x5*uzz2l

      urrzl = c1*(urr1l*ei1z + urr1*ei1zl + urr1zl*ei1 + urr1z*ei1l
     &        -urr2l*ei2z - urr2*ei2zl - urr2zl*ei2 - urr2z*ei2l)
      urzzl = c1*(urz1l*ei1z + urz1*ei1zl + urz1zl*ei1 + urz1z*ei1l
     &        -urz2l*ei2z - urz2*ei2zl + urz2zl*ei2 + urz2z*ei2l)
      uzrzl = -c1*(uzr1l*ei1z + uzr1*ei1zl + uzr1zl*ei1 + uzr1z*ei1l
     &        -uzr2l*ei2z - uzr2*ei2zl - uzr2zl*ei2 - uzr2z*ei2l)
      uzzzl = 2.0*c1*(uzz1l*ei1z + uzz1*ei1zl - uzz1zl*ei1 - uzz1z*ei1l
     &        + uzz2l*ei2z + uzz2*ei2zl + uzz2zl*ei2 + uzz2z*ei2l)

      trr1l = c3*urrrl + c4*(urrl/e - urr*el/e2 + uzrzl)
      trr2l = 0.5*(urrzl + uzrrl)
      tzr1l = c3*uzrzl + c4*(urrl/e - urr*el/e2 + urrrl)
      tzr2l = trr2l
      trz1l = c3*urzrl + c4*(urzl/e - urz*el/e2 + uzzzl)
      trz2l = 0.5*(urzzl + uzzrl)
      tzz1l = c3*uzzzl + c4*(urzl/e - urz*el/e2 + urzrl)
      tzz2l = trz2l

      trrl = 2.*g*(trr1l*xnur + trr1*xnurl + trr2l*xnuz + trr2*xnuzl)
      tzrl = 2.*g*(tzr1l*xnuz + tzr1*xnuzl + tzr2l*xnur + tzr2*xnurl)
      trzl = 2.*g*(trz1l*xnur + trz1*xnurl + trz2l*xnuz + trz2*xnuzl)
      tzzl = 2.*g*(tzz1l*xnuz + tzz1*xnuzl + tzz2l*xnur + tzzz*xnurl)

      end if

 450  continue

      if (ntype.eq.1.or.ntype.eq.2)then

        do k = 1,3

          k2 = 2*k
          k1 = k2-1
          bk(1,k1) = bk(1,k1)+gwpj(k)*urr
          bk(2,k1) = bk(2,k1)+gwpj(k)*urz
          bk(1,k2) = bk(1,k2)+gwpj(k)*uzr
          bk(2,k2) = bk(2,k2)+gwpj(k)*uzz

          if (isens .eq. 2) then
            bkl(1,k1) = bkl(1,k1)+gwpj(k)*urrl+gwpjl(k)*urr
            bkl(2,k1) = bkl(2,k1)+gwpj(k)*urzl+gwpjl(k)*urz
            bkl(1,k2) = bkl(1,k2)+gwpj(k)*uzrl+gwpjl(k)*uzr
            bkl(2,k2) = bkl(2,k2)+gwpj(k)*uzzl+gwpjl(k)*uzz
          end if

          if(k.ne.ising) then

            ak(1,k1) = ak(1,k1)+gwpj(k)*trr
            ak(2,k1) = ak(2,k1)+gwpj(k)*trz
            ak(1,k2) = ak(1,k2)+gwpj(k)*tzr
            ak(2,k2) = ak(2,k2)+gwpj(k)*tzz

            if (isens .eq. 2) then
              akl(1,k1) = akl(1,k1)+gwpj(k)*trrl+gwpjl(k)*trr
              akl(2,k1) = akl(2,k1)+gwpj(k)*trzl+gwpjl(k)*trz
              akl(1,k2) = akl(1,k2)+gwpj(k)*tzrl+gwpjl(k)*tzr
              akl(2,k2) = akl(2,k2)+gwpj(k)*tzzl+gwpjl(k)*tzz
            end if

          end if

        end do

      end if

        end do
      end do

      return
      end
      subroutine bndstr ( ir, sig ) 

c*********************************************************************72
c
cc BNDSTR calculates stresses at a boundary node.
c
c  Modified:
c
c    14 December 2007
c
c  Author:
c
c    Ragevendra Aithal,
c    Jeff Borggaard.
c
      implicit double precision (a-h,o-z)

      double precision dshp(3)
      double precision sig(4)
      double precision u(6)
      double precision x(6)

      common/trnode/ in,rn(300,2) 
      common/geomat/ cord(200,2),ntnode,ntelem,icon(100,3)
      common/temp4/ tt(600),uu(600),tp(600),up(600) 
      common/maprop/ emod,pr,rho,shmod,iplstr,iprob 
      common/inter / nip,neli(200),xip(200),yip(200),nbp,nbn(200)

      data tol / 1.0e-08 /

      nel = neli(ir)
      nn = icon(nel,nbn(ir))

      if ( nbn(ir) .eq. 1 ) then
        eta = 0.0
      else if ( nbn(ir) .eq. 2 ) then
        eta = 0.5
      else if ( nbn(ir) .eq. 3 ) then
        eta = 1.0 
      end if
c
c  Calculate shape function derivatives, coordinate derivatives
c  and local derivatives of displacements at the node of interest
c
      dshp(1) =  4.0 * eta - 3.0
      dshp(2) = -8.0 * eta + 4.0
      dshp(3) =  4.0 * eta - 1.0

      nf1 = icon(nel,1)
      nf2 = icon(nel,2)
      nf3 = icon(nel,3)
      x(1) = cord(nf1,1)
      x(2) = cord(nf1,2)
      x(3) = cord(nf2,1)
      x(4) = cord(nf2,2)
      x(5) = cord(nf3,1)
      x(6) = cord(nf3,2)

      xj1 = dshp(1)*x(1)+dshp(2)*x(3)+dshp(3)*x(5)
      xj2 = dshp(1)*x(2)+dshp(2)*x(4)+dshp(3)*x(6)
      rj = sqrt(xj1*xj1+xj2*xj2)
      rn1 = xj2/rj
      rn2 = -xj1/rj

      u(1) = uu(2*nf1-1)
      u(2) = uu(2*nf1)
      u(3) = uu(2*nf2-1)
      u(4) = uu(2*nf2)	
      u(5) = uu(2*nf3-1)
      u(6) = uu(2*nf3)

      ud1 = dshp(1)*u(1)+dshp(2)*u(3)+dshp(3)*u(5)
      ud2 = dshp(1)*u(2)+dshp(2)*u(4)+dshp(3)*u(6)
      u1 = uu(2*nn-1)
      u2 = uu(2*nn)

      if ( nbn(ir) .eq. 1 ) then
        ng = 3 * ( nel - 1 ) + 1
      else if ( nbn(ir) .eq. 2 ) then
        ng = 3 * ( nel - 1 ) + 2
      else if ( nbn(ir) .eq. 3 ) then
        ng = 3 * ( nel - 1 ) + 3
      end if

      t1 = tt(2*ng-1)
      t2 = tt(2*ng)

      if ( iplstr .ne. 2 ) then
        xlame1 = shmod
        if(pr.eq.0.0)pr = tol
        if(pr.eq.0.5)pr = 0.5-tol
        xlame2 = 2.*shmod*pr/(1.0-2.0*pr)
        c = xlame2+2.0*xlame1
        dela = xlame1*(xj2*xj2-c*xj1*xj1/xlame2)
        a11 = xj2*xj1
        a21 = xlame1*xj1
        a31 = xlame1*xj2
        fac = c/xlame2+a21*xj1*(c*c/xlame2-xlame2)/dela/xlame2
        b11 = -rn2*a11*(c*c/xlame2-xlame2)/dela+rn1*fac
        b12 = rn1
        b13 = -rn2
        b21 = -rn2*fac
        b22 = -rn1*a11*(c*c/xlame2-xlame2)/dela-rn2
        b23 = rn1*fac
        b31 = -rn2*rn2
        b32 = -rn1*rn1
        b33 = rn1*rn2
        delb = rn1*(-rn2*a11*(c*c/xlame2-xlame2)/dela+rn1*fac)-
     &  rn2*rn2
        fac = a21*(xlame2-c*c/xlame2)*ud1/dela+a31*(xlame2-
     &  c*c/xlame2)*ud2/dela
        sig(1) = (b11*t1+b21*t2+b31*fac)/delb
        sig(2) = (b12*t1+b22*t2+b32*fac)/delb
        sig(3) = (b13*t1+b23*t2+b33*fac)/delb
        return
      end if
c
c  Axisymmetric case
c
      xlame1 = shmod

      if(pr.eq.0.0)then
        pr = tol
      end if

      if(pr.eq.0.5) then
        pr = 0.5-tol
      endif

      if(xj1.eq.0.0)then
        xj1 = tol
      end if

      if(xj2.eq.0.0) then
        xj2 = tol
      end if

      xlame2 = 2.0*shmod*pr/(1.-2.*pr)
      c = xlame2+2.0*xlame1
      d = xlame2*xlame2-xlame2*c

      temp1 = xj2*xlame2/d/xj1-xj1*xlame2/d/xj2
      temp2 = xj1*c/xj2/d-xj2*xlame2/xj1/d

      a11 = rn1
      a12 = xlame1*rn2*temp1
      a13 = xlame1*rn2*temp2
      a22 = rn2+xlame1*rn1*temp1
      a23 = xlame1*rn1*temp2
      a31 = 1.0
      a32 = 1.0
      a33 = (c*c-xlame2*xlame2)/d

      dela = a11*(a22*(c*c-xlame2*xlame2)/d-a23)+a12*a23-a22*a13

      b11 = a22*(c*c-xlame2*xlame2)/d-a23
      b12 = a23
      b13 = -a22
      b21 = a13-a12*(c*c-xlame2*xlame2)/d
      b22 = a11*(c*c-xlame2*xlame2)/d-a13
      b23 = a12-a11
      b31 = a12*a23-a22*a13
      b32 = -a11*a23
      b33 = a11*a22

      rr = cord(nn,1)
      if ( rr .eq. 0.0 ) then
        rr = 0.000001
      end if

      r13 = u1*(2.*xlame2+c*(c*c-xlame2*xlame2)/d)/rr
      r21 = ud1-u1*xj1*(-xlame2*xlame2+c*c)/rr/d
      r22 = ud2-u1*xj2/rr

      f1 = t1-r21*xlame1*rn2/xj2-r22*xlame1*rn2/xj1
      f2 = t2-r21*xlame1*rn1/xj2-r22*xlame1*rn1/xj1
      f3 = r13

      x11 = (b11*f1+b21*f2+b31*f3)/dela
      x12 = (b12*f1+b22*f2+b32*f3)/dela
      x13 = (b13*f1+b23*f2+b33*f3)/dela
      x31 = (r21-xj1*xlame2*x12/d+xj1*c*x13/d)/xj2
      x32 = (r22+xj2*xlame2*x12/d-xj2*xlame2*x13/d)/xj1

      sig(1) = x11
      sig(2) = x12
      sig(3) = xlame1*(x31+x32)
      sig(4) = x13

      return
      end
      subroutine bunsnt

c*********************************************************************72
c
cc BUNSNT recovers the unknown sensitivity at the node.
c
c  Modified:
c
c    14 December 2007
c
c  Author:
c
c    Ragevendra Aithal,
c    Jeff Borggaard.
c
      implicit real*8 (a-h,o-z)

      integer*4 nelz, kzone, nkind, nprops

      common/trnode/in,rn(300,2)
      common/geomat/cord(200,2),ntnode,ntelem,icon(100,3)
      common/temp4/tt(600),uu(600),tp(600),up(600)
      common/temp5/tl(600),ul(600)
      common/geomal/cordl(200,2),valuep
      common/maprop/emod,pr,rho,shmod,iplstr,iprob

      e = emod
      pratio  = pr
      tkk  = 1.0
      nkind = iplstr
      g = shmod
      c0 = pratio
      c1 = pratio / ( 1.0 - pratio )
      c2 = e / ( 1.0 - pratio * pratio )
      c3 = 1.0 / sqrt ( 2.0 )

      iout = 6
      write ( iout, '(a)' ) ' '
      write(iout,10) kzone, nelz
   10 format('------- zone ',i2 ,' boundary stress',
     x           ' sensitivities ------- ',/,
     x         1h ,'             calculated for ',i3,' elements',/,
     11h ,'     s11pl        s12pl        s22pl       tmaxl',
     1'         svml',/,
     11h ,'     s11l         s12l         s22l        s1l  ',
     1'         s2l          s3l')
c
c  Here the original cord is read.
c
      rewind 12

      do i = 1,ntnode
        read(12)node,cord(node,1),cord(node,2)
      end do

      do l = 1, ntelem

          lcol3 = (l-1)*3 + 1
          l1 = lcol3
          l2 = l1 + 1
          l3 = l2 + 1
          n1 = icon(l,1)
          n2 = icon(l,2)
          n3 = icon(l,3)
c
c  Here the original coordinate matrix is read in to the mem
c
          x1 = cord(n1,1)
          y1 = cord(n1,2)
          x2 = cord(n2,1)
          y2 = cord(n2,2)
          x3 = cord(n3,1)
          y3 = cord(n3,2)
          ux1 = uu(2*n1-1)
          uy1 = uu(2*n1)
          ux2 = uu(2*n2-1)
          uy2 = uu(2*n2)
          ux3 = uu(2*n3-1)
          uy3 = uu(2*n3)

          x1l = cordl(n1,1)
          y1l = cordl(n1,2)
          x2l = cordl(n2,1)
          y2l = cordl(n2,2)
          x3l = cordl(n3,1)
          y3l = cordl(n3,2)
          ux1l = ul(2*n1-1)
          uy1l = ul(2*n1)
          ux2l = ul(2*n2-1)
          uy2l = ul(2*n2)
          ux3l = ul(2*n3-1)
          uy3l = ul(2*n3)
c
c  Set sign of tractions.
c
          sign = 1.0
          a  = 0.0
          h1a = 4.0*a-3.0
          h2a = -8.0*a+4.0
          h3a = 4.0*a-1.0
          ynu = -(h1a*x1 + h2a*x2 + h3a*x3)
          xnu =   h1a*y1 + h2a*y2 + h3a*y3
          xj = sqrt ( xnu**2 + ynu**2 )
          xj2 = xj * xj
          xnu = xnu / xj
          ynu = ynu / xj
          c = xnu
          s = ynu

          tx  = tt(2*l1-1)
          ty = tt(2*l1)
          s11p = c * tx + s * ty
          s12p = -s * tx + c * ty
          uxa = h1a * ux1 + h2a * ux2 + h3a * ux3
          uya = h1a * uy1 + h2a * uy2 + h3a * uy3
          e22p = ( - s * uxa + c * uya ) / xj
          s22p = c1 * s11p + c2 * e22p

          if ( nkind .eq. 2 .and. x1 .ne. 0 ) then

            s22p = c1 * s11p + c2 *(e22p+c0*ux1/x1)

          end if

          tmax = sqrt ( 0.25*( s11p - s22p )**2 + s12p*s12p )
          s1 = 0.5 * ( s11p + s22p ) + tmax
          s2 = 0.5 * ( s11p + s22p ) - tmax
          s3 = 0.0

          if ( nkind .eq. 1 ) then
            s3 = pratio * ( s11p + s22p )
          end if

          if ( x1 .ne. 0) then
            if ( nkind .eq. 2 ) then
              s3 = pratio*(s11p+s22p)+e*ux1/x1
            end if
          end if

          svm = c3 * sqrt ( (s1-s2)**2 + (s2-s3)**2 + (s3-s1)**2 )

          cc = c * c
          ss = s * s
          sc = s * c
          s11 = s11p * cc + s22p * ss - 2.0 * s12p * sc
          s22 = s11p * ss + s22p * cc + 2.0 * s12p * sc
          s12 = ( s11p - s22p ) * sc + s12p * ( cc - ss )

          dxda = h1a*x1 + h2a*x2 + h3a*x3
          dyda = h1a*y1 + h2a*y2 + h3a*y3
          dxdal = h1a*x1l + h2a*x2l + h3a*x3l
          dydal = h1a*y1l + h2a*y2l + h3a*y3l
          xjl = (dxda*dxdal + dyda*dydal)/xj
          ynul = xjl*dxda/xj2 - dxdal/xj
          xnul = -(xjl*dyda/xj2 - dydal/xj)

          cl = xnul
          sl = ynul

          txl = tl(2*l1-1)
          tyl = tl(2*l1)
          s11pl =  c * txl + s * tyl + cl * tx + sl * ty
          s12pl = -s * txl + c * tyl - sl * tx + cl * ty

          uxal = h1a*ux1l + h2a*ux2l + h3a*ux3l
          uyal = h1a*uy1l + h2a*uy2l + h3a*uy3l
          uxa = h1a*ux1 + h2a*ux2 + h3a*ux3
          uya = h1a*uy1 + h2a*uy2 + h3a*uy3
          u2pa = -s * uxa + c * uya
          u2pal = -s * uxal + c * uyal - sl * uxa + cl * uya
          e22pl = u2pal / xj - xjl * u2pa / (xj*xj)

          if ( x1 .ne. 0.0 ) then
            e33pl = (ux1l*x1 - ux1*x1l)/x1/x1
          end if

          s22pl = c1 * s11pl + c2 * e22pl

          if ( nkind .eq. 2 ) then
            s22pl = c1 * s11pl + c2 *(e22pl+c0*e33pl)
          end if

          tmaxl = ( 0.25*( s11p - s22p )*( s11pl - s22pl )
     &            + s12p*s12pl ) / tmax
          s1l = 0.5 * ( s11pl + s22pl ) + tmaxl
          s2l = 0.5 * ( s11pl + s22pl ) - tmaxl

          if ( nkind .eq. 1 ) then
            s3l = pratio * ( s11pl + s22pl )
          else if ( nkind .eq. 2 ) then
            s3l = pratio*(s11pl+s22pl)+e*e33pl
          else
            s3l = 0.0
          end if

          svml = c3 * ( (s1-s2)*(s1l-s2l) + (s2-s3)*(s2l-s3l)
     &                + (s3-s1)*(s3l-s1l) ) / svm

          a = cc
          b = ss
          d = sc
          al = 2.0*c*cl
          bl = 2.0*s*sl
          dl = sl*c + s*cl
          s11l = s11pl * a  + s22pl * b  - 2.0 * s12pl * d
     &         + s11p  * al + s22p  * bl - 2.0 * s12p  * dl
          s22l = s11pl * b  + s22pl * a  + 2.0 * s12pl * d
     &         + s11p  * bl + s22p  * al + 2.0 * s12p  * dl
          s12l = ( s11pl - s22pl ) * d  + s12pl * ( a  - b  )
     &         + ( s11p  - s22p  ) * dl + s12p  * ( al - bl )
          write ( iout, '(a,i4,a,3i4)' ) 
     &    '  Element = ', i, '  Nodes = ', n1, n2, n3
          write ( iout, '(5g14.6)' ) s11pl, s12pl, s22pl, tmaxl, svml
          write ( iout, '(6g14.6)' ) s11l,  s12l, s22l, s1l, s2l, s3l

          a = 0.50
          h1a = 4.0*a-3.0
          h2a = -8.0*a+4.0
          h3a = 4.0*a-1.
          ynu = -(h1a*x1 + h2a*x2 + h3a*x3)
          xnu =   h1a*y1 + h2a*y2 + h3a*y3
          xj = sqrt (xnu**2 + ynu**2)
          xj2 = xj * xj
          c = xnu / xj
          s = ynu / xj

          tx = tt(2*l2-1)
          ty = tt(2*l2)
          s11p = c * tx + s * ty
          s12p = -s * tx + c * ty
          uxa = h1a * ux1 + h2a * ux2 + h3a * ux3
          uya = h1a * uy1 + h2a * uy2 + h3a * uy3
          e22p = ( - s * uxa + c * uya ) / xj
          s22p = c1 * s11p + c2 * e22p

          if (nkind .eq. 2 .and. x2 .ne. 0.0 ) then
            s22p = c1 * s11p + c2 *(e22p+c0*ux2/x2)
          end if

          tmax = sqrt ( 0.25*( s11p - s22p )**2 + s12p*s12p )
          s1 = 0.5 * ( s11p + s22p ) + tmax
          s2 = 0.5 * ( s11p + s22p ) - tmax
          s3 = 0.0

          if ( nkind .eq. 1 )  then
            s3 = pratio * ( s11p + s22p )
          else if ( nkind .eq. 2 .and. x2.ne.0.0 ) then
            s3 = pratio*(s11p+s22p)+e*ux2/x2
          else
            s3 = 0.0
          end if

          svm = c3 * sqrt ( (s1-s2)**2 + (s2-s3)**2 + (s3-s1)**2 )

          cc = c * c
          ss = s * s
          sc = s * c
          s11 = s11p * cc + s22p * ss - 2.0 * s12p * sc
          s22 = s11p * ss + s22p * cc + 2.0 * s12p * sc
          s12 = ( s11p - s22p ) * sc + s12p * ( cc - ss )

          dxda = h1a*x1 + h2a*x2 + h3a*x3
          dyda = h1a*y1 + h2a*y2 + h3a*y3
          dxdal = h1a*x1l + h2a*x2l + h3a*x3l
          dydal = h1a*y1l + h2a*y2l + h3a*y3l
          xjl = (dxda*dxdal + dyda*dydal)/xj
          ynul = xjl*dxda/xj2 - dxdal/xj
          xnul = -(xjl*dyda/xj2 - dydal/xj)

          cl = xnul
          sl = ynul

          txl  = tl(2*l2-1)
          tyl  = tl(2*l2)
          s11pl =  c * txl + s * tyl + cl * tx + sl * ty
          s12pl = -s * txl + c * tyl - sl * tx + cl * ty

          uxal = h1a*ux1l + h2a*ux2l + h3a*ux3l
          uyal = h1a*uy1l + h2a*uy2l + h3a*uy3l
          uxa = h1a*ux1 + h2a*ux2 + h3a*ux3
          uya = h1a*uy1 + h2a*uy2 + h3a*uy3
          u2pa = -s * uxa + c * uya
          u2pal = -s * uxal + c * uyal - sl * uxa + cl * uya
          e22pl = u2pal / xj - xjl * u2pa / (xj*xj)

          if ( x2 .ne. 0.0E+00 ) then
            e33pl = (ux2l*x2 - ux2*x2l)/x2/x2
          else
            e33pl = 0.0E+00
          end if

          s22pl = c1 * s11pl + c2 * e22pl

          if ( nkind .eq. 2 ) then
            s22pl = c1 * s11pl + c2 *(e22pl+c0*e33pl)
          end if

          tmaxl = ( 0.25*( s11p - s22p )*( s11pl - s22pl )
     &    + s12p*s12pl ) / tmax

          s1l = 0.5 * ( s11pl + s22pl ) + tmaxl
          s2l = 0.5 * ( s11pl + s22pl ) - tmaxl

          if ( nkind .eq. 1 )  then 
            s3l = pratio * ( s11pl + s22pl )
          else if ( nkind .eq. 2 )  then
            s3l = pratio*(s11pl+s22pl)+e*e33pl
          else
            s3l = 0.0
          end if

          svml = c3 * ( (s1-s2)*(s1l-s2l) + (s2-s3)*(s2l-s3l)
     &    + (s3-s1)*(s3l-s1l) ) / svm

          a = cc
          b = ss
          d = sc
          al = 2.0*c*cl
          bl = 2.0*s*sl
          dl = sl*c + s*cl
          s11l = s11pl * a  + s22pl * b  - 2.0 * s12pl * d
     &         + s11p  * al + s22p  * bl - 2.0 * s12p  * dl
          s22l = s11pl * b  + s22pl * a  + 2.0 * s12pl * d
     &         + s11p  * bl + s22p  * al + 2.0 * s12p  * dl
          s12l = ( s11pl - s22pl ) * d  + s12pl * ( a  - b  )
     &         + ( s11p  - s22p  ) * dl + s12p  * ( al - bl )

          write ( iout, '(a,i4,a,3i4)' ) 
     &    '  Element = ', i, '  Nodes = ', n1, n2, n3
          write ( iout, '(5g14.6)' ) s11pl, s12pl, s22pl, tmaxl, svml
          write ( iout, '(6g14.6)' ) s11l, s12l, s22l, s1l, s2l, s3l

          a = 1.0
          h1a = 4.0*a-3.0
          h2a = -8.0*a+4.0
          h3a = 4.0*a-1.
          ynu = -(h1a*x1 + h2a*x2 + h3a*x3)
          xnu =   h1a*y1 + h2a*y2 + h3a*y3
          xj = sqrt ( xnu**2 + ynu**2 )
          xj2 = xj * xj
          c = xnu / xj
          s = ynu / xj

          tx = tt(2*l3-1)
          ty = tt(2*l3)
          s11p = c * tx + s * ty
          s12p = -s * tx + c * ty
          uxa = h1a * ux1 + h2a * ux2 + h3a * ux3
          uya = h1a * uy1 + h2a * uy2 + h3a * uy3
          e22p = ( - s * uxa + c * uya ) / xj

          s22p = c1 * s11p + c2 * e22p
          if (nkind .eq. 2 .and. x3 .ne. 0.0 ) then
            s22p = c1 * s11p + c2 *(e22p+c0*ux3/x3)
          end if

          tmax = sqrt ( 0.25*( s11p - s22p )**2 + s12p*s12p )
          s1 = 0.5 * ( s11p + s22p ) + tmax
          s2 = 0.5 * ( s11p + s22p ) - tmax

          if ( nkind .eq. 1 ) then
            s3 = pratio * ( s11p + s22p )
          else if ( nkind .eq. 2 .and. x3 .ne. 0.0 ) then
            s3 = pratio*(s11p+s22p)+e*ux3/x3
          else
            s3 = 0.0
          end if

          svm = c3 * sqrt ( (s1-s2)**2 + (s2-s3)**2 + (s3-s1)**2 )

          cc = c * c
          ss = s * s
          sc = s * c
          s11 = s11p * cc + s22p * ss - 2.0 * s12p * sc
          s22 = s11p * ss + s22p * cc + 2.0 * s12p * sc
          s12 = ( s11p - s22p ) * sc + s12p * ( cc - ss )

          dxda = h1a*x1 + h2a*x2 + h3a*x3
          dyda = h1a*y1 + h2a*y2 + h3a*y3
          dxdal = h1a*x1l + h2a*x2l + h3a*x3l
          dydal = h1a*y1l + h2a*y2l + h3a*y3l
          xjl = (dxda*dxdal + dyda*dydal)/xj
          ynul = xjl*dxda/xj2 - dxdal/xj
          xnul = -(xjl*dyda/xj2 - dydal/xj)

          cl = xnul
          sl = ynul

          txl = tl(2*l3-1)
          tyl = tl(2*l3)
          s11pl =  c * txl + s * tyl + cl * tx + sl * ty
          s12pl = -s * txl + c * tyl - sl * tx + cl * ty

          uxal = h1a*ux1l + h2a*ux2l + h3a*ux3l
          uyal = h1a*uy1l + h2a*uy2l + h3a*uy3l
          uxa = h1a*ux1 + h2a*ux2 + h3a*ux3
          uya = h1a*uy1 + h2a*uy2 + h3a*uy3
          u2pa = -s * uxa + c * uya
          u2pal = -s * uxal + c * uyal - sl * uxa + cl * uya
          e22pl = u2pal / xj - xjl * u2pa / (xj*xj)

          if (x3.ne.0.0) then
            e33pl = (ux3l*x3 - ux3*x3l)/x3/x3
          end if

          s22pl = c1 * s11pl + c2 * e22pl

          if ( nkind .eq. 2 ) then
            s22pl = c1 * s11pl + c2 *(e22pl+c0*e33pl)
          end if

          tmaxl = ( 0.25*( s11p - s22p )*( s11pl - s22pl )
     &            + s12p*s12pl ) / tmax
          s1l = 0.5 * ( s11pl + s22pl ) + tmaxl
          s2l = 0.5 * ( s11pl + s22pl ) - tmaxl

          if ( nkind .eq. 1 )  then
            s3l = pratio * ( s11pl + s22pl )
          else if ( nkind .eq. 2 )  then
            s3l = pratio*(s11pl+s22pl)+e*e33pl
          else
            s3l = 0.0
          end if

          svml = c3 * ( (s1-s2)*(s1l-s2l) + (s2-s3)*(s2l-s3l)
     &                + (s3-s1)*(s3l-s1l) ) / svm

          a = cc
          b = ss
          d = sc
          al = 2.0D+00*c*cl
          bl = 2.0D+00*s*sl
          dl = sl*c + s*cl
          s11l = s11pl * a  + s22pl * b  - 2.0 * s12pl * d
     &         + s11p  * al + s22p  * bl - 2.0 * s12p  * dl
          s22l = s11pl * b  + s22pl * a  + 2.0 * s12pl * d
     &         + s11p  * bl + s22p  * al + 2.0 * s12p  * dl
          s12l = ( s11pl - s22pl ) * d  + s12pl * ( a  - b  )
     &         + ( s11p  - s22p  ) * dl + s12p  * ( al - bl )

          write ( iout, '(a,i4,a,3i4)' ) 
     &    '  Element = ', i, '  Nodes = ', n1, n2, n3
          write ( iout, '(5g14.6)' ) s11pl, s12pl, s22pl, tmaxl, svml
          write ( iout, '(6g14.6)' ) s11l, s12l, s22l, s1l, s2l, s3l
      end do

      return
      end
      subroutine force ( gmat, fmat, gmatl, fmatl, nsize, ncolg ) 

c*********************************************************************72
c
cc FORCE sets the right hand side for the sensitivity equation.
c
c  Discussion:
c
c    this subroutine finds the righthand side vector needed for the 
c    calculation of the sensitivity.
c
c    if the trangular factorization is used then the previous fmat
c    could be saved and used again.
c
c  Modified:
c
c    14 December 2007
c
c  Author:
c
c    Ragevendra Aithal,
c    Jeff Borggaard.
c
      implicit real*8 (a-h,o-z) 

      dimension gmat(nsize,ncolg),fmat(nsize,nsize+1) 
      dimension gmatl(nsize,ncolg),fmatl(nsize,nsize+1) 

      common/geomat/ cord(200,2),ntnode,ntelem,icon(100,3)
      common/maprop/emod,pr,rho,shmod,iplstr,iprob
      common/geomal/ cordl(200,2),valuep
      common/temp4/ t(600),u(600),tp(600),up(600)
      common/temp6 /tpl(600),upl(600)
c
c  The difference of the g and f matrix is found here.
c
      do i = 1, nsize
        do j = 1,ncolg
          gmatl(i,j) = gmat(i,j)-gmatl(i,j)
        end do
        do k = 1, nsize
          fmatl(i,k) = fmat(i,k)-fmatl(i,k) 
        end do
      end do

      do i = 1, ncolg
        tempt = tpl(i)
        tpl(i) = (tp(i) - tpl(i))/valuep
        tp(i) = tempt
      end do

       do i = 1, nsize
        tempu = upl(i)
        upl(i) = (up(i) - upl(i))/valuep
        up(i) = tempu
      end do
c
c  Actual force vector for the calculation of the sensitivity is found
c
      do i = 1, nsize
        fmat(i,nsize+1) = 0.0D+00
        do j = 1, nsize
          fmat(i,(nsize+1)) = fmatl(i,j)
     &    / valuep*(up(j)-u(j))+fmat(i,(nsize+1))
        end do
        do j = 1, ncolg
          fmat(i,(nsize+1)) = fmat(i,(nsize+1))
     &    + gmatl(i,j)/valuep*(t(j)-tp(j))
        end do
      end do

      rewind 13
      read(13) ((fmat(i,j),j = 1,nsize),i = 1,nsize)
      read(13) ((gmat(i,j),j = 1,ncolg),i = 1,nsize)

      do i = 1, nsize
        do j = 1, ncolg
          fmat(i,(nsize+1)) = fmat(i,(nsize+1))-gmat(i,j)*tpl(j)
        end do
      end do

      do i = 1,nsize
        do j = 1,nsize
          fmat(i,(nsize+1)) = fmat(i,(nsize+1))+fmat(i,j)*upl(j)
        end do
      end do

      read(13)((fmat(i,j),j = 1,nsize),i = 1,nsize)

      return
      end 
      subroutine force1 ( gmat, fmat, gmatl, fmatl, nsize, ncolg ) 

c*********************************************************************72
c
cc FORCE1 finds the right hand side of the sensitivity equation.
c
c  Discussion:
c
c    This subroutine finds the righthand side vector needed for the 
c    calculation of the sensitivity.
c
c    If the triangular factorization is used then the previous FMAT
c    could be saved and used again.
c
c  Modified:
c
c    14 December 2007
c
c  Author:
c
c    Ragevendra Aithal,
c    Jeff Borggaard.
c
      implicit real*8 (a-h,o-z) 

      dimension fmat(nsize,nsize+1) 
      dimension fmatl(nsize,nsize+1) 
      dimension gmat(nsize,ncolg)
      dimension gmatl(nsize,ncolg)

      common/geomat/ cord(200,2),ntnode,ntelem,icon(100,3)
      common/maprop/emod,pr,rho,shmod,iplstr,iprob
      common/geomal/ cordl(200,2),valuep
      common/temp4/ t(600),u(600),tp(600),up(600)
      common/temp6 /tpl(600),upl(600)
c
c    actual force vector for the calculation fo the sensitivity is found
c    here.
c
      rewind 13
      read(13)((fmat(i,j),j = 1,nsize),i = 1,nsize)
      read(13)((gmat(i,j),j = 1,ncolg),i = 1,nsize)

      do i = 1,nsize
        fmat(i,nsize+1) = 0.0D+00
        do j = 1,nsize
          fmat(i,(nsize+1)) = -fmatl(i,j)*(u(j)-up(j))+fmat(i,(nsize+1))
        end do
        do j = 1,ncolg
          fmat(i,(nsize+1)) = fmat(i,(nsize+1))+gmatl(i,j)*(t(j)-tp(j))
        end do
      end do

      do i = 1,nsize
        do j = 1,ncolg
          fmat(i,(nsize+1)) = fmat(i,(nsize+1))-gmat(i,j)*tpl(j)
        end do
      end do

      do i = 1,nsize
        do j = 1,nsize
          fmat(i,(nsize+1)) = fmat(i,(nsize+1))+fmat(i,j)*upl(j)
        end do
      end do

      read(13)((fmat(i,j),j = 1,nsize),i = 1,nsize)

      return
      end 
      subroutine gaussp 

c*********************************************************************72
c
cc GAUSSP sets the Gauss points and weights for the interval [0,1].
c
c  Modified:
c
c    14 December 2007
c
c  Author:
c
c    Ragevendra Aithal,
c    Jeff Borggaard.
c
      implicit real*8 (a-h,o-z) 

      integer i
      integer j

      common/gauswp/gp(12,12),gpw(12,12),gsp(8,8),gspw(8,8) 

      do i = 1,12
        do j = 1,12
          gp(i,j) = 0.0D+00
          gpw(i,j) = 0.0D+00
        end do
      end do

      gp(1,1) = 0.5D+00
      gpw(1,1) = 1.0D+00

      gp(2,1) = 0.2113249D+00 
      gp(2,2) = 0.7886751D+00 
      gpw(2,1) = 0.5D+00
      gpw(2,2) = 0.5D+00

      gp(3,1) = 0.8872983D+00 
      gp(3,2) = 0.5D+00 
      gp(3,3) = 0.1127017D+00
      gpw(3,1) = 0.27777778D+00
      gpw(3,2) = 0.44444444D+00
      gpw(3,3) = 0.27777778D+00

      gp(4,1) = 0.06943184D+00
      gp(4,2) = 0.33000948D+00
      gp(4,3) = 0.66999052D+00
      gp(4,4) = 0.93056816D+00
      gpw(4,1) = 0.17392742D+00
      gpw(4,2) = 0.32607258D+00 
      gpw(4,3) = gpw(4,2) 
      gpw(4,4) = gpw(4,1) 

      gp(5,1) = .0469101D+00
      gp(5,2) = .2307654D+00
      gp(5,3) = .5D+00
      gp(5,4) = .7692347D+00
      gp(5,5) = .9530899D+00
      gpw(5,1) = .1184634D+00
      gpw(5,2) = .2393143D+00 
      gpw(5,3) = .2844444D+00
      gpw(5,4) = gpw(5,2) 
      gpw(5,5) = gpw(5,1) 

      gp(6,1) = .0337653D+00
      gp(6,2) = .1693954D+00
      gp(6,3) = .3806905D+00
      gp(6,4) = .6193096D+00
      gp(6,5) = .8306047D+00
      gp(6,6) = .9662348D+00
      gpw(6,1) = .0856622D+00
      gpw(6,2) = .1803808D+00
      gpw(6,3) = .233957D+00
      gpw(6,4) = gpw(6,3) 
      gpw(6,5) = gpw(6,2) 
      gpw(6,6) = gpw(6,1) 

      gp(7,1) = .0254461D+00
      gp(7,2) = .1292345D+00
      gp(7,3) = .2970775D+00
      gp(7,4) = .5D+00
      gp(7,5) = .7029226D+00
      gp(7,6) = .8707656D+00
      gp(7,7) = .974554D+00
      gpw(7,1) = .0647425D+00
      gpw(7,2) = .1398527D+00
      gpw(7,3) = .190915D+00
      gpw(7,4) = .2089796 D+00
      gpw(7,5) = gpw(7,3) 
      gpw(7,6) = gpw(7,2) 
      gpw(7,7) = gpw(7,1) 

      gp(8,1) = 0.019855D+00
      gp(8,2) = 0.101667D+00
      gp(8,3) = 0.237234D+00
      gp(8,4) = 0.4082825D+00
      gp(8,5) = 0.5917175D+00
      gp(8,6) = 0.7627660D+00
      gp(8,7) = 0.898333D+00
      gp(8,8) = 0.980145D+00
      gpw(8,1) = 0.050614D+00
      gpw(8,2) = 0.1111905D+00
      gpw(8,3) = 0.1568533D+00
      gpw(8,4) = 0.1813419D+00
      gpw(8,5) = gpw(8,4) 
      gpw(8,6) = gpw(8,3) 
      gpw(8,7) = gpw(8,2) 
      gpw(8,8) = gpw(8,1) 

      gp(12,1) = .0092197D+00
      gp(12,2) = .0479414D+00
      gp(12,3) = .1150487D+00
      gp(12,4) = .2063411D+00
      gp(12,5) = .3160843D+00
      gp(12,6) = .4373833D+00
      gp(12,7) = .5626167D+00
      gp(12,8) = .6839157D+00
      gp(12,9) = .793659D+00
      gp(12,10) = .8849513D+00
      gp(12,11) = .9520586D+00
      gp(12,12) = .9907803D+00
      gpw(12,1) = .0235877D+00
      gpw(12,2) = .0534697D+00
      gpw(12,3) = .0800392D+00
      gpw(12,4) = .1015837D+00
      gpw(12,5) = .1167463D+00
      gpw(12,6) = .1245735D+00
      gpw(12,7) = gpw(12,6) 
      gpw(12,8) = gpw(12,5) 
      gpw(12,9) = gpw(12,4) 
      gpw(12,10) = gpw(12,3)
      gpw(12,11) = gpw(12,2)
      gpw(12,12) = gpw(12,1)
c
c   setting gaussian sampling points and weightages
c   for -1 to +1 end limits
c
      gsp(2,2) = 0.5773502691D+00
      gspw(2,2) = 1.0000000000D+00

      gsp(3,3) = 0.7745966692D+00
      gsp(3,2) = 0.0D+00
      gspw(3,3) = 0.5555555555D+00
      gspw(3,2) = 0.8888888888D+00

      gsp(4,4) = 0.8611363115D+00
      gsp(4,3) = 0.3399810435D+00
      gspw(4,4) = 0.3478548451D+00
      gspw(4,3) = 0.6521451548D+00

      gsp(5,5) = 0.9061798459D+00
      gsp(5,4) = 0.5384693101D+00
      gsp(5,3) = 0.0D+00
      gspw(5,5) = 0.2369268850D+00
      gspw(5,4) = 0.4786286704D+00
      gspw(5,3) = 0.5688888888D+00

      do i = 2,8 

        j = 0 

 1020   continue

        j = j + 1 
        jp = i-j+1
        ia = i / 2
        xa = i / 2

        if (ia .lt. xa) then
          ia = ia + 1 
        end if

        if ( ia .lt. jp ) then
          gsp(i,j) = - gsp(i,jp)
          gspw(i,j) = gspw(i,jp)
          go to 1020
        end if

      end do

      return
      end 
      subroutine input ( nsize, ncolg )

c*********************************************************************72
c
cc INPUT reads input defining the problem to be solved.
c
c  Modified:
c
c    14 December 2007
c
c  Author:
c
c    Ragevendra Aithal,
c    Jeff Borggaard.
c
      implicit real*8 (a-h,o-z) 

      common / corner / intype(100)
      common / maprop / emod,pr,rho,shmod,iplstr,iprob
      common / minlth / rmaxl 
      common / geomat / cord(200,2),ntnode,ntelem,icon(100,3)
      common / boundc / bc(100,6),ibtype(100,6) 
      common / segmen/ nseg,isnod(10),isele(10) 
      common / fixdis / ifxd(200,2),fxdp(200,2)
      common / body  / omega,accn
      common / mesh / x1,y1,x2,y2,x3,y3,x4,y4
      common / inter / nip,neli(200),xip(200),yip(200),nbp,nbn(200)

      character*80 title(8)
      data toli / 0.0000001D+00 /
      write(*,99)
c 
c  input number of lines in the title ; maximum of 8 lines permitted;
c  each line can contain up to 80 characters
c
      read(15,*)nline
c
c  read title 
c
      do i = 1,nline
        read(15,105) title(i)
        write(*,106) title(i)
      end do
c 
c  read flags for the type of problem
c  iplstr = 0  plane strain
c  iplstr = 1  plane stress
c  iplstr = 2  axisymmetric elasticity
c  iplstr = 3  axisymmetric potential
c
c  iprob = 0    no body forces specified
c  iprob = 1    gravitational body forces included
c  iprob = 2    cetrifugal body forces included
c
      read(15,*) iplstr,iprob

      if(iplstr .eq. 0 ) then
        write ( *,109) 
      else if(iplstr.eq.1) then
        write ( *,110) 
      else if(iplstr.eq.2) then
        write ( *,117)
      else if(iplstr.eq.3) then
        write ( *,126)
      end if

      if(iprob.eq.0 ) write ( *,112)
      if(iprob.eq.1 ) write ( *,113)
      if(iprob.eq.2 ) write ( *,114)
c 
c  emod = elastic modulus,  pr = poission's ratio,  rho = mass density
c
      read(15,*) emod,pr,rho 
      write(*,101) emod,pr,rho
c
c   read acceleration due to gravity
c
      if (iprob.eq.1) then
        read(15,*)  accn
        write(*,118) accn
c 
c  read rotational speed
c
      else if (iprob.eq.2) then
        read(15,*)  omega
        write(*,119) omega
      end if
c
c  ntnode = total no. of nodes in the entire boundary discretization
c  ntelem = total no. of elements in the entire boundary discretization 
c
      read(15,*) ntnode,ntelem 
      write ( *,102) ntnode,ntelem

      i = 1
      l = 0 
      ntn = 0 
      nte = 0 
c 
c   read no. of surfaces = nseg
c                        = 1 for a singly connected region
c                        > 1 for a multiply connected region
c 
      read(15,*) nseg

      do iseg = 1,nseg
c
c    the following data concern each bounding surface
c 
c    read no. of nodes and elements in each surface 
c    nnode = no. of nodes , nelem = no. of elements 
c 
        read(15,*) nnode,nelem
        isnod(iseg) = nnode
        isele(iseg) = nelem
        nubn  = ntn + nnode 
c
c   irflag = 1   read the coordinates of every node
c
        read(15,*) irflag
c
c  read coordinates of nodes
c  inod = node number
c  cord(inod,1) = x or r-cordinate of inod
c  cord(inod,2) = y or z-cordinate of inod
c
        do ik = 1,nnode
          read(15,*)  inod,cord(inod,1),cord(inod,2)
        end do

        in1 = ntn + 1
        ie1 = nte + 1
        ie2 = nte + nelem

        do ic = ie1,ie2
          do j = 1,3 
            icon(ic,j) = l + j 
          end do
          l = l + 2
        end do

        if (nnode .ne. 2*nelem+1) then
          icon(ie2,3) = in1
        end if

        ntn = ntn + nnode 
        nte = nte + nelem 
        i = i + 1

      end do

      write ( *,*)

      if(iplstr.eq.2)then
        write ( *,127)
      else
        write ( *,107)
      end if

      do i = 1,ntnode
        write ( *,103) i,(cord(i,j),j = 1,2)
      end do
c 
c  print connectivity
c
      write ( *,98)
      write ( *,108)
      do i = 1,ntelem
        write ( *,104)  i,(icon(i,j),j = 1,3)
      end do

      do j = 1,6 
        do i = 1,ntelem
          ibtype(i,j) = 1 
          bc(i,j)  = 0.0D+00
        end do
      end do

      do k = 1,2
        do i = 1,ntnode
          ifxd(i,k) = 0
          fxdp(i,k) = 0.0D+00
        end do
      end do
c
c   read boundary conditions
c   nbcc = total number of specified boundary conditions
c   nel = element no. ; nb = node no. (1,2,or 3) 
c   it1 = type of bc in x1 or r direction  also used for potential
c   it2 = type of bc in x2 or z direction 
c   it1/it2 = 1 , traction specified in (x1/x2 or r/z) direction
c   it1/it2 = 2 , displacement specified in (x1/x2 or r/z) direction 
c   val1 = value of specified bc corresponding to it1 
c   val2 = value of specified bc corresponding to it2
c   default boundary condition is traction free
c 
      read(15,*)nbcc 

      do i = 1,nbcc
         read(15,*) nel,nb,it1,val1,it2,val2 
           nb1 = 2*nb - 1 
           nb2 = nb1 + 1
           nf = icon(nel,nb)
           bc(nel,nb1) = val1 
           bc(nel,nb2) = val2 
           if (it1.eq.2) then 
             ibtype(nel,nb1) = it1 
             fxdp(nf,1) = val1
           end if
           if (it2.eq.2) then
             ibtype(nel,nb2) = it2 
             fxdp(nf,2) = val2
           end if
      end do

      call setup ( iplstr, nsize, ncolg, rmaxl, intype )

      write ( *,115)
      write ( *,116)
      write ( *,120)
      write ( *,121)
      write ( *,122)
      do i = 1,ntelem
        write ( *,123)i,(ibtype(i,j),j = 1,2),(bc(i,j),j = 1,2) 
      end do
      write ( *,124)
      write ( *,120)
      write ( *,121)
      write ( *,122)
      do i = 1,ntelem
        write ( *,123)i,(ibtype(i,j),j = 3,4),(bc(i,j),j = 3,4)
      end do
      write ( *,125)
      write ( *,120)
      write ( *,121)
      write ( *,122)
      do i = 1,ntelem
        write ( *,123)i,(ibtype(i,j),j = 5,6),(bc(i,j),j = 5,6)
      end do
c
c  nip = total number of 'interior points' where stresses and displacements are desired 
c  xip(),yip() = global coordinates of an interior point 
c 
      read(15,*)nip

      do k = 1,nip
        read(15,*)xip(k),yip(k)
      end do
c
c  nbp = total number of nodal points at which stresses are to be calculated
c  neli() = no. of the element on which the nodal point of interest lies 
c  nbn() = order of the above node 1 , 2 or 3 on the element
c
      read(15,*)nbp

      do k = 1,nbp
        read(15,*)neli(k),nbn(k)        
      end do

 98   format(//'details of connectivity :'/)
 99   format('problem particulars :'/)
 101  format(//'elastic modulus,poisson''s ratio,rho :',3(e10.5,3x)
     *)
 102  format(' number of nodes , number of elements :',2(i5,5x)) 
 103  format(9x,i5,5x,2(f10.4,5x))
 104  format(5x,i3,8x,3(i3,5x)) 
 108  format(2x,11helement no.,2x,6hnode a,2x,6hnode b,2x,6hnode c)
 105  format(a80)
 106  format(1x,a80)
 107  format(//10x,' node no.',5x,'x-cord',9x,'y-cord'/) 
 127  format(//10x,' node no.',5x,'r-cord',9x,'z-cord'/) 
 109  format(1x,'specified problem type - plane strain ') 
 110  format(1x,'specified problem type - plane stress ') 
 112  format(1x,'no body forces')
 113  format(1x,'gravity forces present')
 114  format(1x,'centrifugal forces present')
 115  format(//1x,30hboundary conditions provided :)
 116  format(/1x,11helement no.,11x,6hnode a)
 117  format(1x,'specified problem type - axisymmetric') 
 118  format(/1x,'acceleration due to gravity :',f10.5/)
 119  format(/1x,'speed for centrifugal analysis :',f12.5/)
 120  format(1x,11('_'),11x,6('_'))
 121  format(13x,8hbc types)
 122  format(12x,11hdir.1 dir.2,7x,4hval1,8x,4hval2)
 123  format(4x,i3,7x,i1,5x,i1,3x,f10.4,2x,f10.4)
 124  format(/1x,11helement no.,11x,6hnode b)
 125  format(/1x,11helement no.,11x,6hnode c)
 126  format(1x,'specified problem type - axipotential')
      return  
      end 
      subroutine inputl

c*********************************************************************72
c
cc INPUTL reads the sensitivity information.
c
c  Modified:
c
c    14 December 2007
c
c  Author:
c
c    Ragevendra Aithal,
c    Jeff Borggaard.
c
      implicit real*8(a-h,o-z)

      common / maprop / emod,pr,rho,shmod,iplstr,iprob
      common / geomat / cord(200,2),ntnode,ntelem,icon(100,3) 
      common / geomal / cordl(200,2),valuep

      open ( unit = 12, file = 'sensitivity.dat', form = 'unformatted' )

      toli = 0.000001D+00

      do inod = 1,ntnode
        write(12) inod,cord(inod,1),cord(inod,2)
      end do

      read(15,*)valuep

      do ik = 1,ntnode
         read(15,*)  inod,cordl(inod,1),cordl(inod,2)  
         do i1 = 1,2
           temp = (cordl(inod,i1)-cord(inod,i1))/valuep
           cord(inod,i1) = cordl(inod,i1)
           cordl(inod,i1) = temp
         end do
       end do

      do i = 1,ntnode
        write ( *, '(9x,i5,5x,2(f10.4,5x))' ) i,(cordl(i,j),j = 1,2)
      end do

      return  
      end 
      subroutine inputl1

c*********************************************************************72
c
cc INPUTL1 reads the sensitivity information when ISENS = 2.
c
c  Modified:
c
c    14 December 2007
c
c  Author:
c
c    Ragevendra Aithal,
c    Jeff Borggaard.
c
      implicit real*8(a-h,o-z)

      common / maprop / emod,pr,rho,shmod,iplstr,iprob
      common / geomat / cord(200,2),ntnode,ntelem,icon(100,3) 
      common / geomal / cordl(200,2),valuep

      open(unit = 12,file = 'sensitivity.dat',form = 'unformatted')

      toli = 0.000001D+00

      do inod = 1,ntnode
        write(12) inod,cord(inod,1),cord(inod,2)
      end do

      read(15,*)valuep

       do ik = 1,ntnode
         read(15,*)  inod,cordl(inod,1),cordl(inod,2)
         do i1 = 1,2
           cordl(inod,i1) = (cordl(inod,i1)-cord(inod,i1))/valuep
         end do
       end do

      do i = 1,ntnode
        write ( *, '(9x,i5,5x,2(f10.4,5x))' ) i,(cordl(i,j),j = 1,2)
      end do
 
      return  
      end 
      subroutine integra1 ( isens, gmat, fmat, gmatl, fmatl, nsize, 
     &  ncolg ) 

c*********************************************************************72
c
cc INTEGRA1 carries out the integration needed for assembly.
c
c  Modified:
c
c    14 December 2007
c
c  Author:
c
c    Ragevendra Aithal,
c    Jeff Borggaard.
c
      implicit real*8 (a-h,o-z) 

      dimension gmat(nsize,ncolg),fmat(nsize,nsize) 
      dimension gmatl(nsize,ncolg),fmatl(nsize,nsize)
      dimension xy(3,2),gk(2,6),fk(2,6),sgp(4,12),sl(4),tk(4,6) 
      dimension xyl(3,2),gkl(2,6),fkl(2,6)
      dimension rnl(300,2)

      common /geomat/ cord(200,2),ntnode,ntelem,icon(100,3)
      common /geomal/ cordl(200,2),valuep
      common/intcon/cons1,cons2,cons3,cons4 
      common /maprop/ emod,pr,rho,shmod,iplstr,iprob
      common/trnode/in,rn(300,2) 
      common/minlth/ rmaxl 
c 
c  for implicit differentiation sensitivities, normal sensitivities
c  are calculated for inflation mode calculation for a,l b,l
c
      if (isens .eq. 2) then

        do i = 1,ntelem
          do j = 1,3

            in = (i-1)*3 + j
            n1 = icon(i,1)
            n2 = icon(i,2)
            n3 = icon(i,3)
            r1 = cord(n1,1)
            r2 = cord(n2,1)
            r3 = cord(n3,1)
            z1 = cord(n1,2)
            z2 = cord(n2,2)
            z3 = cord(n3,2)
            r1l = cordl(n1,1)
            r2l = cordl(n2,1)
            r3l = cordl(n3,1)
            z1l = cordl(n1,2)
            z2l = cordl(n2,2)
            z3l = cordl(n3,2)

            if (j.eq.1) then
              s = 0.0D+00
            else if (j.eq.2) then
              s = 0.5D+00
            else if (j.eq.3) then
              s = 1.0D+00
            end if

            ra = (4.0*s-3.0)*r1 + (-8.0*s+4.0)*r2 
     &        + (4.0*s-1.0D+00)*r3
            za = (4.0*s-3.0)*z1 + (-8.0*s+4.0)*z2 
     &        + (4.0*s-1.0D+00)*z3
            ral = (4.0*s-3.0)*r1l + (-8.0*s+4.0)*r2l 
     &        + (4.0*s-1.0D+00)*r3l
            zal = (4.0*s-3.0)*z1l + (-8.0*s+4.0)*z2l 
     &        + (4.0*s-1.0D+00)*z3l
            xj = dsqrt(ra*ra + za*za)
            xjl = (ra*ral + za*zal)/xj
            rnl(in,1) = -xjl*za/xj/xj + zal/xj
            rnl(in,2) = xjl*ra/xj/xj - ral/xj

          end do
        end do

      end if

      toli = 1.0D-02

      do i = 1,nsize

        do j = 1,nsize 
          fmat(i,j) = 0.0D+00
          if (isens.eq.2) then
            fmatl(i,j) = 0.0D+00
          end if
        end do

        do k = 1,ncolg 
          gmat(i,k) = 0.0D+00
          if (isens.eq.2) then
            gmatl(i,k) = 0.0D+00
          end if
        end do

      end do
c 
c  integration loop
c 
c  loop over all the nodes
c
      do i = 1,ntnode

        xpt = cord(i,1) 
        ypt = cord(i,2) 
        xpl = cordl(i,1)
        ypl = cordl(i,2)
c
c  loop over all the elements
c
        do j = 1,ntelem

          ising = 0 

          do k = 1,2 
            do km = 1,6
              gk(k,km) = 0.0D+00
              fk(k,km) = 0.0D+00
              gkl(k,km) = 0.0D+00
              fkl(k,km) = 0.0D+00
            end do
          end do

          node1 = icon(j,1) 
          node2 = icon(j,2) 
          node3 = icon(j,3) 

          do k = 1,2 
            xy(1,k) = cord(node1,k)
            xy(2,k) = cord(node2,k)
            xy(3,k) = cord(node3,k)
            xyl(1,k) = cordl(node1,k)
            xyl(2,k) = cordl(node2,k)
            xyl(3,k) = cordl(node3,k)
          end do
c
c  Check for singularity 
c
          r1  = sqrt((xy(1,1)-xpt)**2 + (xy(1,2)-ypt)**2 ) 
          r2  = sqrt((xy(2,1)-xpt)**2 + (xy(2,2)-ypt)**2 ) 
          r3  = sqrt((xy(3,1)-xpt)**2 + (xy(3,2)-ypt)**2 ) 

          if (r1.lt.toli) then
            ising = 1
          else if (r2.lt.toli) then
            ising = 2 
          else if (r3.lt.toli) then
            ising = 3 
          end if
c
c  Order of gauss points
c
          call autint(ngp,xpt,ypt,xy,nsub,rmaxl) 
c
c  Set subsegments
c  some information is given in pages 378-381 of text
c  also see references 1-4 mentioned at the end of chapter
c  fifteen of text

          if ( ising .eq. 0 ) then
            call norset ( r1, r3, ngp, nsub, sgp, sl )
          else if ( 0 .lt. ising ) then
            call sinset ( ngp, nsub, sgp, sl, ising ) 
          end if

          if (iplstr.eq.2) then
            call axkrnl(xpt,ypt,xy,gk,fk,gkl,fkl,ngp,nsub,sgp,sl,
     &      1,ising,isens,tk,tk,xyl,xpl,ypl)
          else
            call plkrnl(xpt,ypt,xy,gk,fk,ngp,nsub,sgp,sl,ising)
          end if
c
c  assembly of g and f matrices
c  see page 381 of text
c
          j1 = 2*i - 2
          k1 = 2*node1 - 2
          k2 = 2*node2 - 2
          k3 = 2*node3 - 2
          l1 = 6*(j-1)
          l2 = l1 + 2
          l3 = l2 + 2

          do k = 1,2 
            do jk = 1,2
              gmat(j1+jk,l1+k) = gk(jk,k)
              gmat(j1+jk,l2+k) = gk(jk,k+2)
              gmat(j1+jk,l3+k) = gk(jk,k+4)
              fmat(j1+jk,k1+k) = fmat(j1+jk,k1+k) + fk(jk,k)
              fmat(j1+jk,k2+k) = fmat(j1+jk,k2+k) + fk(jk,k+2)
              fmat(j1+jk,k3+k) = fmat(j1+jk,k3+k) + fk(jk,k+4)
              if (isens .eq. 2) then
                gmatl(j1+jk,l1+k) = gkl(jk,k)
                gmatl(j1+jk,l2+k) = gkl(jk,k+2)
                gmatl(j1+jk,l3+k) = gkl(jk,k+4)
                fmatl(j1+jk,k1+k) = fmatl(j1+jk,k1+k) + fkl(jk,k)
                fmatl(j1+jk,k2+k) = fmatl(j1+jk,k2+k) + fkl(jk,k+2)
                fmatl(j1+jk,k3+k) = fmatl(j1+jk,k3+k) + fkl(jk,k+4)
              end if
            end do
          end do

        end do
      end do
c
c   diagonal terms of [f] matrix 
c   by using rigid body translations.
c   see page nos. 96-97 of text
c
      do i = 1,nsize,2 

        sum11 = 0.0D+00
        sum12 = 0.0D+00
        suml12 = 0.0D+00
        sum21 = 0.0D+00
        sum22 = 0.0D+00
        suml22 = 0.0D+00

        do j = 1,nsize,2
          sum11 = sum11 + fmat(i,j)
          sum12 = sum12 + fmat(i,j+1)
          suml12 = suml12 + fmatl(i,j+1)
          sum21 = sum21 + fmat(i+1,j)
          sum22 = sum22 + fmat(i+1,j+1)
          suml22 = suml22 + fmatl(i+1,j+1)
        end do

        if(iplstr.eq.2)then
          sum11 = 0.0D+00
          sum21 = 0.0D+00
          suml11 = 0.0D+00
          suml21 = 0.0D+00
        end if

        fmat(i,i) = - sum11 
        fmat(i,i+1) = - sum12
        fmat(i+1,i) = - sum21
        fmat(i+1,i+1) = - sum22

        if (isens .eq. 2) then
          fmatl(i,i) = - suml11 
          fmatl(i,i+1) = - suml12
          fmatl(i+1,i) = - suml21
          fmatl(i+1,i+1) = - suml22
        end if

      end do
c
c  jump terms for axisymmetric problems
c  see paper entitled "new axisymmetric bem formulation for body
c  forces using particular integrals" - d. p. henry,jr.,
c  d. a. pape and p. k. banerjee, in journal of engineering
c  mechanics(asce), vol.113, no.5, may 1987, pp.671-688,
c  and references 12 and 15 therein. 
c
      if (iplstr.eq.2)then

        c1 = emod/((1.0D+00-2.0D+00*pr)*(1.0D+00+pr))
        c2 = 2.0D+00*pr*c1

        do l = 1,ntelem
          do m = 1,2

            n = icon(l,m)
            l2 = n*2
            l1 = l2-1
            ri1 = cord(n,1)
            sumg1 = 0.0D+00
            sumg2 = 0.0D+00
            sumf1 = 0.0D+00
            sumf2 = 0.0D+00
            sumlg1 = 0.0D+00
            sumlg2 = 0.0D+00
            sumlg3 = 0.0D+00
            sumlg4 = 0.0D+00
            sumlf1 = 0.0D+00
            sumlf2 = 0.0D+00
            k = 1

            do i = 1,ntelem
              do j = 1,3
                k2 = k*2
                k1 = k2-1
                in = (i-1)*3+j
                k = k+1
                t1 = c1*rn(in,1)
                t2 = c2*rn(in,2)
                sumg1 = sumg1+gmat(l1,k1)*t1+gmat(l1,k2)*t2
                sumg2 = sumg2+gmat(l2,k1)*t1+gmat(l2,k2)*t2

                if (isens .eq. 2) then
                  sumlg1 = sumlg1+gmatl(l1,k1)*t1+gmatl(l1,k2)*t2
                  sumlg2 = sumlg2+gmatl(l2,k1)*t1+gmatl(l2,k2)*t2
                  t1l = c1*rnl(in,1)
                  t2l = c2*rnl(in,2)
                  sumlg3 = sumlg3+gmat(l1,k1)*t1l+gmat(l1,k2)*t2l
                  sumlg4 = sumlg4+gmat(l2,k1)*t1l+gmat(l2,k2)*t2l
                end if

              end do
            end do

            k = 0
            do nf = 1,ntnode
              k = k+1
              k2 = k*2
              k1 = k2-1
              r1 = cord(nf,1)
              sumf1 = sumf1+fmat(l1,k1)*r1
              sumf2 = sumf2+fmat(l2,k1)*r1
              sumlf1 = sumlf1+fmatl(l1,k1)*r1
              sumlf2 = sumlf2+fmatl(l2,k1)*r1
            end do

            if(abs(ri1).lt.toli)ri1 = toli
            fmat(l1,l1) =  (sumg1-sumf1)/ri1
            fmat(l2,l1) =  (sumg2-sumf2)/ri1

            if (isens .eq. 2) then

              k = 0
              sumlf3 = 0.0D+00
              sumlf4 = 0.0D+00

              do nf = 1,ntnode
                k = k + 1
                k2 = k*2
                k1 = k2 - 1
                rl1 = cordl(nf,1)
                sumlf3 = sumlf3 + fmat(l1,k1)*rl1
                sumlf4 = sumlf4 + fmat(l2,k1)*rl1
              end do

              fmatl(l1,l1) =  (sumlg1+sumlg3-sumlf1-sumlf3)/ri1
              fmatl(l2,l1) =  (sumlg2+sumlg4-sumlf2-sumlf4)/ri1
            end if

          end do
        end do
      end if

      do i = 1,ntelem
        do j = 1,3
          nf = icon(i,j)
          r = cord(nf,1)
          if (r.eq.0.0D+00) then
            fmat(2*nf-1,2*nf-1) = 1.0D+00
            fmatl(2*nf-1,2*nf-1) = 1.0D+00
            gmat(2*nf-1,2*nf-1) = 0.0D+00
            gmatl(2*nf-1,2*nf-1) = 0.0D+00
          end if
        end do
      end do

 708  format(1x,8(1p1e13.4,1x),/) 
      return
      end 
      subroutine intstr

c*********************************************************************72
c
cc INTSTR calculates interior stresses and displacements.
c
c  Modified:
c
c    14 December 2007
c
c  Author:
c
c    Ragevendra Aithal,
c    Jeff Borggaard.
c
      implicit real*8 (a-h,o-z)

      dimension sgp(4,12),sl(4),bk(2,6),ak(2,6)
      dimension xy(3,2),sig(4),t(6),u(6),ekij(4,6),tkij(4,6),uint(2)

      common /temp4/ tt(600),uu(600),tp(600),up(600) 
      common /maprop/ emod,pr,rho,shmod,iplstr,iprob
      common /trnode/in,rn(300,2)
      common /geomat/cord(200,2),ntnode,ntelem,icon(100,3)
      common /parcns/ c1,c2,c3,c4,c5,c6,c7,c8,c9,c10
      common /body/ omega,accn
      common /inter / nip,neli(200),xip(200),yip(200),nbp,nbn(200)
      common/minlth/ rmaxl 

      if(iplstr.eq.2)then
        ks = 4
      else
        ks = 3
      end if

      l = 0

      if(nip.eq. 0 )go to 101

  99  continue

      l = l + 1 
      if(l.eq.1)then
      write ( *,900)
      if(iplstr.eq.2)then
      write ( *,920)
      else
      write ( *,910)
      end if
      end if 
      xpt = xip(l)
      ypt = yip(l)
      x1 = xpt
      x2 = ypt

      do i = 1,4 
        sig(i) = 0.0D+00
      end do

      do i = 1,2
        uint(i) = 0.0D+00
      end do

      toli = .00001D+00
      kount = 0

      do j = 1,ntelem 

        ising = 0
 
        do k = 1,ks 
          do km = 1,6
            tkij(k,km) = 0.0D+00
            ekij(k,km) = 0.0D+00
          end do
        end do

        do k = 1,2
          do km = 1,6
            bk(k,km) = 0.0D+00
            ak(k,km) = 0.0D+00
          end do
        end do

        node1 = icon(j,1) 
        node2 = icon(j,2) 
        node3 = icon(j,3) 

        do k = 1,2 
          xy(1,k) = cord(node1,k)
          xy(2,k) = cord(node2,k)
          xy(3,k) = cord(node3,k)
        end do
c
c  Check for singularity.
c
        r1 = sqrt((xy(1,1)-xpt)**2 + (xy(1,2)-ypt)**2 ) 
        r2 = sqrt((xy(2,1)-xpt)**2 + (xy(2,2)-ypt)**2 ) 
        r3 = sqrt((xy(3,1)-xpt)**2 + (xy(3,2)-ypt)**2 ) 

        if (r1 .lt. toli) then
          ising = 1 
        else if (r2 .lt. toli) then
          ising = 2 
        else if (r3 .lt. toli) then
          ising = 3 
        end if
c
c  order of gauss points
c
        call autint(ngp,xpt,ypt,xy,nsub,rmaxl) 

        if(ising.gt. 0 )then
          write ( *,901)xpt,ypt
          go to 90
        end if

        call norset ( r1, r3, ngp, nsub, sgp, sl )

        if(iplstr.eq.2)then
          call axkrnl(xpt,ypt,xy,bk,ak,ngp,nsub,sgp,
     &    sl,2,0,tkij,ekij)
        else
          call plkrnl(xpt,ypt,xy,bk,ak,ngp,nsub,sgp,
     &    sl,ising)
          call strker(xpt,ypt,xy,tkij,ekij,ngp,nsub) 
        end if

        k1 = 2*node1 - 2
        k2 = 2*node2 - 2
        k3 = 2*node3 - 2
        l1 = 6*(j-1) 
        l2 = l1 + 2
        l3 = l2 + 2 
        t(1) = tt(l1+1)-tp(l1+1) 
        t(2) = tt(l1+2)-tp(l1+2) 
        t(3) = tt(l2+1)-tp(l2+1) 
        t(4) = tt(l2+2)-tp(l2+2) 
        t(5) = tt(l3+1)-tp(l3+1) 
        t(6) = tt(l3+2)-tp(l3+2) 
        u(1) = uu(k1+1)-up(k1+1) 
        u(2) = uu(k1+2)-up(k1+2) 
        u(3) = uu(k2+1)-up(k2+1) 
        u(4) = uu(k2+2)-up(k2+2) 
        u(5) = uu(k3+1)-up(k3+1) 
        u(6) = uu(k3+2)-up(k3+2) 

        do k = 1, ks 
          do i = 1, 6 
            sig(k) = sig(k)+tkij(k,i)*t(i)-ekij(k,i)*u(i) 
          end do
        end do

        do k = 1, 2
          do i = 1, 6
            uint(k) = uint(k)+bk(k,i)*t(i)-ak(k,i)*u(i)
          end do
        end do

      end do

  103 continue

      if ( iprob .eq. 0 ) then
        go to 198
      end if

      x1x1 = x1*x1
      x2x2 = x2*x2
      x1x2 = x1*x2
      xnxn =  x1x1+x2x2
c
c   Add particular displacements if necessary
c
      if (iplstr.eq.2) then
c
c   Axisymmetric problem
c
        if ( iprob .eq. 1 .and. kount .eq. 0 ) then
c
c   Gravity force particular integrals
c
           uint(1) = uint(1)-rho*accn*pr*x1x2/emod
           uint(2) = uint(2)+0.5D+00*rho*accn*(x2x2+ pr*x1x1)/emod

        else if ( iprob .eq. 2 .and. kount .eq. 0 ) then
c
c  Centrifugal force particular integrals
c
           uint(1) = uint(1)+c1*(0.5D+00*(2.0D+00+pr)*x1x1 
     &       + (1.0D+00 - 2.0D+00*pr)*x2x2)*x1
           uint(2) = uint(2)-c1*x1x1*x2

        end if
        else
c
c  Plane problem
c
        if(iprob.eq.1.and.kount.eq. 0 )then
c
c  Gravity force particular integrals  
c
          uint(1) = uint(1) -0.5D+00* rho*accn*x1*x2*pr/shmod
          uint(2) = uint(2) 
     &      +0.25D+00*rho*accn*((1.0D+00-pr)*x2*x2 +pr*x1*x1)/shmod

        elseif (iprob.eq.2.and.kount.eq. 0 )then
c
c  Centrifugal force particular integrals
c
          uint(1) = uint(1)+c1* xnxn *x1
          uint(2) = uint(2)+c1* xnxn *x2
        end if 
        end if
c
c  Add particular stresses
c
       if(iplstr.eq.2)then
c
c  Axisymmetric problem
c
         if(iprob.eq.1 .and. kount.eq. 0 )then
c
c  Gravity force particular integrals
c
           sig(2) = sig(2)+rho*accn*x2
         elseif(iprob.eq.2.and.kount.eq. 0 )then
c
c  Centrifugal force particular integrals
c
           sig(1) = sig(1)+c4*(c5*x1x1+c6*x2x2)
           sig(2) = sig(2)-c4*(c8*x1x1-c9*x2x2)
           sig(3) = sig(3)+c4*(-c9*x1x2)
           sig(4) = sig(4)+c4*(c10*x1x1+c6*x2x2)

         end if

       else
c
c  Plane problem
c
       if(iprob.eq.1 .and. kount.eq. 0 )then
c
c  Gravity force particular integral
c
        sig(2) = sig(2)+rho*accn*x2
      elseif(iprob.eq.2.and.kount.eq. 0 )then
c
c  Centrifugal force particular integrals.
c
        sig(1) = sig(1) + c1 *(c2*xnxn + c3*x1*x1)
        sig(2) = sig(2) + c1 *(c2*xnxn + c3*x2*x2)
        sig(3) = sig(3) + c1 * c3*x1*x2
      end if
      end if

 198  continue

      x = (sig(1)-sig(2))/2.0D+00
      tmax = x**2+sig(3)**2
      tmax = 2.0D+00*dsqrt(tmax)
      x = sig(1)+sig(2)
      s1 = (x+tmax)/2.0D+00
      s2 = (x-tmax)/2.0D+00

      if(kount.eq. 0 )then
        write ( *,100)
     &  xpt,ypt,(sig(i),i = 1,ks),s1,s2,tmax,(uint(k),k = 1,2)
      else
        write ( *,110)x1,x2,(sig(i),i = 1,ks),s1,s2,tmax
      end if

  90  continue

      if(l.lt.nip)go to 99 

  101 continue

      if ( nbp .eq. 0 ) then
        return
      end if

       do kk = 1,4
         sig(kk) = 0.0D+00
       end do  

       kount = kount + 1

      if(kount.gt.nbp) then
        return
      end if

       if(kount.eq.1)then 
         write ( *,905)
         if(iplstr.ne.2)then
           write ( *,930)
         else
           write ( *,940)
        end if
      end if

       n1 = neli(kount)
       n2 = nbn(kount)
       node = icon(n1,n2)
       x1 = cord(node,1)
       x2 = cord(node,2)
       call bndstr(kount,sig)

       if(kount.le.nbp)go to 103

      return

  100 format(1x,2(e10.3,1x),9(g11.4,1x))
  110 format(1x,2(e10.3,1x),7(g11.4,1x))
  901 format(' error ! ','(',g11.4,','g11.4')','too close to the boundar
     *y')
  900 format(//1x,'internal results :'/)
  905 format(//1x,'boundary stresses :'/)
  910 format(1x,  'x-cord',4x,'y-cord',5x,'sig-x',7x, 
     * 'sig-y',7x,'sig-xy',7x,'s1',8x,'s2',12x,'tmax',8x,'ux',8x,'uy'/) 
  920 format(1x,  'r-cord',4x,'z-cord',5x,'sig-r',7x, 
     * 'sig-z',7x,'sig-rz',7x,'sig-t',8x,'s1',8x,'s2',10x,'tmax',8x,'ur'
     *,9x,'uz'/) 
  930 format(1x,  'x-cord',4x,'y-cord',5x,'sig-x',7x, 
     * 'sig-y',7x,'sig-xy',7x,'s1',8x,'s2',12x,'tmax'/) 
  940 format(1x,  'r-cord',4x,'z-cord',5x,'sig-r',7x, 
     * 'sig-z',7x,'sig-rz',7x,'sig-t',8x,'s1',8x,'s2',10x,'tmax'/)
      end 
      subroutine matcon 

c*********************************************************************72
c
cc MATCON calculates material constants for displacement and stress kernels.
c
c     refer to page nos. 80-81 of text.
c
c  Modified:
c
c    14 December 2007
c
c  Author:
c
c    Ragevendra Aithal,
c    Jeff Borggaard.
c
      implicit real*8 (a-h,o-z) 

      common /maprop/ emod, pr , rho ,shmod,iplstr,iprob  
      common /intcon/ cons1, cons2, cons3, cons4 

      pi = 3.141592653589793D+00
c
c  Plane stress.
c
      if  ( iplstr .eq. 1 ) then  
        emod = emod *(1.0D+00+2.0D+00*pr)/ (1.0D+00 + pr)**2 
        pr = pr / (1.0D+00 + pr) 
      end if

      shmod = 0.5D+00 * emod / ( 1.0D+00 + pr ) 
      cons1 = -1.0D+00 / ( 8.0D+00 * shmod * pi * ( 1.0D+00 - pr ) )
      cons2 = 3.0D+00 - 4.0D+00 * pr
      cons3 = -1.0D+00 / ( 4.0D+00 * pi * ( 1.0D+00 - pr ) )
      cons4 = 1.0D+00 - 2.0D+00 * pr
 
      return
      end 
      subroutine meshgn ( )

c*********************************************************************72
c
cc MESHGN generates the mesh.
c
c  Discussion:
c
c    this program generates input data for banerjee's boundary
c    element program dbem.
c
c    the mesh is for an axisymmetric analysis of a hollow sphere
c    value - the magnitude of perturbation of inner radius to be
c    used in a finite difference sensitivity program
c
c  Modified:
c
c    14 December 2007
c
c  Author:
c
c    Ragevendra Aithal,
c    Jeff Borggaard.
c
      implicit real*8 (a-h, o-z)

      rin = 4.0D+00
      rout = 20.0D+00
      nelr = 20
      nelt = 20
      value = 0.001D+00
      press = 1000.0D+00

      open ( unit = 15, file = 'input.txt', form = 'formatted')

      nline = 1
      write(15,*) nline

      write(15,'(a)' ) 'Axisymmetric sensitivity of a hollow sphere'
      iplstr = 2
      iprob = 0
      write(15,*) iplstr, iprob
      e = 0.3D+08
      pratio = 0.3D+00
      rho = 0.0D+00
      pi = 3.141592653589793D+00
      write(15,*) e, pratio, rho
      ntelem = 2*nelt + nelr
      ntnode = 2*ntelem + 1
      write(15,*) ntnode, ntelem
      nseg = 1
      write(15,*) nseg
      write(15,*) ntnode, ntelem
      ifr = 1
      write(15,*) ifr
c
c  specifying node locations on side 1
c  located on inside of sphere
c
10    continue

      tlen = ( pi / 2.0D+00 ) / dble ( nelt )
      delt = tlen / 2.0D+00
      theta = pi / 2.0D+00
      num = 1
      rcord = 0.0D+00
      zcord = rin

      do i = 1, 2 * nelt
        write(15,*) num, rcord, zcord
        theta = theta - delt
        rcord = rin * cos ( theta )
        zcord = rin * sin ( theta )
        num = num + 1
      end do
c
c  Specifying node locations on side 2
c
      elenr = (rout - rin)/ dble ( nelr )
      delr = elenr/2.0D+00
      rcord = rin
      zcord = 0.0D+00

      do i = 1, 2*nelr
        write(15,*) num, rcord, zcord
        rcord = rcord + delr
        num = num + 1
      end do
c
c  specifying node locations on side 3
c
      theta = 0.0D+00

      do i = 1, 2*nelt
        rcord = rout*cos(theta)
        zcord = rout*sin(theta)
        write(15,*) num, rcord, zcord
        theta = theta + delt
        num = num + 1
      end do

      write(15,*) num, 0.0D+00, rout
c
c  specifying boundary conditions
c
      nbcc = 3*(nelr + nelt) + 1
      write(15,*) nbcc
      it1 = 1
      val1 = 0.0D+00
      it2 = 1
      val2 = 0.0D+00
      write(15,*) 1, 1, it1, val1, it2, val2
      it1 = 1
      val1 = 0.0D+00
      it2 = 2
      val2 = 0.0D+00

      do i = nelt+1, nelt+nelr
        do j = 1, 3
          write(15,*) i, j, it1, val1, it2, val2
        end do
      end do

      theta = 0.0D+00
      it1 = 1
      val1 = press
      it2 = 1
      val2 =  0.0D+00

      do i = nelr + nelt + 1, nelr + 2*nelt - 1
        write(15,*) i, 1, it1, val1, it2, val2
        theta = theta + delt
        val1 = press*cos(theta)
        val2 = press*sin(theta)
        write(15,*) i, 2, it1, val1, it2, val2
        theta = theta + delt
        val1 = press*cos(theta)
        val2 = press*sin(theta)
        write(15,*) i, 3, it1, val1, it2, val2
      end do

      write(15,*) 2*nelt+nelr, 1, it1, val1, it2, val2
      theta = theta + delt
      val1 = press*cos(theta)
      val2 = press*sin(theta)
      write(15,*) 2*nelt+nelr, 2, it1, val1, it2, val2
      it1 = 1
      it2 = 1
      val1 = 0.0D+00
      val2 = press
      write(15,*) 2*nelt+nelr, 3, it1, val1, it2, val2
      nip = 0
      write(15,*) nip
      nbp = 2
      icon = 1
      icon1 = 2
      inod = 3
      inod1 = 1
      write(15,*) nbp
      write(15,*) icon, inod
      write(15,*) icon1, inod1
      write(15,*) value
      rin = rin + value

      if (value .ne. 0.0D+00 ) then
        value = 0.0D+00
        go to 10
      end if

      rewind 15

      return
      end
      subroutine normal ( xy, in, rn ) 

c*********************************************************************72
c
cc NORMAL handles the calculation of direction cosines of the node normal.
c
c  Modified:
c
c    28 December 2007
c
c  Author:
c
c    Ragevendra Aithal,
c    Jeff Borggaard.
c
c  Parameters:
c
c    Input, double precision XY(3,2), coordinates of the nodes of the element.
c
c    Input/output, integer IN, the number of nodes for which a normal
c    has been computed.
c
c    Input/output, double precision RN(300,2), contains the updated
c    normal vector information.
c
      implicit none 

      double precision eta(3)
      integer i
      integer in
      double precision rn(300,2)
      double precision xy(3,2)

      save eta

      data eta / 0.0D+00, 0.5D+00, 1.0D+00 / 
 
      do i = 1, 3
        in = in + 1
        call normpt ( eta(i), xy, rn(in,1), rn(in,2) )
      end do

      return
      end 
      subroutine normpt ( eta, xy, rn1, rn2 ) 

c*********************************************************************72
c
cc NORMPT calculates the direction of the normal vector.
c
c  Modified:
c
c    28 December 2007
c
c  Author:
c
c    Ragevendra Aithal,
c    Jeff Borggaard.
c
c  Parameters:
c
c    Input, double precision ETA, the intrinsic coordinate of the point.
c
c    Input, double precision XY(3,2), coordinates of the nodes of the element.
c
c    Output, double precision RN1, RN2, the direction of the normal vector.
c
      implicit none

      double precision dshp(3)
      double precision dxi(2)
      double precision eta
      integer i
      integer j
      double precision rjac
      double precision rn1
      double precision rn2
      double precision xy(3,2)

      dshp(1) =  4.0D+00 * eta - 3.0D+00
      dshp(2) = -8.0D+00 * eta + 4.0D+00
      dshp(3) =  4.0D+00 * eta - 1.0D+00

      do i = 1, 2
        dxi(i) = 0.0D+00
        do j = 1, 3
          dxi(i) = dxi(i) + dshp(j) * xy(j,i)
        end do
      end do

      rjac = sqrt ( dxi(1) * dxi(1) + dxi(2) * dxi(2) ) 

      rn1 =   dxi(2) / rjac 
      rn2 = - dxi(1) / rjac 

      return
      end 
      subroutine norset ( rn1, rn2, ngp, nseg, sgp, sl )

c*********************************************************************72
c
cc NORSET sets up integration kernels for the nonsingular case.
c
c  Discussion:
c
c    This routine sets up integration by subsegmentation of  
c    the kernel-shape function products for the non-singular case. 
c
c  Modified:
c
c    27 December 2007
c
c  Author:
c
c    Ragevendra Aithal,
c    Jeff Borggaard.
c
c  Parameters:
c
c    Input, double precision RN1, ?
c
c    Input, double precision RN2, ?
c
c    Input, integer NGP, the number of Gauss points.
c
c    Input, integer NSEG, the number of segments.
c
c    Output, double precision SGP(NSEG,NGP), ?
c
c    Output, double precision SL(NSEG), ?
c
      implicit real*8(a-h,o-z)

      double precision sgp(4,12)
      double precision sl(4) 
      double precision stepl
      double precision toler
      double precision toli
      double precision xsold

      common/gauswp/gp(12,12),gpw(12,12),gsp(8,8),gspw(8,8) 
c 
c  Nonsingular segment scale factor.
c
      toler = 1.0D-07
      toli = 1.0D-02
      xsold = 0.0D+00 
      stepl = 1.0D+00

      if ( 1 .lt. nseg ) then

        rdif = abs ( rn1 - rn2 ) 

        if ( toler .le. rdif ) then
          val1 = 1.0D+00 / rn1
          val2 = 1.0D+00 / rn2
          valdf = val2 - val1
        end if

      end if

      segln = 1.0D+00 / dble ( nseg )

      do inm = 1, nseg 

        do i = 1, ngp 
          sgp(inm,i) = gp(ngp,i)
        end do
 
        if ( 1 .lt. nseg ) then

          if ( toler .le. rdif ) then
            vali = dble ( inm ) * segln * valdf + val1 
            rnewi = 1.0D+00 / vali
            xsnew = ( rnewi - rn1 ) / ( rn2 - rn1 ) 
            stepl = xsnew - xsold 
          else
            stepl = segln 
            xsold = dble ( inm - 1 ) * segln 
          end if

          do i = 1, ngp 
            sgp(inm,i) = xsold + stepl * gp(ngp,i)
          end do
  
          xsold = xsnew 

        end if

        sl(inm) = stepl

      end do
 
      return
      end 
      subroutine partin ( isens )

c*********************************************************************72
c
cc PARTIN computes particular integrals.
c
c  Modified:
c
c    14 December 2007
c
c  Author:
c
c    Ragevendra Aithal,
c    Jeff Borggaard.
c
      implicit real*8 (a-h,o-z)

      integer i
      integer in
      integer iprob
      integer isens
      integer j

      common / temp4 / t(600),u(600),tp(600),up(600) 
      common / temp6 / tpl(600),upl(600)
      common / geomat / cord(200,2),ntnode,ntelem,icon(100,3)
      common / geomal /cordl(200,2), valuep
      common / maprop / emod,pr,rho,shmod,iplstr,iprob
      common / trnode / in,rn(300,2) 
      common / parcns / c1,c2,c3,c4,c5,c6,c7,c8,c9,c10
      common / body / omega,accn  

      c1 =  -rho*omega*omega*(1.0D+00+pr)*(1.0D+00-2.0D+00*pr)
     &  /(8.0D+00*emod*(1.0D+00-pr)) 
      c2 =  2.0D+00*shmod*(1.0D+00+2.0D+00*pr)/(1.0D+00-2.0D+00*pr)
      c3 =  4.0D+00*shmod 
      c4 =  -rho*omega*omega/8.0D+00
      c5 =  (3.0D+00-1.5D+00*pr - pr*pr)/(1.0D+00-pr)
      c6 =  (1.0D+00-2.0D+00*pr)/(1.0D+00-pr)
      c7 =  0.5D+00*(1.0D+00-4.0D+00*pr)*c6
      c8 =  (1.0D+00-5.0D+00*pr-2.0D+00*pr*pr)/(1.0D+00-pr)
      c9 =  2.0D+00*c6*pr
      c10 =  0.5D+00*(2.0D+00*pr*pr+3.0D+00*pr+2.0D+00)/(1.0D+00-pr)

      in = 0
      do i = 1, ncolg
        tp(i) = 0.0D+00
        up(i) = 0.0D+00
      end do

      if ( iprob .eq. 0 ) then
        return
      end if

      do i = 1, ntelem
        do j = 1, 3

          in = in+1
          nf = icon(i,j) 
          ng = (i-1)*3 + j
          ng1 = 2*ng-1
          ng2 = ng1+1
          nf1 = 2*nf-1
          nf2 = nf1+1
          x1 = cord(nf,1)
          x2 = cord(nf,2)
          x1x1 = x1*x1
          x2x2 = x2*x2
          x1x2 = x1*x2
          r11 =  rn(in,1) 
          r12 =  rn(in,2) 
c
c  Axisymmetric problem.
c
          if ( iplstr .eq. 2 ) then

            z = x2
            r = x1
            xnur = r11
            xnuz = r12
            rr = x1x1
            rz = x1x2
            zz = x2x2 

            if ( isens .eq. 2 ) then
              rl = cordl(nf,1)
              zl = cordl(nf,2)
              n1 = icon(i,1)
              n2 = icon(i,2)
              n3 = icon(i,3)
              r1 = cord(n1,1)
              r2 = cord(n2,1)
              r3 = cord(n3,1)
              z1 = cord(n1,2)
              z2 = cord(n2,2)
              z3 = cord(n3,2)
              r1l = cordl(n1,1)
              r2l = cordl(n2,1)
              r3l = cordl(n3,1)
              z1l = cordl(n1,2)
              z2l = cordl(n2,2)
              z3l = cordl(n3,2)

              if (j.eq.1) then
                s = 0.0D+00
              else if (j.eq.2) then
                s = 0.5D+00
              else if (j.eq.3) then
                s = 1.0D+00
              end if

              ra = (4.0D+00*s-3.0D+00)*r1 + (-8.0D+00*s+4.0D+00)*r2 
     &          + (4.0D+00*s-1.0D+00)*r3
              za = (4.0D+00*s-3.0D+00)*z1 + (-8.0D+00*s+4.0D+00)*z2 
     &          + (4.0D+00*s-1.0D+00)*z3
              ral = (4.0D+00*s-3.0D+00)*r1l + (-8.0D+00*s+4.0D+00)*r2l 
     &          + (4.0D+00*s-1.0D+00)*r3l
              zal = (4.0D+00*s-3.0D+00)*z1l + (-8.0D+00*s+4.0D+00)*z2l 
     &          + (4.0D+00*s-1.0D+00)*z3l
              xj = dsqrt(ra*ra + za*za)
              xjl = (ra*ral + za*zal)/xj
              xnurl = -xjl*za/xj/xj + zal/xj
              xnuzl = xjl*ra/xj/xj - ral/xj
            end if

            if(iprob.eq.1) then
c
c  Gravity force particular integrals
c
              tp(ng1) = 0.0D+00
              tp(ng2) = rho*accn*z*xnuz
              up(nf1) =  -rho*accn*pr*rz/emod
              up(nf2) = 0.5D+00*rho*accn*(zz+ pr*rr)/emod

              if (isens.eq.2) then
                tpl(ng1) = 0.0D+00
                tpl(ng2) = rho*accn*(zl*xnuz + z*xnuzl)
                upl(nf1) = -rho*accn*pr*(rl*z + r*zl)/emod
                upl(nf2) = rho*accn*(z*zl + pr*r*rl)/emod
               end if

             elseif(iprob.eq.2)then
c
c  Centrifugal force particular integrals
c
              tp(ng1) = c4*((c5*rr+c6*zz)*xnur+ (-c9*rz)*xnuz)
              tp(ng2) = c4*((-c9*rz)*xnur+ (-c8*rr+c9*zz)*xnuz)
              up(nf1) = c1*(0.5*(2.0+pr)*rr + (1.0-2.0*pr)*zz)*r
              up(nf2) = -c1*rr*z

              if (isens.eq.2) then

                tpl(ng1) = c4*(2.*(c5*r*rl+c6*z*zl)*xnur
     &          +(c5*rr+c6*zz)*xnurl)
     &          + c4*(-c9*(rl*z+r*zl)*xnuz - c9*(rz)*xnuzl)

                tpl(ng2) = c4*(-c9*(rl*z+r*zl)*xnur-c9*rz*xnurl)
     &          + c4*(2.0*(-c8*r*rl+c9*z*zl)*xnuz+(-c8*rr+c9*zz)*xnuzl)

                upl(nf1) = c1*(1.5D+00*(2.0D+00+pr)*rr*rl 
     &            + (1.0D+00-2.0D+00*pr)
     &            *(2.0D+00*z*zl*r
     &            +zz*rl))

                upl(nf2) = -c1*(2.0D+00*r*rl*z + rr*zl)

              end if

            end if

          else
c
c  Plane problem
c
            if(iprob.eq.1)then

              y = x2
              x = x1
              xnux = r11
              xnuy = r12

              if (isens .eq. 2) then

                xl = cordl(nf,1)
                xl = cordl(nf,2)
                n1 = icon(i,1)
                n2 = icon(i,2)
                n3 = icon(i,3)
                x1 = cord(n1,1)
                x2 = cord(n2,1)
                x3 = cord(n3,1)
                y1 = cord(n1,2)
                y2 = cord(n2,2)
                y3 = cord(n3,2)
                x1l = cordl(n1,1)
                x2l = cordl(n2,1)
                x3l = cordl(n3,1)
                y1l = cordl(n1,2)
                y2l = cordl(n2,2)
                y3l = cordl(n3,2)

                if (j.eq.1) then
                  s = 0.0D+00
                else if (j.eq.2) then
                  s = 0.5D+00
                else if (j.eq.3) then 
                  s = 1.0D+00
                end if

                xa = (4.0*s-3.0)*x1 + (-8.*s+4.)*x2 + (4.*s-1.)*x3
                ya = (4.0*s-3.0)*y1 + (-8.*s+4.)*y2 + (4.*s-1.)*y3
                xal = (4.0*s-3.0)*x1l + (-8.*s+4.)*x2l + (4.*s-1.)*x3l
                yal = (4.0*s-3.0)*y1l + (-8.*s+4.)*y2l + (4.*s-1.)*y3l
                xj = dsqrt(xa*xa + ya*ya)
                xjl = (xa*xal + ya*yal)/xj
                xnuxl = -xjl*ya/xj/xj + yal/xj
                xnuyl = xjl*xa/xj/xj - xal/xj
              end if
c
c  Gravity force particular integrals  
c
              tp(ng1) = 0.0D+00 
              tp(ng2) =  rho *accn* y * xnuy
              up(nf1) = -0.5D+00* rho*accn*x*y*pr/shmod
              up(nf2) = 0.25D+00*rho*accn*((1.0-pr)*y*y +pr*x*x)/shmod

              if (isens .eq. 2) then
                tpl(ng1) = 0.0D+00
                tpl(ng2) = rho*accn*(yl*xnuy + y*xnuyl)
                upl(nf1) = -0.5D+00*rho*accn*pr/shmod*(xl*y + x*yl)
                upl(nf2) = 0.5D+00*rho*accn*((1.0-pr)*y*yl + pr*x*xl)/shmod
              end if

            elseif (iprob.eq.2)then
c
c  Centrifugal force particular integrals
c
              xnxn =  x1x1 + x2x2
              xjnj =  x*xnux + y*xnuy
              tp(ng1) =  c1*(c2*xnxn*xnux + c3*x*xjnj)
              tp(ng2) =  c1*(c2*xnxn*xnuy + c3*y*xjnj)
              up(nf1) =  c1* xnxn *x1
              up(nf2) =  c1* xnxn *x2

              if (isens .eq. 2) then
                xnxnl = 2.0D+00*(x*xl + y*yl)
                xjxjl = xl*xnux + x*xnuxl + yl*xnuy + y*xnuyl
                tpl(ng1) = c1*(c2*(xnxnl*xnux + xnxn*xnuxl) +
     &                 c3*(xl*xjxj + x*xjxjl))
                tpl(ng2) = c1*(c2*(xnxnl*xnuy + xnxn*xnuyl) +
     &                 c3*(yl*xjxj + y*xjxjl))
                upl(nf1) = c1*(xnxnl*x + xnxn*xl)
                upl(nf2) = c1*(xnxnl*y + xnxn*yl)
              end if

            end if 
          end if

        end do
      end do

      return
      end
      subroutine plkrnl ( xpt, ypt, xy, bk, ak, ngp, nsub, sgp, sl, 
     &  ising )

c*********************************************************************72
c
cc PLKRNL integrates BIE kernel shape function products for stress and strain.
c
c  Modified:
c
c    14 December 2007
c
c  Author:
c
c    Ragevendra Aithal,
c    Jeff Borggaard.
c
      implicit real*8(a-h,o-z)

      real*8 ak(2,6)
      real*8 bk(2,6)
      real*8 gwpj(3)
      real*8 sgp(4,12)
      real*8 sl(4) 
      real*8 xy(3,2)

      common/gauswp/gp(12,12),gpw(12,12),gsp(8,8),gspw(8,8) 
      common/intcon/cons1,cons2,cons3,cons4

      shn1(s) =  2.0D+00 *     ( s - 0.5D+00 ) * ( s - 1.0D+00 )
      shn2(s) = -4.0D+00 * s                   * ( s - 1.0D+00 ) 
      shn3(s) =  2.0D+00 * s * ( s - 0.5D+00 )

      dsh1(s) =  4.0D+00 * s - 3.0D+00
      dsh2(s) = -8.0D+00 * s + 4.0 D+00
      dsh3(s) =  4.0D+00 * s - 1.0D+00

      tol = 0.000001D+00
 
      xa = xy(1,1)
      xb = xy(2,1)
      xc = xy(3,1)
      ya = xy(1,2)
      yb = xy(2,2)
      yc = xy(3,2)
c
c  Start Gaussian integration loop.
c 
      do i =  1, nsub

        xl =  sl(i)

        do ig = 1, ngp 
c 
c  Find coordinates of Gauss point.
c 
          s = sgp(i,ig) 
          xgp = shn1(s) * xa + shn2(s) * xb + shn3(s) * xc
          ygp = shn1(s) * ya + shn2(s) * yb + shn3(s) * yc
c
c  Find 1D Jacobian (scale factor) 
c
          dxis =  dsh1(s)*xa +dsh2(s)*xb +dsh3(s)*xc
          dyis =  dsh1(s)*ya +dsh2(s)*yb +dsh3(s)*yc
          rjac = sqrt(dxis*dxis + dyis*dyis)
c
c  Coordinates of vector between Gauss and field points.
c
          xpg = xgp -xpt
          ypg = ygp -ypt
c
c  Outward normal at Gauss point.
c
          xngp =  dyis/rjac 
          yngp = -dxis/rjac 
          rpg  = sqrt(xpg*xpg +ypg*ypg) 

          if(rpg.lt.tol)then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' )'PLKRNL - Warning:'
            write ( *, '(a)' ) '  Numerical error likely.'
            write ( *, '(a)' ) 
     &        '  A point is too close to a quadrature point.'
            write ( *, '(2g14.6)' )xpt,ypt
          end if

          rpg2 = 1.0D+00/rpg/rpg
          xrp = xpg*xpg*rpg2
          yrp = ypg*ypg*rpg2
          xyrp =  xpg*ypg*rpg2
          rnxy =  xngp*xpg + yngp*ypg 
c
c  Integration of kernels 
c
          gpwj = gpw(ngp,ig)*rjac*xl
          gwpj(1) = shn1(s)*gpwj
          gwpj(2) = shn2(s)*gpwj
          gwpj(3) = shn3(s)*gpwj
          rterm = cons1*cons2*dlog(rpg) 
c 
c  Integration of bie kernel gij
c  see page no. 80 of text for details
c 
          do k = 1,3 
            k2 = 2*k
            km = k2-1 
            bk(1,km) = bk(1,km)+gwpj(k)*(rterm-cons1*xrp) 
            bk(2,k2) = bk(2,k2)+gwpj(k)*(rterm-cons1*yrp) 
            bk(1,k2) = bk(1,k2)-gwpj(k)*cons1*xyrp
            bk(2,km) = bk(1,k2) 
          end do
c  
c  Integration of bie kernel fij 
c  see page no. 81 of text for details
c 
          rnyx =  yngp*xpg -xngp*ypg

          do k = 1,3
            if (k.ne.ising) then
              k2 = 2*k
              km = k2-1 
              bterm = gwpj(k)*cons3*rpg2
              ak(1,km) = ak(1,km)+ bterm*(cons4 +2.0D+00*xrp)*rnxy
              ak(2,k2) = ak(2,k2)+ bterm*(cons4 +2.0D+00*yrp)*rnxy
              ak(2,km) = ak(2,km)+ bterm*(cons4*rnyx +2.0D+00*xyrp*rnxy)
              ak(1,k2) = ak (1,k2)+bterm*(-cons4*rnyx+2.0D+00*xyrp*rnxy) 
            end if
          end do

        end do
      end do

      return
      end 
      subroutine setup ( iplstr, nsize, ncolg, rmaxl, intype ) 

c*********************************************************************72
c
cc SETUP sets up a particular problem.
c
c  Modified:
c
c    14 December 2007
c
c  Author:
c
c    Ragevendra Aithal,
c    Jeff Borggaard.
c
c  Parameters:
c
c    Input, integer IPLSTR, chooses the problem type.
c    0, plane strain,
c    1, plane stress,
c    2, axisymmetric elasticity,
c    3, axisymmetric potential.
c
c    Output, integer NSIZE, total number of rows.
c
c    Output, integer NCOLG, total number of columns.
c
c    Output, double precision RMAXL, the maximum length of an element.
c
c    Output, integer INTYPE(NTELEM), ?
c
      implicit real*8 (a-h,o-z) 

      integer i
      integer in
      integer intype(*)
      integer iplstr
      integer it
      integer n1
      integer n2
      integer n3
      integer ncolg
      integer nsize
      integer ntelem
      integer ntnode
      double precision rmaxl
      double precision rn(300,2)
      double precision selth
      double precision xy(3,2) 

      common /geomat/ cord(200,2),ntnode,ntelem,icon(100,3)
      common /trnode/ in,rn 
c
c  For axisymmetric potential analysis, there is only one unknown per node.
c
      if ( iplstr .eq. 3 ) then
        nsize = ntnode
        ncolg = 3 * ntelem
      else
        nsize = 2 * ntnode
        ncolg = 6 * ntelem
      end if

      do i = 1, ntelem 
        intype(i) = 0
      end do

      it = 0
      in = 0

      rmaxl = 0.0D+00

      do i = 1, ntelem

        n1 = icon(i,1)
        n2 = icon(i,2)
        n3 = icon(i,3)

        do j = 1, 2 
          xy(1,j) = cord(n1,j) 
          xy(2,j) = cord(n2,j) 
          xy(3,j) = cord(n3,j) 
        end do
c
c  Compute the normal vector at the point.
c
        call normal ( xy, in, rn ) 
c
c  Keep track of the maximum element length.
c
        selth = sqrt ( ( xy(3,1) - xy(1,1) )**2 
     &               + ( xy(3,2) - xy(1,2) )**2 )

        rmaxl = max ( rmaxl, selth )

      end do

      return
      end 
      subroutine sinset ( ngp, nsub, sgp, sl, ising ) 

c*********************************************************************72
c
cc SINSET sets up integration for the singular case.
c
c  Discussion:
c
c    This routine sets up integration by subsegmentation of 
c    the kernel-shape function products for the singular case.
c
c  Modified:
c
c    14 December 2007
c
c  Author:
c
c    Ragevendra Aithal,
c    Jeff Borggaard.
c
      implicit real*8 (a-h,o-z) 

      dimension ai(4),bi(4) ,sgp(4,12),sl(4) 

      common/gauswp/gp(12,12),gpw(12,12),gsp(8,8),gspw(8,8) 

      nsub = 4
      seg1 = 0.125
c 
c  Set up segment transformation parameters
c 
      if (ising.eq.1)then

        seg2 = (1.0D+00-seg1) / dble (nsub-1) 
        bi(1) = seg1
        i1 = 2
        i2 = nsub 

        do i = i1,i2 
          bi(i) = seg2
        end do

      elseif(ising.eq.2) then

        seg2 = (1.0D+00-2.0D+00*seg1)/(nsub-2)
        i1 = 1
        i2 = nsub/2 - 1 
        do i = i1,i2 
          bi(i) = seg2
        end do
        i1 = i2 + 1 
        i2 = i1 + 1 
        do i = i1, i2 
          bi(i) = seg1
        end do
        i1 = i2 + 1 
        i2 = nsub 
        do i = i1, i2 
          bi(i) = seg2
        end do

      elseif ( ising .eq. 3 ) then

        seg2 = ( 1.0D+00 - seg1 ) / (nsub-1)
        i1  = 1 
        i2 = nsub-1 
        do i = i1, i2 
          bi(i) = seg2
        end do
        bi(nsub) = seg1 

      end if

        ai(1) = 0.0D+00
        do i = 2, nsub
          ai(i) = ai(i-1) + bi(i-1) 
        end do
c
c  Integration subsegment loop 
c
      do inm = 1, nsub 
        sl(inm) = bi(inm)
        do i = 1, ngp 
          sgp(inm,i) = ai(inm) + bi(inm) * gp(ngp,i)
        end do
      end do

      return
      end
      subroutine solver ( amat, nsize )

c*********************************************************************72
c
cc SOLVER solves the linear system.
c
c  Discussion:
c
c    this routine creates a diagonal matrix with 1's on
c    the diagonal, the solutions are in the jn column
c
c  Modified:
c
c    14 December 2007
c
c  Author:
c
c    Ragevendra Aithal,
c    Jeff Borggaard.
c
      implicit real*8 (a-h,o-z)

      integer nsize

      real*8 amat(nsize,nsize+1) 
      integer i
      integer j
      integer jn
      real*8 rmax

      jn = nsize + 1
      rmax = 0.0D+00

      do i = 1, nsize
        rmax = max ( rmax, abs ( amat(i,i) ) )
      end do

      do j = 1, nsize

        z = amat(j,j)

        if ( dabs ( z ) .lt. 1.0D-08 * rmax ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SOLVER - Fatal error!'
          write ( *, '(a,i8)' ) '  Singular matrix, column ', j
          stop
        end if

        do l = j, jn
          amat(j,l) = amat(j,l) / z
        end do

        do i = 1, nsize

          if ( i .ne. j .and. amat(i,j) .ne. 0.0D+00 ) then

            z = amat(i,j)

            do l = j, jn
              amat(i,l) = amat(i,l) - z * amat(j,l)
            end do

          end if

        end do

      end do

      return
      end
      subroutine solver_lu ( amat, nsize )

c*********************************************************************72
c
cc SOLVER_LU solves the system using LU decomposition.
c
c  Modified:
c
c    14 December 2007
c
c  Author:
c
c    Ragevendra Aithal,
c    Jeff Borggaard.
c
      implicit real*8 (a-h,o-z)

      dimension amat(nsize,nsize+1)
c
c  performing an lu decomposition of amat
c
      jn = nsize + 1
      do i = 1, nsize-1
         z = amat(i,i)
         do j = i+1, nsize
            y = amat(j,i)/z
            do k = j, nsize+1
              amat(j,k) = amat(j,k) - y*amat(i,k)
            end do
            amat(j,i) = y
         end do
      end do
c
c  solving system using backsubstitution
c
      amat(nsize,jn) = amat(nsize,jn)/amat(nsize,nsize)

      do i = nsize-1, 1, -1
         sum = amat(i,jn)
         do j = i, nsize
           sum = sum - amat(i,j)*amat(j,jn)
        end do
        amat(i,jn) = sum/amat(i,i)
      end do

      return
      end
      subroutine store ( gmat, fmat, gmatl, fmatl, nsize, ncolg ) 

c*********************************************************************72
c
cc STORE stores the F and G matrices for subsequent usage.
c
c  Discussion:
c
c    gmat and fmat are the original matrix.
c    gmatl and fmatl are used for storing the matrix for sensitivity
c    calculation.
c
c  Modified:
c
c    14 December 2007
c
c  Author:
c
c    Ragevendra Aithal,
c    Jeff Borggaard.
c
      implicit real*8 (a-h,o-z) 

      dimension gmat(nsize,ncolg),fmat(nsize,nsize) 
      dimension gmatl(nsize,ncolg),fmatl(nsize,nsize) 
      common/geomat/ cord(200,2),ntnode,ntelem,icon(100,3)
      common /temp6 / tpl(600), upl(600)
      common /temp4 / t(600),u(600),tp(600),up(600)
c
c  Loop over all the elements
c
      do i = 1, nsize

        do j = 1, ncolg
          gmatl(i,j) = gmat(i,j)
        end do
        do j = 1, nsize
          fmatl(i,j) = fmat(i,j)
        end do
      end do

      do j = 1, ncolg
        tpl(j) = tp(j)
      end do

      do j = 1, nsize
        upl(j) = up(j)
      end do

      return
      end 
      subroutine strker ( xpt, ypt, xy, uk, tk, ngp, nseg ) 

c*********************************************************************72
c
cc STRKER computes interior stress kernels.
c
c  Modified:
c
c    14 December 2007
c
c  Author:
c
c    Ragevendra Aithal,
c    Jeff Borggaard.
c
      implicit real*8 (a-h,o-z)

      dimension xy(3,2),xi(2),dxi(2),sgp(12) 
      dimension uk(4,6),tk(4,6),gwpj(3) 
      common/intcon/cons1,cons2,cons3,cons4
      common/maprop/ emod,pr,rho,shmod,iplstr,iprob 
      common/gauswp/gp(12,12),gpw(12,12),gsp(8,8),gspw(8,8) 

      shn1(s) =  2.0D+00 *     ( s - 0.5D+00 ) * ( s - 1.0D+00 ) 
      shn2(s) = -4.0D+00 * s                   * ( s - 1.0D+00 )
      shn3(s) =  2.0D+00 * s * ( s - 0.5D+00 ) 

      dsh1(s) =  4.0D+00 * s - 3.0D+00 
      dsh2(s) = -8.0D+00 * s + 4.0D+00
      dsh3(s) =  4.0D+00 * s - 1.0 D+00
 
      pi = 3.141592653589793D+00
 
      cons5 = shmod / ( 2.0D+00 * pi * ( 1.0D+00 - pr ) ) 
      cons6 = 1.0D+00 - 4.0D+00 * pr
      cons7 = 0.25D+00 * cons4 / ( shmod * ( 1.0D+00 - pr ) )
c
c  non singular segment scale factor.
c
      toler = 1.0D-06
      xsold = 0.0D+00
      stepl = 1.0D+00

      if ( 1 .lt. nseg ) then

        rn1 = sqrt ( ( xpt - xy(1,1) )**2 + ( ypt - xy(1,2) )**2 ) 
        rn2 = sqrt ( ( xpt - xy(3,1) )**2 + ( ypt - xy(3,2) )**2 ) 
        rdif  = abs ( rn1 - rn2 )

        if ( toler .le. rdif ) then
          val1 = 1.0D+00 / rn1
          val2 = 1.0D+00 / rn2
          valdf = val2 -val1
        end if

      end if

      segln = 1.0D+00 / nseg
c 
c  Segment loop
c 
      do inm = 1,nseg 

        do i = 1,ngp 
          sgp(i) = gp(ngp,i)
        end do

        if ( 1 .lt. nseg ) then
 
          if ( toler .le. rdif ) then
            vali = inm*segln*valdf + val1 
            rnewi = 1.0D+00 / vali
            xsnew = (rnewi-rn1)/(rn2-rn1) 
            stepl = xsnew - xsold 
          else
            stepl = segln 
            xsold = (inm-1)*segln 
          end if

          do i = 1,ngp 
            sgp(i) = xsold + stepl*gp(ngp,i)
          end do

          xsold = xsnew 

        end if
c 
c  start gaussian integration loop 
c 
        do ig = 1,ngp 

          eta = sgp(ig)
          xg = 0.0D+00
          yg = 0.0D+00
          xi(1) = xy(1,1)*shn1(eta)+xy(2,1)*shn2(eta)+xy(3,1)*shn3(eta) 
          xi(2) = xy(1,2)*shn1(eta)+xy(2,2)*shn2(eta)+xy(3,2)*shn3(eta) 
          dxi(1) = xy(1,1)*dsh1(eta)+xy(2,1)*dsh2(eta)+xy(3,1)*dsh3(eta)
          dxi(2) = xy(1,2)*dsh1(eta)+xy(2,2)*dsh2(eta)+xy(3,2)*dsh3(eta)
          xg = xi(1) 
          yg = xi(2) 
          y1 = xi(1) - xpt 
          y2 = xi(2) - ypt 
          r = sqrt(y1*y1 + y2*y2)
          rj = sqrt(dxi(1)*dxi(1) + dxi(2)*dxi(2)) 
          xn = dxi(2)/rj 
          yn = - dxi(1)/rj 
          z1 = y1/r
          z2 = y2/r
          z1s2 = 2.0D+00*z1*z1 
          z2s2 = 2.0D+00*z2*z2 
          z1s4 = 2.0D+00*z1s2
          z2s4 = 2.0D+00*z2s2
          z12 = z1*z2
          a1 = -cons3*z1/r 
          a2 = -cons3*z2/r 
          b1 = cons4 + z1s2
          b2 = cons4 + z2s2
          b3 = -cons4 + z1s2 
          b4 = -cons4 + z2s2 
          gpwj = gpw(ngp,ig)*rj*stepl
          gwpj(1) = shn1(eta)*gpwj
          gwpj(2) = shn2(eta)*gpwj
          gwpj(3) = shn3(eta)*gpwj
c 
c  EJKI kernel.
c 
          do k = 1,3
            k2 = 2*k
            km = k2 - 1 
            rt1 = a1*gwpj(k)
            rt2 = a2*gwpj(k)
            uk(1,km) = uk(1,km) + rt1*b1
            uk(1,k2) = uk(1,k2) + rt2*b3
            uk(2,km) = uk(2,km) + rt1*b4
            uk(2,k2) = uk(2,k2) + rt2*b2
            uk(3,km) = uk(3,km) + rt2*b1
            uk(3,k2) = uk(3,k2) + rt1*b2
          end do

          zn = 2.0D+00*(z1*xn + z2*yn)
          znz1 = zn*z1 
          znz2 = zn*z2 
          xz12 = xn*z12
          yz12 = yn*z12
          const = cons5/(r*r)
c 
c  TJKI kernel.
c 
          do k = 1, 3
            k2 = 2*k
            k1 = k2 - 1 
            rterm = const*gwpj(k) 
            tk(1,k1) = tk(1,k1)+rterm*(znz1*(1.-z1s4)+(1.+z1s2)*xn) 
            tk(1,k2) = tk(1,k2)+rterm*(znz2*(cons4-z1s4)+4.0D+00*pr*xz12
     &                         -(cons6-cons4*z1s2)*yn ) 
            tk(2,k1) = tk(2,k1)+rterm*(znz1*(cons4-z2s4)+4.0D+00*pr*yz12
     &                         -(cons6-cons4*z2s2)*xn ) 
            tk(2,k2) = tk(2,k2)+rterm*(znz2*(1.0D+00-z2s4)
     &        +(1.0D+00+z2s2)*yn) 
            tk(3,k1) = tk(3,k1)+rterm*(znz2*(pr-z1s4)+2.0D+00*xz12*(pr+ 
     &                          cons4)+(cons4+pr*z1s2)*yn ) 
            tk(3,k2) = tk(3,k2)+rterm*(znz1*(pr-z2s4)+2.0D+00*yz12*(pr+ 
     &                          cons4)+(cons4+pr*z2s2)*xn ) 
          end do

        end do

      end do

      return 
      end
