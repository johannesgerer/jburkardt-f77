      program main

c*********************************************************************72
c
cc MAIN is the main program for TCELL_MASS.
c
c  Discussion:
c
c    TCELL_MASS generates mass matrices for the TCELL problem.
c
c    The finite element model uses piecewise continuous quadratic basis 
c    functions on a mesh of 6-node triangles.
c
c  Input files:
c
c    ELENODE.DAT, created by TCELL.F, contains the number of nodes and elements,
c    the node coordinates, and the node indices that form each element.
c
c    GEN_001.TXT through GEN_008.TXT.
c
c    UP1ST.DAT, the steady state flow solution for ALPHA = 1.
c
c    UP3ST.DAT, the steady state flow solution for ALPHA = 3.
c
c  Output files:
c
c    INIT.OUT
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 January 2004
c
c  Author: 
c
c    Hyung-Chun Lee,
c    Department of Mathematics
c    Ajou University, Korea
c
      implicit double precision(a-h,o-z)

      parameter(ncvtb=  8)
      parameter(ncvtb2= ncvtb+2)

      integer nx
      parameter (nx=41)

      integer ny
      parameter (ny=41)

      integer me
      parameter (me=2*(nx-1)*(ny-1))

      integer mn
      parameter (mn=(2*nx-1)*(2*ny-1))

      integer mn2
      parameter (mn2=2*mn)

      double precision ak0(ncvtb)
      double precision area(me)
      double precision ffck(ncvtb)
      double precision ffl2(ncvtb,ncvtb)
      double precision ffnk(ncvtb)
      double precision fgra(ncvtb,ncvtb)
      double precision g(mn2,ncvtb2)
      integer indx(mn)
      integer node(me,6)
      double precision qklm(ncvtb,ncvtb,ncvtb)
      double precision rkl(ncvtb,ncvtb)
      double precision we(13)
      double precision xe(me,13)
      double precision xn(mn)
      double precision ye(me,13)
      double precision yn(mn)

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TCELL_MASS:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Compute mass matrix and other ODE'
      write ( *, '(a)' ) '  coefficients for a finite element solution'
      write ( *, '(a)' ) '  of the Navier Stokes equations in a TCELL'
      write ( *, '(a)' ) '  region.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Maximum number of elements = ', me
      write ( *, '(a,i8)' ) '  Maximum number of nodes =    ', mn

      nln = 6
      xlngth = 1.0D+00
      ylngth = 1.0D+00
c
c  Read the number of nodes and elements, the node coordinates, and the node indices
c  that form each element.
c
      call grid_read ( me, mn, xn, yn, area, node, ne, nn )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Number of nodes = ', nn
      write ( *, '(a,i6)' ) '  Number of elements = ', ne

      call quad13(we,xe,ye,xn,yn,node,ne,me)
c
c  Read in the basis vectors and steady state ALPHA=3 solution.
c
      do igen = 1, ncvtb + 1

        if (igen.eq.1) then
          open ( unit = 1, file = 'gen_001.txt', status='old')
        else if (igen.eq.2) then
          open ( unit = 1, file = 'gen_002.txt', status='old')
        else if (igen.eq.3) then
          open ( unit = 1, file = 'gen_003.txt', status='old')
        else if (igen.eq.4) then
          open ( unit = 1, file = 'gen_004.txt', status='old')
        else if (igen.eq.5) then
          open ( unit = 1, file = 'gen_005.txt', status='old')
        else if (igen.eq.6) then
          open ( unit = 1, file = 'gen_006.txt', status='old')
        else if (igen.eq.7) then
          open ( unit = 1, file = 'gen_007.txt', status='old')
        else if (igen.eq.8) then
          open ( unit = 1, file = 'gen_008.txt', status='old')
        else if (igen.eq.9) then
          open ( unit = 1, file = 'up3st.dat', status='old')
        end if

        do ic = 1, nn
          read ( 1, * ) uu,vv
          g(ic,igen) = uu
          g(ic+nn,igen) = vv
        end do
        close ( unit = 1 )
      end do
c
c  Read the steady state ALPHA = 1 solution and subtract 1/3 of the
c  ALPHA = 3 solution.
c
      open ( unit = 1, file = 'up1st.dat', status = 'old' )

      do ic = 1, nn
        read(1,*) uu,vv
        g(ic,ncvtb2)=uu-1.0D+00/3.0D+00*g(ic,ncvtb+1)
        g(ic+nn,ncvtb2)=vv-1.0D+00/3.0D+00*g(ic+nn,ncvtb+1)
      end do

      close ( unit = 1 )
c
c  Compute ? and write it out.
c
      write(*,*) 'Compute (Phi_k,v(0)-beta(0)*v_r)'

      open ( unit = 1, file = 'init.out', status = 'unknown' )

      do i = 1, ncvtb
        temp1=0.0D+00
        temp2=0.0D+00
        call assem1(g,xn,yn,xe,ye,we,area,node,temp1,temp2,
     &           ne,nln,me,mn2,nn,ncvtb,ncvtb+2,i)
        ak0(i)=temp1
      end do
      write(1,63) (ak0(i),i=1,ncvtb)
      close ( unit = 1 )
c
c  Compute (Phi_i,Phi_j) and (Phi_i',Phi_j')
c
      write(*,*) 'Compute (Phi_i,Phi_j) and (Grad(Phi_i),Grad(Phi_j))'

      do i = 1, ncvtb
        do j = 1, ncvtb
          temp1=0.0D+00
          temp2=0.0D+00
          call assem1(g,xn,yn,xe,ye,we,area,node,temp1,temp2,
     &           ne,nln,me,mn2,nn,ncvtb,i,j)
          ffl2(i,j)=temp1
          fgra(i,j)=temp2
        end do
      end do

      write(*,*) ffl2(1,1),ffl2(1,2),ffl2(1,3)
      write(*,*) ffl2(2,1),ffl2(2,2),ffl2(2,3)
      write(*,*) fgra(1,1),fgra(1,2),fgra(1,3)
      write(*,*) fgra(2,1),fgra(2,2),fgra(2,3)

c      if(ffl2(1,1).ne.1.D+00) stop

      write(*,*) 'Compute (Phi_i,v_r)'

      do i = 1, ncvtb
        temp1=0.0D+00
        temp2=0.0D+00
        call assem1(g,xn,yn,xe,ye,we,area,node,temp1,temp2,
     &           ne,nln,me,mn2,nn,ncvtb,i,ncvtb+1)
        ffnk(i)=temp1
      end do

      write(*,*) 'Compute Q_klm'

      do i = 1, ncvtb
        do j = 1, ncvtb
          do k = 1, ncvtb
               temp=0.0D+00
          call assem2(g,xn,yn,xe,ye,we,area,node,temp,
     &           ne,nln,me,mn2,nn,ncvtb,i,j,k)
            qklm(i,j,k)=temp
          end do
          write(*,*) i,j
        end do
      end do
c
c   Compute R_kl
c
      write(*,*) 'Compute R_kl first'

        kk=ncvtb+1

      do i = 1, ncvtb
        do j = 1, ncvtb
             temp=0.0D+00
           call assem2(g,xn,yn,xe,ye,we,area,node,temp,
     &           ne,nln,me,mn2,nn,ncvtb,i,j,kk)
            rkl(i,j)=temp
            write(*,*) i,j
        end do
      end do

      write(*,*) 'Compute R_kl second'

      do i = 1, ncvtb
        do j = 1, ncvtb
           temp=0.0D+00
           call assem2(g,xn,yn,xe,ye,we,area,node,temp,
     &           ne,nln,me,mn2,nn,ncvtb,i,kk,j)
           rkl(i,j)=rkl(i,j)+temp
           write(*,*) i,j
        end do
      end do
c
c   Compute C_k
c
      write(*,*) 'compute C_k'

      do i = 1, ncvtb
        temp=0.0D+00
        call assem2(g,xn,yn,xe,ye,we,area,node,temp,
     &      ne,nln,me,mn2,nn,ncvtb,i,kk,kk)
        ffck(i)=temp
      end do
c 
c  Writing result
c
      open ( unit = 1, file = 'matrix.out', status = 'unknown' )

      do i = 1, ncvtb
        write(1,63) (ffl2(i,j),j=1,ncvtb)
      end do
      do i = 1, ncvtb
        write(1,63) (fgra(i,j),j=1,ncvtb)
      end do
      do i = 1, ncvtb
        do j = 1, ncvtb
          write(1,63) (qklm(i,j,k),k=1,ncvtb)
        end do
      end do
      write(1,63) (ffnk(i), i=1,ncvtb)
      do i = 1, ncvtb
        write(1,63) (rkl(i,j),j=1,ncvtb)
      end do
      write(1,63) (ffck(i),i=1,ncvtb)

      close ( unit = 1 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TCELL_MASS:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop

61    format(4g15.5)
62    format(g15.5)
63    format(16g25.15)
      end
      subroutine assem1(g,xn,yn,xe,ye,we,area,node,temp1,temp2,
     &                 ne,nln,me,mn2,nn,ncvtb,ncvt1,ncvt2)

c*********************************************************************72
c
cc ASSEM1 determines the initial value of the coefficient vector.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 January 2004
c
c  Author: 
c
c    Hyung-Chun Lee Ph.D
c    Department of Mathematics
c    Ajou University, Korea
c
c  Parameters:
c
      implicit double precision(a-h,o-z)

      double precision area(*)
      double precision g(mn2,*)
      integer node(me,*)
      double precision we(*)
      double precision xe(me,*)
      double precision xn(*)
      double precision ye(me,*)
      double precision yn(*)

      do ie = 1, ne
        do iq = 1, 13
          x=xe(ie,iq)
          y=ye(ie,iq)
          ar=area(ie)*we(iq)
          do iln = 1, nln
            in=node(ie,iln)
            call qbf(bb,bx,by,x,y,ie,iln,xn,yn,node,me)
            do jln = 1, nln

              jn=node(ie,jln)

              call qbf(bbb,bbx,bby,x,y,ie,jln,xn,yn,node,me)

              aij=bb*g(in,ncvt1)*bbb*g(jn,ncvt2)
     &             +bb*g(in+nn,ncvt1)*bbb*g(jn+nn,ncvt2)

              bij=bx*g(in,ncvt1)*bbx*g(jn,ncvt2)
     &                   +by*g(in,ncvt1)*bby*g(jn,ncvt2)
     &             +bx*g(in+nn,ncvt1)*bbx*g(jn+nn,ncvt2)
     &                   +by*g(in+nn,ncvt1)*bby*g(jn+nn,ncvt2)

              temp1=temp1+aij*ar
                  temp2=temp2+bij*ar

            end do
          end do
        end do
      end do

      return
      end
      subroutine assem2(g,xn,yn,xe,ye,we,area,node,temp,
     &          ne,nln,me,mn2,nn,ncvtb,ncvt1,ncvt2,ncvt3)

c*********************************************************************72
c
cc ASSEM2 determines the mass matrix and other coefficients in the ODE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 January 2004
c
c  Author: 
c
c    Hyung-Chun Lee Ph.D
c    Department of Mathematics
c    Ajou University, Korea
c
c  Parameters:
c
      implicit double precision(a-h,o-z)

      double precision area(*)
      double precision g(mn2,*)
      integer node(me,*)
      double precision we(*)
      double precision xe(me,*)
      double precision xn(*)
      double precision ye(me,*)
      double precision yn(*)
c
c  Evaluate
c
      do ie = 1, ne
        do iq = 1, 13
          x=xe(ie,iq)
          y=ye(ie,iq)
          ar=area(ie)*we(iq)
          do iln = 1, nln
            in=node(ie,iln)
            call qbf(bb,bx,by,x,y,ie,iln,xn,yn,node,me)
            do jln = 1, nln
              jn=node(ie,jln)
              call qbf(bbb,bbx,bby,x,y,ie,jln,xn,yn,node,me)
              do kln = 1, nln
                    kn=node(ie,kln)
                call qbf(bbbb,bbbx,bbby,x,y,ie,kln,xn,yn,node,me)

                cij=bb*g(in,ncvt1)*bbb*g(jn,ncvt2)*bbbx*g(kn,ncvt3)
     &            +bb*g(in,ncvt1)*bbb*g(jn+nn,ncvt2)*bbby*g(kn,ncvt3)
     &            +bb*g(in+nn,ncvt1)*bbb*g(jn,ncvt2)*bbbx*g(kn+nn,ncvt3)
     &         +bb*g(in+nn,ncvt1)*bbb*g(jn+nn,ncvt2)*bbby*g(kn+nn,ncvt3)

                temp=temp+cij*ar

              end do
             end do
          end do
        end do
      end do

      return
      end
      subroutine grid_read ( me, mn, xn, yn, area, node, ne, nn)

c*********************************************************************72
c
cc GRID_READ reads the node and element information and computes element areas.
c
c  Discussion:
c
c    The routine reads a file named "ELENODE.DAT".
c
c    The first line of this file contains NN and NE, the number of nodes
c    and the number of elements.
c
c    Each of the next NN lines lists a node index, and an X and Y coordinate.
c
c    Each of the next NE lines lists an element index, and the six nodes
c    that make up that element, in a particular order.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 January 2004
c
c  Author: 
c
c    Hyung-Chun Lee Ph.D
c    Department of Mathematics
c    Ajou University, Korea
c
c  Parameters:
c
c    Input, integer ME, the maximum number of elements.
c
c    Input, integer MN, the maximum number of nodes.
c
c    Output, double precision XN(NN), YN(NN), the coordinates of the nodes.
c
c    Output, double precision AREA(NE), the area of each element.
c
c    Output, integer NODE(ME,6), the indices of the nodes that comprise each element.
c
c    Output, integer NE, the number of elements.
c
c    Output, integer NN, the number of nodes.
c
      implicit double precision (a-h,o-z)

      integer me
      integer mn

      double precision area(me)
      integer ne
      integer nn
      integer node(me,6)
      double precision xn(mn)
      double precision yn(mn)

      open ( unit = 1, file = 'elenode.dat', status = 'old' )
c
c  Read the data 
c
      read ( 1, * ) nn, ne

      do i = 1, nn
        read ( 1, * ) ii, xn(i), yn(i)
      end do

      do it = 1, ne
        read ( 1, '(7i6)' ) itt, ( node(it,i), i = 1, 6 )
      end do

      close ( unit = 1 )

      do ie = 1, ne
        i1=node(ie,1)
        i2=node(ie,2)
        i3=node(ie,3)
        x1=xn(i1)
        x2=xn(i2)
        x3=xn(i3)
        y1=yn(i1)
        y2=yn(i2)
        y3=yn(i3)
        area(ie)=0.5D+00*dabs(y1*(x2-x3)+y2*(x3-x1)+y3*(x1-x2))
      end do

      return
      end
      subroutine qbf ( bb, bx, by, x, y, ie, iln, xn, yn, node, me )

c*********************************************************************72
c
cc QBF evaluates the quadratic basic functions.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 January 2004
c
c  Author: 
c
c    Hyung-Chun Lee Ph.D
c    Department of Mathematics
c    Ajou University, Korea
c
c  Parameters:
c
c    Output, double precision BB, BX, BY, the value of the basis function,
c    and its X and Y derivatives, associated with the given element, and evaluated
c    at the given point.
c
c    Input, double precision X, Y, the coordinates of the point where the evaluation
c    is to take place.
c
c    Input, integer IE, the index of the element within which the evaluation point lies.
c
c    Input, integer ILN, the index, between 1 and 6, of the basis function.  This is also
c    the local index of the associated node in the (6 node) triangular element.
c
c    Input, double precision XN(*), YN(*), the coordinates of all nodes.
c
c    Input, integer NODE(ME,6), lists the nodes that form each element.
c
c    Input, integer ME, the maximum number of elements.
c
      implicit double precision(a-h,o-z)

      double precision bb
      double precision bx
      double precision by
      integer node(me,*)
      double precision xn(*)
      double precision yn(*)

      in1=1+mod(iln-1,3)
      in2=mod(in1,3)+1
      in3=mod(in1+1,3)+1
      i1=node(ie,in1)
      i2=node(ie,in2)
      i3=node(ie,in3)
      x1=xn(i1)
      x2=xn(i2)
      x3=xn(i3)
      y1=yn(i1)
      y2=yn(i2)
      y3=yn(i3)
      d=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)

      if ( d == 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'QBF - Fatal error!'
        write ( *, '(a)' ) '  Determinant value D = 0.'
        write ( *, '(a,i6)' ) '  Element number = ', ie
        write ( *, '(a,i6,2g14.6)' ) '  I1, X1, Y1 = ', i1, x1, y1
        write ( *, '(a,i6,2g14.6)' ) '  I2, X2, Y2 = ', i2, x2, y2
        write ( *, '(a,i6,2g14.6)' ) '  I3, X3, Y3 = ', i3, x3, y3
        stop
      end if

      t=1.0D+00+((y2-y3)*(x-x1)+(x3-x2)*(y-y1))/d

      if ( iln .le. 3 ) then
        bb=t*(2.0D+00*t-1.0D+00)
        bx=(y2-y3)*(4.0D+00*t-1.0D+00)/d
        by=(x3-x2)*(4.0D+00*t-1.0D+00)/d
      else
        j1=i2
        j2=i3
        j3=i1
        a1=xn(j1)
        a2=xn(j2)
        a3=xn(j3)
        b1=yn(j1)
        b2=yn(j2)
        b3=yn(j3)
        c=(a2-a1)*(b3-b1)-(a3-a1)*(b2-b1)
        s=1.0D+00+((b2-b3)*(x-a1)+(a3-a2)*(y-b1))/c
        bb=4.0D+00*s*t
        bx=4.0D+00*(t*(b2-b3)/c+s*(y2-y3)/d)
        by=4.0D+00*(t*(a3-a2)/c+s*(x3-x2)/d)
      end if

      return
      end
      subroutine quad13 ( we, xe, ye, xn, yn, node, ne, me )

c*********************************************************************72
c
cc QUAD13 sets the nodes and weights for 13 point Gauss quadrature.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 January 2004
c
c  Author: 
c
c    Hyung-Chun Lee Ph.D
c    Department of Mathematics
c    Ajou University, Korea
c
c  Parameters:
c
      implicit double precision (a-h,o-z)

      integer node(me,*)
      double precision we(*)
      double precision xe(me,*)
      double precision xn(*)
      double precision ye(me,*)
      double precision yn(*)

      do i = 1, 3
        we(i)=0.175615257433204D+00
        ii=i+3
        we(ii)=0.053347235608839D+00
        ii=i+6
        iii=ii+3
        we(ii)=0.077113760890257D+00
        we(iii)=we(ii)
      end do

      we(13)=-0.14957004446767D+00

      z1=0.479308067841923D+00
      z2=0.260345966079038D+00
      z3=0.869739794195568D+00
      z4=0.065130102902216D+00
      z5=0.638444188569809D+00
      z6=0.312865496004875D+00
      z7=0.048690315425316D+00

      do ie = 1, ne
        i1=node(ie,1)
        i2=node(ie,2)
        i3=node(ie,3)
        x1=xn(i1)
        x2=xn(i2)
        x3=xn(i3)
        y1=yn(i1)
        y2=yn(i2)
        y3=yn(i3)
        xe(ie,1)=z1*x1+z2*x2+z2*x3
        ye(ie,1)=z1*y1+z2*y2+z2*y3
        xe(ie,2)=z2*x1+z1*x2+z2*x3
        ye(ie,2)=z2*y1+z1*y2+z2*y3
        xe(ie,3)=z2*x1+z2*x2+z1*x3
        ye(ie,3)=z2*y1+z2*y2+z1*y3
        xe(ie,4)=z3*x1+z4*x2+z4*x3
        ye(ie,4)=z3*y1+z4*y2+z4*y3
        xe(ie,5)=z4*x1+z3*x2+z4*x3
        ye(ie,5)=z4*y1+z3*y2+z4*y3
        xe(ie,6)=z4*x1+z4*x2+z3*x3
        ye(ie,6)=z4*y1+z4*y2+z3*y3
        xe(ie,7)=z5*x1+z6*x2+z7*x3
        ye(ie,7)=z5*y1+z6*y2+z7*y3
        xe(ie,8)=z5*x1+z7*x2+z6*x3
        ye(ie,8)=z5*y1+z7*y2+z6*y3
        xe(ie,9)=z6*x1+z5*x2+z7*x3
        ye(ie,9)=z6*y1+z5*y2+z7*y3
        xe(ie,10)=z6*x1+z7*x2+z5*x3
        ye(ie,10)=z6*y1+z7*y2+z5*y3
        xe(ie,11)=z7*x1+z5*x2+z6*x3
        ye(ie,11)=z7*y1+z5*y2+z6*y3
        xe(ie,12)=z7*x1+z6*x2+z5*x3
        ye(ie,12)=z7*y1+z6*y2+z5*y3
        xe(ie,13)=(x1+x2+x3)/3.0D+00
        ye(ie,13)=(y1+y2+y3)/3.0D+00
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
