      program main

c*********************************************************************72
c
cc MAIN is the main program for PPPACK_PRB.
c
c  Discussion:
c
c    PPPACK_PRB calls sample problems for the PPPACK library.
c
c  Modified:
c
c    15 February 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PPPACK_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the PPPACK library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
      call test05 ( )
      call test06 ( )
      call test07 ( )
      call test08 ( )
      call test09 ( )

      call test10 ( )
      call test11 ( )
      call test12 ( )
      call test13 ( )
      call test14 ( )
      call test15 ( )
      call test16 ( )
      call test17 ( )
      call test18 ( )
      call test19 ( )

      call test20 ( )
      call test21 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PPPACK_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01
c
chapter ii. runge example
c  from  * a practical guide to splines *  by c. de boor
      implicit none

      integer i,istep,j,k,n,nmk,nm1
      double precision aloger,algerp,d(20),decay,dx,errmax,g,h,
     &  pnatx,step,tau(20),x
      data step, istep /20.0D+00, 20/

      g(x) = 1.0D+00/(1.0D+00+(5.0D+00*x)**2)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) ' '

      print 600
  600 format(28h  n   max.error   decay exp.//)
      decay = 0.0D+00
      do 40 n=2,20,2
c        choose interpolation points  tau(1), ..., tau(n) , equally
c        spaced in (-1,1), and set  d(i) = g(tau(i)), i=1,...,n.
         nm1 = n-1
         h = 2.0D+00/dble(nm1)
         do 10 i=1,n
            tau(i) = dble(i-1)*h - 1.0D+00
   10       d(i) = g(tau(i))
c        calculate the divided differences for the newton form.
c
         do 20 k=1,nm1
            nmk = n-k
            do 20 i=1,nmk
   20          d(i) = (d(i+1)-d(i))/(tau(i+k)-tau(i))
c
c        estimate max.interpolation error on (-1,1).
         errmax = 0.0D+00
         do 30 i=2,n
            dx = (tau(i)-tau(i-1))/step
            do 30 j=1,istep
               x = tau(i-1) + dble(j)*dx
c              evaluate interp.pol. by nested multiplication
c
               pnatx = d(1)
               do 29 k=2,n
   29             pnatx = d(k) + (x-tau(k))*pnatx
c
   30          errmax = dmax1(errmax,dabs(g(x)-pnatx))
         aloger = dlog(errmax)
         if (n .gt. 2)  decay =
     *       (aloger - algerp)/dlog(dble(n)/dble(n-2))
         algerp = aloger
   40    print 640,n,errmax,decay
  640 format(i3,e12.4,f11.2)

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02
c
chapter iv. runge example, with cubic hermite interpolation
c  from  * a practical guide to splines *  by c. de boor
      implicit none

      integer i,istep,j,n,nm1
      double precision aloger,algerp,c(4,20),decay,divdf1,divdf3,
     &  dtau,dx,errmax,g,h,x
     &    ,pnatx,step,tau(20)
      data step, istep /20.0D+00, 20/

      g(x) = 1.0D+00/(1.0D+00+(5.0D+00*x)**2)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) ' '

      print 600
  600 format(28h  n   max.error   decay exp.//)
      decay = 0.0D+00
      do 40 n=2,20,2
c        choose interpolation points  tau(1), ..., tau(n) , equally
c        spaced in (-1,1), and set  c(1,i) = g(tau(i)), c(2,i) =
c        gprime(tau(i)) = -50.*tau(i)*g(tau(i))**2, i=1,...,n.
         nm1 = n-1
         h = 2.0D+00/dble(nm1)
         do 10 i=1,n
            tau(i) = dble(i-1)*h - 1.0D+00
            c(1,i) = g(tau(i))
   10       c(2,i) = -50.0D+00*tau(i)*c(1,i)**2
c        calculate the coefficients of the polynomial pieces
c
         do 20 i=1,nm1
            dtau = tau(i+1) - tau(i)
            divdf1 = (c(1,i+1) - c(1,i))/dtau
            divdf3 = c(2,i) + c(2,i+1) - 2.0D+00*divdf1
            c(3,i) = (divdf1 - c(2,i) - divdf3)/dtau
   20       c(4,i) = (divdf3/dtau)/dtau
c
c        estimate max.interpolation error on (-1,1).
         errmax = 0.0D+00
         do 30 i=2,n
            dx = (tau(i)-tau(i-1))/step
            do 30 j=1,istep
               h = dble(j)*dx
c              evaluate (i-1)st cubic piece
c
               pnatx = c(1,i-1)+h*(c(2,i-1)+h*(c(3,i-1)+h*c(4,i-1)))
c
   30          errmax = dmax1(errmax,dabs(g(tau(i-1)+h)-pnatx))
         aloger = dlog(errmax)
         if (n .gt. 2)  decay =
     *       (aloger - algerp)/dlog(dble(n)/dble(n-2))
         algerp = aloger
   40    print 640,n,errmax,decay
  640 format(i3,e12.4,f11.2)
      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03
c
chapter ix. example comparing the b-representation of a cubic  f  with
c  its values at knot averages.
c  from  * a practical guide to splines *  by c. de boor
c
      implicit none

      integer i,id,j,jj,n,nm4
      double precision bcoef(23),d(4),d0(4),dtip1,dtip2,
     &  f(23),t(27),tave(23),x
c            the taylor coefficients at  0  for the polynomial  f  are
      data d0 /-162.0D+00,99.0D+00,-18.0D+00,1.0D+00/

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) ' '
c
c                 set up knot sequence in the array  t .
      n = 13
      do 5 i=1,4
         t(i) = 0.0D+00
    5    t(n+i) = 10.0D+00
      nm4 = n-4
      do 6 i=1,nm4
    6    t(i+4) = dble(i)
c
      do 50 i=1,n
c        use nested multiplication to get taylor coefficients  d  at
c               t(i+2)  from those at  0 .
         do 20 j=1,4
   20       d(j) = d0(j)
         do 21 j=1,3
            id = 4
            do 21 jj=j,3
               id = id-1
   21          d(id) = d(id) + d(id+1)*t(i+2)
c
c                compute b-spline coefficients by formula (9).
         dtip1 = t(i+2) - t(i+1)
         dtip2 = t(i+3) - t(i+2)
         bcoef(i) = d(1) + (d(2)*(dtip2-dtip1)
     &     -d(3)*dtip1*dtip2)/3.0D+00
c
c                  evaluate  f  at corresp. knot average.
         tave(i) = (t(i+1) + t(i+2) + t(i+3))/3.0D+00
         x = tave(i)
   50    f(i) = d0(1) + x*(d0(2) + x*(d0(3) + x*d0(4)))
c
      print 650, (i,tave(i), f(i), bcoef(i),i=1,n)
  650 format(45h  i   tave(i)      f at tave(i)      bcoef(i)//
     *       (i3,f10.5,2f16.5))
      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04
c
chapter x. example 1. plotting some b-splines
c  from  * a practical guide to splines *  by c. de boor
calls  bsplvb, interv
      implicit none

      integer i,j,k,left,leftmk,mflag,n,npoint
      double precision dx,t(10),values(7),x,xl
c  dimension, order and knot sequence for spline space are specified...
      data n,k /7,3/
      data t /3*0.0D+00,2*1.0D+00,3.0D+00,4.0D+00,3*6.0D+00/
c  b-spline values are initialized to 0., number of evaluation points...
      data values /7*0.0D+00/, npoint /31/

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) ' '

c  set leftmost evaluation point  xl , and spacing  dx  to be used...
      xl = t(k)
      dx = (t(n+1)-t(k))/dble(npoint-1)
c
      print 600,(i,i=1,5)
  600 format('1  x',8x,5('b',i1,'(x)',7x))
c
      do 10 i=1,npoint
         x = xl + dble(i-1)*dx
c                     locate  x  with respect to knot array  t .
         call interv ( t, n+1, x, left, mflag )
         leftmk = left - k
c           get b(i,k)(x)  in  values(i) , i=1,...,n .  k  of these,
c        viz.  b(left-k+1,k)(x), ..., b(left,k)(x),  are supplied by
c        bsplvb . all others are known to be zero a priori.
         call bsplvb ( t, k, 1, x, left, values(leftmk+1) )
c
         print 610, x, (values(j),j=3,7)
  610    format(f7.3,5f12.7)
c
c           zero out the values just computed in preparation for next
c        evalulation point .
         do 10 j=1,k
   10       values(leftmk+j) = 0.0D+00
      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05
c
chapter x. example 2. plotting the pol,s which make up a b-spline
c  from  * a practical guide to splines *  by c. de boor
calls  bsplvb
c
      implicit none

      integer ia,left
      double precision biatx(4),t(11),values(4),x
c                 knot sequence set here....
      data t / 4*0.0D+00,1.0D+00,3.0D+00,4.0D+00,4*6.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) ' '

      do 20 ia=1,40
         x = dble(ia)/5.0D+00 - 1.0D+00
         do 10 left=4,7
            call bsplvb ( t, 4, 1, x, left, biatx )
c
c           according to  bsplvb  listing,  biatx(.) now contains value
c           at  x  of polynomial which agrees on the interval  (t(left),
c           t(left+1) )  with the b-spline  b(left-4 + . ,4,t) . hence,
c           biatx(8-left)  now contains value of that polynomial for
c           b(left-4 +(8-left) ,4,t) = b(4,4,t) . since this b-spline
c           has support  (t(4),t(8)), it makes sense to run  left = 4,
c           ...,7, storing the resulting values in  values(1),...,
c           values(4)  for later printing.
c
   10       values(left-3) = biatx(8-left)
   20    print 620, x, values
  620 format(f10.1,4f16.8)
      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06
c
chapter x. example 3. construction and evaluation of the pp-representat-
c                     ion of a b-spline.
c  from  * a practical guide to splines *  by c. de boor
calls bsplpp(bsplvb),ppvalu(interv)
c
      implicit none

      integer ia,l
      double precision bcoef(7),break(5),coef(4,4),scrtch(4,4),
     &  t(11),value,x,ppvalu
c                set knot sequence  t  and b-coeffs for  b(4,4,t)  ....
      data t / 4*0.0D+00,1.0D+00,3.0D+00,4.0D+00,4*6.0D+00 /
      data bcoef / 3*0.0D+00,1.0D+00,3*0.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' ) ' '

c                          construct pp-representation  ....
      call bsplpp ( t, bcoef, 7, 4, scrtch, break, coef, l )
c
c     as a check, evaluate  b(4,4,t)  from its pp-repr. on a fine mesh.
c             the values should agree with (some of) those generated in
c                        example 2 .
      do 20 ia=1,40
         x = dble(ia)/5.0D+00 - 1.0D+00
         value = ppvalu ( break, coef, l, 4, x, 0 )
   20    print 620, x, value
  620 format(f10.1,f20.8)
      return
      end
      subroutine test07 ( )

c*********************************************************************72
c
cc TEST07
c
chapter x.  example 4. construction of a b-spline via  bvalue
c  from  * a practical guide to splines *  by c. de boor
calls  bvalue(interv)
      implicit none

      integer ia
      double precision bcoef(1),bvalue,t(5),value,x
c                  set knot sequence  t  and  b-coeffs for  b(1,4,t)
      data t / 0.0D+00,1.0D+00,3.0D+00,4.0D+00,6.0D+00 /
      data bcoef / 1.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07'
      write ( *, '(a)' ) ' '

c       evaluate  b(1,4,t)  on a fine mesh.  on (0,6), the values should
c               coincide with those obtained in example 3 .
      do 20 ia=1,40
         x = dble(ia)*.2 - 1.0D+00
         value = bvalue ( t, bcoef, 1, 4, x, 0 )
   20    print 620, x, value
  620 format(f10.1,f20.8)
      return
      end
      subroutine test08 ( )

c*********************************************************************72
c
cc TEST08
c
chapter xii, example 3. cubic spline interpolation with good knots
c  from  * a practical guide to splines *  by c. de boor
calls cubspl
      implicit none

      integer i,irate,istep,j,n,nhigh,nlow,nm1,run
      double precision algerp,aloger,c(4,20),decay,dx,errmax,g,h,
     &  pnatx,step,tau(20)
     &    ,x
      data step, istep /20.0D+00, 20/
      g(x) = dsqrt(x+1.0D+00)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08'
      write ( *, '(a)' ) ' '

      do run = 1, 3

        if ( run == 1 ) then
          irate = 8
          nlow = 4
          nhigh = 10
        else if ( run == 2 ) then
          irate = 6
          nlow = 4
          nhigh = 20
        else if ( run == 3 ) then
          irate = 4
          nlow = 4
          nhigh = 20
        end if

      decay = 0.0D+00
c     read *,irate,nlow,nhigh
      print 600
  600 format(28h  n   max.error   decay exp./)
      do 40 n=nlow,nhigh,2
         nm1 = n-1
         h = 1.0D+00/dble(nm1)
         do 10 i=1,n
            tau(i) = 2.0D+00*(dble(i-1)*h)**irate - 1.0D+00
   10       c(1,i) = g(tau(i))
c        construct cubic spline interpolant.
         call cubspl ( tau, c, n, 0, 0 )
c        estimate max.interpolation error on (-1,1).
         errmax = 0.0D+00
         do 30 i=1,nm1
            dx = (tau(i+1)-tau(i))/step
            do 30 j=1,istep
               h = dble(j)*dx
               pnatx = c(1,i)+h*(c(2,i)
     &           +h*(c(3,i)+h*c(4,i)/3.0D+00)/2.0D+00)
               x = tau(i) + h
   30          errmax = dmax1(errmax,dabs(g(x)-pnatx))
         aloger = dlog(errmax)
         if (n .gt. nlow)  decay =
     &       (aloger - algerp)/dlog(dble(n)/dble(n-2))
         algerp = aloger
   40    print 640,n,errmax,decay
  640 format(i3,e12.4,f11.2)

      end do

      return
      end
      subroutine test09 ( )

c*********************************************************************72
c
cc TEST09
c
chapter xii, example 2. cubic spline interpolation with good knots
c  from  * a practical guide to splines *  by c. de boor
calls  cubspl, newnot
      implicit none

      integer i,istep,iter,itermx,j,n,nhigh,nlow,nmax,nm1,run
      parameter (nmax=20)
      double precision algerp,aloger,decay,dx,errmax,
     &  c(4,nmax),g,h,pnatx
     &    ,scrtch(2,nmax),step,tau(nmax),taunew(nmax),x
c
c     istep and step = dble(istep) specify point density for error det-
c     ermination.
      data step, istep /20.0D+00, 20/
c                               the function  g  is to be interpolated .
      g(x) = dsqrt(x+1.0D+00)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST09'
      write ( *, '(a)' ) ' '

      do run = 1, 2

        if ( run == 1 ) then
          itermx = 0
          nlow = 4
          nhigh = 20
        else
          itermx = 3
          nlow = 4
          nhigh = 20
        end if

      decay = 0.0D+00
c      read in the number of iterations to be carried out and the lower
c      and upper limit for the number  n  of data points to be used.
c     read 500,itermx,nlow,nhigh
  500 format(3i3)
      print 600, itermx
  600 format(i4,22h cycles through newnot//
     &      28h  n   max.error   decay exp./)
c                                 loop over  n = number of data points .
      do 40 n=nlow,nhigh,2
c                                        knots are initially equispaced.
         nm1 = n-1
         h = 2.0D+00/dble(nm1)
         do 10 i=1,n
   10       tau(i) = dble(i-1)*h - 1.0D+00
         iter = 1
c               construct cubic spline interpolant. then, itermx  times,
c                determine new knots from it and find a new interpolant.
   11       do 15 i=1,n
   15          c(1,i) = g(tau(i))
            call cubspl ( tau, c, n, 0, 0 )
            if (iter .gt. itermx)       go to 19
            iter = iter+1
            call newnot(tau,c,nm1,4,taunew,nm1,scrtch)
            do 18 i=1,n
   18          tau(i) = taunew(i)
                                        go to 11
   19    continue
c                            estimate max.interpolation error on (-1,1).
         errmax = 0.0D+00
c                           loop over polynomial pieces of interpolant .
         do 30 i=1,nm1
            dx = (tau(i+1)-tau(i))/step
c                interpolation error is calculated at  istep  points per
c                    polynomial piece .
            do 30 j=1,istep
               h = dble(j)*dx
               pnatx = c(1,i)+h*(c(2,i)
     &           +h*(c(3,i)+h*c(4,i)/3.0D+00)/2.0D+00)
   30       errmax = dmax1(errmax,dabs(g(tau(i)+h)-pnatx))
c                                             calculate decay exponent .
         aloger = dlog(errmax)
         if (n .gt. nlow)  decay =
     &       (aloger - algerp)/dlog(dble(n)/dble(n-2))
         algerp = aloger
c
   40    print 640,n,errmax,decay
  640 format(i3,e12.4,f11.2)

      end do

      return
      end
      subroutine test10 ( )

c*********************************************************************72
c
cc TEST10
c
chapter xii, example 4. quasi-interpolant with good knots.
c  from  * a practical guide to splines *  by c. de boor
calls bsplpp(bsplvb)
c
      implicit none

      integer i,irate,istep,j,l,m,mmk,n,nlow,nhigh,nm1,run
      double precision algerp,aloger,bcoef(22),break(20),c(4,20),
     &  decay,dg,ddg
     &    ,dtip1,dtip2,dx,errmax,g,h,pnatx,scrtch(4,4),step,t(26),
     &  taui,x

c     istep and step = dble(istep) specify point density for error det-
c     ermination.
      data step, istep /20.0D+00, 20/
c        g  is the function to be approximated,  dg  is its first, and
c        ddg  its second derivative .
      g(x) = dsqrt(x+1.0D+00)
      dg(x) = 0.5D+00/g(x)
      ddg(x) = -0.5D+00*dg(x)/(x+1.0D+00)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST10'
      write ( *, '(a)' ) ' '

      do run = 1, 2

        if ( run == 1 ) then
          irate = 8
          nlow = 4
          nhigh = 10
        else if ( run == 2 ) then
          irate = 6
          nlow = 4
          nhigh = 20
        end if

      decay = 0.0D+00
c      read in the exponent  irate  for the knot distribution and the
c      lower and upper limit for the number  n  .
c     read 500,irate,nlow,nhigh
  500 format(3i3)
      print 600
  600 format(28h  n   max.error   decay exp./)
c                              loop over  n = dim( spline(4,t) ) - 2 .
c              n  is chosen as the parameter in order to afford compar-
c              ison with examples 2 and 3  in which cubic spline interp-
c              olation at  n  data points was used .
      do 40 n=nlow,nhigh,2
         nm1 = n-1
         h = 1.0D+00/dble(nm1)
         m = n+2
         mmk = m-4
         do 5 i=1,4
            t(i) = -1.0D+00
    5       t(m+i) = 1.0D+00
c                 interior knots are equidistributed with respect to the
c                 function  (x + 1)**(1/irate) .
         do 6 i=1,mmk
    6       t(i+4) = 2.0D+00*(dble(i)*h)**irate - 1.0D+00
c                                           construct quasi-interpolant.
c                                                bcoef(1) = g(-1.) = 0.
         bcoef(1) = 0.0D+00
         dtip2 = t(5) - t(4)
         taui = t(5)
c                           special choice of  tau(2)  to avoid infinite
c                           derivatives of  g  at left endpoint .
         bcoef(2) = g(taui) - 2.0D+00*dtip2*dg(taui)/3.0D+00
     &             + dtip2**2*ddg(taui)/6.0D+00
         do 15 i=3,m
            taui = t(i+2)
            dtip1 = dtip2
            dtip2 = t(i+3) - t(i+2)
c                                      formula xii(30) of text is used .
   15       bcoef(i) = g(taui) + (dtip2-dtip1)*dg(taui)/3.0D+00
     &               - dtip1*dtip2*ddg(taui)/6.0D+00
c                                         convert to pp-representation .
         call bsplpp(t,bcoef,m,4,scrtch,break,c,l)
c                            estimate max.interpolation error on (-1,1).
         errmax = 0.0D+00
c                                             loop over cubic pieces ...
         do 30 i=1,l
            dx = (break(i+1)-break(i))/step
c                       error is calculated at  istep  points per piece.
            do 30 j=1,istep
               h = dble(j)*dx
               pnatx = c(1,i)+h*(c(2,i)
     &           +h*(c(3,i)+h*c(4,i)/3.0D+00)/2.0D+00)
   30          errmax = dmax1(errmax,dabs(g(break(i)+h)-pnatx))
c                                             calculate decay exponent .
         aloger = dlog(errmax)
         if (n .gt. nlow)  decay =
     &       (aloger - algerp)/dlog(dble(n)/dble(n-2))
         algerp = aloger
c
   40    print 640,n,errmax,decay
  640 format(i3,e12.4,f11.2)

      end do

      return
      end
      subroutine test11 ( )

c*********************************************************************72
c
cc TEST11
c
chapter xiii, example 1. a large norm amplifies noise.
c  from  * a practical guide to splines *  by c. de boor
calls splint(banfac/slv,bsplvb),bsplpp(bsplvb*),ppvalu(interv),round
c     an initially uniform data point distribution of  n  points is
c  changed  i t e r m x  times by moving the  j c l o s e -th data point
c  toward its left neighbor, cutting the distance between the two by a
c  factor of  r a t e  each time . to demonstrate the corresponding
c  increase in the norm of cubic spline interpolation at these data
c  points, the data are taken from a cubic polynomial (which would be
c  interpolated exactly) but with noise of size  s i z e  added. the re-
c  sulting noise or error in the interpolant (compared to the cubic)
c  gives the norm or noise amplification factor and is printed out
c  together with the diminishing distance  h  between the two data
c  points.
      implicit none

      integer i,iflag,istep,iter,itermx,j,jclose,l,n,nmax,nm1,run
      parameter (nmax=200)
      double precision amax,bcoef(nmax),break(nmax),coef(nmax*4),
     &  dx,fltnm1,fx,g
     & ,gtau(nmax),h,rate,scrtch(nmax*7),size,step,t(nmax+4),tau(nmax),x
      double precision ppvalu
      double precision round

      common /rount/ size
      data step, istep / 20.0D+00, 20 /
c                                          function to be interpolated .
      g(x) = 1.0D+00+x*(1.0D+00+x*(1.0D+00+x))

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST11'
      write ( *, '(a)' ) ' '

      do run = 1, 2

        if ( run == 1 ) then
          n = 7
          itermx = 10
          jclose = 4
          size = 0.000001D+00
          rate = 2.0D+00
        else if ( run == 2 ) then
          n = 7
          itermx = 10
          jclose = 4
          size = 0.0D+00
          rate = 2.0D+00
        end if

c     read 500,n,itermx,jclose,size,rate
  500 format(3i3/e10.3/e10.3)
      print 600,size
  600 format(16h size of noise =,e8.2//
     *      25h    h           max.error)
c                                       start with uniform data points .
      nm1 = n - 1
      fltnm1 = dble(nm1)
      do 10 i=1,n
   10    tau(i) = dble(i-1)/fltnm1
c                    set up knot sequence for not-a-knot end condition .
      do 11 i=1,4
         t(i) = tau(1)
   11    t(n+i) = tau(n)
      do 12 i=5,n
   12    t(i) = tau(i-2)
c
      do 100 iter=1,itermx
         do 21 i=1,n
   21       gtau(i) = round(g(tau(i)))
         call splint ( tau, gtau, t, n, 4, scrtch, bcoef, iflag )
                                        go to (24,23),iflag
   23    print 623
  623    format(27h something wrong in splint.)
                                        return
   24    call bsplpp ( t, bcoef, n, 4, scrtch, break, coef, l )
c                                    calculate max.interpolation error .
         amax = 0.0D+00
         do 30 i=4,n
            dx = (break(i-2) - break(i-3))/step
            do 25 j=2,istep
               x = break(i-2) - dx*dble(j-1)
               fx = ppvalu(break,coef,l,4,x,0)
   25          amax = dmax1(amax,dabs(fx-g(x)))
   30       continue
         h = tau(jclose) - tau(jclose-1)
         print 630,h,amax
  630    format(e9.2,e15.3)
c                 move tau(jclose) toward its left neighbor so as to cut
c                 their distance by a factor of  rate .
         tau(jclose) = (tau(jclose)
     &     + (rate-1.0D+00)*tau(jclose-1))/rate
  100    continue

      end do

      return
      end
      subroutine test12 ( )

c*********************************************************************72
c
cc TEST12
c
chapter xiii, example 2. cubic spline interpolant at knot averages
c                        with good knots.
c  from  * a practical guide to splines *  by c. de boor
calls splint(banfac/slv,bsplvb),newnot,bsplpp(bsplvb*)
      implicit none

      integer i,iflag,istep,iter,itermx,j,
     &  l,n,nhigh,nmax,nmk,nlow,run
      parameter (nmax=20)
      double precision algerp,aloger,bcoef(nmax+2),break(nmax),
     &  decay,dx,errmax
     &    ,c(4,nmax),g,gtau(nmax),h,pnatx,scrtch(nmax*7),step
     &    ,t(nmax+6),tau(nmax),tnew(nmax),x
c
c     istep and step = dble(istep) specify point density for error
c     determination.
c
      data step, istep /20.0D+00, 20/
c                               the function  g  is to be interpolated .
      g(x) = dsqrt(x+1.0D+00)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST12'
      write ( *, '(a)' ) ' '

      do run = 1, 3

        if ( run == 1 ) then
          itermx = 0
          nlow = 4
          nhigh = 20
        else if ( run == 2 ) then
          itermx = 3
          nlow = 4
          nhigh = 20
        else if ( run == 3 ) then
          itermx = 6
          nlow = 4
          nhigh = 20
        end if

      decay = 0.0D+00
c      read in the number of iterations to be carried out and the lower
c      and upper limit for the number  n  of data points to be used.
c     read 500,itermx,nlow,nhigh
  500 format(3i3)
      print 600, itermx
  600 format(i4,22h cycles through newnot//
     &      28h  n   max.error   decay exp./)
c                                 loop over  n = number of data points .
      do 40 n=nlow,nhigh,2
    4    nmk = n-4
         h = 2.0D+00/dble(nmk+1)
         do 5 i=1,4
            t(i) = -1.0D+00
    5       t(n+i) = 1.0D+00
         if (nmk .lt. 1)                go to 10
         do 6 i=1,nmk
    6       t(i+4) = dble(i)*h - 1.0D+00
   10    iter = 1
c               construct cubic spline interpolant. then, itermx  times,
c                determine new knots from it and find a new interpolant.
   11       do 12 i=1,n
               tau(i) = (t(i+1)+t(i+2)+t(i+3))/3.0D+00
   12          gtau(i) = g(tau(i))
            call splint ( tau, gtau, t, n, 4, scrtch, bcoef, iflag )
            call bsplpp ( t, bcoef, n, 4, scrtch, break, c, l )
            if (iter .gt. itermx)       go to 19
            iter = iter + 1
            call newnot ( break, c, l, 4, tnew, l, scrtch )
   15       do 18 i=2,l
   18          t(3+i) = tnew(i)
                                        go to 11
c                      estimate max.interpolation error on  (-1,1) .
   19    errmax = 0.0D+00
c                           loop over polynomial pieces of interpolant .
         do 30 i=1,l
            dx = (break(i+1)-break(i))/step
c                interpolation error is calculated at  istep  points per
c                    polynomial piece .
            do 30 j=1,istep
               h = dble(j)*dx
               pnatx = c(1,i)+h*(c(2,i)
     &           +h*(c(3,i)+h*c(4,i)/3.0D+00)/2.0D+00)
   30       errmax = dmax1(errmax,dabs(g(break(i)+h)-pnatx))
c                                             calculate decay exponent .
         aloger = dlog(errmax)
         if (n .gt. nlow)  decay =
     &       (aloger - algerp)/dlog(dble(n)/dble(n-2))
         algerp = aloger
c
   40    print 640,n,errmax,decay
  640 format(i3,e12.4,f11.2)

      end do

      return
      end
      subroutine test13 ( )

c*********************************************************************72
c
cc TEST13
c
chapter xiii, example 2m. cubic spline interpolant at knot averages
c                        with good knots. modified around label 4.
c  from  * a practical guide to splines *  by c. de Boor
calls splint(banfac/slv,bsplvb),newnot,bsplpp(bsplvb*)
      implicit none

      integer i,iflag,istep,iter,itermx,j,n,nhigh,nmax,nmk,nlow,run
      integer l
      parameter (nmax=20)
      double precision algerp,aloger,bcoef(nmax+2),break(nmax),
     &  decay,dx,errmax,
     &    c(4,nmax),g,gtau(nmax),h,pnatx,scrtch(nmax*7),step,t(nmax+6),
     &    tau(nmax),tnew(nmax),x
c
c     istep and step = dble(istep) specify point density for error
c     determination.
c
      data step, istep /20.0D+00, 20/
c                               the function  g  is to be interpolated .
      g(x) = dsqrt(x+1.0D+00)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST13'
      write ( *, '(a)' ) ' '

      do run = 1, 3

        if ( run == 1 ) then
          itermx = 0
          nlow = 4
          nhigh = 20
        else if ( run == 2 ) then
          itermx = 1
          nlow = 4
          nhigh = 20
        else if ( run == 3 ) then
          itermx = 2
          nlow = 4
          nhigh = 20
        end if

      decay = 0.0D+00
c      read in the number of iterations to be carried out and the lower
c      and upper limit for the number  n  of data points to be used.
c     read 500,itermx,nlow,nhigh
  500 format(3i3)
      print 600, itermx
  600 format(i4,22h cycles through newnot//
     &      28h  n   max.error   decay exp./)
c                                 loop over  n = number of data points .
      do 40 n=nlow,nhigh,2
         if (n .eq. nlow)               go to 4
         call newnot ( break, c, l, 4, tnew, l+2, scrtch )
         l = l + 2
         t(5+l) = 1.0D+00
         t(6+l) = 1.0D+00
         iter = 1
                                        go to 15
    4    nmk = n-4
         h = 2.0D+00/dble(nmk+1)
         do 5 i=1,4
            t(i) = -1.0D+00
    5       t(n+i) = 1.0D+00
         if (nmk .lt. 1)                go to 10
         do 6 i=1,nmk
    6       t(i+4) = dble(i)*h - 1.0D+00
   10    iter = 1
c               construct cubic spline interpolant. then, itermx  times,
c                determine new knots from it and find a new interpolant.
   11       do 12 i=1,n
               tau(i) = (t(i+1)+t(i+2)+t(i+3))/3.0D+00
   12          gtau(i) = g(tau(i))
            call splint ( tau, gtau, t, n, 4, scrtch, bcoef, iflag )
            call bsplpp ( t, bcoef, n, 4, scrtch, break, c, l )
            if (iter .gt. itermx)       go to 19
            iter = iter + 1
            call newnot ( break, c, l, 4, tnew, l, scrtch )
   15       do 18 i=2,l
   18          t(3+i) = tnew(i)
                                        go to 11
c                      estimate max.interpolation error on  (-1,1) .
   19    errmax = 0.0D+00
c                           loop over polynomial pieces of interpolant .
         do 30 i=1,l
            dx = (break(i+1)-break(i))/step
c                interpolation error is calculated at  istep  points per
c                    polynomial piece .
            do 30 j=1,istep
               h = dble(j)*dx
               pnatx = c(1,i)+h*(c(2,i)
     &           +h*(c(3,i)+h*c(4,i)/3.0D+00)/2.0D+00)
   30       errmax = dmax1(errmax,dabs(g(break(i)+h)-pnatx))
c                                             calculate decay exponent .
         aloger = dlog(errmax)
         if (n .gt. nlow)  decay =
     &       (aloger - algerp)/dlog(dble(n)/dble(n-2))
         algerp = aloger
c
   40    print 640,n,errmax,decay
  640 format(i3,e12.4,f11.2)

      end do

      return
      end
      subroutine test14 ( )

c*********************************************************************72
c
cc TEST14
c
chapter xiii, example 3 . test of optimal spline interpolation routine
c                         on titanium heat data .
c  from  * a practical guide to splines *  by c. de boor
calls titand,splopt(bsplvb,banfac/slv),splint(*),bvalue(interv)
c
C     data n,k /12,5/
c     lenscr = (n-k)(2k+3)+5k+3  is the length of scrtch required in
c                                splopt .
      implicit none

      integer k,lenscr,n,ntitan
      parameter (n=12,ntitan=49,k=5,lenscr=(n-k)*(2*k+3)+5*k+3)
      integer i,iflag,ipick(n),ipicki,lx,nmk
      double precision a(n),bvalue,gtitan(ntitan),gtau(ntitan),
     &  scrtch(lenscr),t(n+k),tau(n)
     &    ,x(ntitan)

      data ipick /1,5,11,21,27,29,31,33,35,40,45,49/

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST14'
      write ( *, '(a)' ) ' '

      call titand ( x, gtitan, lx )
      do 10 i=1,n
         ipicki = ipick(i)
         tau(i) = x(ipicki)
   10    gtau(i) = gtitan(ipicki)
      call splopt ( tau, n, k, scrtch, t, iflag )
      if (iflag .gt. 1)                 return
      call splint ( tau, gtau, t, n, k, scrtch, a, iflag )
      if (iflag .gt. 1)                 return
      do 20 i=1,lx
         gtau(i) = bvalue ( t, a, n, k, x(i), 0 )
   20    scrtch(i) = gtitan(i) - gtau(i)
      print 620,(i,x(i),gtitan(i),gtau(i),scrtch(i),i=1,lx)
  620 format(41h  i, data point, data, interpolant, error//
     &         (i3,f8.0,f10.4,f9.4,e11.3))
      nmk = n-k
      print 621,(i,t(k+i),i=1,nmk)
  621 format(///16h optimal knots =/(i5,f15.9))
      return
      end
      subroutine test15 ( )

c*********************************************************************72
c
cc TEST15
c
chapter xiv, example 1. cubic smoothing spline
c  from  * a practical guide to splines *  by c. de boor
calls  bsplpp(bsplvb), ppvalu(interv), smooth(setupq,chol1d)
c     values from a cubic b-spline are rounded to  n d i g i t  places
c  after the decimal point, then smoothed via  s m o o t h  for
c  various values of the control parameter  s .
      implicit none

      integer i,is,j,l,lsm,ndigit,npoint,ns
      double precision a(61,4),bcoef(7),break(5),coef(4,4),
     &  coefsm(4,60),dely,dy(61),ppvalu
     &    ,s(20),scrtch(427),sfp,smooth,t(11),tenton,x(61),y(61)
      equivalence (scrtch,coefsm)
      data t /4*0.0D+00,1.0D+00,3.0D+00,4.0D+00,4*6.0D+00/
      data bcoef /3*0.0D+00,1.0D+00,3*0.0D+00/

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST15'
      write ( *, '(a)' ) ' '

      call bsplpp(t,bcoef,7,4,scrtch,break,coef,l)
      npoint = 61

      ndigit = 2
      ns = 7
      do i = 1, 7
        s(i) = 6.0D+00 * 10.0D+00**(6-i)
      end do

c     read 500,ndigit,ns,(s(i),i=1,ns)
  500 format(2i3/(e10.4))
      print 600,ndigit
  600 format(24h exact values rounded to,i2,21h digits after decimal
     &       ,7h point./)
      tenton = 10.0D+00**ndigit
      dely = 0.5D+00/tenton
      do 10 i=1,npoint
         x(i) = 0.1D+00*dble(i-1)
         y(i) = ppvalu(break,coef,l,4,x(i),0)
         y(i) = dble(int(y(i)*tenton + 0.5D+00))/tenton
   10    dy(i) = dely
      do 15 i=1,npoint,5
         do 15 j=1,4
   15       a(i,j) = ppvalu(break,coef,l,4,x(i),j-1)
      print 615,(i,(a(i,j),j=1,4),i=1,npoint,5)
  615 format(52h value and derivatives of noisefree function at some
     &      ,7h points/(i4,4e15.7))
      do 20 is=1,ns
         sfp = smooth ( x, y, dy, npoint, s(is), scrtch, a )
         lsm = npoint - 1
         do 16 i=1,lsm
            do 16 j=1,4
   16          coefsm(j,i) = a(i,j)
         do 18 i=1,npoint,5
            do 18 j=1,4
   18          a(i,j) = ppvalu(x,coefsm,lsm,4,x(i),j-1)
   20    print 620, s(is), sfp, (i,(a(i,j),j=1,4),i=1,npoint,5)
  620 format(15h prescribed s =,e10.3,23h, s(smoothing spline) =,e10.3/
     &       54h value and derivatives of smoothing spline at corresp.
     &      ,7h points/(i4,4e15.7))
      return
      end
      subroutine test16 ( )

c*********************************************************************72
c
cc TEST16
c
c  main program for example in chapter xv.
c  from  * a practical guide to splines *  by c. de boor
c  solution of a second order nonlinear two point boundary value
c  problem  on (0., 1.) , by collocation with pp functions having
c  pieces of order  6 .  2  passes through newnot are to be made,
c  without any knots being added. newton iteration is to be stopped
c  when two iterates agree to  6  decimal places.
      implicit none

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST16'
      write ( *, '(a)' ) ' '

      call colloc ( 0.0D+00, 1.0D+00, 4, 6, 2, 0.0D+00, 1.0D-6 )

      return
      end
      subroutine test17 ( )

c*********************************************************************72
c
cc TEST17
c
chapter xvi, example 2, taut spline interpolation to the
c  titanium heat data of example xiii.3.
c  from  * a practical guide to splines *  by c. de boor
calls titand,tautsp,ppvalu(interv)
      implicit none

      integer i,iflag,ipick(12),ipicki,k,l,lx,n,npoint
      double precision break(122),coef(4,22),gamma,gtau(49),
     &  gtitan(49),plotf(201)
     &    ,plott(201),plotts(201),ppvalu,
     &  scrtch(119),step,tau(12),x(49)
      data n,ipick /12,1,5,11,21,27,29,31,33,35,40,45,49/
      data k,npoint /4,201/

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST17'
      write ( *, '(a)' ) ' '

      call titand ( x, gtitan, lx )
      do 10 i=1,n
         ipicki = ipick(i)
         tau(i) = x(ipicki)
   10    gtau(i) = gtitan(ipicki)
      call tautsp(tau,gtau,n,0.0D+00,scrtch,break,coef,l,k,iflag)
      if (iflag .gt. 1)                 return
      step = (tau(n) - tau(1))/dble(npoint-1)
      do 20 i=1,npoint
         plott(i) = tau(1) + step*dble(i-1)
   20    plotf(i) = ppvalu(break,coef,l,k,plott(i),0)
    1 print 601
  601 format(18h gamma = ? (f10.3))
c     read 500,gamma
      gamma = 2.5D+00
  500 format(f10.3)
      if (gamma .lt. 0.0D+00)                return
      call tautsp(tau,gtau,n,gamma,scrtch,break,coef,l,k,iflag)
      if (iflag .gt. 1)                 return
      do 30 i=1,npoint
   30    plotts(i) = ppvalu(break,coef,l,k,plott(i),0)
      print 602,gamma,(plott(i),plotf(i),plotts(i),i=1,npoint)
  602 format(/42h cubic spline vs. taut spline with gamma =,f6.3//
     &      (f15.7,2e20.8))
c                                       go to 1
      return
      end
      subroutine test18 ( )

c*********************************************************************72
c
cc TEST18
c
chapter xvi, example 3. two parametrizations of some data
c  from  * a practical guide to splines *  by c. de boor
calls  splint(bsplvb,banfac/slv),bsplpp(bsplvb*),banslv*,ppvalu(interv)
c     parameter k=4,kpkm1=7, n=8,npk=12, npiece=6,npoint=21
c     integer i,icount,iflag,kp1,l
c     real bcoef(n),break(npiece),ds,q(n,kpkm1),s(n),scrtch(k,k)
c    *    ,ss,t(npk),x(n),xcoef(k,npiece),xx(npoint),y(n)
c    *    ,ycoef(k,npiece),yy(npoint)
      implicit none

      integer i,icount,iflag,k,kpkm1,kp1,l,n,npoint
      data k,kpkm1,n,npoint /4,7,8,21/
      double precision bcoef(8),break(6),ds,ppvalu,
     &  q(8,7),s(8),scrtch(4,4)
     &    ,ss,t(12),x(8),xcoef(4,6),xx(21),y(8)
     &    ,ycoef(4,6),yy(21)
      data x / 0.0D+00, 0.1D+00, 0.2D+00, 0.3D+00, 0.301D+00,
     &  0.4D+00, 0.5D+00, 0.6D+00/

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST18'
      write ( *, '(a)' ) ' '

c  *** compute y-component and set 'natural' parametrization
      do 1 i=1,n
         y(i) = (x(i)-0.3D+00)**2
    1    s(i) = x(i)
      print 601
  601 format(26h 'natural' parametrization/6x,1hx,11x,1hy)
      icount = 1
c  *** convert data abscissae to knots. note that second and second
c     last data abscissae are not knots.
    5 do 6 i=1,k
         t(i) = s(1)
    6    t(n+i) = s(n)
      kp1 = k+1
      do 7 i=kp1,n
    7    t(i) = s(i+2-k)
c  *** interpolate to x-component
      call splint(s,x,t,n,k,q,bcoef,iflag)
      call bsplpp(t,bcoef,n,k,scrtch,break,xcoef,l)
c  *** interpolate to y-component. since data abscissae and knots are
c     the same for both components, we only need to use backsubstitution
      do 10 i=1,n
   10    bcoef(i) = y(i)
      call banslv(q,kpkm1,n,k-1,k-1,bcoef)
      call bsplpp(t,bcoef,n,k,scrtch,break,ycoef,l)
c  *** evaluate curve at some points near the potential trouble spot,
c     the fourth and fifth data points.
      ss = s(3)
      ds = (s(6)-s(3))/dble(npoint-1)
      do 20 i=1,npoint
         xx(i) = ppvalu(break,xcoef,l,k,ss,0)
         yy(i) = ppvalu(break,ycoef,l,k,ss,0)
   20    ss = ss + ds
      print 620,(xx(i),yy(i),i=1,npoint)
  620 format(2f12.7)
      if (icount .ge. 2)                return
c  *** now repeat the whole process with uniform parametrization
      icount = icount + 1
      do 30 i=1,n
   30    s(i) = dble(i)
      print 630
  630 format(/26h 'uniform' parametrization/6x,1hx,11x,1hy)
                                        go to 5
      end
      subroutine test19 ( )

c*********************************************************************72
c
cc TEST19
c
chapter xvii, example 2. bivariate spline interpolation
c  from  * a practical guide to splines *  by c. de Boor
calls  spli2d(bsplvb,banfac/slv),interv,bvalue(bsplvb*,interv*)
      implicit none

      integer i,iflag,j,jj,kp1,kx,ky,lefty,mflag,nx,ny
      parameter (nx=7,kx=3,ny=6,ky=4)
      double precision bcoef(nx,ny),bvalue,g,taux(nx),tauy(ny),
     &  tx(nx+kx),ty(ny+ky)
     &         ,work1(nx,ny),work2(nx),work3(nx*ny),x,y

      g(x,y) = dmax1(x-3.5D+00,0.0D+00)**2
     &  + dmax1(y-3.0D+00,0.0D+00)**3

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST19'
      write ( *, '(a)' ) ' '

c *** set up data points and knots
c     in x, interpolate between knots by parabolic splines, using
c     not-a-knot end condition
      do 10 i=1,nx
   10    taux(i) = dble(i)
      do 11 i=1,kx
         tx(i) = taux(1)
   11    tx(nx+i) = taux(nx)
      kp1 = kx+1
      do 12 i=kp1,nx
   12    tx(i) = (taux(i-kx+1) + taux(i-kx+2))/2.0D+00
c     in y, interpolate at knots by cubic splines, using not-a-knot
c     end condition
      do 20 i=1,ny
   20    tauy(i) = dble(i)
      do 21 i=1,ky
         ty(i) = tauy(1)
   21    ty(ny+i) = tauy(ny)
      kp1 = ky+1
      do 22 i=kp1,ny
   22    ty(i) = tauy(i-ky+2)
c  *** generate and print out function values
      print 620,(tauy(i),i=1,ny)
  620 format(' given data'//6f12.1)
      do 32 i=1,nx
         do 31 j=1,ny
   31       bcoef(i,j) = g(taux(i),tauy(j))
   32    print 632,taux(i),(bcoef(i,j),j=1,ny)
  632 format(f5.1,6e12.5)
c
c  *** construct b-coefficients of interpolant
c
      call spli2d(taux,bcoef,tx,nx,kx,ny,work2,work3,work1,iflag)
      call spli2d(tauy,work1,ty,ny,ky,nx,work2,work3,bcoef,iflag)
c
c  *** evaluate interpolation error at mesh points and print out
      print 640,(tauy(j),j=1,ny)
  640 format(//' interpolation error'//6f12.1)
      do 45 j=1,ny
         call interv(ty,ny+1,tauy(j),lefty,mflag)
         do 45 i=1,nx
            do 41 jj=1,ky
   41          work2(jj)=bvalue(tx,bcoef(1,lefty-ky+jj),nx,kx,taux(i),0)
   45       work1(i,j) = g(taux(i),tauy(j)) -
     &                  bvalue(ty(lefty-ky+1),work2,ky,ky,tauy(j),0)
      do 46 i=1,nx
   46    print 632,taux(i),(work1(i,j),j=1,ny)
      return
      end
      subroutine test20 ( )

c*********************************************************************72
c
cc TEST20
c
chapter xvii, example 3. bivariate spline interpolation
c  followed by conversion to pp-repr and evaluation
c  from  * a practical guide to splines *  by c. de boor
calls  spli2d(bsplvb,banfac/slv),interv,bspp2d(bsplvb*),ppvalu(interv*)
      implicit none

      integer i,iflag,j,jj,kp1,kx,ky,lefty,lx,ly,mflag,nx,ny
      parameter (nx=7,kx=3,ny=6,ky=4, lx=nx+1-kx, ly=ny+1-ky)
      double precision bcoef(nx,ny),g,ppvalu,taux(nx),
     &  tauy(ny),tx(nx+kx),ty(ny+ky)
     &         ,work1(nx,ny),work2(nx),work3(nx*ny)
     &         ,breakx(lx+1),breaky(ly+1),coef(kx,lx,ky,ly)
     &         ,work4(kx,kx,ny),work5(ny,kx,lx),work6(ky,ky,nx*kx),
     & x,y
c
c     note that , with the above parameters, lx=5, ly = 3
c
      equivalence (work4,work6)
      g(x,y) = dmax1(x-3.5D+00,0.0D+00)**2
     &  + dmax1(y-3.0D+00,0.0D+00)**3

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST20'
      write ( *, '(a)' ) ' '

c *** set up data points and knots
c     in x, interpolate between knots by parabolic splines, using
c     not-a-knot end condition
      do 10 i=1,nx
   10    taux(i) = dble(i)
      do 11 i=1,kx
         tx(i) = taux(1)
   11    tx(nx+i) = taux(nx)
      kp1 = kx+1
      do 12 i=kp1,nx
   12    tx(i) = (taux(i-kx+1) + taux(i-kx+2))/2.0D+00
c     in y, interpolate at knots by cubic splines, using not-a-knot
c     end condition
      do 20 i=1,ny
   20    tauy(i) = dble(i)
      do 21 i=1,ky
         ty(i) = tauy(1)
   21    ty(ny+i) = tauy(ny)
      kp1 = ky+1
      do 22 i=kp1,ny
   22    ty(i) = tauy(i-ky+2)
c  *** generate and print out function values
      print 620,(tauy(i),i=1,ny)
  620 format(11h given data//6f12.1)
      do 32 i=1,nx
         do 31 j=1,ny
   31       bcoef(i,j) = g(taux(i),tauy(j))
   32    print 632,taux(i),(bcoef(i,j),j=1,ny)
  632 format(f5.1,6e12.5)
c
c  *** construct b-coefficients of interpolant
c
      call spli2d(taux,bcoef,tx,nx,kx,ny,work2,work3,work1,iflag)
      call spli2d(tauy,work1,ty,ny,ky,nx,work2,work3,bcoef,iflag)
c
c  *** convert to pp-representation
c
      call bspp2d(tx,bcoef,nx,kx,ny,work4,breakx,work5)
      call bspp2d(ty,work5,ny,ky,lx*kx,work6,breaky,coef)
c
c  *** evaluate interpolation error at mesh points and print out
      print 640,(tauy(j),j=1,ny)
  640 format(//' interpolation error'//6f12.1)
      do 45 j=1,ny
         call interv(breaky,ly,tauy(j),lefty,mflag)
         do 45 i=1,nx
            do 41 jj=1,ky
   41        work2(jj)=ppvalu(breakx,coef(1,1,jj,lefty),lx,kx,taux(i),0)
   45       work1(i,j) = g(taux(i),tauy(j)) -
     &      ppvalu(breaky(lefty),work2,1,ky,tauy(j),0)
      do 46 i=1,nx
   46    print 632,taux(i),(work1(i,j),j=1,ny)
      return
      end
      subroutine test21 ( )

c*********************************************************************72
c
cc TEST21
c
c  main program for least-squares approximation by splines
c  from  * a practical guide to splines *  by c. de Boor (7 may 92)
calls setdat,l2knts,l2appr(bsplvb,bchfac,bchslv),bsplpp(bsplvb*)
c     ,l2err(ppvalu(interv)),ppvalu*,newnot
c
c  the program, though ostensibly written for l2-approximation, is typ-
c  ical for programs constructing a pp approximation to a function gi-
c  ven in some sense. the subprogram  l 2 a p p r , for instance, could
c  easily be replaced by one carrying out interpolation or some other
c  form of approximation.
c
c******  i n p u t  ******
c  is expected in  s e t d a t  (quo vide), specifying both the data to
c  be approximated and the order and breakpoint sequence of the pp ap-
c  proximating function to be used. further,  s e t d a t  is expected
c  to  t e r m i n a t e  the run (for lack of further input or because
c   i c o u n t  has reached a critical value).
c     the number  n t i m e s  is read in in the main program. it speci
c  fies the number of passes through the knot improvement algorithm in
c  n e w n o t  to be made. also,  a d d b r k  is read in to specify
c  that, on the average, addbrk knots are to be added per pass through
c  newnot. for example,  addbrk = .34  would cause a knot to be added
c  every third pass (as long as  ntimes .lt. 50).
c
c******  p r i n t e d  o u t p u t  ******
c  is governed by the three print control hollerith strings
c  p r b c o  = TRUE  gives printout of b-spline coeffs. of approxim.
c  p r p c o  = TRUE  gives printout of pp repr. of approximation.
c  p r f u n  = TRUE  gives printout of approximation and error at
c                     every data point.
c  the order  k , the number of pieces  l, and the interior breakpoints
c  are always printed out as are (in l2err) the mean, mean square, and
c  maximum errors in the approximation.
c
      implicit none

      integer i,icount,ii,j,k,l,lbegin,lnew,ll,lpkmax,ltkmax,n,nt,ntau
     &        ,ntimes,ntmax,run,inot

      logical prbco
      logical prpco
      logical prfun

      parameter (lpkmax=100,ntmax=200,ltkmax=2000)
      double precision addbrk,bcoef(lpkmax),break,coef,
     &  gtau,q(ltkmax),scrtch(ntmax)
     &    ,t(ntmax),tau,totalw,weight
      common / i4data / ntau
      common / r8data / tau(ntmax),gtau(ntmax),weight(ntmax),totalw
      common /approx/ break(lpkmax),coef(ltkmax),l,k

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST21'
      write ( *, '(a)' ) ' '

      do run = 1, 3

      icount = 0

      if ( run == 1 ) then
        ntimes = 1
        addbrk = 0.0D+00
        prbco = .false.
        prpco = .false.
        prfun = .false.
        inot = 0
        call setdat1 ( icount )
      else if ( run == 2 ) then
        ntimes = 20
        addbrk = 1.0D+00
        prbco = .false.
        prpco = .false.
        prfun = .false.
        inot = 1
        call setdat2 ( icount )
      else if ( run == 3 ) then
        ntimes = 4
        addbrk = 2.0D+00
        prbco = .false.
        prpco = .false.
        prfun = .false.
        inot = 0
        call setdat3 ( icount )
      end if

c        i c o u n t  provides communication with the data-input-and-
c     termination routine  s e t d a t . it is initialized to  0  to
c     signal to setdat when it is being called for the first time. after
c     that, it is up to setdat to use icount for keeping track of the
c     passes through setdat .
c
c     information about the function to be approximated and order and
c     breakpoint sequence of the approximating pp functions is gathered
c     by a
c   1 call setdat(icount)
c
c     breakpoints are translated into knots, and the number  n  of
c     b-splines to be used is obtained by a
c
      call l2knts ( break, l, k, t, n )
c
c     the integer  n t i m e s  and the real  a d d b r k  are requested
c     as well as the print controls  p r b c o ,  p r p c o  and
c     p r f u n .  ntimes  passes  are made through the subroutine new-
c     not, with an increase of  addbrk  knots for every pass .
c     print 600
  600 format(' ntimes,addbrk , prbco,prpco,prfun =? (i3,f10.5/3a2)')
c     read 500,ntimes,addbrk,prbco,prpco,prfun
  500 format(i3,f10.5/3a2)
c
      lbegin = l
      nt = 0
c
c        the b-spline coeffs.  b c o e f  of the l2-approx. are obtain-
c        ed by a
c
   10    call l2appr ( t, n, k, q, scrtch, bcoef )

         if ( prbco ) then
           print 609, (bcoef(i),i=1,n)
         end if

  609    format(//' b-spline coefficients'/(4e20.10))
c
c        convert the b-repr. of the approximation to pp repr.
c
         call bsplpp ( t, bcoef, n, k, q, break, coef, l )
         print 610, k, l, (break(ll),ll=2,l)
  610    format(//' approximation by splines of order',i3,' on ',
     &         i3,' intervals. breakpoints -'/(4e20.10))

         if ( prpco ) then
           print 611
  611      format(/' pp-representation for approximation')
           do 12 i=1,l
              ii = (i-1)*k
   12         print 613,break(i),(coef(ii+j),j=1,k)
  613      format(f9.3,4e20.10/(11x,4e20.10))
         end if
c
c        compute and print out various error norms by a
   15    continue

         call l2err ( prfun, scrtch, q )
c
c        if newnot has been applied less than  n t i m e s  times, try
c        it again to obtain, from the current approx. a possibly improv-
c        ed sequence of breakpoints with  addbrk  more breakpoints (on
c        the average) than the current approximation has.
c           if only an increase in breakpoints is wanted, without the
c        adjustment that newnot provides, a fake newnot routine could be
c        used here which merely returns the breakpoints for  l n e w
c        equal intervals .
c
         if ( nt .lt. ntimes ) then
c        if (nt .ge. ntimes)            go to 1
         lnew = lbegin + dble(nt)*addbrk

         if ( inot .eq. 0 ) then
           call newnot (break, coef, l, k, scrtch, lnew, t )
         else if ( inot .eq. 1 ) then
           call evnnot ( break, coef, l, k, scrtch, lnew, t )
         end if

         call l2knts ( scrtch, lnew, k, t, n )
         nt = nt + 1
                                        go to 10
        end if

      end do

      return
      end
      subroutine setdat1(icount)

c*********************************************************************72
c
cc SETDAT1
c
c  from  * a practical guide to splines *  by c. de boor
c  to be called in main program  l 2 m a i n .
c     this routine is set up to provide the specific data for example 2
c     in chapter xiv. for a general purpose l2-approximation program, it
c     would have to be replaced by a subroutine reading in
c           ntau, tau(i), gtau(i), i=1,...,ntau
c     and reading in or setting
c            k, l, break(i),i=1,...,l+1,  and weight(i),i=1,...,ntau,
c           as well as  totalw = sum( weight(i) , i=1,...,ntau).
c    i c o u n t  is equal to zero when setdat is called in  l 2 m a i n
c     for the first time. after that, it is up to setdat to use icount
c     for keeping track of the passes through setdat . this is important
c     since l2main relies on setdat for  t e r m i n a t i o n .
      implicit none

      integer icount,  i,k,l,lp1,ntau,ntaum1
      double precision break,coef,gtau,step,tau,totalw,weight
      integer lpkmax
      parameter (lpkmax=100)
      integer ntmax
      parameter (ntmax=200)
      integer ltkmax
      parameter (ltkmax=2000)
      common / i4data / ntau
      common / r8data / tau(ntmax),gtau(ntmax),weight(ntmax),totalw
      common /approx/ break(lpkmax),coef(ltkmax),l,k

      if (icount .gt. 0)                stop
      icount = icount + 1
      ntau = 10
      ntaum1 = ntau-1
      do 8 i=1,ntaum1
    8    tau(i) = 1.0D+00 - 0.5D+00**(i-1)
      tau(ntau) = 1.0D+00
      do 9 i=1,ntau
    9    gtau(i) = tau(i)**2  + 1.0D+00
      do 10 i=1,ntau
   10    weight(i) = 1.0D+00
      totalw = ntau
      l = 6
      lp1 = l+1
      step = 1.0D+00/dble(l)
      k = 2
      do 11 i=1,lp1
   11    break(i) = (i-1)*step
                                        return
      end
      subroutine setdat2(icount)

c*********************************************************************72
c
cc SETDAT2
c
c  from  * a practical guide to splines *  by c. de Boor (7 may 92)
c  to be called in main program  l 2 m a i n .
c     this routine is set up to provide the specific data for example 3
c     in chapter xiv.
c
      implicit none

      integer icount,  i,k,l,ntau
      double precision break,coef,gtau,roun,step,tau,totalw,weight,x
      integer lpkmax
      parameter (lpkmax=100)
      integer ntmax
      parameter (ntmax=200)
      integer ltkmax
      parameter (ltkmax=2000)
      common / i4data / ntau
      common / r8data / tau(ntmax),gtau(ntmax),weight(ntmax),totalw
      common /approx/ break(lpkmax),coef(ltkmax),l,k

      roun(x) = dble ( int( x * 100.0D+00 ) ) /100.0D+00

      if (icount .gt. 0)                stop
      icount = icount + 1
      ntau = 65
      step = 3.0D+00/dble(ntau-1)

      do 10 i=1,ntau
         tau(i) = (i-1)*step
         gtau(i) = roun(dexp(tau(i)))
   10    weight(i) = 1.0D+00
      totalw = ntau
      l = 1
      break(1) = tau(1)
      break(2) = tau(ntau)
      k = 3
                                        return
      end
      subroutine setdat3(icount)

c*********************************************************************72
c
cc SETDAT3
c
c  from  * a practical guide to splines *  by c. de boor
c  to be called in main program  l 2 m a i n .
c       calls titand
c     this routine is set up to provide the specific data for example 4
c     in chapter xiv.
      implicit none

      integer icount,  i,k,l,n,ntau
      double precision break,brkpic(9),coef,gtau,tau,totalw,weight
      integer lpkmax
      parameter (lpkmax=100)
      integer ntmax
      parameter (ntmax=200)
      integer ltkmax
      parameter (ltkmax=2000)
      common / i4data / ntau
      common / r8data / tau(ntmax),gtau(ntmax),weight(ntmax),totalw
      common /approx/ break(lpkmax),coef(ltkmax),l,k

      data brkpic/
     &  595.0D+00, 730.985D+00, 794.414D+00, 844.476D+00, 880.06D+00,
     &  907.814D+00, 938.001D+00, 976.752D+00, 1075.0D+00/, n/9/
      if (icount .gt. 0)                stop
      icount = icount + 1
      call titand ( tau, gtau, ntau )
      do 10 i=1,ntau
   10    weight(i) = 1.0D+00
      totalw = ntau
      l = n-1
      k = 5
      do 11 i=1,n
   11    break(i) = brkpic(i)
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
