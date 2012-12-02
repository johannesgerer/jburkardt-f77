      program main

c*********************************************************************72
c
cc FISHPRB calls the FISHPACK tests.
c
c     KPRINT = 0  Quick checks - No printing.
c                 Driver       - Short pass or fail message printed.
c              1  Quick checks - No message printed for passed tests,
c                                short message printed for failed tests.
c                 Driver       - Short pass or fail message printed.
c              2  Quick checks - Print short message for passed tests,
c                                fuller information for failed tests.
c                 Driver       - Pass or fail message printed.
c              3  Quick checks - Print complete quick check results.
c                 Driver       - Pass or fail message printed.
c
c  Modified:
c
c    21 September 2010
c
      implicit none

      integer ipass
      integer kprint

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FISHPACK_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the FISHPACK library.'

      kprint=2

      call xermax ( 1000 )

      if (kprint .le. 1) then
        call xsetf(0)
      else
        call xsetf(1)
      end if
c
c  Test HWSCRT
c
      call test01 ( kprint, ipass )
c
c  Test HWSPLR
c
      call test02 ( kprint, ipass )
c
c  Test HWSCYL
c
      call test03 ( kprint, ipass )
c
c  Test HWSSSP
c
      call test04 ( kprint, ipass )
c
c  Test HWSCSP
c
      call test05 ( kprint, ipass )
 
      call test06 ( kprint, ipass )
c
c  Test GENBUN
c
      call test07 ( kprint, ipass )
c
c  Test BLKTRI.
c
      call test08 ( kprint, ipass )
c
c  Test HWSCRT.
c
      call test09 ( )
c
c  Test HWSPLR.
c
      call test10 ( )
c
c  Test BLKTRI.
c
      call test11 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FISHPACK_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( kprint, ipass )

c*********************************************************************72
c
cc TEST01 illustrates the use of HWSCRT.
c
c  Discussion:
c
c    The equation solved is:
c
c      (D/DX)(DU/DX) + (D/DY)(DU/DY) - 4*U
c
c      = (2 - (4 + PI**2/4)*X**2)*COS((Y+1)*PI/2)
c
c    on the rectangle 0 .LT. X .LT. 2, -1 .LT. Y .LT. 3 with the
c    with the boundary conditions
c
c      U(0,Y) = 0
c                                          -1 .LE. Y .LE. 3
c      (DU/DX)(2,Y) = 4*COS((Y+1)*PI/2)
c
c    and with U periodic in Y.
c
c    The X-interval will be divided into 40 intervals and the
c    Y-interval will be divided into 80 intervals.
c
c  Modified:
c
c    21 September 2010
c
      implicit none

      real a
      real b
      real bda
      real bdb(81)
      real bdc
      real bdd
      real c
      real d
      real elmbda
      real ermax
      real err
      real f(45,82)
      integer i
      integer idimf
      integer ierror
      integer ipass
      integer j
      integer kprint
      integer m
      integer mbdcnd
      integer n
      integer nbdcnd
      real pertrb
      real pi
      real piby2
      real pimach
      real pisq
      real w(1200)
      real x(41)
      real y(81)
      real z
c
c  From dimension statement we get value of IDIMF.  Also note that W
c  is dimensioned 6*(N+1) + 8*(M+1).
c
      idimf = 45
      ermax = 0.001
      a = 0.0
      b = 2.0
      m = 40
      mbdcnd = 2
      c = -1.0
      d = 3.0
      n = 80
      nbdcnd = 0
      elmbda = -4.0
c
c  Auxiliary quantities.
c
      pi = pimach()
      piby2 = pi/2.
      pisq = pi**2
c
c  Generate and store grid points for the purpose of computing
c  boundary data and the right side of the Helmholtz equation.
c
      do i=1,m+1
        x(i) = (i-1)/20.0
      end do

      do j=1,n+1
        y(j) = -1.0+(j-1)/20.0
      end do
c
c  Generate boundary data.
c
      do j=1,n+1
        bdb(j) = 4.0*cos((y(j)+1.)*piby2)
      end do
c
c  BDA, BDC, and BDD are dummy variables.
c
      do j=1,n+1
        f(1,j) = 0.0
      end do
c
c  Generate right side of equation.
c
      do i=2,m+1
        do j=1,n+1
          f(i,j) = (2.-(4.+pisq/4.)*x(i)**2)*cos((y(j)+1.)*piby2)
        end do
      end do

      call hwscrt(a,b,m,mbdcnd,bda,bdb,c,d,n,nbdcnd,bdc,bdd,elmbda,f,
     &  idimf,pertrb,ierror,w)
c
c  Compute discretization error.  The exact solution is
c    U(X,Y) = X**2*COS((Y+1)*PIBY2)
c
      err = 0.0
      do i=1,m+1
        do j=1,n+1
          z = abs(f(i,j)-x(i)**2*cos((y(j)+1.)*piby2))
          if (z .gt. err) err = z
        end do
      end do

      ipass = 1
      if (err.gt.ermax) ipass = 0
      if (kprint.eq.0) return
      if (kprint.ge.2 .or. ipass.eq.0) then
        write(*,*)' '
        write(*,1001) ierror,err,int(w(1))
        if (ipass.eq.1) then
          write(*,*)'Test passed.'
        else
          write(*, *)'Test failed.'
        end if
      end if
      return

 1001 format (' TEST01: HWSCRT example'/
     &        ' the output from the NCAR Control Data 7600 was'/,
     &        ' ierror = 0'/,
     &        ' discretization error = 5.36508e-04'/,
     &        ' required length of w array = 880'/,
     &        ' the output from your computer is'/,
     &        ' ierror =',i2,/,
     &        ' discretization error =',1pe12.5,/,
     &        ' required length of w array =',i4)
      end
      subroutine test02 ( kprint, ipass )

c*********************************************************************72
c
cc TEST02 illustrates the use of HWSPLR.
c
c  The equation solved is
c
c     (1/R)(D/DR)(R*(DU/DR)) + (1/R**2)(D/DTHETA)(DU/DTHETA) = 16*R**2
c
c  on the quarter-disk 0 .LT. R .LT. 1, 0 .LT. THETA .LT. PI/2 
c  with the boundary conditions
c
c     U(1,THETA) = 1 - COS(4*THETA), 0 .LE. THETA .LE. 1
c
c  and
c
c     (DU/DTHETA)(R,0) = (DU/DTHETA)(R,PI/2) = 0,  0 .LE. R .LE. 1.
c
c  (Note that the solution U is unspecified at R = 0.)
c  The R-interval will be divided into 50 panels and the
c  THETA-interval will be divided into 48 panels.
c
c  Modified:
c
c    21 September 2010
c
      implicit none

      real a
      real b
      real bda
      real bdb
      real bdc(51)
      real bdd(51)
      real c
      real d
      real elmbda
      real ermax
      real err
      real f(100,50)
      integer i
      integer idimf
      integer ierror
      integer ipass
      integer j
      integer kprint
      integer m
      integer mbdcnd
      integer n
      integer nbdcnd
      real pertrb
      real pi
      real pimach
      real r(51)
      real theta(49)
      real w(1200)
      real z
c
c  From dimension statement we get value of IDIMF.  Also note that W
c  is dimensioned 6*(N+1) + 8*(M+1).
c
      idimf = 100
      ermax=0.001
      a = 0.0
      b = 1.0
      m = 50
      mbdcnd = 5
      c = 0.
      pi = pimach()
      d = pi/2.0
      n = 48
      nbdcnd = 3
      elmbda = 0.0
c
c  Generate and store grid points for the purpose of computing
c  boundary data and the right side of the Poisson equation.
c
      do i=1,m+1
        r(i) = (i-1)/50.0
      end do

      do j=1,n+1
        theta(j) = (j-1)*pi/96.0
      end do
c
c  Generate boundary data.
c
      do i=1,m+1
        bdc(i) = 0.0
        bdd(i) = 0.0
      end do
c
c  BDA and BDB are dummy variables.
c
      do j=1,n+1
        f(m+1,j) = 1.0-cos(4.0*theta(j))
      end do
c
c  Generate right side of equation.
c
      do i=1,m
        do j=1,n+1
          f(i,j) = 16.0*r(i)**2
        end do
      end do

      call hwsplr(a,b,m,mbdcnd,bda,bdb,c,d,n,nbdcnd,bdc,bdd,elmbda,f,
     &  idimf,pertrb,ierror,w)
c
c  Compute discretization error.  The exact solution is
c    U(R,THETA) = R**4*(1 - COS(4*THETA))
c
      err = 0.0
      do i=1,m+1
        do j=1,n+1
          z = abs(f(i,j)-r(i)**4*(1.-cos(4.*theta(j))))
          if (z .gt. err) err = z
        end do
      end do

      ipass = 1
      if (err.gt.ermax) ipass = 0
      if (kprint.eq.0) return
      if (kprint.ge.2 .or. ipass.eq.0) then
         write(*,*)' '
         write(*,1001) ierror,err,int(w(1))
         if (ipass.eq.1) then
            write(*,*)'Test passed.'
         else
            write(*, *)'Test failed.'
         end if
      end if
      return

 1001 format (' TEST02: HWSPLR example'/,
     &       ' the output from the NCAR Control Data 7600 was',/,
     &       ' ierror = 0',/,
     &        ' discretization error = 6.19134e-04',/,
     &        ' required length of w array = 882',/,
     &       ' the output from your computer is',/,
     &       ' ierror =',i2,/,
     &        ' discretization error =',1pe12.5,/,
     &        ' required length of W array =',i4)
      end
      subroutine test03(kprint, ipass)

c*********************************************************************72
c
cc TEST03 illustrates the use of HWSCYL.
c
c  The equation solved is
c
c     (1/R)(D/DR)(R*(DU/DR)) + (D/DZ)(DU/DZ)
c
c     = (2*R*Z)**2*(4*Z**2 + 3*R**2)
c
c  on the rectangle 0 .LT. R .LT. 1, 0 .LT. Z .LT. 1 with the
c  boundary conditions
c
c     U(0,Z) unspecified
c                                            0 .LE. Z .LE. 1
c     (DU/DR)(1,Z) = 4*Z**4
c
c     and
c
c     (DU/DZ)(R,0) = 0
c                                            0 .LE. R .LE. 1
c     (DU/DZ)(R,1) = 4*R**4 .
c
c  The R-interval will be divided into 50 panels and the
c  Z-interval will be divided into 100 panels.
c
c  Modified:
c
c    21 September 2010
c
      implicit none

      real a
      real b
      real bda(101)
      real bdb(101)
      real bdc(51)
      real bdd(51)
      real c
      real d
      real elmbda
      real ermax
      real err
      real f(75,105)
      integer i
      integer idimf
      integer ierror
      integer ipass
      integer j
      integer kprint
      integer m
      integer mbdcnd
      integer n
      integer nbdcnd
      real pertrb
      real r(51)
      real w(1200)
      real x
      real z(101)
c
c  From dimension statement we get value of IDIMF.  Also note that W
c  is dimensioned 6*(N+1) + 8*(M+1).
c
      idimf = 75
      ermax=0.001
      a = 0.0
      b = 1.0
      m = 50
      mbdcnd = 6
      c = 0.0
      d = 1.0
      n = 100
      nbdcnd = 3
      elmbda = 0.0
c
c  Generate and store grid points for the purpose of computing
c  boundary data and the right side of the Poisson equation.
c
      do i=1,m+1
        r(i) = (i-1)/50.0
      end do

      do j=1,n+1
        z(j) = (j-1)/100.0
      end do
c
c  Generate boundary data.
c
      do j=1,n+1
        bdb(j) = 4.0*z(j)**4
      end do

      do i=1,m+1
        bdc(i) = 0.
        bdd(i) = 4.0*r(i)**4
      end do
c
c  BDA is a dummy variable.
c
c
c  Generate right side of equation.
c
      do i=1,m+1
        do j=1,n+1
          f(i,j) = 4.0*r(i)**2*z(j)**2*(4.*z(j)**2+3.*r(i)**2)
        end do
      end do

      call hwscyl(a,b,m,mbdcnd,bda,bdb,c,d,n,nbdcnd,bdc,bdd,elmbda,f,
     &  idimf,pertrb,ierror,w)
c
c  Compute discretization error by minimizing over all A the function
c  NORM(F(I,J) - A*1 - U(R(I),Z(J))).  The exact solution is
c    U(R,Z) = (R*Z)**4 + arbitrary constant.
c
      x = 0.0
      do i=1,m+1
        do j=1,n+1
          x = x+f(i,j)-(r(i)*z(j))**4
        end do
      end do

      x = x/((n+1)*(m+1))
      do i=1,m+1
        do j=1,n+1
          f(i,j) = f(i,j)-x
        end do
      end do

      err = 0.0
      do i=1,m+1
        do j=1,n+1
          x = abs(f(i,j)-(r(i)*z(j))**4)
          if (x .gt. err) err = x
        end do
      end do

      ipass = 1
      if (err.gt.ermax) ipass = 0
      if (kprint.eq.0) return
      if (kprint.ge.2 .or. ipass.eq.0) then
        write(*,*)' '
        write(*,1001) ierror,pertrb,err,int(w(1))
        if (ipass.eq.1) then
          write(*,*)'Test passed.'
        else
          write(*, *)'Test failed.'
        end if
      end if
      return

 1001 format (' TEST03: HWSCYL example'/
     &       ' the output from the NCAR Control Data 7600 was',/,
     &       ' ierror = 0',/,
     &       ' pertrb = 2.26734e-04',/,
     &        ' discretization error = 3.73672e-04',/,
     &        ' required length of w array = 1118',/,
     &       ' the output from your computer is',/,
     &       ' ierror =',i2,/,
     &       ' pertrb =',e12.5,/,
     &        ' discretization error =',1pe12.5,/,
     &        ' required length of w array =',i4)
      end
      subroutine test04(kprint, ipass)

c*********************************************************************72
c
cc TEST04 illustrates the use of HWSSSP.
c
c  Modified:
c
c    21 September 2010
c
      implicit none

      real bdpf
      real bdps
      real bdtf(73)
      real bdts
      real dphi
      real dtheta
      real elmbda
      real ermax
      real err
      real f(19,73)
      integer i
      integer idimf
      integer ierror
      integer ipass
      integer j
      integer kprint
      integer m
      integer mbdcnd
      integer n
      integer nbdcnd
      real pertrb
      real pf
      real pi
      real pimach
      real ps
      real sinp(73)
      real sint(19)
      real tf
      real ts
      real w(1200)
      real z
c
c  The value of IDIMF is the first dimension of F.  W is
c  dimensioned 11*(M+1)+6*(N+1)=647 since M=18 and N=72.
c
      pi = pimach()
      ermax = 0.005
      ts = 0.0
      tf = pi/2.0
      m = 18
      mbdcnd = 6
      ps = 0.0
      pf = pi+pi
      n = 72
      nbdcnd = 0
      elmbda = 0.0
      idimf = 19
c
c  Generate sines for use in subsequent computations
c
      dtheta = tf/m
      do i=1,m+1
        sint(i) = sin((i-1)*dtheta)
      end do

      dphi = (pi+pi)/n
      do j=1,n+1
        sinp(j) = sin((j-1)*dphi)
      end do
c
c  Compute right side of equation and store in F
c
      do j=1,n+1
        do i=1,m+1
           f(i,j) = 2.0-6.0*(sint(i)*sinp(j))**2
        end do
      end do
c
c  Store derivative data at the equator
c
      do j=1,n+1
        bdtf(j) = 0.0
      end do

      call hwsssp(ts,tf,m,mbdcnd,bdts,bdtf,ps,pf,n,nbdcnd,bdps,bdpf,
     &  elmbda,f,idimf,pertrb,ierror,w)
c
c  Compute discretization error.  Since problem is singular, the
c  solution must be normalized.
c
      err = 0.0
      do j=1,n+1
        do i=1,m+1
          z = abs(f(i,j)-(sint(i)*sinp(j))**2-f(1,1))
          if (z .gt. err) err = z
        end do
      end do

      ipass = 1
      if (err.gt.ermax) ipass = 0
      if (kprint.eq.0) return
      if (kprint.ge.2 .or. ipass.eq.0) then
        write(*,*)' '
        write(*,1001) ierror,err,int(w(1))
        if (ipass.eq.1) then
          write(*,*)'Test passed.'
        else
          write(*, *)'Test failed.'
        end if
      end if
      return

 1001 format (' TEST04: HWSSP example'/,
     &        ' the output from the NCAR Control Data 7600 was',/,
     &        ' ierror = 0',/,
     &        ' discretization error = 3.38107e-03',/,
     &        ' required length of w array = 600',/,
     &        ' the output from your computer is',/,
     &        ' ierror =',i2,/,
     &        ' discretization error =',1pe12.5, /,
     &        ' required length of W array =',i4)
      end
      subroutine test05(kprint, ipass)

c*********************************************************************72
c
cc TEST05 illustrates the use of HWSCSP.
c
c  Modified:
c
c    21 September 2010
c
      implicit none

      real bdrf
      real bdrs
      real bdtf(33)
      real bdts
      real ci4
      real dr
      real dtheta
      real elmbda
      real ermax
      real err
      real f(48,33)
      integer i
      integer idimf
      integer ierror
      integer intl
      integer ipass
      integer j
      integer kprint
      integer m
      integer mbdcnd
      integer n
      integer nbdcnd
      real pertrb
      real pi
      real pimach
      real r(33)
      real rf
      real rs
      real tf
      real theta(48)
      real ts
      real w(1200)
      real z
c
c  The value of IDIMF is the first dimension of F.  Since M=36, N=32,
c  L=N therefore K=5 and W is dimensioned 2*(L+1)*(K-1) + 6*(M+N)
c  + MAX(4*N,6*M) + 14 = 902.
c
      ermax=0.001
      pi = pimach()
      intl = 0
      ts = 0.0
      tf = pi/2.0
      m = 36
      mbdcnd = 6
      rs = 0.0
      rf = 1.0
      n = 32
      nbdcnd = 5
      elmbda = 0.0
      idimf = 48
c
c  Generate and store grid points for the purpose of computing the
c  boundary data and the right side of the equation.
c
      dtheta = tf/m
      do i=1,m+1
        theta(i) = (i-1)*dtheta
      end do

      dr = 1.0/n
      do j=1,n+1
        r(j) = (j-1)*dr
      end do
c
c  Generate normal derivative data at equator
c
      do j=1,n+1
        bdtf(j) = 0.0
      end do
c
c  Compute boundary data on the surface of the sphere
c
      do i=1,m+1
        f(i,n+1) = cos(theta(i))**4
      end do
c
c  Compute right side of equation
c
      do i=1,m+1
        ci4 = 12.0*cos(theta(i))**2
        do j=1,n
          f(i,j) = ci4*r(j)**2
        end do
      end do

      call hwscsp(intl,ts,tf,m,mbdcnd,bdts,bdtf,rs,rf,n,nbdcnd,bdrs,
     &  bdrf,elmbda,f,idimf,pertrb,ierror,w)
c
c  Compute discretization error
c
      err = 0.
      do i=1,m+1
        ci4 = cos(theta(i))**4
        do j=1,n
          z = abs(f(i,j)-ci4*r(j)**4)
          if (z .gt. err) err = z
        end do
      end do

      ipass = 1
      if (err.gt.ermax) ipass = 0
      if (kprint.ne.0) then
         if (kprint.ge.2 .or. ipass.eq.0) then
            write(*,*)' '
            write(*,1001) ierror,err,int(w(1))
            if (ipass.eq.1) then
               write(*, *)'Test passed.'
            else
               write(*, *)'Test failed.'
            end if
         end if
      end if

 1001 format (' TEST05: HWSCSP example 1'/
     &       ' the output from the NCAR Control Data 7600 was',/,
     &       ' ierror = 0',/,
     &        ' discretization error = 7.99842e-04',/,
     &        ' required length of w array = 775',/,
     &       ' the output from your computer is',/,
     &       ' ierror =',i2,/,
     &        ' discretization error =',1pe12.5,/,
     &        ' required length of w array =',i4)

      return
      end
      subroutine test06(kprint, ipass)

c*********************************************************************72
c
cc TEST06 illustrates the use of HWSCSP.
c
c  Discussion:
c
c    This program illustrates the use of HWSCSP to solve
c    a three dimensional problem which has longitudinal dependence.
c
c  Modified:
c
c    21 September 2010
c
      implicit none

      real bdrf
      real bdrs
      real bdtf(33)
      real bdts
      real dphi
      real dr
      real dtheta
      real elmbda
      real ermax
      real err
      real f(48,33)
      integer i
      integer idimf
      integer ierror
      integer intl
      integer ipass
      integer j
      integer kprint
      integer m
      integer mbdcnd
      integer n
      integer nbdcnd
      real pertrb
      real pi
      real pimach
      real r(33)
      real rf
      real rs
      real si
      real tf
      real theta(48)
      real ts
      real w(1200)
      real z
c
c  The value of IDIMF is the first dimension of F.  Since M=36, N=32,
c  L=N therefore K=5 and W is dimensioned 2*(L+1)*(K-1) + 6*(M+N)
c  + MAX(4*N,6*M) + 14 = 902.
c
      ermax=0.001
      pi = pimach()
      intl = 0
      ts = 0.0
      tf = pi/2.0
      m = 36
      mbdcnd = 2
      rs = 0.0
      rf = 1.0
      n = 32
      nbdcnd = 1
      idimf = 48
      dphi = pi/72.0
      elmbda = -2.0*(1.0-cos(dphi))/dphi**2
c
c  Generate and store grid points for the purpose of computing the
c  boundary data and the right side of the equation.
c
      dtheta = tf/m
      do i=1,m+1
        theta(i) = (i-1)*dtheta
      end do

      dr = 1.0/n
      do j=1,n+1
        r(j) = (j-1)*dr
      end do
c
c  Generate normal derivative data at equator
c
      do j=1,n+1
        bdtf(j) = 0.0
      end do
c
c  Compute boundary data on the surface of the sphere
c
      do i=1,m+1
        f(i,n+1) = sin(theta(i))
      end do
c
c  Compute right side of the equation
c
      do j=1,n
        do i=1,m+1
          f(i,j) = 0.0
        end do
      end do

      call hwscsp(intl,ts,tf,m,mbdcnd,bdts,bdtf,rs,rf,n,nbdcnd,bdrs,
     &  bdrf,elmbda,f,idimf,pertrb,ierror,w)
c
c  Compute discretization error   (Fourier coefficients)
c
      err = 0.
      do i=1,m+1
        si = sin(theta(i))
        do j=1,n+1
          z = abs(f(i,j)-r(j)*si)
          if (z .gt. err) err = z
        end do
      end do

      if (err.gt.ermax) ipass = 0
      if (kprint.eq.0) return
      if (kprint.ge.2 .or. ipass.eq.0) then
        write(*,*)' '
        write(*,1001) ierror,err,int(w(1))
        if (ipass.eq.1) then
          write(*,*)'Test passed.'
        else
          write(*, *)'Test failed.'
        end if
      end if
      return

 1001 format (' TEST06: HWSCSP example 2'/
     &       ' the output from the NCAR Control Data 7600 was',/
     &       ' ierror = 0',/
     &        ' discretization error = 5.86824e-05',/
     &        ' required length of w array = 775',/
     &       ' the output from your computer is'/,
     &       ' ierror =',i2,/,
     &        ' discretization error =',1pe12.5,/,
     &        ' required length of w array =',i4)
      end
      subroutine test07(kprint, ipass)

c*********************************************************************72
c
cc TEST07 illustrates the use of GENBUN.
c
c  Modified:
c
c    21 September 2010
c
      implicit none

      real a(20)
      real b(20)
      real c(20)
      real deltax
      real deltay
      real dysq
      real ermax
      real err
      real f(25,130)
      integer i
      integer idimy
      integer ierror
      integer ipass
      integer j
      integer kprint
      integer m
      integer mperod
      integer n
      integer nperod
      real pi
      real pimach
      real s
      real t
      real w(1200)
      real x(20)
      real y(120)
      real z
c
c  From dimension statement we get value of IDIMY.  Also note that
c  W(.) is dimensioned 6*N + 5*M.
c
      ermax=0.01
      idimy = 25
      mperod = 1
      m = 20
      deltax = 1.0/m
      nperod = 0
      n = 120
      pi = pimach()
      deltay = 2.0*pi/n
c
c  Set the grid points.
c
      do i=1,m
        x(i) = real(i-1)/real(m)
      end do

      do j=1,n
        y(j) = -pi + 2*(j-1)*pi/real(n)
      end do
c
c  Generate coefficients.
c
      s = (deltay/deltax)**2
      t = s*deltax
      a(1) = 0.0
      b(1) = -2.0*s
      c(1) = 2.0*s
      do i=2,m
        a(i) = (1.0+x(i))**2*s + (1.0+x(i))*t
        c(i) = (1.0+x(i))**2*s - (1.0+x(i))*t
        b(i) = -2.0*(1.0+x(i))**2*s
      end do
      c(m) = 0.0
c
c  Generate right side of equation for I = 1 showing introduction of
c  boundary data.
c
      dysq = deltay**2
      do j=1,n
        f(1,j) = dysq*(11.0 + 8.0/deltax)*sin(y(j))
      end do
c
c  Generate the right hand side.
c
      do i=2,m-1
        do j=1,n
          f(i,j) = dysq*3.0*(1.0+x(i))**4*sin(y(j))
        end do
      end do
c
c  Generate right side for I=M showing introduction of
c  boundary data.
c
      do j=1,n
        f(m,j)=dysq*(3.0*(1.+x(m))**4 - 16.0*((1.0+x(m))/deltax)**2
     &    +16.0*(1.0+x(m))/deltax)*sin(y(j))
      end do
c
c  Call GENBUN to solve the system.
c
      call genbun(nperod,n,mperod,m,a,b,c,idimy,f,ierror,w)
c
c  Compute the discretization error.  The exact solution is
c    U(X,Y) = (1+X)**4*SIN(Y)
c
      err = 0.0
      do i=1,m
        do j=1,n
          z = abs(f(i,j)-(1.0+x(i))**4*sin(y(j)))
          if (z .gt. err) err = z
        end do
      end do

      ipass = 1
      if (err.gt.ermax) ipass = 0
      if (kprint.eq.0) return
      if (kprint.ge.2 .or. ipass.eq.0) then
        write(*,*)' '
        write(*,1001) ierror, err, int(w(1))
        if (ipass.eq.1) then
          write(*,*)'Test passed.'
        else
          write(*, *)'Test failed.'
        end if
      end if
      return

 1001 format (' TEST07: GENBUN example'/
     &       ' the output from the NCAR Control Data 7600 was',/,
     &       ' ierror = 0',/,
     &        ' discretization error = 7.94113e-03',/,
     &        ' required length of w array = 740',/,
     &       ' the output from your computer is'/,
     &       ' ierror =',i2,/,
     &        ' discretization error =',1pe12.5,/,
     &        ' required length of W array =',i4)
      end
      subroutine test08(kprint, ipass)

c*********************************************************************72
c
cc TEST08 illustrates the use of BLKTRI.
c
c  Modified:
c
c    21 September 2010
c
      implicit none

      integer m
      integer n

      parameter ( m = 50 )
      parameter ( n = 63 )

      integer idimy
      parameter ( idimy = m )

      real am(m)
      real an(n)
      real bm(m)
      real bn(n)
      real cm(m)
      real cn(n)
      real deltas
      real deltat
      real ermax
      real err
      integer flag
      real hds
      real hdt
      integer i
      integer ierror
      integer ipass
      integer j
      integer kprint
      integer mp
      integer np
      real s(m)
      real t(n)
      real tds
      real tdt
      real temp1
      real temp2
      real temp3
      real true
      real w(1952)
      real y(idimy,n)

      ermax = 0.001
      np = 1
      mp = 1
c
c  Generate and store the grid points.
c
      deltas = 1.0/(m+1)
      do i = 1, m
        s(i) = i * deltas
      end do

      deltat = 1.0/(n+1)
      do j = 1, n
        t(j) = j * deltat
      end do
c
c  Compute the coefficients AM,BM,CM corresponding to the S direction
c
      hds = deltas/2.0
      tds = deltas+deltas

      do i = 1, m
        temp1 = 1.0/(s(i)*tds)
        temp2 = 1.0/((s(i)-hds)*tds)
        temp3 = 1.0/((s(i)+hds)*tds)
        am(i) = temp1*temp2
        cm(i) = temp1*temp3
        bm(i) = -(am(i)+cm(i))
      end do
c
c  Compute the coefficients AN,BN,CN corresponding to the T direction
c
      hdt = deltat/2.0
      tdt = deltat+deltat

      do j = 1, n
        temp1 = 1.0/(t(j)*tdt)
        temp2 = 1.0/((t(j)-hdt)*tdt)
        temp3 = 1.0/((t(j)+hdt)*tdt)
        an(j) = temp1*temp2
        cn(j) = temp1*temp3
        bn(j) = -(an(j)+cn(j))
      end do
c
c  Compute the right hand side of the equation.
c
      do j = 1, n
        do i = 1, m
          y(i,j) = 3.75 * s(i) * t(j) * ( s(i)**4 + t(j)**4 )
        end do
      end do
c
c  Include nonhomogeneous boundary into right side. 
c
c  Note that the corner Y(M,N) includes contributions 
c  from both boundaries.
c
      do j = 1, n
        y(m,j) = y(m,j) - cm(m) * t(j)**5
      end do

      do i = 1, m
        y(i,n) = y(i,n) - cn(n) * s(i)**5
      end do
c
c  Initialize the code.
c
      flag = 0

      call blktri(flag,np,n,an,bn,cn,mp,m,am,bm,cm,idimy,y,ierror,w)

      if ( ierror .ne. 0 ) then
        ipass = 0
        return
      end if
c
c  Solve the system.
c
      flag = 1

      call blktri(flag,np,n,an,bn,cn,mp,m,am,bm,cm,idimy,y,ierror,w)

      if ( ierror .ne. 0 ) then
        ipass = 0
        return
      end if
c
c  Compute the discretization error.
c
      err = 0.0
      do j = 1, n
        do i = 1, m
          true = ( s(i) * t(j) )**5
          err = max ( err, abs ( y(i,j) - true ) )
        end do
      end do

      if ( err .gt. ermax ) then
        ipass = 0
      else
        ipass = 1
      end if
c
c  Report the results.
c
      if ( 
     &  ( kprint.eq. 1 .and. ipass.eq.0 ) .or.
     &  kprint .ge. 2 ) then

        write(*,*)' '
        write(*,1001) ierror,err,int(w(1))

        if (ipass.eq.1) then
          write(*, *) 'Test passed'
        else
          write(*, *) 'Test failed'
        end if

      end if

      return

 1001 format (' TEST08: BLKTRI example'/
     &        ' the output from the NCAR Control Data 7600 was' /,
     &        ' ierror = 0' /,
     &        ' discretization error = 1.6478e-05' /,
     &        ' required length of w array = 823' /,
     &        ' the output from your computer is' /,
     &        ' ierror =',i2, /,
     &        ' discretization error =',1pe12.5, /,
     &        ' required length of w array =', i4)
      end
      subroutine test09 ( )

c*********************************************************************72
c
cc TEST09 tests HWSCRT, 2-D Cartesian coordinate Helmholtz equation.
c
c  Discussion:
c
c    The equation has the form:
c
c      Uxx + Uyy + Lambda*U = F(X,Y)
c
c    We use the simple test case 
c
c      U = X**2 + Y**2, 
c
c    for X in (0,1), Y in (0,1).
c
c  Modified:
c
c    21 September 2010
c
      implicit none

      integer m
      integer n
      integer nwork

      parameter (m=5)
      parameter (n=5)
c
c  NWORK=4*(N+1)+(M+1)*(13+ALOG2(N+1))
c
      parameter (nwork=4*(n+1)+(m+1)*(13+3))

      real a
      real b
      real bda(n+1)
      real bdb(n+1)
      real bdc(m+1)
      real bdd(m+1)
      real c
      real d
      real dx
      real dy
      real elmbda
      real f(m+1,n+1)
      integer i
      integer idimf
      integer ierror
      integer j
      integer mbdcnd
      integer nbdcnd
      real pertrb
      real w(nwork)
      real x
      real y

      write(*,*)' '
      write(*,*)'TEST09'
      write(*,*)'test hwscrt, helmholtz equation,'
      write(*,*)'2-d cartesian coordinates.'
      write(*,*)' '
c
c  Define corners of region:
c
      a=0.0
      b=1.0
      c=0.0
      d=1.0
c
c  Define boundary conditions codes:
c
      mbdcnd=1
      nbdcnd=1
c
c  Set boundary conditions:
c
      dy=(d-c)/real(n)
      dx=(b-a)/real(m)

      x=a
      do j=1,n+1
        y=(j-1)*dy
        f(1,j)=x*x+y*y
      end do

      x=b
      do j=1,n+1
        y=(j-1)*dy
        f(m+1,j)=x*x+y*y
      end do

      y=c
      j=1
      do i=1,m+1
        x=(i-1)*dx
        f(i,j)=x*x+y*y
      end do

      y=d
      do i=1,m+1
        x=(i-1)*dx
        f(i,n+1)=x*x+y*y
      end do
c
c  Set right hand side:
c
      do i=2,m
        do j=2,n
          f(i,j)=4.0
        end do
      end do
c
c  Specify lambda.
c
      elmbda=0.0

      idimf=m+1

      call hwscrt(a,b,m,mbdcnd,bda,bdb,c,d,n,nbdcnd,bdc,bdd,elmbda,f,
     &  idimf,pertrb,ierror,w)

      write(*,*)' '

      if ( ierror.ne.0)then

        write(*,*)'hwscrt returned error flag ierror=',ierror

      else

        write(*,*)'hwscrt computed approximate solution.'
        write(*,*)' '

        if ( pertrb.ne.0.0)then
          write(*,*)'required perturbation to be subtracted=',pertrb
          write(*,*)' '
        end if

        do i=1,m+1
          write(*,'(1x,6g12.4)')(f(i,j),j=1,n+11)
        end do

      end if

      return
      end
      subroutine test10 ( )

c*********************************************************************72
c
cc TEST10 tests HWSPLR, 2-D polar coordinate system Helmholtz equation.
c
c  Discussion:
c
c    The equation has the form:
c
c      (1/R) D/DR R DU/DR + (1/R*R) D/DTHETA DU/DTHETA + Lambda*U = F(X,Y)
c
c    We use the simple test case U=R*R, for R in (1,2), Y in (0,PI/2).
c
c  Modified:
c
c    21 September 2010
c
      implicit none

      integer m
      integer n
      integer nwork

      parameter (m=5)
      parameter (n=5)
c
c  NWORK=4*(N+1)+(M+1)*(13+ALOG2(N+1))
c
      parameter (nwork=4*(n+1)+(m+1)*(13+3))

      real a
      real b
      real bda(n+1)
      real bdb(n+1)
      real bdc(m+1)
      real bdd(m+1)
      real c
      real d
      real elmbda
      real f(m+1,n+1)
      integer i
      integer idimf
      integer ierror
      integer j
      integer mbdcnd
      integer nbdcnd
      real pertrb
      real r
      real w(nwork)

      write(*,*)' '
      write(*,*)'TEST10'
      write(*,*)'test hwsplr, helmholtz equation,'
      write(*,*)'2-d polar coordinates.'
      write(*,*)' '
c
c  Define corners of region:
c
      a = 1.0
      b = 2.0
      c = 0.0
      d = 3.14159265/2.0
c
c  Define boundary conditions codes:
c
      mbdcnd = 1
      nbdcnd = 3
c
c  Set boundary conditions:
c
      r = a
      do j = 1, n + 1
        f(1,j) = r * r
      end do

      r = b
      do j = 1, n + 1
        f(m+1,j) = r * r
      end do

      do i = 1, m + 1
        bdc(i) = 0.0
      end do

      do i = 1, m + 1
        bdd(i) = 0.0
      end do
c
c  Set the right hand side:
c
      do i = 2, m
        do j = 1, n + 1
          f(i,j) = 4.0
        end do
      end do
c
c  Specify lambda.
c
      elmbda = 0.0
c
c  IDIMF is the first dimension of F as declared here.
c
      idimf = m + 1

      call hwsplr ( a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, 
     &  bdd, elmbda, f, idimf, pertrb, ierror, w )

      write(*,*)' '

      if ( ierror.ne.0)then

        write(*,*)'HWSPLR returned error flag ierror=',ierror

      else

        write(*,*)'HWSPLR computed approximate solution.'
        write(*,*)' '

        if ( pertrb .ne. 0.0 ) then
          write(*,*)'required perturbation to be subtracted=',pertrb
          write(*,*)' '
        end if

        do i=1,m+1
          write(*,'(1x,6g12.4)')(f(i,j),j=1,n+1)
        end do

      end if

      return
      end
      subroutine test11 ( )

c*********************************************************************72
c
cc TEST11 uses BLKTRI to solve a simple PDE.
c
c  The domain is the unit square, 0 <= x,y <= 1.
c
c  There are M vertical grid lines, and N horizontal grid lines,
c  equally spaced in the interior of the square.  For simplicity,
c  we take M = N.  This does not count the extra lines that define the 
c  boundaries.
c
c  At each interior node, we wish to approximate the solution
c  of a partial differential equation:
c
c    d2u/dx2 + d2u/dy2 = 4.
c
c  Along the boundaries, we apply the Dirichlet condition
c    u(x,y) = x**2 + y**2.
c  (In fact, this formula yields the exact solution to the PDE over
c  the whole square).
c
c  Our partial differential equation is approximated by
c
c                     u(i,j+1)
c    + u(i-1,j) - 4.0*u(i,j)   + u(i+1,j)
c                   + u(i,j-1)             = 4.0 * h**2
c
c  Here "H" is the grid spacing.
c
c  This equation is written at every interior grid point.  For
c  grid points next to the boundary, it is necessary to evaluate
c  any references to U on the boundary, and carry those terms to
c  the right hand side.
c
c  The FISHPACK routine BLKTRI is used to handle this system.
c
c  Modified:
c
c    21 September 2010
c
      implicit none

      integer m
      parameter ( m = 5 )
c
c  We choose to set N = M to simplify the equations.
c
      integer n
      parameter ( n = m )

      integer idimy
      parameter ( idimy = m )

      integer nw
      parameter ( nw = 100 )

      real am(m)
      real an(n)
      real bm(m)
      real bn(n)
      real cm(m)
      real cn(n)
      real delta
      real err
      real errmax
      integer flag
      integer i
      integer ierror
      integer imax
      integer j
      integer jmax
      integer mp
      integer np
      real u
      real ubot
      real uleft
      real urite
      real utop
c
c  Replace this dimension by something more reasonable.
c
      real w(nw)
      real x(0:m+1)
      real y(0:n+1)
      real yy(idimy,n)
c
c  We set MP and NP to indicate that the boundary conditions are not periodic.
c
      np = 1
      mp = 1
c
c  This is the spacing between nodes.
c
      delta = 1.0 / real ( n + 1 )
c
c  Generate and store the grid points.
c  This is just for our convenience.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST11'
      write ( *, '(a)' ) '  An example of the direct use of BLKTRI.'
 
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'I, X(I)'
      write ( *, '(a)' ) ' '

      do i = 0, m + 1
        x(i) = real ( i ) / real ( m + 1 )
        write ( *, * ) i, x(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'J, Y(J)'
      write ( *, '(a)' ) ' '
      do j = 0, n + 1
        y(j) = real ( j ) / real ( n+1 )
        write(*,*) j, y(j)
      end do
c
c  Compute the coefficients AM, BM, CM corresponding to the X direction.
c
      do i = 1, m
        am(i) = 1.0
        bm(i) = -2.0
        cm(i) = 1.0
      end do
c
c  Compute the coefficients AN, BN, CN corresponding to the Y direction.
c
      do j = 1, n
        an(j) = 1.0
        bn(j) = -2.0
        cn(j) = 1.0
      end do
c
c  Set the right hand side of the equations.
c
      do j = 1, n
        do i = 1, m
          yy(i,j) = 4.0 * delta**2
        end do
      end do
c
c  Set the boundary conditions implicitly, by carrying terms involving
c  boundary values to the right hand side, and resetting the corresponding
c  coefficients to 0.
c
      i = 1
      do j = 1, n
        uleft = x(i-1)**2 + y(j)**2
        yy(i,j) = yy(i,j) - am(1) * uleft
      end do
      am(1) = 0.0

      i = m
      do j = 1, n
        urite = x(i+1)**2 + y(j)**2
        yy(i,j) = yy(i,j) - cm(m) * urite
      end do
      cm(m) = 0.0

      j = 1
      do i = 1, m
        ubot = x(i)**2 + y(j-1)**2
        yy(i,j) = yy(i,j) - an(1) * ubot
      end do
      an(1) = 0.0

      j = n
      do i = 1, m
        utop = x(i)**2 + y(j+1)**2
        yy(i,j) = yy(i,j) - cn(n) * utop
      end do
      cn(n) = 0.0

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'I,J,YY(I,J)'
      write ( *, '(a)' ) ' '
      do i = 1, m
        do j = 1, n
          write ( *, '(2x,i4,2x,i4,2x,g14.6)' ) i, j, yy(i,j)
        end do
        write ( *, '(a)' ) ' '
      end do
c
c  Initialize the code.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Call with FLAG = 0 to initialize.'

      flag = 0

      call blktri ( flag, np, n, an, bn, cn, mp, m, am, bm, cm, idimy,
     &  yy, ierror, w )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Fatal error!'
        write ( *, '(a,i8)' ) '  BLKTRI returned IERROR = ', ierror
        return
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 
     &  '  BLKTRI requires W dimension at least ', int ( w(1) )

      if ( nw .lt. int ( w(1) ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)') 'Fatal error!'
        write ( *, '(a,i8,a)' ) 
     &    '  User allocated only ', nw, ' entries.'
        stop
      end if
c
c  Solve the system.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Call with FLAG = 1 to solve the system.'

      flag = 1

      call blktri ( flag, np, n, an, bn, cn, mp, m, am, bm, cm, idimy,
     &  yy, ierror, w )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Fatal error!'
        write ( *, '(a)' ) '  BLKTRI returned IERROR = ', ierror
        return
      end if
c
c  Compute the discretization error.
c
      errmax = 0.0
      imax = 0
      jmax = 0

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'I,J,TRUE,YY(I,J),ERR'

      do j = 1, n
        write ( *, '(a)' ) ' '
        do i = 1, m
          u = x(i)**2 + y(j)**2
          err = abs ( yy(i,j) - u ) 
          write ( *, '(2x,i4,2x,i4,g14.6,2x,g14.6,2x,g14.6)' ) 
     &      i, j, u, yy(i,j), err
          if ( errmax .lt. err ) then
            imax = i
            jmax = j
            errmax = err
          end if
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Maximum error is ', errmax
      write ( *, '(a,i4,2x,i4)' ) '  I, J = ', imax, jmax

      return
      end
