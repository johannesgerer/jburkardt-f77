      program main

c*********************************************************************72
c
cc MAIN is the main program for INOUT.
c
c  Discussion:
c
c    INOUT solves the steady incompressible Navier-Stokes equations
c    in the 2D "INOUT" region.
c
c    The fluid flow problem is formulated in terms of
c    primitive variables - u,v, and p.
c
c    u_t - laplacian u + (u.grad)u + grad p = f
c                                div u = 0
c
c    Boundary conditions:  (u,v)=(0,0) on top
c                          (u,v)=0 on left, right and bottom
c
c    This version uses finite element techniques
c    with piecewise linear functions on triangles to approximate
c    the pressure and quadratics on triangles for the velocity
c    (Taylor-Hood element), isoparametric element
c
c  Input files:
c
c    FUINI.DAT contains the initial values of the solution coefficients.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 May 2004
c
c  Author:
c
c    Hyung-Chun Lee,
c    Department of Mathematics,
c    Ajou University, Korea
c
c  Local parameters:
c
c    FILENAME is the name to be used for the first output file.
c    Each subsequent output file created by the program will have an
c    incremented name.
c
      implicit double precision(a-h,o-z)

      parameter ( nx = 21 )
      parameter ( ny = 21 )
      parameter ( mx = 2 * nx - 1 )
      parameter ( my = 2 * ny - 1 )
      parameter ( maxel = 2 * (nx-1) * (ny-1) )
      parameter ( maxnd = mx * my )
      parameter ( maxun = 2 * mx * my + nx * ny )
      parameter ( minun = 27 * ny )

      double precision a(minun,maxun)
      double precision area(maxel)
      double precision f(maxun)
      character*20 filename
      double precision g(maxun)
      double precision gg(maxun)
      double precision gp(maxun)
      integer indx(maxnd,2)
      integer insc(maxnd)
      integer ipivot(maxun)
      integer node(maxnd,6)
      double precision uold(maxun)
      double precision xc(maxnd)
      double precision xm(maxel,3)
      double precision yc(maxnd)
      double precision ym(maxel,3)

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'INOUT:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Solve the Navier Stokes fluid flow'
      write ( *, '(a)' ) '  equations in the INOUT region,'
      write ( *, '(a)' ) '  using finite elements.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Maximum number of nodes =    ', maxnd
      write ( *, '(a,i6)' ) '  Maximum number of elements = ', maxel

      filename = 'up000.txt'
      pi = 4.0d0 * datan ( 1.0d0 )
c
c  Set functions for input profile
c
c   set input data
c
       open(unit=0,file='xy.dat')
       open(unit=7,file='up.dat')

       iwrite=0
       xlngth=1.0D+00
       ylngth=1.0D+00
       visc_inv = 300.0D+00
       mrow1=minun
       nsim=3
       nsteps=10
       tolns=1.0D-06
       tolopt=1.0D-06
c
c  setgrd constructs grid, numbers unknowns, calculates areas,
c  and points for midpoint quadrature rule
c  SETGRD sets bandwidth and neqn1
c
      call setgrd(xc,yc,area,xm,ym,xlngth,ylngth,
     &  node,indx,insc,nlband,nband,
     &  nx,ny,nelemn,np,nnodes,nuk,nquad,neqn1,
     &  iwrite,maxnd,maxel)

      write ( *, '(a,i6)' ) '  Number of nodes =    ', np
      write ( *, '(a,i6)' ) '  Number of elements = ', nelemn

      call grid_write ( maxel, maxnd, xc, yc, node, nelemn, np )

      nuband=nlband
      nrow1=nlband+nlband+nuband+1
      ncol1=neqn1
      ny2=ny+ny-1

      deltat=0.01D+00
      rdel=1.0D+00/deltat

      do i=1,neqn1
        f(i)=0.0D+00
      end do
c
c  Initialize the solution.
c
      call solution_init ( neqn1, uold )
c
c  Carry out the time iteration.
c
      do iter=1,5000

        open(unit=10,file='fu.dat')

        darg=0.1d0*dble(iter)*deltat*pi-0.5d0*pi
        alpha=2.0D+00*dsin(darg)+3.0D+00

        do i=1,neqn1
          g(i)=f(i)
          f(i)=0.0d0
        end do

        call nstoke(xc,yc,area,xm,ym,
     &    a,f,g,uold,visc_inv,tolns,xlngth,ylngth,
     &    node,indx,insc,ipivot,mrow1,
     &    nlband,nuband,nband,nrow1,ncol1,
     &    nelemn,np,nnodes,nuk,nquad,neqn1,
     &    nsteps,nsim,iwrite,maxnd,maxel,rdel,alpha)
c
c  Save u=(gx,gy) to 'ue.dat' for 't.f'
c
        do i=1,neqn1
          uold(i)=f(i)
        end do

        moit=mod(iter,10)

        if (moit.eq.0) then

          call file_name_inc ( filename )
          write ( *, * ) 'Creating ', filename
          open (unit=1,file=filename)

          do ic =1,np

            iuku=indx(ic,1)
            iukv=indx(ic,2)
            iukp=insc(ic)

            if(iuku.eq.0) then
              gx=0.0d0
              if(iukp.eq.0) then
                pp=0.0d0
              else
                pp=f(iukp)
              end if
              go to 54
            else if(iuku.eq.-1) then
              gx=ubdry(1,ic,xc,yc,ylngth,alpha)
              if(iukp.eq.0) then
                pp=0.0d0
              else
                pp=f(iukp)
              end if
              go to 54
            else
              gx=f(iuku)
            end if

            if(iukp.eq.0) then
              pp=0.0d0
            else
              pp=f(iukp)
            end if

54          continue

            if(iukv.le.0) then
              gy=0.0d0
            else
              gy=f(iukv)
            end if

            if(iukp.eq.0) then
              pp=0.0d0
            else
              pp=f(iukp)
            end if

            g(ic)=gx
            gg(ic)=gy
            gp(ic)=pp
            write(1,3031) gx,gy

          end do

          close (unit=1)

        end if

        write(10,*) iter
        do i=1,neqn1
          write(10,*) f(i)
        end do

        close(10)

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'INOUT:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

3031  format(2e25.15)
3033  format(6d10.3)
1040  format(1x)
1050  format(' number of unknowns=',i4,'   number of triangles= ',i4)
1100  format(' nlband=',i4,' nuband=', i4,' nband=',i4/
     &   ' nrow1=',i4,' ncol1=',i4,' # of pts=',i4)
1300  format(' nx=  ',i4)
1310  format(' ny = ',i4)
1410  format(' visc_inv=',f14.7)

      stop
      end
      function ch_is_digit ( c )

c*********************************************************************72
c
cc CH_IS_DIGIT returns TRUE if a character is a decimal digit.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character C, the character to be analyzed.
c
c    Output, logical CH_IS_DIGIT, TRUE if C is a digit, FALSE otherwise.
c
      implicit none

      character c
      logical ch_is_digit

      if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then
        ch_is_digit = .true.
      else
        ch_is_digit = .false.
      end if

      return
      end
      subroutine ch_to_digit ( c, digit )

c*********************************************************************72
c
cc CH_TO_DIGIT returns the integer value of a base 10 digit.
c
c  Example:
c
c     C   DIGIT
c    ---  -----
c    '0'    0
c    '1'    1
c    ...  ...
c    '9'    9
c    ' '    0
c    'X'   -1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character C, the decimal digit, '0' through '9' or blank
c    are legal.
c
c    Output, integer DIGIT, the corresponding integer value.  If C was
c    'illegal', then DIGIT is -1.
c
      implicit none

      character c
      integer digit

      if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

        digit = ichar ( c ) - 48

      else if ( c .eq. ' ' ) then

        digit = 0

      else

        digit = -1

      end if

      return
      end
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)

c*********************************************************************72
c
cc DAXPY adds a multiple of one vector to another.
c
c     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
c     JACK DONGARRA, LINPACK, 3/11/78.
c
      implicit double precision(a-h,o-z)

      DOUBLE PRECISION DX(1),DY(1),DA
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N

      IF(N.LE.0)RETURN
      IF (DA .EQ. 0.0D0) RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
c
c    CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
c    NOT EQUAL TO 1
c
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
      end do
      RETURN
c
c   CODE FOR BOTH INCREMENTS EQUAL TO 1
c
c
c   CLEAN-UP LOOP
c
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I + 1) = DY(I + 1) + DA*DX(I + 1)
        DY(I + 2) = DY(I + 2) + DA*DX(I + 2)
        DY(I + 3) = DY(I + 3) + DA*DX(I + 3)
   50 CONTINUE
      RETURN
      END
      function ddot ( n, dx, INCX, DY, INCY )

c*********************************************************************72
c
cc DDOT forms the dot product of two vectors.
c
c     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
c     JACK DONGARRA, LINPACK, 3/11/78.
c
      implicit none

      double precision ddot
      double precision dtemp
      DOUBLE PRECISION dx(*)
      double precision dy(*)
      INTEGER I
      integer INCX
      integer INCY
      integer IX
      integer IY
      integer M
      integer N

      DDOT = 0.0D0
      DTEMP = 0.0D0

      if ( N .LE. 0 ) then
        RETURN
      end if

      IF ( INCX .EQ. 1 .AND. INCY .EQ. 1 ) GO TO 20
c
c    CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
c    NOT EQUAL TO 1
c
      IX = 1
      IY = 1
      IF ( INCX .LT. 0 ) IX = (-N+1)*INCX + 1
      IF ( INCY .LT. 0 ) IY = (-N+1)*INCY + 1
      do I = 1, N
        DTEMP = DTEMP + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
      end do
      DDOT = DTEMP
      RETURN
c
c    CODE FOR BOTH INCREMENTS EQUAL TO 1
c
c
c    CLEAN-UP LOOP
c
   20 M = MOD(N,5)

      do I = 1, M
        DTEMP = DTEMP + DX(I)*DY(I)
      end do
      IF ( N .LT. 5 ) GO TO 60

      do I = M+1,N,5
      DTEMP = DTEMP + DX(I)*DY(I) + DX(I + 1)*DY(I + 1) +
     &  DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
      end do

   60 DDOT = DTEMP
      RETURN
      END
      subroutine dgbfa(abd,lda,n,ml,mu,ipvt,info)

c*********************************************************************72
c
cc DGBFA factors a double precision band matrix by elimination.
c
c     dgbfa is usually called by dgbco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c
c     on entry
c
c        abd     double precision(lda, n)
c                contains the matrix in band storage.  the columns
c                of the matrix are stored in the columns of  abd  and
c                the diagonals of the matrix are stored in rows
c                ml+1 through 2*ml+mu+1 of  abd .
c                see the comments below for details.
c
c        lda     integer
c                the leading dimension of the array  abd .
c                lda must be .ge. 2*ml + mu + 1 .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c                0 .le. ml .lt. n .
c
c        mu      integer
c                number of diagonals above the main diagonal.
c                0 .le. mu .lt. n .
c                more efficient if  ml .le. mu .
c     on return
c
c        abd     an upper triangular matrix in band storage and
c                the multipliers which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgbsl will divide by zero if
c                     called.  use  rcond  in dgbco for a reliable
c                     indication of singularity.
c
c     band storage
c
c           if  a  is a band matrix, the following program segment
c           will set up the input.
c
c                   ml = (band width below the diagonal)
c                   mu = (band width above the diagonal)
c                   m = ml + mu + 1
c                   do 20 j = 1, n
c                      i1 = max0(1, j-mu)
c                      i2 = min0(n, j+ml)
c                      do i = i1, i2
c                         k = i - j + m
c                         abd(k,j) = a(i,j)
c                      end do
c                20 continue
c
c           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
c           in addition, the first  ml  rows in  abd  are used for
c           elements generated during the triangularization.
c           the total number of rows needed in  abd  is  2*ml+mu+1 .
c           the  ml+mu by ml+mu  upper left triangle and the
c           ml by ml  lower right triangle are not referenced.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c     fortran max0,min0
c
      implicit double precision(a-h,o-z)

      integer lda,n,ml,mu,ipvt(1),info
      double precision abd(lda,1)

      double precision t
      integer i,idamax,i0,j,ju,jz,j0,j1,k,kp1,l,lm,m,mm,nm1

      m = ml + mu + 1
      info = 0
c
c     zero initial fill-in columns
c
      j0 = mu + 2
      j1 = min0(n,m) - 1
      if (j1 .lt. j0) go to 30
      do jz = j0, j1
         i0 = m + 1 - jz
         do i = i0, ml
           abd(i,jz) = 0.0d0
         end do
      end do
   30 continue
      jz = j1
      ju = 0
c
c  gaussian elimination with partial pivoting
c
      nm1 = n - 1
      if (nm1 .lt. 1) go to 130
      do 120 k = 1, nm1
         kp1 = k + 1
c
c        zero next fill-in column
c
         jz = jz + 1
         if (jz .gt. n) go to 50
            do i = 1, ml
               abd(i,jz) = 0.0d0
            end do
   50    continue
c
c        find l = pivot index
c
         lm = min(ml,n-k)
         l = idamax(lm+1,abd(m,k),1) + m - 1
         ipvt(k) = l + k - m
c
c        zero pivot implies this column already triangularized
c
         if (abd(l,k) .eq. 0.0d0) go to 100
c
c           interchange if necessary
c
            if (l .eq. m) go to 60
               t = abd(l,k)
               abd(l,k) = abd(m,k)
               abd(m,k) = t
   60       continue
c
c           compute multipliers
c
            t = -1.0d0/abd(m,k)
            call dscal(lm,t,abd(m+1,k),1)
c
c           row elimination with column indexing
c
            ju = min0(max0(ju,mu+ipvt(k)),n)
            mm = m
            if (ju .lt. kp1) go to 90
            do 80 j = kp1, ju
               l = l - 1
               mm = mm - 1
               t = abd(l,j)
               if (l .eq. mm) go to 70
                  abd(l,j) = abd(mm,j)
                  abd(mm,j) = t
   70          continue
               call daxpy(lm,t,abd(m+1,k),1,abd(mm+1,j),1)
   80       continue
   90       continue
         go to 110
  100    continue
            info = k
  110    continue
  120 continue
  130 continue
      ipvt(n) = n
      if (abd(m,n) .eq. 0.0d0) info = n
      return
      end
      subroutine dgbsl(abd,lda,n,ml,mu,ipvt,b,job)

c*********************************************************************72
c
cc DGBSL solves a double precision band linear system.
c
c    The linear system has the form
c
c     a * x = b  or  trans(a) * x = b
c
c     The matrix has been factored by dgbco or dgbfa.
c
c     on entry
c
c        abd     double precision(lda, n)
c                the output from dgbco or dgbfa.
c
c        lda     integer
c                the leading dimension of the array  abd .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c
c        mu      integer
c                number of diagonals above the main diagonal.
c
c        ipvt    integer(n)
c                the pivot vector from dgbco or dgbfa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b , where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if dgbco has set rcond .gt. 0.0
c        or dgbfa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do j = 1, p
c              call dgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
c           end do
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c     fortran min0
c
      implicit double precision(a-h,o-z)

      integer lda,n,ml,mu,ipvt(1),job
      double precision abd(lda,1),b(1)

      double precision ddot,t
      integer k,kb,l,la,lb,lm,m,nm1

      m = mu + ml + 1
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve l*y = b
c
         if (ml .eq. 0) go to 30
         if (nm1 .lt. 1) go to 30
            do 20 k = 1, nm1
               lm = min0(ml,n-k)
               l = ipvt(k)
               t = b(l)
               if (l .eq. k) go to 10
                  b(l) = b(k)
                  b(k) = t
   10          continue
               call daxpy(lm,t,abd(m+1,k),1,b(k+1),1)
   20       continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/abd(m,k)
            lm = min0(k,m) - 1
            la = m - lm
            lb = k - lm
            t = -b(k)
            call daxpy(lm,t,abd(la,k),1,b(lb),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            lm = min0(k,m) - 1
            la = m - lm
            lb = k - lm
            t = ddot(lm,abd(la,k),1,b(lb),1)
            b(k) = (b(k) - t)/abd(m,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (ml .eq. 0) go to 90
         if (nm1 .lt. 1) go to 90
            do 80 kb = 1, nm1
               k = n - kb
               lm = min0(ml,n-k)
               b(k) = b(k) + ddot(lm,abd(m+1,k),1,b(k+1),1)
               l = ipvt(k)
               if (l .eq. k) go to 70
                  t = b(l)
                  b(l) = b(k)
                  b(k) = t
   70          continue
   80       continue
   90    continue
  100 continue
      return
      end
      subroutine digit_inc ( c )

c*********************************************************************72
c
cc DIGIT_INC increments a decimal digit.
c
c  Example:
c
c    Input  Output
c    -----  ------
c    '0'    '1'
c    '1'    '2'
c    ...
c    '8'    '9'
c    '9'    '0'
c    'A'    'A'
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character C, a digit to be incremented.
c
      implicit none

      character c
      integer digit

      call ch_to_digit ( c, digit )

      if ( digit .eq. -1 ) then
        return
      end if

      digit = digit + 1

      if ( digit .eq. 10 ) then
        digit = 0
      end if

      call digit_to_ch ( digit, c )

      return
      end
      subroutine digit_to_ch ( digit, c )

c*********************************************************************72
c
cc DIGIT_TO_CH returns the character representation of a decimal digit.
c
c  Example:
c
c    DIGIT   C
c    -----  ---
c      0    '0'
c      1    '1'
c    ...    ...
c      9    '9'
c     17    '*'
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIGIT, the digit value between 0 and 9.
c
c    Output, character C, the corresponding character, or '*' if DIGIT
c    was illegal.
c
      implicit none

      character c
      integer digit

      if ( 0 .le. digit .and. digit .le. 9 ) then

        c = char ( digit + 48 )

      else

        c = '*'

      end if

      return
      end
      subroutine dscal ( n, da, dx, incx )

c*********************************************************************72
c
cc DSCAL scales a vector by a constant.
c
c  Author:
c
c    Jack Dongarra
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, double precision DA, the scale factor.
c
c    Input/output, double precision DX(*), the vector to be scaled.
c
c    Input, integer INCX, the increment between successive entries.
c
      implicit none

      double precision da
      double precision dx(*)
      integer i
      integer incx
      integer m
      integer n
      integer nincx

      if ( n .le. 0 ) then
        return
      end if

      if ( incx .eq. 1 ) then

        m = mod ( n, 5 )

        do i = 1, m
          dx(i) = da * dx(i)
        end do

        do i = m+1, n, 5
          dx(i)   = da * dx(i)
          dx(i+1) = da * dx(i+1)
          dx(i+2) = da * dx(i+2)
          dx(i+3) = da * dx(i+3)
          dx(i+4) = da * dx(i+4)
        end do

      else

        nincx = n * incx
        do i = 1, nincx, incx
          dx(i) = da * dx(i)
        end do

      end if

      return
      end
      subroutine file_name_inc ( file_name )

c*********************************************************************72
c
cc FILE_NAME_INC generates the next filename in a series.
c
c  Discussion:
c
c    It is assumed that the digits in the name, whether scattered or
c    connected, represent a number that is to be increased by 1 on
c    each call.  If this number is all 9's on input, the output number
c    is all 0's.  Non-numeric letters of the name are unaffected, and
c    if the name contains no digits, then nothing is done.
c
c  Example:
c
c      Input          Output
c      -----          ------
c      a7to11.txt     a7to12.txt
c      a7to99.txt     a8to00.txt
c      a9to99.txt     a0to00.txt
c      cat.txt        cat.txt
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character ( len = * ) FILE_NAME.
c    On input, a character string to be incremented.
c    On output, the incremented string.
c
      implicit none

      character c
      logical ch_is_digit
      character*(*) file_name
      integer i
      integer lens

      lens = len ( file_name )

      do i = lens, 1, -1

        c = file_name(i:i)

        if ( ch_is_digit ( c ) ) then

          call digit_inc ( c )

          file_name(i:i) = c

          if ( c .ne. '0' ) then
            return
          end if

        end if

      end do

      return
      end
      subroutine grid_write ( me, mn, xn, yn, node, ne, nn )

c*********************************************************************72
c
cc GRID_WRITE writes the node and element information to a file.
c
c  Discussion:
c
c    The routine writes a file named "ELENODE.DAT".
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
c    Input, double precision XN(NN), YN(NN), the coordinates of the nodes.
c
c    Input, integer NODE(ME,6), the indices of the nodes that comprise each element.
c
c    Input, integer NE, the number of elements.
c
c    Input, integer NN, the number of nodes.
c
      implicit none

      integer me
      integer mn

      integer i
      integer j
      integer ne
      integer nn
      integer node(me,6)
      double precision xn(mn)
      double precision yn(mn)

      open ( unit = 10, file = 'elenode.dat', status = 'unknown' )

      write ( 10, * ) nn, ne

      do i = 1, nn
        write ( 10, * ) i, xn(i), yn(i)
      end do

      do i = 1, ne
        write ( 10, '(7i6)' ) i, ( node(i,j), j = 1, 6 )
      end do

      close ( unit = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GRID_WRITE - Wrote the file "elenode.dat".'

      return
      end
      function idamax ( n, dx, incx )

c*********************************************************************72
c
cc IDAMAX finds the vector element of largest magnitude.
c
c  Author:
c
c    Jack Dongarra
c
c  Parameters:
c
c    Input, integer N, the number of elements in the vector.
c
c    Input, double precision DX(*), the vector to examine.
c
c    Input, integer INCX, the increment between successive entries.
c
c    Output, integer IDAMAX, the index of the vector element of
c    largest magnitude.
c
      implicit none

      double precision dmax
      double precision dx(*)
      integer i
      integer idamax
      integer incx
      integer ix
      integer n

      idamax = 0
      if ( n < 1 ) then
        return
      end if

      idamax = 1
      if ( n == 1 ) then
        return
      end if

      if ( incx == 1 ) then

        dmax = dabs ( dx(1) )
        do i = 2, n
          if ( dmax < dabs ( dx(i) ) ) then
            idamax = i
            dmax = dabs ( dx(i) )
          end if
        end do

      else

        ix = 1
        dmax = dabs ( dx(1) )
        ix = ix + incx
        do i = 2, n
          if ( dmax < dabs ( dx(ix) ) ) then
            idamax = i
            dmax = dabs ( dx(ix) )
            ix = ix + incx
          end if
        end do

      end if

      return
      end
      subroutine nstoke(xc,yc,area,xm,ym,
     &  a,f,g,uold,visc_inv,tolns,xlngth,ylngth,
     &  node,indx,insc,ipivot,mrow1,
     &  nlband,nuband,nband,nrow1,ncol1,
     &  nelemn,np,nnodes,nuk,nquad,neqn1,
     &  nsteps,nsim,iwrite,maxnd,maxel,rdel,alpha)

c*********************************************************************72
c
cc NSTOKE solves the Navier-Stokes equations using Taylor-Hood elements.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 June 2002
c
c  Author:
c
c    Hyung-Chun Lee,
c    Department of Mathematics,
c    Ajou University, Korea
c
      implicit double precision(a-h,o-z)

      integer maxel
      integer maxnd
      integer mrow1

      double precision a(mrow1,*)
      double precision area(*)
      double precision f(*)
      double precision g(*)
      integer indx(maxnd,*)
      integer insc(*)
      integer ipivot(*)
      integer node(maxnd,*)
      double precision un(2)
      double precision uny(2)
      double precision unx(2)
      double precision uold(*)
      double precision xc(*)
      double precision xm(maxel,*)
      double precision yc(*)
      double precision ym(maxel,*)
c
c  zero arrays
c  g array contains the previous iterate
c  f array contains the right hand side initially and then the current
c  iterate is overwritten on f
c  T contains the Temperature for RHS
c
      visc=1.0d0/visc_inv
c
c  matrix assembly triangle by triangle
c  nsim is the number of simple iterations performed
c
      csim=0.0D+00
      do 700  iter=1, nsteps

        niter=iter
        if(iter.gt.nsim) csim=1.0D+00

        do i=1,nrow1
          do j=1,ncol1
            a(i,j)=0.0D+00
          end do
        end do

        do 240 it=1,nelemn
          arr=area(it)/3.0D+00
          do 230 iquad=1,nquad
            y=ym(it,iquad)
            x=xm(it,iquad)
            call trans(it,x,y,det,xix,xiy,etax,etay,xc,yc,node,maxnd)
            ar=arr*det

            do kk=1,2
              un(kk)=0.0D+00
              uny(kk)=0.0D+00
              unx(kk)=0.0D+00
            end do

            do iq=1,nnodes

              call refqbf(x,y,iq,bb,tbx,tby)
              bx=tbx*xix+tby*etax
              by=tbx*xiy+tby*etay
              ip=node(it,iq)

              do iuk=1,2
                iun=indx(ip,iuk)
                if(iun.eq.0) then
                  go to 140
                else if(iun.gt.0) then
                  un(iuk)=un(iuk)+bb*g(iun)
                  unx(iuk)=unx(iuk)+bx*g(iun)
                  uny(iuk)=uny(iuk)+by*g(iun)
                else
                  ubc=ubdry(iuk,ip,xc,yc,ylngth,alpha)
                  un(iuk)=un(iuk)+bb*ubc
                  unx(iuk)=unx(iuk)+bx*ubc
                  uny(iuk)=uny(iuk)+by*ubc
                end if
              end do

 140          continue

            end do

            do 220 iq=1,nnodes
              ip=node(it,iq)
              call refqbf(x,y,iq,bb,tbx,tby)
              bx=tbx*xix+tby*etax
              by=tbx*xiy+tby*etay
              bbl=refbsp(x,y,iq)
              do 210 iuk=1,nuk
                i=insc(ip)
                if(iuk.ne.3) i=indx(ip,iuk)
                if(i.le.0) go to 210

                if(iuk.eq.1) then
                  f(i)=f(i)+csim*( (un(1)*unx(1)+un(2)*uny(1))*bb )*ar
     &              +rdel*uold(i)*bb*ar
                else if(iuk.eq.2) then
                  f(i)=f(i)+csim*( (un(1)*unx(2)+un(2)*uny(2))*bb )*ar
     &              +rdel*uold(i)*bb*ar
                end if

                do 200 iqq=1,nnodes

                  ipp=node(it,iqq)
                  call refqbf(x,y,iqq,bbb,tbbx,tbby)
                  bbx=tbbx*xix+tbby*etax
                  bby=tbbx*xiy+tbby*etay
                  bbbl=refbsp(x,y,iqq)

                  do 190 iukk=1,nuk

                    j=insc(ipp)
                    if(iukk.ne.3) j=indx(ipp,iukk)
                    if(j.eq.0) go to 190
                    aij=0.0D+00
                    if(i.eq.neqn1) go to 190

                    if(iuk.eq.1) then

                      if(iukk.eq.1) then
                        aij=visc*(by*bby+bx*bbx)
     &                    +(bbb*unx(1)*bb)*csim
     &                    +bb*bbx*un(1)
     &                    +bb*bby*un(2) + rdel*(bb*bbb)
                      else if(iukk.eq.2) then
                        aij=csim*(bb*bbb*uny(1))
                      else if(iukk.eq.3) then
                        aij=-bx*bbbl
                      end if

                    else if(iuk.eq.2) then

                      if(iukk.eq.1) then
                        aij= csim*(bb*bbb*unx(2))
                      else if(iukk.eq.2) then
                        aij=(visc*(by*bby+bx*bbx)
     &                    +(bb*bbb*uny(2))*csim
     &                    +bb*bby*un(2)
     &                    +bb*bbx*un(1)) + rdel*(bb*bbb)
                      else if(iukk.eq.3) then
                        aij=-by*bbbl
                      end if

                    else

                      if(iukk.eq.1) then
                        aij=bbx*bbl
                      else if(iukk.eq.2) then
                        aij=bby*bbl
                      else
                        aij = 0.0
                      end if

                    end if

 180                continue

                    if(j.lt.0) go to 185
                    iuse=i-j+nband
                    a(iuse,j)=a(iuse,j)+aij*ar
                    go to 190
c
c  add terms to rhs for inhomogeneous boundary condition
c
 185                continue
                    f(i)=f(i)-ar*ubdry(iukk,ipp,xc,yc,ylngth,alpha)*aij
 190              continue
 200            continue
 210          continue
 220        continue
 230      continue
 240    continue

        f(neqn1)=0.0D+00
        do j=neqn1-nlband,neqn1-1
          i=neqn1-j+nband
          a(i,j)=0.0D+00
        end do
        a(nband,neqn1)=1.0D+00
c
c  solve system
c
        job=0
        call dgbfa(a,mrow1,neqn1,nlband,nuband,ipivot,info)
        call dgbsl(a,mrow1,neqn1,nlband,nuband,ipivot,f,job)

        if(info.eq.0) go to 509

 509    continue
c
c  check for convergence
c
        diff=0.0D+00
        do i=1,neqn1
          diff=diff+(g(i)-f(i))**2
        end do
        diff=sqrt(diff)
        write(6,1045) iter,diff

        if ( diff.le.tolns ) then
          go to 750
        end if

        do i=1,neqn1
          g(i)=f(i)
        end do

        do i=1,neqn1
          f(i)=0.0D+00
        end do

700   continue

750   continue

1045  format('  for iteration no.',i3,' difference in iterates is '
     &  ,2d14.8)
1050  format(/,' solutions converged in',i5,' iterations ',//)
1060  format(/,' solution calculated for reynolds number ',f8.2)

      return
      end
      function refbsp ( x, y, iq )

c*********************************************************************72
c
cc REFBSP evaluates a linear basis functions on the reference triangle.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 June 2002
c
c  Author:
c
c    Hyung-Chun Lee,
c    Department of Mathematics,
c    Ajou University, Korea
c
      implicit none

      double precision refbsp
      double precision x
      double precision y
      integer iq

      if ( iq .eq. 1 ) then
        refbsp = 1.0D+00 - x - y
      else  if ( iq .eq. 2 ) then
        refbsp = x
      else if ( iq .eq. 3 ) then
        refbsp = y
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'REFBSP - Fatal error!'
        write ( *, '(a)' ) '  Illegal input value of IQ = ', iq
        stop
      end if

      return
      end
      subroutine refqbf(x,y,in,bb,bx,by)

c*********************************************************************72
c
cc REFQBF evaluates quadratic basis functions on the reference triangle.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 June 2002
c
c  Author:
c
c    Hyung-Chun Lee,
c    Department of Mathematics,
c    Ajou University, Korea
c
      implicit none

      double precision bb
      double precision bx
      double precision by
      integer in
      double precision x
      double precision y

      if(in.eq.1) then
        bb=(1.0D+00-x-y)*(1.0D+00-2.0D+00*x-2.0D+00*y)
        bx=-3.0D+00+4.0D+00*x+4.0D+00*y
        by=-3.0D+00+4.0D+00*x+4.0D+00*y
      else if(in.eq.2) then
        bb=x*(2.d0*x-1.d0)
        bx=4.d0*x-1.d0
        by=0.d0
      else if(in.eq.3) then
        bb=y*(2.d0*y-1.d0)
        bx=0.d0
        by=4.d0*y-1.d0
      else if(in.eq.4) then
        bb=4.d0*x*(1.d0-x-y)
        bx=4.d0*(1.d0-2.d0*x-y)
        by=-4.d0*x
      else if(in.eq.5) then
        bb=4.d0*x*y
        bx=4.d0*y
        by=4.d0*x
      else if(in.eq.6) then
        bb=4.d0*y*(1.d0-x-y)
        bx=-4.d0*y
        by=4.d0*(1.d0-x-2.d0*y)
      end if

      return
      end
      subroutine setgrd(xc,yc,area,xm,ym,xlngth,ylngth,
     &  node,indx,insc,nlband,nband,
     &  nx,ny,nelemn,np,nnodes,nuk,nquad,neqn1,
     &  iwrite,maxnd,maxel)

c*********************************************************************72
c
cc SETGRD sets up the grid for the problem.
c
c  Discussion:
c
c    We are using quadratics for the velocity and linears for the
c    pressure
c
c    Input needed is nx,ny,xlngth,ylngth,write
c    Computes arrays node,area,xc,yc,xm,ym,index,insc
c    Computes neqn1, bandwidth information
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 May 2004
c
c  Author:
c
c    Hyung-Chun Lee,
c    Department of Mathematics,
c    Ajou University, Korea
c
c  Parameters:
c
c    Output, double precision XC(NP), YC(NP), the coordinates of the nodes.
c
c    Output, double precision AREA(NELEMN), the area of each element.
c
c    Output, double precision XM(MAXEL,3), YM(MAXEL,3), the coordinates
c    of quadrature points in each element.
c
c    Input, double precision XLNGTH, YLNGTH, the total width and height
c    of the region.
c
c    Output, integer NODE(MAXEL,6), the nodes that make up each element.
c
c    Output, integer INDX(MAXND,3), lists the indices of the U, V, and P
c    variables associated with the node.
c
c    Output, integer INSC(NP), is zero if a node is not a pressure node.
c    Otherwise, it is the index of the unknown pressure associated with the node.
c
c    Output, integer NLBAND, the half bandwidth for the finite element matrix.
c
c    Output, integer NBAND, the bandwidth for the finite element matrix.
c
c    Input, integer NX, NY, specifies the density of elements in
c    the X and Y directions.
c
c    Input, integer NELEMN, the number of elements.
c
c    Input, integer NP, the number of nodes.
c
c    Input, integer NNODES, the number of nodes per element.
c
c    Output, integer NUK, the maximum number of unknowns associated with one node.
c
c    Output, integer NQUAD, the number of quadrature points.
c
c    Output, integer NEQN1, the total number of unknowns.
c
c    Input, integer IWRITE, an output unit for debugging information.
c
c    Input, integer MAXND, the maximum number of nodes.
c
c    Input, integer MAXEL, the maximum number of elements.
c
      implicit double precision(a-h,o-z)

      integer maxel
      integer maxnd

      double precision area(*)
      integer indx(maxnd,*)
      integer insc(*)
      integer nnodes
      integer node(maxnd,*)
      integer nuk
      double precision xc(*)
      double precision xm(maxel,*)
      double precision yc(*)
      double precision ym(maxel,*)

c
c  Set parameters for Taylor-Hood element
c
      nnodes=6
      nuk=3
c
c  Construct grid
c  Coordinates and ordering of unknowns
c
      nym1=ny-1
      nxm1=nx-1
      nrow=nx+nxm1
      ncol=ny+nym1
      hx=xlngth/dble(nxm1)
      hx1=hx/2.0d0
      hy=ylngth/dble(nym1)
      hy1=hy/2.0d0

      i=0
      ip=0
      it1=-1
      xx=0.0D+00-hx1
      do 40 ic=1,nrow
        xx=xx+hx1
        icnt=mod(ic,2)
        do 35 jc=1,ncol
          jcnt=mod(jc,2)
          yy=hy1*dble(jc-1)
          ip=ip+1
          xc(ip)=xx
          yc(ip)=yy
          if(icnt.eq.1.and.jcnt.eq.1) go to 5
          go to 10
   5      if(jc.eq.ncol.or.ic.eq.nrow) go to 10
          it1=it1+2
          nelemn=it1+1
          ip1=ip+ncol
          ip2=ip+ncol+ncol
          node(it1,1)=ip
          node(it1,2)=ip2+2
          node(it1,3)=ip2
          node(it1,4)=ip1+1
          node(it1,5)=ip2+1
          node(it1,6)=ip1
          node(nelemn,1)=ip
          node(nelemn,2)=ip2+2
          node(nelemn,3)=ip+2
          node(nelemn,4)=ip1+1
          node(nelemn,5)=ip1+2
          node(nelemn,6)=ip+1
   10     continue
c
c Index for Dirichlet B.C.
c
          if(jc.eq.1.or.jc.eq.ncol) go to 20
          if(ic.eq.1) go to 20
          if(ic.eq.nrow.and.jc.le.33) go to 20
          i=i+2
          indx(ip,1)=i-1
          indx(ip,2)=i
          go to 29
   20     continue
          indx(ip,1)=-1
          indx(ip,2)=-1
   29     continue
          if(jcnt.eq.0.or.icnt.eq.0) go to 30
          i=i+1
          insc(ip)=i
          go to 35
   30     insc(ip)=0
   35   continue
   40 continue

      np=ip
      neqn1=i
      iw=2
      if (iw.eq.1) then
        write(iwrite,1089) (i,xc(i),yc(i),indx(i,1),
     &                indx(i,2),insc(i),i=1,np)
        write(iwrite,1099) (it,(node(it,i),i=1,6),it=1,nelemn)
      end if

c  set quadrature information
c  quadrature rule is midpoint
c  coordinates and weights are set by routine QD7PT
c
      nquad=3
      do it=1,nelemn
        xm(it,1)=0.5D+00
        xm(it,2)=0.5D+00
        xm(it,3)=0.0D+00
        ym(it,1)=0.0D+00
        ym(it,2)=0.5D+00
        ym(it,3)=0.5D+00
        area(it)=0.5D+00
      end do
c
c  half band width
c
      nlband=0
      do 90 it=1,nelemn
        do 90 iq=1,nnodes
          ip=node(it,iq)
          do 80 iuk=1,nuk
            if(iuk.eq.3) then
              i=insc(ip)
            else
              i=indx(ip,iuk)
            end if
            if(i.le.0) go to 80
            do 70 iqq=1,nnodes
              ipp=node(it,iqq)
              do 70 iukk=1,nuk
                if(iukk.eq.nuk) then
                  j=insc(ipp)
                else
                  j=indx(ipp,iukk)
                end if
                if(i.gt.j) go to 70
                ij=j-i
                if(ij.gt.nlband) nlband=ij
   70       continue
   80     continue
   90 continue
      nband=nlband+nlband+1

 1099 format(7i6)
 1089 format(i5,2f12.5,3i5)
 2001 format(' unknown numbers along line:', 7i5)

      return
      end
      subroutine solution_init ( neqn1, uold )

c*********************************************************************72
c
cc SOLUTION_INIT returns an initial value for the solution coefficient vector.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 April 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NEQN1, the number of coefficients.
c
c    Output, double precision UOLD(NEQN1), the initial value for the
c    solution coefficients.
c
      implicit none

      integer neqn1

      integer i
      double precision uold(neqn1)

      open ( unit = 10, file = 'fuini.dat' )

      do i = 1, neqn1
        read ( 10, * ) uold(i)
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
      subroutine trans(it,xq,yq,det,pj11,pj21,pj12,pj22,
     &  xc,yc,node,maxnd)

c*********************************************************************72
c
cc TRANS transforms data between the reference and physical elements.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 June 2002
c
c  Author:
c
c    Hyung-Chun Lee,
c    Department of Mathematics,
c    Ajou University, Korea
c
      implicit none

      integer maxnd

      double precision det
      double precision f1x
      double precision f1y
      double precision f2x
      double precision f2y
      integer i1
      integer i2
      integer i3
      integer i4
      integer i5
      integer i6
      integer it
      integer node(maxnd,*)
      double precision pj11
      double precision pj12
      double precision pj21
      double precision pj22
      double precision x1
      double precision x2
      double precision x3
      double precision x4
      double precision x5
      double precision x6
      double precision xc(*)
      double precision xq
      double precision y1
      double precision y2
      double precision y3
      double precision y4
      double precision y5
      double precision y6
      double precision yc(*)
      double precision yq

      i1=node(it,1)
      i2=node(it,2)
      i3=node(it,3)
      i4=node(it,4)
      i5=node(it,5)
      i6=node(it,6)
      x1=xc(i1)
      y1=yc(i1)
      x2=xc(i2)
      y2=yc(i2)
      x3=xc(i3)
      y3=yc(i3)
      x4=xc(i4)
      y4=yc(i4)
      x5=xc(i5)
      y5=yc(i5)
      x6=xc(i6)
      y6=yc(i6)
c
c  Compute partial derivatives at point (xq,yq)
c
      f1x=x1*(-3.0D+00+4.0D+00*xq+4.0D+00*yq)
     &  +x2*(4.0D+00*xq-1.0D+00)
     &  +x4*4.0D+00*(1.0D+00-2.0D+00*xq-yq)
     &  +x5*4.0D+00*yq + x6*4.0D+00*(-yq)

      f1y=x1*(-3.d0+4.d0*xq+4.d0*yq)
     &  +x3*(4.d0*yq-1.d0)
     &  +x4*4.d0*(-xq) + x5*4.d0*xq
     &  +x6*4.d0*(1.d0-xq-2.d0*yq)

      f2x=y1*(-3.d0+4.d0*xq+4.d0*yq)
     &  +y2*(4.d0*xq-1.d0)
     &  +y4*4.d0*(1.d0-2.d0*xq-yq)
     &  +y5*4.d0*yq + y6*4.d0*(-yq)

      f2y=y1*(-3.d0+4.d0*xq+4.d0*yq)
     &  +y3*(4.d0*yq-1.d0)
     &  +y4*4.d0*(-xq) + y5*4.d0*xq
     &  +y6*4.d0*(1.d0-xq-2.d0*yq)
c
c  Compute determinant of transformation evaluated at point (xq,yq)
c
      det = f1x*f2y-f1y*f2x
c
c  Compute j11, j22, j21, j22
c
      pj11 = f2y/det
      pj12 = -f2x/det
      pj21 = -f1y/det
      pj22 = f1x/det

      det = dabs(det)

      return
      end
      function ubdry ( iuk, ip, xc, yc, ylngth, alpha )

c*********************************************************************72
c
cc UBDRY evaluates the boundary conditions.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 June 2002
c
c  Author:
c
c    Hyung-Chun Lee,
c    Department of Mathematics,
c    Ajou University, Korea
c
      implicit none

      double precision alpha
      integer ip
      integer iuk
      double precision ubdry
      double precision x
      double precision xc(*)
      double precision y
      double precision yc(*)
      double precision ylngth

      x = xc(ip)
      y = yc(ip)
      ubdry = 0.0D+00

      if ( x .eq. 0.0D+00 .and. y .le. 0.2D+00 ) then

        if ( iuk .eq. 1 ) then
          ubdry = alpha * 10.0D+00 * y * ( 0.2D+00 - y )
        end if

      end if

      return
      end
