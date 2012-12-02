      program main

c*********************************************************************72
c  
cc MAIN is the main program for TCELL.
c
c  Discussion:
c
c    TCELL solves the steady incompressible Navier-Stokes equations
c    in a 2D "T"-shaped region.
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
c    30 January 2004
c
c  Author:
c
c    Hyung-Chun Lee,
c    Department of Mathematics,
c    Ajou University, Korea
c
c  Local parameters:
c
c    BC_TYPE selects the boundary conditions, by controlling the value of ALPHA.
c    Legal values are 1, for a step function, 2 for a "hat" function, 3 for a sinusoid.
c
c    FILENAME is the name to be used for the first output file.
c    Each subsequent output file created by the program will have an
c    incremented name.
c
      implicit none

      integer nx
      integer ny
      integer mx
      integer my
      integer maxel
      integer maxnd
      integer maxun
      integer minun

      parameter(nx=41)
      parameter(ny=41)
      parameter(mx=2*nx-1)
      parameter(my=2*ny-1)
      parameter(maxel=  2*(nx-1)*(ny-1) )
      parameter(maxnd=  mx*my)
      parameter(maxun=  2*mx*my+nx*ny)
      parameter(minun=  27*ny)

      double precision a(minun,maxun)
      double precision alpha
      double precision area(maxel)
      integer bc_type
      double precision deltat
      double precision f(maxun)
      character*20 filename
      double precision g(maxun)
      double precision gg(maxun)
      double precision gp(maxun)
      double precision gx
      double precision gy
      integer i
      integer ic
      integer indx(maxnd,2)
      integer insc(maxnd)
      integer ipivot(maxun)
      integer iter
      integer iukp
      integer iuku
      integer iukv
      integer iwrite
      integer mrow1
      integer nband
      integer ncol1
      integer nelemn
      integer neqn1
      integer nlband
      integer nnodes
      integer node(maxel,6)
      integer np
      integer nquad
      integer nrow1
      integer nsim
      integer nsteps
      integer nuband
      integer nuk
      integer ny2
      double precision pi
      double precision pp
      double precision rdel
      double precision reynld
      double precision tolns
      double precision tolopt
      double precision ubdry
      double precision uold(maxun)
      double precision xc(maxnd)
      double precision xlngth
      double precision xm(maxel,3)
      double precision yc(maxnd)
      double precision ylngth
      double precision ym(maxel,3)

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TCELL:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Solve the Navier Stokes fluid flow'
      write ( *, '(a)' ) '  equations in a TCELL region,'
      write ( *, '(a)' ) '  using finite elements.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Maximum number of nodes =    ', maxnd
      write ( *, '(a,i6)' ) '  Maximum number of elements = ', maxel

      bc_type = 1
      filename = 'up000.txt'
c
c  Set functions for the input profile.
c
      iwrite = 0
c
c   input
c
      xlngth = 1.0d0
      ylngth = 1.0d0
      reynld = 1.0d0
      mrow1 = minun
      nsim = 3
      nsteps = 10
      tolns = 1.0d-6
      tolopt = 1.0d-06
      pi = 4.0d0 * datan ( 1.0D+00 )
c
c  2.  setgrd constructs grid, numbers unknowns, calculates areas,
c      and points for midpoint quadrature rule
c      SETGRD sets bandwidth and NEQN1
c
      call setgrd(xc,yc,area,xm,ym,xlngth,ylngth,
     &    node,indx,insc,nlband,nband,
     &    nx,ny,nelemn,np,nnodes,nuk,nquad,neqn1,
     &    iwrite,maxnd,maxel)

      write ( *, '(a,i6)' ) '  Number of nodes =    ', np
      write ( *, '(a,i6)' ) '  Number of elements = ', nelemn
      write ( *, '(a,i6)' ) '  Lower bandwidth =    ', nlband
      write ( *, '(a,i6)' ) '  Total bandwidth =    ', nband

      call grid_write ( maxel, maxnd, xc, yc, node, nelemn, np )

      nuband = nlband
      nrow1 = nlband+nlband+nuband+1
      ncol1 = neqn1
      ny2 = ny+ny-1

      deltat = 0.0002d0
      rdel=1.0d0 / deltat

      do i = 1, neqn1
        f(i) = 0.0d0
      end do
c
c  Initialize the solution.
c
      call solution_init ( neqn1, uold )
c
c  Carry out the time iteration.
c
      do iter = 1, 500

        if ( bc_type == 1 ) then

          if ( iter .le. 250 ) then
            alpha = 5.0d0
          else 
            alpha = 1.0d0
          endif

        else if ( bc_type == 2 ) then

          if ( iter .le. 250 ) then
            alpha = 80.0d0 * dble(iter)*deltat + 1.0d0
          else 
            alpha = -80.0d0 * dble(iter)*deltat + 9.0d0
          endif

        else if ( bc_type == 3 ) then

          alpha = 2.0d0 * sin ( dble(iter) * 0.01d0 * pi )   

        end if

        do i = 1, neqn1
          g(i) = f(i)
        end do

        do i = 1, neqn1
          f(i) = 0.0d0
        end do

        call nstoke ( xc, yc, area, xm, ym,
     &     a, f, g, uold, reynld, tolns, xlngth, ylngth,
     &     node, indx, insc, ipivot, mrow1,
     &     nlband, nuband, nband, nrow1, ncol1,
     &     nelemn, np, nnodes, nuk, nquad, neqn1,
     &     nsteps, nsim, iwrite, maxnd, maxel, rdel, alpha )
c
c  Save u=(gx,gy) to 'ue.dat' for 't.f'
c
        do i = 1, neqn1
          uold(i) = f(i)
        end do
c
c  Increment the filename, and 
c  save u=(gx,gy) to 'up???.dat'
c
        call file_name_inc ( filename )

        write ( *, * ) 'Creating ', filename

        open ( unit = 1, file = filename )

        do ic = 1, np

          iuku=indx(ic,1)
          iukv=indx(ic,2)
          iukp=insc(ic)

          if ( iuku .eq. 0 ) then
            gx=0.0d0
            if ( iukp .eq. 0 ) then
              pp=0.0d0
            else
              pp=f(iukp)
            end if
            go to 54
          else if ( iuku .eq. -1 ) then
            gx = alpha * ubdry ( 1, ic, xc, yc )
            if ( iukp .eq. 0 ) then
              pp=0.0d0
            else
              pp=f(iukp)
            end if
            go to 54
          else
            gx=f(iuku)
          end if

          if ( iukp .eq. 0 ) then
            pp=0.0d0
          else
            pp=f(iukp)
          end if

54        continue

          if ( iukv .le. 0 ) then
            gy = 0.0d0
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
          write ( 1, '(e25.15,2x,e25.15)' ) gx, gy

        end do

        close (unit=1)

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TCELL:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

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
c
      DOUBLE PRECISION DX(1),DY(1),DA
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
c
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
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
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

      do I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I + 1) = DY(I + 1) + DA*DX(I + 1)
        DY(I + 2) = DY(I + 2) + DA*DX(I + 2)
        DY(I + 3) = DY(I + 3) + DA*DX(I + 3)
      end do

      return
      end
      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)

c*********************************************************************72
c
cc DDOT forms the dot product of two vectors.
c
c     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
c     JACK DONGARRA, LINPACK, 3/11/78.
c
      implicit double precision(a-h,o-z)

      DOUBLE PRECISION DX(1),DY(1),DTEMP
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N

      DDOT = 0.0D0
      DTEMP = 0.0D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
c
c    CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
c    NOT EQUAL TO 1
c
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = DTEMP + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      DDOT = DTEMP
      RETURN
c
c    CODE FOR BOTH INCREMENTS EQUAL TO 1
c
c
c    CLEAN-UP LOOP
c
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + DX(I)*DY(I)
   30 CONTINUE
      IF( N .LT. 5 ) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
      DTEMP = DTEMP + DX(I)*DY(I) + DX(I + 1)*DY(I + 1) +
     &  DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
   50 CONTINUE
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
c                      do 10 i = i1, i2
c                         k = i - j + m
c                         abd(k,j) = a(i,j)
c                10    continue
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
      do 20 jz = j0, j1
         i0 = m + 1 - jz
         do 10 i = i0, ml
            abd(i,jz) = 0.0d0
   10    continue
   20 continue
   30 continue
      jz = j1
      ju = 0
c
c     gaussian elimination with partial pivoting
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
         if (ml .lt. 1) go to 50
            do 40 i = 1, ml
               abd(i,jz) = 0.0d0
   40       continue
   50    continue
c
c        find l = pivot index
c
         lm = min0(ml,n-k)
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
c           do 10 j = 1, p
c              call dgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
c        10 continue
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
      SUBROUTINE DSCAL(N,DA,DX,INCX)

c*********************************************************************72
c
cc DSCAL scales a vector by a constant.
c
c     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.
c     JACK DONGARRA, LINPACK, 3/11/78.
c
      implicit double precision(a-h,o-z)
c
      DOUBLE PRECISION DA,DX(1)
      INTEGER I,INCX,M,MP1,N,NINCX
c
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GO TO 20
c
c   CODE FOR INCREMENT NOT EQUAL TO 1
c
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
      DX(I) = DA*DX(I)
   10 CONTINUE
      RETURN
c
c   CODE FOR INCREMENT EQUAL TO 1
c
c
c   CLEAN-UP LOOP
c
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DX(I) = DA*DX(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DX(I) = DA*DX(I)
        DX(I + 1) = DA*DX(I + 1)
        DX(I + 2) = DA*DX(I + 2)
        DX(I + 3) = DA*DX(I + 3)
        DX(I + 4) = DA*DX(I + 4)
   50 CONTINUE
      RETURN
      END
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
      INTEGER FUNCTION IDAMAX(N,DX,INCX)

c*********************************************************************72
c
cc IDAMAX finds the vector element of largest magnitude.
c
c     JACK DONGARRA, LINPACK, 3/11/78.
c
      implicit double precision(a-h,o-z)

      DOUBLE PRECISION DX(1),DMAX
      INTEGER I,INCX,IX,N

      IDAMAX = 0
      IF( N .LT. 1 ) RETURN
      IDAMAX = 1
      IF(N.EQ.1)RETURN
      IF(INCX.EQ.1)GO TO 20
c
c   CODE FOR INCREMENT NOT EQUAL TO 1
c
      IX = 1
      DMAX = DABS(DX(1))
      IX = IX + INCX
      DO 10 I = 2,N
       IF(DABS(DX(IX)).LE.DMAX) GO TO 5
        IDAMAX = I
        DMAX = DABS(DX(IX))
    5   IX = IX + INCX
   10 CONTINUE
      RETURN
c
c   CODE FOR INCREMENT EQUAL TO 1
c
   20 DMAX = DABS(DX(1))
      DO 30 I = 2,N
        IF(DABS(DX(I)).LE.DMAX) GO TO 30
        IDAMAX = I
        DMAX = DABS(DX(I))
   30 CONTINUE
      RETURN
      END
      subroutine nstoke ( xc, yc, area, xm, ym,
     &  a, f, g, uold, reynld, tolns, xlngth, ylngth,
     &  node, indx, insc, ipivot, mrow1,
     &  nlband, nuband, nband, nrow1, ncol1,
     &  nelemn, np, nnodes, nuk, nquad, neqn1,
     &  nsteps, nsim, iwrite, maxnd, maxel, rdel, alpha )

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
c    12 May 2005
c
c  Author:
c
c    Hyung-Chun Lee,
c    Department of Mathematics,
c    Ajou University, Korea
c
c  Parameters:
c
      implicit none

      integer maxel
      integer maxnd
      integer mrow1

      double precision a(mrow1,*)
      double precision aij
      double precision alpha
      double precision ar
      double precision area(*)
      double precision arr
      double precision bb
      double precision bbb
      double precision bbbl
      double precision bbl
      double precision bbx
      double precision bby
      double precision bx
      double precision by
      double precision csim
      double precision det
      double precision diff
      double precision etax
      double precision etay
      double precision f(*)
      double precision g(*)
      integer i
      integer indx(maxnd,*)
      integer info
      integer insc(*)
      integer ip
      integer ipivot(*)
      integer ipp
      integer iq
      integer iqq
      integer iquad
      integer it
      integer iter
      integer iuk
      integer iukk
      integer iun
      integer iuse
      integer iwrite
      integer j
      integer job
      integer kk
      integer nband
      integer ncol1
      integer nelemn
      integer neqn1
      integer nlband
      integer nnodes
      integer node(maxel,*)
      integer np
      integer nquad
      integer nrow1
      integer nsim
      integer nsteps
      integer nuband
      integer nuk
      double precision rdel
      double precision refbsp
      double precision reynld
      double precision tbbx
      double precision tbby
      double precision tbx
      double precision tby
      double precision tolns
      double precision ubc
      double precision ubdry
      double precision un(2)
      double precision unx(2)
      double precision uny(2)
      double precision uold(*)
      double precision uold_qp
      double precision visc
      double precision vold_qp
      double precision x
      double precision xc(*)
      double precision xix
      double precision xiy
      double precision xlngth
      double precision xm(maxel,*)
      double precision y
      double precision yc(*)
      double precision ylngth
      double precision ym(maxel,*)
c
c  zero arrays
c  g array contains the previous iterate
c  f array contains the right hand side initially and then the current
c  iterate is overwritten on f
c
      visc = 1.0D+00 / reynld
c
c  matrix assembly triangle by triangle
c  nsim is the number of simple iterations performed
c
      do iter = 1, nsteps

        if ( iter .le. nsim ) then
          csim = 0.0D+00
        else
          csim = 1.0D+00
        end if

        do i = 1, nrow1
          do j = 1, ncol1
            a(i,j) = 0.0D+00
          end do
        end do

        do it = 1, nelemn

          arr = area(it) / 3.0D+00

          do iquad = 1, nquad

            y = ym(it,iquad)
            x = xm(it,iquad)
            call trans(it,x,y,det,xix,xiy,etax,etay,xc,yc,node,maxel)
            ar = arr * det

            do kk = 1, 2
              un(kk)  = 0.0D+00
              uny(kk) = 0.0D+00
              unx(kk) = 0.0D+00
            end do
            uold_qp = 0.0D+00
            vold_qp = 0.0D+00

            do iq = 1, nnodes

              call refqbf ( x, y, iq, bb, tbx, tby )

              bx = tbx * xix + tby * etax
              by = tbx * xiy + tby * etay
              ip = node(it,iq)

              do iuk = 1, 2

                iun = indx(ip,iuk)

                if ( 0 .lt. iun ) then

                  un(iuk)  = un(iuk)  + bb * g(iun)
                  unx(iuk) = unx(iuk) + bx * g(iun)
                  uny(iuk) = uny(iuk) + by * g(iun)

                  if ( iuk .eq. 1 ) then
                    uold_qp = uold_qp + bb * uold(iun)
                  else if ( iuk .eq. 2 ) then
                    vold_qp = vold_qp + bb * uold(iun)
                  end if

                else if ( iun .lt. 0 ) then

                  ubc = alpha * ubdry ( iuk, ip, xc, yc )

                  un(iuk)  = un(iuk)  + bb * ubc
                  unx(iuk) = unx(iuk) + bx * ubc
                  uny(iuk) = uny(iuk) + by * ubc

                  if ( iuk .eq. 1 ) then
                    uold_qp = uold_qp + bb * ubc
                  else if ( iuk .eq. 2 ) then
                    vold_qp = vold_qp + bb * ubc
                  end if

                end if

              end do

            end do
   
            do iq = 1, nnodes

              ip = node(it,iq)
              call refqbf ( x, y, iq, bb, tbx, tby )
              bx = tbx * xix + tby * etax
              by = tbx * xiy + tby * etay
              bbl = refbsp ( x, y, iq )

              do 210 iuk = 1, nuk

                if ( iuk .eq. 3 ) then
                  i = insc(ip)
                else
                  i = indx(ip,iuk)
                end if

                if(i.le.0) go to 210

                if ( iuk .eq. 1 ) then
                  f(i) = f(i) 
     &                + csim *( (un(1)*unx(1)+un(2)*uny(1))*bb ) * ar
     &                + rdel * uold_qp * bb*ar
                else if ( iuk .eq. 2 ) then
                  f(i)=f(i)
     &               + csim *( (un(1)*unx(2)+un(2)*uny(2))*bb ) * ar
     &               + rdel * vold_qp * bb*ar
                end if

                do iqq = 1, nnodes

                  ipp = node(it,iqq)
                  call refqbf ( x, y, iqq, bbb, tbbx, tbby )
                  bbx = tbbx * xix + tbby * etax
                  bby = tbbx * xiy + tbby * etay
                  bbbl = refbsp ( x, y, iqq )

                  do 190 iukk = 1, nuk

                    if ( iukk < 3 ) then
                      j = indx(ipp,iukk)
                    else
                      j = insc(ipp)
                    end if

                    if ( j .eq. 0 ) go to 190
                    aij = 0.0D+00
                    if ( i .eq. neqn1 ) go to 190

                    if ( iuk .eq. 1 ) then

                      if ( iukk .eq. 1 ) then
                        aij = visc * (by*bby+bx*bbx)
     &                     + (bbb*unx(1)*bb)*csim
     &                     + bb*bbx*un(1)
     &                     + bb*bby*un(2) + rdel*(bb*bbb)
                      else if ( iukk .eq. 2 ) then
                        aij = csim * bb * bbb * uny(1)
                      else if ( iukk.eq.3 ) then
                        aij = -bx * bbbl
                      end if

                    else if ( iuk .eq. 2 ) then

                      if ( iukk .eq. 1) then
                        aij = csim * bb * bbb * unx(2)
                      else if ( iukk .eq. 2 ) then
                        aij = (visc*(by*bby+bx*bbx)
     &                     + (bb*bbb*uny(2))*csim
     &                     + bb*bby*un(2)
     &                     + bb*bbx*un(1)) + rdel * bb * bbb
                      else if ( iukk .eq. 3 ) then
                        aij = -by * bbbl
                      end if

                    else if ( iuk .eq. 3 ) then

                      if ( iukk .eq. 1 ) then
                        aij = bbx * bbl
                      else if ( iukk .eq. 2 ) then
                        aij = bby * bbl
                      end if

                    end if
c
c  An inhomogeneous boundary condition results in a term added
c  to the right hand side.
c
                    if ( j .lt. 0 ) then
                      ubc = alpha * ubdry ( iukk, ipp, xc, yc )
                      f(i) = f(i) - ar * aij * ubc
                    else
                      iuse = i - j + nband
                      a(iuse,j) = a(iuse,j) + ar * aij
                    end if

 190              continue
                end do
 210          continue
            end do
          end do
        end do
c
c  Replace the last equation by an equation that forces the pressure to be zero
c  at the last pressure node.
c
        f(neqn1) = 0.0D+00
        do j = neqn1-nlband, neqn1-1
          i = neqn1 - j + nband
          a(i,j) = 0.0D+00
        end do
        a(nband,neqn1) = 1.0D+00
c
c  Factor the matrix.
c
        call dgbfa ( a, mrow1, neqn1, nlband, nuband, ipivot, info )

        if ( info .ne. 0 ) then
          write ( *, * ) ' '
          write ( *, * ) 'NSTOKE - Fatal error!'
          write ( *, * ) '  DGBFA returned INFO = ', info
          stop
        end if
c
c  Solve the system.
c
        job = 0
        call dgbsl ( a, mrow1, neqn1, nlband, nuband, ipivot, f, job )
c
c  Check for convergence
c
        diff = 0.0D+00
        do i = 1, neqn1
          diff = diff + ( g(i) - f(i) )**2
        end do
        diff = sqrt ( diff )
        write(6,1045) iter, diff

        if ( diff .le. tolns ) then
          go to 750
        end if

        do i = 1, neqn1
          g(i) = f(i)
          f(i) = 0.0D+00
        end do

      end do

750    continue

1045  format('  for iteration no.',i3,' difference in iterates is '
     &      ,2d14.8)

      return
      end
      function refbsp ( x, y, iq )

c*********************************************************************72
c
cc REFBSP evaluates a linear basis function on the reference triangle.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 April 2005
c
c  Author:
c
c    Hyung-Chun Lee,
c    Department of Mathematics,
c    Ajou University, Korea
c
c  Parameters:
c
      implicit none

      integer iq
      double precision refbsp
      double precision x
      double precision y

      if ( iq == 1 ) then
        refbsp = 1.0D+00 - x - y
      else if ( iq == 2 ) then
        refbsp = x
      else if ( iq == 3 ) then
        refbsp = y
      else
        refbsp = 0.0D+00
      end if

      return
      end
      subroutine refqbf ( x, y, in, bb, bx, by )

c*********************************************************************72
c
cc REFQBF evaluates a quadratic basis function on the reference triangle.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 April 2005
c
c  Author:
c
c    Hyung-Chun Lee,
c    Department of Mathematics,
c    Ajou University, Korea
c
c  Parameters:
c
c    Input, double precision X, Y, coordinates in the reference triangle
c    of the point where the basis function is to be evaluated.
c
c    Input, integer IN, the index of the basis function.
c
c    Output, double precision BB, BX, BY, the value of the basis function,
c    and its X and Y derivatives, at the given point.
c
      implicit none

      double precision bb
      double precision bx
      double precision by
      integer in
      double precision x
      double precision y

      if (in.eq.1) then
        bb=(1.0d0-x-y)*(1.0d0-2.0d0*x-2.0d0*y)
        bx=-3.0d0+4.0d0*x+4.0d0*y
        by=-3.0d0+4.0d0*x+4.0d0*y
      else if (in.eq.2) then
        bb=x*(2.0d0*x-1.0d0)
        bx=4.0d0*x-1.0d0
        by=0.0d0
      else if (in.eq.3) then
        bb=y*(2.0d0*y-1.0d0)
        bx=0.0d0
        by=4.0d0*y-1.0d0
      else if (in.eq.4) then
        bb=4.0d0*x*(1.0d0-x-y)
        bx=4.0d0*(1.0d0-2.0d0*x-y)
        by=-4.0d0*x
      else if (in.eq.5) then
        bb=4.0d0*x*y
        bx=4.0d0*y
        by=4.0d0*x
      else if (in.eq.6) then
        bb=4.0d0 * y * ( 1.0d0 - x - y )
        bx=-4.0d0*y
        by=4.0d0*(1.0d0-x-2.0d0*y)
      else
        bb = 0.0D+00
        bx = 0.0D+00
        by = 0.0D+00
      end if

      return
      end
      subroutine setgrd ( xc, yc, area, xm, ym, xlngth, ylngth,
     &   node, indx, insc, nlband, nband, nx, ny, nelemn, np, 
     &   nnodes, nuk, nquad, neqn1, iwrite, maxnd, maxel )

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
c    13 April 2004
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
      implicit none

      integer maxel
      integer maxnd
      integer np

      double precision area(*)
      double precision hx
      double precision hx1
      double precision hy
      double precision hy1
      integer i
      integer ic
      integer icnt
      integer ij
      integer indx(maxnd,*)
      integer insc(*)
      integer ip
      integer ip1
      integer ip2
      integer ipp
      integer iq
      integer iqq
      integer iquater1
      integer iquater2
      integer it
      integer it1
      integer iuk
      integer iukk
      integer iw
      integer iwrite
      integer j
      integer jc
      integer jcnt
      integer nband
      integer ncol
      integer nelemn
      integer neqn1
      integer nlband
      integer nnodes
      integer node(maxel,*)
      integer nquad
      integer nrow
      integer nuk
      integer nx
      integer nxm1
      integer ny
      integer nym1
      double precision xc(np)
      double precision xlngth
      double precision xm(maxel,*)
      double precision xx
      double precision yc(np)
      double precision ylngth
      double precision ym(maxel,*)
      double precision yy
c
c  Set parameters for the Taylor-Hood element.
c
      nnodes=6
      nuk=3
c
c  construct grid
c  coordinates and ordering of unknowns
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
      xx=0.0d0-hx1

      iquater1=(nrow-1)/4-2
      iquater2=3*(nrow-1)/4

      do ic=1,iquater1

        xx=xx+hx1
        icnt=mod(ic,2)

        do jc=1,ny

          jcnt=mod(jc,2)
          yy=0.5d0+hy1*dble(jc-1)
          ip=ip+1
          xc(ip)=xx
          yc(ip)=yy

          if ( icnt .eq. 1 .and. jcnt .eq. 1 .and. jc < ny ) then
  
            it1=it1+2
            nelemn=it1+1
            ip1=ip+ny
            ip2=ip+ny+ny
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

          end if
c
c  Index for Dirichlet B.C.
c
          if ( jc .ne. 1 .and. jc .ne. ny .and. ic .ne. 1 ) then
            i=i+2
            indx(ip,1)=i-1
            indx(ip,2)=i
          else
            indx(ip,1)=-1
            indx(ip,2)=-1
          end if

          if ( jcnt .eq. 1 .and. icnt .eq. 1 ) then
            i=i+1
            insc(ip)=i
          else
            insc(ip)=0
          end if

        end do
      end do

      do 140 ic=iquater1+1,iquater1+2
        xx=xx+hx1
        icnt=mod(ic,2)
        do 135 jc=1,ny
          jcnt=mod(jc,2)
          yy=0.5d0+hy1*dble(jc-1)
          ip=ip+1
          xc(ip)=xx
          yc(ip)=yy
          if(icnt.eq.1.and.jcnt.eq.1) go to 15
          go to 110
   15     if(jc.eq.ny) go to 110
          it1=it1+2
          nelemn=it1+1
          ip1=ip+ny
          ip2=ip+ny+ncol
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
  110     continue
c
c Index for Dirichlet B.C.
c
        if(jc.eq.1.or.jc.eq.ny) go to 120
        i=i+2
          indx(ip,1)=i-1
          indx(ip,2)=i
          go to 129
  120     continue
          indx(ip,1)=-1
          indx(ip,2)=-1
  129     continue
          if(jcnt.eq.0.or.icnt.eq.0) go to 130
          i=i+1
          insc(ip)=i
          go to 135
  130     insc(ip)=0
  135   continue
  140 continue

      do 240 ic=iquater1+3,iquater2
        xx=xx+hx1
        icnt=mod(ic,2)
        do 235 jc=1,ncol
          jcnt=mod(jc,2)
          yy=hy1*dble(jc-1)
          ip=ip+1
          xc(ip)=xx
          yc(ip)=yy
          if(icnt.eq.1.and.jcnt.eq.1) go to 25
          go to 210
   25     if(jc.eq.ncol) go to 210
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
  210     continue
c
c Index for Dirichlet B.C.
c
        if(jc.eq.1.or.jc.eq.ncol) go to 220
        if(ic.eq.iquater1+3.and.jc.le.ny) go to 220
        i=i+2
          indx(ip,1)=i-1
          indx(ip,2)=i
          go to 229
  220     continue
          indx(ip,1)=-1
          indx(ip,2)=-1
  229     continue
          if(jcnt.eq.0.or.icnt.eq.0) go to 230
          i=i+1
          insc(ip)=i
          go to 235
  230     insc(ip)=0
  235   continue
  240 continue

        xx=xx+hx1
        icnt=mod(iquater2+1,2)
        do 255 jc=1,ny-1
          jcnt=mod(jc,2)
          yy=hy1*dble(jc-1)
          ip=ip+1
          xc(ip)=xx
          yc(ip)=yy
          indx(ip,1)=-1
          indx(ip,2)=-1
          if(jcnt.eq.0.or.icnt.eq.0) go to 930
          i=i+1
          insc(ip)=i
          go to 255
  930     insc(ip)=0
  255   continue

       xx=xx-hx1
      do 340 ic=iquater2+1,nrow
        xx=xx+hx1
        icnt=mod(ic,2)
        do 335 jc=1,ny
          jcnt=mod(jc,2)
          yy=0.5d0+hy1*dble(jc-1)
          ip=ip+1
          xc(ip)=xx
          yc(ip)=yy
          if(icnt.eq.1.and.jcnt.eq.1) go to 435
          go to 310
  435     if(ic.eq.nrow.or.jc.eq.ny) go to 310
          it1=it1+2
          nelemn=it1+1
          ip1=ip+ny
          ip2=ip+ny+ny
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
  310     continue
c
c Index for Dirichlet B.C.
c
        if(jc.eq.1.or.jc.eq.ny) go to 320
          i=i+2
          indx(ip,1)=i-1
          indx(ip,2)=i
          go to 329
  320     continue
          indx(ip,1)=-1
          indx(ip,2)=-1
  329     continue
          if(jcnt.eq.0.or.icnt.eq.0) go to 330
          i=i+1
          insc(ip)=i
          go to 335
  330     insc(ip)=0
  335   continue
  340 continue

      np=ip
      neqn1=i
      iw=2

      if (iw.eq.1) then
        write(iwrite,1089) (i,xc(i),yc(i),indx(i,1),
     &                  indx(i,2),insc(i),i=1,np)
        write(iwrite,1099) (it,(node(it,i),i=1,6),it=1,nelemn)
      end if
c
c  Set the quadrature information
c  quadrature rule is midpoint
c  coordinates and weights are set by routine QD7PT
c
      nquad=3

      do it=1,nelemn
        xm(it,1)=0.5d0
        xm(it,2)=0.5d0
        xm(it,3)=0.0d0
        ym(it,1)=0.0d0
        ym(it,2)=0.5d0
        ym(it,3)=0.5d0
        area(it)=0.5d0
      end do
c
c  half band width
c
      nlband = 0

      do it = 1, nelemn

        do iq = 1, nnodes

          ip = node(it,iq)

          do iuk = 1, nuk

            if ( iuk .eq. 3 ) then
              i = insc(ip)
            else
              i = indx(ip,iuk)
            end if

            if ( 0 .lt. i ) then

              do iqq = 1, nnodes

                ipp = node(it,iqq)

                do iukk = 1, nuk

                  if ( iukk .eq. nuk ) then
                    j = insc(ipp)
                  else
                    j = indx(ipp,iukk)
                  end if

                  if ( i .le. j ) then
                    ij = j - i
                    if ( ij .gt. nlband ) then
                      nlband = ij
                    end if
                  end if

                end do

              end do

            end if

          end do
        end do
      end do

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
c    Hyung-Chun Lee
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
      subroutine trans ( it, xq, yq, det, pj11, pj21, pj12, pj22,
     &  xc, yc, node, maxel )

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
c  Parameters:
c
      implicit none

      integer maxel

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
      integer node(maxel,*)
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
      f1x=x1*(-3.d0+4.d0*xq+4.d0*yq)
     &       +x2*(4.d0*xq-1.d0)
     &       +x4*4.d0*(1.d0-2.d0*xq-yq)
     &       +x5*4.d0*yq + x6*4.d0*(-yq)

      f1y=x1*(-3.d0+4.d0*xq+4.d0*yq)
     &       +x3*(4.d0*yq-1.d0)
     &       +x4*4.d0*(-xq) + x5*4.d0*xq
     &       +x6*4.d0*(1.d0-xq-2.d0*yq)

      f2x=y1*(-3.d0+4.d0*xq+4.d0*yq)
     &       +y2*(4.d0*xq-1.d0)
     &       +y4*4.d0*(1.d0-2.d0*xq-yq)
     &       +y5*4.d0*yq + y6*4.d0*(-yq)

      f2y=y1*(-3.d0+4.d0*xq+4.d0*yq)
     &       +y3*(4.d0*yq-1.d0)
     &       +y4*4.d0*(-xq) + y5*4.d0*xq
     &       +y6*4.d0*(1.d0-xq-2.d0*yq)
c
c  Compute determinant of transformation evaluated at point (xq,yq).
c
      det = f1x * f2y - f1y * f2x
c
c  Compute j11, j22, j21, j22.
c
      pj11 =  f2y / det
      pj12 = -f2x / det
      pj21 = -f1y / det
      pj22 =  f1x / det

      det = dabs ( det )

      return
      end
      function ubdry ( iuk, ip, xc, yc )

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
c    15 April 2005
c
c  Author:
c
c    Hyung-Chun Lee,
c    Department of Mathematics,
c    Ajou University, Korea
c
c  Parameters:
c
c    Input, integer IUK, indicates the type of the unknown.
c    1, horizontal velocity.
c    2, vertical velocity.
c    3, pressure.
c
c    Input, integer IP, the index of the node.
c
c    Input, double precision XC(*), YC(*), the node coordinates.
c
c    Output, double precision UBDRY, the value of the boundary condition
c    applied at this node, if any.
c
      implicit none

      integer ip
      integer iuk
      double precision ubdry
      double precision x
      double precision xc(*)
      double precision y
      double precision yc(*)

      x = xc(ip)
      y = yc(ip)

      if ( x .eq. 0.0 ) then

        if ( iuk .eq. 1 ) then
          ubdry = 16.0D+00 * ( y - 0.5D+00 ) * ( 1.0D+00 - y )
        else if ( iuk .eq. 2 ) then
          ubdry = 0.0D+00
        else
          ubdry = 0.0D+00
        end if

      else

        ubdry = 0.0D+00

      end if

      return
      end
