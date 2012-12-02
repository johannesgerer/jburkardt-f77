      subroutine daxpy ( n, da, dx, incx, dy, incy )

c*********************************************************************72
c
cc DAXPY computes constant times a vector plus a vector.
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
c
c    This routine uses unrolled loops for increments equal to one.
c
c  Modified:
c
c    18 December 2008
c
c  Author:
c
c    Jack Dongarra
c
c  Reference:
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of elements in DX and DY.
c
c    Input, double precision DA, the multiplier of DX.
c
c    Input, double precision DX(*), the first vector.
c
c    Input, integer INCX, the increment between successive entries of DX.
c
c    Input/output, double precision DY(*), the second vector.
c    On output, DY(*) has been replaced by DY(*) + DA * DX(*).
c
c    Input, integer INCY, the increment between successive entries of DY.
c
      implicit none

      double precision da
      double precision dx(*)
      double precision dy(*)
      integer i
      integer incx
      integer incy
      integer ix
      integer iy
      integer m
      integer n

      if ( n .le. 0 ) then
        return
      end if

      if ( da .eq. 0.0d0 ) then
        return
      end if

      if ( incx .ne. 1 .or. incy .ne. 1 ) then

        if ( incx .lt. 0 ) then
          ix = (-n+1)*incx + 1
        else
          ix = 1
        end if

        if ( incy .lt. 0 ) then
          iy = (-n+1)*incy + 1
        else
          iy = 1
        end if

        do i = 1, n
          dy(iy) = dy(iy) + da*dx(ix)
          ix = ix + incx
          iy = iy + incy
        end do

      else

        m = mod(n,4)

        do i = 1, m
          dy(i) = dy(i) + da*dx(i)
        end do

        do i = m + 1, n, 4
          dy(i) = dy(i) + da*dx(i)
          dy(i + 1) = dy(i + 1) + da*dx(i + 1)
          dy(i + 2) = dy(i + 2) + da*dx(i + 2)
          dy(i + 3) = dy(i + 3) + da*dx(i + 3)
        end do

      end if

      return
      end
      function ddot ( n, dx, incx, dy, incy )

c*********************************************************************72
c
cc DDOT forms the dot product of two vectors.
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
c
c    This routine uses unrolled loops for increments equal to one.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Jack Dongarra
c
c  Reference:
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vectors.
c
c    Input, double precision DX(*), the first vector.
c
c    Input, integer INCX, the increment between successive entries in DX.
c
c    Input, double precision DY(*), the second vector.
c
c    Input, integer INCY, the increment between successive entries in DY.
c
c    Output, double precision DDOT, the sum of the product of the 
c    corresponding entries of DX and DY.
c
      implicit none

      double precision ddot
      double precision dx(*)
      double precision dy(*)
      double precision dtemp
      integer i
      integer incx
      integer incy
      integer ix
      integer iy
      integer m
      integer n

      ddot = 0.0d0
      dtemp = 0.0d0

      if ( n .le. 0 ) then
        return
      end if

      if(incx.eq.1.and.incy.eq.1)go to 20
c
c  code for unequal increments or equal increments not equal to 1
c
        if ( incx .lt. 0 ) then
          ix = (-n+1)*incx + 1
        else
          ix = 1
        end if

        if ( incy .lt. 0 ) then
          iy = (-n+1)*incy + 1
        else
          iy = 1
        end if

      do i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
      end do
      ddot = dtemp
      return
c
c  code for both increments equal to 1
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
      end do
      if( n .lt. 5 ) go to 60
   40 continue
      do i = m+1, n, 5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     &   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
      end do

   60 ddot = dtemp

      return
      end
      function dnrm2 ( n, x, incx )

c*********************************************************************72
c
cc DNRM2 returns the euclidean norm of a vector. 
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
c
c      DNRM2 ( X ) = sqrt ( X' * X )
c
c  Author:
c
c    Sven Hammarling
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, double precision X(*), the vector whose norm is to be computed.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Output, double precision DNRM2, the Euclidean norm of X.
c
      implicit none

      double precision dnrm2
      integer incx
      integer ix
      integer n
      double precision x( * )


      double precision one         , zero
      parameter           ( one = 1.0d+0, zero = 0.0d+0 )

      double precision      absxi, norm, scale, ssq

      if( n.lt.1 .or. incx.lt.1 )then
         norm  = zero
      else if( n.eq.1 )then
         norm  = abs( x( 1 ) )
      else
         scale = zero
         ssq   = one
         do ix = 1, 1 + ( n - 1 )*incx, incx
            if( x( ix ).ne.zero )then
               absxi = abs( x( ix ) )
               if( scale.lt.absxi )then
                  ssq   = one   + ssq*( scale/absxi )**2
                  scale = absxi
               else
                  ssq   = ssq   +     ( absxi/scale )**2
               end if
            end if
         end do
         norm  = scale * sqrt( ssq )
      end if

      dnrm2 = norm

      return
      end
      subroutine drot ( n, dx, incx, dy, incy, c, s )

c*********************************************************************72
c
cc DROT applies a plane rotation.
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Jack Dongarra
c
c  Reference:
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vectors.
c
c    Input/output, double precision X(*), one of the vectors to be rotated.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Input/output, double precision Y(*), one of the vectors to be rotated.
c
c    Input, integer INCY, the increment between successive elements of Y.
c
c    Input, double precision C, S, parameters (presumably the cosine and
c    sine of some angle) that define a plane rotation.
c
      implicit none

      double precision dtemp
      double precision dx(*)
      double precision dy(*),c,s
      integer i,incx,incy,ix,iy,n

      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c  code for unequal increments or equal increments not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        dtemp = c*dx(ix) + s*dy(iy)
        dy(iy) = c*dy(iy) - s*dx(ix)
        dx(ix) = dtemp
        ix = ix + incx
        iy = iy + incy
      end do
      return
c
c  code for both increments equal to 1
c
   20 do i = 1,n
        dtemp = c*dx(i) + s*dy(i)
        dy(i) = c*dy(i) - s*dx(i)
        dx(i) = dtemp
      end do
      return
      end
      subroutine drotg ( da, db, c, s )

c*********************************************************************72
c
cc DROTG constructs a Givens plane rotation.
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
c
c    Given values A and B, this routine computes
c
c    SIGMA = sign ( A ) if abs ( A ) >  abs ( B )
c          = sign ( B ) if abs ( A ) <= abs ( B );
c
c    R     = SIGMA * ( A * A + B * B );
c
c    C = A / R if R is not 0
c      = 1     if R is 0;
c
c    S = B / R if R is not 0,
c        0     if R is 0.
c
c    The computed numbers then satisfy the equation
c
c    (  C  S ) ( A ) = ( R )
c    ( -S  C ) ( B ) = ( 0 )
c
c    The routine also computes
c
c    Z = S     if abs ( A ) > abs ( B ),
c      = 1 / C if abs ( A ) <= abs ( B ) and C is not 0,
c      = 1     if C is 0.
c
c    The single value Z encodes C and S, and hence the rotation:
c
c    If Z = 1, set C = 0 and S = 1;
c    If abs ( Z ) < 1, set C = sqrt ( 1 - Z * Z ) and S = Z;
c    if abs ( Z ) > 1, set C = 1/ Z and S = sqrt ( 1 - C * C );
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Jack Dongarra
c
c  Reference:
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input/output, double precision SA, SB.  On input, SA and SB are the values
c    A and B.  On output, SA is overwritten with R, and SB is
c    overwritten with Z.
c
c    Output, double precision C, S, the cosine and sine of the
c    Givens rotation.
c
      implicit none

      double precision c
      double precision da
      double precision db,s,roe,scale,r,z

      if( dabs(da) .gt. dabs(db) ) then
        roe = da
      else
        roe = db
      end if

      scale = dabs(da) + dabs(db)

      if( scale .eq. 0.0d0 ) then
         c = 1.0d0
         s = 0.0d0
         r = 0.0d0
         z = 0.0d0
      else
        r = scale*dsqrt((da/scale)**2 + (db/scale)**2)
        r = dsign(1.0d0,roe)*r
        c = da/r
        s = db/r
        z = 1.0d0
        if( dabs(da) .gt. dabs(db) ) z = s
        if( dabs(db) .ge. dabs(da) .and. c .ne. 0.0d0 ) z = 1.0d0/c
      end if

      da = r
      db = z

      return
      end
      subroutine dscal ( n, da, dx, incx )

c*********************************************************************72
c
cc DSCAL scales a vector by a constant.
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Jack Dongarra
c
c  Reference:
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, double precision SA, the multiplier.
c
c    Input/output, double precision X(*), the vector to be scaled.
c
c    Input, integer INCX, the increment between successive entries of X.
c
      implicit none

      double precision da
      double precision dx(*)
      integer i
      integer incx,m,n,nincx

      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c  code for increment not equal to 1
c
      nincx = n*incx
      do i = 1,nincx,incx
        dx(i) = da*dx(i)
      end do
      return
c
c  code for increment equal to 1
c
c
c  clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        dx(i) = da*dx(i)
      end do
      if( n .lt. 5 ) return
   40 continue
      do i = m+1, n, 5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
      end do

      return
      end
      subroutine dsvdc ( x, ldx, n, p, s, e, u, ldu, v, ldv, work, 
     &  job, info )

c*********************************************************************72
c
cc DSVDC computes the singular value decomposition of a matrix.
c
c     dsvdc is a subroutine to reduce a double precision nxp matrix x
c     by orthogonal transformations u and v to diagonal form.  the
c     diagonal elements s(i) are the singular values of x.  the
c     columns of u are the corresponding left singular vectors,
c     and the columns of v the right singular vectors.
c
c     on entry
c
c         x         double precision(ldx,p), where ldx.ge.n.
c                   x contains the matrix whose singular value
c                   decomposition is to be computed.  x is
c                   destroyed by dsvdc.
c
c         ldx       integer.
c                   ldx is the leading dimension of the array x.
c
c         n         integer.
c                   n is the number of rows of the matrix x.
c
c         p         integer.
c                   p is the number of columns of the matrix x.
c
c         ldu       integer.
c                   ldu is the leading dimension of the array u.
c                   (see below).
c
c         ldv       integer.
c                   ldv is the leading dimension of the array v.
c                   (see below).
c
c         work      double precision(n).
c                   work is a scratch array.
c
c         job       integer.
c                   job controls the computation of the singular
c                   vectors.  it has the decimal expansion ab
c                   with the following meaning
c
c                        a.eq.0    do not compute the left singular
c                                  vectors.
c                        a.eq.1    return the n left singular vectors
c                                  in u.
c                        a.ge.2    return the first min(n,p) singular
c                                  vectors in u.
c                        b.eq.0    do not compute the right singular
c                                  vectors.
c                        b.eq.1    return the right singular vectors
c                                  in v.
c
c     on return
c
c         s         double precision(mm), where mm=max(n+1,p).
c                   the first min(n,p) entries of s contain the
c                   singular values of x arranged in descending
c                   order of magnitude.
c
c         e         double precision(mm), where mm=max(n+1,p).
c                   e ordinarily contains zeros.  however see the
c                   discussion of info for exceptions.
c
c         u         double precision(ldu,k), where ldu.ge.n.  if
c                                   joba.eq.1 then k.eq.n, if joba.ge.2
c                                   then k.eq.min(n,p).
c                   u contains the matrix of left singular vectors.
c                   u is not referenced if joba.eq.0.  if n.le.p
c                   or if joba.eq.2, then u may be identified with x
c                   in the subroutine call.
c
c         v         double precision(ldv,p), where ldv.ge.p.
c                   v contains the matrix of right singular vectors.
c                   v is not referenced if job.eq.0.  if p.le.n,
c                   then v may be identified with x in the
c                   subroutine call.
c
c         info      integer.
c                   the singular values (and their corresponding
c                   singular vectors) s(info+1),s(info+2),...,s(m)
c                   are correct (here m=min(n,p)).  thus if
c                   info.eq.0, all the singular values and their
c                   vectors are correct.  in any event, the matrix
c                   b = trans(u)*x*v is the bidiagonal matrix
c                   with the elements of s on its diagonal and the
c                   elements of e on its super-diagonal (trans(u)
c                   is the transpose of u).  thus the singular
c                   values of x and b are the same.
c
c     linpack. this version dated 08/14/78 .
c              correction made to shift 2/84.
c     g.w. stewart, university of maryland, argonne national lab.
c
c     dsvdc uses the following functions and subprograms.
c
c     external drot
c     blas daxpy,ddot,dscal,dswap,dnrm2,drotg
c     fortran dabs,dmax1,max0,min0,mod,dsqrt
c
      integer ldx,n,p,ldu,ldv,job,info
      double precision x(ldx,*),s(*),e(*),u(ldu,*),v(ldv,*),work(*)

      integer i,iter,j,jobu,k,kase,kk,l,ll,lls,lm1,lp1,ls,lu,m,maxit,
     *        mm,mm1,mp1,nct,nctp1,ncu,nrt,nrtp1
      double precision ddot,t,r
      double precision b,c,cs,el,emm1,f,g,dnrm2,scale,shift,sl,sm,sn,
     *                 smm1,t1,test,ztest
      logical wantu,wantv
c
c
c     set the maximum number of iterations.
c
      maxit = 30
c
c     determine what is to be computed.
c
      wantu = .false.
      wantv = .false.
      jobu = mod(job,100)/10
      ncu = n
      if (jobu .gt. 1) ncu = min0(n,p)
      if (jobu .ne. 0) wantu = .true.
      if (mod(job,10) .ne. 0) wantv = .true.
c
c     reduce x to bidiagonal form, storing the diagonal elements
c     in s and the super-diagonal elements in e.
c
      info = 0
      nct = min0(n-1,p)
      nrt = max0(0,min0(p-2,n))
      lu = max0(nct,nrt)
      if (lu .lt. 1) go to 170
      do 160 l = 1, lu
         lp1 = l + 1
         if (l .gt. nct) go to 20
c
c           compute the transformation for the l-th column and
c           place the l-th diagonal in s(l).
c
            s(l) = dnrm2(n-l+1,x(l,l),1)
            if (s(l) .eq. 0.0d0) go to 10
               if (x(l,l) .ne. 0.0d0) s(l) = dsign(s(l),x(l,l))
               call dscal(n-l+1,1.0d0/s(l),x(l,l),1)
               x(l,l) = 1.0d0 + x(l,l)
   10       continue
            s(l) = -s(l)
   20    continue
         if (p .lt. lp1) go to 50
         do 40 j = lp1, p
            if (l .gt. nct) go to 30
            if (s(l) .eq. 0.0d0) go to 30
c
c              apply the transformation.
c
               t = -ddot(n-l+1,x(l,l),1,x(l,j),1)/x(l,l)
               call daxpy(n-l+1,t,x(l,l),1,x(l,j),1)
   30       continue
c
c           place the l-th row of x into  e for the
c           subsequent calculation of the row transformation.
c
            e(j) = x(l,j)
   40    continue
   50    continue
         if (.not.wantu .or. l .gt. nct) go to 70
c
c           place the transformation in u for subsequent back
c           multiplication.
c
            do 60 i = l, n
               u(i,l) = x(i,l)
   60       continue
   70    continue
         if (l .gt. nrt) go to 150
c
c           compute the l-th row transformation and place the
c           l-th super-diagonal in e(l).
c
            e(l) = dnrm2(p-l,e(lp1),1)
            if (e(l) .eq. 0.0d0) go to 80
               if (e(lp1) .ne. 0.0d0) e(l) = dsign(e(l),e(lp1))
               call dscal(p-l,1.0d0/e(l),e(lp1),1)
               e(lp1) = 1.0d0 + e(lp1)
   80       continue
            e(l) = -e(l)
            if (lp1 .gt. n .or. e(l) .eq. 0.0d0) go to 120
c
c              apply the transformation.
c
               do 90 i = lp1, n
                  work(i) = 0.0d0
   90          continue
               do 100 j = lp1, p
                  call daxpy(n-l,e(j),x(lp1,j),1,work(lp1),1)
  100          continue
               do 110 j = lp1, p
                  call daxpy(n-l,-e(j)/e(lp1),work(lp1),1,x(lp1,j),1)
  110          continue
  120       continue
            if (.not.wantv) go to 140
c
c              place the transformation in v for subsequent
c              back multiplication.
c
               do 130 i = lp1, p
                  v(i,l) = e(i)
  130          continue
  140       continue
  150    continue
  160 continue
  170 continue
c
c     set up the final bidiagonal matrix or order m.
c
      m = min0(p,n+1)
      nctp1 = nct + 1
      nrtp1 = nrt + 1
      if (nct .lt. p) s(nctp1) = x(nctp1,nctp1)
      if (n .lt. m) s(m) = 0.0d0
      if (nrtp1 .lt. m) e(nrtp1) = x(nrtp1,m)
      e(m) = 0.0d0
c
c     if required, generate u.
c
      if (.not.wantu) go to 300
         if (ncu .lt. nctp1) go to 200
         do 190 j = nctp1, ncu
            do 180 i = 1, n
               u(i,j) = 0.0d0
  180       continue
            u(j,j) = 1.0d0
  190    continue
  200    continue
         if (nct .lt. 1) go to 290
         do 280 ll = 1, nct
            l = nct - ll + 1
            if (s(l) .eq. 0.0d0) go to 250
               lp1 = l + 1
               if (ncu .lt. lp1) go to 220
               do 210 j = lp1, ncu
                  t = -ddot(n-l+1,u(l,l),1,u(l,j),1)/u(l,l)
                  call daxpy(n-l+1,t,u(l,l),1,u(l,j),1)
  210          continue
  220          continue
               call dscal(n-l+1,-1.0d0,u(l,l),1)
               u(l,l) = 1.0d0 + u(l,l)
               lm1 = l - 1
               if (lm1 .lt. 1) go to 240
               do 230 i = 1, lm1
                  u(i,l) = 0.0d0
  230          continue
  240          continue
            go to 270
  250       continue
               do 260 i = 1, n
                  u(i,l) = 0.0d0
  260          continue
               u(l,l) = 1.0d0
  270       continue
  280    continue
  290    continue
  300 continue
c
c     if it is required, generate v.
c
      if (.not.wantv) go to 350
         do 340 ll = 1, p
            l = p - ll + 1
            lp1 = l + 1
            if (l .gt. nrt) go to 320
            if (e(l) .eq. 0.0d0) go to 320
               do 310 j = lp1, p
                  t = -ddot(p-l,v(lp1,l),1,v(lp1,j),1)/v(lp1,l)
                  call daxpy(p-l,t,v(lp1,l),1,v(lp1,j),1)
  310          continue
  320       continue
            do 330 i = 1, p
               v(i,l) = 0.0d0
  330       continue
            v(l,l) = 1.0d0
  340    continue
  350 continue
c
c     main iteration loop for the singular values.
c
      mm = m
      iter = 0
  360 continue
c
c        quit if all the singular values have been found.
c
c     ...exit
         if (m .eq. 0) go to 620
c
c        if too many iterations have been performed, set
c        flag and return.
c
         if (iter .lt. maxit) go to 370
            info = m
c     ......exit
            go to 620
  370    continue
c
c        this section of the program inspects for
c        negligible elements in the s and e arrays.  on
c        completion the variables kase and l are set as follows.
c
c           kase = 1     if s(m) and e(l-1) are negligible and l.lt.m
c           kase = 2     if s(l) is negligible and l.lt.m
c           kase = 3     if e(l-1) is negligible, l.lt.m, and
c                        s(l), ..., s(m) are not negligible (qr step).
c           kase = 4     if e(m-1) is negligible (convergence).
c
         do 390 ll = 1, m
            l = m - ll
c        ...exit
            if (l .eq. 0) go to 400
            test = dabs(s(l)) + dabs(s(l+1))
            ztest = test + dabs(e(l))
            if (ztest .ne. test) go to 380
               e(l) = 0.0d0
c        ......exit
               go to 400
  380       continue
  390    continue
  400    continue
         if (l .ne. m - 1) go to 410
            kase = 4
         go to 480
  410    continue
            lp1 = l + 1
            mp1 = m + 1
            do 430 lls = lp1, mp1
               ls = m - lls + lp1
c           ...exit
               if (ls .eq. l) go to 440
               test = 0.0d0
               if (ls .ne. m) test = test + dabs(e(ls))
               if (ls .ne. l + 1) test = test + dabs(e(ls-1))
               ztest = test + dabs(s(ls))
               if (ztest .ne. test) go to 420
                  s(ls) = 0.0d0
c           ......exit
                  go to 440
  420          continue
  430       continue
  440       continue
            if (ls .ne. l) go to 450
               kase = 3
            go to 470
  450       continue
            if (ls .ne. m) go to 460
               kase = 1
            go to 470
  460       continue
               kase = 2
               l = ls
  470       continue
  480    continue
         l = l + 1
c
c        perform the task indicated by kase.
c
         go to (490,520,540,570), kase
c
c        deflate negligible s(m).
c
  490    continue
            mm1 = m - 1
            f = e(m-1)
            e(m-1) = 0.0d0
            do 510 kk = l, mm1
               k = mm1 - kk + l
               t1 = s(k)
               call drotg(t1,f,cs,sn)
               s(k) = t1
               if (k .eq. l) go to 500
                  f = -sn*e(k-1)
                  e(k-1) = cs*e(k-1)
  500          continue
               if (wantv) call drot(p,v(1,k),1,v(1,m),1,cs,sn)
  510       continue
         go to 610
c
c        split at negligible s(l).
c
  520    continue
            f = e(l-1)
            e(l-1) = 0.0d0
            do 530 k = l, m
               t1 = s(k)
               call drotg(t1,f,cs,sn)
               s(k) = t1
               f = -sn*e(k)
               e(k) = cs*e(k)
               if (wantu) call drot(n,u(1,k),1,u(1,l-1),1,cs,sn)
  530       continue
         go to 610
c
c        perform one qr step.
c
  540    continue
c
c           calculate the shift.
c
            scale = dmax1(dabs(s(m)),dabs(s(m-1)),dabs(e(m-1)),
     *                    dabs(s(l)),dabs(e(l)))
            sm = s(m)/scale
            smm1 = s(m-1)/scale
            emm1 = e(m-1)/scale
            sl = s(l)/scale
            el = e(l)/scale
            b = ((smm1 + sm)*(smm1 - sm) + emm1**2)/2.0d0
            c = (sm*emm1)**2
            shift = 0.0d0
            if (b .eq. 0.0d0 .and. c .eq. 0.0d0) go to 550
               shift = dsqrt(b**2+c)
               if (b .lt. 0.0d0) shift = -shift
               shift = c/(b + shift)
  550       continue
            f = (sl + sm)*(sl - sm) + shift
            g = sl*el
c
c           chase zeros.
c
            mm1 = m - 1
            do 560 k = l, mm1
               call drotg(f,g,cs,sn)
               if (k .ne. l) e(k-1) = f
               f = cs*s(k) + sn*e(k)
               e(k) = cs*e(k) - sn*s(k)
               g = sn*s(k+1)
               s(k+1) = cs*s(k+1)
               if (wantv) call drot(p,v(1,k),1,v(1,k+1),1,cs,sn)
               call drotg(f,g,cs,sn)
               s(k) = f
               f = cs*e(k) + sn*s(k+1)
               s(k+1) = -sn*e(k) + cs*s(k+1)
               g = sn*e(k+1)
               e(k+1) = cs*e(k+1)
               if (wantu .and. k .lt. n)
     *            call drot(n,u(1,k),1,u(1,k+1),1,cs,sn)
  560       continue
            e(m-1) = f
            iter = iter + 1
         go to 610
c
c        convergence.
c
  570    continue
c
c           make the singular value  positive.
c
            if (s(l) .ge. 0.0d0) go to 580
               s(l) = -s(l)
               if (wantv) call dscal(p,-1.0d0,v(1,l),1)
  580       continue
c
c           order the singular value.
c
  590       if (l .eq. mm) go to 600
c           ...exit
               if (s(l) .ge. s(l+1)) go to 600
               t = s(l)
               s(l) = s(l+1)
               s(l+1) = t
               if (wantv .and. l .lt. p)
     *            call dswap(p,v(1,l),1,v(1,l+1),1)
               if (wantu .and. l .lt. n)
     *            call dswap(n,u(1,l),1,u(1,l+1),1)
               l = l + 1
            go to 590
  600       continue
            iter = 0
            m = m - 1
  610    continue
      go to 360
  620 continue
      return
      end
      subroutine dswap ( n, dx, incx, dy, incy )

c*********************************************************************72
c
cc DSWAP interchanges two vectors.
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Jack Dongarra
c
c  Reference:
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vectors.
c
c    Input/output, double precision X(*), one of the vectors to swap.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Input/output, double precision Y(*), one of the vectors to swap.
c
c    Input, integer INCY, the increment between successive elements of Y.
c
      implicit none

      double precision dx(*)
      double precision dy(*),dtemp
      integer i,incx,incy,ix,iy,m,n

      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c  code for unequal increments or equal increments not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
      end do
      return
c
c  code for both increments equal to 1
c
c
c  clean-up loop
c
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
      end do
      if( n .lt. 3 ) return
   40 continue

      do i = m+1, n, 3
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
        dtemp = dx(i + 1)
        dx(i + 1) = dy(i + 1)
        dy(i + 1) = dtemp
        dtemp = dx(i + 2)
        dx(i + 2) = dy(i + 2)
        dy(i + 2) = dtemp
      end do

      return
      end
      subroutine phi1 ( n, r, r0, v )

c*********************************************************************72
c
cc PHI1 evaluates the multiquadric radial basis function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
c    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
c    Third Edition,
c    Cambridge University Press, 2007,
c    ISBN13: 978-0-521-88068-8,
c    LC: QA297.N866.
c
c  Parameters:
c
c    Input, integer N, the number of points.
c
c    Input, double precision R(N), the radial separation.
c    0 < R.
c
c    Input, double precision R0, a scale factor.
c
c    Output, double precision V(N), the value of the radial basis function.
c
      implicit none

      integer n

      integer i
      double precision r(n)
      double precision r0
      double precision v(n)

      do i = 1, n
        v(i) = sqrt ( r(i)**2 + r0**2 )
      end do

      return
      end
      subroutine phi2 ( n, r, r0, v )

c*********************************************************************72
c
cc PHI2 evaluates the inverse multiquadric radial basis function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
c    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
c    Third Edition,
c    Cambridge University Press, 2007,
c    ISBN13: 978-0-521-88068-8,
c    LC: QA297.N866.
c
c  Parameters:
c
c    Input, integer N, the number of points.
c
c    Input, double precision R(N), the radial separation.
c    0 < R.
c
c    Input, double precision R0, a scale factor.
c
c    Output, double precision V(N), the value of the radial basis function.
c
      implicit none

      integer n

      integer i
      double precision r(n)
      double precision r0
      double precision v(n)

      do i = 1, n
        v(i) = 1.0D+00 / sqrt ( r(i)**2 + r0**2 )
      end do

      return
      end
      subroutine phi3 ( n, r, r0, v )

c*********************************************************************72
c
cc PHI3 evaluates the thin-plate spline radial basis function.
c
c  Discussion:
c
c    Note that PHI3(R,R0) is negative if R < R0.  Thus, for this basis function,
c    it may be desirable to choose a value of R0 smaller than any possible R.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
c    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
c    Third Edition,
c    Cambridge University Press, 2007,
c    ISBN13: 978-0-521-88068-8,
c    LC: QA297.N866.
c
c  Parameters:
c
c    Input, integer N, the number of points.
c
c    Input, double precision R(N), the radial separation.
c    0 < R.
c
c    Input, double precision R0, a scale factor.
c
c    Output, double precision V(N), the value of the radial basis function.
c
      implicit none

      integer n

      integer i
      double precision r(n)
      double precision r0
      double precision v(n)

      do i = 1, n
        if ( r(i) .le. 0.0D+00 ) then
          v(i) = 0.0D+00
        else
          v(i) = r(i)**2 * log ( r(i) / r0 )
        end if
      end do

      return
      end
      subroutine phi4 ( n, r, r0, v )

c*********************************************************************72
c
cc PHI4 evaluates the gaussian radial basis function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
c    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
c    Third Edition,
c    Cambridge University Press, 2007,
c    ISBN13: 978-0-521-88068-8,
c    LC: QA297.N866.
c
c  Parameters:
c
c    Input, integer N, the number of points.
c
c    Input, double precision R(N), the radial separation.
c    0 < R.
c
c    Input, double precision R0, a scale factor.
c
c    Output, double precision V(N), the value of the radial basis function.
c
      implicit none

      integer n

      integer i
      double precision r(n)
      double precision r0
      double precision v(n)

      do i = 1, n
        v(i) = exp ( - 0.5D+00 * r(i)**2 / r0**2 )
      end do

      return
      end
      subroutine r8mat_solve_svd ( m, n, a, b, x )

c*********************************************************************72
c
cc R8MAT_SOLVE_SVD solves a linear system A*x=b using the SVD.
c
c  Discussion:
c
c    When the system is determined, the solution is the solution in the
c    ordinary sense, and A*x = b.
c
c    When the system is overdetermined, the solution minimizes the
c    L2 norm of the residual ||A*x-b||.
c
c    When the system is underdetermined, ||A*x-b|| should be zero, and
c    the solution is the solution of minimum L2 norm, ||x||.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns
c    in the matrix A.
c
c    Input, double precision A(M,N), the matrix.
c
c    Input, double precision B(M), the right hand side.
c
c    Output, double precision X(N), the solution.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision a_copy(m,n)
      double precision a_pseudo(n,m)
      double precision b(m)
      double precision e(max(m+1,n))
      integer i
      integer info
      integer j
      integer k
      integer l
      integer lda
      integer ldu
      integer ldv
      integer job
      integer lwork
      double precision s(m,n)
      double precision sp(n,m)
      double precision sdiag(max(m+1,n))
      double precision u(m,m)
      double precision v(n,n)
      double precision work(m)
      double precision x(n)
c
c  Compute the SVD decomposition.
c
      job = 11
      lda = m
      ldu = m
      ldv = n

      do j = 1, n
        do i = 1, m
          a_copy(i,j) = a(i,j)
        end do
      end do

      call dsvdc ( a_copy, lda, m, n, sdiag, e, u, ldu, v, ldv, work, 
     &  job, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8MAT_SOLVE_SVD - Failure!'
        write ( *, '(a)' ) '  The SVD could not be calculated.'
        write ( *, '(a)' ) '  LINPACK routine DSVDC returned a nonzero'
        write ( *, '(a,i8)' ) '  value of the error flag, INFO = ', info
        return
      end if

      do j = 1, n
        do i = 1, m
          s(i,j) = 0.0D+00
        end do
      end do

      do i = 1, min ( m, n )
        s(i,i) = sdiag(i)
      end do
c
c  Compute the pseudo inverse.
c
      do j = 1, m
        do i = 1, n
          sp(i,j) = 0.0D+00
        end do
      end do

      do i = 1, min ( m, n )
        if ( s(i,i) .ne. 0.0D+00 ) then
          sp(i,i) = 1.0D+00 / s(i,i)
        end if
      end do

      do j = 1, m
        do i = 1, n
          a_pseudo(i,j) = 0.0D+00
          do k = 1, n
            do l = 1, m
              a_pseudo(i,j) = a_pseudo(i,j) + v(i,k) * sp(k,l) * u(j,l)
            end do
          end do
        end do
      end do
c
c  Compute x = A_pseudo * b
c
      do i = 1, n
        x(i) = 0.0D+00
        do j = 1, m
          x(i) = x(i) + a_pseudo(i,j) * b(j)
        end do
      end do

      return
      end
      subroutine rbf_interp ( m, nd, xd, r0, phi, w, ni, xi, fi )

c*********************************************************************72
c
cc RBF_INTERP evaluates a radial basis function interpolant.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
c    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
c    Third Edition,
c    Cambridge University Press, 2007,
c    ISBN13: 978-0-521-88068-8,
c    LC: QA297.N866.
c
c  Parameters:
c
c    Input, integer M, the spatial dimension.
c
c    Input, integer ND, the number of data points.
c
c    Input, double precision XD(M,ND), the data points.
c
c    Input, double precision R0, a scale factor.  R0 should be larger than the typical
c    separation between points, but smaller than the maximum separation.
c    The value of R0 has a significant effect on the resulting interpolant.
c
c    Input, subroutine PHI ( N, R, R0, V ), a subroutine to evaluate the radial
c    basis functions.
c
c    Input, double precision W(ND), the weights, as computed by RBF_WEIGHTS.
c
c    Input, integer NI, the number of interpolation points.
c
c    Input, double precision XI(M,NI), the interpolation points.
c
c    Output, double precision FI(NI), the interpolated values.
c
      implicit none

      integer m
      integer nd
      integer ni

      double precision fi(ni)
      integer i
      integer j
      external phi
      double precision r(nd)
      double precision r0
      double precision v(nd)
      double precision w(nd)
      double precision xd(m,nd)
      double precision xi(m,ni)

      do i = 1, ni

        do j = 1, nd
          r(j) = sqrt ( sum ( ( xi(1:m,i) - xd(1:m,j) )**2 ) )
        end do

        call phi ( nd, r, r0, v )

        fi(i) = dot_product ( v, w )

      end do

      return
      end
      subroutine rbf_weight ( m, nd, xd, r0, phi, fd, w )

c*********************************************************************72
c
cc RBF_WEIGHT computes weights for radial basis function interpolation.
c
c  Discussion:
c
c    We assume that there are N (nonsingular) equations in N unknowns.
c
c    However, it should be clear that, if we are willing to do some kind
c    of least squares calculation, we could allow for singularity,
c    inconsistency, or underdetermine systems.  This could be associated
c    with data points that are very close or repeated, a smaller number
c    of data points than function values, or some other ill-conditioning
c    of the system arising from a peculiarity in the point spacing.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
c    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
c    Third Edition,
c    Cambridge University Press, 2007,
c    ISBN13: 978-0-521-88068-8,
c    LC: QA297.N866.
c
c  Parameters:
c
c    Input, integer M, the spatial dimension.
c
c    Input, integer ND, the number of data points.
c
c    Input, double precision XD(M,ND), the data points.
c
c    Input, double precision R0, a scale factor.  R0 should be larger than the typical
c    separation between points, but smaller than the maximum separation.
c    The value of R0 has a significant effect on the resulting interpolant.
c
c    Input, subroutine PHI ( N, R, R0, V ), a subroutine to evaluate the radial
c    basis functions.
c
c    Input, double precision FD(ND), the function values at the data points.
c
c    Output, double precision W(ND), the weights.
c
      implicit none

      integer m
      integer nd

      double precision a(nd,nd)
      double precision fd(nd)
      integer i
      integer j
      external phi
      double precision r(nd)
      double precision r0
      double precision v(nd)
      double precision w(nd)
      double precision xd(m,nd)

      do i = 1, nd

        do j = 1, nd
          r(j) = sqrt ( sum ( ( xd(1:m,i) - xd(1:m,j) )**2 ) )
        end do

        call phi ( nd, r, r0, v )

        do j = 1, nd
          a(i,j) = v(j)
        end do

      end do
c
c  Solve for the weights.
c
      call r8mat_solve_svd ( nd, nd, a, fd, w )

      return
      end
