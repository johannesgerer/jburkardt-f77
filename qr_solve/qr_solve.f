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

      if ( da .eq. 0.0D+00 ) then
        return
      end if

      if ( incx .ne. 1 .or. incy .ne. 1 ) then

        if ( incx .lt. 0 ) then
          ix = (-n+1) * incx + 1
        else
          ix = 1
        end if

        if ( incy .lt. 0 ) then
          iy = (-n+1) * incy + 1
        else
          iy = 1
        end if

        do i = 1, n
          dy(iy) = dy(iy) + da * dx(ix)
          ix = ix + incx
          iy = iy + incy
        end do

      else

        m = mod ( n, 4 )

        do i = 1, m
          dy(i) = dy(i) + da * dx(i)
        end do

        do i = m + 1, n, 4
          dy(i)     = dy(i)     + da * dx(i)
          dy(i + 1) = dy(i + 1) + da * dx(i + 1)
          dy(i + 2) = dy(i + 2) + da * dx(i + 2)
          dy(i + 3) = dy(i + 3) + da * dx(i + 3)
        end do

      end if

      return
      end
      subroutine dcopy ( n, dx, incx, dy, incy )

c*********************************************************************72
c
cc DCOPY copies a vector.
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
c
c    The routine uses unrolled loops for increments equal to one.
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
c    Input, integer N, the number of elements in DX and DY.
c
c    Input, double precision DX(*), the first vector.
c
c    Input, integer INCX, the increment between successive entries of DX.
c
c    Output, double precision DY(*), the second vector.
c
c    Input, integer INCY, the increment between successive entries of DY.
c
      implicit none

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
c
c  code for unequal increments or equal increments not equal to 1
c
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

        do i = 1,n
          dy(iy) = dx(ix)
          ix = ix + incx
          iy = iy + incy
        end do
c
c  code for both increments equal to 1
c
      else

        m = mod(n,7)

        do i = 1,m
          dy(i) = dx(i)
        end do

        do i = m + 1, n, 7
          dy(i) = dx(i)
          dy(i + 1) = dx(i + 1)
          dy(i + 2) = dx(i + 2)
          dy(i + 3) = dx(i + 3)
          dy(i + 4) = dx(i + 4)
          dy(i + 5) = dx(i + 5)
          dy(i + 6) = dx(i + 6)
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

      ddot = 0.0D+00
      dtemp = 0.0D+00

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
          iy = (-n+1) * incy + 1
        else
          iy = 1
        end if

      do i = 1, n
        dtemp = dtemp + dx(ix) * dy(iy)
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
        dtemp = dtemp + dx(i) * dy(i)
      end do
      if( n .lt. 5 ) go to 60
   40 continue
      do i = m+1, n, 5
        dtemp = dtemp + dx(i)     * dy(i) 
     &                + dx(i + 1) * dy(i + 1) 
     &                + dx(i + 2) * dy(i + 2) 
     &                + dx(i + 3) * dy(i + 3) 
     &                + dx(i + 4) * dy(i + 4)
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

      double precision absxi
      double precision dnrm2
      integer incx
      integer ix
      integer n
      double precision norm
      double precision one
      double precision scale
      double precision ssq
      double precision x( * )
      double precision zero


      parameter ( one = 1.0d+0, zero = 0.0d+0 )

      if( n.lt.1 .or. incx.lt.1 )then
        norm  = zero
      else if( n.eq.1 )then
        norm  = abs( x( 1 ) )
      else
        scale = zero
        ssq   = one
        do ix = 1, 1 + ( n - 1 ) * incx, incx
          if( x( ix ).ne.zero )then
            absxi = abs( x( ix ) )
            if( scale.lt.absxi )then
              ssq   = one   + ssq * ( scale/absxi )**2
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
      subroutine dqrank ( a, lda, m, n, tol, kr, jpvt, qraux, work )

c*********************************************************************72
c
cc DQRANK computes the QR factorization of a rectangular matrix.
c
c  Discussion:
c
c    This routine is used in conjunction with sqrlss to solve
c    overdetermined, underdetermined and singular linear systems
c    in a least squares sense.
c
c    DQRANK uses the LINPACK subroutine DQRDC to compute the QR
c    factorization, with column pivoting, of an M by N matrix A.
c    The numerical rank is determined using the tolerance TOL.
c
c    Note that on output, ABS ( A(1,1) ) / ABS ( A(KR,KR) ) is an estimate
c    of the condition number of the matrix of independent columns,
c    and of R.  This estimate will be <= 1/TOL.
c
c  Reference:
c
c    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c  Parameters:
c
c    Input/output, double precision A(LDA,N).  On input, the matrix whose
c    decomposition is to be computed.  On output, the information from DQRDC.
c    The triangular matrix R of the QR factorization is contained in the
c    upper triangle and information needed to recover the orthogonal
c    matrix Q is stored below the diagonal in A and in the vector QRAUX.
c
c    Input, integer LDA, the leading dimension of A, which must
c    be at least M.
c
c    Input, integer M, the number of rows of A.
c
c    Input, integer N, the number of columns of A.
c
c    Input, double precision TOL, a relative tolerance used to determine the
c    numerical rank.  The problem should be scaled so that all the elements
c    of A have roughly the same absolute accuracy, EPS.  Then a reasonable
c    value for TOL is roughly EPS divided by the magnitude of the largest
c    element.
c
c    Output, integer KR, the numerical rank.
c
c    Output, integer JPVT(N), the pivot information from DQRDC.
c    Columns JPVT(1), ..., JPVT(KR) of the original matrix are linearly
c    independent to within the tolerance TOL and the remaining columns
c    are linearly dependent.
c
c    Output, double precision QRAUX(N), will contain extra information defining
c    the QR factorization.
c
c    Workspace, double precision WORK(N).
c
      implicit none

      integer lda
      integer n

      double precision a(lda,n)
      integer j
      integer jpvt(n)
      integer k
      integer kr
      integer m
      double precision qraux(n)
      double precision tol
      double precision work(n)

      do j = 1, n
        jpvt(j) = 0
      end do

      call dqrdc ( a, lda, m, n, qraux, jpvt, work, 1 )

      kr = 0
      k = min ( m, n )

      do j = 1, k
        if ( abs ( a(j,j) ) .le. tol * abs ( a(1,1) ) ) then
          return
        end if
        kr = j
      end do

      return
      end
      subroutine dqrdc(x,ldx,n,p,qraux,jpvt,work,job)

c*********************************************************************72

      integer ldx,n,p,job
      integer jpvt(1)
      double precision x(ldx,1),qraux(1),work(1)
c
c     dqrdc uses householder transformations to compute the qr
c     factorization of an n by p matrix x.  column pivoting
c     based on the 2-norms of the reduced columns may be
c     performed at the users option.
c
c     on entry
c
c        x       double precision(ldx,p), where ldx .ge. n.
c                x contains the matrix whose decomposition is to be
c                computed.
c
c        ldx     integer.
c                ldx is the leading dimension of the array x.
c
c        n       integer.
c                n is the number of rows of the matrix x.
c
c        p       integer.
c                p is the number of columns of the matrix x.
c
c        jpvt    integer(p).
c                jpvt contains integers that control the selection
c                of the pivot columns.  the k-th column x(k) of x
c                is placed in one of three classes according to the
c                value of jpvt(k).
c
c                   if jpvt(k) .gt. 0, then x(k) is an initial
c                                      column.
c
c                   if jpvt(k) .eq. 0, then x(k) is a free column.
c
c                   if jpvt(k) .lt. 0, then x(k) is a final column.
c
c                before the decomposition is computed, initial columns
c                are moved to the beginning of the array x and final
c                columns to the end.  both initial and final columns
c                are frozen in place during the computation and only
c                free columns are moved.  at the k-th stage of the
c                reduction, if x(k) is occupied by a free column
c                it is interchanged with the free column of largest
c                reduced norm.  jpvt is not referenced if
c                job .eq. 0.
c
c        work    double precision(p).
c                work is a work array.  work is not referenced if
c                job .eq. 0.
c
c        job     integer.
c                job is an integer that initiates column pivoting.
c                if job .eq. 0, no pivoting is done.
c                if job .ne. 0, pivoting is done.
c
c     on return
c
c        x       x contains in its upper triangle the upper
c                triangular matrix r of the qr factorization.
c                below its diagonal x contains information from
c                which the orthogonal part of the decomposition
c                can be recovered.  note that if pivoting has
c                been requested, the decomposition is not that
c                of the original matrix x but that of x
c                with its columns permuted as described by jpvt.
c
c        qraux   double precision(p).
c                qraux contains further information required to recover
c                the orthogonal part of the decomposition.
c
c        jpvt    jpvt(k) contains the index of the column of the
c                original matrix that has been interchanged into
c                the k-th column, if pivoting was requested.
c
c     linpack. this version dated 08/14/78 .
c     g.w. stewart, university of maryland, argonne national lab.
c
c     dqrdc uses the following functions and subprograms.
c
c     blas daxpy,ddot,dscal,dswap,dnrm2
c     fortran dabs,dmax1,min0,dsqrt
c
c     internal variables
c
      integer j,jp,l,lp1,lup,maxj,pl,pu
      double precision maxnrm,dnrm2,tt
      double precision ddot,nrmxl,t
      logical negj,swapj
c
c
      pl = 1
      pu = 0
      if (job .eq. 0) go to 60
c
c        pivoting has been requested.  rearrange the columns
c        according to jpvt.
c
         do 20 j = 1, p
            swapj = jpvt(j) .gt. 0
            negj = jpvt(j) .lt. 0
            jpvt(j) = j
            if (negj) jpvt(j) = -j
            if (.not.swapj) go to 10
               if (j .ne. pl) call dswap(n,x(1,pl),1,x(1,j),1)
               jpvt(j) = jpvt(pl)
               jpvt(pl) = j
               pl = pl + 1
   10       continue
   20    continue
         pu = p
         do 50 jj = 1, p
            j = p - jj + 1
            if (jpvt(j) .ge. 0) go to 40
               jpvt(j) = -jpvt(j)
               if (j .eq. pu) go to 30
                  call dswap(n,x(1,pu),1,x(1,j),1)
                  jp = jpvt(pu)
                  jpvt(pu) = jpvt(j)
                  jpvt(j) = jp
   30          continue
               pu = pu - 1
   40       continue
   50    continue
   60 continue
c
c     compute the norms of the free columns.
c
      if (pu .lt. pl) go to 80
      do 70 j = pl, pu
         qraux(j) = dnrm2(n,x(1,j),1)
         work(j) = qraux(j)
   70 continue
   80 continue
c
c     perform the householder reduction of x.
c
      lup = min0(n,p)
      do 200 l = 1, lup
         if (l .lt. pl .or. l .ge. pu) go to 120
c
c           locate the column of largest norm and bring it
c           into the pivot position.
c
            maxnrm = 0.0d0
            maxj = l
            do 100 j = l, pu
               if (qraux(j) .le. maxnrm) go to 90
                  maxnrm = qraux(j)
                  maxj = j
   90          continue
  100       continue
            if (maxj .eq. l) go to 110
               call dswap(n,x(1,l),1,x(1,maxj),1)
               qraux(maxj) = qraux(l)
               work(maxj) = work(l)
               jp = jpvt(maxj)
               jpvt(maxj) = jpvt(l)
               jpvt(l) = jp
  110       continue
  120    continue
         qraux(l) = 0.0d0
         if (l .eq. n) go to 190
c
c           compute the householder transformation for column l.
c
            nrmxl = dnrm2(n-l+1,x(l,l),1)
            if (nrmxl .eq. 0.0d0) go to 180
               if (x(l,l) .ne. 0.0d0) nrmxl = dsign(nrmxl,x(l,l))
               call dscal(n-l+1,1.0d0/nrmxl,x(l,l),1)
               x(l,l) = 1.0d0 + x(l,l)
c
c              apply the transformation to the remaining columns,
c              updating the norms.
c
               lp1 = l + 1
               if (p .lt. lp1) go to 170
               do 160 j = lp1, p
                  t = -ddot(n-l+1,x(l,l),1,x(l,j),1)/x(l,l)
                  call daxpy(n-l+1,t,x(l,l),1,x(l,j),1)
                  if (j .lt. pl .or. j .gt. pu) go to 150
                  if (qraux(j) .eq. 0.0d0) go to 150
                     tt = 1.0d0 - (dabs(x(l,j))/qraux(j))**2
                     tt = dmax1(tt,0.0d0)
                     t = tt
                     tt = 1.0d0 + 0.05d0*tt*(qraux(j)/work(j))**2
                     if (tt .eq. 1.0d0) go to 130
                        qraux(j) = qraux(j)*dsqrt(t)
                     go to 140
  130                continue
                        qraux(j) = dnrm2(n-l,x(l+1,j),1)
                        work(j) = qraux(j)
  140                continue
  150             continue
  160          continue
  170          continue
c
c              save the transformation.
c
               qraux(l) = x(l,l)
               x(l,l) = -nrmxl
  180       continue
  190    continue
  200 continue
      return
      end
      subroutine dqrls ( a, lda, m, n, tol, kr, b, x, rsd, jpvt, 
     &  qraux, work, itask, ind )

c*********************************************************************72
c
cc DQRLS factors and solves a linear system in the least squares sense.
c
c  Discussion:
c
c    The linear system may be overdetermined, underdetermined or singular.
c    The solution is obtained using a QR factorization of the
c    coefficient matrix.
c
c    DQRLS can be efficiently used to solve several least squares
c    problems with the same matrix A.  The first system is solved
c    with ITASK = 1.  The subsequent systems are solved with
c    ITASK = 2, to avoid the recomputation of the matrix factors.
c    The parameters KR, JPVT, and QRAUX must not be modified
c    between calls to DQRLS.
c
c    DQRLS is used to solve in a least squares sense
c    overdetermined, underdetermined and singular linear systems.
c    The system is A*X approximates B where A is M by N.
c    B is a given M-vector, and X is the N-vector to be computed.
c    A solution X is found which minimimzes the sum of squares (2-norm)
c    of the residual,  A*X - B.
c
c    The numerical rank of A is determined using the tolerance TOL.
c
c    DQRLS uses the LINPACK subroutine DQRDC to compute the QR
c    factorization, with column pivoting, of an M by N matrix A.
c
c  Modified:
c
c    15 April 2012
c
c  Reference:
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input/output, double precision A(LDA,N), an M by N matrix.
c    On input, the matrix whose decomposition is to be computed.
c    In a least squares data fitting problem, A(I,J) is the
c    value of the J-th basis (model) function at the I-th data point.
c    On output, A contains the output from DQRDC.  The triangular matrix R
c    of the QR factorization is contained in the upper triangle and
c    information needed to recover the orthogonal matrix Q is stored
c    below the diagonal in A and in the vector QRAUX.
c
c    Input, integer LDA, the leading dimension of A.
c
c    Input, integer M, the number of rows of A.
c
c    Input, integer N, the number of columns of A.
c
c    Input, double precision TOL, a relative tolerance used to determine the
c    numerical rank.  The problem should be scaled so that all the elements
c    of A have roughly the same absolute accuracy EPS.  Then a reasonable
c    value for TOL is roughly EPS divided by the magnitude of the largest
c    element.
c
c    Output, integer KR, the numerical rank.
c
c    Input, double precision B(M), the right hand side of the linear system.
c
c    Output, double precision X(N), a least squares solution to the linear
c    system.
c
c    Output, double precision RSD(M), the residual, B - A*X.  RSD may
c    overwrite B.
c
c    Workspace, integer JPVT(N), required if ITASK = 1.
c    Columns JPVT(1), ..., JPVT(KR) of the original matrix are linearly
c    independent to within the tolerance TOL and the remaining columns
c    are linearly dependent.  ABS ( A(1,1) ) / ABS ( A(KR,KR) ) is an estimate
c    of the condition number of the matrix of independent columns,
c    and of R.  This estimate will be <= 1/TOL.
c
c    Workspace, double precision QRAUX(N), required if ITASK = 1.
c
c    Workspace, double precision WORK(N), required if ITASK = 1.
c
c    Input, integer ITASK.
c    1, DQRLS factors the matrix A and solves the least squares problem.
c    2, DQRLS assumes that the matrix A was factored with an earlier
c       call to DQRLS, and only solves the least squares problem.
c
c    Output, integer IND, error code.
c    0:  no error
c    -1: LDA < M   (fatal error)
c    -2: N < 1     (fatal error)
c    -3: ITASK < 1 (fatal error)
c
      implicit none

      integer lda
      integer m
      integer n

      double precision a(lda,n)
      double precision b(m)
      integer ind
      integer itask
      integer jpvt(n)
      integer kr
      double precision qraux(n)
      double precision rsd(m)
      double precision tol
      double precision work(n)
      double precision x(n)

      if ( lda .lt. m ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DQRLS - Fatal error!'
        write ( *, '(a)' ) '  LDA < M.'
        stop
      end if

      if ( n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DQRLS - Fatal error!'
        write ( *, '(a)' ) '  N <= 0.'
        stop
      end if

      if ( itask .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DQRLS - Fatal error!'
        write ( *, '(a)' ) '  ITASK < 1.'
        stop
      end if

      ind = 0
c
c  Factor the matrix.
c
      if ( itask .eq. 1 ) then
        call dqrank ( a, lda, m, n, tol, kr, jpvt, qraux, work )
      end if
c
c  Solve the least-squares problem.
c
      call dqrlss ( a, lda, m, n, kr, b, x, rsd, jpvt, qraux )

      return
      end
      subroutine dqrlss ( a, lda, m, n, kr, b, x, rsd, jpvt, qraux )

c*****************************************************************************80
c
cc DQRLSS solves a linear system in a least squares sense.
c
c  Discussion:
c
c    DQRLSS must be preceeded by a call to DQRANK.
c
c    The system is to be solved is
c      A * X = B
c    where
c      A is an M by N matrix with rank KR, as determined by DQRANK,
c      B is a given M-vector,
c      X is the N-vector to be computed.
c
c    A solution X, with at most KR nonzero components, is found which
c    minimizes the 2-norm of the residual (A*X-B).
c
c    Once the matrix A has been formed, DQRANK should be
c    called once to decompose it.  Then, for each right hand
c    side B, DQRLSS should be called once to obtain the
c    solution and residual.
c
c  Modified:
c
c    15 April 2012
c
c  Parameters:
c
c    Input, double precision A(LDA,N), the QR factorization information
c    from DQRANK.  The triangular matrix R of the QR factorization is
c    contained in the upper triangle and information needed to recover
c    the orthogonal matrix Q is stored below the diagonal in A and in
c    the vector QRAUX.
c
c    Input, integer LDA, the leading dimension of A, which must
c    be at least M.
c
c    Input, integer M, the number of rows of A.
c
c    Input, integer N, the number of columns of A.
c
c    Input, integer KR, the rank of the matrix, as estimated
c    by DQRANK.
c
c    Input, double precision B(M), the right hand side of the linear system.
c
c    Output, double precision X(N), a least squares solution to the
c    linear system.
c
c    Output, double precision RSD(M), the residual, B - A*X.  RSD may
c    overwite B.
c
c    Input, integer JPVT(N), the pivot information from DQRANK.
c    Columns JPVT(1), ..., JPVT(KR) of the original matrix are linearly
c    independent to within the tolerance TOL and the remaining columns
c    are linearly dependent.
c
c    Input, double precision QRAUX(N), auxiliary information from DQRANK
c    defining the QR factorization.
c
      implicit none

      integer lda
      integer m
      integer n

      double precision a(lda,n)
      double precision b(m)
      integer info
      integer j
      integer jpvt(n)
      integer k
      integer kr
      double precision qraux(n)
      double precision rsd(m)
      double precision t
      double precision x(n)

      if ( kr .ne. 0 ) then
        call dqrsl ( a, lda, m, kr, qraux, b, rsd, rsd, x, rsd, rsd, 
     &    110, info )
      end if

      do j = 1, n
        jpvt(j) = - jpvt(j)
      end do

      do j = kr + 1, n
        x(j) = 0.0D+00
      end do

      do j = 1, n

        if ( jpvt(j) .le. 0 ) then

          k = -jpvt(j)
          jpvt(j) = k

10        continue

          if ( k .ne. j ) then
            t = x(j)
            x(j) = x(k)
            x(k) = t
            jpvt(k) = -jpvt(k)
            k = jpvt(k)
            go to 10
          end if

        end if

      end do

      return
      end
      subroutine dqrsl(x,ldx,n,k,qraux,y,qy,qty,b,rsd,xb,job,info)

c*********************************************************************72

      integer ldx,n,k,job,info
      double precision x(ldx,1),qraux(1),y(1),qy(1),qty(1),b(1),rsd(1),
     *                 xb(1)
c
c     dqrsl applies the output of dqrdc to compute coordinate
c     transformations, projections, and least squares solutions.
c     for k .le. min(n,p), let xk be the matrix
c
c            xk = (x(jpvt(1)),x(jpvt(2)), ... ,x(jpvt(k)))
c
c     formed from columnns jpvt(1), ... ,jpvt(k) of the original
c     n x p matrix x that was input to dqrdc (if no pivoting was
c     done, xk consists of the first k columns of x in their
c     original order).  dqrdc produces a factored orthogonal matrix q
c     and an upper triangular matrix r such that
c
c              xk = q * (r)
c                       (0)
c
c     this information is contained in coded form in the arrays
c     x and qraux.
c
c     on entry
c
c        x      double precision(ldx,p).
c               x contains the output of dqrdc.
c
c        ldx    integer.
c               ldx is the leading dimension of the array x.
c
c        n      integer.
c               n is the number of rows of the matrix xk.  it must
c               have the same value as n in dqrdc.
c
c        k      integer.
c               k is the number of columns of the matrix xk.  k
c               must nnot be greater than min(n,p), where p is the
c               same as in the calling sequence to dqrdc.
c
c        qraux  double precision(p).
c               qraux contains the auxiliary output from dqrdc.
c
c        y      double precision(n)
c               y contains an n-vector that is to be manipulated
c               by dqrsl.
c
c        job    integer.
c               job specifies what is to be computed.  job has
c               the decimal expansion abcde, with the following
c               meaning.
c
c                    if a.ne.0, compute qy.
c                    if b,c,d, or e .ne. 0, compute qty.
c                    if c.ne.0, compute b.
c                    if d.ne.0, compute rsd.
c                    if e.ne.0, compute xb.
c
c               note that a request to compute b, rsd, or xb
c               automatically triggers the computation of qty, for
c               which an array must be provided in the calling
c               sequence.
c
c     on return
c
c        qy     double precision(n).
c               qy conntains q*y, if its computation has been
c               requested.
c
c        qty    double precision(n).
c               qty contains trans(q)*y, if its computation has
c               been requested.  here trans(q) is the
c               transpose of the matrix q.
c
c        b      double precision(k)
c               b contains the solution of the least squares problem
c
c                    minimize norm2(y - xk*b),
c
c               if its computation has been requested.  (note that
c               if pivoting was requested in dqrdc, the j-th
c               component of b will be associated with column jpvt(j)
c               of the original matrix x that was input into dqrdc.)
c
c        rsd    double precision(n).
c               rsd contains the least squares residual y - xk*b,
c               if its computation has been requested.  rsd is
c               also the orthogonal projection of y onto the
c               orthogonal complement of the column space of xk.
c
c        xb     double precision(n).
c               xb contains the least squares approximation xk*b,
c               if its computation has been requested.  xb is also
c               the orthogonal projection of y onto the column space
c               of x.
c
c        info   integer.
c               info is zero unless the computation of b has
c               been requested and r is exactly singular.  in
c               this case, info is the index of the first zero
c               diagonal element of r and b is left unaltered.
c
c     the parameters qy, qty, b, rsd, and xb are not referenced
c     if their computation is not requested and in this case
c     can be replaced by dummy variables in the calling program.
c     to save storage, the user may in some cases use the same
c     array for different parameters in the calling sequence.  a
c     frequently occuring example is when one wishes to compute
c     any of b, rsd, or xb and does not need y or qty.  in this
c     case one may identify y, qty, and one of b, rsd, or xb, while
c     providing separate arrays for anything else that is to be
c     computed.  thus the calling sequence
c
c          call dqrsl(x,ldx,n,k,qraux,y,dum,y,b,y,dum,110,info)
c
c     will result in the computation of b and rsd, with rsd
c     overwriting y.  more generally, each item in the following
c     list contains groups of permissible identifications for
c     a single callinng sequence.
c
c          1. (y,qty,b) (rsd) (xb) (qy)
c
c          2. (y,qty,rsd) (b) (xb) (qy)
c
c          3. (y,qty,xb) (b) (rsd) (qy)
c
c          4. (y,qy) (qty,b) (rsd) (xb)
c
c          5. (y,qy) (qty,rsd) (b) (xb)
c
c          6. (y,qy) (qty,xb) (b) (rsd)
c
c     in any group the value returned in the array allocated to
c     the group corresponds to the last member of the group.
c
c     linpack. this version dated 08/14/78 .
c     g.w. stewart, university of maryland, argonne national lab.
c
c     dqrsl uses the following functions and subprograms.
c
c     blas daxpy,dcopy,ddot
c     fortran dabs,min0,mod
c
c     internal variables
c
      integer i,j,jj,ju,kp1
      double precision ddot,t,temp
      logical cb,cqy,cqty,cr,cxb
c
c
c     set info flag.
c
      info = 0
c
c     determine what is to be computed.
c
      cqy = job/10000 .ne. 0
      cqty = mod(job,10000) .ne. 0
      cb = mod(job,1000)/100 .ne. 0
      cr = mod(job,100)/10 .ne. 0
      cxb = mod(job,10) .ne. 0
      ju = min0(k,n-1)
c
c     special action when n=1.
c
      if (ju .ne. 0) go to 40
         if (cqy) qy(1) = y(1)
         if (cqty) qty(1) = y(1)
         if (cxb) xb(1) = y(1)
         if (.not.cb) go to 30
            if (x(1,1) .ne. 0.0d0) go to 10
               info = 1
            go to 20
   10       continue
               b(1) = y(1)/x(1,1)
   20       continue
   30    continue
         if (cr) rsd(1) = 0.0d0
      go to 250
   40 continue
c
c        set up to compute qy or qty.
c
         if (cqy) call dcopy(n,y,1,qy,1)
         if (cqty) call dcopy(n,y,1,qty,1)
         if (.not.cqy) go to 70
c
c           compute qy.
c
            do 60 jj = 1, ju
               j = ju - jj + 1
               if (qraux(j) .eq. 0.0d0) go to 50
                  temp = x(j,j)
                  x(j,j) = qraux(j)
                  t = -ddot(n-j+1,x(j,j),1,qy(j),1)/x(j,j)
                  call daxpy(n-j+1,t,x(j,j),1,qy(j),1)
                  x(j,j) = temp
   50          continue
   60       continue
   70    continue
         if (.not.cqty) go to 100
c
c           compute trans(q)*y.
c
            do 90 j = 1, ju
               if (qraux(j) .eq. 0.0d0) go to 80
                  temp = x(j,j)
                  x(j,j) = qraux(j)
                  t = -ddot(n-j+1,x(j,j),1,qty(j),1)/x(j,j)
                  call daxpy(n-j+1,t,x(j,j),1,qty(j),1)
                  x(j,j) = temp
   80          continue
   90       continue
  100    continue
c
c        set up to compute b, rsd, or xb.
c
         if (cb) call dcopy(k,qty,1,b,1)
         kp1 = k + 1
         if (cxb) call dcopy(k,qty,1,xb,1)
         if (cr .and. k .lt. n) call dcopy(n-k,qty(kp1),1,rsd(kp1),1)
         if (.not.cxb .or. kp1 .gt. n) go to 120
            do 110 i = kp1, n
               xb(i) = 0.0d0
  110       continue
  120    continue
         if (.not.cr) go to 140
            do 130 i = 1, k
               rsd(i) = 0.0d0
  130       continue
  140    continue
         if (.not.cb) go to 190
c
c           compute b.
c
            do 170 jj = 1, k
               j = k - jj + 1
               if (x(j,j) .ne. 0.0d0) go to 150
                  info = j
c           ......exit
                  go to 180
  150          continue
               b(j) = b(j)/x(j,j)
               if (j .eq. 1) go to 160
                  t = -b(j)
                  call daxpy(j-1,t,x(1,j),1,b,1)
  160          continue
  170       continue
  180       continue
  190    continue
         if (.not.cr .and. .not.cxb) go to 240
c
c           compute rsd or xb as required.
c
            do 230 jj = 1, ju
               j = ju - jj + 1
               if (qraux(j) .eq. 0.0d0) go to 220
                  temp = x(j,j)
                  x(j,j) = qraux(j)
                  if (.not.cr) go to 200
                     t = -ddot(n-j+1,x(j,j),1,rsd(j),1)/x(j,j)
                     call daxpy(n-j+1,t,x(j,j),1,rsd(j),1)
  200             continue
                  if (.not.cxb) go to 210
                     t = -ddot(n-j+1,x(j,j),1,xb(j),1)/x(j,j)
                     call daxpy(n-j+1,t,x(j,j),1,xb(j),1)
  210             continue
                  x(j,j) = temp
  220          continue
  230       continue
  240    continue
  250 continue
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

      double precision c
      double precision dtemp
      double precision dx(*)
      double precision dy(*)
      double precision s
      integer i,incx,incy,ix,iy,n

      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c  code for unequal increments or equal increments not equal to 1
c
      if(incx.lt.0)then
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
      double precision db
      double precision r
      double precision roe
      double precision s
      double precision scale
      double precision z

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
      integer incx
      integer m
      integer n
      integer nincx

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
        dx(i)     = da * dx(i)
        dx(i + 1) = da * dx(i + 1)
        dx(i + 2) = da * dx(i + 2)
        dx(i + 3) = da * dx(i + 3)
        dx(i + 4) = da * dx(i + 4)
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
      if(incx.lt.0) then
        ix = (-n+1)*incx + 1
      end if
      if(incy.lt.0) then
        iy = (-n+1)*incy + 1
      end if
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
      subroutine normal_solve ( m, n, a, b, x, flag )

c*********************************************************************72
c
cc NORMAL_SOLVE solves a linear system using the normal equations.
c
c  Discussion:
c
c    Given a presumably rectangular MxN system of equations A*x=b, this routine
c    sets up the NxN system A'*A*x=A'b.  Assuming N <= M, and that A has full
c    column rank, the system will be solvable, and the vector x that is returned
c    will minimize the Euclidean norm of the residual.
c
c    One drawback to this approach is that the condition number of the linear
c    system A'*A is effectively the square of the condition number of A, 
c    meaning that there is a substantial loss of accuracy.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 April 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, integer M, the number of rows of A.
c
c    Input, integer N, the number of columns of A.
c    It must be the case that N <= M.
c
c    Input, double precision A(M,N), the matrix.
c    The matrix must have full column rank.
c
c    Input, double precision B(M), the right hand side.
c
c    Output, double precision X(N), the least squares solution.
c
c    Output, integer FLAG,
c    0, no error was detected.
c    1, an error occurred.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision ata(n,n)
      double precision ata_c(n,n)
      double precision atb(n)
      double precision b(m)
      integer flag
      double precision x(n)

      flag = 0

      if ( m .lt. n ) then
        flag = 1
        return
      end if

      call r8mat_mtm ( n, m, n, a, a, ata )

      call r8mat_mtv ( m, n, a, b, atb )

      call r8mat_cholesky_factor ( n, ata, ata_c, flag )

      if ( flag .ne. 0 ) then
        return
      end if

      call r8mat_cholesky_solve ( n, ata_c, atb, x )

      return
      end
      subroutine qr_solve ( m, n, a, b, x )

c*********************************************************************72
c
cc QR_SOLVE solves a linear system in the least squares sense.
c
c  Discussion:
c
c    If the matrix A has full column rank, then the solution X should be the
c    unique vector that minimizes the Euclidean norm of the residual.
c
c    If the matrix A does not have full column rank, then the solution is
c    not unique; the vector X will minimize the residual norm, but so will
c    various other vectors.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 April 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, integer M, the number of rows of A.
c
c    Input, integer N, the number of columns of A.
c
c    Input, double precision A(M,N), the matrix.
c
c    Input, double precision B(M), the right hand side.
c
c    Output, double precision X(N), the least squares solution.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision a_qr(m,n)
      double precision amax
      double precision b(m)
      integer i
      integer ind
      integer itask
      integer j
      integer jpvt(n)
      integer kr
      integer lda
      double precision qraux(n)
      double precision r(m)
      double precision r8_epsilon
      double precision tol
      double precision work(n)
      double precision x(n)

      do j = 1, n
        do i = 1, m
          a_qr(i,j) = a(i,j)
        end do
      end do

      lda = m

      call r8mat_amax ( m, n, a_qr, amax )

      tol = r8_epsilon ( ) / amax
      
      itask = 1

      call dqrls ( a_qr, lda, m, n, tol, kr, b, x, r, jpvt, qraux, 
     &  work, itask, ind )

      return
      end
      subroutine svd_solve ( m, n, a, b, x )

c*********************************************************************72
c
cc SVD_SOLVE solves a linear system in the least squares sense.
c
c  Discussion:
c
c    The vector X returned by this routine should always minimize the 
c    Euclidean norm of the residual ||A*x-b||.
c
c    If the matrix A does not have full column rank, then there are multiple
c    vectors that attain the minimum residual.  In that case, the vector
c    X returned by this routine is the unique such minimizer that has the 
c    the minimum possible Euclidean norm, that is, ||A*x-b|| and ||x||
c    are both minimized.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 April 2012
c
c  Reference:
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, integer M, the number of rows of A.
c
c    Input, integer N, the number of columns of A.
c
c    Input, double precision A(M,N), the matrix.
c
c    Input, double precision B(M), the right hand side.
c
c    Output, double precision X(N), the least squares solution.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision a_copy(m,n)
      double precision b(m)
      double precision e(max(m+1,n))
      integer i
      integer info
      integer j
      integer lda
      integer ldu
      integer ldv
      integer job
      double precision r8_epsilon
      double precision sdiag(max(m+1,n))
      double precision smax
      double precision stol
      double precision sub(n)
      double precision u(m,m)
      double precision ub(m)
      double precision v(n,n)
      double precision work(m)
      double precision x(n)
c
c  Get the SVD.
c
      do j = 1, n
        do i = 1, m
          a_copy(i,j) = a(i,j)
        end do
      end do

      lda = m
      ldu = m
      ldv = n
      job = 11

      call dsvdc ( a_copy, lda, m, n, sdiag, e, u, ldu, v, ldv, work, 
     &  job, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SVD_SOLVE - Failure!'
        write ( *, '(a)' ) '  The SVD could not be calculated.'
        write ( *, '(a)' ) '  LINPACK routine DSVDC returned a nonzero'
        write ( *, '(a,i8)' ) '  value of the error flag, INFO = ', info
        stop
      end if

      call r8mat_mtv ( m, m, u, b, ub )

      do i = 1, n
        sub(i) = 0.0D+00
      end do
c
c  For singular problems, there may be tiny but nonzero singular values
c  that should be ignored.  This is a reasonable attempt to avoid such 
c  problems, although in general, the user might wish to control the tolerance.
c
      call r8vec_max ( n, sdiag, smax )

      if ( smax .le. r8_epsilon ( ) ) then
        smax = 1.0D+00
      end if

      stol = r8_epsilon ( ) * smax

      do i = 1, n
        if ( i .le. m ) then
          if ( stol .le. sdiag(i) ) then
            sub(i) = ub(i) / sdiag(i)
          end if
        end if
      end do

      call r8mat_mv ( n, n, v, sub, x )

      return
      end
