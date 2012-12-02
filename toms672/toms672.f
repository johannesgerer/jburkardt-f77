      subroutine assign ( n, pnodes, iwork, work, ldw, recur, t, ierr )

c*********************************************************************72
c
cc ASSIGN generates the polynomial whose roots are the preassigned nodes.
c
c  Discussion:
c
c    This routine generates the initial polynomial T whose roots are the required
c    pre-assigned nodes.
c
c  Modified:
c
c    02 March 2009
c
c  Author:
c
c    Original FORTRAN77 version by Thomas Patterson.
c
c  Reference:
c
c    Thomas Patterson,
c    An algorithm for generating interpolatory quadrature rules of the highest
c    degree of precision with preassigned nodes for general weight functions,
c    Transactions on Mathematical Software,
c    Volume 15, Number 2, June 1989, pages 123-136.
c
c    Thomas Patterson,
c    Algorithm 672:
c    EXTEND: generation of interpolatory quadrature rules of the highest degree
c    of precision with preassigned nodes for general weight functions,
c    Transactions on Mathematical Software,
c    Volume 15, Number 2, June 1989, pages 137-143.
c
c  Parameters:
c
c    Input, integer N, the number of pre-assigned nodes to be used to
c    generate T.
c
c    Input, double precision PNODES(N), the pre-assigned nodes to be be used to
c    generate T.
c
c    Workspace, integer IWORK(N).
c
c    Workspace, double precision WORK(N,N+1).
c
c    Input, integer LDW, the leading dimension of WORK.
c
c    Input, external RECUR ( ), the user supplied function which defines the
c    orthogonal polynomials.  Given K, CALL RECUR(K,C,D,E) gives
c    the coefficients C,D and E such that,
c      P(K+1,X) = (C*X+D)*P(K,X)+E*P(K-1,X)
c    The parameters are defined as follows:
c      K = Index
c      C, D, E = parameters in the recurrence relation (functions of K).
c
c    Output, double precision T(0:N), the coefficients of the polynomial whose
c    roots define the pre-assigned nodes of the quadrature rule and expressed
c    as:
c      H0* SUM (I = 0 to N) T(I)/HI*P(I,X)
c    T(I) holds the value of TI.
c
c    Output, integer IERR
c    * 0, No error detected
c    * 1, The linear system of equations used to generate the polynomial T
c      became singular or very ill-conditioned.
c
      implicit none

      integer ldw
      integer n

      double precision c0
      double precision d0
      double precision e0
      integer ierr
      integer info
      integer iwork(n)
      integer j
      integer k
      double precision p
      double precision p0
      double precision p1
      double precision pnodes(n)
      external recur
      double precision t(0:n)
      double precision work(0:ldw-1,0:n)
      double precision x

      ierr = 0
c
c  Set up the linear system of equations.
c
      do j = 1, n

        x = pnodes(j)
        p0 = 0.0D+00
        p1 = 1.0D+00
        p = p1

        do k = 0, n
          work(j-1,k) = p
          call recur ( k, c0, d0, e0 )
          p = ( c0 * x + d0 ) * p1 + e0 * p0
          p0 = p1
          p1 = p
        end do

      end do
c
c  Solve linear system.
c
      call dgefa ( work, ldw, n, iwork, info )

      if ( info .ne. 0 ) then
        ierr = 1
        return
      end if

      call dgesl ( work, ldw, n, iwork, work(0,n), 0 )

      do j = 0, n - 1
        t(j) = - work(j,n)
      end do

      t(n) = 1.0D+00
c
c  Weight with moments.
c
      call transf ( t, 0, n, recur, 1 )

      return
      end
      subroutine bair ( n, polin, polout, a0, a1, a2, recur, idigit,
     &  errval, ifail )

c*********************************************************************72
c
cc BAIR seeks roots of a polynomial.
c
c  Discussion:
c
c    This function carries out a generalized Bairstow root extraction
c    for the polynomial:
c
c      SUM(I = 0 to N)  POLIN(I) * P(I,X).
c
c    It calculates the root as a quadratic factor:
c
c      A2 * P(2,X) - A1 * P(1,X) - A0 * P(0,X)
c
c    where P(I,X) is a general orthogonal polynomial of degree I.
c
c  Modified:
c
c    27 February 2009
c
c  Author:
c
c    Original FORTRAN77 version by Thomas Patterson.
c
c  Reference:
c
c    Gene Golub, Thomas Robertson,
c    A generalized Bairstow Algorithm,
c    Communications of the ACM,
c    Volume 10, Number 6, June 1967, pages 371-373.
c
c    Thomas Patterson,
c    An algorithm for generating interpolatory quadrature rules of the highest
c    degree of precision with preassigned nodes for general weight functions,
c    Transactions on Mathematical Software,
c    Volume 15, Number 2, June 1989, pages 123-136.
c
c    Thomas Patterson,
c    Algorithm 672:
c    EXTEND: generation of interpolatory quadrature rules of the highest degree
c    of precision with preassigned nodes for general weight functions,
c    Transactions on Mathematical Software,
c    Volume 15, Number 2, June 1989, pages 137-143.
c
c  Parameters:
c
c    Input, integer N, the degree of input polynomial POLIN.
c
c    Input, double precision POLIN(0:N), coefficients of the polynomial whose
c    quadratic factor is to be found, i.e.
c      POLIN = SUM(I = 0 to N) POLIN(I)*P(I,X)
c
c    Output, double precision POLOUT(0:N-2), coefficients of the deflated
c    polynomial of degree N-2 with the quadratic factor removed, i.e.
c      POLOUT = SUM(I = 0 to N-2) POLOUT(I)*P(I,X)
c
c    Input/output, double precision A0, A1, A2, on input, the estimated
c    quadratic factors.  On output, the estimate has been improved.
c
c    Input, external RECUR ( ), the function which defines the orthogonal
c    polynomials.  See EXTEND for full description.
c
c    Input, integer IDIGIT, the node convergence parameter, an integer greater
c    than 0.  An attempt is made to calculate the nodes to the maximum
c    accuracy possible by the machine precision available.
c    IDIGIT controls the assessment procedure to take account of
c    round-off errors and specifies the number of least signific
c    decimal digits that can be ignored (i.e. attributed
c    to round-off) in the computed relative error.  A typical value is 5.
c
c    Output, double precision ERRVAL, the mean value of the correction to
c    the coefficients of the quadratic factor. May be used as a measure of the
c    root accuracy when convergence is not achieved.
c
c    Output, integer IFAIL, error flag.
c    * 0, Quadratic factor found.
c    * 1, Convergence not achieved after 50 iterations.
c
      implicit none

      integer n

      double precision a
      double precision a0
      double precision a1
      double precision a2
      double precision aa
      double precision ab
      double precision alpha
      double precision b
      double precision ba
      double precision bb
      double precision beta
      double precision delta
      double precision eps
      double precision errval
      double precision f1
      double precision f2
      integer idigit
      integer ifail
      integer iter
      double precision polin(0:n)
      double precision polout(0:n-2)
      external recur
      double precision scale
      double precision tol

      ifail = 0
      iter = 50
      errval = 0.0D+00
c
c  Special cases.
c
      if ( n .eq. 1 ) then
        a0 = - polin(0)
        a1 = - polin(1)
        a2 = 0.0D+00
        return
      end if

      if ( n .eq. 2 ) then
        a0 = - polin(0)
        a1 = - polin(1)
        a2 =  polin(2)
        return
      end if
c
c  Estimated ALPHA and BETA.
c
      tol = 10.0D+00**( - max ( 1, idigit ) )
      alpha = a1 / a2
      beta = a0 / a2

10    continue

      iter = iter - 1

      if ( iter .lt. 0 ) then
        ifail = 1
        goto 20
      end if

      call qfact ( n, polin, polout, recur, alpha, beta, a, b, aa,
     &  ab, ba, bb )

      scale = max ( abs ( ab ), abs ( bb ) )
      f1 = ab / scale
      f2 = bb / scale
      delta = ( b * f1 - a * f2 ) / ( aa * f2 - ba * f1 )
      scale = max ( abs ( ba ), abs ( aa ) )
      f1 = ba / scale
      f2 = aa / scale
      eps = ( a * f1 - b * f2 ) / ( bb * f2 - ab * f1 )
      alpha = alpha + delta
      beta = beta + eps
c
c  Test for convergence.
c  Stop if correction is less than 1/TOL times the smallest machine
c  relative error.
c
      if ( abs ( alpha ) + tol * abs ( delta ) .ne. abs ( alpha ) .or.
     &     abs ( beta ) + tol * abs ( eps ) .ne. abs ( beta ) ) then
        go to 10
      end if
c
c  Final iteration to tidy up result.
c
      call qfact ( n, polin, polout, recur, alpha, beta, a, b, aa,
     &  ab, ba, bb )

      scale = max ( abs ( ab ), abs ( bb ) )
      f1 = ab / scale
      f2 = bb / scale
      delta = ( b * f1 - a * f2 ) / ( aa * f2 - ba * f1 )
      scale = max ( abs ( ba ), abs ( aa ) )
      f1 = ba / scale
      f2 = aa / scale
      eps = ( a * f1 - b * f2 ) / ( bb * f2 - ab * f1 )
      alpha = alpha + delta
      beta = beta + eps

20    continue

      a0 = beta
      a1 = alpha
      a2 = 1.0D+00
      errval = 0.5D+00 * ( abs ( eps ) + abs ( delta ) )

      return
      end
      subroutine check ( n, pnodes, wt, k, h0, recur, test, ierr )

c*********************************************************************72
c
cc CHECK tests a computed quadrature rule.
c
c  Purpose:
c
c    This function carries out a test of the given quadrature rule by
c    computing the appropriate integral of
c
c      W(X) * P(K,X) * P(K,X)
c
c    over the region associated with the weight function W(X) and the
c    orthogonal polynomials P(K,X).
c
c  Modified:
c
c    02 March 2009
c
c  Author:
c
c    Original FORTRAN77 version by Thomas Patterson.
c
c  Reference:
c
c    Thomas Patterson,
c    An algorithm for generating interpolatory quadrature rules of the highest
c    degree of precision with preassigned nodes for general weight functions,
c    Transactions on Mathematical Software,
c    Volume 15, Number 2, June 1989, pages 123-136.
c
c    Thomas Patterson,
c    Algorithm 672:
c    EXTEND: generation of interpolatory quadrature rules of the highest degree
c    of precision with preassigned nodes for general weight functions,
c    Transactions on Mathematical Software,
c    Volume 15, Number 2, June 1989, pages 137-143.
c
c  Parameters:
c
c    Input, integer N, the number of nodes in the quadrature rule.
c
c    Input, double precision PNODES(N), the nodes.
c
c    Input, double precision WT(N), the weights.
c
c    Input, integer K, the index of the orthogonal polynomial whose weighted
c    square is to be integrated.
c
c    Input, double precision H0, the integral of the orthogonality weight
c    function over the interval of integration. Zero moment integral.  Note that
c    P(0,X) is arbitrarily taken to be 1.0.
c
c    Input, external RECUR (), the function which defines the
c    orthogonal polynomials.  See EXTEND for a full description.
c
c    Output, double precision TEST, the approximate value of the test integral
c    normalized to unity.  Thus, abs(TEST-1) gives a measure of the
c    quality of the calculated rule.
c
c    Output, integer IERR, error flag.
c    * 0, OK.
c    * 1, Rule quality unsatisfactory.
c    * 2, Invalid values for input arguments.
c
      implicit none

      integer n

      double precision cj
      double precision dj
      double precision ej
      double precision h0
      integer i
      integer ierr
      integer j
      integer k
      double precision p
      double precision p0
      double precision p1
      double precision pnodes(n)
      external recur
      double precision test
      double precision tol
      parameter ( tol = 0.0000001D+00 )
      double precision wt(n)
      double precision x

      ierr = 0

      if ( k .lt. 0 ) then
        ierr = 2
        return
      end if

      if ( n .lt. 1 ) then
        ierr = 2
        return
      end if

      if ( h0 .le. 0.0D+00 ) then
        ierr = 2
        return
      end if

      test = 0.0D+00

      do i = 1, n

        p1 = 1.0D+00

        if ( 0 .lt. k ) then

          p0 = 0.0D+00
          x = pnodes(i)
c
c  Calculate the integrand.
c
          do j = 0, k - 1
            call recur ( j, cj, dj, ej )
            p = ( cj * x + dj ) * p1 + ej * p0
            p0 = p1
            p1 = p
          end do

        end if

        test = test + p1 * p1 * wt(i)

      end do

      test = test / h0

      if ( k .eq. 0 ) then
        return
      end if
c
c  Calculate the exact value.
c
      j = 0
      call recur ( j, p, p0, p1 )

      do j = 1, k
        call recur ( j, cj, dj, ej )
        p = - p * ej
      end do
c
c  Normalize the result to unity.
c
      test = test * cj / p
c
c  Test for rule quality.
c
      if ( tol .lt. abs ( test - 1.0D+00 ) ) then
        ierr = 1
      end if

      return
      end
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

      if ( da .eq. 0.0d0 ) then
        return
      end if

      if ( incx .ne. 1 .or. incy .ne. 1 ) then

        ix = 1
        iy = 1
        if ( incx .lt. 0 ) ix = (-n+1)*incx + 1
        if ( incy .lt. 0 ) iy = (-n+1)*incy + 1

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
      integer i,incx,incy,ix,iy,m,n

      ddot = 0.0d0
      dtemp = 0.0d0
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
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
      end do
      ddot = dtemp
      return
c
c  code for both increments equal to 1
c
c
c  clean-up loop
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
      subroutine dgefa(a,lda,n,ipvt,info)

c*********************************************************************72
c
cc DGEFA factors a double precision matrix by gaussian elimination.
c
c     dgefa is usually called by dgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
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
c                     indicate that dgesl or dgedi will divide by zero
c                     if called.  use  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
      integer lda,n,ipvt(1),info
      double precision a(lda,1)
      double precision t
      integer idamax,j,k,kp1,l,nm1
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end
      subroutine dgesl(a,lda,n,ipvt,b,job)

c*********************************************************************72
c
cc DGESL solves a linear system factored by DGEFA.
c
c  DGESL solves the double precision system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
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
c        called correctly and if dgeco has set rcond .gt. 0.0
c        or dgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call dgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
      integer lda,n,ipvt(1),job
      double precision a(lda,1),b(1)
      double precision ddot,t
      integer k,kb,l,nm1
c
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call daxpy(k-1,t,a(1,k),1,b(1),1)
         end do
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do k = 1, n
            t = ddot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
         end do
c
c        now solve trans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
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
      integer i,incx,m,n,nincx

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
      subroutine eprod ( n, j, coeff, work, lw, recur, ifail )

c*********************************************************************72
c
cc EPROD expands a product of two orthogonal polynomials.
c
c  Discussion:
c
c    This function calculates the expansion of a product of two orthogonal
c    polynomials:
c
c      P(N,X) * P(J,X) = SUM (I = N-J to N+J ) COEFF(I) * P(I,X)
c
c    where J must not exceed N.  The orthogonal polynomials are defined
c    by the recurrence relation calculated by RECUR.
c
c    For proper initialization, the function must first be called
c    with J = 0 and the required value of N.  Subsequent calls must be in
c    the order J = 1,2,,,,,N with the appropriate expansion being
c    generated from previous values and returned in COEFF(*).  The
c    coefficients of P(N-J,X),...., P(N+J,X) are stored in the array
c    COEFF(1),...,COEFF(2*J+1).
c
c  Modified:
c
c    02 March 2009
c
c  Author:
c
c    Original FORTRAN77 version by Thomas Patterson.
c
c  Reference:
c
c    Thomas Patterson,
c    An algorithm for generating interpolatory quadrature rules of the highest
c    degree of precision with preassigned nodes for general weight functions,
c    Transactions on Mathematical Software,
c    Volume 15, Number 2, June 1989, pages 123-136.
c
c    Thomas Patterson,
c    Algorithm 672:
c    EXTEND: generation of interpolatory quadrature rules of the highest degree
c    of precision with preassigned nodes for general weight functions,
c    Transactions on Mathematical Software,
c    Volume 15, Number 2, June 1989, pages 137-143.
c
c  Parameters:
c
c    Input, integer N, the highest polynomial degree. Note that after the initial
c    call with J = 0 the value of N in this argument is ignored.
c
c    Input, integer J, the current product of P(J,X) with P(N,X) to be
c    calculated.  Note that this function must be first called with J = 0 and
c    the required largest N.  Subsequent calls must be
c    in the order J = 1,2,..,N.
c
c    Output, double precision COEFF(2*J+1), the coefficients of the expansion.
c
c    Workspace, double precision WORK(2*J+1,2).
c    The contents of this work area must not be altered between
c    calls by the calling program.
c
c    Input, integer LW, leading dimension of WORK in the calling program
c
c    Input, external RECUR ( ), the function which defines the
c    orthogonal polynomials.  See EXTEND for a full description.
c
c    Output, integer IFAIL, error flag.
c    * 0, Result OK.
c    * 1, J exceeds N.
c    * 2, J has not been called sequentially.
c
      implicit none

      integer j
      integer lw

      double precision ci
      double precision cj1
      double precision coeff(2*j+1)
      double precision di
      double precision dj1
      double precision ei
      double precision ej1
      integer i
      integer ibot
      integer ifail
      integer ii
      integer itop
      integer ix(2)
      integer j2
      integer last
      integer n
      external recur
      integer s
      integer ss
      double precision work(lw,2)

      save ix
      save last
      save ss

      ifail = 0
c
c  Initialize.
c
      if ( j .eq. 0 ) then
        ix(1) = 1
        ix(2) = 2
        coeff(1) = 1.0D+00
        work(1,2) = 1.0D+00
        last = 0
        ss = n
        return
      end if

      s = ss
c
c  Check that J does not exceed S value.
c
      if ( s .lt. j ) then
        ifail = 1
        return
      end if
c
c  Check that J is used sequentially.
c
      if ( last .ne. j - 1 ) then
        ifail = 2
        return
      end if

      last = j
      j2 = j + j
      call recur ( j - 1, cj1, dj1, ej1 )

      if ( j .eq. 1 ) then

        do i = 1, j2 + 1
          coeff(i) = 0.0D+00
        end do

      else

        do i = 1, j2 - 3
          coeff(i+2) = work(i,ix(1)) * ej1
        end do

        coeff(1)    = 0.0D+00
        coeff(2)    = 0.0D+00
        coeff(j2)   = 0.0D+00
        coeff(j2+1) = 0.0D+00

      end if

      ibot = s - j + 1
      itop = s + j - 1

      do ii = ibot, itop
        i = ii - ibot + 1
        call recur ( ii, ci, di, ei )
        coeff(i+2) = coeff(i+2) + ( work(i,ix(2) ) / ci ) * cj1
        coeff(i+1) = coeff(i+1)
     &    + work(i,ix(2)) * ( dj1 - ( cj1 / ci ) * di )
        coeff(i) = coeff(i) - ( work(i,ix(2)) / ci ) * cj1 * ei
      end do

      ii = ix(1)
      ix(1) = ix(2)
      ix(2) = ii

      do i = 1, j2 + 1
        work(i,ix(2)) = coeff(i)
      end do

      return
      end
      subroutine extend ( n, m, m0, t, recur, symmet, start, pnodes,
     &  h0, nexp, idigit, wt, nodes, qr, qi, err, ext, iwork,
     &  worka, lda, workb, ldb, iflag )

c*********************************************************************72
c
cc EXTEND extends a quadrature rule by adding new nodes.
c
c  Discussion:
c
c    This function calculates the N+M node quadrature rule composed of
c    N pre-assigned nodes together with M nodes chosen optimally to
c    achieve algebraic degree of precision of at least N+2*M-1.
c
c    The orthogonal system of polynomials associated with the quadrature
c    weight is defined generally by the recurrence relation specified in the
c    user supplied function RECUR.
c
c  Modified:
c
c    02 March 2009
c
c  Author:
c
c    Original FORTRAN77 version by Thomas Patterson.
c
c  Reference:
c
c    Thomas Patterson,
c    An algorithm for generating interpolatory quadrature rules of the highest
c    degree of precision with preassigned nodes for general weight functions,
c    Transactions on Mathematical Software,
c    Volume 15, Number 2, June 1989, pages 123-136.
c
c    Thomas Patterson,
c    Algorithm 672:
c    EXTEND: generation of interpolatory quadrature rules of the highest degree
c    of precision with preassigned nodes for general weight functions,
c    Transactions on Mathematical Software,
c    Volume 15, Number 2, June 1989, pages 137-143.
c
c  Parameters:
c
c    Input/output, integer N; on input, the number of pre-assigned nodes,
c    and the upper limit to the expansion.  On output, if the computation was
c    successful, N is reset to N+M, which is the appropriate next input value
c    in cases where this routine is called iteratively.
c
c    Input, integer M, the number of nodes to be optimally added.
c
c    Input/output, integer M0; on input, the lower limit to the expansion of T.
c    This is ignored if START is TRUE.  If the computation is successful, the
c    output value of M0 is reset to M, which is the appropriate next input
c    value in cases where this routine is called iteratively.
c
c    Input/output, double precision T(0:max(N-M0,M)); on input, the coefficients
c    TI of the polynomial whose roots define the N pre-assigned nodes of the
c    quadrature rule and expressed as:
c      SUM (I = M0 to N) (TI/HI)*P(I,X)
c    where HI is the integral of W(X)*P(I,X)**2 over the
c    interval for which orthogonality with respect the weight
c    W(X) is defined (moment integrals) and P(I,X) is the
c    orthogonal polynomial of degree I.  Element T(I-M0) holds the
c    value of TI.
c    Note that T is either,
c    (1) provided explicitly,
c    (2) generated automatically from the N pre-assigned nodes
c        given in PNODES(*) (if START is TRUE.)
c    or,
c    (3) generated from a previous call to the function.
c    This array should be declared to have at least
c    max(N-M0+1,M+1) elements in the calling program.
c    The service function TRANSF can be used to transform
c    the expansion to the required input form if desired
c    with the parameter IFLAG set to 1.
c    On output, the coefficients TI of the new orthogonal expansion whose
c    roots are the nodes of the extended quadrature rule (that is, the
c    pre-assigned nodes plus the extended nodes) and expressed as:
c      SUM (I = M to N+M) (TI/HI)*P(I,X)
c    T(I-M) holds the value of TI.
c    This polynomial can be used as input for further extensions.  The service
c    function TRANSF can be used to remove the moment factors from the
c    expansion if desired with the parameter IFLAG set to 0.
c
c    Input, external RECUR ( ), the user supplied function which defines the
c    orthogonal polynomials.  Given K, CALL RECUR(K,C,D,E) gives
c    the coefficients C,D and E such that,
c      P(K+1,X) = (C*X+D)*P(K,X)+E*P(K-1,X)
c    The parameters are defined as follows:
c      K = Index
c      C, D, E = Parameters in the recurrence relation (functions of K)
c
c    Input, logical SYMMET,
c    * FALSE, if no advantage is to be taken of symmetry, if any,
c      about x = 0 in the interval of integration and the
c      orthogonality  weight function. Note that if symmetry in
c      fact does exist setting this parameter to zero will still
c      produce correct results - only efficiency is effected.
c    * TRUE, if the extended rule computations should
c      exploit symmetry about x = 0 in the interval of
c      integration and the orthogonality  weight function.
c      This reduces the size of the system of linear equations
c      determining EXT by a factor of about 2 (see WORKA). If
c      symmetry does not in fact exist erroneous results will be
c      produced.
c
c    Input, logical START,
c    * TRUE, then the polynomial T is generated to have
c      the pre-assigned nodes (PNODES) as its roots.
c    * FALSE. then the supplied values of the coefficients
c      of T are used directly.
c
c    Input/output, double precision PNODES(N+M), on input, the pre-assigned nodes.
c    On output, the nodes of the extended quadrature rule made up from the
c    original pre-assigned nodes and the new optimally extended nodes.  These
c    values can be used in subsequent iterative use.
c
c    Input, double precision H0, the integral of the orthogonality weight
c    function over the interval of integration. Zero moment integral.
c
c    Input, integer NEXP, the largest negative decimal exponent supported on the
c    computer.  (Positive number - typical value 38 for VAX/VMS).
c    Weights less than approximately 10**(-NEXP) are set to zero
c    when the Christoffel-Darboux identity is used (N = M).
c    This may be set to INT(LOG10(X1MACH(2))) where X is set to
c    correspond to the appropriate precision in the PORT library.
c
c    Input, integer IDIGIT, the node convergence parameter, an integer greater
c    than 0.  An attempt is made to calculate the nodes to the maximum
c    accuracy possible by the machine precision available.
c    IDIGIT controls the assessment procedure to take account of
c    round-off errors and specifies the number of least significan
c    decimal digits that can be ignored (i.e. attributed
c    to round-off) in the computed relative error.  A typical
c    value is 5.
c
c    Output, double precision WT(N+M), the quadrature weights for
c    the extended rule associated with the nodes in PNODES.
c
c    Output, integer NODES, the number of extended nodes found.  Normally this
c    equals M, but NODES will be less than M in cases where the computation
c    was terminated because an error condition was encountered.
c
c    Output, double precision QR(M), the real parts of the extended nodes.
c
c    Output, double precision QI(M), the imaginary parts of the extended
c    nodes.
c
c    Output, double precision ERR(M), a measure of the relative error in the
c    nodes.  This may be inspected if the convergence error flag has been raised
c    (IFLAG = 3) to decide if the nodes in question are acceptable.  (ERR(*)
c    actually gives the mean last correction to the quadratic factor in the
c    generalized Bairstow root finder (see BAIR).
c
c    Output, double precision EXT(0:M), the coefficients of the polynomial whose
c    roots are the  extended nodes (QRS(*),QIS(*)) and expressed as:
c      EXT = SUM (I = 0 to M) EXT(I)*P(I,X).
c
c    Output, integer IWORK(max(M,N)), node convergence flags.  Elements 1 to NODES
c    give information on the convergence of the roots of the polynomial EXT
c    corresponding to each extended node.
c    * 0, Convergence of I th root satisfactory;
c    * 1, Convergence of I th root unsatisfactory.
c
c    Workspace, double precision WORKA(max(M+1,N+1),max(M+1,N+1)).
c    If SYMMET = TRUE, the dimensions can be reduced to max(M/2+1,N) by
c    max(M/2+1,N+1).
c
c    Input, integer LDA, the leading dimension of WORKA in the calling program.
c
c    Input, double precision WORKB(2*M+1,3).
c
c    Input, integer LDB, the leading dimension of WORKB.
c
c    Output, integer IFLAG, error flag.
c    * 0, No error detected;
c    * 1, The linear system of equations defining the polynomial
c      whose roots are the extended nodes became singular or
c      very  ill-conditioned.   (FATAL).
c    * 2, The linear system of equations used to generate the
c      polynomial T when START is TRUE became singular
c      or very ill-conditioned. (FATAL).
c    * 3, Poor convergence has been detected in the calculation
c      of the roots of EXT corresponding to the new
c      nodes or all nodes have not been found (M not equal
c      to NODES). See also ERR(*).
c    * 4, Possible imaginary nodes detected.
c    * 5, Value of N and M incompatible for SYMMET = TRUE.
c      Both cannot be odd. (FATAL)
c    * 6, Test of new quadrature rule has failed.
c
      implicit none

      integer lda
      integer ldb
      integer m
      integer n

      double precision err(m)
      double precision ext(0:m)
      double precision h0
      integer i
      integer ideg
      integer idigit
      integer ierr
      integer iflag
      integer iwork(lda)
      integer k
      integer m0
      integer nexp
      integer nlast
      integer nodes
      integer num
      double precision pnodes(n+m)
      double precision qi(m)
      double precision qr(m)
      external recur
      integer s
      logical start
      logical symmet
      double precision t(0:n)
      double precision test
      double precision worka(0:lda-1,0:*)
      double precision workb(0:ldb-1,3)
      double precision wt(n+m)
c
c  Check the size of LDA.
c
      if ( lda < max ( n + 1, m + 1 ) ) then
        iflag = 7
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EXTEND - Fatal error!'
        write ( *, '(a)' ) '  LDA < max ( N + 1, M + 1 ).'
        return
      end if
c
c  Check the size of LDB.
c
      if ( ldb < 2 * m + 1 ) then
        iflag = 8
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EXTEND - Fatal error!'
        write ( *, '(a)' ) '  LDB < 2  M + 1.'
        return
      end if

      iflag = 0
      nodes = 0
      ideg = n + 2 * m - 1
c
c  Look for incompatible values of N and M.
c
      if ( symmet ) then
c
c  Both N and M cannot be odd.
c
        if ( mod ( n, 2 ) .eq. 1 .and. mod ( m, 2 ) .eq. 1 ) then
          iflag = 5
          return
        end if
      end if
c
c  Generate if required the initial T polynomial corresponding to
c  prescribed pre-assigned nodes.
c
      if ( start .and. n .ne. 0 ) then
        call assign ( n, pnodes, iwork, worka, lda, recur, t, ierr )
        m0 = 0
        if ( ierr .ne. 0 ) then
          iflag = 2
          return
        end if
      end if
      nlast = n
c
c  Generate extended expansion coefficients and overwrite T.
c
      call gener ( t, m0, n, m, recur, symmet, ext, iwork, worka, lda,
     &  workb, ldb, ierr )

      if ( ierr .ne. 0 ) then
        iflag = 1
        return
      end if
c
c  Find extended nodes as roots of EXT(*).
c
      call solve ( ext, m, symmet, recur, idigit, qr, qi,
     &  nodes, err, iwork, workb, ldb, ierr )

      if ( ierr .ne. 0 ) then
        iflag = ierr + 2
      end if

      if ( iflag .ne. 0 ) then
        return
      end if
c
c  Accumulate nodes for extended rule.
c
      do i = 1, m
        pnodes(nlast+i) = qr(i)
      end do
c
c  Reorder.
c
      call rsort ( pnodes, n, 1 )
c
c  Compute weights (only for positive nodes if symmetric).
c
      if ( symmet ) then
        num = ( n + 1 ) / 2
      else
        num = n
      end if

      do i = 1, num
        call weight ( t, m0, n, pnodes(i), recur, h0, nexp, wt(i) )
        if ( symmet ) then
          wt(n-i+1) = wt(i)
        end if
      end do
c
c  Test the new rule.
c
      do k = 0, min ( 4, ideg / 2 )
        call check ( n, pnodes, wt, k, h0, recur, test, ierr )
        if ( ierr .eq. 1 ) then
          iflag = 6
          return
        end if
      end do

      return
      end
      subroutine gener ( t, m0, n, m, recur, symmet, ext, iwork,
     &  worka, lda, workb, ldb, iflag )

c*********************************************************************72
c
cc GENER calculates the polynomial defining the optimal new nodes.
c
c  Discussion:
c
c    Given N pre-assigned quadrature nodes defined as the roots of the
c    polynomial expansion:
c
c      SUM (I = M0 to N) (TI/HI)*P(I,X)
c
c    calculate the polynomial expansion:
c
c      SUM (I = 0 to M) SI*P(I,X)
c
c    whose roots are the M optimal nodes and the new expansion:
c
c      SUM (I = M to N+M) (RI/HI)*P(I,X)
c
c    whose roots are to the (N + M) nodes of the full extended quadrature rule.
c
c  Modified:
c
c    02 March 2009
c
c  Author:
c
c    Original FORTRAN77 version by Thomas Patterson.
c
c  Reference:
c
c    Thomas Patterson,
c    An algorithm for generating interpolatory quadrature rules of the highest
c    degree of precision with preassigned nodes for general weight functions,
c    Transactions on Mathematical Software,
c    Volume 15, Number 2, June 1989, pages 123-136.
c
c    Thomas Patterson,
c    Algorithm 672:
c    EXTEND: generation of interpolatory quadrature rules of the highest degree
c    of precision with preassigned nodes for general weight functions,
c    Transactions on Mathematical Software,
c    Volume 15, Number 2, June 1989, pages 137-143.
c
c  Parameters:
c
c    Input/output, double precision T(0:N); on input, the coefficients TI of
c    the polynomial whose roots define the N pre-assigned nodes of the
c    quadrature rule and expressed as:
c      SUM (I = M0 to N) (TI/HI)*P(I,X)
c    where HI is the integral of W(X)*P(I,X)**2 over the
c    interval for which orthogonality with respect the weight
c    W(X) is defined (moment integrals) and P(I,X) is the
c    orthogonal polynomial of degree I. T(I-M0) holds the
c    value of TI. This array should be declared to have at least
c    max(N-M0+1,M+1) elements in the calling program.
c    On output, the coefficients of the new orthogonal
c    expansion whose roots are the nodes of the extended quadrature rule
c    (that is the pre-assigned nodes plus the extended nodes).
c    It is expressed as:
c      SUM (I = M to N+M) (TI/HI)*P(I,X)
c    where N and M have their original values. T(I-M) holds
c    the value of TI. See input argument of T for definitions.
c
c    Input/output, integer M0, the lower limit to the expansion of T.  On
c    output, this is updated to correspond with the output value of T.
c
c    Input/output, integer N, the upper limit to expansion of T.  On output,
c    this is updated to correspond with the output value of T.
c
c    Input, integer M, the number of nodes to be optimally added.
c
c    Input, external RECUR ( ), the user supplied function which defines the
c    orthogonal polynomials.  Given K, CALL RECUR(K,C,D,E) gives
c    the coefficients C,D and E such that,
c      P(K+1,X) = (C*X+D)*P(K,X)+E*P(K-1,X)
c    The parameters are defined as follows:
c      K = Index
c      C, D, E = parameters in the recurrence relation (functions of K).
c
c    Input, SYMMET
c    * FALSE, if no advantage is to be taken of symmetry, if any,
c      about x = 0 in the interval of integration and the
c      orthogonality  weight function. Note that if symmetry in
c      fact does exist setting this parameter to zero will still
c      produce correct results - only efficiency is effected.
c    * TRUE, if the extended rule computations should
c      exploit symmetry about x = 0 in the interval of
c      integration and the orthogonality  weight function.
c      This reduces the size of the system of linear equations
c      determining EXT by a factor of about 2 (see WORKA). If
c      symmetry does not in fact exist erroneous results will be
c      produced.
c
c    Output, double precision EXT(0:M), the coefficients of the polynomial whose
c    roots are the  new extended nodes and expressed as:
c      EXT = SUM (I = 0 to M) EXT(I)*P(I,X)
c
c    Workspace, integer IWORK(M).
c
c    Workspace, double precision WORKA(M+1,max(M+1,N+1)).
c    If SYMMET = TRUE, the dimension can be reduced to
c    M/2+1 by max(M/2+1,N/2+1).
c
c    Input, integer LDA, the number of elements in the leading dimension of
c    WORKA declared in the calling program.
c
c    Input, double precision WORKB(2*M+1,3).
c
c    Input, integer LDB, the number of elements in the leading dimension of
c    WORKB declared in the calling program.
c
c    Output, integer IFLAG, error flag.
c    * 0, No error detected
c    * 1, The linear system of equations defining the polynomial whose roots
c      are the extended nodes became singular or very ill-conditioned.
c
      implicit none

      integer lda
      integer ldb
      integer m
      integer n

      double precision ext(0:m)
      integer i
      integer ibot
      integer ic
      integer ifail
      integer iflag
      integer info
      integer ir
      integer iref
      integer itop
      integer iwork(m)
      integer j
      integer m0
      logical miss
      logical msodd
      integer neq
      logical neven
      integer nm
      double precision pmax
      external recur
      integer s
      logical symmet
      double precision t(0:n)
      double precision total
      double precision worka(0:lda-1,0:*)
      double precision workb(0:ldb-1,*)

      iflag = 0
c
c  Look for trivial case.
c
      if ( n .eq. 0 ) then

        do i = 0, m - 1
          ext(i) = 0.0D+00
        end do

        ext(m) = 1.0D+00
        t(0) = 1.0D+00
        n = m
        m0 = m
        return

      end if
c
c  General case.
c
      neven = mod ( n, 2 ) .eq. 0
      nm = n + m
c
c  Form matrix.
c
      do s = 0, m

        msodd = mod ( m + s, 2 ) .eq. 1

        if ( neven .and. msodd .and. symmet ) then
          go to 60
        end if

        do j = 0, s

          call eprod ( s, j, workb(0,1), workb(0,2), ldb, recur,
     &      ifail )

          if ( mod ( n + s + j, 2 ) .ne. 1 .or. .not. symmet ) then

            iref = s - j
            itop = min ( n, j + s )
            ibot = max ( m0, iref )

            total = 0.0D+00
            do i = ibot, itop
              total = total + t(i-m0) * workb(i-iref,1)
            end do

            if ( .not. symmet ) then

              worka(s,j) = total
              worka(j,s) = total

            else


              if ( neven ) then
                worka(s/2,j/2) = total
                worka(j/2,s/2) = total
              else if ( msodd ) then
                worka(s/2,j/2) = total
              else
                worka(j/2,s/2) = total
              end if

            end if

          end if

        end do

60      continue

      end do

      if ( symmet ) then
        neq = m / 2
      else
        neq = m
      end if
c
c  Solve for expansion coefficients.
c
      if ( 0 .lt. neq ) then

        call dgefa ( worka, lda, neq, iwork, info )

        if ( info .ne. 0 ) then
          iflag = 1
          return
        end if

        call dgesl ( worka, lda, neq, iwork, worka(0,neq), 0 )
c
c  Store expansion coefficients.
c
        do j = 0, neq - 1
          ext(j) = - worka(j,neq)
        end do
        ext(neq) = 1.0D+00

      end if
c
c  Calculate new T polynomial.
c
      if ( .not. symmet ) then
c
c  Non-symmetric case.
c
        do s = m, nm

          if ( s .ne. m ) then

            do j = 0, m

              call eprod  (s, j, workb(0,1), workb(0,2), ldb, recur,
     &          ifail )

              iref = s - j
              itop = min ( n, j + s )
              ibot = max ( m0, iref )

              total = 0.0D+00
              do i = ibot, itop
                ir = i - iref
                total = total + t(i-m0) * workb(i-iref,1)
             end do

              worka(m,j) = total

            end do

          end if

          total = 0.0D+00
          do i = 0, m
            total = total + ext(i) * worka(m,i)
          end do

          worka(m-1,s-m) = total

        end do
c
c  Overwrite old values of T.
c
        do i = 0, n
          t(i) = worka(m-1,i)
        end do
c
c  Symmetric case.
c
      else

        do s = m, nm

          if ( mod ( n + m + s, 2 ) .ne. 1 ) then

            do j = 0, m

              call eprod ( s, j, workb(0,1), workb(0,2), ldb, recur,
     &          ifail )

              if ( mod ( n + s + j, 2 ) .ne. 1 ) then

                iref = s - j
                itop = min ( n, j + s )
                ibot = max ( m0, iref )

                total = 0.0D+00
                do i = ibot, itop
                  ir = i - iref
                  total = total + t(i-m0) * workb(i-iref,1)
                end do

                worka(neq,j/2) = total

             end if

            end do

            total = 0.0D+00
            do i = 0, neq
              total = total + ext(i) * worka(neq,i)
            end do

            if ( 0 .le. neq - 1 ) then
              worka(neq-1,(s-m)/2) = total
            end if

          end if

        end do
c
c  Overwrite old values of T in full unsymmetric form.
c
        ic = n / 2
        miss = .true.

        do j = n, 0, -1
          miss = .not. miss
          if ( miss ) then
            t(j) = 0.0D+00
          else
            if ( 0 .le. neq - 1 ) then
              t(j) = worka(neq-1,ic)
              ic = ic - 1
            end if
          end if
        end do
c
c  Convert EXT to full unsymmetric form.
c
        workb(m,1) = 1.0D+00
        ic = neq - 1
        miss = .false.

        do j = m - 1, 0, -1
          miss = .not. miss
          if ( miss ) then
            workb(j,1) = 0.0D+00
          else
            workb(j,1) = ext(ic)
            ic = ic - 1
          end if
        end do

        do j = 0, m
          ext(j) = workb(j,1)
        end do

      end if
c
c  Scale new T polynomial.
c
      pmax = 0.0D+00
      do i = 0, n
        pmax = max ( pmax, abs ( t(i) ) )
      end do

      do i = 0, n
        t(i) = t(i) / pmax
      end do

      n = nm
      m0 = m

      return
      end
      function idamax ( n, dx, incx )

c*********************************************************************72
c
cc IDAMAX finds the index of element having maximum absolute value.
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
c    Input, double precision X(*), the vector to be examined.
c
c    Input, integer INCX, the increment between successive entries of SX.
c
c    Output, integer IDAMAX, the index of the element of SX of maximum
c    absolute value.
c
      implicit none

      double precision dx(*),dmax
      integer idamax
      integer i,incx,ix,n

      idamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c  code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do  i = 2,n
         if ( dmax .lt. dabs ( dx(ix) ) ) then
           idamax = i
           dmax = dabs(dx(ix))
         end if
         ix = ix + incx
      end do
      return
c
c  code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do i = 2,n
        if( dmax .lt. dabs(dx(i)) ) then
          idamax = i
          dmax = dabs(dx(i))
        end if
      end do

      return
      end
      subroutine lfact ( gamma, delta, n, xnode, recur, r, dr )

c*********************************************************************72
c
cc LFACT removes a linear factor from a polynomial expansion.
c
c  Discussion:
c
c    This function removes the linear factor (X-XNODE) from the polynomial
c    expansion:
c
c      SUM(I = 0 to N) GAMMA(I) * P(I,X)
c
c    to give the quotient:
c
c      SUM (I = 0 to N-1) DELTA(I) * P(I,X).
c
c    and the remainder and its derivative at XNODE.
c
c  Modified:
c
c    02 March 2009
c
c  Author:
c
c    Original FORTRAN77 version by Thomas Patterson.
c
c  Reference:
c
c    Thomas Patterson,
c    An algorithm for generating interpolatory quadrature rules of the highest
c    degree of precision with preassigned nodes for general weight functions,
c    Transactions on Mathematical Software,
c    Volume 15, Number 2, June 1989, pages 123-136.
c
c    Thomas Patterson,
c    Algorithm 672:
c    EXTEND: generation of interpolatory quadrature rules of the highest degree
c    of precision with preassigned nodes for general weight functions,
c    Transactions on Mathematical Software,
c    Volume 15, Number 2, June 1989, pages 137-143.
c
c  Parameters:
c
c    Input, double precision GAMMA(0:N), the polynomial from which the factor
c    is to be removed and expressed as:
c      GAMMA = SUM (I = 0 to N) GAMMA(I)*P(I,X)
c
c    Output, double precision DELTA(0:N-1), the quotient polynomial expressed as:
c      DELTA = SUM (I = 0 to N-1) DELTA(I)*P(I,X)
c
c    Input, integer N, the degree of GAMMA.
c
c    Input, double precision XNODE, the node to be removed.
c
c    Input, external RECUR ( ), the function which defines the orthogonal
c    polynomials.  See EXTEND for a full description.
c
c    Output, double precision R, the remainder from division.
c
c    Output, double precision DR, the derivative of R with respect to XNODE.
c    (-R/DR is the Newton correction).
c
      implicit none

      integer n

      double precision bk1
      double precision bk2
      double precision ck
      double precision ck1
      double precision ckm1
      double precision dbk1
      double precision dbk2
      double precision delta(0:n-1)
      double precision dk
      double precision dk1
      double precision dkm1
      double precision dr
      double precision ek
      double precision ek1
      double precision ekm1
      double precision gamma(0:n)
      integer k
      double precision r
      external recur
      double precision xnode

      bk1 = 0.0D+00
      bk2 = 0.0D+00
      dbk1 = 0.0D+00
      dbk2 = 0.0D+00
      call recur ( n, ck, dk, ek )
      call recur ( n + 1, ck1, dk1, ek1 )

      do k = n, 0, -1

        r = gamma(k) + ( dk + xnode * ck ) * bk1 + ek1 * bk2
        dr = ( dk + xnode * ck ) * dbk1 + ek1 * dbk2 + ck * bk1
        bk2 = bk1
        bk1 = r
        dbk2 = dbk1
        dbk1 = dr

        if ( k .ne. 0 ) then
          call recur ( k - 1, ckm1, dkm1, ekm1 )
          delta(k-1) = r * ckm1
        end if

        ek1 = ek
        ck = ckm1
        dk = dkm1
        ek = ekm1

      end do

      return
      end
      subroutine newton ( t, n, xnode, recur, idigit, delta, errval,
     &  ifail )

c*********************************************************************72
c
cc NEWTON applies Newton's method for a single root of a polynomial.
c
c  Discussion:
c
c    This function applies Newton's method to find a single root of the
c    polynomial T:
c
c      T = SUM (I = 0 to N) T(I) * P(I,X)
c
c    where P(I,X) are the orthogonal polymonials whose recurrence
c    relation is defined by RECUR.
c
c    The value of T is found from the remainder when T is divided
c    by (X-XNODE).  The derivative (of the remainder) is
c    calculated simultaneously.
c
c    The deflated polynomial:
c
c      DELTA = SUM (I = 0 to N-1) DELTA(I) * P(I,X)
c
c    is also computed.
c
c  Modified:
c
c    02 March 2009
c
c  Author:
c
c    Original FORTRAN77 version by Thomas Patterson.
c
c  Reference:
c
c    Thomas Patterson,
c    An algorithm for generating interpolatory quadrature rules of the highest
c    degree of precision with preassigned nodes for general weight functions,
c    Transactions on Mathematical Software,
c    Volume 15, Number 2, June 1989, pages 123-136.
c
c    Thomas Patterson,
c    Algorithm 672:
c    EXTEND: generation of interpolatory quadrature rules of the highest degree
c    of precision with preassigned nodes for general weight functions,
c    Transactions on Mathematical Software,
c    Volume 15, Number 2, June 1989, pages 137-143.
c
c  Parameters:
c
c    Input, double precision T(0:N), the polynomial whose roots define the nodes
c    of the quadrature rule and expressed as:
c      T = SUM (I = 0 to N) T(I)*P(I,X)
c
c    Input, integer N, the degree of the expansion of T.
c
c    Input/output, double precision XNODE, on input, a rough estimate for
c    the root.  On output, this estimate has been improved using Newton's method.
c
c    Input, external RECUR ( ), the function which defines the orthogonal
c    polynomials.  See EXTEND for a full description.
c
c    Input, integer IDIGIT, the node convergence parameter, an integer greater
c    than 0.  An attempt is made to calculate the nodes to the maximum
c    accuracy possible by the machine precision available.
c    IDIGIT controls the assessment procedure to take account of
c    round-off errors and specifies the number of least significan
c    decimal digits that can be ignored, that is, attributed
c    to round-off, in the computed relative error.  A typical value is 5.
c
c    Output, double precision DELTA(0:N-1), the coefficients of the deflated
c    polynomial.
c
c    Output, double precision ERRVAL, the value of the correction.  May be used
c    as a measure of the root accuracy when convergence is not achieved.
c
c    Output, integer IFAIL, error flag.
c    * 0, Convergence OK.
c    * 1, Unsatisfactory convergence after 50 iterations.
c
      implicit none

      integer n

      double precision delta(0:n-1)
      double precision dr
      double precision eps
      double precision errval
      integer idigit
      integer ifail
      integer iter
      double precision r
      external recur
      double precision t(0:n)
      double precision tol
      double precision xnode

      iter = 50
      tol = 10.0D+00**( - max ( 1, idigit ) )

10    continue

      iter = iter - 1

      if ( iter .lt. 0 ) then
        ifail = 1
        errval = abs ( eps )
        return
      end if

      call lfact ( t, delta, n, xnode, recur, r, dr )
      eps = - r / dr
      xnode = xnode + eps

      if ( abs ( xnode ) + tol * abs ( eps ) .ne. abs ( xnode ) ) then
        go to 10
      end if
c
c  Final iteration.
c
      call lfact ( t, delta, n, xnode, recur, r, dr )
      eps = - r / dr
      xnode = xnode + eps
      ifail = 0
      errval = abs ( eps )

      return
      end
      subroutine qfact ( n, gamma, delta, recur, alpha, beta, a, b,
     &  aa, ab, ba, bb )

c*********************************************************************72
c
cc QFACT divides a polynomial by a quadratic factor.
c
c  Discussion:
c
c    This functioin divides the polynomial:
c
c      SUM(I = 0 to N) GAMMA(I) * P(I,X)
c
c    by the quadratic factor:
c
c      P(2,X) - ALPHA * P(1,X) - BETA * P(0,X)
c
c    giving the quotient:
c
c      SUM(I = 0 to N-2) DELTA(I) * P(I,X)
c
c    and remainder:
c
c      A * P(1,X) + B * P(0,X)
c
c    where P(I,X) is the orthogonal polynomial of degree I defined by RECUR.
c
c  Modified:
c
c    02 March 2009
c
c  Author:
c
c    Original FORTRAN77 version by Thomas Patterson.
c
c  Reference:
c
c    Thomas Patterson,
c    An algorithm for generating interpolatory quadrature rules of the highest
c    degree of precision with preassigned nodes for general weight functions,
c    Transactions on Mathematical Software,
c    Volume 15, Number 2, June 1989, pages 123-136.
c
c    Thomas Patterson,
c    Algorithm 672:
c    EXTEND: generation of interpolatory quadrature rules of the highest degree
c    of precision with preassigned nodes for general weight functions,
c    Transactions on Mathematical Software,
c    Volume 15, Number 2, June 1989, pages 137-143.
c
c  Parameters:
c
c    Input, integer N, the degree of GAMMA.
c
c    Input, double precision GAMMA(0:N), a polynomial to be divided by
c    a quadratic factor.
c
c    Output, double precision DELTA(0:N-2), the quotient polynomial of
c    degree N-2.
c
c    Input, external RECUR(), the function which defines the
c    orthogonal polynomials.  See EXTEND for a full description.
c
c    Input, double precision ALPHA, BETA, the coefficients of the quadratic
c    factor.
c
c    Output, double precision A, B, the remainder coefficients.
c
c    Output, double precision AA, the partial of A with respect to ALPHA.
c
c    Output, double precision AB, the partial of A with respect to BETA.
c
c    Output, double precision BA, the partial of B with respect to ALPHA.
c
c    Output, double precision BB, the partial of B with respect to BETA.
c
      implicit none

      integer n

      double precision a
      double precision aa
      double precision ab
      double precision adn
      double precision adnm1
      double precision adnp1
      double precision adnp2
      double precision alpha
      double precision b
      double precision ba
      double precision bb
      double precision bdn
      double precision bdnm1
      double precision bdnp1
      double precision bdnp2
      double precision beta
      double precision c0
      double precision c1
      double precision c2
      double precision c3
      double precision cf1
      double precision cf2
      double precision cf3
      double precision cf4
      double precision cf5
      double precision ck
      double precision ck1
      double precision ck2
      double precision ck3
      double precision ck4
      double precision d
      double precision d0
      double precision d1
      double precision d2
      double precision d3
      double precision da
      double precision db
      double precision delta(0:n-2)
      double precision dk
      double precision dk1
      double precision dk2
      double precision dk3
      double precision dk4
      double precision dn
      double precision dnm1
      double precision dnp1
      double precision dnp2
      double precision e0
      double precision e1
      double precision e2
      double precision e3
      double precision ek
      double precision ek1
      double precision ek2
      double precision ek3
      double precision ek4
      double precision gamma(0:n)
      integer k
      double precision r0
      double precision r1
      double precision r2
      external recur
      double precision rk1
      double precision rk2
      double precision rk3
      double precision rs
      double precision rs0
      double precision rs1
      double precision sk2
      double precision sn1
      double precision sn2
      double precision sn3
      double precision sn4
      double precision tk2
      double precision uk2
      double precision v1
      double precision vk2
      double precision vlk1
      double precision vlk2
      double precision vlk3
      double precision vlk4
      double precision vm1
      double precision vm2
      double precision vmk1
      double precision vmk2
      double precision vmk3
      double precision w0
      double precision w1
      double precision wk2
c
c  Initialize coefficients.
c
      dnp2 = 0.0D+00
      dnp1 = 0.0D+00
      dn   = 0.0D+00
      dnm1 = 0.0D+00
c
c  Partial coefficients wrt ALPHA.
c
      adnp2 = 0.0D+00
      adnp1 = 0.0D+00
      adn   = 0.0D+00
      adnm1 = 0.0D+00
c
c  Partial coefficients wrt BETA.
c
      bdnp2 = 0.0D+00
      bdnp1 = 0.0D+00
      bdn   = 0.0D+00
      bdnm1 = 0.0D+00
c
c  Scaling parameters.
c
      sn1 = 1.0D+00
      sn2 = 1.0D+00
      sn3 = 1.0D+00
      sn4 = 1.0D+00

      call recur ( 0, c0, d0, e0 )
      call recur ( 1, c1, d1, e1 )
      call recur ( 2, c2, d2, e2 )
      call recur ( 3, c3, d3, e3 )

      r0 = - c0 * e1 / c1
      r1 = - c0 * e2 / c2
      r2 = - c0 * e3 / c3

      vm1 = d0 - c0 * d1 / c1
      vm2 = d0 - c0 * d2 / c2

      w0 = - r1 * e1
      w1 = - c1 * r2 * e2 / c2
      v1 = d1 * r1 - c1 * vm2 * e2 / c2 - c1 * r1 * d1 / c1
      k = n - 2

      call recur ( k + 4, ck4, dk4, ek4 )
      call recur ( k + 3, ck3, dk3, ek3 )
      call recur ( k + 2, ck2, dk2, ek2 )
      call recur ( k + 1, ck1, dk1, ek1 )

      vlk4 = c0 / ck3
      vlk3 = c0 / ck2
      vlk2 = c0 / ck1
      rk3 = - c0 * ek4 / ck4
      rk2 = - c0 * ek3 / ck3
      vmk3 = d0 - dk3 * vlk4
      vmk2 = d0 - dk2 * vlk3
c
c  Extract quadratic factor and find partial derivatives
c
      do k = n - 2, 0, - 1

        call recur ( k, ck, dk, ek )
        vlk1 = c0 / ck
        rk1 = - c0 * ek2 / ck2
        vmk1 = d0 - dk1 * vlk2
        sk2 = c1 * vlk1 * vlk2 / c0
        tk2 = vlk2 * ( d1 - c1 * dk2 / ck2 ) + c1 * vmk1 / ck1
        uk2 = d1 * vmk2 + e1 - c1 * vlk3 * ek3 / ck3
     &    - c1 * vmk2 * dk2 / ck2 + c1 * rk1 / ck1
        vk2 = d1 * rk2 - c1 * vmk3 * ek3 / ck3 - c1 * rk2 * dk2 / ck2
        wk2 = - c1 * rk3 * ek3 / ck3
        cf1 = ( alpha * vlk2 - tk2 ) / sn1
        cf2 = ( beta + alpha*  vmk2 - uk2 ) / sn2
        cf3 = ( alpha * rk2 - vk2 ) / sn3
        cf4 = - wk2 / sn4
        rs = gamma(k+2)
        d = rs + cf1 * dnm1 + cf2 * dn + cf3 * dnp1 + cf4 * dnp2
        delta(k) = d / sk2
        da = vlk2 * dnm1 / sn1 + vmk2 * dn / sn2 + rk2 * dnp1 / sn3
     &    + cf1 * adnm1 + cf2 * adn + cf3 * adnp1 + cf4 * adnp2
        db = dn / sn2 + cf1 * bdnm1 + cf2 * bdn + cf3 * bdnp1
     &    + cf4 * bdnp2
c
c  Recycle old values.
c
        sn4 = sn3
        sn3 = sn2
        sn2 = sn1
        sn1 = sk2
        dnp2 = dnp1
        dnp1 = dn
        dn = dnm1
        dnm1 = d
        adnp2 = adnp1
        adnp1 = adn
        adn = adnm1
        adnm1 = da
        bdnp2 = bdnp1
        bdnp1 = bdn
        bdn = bdnm1
        bdnm1 = db
        ck4 = ck3
        ck3 = ck2
        ck2 = ck1
        ck1 = ck
        dk4 = dk3
        dk3 = dk2
        dk2 = dk1
        dk1 = dk
        ek4 = ek3
        ek3 = ek2
        ek2 = ek1
        ek1 = ek
        vlk4 = vlk3
        vlk3 = vlk2
        vlk2 = vlk1
        rk3 = rk2
        rk2 = rk1
        vmk3 = vmk2
        vmk2 = vmk1

      end do

      cf1 = alpha
      cf2 = beta + alpha * vm1 - r1
      cf3 = alpha * r1 - v1
      cf4 = - w1
      cf5 = alpha * r0
      rs0 = gamma(0)
      rs1 = gamma(1)
      dnm1 = dnm1 / sn1
      dn = dn / sn2
      dnp1 = dnp1 / sn3
      dnp2 = dnp2 / sn4
      adnm1 = adnm1 / sn1
      adn = adn / sn2
      adnp1 = adnp1 / sn3
      adnp2 = adnp2 / sn4
      bdnm1 = bdnm1 / sn1
      bdn = bdn / sn2
      bdnp1 = bdnp1 / sn3
      bdnp2 = bdnp2 / sn4
c
c  Remainder.
c
      a = rs1 + cf1 * dnm1 + cf2 * dn + cf3 * dnp1 + cf4 * dnp2
      b = rs0 + beta * dnm1 + cf5 * dn - w0 * dnp1
c
c  Partials.
c
      aa = dnm1 + vm1 * dn + r1 * dnp1 + cf1 * adnm1 + cf2 * adn
     &  + cf3 * adnp1 + cf4 * adnp2
      ab = dn + cf1 * bdnm1 + cf2 * bdn + cf3 * bdnp1 + cf4 * bdnp2
      ba = r0 * dn + beta * adnm1 + cf5 * adn - w0 * adnp1
      bb = dnm1 + beta * bdnm1 + cf5 * bdn - w0 * bdnp1

      return
      end
      subroutine roots ( a0, a1, a2, zreal1, zimag1, zreal2, zimag2,
     &  recur, info )

c*********************************************************************72
c
cc ROOTS calculates roots of a quadratic factor.
c
c  Discussion:
c
c    This function calculates the roots corresponding to the quadratic factor
c
c      A2 * P(2,X) - A1 * P(1,X) - A0 * P(0,X)
c
c    where P(I,X) is a general orthogonal polynomial of degree I
c    defined by the recurrence calculated by RECUR.
c
c  Modified:
c
c    02 March 2009
c
c  Author:
c
c    Original FORTRAN77 version by Thomas Patterson.
c
c  Reference:
c
c    Thomas Patterson,
c    An algorithm for generating interpolatory quadrature rules of the highest
c    degree of precision with preassigned nodes for general weight functions,
c    Transactions on Mathematical Software,
c    Volume 15, Number 2, June 1989, pages 123-136.
c
c    Thomas Patterson,
c    Algorithm 672:
c    EXTEND: generation of interpolatory quadrature rules of the highest degree
c    of precision with preassigned nodes for general weight functions,
c    Transactions on Mathematical Software,
c    Volume 15, Number 2, June 1989, pages 137-143.
c
c  Parameters:
c
c    Input, double precision A0, A1, A2, the coefficients of the quadratic
c    factor.
c
c    Output, double precision ZREAL1, ZIMAG1, the real and imaginary parts
c    of root 1.
c
c    Output, double precision ZREAL2, ZIMAG2, the real and Imaginary parts
c    of root 2.
c
c    Input, external RECUR ( ), the function which defines the orthogonal
c    polynomials.  See EXTEND for full description.
c
c    Output, integer INFO, error flag.
c    * 0, two roots found.
c    * 1, one root only (A2 = 0).
c
      implicit none

      double precision a0
      double precision a1
      double precision a2
      double precision aa
      double precision bb
      double precision c0
      double precision c1
      double precision cc
      double precision d0
      double precision d1
      double precision e0
      double precision e1
      integer info
      integer k
      external recur
      double precision z
      double precision zimag1
      double precision zimag2
      double precision zr
      double precision zreal1
      double precision zreal2

      info = 0

      k = 0
      call recur ( k, c0, d0, e0 )

      if ( a2 .eq. 0.0D+00 ) then
        zreal1 = - ( a0 + a1 * d0 ) / a1 / c0
        zreal2 = 0.0D+00
        zimag1 = 0.0D+00
        zimag2 = 0.0D+00
        info = 1
        return
      end if

      k = 1
      call recur ( k, c1, d1, e1 )

      aa = - c0 * c1 * a2
      bb = - a2 * ( c0 * d1 + d0 * c1 ) + c0 * a1
      cc = - d0 * d1 * a2 - e1 * a2 + a0 + a1 * d0
      z = bb * bb - 4.0D+00 * aa * cc
      zr = sqrt ( abs ( z ) )

      if ( 0.0D+00 .le. z ) then
        zimag1 = 0.0D+00
        zimag2 = 0.0D+00
        zreal1 = 0.5D+00 * ( - bb - sign ( zr, bb ) ) / aa
        zreal2 = cc / aa / zreal1
      else
        zreal1 = - 0.5D+00 * bb / aa
        zreal2 = zreal1
        zimag1 = 0.5D+00 * zr / aa
        zimag2 = - zimag1
      end if

      return
      end
      subroutine rsort ( a, n, iflag )

c*********************************************************************72
c
cc RSORT carries out a simple ripple sort.
c
c  Modified:
c
c    16 February 2009
c
c  Author:
c
c    Original FORTRAN77 version by Thomas Patterson.
c
c  Reference:
c
c    Thomas Patterson,
c    An algorithm for generating interpolatory quadrature rules of the highest
c    degree of precision with preassigned nodes for general weight functions,
c    Transactions on Mathematical Software,
c    Volume 15, Number 2, June 1989, pages 123-136.
c
c    Thomas Patterson,
c    Algorithm 672:
c    EXTEND: generation of interpolatory quadrature rules of the highest degree
c    of precision with preassigned nodes for general weight functions,
c    Transactions on Mathematical Software,
c    Volume 15, Number 2, June 1989, pages 137-143.
c
c  Parameters:
c
c    Input/output, double precision A(N), the array to be sorted.
c
c    Input, integer N, the number of elements to be sorted
c
c    Input, integer IFLAG, determines the sort direction.
c    * 0, for ascending sort;
c    * 1, for descending sort.
c
      implicit none

      integer n

      double precision a(n)
      logical ascend
      logical done
      integer iflag
      integer j
      integer k
      integer k1
      integer k2
      double precision val

      ascend = iflag .eq. 0
c
c  Begin scans.
c
      do j = n - 1, 1, - 1

        done = .true.

        do k = 1, j

          if ( ascend ) then
            k1 = k
            k2 = k + 1
          else
            k1 = k + 1
            k2 = k
          end if
c
c  Exchange elements.
c
          if ( a(k2) .lt. a(k1) ) then
            val = a(k1)
            a(k1) = a(k2)
            a(k2) = val
            done = .false.
          end if

        end do

        if ( done ) then
          return
        end if

      end do

      return
      end
      subroutine solve ( ext, m, symmet, recur, idigit, qr, qi,
     &  nodes, err, icheck, work, ldw, ierr )

c*********************************************************************72
c
cc SOLVE calculates roots of an orthogonal polynomial expansion.
c
c  Discussion:
c
c    This function calculates the roots of the orthogonal polynomial expansion:
c
c      SUM (I = 0 to M) EXT(I)*P(I,X)
c
c  Modified:
c
c    16 February 2009
c
c  Author:
c
c    Original FORTRAN77 version by Thomas Patterson.
c
c  Reference:
c
c    Thomas Patterson,
c    An algorithm for generating interpolatory quadrature rules of the highest
c    degree of precision with preassigned nodes for general weight functions,
c    Transactions on Mathematical Software,
c    Volume 15, Number 2, June 1989, pages 123-136.
c
c    Thomas Patterson,
c    Algorithm 672:
c    EXTEND: generation of interpolatory quadrature rules of the highest degree
c    of precision with preassigned nodes for general weight functions,
c    Transactions on Mathematical Software,
c    Volume 15, Number 2, June 1989, pages 137-143.
c
c  Parameters:
c
c    Input, double precision EXT(0:M), the coefficients of the polynomial whose
c    roots are required (nodes of the quadrature rule) and expressed as:
c      SUM (I = 0 to M) EXT(I)*P(I,X)
c    The recurrence relation for the orthogonal polynomials
c    P(I,X) is defined by RECUR.
c
c    Input, integer M, the upper limit to expansion EXT (polynomial degree).
c
c    Input, logical SYMMET
c    * FALSE, if no advantage can be taken of symmetry
c      about x = 0 in the interval of integration and the
c      orthogonality weight function.
c    * TRUE, if symmetry exists about x = 0 in the interval of
c      integration and the orthogonality weight function.
c
c    Input, external RECUR ( ), the user supplied function which defines the
c    orthogonal polynomials.  Given K, RECUR ( K, C, D, E ) gives
c    the coefficients C,D and E such that,
c      P(K+1,X) = (C*X+D)*P(K,X)+E*P(K-1,X)
c    The parameters are defined as follows:
c      K = Index
c      C, D, E = parameters in the recurrence relation (functions of K).
c
c    Input, integer IDIGIT, the node convergence parameter, an integer greater
c    than 0.  An attempt is made to calculate the nodes to the maximum
c    accuracy possible by the machine precision available.
c    IDIGIT controls the assessment procedure to take account of
c    round-off errors and specifies the number of least significan
c    decimal digits that can be ignored (i.e. attributed
c    to round-off) in the computed relative error.  A typical
c    value is 5.
c
c    Output, double precision QR(M), the real parts of the roots of EXT.
c
c    Output, double precision QI(M), the imaginary parts of the roots
c    of EXT.  (Hopefully these values are zero!).
c
c    Output, integer NODES, the number of extended nodes found.  Normally equals
c    M but see IERR.
c
c    Output, double precision ERR(M), a measure of the relative error in the
c    roots.  This may be inspected if the convergence error flag has been
c    raised (IERR = 2) to decide if the roots in question are acceptable.
c    (ERR(*) actually gives the mean last correction to the quadratic factor
c    in the generalized Bairstow root finder (see BAIR).
c
c    Output, integer ICHECK(M), root convergence flags.  Elements 1 to NODES
c    give information on the convergence of the roots of the polynomial EXT.
c    * 0, Convergence of I th root satisfactory;
c    * 1, Convergence of I th root unsatisfactory.
c
c    Worskpace, double precision WORK(M+1,2).
c
c    Input, integer LDW, the leading dimension of WORK (which must be at
c    least M+1).
c
c    Output, integer IERR
c    * 0, No error detected
c    * 1, Possible imaginary nodes detected.
c    * 2, Poor convergence has been detected in the calculation
c      of the roots of EXT or all roots have not
c      been found (M not equal to NODES).  See also ERR.
c
      implicit none

      integer ldw
      integer m

      double precision a0
      double precision a1
      double precision a2
      double precision c0
      double precision c1
      double precision d0
      double precision d1
      double precision det
      double precision e0
      double precision e1
      double precision err(m)
      double precision errval
      double precision ext(0:m)
      integer i
      integer icheck(m)
      integer idigit
      integer ierr
      integer ifail
      integer info
      integer ip1
      integer ip2
      integer j
      integer nodes
      integer nroot
      double precision p1a
      double precision p1b
      double precision p2a
      double precision p2b
      double precision pmax
      double precision qi(m)
      double precision qr(m)
      external recur
      logical reset
      double precision rt1
      double precision rt2
      double precision sa0
      double precision sa1
      double precision sfa0
      double precision sfa1
      logical symmet
      double precision vrt1
      parameter ( vrt1 = 0.0000001D+00 )
      double precision work(0:ldw-1,2)
      double precision zi1
      double precision zi2
      double precision zr1
      double precision zr2

      nodes = 0
      ierr = 0
c
c  If M is odd, find and remove initial real root using NEWTON iteration.
c  Set WORK(*,1) to polynomial to be processed.
c
      if ( mod ( m, 2 ) .eq. 1 ) then

        zr1 = vrt1
        call newton ( ext, m, zr1, recur, idigit, work(0,1), errval,
     &    ifail )
        nodes = nodes + 1
        icheck(nodes) = ifail
        err(nodes) = errval
        qr(nodes) = zr1
        qi(nodes) = 0.0D+00
        nroot = m - 1

      else

        do i = 0, m
          work(i,1) = ext(i)
        end do
        nroot = m

      end if
c
c  Find the remaining root pairs.
c  Calculate seed approximation for quadratic factor.
c
      if ( nroot .ne. 0 ) then

        call recur ( 0, c0, d0, e0 )
        call recur ( 1, c1, d1, e1 )
        rt1 = vrt1

        if ( symmet ) then
          rt2 = - rt1
        else
          rt2 = 0.0D+00
        end if

        p1a = c0 * rt1 + d0
        p1b = c0 * rt2 + d0
        p2a = ( c1 * rt1 + d1 ) * p1a + e1
        p2b = ( c1 * rt2 + d1 ) * p1b + e1
        det = c0 * ( rt1 - rt2 )
        sa1 = ( p2a - p2b ) / det
        sa0 = ( p1a * p2b - p1b * p2a ) / det
        reset = .true.
c
c  Alternate approximation which introduces small complex component.
c
        rt1 = vrt1
        rt2 = vrt1
        sfa1 = ( c0 * d1 + d0 * c1 ) / c0 + 2.0D+00 * c1 * rt1
        sfa0 = d0 * d1 + e1 - d0 * sfa1
     &    - c0 * c1 * ( rt1 * rt1 + rt2 * rt2 )
c
c  IP1 points to the current deflated polynomial.
c
        ip1 = 1
        ip2 = 2

20      continue

          if ( reset ) then
            a2 = 1.0D+00
            a1 = sa1
            a0 = sa0
            reset = .false.
          end if

          call bair ( nroot, work(0,ip1), work(0,ip2), a0, a1, a2, 
     &      recur, idigit, errval, ifail )
c
c  On failure, try again with complex components introduced.
c
          if ( ifail .ne. 0 ) then
            a2 = 1.0D+00
            a1 = sfa1
            a0 = sfa0
            reset = .true.
            call bair ( nroot, work(0,ip1), work(0,ip2), a0, a1, a2,
     &        recur, idigit, errval, ifail )
          end if
c
c  Apply Bairstow to full expansion to avoid error accumulation.
c
          call bair ( m, ext, work(0,ip2), a0, a1, a2, recur, idigit,
     &      errval, ifail )
c
c  Tidy up the quotient polynomial.
c
          call qfact ( nroot, work(0,ip1), work(0,ip2), recur, a1, a0,
     &      zr1, zr1, zr1, zr1, zr1, zr1 )

          call roots ( a0, a1, a2, zr1, zi1, zr2, zi2, recur, info )
c
c  Record node information.
c
          nodes = nodes + 1
          icheck(nodes) = ifail
          err(nodes) = errval
          qr(nodes) = zr1
          qi(nodes) = zi1

          nodes = nodes + 1
          icheck(nodes) = ifail
          err(nodes) = errval
          qr(nodes) = zr2
          qi(nodes) = zi2

          nroot = nroot - 2
c
c  Make the deflated polynomial current.
c
          i = ip1
          ip1 = ip2
          ip2 = i
c
c  Scale the deflated polynomial.
c
          if ( nroot .le. 0 ) then
            go to 30
          end if

          pmax = 0.0D+00
          do i = 0, nroot
            pmax = max ( pmax, abs ( work(i,ip1) ) )
          end do

          do i = 0, nroot
            work(i,ip1) = work(i,ip1) / pmax
          end do

          go to 20

30      continue

      end if
c
c  Check for poor convergence.
c
      i = 0
      do j = 1, nodes
        i = i + icheck(j)
      end do

      if ( i .ne. 0 ) then
        ierr = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SOLVE - Warning:'
        write ( *, '(a)' ) '  Poor convergence for some roots.'
        return
      end if

      if ( nodes .ne. m ) then
        ierr = 1
        return
      end if
c
c  Look for possible imaginary nodes.
c
      do j = 1, nodes
        if ( qi(j) .ne. 0.0D+00 ) then
          ierr = 2
          return
        end if
      end do

      return
      end
      subroutine transf ( t, m, n, recur, iflag )

c*********************************************************************72
c
cc TRANSF scales a polynomial expansion with respect to the moments.
c
c  Discussion:
c
c    This function scales the polynomial expansion:
c
c      SUM (M to N) TI*P(I,X)
c
c    with respect to the moments HI of the orthogonality weight function
c    giving the expansion:
c
c      H0 * SUM (M to N) (TI/HI)*P(I,X)
c
c    or
c
c      (1/H0) * SUM (M to N) (TI*HI)*P(I,X)
c
c    depending on the value of IFLAG.
c
c  Modified:
c
c    27 February 2009
c
c  Author:
c
c    Original FORTRAN77 version by Thomas Patterson.
c
c  Reference:
c
c    Thomas Patterson,
c    An algorithm for generating interpolatory quadrature rules of the highest
c    degree of precision with preassigned nodes for general weight functions,
c    Transactions on Mathematical Software,
c    Volume 15, Number 2, June 1989, pages 123-136.
c
c    Thomas Patterson,
c    Algorithm 672:
c    EXTEND: generation of interpolatory quadrature rules of the highest degree
c    of precision with preassigned nodes for general weight functions,
c    Transactions on Mathematical Software,
c    Volume 15, Number 2, June 1989, pages 137-143.
c
c  Parameters:
c
c    Input/output, double precision T(0:N), the coefficients TI of the
c    polynomial expansion to be scaled and expressed as:
c      SUM (I = M to N) TI*P(I,X)
c    T(I-M) holds the value of TI.  On output, the polynomial has been scaled.
c
c    Input, integer M, the lower limit to the expansion of T.
c
c    Input, integer N, the upper limit to the expansion of T.
c
c    Input, external RECUR ( ), the function which defines the orthogonal
c    polynomials.  See EXTEND for a full description.
c
c    Input, integer IFLAG, indicates the operation to be carried out.
c    * 0, if coefficient TI is to be replaced by TI*(H0/HI).
c    * 1, if coefficient TI is to be replaced by TI*(HI/H0).
c
      implicit none

      double precision ck
      double precision ckm1
      double precision dk
      double precision ek
      double precision h
      integer iflag
      integer k
      integer m
      integer n
      external recur
      double precision t(0:n)

      h = 1.0D+00

      do k = 0, n

        call recur ( k, ck, dk, ek )

        if ( k .ne. 0 ) then
          h = - ckm1 / ck * ek * h
        end if

        if ( m .le. k ) then
          if ( iflag .eq. 0 ) then
            t(k-m) = t(k-m) / h
          else
            t(k-m) = t(k-m) * h
          end if
        end if

        ckm1 = ck

      end do

      return
      end
      subroutine weight ( t, m, n, xnode, recur, h0, nexp, wt )

c*********************************************************************72
c
cc WEIGHT calculates quadrature weights.
c
c  Discussion:
c
c    This function calculates the quadrature weight associated with the node
c    XNODE in the rule whose nodes are defined by the roots of polynomial T.
c
c    The weight is calculated by dividing T by (X-XNODE) to give:
c
c      S(X) = T(X)/(X-XNODE) = SUM (0 to N-1) G(I)*P(I,X).
c
c    S(X) is then divided by (X-XNODE) to give the remainder R.
c
c    The weight is finally given by H0*G(0)/R. If N = M the
c    Christoffel-Darboux identity result is used to reduce extreme
c    cancellation effects at high degree.
c
c  Modified:
c
c    27 February 2009
c
c  Author:
c
c    Original FORTRAN77 version by Thomas Patterson.
c
c  Reference:
c
c    Thomas Patterson,
c    An algorithm for generating interpolatory quadrature rules of the highest
c    degree of precision with preassigned nodes for general weight functions,
c    Transactions on Mathematical Software,
c    Volume 15, Number 2, June 1989, pages 123-136.
c
c    Thomas Patterson,
c    Algorithm 672:
c    EXTEND: generation of interpolatory quadrature rules of the highest degree
c    of precision with preassigned nodes for general weight functions,
c    Transactions on Mathematical Software,
c    Volume 15, Number 2, June 1989, pages 137-143.
c
c  Parameters:
c
c    Input, double precision T(0:N), the coefficients TI of the polynomial
c    whose roots define the N pre-assigned nodes of the quadrature
c    rule and expressed as:
c      SUM (I = M to N) (TI/HI)*P(I,X)
c    where HI is the integral of W(X)*P(I,X)**2 over the
c    interval for which orthogonality with respect the weight
c    W(X) is defined (moment integrals) and P(I,X) is the
c    orthogonal polynomial of degree I. T(I-M) holds the
c    value of TI. This array should be declared to have at least
c    N-M+1 elements in the calling program.
c
c    Input, integer M, the lower limit to the expansion of T.
c
c    input, integer N, the upper limit to the expansion of T.
c
c    Input, double precision XNODE, the node whose weight is required
c
c    Input, external RECUR ( ), the function which defines the orthogonal
c    polynomials.  See EXTEND for a full description.
c
c    Input, double precision H0, th eintegral of the orthogonality weight
c    function over the interval of integration. Zero moment integral.  Note that
c    P(0,X) is arbitrarily taken to be 1.0
c
c    Input, integer NEXP, the largest negative decimal exponent supported on the
c    computer. (Positive number - typical value 38).
c    Weights less than approximately 10**(-NEXP) are set to zero
c    when the Christoffel-Darboux identity is used (N = M).
c
c    Output, double precision WT, the weight associated with XNODE.
c
      implicit none

      double precision bb
      double precision bk1
      double precision bk2
      double precision ck
      double precision ck1
      double precision ckm1
      double precision d0
      double precision dd
      double precision dk
      double precision dk1
      double precision dk2
      double precision dkm1
      double precision e0
      double precision ek
      double precision ek1
      double precision ekm1
      double precision h
      double precision h0
      integer iscale
      integer itest
      integer j
      integer k
      integer m
      integer n
      integer nexp
      external recur
      double precision rk1
      double precision rk2
      double precision rs
      double precision scale
      double precision t(0:n)
      double precision wt
      double precision xnode
c
c  Check for special case.
c
c  Use Christoffel-Darboux result.
c
      if ( m .eq. n ) then

        bk1 = 0.0D+00
        bk2 = 1.0D+00
        dk1 = 0.0D+00
        dk2 = 0.0D+00
        iscale = 0
        k = 0
        call recur ( k, h, d0, e0 )

        do k = 0, n - 1

          call recur ( k, ck, dk, ek )

          if ( 0 .lt. k ) then
            h = - ek * h
          end if

          bb = ( ck * xnode + dk ) * bk2 + ek * bk1
          dd = ( ck * xnode + dk ) * dk2 + ek * dk1 + ck * bk2
          bk1 = bk2
          bk2 = bb
          dk1 = dk2
          dk2 = dd

          if ( bk2 .ne. 0.0D+00 ) then
            j = int ( log10 ( abs ( bk2 ) ) )
            if ( 2 .lt. abs ( j ) ) then
c
c  Scale to control overflow/underflow.
c
              iscale = iscale - 2 * j
              scale = 10.0D+00**j
              bk2 = bk2 / scale
              bk1 = bk1 / scale
              dk1 = dk1 / scale
              dk2 = dk2 / scale
            end if
          end if

          if ( h .ne. 0.0D+00 ) then
            j = int ( log10 ( abs ( h ) ) )
            if ( 2 .le. abs ( j ) ) then
              iscale = iscale + j
              h = h / 10.0D+00**j
            end if
          end if

        end do

        wt = h0 * h / dk2 / bk1

        if ( wt .ne. 0.0D+00 ) then
          itest = int ( log10 ( abs ( wt ) ) ) + iscale
          if ( - nexp .le. itest ) then
            wt = wt * 10.0D+00**iscale
          else
            wt = 0.0D+00
          end if
        end if
        return
      end if
c
c  General case.
c
      bk2 = 0.0D+00
      bk1 = 0.0D+00
      rk2 = 0.0D+00
      rk1 = 0.0D+00
      call recur ( n, ck, dk, ek )
      call recur ( n + 1, ck1, dk1, ek1 )
      h = 1.0D+00
      iscale = 0

      do k = n, 1, -1

        if ( m .le. k ) then

          rs = t(k-m) / h
c
c  Scale and adjust for possible overflow/underflow.
c
          if ( nexp .lt. iscale ) then
            rs = 0.0D+00
          else
            rs = rs / 10.0D+00**iscale
          end if
        else
          rs = 0.0D+00
        end if

        bb = rs + ( dk + xnode * ck ) * bk1 + ek1 * bk2
        bk2 = bk1
        bk1 = bb

        call recur ( k - 1, ckm1, dkm1, ekm1 )

        if ( n .ne. m ) then
          h = - h * ck / ek / ckm1
        end if

        bb = bb * ckm1
        wt = bb + ( dkm1 + xnode * ckm1 ) * rk1 + ek * rk2
        rk2 = rk1
        rk1 = wt
        ck1 = ck
        dk1 = dk
        ek1 = ek
        ck = ckm1
        dk = dkm1
        ek = ekm1

        if ( bk1 .ne. 0.0D+00 ) then

          j = int ( log10 ( abs ( bk1 ) ) )
c
c  Scale to control overflow/underflow.
c
          if ( 2 .lt. abs ( j ) ) then
            iscale = iscale + j
            scale = 10.0D+00**j
            bk1 = bk1 / scale
            bk2 = bk2 / scale
            rk1 = rk1 / scale
            rk2 = rk2 / scale
          end if
        end if

      end do

      wt = h0 * bb / wt

      return
      end
