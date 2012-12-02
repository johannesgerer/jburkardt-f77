      subroutine addchg ( dif, size, xnew, chg, n )

c*********************************************************************72
c
cc ADDCHG adds an increment vector to a vector.
c
c  Discussion:
c
c    It also computes the norm of the increment, and the updated vector.  
c    It is useful for controlling iterations.
c
c  Modified:
c
c    12 December 2008
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, real DIF, the sum of the absolute values of the entries
c    of CHG.
c   
c    Output, real SIZE, the sum of the absolute values of the entries
c    of the updated vector XNEW.
c   
c    Input/output, real XNEW(N).  On input, the vector to be updated,
c    on output, the updated vector.
c   
c    Input, real CHG(N), the vector of updates to XNEW.
c   
c    Input, integer N, the dimension of CHG and XNEW.
c
      implicit none

      integer n

      real chg(n)
      real dif
      integer i
      real size
      real xnew(n)

      dif = 0.0E+00
      size = 0.0E+00
      do i = 1, n
        dif = dif + abs ( chg(i) )
        xnew(i) = xnew(i) + chg(i)
        size = size + abs ( xnew(i) )
      end do

      return
      end
      subroutine ahess ( a, la, n, w )

c*********************************************************************72
c
cc AHESS balances a real matrix, and reduces it to upper Hessenberg form.
c
c  Discussion:
c
c    After AHESS has been called, a call to SIM will display 
c    the similarity transformation that was used.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c  Input:
c   
c    Input/output, real A(N*N+2*N+1).  On input, the matrix.  On output,
c    the Hessenberg matrix and other information.
c   
c    Input, integer LA, the leading dimension of the array.
c   
c    Input, integer N, the rank of the matrix.
c
c    Workspace, real W(N).
c
      implicit none

      integer n

      real a(n*n+2*n+1)
      integer c
      integer d
      integer e
      integer g
      integer h
      integer i
      integer j
      integer k
      integer l
      integer la
      integer m
      integer o
      integer p
      integer q
      real r
      real s
      real t
      real u
      real v
      real w(n)

      if ( n .lt. la ) then
        call pack ( a, la, n )
      end if

      i = n * n
      m = n + 1
      o = m + 1
      e = i + m
      j = i + 1

      call bal ( a, n, n, w, a(j) )

      do k = 1,n
        a(e+k) = w(k)
      end do
      v = w(1)
      j = m
      k = i

20    continue

      k = k - n

30    continue

      a(i+j) = a(i)
      i = i - 1
      if ( i .gt. k ) then
        go to 30
      end if

      j = j - 1
      if ( k .gt. 0 ) go to 20
      a(1) = 2231
      a(2) = n
      k = 4
      l = o
      d = 1
      c = 2

40    continue

      if ( c .ge. n ) go to 200
      p = k + 1

      do i = p,l
        if ( a(i) .ne. 0.0E+00 ) go to 60
      end do

      a(l+1) = 0.0E+00
      go to 190
60    t = abs ( a(k) )
      if ( t .ne. 0.0E+00 ) then
        u = 1.0E+00 / t
      end if
      r = 1.0E+00

      do j = i,l
        s = abs ( a(j) )
        if ( t .lt. s ) then
          u = 1.0E+00 / s
          r = 1.0E+00 + r * ( t * u )**2
          t = s
        else
          r = r + ( s * u )**2
        end if
      end do

      s = t * sqrt ( r )
      r = a(k)
      u = 1.0E+00 / sqrt ( s * ( s + abs ( r ) ) )
      if ( r .lt. 0.0E+00 ) then
        s = -s
      end if
      i = l
90    a(i+1) = u * a(i)
      i = i - 1
      if ( i .gt. k ) go to 90
      a(k) = - s
      a(p) = u * ( r + s )
      h = l

      do i = 1, n
        w(i) = 0.0E+00
      end do

110   h = h + m
      s = a(p)
      p = p + 1
      q = h - n
      do i = 1, d
        w(i) = w(i) + s*a(i+q)
      end do
      j = k - d
      t = 0.0E+00
      do i = c, n
        r = a(i+q)
        t = t + r * a(i+j)
        w(i) = w(i) + r * s
      end do
      a(h+1) = t

      if ( h .lt. e ) then
        go to 110
      end if

      t = 0.
      h = l + 1
      p = k + 1
      j = c - p

      do i = p, h
        t = t + w(i+j) * a(i)
      end do

      do i = c, n
        w(i) = w(i) - t*a(i-j)
      end do

      h = l
160   g = h + 2
      q = h + m
      h = h + c
      t = a(q+1)
      s = a(p)
      p = p + 1
      j = 1 - g
      do i = g, h
        a(i) = a(i) - w(i+j)*s
      end do
      i = h
      h = q
      q = k - i
      g = i + 1
      do i = g, h
        a(i) = a(i) - a(i+q)*t - w(i+j)*s
      end do
      if ( h .lt. e ) go to 160
190   k = k + o
      l = l + m
      d = c
      c = c + 1
      go to 40
200   a(e+1) = v
      return
      end
      function amag ( c )

c*********************************************************************72
c
cc AMAG returns the L1 norm of a complex number.
c
c  Discussion:
c
c    The L1 norm of the complex number A+Bi is |A| + |B|.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input, complex C, the number.
c
c    Output, real AMAG, the L1 norm of the number.
c
      implicit none

      real amag
      complex c

      amag = abs ( real ( c ) ) + abs ( aimag ( c ) )

      return
      end
      subroutine bal ( a, la, n, d, w )

c*********************************************************************72
c
cc BAL balances (rescales the rows and columns of) a real matrix.
c
c  Discussion:
c
c    This is useful for stability and accuracy.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c  Input:
c
c    Input/output, real A(LA,N), the NxN matrix.
c    On output, output A = D * input A * inverse ( D ).
c  
c    Input, integer LA, the leading dimension of A.
c
c    Input, integer N, the number of rows and columns in A.
c
c    Output, real D(N), the diagonal of the balancing matrix.
c
c    Workspace, real W(2*N).
c
      implicit none

      integer la
      integer n

      real a(la,n)
      real b
      real c
      real d(n)
      integer i
      integer j
      integer k
      integer l
      integer m
      real q
      real r
      real s
      real t
      real w(2*n)

      t = 1.0E+00
10    continue
      t = t + t
      if ( ( 1.0E+00 + t ) - t .eq. 1.0E+00 ) then
        go to 10
      end if

      b = 0.
20    b = b + 1
      if ( t+b .eq. t ) go to 20

      if ( t+2.*b .gt. t+b ) go to 30
      b = b + b
30    q = alog ( b )
      q = 0.5E+00 / q

      do i = 1, n
        d(i) = 1.0E+00
        w(i) = 1.0E+00
        w(i+n) = 0.0E+00
      end do

      m = n + 1
      l = n + n

      do j = 1, n
        do i = m, l
          w(i) = w(i) + abs ( a(i-n,j) )
        end do
      end do

60    continue

      l = 0

      do  110 j = 1, n

        c = 0.
        do i = 1, n
          s = a(i,j)*w(i)
          a(i,j) = s
          c = c + abs ( s )
        end do
        if ( c .eq. 0. ) go to 110
        r = w(j+n)
        if ( r .le. 0.0E+00 ) go to 110
        s = .5 + q * alog ( c / r )
        if ( s .lt. 0. ) go to 80
        i = s
        if ( i .eq. s ) i = i - 1
        go to 90
80      i = s - 1
90      t = b**i
        s = 1.0E+00 / t
        w(j) = 1.0E+00
        if ( t*r+s*c .gt. 0.95*(r+c) ) go to 110
        l = 1
        w(j) = t
        d(j) = d(j)*s

        do i = 1, n
          r = a(i,j)
          c = r*s
          k = i + n
          w(k) = ( w(k) - abs ( r ) ) + abs ( c )
          a(i,j) = c
        end do

        w(j+n) = t*w(j+n)

110   continue

      if ( l .eq. 1 ) go to 60

      do j = 1, n
        do i = 1, n
          a(i,j) = a(i,j) * w(i)
        end do
      end do

      return
      end
      subroutine basis ( b, lb, n, a, c )

c*********************************************************************72
c
cc BASIS computes an orthonormal basis for a collection of vectors
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, real B(LB,N), the orthonormal basis.  B can be identified with A.
c
c    Input, integer LB, the leading dimension of B.
c
c    Input, integer N, the number of vectors in the basis.
c
c    Input, real A(*), the factorization of the vectors computed by QR.
c
c    Input, real C, the cutoff value.
c
      implicit none

      integer lb

      real a(*)
      real b(lb,*)
      real c
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      integer o
      real t

      t = a(1)

      if ( abs ( t ) .ne. 3230 ) then
        write ( *, '(a)' ) ' '
        write(*,*)'basis: error - you must factor the array of,'
        write(*,*)'vectors, using qr, before calling basis.'
      end if

      m = a(2)
      n = a(3)

      if ( lb .lt. m ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: the leading dimension of argument b in'
        write(*,*) 'basis must be greater than or equal'
        write(*,*) 'to the number of components in each basis vector'
      end if

      k = 4
      l = min ( m, n ) - 1

      o = m + 1
      do j = 1, l
        if ( abs ( a(j+k-1) ) .le. c ) then
          go to 110
        end if
        do i = j, m
          b(i,j) = a(i+k)
        end do
        k = k + o
      end do

      j = l + 1

      if ( abs ( a(k+l) ) .le. c ) then
        go to 110
      end if

      if ( n .lt. m ) go to 80

      n = m

      do j = 1, m
        do i = 1, m
          b(i,j) = 0.0E+00
        end do
        b(j,j) = 1.0E+00
      end do

      return

80    do i = n, m
        b(i,n) = a(i+k)
      end do

      call hsr3(b,lb,m,n)

      return

110   continue

      do k = j, n
        do i = 1, m
          b(i,k) = 0.0
        end do
      end do

      n = j - 1
      if ( n .gt. 0 ) then
        call hsr3 ( b, lb, m, n )
      end if

      return
      end
      function bcon ( a, b )

c*********************************************************************72
c
cc BCON estimates the condition number of a band matrix.
c
c  Discussion:
c
c    The matrix must already have been factored by BFACT.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c    William Hager,
c    Condition Estimates,
c    SIAM Journal on Scientific and Statistical Computing,
c    Volume 5, Number 2, June 1984, pages 311-316.
c
c  Parameters:
c   
c    Input, real A(*), the factorization information from BFACT.
c   
c    Workspace, real B(N).
c
c    Output, real BCON, the condition number of the matrix.
c
      implicit none

      real a(*)
      real b(*)
      real bcon
      real c
      real d
      integer i
      integer j
      integer m
      integer n

      c = a(1)

      if ( abs ( c ) .ne. 1231 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor with bfact before',
     &  'estimating condition'
        stop
      end if

      if ( c .lt. 0.0E+00 ) then
        bcon = -1.0E+00
        return
      end if

      m = 0
      n = a(2)
      c = 1.0E+00 / a(2)
      do j = 1, n
        b(j) = c
      end do
      go to 50

30    do j = 1, n
        b(j) = 0.0E+00
      end do
      b(m) = 1.0E+00

50    call bsolve ( b, a, b )

      c = 0.0E+00
      do j = 1, n
        c = c + abs ( b(j) )
        if ( b(j) .lt. 0.0E+00 ) then
          b(j) = -1.0E+00
        else
          b(j) = 1.0E+00
        end if
      end do

      call btrans ( b, a, b )
c
c  I is the index of the largest entry of B.
c
      i = 1
      do j = 1, n
        if ( abs ( b(i) ) .lt. abs ( b(j) ) ) then
          i = j
        end if
      end do

      if ( m .eq. 0 ) go to 90
      if ( m .eq. i ) go to 100
      if ( d .ge. c ) go to 100
90    m = i
      d = c
      go to 30
100   c = c * a(3)
      c = max ( c, 1.0E+00 )
      bcon = c
      return
      end
      function bdet ( iexp, a )

c*********************************************************************72
c
cc BDET computes the determinant of a band matrix factored by BFACT.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, integer IEXP, a power of ten that is part of the
c    determinant.
c   
c    Input, real A(*), factorization information from BFACT.
c
c    Output, real BDET, the mantissa of the determinant.
c    Determinant = BDET * 10.0^IEXP.
c
      implicit none

      real a(*)
      real bdet
      real c
      real d
      real f
      real g
      integer h
      integer i
      integer iexp
      integer j
      integer k
      integer l
      integer m
      integer n

      iexp=0
      bdet=0.0E+00
      d = a(1)

      if ( abs ( d ) .ne. 1231 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor with bfact before',
     &  ' computing determinant'
        return
      end if

      iexp = 0
      if ( d .lt. 0.0E+00 ) go to 70
      n = a(2)
      if ( n .eq. 1 ) go to 80
      d = 1.0E+00
      f = 2.0E+00**64
      g = 1.0E+00 / f
      h = 64
      j = a(4) + a(5) + 1
      m = 1 + j + a(4)
      k = 6
      l = 5 - m + m*n
      n = 0

      do i = k, l, m
        n = n + 1
        if ( a(i) .gt. n ) d = -d
        d = d*a(i+j)
20      if ( abs ( d ) .lt. f ) go to 30
        iexp = iexp + h
        d = d*g
        go to 20
30      if ( abs ( d ) .gt. g ) go to 40
        iexp = iexp - h
        d = d*f
        go to 30
40      continue
      end do

      d = d*a(j+l+1)
      if ( iexp .ne. 0 ) go to 50
      bdet = d
      return
50    if ( d .eq. 0. ) go to 90
      c = alog10 ( abs ( d ) ) + iexp * alog10 ( 2.0 )
      iexp = c
      c = c - iexp
      if ( c .le. 0.0 ) go to 60
      c = c - 1
      iexp = iexp + 1
60    f = 10.0E+00**c
      if ( d .lt. 0.0E+00 ) f = -f
      bdet = f
      return
70    bdet = 0.0E+00
      return
80    bdet = a(7)
      return
90    iexp = 0
      go to 70
      end
      subroutine bfact ( a, la, n, l, u )

c*********************************************************************72
c
cc BFACT factors a band matrix with partial pivoting.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input/output, real A(5+N*(2*L+U+2)).  On input, information defining the
c    matrix.  On output, information defining the factored matrix.
c
c    Input, integer LA, the leading dimension of array A.
c
c    Input, integer N, the number of rows and columns in the matrix.
c
c    Input, integer L, U, the number of bands below and above
c    the diagonal.   
c   
      implicit none

      real a(*)
      integer c
      integer d
      integer e
      integer f
      integer g
      integer h
      integer i
      integer j
      integer k
      integer l
      integer la
      integer m
      integer n
      integer o
      integer p
      integer q
      integer r
      real s
      real t
      integer u
      real v

      p = 1 + l
      c = p + u
      f = la*n + c
      k = l
10    if ( k .eq. 0 ) go to 30
      f = f - la
      i = f
      j = f - k
      k = k - 1
20    a(i) = 0.
      i = i - 1
      if ( i .gt. j ) go to 20
      go to 10
30    f = la*n
      k = u
40    if ( k .eq. 0 ) go to 60
      f = f - la
      i = f
      j = f + k
      k = k - 1
50    i = i + 1
      a(i) = 0.
      if ( i .lt. j ) go to 50
      go to 40
60    if ( c .lt. la ) call rpack(a,la,c,n)
      if ( u .eq. 0 ) go to 100
      q = c + 1
      d = c - u
      j = n*c
      k = n
70    k = k - 1
      if ( k .lt. 0 ) go to 100
      j = j - d
      i = j - c
      e = j - u
      f = j - min ( u, k )
80    if ( j .eq. f ) go to 90
      a(j) = a(i)
      i = i - q
      j = j - 1
      go to 80
90    if ( j .eq. e ) go to 70
      a(j) = 0.
      j = j - 1
      go to 90
100   v = 0.0E+00
      e = c + p
      k = n
      j = 5 + e*n
      i = c*n - j
110   if ( j .eq. 5 ) go to 150
      s = 0.
      f = j - c
120   t = a(i+j)
      a(j) = t
      s = s + abs ( t )
      j = j - 1
      if ( j .gt. f ) go to 120
      if ( v .lt. s ) v = s
      i = i + p
      f = f - l
130   if ( j .eq. f ) go to 140
      a(j) = 0.
      j = j - 1
      go to 130
140   a(j) = k
      k = k - 1
      j = j - 1
      go to 110
150   a(1) = 1231
      a(2) = n
      a(3) = v
      a(4) = l
      a(5) = u
      i = 5 - l
      if ( l .eq. 0 ) go to 230
      c = c + l
      d = l + u
      r = p + u
      k = 0
160   k = k + 1
      i = i + e
      if ( k .eq. n ) go to 260
      m = i + 1
      q = i
      o = min ( l, n - k )
      p = i + o
      do j = m, p
        if ( abs ( a(j) ) .gt. abs ( a(q) ) ) then
          q = j
        end if
      end do
      j = i - r
      h = q - i
      a(j) = k + h
      t = a(q)
      if ( t .eq. 0. ) go to 220
      a(q) = a(i)
      a(i) = t
      do j = m, p
        a(j) = a(j) / t
      end do
      f = i + c * min ( d, n - k )
      g = c - o
190   m = p + g
      p = m + h
      t = a(p)
      a(p) = a(m)
      a(m) = t
      p = m + o
      if ( t .eq. 0. ) go to 210
      q = i - m
      m = m + 1
      do j = m, p
        a(j) = a(j) - t*a(j+q)
      end do

210   if ( p .lt. f ) go to 190
      go to 160
220   a(1) = -1231
      go to 160
230   j = 5 + e*n
240   i = i + e
      if ( a(i) .eq. 0. ) go to 250
      if ( i .lt. j ) go to 240
      return
250   a(1) = -1231
      return
260   if ( a(i) .eq. 0. ) go to 250
      return
      end
      subroutine bidag ( d, b, a, la, m, n )

c*********************************************************************72
c
cc BIDAG reduces a general matrix to bidiagonal form
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c  Input:
c   
c  D      ARRAY WITH AT LEAST N ELEMENTS
c   
c  B      ARRAY WITH AT LEAST M ELEMENTS
c   
c  A      ARRAY CONTAINING COEFFICIENT MATRIX
c   
c  LA     Input, INTEGER LA, leading dimension of array A
c   
c  M      ROW DIMENSION OF MATRIX STORED IN A
c   
c  N      COLUMN DIMENSION OF MATRIX STORED IN A
c   
c  Output:
c   
c  D      DIAGONAL OF BIDIAGONAL FORM
c   
c  B      SUPERDIAGONAL (IF M .GE. N) OR SUBDIAGONAL (IF N .GT. M) OF
c         BIDIAGONAL FORM
c   
c  A      THE HOUSEHOLDER VECTORS USED IN THE REDUCTION PROCESS
c
      implicit none

      integer la

      real a(la,*)
      real b(*)
      real d(*)
      integer h
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      real r
      real s
      real t
      real u

      k = 1
      h = 2
      if ( m .lt. n ) go to 220

      if ( m .le. 1 ) then
        d(1) = a(1,1)
        return
      end if

10    continue

      j = k
      k = h
      do i = k, m
        if ( a(i,j) .ne. 0. ) go to 40
      end do
      d(j) = a(j,j)
      a(j,j) = 0.0E+00
      do i = k, n
         d(i) = a(j,i)
      end do
      go to 110

40    t = abs ( a(j,j) )
      if ( t .ne. 0. ) u = 1.0E+00 / t
      r = 1.0E+00
      do 60 l = i, m
        s = abs ( a(l,j) )
        if ( s .le. t ) go to 50
        u = 1.0E+00 / s
        r = 1.0E+00 + r*(t*u)**2
        t = s
        go to 60
50      r = r + (s*u)**2
60    continue

      s = t * sqrt ( r )
      r = a(j,j)
      u = 1.0E+00 / sqrt ( s* ( s + abs ( r ) ) )
      if ( r .lt. 0. ) s = -s
      d(j) = -s

      a(j,j) = u*(r+s)
      do i = k, m
        a(i,j) = a(i,j)*u
      end do

      if ( k .gt. n ) return

      do l = k, n
        t = 0.0
        do i = j, m
          t = t + a(i,j)*a(i,l)
        end do
        a(j,l) = a(j,l) - t*a(j,j)
        d(l) = a(j,l)
        do i = k, m
          a(i,l) = a(i,l) - t*a(i,j)
        end do
      end do

110   h = k + 1
      if ( k .lt. n ) go to 120
      if ( k .gt. n ) return
      if ( m .eq. n ) go to 210
      b(j) = a(j,n)
      go to 10

120   do i = h, n
        if ( d(i) .ne. 0. ) go to 140
      end do
      b(j) = d(k)
      a(j,k) = 0.0E+00
      go to 10
140   t = abs ( d(k) )
      if ( t .ne. 0. ) u = 1.0E+00 / t
      r = 1.0E+00
      do 160 l = i, n
        s = abs ( d(l) )
        if ( s .le. t ) go to 150
        u = 1.0E+00 / s
        r = 1.0E+00 + r*(t*u)**2
        t = s
        go to 160
150     r = r + (s*u)**2
160   continue
      s = t * sqrt ( r )
      r = d(k)
      u = 1.0E+00 / sqrt ( s * ( s + abs ( r ) ) )
      if ( r .lt. 0. ) s = -s
      d(k) = u*(r+s)
      do i = h, n
        d(i) = d(i)*u
      end do
      b(j) = -s
      do i = k, m
        b(i) = 0.0
      end do

      do l = k, n
        t = d(l)
        a(j,l) = t
        do i = k, m
          b(i) = b(i) + t*a(i,l)
        end do
      end do

      do l = k, n
        t = d(l)
        do i = k, m
          a(i,l) = a(i,l) - t*b(i)
        end do
      end do

      go to 10
210   d(n) = a(n,n)
      b(n-1) = a(n-1,n)
      return

220   do i = k, n
        d(i) = a(k,i)
      end do
240   j = k
      k = h
      do i = k, n
        if ( d(i) .ne. 0. ) go to 260
      end do
      d(j) = a(j,j)
      a(j,j) = 0.0E+00
      go to 330

260   t = abs ( d(j) )
      if ( t .ne. 0. ) u = 1.0E+00 / t
      r = 1.0E+00
      do l = i, n
        s = abs ( d(l) )
        if ( s .gt. t ) then
          u = 1.0E+00 / s
          r = 1.0E+00 + r*(t*u)**2
          t = s
        else
          r = r + (s*u)**2
        end if
      end do

      s = t * sqrt ( r )
      r = d(j)
      u = 1.0E+00 / sqrt ( s * ( s + abs ( r ) ) )
      if ( r .lt. 0. ) s = -s
      d(j) = u*(r+s)
      do i = k, n
        d(i) = d(i)*u
      end do
      u = -s
      if ( k .gt. m ) go to 340
      do i = k, m
        b(i) = 0.
      end do

      do l = j, n
        t = d(l)
        a(j,l) = t
        do i = k, m
          b(i) = b(i) + t*a(i,l)
        end do
      end do

      do l = j, n
        t = d(l)
        do i = k, m
          a(i,l) = a(i,l) - t*b(i)
        end do
      end do

330   h = k + 1
      d(j) = u
      if ( k .lt. m ) go to 360
      if ( k .gt. m ) return
      b(j) = a(m,j)
      go to 220

340   do i = j, n
        a(j,i) = d(i)
      end do
      go to 330

360   do i = h, m
        if ( a(i,j) .ne. 0. ) go to 380
      end do

      b(j) = a(k,j)
      a(k,j) = 0.
      go to 240
380   continue

      t = abs ( a(k,j) )
      if ( t .ne. 0. ) u = 1.0E+00 / t
      r = 1.0E+00

      do l = i, m
        s = abs ( a(l,j) )
        if ( s .gt. t ) then
          u = 1.0E+00 / s
          r = 1.0E+00 + r*(t*u)**2
          t = s
        else
          r = r + (s*u)**2
        end if
      end do

      s = t * sqrt ( r )
      r = a(k,j)
      u = 1.0E+00 / sqrt ( s * ( s + abs ( r ) ) )
      if ( r .lt. 0. ) s = -s
      b(j) = -s
      a(k,j) = u*(r+s)
      do i = h, m
        a(i,j) = a(i,j)*u
      end do

      do l = k, n
        t = 0.
        do i = k, m
          t = t + a(i,j)*a(i,l)
        end do
        a(k,l) = a(k,l) - t*a(k,j)
        d(l) = a(k,l)
        do i = h, m
          a(i,l) = a(i,l) - t*a(i,j)
        end do
      end do

      go to 240
      end
      subroutine bidag2 ( d, b, q, lq, iq, p, lp, ip, a, la, m, n )

c*********************************************************************72
c
cc BIDAG2 reduces a general matrix to bidiagonal form.
c
c  Discussion:
c
c    A = Q TIMES BIDIAGONAL MATRIX TIMES P TRANSPOSE
c     
c    EITHER P OR Q CAN BE IDENTIFIED WITH A BUT NOT BOTH. WHEN EITHER P OR Q ARE
c    IDENTIFIED WITH A, THEN THE HOUSEHOLDER VECTORS IN A ARE DESTROYED
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, real D(N), the diagonal of the bidiagonal form.
c   
c    Output, real B(M), if N <= M, the superdiagonal, or if M < N, 
c    the subdiagonal of the bidiagonal form.
c
c    Output, real Q(LQ,*), the Q factor.
c
c    Input, integer LQ, leading dimension of array Q
c   
c    Input, integer IQ     AN INTEGER WHICH INDICATES WHICH COLUMNS OF Q TO COMPUTE
c         = 0 MEANS NONE,
c         = 1 MEANS FIRST L,
c         = 2 MEANS LAST M-L,
c         = 3 MEANS ALL M WHERE L = MIN(M,N))
c   
c    Output, real P(LP,*), the P factor.
c
c    Input, integer LP, leading dimension of array P
c   
c    Input, integer IP, INDICATES WHICH COLUMNS OF P TO COMPUTE.
c   
c    Input/output, real A(LA,*), on input the coefficient matrix.
c    On output, the Householder vectors used in the reduction process.
c   
c    Input, integer LA, the leading dimension of array A
c
c    Input, integer M, N, the row and column dimension of A.
c
      implicit none

      integer la
      integer lp
      integer lq

      real a(la,*)
      real b(*)
      real d(*)
      integer h
      integer i
      integer ip
      integer iq
      integer j
      integer jp
      integer jq
      integer k
      integer l
      integer m
      integer n
      real p(lp,*)
      real q(lq,*)
      real r
      real s
      real t
      real u

      l = min ( m, n )

      if ( iq .lt. 0 .or. 3 .lt. iq ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: input parameter iq for bidag2'
        write(*,*) 'either less than 0 or greater than 3'
        stop
      end if

      jq = iq
      if ( iq .le. 1 ) go to 30
      if ( iq .eq. 3 ) go to 30
      if ( m .eq. l ) jq = 0
30    if ( ip .ge. 0 ) go to 50
40    continue
      write ( *, '(a)' ) ' '
      write(*,*) 'error: input parameter ip for bidag2'
      write(*,*) 'either less than 0 or greater than 3'
      stop
50    if ( ip .gt. 3 ) go to 40
      jp = ip
      if ( ip .le. 1 ) go to 60
      if ( ip .eq. 3 ) go to 60
      if ( n .eq. l ) jp = 0
60    k = 1
      h = 2
      if ( m .lt. n ) go to 330
      if ( m .gt. 1 ) go to 70
      d(1) = a(1,1)
      if ( iq .gt. 0 ) q(1,1) = 1.0E+00
      if ( ip .gt. 0 ) p(1,1) = 1.0E+00
      return
70    j = k
      k = h
      do i = k, m
        if ( a(i,j) .ne. 0. ) go to 110
      end do
      d(j) = a(j,j)
      a(j,j) = 0.
      do i = k, n
        d(i) = a(j,i)
      end do
      if ( jq .eq. 0 ) go to 200

      do i = j, m
        q(i,j) = 0.0
      end do

      go to 200

110   t = abs ( a(j,j) )
      if ( t .ne. 0. ) u = 1.0E+00 / t
      r = 1.0E+00
      do l = i, m
        s = abs ( a(l,j) )
        if ( s .gt. t ) then
          u = 1.0E+00 / s
          r = 1.0E+00 + r*(t*u)**2
          t = s
        else
          r = r + (s*u)**2
        end if
      end do

      s = t * sqrt ( r )
      r = a(j,j)
      u = 1.0E+00 / sqrt ( s * ( s + abs ( r ) ) )
      if ( r .lt. 0. ) s = -s
      d(j) = -s
      a(j,j) = u*(r+s)
      do i = k, m
        a(i,j) = a(i,j)*u
      end do
      if ( jq .eq. 0 ) go to 160
      do i = j, m
        q(i,j) = a(i,j)
      end do
160   if ( k .gt. n ) go to 620

      do l = k, n
        t = 0.
        do i = j, m
          t = t + a(i,j)*a(i,l)
        end do
        a(j,l) = a(j,l) - t*a(j,j)
        d(l) = a(j,l)
        do i = k, m
          a(i,l) = a(i,l) - t*a(i,j)
        end do
      end do

200   h = k + 1
      if ( k .lt. n ) go to 210
      if ( k .gt. n ) go to 620
      if ( m .eq. n ) go to 610
      b(j) = a(j,n)
      go to 70
210   do i = h, n
        if ( d(i) .ne. 0. ) go to 240
      end do
      b(j) = d(k)
      a(j,k) = 0.
      if ( ip .eq. 0 ) go to 70
      do i = k, n
        p(i,j) = 0.
      end do
      go to 70
240   t = abs ( d(k) )
      if ( t .ne. 0. ) u = 1.0E+00 / t
      r = 1.0E+00
      do 260 l = i, n
        s = abs ( d(l) )
        if ( s .le. t ) go to 250
        u = 1.0E+00 / s
        r = 1.0E+00 + r*(t*u)**2
        t = s
        go to 260
250     r = r + (s*u)**2
260   continue
      s = t * sqrt ( r )
      r = d(k)
      u = 1.0E+00 / sqrt ( s * ( s + abs ( r ) ) )
      if ( r .lt. 0. ) s = -s
      d(k) = u*(r+s)
      do i = h, n
        d(i) = d(i)*u
      end do
      if ( ip .eq. 0 ) go to 290
      do i = k, n
        p(i,j) = d(i)
      end do
290   b(j) = -s
      do i = k, m
        b(i) = 0.
      end do

      do l = k, n
        t = d(l)
        a(j,l) = t
        do i = k, m
          b(i) = b(i) + t*a(i,l)
        end do
      end do

      do l = k, n
        t = d(l)
        do i = k, m
          a(i,l) = a(i,l) - t*b(i)
        end do
      end do

      go to 70

330   do i = k, n
        d(i) = a(k,i)
      end do
350   j = k
      k = h
      do i = k, n
        if ( d(i) .ne. 0. ) go to 370
      end do
      u = d(j)
      d(j) = 0.
      go to 440
370   t = abs ( d(j) )
      if ( t .ne. 0. ) u = 1.0E+00 / t
      r = 1.0E+00
      do 390 l = i, n
        s = abs ( d(l) )
        if ( s .le. t ) go to 380
        u = 1.0E+00 / s
        r = 1.0E+00 + r*(t*u)**2
        t = s
        go to 390
380     r = r + (s*u)**2
390   continue
      s = t * sqrt ( r )
      r = d(j)
      u = 1.0E+00 / sqrt ( s * ( s + abs ( r ) ) )
      if ( r .lt. 0. ) s = -s
      d(j) = u*(r+s)
      do i = k, n
        d(i) = d(i)*u
      end do
      u = -s

      if ( k .gt. m ) go to 470

      do i = k, m
        b(i) = 0.
      end do

      do l = j, n
        t = d(l)
        a(j,l) = t
        do i = k, m
          b(i) = b(i) + t*a(i,l)
        end do
      end do

      do l = j, n
        t = d(l)
        do i = k, m
          a(i,l) = a(i,l) - t*b(i)
        end do
      end do

440   h = k + 1

      if ( ip .eq. 0 ) go to 460

      do i = j, n
        p(i,j) = d(i)
      end do

460   d(j) = u
      if ( k .lt. m ) go to 490
      if ( k .gt. m ) go to 620
      b(j) = a(m,j)
      go to 330

470   do i = j, n
        a(j,i) = d(i)
      end do
      go to 440

490   do i = h, m
        if ( a(i,j) .ne. 0. ) go to 520
      end do
      b(j) = a(k,j)
      a(k,j) = 0.
      if ( iq .eq. 0 ) go to 330
      do i = k, m
        q(i,j) = 0.
      end do
      go to 330

520   t = abs ( a(k,j) )
      if ( t .ne. 0. ) u = 1.0E+00 / t
      r = 1.0E+00
      do 540 l = i, m
        s = abs ( a(l,j) )
        if ( s .le. t ) go to 530
        u = 1.0E+00 / s
        r = 1.0E+00 + r*(t*u)**2
        t = s
        go to 540
530     r = r + (s*u)**2
540   continue 
      s = t * sqrt ( r )
      r = a(k,j)
      u = 1.0E+00 / sqrt ( s * ( s + abs ( r ) ) )
      if ( r .lt. 0. ) s = -s
      b(j) = -s
      a(k,j) = u*(r+s)
      do i = h, m
        a(i,j) = a(i,j)*u
      end do
      if ( iq .eq. 0 ) go to 570
      do i = k, m
        q(i,j) = a(i,j)
      end do

570   do l = k, n
        t = 0.
        do i = k, m
          t = t + a(i,j)*a(i,l)
        end do
        a(k,l) = a(k,l) - t*a(k,j)
        d(l) = a(k,l)
        do i = h, m
          a(i,l) = a(i,l) - t*a(i,j)
        end do
      end do

      go to 350
610   d(n) = a(n,n)
      b(n-1) = a(n-1,n)
620   if ( jq .eq. 0 ) go to 650
      if ( n .gt. m ) go to 640
      if ( n .eq. m ) go to 630
      if ( jq .eq. 1 ) call hsr3(q,lq,m,n)
      if ( jq .eq. 2 ) call hsr4(q,lq,m,n)
      if ( jq .eq. 3 ) call hsr5(q,lq,m,n)
      go to 650
630   call hsr2(q,lq,m)
      go to 650
640   call hsr1(q,lq,m)
650   if ( jp .eq. 0 ) return
      if ( n .le. m ) go to 660
      if ( jp .eq. 1 ) call hsr3(p,lp,n,m)
      if ( jp .eq. 2 ) call hsr4(p,lp,n,m)
      if ( jp .eq. 3 ) call hsr5(p,lp,n,m)
      return
660   call hsr1(p,lp,n)
      return
      end
      subroutine bmult ( y, x, a, la, n, l, u )

c*********************************************************************72
c
cc BMULT multiplies a real band matrix A by a vector X.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c  Input:
c
c    Output, real Y(N), the product of the matrix times X.
c
c    Input, real X(N), the vector to be multiplied.
c
c    Input, real A(LA,N), containing the (L+1+U)xN array of
c    bands.
c   
c    Input, integer LA, the leading dimension of array A.
c
c    Input, integer N, the number of rows and columns of A.
c
c    Input, integer L, U, the number of bands below and above the diagonal.
c
      implicit none

      integer la
      integer n

      real a(la,*)
      integer i
      integer j
      integer k
      integer l
      integer m
      integer o
      integer p
      integer q
      real t
      integer u
      real x(n)
      real y(n)

      m = u + 1

      if ( l .le. 0 ) then

        do j = 1, n
          t = 0.0E+00
          k = min ( n, j + u )
          m = m + 1
          do i = j, k
            t = t + a(m-i,j) * x(i)
          end do
          y(j) = t
        end do

      else

        o = m
        p = n - 1

        do i = 1, p
          y(i) = 0.0E+00
        end do

        y(n) = a(m,n) * x(n)

        do j = 1, p

          t = 0.0E+00
          k = min ( n, j + u )
          m = m + 1
          do i = j, k
            t = t + a(m-i,j) * x(i)
          end do
          y(j) = y(j) + t
          k = j + 1
          q = min ( n, j + l )
          o = o - 1
          t = x(j)
          do i = k, q
            y(i) = y(i) + a(o+i,j) * t
          end do

        end do

      end if

      return
      end
      subroutine bsolve ( x, a, b )

c*********************************************************************72
c
cc BSOLVE solves a factored band system.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, real X(N), the solution.
c
c    Input, real A(*), factorization information from BFACT.
c
c    Input, real B(N), the right hand side.
c
      implicit none

      real a(*)
      real b(*)
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      integer o
      integer p
      integer q
      real t
      real x(*)

      t = a(1)

      if ( abs ( t ) .ne. 1231 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor with bfact before solving'
        stop
      end if

      n = a(2)
      l = a(4)
      m = a(5)
      o = l + m
      m = 1 + l + o
      j = 5 - l
      k = 0
      if ( t .lt. 0.0E+00 ) go to 90

      do i = 1, n
        x(i) = b(i)
      end do

      q = 1
      if ( l .gt. 0 ) go to 40
      k = n
      j = j + m*n
      if ( m .gt. 1 ) go to 60

      do k = 1, n
        x(k) = x(k) / a(5+k+k)
      end do

      return
40    j = j + m
      i = a(j+k-o)
      k = q
      q = k + 1
      if ( k .eq. n ) go to 60
      t = x(i)
      x(i) = x(k)
      x(k) = t
      if ( t .eq. 0 ) go to 40
      p = min ( k + l, n )
      do i = q, p
        x(i) = x(i) - t*a(i+j)
      end do
      go to 40
60    t = x(k) / a(j+k)
70    x(k) = t
      if ( k .eq. 1 ) return
      q = max ( 1, k - o )
      k = k - 1
      do i = q, k
        x(i) = x(i) - t*a(i+j)
      end do
      j = j - m
      go to 60
90    j = j + m
      k = k + 1
      if ( a(j+k) .ne. 0. ) go to 90

      do i = 1, n
        x(i) = 0.0
      end do

      t = 1.0E+00
      if ( m .gt. 1 ) go to 70
      x(k) = 1.0E+00
      return
      end
      subroutine btrans ( x, a, b )

c*********************************************************************72
c
cc BTRANS solves transpose of a factored band system.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, real X(N), the solution.
c
c    Input, real A(*), factorization information from BFACT.
c   
c    Input, real B(N), the right hand side.
c   
      implicit none

      real a(*)
      real b(*)
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      integer o
      integer p
      integer q
      real t
      real x(*)

      t = a(1)

      if ( abs ( t ) .ne. 1231 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor with bfact before solving'
        stop
      end if

      n = a(2)
      l = a(4)
      m = a(5)
      o = l + m - 1
      m = 2 + l + o
      j = 7 + o
      k = 1
      if ( t .lt. 0. ) go to 100
      t = 0.
      if ( o .ge. 0 ) go to 40
      do k = 1, n
        x(k) = b(k) / a(k+k+5)
      end do
      return
30    if ( b(k) .ne. 0. ) go to 40
      x(k) = 0.
      k = k + 1
      if ( k .le. n ) go to 30
      return
40    j = j - m + m*k
50    x(k) = (b(k)-t) / a(j+k)
      if ( k .eq. n ) go to 70
      t = 0.
      j = j + m
      p = max ( 1, k - o )
      do i = p, k
        t = t + x(i)*a(i+j)
      end do
      k = k + 1
      go to 50
70    if ( l .eq. 0 ) return
      o = o + 2
80    if ( k .eq. 1 ) return
      j = j - m
      q = k
      k = k - 1
      p = min ( n, k + l )
      t = x(k)
      do i = q, p
        t = t - x(i)*a(i+j)
      end do
      i = a(j+k-o)
      x(k) = x(i)
      x(i) = t
      go to 80
100   i = 8 + o + n + m*n
      q = n + 1
110   i = i - m - 1
      q = q - 1
      if ( a(i) .ne. 0. ) go to 110
      j = j + m*(q-k)
      k = q
      do i = 1, n
        x(i) = 0.
      end do
      x(k) = 1.0E+00
      if ( o .lt. 0 ) return
130   if ( k .eq. n ) go to 70
      t = 0.
      j = j + m
      p = max ( q, k - o )
      do i = p, k
        t = t - x(i)*a(i+j)
      end do
      k = k + 1
      x(k) = t / a(j+k)
      go to 130
      end
      subroutine bvert ( v, lv, a )

c*********************************************************************72
c
cc BVERT inverts a band matrix
c
c  Discussion:
c
c    The inverse of a band matrix is generally a full matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, real V(LV,N), the NxN inverse of A.
c
c    Input, integer LV, the leading dimension of array V.
c   
c    Input, real A(*), factorization information from BFACT.
c
      implicit none

      integer lv

      real a(*)
      integer g
      integer h
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      integer o
      integer p
      integer q
      real t
      real v(lv,*)

      t = a(1)

      if ( abs ( t ) .ne. 1231 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor with bfact before inverting'
        stop
      end if

      if ( t .le. 0. ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: matrix has no inverse'
        stop
      end if

      n = a(2)
      l = a(4)
      m = a(5)
      o = l + m
      m = 1 + l + o
      g = 5 - l - m

      do 80 h = 1, n
        j = g + m*h
        k = h - 1
        q = h
        do i = 1, n
          v(i,h) = 0.
        end do
        v(h,h) = 1.0E+00
        if ( l .gt. 0 ) go to 40
        k = n
        j = 5 - l + m*n
        if ( m .gt. 1 ) go to 60
        v(h,h) = v(h,h) / a(5+h+h)
        go to 80
40      j = j + m
        i = a(j+k-o)
        k = q
        q = k + 1
        if ( k .eq. n ) go to 60
        t = v(i,h)
        v(i,h) = v(k,h)
        v(k,h) = t
        if ( t .eq. 0 ) go to 40
        p = min ( k + l, n )
        do i = q, p
          v(i,h) = v(i,h) - t*a(i+j)
        end do
        go to 40
60      t = v(k,h) / a(j+k)
        v(k,h) = t
        if ( k .eq. 1 ) return
        q = max ( 1, k - o )
        k = k - 1
        do i = q, k
          v(i,h) = v(i,h) - t*a(i+j)
        end do
        j = j - m
        go to 60
80    continue

      return
      end
      subroutine cahess ( a, la, n, w )

c*********************************************************************72
c
cc CAHESS balances a complex matrix, and reduces it to Upper Hessenberg form.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input/output, complex A(1+N*(N+2)).  On input, the matrix.
c    On output, the Hessenberg matrix.
c
c    Input, INTEGER LA, leading dimension of array A
c
c    Input, integer N, the order of the matrix.
c
c    Workspace, complex W(N).
c
      implicit none

      integer n

      complex a(*)
      real amag
      integer c
      integer d
      integer e
      integer g
      integer h
      integer i
      integer j
      integer k
      integer l
      integer la
      integer m
      integer o
      integer p
      integer q
      real r
      real s
      real sqr
      real t
      real u
      complex v
      complex w(n)
      complex x
      complex y
      complex z

      if ( la .gt. n ) then
        call cpack(a,la,n)
      end if

      i = n*n
      u = 0.
      m = n + 1
      o = m + 1
      e = i + m
      j = i + 1
      call cball(a,n,n,w,a(j))
      j = m
      k = i
      c = 0
      g = e
      h = e + n
10    c = c + 1
      g = g + 1
      a(g) = real ( w(c) )
      if ( g .eq. h ) go to 20
      g = g + 1
      a(g) = aimag(w(c))
      if ( g .lt. h ) go to 10
20    v = a(e+1)
30    k = k - n
40    a(i+j) = a(i)
      i = i - 1
      if ( i .gt. k ) go to 40
      j = j - 1
      if ( k .gt. 0 ) go to 30
      a(1) = 2233
      a(2) = n
      k = 4
      l = o
      d = 1
      c = 2
50    if ( c .ge. n ) go to 210
      p = k + 1
      do i = p, l
        if ( amag ( a(i) ) .ne. 0. ) go to 70
      end do
      a(l+1) = (0.,0.)
      go to 200
70    t = amag ( a(k) )
      if ( t .ne. 0. ) u = 1.0E+00 / t 
      r = sqr(a(k),u)
      do 90 j = i, l
        s = amag ( a(j) )
        if ( s .le. t ) go to 80
        u = 1.0E+00 / s
        r = sqr(a(j),u) + r*(t*u)**2
        t = s
        go to 90
80      r = r + sqr(a(j),u)
90    continue
      s = t * sqrt ( r )
      z = a(k)
      t = cabs ( z )
      u = 1.0E+00 / sqrt ( s * ( s + t ) )
      if ( t .ne. 0. ) z = z / t
      if ( t .eq. 0. ) z = (1.0E+00,0.)
      i = l
100   a(i+1) = u * conjg ( a(i) )
      i = i - 1
      if ( i .gt. k ) go to 100
      a(k) = -z*s
      a(p) = conjg ( z ) * u * ( t + s )
      h = l
      do i = 1, n
        w(i) = (0.,0.)
      end do
120   h = h + m
      y = conjg ( a(p) )
      p = p + 1
      q = h - n
      do i = 1, d
        w(i) = w(i) + y*a(i+q)
      end do
      j = k - d
      z = (0.,0.)

      do i = c, n
        x = a(i+q)
        z = z + x*a(i+j)
        w(i) = w(i) + x*y
      end do

      a(h+1) = z
      if ( h .lt. e ) go to 120
      z = (0.,0.)
      h = l + 1
      p = k + 1
      j = c - p
      do i = p, h
        z = z + w(i+j)*a(i)
        a(i) = conjg ( a(i) )
      end do
      do i = c, n
        w(i) = w(i) - z*a(i-j)
      end do
      h = l
170   g = h + 2
      q = h + m
      h = h + c
      z = a(q+1)
      y = conjg ( a(p) )
      p = p + 1
      j = 1 - g
      do i = g, h
        a(i) = a(i) - w(i+j)*y
      end do
      i = h
      h = q
      q = k - i
      g = i + 1
      do i = g, h
        a(i) = a(i) - a(i+q)*z - w(i+j)*y
      end do
      if ( h .lt. e ) go to 170
200   k = k + o
      l = l + m
      d = c
      c = c + 1
      go to 50
210   a(e+1) = v
      return
      end
      subroutine cball ( a, la, n, d, w )

c*********************************************************************72
c
cc CBALL balances (rescales the rows and columns of) a complex matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    A     --COMPLEX ARRAY CONTAINING MATRIX
c
c    LA    --LEADING (ROW) DIMENSION OF ARRAY A
c
c    N     --DIMENSION OF MATRIX STORED IN A
c
c    W     --WORK ARRAY (AT LEAST 2N REAL ELEMENTS)
c
c  OUTPUT:
c
c    A     --BALANCED ARRAY
c    (NEW A = D TIMES OLD A TIMES D SUP -1)
c
c    D     --REAL ARRAY STORING DIAGONAL OF D MATRIX
c
      implicit none

      integer la

      complex a(la,*)
      real amag
      real b
      real c
      real d(*)
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      real q
      real r
      real s
      real t
      real w(1)
      complex y
      complex z

      t = 1.0E+00
10    t = t + t
      if ( (1.0E+00 + t ) - t .eq. 1.0E+00 ) go to 10
      b = 0.0E+00
20    b = b + 1
      if ( t+b .eq. t ) go to 20
      if ( t+2.0E+00*b .gt. t+b ) go to 30
      b = b + b
30    q = alog ( b )
      q = 0.5E+00 / q

      do i = 1, n
        d(i) = 1.0E+00
        w(i) = 1.0E+00
        w(i+n) = 0.0E+00
      end do

      m = n + 1
      l = n + n

      do j = 1, n
        do i = m,l
          w(i) = w(i) + amag ( a(i-n,j) )
        end do
      end do

60    l = 0
      do  110 j = 1, n
        c = 0.0E+00
        do i = 1, n
          z = a(i,j)*w(i)
          a(i,j) = z
          c = c + amag ( z )
        end do
        if ( c .eq. 0.0E+00 ) go to 110
        r = w(j+n)
        if ( r .le. 0.0E+00 ) go to 110
        s = .5 + q * alog ( c / r )
        if ( s .lt. 0.0E+00 ) go to 80
        i = s
        if ( i .eq. s ) i = i - 1
        go to 90
80      i = s - 1
90      t = b**i
        s = 1.0E+00 / t 
        w(j) = 1.0E+00
        if ( t*r+s*c .gt. .95*(r+c) ) go to 110
        l = 1
        w(j) = t
        d(j) = d(j)*s

        do i = 1, n
          z = a(i,j)
          y = z*s
          k = i + n
          w(k) = ( w(k) - amag ( z ) ) + amag ( y )
          a(i,j) = y
        end do

        w(j+n) = t*w(j+n)
110   continue

      if ( l .eq. 1 ) go to 60

      do j = 1, n
        do i = 1, n
          a(i,j) = a(i,j)*w(i)
        end do
      end do

      return
      end
      subroutine cc ( z, f, e, n, y, y0 )

c*********************************************************************72
c
cc CC is used by CZERO.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      integer n

      real amag
      real e(n)
      complex f(n)
      integer i
      real s
      real t
      complex x
      real y
      real y0
      complex z(n)

      do i = 1, n
        s = 0.0
        x = z(i)
        t = amag ( x )
        if ( t .gt. 1.0E+00 ) go to 20
        if ( t .gt. y0 ) go to 30
        if ( t .eq. 0.0 ) go to 30
10      t = t*y
        x = x*y
        s = s - 1.0E+00
        if ( t .le. y0 ) go to 10
        go to 30
20      t = t*y0
        x = x*y0
        s = s + 1.0E+00
        if ( t .gt. 1.0E+00 ) go to 20
30      e(i) = s
        f(i) = x
      end do

      return
      end
      subroutine cdiag ( e, v, lv, a, la, n )

c*********************************************************************72
c
cc CDIAG diagonalizes a general complex matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
C     |    INPUT:                                              |
C     |                                                        |
c    Input, integer LV, the leading dimension of the array V.
c
C     |                                                        |
C     |         A     --COMPLEX ARRAY CONTAINING COEFFICIENT   |
C     |                 MATRIX (LENGTH AT LEAST 1 + N(N+2))    |
C     |                                                        |
C     |         LA    --LEADING (ROW) DIMENSION OF ARRAY A     |
C     |                                                        |
C     |         N     --DIMENSION OF MATRIX STORED IN A        |
C     |                                                        |
C     |    OUTPUT:                                             |
C     |                                                        |
C     |         E     --COMPLEX ARRAY OF EIGENVALUES           |
C     |                                                        |
C     |         V     --COMPLEX ARRAY OF EIGENVECTORS 
c
      implicit none

      integer lv

      complex a(*)
      complex e(*)
      integer i
      integer j
      integer k
      integer l
      integer la
      integer m
      integer n
      complex v(lv,*)

      call cahess(a,la,n,e)
      call csim(v,lv,a)
      j = 2
      k = 1
      m = 1
      l = 2

10    continue

      do i = m, l
        a(i) = a(i+j)
      end do

      m = l + 1
      j = j + n - k
      k = k + 1
      l = m + k
      if ( k .lt. n ) go to 10
      if ( k .gt. n ) go to 30
      l = l - 1
      go to 10

30    i = (n*(n+3)) / 2
      j = i + n
      k = j + 1 + n / 2
      call dag(e,v,lv,a,n,a(k),a(j),a(i))
      return
      end
      subroutine cediag ( e, v, lv, a, n, w )

c*********************************************************************72
c
cc CEDIAG diagonalizes a complex Hessenberg matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
C     |    INPUT:                                              |
C     |                                                        |
c    Input, integer LV, the leading dimension of the array V.
c
C     |         A     --COEFFICIENTS OF HESSENBERG MATRIX      |
C     |                 PACKED AT START OF COMPLEX ARRAY       |
c
c    Input, integer N, the dimension of the matrix.
c
C     |         W     --REAL WORK ARRAY WITH AT LEAST          |
C     |                 4N ELEMENTS                            |
C     |                                                        |
C     |    OUTPUT:                                             |
C     |                                                        |
C     |         E     --COMPLEX ARRAY OF EIGENVALUES           |
C     |                                                        |
C     |         V     --COMPLEX ARRAY OF EIGENVECTORS
c
      implicit none

      integer lv

      complex a(*)
      complex e(*)
      integer i
      integer j
      integer l
      integer m
      integer n
      complex v(lv,*)
      real w(*)

      do j = 1, n
        do i = 1, n
          v(i,j) = 0.0E+00
        end do
        v(j,j) = 1.0E+00
      end do

      m = n + 1
      l = m + n
      call dag(e,v,lv,a,n,w,w(m),w(l))

      return
      end
      subroutine cef ( w, z, f, e, d, r, x, n, n1, t, t1, y, y0 )

c*********************************************************************72
c
cc CEF
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      integer n

      real a
      real amag
      real b
      complex c
      real d(*)
      real e(*)
      complex f(*)
      integer i
      integer n1
      complex p
      real q
      real r
      real s
      real t
      real t1
      complex w(*)
      complex x
      real y
      real y0
      complex z(n)

      p = 1.0E+00
      q = 0.0E+00
      s = e(1)
      c = x*r
      if ( r .gt. 1.0E+00 ) go to 50
      t = 1.0E+00
      t1 = 1.0E+00

      do 20 i = 1, n
        a = e(i) + q
        if ( a .gt. s ) s = a
        w(i) = p*f(i)
        d(i) = a
        t = r*t + amag ( z(n1-i) )
        p = p*c
        if ( amag ( p ) .gt. y0 ) go to 20
10      p = p*y
        q = q - 1.0E+00
        if ( amag ( p ) .le. y0 ) go to 10
20    continue

30    if ( q .gt. s ) s = q
      do i = 1, n
        a = d(i) - s
        w(i) = w(i)*y**a
      end do
      a = q - s
      w(n1) = p*y**a
      return
50    t1 = r
      b = 1.0E+00 / r
      t = 0.0E+00
      do 70 i = 1, n
        a = e(i) + q
        if ( a .gt. s ) s = a
        w(i) = p*f(i)
        d(i) = a
        t = b*t + amag ( z(i) )
        p = p*c
        if ( amag ( p ) .le. 1.0E+00 ) go to 70
60      p = p*y0
        q = q + 1.0E+00
        if ( amag ( p ) .gt. 1.0E+00 ) go to 60
70    continue
      t = b*t + 1.0E+00
      go to 30
      end
      subroutine ceig ( e1, e2, a, b, c, d ) 

c*********************************************************************72
c
cc CEIG
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      complex a
      real amag
      complex b
      complex c
      complex d
      complex e1
      complex e2
      real p
      real q
      real r
      real s
      real t
      complex u
      complex v
      complex x
      complex y

      q = amag ( b )
      r = amag ( c )
      t = sqrt ( q ) * sqrt ( r )
      u = .5*(a-d)
      s = amag ( u )
      if ( s .lt. t ) go to 20
      if ( t .ne. 0. ) go to  10
      e1 = a
      e2 = d
      return
10    p = 1.0E+00 / t
      s = 1.0E+00 / s
      v = csqrt ( ( s * u )**2 + ((b*p)*(c*p))*(s*t)**2)
      go to 30
20    s = 1.0E+00 / t
      v = csqrt ( (s*u)**2+(b*s)*(c*s))
30    u = s*u
      x = u + v
      y = u - v
      if ( amag ( x ) .gt. amag ( y ) ) go to 50
      if ( q .gt. r ) go to 40
      e1 = a + b*(c*s) / y
      e2 = d - b*(c*s) / y
      return
40    e1 = a + (b*s)*c / y
      e2 = d - (b*s)*c / y
      return
50    if ( q .gt. r ) go to 60
      e1 = a + b*(c*s) / x
      e2 = d - b*(c*s) / x
      return
60    e1 = a + (b*s)*c / x
      e2 = d - (b*s)*c / x

      return
      end
      subroutine cemult ( y, x, a, n ) 

c*********************************************************************72
c
cc CEMULT mltiplies a complex Hessenberg matrix by a complex vector.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, complex Y(N), the product A*X.
c
c    Input, complex X(N), the vector to be multiplied.
c
c    Input, complex A((N*(N+3))/2-1), the matrix.
c
c    Input, integer N, the number of rows and columns of A.
c
      implicit none

      integer n

      complex a(*)
      integer i
      integer j
      integer k
      complex t
      complex x(n)
      complex y(n)

      if ( n .eq. 1 ) then
        y(1) = a(1) * x(1)
        return
      end if

      y(1) = ( 0.0E+00, 0.0E+00 )
      j = 0
      k = 1

10    continue

        t = x(k)
        do i = 1, k
          y(i) = y(i) + t * a(i+j)
        end do

        if ( k .eq. n ) then
          go to 20
        end if

        k = k + 1
        j = j + k
        y(k) = t * a(j)

      go to 10

20    continue

      return
      end
      subroutine cevals ( e, a, n, w )

c*********************************************************************72
c
cc CEVALS computes all eigenvalues of a complex Hessenberg matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
C     |    INPUT:                                              |
C     |                                                        |
C     |         A     --COEFFICIENTS OF HESSENBERG MATRIX      |
C     |                 PACKED AT START OF COMPLEX ARRAY       |
C     |                                                        |
C     |         N     --DIMENSION OF MATRIX STORED IN A        |
C     |                                                        |
C     |         W     --REAL WORK ARRAY WITH AT LEAST          |
C     |                 4N ELEMENTS                            |
C     |                                                        |
C     |    OUTPUT:                                             |
C     |                                                        |
C     |         E     --COMPLEX ARRAY OF EIGENVALUES 
c
      implicit none

      complex a(*)
      complex e(*)
      integer l
      integer m
      integer n
      real w(*)

      m = n + 1
      l = m + n
      call vls(e,a,n,w,w(m),w(l))

      return
      end
      subroutine cevect ( e, x, a, n )

c*********************************************************************72
c
cc CEVECT computes an eigenvector of a complex Hessenberg matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input/output, complex E, an estimated eigenvalue.  On output,
c    the estimate may have been improved.
c
C     |    INPUT:                                              |
C     |                                                        |
C     |         A     --HESSENBERG MATRIX PACKED AT START OF   |
C     |                 COMPLEX ARRAY WITH AT LEAST N(N+7)/2 -2|
C     |                 ELEMENTS (ALGORITHM DESTROYS ELEMENTS) |
C     |                                                        |
C     |    OUTPUT:                                             |
C     |                                                        |
C     |         X     --ARRAY CONTAINING EIGENVECTOR    
c
c    Input, integer N, the order of the matrix.
c
      implicit none

      complex a(*)
      real amag
      complex e
      integer f
      integer g
      integer h
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      integer o
      integer p
      integer q
      real r
      real s
      real sq
      real t
      complex x(*)
      complex y
      complex z

      if ( n .le. 1 ) then
        e = a(1)
        x(1) = (1.0E+00,0.)
        return
      end if

      g = (n*(n+1)) / 2
      f = g
      h = g + n - 2
      m = 1
      l = 2
      j = 0
      k = 1
20    if ( k .gt. n ) go to 40
      p = l - 1
      a(p) = a(p) - e
      if ( k .eq. n ) l = p
      do i = m, l
        a(i-j) = a(i)
      end do
      j = k
      k = k + 1
      a(h+k) = a(l)
      m = l + 1
      l = m + k
      go to 20
40    i = g + 1
      j = h + 2
      k = h + n
50    if ( i .gt. k ) go to 60
      a(i) = a(j)
      i = i + 2
      j = j + 1
      go to 50
60    k = n
      m = g
      l = g - n + 1
      g = h + n - 1
      r = cabs ( a(m) ) + cabs ( a(g) ) 
70    if ( k .eq. 1 ) go to 140
      j = k
      k = k - 1
      y = a(m)
      z = a(g)
      s = cabs ( y )
      t = cabs ( z )
      if ( s .ge. t ) go to 100
      a(m) = z
      if ( r .lt. t ) go to 80
      r = t
      p = j
      o = m
80    z = y / z
      a(g) = z
      a(g+1) = (1.0E+00,0.)
      g = g - 2
      l = l - k
      m = m - j

      do i = l, m
        q = i + k
        y = a(i)
        a(i) = a(q) - y*z
        a(q) = y
      end do

      go to 70
100   if ( r .lt. t ) go to 110
      r = t
      p = j
      o = m
110   if ( s .eq. 0. ) go to 130
      z = z / y
      l = l - k
      m = m - j
      do i = l, m
        a(i) = a(i) - z*a(i+k)
      end do
130   a(g) = z
      a(g+1) = (0.,0.)
      g = g - 2
      go to 70
140   t = cabs ( a(1) )
      if ( t .gt. r ) go to 150
      r = t
      p = k
150   do i = 1, n
        x(i) = (0.,0.)
      end do
      z = (1.0E+00,0.)
      l = f
      j = o - p
      k = p
      go to 180
170   z = x(k) / a(j+k)
180   x(k) = z
      if ( k .eq. 1 ) go to 200
      k = k - 1
      do i = 1, k
        x(i) = x(i) - z*a(i+j)
      end do
      j = j - k
      go to 170
200   if ( k .eq. n ) go to 210
      j = k
      k = k + 1
      l = l + 2
      x(k) = x(k) - a(l-1)*x(j)
      if ( real ( a(l) ) .eq. 0. ) go to 200
      z = x(k)
      x(k) = x(j)
      x(j) = z
      go to 200
210   if ( r .eq. 0. ) go to 310

      s = 0.0E+00
      do i = 1, n
        t = amag ( x(i) )
        if ( t .gt. s ) s = t
      end do

      p = h + n
      r = 0.
      s = 1.0E+00 / s
      do i = 1, n
        z = s*x(i)
        r = r + sq(z)
        x(i) = z
        a(i+p) = z
      end do

      y = ( 0.0E+00, 0.0E+00 )
      l = f
      j = l - n
      k = n
240   z = x(k) / a(j+k)
      x(k) = z
      if ( k .eq. 1 ) go to 260
      k = k - 1
      do i = 1, k
        x(i) = x(i) - z*a(i+j)
      end do
      j = j - k
      go to 240
260   t = amag ( x(1) )
270   if ( k .eq. n ) go to 290
      j = k
      k = k + 1
      l = l + 2
      z = x(j)
      x(k) = x(k) - a(l-1)*z
      s = amag ( x(k) )
      if ( s .gt. t ) t = s
      if ( real ( a(l) ) .eq. 1.0E+00 ) go to 280
      y = y + conjg ( a(p+j) ) * z
      go to 270
280   x(j) = x(k)
      x(k) = z
      y = y + conjg ( a(p+j) ) * x(j)
      go to 270
290   y = y + conjg ( a(p+n) ) * x(n)
      if ( amag ( y ) .ne. 0.0E+00 ) y = r / y
      s = 0.0E+00
      t = 1.0E+00 / t
      do i = 1, n
        s = s + sq(y*x(i))
        x(i) = t*x(i)
      end do
      if ( r+r .ge. s ) e = e + y
      return
310   t = 0.0E+00
      do i = 1, n
        s = amag ( x(i) )
        if ( s .gt. t ) t = s
      end do

      t = 1.0E+00 / t

      do i = 1, n
        x(i) = t*x(i)
      end do

      return
      end
      subroutine chess ( a, la, n, w )

c*********************************************************************72
c
cc CHESS reduces a complex matrix to Hessenberg form.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
C     |    INPUT:                                              |
C     |                                                        |
C     |         A     --COMPLEX ARRAY CONTAINING MATRIX        |
C     |                 (LENGTH AT LEAST 2 + N(N+1))           |
C     |                                                        |
C     |         LA    --LEADING (ROW) DIMENSION OF ARRAY A     |
C     |                                                        |
C     |         N     --DIMENSION OF MATRIX STORED IN A        |
C     |                                                        |
C     |         W     --WORK ARRAY WITH AT LEAST N COMPLEX     |
C     |                 ELEMENTS                               |
C     |                                                        |
C     |    OUTPUT:                                             |
C     |                                                        |
C     |         A     --HESSENBERG MATRIX  
c
      implicit none

      integer n

      complex a(*)
      real amag
      integer c
      integer d
      integer e
      integer g
      integer h
      integer i
      integer j
      integer k
      integer l
      integer la
      integer m
      integer o
      integer p
      integer q
      real r
      real s
      real sqr
      real t
      real u
      complex w(n)
      complex x
      complex y
      complex z

      if ( la .gt. n ) then
        call cpack(a,la,n)
      end if

      i = n*n
      u = 0.
      m = n + 1
      o = m + 1
      e = i + m
      j = m
      k = i
10    k = k - n
20    a(i+j) = a(i)
      i = i - 1
      if ( i .gt. k ) go to 20
      j = j - 1
      if ( k .gt. 0 ) go to 10
      a(1) = 2232
      a(2) = n
      k = 4
      l = o
      d = 1
      c = 2
30    if ( c .ge. n ) return
      p = k + 1
      do i = p, l
        if ( amag ( a(i) ) .ne. 0.0E+00 ) go to 50
      end do
      a(l+1) = (0.,0.)
      go to 180
50    t = amag ( a(k) )
      if ( t .ne. 0.0E+00 ) u = 1.0E+00 / t
      r = sqr(a(k),u)

      do j = i, l
        s = amag ( a(j) )
        if ( s .gt. t ) then
          u = 1.0E+00 / s
          r = sqr(a(j),u) + r*(t*u)**2
          t = s
        else
          r = r + sqr(a(j),u)
        end if
      end do

      s = t * sqrt ( r )
      z = a(k)
      t = cabs ( z )
      u = 1.0E+00 / sqrt ( s * ( s + t ) )
      if ( t .ne. 0.0E+00 ) z = z / t
      if ( t .eq. 0.0E+00 ) z = (1.0E+00,0.)
      i = l
80    a(i+1) = u*conjg ( a(i) )
      i = i - 1
      if ( i .gt. k ) go to 80
      a(k) = -z*s
      a(p) = conjg ( z ) * u * ( t + s )
      h = l

      do i = 1, n
        w(i) = (0.,0.)
      end do

100   h = h + m
      y = conjg ( a(p) )
      p = p + 1
      q = h - n
      do i = 1, d
        w(i) = w(i) + y*a(i+q)
      end do
      j = k - d
      z = (0.0E+00,0.0E+00)
      do i = c, n
        x = a(i+q)
        z = z + x*a(i+j)
        w(i) = w(i) + x*y
      end do

      a(h+1) = z
      if ( h .lt. e ) go to 100
      z = (0.,0.)
      h = l + 1
      p = k + 1
      j = c - p
      do i = p, h
        z = z + w(i+j)*a(i)
        a(i) = conjg ( a(i) )
      end do
      do i = c, n
        w(i) = w(i) - z*a(i-j)
      end do
      h = l
150   g = h + 2
      q = h + m
      h = h + c
      z = a(q+1)
      y = conjg ( a(p) )
      p = p + 1
      j = 1 - g
      do i = g, h
        a(i) = a(i) - w(i+j)*y
      end do
      i = h
      h = q
      q = k - i
      g = i + 1
      do i = g, h
        a(i) = a(i) - a(i+q)*z - w(i+j)*y
      end do
      if ( h .lt. e ) go to 150
180   k = k + o
      l = l + m
      d = c
      c = c + 1
      go to 30
      end
      function con ( a, b )

c*********************************************************************72
c
cc CON estimates the 1-norm condition number of a general matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c    William Hager,
c    Condition Estimates,
c    SIAM Journal on Scientific and Statistical Computing,
c    Volume 5, Number 2, June 1984, pages 311-316.
c
c  Parameters:
c
c    Input, real A(*), the factorization information from FACT.
c
c    Workspace, real B(N).
c
c    Output, real CON, the condition number of the matrix.
c
      implicit none

      real a(*)
      real b(*)
      real c
      real con
      real d
      integer i
      integer j
      integer m
      integer n

      c = a(1)

      if ( abs ( c ) .ne. 1230 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CON - Fatal error!'
        write ( *, '(a)' ) '  Input matrix must be factored first!'
        stop
      end if

      if ( c .lt. 0.0E+00 ) then
        con = -1.0E+00
        return
      end if

      m = 0
      n = int ( a(2) )
      c = 1.0E+00 / a(2)

      do j = 1, n
        b(j) = c
      end do

10    continue

      call solve ( b, a, b )

      c = 0.0E+00
      do j = 1, n
        c = c + abs ( b(j) )
        if ( b(j) .lt. 0.0E+00 ) then
          b(j) = -1.0E+00
        else
          b(j) = 1.0E+00
        end if
      end do

      call trans ( b, a, b )
c
c  I is the index of the largest entry of B.
c
      i = 1
      do j = 1, n
        if ( abs ( b(i) ) .lt. abs ( b(j) ) ) then
          i = j
        end if
      end do

      if ( 0 .lt. m ) then
        if ( m .eq. i .or. c .le. d ) then
          c = c * a(3)
          c = max ( c, 1.0E+00 )
          go to 20
        end if
      end if

      m = i
      d = c

      do j = 1, n
        b(j) = 0.0E+00
      end do
      b(m) = 1.0E+00

      go to 10

20    continue

      con = c

      end
      subroutine cong ( x, e, it, step, t, limit, n, m, value, grad, 
     &  both, pre, h )

c*********************************************************************72
c
cc CONG minimizes a function using the conjugate gradient method.
c 
c  Discussion:
c                                                             
c    This routine minimizes a function using the Fletcher-Reeves form   
c    of the conjugate gradient method with optional preconditioning.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input/output, real X(N).  On input, a starting estimate, on
c    output, the computed minimizer.
c
c    Output, real E, the maximum-norm of the gradient vector at X.
c
c    Output, integer IT, the number of iterations performed.
c
c    Output, real STEP, the final stepsize.
c            
c         input:                                              
c                                                                                                                      
c              step  --starting guess for minimizer in direc- 
c                      tion of negative gradient during first 
c                      iteration (e. g. step=1) when step=0,  
c                      the program selects a starting guess   
c                                                             
c              t     --computing tolerance (iterations stop   
c                      when max-norm of gradient .le. t)      
c                                                             
c              limit --maximum number of iterations           
c                                                             
c              n     --number of unknowns                     
c                                                             
c              m     --number of iterations until the search  
c                      directions are renormalized along the  
c                      negative gradient (typically, m = n)   
c                                                             
c              value --name of cost evaluation func. routine  
c                      (external in main program)             
c                      value(x) is value of cost at x         
c                                                             
c              grad  --name of gradient evaluation routine 
c                      (external in main program)             
c                      grad(g,x) puts in g the gradient at x  
c                                                             
c              both  --name of routine to evaluate both cost  
c                      and its gradient (external in main     
c                      program) both(v,g,x) puts the value in 
c                      v and the gradient in g for the point x
c                                                             
c              pre   --name of preconditioning routine     
c                      (external in main program)             
c                      pre(y,z) applies the preconditioner to 
c                     z, storing the result in y.            
c                      if preconditioning not used set y = z  
c                                                             
c              h     --work array (length at least 3n)        
c                                                                                                                       
      implicit none

      integer n

      real a
      real a1
      real a2
      real a3
      real a4
      real a5
      real a6
      real a7
      real a8
      real b
      external both
      real c
      real c0
      real c1
      real d
      real d0
      real da
      real db
      real e
      real f
      real f0
      real f1
      real fa
      real fb
      real fc
      real fd
      real fv
      real g
      external grad
      real h(n,*)
      integer i
      integer iq
      integer it
      integer j
      integer k
      integer l
      integer limit
      real l3
      integer m
      integer na
      integer nb
      integer nc
      integer nd
      real p
      external pre
      real q
      real r
      real s
      real step
      real t
      real v
      real value
      external value
      real w
      real x(n)
      real y(50)
      real z(50)

      save a1
      save a2
      save a3
      save a4
      save a5
      save a6
      save a7

      data a1 /0.10/
      data a2 /0.90/
      data a3 /5.0/
      data a4 /0.20/
      data a5 /10.0/
      data a6 /0.90/
      data a7 /0.3/

      a8 = a3 + 0.01E+00
      it = 0
      call both ( f, h(1,3), x )
      e = 0.
      do i = 1, n
        if ( abs ( h(i,3) ) .gt. e ) then
          e = abs ( h(i,3) )
        end if
      end do
      if ( e .le. t ) return
      l3 = 1.0E+00 / alog ( a3 )
      call pre(h(1,2),h(1,3))
      a = step
      if ( a .gt. 0.0E+00 ) go to 30
      do i = 1, n
        if ( abs ( x(i) ) .gt. a ) a = abs ( x(i) )
      end do
      a = 0.01E+00 * a / e

      if ( a .eq. 0.0E+00 ) then
        a = 1.0E+00
      end if

30    g = 0.
      do i = 1, n
        g = g + h(i,2)*h(i,3)
      end do
      if ( g .lt. 0.0E+00 ) then
        go to 620
      end if
50    l = 0
      do i = 1, n
        h(i,1) = -h(i,2)
      end do
      d = -g
70    fa = fv(a,x,h,n,value)
      c0 = a
      f0 = fa
      j = 2
      y(1) = 0.0E+00
      z(1) = f
      y(2) = a
      z(2) = fa
      v = a1*d
      w = a2*d
      iq = 0
      if ( fa .le. f ) go to 80
      c = a
      b = 0.
      a = 0.
      fc = fa
      fb = f
      fa = f
      go to 90
80    c = 0.
      b = 0.
      fc = f
      fb = f
      iq = 1
90    na = 0
      nb = 0
      nc = 0
      nd = 0
      q = (d+(f-f0) / c0 ) / c0
      if ( q .lt. 0. ) go to 110
      q = a
100   nd = nd + 1
      if ( nd .gt. 25 ) go to 610
      q = a3*q
      p = fv(q,x,h,n,value)
      call ins(q,p,a,b,c,fa,fb,fc,j,y,z)
      if ( p-f .lt. w*q ) go to 100
      go to 260
110   q = .5*d / q
      if ( q .lt. .01*c0 ) q = .01*c0
      p = fv(q,x,h,n,value)
      if ( p .le. f0 ) go to 120
      f1 = f0
      c1 = c0
      f0 = p
      c0 = q
      go to 130
120   f1 = p
      c1 = q
130   call ins(q,p,a,b,c,fa,fb,fc,j,y,z)
      if ( a .eq. 0. ) go to 140
      if ( fa-f .ge. v*a ) go to 160
      if ( fa-f .lt. w*a ) go to 210
      go to 280
140   q = c0
      if ( c1 .lt. q ) q = c1
150   na = na + 1
      if ( na .gt. 25 ) go to 630
      q = a4*q
      p = fv(q,x,h,n,value)
      call ins(q,p,a,b,c,fa,fb,fc,j,y,z)
      if ( p-f .ge. v*q ) go to 150
      go to 250
160   if ( c0 .gt. c1 ) go to 200
      if ( f0-f .gt. v*c0 ) go to 180
      if ( f0-f .ge. w*c0 ) go to 320
      if ( c1 .le. a5*c0 ) go to 320
      r = alog (c1 / c0 )
      s = - int ( r * l3 + 0.999 )
      r = 0.999 * exp ( r / s )
      q = c1
170   q = q*r
      if ( q .lt. c0 ) go to 320
      p = fv(q,x,h,n,value)
      call ins(q,p,a,b,c,fa,fb,fc,j,y,z)
      na = na + 1
      if ( p-f .gt. v*q ) go to 170
      go to 320
180   q = c0
190   na = na + 1
      if ( na .gt. 25 ) go to 630
      q = a4*q
      p = fv(q,x,h,n,value)
      call ins(q,p,a,b,c,fa,fb,fc,j,y,z)
      if ( p-f .ge. v*q ) go to 190
      go to 250
200   q = a
      go to 190
210   if ( c0 .lt. c1 ) go to 290
      if ( f0-f .ge. v*c0 ) go to 230
      if ( f0-f .ge. w*c0 ) go to 250
      q = c0
220   nd = nd  + 1
      if ( nd .gt. 25 ) go to 610
      q = a3*q
      p = fv(q,x,h,n,value)
      call ins(q,p,a,b,c,fa,fb,fc,j,y,z)
      if ( p-f .lt. w*q ) go to 220
      go to 250
230   if ( c0 .le. a5*c1 ) go to 250
      r = alog ( c0 / c1 )
      s = int(r*l3+0.999)
      r = 1.001E+00 * exp ( r / s )
      q = a
240   q = q*r
      if ( q .gt. c0 ) go to 250
      nd = nd + 1
      p = fv(q,x,h,n,value)
      call ins(q,p,a,b,c,fa,fb,fc,j,y,z)
      if ( p-f .lt. w*q ) go to 240
250   if ( iq .eq. 1 ) go to 320
260   if ( b .eq. 0. ) go to 280
      if ( c .eq. 0. ) go to 270
      v = c - a
      w = a - b
      r = 1.0E+00 / v
      s = 1.0E+00 / w
      p = fc - fa
      q = fb - fa
      e = p*r + q*s
      if ( sign(e,c-b) .ne. e ) go to 320
      if ( e .eq. 0. ) go to 320
      q = (p*r)*w - (q*s)*v
      q = a - .5*q / e
      p = fv(q,x,h,n,value)
      call ins(q,p,a,b,c,fa,fb,fc,j,y,z)
      go to 320
270   r = 1.0E+00 / a
      s = 1.0E+00 / b
      p = r*(fa-f) - d
      q = s*(fb-f) - d
      e = a - b
      v = (r*p-s*q) / e
      w = (a*q*s-b*p*r) / e
      v = w*w-3.*v*d
      if ( v .lt. 0. ) v = 0.
      v = sqrt ( v )
      if ( w+v .eq. 0. ) go to 320
      q = -d / (w+v)
      if ( q .le. 0. ) go to 320
      p = fv(q,x,h,n,value)
      call ins(q,p,a,b,c,fa,fb,fc,j,y,z)
      go to 320
280   if ( iq .eq. 1 ) go to  320
      q = (d+(f-fa) / a ) / a
      if ( q .ge. 0. ) go to 320
      q = .5*d / q
      p = fv(q,x,h,n,value)
      call ins(q,p,a,b,c,fa,fb,fc,j,y,z)
      go to 320
290   if ( f0-f .gt. v*c0 ) go to 300
      if ( f0-f .gt. w*c0 ) go to 320
300   q = a
310   nd = nd + 1
      if ( nd .gt. 25 ) go to 610
      q = a3*q
      p = fv(q,x,h,n,value)
      call ins(q,p,a,b,c,fa,fb,fc,j,y,z)
      if ( p-f .lt. w*q ) go to 310
      go to 250
320   da = fd ( a, x, h, n, grad )
      if ( da .gt. a6*g ) go to 410
      if ( da .ge. 0. ) go to 560
      r = a
      q = 0.
      do 330 i = 1, j
        if ( y(i) .gt. a ) go to 370
        if ( y(i) .le. q ) go to 330
        if ( y(i) .eq. a ) go to 330
        q = y(i)
330   continue
      if ( a .le. a8*q ) go to 560
      q = a
340   nd = nd + 1
      if ( nd .gt. 25 ) go to 610
      q = a3*q
      p = fv(q,x,h,n,value)
      f1 = fa
      call ins(q,p,a,b,c,fa,fb,fc,j,y,z)
      if ( p .lt. f1 ) go to 340
      if ( a .gt. r ) go to 360
      do i = 1, n
        h(i,2) = x(i) + a*h(i,1)
      end do
      go to 560
360   da = fd(a,x,h,n,grad)
      if ( da .gt. a6*g ) go to 410
      go to 560
370   q = y(i)
      do 380 k = i, j
        if ( y(k) .le. a ) go to 380
        if ( y(k) .lt. q ) q = y(k)
380   continue
      if ( q .le. a5*a ) go to 560
      f0 = alog ( q / a )
      s = int(f0*l3+0.999)
      f0 = 1.001E+00 * exp ( f0 / s )
      s = a
390   s = s*f0
      if ( s .ge. q ) go to 320
      p = fv(s,x,h,n,value)
      f1 = fa
      call ins(s,p,a,b,c,fa,fb,fc,j,y,z)
      if ( p .lt. f1 ) go to 390
      if ( a .gt. r ) go to 320
      do i = 1, n
        h(i,2) = x(i) + a*h(i,1)
      end do
      go to 560
410   b = 0.
      k = 1
      i = k
420   i = i + 1
      if ( i .gt. j ) go to 430
      if ( y(i) .ge. a ) go to 420
      if ( y(i) .lt. b ) go to 420
      b = y(i)
      k = i
      go to 420
430   fb = z(k)
      db = d
      if ( b .ne. 0. ) db = fd(b,x,h,n,grad)
440   w = 2.*abs ( b - a )
      call cub(c,a,b,fa,fb,da,db)
      nc = 1
      go to 480
450   w = .5*w
      if ( w .lt. abs ( c0 - c ) ) go to 550
      if ( c0 .lt. c ) go to 460
      if ( d0 .ge. d ) go to 470
      go to 550
460   if ( d0 .gt. d ) go to 550
470   call cub(c,c,c0,f,f0,d,d0)
      nc = nc + 1
      if ( nc .gt. 30 ) go to 600
480   r = max ( a, b )
      s = min ( a, b )
      if ( c .gt. r ) go to 490
      if ( c .gt. s ) go to 500
      c = s + (s-c)
      s = .5*(a+b)
      if ( c .gt. s ) c = s
      go to 500
490   c = r - (c-r)
      s = .5*(a+b)
      if ( c .lt. s ) c = s
500   c0 = a
      f0 = fa
      d0 = da
      call fvd(f,d,c,x,h,n,both)
      if ( f .lt. fa ) go to 510
      b = c
      fb = f
      db = d
      go to 450
510   if ( c .lt. a ) go to 540
      if ( d .lt. 0. ) go to 530
520   b = a
      fb = fa
      db = da
530   a = c
      fa = f
      da = d
      if ( d .gt. a6*g ) go to 450
      go to 560
540   if ( d .lt. 0. ) go to 520
      go to 530
550   c = .5*(a+b)
      nb = nb + 1
      w = abs ( b - a )
      go to 500
560   e = 0.
      do i = 1, n
        if ( abs ( h(i,3) ) .gt. e ) e = abs ( h(i,3) )
        x(i) = h(i,2)
      end do
      it = it + 1
      if ( e .le. t ) go to 660
      if ( it .ge. limit ) go to 660
      f = fa
      d = da
      a = a7*a
      call pre(h(1,2),h(1,3))
      r = 0.
      do i = 1, n
        r = r + h(i,2)*h(i,3)
      end do
      if ( r .lt. 0. ) go to 620
      s = r / g
      g = r
      l = l + 1
      if ( l .ge. m ) go to 50
      d = 0.
      do i = 1, n
        h(i,1) = -h(i,2) + s*h(i,1)
        d = d + h(i,1)*h(i,3)
      end do
      go to 70
600   if ( d .lt. g ) go to 560
      write(*,*) 'unable to obtain descent direction'
      stop
610   write(*,*) 'the function decreases with no minimum'
      stop
620   write(*,*) 'preconditioner not positive definite'
      stop
630   q = q*a3**25
      nd = 0
640   nd = nd + 1
      if ( nd .gt. 25 ) go to 650
      q = a3*q
      p = fv(q,x,h,n,value)
      call ins(q,p,a,b,c,fa,fb,fc,j,y,z)
      if ( p-f .gt. v*q ) go to 640
      go to 130
650   write(*,*) 'unable to satisfy armijo condition'
      stop
660   step = a
      return
      end
      subroutine cp ( t, b, v, l, d, p, n )

c*********************************************************************72
c
cc CP
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real a
      real b
      real c
      real d(*)
      real e
      real f
      integer i
      integer j
      integer l
      integer m
      integer n
      real p(*)
      real t
      real u
      real v
      real z

      l = 0
      z = 0.0E+00
      v = z
      if ( n .gt. 1 ) go to 10
      b = d(1) - t
      if ( b .le. z ) l = 1
      return
10    f = 65536.0E+00**4
      e = 1.0E+00 / f
      u = 16.
      m = 0
      i = 1
20    c = 1.0E+00
      b = d(i) - t
30    if ( b .le. z ) go to 70
      if ( i .ge. n ) go to 130
40    j = i
      i = i + 1
      a = (d(i)-t)*b - p(j)*c
      c = b
      b = a
50    a = abs ( b )
      if ( a .gt. f ) go to 60
      if ( a .gt. e ) go to 30
      if ( a .eq. z ) go to 70
      c = c*f
      b = b*f
      v = v - u
      go to 50
60    c = c*e
      b = b*e
      v = v + u
      go to 50
70    l = l + 1
      if ( i .ge. n ) go to 130
      if ( b .lt. z ) go to 90
      if ( p(i) .gt. 0. ) go to 90
      i = i + 1
      m = 1
      v = z
      go to 20
80    if ( b .ge. z ) go to 120
      if ( i .ge. n ) go to 130
90    j = i
      i = i + 1
      a = (d(i)-t)*b - p(j)*c
      c = b
      b = a
100   a = abs ( b )
      if ( a .gt. f ) go to 110
      if ( a .gt. e ) go to 80
      if ( a .eq. z ) go to 120
      c = c*f
      b = b*f
      v = v - u
      go to 100
110   c = c*e
      b = b*e
      v = v + u
      go to 100
120   l = l + 1
      if ( i .ge. n ) go to 130
      if ( b .gt. z ) go to 40
      if ( p(i) .gt. 0. ) go to 40
      i = i + 1
      m = 1
      v = z
      go to 20
130   if ( m .eq. 1 ) go to 160
      if ( b .eq. 0. ) go to 160
      a = 1.0E+00 / u
      if ( abs ( b ) .lt. a ) go to 150
140   if ( abs ( b ) .lt. 1.0E+00 ) return
      b = b*a
      v = v + 1.0E+00
      go to 140
150   b = b*u
      v = v - 1.0E+00
      if ( abs ( b ) .lt. a ) go to 150
      return
160   v = z
      b = z
      return
      end
      subroutine cpack ( a, la, n )

c*********************************************************************72
c
cc CPACK packs the elements of a complex square array.
c
c  Discussion:
c
c    The packing is done in such a way that the elements of the matrix
c    are stored sequentially.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input/output, complex A(*).
c    On input A is assumed to be an NxN matrix stored as an LAxN array.
c    On output, A has been packed as an NxN array.
c
c    Input, integer LA, the leading dimension of the array.
c    N <= LA.
c
c    Input, integer N, the number of rows and columns in the matrix.
c
      implicit none

      complex a(*)
      integer h
      integer i
      integer j
      integer jhi
      integer jlo
      integer k
      integer l
      integer la
      integer n
      integer o

      h = la - n

      if ( h .eq. 0 ) then
        return
      end if

      if ( h .lt. 0 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: la argument in qpack must be .ge. n argument'
        stop
      end if

      i = 0
      k = 1
      l = n
      o = n * n

10    continue

      if ( l .eq. o ) then
        jlo = n * n + 1
        jhi = la * n
        do j = jlo, jhi
          a(j) = cmplx ( 0.0E+00, 0.0E+00 )
        end do
        return
      end if

      i = i + h
      k = k + n
      l = l + n
      do j = k,l
        a(j) = a(i+j)
      end do

      go to 10

      end
      subroutine cpower ( e, x, n, dif, size, ndigit, limit, mult, w )

c*********************************************************************72
c
cc CPOWER uses the power method for a complex matrix.
c
c  Discussion:
c
c    This routine computes an eigenvalue of largest modulus, and
c    the corresponding eigenvector, using the power method.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, complex E, the computed eigenvalue.
C     |                                                        |
C     |    INPUT:                                              |
C     |                                                        |
C     |         X     --COMPLEX ARRAY CONTAINING STARTING GUESS|
C     |                 FOR THE EIGENVECTOR SPACE              |
C     |                                                        |
C     |         N     --MATRIX DIMENSION                       |
C     |                                                        |
C     |         NDIGIT--DESIRED NUMBER OF CORRECT DIGITS       |
C     |                                                        |
C     |         LIMIT --MAXIMUM NUMBER OF ITERATIONS           |
C     |                                                        |
C     |         MULT  --NAME OF SUBROUTINE TO MULTIPLY MATRIX A|
C     |                 BY VECTOR (NAME EXTERNAL IN MAIN PROG.)|
C     |                 MULT(P,V) STORES IN P THE PRODUCT AV   |
C     |                 WHERE P AND V ARE COMPLEX ARRAYS       |
C     |                                                        |
C     |         W     --COMPLEX WORK ARRAY WITH N ELEMENTS     |
C     |                                                        |
C     |    OUTPUT:                                             |
C     |                                                        |
C     |         X     --EIGENVECTOR                            |
C     |                                                        |                             |
C     |         DIF   --INPUT FOR SUBROUTINE WHATIS            |
C     |                                                        |
C     |         SIZE  --INPUT FOR SUBROUTINE WHATIS            |
c
      implicit none

      real amag
      real dif
      complex e
      integer i
      integer j
      integer l
      integer limit
      integer n
      integer ndigit
      external mult
      real r
      real s
      real size
      real sq
      real t
      complex u
      complex v
      complex w(*)
      complex x(*)
      complex y

      t = 0.

      do i = 1, n
        s = amag ( x(i) )
        if ( s .gt. t ) t = s
      end do

      if ( t .eq. 0. ) then
        write(*,*) 'starting guess for routine cpower cannot be zero'
        stop
      end if

      t = 1.0E+00 / t
      do i = 1, n
        w(i) = t*x(i)
      end do
      r = 10.0E+00**(-ndigit)
      l = 0
40    l = l + 1
      do i = 1, n
        x(i) = w(i)
      end do
      call mult ( w, x )

      t = 0.
      do i = 1, n
        s = amag ( w(i) )
        if ( s .gt. t ) then
          t = s
          j = i
        end if
      end do

      if ( t .eq. 0. ) go to 110
      v = 1.0E+00 / w(j)
      u = (0.,0.)
      if ( amag ( x(j) ) .ne. 0. ) u = 1.0E+00 / x(j)

      dif = 0.
      do i = 1, n
        y = x(i)*u
        w(i) = w(i)*v
        dif = dif + amag ( w(i) - y )
      end do

      if ( dif .le. r ) go to 100
      if ( l .lt. limit ) go to 40
      size = 3*l + 2
80    e = (0.,0.)
      t = 0.0E+00
      do i = 1, n
        e = e + conjg ( x(i) ) * w(i)
        t = t + sq(x(i))
        x(i) = w(i)
      end do
      e = e / ( t * v )
      size = 3*l + 2
      return
100   size = 3*l + 1
      go to 80
110   e = (0.,0.)
      size = 3*l + 1
      do i = 1, n
        x(i) = w(i)
      end do

      return
      end
      subroutine csim ( p, lp, a )

c*********************************************************************72
c
cc CSIM computes the Hessenberg similarity transform for a complex matrix.
c
c  Discussion:
c
C     |COMPUTE THE SIMILARITY TRANSFORMATION USED IN REDUCTION |
C     |TO HESSENBERG FORM: ORIGINAL A = P TIMES HESSENBERG     |
C     |                MATRIX TIMES P INVERSE                  |
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c
C     |                                                        |
C     |    INPUT:                                              |
C     |                                                        |
C     |         A     --OUTPUT FROM CHESS OR CAHESS            |
C     |                                                        |
C     |         LP    --LEADING (ROW) DIMENSION OF ARRAY P     |
C     |                                                        |
C     |    OUTPUT:                                             |
C     |                                                        |
C     |         P     --SIMILARITY TRANSFORMATION              |
C     |                 (COMPLEX ARRAY)      
c
      implicit none

      integer lp

      complex a(*)
      integer h
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      complex p(lp,*)
      real t
      complex z

      t = a(1)

      if ( t .ne. 2232 .and. t .ne. 2233 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must reduce matrix to hessenberg form using'
        write(*,*) 'either routine chess or cahess before using'
        write(*,*) 'routine csim to compute the similarity'
        write(*,*) 'transformation'
        stop
      end if

      n = a(2)
      p(1,1) = (1.0E+00,0.)
      if ( n .eq. 1 ) return
      p(1,2) = (0.,0.)
      if ( n .gt. 2 ) go to 20
      p(2,1) = (0.,0.)
      p(2,2) = (1.0E+00,0.)
      go to 120
20    l = n - 2
      m = n + 1
      h = 3
      do j = 1, l
        k = j + 1
        do i = k, n
          p(i,j) = a(i+h)
        end do
        h = h + m
      end do
      k = n
      m = n - 1
      z = conjg ( p(k,l) )
      p(n,n) = (1.0E+00,0.) - z*p(n,l)
      p(m,n) = -z*p(m,l)
60    j = m
      m = l
      l = l - 1
      if ( l .eq. 0 ) go to 90
      z = (0.,0.)
      do i = j, n
        z = z + conjg ( p(i,l) ) * p(i,k)
      end do
      p(m,k) = -z*p(m,l)
      do i = j, n
        p(i,k) = p(i,k) - z*p(i,l)
      end do
      go to 60
90    p(1,k) = (0.,0.)
      k = k - 1
      m = k
      l = k - 1
      z = -conjg ( p(k,l) )

      do i = k, n
        p(i,k) = z*p(i,l)
      end do

      p(k,k) = (1.0E+00,0.) + p(k,k)
      if ( l .gt. 1 ) go to 60

      do i = 2, n
        p(i,1) = (0.,0.)
      end do

120   if ( t .eq. 2232 ) return
      k = 1 + n + n*n

      do j = 1, n
        do i = 1, n
          p(i,j) = a(i+k)*p(i,j)
        end do
      end do

      return
      end
      subroutine csolve ( x, a, b )

c*********************************************************************72
c
cc CSOLVE solves a general complex factored system.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, complex X(N), the solution.
c
c    Input, complex A(*), the factor information from CFACT.
c
c    Input, complex B(N), the right hand side.
c
      implicit none

      complex a(*)
      complex b(*)
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      complex t
      complex x(*)

      i = int ( a(1) )

      if ( abs ( i ) .ne. 1239 ) then
        write ( *, '(a)' ) ' '
        write(6,*) 'error: must factor before solving'
        stop
      end if
c
c  forward elimination.
c
      n = int(a(2))
      m = n + 1
      j = 4 - m
      if ( i .lt. 0. ) go to 80
      do i = 1, n
        x(i) = b(i)
      end do
      k = 1
30    j = j + m
      if ( cabs ( a(j+k) ) .eq. 0. ) go to 80
      if ( k .eq. n ) go to 50
      l = int(a(j))
      t = x(l)
      x(l) = x(k)
      x(k) = t
      k = k + 1
      if ( cabs ( t ) .eq. 0. ) go to 30
      do i = k, n
        x(i) = x(i) - t*a(i+j)
      end do
      go to 30
c
c back substitution by columns.
c
50    t = x(k) / a(j+k)
60    x(k) = t
      if ( k .eq. 1 ) return
      k = k - 1
      do i = 1, k
        x(i) = x(i) - t*a(i+j)
      end do
      j = j - m
      go to 50
c
c  compute null vector.
c
80    k = 0
90    k = k + 1
      j = j + m
      if ( cabs ( a(j+k) ) .ne. 0. ) go to 90
      do i = 1, n
        x(i) = 0.0E+00
      end do
      t = 1.0E+00
      go to 60
      end
      subroutine ctrans ( x, a, b )

c*********************************************************************72
c
cc CTRANS solves a transposed general complex factored system.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, complex X(N), the solution.
c
c    Input, complex A(*), factorization information from CFACT.
c
c    Input, complex B(N), the right hand side.
c                                                       |
      implicit none

      complex a(*)
      complex b(*)
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      complex t
      complex x(*)

      i = int ( a(1) )

      if ( abs ( i ) .ne. 1239 ) then
        write ( *, '(a)' ) ' '
        write(6,*) 'error: must factor with cfact before solving'
        stop
      end if

      n = int(a(2))
      m = n + 1
      if ( i .lt. 0. ) go to 80
      t = 0.
      j = 4
      k = 1
c
c  skip over zeros.
c
20    if ( b(k) .ne. 0.0E+00 ) go to 30
      x(k) = 0.
      k = k + 1
      if ( k .le. n ) go to 20
      return
c
c  Forward substitution.
c
30    j = j - m + m*k
40    x(k) = (b(k)-t) / a(j+k)
      if ( k .eq. n ) go to 60
      t = 0.
      j = j + m
      do i = 1, k
        t = t + a(i+j)*x(i)
      end do
      k = k + 1
      go to 40
c
c  back substitution.
c
60    if ( k .eq. 1 ) return
      j = j - m
      t = x(k-1)
      do i = k, n
        t = t - x(i)*a(i+j)
      end do
      k = k - 1
      i = int(a(j))
      x(k) = x(i)
      x(i) = t
      go to 60
c
c  compute null vector
c
80    i = 5 + n + m*n
      l = m
90    i = i - m - 1
      l = l - 1
      if ( cabs ( a(i) ) .ne. 0. ) go to 90
      k = l
      j = i - k
      do i = 1, n
        x(i) = 0.
      end do
      x(k) = 1.0E+00
110   if ( k .eq. n ) go to 60
      t = 0.
      j = j + m
      do i = l, k
        t = t - a(i+j)*x(i)
      end do
      k = k + 1
      x(k) = t / a(j+k)
      go to 110
      end
      subroutine cub ( x, a, b, c, d, e, f )

c*********************************************************************72
c
cc CUB
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real a
      real b
      real c
      real d
      real e
      real f
      real g
      real v
      real w
      real x
      real y
      real z

      g = b - a

      if ( g .eq. 0.0E+00 ) then
        x = a
        return
      end if

      v = e + f - 3.0E+00 * ( d - c ) / g
      w = v * v - e * f
      w = max ( w, 0.0E+00 )
      w = sign ( sqrt ( w ), g )
      y = e + v
      z = f + v
      if ( sign(y,g) .ne. y ) go to 30
      if ( sign(z,g) .ne. z ) go to 20
      if ( z .eq. 0. ) go to 20
10    x = b - g*f / (z+w)
      return
20    if ( c .lt. d ) x = a
      if ( c .ge. d ) x = b
      return
30    if ( sign(z,g) .ne. z ) go to 40
      if ( abs ( e ) .gt. abs ( f ) ) go to 10
40    x = a + g*e / (y-w)
      return
      end
      subroutine cubic ( x, a, b, c, d, e, f )

c*********************************************************************72
c
cc CUBIC
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real a
      real b
      real c
      real d
      real e
      real f
      real g
      real v
      real w
      real x
      real y
      real z

      g = b - a
      v = e + f - 3 * ( d - c ) / g
      w = v * v - e * f
      w = max ( w, 0.0E+00 )
      w = sign ( sqrt ( w ), g )
      y = e + v
      z = f + v
      if ( sign ( y, g ) .ne. y ) go to 40
      if ( sign ( z, g ) .ne. z ) go to 20
      if ( z .eq. 0.0E+00 ) go to 20
      if ( g * z .le. 0.0E+00 ) go to 20
10    x = b - g * f / ( z + w )
      return
20    x = y + z
      if ( x .eq. 0.0E+00 ) go to 30
      x = a + g * ( y + w ) / ( y + z )
      return
30    x = 0.5E+00 * ( a + b )
      return
40    if ( sign ( z, g ) .ne. z ) go to 50
      if ( abs ( e ) .gt. abs ( f ) ) go to 10
50    x = a + g * e / ( y - w )
      return
      end
      subroutine cunpack ( a, la, n )

c*********************************************************************72
c
cc CUNPACK reverses the operation of CPACK by unpacking the matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input/output, complex A(*).
c    On input, the NxN entries of A were packed into the front of an LAxN array.
c    On output, the NxN entries of A have been unpacked into the LAxN array.
c
c    Input, integer LA, the leading dimension of the array.
c    N <= LA.
c
c    Input, integer N, the number of rows and columns in the matrix.
c
      implicit none

      complex a(*)
      integer i
      integer ii
      integer j
      integer jj
      integer la
      integer n

      if ( la .eq. n ) then
        return
      end if

      if ( la .lt. n ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: la in cunpack must be .ge. n argument!'
        return
      end if

      do jj = 1, n
        j = n + 1 - jj
        do ii = 1, n
          i = n + 1 - ii
          a((j-1)*la+i) = a((j-1)*n+i)
        end do
      end do

      do j = 1, n
        do i = n + 1, la
          a((j-1)*la+i) = cmplx ( 0.0, 0.0 )
        end do
      end do

      return
      end
      subroutine cvals ( e, a, la, n, v )

c*********************************************************************72
c
cc CVALS computes the eigenvalues of a general complex matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
C     |    INPUT:                                              |
C     |                                                        |
C     |         A     --COMPLEX ARRAY CONTAINING MATRIX        |
C     |                 (LENGTH AT LEAST 1 + N(N+2))           |
C     |                                                        |
C     |         LA    --LEADING (ROW) DIMENSION OF ARRAY A     |
C     |                                                        |
C     |         N     --DIMENSION OF MATRIX STORED IN A        |
C     |                                                        |
C     |    OUTPUT:                                             |
C     |                                                        |
C     |         E     --COMPLEX ARRAY OF EIGENVALUES           |
C     |                                                        |
C     |         V     --COMPLEX ARRAY WITH AT LEAST 1+N(N+2)   |
C     |                 ELEMENTS USED AS INPUT FOR SUBROUTINE  |
C     |                 VECT WHEN EIGENVECTORS ARE COMPUTED    |
C     |                 (IF EIGENVECTORS ARE NOT DESIRED,      |
C     |                 ARGUMENT V CAN BE IDENTIFIED WITH A)   |
c
      implicit none

      complex a(*)
      complex e(*)
      integer i
      integer j
      integer k
      integer l
      integer la
      integer m
      integer n
      complex v(*)

      call cahess ( a, la, n, e )
      l = 1 + n*(n+2)

      do i = 1, l
        v(i) = a(i)
      end do

      if ( n .le. 1 ) then
        e(1) = a(3)
        return
      end if

      j = 2
      k = 1
      m = 1
      l = 2
30    do i = m, l
        a(i) = a(i+j)
      end do
      m = l + 1
      j = j + n - k
      k = k + 1
      l = m + k
      if ( k .lt. n ) go to 30
      if ( k .gt. n ) go to 50
      l = l - 1
      go to 30
50    i = (n*(n+3)) / 2
      j = i + n
      k = j + 1 + n / 2
      call vls ( e, a, n, a(k), a(j), a(i) )
      return
      end
      subroutine cvc ( f, g, w, x, y, z, v, k )

c*********************************************************************72
c
cc CVC
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real a
      real b
      real c
      real d
      real f
      real g
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      real p
      real q
      complex r
      complex s
      real v(*)
      real w(*)
      real x(*)
      real y(*)
      real z(*)

      l = k
      j = (k*(k+3)) / 2 - 1
      if ( l .eq. 0 ) return

10    continue

      m = l - 1
      d = f - v(j)
      if ( m .eq. 0 ) go to 40
      b = v(j-l)
      if ( b .eq. 0. ) go to 40
      a = f - v(j-l-1)
      c = v(j-1)
      s = cmplx(a*d-g*g-b*c,g*(a+d))
      if ( cabs ( s ) .ne. 0. ) go to 20
      write(*,*) 'more than two eigenvalues with the same magnitude'
      write(*,*) 'have been detected. routine mpower must stop'
      stop
20    r = cmplx(c*y(l)+d*y(m)-g*z(m),c*z(l)+g*y(m)+d*z(m)) / s
      s = cmplx(b*y(m)+a*y(l)-g*z(l),g*y(l)+a*z(l)+b*z(m)) / s
      a = real ( r )
      b = real ( s )
      c = aimag(r)
      d = aimag(s)
      w(m) = a
      w(l) = b
      x(m) = c
      x(l) = d
      n = j - l
      m = n - l
      j = m - 1
      l = l - 2

      if ( l .eq. 0 ) then
        return
      end if

      do i = 1, l
        p = v(m+i)
        q = v(n+i)
        y(i) = y(i) + a*p + b*q
        z(i) = z(i) + c*p + d*q
      end do

      go to 10
40    r = cmplx(y(l),z(l)) / cmplx(d,g)
      a = real ( r )
      b = aimag(r)
      w(l) = a
      x(l) = b
      n = j - l
      j = n - 1
      l = l - 1
      if ( l .eq. 0 ) return

      do i = 1, l
        y(i) = y(i) + a*v(n+i)
        z(i) = z(i) + b*v(n+i)
      end do

      go to 10
      end
      subroutine cvect ( e, x, v, w )

c*********************************************************************72
c
cc CVECT computes an eigenvector of a general complex matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
C     |    INPUT:                                              |
C     |                                                        |
C     |         E     --ESTIMATE OF EIGENVALUE                 |
C     |                 (COMPLEX VARIABLE)                     |
C     |                                                        |
C     |         V     --OUTPUT FROM EITHER SUBROUTINE CHESS,   |
C     |                 CAHESS, OR CVALS                       |
C     |                                                        |
C     |         W     --COMPLEX WORK ARRAY WITH AT LEAST       |
C     |                 N(N+9)/2 - 3 ELEMENTS                  |
C     |                                                        |
C     |    OUTPUT:                                             |
C     |                                                        |
C     |         E     --IMPROVED ESTIMATE FOR EIGENVALUE       |
C     |                                                        |
C     |         X     --COMPLEX ARRAY CONTAINING EIGENVECTOR   |
c
      implicit none

      real amag
      complex e
      integer h
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      integer o
      real s
      real t
      complex v(*)
      complex w(*)
      complex x(*)
      complex z

      t = v(1)
      if ( t .eq. 2232 ) go to 10
      if ( t .eq. 2233 ) go to 10
        write ( *, '(a)' ) ' '
      write(*,*) 'error: must process the coefficient matrix using'
      write(*,*) 'one of the following routines before computing'
      write(*,*) 'an eigenvector using cvect: chess, cahess, or cvals'
      stop
10    n = v(2)
      if ( n .gt. 1 ) go to 20
      e = v(3)
      x(1) = (1.0E+00,0.)
      return
20    h = (n*(n+3)) / 2
      o = n + 1
      j = 2
      k = 1
      l = 2
      m = 1
30    do i = m, l
        w(i) = v(i+j)
      end do
      m = l + 1
      j = j + n - k
      k = k + 1
      l = m + k
      if ( l .lt. h ) go to 30
      if ( l .gt. h ) go to 50
      l = l - 1
      go to 30
50    call cevect ( e, x, w, n )
      k = n
      j = 3 + (n-2)*o
60    j = j - o
      k = k - 1
      if ( k .le. 1 ) go to 90
      if ( amag ( v(j+k) ) .eq. 0. ) go to 60
      z = (0.,0.)
      do i = k, n
        z = z + x(i)*conjg ( v(i+j) )
      end do
      do i = k, n
        x(i) = x(i) - z*v(i+j)
      end do
      go to 60
90    if ( real ( v(1) ) .eq. 2232 ) go to 110
      j = 1 + n*o
      do i = 1, n
        x(i) = x(i)*v(i+j)
      end do

110   s = 0.
      do i = 1, n
        t = amag ( x(i) )
        if ( t .gt. s ) s = t
      end do
      if ( s .ne. 0. ) s = 1.0E+00 / s 

      do i = 1, n
        x(i) = s*x(i)
      end do

      return
      end
      subroutine cvert ( v, lv, n, w )

c*********************************************************************72
c
cc CVERT inverts a complex matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input/output, complex V(LV,N).  On input, the NxN matrix.
c    On output, the inverse.
c
c    Input, integer LV, the leading dimension of array V.
c
c    Input, integer N, the number of rows and columns of the matrix.
c
c    Workspace, integer W(N-1).
c
      implicit none

      integer lv
      integer n

      integer i
      integer j
      integer k
      integer l
      integer m
      integer p
      real s
      real t
      complex v(lv,*)
      integer w(n-1)
      complex y
      complex z

      if ( n .eq. 1 ) go to 110
      l = 0
      m = 1
10    if ( l .eq. n ) go to 90
      k = l
      l = m
      m = m + 1
      p = l
      s = cabs ( v(l,l) )
      if ( m .gt. n ) go to 30

      do i = m, n
        t = cabs ( v(i,l) )
        if ( t .gt. s ) then
          p = i
          s = t
        end if
      end do

      w(l) = p
30    if ( s .eq. 0. ) go to 120
      y = v(p,l)
      v(p,l) = v(l,l)
      v(l,l) = -1.0E+00
      y = 1.0E+00 / y
      do i = 1, n
        v(i,l) = -y*v(i,l)
      end do
      j = l
50    j = j + 1
      if ( j .gt. n ) j = 1
      if ( j .eq. l ) go to 10
      z = v(p,j)
      v(p,j) = v(l,j)
      v(l,j) = z
      if ( cabs ( z ) .eq. 0. ) go to 50
      if ( k .eq. 0 ) go to 70
      do i = 1, k
        v(i,j) = v(i,j) + z*v(i,l)
      end do
70    v(l,j) = y*z
      if ( m .gt. n ) go to 50
      do i = m, n
        v(i,j) = v(i,j) + z*v(i,l)
      end do
      go to 50
90    l = w(k)

      do i = 1, n
        y = v(i,l)
        v(i,l) = v(i,k)
        v(i,k) = y
      end do

      k = k - 1
      if ( k .gt. 0 ) go to 90
      return
110   if ( cabs ( v(1,1) ) .eq. 0. ) go to 120
      v(1,1) = 1.0E+00 / v(1,1)
      return
120   continue
      write ( *, '(a)' ) ' '
      write(*,*) 'error: matrix has no inverse'
      stop
      end
      subroutine czero ( z, p, nd, w )

c*********************************************************************72
c
cc CZERO computes zeros of a monic complex polynomial.
c
c  Discussion:
c
c    The polynomial has the form:
c
c      P(1) + P(2)X +...+ P(ND)X**(ND-1) + X**ND 
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, complex Z(ND), the estimated zeros of the polynomial.
c
c    Input, complex P(ND), the polynomial coefficients.
c
c    Input, integer ND, the number of coefficients in the polynomial.
c
c    Workspace, complex W(*).  The size necessary for W depends on the
c    problem.  A dimension between 4*ND+32 and 6*ND is usually enough.
c    On output, W(1) contains the number of computed zeros.
c
      implicit none

      integer nd

      real a
      real amag
      real b
      real c
      real d
      real e
      integer f
      integer g
      integer h
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      integer nc
      integer nz
      integer n1
      integer o
      integer o1
      integer o2
      integer o3
      integer o4
      complex p(nd)
      complex pi
      integer q
      complex r
      complex s
      real s1
      real t
      real t1
      real t2
      complex u
      complex v
      real v1
      complex w(*)
      complex x
      real y
      real y0
      complex z(nd)

      b=0.0
      n=nd

      do i = 1, n
        z(i) = p(i)
      end do

20    continue

      if ( n .le. 1 ) then
        return
      end if

      if ( n .eq. 2 ) then
        call qad ( z(2), z(1) )
        return
      end if
c
c  Remove zeros at the origin.
c
      if ( amag ( z(1) ) .eq. 0.0 ) then
        do i = 2, n
          z(i-1) = z(i)
        end do
        z(n) = cmplx(0.0,0.0)
        n = n - 1
        go to 20
      end if

      nz = n
c
c  Compute the machine precision.
c
      e = 1.0E+00
60    continue
      e = 0.5*e
      d = 1.0E+00 + e
      if ( d .gt. 1.0E+00 ) go to 60
      e = e + e

      m = 5
      pi = cmplx(0.0, 2.0 * acos ( -1.0E+00 ) )
      x = cmplx(1.0E+00,0.0)
      y = 65536.
      y0 = 1.0E+00 / y
70    continue
      n1 = n + 1
      o1 = n1 + n
      o2 = o1 + n
      o3 = o2 + n1 / 2
      o4 = o3 + n1 / 2
      o = o4 - 1
      call cc ( z, w(o1), w(o2), n, y, y0 )
      g = o4 + n1
      nc = 32*n
      k = 10
      l = 8
      q = 0
      j = 0
      f = 0
80    continue
      f = f + 1
      q = q + j
      if ( f .lt. 6 ) go to 90
      if ( l .gt. nc ) go to 370
90    continue
      l = l + l
c
c  Estimate the radius of the smallest circle containing a zero.
c
      call rad ( z, w, w(n1), k, n, a )
      call cef ( w(o4), z, w(o1), w(o2), w(o3), a, x, n, n1, t, t1, 
     &  y, y0 )
      t = n*t*e
      j = 1
      i = n1
      h = i
100   i = i / 2
      j = j + j
      if ( i .gt. 0 ) go to 100
      if ( h+h .eq. j ) go to 120
      h = j
      j = h + o
      do i = g, j
        w(i) = (0.0,0.0)
      end do
120   j = l + o4
      call ff ( w(o4), w(j), h, l )
      c = amag ( w(o4) )
      j = 0
      h = l - 1
c
c  Determine a starting guess.
c
      do i = 1, h
        d = amag ( w(o4+i) )
        if ( d .le. c ) then
          j = i
          c = d
        end if
      end do

      s = a*x * cexp ( j * pi / l )
      s1 = cabs ( s )
      if ( s1 .lt. 1.0E+00 ) s1 = 1.0E+00
      x = x * cexp ( pi / ( 4 * l ) )
      x = x / cabs ( x )
c
c  Jenkins-Traub iterations.
c
      do i = 1, m
        call jt ( z, w, n, v, s, s1 )
      end do

      c = 0.0
      j = 0
150   call jt ( z, w, n, v, s, s1 )
      v1 = amag ( v )
      if ( v1 .eq. 0. ) go to 240
      if ( s1 .gt. 1.0E+00 ) go to 170
      u = (1.0E+00,0.0)
      i = n
160   i = i - 1
      u = w(i) + u*s
      if ( i .gt. 1 ) go to 160
      if ( amag ( u ) .eq. 0.0 ) go to 80
      go to 190
170   u = w(1)
      r = 1.0E+00 / s
      do i = 2, n
        u = w(i) + u*r
      end do
      if ( amag ( u ) .eq. 0.0 ) go to 80
      u = u*r
190   r = s - v / u
      j = j + 1
      d = amag ( r - s )
      s = r
      if ( s1 .le. a ) go to 200
      call err ( z, n, s1, t2 )
      t2 = n*t2*e
      if ( v1 .le. t2 ) go to 240
      go to 210
200   if ( n * alog ( t1 / s1 ) + alog ( t / v1 ) .ge. 0.0E+00 ) then
        go to 230
      end if

210   s1 = cabs ( s )
      if ( s1 .lt. 1.0E+00 ) s1 = 1.0E+00
      if ( j .lt. 3 ) go to 220
      if ( .5*d .gt. b+c ) go to 80
      if ( j .lt. m ) go to 220
      if ( 4.0*d .gt. b+c ) go to 80
220   b = c
      c = d
      go to 150
230   call err ( z, n, s1, t2 )
      t2 = n*e*t2
      if ( v1 .gt. t2 ) go to 210
240   r = z(n) + s
      i = n
250   i = i - 1
      v = z(i) + r*s
      z(i) = r
      r = v
      if ( i .gt. 1 ) go to 250
      z(n) = s
      q = q + j
      n = n - 1
      if ( n .gt. 2 ) go to 70
      call qad ( z(2), z(1) )
      nc = nd
      n = 1
260   l = n
      n = nd
      m = n + n
      j = 1
      k = n - 1
      do i = 2, k
        w(j) = j*p(i)
        w(j+n) = j*i*p(i+1)
        w(j+m) = amag ( p(j) )
        j = i
      end do
      w(k) = k*p(n)
      w(n) = n
      w(m-1) = n*k
      w(m) = 0.0
      w(m+k) = amag ( p(k) )
      w(m+n) = amag ( p(n) )
      o1 = m + 1
      n1 = n + 1
c
c  Refine the zeros.
c
      do 340 k = l, nz
        s = z(k)
        do j = 1, 3
          call rfn ( p, w, w(n1), w(o1), n, s, v, c, t )
          t = n*t*e
          if ( c .le. t ) go to 330
        end do
        a = cabs ( s )
        if ( a .gt. 1.0E+00 ) go to 300
        v = (1.0E+00,0.0)
        t = 1.0E+00
        i = n
290     v = v*s + p(i)
        t = t*a + real ( w(i+m) )
        i = i - 1
        if ( i .gt. 0 ) go to 290
        go to 320
300     r = 1.0E+00 / s
        a = 1.0E+00 / a
        v = (0.0,0.0)
        t = 0.
        do i = 1, n
          v = v*r + p(i)
          t = t*a + real ( w(i+m) )
        end do
        v = v*r + cmplx(1.0E+00,0.0)
        t = t*a + 1.0E+00
320     a = amag ( v )
        if ( a .gt. n*t*e ) go to 340
330     z(k) = s
340   continue

      do i = l, n
        w(n1-i) = z(i)
      end do
      l = n1 - l
      do  i = 1, l
        z(i) = w(i)
      end do
      w(1) = nc
      return
370   nc = nd - n + 1
      write(*,*) 'since the stopping criterion has not been satisfied'
      write(*,*) 'after',q,'iterations, we stop while computing'
      write(*,*) 'root number',nc
      z(n) = s
      go to 260
      end
      subroutine dag ( e, v, lv, a, n, d, p, w )

c*********************************************************************72
c
cc DAG diagonalizes a complex Hessenberg matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input, integer LV, the leading dimension of the array V.
c
      implicit none

      integer lv

      complex a(*)
      real amag
      real b
      real d(*)
      complex e(*)
      integer f
      integer f0
      integer g
      integer g0
      integer h
      integer i
      integer ig
      integer ii
      integer j
      integer k
      integer kl
      integer km
      integer ks
      integer k0
      integer k1
      integer l
      integer ll
      integer l0
      integer l1
      integer m
      integer n
      integer ns
      integer o
      real p(*)
      real q
      real r
      real s
      real sqr
      real t
      real t0
      real t1
      real t2
      real u
      complex v(lv,*)
      complex w(*)
      complex x
      complex y
      complex z
      complex z1
      complex z2
      complex z3
      complex z4

      do i = 1, n
        d(i) = 1.0E+00
      end do

      if ( n .eq. 1 ) go to 440
      b = 65536.0E+00**(-3)
      t = 1.0E+00
20    t = 0.5 * t
      s = 1.0E+00 + t
      if ( s .gt. 1.0E+00 ) go to 20
      t0 = t + t
      t2 = t0*t0
      ns = 50*n
      ll = 0
      l1 = 0
      kl = (n*(n+3)) / 2 - 1
      m = n
      km = kl - m
30    i = km
      j = m
40    if ( amag ( a(i) ) .eq. 0. ) go to 50
      i = i - j
      j = j - 1
      if ( j .gt. 1 ) go to 40
50    k0 = j
      l = j + 1
      g = i + l
      f = g - 1
      f0 = f
      g0 = g
      l0 = l
      if ( j .eq. m ) go to 430
      do i = j, m
        p(i) = sqrt ( d(i) )
      end do
      s = 0.
      k = j
70    h = f - j
      t = p(k)
      do i = f, g
        s = s + t* amag ( a(i) ) * p(i-h)
      end do
      k = l
      l = l + 1
      f = f + k
      g = g + l
      if ( k .lt. m ) go to 70
      if ( k .gt. m ) go to 90
      g = g - 1
      go to 70
90    if ( s .eq. 0. ) s = 1.0E+00
      t1 = 1.0E+00 / ( s * t0 )
      go to 410
100   ll = ll + 1
      if ( ll .gt. ns ) go to 550
      t = 1.0E+00
      do i = k0, m
        if ( d(i) .lt. t ) t = d(i)
      end do
      if ( t .gt. b ) go to 210

      do i = k0, m
        p(i) = sqrt ( d(i) )
        d(i) = 1.0E+00
      end do

      f = f0
      g = g0
      l = l0
      k = k0
130   t = p(k)

      if ( k0 .eq. 1 ) go to 150

      ii = f - k0 + 1
      h = f - 1
      do i = ii, h
        a(i) = a(i)*t
      end do

150   h = f - k0
      do i = f, g
        a(i) = a(i)*t*p(i-h)
      end do
      do i = 1, n
        v(i,k) = v(i,k)*t
      end do
      k = l
      l = l + 1
      f = f + k
      g = g + l
      if ( k .lt. m ) go to 130
      if ( k .gt. m ) go to 180
      g = g - 1
      go to 130
180   g = g - 1
190   if ( k .gt. n ) go to 210
      h = f - k0
      do i = f, g
        a(i) = a(i)*p(i-h)
      end do
      k = k + 1
      f = f + k
      g = g + k
      go to 190
210   s = d(m-1)
      t = d(m)
      z1 = a(km-1)*s
      z2 = a(km)*s
      z3 = a(kl-1)*t
      z4 = a(kl)*t
      call ceig ( z, z, z1, z2, z3, z4 )
      k = k0
      l = l0
      f = f0
      g = g0
      ks = k
      k1 = k - 1
      z1 = d(k)*a(f) - z
      z2 = d(k)*a(g)
220   t = amag ( z1 ) + amag ( z2 )
      if ( t .eq. 0. ) go to 240
      t = 1.0E+00 / t
      q = sqr(z1,t)
      r = sqr(z2,t)
      q = d(k)*q
      r = d(l)*r
      if ( q .gt. r ) go to 230
      z4 = z1 / z2
      z3 = (d(k) / d(l))*conjg ( z4 )
      e(k) = z3
      w(k) = z4
      p(k) = 1.0E+00
      s = r / (q+r)
      r = d(k)
      d(k) = d(l)*s
      q = d(l)
      d(l) = r*s
      go to 250
230   z4 = z2 / z1
      z3 = (d(l) / d(k))*conjg ( z4 )
      e(k) = z3
      w(k) = z4
      p(k) = 0.
      s = q / (q+r)
      d(k) = d(k)*s
      q = d(l)
      d(l) = d(l)*s
      go to 250
240   p(k) = 0.
      z3 = (0.,0.)
      z4 = (0.,0.)
      e(k) = z3
      w(k) = z4
      q = d(l)
250   if ( k .gt. ks ) go to 270
      if ( p(k) .eq. 0. ) go to 260
      y = a(g) + z3*a(f)
      a(g) = z4*a(g) - a(f)
      a(f) = y
      go to 290
260   y = a(f) + z3*a(g)
      a(g) = a(g) - z4*a(f)
      a(f) = y
      go to 290
270   i = g - 1
      if ( p(k) .eq. 0. ) go to 280
      x = a(h)*z3 + x
      a(h) = x
      y = a(g) + z3*a(i)
      a(g) = z4*a(g) - a(i)
      a(i) = y
      if ( amag ( x ) .ne. 0. ) go to 290
      k0 = k
      l0 = l
      f0 = i
      g0 = g
      go to 290
280   x = x*z3 + a(h)
      a(h) = x
      y = a(i) + z3*a(g)
      a(g) = a(g) - z4*a(i)
      a(i) = y
      if ( amag ( x ) .ne. 0. ) go to 290
      k0 = k
      l0 = l
      f0 = i
      g0 = g
290   j = k
      k = l
      l = l + 1
      do 310 i = ks, j
        o = g + i
        h = o + 1
        if ( p(i) .eq. 0. ) go to 300
        y = a(h) + e(i)*a(o)
        a(h) = w(i)*a(h) - a(o)
        a(o) = y
        go to 310
300     y = a(o) + e(i)*a(h)
        a(h) = a(h) - w(i)*a(o)
        a(o) = y
310   continue
      z3 = conjg ( z3 )
      z4 = conjg ( z4 )
      o = f - k1
      if ( p(j) .eq. 0. ) go to 340
      z1 = q*a(g+k) - z*w(j)

      do i = o, g
        h = i + k
        y = a(h) + z3*a(i)
        a(h) = z4*a(h) - a(i)
        a(i) = y
      end do

      do i = 1, n
        y = v(i,k) + z3*v(i,j)
        v(i,k) = z4*v(i,k) - v(i,j)
        v(i,j) = y
      end do

      f = f + k
      h = g
      g = g + l
      if ( l .gt. m ) go to 370
      z2 = q*a(g)
      x = a(g)
      a(g) = z4*a(g)
      go to 220
340   z1 = q*a(g+k) - z

      do i = o, g
        h = i + k
        y = a(i) + z3*a(h)
        a(h) = a(h) - z4*a(i)
        a(i) = y
      end do

      do i = 1, n
        y = v(i,j) + z3*v(i,k)
        v(i,k) = v(i,k) - z4*v(i,j)
        v(i,j) = y
      end do

      f = f + k
      h = g
      g = g + l
      if ( l .gt. m ) go to 370
      z2 = q*a(g)
      x = z3*a(g)
      go to 220
370   ig = g
      ii = l
380   if ( ii .gt. n ) go to 410
      ii = ii + 1

      do i = ks, j
        o = ig + i
        h = o + 1
        if ( p(i) .ne. 0. ) then
          y = a(h) + e(i)*a(o)
          a(h) = w(i)*a(h) - a(o)
          a(o) = y
        else
          y = a(o) + e(i)*a(h)
          a(h) = a(h) - w(i)*a(o)
          a(o) = y
        end if
      end do

      ig = ig + ii
      go to 380
410   t = d(m)
      s = d(m-1)
      q = amag ( a(km) )
      if ( ((t*q)*t1)*((s*q)*t1) .gt. 1.0E+00 ) go to 100
      l1 = l1 + 1
      if ( l1 .gt. 30 ) go to 420
      r = amag ( a(kl) )
      u = max ( q, r )
      if ( u .eq. 0. ) go to 420
      if ( (s*q)*(q / u) .gt. t2*(t*r)*(r / u) ) go to 100
420   l1 = 0
      e(m) = a(kl)*d(m)
      kl = km - 1
      km = km - m
      m = m - 1
      if ( m .gt. k0 ) go to 410
430   e(m) = a(kl)*d(m)
      kl = km - 1
      km = km - m
      m = m - 1
      if ( m .gt. 1 ) go to 30
440   e(1) = a(1)*d(1)
      p(1) = n
      km = (n*(n+1)) / 2 - 1
      k = n
450   x = e(k)
      w(k) = (1.0E+00,0.)
      l = k
      k = k - 1
      if ( k .eq. 0 ) go to 520
      do i = 1, k
        w(i) = -a(i+km)
      end do
      km = km - l
      f = km
      j = k
470   y = a(f+j) - x / d(j)
      z = (0.,0.)
      if ( amag ( y ) .ne. 0. ) z = w(j) / y 
      w(j) = z
      if ( j .eq. 1 ) go to 490
      g = f
      f = f - j
      j = j - 1
      do i = 1, j
        w(i) = w(i) - z*a(i+g)
      end do
      go to 470
490   x = w(l)
      do i = 1, n
        v(i,l) = x*v(i,l)
      end do
      do j = 1, k
        x = w(j)
        do i = 1, n
          v(i,l) = v(i,l) + x*v(i,j)
        end do
      end do
520   t = 0.
      do i = 1, n
        s = amag ( v(i,l) )
        if ( t .lt. s ) t = s
      end do
      if ( t .ne. 0. ) t = 1.0E+00 /t
      do i = 1, n
        v(i,l) = v(i,l)*t
      end do
      if ( k .gt. 0 ) go to 450
      return
550   m = n - m + 1
      write(*,*) 'since the stopping criterion not satisfied'
      write(*,*) 'after',ns,'iterations, we stop while computing'
      write(*,*) 'eigenvalue number',m
      p(1) = m
      return
      end
      subroutine db ( x, g )

c*********************************************************************72
c
cc DB
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real g(*)
      real x

      x = g(1)

      return
      end
      function det ( iexp, a )

c*********************************************************************72
c
cc DET computes the determinant of a matrix factored by FACT.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, integer IEXP, a power of ten that is part of the
c    determinant.
c
c    Input, real A(*), factorization information from FACT.
c
c    Output, real DET, the mantissa of the determinant.
c    Determinant = DET * 10.0^IEXP.
c
      implicit none

      real a(*)
      real c
      real d
      real det
      real f
      real g
      integer h
      integer i
      integer iexp
      integer j
      integer k
      integer l
      integer m
      integer n

      iexp=0
      det=0.0
      d = a(1)

      if ( abs ( d ) .ne. 1230 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor before computing determinant'
        return
      end if

      iexp = 0
      if ( d .lt. 0.0E+00 ) go to 70
      n = a(2)
      if ( n .eq. 1 ) go to 80
      d = 1.0E+00
      f = 2.0E+00**64
      g = 1.0E+00 /f
      h = 64
      m = n + 1
      j = 0
      k = 4
      l = 3 - m + m*n

      do i = k, l, m
        j = j + 1
        if ( a(i) .gt. j ) d = -d
        d = d*a(i+j)
20      if ( abs ( d ) .lt. f ) go to 30
        iexp = iexp + h
        d = d*g
        go to 20
30      if ( abs ( d ) .gt. g ) go to 40
        iexp = iexp - h
        d = d*f
        go to 30
40      continue
      end do

      d = d*a(l+m)
      if ( iexp .ne. 0 ) go to 50
      det = d
      return
50    if ( d .eq. 0. ) go to 90
      c = alog10 ( abs ( d ) ) + iexp * alog10 ( 2.0E+00 )
      iexp = c
      c = c - iexp
      if ( c .le. 0.0 ) go to 60
      c = c - 1
      iexp = iexp + 1
60    f = 10.0E+00**c
      if ( d .lt. 0. ) f = -f
      det = f
      return
70    det = 0.
      return
80    det = a(5)
      return
90    iexp = 0
      go to 70
      end
      subroutine diag ( e, v, lv, a, la, n )

c*********************************************************************72
c
cc DIAG diagonalizes a general real matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input, integer LV, the leading dimension of the array V.
c
      implicit none

      real a(*)
      complex e(*)
      integer i
      integer j
      integer k
      integer l
      integer la
      integer lv
      integer m
      integer n
      integer o
      integer p
      real v(*)

      if ( n .le. 1 ) then
        e(1) = a(1)
        v(1) = 1.0E+00
        return
      end if

      if ( lv .lt. n ) then
        write(*,*) 'lv argument in sub diag must be .ge. n argument'
        stop
      end if

      call ahess ( a, la, n, e )
      call sim ( v, n, a )
      o = n + 1
      p = n + 2
      i = 1 + n*o
      j = n + n - 3
30    a(i+j) = a(i)
      i = i - 1
      if ( i .gt. 0 ) go to 30
      j = 1
      k = n*p - o
      m = n + n
      l = m + 1
40    do i = m, l
        a(j) = a(i)
        a(j+1) = 0.
        j = j + 2
      end do
      m = m + o
      l = l + p
      if ( j .lt. k ) go to 40
      if ( j .gt. k ) go to 60
      l = l - 1
      go to 40
60    i = n*n
      j = i + i
70    v(j) = 0.
      v(j-1) = v(i)
      i = i - 1
      j = j - 2
      if ( i .gt. 0 ) go to 70
      i = k + n + n
      j = i + n + n
      k = j + n
      call dag ( e, v, n, a, n, a(k), a(j), a(i) )
      k = 2*(lv-n)
      i = 2*n*n
      j = n*k - k
80    l = i - n - n
90    v(i+j) = v(i)
      i = i - 1
      if ( i .gt. l ) go to 90
      j = j - k
      if ( i .gt. 0 ) go to 80
      return
      end
      function econ ( a, b )

c*********************************************************************72
c
cc ECON computes the condition number of an upper Hessenberg matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c    William Hager,
c    Condition Estimates,
c    SIAM Journal on Scientific and Statistical Computing,
c    Volume 5, Number 2, June 1984, pages 311-316.
c
c  Parameters:
c
c    Input, real A(*), the factorization information from EFACT.
c
c    Workspace, real B(N).
c
c    Output, real ECON, the condition number of the matrix.
c
      implicit none

      real a(*)
      real b(*)
      real c
      real d
      real econ
      integer i
      integer j
      integer m
      integer n

      c = a(1)

      if ( abs ( c ) .ne. 1237 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor with EFACT before',
     &  'estimating condition'
        stop
      end if

      if ( c .lt. 0.0E+00 ) then
        econ = -1.0E+00
        return
      end if

      m = 0
      n = a(2)
      c = 1.0E+00 / a(2)
      do j = 1, n
        b(j) = c
      end do
      go to 50

30    continue

      do j = 1, n
        b(j) = 0.0E+00
      end do

      b(m) = 1.0E+00

50    continue

      call esolve ( b, a, b )

      c = 0.0E+00
      do j = 1, n
        c = c + abs ( b(j) )
        if ( b(j) .lt. 0.0E+00 ) then
          b(j) = -1.0E+00
        else
          b(j) = 1.0E+00
        end if
      end do

      call etrans ( b, a, b )
c
c  I is the index of the largest entry of B.
c
      i = 1
      do j = 1, n
        if ( abs ( b(i) ) .lt. abs ( b(j) ) ) then
          i = j
        end if
      end do

      if ( m .eq. 0 ) go to 90
      if ( m .eq. i ) go to 100
      if ( d .ge. c ) go to 100
90    m = i
      d = c
      go to 30
100   c = c*a(3)
      c = max ( c, 1.0E+00 )
      econ = c
      return
      end
      function edet ( iexp, a )

c*********************************************************************72
c
cc EDET computes the determinant of an upper Hessenberg matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, integer IEXP, a power of ten that is part of the
c    determinant.
c
c    Input, real A(*), factorization information from EFACT.
c
c    Output, real EDET, the mantissa of the determinant.
c    Determinant = EDET * 10.0^IEXP.
c
      implicit none

      real a(*)
      real c
      real d
      real edet 
      real f
      real g
      integer h
      integer i
      integer iexp
      integer j
      integer n

      iexp=0
      edet=0.0
      d = a(1)

      if ( abs ( d ) .ne. 1237 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor with EFACT before',
     &  ' computing determinant'
        return
      end if

      iexp = 0
      if ( d .lt. 0. ) go to 90
      d = 1.0E+00
      f = 2.0E+00**64
      g = 1.0E+00 /f
      h = 64
      n = a(2)
      i = 3
      j = 1
20    i = i + j
      j = j + 1
      d = d*a(i)
      if ( j .gt. n ) go to 50
30    if ( abs ( d ) .lt. f ) go to 40
      iexp = iexp + h
      d = d*g
      go to 30
40    if ( abs ( d ) .gt. g ) go to 20
      iexp = iexp - h
      d = d*f
      go to 40
50    i = 5 + (n*(n+1))/2
      j = i + n + n - 2
60    if ( a(i) .eq. 1.0E+00 ) d = -d
      i = i + 2
      if ( i .lt. j ) go to 60
      if ( iexp .ne. 0 ) go to 70
      edet = d
      return
70    if ( d .eq. 0. ) go to 100
      c = alog10 ( abs ( d ) ) + iexp * alog10 ( 2.0E+00 )
      iexp = c
      c = c - iexp
      if ( c .le. 0.0 ) go to 80
      c = c - 1
      iexp = iexp + 1
80    f = 10.0E+00**c
      if ( d .lt. 0. ) f = -f
      edet = f
      return
90    edet = 0.
      return
100   iexp = 0
      go to 90
      end
      subroutine ediag ( e, v, lv, a, n, w )

c*********************************************************************72
c
cc EDIAG diagonalizes an upper Hessenberg matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input, integer LV, the leading dimension of the array V.
c
C     |    INPUT:                                              |
C     |                                                        |                            |
C     |         A     --COEFFICIENTS OF HESSENBERG MATRIX      |
C     |                 PACKED AT START OF REAL ARRAY          |
C     |                 (LENGTH AT LEAST (N+1)(N+2) - 4)       |
C     |                                                        |
C     |         N     --DIMENSION OF MATRIX STORED IN A        |
C     |                                                        |
C     |         W     --REAL WORK ARRAY WITH AT LEAST          |
C     |                 4N ELEMENTS                            |
C     |                                                        |
C     |    OUTPUT:                                             |
C     |                                                        |
C     |         E     --COMPLEX ARRAY OF EIGENVALUES           |
C     |                                                        |
C     |         V     --COMPLEX ARRAY OF EIGENVECTORS  
c
      implicit none

      integer lv

      real a(*)
      real amag
      complex e(*)
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      real r
      real s
      real t
      complex v(lv,*)
      real w(*)
      complex z

      j = (n+1)*(n+2) - 4
      i = j/2
10    a(j-1) = a(i)
      a(j) = 0.0E+00
      i = i - 1
      j = j - 2
      if ( i .gt. 0 ) go to 10

      do j = 1, n
        do i = 1, n
          v(i,j) = 0.0E+00
        end do
        v(j,j) = 1.0E+00
      end do

      m = n + 1
      l = m + n
      call dag ( e, v, lv, a, n, w, w(m), w(l) )
      j = m
40    j = j - 1
      if ( j .le. 1 ) go to 90
      z = conjg ( e(j) )
      r = abs ( aimag(z) )
      s = r
      l = j - 1
      do i = 1, l
         t = amag ( e(i) - z )
         if ( t .lt. s ) then
           k = i
           s = t
         end if
      end do
      if ( r .gt. s ) go to 70
      e(j) = real ( z )
      do i = 1, n
        v(i,j) = real ( v(i,j) )
      end do
      go to 40
70    e(k) = e(l)
      e(l) = z

      do i = 1, n
        v(i,k) = v(i,l)
        v(i,l) = conjg ( v(i,j) )
      end do

      j = l
      go to 40
90    if ( j .lt. 1 ) return
      e(j) = real ( e(j) )

      do i = 1, n
        v(i,j) = real ( v(i,j) )
      end do

      return
      end
      subroutine efact ( a, n )

c*********************************************************************72
c
cc EFACT factors an upper Hessenberg matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input/output, real A(1 + (N*(N+5))/2).  On input, the initial entries
c    of the array contain the upper Hessenberg portion of the matrix, 
c    stored by rows.  On output, factorization information. 
c
c    Input, integer N, the number of rows and columns in the matrix.
c
      implicit none

      real a(*)
      integer g
      integer h
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      real r
      real s
      real t

      g = (n*(n+1))/2
      h = g + n + 1
      m = 1
      l = 2
      j = 0
      k = 1
      r = 0.
10    if ( m .lt. g ) go to 20
      if ( m .gt. g ) go to 40
      l = l - 1
20    s = 0.
      do i = m, l
        t = a(i)
        s = s + abs ( t )
        a(i-j) = t
      end do
      if ( s .gt. r ) r = s
      j = k
      k = k + 1
      a(h+j) = t
      m = l + 1
      l = m + k
      go to 10
40    g = g + 3
      i = g + 1
      j = h + 1
      k = h + n
50    if ( i .gt. k ) go to 60
      a(i) = a(j)
      i = i + 2
      j = j + 1
      go to 50
60    i = g
70    a(i) = a(i-3)
      i = i - 1
      if ( i .gt. 3 ) go to 70
      a(1) = 1237
      a(2) = n
      a(3) = r
      k = n
      m = g
      l = g - n + 1
      g = h + n - 1
80    if ( k .eq. 1 ) go to 140
      j = k
      k = k - 1
      s = a(m)
      t = a(g)
      if ( abs ( s ) .ge. abs ( t ) ) go to 100
      a(m) = t
      t = s/t
      a(g) = t
      a(g+1) = 1.0E+00
      g = g - 2
      l = l - k
      m = m - j
      do i = l, m
        h = i + k
        s = a(i)
        a(i) = a(h) - t*s
        a(h) = s
      end do
      go to 80
100   if ( s .eq. 0. ) go to 130
      t = t/s
      l = l - k
      m = m - j
      if ( t .eq. 0. ) go to 120
      do i = l, m
        a(i) = a(i) - t*a(i+k)
      end do
120   a(g) = t
      a(g+1) = 0.
      g = g - 2
      go to 80
130   a(1) = -1237
      go to 120
140   if ( a(4) .eq. 0. ) a(1) = -1237
      return
      end
      subroutine eig2 ( e1, e2, a, b, c, d )

c*********************************************************************72
c
cc EIG2
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real a
      real b
      real c
      real d
      complex e1
      complex e2
      real q
      real r
      real s
      real t
      real u
      real v
      real w

      q = abs ( b )
      r = abs ( c )
      t = sqrt ( q ) * sqrt ( r )
      u = .5*(a-d)
      s = abs ( u )
      if ( s .lt. t ) go to 20

      if ( t .eq. 0.0 ) then
        e1 = a
        e2 = d
        return
      end  if

      w = 1.0E+00 /t
      s = 1.0E+00 /s
      v = (s*u)**2+((b*w)*(c*w))*(s*t)**2
      go to 30
20    s = 1.0E+00 /t
      v = (s*u)**2+(b*s)*(c*s)
30    if ( v .gt. 0. ) go to 40
      e1 = cmplx(.5*(a+d), sqrt ( - v ) * t )
      e2 = conjg ( e1 )
      return
40    v = s * abs ( u ) + sqrt ( v )
      if ( u .gt. 0. ) go to 60
      if ( q .gt. r ) go to 50
      e1 = a - b*(c*s)/v
      e2 = d + b*(c*s)/v
      return
50    e1 = a - (b*s)*c/v
      e2 = d + (b*s)*c/v
      return
60    if ( q .gt. r ) go to 70
      e1 = a + b*(c*s)/v
      e2 = d - b*(c*s)/v
      return
70    e1 = a + (b*s)*c/v
      e2 = d - (b*s)*c/v
      return
      end
      subroutine eig3 ( ea, eb, a, b, y, z )

c*********************************************************************72
c
cc EIG3
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real a
      real b
      real c
      real ea
      real eb
      real s
      real t
      real y
      real z

      t = .5*(b-a)
      c = sqrt ( abs ( y ) ) * sqrt ( abs ( z ) )
  
      if ( abs ( t ) .gt. abs ( c ) ) then
        t = abs ( c ) / t
        s = t * abs ( c ) / ( 1.0E+00 + sqrt ( 1.0E+00 + t * t ) )
        ea = a - s
        eb = b + s
        return
      end if

      if ( c .eq. 0. ) then
        ea = a
        eb = b
        return
      end if

      t = t / abs ( c )
      s = abs ( c ) / ( abs ( t ) + sqrt ( 1.0E+00 + t * t ) )

      if ( t .ge. 0. ) then
        ea = a - s
        eb = b + s
      else
        ea = a + s
        eb = b - s
      end if

      return
      end
      subroutine emult ( y, x, a, n )

c*********************************************************************72
c
cc EMULT multiplies an upper Hessenberg matrix times a vector.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, real Y(N), the product of the matrix times X.
c
c    Input, real X(N), the vector to be multiplied.
c
c    Input, real A((N*(N+3))/2-1), the upper Hessenberg matrix, packed into the array.
c
c    Input, integer N, the number of rows and columns in the matrix.
c
      implicit none

      integer n

      real a(*)
      integer i
      integer j
      integer k
      real t
      real x(n)
      real y(n)

      if ( n .le. 1 ) then
        y(1) = a(1) * x(1)
        return
      end if

      y(1) = 0.0
      j = 0
      k = 1
20    t = x(k)
      do i = 1, k
        y(i) = y(i) + t*a(i+j)
      end do
      if ( k .eq. n ) return
      k = k + 1
      j = j + k
      y(k) = t*a(j)
      go to 20
      end
      subroutine eql ( x, f, g, h, n )

c*********************************************************************72
c
cc EQL
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real f(*)
      real g(*)
      real h(*)
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      integer nm1
      real r
      real s
      real u
      real x(*)
      real z

      z = 0.
      i = 1
      r = x(n)
      s = f(n)
      u = g(n)
      if ( s .eq. z ) go to 30
10    if ( f(i) .eq. z ) go to 20
      if ( u .lt. g(i) ) go to 30
      if ( u .gt. g(i) ) go to 20
      if ( abs ( f(i) ) .ge. abs ( s ) ) go to 30
20    i = i + 1
      if ( i .lt. n ) go to 10
      go to 50
30    m = n + i
      nm1 = n - 1
      do j = i, nm1
        k = m - j
        l = k - 1
        x(k) = x(l)
        f(k) = f(l)
        g(k) = g(l)
      end do
      x(i) = r
      f(i) = s
      g(i) = u
50    u = g(n)

      do 80 i = 1, n
        if ( f(i) .eq. z ) go to 60
        if ( g(i)-u .gt. -99.0 ) go to 70
60      h(i) = z
        go to 80
70      h(i) = f(i)*16.0**(g(i)-u)
80    continue

      return
      end
      subroutine err ( z, n, s, t )

c*********************************************************************72
c
cc ERR
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      integer n

      real amag
      integer i
      real r
      real s
      real t
      complex z(n)

      if ( s .le. 1.0E+00 ) then

        t = 1.0E+00
        do i = n, 1, -1
          t = t * s + amag ( z(i) )
        end do
 
      else

        r = 1.0E+00 / s
 
        t = amag ( z(1) )
        do i = 2, n
          t = t * r + amag ( z(i) )
        end do

        t = t * r + 1.0E+00

      end if

      return
      end
      subroutine esolve ( x, a, b )

c*********************************************************************72
c
cc ESOLVE solves an upper Hessenberg linear system factored by EFACT.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, real X(N), the solution.
c
c    Input, real A(*), factorization information from EFACT.
c
c    Input, real B(N), the right hand side.
c
      implicit none

      real a(*)
      real b(*)
      integer i
      integer j
      integer k
      integer l
      integer n
      real t
      real x(*)

      t = a(1)

      if ( abs ( t ) .ne. 1237 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor with efact before solving'
        stop
      end if

      n = a(2)
      j = 3 + (n*(n-1))/2
      k = n
      l = j + n
      if ( t .lt. 0. ) go to 80
      do i = 1, n
        x(i) = b(i)
      end do
30    t = b(k)/a(j+k)
40    x(k) = t
      if ( k .eq. 1 ) go to 70
      k = k - 1
      if ( t .eq. 0.0E+00 ) go to 60
      do i = 1, k
        x(i) = x(i) - t*a(i+j)
      end do
60    j = j - k
      go to 30
70    if ( k .eq. n ) return
      j = k
      k = k + 1
      l = l + 2
      x(k) = x(k) - a(l-1)*x(j)
      if ( a(l) .eq. 0. ) go to 70
      t = x(k)
      x(k) = x(j)
      x(j) = t
      go to 70
80    j = 3
      k = 0
90    j = j + k
      k = k + 1
      if ( a(j+k) .ne. 0. ) go to 90

      do i = 1, n
        x(i) = 0.0
      end do

      t = 1.0E+00
      go to 40
      end
      subroutine etrans ( x, a, b )

c*********************************************************************72
c
cc ETRANS solves a transposed upper Hessenberg linear system factored by EFACT.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, real X(N), the solution.
c
c    Input, real A(*), factorization information from EFACT.
c
c    Input, real B(N), the right hand side.
c
      implicit none

      real a(*)
      real b(*)
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      real t
      real x(*)

      t = a(1)

      if ( abs ( t ) .ne. 1237 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor with efact before solving'
        stop
      end if

      n = a(2)
      k = n
      i = 3 + (n*(n+5))/2
      if ( t .lt. 0.0E+00 ) go to 80
      do j = 1, n
        x(j) = b(j)
      end do
30    if ( k .lt. 2 ) go to 50
      j = k
      k = k - 1
      i = i - 2
      if ( a(i) .eq. 0.0E+00 ) go to 40
      t = x(j)
      x(j) = x(k)
      x(k) = t
40    x(k) = x(k) - a(i-1)*x(j)
      go to 30
50    j = 4
      t = x(k)
60    x(k) = t/a(j)
      l = k
      k = k + 1
      if ( k .gt. n ) return
      t = x(k)
      do i = 1, l
        t = t - a(i+j)*x(i)
      end do
      j = j + k
      go to 60
80    j = i - n - n
90    if ( a(j) .eq. 0.0E+00 ) go to 100
      j = j - k
      k = k - 1
      go to 90
100   do i = 1, n
        x(i) = 0.0E+00
      end do
      x(k) = 1.0E+00
      l = k
      go to 130
120   x(k) = t/a(j)
130   if ( k .eq. n ) return
      m = k
      k = k + 1
      t = x(k)
      do i = l, m
        t = t - a(i+j)*x(i)
      end do
      j = j + k
      go to 120
      end
      subroutine evals ( e, a, n, w )

c*********************************************************************72
c
cc EVALS computes the eigenvalues of an upper Hessenberg matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
C     |         A     --COEFFICIENTS OF HESSENBERG MATRIX      |
C     |                 PACKED AT START OF REAL ARRAY          |
C     |                 (LENGTH AT LEAST (N+1)(N+2) - 4)       |
C     |                                                        |
C     |         N     --DIMENSION OF MATRIX STORED IN A        |
C     |                                                        |
C     |         W     --REAL WORK ARRAY WITH AT LEAST          |
C     |                 4N ELEMENTS                            |
C     |                                                        |
C     |    OUTPUT:                                             |
C     |                                                        |
C     |         E     --COMPLEX ARRAY OF EIGENVALUES  
c
      implicit none

      real a(*)
      real amag
      complex e(*)
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      real r
      real s
      real t
      real w(*)
      complex z

      j = (n+1)*(n+2) - 4
      i = j/2
10    a(j-1) = a(i)
      a(j) = 0.0E+00
      i = i - 1
      j = j - 2
      if ( i .gt. 0 ) go to 10
      m = n + 1
      l = m + n
      call vls ( e, a, n, w, w(m), w(l) )
      i = m
20    i = i - 1
      if ( i .le. 1 ) go to 50
      z = conjg ( e(i) )
      r = abs ( aimag(z) )
      s = r
      l = i - 1

      do j = 1, l
        t = amag ( e(j) - z )
        if ( t .lt. s ) then
          k = j
          s = t
        end if
      end do

      if ( r .gt. s ) go to 40
      e(i) = real ( z )
      go to 20
40    e(k) = e(l)
      e(l) = z
      i = l
      go to 20
50    if ( i .eq. 1 ) e(i) = real ( e(i) )
      return
      end
      subroutine evect ( e, x, a, n )

c*********************************************************************72
c
cc EVECT computes an eigenvector of an upper Hessenberg matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input/output, real E, the estimated eigenvalue.  On output, the
c    estimate may have been improved.
c
c    Output, real X(N), the eigenvector.
c
C     A     --HESSENBERG MATRIX PACKED AT START OF
C                      ARRAY WITH LENGTH AT LEAST N(N+7)/2 - 2
C                      (ALGORITHM DESTROYS ELEMENTS)
c
c    Input, integer N, the dimension of the matrix.
c
      implicit none

      real a(*)
      real e
      integer f
      integer g
      integer h
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      integer o
      integer p
      integer q
      real r
      real s
      real t
      real u
      real v
      real x(*)

      if ( n .le. 1 ) then
        e = a(1)
        x(1) = 1.0E+00
        return
      end if

      g = (n*(n+1))/2
      f = g
      h = g + n - 2
      m = 1
      l = 2
      j = 0
      k = 1
20    if ( k .gt. n ) go to 40
      p = l - 1
      a(p) = a(p) - e
      if ( k .eq. n ) l = p
      do i = m, l
        a(i-j) = a(i)
      end do
      j = k
      k = k + 1
      a(h+k) = a(l)
      m = l + 1
      l = m + k
      go to 20
40    i = g + 1
      j = h + 2
      k = h + n
50    if ( i .gt. k ) go to 60
      a(i) = a(j)
      i = i + 2
      j = j + 1
      go to 50
60    k = n
      m = g
      l = g - n + 1
      g = h + n - 1
      r = abs ( a(m) ) + abs ( a(g) )
70    if ( k .eq. 1 ) go to 140
      j = k
      k = k - 1
      s = a(m)
      t = a(g)
      if ( abs ( s ) .ge. abs ( t ) ) go to 100
      a(m) = t
      if ( r .lt. abs ( t ) ) go to 80
      r = abs ( t )
      p = j
      o = m
80    t = s/t
      a(g) = t
      a(g+1) = 1.0E+00
      g = g - 2
      l = l - k
      m = m - j
      do i = l, m
        q = i + k
        s = a(i)
        a(i) = a(q) - t*s
        a(q) = s
      end do
      go to 70
100   if ( r .lt. abs ( t ) ) go to 110
      r = abs ( t )
      p = j
      o = m
110   if ( s .eq. 0.0E+00 ) go to 130
      t = t/s
      l = l - k
      m = m - j
      do i = l, m
        a(i) = a(i) - t*a(i+k)
      end do
130   a(g) = t
      a(g+1) = 0.0E+00
      g = g - 2
      go to 70
140   t = abs ( a(1) )
      if ( t .gt. r ) go to 150
      r = t
      p = k
150   do i = 1, n
        x(i) = 0.0E+00
      end do
      t = 1.0E+00
      l = f
      j = o - p
      k = p
      go to 180
170   t = x(k)/a(j+k)
180   x(k) = t
      if ( k .eq. 1 ) go to 200
      k = k - 1
      do i = 1, k
        x(i) = x(i) - t*a(i+j)
      end do
      j = j - k
      go to 170
200   if ( k .eq. n ) go to 210
      j = k
      k = k + 1
      l = l + 2
      x(k) = x(k) - a(l-1)*x(j)
      if ( a(l) .eq. 0.0E+00 ) go to 200
      t = x(k)
      x(k) = x(j)
      x(j) = t
      go to 200
210   if ( r .eq. 0.0E+00 ) go to 310
      s = 0.0E+00
      do i = 1, n
        t = abs ( x(i) )
        if ( t .gt. s ) s = t
      end do
      p = h + n
      r = 0.0E+00
      s = 1.0E+00 /s
      do i = 1, n
        t = s*x(i)
        r = r + t*t
        x(i) = t
        a(i+p) = t
      end do
      v = 0.0E+00
      l = f
      j = l - n
      k = n
240   t = x(k)/a(j+k)
      x(k) = t
      if ( k .eq. 1 ) go to 260
      k = k - 1
      do i = 1, k
        x(i) = x(i) - t*a(i+j)
      end do
      j = j - k
      go to 240
260   u = abs ( x(1) )
270   if ( k .eq. n ) go to 290
      j = k
      k = k + 1
      l = l + 2
      t = x(j)
      s = x(k) - a(l-1)*t
      x(k) = s
      if ( abs ( s ) .gt. u ) u = abs ( s )
      if ( a(l) .eq. 1.0E+00 ) go to 280
      v = v + a(p+j)*t
      go to 270
280   x(j) = s
      x(k) = t
      v = v + a(p+j)*x(j)
      go to 270
290   v = v + a(p+n)*x(n)
      if ( v .ne. 0.0E+00 ) v = r/v
      s = 0.0E+00
      t = 1.0E+00 /u
      do i = 1, n
        s = s + (v*x(i))**2
        x(i) = t*x(i)
      end do
      if ( r+r .ge. s ) e = e + v
      return
310   u = 0.0E+00
      do i = 1, n
        t = abs ( x(i) )
        if ( t .gt. u ) u = t
      end do
      t = 1.0E+00 /u
      do i = 1, n
        x(i) = t*x(i)
      end do

      return
      end
      subroutine evert ( v, lv, a )

c*********************************************************************72
c
cc EVERT inverts an upper Hessenberg matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, real V(LV,N), the NxN inverse of the matrix.
c
c    Input, integer LV, the leading dimension of V.
c
c    Input, real A(*), factorization information from EFACT.
c
      implicit none

      integer lv

      real a(*)
      integer i
      integer j
      integer n
      real t
      real v(lv,*)

      t = a(1)

      if ( abs ( t ) .ne. 1237 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor with efact before inverting'
        stop
      end if

      if ( t .le. 0.0E+00 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: matrix has no inverse'
        stop
      end if

      n = a(2)

      do j = 1, n
        do i = 1, n
          v(i,j) = 0.0E+00
        end do
        v(j,j) = 1.0E+00
        call esolve ( v(1,j), a, v(1,j) )
      end do

      return
      end
      subroutine fact ( a, la, n )

c*********************************************************************72
c
cc FACT factors a general matrix, using partial pivoting.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input/output, real A(3+N*(N+1)).  On input, information defining the
c    matrix, stored as an LAxN matrix.  On output, factorization information.
c    A(1) is 1230 if the factorization was successful, and -1230 if it failed.
c    A(2) contains the value of N.
c    A(3) contains the L1 norm of the matrix.
c
c    Input, integer LA, the leading dimension of A.
c
c    Input, integer N, the number of rows and columns in the matrix.
c
      implicit none

      real a(*)
      integer e
      integer f
      integer g
      integer h
      integer i
      integer j
      integer k
      integer l
      integer la
      integer m
      integer n
      integer o
      integer p
      real r
      real s
      real t

      if ( la .gt. n ) then
        call pack ( a, la, n )
      end if

      r = 0.0E+00
      o = n + 1
      p = o + 1
      l = 5 + n * p
      i = -n - 3

10    continue

      l = l - o

      if ( l .ne. 4 ) then

        s = 0.0E+00

        do k = 1, n
          j = l - k
          t = a(i+j)
          a(j) = t
          s = s + abs ( t )
        end do

        r = max ( r, s )
        i = i + 1
        go to 10

      end if

      a(1) = 1230
      a(2) = n
      a(3) = r

      i = 5 - p
      k = 1

20    continue

      i = i + p
      if ( k .eq. n ) then
        if ( a(i) .eq. 0.0E+00 ) then
          a(1) = -1230
        end if
        return
      end if

      e = n - k
      m = i + 1
      h = i
      l = i + e
      do j = m, l
        if ( abs ( a(j) ) .gt. abs ( a(h) ) ) then
          h = j
        end if
      end do
      g = h - i
      j = i - k
      a(j) = g + k
      t = a(h)
      a(h) = a(i)
      a(i) = t
      k = k + 1

      if ( t .eq. 0.0E+00 ) then
        a(1) = -1230
        go to 40
      end if

      do j = m, l
        a(j) = a(j) / t
      end do
      f = i + e * o

30    continue

      j = k + l
      h = j + g
      t = a(h)
      a(h) = a(j)
      a(j) = t
      l = e + j

      if ( t .ne. 0.0E+00 ) then
        h = i - j
        m = j + 1
        do j = m, l
          a(j) = a(j) - t * a(j+h)
        end do

      end if

      if ( l .lt. f ) then
        go to 30
      end if

      go to 20
      end
      function fd ( a, x, h, n, grad )

c*********************************************************************72
c
cc FD
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      integer n

      real a
      real d
      real fd
      external grad
      real h(n,3)
      integer i
      real x(n)

      do i = 1, n
        h(i,2) = x(i) + a * h(i,1)
      end do

      call grad ( h(1,3), h(1,2) )

      d = 0.0E+00
      do i = 1, n
        d = d + h(i,1) * h(i,3)
      end do

      fd = d

      return
      end
      subroutine ff ( a, b, n, l )

c*********************************************************************72
c
cc FF evaluates a Fast Fourier Transform for CZERO.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      complex a(*)
      complex b(*)
      integer i
      integer il
      integer i1
      integer i2
      integer j
      integer k
      integer l
      integer l1
      integer l2
      integer m
      integer m0
      integer m1
      integer n
      real p
      complex s
      complex t
      complex w

      m = n

      if ( l .lt. m ) then

        m = m - 1

        do j = l, m, l
          do i = 1, l
            a(i) = a(i) + a(i+j)
          end do
        end do

        m = l

      end if

      p = acos ( -1.0E+00 )
      w = cexp ( m * cmplx ( 0.0E+00, p ) / l )
      l1 = l - 1
      l2 = l/2

      do j = m, l1, m
        do i = 1, m
          a(i+j) = a(i)
        end do
      end do

10    continue

      m0 = m
      m = m/2
      m1 = m - 1
      if ( m .eq. 0 ) return
      s = 1
      t = -1
      i1 = 0
      do j = 1, l, m0
        il = j + m1
        i2 = i1 + l2
        do i = j, il
          k = i + m
          b(i+i1) = a(i) + s * a(k)
          b(i+i2) = a(i) + t * a(k)
        end do
        s = w * s
        t = w * t
        i1 = i1 - m
      end do

      do i = 1, l
        a(i) = b(i)
      end do

      w = csqrt ( w )
      go to 10
      end
      subroutine ffc ( a, n, w )

c*********************************************************************72
c
cc FFC computes the conjugate Fast Fourier Transform.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      complex a(*)
      integer d
      integer e
      integer f
      integer g
      integer h
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      integer np
      integer o
      integer p(25)
      complex s
      complex t
      complex u
      complex v
      complex w(*)

      save np
      save p

      data p /2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,
     & 73,79,83,89,97/

      data np /25/

      m = n
      f = 0
      u = cmplx(0.0E+00, -2. * acos ( -1.0E+00 ) / n )
10    if ( m .eq. 1 ) go to 900
      do i = 1, np
        if ( (m/p(i))*p(i) .eq. m ) go to 30
      end do
      l = m
      go to 40
30    l = p(i)
40    o = m
      m = m/l
      v = cexp ( m * u )
      s = (1.0E+00,0.0E+00)
      h = 0
      if ( f .eq. 1 ) go to 470
      if ( l .eq. 2 ) go to 50
      if ( l .eq. 3 ) go to 170
      go to 290
50    go to (150,130,110,90),m
60    j = -h
70    i = h + 1
      h = h + m
      e = j + m
      do k = i, h
        w(k) = a(j+k) + s*a(e+k)
      end do
      j = e
      s = s*v
      if ( j+h .lt. n ) go to 70
      if ( h .lt. n ) go to 60
      f = 1
      go to 10
90    j = -h
100   h = h + 1
      e = j + m
      w(h) = a(j+h) + s*a(e+h)
      h = h + 1
      w(h) = a(j+h) + s*a(e+h)
      h = h + 1
      w(h) = a(j+h) + s*a(e+h)
      h = h + 1
      w(h) = a(j+h) + s*a(e+h)
      j = e
      s = s*v
      if ( j+h .lt. n ) go to 100
      if ( h .lt. n ) go to 90
      f = 1
      go to 10
110   j = -h
120   h = h + 1
      e = j + m
      w(h) = a(j+h) + s*a(e+h)
      h = h + 1
      w(h) = a(j+h) + s*a(e+h)
      h = h + 1
      w(h) = a(j+h) + s*a(e+h)
      j = e
      s = s*v
      if ( j+h .lt. n ) go to 120
      if ( h .lt. n ) go to 110
      f = 1
      go to 10
130   j = -h
140   h = h + 1
      e = j + m
      w(h) = a(j+h) + s*a(e+h)
      h = h + 1
      w(h) = a(j+h) + s*a(e+h)
      j = e
      s = s*v
      if ( j+h .lt. n ) go to 140
      if ( h .lt. n ) go to 130
      f = 1
      go to 10
150   j = -h
160   h = h + 1
      e = j + m
      w(h) = a(j+h) + s*a(e+h)
      j = e
      s = s*v
      if ( j+h .lt. n ) go to 160
      if ( h .lt. n ) go to 150
      f = 1
      go to 10
170   go to (270,250,230,210),m
180   j = -h
190   i = h + 1
      h = h + m
      e = j + m
      d = e + m
      t = s*s
      do k = i, h
        w(k) = a(j+k) + s*a(e+k) + t*a(d+k)
      end do
      j = d
      s = s*v
      if ( j+h .lt. n ) go to 190
      if ( h .lt. n ) go to 180
      f = 1
      go to 10
210   j = -h
220   h = h + 1
      e = j + m
      d = e + m
      t = s*s
      w(h) = a(j+h) + s*a(e+h) + t*a(d+h)
      h = h + 1
      w(h) = a(j+h) + s*a(e+h) + t*a(d+h)
      h = h + 1
      w(h) = a(j+h) + s*a(e+h) + t*a(d+h)
      h = h + 1
      w(h) = a(j+h) + s*a(e+h) + t*a(d+h)
      j = d
      s = s*v
      if ( j+h .lt. n ) go to 220
      if ( h .lt. n ) go to 210
      f = 1
      go to 10
230   j = -h
240   h = h + 1
      e = j + m
      d = e + m
      t = s*s
      w(h) = a(j+h) + s*a(e+h) + t*a(d+h)
      h = h + 1
      w(h) = a(j+h) + s*a(e+h) + t*a(d+h)
      h = h + 1
      w(h) = a(j+h) + s*a(e+h) + t*a(d+h)
      j = d
      s = s*v
      if ( j+h .lt. n ) go to 240
      if ( h .lt. n ) go to 230
      f = 1
      go to 10
250   j = -h
260   h = h + 1
      e = j + m
      d = e + m
      t = s*s
      w(h) = a(j+h) + s*a(e+h) + t*a(d+h)
      h = h + 1
      w(h) = a(j+h) + s*a(e+h) + t*a(d+h)
      j = d
      s = s*v
      if ( j+h .lt. n ) go to 260
      if ( h .lt. n ) go to 250
      f = 1
      go to 10
270   j = -h
280   h = h + 1
      e = j + m
      d = e + m
      t = s*s
      w(h) = a(j+h) + s*a(e+h) + t*a(d+h)
      j = d
      s = s*v
      if ( j+h .lt. n ) go to 280
      if ( h .lt. n ) go to 270
      f = 1
      go to 10
290   go to (440,410,380,350),m
300   j = -h
310   t = (1.0E+00,0.0E+00)
      i = h + 1
      h = h + m
      g = j + o
      do k = i, h
        w(k) =  (0.0E+00,0.0E+00)
      end do
330   do k = i, h
        w(k) = w(k) + t*a(j+k)
      end do
      t = t*s
      j = j + m
      if ( j .lt. g ) go to 330
      j = j - m
      s = s*v
      if ( j+h .lt. n ) go to 310
      if ( h .lt. n ) go to 300
      f = 1
      go to 10
350   j = -h
360   t = (1.0E+00,0.0E+00)
      i = h + 1
      e = i + 1
      d = e + 1
      h = h + m
      g = j + o
      w(i) = (0.0E+00,0.0E+00)
      w(e) = (0.0E+00,0.0E+00)
      w(d) = (0.0E+00,0.0E+00)
      w(h) = (0.0E+00,0.0E+00)
370   w(i) = w(i) + t*a(j+i)
      w(e) = w(e) + t*a(j+e)
      w(d) = w(d) + t*a(j+d)
      w(h) = w(h) + t*a(j+h)
      t = t*s
      j = j + m
      if ( j .lt. g ) go to 370
      j = j - m
      s = s*v
      if ( j+h .lt. n ) go to 360
      if ( h .lt. n ) go to 350
      f = 1
      go to 10
380   j = -h
390   t = (1.0E+00,0.0E+00)
      i = h + 1
      e = i + 1
      h = h + m
      g = j + o
      w(i) = (0.0E+00,0.0E+00)
      w(e) = (0.0E+00,0.0E+00)
      w(h) = (0.0E+00,0.0E+00)
400   w(i) = w(i) + t*a(j+i)
      w(e) = w(e) + t*a(j+e)
      w(h) = w(h) + t*a(j+h)
      t = t*s
      j = j + m
      if ( j .lt. g ) go to 400
      j = j - m
      s = s*v
      if ( j+h .lt. n ) go to 390
      if ( h .lt. n ) go to 380
      f = 1
      go to 10
410   j = -h
420   t = (1.0E+00,0.0E+00)
      i = h + 1
      h = h + m
      g = j + o
      w(i) = (0.0E+00,0.0E+00)
      w(h) = (0.0E+00,0.0E+00)
430   w(i) = w(i) + t*a(j+i)
      w(h) = w(h) + t*a(j+h)
      t = t*s
      j = j + m
      if ( j .lt. g ) go to 430
      j = j - m
      s = s*v
      if ( j+h .lt. n ) go to 420
      if ( h .lt. n ) go to 410
      f = 1
      go to 10
440   j = -h
450   t = (1.0E+00,0.0E+00)
      i = h + 1
      h = h + m
      g = j + o
      w(i) = (0.0E+00,0.0E+00)
460   w(i) = w(i) + t*a(j+i)
      t = t*s
      j = j + m
      if ( j .lt. g ) go to 460
      j = j - m
      s = s*v
      if ( j+h .lt. n ) go to 450
      if ( h .lt. n ) go to 440
      f = 1
      go to 10
470   if ( l .eq. 2 ) go to 480
      if ( l .eq. 3 ) go to 600
      go to 720
480   go to (580,560,540,520),m
490   j = -h
500   i = h + 1
      h = h + m
      e = j + m
      do k = i, h
        a(k) = w(j+k) + s*w(e+k)
      end do
      j = e
      s = s*v
      if ( j+h .lt. n ) go to 500
      if ( h .lt. n ) go to 490
      f = 0
      go to 10
520   j = -h
530   h = h + 1
      e = j + m
      a(h) = w(j+h) + s*w(e+h)
      h = h + 1
      a(h) = w(j+h) + s*w(e+h)
      h = h + 1
      a(h) = w(j+h) + s*w(e+h)
      h = h + 1
      a(h) = w(j+h) + s*w(e+h)
      j = e
      s = s*v
      if ( j+h .lt. n ) go to 530
      if ( h .lt. n ) go to 520
      f = 0
      go to 10
540   j = -h
550   h = h + 1
      e = j + m
      a(h) = w(j+h) + s*w(e+h)
      h = h + 1
      a(h) = w(j+h) + s*w(e+h)
      h = h + 1
      a(h) = w(j+h) + s*w(e+h)
      j = e
      s = s*v
      if ( j+h .lt. n ) go to 550
      if ( h .lt. n ) go to 540
      f = 0
      go to 10
560   j = -h
570   h = h + 1
      e = j + m
      a(h) = w(j+h) + s*w(e+h)
      h = h + 1
      a(h) = w(j+h) + s*w(e+h)
      j = e
      s = s*v
      if ( j+h .lt. n ) go to 570
      if ( h .lt. n ) go to 560
      f = 0
      go to 10
580   j = -h
590   h = h + 1
      e = j + m
      a(h) = w(j+h) + s*w(e+h)
      j = e
      s = s*v
      if ( j+h .lt. n ) go to 590
      if ( h .lt. n ) go to 580
      f = 0
      go to 10
600   go to (700,680,660,640),m
610   j = -h
620   i = h + 1
      h = h + m
      e = j + m
      d = e + m
      t = s*s
      do k = i, h
        a(k) = w(j+k) + s*w(e+k) + t*w(d+k)
      end do
      j = d
      s = s*v
      if ( j+h .lt. n ) go to 620
      if ( h .lt. n ) go to 610
      f = 0
      go to 10
640   j = -h
650   h = h + 1
      e = j + m
      d = e + m
      t = s*s
      a(h) = w(j+h) + s*w(e+h) + t*w(d+h)
      h = h + 1
      a(h) = w(j+h) + s*w(e+h) + t*w(d+h)
      h = h + 1
      a(h) = w(j+h) + s*w(e+h) + t*w(d+h)
      h = h + 1
      a(h) = w(j+h) + s*w(e+h) + t*w(d+h)
      j = d
      s = s*v
      if ( j+h .lt. n ) go to 650
      if ( h .lt. n ) go to 640
      f = 0
      go to 10
660   j = -h
670   h = h + 1
      e = j + m
      d = e + m
      t = s*s
      a(h) = w(j+h) + s*w(e+h) + t*w(d+h)
      h = h + 1
      a(h) = w(j+h) + s*w(e+h) + t*w(d+h)
      h = h + 1
      a(h) = w(j+h) + s*w(e+h) + t*w(d+h)
      j = d
      s = s*v
      if ( j+h .lt. n ) go to 670
      if ( h .lt. n ) go to 660
      f = 0
      go to 10
680   j = -h
690   h = h + 1
      e = j + m
      d = e + m
      t = s*s
      a(h) = w(j+h) + s*w(e+h) + t*w(d+h)
      h = h + 1
      a(h) = w(j+h) + s*w(e+h) + t*w(d+h)
      j = d
      s = s*v
      if ( j+h .lt. n ) go to 690
      if ( h .lt. n ) go to 680
      f = 0
      go to 10
700   j = -h
710   h = h + 1
      e = j + m
      d = e + m
      t = s*s
      a(h) = w(j+h) + s*w(e+h) + t*w(d+h)
      j = d
      s = s*v
      if ( j+h .lt. n ) go to 710
      if ( h .lt. n ) go to 700
      f = 0
      go to 10
720   go to (870,840,810,780),m
730   j = -h
740   t = (1.0E+00,0.0E+00)
      i = h + 1
      h = h + m
      g = j + o
      do k = i, h
        a(k) =  (0.0E+00,0.0E+00)
      end do
760   do k = i, h
        a(k) = a(k) + t*w(j+k)
      end do
      t = t*s
      j = j + m
      if ( j .lt. g ) go to 760
      j = j - m
      s = s*v
      if ( j+h .lt. n ) go to 740
      if ( h .lt. n ) go to 730
      f = 0
      go to 10
780   j = -h
790   t = (1.0E+00,0.0E+00)
      i = h + 1
      e = i + 1
      d = e + 1
      h = h + m
      g = j + o
      a(i) = (0.,0.)
      a(e) = (0.,0.)
      a(d) = (0.,0.)
      a(h) = (0.,0.)
800   a(i) = a(i) + t*w(j+i)
      a(e) = a(e) + t*w(j+e)
      a(d) = a(d) + t*w(j+d)
      a(h) = a(h) + t*w(j+h)
      t = t*s
      j = j + m
      if ( j .lt. g ) go to 800
      j = j - m
      s = s*v
      if ( j+h .lt. n ) go to 790
      if ( h .lt. n ) go to 780
      f = 0
      go to 10
810   j = -h
820   t = (1.0E+00,0.0E+00)
      i = h + 1
      e = i + 1
      h = h + m
      g = j + o
      a(i) = (0.0E+00,0.0E+00)
      a(e) = (0.0E+00,0.0E+00)
      a(h) = (0.0E+00,0.0E+00)
830   a(i) = a(i) + t*w(j+i)
      a(e) = a(e) + t*w(j+e)
      a(h) = a(h) + t*w(j+h)
      t = t*s
      j = j + m
      if ( j .lt. g ) go to 830
      j = j - m
      s = s*v
      if ( j+h .lt. n ) go to 820
      if ( h .lt. n ) go to 810
      f = 0
      go to 10
840   j = -h
850   t = (1.0E+00,0.0E+00)
      i = h + 1
      h = h + m
      g = j + o
      a(i) = (0.0E+00,0.0E+00)
      a(h) = (0.0E+00,0.0E+00)
860   a(i) = a(i) + t*w(j+i)
      a(h) = a(h) + t*w(j+h)
      t = t*s
      j = j + m
      if ( j .lt. g ) go to 860
      j = j - m
      s = s*v
      if ( j+h .lt. n ) go to 850
      if ( h .lt. n ) go to 840
      f = 0
      go to 10
870   j = -h
880   t = (1.0E+00,0.0E+00)
      i = h + 1
      h = h + m
      g = j + o
      a(i) = (0.0E+00,0.0E+00)
890   a(i) = a(i) + t*w(j+i)
      t = t*s
      j = j + m
      if ( j .lt. g ) go to 890
      j = j - m
      s = s*v
      if ( j+h .lt. n ) go to 880
      if ( h .lt. n ) go to 870
      f = 0
      go to 10
900   if ( f .eq. 0 ) return

      do i = 1, n
        a(i) = w(i)
      end do

      return
      end
      subroutine ffc0 ( a, n, w )

c*********************************************************************72
c
cc FFC0
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      complex a(*)
      integer f
      integer g
      integer h
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      integer np
      integer o
      integer p(25)
      complex s
      complex t
      complex u
      complex v
      complex w(*)

      save np
      save p

      data p /2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,
     & 73,79,83,89,97/

      data np /25/

      m = n
      f = 0
      u = cmplx(0.0E+00, -2.0E+00 * acos ( -1.0E+00 ) / n )
10    if ( m .eq. 1 ) go to 150
      do i = 1, np
        if ( (m/p(i))*p(i) .eq. m ) go to 30
      end do
      l = m
      go to 40
30    l = p(i)
40    o = m
      m = m/l
      v = cexp ( m * u )
      s = (1.0E+00,0.0E+00)
      h = 0
      if ( f .eq. 1 ) go to 100
50    j = -h
60    t = (1.0E+00,0.0E+00)
      i = h + 1
      h = h + m
      g = j + o
      do k = i, h
        w(k) =  (0.0E+00,0.0E+00)
      end do
80    do k = i, h
        w(k) = w(k) + t*a(j+k)
      end do
      t = t*s
      j = j + m
      if ( j .lt. g ) go to 80
      j = j - m
      s = s*v
      if ( j+h .lt. n ) go to 60
      if ( h .lt. n ) go to 50
      f = 1
      go to 10
100   j = -h
110   t = (1.0E+00,0.0E+00)
      i = h + 1
      h = h + m
      g = j + o
      do k = i, h
        a(k) =  (0.0E+00,0.0E+00)
      end do

130   do k = i, h
        a(k) = a(k) + t*w(j+k)
      end do
      t = t*s
      j = j + m
      if ( j .lt. g ) go to 130
      j = j - m
      s = s*v
      if ( j+h .lt. n ) go to 110
      if ( h .lt. n ) go to 100
      f = 0
      go to 10
150   if ( f .eq. 0 ) return

      do i = 1, n
        a(i) = w(i)
      end do

      return
      end
      subroutine fft ( a, n, w )

c*********************************************************************72
c
cc FFT computes the Fast Fourier Transform.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      integer np
      parameter ( np = 25 )

      complex a(*)
      integer d
      integer e
      integer f
      integer g
      integer h
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      integer o
      integer p(np)
      complex s
      complex t
      complex u
      complex v
      complex w(*)

      save p

      data p / 
     &   2,  3,  5,  7, 11, 13, 17, 19, 23, 29,
     &  31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
     &  73, 79, 83, 89, 97 /

      m = n
      f = 0
      u = cmplx(0.0E+00, 2.0E+00 * acos ( -1.0E+00 ) / n )

10    continue
c
c  Done when M is divided down to 1.
c
      if ( m .eq. 1 ) then

        if ( f .ne. 0 ) then
          do i = 1, n
            a(i) = w(i)
          end do
        end if

        return

      end if
c
c  Search for small divisors of M.
c
      l = m
      do i = 1, np
        if ( ( m / p(i) ) * p(i) .eq. m ) then
          l = p(i)
          go to 30
        end if
      end do

30    continue

      o = m
      m = m / l
      v = cexp ( m * u )
      s = (1.0E+00,0.0E+00)
      h = 0
      if ( f .eq. 1 ) go to 470
      if ( l .eq. 2 ) go to 50
      if ( l .eq. 3 ) go to 170
      go to 290
50    go to (150,130,110,90),m
60    j = -h
70    i = h + 1
      h = h + m
      e = j + m
      do k = i, h
        w(k) = a(j+k) + s*a(e+k)
      end do
      j = e
      s = s*v
      if ( j+h .lt. n ) go to 70
      if ( h .lt. n ) go to 60
      f = 1
      go to 10
90    j = -h
100   h = h + 1
      e = j + m
      w(h) = a(j+h) + s*a(e+h)
      h = h + 1
      w(h) = a(j+h) + s*a(e+h)
      h = h + 1
      w(h) = a(j+h) + s*a(e+h)
      h = h + 1
      w(h) = a(j+h) + s*a(e+h)
      j = e
      s = s*v
      if ( j+h .lt. n ) go to 100
      if ( h .lt. n ) go to 90
      f = 1
      go to 10
110   j = -h
120   h = h + 1
      e = j + m
      w(h) = a(j+h) + s*a(e+h)
      h = h + 1
      w(h) = a(j+h) + s*a(e+h)
      h = h + 1
      w(h) = a(j+h) + s*a(e+h)
      j = e
      s = s*v
      if ( j+h .lt. n ) go to 120
      if ( h .lt. n ) go to 110
      f = 1
      go to 10
130   j = -h
140   h = h + 1
      e = j + m
      w(h) = a(j+h) + s*a(e+h)
      h = h + 1
      w(h) = a(j+h) + s*a(e+h)
      j = e
      s = s*v
      if ( j+h .lt. n ) go to 140
      if ( h .lt. n ) go to 130
      f = 1
      go to 10
150   j = -h
160   h = h + 1
      e = j + m
      w(h) = a(j+h) + s*a(e+h)
      j = e
      s = s*v
      if ( j+h .lt. n ) go to 160
      if ( h .lt. n ) go to 150
      f = 1
      go to 10
170   go to (270,250,230,210),m
180   j = -h
190   i = h + 1
      h = h + m
      e = j + m
      d = e + m
      t = s*s
      do k = i, h
        w(k) = a(j+k) + s*a(e+k) + t*a(d+k)
      end do
      j = d
      s = s*v
      if ( j+h .lt. n ) go to 190
      if ( h .lt. n ) go to 180
      f = 1
      go to 10
210   j = -h
220   h = h + 1
      e = j + m
      d = e + m
      t = s*s
      w(h) = a(j+h) + s*a(e+h) + t*a(d+h)
      h = h + 1
      w(h) = a(j+h) + s*a(e+h) + t*a(d+h)
      h = h + 1
      w(h) = a(j+h) + s*a(e+h) + t*a(d+h)
      h = h + 1
      w(h) = a(j+h) + s*a(e+h) + t*a(d+h)
      j = d
      s = s*v
      if ( j+h .lt. n ) go to 220
      if ( h .lt. n ) go to 210
      f = 1
      go to 10
230   j = -h
240   h = h + 1
      e = j + m
      d = e + m
      t = s*s
      w(h) = a(j+h) + s*a(e+h) + t*a(d+h)
      h = h + 1
      w(h) = a(j+h) + s*a(e+h) + t*a(d+h)
      h = h + 1
      w(h) = a(j+h) + s*a(e+h) + t*a(d+h)
      j = d
      s = s*v
      if ( j+h .lt. n ) go to 240
      if ( h .lt. n ) go to 230
      f = 1
      go to 10
250   j = -h
260   h = h + 1
      e = j + m
      d = e + m
      t = s*s
      w(h) = a(j+h) + s*a(e+h) + t*a(d+h)
      h = h + 1
      w(h) = a(j+h) + s*a(e+h) + t*a(d+h)
      j = d
      s = s*v
      if ( j+h .lt. n ) go to 260
      if ( h .lt. n ) go to 250
      f = 1
      go to 10
270   j = -h
280   h = h + 1
      e = j + m
      d = e + m
      t = s*s
      w(h) = a(j+h) + s*a(e+h) + t*a(d+h)
      j = d
      s = s*v
      if ( j+h .lt. n ) go to 280
      if ( h .lt. n ) go to 270
      f = 1
      go to 10
290   go to (440,410,380,350),m
300   j = -h
310   t = (1.0E+00,0.0E+00)
      i = h + 1
      h = h + m
      g = j + o
      do k = i, h
        w(k) =  (0.0E+00,0.0E+00)
      end do

330   do k = i, h
        w(k) = w(k) + t*a(j+k)
      end do
      t = t*s
      j = j + m
      if ( j .lt. g ) go to 330
      j = j - m
      s = s*v
      if ( j+h .lt. n ) go to 310
      if ( h .lt. n ) go to 300
      f = 1
      go to 10
350   j = -h
360   t = (1.0E+00,0.0E+00)
      i = h + 1
      e = i + 1
      d = e + 1
      h = h + m
      g = j + o
      w(i) = (0.,0.)
      w(e) = (0.,0.)
      w(d) = (0.,0.)
      w(h) = (0.,0.)
370   w(i) = w(i) + t*a(j+i)
      w(e) = w(e) + t*a(j+e)
      w(d) = w(d) + t*a(j+d)
      w(h) = w(h) + t*a(j+h)
      t = t*s
      j = j + m
      if ( j .lt. g ) go to 370
      j = j - m
      s = s*v
      if ( j+h .lt. n ) go to 360
      if ( h .lt. n ) go to 350
      f = 1
      go to 10
380   j = -h
390   t = (1.0E+00,0.)
      i = h + 1
      e = i + 1
      h = h + m
      g = j + o
      w(i) = (0.,0.)
      w(e) = (0.,0.)
      w(h) = (0.,0.)
400   w(i) = w(i) + t*a(j+i)
      w(e) = w(e) + t*a(j+e)
      w(h) = w(h) + t*a(j+h)
      t = t*s
      j = j + m
      if ( j .lt. g ) go to 400
      j = j - m
      s = s*v
      if ( j+h .lt. n ) go to 390
      if ( h .lt. n ) go to 380
      f = 1
      go to 10
410   j = -h
420   t = (1.0E+00,0.)
      i = h + 1
      h = h + m
      g = j + o
      w(i) = (0.,0.)
      w(h) = (0.,0.)
430   w(i) = w(i) + t*a(j+i)
      w(h) = w(h) + t*a(j+h)
      t = t*s
      j = j + m
      if ( j .lt. g ) go to 430
      j = j - m
      s = s*v
      if ( j+h .lt. n ) go to 420
      if ( h .lt. n ) go to 410
      f = 1
      go to 10
440   j = -h
450   t = (1.0E+00,0.)
      i = h + 1
      h = h + m
      g = j + o
      w(i) = (0.,0.)
460   w(i) = w(i) + t*a(j+i)
      t = t*s
      j = j + m
      if ( j .lt. g ) go to 460
      j = j - m
      s = s*v
      if ( j+h .lt. n ) go to 450
      if ( h .lt. n ) go to 440
      f = 1
      go to 10
470   if ( l .eq. 2 ) go to 480
      if ( l .eq. 3 ) go to 600
      go to 720
480   go to (580,560,540,520),m
490   j = -h
500   i = h + 1
      h = h + m
      e = j + m
      do k = i, h
        a(k) = w(j+k) + s*w(e+k)
      end do
      j = e
      s = s*v
      if ( j+h .lt. n ) go to 500
      if ( h .lt. n ) go to 490
      f = 0
      go to 10
520   j = -h
530   h = h + 1
      e = j + m
      a(h) = w(j+h) + s*w(e+h)
      h = h + 1
      a(h) = w(j+h) + s*w(e+h)
      h = h + 1
      a(h) = w(j+h) + s*w(e+h)
      h = h + 1
      a(h) = w(j+h) + s*w(e+h)
      j = e
      s = s*v
      if ( j+h .lt. n ) go to 530
      if ( h .lt. n ) go to 520
      f = 0
      go to 10
540   j = -h
550   h = h + 1
      e = j + m
      a(h) = w(j+h) + s*w(e+h)
      h = h + 1
      a(h) = w(j+h) + s*w(e+h)
      h = h + 1
      a(h) = w(j+h) + s*w(e+h)
      j = e
      s = s*v
      if ( j+h .lt. n ) go to 550
      if ( h .lt. n ) go to 540
      f = 0
      go to 10
560   j = -h
570   h = h + 1
      e = j + m
      a(h) = w(j+h) + s*w(e+h)
      h = h + 1
      a(h) = w(j+h) + s*w(e+h)
      j = e
      s = s*v
      if ( j+h .lt. n ) go to 570
      if ( h .lt. n ) go to 560
      f = 0
      go to 10
580   j = -h
590   h = h + 1
      e = j + m
      a(h) = w(j+h) + s*w(e+h)
      j = e
      s = s*v
      if ( j+h .lt. n ) go to 590
      if ( h .lt. n ) go to 580
      f = 0
      go to 10
600   go to (700,680,660,640),m
610   j = -h
620   i = h + 1
      h = h + m
      e = j + m
      d = e + m
      t = s*s
      do k = i, h
        a(k) = w(j+k) + s*w(e+k) + t*w(d+k)
      end do
      j = d
      s = s*v
      if ( j+h .lt. n ) go to 620
      if ( h .lt. n ) go to 610
      f = 0
      go to 10
640   j = -h
650   h = h + 1
      e = j + m
      d = e + m
      t = s*s
      a(h) = w(j+h) + s*w(e+h) + t*w(d+h)
      h = h + 1
      a(h) = w(j+h) + s*w(e+h) + t*w(d+h)
      h = h + 1
      a(h) = w(j+h) + s*w(e+h) + t*w(d+h)
      h = h + 1
      a(h) = w(j+h) + s*w(e+h) + t*w(d+h)
      j = d
      s = s*v
      if ( j+h .lt. n ) go to 650
      if ( h .lt. n ) go to 640
      f = 0
      go to 10
660   j = -h
670   h = h + 1
      e = j + m
      d = e + m
      t = s*s
      a(h) = w(j+h) + s*w(e+h) + t*w(d+h)
      h = h + 1
      a(h) = w(j+h) + s*w(e+h) + t*w(d+h)
      h = h + 1
      a(h) = w(j+h) + s*w(e+h) + t*w(d+h)
      j = d
      s = s*v
      if ( j+h .lt. n ) go to 670
      if ( h .lt. n ) go to 660
      f = 0
      go to 10
680   j = -h
690   h = h + 1
      e = j + m
      d = e + m
      t = s*s
      a(h) = w(j+h) + s*w(e+h) + t*w(d+h)
      h = h + 1
      a(h) = w(j+h) + s*w(e+h) + t*w(d+h)
      j = d
      s = s*v
      if ( j+h .lt. n ) go to 690
      if ( h .lt. n ) go to 680
      f = 0
      go to 10
700   j = -h
710   h = h + 1
      e = j + m
      d = e + m
      t = s*s
      a(h) = w(j+h) + s*w(e+h) + t*w(d+h)
      j = d
      s = s*v
      if ( j+h .lt. n ) go to 710
      if ( h .lt. n ) go to 700
      f = 0
      go to 10
720   go to (870,840,810,780),m
730   j = -h
740   t = (1.0E+00,0.)
      i = h + 1
      h = h + m
      g = j + o
      do k = i, h
        a(k) =  (0.,0.)
      end do
760   do k = i, h
        a(k) = a(k) + t*w(j+k)
      end do
      t = t*s
      j = j + m
      if ( j .lt. g ) go to 760
      j = j - m
      s = s*v
      if ( j+h .lt. n ) go to 740
      if ( h .lt. n ) go to 730
      f = 0
      go to 10
780   j = -h
790   t = (1.0E+00,0.)
      i = h + 1
      e = i + 1
      d = e + 1
      h = h + m
      g = j + o
      a(i) = (0.,0.)
      a(e) = (0.,0.)
      a(d) = (0.,0.)
      a(h) = (0.,0.)
800   a(i) = a(i) + t*w(j+i)
      a(e) = a(e) + t*w(j+e)
      a(d) = a(d) + t*w(j+d)
      a(h) = a(h) + t*w(j+h)
      t = t*s
      j = j + m
      if ( j .lt. g ) go to 800
      j = j - m
      s = s*v
      if ( j+h .lt. n ) go to 790
      if ( h .lt. n ) go to 780
      f = 0
      go to 10
810   j = -h
820   t = (1.0E+00,0.)
      i = h + 1
      e = i + 1
      h = h + m
      g = j + o
      a(i) = (0.,0.)
      a(e) = (0.,0.)
      a(h) = (0.,0.)
830   a(i) = a(i) + t*w(j+i)
      a(e) = a(e) + t*w(j+e)
      a(h) = a(h) + t*w(j+h)
      t = t*s
      j = j + m
      if ( j .lt. g ) go to 830
      j = j - m
      s = s*v
      if ( j+h .lt. n ) go to 820
      if ( h .lt. n ) go to 810
      f = 0
      go to 10
840   j = -h
850   t = (1.0E+00,0.)
      i = h + 1
      h = h + m
      g = j + o
      a(i) = (0.,0.)
      a(h) = (0.,0.)
860   a(i) = a(i) + t*w(j+i)
      a(h) = a(h) + t*w(j+h)
      t = t*s
      j = j + m
      if ( j .lt. g ) go to 860
      j = j - m
      s = s*v
      if ( j+h .lt. n ) go to 850
      if ( h .lt. n ) go to 840
      f = 0
      go to 10
870   j = -h
880   t = (1.0E+00,0.)
      i = h + 1
      h = h + m
      g = j + o
      a(i) = (0.,0.)
890   a(i) = a(i) + t*w(j+i)
      t = t*s
      j = j + m
      if ( j .lt. g ) go to 890
      j = j - m
      s = s*v
      if ( j+h .lt. n ) go to 880
      if ( h .lt. n ) go to 870
      f = 0
      go to 10

      end
      subroutine fft0 ( a, n, w )

c*********************************************************************72
c
cc FFT0
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      integer np
      parameter ( np = 25 )

      complex a(*)
      integer f
      integer g
      integer h
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      integer o
      integer p(np)
      complex s
      complex t
      complex u
      complex v
      complex w(*)

      save p

      data p / 2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,
     & 73,79,83,89,97/

      m = n
      f = 0
      u = cmplx(0., 2. * acos ( -1.0E+00 ) / n )

10    continue
c
c  Done when M is divided down to 1.
c
      if ( m .eq. 1 ) then

        if ( f .ne. 0 ) then
          do i = 1, n
            a(i) = w(i)
          end do
        end if

        return

      end if
c
c  Search for small divisors of M.
c
      l = m
      do i = 1, np
        if ( (m/p(i))*p(i) .eq. m ) then
          l = p(i)
          go to 30
        end if
      end do

30    continue

      o = m
      m = m/l
      v = cexp ( m * u )
      s = (1.0E+00,0.0E+00)
      h = 0
      if ( f .eq. 1 ) go to 100
50    j = -h
60    t = (1.0E+00,0.)
      i = h + 1
      h = h + m
      g = j + o
      do k = i, h
        w(k) =  (0.0, 0.0 )
      end do
80    continue

      do k = i, h
        w(k) = w(k) + t*a(j+k)
      end do

      t = t*s
      j = j + m
      if ( j .lt. g ) go to 80
      j = j - m
      s = s*v
      if ( j+h .lt. n ) go to 60
      if ( h .lt. n ) go to 50
      f = 1
      go to 10
100   j = -h
110   t = (1.0E+00,0.0E+00)
      i = h + 1
      h = h + m
      g = j + o

      do k = i, h
        a(k) =  (0.0E+00,0.0E+00)
      end do

130   continue

      do k = i, h
        a(k) = a(k) + t*w(j+k)
      end do

      t = t*s
      j = j + m
      if ( j .lt. g ) go to 130
      j = j - m
      s = s*v
      if ( j+h .lt. n ) go to 110
      if ( h .lt. n ) go to 100
      f = 0
      go to 10

      return
      end
      subroutine fgv ( x, y, s, p, q, a, b )

c*********************************************************************72
c
cc FGV
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real a
      real b
      real c
      real p
      real q
      real r
      real s
      real t
      real x
      real y

      if ( abs ( p ) .gt. abs ( q ) ) go to 10
      if ( q .eq. 0.0E+00 ) go to 110
      r = a / b
      s = p / q
      t = abs ( r ) * s * s
      if ( t .lt. 1 ) go to 70
      t = t / (1.0E+00+t)
      r = b / a
      s = q / p
      go to 20
10    r = b / a
      s = q / p
      t = abs ( r ) * s * s
      if ( t .gt. 1.0E+00 ) go to 60
      t = 1.0E+00 / (1.0E+00+t)
20    continue

      a = sign ( a * t, p )

      if ( sign ( r, p ) .ne. r ) then
        b = - abs ( b * t )
      else
        b = abs ( b * t )
      end if

      y = s
      x = abs ( r ) * s
      s = 0.0E+00
      return
60    t = t / (1.0E+00+t)
      r = a / b
      s = p / q
      go to 80
70    t = 1.0E+00 / (1.0E+00+t)
80    c = a
      a = sign ( b * t, q )
      if ( sign ( r, q ) .eq. r ) go to 90
      b = - abs ( c * t )
      go to 100
90    b = abs ( c * t )
100   y = s
      x = abs ( r ) * s
      s = 1.0E+00
      return
110   x = 0.0E+00
      y = 0.0E+00
      s = 0.0E+00
      return
      end
      function fv ( a, x, h, n, value )

c*********************************************************************72
c
cc FV
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      integer n

      real a
      real fv
      real h(n,*)
      integer i
      external value
      real value
      real x(n)

      do i = 1, n
        h(i,2) = x(i) + a * h(i,1)
      end do

      fv = value ( h(1,2) )

      return
      end
      subroutine fvd ( v, d, a, x, h, n, both )

c*********************************************************************72
c
cc FVD
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      integer n

      real a
      external both
      real d
      real h(n,*)
      integer i
      real v
      real x(n)

      do i = 1, n
        h(i,2) = x(i) + a * h(i,1)
      end do

      call both ( v, h(1,3), h(1,2) )

      d = 0.0E+00
      do i = 1, n
        d = d + h(i,1) * h(i,3)
      end do

      return
      end
      function hcon ( a, b )

c*********************************************************************72
c
cc HCON estimates the condition number of a symmetric band matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c    William Hager,
c    Condition Estimates,
c    SIAM Journal on Scientific and Statistical Computing,
c    Volume 5, Number 2, June 1984, pages 311-316.
c
c  Parameters:
c
c    Input, real A(4+(H+1)*N), the factorization information from HFACT.
c
c    Workspace, real B(N).
c
c    Output, real HCON, the condition number of the matrix.
c
      implicit none

      real a(*)
      real b(*)
      real c
      real d
      real hcon
      integer i
      integer j
      integer m
      integer n

      c = a(1)

      if ( abs ( c ) .ne. 1232 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor with hfact before',
     &  'estimating condition'
        stop
      end if

      if ( c .lt. 0.0E+00 ) then
        hcon = -1.0E+00
        return
      end if

      m = 0
      n = a(2)
      j = 5 + n + n * a(4)
      c = 1.0E+00 /a(2)
      do j = 1, n
        b(j) = c
      end do
      go to 50
30    do j = 1, n
        b(j) = 0.0E+00
      end do
      b(m) = 1.0E+00
50    call hsolve ( b, a, b )

      c = 0.0E+00
      do j = 1, n
        c = c + abs ( b(j) )
        if ( b(j) .lt. 0.0E+00 ) then
          b(j) = -1.0E+00
        else
          b(j) = 1.0E+00
        end if
      end do

      call hsolve ( b, a, b )
c
c  I is the index of the largest entry of B.
c
      i = 1
      do j = 1, n
        if ( abs ( b(i) ) .lt. abs ( b(j) ) ) then
          i = j
        end if
      end do

      if ( m .eq. 0 ) go to 90
      if ( m .eq. i ) go to 100
      if ( d .ge. c ) go to 100
90    m = i
      d = c
      go to 30
100   c = c * a(3)
      c = max ( c, 1.0E+00 )
      hcon = c
      return
      end
      function hdet ( iexp, a )

c*********************************************************************72
c
cc HDET computes the determinant of a symmetric band matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, integer IEXP, a power of ten that is part of the
c    determinant.
c
c    Input, real A(4+(H+1)*N), factorization information from HFACT.
c
c    Output, real HDET, the mantissa of the determinant.
c    Determinant = HDET * 10.0^IEXP.
c
      implicit none

      real a(*)
      real c
      real d
      real f
      real g
      real hdet
      integer h
      integer iexp
      integer i
      integer k
      integer l
      integer m
      integer n

      iexp=0
      hdet=0.0E+00
      d = a(1)

      if ( abs ( d ) .ne. 1232 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor with hfact before',
     &  ' computing determinant'
        return
      end if

      iexp = 0
      if ( d .lt. 0.0E+00 ) go to 70
      d = 1.0E+00
      f = 2.0E+00**64
      g = 1.0E+00 / f
      h = 64
      n = a(2)
      m = a(4) + 1
      k = 5
      l = 4 + m*n

      do i = k, l, m

        d = d * a(i)
20      if ( abs ( d ) .lt. f ) go to 30
        iexp = iexp + h
        d = d*g
        go to 20

30      if ( abs ( d ) .le. g ) then
          iexp = iexp - h
          d = d * f
          go to 30
        end if

      end do

      if ( iexp .eq. 0 ) then
        hdet = d
        return
      end if

      if ( d .eq. 0.0E+00 ) go to 80
      c = alog10 ( abs ( d ) ) + iexp * alog10 ( 2.0E+00 )
      iexp = c
      c = c - iexp
      if ( c .le. 0.0E+00 ) go to 60
      c = c - 1
      iexp = iexp + 1
60    f = 10.0E+00**c
      if ( d .lt. 0.0E+00 ) f = -f
      hdet = f
      return
70    hdet = 0.0E+00
      return
80    iexp = 0
      go to 70
      end
      subroutine hdiag ( e, v, lv, a, la, n, h, w )

c*********************************************************************72
c
cc HDIAG diagonalizes a symmetric band matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input, integer LV, the leading dimension of the array V.
c
      implicit none

      integer la
      integer lv

      real a(la,*)
      real e(*)
      integer h
      integer i
      integer j
      integer m
      integer n
      real v(lv,*)
      real w(*)

      if ( h .le. 0 ) then

        do j = 1, n
          do i = 1, n
            v(i,j) = 0.0E+00
          end do
          v(j,j) = 1.0E+00
          e(j) = a(1,j)
        end do

      else

        call hsim ( v, lv, e, v, a, la, n, h, w )

        do i = 1, n - 1
          w(i) = e(i)
          w(i+n) = v(i,1)
          v(i,1) = 0.0E+00
        end do

        v(n,1) = 0.0E+00
        v(1,1) = 1.0E+00
        w(n) = e(n)
        m = n + 1
        call tdg ( e, v, lv, w, w(m), n )

      end if

      return
      end
      subroutine hess ( a, la, n, w )

c*********************************************************************72
c
cc HESS reduces a real matrix to Hessenberg form.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real a(*)
      integer c
      integer d
      integer e
      integer g
      integer h
      integer i
      integer j
      integer k
      integer l
      integer la
      integer m
      integer n
      integer o
      integer p
      integer q
      real r
      real s
      real t
      real u
      real w(*)

      if ( la .gt. n ) then
        call pack ( a, la, n )
      end if

      i = n*n
      m = n + 1
      o = m + 1
      e = i + m
      j = m
      k = i

10    k = k - n

20    a(i+j) = a(i)
      i = i - 1
      if ( i .gt. k ) go to 20
      j = j - 1
      if ( k .gt. 0 ) go to 10
      a(1) = 2230
      a(2) = n
      k = 4
      l = o
      d = 1
      c = 2

30    if ( c .ge. n ) return
      p = k + 1
      do i = p, l
        if ( a(i) .ne. 0.0E+00 ) go to 50
      end do
      a(l+1) = 0.0E+00
      go to 180

50    t = abs ( a(k) )
      if ( t .ne. 0.0E+00 ) then
        u = 1.0E+00 /t
      end if
      r = 1.0E+00

      do j = i, l
        s = abs ( a(j) )
        if ( s .gt. t ) then
          u = 1.0E+00 /s
          r = 1.0E+00 + r*(t*u)**2
          t = s
        else
          r = r + (s*u)**2
        end if
      end do

      s = t * sqrt ( r )
      r = a(k)
      u = 1.0E+00 / sqrt ( s * ( s + abs ( r ) ) )
      if ( r .lt. 0.0E+00 ) s = -s
      i = l
80    a(i+1) = u*a(i)
      i = i - 1
      if ( i .gt. k ) go to 80
      a(k) = -s
      a(p) = u*(r+s)
      h = l
      do i = 1, n
        w(i) = 0.0E+00
      end do
100   h = h + m
      s = a(p)
      p = p + 1
      q = h - n
      do i = 1, d
        w(i) = w(i) + s*a(i+q)
      end do
      j = k - d

      t = 0.0E+00
      do i = c, n
        r = a(i+q)
        t = t + r*a(i+j)
        w(i) = w(i) + r*s
      end do

      a(h+1) = t
      if ( h .lt. e ) go to 100
      t = 0.0E+00
      h = l + 1
      p = k + 1
      j = c - p
      do i = p, h
        t = t + w(i+j)*a(i)
      end do
      do i = c, n
        w(i) = w(i) - t*a(i-j)
      end do
      h = l
150   g = h + 2
      q = h + m
      h = h + c
      t = a(q+1)
      s = a(p)
      p = p + 1
      j = 1 - g
      do i = g, h
        a(i) = a(i) - w(i+j)*s
      end do
      i = h
      h = q
      q = k - i
      g = i + 1
      do i = g, h
        a(i) = a(i) - a(i+q)*t - w(i+j)*s
      end do
      if ( h .lt. e ) go to 150
180   k = k + o
      l = l + m
      d = c
      c = c + 1
      go to 30
      end
      subroutine hfact ( a, la, n, h )

c*********************************************************************72
c
cc HFACT factors a symmetric band matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input/output, real A(4+(H+1)*N); on input, the diagonal and subdiagonal bands
c    of the matrix.  On output, factorization information.
c
c    Input, integer LA, the leading dimension of A.
c
c    Input, integer N, the number of rows and columns in the matrix.
c
c    Input, integer H, the half bandwidth.
c
      implicit none

      real a(*)
      integer d
      integer e
      integer f
      integer g
      integer h
      integer i
      integer j
      integer k
      integer l
      integer la
      integer m
      integer n
      real q
      real r
      real s
      real t

      g = h + 1
      f = la*n + g
      k = h
10    if ( k .eq. 0 ) go to 30
      f = f - la
      i = f
      j = f - k
      k = k - 1
20    a(i) = 0.0E+00
      i = i - 1
      if ( i .gt. j ) go to 20
      go to 10
30    if ( la .gt. g ) then
        call rpack ( a, la, g, n )
      end if
      f = 2 - g + g*n
      a(f+h) = 0.0E+00
      q = 0.0E+00
      r = q
      i = 0
      j = 0
      l = 0
40    l = l + 1
      if ( l .gt. n ) go to 80
      d = j + g
      j = j + 1
      e = j + min ( h, n - l )
      s = q + abs ( a(j) )
      k = f
50    if ( j .eq. e ) go to 60
      j = j + 1
      t = abs ( a(j) )
      s = s + t
      a(k) = a(k+1) + t
      k = k + 1
      go to 50
60    q = a(f)
      if ( r .lt. s ) r = s
70    if ( j .eq. d ) go to 40
      j = j + 1
      a(j) = 0.0E+00
      go to 70
80    j = 4 + g*n
90    a(j) = a(j-4)
      j = j - 1
      if ( j .gt. 4 ) go to 90
      a(1) = 1232
      a(2) = n
      a(3) = r
      a(4) = h
      i = 5 - g
      k = 0
100   k = k + 1
      i = i + g
      if ( k .eq. n ) go to 140
      m = min ( h, n - k )
      s = a(i)
      if ( s .eq. 0.0E+00 ) go to 130
      d = 0
      e = i
      j = i
110   if ( m .eq. 0 ) go to 100
      j = j + 1
      t = a(j)/s
      d = d - h
      e = e + g
      m = m - 1
      if ( t .eq. 0.0E+00 ) go to 110
      f = e + m
      do l = e, f
        a(l) = a(l) - t*a(d+l)
      end do
      go to 110
130   a(1) = -1232
      go to 100
140   if ( a(i) .ne. 0.0E+00 ) return
      a(1) = -1232
      return
      end
      subroutine hhess ( d, u, a, la, n, h, w )

c*********************************************************************72
c
cc HHESS reduces a real symmetric matrix to tridiagonal form.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      integer la

      real a(la,*)
      real d(*)
      integer g
      integer h
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      integer o
      real u(*)
      real w(*)

      m = n + 1
      k = n - 1

      if ( h .le. 0 ) then

        d(n) = a(1,n)

        if ( n .eq. 1 ) then
          return
        end if

        do j = 1, k
          d(j) = a(1,j)
          u(j) = 0.0E+00
        end do

        return

      end if

      if ( h .le. 1 ) then

        do j = 1, k
          d(j) = a(1,j)
          u(j) = a(2,j)
        end do
        d(n) = a(1,n)

        return

      end if

      g = h + 1

      k = 0
      do j = 1, n
        l = min ( g, m - j )
        do i = 1, l
          w(i+k) = a(i,j)
        end do
        k = k + g
      end do

      l = n * g + 1
      do i =  1, n
        l = l - g
        k = 0
        o = l + min ( n - i, h )
        do j = l, o
          w(j) = w(j-k)
          k = k + g
        end do
      end do

      j = 0
      k = 0
      l = 0
      do i = 1, n
        m = l + 1
        o = min ( i, g )
        l = l + o
        do j = m, l
          w(j) = w(j+k)
        end do
        k = k + g - o
      end do

      g = h - 1
      i = l + 1
      j = i + g
      k = j + g
      l = k + g
      m = l + g
      g = m + g
      call htr ( w, w(i), w(j), w(k), w(l), w(m), w(g), d, u, n, h )

      return
      end
      subroutine hmult ( y, x, a, la, n, h )

c*********************************************************************72
c
cc HMULT multiplies a symmetric band matrix times a vector.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, real Y(N), the product of the matrix times X.
c
c    Input, real X(N), the vector to be multiplied.
c
c    Input, real A(LA,N), the (H+1)*N matrix of bands.
c
c    Input, integer LA, the leading dimension of A.
c
c    Input, integer N, the number of rows and columns in the matrix.
c
c    Input, integer H, the number of bands.
c
      implicit none

      integer la
      integer n

      real a(la,n)
      integer h
      integer i
      integer j
      integer k
      integer l
      integer m
      real r
      real s
      real t
      real x(n)
      real y(n)

      if ( h .eq. 0 ) then

        do i = 1, n
          y(i) = x(i) * a(1,i)
        end do

        return

      end if

      do i = 1, n
        y(i) = 0.0E+00
      end do

      k = 1

10    continue

      j = k
      k = k + 1
      s = x(j)
      t = a(1,j) * s

      if ( k .gt. n ) then
        y(n) = y(n) + t
        return
      end if

      m = 2 - k
      l = min ( j + h, n )

      do i = k, l
        r = a(i+m,j)
        t = t + r * x(i)
        y(i) = y(i) + s * r
      end do

      y(j) = y(j) + t

      go to 10
      end
      subroutine hse ( y, w, k, l, n )

c*********************************************************************72
c
cc HSE
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      integer n

      integer i
      integer j
      integer k
      integer l
      real s
      real t
      real w(*)
      real y(n)

      j = k*n - (k*(1+k))/2

      t = 0.0E+00
      do i = l, n
        t = t + y(i)**2
      end do
      t = sqrt ( t )

      s = y(l)
      if ( s .lt. 0.0E+00 ) s = s - t
      if ( s .ge. 0.0E+00 ) s = s + t
      t = sqrt ( abs ( s * t ) )
      if ( t .ne. 0.0E+00 ) t = 1.0E+00 /t
      do i = l, n
        y(i) = t*y(i)
        w(i+j) = y(i)
      end do
      y(l) = s*t
      w(l+j) = y(l)

      return
      end
      subroutine hsim ( p, lp, d, u, a, la, n, h, w )

c*********************************************************************72
c
cc HSIM reduces a real symmetric band matrix to tridiagonal form.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      integer la
      integer lp

      real a(la,*)
      real d(*)
      integer g
      integer h
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      integer o
      real p(lp,*)
      real u(*)
      real w(*)

      m = n + 1
      k = n - 1
      if ( h .gt. 0 ) go to 50
      if ( n .eq. 1 ) go to 40

      do j = 1, k
        do i = 1, n
          p(i,j) = 0.0E+00
        end do
        p(j,j) = 1.0E+00
        d(j) = a(1,j)
        u(j) = 0.0E+00
      end do

      do i = 1, n
        p(i,n) = 0.0E+00
      end do

40    p(n,n) = 1.0E+00
      d(n) = a(1,n)
      return
50    if ( h .gt. 1 ) go to 90

      do j = 1, k
        do i = 1, n
          p(i,j) = 0.0E+00
        end do
        p(j,j) = 1.0E+00
        d(j) = a(1,j)
        u(j) = a(2,j)
      end do

      do i = 1, n
        p(i,n) = 0.0E+00
      end do

      p(n,n) = 1.0E+00
      d(n) = a(1,n)
      return
90    g = h + 1

      k = 0
      do j = 1, n
        l = min ( g, m - j )
        do i = 1, l 
          w(i+k) = a(i,j)
        end do
        do i = 1, n
          p(i,j) = 0.0E+00
        end do
        p(j,j) = 1.0E+00
        k = k + g
      end do

      l = n*g + 1
      do i =  1, n
        l = l - g
        k = 0
        o = l + min ( n - i, h )
        do j = l, o
          w(j) = w(j-k)
          k = k + g
        end do
      end do

      j = 0
      k = 0
      l = 0
      do i = 1, n
        m = l + 1
        o = min ( i, g )
        l = l + o
        do j = m, l
          w(j) = w(j+k)
        end do
        k = k + g - o
      end do

      g = h - 1
      i = l + 1
      j = i + g
      k = j + g
      l = k + g
      m = l + g
      g = m + g

      call htr2 ( w, w(i), w(j), w(k), w(l), w(m), w(g), d, u, n, h, 
     &  p, lp )

      return
      end
      subroutine hsolve ( x, a, b )

c*********************************************************************72
c
cc HSOLVE solves a symmetric band matrix linear system.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, real X(N), the solution.
c
c    Input, real A(4+(H+1)*N), factorization information from HFACT.
c
c    Input, real B(N), the right hand side.
c
      implicit none

      real a(*)
      real b(*)
      integer h
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      real t
      real x(*)

      t = a(1)

      if ( abs ( t ) .ne. 1232 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor with hfact before solving'
        stop
      end if

      n = a(2)
      h = a(4)
      j = 4 - h
      k = 1
      if ( t .lt. 0.0E+00 ) go to 90

      do i = 1, n
        x(i) = b(i)
      end do

      if ( h .gt. 0 ) go to 40

      do k = 1, n
        x(k) = x(k) / a(4+k)
      end do

      return
40    j = j + h
      if ( k .eq. n ) go to 60
      t = x(k)/a(j+k)
      l = min ( n, k + h )
      k = k + 1
      if ( t .eq. 0.0E+00 ) go to 40
      do i = k, l
        x(i) = x(i) - t*a(i+j)
      end do
      go to 40
60    x(n) = x(n)/a(j+k)
70    if ( k .eq. 1 ) return
      j = j - h
      m = k
      k = k - 1
      l = min ( n, k + h )

      t = x(k)
      do i = m, l
        t = t - x(i)*a(i+j)
      end do
      x(k) = t/a(j+k)

      go to 70
90    k = 0
100   k = k + 1
      j = j + h
      if ( a(j+k) .ne. 0.0E+00 ) go to 100

      do i = 1, n
        x(i) = 0.0E+00
      end do
      x(k) = 1.0E+00
      if ( h .eq. 0 ) return
      go to 70
      end
      subroutine hsr1 ( a, la, n )

c*********************************************************************72
c
cc HSR1
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      integer la

      real a(la,*)
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      real s

      a(1,1) = 1.0E+00
      if ( n .eq. 1 ) then
        return
      end if

      a(1,2) = 0.0E+00

      if ( n .eq. 2 ) then
        a(2,1) = 0.0E+00
        a(2,2) = 1.0E+00
        return
      end if

      l = n - 2
      m = n - 1
      k = n
      s = a(k,l)
      a(n,n) = 1.0E+00 - s*a(n,l)
      a(m,n) = -s*a(m,l)
20    j = m
      m = l
      l = l - 1
      if ( l .eq. 0 ) go to 50
      s = 0.0E+00
      do i = j, n
        s = s + a(i,l)*a(i,k)
      end do
      a(m,k) = -s*a(m,l)
      do i = j, n
        a(i,k) = a(i,k) - s*a(i,l)
      end do
      go to 20
50    a(1,k) = 0.0E+00
      k = k - 1
      m = k
      l = k - 1
      s = -a(m,l)
      do i = k, n
         a(i,k) = s*a(i,l)
      end do
      a(k,k) = 1.0E+00 + a(k,k)
      if ( l .gt. 1 ) go to 20
      do i = 2, n
        a(i,1) = 0.0E+00
      end do

      return
      end
      subroutine hsr2 ( a, la, n )

c*********************************************************************72
c
cc HSR2
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real a(la,*)
      integer i
      integer j
      integer k
      integer l
      integer la
      integer m
      integer n
      real s

      if ( n .eq. 1 ) then
        a(1,1) = 1.0E+00
        return
      end if

      l = n - 2
      m = n - 1
      k = n
      s = a(n,m)
      a(n,n) = 1.0E+00 - s*a(n,m)
      a(m,n) = -s*a(m,m)
20    j = m
      m = m - 1
      if ( m .eq. 0 ) go to 50
      s = 0.0E+00
      do i = j, n
        s = s + a(i,m)*a(i,k)
      end do
      a(m,k) = -s*a(m,m)
      do i = j, n
        a(i,k) = a(i,k) - s*a(i,m)
      end do
      go to 20
50    k = k - 1
      m = k
      s = -a(k,k)
      do i = k, n
         a(i,k) = s*a(i,k)
      end do
      a(k,k) = 1.0E+00 + a(k,k)
      if ( k .gt. 1 ) go to 20
      return
      end
      subroutine hsr3 ( a, la, m, n )

c*********************************************************************72
c
cc HSR3
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      integer la

      real a(la,*)
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      real s

      k = n

      if ( m .lt. n ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: argument m must be .ge. n in routine hsr3'
        stop
      end if

10    continue

      s = -a(k,k)

      do i = k, m
        a(i,k) = s * a(i,k)
      end do

      a(k,k) = 1.0E+00 + a(k,k)
      if ( k .eq. 1 ) return
      l = k
30    j = l
      l = l - 1
      if ( l .eq. 0 ) go to 60
      s = 0.0E+00
      do i = j, m
        s = s + a(i,l)*a(i,k)
      end do
      a(l,k) = -s*a(l,l)
      do i = j, m
        a(i,k) = a(i,k) - s*a(i,l)
      end do
      go to 30
60    k = k - 1
      go to 10
      end
      subroutine hsr4 ( a, la, m, n )

c*********************************************************************72
c
cc HSR4
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      integer la

      real a(la,*)
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      real s

      k = m

      if ( m .le. n ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: argument m must be .gt. n in routine hsr4'
        stop
      end if

10    continue

      s = -a(k,n)
      do i = n, m
        a(i,k) = s*a(i,n)
      end do
      a(k,k) = 1.0E+00 + a(k,k)
      l = n
30    j = l
      l = l - 1
      if ( l .eq. 0 ) go to 60
      s = 0.0E+00
      do i = j, m
        s = s + a(i,l)*a(i,k)
      end do
      a(l,k) = -s*a(l,l)
      do i = j, m
        a(i,k) = a(i,k) - s*a(i,l)
      end do
      go to 30

60    continue

      k = k - 1
      if ( k .gt. n ) go to 10

      k = m - n

      do j = 1, k
        do i = 1, m
          a(i,j) = a(i,j+n)
        end do
      end do

      return
      end
      subroutine hsr5 ( a, la, m, n )

c*********************************************************************72
c
cc HSR5
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      integer la

      real a(la,*)
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      real s

      if ( m .lt. n ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: argument m must be .ge. n in routine hsr5'
        stop
      end if

      if ( m .eq. n ) then
        call hsr3 ( a, la, m, n )
        return
      end if

      k = m

20    continue

      s = -a(k,n)
      do i = n, m
        a(i,k) = s*a(i,n)
      end do
      a(k,k) = 1.0E+00 + a(k,k)
      l = n
40    j = l
      l = l - 1
      if ( l .eq. 0 ) go to 70
      s = 0.0E+00
      do i = j, m
        s = s + a(i,l)*a(i,k)
      end do
      a(l,k) = -s*a(l,l)
      do i = j, m
        a(i,k) = a(i,k) - s*a(i,l)
      end do
      go to 40
70    k = k - 1
      if ( k .gt. n ) go to 20
      call hsr3 ( a, la, m, n )

      return
      end
      subroutine htr ( a, b, c, z, bn, cn, zn, d, uu, n, h )

c*********************************************************************72
c
cc HTR
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real a(*)
      real b(*)
      real bn(*)
      real c(*)
      real cn(*)
      real d(*)
      real e
      integer f
      integer g
      integer h
      integer hg
      integer hw
      integer h1
      integer h2
      integer i
      integer ix
      integer iy
      integer j
      integer jl
      integer jp
      integer k
      integer l
      integer m
      integer n
      integer nh
      integer o
      real p
      real q
      real qq
      real r
      real s
      real t
      integer u
      real uu(*)
      integer v
      integer w
      real x
      real y
      integer z(*)
      integer zn(*)

      e = 65536.0**(-3)

      do i = 1, n
        d(i) = 1.0E+00
      end do

      g = h + 1
      f = g + 1
      hg = h*g
      h1 = h - 1
      h2 = h1 + h1
      o = h2 - n
      nh = n - h
      w = 1

20    continue

      hw = h + w
      v = w
      w = w + 1
      if ( w .ge. n ) go to 410

      t = 1.0E+00
      do i = w, n
        if ( d(i) .lt. t ) t = d(i)
      end do

      if ( t .gt. e ) go to 130

      do 120 i = w, n
        t = sqrt ( d(i) )
        d(i) = 1.0E+00
        if ( i .lt. g ) go to 70
        j = g*i - hg/2
        k = j - h
40      a(j) = t*a(j)
        if ( j .eq. k ) go to 50
        j = j - 1
        go to 40
50      k = j + f * min ( h, n - i )
60      a(j) = t*a(j)
        j = j + f
        if ( j .le. k ) go to 60
        go to 120
70      j = (i*(i+1))/2
        k = j - i + 1
80      a(j) = t*a(j)
        if ( j .eq. k ) go to 90
        j = j - 1
        go to 80
90      k = i
100     a(j) = t*a(j)
        k = k + 1
        j = j + k
        if ( k .le. g ) go to 100
        k = g - min ( i, n - h )
        if ( k .ge. h ) go to 120
110     a(j) = t*a(j)
        j = j + f
        k = k + 1
        if ( k .lt. h ) go to 110
120   continue

130   u = hg/2 + f * min ( v, nh ) - v - f
      j = u
      k = u + f
      ix = min ( hw, n )
      if ( ix .eq. g ) j = j + 1
      m = max ( 0, v - nh )
      go to 330
140   if ( ix .gt. nh ) go to 20
      ix = iy + h2
      u = u + hg
      j = u
      m = max ( 0, ix - n )
      if ( m .eq. 0 ) go to 150
      ix = n
      m = m - 1
      j = j - m*f
      k = j - f
150   m = m + 1
      s = a(j)
      jp = j + 1
      do 170 i = 1, m
        t = a(jp)
        if ( z(i) .gt. 0 ) go to 160
        a(j) = t + s*b(i)
        s = t*c(i) - s
        j = jp
        jp = jp + 1
        go to 170
160     a(j) = s - t*b(i)
        s = t + s*c(i)
        j = jp
        jp = jp + 1
170   continue
      a(j) = s
      if ( ix .lt. n ) go to 190
      if ( m-o .gt. iy ) go to 190
      if ( m .eq. h1 ) go to 20
      j = k
      go to 150
190   k = j + g
      if ( z(m) .gt. 0 ) go to 200
      t = -a(k)
      a(k) = -t*b(m)
      go to 210
200   t = a(k)*c(m)
210   l = j - h
220   iy = ix
      ix = ix - 1
      x = d(ix)
      y = d(iy)
      if ( abs ( s ) .ge. abs ( t ) ) go to 230
      r = x/y
      p = -s/t
      q = r*p*p
      if ( q .le. 1.0E+00 ) go to 250
      p = t/s
      qq = q/(1.0E+00+q)
      q = p/r
      go to 290
230   if ( s .eq. 0.0E+00 ) go to 240
      r = y/x
      p = t/s
      q = r*p*p
      if ( q .le. 1.0E+00 ) go to 280
      p = -s/t
      qq = q/(1.0E+00+q)
      q = p/r
      go to 260
240   p = 0.0E+00
      q = 0.0E+00
      qq = 1.0E+00
      go to 290
250   qq = 1.0E+00 /(1.0E+00+q)
      q = r*p
260   zn(m) = 0
      a(j) = q*s - t
      d(ix) = qq*y
      d(iy) = qq*x
      k = k - j + l
      jl = j - 1
      j = l
      l = k + 1
      jp = j + 1
      s = a(l)
      a(l) = a(j) + p*s
      a(j) = q*a(j) - s
      t = q*s - a(k)
      r = s + p*a(k)
      a(j) = q*a(j) - t
      s = a(l)
      a(k) = s + p*r
      a(l) = q*s - r
      if ( jp .gt. jl ) go to 310
      do i = jp, jl
        l = l + 1
        s = a(l)
        a(l) = a(i) + p*s
        a(i) = q*a(i) - s
      end do
      go to 310
280   qq = 1.0E+00 /(1.0E+00+q)
      q = r*p
290   zn(m) = 1
      a(j) = s + q*t
      d(ix) = qq*x
      d(iy) = qq*y
      k = k - j + l
      jl = j - 1
      j = l
      l = k + 1
      jp = j + 1
      s = a(l)
      a(l) = s - p*a(j)
      a(j) = a(j) + q*s
      t = s + q*a(k)
      r = a(k) - p*s
      a(j) = a(j) + q*t
      s = a(l)
      a(k) = r - p*s
      a(l) = s + q*r

      do i = jp, jl
        l = l + 1
        s = a(l)
        a(l) = s - p*a(i)
        a(i) = a(i) + q*s
      end do

310   continue

      bn(m) = p
      cn(m) = q
      if ( m .eq. h1 ) go to 350
      if ( ix .lt. hw ) go to 320
      j = j - 2 - m
      go to 150
320   if ( ix .gt. g ) go to 340
      j = j - v
      k = j + ix
330   l = j + w - ix
      m = m + 1
      s = a(j)
      t = a(k)
      k = k - 1
      go to 220
340   k = j + ix - v
      j = k - f
      go to 330
350   do i = 1, h1
        b(i) = bn(i)
        c(i) = cn(i)
        z(i) = zn(i)
      end do
      j = j + f
      l = max ( iy - nh, 2 )
      k = g - h
      m = h
370   j = j + k + m
      if ( ix .gt. m ) go to 380
      j = iy + h - m
      j = 2 + (j*j+j)/2
380   m = m - 1
      if ( m .lt. l ) go to 140
      s = a(j)
      jp = j + 1

      do i = m, h1
        t = a(jp)
        if ( z(i) .le. 0 ) then
          a(j) = t + s*b(i)
          s = t*c(i) - s
          j = jp
          jp = jp + 1
        else
          a(j) = s - t*b(i)
          s = t + s*c(i)
          j = jp
          jp = jp + 1
        end if
      end do

      a(j) = s
      go to 370
410   t = 1.0E+00
      d(1) = a(1)
      j = 1
      do i = 2, n
        k = i - 1
        if ( i .gt. g ) go to 420
        j = j + k
        go to 430
420     j = j + g
430     s = d(i)
        d(i) = s*a(j)
        s = sqrt ( s )
        uu(k) = s*t*a(j+1)
        t = s
      end do

      return
      end
      subroutine htr2 ( a, b, c, z, bn, cn, zn, d, uu, n, h, v, lv )

c*********************************************************************72
c
cc HTR2
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input, integer LV, the leading dimension of the array V.
c
      implicit none

      integer lv

      real a(*)
      real b(*)
      real bn(*)
      real c(*)
      real cn(*)
      real d(*)
      real e
      integer f
      integer g
      integer h
      integer h1
      integer h2
      integer hg
      integer hw
      integer i
      integer ix
      integer iy
      integer j
      integer k
      integer l
      integer m
      integer n
      integer nh
      integer o
      real p
      real q
      real qq
      real r
      real s
      real t
      integer u
      real uu(*)
      real v(lv,*)
      integer w
      integer w0
      real x
      real y
      integer z(*)
      integer zn(*)

      e = 65536.0**(-3)

      do i = 1, n
        d(i) = 1.0E+00
      end do

      g = h + 1
      f = g + 1
      hg = h*g
      h1 = h - 1
      h2 = h1 + h1
      o = h2 - n
      nh = n - h
      w = 1

20    continue

      hw = h + w
      w0 = w
      w = w + 1
      if ( w .ge. n ) go to 450

      t = 1.0E+00
      do i = w, n
        t = min ( t, d(i) )
      end do

      if ( t .gt. e ) go to 140

      do 130 j = w, n
        t = sqrt ( d(j) )
        d(j) = 1.0E+00
        do i = 1, n
          v(i,j) = t*v(i,j)
        end do
        if ( j .lt. g ) go to 80
        k = g*j - hg/2
        l = k - h
50      a(k) = t*a(k)
        if ( k .eq. l ) go to 60
        k = k - 1
        go to 50
60      l = k + f * min ( h, n - j )
70      a(k) = t*a(k)
        k = k + f
        if ( k .le. l ) go to 70
        go to 130
80      k = (j*(j+1))/2
        l = k - j + 1
90      a(k) = t*a(k)
        if ( k .eq. l ) go to 100
        k = k - 1
        go to 90
100     l = j
110     a(k) = t*a(k)
        l = l + 1
        k = k + l
        if ( l .le. g ) go to 110
        l = g - min ( j, n - h )
        if ( l .ge. h ) go to 130
120     a(k) = t*a(k)
        k = k + f
        l = l + 1
        if ( l .lt. h ) go to 120
130   continue

140   u = hg/2 + f * min ( w0, nh ) - w0 - f
      j = u
      k = u + f
      ix = min ( hw, n )
      if ( ix .eq. g ) j = j + 1
      m = max ( 0, w0 - nh )
      go to 360
150   if ( ix .gt. nh ) go to 20
      ix = iy + h2
      u = u + hg
      j = u
      m = max ( 0, ix - n )
      if ( m .eq. 0 ) go to 160
      ix = n
      m = m - 1
      j = j - m*f
      k = j - f
160   m = m + 1
      i = 0
      s = a(j)
170   if ( i .eq. m ) go to 190
      i = i + 1
      t = a(j+1)
      if ( z(i) .gt. 0 ) go to 180
      a(j) = t + s*b(i)
      s = t*c(i) - s
      j = j + 1
      go to 170
180   a(j) = s - t*b(i)
      s = t + s*c(i)
      j = j + 1
      go to 170
190   a(j) = s
      if ( ix .lt. n ) go to 200
      if ( m-o .gt. iy ) go to 200
      if ( m .eq. h1 ) go to 20
      j = k
      go to 160
200   k = j + g
      if ( z(i) .gt. 0 ) go to 210
      t = -a(k)
      a(k) = -t*b(i)
      go to 220
210   t = a(k)*c(i)
220   l = j - h
230   iy = ix
      ix = ix - 1
      x = d(ix)
      y = d(iy)
      if ( abs ( s ) .ge. abs ( t ) ) go to 240
      r = x/y
      p = -s/t
      q = r*p*p
      if ( q .le. 1.0E+00 ) go to 260
      p = t/s
      qq = q/(1.0E+00+q)
      q = p/r
      go to 310
240   if ( s .eq. 0. ) go to 250
      r = y/x
      p = t/s
      q = r*p*p
      if ( q .le. 1.0E+00 ) go to 300
      p = -s/t
      qq = q/(1.0E+00+q)
      q = p/r
      go to 270
250   p = 0.0E+00
      q = 0.0E+00
      qq = 1.0E+00
      go to 310
260   qq = 1.0E+00 /(1.0E+00+q)
      q = r*p
270   zn(m) = 0
      a(j) = q*s - t
      d(ix) = qq*y
      d(iy) = qq*x

      do i = 1, n
        t = v(i,ix)
        s = v(i,iy)
        v(i,ix) = q*t - s
        v(i,iy) = t + p*s
      end do

290   j = j - 1
      s = a(k)
      a(k) = a(j) + p*s
      a(j) = q*a(j) - s
      k = k - 1
      if ( j .gt. l ) go to 290
      t = q*s - a(k)
      r = s + p*a(k)
      a(j) = q*a(j) - t
      i = k + 1
      s = a(i)
      a(k) = s + p*r
      a(i) = q*s - r
      go to 340
300   qq = 1.0E+00 /(1.0E+00+q)
      q = r*p
310   zn(m) = 1
      a(j) = s + q*t
      d(ix) = qq*x
      d(iy) = qq*y

      do i = 1, n
        t = v(i,ix)
        s = v(i,iy)
        v(i,ix) = t + q*s
        v(i,iy) = s - p*t
      end do

330   j = j - 1
      s = a(k)
      a(k) = s - p*a(j)
      a(j) = a(j) + q*s
      k = k - 1
      if ( j .gt. l ) go to 330
      t = s + q*a(k)
      r = a(k) - p*s
      a(j) = a(j) + q*t
      i = k + 1
      s = a(i)
      a(k) = r - p*s
      a(i) = s + q*r
340   bn(m) = p
      cn(m) = q
      if ( m .eq. h1 ) go to 380
      if ( ix .lt. hw ) go to 350
      j = j - 2 - m
      go to 160
350   if ( ix .gt. g ) go to 370
      j = j - w0
      k = j + ix
360   l = j + w - ix
      m = m + 1
      s = a(j)
      t = a(k)
      k = k - 1
      go to 230
370   k = j + ix - w0
      j = k - f
      go to 360
380   i = 1
390   b(i) = bn(i)
      c(i) = cn(i)
      z(i) = zn(i)
      i = i + 1
      if ( i .lt. h ) go to 390
      j = j + f
      l = max ( iy - nh, 2 )
      k = g - h
      m = h
400   j = j + k + m
      if ( ix .gt. m ) go to 410
      j = iy + h - m
      j = 2 + (j*j+j)/2
410   m = m - 1
      if ( m .lt. l ) go to 150
      i = m
      s = a(j)
420   t = a(j+1)
      if ( z(i) .gt. 0 ) go to 430
      a(j) = t + s*b(i)
      s = t*c(i) - s
      go to 440
430   a(j) = s - t*b(i)
      s = t + s*c(i)
440   j = j + 1
      i = i + 1
      if ( i .lt. h ) go to 420
      a(j) = s
      go to 400
450   t = 1.0E+00
      d(1) = a(1)
      l = 1
      do j = 2, n
        k = j - 1
        if ( j .gt. g ) go to 460
        l = l + k
        go to 470
460     l = l + g
470     s = d(j)
        d(j) = s*a(l)
        s = sqrt ( s )
        uu(k) = s*t*a(l+1)
        t = s
        do i = 1, n
          v(i,j) = v(i,j)*t
        end do
      end do

      return
      end
      subroutine hvals ( e, d, u, a, la, n, h, w )

c*********************************************************************72
c
cc HVALS computes the eigenvalues of a symmetric band matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real a(*)
      real d(*)
      real e(*)
      integer h
      integer la
      integer n
      real u(*)
      real w(*)

      call hhess ( d, u, a, la, n, h, w )

      call tvals ( e, u, d, u, n, w )

      return
      end
      subroutine hvect ( e, x, a, la, n, h, w )

c*********************************************************************72
c
cc HVECT computes the eigenvectors of a symmetric band matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
      implicit none

      integer la
      integer n

      real a(la,n)
      integer b
      integer c
      integer d
      real e
      integer f
      integer g
      integer h
      integer i
      integer j
      integer k
      integer l
      integer m
      integer o
      integer p
      integer q
      real r
      real s
      real t
      real u
      real v
      real w(*)
      real x(n)
      integer y
      integer z

      if ( h .le. 0 ) then

        s = abs ( a(1,1) - e )

        do j = 1, n
          x(j) = 0.0E+00
          t = abs ( a(1,j) - e )
          if ( t .le. s ) then
            s = t
            k = j
          end if
        end do

        x(k) = 1.0E+00
        e = a(1,k)
        return

      end if

      c = 1 + 3*h
      d = h + h
      p = h + 1
      m = 1 - h

      do 70 j = 1, n

        x(j) = 0.0E+00
        l = m + d
        m = m + p
        do i = m, l
          w(i) = 0.0E+00
        end do
        m = l + p
        w(m) = a(1,j) - e
        f = m
        l = n - j
        if ( l .ge. h ) go to 40
        if ( l .eq. 0 ) go to 70
        l = l + m
        go to 50
40      l = m + h
50      m = m + 1
        k = 2 - m

        do i = m, l
          t = a(k+i,j)
          w(i) = t
          f = f + c
          w(f) = t
        end do

70    continue

      r = 0.0E+00
      q = p + p
      l = c + 1
      do i = q, l
        t = abs ( w(i) )
        if ( t .gt. r ) r = t
      end do

      i = -h
      z = p + h
      k = 0
90    k = k + 1
      i = i + l
      if ( k .eq. n ) go to 160
      m = i + 1
      q = i
      o = min ( h, n - k )
      p = i + o
      s = abs ( w(i) )

      do j = m, p
        t = abs ( w(j) )
        if ( s .le. t ) then
          q = j
          s = t
        end if
      end do

      if ( r .lt. s ) go to 110
      r = s
      y = k
      if ( r .eq. 0 ) go to 170
110   j = i - z
      b = q - i
      w(j) = k + b
      t = w(q)
      w(q) = w(i)
      w(i) = t
      do j = m, p
        w(j) = w(j)/t
      end do
      f = i + c * min ( d, n - k )
      g = c - o
130   m = p + g
      p = m + b
      t = w(p)
      w(p) = w(m)
      w(m) = t
      p = m + o
      if ( t .eq. 0.0E+00 ) go to 150
      q = i - m
      m = m + 1
      do j = m, p
        w(j) = w(j) - t*w(j+q)
      end do
150   if ( p .lt. f ) go to 130
      go to 90
160   if ( abs ( w(i) ) .ge. r ) go to 170
      r = w(i)
      y = n
170   k = y
      j = c*k - h
      t = 1.0E+00
      go to 190
180   t = x(k)/w(j+k)
190   x(k) = t
      if ( k .eq. 1 ) go to 210
      q = max ( 1, k - d )
      k = k - 1
      do i = q, k
        x(i) = x(i) - t*w(i+j)
      end do
      j = j - c
      go to 180
210   if ( r .eq. 0.0E+00 ) go to 360
      s = 0.0E+00
      do i = 1, n
        t = abs ( x(i) )
        if ( t .gt. s ) s = t
      end do

      r = 0.0E+00
      s = 1.0E+00 /s
      do i = 1, n
        t = s*x(i)
        r = r + t*t
        x(i) = t
      end do

      j = -h
      q = 1
      p = h + 1
      do i = 1, p
        w(i+1) = x(i)
      end do
250   j = j + c
      k = q
      if ( k .eq. n ) go to 270
      q = k + 1
      i = w(j+k-z)
      t = x(i)
      x(i) = x(k)
      x(k) = t
      o = min ( k + h, n )
      do i = q, o
        x(i) = x(i) - t*w(i+j)
      end do
      if ( o .eq. n ) go to 250
      w(q+j) = x(o+1)
      go to 250
270   u = 0.0E+00
280   t = x(k)/w(j+k)
      if ( abs ( t ) .gt. u ) u = abs ( t )
      x(k) = t
      if ( k .eq. 1 ) go to 300
      q = max ( 1, k - d )
      k = k - 1
      do i = q, k
        x(i) = x(i) - t*w(i+j)
      end do
      j = j - c
      go to 280
300   v = 0.0E+00
      do i = 1, p
        v = v + x(i)*w(i+1)
      end do
      if ( p .eq. n ) go to 330
      j = z + 1
      k = p + 1
      do i = k, n
        v = v + x(i)*w(j)
        j = j + l
      end do
330   if ( v .ne. 0.0E+00 ) v = r/v
      s = 0.0E+00
      t = 1.0E+00 /u
      u = 0.0E+00
      do i = 1, n
        s = s + (x(i)*v)**2
        u = u + (t*x(i))**2
      end do
      t = t / sqrt ( u )
      do i = 1, n
        x(i) = t*x(i)
      end do
      if ( r+r .ge. s ) e = e + v
      return

360   u = 0.0E+00
      do i = 1, n
        t = abs ( x(i) )
        if ( t .gt. u ) u = t
      end do

      t = 1.0E+00 /u

      u = 0.0E+00
      do i = 1, n
        u = u + (t*x(i))**2
      end do

      t = t / sqrt ( u )

      do i = 1, n
        x(i) = t*x(i)
      end do

      return
      end
      subroutine hvert ( v, a )

c*********************************************************************72
c
cc HVERT inverts a symmetric band matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, real V((N*(N+1))/2), the inverse matrix.
c
c    Input, real A(4+(H+1)*N), the factorization information from HFACT.
c
      implicit none

      real a(*)
      integer e
      integer f
      integer g
      integer h
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      integer o
      integer p
      real t
      real v(*)

      t = a(1)

      if ( abs ( t ) .ne. 1232 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor with hfact before inverting'
        stop
      end if

      if ( t .le. 0.0E+00 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: matrix has no inverse'
        stop
      end if

      n = a(2)
      h = a(4)
      e = h + 1
      f = 4 - h - h
      p = 0

      do 90 g = 1, n
        o = p + 1
        p = o + n - g
        j = f - o + g*e
        k = o
        do i = o, p
          v(i) = 0.0E+00
        end do
        v(o) = 1.0E+00
        if ( h .gt. 0 ) go to 40
        v(o) = v(o)/a(4+g)
        go to 90
40      j = j + h
        if ( k .eq. p ) go to 60
        t = v(k)/a(j+k)
        l = min ( p, k + h )
        k = k + 1
        if ( t .eq. 0.0E+00 ) go to 40
        do i = k, l
          v(i) = v(i) - t*a(i+j)
        end do
        go to 40
60      v(p) = v(p)/a(j+k)
70      if ( k .eq. o ) go to 90
        j = j - h
        m = k
        k = k - 1
        l = min ( p, k + h )
        t = v(k)
        do i = m, l
          t = t - v(i)*a(i+j)
        end do
        v(k) = t/a(j+k)
        go to 70
90    continue

      return
      end
      function icon ( a, b )

c*********************************************************************72
c
cc ICON estimates the condition number of a symmetric matrix.
c
c  Discussion:
c
c    Symmetric pivoting is used.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c    William Hager,
c    Condition Estimates,
c    SIAM Journal on Scientific and Statistical Computing,
c    Volume 5, Number 2, June 1984, pages 311-316.
c
c  Parameters:
c
c    Input, real A(*), the factorization information from IFACT.
c
c    Workspace, real B(N).
c
c    Output, real ICON, the condition number of the matrix.
c
      implicit none

      real a(*)
      real b(*)
      real c
      real d
      integer i
      real icon
      integer j
      integer m
      integer n

      c = a(1)

      if ( abs ( c ) .ne. 1237 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor with ifact before',
     &  'estimating condition'
        stop
      end if

      if ( c .lt. 0.0E+00 ) then
        icon = -1.0E+00
        return
      end if

      m = 0
      n = a(2)
      c = 1.0E+00 / a(2)
      do j = 1, n
        b(j) = c
      end do

      go to 50

30    do j = 1, n
        b(j) = 0.0E+00
      end do
      b(m) = 1.0E+00

50    call isolve ( b, a, b )

      c = 0.0E+00
      do j = 1, n
        c = c + abs ( b(j) )
        if ( b(j) .lt. 0.0E+00 ) then
          b(j) = -1.0E+00
        else
          b(j) = 1.0E+00
        end if
      end do

      call isolve ( b, a, b )
c
c  I is the index of the largest entry of B.
c
      i = 1
      do j = 1, n
        if ( abs ( b(i) ) .lt. abs ( b(j) ) ) then
          i = j
        end if
      end do

      if ( m .eq. 0 ) go to 90
      if ( m .eq. i ) go to 100
      if ( d .ge. c ) go to 100
90    m = i
      d = c
      go to 30

100   c = c * a(3)
      c = max ( c, 1.0E+00 )
      icon = c
      return
      end
      function idet ( iexp, a )

c*********************************************************************72
c
cc IDET computes the determinant of a symmetric matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, integer IEXP, a power of ten that is part of the
c    determinant.
c
c    Input, real A(*), factorization information from IFACT.
c
c    Output, real IDET, the mantissa of the determinant.
c    Determinant = IDET * 10.0^IEXP.
c
      implicit none

      real a(*)
      real c
      real d
      real f
      real g
      integer h
      integer i
      real idet
      integer iexp
      integer k
      integer l
      integer n

      iexp=0
      idet=0.0E+00
      d = a(1)

      if ( abs ( d ) .ne. 1237 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor with ifact before',
     &  ' computing determinant'
        return
      end if

      iexp = 0
      if ( d .lt. 0.0E+00 ) go to 80
      n = a(2)
      if ( n .eq. 1 ) go to 90
      d = 1.0E+00
      f = 2.0E+00**64
      g = 1.0E+00 /f
      h = 64
      k = 10 + (n*(n-3))/2
      l = k - 5 + 4*n

      do 50 i = k, l, 4
        if ( a(i-2) .ne. 0.0E+00 ) d = -d
        d = d*a(i)
20      if ( abs ( d ) .gt. f ) go to 40
30      if ( abs ( d ) .gt. g ) go to 50
        iexp = iexp - h
        d = d*f
        go to 30
40      iexp = iexp + h
        d = d*g
        go to 20
50    continue

      d = d*a(l+1)
      if ( iexp .ne. 0 ) go to 60
      idet = d
      return
60    if ( d .eq. 0.0E+00 ) go to 100
      c = alog10 ( abs ( d ) ) + iexp * alog10 ( 2.0E+00 )
      iexp = c
      c = c - iexp
      if ( c .le. 0.0E+00 ) go to 70
      c = c - 1
      iexp = iexp + 1
70    f = 10.0E+00**c
      if ( d .lt. 0.0E+00 ) f = -f
      idet = f
      return
80    idet = 0.0E+00
      return
90    idet = a(9)
      return
100   iexp = 0
      go to 80
      end
      subroutine ifact ( a, n )

c*********************************************************************72
c
cc IFACT computes the LU factorization of a symmetric matrix.
c
c  Discussion:
c
c    Partial pivoting is used.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input/output, real A(5+(N*(N+7))/2); on input, information defining the
c    rows of the matrix, starting from the diagonal.  On output, 
c    factorization information.
c
c    Input, integer N, the number of rows and columns in the matrix.
c
      implicit none

      real a(*)
      integer b
      integer c
      integer d
      integer e
      integer f
      integer g
      integer h
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      integer o
      integer p
      integer q
      real r
      real s
      real t
      real u
      real v
      real w

      m = (n+n*n)/2
      l = m + n
      i = m
10    i = i + 1
      a(i) = 0.0E+00
      if ( i .lt. l ) go to 10
      i = -l
      j = m
      k = m
      r = 0.0E+00
      s = 0.0E+00
20    i = i + l - k
      k = k + 1
      j = k
      s = abs ( a(i+j) )
30    if ( j .eq. l ) go to 40
      j = j + 1
      t = abs ( a(i+j) )
      s = s + t
      a(j) = a(j) + t
      go to 30
40    s = s + a(k)
      if ( r .lt. s ) r = s
      if ( k .lt. l ) go to 20
      j = m + 3
50    a(j) = a(j-3)
      j = j - 1
      if ( j .gt. 3 ) go to 50
      a(1) = 1237
      a(2) = n
      a(3) = r
      if ( n .gt. 1 ) go to 60
      a(9) = a(4)
      a(4) = 1235
      a(5) = 1
      a(6) = abs ( a(9) )
      if ( a(9) .ne. 0.0E+00 ) return
      a(1) = -1237
      a(4) = -1235
      return
60    if ( n .eq. 2 ) go to 250
      e = 7 + (n*(n+5))/2
      h = n
      k = 4
70    g = h - 1
      d = n - g
      if ( h .gt. 2 ) go to 80
      c = 0
      go to 150
80    l = k + g
      i = k + 1
      p = i
      do j = i, l
        if ( abs ( a(j) ) .gt. abs ( a(p) ) ) p = j
      end do
      s = a(p)
      a(e+d) = d + p - k
      c = p - i
      if ( s .eq. 0.0E+00 ) go to 150
      if ( c .eq. 0 ) go to 130
      a(p) = a(i)
      a(i) = s
      i = k + h + 1
      l = i + c - 2
      p = l + g
      if ( i .gt. l ) go to 110
      o = g + i - 2

      do j = i, l
        t = a(j)
        a(j) = a(p)
        a(p) = t
        p = p + o - j
      end do

110   j = k + h
      t = a(j)
      a(j) = a(p)
      a(p) = t
      i = l + 2
      l = k + g + g
      if ( i .gt. l ) go to 130
      o = (c*(g+g-c-1))/2
      do j = i, l
        t = a(j)
        p = j + o
        a(j) = a(p)
        a(p) = t
      end do
130   i = k + 2
      l = k + g
      do j = i, l
        a(j) = a(j)/s
      end do
150   q = k + g + g
      p = k + h + 1
      if ( d .gt. 1 ) go to 160
      w = a(k+h)
      o = -g
      go to 210
160   m = n - 1
      v = a(k+h)
      t = 0.0E+00
      s = a(5)
      i = 4
      j = 4 + d + c
      u = a(j)
      l = n
      if ( h .eq. 2 ) go to 230
170   i = i + l
      r = a(j)
      b = j - c
      o = b - p + 1
      a(j) = a(b)
      a(b) = r
      j = j + m
      r = a(i+1)
      if ( i .eq. k ) go to 190
      w = s*t + u*a(i) + a(j)*r
      do f = p, q
        a(f) = a(f) - w*a(f+o)
      end do
      v = v - u*w
      s = r
      t = u
      u = a(j)
      l = m
      m = m - 1
      go to 170
190   w = s*t + u*a(i) + r

      do f = p, q
        a(f) = a(f) - w*a(f+o)
      end do

      w = v - u*w
      a(k+h) = w - u*r
      o = o + h
210   do f = p, q
        a(f) = a(f) - w*a(f+o)
      end do
      k = k + h
      h = h - 1
      go to 70
230   i = i + l
      j = j + m
      r = a(i+1)
      if ( i .eq. k ) go to 240
      w = s*t + u*a(i) + a(j)*r
      v = v - u*w
      s = r
      t = u
      u = a(j)
      l = m
      m = m - 1
      go to 230
240   w = s*t + u*a(i) + r
      v = v - u*w
      a(k+h) = v - u*r
250   i = 4
      k = 4
      h = n
      m = 5 + (n*(n+1))/2
260   a(m) = a(k)
      a(m+1) = a(k+1)
      if ( h .eq. 2 ) go to 280
      o = k - i + 2
      l = i - 3 + h
      do j = i, l
        a(j) = a(j+o)
      end do
      i = l + 1
      k = k + h
      h = h - 1
      m = m + 2
      go to 260
280   m = m + 2
      a(m) = a(k+2)
      a(m+1) = 0.0E+00
      i = 6 + (n*(n+1))/2
      l = i + n + n - 2
      k = i - n - n - 1
      m = k

      do j = i, l, 2
        a(k) = a(j)
        a(k+1) = a(j-1)
        a(k+2) = a(j)
        k = k + 3
      end do

      call pfact ( a(m), 3, n )

      if ( a(m) .lt. 0.0E+00 ) then
        a(1) = -1237
      end if

      return
      end
      subroutine inp ( a, x, y, z, u, v, w, l, r )

c*********************************************************************72
c
cc INP
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real a
      real b
      real c
      real l
      real p
      real q
      real r
      real s
      real t
      real u
      real v
      real w
      real x
      real y
      real z

      s = z - x
      t = (y-x)/s
      a = (u-v)/t + (w-v)/(1.0E+00-t)
      if ( a .eq. 0.0E+00 ) go to 40
      b = .5*(w-u)/a - .5
      c = u/a
      t = sqrt ( abs ( c ) )
      if ( abs ( b ) .lt. sign(t,c) ) go to 60
      t = max ( t, abs ( b ) )
      if ( t .eq. 0.0E+00 ) go to 50
      q = 1.0E+00 /t
      p = sqrt ( ( q * b )**2 - q*c*q)
      p = t*p
      if ( abs ( p + b ) .gt. abs ( p - b ) ) go to 10
      q = p - b
      go to 20
10    q = -(b+p)
20    p = c/q
      q = x + s*q
      p = x + s*p
      if ( q .lt. l ) go to 30
      if ( q .gt. r ) go to 30
      a = q
      return
30    a = p
      return
40    if ( u .eq. w ) go to 50
      a = x + s*u/(u-w)
      return
50    a = l
      return
60    a = x - s*b
      return
      end
      subroutine ins ( s, f, a, b, c, fa, fb, fc, j, y, z )

c*********************************************************************72
c
cc INS
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real a
      real b
      real c
      real f
      real fa
      real fb
      real fc
      integer j
      real s
      real y(*)
      real z(*)

      j = j + 1
      y(j) = s
      z(j) = f

      if ( f .le. fa ) then
        c = b
        b = a
        a = s
        fc = fb
        fb = fa
        fa = f
      else if ( f .le. fb ) then
        c = b
        b = s
        fc = fb
        fb = f
      else if ( f .le. fc ) then
        c = s
        fc = f
      end if

      return
      end
      function isig ( x )

c*********************************************************************72
c
cc ISIG returns a sign function of a real number.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input, real X, a number.
c
c    Output, integer ISIG,
c    -1, 0, or +1, depending on whether X is negative, zero, or positive.
c
      implicit none

      integer isig
      real x

      if ( x .lt. 0.0E+00 ) then
        isig = -1
      else if ( x .eq. 0.0E+00 ) then
        isig = 0
      else
        isig = +1
      end if

      return
      end
      subroutine isolve ( x, a, b )

c*********************************************************************72
c
cc ISOLVE solves a symmetric linear system.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, real X(N), the solution.
c
c    Input, real A(*), factorization information from IFACT.
c
c    Input, real B(N), the right hand side.
c
      implicit none

      real a(*)
      real b(*)
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      integer o
      real t
      real x(*)

      t = a(1)

      if ( abs ( t ) .ne. 1237 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor with ifact before solving'
        stop
      end if

      n = a(2)

      if ( n .le. 2 ) then
        call psolve ( x, a(4), b )
        return
      end if

      o = 6 + (n*(n+5))/2
      m = n - 1
      if ( t .lt. 0.0E+00 ) go to 80

      do i = 1, n
        x(i) = b(i)
      end do

      do k = 2, m
        j = a(k+o)
        t = x(j)
        x(j) = x(k)
        x(k) = t
      end do

      k = 2
      i = 1

50    continue

      t = x(k)
      k = k + 1

      if ( t .ne. 0.0E+00 ) then
        do j = k, n
          x(j) = x(j) - t*a(i+j)
        end do
      end if

      i = i + n - k
      if ( k .lt. n ) go to 50

80    continue

      i = 5 + (n*(n-3))/2
      call psolve ( x, a(i), x )
      k = n
      j = m
      l = i - 1 - n

90    continue

      t = x(j)

      do i = k, n
        t = t - x(i)*a(i+l)
      end do

      k = j
      x(k) = t
      j = j - 1
      l = l - n + k
      if ( k .gt. 2 ) go to 90

      k = m

110   continue

      j = a(k+o)
      t = x(j)
      x(j) = x(k)
      x(k) = t
      k = k - 1
      if ( k .gt. 1 ) go to 110

      return
      end
      subroutine ivert ( v, w )

c*********************************************************************72
c
cc IVERT inverts a symmetric matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input/output, real V(3+(N*(N+1))/2).  On input, the factorization
c    information from IFACT.  On output, the inverse matrix.
c
c    Workspace, real W(2*N-2).
c
      implicit none

      integer h
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      integer o
      integer p
      real t
      real v(*)
      real w(*)

      t = v(1)

      if ( abs ( t ) .ne. 1237 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor with ifact before inverting'
        stop
      end if

      if ( t .le. 0.0E+00 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: matrix has no inverse'
        stop
      end if

      n = v(2)

      if ( n .le. 1 ) then
        v(1) = 1.0E+00 /v(9)
        return
      end if

      if ( n .le. 2 ) then
        w(1) = 1.0E+00
        w(2) = 0.0E+00
        call psolve ( w, v(4), w )
        v(1) = w(1)
        v(2) = w(2)
        w(1) = 0.0E+00
        w(2) = 1.0E+00
        call psolve ( w, v(4), w )
        v(3) = w(2)
        return
      end if

      o = 1 + (n*(n+1))/2
      l = n + n - 2
      k = 9 + o + l
      m = n + 1
      j = k - m
      do i = m, l
        w(i) = v(i+j)
      end do
      j = n + n - 4

60    continue

      k = k - 1
      v(k+j) = v(k)
      if ( k .gt. 1 ) go to 60

      m = n - 1
      do i = 1, n
        w(i) = 0.0E+00
      end do
      w(1) = 1.0E+00
      call psolve ( w, v(o), w )
      l = n
      j = m
      p = o - 1 - n

80    continue

      t = w(j)
      do i = l, n
        t = t - w(i)*v(i+p)
      end do
      l = j
      w(l) = t
      j = j - 1
      p = p - n + l
      if ( l .gt. 2 ) go to 80

      do i = 1, n
        v(i) = w(i)
      end do

      do k = 2, n

        do i = 1, n
          w(i) = 0.0E+00
        end do
        w(k) = 1.0E+00
        if ( k .eq. n ) go to 140
        i = (k*(n+n-k-1))/2
        l = k
120     t = w(l)
        l = l + 1
        do j = l, n
          w(j) = w(j) - t*v(i+j)
        end do
        i = i + n - l
        if ( l .lt. n ) go to 120
140     call psolve ( w, v(o), w )
        l = n
        j = m
        p = o - 1 - n
150     t = w(j)
        do i = l, n
          t = t - w(i)*v(i+p)
        end do
        l = j
        w(l) = t
        j = j - 1
        p = p - n + l
        if ( l .gt. k ) go to 150
        j = ((k-1)*(n+n-k))/2
        do i = k, n
          v(i+j) = w(i)
        end do

      end do

      k = n
190   k = k - 1
      if ( k .le. 1 ) return
      l = w(m+k)
      if ( k .eq. l ) go to 190
      h = l
      o = n
      i = k
      j = 1 + ((k-1)*(n+n-k+2))/2
200   t = v(i)
      v(i) = v(l)
      v(l) = t
      o = o - 1
      i = i + o
      l = l + o
      if ( i .lt. j ) go to 200
      p = j
      o = k - j - n - 1
      i = l
210   j = j + 1
      i = i - j - o
      if ( j .eq. l ) go to 220
      t = v(j)
      v(j) = v(i)
      v(i) = t
      go to 210
220   t = v(p)
      v(p) = v(i)
      v(i) = t
      if ( h .eq. n ) go to 190
      j = i + n - h
230   i = i + 1
      l = l + 1
      t = v(i)
      v(i) = v(l)
      v(l) = t
      if ( i .lt. j ) go to 230
      go to 190
      end
      subroutine jt ( p, q, n, v, s, s1 )

c*********************************************************************72
c
cc JT (???Jenkins-Traub?))
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real amag
      integer i
      integer m
      integer n
      complex p(*)
      complex q(*)
      complex r
      complex s
      real s1
      complex t
      complex u
      complex v
      complex w

      if ( s1 .gt. 1.0E+00 ) go to 70
      u = (0.0E+00,0.0E+00)
      v = (1.0E+00,0.0E+00)
      i = n
10    v = p(i) + v*s
      u = q(i) + u*s
      i = i - 1
      if ( i .gt. 0 ) go to 10
      if ( amag ( u ) .ne. 0.0 ) go to 50
      m = n
20    r = q(m)
      q(m) = (0.0,0.0)
      i = m
30    i = i - 1
      u = q(i) + r*s
      q(i) = r
      r = u
      if ( i .gt. 1 ) go to 30
      u = (0.0,0.0)
      i = m
40    u = q(i) + u*s
      i = i - 1
      if ( i .gt. 0 ) go to 40
      m = m - 1
      if ( amag ( u ) .eq. 0.0 ) go to 20
50    r = v/u
      t = p(n) - r*q(n) + s
      q(n) = (1.0E+00,0.0)
      i = n
60    i = i - 1
      u = p(i) - r*q(i) + t*s
      q(i) = t
      t = u
      if ( i .gt. 1 ) go to 60
      return
70    w = 1.0E+00 /s
      u = q(1)
      v = p(1)
      do i = 2, n
        u = u*w + q(i)
        v = v*w + p(i)
      end do
      u = u*w
      v = v*w + 1.0E+00
      if ( amag ( u ) .ne. 0. ) go to 120
      m = n
90    r = q(m)
      q(m) = (0.0,0.0)
      i = m
100   i = i - 1
      u = q(i) + r*s
      q(i) = r
      r = u
      if ( i .gt. 1 ) go to 100
      u = (0.0,0.0)
      m = m - 1
      do i = 1, n
        u = u*w + q(i)
      end do
      u = u*w
      if ( amag ( u ) .eq. 0. ) go to 90
120   r = v/u

      u = p(1) - r*q(1)
      q(1) = -w*u
      do i = 2, n
        u = p(i) - r*q(i) + u*w
        q(i) = -w*u
      end do
      q(n) = (1.0E+00,0.0)

      return
      end
      function kcon ( a, b )

c*********************************************************************72
c
cc KCON estimates the condition number of a general matrix.
c
c  Modified:
c
c    20 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c    William Hager,
c    Condition Estimates,
c    SIAM Journal on Scientific and Statistical Computing,
c    Volume 5, Number 2, June 1984, pages 311-316.
c
c  Parameters:
c
c    Input, real A(*), the factorization information from KFACT.
c
c    Workspace, real B(N).
c
c    Output, real KCON, the condition number of the matrix.
c
      implicit none

      real a(*)
      real b(*)
      real c
      real d
      integer i
      integer j
      real kcon
      integer m
      integer n

      c = a(1)

      if ( abs ( c ) .ne. 1236 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor with KFACT',
     &  'before estimating condition'
        stop
      end if

      if ( c .lt. 0.0E+00 ) then
        kcon = -1.0E+00
        return
      end if

      m = 0
      n = a(2)
      c = 1.0E+00 /a(2)
      do j = 1, n
        b(j) = c
      end do
      go to 50
30    do j = 1, n
        b(j) = 0.0E+00
      end do
      b(m) = 1.0E+00

50    call ksolve ( b, a, b )

      c = 0.0E+00
      do j = 1, n
        c = c + abs ( b(j) )
        if ( b(j) .lt. 0.0E+00 ) then
          b(j) = -1.0E+00
        else
          b(j) = 1.0E+00
        end if
      end do

      call ktrans ( b, a, b )
c
c  I is the index of the largest entry of B.
c
      i = 1
      do j = 1, n
        if ( abs ( b(i) ) .lt. abs ( b(j) ) ) then
          i = j
        end if
      end do

      if ( m .eq. 0 ) go to 90
      if ( m .eq. i ) go to 100
      if ( d .ge. c ) go to 100
90    m = i
      d = c
      go to 30
100   c = c * a(3)
      c = max ( c, 1.0E+00 )
      kcon = c
      return
      end
      function kdet ( iexp, a )

c*********************************************************************72
c
cc KDET computes the determinant of a general matrix.
c
c  Discussion:
c
c    The determinant is KDET * 10^IEXP.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, integer IEXP, the power of 10 of the determinant.
c
c    Input, real A(N*N+2*N+2), factorization information from KFACT.
c
c    Output, real KDET, the mantissa of the determinant.
c    Determinant = KDET * 10.0^IEXP.
c
      implicit none

      real a(*)
      real c
      real d
      real f
      real g
      integer h
      integer i
      integer iexp
      integer j
      integer k
      real kdet
      integer l
      integer m
      integer n

      iexp = 0
      kdet = 0.0
      d = a(1)

      if ( abs ( d ) .ne. 1236 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'kdet   error, matrix must be factored by KFACT',
     &  ' before kdet is called.'
        return
      end if

      iexp = 0

      if ( d .lt. 0.0 ) then
        kdet = 0.0
        return
      end if

      n = a(2)

      if ( n .eq. 1 ) go to 80
      d = 1.0E+00
      f = 2.0E+00**64
      g = 1.0E+00 /f
      h = 64
      m = n + 1
      j = 0
      k = 4
      l = 3 - m + m*n
      n = l + m

      do i = k, l, m
        j = j + 1
        if ( a(i) .gt. j ) d = -d
        if ( a(j+n) .gt. j ) d = -d
        d = d*a(i+j)
20      if ( abs ( d ) .lt. f ) go to 30
        iexp = iexp + h
        d = d*g
        go to 20
30      if ( abs ( d ) .gt. g ) go to 40
        iexp = iexp - h
        d = d*f
        go to 30
40      continue
      end do

      d = d*a(l+m)
      if ( iexp .ne. 0 ) go to 50
      kdet = d
      return

50    if ( d .eq. 0. ) go to 90
      c = alog10 ( abs ( d ) ) + iexp * alog10 ( 2.0E+00 )
      iexp = c
      c = c - iexp
      if ( c .le. 0.0 ) go to 60
      c = c - 1
      iexp = iexp + 1
60    f = 10.0E+00**c
      if ( d .lt. 0. ) f = -f
      kdet = f
      return
80    kdet = a(5)
      return
90    iexp = 0
      kdet = 0.0
      return
      end
      subroutine kfact ( a, la, n )

c*********************************************************************72
c
cc KFACT factors a general matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input/output, real A(N*N+2*N+2).
c    On input, the initial NxN segment contains the matrix.
c    On output, factorization information.
c    A(1) contains 1236 if the factorization was successful, -1236 otherwise.
c    A(2) contains N.
c    A(3) contains the L1 norm of the matrix.
c
c    Input, integer LA, the leading dimension of A.
c
c    Input, integer N, the number of rows and columns in the matrix.
c
      implicit none

      integer la
      integer n

      real a(la*n)
      integer b
      integer c
      integer d
      integer e
      integer f
      integer g
      integer h
      integer i
      integer j
      integer k
      integer l
      integer m
      integer o
      integer p
      integer q
      real r
      real s
      real t

      if ( n .lt. la ) then
        call pack ( a, la, n )
      end if

      r = 0.0E+00
      o = n + 1
      p = o + 1
      l = 5 + n * p
      i = - n - 3

10    continue

      l = l - o

      if ( l .ne. 4 ) then

        s = 0.0E+00
        do k = 1, n
          j = l - k
          t = a(i+j)
          a(j) = t
          s = s + abs ( t )
        end do

        r = max ( r, s )
        i = i + 1
        go to 10

      end if

      a(1) = 1236
      a(2) = n
      a(3) = r
      q = 3 + n * o
      i = 5 - p
      b = 0
      k = 1
40    i = i + p
      if ( k .eq. n ) go to 120
      e = n - k
      h = i

      do m = i, q, o
        l = i + e
        do j = m, l
          if ( abs ( a(j) ) .gt. abs ( a(h) ) ) h = j
        end do
      end do

      c = (h-4)/o
      d = 4 + o*c + k
      g = h - d
      h = d - i
      l = i + e
      f = i - b
      do j = f, l
        t = a(j)
        m = j + h
        a(j) = a(m)
        a(m) = t
      end do
      j = i - k
      a(j) = g + k
      h = g + i
      t = a(h)
      a(h) = a(i)
      a(i) = t
      b = k
      k = k + 1
      if ( t .eq. 0. ) go to 120
      m = i + 1
      do j = m, l
        a(j) = a(j)/t
      end do
      f = i + e*o
80    j = k + l
      h = j + g
      t = a(h)
      a(h) = a(j)
      a(j) = t
      l = e + j
      if ( t .eq. 0. ) go to 100
      h = i - j
      m = j + 1

      do j = m, l
        a(j) = a(j) - t*a(j+h)
      end do

100   if ( l .lt. f ) go to 80
      a(l+b) = c + 1
      go to 40
110   a(1) = -1236
      return
120   if ( a(i) .eq. 0 ) go to 110
      return
      end
      subroutine ksolve ( x, a, b )

c*********************************************************************72
c
cc KSOLVE solves a general linear system.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, real X(N), the solution.
c
c    Input, real A(*), factorization information from KFACT.
c
c    Input, real B(N), the right hand side.
c
      implicit none

      real a(*)
      real b(*)
      integer h
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      real t
      real x(*)

      t = a(1)

      if ( abs ( t ) .ne. 1236 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor with KFACT before solving'
        stop
      end if

      n = a(2)
      h = n
      m = n + 1
      j = 4 - m
      if ( t .lt. 0. ) go to 100
      do i = 1, n
        x(i) = b(i)
      end do
      k = 1
30    j = j + m
      if ( k .eq. n ) go to 50
      l = a(j)
      t = x(l)
      x(l) = x(k)
      x(k) = t
      k = k + 1
      if ( t .eq. 0. ) go to 30
      do i = k, n
        x(i) = x(i) - t*a(i+j)
      end do
      go to 30
50    t = x(k)/a(j+k)
60    x(k) = t
      if ( k .eq. 1 ) go to 80
      k = k - 1
      do i = 1, k
        x(i) = x(i) - t*a(i+j)
      end do
      j = j - m
      go to 50
80    l = 3 + m*n
90    h = h - 1
      if ( h .eq. 0 ) return
      i = a(h+l)
      t = x(h)
      x(h) = x(i)
      x(i) = t
      go to 90
100   k = 0
110   k = k + 1
      j = j + m
      if ( a(j+k) .ne. 0. ) go to 110
      do i = 1, n
        x(i) = 0.
      end do
      t = 1.0E+00
      h = k
      go to 60
      end
      subroutine ktrans ( x, a, b )

c*********************************************************************72
c
cc KTRANS solves a transposed general linear system.
c
c  Modified:
c
c    20 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, real X(N), the solution.
c
c    Input, real A(*), factorization information from KFACT.
c   
c    Input, real B(N), the right hand side.
c   
      implicit none

      real a(*)
      real b(*)
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      real t
      real x(*)

      t = a(1)

      if ( abs ( t ) .ne. 1236 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor with kfact before solving'
        stop
      end if

      n = a(2)
      m = n + 1
      if ( t .lt. 0. ) go to 110

      do i = 1, n
        x(i) = b(i)
      end do

      j = 3 + m*n
      k = n - 1
      do i = 1, k
        l = a(i+j)
        t = x(i)
        x(i) = x(l)
        x(l) = t
      end do

      t = 0.
      j = 4
      k = 1
50    if ( x(k) .ne. 0. ) go to 60
      k = k + 1
      j = j + m
      if ( k .le. n ) go to 50
      return
60    x(k) = (x(k)-t)/a(j+k)
70    if ( k .eq. n ) go to 90
      t = 0.
      j = j + m
      do i = 1, k
        t = t + a(i+j)*x(i)
      end do
      k = k + 1
      go to 60
90    if ( k .eq. 1 ) return
      j = j - m
      t = x(k-1)

      do i = k, n
        t = t - x(i)*a(i+j)
      end do

      k = k - 1
      i = a(j)
      x(k) = x(i)
      x(i) = t
      go to 90
110   i = 5 + n + m*n
      l = m
120   i = i - m - 1
      l = l - 1
      if ( a(i) .ne. 0. ) go to 120
      j = j + m*(l-k)
      k = l
      do i = 1, n
        x(i) = 0.
      end do
      x(k) = 1.0E+00
      go to 70
      end
      subroutine kvert ( v, lv, n, w )

c*********************************************************************72
c
cc KVERT inverts a general matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input/output, real V(LV,N); on input, the NxN matrix.  On output,
c    the NxN inverse matrix.
c
c    Input, integer LV, the leading dimension of the array V.
c
c    Input, integer N, the order of the matrix.
c
c    Workspace, real W(2*N).
c
      implicit none

      integer lv
      integer n

      integer h
      integer i
      integer j
      integer k
      integer l
      integer m
      integer o
      integer p
      integer q
      real s
      real t
      real v(lv,n)
      real w(2*n)

      if ( n .eq. 1 ) then
        if ( v(1,1) .eq. 0.0E+00 ) then
          write(*,*) 'kvert - matrix is not invertible!'
          stop
        end if
        v(1,1) = 1.0E+00 / v(1,1)
        return
      end if

      o = n + 1
      l = 0
      m = 1

10    if ( l .eq. n ) go to 90

      k = l
      l = m
      m = m + 1
      p = l
      q = l
      s = abs ( v(l,l) )

      do h = l, n
        do i = l, n
          t = abs ( v(i,h) )
          if ( t .gt. s ) then
            p = i
            q = h
            s = t
          end if
        end do
      end do

      w(n+l) = p
      w(o-l) = q

      do i = 1, n
        t = v(i,l)
        v(i,l) = v(i,q)
        v(i,q) = t
      end do

      s = v(p,l)
      v(p,l) = v(l,l)

      if ( s .eq. 0. ) then
        write(*,*) 'kvert - matrix is not invertible!'
        stop
      end if

      v(l,l) = -1.0E+00
      s = 1.0E+00 /s
      do i = 1, n
        v(i,l) = -s*v(i,l)
      end do

      j = l

50    j = j + 1

      if ( j .gt. n ) j = 1

      if ( j .eq. l ) go to 10

      t = v(p,j)
      v(p,j) = v(l,j)
      v(l,j) = t

      do i = 1, k
        v(i,j) = v(i,j) + t*v(i,l)
      end do
      v(l,j) = s*t
      do i = m, n
        v(i,j) = v(i,j) + t*v(i,l)
      end do

      go to 50

90    l = w(k+n)

      do i = 1, n
        t = v(i,l)
        v(i,l) = v(i,k)
        v(i,k) = t
      end do

      k = k - 1
      if ( k .gt. 0 ) go to 90

      do j = 1, n
        do i = 2, n
          p = w(i)
          h = o - i
          t = v(p,j)
          v(p,j) = v(h,j)
          v(h,j) = t
        end do
      end do

      return
      end
      subroutine lancz ( q, d, u, p, j, n, mult, w )

c*********************************************************************72
c
cc LANCZ performs an iteration of Lanczos method to reduce a matrix to tridiagonal form.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      integer n

      real d(*)
      integer i
      integer j
      external mult
      real p(*)
      real q(n)
      real u(*)
      real w(n,*)
      real s
      real t

      if ( j .gt. n ) then
        return
      end if

      if ( j .eq. 1 ) then

        t = 0.
        do i = 1, n
          t = t + q(i)**2
        end do

        if ( t .eq. 0. ) then
          return
        end if

        t = 1.0E+00 / sqrt ( t )
        do i = 1, n
          q(i) = t*q(i)
          w(i,1) = q(i)
        end do

        call mult ( w, q )

        t = 0.
        do i = 1, n
          t = t + q(i) * w(i,1)
        end do

        d(1) = t
        do i = 1, n
          w(i,2) = w(i,1) - t*q(i)
          w(i,1) = q(i)
        end do

      else

        t = 0.
        do i = 1, n
          t = t + w(i,2)**2
        end do

        if ( t .eq. 0. ) then
          return
        end if

        p(j-1) = t
        t = sqrt ( t )
        u(j-1) = t

        s = 1.0E+00 /t
        do i = 1, n
          w(i,2) = s*w(i,2)
          q(i) = w(i,2)
        end do

        call mult ( w(1,2), q )

        s = 0.
        do i = 1, n
          w(i,1) = w(i,2) - t*w(i,1)
          s = s + q(i)*w(i,1)
        end do

        d(j) = s

        do i = 1, n
          w(i,2) = w(i,1) - s*q(i)
          w(i,1) = q(i)
        end do

      end if

      return
      end
      subroutine mgrid ( k )

c*********************************************************************72
c                                                    
cc MGRID solves a 1D PDE using the multigrid method.
c
c  Discussion:
c
c    Solves  -x''(t) = 1, 
c    with the boundary condtions x(0) = x(1) = 0  by the multigrid method.  
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input, integer K, the mesh index.  K should be 1 or greater.
c    The number of intervals used is 2^K.
c
      implicit none

      real d0
      real d1
      real difmax
      integer i
      integer it
      integer j
      integer k
      integer l
      integer ll
      integer m
      integer n
      integer nl
      real r(518)
      real s
      real tol
      real u(518)
      real utol

      write ( *, * ) ' '
      write ( *, * ) 'MGRID'
      write ( *, * ) '  example of multigrid method.'
      write ( *, * ) '  solve -u''(t)=1, u(0)=u(1)=0'
      write ( *, * ) '  solution is u(t)=0.5*(-t*t+t)'
      write ( *, * ) ' '

      write ( *, '(a,i4)' ) '  Mesh index K = ', k
      n = 2**k
      write ( *, '(a,i6)' ) '  Number of intervals N=2^K = ', n
      it = 4
      tol = 0.0001
      utol = 0.7
c 
c  Set the right hand side.
c 
      do i = 1, n
        r(i) = 1.0E+00
      end do
c 
c  Initialize.
c 
      s = (1.0E+00/n)**2
      do i = 1, n
        r(i) = s * r(i)
      end do
      m = n
      l = 1
      ll = n
      nl = n + n + k - 2
      do i = 1, nl
        u(i) = 0.0E+00
      end do
      d1 = 0.0E+00
c 
c  Gauss-seidel iteration
c 
10    continue

      j = 0

50    continue

      d0 = d1
      d1 = 0.0E+00
      j = j + 1
      i = l

60    continue

      i = i + 1
      s = 0.5E+00 * ( u(i-1) + u(i+1) + r(i) )
      d1 = d1 + abs ( s - u(i) )
      u(i) = s
      if ( i .lt. ll ) then
        go to 60
      end if

      write(*,70) d1
70    format(' dif:',f20.10)

      if ( j .lt. it ) goto 50

      if ( d1  .lt.  tol ) then
        write ( *, * ) '  Time to go up!'
        go to 100
      end if

      if ( d1 / d0 .lt. utol ) then
        write ( *, * ) '  D1/D0 = ', d1/d0
        goto 50
      end if

      if ( n .eq. 2 ) then
        write ( *, *) 'HEY WHAT THE HELL?'
        go to 50
      end if
c
c  coarser mesh (slow convergence)
c
      i = ll + 2

80    continue

      l = l + 2
      i = i + 1

      if ( l .le. ll ) then
        u(i) = 0.0E+00
        r(i) = 4.0E+00 * ( r(l) + u(l-1) - 2.0E+00 * u(l) + u(l+1) )
        go to 80
      end if

      n = n / 2
      write(*,*) '  Go down to mesh intervals:', n
      ll = ll + n + 1
      l = l + 1
      go to 10
c     
c  finer mesh (fast convergence)
c     
100   if ( n .eq. m ) go to 120
      i = l - 3
      j = ll

110   continue

      u(i) = u(i) + u(j)
      u(i+1) = u(i+1) + 0.5E+00 * ( u(j) + u(j+1) )
      i = i - 2
      j = j - 1

      if ( j .gt. l ) then
        go to 110
      end if

      n = n + n

      write(*,*) '  Go up to mesh intervals:',n

      ll = l - 2
      l = i
      u(i+1) = u(i+1) + 0.5E+00 * u(j+1)

      go to 10
c 
c  Computation completed.
c 
120   continue

      do i = 1, n + 1
        s = real ( i - 1 ) / real ( n )
        write(*,'(i5,f10.5,f15.6)') i, s, u(i)
      end do

      difmax = 0.0E+00
      do i = 1, n + 1
        s = real ( i - 1 ) / real ( n )
        difmax = max ( difmax, abs ( u(i) - 0.5 * ( - s * s + s ) ) )
      end do
      write ( *, * ) 'maximum error=',difmax

      return
      end
      subroutine mpower ( ex, x, ey, y, n, ne, dif, size, ndigit, 
     &  limit, mult, w )

c*********************************************************************72
c
cc MPOWER
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real a
      real b
      real c
      real d
      real dif
      real e
      complex ex
      complex ey
      real f
      real g
      integer h
      integer i
      integer j
      integer k
      integer l
      integer limit
      integer m
      integer m1
      integer m2
      external mult
      integer n
      integer n2
      integer n3
      integer n4
      integer n5
      integer ndigit
      integer ne
      integer np
      integer ns
      integer o
      integer p
      integer p1
      real q
      real r
      real s
      real size
      real t
      real w(*)
      real x(*)
      real y(*)

      save np

      data np /0/

      if ( ne .eq. 0 ) then
        write(*,*) 'since the requested eigenpairs have been computed,'
        write(*,*) 'routine mpower will quit'
        return
      end if

      if ( ne .ne. np ) then
        h = 0
        ns = ne
        n3 = (ns*(3+n+n-ns))/2
      end if

      o = h + 1
      p = o + 1
      p1 = p + 1
      m = n + 1
      m1 = n + h
      m2 = n + o
      n2 = n + n
      n4 = (o*p)/2 + n3 - 2
      n5 = (h*(n+n-o))/2
      q = 10.0E+00**(-2*ndigit)
      l = 0
      j = 1
      do i = 1, n
        x(i) = x(j)
        j = j + 2
      end do

      if ( h .ne. 0 ) then

        call xp ( x, w(m), h, n )

        do i = 1, h
          x(i) = 0.0E+00
        end do

      end if

      t = 0.0E+00
      do i = o, n
        t = t + x(i)**2
      end do

      if ( t .eq. 0. ) then
        write(*,*) 'starting guess for routine mpower cannot be zero'
        write(*,*) 'and it cannot lie in the space spanned by the'
        write(*,*) 'previously computed eigenvectors'
        stop
      end if

      t = 1.0E+00 / sqrt ( t )
      do i = o, n
        x(i) = t*x(i)
      end do

90    continue

      do i = 1, n
        y(i) = x(i)
      end do

      call stp ( mult, y, w, h, m, n, n5 )

      do i = 1, n
        y(i+n) = y(i)
      end do

      call stp ( mult, y(m), w, h, m, n, n5 )
      a = 0.

      do i = o, n
        a = a + y(i+n)**2
      end do

      if ( a .eq. 0. ) go to 670

      a = 1.0E+00 / sqrt ( a )
      s = 0.
      do i = o, n
        s = s + x(i)*y(i)
        j = i + n
        y(j) = y(j)*a
      end do
      a = x(o)
      b = 1.0E+00 + abs ( a )

      if ( a .lt. 0. ) go to 140
      a = a + 1.0E+00
      s = s + y(o)
      go to 150
140   a = a - 1.0E+00
      s = s - y(o)
150   s = s/b

      y(o) = y(o) - a*s
      c = 0.
      d = 0.
      do i = p, n
        t = y(i) - s*x(i)
        y(i) = t
        c = c + t*x(i)
        d = d + t*t
      end do

      if ( d .eq. 0. ) then
        c = x(p)
        d = 1.0E+00
        y(p) = 1.0E+00
      end if

      d = 1.0E+00 / sqrt ( d )
      c = c/b
      y(o) = -c*a*d
      do i = p, n
        y(i) = d*(y(i)-c*x(i))
      end do

      s = 0.
      t = 0.
      do i = o, n
        j = i + n
        s = s + y(j)*x(i)
        t = t + y(j)*y(i)
      end do

      d = 0.
      do i = o, n
        a = y(i+n)
        d = d + (a-s*x(i)-t*y(i))**2
        x(i) = a
      end do

      l = l + 2
      if ( l .ge. limit ) go to 210
      if ( d .gt. q ) go to 90
210   s = 0.
      do i = o, n
        s = s + x(i)*y(i)
      end do
      a = x(o)
      b = 1.0E+00 + abs ( a )
      if ( a .lt. 0. ) go to 230
      a = a + 1.0E+00
      s = s + y(o)
      go to 240
230   a = a - 1.0E+00
      s = s - y(o)
240   s = s/b

      y(o) = y(o) - a*s
      c = 0.
      d = 0.
      do i = p, n
        t = y(i) - s*x(i)
        y(i) = t
        c = c + t*x(i)
        d = d + t*t
      end do

      if ( d .ne. 0. ) go to 260
      c = x(p)
      d = 1.0E+00
      y(p) = 1.0E+00
260   c = c/b
      d = 1.0E+00 / sqrt ( d )
      y(o) = -c*a*d
      y(m2) = y(o)
      x(m2) = x(o)
      do i = p, n
        j = i + n
        t = d*(y(i)-c*x(i))
        y(i) = t
        y(j) = t
        x(j) = x(i)
      end do
      call stp ( mult, x, w, h, m, n, n5 )
      call stp ( mult, y, w, h, m, n, n5 )

      do i = 1, h
        x(i) = 0.
        y(i) = 0.
      end do

      a = 0.
      b = 0.
      c = 0.
      d = 0.
      s = 0.
      t = 0.
      do i = o, n
        j = i + n
        f = x(i)
        g = y(i)
        s = s + f*f
        t = t + g*g
        a = a + f*x(j)
        b = b + g*x(j)
        c = c + f*y(j)
        d = d + g*y(j)
      end do

      if ( t .eq. 0. ) go to 710
      if ( s .eq. 0. ) go to 760

      f = 1.0E+00 / sqrt ( s )
      g = 1.0E+00 / sqrt ( t )
      r = 0.
      s = 0.
      t = 0.
      do i = o, n
        j = i + n
        r = r + (x(i)-a*x(j))**2
        s = s + (x(i)-a*x(j)-c*y(j))**2
        t = t + (y(i)-b*x(j)-d*y(j))**2
        x(i) = f*x(i)
        y(i) = g*y(i)
      end do

      t = f*f*s + g*g*t
      s = f*f*r
      l = l + 2
      if ( l .ge. limit ) go to 400
      if ( t .le. q ) go to 410
      if ( s .gt. q ) go to 210
320   size = 3*l + 1
330   dif = sqrt ( s )
340   ex = cmplx(a,0.)
      ey = (0.,0.)
      do i = 1, n
        x(i+n) = x(i)
      end do
      call stp ( mult, x(m), w, h, m, n, n5 )
      call rvc ( a, x, x(m), w(n3), h )
      call px ( x, w(m), h, n, n5 )
      if ( o .ge. ns ) go to 380
      call hse ( x(m), w(m), h, o, n )
      do i = 1, n
        y(i) = 0.
      end do
      y(o) = 1.0E+00
      n5 = (o*(n+n-p))/2
      call stp ( mult, y, w, o, m, n, n5 )
      y(p) = 0.
      do i = 1, p
        w(n4+i) = y(i)
      end do
380   j = n2
      i = n
390   k = j - 1
      y(k) = 0.
      y(j) = 0.
      x(k) = x(i)
      x(j) = 0.
      j = j - 2
      i = i - 1
      if ( i .gt. 0 ) go to 390
      ne = ne - 1
      np = ne
      h = o
      return
400   size = 3*l + 2
      if ( s+s .lt. t ) go to 330
      go to 420
410   size = 3*l + 1
420   dif = sqrt ( t )
      if ( h .eq. 0 ) go to 440
      do i = m, m1
        x(i) = 0.
        y(i) = 0.
      end do

440   if ( p .ge. ns ) go to 500

      call hse ( x, w(m), h, o, n )
      t = 0.
      do i = o, n
        t = t + x(i)*y(i)
      end do
      do i = o, n
        y(i) = y(i) - t*x(i)
      end do
      call hse ( y, w(m), o, p, n )
      do i = 1, n
        x(i) = 0.
        y(i) = 0.
      end do
      x(o) = 1.0E+00
      n5 = (p*(n+n-p1))/2
      call stp ( mult, x, w, p, m, n, n5 )
      do i = 1, p
        w(i+n4) = x(i)
      end do
      y(p) = 1.0E+00
      call stp ( mult, y, w, p, m, n, n5 )
      n5 = (h*(n+n-o))/2
      y(p1) = 0.
      n4 = n4 + p
      do i = 1, p1
        w(i+n4) = y(i)
      end do

500   s = .5*(d-a)
      t = b*c
      r = s*s + t
      if ( r .lt. 0. ) go to 610
      r = abs ( s ) + sqrt ( r )
      if ( s .lt. 0. ) go to 510
      ex = cmplx(a-t/r,0.)
      ey = cmplx(a+r,0.)
      go to 520
510   ex = cmplx(a+t/r,0.)
      ey = cmplx(a-r,0.)
520   s = cabs ( a - ex ) + abs ( b )
      t = cabs ( d - ex ) + abs ( c )
      if ( s .gt. t ) go to 540
      if ( t .ne. 0. ) go to 530
      e = 0.
      f = 1.0E+00
      a = 1.0E+00
      b = 0.
      go to 580
530   e = c/t
      f = (ex-d)/t
      go to 550
540   e = (ex-a)/s
      f = b/s
550   s = cabs ( a - ey ) + abs ( b )
      t = cabs ( d - ey ) + abs ( c )
      if ( s .gt. t ) go to 570
      if ( t .ne. 0. ) go to 560
      e = 0.
      f = 1.0E+00
      a = 1.0E+00
      b = 0.
      go to 580
560   a = c/t
      b = (ey-d)/t
      go to 580
570   a = (ey-a)/s
      b = b/s

580   do i = o, n
        j = i + n
        s = f*x(j) + e*y(j)
        t = b*x(j) + a*y(j)
        x(i) = s
        y(i) = t
        x(j) = s
        y(j) = t
      end do

      call stp ( mult, x(m), w, h, m, n, n5 )
      t = ex
      call rvc ( t, x, x(m), w(n3), h )
      call px ( x, w(m), h, n, n5 )
      call stp ( mult, y(m), w, h, m, n, n5 )
      t = ey
      call rvc ( t, y, y(m), w(n3), h )
      call px ( y, w(m), h, n, n5 )
      j = n2
      i = n
600   k = j - 1
      x(k) = x(i)
      x(j) = 0.
      y(k) = y(i)
      y(j) = 0.
      i = i - 1
      j = j - 2
      if ( i .gt. 0 ) go to 600
      ne = ne - 2
      ne = max ( ne, 0 )
      np = ne
      h = p
      return
610   q = a + s
      r = sqrt ( - r )
      ex = cmplx(q,r)
      ey = cmplx(q,-r)
      s = cabs ( a - ex ) + abs ( b )
      t = cabs ( d - ex ) + abs ( c )
      if ( s .gt. t ) go to 650
      b = c/t
      c = (q-d)/t
      a = r/t

      do i = o, n
        j = i + n
        t = c*x(j) + b*y(j)
        s = a*x(j)
        x(i) = t
        y(i) = s
        x(j) = t
        y(j) = s
      end do

630   call stp(mult,x(m),w,h,m,n,n5)
      call stp(mult,y(m),w,h,m,n,n5)
      call cvc(q,r,x,y,x(m),y(m),w(n3),h)
      call px(x,w(m),h,n,n5)
      call px(y,w(m),h,n,n5)
      j = n2
      i = n
640   k = j - 1
      x(j) = y(i)
      y(j) = -y(i)
      x(k) = x(i)
      y(k) = x(i)
      j = j - 2
      i = i - 1
      if ( i .gt. 0 ) go to 640
      ne = ne - 2
      ne = max ( ne, 0 )
      np = ne
      h = p
      return
650   c = b/s
      b = (q-a)/s
      d = r/s
      do i = o, n
        j = i + n
        t = c*x(j) + b*y(j)
        s = d*y(j)
        x(i) = t
        y(i) = s
        x(j) = t
        y(j) = s
      end do

      go to 630
670   a = 0.
      dif = 0.
      size = 3*l + 1
      do i = o, n
        if ( y(i) .ne. 0. ) go to 690
      end do
      go to 340
690   do i = o, n
        x(i) = y(i)
      end do
      go to 340
710   if ( s .ne. 0. ) go to 730
      do i = o, n
        x(i) = x(i+n)
      end do
      go to 670
730   s = 1.0E+00 / sqrt ( s )
      a = a*s
      do i = o, n
        x(i) = x(i)*s
      end do

      s = 0.
      do i = o, n
        s = s + (x(i)-a*x(i+n))**2
      end do

      if ( s .gt. q ) go to 90

      a = 0.
      go to 320

760   t = 1.0E+00 / sqrt ( t )
      do i = o, n
        x(i) = t*y(i)
      end do

      go to 90
      end
      subroutine mult ( y, x, a, la, m, n )

c*********************************************************************72
c
cc MULT computes a matrix-vector product.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, real Y(M), the matrix-vector product A*X.
c
c    Input, real X(N), the vector to be multiplied.
c
c    Input, real A(LA,N), the MxN matrix.
c
c    Input, integer LA, the leading dimension of the matrix.
c
c    Input, integer M, N, the number of rows and columns of A.
c
      implicit none

      integer la
      integer m
      integer n

      real a(la,n)
      integer i
      integer j
      real t
      real x(n)
      real y(m)

      do i = 1, m
        y(i) = 0.0
      end do

      do j = 1, n
        do i = 1, m
          y(i) = y(i) + a(i,j) * x(j)
        end do
      end do

      return
      end
      subroutine newton ( t, s, k, d, u, n )

c*********************************************************************72
c
cc NEWTON applies one step of Newton's method to the characteristic polynomial for a tridiagonal matrix
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real a
      real b
      real c
      real d(*)
      real e
      real f
      real g
      real h
      integer i
      integer isig
      integer j
      integer k
      integer l
      integer n
      real p
      real q
      real s
      real t
      real u(*)
      real x
      real y
      real z

      k = 0

      if ( n .eq. 1 ) then
        s = d(1) - t
        t = d(1)
        if ( s .le. 0. ) then
          k = 1
        end if
        return
      end if

      z = 0.
      y = 65536.0E+00**4
      x = 1.0E+00 /y
      c = 2.
      b = c*(d(1)-t)
      e = -c
      f = z
      g = z
      h = z
      i = 1
      j = 1
      l = 1
      if ( isig(b) .eq. l ) go to 20
      k = k + 1
      l = -1
      if ( b .lt. z ) go to 20
      if ( u(i) .eq. z ) go to 90
20    p = u(i)
      i = i + 1
      q = d(i) - t
      a = q*b - p*c
      if ( isig(a) .eq. l ) go to 30
      k = k + 1
      l = -l
      if ( a .ne. 0. ) go to 30
      if ( i .ge. n ) go to 140
      if ( u(i) .eq. 0. ) go to 90
30    c = b
      b = a
      a = q*e - p*f - c
      f = e
      e = a
      a = q*g - p*h - f
      h = g
      g = a
      if ( i .ge. n ) go to 60
      a = abs ( b )
      if ( abs ( e ) .gt. a ) a = abs ( e )
      if ( abs ( g ) .gt. a ) a = abs ( g )
40    if ( a .gt. y ) go to 50
      if ( a .gt. x ) go to 20
      if ( a .eq. z ) go to 20
      a = a*y
      b = b*y
      c = c*y
      e = e*y
      f = f*y
      g = g*y
      h = h*y
      go to 40
50    a = a*x
      b = b*x
      c = c*x
      e = e*x
      f = f*x
      g = g*x
      h = h*x
      go to 40
60    if ( e .eq. z ) go to 150
      x = b/e
      if ( g .eq. z ) go to 70
      y = e/(g+g)
      s = x*y/(x-y)
      go to 80
70    s = -x
80    t = t + s
      return
90    i = i + 1
      b = d(i) - t
      c = 1.0E+00
      l = 1
      if ( isig(b) .eq. l ) go to 100
      k = k + 1
      l = -1
      if ( b .lt. 0. ) go to 100
      if ( i .ge. n ) go to 140
      if ( u(i) .eq. 0. ) go to 90
100   if ( i .ge. n ) go to 140
      j = i
      i = i + 1
      a = (d(i)-t)*b - u(j)*c
      if ( isig(a) .eq. l ) go to 110
      k = k + 1
      l = -l
      if ( a .ne. 0. ) go to 110
      if ( i .ge. n ) go to 140
      if ( u(i) .eq. 0. ) go to 90
110   c = b
      b = a
120   a = abs ( b )
      if ( a .gt. y ) go to 130
      if ( a .gt. x ) go to 100
      if ( a .eq. z ) go to 100
      c = c*y
      b = b*y
      go to 120
130   c = c*x
      b = b*x
      go to 120
140   s = 0.
      return
150   s = t + sign(1.0E+00,t)
      t = t + s
      return
      end
      function norm1 ( n, l, mult, tmult, w )

c*********************************************************************72
c
cc NORM1 estimates the 1-norm of a matrix.
c
c  Discussion:
c
c    The description above cannot be correct, since the 1-norm of a
c    matrix is simply the maximum of the sums of the absolute values
c    of the columns.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real c
      real d
      real e
      integer i
      integer j
      integer k
      integer l
      integer m
      external mult
      integer n
      real norm1
      integer o
      integer p
      integer q
      external tmult
      real w(*)

      k = 0
      m = 0
      o = 0
      c = 1.0E+00/n

      do j = 1, n
        w(j) = c
      end do

      if ( l .eq. 1 ) go to 50

      do j = 1, n
        w(j+n) = 0.
      end do

      go to 50

30    do j = 1, n
        w(j) = 0.
      end do
      w(m) = 1.0E+00

50    call mult ( w, m )

      c = 0.
      do 70 j = 1, n
        c = c + abs ( w(j) )
        if ( w(j) .lt. 0. ) go to 60
        w(j) = 1.0E+00
        go to 70
60      w(j) = -1.0E+00
70    continue

      call tmult ( w )
      i = 1
      do j = 1, n
        if ( abs ( w(i) ) .lt. abs ( w(j) ) ) i = j
      end do
      if ( m .eq. 0 ) go to 90
      if ( d .ge. c ) go to 110
      if ( m .eq. i ) go to 100
90    m = i
      d = c
      o = o + 1
      if ( l .gt. 1 ) w(i+n) = 1.0E+00
      go to 30
100   d = c
110   if ( o .ge. n ) go to 130
120   k = k + 1
      if ( k .lt. l ) go to 140
130   norm1 = d
      return
140   c = 1.0E+00 /(n-o)

      do j = 1, n
        w(j) = 0.
        if ( w(j+n) .le. 0. ) then
          w(j) = c
        end if
      end do

      m = 0
      q = o
      go to 180

160   do j = 1, n
        w(j) = 0.
      end do
      w(m) = 1.0E+00

180   call mult ( w, m )
      if ( m .gt. 0 ) go to 210
      do j = 1, n
        if ( w(j+n) .eq. 0. ) go to 200
      end do
200   m = j
      p = j

210   c = 0.
      do 230 j = 1, n
        c = c + abs ( w(j) )
        if ( w(j) .lt. 0. ) go to 220
        w(j) = 1.0E+00
        go to 230
220     w(j) = -1.0E+00
230   continue

      call tmult ( w )
      i = m
      do j = p, n
        if ( w(j+n) .le. 0. ) then
          if ( abs ( w(i) ) .lt. abs ( w(j) ) ) i = j
        end if
      end do
      if ( o .eq. q ) go to 250
      if ( e .ge. c ) go to 270
      if ( m .eq. i ) go to 260
250   m = i
      e = c
      o = o + 1
      w(i+n) = 1.0E+00
      if ( o .gt. n ) go to 270
      go to 160
260   e = c
270   if ( e .gt. d ) d = e
      go to 110
      end
      subroutine null ( b, lb, n, a, c )

c*********************************************************************72
c 
cc NULL computes an orthonormal basis for the space perpendicular to a given collection of vectors
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      integer lb

      real a(*)
      real b(lb,*)
      real c
      integer g
      integer h
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      integer o
      real s
      real t

      t = a(1)

      if ( abs ( t ) .ne. 3230 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor array of vectors using'
        write(*,*) 'routine qr before using routine null'
        write(*,*) 'to compute a basis for the null space'
        stop
      end if

      m = a(2)
      n = a(3)

      if ( lb .lt. m ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: the leading dimension of argument b in'
        write(*,*) 'routine null must be greater than or equal'
        write(*,*) 'to the number of components in each basis vector'
        stop
      end if

      k = 3
      l = min ( n, m )
      o = m + 1

      do j = 1, l
        if ( abs ( a(j+k) ) .le. c ) go to 110
        k = k + o
      end do

      if ( n .lt. m ) go to 40
      n = 0
      return
40    g = 3 + o*n - m
      k = m - n
50    l = k + n
      s = -a(l+g)
      do i = n, m
        b(i,k) = s*a(i+g)
      end do
      b(l,k) = b(l,k) + 1.0E+00
      h = g
      l = n
70    j = l
      l = l - 1
      if ( l .eq. 0 ) go to 100
      h = h - o
      s = 0.
      do i = j, m
        s = s + a(i+h)*b(i,k)
      end do
      b(l,k) = -s*a(l+h)
      do i = j, m
        b(i,k) = b(i,k) - s*a(i+h)
      end do
      go to 70
100   k = k - 1
      if ( k .gt. 0 ) go to 50
      n = m - n
      return
110   n = j - 1
      if ( n .gt. 0 ) go to 40
      n = m
      do j = 1, m
        do i = 1, m
          b(i,j) = 0.
        end do
        b(j,j) = 1.0E+00
      end do

      return
      end
      subroutine omult ( y, x, a, n, w )

c*********************************************************************72
c
cc OMULT multiplies a circulant matrix times a vector.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, complex Y(N), the product of the matrix times X.
c
c    Input, complex X(N), the vector to be multiplied.
c
c    Input, complex A(N), the matrix.
c
c    Input, integer N, the number of rows and columns in the matrix.
c
c    Workspace, complex W(N).
c
      implicit none

      integer n

      complex a(n)
      integer i
      real t
      complex w(n)
      complex x(n)
      complex y(n)

      t = n
      t = 1.0E+00 / t

      do i = 1, n
        y(i) = t * conjg ( x(i) )
      end do

      call fft ( y, n, w )

      do i = 1, n
        y(i) = conjg ( y(i) ) * a(i)
      end do

      call fft ( y, n, w )

      return
      end
      subroutine osolve ( x, a, n, b, w )

c*********************************************************************72
c
cc OSOLVE solves a circulant linear system.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, complex X(N), the solution.
c
c    Input, complex A(*), factorization information from OFACT.
c
c    Input, integer N, the number of rows and columns in the matrix.
c
c    Input, complex B(N), the right hand side.
c
c    Workspace, complex W(N).
c
      implicit none

      integer n

      complex a(*)
      complex b(n)
      integer i
      real t
      complex w(n)
      complex x(n)

      t = n
      t = 1.0E+00 / t

      do i = 1, n
        x(i) = t * conjg ( b(i) )
      end do

      call fft ( x, n, w )

      do i = 1, n
        x(i) = conjg ( x(i) ) / a(i)
      end do

      call fft ( x, n, w )

      return
      end
      subroutine ovals ( a, n, w )

c*********************************************************************72
c
cc OVALS computes the eigenvalues of a circulant matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      integer n

      complex a(n)
      complex w(n)

      call fft ( a, n, w )

      return
      end
      subroutine over ( x, a, b )

c*********************************************************************72
c
cc OVER computes the least squares solution to an overdetermined linear system.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real a(*)
      real b(*)
      integer h
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      integer o
      real t
      real x(*)

      t = a(1)

      if ( abs ( t ) .ne. 3230 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must qr factor coefficient matrix'
        write(*,*) 'before solving system'
        stop
      end if

      if ( t .le. 0. ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'singular system - compute regularized solution'
        write(*,*) '(see section 6-11)'
        stop
      end if

      m = a(2)
      n = a(3)

      if ( m .lt. n ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: the number of equations is less than'
        write(*,*) 'the number of unknowns. for an overdetermined'
        write(*,*) 'system, there must be more equations than unknowns'
        stop
      end if

      if ( m .eq. 1 ) then
        x(1) = b(1)/a(4)
        return
      end if

      do i = 1, m
        x(i) = b(i)
      end do

      o = m + 1
      l = n
      if ( m .eq. n ) l = n - 1

      k = 4
      do j = 1, l
        t = 0.
        do i = j, m
          t = t + x(i)*a(i+k)
        end do
        do i = j, m
          x(i) = x(i) - t*a(i+k)
        end do
        k = k + o
      end do

      j = n
      k = 3 + o*n
      h = k
90    k = k - o
      t = x(j)/a(j+k)
      x(j) = t
      if ( j .eq. 1 ) go to 110
      j = j - 1

      do i = 1, j
        x(i) = x(i) - t*a(i+k)
      end do

      go to 90

110   if ( n .eq. 1 ) return
120   j = a(l+h)
      t = x(j)
      x(j) = x(l)
      x(l) = t
      l = l - 1
      if ( l .gt. 0 ) go to 120
      return
      end
      subroutine overt ( a, n, w )

c*********************************************************************72
c
cc OVERT inverts a circulant matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input/output, complex A(N).  On input, information defining the 
c    NxN matrix.  On output, information defining the inverse matrix.
c
c    Input, integer N, the order of the matrix.
c
c    Workspace, W(N).
c
      implicit none

      integer n

      complex a(n)
      integer i
      real t
      complex w(n)

      call fft ( a, n, w )

      do i = 1, n
        t = abs ( real ( a(i) ) ) + abs ( aimag ( a(i) ) )
        if ( t .eq. 0. ) then
          write ( *, '(a)' ) ' '
          write(*,*) 'error: matrix has no inverse'
          stop
        end if
        a(i) = 1.0E+00 / conjg ( a(i) )
      end do

      call fft ( a, n, w )

      do i = 1, n
        a(i) = conjg ( a(i) )
      end do

      return
      end
      subroutine pack ( a, la, n )

c*********************************************************************72
c
cc PACK packs the elements of a square matrix into contiguous storage.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input/output, real A(*).
c    On input A is assumed to be an NxN matrix stored as an LAxN array.
c    On output, A has been packed as an NxN array.
c
c    Input, integer LA, the leading dimension of the array.
c    N <= LA.
c
c    Input, integer N, the number of rows and columns in the matrix.
c
      implicit none

      integer la
      integer n

      real a(la*n)
      integer h
      integer i
      integer j
      integer jhi
      integer jlo
      integer k
      integer l
      integer o

      h = la - n

      if ( h .eq. 0 ) then
        return
      end if

      if ( h .le. 0 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: la argument in pack must be .ge. n argument'
        stop
      end if

      i = 0
      k = 1
      l = n
      o = n * n

10    continue

      if ( l .eq. o ) then
        jlo = n * n + 1
        jhi = la * n
        do j = jlo, jhi
          a(j) = 0.0E+00
        end do
        return
      end if

      i = i + h
      k = k + n
      l = l + n
      do j = k, l
        a(j) = a(i+j)
      end do

      go to 10
      end
      function pcon ( a, b )

c*********************************************************************72
c
cc PCON estimates the condition number of a tridiagonal matrix.
c
c  Discussion:
c
c    Partial pivoting is used.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c    William Hager,
c    Condition Estimates,
c    SIAM Journal on Scientific and Statistical Computing,
c    Volume 5, Number 2, June 1984, pages 311-316.
c
c  Parameters:
c
c    Input, real A(*), the factorization information from PFACT.
c
c    Workspace, real B(N).
c
c    Output, real PCON, the condition number of the matrix.
c
      implicit none

      real a(*)
      real b(*)
      real c
      real d
      integer i
      integer j
      integer m
      integer n
      real pcon

      c = a(1)

      if ( abs ( c ) .ne. 1235 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor with pfact before',
     &  'estimating condition'
        stop
      end if

      if ( c .lt. 0.0E+00 ) then
        pcon = -1.0E+00
        return
      end if

      m = 0
      n = a(2)
      c = 1.0E+00 /a(2)
      do j = 1, n
        b(j) = c
      end do
      go to 50
30    do j = 1, n
        b(j) = 0.
      end do
      b(m) = 1.0E+00
50    call psolve ( b, a, b )

      c = 0.0E+00
      do j = 1, n
        c = c + abs ( b(j) )
        if ( b(j) .lt. 0.0E+00 ) then
          b(j) = -1.0E+00
        else
          b(j) = 1.0E+00
        end if
      end do

      call ptrans ( b, a, b )
c
c  I is the index of the largest entry of B.
c
      i = 1
      do j = 1, n
        if ( abs ( b(i) ) .lt. abs ( b(j) ) ) then
          i = j
        end if
      end do

      if ( m .eq. 0 ) go to 90
      if ( m .eq. i ) go to 100
      if ( d .ge. c ) go to 100
90    m = i
      d = c
      go to 30
100   c = c * a(3)
      c = max ( c, 1.0E+00 )
      pcon = c
      return
      end
      function pdet ( iexp, a )

c*********************************************************************72
c
cc PDET computes the determinant of a tridiagonal matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, integer IEXP, a power of ten that is part of the
c    determinant.
c
c    Input, real A(*), factorization information from PFACT.
c
c    Output, real PDET, the mantissa of the determinant.
c    Determinant = PDET * 10.0^IEXP.
c
      implicit none

      real a(*)
      real c
      real d
      real f
      real g
      integer h
      integer i
      integer iexp
      integer l
      integer n
      real pdet

      iexp=0
      pdet=0.0
      d = a(1)

      if ( abs ( d ) .ne. 1235 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor with pfact before',
     &  ' computing determinant'
        return
      end if

      iexp = 0
      if ( d .lt. 0. ) go to 80
      n = a(2)
      if ( n .eq. 1 ) go to 90
      d = 1.0E+00
      f = 2.0E+00**64
      g = 1.0E+00 /f
      h = 64
      l = 1 + 4*n

      do 50 i = 6, l, 4
        if ( a(i-2) .ne. 0. ) d = -d
        d = d*a(i)
20      if ( abs ( d ) .gt. f ) go to 40
30      if ( abs ( d ) .gt. g ) go to 50
        iexp = iexp - h
        d = d*f
        go to 30
40      iexp = iexp + h
        d = d*g
        go to 20
50    continue

      d = d*a(l+1)
      if ( iexp .ne. 0 ) go to 60
      pdet = d
      return
60    if ( d .eq. 0. ) go to 100
      c = alog10 ( abs ( d ) ) + iexp * alog10 ( 2.0E+00 )
      iexp = c
      c = c - iexp
      if ( c .le. 0.0 ) go to 70
      c = c - 1
      iexp = iexp + 1
70    f = 10.0E+00**c
      if ( d .lt. 0. ) f = -f
      pdet = f
      return
80    pdet = 0.
      return
90    pdet = a(6)
      return
100   iexp = 0
      go to 80
      end
      subroutine pfact ( a, la, n )

c*********************************************************************72
c
cc PFACT factors a tridiagonal matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input/output, real A(3+4*N); on input, information defining the matrix.
c    On output, factorization information.
c
c    Input, integer LA, the leading dimension of the array A.
c
c    Input, integer N, the number of rows and columns in the matrix.
c
      implicit none

      real a(*)
      integer g
      integer h
      integer i
      integer j
      integer k
      integer l
      integer la
      integer n
      real o
      real p
      real q
      real r
      real s
      real t

      k = la*n - la + 1
      a(k) = 0.
      k = k + 2
      a(k) = 0.
      if ( la .gt. 3 ) then
        call rpack(a,la,3,n)
      end if
      k = 3*n
      i = 3 + k + n
      j = k - i
      r = 0.
      s = 0.
10    if ( i .eq. 3 ) go to 20
      p = a(i+j)
      a(i) = p
      i = i - 1
      q = a(i+j)
      a(i) = q
      t = abs ( p ) + abs ( q )
      i = i - 1
      q = a(i+j)
      a(i) = q
      r = r + abs ( q )
      s = max ( s, r )
      r = t
      i = i - 2
      j = j + 1
      go to 10
20    continue

      s = max ( s, r )
      a(1) = 1235
      a(2) = n
      a(3) = s
      o = 2.0E+00**(-64)
      t = .5*o
      s = t
30    t = .5*t
      p = s
      s = s + t
      if ( s .ge. o ) go to 40
      if ( s+t .gt. s ) go to 30
40    r = a(6)
      l = k + n
      g = 1
50    g = g + 4
      if ( g .gt. l ) go to 90
      h = g - 1
      i = g + 2
      j = g + 5
      q = a(i)
      if ( abs ( r ) .ge. abs ( q ) ) go to 70
      t = a(g)
      a(g) = a(j)
      a(j) = t
      t = r/q
      k = g + 1
      a(k) = q
      k = j - 1
      s = a(k)
      if ( s .eq. 0. ) go to 60
      if ( s .eq. o ) s = p
      a(k) = -s*t
      a(h) = s
      go to 80
60    a(k) = s
      a(h) = o
      go to 80
70    a(h) = 0.
      if ( r .eq. 0. ) go to 100
      t = q/r
80    r = a(j) - t*a(g)
      a(j) = r
      a(i) = t
      go to 50
90    if ( r .eq. 0. ) a(1) = -1235
      return
100   a(1) = -1235
      go to 50
      end
      subroutine power ( ex, x, ey, y, n, dif, size, ndigit, limit,
     &  mult )

c*********************************************************************72
c
cc POWER
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real a
      real b
      real c
      real d
      real dif
      real e
      complex ex
      complex ey
      real f
      real g
      integer i
      integer j
      integer k
      integer l
      integer limit
      integer m
      external mult
      integer n
      integer n2
      integer ndigit
      real q
      real r
      real s
      real size
      real t
      real x(*)
      real y(*)

      m = n + 1
      n2 = n + n
      q = 10.0E+00**(-2*ndigit)
      l = 0
      t = 0.

      do i = 1, n2, 2
        t = t + x(i)**2
      end do

      if ( t .eq. 0. ) then
        write(*,*) 'starting guess for routine power cannot be zero'
        stop
      end if

      t = 1.0E+00 / sqrt ( t )

      j = 1
      do i = 1, n
        x(i) = t * x(j)
        j = j + 2
      end do

40    do i = 1, n
        y(i) = x(i)
      end do

      call mult ( y, x )

      do i = 1, n
        y(i+n) = y(i)
      end do

      call mult ( y(m), y )
      a = 0.
      do i = 1, n
        a = a + y(i+n)**2
      end do
      if ( a .eq. 0. ) go to 440
      a = 1.0E+00 / sqrt ( a )
      s = 0.
      do i = 1, n
        s = s + x(i)*y(i)
        j = i + n
        y(j) = y(j)*a
      end do
      a = x(1)
      b = 1.0E+00 + abs ( a )
      if ( a .lt. 0. ) go to 90
      a = a + 1.0E+00
      s = s + y(1)
      go to 100
90    a = a - 1.0E+00
      s = s - y(1)
100   s = s/b

      y(1) = y(1) - a*s
      c = 0.
      d = 0.
      do i = 2, n
        t = y(i) - s*x(i)
        y(i) = t
        c = c + t*x(i)
        d = d + t*t
      end do

      if ( d .ne. 0. ) go to 120
      c = x(2)
      d = 1.0E+00
      y(2) = 1.0E+00
120   d = 1.0E+00 / sqrt ( d )
      c = c/b
      y(1) = -c*a*d
      do i = 2, n
        y(i) = d*(y(i)-c*x(i))
      end do

      s = 0.
      t = 0.
      do i = 1, n
        j = i + n
        s = s + y(j)*x(i)
        t = t + y(j)*y(i)
      end do

      d = 0.
      do i = 1, n
        a = y(i+n)
        d = d + (a-s*x(i)-t*y(i))**2
        x(i) = a
      end do

      l = l + 2
      if ( l .ge. limit ) go to 160
      if ( d .gt. q ) go to 40
160   s = 0.
      do i = 1, n
        s = s + x(i)*y(i)
      end do
      a = x(1)
      b = 1.0E+00 + abs ( a )
      if ( a .lt. 0. ) go to 180
      a = a + 1.0E+00
      s = s + y(1)
      go to 190
180   a = a - 1.0E+00
      s = s - y(1)
190   s = s/b
      y(1) = y(1) - a*s

      c = 0.
      d = 0.
      do i = 2, n
        t = y(i) - s*x(i)
        y(i) = t
        c = c + t*x(i)
        d = d + t*t
      end do

      if ( d .ne. 0. ) go to 210
      c = x(2)
      d = 1.0E+00
      y(2) = 1.0E+00
210   c = c/b
      d = 1.0E+00 / sqrt ( d )
      y(1) = -c*a*d
      y(m) = y(1)
      x(m) = x(1)
      do i = 2, n
        j = i + n
        t = d*(y(i)-c*x(i))
        y(i) = t
        y(j) = t
        x(j) = x(i)
      end do
      call mult ( x, x(m) )
      call mult ( y, y(m) )
      a = 0.
      b = 0.
      c = 0.
      d = 0.
      s = 0.
      t = 0.
      do i = 1, n
        j = i + n
        f = x(i)
        g = y(i)
        s = s + f*f
        t = t + g*g
        a = a + f*x(j)
        b = b + g*x(j)
        c = c + f*y(j)
        d = d + g*y(j)
      end do
      if ( t .eq. 0. ) go to 490
      if ( s .eq. 0. ) go to 540

      f = 1.0E+00 / sqrt ( s )
      g = 1.0E+00 / sqrt ( t )
      r = 0.0E+00
      s = 0.0E+00
      t = 0.0E+00
      do i = 1, n
        j = i + n
        r = r + (x(i)-a*x(j))**2
        s = s + (x(i)-a*x(j)-c*y(j))**2
        t = t + (y(i)-b*x(j)-d*y(j))**2
        x(i) = f*x(i)
        y(i) = g*y(i)
      end do

      t = f*f*s + g*g*t
      s = f*f*r
      l = l + 2
      if ( l .ge. limit ) go to 280
      if ( t .le. q ) go to 290
      if ( s .gt. q ) go to 160
250   size = 3*l + 1
260   dif = sqrt ( s )
      ex = cmplx(a,0.)
      ey = (0.,0.)
      j = n2
      i = n
270   k = j - 1
      y(k) = 0.
      y(j) = 0.
      x(k) = x(i)
      x(j) = 0.
      j = j - 2
      i = i - 1
      if ( j .gt. 0 ) go to 270
      return
280   size = 3*l + 2
      if ( s+s .lt. t ) go to 260
      go to 300
290   size = 3*l + 1
300   dif = sqrt ( t )
      s = .5*(d-a)
      t = b*c
      r = s*s + t
      if ( r .lt. 0. ) go to 400
      r = abs ( s ) + sqrt ( r )
      if ( s .lt. 0. ) go to 310
      ex = cmplx(a-t/r,0.)
      ey = cmplx(a+r,0.)
      go to 320
310   ex = cmplx(a+t/r,0.)
      ey = cmplx(a-r,0.)
320   s = cabs ( a - ex ) + abs ( b )
      t = cabs ( d - ex ) + abs ( c )
      if ( s .gt. t ) go to 340
      if ( t .ne. 0. ) go to 330
      e = 0.
      f = 1.0E+00
      a = 1.0E+00
      b = 0.
      go to 380
330   e = c/t
      f = (ex-d)/t
      go to 350
340   e = (ex-a)/s
      f = b/s
350   s = cabs ( a - ey ) + abs ( b )
      t = cabs ( d - ey ) + abs ( c )
      if ( s .gt. t ) go to 370
      if ( t .ne. 0. ) go to 360
      e = 0.
      f = 1.0E+00
      a = 1.0E+00
      b = 0.
      go to 380
360   a = c/t
      b = (ey-d)/t
      go to 380
370   a = (ey-a)/s
      b = b/s
380   j = 0

      do i = m, n2
        j = j + 2
        k = j - 1
        s = f*x(i) + e*y(i)
        t = b*x(i) + a*y(i)
        y(j) = 0.
        x(j) = 0.
        x(k) = s
        y(k) = t
      end do

      return
400   q = a + s
      r = sqrt ( - r )
      ex = cmplx(q,r)
      ey = cmplx(q,-r)
      s = cabs ( a - ex ) + abs ( b )
      t = cabs ( d - ex ) + abs ( c )
      j = 0
      if ( s .gt. t ) go to 420
      b = c/t
      c = (q-d)/t
      a = r/t
      do i = m, n2
        j = j + 2
        k = j - 1
        t = c*x(i) + b*y(i)
        s = a*x(i)
        x(k) = t
        y(k) = t
        x(j) = s
        y(j) = -s
      end do

      return
420   c = b/s
      b = (q-a)/s
      d = r/s
      do i = m, n2
        j = j + 2
        k = j - 1
        t = c*x(i) + b*y(i)
        s = d*y(i)
        x(k) = t
        y(k) = t
        x(j) = s
        y(j) = -s
      end do

      return
440   ex = (0.,0.)
      ey = (0.,0.)
      dif = 0.
      size = 3*l + 1
      do i = 1, n
        if ( y(i) .ne. 0. ) go to 470
      end do
      j = n2
      i = n
460   k = j - 1
      x(k) = x(i)
      x(j) = 0.
      y(k) = 0.
      y(j) = 0.
      j = j - 2
      i = i - 1
      if ( j .gt. 0 ) go to 460
      return
470   j = n2
      i = n
480   k = j - 1
      x(k) = y(i)
      x(j) = 0.
      y(k) = 0.
      y(j) = 0.
      j = j - 2
      i = i - 1
      if ( j .gt. 0 ) go to 480
      return
490   if ( s .ne. 0. ) go to 510
      do i = 1, n
        x(i) = x(i+n)
      end do
      go to 440
510   s = 1.0E+00 / sqrt ( s )
      a = a*s
      do i = 1, n
        x(i) = x(i)*s
      end do
      s = 0.
      do i = 1, n
        s = s + (x(i)-a*x(i+n))**2
      end do
      if ( s .gt. q ) go to 40
      go to 250
540   t = 1.0E+00 / sqrt ( t )
      do i = 1, n
        x(i) = t*y(i)
      end do
      go to 40
      end
      subroutine precg ( x, dif, size, ndigit, limit, b, n, mult, 
     &  pre, w )

c*********************************************************************72
c
cc PRECG solves a linear system using (preconditioned) conjugate gradients.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      integer n

      real b(*)
      real dif
      integer i
      integer limit
      external mult
      integer ndigit
      external pre
      real r
      real s
      real size
      real t
      real w(n,*)
      real x(*)

      call mult ( w, x )

      do i = 1, n
        w(i,1) = b(i) - w(i,1)
      end do

      call pre ( w(1,2), w )

      s = 0.0E+00
      do i = 1, n
        s = s + w(i,1) * w(i,2)
      end do

      if ( s .gt. 0. ) go to 30
      dif = 0.
      size = 0.
      go to 80
30    call mult ( b, w(1,2) )
      t = 0.
      do i = 1, n
        t = t + w(i,2)*b(i)
      end do
      dif = 0.
      size = 0.
      if ( t .eq. 0. ) go to 80

      r = s/t
      do i = 1, n
        t = r*w(i,2)
        dif = dif + abs ( t )
        x(i) = x(i) + t
        size = size + abs ( x(i) )
        w(i,1) = w(i,1) - r*b(i)
      end do

      call pre ( b, w )
      t = 0.
      do i = 1, n
        t = t + b(i)*w(i,1)
      end do
      r = t/s
      do i = 1, n
        w(i,2) = b(i) + r*w(i,2)
      end do

      s = t

80    call stopit ( dif, size, ndigit, limit )

      if ( dif .gt. 0. ) go to 30
      return
      end
      subroutine pseudo ( a, la, q, lq, mq, d, p, lp, mp, r )

c*********************************************************************72
c
cc PSEUDO computes the regularized pseudoinverse,
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      integer la
      integer lp
      integer lq

      real a(la,*)
      real d(*)
      integer i
      integer j
      integer k
      integer mp
      integer mq
      integer n
      real q(lq,*)
      real p(lp,*)
      real r
      real s
      real t
      real y
      real z

      n = min ( mq, mp )
      z = abs ( r )
      y = sqrt ( z )

      do j = 1, mq
        do i = 1, mp
          a(i,j) = 0.0
        end do
      end do

      do 50 k = 1, n

        if ( d(k) .eq. 0. ) go to 50
        if ( d(k) .gt. y ) go to 20
        s = d(k)/(d(k)**2+z)
        go to 30
20      s = 1.0E+00 /(d(k)+z/d(k))
30      do j = 1, mq
          t = s*q(j,k)
          do i = 1, mp
            a(i,j) = a(i,j) + t*p(i,k)
          end do
        end do

50    continue

      return
      end
      subroutine psolve ( x, a, b )

c*********************************************************************72
c
cc PSOLVE solves a tridiagonal linear system.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, real X(N), the solution.
c
c    Input, real A(*), factorization information from PFACT.
c
c    Input, real B(N), the right hand side.
c
      implicit none

      real a(*)
      real b(*)
      integer i
      integer j
      integer k
      integer n
      real s
      real t
      real x(*)

      t = a(1)

      if ( abs ( t ) .ne. 1235 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor with pfact before solving'
        stop
      end if

      n = a(2)
      s = 2.0E+00**(-64)
      if ( t .lt. 0. ) go to 70
      x(1) = b(1)
      j = 1
      k = 3
20    k = k + 4
      if ( j .eq. n ) go to 40
      i = j
      j = j + 1
      x(j) = b(j)
      if ( a(k-3) .eq. 0. ) go to 30
      t = x(j)
      x(j) = x(i)
      x(i) = t
30    x(j) = x(j) - a(k)*x(i)
      go to 20
40    k = k - 1
      x(j) = x(j)/a(k)
50    if ( j .eq. 1 ) return
      j = j - 1
      k = k - 4
      x(j) = (x(j)-a(k-1)*x(j+1))/a(k)
60    if ( j .eq. 1 ) return
      j = j - 1
      k = k - 4
      t = a(k-2)
      if ( t .eq. s ) t = 0.
      x(j) = (x(j)-a(k-1)*x(j+1)-t*x(j+2))/a(k)
      go to 60
70    j = 0
      k = 2
80    k = k + 4
      j = j + 1
      if ( a(k) .ne. 0. ) go to 80
      do i = 1, n
        x(i) = 0.
      end do
      x(j) = 1.0E+00
      go to 50
      end
      subroutine ptrans ( x, a, b )

c*********************************************************************72
c
cc PTRANS solves a transposed tridiagonal linear system.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, real X(N), the solution.
c
c    Input, real A(*), factorization information from PFACT.
c   
c    Input, real B(N), the right hand side.
c   
      implicit none

      real a(*)
      real b(*)
      integer i
      integer j
      integer k
      integer n
      real s
      real t
      real x(*)

      t = a(1)

      if ( abs ( t ) .ne. 1235 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor with pfact before solving'
        stop
      end if

      n = a(2)
      s = 2.0E+00**(-64)
      if ( t .lt. 0. ) go to 70
      k = 1
      x(1) = b(1)/a(6)
      if ( n .eq. 1 ) return
      k = 10
      x(2) = (b(2)-x(1)*a(5))/a(k)
      i = 2
      if ( n .eq. 2 ) go to 40
20    i = i + 1
      k = k + 4
      t = a(k-10)
      if ( t .eq. s ) t = 0.
      x(i) = (b(i)-x(i-1)*a(k-5)-x(i-2)*t)/a(k)
30    if ( i .lt. n ) go to 20
40    k = k + 1
50    if ( i .eq. 1 ) return
      k = k - 4
      j = i
      i = i - 1
      t = x(i) - x(j)*a(k)
      if ( a(k-3) .eq. 0. ) go to 60
      x(i) = x(j)
      x(j) = t
      go to 50
60    x(i) = t
      go to 50
70    i = n
      k = 2 + 4*n
80    if ( a(k) .eq. 0. ) go to 90
      i = i - 1
      k = k - 4
      go to 80
90    continue

      do j = 1, i
        x(j) = 0.
      end do

      x(i) = 1.0E+00
      if ( i .eq. n ) go to 40
      j = i
      i = i + 1
      k = k + 4
      x(i) = -x(j)*a(k-5)/a(k)
110   if ( i .eq. n ) go to 40
      j = i
      i = i + 1
      k = k + 4
      t = a(k-10)
      if ( t .eq. s ) t = 0.
      x(i) = -(x(j)*a(k-5)+x(j-1)*t)/a(k)
      go to 110
      end
      subroutine pvect ( e, x, d, u, n, w )

c*********************************************************************72
c
cc PVECT computes the eigenvectors of a tridiagonal matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real d(*)
      real e
      integer f
      integer g
      integer h
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      real o
      real p
      real q
      real r
      real s
      real t
      real u(*)
      real v
      real w(*)
      real x(*)
      real y
      real z

      if ( n .eq. 1 ) then
        e = d(1)
        x(1) = 1.0E+00
        return
      end if

      m = n - 1
      j = 2

      do i = 1, m
        x(i) = 0.
        w(j) = u(i)
        w(j+1) = d(i) - e
        w(j+2) = u(i)
        j = j + 4
      end do

      w(j) = 0.
      w(j+1) = d(n) - e
      o = 2.0E+00**(-64)
      t = .5*o
      s = t
30    t = .5*t
      p = s
      s = s + t
      if ( s .ge. o ) go to 40
      if ( s+t .gt. s ) go to 30
40    r = w(3)
      v = abs ( r ) + abs ( w(4) )
      f = 4
      l = 4*n - 3
      g = -2
50    g = g + 4
      if ( g .gt. l ) go to 110
      h = g - 1
      i = g + 2
      j = g + 5
      q = w(i)
      y = abs ( q )
      z = abs ( r )
      if ( z .ge. y ) go to 80
      if ( v .ge. y ) go to 60
      v = y
      f = i
60    t = w(g)
      w(g) = w(j)
      w(j) = t
      t = r/q
      k = g + 1
      w(k) = q
      k = j - 1
      s = w(k)
      if ( s .eq. 0. ) go to 70
      if ( s .eq. o ) s = p
      w(k) = -s*t
      w(h) = s
      go to 100
70    w(k) = s
      w(h) = o
      go to 100
80    w(h) = 0.
      if ( v .ge. z ) go to 90
      v = z
      f = i
90    if ( r .eq. 0. ) go to 120
      t = q/r
100   r = w(j) - t*w(g)
      w(j) = r
      w(i) = t
      go to 50
110   if ( abs ( r ) .ge. v ) go to 120
      v = r
      f = j + 1
120   j = f/4
      x(j) = 1.0E+00
      if ( j .eq. 1 ) go to 140
      k = f - 5
      j = j - 1
      x(j) = (x(j)-w(k-1)*x(j+1))/w(k)
130   if ( j .eq. 1 ) go to 140
      j = j - 1
      k = k - 4
      t = w(k-2)
      if ( t .eq. o ) t = 0.
      x(j) = (x(j)-w(k-1)*x(j+1)-t*x(j+2))/w(k)
      go to 130
140   if ( v .eq. 0. ) go to 230

      s = 0.
      do i = 1, n
        t = abs ( x(i) )
        if ( t .gt. s ) s = t
      end do

      r = 0.
      s = 1.0E+00 /s
      do i = 1, n
        t = s*x(i)
        r = r + t*t
        x(i) = t
      end do

      k = 0
      j = 1
      y = x(1)
170   k = k + 4
      i = j
      j = j + 1
      s = w(k)
      w(k) = y
      y = x(j)
      if ( w(k-3) .eq. 0. ) go to 180
      t = x(j)
      x(j) = x(i)
      x(i) = t
180   x(j) = x(j) - s*x(i)
      if ( j .lt. n ) go to 170
      s = x(j)/w(k+3)
      x(j) = s
      t = abs ( s )
      v = s*y
      j = j - 1
      k = k - 1
      s = (x(j)-w(k-1)*s)/w(k)
      x(j) = s
      if ( abs ( s ) .gt. t ) t = abs ( s )
      v = v + s*w(k+1)
190   if ( j .eq. 1 ) go to 200
      j = j - 1
      k = k - 4
      z = w(k-2)
      if ( z .eq. o ) z = 0.
      s = (x(j)-w(k-1)*s-z*x(j+2))/w(k)
      x(j) = s
      if ( abs ( s ) .gt. t ) t = abs ( s )
      v = v + s*w(k+1)
      go to 190
200   if ( v .ne. 0. ) v = r/v
      s = 0.
      t = 1.0E+00 /t
      z = 0.
      do i = 1, n
        s = s + (x(i)*v)**2
        z = z + (t*x(i))**2
      end do
      t = t / sqrt ( z )
      do i = 1, n
        x(i) = t*x(i)
      end do
      if ( r+r .ge. s ) e = e + v
      return
230   t = 0.
      do i = 1, n
        s = abs ( x(i) )
        if ( s .gt. t ) t = s
      end do
      t = 1.0E+00 /t
      z = 0.
      do i = 1, n
        z = z + (t*x(i))**2
      end do
      t = t / sqrt ( z )
      do i = 1, n
        x(i) = t*x(i)
      end do

      return
      end
      subroutine pvert ( v, lv, a )

c*********************************************************************72
c
cc PVERT inverts a tridiagonal matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, real V(LV,N), the NxN inverse matrix.
c
c    Input, integer LV, the leading dimension of the array V.
c
c    Input, real A(*), factorization information from PFACT.
c
      implicit none

      integer lv

      real a(*)
      integer i
      integer j
      integer k
      integer l
      integer n
      real s
      real t
      real v(lv,*)

      t = a(1)

      if ( abs ( t ) .ne. 1235 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor with pfact before inverting'
        stop
      end if

      if ( t .le. 0. ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: matrix has no inverse'
        stop
      end if

      n = a(2)
      s = 2.0E+00**(-64)

      do 110 l = 1, n

        if ( l .gt. 1 ) go to 30
        v(l,l) = 1.0E+00
        j = l
        k = 3
        go to 60
30      j = l - 1
        k = 4*j - 1
        do i = 1, j
          v(i,l) = 0.
        end do
        k = k + 4
        i = j
        j = j + 1
        v(j,l) = 1.0E+00
        if ( a(k-3) .eq. 0. ) go to 50
        t = v(j,l)
        v(j,l) = v(i,l)
        v(i,l) = t
50      v(j,l) = v(j,l) - a(k)*v(i,l)
60      k = k + 4
        if ( j .eq. n ) go to 80
        i = j
        j = j + 1
        v(j,l) = 0.
        if ( a(k-3) .eq. 0. ) go to 70
        t = v(j,l)
        v(j,l) = v(i,l)
        v(i,l) = t
70      v(j,l) = v(j,l) - a(k)*v(i,l)
        go to 60
80      k = k - 1
        v(j,l) = v(j,l)/a(k)
90      if ( j .eq. 1 ) go to 110
        j = j - 1
        k = k - 4
        v(j,l) = (v(j,l)-a(k-1)*v(j+1,l))/a(k)
100     if ( j .eq. 1 ) go to 110
        j = j - 1
        k = k - 4
        t = a(k-2)
        if ( t .eq. s ) t = 0.
        v(j,l) = (v(j,l)-a(k-1)*v(j+1,l)-t*v(j+2,l))/a(k)
        go to 100

110   continue

      return
      end
      subroutine pwk ( e, d, u, n, w )

c*********************************************************************72
c
cc PWK
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real a
      real b
      real c
      real d(*)
      real e(*)
      real g
      real h
      integer i
      integer j
      integer k
      integer kk
      integer l
      integer m
      integer n
      real oc
      real og
      real p
      real r
      real s
      real t
      real u(*)
      real v
      real w(*)

      e(n) = d(n)

      if ( n .eq. 1 ) then
        return
      end if

      k = 1
      m = n - 1
      l = n + 1
      t = 1.0E+00

      do i = 1, m

        e(i) = d(i)
        j = l - i
        w(j) = u(j-1)

        if ( w(j) .lt. 0.0E+00 ) then
          w(1) = 2.0E+00
          return
        else if ( w(j) .eq. 0.0E+00 ) then
          if ( t .eq. 1.0E+00 ) k = j
          t = 0.0E+00
        end if

      end do

      w(1) = 0.0E+00
      t = 1.0E+00
30    t = t / 2.0E+00
      s = 1.0E+00 + t
      if ( s .gt. 1.0E+00 ) go to 30
      l = n
40    m = 0
      kk = 0
      v = 0.0E+00
50    p = w(l)
      s = sqrt ( p )
      h = e(l)
      r = s + abs ( h )
      v = max ( v, r )
      if ( s .gt. t * v ) go to 100
      if ( s .le. t * r ) go to 60
      v = r
      if ( kk .eq. 0 ) go to 90
60    l = l - 1
      if ( l .eq. 1 ) return
      if ( k .le. l ) go to 40
      k = 1
      do i = 2, l
        if ( w(i) .eq. 0.0E+00 ) k = i
      end do
      go to 40
80    w(l) = s * p
      e(l) = g + h
      go to 50
90    kk = 1
100   m = m + 1
      if ( m .gt. 40 ) go to 130
      r = 0.5E+00 * ( e(l-1) - h )
      if ( r .lt. 0.0E+00 ) then
        h = h + p / ( abs ( r ) + sqrt ( r * r + p ) )
      else
        h = h - p / ( abs ( r ) + sqrt ( r * r + p ) )
      end if
      c = 1.0E+00
      s = 0.0E+00
      i = k
      g = e(i) - h
      p = g * g
110   j = i
      i = i + 1
      if ( i .gt. l ) go to 80
      b = w(i)
      r = p + b
      a = s * r
      if ( a .eq. 0.0E+00 ) k = j
      w(j) = a
      oc = c
      c = p/r
      s = b/r
      a = e(i)
      og = g
      g = c * ( a - h ) - s * og
      e(j) = og + ( a - g )
      if ( c .eq. 0.0E+00 ) go to 120
      p = g * g / c
      go to 110
120   p = oc * b
      go to 110
130   w(1) = 1.0E+00
      return
      end
      subroutine px ( y, z, k, n, n5 )

c*********************************************************************72
c
cc PX
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      integer n5
      real t
      real y(*)
      real z(*)

      if ( k .eq. 0 ) then
        return
      end if

      l = n5

      do m = 1, k

        j = k - m + 1
        l = l - n + j

        t = 0.
        do i = j, n
          t = t + z(l+i) * y(i)
        end do

        do i = j, n
           y(i) = y(i) - t * z(l+i)
        end do

      end do

      return
      end
      subroutine qad ( b, c )

c*********************************************************************72
c
cc QAD determines the roots of a quadratic polynomial.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real amag
      complex b
      complex c
      complex r
      complex s
      real t

      t = amag ( b )
      r = cmplx(4.0,0.0)*c

      if ( t .le. sqrt ( amag ( r ) ) )then
        r = csqrt ( b * b - r )
      else
        s = 1.0E+00/t
        r = t * csqrt ( ( s * b )**2 - r * s * s )
      end if

      s = r + b
      r = r - b
      t = amag ( r )

      if ( t .lt. amag ( s ) )then
        c = -(c+c)/s
        b = cmplx(-.5,0.0)*s
      else
        if ( t .eq. 0.0 ) return
        c = (c+c)/r
        b = cmplx(.5,0.0)*r
      end if

      return
      end
      subroutine qr ( a, la, m, n )

c*********************************************************************72
c                            
cc QR computes the QR factorization of a general matrix with column pivoting.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c                                             
c         input:                                              
c                                                             
c              a     --array containing matrix                
c                      (length at least 3+n+mn + min ( m, n ) )    
c                                                             
c              la    --leading (row) dimension of array a     
c                                                             
c              m     --number of rows in coefficient matrix   
c                                                             
c              n     --number of columns in coeff. matrix     
c                                                             
c         output:                                             
c                                                             
c              a     --factored matrix                        
c                
      implicit none
                                      
      real a(*)
      integer b
      integer c
      integer d
      integer e
      integer f
      integer g
      integer h
      integer i
      integer j
      integer k
      integer l
      integer la
      integer m
      integer n
      integer o
      integer p
      integer q
      real r
      real s
      real t
      real u
      real v
      real w
      integer z

      if ( la .gt. m ) then
        call rpack ( a, la, m, n )
      end if

      u = 1.0E+00
10    u = .5*u
      t = 1 + u
      if ( t .gt. 1.0E+00 ) go to 10
      u = sqrt ( 40. * u )
      o = m + 1
      i = 2 + n*o
      f = i
      if ( m .lt. n ) f = 2 + m*o
      d = i - m
      j = 2 + n
      h = i + 1
      c = h
      l = h + n
      g = h
      v = 0.
c
c  Compute 2-norm of each column.
c
20    s = 0.
      r = s
      k = i - m
30    t = a(i-j)
      a(i) = t
      i = i - 1
      if ( t .ne. 0. ) go to 40
      if ( i .gt. k ) go to 30
      go to 70
40    s = abs ( t )
      r = 1.0E+00
      if ( i .eq. k ) go to 70
50    t = a(i-j)
      a(i) = t
      i = i - 1
      if ( abs ( t ) .gt. s ) go to 60
      r = r + (t/s)**2
      if ( i .gt. k ) go to 50
      go to 70
60    r = 1 + r*(s/t)**2
      s = abs ( t )
      if ( i .gt. k ) go to 50
70    w = s * sqrt ( r )
      a(g) = w
      a(l) = u*w
      if ( w .lt. v ) go to 80
      v = w
      p = g
80    l = l - 1
      g = i
      i = i - 1
      j = j - 1
      k = i - m
      if ( k .gt. 2 ) go to 20
      a(1) = 3230
      a(2) = m
      a(3) = n
      k = 4
      g = m
c
c  start factorization
c
90    e = k + g
      l = p - m
      j = p - e
      h = h + 1
      a(h+j/o) = a(h)
      a(h) = (p-3)/o
c
c  interchange columns
c
      do i = l, p
        t = a(i)
        q = i - j
        a(i) = a(q)
        a(q) = t
      end do

      s = 0.
      t = a(e)
      if ( t .eq. 0. ) go to 240
      if ( k .eq. f ) return
      e = e - 1
      do i = k, e
        s = s + (a(i)/t)**2
      end do
      if ( s .eq. 0. ) go to 240
      s = t * sqrt ( s )
      t = a(k)
      a(k) = s
      if ( t .ge. 0. ) a(k) = -s
      r = 1.0E+00 / sqrt ( s * ( s + abs ( t ) ) )
      i = e
c
c  Store Householder matrix.
c
120   a(i+1) = r*a(i)
      i = i - 1
      if ( i .gt. k ) go to 120
      if ( t .lt. 0. ) s = -s
      a(k+1) = r*(t+s)
      if ( k .gt. d ) return
      j = k
      e = -1
      z = h
      g = g - 1
      v = 0.
130   if ( j .gt. d ) go to 230
      j = j + o
      e = e + o
      l = j + g
      z = z + 1
      b = l + 1
      s = a(b)
      t = 0.
      if ( s .eq. 0. ) go to 160
      do i = j, l
        t = t + a(i)*a(i-e)
      end do
c
c  update factorization
c
      do i = j, l
        a(i) = a(i) - t*a(i-e)
      end do

      t = s*sqrt ( abs ( 1.0E+00 - (a(j)/s)**2 ) )
      a(b) = t
      if ( t .lt. a(z) ) go to 170
160   if ( t .lt. v ) go to 130
      v = t
      p = b
      go to 130
170   i = j + 1
c
c  compute column 2-norm
c
      s = 0.
      r = s
      do q = i, l
        if ( a(q) .ne. 0. ) go to 190
      end do
      go to 220

190   s = abs ( a(q) )
      do i = q, l
        t = abs ( a(i) )
        if ( t .le. s ) then
          r = r + (t/s)**2
        else
          r = 1 + r*(s/t)**2
          s = t
        end if
      end do

220   t = s * sqrt ( r )
      a(b) = t
      a(z) = u*t
      if ( t .lt. v ) go to 130
      v = t
      p = b
      go to 130
230   k = k + o + 1
      if ( k .lt. f ) go to 90
      if ( m .ge. n ) go to 90
      return
240   a(h) = 0.
      a(1) = -3230
      return
      end
      subroutine quasi ( x, h, lh, n, dif, size, ndigit, limit, 
     &  func, w )

c*********************************************************************72
c
cc QUASI uses a quasi-Newton method to solve a nonlinear system.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      integer lh
      integer n

      real dif
      external func
      real h(lh,*)
      integer i
      integer j
      integer limit
      integer ndigit
      real s
      real size
      real t
      real u
      real w(n,*)
      real x(*)

      do i = 1, n
        w(i,1) = 0.
      end do

      go to 30

20    call addchg ( dif, size, x, w, n )
      call stopit(dif,size,ndigit,limit)
      if ( dif .le. 0. ) return

30    do i = 1, n
        w(i,2) = w(i,1)
        w(i,1) = 0.
      end do

      call func ( w(1,3), x )

      do j = 1, n
        t = w(j,3)
        do i = 1, n
          w(i,1) = w(i,1) - t*h(i,j)
        end do
      end do

      t = 0.
      u = 0.
      do i = 1, n
        s = w(i,2)
        t = t + s*w(i,1)
        u = u + s*s
      end do

      t = u - t

      if ( t .eq. 0. ) go to 20

      do j = 1, n

        s = 0.0E+00
        do i = 1, n
          s = s + h(i,j)*w(i,2)
        end do

        s = s/t
        do i = 1, n
          h(i,j) = h(i,j) + s*w(i,1)
         end do

      end do

      u = u/t

      do i = 1, n
        w(i,1) = u*w(i,1)
      end do

      go to 20
      end
      function r4mat_norm_l1 ( m, n, a )

c*********************************************************************72
c
cc R4MAT_NORM_L1 returns the matrix L1 norm of an R4MAT.
c
c  Discussion:
c
c    An R4MAT is an array of R4 values.
c
c    The matrix L1 norm is defined as:
c
c      R4MAT_NORM_L1 = max ( 1 <= J <= N )
c        sum ( 1 <= I <= M ) abs ( A(I,J) ).
c
c    The matrix L1 norm is derived from the vector L1 norm, and
c    satisifies:
c
c      r4vec_norm_l1 ( A * x ) <= r4mat_norm_l1 ( A ) * r4vec_norm_l1 ( x ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows in A.
c
c    Input, integer N, the number of columns in A.
c
c    Input, real A(M,N), the matrix whose L1 norm is desired.
c
c    Output, real R4MAT_NORM_L1, the L1 norm of A.
c
      implicit none

      integer m
      integer n

      real a(m,n)
      integer i
      integer j
      real r4mat_norm_l1
      real temp

      r4mat_norm_l1 = 0.0E+00

      do j = 1, n
        temp = 0.0E+00
        do i = 1, m
          temp = temp + abs ( a(i,j) )
        end do
        r4mat_norm_l1 = max ( r4mat_norm_l1, temp )
      end do

      return
      end
      subroutine r4mat_print ( m, n, a, title )

c*********************************************************************72
c
cc R4MAT_PRINT prints an R4MAT.
c
c  Discussion:
c
c    An R4MAT is an array of R4 values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows in A.
c
c    Input, integer N, the number of columns in A.
c
c    Input, real A(M,N), the matrix.
c
c    Input, character*(*) TITLE, a title.
c
      implicit none

      integer m
      integer n

      real a(m,n)
      character*(*) title

      call r4mat_print_some ( m, n, a, 1, 1, m, n, title )

      return
      end
      subroutine r4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi,
     &  title )

c*********************************************************************72
c
cc R4MAT_PRINT_SOME prints some of an R4MAT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, real A(M,N), an M by N matrix to be printed.
c
c    Input, integer ILO, JLO, the first row and column to print.
c
c    Input, integer IHI, JHI, the last row and column to print.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer incx
      parameter ( incx = 5 )
      integer m
      integer n

      real a(m,n)
      character * ( 14 ) ctemp(incx)
      integer i
      integer i2hi
      integer i2lo
      integer ihi
      integer ilo
      integer inc
      integer j
      integer j2
      integer j2hi
      integer j2lo
      integer jhi
      integer jlo
      character * ( * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )

      if ( m .le. 0 .or. n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  (None)'
        return
      end if

      do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

        j2hi = j2lo + incx - 1
        j2hi = min ( j2hi, n )
        j2hi = min ( j2hi, jhi )

        inc = j2hi + 1 - j2lo

        write ( *, '(a)' ) ' '

        do j = j2lo, j2hi
          j2 = j + 1 - j2lo
          write ( ctemp(j2), '(i7,7x)') j
        end do

        write ( *, '(''  Col   '',5a14)' ) ( ctemp(j), j = 1, inc )
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) ' '

        i2lo = max ( ilo, 1 )
        i2hi = min ( ihi, m )

        do i = i2lo, i2hi

          do j2 = 1, inc

            j = j2lo - 1 + j2

            write ( ctemp(j2), '(g14.6)' ) a(i,j)

          end do

          write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

        end do

      end do

      return
      end
      subroutine r4vec_print ( n, a, title )

c*********************************************************************72
c
cc R4VEC_PRINT prints an R4VEC.
c
c  Discussion:
c
c    An R4VEC is an array of R4's.
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
c    Input, integer N, the number of components of the vector.
c
c    Input, real A(N), the vector to be printed.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer n

      real a(n)
      integer i
      character ( len = * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
      end do

      return
      end
      subroutine r4vec_uniform_01 ( n, seed, r )

c*********************************************************************72
c
cc R4VEC_UNIFORM_01 returns a unit pseudorandom R4VEC.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 March 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Second Edition,
c    Springer, 1987,
c    ISBN: 0387964673,
c    LC: QA76.9.C65.B73.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, December 1986, pages 362-376.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley, 1998,
c    ISBN: 0471134031,
c    LC: T57.62.H37.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, Number 2, 1969, pages 136-143.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, real R(N), the vector of pseudorandom values.
c
      implicit none

      integer n

      integer i
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer k
      integer seed
      real r(n)

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R4VEC_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      do i = 1, n

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed .lt. 0 ) then
          seed = seed + i4_huge
        end if

        r(i) = real ( seed ) * 4.656612875E-10

      end do

      return
      end
      subroutine rad ( p, q, r, k, n, a )

c*********************************************************************72
c
cc RAD estimates the radius of smallest circle containing a zero for CZERO.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real a
      real amag
      real c
      real d
      real f
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      complex p(*)
      complex q(*)
      complex r(*)
      real s
      complex t

      save j

      data j /1/

      if ( j .eq. k ) go to 20
      j = 1
      l = 1

      do i = 2, n
        q(l) = l*p(i)/n
        l = i
      end do

      q(n) = 1
      c = n * cabs ( p(1) )
      f = 1.0E+00/n
      a = cabs ( p(1) )**(1.0E+00/ real ( n ) )
      go to 40
20    do i = 1, n
        q(i) = r(i)
      end do

40    l = j + k
50    if ( amag ( q(1) ) .eq. 0.0 ) go to 90
      j = j+1
      t = p(1)/q(1)
      do i = 2, n
        q(i-1) = p(i) - t*q(i)
      end do
      s = 1.0E+00/ real ( j )
      f = f**(1.0E+00-s) * cabs ( t )**s
      d = f*(c/cabs ( q(1) ) )**s
      a = min ( a, d )
70    if ( j .lt. l ) go to 50
      k = j
      do i = 1, n
        r(i) = q(i)
      end do
      return
90    m = j
100   j = j + 1
      do i = 2, n
        q(i-1) = q(i)
      end do
      q(n) = (0.0,0.0)
      if ( amag ( q(1) ) .eq. 0.0 ) go to 100
      s = 1.0E+00/j
      f = f**( real ( m ) * s )
      d = f*(c/ cabs ( q(1) ) )**s
      a = min ( a, d )
      t = p(1)/q(1)
      do i = 2, n
        q(i-1) = p(i) - t*q(i)
      end do
      q(n) = (1.0E+00,0.0)
      f = f * cabs ( t )**s
      go to 70
      end
      subroutine rfn ( p, q, r, z, n, s, v, b, c )

c*********************************************************************72
c
cc RFN carries out a Newton step for CZERO.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real a
      real amag
      real b
      real c
      integer i
      complex m
      integer n
      complex p(*)
      complex q(*)
      complex r(*)
      complex s
      complex t
      complex u
      complex v
      complex w
      complex z(*)

      a = cabs ( s )
      if ( a .gt. 1.0E+00 ) go to 20
      v = (1.0E+00,0.0)
      u = (0.0,0.0)
      w = (0.0,0.0)
      c = 1.0E+00
      i = n

10    v = v*s + p(i)
      u = u*s + q(i)
      w = w*s + r(i)
      c = c*a + real ( z(i) )
      i = i - 1
      if ( i .gt. 0 ) go to 10
      go to 40

20    t = cmplx(1.0E+00,0.0)/s
      a = 1.0E+00 /a
      v = p(1)
      u = q(1)
      w = r(1)
      c = real ( z(1) )

      do i = 2, n
        v = v*t + p(i)
        u = u*t + q(i)
        w = w*t + r(i)
        c = c*a + real ( z(i) )
      end do

      v = v*t + cmplx(1.0E+00,0.0)
      u = u*t
      w = w*t
      c = c*a + 1.0E+00
40    b = amag ( v )
      if ( b .eq. 0.0 ) return
      a = amag ( u )
      if ( a .eq. 0.0 ) return
      t = v/u
      a = amag ( w )

      if ( a .eq. 0.0 ) then
        s = s - t
        return
      end if

      u = u/w
      m = u - t
      a = amag ( m )

      if ( a .eq. 0.0 ) then
        return
      end if

      m = u/m
      s = s - m*t
      return
      end
      function root ( y, z, t, f )

c*********************************************************************72
c
cc ROOT solves a scalar nonlinear equation.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real a
      real b
      real c
      real d
      real d2
      real d23
      real d3
      real d34
      real d4
      real d42
      real e
      external f
      real f
      real fa
      real fb
      real fc
      real fd
      real fl
      real fr
      real l
      real p
      real p2
      real p3
      real p4
      real q
      real r
      real root
      real s
      real t
      real u
      real v
      real w
      real y
      real z

      l = y
      r = z

      if ( l .gt. r ) then
        a = l
        l = r
        r = a
      end if

      fl = f(l)
      fr = f(r)
      if ( sign(fl,fr) .eq. fl ) go to 220
      q = fl / abs ( fl )
      s = t + t
      u = 1.0E+00
20    u = .5*u
      a = 1.0E+00 + u
      if ( a .gt. 1.0E+00 ) go to 20
      u = 5.*u
      v = .5*u
30    e = r - l
      if ( e .le. s ) go to 210
      if ( e .le. u*( abs ( l ) + abs ( r ) ) ) go to 210
      if ( abs ( fl ) .gt. abs ( fr ) ) go to 40
      a = l
      fa = fl
      b = r
      fb = fr
      go to 50
40    a = r
      fa = fr
      b = l
      fb = fl
50    c = a - fa*(a-b)/(fa-fb)
      p = c
      w = max ( t, v * ( abs ( l ) + abs ( r ) ) )
      if ( abs ( c - l ) .lt. w ) c = l + w
      if ( abs ( c - r ) .lt. w ) c = r - w
      fc = f(c)
      if ( sign(fc,q) .eq. fc ) go to 60
      r = c
      fr = fc
      go to 70
60    l = c
      fl = fc
70    w = r - l
      if ( w .le. s ) go to 250
      if ( w .le. u*( abs ( l ) + abs ( r ) ) ) go to 250
      if ( abs ( fc ) .ge. abs ( fb ) ) go to 90
      if ( abs ( fc ) .gt. abs ( fa ) ) go to 80
      w = c
      c = b
      b = a
      a = w
      w = fc
      fc = fb
      fb = fa
      fa = w
      go to 90
80    w = c
      c = b
      b = w
      w = fc
      fc = fb
      fb = w
90    if ( a .lt. l ) go to 190
      if ( a .gt. r ) go to 190

      call inp ( d, a, b, c, fa, fb, fc, l, r )

100   p = d
      w = max ( t, v * ( abs ( l ) + abs ( r ) ) )
      if ( abs ( l - d ) .lt. w ) go to 110
      if ( abs ( r - d ) .gt. w ) go to 130
110   if ( d+d .gt. l+r ) go to 120
      d = l + w
      go to 130
120   d = r - w
130   if ( d .le. l ) go to 190
      if ( d .ge. r ) go to 190
      e = .5*e
      if ( e .lt. abs ( a - d ) ) go to 190
      fd = f(d)
      if ( sign(fd,q) .eq. fd ) go to 140
      r = d
      fr = fd
      go to 150
140   l = d
      fl = fd
150   w = r - l
      if ( w .le. s ) go to 250
      if ( w .le. u*( abs ( l ) + abs ( r ) ) ) go to 250
      w = abs ( fd )
      if ( w .le. abs ( fa ) ) go to 170
      if ( w .le. abs ( fb ) ) go to 160
      if ( w .ge. abs ( fc ) ) go to 180
      w = d
      d = c
      c = w
      w = fd
      fd = fc
      fc = w
      go to 180
160   w = d
      d = c
      c = b
      b = w
      w = fd
      fd = fc
      fc = fb
      fb = w
      go to 180
170   w = d
      d = c
      c = b
      b = a
      a = w
      w = fd
      fd = fc
      fc = fb
      fb = fa
      fa = w
180   d4 = d - a
      d3 = c - a
      d2 = b - a
      d34 = c - d
      d42 = d - b
      d23 = b - c
      p2 = 0.
      p3 = 0.
      p4 = 0.
      if ( d34 .ne. 0. ) p2 = 1.0E+00 /d34
      if ( d42 .ne. 0. ) p3 = 1.0E+00 /d42
      if ( d23 .ne. 0. ) p4 = 1.0E+00 /d23
      p2 = (fb-fa)/(d2*(1.0E+00+(d2/d3)*p2*d42+(d2/d4)*p2*d23))
      p3 = (fc-fa)/(d3*(1.0E+00+(d3/d2)*p3*d34+(d3/d4)*p3*d23))
      p4 = (fd-fa)/(d4*(1.0E+00+(d4/d2)*p4*d34+(d4/d3)*p4*d42))
      p2 = p2 + p3 + p4
      if ( p2 .eq. 0. ) go to 190
      d = a - fa/p2
      go to 100
190   d = .5*(l+r)
      fd = f(d)
      if ( sign(fd,q) .eq. fd ) go to 200
      r = d
      fr = fd
      go to 30
200   l = d
      fl = fd
      go to 30
210   root = .5*(l+r)
      return
220   if ( fl .eq. 0. ) go to 230
      if ( fr .eq. 0. ) go to 240
      write ( *, '(a)' ) ' '
      write(*,*) 'error: function has same sign at both starting points'
      stop
230   root = l
      return
240   root = r
      return
250   if ( p .lt. l ) go to 210
      if ( p .gt. r ) go to 210
      root = p
      return
      end
      subroutine rpack ( a, la, m, n )

c*********************************************************************72
c
cc RPACK packs a rectangular matrix into contiguous storage.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input/output, real A(*).
c    On input A is assumed to be an MxN matrix stored as an LAxN array.
c    On output, A has been packed as an MxN array.
c
c    Input, integer LA, the leading dimension of the array.
c    M <= LA.
c
c    Input, integer M, N, the number of rows and columns in the matrix.
c
      implicit none

      integer la
      integer n

      real a(la*n)
      integer h
      integer i
      integer j
      integer jhi
      integer jlo
      integer k
      integer l
      integer m
      integer o

      h = la - m

      if ( h .eq. 0 ) then
        return
      end if

      if ( h .lt. 0 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: la argument in rpack must be .ge. m argument'
        stop
      end if

      i = 0
      k = 1
      l = m
      o = m * n

20    continue

      if ( l .eq. o ) then
        jlo = m*n+1
        jhi = la*n
        do j = jlo, jhi
          a(j) = 0.0E+00
        end do
        return
      end if

      i = i + h
      k = k + m
      l = l + m
      do j = k, l
        a(j) = a(i+j)
      end do

      go to 20
      end
      subroutine rsolve ( x, b, q, lq, mq, d, p, lp, mp, r )

c*********************************************************************72
c
cc RSOLVE computes the regularized solution to a linear system.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, real X(N), the solution.
c
c    Input, real B(N), the right hand side.
c
      implicit none

      integer lp
      integer lq

      real b(*)
      real d(*)
      integer i
      integer j
      integer mp
      integer mq
      integer n
      real p(lp,*)
      real q(lq,*)
      real r
      real t
      real x(*)
      real y
      real z

      n = min ( mq, mp )
      z = abs ( r )
      y = sqrt ( z )

      do i = 1, mp
        x(i) = 0.0
      end do

      do 60 j = 1, n
        if ( d(j) .eq. 0. ) go to 60
        t = 0.
        do i = 1, mq
          t = t + q(i,j)*b(i)
        end do
        if ( d(j) .gt. y ) go to 30
        t = t*(d(j)/(d(j)**2+z))
        go to 40
30      t = t/(d(j)+z/d(j))
40      do i = 1, mp
          x(i) = x(i) + t*p(i,j)
        end do
60    continue

      return
      end
      subroutine runpack ( a, la, m, n )

c*********************************************************************72
c
cc RUNPACK reverses the operation of RPACK by unpacking a rectangular matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input/output, real A(*).
c    On input, the MxN entries of A were packed into the front of an LAxN array.
c    On output, the MxN entries of A have been unpacked into the LAxN array.
c
c    Input, integer LA, the leading dimension of the array.
c    M <= LA.
c
c    Input, integer M, N, the number of rows and columns in the matrix.
c
      implicit none

      real a(*)
      integer i
      integer ii
      integer j
      integer jj
      integer la
      integer m
      integer n

      if ( la .eq. n ) then
        return
      end if

      if ( la .lt. m ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: la in runpack must be .ge. m argument!'
        return
      end if

      do jj = 1, n
        j = n + 1 - jj
        do ii = 1, m
          i = m + 1 - ii
          a((j-1)*la+i) = a((j-1)*m+i)
        end do
      end do

      do j = 1, n
        do i = m + 1, la
          a((j-1)*la+i) = 0.0
        end do
      end do

      return
      end
      subroutine rvc ( f, x, y, v, k )

c*********************************************************************72
c
cc RVC
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real a
      real b
      real c
      real d
      real f
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      real r
      real s
      real v(*)
      real x(*)
      real y(*)

      l = k
      j = (k*(k+3))/2 - 1
      if ( l .eq. 0 ) return
10    m = l - 1
      d = f - v(j)
      if ( m .eq. 0 ) go to 50
      b = v(j-l)
      if ( b .eq. 0. ) go to 50
      a = f - v(j-l-1)
      c = v(j-1)
      s = a*d - b*c
      if ( s .ne. 0. ) go to 30
20    write(*,*) 'more than two eigenvalues with the same magnitude'
      write(*,*) 'have been detected. routine mpower must stop'
      stop
30    r = (c*y(l)+d*y(m))/s
      s = (b*y(m)+a*y(l))/s
      x(m) = r
      x(l) = s
      n = j - l
      m = n - l
      j = m - 1
      l = l - 2
      if ( l .eq. 0 ) return
      do i = 1, l
        y(i) = y(i) + r*v(m+i) + s*v(n+i)
      end do
      go to 10
50    if ( d .eq. 0. ) go to 20
      r = y(l)/d
      x(l) = r
      n = j - l
      j = n - 1
      l = m
      if ( l .eq. 0 ) return
      do i = 1, l
        y(i) = y(i) + r*v(n+i)
      end do
      go to 10
      end
      subroutine scl ( d, u, n, q, lq, mq, p, lp, mp, e, f, b, j, k, 
     &  jl, jr )

c*********************************************************************72
c
cc SCI
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      integer lp
      integer lq

      real b
      real d(*)
      real e(*)
      real f(*)
      integer h
      integer i
      integer j
      integer jl
      integer jr
      integer k
      integer l
      integer m
      integer mp
      integer mq
      integer n
      real p(lp,*)
      real q(lq,*)
      real r
      real s
      real t
      real u(*)
      real v

      t = 1.0E+00

      do i = j, k
        if ( abs ( e(i) ) .lt. t ) t = abs ( e(i) )
        if ( abs ( f(i) ) .lt. t ) t = abs ( f(i) )
      end do

      if ( t .gt. b ) return
      r = e(j)
      r = sign ( sqrt ( abs ( r ) ), r )
      e(j) = 1.0E+00
      s = f(j)
      s = sign ( sqrt ( abs ( s ) ), s )
      f(j) = 1.0E+00
      d(j) = r*d(j)*s
      t = r
      if ( jl .eq. 1 ) go to 20
      if ( jr .ne. 1 ) go to 40
      t = s
20    do i = 1, mq
        q(i,j) = q(i,j)*t
      end do
40    t = s
      if ( jr .eq. 2 ) go to 50
      if ( jl .ne. 2 ) go to 70
      t = r
50    do i = 1, mp
        p(i,j) = p(i,j)*t
      end do
70    l = j + 1
      h = j

      do 130 m = l, k

        t = e(m)
        t = sign ( sqrt ( abs ( t ) ), t )
        e(m) = 1.0E+00
        s = f(m)
        s = sign ( sqrt ( abs ( s ) ), s )
        f(m) = 1.0E+00
        d(m) = s*d(m)*t
        u(h) = r*u(h)*s
        h = m
        r = t
        v = t
        if ( jl .eq. 1 ) go to 80
        if ( jr .ne. 1 ) go to 100
        v = s
80      do i = 1, mq
          q(i,m) = q(i,m)*v
        end do
100     v = s
        if ( jr .eq. 2 ) go to 110
        if ( jl .ne. 2 ) go to 130
        v = t
110     do i = 1, mp
          p(i,m) = p(i,m)*v
        end do

130   continue

      return
      end
      function scon ( a, b )

c*********************************************************************72
c
cc SCON estimates the condition number of a symmetric matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c    William Hager,
c    Condition Estimates,
c    SIAM Journal on Scientific and Statistical Computing,
c    Volume 5, Number 2, June 1984, pages 311-316.
c
c  Parameters:
c
c    Input, real A(*), the factor information from SFACT.
c
c    Workspace, real B(N).
c
c    Output, real SCON, the condition number of the matrix.
c
      implicit none

      real a(*)
      real b(*)
      real c
      real d
      integer i
      integer j
      integer m
      integer n
      real scon

      c = a(1)

      if ( abs ( c ) .ne. 1233 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor with SFACT before',
     &  'estimating condition'
        stop
      end if

      if ( c .lt. 0.0E+00 ) then
        scon = -1.0E+00
        return
      end if

      m = 0
      n = a(2)
      c = 1.0E+00 /a(2)
      do j = 1, n
        b(j) = c
      end do
      go to 50
30    do j = 1, n
        b(j) = 0.
      end do
      b(m) = 1.0E+00
50    call ssolve ( b, a, b )

      c = 0.0E+00
      do j = 1, n
        c = c + abs ( b(j) )
        if ( b(j) .lt. 0.0E+00 ) then
          b(j) = -1.0E+00
        else
          b(j) = 1.0E+00
        end if
      end do

      call ssolve ( b, a, b )
c
c  I is the index of the largest entry of B.
c
      i = 1
      do j = 1, n
        if ( abs ( b(i) ) .lt. abs ( b(j) ) ) then
          i = j
        end if
      end do

      if ( m .eq. 0 ) go to 90
      if ( m .eq. i ) go to 100
      if ( d .ge. c ) go to 100
90    m = i
      d = c
      go to 30
100   c = c * a(3)
      c = max ( c, 1.0E+00 )
      scon = c
      return
      end
      function sdet ( iexp, a )

c*********************************************************************72
c
cc SDET computes the determinant of a symmetric matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, integer IEXP, a power of ten that is part of the
c    determinant.
c   
c    Input, real A(*), factorization information from SFACT.
c
c    Output, real SDET, the mantissa of the determinant.
c    Determinant = SDET * 10.0^IEXP.
c
      implicit none

      real a(*)
      real c
      real d
      real f
      real g
      integer h
      integer i
      integer iexp
      integer n
      real sdet

      iexp=0
      sdet=0.0
      d = a(1)

      if ( abs ( d ) .ne. 1233 )then
        write ( *, '(a)' ) ' '
        write(*,*)'sdet   - error!  matrix must be factored by SFACT!'
        return
      end if

      iexp = 0

      if ( d .lt. 0. ) then
        sdet = 0.0
        return
      end if

      d = 1.0E+00
      f = 2.0E+00**64
      g = 1.0E+00 /f
      h = 64
      n = a(2)
      i = 4
20    d = d*a(i)
      i = i + n
      n = n - 1
      if ( n .eq. 0 ) go to 50
30    if ( abs ( d ) .lt. f ) go to 40
      iexp = iexp + h
      d = d*g
      go to 30
40    if ( abs ( d ) .gt. g ) go to 20
      iexp = iexp - h
      d = d*f
      go to 40
50    if ( iexp .ne. 0 ) go to 60
      sdet = d
      return
60    if ( d .eq. 0. ) go to 90
      c = alog10 ( abs ( d ) ) + iexp * alog10 ( 2.0E+00 )
      iexp = c
      c = c - iexp
      if ( c .le. 0.0 ) go to 70
      c = c - 1
      iexp = iexp + 1
70    f = 10.0E+00**c
      if ( d .lt. 0. ) f = -f
      sdet = f
      return
90    iexp = 0
      sdet = 0.0
      return
      end
      subroutine sdiag ( e, v, lv, a, n, w ) 

c*********************************************************************72
c
cc SDIAG
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input, integer LV, the leading dimension of the array V.
c
      implicit none

      real a(*)
      real e(*)
      integer lv
      integer m
      integer n
      real v(*)
      real w(*)

      m = n + 1
      call shess ( w, w(m), a, n )
      call ssim ( v, lv, a )
      call tdg ( e, v, lv, w, w(m), n )

      return
      end
      subroutine sdiag2 ( e, v, lv, n, w )

c*********************************************************************72
c
cc SDIAG2 is the same as SDIAG, but the matrix is not stored in compressed format.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input, integer LV, the leading dimension of the array V.
c
      implicit none

      real e(*)
      integer i
      integer j
      integer k
      integer l
      integer lv
      integer m
      integer n
      integer o
      real v(*)
      real w(*)

      if ( n .eq. 1 ) then
        e(1) = v(1)
        v(1) = 1.0E+00
        return
      end if

      j = 0
      k = 1
      l = n
      m = n
      o = lv - n
20    do i = k, l
        v(i) = v(i+j)
      end do
      k = l + 1
      m = m - 1
      l = l + m
      o = o + 1
      j = j + o
      if ( m .gt. 0 ) go to 20
      m = n + 1
      if ( n .gt. 2 ) go to 40
      w(1) = v(1)
      w(2) = v(3)
      w(3) = v(2)
      v(1) = 2234
      v(2) = 2
      go to 50
40    call shess ( w, w(m), v, n )
50    call ssim ( v, lv, v )
      call tdg ( e, v, lv, w, w(m), n )
      return
      end
      subroutine sfact ( a, n, w )

c*********************************************************************72
c
cc SFACT computes the LU factorization of a symmetric matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input/output, real A((N*(N+1))/2); on input, information defining the
c    matrix, that is, a packed list of the row information, starting at the
c    diagonals.  On output, factorization information.
c
c    Input, integer N, the number of rows and columns in the matrix.
c
c    Workspace, real W(N).
c
      implicit none

      integer n

      real a(*)
      integer g
      integer h
      integer i
      integer j
      integer k
      integer l
      integer m
      real r
      real s
      real t
      real w(n)

      do i = 1, n
        w(i) = 0.
      end do

      i = -n
      k = 0
      r = 0.
      s = 0.
20    i = i + n - k
      k = k + 1
      j = k
      s = abs ( a(i+j) )
30    if ( j .eq. n ) go to 40
      j = j + 1
      t = abs ( a(i+j) )
      s = s + t
      w(j) = w(j) + t
      go to 30
40    s = s + w(k)
      r = max ( r, s )
      if ( k .lt. n ) go to 20
      j = 3 + (n+n*n)/2
50    a(j) = a(j-3)
      j = j - 1
      if ( j .gt. 3 ) go to 50
      a(1) = 1233
      a(2) = n
      a(3) = r
      h = n
      k = 4
60    if ( h .eq. 1 ) go to 90
      s = a(k)
      k = k + h
      g = k
      h = h - 1
      m = h
      if ( s .eq. 0. ) go to 100
      j = 0
70    j = j - m
      m = m - 1
      l = g + m
      t = a(g+j)/s
      do i = g, l
        a(i) = a(i) - t*a(i+j)
      end do
      g = l + 1
      if ( m .gt. 0 ) go to 70
      go to 60
90    if ( a(k) .ne. 0. ) return
      a(1) = -1233
      return
100   a(1) = -1233
      go to 60
      end
      subroutine sft ( s, a, b, c, d, e2, e1, e0, f2, f1 )

c*********************************************************************72
c 
cc SFT
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real a
      real b
      real c
      real d
      real e0
      real e1
      real e2
      real f1
      real f2
      real g0
      real g1
      real g2
      real h1
      real h2
      real s
      real w
      real x
      real y
      real z

      g0 = abs ( e0 )
      g1 = abs ( e1 )
      g2 = abs ( e2 )
      h1 = abs ( f1 )
      h2 = abs ( f2 )

      w = (a*g2)*(a*h2) + (c*g1)*(c*h2)
      x = (b*g1)*(b*h1) + (d*g0)*(d*h1)
      y = (b*g1)*(c*h1)
      z = (b*g1)*(c*h2)
      call eig3 ( s, s, x, w, y, z )

      return
      end
      subroutine shess ( d, u, a, n )

c*********************************************************************72
c
cc SHESS
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real a(*)
      real d(*)
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      integer o
      integer p
      integer q
      real r
      real s
      real t
      real u(*)
      real v
      real w

      d(1) = a(1)
      if ( n .gt. 2 ) go to 20
      a(1) = 2234
      if ( n .eq. 1 ) go to 10
      u(1) = a(2)
      d(2) = a(3)
      a(4) = d(2)
      a(5) = u(1)
10    a(2) = n
      a(3) = d(1)
      return
20    m = n - 1
      q = 1
      k = 2
      l = n
30    if ( q .ge. m ) go to 150
      w = a(k)
      p = k + 1
      do i = p, l
        if ( a(i) .ne. 0. ) go to 50
      end do
      u(q) = w
      k = l + 2
      i = q
      q = q + 1
      d(q) = a(l+1)
      l = l + n - i
      go to 30
50    t = abs ( w )
      if ( t .ne. 0. ) v = 1.0E+00 /t
      r = 1.0E+00
      do 70 j = i, l
        s = abs ( a(j) )
        if ( s .le. t ) go to 60
        v = 1.0E+00 /s
        r = 1.0E+00 + r*(t*v)**2
        t = s
        go to 70
60      r = r + (s*v)**2
70    continue

      s = t * sqrt ( r )
      v = 1.0E+00 / sqrt ( s * ( s + abs ( w ) ) )
      if ( w .lt. 0. ) s = -s
      u(q) = -s
      a(k) = v*(w+s)
      q = q + 1
      d(q) = a(k)
      u(q) = 0.
      o = q - k
      do i = p, l
        a(i) = v*a(i)
        j = o + i
        d(j) = a(i)
        u(j) = 0.
      end do
      i = q
      l = l + 1
      o = l - q
      s = 0.

90    v = d(i)
      t = v*a(i+o)
      if ( i .ge. n ) go to 110
      p = i + 1

      do j = p, n
        w = a(j+o)
        t = t + d(j)*w
        u(j) = u(j) + v*w
      end do

      u(i) = u(i) + t
      s = s + u(i)*v
      o = o + n - i
      i = p
      go to 90

110   u(n) = u(n) + t
      s = .5*(s+u(n)*v)
      do i = q, n
        u(i) = s*d(i) - u(i)
      end do
      o = l - q

      do j = q, n
        t = u(j)
        s = d(j)
        do i = j, n
          p = i + o
          a(p) = a(p) + d(i)*t + u(i)*s
        end do
        o = o + n - j
      end do

      d(q) = a(l)
      k = l + 1
      l = l + n - q
      go to 30
150   u(q) = a(l)
      d(q+1) = a(l+1)
      k = n - 1
      m = 1
      l = k
      o = 1
160   do i = m, l
        a(i) = a(i+o)
      end do
      m = l + 1
      k = k - 1
      l = l + k
      o = o + 1
      if ( k .gt. 1 ) go to 160
      j = (n*(n-1))/2 + 1
      i = j
180   a(i) = a(i-2)
      i = i - 1
      if ( i .gt. 2 ) go to 180
      a(1) = 2234
      a(2) = n
      return
      end
      subroutine sim ( p, lp, a )

c*********************************************************************72
c
cc SIM
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      integer lp

      real a(*)
      integer h
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      real p(lp,*)
      real t

      t = a(1)
      if ( t .eq. 2230 ) go to 10
      if ( t .eq. 2231 ) go to 10
      write ( *, '(a)' ) ' '
      write(*,*) 'error: must reduce matrix to hessenberg form using'
      write(*,*) 'routine hess before using routine sim to'
      write(*,*) 'compute the similarity transformation'
      stop
10    continue

      n = a(2)
      l = n - 2
      m = n + 1
      h = 3

      do j = 1, l
        k = j + 1
        do i = k, n
          p(i,j) = a(i+h)
        end do
        h = h + m
      end do

      call hsr1 ( p, lp, n )

      if ( t .eq. 2230 ) then
        return
      end if

      k = 1 + n + n*n

      do j = 1, n
        do i = 1, n
          p(i,j) = a(i+k)*p(i,j)
        end do
      end do

      return
      end
      subroutine sing ( q, lq, iq, s, p, lp, ip, a, la, m, n, w )

c*********************************************************************72
c
cc SING computes the singular value decomposition of a general matrix.
c
c  Discussion:
c     
c    A = Q times diagonal matrix times P transpose     
c                                                  
c              either p or q can be identified with a but not 
c              both. when either p or q are identified with a,
c              the householder vectors in a are destroyed. if 
c              iq = 2, q must have m columns even though the  
c              output array has just m-l columns. similarly if
c              ip = 2, p must have n columns even though the  
c              output array has just n-l columns.             
c     
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c                                  
c         input:                                             
c                                                             
c              s     --array with at least n elements         
c                                                             
c              lq    --leading (row) dimension of array q     
c                                                             
c              iq    --an integer which indicates which col-  
c                      umns of q to compute (= 0 means none,  
c                      = 1 means first l, = 2 means last m-l, 
c                      = 3 means all m where l = min ( m, n ))    
c                                                             
c              lp    --leading (row) dimension of array p     
c                                                             
c              ip    --an integer (like iq) which indicates   
c                      which columns of p to compute          
c                                                             
c              a     --array containing coefficient matrix    
c                                                             
c              la    --leading (row) dimension of array a     
c                                                             
c              m     --row dimension of matrix stored in a    
c                                                             
c              n     --column dimension of matrix stored in a 
c                                                             
c              w     --work array(length at least max(m,3l-1))
c                                                             
c         output:                                             
c                                                            
c              q     --q factor in the singular value decomp. 
c                                                             
c              s     --singular values in decreasing order    
c                                                             
c              p     --p factor in the singular value decomp. 
c                                                             
c              a     --the householder vectors used in the    
c                      reduction process                      
c
      implicit none
                                                                                                               
      integer la
      integer lp
      integer lq

      real a(la,*)
      integer i
      integer ip
      integer iq
      integer iu
      integer jl
      integer jr
      integer l
      integer m
      integer n
      real p(lp,*)
      real q(lq,*)
      real s(*)
      real w(*)

      if ( iq .ge. 0 ) go to 20
10    continue

      write ( *, '(a)' ) ' '
      write(*,*) 'error: input parameter iq for routine sing'
      write(*,*) 'either less than 0 or greater than 3'
      stop
20    continue

      if ( iq .gt. 3 ) then
        go to 10
      end if

      jl = 0
      if ( iq .eq. 0 ) go to 30
      if ( iq .eq. 2 ) go to 30
      jl = 1
30    if ( ip .ge. 0 ) go to 50
40    continue
      write ( *, '(a)' ) ' '
      write(*,*) 'error: input parameter ip for routine sing'
      write(*,*) 'either less than 0 or greater than 3'
      stop
50    if ( ip .gt. 3 ) go to 40
      jr = 0
      if ( ip .eq. 0 ) go to 60
      if ( ip .eq. 2 ) go to 60
      jr = 1
60    call bidag2 ( s, w, q, lq, iq, p, lp, ip, a, la, m, n )
      l = min ( m, n )
      if ( l .gt. 1 ) go to 80
      if ( s(1) .ge. 0. ) return
      s(1) = -s(1)
      if ( jl .eq. 0. ) return

      do i = 1, m
        q(i,1) = -q(i,1)
      end do

      return
80    iu = 0
      if ( m .ge. n ) go to 90
      iu = 1
90    call singb(s,l,w,iu,q,lq,m,jl,p,lp,n,jr,w(l),w(l+l))
      return
      end
      subroutine singb ( d, n, u, iu, q, lq, mq, iq, p, lp, mp, ip, 
     &  e, f )

c*********************************************************************72
c
cc SINGB computes the singular value decomposition of a bidiagonal matrix
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      integer lp
      integer lq

      real b
      real c
      real d(*)
      real e(*)
      real f(*)
      integer g
      integer h
      integer i
      integer id
      integer ip
      integer iq
      integer iu
      integer j
      integer jl
      integer jr
      integer j0
      integer k
      integer k2
      integer l
      integer ll
      integer l1
      integer m
      integer mp
      integer mq
      integer n
      integer ns
      real p(lp,*)
      real q(lq,*)
      real r
      real s
      real t
      real t0
      real t1
      real t2
      real t3
      real u(*)
      real v
      real w
      real x
      real y
      real z

      if ( n .gt. 1 ) go to 10
      if ( iq .eq. 1 ) q(1,1) = 1.0E+00
      if ( ip .eq. 1 ) p(1,1) = 1.0E+00

      if ( d(1) .ge. 0. ) then
        return
      end if

      d(1) = -d(1)
      if ( iq .eq. 1 ) q(1,1) = -1.0E+00
      return
10    jl = iq
      if ( jl .eq. 0 ) jl = 3
      if ( jl .ne. 3 ) jl = 1
      jr = ip
      if ( jr .eq. 0 ) jr = 3
      if ( jr .ne. 3 ) jr = 2
      if ( iu .eq. 0 ) go to 20
      i = jl
      jl = jr
      jr = i
20    j = 0
      l = n - 1
      k2 = n - 2

      do i = 1, l
        e(i) = 1.0E+00
        f(i) = 1.0E+00
        if ( u(i) .eq. 0. ) j = i
        if ( d(i) .eq. 0. ) j = i
      end do

      e(n) = 1.0E+00
      f(n) = 1.0E+00
      b = 65536.0E+00**(-3)
      t = 1.0E+00
40    t = .5*t
      s = 1.0E+00 + t
      if ( s .gt. 1.0E+00 ) go to 40
      t0 = 1.0E+00 /(t+t)
      t2 = (t+t)**2
      ns = 50*n
      l1 = 0
      k = n
      ll = 0
      go to 70
50    j = 0
      do i = 1, l
        if ( u(i) .eq. 0. ) j = i
        if ( d(i) .eq. 0. ) j = i
      end do
70    if ( j .eq. 0 ) go to 140
      if ( u(j) .eq. 0. ) go to 140
      i = j
      v = u(j)
      u(j) = 0.
      s = -v * abs ( e ( j ) )
80    h = i
      i = i + 1
      r = abs ( e(i) ) * d(i)
      call fgv(x,y,t,r,s,e(j),e(i))
      if ( t .eq. 1.0E+00 ) go to 90
      d(i) = d(i) - y*v
      call sng0(q,lq,mq,p,lp,mp,jl,j,i,x,y)
      if ( i .eq. k ) go to 100
      v = x*u(i)
      s = -v * abs ( e(j) )
      go to 80
90    d(i) = d(i)*y - v
      call sng1(q,lq,mq,p,lp,mp,jl,j,i,x,y)
      if ( i .eq. k ) go to 100
      v = u(i)
      u(i) = v*y
      s = -v * abs ( e(j) )
      go to 80
100   if ( j .eq. 1 ) go to 130
      i = j - 1
      s = u(i)
      u(i) = 0.
110   h = i
      i = i - 1
      r = d(h)
      call fgv(x,y,t,r,s,f(h),f(j))
      if ( t .eq. 1.0E+00 ) go to 120
      d(h) = r + x*s
      call sng0(q,lq,mq,p,lp,mp,jr,h,j,x,y)
      if ( h .eq. 1 ) go to 130
      s = -y*u(i)
      go to 110
120   d(h) = x*r + s
      call sng1(q,lq,mq,p,lp,mp,jr,h,j,x,y)
      if ( h .eq. 1 ) go to 130
      s = -u(i)
      u(i) = x*u(i)
      go to 110
130   call scl(d,u,n,q,lq,mq,p,lp,mp,e,f,b,1,k,jl,jr)
140   j = j + 1
      if ( j .eq. k ) go to 320

      s = 0.
      t = 0.
      do i = j, l
        x = abs ( ( d(i) * e(i) ) * ( d(i) * f(i) ) )
        y = abs ( ( u(i) * e(i) ) * ( u(i) * f(i+1) ) )
        s = max ( s, x, y )
        t = t + x + y
      end do

      x = abs ( ( e(k) * d(k) ) * ( f(k) * d(k) ) )
      s = max ( s, x )
      t3 = s*t2
      t = t + x
      if ( t .eq. 0. ) t = 1.0E+00
      t1 = t0/t
      go to 280
160   ll = ll + 1
      if ( ll .gt. ns ) go to 530
      if ( l .gt. j ) go to 170
      s = 0.
      t = 0.
      go to 180
170   s = u(k2)
      t = e(k2)
180   call sft(c,d(k),d(l),u(l),s,e(k),e(l),t,f(k),f(l))
      v = abs ( e(j) )
      w = abs ( f(j) )
      t = d(j)*w
      r = t*(d(j)*v) - c
      s = t*(u(j)*v)
      id = 1
      j0 = j
      i = j
      h = j - 1
190   g = h
      h = i
      i = i + 1
      z = abs ( f(i) )
      call fgv(x,y,t,r,s,f(h),f(i))
      v = d(h)
      if ( t .eq. 1 ) go to 210
      if ( h .eq. j ) go to 200
      t = u(g) + x*s
      u(g) = t
      if ( abs ( ( t *e(g))*(t*f(h)) ) .gt. t3 ) go to 200
      j = h
      id = 1
200   r = v + x*u(h)
      s = x*d(i)
      u(h) = u(h) - y*v
      call sng0(q,lq,mq,p,lp,mp,jr,h,i,x,y)
      go to 230
210   if ( h .eq. j ) go to 220
      t = x*u(g) + s
      u(g) = t
      if ( abs ( (t*e(g))*(t*f(h)) ) .gt. t3 ) go to 220
      j = h
      id = 1
220   r = x*v + u(h)
      s = d(i)
      u(h) = y*u(h) - v
      d(i) = y*s
      call sng1(q,lq,mq,p,lp,mp,jr,h,i,x,y)
230   call fgv(x,y,t,r,s,e(h),e(i))
      if ( t .eq. 1.0E+00 ) go to 250
      t = r + x*s
      d(h) = t
      if ( abs ( (t*e(h))*(t*f(h)) ) .gt. t3 ) go to 240
      id = 0
      j = h
240   r = u(h) + x*d(i)
      d(i) = d(i) - y*u(h)
      u(h) = r
      call sng0(q,lq,mq,p,lp,mp,jl,h,i,x,y)
      if ( i .eq. k ) go to 270
      s = x*u(i)
      go to 190
250   t = s + x*r
      d(h) = t
      if ( abs ( (t*e(h))*(t*f(h)) ) .gt. t3 ) go to 260
      id = 0
      j = h
260   r = d(i) + x*u(h)
      d(i) = y*d(i) - u(h)
      u(h) = r
      call sng1(q,lq,mq,p,lp,mp,jl,h,i,x,y)
      if ( i .eq. k ) go to 270
      s = u(i)
      u(i) = s*y
      go to 190
270   call scl(d,u,n,q,lq,mq,p,lp,mp,e,f,b,j0,k,jl,jr)
      if ( id .eq. 0 ) go to 70
280   w = e(l)
      x = e(k)
      y = f(l)
      z = f(k)
      r = abs ( (x*d(k))*(z*d(k)) )
      s = abs ( (w*d(l))*(y*d(l)) )
      t = abs ( (x*u(l))*(y*u(l)) )
      if ( (s*t1)*(t*t1) .gt. 1.0E+00 ) go to 160
      l1 = l1 + 1
      if ( l1 .gt. 40 ) go to 290
      if ( t .eq. 0. ) go to 290
      r = r + t
      if ( (s/r)*(t/r) .gt. t2 ) go to 160
290   l1 = 0
      if ( s .gt. r ) go to 310
      r = -d(k) * abs ( x )
      s = u(l) * abs ( w )
      call fgv(w,y,t,r,s,e(l),e(k))
      x = e(k)
      if ( t .eq. 1.0E+00 ) go to 300
      d(k) = d(k) - y*u(l)
      call sng0(q,lq,mq,p,lp,mp,jl,l,k,w,y)
      go to 310
300   d(l) = d(l)*w
      d(k) = y*d(k) - u(l)
      call sng1(q,lq,mq,p,lp,mp,jl,l,k,w,y)
310   t = sign ( sqrt ( abs ( x ) ), x ) * d(k) 
     &  * sign ( sqrt ( abs ( z ) ), z )
      if ( t .lt. 0. ) e(k) = -e(k)
      d(k) = abs ( t )
      k = l
      l = k2
      k2 = k2 - 1
      if ( k .gt. j ) go to 280
320   x = e(k)
      z = f(k)
      t = sign ( sqrt ( abs ( x ) ), x ) 
     &  * d(k) * sign ( sqrt ( abs ( z ) ), z )
      if ( t .lt. 0. ) e(k) = -e(k)
      d(k) = abs ( t )
      k = l
      l = k2
      k2 = k2 - 1
      if ( k .gt. 1 ) go to 50
      if ( k .eq. 0 ) go to 330
      go to 320
330   go to (340,360,380),jl
340   continue

      do j = 1, n
        t = e(j)
        t = sign ( sqrt ( abs ( t ) ), t )
        do i = 1, mq
          q(i,j) = q(i,j)*t
        end do
      end do

      go to 380

360   do j = 1, n
        t = e(j)
        t = sign ( sqrt ( abs ( t ) ), t )
        do i = 1, mp
          p(i,j) = p(i,j)*t
        end do
      end do

380   go to (390,410,430),jr

390   do j = 1, n
        t = f(j)
        t = sign ( sqrt ( abs ( t ) ), t )
        do i = 1, mq
          q(i,j) = q(i,j)*t
        end do
      end do

      go to 430

410   do j = 1, n
        t = f(j)
        t = sign ( sqrt ( abs ( t ) ), t )
        do i = 1, mp
          p(i,j) = p(i,j)*t
        end do
      end do

430   call sort2(d,e,f,n)

      do i = 1, n
        j = e(i)
        f(j) = i
      end do

      m = n + 1

      do 520 j = 1, n

        l = m - j
        k = e(l)
        i = f(j)
        e(i) = k
        f(k) = i
        t = d(j)
        d(j) = d(k)
        d(k) = t
        if ( jl .eq. 1 ) go to 450
        if ( jr .ne. 1 ) go to 480
450     s = 0.
        do i = 1, mq
          t = q(i,k)
          q(i,k) = q(i,j)
          q(i,j) = t
          s = s + t*t
        end do
        s = 1.0E+00 / sqrt ( s )

        do i = 1, mq
          q(i,j) = s*q(i,j)
        end do

480     if ( jr .eq. 2 ) go to 490
        if ( jl .ne. 2 ) go to 520
490     s = 0.
        do i = 1, mp
          t = p(i,k)
          p(i,k) = p(i,j)
          p(i,j) = t
          s = s + t*t
        end do
        s = 1.0E+00 / sqrt ( s )
        do i = 1, mp
          p(i,j) = s*p(i,j)
        end do

520   continue

      e(1) = n
      return
530   k = n - k + 1
      write(*,*) 'since the stopping criterion not satisfied'
      write(*,*) 'after',ns,'iterations, we stop while computing'
      write(*,*) 'eigenvalue number',k
      e(1) = k
      return
      end
      subroutine slice ( e, k, y, z, l, d, u, n, w )

c*********************************************************************72
c
cc SLICE computes the eigenvalues in a given interval for a tridiagonal matrix.
c
c  Discussion:
c
c    The cross-diagonal products of the triadiagonal matrix must be nonnegative,
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real a
      real b
      real c
      real d(*)
      real e(*)
      integer f
      real fl
      real fr
      integer g
      real gl
      real gr
      integer h
      integer i
      integer j
      integer k
      integer ki
      integer k0
      integer k2
      real l(*)
      integer m
      integer n
      integer o
      integer p
      integer q
      real r
      real s
      real t
      real u(*)
      real w(*)
      real y
      real z

      a = y
      b = z

      if ( a .gt. b ) then
        a = z
        b = y
      end if

      if ( n .le. 1 ) go to 220

      m = n - 1
      fr = d(1)
      fl = d(1)
      s = 0.

      do i = 1, m
        w(i) = l(i)*u(i)
        if ( w(i) .lt. 0. ) go to 230
        s = s + abs ( u(i) )
        t = d(i) + s
        fl = max ( fl, t )
        t = d(i) - s
        fr = min ( fr, t )
        s = abs ( l(i) )
      end do

      t = d(n) + s
      fl = max ( fl, t )
      t = d(n) - s
      fr = min ( fr, t )
      t = 1.0E+00
30    t = .5*t
      s = 1.0E+00 + t
      if ( s .gt. 1.0E+00 ) go to 30
      t = 8.*n*t* max ( abs ( fl ), abs ( fr ) )
      s = t + t
      call cp(a,fl,gl,i,d,w,n)
      call cp(b,fr,gr,j,d,w,n)
      k = j - i
      k0 = i
      if ( k .le. 0 ) go to 240
      do i = 1, k
        w(i+m) = 0.
      end do
      k2 = k + k
      f = n + k
      g = f + k2
      h = g + k2
      p = f - 1
      q = p + k
      call sto(a,fl,gl,0,k,w(n),w(f),w(g),w(h))
      call sto(b,fr,gr,k,k,w(n),w(f),w(g),w(h))

      do 210 i = 1, k

        r = w(m+i)
        if ( r .eq. 0. ) go to 110
        if ( r .eq. 1.0E+00 ) go to 50
        if ( r .eq. 2. ) go to 80
        a = w(p+i)
        b = w(q+i)
        if ( abs ( a - b ) .le. s ) go to 200
        go to 190
50      a = w(p+i)
        do j = i, k
          if ( w(m+j) .ge. 2. ) go to 70
        end do
70      b = w(q+j)
        go to 170
80      b = w(q+i)
        o = m + i + 1
        do j = 1, i
          r = w(o-j)
          if ( r .eq. 1.0E+00 ) go to 100
          if ( r .eq. 3. ) go to 100
        end do
100     a = w(o-j+k)
        go to 170
110     o = m + i + 1
        do j = 1, i
          r = w(o-j)
          if ( r .eq. 1.0E+00 ) go to 130
          if ( r .eq. 3. ) go to 130
        end do
130     a = w(o-j+k)
        do j = i, k
          if ( w(m+j) .ge. 2. ) go to 150
        end do
150     b = w(q+j)
        go to 170
160     if ( w(m+i) .eq. 3. ) go to 190
170     if ( abs ( a - b ) .le. s ) go to 200
        c = .5*(a+b)
        call cp(c,fl,gl,j,d,w,n)
        j = j - k0
        call sto(c,fl,gl,j,k,w(n),w(f),w(g),w(h))
        if ( j .ge. i ) go to 180
        a = c
        go to 160
180     b = c
        go to 160
190     j = h + i + i - 2
        call db(gl,w(j))
        call db(gr,w(j+k+k))
        ki = i + k0
        call stm2
     &    ( e(i),a,b,w(p+i+k2),w(q+i+k2),gl,gr,ki-1,ki,ki,t,d,w,n)
        go to 210
200     e(i) = .5*(a+b)
210   continue

      return
220   if ( d(1) .lt. a ) go to 240
      if ( d(1) .gt. b ) go to 240
      e(1) = d(1)
      k = 1
      return
230   continue
      write ( *, '(a)' ) ' '
      write(*,*) 'error: routine slice can only be applied',
     &  'when the cross-diagonal products are nonnegative'
      stop
240   k = 0
      return
      end
      subroutine smult ( y, x, a, n )

c*********************************************************************72
c
cc SMULT multiplies a symmetric matrix times a vector.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, real Y(N), the product A*x.
c
c    Input, real X(N), the vector to be multiplied.
c
c    Input, real A((N*(N+1))/2), the symmetric array, in packed storage,
c    containing the elements of successive rows starting at the
c    diagonal element.
c
c    Input, integer N, the number of rows and columns in A.
c
      implicit none

      integer n

      real a(*)
      integer i
      integer j
      integer k
      integer l
      real r
      real s
      real t
      real x(n)
      real y(n)

      do i = 1, n
        y(i) = 0.0
      end do

      k = 1
      l = 0

      do j = 2, n
        t = x(k)
        s = a(k+l) * t
        do i = j, n
          r = a(i+l)
          s = s + r * x(i)
          y(i) = y(i) + r * t
        end do
        y(k) = y(k) + s
        l = l + n - k
        k = j
      end do

      y(n) = y(n) + a(k+l) * x(n)

      return
      end
      subroutine sng0 ( q, lq, m, p, lp, n, l, j, k, x, y )

c*********************************************************************72
c
cc SNG0
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      integer lp
      integer lq

      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      real p(lp,*)
      real q(lq,*)
      real s
      real t
      real x
      real y

      if ( l .eq. 1 ) then

        do i = 1, m
          t = q(i,j)
          s = q(i,k)
          q(i,j) = t + x * s
          q(i,k) = s - y * t
        end do

      else if ( l .eq. 2 ) then

        do i = 1, n
          t = p(i,j)
          s = p(i,k)
          p(i,j) = t + x*s
          p(i,k) = s - y*t
        end do

      end if

      return
      end
      subroutine sng1 ( q, lq, m, p, lp, n, l, j, k, x, y )

c*********************************************************************72
c
cc SNG1
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      integer lp
      integer lq

      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      real p(lp,*)
      real q(lq,*)
      real s
      real t
      real x
      real y

      if ( l .eq. 1 ) then

        do i = 1, m
          t = q(i,j)
          s = q(i,k)
          q(i,j) = x*t + s
          q(i,k) = y*s - t
        end do

      else if ( l .eq. 2 ) then

        do i = 1, n
          t = p(i,j)
          s = p(i,k)
          p(i,j) = x*t + s
          p(i,k) = y*s - t
        end do

      end if

      return
      end
      subroutine solve ( x, a, b )

c*********************************************************************72
c
cc SOLVE solves a general factored linear system.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, real X(N), the solution.
c
c    Input, real A(*), factorization information from FACT.
c
c    Input, real B(N), the right hand side.
c
      implicit none

      real a(*)
      real b(*)
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      real t
      real x(*)

      t = a(1)

      if ( abs ( t ) .ne. 1230 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor before solving'
        stop
      end if

      n = a(2)
      m = n + 1
      j = 4 - m
      if ( t .lt. 0. ) go to 80
      do i = 1, n
        x(i) = b(i)
      end do
      k = 1
30    j = j + m
      if ( a(j+k) .eq. 0. ) go to 80
      if ( k .eq. n ) go to 50
      l = a(j)
      t = x(l)
      x(l) = x(k)
      x(k) = t
      k = k + 1
      if ( t .eq. 0. ) go to 30
      do i = k, n
        x(i) = x(i) - t*a(i+j)
      end do
      go to 30
50    t = x(k)/a(j+k)
60    x(k) = t
      if ( k .eq. 1 ) return
      k = k - 1
      do i = 1, k
        x(i) = x(i) - t*a(i+j)
      end do
      j = j - m
      go to 50
80    k = 0
90    k = k + 1
      j = j + m
      if ( a(j+k) .ne. 0. ) go to 90

      do i = 1, n
        x(i) = 0.0
      end do

      t = 1.0E+00
      go to 60
      end
      subroutine sort ( x, y, n )

c*********************************************************************72
c
cc SORT sorts a real array.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input/output, real X(N), the array to be sorted.
c
c    Workspace, real Y(N).
c
c    Input, integer N, the dimension of X.
c
      implicit none

      integer n

      integer i
      integer j
      integer k
      integer l
      integer m
      real s
      real t
      real x(n)
      real y(n)

      i = 1

10    k = i
20    j = i
      i = i + 1
      if ( j .eq. n ) then
        go to 30
      end if
      if ( x(i) .ge. x(j) ) go to 20
      y(k) = i
      go to 10

30    if ( k .eq. 1 ) then
        return
      end if

      y(k) = n + 1
40    m = 1
      l = 1
50    i = l
      if ( i .gt. n ) go to 120
      s = x(i)
      j = y(i)
      k = j
      if ( j .gt. n ) go to 100
      t = x(j)
      l = y(j)
      x(i) = l
60    if ( s .gt. t ) go to 70
      y(m) = s
      m = m + 1
      i = i + 1
      if ( i .eq. k ) go to 80
      s = x(i)
      go to 60
70    y(m)= t
      m = m + 1
      j = j + 1
      if ( j .eq. l ) go to 110
      t = x(j)
      go to 60
80    y(m) = t
      k = m + l - j
      i = j - m
90    m = m + 1
      if ( m .eq. k ) go to 50
      y(m) = x(m+i)
      go to 90
100   x(i) = j
      l = j
110   y(m) = s
      k = m + k - i
      i = i - m
      go to 90
120   i = 1
130   k = i
      j = x(i)
140   x(i) = y(i)
      i = i + 1
      if ( i .lt. j ) go to 140
      y(k) = i
      if ( i .le. n ) go to 130
      if ( k .eq. 1 ) return
      go to 40
      end
      subroutine sort2 ( x, y, w, n )

c*********************************************************************72
c
cc SORT2 index-sorts a real array.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input, real X(N), the array to be sorted.
c
c    Output, integer Y(N), indexes the entries of X in increasing order.
c
c    Workspace, integer W(N).
c
c    Input, integer N, the dimension of X.
c
      implicit none

      integer n

      integer i
      integer j
      integer k
      integer l
      integer m
      integer p
      integer q
      real s
      real t
      integer w(n)
      real x(n)
      integer y(n)

      i = 1
10    k = i
20    j = i
      y(i) = i
      i = i + 1
      if ( j .eq. n ) go to 30
      if ( x(i) .ge. x(j) ) go to 20
      w(k) = i
      go to 10
30    if ( k .eq. 1 ) return
      w(k) = n + 1
40    m = 1
      l = 1
50    i = l
      if ( i .gt. n ) go to 120
      p = y(i)
      s = x(p)
      j = w(i)
      k = j
      if ( j .gt. n ) go to 100
      q = y(j)
      t = x(q)
      l = w(j)
      y(i) = l
60    if ( s .gt. t ) go to 70
      w(m) = p
      m = m + 1
      i = i + 1
      if ( i .eq. k ) go to 80
      p = y(i)
      s = x(p)
      go to 60
70    w(m)= q
      m = m + 1
      j = j + 1
      if ( j .eq. l ) go to 110
      q = y(j)
      t = x(q)
      go to 60
80    w(m) = q
      k = m + l - j
      i = j - m
90    m = m + 1
      if ( m .eq. k ) go to 50
      w(m) = y(m+i)
      go to 90
100   y(i) = j
      l = j
110   w(m) = p
      k = m + k - i
      i = i - m
      go to 90
120   i = 1
130   k = i
      j = y(i)
140   y(i) = w(i)
      i = i + 1
      if ( i .lt. j ) go to 140
      w(k) = i
      if ( i .le. n ) go to 130
      if ( k .gt. 1 ) go to 40
      return
      end
      function sq ( a )

c*********************************************************************72
c
cc SQ
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      complex a
      complex b
      real c(2)
      real sq

      equivalence (b,c)

      b = a
      sq = c(1)**2 + c(2)**2

      return
      end
      function sqr ( c, t )

c*********************************************************************72
c
cc SQR
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      complex c
      complex d
      real e(2)
      real sqr
      real t

      equivalence (d,e)

      d = c
      sqr = (e(1)*t)**2 + (e(2)*t)**2

      return
      end
      subroutine ssim ( p, lp, a )

c*********************************************************************72
c
cc SSIM computes a similarity transform for reduction to tridiagonal form.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      integer lp

      real a(*)
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      integer o
      real p(lp,*)
      real t

      t = a(1)

      if ( t .ne. 2234 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must reduce matrix to tridiagonal form using'
        write(*,*) 'routine shess before using routine ssim to'
        write(*,*) 'compute the similarity transformation'
        stop
      end if

      n = a(2)

      if ( n .lt. 2 ) then
        call hsr2 ( p, lp, n )
        return
      end if

      k = 1
      l = n - 2
      o = (l*n-l)/2

      do m = 1, n - 2
        j = n - m
        do i = j, n
          p(i,j) = a(i+o)
        end do
        k = k + 1
        o = o - k
      end do

      do i = 1, n
        p(i,1) = 0.
      end do

      call hsr2 ( p, lp, n )

      return
      end
      subroutine ssolve ( x, a, b )

c*********************************************************************72
c
cc SSOLVE
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, real X(N), the solution.
c
c    Input, real A(*), the factorization output from SFACT.
c
c    Input, real B(N), the right hand side.
c
      implicit none

      real a(*)
      real b(*)
      integer i
      integer j
      integer k
      integer l
      integer n
      real t
      real x(*)

      t = a(1)

      if ( abs ( t ) .ne. 1233 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor with sfact before solving'
        stop
      end if

      n = a(2)
      l = 3
      k = 1
      if ( t .lt. 0. ) go to 80
      do i = 1, n
        x(i) = b(i)
      end do

30    if ( k .eq. n ) go to 50
      t = x(k)/a(k+l)
      j = l
      l = l + n - k
      k = k + 1
      if ( t .eq. 0. ) go to 30
      do i = k, n
        x(i) = x(i) - t*a(i+j)
      end do
      go to 30
50    x(n) = x(n)/a(k+l)
60    if ( k .eq. 1 ) return
      j = k
      k = k - 1
      l = l + k - n
      t = x(k)
      do i = j, n
        t = t - x(i)*a(i+l)
      end do
      x(k) = t/a(k+l)
      go to 60
80    if ( a(k+l) .eq. 0. ) go to 90
      l = l + n - k
      k = k + 1
      go to 80
90    continue

      do i = 1, n
        x(i) = 0.0
      end do

      x(k) = 1.0E+00
      go to 60
      end
      subroutine stm ( e, y, z, k, t, d, p, n )

c*********************************************************************72
c
cc STM
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real a
      real b
      real c
      real d(*)
      real d2
      real d23
      real d3
      real d34
      real d4
      real d42
      real e
      real f(4)
      real fl
      real fr
      real g(4)
      real gl
      real gr
      real h(4)
      integer i
      integer j
      integer k
      real l
      integer m
      integer n
      real p(*)
      real p2
      real p3
      real p4
      real q
      real r
      real s
      real t
      real u
      real v
      real w
      real x(4)
      real y
      real z

      l = y
      r = z
      if ( l .lt. r ) go to 10
      l = z
      r = y
10    if ( n .eq. 1 ) go to 210
      s = t + t
      u = 1.0E+00
20    u = .5*u
      a = 1.0E+00 + u
      if ( a .gt. 1.0E+00 ) go to 20
      u = 5.*u
      v = .5*u
      call cp(l,fl,gl,i,d,p,n)
      call cp(r,fr,gr,j,d,p,n)
      if ( i .ge. k ) go to 190
      if ( j .lt. k ) go to 200
30    e = r - l
      if ( e .le. s ) go to 180
      if ( e .le. u*( abs ( l ) + abs ( r ) ) ) go to 180
      if ( j .eq. i+1 ) go to 70
40    a = .5*(l+r)
      call cp(a,b,c,m,d,p,n)
      if ( k .le. m ) go to 50
      l = a
      i = m
      fl = b
      gl = c
      go to 30
50    r = a
      j = m
      fr = b
      gr = c
      go to 30
60    e = r - l
      if ( e .le. s ) go to 180
      if ( e .le. u*( abs ( l ) + abs ( r ) ) ) go to 180
70    x(1) = l
      f(1) = fl
      g(1) = gl
      x(2) = r
      f(2) = fr
      g(2) = gr
      call eql(x,f,g,h,2)
      if ( h(1) .eq. h(2) ) go to 160
      a = x(1) - h(1)*(x(1)-x(2))/(h(1)-h(2))
      q = a
      w = max ( t, v * ( abs ( l ) + abs ( r ) ) )
      if ( abs ( a - l ) .lt. w ) a = l + w
      if ( abs ( a - r ) .lt. w ) a = r - w
      call cp(a,b,c,j,d,p,n)
      if ( i .ge. j ) go to 80
      r = a
      fr = b
      gr = c
      go to 90
80    l = a
      fl = b
      gl = c
90    x(3) = a
      f(3) = b
      g(3) = c
      w = r - l
      if ( w .le. s ) go to 220
      if ( w .le. u*( abs ( l ) + abs ( r ) ) ) go to 220
      call eql(x,f,g,h,3)
      call inp(a,x(1),x(2),x(3),h(1),h(2),h(3),l,r)
      b = l
      if ( abs ( a - l ) .gt. abs ( a - r ) ) b = r
100   q = a
      w = max ( t, v * ( abs ( l ) + abs ( r ) ) )
      if ( abs ( a - l ) .lt. w ) go to 110
      if ( abs ( a - r ) .gt. w ) go to 130
110   if ( a+a .gt. l+r ) go to 120
      a = l + w
      go to 130
120   a = r - w
130   if ( a .le. l ) go to 160
      if ( a .ge. r ) go to 160
      e = .5*e
      if ( e .lt. abs ( b - a ) ) go to 160
      call cp(a,b,c,j,d,p,n)
      if ( i .ge. j ) go to 140
      r = a
      fr = b
      gr = c
      go to 150
140   l = a
      fl = b
      gl = c
150   w = r - l
      if ( w .le. s ) go to 220
      if ( w .le. u*( abs ( l ) + abs ( r ) ) ) go to 220
      x(4) = a
      f(4) = b
      g(4) = c
      call eql ( x, f, g, h, 4 )
      if ( x(1) .lt. l ) go to 160
      if ( x(1) .gt. r ) go to 160
      b = x(1)
      d4 = x(4) - b
      d3 = x(3) - b
      d2 = x(2) - b
      d34 = x(3) - x(4)
      d42 = x(4) - x(2)
      d23 = x(2) - x(3)
      p2 = 1.0E+00 /d34
      p3 = 1.0E+00 /d42
      p4 = 1.0E+00 /d23
      p2 = (h(2)-h(1))/(d2*(1.0E+00+(d2/d3)*p2*d42+(d2/d4)*p2*d23))
      p3 = (h(3)-h(1))/(d3*(1.0E+00+(d3/d2)*p3*d34+(d3/d4)*p3*d23))
      p4 = (h(4)-h(1))/(d4*(1.0E+00+(d4/d2)*p4*d34+(d4/d3)*p4*d42))
      p2 = p2 + p3 + p4
      if ( p2 .eq. 0. ) go to 160
      a = b - h(1)/p2
      go to 100
160   a = .5*(l+r)
      call cp(a,b,c,j,d,p,n)
      if ( i .ge. j ) go to 170
      r = a
      fr = b
      gr = c
      go to 60
170   l = a
      fl = b
      gl = c
      go to 60
180   e = .5*(l+r)
      return
190   e = l
      return
200   e = r
      return
210   e = d(1)
      e = max ( e, l )
      e = min ( e, r )
      return
220   e = q
      return
      end
      subroutine stm2 ( e, y, z, fy, fz, gy, gz, iy, iz, k, t, d, p, n )

c*********************************************************************72
c
cc STM2
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real a
      real b
      real c
      real d(*)
      real d2
      real d23
      real d3
      real d34
      real d4
      real d42
      real e
      real f(4)
      real fl
      real fr
      real fy
      real fz
      real g(4)
      real gl
      real gr
      real gy
      real gz
      real h(4)
      integer i
      integer iy
      integer iz
      integer j
      integer k
      real l
      integer m
      integer n
      real p(*)
      real p2
      real p3
      real p4
      real q
      real r
      real s
      real t
      real u
      real v
      real w
      real x(4)
      real y
      real z

      l = y
      r = z
      fl = fy
      fr = fz
      gl = gy
      gr = gz
      i = iy
      j = iz
      if ( l .lt. r ) go to 10
      l = z
      r = y
      fl = fz
      fr = fy
      gl = gz
      gr = gy
      i = iz
      j = iy
10    if ( n .eq. 1 ) go to 210
      s = t + t
      u = 1.0E+00
20    u = .5*u
      a = 1.0E+00 + u
      if ( a .gt. 1.0E+00 ) go to 20
      u = 5.*u
      v = .5*u
      if ( i .ge. k ) go to 190
      if ( j .lt. k ) go to 200
30    e = r - l
      if ( e .le. s ) go to 180
      if ( e .le. u*( abs ( l ) + abs ( r ) ) ) go to 180
      if ( j .eq. i+1 ) go to 70
40    a = .5*(l+r)
      call cp(a,b,c,m,d,p,n)
      if ( k .le. m ) go to 50
      l = a
      i = m
      fl = b
      gl = c
      go to 30
50    r = a
      j = m
      fr = b
      gr = c
      go to 30
60    e = r - l
      if ( e .le. s ) go to 180
      if ( e .le. u*( abs ( l ) + abs ( r ) ) ) go to 180
70    x(1) = l
      f(1) = fl
      g(1) = gl
      x(2) = r
      f(2) = fr
      g(2) = gr
      call eql(x,f,g,h,2)
      if ( h(1) .eq. h(2) ) go to 160
      a = x(1) - h(1)*(x(1)-x(2))/(h(1)-h(2))
      q = a
      w = max ( t, v * ( abs ( l ) + abs ( r ) ) )
      if ( abs ( a - l ) .lt. w ) a = l + w
      if ( abs ( a - r ) .lt. w ) a = r - w
      call cp(a,b,c,j,d,p,n)
      if ( i .ge. j ) go to 80
      r = a
      fr = b
      gr = c
      go to 90
80    l = a
      fl = b
      gl = c
90    x(3) = a
      f(3) = b
      g(3) = c
      w = r - l
      if ( w .le. s ) go to 220
      if ( w .le. u*( abs ( l ) + abs ( r ) ) ) go to 220
      call eql(x,f,g,h,3)
      call inp(a,x(1),x(2),x(3),h(1),h(2),h(3),l,r)
      b = l
      if ( abs ( a - l ) .gt. abs ( a - r ) ) b = r
100   q = a
      w = max ( t, v * ( abs ( l ) + abs ( r ) ) )
      if ( abs ( a - l ) .lt. w ) go to 110
      if ( abs ( a - r ) .gt. w ) go to 130
110   if ( a+a .gt. l+r ) go to 120
      a = l + w
      go to 130
120   a = r - w
130   if ( a .le. l ) go to 160
      if ( a .ge. r ) go to 160
      e = .5*e
      if ( e .lt. abs ( b - a ) ) go to 160
      call cp(a,b,c,j,d,p,n)
      if ( i .ge. j ) go to 140
      r = a
      fr = b
      gr = c
      go to 150
140   l = a
      fl = b
      gl = c
150   w = r - l
      if ( w .le. s ) go to 220
      if ( w .le. u*( abs ( l ) + abs ( r ) ) ) go to 220
      x(4) = a
      f(4) = b
      g(4) = c
      call eql(x,f,g,h,4)
      if ( x(1) .lt. l ) go to 160
      if ( x(1) .gt. r ) go to 160
      b = x(1)
      d4 = x(4) - b
      d3 = x(3) - b
      d2 = x(2) - b
      d34 = x(3) - x(4)
      d42 = x(4) - x(2)
      d23 = x(2) - x(3)
      p2 = 1.0E+00 /d34
      p3 = 1.0E+00 /d42
      p4 = 1.0E+00 /d23
      p2 = (h(2)-h(1))/(d2*(1.0E+00+(d2/d3)*p2*d42+(d2/d4)*p2*d23))
      p3 = (h(3)-h(1))/(d3*(1.0E+00+(d3/d2)*p3*d34+(d3/d4)*p3*d23))
      p4 = (h(4)-h(1))/(d4*(1.0E+00+(d4/d2)*p4*d34+(d4/d3)*p4*d42))
      p2 = p2 + p3 + p4
      if ( p2 .eq. 0. ) go to 160
      a = b - h(1)/p2
      go to 100
160   a = .5*(l+r)
      call cp(a,b,c,j,d,p,n)
      if ( i .ge. j ) go to 170
      r = a
      fr = b
      gr = c
      go to 60
170   l = a
      fl = b
      gl = c
      go to 60
180   e = .5*(l+r)
      return
190   e = l
      return
200   e = r
      return
210   e = d(1)
      if ( l .gt. e ) e = l
      if ( r .lt. e ) e = r
      return
220   e = q
      return
      end
      subroutine sto ( a, b, c, i, k, p, x, f, g )

c*********************************************************************72
c
cc STO
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real a
      real b
      real c
      real f(*)
      real g(*)
      integer i
      integer j
      integer k
      real p(*)
      real x(*)

      if ( i .eq. k ) go to 30
      j = i + 1
      if ( p(j) .eq. 0. ) go to 20
      if ( p(j) .eq. 2. ) go to 20
      if ( x(j) .ge. a ) go to 30
10    x(j) = a
      f(j) = b
      g(j) = c
      go to 30
20    p(j) = p(j) + 1.0E+00
      go to 10
30    if ( i .eq. 0 ) return
      if ( p(i) .lt. 2. ) go to 50
      j = i + k
      if ( x(j) .le. a ) return
40    x(j) = a
      f(j) = b
      g(j) = c
      return
50    p(i) = p(i) + 2.
      j = i + k
      go to 40
      end
      subroutine stopit ( dif, size, ndigit, limit )

c*********************************************************************72
c
cc STOPIT performs a test for convergence.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real dif
      real e
      integer i
      integer limit
      integer ndigit
      real size
      real t

      save i

      data i / 0 /

      dif = abs ( dif )
      size = abs ( size )

      if ( i .le. 0 ) then
        t = 10.0E+00**(-ndigit)
      end if

      i = i + 1
      e = 3 * i

      if ( dif .gt. t*size ) go to 20
      e = e + 1.0E+00
      go to 30
20    if ( i .lt. limit ) go to 40
      e = e + 2.
30    dif = -dif
      i = 0
40    size = e

      return
      end
      subroutine stp ( mult, x, w, k, m, n, n5 )

c*********************************************************************72
c
cc STP
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      integer i
      integer k
      integer m
      external mult
      integer n
      integer n5
      real w(*)
      real x(*)

      call px ( x, w(m), k, n, n5 )

      do i = 1, n
        w(i) = x(i)
      end do

      call mult ( x, w )

      call xp ( x, w(m), k, n )

      return
      end
      subroutine svals ( e, d, u, a, n )

c*********************************************************************72
c
cc SVALS computes the eigenvalues of a symmetric matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
C     |         A     --COEFFICIENT MATRIX IN COMPRESSED FORMAT|
C     |                                                        |                                                       |
C     |    OUTPUT:                                             |
C     |                                                        |
C     |         E     --EIGENVALUES                            |
C     |                                                        |
C     |         D     --DIAGONAL OF REDUCED TRIDIAGONAL MATRIX |
C     |                                                        |
C     |         U     --SUPERDIAGONAL OF REDUCED MATRIX 
c
c    Input, integer N, the dimension of the matrix.
c
      implicit none

      real a(*)
      real d(*)
      real e(*)
      integer j
      integer n
      real t
      real u(*)

      call shess ( d, u, a, n )
      j = 1 + (n*(n-1))/2
      t = a(j)
      call tvals ( e, u, d, u, n, a(j) )
      a(j) = t

      return
      end
      subroutine svect ( e, x, a, d, u, w )

c*********************************************************************72
c
cc SVECT computes an eigenvector of a symmetric matrix.
c
c  Discussion:
c
c    The eigenvector corresponds to the input eigenvalue.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
C     |    INPUT:                                              |
C     |                                                        |
C     |         E     --EIGENVALUE ESTIMATE                    |
C     |                                                        |
C     |         A     --SHESS'S OUTPUT                         |
C     |                                                        |
C     |         W     --WORK ARRAY (LENGTH AT LEAST 4N)        |
C     |                                                        |
C     |    OUTPUT:                                             |
C     |                                                        |
C     |         E     --IMPROVED ESTIMATE FOR EIGENVALUE       |
C     |                                                        |
C     |         X     --EIGENVECTOR    
c
      implicit none

      real a(*)
      real d(*)
      real e
      integer i
      integer j
      integer k
      integer l
      integer n
      real t
      real u(*)
      real w(*)
      real x(*)

      t = a(1)

      if ( t .ne. 2234 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must process with coefficient matrix using'
        write(*,*) 'routine shess before computing an eigenvector'
        write(*,*) 'using svect'
        return
      end if

      n = a(2)
      call tvect(e,x,u,d,u,n,w)
      l = 1
      k = n
      j = 1 + (n*(n-3))/2
20    k = k - 1
      l = l + 1
      if ( k .le. 1 ) return

      t = 0.
      do i = k, n
        t = t + x(i)*a(i+j)
      end do

      do i = k, n
        x(i) = x(i) - t*a(i+j)
      end do

      j = j - l
      go to 20
      end
      subroutine svert ( v, n, w )

c*********************************************************************72
c
cc SVERT inverts a symmetric matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input/output, real V((N*(N+1))/1); on input, the factorization information
c    from SFACT.  On output, the NxN inverse matrix, using symmetric storage 
c    of each row, starting at the diagonal.
c
c    Input, integer N, the number of rows and columns in the matrix.
c
c    Workspace, real W(N).                  
c
      implicit none

      integer n

      integer g
      integer h
      integer i
      integer j
      integer k
      integer l
      integer m
      real s
      real t
      real v(*)
      real w(n)

      h = n
      k = 1
10    if ( h .eq. 1 ) go to 40
      s = v(k)
      k = k + h
      g = k
      h = h - 1
      m = h
      if ( s .eq. 0. ) go to 50
      j = 0
20    j = j - m
      m = m - 1
      l = g + m
      t = v(g+j)/s
      do i = g, l
        v(i) = v(i) - t*v(i+j)
      end do
      g = l + 1
      if ( m .gt. 0 ) go to 20
      go to 10
40    if ( v(k) .ne. 0. ) go to 60
50    continue
      write ( *, '(a)' ) ' '
      write(*,*) 'error: zero pivot encountered'
      stop
60    g = n + n

      do m = 1, n

        l = ((g-m)*(m-1))/2
        h = l
        k = m
        do i = m, n
          w(i) = 0.
        end do
        w(m) = 1.0E+00
80      if ( k .eq. n ) go to 100
        t = w(k)/v(k+l)
        j = l
        l = l + n - k
        k = k + 1
        do i = k, n
          w(i) = w(i) - t*v(i+j)
        end do
        go to 80
100     w(n) = w(n)/v(k+l)
110     if ( k .eq. m ) go to 130
        j = k
        k = k - 1
        l = l + k - n
        t = w(k)
        do i = j, n
          t = t - w(i)*v(i+l)
        end do
        w(k) = t/v(k+l)
        go to 110

130     do i = m, n
          v(i+h) = w(i)
        end do

      end do

      return
      end
      function tcon ( l, d, u, b )

c*********************************************************************72
c
cc TCON estimates the condition number of a tridiagonal matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c    William Hager,
c    Condition Estimates,
c    SIAM Journal on Scientific and Statistical Computing,
c    Volume 5, Number 2, June 1984, pages 311-316.
c
c  Parameters:
c
c    Input, real L(*), D(*), U(*), factorization information from TFACT.
c
c    Workspace, real B(N).
c
c    Output, real TCON, the condition number of the matrix.
c
      implicit none

      real b(*)
      real c
      real d(*)
      real e
      integer i
      integer j
      real l(*)
      integer m
      integer n
      real tcon
      real u(*)

      c = d(1)

      if ( abs ( c ) .ne. 1234 .and. abs ( c ) .ne. 1238 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor with tfact before',
     &  'estimating condition'
        stop
      end if

      if ( c .lt. 0.0E+00 ) then
        tcon = -1.0E+00
        return
      end if

      m = 0
      n = d(2)
      c = 1.0E+00 /d(2)
      do j = 1, n
        b(j) = c
      end do
      go to 50
30    do j = 1, n
        b(j) = 0.
      end do
      b(m) = 1.0E+00
50    call tsolve(b,l,d,u,b)

      c = 0.0E+00
      do j = 1, n
        c = c + abs ( b(j) )
        if ( b(j) .lt. 0.0E+00 ) then
          b(j) = -1.0E+00
        else
          b(j) = 1.0E+00
        end if
      end do

      call ttrans(b,l,d,u,b)
c
c  I is the index of the largest entry of B.
c
      i = 1
      do j = 1, n
        if ( abs ( b(i) ) .lt. abs ( b(j) ) ) then
          i = j
        end if
      end do

      if ( m .eq. 0 ) go to 90
      if ( m .eq. i ) go to 100
      if ( e .ge. c ) go to 100
90    m = i
      e = c
      go to 30
100   c = c*d(3)
      c = max ( c, 1.0E+00 )
      tcon = c
      return
      end
      function tdet ( iexp, d )

c*********************************************************************72
c
cc TDET computes the determinant of a tridiagonal matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, integer IEXP, a power of ten that is part of the
c    determinant.
c
c    Input, real D(*), factorization information from TFACT.
c
c    Output, real TDET, the mantissa of the determinant.
c    Determinant = TDET * 10.0^IEXP.
c
      implicit none

      real c
      real d(*)
      real f
      real g
      integer h
      integer i
      integer iexp
      integer l
      integer n
      real t
      real tdet

      iexp=0
      tdet=0.0
      t = d(1)

      if ( abs ( t ) .ne. 1234 .and. abs ( t ) .ne. 1238 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor with tfact before',
     &  ' computing determinant'
        return
      end if

      iexp = 0
      if ( t .lt. 0. ) go to 80
      t = 1.0E+00
      f = 65536.0E+00**4
      g = 1.0E+00 /f
      h = 64
      n = d(2)
      l = 3 + n

      do 50 i = 4, l
        t = t*d(i)
20      if ( abs ( t ) .gt. f ) go to 40
30      if ( abs ( t ) .gt. g ) go to 50
        iexp = iexp - h
        t = t*f
        go to 30
40      iexp = iexp + h
        t = t*g
        go to 20
50    continue

      if ( iexp .ne. 0 ) go to 60
      tdet = t
      return
60    if ( t .eq. 0. ) go to 90
      c = alog10 ( abs ( t ) ) + iexp * alog10 ( 2.0E+00 )
      iexp = c
      c = c - iexp
      if ( c .le. 0.0 ) go to 70
      c = c - 1
      iexp = iexp + 1
70    f = 10.0E+00**c
      if ( t .lt. 0. ) f = -f
      tdet = f
      return
80    tdet = 0.
      return
90    iexp = 0
      go to 80
      end
      subroutine tdg ( e, v, lv, d, u, n )

c*********************************************************************72
c
cc TDG
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input, integer LV, the leading dimension of the array V.
c
      implicit none

      integer lv

      real a
      real b
      real d(*)
      real e(*)
      real f
      integer g
      integer h
      integer i
      integer j
      integer k
      integer l
      integer ll
      integer l1
      integer m
      integer n
      integer ns
      real p
      real q
      real r
      real s
      real t
      real t0
      real t1
      real t2
      real u(*)
      real v(lv,*)
      real w
      real x
      real y
      real z

      e(1) = d(1)
      if ( n .le. 1 ) return
      j = 1
      h = 1
      d(1) = 1.0E+00

      t = 0.
      do i = 2, n
        e(i) = d(i)
        d(i) = 1.0E+00
        t = t + abs ( e(h) ) + abs ( u(h) )
        if ( u(h) .eq. 0. ) j = i
        h = i
      end do

      b = 65536.0E+00**(-3)
      r = 1.0E+00
20    r = .5*r
      s = 1.0E+00 + r
      if ( s .gt. 1.0E+00 ) go to 20
      t0 = r + r
      t2 = t0*t0
      ns = 50*n
      k = n
      ll = 0
      l = n - 1
      l1 = 0
      if ( k .eq. j ) go to 250
      go to 60
30    k = l
      l = l - 1
      j = 0
      do i = 1, l
        if ( u(i) .eq. 0. ) j = i
      end do
      j = j + 1
      if ( k .eq. j ) go to 250

      t = 0.
      r = sqrt ( d(j) )
      do i = j, l
        s = sqrt ( d(i+1) )
        t = t + d(i) * abs ( e(i) ) + r * abs ( u(i) ) * s
        r = s
      end do

60    t = t + abs ( e(k) )
      if ( t .eq. 0. ) t = 1.0E+00
      t1 = 1.0E+00 /(t*t0)
      go to 230
70    ll = ll + 1
      if ( ll .gt. ns ) go to 330
      t = 1.0E+00
      do i = j, k
        if ( d(i) .lt. t ) t = d(i)
      end do
      if ( t .gt. b ) go to 120
      e(j) = e(j)*d(j)
      s = sqrt ( d(j) )
      d(j) = 1.0E+00
      do i = 1, n
        v(i,j) = v(i,j)*s
      end do
      if ( j .gt. 1 ) u(j-1) = u(j-1)*s

      h = j
      g = j + 1
      do m = g, k
        e(m) = e(m)*d(m)
        t = sqrt ( d(m) )
        d(m) = 1.0E+00
        u(h) = s*u(h)*t
        s = t
        do i = 1, n
          v(i,m) = v(i,m)*t
        end do
        h = m
      end do

120   t = d(k)
      s = d(l)
      call eig3(f,f,s*e(l),t*e(k),s*u(l),t*u(l))
      h = j - 1
      i = j
      p = d(j)*e(j) - f
      q = d(j)*u(j)
130   g = h
      h = i
      i = i + 1
      y = d(h)
      z = d(i)
      if ( abs ( p ) .gt. abs ( q ) ) go to 140
      r = y/z
      s = p/q
      t = r*s*s
      if ( t .lt. 1.0E+00 ) go to 190
      t = t/(1.0E+00+t)
      r = z/y
      s = q/p
      go to 150
140   r = z/y
      s = q/p
      t = r*s*s
      if ( t .gt. 1.0E+00 ) go to 180
      t = 1.0E+00 /(1.0E+00+t)
150   d(h) = y*t
      d(i) = z*t
      y = s
      x = r*s
      do m = 1, n
        t = v(m,h)
        s = v(m,i)
        v(m,h) = t + x*s
        v(m,i) = s - y*t
      end do

      if ( h .eq. j ) go to 170
      t = u(g) + x*a
      u(g) = t
      if ( t .eq. 0. ) j = h
170   t = e(h)
      p = u(h)
      q = t + x*p
      r = p - y*t
      w = e(i) - y*p
      s = p + x*e(i)
      u(h) = s - y*q
      e(h) = q + x*s
      e(i) = w - r*y
      if ( i .eq. k ) go to 230
      q = u(i)
      a = x*q
      p = z*w - f
      q = q*z
      go to 130
180   t = t/(1.0E+00+t)
      r = y/z
      s = p/q
      go to 200
190   t = 1.0E+00 /(1.0E+00+t)
200   d(h) = z*t
      d(i) = y*t
      y = s
      x = r*s

      do m = 1, n
        t = v(m,h)
        s = v(m,i)
        v(m,h) = x*t + s
        v(m,i) = y*s - t
      end do

      if ( h .eq. j ) go to 220
      t = x*u(g) + a
      u(g) = t
      if ( t .eq. 0. ) j = h
220   t = e(h)
      p = u(h)
      q = x*t + p
      r = y*p - t
      w = y*e(i) - p
      s = x*p + e(i)
      u(h) = y*s - q
      e(h) = x*q + s
      e(i) = y*w - r
      if ( i .eq. k ) go to 230
      a = u(i)
      q = a
      u(i) = y*a
      p = z*w - y*f
      q = q*z
      go to 130
230   t = d(k)
      s = d(l)
      q = e(k)
      p = u(l)
      if ( ((t*p)*t1)*((s*p)*t1) .gt. 1.0E+00 ) go to 70
      l1 = l1 + 1
      if ( l1 .gt. 30 ) go to 240
      r = max ( abs ( p ), abs ( q ) )
      if ( r .eq. 0. ) go to 240

      if ( abs ( ( s * p ) * (p/r)) .gt. t2 * abs ( (t*q)*(q/r)) ) then
        go to 70
      end if

240   l1 = 0
      e(k) = q*t
      k = l
      l = l - 1
      if ( k .gt. j ) go to 230
250   e(j) = e(j)*d(j)
      if ( j .gt. 2 ) go to 30
      if ( j .eq. 1 ) go to 270
      e(1) = e(1)*d(1)

270   do j = 1, n
        t = sqrt ( d(j) )
        do i = 1, n
          v(i,j) = v(i,j)*t
        end do
      end do
c
c  Index sort E.
c
      call sort2(e,d,u,n)

      do i = 1, n
        j = d(i)
        u(j) = i
      end do

      do j = 1, n

        k = d(j)
        i = u(j)
        d(i) = k
        u(k) = i
        t = e(j)
        e(j) = e(k)
        e(k) = t

        s = 0.
        do i = 1, n
          t = v(i,k)
          v(i,k) = v(i,j)
          v(i,j) = t
          s = s + t*t
        end do

        s = 1.0E+00 / sqrt ( s )
        do i = 1, n
          v(i,j) = s*v(i,j)
        end do

      end do

      d(1) = n
      return
330   k = n - k + 1
      write(*,*) 'since the stopping criterion not satisfied'
      write(*,*) 'after',ns,'iterations, we stop while computing'
      write(*,*) 'eigenvalue number',k
      d(1) = k
      return
      end
      subroutine tdiag ( e, v, lv, l, d, u, n )

c*********************************************************************72
c
cc TDIAG diagonalizes a tridiagonal matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input, integer LV, the leading dimension of the array V.
c
      implicit none

      integer lv

      real d(*)
      real e(*)
      integer i
      integer j
      integer k
      real l(*)
      integer n
      real s
      real t
      real u(*)
      real v(lv,*)

      do i = 1, n
        v(i,1) = 0.0
      end do

      v(1,1) = 1.0E+00

      if ( n .le. 1 ) then
        e(1) = d(1)
        return
      end if

      k = 1
      s = l(1)

      do j = 2, n

        do i = 1, n
          v(i,j) = 0.0
        end do

        if ( s .ne. 0. ) go to 40
        if ( u(k) .ne. 0.0 ) go to 70
        v(j,j) = 1.0E+00
        s = l(j)
        go to 49
40      if ( u(k) .eq. 0.0 ) go to 70
        t = s/u(k)
        if ( t .lt. 0. ) go to 60
        t = sqrt ( t )
        v(j,j) = t
        u(k) = t*u(k)
        if ( j .eq. n ) go to 49
        u(j) = u(j)/t
        s = t*l(j)

49      continue

        k = j

      end do

      call tdg(e,v,lv,d,u,n)
      return
60    continue
      write ( *, '(a)' ) ' '
      write(*,*) 'error: routine tdiag can only be used when'
      write(*,*) 'the cross-diagonal products are nonnegative'
      stop
70    continue
      write ( *, '(a)' ) ' '
      write(*,*) 'error: routine tdiag can only be used if the'
      write(*,*) 'cross-diagonal product a sub i+1,i times a sub i,i+1'
      write(*,*) 'vanishes only when both a sub i+1,i and a sub i,i+1'
      write(*,*) 'are zero'
      stop
      end
      subroutine tfact ( l, d, u, n )

c*********************************************************************72
c
cc TFACT computes the LU factorization of a tridiagonal matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input/output, real L(N), D(3+N), U(N); on input, information defining 
c    the matrix.  On output, factorization information.
c
c    Input, integer N, the number of rows and columns in the matrix.
c
      implicit none

      integer n

      real d(n+3)
      integer i
      integer j
      integer k
      real l(n)
      integer m
      real r
      real s
      real t
      real u(n)
      real w
      real x
      real y

      i = 4 + n
      j = 1 + n
      do k = 1, n
        d(i-k) = d(j-k)
      end do
      x = 1234
      t = u(1)
      s = l(1)
      u(1) = 0.0
      l(1) = 1.0E+00
      if ( u(1) .eq. 1.0E+00 ) x = 1238
      u(1) = t
      l(1) = s
      d(1) = x
      d(2) = n
      r = 0.
      s = 0.
      w = d(4)
      t = w
      m = n - 1

      do k = 1, m

        r = r + abs ( l(k) ) + abs ( w )
        s = max ( s, r )
        r = abs ( u(k) )
        j = k + 4
        w = d(j)

        if ( t .eq. 0. ) then
          d(1) = -x
          t = w
        else
          y = l(k)/t
          t = w - y*u(k)
          l(k) = y
          d(j) = t
        end if

      end do

      r = r + abs ( w )
      s = max ( s, r )
      d(3) = s
      if ( t .eq. 0. ) d(1) = -x

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
      subroutine tmult ( y, x, l, d, u, n )

c*********************************************************************72
c
cc TMULT multiplies a tridiagonal matrix times a vector.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, real Y(N), the product of the matrix times X.
c
c    Input, real X(N), the vector to be multiplied.
c
c    Input, real L(N), D(N), U(N), the lower, diagonal, and upper values
c    that define the matrix.
c
c    Input, integer N, the number of rows and columns in the matrix.
c
      implicit none

      integer n

      real d(n)
      integer i
      integer j
      real l(n)
      integer m
      real u(n)
      real x(n)
      real y(n)

      if ( n .le. 1 ) then
        y(1) = d(1) * x(1)
        return
      end if

      y(1) = d(1) * x(1) + u(1) * x(2)
      m = n - 1
      j = 1

      do i = 2, n - 1
        y(i) = l(j)*x(j) + d(i) * x(i) + u(i) * x(i+1)
        j = i
      end do

      y(n) = l(j) * x(j) + d(n) * x(n)

      return
      end
      subroutine trans ( x, a, b )

c*********************************************************************72
c
cc TRANS solves A' * x = b, after A has been factored.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, real X(N), the solution.
c
c    Input, real A(*), factorization information from FACT.
c   
c    Input, real B(N), the right hand side.
c   
      implicit none

      real a(*)
      real b(*)
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      real t
      real x(*)

      t = a(1)

      if ( abs ( t ) .ne. 1230 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor before solving'
        stop
      end if

      n = a(2)
      m = n + 1
      if ( t .lt. 0.0E+00 ) go to 80
      t = 0.0E+00
      j = 4
      k = 1
20    if ( b(k) .ne. 0.0E+00 ) go to 30
      x(k) = 0.0E+00
      k = k + 1
      if ( k .le. n ) go to 20
      return
30    j = j - m + m*k
40    x(k) = (b(k)-t)/a(j+k)
      if ( k .eq. n ) go to 60

      t = 0.0
      j = j + m
      do i = 1, k
        t = t + a(i+j)*x(i)
      end do

      k = k + 1
      go to 40

60    if ( k .eq. 1 ) then
        return
      end if

      j = j - m
      t = x(k-1)
      do i = k, n
        t = t - x(i)*a(i+j)
      end do
      k = k - 1
      i = a(j)
      x(k) = x(i)
      x(i) = t
      go to 60
80    i = 5 + n + m*n
      l = m
90    i = i - m - 1
      l = l - 1
      if ( a(i) .ne. 0.0E+00 ) go to 90
      j = j + m*(l-k)
      k = l

      do i = 1, n
        x(i) = 0.0E+00
      end do

      x(k) = 1.0E+00
110   if ( k .eq. n ) go to 60
      t = 0.
      j = j + m
      do i = l, k
        t = t - a(i+j)*x(i)
      end do

      k = k + 1
      x(k) = t/a(j+k)
      go to 110
      end
      subroutine tsolve ( x, l, d, u, b )

c*********************************************************************72
c
cc TSOLVE solves a tridiagonal linear system.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, real X(N), the solution.
c
c    Input, real L(*), D(*), U(*), factorization information from TFACT.
c
c    Input, real B(N), the right hand side.
c
      implicit none

      real b(*)
      real d(*)
      integer i
      integer j
      integer k
      real l(*)
      integer m
      integer n
      real t
      real u(*)
      real x(*)

      t = d(1)

      if ( abs ( t ) .ne. 1234 .and. abs ( t ) .eq. 1238 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor with tfact before solving'
        stop
      end if

      n = d(2)
      if ( n .gt. 1 ) go to 30
      if ( t .lt. 0.0E+00 ) go to 20
      x(1) = b(1)/d(4)
      return
20    x(1) = 1.0E+00
      return

30    if ( t .lt. 0.0E+00 ) go to 90
      j = 1
      x(1) = b(1)
      do k = 2, n
        x(k) = b(k) - x(j)*l(j)
        j = k
      end do
      x(n) = x(n)/d(n+3)
      m = n - 1

50    if ( abs ( t ) .ne. 1238 ) then

        do j = 1, m
          k = n - j
          x(k) = (x(k)-x(k+1)*u(k))/d(k+3)
        end do

      else

        do j = 1, m
          k = n - j
          x(k) = x(k)/d(k+3) - x(k+1)*u(k)
        end do

      end if
  
      return

90    j = n + 4

      do i = 1, n
        if ( d(j-i) .eq. 0.0E+00 ) k = i
        x(i) = 0.0E+00
      end do

      k = n - k + 1
      x(k) = 1.0E+00
      n = k
      m = k - 1
      if ( k .gt. 1 ) go to 50
      return
      end
      subroutine ttrans ( x, l, d, u, b )

c*********************************************************************72
c
cc TTRANS solves a transposed tridiagonal linear system.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, real X(N), the solution.
c
c    Input, real L(*), D(*), U(*), factorization information from TFACT.
c   
c    Input, real B(N), the right hand side.
c   
      implicit none

      real b(*)
      real d(*)
      integer i
      integer j
      integer k
      real l(*)
      integer m
      integer n
      real s
      real t
      real u(*)
      real x(*)

      t = d(1)

      if ( abs ( t ) .ne. 1234 .and. abs ( t ) .ne. 1238 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor with tfact before solving'
        stop
      end if

      n = d(2)
      if ( n .gt. 1 ) go to 30

      if ( t .lt. 0.0E+00 ) go to 20
      x(1) = b(1)/d(4)
      return
20    x(1) = 1.0E+00
      return
30    if ( t .lt. 0.0E+00 ) go to 90
      j = 1
      x(1) = b(1) / d(4)

      if ( t .eq. 1238 ) go to 50
      do k = 2, n
        x(k) = (b(k) - x(j)*u(j))/d(k+3)
        j = k
      end do
      go to 70

50    s = d(4)

      do k = 2, n
        t = d(k+3)
        x(k) = (b(k)-x(j)*u(j)*s)/t
        s = t
        j = k
      end do

70    m = n - 1
      do j = 1, m
        k = n - j
        x(k) = x(k) - x(k+1)*l(k)
      end do

      return
90    continue

      do i = 1, n
        if ( d(i+3) .eq. 0.0E+00 ) j = i
        x(i) = 0.0E+00
      end do

      x(j) = 1.0E+00
      if ( j .eq. n ) go to 70
      m = j + 1
      if ( t .eq. 1238 ) go to 120
      do k = m, n
        x(k) = -x(j)*u(j)/d(k+3)
        j = k
      end do
      go to 70

120   s = d(m+2)

      do k = m, n
        t = d(k+3)
        x(k) = -x(j)*u(j)*s/t
        s = t
        j = k
      end do

      go to 70
      end
      subroutine tval ( e, k, l, d, u, n, w )

c*********************************************************************72
c
cc TVAL computes the k-th smallest or largest eigenvalue of a tridiagonal matrix whose cross-diagonal products are nonnegative.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real d(*)
      real e
      integer i
      integer k
      real l(*)
      integer m
      integer n
      real s
      real t
      real u(*)
      real w(*)
      real y
      real z

      if ( n .le. 1 ) then
        e = d(1)
        return
      end if

      m = n - 1
      y = d(1)
      z = d(1)
      s = 0.

      do i = 1, m

        w(i) = l(i)*u(i)

        if ( w(i) .lt. 0. ) then
          write ( *, '(a)' ) ' '
          write(*,*) 'error: routine tval can only be applied'
          write(*,*) 'when the cross-diagonal products are nonnegative'
          stop
        end if

        s = s + abs ( u(i) )
        t = d(i) + s
        z = max ( z, t )
        t = d(i) - s
        y = min ( y, t )
        s = abs ( l(i) )

      end do

      t = d(n) + s
      z = max ( z, t )
      t = d(n) - s
      y = min ( y, t )
      m = k
      if ( k .le. 0 ) then
        m = n + k + 1
      end if

      t = 1.0E+00

10    continue

      t = .5*t
      s = 1.0E+00 + t
      if ( s .gt. 1.0E+00 ) go to 10

      t = 8.*n*t* max ( abs ( y ), abs ( z ) )

      if ( t .eq. 0.0E+00 ) then
        e = d(1)
        return
      end if

      call stm(e,y,z,m,t,d,w,n)

      return
      end
      subroutine tvals ( e, l, d, u, n, w )

c*********************************************************************72
c
cc TVALS computes the eigenvalues of a tridiagonal matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real d(*)
      real e(*)
      integer i
      integer j
      integer k
      real l(*)
      integer m
      integer n
      real q
      real r
      real s
      real t
      real u(*)
      real v
      real w(*)

      if ( n .le. 1 ) then
        e(1) = d(1)
        return
      end if

      m = n - 1

      do i = 1, m
        w(i) = l(i)*u(i)
      end do

      call pwk(e,d,w,n,w)

      if ( w(1) .eq. 1.0E+00 ) go to 70
      if ( w(1) .eq. 2. ) go to 80
      t = 1.0E+00
30    t = t/2.
      s = 1.0E+00 + t
      if ( s .gt. 1.0E+00 ) go to 30
      t = sqrt ( t )

      do i = 1, m
        w(i) = l(i)*u(i)
      end do

      do i = 1, n

        q = 0.
        j = 0
        v = e(i)
50      call newton(v,s,k,d,w,n)
        j = j + 1

        if ( j .le. 2 ) then
          s = abs ( s )
          r = abs ( v ) + s
          if ( q .lt. r ) q = r
          if ( s .gt. t*q ) go to 50
          e(i) = v
       end if

      end do
!
!  Sort the values.
!
      call sort ( e, w, n )

      return
70    continue
      write ( *, '(a)' ) ' '
      write(*,*) 'the qr method did not converge'
      return
80    continue
      write ( *, '(a)' ) ' '
      write(*,*) 'error:to use routine pvals, the cross-diagonal'
      write(*,*) 'product l(i)*u(i) must be nonnegative for every i'
      return
      end
      subroutine tvect ( e, x, l, d, u, n, w )

c*********************************************************************72
c
cc TVECT computes the eigenvectors of a tridiagonal matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      integer n

      real d(n)
      real e
      integer f
      integer g
      integer h
      integer i
      integer j
      integer k
      real l(n)
      integer m
      real o
      real p
      real q
      real r
      real s
      real t
      real u(n)
      real v
      real w(*)
      real x(*)
      real y
      real z

      if ( n .le. 1 ) then
        e = d(1)
        x(1) = 1.0E+00
        return
      end if

      m = n - 1
      j = 2
      do i = 1, n - 1
        x(i) = 0.
        w(j) = u(i)
        w(j+1) = d(i) - e
        w(j+2) = l(i)
        j = j + 4
      end do

      w(j) = 0.
      w(j+1) = d(n) - e
      o = 65536.0E+00**(-4)
      t = .5*o
      s = t
30    t = .5*t
      p = s
      s = s + t
      if ( s .ge. o ) go to 40
      if ( s+t .gt. s ) go to 30
40    r = w(3)
      v = abs ( r ) + abs ( w(4) )
      f = 4
      m = 4*n - 3
      g = -2
50    g = g + 4
      if ( g .gt. m ) go to 110
      h = g - 1
      i = g + 2
      j = g + 5
      q = w(i)
      y = abs ( q )
      z = abs ( r )
      if ( z .ge. y ) go to 80
      if ( v .ge. y ) go to 60
      v = y
      f = i
60    t = w(g)
      w(g) = w(j)
      w(j) = t
      t = r/q
      k = g + 1
      w(k) = q
      k = j - 1
      s = w(k)
      if ( s .eq. 0. ) go to 70
      if ( s .eq. o ) s = p
      w(k) = -s*t
      w(h) = s
      go to 100
70    w(k) = s
      w(h) = o
      go to 100
80    w(h) = 0.
      if ( v .ge. z ) go to 90
      v = z
      f = i
90    if ( r .eq. 0. ) go to 120
      t = q/r
100   r = w(j) - t*w(g)
      w(j) = r
      w(i) = t
      go to 50
110   if ( abs ( r ) .ge. v ) go to 120
      v = r
      f = j + 1
120   j = f/4
      x(j) = 1.0E+00
      if ( j .eq. 1 ) go to 140
      k = f - 5
      j = j - 1
      x(j) = (x(j)-w(k-1)*x(j+1))/w(k)
130   if ( j .eq. 1 ) go to 140
      j = j - 1
      k = k - 4
      t = w(k-2)
      if ( t .eq. o ) t = 0.
      x(j) = (x(j)-w(k-1)*x(j+1)-t*x(j+2))/w(k)
      go to 130
140   if ( v .eq. 0. ) go to 230
      s = 0.
      do i = 1, n
        t = abs ( x(i) )
        if ( t .gt. s ) s = t
      end do
      r = 0.
      s = 1.0E+00 /s
      do i = 1, n
        t = s*x(i)
        r = r + t*t
        x(i) = t
      end do
      k = 0
      j = 1
      y = x(1)
170   k = k + 4
      i = j
      j = j + 1
      s = w(k)
      w(k) = y
      y = x(j)
      if ( w(k-3) .eq. 0. ) go to 180
      t = x(j)
      x(j) = x(i)
      x(i) = t
180   x(j) = x(j) - s*x(i)
      if ( j .lt. n ) go to 170
      s = x(j)/w(k+3)
      x(j) = s
      t = abs ( s )
      v = s*y
      j = j - 1
      k = k - 1
      s = (x(j)-w(k-1)*s)/w(k)
      x(j) = s
      if ( abs ( s ) .gt. t ) t = abs ( s )
      v = v + s*w(k+1)
190   if ( j .eq. 1 ) go to 200
      j = j - 1
      k = k - 4
      z = w(k-2)
      if ( z .eq. o ) z = 0.
      s = (x(j)-w(k-1)*s-z*x(j+2))/w(k)
      x(j) = s
      if ( abs ( s ) .gt. t ) t = abs ( s )
      v = v + s*w(k+1)
      go to 190
200   if ( v .ne. 0. ) v = r/v
      s = 0.
      t = 1.0E+00 /t
      z = 0.
      do i = 1, n
        s = s + (x(i)*v)**2
        z = z + (t*x(i))**2
      end do
      t = t / sqrt ( z )
      do i = 1, n
        x(i) = t*x(i)
      end do
      if ( r+r .ge. s ) e = e + v
      return

230   continue

      t = 0.0E+00
      do i = 1, n
        s = abs ( x(i) )
        t = max ( t, s )
      end do
      t = 1.0E+00 / t

      z = 0.0E+00
      do i = 1, n
        z = z + ( t * x(i) )**2
      end do

      t = t / sqrt ( z )
      do i = 1, n
        x(i) = t * x(i)
      end do

      return
      end
      subroutine tvert ( v, lv, l, d, u )

c*********************************************************************72
c
cc TVERT inverts a tridiagonal matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, real V(LV,N), the NxN inverse matrix.
c
c    Input, integer LV, the leading dimension of the array V.
c
c    Input, real L(*), D(*), U(*), factorization information from TFACT.
c
      implicit none

      integer lv

      real d(*)
      integer h
      integer i
      integer j
      integer k
      real l(*)
      integer m
      integer n
      integer o
      real t
      real u(*)
      real v(lv,*)

      t = d(1)

      if ( abs ( t ) .ne. 1234 .and. abs ( t ) .ne. 1238 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must factor with tfact before inverting'
        stop
      end if

      if ( t .le. 0.0E+00 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: matrix has no inverse'
        stop
      end if

      n = d(2)

      if ( n .eq. 1 ) then
        v(1,1) = 1.0E+00 / d(4)
        return
      end if

      m = n - 1

      do h = 1, n

        o = h - 1
        do i = 1, h - 1
          v(i,h) = 0.0E+00
        end do
        v(h,h) = 1.0E+00

        o = h + 1
        j = h
        do k = o, n
          v(k,h) = - v(j,h) * l(j)
          j = k
        end do

        v(n,h) = v(n,h) / d(n+3)

        if ( t .ne. 1238 ) then

          do j = 1, n - 1
            k = n - j
            v(k,h) = ( v(k,h)-v(k+1,h) * u(k) ) / d(k+3)
          end do

        else

          do j = 1, n - 1
            k = n - j
            v(k,h) = v(k,h) / d(k+3) - v(k+1,h) * u(k)
          end do

        end if

      end do

      return
      end
      subroutine under ( x, a, b )

c*********************************************************************72
c
cc UNDER computes the least squares (minimum norm) solution to an underdetermined linear system.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, real X(N), the solution.
c
c    Input, real A(*), the factorization information from QR for the MxN matrix.
c
c    Input, real B(M), the right hand side.
c
      implicit none

      real a(*)
      real b(*)
      integer h
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      integer o
      real t
      real x(*)

      t = a(1)

      if ( abs ( t ) .ne. 3230 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must qr factor coefficient matrix'
        write(*,*) 'before solving system'
        stop
      end if

      if ( t .le. 0. ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'singular system - compute regularized solution'
        write(*,*) '(see section 6-11)'
        stop
      end if

      m = a(2)
      n = a(3)

      if ( m .lt. n ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: the number of equations is greater than'
        write(*,*) 'the number of unknowns. for an underdetermined'
        write(*,*) 'system, there are fewer equations than unknowns'
        stop
      end if

      if ( m .eq. 1 ) then
        x(1) = b(1) / a(4)
        return
      end if

      do i = 1, n
        x(i) = b(i)
      end do

      o = m + 1
      l = n
      if ( m .eq. n ) then
        l = n - 1
      end if

      k = 3 + o * n
      do i = 1, l
        j = a(i+k)
        t = x(j)
        x(j) = x(i)
        x(i) = t
      end do

      k = 3
      j = 1
      t = x(1)

80    continue

      x(j) = t / a(j+k)
      k = k + o
      if ( j .eq. n ) go to 100
      h = j
      j = j + 1
      t = x(j)
      do i = 1, h
        t = t - a(k+i)*x(i)
      end do
      go to 80

100   continue

      j = n + 1
      do i = n + 1, m
        x(i) = 0.
      end do
      k = 4 + l*o

130   continue

      k = k - o
      t = 0.
      do i = l, m
        t = t + a(k+i)*x(i)
      end do

      do i = l, m
        x(i) = x(i) - t*a(k+i)
      end do

      l = l - 1

      if ( l .gt. 0 ) then
        go to 130
      end if

      return
      end
      subroutine unpack ( a, la, n )

c*********************************************************************72
c
cc UNPACK reverses the operation of PACK by unpacking a square matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input/output, real A(*).
c    On input, the NxN entries of A were packed into the front of an LAxN array.
c    On output, the NxN entries of A have been unpacked into the LAxN array.
c
c    Input, integer LA, the leading dimension of the array.
c    N <= LA.
c
c    Input, integer N, the number of rows and columns in the matrix.
c
      implicit none

      real a(*)
      integer i
      integer ii
      integer j
      integer jj
      integer la
      integer n

      if ( la .eq. n ) then
        return
      end if

      if ( la .lt. n ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: la in unpack must be .ge. n argument!'
        return
      end if

      do jj = 1, n
        j = n + 1 - jj
        do ii = 1, n
          i = n + 1 - ii
          a((j-1)*la+i) = a((j-1)*n+i)
        end do
      end do

      do j = 1, n
        do i = n + 1, la
          a((j-1)*la+i) = 0.0E+00
        end do
      end do

      return
      end
      subroutine update ( dif, size, new, old, n )

c*********************************************************************72
c
cc UPDATE updates a vector, and computes the difference and norm.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Output, real DIF, the L1 norm of NEW - OLD.
c
c    Output, real SIZE, the L1 norm of NEW.
c
c    Input, real NEW(N), the new vector.
c
c    Input/output, real OLD(N); on input, the old vector.
c    On output, a copy of NEW.
c
c    Input, integer N, the size of the vector.
c
      implicit none

      integer n

      real dif
      integer i
      real new(n)
      real old(n)
      real size

      dif = 0.0E+00
      size = 0.0E+00
      do i = 1, n
        dif = dif + abs ( new(i) - old(i) )
        size = size + abs ( new(i) )
        old(i) = new(i)
      end do

      return
      end
      subroutine vals ( e, a, la, n, v )

c*********************************************************************72
c
cc VALS computes all eigenvalues of a general real matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      real a(*)
      real amag
      complex e(*)
      integer i
      integer j
      integer k
      integer l
      integer la
      integer m
      integer n
      integer o
      integer p
      real r
      real s
      real t
      real v(*)
      complex z

      call ahess(a,la,n,e)
      o = n + 1
      p = n + 2
      l = 1 + n*p

      do i = 1, l
        v(i) = a(i)
      end do

      if ( n .eq. 1 ) then
        e(1) = a(3)
        return
      end if

      if ( n .gt. 2 ) go to 30
      a(1) = v(3)
      a(2) = 0.0E+00
      a(3) = v(4)
      a(4) = 0.0E+00
      a(5) = v(6)
      a(6) = 0.0E+00
      a(8) = 0.0E+00
      i = 9
      go to 80
30    j = 1
      k = n * n - 1
      m = 3
      l = 4
40    do i = m, l
        a(j) = v(i)
        a(j+1) = 0.0E+00
        j = j + 2
      end do
      m = m + o
      l = l + p
      if ( l .lt. k ) go to 40
      a(j) = a(m)
      a(j+1) = 0.0E+00
      j = j + 4*n
      m = m + o
      l = l + p
      i = j
60    l = l - 1
      j = j - 2
      a(j) = a(l)
      a(j+1) = 0.0E+00
      if ( l .gt. m ) go to 60
      l = l - 1
      m = m - n
70    l = l - 1
      j = j - 2
      a(j) = a(l)
      a(j+1) = 0.0E+00
      if ( l .gt. m ) go to 70
80    j = i + n
      k = j + n
      call vls(e,a,n,a(i),a(j),a(k))
      i = o
90    i = i - 1
      if ( i .le. 1 ) go to 120
      z = conjg ( e(i) )
      r = abs ( aimag ( z ) )
      s = r
      l = i - 1

      do j = 1, l
        t = amag ( e(j) - z )
        if ( t .lt. s ) then
          k = j
          s = t
        end if
      end do

      if ( r .gt. s ) go to 110
      e(i) = real ( z )
      go to 90
110   e(k) = e(l)
      e(l) = z
      i = l
      go to 90
120   if ( i .eq. 1 ) then
        e(i) = real ( e(i) )
      end if

      return
      end
      subroutine vect ( e, x, v, w )

c*********************************************************************72
c
cc VECT computes an eigenvector of a general real matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      complex e
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      integer o
      integer p
      real r
      real s
      real t
      real v(*)
      real w(*)
      real x(*)

      t = v(1)

      if ( t .ne. 2230 .and. t .ne. 2231 ) then
        write ( *, '(a)' ) ' '
        write(*,*) 'error: must process the coefficient matrix using'
        write(*,*) 'one of the following routines before computing'
        write(*,*) 'an eigenvector using vect: hess, ahess, or vals'
        stop
      end if

      n = v(2)

      if ( n .le. 1 ) then
        e = v(3)
        x(1) = (1.0E+00,0.0E+00)
        return
      end if

      o = n + 1
      if ( aimag(e) .ne. 0.0E+00 ) go to 140
      p = (n*(n+3))/2
      j = 2
      k = 1
      l = 2
      m = 1

30    continue

      do i = m, l
        w(i) = v(i+j)
      end do
      m = l + 1
      j = j + n - k
      k = k + 1
      l = m + k
      if ( l .lt. p ) go to 30
      l = l - 1
      if ( l .lt. p ) go to 30
      t = real ( e )
      call evect(t,x,w,n)
      e = t
      k = n
      j = 3 + (n-2)*o
60    j = j - o
      k = k - 1
      if ( k .le. 1 ) go to 90
      if ( v(j+k) .eq. 0. ) go to 60
      t = 0.0E+00
      do i = k, n
        t = t + x(i)*v(i+j)
      end do
      do i = k, n
        x(i) = x(i) - t*v(i+j)
      end do
      go to 60
90    if ( v(1) .eq. 2230 ) go to 110
      j = 1 + n*o

      do i = 1, n
        x(i) = x(i)*v(i+j)
      end do

110   s = 0.0E+00
      do i = 1, n
        t = abs ( x(i) )
        if ( t .gt. s ) s = t
      end do
      if ( s .ne. 0. ) s = 1.0E+00 /s
      i = n + n
      j = n
130   x(i) = 0.0E+00
      x(i-1) = s*x(j)
      i = i - 2
      j = j - 1
      if ( j .gt. 0 ) go to 130
      return
140   p = n + 2
      j = 1
      k = 2 + n*n
      m = 3
      l = 4
150   do i = m, l
        w(j) = v(i)
        w(j+1) = 0.0E+00
        j = j + 2
      end do
      m = m + o
      l = l + p
      if ( m .lt. k ) go to 150
      l = l - 1
      if ( m .eq. k ) go to 150
      call cevect(e,x,w,n)

      j = 1
      do i = 1, n
        w(i) = x(j)
        w(i+n) = x(j+1)
        j = j + 2
      end do

      k = n
      j = 3 + (n-2)*o
180   j = j - o
      k = k - 1
      if ( k .le. 1 ) go to 210
      if ( v(j+k) .eq. 0. ) go to 180

      t = 0.
      s = 0.
      do i = k, n
        r = v(i+j)
        t = t + w(i)*r
        s = s + w(i+n)*r
      end do

      do i = k, n
        r = v(i+j)
        l = i + n
        w(i) = w(i) - t*r
        w(l) = w(l) - s*r
      end do

      go to 180

210   if ( v(1) .eq. 2230 ) go to 230

      j = 1 + n*o
      do i = 1, n
        r = v(i+j)
        l = i + n
        w(i) = w(i)*r
        w(l) = w(l)*r
      end do

230   s = 0.0E+00
      do i = 1, n
        t = abs ( w(i) ) + abs ( w(i+n) )
        if ( t .gt. s ) s = t
      end do
      if ( s .ne. 0. ) s = 1.0E+00 /s
      i = n + n
      j = n
250   x(i) = s*w(j+n)
      x(i-1) = s*w(j)
      i = i - 2
      j = j - 1
      if ( j .gt. 0 ) go to 250
      return
      end
      subroutine vert ( v, lv, n, w )

c*********************************************************************72
c
cc VERT inverts a general matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c  V      ARRAY CONTAINING MATRIX
c
c    Input, integer LV, the leading dimension of the array V.
c
c  N      DIMENSION OF MATRIX STORED IN ARRAY V
c
c    Workspace, real W(N-1).
c
c  Output:
c
c  V      INVERSE
c
      implicit none

      integer lv
      integer n

      integer i
      integer j
      integer k
      integer l
      integer m
      integer p
      real s
      real t
      real v(lv,*)
      real w(n-1)

      if ( n .eq. 1 ) go to 110
      l = 0
      m = 1
10    if ( l .eq. n ) go to 90
      k = l
      l = m
      m = m + 1
      p = l
      if ( m .gt. n ) go to 30

      s = abs ( v(l,l) )
      do i = m, n
        t = abs ( v(i,l) )
        if ( t .gt. s ) then
          p = i
          s = t
        end if
      end do

      w(l) = p
30    s = v(p,l)
      v(p,l) = v(l,l)
      if ( s .eq. 0.0E+00 ) go to 120
      v(l,l) = -1.0E+00
      s = 1.0E+00 /s
      do i = 1, n
        v(i,l) = -s*v(i,l)
      end do
      j = l
50    j = j + 1
      if ( j .gt. n ) j = 1
      if ( j .eq. l ) go to 10
      t = v(p,j)
      v(p,j) = v(l,j)
      v(l,j) = t
      if ( t .eq. 0. ) go to 50
      if ( k .eq. 0 ) go to 70
      do i = 1, k
        v(i,j) = v(i,j) + t*v(i,l)
      end do
70    v(l,j) = s*t
      if ( m .gt. n ) go to 50
      do i = m, n
        v(i,j) = v(i,j) + t*v(i,l)
      end do
      go to 50
90    l = w(k)

      do i = 1, n
        t = v(i,l)
        v(i,l) = v(i,k)
        v(i,k) = t
      end do

      k = k - 1
      if ( k .gt. 0 ) go to 90
      return
110   if ( v(1,1) .eq. 0. ) go to 120
      v(1,1) = 1.0E+00 /v(1,1)
      return
120   continue
      write ( *, '(a)' ) ' '
      write(*,*) 'error: matrix has no inverse'
      stop
      end
      subroutine vls ( e, a, n, d, p, w )

c*********************************************************************72
c
cc VLS computes all eigenvalues of a complex Hessenberg matrix.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      complex a(*)
      real amag
      real b
      real d(*)
      complex e(*)
      integer f
      integer f0
      integer g
      integer g0
      integer h
      integer i
      integer j
      integer k
      integer kl
      integer km
      integer ks
      integer k0
      integer l
      integer ll
      integer l0
      integer l1
      integer m
      integer n
      integer ns
      integer o
      real p(*)
      real q
      real r
      real s
      real sqr
      real t
      real t0
      real t1
      real t2
      real u
      complex w(*)
      complex x
      complex y
      complex z
      complex z1
      complex z2
      complex z3
      complex z4

      do i = 1, n
        d(i) = 1.0E+00
      end do

      if ( n .eq. 1 ) go to 330
      b = 65536.0E+00**(-3)
      t = 1.0E+00
20    t = .5*t
      s = 1.0E+00 + t
      if ( s .gt. 1.0E+00 ) go to 20
      t0 = t + t
      t2 = t0*t0
      ns = 50*n
      ll = 0
      kl = (n*(n+3))/2 - 1
      m = n
      km = kl - m
30    i = km
      j = m
40    if ( amag ( a(i) ) .eq. 0. ) go to 50
      i = i - j
      j = j - 1
      if ( j .gt. 1 ) go to 40
50    k0 = j
      l = j + 1
      g = i + l
      f = g - 1
      f0 = f
      g0 = g
      l0 = l
      if ( j .eq. m ) go to 320
      do i = j, m
        p(i) = sqrt ( d(i) )
      end do
      s = 0.
      k = j
70    h = f - j
      t = p(k)
      do i = f, g
        s = s + t * amag ( a(i) ) * p(i-h)
      end do
      k = l
      l = l + 1
      f = f + k
      g = g + l
      if ( k .lt. m ) go to 70
      if ( k .gt. m ) go to 90
      g = g - 1
      go to 70
90    if ( s .eq. 0. ) s = 1.0E+00
      t1 = 1.0E+00 /(s*t0)
      l1 = 0
      go to 300
100   ll = ll + 1
      if ( ll .gt. ns ) go to 340
      t = 1.0E+00
      do i = k0, m
        if ( d(i) .lt. t ) t = d(i)
      end do
      if ( t .gt. b ) go to 150
      do i = k0, m
        p(i) = sqrt ( d(i) )
        d(i) = 1.0E+00
      end do
      f = f0
      g = g0
      l = l0
      k = k0
130   h = f - k0
      t = p(k)
      do i = f, g
        a(i) = a(i)*t*p(i-h)
      end do
      k = l
      l = l + 1
      f = f + k
      g = g + l
      if ( k .lt. m ) go to 130
      if ( k .gt. m ) go to 150
      g = g - 1
      go to 130
150   s = d(m-1)
      t = d(m)
      z1 = a(km-1)*s
      z2 = a(km)*s
      z3 = a(kl-1)*t
      z4 = a(kl)*t
      call ceig(z,z,z1,z2,z3,z4)
160   k = k0
      l = l0
      f = f0
      g = g0
      ks = k
      z1 = d(k)*a(f) - z
      z2 = d(k)*a(g)
170   t = amag ( z1 ) + amag ( z2 )
      if ( t .eq. 0. ) go to 190
      t = 1.0E+00 /t
      q = sqr(z1,t)
      r = sqr(z2,t)
      q = d(k)*q
      r = d(l)*r
      if ( q .gt. r ) go to 180
      z4 = z1/z2
      z3 = (d(k)/d(l)) * conjg ( z4 )
      e(k) = z3
      w(k) = z4
      p(k) = 1.0E+00
      s = r/(q+r)
      r = d(k)
      d(k) = d(l)*s
      q = d(l)
      d(l) = r*s
      go to 200
180   z4 = z2/z1
      z3 = (d(l)/d(k)) * conjg ( z4 )
      e(k) = z3
      w(k) = z4
      p(k) = 0.
      s = q/(q+r)
      d(k) = d(k)*s
      q = d(l)
      d(l) = d(l)*s
      go to 200
190   p(k) = 0.
      z3 = (0.,0.)
      z4 = (0.,0.)
      e(k) = z3
      w(k) = z4
      q = d(l)
200   if ( k .gt. ks ) go to 220
      if ( p(k) .eq. 0. ) go to 210
      y = a(g) + z3*a(f)
      a(g) = z4*a(g) - a(f)
      a(f) = y
      go to 240
210   y = a(f) + z3*a(g)
      a(g) = a(g) - z4*a(f)
      a(f) = y
      go to 240
220   i = g - 1
      if ( p(k) .eq. 0. ) go to 230
      x = a(h)*z3 + x
      a(h) = x
      y = a(g) + z3*a(i)
      a(g) = z4*a(g) - a(i)
      a(i) = y
      if ( amag ( x ) .ne. 0. ) go to 240
      k0 = k
      l0 = l
      f0 = i
      g0 = g
      go to 240
230   x = x*z3 + a(h)
      a(h) = x
      y = a(i) + z3*a(g)
      a(g) = a(g) - z4*a(i)
      a(i) = y
      if ( amag ( x ) .ne. 0. ) go to 240
      k0 = k
      l0 = l
      f0 = i
      g0 = g
240   j = k
      k = l
      l = l + 1
      do i = ks, j
        o = g + i
        h = o + 1
        if ( p(i) .ne. 0. ) then
          y = a(h) + e(i)*a(o)
          a(h) = w(i)*a(h) - a(o)
          a(o) = y
        else
          y = a(o) + e(i)*a(h)
          a(h) = a(h) - w(i)*a(o)
          a(o) = y
        end if
      end do

      z3 = conjg ( z3 )
      z4 = conjg ( z4 )
      if ( p(j) .eq. 0. ) go to 280
      z1 = q*a(g+k) - z*w(j)
      do i = f, g
        h = i + k
        y = a(h) + z3*a(i)
        a(h) = z4*a(h) - a(i)
        a(i) = y
      end do
      if ( l .gt. m ) go to 300
      f = f + k
      h = g
      g = g + l
      z2 = q*a(g)
      x = a(g)
      a(g) = z4*a(g)
      go to 170
280   z1 = q*a(g+k) - z
      do i = f, g
        h = i + k
        y = a(i) + z3*a(h)
        a(h) = a(h) - z4*a(i)
        a(i) = y
      end do
      if ( l .gt. m ) go to 300
      f = f + k
      h = g
      g = g + l
      z2 = q*a(g)
      x = z3*a(g)
      go to 170
300   t = d(m)
      s = d(m-1)
      q = amag ( a(km) )
      if ( ((t*q)*t1)*((s*q)*t1) .gt. 1.0E+00 ) go to 100
      l1 = l1 + 1
      if ( l1 .gt. 30 ) go to 310
      r = amag ( a(kl) )
      u = max ( q, r )
      if ( u .eq. 0. ) go to 310
      if ( (s*q)*(q/u) .gt. t2*(t*r)*(r/u) ) go to 100
310   l1 = 0
      e(m) = a(kl)*d(m)
      kl = km - 1
      km = km - m
      m = m - 1
      if ( m .gt. k0 ) go to 300
320   e(m) = a(kl)*d(m)
      kl = km - 1
      km = km - m
      m = m - 1
      if ( m .gt. 1 ) go to 30
330   e(1) = a(1)*d(1)
      p(1) = n
      return
340   m = n - m + 1
      write(*,*) 'since the stopping criterion not satisfied'
      write(*,*) 'after',ns,'iterations, we stop while computing'
      write(*,*) 'eigenvalue number',m
      p(1) = m
      return
      end
      subroutine whatis ( dif, t )

c*********************************************************************72
c
cc WHATIS prints the iteration number, iteration difference, and stopping criterion.
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input, real DIF, the iteration difference.
c
c    Input, real T, ?
c
      implicit none

      real dif
      integer i
      real s
      real t

      s = abs ( dif )
      i = int ( abs ( t ) / 3.0E+00 )
      write ( *, * ) 'whatis - iteration:', i, ' difference:', s
      i = abs ( t ) - 3 * i - 1

      if ( i .eq. 0 ) then
        write(*,*) 'whatis - current error is below user tolerance.'
      else if ( i .gt. 0 ) then
        write(*,*) 'whatis - reached maximum number of iterations.'
      end if

      return
      end
      subroutine xp ( y, z, k, n )

c*********************************************************************72
c
cc XP
c
c  Modified:
c
c    16 November 2011
c
c  Author:
c
c    William Hager
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
      implicit none

      integer n

      integer i
      integer j
      integer k
      integer l
      real t
      real y(n)
      real z(*)

      if ( k .eq. 0 ) then
        return
      end if

      l = 0

      do j = 1, k

        t = 0.0E+00
        do i = j, n
          t = t + z(l+i) * y(i)
        end do

        do i = j, n
          y(i) = y(i) - t * z(l+i)
        end do

        l = l + n - j

      end do

      return
      end
