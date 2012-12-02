      subroutine bakvec(nm,n,t,e,m,z,ierr)

c*********************************************************************72
c
cc BAKVEC determines eigenvectors by reversing the FIGI transformation.
c
c     this routine forms the eigenvectors of a nonsymmetric
c     tridiagonal matrix by back transforming those of the
c     corresponding symmetric matrix determined by  figi.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        t contains the nonsymmetric matrix.  its subdiagonal is
c          stored in the last n-1 positions of the first column,
c          its diagonal in the n positions of the second column,
c          and its superdiagonal in the first n-1 positions of
c          the third column.  t(1,1) and t(n,3) are arbitrary.
c
c        e contains the subdiagonal elements of the symmetric
c          matrix in its last n-1 positions.  e(1) is arbitrary.
c
c        m is the number of eigenvectors to be back transformed.
c
c        z contains the eigenvectors to be back transformed
c          in its first m columns.
c
c     on output
c
c        t is unaltered.
c
c        e is destroyed.
c
c        z contains the transformed eigenvectors
c          in its first m columns.
c
c        ierr is set to
c          zero       for normal return,
c          2*n+i      if e(i) is zero with t(i,1) or t(i-1,3) non-zero.
c                     in this case, the symmetric matrix is not similar
c                     to the original matrix, and the eigenvectors
c                     cannot be found by this program.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,m,n,nm,ierr
      double precision t(nm,3),e(n),z(nm,m)

      ierr = 0
      if (m .eq. 0) go to 1001
      e(1) = 1.0d0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
         if (e(i) .ne. 0.0d0) go to 80
         if (t(i,1) .ne. 0.0d0 .or. t(i-1,3) .ne. 0.0d0) go to 1000
         e(i) = 1.0d0
         go to 100
   80    e(i) = e(i-1) * e(i) / t(i-1,3)
  100 continue
c
      do 120 j = 1, m
c
         do 120 i = 2, n
         z(i,j) = z(i,j) * e(i)
  120 continue
c
      go to 1001
c     .......... set error -- eigenvectors cannot be
c                found by this program ..........
 1000 ierr = 2 * n + i
 1001 return
      end
      subroutine balanc ( nm, n, a, low, igh, scale )

c*********************************************************************72
c
cc BALANC balances a real matrix before eigenvalue calculations.
c
c  Discussion:
c
c    This routine is a translation of the ALGOL procedure BALANCE,
c    num. math. 13, 293-304(1969) by Parlett and Reinsch.
c    handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
c
c    This routine balances a real matrix and isolates
c    eigenvalues whenever possible.
c
c    Suppose that the principal submatrix in rows LOW through IGH
c    has been balanced, that P(J) denotes the index interchanged
c    with J during the permutation step, and that the elements
c    of the diagonal matrix used are denoted by D(I,J).  Then
c        SCALE(J) = P(J),    for J = 1,...,LOW-1
c                 = D(J,J),      J = LOW,...,IGH
c                 = P(J)         J = IGH+1,...,N.
c    the order in which the interchanges are made is N to IGH+1,
c    then 1 to LOW-1.
c
c    Note that 1 is returned for IGH if IGH is zero formally.
c
c  Parameters:
c
c    Input, integer NM, the row dimension of two-dimensional
c    array parameters as declared in the calling program
c    dimension statement.
c
c    Input, integer N, the order of the matrix.
c
c    Input/output, double precision A(NM,N).  On input, the matrix to be 
c    balanced.  On output, the balanced matrix.
c
c    Output, integer LOW, IGH, such that A(I,J) is equal to zero if
c    (1) I is greater than J and
c    (2) J=1,...,LOW-1 or I=IGH+1,...,N.
c
c    Output, double precision SCALE(N), information determining the
c    permutations and scaling factors used.
c
      integer i,j,k,l,m,n,jj,nm,igh,low,iexc
      double precision a(nm,n),scale(n)
      double precision c,f,g,r,s,b2,radix
      logical noconv

      radix = 16.0d0

      b2 = radix * radix
      k = 1
      l = n
      go to 100
c
c  In-line procedure for row and column exchange.
c
   20 scale(m) = j
      if (j .eq. m) go to 50

      do i = 1, l
         f = a(i,j)
         a(i,j) = a(i,m)
         a(i,m) = f
      end do

      do i = k, n
         f = a(j,i)
         a(j,i) = a(m,i)
         a(m,i) = f
      end do

   50 go to (80,130), iexc
c
c  Search for rows isolating an eigenvalue and push them down.
c
   80 if (l .eq. 1) go to 280
      l = l - 1
c
c  For J = L step -1 until 1 do.
c
  100 do 120 jj = 1, l
         j = l + 1 - jj

         do 110 i = 1, l
            if (i .eq. j) go to 110
            if (a(j,i) .ne. 0.0d0) go to 120
  110    continue

         m = l
         iexc = 1
         go to 20
  120 continue

      go to 140
c
c  Search for columns isolating an eigenvalue and push them left.
c
  130 k = k + 1

  140 do 170 j = k, l

         do 150 i = k, l
            if (i .eq. j) go to 150
            if (a(i,j) .ne. 0.0d0) go to 170
  150    continue

         m = k
         iexc = 2
         go to 20
  170 continue
c
c  Now balance the submatrix in rows K to L.
c
      do i = k, l
        scale(i) = 1.0d0
      end do
c
c  Iterative loop for norm reduction.
c
  190 noconv = .false.

      do 270 i = k, l
         c = 0.0d0
         r = 0.0d0

         do 200 j = k, l
            if (j .eq. i) go to 200
            c = c + dabs(a(j,i))
            r = r + dabs(a(i,j))
  200    continue
c
c  Guard against zero C or R due to underflow.
c
         if (c .eq. 0.0d0 .or. r .eq. 0.0d0) go to 270
         g = r / radix
         f = 1.0d0
         s = c + r
  210    if (c .ge. g) go to 220
         f = f * radix
         c = c * b2
         go to 210
  220    g = r * radix
  230    if (c .lt. g) go to 240
         f = f / radix
         c = c / b2
         go to 230
c
c  Now balance.
c
  240    continue

         if ((c + r) / f .ge. 0.95d0 * s) go to 270
         g = 1.0d0 / f
         scale(i) = scale(i) * f
         noconv = .true.

         do j = k, n
           a(i,j) = a(i,j) * g
         end do

         do j = 1, l
           a(j,i) = a(j,i) * f
         end do

  270 continue

      if (noconv) go to 190

  280 low = k
      igh = l

      return
      end
      subroutine balbak(nm,n,low,igh,scale,m,z)

c*********************************************************************72
c
cc BALBAK determines eigenvectors by undoing the BALANC transformation.
c
c     this routine is a translation of the algol procedure balbak,
c     num. math. 13, 293-304(1969) by parlett and reinsch.
c     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
c
c     this routine forms the eigenvectors of a real general
c     matrix by back transforming those of the corresponding
c     balanced matrix determined by  balanc.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by  balanc.
c
c        scale contains information determining the permutations
c          and scaling factors used by  balanc.
c
c        m is the number of columns of z to be back transformed.
c
c        z contains the real and imaginary parts of the eigen-
c          vectors to be back transformed in its first m columns.
c
c     on output
c
c        z contains the real and imaginary parts of the
c          transformed eigenvectors in its first m columns.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,k,m,n,ii,nm,igh,low
      double precision scale(n),z(nm,m)
      double precision s

      if (m .eq. 0) go to 200
      if (igh .eq. low) go to 120
c
      do 110 i = low, igh
         s = scale(i)
c     .......... left hand eigenvectors are back transformed
c                if the foregoing statement is replaced by
c                s=1.0d0/scale(i). ..........
         do 100 j = 1, m
  100    z(i,j) = z(i,j) * s
c
  110 continue
c     ......... for i=low-1 step -1 until 1,
c               igh+1 step 1 until n do -- ..........
  120 do 140 ii = 1, n
         i = ii
         if (i .ge. low .and. i .le. igh) go to 140
         if (i .lt. low) i = low - ii
         k = scale(i)
         if (k .eq. i) go to 140
c
         do 130 j = 1, m
            s = z(i,j)
            z(i,j) = z(k,j)
            z(k,j) = s
  130    continue
c
  140 continue
c
  200 return
      end
      subroutine bandr(nm,n,mb,a,d,e,e2,matz,z)

c*********************************************************************72
c
cc BANDR reduces a symmetric band matrix to symmetric tridiagonal form.
c
C  REFORMULATED S2 IN LOOP 500 TO AVOID OVERFLOW. (9/29/89 BSG)
c
c     this routine is a translation of the algol procedure bandrd,
c     num. math. 12, 231-241(1968) by schwarz.
c     handbook for auto. comp., vol.ii-linear algebra, 273-283(1971).
c
c     this subroutine reduces a real symmetric band matrix
c     to a symmetric tridiagonal matrix using and optionally
c     accumulating orthogonal similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        mb is the (half) band width of the matrix, defined as the
c          number of adjacent diagonals, including the principal
c          diagonal, required to specify the non-zero portion of the
c          lower triangle of the matrix.
c
c        a contains the lower triangle of the symmetric band input
c          matrix stored as an n by mb array.  its lowest subdiagonal
c          is stored in the last n+1-mb positions of the first column,
c          its next subdiagonal in the last n+2-mb positions of the
c          second column, further subdiagonals similarly, and finally
c          its principal diagonal in the n positions of the last column.
c          contents of storages not part of the matrix are arbitrary.
c
c        matz should be set to .true. if the transformation matrix is
c          to be accumulated, and to .false. otherwise.
c
c     on output
c
c        a has been destroyed, except for its last two columns which
c          contain a copy of the tridiagonal matrix.
c
c        d contains the diagonal elements of the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        e2 contains the squares of the corresponding elements of e.
c          e2 may coincide with e if the squares are not needed.
c
c        z contains the orthogonal transformation matrix produced in
c          the reduction if matz has been set to .true.  otherwise, z
c          is not referenced.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated september 1989.
c
      integer j,k,l,n,r,i1,i2,j1,j2,kr,mb,mr,m1,nm,n2,r1,ugl,maxl,maxr
      double precision a(nm,mb),d(n),e(n),e2(n),z(nm,n)
      double precision g,u,b1,b2,c2,f1,f2,s2,dmin,dminrt
      logical matz

      dmin = 2.0d0**(-64)
      dminrt = 2.0d0**(-32)
c     .......... initialize diagonal scaling matrix ..........
      do 30 j = 1, n
   30 d(j) = 1.0d0
c
      if (.not. matz) go to 60
c
      do 50 j = 1, n
c
         do 40 k = 1, n
   40    z(j,k) = 0.0d0
c
         z(j,j) = 1.0d0
   50 continue
c
   60 m1 = mb - 1
      if (m1 - 1) 900, 800, 70
   70 n2 = n - 2
c
      do 700 k = 1, n2
         maxr = min0(m1,n-k)
c     .......... for r=maxr step -1 until 2 do -- ..........
         do 600 r1 = 2, maxr
            r = maxr + 2 - r1
            kr = k + r
            mr = mb - r
            g = a(kr,mr)
            a(kr-1,1) = a(kr-1,mr+1)
            ugl = k
c
            do 500 j = kr, n, m1
               j1 = j - 1
               j2 = j1 - 1
               if (g .eq. 0.0d0) go to 600
               b1 = a(j1,1) / g
               b2 = b1 * d(j1) / d(j)
               IF (ABS(B1) .GT. 1.0D0) THEN
                  U = 1.0D0 / B1
                  S2 = U / (U + B2)
               ELSE
                  S2 = 1.0D0 / (1.0D0 + B1 * B2)
               ENDIF
c
               if (s2 .ge. 0.5d0 ) go to 450
               b1 = g / a(j1,1)
               b2 = b1 * d(j) / d(j1)
               c2 = 1.0d0 - s2
               d(j1) = c2 * d(j1)
               d(j) = c2 * d(j)
               f1 = 2.0d0 * a(j,m1)
               f2 = b1 * a(j1,mb)
               a(j,m1) = -b2 * (b1 * a(j,m1) - a(j,mb)) - f2 + a(j,m1)
               a(j1,mb) = b2 * (b2 * a(j,mb) + f1) + a(j1,mb)
               a(j,mb) = b1 * (f2 - f1) + a(j,mb)
c
               do 200 l = ugl, j2
                  i2 = mb - j + l
                  u = a(j1,i2+1) + b2 * a(j,i2)
                  a(j,i2) = -b1 * a(j1,i2+1) + a(j,i2)
                  a(j1,i2+1) = u
  200          continue
c
               ugl = j
               a(j1,1) = a(j1,1) + b2 * g
               if (j .eq. n) go to 350
               maxl = min0(m1,n-j1)
c
               do 300 l = 2, maxl
                  i1 = j1 + l
                  i2 = mb - l
                  u = a(i1,i2) + b2 * a(i1,i2+1)
                  a(i1,i2+1) = -b1 * a(i1,i2) + a(i1,i2+1)
                  a(i1,i2) = u
  300          continue
c
               i1 = j + m1
               if (i1 .gt. n) go to 350
               g = b2 * a(i1,1)
  350          if (.not. matz) go to 500
c
               do 400 l = 1, n
                  u = z(l,j1) + b2 * z(l,j)
                  z(l,j) = -b1 * z(l,j1) + z(l,j)
                  z(l,j1) = u
  400          continue
c
               go to 500
c
  450          u = d(j1)
               d(j1) = s2 * d(j)
               d(j) = s2 * u
               f1 = 2.0d0 * a(j,m1)
               f2 = b1 * a(j,mb)
               u = b1 * (f2 - f1) + a(j1,mb)
               a(j,m1) = b2 * (b1 * a(j,m1) - a(j1,mb)) + f2 - a(j,m1)
               a(j1,mb) = b2 * (b2 * a(j1,mb) + f1) + a(j,mb)
               a(j,mb) = u
c
               do 460 l = ugl, j2
                  i2 = mb - j + l
                  u = b2 * a(j1,i2+1) + a(j,i2)
                  a(j,i2) = -a(j1,i2+1) + b1 * a(j,i2)
                  a(j1,i2+1) = u
  460          continue
c
               ugl = j
               a(j1,1) = b2 * a(j1,1) + g
               if (j .eq. n) go to 480
               maxl = min0(m1,n-j1)
c
               do 470 l = 2, maxl
                  i1 = j1 + l
                  i2 = mb - l
                  u = b2 * a(i1,i2) + a(i1,i2+1)
                  a(i1,i2+1) = -a(i1,i2) + b1 * a(i1,i2+1)
                  a(i1,i2) = u
  470          continue
c
               i1 = j + m1
               if (i1 .gt. n) go to 480
               g = a(i1,1)
               a(i1,1) = b1 * a(i1,1)
  480          if (.not. matz) go to 500
c
               do 490 l = 1, n
                  u = b2 * z(l,j1) + z(l,j)
                  z(l,j) = -z(l,j1) + b1 * z(l,j)
                  z(l,j1) = u
  490          continue
c
  500       continue
c
  600    continue
c
         if (mod(k,64) .ne. 0) go to 700
c     .......... rescale to avoid underflow or overflow ..........
         do 650 j = k, n
            if (d(j) .ge. dmin) go to 650
            maxl = max0(1,mb+1-j)
c
            do 610 l = maxl, m1
  610       a(j,l) = dminrt * a(j,l)
c
            if (j .eq. n) go to 630
            maxl = min0(m1,n-j)
c
            do 620 l = 1, maxl
               i1 = j + l
               i2 = mb - l
               a(i1,i2) = dminrt * a(i1,i2)
  620       continue
c
  630       if (.not. matz) go to 645
c
            do 640 l = 1, n
  640       z(l,j) = dminrt * z(l,j)
c
  645       a(j,mb) = dmin * a(j,mb)
            d(j) = d(j) / dmin
  650    continue
c
  700 continue
c     .......... form square root of scaling matrix ..........
  800 do 810 j = 2, n
  810 e(j) = dsqrt(d(j))
c
      if (.not. matz) go to 840
c
      do 830 j = 1, n
c
         do 820 k = 2, n
  820    z(j,k) = e(k) * z(j,k)
c
  830 continue
c
  840 u = 1.0d0
c
      do 850 j = 2, n
         a(j,m1) = u * e(j) * a(j,m1)
         u = e(j)
         e2(j) = a(j,m1) ** 2
         a(j,mb) = d(j) * a(j,mb)
         d(j) = a(j,mb)
         e(j) = a(j,m1)
  850 continue
c
      d(1) = a(1,mb)
      e(1) = 0.0d0
      e2(1) = 0.0d0
      go to 1001
c
  900 do 950 j = 1, n
         d(j) = a(j,mb)
         e(j) = 0.0d0
         e2(j) = 0.0d0
  950 continue
c
 1001 return
      end
      subroutine bandv(nm,n,mbw,a,e21,m,w,z,ierr,nv,rv,rv6)

c*********************************************************************72
c
cc BANDV finds eigenvectors from eigenvalues, for a real symmetric band matrix.
c
c     this subroutine finds those eigenvectors of a real symmetric
c     band matrix corresponding to specified eigenvalues, using inverse
c     iteration.  the subroutine may also be used to solve systems
c     of linear equations with a symmetric or non-symmetric band
c     coefficient matrix.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        mbw is the number of columns of the array a used to store the
c          band matrix.  if the matrix is symmetric, mbw is its (half)
c          band width, denoted mb and defined as the number of adjacent
c          diagonals, including the principal diagonal, required to
c          specify the non-zero portion of the lower triangle of the
c          matrix.  if the subroutine is being used to solve systems
c          of linear equations and the coefficient matrix is not
c          symmetric, it must however have the same number of adjacent
c          diagonals above the main diagonal as below, and in this
c          case, mbw=2*mb-1.
c
c        a contains the lower triangle of the symmetric band input
c          matrix stored as an n by mb array.  its lowest subdiagonal
c          is stored in the last n+1-mb positions of the first column,
c          its next subdiagonal in the last n+2-mb positions of the
c          second column, further subdiagonals similarly, and finally
c          its principal diagonal in the n positions of column mb.
c          if the subroutine is being used to solve systems of linear
c          equations and the coefficient matrix is not symmetric, a is
c          n by 2*mb-1 instead with lower triangle as above and with
c          its first superdiagonal stored in the first n-1 positions of
c          column mb+1, its second superdiagonal in the first n-2
c          positions of column mb+2, further superdiagonals similarly,
c          and finally its highest superdiagonal in the first n+1-mb
c          positions of the last column.
c          contents of storages not part of the matrix are arbitrary.
c
c        e21 specifies the ordering of the eigenvalues and contains
c            0.0d0 if the eigenvalues are in ascending order, or
c            2.0d0 if the eigenvalues are in descending order.
c          if the subroutine is being used to solve systems of linear
c          equations, e21 should be set to 1.0d0 if the coefficient
c          matrix is symmetric and to -1.0d0 if not.
c
c        m is the number of specified eigenvalues or the number of
c          systems of linear equations.
c
c        w contains the m eigenvalues in ascending or descending order.
c          if the subroutine is being used to solve systems of linear
c          equations (a-w(r)*i)*x(r)=b(r), where i is the identity
c          matrix, w(r) should be set accordingly, for r=1,2,...,m.
c
c        z contains the constant matrix columns (b(r),r=1,2,...,m), if
c          the subroutine is used to solve systems of linear equations.
c
c        nv must be set to the dimension of the array parameter rv
c          as declared in the calling program dimension statement.
c
c     on output
c
c        a and w are unaltered.
c
c        z contains the associated set of orthogonal eigenvectors.
c          any vector which fails to converge is set to zero.  if the
c          subroutine is used to solve systems of linear equations,
c          z contains the solution matrix columns (x(r),r=1,2,...,m).
c
c        ierr is set to
c          zero       for normal return,
c          -r         if the eigenvector corresponding to the r-th
c                     eigenvalue fails to converge, or if the r-th
c                     system of linear equations is nearly singular.
c
c        rv and rv6 are temporary storage arrays.  note that rv is
c          of dimension at least n*(2*mb-1).  if the subroutine
c          is being used to solve systems of linear equations, the
c          determinant (up to sign) of a-w(m)*i is available, upon
c          return, as the product of the first n elements of rv.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,k,m,n,r,ii,ij,jj,kj,mb,m1,nm,nv,ij1,its,kj1,mbw,m21,
     x        ierr,maxj,maxk,group
      double precision a(nm,mbw),w(m),z(nm,m),rv(nv),rv6(n)
      double precision u,v,uk,xu,x0,x1,e21,eps2,eps3,eps4,norm,order,
     x       epslon,pythag

      ierr = 0
      if (m .eq. 0) go to 1001
      mb = mbw
      if (e21 .lt. 0.0d0) mb = (mbw + 1) / 2
      m1 = mb - 1
      m21 = m1 + mb
      order = 1.0d0 - dabs(e21)
c     .......... find vectors by inverse iteration ..........
      do 920 r = 1, m
         its = 1
         x1 = w(r)
         if (r .ne. 1) go to 100
c     .......... compute norm of matrix ..........
         norm = 0.0d0
c
         do 60 j = 1, mb
            jj = mb + 1 - j
            kj = jj + m1
            ij = 1
            v = 0.0d0
c
            do 40 i = jj, n
               v = v + dabs(a(i,j))
               if (e21 .ge. 0.0d0) go to 40
               v = v + dabs(a(ij,kj))
               ij = ij + 1
   40       continue
c
            norm = dmax1(norm,v)
   60    continue
c
         if (e21 .lt. 0.0d0) norm = 0.5d0 * norm
c     .......... eps2 is the criterion for grouping,
c                eps3 replaces zero pivots and equal
c                roots are modified by eps3,
c                eps4 is taken very small to avoid overflow ..........
         if (norm .eq. 0.0d0) norm = 1.0d0
         eps2 = 1.0d-3 * norm * dabs(order)
         eps3 = epslon(norm)
         uk = n
         uk = dsqrt(uk)
         eps4 = uk * eps3
   80    group = 0
         go to 120
c     .......... look for close or coincident roots ..........
  100    if (dabs(x1-x0) .ge. eps2) go to 80
         group = group + 1
         if (order * (x1 - x0) .le. 0.0d0) x1 = x0 + order * eps3
c     .......... expand matrix, subtract eigenvalue,
c                and initialize vector ..........
  120    do 200 i = 1, n
            ij = i + min0(0,i-m1) * n
            kj = ij + mb * n
            ij1 = kj + m1 * n
            if (m1 .eq. 0) go to 180
c
            do 150 j = 1, m1
               if (ij .gt. m1) go to 125
               if (ij .gt. 0) go to 130
               rv(ij1) = 0.0d0
               ij1 = ij1 + n
               go to 130
  125          rv(ij) = a(i,j)
  130          ij = ij + n
               ii = i + j
               if (ii .gt. n) go to 150
               jj = mb - j
               if (e21 .ge. 0.0d0) go to 140
               ii = i
               jj = mb + j
  140          rv(kj) = a(ii,jj)
               kj = kj + n
  150       continue
c
  180       rv(ij) = a(i,mb) - x1
            rv6(i) = eps4
            if (order .eq. 0.0d0) rv6(i) = z(i,r)
  200    continue
c
         if (m1 .eq. 0) go to 600
c     .......... elimination with interchanges ..........
         do 580 i = 1, n
            ii = i + 1
            maxk = min0(i+m1-1,n)
            maxj = min0(n-i,m21-2) * n
c
            do 360 k = i, maxk
               kj1 = k
               j = kj1 + n
               jj = j + maxj
c
               do 340 kj = j, jj, n
                  rv(kj1) = rv(kj)
                  kj1 = kj
  340          continue
c
               rv(kj1) = 0.0d0
  360       continue
c
            if (i .eq. n) go to 580
            u = 0.0d0
            maxk = min0(i+m1,n)
            maxj = min0(n-ii,m21-2) * n
c
            do 450 j = i, maxk
               if (dabs(rv(j)) .lt. dabs(u)) go to 450
               u = rv(j)
               k = j
  450       continue
c
            j = i + n
            jj = j + maxj
            if (k .eq. i) go to 520
            kj = k
c
            do 500 ij = i, jj, n
               v = rv(ij)
               rv(ij) = rv(kj)
               rv(kj) = v
               kj = kj + n
  500       continue
c
            if (order .ne. 0.0d0) go to 520
            v = rv6(i)
            rv6(i) = rv6(k)
            rv6(k) = v
  520       if (u .eq. 0.0d0) go to 580
c
            do 560 k = ii, maxk
               v = rv(k) / u
               kj = k
c
               do 540 ij = j, jj, n
                  kj = kj + n
                  rv(kj) = rv(kj) - v * rv(ij)
  540          continue
c
               if (order .eq. 0.0d0) rv6(k) = rv6(k) - v * rv6(i)
  560       continue
c
  580    continue
c     .......... back substitution
c                for i=n step -1 until 1 do -- ..........
  600    do 630 ii = 1, n
            i = n + 1 - ii
            maxj = min0(ii,m21)
            if (maxj .eq. 1) go to 620
            ij1 = i
            j = ij1 + n
            jj = j + (maxj - 2) * n
c
            do 610 ij = j, jj, n
               ij1 = ij1 + 1
               rv6(i) = rv6(i) - rv(ij) * rv6(ij1)
  610       continue
c
  620       v = rv(i)
            if (dabs(v) .ge. eps3) go to 625
c     .......... set error -- nearly singular linear system ..........
            if (order .eq. 0.0d0) ierr = -r
            v = dsign(eps3,v)
  625       rv6(i) = rv6(i) / v
  630    continue
c
         xu = 1.0d0
         if (order .eq. 0.0d0) go to 870
c     .......... orthogonalize with respect to previous
c                members of group ..........
         if (group .eq. 0) go to 700
c
         do 680 jj = 1, group
            j = r - group - 1 + jj
            xu = 0.0d0
c
            do 640 i = 1, n
  640       xu = xu + rv6(i) * z(i,j)
c
            do 660 i = 1, n
  660       rv6(i) = rv6(i) - xu * z(i,j)
c
  680    continue
c
  700    norm = 0.0d0
c
         do 720 i = 1, n
  720    norm = norm + dabs(rv6(i))
c
         if (norm .ge. 0.1d0) go to 840
c     .......... in-line procedure for choosing
c                a new starting vector ..........
         if (its .ge. n) go to 830
         its = its + 1
         xu = eps4 / (uk + 1.0d0)
         rv6(1) = eps4
c
         do 760 i = 2, n
  760    rv6(i) = xu
c
         rv6(its) = rv6(its) - eps4 * uk
         go to 600
c     .......... set error -- non-converged eigenvector ..........
  830    ierr = -r
         xu = 0.0d0
         go to 870
c     .......... normalize so that sum of squares is
c                1 and expand to full order ..........
  840    u = 0.0d0
c
         do 860 i = 1, n
  860    u = pythag(u,rv6(i))
c
         xu = 1.0d0 / u
c
  870    do 900 i = 1, n
  900    z(i,r) = rv6(i) * xu
c
         x0 = x1
  920 continue
c
 1001 return
      end
      subroutine bisect(n,eps1,d,e,e2,lb,ub,mm,m,w,ind,ierr,rv4,rv5)

c*********************************************************************72
c
cc BISECT computes some eigenvalues of a real symmetric tridiagonal matrix.
c
c     this subroutine is a translation of the bisection technique
c     in the algol procedure tristurm by peters and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 418-439(1971).
c
c     this subroutine finds those eigenvalues of a tridiagonal
c     symmetric matrix which lie in a specified interval,
c     using bisection.
c
c     on input
c
c        n is the order of the matrix.
c
c        eps1 is an absolute error tolerance for the computed
c          eigenvalues.  if the input eps1 is non-positive,
c          it is reset for each submatrix to a default value,
c          namely, minus the product of the relative machine
c          precision and the 1-norm of the submatrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c        e2 contains the squares of the corresponding elements of e.
c          e2(1) is arbitrary.
c
c        lb and ub define the interval to be searched for eigenvalues.
c          if lb is not less than ub, no eigenvalues will be found.
c
c        mm should be set to an upper bound for the number of
c          eigenvalues in the interval.  warning. if more than
c          mm eigenvalues are determined to lie in the interval,
c          an error return is made with no eigenvalues found.
c
c     on output
c
c        eps1 is unaltered unless it has been reset to its
c          (last) default value.
c
c        d and e are unaltered.
c
c        elements of e2, corresponding to elements of e regarded
c          as negligible, have been replaced by zero causing the
c          matrix to split into a direct sum of submatrices.
c          e2(1) is also set to zero.
c
c        m is the number of eigenvalues determined to lie in (lb,ub).
c
c        w contains the m eigenvalues in ascending order.
c
c        ind contains in its first m positions the submatrix indices
c          associated with the corresponding eigenvalues in w --
c          1 for eigenvalues belonging to the first submatrix from
c          the top, 2 for those belonging to the second submatrix, etc..
c
c        ierr is set to
c          zero       for normal return,
c          3*n+1      if m exceeds mm.
c
c        rv4 and rv5 are temporary storage arrays.
c
c     the algol procedure sturmcnt contained in tristurm
c     appears in bisect in-line.
c
c     note that subroutine tql1 or imtql1 is generally faster than
c     bisect, if more than n/4 eigenvalues are to be found.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,k,l,m,n,p,q,r,s,ii,mm,m1,m2,tag,ierr,isturm
      double precision d(n),e(n),e2(n),w(mm),rv4(n),rv5(n)
      double precision u,v,lb,t1,t2,ub,xu,x0,x1,eps1,tst1,tst2,epslon
      integer ind(mm)

      ierr = 0
      tag = 0
      t1 = lb
      t2 = ub
c     .......... look for small sub-diagonal entries ..........
      do 40 i = 1, n
         if (i .eq. 1) go to 20
         tst1 = dabs(d(i)) + dabs(d(i-1))
         tst2 = tst1 + dabs(e(i))
         if (tst2 .gt. tst1) go to 40
   20    e2(i) = 0.0d0
   40 continue
c     .......... determine the number of eigenvalues
c                in the interval ..........
      p = 1
      q = n
      x1 = ub
      isturm = 1
      go to 320
   60 m = s
      x1 = lb
      isturm = 2
      go to 320
   80 m = m - s
      if (m .gt. mm) go to 980
      q = 0
      r = 0
c     .......... establish and process next submatrix, refining
c                interval by the gerschgorin bounds ..........
  100 if (r .eq. m) go to 1001
      tag = tag + 1
      p = q + 1
      xu = d(p)
      x0 = d(p)
      u = 0.0d0
c
      do 120 q = p, n
         x1 = u
         u = 0.0d0
         v = 0.0d0
         if (q .eq. n) go to 110
         u = dabs(e(q+1))
         v = e2(q+1)
  110    xu = dmin1(d(q)-(x1+u),xu)
         x0 = dmax1(d(q)+(x1+u),x0)
         if (v .eq. 0.0d0) go to 140
  120 continue
c
  140 x1 = epslon(dmax1(dabs(xu),dabs(x0)))
      if (eps1 .le. 0.0d0) eps1 = -x1
      if (p .ne. q) go to 180
c     .......... check for isolated root within interval ..........
      if (t1 .gt. d(p) .or. d(p) .ge. t2) go to 940
      m1 = p
      m2 = p
      rv5(p) = d(p)
      go to 900
  180 x1 = x1 * (q - p + 1)
      lb = dmax1(t1,xu-x1)
      ub = dmin1(t2,x0+x1)
      x1 = lb
      isturm = 3
      go to 320
  200 m1 = s + 1
      x1 = ub
      isturm = 4
      go to 320
  220 m2 = s
      if (m1 .gt. m2) go to 940
c     .......... find roots by bisection ..........
      x0 = ub
      isturm = 5
c
      do 240 i = m1, m2
         rv5(i) = ub
         rv4(i) = lb
  240 continue
c     .......... loop for k-th eigenvalue
c                for k=m2 step -1 until m1 do --
c                (-do- not used to legalize -computed go to-) ..........
      k = m2
  250    xu = lb
c     .......... for i=k step -1 until m1 do -- ..........
         do 260 ii = m1, k
            i = m1 + k - ii
            if (xu .ge. rv4(i)) go to 260
            xu = rv4(i)
            go to 280
  260    continue
c
  280    if (x0 .gt. rv5(k)) x0 = rv5(k)
c     .......... next bisection step ..........
  300    x1 = (xu + x0) * 0.5d0
         if ((x0 - xu) .le. dabs(eps1)) go to 420
         tst1 = 2.0d0 * (dabs(xu) + dabs(x0))
         tst2 = tst1 + (x0 - xu)
         if (tst2 .eq. tst1) go to 420
c     .......... in-line procedure for sturm sequence ..........
  320    s = p - 1
         u = 1.0d0
c
         do 340 i = p, q
            if (u .ne. 0.0d0) go to 325
            v = dabs(e(i)) / epslon(1.0d0)
            if (e2(i) .eq. 0.0d0) v = 0.0d0
            go to 330
  325       v = e2(i) / u
  330       u = d(i) - x1 - v
            if (u .lt. 0.0d0) s = s + 1
  340    continue
c
         go to (60,80,200,220,360), isturm
c     .......... refine intervals ..........
  360    if (s .ge. k) go to 400
         xu = x1
         if (s .ge. m1) go to 380
         rv4(m1) = x1
         go to 300
  380    rv4(s+1) = x1
         if (rv5(s) .gt. x1) rv5(s) = x1
         go to 300
  400    x0 = x1
         go to 300
c     .......... k-th eigenvalue found ..........
  420    rv5(k) = x1
      k = k - 1
      if (k .ge. m1) go to 250
c     .......... order eigenvalues tagged with their
c                submatrix associations ..........
  900 s = r
      r = r + m2 - m1 + 1
      j = 1
      k = m1
c
      do 920 l = 1, r
         if (j .gt. s) go to 910
         if (k .gt. m2) go to 940
         if (rv5(k) .ge. w(l)) go to 915
c
         do 905 ii = j, s
            i = l + s - ii
            w(i+1) = w(i)
            ind(i+1) = ind(i)
  905    continue
c
  910    w(l) = rv5(k)
         ind(l) = tag
         k = k + 1
         go to 920
  915    j = j + 1
  920 continue
c
  940 if (q .lt. n) go to 100
      go to 1001
c     .......... set error -- underestimate of number of
c                eigenvalues in interval ..........
  980 ierr = 3 * n + 1
 1001 lb = t1
      ub = t2
      return
      end
      subroutine bqr(nm,n,mb,a,t,r,ierr,nv,rv)

c*********************************************************************72
c
cc BQR finds the smallest eigenvalue of a real symmetric band matrix.
c
c     this subroutine is a translation of the algol procedure bqr,
c     num. math. 16, 85-92(1970) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol ii-linear algebra, 266-272(1971).
c
c     this subroutine finds the eigenvalue of smallest (usually)
c     magnitude of a real symmetric band matrix using the
c     qr algorithm with shifts of origin.  consecutive calls
c     can be made to find further eigenvalues.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        mb is the (half) band width of the matrix, defined as the
c          number of adjacent diagonals, including the principal
c          diagonal, required to specify the non-zero portion of the
c          lower triangle of the matrix.
c
c        a contains the lower triangle of the symmetric band input
c          matrix stored as an n by mb array.  its lowest subdiagonal
c          is stored in the last n+1-mb positions of the first column,
c          its next subdiagonal in the last n+2-mb positions of the
c          second column, further subdiagonals similarly, and finally
c          its principal diagonal in the n positions of the last column.
c          contents of storages not part of the matrix are arbitrary.
c          on a subsequent call, its output contents from the previous
c          call should be passed.
c
c        t specifies the shift (of eigenvalues) applied to the diagonal
c          of a in forming the input matrix. what is actually determined
c          is the eigenvalue of a+ti (i is the identity matrix) nearest
c          to t.  on a subsequent call, the output value of t from the
c          previous call should be passed if the next nearest eigenvalue
c          is sought.
c
c        r should be specified as zero on the first call, and as its
c          output value from the previous call on a subsequent call.
c          it is used to determine when the last row and column of
c          the transformed band matrix can be regarded as negligible.
c
c        nv must be set to the dimension of the array parameter rv
c          as declared in the calling program dimension statement.
c
c     on output
c
c        a contains the transformed band matrix.  the matrix a+ti
c          derived from the output parameters is similar to the
c          input a+ti to within rounding errors.  its last row and
c          column are null (if ierr is zero).
c
c        t contains the computed eigenvalue of a+ti (if ierr is zero).
c
c        r contains the maximum of its input value and the norm of the
c          last column of the input matrix a.
c
c        ierr is set to
c          zero       for normal return,
c          n          if the eigenvalue has not been
c                     determined after 30 iterations.
c
c        rv is a temporary storage array of dimension at least
c          (2*mb**2+4*mb-3).  the first (3*mb-2) locations correspond
c          to the algol array b, the next (2*mb-1) locations correspond
c          to the algol array h, and the final (2*mb**2-mb) locations
c          correspond to the mb by (2*mb-1) algol array u.
c
c     note. for a subsequent call, n should be replaced by n-1, but
c     mb should not be altered even when it exceeds the current n.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,k,l,m,n,ii,ik,jk,jm,kj,kk,km,ll,mb,mk,mn,mz,
     x        m1,m2,m3,m4,ni,nm,nv,its,kj1,m21,m31,ierr,imult
      double precision a(nm,mb),rv(nv)
      double precision f,g,q,r,s,t,tst1,tst2,scale,pythag

      ierr = 0
      m1 = min0(mb,n)
      m = m1 - 1
      m2 = m + m
      m21 = m2 + 1
      m3 = m21 + m
      m31 = m3 + 1
      m4 = m31 + m2
      mn = m + n
      mz = mb - m1
      its = 0
c     .......... test for convergence ..........
   40 g = a(n,mb)
      if (m .eq. 0) go to 360
      f = 0.0d0
c
      do 50 k = 1, m
         mk = k + mz
         f = f + dabs(a(n,mk))
   50 continue
c
      if (its .eq. 0 .and. f .gt. r) r = f
      tst1 = r
      tst2 = tst1 + f
      if (tst2 .le. tst1) go to 360
      if (its .eq. 30) go to 1000
      its = its + 1
c     .......... form shift from bottom 2 by 2 minor ..........
      if (f .gt. 0.25d0 * r .and. its .lt. 5) go to 90
      f = a(n,mb-1)
      if (f .eq. 0.0d0) go to 70
      q = (a(n-1,mb) - g) / (2.0d0 * f)
      s = pythag(q,1.0d0)
      g = g - f / (q + dsign(s,q))
   70 t = t + g
c
      do 80 i = 1, n
   80 a(i,mb) = a(i,mb) - g
c
   90 do 100 k = m31, m4
  100 rv(k) = 0.0d0
c
      do 350 ii = 1, mn
         i = ii - m
         ni = n - ii
         if (ni .lt. 0) go to 230
c     .......... form column of shifted matrix a-g*i ..........
         l = max0(1,2-i)
c
         do 110 k = 1, m3
  110    rv(k) = 0.0d0
c
         do 120 k = l, m1
            km = k + m
            mk = k + mz
            rv(km) = a(ii,mk)
  120    continue
c
         ll = min0(m,ni)
         if (ll .eq. 0) go to 135
c
         do 130 k = 1, ll
            km = k + m21
            ik = ii + k
            mk = mb - k
            rv(km) = a(ik,mk)
  130    continue
c     .......... pre-multiply with householder reflections ..........
  135    ll = m2
         imult = 0
c     .......... multiplication procedure ..........
  140    kj = m4 - m1
c
         do 170 j = 1, ll
            kj = kj + m1
            jm = j + m3
            if (rv(jm) .eq. 0.0d0) go to 170
            f = 0.0d0
c
            do 150 k = 1, m1
               kj = kj + 1
               jk = j + k - 1
               f = f + rv(kj) * rv(jk)
  150       continue
c
            f = f / rv(jm)
            kj = kj - m1
c
            do 160 k = 1, m1
               kj = kj + 1
               jk = j + k - 1
               rv(jk) = rv(jk) - rv(kj) * f
  160       continue
c
            kj = kj - m1
  170    continue
c
         if (imult .ne. 0) go to 280
c     .......... householder reflection ..........
         f = rv(m21)
         s = 0.0d0
         rv(m4) = 0.0d0
         scale = 0.0d0
c
         do 180 k = m21, m3
  180    scale = scale + dabs(rv(k))
c
         if (scale .eq. 0.0d0) go to 210
c
         do 190 k = m21, m3
  190    s = s + (rv(k)/scale)**2
c
         s = scale * scale * s
         g = -dsign(dsqrt(s),f)
         rv(m21) = g
         rv(m4) = s - f * g
         kj = m4 + m2 * m1 + 1
         rv(kj) = f - g
c
         do 200 k = 2, m1
            kj = kj + 1
            km = k + m2
            rv(kj) = rv(km)
  200    continue
c     .......... save column of triangular factor r ..........
  210    do 220 k = l, m1
            km = k + m
            mk = k + mz
            a(ii,mk) = rv(km)
  220    continue
c
  230    l = max0(1,m1+1-i)
         if (i .le. 0) go to 300
c     .......... perform additional steps ..........
         do 240 k = 1, m21
  240    rv(k) = 0.0d0
c
         ll = min0(m1,ni+m1)
c     .......... get row of triangular factor r ..........
         do 250 kk = 1, ll
            k = kk - 1
            km = k + m1
            ik = i + k
            mk = mb - k
            rv(km) = a(ik,mk)
  250    continue
c     .......... post-multiply with householder reflections ..........
         ll = m1
         imult = 1
         go to 140
c     .......... store column of new a matrix ..........
  280    do 290 k = l, m1
            mk = k + mz
            a(i,mk) = rv(k)
  290    continue
c     .......... update householder reflections ..........
  300    if (l .gt. 1) l = l - 1
         kj1 = m4 + l * m1
c
         do 320 j = l, m2
            jm = j + m3
            rv(jm) = rv(jm+1)
c
            do 320 k = 1, m1
               kj1 = kj1 + 1
               kj = kj1 - m1
               rv(kj) = rv(kj1)
  320    continue
c
  350 continue
c
      go to 40
c     .......... convergence ..........
  360 t = t + g
c
      do 380 i = 1, n
  380 a(i,mb) = a(i,mb) - g
c
      do 400 k = 1, m1
         mk = k + mz
         a(n,mk) = 0.0d0
  400 continue
c
      go to 1001
c     .......... set error -- no convergence to
c                eigenvalue after 30 iterations ..........
 1000 ierr = n
 1001 return
      end
      subroutine cbabk2(nm,n,low,igh,scale,m,zr,zi)

c*********************************************************************72
c
cc CBABK2 finds eigenvectors by undoing the CBAL transformation.
c
c     this subroutine is a translation of the algol procedure
c     cbabk2, which is a complex version of balbak,
c     num. math. 13, 293-304(1969) by parlett and reinsch.
c     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
c
c     this subroutine forms the eigenvectors of a complex general
c     matrix by back transforming those of the corresponding
c     balanced matrix determined by  cbal.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by  cbal.
c
c        scale contains information determining the permutations
c          and scaling factors used by  cbal.
c
c        m is the number of eigenvectors to be back transformed.
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the eigenvectors to be
c          back transformed in their first m columns.
c
c     on output
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the transformed eigenvectors
c          in their first m columns.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,k,m,n,ii,nm,igh,low
      double precision scale(n),zr(nm,m),zi(nm,m)
      double precision s

      if (m .eq. 0) go to 200
      if (igh .eq. low) go to 120
c
      do 110 i = low, igh
         s = scale(i)
c     .......... left hand eigenvectors are back transformed
c                if the foregoing statement is replaced by
c                s=1.0d0/scale(i). ..........
         do 100 j = 1, m
            zr(i,j) = zr(i,j) * s
            zi(i,j) = zi(i,j) * s
  100    continue
c
  110 continue
c     .......... for i=low-1 step -1 until 1,
c                igh+1 step 1 until n do -- ..........
  120 do 140 ii = 1, n
         i = ii
         if (i .ge. low .and. i .le. igh) go to 140
         if (i .lt. low) i = low - ii
         k = scale(i)
         if (k .eq. i) go to 140
c
         do 130 j = 1, m
            s = zr(i,j)
            zr(i,j) = zr(k,j)
            zr(k,j) = s
            s = zi(i,j)
            zi(i,j) = zi(k,j)
            zi(k,j) = s
  130    continue
c
  140 continue
c
  200 return
      end
      subroutine cbal(nm,n,ar,ai,low,igh,scale)

c*********************************************************************72
c
cc CBAL balances a complex matrix before eigenvalue calculations.
c
c     this subroutine is a translation of the algol procedure
c     cbalance, which is a complex version of balance,
c     num. math. 13, 293-304(1969) by parlett and reinsch.
c     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
c
c     this subroutine balances a complex matrix and isolates
c     eigenvalues whenever possible.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the complex matrix to be balanced.
c
c     on output
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the balanced matrix.
c
c        low and igh are two integers such that ar(i,j) and ai(i,j)
c          are equal to zero if
c           (1) i is greater than j and
c           (2) j=1,...,low-1 or i=igh+1,...,n.
c
c        scale contains information determining the
c           permutations and scaling factors used.
c
c     suppose that the principal submatrix in rows low through igh
c     has been balanced, that p(j) denotes the index interchanged
c     with j during the permutation step, and that the elements
c     of the diagonal matrix used are denoted by d(i,j).  then
c        scale(j) = p(j),    for j = 1,...,low-1
c                 = d(j,j)       j = low,...,igh
c                 = p(j)         j = igh+1,...,n.
c     the order in which the interchanges are made is n to igh+1,
c     then 1 to low-1.
c
c     note that 1 is returned for igh if igh is zero formally.
c
c     the algol procedure exc contained in cbalance appears in
c     cbal  in line.  (note that the algol roles of identifiers
c     k,l have been reversed.)
c
c     arithmetic is real throughout.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,k,l,m,n,jj,nm,igh,low,iexc
      double precision ar(nm,n),ai(nm,n),scale(n)
      double precision c,f,g,r,s,b2,radix
      logical noconv

      radix = 16.0d0
c
      b2 = radix * radix
      k = 1
      l = n
      go to 100
c     .......... in-line procedure for row and
c                column exchange ..........
   20 scale(m) = j
      if (j .eq. m) go to 50
c
      do 30 i = 1, l
         f = ar(i,j)
         ar(i,j) = ar(i,m)
         ar(i,m) = f
         f = ai(i,j)
         ai(i,j) = ai(i,m)
         ai(i,m) = f
   30 continue
c
      do 40 i = k, n
         f = ar(j,i)
         ar(j,i) = ar(m,i)
         ar(m,i) = f
         f = ai(j,i)
         ai(j,i) = ai(m,i)
         ai(m,i) = f
   40 continue
c
   50 go to (80,130), iexc
c     .......... search for rows isolating an eigenvalue
c                and push them down ..........
   80 if (l .eq. 1) go to 280
      l = l - 1
c     .......... for j=l step -1 until 1 do -- ..........
  100 do 120 jj = 1, l
         j = l + 1 - jj
c
         do 110 i = 1, l
            if (i .eq. j) go to 110
            if (ar(j,i) .ne. 0.0d0 .or. ai(j,i) .ne. 0.0d0) go to 120
  110    continue
c
         m = l
         iexc = 1
         go to 20
  120 continue
c
      go to 140
c     .......... search for columns isolating an eigenvalue
c                and push them left ..........
  130 k = k + 1
c
  140 do 170 j = k, l
c
         do 150 i = k, l
            if (i .eq. j) go to 150
            if (ar(i,j) .ne. 0.0d0 .or. ai(i,j) .ne. 0.0d0) go to 170
  150    continue
c
         m = k
         iexc = 2
         go to 20
  170 continue
c     .......... now balance the submatrix in rows k to l ..........
      do 180 i = k, l
  180 scale(i) = 1.0d0
c     .......... iterative loop for norm reduction ..........
  190 noconv = .false.
c
      do 270 i = k, l
         c = 0.0d0
         r = 0.0d0
c
         do 200 j = k, l
            if (j .eq. i) go to 200
            c = c + dabs(ar(j,i)) + dabs(ai(j,i))
            r = r + dabs(ar(i,j)) + dabs(ai(i,j))
  200    continue
c     .......... guard against zero c or r due to underflow ..........
         if (c .eq. 0.0d0 .or. r .eq. 0.0d0) go to 270
         g = r / radix
         f = 1.0d0
         s = c + r
  210    if (c .ge. g) go to 220
         f = f * radix
         c = c * b2
         go to 210
  220    g = r * radix
  230    if (c .lt. g) go to 240
         f = f / radix
         c = c / b2
         go to 230
c     .......... now balance ..........
  240    if ((c + r) / f .ge. 0.95d0 * s) go to 270
         g = 1.0d0 / f
         scale(i) = scale(i) * f
         noconv = .true.
c
         do 250 j = k, n
            ar(i,j) = ar(i,j) * g
            ai(i,j) = ai(i,j) * g
  250    continue
c
         do 260 j = 1, l
            ar(j,i) = ar(j,i) * f
            ai(j,i) = ai(j,i) * f
  260    continue
c
  270 continue
c
      if (noconv) go to 190
c
  280 low = k
      igh = l
      return
      end
      subroutine cdiv(ar,ai,br,bi,cr,ci)

c*********************************************************************72
c
cc CDIV emulates complex division, using real arithmetic.
c
c     complex division, (cr,ci) = (ar,ai)/(br,bi)
c
      double precision ar,ai,br,bi,cr,ci
      double precision s,ars,ais,brs,bis

      s = dabs(br) + dabs(bi)
      ars = ar/s
      ais = ai/s
      brs = br/s
      bis = bi/s
      s = brs**2 + bis**2
      cr = (ars*brs + ais*bis)/s
      ci = (ais*brs - ars*bis)/s
      return
      end
      subroutine cg(nm,n,ar,ai,wr,wi,matz,zr,zi,fv1,fv2,fv3,ierr)

c*********************************************************************72
c
cc CG gets eigenvalues and eigenvectors of a complex general matrix.
c
c     this subroutine calls the recommended sequence of
c     routines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a complex general matrix.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrix  a=(ar,ai).
c
c        ar  and  ai  contain the real and imaginary parts,
c        respectively, of the complex general matrix.
c
c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        wr  and  wi  contain the real and imaginary parts,
c        respectively, of the eigenvalues.
c
c        zr  and  zi  contain the real and imaginary parts,
c        respectively, of the eigenvectors if matz is not zero.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for comqr
c           and comqr2.  the normal completion code is zero.
c
c        fv1, fv2, and  fv3  are temporary storage arrays.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c

      integer n,nm,is1,is2,ierr,matz
      double precision ar(nm,n),ai(nm,n),wr(n),wi(n),zr(nm,n),zi(nm,n),
     x       fv1(n),fv2(n),fv3(n)

      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 call  cbal(nm,n,ar,ai,is1,is2,fv1)
      call  corth(nm,n,is1,is2,ar,ai,fv2,fv3)
      if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  comqr(nm,n,is1,is2,ar,ai,wr,wi,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 call  comqr2(nm,n,is1,is2,fv2,fv3,ar,ai,wr,wi,zr,zi,ierr)
      if (ierr .ne. 0) go to 50
      call  cbabk2(nm,n,is1,is2,fv1,n,zr,zi)
   50 return
      end
      subroutine ch(nm,n,ar,ai,w,matz,zr,zi,fv1,fv2,fm1,ierr)

c*********************************************************************72
c
cc CH gets eigenvalues and eigenvectors of a complex Hermitian matrix.
c
c     this subroutine calls the recommended sequence of
c     routines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a complex hermitian matrix.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrix  a=(ar,ai).
c
c        ar  and  ai  contain the real and imaginary parts,
c        respectively, of the complex hermitian matrix.
c
c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        w  contains the eigenvalues in ascending order.
c
c        zr  and  zi  contain the real and imaginary parts,
c        respectively, of the eigenvectors if matz is not zero.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for tqlrat
c           and tql2.  the normal completion code is zero.
c
c        fv1, fv2, and  fm1  are temporary storage arrays.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,n,nm,ierr,matz
      double precision ar(nm,n),ai(nm,n),w(n),zr(nm,n),zi(nm,n),
     x       fv1(n),fv2(n),fm1(2,n)

      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 call  htridi(nm,n,ar,ai,w,fv1,fv2,fm1)
      if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  tqlrat(n,w,fv2,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 do 40 i = 1, n
c
         do 30 j = 1, n
            zr(j,i) = 0.0d0
   30    continue
c
         zr(i,i) = 1.0d0
   40 continue
c
      call  tql2(nm,n,w,fv1,zr,ierr)
      if (ierr .ne. 0) go to 50
      call  htribk(nm,n,ar,ai,fm1,n,zr,zi)
   50 return
      end
      subroutine cinvit(nm,n,ar,ai,wr,wi,select,mm,m,zr,zi,
     x                  ierr,rm1,rm2,rv1,rv2)

c*********************************************************************72
c
cc CINVIT gets eigenvectors from eigenvalues, for a complex Hessenberg matrix.
c
c     this subroutine is a translation of the algol procedure cx invit
c     by peters and wilkinson.
c     handbook for auto. comp. vol.ii-linear algebra, 418-439(1971).
c
c     this subroutine finds those eigenvectors of a complex upper
c     hessenberg matrix corresponding to specified eigenvalues,
c     using inverse iteration.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the hessenberg matrix.
c
c        wr and wi contain the real and imaginary parts, respectively,
c          of the eigenvalues of the matrix.  the eigenvalues must be
c          stored in a manner identical to that of subroutine  comlr,
c          which recognizes possible splitting of the matrix.
c
c        select specifies the eigenvectors to be found.  the
c          eigenvector corresponding to the j-th eigenvalue is
c          specified by setting select(j) to .true..
c
c        mm should be set to an upper bound for the number of
c          eigenvectors to be found.
c
c     on output
c
c        ar, ai, wi, and select are unaltered.
c
c        wr may have been altered since close eigenvalues are perturbed
c          slightly in searching for independent eigenvectors.
c
c        m is the number of eigenvectors actually found.
c
c        zr and zi contain the real and imaginary parts, respectively,
c          of the eigenvectors.  the eigenvectors are normalized
c          so that the component of largest magnitude is 1.
c          any vector which fails the acceptance test is set to zero.
c
c        ierr is set to
c          zero       for normal return,
c          -(2*n+1)   if more than mm eigenvectors have been specified,
c          -k         if the iteration corresponding to the k-th
c                     value fails,
c          -(n+k)     if both error situations occur.
c
c        rm1, rm2, rv1, and rv2 are temporary storage arrays.
c
c     the algol procedure guessvec appears in cinvit in line.
c
c     calls cdiv for complex division.
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,k,m,n,s,ii,mm,mp,nm,uk,ip1,its,km1,ierr
      double precision ar(nm,n),ai(nm,n),wr(n),wi(n),zr(nm,mm),
     x       zi(nm,mm),rm1(n,n),rm2(n,n),rv1(n),rv2(n)
      double precision x,y,eps3,norm,normv,epslon,growto,ilambd,pythag,
     x       rlambd,ukroot
      logical select(n)

      ierr = 0
      uk = 0
      s = 1
c
      do 980 k = 1, n
         if (.not. select(k)) go to 980
         if (s .gt. mm) go to 1000
         if (uk .ge. k) go to 200
c     .......... check for possible splitting ..........
         do 120 uk = k, n
            if (uk .eq. n) go to 140
            if (ar(uk+1,uk) .eq. 0.0d0 .and. ai(uk+1,uk) .eq. 0.0d0)
     x         go to 140
  120    continue
c     .......... compute infinity norm of leading uk by uk
c                (hessenberg) matrix ..........
  140    norm = 0.0d0
         mp = 1
c
         do 180 i = 1, uk
            x = 0.0d0
c
            do 160 j = mp, uk
  160       x = x + pythag(ar(i,j),ai(i,j))
c
            if (x .gt. norm) norm = x
            mp = i
  180    continue
c     .......... eps3 replaces zero pivot in decomposition
c                and close roots are modified by eps3 ..........
         if (norm .eq. 0.0d0) norm = 1.0d0
         eps3 = epslon(norm)
c     .......... growto is the criterion for growth ..........
         ukroot = uk
         ukroot = dsqrt(ukroot)
         growto = 0.1d0 / ukroot
  200    rlambd = wr(k)
         ilambd = wi(k)
         if (k .eq. 1) go to 280
         km1 = k - 1
         go to 240
c     .......... perturb eigenvalue if it is close
c                to any previous eigenvalue ..........
  220    rlambd = rlambd + eps3
c     .......... for i=k-1 step -1 until 1 do -- ..........
  240    do 260 ii = 1, km1
            i = k - ii
            if (select(i) .and. dabs(wr(i)-rlambd) .lt. eps3 .and.
     x         dabs(wi(i)-ilambd) .lt. eps3) go to 220
  260    continue
c
         wr(k) = rlambd
c     .......... form upper hessenberg (ar,ai)-(rlambd,ilambd)*i
c                and initial complex vector ..........
  280    mp = 1
c
         do 320 i = 1, uk
c
            do 300 j = mp, uk
               rm1(i,j) = ar(i,j)
               rm2(i,j) = ai(i,j)
  300       continue
c
            rm1(i,i) = rm1(i,i) - rlambd
            rm2(i,i) = rm2(i,i) - ilambd
            mp = i
            rv1(i) = eps3
  320    continue
c     .......... triangular decomposition with interchanges,
c                replacing zero pivots by eps3 ..........
         if (uk .eq. 1) go to 420
c
         do 400 i = 2, uk
            mp = i - 1
            if (pythag(rm1(i,mp),rm2(i,mp)) .le.
     x          pythag(rm1(mp,mp),rm2(mp,mp))) go to 360
c
            do 340 j = mp, uk
               y = rm1(i,j)
               rm1(i,j) = rm1(mp,j)
               rm1(mp,j) = y
               y = rm2(i,j)
               rm2(i,j) = rm2(mp,j)
               rm2(mp,j) = y
  340       continue
c
  360       if (rm1(mp,mp) .eq. 0.0d0 .and. rm2(mp,mp) .eq. 0.0d0)
     x         rm1(mp,mp) = eps3
            call cdiv(rm1(i,mp),rm2(i,mp),rm1(mp,mp),rm2(mp,mp),x,y)
            if (x .eq. 0.0d0 .and. y .eq. 0.0d0) go to 400
c
            do 380 j = i, uk
               rm1(i,j) = rm1(i,j) - x * rm1(mp,j) + y * rm2(mp,j)
               rm2(i,j) = rm2(i,j) - x * rm2(mp,j) - y * rm1(mp,j)
  380       continue
c
  400    continue
c
  420    if (rm1(uk,uk) .eq. 0.0d0 .and. rm2(uk,uk) .eq. 0.0d0)
     x      rm1(uk,uk) = eps3
         its = 0
c     .......... back substitution
c                for i=uk step -1 until 1 do -- ..........
  660    do 720 ii = 1, uk
            i = uk + 1 - ii
            x = rv1(i)
            y = 0.0d0
            if (i .eq. uk) go to 700
            ip1 = i + 1
c
            do 680 j = ip1, uk
               x = x - rm1(i,j) * rv1(j) + rm2(i,j) * rv2(j)
               y = y - rm1(i,j) * rv2(j) - rm2(i,j) * rv1(j)
  680       continue
c
  700       call cdiv(x,y,rm1(i,i),rm2(i,i),rv1(i),rv2(i))
  720    continue
c     .......... acceptance test for eigenvector
c                and normalization ..........
         its = its + 1
         norm = 0.0d0
         normv = 0.0d0
c
         do 780 i = 1, uk
            x = pythag(rv1(i),rv2(i))
            if (normv .ge. x) go to 760
            normv = x
            j = i
  760       norm = norm + x
  780    continue
c
         if (norm .lt. growto) go to 840
c     .......... accept vector ..........
         x = rv1(j)
         y = rv2(j)
c
         do 820 i = 1, uk
            call cdiv(rv1(i),rv2(i),x,y,zr(i,s),zi(i,s))
  820    continue
c
         if (uk .eq. n) go to 940
         j = uk + 1
         go to 900
c     .......... in-line procedure for choosing
c                a new starting vector ..........
  840    if (its .ge. uk) go to 880
         x = ukroot
         y = eps3 / (x + 1.0d0)
         rv1(1) = eps3
c
         do 860 i = 2, uk
  860    rv1(i) = y
c
         j = uk - its + 1
         rv1(j) = rv1(j) - eps3 * x
         go to 660
c     .......... set error -- unaccepted eigenvector ..........
  880    j = 1
         ierr = -k
c     .......... set remaining vector components to zero ..........
  900    do 920 i = j, n
            zr(i,s) = 0.0d0
            zi(i,s) = 0.0d0
  920    continue
c
  940    s = s + 1
  980 continue
c
      go to 1001
c     .......... set error -- underestimate of eigenvector
c                space required ..........
 1000 if (ierr .ne. 0) ierr = ierr - n
      if (ierr .eq. 0) ierr = -(2 * n + 1)
 1001 m = s - 1
      return
      end
      subroutine combak(nm,low,igh,ar,ai,int,m,zr,zi)

c*********************************************************************72
c
cc COMBAK determines eigenvectors by undoing the COMHES transformation.
c
c     this subroutine is a translation of the algol procedure combak,
c     num. math. 12, 349-368(1968) by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
c
c     this subroutine forms the eigenvectors of a complex general
c     matrix by back transforming those of the corresponding
c     upper hessenberg matrix determined by  comhes.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        low and igh are integers determined by the balancing
c        routine  cbal.  if  cbal  has not been used,
c          set low=1 and igh equal to the order of the matrix.
c
c        ar and ai contain the multipliers which were used in the
c          reduction by  comhes  in their lower triangles
c          below the subdiagonal.
c
c        int contains information on the rows and columns
c          interchanged in the reduction by  comhes.
c          only elements low through igh are used.
c
c        m is the number of eigenvectors to be back transformed.
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the eigenvectors to be
c          back transformed in their first m columns.
c
c     on output
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the transformed eigenvectors
c          in their first m columns.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c

      integer i,j,m,la,mm,mp,nm,igh,kp1,low,mp1
      double precision ar(nm,igh),ai(nm,igh),zr(nm,m),zi(nm,m)
      double precision xr,xi
      integer int(igh)

      if (m .eq. 0) go to 200
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
c     .......... for mp=igh-1 step -1 until low+1 do -- ..........
      do 140 mm = kp1, la
         mp = low + igh - mm
         mp1 = mp + 1
c
         do 110 i = mp1, igh
            xr = ar(i,mp-1)
            xi = ai(i,mp-1)
            if (xr .eq. 0.0d0 .and. xi .eq. 0.0d0) go to 110
c
            do 100 j = 1, m
               zr(i,j) = zr(i,j) + xr * zr(mp,j) - xi * zi(mp,j)
               zi(i,j) = zi(i,j) + xr * zi(mp,j) + xi * zr(mp,j)
  100       continue
c
  110    continue
c
         i = int(mp)
         if (i .eq. mp) go to 140
c
         do 130 j = 1, m
            xr = zr(i,j)
            zr(i,j) = zr(mp,j)
            zr(mp,j) = xr
            xi = zi(i,j)
            zi(i,j) = zi(mp,j)
            zi(mp,j) = xi
  130    continue
c
  140 continue
c
  200 return
      end
      subroutine comhes(nm,n,low,igh,ar,ai,int)

c*********************************************************************72
c
cc COMHES transforms a complex general matrix to upper Hessenberg form.
c
c     this subroutine is a translation of the algol procedure comhes,
c     num. math. 12, 349-368(1968) by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
c
c     given a complex general matrix, this subroutine
c     reduces a submatrix situated in rows and columns
c     low through igh to upper hessenberg form by
c     stabilized elementary similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c        routine  cbal.  if  cbal  has not been used,
c          set low=1, igh=n.
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the complex input matrix.
c
c     on output
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the hessenberg matrix.  the
c          multipliers which were used in the reduction
c          are stored in the remaining triangles under the
c          hessenberg matrix.
c
c        int contains information on the rows and columns
c          interchanged in the reduction.
c          only elements low through igh are used.
c
c     calls cdiv for complex division.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c

      integer i,j,m,n,la,nm,igh,kp1,low,mm1,mp1
      double precision ar(nm,n),ai(nm,n)
      double precision xr,xi,yr,yi
      integer int(igh)

      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
c
      do 180 m = kp1, la
         mm1 = m - 1
         xr = 0.0d0
         xi = 0.0d0
         i = m
c
         do 100 j = m, igh
            if (dabs(ar(j,mm1)) + dabs(ai(j,mm1))
     x         .le. dabs(xr) + dabs(xi)) go to 100
            xr = ar(j,mm1)
            xi = ai(j,mm1)
            i = j
  100    continue
c
         int(m) = i
         if (i .eq. m) go to 130
c     .......... interchange rows and columns of ar and ai ..........
         do 110 j = mm1, n
            yr = ar(i,j)
            ar(i,j) = ar(m,j)
            ar(m,j) = yr
            yi = ai(i,j)
            ai(i,j) = ai(m,j)
            ai(m,j) = yi
  110    continue
c
         do 120 j = 1, igh
            yr = ar(j,i)
            ar(j,i) = ar(j,m)
            ar(j,m) = yr
            yi = ai(j,i)
            ai(j,i) = ai(j,m)
            ai(j,m) = yi
  120    continue
c     .......... end interchange ..........
  130    if (xr .eq. 0.0d0 .and. xi .eq. 0.0d0) go to 180
         mp1 = m + 1
c
         do 160 i = mp1, igh
            yr = ar(i,mm1)
            yi = ai(i,mm1)
            if (yr .eq. 0.0d0 .and. yi .eq. 0.0d0) go to 160
            call cdiv(yr,yi,xr,xi,yr,yi)
            ar(i,mm1) = yr
            ai(i,mm1) = yi
c
            do 140 j = m, n
               ar(i,j) = ar(i,j) - yr * ar(m,j) + yi * ai(m,j)
               ai(i,j) = ai(i,j) - yr * ai(m,j) - yi * ar(m,j)
  140       continue
c
            do 150 j = 1, igh
               ar(j,m) = ar(j,m) + yr * ar(j,i) - yi * ai(j,i)
               ai(j,m) = ai(j,m) + yr * ai(j,i) + yi * ar(j,i)
  150       continue
c
  160    continue
c
  180 continue
c
  200 return
      end
      subroutine comlr(nm,n,low,igh,hr,hi,wr,wi,ierr)

c*********************************************************************72
c
cc COMLR gets all eigenvalues of a complex upper Hessenberg matrix.
c
c     this subroutine is a translation of the algol procedure comlr,
c     num. math. 12, 369-376(1968) by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 396-403(1971).
c
c     this subroutine finds the eigenvalues of a complex
c     upper hessenberg matrix by the modified lr method.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c        routine  cbal.  if  cbal  has not been used,
c          set low=1, igh=n.
c
c        hr and hi contain the real and imaginary parts,
c          respectively, of the complex upper hessenberg matrix.
c          their lower triangles below the subdiagonal contain the
c          multipliers which were used in the reduction by  comhes,
c          if performed.
c
c     on output
c
c        the upper hessenberg portions of hr and hi have been
c          destroyed.  therefore, they must be saved before
c          calling  comlr  if subsequent calculation of
c          eigenvectors is to be performed.
c
c        wr and wi contain the real and imaginary parts,
c          respectively, of the eigenvalues.  if an error
c          exit is made, the eigenvalues should be correct
c          for indices ierr+1,...,n.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the limit of 30*n iterations is exhausted
c                     while the j-th eigenvalue is being sought.
c
c     calls cdiv for complex division.
c     calls csroot for complex square root.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,l,m,n,en,ll,mm,nm,igh,im1,itn,its,low,mp1,enm1,ierr
      double precision hr(nm,n),hi(nm,n),wr(n),wi(n)
      double precision si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,tst1,tst2

      ierr = 0
c     .......... store roots isolated by cbal ..........
      do 200 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 200
         wr(i) = hr(i,i)
         wi(i) = hi(i,i)
  200 continue
c
      en = igh
      tr = 0.0d0
      ti = 0.0d0
      itn = 30*n
c     .......... search for next eigenvalue ..........
  220 if (en .lt. low) go to 1001
      its = 0
      enm1 = en - 1
c     .......... look for single small sub-diagonal element
c                for l=en step -1 until low d0 -- ..........
  240 do 260 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 300
         tst1 = dabs(hr(l-1,l-1)) + dabs(hi(l-1,l-1))
     x            + dabs(hr(l,l)) + dabs(hi(l,l))
         tst2 = tst1 + dabs(hr(l,l-1)) + dabs(hi(l,l-1))
         if (tst2 .eq. tst1) go to 300
  260 continue
c     .......... form shift ..........
  300 if (l .eq. en) go to 660
      if (itn .eq. 0) go to 1000
      if (its .eq. 10 .or. its .eq. 20) go to 320
      sr = hr(en,en)
      si = hi(en,en)
      xr = hr(enm1,en) * hr(en,enm1) - hi(enm1,en) * hi(en,enm1)
      xi = hr(enm1,en) * hi(en,enm1) + hi(enm1,en) * hr(en,enm1)
      if (xr .eq. 0.0d0 .and. xi .eq. 0.0d0) go to 340
      yr = (hr(enm1,enm1) - sr) / 2.0d0
      yi = (hi(enm1,enm1) - si) / 2.0d0
      call csroot(yr**2-yi**2+xr,2.0d0*yr*yi+xi,zzr,zzi)
      if (yr * zzr + yi * zzi .ge. 0.0d0) go to 310
      zzr = -zzr
      zzi = -zzi
  310 call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)
      sr = sr - xr
      si = si - xi
      go to 340
c     .......... form exceptional shift ..........
  320 sr = dabs(hr(en,enm1)) + dabs(hr(enm1,en-2))
      si = dabs(hi(en,enm1)) + dabs(hi(enm1,en-2))
c
  340 do 360 i = low, en
         hr(i,i) = hr(i,i) - sr
         hi(i,i) = hi(i,i) - si
  360 continue
c
      tr = tr + sr
      ti = ti + si
      its = its + 1
      itn = itn - 1
c     .......... look for two consecutive small
c                sub-diagonal elements ..........
      xr = dabs(hr(enm1,enm1)) + dabs(hi(enm1,enm1))
      yr = dabs(hr(en,enm1)) + dabs(hi(en,enm1))
      zzr = dabs(hr(en,en)) + dabs(hi(en,en))
c     .......... for m=en-1 step -1 until l do -- ..........
      do 380 mm = l, enm1
         m = enm1 + l - mm
         if (m .eq. l) go to 420
         yi = yr
         yr = dabs(hr(m,m-1)) + dabs(hi(m,m-1))
         xi = zzr
         zzr = xr
         xr = dabs(hr(m-1,m-1)) + dabs(hi(m-1,m-1))
         tst1 = zzr / yi * (zzr + xr + xi)
         tst2 = tst1 + yr
         if (tst2 .eq. tst1) go to 420
  380 continue
c     .......... triangular decomposition h=l*r ..........
  420 mp1 = m + 1
c
      do 520 i = mp1, en
         im1 = i - 1
         xr = hr(im1,im1)
         xi = hi(im1,im1)
         yr = hr(i,im1)
         yi = hi(i,im1)
         if (dabs(xr) + dabs(xi) .ge. dabs(yr) + dabs(yi)) go to 460
c     .......... interchange rows of hr and hi ..........
         do 440 j = im1, en
            zzr = hr(im1,j)
            hr(im1,j) = hr(i,j)
            hr(i,j) = zzr
            zzi = hi(im1,j)
            hi(im1,j) = hi(i,j)
            hi(i,j) = zzi
  440    continue
c
         call cdiv(xr,xi,yr,yi,zzr,zzi)
         wr(i) = 1.0d0
         go to 480
  460    call cdiv(yr,yi,xr,xi,zzr,zzi)
         wr(i) = -1.0d0
  480    hr(i,im1) = zzr
         hi(i,im1) = zzi
c
         do 500 j = i, en
            hr(i,j) = hr(i,j) - zzr * hr(im1,j) + zzi * hi(im1,j)
            hi(i,j) = hi(i,j) - zzr * hi(im1,j) - zzi * hr(im1,j)
  500    continue
c
  520 continue
c     .......... composition r*l=h ..........
      do 640 j = mp1, en
         xr = hr(j,j-1)
         xi = hi(j,j-1)
         hr(j,j-1) = 0.0d0
         hi(j,j-1) = 0.0d0
c     .......... interchange columns of hr and hi,
c                if necessary ..........
         if (wr(j) .le. 0.0d0) go to 580
c
         do 540 i = l, j
            zzr = hr(i,j-1)
            hr(i,j-1) = hr(i,j)
            hr(i,j) = zzr
            zzi = hi(i,j-1)
            hi(i,j-1) = hi(i,j)
            hi(i,j) = zzi
  540    continue
c
  580    do 600 i = l, j
            hr(i,j-1) = hr(i,j-1) + xr * hr(i,j) - xi * hi(i,j)
            hi(i,j-1) = hi(i,j-1) + xr * hi(i,j) + xi * hr(i,j)
  600    continue
c
  640 continue
c
      go to 240
c     .......... a root found ..........
  660 wr(en) = hr(en,en) + tr
      wi(en) = hi(en,en) + ti
      en = enm1
      go to 220
c     .......... set error -- all eigenvalues have not
c                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end
      subroutine comlr2(nm,n,low,igh,int,hr,hi,wr,wi,zr,zi,ierr)

c*********************************************************************72
c
cc COMLR2 gets eigenvalues/vectors of a complex upper Hessenberg matrix.
c
C  MESHED overflow control WITH vectors of isolated roots (10/19/89 BSG)
C  MESHED overflow control WITH triangular multiply (10/30/89 BSG)
c
c
c     this subroutine is a translation of the algol procedure comlr2,
c     num. math. 16, 181-204(1970) by peters and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a complex upper hessenberg matrix by the modified lr
c     method.  the eigenvectors of a complex general matrix
c     can also be found if  comhes  has been used to reduce
c     this general matrix to hessenberg form.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c        routine  cbal.  if  cbal  has not been used,
c          set low=1, igh=n.
c
c        int contains information on the rows and columns interchanged
c          in the reduction by  comhes, if performed.  only elements
c          low through igh are used.  if the eigenvectors of the hessen-
c          berg matrix are desired, set int(j)=j for these elements.
c
c        hr and hi contain the real and imaginary parts,
c          respectively, of the complex upper hessenberg matrix.
c          their lower triangles below the subdiagonal contain the
c          multipliers which were used in the reduction by  comhes,
c          if performed.  if the eigenvectors of the hessenberg
c          matrix are desired, these elements must be set to zero.
c
c     on output
c
c        the upper hessenberg portions of hr and hi have been
c          destroyed, but the location hr(1,1) contains the norm
c          of the triangularized matrix.
c
c        wr and wi contain the real and imaginary parts,
c          respectively, of the eigenvalues.  if an error
c          exit is made, the eigenvalues should be correct
c          for indices ierr+1,...,n.
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the eigenvectors.  the eigenvectors
c          are unnormalized.  if an error exit is made, none of
c          the eigenvectors has been found.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the limit of 30*n iterations is exhausted
c                     while the j-th eigenvalue is being sought.
c
c
c     calls cdiv for complex division.
c     calls csroot for complex square root.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated october 1989.
c
      integer i,j,k,l,m,n,en,ii,jj,ll,mm,nm,nn,igh,im1,ip1,
     x        itn,its,low,mp1,enm1,iend,ierr
      double precision hr(nm,n),hi(nm,n),wr(n),wi(n),zr(nm,n),zi(nm,n)
      double precision si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,norm,tst1,tst2
      integer int(igh)

      ierr = 0
c     .......... initialize eigenvector matrix ..........
      do 100 i = 1, n
c
         do 100 j = 1, n
            zr(i,j) = 0.0d0
            zi(i,j) = 0.0d0
            if (i .eq. j) zr(i,j) = 1.0d0
  100 continue
c     .......... form the matrix of accumulated transformations
c                from the information left by comhes ..........
      iend = igh - low - 1
      if (iend .le. 0) go to 180
c     .......... for i=igh-1 step -1 until low+1 do -- ..........
      do 160 ii = 1, iend
         i = igh - ii
         ip1 = i + 1
c
         do 120 k = ip1, igh
            zr(k,i) = hr(k,i-1)
            zi(k,i) = hi(k,i-1)
  120    continue
c
         j = int(i)
         if (i .eq. j) go to 160
c
         do 140 k = i, igh
            zr(i,k) = zr(j,k)
            zi(i,k) = zi(j,k)
            zr(j,k) = 0.0d0
            zi(j,k) = 0.0d0
  140    continue
c
         zr(j,i) = 1.0d0
  160 continue
c     .......... store roots isolated by cbal ..........
  180 do 200 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 200
         wr(i) = hr(i,i)
         wi(i) = hi(i,i)
  200 continue
c
      en = igh
      tr = 0.0d0
      ti = 0.0d0
      itn = 30*n
c     .......... search for next eigenvalue ..........
  220 if (en .lt. low) go to 680
      its = 0
      enm1 = en - 1
c     .......... look for single small sub-diagonal element
c                for l=en step -1 until low do -- ..........
  240 do 260 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 300
         tst1 = dabs(hr(l-1,l-1)) + dabs(hi(l-1,l-1))
     x            + dabs(hr(l,l)) + dabs(hi(l,l))
         tst2 = tst1 + dabs(hr(l,l-1)) + dabs(hi(l,l-1))
         if (tst2 .eq. tst1) go to 300
  260 continue
c     .......... form shift ..........
  300 if (l .eq. en) go to 660
      if (itn .eq. 0) go to 1000
      if (its .eq. 10 .or. its .eq. 20) go to 320
      sr = hr(en,en)
      si = hi(en,en)
      xr = hr(enm1,en) * hr(en,enm1) - hi(enm1,en) * hi(en,enm1)
      xi = hr(enm1,en) * hi(en,enm1) + hi(enm1,en) * hr(en,enm1)
      if (xr .eq. 0.0d0 .and. xi .eq. 0.0d0) go to 340
      yr = (hr(enm1,enm1) - sr) / 2.0d0
      yi = (hi(enm1,enm1) - si) / 2.0d0
      call csroot(yr**2-yi**2+xr,2.0d0*yr*yi+xi,zzr,zzi)
      if (yr * zzr + yi * zzi .ge. 0.0d0) go to 310
      zzr = -zzr
      zzi = -zzi
  310 call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)
      sr = sr - xr
      si = si - xi
      go to 340
c     .......... form exceptional shift ..........
  320 sr = dabs(hr(en,enm1)) + dabs(hr(enm1,en-2))
      si = dabs(hi(en,enm1)) + dabs(hi(enm1,en-2))
c
  340 do 360 i = low, en
         hr(i,i) = hr(i,i) - sr
         hi(i,i) = hi(i,i) - si
  360 continue
c
      tr = tr + sr
      ti = ti + si
      its = its + 1
      itn = itn - 1
c     .......... look for two consecutive small
c                sub-diagonal elements ..........
      xr = dabs(hr(enm1,enm1)) + dabs(hi(enm1,enm1))
      yr = dabs(hr(en,enm1)) + dabs(hi(en,enm1))
      zzr = dabs(hr(en,en)) + dabs(hi(en,en))
c     .......... for m=en-1 step -1 until l do -- ..........
      do 380 mm = l, enm1
         m = enm1 + l - mm
         if (m .eq. l) go to 420
         yi = yr
         yr = dabs(hr(m,m-1)) + dabs(hi(m,m-1))
         xi = zzr
         zzr = xr
         xr = dabs(hr(m-1,m-1)) + dabs(hi(m-1,m-1))
         tst1 = zzr / yi * (zzr + xr + xi)
         tst2 = tst1 + yr
         if (tst2 .eq. tst1) go to 420
  380 continue
c     .......... triangular decomposition h=l*r ..........
  420 mp1 = m + 1
c
      do 520 i = mp1, en
         im1 = i - 1
         xr = hr(im1,im1)
         xi = hi(im1,im1)
         yr = hr(i,im1)
         yi = hi(i,im1)
         if (dabs(xr) + dabs(xi) .ge. dabs(yr) + dabs(yi)) go to 460
c     .......... interchange rows of hr and hi ..........
         do 440 j = im1, n
            zzr = hr(im1,j)
            hr(im1,j) = hr(i,j)
            hr(i,j) = zzr
            zzi = hi(im1,j)
            hi(im1,j) = hi(i,j)
            hi(i,j) = zzi
  440    continue
c
         call cdiv(xr,xi,yr,yi,zzr,zzi)
         wr(i) = 1.0d0
         go to 480
  460    call cdiv(yr,yi,xr,xi,zzr,zzi)
         wr(i) = -1.0d0
  480    hr(i,im1) = zzr
         hi(i,im1) = zzi
c
         do 500 j = i, n
            hr(i,j) = hr(i,j) - zzr * hr(im1,j) + zzi * hi(im1,j)
            hi(i,j) = hi(i,j) - zzr * hi(im1,j) - zzi * hr(im1,j)
  500    continue
c
  520 continue
c     .......... composition r*l=h ..........
      do 640 j = mp1, en
         xr = hr(j,j-1)
         xi = hi(j,j-1)
         hr(j,j-1) = 0.0d0
         hi(j,j-1) = 0.0d0
c     .......... interchange columns of hr, hi, zr, and zi,
c                if necessary ..........
         if (wr(j) .le. 0.0d0) go to 580
c
         do 540 i = 1, j
            zzr = hr(i,j-1)
            hr(i,j-1) = hr(i,j)
            hr(i,j) = zzr
            zzi = hi(i,j-1)
            hi(i,j-1) = hi(i,j)
            hi(i,j) = zzi
  540    continue
c
         do 560 i = low, igh
            zzr = zr(i,j-1)
            zr(i,j-1) = zr(i,j)
            zr(i,j) = zzr
            zzi = zi(i,j-1)
            zi(i,j-1) = zi(i,j)
            zi(i,j) = zzi
  560    continue
c
  580    do 600 i = 1, j
            hr(i,j-1) = hr(i,j-1) + xr * hr(i,j) - xi * hi(i,j)
            hi(i,j-1) = hi(i,j-1) + xr * hi(i,j) + xi * hr(i,j)
  600    continue
c     .......... accumulate transformations ..........
         do 620 i = low, igh
            zr(i,j-1) = zr(i,j-1) + xr * zr(i,j) - xi * zi(i,j)
            zi(i,j-1) = zi(i,j-1) + xr * zi(i,j) + xi * zr(i,j)
  620    continue
c
  640 continue
c
      go to 240
c     .......... a root found ..........
  660 hr(en,en) = hr(en,en) + tr
      wr(en) = hr(en,en)
      hi(en,en) = hi(en,en) + ti
      wi(en) = hi(en,en)
      en = enm1
      go to 220
c     .......... all roots found.  backsubstitute to find
c                vectors of upper triangular form ..........
  680 norm = 0.0d0
c
      do 720 i = 1, n
c
         do 720 j = i, n
            tr = dabs(hr(i,j)) + dabs(hi(i,j))
            if (tr .gt. norm) norm = tr
  720 continue
c
      hr(1,1) = norm
      if (n .eq. 1 .or. norm .eq. 0.0d0) go to 1001
c     .......... for en=n step -1 until 2 do -- ..........
      do 800 nn = 2, n
         en = n + 2 - nn
         xr = wr(en)
         xi = wi(en)
         hr(en,en) = 1.0d0
         hi(en,en) = 0.0d0
         enm1 = en - 1
c     .......... for i=en-1 step -1 until 1 do -- ..........
         do 780 ii = 1, enm1
            i = en - ii
            zzr = 0.0d0
            zzi = 0.0d0
            ip1 = i + 1
c
            do 740 j = ip1, en
               zzr = zzr + hr(i,j) * hr(j,en) - hi(i,j) * hi(j,en)
               zzi = zzi + hr(i,j) * hi(j,en) + hi(i,j) * hr(j,en)
  740       continue
c
            yr = xr - wr(i)
            yi = xi - wi(i)
            if (yr .ne. 0.0d0 .or. yi .ne. 0.0d0) go to 765
               tst1 = norm
               yr = tst1
  760          yr = 0.01d0 * yr
               tst2 = norm + yr
               if (tst2 .gt. tst1) go to 760
  765       continue
            call cdiv(zzr,zzi,yr,yi,hr(i,en),hi(i,en))
c     .......... overflow control ..........
            tr = dabs(hr(i,en)) + dabs(hi(i,en))
            if (tr .eq. 0.0d0) go to 780
            tst1 = tr
            tst2 = tst1 + 1.0d0/tst1
            if (tst2 .gt. tst1) go to 780
            do 770 j = i, en
               hr(j,en) = hr(j,en)/tr
               hi(j,en) = hi(j,en)/tr
  770       continue
c
  780    continue
c
  800 continue
c     .......... end backsubstitution ..........
c     .......... vectors of isolated roots ..........
      do  840 i = 1, N
         if (i .ge. low .and. i .le. igh) go to 840
c
         do 820 j = I, n
            zr(i,j) = hr(i,j)
            zi(i,j) = hi(i,j)
  820    continue
c
  840 continue
c     .......... multiply by transformation matrix to give
c                vectors of original full matrix.
c                for j=n step -1 until low do -- ..........
      do 880 jj = low, N
         j = n + low - jj
         m = min0(j,igh)
c
         do 880 i = low, igh
            zzr = 0.0d0
            zzi = 0.0d0
c
            do 860 k = low, m
               zzr = zzr + zr(i,k) * hr(k,j) - zi(i,k) * hi(k,j)
               zzi = zzi + zr(i,k) * hi(k,j) + zi(i,k) * hr(k,j)
  860       continue
c
            zr(i,j) = zzr
            zi(i,j) = zzi
  880 continue
c
      go to 1001
c     .......... set error -- all eigenvalues have not
c                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end
      subroutine comqr(nm,n,low,igh,hr,hi,wr,wi,ierr)

c*********************************************************************72
c
cc COMQR gets eigenvalues of a complex upper Hessenberg matrix.
c
c     this subroutine is a translation of a unitary analogue of the
c     algol procedure  comlr, num. math. 12, 369-376(1968) by martin
c     and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 396-403(1971).
c     the unitary analogue substitutes the qr algorithm of francis
c     (comp. jour. 4, 332-345(1962)) for the lr algorithm.
c
c     this subroutine finds the eigenvalues of a complex
c     upper hessenberg matrix by the qr method.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c        routine  cbal.  if  cbal  has not been used,
c          set low=1, igh=n.
c
c        hr and hi contain the real and imaginary parts,
c          respectively, of the complex upper hessenberg matrix.
c          their lower triangles below the subdiagonal contain
c          information about the unitary transformations used in
c          the reduction by  corth, if performed.
c
c     on output
c
c        the upper hessenberg portions of hr and hi have been
c          destroyed.  therefore, they must be saved before
c          calling  comqr  if subsequent calculation of
c          eigenvectors is to be performed.
c
c        wr and wi contain the real and imaginary parts,
c          respectively, of the eigenvalues.  if an error
c          exit is made, the eigenvalues should be correct
c          for indices ierr+1,...,n.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the limit of 30*n iterations is exhausted
c                     while the j-th eigenvalue is being sought.
c
c     calls cdiv for complex division.
c     calls csroot for complex square root.
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,l,n,en,ll,nm,igh,itn,its,low,lp1,enm1,ierr
      double precision hr(nm,n),hi(nm,n),wr(n),wi(n)
      double precision si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,norm,tst1,tst2,
     x       pythag

      ierr = 0
      if (low .eq. igh) go to 180
c     .......... create real subdiagonal elements ..........
      l = low + 1
c
      do 170 i = l, igh
         ll = min0(i+1,igh)
         if (hi(i,i-1) .eq. 0.0d0) go to 170
         norm = pythag(hr(i,i-1),hi(i,i-1))
         yr = hr(i,i-1) / norm
         yi = hi(i,i-1) / norm
         hr(i,i-1) = norm
         hi(i,i-1) = 0.0d0
c
         do 155 j = i, igh
            si = yr * hi(i,j) - yi * hr(i,j)
            hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
            hi(i,j) = si
  155    continue
c
         do 160 j = low, ll
            si = yr * hi(j,i) + yi * hr(j,i)
            hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
            hi(j,i) = si
  160    continue
c
  170 continue
c     .......... store roots isolated by cbal ..........
  180 do 200 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 200
         wr(i) = hr(i,i)
         wi(i) = hi(i,i)
  200 continue
c
      en = igh
      tr = 0.0d0
      ti = 0.0d0
      itn = 30*n
c     .......... search for next eigenvalue ..........
  220 if (en .lt. low) go to 1001
      its = 0
      enm1 = en - 1
c     .......... look for single small sub-diagonal element
c                for l=en step -1 until low d0 -- ..........
  240 do 260 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 300
         tst1 = dabs(hr(l-1,l-1)) + dabs(hi(l-1,l-1))
     x            + dabs(hr(l,l)) + dabs(hi(l,l))
         tst2 = tst1 + dabs(hr(l,l-1))
         if (tst2 .eq. tst1) go to 300
  260 continue
c     .......... form shift ..........
  300 if (l .eq. en) go to 660
      if (itn .eq. 0) go to 1000
      if (its .eq. 10 .or. its .eq. 20) go to 320
      sr = hr(en,en)
      si = hi(en,en)
      xr = hr(enm1,en) * hr(en,enm1)
      xi = hi(enm1,en) * hr(en,enm1)
      if (xr .eq. 0.0d0 .and. xi .eq. 0.0d0) go to 340
      yr = (hr(enm1,enm1) - sr) / 2.0d0
      yi = (hi(enm1,enm1) - si) / 2.0d0
      call csroot(yr**2-yi**2+xr,2.0d0*yr*yi+xi,zzr,zzi)
      if (yr * zzr + yi * zzi .ge. 0.0d0) go to 310
      zzr = -zzr
      zzi = -zzi
  310 call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)
      sr = sr - xr
      si = si - xi
      go to 340
c     .......... form exceptional shift ..........
  320 sr = dabs(hr(en,enm1)) + dabs(hr(enm1,en-2))
      si = 0.0d0
c
  340 do 360 i = low, en
         hr(i,i) = hr(i,i) - sr
         hi(i,i) = hi(i,i) - si
  360 continue
c
      tr = tr + sr
      ti = ti + si
      its = its + 1
      itn = itn - 1
c     .......... reduce to triangle (rows) ..........
      lp1 = l + 1
c
      do 500 i = lp1, en
         sr = hr(i,i-1)
         hr(i,i-1) = 0.0d0
         norm = pythag(pythag(hr(i-1,i-1),hi(i-1,i-1)),sr)
         xr = hr(i-1,i-1) / norm
         wr(i-1) = xr
         xi = hi(i-1,i-1) / norm
         wi(i-1) = xi
         hr(i-1,i-1) = norm
         hi(i-1,i-1) = 0.0d0
         hi(i,i-1) = sr / norm
c
         do 490 j = i, en
            yr = hr(i-1,j)
            yi = hi(i-1,j)
            zzr = hr(i,j)
            zzi = hi(i,j)
            hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
            hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
            hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
            hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
  490    continue
c
  500 continue
c
      si = hi(en,en)
      if (si .eq. 0.0d0) go to 540
      norm = pythag(hr(en,en),si)
      sr = hr(en,en) / norm
      si = si / norm
      hr(en,en) = norm
      hi(en,en) = 0.0d0
c     .......... inverse operation (columns) ..........
  540 do 600 j = lp1, en
         xr = wr(j-1)
         xi = wi(j-1)
c
         do 580 i = l, j
            yr = hr(i,j-1)
            yi = 0.0d0
            zzr = hr(i,j)
            zzi = hi(i,j)
            if (i .eq. j) go to 560
            yi = hi(i,j-1)
            hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
  560       hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  580    continue
c
  600 continue
c
      if (si .eq. 0.0d0) go to 240
c
      do 630 i = l, en
         yr = hr(i,en)
         yi = hi(i,en)
         hr(i,en) = sr * yr - si * yi
         hi(i,en) = sr * yi + si * yr
  630 continue
c
      go to 240
c     .......... a root found ..........
  660 wr(en) = hr(en,en) + tr
      wi(en) = hi(en,en) + ti
      en = enm1
      go to 220
c     .......... set error -- all eigenvalues have not
c                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end
      subroutine comqr2(nm,n,low,igh,ortr,orti,hr,hi,wr,wi,zr,zi,ierr)

c*********************************************************************72
c
cc COMQR2 gets eigenvalues/vectors of a complex upper Hessenberg matrix.
c
C  MESHED overflow control WITH vectors of isolated roots (10/19/89 BSG)
C  MESHED overflow control WITH triangular multiply (10/30/89 BSG)
c
c
c     this subroutine is a translation of a unitary analogue of the
c     algol procedure  comlr2, num. math. 16, 181-204(1970) by peters
c     and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
c     the unitary analogue substitutes the qr algorithm of francis
c     (comp. jour. 4, 332-345(1962)) for the lr algorithm.
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a complex upper hessenberg matrix by the qr
c     method.  the eigenvectors of a complex general matrix
c     can also be found if  corth  has been used to reduce
c     this general matrix to hessenberg form.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c         routine  cbal.  if  cbal  has not been used,
c          set low=1, igh=n.
c
c        ortr and orti contain information about the unitary trans-
c          formations used in the reduction by  corth, if performed.
c          only elements low through igh are used.  if the eigenvectors
c          of the hessenberg matrix are desired, set ortr(j) and
c          orti(j) to 0.0d0 for these elements.
c
c        hr and hi contain the real and imaginary parts,
c          respectively, of the complex upper hessenberg matrix.
c          their lower triangles below the subdiagonal contain further
c          information about the transformations which were used in the
c          reduction by  corth, if performed.  if the eigenvectors of
c          the hessenberg matrix are desired, these elements may be
c          arbitrary.
c
c     on output
c
c        ortr, orti, and the upper hessenberg portions of hr and hi
c          have been destroyed.
c
c        wr and wi contain the real and imaginary parts,
c          respectively, of the eigenvalues.  if an error
c          exit is made, the eigenvalues should be correct
c          for indices ierr+1,...,n.
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the eigenvectors.  the eigenvectors
c          are unnormalized.  if an error exit is made, none of
c          the eigenvectors has been found.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the limit of 30*n iterations is exhausted
c                     while the j-th eigenvalue is being sought.
c
c     calls cdiv for complex division.
c     calls csroot for complex square root.
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated october 1989.
c
      integer i,j,k,l,m,n,en,ii,jj,ll,nm,nn,igh,ip1,
     x        itn,its,low,lp1,enm1,iend,ierr
      double precision hr(nm,n),hi(nm,n),wr(n),wi(n),zr(nm,n),zi(nm,n),
     x       ortr(igh),orti(igh)
      double precision si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,norm,tst1,tst2,
     x       pythag

      ierr = 0
c     .......... initialize eigenvector matrix ..........
      do 101 j = 1, n
c
         do 100 i = 1, n
            zr(i,j) = 0.0d0
            zi(i,j) = 0.0d0
  100    continue
         zr(j,j) = 1.0d0
  101 continue
c     .......... form the matrix of accumulated transformations
c                from the information left by corth ..........
      iend = igh - low - 1
      if (iend) 180, 150, 105
c     .......... for i=igh-1 step -1 until low+1 do -- ..........
  105 do 140 ii = 1, iend
         i = igh - ii
         if (ortr(i) .eq. 0.0d0 .and. orti(i) .eq. 0.0d0) go to 140
         if (hr(i,i-1) .eq. 0.0d0 .and. hi(i,i-1) .eq. 0.0d0) go to 140
c     .......... norm below is negative of h formed in corth ..........
         norm = hr(i,i-1) * ortr(i) + hi(i,i-1) * orti(i)
         ip1 = i + 1
c
         do 110 k = ip1, igh
            ortr(k) = hr(k,i-1)
            orti(k) = hi(k,i-1)
  110    continue
c
         do 130 j = i, igh
            sr = 0.0d0
            si = 0.0d0
c
            do 115 k = i, igh
               sr = sr + ortr(k) * zr(k,j) + orti(k) * zi(k,j)
               si = si + ortr(k) * zi(k,j) - orti(k) * zr(k,j)
  115       continue
c
            sr = sr / norm
            si = si / norm
c
            do 120 k = i, igh
               zr(k,j) = zr(k,j) + sr * ortr(k) - si * orti(k)
               zi(k,j) = zi(k,j) + sr * orti(k) + si * ortr(k)
  120       continue
c
  130    continue
c
  140 continue
c     .......... create real subdiagonal elements ..........
  150 l = low + 1
c
      do 170 i = l, igh
         ll = min0(i+1,igh)
         if (hi(i,i-1) .eq. 0.0d0) go to 170
         norm = pythag(hr(i,i-1),hi(i,i-1))
         yr = hr(i,i-1) / norm
         yi = hi(i,i-1) / norm
         hr(i,i-1) = norm
         hi(i,i-1) = 0.0d0
c
         do 155 j = i, n
            si = yr * hi(i,j) - yi * hr(i,j)
            hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
            hi(i,j) = si
  155    continue
c
         do 160 j = 1, ll
            si = yr * hi(j,i) + yi * hr(j,i)
            hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
            hi(j,i) = si
  160    continue
c
         do 165 j = low, igh
            si = yr * zi(j,i) + yi * zr(j,i)
            zr(j,i) = yr * zr(j,i) - yi * zi(j,i)
            zi(j,i) = si
  165    continue
c
  170 continue
c     .......... store roots isolated by cbal ..........
  180 do 200 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 200
         wr(i) = hr(i,i)
         wi(i) = hi(i,i)
  200 continue
c
      en = igh
      tr = 0.0d0
      ti = 0.0d0
      itn = 30*n
c     .......... search for next eigenvalue ..........
  220 if (en .lt. low) go to 680
      its = 0
      enm1 = en - 1
c     .......... look for single small sub-diagonal element
c                for l=en step -1 until low do -- ..........
  240 do 260 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 300
         tst1 = dabs(hr(l-1,l-1)) + dabs(hi(l-1,l-1))
     x            + dabs(hr(l,l)) + dabs(hi(l,l))
         tst2 = tst1 + dabs(hr(l,l-1))
         if (tst2 .eq. tst1) go to 300
  260 continue
c     .......... form shift ..........
  300 if (l .eq. en) go to 660
      if (itn .eq. 0) go to 1000
      if (its .eq. 10 .or. its .eq. 20) go to 320
      sr = hr(en,en)
      si = hi(en,en)
      xr = hr(enm1,en) * hr(en,enm1)
      xi = hi(enm1,en) * hr(en,enm1)
      if (xr .eq. 0.0d0 .and. xi .eq. 0.0d0) go to 340
      yr = (hr(enm1,enm1) - sr) / 2.0d0
      yi = (hi(enm1,enm1) - si) / 2.0d0
      call csroot(yr**2-yi**2+xr,2.0d0*yr*yi+xi,zzr,zzi)
      if (yr * zzr + yi * zzi .ge. 0.0d0) go to 310
      zzr = -zzr
      zzi = -zzi
  310 call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)
      sr = sr - xr
      si = si - xi
      go to 340
c     .......... form exceptional shift ..........
  320 sr = dabs(hr(en,enm1)) + dabs(hr(enm1,en-2))
      si = 0.0d0
c
  340 do 360 i = low, en
         hr(i,i) = hr(i,i) - sr
         hi(i,i) = hi(i,i) - si
  360 continue
c
      tr = tr + sr
      ti = ti + si
      its = its + 1
      itn = itn - 1
c     .......... reduce to triangle (rows) ..........
      lp1 = l + 1
c
      do 500 i = lp1, en
         sr = hr(i,i-1)
         hr(i,i-1) = 0.0d0
         norm = pythag(pythag(hr(i-1,i-1),hi(i-1,i-1)),sr)
         xr = hr(i-1,i-1) / norm
         wr(i-1) = xr
         xi = hi(i-1,i-1) / norm
         wi(i-1) = xi
         hr(i-1,i-1) = norm
         hi(i-1,i-1) = 0.0d0
         hi(i,i-1) = sr / norm
c
         do 490 j = i, n
            yr = hr(i-1,j)
            yi = hi(i-1,j)
            zzr = hr(i,j)
            zzi = hi(i,j)
            hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
            hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
            hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
            hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
  490    continue
c
  500 continue
c
      si = hi(en,en)
      if (si .eq. 0.0d0) go to 540
      norm = pythag(hr(en,en),si)
      sr = hr(en,en) / norm
      si = si / norm
      hr(en,en) = norm
      hi(en,en) = 0.0d0
      if (en .eq. n) go to 540
      ip1 = en + 1
c
      do 520 j = ip1, n
         yr = hr(en,j)
         yi = hi(en,j)
         hr(en,j) = sr * yr + si * yi
         hi(en,j) = sr * yi - si * yr
  520 continue
c     .......... inverse operation (columns) ..........
  540 do 600 j = lp1, en
         xr = wr(j-1)
         xi = wi(j-1)
c
         do 580 i = 1, j
            yr = hr(i,j-1)
            yi = 0.0d0
            zzr = hr(i,j)
            zzi = hi(i,j)
            if (i .eq. j) go to 560
            yi = hi(i,j-1)
            hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
  560       hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  580    continue
c
         do 590 i = low, igh
            yr = zr(i,j-1)
            yi = zi(i,j-1)
            zzr = zr(i,j)
            zzi = zi(i,j)
            zr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            zi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
            zr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            zi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  590    continue
c
  600 continue
c
      if (si .eq. 0.0d0) go to 240
c
      do 630 i = 1, en
         yr = hr(i,en)
         yi = hi(i,en)
         hr(i,en) = sr * yr - si * yi
         hi(i,en) = sr * yi + si * yr
  630 continue
c
      do 640 i = low, igh
         yr = zr(i,en)
         yi = zi(i,en)
         zr(i,en) = sr * yr - si * yi
         zi(i,en) = sr * yi + si * yr
  640 continue
c
      go to 240
c     .......... a root found ..........
  660 hr(en,en) = hr(en,en) + tr
      wr(en) = hr(en,en)
      hi(en,en) = hi(en,en) + ti
      wi(en) = hi(en,en)
      en = enm1
      go to 220
c     .......... all roots found.  backsubstitute to find
c                vectors of upper triangular form ..........
  680 norm = 0.0d0
c
      do 720 i = 1, n
c
         do 720 j = i, n
            tr = dabs(hr(i,j)) + dabs(hi(i,j))
            if (tr .gt. norm) norm = tr
  720 continue
c
      if (n .eq. 1 .or. norm .eq. 0.0d0) go to 1001
c     .......... for en=n step -1 until 2 do -- ..........
      do 800 nn = 2, n
         en = n + 2 - nn
         xr = wr(en)
         xi = wi(en)
         hr(en,en) = 1.0d0
         hi(en,en) = 0.0d0
         enm1 = en - 1
c     .......... for i=en-1 step -1 until 1 do -- ..........
         do 780 ii = 1, enm1
            i = en - ii
            zzr = 0.0d0
            zzi = 0.0d0
            ip1 = i + 1
c
            do 740 j = ip1, en
               zzr = zzr + hr(i,j) * hr(j,en) - hi(i,j) * hi(j,en)
               zzi = zzi + hr(i,j) * hi(j,en) + hi(i,j) * hr(j,en)
  740       continue
c
            yr = xr - wr(i)
            yi = xi - wi(i)
            if (yr .ne. 0.0d0 .or. yi .ne. 0.0d0) go to 765
               tst1 = norm
               yr = tst1
  760          yr = 0.01d0 * yr
               tst2 = norm + yr
               if (tst2 .gt. tst1) go to 760
  765       continue
            call cdiv(zzr,zzi,yr,yi,hr(i,en),hi(i,en))
c     .......... overflow control ..........
            tr = dabs(hr(i,en)) + dabs(hi(i,en))
            if (tr .eq. 0.0d0) go to 780
            tst1 = tr
            tst2 = tst1 + 1.0d0/tst1
            if (tst2 .gt. tst1) go to 780
            do 770 j = i, en
               hr(j,en) = hr(j,en)/tr
               hi(j,en) = hi(j,en)/tr
  770       continue
c
  780    continue
c
  800 continue
c     .......... end backsubstitution ..........
c     .......... vectors of isolated roots ..........
      do  840 i = 1, N
         if (i .ge. low .and. i .le. igh) go to 840
c
         do 820 j = I, n
            zr(i,j) = hr(i,j)
            zi(i,j) = hi(i,j)
  820    continue
c
  840 continue
c     .......... multiply by transformation matrix to give
c                vectors of original full matrix.
c                for j=n step -1 until low do -- ..........
      do 880 jj = low, N
         j = n + low - jj
         m = min0(j,igh)
c
         do 880 i = low, igh
            zzr = 0.0d0
            zzi = 0.0d0
c
            do 860 k = low, m
               zzr = zzr + zr(i,k) * hr(k,j) - zi(i,k) * hi(k,j)
               zzi = zzi + zr(i,k) * hi(k,j) + zi(i,k) * hr(k,j)
  860       continue
c
            zr(i,j) = zzr
            zi(i,j) = zzi
  880 continue
c
      go to 1001
c     .......... set error -- all eigenvalues have not
c                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end
      subroutine cortb(nm,low,igh,ar,ai,ortr,orti,m,zr,zi)

c*********************************************************************72
c
cc CORTB determines eigenvectors by undoing the CORTH transformation.
c
c     this subroutine is a translation of a complex analogue of
c     the algol procedure ortbak, num. math. 12, 349-368(1968)
c     by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
c
c     this subroutine forms the eigenvectors of a complex general
c     matrix by back transforming those of the corresponding
c     upper hessenberg matrix determined by  corth.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        low and igh are integers determined by the balancing
c        routine  cbal.  if  cbal  has not been used,
c          set low=1 and igh equal to the order of the matrix.
c
c        ar and ai contain information about the unitary
c          transformations used in the reduction by  corth
c          in their strict lower triangles.
c
c        ortr and orti contain further information about the
c          transformations used in the reduction by  corth.
c          only elements low through igh are used.
c
c        m is the number of columns of zr and zi to be back transformed.
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the eigenvectors to be
c          back transformed in their first m columns.
c
c     on output
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the transformed eigenvectors
c          in their first m columns.
c
c        ortr and orti have been altered.
c
c     note that cortb preserves vector euclidean norms.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,m,la,mm,mp,nm,igh,kp1,low,mp1
      double precision ar(nm,igh),ai(nm,igh),ortr(igh),orti(igh),
     x       zr(nm,m),zi(nm,m)
      double precision h,gi,gr

      if (m .eq. 0) go to 200
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
c     .......... for mp=igh-1 step -1 until low+1 do -- ..........
      do 140 mm = kp1, la
         mp = low + igh - mm
         if (ar(mp,mp-1) .eq. 0.0d0 .and. ai(mp,mp-1) .eq. 0.0d0)
     x      go to 140
c     .......... h below is negative of h formed in corth ..........
         h = ar(mp,mp-1) * ortr(mp) + ai(mp,mp-1) * orti(mp)
         mp1 = mp + 1
c
         do 100 i = mp1, igh
            ortr(i) = ar(i,mp-1)
            orti(i) = ai(i,mp-1)
  100    continue
c
         do 130 j = 1, m
            gr = 0.0d0
            gi = 0.0d0
c
            do 110 i = mp, igh
               gr = gr + ortr(i) * zr(i,j) + orti(i) * zi(i,j)
               gi = gi + ortr(i) * zi(i,j) - orti(i) * zr(i,j)
  110       continue
c
            gr = gr / h
            gi = gi / h
c
            do 120 i = mp, igh
               zr(i,j) = zr(i,j) + gr * ortr(i) - gi * orti(i)
               zi(i,j) = zi(i,j) + gr * orti(i) + gi * ortr(i)
  120       continue
c
  130    continue
c
  140 continue
c
  200 return
      end
      subroutine corth(nm,n,low,igh,ar,ai,ortr,orti)

c*********************************************************************72
c
cc CORTH transforms a complex general matrix to upper Hessenberg form.
c
c     this subroutine is a translation of a complex analogue of
c     the algol procedure orthes, num. math. 12, 349-368(1968)
c     by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
c
c     given a complex general matrix, this subroutine
c     reduces a submatrix situated in rows and columns
c     low through igh to upper hessenberg form by
c     unitary similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c        routine  cbal.  if  cbal  has not been used,
c          set low=1, igh=n.
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the complex input matrix.
c
c     on output
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the hessenberg matrix.  information
c          about the unitary transformations used in the reduction
c          is stored in the remaining triangles under the
c          hessenberg matrix.
c
c        ortr and orti contain further information about the
c          transformations.  only elements low through igh are used.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,m,n,ii,jj,la,mp,nm,igh,kp1,low
      double precision ar(nm,n),ai(nm,n),ortr(igh),orti(igh)
      double precision f,g,h,fi,fr,scale,pythag

      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
c
      do 180 m = kp1, la
         h = 0.0d0
         ortr(m) = 0.0d0
         orti(m) = 0.0d0
         scale = 0.0d0
c     .......... scale column (algol tol then not needed) ..........
         do 90 i = m, igh
   90    scale = scale + dabs(ar(i,m-1)) + dabs(ai(i,m-1))
c
         if (scale .eq. 0.0d0) go to 180
         mp = m + igh
c     .......... for i=igh step -1 until m do -- ..........
         do 100 ii = m, igh
            i = mp - ii
            ortr(i) = ar(i,m-1) / scale
            orti(i) = ai(i,m-1) / scale
            h = h + ortr(i) * ortr(i) + orti(i) * orti(i)
  100    continue
c
         g = dsqrt(h)
         f = pythag(ortr(m),orti(m))
         if (f .eq. 0.0d0) go to 103
         h = h + f * g
         g = g / f
         ortr(m) = (1.0d0 + g) * ortr(m)
         orti(m) = (1.0d0 + g) * orti(m)
         go to 105
c
  103    ortr(m) = g
         ar(m,m-1) = scale
c     .......... form (i-(u*ut)/h) * a ..........
  105    do 130 j = m, n
            fr = 0.0d0
            fi = 0.0d0
c     .......... for i=igh step -1 until m do -- ..........
            do 110 ii = m, igh
               i = mp - ii
               fr = fr + ortr(i) * ar(i,j) + orti(i) * ai(i,j)
               fi = fi + ortr(i) * ai(i,j) - orti(i) * ar(i,j)
  110       continue
c
            fr = fr / h
            fi = fi / h
c
            do 120 i = m, igh
               ar(i,j) = ar(i,j) - fr * ortr(i) + fi * orti(i)
               ai(i,j) = ai(i,j) - fr * orti(i) - fi * ortr(i)
  120       continue
c
  130    continue
c     .......... form (i-(u*ut)/h)*a*(i-(u*ut)/h) ..........
         do 160 i = 1, igh
            fr = 0.0d0
            fi = 0.0d0
c     .......... for j=igh step -1 until m do -- ..........
            do 140 jj = m, igh
               j = mp - jj
               fr = fr + ortr(j) * ar(i,j) - orti(j) * ai(i,j)
               fi = fi + ortr(j) * ai(i,j) + orti(j) * ar(i,j)
  140       continue
c
            fr = fr / h
            fi = fi / h
c
            do 150 j = m, igh
               ar(i,j) = ar(i,j) - fr * ortr(j) - fi * orti(j)
               ai(i,j) = ai(i,j) + fr * orti(j) - fi * ortr(j)
  150       continue
c
  160    continue
c
         ortr(m) = scale * ortr(m)
         orti(m) = scale * orti(m)
         ar(m,m-1) = -g * ar(m,m-1)
         ai(m,m-1) = -g * ai(m,m-1)
  180 continue
c
  200 return
      end
      subroutine csroot(xr,xi,yr,yi)

c*********************************************************************72
c
cc CSROOT computes the complex square root of a complex quantity.
c
c     (yr,yi) = complex dsqrt(xr,xi)
c     branch chosen so that yr .ge. 0.0 and sign(yi) .eq. sign(xi)
c
      double precision xr,xi,yr,yi
      double precision s,tr,ti,pythag

      tr = xr
      ti = xi
      s = dsqrt(0.5d0*(pythag(tr,ti) + dabs(tr)))
      if (tr .ge. 0.0d0) yr = s
      if (ti .lt. 0.0d0) s = -s
      if (tr .le. 0.0d0) yi = s
      if (tr .lt. 0.0d0) yr = 0.5d0*(ti/yi)
      if (tr .gt. 0.0d0) yi = 0.5d0*(ti/yr)
      return
      end
      subroutine elmbak(nm,low,igh,a,int,m,z)

c*********************************************************************72
c
cc ELMBAK determines eigenvectors by undoing the ELMHES transformation.
c
c     this subroutine is a translation of the algol procedure elmbak,
c     num. math. 12, 349-368(1968) by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
c
c     this subroutine forms the eigenvectors of a real general
c     matrix by back transforming those of the corresponding
c     upper hessenberg matrix determined by  elmhes.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        low and igh are integers determined by the balancing
c        routine  balanc.  if  balanc  has not been used,
c          set low=1 and igh equal to the order of the matrix.
c
c        a contains the multipliers which were used in the
c          reduction by  elmhes  in its lower triangle
c          below the subdiagonal.
c
c        int contains information on the rows and columns
c          interchanged in the reduction by  elmhes.
c          only elements low through igh are used.
c
c        m is the number of columns of z to be back transformed.
c
c        z contains the real and imaginary parts of the eigen-
c          vectors to be back transformed in its first m columns.
c
c     on output
c
c        z contains the real and imaginary parts of the
c          transformed eigenvectors in its first m columns.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c

      integer i,j,m,la,mm,mp,nm,igh,kp1,low,mp1
      double precision a(nm,igh),z(nm,m)
      double precision x
      integer int(igh)

      if (m .eq. 0) go to 200
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
c     .......... for mp=igh-1 step -1 until low+1 do -- ..........
      do 140 mm = kp1, la
         mp = low + igh - mm
         mp1 = mp + 1
c
         do 110 i = mp1, igh
            x = a(i,mp-1)
            if (x .eq. 0.0d0) go to 110
c
            do 100 j = 1, m
  100       z(i,j) = z(i,j) + x * z(mp,j)
c
  110    continue
c
         i = int(mp)
         if (i .eq. mp) go to 140
c
         do 130 j = 1, m
            x = z(i,j)
            z(i,j) = z(mp,j)
            z(mp,j) = x
  130    continue
c
  140 continue
c
  200 return
      end
      subroutine elmhes ( nm, n, low, igh, a, int )

c*********************************************************************72
c
cc ELMHES transforms a real general matrix to upper Hessenberg form.
c
c     this subroutine is a translation of the algol procedure elmhes,
c     num. math. 12, 349-368(1968) by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
c
c     given a real general matrix, this subroutine
c     reduces a submatrix situated in rows and columns
c     low through igh to upper hessenberg form by
c     stabilized elementary similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c        routine  balanc.  if  balanc  has not been used,
c          set low=1, igh=n.
c
c        a contains the input matrix.
c
c     on output
c
c        a contains the hessenberg matrix.  the multipliers
c          which were used in the reduction are stored in the
c          remaining triangle under the hessenberg matrix.
c
c        int contains information on the rows and columns
c          interchanged in the reduction.
c          only elements low through igh are used.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c

      integer i,j,m,n,la,nm,igh,kp1,low,mm1,mp1
      double precision a(nm,n)
      double precision x,y
      integer int(igh)

      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200

      do 180 m = kp1, la
         mm1 = m - 1
         x = 0.0d0
         i = m

         do 100 j = m, igh
            if (dabs(a(j,mm1)) .le. dabs(x)) go to 100
            x = a(j,mm1)
            i = j
  100    continue

         int(m) = i
         if (i .eq. m) go to 130
c
c  interchange rows and columns of a.
c
         do 110 j = mm1, n
            y = a(i,j)
            a(i,j) = a(m,j)
            a(m,j) = y
  110    continue

         do 120 j = 1, igh
            y = a(j,i)
            a(j,i) = a(j,m)
            a(j,m) = y
  120    continue
c
c  end interchange.
c
  130    if (x .eq. 0.0d0) go to 180
         mp1 = m + 1

         do 160 i = mp1, igh
            y = a(i,mm1)
            if (y .eq. 0.0d0) go to 160
            y = y / x
            a(i,mm1) = y

            do 140 j = m, n
  140       a(i,j) = a(i,j) - y * a(m,j)

            do 150 j = 1, igh
  150       a(j,m) = a(j,m) + y * a(j,i)

  160    continue

  180 continue

  200 continue

      return
      end
      subroutine eltran(nm,n,low,igh,a,int,z)

c*********************************************************************72
c
cc ELTRAN accumulates similarity transformations used by ELMHES.
c
c     this subroutine is a translation of the algol procedure elmtrans,
c     num. math. 16, 181-204(1970) by peters and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
c
c     this subroutine accumulates the stabilized elementary
c     similarity transformations used in the reduction of a
c     real general matrix to upper hessenberg form by  elmhes.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c        routine  balanc.  if  balanc  has not been used,
c          set low=1, igh=n.
c
c        a contains the multipliers which were used in the
c          reduction by  elmhes  in its lower triangle
c          below the subdiagonal.
c
c        int contains information on the rows and columns
c          interchanged in the reduction by  elmhes.
c          only elements low through igh are used.
c
c     on output
c
c        z contains the transformation matrix produced in the
c          reduction by  elmhes.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,n,kl,mm,mp,nm,igh,low,mp1
      double precision a(nm,igh),z(nm,n)
      integer int(igh)
c
c     .......... initialize z to identity matrix ..........
      do 80 j = 1, n
c
         do 60 i = 1, n
   60    z(i,j) = 0.0d0
c
         z(j,j) = 1.0d0
   80 continue
c
      kl = igh - low - 1
      if (kl .lt. 1) go to 200
c     .......... for mp=igh-1 step -1 until low+1 do -- ..........
      do 140 mm = 1, kl
         mp = igh - mm
         mp1 = mp + 1
c
         do 100 i = mp1, igh
  100    z(i,mp) = a(i,mp-1)
c
         i = int(mp)
         if (i .eq. mp) go to 140
c
         do 130 j = mp, igh
            z(mp,j) = z(i,j)
            z(i,j) = 0.0d0
  130    continue
c
         z(i,mp) = 1.0d0
  140 continue
c
  200 return
      end
      double precision function epslon (x)

c*********************************************************************72
c
cc EPSLON estimate unit roundoff in quantities of size X.
c
c     this program should function properly on all systems
c     satisfying the following two assumptions,
c        1.  the base used in representing floating point
c            numbers is not a power of three.
c        2.  the quantity  a  in statement 10 is represented to
c            the accuracy used in floating point variables
c            that are stored in memory.
c     the statement number 10 and the go to 10 are intended to
c     force optimizing compilers to generate code satisfying
c     assumption 2.
c     under these assumptions, it should be true that,
c            a  is not exactly equal to four-thirds,
c            b  has a zero for its last bit or digit,
c            c  is not exactly equal to one,
c            eps  measures the separation of 1.0 from
c                 the next larger floating point number.
c     the developers of eispack would appreciate being informed
c     about any systems where these assumptions do not hold.
c
c     this version dated 4/6/83.
c
      double precision a,b,c,eps
      double precision x

      a = 4.0d0/3.0d0
   10 b = a - 1.0d0
      c = b + b + b
      eps = dabs(c-1.0d0)
      if (eps .eq. 0.0d0) go to 10
      epslon = eps*dabs(x)
      return
      end
      subroutine figi(nm,n,t,d,e,e2,ierr)

c*********************************************************************72
c
cc FIGI transforms a real nonsymmetric tridiagonal matrix to symmetric form.
c
c     given a nonsymmetric tridiagonal matrix such that the products
c     of corresponding pairs of off-diagonal elements are all
c     non-negative, this subroutine reduces it to a symmetric
c     tridiagonal matrix with the same eigenvalues.  if, further,
c     a zero product only occurs when both factors are zero,
c     the reduced matrix is similar to the original matrix.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        t contains the input matrix.  its subdiagonal is
c          stored in the last n-1 positions of the first column,
c          its diagonal in the n positions of the second column,
c          and its superdiagonal in the first n-1 positions of
c          the third column.  t(1,1) and t(n,3) are arbitrary.
c
c     on output
c
c        t is unaltered.
c
c        d contains the diagonal elements of the symmetric matrix.
c
c        e contains the subdiagonal elements of the symmetric
c          matrix in its last n-1 positions.  e(1) is not set.
c
c        e2 contains the squares of the corresponding elements of e.
c          e2 may coincide with e if the squares are not needed.
c
c        ierr is set to
c          zero       for normal return,
c          n+i        if t(i,1)*t(i-1,3) is negative,
c          -(3*n+i)   if t(i,1)*t(i-1,3) is zero with one factor
c                     non-zero.  in this case, the eigenvectors of
c                     the symmetric matrix are not simply related
c                     to those of  t  and should not be sought.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c

      integer i,n,nm,ierr
      double precision t(nm,3),d(n),e(n),e2(n)

      ierr = 0
c
      do 100 i = 1, n
         if (i .eq. 1) go to 90
         e2(i) = t(i,1) * t(i-1,3)
         if (e2(i)) 1000, 60, 80
   60    if (t(i,1) .eq. 0.0d0 .and. t(i-1,3) .eq. 0.0d0) go to 80
c     .......... set error -- product of some pair of off-diagonal
c                elements is zero with one member non-zero ..........
         ierr = -(3 * n + i)
   80    e(i) = dsqrt(e2(i))
   90    d(i) = t(i,2)
  100 continue
c
      go to 1001
c     .......... set error -- product of some pair of off-diagonal
c                elements is negative ..........
 1000 ierr = n + i
 1001 return
      end
      subroutine figi2(nm,n,t,d,e,z,ierr)

c*********************************************************************72
c
cc FIGI2 transforms a real nonsymmetric tridiagonal matrix to symmetric form.
c
c     given a nonsymmetric tridiagonal matrix such that the products
c     of corresponding pairs of off-diagonal elements are all
c     non-negative, and zero only when both factors are zero, this
c     routine reduces it to a symmetric tridiagonal matrix
c     using and accumulating diagonal similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        t contains the input matrix.  its subdiagonal is
c          stored in the last n-1 positions of the first column,
c          its diagonal in the n positions of the second column,
c          and its superdiagonal in the first n-1 positions of
c          the third column.  t(1,1) and t(n,3) are arbitrary.
c
c     on output
c
c        t is unaltered.
c
c        d contains the diagonal elements of the symmetric matrix.
c
c        e contains the subdiagonal elements of the symmetric
c          matrix in its last n-1 positions.  e(1) is not set.
c
c        z contains the transformation matrix produced in
c          the reduction.
c
c        ierr is set to
c          zero       for normal return,
c          n+i        if t(i,1)*t(i-1,3) is negative,
c          2*n+i      if t(i,1)*t(i-1,3) is zero with
c                     one factor non-zero.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,n,nm,ierr
      double precision t(nm,3),d(n),e(n),z(nm,n)
      double precision h

      ierr = 0
c
      do 100 i = 1, n
c
         do 50 j = 1, n
   50    z(i,j) = 0.0d0
c
         if (i .eq. 1) go to 70
         h = t(i,1) * t(i-1,3)
         if (h) 900, 60, 80
   60    if (t(i,1) .ne. 0.0d0 .or. t(i-1,3) .ne. 0.0d0) go to 1000
         e(i) = 0.0d0
   70    z(i,i) = 1.0d0
         go to 90
   80    e(i) = dsqrt(h)
         z(i,i) = z(i-1,i-1) * e(i) / t(i-1,3)
   90    d(i) = t(i,2)
  100 continue
c
      go to 1001
c     .......... set error -- product of some pair of off-diagonal
c                elements is negative ..........
  900 ierr = n + i
      go to 1001
c     .......... set error -- product of some pair of off-diagonal
c                elements is zero with one member non-zero ..........
 1000 ierr = 2 * n + i
 1001 return
      end
      subroutine hqr(nm,n,low,igh,h,wr,wi,ierr)

c*********************************************************************72
c
cc HQR computes all eigenvalues of a real upper Hessenberg matrix.
c
C  RESTORED CORRECT INDICES OF LOOPS (200,210,230,240). (9/29/89 BSG)
c
c     this subroutine is a translation of the algol procedure hqr,
c     num. math. 14, 219-231(1970) by martin, peters, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 359-371(1971).
c
c     this subroutine finds the eigenvalues of a real
c     upper hessenberg matrix by the qr method.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c        routine  balanc.  if  balanc  has not been used,
c          set low=1, igh=n.
c
c        h contains the upper hessenberg matrix.  information about
c          the transformations used in the reduction to hessenberg
c          form by  elmhes  or  orthes, if performed, is stored
c          in the remaining triangle under the hessenberg matrix.
c
c     on output
c
c        h has been destroyed.  therefore, it must be saved
c          before calling  hqr  if subsequent calculation and
c          back transformation of eigenvectors is to be performed.
c
c        wr and wi contain the real and imaginary parts,
c          respectively, of the eigenvalues.  the eigenvalues
c          are unordered except that complex conjugate pairs
c          of values appear consecutively with the eigenvalue
c          having the positive imaginary part first.  if an
c          error exit is made, the eigenvalues should be correct
c          for indices ierr+1,...,n.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the limit of 30*n iterations is exhausted
c                     while the j-th eigenvalue is being sought.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated september 1989.
c
      integer i,j,k,l,m,n,en,ll,mm,na,nm,igh,itn,its,low,mp2,enm2,ierr
      double precision h(nm,n),wr(n),wi(n)
      double precision p,q,r,s,t,w,x,y,zz,norm,tst1,tst2
      logical notlas

      ierr = 0
      norm = 0.0d0
      k = 1
c     .......... store roots isolated by balanc
c                and compute matrix norm ..........
      do 50 i = 1, n
c
         do 40 j = k, n
   40    norm = norm + dabs(h(i,j))
c
         k = i
         if (i .ge. low .and. i .le. igh) go to 50
         wr(i) = h(i,i)
         wi(i) = 0.0d0
   50 continue
c
      en = igh
      t = 0.0d0
      itn = 30*n
c     .......... search for next eigenvalues ..........
   60 if (en .lt. low) go to 1001
      its = 0
      na = en - 1
      enm2 = na - 1
c     .......... look for single small sub-diagonal element
c                for l=en step -1 until low do -- ..........
   70 do 80 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 100
         s = dabs(h(l-1,l-1)) + dabs(h(l,l))
         if (s .eq. 0.0d0) s = norm
         tst1 = s
         tst2 = tst1 + dabs(h(l,l-1))
         if (tst2 .eq. tst1) go to 100
   80 continue
c     .......... form shift ..........
  100 x = h(en,en)
      if (l .eq. en) go to 270
      y = h(na,na)
      w = h(en,na) * h(na,en)
      if (l .eq. na) go to 280
      if (itn .eq. 0) go to 1000
      if (its .ne. 10 .and. its .ne. 20) go to 130
c     .......... form exceptional shift ..........
      t = t + x
c
      do 120 i = low, en
  120 h(i,i) = h(i,i) - x
c
      s = dabs(h(en,na)) + dabs(h(na,enm2))
      x = 0.75d0 * s
      y = x
      w = -0.4375d0 * s * s
  130 its = its + 1
      itn = itn - 1
c     .......... look for two consecutive small
c                sub-diagonal elements.
c                for m=en-2 step -1 until l do -- ..........
      do 140 mm = l, enm2
         m = enm2 + l - mm
         zz = h(m,m)
         r = x - zz
         s = y - zz
         p = (r * s - w) / h(m+1,m) + h(m,m+1)
         q = h(m+1,m+1) - zz - r - s
         r = h(m+2,m+1)
         s = dabs(p) + dabs(q) + dabs(r)
         p = p / s
         q = q / s
         r = r / s
         if (m .eq. l) go to 150
         tst1 = dabs(p)*(dabs(h(m-1,m-1)) + dabs(zz) + dabs(h(m+1,m+1)))
         tst2 = tst1 + dabs(h(m,m-1))*(dabs(q) + dabs(r))
         if (tst2 .eq. tst1) go to 150
  140 continue
c
  150 mp2 = m + 2
c
      do 160 i = mp2, en
         h(i,i-2) = 0.0d0
         if (i .eq. mp2) go to 160
         h(i,i-3) = 0.0d0
  160 continue
c     .......... double qr step involving rows l to en and
c                columns m to en ..........
      do 260 k = m, na
         notlas = k .ne. na
         if (k .eq. m) go to 170
         p = h(k,k-1)
         q = h(k+1,k-1)
         r = 0.0d0
         if (notlas) r = h(k+2,k-1)
         x = dabs(p) + dabs(q) + dabs(r)
         if (x .eq. 0.0d0) go to 260
         p = p / x
         q = q / x
         r = r / x
  170    s = dsign(dsqrt(p*p+q*q+r*r),p)
         if (k .eq. m) go to 180
         h(k,k-1) = -s * x
         go to 190
  180    if (l .ne. m) h(k,k-1) = -h(k,k-1)
  190    p = p + s
         x = p / s
         y = q / s
         zz = r / s
         q = q / p
         r = r / p
         if (notlas) go to 225
c     .......... row modification ..........
         do 200 j = k, EN
            p = h(k,j) + q * h(k+1,j)
            h(k,j) = h(k,j) - p * x
            h(k+1,j) = h(k+1,j) - p * y
  200    continue
c
         j = min0(en,k+3)
c     .......... column modification ..........
         do 210 i = L, j
            p = x * h(i,k) + y * h(i,k+1)
            h(i,k) = h(i,k) - p
            h(i,k+1) = h(i,k+1) - p * q
  210    continue
         go to 255
  225    continue
c     .......... row modification ..........
         do 230 j = k, EN
            p = h(k,j) + q * h(k+1,j) + r * h(k+2,j)
            h(k,j) = h(k,j) - p * x
            h(k+1,j) = h(k+1,j) - p * y
            h(k+2,j) = h(k+2,j) - p * zz
  230    continue
c
         j = min0(en,k+3)
c     .......... column modification ..........
         do 240 i = L, j
            p = x * h(i,k) + y * h(i,k+1) + zz * h(i,k+2)
            h(i,k) = h(i,k) - p
            h(i,k+1) = h(i,k+1) - p * q
            h(i,k+2) = h(i,k+2) - p * r
  240    continue
  255    continue
c
  260 continue
c
      go to 70
c     .......... one root found ..........
  270 wr(en) = x + t
      wi(en) = 0.0d0
      en = na
      go to 60
c     .......... two roots found ..........
  280 p = (y - x) / 2.0d0
      q = p * p + w
      zz = dsqrt(dabs(q))
      x = x + t
      if (q .lt. 0.0d0) go to 320
c     .......... real pair ..........
      zz = p + dsign(zz,p)
      wr(na) = x + zz
      wr(en) = wr(na)
      if (zz .ne. 0.0d0) wr(en) = x - w / zz
      wi(na) = 0.0d0
      wi(en) = 0.0d0
      go to 330
c     .......... complex pair ..........
  320 wr(na) = x + p
      wr(en) = x + p
      wi(na) = zz
      wi(en) = -zz
  330 en = enm2
      go to 60
c     .......... set error -- all eigenvalues have not
c                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end
      subroutine hqr2(nm,n,low,igh,h,wr,wi,z,ierr)

c*********************************************************************72
c
cc HQR2 computes eigenvalues and eigenvectors of a real upper Hessenberg matrix.
c
c     this subroutine is a translation of the algol procedure hqr2,
c     num. math. 16, 181-204(1970) by peters and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a real upper hessenberg matrix by the qr method.  the
c     eigenvectors of a real general matrix can also be found
c     if  elmhes  and  eltran  or  orthes  and  ortran  have
c     been used to reduce this general matrix to hessenberg form
c     and to accumulate the similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c        routine  balanc.  if  balanc  has not been used,
c          set low=1, igh=n.
c
c        h contains the upper hessenberg matrix.
c
c        z contains the transformation matrix produced by  eltran
c          after the reduction by  elmhes, or by  ortran  after the
c          reduction by  orthes, if performed.  if the eigenvectors
c          of the hessenberg matrix are desired, z must contain the
c          identity matrix.
c
c     on output
c
c        h has been destroyed.
c
c        wr and wi contain the real and imaginary parts,
c          respectively, of the eigenvalues.  the eigenvalues
c          are unordered except that complex conjugate pairs
c          of values appear consecutively with the eigenvalue
c          having the positive imaginary part first.  if an
c          error exit is made, the eigenvalues should be correct
c          for indices ierr+1,...,n.
c
c        z contains the real and imaginary parts of the eigenvectors.
c          if the i-th eigenvalue is real, the i-th column of z
c          contains its eigenvector.  if the i-th eigenvalue is complex
c          with positive imaginary part, the i-th and (i+1)-th
c          columns of z contain the real and imaginary parts of its
c          eigenvector.  the eigenvectors are unnormalized.  if an
c          error exit is made, none of the eigenvectors has been found.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the limit of 30*n iterations is exhausted
c                     while the j-th eigenvalue is being sought.
c
c     calls cdiv for complex division.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c

      integer i,j,k,l,m,n,en,ii,jj,ll,mm,na,nm,nn,
     x        igh,itn,its,low,mp2,enm2,ierr
      double precision h(nm,n),wr(n),wi(n),z(nm,n)
      double precision p,q,r,s,t,w,x,y,ra,sa,vi,vr,zz,norm,tst1,tst2
      logical notlas

      ierr = 0
      norm = 0.0d0
      k = 1
c     .......... store roots isolated by balanc
c                and compute matrix norm ..........
      do 50 i = 1, n
c
         do 40 j = k, n
   40    norm = norm + dabs(h(i,j))
c
         k = i
         if (i .ge. low .and. i .le. igh) go to 50
         wr(i) = h(i,i)
         wi(i) = 0.0d0
   50 continue
c
      en = igh
      t = 0.0d0
      itn = 30*n
c     .......... search for next eigenvalues ..........
   60 if (en .lt. low) go to 340
      its = 0
      na = en - 1
      enm2 = na - 1
c     .......... look for single small sub-diagonal element
c                for l=en step -1 until low do -- ..........
   70 do 80 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 100
         s = dabs(h(l-1,l-1)) + dabs(h(l,l))
         if (s .eq. 0.0d0) s = norm
         tst1 = s
         tst2 = tst1 + dabs(h(l,l-1))
         if (tst2 .eq. tst1) go to 100
   80 continue
c     .......... form shift ..........
  100 x = h(en,en)
      if (l .eq. en) go to 270
      y = h(na,na)
      w = h(en,na) * h(na,en)
      if (l .eq. na) go to 280
      if (itn .eq. 0) go to 1000
      if (its .ne. 10 .and. its .ne. 20) go to 130
c     .......... form exceptional shift ..........
      t = t + x
c
      do 120 i = low, en
  120 h(i,i) = h(i,i) - x
c
      s = dabs(h(en,na)) + dabs(h(na,enm2))
      x = 0.75d0 * s
      y = x
      w = -0.4375d0 * s * s
  130 its = its + 1
      itn = itn - 1
c     .......... look for two consecutive small
c                sub-diagonal elements.
c                for m=en-2 step -1 until l do -- ..........
      do 140 mm = l, enm2
         m = enm2 + l - mm
         zz = h(m,m)
         r = x - zz
         s = y - zz
         p = (r * s - w) / h(m+1,m) + h(m,m+1)
         q = h(m+1,m+1) - zz - r - s
         r = h(m+2,m+1)
         s = dabs(p) + dabs(q) + dabs(r)
         p = p / s
         q = q / s
         r = r / s
         if (m .eq. l) go to 150
         tst1 = dabs(p)*(dabs(h(m-1,m-1)) + dabs(zz) + dabs(h(m+1,m+1)))
         tst2 = tst1 + dabs(h(m,m-1))*(dabs(q) + dabs(r))
         if (tst2 .eq. tst1) go to 150
  140 continue
c
  150 mp2 = m + 2
c
      do 160 i = mp2, en
         h(i,i-2) = 0.0d0
         if (i .eq. mp2) go to 160
         h(i,i-3) = 0.0d0
  160 continue
c     .......... double qr step involving rows l to en and
c                columns m to en ..........
      do 260 k = m, na
         notlas = k .ne. na
         if (k .eq. m) go to 170
         p = h(k,k-1)
         q = h(k+1,k-1)
         r = 0.0d0
         if (notlas) r = h(k+2,k-1)
         x = dabs(p) + dabs(q) + dabs(r)
         if (x .eq. 0.0d0) go to 260
         p = p / x
         q = q / x
         r = r / x
  170    s = dsign(dsqrt(p*p+q*q+r*r),p)
         if (k .eq. m) go to 180
         h(k,k-1) = -s * x
         go to 190
  180    if (l .ne. m) h(k,k-1) = -h(k,k-1)
  190    p = p + s
         x = p / s
         y = q / s
         zz = r / s
         q = q / p
         r = r / p
         if (notlas) go to 225
c     .......... row modification ..........
         do 200 j = k, n
            p = h(k,j) + q * h(k+1,j)
            h(k,j) = h(k,j) - p * x
            h(k+1,j) = h(k+1,j) - p * y
  200    continue
c
         j = min0(en,k+3)
c     .......... column modification ..........
         do 210 i = 1, j
            p = x * h(i,k) + y * h(i,k+1)
            h(i,k) = h(i,k) - p
            h(i,k+1) = h(i,k+1) - p * q
  210    continue
c     .......... accumulate transformations ..........
         do 220 i = low, igh
            p = x * z(i,k) + y * z(i,k+1)
            z(i,k) = z(i,k) - p
            z(i,k+1) = z(i,k+1) - p * q
  220    continue
         go to 255
  225    continue
c     .......... row modification ..........
         do 230 j = k, n
            p = h(k,j) + q * h(k+1,j) + r * h(k+2,j)
            h(k,j) = h(k,j) - p * x
            h(k+1,j) = h(k+1,j) - p * y
            h(k+2,j) = h(k+2,j) - p * zz
  230    continue
c
         j = min0(en,k+3)
c     .......... column modification ..........
         do 240 i = 1, j
            p = x * h(i,k) + y * h(i,k+1) + zz * h(i,k+2)
            h(i,k) = h(i,k) - p
            h(i,k+1) = h(i,k+1) - p * q
            h(i,k+2) = h(i,k+2) - p * r
  240    continue
c     .......... accumulate transformations ..........
         do 250 i = low, igh
            p = x * z(i,k) + y * z(i,k+1) + zz * z(i,k+2)
            z(i,k) = z(i,k) - p
            z(i,k+1) = z(i,k+1) - p * q
            z(i,k+2) = z(i,k+2) - p * r
  250    continue
  255    continue
c
  260 continue
c
      go to 70
c     .......... one root found ..........
  270 h(en,en) = x + t
      wr(en) = h(en,en)
      wi(en) = 0.0d0
      en = na
      go to 60
c     .......... two roots found ..........
  280 p = (y - x) / 2.0d0
      q = p * p + w
      zz = dsqrt(dabs(q))
      h(en,en) = x + t
      x = h(en,en)
      h(na,na) = y + t
      if (q .lt. 0.0d0) go to 320
c     .......... real pair ..........
      zz = p + dsign(zz,p)
      wr(na) = x + zz
      wr(en) = wr(na)
      if (zz .ne. 0.0d0) wr(en) = x - w / zz
      wi(na) = 0.0d0
      wi(en) = 0.0d0
      x = h(en,na)
      s = dabs(x) + dabs(zz)
      p = x / s
      q = zz / s
      r = dsqrt(p*p+q*q)
      p = p / r
      q = q / r
c     .......... row modification ..........
      do 290 j = na, n
         zz = h(na,j)
         h(na,j) = q * zz + p * h(en,j)
         h(en,j) = q * h(en,j) - p * zz
  290 continue
c     .......... column modification ..........
      do 300 i = 1, en
         zz = h(i,na)
         h(i,na) = q * zz + p * h(i,en)
         h(i,en) = q * h(i,en) - p * zz
  300 continue
c     .......... accumulate transformations ..........
      do 310 i = low, igh
         zz = z(i,na)
         z(i,na) = q * zz + p * z(i,en)
         z(i,en) = q * z(i,en) - p * zz
  310 continue
c
      go to 330
c     .......... complex pair ..........
  320 wr(na) = x + p
      wr(en) = x + p
      wi(na) = zz
      wi(en) = -zz
  330 en = enm2
      go to 60
c     .......... all roots found.  backsubstitute to find
c                vectors of upper triangular form ..........
  340 if (norm .eq. 0.0d0) go to 1001
c     .......... for en=n step -1 until 1 do -- ..........
      do 800 nn = 1, n
         en = n + 1 - nn
         p = wr(en)
         q = wi(en)
         na = en - 1
         if (q) 710, 600, 800
c     .......... real vector ..........
  600    m = en
         h(en,en) = 1.0d0
         if (na .eq. 0) go to 800
c     .......... for i=en-1 step -1 until 1 do -- ..........
         do 700 ii = 1, na
            i = en - ii
            w = h(i,i) - p
            r = 0.0d0
c
            do 610 j = m, en
  610       r = r + h(i,j) * h(j,en)
c
            if (wi(i) .ge. 0.0d0) go to 630
            zz = w
            s = r
            go to 700
  630       m = i
            if (wi(i) .ne. 0.0d0) go to 640
            t = w
            if (t .ne. 0.0d0) go to 635
               tst1 = norm
               t = tst1
  632          t = 0.01d0 * t
               tst2 = norm + t
               if (tst2 .gt. tst1) go to 632
  635       h(i,en) = -r / t
            go to 680
c     .......... solve real equations ..........
  640       x = h(i,i+1)
            y = h(i+1,i)
            q = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i)
            t = (x * s - zz * r) / q
            h(i,en) = t
            if (dabs(x) .le. dabs(zz)) go to 650
            h(i+1,en) = (-r - w * t) / x
            go to 680
  650       h(i+1,en) = (-s - y * t) / zz
c
c     .......... overflow control ..........
  680       t = dabs(h(i,en))
            if (t .eq. 0.0d0) go to 700
            tst1 = t
            tst2 = tst1 + 1.0d0/tst1
            if (tst2 .gt. tst1) go to 700
            do 690 j = i, en
               h(j,en) = h(j,en)/t
  690       continue
c
  700    continue
c     .......... end real vector ..........
         go to 800
c     .......... complex vector ..........
  710    m = na
c     .......... last vector component chosen imaginary so that
c                eigenvector matrix is triangular ..........
         if (dabs(h(en,na)) .le. dabs(h(na,en))) go to 720
         h(na,na) = q / h(en,na)
         h(na,en) = -(h(en,en) - p) / h(en,na)
         go to 730
  720    call cdiv(0.0d0,-h(na,en),h(na,na)-p,q,h(na,na),h(na,en))
  730    h(en,na) = 0.0d0
         h(en,en) = 1.0d0
         enm2 = na - 1
         if (enm2 .eq. 0) go to 800
c     .......... for i=en-2 step -1 until 1 do -- ..........
         do 795 ii = 1, enm2
            i = na - ii
            w = h(i,i) - p
            ra = 0.0d0
            sa = 0.0d0
c
            do 760 j = m, en
               ra = ra + h(i,j) * h(j,na)
               sa = sa + h(i,j) * h(j,en)
  760       continue
c
            if (wi(i) .ge. 0.0d0) go to 770
            zz = w
            r = ra
            s = sa
            go to 795
  770       m = i
            if (wi(i) .ne. 0.0d0) go to 780
            call cdiv(-ra,-sa,w,q,h(i,na),h(i,en))
            go to 790
c     .......... solve complex equations ..........
  780       x = h(i,i+1)
            y = h(i+1,i)
            vr = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i) - q * q
            vi = (wr(i) - p) * 2.0d0 * q
            if (vr .ne. 0.0d0 .or. vi .ne. 0.0d0) go to 784
               tst1 = norm * (dabs(w) + dabs(q) + dabs(x)
     x                      + dabs(y) + dabs(zz))
               vr = tst1
  783          vr = 0.01d0 * vr
               tst2 = tst1 + vr
               if (tst2 .gt. tst1) go to 783
  784       call cdiv(x*r-zz*ra+q*sa,x*s-zz*sa-q*ra,vr,vi,
     x                h(i,na),h(i,en))
            if (dabs(x) .le. dabs(zz) + dabs(q)) go to 785
            h(i+1,na) = (-ra - w * h(i,na) + q * h(i,en)) / x
            h(i+1,en) = (-sa - w * h(i,en) - q * h(i,na)) / x
            go to 790
  785       call cdiv(-r-y*h(i,na),-s-y*h(i,en),zz,q,
     x                h(i+1,na),h(i+1,en))
c
c     .......... overflow control ..........
  790       t = dmax1(dabs(h(i,na)), dabs(h(i,en)))
            if (t .eq. 0.0d0) go to 795
            tst1 = t
            tst2 = tst1 + 1.0d0/tst1
            if (tst2 .gt. tst1) go to 795
            do 792 j = i, en
               h(j,na) = h(j,na)/t
               h(j,en) = h(j,en)/t
  792       continue
c
  795    continue
c     .......... end complex vector ..........
  800 continue
c     .......... end back substitution.
c                vectors of isolated roots ..........
      do 840 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 840
c
         do 820 j = i, n
  820    z(i,j) = h(i,j)
c
  840 continue
c     .......... multiply by transformation matrix to give
c                vectors of original full matrix.
c                for j=n step -1 until low do -- ..........
      do 880 jj = low, n
         j = n + low - jj
         m = min0(j,igh)
c
         do 880 i = low, igh
            zz = 0.0d0
c
            do 860 k = low, m
  860       zz = zz + z(i,k) * h(k,j)
c
            z(i,j) = zz
  880 continue
c
      go to 1001
c     .......... set error -- all eigenvalues have not
c                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end
      subroutine htrib3(nm,n,a,tau,m,zr,zi)

c*********************************************************************72
c
cc HTRIB3 determines eigenvectors by undoing the HTRID3 transformation.
c
c     this subroutine is a translation of a complex analogue of
c     the algol procedure trbak3, num. math. 11, 181-195(1968)
c     by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine forms the eigenvectors of a complex hermitian
c     matrix by back transforming those of the corresponding
c     real symmetric tridiagonal matrix determined by  htrid3.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        a contains information about the unitary transformations
c          used in the reduction by  htrid3.
c
c        tau contains further information about the transformations.
c
c        m is the number of eigenvectors to be back transformed.
c
c        zr contains the eigenvectors to be back transformed
c          in its first m columns.
c
c     on output
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the transformed eigenvectors
c          in their first m columns.
c
c     note that the last component of each returned vector
c     is real and that vector euclidean norms are preserved.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,k,l,m,n,nm
      double precision a(nm,n),tau(2,n),zr(nm,m),zi(nm,m)
      double precision h,s,si

      if (m .eq. 0) go to 200
c     .......... transform the eigenvectors of the real symmetric
c                tridiagonal matrix to those of the hermitian
c                tridiagonal matrix. ..........
      do 50 k = 1, n
c
         do 50 j = 1, m
            zi(k,j) = -zr(k,j) * tau(2,k)
            zr(k,j) = zr(k,j) * tau(1,k)
   50 continue
c
      if (n .eq. 1) go to 200
c     .......... recover and apply the householder matrices ..........
      do 140 i = 2, n
         l = i - 1
         h = a(i,i)
         if (h .eq. 0.0d0) go to 140
c
         do 130 j = 1, m
            s = 0.0d0
            si = 0.0d0
c
            do 110 k = 1, l
               s = s + a(i,k) * zr(k,j) - a(k,i) * zi(k,j)
               si = si + a(i,k) * zi(k,j) + a(k,i) * zr(k,j)
  110       continue
c     .......... double divisions avoid possible underflow ..........
            s = (s / h) / h
            si = (si / h) / h
c
            do 120 k = 1, l
               zr(k,j) = zr(k,j) - s * a(i,k) - si * a(k,i)
               zi(k,j) = zi(k,j) - si * a(i,k) + s * a(k,i)
  120       continue
c
  130    continue
c
  140 continue
c
  200 return
      end
      subroutine htribk(nm,n,ar,ai,tau,m,zr,zi)

c*********************************************************************72
c
cc HTRIBK determines eigenvectors by undoing the HTRIDI transformation.
c
c     this subroutine is a translation of a complex analogue of
c     the algol procedure trbak1, num. math. 11, 181-195(1968)
c     by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine forms the eigenvectors of a complex hermitian
c     matrix by back transforming those of the corresponding
c     real symmetric tridiagonal matrix determined by  htridi.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        ar and ai contain information about the unitary trans-
c          formations used in the reduction by  htridi  in their
c          full lower triangles except for the diagonal of ar.
c
c        tau contains further information about the transformations.
c
c        m is the number of eigenvectors to be back transformed.
c
c        zr contains the eigenvectors to be back transformed
c          in its first m columns.
c
c     on output
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the transformed eigenvectors
c          in their first m columns.
c
c     note that the last component of each returned vector
c     is real and that vector euclidean norms are preserved.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,k,l,m,n,nm
      double precision ar(nm,n),ai(nm,n),tau(2,n),zr(nm,m),zi(nm,m)
      double precision h,s,si

      if (m .eq. 0) go to 200
c     .......... transform the eigenvectors of the real symmetric
c                tridiagonal matrix to those of the hermitian
c                tridiagonal matrix. ..........
      do 50 k = 1, n
c
         do 50 j = 1, m
            zi(k,j) = -zr(k,j) * tau(2,k)
            zr(k,j) = zr(k,j) * tau(1,k)
   50 continue
c
      if (n .eq. 1) go to 200
c     .......... recover and apply the householder matrices ..........
      do 140 i = 2, n
         l = i - 1
         h = ai(i,i)
         if (h .eq. 0.0d0) go to 140
c
         do 130 j = 1, m
            s = 0.0d0
            si = 0.0d0
c
            do 110 k = 1, l
               s = s + ar(i,k) * zr(k,j) - ai(i,k) * zi(k,j)
               si = si + ar(i,k) * zi(k,j) + ai(i,k) * zr(k,j)
  110       continue
c     .......... double divisions avoid possible underflow ..........
            s = (s / h) / h
            si = (si / h) / h
c
            do 120 k = 1, l
               zr(k,j) = zr(k,j) - s * ar(i,k) - si * ai(i,k)
               zi(k,j) = zi(k,j) - si * ar(i,k) + s * ai(i,k)
  120       continue
c
  130    continue
c
  140 continue
c
  200 return
      end
      subroutine htrid3(nm,n,a,d,e,e2,tau)

c*********************************************************************72
c
cc HTRID3 tridiagonalizes a complex hermitian packed matrix.
c
c     this subroutine is a translation of a complex analogue of
c     the algol procedure tred3, num. math. 11, 181-195(1968)
c     by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a complex hermitian matrix, stored as
c     a single square array, to a real symmetric tridiagonal matrix
c     using unitary similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        a contains the lower triangle of the complex hermitian input
c          matrix.  the real parts of the matrix elements are stored
c          in the full lower triangle of a, and the imaginary parts
c          are stored in the transposed positions of the strict upper
c          triangle of a.  no storage is required for the zero
c          imaginary parts of the diagonal elements.
c
c     on output
c
c        a contains information about the unitary transformations
c          used in the reduction.
c
c        d contains the diagonal elements of the the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        e2 contains the squares of the corresponding elements of e.
c          e2 may coincide with e if the squares are not needed.
c
c        tau contains further information about the transformations.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,k,l,n,ii,nm,jm1,jp1
      double precision a(nm,n),d(n),e(n),e2(n),tau(2,n)
      double precision f,g,h,fi,gi,hh,si,scale,pythag

      tau(1,n) = 1.0d0
      tau(2,n) = 0.0d0
c     .......... for i=n step -1 until 1 do -- ..........
      do 300 ii = 1, n
         i = n + 1 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 1) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + dabs(a(i,k)) + dabs(a(k,i))
c
         if (scale .ne. 0.0d0) go to 140
         tau(1,l) = 1.0d0
         tau(2,l) = 0.0d0
  130    e(i) = 0.0d0
         e2(i) = 0.0d0
         go to 290
c
  140    do 150 k = 1, l
            a(i,k) = a(i,k) / scale
            a(k,i) = a(k,i) / scale
            h = h + a(i,k) * a(i,k) + a(k,i) * a(k,i)
  150    continue
c
         e2(i) = scale * scale * h
         g = dsqrt(h)
         e(i) = scale * g
         f = pythag(a(i,l),a(l,i))
c     .......... form next diagonal element of matrix t ..........
         if (f .eq. 0.0d0) go to 160
         tau(1,l) = (a(l,i) * tau(2,i) - a(i,l) * tau(1,i)) / f
         si = (a(i,l) * tau(2,i) + a(l,i) * tau(1,i)) / f
         h = h + f * g
         g = 1.0d0 + g / f
         a(i,l) = g * a(i,l)
         a(l,i) = g * a(l,i)
         if (l .eq. 1) go to 270
         go to 170
  160    tau(1,l) = -tau(1,i)
         si = tau(2,i)
         a(i,l) = g
  170    f = 0.0d0
c
         do 240 j = 1, l
            g = 0.0d0
            gi = 0.0d0
            if (j .eq. 1) go to 190
            jm1 = j - 1
c     .......... form element of a*u ..........
            do 180 k = 1, jm1
               g = g + a(j,k) * a(i,k) + a(k,j) * a(k,i)
               gi = gi - a(j,k) * a(k,i) + a(k,j) * a(i,k)
  180       continue
c
  190       g = g + a(j,j) * a(i,j)
            gi = gi - a(j,j) * a(j,i)
            jp1 = j + 1
            if (l .lt. jp1) go to 220
c
            do 200 k = jp1, l
               g = g + a(k,j) * a(i,k) - a(j,k) * a(k,i)
               gi = gi - a(k,j) * a(k,i) - a(j,k) * a(i,k)
  200       continue
c     .......... form element of p ..........
  220       e(j) = g / h
            tau(2,j) = gi / h
            f = f + e(j) * a(i,j) - tau(2,j) * a(j,i)
  240    continue
c
         hh = f / (h + h)
c     .......... form reduced a ..........
         do 260 j = 1, l
            f = a(i,j)
            g = e(j) - hh * f
            e(j) = g
            fi = -a(j,i)
            gi = tau(2,j) - hh * fi
            tau(2,j) = -gi
            a(j,j) = a(j,j) - 2.0d0 * (f * g + fi * gi)
            if (j .eq. 1) go to 260
            jm1 = j - 1
c
            do 250 k = 1, jm1
               a(j,k) = a(j,k) - f * e(k) - g * a(i,k)
     x                         + fi * tau(2,k) + gi * a(k,i)
               a(k,j) = a(k,j) - f * tau(2,k) - g * a(k,i)
     x                         - fi * e(k) - gi * a(i,k)
  250       continue
c
  260    continue
c
  270    do 280 k = 1, l
            a(i,k) = scale * a(i,k)
            a(k,i) = scale * a(k,i)
  280    continue
c
         tau(2,l) = -si
  290    d(i) = a(i,i)
         a(i,i) = scale * dsqrt(h)
  300 continue
c
      return
      end
      subroutine htridi(nm,n,ar,ai,d,e,e2,tau)

c*********************************************************************72
c
cc HTRIDI tridiagonalizes a complex hermitian matrix.
c
c     this subroutine is a translation of a complex analogue of
c     the algol procedure tred1, num. math. 11, 181-195(1968)
c     by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a complex hermitian matrix
c     to a real symmetric tridiagonal matrix using
c     unitary similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the complex hermitian input matrix.
c          only the lower triangle of the matrix need be supplied.
c
c     on output
c
c        ar and ai contain information about the unitary trans-
c          formations used in the reduction in their full lower
c          triangles.  their strict upper triangles and the
c          diagonal of ar are unaltered.
c
c        d contains the diagonal elements of the the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        e2 contains the squares of the corresponding elements of e.
c          e2 may coincide with e if the squares are not needed.
c
c        tau contains further information about the transformations.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,k,l,n,ii,nm,jp1
      double precision ar(nm,n),ai(nm,n),d(n),e(n),e2(n),tau(2,n)
      double precision f,g,h,fi,gi,hh,si,scale,pythag

      tau(1,n) = 1.0d0
      tau(2,n) = 0.0d0
c
      do 100 i = 1, n
  100 d(i) = ar(i,i)
c     .......... for i=n step -1 until 1 do -- ..........
      do 300 ii = 1, n
         i = n + 1 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 1) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + dabs(ar(i,k)) + dabs(ai(i,k))
c
         if (scale .ne. 0.0d0) go to 140
         tau(1,l) = 1.0d0
         tau(2,l) = 0.0d0
  130    e(i) = 0.0d0
         e2(i) = 0.0d0
         go to 290
c
  140    do 150 k = 1, l
            ar(i,k) = ar(i,k) / scale
            ai(i,k) = ai(i,k) / scale
            h = h + ar(i,k) * ar(i,k) + ai(i,k) * ai(i,k)
  150    continue
c
         e2(i) = scale * scale * h
         g = dsqrt(h)
         e(i) = scale * g
         f = pythag(ar(i,l),ai(i,l))
c     .......... form next diagonal element of matrix t ..........
         if (f .eq. 0.0d0) go to 160
         tau(1,l) = (ai(i,l) * tau(2,i) - ar(i,l) * tau(1,i)) / f
         si = (ar(i,l) * tau(2,i) + ai(i,l) * tau(1,i)) / f
         h = h + f * g
         g = 1.0d0 + g / f
         ar(i,l) = g * ar(i,l)
         ai(i,l) = g * ai(i,l)
         if (l .eq. 1) go to 270
         go to 170
  160    tau(1,l) = -tau(1,i)
         si = tau(2,i)
         ar(i,l) = g
  170    f = 0.0d0
c
         do 240 j = 1, l
            g = 0.0d0
            gi = 0.0d0
c     .......... form element of a*u ..........
            do 180 k = 1, j
               g = g + ar(j,k) * ar(i,k) + ai(j,k) * ai(i,k)
               gi = gi - ar(j,k) * ai(i,k) + ai(j,k) * ar(i,k)
  180       continue
c
            jp1 = j + 1
            if (l .lt. jp1) go to 220
c
            do 200 k = jp1, l
               g = g + ar(k,j) * ar(i,k) - ai(k,j) * ai(i,k)
               gi = gi - ar(k,j) * ai(i,k) - ai(k,j) * ar(i,k)
  200       continue
c     .......... form element of p ..........
  220       e(j) = g / h
            tau(2,j) = gi / h
            f = f + e(j) * ar(i,j) - tau(2,j) * ai(i,j)
  240    continue
c
         hh = f / (h + h)
c     .......... form reduced a ..........
         do 260 j = 1, l
            f = ar(i,j)
            g = e(j) - hh * f
            e(j) = g
            fi = -ai(i,j)
            gi = tau(2,j) - hh * fi
            tau(2,j) = -gi
c
            do 260 k = 1, j
               ar(j,k) = ar(j,k) - f * e(k) - g * ar(i,k)
     x                           + fi * tau(2,k) + gi * ai(i,k)
               ai(j,k) = ai(j,k) - f * tau(2,k) - g * ai(i,k)
     x                           - fi * e(k) - gi * ar(i,k)
  260    continue
c
  270    do 280 k = 1, l
            ar(i,k) = scale * ar(i,k)
            ai(i,k) = scale * ai(i,k)
  280    continue
c
         tau(2,l) = -si
  290    hh = d(i)
         d(i) = ar(i,i)
         ar(i,i) = hh
         ai(i,i) = scale * dsqrt(h)
  300 continue
c
      return
      end
      subroutine imtql1(n,d,e,ierr)

c*********************************************************************72
c
cc IMTQL1 computes all eigenvalues of a symmetric tridiagonal matrix.
c
c     this subroutine is a translation of the algol procedure imtql1,
c     num. math. 12, 377-383(1968) by martin and wilkinson,
c     as modified in num. math. 15, 450(1970) by dubrulle.
c     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
c
c     this subroutine finds the eigenvalues of a symmetric
c     tridiagonal matrix by the implicit ql method.
c
c     on input
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct and
c          ordered for indices 1,2,...ierr-1, but may not be
c          the smallest eigenvalues.
c
c        e has been destroyed.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,l,m,n,ii,mml,ierr
      double precision d(n),e(n)
      double precision b,c,f,g,p,r,s,tst1,tst2,pythag

      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
c
      e(n) = 0.0d0
c
      do 290 l = 1, n
         j = 0
c     .......... look for small sub-diagonal element ..........
  105    do 110 m = l, n
            if (m .eq. n) go to 120
            tst1 = dabs(d(m)) + dabs(d(m+1))
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
  110    continue
c
  120    p = d(l)
         if (m .eq. l) go to 215
         if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         g = (d(l+1) - p) / (2.0d0 * e(l))
         r = pythag(g,1.0d0)
         g = d(m) - p + e(l) / (g + dsign(r,g))
         s = 1.0d0
         c = 1.0d0
         p = 0.0d0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            i = m - ii
            f = s * e(i)
            b = c * e(i)
            r = pythag(f,g)
            e(i+1) = r
            if (r .eq. 0.0d0) go to 210
            s = f / r
            c = g / r
            g = d(i+1) - p
            r = (d(i) - g) * s + 2.0d0 * c * b
            p = s * r
            d(i+1) = g + p
            g = c * r - b
  200    continue
c
         d(l) = d(l) - p
         e(l) = g
         e(m) = 0.0d0
         go to 105
c     .......... recover from underflow ..........
  210    d(i+1) = d(i+1) - p
         e(m) = 0.0d0
         go to 105
c     .......... order eigenvalues ..........
  215    if (l .eq. 1) go to 250
c     .......... for i=l step -1 until 2 do -- ..........
         do 230 ii = 2, l
            i = l + 2 - ii
            if (p .ge. d(i-1)) go to 270
            d(i) = d(i-1)
  230    continue
c
  250    i = 1
  270    d(i) = p
  290 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end
      subroutine imtql2(nm,n,d,e,z,ierr)

c*********************************************************************72
c
cc IMTQL2 computes all eigenvalues/vectors of a symmetric tridiagonal matrix.
c
c     this subroutine is a translation of the algol procedure imtql2,
c     num. math. 12, 377-383(1968) by martin and wilkinson,
c     as modified in num. math. 15, 450(1970) by dubrulle.
c     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a symmetric tridiagonal matrix by the implicit ql method.
c     the eigenvectors of a full symmetric matrix can also
c     be found if  tred2  has been used to reduce this
c     full matrix to tridiagonal form.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c        z contains the transformation matrix produced in the
c          reduction by  tred2, if performed.  if the eigenvectors
c          of the tridiagonal matrix are desired, z must contain
c          the identity matrix.
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct but
c          unordered for indices 1,2,...,ierr-1.
c
c        e has been destroyed.
c
c        z contains orthonormal eigenvectors of the symmetric
c          tridiagonal (or full) matrix.  if an error exit is made,
c          z contains the eigenvectors associated with the stored
c          eigenvalues.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,k,l,m,n,ii,nm,mml,ierr
      double precision d(n),e(n),z(nm,n)
      double precision b,c,f,g,p,r,s,tst1,tst2,pythag

      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
c
      e(n) = 0.0d0
c
      do 240 l = 1, n
         j = 0
c     .......... look for small sub-diagonal element ..........
  105    do 110 m = l, n
            if (m .eq. n) go to 120
            tst1 = dabs(d(m)) + dabs(d(m+1))
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
  110    continue
c
  120    p = d(l)
         if (m .eq. l) go to 240
         if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         g = (d(l+1) - p) / (2.0d0 * e(l))
         r = pythag(g,1.0d0)
         g = d(m) - p + e(l) / (g + dsign(r,g))
         s = 1.0d0
         c = 1.0d0
         p = 0.0d0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            i = m - ii
            f = s * e(i)
            b = c * e(i)
            r = pythag(f,g)
            e(i+1) = r
            if (r .eq. 0.0d0) go to 210
            s = f / r
            c = g / r
            g = d(i+1) - p
            r = (d(i) - g) * s + 2.0d0 * c * b
            p = s * r
            d(i+1) = g + p
            g = c * r - b
c     .......... form vector ..........
            do 180 k = 1, n
               f = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * f
               z(k,i) = c * z(k,i) - s * f
  180       continue
c
  200    continue
c
         d(l) = d(l) - p
         e(l) = g
         e(m) = 0.0d0
         go to 105
c     .......... recover from underflow ..........
  210    d(i+1) = d(i+1) - p
         e(m) = 0.0d0
         go to 105
  240 continue
c     .......... order eigenvalues and eigenvectors ..........
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
c
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
c
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
c
         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue
c
  300 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end
      subroutine imtqlv(n,d,e,e2,w,ind,ierr,rv1)

c*********************************************************************72
c
cc IMTQLV computes all eigenvalues of a real symmetric tridiagonal matrix.
c
c     this subroutine is a variant of  imtql1  which is a translation of
c     algol procedure imtql1, num. math. 12, 377-383(1968) by martin and
c     wilkinson, as modified in num. math. 15, 450(1970) by dubrulle.
c     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
c
c     this subroutine finds the eigenvalues of a symmetric tridiagonal
c     matrix by the implicit ql method and associates with them
c     their corresponding submatrix indices.
c
c     on input
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c        e2 contains the squares of the corresponding elements of e.
c          e2(1) is arbitrary.
c
c     on output
c
c        d and e are unaltered.
c
c        elements of e2, corresponding to elements of e regarded
c          as negligible, have been replaced by zero causing the
c          matrix to split into a direct sum of submatrices.
c          e2(1) is also set to zero.
c
c        w contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct and
c          ordered for indices 1,2,...ierr-1, but may not be
c          the smallest eigenvalues.
c
c        ind contains the submatrix indices associated with the
c          corresponding eigenvalues in w -- 1 for eigenvalues
c          belonging to the first submatrix from the top,
c          2 for those belonging to the second submatrix, etc..
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c        rv1 is a temporary storage array.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,k,l,m,n,ii,mml,tag,ierr
      double precision d(n),e(n),e2(n),w(n),rv1(n)
      double precision b,c,f,g,p,r,s,tst1,tst2,pythag
      integer ind(n)

      ierr = 0
      k = 0
      tag = 0
c
      do 100 i = 1, n
         w(i) = d(i)
         if (i .ne. 1) rv1(i-1) = e(i)
  100 continue
c
      e2(1) = 0.0d0
      rv1(n) = 0.0d0
c
      do 290 l = 1, n
         j = 0
c     .......... look for small sub-diagonal element ..........
  105    do 110 m = l, n
            if (m .eq. n) go to 120
            tst1 = dabs(w(m)) + dabs(w(m+1))
            tst2 = tst1 + dabs(rv1(m))
            if (tst2 .eq. tst1) go to 120
c     .......... guard against underflowed element of e2 ..........
            if (e2(m+1) .eq. 0.0d0) go to 125
  110    continue
c
  120    if (m .le. k) go to 130
         if (m .ne. n) e2(m+1) = 0.0d0
  125    k = m
         tag = tag + 1
  130    p = w(l)
         if (m .eq. l) go to 215
         if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         g = (w(l+1) - p) / (2.0d0 * rv1(l))
         r = pythag(g,1.0d0)
         g = w(m) - p + rv1(l) / (g + dsign(r,g))
         s = 1.0d0
         c = 1.0d0
         p = 0.0d0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            i = m - ii
            f = s * rv1(i)
            b = c * rv1(i)
            r = pythag(f,g)
            rv1(i+1) = r
            if (r .eq. 0.0d0) go to 210
            s = f / r
            c = g / r
            g = w(i+1) - p
            r = (w(i) - g) * s + 2.0d0 * c * b
            p = s * r
            w(i+1) = g + p
            g = c * r - b
  200    continue
c
         w(l) = w(l) - p
         rv1(l) = g
         rv1(m) = 0.0d0
         go to 105
c     .......... recover from underflow ..........
  210    w(i+1) = w(i+1) - p
         rv1(m) = 0.0d0
         go to 105
c     .......... order eigenvalues ..........
  215    if (l .eq. 1) go to 250
c     .......... for i=l step -1 until 2 do -- ..........
         do 230 ii = 2, l
            i = l + 2 - ii
            if (p .ge. w(i-1)) go to 270
            w(i) = w(i-1)
            ind(i) = ind(i-1)
  230    continue
c
  250    i = 1
  270    w(i) = p
         ind(i) = tag
  290 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end
      subroutine invit(nm,n,a,wr,wi,select,mm,m,z,ierr,rm1,rv1,rv2)

c*********************************************************************72
c
cc INVIT computes eigenvectors given eigenvalues, for a real upper Hessenberg matrix.
c
c     this subroutine is a translation of the algol procedure invit
c     by peters and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 418-439(1971).
c
c     this subroutine finds those eigenvectors of a real upper
c     hessenberg matrix corresponding to specified eigenvalues,
c     using inverse iteration.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        a contains the hessenberg matrix.
c
c        wr and wi contain the real and imaginary parts, respectively,
c          of the eigenvalues of the matrix.  the eigenvalues must be
c          stored in a manner identical to that of subroutine  hqr,
c          which recognizes possible splitting of the matrix.
c
c        select specifies the eigenvectors to be found. the
c          eigenvector corresponding to the j-th eigenvalue is
c          specified by setting select(j) to .true..
c
c        mm should be set to an upper bound for the number of
c          columns required to store the eigenvectors to be found.
c          note that two columns are required to store the
c          eigenvector corresponding to a complex eigenvalue.
c
c     on output
c
c        a and wi are unaltered.
c
c        wr may have been altered since close eigenvalues are perturbed
c          slightly in searching for independent eigenvectors.
c
c        select may have been altered.  if the elements corresponding
c          to a pair of conjugate complex eigenvalues were each
c          initially set to .true., the program resets the second of
c          the two elements to .false..
c
c        m is the number of columns actually used to store
c          the eigenvectors.
c
c        z contains the real and imaginary parts of the eigenvectors.
c          if the next selected eigenvalue is real, the next column
c          of z contains its eigenvector.  if the eigenvalue is
c          complex, the next two columns of z contain the real and
c          imaginary parts of its eigenvector.  the eigenvectors are
c          normalized so that the component of largest magnitude is 1.
c          any vector which fails the acceptance test is set to zero.
c
c        ierr is set to
c          zero       for normal return,
c          -(2*n+1)   if more than mm columns of z are necessary
c                     to store the eigenvectors corresponding to
c                     the specified eigenvalues.
c          -k         if the iteration corresponding to the k-th
c                     value fails,
c          -(n+k)     if both error situations occur.
c
c        rm1, rv1, and rv2 are temporary storage arrays.  note that rm1
c          is square of dimension n by n and, augmented by two columns
c          of z, is the transpose of the corresponding algol b array.
c
c     the algol procedure guessvec appears in invit in line.
c
c     calls cdiv for complex division.
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,k,l,m,n,s,ii,ip,mm,mp,nm,ns,n1,uk,ip1,its,km1,ierr
      double precision a(nm,n),wr(n),wi(n),z(nm,mm),rm1(n,n),
     x       rv1(n),rv2(n)
      double precision t,w,x,y,eps3,norm,normv,epslon,growto,ilambd,
     x       pythag,rlambd,ukroot
      logical select(n)

      ierr = 0
      uk = 0
      s = 1
c     .......... ip = 0, real eigenvalue
c                     1, first of conjugate complex pair
c                    -1, second of conjugate complex pair ..........
      ip = 0
      n1 = n - 1
c
      do 980 k = 1, n
         if (wi(k) .eq. 0.0d0 .or. ip .lt. 0) go to 100
         ip = 1
         if (select(k) .and. select(k+1)) select(k+1) = .false.
  100    if (.not. select(k)) go to 960
         if (wi(k) .ne. 0.0d0) s = s + 1
         if (s .gt. mm) go to 1000
         if (uk .ge. k) go to 200
c     .......... check for possible splitting ..........
         do 120 uk = k, n
            if (uk .eq. n) go to 140
            if (a(uk+1,uk) .eq. 0.0d0) go to 140
  120    continue
c     .......... compute infinity norm of leading uk by uk
c                (hessenberg) matrix ..........
  140    norm = 0.0d0
         mp = 1
c
         do 180 i = 1, uk
            x = 0.0d0
c
            do 160 j = mp, uk
  160       x = x + dabs(a(i,j))
c
            if (x .gt. norm) norm = x
            mp = i
  180    continue
c     .......... eps3 replaces zero pivot in decomposition
c                and close roots are modified by eps3 ..........
         if (norm .eq. 0.0d0) norm = 1.0d0
         eps3 = epslon(norm)
c     .......... growto is the criterion for the growth ..........
         ukroot = uk
         ukroot = dsqrt(ukroot)
         growto = 0.1d0 / ukroot
  200    rlambd = wr(k)
         ilambd = wi(k)
         if (k .eq. 1) go to 280
         km1 = k - 1
         go to 240
c     .......... perturb eigenvalue if it is close
c                to any previous eigenvalue ..........
  220    rlambd = rlambd + eps3
c     .......... for i=k-1 step -1 until 1 do -- ..........
  240    do 260 ii = 1, km1
            i = k - ii
            if (select(i) .and. dabs(wr(i)-rlambd) .lt. eps3 .and.
     x         dabs(wi(i)-ilambd) .lt. eps3) go to 220
  260    continue
c
         wr(k) = rlambd
c     .......... perturb conjugate eigenvalue to match ..........
         ip1 = k + ip
         wr(ip1) = rlambd
c     .......... form upper hessenberg a-rlambd*i (transposed)
c                and initial real vector ..........
  280    mp = 1
c
         do 320 i = 1, uk
c
            do 300 j = mp, uk
  300       rm1(j,i) = a(i,j)
c
            rm1(i,i) = rm1(i,i) - rlambd
            mp = i
            rv1(i) = eps3
  320    continue
c
         its = 0
         if (ilambd .ne. 0.0d0) go to 520
c     .......... real eigenvalue.
c                triangular decomposition with interchanges,
c                replacing zero pivots by eps3 ..........
         if (uk .eq. 1) go to 420
c
         do 400 i = 2, uk
            mp = i - 1
            if (dabs(rm1(mp,i)) .le. dabs(rm1(mp,mp))) go to 360
c
            do 340 j = mp, uk
               y = rm1(j,i)
               rm1(j,i) = rm1(j,mp)
               rm1(j,mp) = y
  340       continue
c
  360       if (rm1(mp,mp) .eq. 0.0d0) rm1(mp,mp) = eps3
            x = rm1(mp,i) / rm1(mp,mp)
            if (x .eq. 0.0d0) go to 400
c
            do 380 j = i, uk
  380       rm1(j,i) = rm1(j,i) - x * rm1(j,mp)
c
  400    continue
c
  420    if (rm1(uk,uk) .eq. 0.0d0) rm1(uk,uk) = eps3
c     .......... back substitution for real vector
c                for i=uk step -1 until 1 do -- ..........
  440    do 500 ii = 1, uk
            i = uk + 1 - ii
            y = rv1(i)
            if (i .eq. uk) go to 480
            ip1 = i + 1
c
            do 460 j = ip1, uk
  460       y = y - rm1(j,i) * rv1(j)
c
  480       rv1(i) = y / rm1(i,i)
  500    continue
c
         go to 740
c     .......... complex eigenvalue.
c                triangular decomposition with interchanges,
c                replacing zero pivots by eps3.  store imaginary
c                parts in upper triangle starting at (1,3) ..........
  520    ns = n - s
         z(1,s-1) = -ilambd
         z(1,s) = 0.0d0
         if (n .eq. 2) go to 550
         rm1(1,3) = -ilambd
         z(1,s-1) = 0.0d0
         if (n .eq. 3) go to 550
c
         do 540 i = 4, n
  540    rm1(1,i) = 0.0d0
c
  550    do 640 i = 2, uk
            mp = i - 1
            w = rm1(mp,i)
            if (i .lt. n) t = rm1(mp,i+1)
            if (i .eq. n) t = z(mp,s-1)
            x = rm1(mp,mp) * rm1(mp,mp) + t * t
            if (w * w .le. x) go to 580
            x = rm1(mp,mp) / w
            y = t / w
            rm1(mp,mp) = w
            if (i .lt. n) rm1(mp,i+1) = 0.0d0
            if (i .eq. n) z(mp,s-1) = 0.0d0
c
            do 560 j = i, uk
               w = rm1(j,i)
               rm1(j,i) = rm1(j,mp) - x * w
               rm1(j,mp) = w
               if (j .lt. n1) go to 555
               l = j - ns
               z(i,l) = z(mp,l) - y * w
               z(mp,l) = 0.0d0
               go to 560
  555          rm1(i,j+2) = rm1(mp,j+2) - y * w
               rm1(mp,j+2) = 0.0d0
  560       continue
c
            rm1(i,i) = rm1(i,i) - y * ilambd
            if (i .lt. n1) go to 570
            l = i - ns
            z(mp,l) = -ilambd
            z(i,l) = z(i,l) + x * ilambd
            go to 640
  570       rm1(mp,i+2) = -ilambd
            rm1(i,i+2) = rm1(i,i+2) + x * ilambd
            go to 640
  580       if (x .ne. 0.0d0) go to 600
            rm1(mp,mp) = eps3
            if (i .lt. n) rm1(mp,i+1) = 0.0d0
            if (i .eq. n) z(mp,s-1) = 0.0d0
            t = 0.0d0
            x = eps3 * eps3
  600       w = w / x
            x = rm1(mp,mp) * w
            y = -t * w
c
            do 620 j = i, uk
               if (j .lt. n1) go to 610
               l = j - ns
               t = z(mp,l)
               z(i,l) = -x * t - y * rm1(j,mp)
               go to 615
  610          t = rm1(mp,j+2)
               rm1(i,j+2) = -x * t - y * rm1(j,mp)
  615          rm1(j,i) = rm1(j,i) - x * rm1(j,mp) + y * t
  620       continue
c
            if (i .lt. n1) go to 630
            l = i - ns
            z(i,l) = z(i,l) - ilambd
            go to 640
  630       rm1(i,i+2) = rm1(i,i+2) - ilambd
  640    continue
c
         if (uk .lt. n1) go to 650
         l = uk - ns
         t = z(uk,l)
         go to 655
  650    t = rm1(uk,uk+2)
  655    if (rm1(uk,uk) .eq. 0.0d0 .and. t .eq. 0.0d0) rm1(uk,uk) = eps3
c     .......... back substitution for complex vector
c                for i=uk step -1 until 1 do -- ..........
  660    do 720 ii = 1, uk
            i = uk + 1 - ii
            x = rv1(i)
            y = 0.0d0
            if (i .eq. uk) go to 700
            ip1 = i + 1
c
            do 680 j = ip1, uk
               if (j .lt. n1) go to 670
               l = j - ns
               t = z(i,l)
               go to 675
  670          t = rm1(i,j+2)
  675          x = x - rm1(j,i) * rv1(j) + t * rv2(j)
               y = y - rm1(j,i) * rv2(j) - t * rv1(j)
  680       continue
c
  700       if (i .lt. n1) go to 710
            l = i - ns
            t = z(i,l)
            go to 715
  710       t = rm1(i,i+2)
  715       call cdiv(x,y,rm1(i,i),t,rv1(i),rv2(i))
  720    continue
c     .......... acceptance test for real or complex
c                eigenvector and normalization ..........
  740    its = its + 1
         norm = 0.0d0
         normv = 0.0d0
c
         do 780 i = 1, uk
            if (ilambd .eq. 0.0d0) x = dabs(rv1(i))
            if (ilambd .ne. 0.0d0) x = pythag(rv1(i),rv2(i))
            if (normv .ge. x) go to 760
            normv = x
            j = i
  760       norm = norm + x
  780    continue
c
         if (norm .lt. growto) go to 840
c     .......... accept vector ..........
         x = rv1(j)
         if (ilambd .eq. 0.0d0) x = 1.0d0 / x
         if (ilambd .ne. 0.0d0) y = rv2(j)
c
         do 820 i = 1, uk
            if (ilambd .ne. 0.0d0) go to 800
            z(i,s) = rv1(i) * x
            go to 820
  800       call cdiv(rv1(i),rv2(i),x,y,z(i,s-1),z(i,s))
  820    continue
c
         if (uk .eq. n) go to 940
         j = uk + 1
         go to 900
c     .......... in-line procedure for choosing
c                a new starting vector ..........
  840    if (its .ge. uk) go to 880
         x = ukroot
         y = eps3 / (x + 1.0d0)
         rv1(1) = eps3
c
         do 860 i = 2, uk
  860    rv1(i) = y
c
         j = uk - its + 1
         rv1(j) = rv1(j) - eps3 * x
         if (ilambd .eq. 0.0d0) go to 440
         go to 660
c     .......... set error -- unaccepted eigenvector ..........
  880    j = 1
         ierr = -k
c     .......... set remaining vector components to zero ..........
  900    do 920 i = j, n
            z(i,s) = 0.0d0
            if (ilambd .ne. 0.0d0) z(i,s-1) = 0.0d0
  920    continue
c
  940    s = s + 1
  960    if (ip .eq. (-1)) ip = 0
         if (ip .eq. 1) ip = -1
  980 continue
c
      go to 1001
c     .......... set error -- underestimate of eigenvector
c                space required ..........
 1000 if (ierr .ne. 0) ierr = ierr - n
      if (ierr .eq. 0) ierr = -(2 * n + 1)
 1001 m = s - 1 - iabs(ip)
      return
      end
      subroutine minfit(nm,m,n,a,w,ip,b,ierr,rv1)

c*********************************************************************72
c
cc MINFIT solves the least squares problem, for a real overdetermined linear system.
c
c     this subroutine is a translation of the algol procedure minfit,
c     num. math. 14, 403-420(1970) by golub and reinsch.
c     handbook for auto. comp., vol ii-linear algebra, 134-151(1971).
c
c     this subroutine determines, towards the solution of the linear
c                                                        t
c     system ax=b, the singular value decomposition a=usv  of a real
c                                         t
c     m by n rectangular matrix, forming u b rather than u.  householder
c     bidiagonalization and a variant of the qr algorithm are used.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.  note that nm must be at least
c          as large as the maximum of m and n.
c
c        m is the number of rows of a and b.
c
c        n is the number of columns of a and the order of v.
c
c        a contains the rectangular coefficient matrix of the system.
c
c        ip is the number of columns of b.  ip can be zero.
c
c        b contains the constant column matrix of the system
c          if ip is not zero.  otherwise b is not referenced.
c
c     on output
c
c        a has been overwritten by the matrix v (orthogonal) of the
c          decomposition in its first n rows and columns.  if an
c          error exit is made, the columns of v corresponding to
c          indices of correct singular values should be correct.
c
c        w contains the n (non-negative) singular values of a (the
c          diagonal elements of s).  they are unordered.  if an
c          error exit is made, the singular values should be correct
c          for indices ierr+1,ierr+2,...,n.
c
c                                   t
c        b has been overwritten by u b.  if an error exit is made,
c                       t
c          the rows of u b corresponding to indices of correct
c          singular values should be correct.
c
c        ierr is set to
c          zero       for normal return,
c          k          if the k-th singular value has not been
c                     determined after 30 iterations.
c
c        rv1 is a temporary storage array.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,k,l,m,n,ii,ip,i1,kk,k1,ll,l1,m1,nm,its,ierr
      double precision a(nm,n),w(n),b(nm,ip),rv1(n)
      double precision c,f,g,h,s,x,y,z,tst1,tst2,scale,pythag

      ierr = 0
c     .......... householder reduction to bidiagonal form ..........
      g = 0.0d0
      scale = 0.0d0
      x = 0.0d0
c
      do 300 i = 1, n
         l = i + 1
         rv1(i) = scale * g
         g = 0.0d0
         s = 0.0d0
         scale = 0.0d0
         if (i .gt. m) go to 210
c
         do 120 k = i, m
  120    scale = scale + dabs(a(k,i))
c
         if (scale .eq. 0.0d0) go to 210
c
         do 130 k = i, m
            a(k,i) = a(k,i) / scale
            s = s + a(k,i)**2
  130    continue
c
         f = a(i,i)
         g = -dsign(dsqrt(s),f)
         h = f * g - s
         a(i,i) = f - g
         if (i .eq. n) go to 160
c
         do 150 j = l, n
            s = 0.0d0
c
            do 140 k = i, m
  140       s = s + a(k,i) * a(k,j)
c
            f = s / h
c
            do 150 k = i, m
               a(k,j) = a(k,j) + f * a(k,i)
  150    continue
c
  160    if (ip .eq. 0) go to 190
c
         do 180 j = 1, ip
            s = 0.0d0
c
            do 170 k = i, m
  170       s = s + a(k,i) * b(k,j)
c
            f = s / h
c
            do 180 k = i, m
               b(k,j) = b(k,j) + f * a(k,i)
  180    continue
c
  190    do 200 k = i, m
  200    a(k,i) = scale * a(k,i)
c
  210    w(i) = scale * g
         g = 0.0d0
         s = 0.0d0
         scale = 0.0d0
         if (i .gt. m .or. i .eq. n) go to 290
c
         do 220 k = l, n
  220    scale = scale + dabs(a(i,k))
c
         if (scale .eq. 0.0d0) go to 290
c
         do 230 k = l, n
            a(i,k) = a(i,k) / scale
            s = s + a(i,k)**2
  230    continue
c
         f = a(i,l)
         g = -dsign(dsqrt(s),f)
         h = f * g - s
         a(i,l) = f - g
c
         do 240 k = l, n
  240    rv1(k) = a(i,k) / h
c
         if (i .eq. m) go to 270
c
         do 260 j = l, m
            s = 0.0d0
c
            do 250 k = l, n
  250       s = s + a(j,k) * a(i,k)
c
            do 260 k = l, n
               a(j,k) = a(j,k) + s * rv1(k)
  260    continue
c
  270    do 280 k = l, n
  280    a(i,k) = scale * a(i,k)
c
  290    x = dmax1(x,dabs(w(i))+dabs(rv1(i)))
  300 continue
c     .......... accumulation of right-hand transformations.
c                for i=n step -1 until 1 do -- ..........
      do 400 ii = 1, n
         i = n + 1 - ii
         if (i .eq. n) go to 390
         if (g .eq. 0.0d0) go to 360
c
         do 320 j = l, n
c     .......... double division avoids possible underflow ..........
  320    a(j,i) = (a(i,j) / a(i,l)) / g
c
         do 350 j = l, n
            s = 0.0d0
c
            do 340 k = l, n
  340       s = s + a(i,k) * a(k,j)
c
            do 350 k = l, n
               a(k,j) = a(k,j) + s * a(k,i)
  350    continue
c
  360    do 380 j = l, n
            a(i,j) = 0.0d0
            a(j,i) = 0.0d0
  380    continue
c
  390    a(i,i) = 1.0d0
         g = rv1(i)
         l = i
  400 continue
c
      if (m .ge. n .or. ip .eq. 0) go to 510
      m1 = m + 1
c
      do 500 i = m1, n
c
         do 500 j = 1, ip
            b(i,j) = 0.0d0
  500 continue
c     .......... diagonalization of the bidiagonal form ..........
  510 tst1 = x
c     .......... for k=n step -1 until 1 do -- ..........
      do 700 kk = 1, n
         k1 = n - kk
         k = k1 + 1
         its = 0
c     .......... test for splitting.
c                for l=k step -1 until 1 do -- ..........
  520    do 530 ll = 1, k
            l1 = k - ll
            l = l1 + 1
            tst2 = tst1 + dabs(rv1(l))
            if (tst2 .eq. tst1) go to 565
c     .......... rv1(1) is always zero, so there is no exit
c                through the bottom of the loop ..........
            tst2 = tst1 + dabs(w(l1))
            if (tst2 .eq. tst1) go to 540
  530    continue
c     .......... cancellation of rv1(l) if l greater than 1 ..........
  540    c = 0.0d0
         s = 1.0d0
c
         do 560 i = l, k
            f = s * rv1(i)
            rv1(i) = c * rv1(i)
            tst2 = tst1 + dabs(f)
            if (tst2 .eq. tst1) go to 565
            g = w(i)
            h = pythag(f,g)
            w(i) = h
            c = g / h
            s = -f / h
            if (ip .eq. 0) go to 560
c
            do 550 j = 1, ip
               y = b(l1,j)
               z = b(i,j)
               b(l1,j) = y * c + z * s
               b(i,j) = -y * s + z * c
  550       continue
c
  560    continue
c     .......... test for convergence ..........
  565    z = w(k)
         if (l .eq. k) go to 650
c     .......... shift from bottom 2 by 2 minor ..........
         if (its .eq. 30) go to 1000
         its = its + 1
         x = w(l)
         y = w(k1)
         g = rv1(k1)
         h = rv1(k)
         f = 0.5d0 * (((g + z) / h) * ((g - z) / y) + y / h - h / y)
         g = pythag(f,1.0d0)
         f = x - (z / x) * z + (h / x) * (y / (f + dsign(g,f)) - h)
c     .......... next qr transformation ..........
         c = 1.0d0
         s = 1.0d0
c
         do 600 i1 = l, k1
            i = i1 + 1
            g = rv1(i)
            y = w(i)
            h = s * g
            g = c * g
            z = pythag(f,h)
            rv1(i1) = z
            c = f / z
            s = h / z
            f = x * c + g * s
            g = -x * s + g * c
            h = y * s
            y = y * c
c
            do 570 j = 1, n
               x = a(j,i1)
               z = a(j,i)
               a(j,i1) = x * c + z * s
               a(j,i) = -x * s + z * c
  570       continue
c
            z = pythag(f,h)
            w(i1) = z
c     .......... rotation can be arbitrary if z is zero ..........
            if (z .eq. 0.0d0) go to 580
            c = f / z
            s = h / z
  580       f = c * g + s * y
            x = -s * g + c * y
            if (ip .eq. 0) go to 600
c
            do 590 j = 1, ip
               y = b(i1,j)
               z = b(i,j)
               b(i1,j) = y * c + z * s
               b(i,j) = -y * s + z * c
  590       continue
c
  600    continue
c
         rv1(l) = 0.0d0
         rv1(k) = f
         w(k) = x
         go to 520
c     .......... convergence ..........
  650    if (z .ge. 0.0d0) go to 700
c     .......... w(k) is made non-negative ..........
         w(k) = -z
c
         do 690 j = 1, n
  690    a(j,k) = -a(j,k)
c
  700 continue
c
      go to 1001
c     .......... set error -- no convergence to a
c                singular value after 30 iterations ..........
 1000 ierr = k
 1001 return
      end
      subroutine ortbak(nm,low,igh,a,ort,m,z)

c*********************************************************************72
c
cc ORTBAK determines eigenvectors by undoing the ORTHES transformation.
c
c     this subroutine is a translation of the algol procedure ortbak,
c     num. math. 12, 349-368(1968) by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
c
c     this subroutine forms the eigenvectors of a real general
c     matrix by back transforming those of the corresponding
c     upper hessenberg matrix determined by  orthes.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        low and igh are integers determined by the balancing
c        routine  balanc.  if  balanc  has not been used,
c          set low=1 and igh equal to the order of the matrix.
c
c        a contains information about the orthogonal trans-
c          formations used in the reduction by  orthes
c          in its strict lower triangle.
c
c        ort contains further information about the trans-
c          formations used in the reduction by  orthes.
c          only elements low through igh are used.
c
c        m is the number of columns of z to be back transformed.
c
c        z contains the real and imaginary parts of the eigen-
c          vectors to be back transformed in its first m columns.
c
c     on output
c
c        z contains the real and imaginary parts of the
c          transformed eigenvectors in its first m columns.
c
c        ort has been altered.
c
c     note that ortbak preserves vector euclidean norms.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,m,la,mm,mp,nm,igh,kp1,low,mp1
      double precision a(nm,igh),ort(igh),z(nm,m)
      double precision g

      if (m .eq. 0) go to 200
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
c     .......... for mp=igh-1 step -1 until low+1 do -- ..........
      do 140 mm = kp1, la
         mp = low + igh - mm
         if (a(mp,mp-1) .eq. 0.0d0) go to 140
         mp1 = mp + 1
c
         do 100 i = mp1, igh
  100    ort(i) = a(i,mp-1)
c
         do 130 j = 1, m
            g = 0.0d0
c
            do 110 i = mp, igh
  110       g = g + ort(i) * z(i,j)
c     .......... divisor below is negative of h formed in orthes.
c                double division avoids possible underflow ..........
            g = (g / ort(mp)) / a(mp,mp-1)
c
            do 120 i = mp, igh
  120       z(i,j) = z(i,j) + g * ort(i)
c
  130    continue
c
  140 continue
c
  200 return
      end
      subroutine orthes(nm,n,low,igh,a,ort)

c*********************************************************************72
c
cc ORTHES transforms a real general matrix to upper Hessenberg form.
c
c     this subroutine is a translation of the algol procedure orthes,
c     num. math. 12, 349-368(1968) by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
c
c     given a real general matrix, this subroutine
c     reduces a submatrix situated in rows and columns
c     low through igh to upper hessenberg form by
c     orthogonal similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c        routine  balanc.  if  balanc  has not been used,
c          set low=1, igh=n.
c
c        a contains the input matrix.
c
c     on output
c
c        a contains the hessenberg matrix.  information about
c          the orthogonal transformations used in the reduction
c          is stored in the remaining triangle under the
c          hessenberg matrix.
c
c        ort contains further information about the transformations.
c          only elements low through igh are used.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,m,n,ii,jj,la,mp,nm,igh,kp1,low
      double precision a(nm,n),ort(igh)
      double precision f,g,h,scale

      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
c
      do 180 m = kp1, la
         h = 0.0d0
         ort(m) = 0.0d0
         scale = 0.0d0
c     .......... scale column (algol tol then not needed) ..........
         do 90 i = m, igh
   90    scale = scale + dabs(a(i,m-1))
c
         if (scale .eq. 0.0d0) go to 180
         mp = m + igh
c     .......... for i=igh step -1 until m do -- ..........
         do 100 ii = m, igh
            i = mp - ii
            ort(i) = a(i,m-1) / scale
            h = h + ort(i) * ort(i)
  100    continue
c
         g = -dsign(dsqrt(h),ort(m))
         h = h - ort(m) * g
         ort(m) = ort(m) - g
c     .......... form (i-(u*ut)/h) * a ..........
         do 130 j = m, n
            f = 0.0d0
c     .......... for i=igh step -1 until m do -- ..........
            do 110 ii = m, igh
               i = mp - ii
               f = f + ort(i) * a(i,j)
  110       continue
c
            f = f / h
c
            do 120 i = m, igh
  120       a(i,j) = a(i,j) - f * ort(i)
c
  130    continue
c     .......... form (i-(u*ut)/h)*a*(i-(u*ut)/h) ..........
         do 160 i = 1, igh
            f = 0.0d0
c     .......... for j=igh step -1 until m do -- ..........
            do 140 jj = m, igh
               j = mp - jj
               f = f + ort(j) * a(i,j)
  140       continue
c
            f = f / h
c
            do 150 j = m, igh
  150       a(i,j) = a(i,j) - f * ort(j)
c
  160    continue
c
         ort(m) = scale * ort(m)
         a(m,m-1) = scale * g
  180 continue
c
  200 return
      end
      subroutine ortran(nm,n,low,igh,a,ort,z)

c*********************************************************************72
c
cc ORTRAN accumulates similarity transformations generated by ORTHES.
c
c     this subroutine is a translation of the algol procedure ortrans,
c     num. math. 16, 181-204(1970) by peters and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
c
c     this subroutine accumulates the orthogonal similarity
c     transformations used in the reduction of a real general
c     matrix to upper hessenberg form by  orthes.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c        routine  balanc.  if  balanc  has not been used,
c          set low=1, igh=n.
c
c        a contains information about the orthogonal trans-
c          formations used in the reduction by  orthes
c          in its strict lower triangle.
c
c        ort contains further information about the trans-
c          formations used in the reduction by  orthes.
c          only elements low through igh are used.
c
c     on output
c
c        z contains the transformation matrix produced in the
c          reduction by  orthes.
c
c        ort has been altered.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,n,kl,mm,mp,nm,igh,low,mp1
      double precision a(nm,igh),ort(igh),z(nm,n)
      double precision g

c     .......... initialize z to identity matrix ..........
      do 80 j = 1, n
c
         do 60 i = 1, n
   60    z(i,j) = 0.0d0
c
         z(j,j) = 1.0d0
   80 continue
c
      kl = igh - low - 1
      if (kl .lt. 1) go to 200
c     .......... for mp=igh-1 step -1 until low+1 do -- ..........
      do 140 mm = 1, kl
         mp = igh - mm
         if (a(mp,mp-1) .eq. 0.0d0) go to 140
         mp1 = mp + 1
c
         do 100 i = mp1, igh
  100    ort(i) = a(i,mp-1)
c
         do 130 j = mp, igh
            g = 0.0d0
c
            do 110 i = mp, igh
  110       g = g + ort(i) * z(i,j)
c     .......... divisor below is negative of h formed in orthes.
c                double division avoids possible underflow ..........
            g = (g / ort(mp)) / a(mp,mp-1)
c
            do 120 i = mp, igh
  120       z(i,j) = z(i,j) + g * ort(i)
c
  130    continue
c
  140 continue
c
  200 return
      end
      double precision function pythag(a,b)

c*********************************************************************72
c
cc PYTHAG computes SQRT ( A**2 + B**2 ) carefully.
c
c     finds dsqrt(a**2+b**2) without overflow or destructive underflow
c
      double precision a,b
      double precision p,r,s,t,u

      p = dmax1(dabs(a),dabs(b))
      if (p .eq. 0.0d0) go to 20
      r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end
      subroutine qzhes(nm,n,a,b,matz,z)

c*********************************************************************72
c
cc QZHES carries out transformations for a generalized eigenvalue problem.
c
c     this subroutine is the first step of the qz algorithm
c     for solving generalized matrix eigenvalue problems,
c     siam j. numer. anal. 10, 241-256(1973) by moler and stewart.
c
c     this subroutine accepts a pair of real general matrices and
c     reduces one of them to upper hessenberg form and the other
c     to upper triangular form using orthogonal transformations.
c     it is usually followed by  qzit,  qzval  and, possibly,  qzvec.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrices.
c
c        a contains a real general matrix.
c
c        b contains a real general matrix.
c
c        matz should be set to .true. if the right hand transformations
c          are to be accumulated for later use in computing
c          eigenvectors, and to .false. otherwise.
c
c     on output
c
c        a has been reduced to upper hessenberg form.  the elements
c          below the first subdiagonal have been set to zero.
c
c        b has been reduced to upper triangular form.  the elements
c          below the main diagonal have been set to zero.
c
c        z contains the product of the right hand transformations if
c          matz has been set to .true.  otherwise, z is not referenced.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,k,l,n,lb,l1,nm,nk1,nm1,nm2
      double precision a(nm,n),b(nm,n),z(nm,n)
      double precision r,s,t,u1,u2,v1,v2,rho
      logical matz
c
c     .......... initialize z ..........
      if (.not. matz) go to 10
c
      do 3 j = 1, n
c
         do 2 i = 1, n
            z(i,j) = 0.0d0
    2    continue
c
         z(j,j) = 1.0d0
    3 continue
c     .......... reduce b to upper triangular form ..........
   10 if (n .le. 1) go to 170
      nm1 = n - 1
c
      do 100 l = 1, nm1
         l1 = l + 1
         s = 0.0d0
c
         do 20 i = l1, n
            s = s + dabs(b(i,l))
   20    continue
c
         if (s .eq. 0.0d0) go to 100
         s = s + dabs(b(l,l))
         r = 0.0d0
c
         do 25 i = l, n
            b(i,l) = b(i,l) / s
            r = r + b(i,l)**2
   25    continue
c
         r = dsign(dsqrt(r),b(l,l))
         b(l,l) = b(l,l) + r
         rho = r * b(l,l)
c
         do 50 j = l1, n
            t = 0.0d0
c
            do 30 i = l, n
               t = t + b(i,l) * b(i,j)
   30       continue
c
            t = -t / rho
c
            do 40 i = l, n
               b(i,j) = b(i,j) + t * b(i,l)
   40       continue
c
   50    continue
c
         do 80 j = 1, n
            t = 0.0d0
c
            do 60 i = l, n
               t = t + b(i,l) * a(i,j)
   60       continue
c
            t = -t / rho
c
            do 70 i = l, n
               a(i,j) = a(i,j) + t * b(i,l)
   70       continue
c
   80    continue
c
         b(l,l) = -s * r
c
         do 90 i = l1, n
            b(i,l) = 0.0d0
   90    continue
c
  100 continue
c     .......... reduce a to upper hessenberg form, while
c                keeping b triangular ..........
      if (n .eq. 2) go to 170
      nm2 = n - 2
c
      do 160 k = 1, nm2
         nk1 = nm1 - k
c     .......... for l=n-1 step -1 until k+1 do -- ..........
         do 150 lb = 1, nk1
            l = n - lb
            l1 = l + 1
c     .......... zero a(l+1,k) ..........
            s = dabs(a(l,k)) + dabs(a(l1,k))
            if (s .eq. 0.0d0) go to 150
            u1 = a(l,k) / s
            u2 = a(l1,k) / s
            r = dsign(dsqrt(u1*u1+u2*u2),u1)
            v1 =  -(u1 + r) / r
            v2 = -u2 / r
            u2 = v2 / v1
c
            do 110 j = k, n
               t = a(l,j) + u2 * a(l1,j)
               a(l,j) = a(l,j) + t * v1
               a(l1,j) = a(l1,j) + t * v2
  110       continue
c
            a(l1,k) = 0.0d0
c
            do 120 j = l, n
               t = b(l,j) + u2 * b(l1,j)
               b(l,j) = b(l,j) + t * v1
               b(l1,j) = b(l1,j) + t * v2
  120       continue
c     .......... zero b(l+1,l) ..........
            s = dabs(b(l1,l1)) + dabs(b(l1,l))
            if (s .eq. 0.0d0) go to 150
            u1 = b(l1,l1) / s
            u2 = b(l1,l) / s
            r = dsign(dsqrt(u1*u1+u2*u2),u1)
            v1 =  -(u1 + r) / r
            v2 = -u2 / r
            u2 = v2 / v1
c
            do 130 i = 1, l1
               t = b(i,l1) + u2 * b(i,l)
               b(i,l1) = b(i,l1) + t * v1
               b(i,l) = b(i,l) + t * v2
  130       continue
c
            b(l1,l) = 0.0d0
c
            do 140 i = 1, n
               t = a(i,l1) + u2 * a(i,l)
               a(i,l1) = a(i,l1) + t * v1
               a(i,l) = a(i,l) + t * v2
  140       continue
c
            if (.not. matz) go to 150
c
            do 145 i = 1, n
               t = z(i,l1) + u2 * z(i,l)
               z(i,l1) = z(i,l1) + t * v1
               z(i,l) = z(i,l) + t * v2
  145       continue
c
  150    continue
c
  160 continue
c
  170 return
      end
      subroutine qzit(nm,n,a,b,eps1,matz,z,ierr)

c*********************************************************************72
c
cc QZIT carries out iterations to solve a generalized eigenvalue problem.
c
c     this subroutine is the second step of the qz algorithm
c     for solving generalized matrix eigenvalue problems,
c     siam j. numer. anal. 10, 241-256(1973) by moler and stewart,
c     as modified in technical note nasa tn d-7305(1973) by ward.
c
c     this subroutine accepts a pair of real matrices, one of them
c     in upper hessenberg form and the other in upper triangular form.
c     it reduces the hessenberg matrix to quasi-triangular form using
c     orthogonal transformations while maintaining the triangular form
c     of the other matrix.  it is usually preceded by  qzhes  and
c     followed by  qzval  and, possibly,  qzvec.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrices.
c
c        a contains a real upper hessenberg matrix.
c
c        b contains a real upper triangular matrix.
c
c        eps1 is a tolerance used to determine negligible elements.
c          eps1 = 0.0 (or negative) may be input, in which case an
c          element will be neglected only if it is less than roundoff
c          error times the norm of its matrix.  if the input eps1 is
c          positive, then an element will be considered negligible
c          if it is less than eps1 times the norm of its matrix.  a
c          positive value of eps1 may result in faster execution,
c          but less accurate results.
c
c        matz should be set to .true. if the right hand transformations
c          are to be accumulated for later use in computing
c          eigenvectors, and to .false. otherwise.
c
c        z contains, if matz has been set to .true., the
c          transformation matrix produced in the reduction
c          by  qzhes, if performed, or else the identity matrix.
c          if matz has been set to .false., z is not referenced.
c
c     on output
c
c        a has been reduced to quasi-triangular form.  the elements
c          below the first subdiagonal are still zero and no two
c          consecutive subdiagonal elements are nonzero.
c
c        b is still in upper triangular form, although its elements
c          have been altered.  the location b(n,1) is used to store
c          eps1 times the norm of b for later use by  qzval  and  qzvec.
c
c        z contains the product of the right hand transformations
c          (for both steps) if matz has been set to .true..
c
c        ierr is set to
c          zero       for normal return,
c          j          if the limit of 30*n iterations is exhausted
c                     while the j-th eigenvalue is being sought.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,k,l,n,en,k1,k2,ld,ll,l1,na,nm,ish,itn,its,km1,lm1,
     x        enm2,ierr,lor1,enorn
      double precision a(nm,n),b(nm,n),z(nm,n)
      double precision r,s,t,a1,a2,a3,ep,sh,u1,u2,u3,v1,v2,v3,ani,a11,
     x       a12,a21,a22,a33,a34,a43,a44,bni,b11,b12,b22,b33,b34,
     x       b44,epsa,epsb,eps1,anorm,bnorm,epslon
      logical matz,notlas

      ierr = 0
c     .......... compute epsa,epsb ..........
      anorm = 0.0d0
      bnorm = 0.0d0
c
      do 30 i = 1, n
         ani = 0.0d0
         if (i .ne. 1) ani = dabs(a(i,i-1))
         bni = 0.0d0
c
         do 20 j = i, n
            ani = ani + dabs(a(i,j))
            bni = bni + dabs(b(i,j))
   20    continue
c
         if (ani .gt. anorm) anorm = ani
         if (bni .gt. bnorm) bnorm = bni
   30 continue
c
      if (anorm .eq. 0.0d0) anorm = 1.0d0
      if (bnorm .eq. 0.0d0) bnorm = 1.0d0
      ep = eps1
      if (ep .gt. 0.0d0) go to 50
c     .......... use roundoff level if eps1 is zero ..........
      ep = epslon(1.0d0)
   50 epsa = ep * anorm
      epsb = ep * bnorm
c     .......... reduce a to quasi-triangular form, while
c                keeping b triangular ..........
      lor1 = 1
      enorn = n
      en = n
      itn = 30*n
c     .......... begin qz step ..........
   60 if (en .le. 2) go to 1001
      if (.not. matz) enorn = en
      its = 0
      na = en - 1
      enm2 = na - 1
   70 ish = 2
c     .......... check for convergence or reducibility.
c                for l=en step -1 until 1 do -- ..........
      do 80 ll = 1, en
         lm1 = en - ll
         l = lm1 + 1
         if (l .eq. 1) go to 95
         if (dabs(a(l,lm1)) .le. epsa) go to 90
   80 continue
c
   90 a(l,lm1) = 0.0d0
      if (l .lt. na) go to 95
c     .......... 1-by-1 or 2-by-2 block isolated ..........
      en = lm1
      go to 60
c     .......... check for small top of b ..........
   95 ld = l
  100 l1 = l + 1
      b11 = b(l,l)
      if (dabs(b11) .gt. epsb) go to 120
      b(l,l) = 0.0d0
      s = dabs(a(l,l)) + dabs(a(l1,l))
      u1 = a(l,l) / s
      u2 = a(l1,l) / s
      r = dsign(dsqrt(u1*u1+u2*u2),u1)
      v1 = -(u1 + r) / r
      v2 = -u2 / r
      u2 = v2 / v1
c
      do 110 j = l, enorn
         t = a(l,j) + u2 * a(l1,j)
         a(l,j) = a(l,j) + t * v1
         a(l1,j) = a(l1,j) + t * v2
         t = b(l,j) + u2 * b(l1,j)
         b(l,j) = b(l,j) + t * v1
         b(l1,j) = b(l1,j) + t * v2
  110 continue
c
      if (l .ne. 1) a(l,lm1) = -a(l,lm1)
      lm1 = l
      l = l1
      go to 90
  120 a11 = a(l,l) / b11
      a21 = a(l1,l) / b11
      if (ish .eq. 1) go to 140
c     .......... iteration strategy ..........
      if (itn .eq. 0) go to 1000
      if (its .eq. 10) go to 155
c     .......... determine type of shift ..........
      b22 = b(l1,l1)
      if (dabs(b22) .lt. epsb) b22 = epsb
      b33 = b(na,na)
      if (dabs(b33) .lt. epsb) b33 = epsb
      b44 = b(en,en)
      if (dabs(b44) .lt. epsb) b44 = epsb
      a33 = a(na,na) / b33
      a34 = a(na,en) / b44
      a43 = a(en,na) / b33
      a44 = a(en,en) / b44
      b34 = b(na,en) / b44
      t = 0.5d0 * (a43 * b34 - a33 - a44)
      r = t * t + a34 * a43 - a33 * a44
      if (r .lt. 0.0d0) go to 150
c     .......... determine single shift zeroth column of a ..........
      ish = 1
      r = dsqrt(r)
      sh = -t + r
      s = -t - r
      if (dabs(s-a44) .lt. dabs(sh-a44)) sh = s
c     .......... look for two consecutive small
c                sub-diagonal elements of a.
c                for l=en-2 step -1 until ld do -- ..........
      do 130 ll = ld, enm2
         l = enm2 + ld - ll
         if (l .eq. ld) go to 140
         lm1 = l - 1
         l1 = l + 1
         t = a(l,l)
         if (dabs(b(l,l)) .gt. epsb) t = t - sh * b(l,l)
         if (dabs(a(l,lm1)) .le. dabs(t/a(l1,l)) * epsa) go to 100
  130 continue
c
  140 a1 = a11 - sh
      a2 = a21
      if (l .ne. ld) a(l,lm1) = -a(l,lm1)
      go to 160
c     .......... determine double shift zeroth column of a ..........
  150 a12 = a(l,l1) / b22
      a22 = a(l1,l1) / b22
      b12 = b(l,l1) / b22
      a1 = ((a33 - a11) * (a44 - a11) - a34 * a43 + a43 * b34 * a11)
     x     / a21 + a12 - a11 * b12
      a2 = (a22 - a11) - a21 * b12 - (a33 - a11) - (a44 - a11)
     x     + a43 * b34
      a3 = a(l1+1,l1) / b22
      go to 160
c     .......... ad hoc shift ..........
  155 a1 = 0.0d0
      a2 = 1.0d0
      a3 = 1.1605d0
  160 its = its + 1
      itn = itn - 1
      if (.not. matz) lor1 = ld
c     .......... main loop ..........
      do 260 k = l, na
         notlas = k .ne. na .and. ish .eq. 2
         k1 = k + 1
         k2 = k + 2
         km1 = max0(k-1,l)
         ll = min0(en,k1+ish)
         if (notlas) go to 190
c     .......... zero a(k+1,k-1) ..........
         if (k .eq. l) go to 170
         a1 = a(k,km1)
         a2 = a(k1,km1)
  170    s = dabs(a1) + dabs(a2)
         if (s .eq. 0.0d0) go to 70
         u1 = a1 / s
         u2 = a2 / s
         r = dsign(dsqrt(u1*u1+u2*u2),u1)
         v1 = -(u1 + r) / r
         v2 = -u2 / r
         u2 = v2 / v1
c
         do 180 j = km1, enorn
            t = a(k,j) + u2 * a(k1,j)
            a(k,j) = a(k,j) + t * v1
            a(k1,j) = a(k1,j) + t * v2
            t = b(k,j) + u2 * b(k1,j)
            b(k,j) = b(k,j) + t * v1
            b(k1,j) = b(k1,j) + t * v2
  180    continue
c
         if (k .ne. l) a(k1,km1) = 0.0d0
         go to 240
c     .......... zero a(k+1,k-1) and a(k+2,k-1) ..........
  190    if (k .eq. l) go to 200
         a1 = a(k,km1)
         a2 = a(k1,km1)
         a3 = a(k2,km1)
  200    s = dabs(a1) + dabs(a2) + dabs(a3)
         if (s .eq. 0.0d0) go to 260
         u1 = a1 / s
         u2 = a2 / s
         u3 = a3 / s
         r = dsign(dsqrt(u1*u1+u2*u2+u3*u3),u1)
         v1 = -(u1 + r) / r
         v2 = -u2 / r
         v3 = -u3 / r
         u2 = v2 / v1
         u3 = v3 / v1
c
         do 210 j = km1, enorn
            t = a(k,j) + u2 * a(k1,j) + u3 * a(k2,j)
            a(k,j) = a(k,j) + t * v1
            a(k1,j) = a(k1,j) + t * v2
            a(k2,j) = a(k2,j) + t * v3
            t = b(k,j) + u2 * b(k1,j) + u3 * b(k2,j)
            b(k,j) = b(k,j) + t * v1
            b(k1,j) = b(k1,j) + t * v2
            b(k2,j) = b(k2,j) + t * v3
  210    continue
c
         if (k .eq. l) go to 220
         a(k1,km1) = 0.0d0
         a(k2,km1) = 0.0d0
c     .......... zero b(k+2,k+1) and b(k+2,k) ..........
  220    s = dabs(b(k2,k2)) + dabs(b(k2,k1)) + dabs(b(k2,k))
         if (s .eq. 0.0d0) go to 240
         u1 = b(k2,k2) / s
         u2 = b(k2,k1) / s
         u3 = b(k2,k) / s
         r = dsign(dsqrt(u1*u1+u2*u2+u3*u3),u1)
         v1 = -(u1 + r) / r
         v2 = -u2 / r
         v3 = -u3 / r
         u2 = v2 / v1
         u3 = v3 / v1
c
         do 230 i = lor1, ll
            t = a(i,k2) + u2 * a(i,k1) + u3 * a(i,k)
            a(i,k2) = a(i,k2) + t * v1
            a(i,k1) = a(i,k1) + t * v2
            a(i,k) = a(i,k) + t * v3
            t = b(i,k2) + u2 * b(i,k1) + u3 * b(i,k)
            b(i,k2) = b(i,k2) + t * v1
            b(i,k1) = b(i,k1) + t * v2
            b(i,k) = b(i,k) + t * v3
  230    continue
c
         b(k2,k) = 0.0d0
         b(k2,k1) = 0.0d0
         if (.not. matz) go to 240
c
         do 235 i = 1, n
            t = z(i,k2) + u2 * z(i,k1) + u3 * z(i,k)
            z(i,k2) = z(i,k2) + t * v1
            z(i,k1) = z(i,k1) + t * v2
            z(i,k) = z(i,k) + t * v3
  235    continue
c     .......... zero b(k+1,k) ..........
  240    s = dabs(b(k1,k1)) + dabs(b(k1,k))
         if (s .eq. 0.0d0) go to 260
         u1 = b(k1,k1) / s
         u2 = b(k1,k) / s
         r = dsign(dsqrt(u1*u1+u2*u2),u1)
         v1 = -(u1 + r) / r
         v2 = -u2 / r
         u2 = v2 / v1
c
         do 250 i = lor1, ll
            t = a(i,k1) + u2 * a(i,k)
            a(i,k1) = a(i,k1) + t * v1
            a(i,k) = a(i,k) + t * v2
            t = b(i,k1) + u2 * b(i,k)
            b(i,k1) = b(i,k1) + t * v1
            b(i,k) = b(i,k) + t * v2
  250    continue
c
         b(k1,k) = 0.0d0
         if (.not. matz) go to 260
c
         do 255 i = 1, n
            t = z(i,k1) + u2 * z(i,k)
            z(i,k1) = z(i,k1) + t * v1
            z(i,k) = z(i,k) + t * v2
  255    continue
c
  260 continue
c     .......... end qz step ..........
      go to 70
c     .......... set error -- all eigenvalues have not
c                converged after 30*n iterations ..........
 1000 ierr = en
c     .......... save epsb for use by qzval and qzvec ..........
 1001 if (n .gt. 1) b(n,1) = epsb
      return
      end
      subroutine qzval(nm,n,a,b,alfr,alfi,beta,matz,z)

c*********************************************************************72
c
cc QZVAL computes eigenvalues for a generalized eigenvalue problem.
c
c     this subroutine is the third step of the qz algorithm
c     for solving generalized matrix eigenvalue problems,
c     siam j. numer. anal. 10, 241-256(1973) by moler and stewart.
c
c     this subroutine accepts a pair of real matrices, one of them
c     in quasi-triangular form and the other in upper triangular form.
c     it reduces the quasi-triangular matrix further, so that any
c     remaining 2-by-2 blocks correspond to pairs of complex
c     eigenvalues, and returns quantities whose ratios give the
c     generalized eigenvalues.  it is usually preceded by  qzhes
c     and  qzit  and may be followed by  qzvec.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrices.
c
c        a contains a real upper quasi-triangular matrix.
c
c        b contains a real upper triangular matrix.  in addition,
c          location b(n,1) contains the tolerance quantity (epsb)
c          computed and saved in  qzit.
c
c        matz should be set to .true. if the right hand transformations
c          are to be accumulated for later use in computing
c          eigenvectors, and to .false. otherwise.
c
c        z contains, if matz has been set to .true., the
c          transformation matrix produced in the reductions by qzhes
c          and qzit, if performed, or else the identity matrix.
c          if matz has been set to .false., z is not referenced.
c
c     on output
c
c        a has been reduced further to a quasi-triangular matrix
c          in which all nonzero subdiagonal elements correspond to
c          pairs of complex eigenvalues.
c
c        b is still in upper triangular form, although its elements
c          have been altered.  b(n,1) is unaltered.
c
c        alfr and alfi contain the real and imaginary parts of the
c          diagonal elements of the triangular matrix that would be
c          obtained if a were reduced completely to triangular form
c          by unitary transformations.  non-zero values of alfi occur
c          in pairs, the first member positive and the second negative.
c
c        beta contains the diagonal elements of the corresponding b,
c          normalized to be real and non-negative.  the generalized
c          eigenvalues are then the ratios ((alfr+i*alfi)/beta).
c
c        z contains the product of the right hand transformations
c          (for all three steps) if matz has been set to .true.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,n,en,na,nm,nn,isw
      double precision a(nm,n),b(nm,n),alfr(n),alfi(n),beta(n),z(nm,n)
      double precision c,d,e,r,s,t,an,a1,a2,bn,cq,cz,di,dr,ei,ti,tr,u1,
     x       u2,v1,v2,a1i,a11,a12,a2i,a21,a22,b11,b12,b22,sqi,sqr,
     x       ssi,ssr,szi,szr,a11i,a11r,a12i,a12r,a22i,a22r,epsb
      logical matz

      epsb = b(n,1)
      isw = 1
c     .......... find eigenvalues of quasi-triangular matrices.
c                for en=n step -1 until 1 do -- ..........
      do 510 nn = 1, n
         en = n + 1 - nn
         na = en - 1
         if (isw .eq. 2) go to 505
         if (en .eq. 1) go to 410
         if (a(en,na) .ne. 0.0d0) go to 420
c     .......... 1-by-1 block, one real root ..........
  410    alfr(en) = a(en,en)
         if (b(en,en) .lt. 0.0d0) alfr(en) = -alfr(en)
         beta(en) = dabs(b(en,en))
         alfi(en) = 0.0d0
         go to 510
c     .......... 2-by-2 block ..........
  420    if (dabs(b(na,na)) .le. epsb) go to 455
         if (dabs(b(en,en)) .gt. epsb) go to 430
         a1 = a(en,en)
         a2 = a(en,na)
         bn = 0.0d0
         go to 435
  430    an = dabs(a(na,na)) + dabs(a(na,en)) + dabs(a(en,na))
     x      + dabs(a(en,en))
         bn = dabs(b(na,na)) + dabs(b(na,en)) + dabs(b(en,en))
         a11 = a(na,na) / an
         a12 = a(na,en) / an
         a21 = a(en,na) / an
         a22 = a(en,en) / an
         b11 = b(na,na) / bn
         b12 = b(na,en) / bn
         b22 = b(en,en) / bn
         e = a11 / b11
         ei = a22 / b22
         s = a21 / (b11 * b22)
         t = (a22 - e * b22) / b22
         if (dabs(e) .le. dabs(ei)) go to 431
         e = ei
         t = (a11 - e * b11) / b11
  431    c = 0.5d0 * (t - s * b12)
         d = c * c + s * (a12 - e * b12)
         if (d .lt. 0.0d0) go to 480
c     .......... two real roots.
c                zero both a(en,na) and b(en,na) ..........
         e = e + (c + dsign(dsqrt(d),c))
         a11 = a11 - e * b11
         a12 = a12 - e * b12
         a22 = a22 - e * b22
         if (dabs(a11) + dabs(a12) .lt.
     x       dabs(a21) + dabs(a22)) go to 432
         a1 = a12
         a2 = a11
         go to 435
  432    a1 = a22
         a2 = a21
c     .......... choose and apply real z ..........
  435    s = dabs(a1) + dabs(a2)
         u1 = a1 / s
         u2 = a2 / s
         r = dsign(dsqrt(u1*u1+u2*u2),u1)
         v1 = -(u1 + r) / r
         v2 = -u2 / r
         u2 = v2 / v1
c
         do 440 i = 1, en
            t = a(i,en) + u2 * a(i,na)
            a(i,en) = a(i,en) + t * v1
            a(i,na) = a(i,na) + t * v2
            t = b(i,en) + u2 * b(i,na)
            b(i,en) = b(i,en) + t * v1
            b(i,na) = b(i,na) + t * v2
  440    continue
c
         if (.not. matz) go to 450
c
         do 445 i = 1, n
            t = z(i,en) + u2 * z(i,na)
            z(i,en) = z(i,en) + t * v1
            z(i,na) = z(i,na) + t * v2
  445    continue
c
  450    if (bn .eq. 0.0d0) go to 475
         if (an .lt. dabs(e) * bn) go to 455
         a1 = b(na,na)
         a2 = b(en,na)
         go to 460
  455    a1 = a(na,na)
         a2 = a(en,na)
c     .......... choose and apply real q ..........
  460    s = dabs(a1) + dabs(a2)
         if (s .eq. 0.0d0) go to 475
         u1 = a1 / s
         u2 = a2 / s
         r = dsign(dsqrt(u1*u1+u2*u2),u1)
         v1 = -(u1 + r) / r
         v2 = -u2 / r
         u2 = v2 / v1
c
         do 470 j = na, n
            t = a(na,j) + u2 * a(en,j)
            a(na,j) = a(na,j) + t * v1
            a(en,j) = a(en,j) + t * v2
            t = b(na,j) + u2 * b(en,j)
            b(na,j) = b(na,j) + t * v1
            b(en,j) = b(en,j) + t * v2
  470    continue
c
  475    a(en,na) = 0.0d0
         b(en,na) = 0.0d0
         alfr(na) = a(na,na)
         alfr(en) = a(en,en)
         if (b(na,na) .lt. 0.0d0) alfr(na) = -alfr(na)
         if (b(en,en) .lt. 0.0d0) alfr(en) = -alfr(en)
         beta(na) = dabs(b(na,na))
         beta(en) = dabs(b(en,en))
         alfi(en) = 0.0d0
         alfi(na) = 0.0d0
         go to 505
c     .......... two complex roots ..........
  480    e = e + c
         ei = dsqrt(-d)
         a11r = a11 - e * b11
         a11i = ei * b11
         a12r = a12 - e * b12
         a12i = ei * b12
         a22r = a22 - e * b22
         a22i = ei * b22
         if (dabs(a11r) + dabs(a11i) + dabs(a12r) + dabs(a12i) .lt.
     x       dabs(a21) + dabs(a22r) + dabs(a22i)) go to 482
         a1 = a12r
         a1i = a12i
         a2 = -a11r
         a2i = -a11i
         go to 485
  482    a1 = a22r
         a1i = a22i
         a2 = -a21
         a2i = 0.0d0
c     .......... choose complex z ..........
  485    cz = dsqrt(a1*a1+a1i*a1i)
         if (cz .eq. 0.0d0) go to 487
         szr = (a1 * a2 + a1i * a2i) / cz
         szi = (a1 * a2i - a1i * a2) / cz
         r = dsqrt(cz*cz+szr*szr+szi*szi)
         cz = cz / r
         szr = szr / r
         szi = szi / r
         go to 490
  487    szr = 1.0d0
         szi = 0.0d0
  490    if (an .lt. (dabs(e) + ei) * bn) go to 492
         a1 = cz * b11 + szr * b12
         a1i = szi * b12
         a2 = szr * b22
         a2i = szi * b22
         go to 495
  492    a1 = cz * a11 + szr * a12
         a1i = szi * a12
         a2 = cz * a21 + szr * a22
         a2i = szi * a22
c     .......... choose complex q ..........
  495    cq = dsqrt(a1*a1+a1i*a1i)
         if (cq .eq. 0.0d0) go to 497
         sqr = (a1 * a2 + a1i * a2i) / cq
         sqi = (a1 * a2i - a1i * a2) / cq
         r = dsqrt(cq*cq+sqr*sqr+sqi*sqi)
         cq = cq / r
         sqr = sqr / r
         sqi = sqi / r
         go to 500
  497    sqr = 1.0d0
         sqi = 0.0d0
c     .......... compute diagonal elements that would result
c                if transformations were applied ..........
  500    ssr = sqr * szr + sqi * szi
         ssi = sqr * szi - sqi * szr
         i = 1
         tr = cq * cz * a11 + cq * szr * a12 + sqr * cz * a21
     x      + ssr * a22
         ti = cq * szi * a12 - sqi * cz * a21 + ssi * a22
         dr = cq * cz * b11 + cq * szr * b12 + ssr * b22
         di = cq * szi * b12 + ssi * b22
         go to 503
  502    i = 2
         tr = ssr * a11 - sqr * cz * a12 - cq * szr * a21
     x      + cq * cz * a22
         ti = -ssi * a11 - sqi * cz * a12 + cq * szi * a21
         dr = ssr * b11 - sqr * cz * b12 + cq * cz * b22
         di = -ssi * b11 - sqi * cz * b12
  503    t = ti * dr - tr * di
         j = na
         if (t .lt. 0.0d0) j = en
         r = dsqrt(dr*dr+di*di)
         beta(j) = bn * r
         alfr(j) = an * (tr * dr + ti * di) / r
         alfi(j) = an * t / r
         if (i .eq. 1) go to 502
  505    isw = 3 - isw
  510 continue
      b(n,1) = epsb
c
      return
      end
      subroutine qzvec(nm,n,a,b,alfr,alfi,beta,z)

c*********************************************************************72
c
cc QZVEC computes eigenvectors for a generalized eigenvalue problem.
c
c     this subroutine is the optional fourth step of the qz algorithm
c     for solving generalized matrix eigenvalue problems,
c     siam j. numer. anal. 10, 241-256(1973) by moler and stewart.
c
c     this subroutine accepts a pair of real matrices, one of them in
c     quasi-triangular form (in which each 2-by-2 block corresponds to
c     a pair of complex eigenvalues) and the other in upper triangular
c     form.  it computes the eigenvectors of the triangular problem and
c     transforms the results back to the original coordinate system.
c     it is usually preceded by  qzhes,  qzit, and  qzval.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrices.
c
c        a contains a real upper quasi-triangular matrix.
c
c        b contains a real upper triangular matrix.  in addition,
c          location b(n,1) contains the tolerance quantity (epsb)
c          computed and saved in  qzit.
c
c        alfr, alfi, and beta  are vectors with components whose
c          ratios ((alfr+i*alfi)/beta) are the generalized
c          eigenvalues.  they are usually obtained from  qzval.
c
c        z contains the transformation matrix produced in the
c          reductions by  qzhes,  qzit, and  qzval, if performed.
c          if the eigenvectors of the triangular problem are
c          desired, z must contain the identity matrix.
c
c     on output
c
c        a is unaltered.  its subdiagonal elements provide information
c           about the storage of the complex eigenvectors.
c
c        b has been destroyed.
c
c        alfr, alfi, and beta are unaltered.
c
c        z contains the real and imaginary parts of the eigenvectors.
c          if alfi(i) .eq. 0.0, the i-th eigenvalue is real and
c            the i-th column of z contains its eigenvector.
c          if alfi(i) .ne. 0.0, the i-th eigenvalue is complex.
c            if alfi(i) .gt. 0.0, the eigenvalue is the first of
c              a complex pair and the i-th and (i+1)-th columns
c              of z contain its eigenvector.
c            if alfi(i) .lt. 0.0, the eigenvalue is the second of
c              a complex pair and the (i-1)-th and i-th columns
c              of z contain the conjugate of its eigenvector.
c          each eigenvector is normalized so that the modulus
c          of its largest component is 1.0 .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,k,m,n,en,ii,jj,na,nm,nn,isw,enm2
      double precision a(nm,n),b(nm,n),alfr(n),alfi(n),beta(n),z(nm,n)
      double precision d,q,r,s,t,w,x,y,di,dr,ra,rr,sa,ti,tr,t1,t2,w1,x1,
     x       zz,z1,alfm,almi,almr,betm,epsb

      epsb = b(n,1)
      isw = 1
c     .......... for en=n step -1 until 1 do -- ..........
      do 800 nn = 1, n
         en = n + 1 - nn
         na = en - 1
         if (isw .eq. 2) go to 795
         if (alfi(en) .ne. 0.0d0) go to 710
c     .......... real vector ..........
         m = en
         b(en,en) = 1.0d0
         if (na .eq. 0) go to 800
         alfm = alfr(m)
         betm = beta(m)
c     .......... for i=en-1 step -1 until 1 do -- ..........
         do 700 ii = 1, na
            i = en - ii
            w = betm * a(i,i) - alfm * b(i,i)
            r = 0.0d0
c
            do 610 j = m, en
  610       r = r + (betm * a(i,j) - alfm * b(i,j)) * b(j,en)
c
            if (i .eq. 1 .or. isw .eq. 2) go to 630
            if (betm * a(i,i-1) .eq. 0.0d0) go to 630
            zz = w
            s = r
            go to 690
  630       m = i
            if (isw .eq. 2) go to 640
c     .......... real 1-by-1 block ..........
            t = w
            if (w .eq. 0.0d0) t = epsb
            b(i,en) = -r / t
            go to 700
c     .......... real 2-by-2 block ..........
  640       x = betm * a(i,i+1) - alfm * b(i,i+1)
            y = betm * a(i+1,i)
            q = w * zz - x * y
            t = (x * s - zz * r) / q
            b(i,en) = t
            if (dabs(x) .le. dabs(zz)) go to 650
            b(i+1,en) = (-r - w * t) / x
            go to 690
  650       b(i+1,en) = (-s - y * t) / zz
  690       isw = 3 - isw
  700    continue
c     .......... end real vector ..........
         go to 800
c     .......... complex vector ..........
  710    m = na
         almr = alfr(m)
         almi = alfi(m)
         betm = beta(m)
c     .......... last vector component chosen imaginary so that
c                eigenvector matrix is triangular ..........
         y = betm * a(en,na)
         b(na,na) = -almi * b(en,en) / y
         b(na,en) = (almr * b(en,en) - betm * a(en,en)) / y
         b(en,na) = 0.0d0
         b(en,en) = 1.0d0
         enm2 = na - 1
         if (enm2 .eq. 0) go to 795
c     .......... for i=en-2 step -1 until 1 do -- ..........
         do 790 ii = 1, enm2
            i = na - ii
            w = betm * a(i,i) - almr * b(i,i)
            w1 = -almi * b(i,i)
            ra = 0.0d0
            sa = 0.0d0
c
            do 760 j = m, en
               x = betm * a(i,j) - almr * b(i,j)
               x1 = -almi * b(i,j)
               ra = ra + x * b(j,na) - x1 * b(j,en)
               sa = sa + x * b(j,en) + x1 * b(j,na)
  760       continue
c
            if (i .eq. 1 .or. isw .eq. 2) go to 770
            if (betm * a(i,i-1) .eq. 0.0d0) go to 770
            zz = w
            z1 = w1
            r = ra
            s = sa
            isw = 2
            go to 790
  770       m = i
            if (isw .eq. 2) go to 780
c     .......... complex 1-by-1 block ..........
            tr = -ra
            ti = -sa
  773       dr = w
            di = w1
c     .......... complex divide (t1,t2) = (tr,ti) / (dr,di) ..........
  775       if (dabs(di) .gt. dabs(dr)) go to 777
            rr = di / dr
            d = dr + di * rr
            t1 = (tr + ti * rr) / d
            t2 = (ti - tr * rr) / d
            go to (787,782), isw
  777       rr = dr / di
            d = dr * rr + di
            t1 = (tr * rr + ti) / d
            t2 = (ti * rr - tr) / d
            go to (787,782), isw
c     .......... complex 2-by-2 block ..........
  780       x = betm * a(i,i+1) - almr * b(i,i+1)
            x1 = -almi * b(i,i+1)
            y = betm * a(i+1,i)
            tr = y * ra - w * r + w1 * s
            ti = y * sa - w * s - w1 * r
            dr = w * zz - w1 * z1 - x * y
            di = w * z1 + w1 * zz - x1 * y
            if (dr .eq. 0.0d0 .and. di .eq. 0.0d0) dr = epsb
            go to 775
  782       b(i+1,na) = t1
            b(i+1,en) = t2
            isw = 1
            if (dabs(y) .gt. dabs(w) + dabs(w1)) go to 785
            tr = -ra - x * b(i+1,na) + x1 * b(i+1,en)
            ti = -sa - x * b(i+1,en) - x1 * b(i+1,na)
            go to 773
  785       t1 = (-r - zz * b(i+1,na) + z1 * b(i+1,en)) / y
            t2 = (-s - zz * b(i+1,en) - z1 * b(i+1,na)) / y
  787       b(i,na) = t1
            b(i,en) = t2
  790    continue
c     .......... end complex vector ..........
  795    isw = 3 - isw
  800 continue
c     .......... end back substitution.
c                transform to original coordinate system.
c                for j=n step -1 until 1 do -- ..........
      do 880 jj = 1, n
         j = n + 1 - jj
c
         do 880 i = 1, n
            zz = 0.0d0
c
            do 860 k = 1, j
  860       zz = zz + z(i,k) * b(k,j)
c
            z(i,j) = zz
  880 continue
c     .......... normalize so that modulus of largest
c                component of each vector is 1.
c                (isw is 1 initially from before) ..........
      do 950 j = 1, n
         d = 0.0d0
         if (isw .eq. 2) go to 920
         if (alfi(j) .ne. 0.0d0) go to 945
c
         do 890 i = 1, n
            if (dabs(z(i,j)) .gt. d) d = dabs(z(i,j))
  890    continue
c
         do 900 i = 1, n
  900    z(i,j) = z(i,j) / d
c
         go to 950
c
  920    do 930 i = 1, n
            r = dabs(z(i,j-1)) + dabs(z(i,j))
            if (r .ne. 0.0d0) r = r * dsqrt((z(i,j-1)/r)**2
     x                                     +(z(i,j)/r)**2)
            if (r .gt. d) d = r
  930    continue
c
         do 940 i = 1, n
            z(i,j-1) = z(i,j-1) / d
            z(i,j) = z(i,j) / d
  940    continue
c
  945    isw = 3 - isw
  950 continue
c
      return
      end
      subroutine r8mat_uniform_01 ( m, n, seed, r )

c*********************************************************************72
c
cc R8MAT_UNIFORM_01 returns a unit pseudorandom R8MAT.
c
c  Discussion:
c
c    An R8MAT is an array of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 August 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Springer Verlag, pages 201-202, 1983.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, pages 362-376, 1986.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, pages 136-143, 1969.
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns in the array.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R(M,N), the array of pseudorandom values.
c
      implicit none

      integer m
      integer n

      integer i
      integer j
      integer k
      integer seed
      double precision r(m,n)

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8MAT_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      do j = 1, n

        do i = 1, m

          k = seed / 127773

          seed = 16807 * ( seed - k * 127773 ) - k * 2836

          if ( seed .lt. 0 ) then
            seed = seed + 2147483647
          end if

          r(i,j) = dble ( seed ) * 4.656612875D-10

        end do
      end do

      return
      end
      subroutine ratqr(n,eps1,d,e,e2,m,w,ind,bd,type,idef,ierr)

c*********************************************************************72
c
cc RATQR computes selected eigenvalues of a real symmetric tridiagonal matrix.
c
c     this subroutine is a translation of the algol procedure ratqr,
c     num. math. 11, 264-272(1968) by reinsch and bauer.
c     handbook for auto. comp., vol.ii-linear algebra, 257-265(1971).
c
c     this subroutine finds the algebraically smallest or largest
c     eigenvalues of a symmetric tridiagonal matrix by the
c     rational qr method with newton corrections.
c
c     on input
c
c        n is the order of the matrix.
c
c        eps1 is a theoretical absolute error tolerance for the
c          computed eigenvalues.  if the input eps1 is non-positive,
c          or indeed smaller than its default value, it is reset
c          at each iteration to the respective default value,
c          namely, the product of the relative machine precision
c          and the magnitude of the current eigenvalue iterate.
c          the theoretical absolute error in the k-th eigenvalue
c          is usually not greater than k times eps1.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c        e2 contains the squares of the corresponding elements of e.
c          e2(1) is arbitrary.
c
c        m is the number of eigenvalues to be found.
c
c        idef should be set to 1 if the input matrix is known to be
c          positive definite, to -1 if the input matrix is known to
c          be negative definite, and to 0 otherwise.
c
c        type should be set to .true. if the smallest eigenvalues
c          are to be found, and to .false. if the largest eigenvalues
c          are to be found.
c
c     on output
c
c        eps1 is unaltered unless it has been reset to its
c          (last) default value.
c
c        d and e are unaltered (unless w overwrites d).
c
c        elements of e2, corresponding to elements of e regarded
c          as negligible, have been replaced by zero causing the
c          matrix to split into a direct sum of submatrices.
c          e2(1) is set to 0.0d0 if the smallest eigenvalues have been
c          found, and to 2.0d0 if the largest eigenvalues have been
c          found.  e2 is otherwise unaltered (unless overwritten by bd).
c
c        w contains the m algebraically smallest eigenvalues in
c          ascending order, or the m largest eigenvalues in
c          descending order.  if an error exit is made because of
c          an incorrect specification of idef, no eigenvalues
c          are found.  if the newton iterates for a particular
c          eigenvalue are not monotone, the best estimate obtained
c          is returned and ierr is set.  w may coincide with d.
c
c        ind contains in its first m positions the submatrix indices
c          associated with the corresponding eigenvalues in w --
c          1 for eigenvalues belonging to the first submatrix from
c          the top, 2 for those belonging to the second submatrix, etc..
c
c        bd contains refined bounds for the theoretical errors of the
c          corresponding eigenvalues in w.  these bounds are usually
c          within the tolerance specified by eps1.  bd may coincide
c          with e2.
c
c        ierr is set to
c          zero       for normal return,
c          6*n+1      if  idef  is set to 1 and  type  to .true.
c                     when the matrix is not positive definite, or
c                     if  idef  is set to -1 and  type  to .false.
c                     when the matrix is not negative definite,
c          5*n+k      if successive iterates to the k-th eigenvalue
c                     are not monotone increasing, where k refers
c                     to the last such occurrence.
c
c     note that subroutine tridib is generally faster and more
c     accurate than ratqr if the eigenvalues are clustered.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,k,m,n,ii,jj,k1,idef,ierr,jdef
      double precision d(n),e(n),e2(n),w(n),bd(n)
      double precision f,p,q,r,s,ep,qp,err,tot,eps1,delta,epslon
      integer ind(n)
      logical type

      ierr = 0
      jdef = idef
c     .......... copy d array into w ..........
      do 20 i = 1, n
   20 w(i) = d(i)
c
      if (type) go to 40
      j = 1
      go to 400
   40 err = 0.0d0
      s = 0.0d0
c     .......... look for small sub-diagonal entries and define
c                initial shift from lower gerschgorin bound.
c                copy e2 array into bd ..........
      tot = w(1)
      q = 0.0d0
      j = 0
c
      do 100 i = 1, n
         p = q
         if (i .eq. 1) go to 60
         if (p .gt. epslon(dabs(d(i)) + dabs(d(i-1)))) go to 80
   60    e2(i) = 0.0d0
   80    bd(i) = e2(i)
c     .......... count also if element of e2 has underflowed ..........
         if (e2(i) .eq. 0.0d0) j = j + 1
         ind(i) = j
         q = 0.0d0
         if (i .ne. n) q = dabs(e(i+1))
         tot = dmin1(w(i)-p-q,tot)
  100 continue
c
      if (jdef .eq. 1 .and. tot .lt. 0.0d0) go to 140
c
      do 110 i = 1, n
  110 w(i) = w(i) - tot
c
      go to 160
  140 tot = 0.0d0
c
  160 do 360 k = 1, m
c     .......... next qr transformation ..........
  180    tot = tot + s
         delta = w(n) - s
         i = n
         f = dabs(epslon(tot))
         if (eps1 .lt. f) eps1 = f
         if (delta .gt. eps1) go to 190
         if (delta .lt. (-eps1)) go to 1000
         go to 300
c     .......... replace small sub-diagonal squares by zero
c                to reduce the incidence of underflows ..........
  190    if (k .eq. n) go to 210
         k1 = k + 1
         do 200 j = k1, n
            if (bd(j) .le. (epslon(w(j)+w(j-1))) ** 2) bd(j) = 0.0d0
  200    continue
c
  210    f = bd(n) / delta
         qp = delta + f
         p = 1.0d0
         if (k .eq. n) go to 260
         k1 = n - k
c     .......... for i=n-1 step -1 until k do -- ..........
         do 240 ii = 1, k1
            i = n - ii
            q = w(i) - s - f
            r = q / qp
            p = p * r + 1.0d0
            ep = f * r
            w(i+1) = qp + ep
            delta = q - ep
            if (delta .gt. eps1) go to 220
            if (delta .lt. (-eps1)) go to 1000
            go to 300
  220       f = bd(i) / q
            qp = delta + f
            bd(i+1) = qp * ep
  240    continue
c
  260    w(k) = qp
         s = qp / p
         if (tot + s .gt. tot) go to 180
c     .......... set error -- irregular end of iteration.
c                deflate minimum diagonal element ..........
         ierr = 5 * n + k
         s = 0.0d0
         delta = qp
c
         do 280 j = k, n
            if (w(j) .gt. delta) go to 280
            i = j
            delta = w(j)
  280    continue
c     .......... convergence ..........
  300    if (i .lt. n) bd(i+1) = bd(i) * f / qp
         ii = ind(i)
         if (i .eq. k) go to 340
         k1 = i - k
c     .......... for j=i-1 step -1 until k do -- ..........
         do 320 jj = 1, k1
            j = i - jj
            w(j+1) = w(j) - s
            bd(j+1) = bd(j)
            ind(j+1) = ind(j)
  320    continue
c
  340    w(k) = tot
         err = err + dabs(delta)
         bd(k) = err
         ind(k) = ii
  360 continue
c
      if (type) go to 1001
      f = bd(1)
      e2(1) = 2.0d0
      bd(1) = f
      j = 2
c     .......... negate elements of w for largest values ..........
  400 do 500 i = 1, n
  500 w(i) = -w(i)
c
      jdef = -jdef
      go to (40,1001), j
c     .......... set error -- idef specified incorrectly ..........
 1000 ierr = 6 * n + 1
 1001 return
      end
      subroutine rebak(nm,n,b,dl,m,z)

c*********************************************************************72
c
cc REBAK determines eigenvectors by undoing the REDUC transformation.
c
c     this subroutine is a translation of the algol procedure rebaka,
c     num. math. 11, 99-110(1968) by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 303-314(1971).
c
c     this subroutine forms the eigenvectors of a generalized
c     symmetric eigensystem by back transforming those of the
c     derived symmetric matrix determined by  reduc.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix system.
c
c        b contains information about the similarity transformation
c          (cholesky decomposition) used in the reduction by  reduc
c          in its strict lower triangle.
c
c        dl contains further information about the transformation.
c
c        m is the number of eigenvectors to be back transformed.
c
c        z contains the eigenvectors to be back transformed
c          in its first m columns.
c
c     on output
c
c        z contains the transformed eigenvectors
c          in its first m columns.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,k,m,n,i1,ii,nm
      double precision b(nm,n),dl(n),z(nm,m)
      double precision x
c
      if (m .eq. 0) go to 200
c
      do 100 j = 1, m
c     .......... for i=n step -1 until 1 do -- ..........
         do 100 ii = 1, n
            i = n + 1 - ii
            i1 = i + 1
            x = z(i,j)
            if (i .eq. n) go to 80
c
            do 60 k = i1, n
   60       x = x - b(k,i) * z(k,j)
c
   80       z(i,j) = x / dl(i)
  100 continue
c
  200 return
      end
      subroutine rebakb(nm,n,b,dl,m,z)

c*********************************************************************72
c
cc REBAKB determines eigenvectors by undoing the REDUC2 transformation.
c
c     this subroutine is a translation of the algol procedure rebakb,
c     num. math. 11, 99-110(1968) by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 303-314(1971).
c
c     this subroutine forms the eigenvectors of a generalized
c     symmetric eigensystem by back transforming those of the
c     derived symmetric matrix determined by REDUC2.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix system.
c
c        b contains information about the similarity transformation
c          (cholesky decomposition) used in the reduction by REDUC2
c          in its strict lower triangle.
c
c        dl contains further information about the transformation.
c
c        m is the number of eigenvectors to be back transformed.
c
c        z contains the eigenvectors to be back transformed
c          in its first m columns.
c
c     on output
c
c        z contains the transformed eigenvectors
c          in its first m columns.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,k,m,n,i1,ii,nm
      double precision b(nm,n),dl(n),z(nm,m)
      double precision x
c
      if (m .eq. 0) go to 200
c
      do 100 j = 1, m
c     .......... for i=n step -1 until 1 do -- ..........
         do 100 ii = 1, n
            i1 = n - ii
            i = i1 + 1
            x = dl(i) * z(i,j)
            if (i .eq. 1) go to 80
c
            do 60 k = 1, i1
   60       x = x + b(i,k) * z(k,j)
c
   80       z(i,j) = x
  100 continue
c
  200 return
      end
      subroutine reduc(nm,n,a,b,dl,ierr)

c*********************************************************************72
c
cc REDUC reduces the eigenvalue problem A*x=lambda*B*x to A*x=lambda*x.
c
c     this subroutine is a translation of the algol procedure reduc1,
c     num. math. 11, 99-110(1968) by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 303-314(1971).
c
c     this subroutine reduces the generalized symmetric eigenproblem
c     ax=(lambda)bx, where b is positive definite, to the standard
c     symmetric eigenproblem using the cholesky factorization of b.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrices a and b.  if the cholesky
c          factor l of b is already available, n should be prefixed
c          with a minus sign.
c
c        a and b contain the real symmetric input matrices.  only the
c          full upper triangles of the matrices need be supplied.  if
c          n is negative, the strict lower triangle of b contains,
c          instead, the strict lower triangle of its cholesky factor l.
c
c        dl contains, if n is negative, the diagonal elements of l.
c
c     on output
c
c        a contains in its full lower triangle the full lower triangle
c          of the symmetric matrix derived from the reduction to the
c          standard form.  the strict upper triangle of a is unaltered.
c
c        b contains in its strict lower triangle the strict lower
c          triangle of its cholesky factor l.  the full upper
c          triangle of b is unaltered.
c
c        dl contains the diagonal elements of l.
c
c        ierr is set to
c          zero       for normal return,
c          7*n+1      if b is not positive definite.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,k,n,i1,j1,nm,nn,ierr
      double precision a(nm,n),b(nm,n),dl(n)
      double precision x,y
c
      ierr = 0
      nn = iabs(n)
      if (n .lt. 0) go to 100
c     .......... form l in the arrays b and dl ..........
      do 80 i = 1, n
         i1 = i - 1
c
         do 80 j = i, n
            x = b(i,j)
            if (i .eq. 1) go to 40
c
            do 20 k = 1, i1
   20       x = x - b(i,k) * b(j,k)
c
   40       if (j .ne. i) go to 60
            if (x .le. 0.0d0) go to 1000
            y = dsqrt(x)
            dl(i) = y
            go to 80
   60       b(j,i) = x / y
   80 continue
c     .......... form the transpose of the upper triangle of inv(l)*a
c                in the lower triangle of the array a ..........
  100 do 200 i = 1, nn
         i1 = i - 1
         y = dl(i)
c
         do 200 j = i, nn
            x = a(i,j)
            if (i .eq. 1) go to 180
c
            do 160 k = 1, i1
  160       x = x - b(i,k) * a(j,k)
c
  180       a(j,i) = x / y
  200 continue
c     .......... pre-multiply by inv(l) and overwrite ..........
      do 300 j = 1, nn
         j1 = j - 1
c
         do 300 i = j, nn
            x = a(i,j)
            if (i .eq. j) go to 240
            i1 = i - 1
c
            do 220 k = j, i1
  220       x = x - a(k,j) * b(i,k)
c
  240       if (j .eq. 1) go to 280
c
            do 260 k = 1, j1
  260       x = x - a(j,k) * b(i,k)
c
  280       a(i,j) = x / dl(i)
  300 continue
c
      go to 1001
c     .......... set error -- b is not positive definite ..........
 1000 ierr = 7 * n + 1
 1001 return
      end
      subroutine reduc2(nm,n,a,b,dl,ierr)

c*********************************************************************72
c
cc REDUC2 reduces the eigenvalue problem A*B*x=lamdba*x to A*x=lambda*x.
c
c     this subroutine is a translation of the algol procedure reduc2,
c     num. math. 11, 99-110(1968) by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 303-314(1971).
c
c     this subroutine reduces the generalized symmetric eigenproblems
c     abx=(lambda)x or bay=(lambda)y, where b is positive definite,
c     to the standard symmetric eigenproblem using the cholesky
c     factorization of b.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrices a and b.  if the cholesky
c          factor l of b is already available, n should be prefixed
c          with a minus sign.
c
c        a and b contain the real symmetric input matrices.  only the
c          full upper triangles of the matrices need be supplied.  if
c          n is negative, the strict lower triangle of b contains,
c          instead, the strict lower triangle of its cholesky factor l.
c
c        dl contains, if n is negative, the diagonal elements of l.
c
c     on output
c
c        a contains in its full lower triangle the full lower triangle
c          of the symmetric matrix derived from the reduction to the
c          standard form.  the strict upper triangle of a is unaltered.
c
c        b contains in its strict lower triangle the strict lower
c          triangle of its cholesky factor l.  the full upper
c          triangle of b is unaltered.
c
c        dl contains the diagonal elements of l.
c
c        ierr is set to
c          zero       for normal return,
c          7*n+1      if b is not positive definite.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,k,n,i1,j1,nm,nn,ierr
      double precision a(nm,n),b(nm,n),dl(n)
      double precision x,y
c
      ierr = 0
      nn = iabs(n)
      if (n .lt. 0) go to 100
c     .......... form l in the arrays b and dl ..........
      do 80 i = 1, n
         i1 = i - 1
c
         do 80 j = i, n
            x = b(i,j)
            if (i .eq. 1) go to 40
c
            do 20 k = 1, i1
   20       x = x - b(i,k) * b(j,k)
c
   40       if (j .ne. i) go to 60
            if (x .le. 0.0d0) go to 1000
            y = dsqrt(x)
            dl(i) = y
            go to 80
   60       b(j,i) = x / y
   80 continue
c     .......... form the lower triangle of a*l
c                in the lower triangle of the array a ..........
  100 do 200 i = 1, nn
         i1 = i + 1
c
         do 200 j = 1, i
            x = a(j,i) * dl(j)
            if (j .eq. i) go to 140
            j1 = j + 1
c
            do 120 k = j1, i
  120       x = x + a(k,i) * b(k,j)
c
  140       if (i .eq. nn) go to 180
c
            do 160 k = i1, nn
  160       x = x + a(i,k) * b(k,j)
c
  180       a(i,j) = x
  200 continue
c     .......... pre-multiply by transpose(l) and overwrite ..........
      do 300 i = 1, nn
         i1 = i + 1
         y = dl(i)
c
         do 300 j = 1, i
            x = y * a(i,j)
            if (i .eq. nn) go to 280
c
            do 260 k = i1, nn
  260       x = x + a(k,j) * b(k,i)
c
  280       a(i,j) = x
  300 continue
c
      go to 1001
c     .......... set error -- b is not positive definite ..........
 1000 ierr = 7 * n + 1
 1001 return
      end
      subroutine rg(nm,n,a,wr,wi,matz,z,iv1,fv1,ierr)

c*********************************************************************72
c
cc RG computes eigenvalues and eigenvectors of a real general matrix.
c
c     this subroutine calls the recommended sequence of
c     routines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a real general matrix.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrix  a.
c
c        a  contains the real general matrix.
c
c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        wr  and  wi  contain the real and imaginary parts,
c        respectively, of the eigenvalues.  complex conjugate
c        pairs of eigenvalues appear consecutively with the
c        eigenvalue having the positive imaginary part first.
c
c        z  contains the real and imaginary parts of the eigenvectors
c        if matz is not zero.  if the j-th eigenvalue is real, the
c        j-th column of  z  contains its eigenvector.  if the j-th
c        eigenvalue is complex with positive imaginary part, the
c        j-th and (j+1)-th columns of  z  contain the real and
c        imaginary parts of its eigenvector.  the conjugate of this
c        vector is the eigenvector for the conjugate eigenvalue.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for hqr
c           and hqr2.  the normal completion code is zero.
c
c        iv1  and  fv1  are temporary storage arrays.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer n,nm,is1,is2,ierr,matz
      double precision a(nm,n),wr(n),wi(n),z(nm,n),fv1(n)
      integer iv1(n)

      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 call  balanc(nm,n,a,is1,is2,fv1)
      call  elmhes(nm,n,is1,is2,a,iv1)
      if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  hqr(nm,n,is1,is2,a,wr,wi,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 call  eltran(nm,n,is1,is2,a,iv1,z)
      call  hqr2(nm,n,is1,is2,a,wr,wi,z,ierr)
      if (ierr .ne. 0) go to 50
      call  balbak(nm,n,is1,is2,fv1,n,z)
   50 return
      end
      subroutine rgg(nm,n,a,b,alfr,alfi,beta,matz,z,ierr)

c*********************************************************************72
c
cc RGG computes eigenvalues/vectors for the generalized problem A*x = lambda*B*x.
c
c     this subroutine calls the recommended sequence of
c     routines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     for the real general generalized eigenproblem  ax = (lambda)bx.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrices  a  and  b.
c
c        a  contains a real general matrix.
c
c        b  contains a real general matrix.
c
c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        alfr  and  alfi  contain the real and imaginary parts,
c        respectively, of the numerators of the eigenvalues.
c
c        beta  contains the denominators of the eigenvalues,
c        which are thus given by the ratios  (alfr+i*alfi)/beta.
c        complex conjugate pairs of eigenvalues appear consecutively
c        with the eigenvalue having the positive imaginary part first.
c
c        z  contains the real and imaginary parts of the eigenvectors
c        if matz is not zero.  if the j-th eigenvalue is real, the
c        j-th column of  z  contains its eigenvector.  if the j-th
c        eigenvalue is complex with positive imaginary part, the
c        j-th and (j+1)-th columns of  z  contain the real and
c        imaginary parts of its eigenvector.  the conjugate of this
c        vector is the eigenvector for the conjugate eigenvalue.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for qzit.
c           the normal completion code is zero.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer n,nm,ierr,matz
      double precision a(nm,n),b(nm,n),alfr(n),alfi(n),beta(n),z(nm,n)
      logical tf

      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      tf = .false.
      call  qzhes(nm,n,a,b,tf,z)
      call  qzit(nm,n,a,b,0.0d0,tf,z,ierr)
      call  qzval(nm,n,a,b,alfr,alfi,beta,tf,z)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 tf = .true.
      call  qzhes(nm,n,a,b,tf,z)
      call  qzit(nm,n,a,b,0.0d0,tf,z,ierr)
      call  qzval(nm,n,a,b,alfr,alfi,beta,tf,z)
      if (ierr .ne. 0) go to 50
      call  qzvec(nm,n,a,b,alfr,alfi,beta,z)
   50 return
      end
      subroutine rs(nm,n,a,w,matz,z,fv1,fv2,ierr)

c*********************************************************************72
c
cc RS computes eigenvalues and eigenvectors of real symmetric matrix.
c
c     this subroutine calls the recommended sequence of
c     routines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a real symmetric matrix.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrix  a.
c
c        a  contains the real symmetric matrix.
c
c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        w  contains the eigenvalues in ascending order.
c
c        z  contains the eigenvectors if matz is not zero.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for tqlrat
c           and tql2.  the normal completion code is zero.
c
c        fv1  and  fv2  are temporary storage arrays.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer n,nm,ierr,matz
      double precision a(nm,n),w(n),z(nm,n),fv1(n),fv2(n)

      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  tred1(nm,n,a,w,fv1,fv2)
*  tqlrat encounters catastrophic underflow on the Vax
*     call  tqlrat(n,w,fv2,ierr)
      call  tql1(n,w,fv1,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 call  tred2(nm,n,a,w,fv1,z)
      call  tql2(nm,n,w,fv1,z,ierr)
   50 return
      end
      subroutine rsb(nm,n,mb,a,w,matz,z,fv1,fv2,ierr)

c*********************************************************************72
c
cc RSB computes eigenvalues and eigenvectors of a real symmetric band matrix.
c
c     this subroutine calls the recommended sequence of
c     routines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a real symmetric band matrix.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrix  a.
c
c        mb  is the half band width of the matrix, defined as the
c        number of adjacent diagonals, including the principal
c        diagonal, required to specify the non-zero portion of the
c        lower triangle of the matrix.
c
c        a  contains the lower triangle of the real symmetric
c        band matrix.  its lowest subdiagonal is stored in the
c        last  n+1-mb  positions of the first column, its next
c        subdiagonal in the last  n+2-mb  positions of the
c        second column, further subdiagonals similarly, and
c        finally its principal diagonal in the  n  positions
c        of the last column.  contents of storages not part
c        of the matrix are arbitrary.
c
c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        w  contains the eigenvalues in ascending order.
c
c        z  contains the eigenvectors if matz is not zero.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for tqlrat
c           and tql2.  the normal completion code is zero.
c
c        fv1  and  fv2  are temporary storage arrays.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer n,mb,nm,ierr,matz
      double precision a(nm,mb),w(n),z(nm,n),fv1(n),fv2(n)
      logical tf

      if (n .le. nm) go to 5
      ierr = 10 * n
      go to 50
    5 if (mb .gt. 0) go to 10
      ierr = 12 * n
      go to 50
   10 if (mb .le. n) go to 15
      ierr = 12 * n
      go to 50
c
   15 if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      tf = .false.
      call  bandr(nm,n,mb,a,w,fv1,fv2,tf,z)
      call  tqlrat(n,w,fv2,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 tf = .true.
      call  bandr(nm,n,mb,a,w,fv1,fv1,tf,z)
      call  tql2(nm,n,w,fv1,z,ierr)
   50 return
      end
      subroutine rsg(nm,n,a,b,w,matz,z,fv1,fv2,ierr)

c*********************************************************************72
c
cc RSG computes eigenvalues/vectors, A*x=lambda*B*x, A symmetric, B pos-def.
c
c     this subroutine calls the recommended sequence of
c     routines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     for the real symmetric generalized eigenproblem  ax = (lambda)bx.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrices  a  and  b.
c
c        a  contains a real symmetric matrix.
c
c        b  contains a positive definite real symmetric matrix.
c
c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        w  contains the eigenvalues in ascending order.
c
c        z  contains the eigenvectors if matz is not zero.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for tqlrat
c           and tql2.  the normal completion code is zero.
c
c        fv1  and  fv2  are temporary storage arrays.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer n,nm,ierr,matz
      double precision a(nm,n),b(nm,n),w(n),z(nm,n),fv1(n),fv2(n)

      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 call  reduc(nm,n,a,b,fv2,ierr)
      if (ierr .ne. 0) go to 50
      if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  tred1(nm,n,a,w,fv1,fv2)
      call  tqlrat(n,w,fv2,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 call  tred2(nm,n,a,w,fv1,z)
      call  tql2(nm,n,w,fv1,z,ierr)
      if (ierr .ne. 0) go to 50
      call  rebak(nm,n,b,fv2,n,z)
   50 return
      end
      subroutine rsgab(nm,n,a,b,w,matz,z,fv1,fv2,ierr)

c*********************************************************************72
c
cc RSGAB computes eigenvalues/vectors, A*B*x=lambda*x, A symmetric, B pos-def.
c
c     this subroutine calls the recommended sequence of
c     routines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     for the real symmetric generalized eigenproblem  abx = (lambda)x.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrices  a  and  b.
c
c        a  contains a real symmetric matrix.
c
c        b  contains a positive definite real symmetric matrix.
c
c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        w  contains the eigenvalues in ascending order.
c
c        z  contains the eigenvectors if matz is not zero.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for tqlrat
c           and tql2.  the normal completion code is zero.
c
c        fv1  and  fv2  are temporary storage arrays.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer n,nm,ierr,matz
      double precision a(nm,n),b(nm,n),w(n),z(nm,n),fv1(n),fv2(n)

      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 call  reduc2(nm,n,a,b,fv2,ierr)
      if (ierr .ne. 0) go to 50
      if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  tred1(nm,n,a,w,fv1,fv2)
      call  tqlrat(n,w,fv2,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 call  tred2(nm,n,a,w,fv1,z)
      call  tql2(nm,n,w,fv1,z,ierr)
      if (ierr .ne. 0) go to 50
      call  rebak(nm,n,b,fv2,n,z)
   50 return
      end
      subroutine rsgba(nm,n,a,b,w,matz,z,fv1,fv2,ierr)

c*********************************************************************72
c
cc RSGBA computes eigenvalues/vectors, B*A*x=lambda*x, A symmetric, B pos-def.
c
c     this subroutine calls the recommended sequence of
c     routines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     for the real symmetric generalized eigenproblem  bax = (lambda)x.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrices  a  and  b.
c
c        a  contains a real symmetric matrix.
c
c        b  contains a positive definite real symmetric matrix.
c
c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        w  contains the eigenvalues in ascending order.
c
c        z  contains the eigenvectors if matz is not zero.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for tqlrat
c           and tql2.  the normal completion code is zero.
c
c        fv1  and  fv2  are temporary storage arrays.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer n,nm,ierr,matz
      double precision a(nm,n),b(nm,n),w(n),z(nm,n),fv1(n),fv2(n)

      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 call  reduc2(nm,n,a,b,fv2,ierr)
      if (ierr .ne. 0) go to 50
      if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  tred1(nm,n,a,w,fv1,fv2)
      call  tqlrat(n,w,fv2,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 call  tred2(nm,n,a,w,fv1,z)
      call  tql2(nm,n,w,fv1,z,ierr)
      if (ierr .ne. 0) go to 50
      call  rebakb(nm,n,b,fv2,n,z)
   50 return
      end
      subroutine rsm(nm,n,a,w,m,z,fwork,iwork,ierr)

c*********************************************************************72
c
cc RSM computes eigenvalues, some eigenvectors, real symmetric matrix.
c
c     this subroutine calls the recommended sequence of
c     routines from the eigensystem subroutine package (eispack)
c     to find all of the eigenvalues and some of the eigenvectors
c     of a real symmetric matrix.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrix  a.
c
c        a  contains the real symmetric matrix.
c
c        m  the eigenvectors corresponding to the first m eigenvalues
c           are to be computed.
c           if m = 0 then no eigenvectors are computed.
c           if m = n then all of the eigenvectors are computed.
c
c     on output
c
c        w  contains all n eigenvalues in ascending order.
c
c        z  contains the orthonormal eigenvectors associated with
c           the first m eigenvalues.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for tqlrat,
c           imtqlv and tinvit.  the normal completion code is zero.
c
c        fwork  is a temporary storage array of dimension 8*n.
c
c        iwork  is an integer temporary storage array of dimension n.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer n,nm,m,iwork(n),ierr
      integer k1,k2,k3,k4,k5,k6,k7
      double precision a(nm,n),w(n),z(nm,m),fwork(1)

      ierr = 10 * n
      if (n .gt. nm .or. m .gt. nm) go to 50
      k1 = 1
      k2 = k1 + n
      k3 = k2 + n
      k4 = k3 + n
      k5 = k4 + n
      k6 = k5 + n
      k7 = k6 + n
      k8 = k7 + n
      if (m .gt. 0) go to 10
c     .......... find eigenvalues only ..........
      call  tred1(nm,n,a,w,fwork(k1),fwork(k2))
      call  tqlrat(n,w,fwork(k2),ierr)
      go to 50
c     .......... find all eigenvalues and m eigenvectors ..........
   10 call  tred1(nm,n,a,fwork(k1),fwork(k2),fwork(k3))
      call  imtqlv(n,fwork(k1),fwork(k2),fwork(k3),w,iwork,
     x             ierr,fwork(k4))
      call  tinvit(nm,n,fwork(k1),fwork(k2),fwork(k3),m,w,iwork,z,ierr,
     x             fwork(k4),fwork(k5),fwork(k6),fwork(k7),fwork(k8))
      call  trbak1(nm,n,a,fwork(k2),m,z)
   50 return
      end
      subroutine rsp(nm,n,nv,a,w,matz,z,fv1,fv2,ierr)

c*********************************************************************72
c
cc RSP computes eigenvalues and eigenvectors of real symmetric packed matrix.
c
c     this subroutine calls the recommended sequence of
c     routines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a real symmetric packed matrix.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrix  a.
c
c        nv  is an integer variable set equal to the
c        dimension of the array  a  as specified for
c        a  in the calling program.  nv  must not be
c        less than  n*(n+1)/2.
c
c        a  contains the lower triangle of the real symmetric
c        packed matrix stored row-wise.
c
c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        w  contains the eigenvalues in ascending order.
c
c        z  contains the eigenvectors if matz is not zero.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for tqlrat
c           and tql2.  the normal completion code is zero.
c
c        fv1  and  fv2  are temporary storage arrays.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,n,nm,nv,ierr,matz
      double precision a(nv),w(n),z(nm,n),fv1(n),fv2(n)
c

      if (n .le. nm) go to 5
      ierr = 10 * n
      go to 50
    5 if (nv .ge. (n * (n + 1)) / 2) go to 10
      ierr = 20 * n
      go to 50
c
   10 call  tred3(n,nv,a,w,fv1,fv2)
      if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  tqlrat(n,w,fv2,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 do 40 i = 1, n
c
         do 30 j = 1, n
            z(j,i) = 0.0d0
   30    continue
c
         z(i,i) = 1.0d0
   40 continue
c
      call  tql2(nm,n,w,fv1,z,ierr)
      if (ierr .ne. 0) go to 50
      call  trbak3(nm,n,nv,a,n,z)
   50 return
      end
      subroutine rspp ( n, nv, a, w, matz, z, ierr, m, type )

c*********************************************************************72
c
cc RSPP computes some eigenvalues/vectors, real symmetric packed matrix.
c
c  Discussion:
c
c    This routine calls the appropriate routines for the following problem:
c
c    Given a symmetric matrix A, which is stored in a packed mode, find
c    the M smallest or largest eigenvalues, and corresponding eigenvectors.
c
c    The routine RSP returns all eigenvalues and eigenvectors.
c
c  Modified:
c
c    09 February 2008
c
c  Reference:
c
c    James Wilkinson, Christian Reinsch,
c    Handbook for Automatic Computation,
c    Volume II, Linear Algebra, Part 2,
c    Springer, 1971,
c    ISBN: 0387054146,
c    LC: QA251.W67.
c
c    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
c    Matrix Eigensystem Routines, EISPACK Guide,
c    Lecture Notes in Computer Science, Volume 6,
c    Springer Verlag, 1976.
c
c  Parameters:
c
c    Input, integer N, the order of A, the number of rows and columns in the
c    original matrix.
c
c    Input, integer NV, is the of the array A as specified in the calling
c    program.  NV must not be less than N*(N+1)/2.
c
c    Input, double precision A((N*(N+1))/2), on input the lower triangle of the
c    real symmetric matrix, stored row-wise in the vector,
c    in the order A(1,1), / A(2,1), A(2,2), / A(3,1), A(3,2), A(3,3)/
c    and so on.
c
c    Output, double precision W(M), the eigenvalues requested.
c
c    Input, integer MATZ, is set to 0 if only eigenvalues are desired.
c    Otherwise it is set to any non-zero integer for both eigenvalues
c    and eigenvectors.
c
c    Output, double precision Z(N,M), the eigenvectors.
c
c    Output, integer IERR, error flag from RATQR.  IERR=0 on normal return.
c    IERR nonzero, in this case, means that the algorithm broke
c    down while computing an eigenvalue.
c
c    Input, integer M, the number of eigenvalues to be found.
c
c    Input, logical TYPE, set to .true. if the smallest eigenvalues
c    are to be found, or .false. if the largest ones are sought.
c
      implicit none

      integer m
      integer n
      integer nv

      double precision a(nv)
      double precision bd(n)
      double precision d(n)
      double precision e(n)
      double precision e2(n)
      double precision eps1
      integer idef
      integer ierr
      integer iwork(n)
      integer matz
      logical type
      double precision w(m)
      double precision work1(n)
      double precision work2(n)
      double precision work3(n)
      double precision work4(n)
      double precision work6(n)
      double precision z(n,m)
c
c  IDEF =
c    -1 if the matrix is known to be negative definite,
c    +1 if the matrix is known to be positive definite, or
c    0 otherwise.
c
      idef = 0
c
c  Reduce to symmetric tridiagonal form.
c
      call tred3 ( n, nv, a, d, e, e2 )
c
c  Find the eigenvalues.
c
      eps1 = 0.0D+00

      call ratqr ( n, eps1, d, e, e2, m, w, iwork,
     &  bd, type, idef, ierr )

      if ( ierr .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RSPP - Fatal error!'
        write ( *, '(a)' ) '  Error return from RATQR.'
        return
      end if
c
c  Find eigenvectors for the first M eigenvalues.
c
      if ( matz .ne. 0 ) then

        call tinvit ( n, n, d, e, e2, m, w, iwork, z, ierr,
     &  work1, work2, work3, work4, work6 )

        if ( ierr .ne. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'RSPP - Fatal error!'
          write ( *, '(a)' ) '  Error return from TINVIT.'
          return
        end if
c
c  Reverse the transformation.
c
        call trbak3 ( n, n, nv, a, m, z )

      end if

      return
      end
      subroutine rst(nm,n,w,e,matz,z,ierr)

c*********************************************************************72
c
cc RST computes eigenvalues/vectors, real symmetric tridiagonal matrix.
c
c     this subroutine calls the recommended sequence of
c     routines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a real symmetric tridiagonal matrix.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrix.
c
c        w  contains the diagonal elements of the real
c        symmetric tridiagonal matrix.
c
c        e  contains the subdiagonal elements of the matrix in
c        its last n-1 positions.  e(1) is arbitrary.
c
c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        w  contains the eigenvalues in ascending order.
c
c        z  contains the eigenvectors if matz is not zero.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for imtql1
c           and imtql2.  the normal completion code is zero.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,n,nm,ierr,matz
      double precision w(n),e(n),z(nm,n)
c
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  imtql1(n,w,e,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 do 40 i = 1, n
c
         do 30 j = 1, n
            z(j,i) = 0.0d0
   30    continue
c
         z(i,i) = 1.0d0
   40 continue
c
      call  imtql2(nm,n,w,e,z,ierr)
   50 return
      end
      subroutine rt(nm,n,a,w,matz,z,fv1,ierr)

c*********************************************************************72
c
cc RT computes eigenvalues/vectors, real sign-symmetric tridiagonal matrix.
c
c     this subroutine calls the recommended sequence of
c     routines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a special real tridiagonal matrix.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrix  a.
c
c        a  contains the special real tridiagonal matrix in its
c        first three columns.  the subdiagonal elements are stored
c        in the last  n-1  positions of the first column, the
c        diagonal elements in the second column, and the superdiagonal
c        elements in the first  n-1  positions of the third column.
c        elements  a(1,1)  and  a(n,3)  are arbitrary.
c
c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        w  contains the eigenvalues in ascending order.
c
c        z  contains the eigenvectors if matz is not zero.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for imtql1
c           and imtql2.  the normal completion code is zero.
c
c        fv1  is a temporary storage array.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer n,nm,ierr,matz
      double precision a(nm,3),w(n),z(nm,n),fv1(n)
c
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  figi(nm,n,a,w,fv1,fv1,ierr)
      if (ierr .gt. 0) go to 50
      call  imtql1(n,w,fv1,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 call  figi2(nm,n,a,w,fv1,z,ierr)
      if (ierr .ne. 0) go to 50
      call  imtql2(nm,n,w,fv1,z,ierr)
   50 return
      end
      subroutine svd(nm,m,n,a,w,matu,u,matv,v,ierr,rv1)

c*********************************************************************72
c
cc SVD computes the singular value decomposition for a real matrix.
c
c     this subroutine is a translation of the algol procedure svd,
c     num. math. 14, 403-420(1970) by golub and reinsch.
c     handbook for auto. comp., vol ii-linear algebra, 134-151(1971).
c
c     this subroutine determines the singular value decomposition
c          t
c     a=usv  of a real m by n rectangular matrix.  householder
c     bidiagonalization and a variant of the qr algorithm are used.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.  note that nm must be at least
c          as large as the maximum of m and n.
c
c        m is the number of rows of a (and u).
c
c        n is the number of columns of a (and u) and the order of v.
c
c        a contains the rectangular input matrix to be decomposed.
c
c        matu should be set to .true. if the u matrix in the
c          decomposition is desired, and to .false. otherwise.
c
c        matv should be set to .true. if the v matrix in the
c          decomposition is desired, and to .false. otherwise.
c
c     on output
c
c        a is unaltered (unless overwritten by u or v).
c
c        w contains the n (non-negative) singular values of a (the
c          diagonal elements of s).  they are unordered.  if an
c          error exit is made, the singular values should be correct
c          for indices ierr+1,ierr+2,...,n.
c
c        u contains the matrix u (orthogonal column vectors) of the
c          decomposition if matu has been set to .true.  otherwise
c          u is used as a temporary array.  u may coincide with a.
c          if an error exit is made, the columns of u corresponding
c          to indices of correct singular values should be correct.
c
c        v contains the matrix v (orthogonal) of the decomposition if
c          matv has been set to .true.  otherwise v is not referenced.
c          v may also coincide with a if u is not needed.  if an error
c          exit is made, the columns of v corresponding to indices of
c          correct singular values should be correct.
c
c        ierr is set to
c          zero       for normal return,
c          k          if the k-th singular value has not been
c                     determined after 30 iterations.
c
c        rv1 is a temporary storage array.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,k,l,m,n,ii,i1,kk,k1,ll,l1,mn,nm,its,ierr
      double precision a(nm,n),w(n),u(nm,n),v(nm,n),rv1(n)
      double precision c,f,g,h,s,x,y,z,tst1,tst2,scale,pythag
      logical matu,matv
c
      ierr = 0
c
      do 100 i = 1, m
c
         do 100 j = 1, n
            u(i,j) = a(i,j)
  100 continue
c     .......... householder reduction to bidiagonal form ..........
      g = 0.0d0
      scale = 0.0d0
      x = 0.0d0
c
      do 300 i = 1, n
         l = i + 1
         rv1(i) = scale * g
         g = 0.0d0
         s = 0.0d0
         scale = 0.0d0
         if (i .gt. m) go to 210
c
         do 120 k = i, m
  120    scale = scale + dabs(u(k,i))
c
         if (scale .eq. 0.0d0) go to 210
c
         do 130 k = i, m
            u(k,i) = u(k,i) / scale
            s = s + u(k,i)**2
  130    continue
c
         f = u(i,i)
         g = -dsign(dsqrt(s),f)
         h = f * g - s
         u(i,i) = f - g
         if (i .eq. n) go to 190
c
         do 150 j = l, n
            s = 0.0d0
c
            do 140 k = i, m
  140       s = s + u(k,i) * u(k,j)
c
            f = s / h
c
            do 150 k = i, m
               u(k,j) = u(k,j) + f * u(k,i)
  150    continue
c
  190    do 200 k = i, m
  200    u(k,i) = scale * u(k,i)
c
  210    w(i) = scale * g
         g = 0.0d0
         s = 0.0d0
         scale = 0.0d0
         if (i .gt. m .or. i .eq. n) go to 290
c
         do 220 k = l, n
  220    scale = scale + dabs(u(i,k))
c
         if (scale .eq. 0.0d0) go to 290
c
         do 230 k = l, n
            u(i,k) = u(i,k) / scale
            s = s + u(i,k)**2
  230    continue
c
         f = u(i,l)
         g = -dsign(dsqrt(s),f)
         h = f * g - s
         u(i,l) = f - g
c
         do 240 k = l, n
  240    rv1(k) = u(i,k) / h
c
         if (i .eq. m) go to 270
c
         do 260 j = l, m
            s = 0.0d0
c
            do 250 k = l, n
  250       s = s + u(j,k) * u(i,k)
c
            do 260 k = l, n
               u(j,k) = u(j,k) + s * rv1(k)
  260    continue
c
  270    do 280 k = l, n
  280    u(i,k) = scale * u(i,k)
c
  290    x = dmax1(x,dabs(w(i))+dabs(rv1(i)))
  300 continue
c     .......... accumulation of right-hand transformations ..........
      if (.not. matv) go to 410
c     .......... for i=n step -1 until 1 do -- ..........
      do 400 ii = 1, n
         i = n + 1 - ii
         if (i .eq. n) go to 390
         if (g .eq. 0.0d0) go to 360
c
         do 320 j = l, n
c     .......... double division avoids possible underflow ..........
  320    v(j,i) = (u(i,j) / u(i,l)) / g
c
         do 350 j = l, n
            s = 0.0d0
c
            do 340 k = l, n
  340       s = s + u(i,k) * v(k,j)
c
            do 350 k = l, n
               v(k,j) = v(k,j) + s * v(k,i)
  350    continue
c
  360    do 380 j = l, n
            v(i,j) = 0.0d0
            v(j,i) = 0.0d0
  380    continue
c
  390    v(i,i) = 1.0d0
         g = rv1(i)
         l = i
  400 continue
c     .......... accumulation of left-hand transformations ..........
  410 if (.not. matu) go to 510
c     ..........for i=min(m,n) step -1 until 1 do -- ..........
      mn = n
      if (m .lt. n) mn = m
c
      do 500 ii = 1, mn
         i = mn + 1 - ii
         l = i + 1
         g = w(i)
         if (i .eq. n) go to 430
c
         do 420 j = l, n
  420    u(i,j) = 0.0d0
c
  430    if (g .eq. 0.0d0) go to 475
         if (i .eq. mn) go to 460
c
         do 450 j = l, n
            s = 0.0d0
c
            do 440 k = l, m
  440       s = s + u(k,i) * u(k,j)
c     .......... double division avoids possible underflow ..........
            f = (s / u(i,i)) / g
c
            do 450 k = i, m
               u(k,j) = u(k,j) + f * u(k,i)
  450    continue
c
  460    do 470 j = i, m
  470    u(j,i) = u(j,i) / g
c
         go to 490
c
  475    do 480 j = i, m
  480    u(j,i) = 0.0d0
c
  490    u(i,i) = u(i,i) + 1.0d0
  500 continue
c     .......... diagonalization of the bidiagonal form ..........
  510 tst1 = x
c     .......... for k=n step -1 until 1 do -- ..........
      do 700 kk = 1, n
         k1 = n - kk
         k = k1 + 1
         its = 0
c     .......... test for splitting.
c                for l=k step -1 until 1 do -- ..........
  520    do 530 ll = 1, k
            l1 = k - ll
            l = l1 + 1
            tst2 = tst1 + dabs(rv1(l))
            if (tst2 .eq. tst1) go to 565
c     .......... rv1(1) is always zero, so there is no exit
c                through the bottom of the loop ..........
            tst2 = tst1 + dabs(w(l1))
            if (tst2 .eq. tst1) go to 540
  530    continue
c     .......... cancellation of rv1(l) if l greater than 1 ..........
  540    c = 0.0d0
         s = 1.0d0
c
         do 560 i = l, k
            f = s * rv1(i)
            rv1(i) = c * rv1(i)
            tst2 = tst1 + dabs(f)
            if (tst2 .eq. tst1) go to 565
            g = w(i)
            h = pythag(f,g)
            w(i) = h
            c = g / h
            s = -f / h
            if (.not. matu) go to 560
c
            do 550 j = 1, m
               y = u(j,l1)
               z = u(j,i)
               u(j,l1) = y * c + z * s
               u(j,i) = -y * s + z * c
  550       continue
c
  560    continue
c     .......... test for convergence ..........
  565    z = w(k)
         if (l .eq. k) go to 650
c     .......... shift from bottom 2 by 2 minor ..........
         if (its .eq. 30) go to 1000
         its = its + 1
         x = w(l)
         y = w(k1)
         g = rv1(k1)
         h = rv1(k)
         f = 0.5d0 * (((g + z) / h) * ((g - z) / y) + y / h - h / y)
         g = pythag(f,1.0d0)
         f = x - (z / x) * z + (h / x) * (y / (f + dsign(g,f)) - h)
c     .......... next qr transformation ..........
         c = 1.0d0
         s = 1.0d0
c
         do 600 i1 = l, k1
            i = i1 + 1
            g = rv1(i)
            y = w(i)
            h = s * g
            g = c * g
            z = pythag(f,h)
            rv1(i1) = z
            c = f / z
            s = h / z
            f = x * c + g * s
            g = -x * s + g * c
            h = y * s
            y = y * c
            if (.not. matv) go to 575
c
            do 570 j = 1, n
               x = v(j,i1)
               z = v(j,i)
               v(j,i1) = x * c + z * s
               v(j,i) = -x * s + z * c
  570       continue
c
  575       z = pythag(f,h)
            w(i1) = z
c     .......... rotation can be arbitrary if z is zero ..........
            if (z .eq. 0.0d0) go to 580
            c = f / z
            s = h / z
  580       f = c * g + s * y
            x = -s * g + c * y
            if (.not. matu) go to 600
c
            do 590 j = 1, m
               y = u(j,i1)
               z = u(j,i)
               u(j,i1) = y * c + z * s
               u(j,i) = -y * s + z * c
  590       continue
c
  600    continue
c
         rv1(l) = 0.0d0
         rv1(k) = f
         w(k) = x
         go to 520
c     .......... convergence ..........
  650    if (z .ge. 0.0d0) go to 700
c     .......... w(k) is made non-negative ..........
         w(k) = -z
         if (.not. matv) go to 700
c
         do 690 j = 1, n
  690    v(j,k) = -v(j,k)
c
  700 continue
c
      go to 1001
c     .......... set error -- no convergence to a
c                singular value after 30 iterations ..........
 1000 ierr = k
 1001 return
      end
      subroutine tinvit(nm,n,d,e,e2,m,w,ind,z,
     x                  ierr,rv1,rv2,rv3,rv4,rv6)

c*********************************************************************72
c
cc TINVIT computes eigenvectors from eigenvalues, real tridiagonal symmetric.
c
c     this subroutine is a translation of the inverse iteration tech-
c     nique in the algol procedure tristurm by peters and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 418-439(1971).
c
c     this subroutine finds those eigenvectors of a tridiagonal
c     symmetric matrix corresponding to specified eigenvalues,
c     using inverse iteration.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c        e2 contains the squares of the corresponding elements of e,
c          with zeros corresponding to negligible elements of e.
c          e(i) is considered negligible if it is not larger than
c          the product of the relative machine precision and the sum
c          of the magnitudes of d(i) and d(i-1).  e2(1) must contain
c          0.0d0 if the eigenvalues are in ascending order, or 2.0d0
c          if the eigenvalues are in descending order.  if  bisect,
c          tridib, or  imtqlv  has been used to find the eigenvalues,
c          their output e2 array is exactly what is expected here.
c
c        m is the number of specified eigenvalues.
c
c        w contains the m eigenvalues in ascending or descending order.
c
c        ind contains in its first m positions the submatrix indices
c          associated with the corresponding eigenvalues in w --
c          1 for eigenvalues belonging to the first submatrix from
c          the top, 2 for those belonging to the second submatrix, etc.
c
c     on output
c
c        all input arrays are unaltered.
c
c        z contains the associated set of orthonormal eigenvectors.
c          any vector which fails to converge is set to zero.
c
c        ierr is set to
c          zero       for normal return,
c          -r         if the eigenvector corresponding to the r-th
c                     eigenvalue fails to converge in 5 iterations.
c
c        rv1, rv2, rv3, rv4, and rv6 are temporary storage arrays.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,m,n,p,q,r,s,ii,ip,jj,nm,its,tag,ierr,group
      double precision d(n),e(n),e2(n),w(m),z(nm,m),
     x       rv1(n),rv2(n),rv3(n),rv4(n),rv6(n)
      double precision u,v,uk,xu,x0,x1,eps2,eps3,eps4,norm,order,epslon,
     x       pythag
      integer ind(m)

      ierr = 0
      if (m .eq. 0) go to 1001
      tag = 0
      order = 1.0d0 - e2(1)
      q = 0
c     .......... establish and process next submatrix ..........
  100 p = q + 1
c
      do 120 q = p, n
         if (q .eq. n) go to 140
         if (e2(q+1) .eq. 0.0d0) go to 140
  120 continue
c     .......... find vectors by inverse iteration ..........
  140 tag = tag + 1
      s = 0
c
      do 920 r = 1, m
         if (ind(r) .ne. tag) go to 920
         its = 1
         x1 = w(r)
         if (s .ne. 0) go to 510
c     .......... check for isolated root ..........
         xu = 1.0d0
         if (p .ne. q) go to 490
         rv6(p) = 1.0d0
         go to 870
  490    norm = dabs(d(p))
         ip = p + 1
c
         do 500 i = ip, q
  500    norm = dmax1(norm, dabs(d(i))+dabs(e(i)))
c     .......... eps2 is the criterion for grouping,
c                eps3 replaces zero pivots and equal
c                roots are modified by eps3,
c                eps4 is taken very small to avoid overflow ..........
         eps2 = 1.0d-3 * norm
         eps3 = epslon(norm)
         uk = q - p + 1
         eps4 = uk * eps3
         uk = eps4 / dsqrt(uk)
         s = p
  505    group = 0
         go to 520
c     .......... look for close or coincident roots ..........
  510    if (dabs(x1-x0) .ge. eps2) go to 505
         group = group + 1
         if (order * (x1 - x0) .le. 0.0d0) x1 = x0 + order * eps3
c     .......... elimination with interchanges and
c                initialization of vector ..........
  520    v = 0.0d0
c
         do 580 i = p, q
            rv6(i) = uk
            if (i .eq. p) go to 560
            if (dabs(e(i)) .lt. dabs(u)) go to 540
c     .......... warning -- a divide check may occur here if
c                e2 array has not been specified correctly ..........
            xu = u / e(i)
            rv4(i) = xu
            rv1(i-1) = e(i)
            rv2(i-1) = d(i) - x1
            rv3(i-1) = 0.0d0
            if (i .ne. q) rv3(i-1) = e(i+1)
            u = v - xu * rv2(i-1)
            v = -xu * rv3(i-1)
            go to 580
  540       xu = e(i) / u
            rv4(i) = xu
            rv1(i-1) = u
            rv2(i-1) = v
            rv3(i-1) = 0.0d0
  560       u = d(i) - x1 - xu * v
            if (i .ne. q) v = e(i+1)
  580    continue
c
         if (u .eq. 0.0d0) u = eps3
         rv1(q) = u
         rv2(q) = 0.0d0
         rv3(q) = 0.0d0
c     .......... back substitution
c                for i=q step -1 until p do -- ..........
  600    do 620 ii = p, q
            i = p + q - ii
            rv6(i) = (rv6(i) - u * rv2(i) - v * rv3(i)) / rv1(i)
            v = u
            u = rv6(i)
  620    continue
c     .......... orthogonalize with respect to previous
c                members of group ..........
         if (group .eq. 0) go to 700
         j = r
c
         do 680 jj = 1, group
  630       j = j - 1
            if (ind(j) .ne. tag) go to 630
            xu = 0.0d0
c
            do 640 i = p, q
  640       xu = xu + rv6(i) * z(i,j)
c
            do 660 i = p, q
  660       rv6(i) = rv6(i) - xu * z(i,j)
c
  680    continue
c
  700    norm = 0.0d0
c
         do 720 i = p, q
  720    norm = norm + dabs(rv6(i))
c
         if (norm .ge. 1.0d0) go to 840
c     .......... forward substitution ..........
         if (its .eq. 5) go to 830
         if (norm .ne. 0.0d0) go to 740
         rv6(s) = eps4
         s = s + 1
         if (s .gt. q) s = p
         go to 780
  740    xu = eps4 / norm
c
         do 760 i = p, q
  760    rv6(i) = rv6(i) * xu
c     .......... elimination operations on next vector
c                iterate ..........
  780    do 820 i = ip, q
            u = rv6(i)
c     .......... if rv1(i-1) .eq. e(i), a row interchange
c                was performed earlier in the
c                triangularization process ..........
            if (rv1(i-1) .ne. e(i)) go to 800
            u = rv6(i-1)
            rv6(i-1) = rv6(i)
  800       rv6(i) = u - rv4(i) * rv6(i-1)
  820    continue
c
         its = its + 1
         go to 600
c     .......... set error -- non-converged eigenvector ..........
  830    ierr = -r
         xu = 0.0d0
         go to 870
c     .......... normalize so that sum of squares is
c                1 and expand to full order ..........
  840    u = 0.0d0
c
         do 860 i = p, q
  860    u = pythag(u,rv6(i))
c
         xu = 1.0d0 / u
c
  870    do 880 i = 1, n
  880    z(i,r) = 0.0d0
c
         do 900 i = p, q
  900    z(i,r) = rv6(i) * xu
c
         x0 = x1
  920 continue
c
      if (q .lt. n) go to 100
 1001 return
      end
      subroutine tql1(n,d,e,ierr)

c*********************************************************************72
c
cc TQL1 computes all eigenvalues of a real symmetric tridiagonal matrix.
c
c     this subroutine is a translation of the algol procedure tql1,
c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
c     wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
c
c     this subroutine finds the eigenvalues of a symmetric
c     tridiagonal matrix by the ql method.
c
c     on input
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct and
c          ordered for indices 1,2,...ierr-1, but may not be
c          the smallest eigenvalues.
c
c        e has been destroyed.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,l,m,n,ii,l1,l2,mml,ierr
      double precision d(n),e(n)
      double precision c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
c

      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
c
      f = 0.0d0
      tst1 = 0.0d0
      e(n) = 0.0d0
c
      do 290 l = 1, n
         j = 0
         h = dabs(d(l)) + dabs(e(l))
         if (tst1 .lt. h) tst1 = h
c     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
c     .......... e(n) is always zero, so there is no exit
c                through the bottom of the loop ..........
  110    continue
c
  120    if (m .eq. l) go to 210
  130    if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = pythag(p,1.0d0)
         d(l) = e(l) / (p + dsign(r,p))
         d(l1) = e(l) * (p + dsign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
c
         do 140 i = l2, n
  140    d(i) = d(i) - h
c
  145    f = f + h
c     .......... ql transformation ..........
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
  200    continue
c
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + dabs(e(l))
         if (tst2 .gt. tst1) go to 130
  210    p = d(l) + f
c     .......... order eigenvalues ..........
         if (l .eq. 1) go to 250
c     .......... for i=l step -1 until 2 do -- ..........
         do 230 ii = 2, l
            i = l + 2 - ii
            if (p .ge. d(i-1)) go to 270
            d(i) = d(i-1)
  230    continue
c
  250    i = 1
  270    d(i) = p
  290 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end
      subroutine tql2(nm,n,d,e,z,ierr)

c*********************************************************************72
c
cc TQL2 computes all eigenvalues/vectors, real symmetric tridiagonal matrix.
c
c     this subroutine is a translation of the algol procedure tql2,
c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
c     wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a symmetric tridiagonal matrix by the ql method.
c     the eigenvectors of a full symmetric matrix can also
c     be found if  tred2  has been used to reduce this
c     full matrix to tridiagonal form.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c        z contains the transformation matrix produced in the
c          reduction by  tred2, if performed.  if the eigenvectors
c          of the tridiagonal matrix are desired, z must contain
c          the identity matrix.
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct but
c          unordered for indices 1,2,...,ierr-1.
c
c        e has been destroyed.
c
c        z contains orthonormal eigenvectors of the symmetric
c          tridiagonal (or full) matrix.  if an error exit is made,
c          z contains the eigenvectors associated with the stored
c          eigenvalues.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
      double precision d(n),e(n),z(nm,n)
      double precision c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag

      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
c
      f = 0.0d0
      tst1 = 0.0d0
      e(n) = 0.0d0
c
      do 240 l = 1, n
         j = 0
         h = dabs(d(l)) + dabs(e(l))
         if (tst1 .lt. h) tst1 = h
c     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
c     .......... e(n) is always zero, so there is no exit
c                through the bottom of the loop ..........
  110    continue
c
  120    if (m .eq. l) go to 220
  130    if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = pythag(p,1.0d0)
         d(l) = e(l) / (p + dsign(r,p))
         d(l1) = e(l) * (p + dsign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
c
         do 140 i = l2, n
  140    d(i) = d(i) - h
c
  145    f = f + h
c     .......... ql transformation ..........
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
c     .......... form vector ..........
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       continue
c
  200    continue
c
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + dabs(e(l))
         if (tst2 .gt. tst1) go to 130
  220    d(l) = d(l) + f
  240 continue

c     .......... order eigenvalues and eigenvectors ..........
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
c
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
c
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
c
         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue
c
  300 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end
      subroutine tqlrat(n,d,e2,ierr)

c*********************************************************************72
c
cc TQLRAT computes all eigenvalues of a real symmetric tridiagonal matrix.
c
C     This subroutine is a translation of the Algol procedure tqlrat,
C     Algorithm 464, Comm. ACM 16, 689(1973) by Reinsch.
C
C     This subroutine finds the eigenvalues of a symmetric
C     tridiagonal matrix by the rational QL method.
C
C     On input
C
C        N is the order of the matrix.
C
C        D contains the diagonal elements of the input matrix.
C
C        E2 contains the squares of the subdiagonal elements of the
C          input matrix in its last N-1 positions.  E2(1) is arbitrary.
C
C      On output
C
C        D contains the eigenvalues in ascending order.  If an
C          error exit is made, the eigenvalues are correct and
C          ordered for indices 1,2,...IERR-1, but may not be
C          the smallest eigenvalues.
C
C        E2 has been destroyed.
C
C        IERR is set to
C          zero       for normal return,
C          J          if the J-th eigenvalue has not been
C                     determined after 30 iterations.
C
C     Calls PYTHAG for  DSQRT(A*A + B*B) .
C
C     Questions and comments should be directed to Burton S. Garbow,
C     Mathematics and Computer Science Div, Argonne National Laboratory
C
C     This version dated August 1987.
C     Modified by C. Moler to fix underflow/overflow difficulties,
C     especially on the VAX and other machines where epslon(1.0d0)**2
C     nearly underflows.  See the loop involving statement 102 and
C     the two statements just before statement 200.
C
      integer i,j,l,m,n,ii,l1,mml,ierr
      double precision d(n),e2(n)
      double precision b,c,f,g,h,p,r,s,t,epslon,pythag

      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e2(i-1) = e2(i)
c
      f = 0.0d0
      t = 0.0d0
      e2(n) = 0.0d0
c
      do 290 l = 1, n
         j = 0
         h = dabs(d(l)) + dsqrt(e2(l))
         if (t .gt. h) go to 105
         t = h
         b = epslon(t)
         c = b * b
         if (c .ne. 0.0d0) go to 105
c        spliting tolerance underflowed.  look for larger value.
         do 102 i = l, n
            h = dabs(d(i)) + dsqrt(e2(i))
            if (h .gt. t) t = h
  102    continue
         b = epslon(t)
         c = b * b
c     .......... look for small squared sub-diagonal element ..........
  105    do 110 m = l, n
            if (e2(m) .le. c) go to 120
c     .......... e2(n) is always zero, so there is no exit
c                through the bottom of the loop ..........
  110    continue
c
  120    if (m .eq. l) go to 210
  130    if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         l1 = l + 1
         s = dsqrt(e2(l))
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * s)
         r = pythag(p,1.0d0)
         d(l) = s / (p + dsign(r,p))
         h = g - d(l)
c
         do 140 i = l1, n
  140    d(i) = d(i) - h
c
         f = f + h
c     .......... rational ql transformation ..........
         g = d(m)
         if (g .eq. 0.0d0) g = b
         h = g
         s = 0.0d0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            i = m - ii
            p = g * h
            r = p + e2(i)
            e2(i+1) = s * r
            s = e2(i) / r
            d(i+1) = h + s * (h + d(i))
            g = d(i) - e2(i) / g
c           avoid division by zero on next pass
            if (g .eq. 0.0d0) g = epslon(d(i))
            h = g * (p / r)
  200    continue
c
         e2(l) = s * g
         d(l) = h
c     .......... guard against underflow in convergence test ..........
         if (h .eq. 0.0d0) go to 210
         if (dabs(e2(l)) .le. dabs(c/h)) go to 210
         e2(l) = h * e2(l)
         if (e2(l) .ne. 0.0d0) go to 130
  210    p = d(l) + f
c     .......... order eigenvalues ..........
         if (l .eq. 1) go to 250
c     .......... for i=l step -1 until 2 do -- ..........
         do 230 ii = 2, l
            i = l + 2 - ii
            if (p .ge. d(i-1)) go to 270
            d(i) = d(i-1)
  230    continue
c
  250    i = 1
  270    d(i) = p
  290 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end
      subroutine trbak1(nm,n,a,e,m,z)

c*********************************************************************72
c
cc TRBAK1 determines eigenvectors by undoing the TRED1 transformation.
c
c     this subroutine is a translation of the algol procedure trbak1,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine forms the eigenvectors of a real symmetric
c     matrix by back transforming those of the corresponding
c     symmetric tridiagonal matrix determined by  tred1.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        a contains information about the orthogonal trans-
c          formations used in the reduction by  tred1
c          in its strict lower triangle.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is arbitrary.
c
c        m is the number of eigenvectors to be back transformed.
c
c        z contains the eigenvectors to be back transformed
c          in its first m columns.
c
c     on output
c
c        z contains the transformed eigenvectors
c          in its first m columns.
c
c     note that trbak1 preserves vector euclidean norms.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,k,l,m,n,nm
      double precision a(nm,n),e(n),z(nm,m)
      double precision s

      if (m .eq. 0) go to 200
      if (n .eq. 1) go to 200
c
      do 140 i = 2, n
         l = i - 1
         if (e(i) .eq. 0.0d0) go to 140
c
         do 130 j = 1, m
            s = 0.0d0
c
            do 110 k = 1, l
  110       s = s + a(i,k) * z(k,j)
c     .......... divisor below is negative of h formed in tred1.
c                double division avoids possible underflow ..........
            s = (s / a(i,l)) / e(i)
c
            do 120 k = 1, l
  120       z(k,j) = z(k,j) + s * a(i,k)
c
  130    continue
c
  140 continue
c
  200 return
      end
      subroutine trbak3(nm,n,nv,a,m,z)

c*********************************************************************72
c
cc TRBAK3 determines eigenvectors by undoing the TRED3 transformation.
c
c     this subroutine is a translation of the algol procedure trbak3,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine forms the eigenvectors of a real symmetric
c     matrix by back transforming those of the corresponding
c     symmetric tridiagonal matrix determined by  tred3.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        nv must be set to the dimension of the array parameter a
c          as declared in the calling program dimension statement.
c
c        a contains information about the orthogonal transformations
c          used in the reduction by  tred3  in its first
c          n*(n+1)/2 positions.
c
c        m is the number of eigenvectors to be back transformed.
c
c        z contains the eigenvectors to be back transformed
c          in its first m columns.
c
c     on output
c
c        z contains the transformed eigenvectors
c          in its first m columns.
c
c     note that trbak3 preserves vector euclidean norms.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,k,l,m,n,ik,iz,nm,nv
      double precision a(nv),z(nm,m)
      double precision h,s

      if (m .eq. 0) go to 200
      if (n .eq. 1) go to 200
c
      do 140 i = 2, n
         l = i - 1
         iz = (i * l) / 2
         ik = iz + i
         h = a(ik)
         if (h .eq. 0.0d0) go to 140
c
         do 130 j = 1, m
            s = 0.0d0
            ik = iz
c
            do 110 k = 1, l
               ik = ik + 1
               s = s + a(ik) * z(k,j)
  110       continue
c     .......... double division avoids possible underflow ..........
            s = (s / h) / h
            ik = iz
c
            do 120 k = 1, l
               ik = ik + 1
               z(k,j) = z(k,j) - s * a(ik)
  120       continue
c
  130    continue
c
  140 continue
c
  200 return
      end
      subroutine tred1(nm,n,a,d,e,e2)

c*********************************************************************72
c
cc TRED1 transforms a real symmetric matrix to symmetric tridiagonal form.
c
c     this subroutine is a translation of the algol procedure tred1,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a real symmetric matrix
c     to a symmetric tridiagonal matrix using
c     orthogonal similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        a contains the real symmetric input matrix.  only the
c          lower triangle of the matrix need be supplied.
c
c     on output
c
c        a contains information about the orthogonal trans-
c          formations used in the reduction in its strict lower
c          triangle.  the full upper triangle of a is unaltered.
c
c        d contains the diagonal elements of the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        e2 contains the squares of the corresponding elements of e.
c          e2 may coincide with e if the squares are not needed.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,k,l,n,ii,nm,jp1
      double precision a(nm,n),d(n),e(n),e2(n)
      double precision f,g,h,scale

      do 100 i = 1, n
         d(i) = a(n,i)
         a(n,i) = a(i,i)
  100 continue
c     .......... for i=n step -1 until 1 do -- ..........
      do 300 ii = 1, n
         i = n + 1 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 1) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + dabs(d(k))
c
         if (scale .ne. 0.0d0) go to 140
c
         do 125 j = 1, l
            d(j) = a(l,j)
            a(l,j) = a(i,j)
            a(i,j) = 0.0d0
  125    continue
c
  130    e(i) = 0.0d0
         e2(i) = 0.0d0
         go to 300
c
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
c
         e2(i) = scale * scale * h
         f = d(l)
         g = -dsign(dsqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
         if (l .eq. 1) go to 285
c     .......... form a*u ..........
         do 170 j = 1, l
  170    e(j) = 0.0d0
c
         do 240 j = 1, l
            f = d(j)
            g = e(j) + a(j,j) * f
            jp1 = j + 1
            if (l .lt. jp1) go to 220
c
            do 200 k = jp1, l
               g = g + a(k,j) * d(k)
               e(k) = e(k) + a(k,j) * f
  200       continue
c
  220       e(j) = g
  240    continue
c     .......... form p ..........
         f = 0.0d0
c
         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue
c
         h = f / (h + h)
c     .......... form q ..........
         do 250 j = 1, l
  250    e(j) = e(j) - h * d(j)
c     .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)
c
            do 260 k = j, l
  260       a(k,j) = a(k,j) - f * e(k) - g * d(k)
c
  280    continue
c
  285    do 290 j = 1, l
            f = d(j)
            d(j) = a(l,j)
            a(l,j) = a(i,j)
            a(i,j) = f * scale
  290    continue
c
  300 continue
c
      return
      end
      subroutine tred2(nm,n,a,d,e,z)

c*********************************************************************72
c
cc TRED2 transforms a real symmetric matrix to symmetric tridiagonal form.
c
c     this subroutine is a translation of the algol procedure tred2,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a real symmetric matrix to a
c     symmetric tridiagonal matrix using and accumulating
c     orthogonal similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        a contains the real symmetric input matrix.  only the
c          lower triangle of the matrix need be supplied.
c
c     on output
c
c        d contains the diagonal elements of the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        z contains the orthogonal transformation matrix
c          produced in the reduction.
c
c        a and z may coincide.  if distinct, a is unaltered.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,k,l,n,ii,nm,jp1
      double precision a(nm,n),d(n),e(n),z(nm,n)
      double precision f,g,h,hh,scale

      do 100 i = 1, n
c
         do 80 j = i, n
   80    z(j,i) = a(j,i)
c
         d(i) = a(n,i)
  100 continue
c
      if (n .eq. 1) go to 510
c     .......... for i=n step -1 until 2 do -- ..........
      do 300 ii = 2, n
         i = n + 2 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 2) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + dabs(d(k))
c
         if (scale .ne. 0.0d0) go to 140
  130    e(i) = d(l)
c
         do 135 j = 1, l
            d(j) = z(l,j)
            z(i,j) = 0.0d0
            z(j,i) = 0.0d0
  135    continue
c
         go to 290
c
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
c
         f = d(l)
         g = -dsign(dsqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
c     .......... form a*u ..........
         do 170 j = 1, l
  170    e(j) = 0.0d0
c
         do 240 j = 1, l
            f = d(j)
            z(j,i) = f
            g = e(j) + z(j,j) * f
            jp1 = j + 1
            if (l .lt. jp1) go to 220
c
            do 200 k = jp1, l
               g = g + z(k,j) * d(k)
               e(k) = e(k) + z(k,j) * f
  200       continue
c
  220       e(j) = g
  240    continue
c     .......... form p ..........
         f = 0.0d0
c
         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue
c
         hh = f / (h + h)
c     .......... form q ..........
         do 250 j = 1, l
  250    e(j) = e(j) - hh * d(j)
c     .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)
c
            do 260 k = j, l
  260       z(k,j) = z(k,j) - f * e(k) - g * d(k)
c
            d(j) = z(l,j)
            z(i,j) = 0.0d0
  280    continue
c
  290    d(i) = h
  300 continue
c     .......... accumulation of transformation matrices ..........
      do 500 i = 2, n
         l = i - 1
         z(n,l) = z(l,l)
         z(l,l) = 1.0d0
         h = d(i)
         if (h .eq. 0.0d0) go to 380
c
         do 330 k = 1, l
  330    d(k) = z(k,i) / h
c
         do 360 j = 1, l
            g = 0.0d0
c
            do 340 k = 1, l
  340       g = g + z(k,i) * z(k,j)
c
            do 360 k = 1, l
               z(k,j) = z(k,j) - g * d(k)
  360    continue
c
  380    do 400 k = 1, l
  400    z(k,i) = 0.0d0
c
  500 continue
c
  510 do 520 i = 1, n
         d(i) = z(n,i)
         z(n,i) = 0.0d0
  520 continue
c
      z(n,n) = 1.0d0
      e(1) = 0.0d0
      return
      end
      subroutine tred3(n,nv,a,d,e,e2)

c*********************************************************************72
c
cc TRED3 transforms a real symmetric packed matrix to symmetric tridiagonal form.
c
c     this subroutine is a translation of the algol procedure tred3,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a real symmetric matrix, stored as
c     a one-dimensional array, to a symmetric tridiagonal matrix
c     using orthogonal similarity transformations.
c
c     on input
c
c        n is the order of the matrix.
c
c        nv must be set to the dimension of the array parameter a
c          as declared in the calling program dimension statement.
c
c        a contains the lower triangle of the real symmetric
c          input matrix, stored row-wise as a one-dimensional
c          array, in its first n*(n+1)/2 positions.
c
c     on output
c
c        a contains information about the orthogonal
c          transformations used in the reduction.
c
c        d contains the diagonal elements of the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        e2 contains the squares of the corresponding elements of e.
c          e2 may coincide with e if the squares are not needed.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,k,l,n,ii,iz,jk,nv,jm1
      double precision a(nv),d(n),e(n),e2(n)
      double precision f,g,h,hh,scale
c
c     .......... for i=n step -1 until 1 do -- ..........
      do 300 ii = 1, n
         i = n + 1 - ii
         l = i - 1
         iz = (i * l) / 2
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 1) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
            iz = iz + 1
            d(k) = a(iz)
            scale = scale + dabs(d(k))
  120    continue
c
         if (scale .ne. 0.0d0) go to 140
  130    e(i) = 0.0d0
         e2(i) = 0.0d0
         go to 290
c
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
c
         e2(i) = scale * scale * h
         f = d(l)
         g = -dsign(dsqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
         a(iz) = scale * d(l)
         if (l .eq. 1) go to 290
         jk = 1
c
         do 240 j = 1, l
            f = d(j)
            g = 0.0d0
            jm1 = j - 1
            if (jm1 .lt. 1) go to 220
c
            do 200 k = 1, jm1
               g = g + a(jk) * d(k)
               e(k) = e(k) + a(jk) * f
               jk = jk + 1
  200       continue
c
  220       e(j) = g + a(jk) * f
            jk = jk + 1
  240    continue
c     .......... form p ..........
         f = 0.0d0
c
         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue
c
         hh = f / (h + h)
c     .......... form q ..........
         do 250 j = 1, l
  250    e(j) = e(j) - hh * d(j)
c
         jk = 1
c     .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)
c
            do 260 k = 1, j
               a(jk) = a(jk) - f * e(k) - g * d(k)
               jk = jk + 1
  260       continue
c
  280    continue
c
  290    d(i) = a(iz+1)
         a(iz+1) = scale * dsqrt(h)
  300 continue
c
      return
      end
      subroutine tridib(n,eps1,d,e,e2,lb,ub,m11,m,w,ind,ierr,rv4,rv5)

c*********************************************************************72
c
cc TRIDIB computes some eigenvalues of a real symmetric tridiagonal matrix.
c
c     this subroutine is a translation of the algol procedure bisect,
c     num. math. 9, 386-393(1967) by barth, martin, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 249-256(1971).
c
c     this subroutine finds those eigenvalues of a tridiagonal
c     symmetric matrix between specified boundary indices,
c     using bisection.
c
c     on input
c
c        n is the order of the matrix.
c
c        eps1 is an absolute error tolerance for the computed
c          eigenvalues.  if the input eps1 is non-positive,
c          it is reset for each submatrix to a default value,
c          namely, minus the product of the relative machine
c          precision and the 1-norm of the submatrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c        e2 contains the squares of the corresponding elements of e.
c          e2(1) is arbitrary.
c
c        m11 specifies the lower boundary index for the desired
c          eigenvalues.
c
c        m specifies the number of eigenvalues desired.  the upper
c          boundary index m22 is then obtained as m22=m11+m-1.
c
c     on output
c
c        eps1 is unaltered unless it has been reset to its
c          (last) default value.
c
c        d and e are unaltered.
c
c        elements of e2, corresponding to elements of e regarded
c          as negligible, have been replaced by zero causing the
c          matrix to split into a direct sum of submatrices.
c          e2(1) is also set to zero.
c
c        lb and ub define an interval containing exactly the desired
c          eigenvalues.
c
c        w contains, in its first m positions, the eigenvalues
c          between indices m11 and m22 in ascending order.
c
c        ind contains in its first m positions the submatrix indices
c          associated with the corresponding eigenvalues in w --
c          1 for eigenvalues belonging to the first submatrix from
c          the top, 2 for those belonging to the second submatrix, etc..
c
c        ierr is set to
c          zero       for normal return,
c          3*n+1      if multiple eigenvalues at index m11 make
c                     unique selection impossible,
c          3*n+2      if multiple eigenvalues at index m22 make
c                     unique selection impossible.
c
c        rv4 and rv5 are temporary storage arrays.
c
c     note that subroutine tql1, imtql1, or tqlrat is generally faster
c     than tridib, if more than n/4 eigenvalues are to be found.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,k,l,m,n,p,q,r,s,ii,m1,m2,m11,m22,tag,ierr,isturm
      double precision d(n),e(n),e2(n),w(m),rv4(n),rv5(n)
      double precision u,v,lb,t1,t2,ub,xu,x0,x1,eps1,tst1,tst2,epslon
      integer ind(m)

      ierr = 0
      tag = 0
      xu = d(1)
      x0 = d(1)
      u = 0.0d0
c     .......... look for small sub-diagonal entries and determine an
c                interval containing all the eigenvalues ..........
      do 40 i = 1, n
         x1 = u
         u = 0.0d0
         if (i .ne. n) u = dabs(e(i+1))
         xu = dmin1(d(i)-(x1+u),xu)
         x0 = dmax1(d(i)+(x1+u),x0)
         if (i .eq. 1) go to 20
         tst1 = dabs(d(i)) + dabs(d(i-1))
         tst2 = tst1 + dabs(e(i))
         if (tst2 .gt. tst1) go to 40
   20    e2(i) = 0.0d0
   40 continue
c
      x1 = n
      x1 = x1 * epslon(dmax1(dabs(xu),dabs(x0)))
      xu = xu - x1
      t1 = xu
      x0 = x0 + x1
      t2 = x0
c     .......... determine an interval containing exactly
c                the desired eigenvalues ..........
      p = 1
      q = n
      m1 = m11 - 1
      if (m1 .eq. 0) go to 75
      isturm = 1
   50 v = x1
      x1 = xu + (x0 - xu) * 0.5d0
      if (x1 .eq. v) go to 980
      go to 320
   60 if (s - m1) 65, 73, 70
   65 xu = x1
      go to 50
   70 x0 = x1
      go to 50
   73 xu = x1
      t1 = x1
   75 m22 = m1 + m
      if (m22 .eq. n) go to 90
      x0 = t2
      isturm = 2
      go to 50
   80 if (s - m22) 65, 85, 70
   85 t2 = x1
   90 q = 0
      r = 0
c     .......... establish and process next submatrix, refining
c                interval by the gerschgorin bounds ..........
  100 if (r .eq. m) go to 1001
      tag = tag + 1
      p = q + 1
      xu = d(p)
      x0 = d(p)
      u = 0.0d0
c
      do 120 q = p, n
         x1 = u
         u = 0.0d0
         v = 0.0d0
         if (q .eq. n) go to 110
         u = dabs(e(q+1))
         v = e2(q+1)
  110    xu = dmin1(d(q)-(x1+u),xu)
         x0 = dmax1(d(q)+(x1+u),x0)
         if (v .eq. 0.0d0) go to 140
  120 continue
c
  140 x1 = epslon(dmax1(dabs(xu),dabs(x0)))
      if (eps1 .le. 0.0d0) eps1 = -x1
      if (p .ne. q) go to 180
c     .......... check for isolated root within interval ..........
      if (t1 .gt. d(p) .or. d(p) .ge. t2) go to 940
      m1 = p
      m2 = p
      rv5(p) = d(p)
      go to 900
  180 x1 = x1 * (q - p + 1)
      lb = dmax1(t1,xu-x1)
      ub = dmin1(t2,x0+x1)
      x1 = lb
      isturm = 3
      go to 320
  200 m1 = s + 1
      x1 = ub
      isturm = 4
      go to 320
  220 m2 = s
      if (m1 .gt. m2) go to 940
c     .......... find roots by bisection ..........
      x0 = ub
      isturm = 5
c
      do 240 i = m1, m2
         rv5(i) = ub
         rv4(i) = lb
  240 continue
c     .......... loop for k-th eigenvalue
c                for k=m2 step -1 until m1 do --
c                (-do- not used to legalize -computed go to-) ..........
      k = m2
  250    xu = lb
c     .......... for i=k step -1 until m1 do -- ..........
         do 260 ii = m1, k
            i = m1 + k - ii
            if (xu .ge. rv4(i)) go to 260
            xu = rv4(i)
            go to 280
  260    continue
c
  280    if (x0 .gt. rv5(k)) x0 = rv5(k)
c     .......... next bisection step ..........
  300    x1 = (xu + x0) * 0.5d0
         if ((x0 - xu) .le. dabs(eps1)) go to 420
         tst1 = 2.0d0 * (dabs(xu) + dabs(x0))
         tst2 = tst1 + (x0 - xu)
         if (tst2 .eq. tst1) go to 420
c     .......... in-line procedure for sturm sequence ..........
  320    s = p - 1
         u = 1.0d0
c
         do 340 i = p, q
            if (u .ne. 0.0d0) go to 325
            v = dabs(e(i)) / epslon(1.0d0)
            if (e2(i) .eq. 0.0d0) v = 0.0d0
            go to 330
  325       v = e2(i) / u
  330       u = d(i) - x1 - v
            if (u .lt. 0.0d0) s = s + 1
  340    continue
c
         go to (60,80,200,220,360), isturm
c     .......... refine intervals ..........
  360    if (s .ge. k) go to 400
         xu = x1
         if (s .ge. m1) go to 380
         rv4(m1) = x1
         go to 300
  380    rv4(s+1) = x1
         if (rv5(s) .gt. x1) rv5(s) = x1
         go to 300
  400    x0 = x1
         go to 300
c     .......... k-th eigenvalue found ..........
  420    rv5(k) = x1
      k = k - 1
      if (k .ge. m1) go to 250
c     .......... order eigenvalues tagged with their
c                submatrix associations ..........
  900 s = r
      r = r + m2 - m1 + 1
      j = 1
      k = m1
c
      do 920 l = 1, r
         if (j .gt. s) go to 910
         if (k .gt. m2) go to 940
         if (rv5(k) .ge. w(l)) go to 915
c
         do 905 ii = j, s
            i = l + s - ii
            w(i+1) = w(i)
            ind(i+1) = ind(i)
  905    continue
c
  910    w(l) = rv5(k)
         ind(l) = tag
         k = k + 1
         go to 920
  915    j = j + 1
  920 continue
c
  940 if (q .lt. n) go to 100
      go to 1001
c     .......... set error -- interval cannot be found containing
c                exactly the desired eigenvalues ..........
  980 ierr = 3 * n + isturm
 1001 lb = t1
      ub = t2
      return
      end
      subroutine tsturm(nm,n,eps1,d,e,e2,lb,ub,mm,m,w,z,
     x                  ierr,rv1,rv2,rv3,rv4,rv5,rv6)

c*********************************************************************72
c
cc TSTURM computes some eigenvalues/vectors, real symmetric tridiagonal matrix.
c
c     this subroutine is a translation of the algol procedure tristurm
c     by peters and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 418-439(1971).
c
c     this subroutine finds those eigenvalues of a tridiagonal
c     symmetric matrix which lie in a specified interval and their
c     associated eigenvectors, using bisection and inverse iteration.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        eps1 is an absolute error tolerance for the computed
c          eigenvalues.  it should be chosen commensurate with
c          relative perturbations in the matrix elements of the
c          order of the relative machine precision.  if the
c          input eps1 is non-positive, it is reset for each
c          submatrix to a default value, namely, minus the
c          product of the relative machine precision and the
c          1-norm of the submatrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c        e2 contains the squares of the corresponding elements of e.
c          e2(1) is arbitrary.
c
c        lb and ub define the interval to be searched for eigenvalues.
c          if lb is not less than ub, no eigenvalues will be found.
c
c        mm should be set to an upper bound for the number of
c          eigenvalues in the interval.  warning. if more than
c          mm eigenvalues are determined to lie in the interval,
c          an error return is made with no values or vectors found.
c
c     on output
c
c        eps1 is unaltered unless it has been reset to its
c          (last) default value.
c
c        d and e are unaltered.
c
c        elements of e2, corresponding to elements of e regarded
c          as negligible, have been replaced by zero causing the
c          matrix to split into a direct sum of submatrices.
c          e2(1) is also set to zero.
c
c        m is the number of eigenvalues determined to lie in (lb,ub).
c
c        w contains the m eigenvalues in ascending order if the matrix
c          does not split.  if the matrix splits, the eigenvalues are
c          in ascending order for each submatrix.  if a vector error
c          exit is made, w contains those values already found.
c
c        z contains the associated set of orthonormal eigenvectors.
c          if an error exit is made, z contains those vectors
c          already found.
c
c        ierr is set to
c          zero       for normal return,
c          3*n+1      if m exceeds mm.
c          4*n+r      if the eigenvector corresponding to the r-th
c                     eigenvalue fails to converge in 5 iterations.
c
c        rv1, rv2, rv3, rv4, rv5, and rv6 are temporary storage arrays.
c
c     the algol procedure sturmcnt contained in tristurm
c     appears in tsturm in-line.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,k,m,n,p,q,r,s,ii,ip,jj,mm,m1,m2,nm,its,
     x        ierr,group,isturm
      double precision d(n),e(n),e2(n),w(mm),z(nm,mm),
     x       rv1(n),rv2(n),rv3(n),rv4(n),rv5(n),rv6(n)
      double precision u,v,lb,t1,t2,ub,uk,xu,x0,x1,eps1,eps2,eps3,eps4,
     x       norm,tst1,tst2,epslon,pythag
c
      ierr = 0
      t1 = lb
      t2 = ub
c     .......... look for small sub-diagonal entries ..........
      do 40 i = 1, n
         if (i .eq. 1) go to 20
         tst1 = dabs(d(i)) + dabs(d(i-1))
         tst2 = tst1 + dabs(e(i))
         if (tst2 .gt. tst1) go to 40
   20    e2(i) = 0.0d0
   40 continue
c     .......... determine the number of eigenvalues
c                in the interval ..........
      p = 1
      q = n
      x1 = ub
      isturm = 1
      go to 320
   60 m = s
      x1 = lb
      isturm = 2
      go to 320
   80 m = m - s
      if (m .gt. mm) go to 980
      q = 0
      r = 0
c     .......... establish and process next submatrix, refining
c                interval by the gerschgorin bounds ..........
  100 if (r .eq. m) go to 1001
      p = q + 1
      xu = d(p)
      x0 = d(p)
      u = 0.0d0
c
      do 120 q = p, n
         x1 = u
         u = 0.0d0
         v = 0.0d0
         if (q .eq. n) go to 110
         u = dabs(e(q+1))
         v = e2(q+1)
  110    xu = dmin1(d(q)-(x1+u),xu)
         x0 = dmax1(d(q)+(x1+u),x0)
         if (v .eq. 0.0d0) go to 140
  120 continue
c
  140 x1 = epslon(dmax1(dabs(xu),dabs(x0)))
      if (eps1 .le. 0.0d0) eps1 = -x1
      if (p .ne. q) go to 180
c     .......... check for isolated root within interval ..........
      if (t1 .gt. d(p) .or. d(p) .ge. t2) go to 940
      r = r + 1
c
      do 160 i = 1, n
  160 z(i,r) = 0.0d0
c
      w(r) = d(p)
      z(p,r) = 1.0d0
      go to 940
  180 u = q-p+1
      x1 = u * x1
      lb = dmax1(t1,xu-x1)
      ub = dmin1(t2,x0+x1)
      x1 = lb
      isturm = 3
      go to 320
  200 m1 = s + 1
      x1 = ub
      isturm = 4
      go to 320
  220 m2 = s
      if (m1 .gt. m2) go to 940
c     .......... find roots by bisection ..........
      x0 = ub
      isturm = 5
c
      do 240 i = m1, m2
         rv5(i) = ub
         rv4(i) = lb
  240 continue
c     .......... loop for k-th eigenvalue
c                for k=m2 step -1 until m1 do --
c                (-do- not used to legalize -computed go to-) ..........
      k = m2
  250    xu = lb
c     .......... for i=k step -1 until m1 do -- ..........
         do 260 ii = m1, k
            i = m1 + k - ii
            if (xu .ge. rv4(i)) go to 260
            xu = rv4(i)
            go to 280
  260    continue
c
  280    if (x0 .gt. rv5(k)) x0 = rv5(k)
c     .......... next bisection step ..........
  300    x1 = (xu + x0) * 0.5d0
         if ((x0 - xu) .le. dabs(eps1)) go to 420
         tst1 = 2.0d0 * (dabs(xu) + dabs(x0))
         tst2 = tst1 + (x0 - xu)
         if (tst2 .eq. tst1) go to 420
c     .......... in-line procedure for sturm sequence ..........
  320    s = p - 1
         u = 1.0d0
c
         do 340 i = p, q
            if (u .ne. 0.0d0) go to 325
            v = dabs(e(i)) / epslon(1.0d0)
            if (e2(i) .eq. 0.0d0) v = 0.0d0
            go to 330
  325       v = e2(i) / u
  330       u = d(i) - x1 - v
            if (u .lt. 0.0d0) s = s + 1
  340    continue
c
         go to (60,80,200,220,360), isturm
c     .......... refine intervals ..........
  360    if (s .ge. k) go to 400
         xu = x1
         if (s .ge. m1) go to 380
         rv4(m1) = x1
         go to 300
  380    rv4(s+1) = x1
         if (rv5(s) .gt. x1) rv5(s) = x1
         go to 300
  400    x0 = x1
         go to 300
c     .......... k-th eigenvalue found ..........
  420    rv5(k) = x1
      k = k - 1
      if (k .ge. m1) go to 250
c     .......... find vectors by inverse iteration ..........
      norm = dabs(d(p))
      ip = p + 1
c
      do 500 i = ip, q
  500 norm = dmax1(norm, dabs(d(i)) + dabs(e(i)))
c     .......... eps2 is the criterion for grouping,
c                eps3 replaces zero pivots and equal
c                roots are modified by eps3,
c                eps4 is taken very small to avoid overflow ..........
      eps2 = 1.0d-3 * norm
      eps3 = epslon(norm)
      uk = q - p + 1
      eps4 = uk * eps3
      uk = eps4 / dsqrt(uk)
      group = 0
      s = p
c
      do 920 k = m1, m2
         r = r + 1
         its = 1
         w(r) = rv5(k)
         x1 = rv5(k)
c     .......... look for close or coincident roots ..........
         if (k .eq. m1) go to 520
         if (x1 - x0 .ge. eps2) group = -1
         group = group + 1
         if (x1 .le. x0) x1 = x0 + eps3
c     .......... elimination with interchanges and
c                initialization of vector ..........
  520    v = 0.0d0
c
         do 580 i = p, q
            rv6(i) = uk
            if (i .eq. p) go to 560
            if (dabs(e(i)) .lt. dabs(u)) go to 540
            xu = u / e(i)
            rv4(i) = xu
            rv1(i-1) = e(i)
            rv2(i-1) = d(i) - x1
            rv3(i-1) = 0.0d0
            if (i .ne. q) rv3(i-1) = e(i+1)
            u = v - xu * rv2(i-1)
            v = -xu * rv3(i-1)
            go to 580
  540       xu = e(i) / u
            rv4(i) = xu
            rv1(i-1) = u
            rv2(i-1) = v
            rv3(i-1) = 0.0d0
  560       u = d(i) - x1 - xu * v
            if (i .ne. q) v = e(i+1)
  580    continue
c
         if (u .eq. 0.0d0) u = eps3
         rv1(q) = u
         rv2(q) = 0.0d0
         rv3(q) = 0.0d0
c     .......... back substitution
c                for i=q step -1 until p do -- ..........
  600    do 620 ii = p, q
            i = p + q - ii
            rv6(i) = (rv6(i) - u * rv2(i) - v * rv3(i)) / rv1(i)
            v = u
            u = rv6(i)
  620    continue
c     .......... orthogonalize with respect to previous
c                members of group ..........
         if (group .eq. 0) go to 700
c
         do 680 jj = 1, group
            j = r - group - 1 + jj
            xu = 0.0d0
c
            do 640 i = p, q
  640       xu = xu + rv6(i) * z(i,j)
c
            do 660 i = p, q
  660       rv6(i) = rv6(i) - xu * z(i,j)
c
  680    continue
c
  700    norm = 0.0d0
c
         do 720 i = p, q
  720    norm = norm + dabs(rv6(i))
c
         if (norm .ge. 1.0d0) go to 840
c     .......... forward substitution ..........
         if (its .eq. 5) go to 960
         if (norm .ne. 0.0d0) go to 740
         rv6(s) = eps4
         s = s + 1
         if (s .gt. q) s = p
         go to 780
  740    xu = eps4 / norm
c
         do 760 i = p, q
  760    rv6(i) = rv6(i) * xu
c     .......... elimination operations on next vector
c                iterate ..........
  780    do 820 i = ip, q
            u = rv6(i)
c     .......... if rv1(i-1) .eq. e(i), a row interchange
c                was performed earlier in the
c                triangularization process ..........
            if (rv1(i-1) .ne. e(i)) go to 800
            u = rv6(i-1)
            rv6(i-1) = rv6(i)
  800       rv6(i) = u - rv4(i) * rv6(i-1)
  820    continue
c
         its = its + 1
         go to 600
c     .......... normalize so that sum of squares is
c                1 and expand to full order ..........
  840    u = 0.0d0
c
         do 860 i = p, q
  860    u = pythag(u,rv6(i))
c
         xu = 1.0d0 / u
c
         do 880 i = 1, n
  880    z(i,r) = 0.0d0
c
         do 900 i = p, q
  900    z(i,r) = rv6(i) * xu
c
         x0 = x1
  920 continue
c
  940 if (q .lt. n) go to 100
      go to 1001
c     .......... set error -- non-converged eigenvector ..........
  960 ierr = 4 * n + r
      go to 1001
c     .......... set error -- underestimate of number of
c                eigenvalues in interval ..........
  980 ierr = 3 * n + 1
 1001 lb = t1
      ub = t2
      return
      end
