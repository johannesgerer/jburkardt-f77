      subroutine dgbbqr2 ( m, n, ml, mu, a, lda, tau, work, info )

c*********************************************************************72
c
cc DGBBQR2 QR factors an M by N band matrix in GB format, without blocking.
c
c  Discussion:
c
c    DGBBQR2 computes a QR factorization of a real M by N band matrix A
c    with lower band ML and upper band MU: A = Q * R.
c
c    A is stored as a packed band matrix.
c    Input matrix A has ML subdiagonals and MU superdiagonals.
c    Output matrix R has ML + MU superdiagonals.
c
c    Example of A with M = 8, N = 6, ML = 2, and MU = 1.
c    Left, input matrix; Right, output matrix.
c
c      x  x  0  0  0  0       r   r   r   r   0   0
c      x  x  x  0  0  0       v1  r   r   r   r   0
c      x  x  x  x  0  0       v1  v2  r   r   r   r
c      0  x  x  x  x  0       0   v2  v3  r   r   r
c      0  0  x  x  x  x       0   0   v3  v4  r   r
c      0  0  0  x  x  x       0   0   0   v4  v5  r
c      0  0  0  0  x  x       0   0   0   0   v5  v6
c      0  0  0  0  0  x       0   0   0   0   0   v6
c
c  Licensing:
c
c    Copyright (c) 2007, Universidad Jaume I de Castellon
c    All rights reserved.
c
c    Redistribution and use in source and binary forms, with or without
c    modification, are permitted provided that the following conditions 
c    are met:
c      * Redistributions of source code must retain the above copyright
c        notice, this list of conditions and the following disclaimer.
c      * Redistributions in binary form must reproduce the above copyright
c        notice, this list of conditions and the following disclaimer in the
c        documentation and/or other materials provided with the distribution.
c      * Neither the name of the <organization> nor the
c        names of its contributors may be used to endorse or promote 
c        products derived from this software without specific prior written 
c        permission.
c
c    This software is provided by <copyright holder> ''as is'' and any
c    express or implied warranties, including, but not limited to, the 
c    implied warranties of merchantability and fitness for a particular 
c    purpose are disclaimed.  In no event shall <copyright holder> be liable 
c    for any direct, indirect, incidental, special, exemplary, or 
c    consequential damages (including, but not limited to, procurement of 
c    substitute goods or services; loss of use, data, or profits; or 
c    business interruption) however caused and on any theory of liability, 
c    whether in contract, strict liability, or tort (including negligence 
c    or otherwise) arising in any way out of the use of this software, even 
c    if advised of the possibility of such damage.
c
c  Modified:
c
c    04 April 2010
c
c  Author:
c
c    Alfredo Remon, Enrique Quintana-Orti, Gregorio Quintana-Orti,
c    Dept. de Ingenieria y Ciencia de Computadores, Univ. Jaume I
c    12.080 Castellon, Spain
c    {gquintan,remon,quintana}@icc.uji.es
c
c  Reference:
c
c    Alfredo Remon, Enrique Quintana-Orti, Gregorio Quintana-Orti,
c    LAPACK-Style Codes for the QR Factorization of Banded Matrices,
c    To Appear.
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c    0 <= M.
c
c    Input, integer N, the number of columns of the matrix.
c    0 <= N.
c
c    Input, integer ML, the number of nonzero subdiagonals.
c    0 <= ML.
c
c    Input, integer MU, the number of nonzero superdiagonals.
c    0 <= MU.
c
c    Input/output, double precision A(LDA,N).
c    A is stored as a packed matrix: dense A is M by N, whereas
c    band A is (2*ML+MU+1) by N.
c    On entry, the M by N matrix AB in band storage in rows ML+1 to
c    2*ML+MU+1; rows 1 to ML of the array need not be set.
c    This matrix has ML subdiagonals and MU superdiagonals with data.
c    On exit, the elements on and above the diagonal of the array
c    contain the upper band matrix R.
c    The output matrix R has ML+MU superdiagonals.
c    The elements within the band below the diagonal, with the array
c    TAU, represent the orthogonal matrix Q as a product of
c    elementary reflectors.
c
c    Input, integer LDA, the leading dimension of A.
c    2*ML+MU+1 <= LDA.
c
c    Output, double precision TAU(min(M,N)), the scalar factors of the 
c    elementary reflectors.
c
c    Workspace, double precision WORK(min(N,MU+ML)).
c
c    Output, integer INFO, error flag.
c    0: successful exit
c    nonzero, if INFO = -I, the I-th argument had an illegal value.
c
      implicit none

      integer lda
      integer n

      double precision a(lda,n)
      double precision diag
      integer info
      integer j
      integer m
      integer mh
      integer ml
      integer mn
      integer mu
      integer nh
      double precision one
      parameter ( one = 1.0D+00 )
      double precision tau(*)
      double precision work(*)
      double precision zero
      parameter ( zero = 0.0D+00 )
c
c  Test the input arguments.
c
      info = 0

      if ( m .lt. 0 ) then
        info = - 1
      else if ( n .lt. 0 ) then
        info = - 2
      else if ( ml .lt. 0 ) then
        info = - 3
      else if ( mu .lt. 0 ) then
        info = - 4
      else if ( lda .lt. 2 * ml + mu + 1 ) then
        info = - 6
      end if

      if ( info .ne. 0 ) then
        call xerbla ( 'dgbbqr2', -info )
        return
      end if
c
c  Quick return if possible.
c
      mn = min ( m, n )
      if ( mn .eq. 0 ) then
        return
      end if
c
c  Set rows 1:ML to zero.
c
      call dlaset ( 'all', ml, n, zero, zero, a, lda )

      do j = 1, mn

        mh = min ( ml + 1, m - j + 1 )
c
c  Generate reflector H(j) to annihilate A(j+1:j+mh-1,j).
c
        call dlarfg ( mh, a(ml+mu+1,j), a(ml+mu+1+min(1,ml),j), 1, 
     &    tau(j) )

        nh = min ( n - j, mu + ml )
c
c  Apply reflector H(j) to rest of matrix: A(j:j+mh-1,j+1:j+nh).
c
        if ( 0 .lt. nh ) then

          diag = a(ml+mu+1,j)

          a(ml+mu+1,j) = one

          call dlarf ( 'left', mh, nh, a(ml+mu+1,j), 1, tau(j),
     &      a(ml+mu,j+1), lda - 1, work )

          a(ml+mu+1,j) = diag

        end if

      end do

      return
      end
      subroutine dgbbqrf ( nb, m, n, ml, mu, a, lda, tau, work, info )

c*********************************************************************72
c
cc DGBBQRF QR factors an M by N band matrix in GB format, using blocking.
c
c  Discussion:
c
c    DGBBQRF computes a QR factorization of a real M by N band matrix A
c    with lower band ML and upper band MU: A = Q * R.
c
c    A is stored as a packed matrix: dense A is M by N, whereas
c    band A is (2*ML+MU+NB) by N.
c    Input matrix AB has ML subdiagonals and MU superdiagonals.
c    Output matrix R has ML+MU+NB-1 superdiagonals.
c
c    Example of A with M = 8, N = 6, ML = 2, MU = 1, and NB = 2.
c    Left, input matrix; Right, output matrix.
c
c      x  x  0  0  0  0       r   r   r   r   r   0
c      x  x  x  0  0  0       v1  r   r   r   r   r
c      x  x  x  x  0  0       v1  v2  r   r   r   r
c      0  x  x  x  x  0       0   v2  v3  r   r   r
c      0  0  x  x  x  x       0   0   v3  v4  r   r
c      0  0  0  x  x  x       0   0   0   v4  v5  r
c      0  0  0  0  x  x       0   0   0   0   v5  v6
c      0  0  0  0  0  x       0   0   0   0   0   v6
c
c  Licensing:
c
c    Copyright (c) 2007, Universidad Jaume I de Castellon
c    All rights reserved.
c
c    Redistribution and use in source and binary forms, with or without
c    modification, are permitted provided that the following conditions 
c    are met:
c      * Redistributions of source code must retain the above copyright
c        notice, this list of conditions and the following disclaimer.
c      * Redistributions in binary form must reproduce the above copyright
c        notice, this list of conditions and the following disclaimer in the
c        documentation and/or other materials provided with the distribution.
c      * Neither the name of the <organization> nor the
c        names of its contributors may be used to endorse or promote 
c        products derived from this software without specific prior written 
c        permission.
c
c    This software is provided by <copyright holder> ''as is'' and any
c    express or implied warranties, including, but not limited to, the 
c    implied warranties of merchantability and fitness for a particular 
c    purpose are disclaimed.  In no event shall <copyright holder> be liable 
c    for any direct, indirect, incidental, special, exemplary, or 
c    consequential damages (including, but not limited to, procurement of 
c    substitute goods or services; loss of use, data, or profits; or 
c    business interruption) however caused and on any theory of liability, 
c    whether in contract, strict liability, or tort (including negligence 
c    or otherwise) arising in any way out of the use of this software, even 
c    if advised of the possibility of such damage.
c
c  Modified:
c
c    03 April 2010
c
c  Author:
c
c    Alfredo Remon, Enrique Quintana-Orti, Gregorio Quintana-Orti,
c    Dept. de Ingenieria y Ciencia de Computadores, Univ. Jaume I
c    12.080 Castellon, Spain
c    {gquintan,remon,quintana}@icc.uji.es
c
c  Reference:
c
c    Alfredo Remon, Enrique Quintana-Orti, Gregorio Quintana-Orti,
c    LAPACK-Style Codes for the QR Factorization of Banded Matrices,
c    To Appear.
c
c  Parameters:
c
c    Input, integer NB, the block size to use.
c    1 <= NB.
c
c    Input, integer M, the number of rows in the matrix.
c    0 <= M.
c
c    Input, integer N, the number of columns in the matrix.
c    0 <= N.
c
c    Input, integer ML, the number of nonzero subdiagonals.
c    0 <= ML.
c
c    Input, integer MU, the number of nonzero superdiagonals.
c    0 <= MU.
c
c    Input/output, double precision A(LDA,N).
c    A is stored as a packed matrix: dense A is M by N, whereas
c    band A is (2*ML+MU+NB) by N.
c    On entry, the M by N matrix A in band storage in rows ML+NB to
c    2*ML+MU+NB; rows 1 to ML+NB-1 of the array need not be set.
c    This matrix has ML subdiagonals and MU superdiagonals with data.
c    On exit, the elements on and above the diagonal of the array
c    contain the upper band matrix R.
c    The output matrix R has ML+MU+NB-1 superdiagonals.
c    The elements within the band below the diagonal, with the array
c    TAU, represent the orthogonal matrix Q as a product of
c    elementary reflectors.
c
c    Input, integer LDA, the leading dimension of A.
c    max ( 1, M ) <= LDA.
c
c    Output, double precision TAU(min(M,N)), the scalar factors of the 
c    elementary reflectors.
c
c    Workspace, double precision WORK(DIMWORK).
c    NB*NB + MIN(N,ML+MU)*NB + MIN(M,ML+NB)*NB <= DIMWORK.
c
c    Output, integer INFO, error flag.
c    0, no errors.
c    nonzero, if INFO = -I, the I-th argument had an illegal value.
c
      implicit none

      integer lda

      double precision a(lda,*)
      integer am
      integer aml
      integer amu
      integer an
      integer anb
      integer ii
      integer info
      integer irwk
      integer it
      integer iv
      integer j
      integer jb
      integer jj
      integer jlc
      integer m
      integer ml
      integer mn
      integer mu
      integer n
      integer nb
      integer ncol
      double precision tau(*)
      double precision work(*)
      double precision zero
      parameter ( zero = 0.0D+00 )
c
c  Test the input arguments.
c
      info = 0

      if ( nb .lt. 1 ) then
        info = -1
      else if ( m .lt. 0 ) then
        info = -2
      else if ( n .lt. 0 ) then
        info = -3
      else if ( ml .lt. 0 ) then
        info = -4
      else if ( mu .lt. 0 ) then
        info = -5
      else if ( lda .lt. max ( 1, 2 * ml + mu + nb ) ) then
        info = -7
      end if

      if ( info .ne. 0 ) then
        call xerbla ( 'dgbbqrf', -info )
        return
      end if
c
c  Quick return if possible.
c
      if ( m .eq. 0 .or. n .eq. 0 .or. ml .eq. 0 ) then
        return
      end if
c
c  Adjust matrix sizes: AM, AN, AML, AMU, and ANB.
c
      am = min ( m, n + ml )
      an = min ( n, m + mu )
      aml = min ( m - 1, ml )
      amu = min ( n - 1, mu )
      anb = nb
      if ( aml .lt. anb ) then
        anb = aml
      end if

      mn = min ( am, an )

      it = 1
      iv = it + anb * anb
      irwk = iv + min ( aml + anb, am ) * anb

      ncol = min ( an, am - aml )
c
c  Factorization of full band matrix A(:,1:ncol).
c
      do j = 1, ncol, anb

        jb = min ( anb, ncol - j + 1 )
c
c  Factorize block A(j:j+aml+jb-1,j:j+jb-1).
c
        call dgebqr2 ( aml + jb, jb, aml, amu, a(ml+mu+nb,j), lda - 1, 
     &    tau(j), work(irwk), info )

        if ( info .ne. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'DGBBQRF - Fatal error!'
          write ( *, '(a)' ) '  Error return from DGEBQR2.'
          write ( *, '(a,i8)' ) '  INFO = ', info
          stop
        end if

        if ( j + jb .le. an ) then
c
c  WORK(IV) := Y (the lower part of Y is padded with zeros).
c
          do jj = 1, jb

            do ii = 1, aml + jj
              work(iv-1+(aml+jb)*(jj-1)+ii) =
     &          a(ml+mu+nb+ii-jj,j+jj-1)
            end do

            do ii = aml + jj + 1, aml + jb
              work(iv-1+(aml+jb)*(jj-1)+ii) = zero
            end do

          end do
c
c  Form the triangular factor T of the block reflector.
c
          call dlarft ( 'forward', 'columnwise', aml + jb, jb,
     &      work(iv), aml + jb, tau(j), work(it), jb )
c
c  Apply block reflector to A(j:j+aml+jb-1,j+jb:an) from the left.
c
          jlc = min ( an, j + jb - 1 + aml + amu )

          call dlarfb ( 'left', 'transpose', 'forward',
     &      'columnwise', aml + jb, jlc - j - jb + 1, jb,
     &      work(iv), aml + jb, work(it), jb,
     &      a(ml+mu+nb-jb,j+jb), lda - 1,
     &      work(irwk), jlc - j - jb + 1 )
        end if

      end do
c
c  Factorization of rectangular matrix A(:,ncol+1:mn).
c
      do j = ncol + 1, mn, anb

        jb = min ( anb, mn - j + 1 )
c
c  Factorize block A(j:am,j:j+jb-1).
c
        call dgeqr2 ( am - j + 1, jb, a(ml+mu+nb,j), lda - 1, tau(j),
     &    work(irwk), info )

        if ( info .ne. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'DGBBQRF - Fatal error!'
          write ( *, '(a)' ) '  Error return from DGEQR2.'
          write ( *, '(a,i8)' ) '  INFO = ', info
          stop
        end if

        if ( j + jb .le. an ) then
c
c  Form the triangular factor T of the block reflector.
c
          call dlarft ( 'forward', 'columnwise', am - j + 1, jb,
     &      a(ml+mu+nb,j), lda - 1, tau(j), work(it), jb )
c
c  Apply block reflector to A(j:am,j+jb:an) from the left.
c
          call dlarfb ( 'left', 'transpose', 'forward',
     &      'columnwise', am - j + 1, an - j - jb + 1, jb,
     &      a(ml+mu+nb,j), lda - 1, work(it), jb,
     &      a(ml+mu+nb-jb,j+jb), lda - 1,
     &      work(irwk), an - j - jb + 1 )

        end if

      end do

      return
      end
      subroutine dgbbqrs ( m, n, ml, mu, nrhs, a, lda, tau, b, ldb, 
     &  work, lwork, info )

c*********************************************************************72
c
cc DGBBQRS solves A*X = B when A has been factored by DGBBQR2.
c
c  Discussion:
c
c    Solve the least squares problem
c      min || A*X - B ||
c    using the QR factorization
c      A = Q*R
c    computed by DGBBQR2.
c
c    Here, the matrix A is an M by N band matrix, which was stored in the
c    standard LINPACK/LAPACK "GB" format.
c
c    This matrix was QR-factored by DGBBQR2, which is able to take advantage
c    of the GB format, and the results of the factorization overwrote
c    A, with some extra information stored in TAU.
c
c    Now one or more linear systems A * x = b are to solved using
c    the QR procedure.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 April 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix A.
c    0 <= M.
c
c    Input, integer N, the number of columns of the matrix A.
c    0 <= N <= M.
c
c    Input, integer ML, the number of nonzero subdiagonals.
c    0 <= ML.
c
c    Input, integer MU, the number of nonzero superdiagonals.
c    0 <= MU.
c
c    Input, integer NRHS, the number of columns of B.  
c    0 <= NRHS.
c
c    Input, double precision A(LDA,N), part of the QR factorization
c    computed by DGBBQRF.
c
c    Input, integer LDA, the leading dimension of A.
c    2*ML+MU+1 <= LDA.
c
c    Input, double precision TAU(N), the scalar factor of the elementary
c    reflectors H, as computed by DGBBQRF.
c
c    Input/output, double precision B(LDB,NRHS).
c    On entry, the M by NRHS right hand side matrix B.
c    On exit, the N by NRHS solution matrix X.
c
c    Input, integer LDA, the leading dimension of B.
c    M <= LDB.
c
c    Workspace, double precision WORK(LWORK).
c
c    Input, integer LWORK, the length of the array WORK.  LWORK must be at 
c    least NRHS, and should be at least NRHS*NB, where NB is the block size
c    for this environment.
c
c    Output, integer INFO, error flag.
c    0, no errors.
c    nonzero, if INFO = -I, the I-th argument had an illegal value.
c
      implicit none

      integer lda
      integer ldb
      integer lwork
      integer n
      integer nrhs

      double precision a(lda,n)
      double precision aii
      double precision b(ldb,nrhs)
      integer i
      integer i2
      integer info
      integer j
      integer k
      integer ml
      integer mu
      integer m
      double precision tau(n)
      double precision temp
      double precision work(lwork)
c
c  B := Q' * B.
c
c  Q is a real orthogonal matrix of order M, 
c  the product of N elementary reflectors.
c
c  B is an M by NRHS vector.
c 
      do i = 1, n
c
c  H(I) is applied to B(I:M,1:NRHS).
c
        if ( tau(i) .ne. 0.0D+00 ) then

          aii = a(i-i+ml+mu+1,i)
          a(i-i+ml+mu+1,i) = 1.0D+00
c
c  work := B' * v.
c
          do j = 1, nrhs
            temp = 0.0D+00
            do i2 = i, min ( m, i + ml )
              temp = temp + b(i2,j) * a(i2-i+ml+mu+1,i)
            end do
            work(j) = temp
          end do
c
c  B := B - tau * v * work'.
c
          do j = 1, nrhs
            if ( work(j) .ne. 0.0D+00 )then
              temp = work(j)
              do i2 = i, min ( m, i + ml )
                b(i2,j) = b(i2,j) - tau(i) * a(i2-i+ml+mu+1,i) * temp
              end do
            end if
          end do

          a(i-i+ml+mu+1,i) = aii

        end if

      end do
c
c  Solve R*X = B(1:n,:).
c
      do j = 1, nrhs
        do k = n, 1, -1
          if ( b(k,j) .ne. 0.0D+00 )then
            b(k,j) = b(k,j) / a(k-k+ml+mu+1,k)
            do i = max ( 1, k - ml - mu ), k - 1
              b(i,j) = b(i,j) - b(k,j) * a(i-k+ml+mu+1,k) 
           end do
          end if
        end do
      end do

      return
      end
      subroutine dgebqr2 ( m, n, ml, mu, a, lda, tau, work, info )

c*********************************************************************72
c
cc DGEBQR2 QR factors an M by N band matrix in GE format, with no blocking.
c
c  Discussion:
c
c    DGEBQR2 computes a QR factorization of a real M by N band matrix A
c    with lower band ML and upper band MU: A = Q * R.
c
c    A is stored as a general matrix.
c    Input matrix A has ML subdiagonals and MU superdiagonals.
c    Output matrix R has ML+MU superdiagonals.
c
c    Example of A with M = 8, N = 6, ML = 2, and MU = 1.
c    Left, input matrix; Right, output matrix.
c
c      x  x  0  0  0  0       r   r   r   r   0   0
c      x  x  x  0  0  0       v1  r   r   r   r   0
c      x  x  x  x  0  0       v1  v2  r   r   r   r
c      0  x  x  x  x  0       0   v2  v3  r   r   r
c      0  0  x  x  x  x       0   0   v3  v4  r   r
c      0  0  0  x  x  x       0   0   0   v4  v5  r
c      0  0  0  0  x  x       0   0   0   0   v5  v6
c      0  0  0  0  0  x       0   0   0   0   0   v6
c
c    The matrix Q is represented as a product of elementary reflectors
c
c      Q = H(1) H(2) . . . H(k), where k = min(m,n).
c
c    Each H(i) has the form
c
c      H(i) = I - tau * v * v'
c
c    where tau is a real scalar, and v is a real vector with
c    v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
c    and tau in TAU(i).
c
c  Licensing:
c
c    Copyright (c) 2007, Universidad Jaume I de Castellon
c    All rights reserved.
c
c    Redistribution and use in source and binary forms, with or without
c    modification, are permitted provided that the following conditions 
c    are met:
c      * Redistributions of source code must retain the above copyright
c        notice, this list of conditions and the following disclaimer.
c      * Redistributions in binary form must reproduce the above copyright
c        notice, this list of conditions and the following disclaimer in the
c        documentation and/or other materials provided with the distribution.
c      * Neither the name of the <organization> nor the
c        names of its contributors may be used to endorse or promote 
c        products derived from this software without specific prior written 
c        permission.
c
c    This software is provided by <copyright holder> ''as is'' and any
c    express or implied warranties, including, but not limited to, the 
c    implied warranties of merchantability and fitness for a particular 
c    purpose are disclaimed.  In no event shall <copyright holder> be liable 
c    for any direct, indirect, incidental, special, exemplary, or 
c    consequential damages (including, but not limited to, procurement of 
c    substitute goods or services; loss of use, data, or profits; or 
c    business interruption) however caused and on any theory of liability, 
c    whether in contract, strict liability, or tort (including negligence 
c    or otherwise) arising in any way out of the use of this software, even 
c    if advised of the possibility of such damage.
c
c  Modified:
c
c    04 April 2010
c
c  Author:
c
c    Alfredo Remon, Enrique Quintana-Orti, Gregorio Quintana-Orti,
c    Dept. de Ingenieria y Ciencia de Computadores, Univ. Jaume I
c    12.080 Castellon, Spain
c    {gquintan,remon,quintana}@icc.uji.es
c
c  Reference:
c
c    Alfredo Remon, Enrique Quintana-Orti, Gregorio Quintana-Orti,
c    LAPACK-Style Codes for the QR Factorization of Banded Matrices,
c    To Appear.
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c    0 <= M.
c
c    Input, integer N, the number of columns of the matrix.
c    0 <= N.
c
c    Input, integer ML, the number of nonzero subdiagonals.
c    0 <= ML.
c
c    Input, integer MU, the number of nonzero superdiagonals.
c    0 <= MU.
c
c    Input/output, double precision(LDA,N).
c    On entry, the M by N matrix A. It has ML subdiagonals and MU
c    superdiagonals.
c    On exit, the elements on and above the diagonal of the array
c    contain the min(M,N) by N upper band matrix R.
c    The output matrix R has ML+MU superdiagonals.
c    The elements within the band below the diagonal, with the array
c    TAU, represent the orthogonal matrix Q as a product of
c    elementary reflectors.
c
c    Input, integer LDA, the leading dimension of the array A.  
c    max ( 1, M ) <= LDA.
c
c    Output, double precision TAU(min(M,N)),
c    the scalar factors of the elementary reflectors.
c
c    Workspace, double precision WORK(min(N,MU+ML)).
c
c    Output, integer INFO, error flag.
c    0, no errors.
c    nonzero, if INFO = -I, the I-th argument had an illegal value.
c
      implicit none

      integer lda
      integer n

      double precision a(lda,n)
      double precision diag
      integer info
      integer j
      integer m
      integer mh
      integer ml
      integer mn
      integer mu
      integer nh
      double precision one
      parameter ( one = 1.0D+00 )
      double precision tau(*)
      double precision work(*)
c
c  Test the input arguments.
c
      info = 0

      if ( m .lt. 0 ) then
        info = -1
      else if ( n .lt. 0 ) then
        info = -2
      else if ( ml .lt. 0 ) then
        info = -3
      else if ( mu .lt. 0 ) then
        info = -4
      else if ( lda .lt. max ( 1, m ) ) then
        info = -6
      end if

      if ( info .ne. 0 ) then
        call xerbla ( 'dgebqr2', -info )
        return
      end if
c
c  Quick return if possible.
c
      mn = min ( m, n )
      if ( mn .eq. 0 ) then
        return
      end if

      do j = 1, mn

        mh = min ( ml + 1, m - j + 1 )
c
c  Generate reflector H(j) to annihilate A(j+1:j+mh-1,j).
c
        call dlarfg ( mh, a(j,j), a(min(j+1,m),j), 1, tau(j) )

        nh = min ( n - j, mu + ml )

        if ( 0 .lt. nh ) then
c
c  Apply reflector H(j) to rest of matrix: A(j:j+mh-1,j+1:j+nh).
c
          diag = a(j,j)

          a(j,j) = one

          call dlarf ( 'left', mh, nh, a(j,j), 1, tau(j), a(j,j+1), 
     &      lda, work )

          a(j,j) = diag

        end if

      end do

      return
      end
      subroutine dgebqrf ( nb, m, n, ml, mu, a, lda, tau, work, info )

c*********************************************************************72
c
cc DGEBQRF QR factors an M by N band matrix stored in GE format, with blocking.
c
c  Discussion:
c
c    DGEBQRF computes a QR factorization of a real M by N band matrix A
c    with lower band ML and upper band MU: A = Q * R.
c
c    A is stored as a general matrix.
c    Input matrix A has ML subdiagonals and MU superdiagonals.
c    Output matrix R has ML+MU+NB-1 superdiagonals.
c
c    Example of A with M = 8, N = 6, ML = 2, MU = 1, and NB = 2.
c    Left, input matrix; Right, output matrix.
c
c      x  x  0  0  0  0       r   r   r   r   r   0
c      x  x  x  0  0  0       v1  r   r   r   r   r
c      x  x  x  x  0  0       v1  v2  r   r   r   r
c      0  x  x  x  x  0       0   v2  v3  r   r   r
c      0  0  x  x  x  x       0   0   v3  v4  r   r
c      0  0  0  x  x  x       0   0   0   v4  v5  r
c      0  0  0  0  x  x       0   0   0   0   v5  v6
c      0  0  0  0  0  x       0   0   0   0   0   v6
c
c  Licensing:
c
c    Copyright (c) 2007, Universidad Jaume I de Castellon
c    All rights reserved.
c
c    Redistribution and use in source and binary forms, with or without
c    modification, are permitted provided that the following conditions 
c    are met:
c      * Redistributions of source code must retain the above copyright
c        notice, this list of conditions and the following disclaimer.
c      * Redistributions in binary form must reproduce the above copyright
c        notice, this list of conditions and the following disclaimer in the
c        documentation and/or other materials provided with the distribution.
c      * Neither the name of the <organization> nor the
c        names of its contributors may be used to endorse or promote 
c        products derived from this software without specific prior written 
c        permission.
c
c    This software is provided by <copyright holder> ''as is'' and any
c    express or implied warranties, including, but not limited to, the 
c    implied warranties of merchantability and fitness for a particular 
c    purpose are disclaimed.  In no event shall <copyright holder> be liable 
c    for any direct, indirect, incidental, special, exemplary, or 
c    consequential damages (including, but not limited to, procurement of 
c    substitute goods or services; loss of use, data, or profits; or 
c    business interruption) however caused and on any theory of liability, 
c    whether in contract, strict liability, or tort (including negligence 
c    or otherwise) arising in any way out of the use of this software, even 
c    if advised of the possibility of such damage.
c
c  Modified:
c
c    04 April 2010
c
c  Author:
c
c    Alfredo Remon, Enrique Quintana-Orti, Gregorio Quintana-Orti,
c    Dept. de Ingenieria y Ciencia de Computadores, Univ. Jaume I
c    12.080 Castellon, Spain
c    {gquintan,remon,quintana}@icc.uji.es
c
c  Reference:
c
c    Alfredo Remon, Enrique Quintana-Orti, Gregorio Quintana-Orti,
c    LAPACK-Style Codes for the QR Factorization of Banded Matrices,
c    To Appear.
c
c  Parameters:
c
c    Input, integer NB, the block size.
c    1 <= NB.
c
c    Input, integer M, the number of rows of the matrix.
c    0 <= M.
c
c    Input, integer N, the number of columns of the matrix.
c    0 <= N.
c
c    Input, integer ML, the number of nonzero subdiagonals.
c    0 <= ML.
c
c    Input, integer MU, the number of nonzero superdiagonals.
c    0 <= MU.
c
c    Input/output, double precision A(LDA,N).
c    On entry, the M by N matrix A. It has ML subdiagonals and MU
c    superdiagonals.
c    On exit, the elements on and above the diagonal of the array
c    contain the min(M,N) by N upper band matrix R.
c    The output matrix R has ML+MU+NB-1 superdiagonals.
c    The elements within the band below the diagonal, with the array
c    TAU, represent the orthogonal matrix Q as a product of
c    elementary reflectors.
c
c    Input, integer LDA, the leading dimension of the array A.  
c    max ( 1, M ) <= LDA.
c
c    Output, double precision TAU(min(M,N)), the scalar factors of the 
c    elementary reflectors.
c
c    Workspace, double precision WORK(DIMWORK).
c    NB*NB + min(N,ML+MU)*NB + min(M,ML+NB)*NB <= DIMWORK.
c
c    Output, integer INFO, error flag.
c    0, no errors.
c    nonzero, if INFO = -I, the I-th argument had an illegal value.
c
      implicit none

      integer lda
      integer n

      double precision a(lda,n)
      integer am
      integer aml
      integer amu
      integer an
      integer anb
      integer ii
      integer info
      integer irwk
      integer it
      integer iv
      integer j
      integer jb
      integer jj
      integer jlc
      integer m
      integer ml
      integer mn
      integer mu
      integer nb
      integer ncol
      double precision tau(*)
      double precision work(*)
      double precision zero
      parameter ( zero = 0.0D+00 )
c
c  Test the input arguments.
c
      info = 0

      if ( nb .lt. 1 ) then
        info = -1
      else if ( m .lt. 0 ) then
        info = -2
      else if ( n .lt. 0 ) then
        info = -3
      else if ( ml .lt. 0 ) then
        info = -4
      else if ( mu .lt. 0 ) then
        info = -5
      else if ( lda .lt. max ( 1, m ) ) then
        info = -7
      end if

      if ( info .ne. 0 ) then
        call xerbla ( 'dgebqrf', -info )
        return
      end if
c
c  Quick return if possible.
c
      if ( m .eq. 0 .or. n .eq. 0 .or. ml .eq. 0 ) then
        return
      end if
c
c  Adjust matrix sizes: AM, AN, AML, AMU, and ANB.
c
      am = min ( m, n + ml )
      an = min ( n, m + mu )
      aml = min ( m - 1, ml )
      amu = min ( n - 1, mu )
      anb = nb
      if ( aml .lt. anb ) then
        anb = aml
      end if

      mn = min ( am, an )

      it = 1
      iv = it + anb * anb
      irwk = iv + min ( aml + anb, am ) * anb

      ncol = min ( an, am - aml )
c
c  Factorization of full band matrix A(:,1:ncol).
c
      do j = 1, ncol, anb

        jb = min ( anb, ncol - j + 1 )
c
c  Factorize block A(j:j+aml+jb-1,j:j+jb-1).
c
        call dgebqr2 ( aml + jb, jb, aml, amu, a(j,j), lda, tau(j),
     &    work(irwk), info )

        if ( info .ne. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'DGEBQRF - Fatal error!'
          write ( *, '(a)' ) '  Error return from DGEBQR2.'
          write ( *, '(a,i8)' ) '  INFO = ', info
          stop
        end if

        if ( j + jb .le. an ) then
c
c  WORK(IV) := Y (the lower part of Y is padded with zeros).
c
          do jj = 1, jb

            do ii = 1, aml + jj
              work(iv-1+(aml+jb)*(jj-1)+ii) =
     &          a(j+ii-1,j+jj-1)
            end do

            do ii = aml + jj + 1, aml + jb
              work(iv-1+(aml+jb)*(jj-1)+ii) = zero
            end do
          end do
c
c  Form the triangular factor T of the block reflector.
c
          call dlarft ( 'forward', 'columnwise', aml + jb, jb,
     &      work(iv), aml + jb, tau(j), work(it), jb )
c
c  Apply block reflector to A(j:j+aml+jb-1,j+jb:an) from the left.
c
          jlc = min ( an, j + jb - 1 + aml + amu )

          call dlarfb ( 'left', 'transpose', 'forward',
     &      'columnwise', aml + jb, jlc - j - jb + 1, jb,
     &      work(iv), aml + jb, work(it), jb,
     &      a(j,j+jb), lda, work(irwk), jlc - j - jb + 1 )

        end if

      end do
c
c  Factorization of rectangular matrix A(:,ncol+1:mn).
c
      do j = ncol + 1, mn, anb

        jb = min ( anb, mn - j + 1 )
c
c  Factorize block A(j:am,j:j+jb-1).
c
        call dgeqr2 ( am - j + 1, jb, a(j,j), lda, tau(j), work(irwk),
     &    info )

        if ( info .ne. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'DGEBQRF - Fatal error!'
          write ( *, '(a)' ) '  Error return from DGEQR2.'
          write ( *, '(a,i8)' ) '  INFO = ', info
          stop
        end if

        if ( j + jb .le. an ) then
c
c  Form the triangular factor T of the block reflector.
c
          call dlarft ( 'forward', 'columnwise', am - j + 1, jb,
     &      a(j,j), lda, tau(j), work(it), jb )
c
c  Apply block reflector to A(j:am,j+jb:an) from the left.
c
          call dlarfb ( 'left', 'transpose', 'forward',
     &      'columnwise', am - j + 1, an - j - jb + 1, jb,
     &      a(j,j), lda, work(it), jb,
     &      a(j,j+jb), lda, work(irwk), an - j - jb + 1 )

        end if

      end do

      return
      end
      subroutine dgeqrs ( m, n, nrhs, a, lda, tau, b, ldb, work, lwork,
     &  info )

c*********************************************************************72
c
cc DGEQRS solves a linear system factored by DGEQRF.
c
c  Discussion:
c
c    This routine is not part of the main LAPACK release.  It only has
c    the status of a "testing" routine.  We need direct access to this 
c    routine to make comparisons in source code and results.
c
c    This routine solves the least squares problem
c      min || A*X - B ||
c    using the QR factorization
c      A = Q*R
c    computed by DGEQRF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 March 2010
c
c  Author:
c
c    This version by John Burkardt.
c
c  Reference:
c
c    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
c    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
c    Sven Hammarling, Alan McKenney, Danny Sorensen,
c    LAPACK User's Guide,
c    Third Edition,
c    SIAM, 1999,
c    ISBN: 0898714478,
c    LC: QA76.73.F25L36
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix A.  
c    0 <= M.
c
c    Input, integer N, the number of columns of the matrix A.  
c    0 <= N <= M.
c
c    Input, integer NRHS, the number of columns of B.  
c    0 <= NRHS.
c
c    Input, double precision A(LDA,N)
c    Details of the QR factorization of the original matrix A as
c    returned by DGEQRF.
c
c    Input, integer LDA, the leading dimension of the array A.  
c    M <= LDA.
c
c    Input, double precision TAU(N)
c    Details of the orthogonal matrix Q.
c
c    Input/output, double precision B(LDB,NRHS),
c    On entry, the M by NRHS right hand side matrix B.
c    On exit, the N by NRHS solution matrix X.
c
c    Input, integer LDB, the leading dimension of the array B. 
c    M <= LDB.
c
c    Workspace, double precision WORK(LWORK).
c
c    Input, integer LWORK,
c    The length of the array WORK.  LWORK must be at least NRHS,
c    and should be at least NRHS*NB, where NB is the block size
c    for this environment.
c
c    Output, integer INFO,
c    = 0: successful exit
c    < 0: if INFO = -I, the I-th argument had an illegal value
c
      implicit none

      integer lda
      integer ldb
      integer lwork
      integer n
      integer nrhs

      double precision a(lda,n)
      double precision b(ldb,nrhs)
      integer info
      integer m
      double precision one
      parameter ( one = 1.0D+00 )
      double precision tau(n)
      double precision work(lwork)

      info = 0

      if ( m .lt. 0 ) then
        info = - 1
      else if ( n .lt. 0 .or. m .lt. n ) then
        info = - 2
      else if ( nrhs .lt. 0 ) then
        info = - 3
      else if ( lda .lt. max ( 1, m ) ) then
        info = - 5
      else if ( ldb .lt. max ( 1, m ) ) then
        info = - 8
      else if ( ( lwork .lt. 1 .or. lwork .lt. nrhs ) .and. 
     &  0 .lt. m .and. 0 .lt. n ) then
         info = - 10
      end if

      if ( info .ne. 0 ) then
         call xerbla ( 'DGEQRS', - info )
         return
      end if
c
c  Quick return if possible.
c
      if ( n .eq. 0 .or. nrhs .eq. 0 .or. m .eq. 0 ) then
        return
      end if
c
c  Compute B := Q' * B.
c
      call dormqr ( 'Left', 'Transpose', m, nrhs, n, a, lda, tau, b, 
     &  ldb, work, lwork, info )
c
c  Solve R * X = B(1:n,:).
c
      call dtrsm ( 'Left', 'Upper', 'No transpose', 'Non-unit', n, nrhs,
     &  one, a, lda, b, ldb )

      return
      end
      subroutine dgeqrs_two ( m, n, nrhs, a, lda, tau, b, ldb, work, 
     &  lwork, info )

c*********************************************************************72
c
cc DGEQRS_TWO solves A*X = B when A has been factored by DGEQRF or DGEBQR2.
c
c  Discussion:
c
c    This routine solve the least squares problem
c      min || A*X - B ||
c    using the QR factorization
c      A = Q*R
c    computed by DGEQRF.
c
c    This routine is a revised version of the LAPACK testing code
c    DGEQRS.  It replaces all the calls to LAPACK subroutines by
c    explicit code.
c
c    It is intended to be a guide for devising a corresponding
c    routine to handle band matrices which have been QR factored
c    by DGBBQR2 or by DGBBQRF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 March 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
c    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
c    Sven Hammarling, Alan McKenney, Danny Sorensen,
c    LAPACK User's Guide,
c    Third Edition,
c    SIAM, 1999,
c    ISBN: 0898714478,
c    LC: QA76.73.F25L36
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix A.  
c    0 <= M.
c
c    Input, integer N, the number of columns of the matrix A.
c    0 <= N <= M.
c
c    Input, integer NRHS, the number of columns of B.  
c    0 <= NRHS.
c
c    Input, double precision A(LDA,N).
c    Details of the QR factorization of the original matrix A as
c    returned by DGEQRF.
c
c    Input, integer LDA, the leading dimension of the array A.  
c    M <= LDA.
c
c    Input, double precision TAU(N), the scalar factor of the elementary
c    reflectors.
c
c    Input/output, double precision B(LDB,NRHS).
c    On entry, the M by NRHS right hand side matrix B.
c    On exit, the N by NRHS solution matrix X.
c
c    Input, integer LDB, the leading dimension of the array B. 
c    M <= LDB.
c
c    Workspace, double precision WORK(LWORK).
c
c    Input, integer LWORK.
c    The length of the array WORK.  LWORK must be at least NRHS,
c    and should be at least NRHS*NB, where NB is the block size
c    for this environment.
c
c    Output, integer INFO.
c    = 0: successful exit
c    < 0: if INFO = -i, the i-th argument had an illegal value
c
      implicit none

      integer lda
      integer ldb
      integer lwork
      integer n
      integer nrhs

      double precision a(lda,n)
      double precision aii
      double precision b(ldb,nrhs)
      integer i
      integer i2
      integer info
      integer j
      integer k
      integer m
      double precision tau(n)
      double precision temp
      double precision work(lwork)
c
c  B := Q' * B.
c
c  Q is an orthogonal matrix of order M, the product of N elementary reflectors.
c
c  B is an M by NRHS vector.
c 
      do i = 1, n
c
c  H(I) is applied to B(I:M,1:NRHS).
c
        if ( tau(i) .ne. 0.0D+00 ) then

          aii = a(i,i)
          a(i,i) = 1.0D+00
c
c  work := B' * v.
c
          do j = 1, nrhs
            temp = 0.0D+00
            do i2 = i, m
              temp = temp + b(i2,j) * a(i2,i)
            end do
            work(j) = temp
          end do
c
c  B := B - tau * v * work'.
c
          do j = 1, nrhs
            if ( work(j) .ne. 0.0D+00 )then
              temp = work(j)
              do i2 = i, m
                b(i2,j) = b(i2,j) - tau(i) * a(i2,i) * temp
              end do
            end if
          end do

          a(i,i) = aii

        end if

      end do
c
c  Solve R*X = B(1:n,:).
c
      do j = 1, nrhs
        do k = n, 1, -1
          if ( b(k,j) .ne. 0.0D+00 )then
            b(k,j) = b(k,j) / a(k,k)
            do i = 1, k - 1
              b(i,j) = b(i,j) - b(k,j) * a(i,k)
            end do
          end if
        end do
      end do

      return
      end
      subroutine r8gb_print ( m, n, ml, mu, a, title )

c*********************************************************************72
c
cc R8GB_PRINT prints an R8GB matrix.
c
c  Discussion:
c
c    The R8GB storage format is for an M by N banded matrix, with lower 
c    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
c    extra superdiagonals, which may be required to store nonzero entries 
c    generated during Gaussian elimination.
c
c    The original M by N matrix is "collapsed" downward, so that diagonals
c    become rows of the storage array, while columns are preserved.  The
c    collapsed array is logically 2*ML+MU+1 by N.  
c
c    R8GB storage is used by LINPACK and LAPACK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c    M must be positive.
c
c    Input, integer N, the number of columns of the matrix.
c    N must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than min(M,N)-1..
c
c    Input, double precision A(2*ML+MU+1,N), the R8GB matrix.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer ml
      integer mu
      integer n

      double precision a(2*ml+mu+1,n)
      integer m
      character ( len = * )  title

      call r8gb_print_some ( m, n, ml, mu, a, 1, 1, m, n, title )

      return
      end
      subroutine r8gb_print_some ( m, n, ml, mu, a, ilo, jlo, ihi, 
     &  jhi, title )

c*********************************************************************72
c
cc R8GB_PRINT_SOME prints some of an R8GB matrix.
c
c  Discussion:
c
c    The R8GB storage format is for an M by N banded matrix, with lower 
c    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
c    extra superdiagonals, which may be required to store nonzero entries 
c    generated during Gaussian elimination.
c
c    The original M by N matrix is "collapsed" downward, so that diagonals
c    become rows of the storage array, while columns are preserved.  The
c    collapsed array is logically 2*ML+MU+1 by N.  
c
c    R8GB storage is used by LINPACK and LAPACK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c    M must be positive.
c
c    Input, integer N, the number of columns of the matrix.
c    N must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than min(M,N)-1.
c
c    Input, double precision A(2*ML+MU+1,N), the R8GB matrix.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column to be printed.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer incx
      parameter ( incx = 5 )
      integer ml
      integer mu
      integer n

      double precision a(2*ml+mu+1,n)
      character ( len = 14 ) ctemp(incx)
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
      integer m
      character ( len = * )  title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
c
c  Print the columns of the matrix, in strips of 5.
c
      do j2lo = jlo, jhi, incx

        j2hi = j2lo + incx - 1
        j2hi = min ( j2hi, n )
        j2hi = min ( j2hi, jhi )

        inc = j2hi + 1 - j2lo

        write ( *, '(a)' ) ' '

        do j = j2lo, j2hi
          j2 = j + 1 - j2lo
          write ( ctemp(j2), '(i7,7x)' ) j
        end do

        write ( *, '(a,5a14)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) '  ---'
c
c  Determine the range of the rows in this strip.
c
        i2lo = max ( ilo, 1 )
        i2lo = max ( i2lo, j2lo - mu - ml )
        i2hi = min ( ihi, m )
        i2hi = min ( i2hi, j2hi + ml )

        do i = i2lo, i2hi
c
c  Print out (up to) 5 entries in row I, that lie in the current strip.
c
          do j2 = 1, inc

            j = j2lo - 1 + j2

            if ( i .lt. j - ml - mu  .or. j + ml .lt. i ) then
              ctemp(j2) = '              '
            else
              write ( ctemp(j2), '(g14.6)' ) a(i-j+ml+mu+1,j)
            end if

          end do

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

        end do

      end do

      return
      end
      subroutine r8ge_print ( m, n, a, title )

c*********************************************************************72
c
cc R8GE_PRINT prints an R8 GE matrix
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
c    20 May 2004
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
c    Input, double precision A(M,N), the matrix.
c
c    Input, character ( len = * ) TITLE, a title to be printed.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      character ( len = * ) title

      call r8ge_print_some ( m, n, a, 1, 1, m, n, title )

      return
      end
      subroutine r8ge_print_some ( m, n, a, ilo, jlo, ihi, jhi,
     &  title )

c*********************************************************************72
c
cc R8GE_PRINT_SOME prints some of an R8 GE matrix.
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
c    Input, double precision A(M,N), an M by N matrix to be printed.
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

      double precision a(m,n)
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

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

        end do

      end do

      return
      end
      subroutine r8vec_print ( n, a, title )

c*********************************************************************72
c
cc R8VEC_PRINT prints an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
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
c    Input, double precision A(N), the vector to be printed.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      character ( len = * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,2x,g16.8)' ) i, a(i)
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
