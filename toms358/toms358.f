      subroutine csvd ( a, mmax, nmax, m, n, p, nu, nv, s, u, v )

c*********************************************************************72
c
cc CSVD computes the singular value decomposition of an M by N complex matrix.
c
c  Discussion:
c
c    This routine requires that N <= M.
c
c    The singular value decomposition of a complex M by N matrix A
c    has the form
c
c      A = U S V*
c
c    where 
c
c      U is an M by M unitary matrix,
c      S is an M by N diagonal matrix,
c      V is an N by N unitary matrix.
c
c    Moreover, the entries of S are nonnegative and occur on the diagonal
c    in descending order.
c
c    Several internal arrays are dimensioned under the assumption
c    that N <= 100.
c
c  Modified:
c
c    03 May 2007
c
c  Author:
c
c    Peter Businger, Gene Golub
c
c  Reference:
c
c    Peter Businger, Gene Golub,
c    Algorithm 358:
c    Singular Value Decomposition of a Complex Matrix,
c    Communications of the ACM,
c    Volume 12, Number 10, October 1969, pages 564-565.
c
c  Parameters:
c
c    Input/output, complex A(MMAX,*), the M by N matrix, which may be
c    augmented by P extra columns to which the transformation U*
c    is to be applied.  On output, A has been overwritten, and
c    if 0 < P, columns N+1 through N+P have been premultiplied by U*.
c
c    Input, integer MMAX, the leading dimension of the arrays A
c    and U.
c
c    Input, integer NMAX, the leading dimension of V, and perhaps
c    the second dimension of A and U.
c
c    Input, integer M, N, the number of rows and columns in A.
c    It must be the case that 1 <= N <= M.  Several internal arrays are
c    dimensioned under the assumption that N <= NBIG, where NBIG
c    is an internal parameter, currently set to 100.
c
c    Input, integer P, the number of vectors, stored in A(*,N+1:N+P),
c    to which the transformation U* should be applied.
c
c    Input, integer NU, the number of columns of U to compute.
c
c    Input, integer NV, the number of columns of V to compute.
c
c    Output, real S(N), the computed singular values.
c
c    Output, complex U(MMAX,NU), the first NU columns of U.
c
c    Output, complex V(NMAX,NV), the first NV columns of V.
c
c  Local Parameters:
c
c    Local, real ETA, the relative machine precision.
c    The original text uses ETA = 1.5E-8.
c
c    Local, integer NBIG, is a parameter used to dimension work arrays.
c    The size of NBIG limits the maximum possible size of N that can
c    be handled.  If you want to work with values of N that are larger,
c    simply increase the value assigned to NBIG below.
c
c    Local, real TOL, the smallest normalized positive number, divided by ETA.
c    The original test uses TOL = 1.E-31.
c
      implicit none

      integer mmax
      integer nmax

      integer nbig
      parameter ( nbig = 100 )

      complex a(mmax,*)
      real b(nbig)
      real c(nbig)
      real cs
      real eps
      real eta
      real f
      real g
      real h
      integer i
      integer j
      integer k
      integer kk
      integer k1
      integer l
      integer l1
      integer ll
      integer m
      integer n
      integer nu
      integer nv
      integer p
      complex q
      complex r
      real s(*)
      real sn
      real t(nbig)
      real tol
      complex u(mmax,*)
      complex v(nmax,*)
      real w
      real x
      real y
      real z

      save eta
      save tol

      data eta / 1.1920929e-07 /
      data tol / 1.5e-31 /
c
c  Check N.
c
      if ( n < 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CSVD - Fatal error!'
        write ( *, '(a)' ) '  Input N < 1.'
        stop
      else if ( nbig < n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CSVD - Fatal error!'
        write ( *, '(a)' ) '  NBIG < N.'
        stop
      end if
c
c  Check M.
c
      if ( m < 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CSVD - Fatal error!'
        write ( *, '(a)' ) '  Input M < 1.'
        stop
      else if ( m < n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CSVD - Fatal error!'
        write ( *, '(a)' ) '  M < N.'
        stop
      end if
c
c  Householder reduction.
c
      c(1) = 0.0e+00
      k = 1

10    continue

      k1 = k + 1
c
c  Elimination of A(I,K), I = K+1, ..., M.
c
      z = 0.0e+00
      do i = k, m
        z = z + ( real ( a(i,k) ) )**2 + ( aimag ( a(i,k) ) )**2
      end do

      b(k) = 0.0e+00

      if ( tol .lt. z ) then

        z = sqrt ( z )
        b(k) = z
        w = cabs ( a(k,k) )

        if ( w .eq. 0.0e+00 ) then
          q = cmplx ( 1.0e+00, 0.0e+00 )
        else
          q = a(k,k) / w
        end if

        a(k,k) = q * ( z + w )

        if ( k .ne. n + p ) then

          do j = k1, n + p

            q = cmplx ( 0.0e+00, 0.0e+00 )
            do i = k, m
              q = q + conjg ( a(i,k) ) * a(i,j)
            end do
            q = q / ( z * ( z + w ) )

            do i = k, m
              a(i,j) = a(i,j) - q * a(i,k)
            end do

          end do
c
c  Phase transformation.
c
          q = -conjg ( a(k,k) ) / cabs ( a(k,k) )

          do j = k1, n + p
            a(k,j) = q * a(k,j)
          end do

        end if

      end if
c
c  Elimination of A(K,J), J = K+2, ..., N
c
      if ( k .eq. n ) then
        go to 140
      end if

      z = 0.0e+00
      do j = k1, n
        z = z + ( real ( a(k,j) ) )**2 + ( aimag ( a(k,j) ) )**2
      end do

      c(k1) = 0.0e+00

      if ( tol .lt. z ) then

        z = sqrt ( z )
        c(k1) = z
        w = cabs ( a(k,k1) )

        if ( w .eq. 0.0e+00 ) then
          q = cmplx ( 1.0e+00, 0.0e+00 )
        else
          q = a(k,k1) / w
        end if

        a(k,k1) = q * ( z + w )

        do i = k1, m

          q = cmplx ( 0.0e+00, 0.0e+00 )

          do j = k1, n
            q = q + conjg ( a(k,j) ) * a(i,j)
          end do

          q = q / ( z * ( z + w ) )

          do j = k1, n
            a(i,j) = a(i,j) - q * a(k,j)
          end do

        end do
c
c  Phase transformation.
c
        q = -conjg ( a(k,k1) ) / cabs ( a(k,k1) )
        do i = k1, m
          a(i,k1) = a(i,k1) * q
        end do

      end if

      k = k1
      go to 10
c
c  Tolerance for negligible elements.
c
140   continue

      eps = 0.0e+00
      do k = 1, n
        s(k) = b(k)
        t(k) = c(k)
        eps = amax1 ( eps, s(k) + t(k) )
      end do

      eps = eps * eta
c
c  Initialization of U and V.
c
      if ( 0 .lt. nu ) then

        do j = 1, nu
          do i = 1, m
            u(i,j) = cmplx ( 0.0e+00, 0.0e+00 )
          end do
          u(j,j) = cmplx ( 1.0e+00, 0.0e+00 )
        end do

      end if

      if ( 0 .lt. nv ) then

        do j = 1, nv
          do i = 1, n
            v(i,j) = cmplx ( 0.0e+00, 0.0e+00 )
          end do
          v(j,j) = cmplx ( 1.0e+00, 0.0e+00 )
        end do

      end if
c
c  QR diagonalization.
c
      do kk = 1, n

        k = n + 1 - kk
c
c  Test for split.
c
220     continue

        do ll = 1, k

          l = k + 1 - ll

          if ( abs ( t(l) ) .le. eps ) then
            go to 290
          end if

          if ( abs ( s(l-1) ) .le. eps ) then
            go to 240
          end if

        end do
c
c  Cancellation of E(L).
c
240     continue

        cs = 0.0e+00
        sn = 1.0e+00
        l1 = l - 1

        do i = l, k

          f = sn * t(i)
          t(i) = cs * t(i)

          if ( abs ( f ) .le. eps ) then
            go to 290
          end if

          h = s(i)
          w = sqrt ( f * f + h * h )
          s(i) = w
          cs = h / w
          sn = - f / w

          if ( 0 .lt. nu ) then

            do j = 1, n
              x = real ( u(j,l1) )
              y = real ( u(j,i) )
              u(j,l1) = cmplx ( x * cs + y * sn, 0.0e+00 )
              u(j,i)  = cmplx ( y * cs - x * sn, 0.0e+00 )
            end do

          end if

          if ( p .ne. 0 ) then

            do j = n + 1, n + p
              q = a(l1,j)
              r = a(i,j)
              a(l1,j) = q * cs + r * sn
              a(i,j)  = r * cs - q * sn
            end do

          end if

        end do
c
c  Test for convergence.
c
290     continue

        w = s(k)

        if ( l .eq. k ) then
          go to 360
        end if
c
c  Origin shift.
c
        x = s(l)
        y = s(k-1)
        g = t(k-1)
        h = t(k)
        f = ( ( y - w ) * ( y + w ) + ( g - h ) * ( g + h ) ) 
     &    / ( 2.0e+00 * h * y )
        g = sqrt ( f * f + 1.0e+00 )
        if ( f .lt. 0.0e+00 ) then
          g = -g
        end if
        f = ( ( x - w ) * ( x + w ) + ( y / ( f + g ) - h ) * h ) / x
c
c  QR step.
c
        cs = 1.0e+00
        sn = 1.0e+00
        l1 = l + 1

        do i = l1, k

          g = t(i)
          y = s(i)
          h = sn * g
          g = cs * g
          w = sqrt ( h * h + f * f )
          t(i-1) = w
          cs = f / w
          sn = h / w
          f = x * cs + g * sn
          g = g * cs - x * sn
          h = y * sn
          y = y * cs

          if ( 0 .lt. nv ) then

            do j = 1, n
              x = real ( v(j,i-1) )
              w = real ( v(j,i) )
              v(j,i-1) = cmplx ( x * cs + w * sn, 0.0e+00 )
              v(j,i)   = cmplx ( w * cs - x * sn, 0.0e+00 )
            end do

          end if

          w = sqrt ( h * h + f * f )
          s(i-1) = w
          cs = f / w
          sn = h / w
          f = cs * g + sn * y
          x = cs * y - sn * g

          if ( 0 .lt. nu ) then

            do j = 1, n
              y = real ( u(j,i-1) )
              w = real ( u(j,i) )
              u(j,i-1) = cmplx ( y * cs + w * sn, 0.0e+00 )
              u(j,i)   = cmplx ( w * cs - y * sn, 0.0e+00 )
            end do

          end if

          if ( p .ne. 0 ) then

            do j = n + 1, n + p
              q = a(i-1,j)
              r = a(i,j)
              a(i-1,j) = q * cs + r * sn
              a(i,j)   = r * cs - q * sn
            end do

          end if

        end do

        t(l) = 0.0e+00
        t(k) = f
        s(k) = x
        go to 220
c
c  Convergence.
c
360     continue

        if ( w .lt. 0.0e+00 ) then

          s(k) = - w

          if ( 0 .lt. nv ) then

            do j = 1, n
              v(j,k) = - v(j,k)
            end do

          end if
 
        end if

      end do
c
c  Sort the singular values.
c
      do k = 1, n

        g = - 1.0e+00
        j = k

        do i = k, n
          if ( g .lt. s(i) ) then
            g = s(i)
            j = i
          end if
        end do

        if ( j .ne. k ) then

          s(j) = s(k)
          s(k) = g
c
c  Interchange V(1:N,J) and V(1:N,K).
c
          if ( 0 .lt. nv ) then

            do i = 1, n
              q      = v(i,j)
              v(i,j) = v(i,k)
              v(i,k) = q
            end do

          end if
c
c  Interchange U(1:N,J) and U(1:N,K).
c
          if ( 0 .lt. nu ) then
  
            do i = 1, n
              q      = u(i,j)
              u(i,j) = u(i,k)
              u(i,k) = q
            end do

          end if
c
c  Interchange A(J,N1:NP) and A(K,N1:NP).
c
          if ( p .ne. 0 ) then

            do i = n + 1, n + p
              q      = a(j,i)
              a(j,i) = a(k,i)
              a(k,i) = q
            end do

          end if

        end if

      end do
c
c  Back transformation.
c
      if ( 0 .lt. nu ) then

        do kk = 1, n

          k = n + 1 - kk

          if ( b(k) .ne. 0.0e+00 ) then

            q = -a(k,k) / cabs ( a(k,k) )

            do j = 1, nu
              u(k,j) = q * u(k,j)
            end do

            do j = 1, nu

              q = cmplx ( 0.0e+00, 0.0e+00 )

              do i = k, m
                q = q + conjg ( a(i,k) ) * u(i,j)
              end do

              q = q / ( cabs ( a(k,k) ) * b(k) )

              do i = k, m
                u(i,j) = u(i,j) - q * a(i,k)
              end do

            end do

          end if

        end do

      end if

      if ( 0 .lt. nv ) then

        if ( 2 .lt. n ) then

          do kk = 2, n

            k = n + 1 - kk
            k1 = k + 1

            if ( c(k1) .ne. 0.0e+00 ) then

              q = -conjg ( a(k,k1) ) / cabs ( a(k,k1) )

              do j = 1, nv
                v(k1,j) = q * v(k1,j)
              end do

              do j = 1, nv

                q = cmplx ( 0.0e+00, 0.0e+00 )

                do i = k1, n
                  q = q + a(k,i) * v(i,j)
                end do

                q = q / ( cabs ( a(k,k1) ) * c(k1) )

                do i = k1, n
                  v(i,j) = v(i,j) - q * conjg ( a(k,i) )
                end do

              end do

            end if

          end do

        end if

      end if

      return
      end
