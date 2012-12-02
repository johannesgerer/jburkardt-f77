      subroutine cholesky ( a, n, nn, u, nullty, ifault )

c*********************************************************************72
c
cc CHOLESKY computes the Cholesky factorization of a PDS matrix.
c
c  Discussion:
c
c    For a positive definite symmetric matrix A, the Cholesky factor U
c    is an upper triangular matrix such that A = U' * U.
c
c    This routine was originally named "CHOL", but that conflicted with
c    a built in MATLAB routine name.
c
c    The missing initialization "II = 0" has been added to the code.
c
c  Modified:
c
c    12 February 2008
c
c  Author:
c
c    Original FORTRAN77 version by Michael Healy.
c    Modifications by AJ Miller.
c    Modifications by John Burkardt.
c
c  Reference:
c
c    Michael Healy,
c    Algorithm AS 6:
c    Triangular decomposition of a symmetric matrix,
c    Applied Statistics,
c    Volume 17, Number 2, 1968, pages 195-197.
c
c  Parameters:
c
c    Input, double precision A((N*(N+1))/2), a positive definite matrix 
c    stored by rows in lower triangular form as a one dimensional array, 
c    in the sequence
c    A(1,1),
c    A(2,1), A(2,2),
c    A(3,1), A(3,2), A(3,3), and so on.
c
c    Input, integer N, the order of A.
c
c    Input, integer NN, the dimension of A, (N*(N+1))/2.
c
c    Output, double precision U((N*(N+1))/2), an upper triangular matrix,
c    stored by columns, which is the Cholesky factor of A.  The program is
c    written in such a way that A and U can share storage.
c
c    Output, integer NULLTY, the rank deficiency of A.  If NULLTY is zero,
c    the matrix is judged to have full rank.
c
c    Output, integer IFAULT, an error indicator.
c    0, no error was detected;
c    1, if N < 1;
c    2, if A is not positive semi-definite.
c    3, if NN < (N*(N+1))/2.
c
c  Local Parameters:
c
c    Local, double precision ETA, should be set equal to the smallest positive 
c    value such that 1.0 + ETA is calculated as being greater than 1.0 in the 
c    accuracy being used.
c
      implicit none

      integer nn

      double precision a(nn)
      double precision eta
      parameter ( eta = 1.0D-09 )
      integer i
      integer icol
      integer ifault
      integer ii
      integer irow
      integer j
      integer k
      integer kk
      integer l
      integer m
      integer n
      integer nullty
      double precision rsq
      double precision u(nn)
      double precision w
      double precision x

      ifault = 0
      nullty = 0

      if ( n .le. 0 ) then
        ifault = 1
        return
      end if

      if ( nn .lt. ( n * ( n + 1 ) ) / 2 ) then
        ifault = 3
        return
      end if

      j = 1
      k = 0
      ii = 0
c
c  Factorize column by column, ICOL = column number.
c
      do icol = 1, n

        ii = ii + icol
        x = eta * eta * a(ii)
        l = 0
        kk = 0
c
c  IROW = row number within column ICOL.
c
        do irow = 1, icol

          kk = kk + irow
          k = k + 1
          w = a(k)
          m = j

          do i = 1, irow - 1
            l = l + 1
            w = w - u(l) * u(m)
            m = m + 1
          end do

          l = l + 1

          if ( irow .eq. icol ) then
            go to 50
          end if

          if ( u(l) .ne. 0.0D+00 ) then

            u(k) = w / u(l)

          else

            u(k) = 0.0D+00

            if ( abs ( x * a(k) ) .lt. w * w ) then
              ifault = 2
              return
            end if

          end if

        end do
c
c  End of row, estimate relative accuracy of diagonal element.
c
   50   continue

        if ( abs ( w ) .le. abs ( eta * a(k) ) ) then

          u(k) = 0.0D+00
          nullty = nullty + 1
       
        else

          if ( w .lt. 0.0D+00 ) then
            ifault = 2
            return
          end if

          u(k) = sqrt ( w )

        end if

        j = j + icol

      end do

      return
      end
      subroutine syminv ( a, n, c, w, nullty, ifault )

c*********************************************************************72
c
cc SYMINV computes the inverse of a symmetric matrix.
c
c  Modified:
c
c    01 February 2008
c
c  Author:
c
c    Michael Healy
c    Modifications by John Burkardt
c
c  Reference:
c
c    Michael Healy,
c    Algorithm AS 7:
c    Inversion of a Positive Semi-Definite Symmetric Matrix,
c    Applied Statistics,
c    Volume 17, Number 2, 1968, pages 198-199.
c
c  Parameters:
c
c    Input, double precision A((N*(N+1))/2), a positive definite matrix stored 
c    by rows in lower triangular form as a one dimensional array, in the sequence
c    A(1,1),
c    A(2,1), A(2,2),
c    A(3,1), A(3,2), A(3,3), and so on.
c
c    Input, integer N, the order of A.
c
c    Output, double precision C((N*(N+1))/2), the inverse of A, or generalized
c    inverse if A is singular, stored using the same storage scheme employed
c    for A.  The program is written in such a way that A and U can share storage.
c
c    Workspace, double precision W(N).
c
c    Output, integer NULLTY, the rank deficiency of A.  If NULLTY is zero,
c    the matrix is judged to have full rank.
c
c    Output, integer IFAULT, error indicator.
c    0, no error detected.
c    1, N < 1.
c    2, A is not positive semi-definite.
c
      implicit none

      integer n

      double precision a((n*(n+1))/2)
      double precision c((n*(n+1))/2)
      integer i
      integer icol
      integer ifault
      integer irow
      integer j
      integer jcol
      integer k
      integer l
      integer mdiag
      integer ndiag
      integer nn
      integer nrow
      integer nullty
      double precision w(n)
      double precision x

      ifault = 0

      if ( n .le. 0 ) then
        ifault = 1
        return
      end if

      nrow = n
c
c  Compute the Cholesky factorization of A.
c  The result is stored in C.
c
      nn = ( n * ( n + 1 ) ) / 2

      call cholesky ( a, n, nn, c, nullty, ifault )

      if ( ifault .ne. 0 ) then
        return
      end if
c
c  Invert C and form the product (Cinv)' * Cinv, where Cinv is the inverse
c  of C, row by row starting with the last row.
c  IROW = the row number, 
c  NDIAG = location of last element in the row.
c
      irow = nrow
      ndiag = nn

10    continue
c
c  Special case, zero diagonal element.
c
      if ( c(ndiag) .eq. 0.0D+00 ) then

        l = ndiag
        do j = irow, nrow
          c(l) = 0.0D+00
          l = l + j
        end do

      else

        l = ndiag
        do i = irow, nrow
          w(i) = c(l)
          l = l + i
        end do

        icol = nrow
        jcol = nn
        mdiag = nn

   30   continue

        l = jcol
 
        if ( icol .eq. irow ) then
          x = 1.0D+00 / w(irow)
        else
          x = 0.0D+00
        end if

        k = nrow

   40   continue

        if ( irow .lt. k ) then
          x = x - w(k) * c(l)
          k = k - 1
          l = l - 1
          if ( l .gt. mdiag ) then
            l = l - k + 1
          end if
          go to 40
        end if

        c(l) = x / w(irow)

        if ( irow .lt. icol ) then
          mdiag = mdiag - icol
          icol = icol - 1
          jcol = jcol - 1
          go to 30
        end if

      end if

      ndiag = ndiag - irow
      irow = irow - 1

      if ( 0 .lt. irow ) then
        go to 10
      end if

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
