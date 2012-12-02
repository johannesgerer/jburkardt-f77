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
c    Michael Healy
c    Modifications by AJ Miller.
c    Modifications by John Burkardt
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
c    Input, integer NN, the dimension of the array used to store A, 
c    which should be at least (N*(N+1))/2.
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
      subroutine subchl ( a, b, n, u, nullty, ifault, ndim, det )

c*********************************************************************72
c 
cc SUBCHL computes the Cholesky factorization of a (subset of a ) PDS matrix.
c
c  Modified:
c
c    11 February 2008
c
c  Author:
c
c    Michael Healy, PR Freeman
c    Modifications by John Burkardt
c
c  Reference:
c
c    PR Freeman,
c    Remark AS R44:
c    A Remark on AS 6 and AS7: Triangular decomposition of a symmetric matrix
c    and Inversion of a positive semi-definite symmetric matrix,
c    Applied Statistics,
c    Volume 31, Number 3, 1982, pages 336-339.
c
c    Michael Healy,
c    Algorithm AS 6:
c    Triangular decomposition of a symmetric matrix,
c    Applied Statistics,
c    Volume 17, Number 2, 1968, pages 195-197.
c
c  Parameters:
c
c    Input, double precision A((M*(M+1))/2), a positive definite matrix 
c    stored by rows in lower triangular form as a one dimensional array, 
c    in the sequence
c    A(1,1),
c    A(2,1), A(2,2),
c    A(3,1), A(3,2), A(3,3), and so on.  
c    In the simplest case, M, the order of A, is equal to N.
c
c    Input, integer B(N), indicates the order in which the rows and columns
c    of A are to be used.  In the simplest case, B = (1,2,3...,N).
c
c    Input, integer N, the order of the matrix, that is, the matrix formed
c    by using B to select N rows and columns of A.
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
c
c    Input, integer NDIM, the dimension of A and U, which might be presumed
c    to be (N*(N+1))/2.
c
c    Output, double precision DET, the determinant of the matrix.
c
      implicit none

      integer n
      integer ndim

      double precision a(ndim)
      integer b(n)
      double precision det
      double precision eta
      parameter ( eta = 1.0D-09 )
      integer i
      integer icol
      integer ifault
      integer ii
      integer ij
      integer irow
      integer j
      integer jj
      integer k
      integer kk
      integer l
      integer m
      integer nullty
      double precision u(ndim)
      double precision w
      double precision x

      ifault = 0
      nullty = 0
      det = 1.0D+00

      if ( n .le. 0 ) then
        ifault = 1
        return
      end if

      ifault = 2
      j = 1
      k = 0

      do icol = 1, n

        ij = ( b(icol) * ( b(icol) - 1 ) ) / 2
        ii = ij + b(icol)
        x = eta * eta * a(ii)
        l = 0

        do irow = 1, icol

          kk = ( b(irow) * ( b(irow) + 1 ) ) / 2
          k = k + 1
          jj = ij + b(irow)
          w = a(jj)
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

            if ( abs ( x * a(kk) ) .lt. w * w ) then
              ifault = 2
              return
            end if

            u(k) = 0.0D+00

          end if

        end do

 50     continue

        if ( abs ( eta * a(kk) ) .le. abs ( w ) ) then

          if ( w .lt. 0.0D+00 ) then
            ifault = 2
            return
          end if

          u(k) = sqrt ( w )

        else

          u(k) = 0.0D+00
          nullty = nullty + 1

        end if

        j = j + icol
        det = det * u(k) * u(k)

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
