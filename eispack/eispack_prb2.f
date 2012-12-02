      program main

c*********************************************************************72
c
cc MAIN is the main program for EISPACK_PRB2.
c
c  Discussion:
c
c    EISPACK_PRB2 does some symmetric eigenproblem tests on EISPACK.
c
c  Modified:
c
c    05 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      integer n

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EISPACK_PRB2'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the EISPACK library.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Do some symmetric eigenproblem tests.'

      n = 1

      do i = 1, 4
        n = n * 4
        call test01 ( n )
      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EISPACK_PRB2'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)') ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( n )

c*********************************************************************72
c
cc TEST01 tests RS.
c
c  Modified:
c
c    05 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n

      double precision a(n,n)
      double precision aq(n,n)
      double precision a2(n,n)
      double precision fv1(n)
      double precision fv2(n)
      integer i
      integer ierror
      integer j
      integer k
      double precision lambda(n)
      double precision lambda2(n)
      double precision lambda_dev
      parameter ( lambda_dev = 1.0D+00 )
      double precision lambda_max
      double precision lambda_mean
      parameter ( lambda_mean = 1.0D+00 )
      double precision lambda_min
      integer matz
      double precision q(n,n)
      double precision q2(n,n)
      double precision r(n,n)
      integer seed
      double precision t1
      double precision t2
      double precision time_setup
      double precision time_solve

      seed = 12345

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  RS computes the eigenvalues and'
      write ( *, '(a)' ) '  eigenvectors of a real symmetric matrix.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix order = ', n
      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) 
     &  '  Initialize random number generator using SEED = ', seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  SYMM_TEST will give us a symmetric matrix'
      write ( *, '(a)' ) '  with known eigenstructure.'

      call cpu_time ( t1 )

      call r8symm_test ( n, lambda_mean, lambda_dev, seed, a, q, 
     &  lambda )

      call cpu_time ( t2 )

      time_setup = t2 - t1

      if ( n <= 5 ) then

        call r8mat_print ( n, n, a, '  The matrix A:' )

        call r8mat_print ( n, n, q, '  The eigenvector matrix Q:' )

      end if

      lambda_min = lambda(1)
      lambda_max = lambda(1)
      do i = 2, n
        lambda_min = min ( lambda_min, lambda(i) )
        lambda_max = max ( lambda_max, lambda(i) )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  LAMBDA_MIN = ', lambda_min
      write ( *, '(a,g14.6)' ) '  LAMBDA_MAX = ', lambda_max

      if ( n .le. 10 ) then

        call r8vec_print ( n, lambda, 
     &  '  The eigenvalue vector LAMBDA:' )

      end if
c
c  Verify the claim that A*Q = Q * LAMBDA.
c
      if ( n .le. 5 ) then

        do i = 1, n
          do j = 1, n
            aq(i,j) = 0.0D+00
            do k = 1, n
              aq(i,j) = aq(i,j) + a(i,k) * q(k,j)
            end do
          end do
        end do

        do j = 1, n
          lambda2(j) = 0.0D+00
          do i = 1, n
            lambda2(j) = lambda2(j) + aq(i,j)**2
          end do
          lambda2(j) = sqrt ( lambda2(j) )
        end do

        call r8vec_print ( n, lambda2, '  The column norms of A*Q:' )

      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Now call EISPACK routine RS'
      write ( *, '(a)' ) '  and see if it can recover Q and LAMBDA.'
c
c  Copy the matrix.
c
      do i = 1, n
        do j = 1, n
          a2(i,j) = a(i,j)
        end do
      end do

      matz = 1

      call cpu_time ( t1 )

      call rs ( n, n, a2, lambda2, matz, q2, fv1, fv2, ierror )

      call cpu_time ( t2 )

      time_solve = t2 - t1 

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST01 - Warning!'
        write ( *, '(a,i8)' ) '  RS returned  error flag = ', ierror
        return
      end if

      lambda_min = lambda2(1)
      lambda_max = lambda2(1)
      do i = 2, n
        lambda_min = min ( lambda_min, lambda2(i) )
        lambda_max = max ( lambda_max, lambda2(i) )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  LAMBDA_MIN = ', lambda_min
      write ( *, '(a,g14.6)' ) '  LAMBDA_MAX = ', lambda_max

      if ( n .le. 10 ) then
        call r8vec_print ( n, lambda2, 
     &  '  The computed eigenvalues Lambda:' )
      end if

      if ( matz .ne. 0 ) then

        if ( n .le. 5 ) then
          call r8mat_print ( n, n, q2, '  The eigenvector matrix:' )

          do i = 1, n
            do j = 1, n
              r(i,j) = 0.0D+00
              do k = 1, n
                r(i,j) = r(i,j) + a(i,k) * q2(k,j)
              end do
            end do
          end do

          do j = 1, n
            do i = 1, n
              r(i,j) = r(i,j) - lambda2(j) * q2(i,j)
            end do
          end do

          call r8mat_print ( n, n, r, 
     &    '  The residual (A-Lambda*I)*X:' )

        end if

      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Setup time = ', time_setup
      write ( *, '(a,g14.6)' ) '  Solve time = ', time_solve

      return
      end
      subroutine r8mat_print ( m, n, a, title )

c*********************************************************************72
c
cc R8MAT_PRINT prints an R8MAT.
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

      call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

      return
      end
      subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi,
     &  title )

c*********************************************************************72
c
cc R8MAT_PRINT_SOME prints some of an R8MAT.
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
c    Input, character ( len = * ) TITLE, an optional title.
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
      integer title_length

      title_length = len_trim ( title )

      if ( 0 .lt. title_length ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) title(1:title_length)
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

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

        end do

      end do

      write ( *, '(a)' ) ' '

      return
      end
      subroutine r8vec_print ( n, a, title )

c*********************************************************************72
c
cc R8VEC_PRINT prints an R8VEC.
c
c  Discussion:
c
c    An R8VEC is an array of double precision real values.
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
c    Input, character * ( * ) TITLE, an optional title.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      integer s_len_trim
      character ( len = * ) title
      integer title_length

      title_length = s_len_trim ( title )
      if ( 0 .lt. title_length ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) title(1:title_length)
      end if

      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,2x,g16.8)' ) i, a(i)
      end do

      return
      end
      function s_len_trim ( s )

c*********************************************************************72
c
cc S_LEN_TRIM returns the length of a string to the last nonblank.
c
c  Modified:
c
c    05 March 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*(*) S, a string.
c
c    Output, integer S_LEN_TRIM, the length of the string to the last nonblank.
c
      implicit none

      integer i
      character*(*) s
      integer s_len_trim

      do i = len ( s ), 1, -1

        if ( s(i:i) .ne. ' ' ) then
          s_len_trim = i
          return
        end if

      end do

      s_len_trim = 0

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
