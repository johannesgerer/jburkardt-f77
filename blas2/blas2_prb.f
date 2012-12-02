      program main

c*********************************************************************72
c
cc MAIN is the main program for BLAS2_PRB.
c
c  Discussion:
c
c    BLAS2_PRB tests the BLAS2 routines.
c
c  Modified:
c
c    15 February 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BLAS2_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the BLAS2 library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
c
c  Terminate
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BLAS2_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests DGEMV.
c
      implicit none

      integer m
      integer n

      parameter ( m = 5 )
      parameter ( n = 5 )

      double precision a(m,n)
      double precision alpha
      double precision beta
      integer i
      integer incx
      integer incy
      integer j
      integer lda
      character trans
      double precision x(n)
      double precision y(m)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  For a general matrix A,'
      write ( *, '(a)' ) '  DGEMV computes '
      write ( *, '(a)' ) '  y := alpha * A * x + beta * y'

      trans = 'N'
      alpha = 2.0D+00
      lda = m
      incx = 1
      beta = 3.0D+00
      incy = 1

      do i = 1, m
        do j = 1, n
          if ( i == j ) then
            a(i,j) = 2.0D+00
          else if ( i == j - 1 .or. i == j + 1 ) then
            a(i,j) = -1.0D+00
          else
            a(i,j) = 0.0D+00
          end if
        end do
      end do

      do i = 1, n
        x(i) = dble ( i )
      end do

      do i = 1, m
        y(i) = dble ( 10 * i )
      end do

      call dgemv ( trans, m, n, alpha, a, lda, x, incx, beta, y, incy )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Result vector Y = '
      write ( *, '(a)' ) ' '

      do i = 1, m
        write ( *, '(2x,g14.6)' ) y(i)
      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests DGBMV.
c
      implicit none

      integer m
      integer n
      integer kl
      integer ku
      integer lda

      parameter ( m = 5 )
      parameter ( n = 5 )
      parameter ( kl = 1 )
      parameter ( ku = 1 )
      parameter ( lda = kl + 1 + ku )

      double precision a(lda,n)
      double precision alpha
      double precision beta
      integer i
      integer incx
      integer incy
      integer j
      integer jhi
      integer jlo
      character trans
      double precision x(n)
      double precision y(m)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  For a general band matrix A,'
      write ( *, '(a)' ) '  DGBMV computes '
      write ( *, '(a)' ) '  y := alpha * A * x + beta * y'

      trans = 'N'
      alpha = 2.0D+00
      incx = 1
      beta = 3.0D+00
      incy = 1

      do i = 1, m

        jlo = max ( 1, i - kl )
        jhi = min ( n, i + ku )

        do j = jlo, jhi

          if ( i == j ) then
            a(ku+1+i-j,j) = 2.0D+00
          else if ( i == j - 1 .or. i == j + 1 ) then
            a(ku+1+i-j,j) = -1.0D+00
          else
            a(ku+1+i-j,j) = 0.0D+00
          end if

        end do
      end do

      do i = 1, n
        x(i) = dble ( i )
      end do

      do i = 1, m
        y(i) = dble ( 10 * i )
      end do

      call dgbmv ( trans, m, n, kl, ku, alpha, a, lda, x,
     &  incx, beta, y, incy )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Result vector Y = '
      write ( *, '(a)' ) ' '

      do i = 1, m
        write ( *, '(2x,g14.6)' ) y(i)
      end do

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests DSYMV.
c
      implicit none

      integer n
      integer lda

      parameter ( n = 5 )
      parameter ( lda = n )

      double precision a(lda,n)
      double precision alpha
      double precision beta
      integer i
      integer incx
      integer incy
      integer j
      character uplo
      double precision x(n)
      double precision y(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  For a general symmetric matrix A,'
      write ( *, '(a)' ) '  DSYMV computes '
      write ( *, '(a)' ) '  y := alpha * A * x + beta * y'

      uplo = 'U'
      alpha = 2.0D+00
      incx = 1
      beta = 3.0D+00
      incy = 1

      do i = 1, n
        do j = 1, n
          if ( i == j ) then
            a(i,j) = 2.0D+00
          else if ( i == j - 1 ) then
            a(i,j) = -1.0D+00
          else
            a(i,j) = 0.0D+00
          end if
        end do
      end do

      do i = 1, n
        x(i) = dble ( i )
      end do

      do i = 1, n
        y(i) = dble ( 10 * i )
      end do

      call dsymv ( uplo, n, alpha, a, lda, x, incx, beta, y, incy )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Result vector Y = '
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,g14.6)' ) y(i)
      end do

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests DSBMV.
c
      implicit none

      integer m
      integer n
      integer k
      integer lda

      parameter ( m = 5 )
      parameter ( n = 5 )
      parameter ( k = 1 )
      parameter ( lda = k + 1 )

      double precision a(lda,n)
      double precision alpha
      double precision beta
      integer i
      integer incx
      integer incy
      integer j
      integer jhi
      integer jlo
      character uplo
      double precision x(n)
      double precision y(m)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) '  For a symmetric band matrix A,'
      write ( *, '(a)' ) '  DSBMV computes '
      write ( *, '(a)' ) '  y := alpha * A * x + beta * y'

      uplo = 'U'
      alpha = 2.0D+00
      incx = 1
      beta = 3.0D+00
      incy = 1

      do i = 1, m

        jhi = min ( n, i + k )

        do j = i, jhi

          if ( i == j ) then
            a(k+1+i-j,j) = 2.0D+00
          else if ( i == j - 1 ) then
            a(k+1+i-j,j) = -1.0D+00
          else
            a(k+1+i-j,j) = 0.0D+00
          end if

        end do
      end do

      do i = 1, n
        x(i) = dble ( i )
      end do

      do i = 1, m
        y(i) = dble ( 10 * i )
      end do

      call dsbmv ( uplo, n, k, alpha, a, lda, x, incx, beta, y, incy )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Result vector Y = '
      write ( *, '(a)' ) ' '

      do i = 1, m
        write ( *, '(2x,g14.6)' ) y(i)
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
c    16 September 2005
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

      character ( len = 8 ) date
      character ( len = 10 ) time

      call date_and_time ( date, time )

      write ( *, '(a8,2x,a10)' ) date, time

      return
      end
