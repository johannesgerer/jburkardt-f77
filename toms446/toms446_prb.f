      program main

c*********************************************************************72
c
cc TOMS446_PRB tests algorithm 446.
c
c  Modified:
c
c    12 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS446_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS446 library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
      call test05 ( )
      call test06 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS446_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests CHEBY, which computes Chebyshev series.
c
      implicit none

      integer n2
      integer nf
      integer npl
      integer nplmax

      parameter ( nf = 5 )
      parameter ( npl = 10 )

      parameter ( n2 = 2 * ( npl - 1 ) )
      parameter ( nplmax = npl )

      external functn
      double precision fxj(nf)
      double precision gc(n2)
      integer i
      integer j
      double precision x(nplmax,nf)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Test CHEBY, which computes the '
      write ( *, '(a)' ) '  Chebyshev series for several functions.'

      call cheby ( nf, npl, nplmax, n2, functn, x, fxj, gc )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '        Sin(x)      Cos(x)    Sin(2x)     Cos(2x)       X^5'
      write ( *, '(a)' ) ' '

      do i = 1, npl
        write ( *, '(5(2x,f10.4))' ) ( x(i,j), j = 1, nf )
      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests MULTPLY, which multiplies two Chebyshev series.
c
      implicit none

      integer n2
      integer nf
      integer npl
      integer nplmax

      parameter ( nf = 5 )
      parameter ( npl = 10 )

      parameter ( n2 = 2 * ( npl - 1 ) )
      parameter ( nplmax = npl )

      external functn
      double precision fxj(nf)
      double precision gc(n2)
      integer i
      integer j
      double precision x(nplmax,nf)
      double precision x1(npl)
      double precision x2(npl)
      double precision x3(npl)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  Test MLTPLY, which computes the '
      write ( *, '(a)' ) '  product of two Chebyshev series.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Multiply series for SIN(X) and COS(X)'
      write ( *, '(a)' ) '  and compare with series for 1/2*SIN(2X).'

      call cheby ( nf, npl, nplmax, n2, functn, x, fxj, gc )

      do i = 1, npl
        x1(i) = x(i,1)
        x2(i) = x(i,2)
        x(i,3) = 0.5D+00 * x(i,3)
      end do

      call mltply ( x1, x2, npl, x3 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '        Sin(x)      Cos(x)   1/2*Sin(2x)     RESULT'
      write ( *, '(a)' ) ' '

      do i = 1, npl
        write ( *, '(5(2x,f10.4))' ) ( x(i,j), j = 1, 3 ), x3(i)
      end do

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests ECHEB, which evaluates a Chebyshev series.
c
      implicit none

      integer n2
      integer nf
      integer npl
      integer nplmax

      parameter ( nf = 5 )
      parameter ( npl = 10 )

      parameter ( n2 = 2 * ( npl - 1 ) )
      parameter ( nplmax = npl )

      external functn
      double precision fval
      double precision fxj(nf)
      double precision gc(n2)
      integer i
      integer j
      integer k
      integer nx
      double precision x(nplmax,nf)
      double precision x2(npl)
      double precision xval

      nx = 6

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  Test ECHEB, which evaluates a '
      write ( *, '(a)' ) '  Chebyshev series.'

      call cheby ( nf, npl, nplmax, n2, functn, x, fxj, gc )

      do j = 1, nf

        do i = 1, npl
          x2(i) = x(i,j)
        end do

        write ( *, '(a)' ) ' '
        if ( j .eq. 1 ) then
          write ( *, '(a)' ) '  Sin(x)'
        else if ( j .eq. 2 ) then
          write ( *, '(a)' ) '  Cos(x)'
        else if ( j .eq. 3 ) then
          write ( *, '(a)' ) '  Sin(2x)'
        else if ( j .eq. 4 ) then
          write ( *, '(a)' ) '  Cos(2x)'
        else if ( j .eq. 5 ) then
          write ( *, '(a)' ) '  x^5'
        end if

        write ( *, '(a)' ) ' '

        do k = 1, nx

          xval = 2.0D+00 * dble ( k - 1 ) / dble ( nx - 1 ) - 1.0D+00

          call functn ( xval, fxj )

          call echeb ( xval, x2, npl, fval )

          write ( *, '(5(2x,f10.4))' ) xval, fxj(j), fval

        end do

      end do

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests EDCHEB, which evaluates the derivative of a Chebyshev series.
c
      implicit none

      integer n2
      integer nf
      integer npl
      integer nplmax

      parameter ( nf = 5 )
      parameter ( npl = 10 )

      parameter ( n2 = 2 * ( npl - 1 ) )
      parameter ( nplmax = npl )

      external functn
      double precision fval
      double precision fxj(nf)
      double precision gc(n2)
      integer i
      integer j
      integer k
      integer nx
      double precision x(nplmax,nf)
      double precision x2(npl)
      double precision xval

      nx = 6

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) '  Test EDCHEB, which evaluates the  '
      write ( *, '(a)' ) '  derivative of a Chebyshev series.'

      call cheby ( nf, npl, nplmax, n2, functn, x, fxj, gc )

      do j = 1, nf

        do i = 1, npl
          x2(i) = x(i,j)
        end do

        write ( *, '(a)' ) ' '
        if ( j .eq. 1 ) then
          write ( *, '(a)' ) '  d/dx Sin(x)'
        else if ( j .eq. 2 ) then
          write ( *, '(a)' ) '  d/dx Cos(x)'
        else if ( j .eq. 3 ) then
          write ( *, '(a)' ) '  d/dx Sin(2x)'
        else if ( j .eq. 4 ) then
          write ( *, '(a)' ) '  d/dx Cos(2x)'
        else if ( j .eq. 5 ) then
          write ( *, '(a)' ) '  d/dx x^5'
        end if

        write ( *, '(a)' ) ' '

        do k = 1, nx

          xval = 2.0D+00 * dble ( k - 1 ) / dble ( nx - 1 ) - 1.0D+00

          call functn_d ( xval, fxj )

          call edcheb ( xval, x2, npl, fval )

          write ( *, '(5(2x,f10.4))' ) xval, fxj(j), fval

        end do

      end do

      return
      end
      subroutine test05

c*********************************************************************72
c
cc TEST05 tests DFRNT, which computes the Chebyshev series of a derivative.
c
      implicit none

      integer n2
      integer nf
      integer npl
      integer nplmax

      parameter ( nf = 5 )
      parameter ( npl = 10 )

      parameter ( n2 = 2 * ( npl - 1 ) )
      parameter ( nplmax = npl )

      external functn
      double precision fxj(nf)
      double precision gc(n2)
      integer i
      integer j
      double precision x(nplmax,nf)
      double precision x2(npl)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) '  Test DFRNT, which computes the '
      write ( *, '(a)' ) '  Chebyshev series for the derivative'
      write ( *, '(a)' ) '  of several functions.'

      call cheby ( nf, npl, nplmax, n2, functn, x, fxj, gc )

      do j = 1, nf
        do i = 1, npl
          x2(i) = x(i,j)
        end do
        call dfrnt ( x2, npl, x2 )
        do i = 1, npl
          x(i,j) = x2(i)
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Chebyshev series for d/dx of:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '        Sin(x)      Cos(x)    Sin(2x)     Cos(2x)       X^5'
      write ( *, '(a)' ) ' '

      do i = 1, npl
        write ( *, '(5(2x,f10.4))' ) ( x(i,j), j = 1, nf )
      end do

      return
      end
      subroutine test06

c*********************************************************************72
c
cc TEST06 tests NTGRT, which computes the Chebyshev series of an indefinite integral.
c
      implicit none

      integer n2
      integer nf
      integer npl
      integer nplmax

      parameter ( nf = 5 )
      parameter ( npl = 10 )

      parameter ( n2 = 2 * ( npl - 1 ) )
      parameter ( nplmax = npl )

      external functn
      double precision fxj(nf)
      double precision gc(n2)
      integer i
      integer j
      double precision x(nplmax,nf)
      double precision x2(npl)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' ) '  Test NTGRT, which computes the '
      write ( *, '(a)' ) '  Chebyshev series for the indefinite'
      write ( *, '(a)' ) '  integral of several functions.'

      call cheby ( nf, npl, nplmax, n2, functn, x, fxj, gc )

      do j = 1, nf
        do i = 1, npl
          x2(i) = x(i,j)
        end do
        call ntgrt ( x2, npl, x2 )
        do i = 1, npl
          x(i,j) = x2(i)
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Series for indefinite integral of:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '        Sin(x)      Cos(x)    Sin(2x)     Cos(2x)       X^5'
      write ( *, '(a)' ) ' '

      do i = 1, npl
        write ( *, '(5(2x,f10.4))' ) ( x(i,j), j = 1, nf )
      end do

      return
      end
      subroutine functn ( x, fxj )

c*********************************************************************72
c
cc FUNCTN evaluates several functions at X.
c
      implicit none

      integer nf

      parameter ( nf = 5 )

      double precision fxj(nf)
      double precision x

      fxj(1) = sin ( x )
      fxj(2) = cos ( x )
      fxj(3) = sin ( 2.0D+00 * x )
      fxj(4) = cos ( 2.0D+00 * x )
      fxj(5) = x**5

      return
      end
      subroutine functn_d ( x, fxj )

c*********************************************************************72
c
cc FUNCTN_D evaluates the derivatives of several functions at X.
c
      implicit none

      integer nf

      parameter ( nf = 5 )

      double precision fxj(nf)
      double precision x

      fxj(1) =  cos ( x )
      fxj(2) = -sin ( x )
      fxj(3) =  2.0D+00 * cos ( 2.0D+00 * x )
      fxj(4) = -2.0D+00 * sin ( 2.0D+00 * x )
      fxj(5) =  5.0D+00 * x**4

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
