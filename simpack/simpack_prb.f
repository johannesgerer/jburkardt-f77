      program main

c*********************************************************************72
c
cc MAIN is the main program for SIMPACK_PRB.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Alan Genz
c
c  Local Parameters:
c
c    Local, integer NDIM, the spatial dimension.
c
      implicit none

      integer ndim
      parameter ( ndim = 5 )
      integer nf
      parameter ( nf = 9 )
      integer lwk
      parameter ( lwk = 20000 )

      double precision absest(nf)
      double precision absreq
      double precision finest(nf)
      external ftest
      integer inform
      integer key
      integer m
      integer maxcls
      integer mincls
      integer n
      integer neval
      integer num
      double precision relreq
      integer subs
      double precision vertxs(lwk)

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SIMPACK_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the SIMPACK library.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Estimate the integral of a vector function'
      write ( *, '(a)' ) '  over an NDIM-dimensional region which is'
      write ( *, '(a)' ) '  a collection of simplexes.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Spatial dimension NDIM = ', ndim
      write ( *, '(a,i8)' ) '  Vector dimension NF =    ', nf
c
c  Initialize vertices for simplex 0 < x_1 < x_2 ... x_5 < 1.
c
      do m = 0, ndim
        do n = 1, ndim-m
          vertxs(n+m*ndim) = 0.0D+00
        end do
        do n = ndim-m+1, ndim
          vertxs(n+m*ndim) = 1.0D+00
        end do
      end do
c
c  Initialize other parameters.
c
      mincls = 0
      num = 20000
      maxcls = num
      key = 4
      absreq = 0.0D+00
      relreq = 1.0D-04
      subs = 1
c
c  Call SMPINT and print results
c
      call smpint ( ndim, nf, mincls, maxcls, ftest, absreq, relreq,
     &  key, subs, lwk, vertxs, 0, finest, absest, neval, inform )
c
c  Print results
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Results:'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Function evaluations NEVAL =    ', neval
      write ( *, '(a,i8)' ) '  Integration rule selector KEY = ', key
      write ( *, '(a,i8)' ) '  Error flag INFORM =             ', inform
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    N  Estimated error    Integral'
      write ( *, '(a)' ) ' '

      do n = 1, nf
        print 99998, n, absest(n), finest(n)
99998   format (3x, i2, 2f15.10)
      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SMPINT_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine ftest ( ndim, z, nfun, f )

c*********************************************************************72
c
cc FTEST evaluates the integrand function.
c
      implicit none

      integer ndim
      integer nfun

      integer i
      double precision f(nfun)
      double precision value
      double precision z(ndim)

      value = 0.0D+00
      do i = 1, ndim
         value = value - z(i) * z(i) / 2.0D+00
      end do
      value = exp ( value )

      f(1) = value
      f(2) = value * z(1)
      f(3) = value * z(1) * z(2)
      f(4) = value * z(1)**2
      f(5) = value * z(1)**3
      f(6) = value * z(1)**4
      f(7) = value * z(1)**5
      f(8) = value * z(1)**6
      f(9) = value * z(1)**7

      return
      end
