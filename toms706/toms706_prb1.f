      program main

c*********************************************************************72
c
cc MAIN is the main program for TOMS706_PRB1.
c
c  Discussion:
c
c    The following program will approximate the integral of 
c      exp(x*x+y*y)
c    over the triangle (0.,0.),(2.,0.),(0.,2.) and print out the
c    values of the estimated integral, the estimated error, the number
c    of function evaluations, and IFAIL.
c
c    We need to choose LENVER  such that:
c    LENVER >= 3*((MAXPTS-37*NUMTRI)/(4*37)) + NUMTRI 
c    This simply means that we have enough space to achieve MAXPTS
c    function evaluations. By choosing LENVER bigger than 
c    this lower bound we may later change MAXPTS and RESTAR and restart
c    the computations from the point where the last run stopped.
c    The reason for stopping may have been that MAXPTS was too small
c    to achieve the requested error.
c
c    With our choice LENVER = 100 any value of MAXPTS
c    between 37 and 5068 will be accepted by the code. Choosing
c    MAXPTS = 2000 allows us to restart with a greater value
c    if necessary.
c
c  Modified:
c
c    09 December 2006
c
c  Author:
c
c    Jarle Berntsen, Terje Espelid
c
      implicit none

      integer lenver
      parameter ( lenver = 100 )
      integer mdiv
      parameter ( mdiv = 1 )
      integer minpts
      parameter ( minpts = 37 )
      integer numfun
      parameter ( numfun = 1 )
      integer numtri
      parameter ( numtri = 1 )

      integer nw
c     parameter ( nw = lenver*(2*numfun+1) +
c    &                max(32*mdiv,8*numtri)*numfun + 1 )
      parameter ( nw = lenver*(2*numfun+1) + 32*numfun + 1 )

      double precision abserr(numfun)
      double precision epsabs
      double precision epsrel
      external f
      integer ifail
      integer iwork(lenver+mdiv)
      integer maxpts
      integer neval
      integer restar
      double precision result(numfun)
      double precision ver(2,3,lenver)
      double precision work(nw)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS706_PRB1:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS706 library.'

      ver(1,1,1) = 0.0D+00
      ver(1,2,1) = 2.0D+00
      ver(1,3,1) = 0.0D+00
      ver(2,1,1) = 0.0D+00
      ver(2,2,1) = 0.0D+00
      ver(2,3,1) = 2.0D+00

      epsabs = 0.0D+00
      epsrel = 1.0D-05

      restar = 0
      maxpts = 2000

      call dcutri ( f, numfun, ver, numtri, minpts, maxpts, epsabs,
     &  epsrel, lenver, nw, restar, result, abserr, neval,
     &  ifail, work, iwork )

      write ( *, '(a)' ) ' '
      write ( *, '(a,d25.16)' ) '  RESULT= ', result(1)
      write ( *, '(a,d25.16)' ) '  ERROR = ', abserr(1)
      write ( *, '(a,i8)' ) '  NEVAL = ', neval
      write ( *, '(a,i8)' ) '  IFAIL = ', ifail

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DCUTRI_PRB1:'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
      subroutine f ( x, numfun, funvls )

c*********************************************************************72
c
cc F evaluates the integrand function.
c
c  Modified:
c
c    09 December 2006
c
c  Parameters:
c
c    Input, double precision X(2), the evaluation point.
c
c    Input, integer NUMFUN, the number of functions.
c
c    Output, double precision FUNVLS(NUMFUN), the value of the
c    functions at X.
c
      implicit none

      integer numfun

      double precision funvls(numfun)
      double precision x(2)

      funvls(1) = exp ( x(1) * x(1) + x(2) * x(2) )

      return
      end

