      program main

c*********************************************************************72
c
cc MAIN is the main program for UNCMIN_PRB.
c
c  Discussion:
c
c    UNCMIN_PRB tests the routines in UNCMIN.
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82,
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
      implicit none

      write ( *, '(a)' ) ' '
      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'UNCMIN_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version.'
      write ( *, '(a)' ) '  Test the UNCMIN library.'

      call test01 ( )
      call test02 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'UNCMIN_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 demonstrates the use of the simple user interface OPTIF0.
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter (n=2)

      integer nr
      parameter (nr=n)

      double precision a(nr,n)
      double precision f
      external f1
      double precision g(n)
      integer i
      integer itrmcd
      double precision work(n,9)
      double precision x(n)
      double precision x0(n)

      write ( *, * ) ' '
      write ( *, * ) 'TEST01'
      write ( *, * ) '  Test OPTIF0, the simple interface to UNCMIN.'
c
c  Initial estimate of solution.
c
      x0(1) = 1.0D+00
      x0(2) = 1.0D+00

      call optif0 ( nr, n, x0, f1, x, f, g, itrmcd, a, work )

      write ( *, * ) ' '
      write ( *, * ) '  Output from UNCMIN:'
      write ( *, * ) '  Error code=', itrmcd
      call explain(itrmcd)
      write ( *, * ) '  F(X*) =', f
      write ( *, * ) '  X* =', (x(i), i = 1,n)

      write ( *, * ) ' '
      write ( *, * ) '  (Partial reference results for comparison:)'
      write ( *, * ) '      19.9145       -20.6011       -5.26250'
      write ( *, * ) '      19.9900       -20.6230        19.9145'
      write ( *, * ) '      20.0100       -20.6230        19.9145'
      write ( *, * ) '      19.9900       -20.6023        19.9145'
      write ( *, * ) '  Error code =           1'
      write ( *, * ) '  F(X*) =    91.0001'
      write ( *, * ) '  X* =    19.9900       -20.6230'

      return
      end
      subroutine f1 ( n, x, f )

c*********************************************************************72
c
cc F1 is the objective function for problem 1.
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 4 )

      integer n

      double precision b(m)
      double precision f
      integer j
      double precision r
      double precision t(m)
      double precision x(n)
c
c  The data to be fitted, (t(i), y(i)).
c
      t(1) =  0.0D+00
      t(2) =  1.0D+00
      t(3) =  2.0D+00
      t(4) =  3.0D+00
      b(1) = 20.0D+00
      b(2) =  9.0D+00
      b(3) =  3.0D+00
      b(4) =  1.0D+00

      f = 0.0D+00
      do j = 1, m
        r = b(j) - x(1) * exp ( x(2) * t(j) )
        f = f + r * r
      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 demonstrates the full user interface routine OPTIF9.
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nr
      parameter (nr=20)

      double precision a(nr,nr)
      external d2fcn
      double precision dlt
      external f2
      double precision fpls
      double precision fscale
      external g2
      double precision gpls(nr)
      double precision gradtl
      integer i
      integer iagflg
      integer iahflg
      integer iexp
      integer iprint
      integer itnlim
      integer itrmcd
      integer itry
      integer lounit
      integer method
      integer msg
      integer n
      integer ndigit
      double precision stepmx
      double precision steptl
      double precision typsiz(nr)
      double precision wrk(nr,8)
      double precision x(nr)
      double precision xpls(nr)

      do itry = 1, 3

        write(*,*)' '
        write(*,*)' '

        n = 3

        x(1) = -1.0D+00
        x(2) = 0.0D+00
        x(3) = 0.0D+00
c
c  Set the variables to default values.
c
        call dfault(n,typsiz,fscale,method,iexp,msg,ndigit,itnlim,
     &  iagflg,iahflg,iprint,lounit,dlt,gradtl,stepmx,steptl)
c
c  Change some values from their defaults.
c
        iprint = 1
        iexp = 0
        gradtl = 0.000001D+00
        iagflg = 1
        iahflg = 0
        itnlim = 25
        method = itry
        stepmx = 8.0D+00
        steptl = 0.0001D+00
        typsiz(1) = 1.0D+00
        typsiz(2) = 1.0D+00
        typsiz(3) = 1.0D+00
        fscale = 50.0D+00
c
c  Call OPTIF9, the full interface to UNCMIN.
c
        write ( *, * ) ' '
        write ( *, * ) 'uncprb  function of 3 variables.'
        write ( *, * ) ' '

        if ( method .eq. 1 ) then
          write ( *, * ) '  Solve problem using line search.'
        else if ( method .eq. 2 ) then
          write ( *, * ) '  Solve problem using double dogleg method.'
        else if ( method .eq. 3 ) then
          write ( *, * ) '  Solve problem using more-hebdon method.'
        end if

        call optif9(nr,n,x,f2,g2,d2fcn,typsiz,fscale,method,
     &  iexp,msg,ndigit,itnlim,iagflg,iahflg,iprint,lounit,dlt,gradtl,
     &  stepmx,steptl,xpls,fpls,gpls,itrmcd,a,wrk)
c
c  Print results.
c
        write ( *, * ) ' '
        write ( *, * ) '  termination code itrmcd=',itrmcd
        call explain(itrmcd)
        write ( *, * ) '  return code msg=',msg
        write ( *, * ) ' '
        write ( *, * ) '  x* =     ',(xpls(i),i =1,n)
        write ( *, * ) ' '
        write ( *, * ) '  f(x*) =  ',fpls
        write ( *, * ) ' '
        write ( *, * ) '  gradient=',(gpls(i),i=1,n)
        write ( *, * ) ' '

      end do

      return
      end
      subroutine f2 ( n, x, f )

c*********************************************************************72
c
cc F2 is the optimization function for TEST02.
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n

      double precision del
      double precision f
      double precision r1
      double precision r2
      double precision r3
      double precision x(n)

      if ( x(1) .lt. 0.0D+00 ) then
        del = 0.5D+00
      else
        del = 0.0D+00
      end if

      r1 = 10.0D+00 * ( x(3) - 10.0D+00 * atan2 ( x(2), x(1) ) + del )
      r2 = 10.0D+00 * ( sqrt ( x(1)**2 + x(2)**2 ) - 1.0D+00 )
      r3 = x(3)

      f = r1 * r1
     &  + r2 * r2
     &  + r3 * r3

      return
      end
      subroutine g2 ( n, x, g )

c*********************************************************************72
c
cc G2 is the gradient function for TEST02
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n

      double precision del
      double precision g(n)
      double precision r1
      double precision r2
      double precision r3
      double precision x(n)

      if ( x(1) .lt. 0.0D+00 ) then
        del = 0.5D+00
      else
        del = 0.0D+00
      end if

      r1 = 10.0D+00 * ( x(3) - 10.0D+00 * atan2 ( x(2), x(1) ) + del )
      r2 = 10.0D+00 * ( sqrt ( x(1)**2 + x(2)**2 ) - 1.0D+00 )
      r3 = x(3)

      g(1) = 2.0D+00 * r1 * ( 100.0D+00 * x(2)
     &  / ( x(1)**2 + x(2)**2 ) )
     &  + 20.0D+00 * r2 * ( x(1)/ sqrt(x(1)**2 + x(2)**2 ) )

      g(2) = 2.0D+00 * r1 * ( -100.0D+00 * x(1)
     &  / ( x(1)**2 + x(2)**2 ) )
     &    + 20.0D+00 * r2 * ( x(2 ) / sqrt ( x(1)**2 + x(2)**2 ) )

      g(3) = 20.0D+00 * r1 + 2.0D+00 * x(3)

      return
      end
