      program main

c*********************************************************************72
c
cc NL2SOL_PRB1 tests NL2SOL and NL2SNO on the Madsen example.
c
c  Modified:
c
c    03 April 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'NL2SOL_PRB1:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the NL2SOL library.'

      call test01 ( )
      call test02 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'NL2SOL_PRB1:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests NL2SOL on the Madsen example.
c
c  Modified:
c
c    27 March 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer meqn
      integer nvar

      parameter ( meqn = 3 )
      parameter ( nvar = 2 )

      integer iv(60+nvar)
      external madj
      external madr
      external ufparm
      integer uiparm(1)
      real urparm(1)
      real v(93+meqn*nvar+3*meqn+nvar*(3*nvar+33)/2)
      real x(nvar)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' ) '  Test the NL2SOL routine,'
      write ( *, '(a)' ) '  which requires a user residual '
      write ( *, '(a)' ) '  and jacobian routines.'
c
c  Set the initial solution estimate.
c
      x(1) = 3.0E+00
      x(2) = 1.0E+00

      iv(1) = 0

      call nl2sol ( meqn, nvar, x, madr, madj, iv, v, uiparm, 
     &  urparm, ufparm )

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests NL2SNO on the Madsen example.
c
c  Modified:
c
c    27 March 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer meqn
      integer nvar

      parameter ( meqn = 3 )
      parameter ( nvar = 2 )

      integer iv(60+nvar)
      external madr
      external ufparm
      integer uiparm(1)
      real urparm(1)
      real v(147)
      real x(nvar)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02:'
      write ( *, '(a)' ) '  Test the NL2SNO routine,'
      write ( *, '(a)' ) '  which requires only a user residual.'
      write ( *, '(a)' ) '  The jacobian is approximated internally.'
c
c  Set the initial solution estimate.
c
      x(1) = 3.0E+00
      x(2) = 1.0E+00

      iv(1) = 0

      call nl2sno ( meqn, nvar, x, madr, iv, v, uiparm, urparm, ufparm )

      return
      end
      subroutine madr ( meqn, nvar, x, nf, r, uiparm, urparm, ufparm )

c*********************************************************************72
c
cc MADR computes the Madsen residual.
c
c  Discussion:
c
c    Given the value of the vector X, this routine computes the 
c    value of F(X), the vector function whose norm we are trying
c    to minimize.
c
c  Modified:
c
c    27 March 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer MEQN, the number of functions.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, real X(NVAR), the current value of the variables.
c
c    Input, integer NF, the number of times the residual routine
c    has been called so far.
c
c    Output, real R(MEQN), the residual vector, that is, the
c    value of the functions for the given input value of the variables.
c
c    Input, integer UIPARM(*), a user array.
c
c    Input, real URPARM(*), a user array.
c
c    Input, external UFPARM, an external reference to a user subroutine
c    or function.
c
      implicit none

      integer meqn
      integer nvar

      integer nf
      real r(meqn)
      external ufparm
      integer uiparm(*)
      real urparm(*)
      real x(nvar)

      r(1) = x(1)**2 + x(2)**2 + x(1) * x(2)
      r(2) = sin ( x(1) )
      r(3) = cos ( x(2) )

      return
      end
      subroutine madj ( meqn, nvar, x, nf, jac, uiparm, urparm, ufparm )

c*********************************************************************72
c
cc MADJ computes the Madsen jacobian.
c
c  Discussion:
c
c    Given the value of the vector X, this routine computes the 
c    jacobian matrix of the vector function F(X), which has the
c    form
c
c      Jac(I,J) = d Fi / d Xj
c
c  Modified:
c
c    27 March 2006
c
c  Parameters:
c
c    Input, integer MEQN, the number of functions.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, real X(NVAR), the current value of the variables.
c
c    Input, integer NF, the number of times the residual routine
c    has been called so far.
c
c    Output, real JAC(MEQN,NVAR), the jacobian matrix.  JAC(I,J) is
c    the derivative of function I with respect to variable J.
c
c    Input, integer UIPARM(*), a user array.
c
c    Input, real URPARM(*), a user array.
c
c    Input, external UFPARM, an external reference to a user subroutine
c    or function.
c
      implicit none

      integer meqn
      integer nvar

      real jac(meqn,nvar)
      integer nf
      external ufparm
      integer uiparm(*)
      real urparm(*)
      real x(nvar)

      jac(1,1) = 2.0E+00 * x(1) + x(2)
      jac(1,2) = 2.0E+00 * x(2) + x(1)

      jac(2,1) = cos ( x(1) )
      jac(2,2) = 0.0E+00

      jac(3,1) = 0.0E+00
      jac(3,2) = -sin ( x(2) )

      return
      end
      subroutine ufparm ( meqn, nvar, x )

c*********************************************************************72
c
cc UFPARM is a user-supplied external routine.
c
c  Discussion:
c
c    The name of the routine, the argument list, and even whether
c    it is a function or subroutine, are left to the user.
c
c    NL2SOL simply passes the external reference from the calling
c    program through to the residual and jacobian routines.
c
c    If the user has no need for this facility, then a dummy
c    routine like this one may be used.
c
c  Modified:
c
c    27 March 2006
c
c  Parameters:
c
c    Input, integer MEQN, the number of functions.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, real X(NVAR), the current value of the variables.
c
      implicit none

      integer meqn
      integer nvar

      real x(nvar)

      return
      end
