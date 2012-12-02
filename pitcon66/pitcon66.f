      subroutine pitcon ( df, fpar, fx, ierror, ipar, iwork, liw, 
     &  nvar, rwork, lrw, xr, slname )

c*********************************************************************72
c
cc PITCON is the driver routine for the continuation algorithm.
c
c  Discussion:
c
c    This is version 6.6 of PITCON, the University of Pittsburgh
c    continuation algorithm.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 May 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      external checkw
      external corect
      external df
      external fx
      external limit
      external setstp
      intrinsic sign
      external slname
      external start
      external tanpar
      external target
      external update

      double precision one
      parameter (one=1.0D+00)

      double precision zero
      parameter (zero=0.0D+00)

      integer liw
      integer lrw
      integer nvar

      double precision dets
      double precision fpar(*)
      integer i
      integer ierror
      integer ifound
      integer ipar(*)
      integer iwork(liw)
      integer iwrite
      integer job
      integer ltc
      integer lwk
      integer lxc
      integer lxf
      double precision rwork(lrw)
      double precision xr(nvar)
c
c  1.  Set a few parameters.
c
      ierror=0
      iwrite=iwork(7)
      lxc=30
      lxf=29+nvar+1
      ltc=29+2*nvar+1
      lwk=29+3*nvar+1
      if(iwork(1).lt.(-10))iwork(1)=-2
      if(iwork(1).gt.5)iwork(1)=0
c
c  2.  Preparations.
c
c  Check entries of IWORK and RWORK, compute some constants.
c
      if(iwork(1).eq.0)then

        call checkw(ierror,iwork,liw,lrw,nvar,rwork)

        if(ierror.ne.0)then
          write(*,*)' '
          write(*,*)'PITCON - Fatal error!'
          write(*,*)'  An error was detected in the user input.'
          write(*,*)'  The program can not proceed.'
          stop
        end if

      end if
c
c  3.  Check user Jacobian routine versus finite difference calculation.
c
      if(iwork(1).lt.0)then

        job=3
        call slname(dets,fx,df,fpar,ierror,iwork(2),ipar,iwork,liw,
     &    job,nvar,rwork,lrw,xr,rwork(lwk))

        if(ierror.ne.0)then

          if(iwrite.gt.0)then
            write(*,*)' '
            write(*,*)'PITCON - Warning!'
            write(*,*)'  An error occurred during the jacobian check.'
          end if

        end if

        return
      end if
c
c  4.  Starting point check
c
c  On first call for a given problem, check that F(XR) is small 
c  enough so that the starting point may be considered to lie on 
c  the curve.
c
c  If this is not the case, call the corrector to try to enforce it.
c
      if(iwork(1).eq.0)then
c
c  See if the initial jacobian should be set up.
c
        if(iwork(4).eq.2)then

          call stajac(df,fpar,fx,ierror,ipar,iwork(2),iwork,liw,
     &      lrw,nvar,rwork,rwork(lwk),xr,slname)

          if(ierror.ne.0)then

            if(iwrite.gt.0)then
              write(*,*)' '
              write(*,*)'PITCON - Fatal error!'
              write(*,*)'  An error occurred during the'
              write(*,*)'  initial jacobian setup.'
              write(*,*)' '
              write(*,*)'  The program can not continue!'
              stop
            end if

            return

          end if

        end if
c
c  Check that the starting point satisfies the equations.
c
        call start(df,fpar,fx,ierror,ipar,iwork(2),iwork,liw,
     &    lrw,nvar,rwork,rwork(lwk),rwork(lxc),rwork(lxf),xr,slname)

        if(ierror.ne.0)then

          if(iwrite.gt.0)then
            write(*,*)' '
            write(*,*)'PITCON - Fatal error!'
            write(*,*)'  An error occurred during the starting point'
            write(*,*)'  check.'
            write(*,*)' '
            write(*,*)'  The starting point does not satisfy the'
            write(*,*)'  accuracy requirements, and PITCON could not'
            write(*,*)'  correct it.'
            write(*,*)' '
            write(*,*)'  The program can not continue!'
          end if

        end if

        do i=1,nvar
          rwork(ltc+i-1)=zero
        enddo
        rwork(ltc+iwork(2)-1)=one

        return
      end if
c
c  5.  Target point
c
c  If IWORK(5) is nonzero, target points are sought.  Check to see if 
c  target component IWORK(5), also called "IT", has value lying between 
c  XC(IT) and XF(IT).  If so, get linearly interpolated starting point, 
c  and use Newton's method to get target point.
c
      call target(df,fpar,fx,ierror,ifound,ipar,iwork,liw,lrw,
     &  nvar,rwork,slname,rwork(lwk),rwork(lxc),rwork(lxf),xr)

      if(ifound.eq.1)then
        if(ierror.ne.0)then
          ierror=9
          if(iwrite.gt.0)then
            write(*,*)' '
            write(*,*)'PITCON - Warning!'
            write(*,*)'  An error occurred during the'
            write(*,*)'  target point computation.'
            write(*,*)' '
            write(*,*)'  The target point returned does'
            write(*,*)'  not satisfy the accuracy requirements.'
            write(*,*)' '
            write(*,*)'  However, the code can continue.'
          end if
        end if

        return
      end if
c
c  6.  Tangent and local continuation parameter calculation.
c
c  Unless the tangent and limit point calculations were already 
c  performed, (because the loop was interrupted for a limit point 
c  calculation), set up and solve the equation for the tangent vector.
c
c  Force the tangent vector to be of unit length, and try to preserve
c  the "sense" or "direction" of the curve by forcing the IPL-th 
c  component of the tangent vector to agree in sign with the IPL-th 
c  component of the previous secant vector.  (On the first step, however, 
c  we have to use the user's input direction to choose the sign).
c
c  Set the local continuation parameter IPC.
c
c  If IWORK(3) is 0, the program is free to vary IPC from step to step.
c  In that case, IPC is normally set to the index of the component of
c  greatest magnitude in the tangent vector.  However, if a limit point
c  appears to be coming in that direction, the index of the second 
c  greatest magnitude component might be chosen instead.
c
      if(iwork(10).ne.3)then

        call tanpar(df,fpar,fx,ierror,ipar,iwork,liw,
     &    lrw,nvar,rwork,slname,rwork(ltc),rwork(lwk),rwork(lxc),
     &    rwork(lxf),xr)

        if(ierror.ne.0)then

          if(iwrite.gt.0)then
            write(*,*)' '
            write(*,*)'PITCON - Fatal error.'
            write(*,*)'  The computation failed while computing the'
            write(*,*)'  parameter and the tangent vector.'
            write(*,*)' '
            write(*,*)'  The program can not proceed!'
          end if

          return
        end if
c
c  7.  Limit point check.
c
c  Skip this section if IWORK(6)=0.
c
c  Otherwise, user has requested a search for limit points in a given
c  index by setting IWORK(6), also called "LIM", to a nonzero value.
c
c  Compare LIM-th components of previous and current tangent vectors.
c  If a sign difference occurs, we assume a limit point has been passed.
c  Attempt to compute a point XR between the previous and current points,
c  for which the LIM-th component of the tangent is zero.
c
c  This search will be guided by a rootfinder.  The scalar function
c  to be zeroed out is the LIM-th tangent vector component.
c
        if((iwork(6).ne.0).and.
     &    (iwork(1).ne.4).and.
     &    (iwork(10).eq.3).and.
     &    (sign(one,rwork(26)).ne.sign(one,rwork(27))))then
 
          call limit(df,fpar,fx,ierror,ipar,iwork,liw,lrw,nvar,rwork,
     &      slname,rwork(ltc),rwork(lwk),rwork(lxc),rwork(lxf),xr)
 
          if(ierror.ne.0)then

            if(iwrite.gt.0)then
              write(*,*)' '
              write(*,*)'PITCON - Warning!'
              write(*,*)'  An error occurred during the'
              write(*,*)'  limit point computation.'
              write(*,*)' '
              write(*,*)'  The computed limit point does not'
              write(*,*)'  satisfy the accuracy requirements.'
              write(*,*)' '
              write(*,*)'  However, the code can continue.'
            end if

          end if
          return
        end if
      end if
c
c  8.  Compute next predictor step length, HTAN.
c
      if(iwork(10).gt.1)then
        call setstp(iwork,liw,lrw,rwork)
      end if
c
c  9.  Continuation step
c
c  Our current data is the current point XC, its tangent vector TC, and
c  a steplength HTAN.  We predict the location of the next point on the
c  curve using the Euler approximation XR=XC+HTAN*TC.
c
c  Newton iteration is applied to this point, to force it to lie on the
c  curve.  In order to make the system square, an augmenting equation
c  is added to the system, specifying that XR(IPC)=XC(IPC)+HTAN*TC(IPC).
c  (The right hand side is a constant.)
c
c  If the Newton correction process fails, the stepsize is reduced and
c  prediction and correction retried.  Failure will most likely be
c  signaled by repeated step reductions, until the minimum allowable
c  stepsize is reached.  If this occurs, PITCON has failed, and cannot
c  proceed along the curve any more.
c
      call trystp(df,fpar,fx,ierror,ipar,iwork,liw,lrw,nvar,rwork,
     &  slname,rwork(ltc),rwork(lwk),rwork(lxf),xr)

      if(ierror.ne.0)then
        if(iwrite.gt.0)then
          write(*,*)' '
          write(*,*)'PITCON - Fatal error.'
          write(*,*)'The computation failed while trying'
          write(*,*)'to compute the next point.'
          write(*,*)' '
          write(*,*)'The program can not proceed!'
        end if
        return
      end if
c
c  10.  Successful step.  Update information.
c
      call update(iwork,liw,lrw,nvar,rwork,rwork(ltc),rwork(lwk),
     &  rwork(lxc),rwork(lxf),xr)
c
c  Compute the convergence "quality", a factor between 1/8 and 8, which
c  tries to estimate how deeply inside the Newton attraction region we
c  are, and hence how much bolder or more timid we could be on the next
c  prediction.
c
      call coqual(iwork,liw,rwork,lrw)

      if(iwrite.ge.3)then
        write(*,*)
     &    'PITCON - Corrector convergence quality factor=',rwork(23)
      end if

      return
      end
      subroutine banjac ( eps, fcol, fpar, fprime, frow, fx, ierror,
     &  ipar, ipc, iwork, jac, liw, nband, neqn, nvar, x, xtemp, 
     &  work1, work2 )

c*********************************************************************72
c
cc BANJAC estimates a banded jacobian matrix.
c
c  Discussion:
c
c    This routine estimates the jacobian matrix FPRIME of the function FX,
c    using forward or central finite differences.  BANJAC is called by
c    BANSLV when the user has specified the jacobian option as 1 or 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c  EPS    Input, REAL EPS, a tolerance used for shifting the X values. 
c         A value of the square root of the machine precision is 
c         usually appropriate.
c
c  FCOL   Output, REAL FCOL(NEQN), the last column of the approximate
c         jacobian, which is allowed to be "full".  This comprises
c         matrix entries FPRIME(1,NVAR) through FPRIME(NEQN,NVAR).
c
c  FPAR   Input/output, REAL FPAR(*), user parameter vector, not
c         touched by this routine, but passed on to user routines.
c
c  FPRIME Output, REAL FPRIME(NBAND,NEQN), is the array into which the
c         the banded portion of the computed jacobian will be stored.
c         The LAPACK general band format is used, assigning entry (I,J)
c         to FPRIME(I-J+ML+MU+1,J), where ML and MU are the lower and
c         upper half bandwidths respectively.
c
c  FROW   Output, REAL FROW(NVAR), storage for the last (augmenting) row
c         of the jacobian, which will be all zero except for a 1 in
c         location IPC.
c
c  FX     Input, EXTERNAL FX, the name of the routine which evaluates the
c         function.
c
c         FX computes the value of the nonlinear function.  This name 
c         must be declared EXTERNAL in the calling program.  FX should 
c         evaluate the NVAR-1 function components at the input point X, 
c         and store the result in the vector FVEC.  An augmenting 
c         equation will be stored in entry NVAR of FVEC by the PITCON 
c         program.
c
c         FX should have the following form:
c
c           subroutine fx(nvar,fpar,ipar,x,fvec,ierror)
c
c           NVAR   Input, INTEGER NVAR, number of variables.
c
c           FPAR   Input/output, REAL FPAR(*), array of user parameters.
c
c           IPAR   Input/output, INTEGER IPAR(*), array of user parameters.
c
c           X      Input, REAL X(NVAR), the point at which function evaluation
c                  is required.
c
c           FVEC   Output, REAL FVEC(NVAR), the value of the function at point
c                  X.  Only the first NVAR-1 entries of FVEC are to be set by
c                  the routine.  PITCON sets the final value itself.
c
c           IERROR Output, INTEGER IERROR, error flag.  FX should set this to 0
c                  if there are no problems, or to 2 if there is a problem.
c
c  IERROR Output, INTEGER IERROR, error flag.  A nonzero value means that
c         there was an error in the user routine FX, or in BANJAC itself.
c         In either case, the jacobian has not been computed.
c
c  IPAR   Input, INTEGER IPAR(*), a user parameter vector passed to FX.
c         However, because this is a problem with a banded jacobian, entries
c         IPAR(1) and IPAR(2) are read by this routine.  IPAR(1) contains
c         ML, the lower half bandwidth of the jacobian, and IPAR(2) contains
c         MU, the upper half bandwidth of the jacobian.
c
c  IPC    Input, INTEGER IPC, the index of the current continuation parameter,
c         which is needed to determine the form of FROW.
c
c  IWORK  Input, INTEGER IWORK(LIW), work and statistics vector.  Only
c         required here so that we can count the number of function
c         evaluations.
c
c  JAC    Input, INTEGER JAC, the user requested jacobian option.  For
c         our purposes, the only two values of interest are:
c
c           1 = estimate jacobian with forward differences,
c           2 = estimate jacobian with central differences (twice the work)
c
c  LIW    Input, INTEGER LIW, the dimension of IWORK.
c
c  NBAND  Input, INTEGER NBAND, the first dimension of the jacobian matrix
c         FPRIME, NBAND=ML+MU+1.
c
c  NEQN   Input, INTEGER NEQN, the number of equations, equal to NVAR-1.
c
c  NVAR   Input, INTEGER NVAR, the number of variables.
c
c  X      Input, REAL X(NVAR), the point at which the jacobian is desired.
c
c  XTEMP,
c  WORK1,
c  WORK2  Work arrays, REAL XTEMP(NVAR), WORK1(NVAR), WORK2(NVAR).
c
      implicit none

      double precision one
      double precision two

      parameter (one=1.0D+00)
      parameter (two=2.0D+00)

      external fx
      external daxpy
      external dcopy
      external dscal

      intrinsic abs
      intrinsic max
      intrinsic min

      integer liw
      integer nband
      integer neqn
      integer nvar

      double precision eps
      double precision fcol(neqn)
      double precision fpar(*)
      double precision fprime(nband,neqn)
      double precision frow(nvar)
      integer iband
      integer ierror
      integer ihi
      integer ilo
      integer ipar(2)
      integer ipc
      integer irow
      integer iwork(liw)
      integer iwrite
      integer j
      integer jac
      integer kcall
      integer mband
      integer ml
      integer mu
      double precision skale
      double precision x(nvar)
      double precision xjac
      double precision xtemp(nvar)
      double precision work1(nvar)
      double precision work2(nvar)

      iwrite=iwork(7)
      ml=ipar(1)
      mu=ipar(2)
      mband=ml+mu+1

      if(jac.eq.1)then

        call fx(nvar,fpar,ipar,x,work2,ierror)
        iwork(22)=iwork(22)+1

        if(ierror.ne.0)then

          if(iwrite.gt.0)then
            write(*,*)'BANJAC - Fatal error!'
            write(*,*)'  The user function returns IERROR=',ierror
          end if

          return
        end if
      end if

      xjac=one
      if(jac.eq.2)xjac=two

      do kcall=1,mband

        call dcopy(nvar,x,1,xtemp,1)

        do j=kcall,neqn,mband
          xtemp(j)=x(j)+eps*(one+abs(x(j)))
        enddo

        call fx(nvar,fpar,ipar,xtemp,work1,ierror)
        iwork(22)=iwork(22)+1

        if(ierror.ne.0)then
          if(iwrite.gt.0)then
            write(*,*)'BANJAC - Fatal error!'
            write(*,*)'  User function returns IERROR=',ierror
          end if
 
          return
        end if

        if(jac.eq.2)then
          call dcopy(nvar,x,1,xtemp,1)

          do j=kcall,neqn,mband
            xtemp(j)=x(j)-eps*(one+abs(x(j)))
          enddo

          call fx(nvar,fpar,ipar,xtemp,work2,ierror)
          iwork(22)=iwork(22)+1

          if(ierror.ne.0)then

            if(iwrite.gt.0)then
              write(*,*)'BANJAC - Fatal error!'
              write(*,*)'  User function returns IERROR=',ierror
            end if

            return
          end if

        end if

        do j=kcall,neqn,mband
          ilo=max(1,j-mu)
          ihi=min(neqn,j+ml)
          irow=ilo-j+ml+mu+1
          iband=ihi-ilo+1
          call daxpy(iband,-one,work2(ilo),1,work1(ilo),1)
          skale=one/(xjac*eps*(one+abs(x(j))))
          call dscal(iband,skale,work1(ilo),1)
          call daxpy(iband,one,work1(ilo),1,fprime(irow,j),1)
        enddo

      enddo
c
c  Compute last column of jacobian, rows 1 to NEQN
c
      call dcopy(nvar,x,1,xtemp,1)
      xtemp(nvar)=x(nvar)+eps*(one+abs(x(nvar)))
      call fx(nvar,fpar,ipar,xtemp,work1,ierror)
      iwork(22)=iwork(22)+1

      if(ierror.ne.0)then
        if(iwrite.gt.0)then
          write(*,*)'BANJAC - Fatal error!'
          write(*,*)'  User function returns IERROR=',ierror
          end if
        return
      end if

      if(jac.eq.2)then

        xtemp(nvar)=x(nvar)-eps*(one+abs(x(nvar)))
        call fx(nvar,fpar,ipar,xtemp,work2,ierror)
        iwork(22)=iwork(22)+1

        if(ierror.ne.0)then
          if(iwrite.gt.0)then
            write(*,*)'BANJAC - Fatal error!'
            write(*,*)'  User function returns IERROR=',ierror
          end if
          return
        end if

      end if

      call daxpy(neqn,-one,work2,1,work1,1)
      skale=one/(xjac*eps*(one+abs(x(nvar))))
      call dscal(neqn,skale,work1,1)
      call daxpy(neqn,one,work1,1,fcol,1)
c
c  Do last row, J=1,NVAR
c
      frow(ipc)=frow(ipc)+one

      return
      end
      subroutine banslv ( dets, fx, df, fpar, ierror, ipc, ipar, 
     &  iwork, liw, job, nvar, rwork, lrw, x, y )

c*********************************************************************72
c
cc BANSLV solves a banded linear system.
c
c  Discussion:
c
c    The linear system has the form
c
c          ( DFDY | DFDZ )
c    A*X = (-------------) * X = B
c          (   E(IPC)    )
c
c  where
c
c    B is a given vector of length NVAR,
c    DFDY is "logically" an NVAR-1 by NVAR-1 matrix, with band 
c      structure,
c    DFDZ is an NVAR-1 vector, 
c    E(IPC) is an NVAR vector whose only nonzero entry is a 1 in
c      position IPC.
c
c
c  DFDY is actually stored compactly, using LAPACK general band storage,
c
c
c  DFDY and DFDZ represent the jacobian of an NVAR-1 dimensional function
c  of NVAR variables, and E(IPC) is the augmenting row, determined by
c  the choice of "continuation parameter" IPC.
c
c
c    BANSLV factors and solves the linear system A*x=b, taking advantage
c    of the bandedness of the DFDY subsystem, which would be lost if the 
c    full system was factored and solved directly.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c  DETS   Output, REAL DETS, the sign of the determinant of the full 
c         matrix A.
c
c  FX     Input, EXTERNAL FX, the name of the routine which evaluates the
c         function.
c
c         FX computes the value of the nonlinear function.  This name must be
c         declared EXTERNAL in the calling program.  FX should evaluate the
c         NVAR-1 function components at the input point X, and store the result
c         in the vector FVEC.  An augmenting equation will be stored in entry
c         NVAR of FVEC by the PITCON program.
c
c         FX should have the following form:
c
c           subroutine fx(nvar,fpar,ipar,x,fvec,ierror)
c
c           NVAR   Input, INTEGER NVAR, number of variables.
c
c           FPAR   Input/output, REAL FPAR(*), user parameters.
c
c           IPAR   Input/output, INTEGER IPAR(*), user parameters.
c
c           X      Input, REAL X(NVAR), the point at which function
c                  evaluation is required.
c
c           FVEC   Output, REAL FVEC(NVAR), value of the function at
c                  X.  Only the first NVAR-1 entries of FVEC are to be set by
c                  the routine.  PITCON sets the final value itself.
c
c           IERROR Output, INTEGER IERROR, error flag.  FX should set this to 0
c                  if there are no problems, or to 2 if there is a problem.
c
c  DF     Input, EXTERNAL FP, the name of the user supplied routine which
c         computes the jacobian (DFDY, DFDZ) of the NVAR-1 nonlinear 
c         functions.  
c
c         Jacobian entries for the NVAR-th, augmenting, function are
c         inserted by this routine.  
c
c         Because of the special banded storage arrangement, great care
c         must be taken to store the information properly.
c
c         Let MU and ML be the upper and lower bandwidths of DFDY.
c         Set NBAND=2*ML+MU+1.
c
c         Then the index K of RWORK into which we want to store entry
c         (I,J) of DFDY is determined as follows:
c
c           K=I+J*(NBAND-1)-ML
c
c         Here we are assuming that I-ML <= J <= I+MU.  For any other values 
c         of I and J, the corresponding entry of DFDY is assumed to be zero.
c
c         There are NVAR-1 entries in DFDZ, and the I-th entry of this
c         vector is stored in RWORK(K), where
c
c           K=(NVAR-1)*NBAND+I
c
c         DF must be a subroutine of the form:
c
c           subroutine df(nvar,fpar,ipar,x,rwork,ierror)
c
c           integer nvar
c
c           real fpar(*)
c           integer ierror
c           integer ipar(*)
c           real rwork(*)
c           real x(nvar)
c
c           ml=ipar(1)
c           mu=ipar(2)
c           nband=2*ml+mu+1
c
c           do i=1,nvar-1
c             do j=max(i-ml,1), min(i+mu,nvar-1)
c               k=i+j*(nband-1)-ml
c               rwork(k)=d f(i)/d x(j)
c             enddo
c           enddo
c
c           do i=1,nvar-1
c             k=(nvar-1)*nband+i
c             rwork(k)=d f(i)/d x(nvar)
c           enddo
c
c           return
c           end
c
c         If a serious error occurs during the execution of DF, and
c         you wish to signal that the continuation code should halt,
c         set IERROR nonzero before return.
c
c  FPAR   Input/output, REAL FPAR(*), a user defined parameter array.
c
c  IERROR Output, INTEGER IERROR, the error return flag.
c
c         0, No errors were detected.
c
c         1, data or storage error, including illegal values for NVAR,
c            IPC, MU, ML, or insufficient storage in RWORK or IWORK.
c
c         2, The user set a nonzero error return in DF.
c
c         3, The matrix DFDY or A is numerically singular.
c            In some cases, a different choice of IPC could rectify
c            this problem.
c
c  IPC    Input, INTEGER IPC, the continuation parameter.
c
c         IPC determines the form of the final, augmenting row of the
c         jacobian matrix.  
c
c         If IPC is equal to NVAR, then the linear system can easily
c         be solved using a standard band matrix solver, once we
c         have solved for X(NVAR).
c
c         But in the general case when IPC is not NVAR, we have to
c         do some work to modify the system, and still be able to
c         use a band solver.
c
c  IPAR   Input, INTEGER IPAR(*).  
c
c         IPAR(1)=ML, the lower bandwidth of DFDY,
c         IPAR(2)=MU, the upper bandwidth of DFDY.
c
c         The other entries of IPAR are not referenced by PITCON, and
c         may be used to pass information between the calling program
c         and the user routines FX and FP.
c
c  IWORK  Workspace, IWORK(*), a work array used by the continuation code
c         to store statistics, pointers, and the pivot vector for the
c         linear equation solver.
c
c         IWORK(13) stores the address of the first entry in IWORK
c         used for pivoting.
c
c         IWORK(15) stores the address of the first entry in RWORK
c         used to store the Jacobian.
c 
c         IWORK(20) counts the number of matrix factorizations.
c
c         IWORK(21) counts the number of linear system back-solves.
c
c  LIW    Input, INTEGER LIW, the dimension of IWORK.
c
c  JOB    Input, INTEGER JOB, controls the action of the routine:
c
c         0, evaluate jacobian, decompose jacobian, compute determinant,
c         and solve linear system.
c
c         1, solve a system with a new right hand side, and a previously
c         factored jacobian.
c
c         2, evaluate jacobian, decompose jacobian, and compute determinant.
c
c         3, Check jacobian matrix.  Call user jacobian routine,
c         multiply by -1.0, add finite difference jacobian,
c         print largest entry.
c
c  NVAR   Input, INTEGER NVAR, the number of variables.  The dimension of X.
c
c  RWORK  Workspace, REAL RWORK(LRW).  RWORK contains workspace used
c         for storing the jacobian, as well as various vectors and
c         scalars.  The IWORK array contains pointers to the beginning
c         locations of some of these objects.
c
c  LRW    Input, INTEGER LRW, the dimension of RWORK.
c
c  X      Input, REAL X(NVAR), the point at which the jacobian is to
c         be evaluated.
c
c  Y      Input/output, REAL Y(NVAR).  Right hand side/solution vector.
c
c         On input with JOB=1, Y contains the right hand side of a linear
c         system to be solved, and on output, contains the solution to that
c         same linear system.
c
      implicit none

      double precision one
      double precision zero

      parameter (one=1.0D+00)
      parameter (zero=0.0D+00)

      external banjac
      external df
      external fx
      external idamax
      intrinsic max
      intrinsic min
      external daxpy
      external ddot
      external dgbdet
      external dgbtrf
      external dgbtrs
      intrinsic sqrt
      external dscal
      external dswap

      integer liw
      integer lrw
      integer nvar

      double precision ak
      double precision det
      double precision dets
      double precision fpar(*)
      integer i
      integer ierror
      integer info
      integer ipar(*)
      integer ipc
      integer irl
      integer iru
      integer idamax
      integer itemp
      integer iwork(liw)
      integer iwrite
      integer j
      integer jac
      integer jack
      integer job
      integer jtemp
      integer k
      integer lda
      integer ldfl
      integer ldfx
      integer ldx
      integer lfxm
      integer lfxp
      integer lilst
      integer lpiv
      integer lrlst
      integer lrowip
      integer mband
      integer ml
      integer mu
      integer nband
      integer ndim
      integer neqn
      integer nswap
      double precision rwork(lrw)
      double precision ddot
      double precision skale
      double precision temp
      double precision x(nvar)
      double precision y(nvar)

      ierror=0
      iwrite=iwork(7)
      ml=ipar(1)
      mu=ipar(2)

      neqn=nvar-1

      if(ml.lt.0)then
        ierror=1
        write(*,*)'BANSLV - Fatal error!'
        write(*,*)'  Illegal lower bandwidth ML=',ml
        write(*,*)'  ML must be at least 0.'
        stop
      elseif(ml.gt.neqn-1)then
        ierror=1
        write(*,*)'BANSLV - Fatal error!'
        write(*,*)'  Illegal lower bandwidth ML=',ml
        write(*,*)'  ML must be no more than ',neqn-1
        stop
      elseif(mu.lt.0)then
        ierror=1
        write(*,*)'BANSLV - Fatal error!'
        write(*,*)'  Illegal upper bandwidth MU=',mu
        write(*,*)'  MU must be at least 0.'
        stop
      elseif(mu.gt.neqn-1)then
        ierror=1
        write(*,*)'BANSLV - Fatal error!'
        write(*,*)'  Illegal upper bandwidth MU=',mu
        write(*,*)'  MU must be no more than ',neqn-1
        stop
      end if

      lda=2*ml+mu+1

      mband=ml+mu+1
      nband=mband+ml
      lpiv=iwork(13)
      lilst=lpiv+neqn-1
      ldfx=iwork(15)
      ldfl=ldfx+nband*neqn
      lrowip=ldfl+nvar
      jac=iwork(9)
c
c  Make sure that the jacobian routine is available if the user
c  has specified an option that requires it!
c
      if(jac.ne.0.and.
     &  (iwork(1).eq.(-1).or.iwork(1).eq.(-2).or.
     &   iwork(1).eq.(-3).or.iwork(1).eq.(-4).or.
     &   iwork(1).eq.(-9).or.iwork(1).eq.(-10)))then
        ierror=4

        if(iwrite.ge.1)then
          write(*,*)'BANSLV - Fatal error!'
          write(*,*)'  An option was selected requiring'
          write(*,*)'  a user jacobian, but the value of'
          write(*,*)'  IWORK(9) indicates no such routine'
          write(*,*)'  is available!'
        end if

        return
      end if
c
c  Check that enough storage is available.
c
      if(jac.eq.0.and.job.ne.3)then
        lrlst=lrowip+nvar-1
      else
        lfxp=lrowip+nvar
        lfxm=lfxp+nvar
        ldx=lfxm+nvar
        lrlst=ldx+nvar-1
      end if

      ndim=neqn*nband+nvar+nvar-1

      if(lilst.gt.liw)then
        ierror=1

        if(iwrite.ge.1)then
          write(*,*)'BANSLV - Fatal error!'
          write(*,*)'  Insufficient integer workspace in IWORK!'
          write(*,*)'  Need workspace size    =',lilst
          write(*,*)'  Available workspace LIW=',liw
        end if

        return
      end if
 
      if(lrlst.gt.lrw)then
        ierror=1
        if(iwrite.ge.1)then
          write(*,*)'BANSLV - Fatal error!'
          write(*,*)'  Insufficient real workspace in RWORK!'
          write(*,*)'  Needed workspace size    =',lrlst
          write(*,*)'  Available workspace LRW=  ',lrw
        end if
        return
      end if
c
c  If JOB=1, we need to solve a linear system.
c  Two entirely different methods are used to solve a linear system,
c  depending on whether the parameter is the last variable or not.
c
      if(job.eq.1)then
        if(ipc.eq.nvar)then
          go to 70
        else
          go to 40
        end if
      end if
c
c  JOB=0 OR 2 means we must 
c    evaluate the jacobian, 
c    factor it,
c    get the sign of the determinant.
c
      do i=1,ndim
        rwork(ldfx+i-1)=zero
      enddo
c
c  If a user jacobian routine is available, invoke it.
c
      if(jac.eq.0)then

        if(iwork(1).gt.(-5).or.iwork(1).lt.(-8))then
          call df(nvar,fpar,ipar,x,rwork(ldfx),ierror)
          iwork(19)=iwork(19)+1

          if(ierror.ne.0)then

            if(iwrite.gt.0)then
              write(*,*)'BANSLV - Fatal error!'
              write(*,*)'  The user jacobian returned IERROR=',ierror
            end if
            
            return
          end if

          rwork(lrowip-1+ipc)=one
        end if

      end if
c
c  If we're going to compare the user jacobian to an approximate 
c  one, negate the user jacobian.
c
      if(job.eq.3)then
        if(iwork(1).le.(-1).and.iwork(1).ge.(-4))then
          call dscal(ndim,-one,rwork(ldfx),1)
        end if
      end if
c
c  If we don't have a user jacobian, or we need to compare the
c  user jacobian to an approximate one, get the approximation.
c
      if((iwork(1).ge.0.and.(jac.eq.1.or.jac.eq.2)).or.
     &   (iwork(1).le.(-1).and.iwork(1).ge.(-8)))then

        jack=jac

        if(job.eq.3)then
          if(mod(iabs(iwork(1)),2).eq.1)then
            jack=1
          else
            jack=2
          end if
        end if

        call banjac(rwork(18),rwork(ldfl),fpar,rwork(ldfx),
     &    rwork(lrowip),fx,ierror,ipar,ipc,iwork,jack,liw,nband,
     &    neqn,nvar,x,rwork(ldx),rwork(lfxp),rwork(lfxm))

        if(ierror.ne.0)then
          if(iwrite.gt.0)then
            write(*,*)'BANSLV - Fatal error!'
            write(*,*)'  BANJAC returns IERROR==',ierror
          end if
          return
        end if

      end if
c
c  If JOB=3, print out information:
c    IWORK(1)=-1, -2, -5, -6, -9, print largest entry of matrix.
c    IWORK(1)=-3, -4, -7, -8, -10, also print out entire matrix.
c
      if(job.eq.3)then
        k=idamax(ndim,rwork(ldfx),1)
        ak=rwork(ldfx+k-1)

        if(k.le.neqn*nband)then
          j=((k-1)/nband)+1
          i=k-(j-1)*nband+j-ml-mu-1
        elseif(k.le.neqn*nband+neqn)then
          i=k-neqn*nband
          j=nvar
        else
          i=nvar
          j=k-neqn*nband-neqn
        end if

        if(iwork(1).le.(-1).and.iwork(1).ge.(-4))then
          write(*,*)' '
          write(*,*)
     &      'BANSLV - Maximum value of FP_Approx(I,J)-FP_User(I,J)'
        elseif(iwork(1).le.(-5).and.iwork(1).ge.(-8))then
          write(*,*)' '
          write(*,*)
     &      'BANSLV - Maximum value of finite difference jacobian:'
        elseif(iwork(1).le.(-9).and.iwork(1).ge.(-10))then
          write(*,*)' '
          write(*,*)
     &      'BANSLV - Maximum value of user supplied jacobian:'
        end if
 
        write(*,*)ak,' I, J=',i,j
        write(*,*)' '

        if(iwork(1).eq.(-3).or.iwork(1).eq.(-4).or.
     &     iwork(1).eq.(-7).or.iwork(1).eq.(-8).or.
     &     iwork(1).eq.(-10))then

          if(iwork(1).eq.(-3).or.iwork(1).eq.(-4))then
            write(*,*)' '
            write(*,*)'BANSLV - Entire difference matrix:'
            write(*,*)'FP_Approx(I,J)-FP_User(I,J)'
            write(*,*)' '
          elseif(iwork(1).eq.(-7).or.iwork(1).eq.(-8))then
            write(*,*)' '
            write(*,*)'BANSLV - Finite difference jacobian:'
            write(*,*)' '
          elseif(iwork(1).eq.(-10))then
            write(*,*)' '
            write(*,*)'BANSLV - User supplied jacobian:'
            write(*,*)' '
          end if

          do i=1,nvar
            do j=1,nvar
              if(j.eq.nvar)then
                k=ldfx-1+(nvar-1)*nband+i
                write(*,*)rwork(k),' I, J=',i,j
              elseif((j-i.le.ml).and.(i-j.le.mu))then
                k=ldfx-1+i+j*(nband-1)-ml
                write(*,*)rwork(k),' I, J=',i,j
              end if
            enddo
            write(*,*)' '
          enddo

        end if
        return
      end if

      if(ipc.eq.nvar)go to 60
c
c  Switch the NVAR-th and IPC-th rows.
c
      irl=max(1,ipc-ml)
      iru=min(neqn,ipc+mu)
      nswap=iru+1-irl
      itemp=ldfx-ml-1+ipc+irl*(nband-1)
      jtemp=nband-1
      call dswap(nswap,rwork(lrowip-1+irl),1,rwork(itemp),jtemp)
c
c  Decompose the submatrix and obtain determinant sign.
c
      call dgbtrf(neqn,neqn,ml,mu,rwork(ldfx),lda,iwork(lpiv),info)

      iwork(20)=iwork(20)+1
      if(info.ne.0)then

        if(iwrite.ge.1)then
          write(*,*)'BANSLV - Fatal error!'
          write(*,*)'  LU factor routine DGBTRF returns INFO=',info
          write(*,*)'  Current index IPC= ',ipc
        end if

        ierror=3
        return
      end if

      call dgbdet(neqn,ml,mu,rwork(ldfx),lda,iwork(lpiv),det)
      dets=zero

      if(det.gt.zero)then
        dets=-one
      elseif(det.lt.zero)then
        dets=one
      end if
c
c  Set the right hand side of the auxilliary system to the last
c  column of the jacobian, minus E(IPC).
c
c  Shuffle the IPC-th and NVAR-th entries of this right hand
c  side, to reflect the pivoting of the equations.
c
c  Then solve the system.
c
      temp=rwork(ldfl+ipc-1)-one
      rwork(ldfl+ipc-1)=rwork(ldfl+nvar-1)
      rwork(ldfl+nvar-1)=temp

      call dgbtrs('n',neqn,ml,mu,1,rwork(ldfx),lda,iwork(lpiv),
     &  rwork(ldfl),neqn,info)

      iwork(21)=iwork(21)+1

      if(info.ne.0)then

        if(iwrite.ge.1)then
          write(*,*)'BANSLV - Fatal error!'
          write(*,*)'  LU backsolve routine DGBTRS returns INFO=',info
        end if

        ierror=1
        return
      end if
c
c  Solve for last entry of auxilliary solution.
c
      rwork(ldfl+nvar-1)=rwork(ldfl+nvar-1)
     & -ddot(neqn,rwork(lrowip),1,rwork(ldfl),1)
c
c  Adjust the sign of the determinant.
c
      if(one+rwork(ldfl+nvar-1).eq.zero)then
        ierror=3

        if(iwrite.ge.1)then
          write(*,*)'BANSLV - Algorithm fails, DENOM=0.0'
        end if

        return
      end if

      if(one+rwork(ldfl+nvar-1).lt.zero)dets=-dets
      if(job.eq.2)return
c
c  Solve the system.
c
40    continue
      if(ipc.eq.nvar)go to 70
c
c  Modify right hand side of main system.
c
      temp=y(ipc)
      y(ipc)=y(nvar)
      y(nvar)=temp
c
c  Solve subsystem.
c
      call dgbtrs('n',neqn,ml,mu,1,rwork(ldfx),lda,iwork(lpiv),
     &  y,neqn,info)
      iwork(21)=iwork(21)+1

      if(info.ne.0)then

        if(iwrite.ge.1)then
          write(*,*)'BANSLV - Fatal error!'
          write(*,*)'  DGBTRS returns INFO=',info
        end if

        ierror=1
        return
      end if
c
c  Solve for last entry of main solution.
c
      y(nvar)=y(nvar)-ddot(neqn,rwork(lrowip),1,y,1)
c
c  Correct the main solution, using a multiple of the "subsolution".
c
      skale=-y(nvar)/(one+rwork(ldfl+nvar-1))
      call daxpy(nvar,skale,rwork(ldfl),1,y,1)
      return
c
c  Factor the matrix for the special case of IPC=NVAR.
c
   60 continue

      call dgbtrf(neqn,neqn,ml,mu,rwork(ldfx),lda,iwork(lpiv),info)
      iwork(20)=iwork(20)+1

      if(info.ne.0)then

        if(iwrite.ge.1)then
          write(*,*)'BANSLV - Fatal error!!'
          write(*,*)'  DGBTRF returns INFO=',info
          write(*,*)'  Current index IPC= ',ipc
        end if

        ierror=3
        return

      end if

      call dgbdet(neqn,ml,mu,rwork(ldfx),lda,iwork(lpiv),det)

      dets=zero

      if(det.lt.zero)then
        dets=-one
      elseif(det.gt.zero)then
        dets=one
      end if

      if(job.eq.2)return
c
c  Solve the linear system for the special case where IPC=NVAR.
c
   70 continue

      call daxpy(neqn,-y(nvar),rwork(ldfl),1,y,1)

      call dgbtrs('n',neqn,ml,mu,1,rwork(ldfx),lda,iwork(lpiv),
     &  y,neqn,info)

      iwork(21)=iwork(21)+1

      return
      end
      subroutine checkw ( ierror, iwork, liw, lrw, nvar, rwork )

c*********************************************************************72
c
cc CHECKW checks the input arrays on the first call.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 May 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision hundrd
      double precision one
      double precision three
      double precision two
      double precision zero

      parameter (hundrd=100.0D+00)
      parameter (one=1.0D+00)
      parameter (three=3.0D+00)
      parameter (two=2.0D+00)
      parameter (zero=0.0D+00)

      integer liw
      integer lrw

      integer i
      integer ierror
      integer iwork(liw)
      integer iwrite
      integer nvar
      double precision r8_epsilon
      double precision rwork(lrw)
      double precision tcos
      double precision temp
      double precision tsin

      iwrite=iwork(7)
 
      if(iwrite.ge.1)then
        write(*,*)' '
        write(*,*)'PITCON 6.6'
        write(*,*)'University of Pittsburgh continuation code'
        write(*,*)' '
        write(*,*)'Last modified on 09 September 1994'
        write(*,*)'This version uses LAPACK for linear algebra.'
        write(*,*)'This version uses double precision arithmetic.'
        write(*,*)' '
      end if
 
      do i=11,17
        rwork(i)=zero
      enddo

      do i=21,29
        rwork(i)=zero
      enddo

      iwork(10)=0
      iwork(11)=0
      iwork(12)=0

      do i=18,28
        iwork(i)=0
      enddo

      rwork(8) = r8_epsilon ( )

      if(iwrite.ge.2)then
        write(*,*)'CHECKW - Machine epsilon=',rwork(8)
      end if
c
c  Set the value of the parameter used in estimating the jacobian
c  via finite differences.
c
      if(rwork(18).eq.0.0D+00)rwork(18)=sqrt(sqrt(rwork(8)))
      if(iwrite.ge.3)write(*,*)
     &  'CHECKW - Jacobian finite difference increment=',rwork(18)
c
c  Set the value of an angle which is "almost" a zero angle.
c
      tcos=sqrt(one-rwork(8))
      tsin=sqrt(rwork(8))
      rwork(10)=two*atan2(tsin,tcos)

      if(nvar.le.1)then
        ierror=1

        if(iwrite.ge.1)then
          write(*,*)'CHECKW - Fatal error!'
          write(*,*)'  The number of variables, NVAR, must be at '
     &      //'least 2.'
          write(*,*)'  The input value is ',nvar
        end if

        return
      end if
c
c  Set entries of IWORK which point to the next free location, and
c  total available space in IWORK and RWORK.
c
      iwork(13)=30
      iwork(14)=liw
      iwork(15)=29+4*nvar+1
      iwork(16)=lrw
c
c  Check allocated sizes of IWORK and RWORK versus minimal requirements.
c  The actual need for RWORK will generally be even greater than is checked
c  for here, depending on the solver and storage method chosen.
c
      if(liw.lt.iwork(13))then

        if(iwrite.ge.1)then
          write(*,*)'CHECKW - Fatal error!'
          write(*,*)'  LIW is too small!  Input value of LIW=',liw
          write(*,*)'  The minimum acceptable value=',iwork(13)
        end if

        ierror=1
        return
      end if
 
      if(lrw.lt.iwork(15))then
        if(iwrite.ge.1)then
          write(*,*)'CHECKW - Fatal error!'
          write(*,*)'  LRW is too small!  Input value of LRW=',lrw
          write(*,*)'  The minimum acceptable value=',iwork(15)
        end if
        ierror=1
        return
      end if
c
c  Check entries of IWORK
c
      if(iwork(2).lt.1.or.iwork(2).gt.nvar)then
        iwork(2)=nvar
        if(iwrite.ge.1)then
          write(*,*)'CHECKW - Note:'
          write(*,*)'  Starting continuation component IWORK(2) at ',
     &      iwork(2)
        end if
      end if

      if(iwork(3).ne.1)iwork(3)=0
      if(iwork(4).lt.0.or.iwork(4).gt.2)iwork(4)=0
      if(iwork(5).lt.1.or.iwork(5).gt.nvar)iwork(5)=0
      if(iwork(6).lt.1.or.iwork(6).gt.nvar)iwork(6)=0
      if(iwork(9).lt.0.or.iwork(9).gt.2)iwork(9)=0
      if(iwork(17).lt.1)iwork(17)=10
c
c  Check entries of RWORK.
c
      if(rwork(1).le.zero)rwork(1)=sqrt(rwork(8))
      if(rwork(2).le.zero)rwork(2)=sqrt(rwork(8))

      if(rwork(3).le.0.0D+00)then
        rwork(3)=sqrt(rwork(8))
        if(iwork(7).ge.1)then
          write(*,*)'CHECKW - The minimum stepsize is too small.'
          write(*,*)'  Increasing RWORK(3) to ',rwork(3)
        end if
      end if

      if(rwork(4).lt.rwork(3))then
        temp=nvar
        rwork(4)=max(rwork(3),sqrt(temp))
      end if

      if(rwork(5).lt.zero)then
        rwork(5)=-rwork(5)
        rwork(6)=-rwork(6)
      end if

      if(rwork(5).lt.rwork(3).or.rwork(5).gt.rwork(4))then
        rwork(5)=(rwork(3)+rwork(4))/two
      end if

      if(rwork(6).ne.(-one))rwork(6)=one
      if(rwork(20).lt.one.or.rwork(20).gt.hundrd)rwork(20)=three

      return
      end
      subroutine coqual ( iwork, liw, rwork, lrw )

c*********************************************************************72
c
cc COQUAL computes the "correction quality".
c
c  Discussion:
c
c    COQUAL computes the factor QUAL which is based on the 'quality' of
c    the Newton correction iteration.  Considerations used include the
c    number of steps taken versus the maximum allowed and the ratio of 
c    the last corrector step to the total correction.
c
c    See the paper "On Steplength Algorithms for a Class of Continuation
c    Methods", listed in the documentation.
c
c    The quality factor, locally called QUAL, is stored in RWORK(23) 
c    on return.  Its value is between 1/8 and 8, with 1/8 signifying 
c    a poor correction process, 1 an average value, and 8 superior.  
c    The value of QUAL is a factor in the size of the next step.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Cor den Heijer, Werner Rheinboldt,
c    On Steplength Algorithms for a Class of Continuation Methods,
c    SIAM Journal on Numerical Analysis,
c    Volume 18, Number 5, October 1981, pages 925-947.
c
c  Parameters:
c
c    Input, INTEGER IWORK(LIW), the integer work array.
c
c    Input, INTEGER LIW, the size of IWORK.
c
c    Input/output, REAL RWORK(LRW), the real work array.
c
c    Input, INTEGER LRW, the size of RWORK.
c
      implicit none

      double precision eight
      double precision one
      double precision thgie

      parameter (eight=8.0D+00)
      parameter (one=1.0D+00)
      parameter (thgie=0.125D+00)

      integer liw
      integer lrw

      double precision base
      double precision bot
      double precision cordis
      double precision esab
      double precision expo
      integer iwork(liw)
      integer iwrite
      integer maxcor
      integer modnew
      integer nave
      integer nmax
      double precision qual
      double precision rwork(lrw)
      double precision stepx
      double precision term
      double precision test
      double precision top

      iwrite=iwork(7)
      stepx=rwork(9)
c
c  CORDIS, RWORK(15), is the Euclidean distance between the first and last
c  points in the Newton correction process.
c
c  STEPX,  RWORK(9), is the size of the last correction.
c
c  MAXCOR, IWORK(17), is the maximum number of Newton corrections allowed.
c
c  MODNEW, IWORK(4), defines the type of Newton corrector process used.
c
c  IWORK(28) contains the actual number of Newton correction steps used.
c
      cordis=rwork(15)
      modnew=iwork(4)
      maxcor=iwork(17)
c
c  Was the correction "minimal"?  Few steps, or a small total size?
c
      if((iwork(28).le.1).or.(cordis.le.eight*rwork(8))) then
        qual=eight
        go to 10
      end if
c
c  Were an "average" number of steps take?
c
      if(modnew.eq.0)then
        nave=(maxcor-1)/2
      else
        nave=maxcor
      end if

      if(iwork(28).eq.nave)then
        qual=one
        go to 10
      end if
c
c  Were many steps taken?
c
      if(modnew.eq.0)then
        nmax=maxcor
      else
        nmax=2*maxcor
      end if

      if(iwork(28).ge.nmax)then
        qual=thgie
        go to 10
      end if
c
c  For MODNEW=0, 
c    W=(STEPX/CORDIS),
c    IEXP=1/(2**(NCOR-1)-1)
c    U=W**IEXP
c    JEXP=2**(NCOR-NAVE)
c    QUAL=(U+1+(1/U)) / (U**JEXP+1+(1/U**JEXP))
c
      if ( modnew .eq. 0 ) then

        top = 2.0D+00**( iwork(28) - nave )
        bot = 2.0D+00**( iwork(28) - 1 ) - 1
        expo = 1.0D+00 / bot
        base = ( stepx / cordis )**expo

        if ( base .le. rwork(8) )then
          qual = 8.0D+00
        else
          top = base + 1.0D+00 + ( 1.0D+00 / base )
          term = base**top
          if ( term .le. rwork(8) ) then
            qual = 8.0D+00
          else
            bot = term + 1.0D+00 + ( 1.0D+00 / term )
            qual = top / bot
          end if
        end if
c
c  For MODNEW nonzero,
c    EXP=(NCOR-NAVE)/(NCOR-1)
c    W=(STEPX/CORDIS)
c    QUAL=W**EXP
c
      else

        top=iwork(28)-nave
        bot=iwork(28)-1
        expo=bot/top
        test=eight**expo
        base=stepx/cordis
        esab=cordis/stepx

        if((iwork(28).lt.nave.and.test.gt.base)
     &    .or.(iwork(28).gt.nave.and.test.lt.base))then
          qual=eight
          go to 10
        end if

        if((iwork(28).lt.nave.and.test.gt.esab)
     &    .or.(iwork(28).gt.nave.and.test.lt.esab))then
          qual=thgie
          go to 10
        end if

        expo=top/bot
        qual=base**expo

      end if
c
c  Store QUAL in RWORK
c
10    continue

      if(qual.gt.eight)qual=eight
      if(qual.lt.thgie)qual=thgie

      rwork(23)=qual

      return
      end
      subroutine corect ( df, fpar, fx, ierror, ihold, ipar, iwork,
     &  nvar, rwork, wk, xr, lrw, liw, icrit, slname )

c*********************************************************************72
c
cc CORECT applies a correction method to an approximate solution of the system.
c
c  Discussion:
c
c    The routine performs the Newton correction of an approximate solution X of
c    the vector equation F(X)=0.  It is required that the output value of
c    X satisfy this equation to within a certain tolerance.
c
c    Either Newton's method, or the chord Newton method may be used, depending
c    on a user chosen parameter.  In the latter case, the jacobian is only
c    evaluated at the starting point.
c
c    The system of NVAR-1 equations is temporarily augmented by an NVAR-th
c    equation which makes the system square and (presumably) nonsingular.
c    The equation is particularly simple:
c
c      X(IHOLD)=B
c
c    were B is some fixed value. In fact, the value B is set to the input value
c    of X(IHOLD).  This corresponds to simply holding the IHOLD-th entry of
c    X fixed during the iteration, which in turn, corresponds to treating
c    the IHOLD-th entry of X as a "parameter" in terms of which the other
c    entries of X may be solved for.
c
c    The linear system to be solved has the form
c
c     DFA(X,IHOLD) * (-DELX) = FA(X)
c
c    where the Jacobian DFA is augmented with an NVAR-th row containing a
c    1 in the IHOLD-th column, and the vector FA is augmented with an NVAR-th
c    value of 0.
c
c
c    After each Newton step, a decision is made whether to continue the iteration,
c    or to accept the current point, or to reject the entire Newton correction
c    process.  The criteria use the following parameters:
c
c      ABSERR is the user's absolute error tolerance.
c      RWORK(8) is the machine epsilon, also called EPMACH.
c      FMP is an adjustment factor, 2.0 if NCOR=1, 1.05 otherwise.
c      FNRM is the maximum-norm of the function value of X.
c      FNRML is the maximum-norm of the function value of the previous iterate.
c      MAXCOR is the maximum number of Newton iterations allowed.
c      NCOR is the iteration counter, with NCOR=0 for the starting point.
c      RELERR is the user's relative error tolerance.
c      STEPX is the maximum-norm of the difference of X and the previous iterate.
c      STEPXL is the previous value of STEPX.
c      X is the current Newton iterate.
c      XNRM is the maximum-norm of X.
c
c    At each step, the current point X will be accepted as the solution to
c    the nonlinear system if any of the following conditions hold:
c
c    Strong acceptance criterion:
c
c    1.  (FNRM.LE.ABSERR) and (STEPX.LE.(ABSERR+RELERR*XNRM))
c
c    Weak acceptance criteria:
c
c    2.  (NCOR.EQ.0) and (FNRM.LE.(.5*ABSERR))
c    3.  (FNRM.LE.8*EPMACH) or (STEPX.LE.8*EPMACH)
c    4.  (NCOR.GE.2) and
c        (FNRM+FNRML).LE.ABSERR and STEPX.LE.8*(ABSERR+RELERR*XNRM)
c    5.  (NCOR.GE.2) and
c        (FNRM.LE.8.0*ABSERR) and (STEPX+STEPXL).LE.(ABSERR+RELERR*XNRM)
c
c    The Newton iteration process is to be aborted if any of the following
c    criteria hold:
c
c    1.  FNRM.GT.(FMP*FNRML+ABSERR)
c    2.  (NCOR.GE.2) AND (ICRIT.EQ.0) and (STEPX.GT.(FMP*STEPXL+ABSERR))
c    3.  NCOR.GE.MAXCOR
c
c    Error conditions returned in IERROR:
c
c    0, no errors detected.
c    1, data or storage error.
c    2, error condition returned from user function or jacobian routine.
c    3, error condition returned from the solver routine.
c    4, a correction step was rejected (function norm increased or
c       size of correction increased)
c    5, too many correction steps were taken without convergence.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 May 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision eight
      double precision one
      double precision onefiv
      double precision two
      double precision zero

      parameter (eight=8.0D+00)
      parameter (one=1.0D+00)
      parameter (onefiv=1.05D+00)
      parameter (two=2.0D+00)
      parameter (zero=0.0D+00)

      external df
      external fx
      external idamax
      external daxpy
      external slname
      external dnrm2

      intrinsic abs

      integer liw
      integer lrw
      integer nvar

      double precision abserr
      double precision dets
      double precision fmp
      double precision fnrm
      double precision fnrml
      double precision fpar(*)
      integer i
      integer icrit
      integer ierror
      integer ifmax
      integer ihold
      integer ipar(*)
      integer idamax
      integer iwork(liw)
      integer iwrite
      integer ixmax
      integer job
      integer ksmax
      integer maxcor
      integer maxnew
      integer modnew
      double precision relerr
      double precision rwork(lrw)
      double precision dnrm2
      double precision stepxl
      double precision tlstep
      double precision wk(nvar)
      double precision xnrm
      double precision xr(nvar)
      double precision xvalue
c
c  Initialize.
c
      abserr=rwork(1)
      relerr=rwork(2)
      modnew=iwork(4)
      iwrite=iwork(7)
      maxcor=iwork(17)
      ierror=0
c
      iwork(28)=0

      if(modnew.eq.0)then
        maxnew=maxcor
      else
        maxnew=2*maxcor
      end if

      fmp=two
      rwork(9)=zero
      xvalue=xr(ihold)
c
c  Get the function value at the starting point.
c
      call fx(nvar,fpar,ipar,xr,wk,ierror)
      iwork(22)=iwork(22)+1

      if(ierror.ne.0)then
        if(iwrite.ge.1)then
          write(*,*)'CORECT - Fatal error!'
          write(*,*)'  User function returned IERROR=',ierror
        end if
        return
      end if

      wk(nvar)=xr(ihold)-xvalue

      ifmax=idamax(nvar,wk,1)
      fnrm=dnrm2(nvar,wk,1)

      ixmax=idamax(nvar,xr,1)
      xnrm=dnrm2(nvar,xr,1)

      if(iwrite.ge.2)then
        write(*,*)' '
        write(*,*)'Step      FX             X             DX'
        write(*,*)' '
        write(*,'(1x,i4,2g14.6)')iwork(28),fnrm,xnrm
        write(*,'(12x,i3,11x,i3)')ifmax,ixmax
      end if

      if(two*fnrm.le.abserr)return
c
c  Carry out up to MAXNEW Newton corrections.
c
      do i=1,maxnew

        iwork(28)=i
c
c  Solve FPRIME*WK=FX, where FX is the current residual or function value.
c
        job=0
        if(i.ne.1.and.i.ne.maxcor.and.modnew.eq.1)job=1
        if(modnew.eq.2)job=1

        call slname(dets,fx,df,fpar,ierror,ihold,ipar,iwork,liw,
     &    job,nvar,rwork,lrw,xr,wk)

        if(ierror.ne.0)then
          if(iwrite.ge.1)then
            write(*,*)'CORECT - Fatal error!'
            write(*,*)'  Solver returned IERROR=',ierror
          end if
          return
        end if
c
c  Subtract WK from XR to get the next iterate.
c
        call daxpy(nvar,-one,wk,1,xr,1)
        stepxl=rwork(9)
        ksmax=idamax(nvar,wk,1)
        rwork(9)=abs(wk(ksmax))

        ixmax=idamax(nvar,xr,1)
        xnrm=dnrm2(nvar,xr,1)
c
c  Compute function value at new iterate and take its norm.
c
        call fx(nvar,fpar,ipar,xr,wk,ierror)
        iwork(22)=iwork(22)+1

        if(ierror.ne.0) then
          if(iwrite.ge.1)then
            write(*,*)'CORECT - Fatal error!'
            write(*,*)'  User function returned IERROR=',ierror
          end if
          return
        end if

        wk(nvar)=xr(ihold)-xvalue
        fnrml=fnrm

        ifmax=idamax(nvar,wk,1)
        fnrm=dnrm2(nvar,wk,1)

        if(iwrite.ge.2) then
          write(*,'(1x,4x,28x,g14.6)')rwork(9)
          write(*,'(1x,4x,28x,7x,i3)')ksmax
          write(*,'(1x,i4,2g14.6)')iwork(28),fnrm,xnrm
          write(*,'(12x,i3,11x,i3)')ifmax,ixmax
        end if
c
c  Check for strong acceptance of function and stepsize.
c
        tlstep=abserr+relerr*xnrm
        if(fnrm.le.abserr.and.rwork(9).le.tlstep)return
c
c  TEMPORARY
c
        if(fnrm.le.abserr)return
c
c  Check for weak acceptance of function and stepsize.
c
        if(fnrm.le.eight*rwork(8).or.rwork(9).le.eight*rwork(8))return

        if(iwork(28).gt.1)then
          if((fnrm+fnrml).le.abserr.and.rwork(9).le.eight*tlstep)return
          if(fnrm.le.eight*abserr.and.(rwork(9)+stepxl).le.tlstep)return
        end if
c
c  Decide if iteration should be aborted
c
        if(iwork(28).gt.1)then

          if(icrit.lt.1.and.rwork(9).gt.(fmp*stepxl+abserr))then

            ierror=4

            if(iwrite.ge.2)then
              write(*,*)'CORECT - Warning!'
              write(*,*)'  The correction DX is not decreasing.'
            end if

            return

          end if

        end if

        if(icrit.lt.2.and.fnrm.gt.(fmp*fnrml+abserr)) then

           ierror=4

           if(iwrite.ge.2)then
             write(*,*)'CORECT - Warning!'
             write(*,*)'  The residual FX is not decreasing.'
           end if

           return
         end if

        fmp=onefiv

      enddo
c
c  Reached maximum number of steps without acceptance or rejection.
c
      ierror=5

      if(iwrite.ge.2)then
        write(*,*)'CORECT - Warning!'
        write(*,*)'  Convergence is too slow.'
      end if

      return
      end
      subroutine denjac ( eps, fpar, fprime, fx, ierror, ipar, ipc, 
     &  iwork, jac, liw, nvar, x, work1, work2 )

c*********************************************************************72
c
cc DENJAC estimates a dense jacobian.
c
c  Discussion:
c
c    The routine computes a matrix FPRIME which estimates the jacobin
c    of the function FX, assuming that full or "dense" storage is to be used. 
c    DENJAC is called  by DENSLV when the user jacobian option JAC is set to 1 or 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c  EPS    Input, REAL EPS, a tolerance to be used for shifting the X 
c         values during the finite differencing.  No single value of EPS
c         will be reliable for all vectors X and functions FX.  Values of
c         EPS have typically been chosen between SQRT(EPSMCH) and
c         SQRT(SQRT(EPSMCH)) where EPSMCH is the machine tolerance.
c
c  FPAR   Input, REAL FPAR(*), a real vector, available for the user
c         to communicate with the FX routine.
c
c  FPRIME Output, REAL FPRIME(NVAR,NVAR), the array into which the
c         jacobian will be stored, including the augmenting last row.
c
c  FX     Input, EXTERNAL FX, the name of the user supplied routine
c         that defines the nonlinear function, and which has the form:
c
c           subroutine fx(nvar,fpar,ipar,x,f,ierror)
c
c           integer nvar
c
c           real f(nvar)
c           real fpar(*)
c           integer ierror
c           integer ipar(*)
c           real x(nvar)
c
c           do i=1,nvar-1
c             f(i)=function(i)(x)
c           enddo
c
c           return
c           end
c
c  IERROR Output, INTEGER IERROR, error return flag.  A nonzero value
c         indicates that an error occurred.
c
c  IPAR   Input, INTEGER IPAR(*), an integer vector, available for the user
c         to communicate with the FX routine.
c
c  IPC    Input, INTEGER IPC, the index of the current continuation parameter,
c         which determines the location of the "1" entry in the last row
c         of the jacobian.
c
c  IWORK  Input/output, INTEGER IWORK(LIW), an array containing pointers
c         and counters used by the continuation code.  In particular,
c         IWORK(22) contains a count of the number of function evaluations.
c
c  JAC    Input, INTEGER JAC, the jacobian option.
c
c         0, the user has supplied the jacobian routine DF.
c
c         1, the program must estimate the jacobian using forward differences.
c
c         2, the program must estimate the jacobian using central differences.
c
c  LIW    Input, INTEGER LIW, the dimension of IWORK.
c
c  NVAR   Input, INTEGER NVAR, the number of variables.
c
c  X      Input, REAL X(NVAR), the point at which the jacobian is to be
c         estimated.
c
c  WORK1,
c  WORK2  Workspace, REAL WORK1(NVAR), WORK2(NVAR).
c
      implicit none

      double precision one
      parameter (one=1.0D+00)

      double precision zero
      parameter (zero=0.0D+00)

      intrinsic abs
      external fx
      external daxpy
      external dscal

      integer liw
      integer nvar

      double precision delm
      double precision delp
      double precision eps
      double precision fpar(*)
      double precision fprime(nvar,nvar)
      integer ipar(*)
      integer ierror
      integer ipc
      integer iwork(liw)
      integer iwrite
      integer j
      integer jac
      double precision skale
      double precision x(nvar)
      double precision xsave
      double precision work1(nvar)
      double precision work2(nvar)

      iwrite=iwork(7)
c
c  If we are using forward differences, then evaluate the function F
c  at the base point X, and save the value for use in the difference
c  quotient.
c
      if(jac.eq.1)then
        call fx(nvar,fpar,ipar,x,work2,ierror)
        delm=zero
        iwork(22)=iwork(22)+1
        if(ierror.ne.0)then
          if(iwrite.gt.0)then
            write(*,*)
     &        'DENJAC - User function returns IERROR=',ierror
          end if

          return
        end if
      end if
c
c  Increment each variable X(J) by a small amount, and evaluate the
c  function there.
c
      do j=1,nvar

        xsave=x(j)
        delp=eps*(one+abs(x(j)))
        x(j)=x(j)+delp
        call fx(nvar,fpar,ipar,x,work1,ierror)
        iwork(22)=iwork(22)+1
        if(ierror.ne.0)then
          if(iwrite.gt.0)then
            write(*,*)
     &        'DENJAC - User function returns IERROR=',ierror
          end if
          return
        end if
c
c  For central difference approximations, decrement each variable X(J)
c  by a small amount, and evaluate the function there.
c
        if(jac.eq.2)then
          delm=-delp
          x(j)=xsave+delm
          call fx(nvar,fpar,ipar,x,work2,ierror)
          iwork(22)=iwork(22)+1

          if(ierror.ne.0)then

            if(iwrite.gt.0)then
              write(*,*)
     &          'DENJAC - User function returns IERROR=',ierror
            end if

            return
          end if

        end if

        x(j)=xsave
c
c  Compute DFDX(*,J)= (F(X+)-F(X-))/DELX.
c
c  Note that we contrive to ADD this quantity to DFDX, without
c  overwriting anything that may already be in DFDX.
c
c  This makes it possible to use this routine also to check
c  the accuracy of the user's jacobian routine.
c
        call daxpy(nvar-1,-one,work2,1,work1,1)
        skale=one/(delp-delm)
        call dscal(nvar-1,skale,work1,1)
        call daxpy(nvar-1,one,work1,1,fprime(1,j),1)

      enddo

      fprime(nvar,ipc)=fprime(nvar,ipc)+one

      return
      end
      subroutine denslv ( dets, fx, df, fpar, ierror, ipc, ipar, iwork, 
     &  liw, job, nvar, rwork, lrw, x, y ) 

c*********************************************************************72
c
cc DENSLV solves a dense linear system
c
c  Discussion:
c
c    This routine manages the solution of the linear system
c
c           (  DF(X)  )
c    (1)    (---------) * Y   = Y
c           (  E(IPC) )
c
c    where the NVAR-1 by NVAR submatrix DF(X) is the jacobian of the
c    nonlinear function F(X), or an approximation thereto, and E(IPC)
c    is a row vector of NVAR entries, consisting of all zeroes except
c    for a 1 in column IPC.
c
c    As a special application, DENSLV may also be called to compare
c    the user supplied jacobian to a finite difference approximation.
c
c    DENSLV is used when the jacobian DF(X) is assumed to be dense.
c    In this case, the corresponding matrix is stored in NVAR*NVAR
c    consecutive entries in the RWORK array, RWORK(LRBEG) through
c    RWORK(LRBEG+NVAR*NVAR-1), with the value of LRBEG stored in
c    IWORK(15).
c
c    If DENSLV is required to factor the matrix, the factored version will
c    overwrite the original matrix, using the same storage area.  As part of
c    the factorization, DENSLV will also need a vector of length NVAR
c    to store the pivot information.  This is done in the tail end of the
c    IWORK array.  The first entry of IWORK devoted to this is recorded
c    in IWORK(13).  In other words, the pivot information starts in
c    IWORK(IWORK(13)).
c
c    If the jacobian of the nonlinear function is dense, the user should
c    specify DENSLV as the external solver argument SLNAME.  The user
c    should either request a finite difference jacobian be calculated,
c    or provide a jacobian routine of the following form:
c
c        subroutine df(nvar,fpar,ipar,x,a,ierr)
c
c        integer nvar
c
c        real a(nvar,nvar)
c        real fpar(*)
c        integer ierr
c        integer ipar(*)
c        real x(nvar)
c
c        do i=1,nvar-1
c          do j=1,nvar
c            a(i,j)=df(i)/dx(j)
c          enddo
c        enddo
c 
c        return
c        end
c
c    where "DF(I)/DX(J)" denotes the derivative of the I-th component of the
c    nonlinear function F(X) with respect to the J-th component of X.
c    The user need not supply the NVAR-th, augmenting, row of the matrix, since
c    this is taken care of automatically.
c
c    DENSLV will directly call the LAPACK routines SGETRF and SGETRS, and
c    the routine SGEDET.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, REAL DETS, the sign of the determinant of the matrix
c    in equation (1).
c
c    Input, EXTERNAL FX, the name of the user supplied routine
c    which evaluates the nonlinear function.
c
c    Input, EXTERNAL DF, the name of the user supplied routine
c    which evaluates the jacobian matrix.
c
c    Input, REAL FPAR(*), a real array which allows the user to pass
c    double precision parameters "through" PITCON to FX and DF.
c
c    Output, INTEGER IERROR, error flag.
c    * 0, normal return.
c    * 1, data or storage error.
c    * 2 , error returned by the derivative routine DF.
c    * 3, the matrix containing the augmented jacobian is singular.
c
c    Input, INTEGER IPC, the index of the continuation parameter,
c    which determines the location of the "1" entry in the final
c    row of the jacobian.
c
c    Input, INTEGER IPAR(*), an integer array available for the user
c    to transmit information to DF and FX.
c
c    Workspace, INTEGER IWORK(LIW), an array containing pointers,
c    counters, and space for a pivot array.
c
c    Input, INTEGER LIW, the dimension of IWORK.
c
c    Input, INTEGER JOB, the action switch.
c     * 0, Evaluate jacobian, factor jacobian, compute determinant,
c     and solve linear system.
c     * 1, solve linear system with previously factored jacobian.
c     * 2, Evaluate jacobian, factor jacobian, compute determinant.
c     * 3, Check jacobian matrix.  Call user jacobian routine,
c     multiply by -1.0, add finite difference jacobian,
c     print largest entry.
c
c    Input, INTEGER NVAR, the number of variables.
c
c    Input/output, REAL RWORK(*).  RWORK contains various scalars,
c    vectors, and the jacobian matrix.  On input, RWORK may contain
c     a previously factored jacobian matrix.  On output, RWORK
c     will generally contain the factored jacobian matrix.
c
c    Input, INTEGER LRW, the dimension of RWORK.
c
c    Input, REAL X(NVAR), the point at which the jacobian is to be
c     evaluated.
c
c    Input/output, REAL Y(NVAR).  On input, Y may contain the right
c     hand side of a linear system to be solved, in which case Y
c     will contain the solution of that linear system on output.
c
      implicit none

      double precision one
      double precision zero

      parameter (one=1.0D+00)
      parameter (zero=0.0D+00)

      external denjac
      external df
      external fx
      external idamax
      intrinsic mod
      external dgedet
      external dgetrf
      external dgetrs
      intrinsic sqrt
      external dscal

      integer liw
      integer lrw
      integer nvar

      double precision det
      double precision dets
      double precision fpar(*)
      integer i
      integer ierror
      integer info
      integer ipar(*)
      integer ipc
      integer idamax
      integer iwork(liw)
      integer iwrite
      integer j
      integer jac
      integer jack
      integer job
      integer k
      integer ldf
      integer lfxm
      integer lfxp
      integer lilst
      integer lpiv
      integer lrlst
      integer ndim
      double precision rwork(lrw)
      double precision x(nvar)
      double precision y(nvar)

      ierror=0
      iwrite=iwork(7)
      lpiv=iwork(13)
      ldf=iwork(15)
      jac=iwork(9)
c
c  Make sure that the jacobian routine is available if the user
c  has specified an option that requires it!
c
      if(jac.ne.0.and.
     &  (iwork(1).eq.(-1).or.iwork(1).eq.(-2).or.
     &   iwork(1).eq.(-3).or.iwork(1).eq.(-4).or.
     &   iwork(1).eq.(-9).or.iwork(1).eq.(-10)))then
        ierror=4

        if(iwrite.ge.1)then
          write(*,*)'DENSLV - Fatal error!'
          write(*,*)'An option was selected requiring'
          write(*,*)'a user jacobian, but the value of'
          write(*,*)'IWORK(9) indicates no such routine'
          write(*,*)'is available!'
        end if

        return
      end if
c
c  Check the amount of storage available
c
      lilst=lpiv+nvar-1

      if(jac.eq.0.and.job.ne.3)then
        lrlst=ldf-1+nvar*nvar
      else
        lrlst=ldf+nvar*nvar+2*nvar-1
      end if
 
      if(lilst.gt.liw)then

        ierror=1

        if(iwrite.ge.1)then
          write(*,*)'DENSLV - Fatal error!'
          write(*,*)
     &      'Insufficient integer workspace in IWORK!'
          write(*,*)'Need workspace size    =',lilst
          write(*,*)'Available workspace LIW=',liw
        end if

        return
      end if
 
      if(lrlst.gt.lrw)then
        ierror=1

        if(iwrite.ge.1)then
          write(*,*)'DENSLV - Fatal error!'
          write(*,*)
     &      'DENSLV - Insufficient real workspace in RWORK!'
          write(*,*)'Need workspace size    =',lrlst
          write(*,*)'Available workspace LRW=',lrw
        end if

        return
      end if
c
c  JOB=1.
c  A linear system is to be solved.  The matrix has already been factored,
c  and the right hand side is supplied.  Return the solution.
c
      if(job.eq.1)then
        call dgetrs('n',nvar,1,rwork(ldf),nvar,iwork(lpiv),y,nvar,info)
        iwork(21)=iwork(21)+1
        if(info.ne.0)then
          if(iwrite.ge.1)then
            write(*,*)'DENSLV - DGETRS returns INFO=',info
          end if
          ierror=3
        end if
        return
      end if

      ndim=nvar*nvar

      do i=1,ndim
        rwork(ldf+i-1)=zero
      enddo
c
c  If a user jacobian routine is available, invoke it.
c
      if(jac.eq.0)then
        if(iwork(1).gt.(-5).or.iwork(1).lt.(-8))then
          call df(nvar,fpar,ipar,x,rwork(ldf),ierror)
          iwork(19)=iwork(19)+1
          rwork(ldf+ipc*nvar-1)=one
        end if
      end if
c
c  If we're going to compare the user jacobian to an approximate one,
c  negate the user jacobian now.
c
      if(job.eq.3)then
        if(iwork(1).le.(-1).and.iwork(1).ge.(-4))then
          call dscal(ndim,-one,rwork(ldf),1)
        end if
      end if
c
c  If we don't have a user jacobian, or we need to compare the
c  user jacobian to an approximate one, get the approximation.
c
      if((iwork(1).ge.0.and.(jac.eq.1.or.jac.eq.2)).or.
     &   (iwork(1).le.(-1).and.iwork(1).ge.(-8)))then

        jack=jac

        if(job.eq.3)then
          if(mod(iabs(iwork(1)),2).eq.1)then
            jack=1
          else
            jack=2
          end if
        end if

        lfxp=ldf+nvar*nvar
        lfxm=ldf+nvar*nvar+nvar

        call denjac(rwork(18),fpar,rwork(ldf),fx,ierror,ipar,ipc,iwork,
     &    jack,liw,nvar,x,rwork(lfxp),rwork(lfxm))

        if(ierror.ne.0)then

          if(iwrite.gt.0)then
            write(*,*)'DENSLV - Fatal error!'
            write(*,*)'  DENJAC returned IERROR=',ierror
          end if

          return
        end if
      end if
c
c  If JOB=3, print out information:
c    IWORK(1)=-1, -2, -5, -6, -9, print largest entry of matrix.
c    IWORK(1)=-3, -4, -7, -8, -10, also print out entire matrix.
c
      if(job.eq.3)then
        k=idamax(ndim,rwork(ldf),1)
        i=mod(k-1,nvar)+1
        j=(k-i)/nvar+1

        if(iwork(1).le.(-1).and.iwork(1).ge.(-4))then
          write(*,*)' '
          write(*,*)
     &      'DENSLV - Maximum value of FP_Approx(I,J)-FP_User(I,J)'
        elseif(iwork(1).le.(-5).and.iwork(1).ge.(-8))then
          write(*,*)' '
          write(*,*)
     &      'DENSLV - Maximum value of finite difference jacobian:'
        elseif(iwork(1).le.(-9).and.iwork(1).ge.(-10))then
          write(*,*)' '
          write(*,*)
     &      'DENSLV - Maximum value of user supplied jacobian:'
        end if

        write(*,*)rwork(ldf+k-1),' I, J=',i,j
        write(*,*)' '

        if(iwork(1).eq.(-3).or.iwork(1).eq.(-4).or.
     &     iwork(1).eq.(-7).or.iwork(1).eq.(-8).or.
     &     iwork(1).eq.(-10))then
          if(iwork(1).eq.(-3).or.iwork(1).eq.(-4))then
            write(*,*)' '
            write(*,*)'DENSLV - Entire difference matrix:'
            write(*,*)'         FP_Approx(I,J)-FP_User(I,J)'
            write(*,*)' '
          elseif(iwork(1).eq.(-7).or.iwork(1).eq.(-8))then
            write(*,*)' '
            write(*,*)'DENSLV - Finite difference jacobian:'
            write(*,*)' '
          elseif(iwork(1).eq.(-10))then
            write(*,*)' '
            write(*,*)'DENSLV - User supplied jacobian:'
            write(*,*)' '
          end if

          do i=1,nvar
            do j=1,nvar
              k=ldf+(j-1)*nvar+i-1
              write(*,*)rwork(k),' i, j=',i,j
            enddo
            write(*,*)' '
          enddo

        end if
        return
      end if
c
c  Decompose matrix into LU factors.
c
      call dgetrf(nvar,nvar,rwork(ldf),nvar,iwork(lpiv),info)

      iwork(20)=iwork(20)+1
      if(info.ne.0)then
        if(iwrite.ge.1)then
          write(*,*)'DENSLV - Fatal error!'
          write(*,*)'  Zero pivot, DGETRF returns INFO=',info
        end if
        ierror=3
        return
      end if
c
c  Compute matrix determinant, and record its sign.
c
      call dgedet(nvar,rwork(ldf),nvar,iwork(lpiv),det)
      dets=zero

      if(det.gt.zero)then
        dets=one
      elseif(det.lt.zero)then
        dets=-one
      end if

      if(job.eq.2)return
c
c  Solve linear system, overwriting right hand side with solution.
c
      call dgetrs('n',nvar,1,rwork(ldf),nvar,iwork(lpiv),y,nvar,info)
      iwork(21)=iwork(21)+1
      if(info.ne.0)then
        if(iwrite.ge.1)then
          write(*,*)'DENSLV - Fatal error!'
          write(*,*)'  DGETRS returns INFO=',info
        end if
        ierror=3
        return
      end if

      return
      end
      subroutine dgbdet ( m, ml, mu, a, ldab, ipivot, det )

c*********************************************************************72
c
cc DGBDET computes the determinant of a band matrix factored by DGBTRF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, INTEGER M, the number of rows and columns in the matrix.
c
c    Input, INTEGER ML, the number of nonzero subdiagonals in the matrix.
c
c    Input, INTEGER MU, the number of nonzero superdiagonals in the matrix.
c
c    Input, double precision A(LDAB,M), the band matrix, as factored by SGBTRF.
c
c    Input, INTEGER IPIVOT(M), the pivot vector, as computed by SGBTRF.
c
c    Output, double precision DET, the determinant of the matrix.
c
      implicit none

      integer ldab
      integer m

      double precision a(ldab,m)
      double precision det
      integer i
      integer ipivot(m)
      integer ml
      integer mu

      det=1.0D+00
 
      do i=1,m
        if(ipivot(i).ne.i)det=-det
      enddo
 
      do i=1,m
        det=a(ml+mu+1,i)*det
      enddo
 
      return
      end
      subroutine dgedet ( m, a, lda, ipivot, det )

c*********************************************************************72
c
cc DGEDET computes the determinant of a matrix factored by DGETRF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c  M      Input, INTEGER M.
c         The number of rows and columns of the matrix A.
c
c  A      Input, double precision A(LDA,M)
c         The factors L and U as computed by DGETRF.
c
c  LDA    Input, INTEGER LDA,
c         The leading dimension of the array A.  LDA >= max(1,M).
c
c  IPIVOT Input, INTEGER IPIVOT(M), as computed by DGETRF.
c
c  DET    Output, double precision DET, the determinant of the 
c         matrix A.
c
      implicit none

      integer lda
      integer m

      double precision a(lda,m)
      double precision det
      integer i
      integer ipivot(m)

      det=1.0D+00
 
      do i=1,m
        det=det*a(i,i)
      enddo
 
      do i=1,m
        if(ipivot(i).ne.i)det=-det
      enddo
 
      return
      end
      subroutine limit ( df, fpar, fx, ierror, ipar, iwork, liw, 
     &  lrw, nvar, rwork, slname, tc, wk, xc, xf, xr )

c*********************************************************************72
c
cc LIMIT is used to seek a limit point between two continuation points.
c
c  Discussion:
c
c    The continuation points are XF and XC, and there is also a tangent
c    vector at XF in TC, and a tangent vector at XC if WK.  It is assumed
c    that the LIM-th components of these tangent vectors differ in sign,
c    indicating the presence of a limit point having a LIM-th component
c    which is exactly zero.
c
c    We solve this problem using a one-dimensional zero finder.  We
c    set up a variable SN, which measures the proportional length along
c    the secant between XF and XC.  XF can now be thought of as the
c    point on the curve corresponding to SN=0, and XC as the point
c    corresponding to SN=1.  For other values of SN, we compute a linear
c    combination of XF and XC, and use Newton iteration.
c
c    During this Newton iteration, we must fix a component of the solution.
c    The component is chosen as the index of the entry of largest magnitude
c    in the secant.  However, should that index be LIM, then the second
c    largest magnitude will be chosen instead.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 May 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      intrinsic abs
      external corect
      external df
      external fx
      external idamax
      intrinsic max
      external root
      external daxpy
      external dcopy
      intrinsic sign
      external slname
      external dnrm2
      external tangnt

      double precision eight
      parameter (eight=8.0D+00)

      integer maxit
      parameter (maxit=25)

      double precision one
      parameter (one=1.0D+00)

      double precision two
      parameter (two=2.0D+00)

      double precision zero
      parameter (zero=0.0D+00)

      integer liw
      integer lrw
      integer nvar

      double precision a
      double precision b
      double precision dirlpc
      double precision fa
      double precision fb
      double precision fpar(*)
      integer icrit
      integer ierror
      integer iflag
      integer imitl
      integer ipar(*)
      integer idamax
      integer iwork(liw)
      integer iwrite
      integer lim
      integer lpc
      integer modsav
      double precision rwork(lrw)
      double precision skale
      double precision sn
      double precision snl
      double precision dnrm2
      double precision tc(nvar)
      double precision temp
      double precision tsn
      double precision wk(nvar)
      double precision xabs
      double precision xc(nvar)
      double precision xdif
      double precision xf(nvar)
      double precision xr(nvar)

      iwrite=iwork(7)
      lim=iwork(6)

      if(iwrite.ge.2)write(*,*)
     &  'LIMIT  - Attempt correction of approximate limit point.'
      modsav=iwork(4)
c
c  The limit point occurs somewhere between the current and previous
c  points.  See if either already computed point can be accepted as
c  the limit point.  If the current point XC is the limit point, then
c  its tangent is already in WK.
c
      if(abs(rwork(27)).le.rwork(1)/two) then
        call dcopy(nvar,xc,1,xr,1)
        go to 30
      end if

      if(abs(rwork(26)).le.rwork(1)/two) then
        call dcopy(nvar,xf,1,xr,1)
        call dcopy(nvar,tc,1,wk,1)
        go to 30
      end if
c
c  If interval is extremely small, simply assign one
c  endpoint of the interval as the answer.
c
      xdif=abs(xc(lim)-xf(lim))
      xabs=max(abs(xc(lim)),abs(xf(lim)))

      if(xdif.le.eight*rwork(8)*(one+xabs)) then
        if(abs(rwork(27)).gt.abs(rwork(26))) then
          call dcopy(nvar,xf,1,xr,1)
          call dcopy(nvar,tc,1,wk,1)
        else
          call dcopy(nvar,xc,1,xr,1)
        end if
        go to 30
      end if
c
c  Begin root-finding iteration on interval (0,1), with function values
c  TLLIM and TCLIM.
c
      a=zero
      fa=rwork(27)
      b=one
      fb=rwork(26)
c
c  Find LPC, the index of the entry of maximum absolute value in the
c  secant between the two continuation points.  However, we will not
c  allow LPC to equal LIM.  We save the sign of the LPC-th entry of
c  the secant so that new tangents may be properly signed.
c
c  Note: If the value of IWORK(3) is 1, then the user has requested that
c  the parameterization index always be held fixed at the value in IWORK(2).
c  The user's choice overrides all considerations, even in this routine.
c  Certain doom would occur if the user chooses unwisely, but that is
c  not our concern!
c
      if(iwork(3).ne.1)then
        temp=xc(lim)
        xc(lim)=xf(lim)
        call dcopy(nvar,xf,1,xr,1)
        call daxpy(nvar,-one,xc,1,xr,1)
        lpc=idamax(nvar,xr,1)
        xc(lim)=temp
      else
        lpc=iwork(2)
      end if

      dirlpc=sign(one,xf(lpc)-xc(lpc))
c
c  The first approximation to the limit point will be whichever endpoint
c  has the smallest LIM-th component of the tangent vector.
c
      if(abs(rwork(26)).ge.abs(rwork(27))) then
        sn=zero
        tsn=rwork(27)
        call dcopy(nvar,xc,1,xr,1)
      else
        sn=one
        tsn=rwork(26)
        call dcopy(nvar,xf,1,xr,1)
      end if

      imitl=0
      if(iwrite.ge.2)then
        write(*,*)
     &    'LIMIT  - For S=',zero,' Tan(X(S))(Lim)=',rwork(27)
        write(*,*)
     &    'LIMIT  - For S=',one,' Tan(X(S))(Lim)=',rwork(26)
      end if
c
c  Call the rootfinder repeatedly for the approximate root SN.
c  Use a linear combination of the points X(SNL) and X(0.0) or X(1.0)
c  to get a starting point for correction to the curve, returning
c  to the curve along the line X(LPC)=constant.  Compute the
c  tangent vector there, and return the LIM-th component of the
c  tangent as the function value whose zero we seek.
c
10    continue

      snl=sn
      call root(a,fa,b,fb,sn,tsn,imitl,iflag,ierror,rwork(8))
      iwork(23)=iwork(23)+1

      if(ierror.ne.0)then
        if(iwrite.ge.1)then
          write(*,*)'LIMIT - Fatal error!'
          write(*,*)'  Rootfinder returns IERROR=',ierror
        end if
        ierror=7
        return
      end if

      if(iflag.eq.(-1).or.iflag.eq.0)go to 30
c
c  Find whether SN lies in (0.0,SNL) or (SNL,1.0).  This will determine
c  how we construct our linear combination of the current points to
c  get a starting point for X(SN).  This somewhat cumbersome procedure
c  occurs because we limit the number of vectors we use.
c
      if(sn.le.snl) then
c
c  If SN lies between 0.0 and SNL, then set the approximation to the point
c  on the curve parameterized by the value SN as:
c
c  X(SN) = (SNL-SN)/(SNL-0.0) * X(0.0) + (SN-0.0)/(SNL-0.0) * X(SNL)
c
        if(snl.le.zero)then
          skale=zero
        else
          skale=sn/snl
          if(skale.le.zero)skale=zero
          if(skale.ge.one)skale=one
        end if

        call dscal(nvar,skale,xr,1)
        skale=one-skale
        call daxpy(nvar,skale,xc,1,xr,1)
c
c  Otherwise, if SN lies between SNL and 1.0, set
c
c  X(SN)=(SN-SNL)/(1.0-SNL)*X(1.0)+(1.0-SN)/(1.0-SNL)*X(SNL)
c
      else

        if(snl.ge.one)then
          skale=zero
        else
          skale=(one-sn)/(one-snl)
          if(skale.le.zero)skale=zero
          if(skale.ge.one)skale=one
        end if

        call dscal(nvar,skale,xr,1)
        skale=one-skale
        call daxpy(nvar,skale,xf,1,xr,1)

      end if
c
c  Try to correct the approximate point so that it lies on the curve.
c  If the user is trying to economize, by using a nonzero value of
c  IWORK(4), then we may try to recover from a failed correction by
c  retrying it with a smaller value of IWORK(4).
c
20    continue
      icrit=0
   
      call corect(df,fpar,fx,ierror,lpc,ipar,iwork,
     &  nvar,rwork,wk,xr,lrw,liw,icrit,slname)
c
c  If the correction fails, see if we can retry.
c
      if(ierror.ne.0.and.iwork(4).gt.0)then
        ierror=0
        iwork(4)=iwork(4)-1
        call dcopy(nvar,xc,1,xr,1)
        call dscal(nvar,one-sn,xr,1)
        call daxpy(nvar,sn,xf,1,xr,1)
        if(iwrite.ge.1)write(*,*)
     &    'LIMIT  - Retry limit computation with IWORK(4)=',iwork(4)
        go to 20
      end if

      iwork(4)=modsav
      if(ierror.ne.0)then
        if(iwrite.ge.1)then
          write(*,*)'LIMIT - Fatal error!'
          write(*,*)'  Corrector returned IERROR=',ierror
        end if
        ierror=7
        return
      end if
c
c  Compute the tangent at the new point.
c
      call tangnt(temp,fx,df,fpar,ierror,lpc,ipar,iwork,nvar,
     &  rwork,wk,xr,liw,lrw,slname)

      if(ierror.ne.0)then
        if(iwrite.ge.1)then
          write(*,*)'LIMIT - Fatal error!'
          write(*,*)'  TANGNT returned IERROR=',ierror
        end if
        ierror=7
        return
      end if
c
c  Adjust the sign of the tangent vector so the LPC-th component
c  has the same sign as the LPC-th component of the secant.
c
      if(dirlpc.ne.sign(one,wk(lpc))) then
        call dscal(nvar,-one,wk,1)
      end if
c
c  See if we can accept the new point as the actual limit point, because
c  the LIM-th component of the tangent is acceptably small.
c
      tsn=wk(lim)
      if(iwrite.ge.2)write(*,*)
     &  'LIMIT  - For S=',sn,' Tan(X(S))(Lim)=',tsn
      if(abs(tsn).le.rwork(1))go to 30
c
c  See if we have reached MAXIT iterations
c
      if(imitl.lt.maxit)go to 10
c
c  The limit point iteration has not produced a point which satisfies
c  our requirement.  We set an error flag, but we return the partial results
c  of the abortive computation anyway.
c
      ierror=8
      if(iwrite.ge.1)then
        write(*,*)'LIMIT - Warning!'
        write(*,*)'  Iteration did not reach limit point after '
        write(*,*)'  taking ',maxit,' steps.'
      end if
c
c  The limit point iteration is over.
c  Compute and store information.
c
30    continue
      call dcopy(nvar,xr,1,wk,1)
      call daxpy(nvar,-one,xc,1,wk,1)
      rwork(14)=rwork(12)+dnrm2(nvar,wk,1)
      iwork(27)=iwork(27)+1
      iwork(1)=4

      return
      end
      function r8_epsilon ( )

c*********************************************************************72
c
cc R8_EPSILON returns the R8 roundoff unit.
c
c  Discussion:
c
c    The roundoff unit is a number R which is a power of 2 with the
c    property that, to the precision of the computer's arithmetic,
c      1 .lt. 1 + R
c    but
c      1 = ( 1 + R / 2 )
c
c    FORTRAN90 provides the superior library routine
c
c      EPSILON ( X )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 March 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision R8_EPSILON, the R8 roundoff unit.
c
      implicit none

      double precision r8
      double precision r8_epsilon
      double precision r8_test

      r8 = 1.0D+00
      r8_test = 1.0D+00 + ( r8 / 2.0D+00 )

10    continue

      if ( 1.0D+00 .lt. r8_test ) then
        r8 = r8 / 2.0D+00
        r8_test = 1.0D+00 + ( r8 / 2.0D+00 )
        go to 10
      end if

      r8_epsilon = r8

      return
      end
      subroutine root ( a, fa, b, fb, u, fu, kount, iflag, ierror, 
     &  epmach )

c*********************************************************************72
c
cc ROOT seeks a root of the scalar equation F(X)=0.0.
c
c  Discussion:
c
c    ROOT must be called repeatedly to find the root.  On the first call,
c    ROOT is given a starting interval (A,B), on which F changes sign,
c    and the function values FA and FB.  ROOT returns a new point U
c    at which the value of F is to be computed.
c
c    The user may accept U as the root, or more likely, return FU
c    to ROOT, allowing it to make a new and better guess for the root.
c
c    See the book by Brent for details.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c  A      Input/output, REAL A.  A is one endpoint of the current 
c         interval in which the root is sought.  The user must set this 
c         before the first call only.  Thereafter, ROOT adjusts A as 
c         the interval shrinks.
c
c  FA     Input/output, REAL FA.  On the very first call only, the 
c         user must set FA to the value of F(A).  Thereafter, the 
c         program resets FA as the value of A changes, and the user 
c         should not alter the value of FA.
c
c  B      Input/output, REAL B.  B is one endpoint of the current 
c         interval in which the root is sought.  The user must set 
c         this before the first call only.  Thereafter, ROOT adjusts
c         B as the interval shrinks.
c
c  FB     Input/output, REAL FB.  On the very first call only, the 
c         user must set FB to the value of F(B).  Thereafter, the 
c         program resets FB as the value of B changes, and the user 
c         should not alter the value of FB.
c
c  U      Output, REAL U, the current approximation to the root.  After every
c         call to ROOT, U will contain the routine's best approximation to
c         the location of the root.
c
c  FU     Input, REAL FU.  On the first call, the user should not set FU.
c         Before the second call, the user should evaluate the function at
c         the point U, returned on the first call, and send that value as FU.
c         This behavior should then be repeated on each subsequent call.
c         The output value of U should be used to evaluate the function, with
c         the result sent back as FU on the next call.
c
c  KOUNT  Input/output, INTEGER KOUNT.  On the first call, set KOUNT to zero.
c         Thereafter, the program will update KOUNT, which counts the number
c         of calls made to ROOT.
c
c  IFLAG  Output, INTEGER IFLAG, reports the status of the search.
c
c         -1, the current bracketing interval (A,B) or (B,A) is smaller than
c         4*EPMACH*ABS(U)+EPMACH, and so U should be accepted as the root.
c
c         0, the input value FU is exactly zero, so U should be accepted
c         as the root.
c
c         Positive values of IFLAG indicated that the program has found a
c         new approximation U to the root.  If a better approximation is
c         desired, return the value of the function in FU, and the program
c         will proceed with the search.  The actual value of IFLAG tells
c         what method was used to produce the current U:
c
c         1, bisection
c         2, linear interpolation
c         3, inverse quadratic interpolation.
c
c  IERROR Output, INTEGER IERROR, error flag.
c         0, no error occurred.
c         7, on first call, FA*FB is greater than 0, and so the given
c         interval is unacceptable.
c
c  EPMACH Input, REAL EPMACH, the relative machine precision.
c
      implicit none

      double precision eight
      double precision one
      double precision onep5
      double precision two
      double precision zero

      parameter (eight=8.0D+00)
      parameter (one=1.0D+00)
      parameter (onep5=1.5D+00)
      parameter (two=2.0D+00)
      parameter (zero=0.0D+00)

      intrinsic abs
      intrinsic sign

      double precision a
      double precision b
      double precision epmach
      double precision fa
      double precision fb
      double precision fu
      double precision halfub
      integer ierror
      integer iflag
      integer kount
      double precision p
      double precision q
      double precision r
      double precision s
      double precision sdel1
      double precision sdel2
      double precision sdel3
      double precision sdel4
      double precision step
      double precision toler
      double precision u
c
c  Segment 1.  The first call is handled specially.
c
      if(kount.le.0)then

        if((fa.gt.zero.and.fb.gt.zero)
     &    .or.(fa.lt.zero.and.fb.lt.zero))then

          ierror=7
          kount=0
          return

        end if

        kount=1
        sdel1=two*abs(b-a)
        sdel2=two*sdel1
        sdel3=two*sdel2
        u=b
        b=a
        fu=fb
        fb=fa

      else
c
c  On calls after the first call, increment the counter, and check
c  whether F(U) is zero.
c
        kount=kount+1

        if(fu.eq.zero)then
          iflag=0
          return
        end if
c
c  If FU has the same sign as FB, then store the value of A in B.
c
        if(sign(one,fu).eq.sign(one,fb))then
          b=a
          fb=fa
        end if

      end if
c
c  Segment 2.  Rearrange points if necessary to ensure that ABS(FU)<ABS(FB).
c
      if(abs(fb).lt.abs(fu)) then
        a=u
        u=b
        b=a
        fa=fu
        fu=fb
        fb=fa
      end if
c
c  Segment 3.  Check to see if we can accept the current estimate because
c  the current change-in-sign interval (B,U) or (U,B) is very small.
c
      toler=two*epmach*abs(u)+epmach
      halfub=(b-u)/two
      sdel4=sdel3
      sdel3=sdel2
      sdel2=sdel1
      sdel1=abs(b-u)

      if(abs(halfub).le.toler) then
        iflag=-1
        a=u
        fa=fu
        return
      end if
c
c  Segment 4.  Compute a new approximate root, of the form U(new)=U(old)+STEP.
c  Methods availabe are linear interpolation, inverse quadratic interpolation,
c  and bisection.
c
      if(abs(fu).ge.abs(fa))then
        iflag=1
        step=halfub
        go to 10
      end if
c
c  If only two points are available, use linear interpolation.
c
      if(a.eq.b)then
        iflag=2
        s=fu/fa
        p=two*halfub*s
        q=one-s
c
c  If three points are available, try inverse quadratic 
c  interpolation.
c
      else

        iflag=3
        s=fu/fa
        q=fa/fb
        r=fu/fb
        p=s*(two*halfub*q*(q-r)-(u-a)*(r-one))
        q=(q-one)*(r-one)*(s-one)

      end if
c
c  Correct the signs of P and Q.
c
      if(p.gt.zero)then
        q=-q
      else
        p=-p
      end if
c
c  If P/Q is too large, use bisection instead.
c
      if((eight*sdel1.gt.sdel4)
     &  .or.(p.ge.onep5*abs(halfub*q)-abs(toler*q)))then
        iflag=1
        step=halfub
        go to 10
      end if

      step=p/q
c
c  Segment 5.  The value of STEP is known.  Update information.
c  The change in sign intevarl is now (A,B) or (B,A).
c
10    continue
      a=u
      fa=fu
      if(abs(step).le.toler)step=sign(toler,halfub)

      u=u+step

      return
      end
      subroutine setstp ( iwork, liw, lrw, rwork )

c*********************************************************************72
c
cc SETSTP computes the stepsize to be used by the Euler prediction step.
c
c  Discussion:
c
c    The formulas underlying the algorithm are:
c
c    ALFMIN = A minimal angle, ALFMIN=2*ARCCOS(1-EPMACH).
c
c    ALPHLC = Angle between last two tangents, value of ARCCOS(TL dot TC),
c             except that ALPHLC must be at least equal to ALFMIN.
c
c    HSEC   = Euclidean norm of the secant step, NORM2(XC-XF).
c
c    HSECL  = Euclidean norm of previous secant step, NORM2(XL-XC).
c
c    ABSNLC = ABS(SIN(.5*ALPHLC))
c
c    CURV   = Previous value of the curvature.
c
c    CURVN  = 2*ABSNLC/HSEC
c
c    CORDIS = Distance between predicted and corrected points.
c             Adjust CORDIS to lie between 0.01*HSEC and HSEC,
c           
c             But if CORDIS=0, meaning the predicted point was accepted,
c             set HTAN=HFACT*HSEC instead of using the first estimate
c             for HTAN.
c
c    Then
c
c      CURVX=CURVN+HSEC*(CURVN-CURV)/(HSEC+HSECL)
c
c    A simpler formula is used if we do not have data at two old points.
c
c    If (IWORK(10).GE.2) CURVX must be at least as large as the maximum
c    of 0.001 and 0.01/HSEC.
c
c    First estimate for stepsize (unless CORDIS=0):
c
c      HTAN=SQRT(2*QUAL*CORDIS/CURVX)
c
c    Adjusted value:
c  
c      HTAN=HTAN*(1.0+HTAN*(TC(IPC)-TL(IPC))/(2*HSEC*TC(IPC)))
c
c    Readjustments to the calculated value:
c
c    If stepsize reduction occurred during the correction of the previous
c    continuation point, HTAN is forced to be less than (HFACT-1)*HSEC/2.
c    The calculated HTAN is forced to lie between (HSEC/HFACT) and (HSEC*HFACT).
c    The calculated HTAN is also forced to lie between HMIN and HMAX.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 May 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision hundth
      double precision one
      double precision thouth
      double precision two
      double precision zero

      parameter (hundth=0.01D+00)
      parameter (one=1.0D+00)
      parameter (thouth=0.001D+00)
      parameter (two=2.0D+00)
      parameter (zero=0.0D+00)

      intrinsic abs
      intrinsic max
      intrinsic min
      intrinsic sin
      intrinsic sqrt

      integer liw
      integer lrw

      double precision curv
      double precision curvx
      double precision htan
      integer iwork(liw)
      double precision rwork(lrw)
      double precision temp

      rwork(11)=max(rwork(10),rwork(11))
c
c  Update estimates of curvature.
c
      curv=rwork(16)
      rwork(16)=two*abs(sin(rwork(11)/two))/rwork(21)
      if(curv.eq.zero)curv=rwork(16)

      if(rwork(22).eq.zero)then
        curvx=rwork(16) 
      else
        curvx=rwork(16)
     &    +rwork(21)*(rwork(16)-curv)/(rwork(21)+rwork(22))
      end if

      curvx=max( curvx, thouth)
      curvx=max( curvx, hundth/rwork(21) )
c
c  If the convergence distance was zero, then set the tentative next
c  stepsize to the maximum growth factor times the size of the secant step.
c
c  Otherwise, use the curvature estimate and other information to
c  compute an optimal step.
c
      if(rwork(15).eq.zero)then
        htan=rwork(20)*rwork(21)
      else
        temp=rwork(23)*rwork(15)
        temp=max(temp,hundth*rwork(21))
        temp=min(temp,rwork(21))
        htan=sqrt(two*temp/curvx)
      end if
c
c  Adjust the step to account for estimated curvature in the direction
c  of the parameter.
c
      if(iwork(18).gt.0)
     &  htan=min(htan,(rwork(20)-one)*rwork(21)/two)

      if(iwork(3).ne.1)then
        temp=one+(one-rwork(25)/rwork(24))*(htan/two)/rwork(21)
        htan=htan*temp
      end if
c
c  Enforce restrictions on growth and size of the step.
c
      htan=max(htan,rwork(21)/rwork(20))
      htan=min(htan,rwork(21)*rwork(20))
      htan=max(htan,rwork(3))
      htan=min(htan,rwork(4))
      rwork(5)=htan

      if(iwork(7).ge.2)then
        write(*,*)'SETSTP - Next stepsize HTAN = ',htan
      end if

      return
      end
      subroutine stajac ( df, fpar, fx, ierror, ipar, ipc, iwork, liw,
     &  lrw, nvar, rwork, wk, xr, slname )

c*********************************************************************72
c
cc STAJAC generates and factors the jacobian the very first time.
c
c  Discussion:
c
c    This routine is used if the user has requested the option to compute
c    the jacobian as rarely as possible.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 May 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      external df
      external fx
      external slname

      integer liw
      integer lrw
      integer nvar

      double precision dets
      double precision fpar(*)
      integer ierror
      integer ipar(*)
      integer ipc
      integer iwork(liw)
      integer iwrite
      integer job
      double precision rwork(lrw)
      double precision wk(nvar)
      double precision xr(nvar)

      iwrite=iwork(7)
c
c  If user is requesting that Jacobian be used as long as possible,
c  then go ahead, generate and factor the first one now.
c
      if(iwork(4).eq.2)then

        if(iwrite.ge.2)then
          write(*,*)'STAJAC - Generating initial jacobian.'
        end if

        job=2
        call slname(dets,fx,df,fpar,ierror,ipc,ipar,iwork,liw,
     &    job,nvar,rwork,lrw,xr,wk)
        rwork(17)=dets

        if(ierror.ne.0)then

          if(iwrite.ge.1)then
            write(*,*)'STAJAC - Serious error!'
            write(*,*)'  Could not factor initial jacobian.'
          end if

          return
        end if

      end if

      return
      end
      subroutine start ( df, fpar, fx, ierror, ipar, ipc, iwork, liw,
     &  lrw, nvar, rwork, wk, xc, xf, xr, slname )

c*********************************************************************72
c
cc START forces the starting point to satisfy the nonlinear system.
c
c  Discussion:
c
c    START makes sure that, on the first call, the input point XR 
c    satisfies the nonlinear equations.  If not, it calls call CORECT 
c    to use Newton iteration to improve the point.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 May 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision one
      parameter (one=1.0D+00)

      intrinsic abs
      external corect
      external df
      external fx
      external idamax
      external daxpy
      external dcopy
      external slname

      integer liw
      integer lrw
      integer nvar

      double precision fpar(*)
      integer icrit
      integer ierror
      integer imax
      integer ipar(*)
      integer ipc
      integer idamax
      integer iwork(liw)
      integer iwrite
      integer modsav
      double precision rwork(lrw)
      double precision wk(nvar)
      double precision xc(nvar)
      double precision xf(nvar)
      double precision xr(nvar)

      iwrite=iwork(7)

      if(iwrite.ge.2)then
        write(*,*)'START  - Checking the initial point.'
        write(*,*)'  Fixing variable number',ipc
      end if

      call dcopy(nvar,xr,1,xc,1)

      modsav=iwork(4)
      icrit=1

20    continue

      call dcopy(nvar,xc,1,xr,1)

      call corect(df,fpar,fx,ierror,ipc,ipar,iwork,
     &  nvar,rwork,wk,xr,lrw,liw,icrit,slname)

      iwork(25)=iwork(25)+iwork(28)
c
c  If an error occurred, then, if possible, retry with ICRIT=2.
c
      if(ierror.ne.0)then

        if(icrit.eq.1)then

          if(iwrite.ge.1)then
            write(*,*)'START -  Warning!'
            write(*,*)'  The first try to correct the starting point'
            write(*,*)'  has failed.  Correction of the starting '
            write(*,*)'  point will be retried.'
          end if

          icrit=2
          go to 20

        else

          icrit=1

        end if

        if(iwork(4).gt.0)then
          iwork(4)=iwork(4)-1
          ierror=0

          if(iwrite.ge.1)then
            write(*,*)'START  - Retry starting point correction.'
            write(*,*)'  Set Newton option IWORK(4)=',iwork(4)
          end if

          go to 20
        end if

      end if

      iwork(4)=modsav

      if(ierror.ne.0)then

        if(iwrite.ge.1)then
          write(*,*)'START - Fatal error!'
          write(*,*)'  Correction of the starting point failed.'
        end if

        return

      end if
c
c  If we were able to correct the starting point, then
c  record necessary data.
c
      call daxpy(nvar,-one,xr,1,xc,1)
      imax=idamax(nvar,xc,1)
      rwork(15)=abs(xc(imax))
      call dcopy(nvar,xr,1,xc,1)
      call dcopy(nvar,xr,1,xf,1)
      call coqual(iwork,liw,rwork,lrw)

      if(iwrite.ge.3)then
        write(*,*)
     &    'START - Corrector convergence quality factor=',rwork(23)
      end if

      rwork(14)=rwork(13)
      iwork(27)=iwork(27)+1
      iwork(10)=1
      iwork(1)=1

      return
      end
      subroutine tangnt ( detsn, fx, df, fpar, ierror, ip, ipar, iwork,
     &  nvar, rwork, tan, xr, liw, lrw, slname )

c*********************************************************************72
c
cc TANGNT computes a tangent vector to the solution curve F(X)=0.  
c
c  Discussion:
c
c    The length of the tangent vector in the Euclidean norm will be 1.  
c    There are two such tangent vectors at each point, one the negative 
c    of the other.  This routine produces one such vector.  The choice of 
c    which is appropriate must be made outside this routine.
c
c    The tangent vector TAN is the solution of the linear system
c
c      DFA(X,IP)*TAN = E(NVAR)
c
c    Here E(I) denotes the I-th basis vector, that is, the vector with 
c    a 1 in the I-th position and 0 elsewhere.  DFA(X,IP) is the NVAR by 
c    NVAR matrix whose first NVAR-1 rows are the jacobian of FX evaluated 
c    at X, and whose last row is the transpose of E(IP).
c
c    After computation, the tangent vector is normalized so that it has
c    a Euclidean norm of 1.  The adjustment of sign is performed outside
c    this routine.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, INTEGER IERROR, error flag.
c    0, no error occurred.
c    1, data or storage error.
c    2, error occurred in user jacobian routine.
c    3, an error occurred in the solver.
c    6, the computed tangent vector is entirely zero.
c
      implicit none

      double precision one
      double precision zero

      parameter (one=1.0D+00)
      parameter (zero=0.0D+00)

      external fx
      external df
      external slname
      external dnrm2

      integer liw
      integer lrw
      integer nvar

      double precision detsn
      double precision fpar(*)
      integer i
      integer ierror
      integer ip
      integer ipar(*)
      integer iwork(liw)
      integer iwrite
      integer job
      double precision rwork(lrw)
      double precision dnrm2
      double precision tan(nvar)
      double precision tnorm
      double precision xr(nvar)

      iwrite=iwork(7)
c
c
c  Set the right hand side of the linear system to be solved.
c
      do i=1,nvar
        tan(i)=zero
      enddo

      tan(nvar)=one
c
c  Call the user-specified solver.
c
      job=0
      if(iwork(4).eq.2)job=1
      call slname(detsn,fx,df,fpar,ierror,ip,ipar,iwork,liw,
     &  job,nvar,rwork,lrw,xr,tan)
 
      if(ierror.ne.0)then
        if(iwrite.ge.1)then
          write(*,*)'TANGNT - Warning!'
          write(*,*)
     &      'The linear solver returned error flag IERROR=',ierror
        end if
        return
      end if
c
c  Normalize the tangent vector.
c
      tnorm=dnrm2(nvar,tan,1)

      if(tnorm.eq.zero)then

        ierror=6

        if(iwrite.ge.1)then
          write(*,*)'TANGNT - Warning!'
          write(*,*)'The computed tangent has zero norm!'
        end if

      else

        do i=1,nvar
          tan(i)=tan(i)/tnorm
        enddo

      end if
 
      return
      end
      subroutine tanpar ( df, fpar, fx, ierror, ipar, iwork, liw,
     &  lrw, nvar, rwork, slname, tc, wk, xc, xf, xr )

c*********************************************************************72
c
cc TANPAR controls the computation of the next tangent vector.
c
c  Discussion:
c
c    This computation also is used to determine the next continuation 
c    parameter index.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 May 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision one
      double precision zero

      parameter (one=1.0D+00)
      parameter (zero=0.0D+00)

      intrinsic abs
      intrinsic atan2
      external df
      external fx
      external idamax
      external dcopy
      external ddot
      intrinsic sign
      external slname
      intrinsic sqrt
      external dscal
      external tangnt

      integer liw
      integer lrw
      integer nvar

      double precision atcipc
      double precision atcjpc
      double precision fpar(*)
      integer idamax
      integer ierror
      integer ipar(*)
      integer ipc
      integer ipl
      integer jpc
      integer itemp
      integer iwork(liw)
      integer iwrite
      double precision rwork(lrw)
      double precision ddot
      double precision tc(nvar)
      double precision tcipc
      double precision tcos
      double precision temp
      double precision tlipc
      double precision tsin
      double precision wk(nvar)
      double precision xc(nvar)
      double precision xf(nvar)
      double precision xr(nvar)

      iwrite=iwork(7)
c
c  Retrieve the previous value of the continuation parameter as IPL.
c  Move the old tangent to WK.
c  Move the sign of the previous determinant to RWORK(29).
c
      ipl=iwork(2)
      call dcopy(nvar,tc,1,wk,1)
      rwork(29)=rwork(17)
c
c  Compute the tangent for the current point and store in TC.
c
      call tangnt(rwork(17),fx,df,fpar,ierror,iwork(2),ipar,iwork,nvar,
     &  rwork,tc,xf,liw,lrw,slname)
 
      if(ierror.ne.0)then

        if(iwrite.ge.1)then
          write(*,*)'TANPAR - Serious error!'
          write(*,*)'  The tangent calculation failed.'
        end if

        return
      end if
c
c  Find the two largest components of tangent vector, which are our
c  candidates for continuation directions, unless the user demands
c  a fixed parameterization.
c
      if(iwork(3).ne.1)then
        ipc=idamax(nvar,tc,1)
        temp=tc(ipc)
        tc(ipc)=zero
        jpc=idamax(nvar,tc,1)
        if(jpc.eq.0)jpc=ipc
        tc(ipc)=temp
      else
        ipc=iwork(2)
        jpc=ipc
      end if
c
c  Adjust the sign of the tangent vector TC.  To do so, compare the sign
c  of the IPL-th component with the sign of the IPL-th component of
c  the secant vector XF-XC.  If a target or limit point, XR, has been
c  computed, use XF-XR instead of XF-XC.  And on the first step, we'll
c  have to use the user's input direction for our comparison.
c
      if(iwork(10).le.1)then
        temp=rwork(6)
      elseif(iwork(1).eq.3.or.iwork(1).eq.4)then
        temp=xf(ipl)-xr(ipl)
      else
        temp=xf(ipl)-xc(ipl)
      end if

      if(sign(one,tc(ipl)).ne.sign(one,temp))then
        call dscal(nvar,-one,tc,1)
        rwork(17)=-rwork(17)
      end if
c
c  Unless we are computing a starting point, record the new state.
c
      if(iwork(10).gt.1)iwork(10)=3
c
c  The index of the entry of largest magnitude in the tangent vector will
c  be used for the local parameterization, unless the user has ordered
c  us to stick to a given parameter, or if a limit point in the first
c  choice appears to be approaching.
c
c  To check this, we compare TC(IPC) and the second largest component, TC(JPC).
c  If TC(JPC) is no less than 0.1 of TC(IPC), and
c     TC(JPC) is larger than TL(JPC), and
c     TC(IPC) is smaller than TL(IPC), then
c  we suspect a limit point may be coming in the IPC direction and we switch
c  our choice to JPC.
c
      if(iwork(3).ne.1)then

        atcipc=abs(tc(ipc))
        atcjpc=abs(tc(jpc))

        if(jpc.ne.ipc.and.iwork(10).gt.1)then
          tlipc=wk(ipc)
          tcipc=tc(ipc)
          temp=abs(wk(jpc))

          if((sign(one,tcipc).eq.sign(one,tlipc))
     &      .and.(atcipc.lt.abs(tlipc))
     &      .and.(atcjpc.ge.max(0.1D+00*atcipc,temp)))then

            if(iwrite.ge.3)then
              write(*,*)'TANPAR - A limit point may be coming'
              write(*,*)'in the preferable index ',ipc
              write(*,*)'We''ll try second-best index.'
            end if

            itemp=ipc
            ipc=jpc
            jpc=itemp
          end if

        end if

        if(iwrite.ge.3)then
          write(*,*)
     &      'TANPAR - Continuation index: First choice= ',ipc
          if(jpc.ne.ipc)write(*,*)
     &      '                             Second choice=',jpc
        end if
      end if
c
c  Record the values of the IPC-th component of the new tangent vector TC,
c  as well as the IPC-th component of the old tangent vector TL.
c  Set the sign of the determinant.
c  Record the value of the LIM-th compoment of the new tangent vector,
c  if a limit point check is being done.
c
      iwork(2)=ipc
      iwork(12)=jpc
      rwork(24)=tc(ipc)
      rwork(25)=wk(ipc)
      rwork(6)=rwork(17)
      if(iwork(6).gt.0)then
        rwork(27)=rwork(26)
        rwork(26)=tc(iwork(6))
        if(iwrite.ge.2)write(*,*)
     &    'TANPAR - Tangent vector has limit component =',rwork(26)
      end if
c
c  Compute the angle between the old and new tangents.
c
      if(iwork(10).le.1)return
      tcos=ddot(nvar,wk,1,tc,1)
      if(tcos.gt.one)tcos=one
      if(tcos.lt.(-one))tcos=-one
      tsin=sqrt(one-tcos*tcos)
      rwork(11)=atan2(tsin,tcos)
c
c  Note possible bifurcation point, if determinant has changed sign.
c
c     if(rwork(17).ne.rwork(29))then
c       if(iwrite.ge.1)write(*,*)
c    &    'TANPAR - Possible bifurcation point detected.'
c     end if

      return
      end
      subroutine target ( df, fpar, fx, ierror, ifound, ipar, iwork, 
     &  liw, lrw, nvar, rwork, slname, wk, xc, xf, xr )

c*********************************************************************72
c
cc TARGET controls the computation of a target point.
c
c  Discussion:
c
c    The calling code has detected that between the points XC and XF must lie a
c    point XR whose IT-th component has the value XIT desired.  TARGET
c    must compute that point.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, INTEGER IFOUND.  Reports whether a target point was detected.
c    0, no target point was detected.
c    1, A target point was detected.
c
      implicit none

      external corect
      external df
      external fx
      external daxpy
      external dcopy
      external slname
      external dnrm2
      external dscal

      double precision one

      parameter (one=1.0D+00)

      integer liw
      integer lrw
      integer nvar

      double precision fpar(*)
      integer icrit
      integer ierror
      integer ifound
      integer ipar(*)
      integer it
      integer iwork(liw)
      integer iwrite
      integer modsav
      double precision rwork(lrw)
      double precision dnrm2
      double precision skale
      double precision wk(nvar)
      double precision xlow
      double precision xup
      double precision xc(nvar)
      double precision xf(nvar)
      double precision xr(nvar)

      iwrite=iwork(7)
      ifound=0
      it=iwork(5)
      if(it.eq.0)return
      if(iwork(10).le.1)return
      if( (iwork(1).eq.3) .and.
     &    (rwork(7).eq.rwork(28)) .and.
     &    (it.eq.iwork(11)) )return

      xlow=xc(it)
      xup=xf(it)
      if(rwork(7).lt.xlow.and.rwork(7).lt.xup)return
      if(rwork(7).gt.xlow.and.rwork(7).gt.xup)return
      ifound=1
      if(iwrite.ge.2)write(*,*)
     &  'TARGET - A target point has been detected.'
c
c  Approximate the target point using the bracketing solutions.
c
      modsav=iwork(4)
10    continue

      if(xlow.ne.xup)then
        skale=(rwork(7)-xlow)/(xup-xlow)
      else
        skale=one
      end if

      call dcopy(nvar,xf,1,xr,1)
      call dscal(nvar,skale,xr,1)
      skale=one-skale
      call daxpy(nvar,skale,xc,1,xr,1)
      xr(it)=rwork(7)
c
c  Call CORECT to compute the exact target point, holding index IT fixed.
c
      icrit=0

      call corect(df,fpar,fx,ierror,it,ipar,iwork,
     &  nvar,rwork,wk,xr,lrw,liw,icrit,slname)

      iwork(24)=iwork(24)+iwork(28)
      if(ierror.ne.0.and.iwork(4).gt.0)then
        ierror=0
        iwork(4)=iwork(4)-1
        if(iwrite.ge.1)write(*,*)
     &    'TARGET - Retry target computation with IWORK(4)=',iwork(4)
        go to 10
      end if

      iwork(4)=modsav
      if(ierror.ne.0)then
        if(iwrite.ge.1)then
          write(*,*)'TARGET - Target point calculation failed.'
        end if
        return
      end if
c
c  Record the values of IT and XIT, and compute the arclength to the target
c  point.
c
      iwork(1)=3
      iwork(11)=it
      rwork(28)=rwork(7)
      call dcopy(nvar,xr,1,wk,1)
      call daxpy(nvar,-one,xc,1,wk,1)
      rwork(14)=rwork(12)+dnrm2(nvar,wk,1)
      iwork(27)=iwork(27)+1
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
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
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
      subroutine trystp ( df, fpar, fx, ierror, ipar, iwork, liw, 
     &  lrw, nvar, rwork, slname, tc, wk, xf, xr )

c*********************************************************************72
c
cc TRYSTP tries to carry out a continuation step.  
c
c  Discussion:
c
c    Given a point X, a tangent vector TAN, and a stepsize H, TRYSTP 
c    estimates the value of the next point on the curve as X + H*TAN, 
c    and then uses Newton iteration to correct this estimate.  
c
c    The routine tries various fixes if the Newton iteration fails to converge, 
c    or does not converge rapidly enough.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 May 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision one
      parameter (one=1.0D+00)

      external corect
      external df
      external fx
      external daxpy
      external dcopy
      external slname

      integer liw
      integer lrw
      integer nvar

      double precision fpar(*)
      integer i
      integer icrit
      integer ierror
      integer ihold
      integer ipar(*)
      integer ipc
      integer iwork(liw)
      integer iwrite
      integer jpc
      integer modsav
      double precision rwork(lrw)
      double precision tc(nvar)
      double precision temp
      double precision wk(nvar)
      double precision xf(nvar)
      double precision xr(nvar)

      iwork(18)=0
      ihold=iwork(2)
      ipc=iwork(2)
      iwrite=iwork(7)
      jpc=iwork(12)
      modsav=iwork(4)
c
c  Set initial guess for next point to XR = XF + RWORK(5) * TC
c
10    continue

      if(iwrite.ge.2)then
        write(*,*)'TRYSTP - Predictor using stepsize ',rwork(5)
      end if

      call dcopy(nvar,xf,1,xr,1)
      call daxpy(nvar,rwork(5),tc,1,xr,1)
c
c  If you return to this statement from a later statement, we will
c  NOT be changing the current iterate, but rather, altering the
c  parameter index, trying to "save" the iteration.
c
20    continue
 
      if(iwrite.ge.2.and.iwork(3).eq.0)then
        write(*,*)'TRYSTP - Corrector is fixing index ',ihold
      end if
 
      icrit=0

      call corect(df,fpar,fx,ierror,ihold,ipar,iwork,
     &  nvar,rwork,wk,xr,lrw,liw,icrit,slname)

      iwork(25)=iwork(25)+iwork(28)
c
c  If the corrector succeeds, we're done.  However, before we 
c  return, we may have to restore the default values of the 
c  continuation parameter and the Newton correction algorithm, 
c  because these may have been temporarily altered in a desperate 
c  situation!
c
      if(ierror.eq.0)then
        if(iwork(3).ne.1)iwork(2)=ihold
        iwork(4)=modsav
        return
      end if
c
c  Only VERY fatal errors should abort the correction process.
c  Abort this process only after kicking and screaming.
c
      if(ierror.eq.2)then

        if(iwrite.ge.1)then
          write(*,*)'TRYSTP - Fatal error during Newton correction'
          write(*,*)'of a continutation point!'
        end if

        return
      end if
c
c  We reach this point if the corrector has failed to converge, but
c  not in a disastrous way.  We have two options: vary the parameter
c  held fixed, or reduce the stepsize.  We prefer to vary the parameter,
c  since that would allow us to take the same "healthy" step.
c
c  We will try varying the parameter if:
c
c    There is a nonzero second choice parameter (JPC); and
c    The corrector was not already using the second choice parameter; and
c    The user did not request via IWORK(3) that we always use a particular
c      parameter.
c
      if(iwork(3).ne.1.and.ihold.ne.jpc.and.jpc.ne.0)then
        ihold=jpc
        if(ierror.eq.5)go to 20
        go to 10
      end if
c
c  If JPC fails as the index held fixed during the corrector iteration,
c  return to using IPC.
c
      ihold=ipc
c
c  If we are using some modified form of Newton's method, and
c  we have reached the minimum stepsize, or
c  we have had two failures in a row on this step, then
c  we will try a better form of Newton's method.
c
      if(iwork(4).gt.0.and.
     &  (rwork(5).lt.rwork(20)*rwork(3).or.iwork(18).ge.2) )then
        iwork(4)=iwork(4)-1

        if(iwrite.ge.1)then
          write(*,*)
     &      'TRYSTP - Retrying step with IWORK(4)=',iwork(4)
        end if

        if(ierror.ne.5)go to 10
        go to 20

      end if
c
c  No convergence, so reduce the stepsize.  At this point, if the
c  stepsize falls below the user-specified minimum, we have to quit.
c
      if(rwork(5).lt.rwork(20)*rwork(3))then

        ierror=4

        if(iwrite.ge.1)then
          temp=rwork(5)/rwork(20)
          write(*,*)'TRYSTP - Warning!'
          write(*,*)'  The predictor stepsize fell below minimum.'
          write(*,*)'  Current step is now ',temp
          write(*,*)'  Minimum step is     ',rwork(3)
        end if

        return
      end if

      rwork(5)=rwork(5)/rwork(20)
      iwork(18)=iwork(18)+1
      if(ierror.ne.5)go to 10
c
c  We're reducing the stepsize, but the corrector iteration was converging,
c  though slowly.  We'll try to salvage the work we had done by using
c  as our new predicted point the linear interpolant between our current
c  accepted point and the corrector iterate that we gave up on.
c
      temp=one/rwork(20)

      do i=1,nvar
        xr(i)=temp*xr(i)+(one-temp)*xf(i)
      enddo

      go to 20
      end
      subroutine update ( iwork, liw, lrw, nvar, rwork, tc, wk, xc, 
     &  xf, xr )

c*********************************************************************72
c
cc UPDATE updates information after a successful continuation step.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 May 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision one
      double precision zero

      parameter (one=1.0D+00)
      parameter (zero=0.0D+00)

      integer liw
      integer lrw
      integer nvar

      integer i
      integer iwork(liw)
      integer iwrite
      double precision rwork(lrw)
      double precision dnrm2
      double precision tc(nvar)
      double precision temp
      double precision wk(nvar)
      double precision xc(nvar)
      double precision xf(nvar)
      double precision xr(nvar)

      iwrite=iwork(7)

      iwork(1)=2
c
c  Note that IWORK(2) may be set to IHOLD here.  This simply reflects
c  the fact that the corrector iteration may have failed with the 
c  preferable index, and had to try the "second-best" parameter.   
c  Hence, our initial expectation that the parameter would be 
c  determined by the maximum entry of the tangent is tempered by our
c  experience with the corrector iteration.
c
      iwork(10)=2
      iwork(26)=iwork(26)+iwork(18)
      iwork(27)=iwork(27)+1
      rwork(22)=rwork(21)
      if(iwrite.ge.1.and.iwork(18).gt.0)write(*,*)
     &  'UPDATE - Predictor stepsize was reduced ',iwork(18),' times.'
      call dcopy(nvar,xr,1,wk,1)
      call daxpy(nvar,-one,xf,1,wk,1)
      rwork(21)=dnrm2(nvar,wk,1)
      rwork(12)=rwork(13)
      rwork(13)=rwork(12)+rwork(21)
      rwork(14)=rwork(13)
c
c  Compute and store "CORDIS", the maximum-norm of the total 
c  correction, MaxNorm(XF+HTAN*TC - XR).
c
      rwork(15)=zero

      if(iwork(28).ne.0)then

        do i=1,nvar
          temp=xr(i)-(xf(i)+rwork(5)*tc(i))
          if(abs(temp).gt.rwork(15))rwork(15)=abs(temp)
        enddo

      end if
c
c  Update XC and XF, so that XC has the older continuation point, and
c  XF has the most recent one.
c
      call dcopy(nvar,xf,1,xc,1)
      call dcopy(nvar,xr,1,xf,1)

      return
      end
