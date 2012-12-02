      program main

c*********************************************************************72
c
cc MAIN tests the DRIV library.
c
c  Discussion:
c
C    A sample problem for DRIV.  The three sample programs below
C    call SDRIV1, SDRIV2, and SDRIV3, which are the easy, moderately
C    involved, and advanced interfaces to SDRIV.  Use SDRIV1 first.
C    If you get unsatisfactory results, advance to a more complicated
C    interface.
C
C    These examples do not attempt to demonstrate all possible uses
C    of DRIV.  They are just examples that run properly.
C
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DRIV_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the DRIV library.' 

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DRIV_PRB'
      write ( *, '(a)' ) '  Normal end of execution.' 

      STOP
      END
      SUBROUTINE TEST01

c*********************************************************************72
c
C  Sample calling for SDRIV1, the simplest interface to SDRIV
C
C  N and LENW are set up as parameters, because other quantities
C  require them as DIMENSIONS.  Moreover, LENW is based on N.
C
      PARAMETER (N=2)
      PARAMETER (LENW=N*N+10*N+225)
      DIMENSION WORK(LENW)
      DIMENSION Y(N)
C
C  Set the initial conditions
C
      T=0.0
      Y(1)=2.0
      Y(2)=0.0
C
C  Initialize MSTATE so that SDRIV1 knows we are beginning a
C  new problem.  On output, MSTATE is an error flag we must check.
C
      MSTATE=1
C
C  Set the error tolerance
C
      EPS=0.0001
      WRITE(*,*)' '
      WRITE(*,*)'TEST01'
      WRITE(*,*)'Test SDRIV1, the simple interface.'
      WRITE(*,*)'Solving Van der Pol equation.'
      WRITE(*,*)'Exact solution has Y(1)=0 for points 3, 5, 7,...'
      WRITE(*,*)' '
      DO 20 I=0,12
        IF(I.EQ.0)THEN
          TOUT=T
        ELSEIF(I.EQ.1)THEN
          TOUT=1.39283880203
        ELSE
          TOUT=TOUT+2.214773875
          ENDIF
        CALL SDRIV1(N,T,Y,TOUT,MSTATE,EPS,WORK,LENW)
        IF(MSTATE.GT.2)THEN
          WRITE(*,*)'Error return MSTATE=',MSTATE
          RETURN
          ENDIF
        WRITE(*,'(1X,3G14.6)')TOUT,(Y(J),J=1,N)
20      CONTINUE
      RETURN
      END
      SUBROUTINE TEST02

c*********************************************************************72
c
C  Sample calling for SDRIV2, the moderately-involved interface.
C
C  N, LENIW and LENW are set up as parameters, because other quantities
C  require them as DIMENSIONS.  Moreover, LENIW and LENW are based on N.
C
      PARAMETER (N=2)
      PARAMETER (LENIW=N+21)
      PARAMETER (LENW=N*N+9*N+204)
C
      EXTERNAL F
      EXTERNAL G
C
      DIMENSION IWORK(LENIW)
      DIMENSION WORK(LENW)
      DIMENSION Y(N)
C
C  Set the initial conditions
C
      T=0.0
      Y(1)=2.0
      Y(2)=0.0
C
C  Initialize MSTATE so that SDRIV2 knows we are beginning a
C  new problem.  On output, MSTATE is an error flag we must check.
C
      MSTATE=1
C
C  Set the error tolerance
C
      EPS=0.0001
C
C  No roots of equations are desired
C
      NROOT=0
C
C  Set a sort of scale value, roughly the smallest physically meaningful
C  value of the solution.  EWT=0 sets pure relative error tolerance,
C  EWT>0 sets a mixed absolute/relative error tolerance.
C
      EWT=EPS
C
C  Set integration method flag.
C  1=Adams method, 2=Gear method, 3=program switches dynamically between
C  both Adams and Gear methods.
C
      MINT=2
      WRITE(*,*)' '
      WRITE(*,*)'TEST02'
      WRITE(*,*)'Test SDRIV2, moderately-involved interface.'
      WRITE(*,*)'Solving Van der Pol equation.'
      WRITE(*,*)'Exact solution has Y(1)=0 for points 3, 5, 7,...'
      WRITE(*,*)' '
      DO 20 I=0,12
        IF(I.EQ.0)THEN
          TOUT=T
        ELSEIF(I.EQ.1)THEN
          TOUT=1.39283880203
        ELSE
          TOUT=TOUT+2.214773875
          ENDIF
        CALL SDRIV2(N,T,Y,F,TOUT,MSTATE,NROOT,EPS,EWT,
     *  MINT,WORK,LENW,IWORK,LENIW,G)
        IF(MSTATE.GT.2)THEN
          WRITE(*,*)'Error return, MSTATE=',MSTATE
          RETURN
          ENDIF
        WRITE(*,'(1X,3G14.6)')TOUT,(Y(J),J=1,N)
20      CONTINUE
      RETURN
      END
      SUBROUTINE TEST03

c*********************************************************************72
c
C  Sample calling for SDRIV3, the advanced interface to SDRIV
C
C  N, LENIW and LENW are set up as parameters, because other quantities
C  require them as DIMENSIONS.  Moreover, LENIW and LENW are based on N.
C
      PARAMETER (N=2)
      PARAMETER (LENIW=N+21)
      PARAMETER (LENW=N*N+9*N+204)
C
      EXTERNAL F
      EXTERNAL G
C
      DIMENSION IWORK(LENIW)
      DIMENSION WORK(LENW)
      DIMENSION Y(N)
C
C  Set the initial conditions
C
      T=0.0
      Y(1)=2.0
      Y(2)=0.0
C
C  Initialize MSTATE so that SDRIV3 knows we are beginning a
C  new problem.  On output, MSTATE is an error flag we must check.
C
      MSTATE=1
C
C  Set the error tolerance
C
      EPS=0.0001
C
C  No roots of equations are desired
C
      NROOT=0
C
C  Set a sort of scale value, roughly the smallest physically meaningful
C  value of the solution.  EWT=0 sets pure relative error tolerance,
C  EWT>0 sets a mixed absolute/relative error tolerance.
C
      EWT=EPS
C
C  Set integration method flag.
C  1=Adams method, 2=Gear method, 3=program switches dynamically between
C  both Adams and Gear methods.
C
      MINT=2
C
C  Set NTASK to allow program the freedom to interpolate past the
C  target value.
C
      NTASK=1
C
C  Set IERROR, the error control indicator
C
      IERROR=3
C
C  Set MITER, the iteration method indicator
C  If MITER is not 0 or 2, we have to supply some routines like
C  USERS and JACOBN.
C
      MITER=2
C 
C  Set IMPL, the implicit method indicator
C  Setting IMPL=0 means we are solving an explicit system, Y'=F.
C  IMPL=1 means we are solving A*Y'=F.  IMPL=2 means we are solving
C  a mixed differential/algebraic system.
C
      IMPL=0
C
C  ML, MU are the lower and upper bandwidths of the jacobian matrix
C  for options MITER=4, 5.  We can ignore these.
C
      ML=0
      MU=0
C
C  MXORD is the maximum order desired.  Can be up to 12 for Adams methods,
C  and up to 5 for Gear.
C
      MXORD=5
C
C  HMAX is the largest single stepsize we will allow
C
      HMAX=0.2
C
C  NDE is the number of differential equations, which is equal to N
C  except for the special problem IMPL=2
C
      NDE=N
C
C  MXSTEP is the maximum number of internal steps allowed on one call.
C
      MXSTEP=200
      WRITE(*,*)' '
      WRITE(*,*)'TEST03'
      WRITE(*,*)'Test SDRIV3, advanced interface.'
      WRITE(*,*)'Solving Van der Pol equation.'
      WRITE(*,*)'Exact solution has Y(1)=0 for points 3, 5, 7,...'
      WRITE(*,*)' '
      DO 20 I=0,12
        IF(I.EQ.0)THEN
          TOUT=T
        ELSEIF(I.EQ.1)THEN
          TOUT=1.39283880203
        ELSE
          TOUT=TOUT+2.214773875
          ENDIF
        CALL SDRIV3(N,T,Y,F,MSTATE,TOUT,NTASK,NROOT,EPS,EWT,IERROR,
     *  MINT,MITER,IMPL,ML,MU,MXORD,HMAX,WORK,LENW,IWORK,LENIW,JACOBN,
     *  FA,NDE,MXSTEP,G)
        IF(MSTATE.GT.2)THEN
          WRITE(*,*)'Error return, MSTATE=',MSTATE
          RETURN
          ENDIF
        WRITE(*,'(1X,3G14.6)')TOUT,(Y(J),J=1,N)
20      CONTINUE
      RETURN
      END
      SUBROUTINE F(N,T,Y,YDOT)

c*********************************************************************72
c
C  Given input of T and Y, compute the right hand sides of the
C  ODE system in YDOT.
C
      REAL Y(*),YDOT(*)
      YDOT(1)=Y(2)
      YDOT(2)=3.0*(1.0-Y(1)*Y(1))*Y(2) - Y(1)
      END
      REAL FUNCTION G(N,T,Y,IROOT)

c*********************************************************************72
c
C  Dummy function.  You only need to put a real routine here for
C  uses of SDRIV2 and SDRIV3 where rootfinding is required.
C
      REAL Y(*)
      G=0.0
      RETURN
      END
      SUBROUTINE USERS

c*********************************************************************72
c
C  Dummy routine.  You only need to put a real routine here for
C  certain advanced uses of the interface routine SDRIV3
C
      RETURN
      END
      SUBROUTINE FA(N,T,Y,A,MATDIM,ML,MU,NDE)

c*********************************************************************72
c
C  Dummy routine, only required for advanced uses of SDRIV3.
C
      REAL Y(*),A(MATDIM,*)
      RETURN
      END
      SUBROUTINE JACOBN(N,T,Y,DFDY,MATDIM,ML,MU)

c*********************************************************************72
c
C  Dummy routine, only required for advanced uses of SDRIV3.
C
      REAL Y(*),DFDY(MATDIM,*)
      RETURN
      END
      SUBROUTINE TEST04

c*********************************************************************72
c
C  Sample calling for CDRIV2, the moderately-involved interface.
C
C  N, LENIW and LENW are set up as parameters, because other quantities
C  require them as DIMENSIONS.  Moreover, LENIW and LENW are based on N.
C
      PARAMETER (N=2)
      PARAMETER (LENIW=N+21)
      PARAMETER (LENW=N*N+16*N+204)
      PARAMETER (PI=3.14159265)

      EXTERNAL FC
      EXTERNAL GC

      DIMENSION IWORK(LENIW)
      COMPLEX   WORK(LENW)
      COMPLEX   Y(N)
C
C  Set the initial conditions
C
      T=0.0
      Y(1)=CMPLX(1.0,0.0)
      Y(2)=CMPLX(0.0,1.0)
C
C  Initialize MSTATE so that SDRIV2 knows we are beginning a
C  new problem.  On output, MSTATE is an error flag we must check.
C
      MSTATE=1
C
C  Set the error tolerance
C
      EPS=0.00001
C
C  No roots of equations are desired
C
      NROOT=0
C
C  Set a sort of scale value, roughly the smallest physically meaningful
C  value of the solution.  EWT=0 sets pure relative error tolerance,
C  EWT>0 sets a mixed absolute/relative error tolerance.
C
      EWT=EPS
C
C  Set integration method flag.
C  1=Adams method, 2=Gear method, 3=program switches dynamically between
C  both Adams and Gear methods.
C
      MINT=2
      WRITE(*,*)' '
      WRITE(*,*)'TEST04'
      WRITE(*,*)'Test CDRIV2, moderately-involved interface.'
      WRITE(*,*)'Solving sine/cosine equation.'
      WRITE(*,*)'Y(1)''= Y(2)'
      WRITE(*,*)'Y(2)''=-Y(1)'
      WRITE(*,*)' '

      NSTEP=12
      WRITE(*,*)'       T      '//
     *          '   REAL(Y(1)) '//
     *          '  AIMAG(Y(1)) '//
     *          '   REAL(Y(2)) '//
     *          '  AIMAG(Y(2)) '
      WRITE(*,*)' '

      DO 20 I=0,NSTEP
        TOUT=(I*2.0*PI)/REAL(NSTEP)

        CALL CDRIV2(N,T,Y,FC,TOUT,MSTATE,NROOT,EPS,EWT,
     *  MINT,WORK,LENW,IWORK,LENIW,GC)
        IF(MSTATE.GT.2)THEN
          WRITE(*,*)'Error return, MSTATE=',MSTATE
          RETURN
          ENDIF
        WRITE(*,'(1X,5G14.6)')TOUT,(Y(J),J=1,N)
20      CONTINUE
      RETURN
      END
      SUBROUTINE FC(N,T,Y,YDOT)

c*********************************************************************72
c
C  Given input of T and Y, compute the right hand sides of the
C  ODE system in YDOT.
C
      REAL T
      COMPLEX Y(*),YDOT(*)
C
      YDOT(1)=Y(2)
      YDOT(2)=-Y(1)
      END
      REAL FUNCTION GC(N,T,Y,IROOT)

c*********************************************************************72
c
C  Dummy function.  You only need to put a real routine here for
C  uses of CDRIV2 and CDRIV3 where rootfinding is required.
C
      COMPLEX Y(*)
      GC=0.0
      RETURN
      END
