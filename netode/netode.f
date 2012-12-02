c
c  Netode.f  21 April 1995
c
      program main

c*********************************************************************72
c
cc MAIN is the main program for NETODE.
c
c  Discussion:
c
c    NETODE is a two fluid network ODE problem solved using LSODI.
c
c    We assume a network of nodes and links in which a
c    two-phase (gas and liquid) fluid is flowing, with the following 
c    properties:
c
c    At a node:
c
c    Q - heat source 
c    X - x coordinate
c    Y - y coordinate
c    Z - z coordinate
c    VOL - volume
c
c    At a link:
c
c    CRS - cross sectional area of link
c    VOL - volume
c
c    For the fluid:
c
c    FRICGW - gas/wall friction coefficient
c    FRICGL - gas/liquid friction coefficient
c    FRICLW - liquid/wall friction coefficient
c    GRAVIT - gravitational density
c    RHCRIT - critical density
c    SCRIT -  critical Sigma
c    WCRIT -  critical W
c
c    We are interested in the following variables
c
c    At a node:
c
c    ALPHA -  the void fraction
c    PRESS -  the pressure
c    HGAS  -  the enthalpy stored in the gas state
c    HLIQ  -  the enthalpy stored in the liquid state
c
c    At a link:
c
c    VGAS  -  the gas velocity
c    VLIQ  -  the liquid velocity
c
c    Derived quantities include:
c
c    RHGAS -  gas density, determined as a function of HGAS and PRESS
c    RHLIQ -  liquid density, a function of HLIQ and PRESS
c    ALPHA -  void fraction evaluated on links
c    GRAV  -  potential gravitational energy at a node
c             GRAV(node) = G*Y(node)
c    HGAS  -  HGAS evaluated at a link.
c    HLIQ  -  HLIQ evaluated at a link
c    PRSLNK - PRESS on the link 
c
c    This program can be run interactively or as a batch job. It is
c    recommended that at the very first run of a problem, run
c    this program interactively to generate the Definition file of the
c    problem.
c
c    Step by step running instructions and sample input:
c
c    Model problem: 
c    Two Node Separation Problem 
c    (II-4-2, EPRI NP-3657, Project 1741-2, 1984, T.A.Porsching)
c---------------------------------------------------------------------
c  the very first run
c
c  $run netode
c  $ IS THE PROBLEM DEFINITION IN A FILE?
c  $n
c  $ GONET:
c  $ GET NETODE PROBLEM DEFINITION.
c  $ ENTER A TITLE FOR THE PROBLEM.
c  $two node separation
c  $ ENTER THE NUMBER OF NODES.
c  $ BETWEEN 2 AND 15
c  $2
c  $ ENTER x, y, z, HEAT SOURCE AND VOLUME'
c  $ FOR EACH NODE.
c  $0.0 0.0 0.0 0.0 1.0
c  $0.0 0.0 1.0 0.0 1.0
c  $ DEFINE LINKS BETWEEN PAIRS OF NODES, BY SPECIFYING
c  $ THE TWO NODES.  TO FINISH, USE 0, 0 FOR THE NODES.
c  $1 2
c  $0 0
c  $ ENTER VOLUME, CROSS-SECTIONAL AREA, AND 2 BOX VALUES
c  $ FOR EACH LINK.
c  $2.0 1.0 1.0 1.0
c  $ ENTER GAS-LIQUID FRICTION COEFFICIENT.
c  $3.43
c  $ ENTER GAS-WALL FRICTION COEFFICIENT.
c  $0.06
c  $ ENTER LIQUID-WALL FRICTION COEFFICIENT.
c  $0.021
c  $ ENTER 3 COMPONENTS OF GRAVITY VECTOR:
c  $0.0 0.0 32.17
c  $ ENTER 0 FOR DONOR CELL, 1 FOR CELL AVERAGE.
c  $0
c  $ ENTER CRITICAL DENSITY, RHO-CRIT
c  $60.0
c  $ ENTER CRITICAL Sigma
c  $0.161
c  $ ENTER CRITICAL W
c  $1.24
c  $ ENTER INITIAL VALUES OF Y
c  $0.25
c  $0.25
c  $100.0
c  $100.0
c  $1190.0
c  $1190.0
c  $295.0
c  $295.0
c  $0.0
c  $0.0
c  $ ENTER INITIAL VALUES OF YDOT
c  $0.0
c  $0.0
c  $0.0
c  $0.0
c  $0.0
c  $0.0
c  $0.0
c  $0.0
c  $-16.0850
c  $-16.0850
c  $ ENTER 0 IF YOU BELIEVE YDOT IS NOT CONSISTENT WITH Y,
c  $      1 OTHERWISE.
c  $0
c  $ ENTER NAME OF RESTART FILE TO CREATE,
c  $ OR return FOR NO RESTART FILE.
c  $t2node.in
c  $ WANT TO GO ON AND SOLVE THE ode?
c  $n
c  $ NETODE IS STOPPING NOW.
c
c  Now, the problem definition is written in the file T2NODE.IN
c---------------------------------------------------------------   
c  On subsequent runs:
c
c  $run netode
c  $ IS THE PROBLEM DEFINITION IN A FILE?
c  $y  
c  $ ENTER NAME OF RESTART FILE TO BE READ:
c  $t2node.in
c  $ WANT TO GO ON AND SOLVE THE ode?
c  $y
c  $ ENTER ABSOLUTE AND RELATIVE ERROR TOLERANCES.
c  $0.005 0.0005
c  $ ENTER THE TARGET TIME
c  $5.0
c  $ ENTER NUMBER OF STEPS TO TAKE.
c  $5
c  $ ENTER LSODI METHOD mf:
c  $ 12 FOR ADAMS FORMULAS, OR 
c  $ 22 FOR BDF FORMULAS.
c  $22
c  $ ENTER itask
c  $  1 TO SEE FINAL SOLUTION ONLY,
c  $  2 TO SEE INTERMEDIATE SOLUTIONS.
c  $1
c  $ (running, output flying on the screen)
c  $ ENTER NAME OF RESTART FILE TO CREATE,
c  $ OR return FOR NO RESTART FILE.
c  $t2node.out
c
c  Now, by loading T2NODE.OUT, you can restart the NETODE at the target
c  time you specified before.
c-----------------------------------------------------------------------
c  BATCH JOB:
c
c   edit a file NETODE.CTL contain the following lines
c ______________
c   y
c   t2node.in
c   y
c   0.005 0.0005
c   5.0
c   5
c   22
c   1
c   t2node.out
c ______________
c
c   edit a file NETODE.COM contain the following lines
c --------------------------------
c   $ASSIGN NETODE.CTL SYS$INPUT
c   $ASSIGN T2NODE.DAT SYS$OUTPUT
c   $RUN NETODE
c ---------------------------------
c
c   SUBMIT NETODE.COM duplicates the interactive run with output 
c                     written in T2NODE.DAT
c
c******************************************************************
c     DRIVER
c****************************************************************
c
c  ALGLNK  The value of the gas fraction on link I.
c
c  ALLLNK  The value of the liquid fraction on link I.
c
c  ATOL    The absolute error tolerance used by LSODI.
c
c  CRSLNK  The cross-sectional area of a link.
c
c  FRICGW  The gas/wall friction coefficient.
c
c  FRICGL  The gas/liquid friction coefficient.
c
c  FRICLW  The liquid/wall friction coefficient.
c
c  GRAVIT  The gravitational vector.
c
c  HOTNOD  The heat source at a node.
c
c  ITASK   1, LSODI is to return only upon reaching TOUT.
c          2, LSODI returns on each step.
c
c  ITOL    1, ATOL and RTOL are both scalars.
c
c  LIW     The storage used for IWORK.
c
c  LNKNOD  If we implicitly order the links by their starting node,
c          and, for equal starting nodes, by link number, then
c          LNKNOD(1,*) lists the links in this order.
c
c          Similarly, LNKNOD(2,*) works for end nodes.
c
c          This is the only variable whose name ends in "NOD" which
c          is actually a link-based quantity.
c
c  LNKOPT  0, use the donor cell strategy for link quantities.
c          1, use cell averaging for link quantities.
c
c  LRW     The storage used for RWORK.
c
c  MAXIW   The maximum storage for IWORK.
c
c  MAXRW   The maximum storage for RWORK.
c
c  MF      12, Adams formulas, approx jacobian, full storage.
c          22, BDF formulas, approx jacobian, full storage.
c
c  NABNOD  NABNOD(1,I) is the number of links for which node
c          I is the beginning node.  Similarly, NABNOD(2,I)
c          is the number of links for which node I is the end node.
c
c  NODLNK  NODLNK(1,I) and NODLNK(2,I) are the numbers of the
c          beginning and end nodes for link I.
c
c  NSTEP   The number of intermediate solutions that LSODI is
c          requested to compute between the current time and the
c          target time.
c
c  NUMEQN  The number of equations.
c
c  NUMLNK  The number of links.
c
c  NUMNOD  The number of nodes.
c
c  POSNOD  The X, Y and Z coordinates of each node.
c
c  RHCRIT  The critical density.
c
c  RTOL    The relative error tolerance used by LSODI.
c
c  SCRIT   The critical value of sigma.
c
c  VOLLNK  The volume associated with a link.
c
c  VOLNOD  The volume associated with a node.
c
c  WCRIT   The critical value of W.
c
c  Y       Real Y(NUMEQN), the unknowns, as used by LSODI.
c
c          The first NUMNOD quantities are ALPHA(I), the gas
c          fraction at node I.
c
c          The next NUMNOD quantities are PRESS(I), the pressure
c          at node I.
c
c          The next NUMNOD quantities are HGAS(I), the gas
c          enthalpy at node I.
c
c          The next NUMNOD quantities are HLIQ(I), the liquid
c          enthalpy at node I.
c
c          The next NUMLNK quantities are VGAS(I), the gas
c          velocity at link I.
c
c          The next NUMLNK quantities are VLIQ(I), the liquid
c          velocity at link I.
c
      INTEGER MAXEQN
      INTEGER MAXIW
      INTEGER MAXLNK
      INTEGER MAXNOD
      INTEGER MAXRW

      PARAMETER (MAXIW=125)
      PARAMETER (MAXLNK=15)
      PARAMETER (MAXNOD=15)
      PARAMETER (MAXRW=2500)

      PARAMETER (MAXEQN=4*MAXNOD+2*MAXLNK)

      REAL ALGLNK(MAXLNK)
      REAL ALLLNK(MAXLNK)
      REAL AMGNOD(MAXNOD)
      REAL AMLNOD(MAXNOD)
      REAL ATOL(1)
      REAL BOXLNK(2,MAXLNK)
      REAL CRSLNK(MAXLNK)
      REAL FRICGL
      REAL FRICGW
      REAL FRICLW
      REAL GAMNOD(MAXNOD)
      REAL GRAVIT(3)
      REAL HGLNK(MAXLNK)
      REAL HLLNK(MAXLNK)
      REAL HOTNOD(MAXNOD)
      INTEGER I
      INTEGER IOPT
      INTEGER ISAVE(36)
      CHARACTER*1 ISAY
      INTEGER ISTATE
      INTEGER ITASK
      INTEGER ITOL
      INTEGER IWORK(MAXIW)
      INTEGER LIW
      INTEGER LNKNOD(2,MAXLNK)
      INTEGER LNKOPT
      INTEGER LRW
      INTEGER MF
      INTEGER NABNOD(2,MAXNOD)
      INTEGER NODLNK(2,MAXLNK)
      INTEGER NSTEP
      INTEGER NUMEQN
      INTEGER NUMLNK
      INTEGER NUMNOD
      REAL POSNOD(3,MAXNOD)
      REAL POTNOD(MAXNOD)
      REAL PR1LNK(MAXLNK)
      REAL PR2LNK(MAXLNK)
      REAL PRSLNK(MAXLNK)
      REAL RGHLNK(MAXLNK)
      REAL RGHNOD(MAXNOD)
      REAL RGPLNK(MAXLNK)
      REAL RGPNOD(MAXNOD)
      REAL RHCRIT
      REAL RHGLNK(MAXLNK)
      REAL RHGNOD(MAXNOD)
      REAL RHLLNK(MAXLNK)
      REAL RHLNOD(MAXNOD)
      REAL RLHLNK(MAXLNK)
      REAL RLHNOD(MAXNOD)
      REAL RLPLNK(MAXLNK)
      REAL RLPNOD(MAXNOD)
      REAL RSAVE(218)
      REAL RTOL(1)
      REAL RWORK(MAXRW)
      REAL SCRIT
      REAL TDELT
      REAL TIME
      CHARACTER*80 TITLE
      REAL TOUT
      REAL TTARG
      REAL VOLLNK(MAXLNK)
      REAL VOLNOD(MAXNOD)
      REAL WCRIT
      REAL Y(MAXEQN)
      REAL YDOT(MAXEQN)

      EXTERNAL ADDNET
      EXTERNAL JACNET
      EXTERNAL RESNET

      COMMON /PROBLM/
     &  ALGLNK,ALLLNK,AMGNOD,AMLNOD,BOXLNK,CRSLNK,FRICGL,FRICGW,FRICLW,
     &  GAMNOD,GRAVIT,HGLNK,HLLNK,HOTNOD,LNKNOD,LNKOPT,NABNOD,NODLNK,
     &  NUMEQN,NUMLNK,NUMNOD,POSNOD,POTNOD,PRSLNK,RGHLNK,RGHNOD,
     &  RGPLNK,RGPNOD,RHCRIT,RHGLNK,RHGNOD,RHLLNK,RHLNOD,RLHLNK,
     &  RLHNOD,RLPLNK,RLPNOD,SCRIT,VOLLNK,VOLNOD,WCRIT,PR1LNK,PR2LNK
c
c  Is the problem data to be read from a file?
c
      WRITE(*,*)' '
      WRITE(*,*)'  IS THE PROBLEM DEFINITION IN A FILE?'
      READ(*,'(A1)')ISAY

      IF (ISAY.EQ.'Y'.OR.ISAY.EQ.'y') THEN

        CALL READER(BOXLNK,CRSLNK,FRICGL,FRICGW,FRICLW,GRAVIT,
     &    HOTNOD,ISAVE,ISTATE,IWORK,LIW,LNKOPT,LRW,
     &    MAXEQN,MAXIW,MAXLNK,MAXNOD,MAXRW,NODLNK,NUMEQN,NUMLNK,
     &    NUMNOD,
     &    POSNOD,RHCRIT,RSAVE,RWORK,SCRIT,TIME,TITLE,VOLLNK,VOLNOD,
     &    WCRIT,Y,YDOT)

      ELSE

        CALL GONET(BOXLNK,CRSLNK,FRICGL,FRICGW,FRICLW,GRAVIT,HOTNOD,
     &    ISAVE,ISTATE,IWORK,LIW,LNKOPT,LRW,MAXEQN,MAXIW,MAXLNK,MAXNOD,
     &    MAXRW,NODLNK,NUMEQN,NUMLNK,NUMNOD,POSNOD,RHCRIT,RSAVE,RWORK,
     &    SCRIT,TIME,TITLE,VOLLNK,VOLNOD,WCRIT,Y,YDOT)

        CALL WRITER(BOXLNK,CRSLNK,FRICGL,FRICGW,FRICLW,GRAVIT,
     &    HOTNOD,ISAVE,ISTATE,IWORK,LIW,LNKOPT,
     &    LRW,MAXEQN,MAXLNK,MAXNOD,NODLNK,NUMEQN,NUMLNK,NUMNOD,
     &    POSNOD,RHCRIT,RSAVE,RWORK,SCRIT,TIME,TITLE,VOLLNK,VOLNOD,
     &    WCRIT,Y,YDOT)

      end if
c
c  Construct the inverse node-link maps.
c
      CALL NLMAP(LNKNOD,NABNOD,NODLNK,NUMLNK,NUMNOD)
c
c  Edit the starting point
c
      CALL EDINET(NUMEQN,TIME,TITLE,Y,YDOT)
c
c  Stop or go on?
c
      WRITE(*,*)' '
      WRITE(*,*)'WANT TO GO ON AND SOLVE THE ode?'
      READ(*,'(A1)')ISAY

      IF(ISAY.NE.'Y'.AND.ISAY.NE.'y')THEN
        WRITE(*,*)' '
        WRITE(*,*)'NETODE IS STOPPING NOW.'
        STOP
      end if
c
c  Error tolerances
c
      ITOL=1
      WRITE(*,*)'ENTER ABSOLUTE AND RELATIVE ERROR TOLERANCES.'
      READ(*,*)ATOL(1),RTOL(1)
c
c  Get target time.
c
      WRITE(*,*)' '
      WRITE(*,*)'ENTER THE TARGET TIME'
      READ(*,*)TTARG

      WRITE(*,*)' '
      WRITE(*,*)'ENTER NUMBER OF STEPS TO TAKE.'
      READ(*,*)NSTEP

      WRITE(*,*)' '
      WRITE(*,*)'ENTER LSODI METHOD mf:'
      WRITE(*,*)'  12 FOR ADAMS FORMULAS, OR '
      WRITE(*,*)'  22 FOR BDF FORMULAS.'
      READ(*,*)MF

      IF(MF.NE.12)THEN
        MF=22
      end if

      WRITE(*,*)'ENTER itask'
      WRITE(*,*)'  1 TO SEE FINAL SOLUTION ONLY,'
      WRITE(*,*)'  2 TO SEE INTERMEDIATE SOLUTIONS.'
      READ(*,*)ITASK

      IF(ITASK.NE.2)THEN
        ITASK=1
      end if
c
c  No optional inputs to LSODI.
c
      IOPT=0
c
c  Call LSODI.
c
      TDELT = (TTARG - TIME)/REAL(NSTEP)

      DO I=1,NSTEP

        TOUT = TIME + TDELT

        CALL LSODI(RESNET,ADDNET,JACNET,NUMEQN,Y,YDOT,TIME,TOUT,ITOL,
     &    RTOL,ATOL,ITASK,ISTATE,IOPT,RWORK,LRW,IWORK,LIW,MF)

        IF(ISTATE.NE.2)THEN
          WRITE(*,*)' '
          WRITE(*,*)'NETODE - FATAL ERROR!'
          WRITE(*,*)'  LSODI INTEGRATION HALTED.'
          WRITE(*,*)'  VALUE OF istate=',ISTATE
          STOP
        end if
c
c  Print out information about current solution.
c
        CALL EDINET(NUMEQN,TIME,TITLE,Y,YDOT)

      end do
c
c  Save problem information on disk for possible restart?
c
      CALL WRITER(BOXLNK,CRSLNK,FRICGL,FRICGW,FRICLW,GRAVIT,
     &  HOTNOD,ISAVE,ISTATE,IWORK,LIW,LNKOPT,
     &  LRW,MAXEQN,MAXLNK,MAXNOD,NODLNK,NUMEQN,NUMLNK,NUMNOD,
     &  POSNOD,RHCRIT,RSAVE,RWORK,SCRIT,TIME,TITLE,VOLLNK,VOLNOD,
     &  WCRIT,Y,YDOT)

      WRITE(*,*)' '
      WRITE(*,*)'LSODI STATISTICS:'
      WRITE(*,*)' '
      WRITE(*,*)'LAST LSODI STEPSIZE=',RWORK(11)
      WRITE(*,*)'NEXT LSODI STEPSIZE=',RWORK(12)
      WRITE(*,*)'ACTUAL VALUE OF t = ',RWORK(13)
      WRITE(*,*)'NUMBER OF lsodi STEPS =',IWORK(11)
      WRITE(*,*)'NUMBER OF RESIDUALS =  ',IWORK(12)
      WRITE(*,*)'NUMBER OF JACOBIANS =  ',IWORK(13)
      WRITE(*,*)'LAST ORDER USED =      ',IWORK(14)
      WRITE(*,*)'NEXT ORDER TO USE =    ',IWORK(15)
      WRITE(*,*)'NEEDED REAL WORKSPACE =',IWORK(17)
      WRITE(*,*)'NEEDED INTEGER SPACE = ',IWORK(18)

      STOP
      END
      SUBROUTINE ADDNET(N,TIME,Y,ML,MU,PA,M0)

c*********************************************************************72
c
cc ADDNET adds the matrix A(T,Y) to the input matrix PA.
c
c  Discussion:
c
c    Here, we assume that the differential system being solved has
c    the form
c
c      A(T,Y)*YDOT = G(T,Y)
c
      INTEGER MAXLNK
      INTEGER MAXNOD
      REAL PCON1

      PARAMETER (MAXLNK=15)
      PARAMETER (MAXNOD=15)
      PARAMETER (PCON1=0.18509)

      INTEGER M0
      INTEGER N

      REAL ALGAS
      REAL ALGAS1
      REAL ALGAS2
      REAL ALGLNK(MAXLNK)
      REAL ALLIQ
      REAL ALLLNK(MAXLNK)
      REAL AMGNOD(MAXNOD)
      REAL AMLNOD(MAXNOD)
      REAL BOX1
      REAL BOX2
      REAL BOXLNK(2,MAXLNK)
      REAL CRSLNK(MAXLNK)
      REAL DRGDH
      REAL DRGDH1
      REAL DRGDH2
      REAL DRGDP
      REAL DRGDP1
      REAL DRGDP2
      REAL DRLDH
      REAL DRLDH1
      REAL DRLDH2
      REAL DRLDP
      REAL DRLDP1
      REAL DRLDP2
      REAL FRICGL
      REAL FRICGW
      REAL FRICLW
      REAL GAMNOD(MAXNOD)
      REAL GRAVIT(3)
      REAL HGAS
      REAL HGAS1
      REAL HGAS2
      REAL HGASS
      REAL HGLNK(MAXLNK)
      REAL HLIQ
      REAL HLIQ1
      REAL HLIQ2
      REAL HLIQS
      REAL HLLNK(MAXLNK)
      REAL HOTNOD(MAXNOD)
      INTEGER I
      INTEGER IALGAS
      INTEGER IEQN
      INTEGER IGDONR
      INTEGER IHGAS
      INTEGER IHLIQ
      INTEGER ILDONR
      INTEGER IPRESS
      INTEGER IVAR
      INTEGER IVGAS
      INTEGER IVLIQ
      INTEGER J
      INTEGER JALGAS
      INTEGER JHGAS
      INTEGER JHLIQ
      INTEGER JPRESS
      INTEGER KALGAS
      INTEGER KHGAS
      INTEGER KHLIQ
      INTEGER KPRESS
      INTEGER LNKNOD(2,MAXLNK)
      INTEGER LNKOPT
      INTEGER ML
      INTEGER MU
      INTEGER NABNOD(2,MAXNOD)
      INTEGER NODE1
      INTEGER NODE2
      INTEGER NODLNK(2,MAXLNK)
      INTEGER NUMEQN
      INTEGER NUMLNK
      INTEGER NUMNOD
      REAL PA(M0,N)
      REAL POSNOD(3,MAXNOD)
      REAL POTNOD(MAXNOD)
      REAL PR1LNK(MAXLNK)
      REAL PR2LNK(MAXLNK)
      REAL PRESS
      REAL PRESS1
      REAL PRESS2
      REAL PRSLNK(MAXLNK)
      REAL RGHLNK(MAXLNK)
      REAL RGHNOD(MAXNOD)
      REAL RGPLNK(MAXLNK)
      REAL RGPNOD(MAXNOD)
      REAL RHCRIT
      REAL RHGAS
      REAL RHGAS1
      REAL RHGAS2
      REAL RHGLNK(MAXLNK)
      REAL RHGNOD(MAXNOD)
      REAL RHLIQ
      REAL RHLIQ1
      REAL RHLIQ2
      REAL RHLLNK(MAXLNK)
      REAL RHLNOD(MAXNOD)
      REAL RLHLNK(MAXLNK)
      REAL RLHNOD(MAXNOD)
      REAL RLPLNK(MAXLNK)
      REAL RLPNOD(MAXNOD)
      REAL SCRIT
      REAL TIME
      REAL VGAS
      REAL VLIQ
      REAL VOLLNK(MAXLNK)
      REAL VOLNOD(MAXNOD)
      REAL WCRIT
      REAL Y(N)

      COMMON /PROBLM/
     &  ALGLNK,ALLLNK,AMGNOD,AMLNOD,BOXLNK,CRSLNK,FRICGL,FRICGW,FRICLW,
     &  GAMNOD,GRAVIT,HGLNK,HLLNK,HOTNOD,LNKNOD,LNKOPT,NABNOD,NODLNK,
     &  NUMEQN,NUMLNK,NUMNOD,POSNOD,POTNOD,PRSLNK,RGHLNK,RGHNOD,
     &  RGPLNK,RGPNOD,RHCRIT,RHGLNK,RHGNOD,RHLLNK,RHLNOD,RLHLNK,
     &  RLHNOD,RLPLNK,RLPNOD,SCRIT,VOLLNK,VOLNOD,WCRIT,PR1LNK,PR2LNK

      DO I=1,NUMNOD

        IALGAS=I
        IPRESS=IALGAS+NUMNOD
        IHGAS=IPRESS+NUMNOD
        IHLIQ=IHGAS+NUMNOD

        ALGAS=Y(IALGAS)
        ALLIQ=1.0-ALGAS
        PRESS=Y(IPRESS)
        HGAS=Y(IHGAS)
        HLIQ=Y(IHLIQ)
c
c  Evaluate RHGAS, DRGDP, DRGDH, HGASS,
c           RHLIQ, DRLDP, DRLDH, HLIQS.
c
        CALL DENSE(DRGDH,DRGDP,DRLDH,DRLDP,HGAS,HGASS,HLIQ,HLIQS,PRESS,
     &    RHGAS,RHLIQ)

        RHGNOD(I)=RHGAS
        RGPNOD(I)=DRGDP
        RGHNOD(I)=DRGDH

        RHLNOD(I)=RHLIQ
        RLPNOD(I)=DRLDP
        RLHNOD(I)=DRLDH

        IF(HOTNOD(I).GT.0.0.AND.HLIQ.GE.HLIQS ) THEN
          GAMNOD(I)=ALLIQ*HOTNOD(I)/(HGASS-HLIQS)
        ELSE IF (HOTNOD(I).LT.0.0.AND.HLIQ.LE.HGASS) THEN
          GAMNOD(I)=ALGAS*HOTNOD(I)/(HGASS-HLIQS)
        ELSE
          GAMNOD(I)=0.0
        end if

        POTNOD(I)=0.0

        DO J=1,3
          POTNOD(I)=POTNOD(I)+GRAVIT(J)*POSNOD(J,I)
        end do

      end do
c
c  Get pressure and other variables at every link.
c
      DO I=1,NUMLNK
        IVGAS=4*NUMNOD+I
        IVLIQ=IVGAS+NUMLNK
        VGAS=Y(IVGAS)
        VLIQ=Y(IVLIQ)
c
c  Get the two end nodes of link I.
c
        NODE1=NODLNK(1,I)
        NODE2=NODLNK(2,I)
        BOX1=BOXLNK(1,I)
        BOX2=BOXLNK(2,I)
c
c  Get the indices and values of node based quantities.
c
        JALGAS=NODE1
        JPRESS=JALGAS+NUMNOD
        JHGAS=JPRESS+NUMNOD
        JHLIQ=JHGAS+NUMNOD
        ALGAS1=Y(JALGAS)
        PRESS1=Y(JPRESS)
        HGAS1=Y(JHGAS)
        HLIQ1=Y(JHLIQ)
        RHGAS1=RHGNOD(NODE1)
        RHLIQ1=RHLNOD(NODE1)
        DRGDH1=RGHNOD(NODE1)
        DRGDP1=RGPNOD(NODE1)
        DRLDH1=RLHNOD(NODE1)
        DRLDP1=RLPNOD(NODE1)
c
c  Same for node 2.
c
        KALGAS=NODE2
        KPRESS=KALGAS+NUMNOD
        KHGAS=KPRESS+NUMNOD
        KHLIQ=KHGAS+NUMNOD
        ALGAS2=Y(KALGAS)
        PRESS2=Y(KPRESS)
        HGAS2=Y(KHGAS)
        HLIQ2=Y(KHLIQ)
        RHGAS2=RHGNOD(NODE2)
        RHLIQ2=RHLNOD(NODE2)
        DRGDH2=RGHNOD(NODE2)
        DRGDP2=RGPNOD(NODE2)
        DRLDH2=RLHNOD(NODE2)
        DRLDP2=RLPNOD(NODE2)
c
c  Evaluate link quantities which depend on nodal values.
c
        CALL DOLINK(ALGAS,ALGAS1,ALGAS2,ALLIQ,BOX1,BOX2,DRGDH,
     &    DRGDH1,DRGDH2,DRGDP,DRGDP1,DRGDP2,DRLDH,DRLDH1,DRLDH2,DRLDP,
     &    DRLDP1,DRLDP2,HGAS,HGAS1,HGAS2,HLIQ,HLIQ1,HLIQ2,LNKOPT,PRESS,
     &    PRESS1,PRESS2,RHGAS,RHGAS1,RHGAS2,RHLIQ,RHLIQ1,RHLIQ2,VGAS,
     &    VLIQ)

        ALGLNK(I)=ALGAS
        ALLLNK(I)=ALLIQ
        HGLNK(I)=HGAS
        HLLNK(I)=HLIQ
        PR1LNK(I)=PRESS1
        PR2LNK(I)=PRESS2
        PRSLNK(I)=PRESS
        RHGLNK(I)=RHGAS
        RHLLNK(I)=RHLIQ
        RGHLNK(I)=DRGDH
        RGPLNK(I)=DRGDP
        RLHLNK(I)=DRLDH
        RLPLNK(I)=DRLDP
      end do
c
c  Add the jacobian terms to the input matrix.
c
      IEQN=0

      DO I=1,NUMNOD
c
c  Read off the nodal values.
c
        IALGAS=I
        IPRESS=IALGAS+NUMNOD
        IHGAS=IPRESS+NUMNOD
        IHLIQ=IHGAS+NUMNOD

        ALGAS=Y(IALGAS)
        ALLIQ=1.0-ALGAS
        PRESS=Y(IPRESS)
        HGAS=Y(IHGAS)
        HLIQ=Y(IHLIQ)
        RHGAS=RHGNOD(I)
        RHLIQ=RHLNOD(I)
        DRGDP=RGPNOD(I)
        DRGDH=RGHNOD(I)
        DRLDP=RLPNOD(I)
        DRLDH=RHLNOD(I)
c
c  1.  Gas continuity equation at node I.
c
        IEQN=IEQN+1

        IVAR=IALGAS
        PA(IEQN,IVAR)=PA(IEQN,IVAR)+RHGAS

        IVAR=IPRESS
        PA(IEQN,IVAR)=PA(IEQN,IVAR)+ALGAS*DRGDP

        IVAR=IHGAS
        PA(IEQN,IVAR)=PA(IEQN,IVAR)+ALGAS*DRGDH
c
c  2.  Liquid continuity equation at node I.
c
        IEQN=IEQN+1

        IVAR=IALGAS
        PA(IEQN,IVAR)=PA(IEQN,IVAR)-RHLIQ

        IVAR=IPRESS
        PA(IEQN,IVAR)=PA(IEQN,IVAR)+ALLIQ*DRLDP

        IVAR=IHLIQ
        PA(IEQN,IVAR)=PA(IEQN,IVAR)+ALLIQ*DRLDH
c
c  3.  Gas energy equation at node I.
c
        IEQN=IEQN+1

        IVAR=IALGAS
        PA(IEQN,IVAR)=PA(IEQN,IVAR)-PCON1*PRESS

        IVAR=IPRESS
        PA(IEQN,IVAR)=PA(IEQN,IVAR)-PCON1*ALGAS

        IVAR=IHGAS
        PA(IEQN,IVAR)=PA(IEQN,IVAR)+ALGAS*RHGAS
c
c  4.  Liquid energy equation at node I.
c
        IEQN=IEQN+1

        IVAR=IALGAS
        PA(IEQN,IVAR)=PA(IEQN,IVAR)+PCON1*PRESS

        IVAR=IPRESS
        PA(IEQN,IVAR)=PA(IEQN,IVAR)-PCON1*ALLIQ

        IVAR=IHLIQ
        PA(IEQN,IVAR)=PA(IEQN,IVAR)+ALLIQ*RHLIQ

      end do
c
c  Now do equations at links.
c
      DO I=1,NUMLNK

        IVGAS=4*NUMNOD+I
        IVLIQ=IVGAS+NUMLNK
        VGAS=Y(IVGAS)
        VLIQ=Y(IVLIQ)
        NODE1=NODLNK(1,I)
        NODE2=NODLNK(2,I)
        BOX1=BOXLNK(1,I)
        BOX2=BOXLNK(2,I)
        IF(VGAS.GE.0.0)THEN
          IGDONR=NODE1
        ELSE
          IGDONR=NODE2
        end if
        IF(VLIQ.GE.0.0)THEN
          ILDONR=NODE1
        ELSE
          ILDONR=NODE2
        end if
        JALGAS=NODE1
        JPRESS=JALGAS+NUMNOD
        JHGAS=JPRESS+NUMNOD
        JHLIQ=JHGAS+NUMNOD
        KALGAS=NODE2
        KPRESS=KALGAS+NUMNOD
        KHGAS=KPRESS+NUMNOD
        KHLIQ=KHGAS+NUMNOD
c
c  5.  Gas momentum equation at link I.
c
        IEQN=IEQN+1

        IVAR=IVGAS
        PA(IEQN,IVAR)=PA(IEQN,IVAR)+ALGLNK(I)*RHGLNK(I)

        IF(LNKOPT.EQ.1)THEN

          IVAR=JALGAS
          PA(IEQN,IVAR)=PA(IEQN,IVAR)+0.5*VGAS*RHGLNK(I)*BOX1
     &      /(BOX1+BOX2)

          IVAR=KALGAS
          PA(IEQN,IVAR)=PA(IEQN,IVAR)+0.5*VGAS*RHGLNK(I)*BOX2
     &      /(BOX1+BOX2)

          IVAR=JPRESS
          PA(IEQN,IVAR)=PA(IEQN,IVAR)+0.25*VGAS*ALGLNK(I)*RGPLNK(I)

          IVAR=KPRESS
          PA(IEQN,IVAR)=PA(IEQN,IVAR)+0.25*VGAS*ALGLNK(I)*RGPLNK(I)

          IF(IGDONR.EQ.NODE1)THEN
            IVAR=JHGAS
          ELSE
            IVAR=KHGAS
          end if
          PA(IEQN,IVAR)=PA(IEQN,IVAR)+0.5*VGAS*ALGLNK(I)*RGHLNK(I)

        end if
c
c  6.  Liquid momentum equation at link I.
c
        IEQN=IEQN+1

        IVAR=IVLIQ

        PA(IEQN,IVAR)=PA(IEQN,IVAR)+ALLLNK(I)*RHLLNK(I)

        IF(LNKOPT.EQ.1)THEN

          IVAR=JALGAS
          PA(IEQN,IVAR)=PA(IEQN,IVAR)-0.5*VLIQ*RHLLNK(I)*BOX1
     &      /(BOX1+BOX2)

          IVAR=KALGAS
          PA(IEQN,IVAR)=PA(IEQN,IVAR)-0.5*VLIQ*RHLLNK(I)*BOX2
     &      /(BOX1+BOX2)

          IVAR=JPRESS
          PA(IEQN,IVAR)=PA(IEQN,IVAR)+0.25*VLIQ*ALLLNK(I)*RLPLNK(I)

          IVAR=KPRESS
          PA(IEQN,IVAR)=PA(IEQN,IVAR)+0.25*VLIQ*ALLLNK(I)*RLPLNK(I)

          IF(ILDONR.EQ.NODE1)THEN
            IVAR=JHLIQ
          ELSE
            IVAR=KHLIQ
          end if
          PA(IEQN,IVAR)=PA(IEQN,IVAR)+0.5*VLIQ*ALLLNK(I)*RLHLNK(I)

        end if

      end do

      RETURN
      END
      SUBROUTINE COOKER(ALGASD,ALLIQD,DRGDH,DRGDP,DRLDH,DRLDP,HGAS,
     &  HLIQ,IERROR,PRESS,PRESS1,PRESS2,RHGAS,RHGASD,RHLIQ,RHLIQD)

c*********************************************************************72
c
cc COOKER solves a nonlinear equation for pressure.
c
c  Discussion:
c
c    The routine computes the Newton iterates of pressure for fixed values
c    of HGAS and HLIQ in an attempt to solve the nonlinear equation:
c
c      ALGASD * RHGASD * RHOL(PRESS,HLIQ) 
c    + ALLIQD * RHLIQD * RHOG(PRESS,HGAS)
c    = RHOL(PRESS,HLIQ) * RHOG(PRESS,HGAS)
c
      REAL ALGASD
      REAL ALLIQD
      REAL DRGDH
      REAL DRGDP
      REAL DRLDH
      REAL DRLDP
      REAL FN
      REAL FP
      REAL HGAS
      REAL HGASS
      REAL HLIQ
      REAL HLIQS
      INTEGER IERROR
      INTEGER ISTEP
      REAL PRESS
      REAL PRESS1
      REAL PRESS2
      REAL RHGAS
      REAL RHGASD
      REAL RHLIQ
      REAL RHLIQD

      IERROR=0
      ISTEP=0
c
c  Initial estimate for pressure.
c
      PRESS=0.5*(PRESS1+PRESS2)

10    CONTINUE
c
c  Get the densities RHGAS and RHLIQ for the current pressure PRESS and
c  fixed enthalpies HGAS and HLIQ, and the partial derivatives
c  DRLDP and DRGDP.
c
      CALL DENSE(DRGDH,DRGDP,DRLDH,DRLDP,HGAS,HGASS,HLIQ,HLIQS,PRESS,
     &  RHGAS,RHLIQ)
c
c  Evaluate the function.
c
      FN=ALGASD*RHGASD*RHLIQ+ALLIQD*RHLIQD*RHGAS-RHGAS*RHLIQ

      IF(ABS(FN).LE.0.005.AND.ISTEP.GT.3)THEN
        RETURN
      end if
c
c  Evaluate the derivative of the function.
c
      FP=(ALGASD*RHGASD-RHGAS)*DRLDP+(ALLIQD*RHLIQD-RHLIQ)*DRGDP

      PRESS=PRESS-FN/FP

      IF(ABS(FN/FP).LE.(0.001+0.0001*ABS(PRESS)).AND.ISTEP.GT.3)THEN
        CALL DENSE(DRGDH,DRGDP,DRLDH,DRLDP,HGAS,HGASS,HLIQ,HLIQS,PRESS,
     &    RHGAS,RHLIQ)
        RETURN
      end if

      IF(ISTEP.LE.5)GO TO 10
c
c  No convergence
c
      IERROR=1
      CALL DENSE(DRGDH,DRGDP,DRLDH,DRLDP,HGAS,HGASS,HLIQ,HLIQS,PRESS,
     &  RHGAS,RHLIQ)

      WRITE(*,*)' '
      WRITE(*,*)'COOKER - NEWTON NOT CONVERGED AFTER ',ISTEP,' STEPS.'
      WRITE(*,*)'  NEWTON RESIDUAL fn=',FN
      WRITE(*,*)'  NEWTON DERIVATIVE fp=',FP
      WRITE(*,*)'  PRESSURE SET TO ',PRESS

      RETURN
      END
      SUBROUTINE DENSE(DRGDH,DRGDP,DRLDH,DRLDP,HGAS,HGASS,HLIQ,HLIQS,
     &  PRESS,RHGAS,RHLIQ)

c*********************************************************************72
c
cc DENSE computes density and saturation enthalpy.
c
c  Discussion:
c
c    The routine is given a press PRESS and an enthalpy H, and computes the
c    density RH, and derivatives DRDH, DRDP, and the saturation
c    enthalpy HS.
c
c    The functional form used for density is:
c
c    Let
c
c      H=MAX(HGASS,HGAS),
c      P=PRESS
c
c    Then
c
c      RHGAS = P / ( (D2+D5*H)*P**2 + (D1+D4*H)*P + (D3+H*D6))
c
c      RHLIQ = EXP(-Z), for HLIQ <= HLIQS (Z depends on PRESS & HLIQS)
c      RHLIQ = RHLIQ(HLIQS,PRESS), for HLIQ > HLIQS.
c
      REAL ATAB(9)
      REAL BOT
      REAL BTAB(9)
      REAL CTAB(3,5)
      REAL DRGDH
      REAL DRGDP
      REAL DRLDH
      REAL DRLDP
      REAL DTAB(6)
      REAL DZDH
      REAL DZDP
      REAL ETAB(5)
      REAL ETABP(5)
      REAL HCOPY
      REAL HGAS
      REAL HGASS
      REAL HLIQ
      REAL HLIQS
      INTEGER J
      REAL PCOPY
      REAL PRESS
      REAL RHGAS
      REAL RHLIQ
      REAL Z

      DATA ATAB /
     &   0.182609e+03, 0.144140e+01, -0.387216e-02, 0.651417e-05,
     &  -0.638144e-08, 0.369701e-11, -0.124626e-14, 0.225589e-18,
     &  -0.169253e-22/

      DATA BTAB /
     &   0.115216e+04, 0.460395e+00, -0.159024e-02, 0.286502e-05,
     &  -0.299850e-08, 0.185137e-11, -0.664224e-15, 0.127776e-18,
     &  -0.101790e-22 /

      DATA CTAB /
     &  -0.41345e+01, -0.59428e-05,  0.15681e-08,
     &   0.13252e-04,  0.63377e-07, -0.40711e-10,
     &   0.15812e-05, -0.39974e-09,  0.25401e-12,
     &  -0.21959e-08,  0.69391e-12, -0.52372e-15,
     &   0.21683e-11, -0.36159e-15,  0.32503e-18 /

      DATA DTAB /
     &  -0.81735849e-03, 0.12378514e-04, -0.10339904e+04,
     &  -0.62941689e-05, -0.87292160e-08, 0.12460225e+01 /
c
c  For the given pressure PRESS, compute the saturated gas and
c  liquid enthalpies HGASS and HLIQS.
c
      IF(PRESS.GT.3200.0)THEN
        PCOPY=3200.0
      ELSEIF(PRESS.LT.0.0)THEN
        PCOPY=0.0
      ELSE
        PCOPY=PRESS
      end if

      CALL POLY(PCOPY,ATAB,9,HLIQS)
      CALL POLY(PCOPY,BTAB,9,HGASS)

      DO J=1,5
        CALL POLY(PCOPY,CTAB(1,J),3,ETAB(J))
        CALL DPOLY(PCOPY,CTAB(1,J),3,ETABP(J))
      end do
c
c  Compute RHGAS, and its H and P partial derivatives.
c
      IF(HGAS.GT.1800.0)THEN
        HCOPY=1800.0
      ELSEIF(HGAS.LT.10.0)THEN
        HCOPY=10.0
      ELSE
        HCOPY=HGAS
      end if

      HCOPY=MAX(HCOPY,HGASS)

      BOT=((DTAB(2)+DTAB(5)*HCOPY)*PRESS+
     &  DTAB(1)+DTAB(4)*HCOPY)*PRESS+DTAB(3)+HCOPY*DTAB(6)

      RHGAS=PCOPY/BOT
      DRGDP=(DTAB(3)+HCOPY*DTAB(6)-PCOPY*(DTAB(2)+HCOPY*DTAB(5)))/BOT**2
      DRGDH=-((DTAB(5)*PCOPY+DTAB(4))*PCOPY+DTAB(6))*PCOPY/BOT**2
c
c  Compute RHLIQ, and its H and P partial derivatives.
c
      IF(HLIQ.GT.1800.0)THEN
        HCOPY=1800.0
      ELSEIF(HLIQ.LT.10.0)THEN
        HCOPY=10.0
      ELSE
        HCOPY=HLIQ
      end if

      HCOPY=MIN(HCOPY,HLIQS)

      CALL POLY(HCOPY,ETAB,5,Z)
      CALL DPOLY(HCOPY,ETAB,5,DZDH)
      CALL POLY(HCOPY,ETABP,5,DZDP)

      RHLIQ=EXP(-Z)
      DRLDH=-RHLIQ*DZDH
      DRLDP=-RHLIQ*DZDP

      RETURN
      END
      SUBROUTINE DOLINK(ALGAS,ALGAS1,ALGAS2,ALLIQ,BOX1,BOX2,DRGDH,
     &  DRGDH1,DRGDH2,DRGDP,DRGDP1,DRGDP2,DRLDH,DRLDH1,DRLDH2,DRLDP,
     &  DRLDP1,DRLDP2,HGAS,HGAS1,HGAS2,HLIQ,HLIQ1,HLIQ2,LNKOPT,PRESS,
     &  PRESS1,PRESS2,RHGAS,RHGAS1,RHGAS2,RHLIQ,RHLIQ1,RHLIQ2,VGAS,
     &  VLIQ)

c*********************************************************************72
c
cc DOLINK defines link quantities.
c
c  Discussion:
c
c    The link quantities are computed as weighted values 
c    of the quantities at the two nodes joined by the link.
c
      REAL ALGAS
      REAL ALGAS1
      REAL ALGAS2
      REAL ALLIQ
      REAL BOX1
      REAL BOX2
      REAL DRGDH
      REAL DRGDH1
      REAL DRGDH2
      REAL DRGDP
      REAL DRGDP1
      REAL DRGDP2
      REAL DRLDH
      REAL DRLDH1
      REAL DRLDH2
      REAL DRLDP
      REAL DRLDP1
      REAL DRLDP2
      REAL HGAS
      REAL HGAS1
      REAL HGAS2
      REAL HGASS
      REAL HLIQ
      REAL HLIQ1
      REAL HLIQ2
      REAL HLIQS
      INTEGER LNKOPT
      REAL PRESS
      REAL PRESS1
      REAL PRESS2
      REAL RHGAS
      REAL RHGAS1
      REAL RHGAS2
      REAL RHLIQ
      REAL RHLIQ1
      REAL RHLIQ2
      REAL VGAS
      REAL VLIQ
c
c  Node 1 is gas and liquid donor.
c
      IF(VGAS.GE.0.0.AND.VLIQ.GE.0.0)THEN

        HGAS=HGAS1
        HLIQ=HLIQ1

        IF(LNKOPT.EQ.0)THEN
          ALGAS=ALGAS1
          ALLIQ=1.0-ALGAS
          DRGDH=DRGDH1
          DRGDP=DRGDP1
          DRLDH=DRLDH1
          DRLDP=DRLDP1
          PRESS=PRESS1
          RHGAS=RHGAS1
          RHLIQ=RHLIQ1
        ELSE
          ALGAS=(ALGAS1*BOX1+ALGAS2*BOX2)/(BOX1+BOX2)
          ALLIQ=1.0-ALGAS
          PRESS=(PRESS1+PRESS2)/2.0
          CALL DENSE(DRGDH,DRGDP,DRLDH,DRLDP,HGAS,HGASS,HLIQ,HLIQS,
     &      PRESS,RHGAS,RHLIQ)
        end if
c
c  Node 2 is gas and liquid donor.
c
      ELSEIF(VGAS.LT.0.0.AND.VLIQ.LT.0.0)THEN

        HGAS=HGAS2
        HLIQ=HLIQ2

        IF(LNKOPT.EQ.0)THEN
          ALGAS=ALGAS2
          ALLIQ=1.0-ALGAS
          DRGDH=DRGDH2
          DRGDP=DRGDP2
          DRLDH=DRLDH2
          DRLDP=DRLDP2
          PRESS=PRESS2
          RHGAS=RHGAS2
          RHLIQ=RHLIQ2
        ELSE
          ALGAS=(ALGAS1*BOX1+ALGAS2*BOX2)/(BOX1+BOX2)
          ALLIQ=1.0-ALGAS
          PRESS=(PRESS1+PRESS2)/2.0
          CALL DENSE(DRGDH,DRGDP,DRLDH,DRLDP,HGAS,HGASS,HLIQ,HLIQS,
     &      PRESS,RHGAS,RHLIQ)
        end if
c
c  Node 1 is gas donor, node 2 is liquid donor.
c
      ELSEIF(VGAS.GE.0.0.AND.VLIQ.LT.0.0)THEN

        HGAS=HGAS1
        HLIQ=HLIQ2
        PRESS=(PRESS1+PRESS2)/2.0

        IF(LNKOPT.EQ.0)THEN
          ALGAS=ALGAS1
          ALLIQ=1.0-ALGAS2
          DRGDH=DRGDH1
          DRGDP=DRGDP1
          DRLDH=DRLDH2
          DRLDP=DRLDP2
          RHGAS=RHGAS1
          RHLIQ=RHLIQ2
        ELSE
          ALGAS=(ALGAS1*BOX1+ALGAS2*BOX2)/(BOX1+BOX2)
          ALLIQ=1.0-ALGAS
          CALL DENSE(DRGDH,DRGDP,DRLDH,DRLDP,HGAS,HGASS,HLIQ,HLIQS,
     &      PRESS,RHGAS,RHLIQ)
        end if
c
c  Node 1 is liquid donor, node 2 is gas donor.
c
      ELSEIF(VGAS.LT.0.0.AND.VLIQ.GE.0.0)THEN

        HGAS=HGAS2
        HLIQ=HLIQ1
        PRESS=(PRESS1+PRESS2)/2.0

        IF(LNKOPT.EQ.0)THEN
          ALGAS=ALGAS2
          ALLIQ=1.0-ALGAS1
          DRGDH=DRGDH2
          DRGDP=DRGDP2
          DRLDH=DRLDH1
          DRLDP=DRLDP1
          RHGAS=RHGAS2
          RHLIQ=RHLIQ1
        ELSE
          ALGAS=(ALGAS1*BOX1+ALGAS2*BOX2)/(BOX1+BOX2)
          ALLIQ=1.0-ALGAS
          CALL DENSE(DRGDH,DRGDP,DRLDH,DRLDP,HGAS,HGASS,HLIQ,HLIQS,
     &      PRESS,RHGAS,RHLIQ)
        end if

      end if

      RETURN
      END
      SUBROUTINE DPOLY(BASE,COEF,NCOEF,PRIME)

c*********************************************************************72
c
cc DPOLY computes the derivative of a polynomial.
c
c  Discusssion:
c
c    The polynomial has the form
c
c    POLY(X) = Sum (I=1 to NCOEF) COEF(I)*BASE**(I-1)
c
      INTEGER NCOEF

      REAL BASE
      REAL COEF(NCOEF)
      INTEGER I
      REAL PRIME

      PRIME=0.0

      IF(NCOEF.LE.1)RETURN

      PRIME=(NCOEF-1)*COEF(NCOEF)

      DO I=1,NCOEF-2
        PRIME=(NCOEF-1-I)*COEF(NCOEF-I)+PRIME*BASE
      end do

      RETURN
      END
      SUBROUTINE EDINET(NUMEQN,TIME,TITLE,Y,YDOT)

c*********************************************************************72
c
cc EDINET prints out information about the current solution.
c
      INTEGER NUMEQN

      INTEGER I
      REAL TIME
      CHARACTER*80 TITLE
      REAL Y(NUMEQN)
      REAL YDOT(NUMEQN)

      WRITE(*,*)' '
      WRITE(*,'(1X,A)')TITLE
      WRITE(*,*)' '
      WRITE(*,*)'CURRENT TIME IS ',TIME
      WRITE(*,*)' '
      WRITE(*,*)'SOLUTION y AND TIME DERIVATIVE ydot:'
      WRITE(*,*)' '
      DO I=1,NUMEQN
        WRITE(*,*)I,Y(I),YDOT(I)
      end do

      RETURN
      END
      SUBROUTINE GONET(BOXLNK,CRSLNK,FRICGL,FRICGW,FRICLW,GRAVIT,HOTNOD,
     &  ISAVE,ISTATE,IWORK,LIW,LNKOPT,LRW,MAXEQN,MAXIW,MAXLNK,MAXNOD,
     &  MAXRW,NODLNK,NUMEQN,NUMLNK,NUMNOD,POSNOD,RHCRIT,RSAVE,RWORK,
     &  SCRIT,TIME,TITLE,VOLLNK,VOLNOD,WCRIT,Y,YDOT)

c*********************************************************************72
c
cc GONET gets the information that initializes the problem.
c
c  Discussion:
c
c    This information is gotten either interactively from the user,
c    or from a restart file.
c
      INTEGER MAXIW
      INTEGER MAXEQN
      INTEGER MAXLNK
      INTEGER MAXNOD
      INTEGER MAXRW

      REAL BOXLNK(2,MAXLNK)
      REAL CRSLNK(MAXLNK)
      REAL FRICGL
      REAL FRICGW
      REAL FRICLW
      REAL GRAVIT(3)
      REAL HOTNOD(MAXNOD)
      INTEGER I
      INTEGER INODE
      INTEGER ISAVE(36)
      INTEGER ISTATE
      INTEGER IWORK(MAXIW)
      INTEGER JNODE
      INTEGER LIW
      INTEGER LNKOPT
      INTEGER LRW
      INTEGER NODLNK(2,MAXLNK)
      INTEGER NUMEQN
      INTEGER NUMLNK
      INTEGER NUMNOD
      REAL POSNOD(3,MAXNOD)
      REAL RHCRIT
      REAL RSAVE(218)
      REAL RWORK(MAXRW)
      REAL SCRIT
      REAL TIME
      CHARACTER*80 TITLE
      REAL VOLLNK(MAXLNK)
      REAL VOLNOD(MAXNOD)
      REAL WCRIT
      REAL Y(MAXEQN)
      REAL YDOT(MAXEQN)

      WRITE(*,*)' '
      WRITE(*,*)'GONET:'
      WRITE(*,*)'  GET NETODE PROBLEM DEFINITION.'
c
c  Problem title.
c
      WRITE(*,*)' '
      WRITE(*,*)'ENTER A TITLE FOR THE PROBLEM.'
      READ(*,'(A)')TITLE
c
c  Get the number of nodes, and the values of the node
c  based quantities.
c
      WRITE(*,*)' '
      WRITE(*,*)'ENTER THE NUMBER OF NODES.'
      WRITE(*,*)'BETWEEN 2 AND ',MAXNOD
      READ(*,*)NUMNOD

      WRITE(*,*)' '
      WRITE(*,*)'ENTER x, y, z, HEAT SOURCE AND VOLUME'
      WRITE(*,*)'FOR EACH NODE.'

      DO I=1,NUMNOD
        READ(*,*)POSNOD(1,I),POSNOD(2,I),POSNOD(3,I),HOTNOD(I),VOLNOD(I)
      end do
c
c  Get links.
c
      WRITE(*,*)' '
      WRITE(*,*)'DEFINE LINKS BETWEEN PAIRS OF NODES, BY SPECIFYING'
      WRITE(*,*)'THE TWO NODES.  TO FINISH, USE 0, 0 FOR THE NODES.'

      NUMLNK=0

10    CONTINUE

      READ(*,*)INODE,JNODE

      IF(INODE.NE.0.AND.JNODE.NE.0)THEN
        NUMLNK=NUMLNK+1
        NODLNK(1,NUMLNK)=INODE
        NODLNK(2,NUMLNK)=JNODE
        GO TO 10
      end if
c
c  Get the value of link quantities.
c
      WRITE(*,*)' '
      WRITE(*,*)'ENTER VOLUME, CROSS-SECTIONAL AREA, AND 2 BOX VALUES'
      WRITE(*,*)'FOR EACH LINK.'
      DO I=1,NUMLNK
        READ(*,*)VOLLNK(I),CRSLNK(I),BOXLNK(1,I),BOXLNK(2,I)
      end do
c
c  Get constitutive quantities.
c
      WRITE(*,*)' '
      WRITE(*,*)'ENTER GAS-LIQUID FRICTION COEFFICIENT.'
      READ(*,*)FRICGL
      WRITE(*,*)' '
      WRITE(*,*)'ENTER GAS-WALL FRICTION COEFFICIENT.'
      READ(*,*)FRICGW
      WRITE(*,*)' '
      WRITE(*,*)'ENTER LIQUID-WALL FRICTION COEFFICIENT.'
      READ(*,*)FRICLW

      WRITE(*,*)' '
      WRITE(*,*)'ENTER 3 COMPONENTS OF GRAVITY VECTOR:'
      READ(*,*)GRAVIT(1),GRAVIT(2),GRAVIT(3)

      WRITE(*,*)'ENTER 0 FOR DONOR CELL, 1 FOR CELL AVERAGE.'
      READ(*,*)LNKOPT

      WRITE(*,*)'ENTER CRITICAL DENSITY, RHO-CRIT'
      READ(*,*)RHCRIT

      WRITE(*,*)'ENTER CRITICAL Sigma'
      READ(*,*)SCRIT

      WRITE(*,*)'ENTER CRITICAL W'
      READ(*,*)WCRIT
c
c  Compute number of equations
c
      NUMEQN=4*NUMNOD+2*NUMLNK
      TIME=0.0
      LIW=MAXIW
      LRW=MAXRW
      DO I=1,LIW
        IWORK(I)=0
      end do
      DO I=1,LRW
        RWORK(I)=0.0
      end do
      DO I=1,36
        ISAVE(I)=0
      end do
      DO I=1,218
        RSAVE(I)=0.0
      end do

      WRITE(*,*)' '
      WRITE(*,*)'ENTER INITIAL VALUES OF Y'
      WRITE(*,*)' '
      DO I=1,NUMEQN
        READ(*,*)Y(I)
      end do

      WRITE(*,*)' '
      WRITE(*,*)'ENTER INITIAL VALUES OF YDOT'
      WRITE(*,*)' '
      DO I=1,NUMEQN
        READ(*,*)YDOT(I)
      end do

      WRITE(*,*)' '
      WRITE(*,*)'ENTER 0 IF YOU BELIEVE YDOT IS NOT CONSISTENT WITH Y,'
      WRITE(*,*)'      1 OTHERWISE.'
      READ(*,*)ISTATE
      IF(ISTATE.NE.1)THEN
        ISTATE=0
      end if

      RETURN
      END
      SUBROUTINE JACNET(NUMEQN,TIME,Y,YDOT,ML,MU,PA,M0)

c*********************************************************************72
c
cc JACNET computes the jacobian of the residual.
c
c  Discussion:
c
c    This routine is used by LSODI.
c
c    It is optional, and LSODI can be told to expect it, or to 
c    ignore it.  In the current implementation, this routine is
c    simply a dummy entry, and LSODI approximates the necessary
c    jacobian information.
c
c    In general, the routine is supposed to compute the jacobian of the
c    residual, which has the form:
c
c      D RES/D Y= d ( G(Y,TIME)-A(Y,TIME)*YDOT) /D Y
c
      INTEGER M0
      INTEGER NUMEQN

      INTEGER ML
      INTEGER MU
      REAL PA(M0,NUMEQN)
      REAL TIME
      REAL Y(NUMEQN)
      REAL YDOT(NUMEQN)

      WRITE(*,*)' '
      WRITE(*,*)'JACNET - FATAL ERROR!'
      WRITE(*,*)'  THIS DUMMY ROUTINE SHOULD never BE CALLED!'

      STOP
      END
      SUBROUTINE NLMAP(LNKNOD,NABNOD,NODLNK,NUMLNK,NUMNOD)

c*********************************************************************72
c
cc NLMAP constructs the neighbor arrays NABNOD and LNKNOD.
c
      INTEGER NUMLNK
      INTEGER NUMNOD

      INTEGER I
      INTEGER J
      INTEGER K
      INTEGER LNKNOD(2,NUMLNK)
      INTEGER N1
      INTEGER N2
      INTEGER NABNOD(2,NUMNOD)
      INTEGER NLINK
      INTEGER NODLNK(2,NUMLNK)
c
c  NABNOD records how many links a given node starts or ends.
c
      DO I=1,NUMNOD
        NABNOD(1,I)=0
        NABNOD(2,I)=0
      end do

      DO I=1,NUMLNK
        N1=NODLNK(1,I)
        NABNOD(1,N1)=NABNOD(1,N1)+1
        N2=NODLNK(2,I)
        NABNOD(2,N2)=NABNOD(2,N2)+1
      end do
c
c  If we implicitly order the links by starting node, and then
c  by link number, LNKNOD(1,*) records this ordering.
c
      DO K=1,2

        NLINK=0

        DO I=1,NUMNOD

          DO J=1,NUMLNK
            IF(NODLNK(K,J).EQ.I)THEN
              NLINK=NLINK+1
              LNKNOD(K,NLINK)=J
            end if
          end do

        end do

      end do

      RETURN
      END
      SUBROUTINE POLY(BASE,COEF,NCOEF,PVALU)

c*********************************************************************72
c
cc POLY computes the value of a polynomial.
c
c  Discussion:
c
c    The polynomial has the form
c
c    POLY(X) = Sum (I=1 to NCOEF) COEF(I)*BASE**(I-1)
c
      INTEGER NCOEF

      REAL BASE
      REAL COEF(NCOEF)
      INTEGER I
      REAL PVALU

      PVALU=COEF(NCOEF)
      DO I=1,NCOEF-1
        PVALU=COEF(NCOEF-I)+BASE*PVALU
      end do

      RETURN
      END
      SUBROUTINE READER(BOXLNK,CRSLNK,FRICGL,FRICGW,FRICLW,GRAVIT,
     &  HOTNOD,ISAVE,ISTATE,IWORK,LIW,LNKOPT,LRW,
     &  MAXEQN,MAXIW,MAXLNK,MAXNOD,MAXRW,NODLNK,NUMEQN,NUMLNK,
     &  NUMNOD,POSNOD,
     &  RHCRIT,RSAVE,RWORK,SCRIT,TIME,TITLE,VOLLNK,VOLNOD,WCRIT,Y,YDOT)

c*********************************************************************72
c
cc READER reads a restart data defining a problem.
c
      INTEGER MAXIW
      INTEGER MAXEQN
      INTEGER MAXLNK
      INTEGER MAXNOD
      INTEGER MAXRW

      REAL BOXLNK(2,MAXLNK)
      REAL CRSLNK(MAXLNK)
      CHARACTER*80 FILER
      REAL FRICGL
      REAL FRICGW
      REAL FRICLW
      REAL GRAVIT(3)
      REAL HOTNOD(MAXNOD)
      INTEGER I
      INTEGER ISAVE(36)
      INTEGER ISTATE
      INTEGER IWORK(MAXIW)
      INTEGER JOB
      INTEGER LIW
      INTEGER LNKOPT
      INTEGER LRW
      INTEGER NODLNK(2,MAXLNK)
      INTEGER NUMEQN
      INTEGER NUMLNK
      INTEGER NUMNOD
      REAL POSNOD(3,MAXNOD)
      REAL RHCRIT
      REAL RSAVE(218)
      REAL RWORK(MAXRW)
      REAL SCRIT
      REAL TIME
      CHARACTER*80 TITLE
      REAL VOLLNK(MAXLNK)
      REAL VOLNOD(MAXNOD)
      REAL WCRIT
      REAL Y(MAXEQN)
      REAL YDOT(MAXEQN)

      WRITE(*,*)' '
      WRITE(*,*)'ENTER NAME OF RESTART FILE TO BE READ:'
      READ(*,'(A)')FILER

      OPEN(UNIT=2,FILE=FILER,STATUS='OLD',ERR=10)

      READ(2,'(A80)')TITLE

      READ(2,*)NUMEQN
      READ(2,*)NUMLNK
      READ(2,*)NUMNOD
      READ(2,*)LIW
      READ(2,*)LRW

      DO I=1,NUMLNK
        READ(2,*)BOXLNK(1,I)
      end do

      DO I=1,NUMLNK
        READ(2,*)BOXLNK(2,I)
      end do

      DO I=1,NUMLNK
        READ(2,*)CRSLNK(I)
      end do

      READ(2,*)FRICGL

      READ(2,*)FRICGW
      READ(2,*)FRICLW

      DO I=1,3
        READ(2,*)GRAVIT(I)
      end do

      DO I=1,NUMNOD
        READ(2,*)HOTNOD(I)
      end do

      DO I=1,36
        READ(2,*)ISAVE(I)
      end do

      READ(2,*)ISTATE

      DO I=1,LIW
        READ(2,*)IWORK(I)
      end do

      READ(2,*)LNKOPT

      DO I=1,NUMLNK
        READ(2,*)NODLNK(1,I)
      end do

      DO I=1,NUMLNK
        READ(2,*)NODLNK(2,I)
      end do

      DO I=1,NUMNOD
        READ(2,*)POSNOD(1,I)
      end do

      DO I=1,NUMNOD
        READ(2,*)POSNOD(2,I)
      end do

      DO I=1,NUMNOD
        READ(2,*)POSNOD(3,I)
      end do

      READ(2,*)RHCRIT

      DO I=1,218
        READ(2,*)RSAVE(I)
      end do

      DO I=1,LRW
        READ(2,*)RWORK(I)
      end do

      READ(2,*)SCRIT
      READ(2,*)TIME

      DO I=1,NUMLNK
        READ(2,*)VOLLNK(I)
      end do

      DO I=1,NUMNOD
        READ(2,*)VOLNOD(I)
      end do

      READ(2,*)WCRIT

      DO I=1,NUMEQN
        READ(2,*)Y(I)
      end do

      DO I=1,NUMEQN
        READ(2,*)YDOT(I)
      end do

      CLOSE(UNIT=2)
c
c  Pass new data to LSODI.
c
       JOB = 2
       CALL SRCOM(RSAVE,ISAVE,JOB)

c     CALL RSCOMI(RSAVE,ISAVE)

      RETURN
c
c  Error opening the file.
c
10    CONTINUE
      WRITE(*,*)' '
      WRITE(*,*)'READER - FATAL ERROR!'
      WRITE(*,*)'  THE INPUT FILE COULD NOT BE READ!'
      STOP
      END
      SUBROUTINE RESNET(N,TIME,Y,YDOT,RES,IRES)

c*********************************************************************72
c
cc RESNET computes the residual.
c
c  Discussion:
c
c    The residual has the form:
c
c    RES = G(Y,T) - A(Y,T) * dY/dT
c
      INTEGER MAXLNK
      INTEGER MAXNOD
      REAL PCON1
      REAL PCON2
      REAL SQRTPI

      PARAMETER (MAXLNK=15)
      PARAMETER (MAXNOD=15)
      PARAMETER (PCON1=0.18509)
      PARAMETER (PCON2=4629.168)
      PARAMETER (SQRTPI=1.77245)

      INTEGER N

      REAL RES(N)
      REAL Y(N)
      REAL YDOT(N)

      REAL ALGAS
      REAL ALGAS1
      REAL ALGAS2
      REAL ALGASL
      REAL ALGLNK(MAXLNK)
      REAL ALLIQ
      REAL ALLIQL
      REAL ALLLNK(MAXLNK)
      REAL AMGNOD(MAXNOD)
      REAL AMLNOD(MAXNOD)
      REAL BOX1
      REAL BOX2
      REAL BOXLNK(2,MAXLNK)
      REAL CRSLNK(MAXLNK)
      REAL DRGDH
      REAL DRGDH1
      REAL DRGDH2
      REAL DRGDP
      REAL DRGDP1
      REAL DRGDP2
      REAL DRLDH
      REAL DRLDH1
      REAL DRLDH2
      REAL DRLDP
      REAL DRLDP1
      REAL DRLDP2
      REAL FGL
      REAL FGW
      REAL FLW
      REAL FRICGL
      REAL FRICGW
      REAL FRICLW
      REAL GAMMA1
      REAL GAMMA2
      REAL GAMNOD(MAXNOD)
      REAL GRAVIT(3)
      REAL HGAS
      REAL HGAS1
      REAL HGAS2
      REAL HGASL
      REAL HGASS
      REAL HGLNK(MAXLNK)
      REAL HLIQ
      REAL HLIQ1
      REAL HLIQ2
      REAL HLIQL
      REAL HLIQS
      REAL HLLNK(MAXLNK)
      REAL HOTNOD(MAXNOD)
      INTEGER I
      INTEGER IALGAS
      INTEGER IEQN
      INTEGER IGDONR
      INTEGER IHGAS
      INTEGER IHLIQ
      INTEGER ILDONR
      INTEGER IPRESS
      INTEGER IRES
      INTEGER IVGAS
      INTEGER IVLIQ
      INTEGER IX1
      INTEGER IX2
      INTEGER J
      INTEGER JALGAS
      INTEGER JHGAS
      INTEGER JHLIQ
      INTEGER JPRESS
      INTEGER KALGAS
      INTEGER KHGAS
      INTEGER KHLIQ
      INTEGER KPRESS
      INTEGER LINK
      INTEGER LNKNOD(2,MAXLNK)
      INTEGER LNKOPT
      INTEGER NABNOD(2,MAXNOD)
      INTEGER NODE1
      INTEGER NODE2
      INTEGER NODLNK(2,MAXLNK)
      INTEGER NUMEQN
      INTEGER NUMLNK
      INTEGER NUMNOD
      REAL PLINK
      REAL POSNOD(3,MAXNOD)
      REAL POTEN1
      REAL POTEN2
      REAL POTNOD(MAXNOD)
      REAL PR1LNK(MAXLNK)
      REAL PR2LNK(MAXLNK)
      REAL PRESS
      REAL PRESS1
      REAL PRESS2
      REAL PRSLNK(MAXLNK)
      REAL RGHLNK(MAXLNK)
      REAL RGHNOD(MAXNOD)
      REAL RGPLNK(MAXLNK)
      REAL RGPNOD(MAXNOD)
      REAL RHCRIT
      REAL RHGAS
      REAL RHGAS1
      REAL RHGAS2
      REAL RHGASL
      REAL RHGLNK(MAXLNK)
      REAL RHGNOD(MAXNOD)
      REAL RHLIQ
      REAL RHLIQ1
      REAL RHLIQ2
      REAL RHLIQL
      REAL RHLLNK(MAXLNK)
      REAL RHLNOD(MAXNOD)
      REAL RLHLNK(MAXLNK)
      REAL RLHNOD(MAXNOD)
      REAL RLPLNK(MAXLNK)
      REAL RLPNOD(MAXNOD)
      REAL SCRIT
      REAL SUM1
      REAL SUM2
      REAL TIME
      REAL VGAS
      REAL VLIQ
      REAL VOLLNK(MAXLNK)
      REAL VOLNOD(MAXNOD)
      REAL WCRIT

      COMMON /PROBLM/
     &  ALGLNK,ALLLNK,AMGNOD,AMLNOD,BOXLNK,CRSLNK,FRICGL,FRICGW,FRICLW,
     &  GAMNOD,GRAVIT,HGLNK,HLLNK,HOTNOD,LNKNOD,LNKOPT,NABNOD,NODLNK,
     &  NUMEQN,NUMLNK,NUMNOD,POSNOD,POTNOD,PRSLNK,RGHLNK,RGHNOD,
     &  RGPLNK,RGPNOD,RHCRIT,RHGLNK,RHGNOD,RHLLNK,RHLNOD,RLHLNK,
     &  RLHNOD,RLPLNK,RLPNOD,SCRIT,VOLLNK,VOLNOD,WCRIT,PR1LNK,PR2LNK
c
c   node based physical quantities
c
      DO I=1,NUMNOD

        IALGAS=I
        IPRESS=IALGAS+NUMNOD
        IHGAS=IPRESS+NUMNOD
        IHLIQ=IHGAS+NUMNOD

        ALGAS=Y(IALGAS)
        ALLIQ=1.0-ALGAS
        PRESS=Y(IPRESS)
        HGAS=Y(IHGAS)
        HLIQ=Y(IHLIQ)
c
c  Evaluate RHGAS, DRGDP, DRGDH, HGASS,
c           RHLIQ, DRLDP, DRLDH, HLIQS.
c
        CALL DENSE(DRGDH,DRGDP,DRLDH,DRLDP,HGAS,HGASS,HLIQ,HLIQS,PRESS,
     &    RHGAS,RHLIQ)

        RHGNOD(I)=RHGAS
        RHLNOD(I)=RHLIQ
        RGPNOD(I)=DRGDP
        RGHNOD(I)=DRGDH
        RLPNOD(I)=DRLDP
        RLHNOD(I)=DRLDH
c
c    should the 2nd HLIQ be HGAS (elseif) ???
c
        IF(HOTNOD(I).GT.0.0.AND.HLIQ.GE.HLIQS)THEN
          GAMNOD(I)=ALLIQ*HOTNOD(I)/(HGASS-HLIQS)
        ELSEIF(HOTNOD(I).LT.0.0.AND.HLIQ.LE.HGASS)THEN
          GAMNOD(I)=ALGAS*HOTNOD(I)/(HGASS-HLIQS)
        ELSE
          GAMNOD(I)=0.0
        end if

        POTNOD(I)=0.0
        DO J=1,3
          POTNOD(I)=POTNOD(I)+GRAVIT(J)*POSNOD(J,I)
        end do

      end do
c
c  Get pressure and other variables at every link.
c
      DO I=1,NUMLNK
        IVGAS=4*NUMNOD+I
        IVLIQ=IVGAS+NUMLNK
        VGAS=Y(IVGAS)
        VLIQ=Y(IVLIQ)
c
c  Get the two end nodes of link I.
c
        NODE1=NODLNK(1,I)
        NODE2=NODLNK(2,I)
        BOX1=BOXLNK(1,I)
        BOX2=BOXLNK(2,I)
c
c  Get the indices and values of node based quantities.
c
        JALGAS=NODE1
        JPRESS=JALGAS+NUMNOD
        JHGAS=JPRESS+NUMNOD
        JHLIQ=JHGAS+NUMNOD
        ALGAS1=Y(JALGAS)
        PRESS1=Y(JPRESS)
        HGAS1=Y(JHGAS)
        HLIQ1=Y(JHLIQ)
        RHGAS1=RHGNOD(NODE1)
        RHLIQ1=RHLNOD(NODE1)
        DRGDH1=RGHNOD(NODE1)
        DRGDP1=RGPNOD(NODE1)
        DRLDH1=RLHNOD(NODE1)
        DRLDP1=RLPNOD(NODE1)
c
c  Same for node 2.
c
        KALGAS=NODE2
        KPRESS=KALGAS+NUMNOD
        KHGAS=KPRESS+NUMNOD
        KHLIQ=KHGAS+NUMNOD
        ALGAS2=Y(KALGAS)
        PRESS2=Y(KPRESS)
        HGAS2=Y(KHGAS)
        HLIQ2=Y(KHLIQ)
        RHGAS2=RHGNOD(NODE2)
        RHLIQ2=RHLNOD(NODE2)
        DRGDH2=RGHNOD(NODE2)
        DRGDP2=RGPNOD(NODE2)
        DRLDH2=RLHNOD(NODE2)
        DRLDP2=RLPNOD(NODE2)
c
c  Evaluate link quantities which depend on nodal values.
c
        CALL DOLINK(ALGAS,ALGAS1,ALGAS2,ALLIQ,BOX1,BOX2,DRGDH,
     &    DRGDH1,DRGDH2,DRGDP,DRGDP1,DRGDP2,DRLDH,DRLDH1,DRLDH2,DRLDP,
     &    DRLDP1,DRLDP2,HGAS,HGAS1,HGAS2,HLIQ,HLIQ1,HLIQ2,LNKOPT,PRESS,
     &    PRESS1,PRESS2,RHGAS,RHGAS1,RHGAS2,RHLIQ,RHLIQ1,RHLIQ2,VGAS,
     &    VLIQ)

        ALGLNK(I)=ALGAS
        ALLLNK(I)=ALLIQ
        HGLNK(I)=HGAS
        HLLNK(I)=HLIQ
        PR1LNK(I)=PRESS1
        PR2LNK(I)=PRESS2
        PRSLNK(I)=PRESS
        RHGLNK(I)=RHGAS
        RHLLNK(I)=RHLIQ
        RGHLNK(I)=DRGDH
        RGPLNK(I)=DRGDP
        RLHLNK(I)=DRLDH
        RLPLNK(I)=DRLDP
      end do
c
c  Ready to compute the residual.
c
      IEQN=0
      IX1=0
      IX2=0

      DO I=1,NUMNOD
c
c  Read off the nodal values.
c
        IALGAS=I
        IPRESS=IALGAS+NUMNOD
        IHGAS=IPRESS+NUMNOD
        IHLIQ=IHGAS+NUMNOD

        ALGAS=Y(IALGAS)
        ALLIQ=1.0-ALGAS
        PRESS=Y(IPRESS)
        HGAS=Y(IHGAS)
        HLIQ=Y(IHLIQ)
        RHGAS=RHGNOD(I)
        RHLIQ=RHLNOD(I)
        DRGDP=RGPNOD(I)
        DRGDH=RGHNOD(I)
        DRLDP=RLPNOD(I)
        DRLDH=RLHNOD(I)
c
c  1.  Gas continuity equation at node I.
c
        SUM1=RHGAS*YDOT(IALGAS)+ALGAS*DRGDP*YDOT(IPRESS)
     &    +ALGAS*DRGDH*YDOT(IHGAS)

        SUM2=GAMNOD(I)

        DO J=1,NABNOD(1,I)

          LINK=LNKNOD(1,IX1+J)
          VGAS=Y(4*NUMNOD+LINK)

          IF(VGAS.NE.0.0)THEN
            ALGASL=ALGLNK(LINK)
            RHGASL=RHGLNK(LINK)
            SUM2=SUM2-CRSLNK(LINK)*ALGASL*RHGASL*VGAS/VOLNOD(I)
          end if

        end do

        DO J=1,NABNOD(2,I)

          LINK=LNKNOD(2,IX2+J)
          VGAS=Y(4*NUMNOD+LINK)

          IF(VGAS.NE.0.0)THEN
            ALGASL=ALGLNK(LINK)
            RHGASL=RHGLNK(LINK)
            SUM2=SUM2+CRSLNK(LINK)*ALGASL*RHGASL*VGAS/VOLNOD(I)
          end if

        end do

        AMGNOD(I)=SUM2

        IEQN=IEQN+1
        RES(IEQN)=SUM2-SUM1
c
c  2.  Liquid continuity equation at node I.
c
        SUM1=-RHLIQ*YDOT(IALGAS)+ALLIQ*DRLDP*YDOT(IPRESS)
     &    +ALLIQ*DRLDH*YDOT(IHLIQ)

        SUM2=-GAMNOD(I)

        DO J=1,NABNOD(1,I)

          LINK=LNKNOD(1,IX1+J)
          VLIQ=Y(4*NUMNOD+NUMLNK+LINK)

          IF(VLIQ.NE.0.0)THEN
            ALLIQL=ALLLNK(LINK)
            RHLIQL=RHLLNK(LINK)
            SUM2=SUM2-CRSLNK(LINK)*ALLIQL*RHLIQL*VLIQ/VOLNOD(I)
          end if

        end do

        DO J=1,NABNOD(2,I)

          LINK=LNKNOD(2,IX2+J)
          VGAS=Y(4*NUMNOD+NUMLNK+LINK)

          IF(VLIQ.NE.0.0)THEN
            ALLIQL=ALLLNK(LINK)
            RHLIQL=RHLLNK(LINK)
            SUM2=SUM2+CRSLNK(LINK)*ALLIQL*RHLIQL*VLIQ/VOLNOD(I)
          end if

        end do

        AMLNOD(I)=SUM2

        IEQN=IEQN+1
        RES(IEQN)=SUM2-SUM1
c
c  3.  Gas energy equation at node I.
c
        SUM1=-PCON1*PRESS*YDOT(IALGAS)-ALGAS*PCON1*YDOT(IPRESS)
     &    +ALGAS*RHGAS*YDOT(IHGAS)
        SUM2=ALGAS*HOTNOD(I)+0.5*(HGAS+HLIQ-PCON1*PRESS*
     &    (RHGAS+RHLIQ)/(RHGAS*RHLIQ))*GAMNOD(I) - HGAS*AMGNOD(I)

        DO J=1,NABNOD(1,I)

          LINK=LNKNOD(1,IX1+J)
          VGAS=Y(4*NUMNOD+LINK)

          IF(VGAS.NE.0.0)THEN

            IF(LNKOPT.EQ.0)THEN
              IF(VGAS.GT.0.0)THEN
                PLINK=PR1LNK(LINK)
              ELSE
                PLINK=PR2LNK(LINK)
              end if
            ELSE
              PLINK=PRSLNK(LINK)
            end if

            ALGASL=ALGLNK(LINK)
            RHGASL=RHGLNK(LINK)
            HGASL=HGLNK(LINK)
            SUM2=SUM2-CRSLNK(LINK)*ALGASL*
     &        (RHGASL*HGASL+PCON1*(PRESS-PLINK))*VGAS/VOLNOD(I)
          end if

        end do

        DO J=1,NABNOD(2,I)

          LINK=LNKNOD(2,IX2+J)
          VGAS=Y(4*NUMNOD+LINK)

          IF(VGAS.NE.0.0)THEN

            IF(LNKOPT.EQ.0)THEN
              IF(VGAS.GT.0.0)THEN
                PLINK=PR1LNK(LINK)
              ELSE
                PLINK=PR2LNK(LINK)
              end if
            ELSE
              PLINK=PRSLNK(LINK)
            end if

            ALGASL=ALGLNK(LINK)
            RHGASL=RHGLNK(LINK)
            HGASL=HGLNK(LINK)
            SUM2=SUM2+CRSLNK(LINK)*ALGASL*
     &        (RHGASL*HGASL+PCON1*(PRESS-PLINK))*VGAS/VOLNOD(I)
          end if

        end do

        IEQN=IEQN+1
        RES(IEQN)=SUM2-SUM1
c
c  4.  Liquid energy equation at node I.
c
        SUM1=PCON1*PRESS*YDOT(IALGAS)-ALLIQ*PCON1*YDOT(IPRESS)
     &    +ALLIQ*RHLIQ*YDOT(IHLIQ)
        SUM2=ALLIQ*HOTNOD(I)-0.5*(HGAS+HLIQ-PCON1*PRESS*
     &    (RHGAS+RHLIQ)/(RHGAS*RHLIQ))*GAMNOD(I) - HLIQ*AMLNOD(I)

        DO J=1,NABNOD(1,I)

          LINK=LNKNOD(1,IX1+J)
          VLIQ=Y(4*NUMNOD+NUMLNK+LINK)

          IF(VLIQ.NE.0.0)THEN

            IF(LNKOPT.EQ.0)THEN
              IF(VLIQ.GT.0.0)THEN
                PLINK=PR1LNK(LINK)
              ELSE
                PLINK=PR2LNK(LINK)
              end if
            ELSE
              PLINK=PRSLNK(LINK)
            end if

            ALLIQL=ALLLNK(LINK)
            RHLIQL=RHLLNK(LINK)
            HLIQL=HLLNK(LINK)
            SUM2=SUM2-CRSLNK(LINK)*ALLIQL*
     &        (RHLIQL*HLIQL+PCON1*(PRESS-PLINK))*VLIQ/VOLNOD(I)
          end if

        end do

        DO J=1,NABNOD(2,I)

          LINK=LNKNOD(2,IX2+J)
          VLIQ=Y(4*NUMNOD+NUMLNK+LINK)

          IF(VLIQ.NE.0.0)THEN

            IF(LNKOPT.EQ.0)THEN
              IF(VLIQ.GT.0.0)THEN
                PLINK=PR1LNK(LINK)
              ELSE
                PLINK=PR2LNK(LINK)
              end if
            ELSE
              PLINK=PRSLNK(LINK)
            end if

            ALLIQL=ALLLNK(LINK)
            RHLIQL=RHLLNK(LINK)
            HLIQL=HLLNK(LINK)
            SUM2=SUM2+CRSLNK(LINK)*ALLIQL*
     &        (RHLIQL*HLIQL+PCON1*(PRESS-PLINK))*VLIQ/VOLNOD(I)
          end if

        end do

        IEQN=IEQN+1
        RES(IEQN)=SUM2-SUM1
        IX1 = IX1 + NABNOD(1,I)
        IX2 = IX2 + NABNOD(2,I)

      end do
c
c  Now do equations at links.
c
      DO I=1,NUMLNK

        IVGAS=4*NUMNOD+I
        IVLIQ=IVGAS+NUMLNK
        VGAS=Y(IVGAS)
        VLIQ=Y(IVLIQ)

        NODE1=NODLNK(1,I)
        NODE2=NODLNK(2,I)
        BOX1=BOXLNK(1,I)
        BOX2=BOXLNK(2,I)
        IF(VGAS.GE.0.0)THEN
          IGDONR=NODE1
        ELSE
          IGDONR=NODE2
        end if
        IF(VLIQ.GE.0.0)THEN
          ILDONR=NODE1
        ELSE
          ILDONR=NODE2
        end if
        JALGAS=NODE1
        JPRESS=JALGAS+NUMNOD
        JHGAS=JPRESS+NUMNOD
        JHLIQ=JHGAS+NUMNOD
        ALGAS1=Y(JALGAS)
        PRESS1=Y(JPRESS)
        HGAS1=Y(JHGAS)
        HLIQ1=Y(JHLIQ)
        RHGAS1=RHGNOD(NODE1)
        RHLIQ1=RHLNOD(NODE1)
        GAMMA1=GAMNOD(NODE1)
        POTEN1=POTNOD(NODE1)

        KALGAS=NODE2
        KPRESS=KALGAS+NUMNOD
        KHGAS=KPRESS+NUMNOD
        KHLIQ=KHGAS+NUMNOD
        ALGAS2=Y(KALGAS)
        PRESS2=Y(KPRESS)
        HGAS2=Y(KHGAS)
        HLIQ2=Y(KHLIQ)
        RHGAS2=RHGNOD(NODE2)
        RHLIQ2=RHLNOD(NODE2)
        GAMMA2=GAMNOD(NODE2)
        POTEN2=POTNOD(NODE2)

        ALGAS=ALGLNK(I)
        ALLIQ=ALLLNK(I)
        RHGAS=RHGLNK(I)
        RHLIQ=RHLLNK(I)
c
c  Evaluate friction terms.
c
        FGW=SQRTPI*RHGAS*FRICGW*ALGAS/(4.0*CRSLNK(I))
        FGL=3.0*RHCRIT*RHGAS*FRICGL*ALGAS*(VGAS-VLIQ)**2 /
     &    (4.0*0.06147*WCRIT*SCRIT)
        FLW=SQRTPI*RHLIQ*FRICLW*ALLIQ/(4.0*CRSLNK(I))
c
c  5.  Gas momentum equation at link I.
c
        SUM1=ALGAS*RHGAS*YDOT(IVGAS)

        SUM2=0.25*VLIQ*(GAMMA1+GAMMA2)-FGW*ABS(VGAS)*VGAS
     &    -FGL*ABS(VGAS-VLIQ)*(VGAS-VLIQ)

        SUM2=SUM2-ALGAS*CRSLNK(I)*(PCON2*(PRESS2-PRESS1)
     &    +RHGAS*(POTEN2-POTEN1))/VOLLNK(I)

        IF(LNKOPT.EQ.0)THEN
          SUM2=SUM2-0.5*VGAS*AMGNOD(IGDONR)
        ELSE
          SUM1=SUM1+0.5*VGAS*RHGAS*
     &      (BOX1*YDOT(JALGAS)+BOX2*YDOT(KALGAS))/(BOX1+BOX2)
          SUM1=SUM1+0.25*VGAS*ALGAS*RGPLNK(I)*
     &      (YDOT(JPRESS)+YDOT(KPRESS))
          IF(IGDONR.EQ.NODE1)THEN
            SUM1=SUM1+0.5*VGAS*ALGAS*RGHLNK(I)*YDOT(JHGAS)
          ELSE
            SUM1=SUM1+0.5*VGAS*ALGAS*RGHLNK(I)*YDOT(KHGAS)
          end if
        end if

        IEQN=IEQN+1
        RES(IEQN)=SUM2-SUM1
c
c  6.  Liquid momentum equation at link I.
c
        SUM1=ALLIQ*RHLIQ*YDOT(IVLIQ)

        SUM2=0.25*VGAS*(GAMMA1+GAMMA2)-FLW*ABS(VLIQ)*VLIQ
     &    +FGL*ABS(VGAS-VLIQ)*(VGAS-VLIQ)

        SUM2=SUM2-ALLIQ*CRSLNK(I)*(PCON2*(PRESS2-PRESS1)
     &    +RHLIQ*(POTEN2-POTEN1))/VOLLNK(I)

        IF(LNKOPT.EQ.0)THEN
          SUM2=SUM2-0.5*VLIQ*AMLNOD(ILDONR)
        ELSE
          SUM1=SUM1-0.5*VLIQ*RHLIQ*
     &      (BOX1*YDOT(JALGAS)+BOX2*YDOT(KALGAS))/(BOX1+BOX2)
          SUM1=SUM1+0.25*VLIQ*ALLIQ*RLPLNK(I)*
     &      (YDOT(JPRESS)+YDOT(KPRESS))
          IF(ILDONR.EQ.NODE1)THEN
            SUM1=SUM1+0.5*VLIQ*ALLIQ*RLHLNK(I)*YDOT(JHLIQ)
          ELSE
            SUM1=SUM1+0.5*VLIQ*ALLIQ*RLHLNK(I)*YDOT(KHLIQ)
          end if
        end if

        IEQN=IEQN+1
        RES(IEQN)=SUM2-SUM1

      end do

      RETURN
      END
      SUBROUTINE WRITER(BOXLNK,CRSLNK,FRICGL,FRICGW,FRICLW,GRAVIT,
     &  HOTNOD,ISAVE,ISTATE,IWORK,LIW,LNKOPT,LRW,
     &  MAXEQN,MAXLNK,MAXNOD,NODLNK,NUMEQN,NUMLNK,NUMNOD,POSNOD,
     &  RHCRIT,RSAVE,RWORK,SCRIT,TIME,TITLE,VOLLNK,VOLNOD,WCRIT,
     &  Y,YDOT)

c*********************************************************************72
c
cc WRITER writes restart data to a file.
c
      INTEGER LIW
      INTEGER LRW
      INTEGER MAXEQN
      INTEGER MAXLNK
      INTEGER MAXNOD

      REAL BOXLNK(2,MAXLNK)
      REAL CRSLNK(MAXLNK)
      CHARACTER*80 FILEW
      REAL FRICGL
      REAL FRICGW
      REAL FRICLW
      REAL GRAVIT(3)
      REAL HOTNOD(MAXNOD)
      INTEGER I
      INTEGER ISAVE(36)
      INTEGER ISTATE
      INTEGER IWORK(LIW)
      INTEGER JOB
      INTEGER LNKOPT
      INTEGER NODLNK(2,MAXLNK)
      INTEGER NUMEQN
      INTEGER NUMLNK
      INTEGER NUMNOD
      REAL POSNOD(3,MAXNOD)
      REAL RHCRIT
      REAL RSAVE(218)
      REAL RWORK(LRW)
      REAL SCRIT
      REAL TIME
      CHARACTER*80 TITLE
      REAL VOLLNK(MAXLNK)
      REAL VOLNOD(MAXNOD)
      REAL WCRIT
      REAL Y(MAXEQN)
      REAL YDOT(MAXEQN)

      WRITE(*,*)' '
      WRITE(*,*)'ENTER NAME OF RESTART FILE TO CREATE,'
      WRITE(*,*)'OR return FOR NO RESTART FILE.'
      READ(*,'(A)')FILEW
c
c  Get data from LSODI.
c
       JOB=1
       CALL SRCOM(RSAVE,ISAVE,JOB)

c     CALL SVCOMI(RSAVE,ISAVE)

      OPEN(UNIT=2,FILE=FILEW,STATUS='UNKNOWN',ERR=10)
      CLOSE(UNIT=2,STATUS='DELETE')

      OPEN(UNIT=2,FILE=FILEW,STATUS='NEW',ERR=10)

      WRITE(2,'(A80)')TITLE

      WRITE(2,*)NUMEQN
      WRITE(2,*)NUMLNK
      WRITE(2,*)NUMNOD
      WRITE(2,*)LIW
      WRITE(2,*)LRW

      DO I=1,NUMLNK
        WRITE(2,*)BOXLNK(1,I)
      end do

      DO I=1,NUMLNK
        WRITE(2,*)BOXLNK(2,I)
      end do

      DO I=1,NUMLNK
        WRITE(2,*)CRSLNK(I)
      end do

      WRITE(2,*)FRICGL
      WRITE(2,*)FRICGW
      WRITE(2,*)FRICLW

      DO I=1,3
        WRITE(2,*)GRAVIT(I)
      end do

      DO I=1,NUMNOD
        WRITE(2,*)HOTNOD(I)
      end do

      DO I=1,36
        WRITE(2,*)ISAVE(I)
      end do

      WRITE(2,*)ISTATE

      DO I=1,LIW
        WRITE(2,*)IWORK(I)
      end do

      WRITE(2,*)LNKOPT

      DO I=1,NUMLNK
        WRITE(2,*)NODLNK(1,I)
      end do

      DO I=1,NUMLNK
        WRITE(2,*)NODLNK(2,I)
      end do

      DO I=1,NUMNOD
        WRITE(2,*)POSNOD(1,I)
      end do

      DO I=1,NUMNOD
        WRITE(2,*)POSNOD(2,I)
      end do

      DO I=1,NUMNOD
        WRITE(2,*)POSNOD(3,I)
      end do

      WRITE(2,*)RHCRIT

      DO I=1,218
        WRITE(2,*)RSAVE(I)
      end do

      DO I=1,LRW
        WRITE(2,*)RWORK(I)
      end do

      WRITE(2,*)SCRIT

      WRITE(2,*)TIME

      DO I=1,NUMLNK
        WRITE(2,*)VOLLNK(I)
      end do

      DO I=1,NUMNOD
        WRITE(2,*)VOLNOD(I)
      end do

      WRITE(2,*)WCRIT

      DO I=1,NUMEQN
        WRITE(2,*)Y(I)
      end do

      DO I=1,NUMEQN
        WRITE(2,*)YDOT(I)
      end do

      CLOSE(UNIT=2)

      RETURN
c
c  Error opening the file.
c
10    CONTINUE
      WRITE(*,*)' '
      WRITE(*,*)'WRITER - WARNING!'
      WRITE(*,*)'  THE OUTPUT FILE COULD NOT BE CREATED.'
      WRITE(*,*)'  YOUR DATA WAS not SAVED!'
      RETURN
      END


