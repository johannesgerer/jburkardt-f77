CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     react_diffP1_pbc.f                                                      C  
C     Fortran program to calculate the 1D solution to the Finite Element      C
C     Method P1 for the reaction-diffusion equation of lambda-omega type.     C
C     Basically like react_diffP1.f, but with time dependent Neumann b.c.'s   C
C     corresponding to the ppw solution, or the homogenous Neumann            C
C     b.c.'s (initial data either ppw or exponentially decaying).             C 
C     10 June 2003 (NB P1 is P2 or P2* in my thesis)                          C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC         
C     Notes: For consistency with the Matlab plotting routine                 C
C     plot_P1.m and other Matlab programs node numbering                      C
C     (J=0 to BIGJ) and profile numbering (N=0 to BIGN) are shifted up one    C
C     unit. The user is prompted to enter the parameter values including      C
C     whether Gaussian or periodic plane wave initial data is to be           C
C     used, and the type of b.c.'s. Data from every 'S' th profile is kept    C
C     and written out to files U_P1.dat, V_P1.dat (corresponding to soln's u  C
C     and v respectively). If TMIN=0 then the 1st profile is the init. data.  C
C     Also, we only keep approx <= 1000 nodal data points per                 C
C     profile as more are not needed for viewing purposes.                    C
C     For this version of GNU FORTRAN there appears to be a limitation        C
C     on the 2-D complex array of 8191-by-8191 so a sparse tri-diagonal       C 
C     solver is used here (just stores the 3 diagonals of the coefficient     C
C     matrix A in 1-dimensional arrays SUB, DIAG and SUP).                    C
C     As arrays of type 'complex' must have components 'real', I could        C
C     not do the calculations in double precision.                            C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     There are two ways to use this program. If we use it to for plotting    C
C     purposes, then choose S appropriately, e.g. choose S such that e.g.     C
C     every 50th profile plotted. If the program is to be used to calculate   C
C     the normed space-time errors (see the Part II paper) then take          C
C     S = R = 1, TMIN = LMIN = 0, (actually doesn't matter) PPW initial       C 
C     data chosen.                                                            C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Variable declarations
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT NONE

      DOUBLE PRECISION P,LAM0,LAM1,OM0,OM1,TMAX,TMIN,DELX,DELT,MU,X
      DOUBLE PRECISION RMIN,RMAX,LMAX,LMIN,XI,SPEED
      DOUBLE PRECISION X1,X2,T,MAXE0,STOREE0,ERROR0,ERROR1,AA,BB
      DOUBLE PRECISION E0,E1,LAM,OM,RC
      INTEGER BIGJ,SMALLJ,J,BIGN,SMALLN,N,DIMJ,DIMN,COUNT,COL,ROW,MAXD
      INTEGER CHOICE,R,S,NR,ERR,BCS
      COMPLEX CI 
      CHARACTER*20 FRMT1,FRMT2

      PARAMETER (CI  = (0.0,1.0), MAXD=120000)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Array declarations
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Obviously for very fine discretisation in the x direction we can
C     increase MAXD. Notice that MAXD has 6 digits, thus on lines 5 and
C     50 in the format statement we have I6.6.

      COMPLEX C(1:MAXD),B(1:MAXD),DIAG(1:MAXD),SUB(1:MAXD),SUP(1:MAXD)
      DOUBLE PRECISION U0(1:MAXD),V0(1:MAXD),DATA(16),U(1:MAXD)

      EXTERNAL TRI

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     User inputs
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c$$$      PRINT *, 'Enter parameter p (real)  '
c$$$      READ *, P
c$$$      PRINT *, 'Enter the positive parameter lambda0 (real)   '
c$$$      READ *, LAM0
c$$$      PRINT *, 'Enter the positive parameter lambda1 (real)   '
c$$$      READ *, LAM1
c$$$      PRINT *, 'Enter the non-zero parameter omega0 (real)   '
c$$$      READ *, OM0
c$$$      PRINT *, 'Enter the non-zero parameter omega1 (real)   '
c$$$      READ *, OM1
c$$$      PRINT *, 'Enter max x value Lmax (real)   '
c$$$      READ *, LMAX
c$$$      PRINT *, 'Enter min x value Lmin (real)   '
c$$$      READ *, LMIN
c$$$      PRINT *, 'Enter max time value Tmax (real)   '
c$$$      READ *, TMAX
c$$$      PRINT *, 'Enter min time value Tmin (real)   '
c$$$      READ *, TMIN

      P = 2.0
      LAM0 = 1.0
      LAM1 = 1.0
      OM0 = 3.0
      OM1 = -3.0
      LMAX = 250.0
      LMIN = 0.0
      TMAX = 110.0
      TMIN = 75.0

c$$$      PRINT *, 'Enter delta x >= LMAX/120000  (real)   '
c$$$C       So if LMAX = 120 then smallest delta x is 0.001 (this can be changed)
c$$$      READ *, DELX
c$$$
c$$$      PRINT *, 'Enter delta t <= 1/(LAM0*(2*P+2)) (real)   '
c$$$      READ *, DELT
c$$$
c$$$      PRINT *, 'Enter 0 for PPW initial data & 1 otherwise  '
c$$$      READ *, CHOICE
c$$$
c$$$
c$$$      IF (CHOICE.EQ.0) THEN
c$$$
c$$$        PRINT *, 'Enter 1 for time dependent b.c.s and 0 otherwise    '
c$$$        READ *, BCS
c$$$
c$$$      END IF
c$$$
c$$$
c$$$      PRINT *, 'Enter S so we plot every Sth profile   '
c$$$      READ *, S 
c$$$      PRINT *, 'Enter 1 to calculate errors and 0 not to    '
c$$$      READ *, ERR  

C       Only enter 1 if S=R=1, TMIN=LMIN=0

      DELX = 1.0
      DELT = 0.0005
      CHOICE = 1
C      BCS = 0
      S = 1167
      ERR = 0

   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Some data that may change
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      MU  = DELT/(DELX**2)
      BIGJ = NINT(LMAX/DELX)
      SMALLJ = NINT(LMIN/DELX)
      BIGN = NINT(TMAX/DELT)
      SMALLN = NINT(TMIN/DELT)
      DIMJ  = BIGJ+1
      DIMN  = BIGN+1
      ERROR0 = 0.0
      ERROR1 = 0.0
      STOREE0 = 0.0

      NR = 1000

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Set up the 2 parameters for controlling how much data to keep
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Want to use for plotting purposes approx NR nodes from (SMALLJ + 1)
C     to DIMJ, step R (recall if LMIN = 0, then SMALLJ = 0)

      IF (BIGJ.LE.NR) THEN

         R = 1

         ELSE

         R  =  (BIGJ-SMALLJ)/NR

      END IF
         
C     PRINT *, 'S =', S
C     PRINT *,'(plot every S th profile from TMIN to TMAX)'
C     PRINT *, 'R =', R
C     PRINT *,'(plot every R th node per profile from LMIN to LMAX)'

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Open all files
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      OPEN(UNIT=9,FILE='u1d.txt',STATUS='UNKNOWN')
      OPEN(UNIT=15,FILE='v1d.txt',STATUS='UNKNOWN')
      OPEN(UNIT=13,FILE='data.txt',
     *                                                STATUS='UNKNOWN')

C     (Recall /extra/tmp/dma0mrg/ corresponds to my hard disc drive)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Initialization including initial data
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


C     Initialize array 'A' to zero

      DO 1 COUNT=1,DIMJ
            SUB(COUNT)=(0.0,0.0)
            SUP(COUNT)=(0.0,0.0)
            DIAG(COUNT)=(0.0,0.0)
 1    CONTINUE

C     Initialize profile vector

      IF (CHOICE.EQ.1) THEN

        DO 3 J = 1,DIMJ
        X = DBLE(J-1)*DELX
c$$$        U0(J) = 0.1*EXP(-0.8*(X**2))
c$$$        U0(J) = 0.1*EXP(-0.8*X)
c$$$        V0(J) = U0(J)
c$$$        C(J) = U0(J) + CI*V0(J)
              
           IF (J.EQ.1)  THEN
              U0(J) = 0.01
              V0(J) = U0(J)
              C(J) = U0(J) + CI*V0(J)
           ELSE
              U0(J) = 0.0
              V0(J) = U0(J)
              C(J) = U0(J) + CI*V0(J)
           END IF
 3      CONTINUE

      ELSE

C     Calculate amplitude of periodic plane waves (CHOICE = 0)

        RMIN = SQRT((2.0*LAM0*((OM1**2)+(LAM1**2)))/(LAM1*(2.0*(OM1**2)+
     *         3.0*(LAM1**2))))

        PRINT *, 'Rmin =', RMIN

        RMAX = SQRT(LAM0/LAM1)

        PRINT *, 'Rmax =', RMAX

        RC = (RMIN + RMAX)/2.0



        LAM = LAM0-LAM1*(RC**2)
        OM = OM0+OM1*(RC**2)
         
        PRINT *, 'Amplitude of PPW =  ',RC

        SPEED = -(OM0 + OM1*(RC**2))/SQRT((LAM0-LAM1*(RC**2)))

        PRINT *, 'Speed of PPW =  ',SPEED

        DO  4 J = 1,DIMJ
          X = DBLE(J-1)*DELX
          U0(J) = RC*COS(SQRT(LAM)*X)
          V0(J) = RC*SIN(SQRT(LAM)*X)
          C(J) = U0(J) + CI*V0(J)
 4      CONTINUE
      END IF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Write to file the initial profiles
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      IF (TMIN.EQ.0.0) THEN


C        Remember this is Dad's special code for writing
C        out an arbitrary number of elements in a row. The
C        (E15.8,1X) refers to the format of this data (followed
C        by a space) and the I6.6 is the format of DIMJ (we
C        can have a maximum of 6 digits and if less than 6 
C        digits the spaces are filled up from the left with
C        zeros)

         WRITE(FRMT1,5)DIMJ

 5       FORMAT("(",I6.6,"(E15.8,1X))")

         WRITE(9,FRMT1), (U0(J),J=(SMALLJ+1),DIMJ,R)
C           Ouputs initial profile for U
         WRITE(15,FRMT1), (V0(J),J=(SMALLJ+1),DIMJ,R)
C           Ouputs initial profile for V

      END IF

C     A simpler way to the above is the following. But have to keep 
C     changing the DIMJ each time (i.e. the 121 on line 5)

c$$$      IF (TMIN.EQ.0.0) THEN
c$$$
c$$$C         Write to file the initial profiles
c$$$
c$$$          WRITE(9,5), (U0(J),J=1,DIMJ)
c$$$C           Ouputs initial profile for U
c$$$          WRITE(15,5), (V0(J),J=1,DIMJ)
c$$$C           Ouputs approx. soln profiles for V
c$$$
c$$$ 5        FORMAT(121(E15.8,1X))
c$$$
c$$$      END IF
         
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Main part of program starts
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Now we loop forward from one time level to the next
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      DO 10 N=1,DIMN
C     To plot the very last profile change BIGN to DIMN


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Assemble coefficient matrix A
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Assemble parts of A that remain fixed

         DO 20 COUNT = 1,BIGJ-1
            SUP(1+COUNT) = -MU
            SUB(1+COUNT)   = -MU
 20      CONTINUE

         DO 30 COUNT = 1,DIMJ
            DIAG(COUNT) = 2.0*MU + 1.0
 30      CONTINUE

         SUP(1)         = -2.0*MU
         SUB(DIMJ) = -2.0*MU

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Update RHS and parts of A changing from level to level
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


         DO 40 J=1,DIMJ
            DIAG(J)=DIAG(J)-DELT*(LAM0-LAM1*(CABS(C(J)))**P
     *             + CI*(OM0+OM1*(CABS(C(J)))**P))

         IF (BCS.EQ.1) THEN
            X=DBLE(J-1)*DELX
            T=DBLE(N-1)*DELT
               IF (J.EQ.1) THEN
                 B(J)=C(J)-(2.0*DELT/DELX)*CI*RC*SQRT(LAM)*EXP(CI*OM*T)
               ELSE IF (J.EQ.DIMJ) THEN
                 B(J)=C(J)+(2.0*DELT/DELX)*CI*RC*SQRT(LAM)*EXP(CI*(OM*T+
     *              SQRT(LAM)*LMAX))
               ELSE
                 B(J) = C(J)
               END IF

         ELSE
                 B(J) = C(J)     
         END IF
     

 40      CONTINUE

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Solve tridiagonal system A*C=B for next time profile
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


         CALL TRI(BIGJ,DIMJ,BIGN,SUB,DIAG,SUP,B,C,MAXD)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Output profiles to file
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Output 1st profile and then every 's th' profile after time 'TMIN'
C     (recall we approx TMIN via DELT*(SMALLN+1))

         IF (((TMIN.NE.0.0).AND.(N.EQ.(SMALLN+1))).OR.
     *   ((N.GT.(SMALLN+1)).AND.(MOD((N-SMALLN-1),S).EQ.0))) THEN

              WRITE(FRMT2,50)DIMJ

 50           FORMAT("(",I6.6,"(E15.8,1X))")

              WRITE(9,FRMT2), (REAL(C(J)),J=(SMALLJ+1),DIMJ,R)
C                Ouputs approx. soln profiles for U
              WRITE(15,FRMT2), (AIMAG(C(J)),J=(SMALLJ+1),DIMJ,R)
C                Ouputs approx. soln profiles for V

                 ELSE

          END IF

c$$$ 45          WRITE(9,50), (REAL(C(J)),J=1,DIMJ)
c$$$C              Ouputs approx. soln profiles for U
c$$$             WRITE(15,50), (AIMAG(C(J)),J=1,DIMJ)
c$$$C              Ouputs approx. soln profiles for V
c$$$
c$$$ 50          FORMAT(121(E15.8,1X))
c$$$               
c$$$               ELSE
c$$$
c$$$         END IF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Calculate normed space-time errors
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        IF (ERR.EQ.1) THEN

          E0 = 0.0
          E1 = 0.0
          DO 60 J=1,DIMJ
            X1 = DBLE(J-1)*DELX  
C           i.e. x_j
            X2 = DBLE(J)*DELX    
C           i.e. x_{j+1}
            T = DBLE(N-1)*DELT
C           i.e. t_n
            U(J) = RC*COS((OM)*T + SQRT(LAM)*X1)
            U(J+1) = RC*COS((OM)*T + SQRT(LAM)*X2)
            AA = U(J+1)-REAL(C(J+1))
            BB = U(J)-REAL(C(J))
            E0 = (AA**2 + AA*BB + BB**2) + E0
C           recalling that (A^3-B^3)/(A-B) = A^2 + A*B + B^2
            E1 = (AA-BB)**2 + E1
         

 60       CONTINUE

          E0 = (DELX/3.0)*E0
          E1 = E1/DELX

          MAXE0 = MAX(E0,STOREE0)
          STOREE0 = MAXE0
          ERROR0 = E0 + ERROR0
          ERROR1 = E1 + ERROR1

        END IF

 10   CONTINUE

      IF (ERR.EQ.1) THEN

        PRINT *, 'L2L2_error =', DELT*ERROR0
        PRINT *, 'L2H1semi_error =', DELT*ERROR1
        PRINT *, 'L2H1error =', DELT*(ERROR1 + ERROR0) 
        PRINT *, 'LinfL2_error = ', MAXE0

      END IF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Output plotiing parameters
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Assign data (plotting parameters) to the array DATA. This will then be 
C     output to the file data_P1.dat, to be imported into the matlab plotting
C     routine plot_P1.m. This is really so we don't have to have the user
C     prompted twice to enter the parameters. Now because the array DATA is 'real'
C     to store integers we must convert them to real (and when importing them into
C     plot_P1.m we convert them back into integers again).
      
      DATA(1)=TMAX
      DATA(2)=TMIN
      DATA(3)=LMAX
      DATA(4)=LMIN
      DATA(5)=DELT
      DATA(6)=DELX
      DATA(7)=REAL(DIMJ)
      DATA(8)=REAL(DIMN)
      DATA(9)=REAL(S)
      DATA(10)=REAL(R)
      DATA(11)=RC
      DATA(12)=LAM0
      DATA(13)=LAM1
      DATA(14)=OM0
      DATA(15)=OM1
      DATA(16)=REAL(CHOICE)

      WRITE(13,1000), (DATA(J),J=1,16)
C        Outputs plotting parameters for the routine plot_P1.m
 1000 FORMAT(E15.8,1X)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Close all units
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      CLOSE(UNIT=9,STATUS='KEEP')
      CLOSE(UNIT=15,STATUS='KEEP')
      CLOSE(UNIT=13,STATUS='KEEP')

      STOP
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Main part of program ends
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Subroutine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      SUBROUTINE TRI(BIGJ,DIMJ,BIGN,SUB,DIAG,SUP,B,C,MAXD)
      
      IMPLICIT NONE

      INTEGER BIGJ,DIMJ,BIGN,COUNT,MAXD
      COMPLEX C(1:MAXD),B(1:MAXD),M,DIAG(1:MAXD),SUB(1:MAXD),SUP(1:MAXD)

C     Forward substitution

      DO 200 COUNT = 2,DIMJ
         M=SUB(COUNT)/DIAG(COUNT-1)
         DIAG(COUNT)=DIAG(COUNT)-M*SUP(COUNT-1)
         B(COUNT)=B(COUNT)-M*B(COUNT-1)
         SUB(COUNT)=0.0
 200  CONTINUE 

C     Back-substitution
      
      C(DIMJ)=B(DIMJ)/DIAG(DIMJ)

      DO 300 COUNT = BIGJ,1,-1
         C(COUNT)=(B(COUNT)-SUP(COUNT)*C(COUNT+1))
     *                 /DIAG(COUNT)

 300  CONTINUE   

      RETURN

      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
