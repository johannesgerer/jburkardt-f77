      program main

c*********************************************************************72
c
cc MAIN is the main program for LAWSON_PRB05.
c
C  EXAMPLE OF THE USE OF SUBROUTINES BNDACC AND BNDSOL TO SOLVE 
C  SEQUENTIALLY THE BANDED LEAST SQUARES PROBLEM THAT ARISES IN     
C  SPLINE CURVE FITTING. 
C   
c  Modified:
c
c    19 October 2008
c
c  Author:
c
c    Charles Lawson, Richard Hanson
c
c  Reference:
c
c    Charles Lawson, Richard Hanson,
c    Solving Least Squares Problems,
c    SIAM, 1995,
c    ISBN: 0898713560,
c    LC: QA275.L38.
c
      integer mdg
      parameter ( mdg = 12 )
      integer mxy
      parameter ( mxy = 12 )
      integer nband
      parameter ( nband = 4 )

      double precision b(10), c(mxy), cov(mxy), g(mdg,5), h
      integer i
      integer ic
      integer ig
      integer ip
      integer ir
      integer j
      integer jt
      integer l
      integer m
      integer mt
      integer nbp
      integer nc
      double precision one, p1, p2, q(4), r, rdummy, rnorm
      double precision sigfac, sigsq, u, x(mxy), y(mxy), yfit, zero
      parameter(one = 1.0d0, zero = 0.0d0)
      data y /2.2d0, 4.0d0, 5.0d0, 4.6d0, 2.8d0, 2.7d0,
     *        3.8d0, 5.1d0, 6.1d0, 6.3d0, 5.0d0, 2.0d0/
 
  160 format (/'   I',8X,'X',10X,'Y',6X,'YFIT',4X,'R=Y-YFIT/1X') 
  170 format (1X,I3,4X,F6.0,4X,F6.2,4X,F6.2,4X,F8.4)
  180 format (/' C =',6F10.5/(4X,6F10.5))  
  190 format (3(2X,2I4,g15.7))     
  200 format (/' COVARIANCE MATRIX OF THE SPLINE COEFFICIENTS.')
  210 format (/' RNORM  =',g15.8)    
  220 format (/' SIGFAC =',g15.8)    

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LAWSON_PRB5'
      write ( *, '(a)' ) '  BNDACC accumulates a banded matrix.'
      write ( *, '(a)' ) '  BNDSOL solves an associated banded least'
      write ( *, '(a)' ) '    squares problem.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Execute a sequence of cubic spline fits'
      write ( *, '(a)' ) '  to a discrete set of data.'
    
      M=MDG 
c
c  Set the data points.
c
      do i = 1, m
        x(i) = dble ( 2 * i )
      end do
C   
C  BEGIN LOOP THRU CASES USING INCREASING NOS OF BREAKPOINTS.
C   
      DO 150 NBP = 5, 10  

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  The number of breakpoints is ', nbp

        NC = NBP + 2     
c
C  SET BREAKPOINTS  
c
        B(1)=X(1)    
        B(NBP)=X(M)  
        H=(B(NBP)-B(1))/(NBP-1) 
        DO I=3,NBP   
          B(I-1)=B(I-2)+H     
        end do

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Breakpoints:'
        write ( *, '(a)' ) ' '
        write (*,'(5(2x,g14.6))' ) ( b(i), i = 1, nbp )     
C
C  INITIALIZE IR AND IP BEFORE FIRST CALL TO BNDACC.     
C   
        IR=1 
        IP=1 
        I=1  
        JT=1
 
   40   continue

        MT=0 

   50   CONTINUE     

        IF (X(I).GT.B(JT+1)) GO TO 60
C   
C  SET ROW  FOR ITH DATA POINT
C   
           U=(X(I)-B(JT))/H 
           IG=IR+MT     
           G(IG,1)=P1(ONE - U) 
           G(IG,2)=P2(ONE - U) 
           G(IG,3)=P2(U)
           G(IG,4)=P1(U)
           G(IG,5)=Y(I) 
           MT=MT+1  
           IF (I.EQ.M) GO TO 60     
           I=I+1
           GO TO 50     
C   
C  SEND BLOCK OF DATA TO PROCESSOR 
C   
   60      CONTINUE

           CALL BNDACC (G,MDG,NBAND,IP,IR,MT,JT)
           IF (I.EQ.M) GO TO 70     
           JT=JT+1  
           GO TO 40     
c
C  COMPUTE SOLUTION C()
c
   70      CONTINUE     
           CALL BNDSOL (1,G,MDG,NBAND,IP,IR,C,NC,RNORM)  
c   
C  WRITE SOLUTION COEFFICIENTS  
c
           write ( *, '(a)' ) ' '
           write ( *, '(a)' ) 'C:'
           write ( *, '(a)' ) ' '
           write ( *, '(5(2x,g14.6))' ) ( c(l), l = 1, nc )

           write (*,210) RNORM  
C   
C  COMPUTE AND PRINT X,Y,YFIT,R=Y-YFIT  
C   
           write (*,160)
           JT=1 
                DO 110 I=1,M
   80           IF (X(I).LE.B(JT+1)) GO TO 90   
                JT=JT+1 
                GO TO 80

   90           U=(X(I)-B(JT))/H    
                Q(1)=P1(ONE - U)   
                Q(2)=P2(ONE - U)   
                Q(3)=P2(U)  
                Q(4)=P1(U)  
                YFIT=ZERO   
                     DO 100 L=1,4   
                     IC=JT-1+L  
  100                YFIT=YFIT+C(IC)*Q(L)   
                R=Y(I)-YFIT 
                write (*,170) I,X(I),Y(I),YFIT,R
  110           CONTINUE
C   
C  COMPUTE RESIDUAL VECTOR NORM. 
C   
           IF (M.LE.NC) GO TO 150   
           SIGSQ=(RNORM**2)/(M-NC) 
           SIGFAC=sqrt(SIGSQ)   
           write (*,220) SIGFAC     
           write (*,200)
C   
C  COMPUTE AND PRINT COLS. OF COVARIANCE.    
C   
                DO 140 J=1,NC   
                     DO 120 I=1,NC  
  120                COV(I)=ZERO    
                COV(J)= ONE
                CALL BNDSOL (2,G,MDG,NBAND,IP,IR,COV,NC,RDUMMY) 
                CALL BNDSOL (3,G,MDG,NBAND,IP,IR,COV,NC,RDUMMY) 
C   
C  COMPUTE THE JTH COL. OF THE COVARIANCE MATRIX. 
c
                     DO 130 I=1,NC  
  130                COV(I)=COV(I)*SIGSQ
  140           write (*,190) (L,J,COV(L),L=1,NC)   
  150      CONTINUE 

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LAWSON_PRB5:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )
    
      STOP  
      END   
      double precision function P1(T)
c*********************************************************************72
c
cc P1
c
      double precision T
      P1 = 0.25d0 * T**2 * T  
      return
      end
      double precision function P2(T)
c*********************************************************************72
c
cc P2
c
      double precision T
      P2 = -(1.0d0 - T)**2 * (1.0d0 + T) * 0.75d0 + 1.0d0
      return
      end
