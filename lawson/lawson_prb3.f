      program main

c*********************************************************************72
c
cc MAIN is the main program for LAWSON_PRB3.
c
C  DEMONSTRATE THE USE OF SUBROUTINE   SVDRS  TO COMPUTE THE    
C  SINGULAR VALUE DECOMPOSITION OF A MATRIX, A, AND SOLVE A LEAST    
C  SQUARES PROBLEM,  A*X=B.  
c
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
C     THE S.V.D.  A= U*(S,0)**T*V**T IS 
C     COMPUTED SO THAT..
C     (1)  U**T*B REPLACES THE M+1 ST. COL. OF  B.  
C   
C     (2)  U**T   REPLACES THE M BY M IDENTITY IN   
C     THE FIRST M COLS. OF  B.  
C   
C     (3) V REPLACES THE FIRST N ROWS AND COLS. 
C     OF A. 
C   
C     (4) THE DIAGONAL ENTRIES OF THE S MATRIX  
C     REPLACE THE FIRST N ENTRIES OF THE ARRAY S.   
C   
C     THE ARRAY S( ) MUST BE DIMENSIONED AT LEAST 3*N .     
C  
      integer MDA, MDB, MX
      parameter(MDA = 8, MDB = 8, MX = 8)
      integer I, J, KP1, KRANK, L, M, MINMN, MN1, MN2, N, NOISE
      double precision A(MDA, MX), AA(MDA, MX), ANOISE, B(MDB, MDB+1)
      double precision DN, DUMMY, FLOATN, GEN, ONE, RHO
      double precision S(MX), SM, SMALL3, SMALL4, SRSMSQ
      double precision T, TAU, TEN, WORK(2*MX), X(MX), ZERO
      parameter(ONE = 1.0d0, SMALL3 = 1.0d-3, SMALL4 = 1.0d-4)
      parameter(TEN = 10.0d0)
      parameter(ZERO = 0.0d0)

  170 format (/////'    M   N'/1X,2I4)
  180 format (/'  SINGULAR VALUES OF MATRIX')  
  190 format (/'  TRANSFORMED RIGHT SIDE, U**T*B')     
  200 format (/
     * '  ESTIMATED PARAMETERS, X=A**(+)*B, COMPUTED BY ''SVDRS''')
  210 format (/'  RESIDUAL VECTOR LENGTH = ',E12.4)     
  220 format (/
     * '  FROBENIUS NORM(A-U*(S,0)**T*V**T)'/
     * '  ---------------------------------  = ',g12.4/
     * '    SQRT(N) * SPECTRAL NORM OF A')
  230 format (2X,I5,E16.8,I5,E16.8,I5,E16.8)  
  240 format (/
     * ' PROG3.  THIS PROGRAM DEMONSTRATES THE ALGORITHM SVDRS.'/) 
  250 format (/
     * '  THE RELATIVE NOISE LEVEL OF THE GENERATED DATA WILL BE',
     * g12.4/
     * '  THE RELATIVE TOLERANCE FOR PSEUDORANK DETERMINATION IS RHO =',
     * g12.4)  
  260 format (/'  ABSOLUTE PSEUDORANK TOLERANCE, TAU = ',g12.4/
     * '  PSEUDORANK = ',I5)

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LAWSON_PRB3'

      DO 166 NOISE=1,2  
         if(NOISE .eq. 1) then
            ANOISE= ZERO
            RHO= SMALL3
         else
            ANOISE= SMALL4
            RHO= TEN * ANOISE
         endif

         write (*,240) 
         write (*,250) ANOISE,RHO  
c
C  INITIALIZE DATA GENERATION FUNCTION   
C   
         DUMMY=GEN(-ONE)

         DO 164 MN1=1,6,5  
            MN2=MN1+2 
            DO 162 M=MN1,MN2  
               DO 160 N=MN1,MN2  
                  write (*,170) M,N 
                  DO 35 I=1,M   
                     DO 20 J=1,M   
   20                   B(I,J)= ZERO
                     B(I,I)= ONE
                     DO 30 J=1,N   
                        A(I,J)=GEN(ANOISE)
                        AA(I,J)=A(I,J)    
   30                continue
   35             continue
                  DO 40 I=1,M   
   40                B(I,M+1)=GEN(ANOISE)  
C   
C  The arrays are now filled in..
C  Compute the singular value decomposition.
C
                  CALL SVDRS (A,MDA,M,N,B,MDB,M+1,S, WORK)    

                  write (*,180) 
                  write (*,230) (I,S(I),I=1,N)  
                  write (*,190) 
                  write (*,230) (I,B(I,M+1),I=1,M)  
C 
C  TEST FOR DISPARITY OF RATIO OF SINGULAR VALUES.   
C  
                  KRANK=N   
                  TAU=RHO*S(1)  
                  DO 50 I=1,N   
                     IF (S(I).LE.TAU) GO TO 60 
   50             CONTINUE  
                  GO TO 70  
   60             KRANK=I-1     
   70             continue
                  write (*,260) TAU,KRANK   
c
C  COMPUTE SOLUTION VECTOR ASSUMING PSEUDORANK IS KRANK  
C
                  DO 80 I=1,KRANK   
   80                B(I,M+1)=B(I,M+1)/S(I)
                  DO 100 I=1,N  
                     SM= ZERO
                     DO 90 J=1,KRANK   
   90                   SM=SM+A(I,J)*B(J,M+1)
  100                X(I)=SM   
c
C  COMPUTE PREDICTED NORM OF RESIDUAL VECTOR.
C
                  SRSMSQ= ZERO     
                  IF (KRANK .ne. M) then
                     KP1=KRANK+1   
                     DO 110 I=KP1,M
  110                   SRSMSQ=SRSMSQ+B(I,M+1)**2 
                     SRSMSQ=sqrt(SRSMSQ)   
                  endif
                  write (*,200) 
                  write (*,230) (I,X(I),I=1,N)  
                  write (*,210) SRSMSQ  
C  COMPUTE THE FROBENIUS NORM OF  A**T- V*(S,0)*U**T.    
C   
C  COMPUTE  V*S  FIRST.  
C   
                  MINMN=min(M,N)   
                  DO 130 J=1,MINMN  
                     DO 130 I=1,N  
                        A(I,J)=A(I,J)*S(J)
  130                continue
  135             continue
                  DN= ZERO
                  DO 150 J=1,M  
                     DO 145 I=1,N  
                        SM= ZERO
                        DO 140 L=1,MINMN  
  140                      SM=SM+A(I,L)*B(L,J)
c
C  COMPUTED DIFFERENCE OF (I,J) TH ENTRY     
C  OF  A**T-V*(S,0)*U**T.
C 
                        T=AA(J,I)-SM  
                        DN=DN+T**2
  145                continue
  150             continue
                  FLOATN = N                 
                  DN=sqrt(DN)/(sqrt(FLOATN)*S(1))     
                  write (*,220) DN  
  160          continue  
  162       continue  
  164    continue  
  166 continue

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LAWSON_PRB3:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop  
      end   
