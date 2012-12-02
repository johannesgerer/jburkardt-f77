C------------------------------------------------------------------------
C spacer         Des Higgins, EMBL Data Library, April 1992             I
C------------------------------------------------------------------------
C Program to read a euclidean distance matrix file and output
C a PCOORD (principal coordinates) analysis.
C The method was first described by John Gower in:
C  Gower,J.C. (1966) "Some distance properties of latent root and vector methods
C  used in multivariate analysis". Biometrika, vol. 53, pp. 325-328.
C
C This software was described in:
C  Higgins, D.G. (1992) "Sequence ordinations: a multivariate analysis approach
C  to analysing large sequence data sets". CABIOS, vol. 8, pp. 15-22.
C
       PROGRAM spacer
       IMPLICIT INTEGER(A-Z)
C MAXN is the maximum number of things that can be ordinated.
C Change this in EVERY subroutine that mentions it!  (Sorry about that, Des.)
       PARAMETER (MAXN=650)

       CHARACTER FNAME*50,FNAME2*50,ALPH(MAXN),LINE*120
       REAL SIMM(MAXN,MAXN)

C unit numbers for file I/O
       LIN   = 66
       LPCO  = 68
       LSYM  = 69

C set the symbols to a default of *
       DO 5 I = 1,MAXN
       ALPH(I) = '*'
5      CONTINUE

1      PRINT*
       FNAME = ' '
       CALL GETSTR('Name of file containing dist. matrix ? >> ',FNAME)
       OPEN(UNIT=LIN,FILE=FNAME,STATUS='OLD',ERR=1)
       READ(LIN,*) NSEQS       
       PRINT*
       PRINT*,'Number of objects  = ',NSEQS
       PRINT*
       READ(LIN,*) ((SIMM(I,J),J=1,I),I=1,NSEQS)

       DO 10 I = 1,NSEQS
       DO 10 J = 1,I
10     SIMM(J,I) = SIMM(I,J)

3      PRINT*
       FNAME2 = FNAME
C***       CALL CHFEXT(FNAME2,'PCOORD',NMLEN)
       CALL CHFEXT(FNAME2,'pcoord',NMLEN)

       CALL GETSTR('Name for P-Coord out. file ? [  '
     &   //FNAME2(1:NMLEN)//' ] >> ',FNAME2)
       OPEN(UNIT=LPCO,FILE=FNAME2,ERR=3,STATUS='NEW')

4      PRINT*
       print*,'You can use a file with one symbol per seq. '
       print*,'Ignore this option by pressing <RETURN>'
       print*
       FNAME2 = FNAME
C***       CALL CHFEXT(FNAME2,'symbs',NMLEN)
       CALL CHFEXT(FNAME2,'sym',NMLEN)

       CALL GETSTR('Name for symbol input file ? [  '
     &   //FNAME2(1:NMLEN)//' ] >> ',FNAME2)
       OPEN(UNIT=Lsym,FILE=FNAME2,ERR=20,STATUS='old')

C read the symbols from the symbol file if one was opened.
7      read(lsym,'(a)',end=20) line
       do 6 i = 1,len(line)
       if(line(i:i).ne.' ') then
         nsymbs = nsymbs + 1
         if(nsymbs.gt.nseqs) then
            print*
            print*,' Too many symbols in symbols file'
            print*
            go to 20
         end if
         alph(nsymbs) = line(i:i)
       end if
6      continue
       go to 7

C the number of eigen vectors to be calculated is NEIG
20     continue
       NEIG = MIN(10,NSEQS-1)
       
       PRINT*
       PRINT*,'Principal coordinates analysis'
       PRINT*

C set each euclidean distance (SIMM(i,j)) to be -1/2( simm(i,j) squared)
       DO 12 I = 1,NSEQS
       DO 12 J = 1,NSEQS
12     SIMM(I,J) = -0.5*(SIMM(I,J)*SIMM(I,J))

C get the eigen vectors and print the results
       CALL PCOORD(SIMM,ALPH,NSEQS,NEIG,LPCO)

       PRINT*
       PRINT*,'PCOORD done. using ',NEIG,' axes'
       PRINT*

       STOP
       END

C--------------------------------------------------------------------------
C CHFEXT                                                                  I
C--------------------------------------------------------------------------
C Take a file name (string*(*) ) and change the filename extension to '.ext'
       SUBROUTINE CHFEXT(FNAME,EXT,NEWLEN)
       IMPLICIT INTEGER(A-Z)
       CHARACTER FNAME*(*),EXT*(*)

       STRLEN = LEN(FNAME)
       DO 10 I = STRLEN,1,-1
       IF(FNAME(I:I).NE.' ') THEN      
         NAMLEN = I
         GO TO 11
       END IF
10     CONTINUE
11     CONTINUE

       STRLEN = LEN(EXT)
       DO 20 I = STRLEN,1,-1
       IF(EXT(I:I).NE.' ') THEN      
         EXTLEN = I
         GO TO 21
       END IF
20     CONTINUE
21     CONTINUE

       BRPOS  = INDEX(FNAME,']')
       NAMBEG = BRPOS + 1
       DOTPOS = INDEX(FNAME(NAMBEG:NAMLEN),'.')
       DOTPOS = DOTPOS + BRPOS
       IF(DOTPOS.EQ.0) THEN
         FNAME  = FNAME(1:NAMLEN)//'.'//EXT(1:EXTLEN)
         NEWLEN = NAMLEN + EXTLEN + 1
       ELSE
         FNAME  = FNAME(1:DOTPOS)//EXT(1:EXTLEN)
         NEWLEN = DOTPOS + EXTLEN
       END IF
   
       RETURN
       END

C------------------------------------------------------------------
C  GETSTR                             JAN. 1989                   I
C------------------------------------------------------------------
C SubROUTINE to prompt on the screen for a string using
C ANSTR to return the value if it is changed.
        SUBROUTINE GETSTR(PROMPT,ANSTR)
        IMPLICIT INTEGER(A-Z)
        CHARACTER*(*) PROMPT,ANSTR
        CHARACTER*70 TEMPSTR
        ANLEN= LEN(ANSTR)
        PLEN = LEN(PROMPT)
        DO 10 I = PLEN,1,-1
        IF(PROMPT(I:I).NE.' ') THEN
           PLEN = I
           GO TO 11
        END IF
10      CONTINUE

11      WRITE(*,'(1X,A,A$)') PROMPT(1:PLEN),' '
        TEMPSTR = ' '
        READ(*,'(A)',ERR=11) TEMPSTR(1:ANLEN)
        IF(TEMPSTR.EQ.' ') RETURN
        ANSTR = TEMPSTR(1:ANLEN)
        RETURN
        END

C----------------------------------------------------------------------------
C PCOORD     
C----------------------------------------------------------------------------
C take the maxnXmaxn array of distances (A) and get the NEIG eigen vectors
        SUBROUTINE PCOORD(A,ALPH,N,NEIG,LOUT)
        PARAMETER (MAXN=650,MAXEIG=12)
        REAL A(MAXN,MAXN),B(MAXN,MAXEIG)
        REAL E(MAXN),R(MAXN),xaxis(MAXN),yaxis(MAXN)
        REAL SUMALL,SUMPOS,SUMNEG
        CHARACTER ALPH(MAXN)

        DO 15 I = 1,N                             
15      E(I) = 0.0
        SUM =0.0

      DO 25 K = 1,N                                                             
      DO 20 J = 1,N                                                             
20    E(K)=E(K)+A(K,J)                                                          
25    SUM = SUM + E(K)                                                          

C     SUM IS MEAN FOR ALL VALUES                                                
C     VECTOR E IS MEAN OF ALL COLUMNS                                           
      DIV = 1./FLOAT(N)                                                         

      DO 26 J=1,N                                                               
26    E(J)=E(J)*DIV                                                             

      SUM = SUM/FLOAT(N**2)                                                     
      WRITE (LOUT,'(A)') ' '         

      DO 35 I = 1,N                                                             
      DO 35 J = 1,N                                                             
35    A(I,J)= A(I,J)-E(I)-E(J)+SUM                                              

      WRITE(LOUT,'(A)') ' PRINCIPAL COORDINATES ANALYSIS'
      WRITE(LOUT,'(A)') ' '
      WRITE(LOUT,'(A)')  
     & ' The method was first described in:'
      WRITE(LOUT,'(A)')  ' '
      WRITE(LOUT,'(A)')  
     & '   Gower,J.C. (1966)'
      WRITE(LOUT,'(A)')  
     & '   Some distance properties of latent root and vector'
      WRITE(LOUT,'(A)')  
     & '   methods used in multivariate analysis.'
      WRITE(LOUT,'(A)')  
     & '   Biometrika, vol. 53, pp. 325-328.'

      WRITE(LOUT,'(A)')  ' '
      WRITE(LOUT,'(A)')  
     & ' This software was described in:'
      WRITE(LOUT,'(A)')  ' '
      WRITE(LOUT,'(A)')  
     & '   Higgins, D.G. (1992)'
      WRITE(LOUT,'(A)')  
     & '   Sequence ordinations: a multivariate analysis approach'
      WRITE(LOUT,'(A)')  
     & '   to analysing large sequence data sets.'
      WRITE(LOUT,'(A)')  
     & '   CABIOS, vol. 8, pp. 15-22.'
      WRITE(LOUT,'(A)')  ' '
      WRITE(LOUT,'(A)')  ' '

      WRITE(LOUT,'(A)') '  DIAGONAL ELEMENTS OF TRANSFORMED MATRIX'
      WRITE(LOUT,'(A)') 
     &      '  Squared distance of each point from centroid'
      WRITE(LOUT,'(1X,5F14.6)') (A(I,I),I=1,N)         
      SUM = 0.0                                                                 

      DO 72 I=1,N                                                               
   72 SUM = SUM+A(I,I)                                                          

      WRITE(LOUT,'(A)') ' '
      WRITE(LOUT,'(A,F14.6)') '  TRACE (total variation ) = ',SUM   

C extract the eigen values (E) and eigen vectors (A after the call) of A
      CALL EIGEN(A,N,LOUT,E)

C sort the eigen values (E) in increasing order and the columns of A accordingly
      call eigensort(a,n,e)

      WRITE(LOUT,'(A)') ' '
      WRITE(LOUT,'(A)') '1 EIGENVALUES OF TRANSFORMED DISTANCE MATRIX'
      WRITE(LOUT,'(1X,5F14.6)') (E(J),J=N,1,-1)

        SUMALL = 0.0
        SUMPOS = 0.0
        SUMNEG = 0.0
        DO 10 I = 1,N
        SUMALL = SUMALL + E(I)
        IF(E(I).GT.0.0) SUMPOS = SUMPOS + E(I)
        if(e(i).lt.0.0) sumneg = sumneg + e(i)
10      CONTINUE

        WRITE(LOUT,'(A)') ' '
        WRITE(LOUT,'(A,F14.6)') 
     &  '  Sum of all eigen values      = ',SUMALL
        WRITE(LOUT,'(A,F14.6)') 
     &  '  Sum of positive eigen values = ',SUMPOS
        WRITE(LOUT,'(A,F14.6)') 
     &  '  Sum of negative eigen values = ',SUMNEG

      DO 80 I=n,n-neig+1,-1
      IF(E(I).LT.0.0) THEN
        neig = n-I
        GO TO 82
      ELSE IF(E(I).EQ.0.0) THEN
        R(I) = 0.0
      ELSE
        R(I) = SQRT(E(I))
      END IF
80    CONTINUE
82    CONTINUE

      DO 81 J=1,neig
      DO 81 I=1,N                         
81    B(I,J) = A(I,n-J+1)*R(n-J+1) 

      WRITE(LOUT,'(A)') ' '      
      WRITE(LOUT,'(A)') '1 SPECIMEN COORDINATES'

      DO 94 I=1,N                                                               
94    WRITE(LOUT,'(1X,I4,2X,12F10.5)') I, (B(I,J),J=1,neig)

      if(neig.ge.2) then
        write(lOut,'(a)') ' '
        WRITE(LOUT,'(A)') 
     &  '1 PLOT OF FIRST (horizontal) AND SECOND (vertical) AXES' 
        do 95 i = 1,n
        xaxis(i) = b(i,1)
        yaxis(i) = b(i,2)
95      continue
        call plotxy(xaxis,yaxis,n,alph,lOut)
      end if

      if(neig.ge.3) then
        write(lOut,'(a)') ' '
        WRITE(LOUT,'(A)') 
     &  '1 PLOT OF FIRST (horizontal) AND THIRD (vertical) AXES'  
        do 96 i = 1,n
        yaxis(i) = b(i,3)
96      continue
        call plotxy(xaxis,yaxis,n,alph,lOut)
      end if

      WRITE(LOUT,'(A)') ' '
      WRITE(LOUT,'(A)') ' ' 

        SEI=0.0
        write(lOut,'(a)') ' '
        WRITE(LOUT,'(A)') 
     &  '1        AXIS     PERCENTAGE VARIATION      CUMULATIVE'
        write(lOut,'(a)') ' '
        DO 97 I = N,N-neig+1,-1
        PERCY = E(I)*100.0/SUM
        SEI = SEI + PERCY
97      WRITE(LOUT,'(A,I4,7X,F10.2,15X,F10.2)')
     &  '       ',N-I+1,PERCY,SEI

      RETURN
      END


C-------------------------------------------------------------------------
C  EIGEN                                                                 I
C-------------------------------------------------------------------------
Calls 2 routines taken from 'Numerical Recipes'.  Input a real symmetric
C array (A) and get back the eigen vectors.  The array D contains all the
C eigen values of A on exit.  The nth element of D is the eigen value of
C the nth eigen vector which is the nth column of A.
        SUBROUTINE EIGEN(A,N,LOUT,D)
        PARAMETER(MAXN=650)
        REAL A(MAXN,MAXN),D(MAXN),E(MAXN)

        CALL TRED2(A,N,D,E)
        CALL TQLI(D,E,N,A)

        RETURN
        END

C-------------------------------------------------------------------------
C eigensort
C-------------------------------------------------------------------------
C sort the N eigenvalues in EVALUES in ascending order and sort the 
C eigen vectors
C in EVECTORS accordingly.  Uses a QUICKSORT adapted from Software Tools by 
C Kernighan and Plauger (1986), page 115.
        subroutine eigensort(evectors,n,evalues)
        implicit integer(a-z)
        parameter(MAXN=650)
        integer n,i,j,p,m
        real evectors(MAXN,MAXN),evalues(MAXN)
        real tempvector(MAXN),tempvalue,pivlin
        integer lst(50),ust(50)

        lst(1) = 1
        ust(1) = n
        p = 1

1       continue
        if(p.le.0) RETURN
        if(lst(p).ge.ust(p)) then
          p = p - 1
        else
          i = lst(p) - 1
          j = ust(p)
          pivlin = evalues(j)
2         continue
          if(i.lt.j) then
3           i = i + 1
            if(evalues(i).lt.pivlin) go to 3
4           j = j - 1
            if(evalues(j).le.pivlin) go to 5
            if(j.gt.i) go to 4
5           if(i.lt.j) then

C swap the 2 eigen values
            tempvalue = evalues(i)
            evalues(i) = evalues(j)
            evalues(j) = tempvalue
C swap the 2 vectors here
            do 20 m = 1,n 
            tempvalue     = evectors(m,i)
            evectors(m,i) = evectors(m,j)
            evectors(m,j) = tempvalue
20          continue

          end if
          go to 2
        end if
        j = ust(p)

            tempvalue  = evalues(i)
            evalues(i) = evalues(j)
            evalues(j) = tempvalue
C swap the two vectors here
            do 21 m = 1,n 
            tempvalue     = evectors(m,i)
            evectors(m,i) = evectors(m,j)
            evectors(m,j) = tempvalue
21          continue

          if(i-lst(p).lt.ust(p)-i) then
            lst(p+1) = lst(p)
            ust(p+1) = i - 1
            lst(p) = i + 1
          else
            lst(p+1) = i + 1
            ust(p+1) = ust(p)
            ust(p) = i - 1
          end if
          p = p + 1
        end if
        go to 1

        end

C-------------------------------------------------------------------------
C  tred2                                                                 I
C-------------------------------------------------------------------------
C Take the real symmetric array A and tridiagonalise it.   Translated from
C a C routine in 'Numerical Recipes'.
        SUBROUTINE TRED2(A,N,D,E)
        PARAMETER(MAXN=650)
        INTEGER I,J,K,L,N
        REAL A(MAXN,MAXN),D(MAXN),E(MAXN)
        REAL SCALE,HH,H,G,F

        DO 10 I = N,2,-1
        L = I - 1
        SCALE = 0.0
        H     = 0.0
        IF(L.GT.1) THEN
          DO 12 K = 1,L
          SCALE = SCALE + ABS(A(I,K))
12        CONTINUE
          IF(SCALE.EQ.0) THEN
            E(I) = A(I,L)
          ELSE
            DO 14 K = 1,L
            A(I,K) = A(I,K)/SCALE
            H = H + A(I,K)*A(I,K)
14          CONTINUE
            F = A(I,L)
            IF(F.GT.0) THEN
               G = -SQRT(H)
            ELSE
               G = SQRT(H)
            END IF
            E(I) = SCALE*G
            H = H - F*G
            A(I,L) = F-G
            F = 0.0
            DO 15 J = 1,L
            A(J,I) = A(I,J)/H
            G = 0.0
            DO 16 K = 1,J
16          G = G + A(J,K)*A(I,K)
            DO 17 K = J+1,L
17          G = G + A(K,J)*A(I,K)
            E(J) = G/H
            F = F + E(J)*A(I,J)
15          CONTINUE
            HH = F/(H+H)
            DO 18 J = 1,L
            F = A(I,J)
            G = E(J) - HH*F
            E(J) = G
            DO 19 K = 1,J
19          A(J,K) = A(J,K) - (F*E(K)+G*A(I,K))
18          CONTINUE            
          END IF
        ELSE
          E(I) = A(I,L)
        END IF
        D(I) = H
10      CONTINUE

        D(1) = 0.0
        E(1) = 0.0
        DO 30 I = 1,N
        L = I - 1
        IF(D(I).GT.0.0) THEN
          DO 32 J = 1,L
          G = 0.0
          DO 34 K = 1,L
34        G = G + A(I,K)*A(K,J)
          DO 36 K = 1,L
36        A(K,J) = A(K,J) - G*A(K,I)
32        CONTINUE
        END IF
        D(I) = A(I,I)
        A(I,I) = 1.0
        DO 38 J = 1,L
        A(I,J) = 0.0
        A(J,I) = 0.0
38      CONTINUE
30      CONTINUE

        RETURN
        END

C-------------------------------------------------------------------------
C  tqli                                                                  I
C-------------------------------------------------------------------------
C Take a real symmetric array that has been tridiagonalised by TRED2 and
C find all eigen values and eigen vectors.   Translated from a C routine in
C 'Numerical Recipes'.
        SUBROUTINE TQLI(D,E,N,Z)
        PARAMETER(MAXN=650)
        REAL D(MAXN),E(MAXN),Z(MAXN,MAXN)
        REAL S,R,P,G,F,DD,C,B,SIGN
        INTEGER N,M,L,ITER,I,K,IM

        DO 1 I = 2,N
1       E(I-1) = E(I)

        E(N) = 0.0

        DO 10 L = 1,N
        ITER = 0
9       CONTINUE
        DO 11 IM = L,N-1
        DD = ABS(D(IM)) + ABS(D(IM+1))
        IF(ABS(E(IM))+DD.EQ.DD) THEN
          M = IM
          GO TO 12
        END IF
11      CONTINUE
        M = N
12      CONTINUE
        IF(M.NE.L) THEN
          ITER = ITER + 1
          IF(ITER.GE.30) PRINT*,'ITER = ',ITER
C          IF(E(L).EQ.0.0) E(L) = 0.000001
          G = (D(L+1)-D(L))/(2.0*E(L))
          R = SQRT((G*G)+1.0)
          G = D(M) - D(L) + E(L)/(G+SIGN(R,G))
          C = 1.0
          S = 1.0
          P = 0.0
          DO 13 I = M-1,L,-1
          F = S*E(I)
          B = C*E(I)
          IF(ABS(F).GE.ABS(G)) THEN
            C = G/F
            R = SQRT((C*C)+1.0)
            E(I+1) = F*R
            S = 1.0/R
            C = C * S
          ELSE
            S = F/G
            R = SQRT((S*S)+1.0)
            E(I+1) = G*R
            C = 1.0/R
            S = S * C
          END IF
          G = D(I+1) - P
          R = (D(I)-G)*S + 2.0*C*B
          P = S*R
          D(I+1) = G + P
          G = C*R-B
          DO 14 K = 1,N
          F = Z(K,I+1)
          Z(K,I+1) = S*Z(K,I) + C*F
          Z(K,I) = C*Z(K,I) - S*F
14        CONTINUE
13        CONTINUE
          D(L) = D(L) - P
          E(L) = G
          E(M) = 0.0
        END IF
        IF(M.NE.L) GO TO 9
10      CONTINUE

        RETURN
        END

C-------------------------------------------------------------------------
C  SIGN                                                                  I
C-------------------------------------------------------------------------
C Used by routine TQLI above.   
        REAL FUNCTION SIGN(A,B)
        REAL A,B

        IF(B.LT.0.0) THEN
          SIGN = - ABS(A)
        ELSE
          SIGN = ABS(A)
        END IF

        RETURN
        END

C-------------------------------------------------------------------------
C PLOTXY
C-------------------------------------------------------------------------
C Scatterplot of the real numbers in the 2 arrays: 
C  X(hor.) and Y(vert.).
        subroutine plotxy(x,y,n,alph,lOut)
        implicit integer (a-z)
        parameter (LINELEN=100,PAGELEN=50)
        real x(*),y(*)
        real xmin,ymin,xmax,ymax,xrange,yrange,xscale,yscale
        character alph(*),graph(0:LINELEN,0:PAGELEN)
        character blank,dash,up,star

        BLANK = ' '
        DASH  = '-'
        UP    = '|'
        STAR  = '*'

        xmax = x(1)
        ymax = y(1)
        xmin = x(1)
        ymin = y(1)

        do 1 i = 0,LINELEN
        do 1 j = 0,PAGELEN
1       graph(i,j) = blank

        do 2 i = 1,n
        IF(X(I).GT.XMAX) XMAX = X(I)
        IF(X(I).LT.XMIN) XMIN = X(I)
        IF(Y(I).GT.YMAX) YMAX = Y(I)
        IF(Y(I).LT.YMIN) YMIN = Y(I)
2       CONTINUE

        xrange = xmax - xmin  
        yrange = ymax - ymin

        xscale = (linelen-1)/xrange
        yscale = (pagelen-1)/yrange

        yzero = 0
        if(ymin.le.0.0) yzero = int((- ymin) * yscale)
        xzero = 0
        if(xmin.le.0.0) xzero = int((- xmin) * xscale)

        do 4 i = 0,linelen
        graph(i,0)       = dash
        graph(i,yzero)   = dash
        graph(i,pagelen) = dash
4       continue

        do 3 i = 0,pagelen
        graph(0,i)       = up
        graph(xzero,i)   = up
        graph(linelen,i) = up
3       continue

        do 10 i = 1,n
        xc = 1 + int((x(i) - xmin) * xscale)
        yc = 1 + int((y(i) - ymin) * yscale)
        if(graph(xc,yc).eq.blank.or.
     &     graph(xc,yc).eq.dash.or.
     &     graph(xc,yc).eq.up) then
             graph(xc,yc) = alph(i)
        else
             graph(xc,yc) = star
        end if
10      continue

        write(lOut,'(a)') ' '
        do 20 i = pagelen,0,-1
        write(lOut,'(a,200a1)') ' ',(graph(j,i),j=0,linelen)
20      continue

        return
        end
