C+FORMATB=1000
      Program X(OUTPUT,TAPE6=OUTPUT)
C
C  This is a comment which is (or was, anyway) mostly in lower case.
C
      CALL Y(I,J)
5     IF(I.GT.7)GoTo 2
      I=I+1
      CALL Z(J)
      GOTO 5
C
  2   Write(6,999)J
  999 FORMAT(1X,'The value of J is',I10)
      STOP
      END
C+LABELB=100
C+INDENTI=4
C+EXEMPTN
      SUBROUTINE   Y   (I,J)
      I=10.0*RANDOM(1.5)
      DO 1 K=1,I
      J=7.5*RANDOM(0.5)
      IF(J.GT.5)GOTO 2
      IF(J.LE.2)GOTO 1
      CALL Z(J)
1     CONTINUE
   2  CALL Z(J)
      RETURN
      END
