      PROGRAM X (OUTPUT,TAPE6=OUTPUT)
C
C  THIS IS A COMMENT WHICH IS (OR WAS, ANYWAY) MOSTLY IN LOWER CASE.
C
      CALL Y (I,J)
   10 IF (I.GT.7) GO TO 20
      I = I+1
      CALL Z (J)
      GO TO 10
C
   20 WRITE (6,1000) J
 1000 FORMAT (1X,'THE VALUE OF J IS',I10)
      STOP
      END
      SUBROUTINE   Y   (I,J)
      I = 10.0*RANDOM(1.5)
      DO 100 K = 1, I
          J = 7.5*RANDOM(0.5)
          IF (J.GT.5) GO TO 110
          IF (J.LE.2) GO TO 100
          CALL Z (J)
  100 CONTINUE
  110 CALL Z (J)
      RETURN
      END
