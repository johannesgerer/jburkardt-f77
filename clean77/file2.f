C+UNCOND=C
      SUBROUTINE SETC (A,B,C)
C
      IF(A.GT.0.0)THEN
      IF(A.GT.B)THEN
      C=A
      ELSE
      C=B/A
      ENDIF
      ELSE
      C=-A/B
      ENDIF
      RETURN
      END
