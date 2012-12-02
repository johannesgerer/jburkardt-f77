      SUBROUTINE COUNT ( C, K, P, N )

c*********************************************************************72
c
C  COUNT COMPUTES THE NUMBER OF PARTITIONS OF AN INTEGER
C  RESTRICTED TO C FOR INTEGERS IN THE RANGE 1 TO N.
C  INPUT:  K -- A POSITIVE INTEGER.
C          C -- AN ARRAY OF K POSITIVE INTEGERS.
C          N -- AN INTEGER LARGER THAN THE MAXIMUM VALUE IN C.
C  OUTPUT: P -- AN ARRAY OF N INTEGERS, WHERE P(M) IS THE
C               NUMBER OF PARTITIONS OF M RESTRICTED TO C.
C
      INTEGER C, P
      DIMENSION C(K), P(N)
C
C  INITIALIZE P
C
      DO 10 I = 1, N
        P(I) = 0
10    CONTINUE
C
C  EACH PASS THROUGH THE OUTER LOOP BELOW TRANSFORMS P FROM
C  PARTITIONS RESTRICTED TO C(1), ..., C(I-1) TO
C  PARTITIONS RESTRICTED TO C(1), ..., C(I).
C
      DO 30 I = 1, K
        J = C(I)
        JP1 = J + 1
        P(J) = P(J) + 1
        DO 20 M = JP1, N
          MMJ = M - J
          P(M) = P(M) + P(MMJ)
20      CONTINUE
30    CONTINUE

      RETURN
      END
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
c  Modified:
c
c    16 September 2005
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

      character ( len = 8 ) date
      character ( len = 10 ) time

      call date_and_time ( date, time )

      write ( *, '(a8,2x,a10)' ) date, time

      return
      end
