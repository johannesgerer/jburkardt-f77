      program main

c*********************************************************************72
c
cc MAIN is the main program for TOMS647_PRB.
c
c  Discussion:
c
c    TOMS647_PRB tests the TOMS647 routines.
c
c    Note that the SOBOL routine relies on a block data routine to
c    initialize common block memory.  This means that the SOBOL
c    routine cannot be properly reinitialized for a differente
c    dimension during a given run.  The code has to be restarted!
c    I started to fix this by replacing the block data by assignment
c    statements, but got disgusted.  For a corrected copy of the code,
c    see the FORTRAN90 version!
c
      implicit none

      integer atmost
      integer dimen
      double precision r8_epsilon
      integer seed
      double precision tiny

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS647_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS647 library.'

      atmost = 7000
      dimen = 5
      seed = 12345
      tiny = r8_epsilon ( )

      call test01 ( dimen, atmost )

      call test02 ( dimen, atmost, tiny )

      call test03 ( dimen, atmost )

      call test04 ( dimen, atmost, seed )

      atmost = 7000
      dimen = 10
      seed = 12345
      tiny = r8_epsilon ( )

      call test01 ( dimen, atmost )

      call test02 ( dimen, atmost, tiny )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  SKIPPING call to broken Sobol code.'

c     call test03 ( dimen, atmost )

      call test04 ( dimen, atmost, seed )

      atmost = 7000
      dimen = 20
      seed = 12345
      tiny = r8_epsilon ( )

      call test01 ( dimen, atmost )

      call test02 ( dimen, atmost, tiny )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  SKIPPING call to broken Sobol code.'

c     call test03 ( dimen, atmost )

      call test04 ( dimen, atmost, seed )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS647_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( dimen, atmost )

c*********************************************************************72
c
cc TEST01 tests GOFAUR.
c
C       THIS PROGRAM TESTS ACCURACY OF
C       NUMERICAL INTEGRATION USING "GOFAUR"
C       AND INTEGRAND (2) OF DAVIS AND
C       RABINOWITZ, PAGE 406
C
      LOGICAL FLAG(2)
      INTEGER DIMEN,ATMOST,I,J
      REAL QUASI(40),F,SUM

      write ( *, '(a)' )
      write ( *, '(a)' ) 'TEST FAURE'
      write ( *, '(a,i8)' ) 'DIMENSION = ',DIMEN
      write ( *, '(a,i8)' ) 'ATMOST = ',ATMOST

      CALL INFAUR(FLAG,DIMEN,ATMOST)

      IF ( .NOT. FLAG(1)) THEN
          write ( *, '(a,i8)' ) 'DIMENSION = ',DIMEN
          write ( *, '(a)' ) 'DIMEN IS NOT OK'
          return
      END IF

      IF ( .NOT. FLAG(2)) THEN
          write ( *,'(a)' ) 'ATMOST = ',ATMOST
          write ( *,'(a)' ) 'ATMOST IS NOT OK'
          return
      END IF

      write ( *,'(a)' ) ' '
      write ( *,'(a)' ) '  Iteration   Estimate'
      write ( *,'(a)' ) ' '

      SUM = 0.0
      DO I = 1, ATMOST
         CALL GOFAUR ( QUASI )
         F = 1.0
         DO J = 1,DIMEN
            F = F*ABS(4.0*QUASI(J)-2.0)
         end do
         SUM = SUM + F
         IF (MOD(I,500).EQ.0) THEN
           write ( *, '(2x,i8,2x,g14.6)' ) i, sum / i
         END IF
      end do

      return
      END
      subroutine test02 ( dimen, atmost, tiny )

c*********************************************************************72
c
cc TEST02 tests GOHALT.
c
C       THIS PROGRAM TESTS ACCURACY OF
C       NUMERICAL INTEGRATION USING "GOHALT"
C       AND INTEGRAND (2) OF DAVIS AND
C       RABINOWITZ, PAGE 406
C
      LOGICAL FLAG(2)
      INTEGER DIMEN,ATMOST,I,J
      double precision F,SUM
      DOUBLE PRECISION QUASI(40),TINY

      write ( *, '(a)' ) ' '
      write ( *,*) 'TEST HALTON'
      write ( *,*) 'DIMENSION = ',DIMEN
      write ( *,*) 'ATMOST = ',ATMOST
      write ( *,*) 'TINY = ',TINY

      CALL INHALT(FLAG,DIMEN,ATMOST,TINY,QUASI)

      IF ( .NOT. FLAG(1)) THEN
          write ( *,*) 'DIMENSION = ',DIMEN
          write ( *,*) 'DIMEN IS NOT OK'
          return
      END IF

      IF ( .NOT. FLAG(2)) THEN
          write ( *,*) 'ATMOST = ',ATMOST
          write ( *,*) 'ATMOST IS NOT OK'
          return
      END IF

      write ( *,'(a)' ) ' '
      write ( *,'(a)' ) '  Iteration   Estimate'
      write ( *,'(a)' ) ' '

      SUM = 0.0

      DO I = 1, ATMOST

         CALL GOHALT ( QUASI )

         F = 1.0
         DO J = 1, DIMEN
            F = F * ABS ( 4.0 * QUASI(J) - 2.0 )
         end do
         SUM = SUM + F

         IF (MOD(I,500).EQ.0) THEN
           write ( *, '(2x,i8,2x,g14.6)' ) i, sum / i
         END IF

      end do

      return
      END
      subroutine test03 ( dimen, atmost )

c*********************************************************************72
c
cc TEST03 tests GOSOBL.
c
C       THIS PROGRAM TESTS ACCURACY OF
C       NUMERICAL INTEGRATION USING "GOSOBL"
C       AND INTEGRAND (2) OF DAVIS AND
C       RABINOWITZ, PAGE 406
C
      LOGICAL FLAG(2)
      INTEGER DIMEN,ATMOST,I,J,TAUS
      REAL QUASI(40),F,SUM

      write ( *, '(a)' ) ' '
      write ( *, * ) 'TEST SOBOL'
      write ( *, * ) 'DIMENSION = ',DIMEN
      write ( *, * ) 'ATMOST = ',ATMOST

      call bdsobl

      call insobl ( flag, dimen, atmost, taus )

      write ( *, '(a,i8)' ) 'TAUS = ', TAUS

      IF ( .NOT. FLAG(1) ) THEN
          write ( *,*) 'DIMENSION = ',DIMEN
          write ( *,*) 'DIMEN IS NOT OK'
          return
      END IF

      IF ( .NOT. FLAG(2) ) THEN
          write ( *,*) 'ATMOST = ',ATMOST
          write ( *,*) 'ATMOST IS NOT OK'
          return
      END IF

      write ( *,'(a)' ) ' '
      write ( *,'(a)' ) '  Iteration   Estimate'
      write ( *,'(a)' ) ' '

      SUM = 0.0
      DO I = 1,ATMOST
         CALL GOSOBL(QUASI)
         F = 1.0
         DO J = 1,DIMEN
            F = F*ABS(4.0*QUASI(J)-2.0)
         end do
         SUM = SUM + F
         IF (MOD(I,500).EQ.0) THEN
            write ( *, '(2x,i8,2x,g14.6)' ) i, sum / i
         END IF
      end do

      return
      END
      subroutine test04 ( dimen, atmost, seed )

c*********************************************************************72
c
cc TEST04 tests UNIF.
c
C       THIS PROGRAM TESTS ACCURACY OF
C       NUMERICAL INTEGRATION USING "UNIF"
C       AND INTEGRAND (2) OF DAVIS AND
C       RABINOWITZ, PAGE 406
C
      INTEGER DIMEN,ATMOST,I,K
      integer seed
      REAL QUASI(40),F,SUM

      write ( *,'(a)' ) ' '
      write ( *,*) 'TEST UNIF'
      write ( *,*) 'DIMENSION = ',DIMEN
      write ( *,*) 'ATMOST = ',ATMOST
      write ( *,*) 'SEED = ',SEED

      write ( *,'(a)' ) ' '
      write ( *,'(a)' ) '  Iteration   Estimate'
      write ( *,'(a)' ) ' '

      SUM = 0.0
      DO I = 1,ATMOST
         F = 1.0
         DO K = 1,DIMEN
            QUASI(K) = UNIF ( seed )
            F = F*ABS(4.0*QUASI(K)-2.0)
         end do
         SUM = SUM + F
         IF (MOD(I,500).EQ.0) THEN
           write ( *, '(2x,i8,2x,g14.6)' ) i, sum / i
         END IF
      end do

      return
      end
