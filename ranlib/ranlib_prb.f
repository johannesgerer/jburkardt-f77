      program main

c*********************************************************************72
c
cc MAIN is the main program for RANLIB_PRB.
c
c  Modified:
c
c    10 December 2007
c
      character*100 phrase

      call test_bot

      phrase = 'randomizer'

      call test_chi ( phrase )

      stop
      end
      subroutine test_bot

c*********************************************************************72
C
cc TEST_BOT is a test program for the bottom level routines
C
c  Modified:
c
c    10 December 2007
c
C     Set up the random number generator
C     .. Local Scalars ..
      INTEGER ians,iblock,igen,iseed1,iseed2,itmp,ix,ixgen,nbad
C     ..
C     .. Local Arrays ..
      INTEGER answer(10000),genlst(5)
C     ..
C     .. External Functions ..
      INTEGER ignlgi
      EXTERNAL ignlgi
C     ..
C     .. External Subroutines ..
      EXTERNAL getsd,initgn,setall,setcgn
C     ..
C     .. Data statements ..
      DATA genlst/1,5,10,20,32/
C     ..
C     .. Executable Statements ..
      nbad = 0
      WRITE (*,9000)

 9000 FORMAT (' For five virual generators of the 32'/
     +       ' This test generates 10000 numbers then resets the block'/
     +       '      and does it again'/
     +       ' Any disagreements are reported -- there should be none'/)
C
C     Set up Generators
C
      CALL setall(12345,54321)
C
C     For a selected set of generators
C
      DO 60,ixgen = 1,5
          igen = genlst(ixgen)
          CALL setcgn(igen)
          WRITE (*,*) ' Testing generator ',igen
C
C     Use 10 blocks
C
          CALL initgn(-1)
          CALL getsd(iseed1,iseed2)
          DO 20,iblock = 1,10
C
C     Generate 1000 numbers
C
              DO 10,ians = 1,1000
                  ix = ians + (iblock-1)*1000
                  answer(ix) = ignlgi()
   10         CONTINUE
              CALL initgn(+1)
   20     CONTINUE
          CALL initgn(-1)
C
C     Do it again and compare answers
C
          CALL getsd(iseed1,iseed2)
C
C     Use 10 blocks
C
          DO 50,iblock = 1,10
C
C     Generate 1000 numbers
C
              DO 40,ians = 1,1000
                  ix = ians + (iblock-1)*1000
C      ANSWER( IX ) = IGNLGI()
                  itmp = ignlgi()
                  IF (.NOT. (itmp.NE.answer(ix))) GO TO 30
                  WRITE (*,9010) iblock,ians,ix,answer(ix),itmp

 9010             FORMAT (' Disagreement on regeneration of numbers'/
     +                   ' Block ',I2,' N within Block ',I2,
     +                   ' Index in answer ',I5/
     +                   ' Originally Generated ',I10,' Regenerated ',
     +                   I10)

                  nbad = nbad + 1
                  IF (nbad.GT.10) STOP ' More than 10 mismatches'
   30             CONTINUE
   40         CONTINUE
              CALL initgn(+1)
   50     CONTINUE
          WRITE (*,*) ' Finished testing generator ',igen
          WRITE (*,*) ' Test completed successfully'
   60 CONTINUE
      return
      END
      subroutine test_chi ( phrase )

c*********************************************************************72
c
cc TEST_CHI tests GENCHI, which generates Chi-Square deviates.
c
c  Modified:
c
c    10 December 2007
c
      real array(1000)
      real av
      real avtr
      real genchi
      real genunf
      real high
      integer i
      integer is1
      integer is2
      real low
      real param(1)
      character*(*) phrase
      character*4 type
      real var
      real vartr
      real xmax
      real xmin

      call phrtsd(phrase,is1,is2)

      call setall(is1,is2)

      write ( * , '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_CHI'
      write ( *, '(a)' ) '  Generate Chi-Square deviates.'

      low = 1.0
      high = 10.0
      param(1) = genunf ( low, high )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Parameters:'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  DF = ', param(1)
      write ( *, '(a)' ) ' '

      DO 30,i = 1,1000
          array(i) = genchi(param(1))
   30 CONTINUE
      CALL stats(array,1000,av,var,xmin,xmax)

      type = 'chis'
      CALL trstat(type,param,avtr,vartr)
      WRITE (*,9020) av,avtr,var,vartr,xmin,xmax

 9020 FORMAT (' Mean Generated: ',T30,G15.7,5X,'True:',T60,
     +       G15.7/' Variance Generated:',T30,G15.7,5X,'True:',T60,
     +       G15.7/' Minimun: ',T30,G15.7,5X,'Maximum:',T60,G15.7)

      return
      end
