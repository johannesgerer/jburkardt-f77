      program main

c*********************************************************************72
c
cc MAIN is the main program for TOMS708_PRB1.
c
c      ALGORITHM 708, COLLECTED ALGORITHMS FROM ACM.
c      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
c      VOL. 18, NO. 3, SEPTEMBER, 1992, PP. 360-373z.
c
c     SAMPLE PROGRAM USING BRATIO. GIVEN THE NONNEGATIVE VALUES
c     A, B, X, Y WHERE A AND B ARE NOT BOTH 0 AND X + Y = 1. THEN
c
c              CALL BRATIO (A, B, X, Y, W, W1, IERR)
c
c     COMPUTES THE VALUES
c
c                W = I (A,B)  AND W1 = 1 - I (A,B).
c                     X                     X
c
c     IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.
c     IF NO INPUT ERRORS ARE DETECTED THEN IERR IS SET TO 0 AND
c     W AND W1 ARE COMPUTED. FOR MORE DETAILS SEE THE IN-LINE
c     DOCUMENTATION OF BRATIO.
c
      implicit none

      real a
      real b
      integer i
      integer ierr
      real w
      real w1
      real x
      real y

      write ( *, '(a)' ) ' '
      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS708_PRB1'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS708 library.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    X     Y           W              W1'
      write ( *, '(a)' ) ' '

      a = 5.3E+00
      b = 10.1E+00
      x = 0.01E+00

      do i = 1, 50

        y = 0.5E+00 + ( 0.5E+00 - x )

        call bratio ( a, b, x, y, w, w1, ierr )

        if ( ierr .ne. 0 ) then
          write ( *, '(f6.2,f6.2,a)' ) x, y, '  Error occurred.'
        else
          write ( *, '(f6.2,f6.2,e16.6,e16.6)') x, y, w, w1
        end if

        x = x + 0.01E+00

      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS708_PRB1'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
