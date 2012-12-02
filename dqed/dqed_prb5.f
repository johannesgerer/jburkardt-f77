c Test program for SLATEC routine DBOLS
c by Steve White

      program main
      implicit none
      integer   NCOLS, MROWS, MDW
c Test problem size
      parameter (NCOLS=3, MROWS=5) 
c
c  Row dimension of W.
c
      parameter (MDW=MROWS)
      integer   NX, NI   
c
c  No options in use.
c
      parameter (NX=0, NI=1)
      double precision W(MDW,NCOLS+1)
      double precision BL(NCOLS), BU(NCOLS)
      double precision X(NCOLS+NX)
      double precision RW(5*NCOLS)
      double precision RNORM
      integer   IND(NCOLS), IOPT(1+NI)
      integer   IW(2*NCOLS)
      integer   MODE
c
c  No special options.
c
      data      IOPT / 7, 99 /
c
c Test problem data
c
      data      W(1,1), W(1,2), W(1,3) /  4.0, -7.0, -7.0 /
      data      W(2,1), W(2,2), W(2,3) /  3.0, -4.0,  0.0 /
      data      W(3,1), W(3,2), W(3,3) / -7.0,  3.0, -3.0 /
      data      W(4,1), W(4,2), W(4,3) /  7.0,  7.0,  1.0 /
      data      W(5,1), W(5,2), W(5,3) / -1.0,  7.0,  6.0 /

      data      W(1,4), W(2,4), W(3,4), W(4,4), W(5,4) 
     &  / -22, -6, -2, 34, 13 /

      data      BL / NCOLS * 0.0 /
c
c  Use lower bound BL only.  BU not used.
c
      data      IND / NCOLS * 1 /
c
c  Ignore all bounds.
c
c      data      IND / NCOLS * 4 /

      print *, "Matrix ============================================="
      call colprint( W(1,1), MROWS )
      call colprint( W(1,2), MROWS )
      call colprint( W(1,3), MROWS )
      print *, "RHS ============================================="
      call colprint( W(1,4), MROWS )
      print *, "Bounds ============================================="

      write ( *, '(10g12.3)' ) bl  

      call DBOLS( W, MDW, MROWS, NCOLS, BL, BU, IND, IOPT, X, 
     &                RNORM, MODE, RW, IW )

      print *, "Solution ============================================="
      write ( *, '(10g12.3)' ) x

      print *, "MODE", MODE
      print *, "Residual", RNORM

      stop
      end
      subroutine colprint( A, MROWS )

      double precision A(MROWS)

      write ( *, '(10g12.3)' ) a

      return
      end
