      program main

c*********************************************************************72
c
cc MAIN tests the SPECIAL_FUNCTIONS library.
c
c  Modified:
c
c    12 April 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPECIAL_FUNCTIONS_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version.'
      write ( *, '(a)' ) '  Test the SPECIAL_FUNCTIONS library.'

      call mairya ( )
      call mairyb ( )
      call mairyzo ( )
      call maswfa ( )
      call maswfb ( )

      call mbernoa ( )
      call mbernob ( )
      call mbeta ( )

      call mcchg ( )
      call mcerror ( )
      call mcerzo ( )
      call mcgama ( )
      call mch12n ( )
      call mchgm ( )
      call mchgu ( )
      call mcik01 ( )
      call mciklv ( )
      call mcikna ( )
      call mciknb ( )
      call mcikva ( )
      call mcikvb ( )
      call mcisia ( )
      call mcisib ( )
      call mcjk ( )
      call mcjy01 ( )
      call mcjylv ( )
      call mcjyna ( )
      call mcjynb ( )
      call mcjyva ( )
      call mcjyvb ( )
      call mclpmn ( )
      call mclpn ( )
      call mclqmn ( )
      call mclqn ( )
      call mcomelp ( )
      call mcpbdn ( )
      call mcpsi ( )
      call mcsphik ( )
      call mcsphjy ( )
      call mcva1 ( )
      call mcva2 ( )
      call mcyzo ( )

      call me1xa ( )
      call me1xb ( )
      call me1z ( )
      call meix ( )
      call melit ( )
      call melit3 ( )
      call menxa ( )
      call menxb ( )
      call merror ( )
      call meulera ( )
      call meulerb ( )

      call mfcoef ( )
      call mfcs ( )
      call mfcszo ( )
      call mffk ( )

      call mgamma ( )

      call mherzo ( )
      call mhygfx ( )
      call mhygfz ( )

      call mik01a ( )
      call mik01b ( )
      call mikna ( )
      call miknb ( )
      call mikv ( )
      call mincob ( )
      call mincog ( )
      call mitairy ( )
      call mitika ( )
      call mitikb ( )
      call mitjya ( )
      call mitjyb ( )
      call mitsh0 ( )
      call mitsl0 ( )
      call mitth0 ( )
      call mittika ( )
      call mittikb ( )
      call mittjya ( )
      call mittjyb ( )

      call mjdzo ( )
      call mjelp ( )
      call mjy01a ( )
      call mjy01b ( )
      call mjyna ( )
      call mjynb ( )
      call mjyv ( )
      call mjyzo ( )

      call mklvna ( )
      call mklvnb ( )
      call mklvnzo ( )

      call mlagzo ( )
      call mlegzo ( )
      call mlgama ( )
      call mlamn ( )
      call mlamv ( )
      call mlpmn ( )
      call mlpmns ( )
      call mlpmv ( )
      call mlpn ( )
      call mlpni ( )
      call mlqmn ( )
      call mlqmns ( )
      call mlqna ( )
      call mlqnb ( )

      call mmtu0 ( )
      call mmtu12 ( )

      call mothpl ( )

      call mpbdv ( )
      call mpbvv ( )
      call mpbwa ( )
      call mpsi ( )

      call mrctj ( )
      call mrcty ( )
      call mrswfo ( )
      call mrswfp ( )

      call mscka ( )
      call msckb ( )
      call msegv ( )
      call msdmn ( )
      call msphi ( )
      call msphj ( )
      call msphk ( )
      call msphy ( )
      call mstvh0 ( )
      call mstvh1 ( )
      call mstvhv ( )
      call mstvl0 ( )
      call mstvl1 ( )
      call mstvlv ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPECIAL_FUNCTIONS_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine mairya ( )

c*********************************************************************72
c
cc MAIRYA tests AIRYA.
c
c  Discussion:
c
c    This program computes Airy functions and their 
c    derivatives using subroutine AIRYA.
c
c  Modified:
c
c    07 February 2012
c
c       Input:   x  --- Argument of Airy function
c       Output:  AI --- Ai(x)
c                BI --- Bi(x)
c                AD --- Ai'(x)
c                BD --- Bi'(x)
c       Example:
c
c   x       Ai(x)          Bi(x)          Ai'(x)         Bi'(x)
c  ----------------------------------------------------------------
c   0   .35502805D+00  .61492663D+00 -.25881940D+00  .44828836D+00
c  10   .11047533D-09  .45564115D+09 -.35206337D-09  .14292361D+10
c  20   .16916729D-26  .21037650D+26 -.75863916D-26  .93818393D+26
c  30   .32082176D-48  .90572885D+47 -.17598766D-47  .49533045D+48
c
c   x       Ai(-x)         Bi(-x)         Ai'(-x)        Bi'(-x)
c  ----------------------------------------------------------------
c   0       .35502805      .61492663     -.25881940      .44828836
c  10       .04024124     -.31467983      .99626504      .11941411
c  20      -.17640613     -.20013931      .89286286     -.79142903
c  30      -.08796819     -.22444694     1.22862060     -.48369473
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      integer test_num
      parameter ( test_num = 4 )

      double precision x_test(test_num)

      save x_test

      data x_test /
     &  0.0D+00, 10.0D+00, 20.0D+00, 30.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MAIRYA'
      write ( *, '(a)' ) '  Test AIRYA'

      write(*,30)
      write(*,40)

      do i = 1, test_num
        x = x_test(i)
        call airya(x,ai,bi,ad,bd)
        write(*,10)x,ai,bi,ad,bd
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) ' '
      write(*,50)
      write(*,40)

      do i = 1, test_num
        x = x_test(i)
        call airya(-x,ai,bi,ad,bd)
        write(*,20)x,ai,bi,ad,bd
      end do

      return
10      FORMAT(1X,F5.1,4D16.8)
20      FORMAT(1X,F5.1,4D16.8)
30      FORMAT(4X,'x',8X,'Ai(x)',11X,'Bi(x)',11X,'Ai''(x)',
     &         10X,'Bi''(x)')
40      FORMAT(2X,'----------------------------------',
     &        '-----------------------------------')
50      FORMAT(4X,'x',8X,'Ai(-x)',10X,'Bi(-x)',10X,
     &        'Ai''(-x)',9X,'Bi''(-x)')
        END
      subroutine mairyb ( )

c*********************************************************************72
c
cc MAIRYB tests AIRYB.
c
c  Discussion:
c
c    This program computes Airy functions and their 
c    derivatives using subroutine AIRYB
c
c  Modified:
c
c    07 February 2012
c
c       Input:   x  --- Argument of Airy function
c       Output:  AI --- Ai(x)
c                BI --- Bi(x)
c                AD --- Ai'(x)
c                BD --- Bi'(x)
c       Example:
c
c   x       Ai(x)          Bi(x)          Ai'(x)         Bi'(x)
c  ----------------------------------------------------------------
c   0   .35502805D+00  .61492663D+00 -.25881940D+00  .44828836D+00
c  10   .11047533D-09  .45564115D+09 -.35206337D-09  .14292361D+10
c  20   .16916729D-26  .21037650D+26 -.75863916D-26  .93818393D+26
c  30   .32082176D-48  .90572885D+47 -.17598766D-47  .49533045D+48
c
c   x       Ai(-x)         Bi(-x)         Ai'(-x)        Bi'(-x)
c  ----------------------------------------------------------------
c   0       .35502805      .61492663     -.25881940      .44828836
c  10       .04024124     -.31467983      .99626504      .11941411
c  20      -.17640613     -.20013931      .89286286     -.79142903
c  30      -.08796819     -.22444694     1.22862060     -.48369473
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      integer test_num
      parameter ( test_num = 4 )

      double precision x_test(test_num)

      save x_test

      data x_test /
     &  0.0D+00, 10.0D+00, 20.0D+00, 30.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MAIRYB'
      write ( *, '(a)' ) '  Test AIRYB'

      WRITE(*,30)
      WRITE(*,40)

      do i = 1, test_num
        x = x_test(i)
        CALL AIRYB(X,AI,BI,AD,BD)
        WRITE(*,10)X,AI,BI,AD,BD
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) ' '
      WRITE(*,50)
      WRITE(*,40)

      do i = 1, test_num
        x = x_test(i)
        CALL AIRYB(-X,AI,BI,AD,BD)
        WRITE(*,20)X,AI,BI,AD,BD
      end do

      return
10      FORMAT(1X,F5.1,4D16.8)
20      FORMAT(1X,F5.1,4D16.8)
30      FORMAT(4X,'x',8X,'Ai(x)',11X,'Bi(x)',11X,'Ai''(x)',
     &         10X,'Bi''(x)')
40      FORMAT(2X,'----------------------------------',
     &        '-----------------------------------')
50      FORMAT(4X,'x',8X,'Ai(-x)',10X,'Bi(-x)',10X,
     &        'Ai''(-x)',9X,'Bi''(-x)')
        END
      subroutine MAIRYZO

c*********************************************************************72
c
cc MAIRYZO tests AIRYZO.
c
c  Discussion:
c
c    This program computes the first NT zeros of Airy functions Ai(x) and Ai'(x), 
c    and the associated values of Ai(a') and Ai'(a), and the first NT zeros of 
c    Airy functions Bi(x) and Bi'(x), and the associated values of Bi(b') and 
c    Bi'(b) using subroutine AIRYZO.
c
c  Modified:
c
c    07 February 2012
c
c       Input :  NT    --- Total number of zeros
c                KF    --- Function code
c                          KF=1 for Ai(x) and Ai'(x)
c                          KF=2 for Bi(x) and Bi'(x)
c       Output:  XA(m) --- a, the m-th zero of Ai(x) or
c                          b, the m-th zero of Bi(x) 
c                XB(m) --- a', the m-th zero of Ai'(x) or
c                          b', the m-th zero of Bi'(x)
c                XC(m) --- Ai(a') or Bi(b')
c                XD(m) --- Ai'(a) or Bi'(b)
c                          ( m --- Serial number of zeros )
c       Example: NT=5
c
c       m         a            Ai'(a)         a'          Ai(a')
c      -----------------------------------------------------------
c       1    -2.33810741     .70121082   -1.01879297    .53565666
c       2    -4.08794944    -.80311137   -3.24819758   -.41901548
c       3    -5.52055983     .86520403   -4.82009921    .38040647
c       4    -6.78670809    -.91085074   -6.16330736   -.35790794
c       5    -7.94413359     .94733571   -7.37217726    .34230124
c
c       m         b            Bi'(b)         b'          Bi(b')
c      -----------------------------------------------------------
c       1    -1.17371322     .60195789   -2.29443968   -.45494438
c       2    -3.27109330    -.76031014   -4.07315509    .39652284
c       3    -4.83073784     .83699101   -5.51239573   -.36796916
c       4    -6.16985213    -.88947990   -6.78129445    .34949912
c       5    -7.37676208     .92998364   -7.94017869   -.33602624
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION XA(50),XB(50),XC(50),XD(50)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MAIRYZO'
      write ( *, '(a)' ) '  Test AIRYZO'

        kf = 1
        nt = 5
      write ( *, '(a)' ) ' '
        WRITE(*,30)KF,NT
           WRITE(*,*)'  m        a             Ai''(a)        a''',
     &               '           Ai(a'')'
        WRITE(*,*)'---------------------------------',
     & '---------------------------'
        CALL AIRYZO(NT,KF,XA,XB,XC,XD)         
      DO K=1,NT
        WRITE(*,20)K,XA(K),XD(K),XB(K),XC(K)
      end do
        write ( *, '(a)' ) ' '
        kf = 2
        nt = 5
        WRITE(*,30)KF,NT
           WRITE(*,*)'  m        b             Bi''(b)        b''',
     &               '           Bi(b'')'
        WRITE(*,*)'---------------------------------',
     & '---------------------------'
        CALL AIRYZO(NT,KF,XA,XB,XC,XD)         
      DO K=1,NT
        WRITE(*,20)K,XA(K),XD(K),XB(K),XC(K)
      end do

      return
20      FORMAT(1X,I3,1X,3F14.8,F13.8)
30      FORMAT(2X,3HKF=,I2,',     ',3HNT=,I3)
35      FORMAT(10X,'KF=1 for Ai(x) and Ai''(x); KF=2 for Bi(x)',
     &          ' and Bi''(x)')
40      FORMAT(10X,'NT is the number of the zeros')
      END
      subroutine MASWFA ( )

c*********************************************************************72
c
cc MASWFA tests ASWFA.
c
c  Modified:
c
c    07 February 2012
c
c
c       Purpose: This program computes the prolate and oblate 
c                spheroidal angular functions of the first
c                kind and their derivatives using subroutine ASWFA
c       Input :  m  --- Mode parameter,  m = 0,1,2,...
c                n  --- Mode parameter,  n = m,m+1,...
c                c  --- Spheroidal parameter
c                x  --- Argument of angular function, |x| < 1.0
c                KD --- Function code
c                       KD=1 for prolate;  KD=-1 for oblate
c                cv --- Characteristic value
c       Output:  S1F --- Angular function of the first kind
c                S1D --- Derivative of the angular function of
c                        the first kind
c       Examples:
c               KD = 1, m = 2, n = 3, c = 3.0 and cv = 14.8277782138
c                  x         Smn(c,x)            Smn'(c,x)
c                --------------------------------------------
c                 0.2      .28261309D+01       .12418631D+02
c                 0.5      .49938554D+01       .92761604D+00
c                 0.8      .31693975D+01      -.12646552D+02
c
c               KD =-1, m = 2, n = 3, c = 3.0 and cv = 8.8093939208
c                  x         Smn(-ic,x)         Smn'(-ic,x)
c                --------------------------------------------
c                 0.2      .29417848D+01       .14106305D+02
c                 0.5      .64138827D+01       .76007194D+01
c                 0.8      .60069873D+01      -.14387479D+02
c
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      integer test_num
      parameter ( test_num = 2 )

      double precision c_test(test_num)
      DIMENSION EG(200)
      integer i
      integer j
      integer kd_test(test_num)
      integer m_test(test_num)
      integer n_test(test_num)

      save c_test
      save kd_test
      save m_test
      save n_test

      data c_test / 3.0D+00, 3.0D+00 /
      data kd_test / 1, -1 /
      data m_test / 2, 2 /
      data n_test / 3, 3 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MASWFA'
      write ( *, '(a)' ) '  Test ASWFA'

      do j = 1, test_num

        kd = kd_test(j)
        m = m_test(j)
        n = n_test(j)
        c = c_test(j)

        write ( *, '(a)' ) ' '
        WRITE(*,10)KD,M,N,C
        CALL SEGV(M,N,C,KD,CV,EG)
        WRITE(*,20)CV
        WRITE(*,*)
        IF (KD.EQ.1 ) THEN
          WRITE(*,*)'    x         Smn(c,x)            Smn''(c,x)'
        ELSE IF (KD.EQ.-1) THEN
          WRITE(*,*)'    x         Smn(-ic,x)         Smn''(-ic,x)'
        ENDIF
        WRITE(*,*)'  --------------------------------------------'
        DO I=0,20
          X=-1.0D0+0.1D0*I
          CALL ASWFA(M,N,C,X,KD,CV,S1F,S1D)  
          WRITE(*,30)X,S1F,S1D
        end do

      end do

      return
10    FORMAT(1X,'KD ='I2,', ','m =',I2,', ','n =',I2,', ','c =',F5.1)
20    FORMAT(1X,' cv =',F18.10)
30    FORMAT(1X,F5.1,2D20.8)
      END
      subroutine maswfb ( )

c*********************************************************************72
c
cc MASWFB tests ASWFB.
c
c  Modified:
c
c    12 March 2012
c
c
c       Purpose: This program computes the prolate and oblate 
c                spheroidal angular functions of the first kind
c                and their derivatives using subroutine ASWFB
c       Input :  m  --- Mode parameter,  m = 0,1,2,...
c                n  --- Mode parameter,  n = m,m+1,...
c                c  --- Spheroidal parameter
c                x  --- Argument of angular function, |x| ó 1.0
c                KD --- Function code
c                       KD=1 for prolate;  KD=-1 for oblate
c                cv --- Characteristic value
c       Output:  S1F --- Angular function of the first kind
c                S1D --- Derivative of the angular function of
c                        the first kind
c       Examples:
c               KD = 1, m = 2, n = 3, c = 3.0 and cv = 14.8277782138
c                  x         Smn(c,x)            Smn'(c,x)
c                --------------------------------------------
c                 0.2      .28261309D+01       .12418631D+02
c                 0.5      .49938554D+01       .92761604D+00
c                 0.8      .31693975D+01      -.12646552D+02
c
c               KD =-1, m = 2, n = 3, c = 3.0 and cv = 8.8093939208
c                  x         Smn(-ic,x)         Smn'(-ic,x)
c                --------------------------------------------
c                 0.2      .29417848D+01       .14106305D+02
c                 0.5      .64138827D+01       .76007194D+01
c                 0.8      .60069873D+01      -.14387479D+02
c
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      integer test_num
      parameter ( test_num = 2 )

      double precision c
      double precision c_test(test_num)
      double precision EG(200)
      integer i
      integer j
      integer kd_test(test_num)
      integer m_test(test_num)
      integer n_test(test_num)

      save c_test
      save kd_test
      save m_test
      save n_test

      data c_test / 3.0D+00, 3.0D+00 /
      data kd_test / 1, -1 /
      data m_test / 2, 2 /
      data n_test / 3, 3 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MASWFB'
      write ( *, '(a)' ) '  Test ASWFB'

      do j = 1, test_num

        kd = kd_test(j)
        m = m_test(j)
        n = n_test(j)
        c = c_test(j)

        write ( *, '(a)' ) ' '
        WRITE(*,10)KD,M,N,C
        CALL SEGV(M,N,C,KD,CV,EG)
        WRITE(*,20)CV
        WRITE(*,*)

        IF (KD.EQ.1 ) THEN
          WRITE(*,*)'    x         Smn(c,x)            Smn''(c,x)'
        ELSE IF (KD.EQ.-1) THEN
          WRITE(*,*)'    x         Smn(-ic,x)         Smn''(-ic,x)'
        ENDIF
        WRITE(*,*)'  --------------------------------------------'
        DO I=0,20
          X=-1.0D0+0.1D0*I
          CALL ASWFB(M,N,C,X,KD,CV,S1F,S1D)  
          WRITE(*,30)X,S1F,S1D
        end do

      end do

10      FORMAT(1X,'KD ='I2,', ','m =',I2,', ','n =',I2,', ','c =',F5.1)
20      FORMAT(1X,' cv =',F18.10)
30      FORMAT(1X,F5.1,2D20.8)
      return
      END
      subroutine mbernoa ( )

c*********************************************************************72
c
cc MBERNOA tests BERNOA.
c
c  Modified:
c
c    12 March 2012
c
c
c       Purpose: This program computes Bernoulli number Bn using 
c              subroutine BERNOA
c       Example: Compute Bernouli number Bn for n = 0,1,...,10
c                Computed results:
c
c                   n            Bn
c                 --------------------------
c                   0     .100000000000D+01
c                   1    -.500000000000D+00
c                   2     .166666666667D+00
c                   4    -.333333333333D-01
c                   6     .238095238095D-01
c                   8    -.333333333333D-01
c                  10     .757575757576D-01
c
      integer n
      parameter ( n = 10 )

      double precision b(0:n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MBERNOA:'
      write ( *, '(a)' ) '  Test BERNOA.'

        CALL BERNOA(N,B)
        WRITE(*,*)'   n            Bn'
        WRITE(*,*)' --------------------------'
        WRITE(*,20)0,B(0)
        WRITE(*,20)1,B(1)
        DO K=2,N,2
           WRITE(*,20)K,B(K)
      end do
20      FORMAT(2X,I3,D22.12)

      return
      END
      subroutine mbernob ( )

c*********************************************************************72
c
cc MBERNOB tests BERNOB.
c
c  Modified:
c
c    12 March 2012
c
c
c       Purpose: This program computes Bernoulli number Bn using 
c              subroutine BERNOB
c       Example: Compute Bernouli number Bn for n = 0,1,...,10
c                Computed results:
c
c                   n            Bn
c                 --------------------------
c                   0     .100000000000D+01
c                   1    -.500000000000D+00
c                   2     .166666666667D+00
c                   4    -.333333333333D-01
c                   6     .238095238095D-01
c                   8    -.333333333333D-01
c                  10     .757575757576D-01
c
c
      integer n
      parameter ( n = 10 )

      double precision b(0:n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MBERNOB:'
      write ( *, '(a)' ) '  Test BERNOB.'

        CALL BERNOB(N,B)
        WRITE(*,*)'   n            Bn'
        WRITE(*,*)' --------------------------'
        WRITE(*,20)0,B(0)
        WRITE(*,20)1,B(1)
        DO 10 K=2,N,2
10         WRITE(*,20)K,B(K)
20      FORMAT(2X,I3,D22.12)

      return
      END
      subroutine mbeta ( )

c*********************************************************************72
c
cc MBETA tests BETA.
c
c  Modified:
c
c    12 March 2012
c
c
c       Purpose: This program computes the beta function 
c                B(p,q) for p > 0 and q > 0 using 
c              subroutine BETA
c       Input :  p  --- Parameter  ( p > 0 )
c                q  --- Parameter  ( q > 0 )
c       Output:  BT --- B(p,q)
c       Examples:
c                 p       q           B(p,q)
c               ---------------------------------
c                1.5     2.0     .2666666667D+00
c                2.5     2.0     .1142857143D+00
c                1.5     3.0     .1523809524D+00
c
c
      implicit none

      integer test_num
      parameter ( test_num = 3 )

      double precision bt
      double precision p
      double precision p_test(test_num)
      double precision q
      double precision q_test(test_num)
      integer test

      save p_test
      save q_test

      data p_test / 1.5D+00, 2.5D+00, 1.5D+00 /
      data q_test / 2.0D+00, 2.0D+00, 3.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MBETA:'
      write ( *, '(a)' ) '  Test BETA.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    p       q           B(p,q)'
      write ( *, '(a)' ) '  ---------------------------------'

      do test = 1, test_num

        p = p_test(test)
        q = q_test(test)

        call beta ( p, q, bt )
        write ( *, '(2x,f5.1,3x,f5.1,d20.10)' ) p, q, bt

      end do

      return
      end
      subroutine mcchg ( )

c*********************************************************************72
c
cc MCCHG tests CCHG.
c
c  Modified:
c
c    13 March 2012
c
c
c       Purpose: This program computes confluent hypergeometric 
c                function M(a,b,z) with real parameters a, b, and 
c                a complex argument z using subroutine CCHG
c       Input :  a --- Parameter
c                b --- Parameter
c                z --- Complex argument
c       Output:  CHG --- M(a,b,z)
c       Examples:
c          a      b        z        Re[M(a,b,z)]   Im[M(a,b,z)]
c         -------------------------------------------------------
c         3.3   4.25    10 + 0i    .61677489D+04    0
c         3.3   4.25    25 + 0i    .95781835D+10  -.15738228D-03
c         3.3   4.25     3 -  i    .75828716D+01  -.86815474D+01
c         3.3   4.25    15 +10i   -.58313765D+06  -.48195426D+05
c
c
        IMPLICIT DOUBLE PRECISION (A,B,X,Y)
        IMPLICIT COMPLEX *16 (C,Z)

      integer test_num
      parameter ( test_num = 4 )

      double precision a
      double precision a_test(test_num)
      double precision b
      double precision b_test(test_num)
      integer test
      double precision x
      double precision x_test(test_num)
      double precision y
      double precision y_test(test_num)

      save a_test
      save b_test
      save x_test
      save y_test

      data a_test / 3.3D+00, 3.3D+00, 3.3D+00, 3.3D+00 /
      data b_test / 4.25D+00, 4.25D+00, 4.25D+00, 4.25D+00 /
      data x_test / 10.0D+00, 25.0D+00, 3.0D+00, 15.0D+00 /
      data y_test /  0.0D+00,  0.0D+00, -1.0D+00, 10.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MCCHG'
      write ( *, '(a)' ) '  Test CCHG'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     A     B     re(Z)     im(Z)     M(a,b,z)'
      write ( *, '(a)' ) ' '

      do test = 1, test_num

        a = a_test(test)
        b = b_test(test)
        x = x_test(test)
        y = y_test(test)

        Z=CMPLX(X,Y)
        CALL CCHG(A,B,Z,CHG)

        write ( *, 
     &    '(2x,f5.1,2x,f5.2,2x,f5.1,2x,f5.1,2x,d18.8,2x,d18.8)' ) 
     &    a, b, x, y, chg

      end do

      return
      END
      subroutine mcerror ( )

c*********************************************************************72
c
cc MCERROR tests CERROR.
c
c  Modified:
c
c    13 March 2012
c
c
c       Purpose: This program computes the error function erf(z) 
c                for a complex argument using subroutine CERROR
c       Input :  x   --- Real part of z
c                y   --- Imaginary part of z  ( y ó 3.0 )
c       Output:  ERR --- Real part of erf(z)
c                ERI --- Imaginary part of erf(z)
c       Example:
c                   x       y       Re[erf(z)]      Im[erf(z)]
c                 ---------------------------------------------
c                  1.0     2.0      -.53664357     -5.04914370
c                  2.0     2.0      1.15131087       .12729163
c                  3.0     2.0       .99896328      -.00001155
c                  4.0     2.0      1.00000057      -.00000051
c                  5.0     2.0      1.00000000       .00000000
c
c
      IMPLICIT COMPLEX *16 (C,Z)  

      integer test_num
      parameter ( test_num = 5 )

      integer test
      double precision x
      double precision x_test(test_num)
      double precision y
      double precision y_test(test_num)

      save x_test
      save y_test

      data x_test / 1.0D+00, 2.0D+00, 3.0D+00, 4.0D+00, 5.0D+00 /
      data y_test / 2.0D+00,  2.0D+00, 2.0D+00, 2.0D+00, 2.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MCERROR'
      write ( *, '(a)' ) '  Test CERROR.'
      write ( *, '(a)' ) ' '
      WRITE(*,*)'   x      y      Re[erf(z)]      Im[erf(z)]'
      WRITE(*,*)' ---------------------------------------------'
      do test = 1, test_num
        x = x_test(test)
        y = y_test(test)
        Z=CMPLX(X,Y)
        CALL CERROR(Z,CER)
      WRITE(*,10) Z,CER
10      FORMAT(1X,F5.1,2X,F5.1,1X,2E16.8)
      end do

      return
      end
      subroutine mcerzo ( )

c*********************************************************************72
c
cc MCERZO tests CERZO
c
c  Modified:
c
c    13 March 2012
c
c
c       Purpose : This program evaluates the complex zeros of error 
c                 function erf(z) using subroutine CERZO
c       Input:    NT --- Total number of zeros
c       Example:  NT = 10
c
c    n     complex zeros of erf(z)     n     complex zeros of erf(z)
c   -------------------------------------------------------------------
c    1   1.450616163 + i 1.880943000   6   4.158998400 + i 4.435571444
c    2   2.244659274 + i 2.616575141   7   4.516319400 + i 4.780447644
c    3   2.839741047 + i 3.175628100   8   4.847970309 + i 5.101588043
c    4   3.335460735 + i 3.646174376   9   5.158767908 + i 5.403332643
c    5   3.769005567 + i 4.060697234  10   5.452192201 + i 5.688837437
c
c
        IMPLICIT DOUBLE PRECISION (E,P,W)
        IMPLICIT COMPLEX *16 (C,Z)

      integer nt
      parameter ( nt = 10 )
      complex*16 ZO(nt)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MCERZO'
      write ( *, '(a)' ) '  Test CERZO.'
      write ( *, '(a)' ) ' '
        WRITE(*,20)NT
        CALL CERZO(NT,ZO)
        WRITE(*,*)
        WRITE(*,*)'  n        Complex zeros of erf(z)'
        WRITE(*,*)'-------------------------------------'
        DO 10 I=1,NT
10         WRITE(*,30) I,ZO(I)
20      FORMAT(2X,'NT=',I3)
30      FORMAT(1X,I3,2X,F13.8,2X,2H+i,F13.8)
      return
        END
      subroutine mcgama ( )

c*********************************************************************72
c
cc MCGAMA tests CGAMA.
c
c  Modified:
c
c    13 March 2012
c
c
c       Purpose: This program computes the gamma function â(z)  
c                or ln[â(z)] for a complex argument using 
c              subroutine CGAMA
c       Input :  x  --- Real part of z
c                y  --- Imaginary part of z
c                KF --- Function code
c                       KF=0 for ln[â(z)]
c                       KF=1 for â(z)
c       Output:  GR --- Real part of ln[â(z)] or â(z)
c                GI --- Imaginary part of ln[â(z)] or â(z)
c       Examples:
c
c         x         y           Re[â(z)]           Im[â(z)]
c       --------------------------------------------------------
c        2.50      5.00     .2267360319D-01    -.1172284404D-01
c        5.00     10.00     .1327696517D-01     .3639011746D-02
c        2.50     -5.00     .2267360319D-01     .1172284404D-01
c        5.00    -10.00     .1327696517D-01    -.3639011746D-02
c
c         x         y          Re[lnâ(z)]         Im[lnâ(z)]
c      ---------------------------------------------------------
c        2.50      5.00    -.3668103262D+01     .5806009801D+01
c        5.00     10.00    -.4285507444D+01     .1911707090D+02
c        2.50     -5.00    -.3668103262D+01    -.5806009801D+01
c        5.00    -10.00    -.4285507444D+01    -.1911707090D+02
c
c
        DOUBLE PRECISION GR,GI

      integer test_num
      parameter ( test_num = 4 )

      integer kf
      integer test
      double precision x
      double precision x_test(test_num)
      double precision y
      double precision y_test(test_num)

      save x_test
      save y_test

      data x_test / 2.5D+00,  5.0D+00,  2.5D+00,   5.0D+00 /
      data y_test / 5.0D+00, 10.0D+00, -5.0D+00, -10.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MCGAMA'
      write ( *, '(a)' ) '  Test CGAMA'

      do kf = 1, 2

        write ( *, '(a)' ) ' '
        IF (KF.EQ.1) THEN
            WRITE(*,*)'       x         y           Re[â(z)]',
     &                '           Im[â(z)]'
        ELSE
            WRITE(*,*)'       x         y          Re[lnâ(z)]',
     &                '         Im[lnâ(z)]'
        ENDIF
        WRITE(*,*)'    ------------------------------------',
     &            '---------------------'
        do test = 1, test_num
          x = x_test(test)
          y = y_test(test)
          CALL CGAMA(X,Y,KF,GR,GI)
          WRITE(*,10)X,Y,GR,GI
        end do

10      FORMAT(1X,2F10.2,2D20.10)

      end do

      return
      end
      subroutine mch12n ( )

c*********************************************************************72
c
cc MCH12N tests CH12N.
c
c  Modified:
c
c    13 March 2012
c
c
c       Purpose: This program computes Hankel functions of
c                the first and second kinds and their 
c                derivatives for a complex argument using
c              subroutine CH12N
c       Input :  z --- Complex argument
c                n --- Order of Hn(1)(z) and Hn(2)(z)
c                      ( n = 0,1,úúú, n ó 250 )
c       Output:  CHF1(n) --- Hn(1)(z)
c                CHD1(n) --- Hn(1)'(z)
c                CHF2(n) --- Hn(2)(z)
c                CHD2(n) --- Hn(2)'(z)
c
c
        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)

      integer test_num
      parameter ( test_num = 4 )

      DIMENSION CHF1(0:250),CHD1(0:250),CHF2(0:250),CHD2(0:250)
      integer n
      integer n_test(test_num)
      integer ns
      integer ns_test(test_num)
      integer test
      double precision x
      double precision x_test(test_num)
      double precision y
      double precision y_test(test_num)

      save n_test
      save x_test
      save y_test

      data n_test / 1, 2, 3, 4 /
      data ns_test / 1, 1, 1, 1 /
      data x_test / 2.5D+00,  5.0D+00,  2.5D+00,   5.0D+00 /
      data y_test / 5.0D+00, 10.0D+00, -5.0D+00, -10.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MCH12N'
      write ( *, '(a)' ) '  Test CH12N'

      do test = 1, test_num

        Z=CMPLX(X,Y)

        IF (N.LE.8) THEN
           NS=1
        ELSE
           ns = ns_test(test)
        ENDIF

        CALL CH12N(N,Z,NM,CHF1,CHD1,CHF2,CHD2)

        WRITE(*,*)
        WRITE(*,*)'   n     Re[Hn(1)(z)]     Im[Hn(1)(z)]',
     &            '      Re[Hn(1)''(z)]     Im[Hn(1)''(z)]'
        WRITE(*,*)' -------------------------------------',
     &               '---------------------------------------'
        DO K=0,NM,NS
           WRITE(*,40)K,CHF1(K),CHD1(K)
        end do

        WRITE(*,*)
        WRITE(*,*)'   n     Re[Hn(2)(z)]     Im[Hn(2)(z)]',
     &            '      Re[Hn(2)''(z)]     Im[Hn(2)''(z)]'
        WRITE(*,*)' -------------------------------------',
     &            '---------------------------------------'
        DO K=0,NM,NS
          WRITE(*,40)K,CHF2(K),CHD2(K)
        end do

      end do
40      FORMAT(1X,I4,4D18.10)
45      FORMAT(3X,3Hz =,F8.3,' + i ',F8.3,' ,',6X,6HNmax =,I4)

      return
        END
      subroutine mchgm ( )

c*********************************************************************72
c
cc MCHGM tests CHGM.
c
c  Modified:
c
c    15 March 2012
c
c
c       Purpose: This program computes the confluent
c                hypergeometric function M(a,b,x) using
c              subroutine CHGM
c       Input  : a  --- Parameter
c                b  --- Parameter ( b <> 0,-1,-2,... )
c                x  --- Argument
c       Output:  HG --- M(a,b,x)
c       Example:
c                   a       b       x          M(a,b,x)
c                 -----------------------------------------
c                  1.5     2.0    20.0     .1208527185D+09
c                  4.5     2.0    20.0     .1103561117D+12
c                 -1.5     2.0    20.0     .1004836854D+05
c                 -4.5     2.0    20.0    -.3936045244D+03
c                  1.5     2.0    50.0     .8231906643D+21
c                  4.5     2.0    50.0     .9310512715D+25
c                 -1.5     2.0    50.0     .2998660728D+16
c                 -4.5     2.0    50.0    -.1806547113D+13
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      double precision a
      double precision a_test(4)
      double precision b
      double precision b_test(1)
      integer i
      integer j
      integer k
      double precision x
      double precision x_test(2)

      save a_test
      save b_test
      save x_test

      data a_test / 1.5D+00, 4.5D+00, -1.5D+00, -4.5D+00 /
      data b_test / 2.0D+00 /
      data x_test / 20.0D+00, 50.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MCHGM'
      write ( *, '(a)' ) '  Test CHGM'
      write ( *, '(a)' ) ' '
        WRITE(*,*)'   a       b       x          M(a,b,x)'
        WRITE(*,*)' -----------------------------------------'

      do k = 1, 2
        x = x_test(k)
        do j = 1, 1
          b = b_test(j)
          do i = 1, 4
            a = a_test(i)

        CALL CHGM(A,B,X,HG)
        WRITE(*,10)A,B,X,HG
10      FORMAT(1X,F5.1,3X,F5.1,3X,F5.1,D20.10)

          end do
        end do
      end do

      return
        END
      subroutine mchgu ( )

c*********************************************************************72
c
cc MCHGU tests CHGU.
c
c  Modified:
c
c    15 March 2012
c
c
c       Purpose: This program computes the confluent
c                hypergeometric function U(a,b,x) using
c              subroutine CHGU
c       Input  : a  --- Parameter
c                b  --- Parameter
c                x  --- Argument  ( x ò 0 )
c       Output:  HU --- U(a,b,x)
c                MD --- Method code
c       Example:
c                a       b       x        U(a,b,x)
c             --------------------------------------
c              -2.5     2.5     5.0     -9.02812446
c              -1.5     2.5     5.0      2.15780560
c               -.5     2.5     5.0      1.76649370
c                .0     2.5     5.0      1.00000000
c                .5     2.5     5.0       .49193496
c               1.5     2.5     5.0       .08944272
c               2.5     2.5     5.0       .01239387
c
c                a       b       x        U(a,b,x)
c             --------------------------------------
c              -2.5     5.0    10.0     -2.31982196
c              -1.5     5.0    10.0      8.65747115
c               -.5     5.0    10.0      2.37997143
c                .0     5.0    10.0      1.00000000
c                .5     5.0    10.0       .38329536
c               1.5     5.0    10.0       .04582817
c               2.5     5.0    10.0       .00444535
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      double precision a
      double precision a_test(7)
      double precision b
      double precision b_test(2)
      integer i
      integer j
      double precision x
      double precision x_test(2)

      save a_test
      save b_test
      save x_test

      data a_test / -2.5D+00, -1.5D+00, -0.5D+00, 0.0D+00, 0.5D+00, 
     &  1.5D+00, 2.5D+00 /
      data b_test / 2.5D+00, 5.0D+00 /
      data x_test / 5.0D+00, 10.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MCHGU'
      write ( *, '(a)' ) '  Test CHGU'
      write ( *, '(a)' ) ' '
        WRITE(*,*)'   a       b       x        U(a,b,x)'
        WRITE(*,*)'--------------------------------------'

      do j = 1, 2
        x = x_test(j)
        b = b_test(j)
        do i = 1, 7
          a = a_test(i)
          CALL CHGU(A,B,X,HU,MD)
          WRITE(*,10)A,B,X,HU
        end do
      end do

10      FORMAT(1X,F5.1,3X,F5.1,3X,F5.1,5x,E15.8)

      return
        END
      subroutine mcik01 ( )

c*********************************************************************72
c
cc MCIK01 tests CIK01.
c
c  Modified:
c
c    16 March 2012
c
c
c       Purpose: This program computes the modified Bessel functions  
c                I0(z), I1(z), K0(z), K1(z), and their derivatives 
c                for a complex argument using subroutine CIK01
c       Input :  z --- Complex argument
c       Output:  CBI0 --- I0(z)
c                CDI0 --- I0'(z)
c                CBI1 --- I1(z)
c                CDI1 --- I1'(z)
c                CBK0 --- K0(z)
c                CDK0 --- K0'(z)
c                CBK1 --- K1(z)
c                CDK1 --- K1'(z)
c       Example: z = 20.0 + i 10.0
c
c     n      Re[In(z)]      Im[In(z)]      Re[In'(z)]     Im[In'(z)]
c    -----------------------------------------------------------------
c     0   -.38773811D+08 -.13750292D+08 -.37852037D+08 -.13869150D+08
c     1   -.37852037D+08 -.13869150D+08 -.36982347D+08 -.13952566D+08
c
c     n      Re[Kn(z)]      Im[Kn(z)]      Re[Kn'(z)]     Im[Kn'(z)]
c    -----------------------------------------------------------------
c     0   -.37692389D-09  .39171613D-09  .38056380D-09 -.40319029D-09
c     1   -.38056380D-09  .40319029D-09  .38408264D-09 -.41545502D-09
c
c
        IMPLICIT DOUBLE PRECISION (X,Y)
        IMPLICIT COMPLEX*16 (C,Z)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MCIK01'
      write ( *, '(a)' ) '  Test CIK01'
      write ( *, '(a)' ) ' '

        x = 20.0D+00
        y = 10.0D+00
        Z=CMPLX(X,Y)
        WRITE(*,30)X,Y
        CALL CIK01(Z,CBI0,CDI0,CBI1,CDI1,CBK0,CDK0,CBK1,CDK1)
        WRITE(*,*)
        WRITE(*,*)'  n      Re[In(z)]      Im[In(z)]',
     &            '      Re[In''(z)]     Im[In''(z)]'
        WRITE(*,*)' -------------------------------',
     &            '----------------------------------'
        WRITE(*,10)CBI0,CDI0
        WRITE(*,20)CBI1,CDI1
        WRITE(*,*)
        WRITE(*,*)'  n      Re[Kn(z)]      Im[Kn(z)]',
     &            '      Re[Kn''(z)]     Im[Kn''(z)]'
        WRITE(*,*)' -------------------------------',
     &            '----------------------------------'
        WRITE(*,10)CBK0,CDK0
        WRITE(*,20)CBK1,CDK1
10      FORMAT(3X,'0',2X,4D15.7)
20      FORMAT(3X,'1',2X,4D15.7)
30      FORMAT(3X,3Hz =,F7.2,' + i',F7.2)

      return
        END
      subroutine mciklv ( )

c*********************************************************************72
c
cc MCIKLV tests CIKLV.
c
c  Modified:
c
c    16 March 2012
c
c
c       Purpose: This program computes modified Bessel functions 
c                Iv(z) and Kv(z) and their derivatives for a  
c                large order and a complex argument using
c              subroutine CIKLV
c       Input:   v --- Order of Iv(z) and Kv(z)
c                z --- Complex argument
c       Output:  CBIV --- Iv(z)
c                CDIV --- Iv'(z)
c                CBKV --- Kv(z)
c                CDKV --- Kv'(z)
c       Examples:
c                v =100.00,    z =   4.00 + i   2.00
c
c       Iv(z) = -.7373606617-123 + .6461109082-123 i
c       Iv'(z)= -.8307094243-122 + .2030132500-121 i
c       Kv(z) = -.3836166007+121 - .3356017795+121 i
c       Kv'(z)=  .1103271276+123 + .2886519240+122 i
c
c                v =100.50,    z =   4.00 + i   2.00
c       Iv(z) = -.1289940051-123 + .6845756182-124 i
c       Iv'(z)= -.1907996261-122 + .2672465997-122 i
c       Kv(z) = -.3008779281+122 - .1593719779+122 i
c       Kv'(z)=  .7653781978+123 + .1857772148+122 i
c
c
        IMPLICIT DOUBLE PRECISION (V,X,Y)
        IMPLICIT COMPLEX*16 (C,Z)

      integer test_num
      parameter ( test_num = 2 )

      integer test
      double precision v
      double precision v_test(test_num)
      double precision x
      double precision x_test(test_num)
      double precision y
      double precision y_test(test_num)

      save v_test
      save x_test
      save y_test

      data v_test / 100.0D+00, 100.5D+00 /
      data x_test / 4.0D+00, 4.0D+00 /
      data y_test / 2.0D+00, 2.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MCIKLV'
      write ( *, '(a)' ) '  Test CIKLV'
      write ( *, '(a)' ) ' '

      do test = 1, test_num

        v = v_test(test)
        x = x_test(test)
        y = y_test(test)

        WRITE(*,10)V,X,Y
        Z=CMPLX(X,Y)
        CALL CIKLV(V,Z,CBIV,CDIV,CBKV,CDKV)
        WRITE(*,*)
        WRITE(*,20)CBIV
        WRITE(*,30)CDIV
        WRITE(*,*)
        WRITE(*,40)CBKV
        WRITE(*,50)CDKV
10      FORMAT(8X,'v =',F6.2,',    ','z =',F7.2,' + i',F7.2)
20      FORMAT(8X,'Iv(z) =',D17.10,' + i ',D17.10)
30      FORMAT(8X,'Iv''(z)=',D17.10,' + i ',D17.10)
40      FORMAT(8X,'Kv(z) =',D17.10,' + i ',D17.10)
50      FORMAT(8X,'Kv''(z)=',D17.10,' + i ',D17.10)

      end do

      return
      END
      subroutine mcikna ( )

c*********************************************************************72
c
cc MCIKNA tests CIKNA.
c
c  Modified:
c
c    17 March 2012
c
c
c       Purpose: This program computes the modified Bessel functions 
c                In(z) and Kn(z), and their derivatives for a  
c                complex argument using subroutine CIKNA
c       Input :  z --- Complex argument of In(z) and Kn(z)
c                n --- Order of In(z) and Kn(z)
c                      ( n = 0,1,úúú, n ó 250 )
c       Output:  CBI(n) --- In(z)
c                CDI(n) --- In'(z)
c                CBK(n) --- Kn(z)
c                CDK(n) --- Kn'(z)
c       Example: z = 4.0 + i 2.0 ,      Nmax = 5
c
c     n     Re[In(z)]      Im[In(z)]      Re[In'(z)]     Im[In'(z)]
c   -----------------------------------------------------------------
c     0  -.19056142D+01  .10403505D+02 -.23059657D+01  .92222463D+01
c     1  -.23059657D+01  .92222463D+01 -.23666457D+01  .83284588D+01
c     2  -.28276772D+01  .62534130D+01 -.24255774D+01  .61553456D+01
c     3  -.25451891D+01  .30884450D+01 -.22270972D+01  .36367893D+01
c     4  -.16265172D+01  .10201656D+01 -.16520416D+01  .16217056D+01
c     5  -.75889410D+00  .15496632D+00 -.94510625D+00  .48575220D+00
c
c     n     Re[Kn(z)]      Im[Kn(z)]      Re[Kn'(z)]     Im[Kn'(z)]
c   -----------------------------------------------------------------
c     0  -.64221754D-02 -.84393648D-02  .74307276D-02  .89585853D-02
c     1  -.74307276D-02 -.89585853D-02  .88041795D-02  .94880091D-02
c     2  -.11186184D-01 -.10536653D-01  .14012532D-01  .10936010D-01
c     3  -.20594336D-01 -.12913435D-01  .27416815D-01  .12106413D-01
c     4  -.43647447D-01 -.13676173D-01  .60982763D-01  .63953943D-02
c     5  -.10137119D+00  .12264588D-03  .14495731D+00 -.37132068D-01
c
c
        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        COMMON CBI(0:250),CDI(0:250),CBK(0:250),CDK(0:250)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MCIKNA'
      write ( *, '(a)' ) '  Test CIKNA'
      write ( *, '(a)' ) ' '

      n = 5
      x = 4.0D+00
      y = 2.0D+00
      ns = 5

        Z=CMPLX(X,Y)
        WRITE(*,40)X,Y,N
        IF (N.LE.8) THEN
           NS=1
        ENDIF
        CALL CIKNA(N,Z,NM,CBI,CDI,CBK,CDK)
        WRITE(*,*)
        WRITE(*,*)'   n      Re[In(z)]       Im[In(z)]',
     &            '       Re[In''(z)]      Im[In''(z)]'
        WRITE(*,*)' -----------------------------------',
     &            '----------------------------------'
        DO K=0,NM,NS
           WRITE(*,30)K,CBI(K),CDI(K)
        end do
        WRITE(*,*)
        WRITE(*,*)'   n      Re[Kn(z)]       Im[Kn(z)]',
     &            '       Re[Kn''(z)]      Im[Kn''(z)]'
        WRITE(*,*)' -----------------------------------',
     &            '----------------------------------'
        DO K=0,NM,NS
           WRITE(*,30)K,CBK(K),CDK(K)
        end do
30      FORMAT(1X,I4,1X,4D16.8)
40      FORMAT(3X,3Hz =,F7.1,' + i',F7.1,' ,',6X,6HNmax =,I4)

        return
        END
      subroutine mciknb ( )

c*********************************************************************72
c
cc MCIKNB tests CIKNB.
c
c  Modified:
c
c    17 March 2012
c
c
c       Purpose: This program computes the modified Bessel functions 
c                In(z) and Kn(z), and their derivatives for a
c                complex argument using subroutine CIKNB
c       Input:   z --- Complex argument
c                n --- Order of In(z) and Kn(z)
c                      ( n = 0,1,úúú, n ó 250 )
c       Output:  CBI(n) --- In(z)
c                CDI(n) --- In'(z)
c                CBK(n) --- Kn(z)
c                CDK(n) --- Kn'(z)
c       Example: Nmax = 5,   z = 4.0 + i 2.0
c
c     n     Re[In(z)]      Im[In(z)]      Re[In'(z)]     Im[In'(z)]
c   -----------------------------------------------------------------
c     0  -.19056142D+01  .10403505D+02 -.23059657D+01  .92222463D+01
c     1  -.23059657D+01  .92222463D+01 -.23666457D+01  .83284588D+01
c     2  -.28276772D+01  .62534130D+01 -.24255774D+01  .61553456D+01
c     3  -.25451891D+01  .30884450D+01 -.22270972D+01  .36367893D+01
c     4  -.16265172D+01  .10201656D+01 -.16520416D+01  .16217056D+01
c     5  -.75889410D+00  .15496632D+00 -.94510625D+00  .48575220D+00
c
c     n     Re[Kn(z)]      Im[Kn(z)]      Re[Kn'(z)]     Im[Kn'(z)]
c   -----------------------------------------------------------------
c     0  -.64221754D-02 -.84393648D-02  .74307276D-02  .89585853D-02
c     1  -.74307276D-02 -.89585853D-02  .88041795D-02  .94880091D-02
c     2  -.11186184D-01 -.10536653D-01  .14012532D-01  .10936010D-01
c     3  -.20594336D-01 -.12913435D-01  .27416815D-01  .12106413D-01
c     4  -.43647447D-01 -.13676173D-01  .60982763D-01  .63953943D-02
c     5  -.10137119D+00  .12264588D-03  .14495731D+00 -.37132068D-01
c
c
        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)

        DIMENSION CBI(0:250),CDI(0:250),CBK(0:250),CDK(0:250)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MCIKNB'
      write ( *, '(a)' ) '  Test CIKNB'
      write ( *, '(a)' ) ' '

      n = 5
      x = 4.0D+00
      y = 2.0D+00
      ns = 5

        WRITE(*,25)N,X,Y
        Z=CMPLX(X,Y)
        WRITE(*,*)
        IF (N.LE.8) THEN
           NS=1
        ENDIF
        CALL CIKNB(N,Z,NM,CBI,CDI,CBK,CDK)
        WRITE(*,*)'   n      Re[In(z)]       Im[In(z)]',
     &            '      Re[In''(z)]       Im[In''(z)]'
        WRITE(*,*)' ---------------------------------',
     &            '------------------------------------'
        DO K=0,NM,NS
           WRITE(*,20)K,CBI(K),CDI(K)
        end do
        WRITE(*,*)
        WRITE(*,*)'   n      Re[Kn(z)]       Im[Kn(z)]',
     &            '      Re[Kn''(z)]       Im[Kn''(z)]'
        WRITE(*,*)' ---------------------------------',
     &            '------------------------------------'
        DO K=0,NM,NS
           WRITE(*,20)K,CBK(K),CDK(K)
        end do
20      FORMAT(1X,1X,I3,1X,4D16.8)
25      FORMAT(3X,6HNmaz =,I3,',    ','z =', F6.1,' + i',F6.1)

      return
        END
      subroutine mcikva ( )

c*********************************************************************72
c
cc MCIKVA tests CIKVA.
c
c  Modified:
c
c    17 March 2012
c
c
c       Purpose: This program computes the modified Bessel functions 
c                Iv(z), Kv(z) and their derivatives for an arbitrary
c                order and complex argument using subroutine CIKVA 
c       Input :  z --- Complex argument z
c                v --- Real order of Iv(z) and Kv(z)
c                      ( v =n+v0,  0 ó n ó 250, 0 ó v0 < 1 )
c       Output:  CBI(n) --- In+v0(z)
c                CDI(n) --- In+v0'(z)
c                CBK(n) --- Kn+v0(z)
c                CDK(n) --- Kn+v0'(z)
c       Example: Compute Iv(z), Kv(z) and their derivatives for
c                v =n+v0, v0=0.25, n =0(1)5, and z =4.0 +i 2.0
c                Computation results:
c
c                v= n+v0,   v0 = .25,   z =  4.0+ i  2.0
c
c      n     Re[Iv(z)]      Im[Iv(z)]     Re[Iv'(z)]     Im[Iv'(z)]
c    -----------------------------------------------------------------
c      0  -.19336550D+01  .10328998D+02 -.23119621D+01  .91612230D+01
c      1  -.24735044D+01  .85964317D+01 -.23898329D+01  .78707023D+01
c      2  -.28460107D+01  .54124063D+01 -.24105909D+01  .55204965D+01
c      3  -.23476775D+01  .24445612D+01 -.21145027D+01  .30604463D+01
c      4  -.13829947D+01  .70848630D+00 -.14732387D+01  .12545751D+01
c      5  -.59879982D+00  .64588999D-01 -.78816416D+00  .32629794D+00
c
c      n     Re[Kv(z)]      Im[Kv(z)]     Re[Kv'(z)]     Im[Kv'(z)]
c     ----------------------------------------------------------------
c      0  -.64820386D-02 -.84715754D-02  .75118612D-02  .89920077D-02
c      1  -.80477525D-02 -.92535355D-02  .96506687D-02  .97789903D-02
c      2  -.12819299D-01 -.11086405D-01  .16310878D-01  .11358076D-01
c      3  -.24574004D-01 -.13462616D-01  .33167751D-01  .11850554D-01
c      4  -.53516204D-01 -.12614703D-01  .75424026D-01  .14407268D-02
c      5  -.12627405D+00  .10581162D-01  .18054884D+00 -.64789392D-01
c
c
        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        COMMON  CBI(0:250),CDI(0:250),CBK(0:250),CDK(0:250)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MCIKVA'
      write ( *, '(a)' ) '  Test CIKVA'
      write ( *, '(a)' ) ' '

      v = 0.25D+00
      x = 4.0D+00
      y = 2.0D+00
      ns = 5

        Z=CMPLX(X,Y)
        N=INT(V)
        V0=V-N
        WRITE(*,25)V0,X,Y
        IF (N.LE.8) THEN
           NS=1
        ENDIF
        CALL CIKVA(V,Z,VM,CBI,CDI,CBK,CDK)
        NM=INT(VM)
        WRITE(*,*)
        WRITE(*,*)'   n      Re[Iv(z)]       Im[Iv(z)] ',
     &            '     Re[Iv''(z)]      Im[Iv''(z)] '
        WRITE(*,*)' ---------------------------------------',
     &            '------------------------------'
        DO K=0,NM,NS
           WRITE(*,20)K,CBI(K),CDI(K)
      end do
        WRITE(*,*)
        WRITE(*,*)'   n      Re[Kv(z)]       Im[Kv(z)] ',
     &            '     Re[Kv''(z)]      Im[Kv''(z)] '
        WRITE(*,*)' ---------------------------------------',
     &            '------------------------------'
        DO K=0,NM,NS
           WRITE(*,20)K,CBK(K),CDK(K)
      end do
20      FORMAT(1X,I4,1X,4D16.8)
25      FORMAT(8X,'v= n+v0',',   ','v0 =',F7.2,',   ','z =',F6.1,
     &        '+ i',F6.1)

      return
        END
      subroutine mcikvb ( )

c*********************************************************************72
c
cc MCIKVB tests CIKVB.
c
c  Modified:
c
c    19 March 2012
c
c
c       Purpose: This program computes the modified Bessel functions 
c                Iv(z), Kv(z) and their derivatives for an arbitrary
c                order and complex argument using subroutine CIKVB 
c       Input :  z --- Complex argument z
c                v --- Real order of Iv(z) and Kv(z)
c                      ( v = n+v0,  0 ó n ó 250, 0 ó v0 < 1 )
c       Output:  CBI(n) --- In+v0(z)
c                CDI(n) --- In+v0'(z)
c                CBK(n) --- Kn+v0(z)
c                CDK(n) --- Kn+v0'(z)
c       Example: v= n+v0,   v0 = .25,   z = 4.0+ i 2.0
c
c    n     Re[Iv(z)]      Im[Iv(z)]     Re[Iv'(z)]     Im[Iv'(z)]
c  -----------------------------------------------------------------
c    0  -.19336550D+01  .10328998D+02 -.23119621D+01  .91612230D+01
c    1  -.24735044D+01  .85964317D+01 -.23898329D+01  .78707023D+01
c    2  -.28460107D+01  .54124063D+01 -.24105909D+01  .55204965D+01
c    3  -.23476775D+01  .24445612D+01 -.21145027D+01  .30604463D+01
c    4  -.13829947D+01  .70848630D+00 -.14732387D+01  .12545751D+01
c    5  -.59879982D+00  .64588999D-01 -.78816416D+00  .32629794D+00
c
c    n     Re[Kv(z)]      Im[Kv(z)]     Re[Kv'(z)]     Im[Kv'(z)]
c   ----------------------------------------------------------------
c    0  -.64820386D-02 -.84715754D-02  .75118612D-02  .89920077D-02
c    1  -.80477525D-02 -.92535355D-02  .96506687D-02  .97789903D-02
c    2  -.12819299D-01 -.11086405D-01  .16310878D-01  .11358076D-01
c    3  -.24574004D-01 -.13462616D-01  .33167751D-01  .11850554D-01
c    4  -.53516204D-01 -.12614703D-01  .75424026D-01  .14407268D-02
c    5  -.12627405D+00  .10581162D-01  .18054884D+00 -.64789392D-01
c
c
        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CBI(0:250),CDI(0:250),CBK(0:250),CDK(0:250)
 
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MCIKVB'
      write ( *, '(a)' ) '  Test CIKVB'
      write ( *, '(a)' ) ' '

      v = 0.25D+00
      x = 4.0D+00
      y = 2.0D+00
      ns = 5

        Z=CMPLX(X,Y)
        N=INT(V)
        V0=V-N
        WRITE(*,25)V0,X,Y
        IF (N.LE.8) THEN
           NS=1
        ENDIF
        CALL CIKVB(V,Z,VM,CBI,CDI,CBK,CDK)
        NM=INT(VM)
        WRITE(*,*)
        WRITE(*,*)'   n      Re[Iv(z)]       Im[Iv(z)] ',
     &            '     Re[Iv''(z)]      Im[Iv''(z)] '
        WRITE(*,*)' ---------------------------------------',
     &            '------------------------------'
        DO K=0,NM,NS
           WRITE(*,20) K,CBI(K),CDI(K)
        end do
        WRITE(*,*)
        WRITE(*,*)'   n      Re[Kv(z)]       Im[Kv(z)] ',
     &            '     Re[Kv''(z)]      Im[Kv''(z)] '
        WRITE(*,*)' ---------------------------------------',
     &            '------------------------------'
        DO K=0,NM,NS
           WRITE(*,20)K,CBK(K),CDK(K)
        end do
20      FORMAT(1X,I4,1X,4D16.8)
25      FORMAT(8X,'v= n+v0',',   ','v0 =',F6.2,',   ','z =',F5.1,
     &        '+ i',F5.1)
      return
        END
      subroutine mcisia ( )

c*********************************************************************72
c
cc MCISIA tests CISIA.
c
c  Modified:
c
c    19 March 2012
c
c
c       Purpose: This program computes the cosine and sine 
c                integrals using subroutine CISIA
c       Input :  x  --- Argument of Ci(x) and Si(x)
c       Output:  CI --- Ci(x)
c                SI --- Si(x)
c       Example:
c                    x         Ci(x)          Si(x)
c                 ------------------------------------
c                   0.0     - ì             .00000000
c                   5.0     -.19002975     1.54993124
c                  10.0     -.04545643     1.65834759
c                  20.0      .04441982     1.54824170
c                  30.0     -.03303242     1.56675654
c                  40.0      .01902001     1.58698512
c
c
        DOUBLE PRECISION CI,SI,X

      integer test_num
      parameter ( test_num = 6 )

      integer test
      double precision x_test(test_num)

      save x_test

      data x_test / 
     &  0.0D+00, 5.0D+00, 10.0D+00, 20.0D+00, 30.0D+00, 40.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MCISIA'
      write ( *, '(a)' ) '  Test CISIA'
      write ( *, '(a)' ) ' '

        WRITE(*,*)'   x         Ci(x)          Si(x)'
        WRITE(*,*)'------------------------------------'

      do test = 1, test_num

        x = x_test(test)

        CALL CISIA(X,CI,SI)

        IF (X.NE.0.0D0) then
          WRITE(*,10)X,CI,SI
        else
          WRITE(*,20)
        end if

      end do

10      FORMAT(1X,F5.1,2F15.8)
20      FORMAT(3X,' .0',4X,' - ì',13X,'.00000000')
       return
        END
      subroutine mcisib ( )

c*********************************************************************72
c
cc MCISIB tests CISIB.
c
c  Modified:
c
c    01 July 2012
c
c  Example:
c
c      x        Ci(x)           Si(x)
c    ------------------------------------
c     0.0    - oo                 0
c     5.0    -.190030D+00      1.549931
c    10.0    -.454563D-01      1.658348
c    20.0     .444201D-01      1.548241
c    30.0    -.330326D-01      1.566757
c    40.0     .190201D-01      1.586985
c
      implicit none

      integer test_num
      parameter ( test_num = 6 )

      double precision ci
      double precision si
      integer test
      double precision x
      double precision x_test(test_num)

      save x_test

      data x_test / 
     &  0.0D+00, 5.0D+00, 10.0D+00, 20.0D+00, 30.0D+00, 40.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MCISIB'
      write ( *, '(a)' ) '  Test CISIB, which computes the'
      write ( *, '(a)' ) '  cosine and sine integrals.'
      write ( *, '(a)' ) ' '

      write ( *, '(a)' ) '   x        ci(x)           si(x)'
      write ( *, '(a)' ) '------------------------------------'

      do test = 1, test_num

        x = x_test(test)

        call cisib ( x, ci, si )

        write ( *, '(1x,f5.1,g16.8,g16.8)' ) x, ci, si

      end do

      return
      end
      subroutine mcjk ( )

c*********************************************************************72
c
cc MCJK tests CJK.
c
c  Modified:
c
c    20 March 2012
c
c
c       Purpose: This program computes the expansion coefficients  
c                for the asymptotic expansion of Bessel functions  
c                with large orders using subroutine CJK
c       Input :  Km   --- Maximum k
c       Output:  A(L) --- Cj(k) where j and k are related to L by
c                         L=j+1+[k*(k+1)]/2; j,k=0,1,2,...,Km
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

        DIMENSION A(231)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MCJK'
      write ( *, '(a)' ) '  Test CJK'
      write ( *, '(a)' ) ' '

      km = 15
      write ( *, '(a,i8)' ) '  Maximum K = ', km

        LM=KM+1+(KM*(KM+1))/2
        CALL CJK(KM,A)
        DO K=1,LM
           WRITE(*,15)K,A(K)
      end do
15      FORMAT(1X,I3,D25.14)
      return
        END
      subroutine mcjy01 ( )

c*********************************************************************72
c
cc MCJY01 tests CJY01.
c
c  Modified:
c
c    20 March 2012
c
c
c       Purpose: This program computes Bessel functions J0(z), J1(z),
c                Y0(z), Y1(z), and their derivatives for a complex
c                argument using subroutine CJY01
c       Input :  z --- Complex argument
c       Output:  CBJ0 --- J0(z)
c                CDJ0 --- J0'(z)
c                CBJ1 --- J1(z)
c                CDJ1 --- J1'(z)
c                CBY0 --- Y0(z)
c                CDY0 --- Y0'(z)
c                CBY1 --- Y1(z)
c                CDY1 --- Y1'(z)
c       Example: z =  4.0 + i  2.0
c
c     n     Re[Jn(z)]       Im[Jn(z)]       Re[Jn'(z)]      Im[Jn'(z)]
c   --------------------------------------------------------------------
c     0  -.13787022D+01   .39054236D+00   .50735255D+00   .12263041D+01
c     1  -.50735255D+00  -.12263041D+01  -.11546013D+01   .58506793D+00
c
c     n     Re[Yn(z)]       Im[Yn(z)]       Re[Yn'(z)]      Im[Yn'(z)]
c   --------------------------------------------------------------------
c     0  -.38145893D+00  -.13291649D+01  -.12793101D+01   .51220420D+00
c     1   .12793101D+01  -.51220420D+00  -.58610052D+00  -.10987930D+01
c
c
        IMPLICIT DOUBLE PRECISION (X,Y)
        IMPLICIT COMPLEX*16 (C,Z)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MCJY01'
      write ( *, '(a)' ) '  Test CJY01'
      write ( *, '(a)' ) ' '

      x = 4.0D+00
      y = 2.0D+00

        Z=CMPLX(X,Y)
        WRITE(*,30)X,Y
        CALL CJY01(Z,CBJ0,CDJ0,CBJ1,CDJ1,CBY0,CDY0,CBY1,CDY1)
        WRITE(*,*)
        WRITE(*,*)'  n      Re[Jn(z)]       Im[Jn(z)]',
     &            '       Re[Jn''(z)]      Im[Jn''(z)]'
        WRITE(*,*)' -------------------------------------',
     &            '-------------------------------'
        WRITE(*,10)CBJ0,CDJ0
        WRITE(*,20)CBJ1,CDJ1
        WRITE(*,*)
        WRITE(*,*)'  n      Re[Yn(z)]       Im[Yn(z)]',
     &            '       Re[Yn''(z)]      Im[Yn''(z)]'
        WRITE(*,*)' -------------------------------------',
     &            '-------------------------------'
        WRITE(*,10)CBY0,CDY0
        WRITE(*,20)CBY1,CDY1
10      FORMAT(3X,'0',2X,4D16.8)
20      FORMAT(3X,'1',2X,4D16.8)
30      FORMAT(3X,3Hz =,F6.2,' + i',F6.2)
      return
        END
      subroutine MCJYLV ( )

c*********************************************************************72
c
cc MCJYLV tests CJYLV
c
c  Modified:
c
c    22 March 2012
c
c
c       Purpose: This program computes Bessel functions Jv(z) 
c                and Yv(z) and their derivatives with a large 
c                order and complex argument using subroutine
c                CJYLV
c       Input:   v --- Order of Jv(z) and Yv(z)
c                z --- Complex argument
c       Output:  CBJV --- Jv(z)
c                CDJV --- Jv'(z)
c                CBYV --- Yv(z)
c                CDYV --- Yv'(z)
c       Examples:
c                v = 100.00,    z = 4.00 + 2.00 i
c
c                Jv(z) = -.6444792518-123 + .6619157435-123 i
c                Jv'(z)= -.6251103777-122 + .1967638668-121 i
c                Yv(z) =  .2403065353+121 + .2472039414+121 i
c                Yv'(z)= -.7275814786+122 - .2533588851+122 i
c
c                v =100.5,     z = 4.00 + 2.00 i
c
c                Jv(z) = -.1161315754-123 + .7390127781-124 i
c                Jv'(z)= -.1588519437-122 + .2652227059-122 i
c                Yv(z) =  .1941381412+122 + .1237578195+122 i
c                Yv'(z)= -.5143285247+123 - .5320026773+122 i
c
c
        IMPLICIT DOUBLE PRECISION (V,X,Y)
        IMPLICIT COMPLEX*16 (C,Z)

      integer test_num
      parameter ( test_num = 2 )

      integer test
      double precision v
      double precision v_test(test_num)
      double precision x
      double precision x_test(test_num)
      double precision y
      double precision y_test(test_num)

      save v_test
      save x_test
      save y_test

      data v_test / 100.0D+00, 100.5D+00 /
      data x_test / 4.0D+00, 4.0D+00 /
      data y_test / 2.0D+00, 2.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MCJYLV'
      write ( *, '(a)' ) '  Test CJYLV'
      write ( *, '(a)' ) ' '

      do test = 1, test_num

        v = v_test(test)
        x = x_test(test)
        y = y_test(test)
        WRITE(*,10)V,X,Y
        Z=CMPLX(X,Y)
        CALL CJYLV(V,Z,CBJV,CDJV,CBYV,CDYV)
        WRITE(*,*)
        WRITE(*,20)CBJV
        WRITE(*,30)CDJV
        WRITE(*,*)
        WRITE(*,40)CBYV
        WRITE(*,50)CDYV
      end do

10      FORMAT(8X,'v = ',F6.2,',    ','z =',F7.2,' + i ',F7.2)
20      FORMAT(8X,'Jv(z) =',D17.10,' + i',D17.10)
30      FORMAT(8X,'Jv''(z)=',D17.10,' + i',D17.10)
40      FORMAT(8X,'Yv(z) =',D17.10,' + i',D17.10)
50      FORMAT(8X,'Yv''(z)=',D17.10,' + i',D17.10)
      return
        END
      subroutine MCJYNA ( )

c*********************************************************************72
c
cc MCJYNA tests CJYNA
c
c  Modified:
c
c    22 March 2012
c
c
c       Purpose: This program computes Bessel functions Jn(z), Yn(z)  
c                and their derivatives for a complex argument using
c              subroutine CJYNA
c       Input :  z --- Complex argument of Jn(z) and Yn(z)
c                n --- Order of Jn(z) and Yn(z)
c                      ( n = 0,1,úúú, n ó 250 )
c       Output:  CBJ(n) --- Jn(z)
c                CDJ(n) --- Jn'(z)
c                CBY(n) --- Yn(z)
c                CDY(n) --- Yn'(z)
c       Eaxmple: z = 4.0 + i 2.0
c
c     n     Re[Jn(z)]       Im[Jn(z)]       Re[Jn'(z)]      Im[Jn'(z)]
c    -------------------------------------------------------------------
c     0  -.13787022D+01   .39054236D+00   .50735255D+00   .12263041D+01
c     1  -.50735255D+00  -.12263041D+01  -.11546013D+01   .58506793D+00
c     2   .93050039D+00  -.77959350D+00  -.72363400D+00  -.72836666D+00
c     3   .93991546D+00   .23042918D+00   .29742236D+00  -.63587637D+00
c     4   .33565567D+00   .49215925D+00   .47452722D+00  -.29035945D-01
c     5  -.91389835D-02   .28850107D+00   .20054412D+00   .19908868D+00
c
c     n     Re[Yn(z)]       Im[Yn(z)]       Re[Yn'(z)]      Im[Yn'(z)]
c   --------------------------------------------------------------------
c     0  -.38145893D+00  -.13291649D+01  -.12793101D+01   .51220420D+00
c     1   .12793101D+01  -.51220420D+00  -.58610052D+00  -.10987930D+01
c     2   .79074211D+00   .86842120D+00   .78932897D+00  -.70142425D+00
c     3  -.29934789D+00   .89064431D+00   .70315755D+00   .24423024D+00
c     4  -.61557299D+00   .37996071D+00   .41126221D-01   .34044655D+00
c     5  -.38160033D+00   .20975121D+00  -.33884827D+00  -.20590670D-01
c
c                z = 20.0 + i 10.0 ,      Nmax = 5
c
c     n     Re[Jn(z)]       Im[Jn(z)]       Re[Jn'(z)]      Im[Jn'(z)]
c   --------------------------------------------------------------------
c     0   .15460268D+04  -.10391216D+04  -.10601232D+04  -.15098284D+04
c     1   .10601232D+04   .15098284D+04   .14734253D+04  -.10783122D+04
c     2  -.14008238D+04   .11175029D+04   .11274890D+04   .13643952D+04
c     3  -.11948548D+04  -.12189620D+04  -.11843035D+04   .11920871D+04
c     4   .96778325D+03  -.12666712D+04  -.12483664D+04  -.93887194D+03
c     5   .13018781D+04   .65878188D+03   .64152944D+03  -.12682398D+04
c
c     n     Re[Yn(z)]       Im[Yn(z)]       Re[Yn'(z)]      Im[Yn'(z)]
c   --------------------------------------------------------------------
c     0   .10391216D+04   .15460268D+04   .15098284D+04  -.10601232D+04
c     1  -.15098284D+04   .10601232D+04   .10783122D+04   .14734253D+04
c     2  -.11175029D+04  -.14008238D+04  -.13643952D+04   .11274890D+04
c     3   .12189620D+04  -.11948548D+04  -.11920871D+04  -.11843035D+04
c     4   .12666712D+04   .96778324D+03   .93887194D+03  -.12483664D+04
c     5  -.65878189D+03   .13018781D+04   .12682398D+04   .64152944D+03
c
c
        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        COMMON CBJ(0:251),CDJ(0:251),CBY(0:251),CDY(0:251)

      integer test_num
      parameter ( test_num = 2 )

      integer test
      integer n
      integer n_test(test_num)
      double precision x
      double precision x_test(test_num)
      double precision y
      double precision y_test(test_num)

      save n_test
      save x_test
      save y_test

      data n_test / 5, 5 /
      data x_test / 4.0D+00, 20.0D+00 /
      data y_test / 2.0D+00, 10.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MCJYNA'
      write ( *, '(a)' ) '  Test CJYNA'
      write ( *, '(a)' ) ' '

      do test = 1, test_num

        n = n_test(test)
        x = x_test(test)
        y = y_test(test)
        ns = 5

        Z=CMPLX(X,Y)
        WRITE(*,40)X,Y,N
        IF (N.LE.8) THEN
           NS=1
        ENDIF
        CALL CJYNA(N,Z,NM,CBJ,CDJ,CBY,CDY)
        WRITE(*,*)
        WRITE(*,*)'   n     Re[Jn(z)]       Im[Jn(z)]',
     &            '       Re[Jn''(z)]      Im[Jn''(z)]'
        WRITE(*,*)' -------------------------------------',
     &            '-------------------------------'
        DO 10 K=0,NM,NS
10         WRITE(*,30)K,CBJ(K),CDJ(K)
        WRITE(*,*)
        WRITE(*,*)'   n     Re[Yn(z)]       Im[Yn(z)]',
     &            '       Re[Yn''(z)]      Im[Yn''(z)]'
        WRITE(*,*)' -------------------------------------',
     &            '-------------------------------'
        DO 20 K=0,NM,NS
20         WRITE(*,30)K,CBY(K),CDY(K)

      end do

30      FORMAT(1X,I4,4D16.8)
40      FORMAT(3X,3Hz =,F5.1,' + i ',F5.1,' ,',6X,6HNmax =,I3)
      return
        END
      subroutine MCJYNB ( )

c*********************************************************************72
c
cc MCJYNB tests CJYNB
c
c  Modified:
c
c    23 March 2012
c
c
c       Purpose: This program computes Bessel functions Jn(z), Yn(z)  
c                and their derivatives for a complex argument using
c              subroutine CJYNB
c       Input :  z --- Complex argument of Jn(z) and Yn(z)
c                n --- Order of Jn(z) and Yn(z)
c                      ( n = 0,1,úúú, n ó 250 )
c       Output:  CBJ(n) --- Jn(z)
c                CDJ(n) --- Jn'(z)
c                CBY(n) --- Yn(z)
c                CDY(n) --- Yn'(z)
c       Eaxmple: z = 4.0 + i 2.0
c
c     n     Re[Jn(z)]       Im[Jn(z)]       Re[Jn'(z)]      Im[Jn'(z)]
c    -------------------------------------------------------------------
c     0  -.13787022D+01   .39054236D+00   .50735255D+00   .12263041D+01
c     1  -.50735255D+00  -.12263041D+01  -.11546013D+01   .58506793D+00
c     2   .93050039D+00  -.77959350D+00  -.72363400D+00  -.72836666D+00
c     3   .93991546D+00   .23042918D+00   .29742236D+00  -.63587637D+00
c     4   .33565567D+00   .49215925D+00   .47452722D+00  -.29035945D-01
c     5  -.91389835D-02   .28850107D+00   .20054412D+00   .19908868D+00
c
c     n     Re[Yn(z)]       Im[Yn(z)]       Re[Yn'(z)]      Im[Yn'(z)]
c   --------------------------------------------------------------------
c     0  -.38145893D+00  -.13291649D+01  -.12793101D+01   .51220420D+00
c     1   .12793101D+01  -.51220420D+00  -.58610052D+00  -.10987930D+01
c     2   .79074211D+00   .86842120D+00   .78932897D+00  -.70142425D+00
c     3  -.29934789D+00   .89064431D+00   .70315755D+00   .24423024D+00
c     4  -.61557299D+00   .37996071D+00   .41126221D-01   .34044655D+00
c     5  -.38160033D+00   .20975121D+00  -.33884827D+00  -.20590670D-01
c
c                z = 20.0 + i  10.0 ,      Nmax =  5
c
c     n     Re[Jn(z)]       Im[Jn(z)]       Re[Jn'(z)]      Im[Jn'(z)]
c   --------------------------------------------------------------------
c     0   .15460268D+04  -.10391216D+04  -.10601232D+04  -.15098284D+04
c     1   .10601232D+04   .15098284D+04   .14734253D+04  -.10783122D+04
c     2  -.14008238D+04   .11175029D+04   .11274890D+04   .13643952D+04
c     3  -.11948548D+04  -.12189620D+04  -.11843035D+04   .11920871D+04
c     4   .96778325D+03  -.12666712D+04  -.12483664D+04  -.93887194D+03
c     5   .13018781D+04   .65878188D+03   .64152944D+03  -.12682398D+04
c
c     n     Re[Yn(z)]       Im[Yn(z)]       Re[Yn'(z)]      Im[Yn'(z)]
c   --------------------------------------------------------------------
c     0   .10391216D+04   .15460268D+04   .15098284D+04  -.10601232D+04
c     1  -.15098284D+04   .10601232D+04   .10783122D+04   .14734253D+04
c     2  -.11175029D+04  -.14008238D+04  -.13643952D+04   .11274890D+04
c     3   .12189620D+04  -.11948548D+04  -.11920871D+04  -.11843035D+04
c     4   .12666712D+04   .96778324D+03   .93887194D+03  -.12483664D+04
c     5  -.65878189D+03   .13018781D+04   .12682398D+04   .64152944D+03
c
c
        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        COMMON CBJ(0:250),CDJ(0:250),CBY(0:250),CDY(0:250)

      integer test_num
      parameter ( test_num = 2 )

      integer test
      integer n
      integer n_test(test_num)
      double precision x
      double precision x_test(test_num)
      double precision y
      double precision y_test(test_num)

      save n_test
      save x_test
      save y_test

      data n_test / 5, 5 /
      data x_test / 4.0D+00, 20.0D+00 /
      data y_test / 2.0D+00, 10.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MCJYNB'
      write ( *, '(a)' ) '  Test CJYNB'
      write ( *, '(a)' ) ' '

      do test = 1, test_num

        n = n_test(test)
        x = x_test(test)
        y = y_test(test)
        ns = 5

        Z=CMPLX(X,Y)
        WRITE(*,40)X,Y,N
        IF (N.LE.8) THEN
           NS=1
        ELSE
           NS = 5
        ENDIF
        CALL CJYNB(N,Z,NM,CBJ,CDJ,CBY,CDY)
        WRITE(*,*)
        WRITE(*,*)'   n     Re[Jn(z)]       Im[Jn(z)]',
     &            '       Re[Jn''(z)]      Im[Jn''(z)]'
        WRITE(*,*)' -------------------------------------',
     &            '-------------------------------'
        DO 10 K=0,NM,NS
10         WRITE(*,30)K,CBJ(K),CDJ(K)
        WRITE(*,*)
        WRITE(*,*)'   n     Re[Yn(z)]       Im[Yn(z)]',
     &            '       Re[Yn''(z)]      Im[Yn''(z)]'
        WRITE(*,*)' -------------------------------------',
     &            '-------------------------------'
        DO 20 K=0,NM,NS
20         WRITE(*,30)K,CBY(K),CDY(K)
      end do

30      FORMAT(1X,I4,4D16.8)
40      FORMAT(3X,3Hz =,F5.1,' + i ',F5.1,' ,',6X,6HNmax =,I3)
      return
        END
      subroutine MCJYVA ( )

c*********************************************************************72
c
cc MCJYVA tests CJYVA
c
c  Modified:
c
c    23 March 2012
c
c
c       Purpose: This program computes Bessel functions Jv(z), Yv(z),
c                and their derivatives for a complex argument using
c              subroutine CJYVA
c       Input :  z --- Complex argument
c                v --- Order of Jv(z) and Yv(z)
c                      ( v = n+v0, 0 ó n ó 250, 0 ó v0 < 1 )
c       Output:  CBJ(n) --- Jn+v0(z)
c                CDJ(n) --- Jn+v0'(z)
c                CBY(n) --- Yn+v0(z)
c                CDY(n) --- Yn+v0'(z)
c       Example:
c                v = n +v0,  v0 = 1/3,   z = 4.0 + i 2.0
c
c     n     Re[Jv(z)]       Im[Jv(z)]      Re[Jv'(z)]      Im[Jv'(z)]
c    ------------------------------------------------------------------
c     0  -.13829878D+01  -.30855145D+00  -.18503756D+00   .13103689D+01
c     1   .82553327D-01  -.12848394D+01  -.12336901D+01   .45079506D-01
c     2   .10843924D+01  -.39871046D+00  -.33046401D+00  -.84574964D+00
c     3   .74348135D+00   .40665987D+00   .45318486D+00  -.42198992D+00
c     4   .17802266D+00   .44526939D+00   .39624497D+00   .97902890D-01
c     5  -.49008598D-01   .21085409D+00   .11784299D+00   .19422044D+00
c
c     n     Re[Yv(z)]      Im[Yv(z)]       Re[Yv'(z)]      Im[Yv'(z)]
c    ------------------------------------------------------------------
c     0   .34099851D+00  -.13440666D+01  -.13544477D+01  -.15470699D+00
c     1   .13323787D+01   .53735934D-01  -.21467271D-01  -.11807457D+01
c     2   .38393305D+00   .10174248D+01   .91581083D+00  -.33147794D+00
c     3  -.49924295D+00   .71669181D+00   .47786442D+00   .37321597D+00
c     4  -.57179578D+00   .27099289D+00  -.12111686D+00   .23405313D+00
c     5  -.25700924D+00   .24858555D+00  -.43023156D+00  -.13123662D+00
c
c
        IMPLICIT DOUBLE PRECISION (V,X,Y)
        IMPLICIT COMPLEX*16 (C,Z)
        COMMON CBJ(0:251),CDJ(0:251),CBY(0:251),CDY(0:251)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MCJYVA'
      write ( *, '(a)' ) '  Test CJYVA'
      write ( *, '(a)' ) ' '

      v = 1.0D+00 / 3.0D+00
      x = 4.0D+00
      y = 2.0D+00

        Z=CMPLX(X,Y)
        N=INT(V)
        V0=V-N
        WRITE(*,25)V0,X,Y
        IF (N.LE.8) THEN
           NS=1
        ELSE
           ns = 5
        ENDIF
        CALL CJYVA(V,Z,VM,CBJ,CDJ,CBY,CDY)
        NM=INT(VM)
        WRITE(*,*)
        WRITE(*,*)'  n       Re[Jv(z)]       Im[Jv(z)]',
     &            '       Re[Jv''(z)]      Im[Jv''(z)]'
        WRITE(*,*)' ----------------------------------',
     &            '-----------------------------------'
        DO K=0,NM,NS
          WRITE(*,20) K,CBJ(K),CDJ(K)
        end do
        WRITE(*,*)
        WRITE(*,*)'  n       Re[Yv(z)]       Im[Yv(z)]',
     &            '       Re[Yv''(z)]      Im[Yv''(z)]'
        WRITE(*,*)' ----------------------------------',
     &            '-----------------------------------'
        DO K=0,NM,NS
          WRITE(*,20) K,CBY(K),CDY(K)
        end do
20      FORMAT(1X,I3,2X,4D16.8)
25      FORMAT(8X,'v = n+v0',',  v0 =',F5.2,',  z =',F7.2,' +',F7.2,'i')

      return
        END
      subroutine MCJYVB ( )

c*********************************************************************72
c
cc MCJYVB tests CJYVB
c
c  Modified:
c
c    24 March 2012
c
c
c       Purpose: This program computes Bessel functions Jv(z), Yv(z),
c                and their derivatives for a complex argument using
c              subroutine CJYVB
c       Input :  z --- Complex argument
c                v --- Order of Jv(z) and Yv(z)
c                      ( v = n+v0, 0 ó n ó 250, 0 ó v0 < 1 )
c       Output:  CBJ(n) --- Jn+v0(z)
c                CDJ(n) --- Jn+v0'(z)
c                CBY(n) --- Yn+v0(z)
c                CDY(n) --- Yn+v0'(z)
c       Example:
c                v = n +v0,  v0 = 1/3,   z = 4.0 + i 2.0
c
c     n     Re[Jv(z)]       Im[Jv(z)]      Re[Jv'(z)]      Im[Jv'(z)]
c    -------------------------------------------------------------------
c     0  -.13829878D+01  -.30855145D+00  -.18503756D+00   .13103689D+01
c     1   .82553327D-01  -.12848394D+01  -.12336901D+01   .45079506D-01
c     2   .10843924D+01  -.39871046D+00  -.33046401D+00  -.84574964D+00
c     3   .74348135D+00   .40665987D+00   .45318486D+00  -.42198992D+00
c     4   .17802266D+00   .44526939D+00   .39624497D+00   .97902890D-01
c     5  -.49008598D-01   .21085409D+00   .11784299D+00   .19422044D+00
c
c     n     Re[Yv(z)]      Im[Yv(z)]       Re[Yv'(z)]      Im[Yv'(z)]
c    -------------------------------------------------------------------
c     0   .34099851D+00  -.13440666D+01  -.13544477D+01  -.15470699D+00
c     1   .13323787D+01   .53735934D-01  -.21467271D-01  -.11807457D+01
c     2   .38393305D+00   .10174248D+01   .91581083D+00  -.33147794D+00
c     3  -.49924295D+00   .71669181D+00   .47786442D+00   .37321597D+00
c     4  -.57179578D+00   .27099289D+00  -.12111686D+00   .23405313D+00
c     5  -.25700924D+00   .24858555D+00  -.43023156D+00  -.13123662D+00
c
c
        IMPLICIT DOUBLE PRECISION (V,X,Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CBJ(0:250),CDJ(0:250),CBY(0:250),CDY(0:250)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MCJYVB'
      write ( *, '(a)' ) '  Test CJYVB'
      write ( *, '(a)' ) ' '

      v = 1.0D+00 / 3.0D+00
      x = 4.0D+00
      y = 2.0D+00

        Z=CMPLX(X,Y)
        N=INT(V)
        V0=V-N
        WRITE(*,25)V0,X,Y
        IF (N.LE.8) THEN
           NS=1
        ELSE
           NS = 5
        ENDIF
        CALL CJYVB(V,Z,VM,CBJ,CDJ,CBY,CDY)
        NM=INT(VM)
        WRITE(*,*)
        WRITE(*,*)'  n       Re[Jv(z)]       Im[Jv(z)]',
     &            '       Re[Jv''(z)]      Im[Jv''(z)]'
        WRITE(*,*)' ----------------------------------',
     &            '-----------------------------------'
        DO 10 K=0,NM,NS
10         WRITE(*,20) K,CBJ(K),CDJ(K)
        WRITE(*,*)
        WRITE(*,*)'  n       Re[Yv(z)]       Im[Yv(z)]',
     &            '       Re[Yv''(z)]      Im[Yv''(z)]'
        WRITE(*,*)' ----------------------------------',
     &            '-----------------------------------'
        DO 15 K=0,NM,NS
15         WRITE(*,20) K,CBY(K),CDY(K)
20      FORMAT(1X,I3,2X,4D16.8)
25      FORMAT(8X,'v = n+v0',',  v0 =',F5.2,',  z =',F7.2,' +',F7.2,'i')
      return
        END
      subroutine MCLPMN ( )

c*********************************************************************72
c
cc MCLPMN tests CLPMN.
c
c  Modified:
c
c    24 March 2012
c
c
c       Purpose: This program computes the associated Legendre 
c                functions Pmn(z) and their derivatives Pmn'(z) for
c                a complex argument using subroutine CLPMN
c       Input :  x --- Real part of z
c                y --- Imaginary part of z
c                m --- Order of Pmn(z),  m = 0,1,2,...,n
c                n --- Degree of Pmn(z), n = 0,1,2,...,N
c       Output:  CPM(m,n) --- Pmn(z)
c                CPD(m,n) --- Pmn'(z)
c       Examples:
c                n = 5, x = 0.5, y = 0.2
c
c       m     Re[Pmn(z)]    Im[Pmn(z)]    Re[Pmn'(z)]   Im[Pmn'(z)]
c      -------------------------------------------------------------
c       0    .252594D+00  -.530293D+00   -.347606D+01  -.194250D+01
c       1    .333071D+01   .135206D+01    .117643D+02  -.144329D+02
c       2   -.102769D+02   .125759D+02    .765713D+02   .598500D+02
c       3   -.572879D+02  -.522744D+02   -.343414D+03   .147389D+03
c       4    .335711D+03  -.389151D+02   -.226328D+03  -.737100D+03
c       5   -.461125D+03   .329122D+03    .187180D+04   .160494D+02
c
c                n = 5, x = 2.5, y = 1.0
c
c       m     Re[Pmn(z)]    Im[Pmn(z)]    Re[Pmn'(z)]   Im[Pmn'(z)]
c      -------------------------------------------------------------
c       0   -.429395D+03   .900336D+03   -.350391D+02   .193594D+04
c       1   -.216303D+04   .446358D+04   -.208935D+03   .964685D+04
c       2   -.883477D+04   .174005D+05   -.123703D+04   .381938D+05
c       3   -.273211D+05   .499684D+05   -.568080D+04   .112614D+06
c       4   -.565523D+05   .938503D+05   -.167147D+05   .219713D+06
c       5   -.584268D+05   .863328D+05   -.233002D+05   .212595D+06
c
c
        IMPLICIT DOUBLE PRECISION (X,Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CPM(0:40,0:40),CPD(0:40,0:40)

      integer test_num
      parameter ( test_num = 12 )

      integer test
      integer m
      integer m_test(test_num)
      integer n
      integer n_test(test_num)
      double precision x
      double precision x_test(test_num)
      double precision y
      double precision y_test(test_num)

      save m_test
      save n_test
      save x_test
      save y_test

      data m_test / 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5 /
      data n_test / 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5 /
      data x_test / 
     &  0.5D+00, 0.5D+00, 0.5D+00, 0.5D+00, 0.5D+00, 0.5D+00,
     &  2.5D+00, 2.5D+00, 2.5D+00, 2.5D+00, 2.5D+00, 2.5D+00 /
      data y_test / 
     &  0.2D+00, 0.2D+00, 0.2D+00, 0.2D+00, 0.2D+00, 0.2D+00,
     &  1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MCLPMN'
      write ( *, '(a)' ) '  Test CLPMN'
      write ( *, '(a)' ) ' '

      do test = 1, test_num

        m = m_test(test)
        n = n_test(test)
        x = x_test(test)
        y = y_test(test)

        WRITE(*,30) M,N,X,Y
        CALL CLPMN(40,M,N,X,Y,CPM,CPD)
        WRITE(*,*)'   m   n    Re[Pmn(z)]    Im[Pmn(z)]    ',
     &            'Re[Pmn''(z)]   Im[Pmn''(z)]'
        WRITE(*,*)' -----------------------------------',
     &            '-------------------------------'
        DO 10 J=0,N
10         WRITE(*,20)M,J,CPM(M,J),CPD(M,J)
20      FORMAT(1X,2I4,1X,2D14.6,1X,2D14.6)
30      FORMAT(1X,'m =',I2,', ','n =',I2,', ','x =',F5.1,
     &         ', ','y =',F5.1)

      end do

      return
        END
      subroutine MCLPN ( )

c*********************************************************************72
c
cc MCLPN tests CLPN.
c
c  Modified:
c
c    24 March 2012
c
c
c       Purpose: This program computes the Legendre polynomials 
c                Pn(z) and Pn'(z) for a complex argument using
c              subroutine CLPN
c       Input :  x --- Real part of z
c                y --- Imaginary part of z
c                n --- Degree of Pn(z), n = 0,1,...,N
c       Output:  CPN(n) --- Pn(z)
c                CPD(n) --- Pn'(z)
c       Example: z = 3.0 +2.0 i
c
c       n    Re[Pn(z)]     Im[Pn(z)]     Re[Pn'(z)]   Im[Pn'(z)]
c      -----------------------------------------------------------
c       0   .100000D+01   .000000D+00   .000000D+00   .000000D+00
c       1   .300000D+01   .200000D+01   .100000D+01   .000000D+00
c       2   .700000D+01   .180000D+02   .900000D+01   .600000D+01
c       3  -.270000D+02   .112000D+03   .360000D+02   .900000D+02
c       4  -.539000D+03   .480000D+03  -.180000D+03   .790000D+03
c       5  -.461700D+04   .562000D+03  -.481500D+04   .441000D+04
c
c
        IMPLICIT DOUBLE PRECISION (X,Y)
        IMPLICIT COMPLEX *16 (C,Z)
        DIMENSION CPN(0:100),CPD(0:100)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MCLPN'
      write ( *, '(a)' ) '  Test CLPN'
      write ( *, '(a)' ) ' '

      n = 5
      x = 3.0D+00
      y = 2.0D+00

        WRITE(*,30)X,Y
        WRITE(*,*)
        CALL CLPN(N,X,Y,CPN,CPD)
        WRITE(*,*)'  n    Re[Pn(z)]     Im[Pn(z)]     Re[Pn''(z)]',
     &            '   Im[Pn''(z)]'
        WRITE(*,*)' ---------------------------------------------',
     &            '--------------'
        DO 10 K=0,N
10         WRITE(*,20)K,CPN(K),CPD(K)
20      FORMAT(1X,I3,4D14.6)
30      FORMAT(3X,'x =',F5.1,',  ','y =',F5.1)

      return
        END
      subroutine MCLQMN ( )

c*********************************************************************72
c
cc MCLQMN tests CLQMN.
c
c  Modified:
c
c    24 March 2012
c
c
c       Purpose: This program computes the associated Legendre 
c                functions Qmn(z) and their derivatives Qmn'(z) for
c                a complex argument using subroutine CLQMN
c       Definition: Qmn(z)=(-1)**m*(1-z*z)**(m/2)*dm/dzm[Qn(z)]
c                   Q0(z)=1/2*LOG[(1+z)/(1-z)]     ( for |z|<1 )
c                   Qmn(z)=(z*z-1)**(m/2)*dm/dzm[Qn(z)]
c                   Q0(z)=1/2*LOG[(z+1)/(z-1)]     ( for |z|>1 )
c       Input :  x --- Real part of z
c                y --- Imaginary part of z
c                m --- Order of Qmn(z)  ( m = 0,1,2,úúú )
c                n --- Degree of Qmn(z) ( n = 0,1,2,úúú )
c       Output:  CQM(m,n) --- Qmn(z)
c                CQD(m,n) --- Qmn'(z)
c       Examples:
c                n = 5, x = 0.5, y = 0.2
c
c       m     Re[Qmn(z)]    Im[Qmn(z)]    Re[Qmn'(z)]   Im[Qmn'(z)]
c      -------------------------------------------------------------
c       0    .987156D+00   .354345D+00    .324023D+01  -.447297D+01
c       1   -.240328D+01   .436861D+01    .281158D+02   .171437D+02
c       2   -.245853D+02  -.138072D+02   -.106283D+03   .913792D+02
c       3    .102723D+03  -.651233D+02   -.362578D+03  -.429802D+03
c       4    .155510D+03   .357712D+03    .196975D+04  -.287414D+02
c       5   -.167357D+04  -.680954D+03   -.193093D+04  -.925757D+03
c
c                n = 5, x = 2.5, y = 1.0
c
c       m     Re[Qmn(z)]    Im[Qmn(z)]    Re[Qmn'(z)]   Im[Qmn'(z)]
c      -------------------------------------------------------------
c       0   -.274023D-04  -.227141D-04    .809834D-04   .210884D-04
c       1    .165620D-03   .136108D-03   -.489095D-03  -.124400D-03
c       2   -.118481D-02  -.948832D-03    .349090D-02   .825057D-03
c       3    .982179D-02   .753264D-02   -.288271D-01  -.596384D-02
c       4   -.927915D-01  -.669521D-01    .270840D+00   .451376D-01
c       5    .985601D+00   .656737D+00   -.285567D+01  -.332533D+00
c
c
        IMPLICIT DOUBLE PRECISION (X,Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CQM(0:40,0:40),CQD(0:40,0:40)

      integer test_num
      parameter ( test_num = 12 )

      integer test
      integer m
      integer m_test(test_num)
      integer n
      integer n_test(test_num)
      double precision x
      double precision x_test(test_num)
      double precision y
      double precision y_test(test_num)

      save m_test
      save n_test
      save x_test
      save y_test

      data m_test / 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5 /
      data n_test / 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5 /
      data x_test / 
     &  0.5D+00, 0.5D+00, 0.5D+00, 0.5D+00, 0.5D+00, 0.5D+00,
     &  2.5D+00, 2.5D+00, 2.5D+00, 2.5D+00, 2.5D+00, 2.5D+00 /
      data y_test / 
     &  0.2D+00, 0.2D+00, 0.2D+00, 0.2D+00, 0.2D+00, 0.2D+00,
     &  1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MCLQMN'
      write ( *, '(a)' ) '  Test CLQMN'
      write ( *, '(a)' ) ' '

      do test = 1, test_num

        m = m_test(test)
        n = n_test(test)
        x = x_test(test)
        y = y_test(test)

        WRITE(*,30)M,N,X,Y
        CALL CLQMN(40,M,N,X,Y,CQM,CQD)
        WRITE(*,*)'   m   n   Re[Qmn(z)]    Im[Qmn(z)]    ',
     &            'Re[Qmn''(z)]   Im[Qmn''(z)]'
        WRITE(*,*)' -----------------------------------',
     &            '------------------------------'
        DO 10 J=0,N
10         WRITE(*,20)M,J,CQM(M,J),CQD(M,J)
20      FORMAT(1X,2I4,2D14.6,1X,2D14.6)
30      FORMAT(1X,'m =',I2,', ','n =',I2,', ','x =',F4.1,
     &         ', ','y =',F4.1)

      end do

      return
        END
      subroutine MCLQN ( )

c*********************************************************************72
c
cc MCLQN tests CLQN.
c
c  Modified:
c
c    24 March 2012
c
c
c       Purpose: This program computes the Legendre polynomials 
c                Qn(z) and Qn'(z) for a complex argument using
c              subroutine CLQN
c       Input :  x --- Real part of z
c                y --- Imaginary part of z
c                n --- Degree of Qn(z), n = 0,1,...
c       Output:  CQN(n) --- Qn(z)
c                CQD(n) --- Qn'(z)
c       Examples:
c
c       z = 0.5 + 0.5 i
c       n    Re[Qn(z)]     Im[Qn(z)]     Re[Qn'(z)]    Im[Qn'(z)]
c      -----------------------------------------------------------
c       0   .402359D+00   .553574D+00   .800000D+00   .400000D+00
c       1  -.107561D+01   .477967D+00   .602359D+00   .115357D+01
c       2  -.136636D+01  -.725018D+00  -.242682D+01   .183390D+01
c       3   .182619D+00  -.206146D+01  -.622944D+01  -.247151D+01
c       4   .298834D+01  -.110022D+01  -.114849D+01  -.125963D+02
c       5   .353361D+01   .334847D+01   .206656D+02  -.123735D+02
c
c       z = 3.0 + 2.0 i
c       n    Re[Qn(z)]     Im[Qn(z)]     Re[Qn'(z)]    Im[Qn'(z)]
c      -----------------------------------------------------------
c       0   .229073D+00  -.160875D+00  -.250000D-01   .750000D-01
c       1   .896860D-02  -.244805D-01   .407268D-02   .141247D-01
c       2  -.736230D-03  -.281865D-02   .190581D-02   .155860D-02
c       3  -.264727D-03  -.227023D-03   .391535D-03   .314880D-04
c       4  -.430648D-04  -.443187D-05   .527190D-04  -.305592D-04
c       5  -.481362D-05   .265297D-05   .395108D-05  -.839883D-05
c
c
        IMPLICIT DOUBLE PRECISION (X,Y)
        IMPLICIT COMPLEX*16 (C,Z)

      integer test_num
      parameter ( test_num = 2 )

        DIMENSION CQN(0:100),CQD(0:100)

      integer test
      integer n
      integer n_test(test_num)
      double precision x
      double precision x_test(test_num)
      double precision y
      double precision y_test(test_num)

      save n_test
      save x_test
      save y_test

      data n_test / 5, 5 /
      data x_test / 0.5D+00, 3.0D+00 /
      data y_test / 0.5D+00, 2.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MCLQN'
      write ( *, '(a)' ) '  Test CLQN'
      write ( *, '(a)' ) ' '

      do test = 1, test_num

        n = n_test(test)
        x = x_test(test)
        y = y_test(test)

        WRITE(*,30)X,Y
        WRITE(*,*)
        CALL CLQN(N,X,Y,CQN,CQD)
        WRITE(*,*)'  n    Re[Qn(z)]     Im[Qn(z)]     Re[Qn''(z)]',
     &            '    Im[Qn''(z)]'
        WRITE(*,*)' ---------------------------------------------',
     &            '--------------'
        DO 10 K=0,N
10         WRITE(*,20)K,CQN(K),CQD(K)

      end do

20      FORMAT(1X,I3,4D14.6)
30      FORMAT(3X,'x =',F5.1,',  ','y =',F5.1)

      return
        END
      subroutine MCOMELP ( )

c*********************************************************************72
c
cc MCOMELP tests COMELP.
c
c  Modified:
c
c    24 March 2012
c
c
c       Purpose: This program computes complete elliptic 
c                integrals K(k) and E(k) using subroutine 
c                COMELP
c       Input  : K  --- Modulus k ( 0 ó k ó 1 )
c       Output : CK --- K(k)
c                CE --- E(k)
c       Example:
c                  k         K(k)          E(K)
c                ---------------------------------
c                 .00      1.570796      1.570796
c                 .25      1.596242      1.545957
c                 .50      1.685750      1.467462
c                 .75      1.910990      1.318472
c                1.00       ì            1.000000
c
c
      integer test_num
      parameter ( test_num = 5 )

        DOUBLE PRECISION HK,CK,CE

      double precision hk_test(test_num)
      integer test

      save hk_test

      data hk_test / 0.00D+00, 0.25D+00, 0.50D+00, 0.75D+00, 1.00D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MCOMELP'
      write ( *, '(a)' ) '  Test COMELP'
      write ( *, '(a)' ) ' '
        WRITE(*,*)'    k         K(k)          E(K)'
        WRITE(*,*)'  ---------------------------------'

      do test = 1, test_num
        hk = hk_test(test)

        CALL COMELP(HK,CK,CE)
        IF (HK.NE.1.0) WRITE(*,10) HK,CK,CE
        IF (HK.EQ.1.0) WRITE(*,20) HK,CE
      end do

10      FORMAT(2X,F5.2,2F14.6)
20      FORMAT(2X,F5.2,7X,'ì',6X,F14.6)

      return
        END
      subroutine MCPBDN ( )

c*********************************************************************72
c
cc MCPBDN tests CPBDN.
c
c  Modified:
c
c    25 March 2012
c
c
c       Purpose: This program computes parabolic cylinder functions 
c                Dn(z) for an integer order and a complex argument
c                using subroutine CPBDN
c       Input :  x --- Real part of z
c                y --- Imaginary part of z
c                n --- Order of Dn(z)
c       Output:  CPB(|n|) --- Dn(z)
c                CPD(|n|) --- Dn'(z)
c       Example:
c                z = 5.0+ 5.0 i
c
c     n     Re[Dn(z)]      Im[Dn(z)]      Re[Dn'(z)]     Im[Dn'(z)]
c   -----------------------------------------------------------------
c     0   .99779828D+00  .66321897D-01 -.23286910D+01 -.26603004D+01
c     1   .46573819D+01  .53206009D+01  .26558457D+01 -.24878635D+02
c     2  -.43138931D+01  .49823592D+02  .14465848D+03 -.10313305D+03
c     3  -.28000219D+03  .21690729D+03  .12293320D+04  .30720802D+03
c     4  -.24716057D+04 -.46494526D+03  .38966424D+04  .82090067D+04
c    -1   .10813809D+00 -.90921592D-01 -.50014908D+00 -.23280660D-01
c    -2   .24998820D-02 -.19760577D-01 -.52486940D-01  .47769856D-01
c    -3  -.15821033D-02 -.23090595D-02 -.68249161D-03  .10032670D-01
c    -4  -.37829961D-03 -.10158757D-03  .89032322D-03  .11093416D-02
c
c
      IMPLICIT DOUBLE PRECISION (X,Y)
      IMPLICIT COMPLEX*16 (C,Z)
      DIMENSION CPB(0:100),CPD(0:100)

      integer j
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MCPBDN'
      write ( *, '(a)' ) '  Test CPBDN'

      x = 5.0D+00
      y = 5.0D+00

      do j = 1, 2

        if ( j == 1 ) then
          n = 5
        else 
          n = - 5
        end if

      write ( *, '(a)' ) ' '
      WRITE(*,20)N,X,Y
      Z=CMPLX(X,Y)
      N0=ABS(N)
      CALL CPBDN(N,Z,CPB,CPD)
      WRITE(*,*)
      IF (N.GE.0) THEN
       WRITE(*,*)'  n     Re[Dn(z)]       Im[Dn(z)]       ',
     &               'Re[Dn''(z)]      Im[Dn''(z)]'
      ELSE
       WRITE(*,*)' -n     Re[Dn(z)]       Im[Dn(z)]       ',
     &               'Re[Dn''(z)]      Im[Dn''(z)]'
      ENDIF
      WRITE(*,*)'-------------------------------------------',
     &            '-------------------------'
      DO 5 I=0,N0
5          WRITE(*,10)I,CPB(I),CPD(I)
10      FORMAT(1X,I3,4D16.8)
20      FORMAT(1X,'N =',I3,',   z =x+iy :',F6.2,'+',F6.2,' i')
      end do

      return
      END
      subroutine MCPSI ( )

c*********************************************************************72
c
cc MCPSI tests CPSI.
c
c  Modified:
c
c    25 March 2012
c
c
c       Purpose: This program computes the psi function psi(z)
c                for a complex argument using subroutine CPSI
c       Input :  x   --- Real part of z
c                y   --- Imaginary part of z
c       Output:  PSR --- Real part of psi(z)
c                PSI --- Imaginary part of psi(z)
c       Examples:
c                   x       y      Re[psi(z)]     Im[psi(z)]
c                 -------------------------------------------
c                  3.0     2.0     1.16459152      .67080728
c                  3.0    -2.0     1.16459152     -.67080728
c                 -3.0     2.0     1.39536075     2.62465344
c                 -3.0    -2.0     1.39536075    -2.62465344
c
c
        DOUBLE PRECISION PSR,PSI

      integer test_num
      parameter ( test_num = 4 )

      integer test
      double precision x
      double precision x_test(test_num)
      double precision y
      double precision y_test(test_num)

      save x_test
      save y_test

      data x_test / 3.0D+00, 3.0D+00, -3.0D+00, -3.0D+00 /
      data y_test / 2.0D+00, -2.0D+00, 2.0D+00, -2.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MCPSI'
      write ( *, '(a)' ) '  Test CPSI'
      write ( *, '(a)' ) ' '
        WRITE(*,*)'   x       y      Re[Psi(z)]      Im[Psi(z)]'
        WRITE(*,*)' ----------------------------------------------'

      do test = 1, test_num

        x = x_test(test)
        y = y_test(test)

        CALL CPSI(X,Y,PSR,PSI)
        WRITE(*,10)X,Y,PSR,PSI
10      FORMAT(1X,F5.1,3X,F5.1,2X,2E16.8)
20      FORMAT(1X,2Hx=,F6.2,6X,2Hy=,F6.2)

      end do

      return
        END
      subroutine MCSPHIK ( )

c*********************************************************************72
c
cc MCSPHIK tests CSPHIK.
c
c  Modified:
c
c    26 March 2012
c
c
c       Purpose: This program computes the modified spherical Bessel 
c                functions and their derivatives for a complex
c                argument using subroutine CSPHIK
c       Input :  z --- Complex argument
c                n --- Order of in(z) & kn(z) ( 0 ó n ó 250 )
c       Output:  CSI(n) --- in(z)
c                CDI(n) --- in'(z)
c                CSK(n) --- kn(z)
c                CDK(n) --- kn'(z)
c       Example: z =4.0+i 2.0
c
c     n     Re[in(z)]      Im[in(z)]     Re[in'(z)]     Im[in'(z)]
c    ---------------------------------------------------------------
c     0   .2118080D+00   .6101922D+01  -.4439356D+00   .4900150D+01
c     1  -.4439356D+00   .4900150D+01  -.5906477D+00   .4053075D+01
c     2  -.9918756D+00   .3028652D+01  -.7574058D+00   .2785396D+01
c     3  -.9663859D+00   .1375561D+01  -.7689911D+00   .1541649D+01
c     4  -.6018277D+00   .4263967D+00  -.5777565D+00   .6482500D+00
c     5  -.2668530D+00   .6640148D-01  -.3214450D+00   .1866032D+00
c
c     n     Re[kn(z)]      Im[kn(z)]     Re[kn'(z)]     Im[kn'(z)]
c    ---------------------------------------------------------------
c     0  -.5010582D-02  -.4034862D-02   .6416184D-02   .4340777D-02
c     1  -.6416184D-02  -.4340777D-02   .8445211D-02   .4487936D-02
c     2  -.1016253D-01  -.4714473D-02   .1392804D-01   .4120703D-02
c     3  -.1893595D-01  -.3973987D-02   .2690088D-01   .3192843D-03
c     4  -.3945464D-01   .2977107D-02   .5690203D-01  -.1873044D-01
c     5  -.8727490D-01   .3689398D-01   .1220481D+00  -.9961483D-01
c
c
        IMPLICIT COMPLEX*16 (C,Z)
        DOUBLE PRECISION X,Y
        DIMENSION CSI(0:250),CDI(0:250),CSK(0:250),CDK(0:250)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MCSPHIK'
      write ( *, '(a)' ) '  Test CSPHIK'
      write ( *, '(a)' ) ' '

      n = 5
      x = 4.0D+00
      y = 2.0D+00
      ns = 1

        WRITE(*,30)N,X,Y
        Z=CMPLX(X,Y)
        IF (N.LE.8) THEN
           NS=1
        ELSE
           ns = 5
        ENDIF
        CALL CSPHIK(N,Z,NM,CSI,CDI,CSK,CDK)
        WRITE(*,*)
        WRITE(*,*)'  n      Re[in(z)]        Im[in(z)]',
     &  '        Re[in''(z)]       Im[in''(z)]'
        WRITE(*,*)'--------------------------------------------',
     &  '----------------------------'
        DO 10 K=0,NM,NS
10         WRITE(*,20)K,CSI(K),CDI(K)
        WRITE(*,*)
        WRITE(*,*)'  n      Re[kn(z)]        Im[kn(z)]',
     &  '        Re[kn''(z)]       Im[kn''(z)]'
        WRITE(*,*)'--------------------------------------------',
     &  '----------------------------'
        DO 15 K=0,NM,NS
15         WRITE(*,20)K,CSK(K),CDK(K)
20      FORMAT(1X,I3,4D17.8)
30      FORMAT(3X,'Nmaz =',I3,',     ','z = ',F8.1,'+ i',F8.1)

      return
        END
      subroutine MCSPHJY ( )

c*********************************************************************72
c
cc MCSPHJY tests CSPHJY.
c
c  Modified:
c
c    26 March 2012
c
c
c       Purpose: This program computes the spherical Bessel functions 
c                jn(z), yn(z), and their derivatives for a complex
c                argument using subroutine CSPHJY
c       Input :  z --- Complex argument
c                n --- Order of jn(z) & yn(z) ( 0 ó n ó 250 )
c       Output:  CSJ(n) --- jn(z)
c                CDJ(n) --- jn'(z)
c                CSY(n) --- yn(z)
c                CDY(n) --- yn'(z)
c       Example: z = 4.0+i 2.0
c
c     n     Re[jn(z)]       Im[jn(z)]       Re[jn'(z)]      Im[jn'(z)]
c   --------------------------------------------------------------------
c     0  -.80651523D+00  -.18941093D+00  -.37101203D-01   .75210758D+00
c     1   .37101203D-01  -.75210758D+00  -.67093420D+00   .11885235D+00
c     2   .60314368D+00  -.27298399D+00  -.24288981D+00  -.40737409D+00
c     3   .42955048D+00   .17755176D+00   .18848259D+00  -.24320520D+00
c     4   .12251323D+00   .22087111D+00   .19660170D+00   .17937264D-01
c     5  -.10242676D-01   .10975433D+00   .68951842D-01   .83020305D-01
c
c     n     Re[yn(z)]       Im[yn(z)]       Re[yn'(z)]      Im[yn'(z)]
c   --------------------------------------------------------------------
c     0   .21734534D+00  -.79487692D+00  -.77049661D+00  -.87010064D-02
c     1   .77049661D+00   .87010064D-02  -.92593503D-01  -.64425800D+00
c     2   .24756293D+00   .56894854D+00   .45127429D+00  -.25839924D+00
c     3  -.23845941D+00   .43646607D+00   .26374403D+00   .12439192D+00
c     4  -.27587985D+00   .20902555D+00  -.67092335D-01   .89500599D-01
c     5  -.70001327D-01   .18807178D+00  -.30472133D+00  -.58661384D-01
c
c
        IMPLICIT COMPLEX*16 (C,Z)
        DOUBLE PRECISION X,Y
        DIMENSION CSJ(0:250),CDJ(0:250),CSY(0:250),CDY(0:250)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MSPHJY'
      write ( *, '(a)' ) '  Test CSPHJY'

      n = 5
      x = 4.0D+00
      y = 2.0D+00
      ns = 1

        Z=CMPLX(X,Y)
        IF (N.LE.8) THEN
           NS=1
        ELSE
           ns = 5
        ENDIF
        CALL CSPHJY(N,Z,NM,CSJ,CDJ,CSY,CDY)
        WRITE(*,*)
        WRITE(*,*)'  n      Re[jn(z)]        Im[jn(z)]',
     &  '        Re[jn''(z)]       Im[jn''(z)]'
        WRITE(*,*)'--------------------------------------------',
     &  '----------------------------'
        DO 10 K=0,NM,NS
10         WRITE(*,20)K,CSJ(K),CDJ(K)
        WRITE(*,*)
        WRITE(*,*)'  n      Re[yn(z)]        Im[yn(z)]',
     &  '        Re[yn''(z)]       Im[yn''(z)]'
        WRITE(*,*)'--------------------------------------------',
     &  '----------------------------'
        DO 15 K=0,NM,NS
15         WRITE(*,20)K,CSY(K),CDY(K)
20      FORMAT(1X,I3,4D17.8)
30      FORMAT(3X,6HNmaz =,I3,',     ','z = ',F8.1,'+ i',F8.1)
      return
        END
      subroutine MCVA1 ( )

c*********************************************************************72
c
cc MCVA1 tests CVA1.
c
c  Modified:
c
c    27 March 2012
c
c
c       Purpose: This program computes a sequence of characteristic 
c                values of Mathieu functions using subroutine CVA1
c       Input :  m  --- Order of Mathieu functions
c                q  --- Parameter of Mathieu functions
c                KD --- Case code
c                       KD=1 for cem(x,q)  ( m = 0,2,4,...)
c                       KD=2 for cem(x,q)  ( m = 1,3,5,...)
c                       KD=3 for sem(x,q)  ( m = 1,3,5,...)
c                       KD=4 for sem(x,q)  ( m = 2,4,6,...)
c       Output:  CV(I) --- Characteristic values; I = 1,2,3,...
c                For KD=1, CV(1), CV(2), CV(3),..., correspond to
c                the characteristic values of cem for m = 0,2,4,...
c                For KD=2, CV(1), CV(2), CV(3),..., correspond to
c                the characteristic values of cem for m = 1,3,5,...
c                For KD=3, CV(1), CV(2), CV(3),..., correspond to
c                the characteristic values of sem for m = 1,3,5,...
c                For KD=4, CV(1), CV(2), CV(3),..., correspond to
c                the characteristic values of sem for m = 0,2,4,...
c
c       Example: Mmax = 12,    q = 25.00
c
c                Characteristic values of Mathieu functions
c
c                  m            a                  b
c                ------------------------------------------
c                  0      -40.256779547
c                  1      -21.314899691      -40.256778985
c                  2       -3.522164727      -21.314860622
c                  3       12.964079444       -3.520941527
c                  4       27.805240581       12.986489953
c                  5       40.050190986       28.062765899
c                  6       48.975786716       41.801071292
c                  7       57.534689001       55.002957151
c                  8       69.524065166       69.057988351
c                  9       85.076999882       85.023356505
c                 10      103.230204804      103.225680042
c                 11      123.643012376      123.642713667
c                 12      146.207690643      146.207674647
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION CV1(200),CV2(200),CVE(200),CVS(200)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MCVA1'
      write ( *, '(a)' ) '  Test CVA1'
      write ( *, '(a)' ) ' '

      mmax = 12
      q = 25.0D+00

        WRITE(*,25)MMAX,Q
        WRITE(*,*)
        CALL CVA1(1,MMAX,Q,CV1)
        CALL CVA1(2,MMAX,Q,CV2)
        DO 10 J=1,MMAX/2+1
           CVE(2*J-1)=CV1(J)
10         CVE(2*J)=CV2(J)
        CALL CVA1(3,MMAX,Q,CV1)
        CALL CVA1(4,MMAX,Q,CV2)
        DO 15 J=1,MMAX/2+1
           CVS(2*J)=CV1(J)
15         CVS(2*J+1)=CV2(J)
        WRITE(*,35)
        WRITE(*,*)
        WRITE(*,*)'  m            a                  b'
        WRITE(*,*)'------------------------------------------'
        DO 20 J=0,MMAX
           IF (J.EQ.0) WRITE(*,30)J,CVE(J+1)
           IF (J.NE.0) WRITE(*,30)J,CVE(J+1),CVS(J+1)
20      CONTINUE
25      FORMAT(3X,6HMmax =,I3,',    ',3Hq =,F6.2)
30      FORMAT(1X,I3,2F19.9)
35      FORMAT(1X,'Characteristic values of Mathieu functions')
      return
        END
      subroutine MCVA2 ( )

c*********************************************************************72
c
cc MCVA2 tests CVA2.
c
c  Modified:
c
c    27 March 2012
c
c
c       Purpose: This program calculates a specific characteristic 
c                value of Mathieu functions using subroutine CVA2
c       Input :  m  --- Order of Mathieu functions
c                q  --- Parameter of Mathieu functions
c                KD --- Case code
c                       KD=1 for cem(x,q)  ( m = 0,2,4,...)
c                       KD=2 for cem(x,q)  ( m = 1,3,5,...)
c                       KD=3 for sem(x,q)  ( m = 1,3,5,...)
c                       KD=4 for sem(x,q)  ( m = 2,4,6,...)
c       Output:  A  --- Characteristic value
c       Example: q = 25.0, m = 0,1,2,...,12
c
c                Characteristic values of Mathieu functions
c
c                  m            a                  b
c                ------------------------------------------
c                  0      -40.256779547
c                  1      -21.314899691      -40.256778985
c                  2       -3.522164727      -21.314860622
c                  3       12.964079444       -3.520941527
c                  4       27.805240581       12.986489953
c                  5       40.050190986       28.062765899
c                  6       48.975786716       41.801071292
c                  7       57.534689001       55.002957151
c                  8       69.524065166       69.057988351
c                  9       85.076999882       85.023356505
c                 10      103.230204804      103.225680042
c                 11      123.643012376      123.642713667
c                 12      146.207690643      146.207674647
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MCVA2'
      write ( *, '(a)' ) '  Test CVA2'
      write ( *, '(a)' ) ' '

      q = 25.0D+00

        WRITE(*,*)'  m            a                  b'
        WRITE(*,*)'-----------------------------------'

      do m = 0, 12

        if ( m == 0 ) then
          kd = 1
          CALL CVA2(KD,M,Q,A)
          WRITE(*,10)M,A
        else if ( mod ( m, 2 ) == 0 ) then
          kd = 1
          CALL CVA2(KD,M,Q,A)
          kd = 3
          CALL CVA2(KD,M,Q,b)
          WRITE(*,10)M,A,b
        else
          kd = 2
          CALL CVA2(KD,M,Q,A)
          kd = 4
          CALL CVA2(KD,M,Q,b)
          WRITE(*,10)M,A,b
        end if

      end do

10      FORMAT(1X,I3,F18.8,F18.8)
20      FORMAT(1X,'Characteristic value of even Mathieu function, a')
30      FORMAT(1X,'Characteristic value of odd Mathieu function, b')
      return
        END
      subroutine MCYZO ( )

c*********************************************************************72
c
cc MCYZO tests CYZO.
c
c  Modified:
c
c    28 March 2012
c
c
c       Purpose : This program evaluates the complex zeros of 
c                 Y0(z), Y0'(z), Y1(z) and Y1'(z), and their 
c                 associated values at the zeros using the
c                 modified Newton's iteration method
c       Input:    NT --- Total number of roots/zeros
c                 KF --- Function choice code
c                        KF=0 for  Y0(z) & Y1(z0)
c                        KF=1 for  Y1(z) & Y0(z1)
c                        KF=2 for  Y1'(z) & Y1(z1')
c                 KC --- Choice code
c                        KC=0 for complex roots
c                        KC=1 for real roots
c       Output:   ZO(L) --- L-th zero of Y0(z) or Y1(z) or Y1'(z)
c                 ZV(L) --- Value of Y0'(z) or Y1'(z) or Y1(z)
c                           at the L-th zero
c       Examples: NT = 5
c
c   No.      z0, Zeros of Y0(z)                Y1(z0)
c  -----------------------------------------------------------------
c    1   -2.403016632 + i .5398823130   .1007476893 - i .8819677101
c    2   -5.519876702 + i .5471800106  -.0292464182 + i .5871695027
c    3   -8.653672403 + i .5484120673   .0149080637 - i .4694587524
c    4  -11.791512030 + i .5488191184  -.0093736817 + i .4023045429
c    5  -14.930906564 + i .5490008289   .0065788031 - i .3575673214
c
c   No.      z1, Zeros of Y1(z)                 Y0(z1)
c  -----------------------------------------------------------------
c    1    -.502743273 + i .7862437145  -.4595276847 + i1.3171019361
c    2   -3.833535193 + i .5623565382   .0483019087 - i .6925128842
c    3   -7.015903683 + i .5533930459  -.0201269494 + i .5186425332
c    4  -10.173573834 + i .5512733877   .0116140017 - i .4320329636
c    5  -13.323739307 + i .5504585830  -.0077719300 + i .3779698048
c
c   No.      z1', Zeros of Y1'(z)                Y1(z1')
c   ----------------------------------------------------------------
c    1     .576785129 + i .9039847922  -.7634970879 + i .5892448647
c    2   -1.940477342 - i .7211859189   .1620640057 + i .9520278864
c    3   -5.333478617 - i .5672196368  -.0317940081 - i .5968536736
c    4   -8.536768577 - i .5560607040   .0154177166 + i .4726011652
c    5  -11.706175219 - i .5528590607  -.0095443768 - i .4037533396
c
c
        IMPLICIT COMPLEX *16 (Z)
        DIMENSION ZO(50),ZV(50)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MCYZO'
      write ( *, '(a)' ) '  Test CYZO'

      nt = 5
      do kc = 0, 1
        do kf = 0, 2

        WRITE(*,*)
        WRITE(*,20)NT,KF,KC
        CALL CYZO(NT,KF,KC,ZO,ZV)

        IF (KF.EQ.0) THEN
           WRITE(*,*)' No.          z0, Zeros of Y0(z)',
     &               '                 Y1(z0)'
        ELSE IF (KF.EQ.1) THEN
           WRITE(*,*)' No.          z1, Zeros of Y1(z)',
     &               '                 Y0(z1)'
        ELSE IF (KF.EQ.2) THEN
           WRITE(*,*)' No.        z1'', Zeros of Y1''(z)',
     &               '                Y1(z1'')'
        ENDIF
        WRITE(*,*)'--------------------------------------',
     &            '----------------------------'
        DO 10 I=1,NT
10         WRITE(*,25)I,ZO(I),ZV(I)

        end do
      end do

15      FORMAT(20X,'*****    Please Wait !    *****')
20      FORMAT(2X,'NT=',I3,',  ','KF=',I3,',  ','KC=',I3)
25      FORMAT(1X,I3,2X,F15.9,3F15.10)
      return
        END
      subroutine ME1XA ( )

c*********************************************************************72
c
cc ME1XA tests E1XA.
c
c  Modified:
c
c    28 March 2012
c
c
c       Purpose: This program computes the exponential integral 
c                E1(x) using subroutine E1XA
c       Input :  x  --- Argument of E1(x)  ( x > 0 )
c       Output:  E1 --- E1(x)
c       Example:
c                  x        E1(x)
c                ----------------------
c                 0.0     .1000000+301
c                 1.0     .2193839E+00
c                 2.0     .4890051E-01
c                 3.0     .1304838E-01
c                 4.0     .3779352E-02
c                 5.0     .1148296E-02
c
c
        DOUBLE PRECISION E1,X

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ME1XA'
      write ( *, '(a)' ) '  Test E1XA'
      write ( *, '(a)' ) ' '

        WRITE(*,*)'   x        E1(x)'
        WRITE(*,*)' ----------------------'
        do i = 0, 5
          x = real ( i )
          CALL E1XA(X,E1)
          WRITE(*,10)X,E1
        end do

10      FORMAT(1X,F5.1,E17.7)

      return
        END
      subroutine ME1XB ( )

c*********************************************************************72
c
cc ME1XB tests E1XB.
c
c  Modified:
c
c    29 March 2012
c
c
c       Purpose: This program computes the exponential integral 
c                E1(x) using subroutine E1XB
c       Input :  x  --- Argument of E1(x)  ( x > 0 )
c       Output:  E1 --- E1(x)
c       Example:
c                  x          E1(x)
c                -------------------------
c                 0.0     .1000000000+301
c                 1.0     .2193839344E+00
c                 2.0     .4890051071E-01
c                 3.0     .1304838109E-01
c                 4.0     .3779352410E-02
c                 5.0     .1148295591E-02
c
c
        DOUBLE PRECISION E1,X

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ME1XB'
      write ( *, '(a)' ) '  Test E1XB'
      write ( *, '(a)' ) ' '

        WRITE(*,*)'   x          E1(x)'
        WRITE(*,*)' -------------------------'
        do i = 0, 5
          x = real ( i )
        CALL E1XB(X,E1)
        WRITE(*,10)X,E1
      end do

10      FORMAT(1X,F5.1,E20.10)
      return
        END
      subroutine ME1Z ( )

c*********************************************************************72
c
cc ME1Z tests E1Z.
c
c  Modified:
c
c    29 March 2012
c
c
c       Purpose: This program computes the complex exponential 
c                integral E1(z) using subroutine E1Z
c       Example:
c                     z            Re[E1(z)]       Im[E1(z)]
c                -----------------------------------------------
c                 3.0    2.0    -.90959209D-02   -.69001793D-02
c                 3.0   -2.0    -.90959209D-02    .69001793D-02
c                -3.0    2.0    -.28074891D+01    .59603353D+01
c                -3.0   -2.0    -.28074891D+01   -.59603353D+01
c                25.0   10.0    -.29302080D-12    .40391222D-12
c                25.0  -10.0    -.29302080D-12   -.40391222D-12
c               -25.0   10.0     .27279957D+10   -.49430610D+09
c               -25.0  -10.0     .27279957D+10    .49430610D+09
c
c
        IMPLICIT COMPLEX*16 (C,Z)
        IMPLICIT DOUBLE PRECISION (D-H,O-Y)

      integer test_num
      parameter ( test_num = 8 )

      integer test
      real x
      real x_test(test_num)
      real y
      real y_test(test_num)

      save x_test
      save y_test

      data x_test /
     &   3.0D+00,   3.0D+00,  -3.0D+00,  -3.0D+00,
     &  25.0D+00,  25.0D+00, -25.0D+00, -25.0D+00 /
      data y_test /
     &   2.0D+00,  -2.0D+00,   2.0D+00,  -2.0D+00,
     &  10.0D+00, -10.0D+00,  10.0D+00, -10.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ME1Z'
      write ( *, '(a)' ) '  Test E1Z'
      write ( *, '(a)' ) ' '
        WRITE(*,*)'       z           Re[E1(z)]        Im[E1(z)]'
        WRITE(*,*)' -----------------------------------------------'

      do test = 1, test_num

        x = x_test(test)
        y = y_test(test)

        Z=CMPLX(X,Y)
        CALL E1Z(Z,CE1)

        WRITE(*,10)X,Y,CE1

      end do

10      FORMAT(1X,F5.1,2X,F5.1,1X,2D17.8)

      return
        END
      subroutine MEIX ( )

c*********************************************************************72
c
cc MEIX tests EIX.
c
c  Modified:
c
c    29 March 2012
c
c
c       Purpose: This program computes the exponential integral 
c                Ei(x) using subroutine EIX
c       Example:
c                  x        Ei(x)
c                -----------------------
c                  0    -.10000000+301
c                  1     .18951178E+01
c                  2     .49542344E+01
c                  3     .99338326E+01
c                  4     .19630874E+02
c                  5     .40185275E+02
c
c
        DOUBLE PRECISION EI,X

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MEIX'
      write ( *, '(a)' ) '  Test EIX'
      write ( *, '(a)' ) ' '

        WRITE(*,*)'   x         Ei(x)'
        WRITE(*,*)'------------------------'
      do i = 0, 5
        x = dble ( i )
        CALL EIX(X,EI)
        WRITE(*,10)X,EI
      end do

10      FORMAT(1X,F5.1,E18.8)
      return
        END
      subroutine MELIT ( )

c*********************************************************************72
c
cc MELIT tests ELIT.
c
c  Modified:
c
c    29 March 2012
c
c
c       Purpose: This program computes complete and incomplete 
c                elliptic integrals F(k,phi) and E(k,phi) using
c              subroutine ELIT
c       Input  : HK  --- Modulus k ( 0 ó k ó 1 )
c                Phi --- Argument ( in degrees )
c       Output : FE  --- F(k,phi)
c                EE  --- E(k,phi)
c       Example:
c                k = .5
c
c                 phi     F(k,phi)       E(k,phi)
c                -----------------------------------
c                   0      .00000000      .00000000
c                  15      .26254249      .26106005
c                  30      .52942863      .51788193
c                  45      .80436610      .76719599
c                  60     1.08955067     1.00755556
c                  75     1.38457455     1.23988858
c                  90     1.68575035     1.46746221
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       double precision k

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MELIT'
      write ( *, '(a)' ) '  Test ELIT'
      write ( *, '(a)' ) ' '
        WRITE(*,*)'   k  phi       F(k,phi)       E(k,phi)'
        WRITE(*,*)' ---------------------------------------'

      k = 0.5D+00

      do i = 0, 6
        phi = dble ( i ) * 15.0D+00
        CALL ELIT(K,PHI,FE,EE)
        WRITE(*,10)K,PHI,FE,EE
      end do

10      FORMAT(1X,F4.2,2x,F4.0,2F15.8)

      return
        END
      subroutine MELIT3 ( )

c*********************************************************************72
c
cc MELIT3 tests ELIT3.
c
c  Modified:
c
c    31 March 2012
c
c
c       Purpose: This program computes the elliptic integral of 
c                the third kind using subroutine ELIT3
c       Input :  Phi --- Argument ( in degrees )
c                 k  --- Modulus   ( 0 ó k ó 1 )
c                 c  --- Parameter ( 0 ó c ó 1 )
c       Output:  EL3 ÄÄÄ Value of the elliptic integral of the
c                        third kind
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      double precision k
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MELIT3'
      write ( *, '(a)' ) '  Test ELIT3'
      write ( *, '(a)' ) ' '
      phi = pi / 4.0D+00
      c = 0.75D+00
      k = 0.50D+00

        write(*,'(a,g14.6,a,g14.6,a,g14.6)')
     &  '  Phi = ', PHI, '  K = ', K, '  C = ', C
        CALL ELIT3(PHI,K,C,EL3)
        WRITE(*,10)EL3
10      FORMAT(1X,'EL3=',F12.8)
      return
        END
      subroutine MENXA ( )

c*********************************************************************72
c
cc MENXA tests ENXA.
c
c  Modified:
c
c    31 March 2012
c
c
c       Purpose: This program computes the exponential integral 
c                En(x) using subroutine ENXA
c       Example: x = 10.0
c                   n         En(x)
c                 ----------------------
c                   0     .45399930D-05
c                   1     .41569689D-05
c                   2     .38302405D-05
c                   3     .35487626D-05
c                   4     .33041014D-05
c                   5     .30897289D-05
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION EN(0:100)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MENXA'
      write ( *, '(a)' ) '  Test ENXA'
      n = 5
      x = 10.0D+00
      write ( *, '(a)' ) ' '
        WRITE(*,'(a,g14.6)' ) '  X = ', x
      write ( *, '(a)') ' '
        WRITE(*,*)'   n         En(x)'
        WRITE(*,*)' ----------------------'

        CALL ENXA(N,X,EN)
        DO 10 K=0,N
           WRITE(*,30)K,EN(K)
10      CONTINUE

20      FORMAT(5X,I3,',   ','x=',F5.1)
30      FORMAT(2X,I3,D18.8)
      return
        END
      subroutine MENXB ( )

c*********************************************************************72
c
cc MENXB tests ENXB.
c
c  Modified:
c
c    31 March 2012
c
c
c       Purpose: This program computes the exponential integral 
c                En(x) using subroutine ENXB
c       Example: x = 10.0
c                   n         En(x)
c                 ----------------------
c                   0     .45399930D-05
c                   1     .41569689D-05
c                   2     .38302405D-05
c                   3     .35487626D-05
c                   4     .33041014D-05
c                   5     .30897289D-05
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION EN(0:100)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MENXB'
      write ( *, '(a)' ) '  Test ENXB'

      n = 5
      x = 10.0D+00
      write ( *, '(a)' ) ' '
        WRITE(*,'(a,g14.6)' ) '  X = ', x
        WRITE(*,*)
        WRITE(*,*)'   n         En(x)'
        WRITE(*,*)' ----------------------'
        CALL ENXB(N,X,EN)
        DO 10 K=0,N
           WRITE(*,30)K,EN(K)
10      CONTINUE
20      FORMAT(5X,I3,',   ','x=',F5.1)
30      FORMAT(2X,I3,D18.8)
      return
        END
      subroutine MERROR ( )

c*********************************************************************72
c
cc MERROR tests ERROR.
c
c  Modified:
c
c    31 March 2012
c
c
c       Purpose: This program computes the error function 
c                erf(x) using subroutine ERROR
c       Input:   x   --- Argument of erf(x)
c       Output:  ERR --- erf(x)
c       Example:
c                  x         erf(x)
c                ---------------------
c                 1.0       .84270079
c                 2.0       .99532227
c                 3.0       .99997791
c                 4.0       .99999998
c                 5.0      1.00000000
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MERROR'
      write ( *, '(a)' ) '  Test ERROR'
      write ( *, '(a)' ) ' '

        WRITE(*,*)'   x         erf(x)'
        WRITE(*,*)'---------------------'

      do i = 1, 5
        x = dble ( i )

        CALL ERROR(X,ER)
        WRITE(*,10)X,ER

      end do

10      FORMAT(1X,F5.2,F15.8)
      return
        END
      subroutine MEULERA ( )

c*********************************************************************72
c
cc MEULERA tests EULERA.
c
c  Modified:
c
c    03 April 2012
c
c
c       Purpose: This program computes Euler number En using
c              subroutine EULERA
c       Example: Compute Euler number En for n = 0,2,...,10
c                Computed results:
c
c                   n            En
c                 --------------------------
c                   0     .100000000000D+01
c                   2    -.100000000000D+01
c                   4     .500000000000D+01
c                   6    -.610000000000D+02
c                   8     .138500000000D+04
c                  10    -.505210000000D+05
c
c
        DOUBLE PRECISION E
        DIMENSION E(0:200)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MEULERA'
      write ( *, '(a)' ) '  Test EULERA'
      write ( *, '(a)' ) ' '

      n = 10

        CALL EULERA(N,E)
        WRITE(*,*)'   n            En'
        WRITE(*,*)' --------------------------'
        DO 10 K=0,N,2
10         WRITE(*,20)K,E(K)
20      FORMAT(2X,I3,D22.12)
      return
        END
      subroutine MEULERB ( )

c*********************************************************************72
c
cc MEULERB tests EULERB.
c
c  Modified:
c
c    03 April 2012
c
c
c       Purpose: This program computes Euler number En using
c              subroutine EULERB
c       Example: Compute Euler number En for n = 0,2,...,10
c                Computed results:
c
c                   n            En
c                 --------------------------
c                   0     .100000000000D+01
c                   2    -.100000000000D+01
c                   4     .500000000000D+01
c                   6    -.610000000000D+02
c                   8     .138500000000D+04
c                  10    -.505210000000D+05
c
c
        DOUBLE PRECISION E
        DIMENSION E(0:200)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MEULERB'
      write ( *, '(a)' ) '  Test EULERB'
      write ( *, '(a)' ) ' '

        n = 10
        CALL EULERB(N,E)
        WRITE(*,*)'   n            En'
        WRITE(*,*)' --------------------------'
        DO 10 K=0,N,2
10         WRITE(*,20)K,E(K)
20      FORMAT(2X,I3,D22.12)
      return
        END
      subroutine MFCOEF ( )

c*********************************************************************72
c
cc MFCOEF tests FCOEF.
c
c  Modified:
c
c    03 April 2012
c
c
c       Purpose: This program computes the expansion coefficients 
c                for Mathieu and modified Mathieu functions using
c              subroutine FCOEF
c       Input :  m  --- Order of Mathieu functions
c                q  --- Parameter of Mathieu functions
c                KD --- Case code
c                       KD=1 for cem(x,q)  ( m = 0,2,4,...)
c                       KD=2 for cem(x,q)  ( m = 1,3,5,...)
c                       KD=3 for sem(x,q)  ( m = 1,3,5,...)
c                       KD=4 for sem(x,q)  ( m = 2,4,6,...)
c                A  --- Characteristic value of Mathieu
c                       functions for given m and q
c       Output:  FC(k) --- Expansion coefficients of Mathieu
c                       functions ( k= 1,2,...,KM )
c                       FC(1),FC(2),FC(3),... correspond to
c                       A0,A2,A4,... for KD=1 case, A1,A3,
c                       A5,... for KD=2 case, B1,B3,B5,...
c                       for KD=3 case and B2,B4,B6,... for
c                       KD=4 case
c       Example: Compute the first 12 expansion coefficients of
c                even and odd Mathieu functions for given
c                m = 10 And q= 5.0
c                Expansion coefficients of Mathieu functions
c
c                   n         Amn(q)            Bmn(q)
c                 ----------------------------------------
c                   0     .16788542D-05
c                   2     .33619515D-04     .33444320D-04
c                   4     .64298667D-03     .64297621D-03
c                   6     .10784807D-01     .10784806D-01
c                   8     .13767512D+00     .13767512D+00
c                  10     .98395564D+00     .98395564D+00
c                  12    -.11280678D+00    -.11280678D+00
c                  14     .58929627D-02     .58929627D-02
c                  16    -.18916571D-03    -.18916571D-03
c                  18     .42264064D-05     .42264064D-05
c                  20    -.70485101D-07    -.70485101D-07
c                  22     .91820256D-09     .91820256D-09
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION FC(251),FB(251)

      m = 10
      q = 5.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MFCOEF'
      write ( *, '(a)' ) '  Test FCOEF'
      write ( *, '(a)' ) ' '

        WRITE(*,30)M,Q
        WRITE(*,*)
        IF (M.EQ.2*INT(M/2)) THEN
           KI=1
           KF=4
           KS=3
        ELSE IF (M.NE.2*INT(M/2)) THEN
           KI=2
           KF=3
           KS=1
        ENDIF
        DO 20 KD=KI,KF,KS
           CALL CVA2(KD,M,Q,A)
           CALL FCOEF(KD,M,Q,A,FC)
           IF (KD.EQ.KI) THEN
              DO 10 L=1,M/2+8
                 FB(L)=FC(L)
10            CONTINUE
           ENDIF
20      CONTINUE
        WRITE(*,40)
        WRITE(*,*)
        WRITE(*,45)
        WRITE(*,*)'----------------------------------------'
        DO 25 L=1,M/2+8
           IF (M.EQ.2*INT(M/2)) THEN
              IF (M.EQ.0) WRITE(*,35)2*L-2,FB(L)
              IF (M.GT.0.AND.L.EQ.1) WRITE(*,35)2*L-2,FB(L)
              IF (M.GT.0.AND.L.NE.1) WRITE(*,35)2*L-2,FB(L),FC(L-1)
           ELSE
              WRITE(*,35)2*L-1,FB(L),FC(L)
           ENDIF
25      CONTINUE
30      FORMAT(1X,2Hm=,I3,',    ',2Hq=,F6.2)
35      FORMAT(1X,I3,2D18.8)
40      FORMAT(1X,'Expansion coefficients of Mathieu functions')
45      FORMAT(1X,'  n         Amn(q)        ','    Bmn(q)')
      return
        END
      subroutine MFCS ( )

c*********************************************************************72
c
cc MFCS tests FCS.
c
c  Modified:
c
c    03 April 2012
c
c
c       Purpose: This program computes the Fresnel integrals 
c                C(x) and S(x) using subroutine FCS
c       Input :  x --- Argument of C(x) and S(x)
c       Output:  C --- C(x)
c                S --- S(x)
c       Example:
c                  x          C(x)          S(x)
c                -----------------------------------
c                 0.0      .00000000      .00000000
c                 0.5      .49234423      .06473243
c                 1.0      .77989340      .43825915
c                 1.5      .44526118      .69750496
c                 2.0      .48825341      .34341568
c                 2.5      .45741301      .61918176
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MFCS'
      write ( *, '(a)' ) '  Test FCS'
      write ( *, '(a)' ) ' '

        WRITE(*,*)'   x          C(x)          S(x)'
        WRITE(*,*)' -----------------------------------'
      do i = 0, 5
        x = dble ( i ) / 2.0
        CALL FCS(X,C,S)
        WRITE(*,10)X,C,S
      end do
10      FORMAT(1X,F5.1,2F15.8)

      return
        END
      subroutine MFCSZO ( )

c*********************************************************************72
c
cc MFCSZO tests FCSZO.
c
c  Modified:
c
c    03 April 2012
c
c
c       Purpose : This program computes the complex zeros of the
c                 Fresnel integral C(z) or S(z) using subroutine
c                 FCSZO
c       Input :   KF  --- Function code
c                         KF=1 for C(z) or KF=2 for S(z)
c                 NT  --- Total number of zeros
c       Output:   ZO(L) --- L-th zero of C(z) or S(z)
c       Example:  NT=10
c
c       n     Complex zeros of C(z)        Complex zeros of S(z)
c     ------------------------------------------------------------
c       1    1.7436675 + i .30573506      2.0092570 + i .28854790
c       2    2.6514596 + i .25290396      2.8334772 + i .24428524
c       3    3.3203593 + i .22395346      3.4675331 + i .21849268
c       4    3.8757345 + i .20474747      4.0025782 + i .20085103
c       5    4.3610635 + i .19066973      4.4741893 + i .18768859
c       6    4.7976077 + i .17970801      4.9006784 + i .17732036
c       7    5.1976532 + i .17081930      5.2929467 + i .16884418
c       8    5.5690602 + i .16339854      5.6581068 + i .16172492
c       9    5.9172173 + i .15706585      6.0011034 + i .15562108
c      10    6.2460098 + i .15156826      6.3255396 + i .15030246
c
c
        IMPLICIT DOUBLE PRECISION (E,P,W)
        IMPLICIT COMPLEX *16 (C,Z)
        DIMENSION ZCO(100), ZSO(100)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MFCS'
      write ( *, '(a)' ) '  Test FCS'
      write ( *, '(a)' ) ' '

      nt = 10
      kf = 1
        CALL FCSZO(KF,NT,ZCO)
      kf = 2
      call fcszo ( kf, nt, zso )
        WRITE(*,*)
     &  '  n        Complex zeros of C(z)   Complex zeros of S(z)'
        WRITE(*,*)
     &  '--------------------------------------------------------'
        DO I=1,NT
           WRITE(*,30) I,ZCO(I), zso(i)
        end do

30      FORMAT(1X,I3,F13.8,2X,2H+i,F13.8, F13.8,2X,2H+i,F13.8)
      return
        END
      subroutine MFFK ( )

c*********************************************************************72
c
cc MFFK tests FFK.
c
c  Modified:
c
c    06 April 2012
c
c
c       Purpose: This program computes the modified Fresnel integrals 
c                F+/-(x) and K+/-(x) using subroutine FFK
c       Input :  x   --- Argument of F+/-(x) and K+/-(x)
c                KS  --- Sign code
c                        KS=0 for calculating F+(x) and K+(x)
c                        KS=1 for calculating F_(x) and K_(x)
c       Output:  FR  --- Re[F+/-(x)]
c                FI  --- Im[F+/-(x)]
c                FM  --- |F+/-(x)|
c                FA  --- Arg[F+/-(x)]  (Degs.)
c                GR  --- Re[K+/-(x)]
c                GI  --- Im[K+/-(x)]
c                GM  --- |K+/-(x)|
c                GA  --- Arg[K+/-(x)]  (Degs.)
c       Example:
c
c         x    Re[F+/-(x)] +/-Im[F+/-(x)]  Mod[F+/-(x)] +/-Arg[F+/-(x)]
c       ----------------------------------------------------------
c        0.0    .62665707    .62665707    .88622693    45.000000
c        2.0    .16519561   -.17811942    .24293233   -47.155835
c        4.0    .03219674   -.12047678    .12470479   -75.037684
c        6.0    .08245304   -.01180212    .08329342    -8.145843
c        8.0   -.05729996    .02493542    .06249048   156.482601
c       10.0    .02553188    .04298617    .04999688    59.291561
c
c         x    Re[K+/-(x)] +/-Im[K+/-(x)] Mod[K+/-ñ(x)] +/-Arg[K+/-(x)]
c       ----------------------------------------------------------
c        0.0    .50000000    .00000000    .50000000     0.000000
c        2.0    .10702394    .08562295    .13705989    38.661047
c        4.0    .05126306    .04818949    .07035714    43.229843
c        6.0    .03368650    .03276566    .04699328    44.206095
c        8.0    .02512396    .02473472    .03525648    44.552712
c       10.0    .02004532    .01984592    .02820772    44.713609
c
c
        IMPLICIT DOUBLE PRECISION (F,G,X)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MFFK'
      write ( *, '(a)' ) '  Test FFK'
      write ( *, '(a)' ) ' '

      do i = 0, 5
        x = dble ( 2 * i )
        WRITE(*,*)'   x      Re[F+/-(x)]    +/-Im[F+/-(x)]     ',
     &             'Mod[F+/-(x)]   +/-Arg[F+/-(x)]'
        WRITE(*,*)' ---------------------------------------',
     &            '-----------------------'
        CALL FFK(0,X,FR,FI,FM,FA,GR,GI,GM,GA)
        WRITE(*,10)X,FR,FI,FM,FA
        WRITE(*,*)
        WRITE(*,*)'   x      Re[K+/-(x)]    +/-Im[K+/-(x)]     ',
     &             'Mod[K+/-(x)]   +/-Arg[K+/-(x)]'
        WRITE(*,*)' ---------------------------------------',
     &            '-----------------------'
        WRITE(*,10)X,GR,GI,GM,GA
      end do

10      FORMAT(1X,F5.1,3F14.8,F14.6)
      return
        END
      subroutine MGAMMA ( )

c*********************************************************************72
c
cc MGAMMA tests GAMMA.
c
c  Modified:
c
c    03 April 2012
c
c
c       Purpose: This program computes the gamma function
c                â(x) using subroutine GAMMA
c       Examples:
c                   x            â(x)
c                ----------------------------
c                  1/3       2.678938534708
c                  0.5       1.772453850906
c                 -0.5      -3.544907701811
c                 -1.5       2.363271801207
c                  5.0      24.000000000000
c
c
        DOUBLE PRECISION A,X,GA
        DIMENSION A(5)
        DATA A/.333333333333333333D0,0.5D0,-0.5D0,-1.5,5.0D0/

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MGAMMA'
      write ( *, '(a)' ) '  Test GAMMA'
      write ( *, '(a)' ) ' '

        WRITE(*,*)'     x            â(x)'
        WRITE(*,*)' ----------------------------'
        DO 10 K=1,5
           X=A(K)
           CALL GAMMA(X,GA)
10         WRITE(*,20)X,GA
20      FORMAT(1X,F8.4,E20.12)
      return
        END
      subroutine MHERZO ( )

c*********************************************************************72
c
cc MHERZO tests HERZO.
c
c  Modified:
c
c    04 April 2012
c
c
c       Purpose : This program computes the zeros of Hermite 
c                 polynomial Ln(x) in the interval [-ì,ì] and the
c                 corresponding weighting coefficients for Gauss-
c                 Hermite integration using subroutine HERZO
c       Input :   n    --- Order of the Hermite polynomial
c                 X(n) --- Zeros of the Hermite polynomial
c                 W(n) --- Corresponding weighting coefficients
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION X(100),W(100)
        n = 13

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MHERZO'
      write ( *, '(a)' ) '  Test HERZO'
      write ( *, '(a)' ) ' '

        WRITE(*,20)N
        CALL HERZO(N,X,W)
        WRITE(*,*)'  Nodes and weights for Gauss-Hermite integration'
        WRITE(*,*)
        WRITE(*,*)'  i             xi                      Wi'
        WRITE(*,*)' -----------------------------------------',
     &            '------------'
        DO 10 J=1,N
10         WRITE(*,30)J,X(J),W(J)
20      FORMAT(1X,'n =',I3)
30      FORMAT(1X,I3,3X,D22.13,3X,D22.13)
      return
        END 
      subroutine MHYGFX ( )

c*********************************************************************72
c
cc MHYGFX tests HYGFX.
c
c  Modified:
c
c    04 April 2012
c
c
c       Purpose: This program computes the hypergeometric function 
c                F(a,b,c,x) using subroutine HYGFX
c       Input :  a --- Parameter
c                b --- Parameter
c                c --- Parameter, c <> 0,-1,-2,...
c                x --- Argument ( x ó 1 )
c       Output:  HF --- F(a,b,c,x)
c       Example:
c              b = 3.30,  c = 6.70
c              a     F(a,b,c,.25)     F(a,b,c,.55)    F(a,b,c,.85)
c            ------------------------------------------------------
c            -2.5   .72356129D+00    .46961432D+00   .29106096D+00
c            -0.5   .93610145D+00    .85187390D+00   .75543187D+00
c             0.5   .10689695D+01    .11795358D+01   .13510497D+01
c             2.5   .14051563D+01    .23999063D+01   .57381566D+01
c
c              a = 3.30,  b = 6.70
c              c     F(a,b,c,.25)     F(a,b,c,.55)    F(a,b,c,.85)
c            ------------------------------------------------------
c            -5.5   .15090670D+05    .10170778D+11   .58682088D+19
c            -0.5  -.21631479D+04   -.30854772D+07  -.10217370D+13
c             0.5   .26451677D+03    .11967860D+06   .92370648D+10
c             4.5   .41946916D+01    .58092729D+02   .20396914D+05
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      integer a_test_num
      parameter ( a_test_num = 4 )
      integer c_test_num
      parameter ( c_test_num = 4 )
      integer x_test_num
      parameter ( x_test_num = 3 )

      double precision a
      double precision a_test(a_test_num)
      double precision b
      double precision c
      double precision c_test(c_test_num)
      integer i
      integer j
      integer l
      double precision x
      double precision x_test(x_test_num)

      save a_test
      save c_test
      save x_test

      data a_test / -2.5D+00, -0.5D+00, 0.5D+00, 2.5D+00 /
      data c_test / -5.5D+00, -0.5D+00, 0.5D+00, 4.5D+00 /
      data x_test / 0.25D+00, 0.55D+00, 0.85D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MHYGFX'
      write ( *, '(a)' ) '  Test HYGFX'
      write ( *, '(a)' ) ' '

      do i = 1, a_test_num
        a = a_test(i)
        b = 3.30D+00
        c = 6.70D+00
        do l = 1, x_test_num
          x = x_test(l)
          CALL HYGFX(A,B,C,X,HF)
          WRITE(*,10)A,B,C,X, HF

        end do
      end do

      a = 3.30D+00
      b = 6.70D+00
      do j = 1, c_test_num
        c = c_test(j)
        do l = 1, x_test_num
          x = x_test(l)
          CALL HYGFX(A,B,C,X,HF)
          WRITE(*,10)A,B,C,X, HF
        end do
      end do

10      FORMAT(1X,'a =',F5.2,',    ','b =',F5.2,',   ','c =',F5.2,
     &         ',    ','x =',F5.2,'F(a,b,c,x)=',D16.8)

      return
        END
      subroutine MHYGFZ ( )

c*********************************************************************72
c
cc MHYGFZ tests HYGFZ.
c
c  Modified:
c
c    04 April 2012
c
c
c       Purpose: This program computes hypergeometric function for
c                a complex argument, F(a,b,c,z), using subroutine
c                HYGFZ
c       Input :  a --- Parameter
c                b --- Parameter
c                c --- Parameter,  c <> 0,-1,-2,...
c                x --- Real part of complex argument z
c                y --- Imaginary part of complex argument z
c                      ( z = x+iy )
c       Output:  ZHF --- F(a,b,c,z)
c       Examples:
c     a     b     c    z = x+ iy             F(a,b,c,z)
c   --------------------------------------------------------------
c    3.2   1.8   6.7   1.0+0.0 i    .54689992D+01+.00000000D+00 i
c    3.2  -1.8   6.7   1.0+0.0 i    .33750635D+00+.00000000D+00 i
c   -5.0   3.3   6.7   5.2+4.8 i    .11682745D+03+.60389104D+03 i
c    3.3  -6.0   3.7   5.2-4.8 i    .17620425D+05+.38293812D+05 i
c   -7.0   3.3  -3.7   5.2-4.8 i   -.11772779D+11-.14382286D+11 i
c    4.3  -8.0  -3.7   5.2+4.8 i    .13161188D+13-.10129870D+12 i
c    3.3   5.8   6.7   0.2+0.1 i    .17330557D+01+.63401030D+00 i
c    3.5  -2.4   6.7   0.2+0.5 i    .64762241D+00-.52110507D+00 i
c    3.3   4.3   6.7   0.8+0.3 i   -.14830086D+01+.83744258D+01 i
c    7.0   5.0   4.1   3.0-1.0 i   -.40376095D-02-.29566326D-02 i
c    5.0   7.0   4.1   3.0-1.0 i   -.40376095D-02-.29566326D-02 i
c    3.5   1.2   9.7   0.6+0.9 i    .10343044D+01+.54473814D+00 i
c    2.1   5.4   9.7   0.5+0.7 i    .68850442D+00+.12274187D+01 i
c    8.7   3.2   6.7   0.5+0.7 i   -.90046505D+00-.11198900D+01 i
c    8.7   2.7   6.7   0.6+0.9 i   -.46083890D+00-.54575701D+00 i
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Y)
        IMPLICIT COMPLEX *16 (Z)

      integer test_num
      parameter ( test_num = 15 )

      double precision a
      double precision a_test(test_num)
      double precision b
      double precision b_test(test_num)
      double precision c
      double precision c_test(test_num)
      integer i
      double precision x
      double precision x_test(test_num)
      double precision y
      double precision y_test(test_num)

      save a_test
      save b_test
      save c_test
      save x_test
      save y_test

      data a_test /
     &  3.2D+00, 3.2D+00, -5.0D+00, 3.3D+00, -7.0D+00,
     &  4.3D+00, 3.3D+00,  3.5D+00, 3.3D+00,  7.0D+00, 
     &  5.0D+00, 3.5D+00,  2.1D+00, 8.7D+00,  8.7D+00 /
      data b_test /
     &   1.8D+00, -1.8D+00,  3.3D+00, -6.0D+00, 3.3D+00,
     &  -8.0D+00,  5.8D+00, -2.4D+00,  4.3D+00, 5.0D+00,
     &   7.0D+00,  1.2D+00,  5.4D+00,  3.2D+00, 2.7D+00 /
      data c_test /
     &   6.7D+00, 6.7D+00, 6.7D+00, 3.7D+00, -3.7D+00,
     &  -3.7D+00, 6.7D+00, 6.7D+00, 6.7D+00,  4.1D+00,
     &   4.1D+00, 9.7D+00, 9.7D+00, 6.7D+00,  6.7D+00 /
      data x_test / 
     &  1.0D+00, 1.0D+00, 5.2D+00, 5.2D+00, 5.2D+00,
     &  5.2D+00, 0.2D+00, 0.2D+00, 0.8D+00, 3.0D+00,
     &  3.0D+00, 0.6D+00, 0.5D+00, 0.5D+00, 0.5D+00 /
      data y_test /
     &   0.0D+00, 0.0D+00, 4.8D+00, -4.8D+00, -4.8D+00,
     &   4.8D+00, 0.1D+00, 0.5D+00,  0.3D+00, -1.0D+00,
     &  -1.0D+00, 0.9D+00, 0.7D+00,  0.7D+00,  0.9D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MHYGFZ'
      write ( *, '(a)' ) '  Test HYGFZ'
      write ( *, '(a)' ) ' '

        WRITE(*,*)'     a      b      c      x      y',
     &            '          Re[F]           Im[F]'
        WRITE(*,*)'   --------------------------------',
     &            '-----------------------------------'

      do i = 1, test_num
        a = a_test(i)
        b = b_test(i)
        c = c_test(i)
        x = x_test(i)
        y = y_test(i)
        Z=CMPLX(X,Y)
        CALL HYGFZ(A,B,C,Z,ZHF)

        WRITE(*,10)A,B,C,X,Y,ZHF
      end do

10      FORMAT(1X,5F7.1,2X,2D16.8)
      return
        END
      subroutine MIK01A ( )

c*********************************************************************72
c
cc MIK01A tests IK01A.
c
c  Modified:
c
c    04 April 2012
c
c
c       Purpose: This program computes the modified Bessel functions 
c                I0(x), I1(x), K0(x), K1(x), and their derivatives 
c                using subroutine IK01A
c       Input :  x   --- Argument ( x ò 0 )
c       Output:  BI0 --- I0(x)
c                DI0 --- I0'(x)
c                BI1 --- I1(x)
c                DI1 --- I1'(x)
c                BK0 --- K0(x)
c                DK0 --- K0'(x)
c                BK1 --- K1(x)
c                DK1 --- K1'(x)
c       Example:
c
c         x      I0(x)         I0'(x)        I1(x)         I1'(x)
c       -------------------------------------------------------------
c        1.0  .1266066D+01  .5651591D+00  .5651591D+00  .7009068D+00
c       10.0  .2815717D+04  .2670988D+04  .2670988D+04  .2548618D+04
c       20.0  .4355828D+08  .4245497D+08  .4245497D+08  .4143553D+08
c       30.0  .7816723D+12  .7685320D+12  .7685320D+12  .7560546D+12
c       40.0  .1489477D+17  .1470740D+17  .1470740D+17  .1452709D+17
c       50.0  .2932554D+21  .2903079D+21  .2903079D+21  .2874492D+21
c
c         x      K0(x)         K0'(x)        K1(x)         K1'(x)
c       -------------------------------------------------------------
c        1.0  .4210244D+00 -.6019072D+00  .6019072D+00 -.1022932D+01
c       10.0  .1778006D-04 -.1864877D-04  .1864877D-04 -.1964494D-04
c       20.0  .5741238D-09 -.5883058D-09  .5883058D-09 -.6035391D-09
c       30.0  .2132477D-13 -.2167732D-13  .2167732D-13 -.2204735D-13
c       40.0  .8392861D-18 -.8497132D-18  .8497132D-18 -.8605289D-18
c       50.0  .3410168D-22 -.3444102D-22  .3444102D-22 -.3479050D-22
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MIK01A'
      write ( *, '(a)' ) '  Test IK01A'
      write ( *, '(a)' ) ' '
        WRITE(*,*)'  x       I0(x)          I0''(x)         I1(x)',
     &            '          I1''(x)'
        WRITE(*,*)'-------------------------------------------',
     &            '----------------------'
      do i = 0, 5
        if ( i == 0 ) then
          x = 1.0D+00
        else
          x = dble ( i * 10 )
        end if

        CALL IK01A(X,BI0,DI0,BI1,DI1,BK0,DK0,BK1,DK1)
        WRITE(*,20)X,BI0,DI0,BI1,DI1
      end do

        WRITE(*,*)
        WRITE(*,*)'  x       K0(x)          K0''(x)         K1(x)',
     &            '          K1''(x)'
        WRITE(*,*)'-------------------------------------------',
     &            '----------------------'
      do i = 0, 5
        if ( i == 0 ) then
          x = 1.0D+00
        else
          x = dble ( i * 10 )
        end if
        CALL IK01A(X,BI0,DI0,BI1,DI1,BK0,DK0,BK1,DK1)
        WRITE(*,20)X,BK0,DK0,BK1,DK1
      end do
10      FORMAT(3x 'x =',F5.1)
20      FORMAT(1X,F4.1,4D15.7)
      return
        END
      subroutine MIK01B ( )

c*********************************************************************72
c
cc MIK01B tests IK01B.
c
c  Modified:
c
c    04 April 2012
c
c
c       Purpose: This program computes the modified Bessel functions 
c                I0(x), I1(x), K0(x), K1(x), and their derivatives 
c                using subroutine IK01B
c       Input :  x   --- Argument ( x ò 0 )
c       Output:  BI0 --- I0(x)
c                DI0 --- I0'(x)
c                BI1 --- I1(x)
c                DI1 --- I1'(x)
c                BK0 --- K0(x)
c                DK0 --- K0'(x)
c                BK1 --- K1(x)
c                DK1 --- K1'(x)
c       Example:
c
c         x      I0(x)         I0'(x)        I1(x)         I1'(x)
c       -------------------------------------------------------------
c        1.0   .126607D+01   .565159D+00   .565159D+00   .700907D+00
c       10.0   .281572D+04   .267099D+04   .267099D+04   .254862D+04
c       20.0   .435583D+08   .424550D+08   .424550D+08   .414355D+08
c       30.0   .781672D+12   .768532D+12   .768532D+12   .756055D+12
c       40.0   .148948D+17   .147074D+17   .147074D+17   .145271D+17
c       50.0   .293255D+21   .290308D+21   .290308D+21   .287449D+21
c
c         x      K0(x)         K0'(x)        K1(x)         K1'(x)
c       -------------------------------------------------------------
c        1.0   .421024D+00  -.601907D+00   .601907D+00  -.102293D+01
c       10.0   .177801D-04  -.186488D-04   .186488D-04  -.196449D-04
c       20.0   .574124D-09  -.588306D-09   .588306D-09  -.603539D-09
c       30.0   .213248D-13  -.216773D-13   .216773D-13  -.220474D-13
c       40.0   .839286D-18  -.849713D-18   .849713D-18  -.860529D-18
c       50.0   .341017D-22  -.344410D-22   .344410D-22  -.347905D-22
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MIK01B'
      write ( *, '(a)' ) '  Test IK01B'
      write ( *, '(a)' ) ' '
        WRITE(*,*)'  x       I0(x)          I0''(x)         I1(x)',
     &            '          I1''(x)'
        WRITE(*,*)'-------------------------------------------',
     &            '----------------------'
      do i = 0, 5

        if ( i == 0 ) then
          x = 1.0D+00
        else
          x = dble ( i * 10 )
        end if

        CALL IK01B(X,BI0,DI0,BI1,DI1,BK0,DK0,BK1,DK1)
        WRITE(*,20)X,BI0,DI0,BI1,DI1
      end do

        WRITE(*,*)
        WRITE(*,*)'  x       K0(x)          K0''(x)         K1(x)',
     &            '          K1''(x)'
        WRITE(*,*)'-------------------------------------------',
     &            '----------------------'
      do i = 0, 5
        if ( i == 0 ) then
          x = 1.0D+00
        else
          x = dble ( i * 10 )
        end if
        CALL IK01B(X,BI0,DI0,BI1,DI1,BK0,DK0,BK1,DK1)
        WRITE(*,20)X,BK0,DK0,BK1,DK1
      end do
10      FORMAT(3x 'x =',F5.1)
20      FORMAT(1X,F4.1,4D15.7)
      return
        END
      subroutine MIKNA ( )

c*********************************************************************72
c
cc MIKNA tests IKNA.
c
c  Modified:
c
c    05 April 2012
c
c
c       Purpose: This program computes modified Bessel functions 
c                In(x) and Kn(x), and their derivatives using
c              subroutine IKNA
c       Input:   x --- Argument of In(x) and Kn(x) ( x ò 0 )
c                n --- Order of In(x) and Kn(x)
c                      ( n = 0,1,úúú, n ó 250 )
c       Output:  BI(n) --- In(x)
c                DI(n) --- In'(x)
c                BK(n) --- Kn(x)
c                DK(n) --- Kn'(x)
c       Example: Nmax = 5,    x = 10.0
c
c     n      In(x)          In'(x)         Kn(x)         Kn'(x)
c    ---------------------------------------------------------------
c     0   .2815717D+04   .2670988D+04   .1778006D-04  -.1864877D-04
c     1   .2670988D+04   .2548618D+04   .1864877D-04  -.1964494D-04
c     2   .2281519D+04   .2214685D+04   .2150982D-04  -.2295074D-04
c     3   .1758381D+04   .1754005D+04   .2725270D-04  -.2968563D-04
c     4   .1226491D+04   .1267785D+04   .3786144D-04  -.4239728D-04
c     5   .7771883D+03   .8378964D+03   .5754185D-04  -.6663236D-04
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION BI(0:250),DI(0:250),BK(0:250),DK(0:250)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MIKNA'
      write ( *, '(a)' ) '  Test IKNA'
      write ( *, '(a)' ) ' '

      n = 5
      x = 10.0D+00

        IF (N.LE.10) THEN
           NS=1
        ELSE
           ns = 5
        ENDIF
        CALL IKNA(N,X,NM,BI,DI,BK,DK)
        WRITE(*,*)'  n      In(x)          In''(x) ',
     &            '        Kn(x)         Kn''(x) '
        WRITE(*,*)' -------------------------------',
     &            '--------------------------------'
        DO 10 K=0,NM,NS
           WRITE(*,20)K,BI(K),DI(K),BK(K),DK(K)
10      CONTINUE
20      FORMAT(1X,I3,4D15.7)
25      FORMAT(3X,6HNmax =,I3,',    ',3Hx =,F5.1)

      return
        END
      subroutine MIKNB ( )

c*********************************************************************72
c
cc MIKNB tests IKNB.
c
c  Modified:
c
c    05 April 2012
c
c
c       Purpose: This program computes modified Bessel functions 
c                In(x) and Kn(x), and their derivatives using
c              subroutine IKNA
c       Input:   x --- Argument of In(x) and Kn(x) ( x ò 0 )
c                n --- Order of In(x) and Kn(x)
c                      ( n = 0,1,úúú, n ó 250 )
c       Output:  BI(n) --- In(x)
c                DI(n) --- In'(x)
c                BK(n) --- Kn(x)
c                DK(n) --- Kn'(x)
c       Example: Nmax = 5,    x = 10.0
c
c     n      In(x)          In'(x)         Kn(x)         Kn'(x)
c    ---------------------------------------------------------------
c     0   .2815717D+04   .2670988D+04   .1778006D-04  -.1864877D-04
c     1   .2670988D+04   .2548618D+04   .1864877D-04  -.1964494D-04
c     2   .2281519D+04   .2214685D+04   .2150982D-04  -.2295074D-04
c     3   .1758381D+04   .1754005D+04   .2725270D-04  -.2968563D-04
c     4   .1226491D+04   .1267785D+04   .3786144D-04  -.4239728D-04
c     5   .7771883D+03   .8378964D+03   .5754185D-04  -.6663236D-04
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION BI(0:250),DI(0:250),BK(0:250),DK(0:250)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MIKNB'
      write ( *, '(a)' ) '  Test IKNB'
      write ( *, '(a)' ) ' '

      n = 5
      x = 10.0D+00

        IF (N.LE.10) THEN
           NS=1
        ELSE
           ns = 5
        ENDIF
        CALL IKNB(N,X,NM,BI,DI,BK,DK)
        WRITE(*,*)'  n      In(x)          In''(x) ',
     &            '        Kn(x)         Kn''(x) '
        WRITE(*,*)' -------------------------------',
     &            '--------------------------------'
        DO 10 K=0,NM,NS
           WRITE(*,20)K,BI(K),DI(K),BK(K),DK(K)
10      CONTINUE
20      FORMAT(1X,I3,4D15.7)
25      FORMAT(3X,6HNmax =,I3,',    ',3Hx =,F5.1)

      return
        END
      subroutine MIKV ( )

c*********************************************************************72
c
cc MIKV tests IKV.
c
c  Modified:
c
c    05 April 2012
c
c
c       Purpose: This program computes modified Bessel functions 
c                Iv(x) and Kv(x) with an arbitrary order, and
c                their derivatives using subroutine IKV
c       Input :  x --- Argument ( x ò 0 )
c                v --- Order of Iv(x) and Kv(x)
c                      ( v = n+v0, 0 ó n ó 250 , 0 ó v0 < 1 )
c       Output:  BI(n) --- In+v0(x)
c                DI(n) --- In+v0'(x)
c                BK(n) --- Kn+v0(x)
c                DK(n) --- Kn+v0'(x)
c       Example: v = n+v0,   v0 = .25,   x = 10.0
c
c     n       Iv(x)          Iv'(x)         Kv(x)          Kv'(x)
c    ----------------------------------------------------------------
c     0   .28064359D+04  .26631677D+04  .17833184D-04 -.18709581D-04
c     1   .25930068D+04  .24823101D+04  .19155411D-04 -.20227611D-04
c     2   .21581842D+04  .21074153D+04  .22622037D-04 -.24245369D-04
c     3   .16218239D+04  .16310915D+04  .29335327D-04 -.32156018D-04
c     4   .11039987D+04  .11526244D+04  .41690000D-04 -.47053577D-04
c     5   .68342498D+03  .74520058D+03  .64771827D-04 -.75695209D-04
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        COMMON BI(0:250),DI(0:250),BK(0:250),DK(0:250)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MIKV'
      write ( *, '(a)' ) '  Test IKV'
      write ( *, '(a)' ) ' '

      v = 5.25D+00
      x = 10.0D+00

        N=INT(V)
        V0=V-N
        WRITE(*,30)V0,X
        IF (N.LE.10) THEN
           NS=1
        ELSE
           NS = 5
        ENDIF
        CALL IKV(V,X,VM,BI,DI,BK,DK)
        NM=INT(VM)
        WRITE(*,*)
        WRITE(*,*)'    n       Iv(x)           Iv''(x)          ',
     &            'Kv(x)           Kv''(x)'
        WRITE(*,*)'  -------------------------------------',
     &            '--------------------------------'
        DO 10 K=0,NM,NS
10         WRITE(*,20) K,BI(K),DI(K),BK(K),DK(K)
20      FORMAT(2X,I4,1X,4D16.8)
30      FORMAT(8X,'v = n+v0',',   ','v0 =',F7.5,',   ','x =',F5.1)
      return
        END
      subroutine MINCOB ( )

c*********************************************************************72
c
cc MINCOB tests INCOB.
c
c  Modified:
c
c    05 April 2012
c
c
c       Purpose: This program computes the incomplete beta
c                function Ix(a,b) using subroutine INCOB
c       Input :  a --- Parameter
c                b --- Parameter
c                x --- Argument ( 0 ó x ó 1 )
c       Output:  BIX --- Ix(a,b)
c       Example:
c                  a       b       x       Ix(a,b)
c                -----------------------------------
c                 1.0     3.0     .25     .57812500
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MINCOB'
      write ( *, '(a)' ) '  Test INCOB'
      write ( *, '(a)' ) ' '

      a = 1.0D+00
      b = 3.0D+00
      x = 0.25D+00
        WRITE(*,*)'   a       b       x       Ix(a,b)'
        WRITE(*,*)' -----------------------------------'
        CALL INCOB(A,B,X,BIX)
        WRITE(*,10)A,B,X,BIX
10      FORMAT(1X,F5.1,3X,F5.1,3X,F5.2,F14.8)
      return
        END
      subroutine MINCOG ( )

c*********************************************************************72
c
cc MINCOG tests INCOG.
c
c  Modified:
c
c    05 April 2012
c
c
c       Purpose: This program computes the incomplete gamma
c                function r(a,x), â(a,x) and P(a,x) using
c              subroutine INCOG
c       Input :  a   --- Parameter
c                x   --- Argument
c       Output:  GIN --- r(a,x)
c                GIM --- â(a,x)
c                GIP --- P(a,x)
c       Example:
c            a     x      r(a,x)         â(a,x)         P(a,x)
c           -------------------------------------------------------
c           3.0   5.0  .17506960D+01  .24930404D+00  .87534798D+00
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MINCOG'
      write ( *, '(a)' ) '  Test INCOG'
      write ( *, '(a)' ) ' '

      a = 3.0D+00
      x = 0.25D+00
        WRITE(*,*)'   a     x      r(a,x)         â(a,x)         P(a,x)'
        WRITE(*,*)' --------------------------------------------',
     &            '------------'
        CALL INCOG(A,X,GIN,GIM,GIP)
        WRITE(*,10)A,X,GIN,GIM,GIP
10      FORMAT(1X,F5.1,1X,F5.1,3D15.8)

      return
        END
      subroutine MITAIRY ( )

c*********************************************************************72
c
cc MITAIRY tests ITAIRY.
c
c  Modified:
c
c    06 April 2012
c
c
c       Purpose: This program computes the integrals of Airy 
c                functions using subroutine ITAIRY
c       Input  : x   --- Upper limit of the integral
c       Output : APT --- Integration of Ai(t) from 0 and x
c                BPT --- Integration of Bi(t) from 0 and x
c                ANT --- Integration of Ai(-t) from 0 and x
c                BNT --- Integration of Bi(-t) from 0 and x
c       Example:
c
c         x      Ai(t)dt       Bi(t)dt       Ai(-t)dt     Bi(-t)dt
c        ----------------------------------------------------------
c         5    .33328759   .32147832D+03    .71788220    .15873094
c        10    .33333333   .14780980D+09    .76569840    .01504043
c        15    .33333333   .49673090D+16    .68358063    .07202621
c        20    .33333333   .47447423D+25    .71173925   -.03906173
c        25    .33333333   .78920820D+35    .70489539    .03293190
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MITAIRY'
      write ( *, '(a)' ) '  Test ITAIRY'
      write ( *, '(a)' ) ' '
        WRITE(*,20)
        WRITE(*,30)

      do i = 5, 25, 5
        x = dble ( i )

        CALL ITAIRY(X,APT,BPT,ANT,BNT)
        WRITE(*,10)X,APT,BPT,ANT,BNT
      end do

10      FORMAT(1X,F5.1,F14.8,2X,D15.8,2F14.8)
20      FORMAT(3X,'x',8X,'Ai(t)dt',7X,'Bi(t)dt',9X,
     &        'Ai(-t)dt',6X,'Bi(-t)dt')
30      FORMAT(2X,'----------------------------------',
     &     '------------------------------')
        return
        END
      subroutine MITIKA ( )

c*********************************************************************72
c
cc MITIKA tests ITIKA.
c
c  Modified:
c
c    06 April 2012
c
c
c       Purpose: This program evaluates the integral of modified 
c                Bessel functions I0(t) and K0(t) with respect to t
c                from 0 to x using subroutine ITIKA
c       Input :  x  --- Upper limit of the integral  ( x ò 0 )
c       Output:  TI --- Integration of I0(t) from 0 to x
c                TK --- Integration of K0(t) from 0 to x
c       Example:
c                    x         I0(t)dt         K0(t)dt
c                 --------------------------------------
c                   5.0    .31848668D+02     1.56738739
c                  10.0    .29930445D+04     1.57077931
c                  15.0    .35262048D+06     1.57079623
c                  20.0    .44758593D+08     1.57079633
c                  25.0    .58991731D+10     1.57079633
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MITIKA'
      write ( *, '(a)' ) '  Test ITIKA'
      write ( *, '(a)' ) ' '

        WRITE(*,*)'   x         I0(t)dt         K0(t)dt'
        WRITE(*,*)' --------------------------------------'
      do i = 1, 5
        x = dble ( 5 * i )
        CALL ITIKA(X,TI,TK)
        WRITE(*,10)X,TI,TK
      end do
10      FORMAT(1X,F5.1,D17.8,F15.8)
      return
        END
      subroutine MITIKB ( )

c*********************************************************************72
c
cc MITIKB tests ITIKB.
c
c  Modified:
c
c    08 April 2012
c
c
c       Purpose: This program evaluates the integral of modified 
c                Bessel functions I0(t) and K0(t) with respect to t
c                from 0 to x using subroutine ITIKB
c       Input :  x  --- Upper limit of the integral  ( x ò 0 )
c       Output:  TI --- Integration of I0(t) from 0 to x
c                TK --- Integration of K0(t) from 0 to x
c       Example:
c                    x         I0(t)dt         K0(t)dt
c                 -------------------------------------
c                   5.0     .318487D+02       1.567387
c                  10.0     .299305D+04       1.570779
c                  15.0     .352619D+06       1.570796
c                  20.0     .447586D+08       1.570796
c                  25.0     .589919D+10       1.570796
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MITIKB'
      write ( *, '(a)' ) '  Test ITIKB'
      write ( *, '(a)' ) ' '
        WRITE(*,*)'   x         I0(t)dt         K0(t)dt'
        WRITE(*,*)'--------------------------------------'

      do i = 1, 5
        x = dble ( 5 * i )

        CALL ITIKB(X,TI,TK)
        WRITE(*,10)X,TI,TK
      end do

10      FORMAT(1X,F5.1,D16.6,F15.6)
      return
        END
      subroutine MITJYA ( )

c*********************************************************************72
c
cc MITJYA tests ITJYA.
c
c  Modified:
c
c    10 April 2012
c
c
c       Purpose: This program evaluates the integral of Bessel
c                functions J0(t) and Y0(t) with respect to t 
c                from 0 to x using subroutine ITJYA
c       Input :  x  --- Upper limit of the integral ( x ò 0 )
c       Output:  TJ --- Integration of J0(t) from 0 to x
c                TY --- Integration of Y0(t) from 0 to x
c       Example:
c                   x         J0(t)dt          Y0(t)dt
c                ---------------------------------------
c                  5.0       .71531192       .19971938
c                 10.0      1.06701130       .24129032
c                 15.0      1.20516194       .00745772
c                 20.0      1.05837882      -.16821598
c                 25.0       .87101492      -.09360793
c                 30.0       .88424909       .08822971
c
c
        DOUBLE PRECISION X,TJ,TY

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MITJYA'
      write ( *, '(a)' ) '  Test ITJYA'
      write ( *, '(a)' ) ' '

        WRITE(*,*)'   x         J0(t)dt          Y0(t)dt'
        WRITE(*,*)'---------------------------------------'
      do i = 1, 6
        x = dble ( 5 * i )
        CALL ITJYA(X,TJ,TY)
        WRITE(*,10)X,TJ,TY
      end do

10      FORMAT(1X,F5.1,2F16.8)
      return
        END
      subroutine MITJYB ( )

c*********************************************************************72
c
cc MITJYB tests ITJYB.
c
c  Modified:
c
c    10 April 2012
c
c
c       Purpose: This program evaluates the integral of Bessel
c                functions J0(t) and Y0(t) with respect to t 
c                from 0 to x using subroutine ITJYB
c       Input :  x  --- Upper limit of the integral ( x ò 0 )
c       Output:  TJ --- Integration of J0(t) from 0 to x
c                TY --- Integration of Y0(t) from 0 to x
c       Example:
c                   x         J0(t)dt          Y0(t)dt
c                ---------------------------------------
c                  5.0       .71531192       .19971938
c                 10.0      1.06701130       .24129032
c                 15.0      1.20516194       .00745772
c                 20.0      1.05837882      -.16821598
c                 25.0       .87101492      -.09360793
c                 30.0       .88424909       .08822971
c
c
        DOUBLE PRECISION X,TJ,TY

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MITJYB'
      write ( *, '(a)' ) '  Test ITJYB'
      write ( *, '(a)' ) ' '
        WRITE(*,*)'   x         J0(t)dt          Y0(t)dt'
        WRITE(*,*)'---------------------------------------'
      do i = 1, 6
        x = dble ( 5 * i )
        CALL ITJYB(X,TJ,TY)
        WRITE(*,10)X,TJ,TY
      end do
10      FORMAT(1X,F5.1,2F16.8)
      return
        END
      subroutine MITSH0 ( )

c*********************************************************************72
c
cc MITSH0 tests ITSH0.
c
c  Modified:
c
c    12 April 2012
c
c
c       Purpose: This program evaluates the integral of 
c                Struve function H0(t) with respect to t 
c                from 0 and x using subroutine ITSH0
c       Input :  x   --- Upper limit  ( x ò 0 )
c       Output:  TH0 --- Integration of H0(t) from 0 and x
c       Example:
c                    x        H0(t)dt
c                 ----------------------
c                   0.0       .0000000
c                   5.0      2.0442437
c                  10.0      2.5189577
c                  15.0      2.5415824
c                  20.0      2.5484517
c                  30.0      3.0625848
c                  40.0      3.1484123
c                  50.0      3.2445168
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MITSH0'
      write ( *, '(a)' ) '  Test ITSH0'
      write ( *, '(a)' ) ' '
        WRITE(*,*)'   x         H0(t)dt'
        WRITE(*,*)'----------------------'
      do i = 0, 10
        x = dble ( 5 * i )
        CALL ITSH0(X,TH0)
        WRITE(*,10)X,TH0
      end do
10      FORMAT(1X,F5.1,E16.7)

      return
        END
      subroutine MITSL0 ( )

c*********************************************************************72
c
cc MITSL0 tests ITSL0.
c
c  Modified:
c
c    12 April 2012
c
c
c       Purpose: This program evaluates the integral of modified 
c                Struve function L0(t) with respect to t from 0  
c                to x using subroutine ITSL0
c       Input :  x   --- Upper limit  ( x ò 0 )
c       Output:  TL0 --- Integration of L0(t) from 0 to x
c       Example:
c                      x        L0(t)dt
c                   -----------------------
c                     0.0    .0000000D+00
c                     5.0    .3003079D+02
c                    10.0    .2990773D+04
c                    15.0    .3526179D+06
c                    20.0    .4475860D+08
c                    30.0    .7955389D+12
c                    40.0    .1508972D+17
c                    50.0    .2962966D+21
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MITSL0'
      write ( *, '(a)' ) '  Test ITSL0'
      write ( *, '(a)' ) ' '

        WRITE(*,*)'   x        L0(t)dt'
        WRITE(*,*)'-----------------------'
      do i = 0, 10
        x = dble ( 5 * i )
        CALL ITSL0(X,TL0)
        WRITE(*,10)X,TL0
      end do

10      FORMAT(1X,F5.1,D16.7)

      return
        END
      subroutine MITTH0 ( )

c*********************************************************************72
c
cc MITTH0 tests ITTH0.
c
c  Modified:
c
c    12 April 2012
c
c
c       Purpose: This program evaluates the integral of H0(t)/t 
c                with respect to t from x to infinity using
c              subroutine ITTH0
c       Input :  x   --- Lower limit  ( x ò 0 )
c       Output:  TTH --- Integration of H0(t)/t from x to infinity
c       Example:
c                    x        H0(t)/t dt
c                 -----------------------
c                   0.0      1.57079633
c                   5.0       .07954575
c                  10.0       .04047175
c                  15.0       .04276558
c                  20.0       .04030796
c                  30.0       .01815256
c                  40.0       .01621331
c                  50.0       .01378661
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MITTH0'
      write ( *, '(a)' ) '  Test ITTH0'
      write ( *, '(a)' ) ' '

        WRITE(*,*) '   x       H0(t)/t dt'
        WRITE(*,*)'-----------------------'
      do i = 0, 10
        x = dble ( 5 * i )
        CALL ITTH0(X,TTH)
        WRITE(*,10)X,TTH
      end do

10      FORMAT(1X,F5.1,1X,E16.8)

      return
        END
      subroutine MITTIKA ( )

c*********************************************************************72
c
cc MITTIKA tests ITTIKA.
c
c  Modified:
c
c    12 April 2012
c
c
c       Purpose: This program computes the integral of [I0(t)-1]/t
c                with respect to t from 0 to x and K0(t)/t with 
c                respect to t from x to ì using subroutine ITTIKA
c       Input :  x   --- Variable in the limits  ( x ò 0 )
c       Output:  TTI --- Integration of [I0(t)-1]/t from 0 to x
c                TTK --- Integration of K0(t)/t from x to ì
c       Example:
c                   x     [1-I0(t)]/tdt     K0(t)/tdt
c                ---------------------------------------
c                  5.0   .71047763D+01   .58635626D-03
c                 10.0   .34081537D+03   .15629282D-05
c                 15.0   .25437619D+05   .59837472D-08
c                 20.0   .23673661D+07   .26790545D-10
c                 25.0   .24652751D+09   .13100706D-12
c
c
        DOUBLE PRECISION X,TTI,TTK

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MITTIKA'
      write ( *, '(a)' ) '  Test ITTIKA'
      write ( *, '(a)' ) ' '
        WRITE(*,*)'   x     [1-I0(t)]/tdt     K0(t)/tdt'
        WRITE(*,*)'---------------------------------------'
      do i = 1, 5
        x = dble ( 5 * i )
        CALL ITTIKA(X,TTI,TTK)
        WRITE(*,10)X,TTI,TTK
      end do
10      FORMAT(1X,F5.1,2D16.8)

      return
        END
       subroutine MITTIKB ( )

c*********************************************************************72
c
cc MITTIKB tests ITTIKB.
c
c  Modified:
c
c    12 April 2012
c
c
c       Purpose: This program computes the integral of [I0(t)-1]/t
c                with respect to t from 0 to x and K0(t)/t with 
c                respect to t from x to ì using subroutine ITTIKB
c       Input :  x   --- Upper limit of the integral
c       Output:  TTI --- Integration of [I0(t)-1]/t from 0 to x
c                TTK --- Integration of K0(t)/t from x to ì
c       Example:
c                   x     [1-I0(t)]/tdt      K0(t)/tdt
c                ---------------------------------------
c                  5.0     .710478D+01     .586361D-03
c                 10.0     .340811D+03     .156293D-05
c                 15.0     .254373D+05     .598363D-08
c                 20.0     .236735D+07     .267906D-10
c                 25.0     .246534D+09     .131007D-12
c
c
        DOUBLE PRECISION X,TTI,TTK

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MITTIKB'
      write ( *, '(a)' ) '  Test ITTIKB'
      write ( *, '(a)' ) ' '
        WRITE(*,*)'   x     [1-I0(t)]/tdt      K0(t)/tdt'
        WRITE(*,*)'---------------------------------------'
      do i = 1, 5
        x = dble ( 5 * i )
        CALL ITTIKB(X,TTI,TTK)
        WRITE(*,10)X,TTI,TTK
      end do
10      FORMAT(1X,F5.1,2D16.6)

      return
        END
      subroutine MITTJYA ( )

c*********************************************************************72
c
cc MITTJYA tests ITTJYA.
c
c  Modified:
c
c    12 April 2012
c
c
c       Purpose: This program computes the integral of [1-J0(t)]/t 
c                with respect to t from 0 to x and Y0(t)/t with 
c                respect to t from x to ì using subroutine ITTJYA
c       Input :  x   --- Variable in the limits  ( x ò 0 )
c       Output:  TTJ --- Integration of [1-J0(t)]/t from 0 to x
c                TTY --- Integration of Y0(t)/t from x to ì
c       Example:
c                  x       [1-J0(t)]/tdt       Y0(t)/tdt
c               -------------------------------------------
c                 5.0     .15403472D+01    -.46322055D-01
c                10.0     .21778664D+01    -.22987934D-01
c                15.0     .25785507D+01     .38573574D-03
c                20.0     .28773106D+01     .85031527D-02
c                25.0     .31082313D+01     .35263393D-02
c
c
        DOUBLE PRECISION X,TTJ,TTY
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MITTJYA'
      write ( *, '(a)' ) '  Test ITTJYA'
      write ( *, '(a)' ) ' '

        WRITE(*,*)'   x      [1-J0(t)]/tdt       Y0(t)/tdt'
        WRITE(*,*)'-------------------------------------------'
      do i = 1, 6
        x = dble ( 5 * i )
        CALL ITTJYA(X,TTJ,TTY)
        WRITE(*,10)X,TTJ,TTY
      end do
10      FORMAT(1X,F5.1,2D18.8)
      return
        END
      subroutine MITTJYB ( )

c*********************************************************************72
c
cc MITTJYB tests ITTJYB.
c
c  Modified:
c
c    11 April 2012
c
c
c       Purpose: This program computes the integral of [1-J0(t)]/t 
c                with respect to t from 0 to x and Y0(t)/t with 
c                respect to t from x to ì using subroutine ITTJYB
c       Input :  x   --- Variable in the limits  ( x ò 0 )
c       Output:  TTJ --- Integration of [1-J0(t)]/t from 0 to x
c                TTY --- Integration of Y0(t)/t from x to ì
c       Example:
c                  x      [1-J0(t)]/tdt       Y0(t)/tdt
c                ----------------------------------------
c                 5.0     .1540347D+01    -.4632208D-01
c                10.0     .2177866D+01    -.2298791D-01
c                15.0     .2578551D+01     .3857453D-03
c                20.0     .2877311D+01     .8503154D-02
c                25.0     .3108231D+01     .3526339D-02
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MITTJYB'
      write ( *, '(a)' ) '  Test ITTJYB'
      write ( *, '(a)' ) ' '

        WRITE(*,*)'   x     [1-J0(t)]/tdt      Y0(t)/tdt'
        WRITE(*,*)'----------------------------------------'
      do i = 1, 5
        x = dble ( 5 * i )
        CALL ITTJYB(X,TTJ,TTY)
        WRITE(*,10)X,TTJ,TTY
      end do
10      FORMAT(1X,F5.1,2D17.7)
       return
        END
      subroutine MJDZO ( )

c*********************************************************************72
c
cc MJDZO tests JDZO.
c
c  Modified:
c
c    06 April 2012
c
c
c       Purpose: This program computes the zeros of Bessel functions 
c                Jn(x) and Jn'(x), and arranges them in the order 
c                of their values 
c       Input :  NT    --- Number of total zeros ( NT ó 1200 )
c       Output:  ZO(L) --- Value of the L-th zero of Jn(x) and 
c                          Jn'(x)
c                N(L)  --- n, order of Jn(x) or Jn'(x) associated
c                          with the L-th zero
c                M(L)  --- m, serial number of the zeros of Jn(x)
c                          or Jn'(x) associated with the L-th zero
c                          ( L is the serial number of all the
c                            zeros of Jn(x) and Jn'(x) )
c                P(L)  --- TM or TE, a code for designating the
c                          zeros of Jn(x) or Jn'(x)
c                          In the waveguide applications, the zeros
c                          of Jn(x) correspond to TM modes and those
c                          of Jn'(x) correspond to TE modes.
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        CHARACTER P(1400)*4
        DIMENSION N(1400),M(1400),ZO(1400)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MJDZO'
      write ( *, '(a)' ) '  Test JDZO'
      write ( *, '(a)' ) ' '

      nt = 9

        WRITE(*,60)NT
        CALL JDZO(NT,N,M,P,ZO)
        WRITE(*,*)
        KS=NT/101+1
        DO 55 K0=1,KS
           WRITE(*,*)' Table           Zeros of Bessel',
     &               ' functions Jn(x) and Jn''(x)'
           WRITE(*,*)
           WRITE(*,*)' ----------------------------------',
     &               '----------------------------------'
           DO 50 K=1,50
              J1=100*(K0-1)+K+1
              J2=J1+50
              IF (J1.LE.NT+1.AND.J2.LE.NT+1) THEN
                  WRITE(*,65)J1-1,P(J1),N(J1),M(J1),ZO(J1),
     &                 J2-1,P(J2),N(J2),M(J2),ZO(J2)
              ELSE IF (J1.LE.NT+1.AND.J2.GT.NT+1) THEN
                  WRITE(*,65)J1-1,P(J1),N(J1),M(J1),ZO(J1)
              ENDIF
50         CONTINUE
           WRITE(*,*)' ----------------------------------',
     &               '----------------------------------'
55      CONTINUE
60      FORMAT(1X,'Total number of the zeros:',I5)
65      FORMAT(1X,I4,3X,A2,I4,2H -,I2,F14.8,3X,1H|,2X,I4,
     &         3X,A2,I4,2H -,I2,F14.8)
75      FORMAT(/)
      return
        END
      subroutine MJELP ( )

c*********************************************************************72
c
cc MJELP tests JELP.
c
c  Modified:
c
c    06 April 2012
c
c
c       Purpose: This program computes Jacobian elliptic functions 
c                sn u, cn u and dn u using subroutine JELP
c       Input  : u   --- Argument of Jacobian elliptic fuctions
c                Hk  --- Modulus k ( 0 ó k ó 1 )
c       Output : ESN --- sn u
c                ECN --- cn u
c                EDN --- dn u
c                EPH --- phi ( in degrees )
c       Example:
c                k = .5, ( K(k) = 1.68575035 ), and u = u0*K
c
c                u0       phi       sn u        cn u        dn u
c              ----------------------------------------------------
c               0.0      .0000    .0000000   1.0000000   1.0000000
c               0.5    47.0586    .7320508    .6812500    .9306049
c               1.0    90.0000   1.0000000    .0000000    .8660254
c               1.5   132.9414    .7320508   -.6812500    .9306049
c               2.0   180.0000    .0000000  -1.0000000   1.0000000
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      double precision k

      k = 0.5D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MJELP'
      write ( *, '(a)' ) '  Test JELP'
      write ( *, '(a)' ) ' '

        WRITE(*,*)'   k        u          phi        sn u',
     &            '        cn u        dn u'
        WRITE(*,*)' -------------------------------------',
     &            '---------------------------'
      do i = 0, 4
        u = dble ( i ) / 2.0D+00
        CALL JELP(U,K,ESN,ECN,EDN,EPH)
        WRITE(*,10)K,U,EPH,ESN,ECN,EDN
      end do

10      FORMAT(1X,F5.3,F12.7,2X,F9.5,3F12.7)
      return
        END
      subroutine MJY01A ( )

c*********************************************************************72
c
cc MJY01A tests JY01A.
c
c  Modified:
c
c    11 April 2012
c
c
c       Purpose: This program computes the Bessel functions  
c                Jn(x) and Yn(x) ( n=0,1 ) and their derivatives 
c                using subroutine JY01A
c       Input :  x   --- Argument of Jn(x) & Yn(x) ( x ò 0 )
c       Output:  BJ0 --- J0(x)
c                DJ0 --- J0'(x)
c                BJ1 --- J1(x)
c                DJ1 --- J1'(x)
c                BY0 --- Y0(x)
c                DY0 --- Y0'(x)
c                BY1 --- Y1(x)
c                DY1 --- Y1'(x)
c       Example:
c
c        x       J0(x)        J0'(x)       J1(x)        J1'(x)
c       ---------------------------------------------------------
c        1     .76519769   -.44005059    .44005059    .32514710
c        5    -.17759677    .32757914   -.32757914   -.11208094
c       10    -.24593576   -.04347275    .04347275   -.25028304
c       20     .16702466   -.06683312    .06683312    .16368301
c       30    -.08636798    .11875106   -.11875106   -.08240961
c       40     .00736689   -.12603832    .12603832    .00421593
c       50     .05581233    .09751183   -.09751183    .05776256
c
c        x       Y0(x)        Y0'(x)       Y1(x)        Y1'(x)
c      ---------------------------------------------------------
c        1     .08825696    .78121282   -.78121282    .86946979
c        5    -.30851763   -.14786314    .14786314   -.33809025
c       10     .05567117   -.24901542    .24901542    .03076962
c       20     .06264060    .16551161   -.16551161    .07091618
c       30    -.11729573   -.08442557    .08442557   -.12010992
c       40     .12593642    .00579351   -.00579351    .12608125
c       50    -.09806500    .05679567   -.05679567   -.09692908
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      double precision x_test(7)

      save x_test 

      data x_test / 
     &  1.0D+00, 5.0D+00, 10.0D+00, 20.0D+00, 30.0D+00,
     &  40.0D+00, 50.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MJY01A'
      write ( *, '(a)' ) '  Test JY01A'
      write ( *, '(a)' ) ' '

        WRITE(*,*)
        WRITE(*,*)'  x          Y0(x)          Y0''(x)         Y1(x)',
     &            '          Y1''(x)'
        WRITE(*,*)'------------------------------------------',
     &            '----------------------------'
      do i = 1, 7
        x = x_test(i)
        CALL JY01A(X,BJ0,DJ0,BJ1,DJ1,BY0,DY0,BY1,DY1)
        WRITE(*,10)X,BJ0,DJ0,BJ1,DJ1

        WRITE(*,10)X,BY0,DY0,BY1,DY1
      end do
10      FORMAT(1X,F5.1,4E16.8)
20      FORMAT(3X,'x =',F5.1)
        return
        END
      subroutine MJY01B ( )

c*********************************************************************72
c
cc MJY01B tests JY01B.
c
c  Modified:
c
c    11 April 2012
c
c
c       Purpose: This program computes Bessel functions Jn(x) 
c                and Yn(x) ( n=0,1 ) and their derivatives 
c                using subroutine JY01B
c       Input :  x   --- Argument of Jn(x) & Yn(x) ( x ò 0 )
c       Output:  BJ0 --- J0(x)
c                DJ0 --- J0'(x)
c                BJ1 --- J1(x)
c                DJ1 --- J1'(x)
c                BY0 --- Y0(x)
c                DY0 --- Y0'(x)
c                BY1 --- Y1(x)
c                DY1 --- Y1'(x)
c       Example:
c
c        x       J0(x)        J0'(x)       J1(x)        J1'(x)
c       ---------------------------------------------------------
c        1     .76519769   -.44005059    .44005059    .32514710
c        5    -.17759677    .32757914   -.32757914   -.11208094
c       10    -.24593576   -.04347275    .04347275   -.25028304
c       20     .16702466   -.06683312    .06683312    .16368301
c       30    -.08636798    .11875106   -.11875106   -.08240961
c       40     .00736689   -.12603832    .12603832    .00421593
c       50     .05581233    .09751183   -.09751183    .05776256
c
c        x       Y0(x)        Y0'(x)       Y1(x)        Y1'(x)
c      ---------------------------------------------------------
c        1     .08825696    .78121282   -.78121282    .86946979
c        5    -.30851763   -.14786314    .14786314   -.33809025
c       10     .05567117   -.24901542    .24901542    .03076962
c       20     .06264060    .16551161   -.16551161    .07091618
c       30    -.11729573   -.08442557    .08442557   -.12010992
c       40     .12593642    .00579351   -.00579351    .12608125
c       50    -.09806500    .05679567   -.05679567   -.09692908
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      double precision x_test(7)

      save x_test 

      data x_test / 
     &  1.0D+00, 5.0D+00, 10.0D+00, 20.0D+00, 30.0D+00,
     &  40.0D+00, 50.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MJY01B'
      write ( *, '(a)' ) '  Test JY01B'
      write ( *, '(a)' ) ' '
        WRITE(*,*)'  x          J0(x)          J0''(x)         J1(x)',
     &            '          J1''(x)'
        WRITE(*,*)'------------------------------------------',
     &            '----------------------------'
      do i = 1, 7
        x = x_test(i)
        CALL JY01B(X,BJ0,DJ0,BJ1,DJ1,BY0,DY0,BY1,DY1)
        WRITE(*,10)X,BJ0,DJ0,BJ1,DJ1
        WRITE(*,*)
        WRITE(*,*)'  x          Y0(x)          Y0''(x)         Y1(x)',
     &            '          Y1''(x)'
        WRITE(*,*)'------------------------------------------',
     &            '----------------------------'
        WRITE(*,10)X,BY0,DY0,BY1,DY1
      end do

10      FORMAT(1X,F5.1,4E16.8)
20      FORMAT(3X,'x =',F5.1)
      return
        END
      subroutine MJYNA ( )

c*********************************************************************72
c
cc MJYNA tests JYNA.
c
c  Modified:
c
c    11 April 2012
c
c
c       Purpose: This program computes Bessel functions  
c                Jn(x) and Yn(x), and their derivatives 
c                using subroutine JYNA
c       Input :  x --- Argument of Jn(x) & Yn(x)  ( x ò 0 )
c                n --- Order of Jn(x) & Yn(x)
c                      ( n = 0,1,2,úúú, n ó 250 )
c       Output:  BJ(n) --- Jn(x)
c                DJ(n) --- Jn'(x)
c                BY(n) --- Yn(x)
c                DY(n) --- Yn'(x)
c       Example:
c                x = 10.0
c
c                n        Jn(x)           Jn'(x)
c              -------------------------------------
c                0    -.2459358D+00   -.4347275D-01
c               10     .2074861D+00    .8436958D-01
c               20     .1151337D-04    .2011954D-04
c               30     .1551096D-11    .4396479D-11
c
c                n        Yn(x)           Yn'(x)
c              -------------------------------------
c                0     .5567117D-01   -.2490154D+00
c               10    -.3598142D+00    .1605149D+00
c               20    -.1597484D+04    .2737803D+04
c               30    -.7256142D+10    .2047617D+11
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION BJ(0:250),BY(0:250),DJ(0:250),DY(0:250)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MJYNA'
      write ( *, '(a)' ) '  Test JYNA'
      write ( *, '(a)' ) ' '
      WRITE(*,*)
     &  '  n        Jn(x)           Jn''(x)' // 
     &  '        Yn(x)           Yn''(x)'
      WRITE(*,*)
     &  '--------------------------------------------------------------'

      x = 10.0D+00
      n = 30
      ns = 10

        CALL JYNA(N,X,NM,BJ,DJ,BY,DY)

        DO 10 K=0,NM,NS
10         WRITE(*,40)K,BJ(K),DJ(K),BY(K),DY(K)

40      FORMAT(1X,I3,1X,2D16.7,2D16.7)

      return
      END
      subroutine MJYNB ( )

c*********************************************************************72
c
cc MJYNB tests JYNB.
c
c  Modified:
c
c    11 April 2012
c
c
c       Purpose: This program computes Bessel functions 
c                Jn(x) and Yn(x), and their derivatives 
c                using subroutine JYNB
c       Input :  x --- Argument of Jn(x) & Yn(x)  ( x ò 0 )
c                n --- Order of Jn(x) & Yn(x)
c                      ( n = 0,1,2,úúú, n ó 250 )
c       Output:  BJ(n) --- Jn(x)
c                DJ(n) --- Jn'(x)
c                BY(n) --- Yn(x)
c                DY(n) --- Yn'(x)
c       Example:
c                x = 10.0
c
c                n        Jn(x)           Jn'(x)
c              -------------------------------------
c                0    -.2459358D+00   -.4347275D-01
c               10     .2074861D+00    .8436958D-01
c               20     .1151337D-04    .2011954D-04
c               30     .1551096D-11    .4396479D-11
c
c                n        Yn(x)           Yn'(x)
c              -------------------------------------
c                0     .5567117D-01   -.2490154D+00
c               10    -.3598142D+00    .1605149D+00
c               20    -.1597484D+04    .2737803D+04
c               30    -.7256142D+10    .2047617D+11
c
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION BJ(0:250),BY(0:250),DJ(0:250),DY(0:250)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MJYNB'
      write ( *, '(a)' ) '  Test JYNB'
      write ( *, '(a)' ) ' '
      WRITE(*,*)
     &  '  n        Jn(x)           Jn''(x)' // 
     &  '        Yn(x)           Yn''(x)'
      WRITE(*,*)
     &  '--------------------------------------------------------------'

      x = 10.0D+00
      n = 30
      ns = 10
      CALL JYNB(N,X,NM,BJ,DJ,BY,DY)

      DO 10 K=0,NM,NS
10         WRITE(*,40)K,BJ(K),DJ(K),BY(K),DY(K)

40      FORMAT(1X,I3,1X,2D16.7,2D16.7)

      return
      END
      subroutine MJYV ( )

c*********************************************************************72
c
cc MJYV tests JYV.
c
c  Modified:
c
c    06 April 2012
c
c
c       Purpose: This program computes Bessel functions Jv(x) and
c                Yv(x) and their derivatives using subroutine JYV
c       Input :  x --- Argument of Jv(x) and Yv(x)
c                v --- Order of Jv(x) and Yv(x)
c                      ( v = n+v0,  0 ó n ó 250, 0 ó v0 < 1 )
c       Output:  BJ(n) --- Jn+v0(x)
c                DJ(n) --- Jn+v0'(x)
c                BY(n) --- Yn+v0(x)
c                DY(n) --- Yn+v0'(x)
c       Example: Compute Jv(x) and Yv(x) and their derivatives
c                for v = 0.25(1.0)5.25 and x = 10.0
c                Computation results:
c
c                v =  5.25,      x = 10.00
c
c        v        Jv(x)         Jv'(x)        Yv(x)         Yv'(x)
c       ------------------------------------------------------------
c        .25   -.20639379    -.13476340     .14493044    -.21381777
c       1.25    .12960355    -.22259423     .21744103     .11775031
c       2.25    .23879467     .07587475    -.09057018     .23781932
c       3.25   -.02214595     .24599211    -.25819761    -.00665596
c       4.25   -.25318954     .08545961    -.07725827    -.22536285
c       5.25   -.19306516    -.15183033     .19252809    -.17833551
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        COMMON BJ(0:250),DJ(0:250),BY(0:250),DY(0:250)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MJYV'
      write ( *, '(a)' ) '  Test JYV'
      write ( *, '(a)' ) ' '

      v = 5.25D+00
      x = 10.0D+00

        IF (V.LE.8) THEN
           NS=1
        ELSE
           NS = 5
        ENDIF
        CALL JYV(V,X,VM,BJ,DJ,BY,DY)  
        NM=INT(VM)
        V0=VM-NM
        WRITE(*,*)'   v         Jv(x)           Jv''(x)',
     &            '          Yv(x)           Yv''(x)'
        WRITE(*,*)' ---------------------------------------------',
     &            '------------------------'
        DO 10 K=0,NM,NS
           VK=K+V0
10         WRITE(*,15)VK,BJ(K),DJ(K),BY(K),DY(K)
15      FORMAT(1X,F6.2,4D16.8)
20      FORMAT(8X,3Hv =,F6.2,',    ',3Hx =,F6.2)

      return
        END
      subroutine MJYZO ( )

c*********************************************************************72
c
cc MJYZO tests JYZO.
c
c  Modified:
c
c    07 April 2012
c
c
c       Purpose: This program computes the zeros of Bessel 
c                functions Jn(x), Yn(x), and their derivatives 
c                using subroutine JYZO
c       Input :  n --- Order of Bessel functions ( n ó 100 )
c                NT --- Number of zeros
c       Output:  RJ0(m) --- m-th zero of Jn(x),  m=1,2,...,NT
c                RJ1(m) --- m-th zero of Jn'(x), m=1,2,...,NT
c                RY0(m) --- m-th zero of Yn(x),  m=1,2,...,NT
c                RY1(m) --- m-th zero of Yn'(x), m=1,2,...,NT
c       Example: n = 1, NT =5
c
c      Zeros of Bessel funcions Jn(x), Yn(x) and their derivatives
c                                 ( n = 1 )
c       m       jnm           j'nm          ynm           y'nm
c      -----------------------------------------------------------
c       1     3.8317060     1.8411838     2.1971413     3.6830229
c       2     7.0155867     5.3314428     5.4296810     6.9415000
c       3    10.1734681     8.5363164     8.5960059    10.1234047
c       4    13.3236919    11.7060049    11.7491548    13.2857582
c       5    16.4706301    14.8635886    14.8974421    16.4400580
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION RJ0(101),RJ1(101),RY0(101),RY1(101)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MJYZO'
      write ( *, '(a)' ) '  Test JYZO'
      write ( *, '(a)' ) ' '

      n = 1
      nt = 5

        CALL JYZO(N,NT,RJ0,RJ1,RY0,RY1)
        WRITE(*,30)
        WRITE(*,40)N
        WRITE(*,*)'  m       jnm           j''nm          ynm',
     &            '           y''nm'
        WRITE(*,*)' ----------------------------------------',
     &            '-------------------'
        DO 10 M=1,NT
10         WRITE(*,50)M,RJ0(M),RJ1(M),RY0(M),RY1(M)
30      FORMAT(2X,'Zeros of Bessel funcions Jn(x), Yn(x)',
     &         ' and their derivatives')
40      FORMAT(30X,'( n =',I2,' )')
50      FORMAT(1X,I3,4F14.7)
      return
        END
      subroutine MKLVNA ( )

c*********************************************************************72
c
cc MKLVNA tests KLVNA.
c
c  Modified:
c
c    06 April 2012
c
c
c       Purpose: This program computes Kelvin functions ber x, 
c                bei x, ker x and kei x, and their derivatives 
c                using subroutine KLVNA
c       Input :  x   --- Argument of Kelvin functions
c       Output:  BER --- ber x
c                BEI --- bei x
c                GER --- ker x
c                GEI --- kei x
c                DER --- ber'x
c                DEI --- bei'x
c                HER --- ker'x
c                HEI --- kei'x
c       Example:
c
c      x       ber x          bei x          ker x          kei x
c    -----------------------------------------------------------------
c      0    .1000000D+01    0              ì            -.7853982D+00
c      5   -.6230082D+01   .1160344D+00  -.1151173D-01   .1118759D-01
c     10    .1388405D+03   .5637046D+02   .1294663D-03  -.3075246D-03
c     15   -.2967255D+04  -.2952708D+04  -.1514347D-07   .7962894D-05
c     20    .4748937D+05   .1147752D+06  -.7715233D-07  -.1858942D-06
c
c      x       ber'x          bei'x          ker'x          kei'x
c    -----------------------------------------------------------------
c      0     0              0            - ì              0
c      5   -.3845339D+01  -.4354141D+01   .1719340D-01  -.8199865D-03
c     10    .5119526D+02   .1353093D+03  -.3155969D-03   .1409138D-03
c     15    .9105533D+02  -.4087755D+04   .5644678D-05  -.5882223D-05
c     20   -.4880320D+05   .1118550D+06  -.7501859D-07   .1906243D-06
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MKLVNA'
      write ( *, '(a)' ) '  Test KLVNA'
      write ( *, '(a)' ) ' '

        WRITE(*,*)'   x        ber x           bei x',
     &            '           ker x           kei x'
        WRITE(*,*)'--------------------------------',
     &            '--------------------------------------'
      do i = 0, 4
        x = dble ( 5 * i )
        CALL KLVNA(X,BER,BEI,GER,GEI,DER,DEI,HER,HEI)
        WRITE(*,10)X,BER,BEI,GER,GEI
      end do

        WRITE(*,*)
        WRITE(*,*)'   x        ber''x           bei''x',
     &            '           ker''x           kei''x'
        WRITE(*,*)'--------------------------------',
     &            '--------------------------------------'
      do i = 0, 4
        x = dble ( 5 * i )
        CALL KLVNA(X,BER,BEI,GER,GEI,DER,DEI,HER,HEI)
        WRITE(*,10)X,DER,DEI,HER,HEI
      end do

10      FORMAT(1X,F5.1,4D16.8)
      return
        END
      subroutine MKLVNB ( )

c*********************************************************************72
c
cc MKLVNB tests KLVNB.
c
c  Modified:
c
c    11 April 2012
c
c
c       Purpose: This program computes Kelvin functions ber x, 
c                bei x, ker x and kei x, and their derivatives 
c                using subroutine KLVNB
c       Input :  x   --- Argument of Kelvin functions
c       Output:  BER --- ber x
c                BEI --- bei x
c                GER --- ker x
c                GEI --- kei x
c                DER --- ber'x
c                DEI --- bei'x
c                HER --- ker'x
c                HEI --- kei'x
c       Example:
c     x       ber x         bei x         ker x         kei x
c   -------------------------------------------------------------
c     0    .100000D+01   .000000D+00    ì           -.785398D+00
c     5   -.623008D+01   .116034D+00  -.115117D-01   .111876D-01
c    10    .138840D+03   .563705D+02   .129466D-03  -.307525D-03
c    15   -.296725D+04  -.295271D+04  -.151433D-07   .796289D-05
c    20    .474894D+05   .114775D+06  -.771523D-07  -.185894D-06
c
c     x       ber'x         bei'x         ker'x         kei'x
c   -------------------------------------------------------------
c     0    .000000D+00   .000000D+00  - ì            .000000D+00
c     5   -.384534D+01  -.435414D+01   .171934D-01  -.819979D-03
c    10    .511952D+02   .135309D+03  -.315597D-03   .140914D-03
c    15    .910555D+02  -.408776D+04   .564468D-05  -.588222D-05
c    20   -.488032D+05   .111855D+06  -.750186D-07   .190624D-06
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MKLVNB'
      write ( *, '(a)' ) '  Test KLVNB'
      write ( *, '(a)' ) ' '
        WRITE(*,*)'   x        ber x           bei x',
     &            '           ker x           kei x'
        WRITE(*,*)'--------------------------------',
     &            '--------------------------------------'

      do i = 0, 4
        x = dble ( 5 * i )
        CALL KLVNB(X,BER,BEI,GER,GEI,DER,DEI,HER,HEI)
        WRITE(*,10)X,BER,BEI,GER,GEI
        WRITE(*,*)
        WRITE(*,*)'   x        ber''x           bei''x',
     &            '           ker''x           kei''x'
        WRITE(*,*)'--------------------------------',
     &            '--------------------------------------'
        WRITE(*,10)X,DER,DEI,HER,HEI
      end do

10      FORMAT(1X,F5.1,4D16.6)

      return
        END
      subroutine MKLVNZO ( )

c*********************************************************************72
c
cc MKLVNZO tests KLVNZO.
c
c  Modified:
c
c    11 April 2012
c
c
c       Purpose: This program computes the first NT zeros of Kelvin 
c                functions and their derivatives using subroutine
c                KLVNZO
c       Input :  NT --- Total number of zeros
c       Example: NT = 5
c
c           Zeros of Kelvin functions ber x, bei x, ker x and kei x
c
c        m       ber x          bei x          ker x          kei x
c       ---------------------------------------------------------------
c        1     2.84891782     5.02622395     1.71854296     3.91466761
c        2     7.23882945     9.45540630     6.12727913     8.34422506
c        3    11.67396355    13.89348785    10.56294271    12.78255715
c        4    16.11356383    18.33398346    15.00268812    17.22314372
c        5    20.55463158    22.77543929    19.44381663    21.66464214
c
c          Zeros of Kelvin Functions ber'x, bei'x, ker'x and kei'x
c
c        m       ber'x          bei'x          ker'x          kei'x
c       ---------------------------------------------------------------
c        1     6.03871081     3.77267330     2.66583979     4.93181194
c        2    10.51364251     8.28098785     7.17212212     9.40405458
c        3    14.96844542    12.74214752    11.63218639    13.85826916
c        4    19.41757493    17.19343175    16.08312025    18.30717294
c        5    23.86430432    21.64114394    20.53067845    22.75379258
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION R1(50),R2(50),R3(50),R4(50),R5(50),R6(50),
     &            R7(50),R8(50)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MKLVNZO'
      write ( *, '(a)' ) '  Test KLVNZO'
      write ( *, '(a)' ) ' '

      nt = 5
        WRITE(*,*)
        WRITE(*,*)'   m       ber x          bei x          ker x',
     &            '          kei x'
        WRITE(*,*)' ------------------------------------------------',
     &            '---------------'
        CALL KLVNZO(NT,1,R1)
        CALL KLVNZO(NT,2,R2)
        CALL KLVNZO(NT,3,R3)
        CALL KLVNZO(NT,4,R4)
        DO 10 L=1,NT
10         WRITE(*,20)L,R1(L),R2(L),R3(L),R4(L)
        CALL KLVNZO(NT,5,R5)
        CALL KLVNZO(NT,6,R6)
        CALL KLVNZO(NT,7,R7)
        CALL KLVNZO(NT,8,R8)
        WRITE(*,*)
        WRITE(*,30)
        WRITE(*,*)
        WRITE(*,*)'   m       ber''x          bei''x          ker''x',
     &            '          kei''x'
        WRITE(*,*)' ------------------------------------------------',
     &            '---------------'
        DO 15 L=1,NT
15         WRITE(*,20)L,R5(L),R6(L),R7(L),R8(L)
        CLOSE (01)
20      FORMAT(1X,I3,1X,F14.8,1X,F14.8,1X,F14.8,1X,F14.8)
25      FORMAT(4X,'Zeros of Kelvin functions ber x, bei x,'
     &        ,' ker x and kei x')
30      FORMAT(4X,'Zeros of Kelvin functions ber''x, bei''x,'
     &        ,' ker''x and kei''x')
35      FORMAT(1X,'NT is the number of the zeros')

      return
        END
      subroutine MLAGZO ( )

c*********************************************************************72
c
cc MLAGZO tests LAGZO.
c
c  Modified:
c
c    07 April 2012
c
c
c       Purpose : This program computes the zeros of Laguerre
c                 polynomial Ln(x) in the interval [0,ì] and the
c                 corresponding weighting coefficients for Gauss-
c                 Laguerre integration using subroutine LAGZO
c       Input :   n    --- Order of the Laguerre polynomial
c                 X(n) --- Zeros of the Laguerre polynomial
c                 W(n) --- Corresponding weighting coefficients
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION X(100),W(100)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MLAGZO'
      write ( *, '(a)' ) '  Test LAGZO'
      write ( *, '(a)' ) ' '

      n = 9

        CALL LAGZO(N,X,W)
        WRITE(*,*)'  Nodes and weights for Gauss-Lagurre integration'
        WRITE(*,*)
        WRITE(*,*)'  i             xi                      Wi'
        WRITE(*,*)' -----------------------------------------',
     &            '------------'
        DO 10 J=1,N
10         WRITE(*,30)J,X(J),W(J)
20      FORMAT(1X,'n =',I3)
30      FORMAT(1X,I3,3X,D22.13,3X,D22.13)

      return
        END
      subroutine MLAMN ( )

c*********************************************************************72
c
cc MLAMN tests LAMN.
c
c  Modified:
c
c    06 April 2012
c
c
c       Purpose: This program computes the lambda functions 
c                and their derivatives using subroutine
c                LAMN
c       Input:   x --- Argument of lambda function
c                n --- Order of lambda function
c                      ( n = 0,1,..., n ó 250 )
c       Output:  BL(n) --- Lambda function of order n
c                DL(n) --- Derivative of lambda function
c       Example: Nmax = 5,  x = 10.00
c
c                 n       lambda(x)        lambda'(x)
c                ---------------------------------------
c                 0    -.24593576D+00    -.43472746D-01
c                 1     .86945492D-02    -.50926063D-01
c                 2     .20370425D-01    -.46703503D-02
c                 3     .28022102D-02     .10540929D-01
c                 4    -.84327431D-02     .89879627D-02
c                 5    -.89879627D-02     .55521954D-03
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION BL(0:250),DL(0:250)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MLAMN'
      write ( *, '(a)' ) '  Test LAMN'
      write ( *, '(a)' ) ' '

      n = 5
      x = 10.0D+00

        IF (N.LE.10) THEN
           NS=1
        ELSE
           ns = 5
        ENDIF
        CALL LAMN(N,X,NM,BL,DL)
        WRITE(*,*) '  n       lambda(x)        lambda''(x)'
        WRITE(*,*)' ---------------------------------------'
        DO K=0,NM,NS
           WRITE(*,20)K,BL(K),DL(K)
        end do
15      FORMAT(1X,3HN =,I4,6X,3Hx =,F8.2)
20      FORMAT(1X,I3,2D18.8)

        return
        END
        subroutine MLAMV ( )

c*********************************************************************72
c
cc MLAMV tests LAMV.
c
c  Modified:
c
c    07 April 2012
c
c
c       Purpose: This program computes the lambda functions 
c                for an arbitrary order, and their derivative 
c                using subroutine LAMV
c       Input :  x --- Argument of lambda function
c                v --- Order of lambda function
c                      ( v = n+v0, 0 ó n ó 250, 0 ó v0 < 1 )
c       Output:  VL(n) --- Lambda function of order n+v0
c                DL(n) --- Derivative of lambda function 
c       Example: x = 10.0
c
c                   v         Lambda(x)        Lambda'(x)
c                 ------------------------------------------
c                  0.25    -.12510515D+00    -.78558916D-01
c                  0.50    -.54402111D-01    -.78466942D-01
c                  0.75    -.13657787D-01    -.66234027D-01
c                  1.00     .86945492D-02    -.50926063D-01
c                  1.25     .19639729D-01    -.36186221D-01
c                  1.50     .23540083D-01    -.23382658D-01
c                  1.75     .23181910D-01    -.12893894D-01
c                  2.00     .20370425D-01    -.46703503D-02
c                  2.25     .16283799D-01     .15101684D-02
c                  2.50     .11691329D-01     .59243767D-02
c
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VL(0:250),DL(0:250)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MLAMV'
      write ( *, '(a)' ) '  Test LAMV'
      write ( *, '(a)' ) ' '

      x = 10.0D+00

      IF (V.LE.8) THEN
         NS=1
      ELSE
         ns = 5
      ENDIF

      WRITE(*,*) '   v         Lambda(x)        Lambda''(x)'
      WRITE(*,*)'-------------------------------------------'

      do i = 7, 10
        v = dble ( i ) / 4.0D+00
      CALL LAMV(V,X,VM,VL,DL)
      NM=INT(VM)
      V0=VM-NM
      DO K=0,NM,NS
       VK=K+V0
           WRITE(*,15)VK,VL(K),DL(K)
      end do
      end do

15      FORMAT(1X,F6.2,2D18.8)
20      FORMAT(1X,'v =',F6.2,'    ','x =',F8.2)

      return
      END
      subroutine MLEGZO ( )

c*********************************************************************72
c
cc MLEGZO tests LEGZO.
c
c  Modified:
c
c    07 April 2012
c
c
c       Purpose : This program computes the zeros of Legendre 
c                 polynomial Pn(x) in the interval [-1,1] and the
c                 corresponding weighting coefficients for Gauss-
c                 Legendre integration using subroutine LEGZO
c       Input :   n    --- Order of the Legendre polynomial
c       Output:   X(n) --- Zeros of the Legendre polynomial
c                 W(n) --- Corresponding weighting coefficients
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION X(120),W(120)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MLEGZO'
      write ( *, '(a)' ) '  Test LEGZO'
      write ( *, '(a)' ) ' '

      n = 9

        CALL LEGZO(N,X,W)
        WRITE(*,*)'  Nodes and weights for Gauss-Legendre integration'
        WRITE(*,*)
        WRITE(*,*)'  i              xi                   Wi'
        WRITE(*,*)' ------------------------------------------------'
        DO 10 I=1,N
10         WRITE(*,20)I,X(I),W(I)
15      FORMAT(1X,'n =',I3)
20      FORMAT(1X,I3,1X,F22.13,D22.13)
      return
        END
      subroutine MLGAMA ( )

c*********************************************************************72
c
cc MLGAMA tests LGAMA.
c
c  Modified:
c
c    06 April 2012
c
c
c       Purpose: This program computes the gamma function
c                â(x) for x > 0 using subroutine LGAMA
c       Examples:
c                  x           â(x)
c                -------------------------
c                 0.5     .1772453851D+01
c                 2.5     .1329340388D+01
c                 5.0     .2400000000D+02
c                 7.5     .1871254306D+04
c                10.0     .3628800000D+06
c
c
        IMPLICIT DOUBLE PRECISION (G,X)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MLGAMA'
      write ( *, '(a)' ) '  Test LGAMA'
      write ( *, '(a)' ) ' '

        WRITE(*,*)'   x           â(x)'
        WRITE(*,*)' -------------------------'
        DO 10 L=0,20,5
           X=0.5D0*L
           IF (L.EQ.0) X=0.5
           CALL LGAMA(1,X,GL)
           WRITE(*,20)X,GL
10      CONTINUE
20      FORMAT(1X,F5.1,D20.10)

      return
        END
      subroutine MLPMN ( )

c*********************************************************************72
c
cc MLPMN tests LPMN.
c
c  Modified:
c
c    07 April 2012
c
c
c       Purpose: This program computes the associated Legendre 
c                functions Pmn(x) and their derivatives Pmn'(x) 
c                using subroutine LPMN
c       Input :  x --- Argument of Pmn(x)
c                m --- Order of Pmn(x),  m = 0,1,2,...,n
c                n --- Degree of Pmn(x), n = 0,1,2,...,N
c       Output:  PM(m,n) --- Pmn(x)
c                PD(m,n) --- Pmn'(x)
c       Example: x = 0.50
c          Pmn(x):
c          m\n        1            2            3            4
c         --------------------------------------------------------
c           0      .500000     -.125000     -.437500     -.289063
c           1     -.866025    -1.299038     -.324760     1.353165
c           2      .000000     2.250000     5.625000     4.218750
c           3      .000000      .000000    -9.742786   -34.099750
c           4      .000000      .000000      .000000    59.062500
c
c          Pmn'(x):
c          m\n        1            2            3            4
c         --------------------------------------------------------
c           0     1.000000     1.500000      .375000    -1.562500
c           1      .577350    -1.732051    -6.278684    -5.773503
c           2      .000000    -3.000000     3.750000    33.750000
c           3      .000000      .000000    19.485572      .000000
c           4      .000000      .000000      .000000  -157.500000
c
c
        IMPLICIT DOUBLE PRECISION (P,X)
        DIMENSION PM(0:100,0:100),PD(0:100,0:100)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MLPMN'
      write ( *, '(a)' ) '  Test LPMN'
      write ( *, '(a)' ) ' '

      x = 0.50D+00
      n = 4

        WRITE(*,*)'  m     n      x          Pmn(x)         Pmn''(x)'
        WRITE(*,*)' ---------------------------------------------------'
        do m = 0, 4
          CALL LPMN(100,M,N,X,PM,PD)
          DO 15 J=0,N
             WRITE(*,10)M,J,X,PM(M,J),PD(M,J)
15        CONTINUE
        end do
10      FORMAT(1X,I3,3X,I3,3X,F5.1,2E17.8)
      return
        END
      subroutine MLPMNS ( )

c*********************************************************************72
c
cc MLPMNS tests LPMNS.
c
c  Modified:
c
c    11 April 2012
c
c
c       Purpose: This program computes the associated Legendre 
c                functions Pmn(x) and their derivatives Pmn'(x) 
c                for a given order using subroutine LPMNS
c       Input :  x --- Argument of Pmn(x)
c                m --- Order of Pmn(x),  m = 0,1,2,...,n
c                n --- Degree of Pmn(x), n = 0,1,2,...,N
c       Output:  PM(n) --- Pmn(x)
c                PD(n) --- Pmn'(x)
c       Examples:
c                m = 1,  N = 5,  x = .5
c                n        Pmn(x)           Pmn'(x)
c               -------------------------------------
c                0    .00000000D+00    .00000000D+00
c                1    .86602540D+00   -.57735027D+00
c                2    .12990381D+01    .17320508D+01
c                3    .32475953D+00    .62786842D+01
c                4   -.13531647D+01    .57735027D+01
c                5   -.19282597D+01   -.43977853D+01
c
c                m = 2,  N = 6,  x = 2.5
c                n        Pmn(x)           Pmn'(x)
c               -------------------------------------
c                0    .00000000D+00    .00000000D+00
c                1    .00000000D+00    .00000000D+00
c                2    .15750000D+02    .15000000D+02
c                3    .19687500D+03    .26625000D+03
c                4    .16832813D+04    .29812500D+04
c                5    .12230859D+05    .26876719D+05
c                6    .81141416D+05    .21319512D+06
c
c
        IMPLICIT DOUBLE PRECISION (P,X,Y)
        DIMENSION PM(0:200),PD(0:200)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MLPMNS'
      write ( *, '(a)' ) '  Test LPMNS'
      write ( *, '(a)' ) ' '
        WRITE(*,*)'  n        Pmn(x)           Pmn''(x)    '
        WRITE(*,*)' -------------------------------------'

      m = 1
      n = 5
      x = 0.5D+00

        CALL LPMNS(M,N,X,PM,PD)

        DO J=0,N
        WRITE(*,20)J,PM(J),PD(J)
        end do

      m = 2
      n = 6
      x = 2.5D+00

        CALL LPMNS(M,N,X,PM,PD)

        DO J=0,N
        WRITE(*,20)J,PM(J),PD(J)
        end do

20      FORMAT(1X,I3,2D17.8)
30      FORMAT(1X,'m =',I2,',  ','n =',I2,',  ','x =',F5.1)
      return
        END
      subroutine MLPMV ( )

c*********************************************************************72
c
cc MLPMV tests LPMV.
c
c  Modified:
c
c    08 April 2012
c
c
c       Purpose: This program computes the associated Legendre 
c                function Pmv(x) with an integer order and an
c                arbitrary nonnegative degree using subroutine 
c                LPMV
c       Input :  x   --- Argument of Pm(x)  ( -1 ó x ó 1 )
c                m   --- Order of Pmv(x)
c                v   --- Degree of Pmv(x)
c       Output:  PMV --- Pmv(x)
c       Example:    m = 4,  x = 0.5
c                    v          Pmv(x)
c                 -----------------------
c                   1.5       .46218726
c                   1.6       .48103143
c                   1.7       .45031429
c                   1.8       .36216902
c                   1.9       .21206446
c                   2.0       .00000000
c                   2.5     -1.51996235
c
c
        IMPLICIT DOUBLE PRECISION (P,V,X)

      integer test_num
      parameter ( test_num = 7 )

      integer test
      double precision v
      double precision v_test(test_num)

      save v_test

      data v_test / 
     &  1.5D+00, 1.6D+00, 1.7D+00, 1.8D+00, 1.9D+00,
     &  2.0D+00, 2.5D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MLPMN'
      write ( *, '(a)' ) '  Test LPMN'
      write ( *, '(a)' ) ' '

      m = 4
      x = 0.5D+00

        WRITE(*,*)'     v        Pmv(x)'
        WRITE(*,*)'  -----------------------'
      do test = 1, test_num
        v = v_test(test)
        CALL LPMV(V,M,X,PMV)
        WRITE(*,10)V,PMV
      end do
10      FORMAT(3X,F5.1,E16.8)
20      FORMAT(3X,'m =',I2,',    ','x =',F6.2)
      return
        END
      subroutine MLPN ( )

c*********************************************************************72
c
cc MLPN tests LPN.
c
c  Modified:
c
c    06 April 2012
c
c
c       Purpose: This program computes the Legendre polynomials 
c                Pn(x) and their derivatives Pn'(x) using
c              subroutine LPN
c       Input :  x --- Argument of Pn(x)
c                n --- Degree of Pn(x) ( n = 0,1,...)
c       Output:  PN(n) --- Pn(x)
c                PD(n) --- Pn'(x)
c       Example:    x = 0.5
c                  n          Pn(x)            Pn'(x)
c                ---------------------------------------
c                  0       1.00000000        .00000000
c                  1        .50000000       1.00000000
c                  2       -.12500000       1.50000000
c                  3       -.43750000        .37500000
c                  4       -.28906250      -1.56250000
c                  5        .08984375      -2.22656250
c
c
        DOUBLE PRECISION PN,PD,X
        DIMENSION PN(0:100),PD(0:100)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MLPN'
      write ( *, '(a)' ) '  Test LPN'
      write ( *, '(a)' ) ' '

      n = 5
      x = 0.5D+00

        CALL LPN(N,X,PN,PD)
        WRITE(*,*)'  n         Pn(x)           Pn''(x)'
        WRITE(*,*)'---------------------------------------'
        DO 10 K=0,N
10         WRITE(*,20)K,PN(K),PD(K)
20      FORMAT(1X,I3,2E17.8)
30      FORMAT(3X,'x =',F5.1)

      return
        END
      subroutine MLPNI ( )

c*********************************************************************72
c
cc MLPNI tests LPNI.
c
c  Modified:
c
c    10 April 2012
c
c
c       Purpose: This program computes the Legendre polynomials 
c                Pn(x), Pn'(x) and the integral of Pn(t) from 0 
c                to x using subroutine LPNI
c       Input :  x --- Argument of Pn(x)
c                n --- Degree of Pn(x) ( n = 0,1,... )
c       Output:  PN(n) --- Pn(x)
c                PD(n) --- Pn'(x)
c                PL(n) --- Integral of Pn(t) from 0 to x
c       Example: x = 0.50
c                n       Pn(x)         Pn'(x)        Pn(t)dt
c               ---------------------------------------------
c                0    1.00000000     .00000000     .50000000
c                1     .50000000    1.00000000     .12500000
c                2    -.12500000    1.50000000    -.18750000
c                3    -.43750000     .37500000    -.14843750
c                4    -.28906250   -1.56250000     .05859375
c                5     .08984375   -2.22656250     .11816406
c
c
        DOUBLE PRECISION PN,PD,PL,X
        DIMENSION PN(0:100),PD(0:100),PL(0:100)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MLPNI'
      write ( *, '(a)' ) '  Test LPNI'
      write ( *, '(a)' ) ' '

      x = 0.5D+00
      n = 5

        WRITE(*,*)'  n        Pn(x)          Pn''(x)         Pn(t)dt'
        WRITE(*,*)' ---------------------------------------------------'
        CALL LPNI(N,X,PN,PD,PL)
        DO 10 K=0,N
10         WRITE(*,20)K,PN(K),PD(K),PL(K)
20      FORMAT(1X,I3,3E16.8)
30      FORMAT(3X,'x =',F5.2)
      return
        END
      subroutine MLQMN ( )

c*********************************************************************72
c
cc MLQMN tests LQMN.
c
c  Modified:
c
c    10 April 2012
c
c
c       Purpose: This program computes the associated Legendre  
c                functions Qmn(x) and their derivatives Qmn'(x) using
c              subroutine LQMN
c       Input :  x --- Argument of Qmn(x) 
c                m --- Order of Qmn(x)  ( m = 0,1,2,úúú )
c                n --- Degree of Qmn(x) ( n = 0,1,2,úúú )
c       Output:  QM(m,n) --- Qmn(x)
c                QD(m,n) --- Qmn'(x)
c       Examples:
c
c       Qmn(x):  x = 0.5
c       n\m      0           1           2           3           4
c       ---------------------------------------------------------------
c        0     .549306   -1.154701    1.333333   -5.388603   26.666667
c        1    -.725347   -1.053063    2.666667   -6.158403   32.000000
c        2    -.818663     .729806    4.069272  -12.316806   42.666667
c        3    -.198655    2.491853    -.493486  -23.778868   85.333333
c        4     .440175    1.934087  -11.036781   -9.325204  186.818394
c
c       Qmn'(x): x = 0.5
c       n\m      0           1           2           3           4
c       ---------------------------------------------------------------
c        0    1.333333    -.769800    4.444444  -20.014809  145.777778
c        1    1.215973   -2.377159    3.555556  -24.633611  156.444444
c        2    -.842707   -5.185328    8.796526  -24.633611  199.111111
c        3   -2.877344   -1.091406   28.115454  -50.976710  227.555556
c        4   -2.233291   11.454786   25.483527 -197.068892  412.039838
c
c       Qmn(x): x = 2.0
c       n\m      0           1           2           3           4
c       ---------------------------------------------------------------
c        0     .549306    -.577350    1.333333   -5.003702   26.666667
c        1     .098612    -.203274     .666667   -3.079201   18.666667
c        2     .021184    -.064946     .277089   -1.539601   10.666667
c        3     .004871    -.019817     .104220    -.679543    5.333333
c        4     .001161    -.005887     .036816    -.276005    2.427640
c
c       Qmn'(x): x = 2.0
c       n\m      0           1           2           3           4
c       ---------------------------------------------------------------
c        0    -.333333     .384900   -1.111111    5.388603  -36.444444
c        1    -.117361     .249384    -.888889    4.618802  -32.000000
c        2    -.037496     .116680    -.519437    3.079201  -23.111111
c        3    -.011442     .046960    -.253375    1.720114  -14.222222
c        4    -.003399     .017331    -.110263     .849589   -7.748516
c
c
        IMPLICIT DOUBLE PRECISION (Q,X)
        DIMENSION QM(0:100,0:100),QD(0:100,0:100)
      double precision x_test(2)

      save x_test

      data x_test / 0.5D+00, 2.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MLQMN'
      write ( *, '(a)' ) '  Test LQMN'

        WRITE(*,*)
        WRITE(*,*)'  m     n      x          Qmn(x)         Qmn''(x)'
        WRITE(*,*)' ---------------------------------------------------'

      n = 4
      do i = 1, 2
        x = x_test(i)
        do m = 0, 4
          CALL LQMN(100,M,N,X,QM,QD)
          DO J=0,N
             WRITE(*,10)M,J,X,QM(M,J),QD(M,J)
          end do
        end do
      end do

10      FORMAT(1X,I3,3X,I3,3X,F5.1,2D17.8)

      return
        END
      subroutine MLQMNS ( )

c*********************************************************************72
c
cc MLQMNS tests LQMNS.
c
c  Modified:
c
c    10 April 2012
c
c
c       Purpose: This program computes the associated Legendre 
c                functions Qmn(x) and their derivatives Qmn'(x) 
c                for a given order using subroutine LQMNS
c       Input :  x --- Argument of Qmn(x)
c                m --- Order of Qmn(x),  m = 0,1,2,...
c                n --- Degree of Qmn(x), n = 0,1,2,...
c       Output:  QM(n) --- Qmn(x)
c                QD(n) --- Qmn'(x)
c       Examples:
c                m = 1,  N = 5,  x = .5
c                n        Qmn(x)           Qmn'(x)
c               -------------------------------------
c                0    .11547005D+01    .76980036D+00
c                1    .10530633D+01    .23771592D+01
c                2   -.72980606D+00    .51853281D+01
c                3   -.24918526D+01    .10914062D+01
c                4   -.19340866D+01   -.11454786D+02
c                5    .93896830D+00   -.18602587D+02
c
c                m = 2,  N = 5,  x = 2.5
c                n        Qmn(x)           Qmn'(x)
c               -------------------------------------
c                0    .95238095D+00   -.52607710D+00
c                1    .38095238D+00   -.36281179D+00
c                2    .12485160D+00   -.17134314D+00
c                3    .36835513D-01   -.66284127D-01
c                4    .10181730D-01   -.22703958D-01
c                5    .26919481D-02   -.71662396D-02
c
c
        IMPLICIT DOUBLE PRECISION (Q,X,Y)
        DIMENSION QM(0:200),QD(0:200)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MLQMNS'
      write ( *, '(a)' ) '  Test LQMNS'

        WRITE(*,*)' '
        WRITE(*,*)'  n        Qmn(x)           Qmn''(x)'
        WRITE(*,*)' -------------------------------------'
        m = 1
        n = 5
        x = 0.5D+00
        WRITE(*,30)M,N,X
        CALL LQMNS(M,N,X,QM,QD)
        DO J=0,N
          WRITE(*,20)J,QM(J),QD(J)
      end do

        m = 2
        n = 5
        x = 2.5D+00
        WRITE(*,30)M,N,X
        CALL LQMNS(M,N,X,QM,QD)
        DO J=0,N
          WRITE(*,20)J,QM(J),QD(J)
      end do

20      FORMAT(1X,I3,2D17.8)
30      FORMAT(1X,'m =',I2,',  ','n =',I2,',  ','x =',F5.1)

      return
        END
      subroutine MLQNA ( )

c*********************************************************************72
c
cc MLQNA tests LQNA.
c
c  Modified:
c
c    10 April 2012
c
c
c       Purpose: This program computes the Legendre functions
c                Qn(x) and Qn'(x) using subroutine LQNA
c       Input :  x  --- Argument of Qn(x)  ( -1 ó x ó 1 )
c                n  --- Degree of Qn(x)  ( n = 0,1,úúú )
c       Output:  QN(n) --- Qn(x)
c                QD(n) --- Qn'(x)
c       Example:  x = 0.50
c                 n        Qn(x)         Qn'(x)
c                ---------------------------------
c                 0      .54930614     1.33333333
c                 1     -.72534693     1.21597281
c                 2     -.81866327     -.84270745
c                 3     -.19865477    -2.87734353
c                 4      .44017453    -2.23329085
c                 5      .55508089     1.08422720
c
c
        DOUBLE PRECISION QN,QD,X
        DIMENSION QN(0:100),QD(0:100)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MLQNA'
      write ( *, '(a)' ) '  Test LQNA'

      x = 0.5D+00
      n = 5

        CALL LQNA(N,X,QN,QD)
        WRITE(*,*)'  n        Qn(x)         Qn''(x)'
        WRITE(*,*)' ---------------------------------'
        DO 10 K=0,N
10         WRITE(*,20)K,QN(K),QD(K)
20      FORMAT(1X,I3,2F15.8)
30      FORMAT(3X,'x =',F5.2)
      return
        END
      subroutine MLQNB ( )

c*********************************************************************72
c
cc MLQNB tests LQNB.
c
c  Modified:
c
c    09 April 2012
c
c
c       Purpose: This program computes the Legendre functions Qn(x) 
c                and Qn'(x) using subroutine LQNB
c       Input :  x  --- Argument of Qn(x)
c                n  --- Degree of Qn(x)  ( n = 0,1,úúú)
c       Output:  QN(n) --- Qn(x)
c                QD(n) --- Qn'(x)
c       Examples:     x1 = 0.50,    x2 = 2.50
c
c       n      Qn(x1)        Qn'(x1)       Qn(x2)          Qn'(x2)
c     ----------------------------------------------------------------
c       0     .54930614    1.33333333   .42364893D+00  -.19047619D+00
c       1    -.72534693    1.21597281   .59122325D-01  -.52541546D-01
c       2    -.81866327    -.84270745   .98842555D-02  -.13109214D-01
c       3    -.19865477   -2.87734353   .17695141D-02  -.31202687D-02
c       4     .44017453   -2.23329085   .32843271D-03  -.72261513D-03
c       5     .55508089    1.08422720   .62335892D-04  -.16437427D-03
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION QN(0:100),QD(0:100)
      double precision x_test(2)

      save x_test

      data x_test / 0.5D+00, 2.5D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MLQNB'
      write ( *, '(a)' ) '  Test LQNB'

      n = 5
      do i = 1, 2
        x = x_test(i)
        WRITE(*,*)
        WRITE(*,*)'  n      x          Qn(x)           Qn''(x)'
        WRITE(*,*)'-------------------------------------------'
        CALL LQNB(N,X,QN,QD)
        DO K=0,N
           IF (X.LE.1.0) THEN
              WRITE(*,20)K,X,QN(K),QD(K)
           ELSE
              WRITE(*,30)K,X,QN(K),QD(K)
           ENDIF
        end do

      end do

20      FORMAT(1X,I3,F8.5,2F17.8)
30      FORMAT(1X,I3,F8.5,2D17.8)

      return
        END
      subroutine MMTU0 ( )

c*********************************************************************72
c
cc MMTU0 tests MTU0.
c
c  Modified:
c
c    06 April 2012
c
c
c       Purpose: This program computes Mathieu functions cem(x,q), 
c                sem(x,q) and their derivatives using subroutine
c                MTU0 ( q ò 0 )
c       Input :  KF  --- Function code
c                        KF=1 for computing cem(x,q) and cem'(x,q)
c                        KF=2 for computing sem(x,q) and sem'(x,q)
c                m   --- Order of Mathieu functions
c                q   --- Parameter of Mathieu functions
c                x   --- Argument of Mathieu functions (in degrees)
c       Output:  CSF --- cem(x,q) or sem(x,q)
c                CSD --- cem'x,q) or sem'x,q)
c       Example: x = 40
c           m     q    cem(x,q)   cem'(x,q)    sem(x,q)  sem'(x,q)
c          --------------------------------------------------------
c           0    5.0   .3025683    .9470247
c           1    5.0   .7669652   1.2873097    .2988052   .9606824
c           2    5.0   .9102723   -.3463855    .7549264  1.4743128
c           5    5.0  -.9810931   -.6328576    .1694850 -4.8676455
c           0   25.0   .0515371    .3823737
c           1   25.0   .2074402   1.2646301    .0515365   .3823777
c           2   25.0  -.5297051  -2.4292679    .2074275  1.2646996
c           5   25.0   .7507159  -3.9047012   1.1881232   .3258081
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      integer i
      integer j
      integer m
      integer m_test(4)
      double precision q
      double precision q_test(2)

      save m_test
      save q_test

      data m_test / 0, 1, 2, 5 /
      data q_test / 5.0D+00, 25.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MMTU0'
      write ( *, '(a)' ) '  Test MTU0'
      write ( *, '(a)' ) ' '

      x = 40.0D+00

      WRITE(*,*) '  m   q       x       cem(x,q)        cem''(x,q)' 
     &  // '       sem(x,q)        sem''(x,q)'
      WRITE(*,*) '--------------------------------------'
     &  // '--------------------------------------'

      do j = 1, 2
        q = q_test(j)
        do i = 1, 4
          m = m_test(i)
          kf = 1
          CALL MTU0(KF,M,Q,X,CSF,CSD)
          kf = 2
          CALL MTU0(KF,M,Q,X,SSF,SSD)
          WRITE(*,20) m, q, X, CSF, CSD, SSF, SSD
        end do
      end do

20      FORMAT(2X,i2,2x,f5.2,2x,F5.1,2F16.9,2F16.9)
      return
        END
      subroutine MMTU12 ( )

c*********************************************************************72
c
cc MMTU12 tests MTU12.
c
c  Modified:
c
c    10 April 2012
c
c
c       Purpose: This program computes the modified Mathieu functions 
c                of the first and second kinds, Mcm(1)(2)(x,q) and 
c                Msm(1)(2)(x,q), and their derivatives using 
c              subroutine MTU12
c       Input:   KF --- Function code
c                       KF=1 for computing Mcm(x,q)
c                       KF=2 for computing Msm(x,q)
c                KC --- Function Code
c                       KC=1 for computing Mcm(1)(x,q) and Mcm(1)'(x,q)
c                            or Msm(1)(x,q) and Msm(1)'(x,q)
c                       KC=2 for computing Mcm(2)(x,q) and Mcm(2)'(x,q)
c                            or Msm(2)(x,q) and Msm(2)'(x,q)
c                       KC=3 for both modified Mathieu functions of the
c                            first and second kinds, and their
c                            derivatives
c                m  --- Order of Mathieu functions
c                q  --- Parameter of Mathieu functions
c                x  --- Argument of Mathieu functions
c       Output:  F1R --- Mcm(1)(x,q) or Msm(1)(x,q)
c                D1R --- Derivative of Mcm(1)(x,q) or Msm(1)(x,q)
c                F2R --- Mcm(2)(x,q) or Msm(2)(x,q)
c                D2R --- Derivative of Mcm(2)(x,q) or Msm(2)(x,q)
c
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MMTU12'
      write ( *, '(a)' ) '  Test MTU12'
      write ( *, '(a)' ) ' '

      kc = 3

      kf = 1
      m = 0
      q = 5.0D+00
      x = 40.0D+00

      CALL MTU12(KF,KC,M,Q,X,F1R,D1R,F2R,D2R)
      WRITE(*,*)
       WRITE(*,*)'   x      Mcm(1)(x,q)    Mcm(1)''(x,q)',
     &               '    Mcm(2)(x,q)     Mcm(2)''(x,q)'
      WRITE(*,*)' --------------------------------------',
     &            '-------------------------------'
      WRITE(*,20)X,F1R,D1R,F2R,D2R
      WRITE(*,*)
      WRITE(*,30)F1R*D2R-F2R*D1R,.63661977236758D0
      WRITE(*,40)

      kf = 2
      m = 0
      q = 5.0D+00
      x = 40.0D+00

      CALL MTU12(KF,KC,M,Q,X,F1R,D1R,F2R,D2R)
      WRITE(*,*)
       WRITE(*,*)'   x      Msm(1)(x,q)    Msm(1)''(x,q)',
     &               '    Msm(2)(x,q)     Msm(2)''(x,q)'
      WRITE(*,*)' --------------------------------------',
     &            '-------------------------------'
      WRITE(*,20)X,F1R,D1R,F2R,D2R
      WRITE(*,*)
      WRITE(*,30)F1R*D2R-F2R*D1R,.63661977236758D0
      WRITE(*,40)

10      FORMAT(1X,4HKF =,I2,',  ',3Hm =,I3,',  ',
     &         3Hq =,F5.1,',  ',3Hx =,F5.1)
20      FORMAT(1X,F5.1,4D16.8)
30      FORMAT(1X,'WRONSKIAN=',E16.8,3X,'should equal   2/PI=',E16.8)
40      FORMAT(1X,/1X,'Caution: This check is not accurate if it ',
     &        'involves',/1X,'         the subtraction of two ',
     &        'similar numbers')
      return
      END
      subroutine MOTHPL ( )

c*********************************************************************72
c
cc MOTHPL tests OTHPL.
c
c  Modified:
c
c    06 April 2012
c
c
c       Purpose: This program computes orthogonal polynomials: 
c                Tn(x) or Un(x) or Ln(x) or Hn(x), and their
c                derivatives using subroutine OTHPL
c       Input :  KF --- Function code
c                       KF=1 for Chebyshev polynomial Tn(x)
c                       KF=2 for Chebyshev polynomial Un(x)
c                       KF=3 for Laguerre polynomial Ln(x)
c                       KF=4 for Hermite polynomial Hn(x)
c                n ---  Order of orthogonal polynomials
c                x ---  Argument
c       Output:  PL(n) --- Tn(x) or Un(x) or Ln(x) or Hn(x)
c                DPL(n)--- Tn'(x) or Un'(x) or Ln'(x) or Hn'(x)
c                          n = 0,1,2,...,N ( N ó 100 )
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION PL(0:100),DPL(0:100)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MPSI'
      write ( *, '(a)' ) '  Test PSI'
      write ( *, '(a)' ) ' '

      x = 0.25D+00
      n = 5
      do kf = 1, 4

        CALL OTHPL(KF,N,X,PL,DPL)
        IF (KF.EQ.1) WRITE(*,*)'  n          Tn(x)            Tn''(x)'
        IF (KF.EQ.2) WRITE(*,*)'  n          Un(x)            Un''(x)'
        IF (KF.EQ.3) WRITE(*,*)'  n          Ln(x)            Ln''(x)'
        IF (KF.EQ.4) WRITE(*,*)'  n          Hn(x)            Hn''(x)'
        WRITE(*,*)'-----------------------------------------'
        DO 10 K=0,N
10         WRITE(*,30)K,PL(K),DPL(K)
20      FORMAT(1X,3HKF=,I3,5X,5HNmax=,I3,5X,2Hx=,F6.3)
30      FORMAT(1X,I3,2D18.8)
      end do

      return
        END
      subroutine MPBDV ( )

c*********************************************************************72
c
cc MPBDV tests PBDV.
c
c  Modified:
c
c    06 April 2012
c
c
c       Purpose: This program computes the parabolic cylinder 
c                functions Dv(x) and their derivatives using
c              subroutine PBDV
c       Input:   x --- Argument of Dv(x)
c                v --- Order of Dv(x)
c       Output:  DV(na) --- Dn+v0(x)
c                DP(na) --- Dn+v0'(x)
c                ( na = |n|, n = int(v), v0 = v-n, |v0| < 1
c                  n = 0,ñ1,ñ2,úúú, |n| ó 100 )
c                PDF --- Dv(x)
c                PDD --- Dv'(x)
c       Example: v = 5.5,  x =10.0,  v0 = 0.5,  n = 0,1,...,5
c
c                  n+v0      Dv(x)           Dv'(x)
c                ---------------------------------------
c                  0.5   .43971930D-10  -.21767183D-09
c                  1.5   .43753148D-09  -.21216995D-08
c                  2.5   .43093569D-08  -.20452956D-07
c                  3.5   .41999741D-07  -.19491595D-06
c                  4.5   .40491466D-06  -.18355745D-05
c                  5.5   .38601477D-05  -.17073708D-04
c
c                Dv(x)= .38601477D-05,  Dv'(x)=-.17073708D-04
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION DV(0:100),DP(0:100)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MPBDV'
      write ( *, '(a)' ) '  Test PBDV'
      write ( *, '(a)' ) ' '

      v = 5.5D+00
      x = 10.0D+00

        NV=INT(V)
        V0=V-NV
        NA=ABS(NV)
        CALL PBDV(V,X,DV,DP,PDF,PDD)

        WRITE(*,*)'   v       Dv(x)           Dv''(x)'
        WRITE(*,*)'---------------------------------------'
        DO 10 K=0,NA
           VK=K*ISIGN(1,NV)+V0
10         WRITE(*,30)VK,DV(K),DP(K)
        WRITE(*,*)
        WRITE(*,40)V,PDF,PDD
20      FORMAT(1X,'v =',F6.2,',    ','x =',F6.2)
30      FORMAT(1X,F5.1,2D16.8)
40      FORMAT(1X,'v =',F5.1,',  Dv(x)=',D14.8,',   Dv''(x)=',D14.8)

      return
        END
      subroutine MPBVV ( )

c*********************************************************************72
c
cc MPBVV tests PBVV.
c
c  Modified:
c
c    06 April 2012
c
c
c       Purpose: This program computes the parabolic cylinder 
c                functions Vv(x) and Vv'(x) using subroutine
c                PBVV
c       Input:   x --- Argument of Vv(x)
c                v --- Order of Vv(x)
c       Output:  VV(na) --- Vv(x)
c                VP(na) --- Vv'(x)
c                ( na = |n|, v = n+v0, n = int(v), |v0| < 1
c                  n = 0,ñ1,ñ2,úúú, |n| ó 100 )
c                PVF --- Vv(x)
c                PVD --- Vv'(x)
c       Example: v = 5.5,  x =10.0,  v0 = 0.5,  n = 0,1,2,...,5
c
c                  n+v0      Vv(x)           Vv'(x)
c                ---------------------------------------
c                  0.5   .18522719D+10   .89761157D+10
c                  1.5   .19016268D+09   .90145854D+09
c                  2.5   .19741946D+08   .91452949D+08
c                  3.5   .20733667D+07   .93751130D+07
c                  4.5   .22038231D+06   .97145511D+06
c                  5.5   .23719356D+05   .10178553D+06
c
c                Vv(x)= .23719356D+05,  Vv'(x)= .10178553D+06
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION VV(0:100),VP(0:100)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MPBVV'
      write ( *, '(a)' ) '  Test PBVV'
      write ( *, '(a)' ) ' '

      v = 5.5D+00
      x = 10.0D+00

        NV=INT(V)
        V0=V-NV
        NA=ABS(NV)
        CALL PBVV(V,X,VV,VP,PVF,PVD)
        WRITE(*,*)'   v       Vv(x)           Vv''(x)'
        WRITE(*,*)'---------------------------------------'
        DO 10 K=0,NA
           VK=K*ISIGN(1,NV)+V0
10         WRITE(*,30)VK,VV(K),VP(K)
        WRITE(*,*)
        WRITE(*,40)V,PVF,PVD
20      FORMAT(1X,'v =',F6.2,',    ','x =',F6.2)
30      FORMAT(1X,F5.1,2D16.8)
40      FORMAT(1X,'v =',F5.1,',  Vv(x)=',D14.8,',   Vv''(x)=',D14.8)
      return
        END
      subroutine MPBWA ( )

c*********************************************************************72
c
cc MPBWA tests PBWA.
c
c  Modified:
c
c    06 April 2012
c
c
c       Purpose: This program computes the parabolic cylinder 
c                functions W(a,ñx) and their derivatives using
c              subroutine PBWA
c       Input  : a --- Parameter  ( 0 ó |a| ó 5 )
c                x --- Argument of W(a,ñx)  ( 0 ó |x| ó 5 )
c       Output : W1F --- W(a,x)
c                W1D --- W'(a,x)
c                W2F --- W(a,-x)
c                W2D --- W'(a,-x)
c       Example: x = 5.0
c                 a      W(a,x)     W'(a,x)    W(a,-x)   W'(a,-x)
c              ----------------------------------------------------
c                0.5   .1871153    .1915744  -.8556585   4.4682493
c                1.5  -.0215853    .0899870 -8.8586002  -9.3971967
c                0.0   .3009549   -.7148233   .6599634   1.7552224
c               -0.5  -.1934088  -1.3474400   .6448148   -.6781011
c               -1.5  -.5266539    .8219516  -.2822774  -1.4582283
c               -5.0   .0893618  -1.8118641   .5386084    .2698553
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      integer test_num
      parameter ( test_num = 6 )

      double precision a
      double precision a_test(test_num)
      integer test
      double precision x

      save a_test

      data a_test / 
     &  0.5D+00, 1.5D+00, 0.0D+00, -0.5D+00, -1.5D+00, -5.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MPBWA'
      write ( *, '(a)' ) '  Test PBWA'
      write ( *, '(a)' ) ' '
        WRITE(*,*)'   a       W(a,x)          W''(a,x)',
     &            '         W(a,-x)         W''(a,-x)'
        WRITE(*,*)' -----------------------------------------',
     &            '----------------------------'

      x = 5.0D+00
      do test = 1, test_num
        a = a_test(test)
        CALL PBWA(A,X,W1F,W1D,W2F,W2D)
        WRITE(*,20)A,W1F,W1D,W2F,W2D
      end do

10      FORMAT(1X,'a=',F5.1,3X,'x=',F5.1)
20      FORMAT(1X,F5.1,4D16.8)
      RETURN
        END
      subroutine MPSI ( )

c*********************************************************************72
c
cc MPSI tests PSI.
c
c  Modified:
c
c    06 April 2012
c
c
c       Purpose: This program computes the psi function
c                using subroutine PSI
c       Input :  x  --- Argument of psi(x)
c       Output:  PS --- psi(x)
c       Examples:
c                   x          Psi(x)
c                ------------------------
c                  .25      -4.227453533
c                  .50      -1.963510026
c                  .75      -1.085860880
c                 1.00       -.577215665
c                 1.25       -.227453533
c                 1.50        .036489974
c                 1.75        .247472454
c                 2.00        .422784335
c
c
        DOUBLE PRECISION X,PS

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MPSI'
      write ( *, '(a)' ) '  Test PSI'
      write ( *, '(a)' ) ' '

        WRITE(*,*)'    x          Psi(x)'
        WRITE(*,*)' ------------------------'

      do i = 1, 8
        x = dble ( i ) / 4.0D+00
        CALL PSI(X,PS)
        WRITE(*,10)X,PS
      end do

10      FORMAT(1X,F6.2,F18.9)
      return
        END
      subroutine MRCTJ ( )

c*********************************************************************72
c
cc MRCTJ tests RCTJ.
c
c  Modified:
c
c    06 April 2012
c
c
c       Purpose: This program computes the Riccati-Bessel 
c                functions of the first kind, and their
c                derivatives using subroutine RCTJ
c       Input:   x --- Argument of Riccati-Bessel function
c                n --- Order of jn(x)  ( 0 ó n ó 250 )
c       Output:  RJ(n) --- xújn(x)
c                DJ(n) --- [xújn(x)]'
c       Example: x = 10.0
c                  n        xújn(x)             [xújn(x)]'
c                --------------------------------------------
c                  0    -.5440211109D+00    -.8390715291D+00
c                  1     .7846694180D+00    -.6224880527D+00
c                  2     .7794219363D+00     .6287850307D+00
c                  3    -.3949584498D+00     .8979094712D+00
c                  4    -.1055892851D+01     .2739869063D-01
c                  5    -.5553451162D+00    -.7782202931D+00
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION RJ(0:250),DJ(0:250)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MRCTJ'
      write ( *, '(a)' ) '  Test RCTJ'
      write ( *, '(a)' ) ' '

      x = 10.0D+00
      n = 5

        IF (N.LE.10) THEN
           NS=1
        ELSE
           ns = 5
        ENDIF

        CALL RCTJ(N,X,NM,RJ,DJ)
        WRITE(*,*)'  n        xújn(x)             [xújn(x)]'''
        WRITE(*,*)'--------------------------------------------'
        DO 10 K=0,NM,NS
10         WRITE(*,20)K,RJ(K),DJ(K)
20      FORMAT(1X,I3,2D20.10)
30      FORMAT(3X,6HNmax =,I3,',    ',3Hx =,F7.2)
      return
        END
      subroutine MRCTY ( )

c*********************************************************************72
c
cc MRCTY tests RCTY.
c
c  Modified:
c
c    10 April 2012
c
c
c       Purpose: This program computes the Riccati-Bessel 
c                functions of the second kind and their
c                derivatives using subroutine RCTY
c       Input:   x --- Argument of Riccati-Bessel function
c                n --- Order of yn(x)
c       Output:  RY(n) --- xúyn(x)
c                DY(n) --- [xúyn(x)]'
c       Example: x = 10.0
c                  n        xúyn(x)             [xúyn(x)]'
c                --------------------------------------------
c                  0     .8390715291D+00    -.5440211109D+00
c                  1     .6279282638D+00     .7762787027D+00
c                  2    -.6506930499D+00     .7580668738D+00
c                  3    -.9532747888D+00    -.3647106133D+00
c                  4    -.1659930220D-01    -.9466350679D+00
c                  5     .9383354168D+00    -.4857670106D+00
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION RY(0:250),DY(0:250)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MRCTY'
      write ( *, '(a)' ) '  Test RCTY'
      write ( *, '(a)' ) ' '

      n = 5
      x = 10.0D+00

        IF (N.LE.10) THEN
           NS=1
        ELSE
           NS = 5
        ENDIF
        WRITE(*,*)
        CALL RCTY(N,X,NM,RY,DY)
        WRITE(*,*)
        WRITE(*,*)'  n        xúyn(x)             [xúyn(x)]'''
        WRITE(*,*)'--------------------------------------------'
        DO 10 K=0,NM,NS
           WRITE(*,20)K,RY(K),DY(K)
10      CONTINUE
20      FORMAT(1X,I3,2D20.10)
30      FORMAT(3X,6HNmax =,I3,',    ',3Hx =,F6.2)
      return
      END
      subroutine MRSWFO ( )

c*********************************************************************72
c
cc MRSWFO tests RSWFO.
c
c  Modified:
c
c    10 April 2012
c
c
c       Purpose: This program computes the radial oblate spheriodal 
c                functions of the first and second kinds, and their 
c                derivatives using subroutine RSWFO
c       Input :  m  --- Mode parameter, m = 0,1,2,...
c                n  --- Mode parameter, n = m,m+1,m+2,...
c                c  --- Spheroidal parameter
c                x  --- Argument (x ò 0)
c                cv --- Characteristic value
c                KF --- Function code
c                       KF=1 for the first kind
c                       KF=2 for the second kind
c                       KF=3 for both the first and second kinds
c       Output:  R1F --- Radial function of the first kind
c                R1D --- Derivative of the radial function of
c                        the first kind
c                R2F --- Radial function of the second kind
c                R2D --- Derivative of the radial function of
c                        the second kind
c       Example:
c                KD=-1, m = 2, n = 3, c = 5.0 and cv = 2.10980581604
c
c    x  R23(1)(-ic,ix) R23(1)'(-ic,ix) R23(2)(-ic,ix) R23(2)'(-ic,ix)
c  ------------------------------------------------------------------
c   0.0      0.0000000 (-1) 4.9911346  (-1)-4.0071049  (-1) 3.3065682
c   0.5 (-1) 2.0215966 (-1) 1.9559661  (-1)-1.5154835  (-1) 6.4482526
c   1.0 (-1) 1.0526737 (-1)-5.4010624  (-1) 1.3065731  (-1) 2.7958491
c   1.5 (-1)-1.1611246 (-2)-8.9414690  (-2) 3.7405936  (-1)-5.0118500
c   5.0 (-2) 3.6463251 (-2)-8.3283065  (-2) 1.5549615  (-1) 1.7544481
c
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EG(200)

      double precision x_test(5)

      save x_test

      data x_test / 
     &  0.0D+00, 0.5D+00, 1.0D+00, 1.5D+00, 5.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MRSWFO'
      write ( *, '(a)' ) '  Test RSWFO'
      write ( *, '(a)' ) ' '
      WRITE(*,*)'   x    Rmn(1)(-ic,ix)  Rmn''(1)(-ic,ix)  ',
     &            ' Rmn(2)(-ic,ix)  Rmn''(2)(-ic,ix)'
      WRITE(*,*)'  ---------------------------------------------',
     &            '---------------------------'

      kf = 3
      m = 2
      n = 3
      c = 5.0D+00

      do i = 1, 5
        x = x_test(i)

      CALL SEGV(M,N,C,-1,CV,EG)
      WRITE(*,20)M,N,C,CV,X
      WRITE(*,*)
      CALL RSWFO(M,N,C,X,CV,KF,R1F,R1D,R2F,R2D)


      WRITE(*,30)X,R1F,R1D,R2F,R2D
       WRITE(*,50)R1F*R2D-R2F*R1D,1.0D0/(C*(X*X+1.0D0))

      end do

10      FORMAT(1X,3HKF=,I3)
20      FORMAT(1X,2Hm=,I2,',   ',2Hn=,I2,',   ',2Hc=,F5.1,
     &         ',    ',4Hcv =,F18.10,',  ',2Hx=,F5.2)
30      FORMAT(1X,F5.1,4D17.8)
40      FORMAT(1X,F5.1,34X,4D17.8)
50      FORMAT(1X,/1X,'Wronskian check:',/1X,'Computed value =',
     &         D17.8,5X,'Exact value =',D17.8)
60      FORMAT(1X,/1X,'Caution: This check is not accurate if it ',
     &        'involves',/1X,'         the subtraction of two ',
     &        'similar numbers')
      END
      subroutine MRSWFP ( )

c*********************************************************************72
c
cc MRSWFP tests RSWFP.
c
c  Modified:
c
c    10 April 2012
c
c       Purpose: This program computes the radial prolate spheriodal 
c                functions of the first and second kinds, and their 
c                derivatives using subroutine RSWFP
c       Input :  m  --- Mode parameter, m = 0,1,2,...
c                n  --- Mode parameter, n = m,m+1,m+2,...
c                c  --- Spheroidal parameter
c                cv --- Characteristic value
c                x  --- Argument of radial function ( x > 1.0 )
c                KF --- Function code
c                       KF=1 for the first kind
c                       KF=2 for the second kind
c                       KF=3 for both the first and second kinds
c       Output:  R1F --- Radial function of the first kind
c                R1D --- Derivative of the radial function of
c                        the first kind
c                R2F --- Radial function of the second kind
c                R2D --- Derivative of the radial function of
c                        the second kind
c       Example:
c                KD= 1, m = 2, n = 3, c = 5.0 and cv =19.1359819110
c
c    x      R23(1)(c,x)   R23(1)'(c,x)    R23(2)(c,x)   R23(2)'(c,x)
c -------------------------------------------------------------------
c 1.(7).1 (-8) 1.6240735      1.6240735 ( 6)-3.0786785 (14) 3.0786784
c  1.005  (-3) 8.0600009      1.5998506     -6.2737713 ( 3) 1.2299041
c   1.1   (-1) 1.3578875      1.0702727 (-1)-4.3693218      3.5698419
c   1.5   (-1) 1.2002958 (-1)-8.0657929 (-1) 1.2623910 (-1) 4.8469847
c   5.0   (-2) 3.9455888 (-2) 4.2556949 (-2)-1.0099734 (-1) 2.0031280
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EG(200)

      double precision x_test(5)

      save x_test

      data x_test / 
     &  1.00000001D+00, 1.005D+00, 1.1D+00, 1.5D+00, 5.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MRSWFP'
      write ( *, '(a)' ) '  Test RSWFP'
      write ( *, '(a)' ) ' '

      kf = 3
      m = 2
      n = 3
      c = 5.0D+00

      WRITE(*,*)'   x      Rmn(1)(c,x)      Rmn''(1)(c,x)     ',
     &            'Rmn(2)(c,x)      Rmn''(2)(c,x)'
      WRITE(*,*)'-----------------------------------------------',
     &            '---------------------------'

      do i = 1, 5
        x = x_test(i)
        CALL SEGV(M,N,C,1,CV,EG)
        WRITE(*,20)M,N,C,CV,X
        WRITE(*,*)
        CALL RSWFP(M,N,C,X,CV,KF,R1F,R1D,R2F,R2D)

        WRITE(*,30)X,R1F,R1D,R2F,R2D
        WRITE(*,50)R1F*R2D-R2F*R1D,1.0D0/(C*(X*X-1.0D0))
        WRITE(*,60)
      end do

10      FORMAT(1X,3HKF=,I3)
20      FORMAT(1X,2Hm=,I2,',   ',2Hn=,I2,',   ',2Hc=,F5.1,
     &         ',   ',4Hcv =,D20.12,',  ',2Hx=,F5.2)
30      FORMAT(1X,F5.2,4D17.8)
40      FORMAT(1X,F5.2,34X,4D17.8)
50      FORMAT(1X,/1X,'Wronskian check:',/1X,'Computed value =',
     &         D17.8,5X,'Exact value =',D17.8)
60      FORMAT(1X,/1X,'Caution: This check is not accurate if it ',
     &        'involves',/1X,'         the subtraction of two ',
     &        'similar numbers')
      return
      END
      subroutine MSCKA ( )

c*********************************************************************72
c
cc MSCKA tests SCKA.
c
c  Modified:
c
c    10 April 2012
c
c       Purpose: This program computes the expansion coefficients 
c                of the prolate and oblate spheroidal functions, 
c                c2k, using subroutine SCKA
c       Input :  m  --- Mode parameter
c                n  --- Mode parameter
c                c  --- Spheroidal parameter
c                cv --- Characteristic value
c                KD --- Function code
c                       KD=1 for prolate; KD=-1 for oblate
c       Output:  CK(k) --- Expansion coefficients ck;
c                          CK(1), CK(2),... correspond to
c                          c0, c2,...
c       Example: Compute the first 13 expansion coefficients C2k for
c                KD= 1, m=2, n=3, c=3.0 and cv=14.8277782138; and
c                KD=-1, m=2, n=3, c=3.0 and cv=8.80939392077
c
c                Coefficients of Prolate and oblate functions
c
c                  k         C2k(c)              C2k(-ic)
c                ---------------------------------------------
c                  0     .9173213327D+01     .2489664942D+02
c                  1     .4718258929D+01    -.1205287032D+02
c                  2     .9841212916D+00     .2410564082D+01
c                  3     .1151870224D+00    -.2735821590D+00
c                  4     .8733916403D-02     .2026057157D-01
c                  5     .4663888254D-03    -.1061946315D-02
c                  6     .1853910398D-04     .4158091152D-04
c                  7     .5708084895D-06    -.1264400411D-05
c                  8     .1402786472D-07     .3074963448D-07
c                  9     .2817194508D-09    -.6120579463D-09
c                 10     .4712094447D-11     .1015900041D-10
c                 11     .6667838485D-13    -.1427953361D-12
c                 12     .8087995432D-15     .1721924955D-14
c
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CK(200),EG(200)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MSCKA'
      write ( *, '(a)' ) '  Test SCKA'
      write ( *, '(a)' ) ' '

      kd = 1
      m = 2
      n = 2
      c = 3.0D+00

      CALL SEGV(M,N,C,KD,CV,EG)
      WRITE(*,30)KD,M,N,C,CV
      CALL SCKA(M,N,C,CV,KD,CK)
       WRITE(*,*)'Coefficients of Prolate function'
       WRITE(*,*)
       WRITE(*,*)'   k            C2k(c)'
      WRITE(*,*)'----------------------------'
      NM=25+INT((N-M)/2+C)
      DO K=1,NM
           WRITE(*,20)K-1,CK(K)
      end do

      kd = -1
      m = 2
      n = 2
      c = 3.0D+00

      CALL SEGV(M,N,C,KD,CV,EG)
      WRITE(*,30)KD,M,N,C,CV
      CALL SCKA(M,N,C,CV,KD,CK)
      WRITE(*,*)
       WRITE(*,*)'Coefficients of Oblate function'
       WRITE(*,*)
       WRITE(*,*)'   k           C2k(-ic)'
      WRITE(*,*)'----------------------------'
      NM=25+INT((N-M)/2+C)
      DO K=1,NM
           WRITE(*,20)K-1,CK(K)
      end do

20      FORMAT(2X,I3,4X,D18.10)
30      FORMAT(1X,3HKD=,I3,',  ',2Hm=,I3,',  ',2Hn=,I3,',  ',2Hc=,
     &         F5.1,',  ',4Hcv =,F18.10)
      return
      END
      subroutine MSCKB ( )

c*********************************************************************72
c
cc MSCKB tests SCKB.
c
c  Modified:
c
c    09 April 2012
c
c       Purpose: This program computes the expansion coefficients 
c                of the prolate and oblate spheroidal functions, 
c                c2k, using subroutine SCKB
c       Input :  m  --- Mode parameter
c                n  --- Mode parameter
c                c  --- Spheroidal parameter
c                cv --- Characteristic value
c                KD --- Function code
c                       KD=1 for prolate; KD=-1 for oblate
c       Output:  CK(k) --- Expansion coefficients ck;
c                          CK(1), CK(2), ... correspond to
c                          c0, c2, ...
c       Example: Compute the first 13 expansion coefficients C2k for
c                KD= 1, m=2, n=3, c=3.0 and cv=14.8277782138; and
c                KD=-1, m=2, n=3, c=3.0 and cv=8.80939392077
c
c                Coefficients of Prolate and oblate functions
c
c                  k         C2k(c)              C2k(-ic)
c                ---------------------------------------------
c                  0     .9173213327D+01     .2489664942D+02
c                  1     .4718258929D+01    -.1205287032D+02
c                  2     .9841212916D+00     .2410564082D+01
c                  3     .1151870224D+00    -.2735821590D+00
c                  4     .8733916403D-02     .2026057157D-01
c                  5     .4663888254D-03    -.1061946315D-02
c                  6     .1853910398D-04     .4158091152D-04
c                  7     .5708084895D-06    -.1264400411D-05
c                  8     .1402786472D-07     .3074963448D-07
c                  9     .2817194508D-09    -.6120579463D-09
c                 10     .4712094447D-11     .1015900041D-10
c                 11     .6667838485D-13    -.1427953361D-12
c                 12     .8087995432D-15     .1721924955D-14
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CK(200),DF(200),EG(200)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MSCKB'
      write ( *, '(a)' ) '  Test SCKB'
      write ( *, '(a)' ) ' '

      kd = 1
      m = 2
      n = 2
      c = 3.0D+00
      CALL SEGV(M,N,C,KD,CV,EG)
      WRITE(*,30)KD,M,N,C,CV
      CALL SDMN(M,N,C,CV,KD,DF)
      CALL SCKB(M,N,C,DF,CK)
      WRITE(*,*)
      WRITE(*,*)'Coefficients of Prolate function'
      WRITE(*,*)
      WRITE(*,*)'   k            C2k(c)'
      WRITE(*,*)'----------------------------'
      NM=25+INT((N-M)/2+C)
      DO K=1,NM
        WRITE(*,20)K-1,CK(K)
      end do

      kd = -1
      m = 2
      n = 2
      c = 3.0D+00
      CALL SEGV(M,N,C,KD,CV,EG)
      WRITE(*,30)KD,M,N,C,CV
      CALL SDMN(M,N,C,CV,KD,DF)
      CALL SCKB(M,N,C,DF,CK)
      WRITE(*,*)
       WRITE(*,*)'Coefficients of Oblate function'
       WRITE(*,*)
       WRITE(*,*)'   k           C2k(-ic)'
      WRITE(*,*)'----------------------------'
      NM=25+INT((N-M)/2+C)
      DO K=1,NM
        WRITE(*,20)K-1,CK(K)
      end do

20      FORMAT(2X,I3,4X,D18.10)
30      FORMAT(1X,3HKD=,I3,',  ',2Hm=,I3,',  ',2Hn=,I3,
     &         ',  ',2Hc=,F5.1,',  ',4Hcv =,F18.10)
      return
      END
      subroutine MSDMN ( )

c*********************************************************************72
c
cc MSDMN tests SDMN.
c
c  Modified:
c
c    09 April 2012
c
c       Purpose: This program computes the expansion coefficients 
c                of the prolate and oblate spheroidal functions, 
c                dk, using subroutine SDMN
c       Input :  m  --- Mode parameter
c                n  --- Mode parameter
c                c  --- Spheroidal parameter
c                cv --- Characteristic value
c                KD --- Function code
c                       KD=1 for prolate; KD=-1 for oblate
c       Output:  DF(k) --- Expansion coefficients dk;
c                          DF(1), DF(2),... correspond to
c                          d0, d2,... for even n-m and d1,
c                          d3,... for odd n-m
c       Example: Compute the first 12 expansion coefficients for
c                KD= 1, m=2, n=2, c=3.0 and cv=7.1511005241; and
c                KD=-1, m=2, n=2, c=3.0 and cv=4.5264604622
c
c                Coefficients of Prolate and oblate functions
c
c                  r          dr(c)             dr(-ic)
c                -------------------------------------------
c                  0     .9237882817D+00    .1115434000D+01
c                  2    -.2901607696D-01    .4888489020D-01
c                  4     .8142246173D-03    .1600845667D-02
c                  6    -.1632270292D-04    .3509183384D-04
c                  8     .2376699010D-06    .5416293446D-06
c                 10    -.2601391701D-08    .6176624069D-08
c                 12     .2209142844D-10    .5407431236D-10
c                 14    -.1494812074D-12    .3745889118D-12
c                 16     .8239302207D-15    .2103624480D-14
c                 18    -.3768260778D-17    .9768323113D-17
c                 20     .1452384658D-19    .3812753620D-19
c                 22    -.4780280430D-22    .1268321726D-21
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION DF(200),EG(200)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MSDMN'
      write ( *, '(a)' ) '  Test SDMN'
      write ( *, '(a)' ) ' '

      kd = 1
      m = 2
      n = 2
      c = 3.0D+00
      CALL SEGV(M,N,C,KD,CV,EG)
      WRITE(*,30)KD,M,N,C,CV
      CALL SDMN(M,N,C,CV,KD,DF)
       WRITE(*,*)'Coefficients of Prolate function'
       WRITE(*,*)
       WRITE(*,*)'   r             dr(c)'
      WRITE(*,*)'----------------------------'
      NM=25+INT(0.5*(N-M)+C)         
      DO K=1,NM 
       IF (N-M.EQ.2*INT((N-M)/2)) THEN
          J=2*(K-1)
       ELSE
          J=2*K-1
        ENDIF
        WRITE(*,20)J,DF(K)
      end do

      kd = -1
      m = 2
      n = 2
      c = 3.0D+00
      CALL SEGV(M,N,C,KD,CV,EG)
      WRITE(*,30)KD,M,N,C,CV
      CALL SDMN(M,N,C,CV,KD,DF)
       WRITE(*,*)'Coefficients of Oblate function'
       WRITE(*,*)
       WRITE(*,*)'   r             dr(c)'
      WRITE(*,*)'----------------------------'
      NM=25+INT(0.5*(N-M)+C)         
      DO K=1,NM 
       IF (N-M.EQ.2*INT((N-M)/2)) THEN
          J=2*(K-1)
       ELSE
          J=2*K-1
        ENDIF
        WRITE(*,20)J,DF(K)
      end do

20      FORMAT(2X,I3,4X,D18.10)
30      FORMAT(1X,3HKD=,I3,',  ',2Hm=,I3,',  ',2Hn=,I3,',  ',2Hc=,
     &         F5.1,',  ',4Hcv =,F18.10)

      return
      END
      subroutine MSEGV ( )

c*********************************************************************72
c
cc MSEGV tests SEGV.
c
c  Modified:
c
c    06 April 2012
c
c       Purpose: This program computes a sequence of characteristic 
c                values of spheroidal prolate and oblate functions 
c                using subroutine SEGV
c       Input :  m  --- Mode parameter
c                n  --- Mode parameter
c                c  --- Spheroidal parameter
c                KD --- Function code
c                       KD=1 for Prolate; KD=-1 for Oblate
c       Output:  CV --- Characteristic value for given m, n and c
c                EG(L) --- Characteristic value for mode m and n'
c                          ( L = n' - m + 1 )
c       Examples:
c                Prolate: ( KD = 1 , m = 1, n = 5, c = 5.0 )
c
c                m      n       c        Lambda mn(c)
c              ---------------------------------------
c                1      1      5.0        5.35042230
c                1      2      5.0       14.64295624
c                1      3      5.0       23.39761312
c                1      4      5.0       32.42194359
c                1      5      5.0       42.65818215
c
c                Oblate: ( KD = -1 , m = 1, n = 5, c = 5.0 )
c
c                m      n       c      Lambda mn(-ic)
c               --------------------------------------
c                1      1      5.0       -7.49338828
c                1      2      5.0       -7.12783752
c                1      3      5.0        2.75036721
c                1      4      5.0        8.69495925
c                1      5      5.0       18.43931577
c
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EG(100)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MSEGV'
      write ( *, '(a)' ) '  Test SEGV'
      write ( *, '(a)' ) ' '

      kd = 1
      m = 1
      n = 5
      c = 5.0D+00
      CALL SEGV(M,N,C,KD,CV,EG)
      WRITE(*,*)'  m      n       c       Lambda mn(c)'
        WRITE(*,*)'---------------------------------------'
        DO L=1,N-M+1
           N1=M+L-1
           WRITE(*,20)M,N1,C,EG(L)
        end do

      kd = -1
      m = 1
      n = 5
      c = 5.0D+00
      CALL SEGV(M,N,C,KD,CV,EG)
      write ( *, '(a)' ) ' '
      WRITE(*,*)'  m      n       c      Lambda mn(-ic)'
        WRITE(*,*)'---------------------------------------'
        DO L=1,N-M+1
           N1=M+L-1
           WRITE(*,20)M,N1,C,EG(L)
      end do

15      FORMAT(1X,'KD =',I2,',  m =',I3,',   n =',I3,',  c =',F5.1)
20      FORMAT(1X,I3,4X,I3,4X,F5.1,F18.8)
      return
      END
      subroutine MSPHI ( )

c*********************************************************************72
c
cc MSPHI tests SPHI.
c
c  Modified:
c
c    08 April 2012
c
c
c       Purpose: This program computes the modified spherical 
c                Bessel functions of the first kind in(x) and 
c                in'(x) using subroutine SPHI
c       Input :  x --- Argument of in(x)
c                n --- Order of in(x) ( 0 ó n ó 250 )
c       Output:  SI(n) --- in(x)
c                DI(n) --- in'(x)
c       Example: x = 10.0
c                  n          in(x)               in'(x)
c                --------------------------------------------
c                  0     .1101323287D+04     .9911909633D+03
c                  1     .9911909633D+03     .9030850948D+03
c                  2     .8039659985D+03     .7500011637D+03
c                  3     .5892079640D+03     .5682828129D+03
c                  4     .3915204237D+03     .3934477522D+03
c                  5     .2368395827D+03     .2494166741D+03
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION SI(0:250),DI(0:250)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MSPHI'
      write ( *, '(a)' ) '  Test SPHI'
      write ( *, '(a)' ) ' '

      n = 5
      x = 10.0D+00

        IF (N.LE.10) THEN
           NS=1
        ELSE
           NS = 5
        ENDIF
        CALL SPHI(N,X,NM,SI,DI)
        WRITE(*,*)'  n          in(x)               in''(x)'
        WRITE(*,*)'--------------------------------------------'
        DO 10 K=0,NM,NS
10         WRITE(*,20)K,SI(K),DI(K)
20      FORMAT(1X,I3,2D20.10)
30      FORMAT(3X,'Nmax =',I3,',     ','x =',F6.1)
      return
      END
      subroutine MSPHJ ( )

c*********************************************************************72
c
cc MSPHJ tests SPHJ.
c
c  Modified:
c
c    06 April 2012
c
c
c       Purpose: This program computes the spherical Bessel 
c                functions jn(x) and jn'(x) using subroutine
c                SPHJ
c       Input :  x --- Argument of jn(x)
c                n --- Order of jn(x)  ( n = 0,1,úúú,ó 250 )
c       Output:  SJ(n) --- jn(x)
c                DJ(n) --- jn'(x)
c       Example:   x =10.0
c                  n          jn(x)               jn(x)
c                --------------------------------------------
c                  0    -.5440211109D-01    -.7846694180D-01
c                  1     .7846694180D-01    -.7009549945D-01
c                  2     .7794219363D-01     .5508428371D-01
c                  3    -.3949584498D-01     .9374053162D-01
c                  4    -.1055892851D+00     .1329879757D-01
c                  5    -.5553451162D-01    -.7226857814D-01
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION SJ(0:250),DJ(0:250)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MSPHJ'
      write ( *, '(a)' ) '  Test SPHJ'
      write ( *, '(a)' ) ' '

      n = 5
      x = 10.0D+00

        IF (N.LE.10) THEN
           NS=1
        ELSE
           ns = 5
        ENDIF
        CALL SPHJ(N,X,NM,SJ,DJ)

        WRITE(*,*)'  n          jn(x)               jn''(x)'
        WRITE(*,*)'--------------------------------------------'
        DO 10 K=0,NM,NS
10         WRITE(*,20)K,SJ(K),DJ(K)
20      FORMAT(1X,I3,2D20.10)
30      FORMAT(3X,6HNmax =,I3,',     ',2Hx=,F5.1)
      return
      END
      subroutine MSPHK ( )

c*********************************************************************72
c
cc MSPHK tests SPHK.
c
c  Modified:
c
c    06 April 2012
c
c
c       Purpose: This program computes the modified spherical 
c                Bessel functions kn(x) and kn'(x) using
c              subroutine SPHK
c       Input :  x --- Argument of kn(x)  ( x ò 0 )
c                n --- Order of kn(x) ( n ó 250 )
c       Output:  SK(n) --- kn(x)
c                DK(n) --- kn'(x)
c       Example: x= 10.0
c                  n          kn(x)               kn'(x)
c                --------------------------------------------
c                  0     .7131404291D-05    -.7844544720D-05
c                  1     .7844544720D-05    -.8700313235D-05
c                  2     .9484767707D-05    -.1068997503D-04
c                  3     .1258692857D-04    -.1451953914D-04
c                  4     .1829561771D-04    -.2173473743D-04
c                  5     .2905298451D-04    -.3572740841D-04
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION SK(0:250),DK(0:250)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MSPHK'
      write ( *, '(a)' ) '  Test SPHK'
      write ( *, '(a)' ) ' '

      n = 5
      x = 10.0D+00

        IF (N.LE.10) THEN
           NS=1
        ELSE
           NS = 5
        ENDIF
        CALL SPHK(N,X,NM,SK,DK)
        WRITE(*,*)
        WRITE(*,*)'  n          kn(x)               kn''(x)'
        WRITE(*,*)'--------------------------------------------'
        DO 10 K=0,NM,NS
10         WRITE(*,20)K,SK(K),DK(K)
20      FORMAT(1X,I3,2D20.10)
30      FORMAT(3X,'Nmax =',I3,',     ','x =',F6.1)
      return
      END
      subroutine MSPHY ( )

c*********************************************************************72
c
cc MSPHY tests SPHY.
c
c  Modified:
c
c    08 April 2012
c
c
c       Purpose: This program computes the spherical Bessel 
c                functions yn(x) and yn'(x) using subroutine
c                SPHY
c       Input :  x --- Argument of yn(x) ( x ò 0 )
c                n --- Order of yn(x) ( n = 0,1,úúú, ó 250 )
c       Output:  SY(n) --- yn(x)
c                DY(n) --- yn'(x)
c       Example:   x = 10.0
c                  n          yn(x)               yn'(x)
c                --------------------------------------------
c                  0     .8390715291D-01    -.6279282638D-01
c                  1     .6279282638D-01     .7134858763D-01
c                  2    -.6506930499D-01     .8231361788D-01
c                  3    -.9532747888D-01    -.2693831344D-01
c                  4    -.1659930220D-02    -.9449751377D-01
c                  5     .9383354168D-01    -.5796005523D-01
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION SY(0:250),DY(0:250)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MSPHY'
      write ( *, '(a)' ) '  Test SPHY'
      write ( *, '(a)' ) ' '

      x = 10.0D+00
      n = 5

        IF (N.LE.10) THEN
           NS=1
        ELSE
           NS = 5
        ENDIF
        CALL SPHY(N,X,NM,SY,DY)
        WRITE(*,*)'  n          yn(x)               yn''(x)'
        WRITE(*,*)'--------------------------------------------'
        DO 10 K=0,NM,NS
10         WRITE(*,20)K,SY(K),DY(K)
20    FORMAT(1X,I3,2D20.10)
30    FORMAT(3X,6HNmax =,I3,',     ',2Hx=,F6.1)
      return
      END
      subroutine MSTVH0 ( )

c*********************************************************************72
c
cc MSTVH0 tests STVH0.
c
c  Modified:
c
c    09 April 2012
c
c
c       Purpose: This program computes Struve function 
c                H0(x) using subroutine STVH0
c       Input :  x   --- Argument of H0(x) ( x ò 0 )
c       Output:  SH0 --- H0(x)
c       Example:
c                   x          H0(x)
c                ----------------------
c                  0.0       .00000000
c                  5.0      -.18521682
c                 10.0       .11874368
c                 15.0       .24772383
c                 20.0       .09439370
c                 25.0      -.10182519
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MSTVH0'
      write ( *, '(a)' ) '  Test STVH0'
      write ( *, '(a)' ) ' '
        WRITE(*,*)'   x          H0(x)'
        WRITE(*,*)'----------------------'
      do i = 0, 5
        x = dble ( 5 * i )
        CALL STVH0(X,SH0)
        WRITE(*,10)X,SH0
      end do
10    FORMAT(1X,F5.1,E16.8)
      return
      END
      subroutine MSTVH1 ( )

c*********************************************************************72
c
cc MSTVH1 tests STVH1.
c
c  Modified:
c
c    09 April 2012
c
c
c       Purpose: This program computes Struve function 
c                H1(x) using subroutine STVH1
c       Input :  x   --- Argument of H1(x) ( x ò 0 )
c       Output:  SH1 --- H1(x)
c       Example:
c                   x          H1(x)
c                -----------------------
c                  0.0       .00000000
c                  5.0       .80781195
c                 10.0       .89183249
c                 15.0       .66048730
c                 20.0       .47268818
c                 25.0       .53880362
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MSTVH1'
      write ( *, '(a)' ) '  Test STVH1'
      write ( *, '(a)' ) ' '
        WRITE(*,*)'   x          H1(x)'
        WRITE(*,*)'-----------------------'
      do i = 0, 5
        x = dble ( 5 * i )
        CALL STVH1(X,SH1)
        WRITE(*,10)X,SH1
      end do
10      FORMAT(1X,F5.1,E16.8)
      return
      END
      subroutine MSTVHV ( )

c*********************************************************************72
c
cc MSTVHV tests STVHV.
c
c  Modified:
c
c    08 April 2012
c
c
c       Purpose:  This program computes Struve function Hv(x) 
c                 for an arbitrary order using subroutine
c                 STVHV
c       Input :   v  --- Order of Hv(x)  ( -8.0 ó v ó 12.5 )
c                 x  --- Argument of Hv(x) ( x ò 0 )
c       Output:   HV --- Hv(x)
c       Example:  x = 10.0
c                   v           Hv(x)
c                 -----------------------
c                   .5     .46402212D+00
c                  1.5     .14452322D+01
c                  2.5     .31234632D+01
c                  3.5     .53730255D+01
c                  4.5     .72083122D+01
c                  5.5     .76851132D+01
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MSTVHV'
      write ( *, '(a)' ) '  Test STVHV'
      write ( *, '(a)' ) ' '

      x = 10.0D+00
        WRITE(*,*)'   v           Hv(x)'
        WRITE(*,*)' -----------------------'
      do i = 1, 11, 2
        x = dble ( i ) / 2.0D+00
        CAll STVHV(V,X,HV)
        WRITE(*,20)V,HV
      end do
20    FORMAT(1X,F5.1,D18.8)
30    FORMAT(1X,'v =',F5.1,6X,'x =',F5.1)
      return
      END
      subroutine MSTVL0 ( )

c*********************************************************************72
c
cc MSTVL0 tests STVL0.
c
c  Modified:
c
c    06 April 2012
c
c
c       Purpose: This program computes modified Struve 
c                function L0(x) using subroutine STVL0
c       Input :  x   --- Argument of L0(x) ( x ò 0 )
c       Output:  SL0 --- L0(x)
c       Example:
c                   x        L0(x)
c               ------------------------
c                  0.0   .00000000D+00
c                  5.0   .27105917D+02
c                 10.0   .28156522D+04
c                 15.0   .33964933D+06
c                 20.0   .43558283D+08
c                 30.0   .78167230D+12
c                 40.0   .14894775D+17
c                 50.0   .29325538D+21
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MSTVL0'
      write ( *, '(a)' ) '  Test STVL0'
      write ( *, '(a)' ) ' '
      WRITE(*,*)'   x        L0(x)'
      WRITE(*,*)'-----------------------'
      do i = 0, 10
        x = dble ( 5 * i )
        CALL STVL0(X,SL0)
        WRITE(*,10)X,SL0
      end do
10    FORMAT(1X,F5.1,D16.8)
      return
      END
      subroutine MSTVL1 ( )

c*********************************************************************72
c
cc MSTVL1 tests STVL1.
c
c  Modified:
c
c    06 April 2012
c
c
c       Purpose: This program computes the modified Struve 
c                function L1(x) using subroutine STVL1
c       Input :  x   --- Argument of L1(x) ( x ò 0 )
c       Output:  SL1 --- L1(x)
c       Example:
c                     x        L1(x)
c                 -----------------------
c                   0.0   .00000000D+00
c                   5.0   .23728216D+02
c                  10.0   .26703583D+04
c                  15.0   .32812429D+06
c                  20.0   .42454973D+08
c                  30.0   .76853204D+12
c                  40.0   .14707396D+17
c                  50.0   .29030786D+21
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MSTVL1'
      write ( *, '(a)' ) '  Test STVL1'
      write ( *, '(a)' ) ' '
      WRITE(*,*)'   x        L1(x)'
      WRITE(*,*)'-----------------------'
      do i = 0, 10
        x = dble ( 5 * i )
        CALL STVL1(X,SL1)
        WRITE(*,10)X,SL1
      end do
10    FORMAT(1X,F5.1,D16.8)
      return
      END
      subroutine MSTVLV ( )

c*********************************************************************72
c
cc MSTVLV tests STVLV.
c
c  Modified:
c
c    06 April 2012
c
c
c       Purpose:  This program computes the modified Struve 
c                 function Lv(x) for an arbitrary order v
c                 using subroutine STVLV
c       Input :   v   --- Order of Lv(x)  ( |v| ó 20 )
c                 x   --- Argument of Lv(x) ( x ò 0 )
c       Output:   SLV --- Lv(x)
c       Example:  x = 10.0
c                   v          Lv(x)
c                 ------------------------
c                  0.5     .27785323D+04
c                  1.5     .24996698D+04
c                  2.5     .20254774D+04
c                  3.5     .14816746D+04
c                  4.5     .98173460D+03
c                  5.5     .59154277D+03
c
c
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MSTVLV'
      write ( *, '(a)' ) '  Test STVLV'
      write ( *, '(a)' ) ' '

      WRITE(*,*)'   v       x            Lv(x)'
      WRITE(*,*)' ---------------------------------'

      x = 10.0D+00

      do i = 1, 11, 2
        v = dble ( i ) / 2.0D+00
        CAll STVLV(V,X,SLV)
        WRITE(*,20)V,X,SLV
      end do

20    FORMAT(1X,F5.1,2x,f5.1,2x,D18.8)
 
      return
      END
