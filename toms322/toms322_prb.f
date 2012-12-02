      program main

c*********************************************************************72
c
cc MAIN is the main program for TOMS322_PRB.
c
c  Discussion:
c
c    Test function prob() that generates critical values of standard
c    probability functions
c
c    Larry Weissman 29-Mar-94  larryw@bioeng.washington.edu
c
c  Modified:
c
c    28 December 2007
c
      implicit none

      write ( *, '(a)' ) ' '
      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS322_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS322 library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS322_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests PROB for evaluating the two-tailed normal distribution.
c
c  Modified:
c
c    11 February 2008
c
      implicit none

      integer nz
      parameter ( nz = 13 )

      integer i
      double precision p
      double precision prob
      double precision pz(nz)
      double precision z(nz)

      data  z/4.417, 3.891, 3.291, 2.576, 2.241,
     1        1.960, 1.645, 1.150, 0.674, 0.319,
     2        0.126, 0.063, 0.0125/

      data pz/1.e-5, 1.e-4, 1.e-3, 1.e-2, 0.025,
     1        0.050, 0.100, 0.250, 0.500, 0.750,
     2        0.900, 0.950, 0.990/

10    format(/15x,a//4x,
     1'N1   N2       X      P-Exp       P-Calc      E-Diff     E-Rel'/
     22x,
     3'---- ----    ------  ----------  ----------  ----------   -----')

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write(*,10)'Testing two-tailed normal distribution'
      do i=1,nz
        p=1.0-prob(1,5000,z(i)**2)
        call report(1,5000,z(i),pz(i),p)
      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests PROB for evaluating the two-tailed Student's T distribution.
c
c  Modified:
c
c    11 February 2008
c
      implicit none

      integer nt
      parameter ( nt = 25 )

      integer i
      integer nut(nt)
      double precision p
      double precision prob
      double precision pt(nt)
      double precision t(nt)

      data   t/1.000, 0.817, 0.700, 0.687, 0.677,
     1         6.314, 2.920, 1.813, 1.725, 1.658,
     2        12.706, 4.303, 2.228, 2.086, 1.980,
     3        63.657, 9.925, 3.169, 2.845, 2.617,
     4        127.32,14.089, 3.581, 3.153, 2.860/
      data nut/    1,     2,    10,    20,   120,
     1             1,     2,    10,    20,   120,
     2             1,     2,    10,    20,   120,
     3             1,     2,    10,    20,   120,
     4             1,     2,    10,    20,   120/
      data  pt/0.500, 0.500, 0.500, 0.500, 0.500,
     1         0.100, 0.100, 0.100, 0.100, 0.100,
     2         0.050, 0.050, 0.050, 0.050, 0.050,
     3         0.010, 0.010, 0.010, 0.010, 0.010,
     4         0.005, 0.005, 0.005, 0.005, 0.005/

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'

      write(*,10)'Testing two-tailed Student''s T'
      do i=1,nt
        p=1.0-prob(1,nut(i),t(i)**2)
        call report(1,nut(i),t(i),pt(i),p)
      end do

10    format(/15x,a//4x,
     1'N1   N2       X      P-Exp       P-Calc      E-Diff     E-Rel'/
     22x,
     3'---- ----    ------  ----------  ----------  ----------   -----')

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests PROB for evaluating the Chi-square distribution.
c
c  Modified:
c
c    11 February 2008
c
      implicit none

      integer nx
      parameter ( nx = 25 )

      integer i
      integer nux(nx)
      double precision p
      double precision prob
      double precision px(nx)
      double precision x(nx)

      data   x/ 0.45,  1.39,  9.34, 19.34, 99.33,
     1          3.84,  5.99, 18.31, 31.41,124.34,
     2          5.02,  7.38, 20.48, 34.17,129.56,
     3          6.63,  9.21, 23.21, 37.57,135.81,
     4          7.88, 10.60, 25.19, 40.00,140.17/
      data nux/    1,     2,    10,    20,   100,
     1             1,     2,    10,    20,   100,
     2             1,     2,    10,    20,   100,
     3             1,     2,    10,    20,   100,
     4             1,     2,    10,    20,   100/
      data  px/0.500, 0.500, 0.500, 0.500, 0.500,
     1         0.050, 0.050, 0.050, 0.050, 0.050,
     2         0.025, 0.025, 0.025, 0.025, 0.025,
     3         0.010, 0.010, 0.010, 0.010, 0.010,
     4         0.005, 0.005, 0.005, 0.005, 0.005/

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'

      write(*,10)'Testing chi-square distribution'
      do i=1,nx
        p=1.0-prob(nux(i),5000,x(i)/nux(i))
        call report(nux(i),5000,x(i),px(i),p)
      end do

10    format(/15x,a//4x,
     1'N1   N2       X      P-Exp       P-Calc      E-Diff     E-Rel'/
     22x,
     3'---- ----    ------  ----------  ----------  ----------   -----')

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests PROB for evaluating the F distribution.
c
c  Modified:
c
c    11 February 2008
c
      implicit none

      integer nf
      parameter ( nf = 100 )

      double precision f(nf)
      integer i
      integer nf1(nf)
      integer nf2(nf)
      double precision p
      double precision pf(nf)
      double precision prob

      data   f/ 16211.,198.50,12.826,9.9439,8.1790,
     1          20000.,199.00,9.4270,6.9865,5.5393,
     2          24224.,199.40,5.8467,3.8470,2.7052,
     3          24836.,199.45,5.2740,3.3178,2.1881,
     4          25359.,199.49,4.7501,2.8058,1.6055,

     5          4052.2,98.503,10.044,8.0960,6.8510,
     6          4999.5,99.000,7.5594,5.8489,4.7865,
     7          6055.8,99.399,4.8492,3.3682,2.4721,
     8          6208.7,99.490,4.4054,2.9377,2.0346,
     9          6339.4,99.491,3.9965,2.5168,1.5330,

     a          161.45,18.513,4.9646,4.3513,3.9201,
     b          199.50,19.000,4.1028,3.4928,3.0718,
     c          241.88,19.396,2.9782,2.3479,1.9105,
     d          248.01,19.446,2.7740,2.1242,1.6587,
     e          253.25,19.487,2.5801,1.8963,1.3519,

     f          1.0000,.66667,.48973,.47192,.45774,
     g          1.5000,1.0000,.74349,.71773,.69717,
     h          2.0419,1.3450,1.0000,.96626,.93943,
     i          2.1190,1.3933,1.0349,1.0000,.97228,
     j          2.1848,1.4344,1.0645,1.0285,1.0000/
      data nf1/ 5* 1, 5  *2, 5* 10, 5* 20, 5*120,
     1          5* 1, 5  *2, 5* 10, 5* 20, 5*120,
     2          5* 1, 5  *2, 5* 10, 5* 20, 5*120,
     3          5* 1, 5  *2, 5* 10, 5* 20, 5*120/
      data nf2/1,2,10,20,120,1,2,10,20,120,1,2,10,20,120,1,2,10,20,120,
     1         1,2,10,20,120,1,2,10,20,120,1,2,10,20,120,1,2,10,20,120,
     2         1,2,10,20,120,1,2,10,20,120,1,2,10,20,120,1,2,10,20,120,
     3         1,2,10,20,120,1,2,10,20,120,1,2,10,20,120,1,2,10,20,120,
     4         1,2,10,20,120,1,2,10,20,120,1,2,10,20,120,1,2,10,20,120/
      data  pf/25*0.005, 25*0.01, 25*0.05, 25*0.50/

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write(*,10)'Testing F distribution'

      do i=1,nf
        p=1.0-prob(nf1(i),nf2(i),f(i))
        call report(nf1(i),nf2(i),f(i),pf(i),p)
      end do

10    format(/15x,a//4x,
     1'N1   N2       X      P-Exp       P-Calc      E-Diff     E-Rel'/
     22x,
     3'---- ----    ------  ----------  ----------  ----------   -----')

      return
      end
      subroutine report ( n1, n2, x, pexp, pcalc )

c*********************************************************************72
c
cc REPORT reports the results of a calculation.
c
c  Discussion:
c
c    This routine reports the difference between the expected probability
c    PEXP and the calculated probability PCALC.   It also reports the
c    relative error.
c
c  Modified:
c
c    11 February 2008
c
      implicit none

      integer n1
      integer n2
      double precision pcalc
      double precision pdiff
      double precision pexp
      double precision prel
      double precision x

      pdiff = abs ( pexp - pcalc )

      if ( abs ( pexp ) .gt. 1.0D-30 ) then
        prel = pdiff / pexp
      else
        prel = 0.0D+00
      end if

      write(*,10)n1,n2,x,pexp,pcalc,pdiff,prel
10    format(1x,i5,i5,f10.4,3g12.4,f8.3)

      return
      end

