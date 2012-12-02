      program main

c*********************************************************************72
c
cc MAIN is the main program of TOMS706_PRB2.
c
c  Discussion:
c
c    TOMS706_PRB2 tests some of the features of DCUTRI.
c
C    It checks that DCUTRI integrates to machine
C    precision all monomials of degree less or equal to 13.
c
C    It checks that DCUTRI integrates correctly a vector
C    of integrands.
c
C    It checks that the restart feature of DCUTRI works.
C
c  Modified:
c
c    09 December 2006
c
c  Author:
c
c    Jarle Berntsen, Terje Espelid
c
      implicit none

      external ftestp,ftestq,ftesto

      integer nw,iwork(501)
      parameter (nw=50000)
      double precision vertce(2,3,500),wrkstr(nw),abserr,
     &                 vertcq(2,3,500)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS706_PRB2:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  A test problem for DCUTRI.'

      vertce(1,1,1) = 0
      vertce(2,1,1) = 0
      vertce(1,2,1) = 1
      vertce(2,2,1) = 0
      vertce(1,3,1) = 0
      vertce(2,3,1) = 1
      vertce(1,1,2) = 1
      vertce(2,1,2) = 1
      vertce(1,2,2) = 1
      vertce(2,2,2) = 0
      vertce(1,3,2) = 0
      vertce(2,3,2) = 1
      vertcq(1,1,1) = -1
      vertcq(2,1,1) = 0
      vertcq(1,2,1) = 0.5d0
      vertcq(2,2,1) = -sqrt(3.d0)/2
      vertcq(1,3,1) = 0.5d0
      vertcq(2,3,1) = sqrt(3.d0)/2
      abserr = 0.6
C
C    TEST FOR INTEGRATING MONOMIALS
C     Selected monomials, degrees 0-13
C
C    We first check the 105 monomials in cartesian
C    coordinates which DCUTRI should integrate to machine
C    precision.
C
      call atest(vertce,1,500,1000,105,ftestp,abserr,nw,0,wrkstr,
     &           iwork)
C
C    Then the same for the 21 monomials in polar coordinates.
C
      call atest(vertcq,1,500,1000,21,ftestq,abserr,nw,0,wrkstr,
     &           iwork)
C
C    Test restart on a peak problem.
C
      vertce(1,1,1) = 0
      vertce(2,1,1) = 0
      vertce(1,2,1) = 1
      vertce(2,2,1) = 0
      vertce(1,3,1) = 0
      vertce(2,3,1) = 1
      vertce(1,1,2) = 1
      vertce(2,1,2) = 1
      vertce(1,2,2) = 1
      vertce(2,2,2) = 0
      vertce(1,3,2) = 0
      vertce(2,3,2) = 1
      abserr = 0.1

      call atest(vertce,2,500,74,1,ftesto,abserr,nw,0,wrkstr,
     &           iwork)
      abserr = 1d-3
      call atest(vertce,2,500,1000,1,ftesto,abserr,nw,1,wrkstr,
     &           iwork)
      call atest(vertce,2,500,2000,1,ftesto,abserr,nw,1,wrkstr,
     &           iwork)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DCUTRI_PRB2:'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
      subroutine atest(vertce,numtri,lenver,maxcls,nfun,tstsub,abserr,
     &                 lenwrk,irest,wrkstr,iwork)

c*********************************************************************72
c
cc ATEST
c
      implicit none

      external tstsub
      integer lenwrk,irest,neval,numtri,lenver,iwork(lenver+1)
      double precision vertce(2,3,lenver),absest(105),finest(105),
     &                 wrkstr(lenwrk),abserr,rel
      save neval,absest,finest
      integer n,maxcls,nfun,ifail

      rel = 0
      call dcutri(tstsub,nfun,vertce,numtri,0,maxcls,abserr,rel,lenver,
     &           lenwrk,irest,finest,absest,neval,ifail,wrkstr,iwork)
      write (*,99999)
99999 FORMAT (' '/5X,'DCUTRI TEST ')
      WRITE (*,99998) NEVAL,IFAIL
99998 FORMAT (5X,'SUBROUTINE CALLS = ',I6,', IFAIL = ',I2)
      WRITE (*,99997)
99997 FORMAT (7X,'N',3X,'ABSOLUTE ERROR',4X,'INTEGRAL')
      do n = 1,nfun
          write (*,99996) n,abs(finest(n)-1),finest(n)
99996     format (6x,i3,e14.2,f21.16)
      end do

      return
      end
      subroutine ftestp(z,nfun,f)

c*********************************************************************72
c
cc FTESTP: selected monomials, degree 0-13
C
      implicit none

      integer isave,nfun,m,k,l,index
      double precision z(2),f(nfun),ex(120),exact

      save isave,ex

      data isave/0/

      if (isave.eq.0) then
          do m = 0,13
              do k = 0,m
                  l = m - k
                  index = (m* (m+1))/2 + k + 1
                  ex(index) = exact(k,l)
              end do
          end do
          isave = 1
      end if

      do m = 0,13
          do k = 0,m
              l = m - k
              index = (m* (m+1))/2 + k + 1
              f(index) = z(1)**k*z(2)**l/ex(index)
          end do
      end do

      return
      end
      subroutine ftestq(z,nfun,f)

c*********************************************************************72
c
cc FTESTQ: selected monomials in polar coordinates, degree 0-13
C
      implicit none

      integer nfun,j,k,i,l
      double precision z(2),f(nfun),m(0:15,0:15),r,alpha,pi
      pi = 3.1415926535897932384d0
      m(0,0) = 1.000000000000000000000000000000d0
      m(2,0) = 0.250000000000000000000000000000d0
      m(4,0) = 0.099999999999999999999999999997d0
      m(6,0) = 0.051785714285714285714285714284d0
      m(8,0) = 0.031428571428571428571428571427d0
      m(10,0) = 0.021103896103896103896103896103d0
      m(12,0) = 0.015163408020550877693734836591d0
      m(14,0) = 0.011429195804195804195804195802d0
      m(3,3) = -0.099999999999999999999999999997d0
      m(5,3) = -0.057142857142857142857142857141d0
      m(7,3) = -0.035714285714285714285714285713d0
      m(9,3) = -0.024025974025974025974025974024d0
      m(11,3) = -0.017132867132867132867132867132d0
      m(13,3) = -0.012787212787212787212787212785d0
      m(15,3) = -0.009891579009226068049597461361d0
      m(6,6) = 0.035714285714285714285714285713d0
      m(8,6) = 0.024999999999999999999999999997d0
      m(10,6) = 0.018181818181818181818181818181d0
      m(12,6) = 0.013686313686313686313686313685d0
      m(14,6) = 0.010614385614385614385614385614d0
      m(9,9) = -0.018181818181818181818181818181d0
      m(11,9) = -0.013986013986013986013986013985d0
      m(13,9) = -0.010989010989010989010989010988d0
      m(15,9) = -0.008807369101486748545572074984d0
      m(12,12) = 0.010989010989010989010989010988d0
      m(14,12) = 0.008928571428571428571428571429d0
      m(15,15) = -0.007352941176470588235294117647d0
      i = 0
      r = sqrt(abs(z(1))**2+abs(z(2))**2)

      if (r.eq.0) then
          alpha = 0
      else
          alpha = acos(z(1)/r)
      end if

      if (z(2).lt.0) then
          alpha = 2*pi - alpha
      end if

      do 220 j = 0,13
          do 210 k = 0,j/3
              if (mod(j+3*k,2).eq.0) then
                  i = i + 1
                  if (r.eq.0) then
                      if (j.eq.0) then
                          f(i) = 4/ (3*sqrt(3.d0))
                      else
                          f(i) = 0
                      end if
                  else
                      if (j.eq.0) then
                          f(i) = 1
                      else
                          f(i) = 1
                          do l = 1,j
                              f(i) = f(i)*r
                          end do
                      end if
                      f(i) = f(i)*cos(3*k*alpha)
                      f(i) = 4*f(i)/ (3*sqrt(3.d0)*m(j,3*k))
                  end if
              end if
  210     continue
  220 continue

      return
      end
      subroutine ftesto(z,nfun,f)

c*********************************************************************72
c
cc FTEST0: Product Peak
C
      implicit none

      integer nfun
      double precision z(2),f(nfun)
      double precision alpha1,alpha2,beta1,beta2,tint

      f(1) = 1.d0
      alpha1 = 10
      alpha2 = 10
      beta1 = 0.2d0
      beta2 = 0.3d0
      tint = 1
      tint = tint* (atan((1-beta1)*alpha1)-atan( - beta1*alpha1))*
     &       alpha1
      tint = tint* (atan((1-beta2)*alpha2)-atan( - beta2*alpha2))*
     &       alpha2
      f(1) = f(1)/ (1.d0/alpha1**2+ (z(1)-beta1)**2)
      f(1) = f(1)/ (1.d0/alpha2**2+ (z(2)-beta2)**2)
      f(1) = f(1)/tint

      return
      end
      double precision function exact(k,l)

c*********************************************************************72
c
cc EXACT
C
      implicit none

      integer k,l,i
      double precision bin,binom

      exact = 0
      do i = 0,l + 1
          bin = binom(l+1,i)
          exact = exact + (-1)** (i)*bin/dble(i+k+1)
      end do
      exact = exact/dble(l+1)

      return
      end
      double precision function binom(n,r)

c*********************************************************************72
c
cc BINOM
c
      implicit none

      integer n,r,i

      binom = 1
      if (r.eq.0 .or. r.eq.n) then
          go to 999
      end if
      do 10 i = n,n - r + 1,-1
          binom = binom*i
   10 continue
      do 20 i = 1,r
          binom = binom/i
   20 continue

  999 return
      end
