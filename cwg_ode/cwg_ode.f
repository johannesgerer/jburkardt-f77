      subroutine srkode ( n, t, tout, y, hmin, eps, error, kflag, 
     &  jstart, yprime )

c*********************************************************************72
c
cc SRKODE applies a Runge-Kutta method to differential equations.
c
c  Discussion:
c
c    The routine applies a Runge-Kutta method of order four to integrate
c    a system of first order differential equations.
c
c    You must have an external statement in your program which calls 
c    this routine, declaring the name of your derivative routine which
c    you will be using.
c
c    You must have dimension statements in your program for ERROR, and 
c    Y, giving them at least as much space as required.
c
c    You must set the values of T and TOUT and the corresponding values of
c    Y before calling.
c
c    You should set a positive value to EPS, HMIN, and N.
c    You should set JSTART = 0 on first call.
c
c    Thereafter, if the results of the current step are satisfactory,
c    don't change anything except TOUT, which you should increase
c    (or decrease) to the next value at which a solution is desired.
c    Then just call the program again to get the next solution value.
c
c  Reference:
c
c    Charles Gear,
c    Numerical Initial Value Problems in Ordinary Differential Equations,
c    Prentice-Hall, 1971,
c    ISBN: 0136266061,
c    LC: QA372.G4.
c
c  Parameters:
c
c    Input, integer N, the number of first order differential equations 
c    to be solved.  N must be greater than 0 and less than or equal
c    to 20.
c
c    Input/output, real T, the current value of the independent variable.
c
c    Input, real TOUT, a point at which the solution Y is desired.  SRKODE will
c    probably compute several intermediate points before it returns with the 
c    value of the solution at this point.  TOUT must be at least HMIN away 
c    from T.
c
c  y      - a vector of dimension n which contains the values of the
c           dependent variables.
c
c  hmin   - the magnitude (ignoring sign) of the smallest stepsize
c           that the user will allow the code to take.
c           hmin should be positive.
c
c  eps    - the error test constant.  the code estimates the error
c           that it makes.  for each component of the solution, the
c           estimated error is required to be less than eps*ymax(i),
c           where ymax(i) is the larger of 1.0 and the largest
c           magnitude of y(i) seen so far.  eps should be positive.
c
c  error  - a vector of dimension n, which contains the estimated
c           single step error for each component.
c
c  kflag  - a completion code with the following meanings-
c
c           if kflag= -1 on return, the smallest possible
c           step (abs(h)=hmin) was taken, though the requested error
c           was not achieved.  integration may proceed, but the results
c           are not guaranteed.
c
c           if kflag= +1 on return, the requested error was achieved,
c           with a stepsize h of magnitude at least hmin.
c
c  jstart - an initialization indicator set by the user with the
c           following meanings-
c
c           jstart=  0, this is the first call to srkode for this
c           problem.  always call first with jstart=0 so that srkode
c           can initialize itself.
c
c           jstart= +1, the results of the last step are acceptable,
c           please take the next step.
c
c  yprime - an external quantity (that is, the name of a subroutine
c           which is passed as an argument to this subroutine).
c           it is the name of the user-provided subroutine which
c           evaluates the vector of derivatives that define the ode.
c           the subroutine may have any name the user likes.  that
c           name must appear in an external statement in the
c           calling program, and in the yprime position in the
c           call to srkode.  the subroutine must have the following
c           form (for this example, we assume the name of the subroutine
c           is diffun and that there are 10 ode's.)
c
c           subroutine diffun(t,y,dy,n)
c           dimension y(n)
c           dimension dy(n)
c
c           (for the given t and y that are input, evaluate dy)
c           dy(1)=...
c           dy(2)=...
c
c           dy(10)=...
c           return
c           end
c
      integer n

      real dy(20)
      real error(n)
      real h
      real hmin
      real t
      real tout
      real y(n)
      real ymax(20)
      real y1(20)
      real y2(20)
      real y3(20)
      real y4(20)
      real y5(20)
      real y6(20)
      real y7(20)
      real y8(20)
      real y9(20)
      external yprime

      save h

      data h / 0.0E+00 /
c
c  Do a few checks
c
      if ( n .le. 0 .or. 20 .lt. n ) then
        write(6,1000)
        write(6,1010)
        stop
      end if

      if ( hmin .le. 0.0e+00 ) then
        write(6,1030)
        write(6,1010)
        stop
      end if

      deltat = tout - t

      if ( jstart .eq. 0 ) then
        h = sign ( hmin, deltat )
      end if

      if ( abs ( deltat ) .lt. hmin ) then
        write(6,1020)
        write(6,1010)
        stop
      end if

      if ( eps .le. 0.0e+00 ) then
        write(6,1040)
        write(6,1010)
        stop
      end if
c
c  Save sign of h as a direction
c
      direct = sign ( 1.0e+00, h )

      if ( jstart .lt. 0 ) then

        do i = 1, n
          y(i) = y1(i)
          ymax(i) = y3(i)
        end do

      else 

        if ( 0 .eq. jstart) then

          do i = 1, n
            ymax(i) = 1.0e+00
          end do

          jstart = 1

        end if

        do i = 1, n
          y1(i) = y(i)
          y3(i) = ymax(i)
        end do

        call yprime ( t, y, y2, n )

      end if
c
c  Take next step towards TOUT.
c
40    continue

      if ( h .gt. 0.0e+00 .and. ( t + h ) .gt. tout ) then
        h = tout - t
      else if ( h .lt. 0.0e+00 .and. ( t + h ) .lt. tout ) then
        h = tout - t
      end if

      kflag = 1
c
c  save the final value of t and calculate the half step
c
  110 continue

      a = h + t
      hhalf = h * 0.5e+00
c
c  Perform one full runge kutta step
c
      call srk1(n,t,y1,y2,y7,h,y4,y8,y9,yprime)
c
c  now perform two half interval runge kutta steps
c
      call srk1(n,t,y1,y2,y7,hhalf,y5,y8,y9,yprime)
      thalf=t+hhalf
      call yprime(thalf,y5,dy,n)
      call srk1(n,thalf,y5,dy,y7,hhalf,y6,y8,y9,yprime)
c
c  Calculate the new maximum y's, the errors, and the maximum
c  relative errors.
c
      errmax = 0.0e+00

      do i = 1, n
        ymax(i) = max ( ymax(i), abs ( y4(i) ) )
        ymax(i) = max ( ymax(i), abs ( y5(i) ) )
        ymax(i) = max ( ymax(i), abs ( y6(i) ) )
        error(i) = abs ( ( y6(i) - y4(i) ) / 15.0e+00 )
        errmax = max ( errmax, error(i) / ( eps * ymax(i) ) )
        y(i) = ( 16.0e+00 * y6(i) - y4(i) ) / 15.0e+00
      end do

      if ( errmax .eq. 0.0e+00 ) then
        h = h * 2.0e+00
      end if

      if ( 0.0 .lt. errmax ) then
        h = h * ( errmax**(-0.2e+00) ) * 0.99e+00
      endif

      if ( 1.0e+00 .lt. errmax ) then

        if ( hmin .lt. abs(h) ) then
          go to 110
        end if

        if ( kflag .lt. 0 ) then
          kflag = -1
          return
        end if

        h = direct * hmin
        kflag = -1
        go to 110

      end if

      kflag = 1

      t=a
      if(h.gt.0.0e+00.and.t.ge.tout)return
      if(h.lt.0.0e+00.and.t.le.tout)return

      go to 40

 1000 format(' fatal error - n.le.0 or n.gt.20')
 1010 format(' this program is forcing a stop')
 1020 format(' fatal error - (tout-t).lt.hmin')
 1030 format(' fatal error - hmin.le.0.0')
 1040 format(' fatal error - eps.le.0.0')
      end
      subroutine srk1 ( n, t, y, dy, dy1, h, y1, y2, y3, yprime )

c*********************************************************************72
c
cc SRK1 performs one Runge Kutta step.  
c
c  Reference:
c
c    Charles Gear,
c    Numerical Initial Value Problems in Ordinary Differential Equations,
c    Prentice-Hall, 1971,
c    ISBN: 0136266061,
c    LC: QA372.G4.
c
c  Parameters:
c
c    Input, integer N, the number of equations.
c
c    Input, real T, the initial value of the independent variable.
c
c    y      initial value of the dependent variables.
c
c    dy     initial value of the derivatives.
c
c    dy1    work space
c
c    h      step size.
c
c    y1     the value of y at t+h as computed by this program.
c
c    y2     work space
c
c    y3     work space
c
c    yprime the external name of the derivative-computing routine.
c
      implicit none

      integer n

      real dy(n)
      real dy1(n)
      real h
      real hhalf
      integer i
      real t
      real y(n)
      real y1(n)
      real y2(n)
      real y3(n)
      external yprime

      hhalf = h * 0.5e+00

      do i = 1, n
        y2(i) = y(i) + hhalf * dy(i)
      end do

      call yprime ( t + hhalf, y2, dy1, n )

      do i = 1, n
        y3(i) = y(i) + hhalf * dy1(i)
        y2(i) = y2(i) + 2.0e+00 * y3(i)
      end do

      call yprime ( t + hhalf, y3, dy1, n )

      do i = 1, n
        y3(i) = y(i) + h * dy1(i)
        y2(i) = y2(i) + y3(i)
      end do

      call yprime ( t + h, y3, dy1, n )

      do i = 1, n
        y1(i) = ( y2(i) - y(i) + hhalf * dy1(i) ) / 3.0e+00
      end do

      return
      end
      subroutine sprode ( n, t, tout, y, hmin, eps, mf, error, kflag,
     & jstart, yprime )

c*********************************************************************72
c
cc SPRODE uses polynomial or rational approximation to solve ODE's.
c
c  Discussion:
c
c    This routine uses polynomial or rational function approximation
c    to solve a system of first order differential equations.
c
c  Reference:
c
c    Charles Gear,
c    Numerical Initial Value Problems in Ordinary Differential Equations,
c    Prentice-Hall, 1971,
c    ISBN: 0136266061,
c    LC: QA372.G4.
c
c  Parameters:
c
c    Input, integer N, the number of first order differential equations.
c    N must be greater than 0 and less than or equal to 20.
c
c  t      - the current value of the independent variable.
c
c  tout   - the next value of t at which a solution is desired.
c           sprode will probably compute several points on its way
c           to tout.  abs(tout-t) must be greater than hmin.
c
c  y      - a vector of dimension n containing the current values
c           of the dependent variables.
c
c  hmin   - the minimum stepsize that the user will allow
c           the subroutine to take.  hmin must be greater than zero.
c
c  eps    - the error test constant.  the estimated errors
c           are required to be less than eps*ymax(i) in each
c           component where ymax(i) is the larger of 1.0 and
c           the maximum magnitude of y(i) seen so far.
c           eps must be positive.
c
c  mf     - the method indicator which is to be selected by
c           the user.  the following are available-
c
c           mf=0 bulirsch-stoer rational extrapolation
c           mf=1 polynomial extrapolation
c
c  error  - a vector of size n, containing an estimate of the
c           error in the i-th solution component.
c
c  kflag  - a completion code set by the program just
c           before it returns.  it has the following meanings-
c
c           kflag=-1  the step was taken with abs(h)=hmin,
c           but the requested error was not achieved.
c           integration may proceed but the results are not
c           guaranteed.
c
c           kflag=+1  the step was successful.
c
c  jstart - an initialization flag set by the user on first call,
c           and on any call where the user felt the last step
c           was unsatisfactory -
c
c           jstart=0 this is the first call to sprode for
c           this problem.
c
c           jstart=+1 this is not the first call.  please take a new
c           step, since the last one was acceptable.
c
c  yprime - an external quantity (that is, the name of a subroutine
c           which is passed as an argument to this subroutine).
c           it is the name of the user-provided subroutine which
c           evaluates the vector of derivatives that define the ode.
c           the subroutine may have any name the user likes.  that
c           name must appear in an external statement in the
c           calling program, and in the yprime position in the
c           call to sprode.  the subroutine must have the following
c           form (for this example, we assume the name of the subroutine
c           is diffun and that there are 10 ode's.)
c
c           subroutine diffun(t,y,dy,n)
c           dimension y(n)
c           dimension dy(n)
c
c           (for the given t and y that are input, evaluate dy)
c           dy(1)=...
c           dy(2)=...
c
c           dy(10)=...
c           return
c           end
c
c
c  calling the program
c
c
c  you must have an external statement in the program which calls
c  this program, declaring the name of your derivative routine which
c  you will be using.
c
c  you must have dimension statements in your program for
c  error, and y, giving them at least as much space
c  as described above.
c
c  you must set the values t and tout, and the corresponding values of
c  y before calling.
c  you should choose mf=0 or mf=1.
c
c  you should set a positive value to eps, hmin, and n.
c  you should set jstart=0 on first call.
c
c  thereafter, if the results of the current step are satisfactory,
c  don't change anything except the value of tout, which you should
c  increase (or decrease) to the next value at which a solution
c  is desired. then just call the program again to get the
c  next solution value.
c
      integer n

      real dy(20)
      real dyn(20)
      real error(n)
      real extrap(20,11)
      integer maxord
      integer maxpts
      real quot(11,2)
      real y(n)
      real ymax(20)
      real ymaxhv(20,12)
      real ymaxsv(20)
      real yn(20)
      real ynhv(20,12)
      real ynm1(20)
      real ynm1hv(20,12)
      external yprime
      real ysave(20)
c
c  the arrays are used for the following data
c
c  ysave  - the initial values of y are saved for a restart.
c  ynm1   - y(n-1), the previous value of y in the midpoint
c           method.
c  yn     - y(n), the current value of y in the midpoint
c           integration.
c  dyn    - the initial value of the derivative of y.
c  ymaxsv - the saved values of ymax at the initial point.
c  quot   - the quotients (h(i)/h(i+m))**2 used in
c           the extrapolation.
c  extrap - the most recent extrapolated values of y in the
c           case of polynomial extrapolation, or of the
c           differences in the case of rational function
c           extrapolation.
c  ynm1hv - the values of ynm1 at the midpoint of the basic
c           interval if the number of substeps is divisible
c           by 4.  this information is used to avoid redoing
c           the integration in case the step is halved.
c  ynhv   - the similar values of yn.
c  ymaxhv - and the same for ymax.
c  error  - the estimates of the single step error are saved here.
c
      if(jstart.ne.0)go to 70
      maxord=11
      maxpts=11
      quot(1,1)=1.0e+00
      quot(2,1)=2.25e+00
      quot(3,1)=4.0e+00
      quot(4,1)=9.0e+00
      quot(5,1)=16.0e+00
      quot(6,1)=36.0e+00
      quot(7,1)=64.0e+00
      quot(8,1)=144.0e+00
      quot(9,1)=256.0e+00
      quot(10,1)=576.0e+00
      quot(11,1)=1024.0e+00
      quot(1,2)=1.0e+00
      quot(2,2)=1.777777777777777e+00
      quot(3,2)=4.0e+00
      quot(4,2)=7.111111111111111e+00
      quot(5,2)=16.0e+00
      quot(6,2)=28.444444444444444e+00
      quot(7,2)=64.0e+00
      quot(8,2)=113.777777777777777e+00
      quot(9,2)=256.0e+00
      quot(10,2)=455.111111111111111e+00
      quot(11,2)=1024.0e+00
c
c  fmax is a number smaller than the first integer that cannot
c  be represented exactly in floating point.
c
      fmax=10000000.0e+00

      do i=1,n
        ymax(i)=1.0e+00
      end do
c
c  Do a few checks
c
      if ( n .le. 0 .or. n .gt. 20 ) then
        write(6,1010)
        write(6,1020)
        stop
      end if

      if(hmin.gt.0.0e+00)go to 30
      write(6,1040)
      write(6,1020)
      stop
   30 continue
      deltat=tout-t
      if(jstart.eq.0)h=sign(hmin,deltat)
      if(abs(deltat).ge.hmin)go to 40
      write(6,1030)
      write(6,1020)
      stop
   40 continue
      if(eps.gt.0.0e+00)go to 50
      write(6,1050)
      write(6,1020)
      stop
   50 continue
      if(mf.eq.0.or.mf.eq.1)go to 60
      write(6,1000)
      write(6,1020)
      stop
   60 continue
      jstart=1
   70 continue
   80 continue
      if(jstart.lt.0)go to 100
c
c  save the values of y and ymax in case a restart is necessary
c
      do 90 i=1,n
        ysave(i)=y(i)
        ymaxsv(i)=ymax(i)
   90   continue
      call yprime(t,y,dyn,n)
      go to 120
c
c  restore the values of y and ymax for a restart
c
  100 continue
      do 110 i=1,n
        y(i)=ysave(i)
        ymax(i)=ymaxsv(i)
  110   continue
  120 continue
c
c  the following counters and switches are used
c
c  j      - is the count through the different sub steps g used.
c  jodd   - is 1 if j is odd, 2 if j is even
c  jhvsv  - is the number of substep sizes for which half way
c           information has been saved.
c  jhvsv1 - the value of jhvsv from the previous cycle
c  m      - the number of pairs of sub steps which make up the step h.
c           m takes the sequence 1,2,3,4,6,8,12,16, etc.
c  mnext  - the next value of m.
c  mtwo   - the next but one value of m.
c  quotsv - the last value of quot is irregular due to the fact that
c           the sequence by the multiples 9/4, 16/9 (odd) or
c           16/9, 9/4 (even) until the final multiple of 4.  however,
c           (h(0)/h(m))**2 is always m**2.  the regular value of
c           quot is saved in quotsv, and replaced by m**2.
c  konv   - is set to +1 initially, and reset to -1 if the error
c           test fails.
c
c  use next step to try to reach tout.
c
  130 continue
      if(h.gt.0.0e+00.and.(t+h).gt.tout)h=tout-t
      if(h.lt.0.0e+00.and.(t+h).lt.tout)h=tout-t
      jhvsv1=0
      kflag=1
  140 continue
      jhvsv=0
  150 continue
      a=h+t
      jodd=1
      m=1
      mnext=2
      mtwo=3
      do 310 j=1,maxpts
        quotsv=quot(j,jodd)
        quot(j,jodd)=m*m
        konv=1
        if(j.le.(maxord/2))konv=-1
        if(j.le.(maxord+1))go to 160
        l=maxord+1
        hchnge=0.7071068e+00*hchnge
        go to 170
  160   continue
        l=j
        hchnge=1.0+real(maxord+1-j)/6.0
  170   continue
        b=h/real(m)
        g=b*0.5e+00
        if(j.gt.jhvsv1)go to 190
c
c  the values of the midpoint integration were saved at the
c  half way point in the previous integration.  use them.
c
        do 180 i=1,n
          yn(i)=ynhv(i,j)
          ynm1(i)=ynm1hv(i,j)
          ymax(i)=ymaxhv(i,j)
  180     continue
        go to 240
c
c  integrate over the range h by 2*m steps of a midpoint method
c
  190   continue
        do 200 i=1,n
          ynm1(i)=ysave(i)
          yn(i)=ysave(i)+g*dyn(i)
          ymax(i)=ymaxsv(i)
  200     continue
        m2=m+m
        tu=t
        do 230 k=2,m2
          tu=tu+g
          call yprime(tu,yn,dy,n)
          do 210 i=1,n
            u=ynm1(i)+b*dy(i)
            ynm1(i)=yn(i)
            yn(i)=u
            u=abs(u)
            if(u.gt.ymax(i))ymax(i)=u
  210       continue
          if((k.ne.m).or.(jhvsv1.ne.0).or.(k.eq.3))go to 230
          jhvsv=jhvsv+1
          do 220 i=1,n
            ynhv(i,jhvsv)=yn(i)
            ynm1hv(i,jhvsv)=ynm1(i)
            ymaxhv(i,jhvsv)=ymax(i)
  220       continue
  230     continue
  240   continue
        call yprime(a,yn,dy,n)
        do 300 i=1,n
          v=extrap(i,1)
c
c  calculate the final value to be used in the extrapolation process
c
          ta=(yn(i)+ynm1(i)+g*dy(i))*0.5e+00
          c=ta
c
c  insert the integral as the first extrapolated value
c
          extrap(i,1)=ta
          if(l.lt.2)go to 290
          if(abs(v)*fmax.lt.abs(c))go to 350
          if(mf.gt.0)go to 270
c
c  perform the extrapolation by rational functions on the second
c  and subsequent intervals
c
          do 260 k=2,l
            b1=quot(k,jodd)*v
            b=b1-c
            u=v
            if(b.eq.0.0e+00)go to 250
            b=(c-v)/b
            u=c*b
            c=b1*b
  250       continue
            v=extrap(i,k)
            extrap(i,k)=u
            ta=ta+u
  260       continue
          go to 290
  270     continue
          do 280 k=2,l
            ta=ta+(ta-v)/(quot(k,jodd)-1.0e+00)
            v=extrap(i,k)
            extrap(i,k)=ta
  280       continue
          go to 290
  290     continue
          u=abs(ta)
          if(u.gt.ymax(i))ymax(i)=u
          error(i)=abs(y(i)-ta)
          y(i)=ta
          if(error(i).gt.eps*ymax(i))konv=-1
  300     continue
        quot(j,jodd)=quotsv
        if(konv.gt.0)go to 330
        jodd=3-jodd
        m=mnext
        mnext=mtwo
        mtwo=m+m
  310   continue
      jhvsv1=jhvsv
  320 continue
      if(abs(h).le.hmin)go to 340
      h=h*0.5e+00
      if(abs(h).ge.hmin)go to 140
      h=sign(hmin,h)
      go to 130
  330 continue
      h=h*hchnge
      t=a
      if(h.gt.0.0e+00.and.t.ge.tout)return
      if(h.lt.0.0e+00.and.t.le.tout)return
      go to 80
  340 continue
      kflag=-1
      go to 330
  350 continue
      quot(j,jodd)=quotsv
      go to 320
 1000 format(' fatal error - mf.ne.0 and mf.ne.1')
 1010 format(' fatal error - n.le.0 or n.gt.20')
 1020 format(' this program is forcing a stop')
 1030 format(' fatal error - (tout-t).lt.hmin')
 1040 format(' fatal error - hmin.le.0.0')
 1050 format(' fatal error - eps.le.0.0')
      end
      subroutine smvode ( n, t, tout, y, hmin, eps, mf, error, kflag,
     & jstart, yprime, party )

c*********************************************************************72
c
cc SMVODE solves an ODE with a multi-value method.
c
c  Discussion:
c
c    This routine is a multi-value method for solving a system of
c    first order ordinary differential equations.  this is
c    the single precision version.
c
c    this subroutine integrates a set of n ordinary differential first
c    order equations using a multivalue method.
c
c  Reference:
c
c    Charles Gear,
c    Numerical Initial Value Problems in Ordinary Differential Equations,
c    Prentice-Hall, 1971,
c    ISBN: 0136266061,
c    LC: QA372.G4.
c
c  Parameters:
c
c  n      - the number of first order differential equations.
c           n must be greater than 0 and less than or equal to 20.
c
c  t      - the current value of the independent variable.
c
c  tout   - the next value of the independent variable at which
c           a solution is desired.
c           abs(tout-t) must be greater than hmin.
c
c  y      - an array of n entries containing the dependent variables
c
c  hmin   - the minimum step size that will be used for the
c           integration.  note that on starting this must be
c           much smaller than the average h expected since
c           a first order method is used initially.
c
c  mf     - the method indicator.  the following are allowed..
c
c           0   an adams predictor corrector is used.
c
c           1   a multi-step method suitable for stiff
c           systems is used.  it will also work for
c           non stiff systems.  however the user
c           must provide a subroutine which
c           evaluates the partial derivatives of
c           the differential equations with respect
c           to the y's.  see description of party below.
c
c           2   the same as case 1, except that the
c           subroutine computes the partial
c           derivatives by numerical differencing
c           of the derivatives.
c
c  error  - an array of n elements which contains the estimated
c           one step error in each component.
c
c  kflag  - a completion code with the following meanings..
c
c           +1   the step was successful.
c           -1   the step was taken with h + hmin, but the
c           requested error was not achieved.
c           -2  the maximum order specified was found to
c           be too large.
c           -3   corrector convergence could not be
c           achieved for h .gt. hmin.
c           -4   the requested error is smaller than can
c           be handled for this problem.
c
c  jstart - an input indicator with the following meanings..
c
c           0   perform the first step.  the first step
c           must be done with this value of jstart
c           so that the subroutine can initialize itself.
c
c           +1 (or any positive value)   take a new step continuing from
c           the last.  note that on return, jstart is set to nq,
c           the current order of the method
c           at exit.  nq is also the order of the maximum
c           derivative available.
c
c  yprime - an external quantity (that is, the name of a subroutine
c           which is passed as an argument to this subroutine).
c           it is the name of the user-provided subroutine which
c           evaluates the vector of derivatives that define the ode.
c           the subroutine may have any name the user likes.  that
c           name must appear in an external statement in the
c           calling program, and in the yprime position in the
c           call to smvode.  the subroutine must have the following
c           form (for this example, we assume the name of the subroutine
c           is diffun and that there are 10 ode's.)
c
c           subroutine diffun(t,y,dy,n)
c           dimension y(n)
c           dimension dy(n)
c
c           (for the given t and y that are input, evaluate dy)
c           dy(1)=...
c           dy(2)=...
c
c           dy(10)=...
c           return
c           end
c
c  party  - an external quantity (that is, the name of a subroutine
c           which is passed as an argument to this subroutine).
c           if mf=1, you must write a subroutine as described below.
c           however, if mf=0 or mf=2, you don't have to write the
c           subroutine, but you still must supply a name in the
c           party position in the call to smvode.  use the
c           name 'pdummy' in that case.  in the case mf=1, party will be
c           the name of the user-provided subroutine which
c           evaluates the partial derivatives of the right hand
c           sides rhs(i) of the ode's with respect to each component
c           y(j) and stores them in dermat(i,j).
c           the subroutine may have any name the user likes.  that
c           name must appear in an external statement in the
c           calling program, and in the party position in the
c           call to smvode.  the subroutine must have the following
c           form (for this example, we assume the name of the subroutine
c           is parfun.)
c
c           subroutine parfun(t,y,dermat,n)
c           dimension dermat(n,n)
c           dimension y(n)
c
c           (for the given t and y, evaluate d rhs(i)/d y(j) )
c
c           do 20 i=1,n
c             do 10 j=1,n
c               dermat(i,j)=.....  (=d rhs(i)/d y(j) )
c           10 continue
c           20 continue
c           return
c           end
c
c
c  calling the program
c
c
c  you must have an external statement in the program which calls
c  this program, declaring the names of your derivative routine which
c  you will be using, and the partial derivative routine (or the name
c  of the dummy routine pdummy).
c
c  you must have dimension statements in your program for
c  error, and y, giving them at least as much space
c  as described above.
c
c  you must set the values of t and tout, and the corresponding values
c  of y before calling.
c
c  you should set a positive value to eps, hmin, and n.
c  you should set jstart=0 on first call.
c
c  thereafter, if the results of the current step are satisfactory,
c  don't change anything except tout, which you should increase
c  (or decrease) to the next value of t at which a solution is
c  desired.  then just call the program again to get the
c  next solution value.
c
      real a(8)
      real dermat(400)
      real error(n)
      real h
      real hnew
      integer ipivot(20)
      integer nq
      external party
      real pertst(7,2,3)
      real save(20,8)
      real save9(20)
      real save10(20)
      real save11(20)
      real save12(20)
      real y(n)
      real ytable(20,8)
      real ymax(20)
      external yprime

      save h
      save hnew
      save nq
      save pertst
      save save
      save save9
      save save10
      save save11
      save save12
      save ytable
      save ymax

      data h / 0.0E+00 /
      data hnew / 0.0E+00 /
      data nq / 0 /
c
c  do a few checks
c
      if(n.gt.0.and.n.le.20)go to 10
      write(6,1000)
      write(6,1010)
      stop
   10 continue
      if(hmin.gt.0.0e+00)go to 20
      write(6,1030)
      write(6,1010)
      stop
   20 continue
      deltat=tout-t
      if(jstart.eq.0)h=sign(hmin,deltat)
      if(abs(deltat).ge.hmin)go to 30
      write(6,1020)
      write(6,1010)
      stop
   30 continue
      if(eps.gt.0.0e+00)go to 40
      write(6,1040)
      write(6,1010)
      stop
   40 continue
      if(mf.ge.0.and.mf.le.2)go to 50
      write(6,1050)
      write(6,1010)
      stop
   50 continue
      if(mf.eq.0)maxder=7
      if(mf.eq.1.or.mf.eq.2)maxder=6
      do i = 1, n
        ytable(i,1) = y(i)
      end do
   70 continue
      iret=1
      kflag=1
      if(h.gt.0.0e+00.and.(t+h).gt.tout)h=tout-t
      if(h.lt.0.0e+00.and.(t+h).lt.tout)h=tout-t
      if (jstart.le.0) go to 120
c
c  begin by saving information for possible restart and changing
c  h by the factor r if the caller has changed h.  all variables
c  dependent on h must also be changed.
c  e is a comparison for errors of the current order nq.  eup is
c  to test for increasing the order, edwn for decreasing the order.
c  hnew is the step size that was used on the last call.
c
   80 do i = 1,n
        do j = 1,k
          save(i,j) = ytable(i,j)
        end do
      end do
      hold = hnew
      if (h.eq.hold) go to 110
  100 racum = h/hold
      iret1 = 1
      go to 820
  110 nqold = nq
      told = t
      racum = 1.0e+00
      if (jstart.gt.0) go to 330
      go to 150
  120 if (jstart.eq.-1) go to 140
c
c  on the first call, the order is set to 1 and the initial
c  derivatives are calculated.
c
      pertst(1,1,1)=2.0e+00
      pertst(2,1,1)=4.5e+00
      pertst(3,1,1)=7.333e+00
      pertst(4,1,1)=10.42e+00
      pertst(5,1,1)=13.7e+00
      pertst(6,1,1)=17.15e+00
      pertst(7,1,1)=1.0e+00
      pertst(1,2,1)=2.0e+00
      pertst(2,2,1)=12.0e+00
      pertst(3,2,1)=24.0e+00
      pertst(4,2,1)=37.89e+00
      pertst(5,2,1)=53.33e+00
      pertst(6,2,1)=70.08e+00
      pertst(7,2,1)=87.97e+00
      pertst(1,1,2)=3.0e+00
      pertst(2,1,2)=6.0e+00
      pertst(3,1,2)=9.167e+00
      pertst(4,1,2)=12.5e+00
      pertst(5,1,2)=15.98e+00
      pertst(6,1,2)=1.0e+00
      pertst(7,1,2)=1.0e+00
      pertst(1,2,2)=12.0e+00
      pertst(2,2,2)=24.0e+00
      pertst(3,2,2)=37.89e+00
      pertst(4,2,2)=53.33e+00
      pertst(5,2,2)=70.08e+00
      pertst(6,2,2)=87.97e+00
      pertst(7,2,2)=1.0e+00
      pertst(1,1,3)=1.0e+00
      pertst(2,1,3)=1.0e+00
      pertst(3,1,3)=0.5e+00
      pertst(4,1,3)=0.1667e+00
      pertst(5,1,3)=0.04133e+00
      pertst(6,1,3)=0.008267e+00
      pertst(7,1,3)=1.0e+00
      pertst(1,2,3)=1.0e+00
      pertst(2,2,3)=1.0e+00
      pertst(3,2,3)=2.0e+00
      pertst(4,2,3)=1.0e+00
      pertst(5,2,3)=0.3157e+00
      pertst(6,2,3)=0.07407e+00
      pertst(7,2,3)=0.0139e+00

      nq = 1

      call yprime(t,ytable,save11,n)

      do i = 1,n
        ytable(i,2) = save11(i)*h
        ymax(i)=1.0e+00
      end do
      hnew = h
      k = 2
      go to 80
c
c  repeat last step by restoring saved information.
c
  140 if (nq.eq.nqold) jstart = 1
      t = told
      nq = nqold
      k = nq + 1
      go to 100
c
c  set the coefficients that determine the order and the method
c  type.  check for excessive order.  the last two statements of
c  this section set iweval  .gt.0 if dermat is to be re-evaluated
c  because of the order change, and then repeat the integration
c  step if it has not yet been done (iret = 1) or skip to a final
c  scaling before exit if it has been completed (iret = 1).
c
  150 if (mf.eq.0) go to 160
      if (nq.gt.6) go to 170
      go to (250,260,270,280,290,300),nq
  160 if (nq.gt.7) go to 170
      go to (180,190,200,210,220,230,240),nq
  170 kflag = -2
      go to 860
c
c  the following coeficients should be defined to the maximum
c  accuracy permitted by the machine.  they are, in the order used..
c
c  -1
c  -1/2,-1/2
c  -5/12,-3/4,-1/6
c  -3/8,-11/12,-1/3,-1/24
c  -251/720,-25/24,-35/72,-5/48,-1/120
c  -95/288,-137/120,-5/8,-17/96,-1/4,-1/720
c  -19087/60480,-49/40,-203/270,-49/192,-7/144,-7/1440,-1/5040
c
c  -1
c  -2/3,-1/3
c  -6/11,-6/11,-1/11
c  -12/25,-7/10,-1/5,-1/50
c  -120/274,-225/274,-85/274,-15/274,-1/274
c  -180/441,-58/63,-15/36,-25/252,-3/252,-1/1764
c
  180 a(1) = -1.0e+00
      a(2)=-1.0e+00
      go to 310
  190 a(1) = -0.500000000e+00
      a(2)=-1.0e+00
      a(3) = -0.500000000e+00
      go to 310
  200 a(1) = -0.4166666666666667e+00
      a(2)=-1.0e+00
      a(3) = -0.750000000e+00
      a(4) = -0.1666666666666667e+00
      go to 310
  210 a(1) = -0.375000000e+00
      a(2)=-1.0e+00
      a(3) = -0.9166666666666667e+00
      a(4) = -0.3333333333333333e+00
      a(5) = -0.0416666666666667e+00
      go to 310
  220 a(1) = -0.3486111111111111e+00
      a(2)=-1.0e+00
      a(3) = -1.0416666666666667e+00
      a(4) = -0.4861111111111111e+00
      a(5) = -0.1041666666666667e+00
      a(6) = -0.008333333333333333e+00
      go to 310
  230 a(1) = -0.3298611111111111e+00
      a(2)=-1.0e+00
      a(3) = -1.141666666666667e+00
      a(4) = -0.625000000e+00
      a(5) = -0.1770833333333333e+00
      a(6) = -0.0250000000e+00
      a(7) = -0.001388888888888889e+00
      go to 310
  240 a(1) = -0.3155919312169312e+00
      a(2)=-1.0e+00
      a(3) = -1.225000000e+00
      a(4) = -0.7518518518518519e+00
      a(5) = -0.2552083333333333e+00
      a(6) = -0.04861111111111111e+00
      a(7) = -0.004861111111111111e+00
      a(8) = -0.0001984126984126984e+00
      go to 310
  250 a(1) = -1.000000000e+00
      a(2)=-1.0e+00
      go to 310
  260 a(1) = -0.6666666666666667e+00
      a(2)=-1.0e+00
      a(3) = -0.3333333333333333e+00
      go to 310
  270 a(1) = -0.5454545454545455e+00
      a(2)=-1.0e+00
      a(3) = a(1)
      a(4) = -0.09090909090909091e+00
      go to 310
  280 a(1) = -0.480000000e+00
      a(2)=-1.0e+00
      a(3) = -0.700000000e+00
      a(4) = -0.200000000e+00
      a(5) = -0.020000000e+00
      go to 310
  290 a(1) = -0.437956204379562e+00
      a(2)=-1.0e+00
      a(3) = -0.8211678832116788e+00
      a(4) = -0.3102189781021898e+00
      a(5) = -0.05474452554744526e+00
      a(6) = -0.0036496350364963504e+00
      go to 310
  300 a(1) = -0.4081632653061225e+00
      a(2)=-1.0e+00
      a(3) = -0.9206349206349206e+00
      a(4) = -0.4166666666666667e+00
      a(5) = -0.0992063492063492e+00
      a(6) = -0.0119047619047619e+00
      a(7) = -0.000566893424036282e+00
  310 k= nq+1
      idoub = k
      mtyp = (4 - mf)/2
      enq2 = 0.5/real(nq+1)
      enq3 = 0.5/real(nq+2)
      enq1 = 0.5/real(nq)
      pepsh = eps
      eup = (pertst(nq,mtyp,2)*pepsh)**2
      e = (pertst(nq,mtyp,1)*pepsh)**2
      edwn = (pertst(nq,mtyp,3)*pepsh)**2
      if (edwn.eq.0.0e+00) go to 850
      bnd = eps*enq3/real(n)
  320 iweval = mf
      go to (330,750),iret
c
c  this section computes the predicted values by effectively
c  multiplying the saved information by the pascal triangle
c  matrix.
c
  330 t = t + h
      do 340 j = 2,k
        do 340 j1 = j,k
          j2 = k - j1 + j -1
          do 340 i = 1,n
  340       ytable(i,j2) = ytable(i,j2) + ytable(i,j2+1)
c
c  Up to 3 corrector iterations are taken.  convergence is tested
c  by requiring changes to be less than bnd which is dependent on
c  the error test constant.
c  the sum of the corrections is accumulated in the array
c  error(1).  it is equal to the k-th derivative of y multiplied
c  by  h**k/(factorial(k-1)*a(k)), and is therefore proportional
c  to the actual errors to the lowest power of h present.  (h**k)
c
      do i = 1,n
        error(i) = 0.0e+00
      end do

      do 510 l = 1,3
        call yprime(t,ytable,save11,n)
c
c  If there has been a change of order or there has been trouble
c  with convergence, dermat is re-evaluated prior to starting the
c  corrector iteration in the case of stiff methods.  iweval is
c  then set to -1 as an indicator that it has been done.
c
        if (iweval.lt.1) go to 440
        if (mf.eq.2) go to 400
        call party(t,ytable,dermat,n)
        r = a(1)*h
        do i = 1,n
          do j=1,n
            index=n*(j-1)+i
            dermat(index)=dermat(index)*r
          end do
        end do
c
c  Add the identity matrix to the jacobian and invert to get dermat.
c
  380   doi = 1,n
          index=n*(i-1)+i
          dermat(index)=1.0e+00+dermat(index)
        end do

        iweval = -1
        call sgefa(dermat,n,n,ipivot,info)
        if(info.eq.0)j1=1
        if(info.ne.0)j1=-1
        if (j1.gt.0) go to 440
        go to 520
c
c  evaluate the jacobian into dermat by numerical differencing.  r
c  is the change made to the element of y.  it is eps relative to y with
c  a minimum of eps**2.
c
  400   do i = 1,n
          save9(i) = ytable(i,1)
        end do

        do j = 1,n
          r = eps*max(eps,abs(save9(j)))
          ytable(j,1) = ytable(j,1) + r
          d = a(1)*h/r
          call yprime(t,ytable,save12,n)
          do i = 1,n
            index=n*(j-1)+i
            dermat(index) = (save12(i) - save11(i))*d
          end do
          ytable(j,1) = save9(j)
        end do

        go to 380
  440   if (mf.ne.0) go to 460
        do 450 i = 1,n
  450     save9(i) = ytable(i,2) - save11(i)*h
        go to 490
  460   do 470 i = 1,n
  470     save12(i) = ytable(i,2) - save11(i)*h
        do 480 i = 1,n
          save9(i)=save12(i)
  480     continue
        job=0
        call sgesl(dermat,n,n,ipivot,save9,job)
  490   nt = n
c
c  correct and see if all changes are less than bnd relative to ymax.
c  if so, the correction is said to have converged.
c
        do 500 i = 1,n
          ytable(i,1) = ytable(i,1) + a(1)*save9(i)
          ytable(i,2) = ytable(i,2) - save9(i)
          error(i) = error(i) + save9(i)
          if (abs(save9(i)).le.(bnd*ymax(i))) nt = nt-1
  500     continue
        if (nt.le.0) go to 560
  510   continue
c
c  the corrector iteration failed to converge in 3 tries.  various
c  possibilities are checked for.  if h is already hmin and
c  this is either adams method or the stiff method in which the
c  matrix dermat has already been re-evaluated, a no convergence exit
c  is taken.  otherwise the matrix dermat is re-evaluated and/or the
c  step is reduced to try and get convergence.
c
  520 t = told
      if ((h.le.(hmin*1.00001e+00)).and.((iweval-mtyp).lt.-1)) go to 530
      if ((mf.eq.0).or.(iweval.ne.0)) racum = racum*0.25e+00
      iweval = mf
      iret1 = 2
      go to 820
  530 kflag = -3
  540 do 550 i = 1,n
        do 550 j = 1,k
  550     ytable(i,j) = save(i,j)
      h = hold
      nq = nqold
      jstart = nq
      go to 860
c
c  the corrector converged and control is passed to statement 520
c  if the error test is o.k., and to 540 otherwise.
c  if the step is o.k. it is accepted.  if idoub has been reduced
c  to one, a test is made to see if the step can be increased
c  at the current order or by going to one higher or one lower.
c  such a change is only made if the step can be increased by at
c  least 1.1.  if no change is possible idoub is set to 10 to
c  prevent further testing for 10 steps.
c  if a change is possible, it is made and idoub is set to
c  nq + 1    to prevent further testing for that number of steps.
c  if the error was too large, the optimum step size for this or
c  lower order is computed, and the step retried.  if it should
c  fail twice more it is an indication that the derivatives that
c  have accumulated in the y array have errors of the wrong order
c  so the first derivatives are recomputed and the order is set
c  to 1.
c
  560 d = 0.0e+00
      do 570 i = 1,n
  570   d = d + (error(i)/ymax(i))**2
      iweval = 0
      if (d.gt.e) go to 610
      if (k.lt.3) go to 590
c
c  complete the correction of the higher order derivatives after a
c  successful step.
c
      do 580 j = 3,k
        do 580 i = 1,n
  580     ytable(i,j) = ytable(i,j) + a(j)*error(i)
  590 kflag = +1
      hnew = h
      if (idoub.le.1) go to 620
      idoub = idoub - 1
      if (idoub.gt.1) go to 770
      do 600 i = 1,n
  600   save10(i) = error(i)
      go to 770
c
c  reduce the failure flag count to check for multiple failures.
c  restore t to its original value and try again unless there have
c  three failures.  in that case the derivatives are assumed to have
c  accumulated errors so a restart from the current values of y is
c  tried.
c
  610 kflag = kflag - 2
      if (h.le.(hmin*1.00001e+00)) go to 810
      t = told
      if (kflag.le.-5) go to 790
c
c  pr1, pr2, and pr3  will contain the amounts by which the step size
c  should be divided at order one lower, at this order, and at order
c  one higher respectively.
c
  620 pr2 = (d/e)**enq2*1.2e+00
      pr3 = 1.0e+20
      if ((nq.ge.maxder).or.(kflag.le.-1)) go to 640
      d = 0.0e+00
      do 630 i = 1,n
  630   d = d + ((error(i) - save10(i))/ymax(i))**2
      pr3 = (d/eup)**enq3*1.4e+00
  640 pr1 = 1.0e+20
      if (nq.le.1) go to 660
      d = 0.0e+00
      do 650 i = 1,n
  650   d = d + (ytable(i,k)/ymax(i))**2
      pr1 = (d/edwn)**enq1*1.3e+00
  660 continue
      if (pr2.le.pr3) go to 720
      if (pr3.lt.pr1) go to 730
  670 r = 1.0e+00/max(pr1,0.0001e+00)
      newq = nq -1
  680 idoub = 10
      if ((kflag.eq.1).and.(r.lt.(1.1e+00))) go to 770
      if (newq.le.nq) go to 700
c
c  compute one additional scaled derivative if order is increased.
c
      do 690 i = 1,n
  690   ytable(i,newq+1) = error(i)*a(k)/real(k)
  700 k = newq + 1
      if (kflag.eq.1) go to 740
      racum = racum*r
      iret1 = 3
      go to 820
  710 if (newq.eq.nq) go to 330
      nq = newq
      go to 150
  720 if (pr2.gt.pr1) go to 670
      newq = nq
      r = 1.0e+00/max(pr2,0.0001e+00)
      go to 680
  730 r = 1.0e+00/max(pr3,0.0001e+00)
      newq = nq + 1
      go to 680
  740 iret = 2
      h = h*r
      hnew = h
      if (nq.eq.newq) go to 750
      nq = newq
      go to 150
  750 r1= 1.0e+00
      do 760 j = 2,k
        r1 = r1*r
        do 760 i = 1,n
  760     ytable(i,j) = ytable(i,j)*r1
      idoub = k
  770 do 780 i = 1,n
  780   ymax(i) = max(ymax(i),abs(ytable(i,1)))
      jstart = nq
      if(h.gt.0.0e+00.and.t.ge.tout)go to 860
      if(h.lt.0.0e+00.and.t.le.tout)go to 860
      go to 70
  790 if (nq.eq.1) go to 850
      call yprime(t,ytable,save11,n)
      r = h/hold
      do 800 i = 1,n
        ytable(1,i) = save(i,1)
        save(i,2) = hold*save11(i)
  800   ytable(i,2) = save(i,2)*r
      nq = 1
      kflag = 1
      go to 150
  810 kflag = -1
      hnew = h
      jstart = nq
      go to 860
c
c  this section scales all variables connected with h and returns
c  to the entering section.
c
  820 racum = max(abs(hmin/hold),racum)
      r1 = 1.0e+00
      do 830 j = 2,k
        r1 = r1*racum
        do 830 i = 1,n
  830     ytable(i,j) = save(i,j)*r1
      h = hold*racum
      do 840 i = 1,n
  840   ytable(i,1) = save(i,1)
      idoub = k
      go to (110,330,710),iret1
  850 kflag = -4
      go to 540
c
c  copy y from ytable before return
c
  860 continue
      do 870 i=1,n
        y(i)=ytable(i,1)
  870   continue
      return
 1000 format(' fatal error - n.le.0 or n.gt.20')
 1010 format(' this program is forcing a stop')
 1020 format(' fatal error - (tout-t).lt.hmin')
 1030 format(' fatal error - hmin.le.0.0')
 1040 format(' fatal error - eps.le.0.0')
 1050 format(' fatal error - mf.lt.0 or mf.gt.2')
      end
      subroutine pdummy ( t, y, dermat, n )

c*********************************************************************72
c
cc PDUMMY is a dummy derivative routine for calls to SMVODE with mf.ne.1.
c
      integer n

      real dermat(n,n)
      real t
      real y(n)

      write(6,1000)
      write(6,1010)
      write(6,1020)
      write(6,1030)
      write(6,1040)
      write(6,1050)
      write(6,1060)
      stop
 1000 format(1x)
 1010 format(' error message from subroutine pdummy')
 1020 format(' apparently, although pdummy was specified as')
 1030 format(' the dummy jacobian routine, which means that')
 1040 format(' the mf option must be 0 or 2, the mf option')
 1050 format(' was set to 1.  this is illegal.  fix it.')
 1060 format(' the program will stop now.')
      end
      subroutine sgefa(a,lda,n,ipvt,info)

c*********************************************************************72
c
cc SGEFA factors a matrix by gaussian elimination.
c
c     sgefa is usually called by sgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for sgeco) = (1 + 9/n)*(time for sgefa) .
c
c     on entry
c
c        a       dimension(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that sgesl or sgedi will divide by zero
c                     if called.  use  rcond  in sgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas saxpy,sscal,isamax
c
c     internal variables
c
      integer lda,n,ipvt(1),info
      dimension a(lda,1)
      integer isamax,j,k,kp1,l,nm1
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = isamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0e0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0e0/a(k,k)
            call sscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call saxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0e0) info = n
      return
      end
      subroutine sgesl(a,lda,n,ipvt,b,job)

c*********************************************************************72
c
cc SGESL solves a linear system that was factored by SGEFA.
c
c     sgesl solves the system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by sgeco or sgefa.
c
c     on entry
c
c        a       dimension(lda, n)
c                the output from sgeco or sgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from sgeco or sgefa.
c
c        b       dimension(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if sgeco has set rcond .gt. 0.0
c        or sgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call sgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call sgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas saxpy,sdot
c
      integer lda,n,ipvt(1),job
      dimension a(lda,1),b(1)
      integer k,kb,l,nm1
c
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call saxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call saxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            t = sdot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + sdot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end
      function isamax ( n, sx, incx )

c*********************************************************************72
c
cc ISAMAX finds the index of element having maximum absolute value.
c
c  Discussion:
c
c    This routine uses single precision real arithmetic.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Jack Dongarra
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, real X(*), the vector to be examined.
c
c    Input, integer INCX, the increment between successive entries of SX.
c
c    Output, integer ISAMAX, the index of the element of SX of maximum
c    absolute value.
c
      implicit none

      integer i
      integer incx
      integer isamax
      integer ix
      integer n
      real sx(*)
      real smax

      isamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      isamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c  code for increment not equal to 1
c
      ix = 1
      smax = abs(sx(1))
      ix = ix + incx
      do i = 2,n
         if(abs(sx(ix)).le.smax) go to 5
         isamax = i
         smax = abs(sx(ix))
    5    ix = ix + incx
      end do

      return
c
c  code for increment equal to 1
c
   20 smax = abs(sx(1))
      do i = 2,n
        if( smax .lt. abs(sx(i)) ) then
          isamax = i
          smax = abs(sx(i))
        end if
      end do

      return
      end
      subroutine saxpy ( n, sa, sx, incx, sy, incy )

c*********************************************************************72
c
cc SAXPY computes constant times a vector plus a vector.
c
c  Discussion:
c
c    This routine uses single precision real arithmetic.
c
c    This routine uses unrolled loop for increments equal to one.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Jack Dongarra
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, real SA, the multiplier.
c
c    Input, real X(*), the vector to be scaled and added to Y.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Input/output, real Y(*), the vector to which a multiple of X is to
c    be added.
c
c    Input, integer INCY, the increment between successive entries of Y.
c
      implicit none

      real sx(*),sy(*),sa
      integer i,incx,incy,ix,iy,m,n

      if(n.le.0)return
      if (sa .eq. 0.0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        sy(iy) = sy(iy) + sa*sx(ix)
        ix = ix + incx
        iy = iy + incy
      end do
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        sy(i) = sy(i) + sa*sx(i)
      end do
      if( n .lt. 4 ) return
   40 continue

      do i = m+1, n, 4
        sy(i) = sy(i) + sa*sx(i)
        sy(i + 1) = sy(i + 1) + sa*sx(i + 1)
        sy(i + 2) = sy(i + 2) + sa*sx(i + 2)
        sy(i + 3) = sy(i + 3) + sa*sx(i + 3)
      end do

      return
      end
      function sdot ( n, sx, incx, sy, incy )

c*********************************************************************72
c
cc SDOT forms the dot product of two vectors.
c
c  Discussion:
c
c    This routine uses single precision real arithmetic.
c
c    This routine uses unrolled loops for increments equal to one.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Jack Dongarra
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vectors.
c
c    Input, real X(*), one of the vectors to be multiplied.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Input, real Y(*), one of the vectors to be multiplied.
c
c    Input, integer INCY, the increment between successive elements of Y.
c
c    Output, real SDOT, the dot product of X and Y.
c
      implicit none

      real sdot
      real sx(*),sy(*),stemp
      integer i,incx,incy,ix,iy,m,n

      stemp = 0.0e0
      sdot = 0.0e0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        stemp = stemp + sx(ix)*sy(iy)
        ix = ix + incx
        iy = iy + incy
      end do
      sdot = stemp
      return
c
c  code for both increments equal to 1
c
c
c  clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        stemp = stemp + sx(i)*sy(i)
      end do

      if( n .lt. 5 ) go to 60
   40 continue

      do i = m+1, n, 5
        stemp = stemp + sx(i)*sy(i) + sx(i + 1)*sy(i + 1) +
     &   sx(i + 2)*sy(i + 2) + sx(i + 3)*sy(i + 3) + sx(i + 4)*sy(i + 4)
      end do

   60 sdot = stemp

      return
      end
      subroutine sscal ( n, sa, sx, incx )

c*********************************************************************72
c
cc SSCAL scales a vector by a constant.
c
c  Discussion:
c
c    This routine uses single precision real arithmetic.
c
c    This routine uses unrolled loops for increment equal to 1.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Jack Dongarra
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, real SA, the multiplier.
c
c    Input/output, real X(*), the vector to be scaled.
c
c    Input, integer INCX, the increment between successive entries of X.
c
      implicit none

      real sa,sx(*)
      integer i,incx,m,n,nincx

      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do i = 1,nincx,incx
        sx(i) = sa*sx(i)
      end do
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        sx(i) = sa*sx(i)
      end do
      if( n .lt. 5 ) return
   40 continue

      do i = m+1, n, 5
        sx(i) = sa*sx(i)
        sx(i + 1) = sa*sx(i + 1)
        sx(i + 2) = sa*sx(i + 2)
        sx(i + 3) = sa*sx(i + 3)
        sx(i + 4) = sa*sx(i + 4)
      end do

      return
      end
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
c    12 January 2007
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

      character * ( 8 ) ampm
      integer d
      character * ( 8 ) date
      integer h
      integer m
      integer mm
      character * ( 9 ) month(12)
      integer n
      integer s
      character * ( 10 ) time
      integer y

      save month

      data month /
     &  'January  ', 'February ', 'March    ', 'April    ',
     &  'May      ', 'June     ', 'July     ', 'August   ',
     &  'September', 'October  ', 'November ', 'December ' /

      call date_and_time ( date, time )

      read ( date, '(i4,i2,i2)' ) y, m, d
      read ( time, '(i2,i2,i2,1x,i3)' ) h, n, s, mm

      if ( h .lt. 12 ) then
        ampm = 'AM'
      else if ( h .eq. 12 ) then
        if ( n .eq. 0 .and. s .eq. 0 ) then
          ampm = 'Noon'
        else
          ampm = 'PM'
        end if
      else
        h = h - 12
        if ( h .lt. 12 ) then
          ampm = 'PM'
        else if ( h .eq. 12 ) then
          if ( n .eq. 0 .and. s .eq. 0 ) then
            ampm = 'Midnight'
          else
            ampm = 'AM'
          end if
        end if
      end if

      write ( *,
     &  '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' )
     &  d, month(m), y, h, ':', n, ':', s, '.', mm, ampm

      return
      end
