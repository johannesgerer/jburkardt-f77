      subroutine cgbd ( iprint, iter, nfun, gnorm, n, x, f, g, stp,
     &  finish, ndes, im, betafr, betapr, beta )

c*********************************************************************72
C
cc CGBD prints monitoring information.
c
c  Discussion:
c
c    The frequency and amount of output are controlled by IPRINT.
c
c  Modified:
c
c    18 December 2008
c
c  Author:
c
c    Jean Charles Gilbert, Jorge Nocedal.
c
c  Reference:
c
c    Jean Charles Gilbert, Jorge Nocedal,
c    Global Convergence Properties of Conjugate Gradient Methods,
c    SIAM Journal on Optimization,
c    Volume 2, Number 1, 1992, pages 21-42. 
c
c  Parameters
c
      implicit none

      double precision beta
      double precision betafr
      double precision betapr
      double precision f
      logical finish
      double precision g(n)
      double precision gnorm
      integer i
      integer im
      integer iprint(2)
      integer iter
      integer n
      integer ndes
      integer nfun
      double precision stp
      double precision x(n)

      if (iter.eq.0)then
           write(*,'(a)' ) ' '
           write(*,10)
           write(*,20) n
           write(*,30) f,gnorm
                 if (iprint(2).ge.1)then
                     write(*,40)
                     write(*,50) (x(i),i=1,n)
                     write(*,60)
                     write(*,50) (g(i),i=1,n)
                 endif
           write(*,10)
           write(*,70)
      else
          if ((iprint(1).eq.0).and.(iter.ne.1.and..not.finish))return
          if (iprint(1).ne.0)then
               if(mod(iter-1,iprint(1)).eq.0.or.finish)then
                     if(iprint(2).gt.1.and.iter.gt.1) write(*,70)
                     write(*,80)iter,nfun,f,gnorm,stp,beta
               else
                     return
               endif
          else
               if( iprint(2).gt.1.and.finish) write(*,70)
               write(*,80)iter,nfun,f,gnorm,stp,beta
          endif
          if (iprint(2).eq.2.or.iprint(2).eq.3)then
                  write(*,40)
                  write(*,50)(x(i),i=1,n)
              if (iprint(2).eq.3)then
                  write(*,60)
                  write(*,50)(g(i),i=1,n)
              endif
          endif
          if (finish) write(*,100)
      endif

 10   FORMAT('*************************************************')
 20   FORMAT(' N=',I5,//,'INITIAL VALUES:')
 30   FORMAT(' F= ',1PD10.3,'   GNORM= ',1PD10.3)
 40   FORMAT(/,' VECTOR X= ')
 50   FORMAT(6(2X,1PD10.3/))
 60   FORMAT(' GRADIENT VECTOR G= ')
 70   FORMAT(/'   I  NFN',4X,'FUNC',7X,'GNORM',6X,
     *   'STEPLEN',4x,'BETA',/,
     *   ' ----------------------------------------------------')
 80   FORMAT(I4,1X,I3,2X,2(1PD10.3,2X),1PD8.1,2x,1PD8.1)
100   FORMAT(/' SUCCESSFUL CONVERGENCE (NO ERRORS).'
     *          ,/,' IFLAG = 0')
c
      return
      end
      subroutine cgfam ( n, x, f, g, d, gold, iprint, eps, w, iflag, 
     &  irest, method, finish )

c*********************************************************************72
c
cc CGFAM implements conjugate gradient methods for unconstrained nonlinear optimization.
c
c  Modified:
c
c    18 December 2008
c
c  Author:
c
c    Jean Charles Gilbert, Jorge Nocedal.
c
c  Reference:
c
c    Jean Charles Gilbert, Jorge Nocedal,
c    Global Convergence Properties of Conjugate Gradient Methods,
c    SIAM Journal on Optimization,
c    Volume 2, Number 1, 1992, pages 21-42. 
C
c  Parameters:
c
c     n      =  number of variables
c     x      =  iterate
c     f      =  function value
c     g      =  gradient value
c     gold   =  previous gradient value
c
c    input, integer iprint(2), controls printing.
c    iprint(1) < 0 : no output is generated
c    iprint(1) = 0 : output only at first and last iteration
c    iprint(1) > 0 : output every iprint(1) iterations
c    iprint(2)     : specifies the type of output generated;
c                    the larger the value (between 0 and 3),
c                    the more information
c    iprint(2) = 0 : no additional information printed
c    iprint(2) = 1 : initial x and gradient vectors printed
c    iprint(2) = 2 : x vector printed every iteration
c    iprint(2) = 3 : x vector and gradient vector printed 
c                       every iteration 
c     eps    =  convergence constant
c     w      =  working array
c     iflag  =  controls termination of code, and return to main
c               program to evaluate function and gradient
c               iflag = -3 : improper input parameters
c               iflag = -2 : descent was not obtained
c               iflag = -1 : line search failure
c               iflag =  0 : initial entry or 
c                            successful termination without error   
c               iflag =  1 : indicates a re-entry with new function values
c               iflag =  2 : indicates a re-entry with a new iterate
c     irest  =  0 (no restarts); 1 (restart every n steps)
c     method =  1 : fletcher-reeves 
c               2 : polak-ribiere
c               3 : positive polak-ribiere ( beta=max{beta,0} )
c
      implicit none

      integer n

      double precision eps
      double precision f
      double precision g(n)

      double precision x(n),d(n),gold(n),w(n)
      integer iprint(2),iflag,irest,method,im,ndes
      double precision gtol,one,zero,gnorm,ddot,stp1,ftol,xtol,stpmin,
     .       stpmax,stp,beta,betafr,betapr,dg0,gg,gg0,dgold,
     .       dgout,dg,dg1
      integer iter,nfun,maxfev,info,i,nfev,nrst,ides
      logical new,finish
c
c     iter: keeps track of the number of iterations
c     nfun: keeps track of the number of function/gradient evaluations
      common /runinf/iter,nfun
      save
      data one,zero/1.0d+0,0.0d+0/
c
c iflag = 1 indicates a re-entry with new function values
      if(iflag.eq.1) go to 72
c
c iflag = 2 indicates a re-entry with a new iterate
      if(iflag.eq.2) go to 80
c
c     initialize
c
c     im =   number of times betapr was negative for method 2 or
c            number of times betapr was 0 for method 3
c
c     ndes = number of line search iterations after wolfe conditions
c            were satisfied
c
      iter= 0
      if(n.le.0) go to 96
      nfun= 1
      new=.true.
      nrst= 0
      im=0
      ndes=0

      do i = 1, n
        d(i) = - g(i)
      end do

      gnorm= dsqrt(ddot(n,g,1,g,1))
      stp1= one/gnorm
c
c     parameters for line search routine
c
c     ftol and gtol are nonnegative input variables. termination
c       occurs when the sufficient decrease condition and the
c       directional derivative condition are satisfied.
c
c     xtol is a nonnegative input variable. termination occurs
c       when the relative width of the interval of uncertainty
c       is at most xtol.
c
c     stpmin and stpmax are nonnegative input variables which
c       specify lower and upper bounds for the step.
c
c     maxfev is a positive integer input variable. termination
c       occurs when the number of calls to fcn is at least
c       maxfev by the end of an iteration.

      ftol= 1.0d-4
      gtol= 1.0d-1
      if(gtol.le.1.d-04) then
        write(*,145)
        gtol=1.d-02
      endif
      xtol= 1.0d-17
      stpmin= 1.0d-20
      stpmax= 1.0d+20
      maxfev= 40

      if(iprint(1).ge.0) call cgbd(iprint,iter,nfun,
     *   gnorm,n,x,f,g,stp,finish,ndes,im,betafr,betapr,beta)
c
c     main iteration loop
c
 8    iter= iter+1
c
c     when nrst>n and irest=1 then restart
c
      nrst= nrst+1
      info=0
c
c
c     call the line search routine of mor'e and thuente
c     (modified for our cg method)
c
c    Jorge More, David Thuente, 
c    Linesearch Algorithms with Guaranteed Sufficient Decrease,
c    ACM Transactions on Mathematical Software,
c    Volume 20, Number 3, September 1994, pages 286-307.
c
      nfev=0

      do i = 1, n
        gold(i) = g(i)
      end do

      dg= ddot(n,d,1,g,1)
      dgold=dg
      stp=one
c
c shanno-phua's formula for trial step
c
      if(.not.new) stp= dg0/dg
      if (iter.eq.1) stp=stp1
      ides=0
      new=.false.
  72  continue
c
c     write(6,*) 'step= ', stp
c
c call to the line search subroutine
c
      call cvsmod(n,x,f,g,d,stp,ftol,gtol,
     *            xtol,stpmin,stpmax,maxfev,info,nfev,w,dg,dgout)

c       info is an integer output variable set as follows:
c         info = 0  improper input parameters.
c         info =-1  a return is made to compute the function and gradient.
c         info = 1  the sufficient decrease condition and the
c                   directional derivative condition hold.
c         info = 2  relative width of the interval of uncertainty
c                   is at most xtol.
c         info = 3  number of calls to fcn has reached maxfev.
c         info = 4  the step is at the lower bound stpmin.
c         info = 5  the step is at the upper bound stpmax.
c         info = 6  rounding errors prevent further progress.
c                   there may not be a step which satisfies the
c                   sufficient decrease and curvature conditions.
c                   tolerances may be too small.

      if (info .eq. -1) then
c       return to fetch function and gradient
        iflag=1
        return
      endif
      if (info .ne. 1) go to 90
c
c     test if descent direction is obtained for methods 2 and 3
c
      gg= ddot(n,g,1,g,1)
      gg0= ddot(n,g,1,gold,1)
      betapr= (gg-gg0)/gnorm**2

      if (irest.eq.1.and.nrst.gt.n) then
        nrst=0
        new=.true.
        go to 75
      endif 
c
      if (method.eq.1) then
        go to 75
      else
        dg1=-gg + betapr*dgout
        if (dg1.lt. 0.0d0 ) go to 75
        if (iprint(1).ge.0) write(6,*) 'no descent'
        ides= ides + 1
        if(ides.gt.5) go to 95
        go to 72
      endif
c
c     determine correct beta value for method chosen
c
c     im =   number of times betapr was negative for method 2 or
c            number of times betapr was 0 for method 3
c
c     ndes = number of line search iterations after wolfe conditions
c            were satisfied
c
  75  nfun= nfun + nfev
      ndes= ndes + ides
      betafr= gg/gnorm**2

      if (nrst.eq.0) then
        beta= zero
      else
        if (method.eq.1) beta=betafr
        if (method.eq.2) beta=betapr
        if ((method.eq.2.or.method.eq.3).and.betapr.lt.0) im=im+1
        if (method.eq.3) beta=max(zero,betapr)
      endif
c
c  Compute the new direction.
c
      do i = 1, n
        d(i) = - g(i) + beta * d(i)
      end do

      dg0= dgold*stp
c
c     return to driver for termination test
c
      gnorm=dsqrt(ddot(n,g,1,g,1))
      iflag=2
      return

  80  continue
c
c call subroutine for printing output
c
      if(iprint(1).ge.0) call cgbd(iprint,iter,nfun,
     *     gnorm,n,x,f,g,stp,finish,ndes,im,betafr,betapr,beta)
      if (finish) then
         iflag = 0
         return
      end if
      go to 8
c
c     end of main iteration loop. error exits.
c
  90  iflag=-1
      write(*,100) info
      return
  95  iflag=-2
      write(*,135) i
      return
  96  iflag= -3
      write(*,140)
c
c     formats
c     -------
c
 100  format(/' iflag= -1 ',/' line search failed. see'
     .          ' documentation of routine cvsmod',/' error return'
     .          ' of line search: info= ',i2,/
     .          ' possible cause: function or gradient are incorrect')
 135  format(/' iflag= -2',/' descent was not obtained')
 140  format(/' iflag= -3',/' improper input parameters (n',
     .       ' is not positive)')
 145  format(/'  gtol is less than or equal to 1.d-04',
     .       / ' it has been reset to 1.d-02')
      return
      end
      subroutine cstepm ( stx, fx, dx, sty, fy, dy, stp, fp, dp,
     &  brackt, stpmin, stpmax, info )

c*********************************************************************72
c
cc CSTEPM computes a safeguarded step for a line search.
c
c  Discussion:
c
c    The routine computes a safeguarded step for a line search, and 
c    updates an interval of uncertainty for a minimizer of the function.
c
c     the parameter stx contains the step with the least function
c     value. the parameter stp contains the current step. it is
c     assumed that the derivative at stx is negative in the
c     direction of the step. if brackt is set true then a
c     minimizer has been bracketed in an interval of uncertainty
c     with endpoints stx and sty.
c
c     the subroutine statement is
c
c       subroutine cstepm(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,
c                        stpmin,stpmax,info)
c
c     where
c
c       stx, fx, and dx are variables which specify the step,
c         the function, and the derivative at the best step obtained
c         so far. the derivative must be negative in the direction
c         of the step, that is, dx and stp-stx must have opposite
c         signs. on output these parameters are updated appropriately.
c
c       sty, fy, and dy are variables which specify the step,
c         the function, and the derivative at the other endpoint of
c         the interval of uncertainty. on output these parameters are
c         updated appropriately.
c
c       stp, fp, and dp are variables which specify the step,
c         the function, and the derivative at the current step.
c         if brackt is set true then on input stp must be
c         between stx and sty. on output stp is set to the new step.
c
c       brackt is a logical variable which specifies if a minimizer
c         has been bracketed. if the minimizer has not been bracketed
c         then on input brackt must be set false. if the minimizer
c         is bracketed then on output brackt is set true.
c
c       stpmin and stpmax are input variables which specify lower
c         and upper bounds for the step.
c
c       info is an integer output variable set as follows:
c         if info = 1,2,3,4,5, then the step has been computed
c         according to one of the five cases below. otherwise
c         info = 0, and this indicates improper input parameters.
c
c
c     argonne national laboratory. minpack project. june 1983
c     jorge j. more', david j. thuente
c
      implicit none

      integer info
      double precision stx,fx,dx,sty,fy,dy,stp,fp,dp,stpmin,stpmax
      logical brackt,bound
      double precision gamma,p,q,r,s,sgnd,stpc,stpf,stpq,theta
      info = 0
c
c     check the input parameters for errors.
c
      if ((brackt .and. (stp .le. min(stx,sty) .or.
     *    stp .ge. max(stx,sty))) .or.
     *    dx*(stp-stx) .ge. 0.0 .or. stpmax .lt. stpmin) return
c
c     determine if the derivatives have opposite sign.
c
      sgnd = dp*(dx/abs(dx))
c
c     first case. a higher function value.
c     the minimum is bracketed. if the cubic step is closer
c     to stx than the quadratic step, the cubic step is taken,
c     else the average of the cubic and quadratic steps is taken.
c
      if (fp .gt. fx) then
         info = 1
         bound = .true.
         theta = 3*(fx - fp)/(stp - stx) + dx + dp
         s = max(abs(theta),abs(dx),abs(dp))
         gamma = s*sqrt((theta/s)**2 - (dx/s)*(dp/s))
         if (stp .lt. stx) gamma = -gamma
         p = (gamma - dx) + theta
         q = ((gamma - dx) + gamma) + dp
         r = p/q
         stpc = stx + r*(stp - stx)
         stpq = stx + ((dx/((fx-fp)/(stp-stx)+dx))/2)*(stp - stx)
         if (abs(stpc-stx) .lt. abs(stpq-stx)) then
            stpf = stpc
         else
           stpf = stpc + (stpq - stpc)/2
           end if
         brackt = .true.
c
c     second case. a lower function value and derivatives of
c     opposite sign. the minimum is bracketed. if the cubic
c     step is closer to stx than the quadratic (secant) step,
c     the cubic step is taken, else the quadratic step is taken.
c
      else if (sgnd .lt. 0.0) then
         info = 2
         bound = .false.
         theta = 3*(fx - fp)/(stp - stx) + dx + dp
         s = max(abs(theta),abs(dx),abs(dp))
         gamma = s*sqrt((theta/s)**2 - (dx/s)*(dp/s))
         if (stp .gt. stx) gamma = -gamma
         p = (gamma - dp) + theta
         q = ((gamma - dp) + gamma) + dx
         r = p/q
         stpc = stp + r*(stx - stp)
         stpq = stp + (dp/(dp-dx))*(stx - stp)
         if (abs(stpc-stp) .gt. abs(stpq-stp)) then
            stpf = stpc
         else
            stpf = stpq
            end if
         brackt = .true.
c
c     third case. a lower function value, derivatives of the
c     same sign, and the magnitude of the derivative decreases.
c     the cubic step is only used if the cubic tends to infinity
c     in the direction of the step or if the minimum of the cubic
c     is beyond stp. otherwise the cubic step is defined to be
c     either stpmin or stpmax. the quadratic (secant) step is also
c     computed and if the minimum is bracketed then the the step
c     closest to stx is taken, else the step farthest away is taken.
c
      else if (abs(dp) .lt. abs(dx)) then
         info = 3
         bound = .true.
         theta = 3*(fx - fp)/(stp - stx) + dx + dp
         s = max(abs(theta),abs(dx),abs(dp))
c
c        the case gamma = 0 only arises if the cubic does not tend
c        to infinity in the direction of the step.
c
         gamma = s*sqrt(max(0.0d0,(theta/s)**2 - (dx/s)*(dp/s)))
         if (stp .gt. stx) gamma = -gamma
         p = (gamma - dp) + theta
         q = (gamma + (dx - dp)) + gamma
         r = p/q
         if (r .lt. 0.0 .and. gamma .ne. 0.0) then
            stpc = stp + r*(stx - stp)
         else if (stp .gt. stx) then
            stpc = stpmax
         else
            stpc = stpmin
            end if
         stpq = stp + (dp/(dp-dx))*(stx - stp)
         if (brackt) then
            if (abs(stp-stpc) .lt. abs(stp-stpq)) then
               stpf = stpc
            else
               stpf = stpq
               end if
         else
            if (abs(stp-stpc) .gt. abs(stp-stpq)) then
               stpf = stpc
            else
               stpf = stpq
               end if
            end if
c
c     fourth case. a lower function value, derivatives of the
c     same sign, and the magnitude of the derivative does
c     not decrease. if the minimum is not bracketed, the step
c     is either stpmin or stpmax, else the cubic step is taken.
c
      else
         info = 4
         bound = .false.
         if (brackt) then
            theta = 3*(fp - fy)/(sty - stp) + dy + dp
            s = max(abs(theta),abs(dy),abs(dp))
            gamma = s*sqrt((theta/s)**2 - (dy/s)*(dp/s))
            if (stp .gt. sty) gamma = -gamma
            p = (gamma - dp) + theta
            q = ((gamma - dp) + gamma) + dy
            r = p/q
            stpc = stp + r*(sty - stp)
            stpf = stpc
         else if (stp .gt. stx) then
            stpf = stpmax
         else
            stpf = stpmin
            end if
         end if
c
c     update the interval of uncertainty. this update does not
c     depend on the new step or the case analysis above.
c
      if (fp .gt. fx) then
         sty = stp
         fy = fp
         dy = dp
      else
         if (sgnd .lt. 0.0) then
            sty = stx
            fy = fx
            dy = dx
            end if
         stx = stp
         fx = fp
         dx = dp
         end if
c
c     compute the new step and safeguard it.
c
      stpf = min(stpmax,stpf)
      stpf = max(stpmin,stpf)
      stp = stpf
      if (brackt .and. bound) then
         if (sty .gt. stx) then
            stp = min(stx+0.66*(sty-stx),stp)
         else
            stp = max(stx+0.66*(sty-stx),stp)
            end if
         end if
      return
      end
      subroutine cvsmod(n,x,f,g,s,stp,ftol,gtol,xtol,
     *           stpmin,stpmax,maxfev,info,nfev,wa,dginit,dgout)

c*********************************************************************72
c
cc CSVMOD finds a step which satisfies a sufficient decrease condition.
c
c  Discussion:
c
c    The routine finds a step which satisfies a sufficient decrease condition
c    and a curvature condition.
c
c     the user must provide a subroutine which calculates the
c     function and the gradient.
c
c     at each stage the subroutine updates an interval of
c     uncertainty with endpoints stx and sty. the interval of
c     uncertainty is initially chosen so that it contains a
c     minimizer of the modified function
c
c          f(x+stp*s) - f(x) - ftol*stp*(gradf(x)'s).
c
c     if a step is obtained for which the modified function
c     has a nonpositive function value and nonnegative derivative,
c     then the interval of uncertainty is chosen so that it
c     contains a minimizer of f(x+stp*s).
c
c     the algorithm is designed to find a step which satisfies
c     the sufficient decrease condition
c
c           f(x+stp*s) .le. f(x) + ftol*stp*(gradf(x)'s),
c
c     and the curvature condition
c
c           abs(gradf(x+stp*s)'s)) .le. gtol*abs(gradf(x)'s).
c
c     if ftol is less than gtol and if, for example, the function
c     is bounded below, then there is always a step which satisfies
c     both conditions. if no step can be found which satisfies both
c     conditions, then the algorithm usually stops when rounding
c     errors prevent further progress. in this case stp only
c     satisfies the sufficient decrease condition.
c
c     the subroutine statement is
c
c        subroutine cvsmod(n,x,f,g,s,stp,ftol,gtol,xtol,
c                   stpmin,stpmax,maxfev,info,nfev,wa,dg,dgout)
c     where
c
c       n is a positive integer input variable set to the number
c         of variables.
c
c       x is an array of length n. on input it must contain the
c         base point for the line search. on output it contains
c         x + stp*s.
c
c       f is a variable. on input it must contain the value of f
c         at x. on output it contains the value of f at x + stp*s.
c
c       g is an array of length n. on input it must contain the
c         gradient of f at x. on output it contains the gradient
c         of f at x + stp*s.
c
c       s is an input array of length n which specifies the
c         search direction.
c
c       stp is a nonnegative variable. on input stp contains an
c         initial estimate of a satisfactory step. on output
c         stp contains the final estimate.
c
c       ftol and gtol are nonnegative input variables. termination
c         occurs when the sufficient decrease condition and the
c         directional derivative condition are satisfied.
c
c       xtol is a nonnegative input variable. termination occurs
c         when the relative width of the interval of uncertainty
c         is at most xtol.
c
c       stpmin and stpmax are nonnegative input variables which
c         specify lower and upper bounds for the step.
c
c       maxfev is a positive integer input variable. termination
c         occurs when the number of calls to fcn is at least
c         maxfev by the end of an iteration.
c
c       info is an integer output variable set as follows:
c
c         info = 0  improper input parameters.
c
c         info =-1  a return is made to compute the function and gradient.
c
c         info = 1  the sufficient decrease condition and the
c                   directional derivative condition hold.
c
c         info = 2  relative width of the interval of uncertainty
c                   is at most xtol.
c
c         info = 3  number of calls to fcn has reached maxfev.
c
c         info = 4  the step is at the lower bound stpmin.
c
c         info = 5  the step is at the upper bound stpmax.
c
c         info = 6  rounding errors prevent further progress.
c                   there may not be a step which satisfies the
c                   sufficient decrease and curvature conditions.
c                   tolerances may be too small.
c
c       nfev is an integer output variable set to the number of
c         calls to fcn.
c
c       wa is a work array of length n.
c
c       *** the following two parameters are a modification to the code
c
c       dg is the initial directional derivative (in the original code
c                 it was computed in this routine0
c
c       dgout is the value of the directional derivative when the wolfe
c             conditions hold, and an exit is made to check descent.
c
c     subprograms called
c
c       cstepm
c
c       fortran-supplied...abs,max,min
c
c     argonne national laboratory. minpack project. june 1983
c     jorge j. more', david j. thuente
c
      implicit none

      integer n,maxfev,info,nfev
      double precision f,stp,ftol,gtol,xtol,stpmin,stpmax
      double precision x(n),g(n),s(n),wa(n)

      save

      integer infoc,j
      logical brackt,stage1
      double precision dg,dgm,dginit,dgtest,dgx,dgxm,dgy,dgym,
     *       finit,ftest1,fm,fx,fxm,fy,fym,p5,p66,stx,sty,
     *       stmin,stmax,width,width1,xtrapf,zero,dgout

      data p5,p66,xtrapf,zero /0.5d0,0.66d0,4.0d0,0.0d0/

      if(info.eq.-1) go to 45
      if(info.eq.1) go to 321
      infoc = 1
c
c     check the input parameters for errors.
c
      if (n .le. 0 .or. stp .le. zero .or. ftol .lt. zero .or.
     *    gtol .lt. zero .or. xtol .lt. zero .or. stpmin .lt. zero
     *    .or. stpmax .lt. stpmin .or. maxfev .le. 0) return
c
c     compute the initial gradient in the search direction
c     and check that s is a descent direction.
c
      if (dginit .ge. zero) return
c
c     initialize local variables.
c
      brackt = .false.
      stage1 = .true.
      nfev = 0
      finit = f
      dgtest = ftol*dginit
      width = stpmax - stpmin
      width1 = width/p5
      do j = 1, n
         wa(j) = x(j)
      end do
c
c     the variables stx, fx, dgx contain the values of the step,
c     function, and directional derivative at the best step.
c     the variables sty, fy, dgy contain the value of the step,
c     function, and derivative at the other endpoint of
c     the interval of uncertainty.
c     the variables stp, f, dg contain the values of the step,
c     function, and derivative at the current step.
c
      stx = zero
      fx = finit
      dgx = dginit
      sty = zero
      fy = finit
      dgy = dginit
c
c     start of iteration.
c
   30 continue
c
c        set the minimum and maximum steps to correspond
c        to the present interval of uncertainty.
c
         if (brackt) then
            stmin = min(stx,sty)
            stmax = max(stx,sty)
         else
            stmin = stx
            stmax = stp + xtrapf*(stp - stx)
         end if
c
c        force the step to be within the bounds stpmax and stpmin.
c
         stp = max(stp,stpmin)
         stp = min(stp,stpmax)
c
c        if an unusual termination is to occur then let
c        stp be the lowest point obtained so far.
c
         if ((brackt .and. (stp .le. stmin .or. stp .ge. stmax))
     *      .or. nfev .ge. maxfev-1 .or. infoc .eq. 0
     *      .or. (brackt .and. stmax-stmin .le. xtol*stmax)) stp = stx
c
c        evaluate the function and gradient at stp
c        and compute the directional derivative.
c
         do j = 1, n
            x(j) = wa(j) + stp*s(j)
         end do
c
c        return to compute function value
c
         info=-1
         return

   45    info=0
         nfev = nfev + 1
         dg = zero
         do j = 1, n
            dg = dg + g(j)*s(j)
         end do
         ftest1 = finit + stp*dgtest
c
c        test for convergence.
c
         if ((brackt .and. (stp .le. stmin .or. stp .ge. stmax))
     *      .or. infoc .eq. 0) info = 6
         if (stp .eq. stpmax .and.
     *       f .le. ftest1 .and. dg .le. dgtest) info = 5
         if (stp .eq. stpmin .and.
     *       (f .gt. ftest1 .or. dg .ge. dgtest)) info = 4
         if (nfev .ge. maxfev) info = 3
         if (brackt .and. stmax-stmin .le. xtol*stmax) info = 2
c        more's code has been modified so that at least one new
c        function value is computed during the line search (enforcing
c        at least one interpolation is not easy, since the code may
c        override an interpolation)
         if (f .le. ftest1 .and. abs(dg) .le. gtol*(-dginit).
     *       and.nfev.gt.1) info = 1
c
c        check for termination.
c
         if (info .ne. 0)then
            dgout=dg
            return
         endif
 321     continue
c
c        in the first stage we seek a step for which the modified
c        function has a nonpositive value and nonnegative derivative.
c
         if (stage1 .and. f .le. ftest1 .and.
     *       dg .ge. min(ftol,gtol)*dginit) stage1 = .false.
c
c  a modified function is used to predict the step only if
c  we have not obtained a step for which the modified
c  function has a nonpositive function value and nonnegative
c  derivative, and if a lower function value has been
c  obtained but the decrease is not sufficient.
c
         if (stage1 .and. f .le. fx .and. f .gt. ftest1) then
c
c  define the modified function and derivative values.
c
            fm = f - stp*dgtest
            fxm = fx - stx*dgtest
            fym = fy - sty*dgtest
            dgm = dg - dgtest
            dgxm = dgx - dgtest
            dgym = dgy - dgtest
c
c  call cstepm to update the interval of uncertainty
c  and to compute the new step.
c
            call cstepm(stx,fxm,dgxm,sty,fym,dgym,stp,fm,dgm,
     *                 brackt,stmin,stmax,infoc)
c
c  reset the function and gradient values for f.
c
            fx = fxm + stx*dgtest
            fy = fym + sty*dgtest
            dgx = dgxm + dgtest
            dgy = dgym + dgtest
         else
c
c  call cstepm to update the interval of uncertainty
c  and to compute the new step.
c
            call cstepm(stx,fx,dgx,sty,fy,dgy,stp,f,dg,
     *                 brackt,stmin,stmax,infoc)
            end if
c
c  force a sufficient decrease in the size of the
c  interval of uncertainty.
c
         if (brackt) then
            if (abs(sty-stx) .ge. p66*width1)
     *         stp = stx + p5*(sty - stx)
            width1 = width
            width = abs(sty-stx)
            end if

         go to 30
      end
      function dasum ( n, dx, incx )

c*********************************************************************72
c
cc DASUM takes the sum of the absolute values.
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
c
c  Modified:
c
c    18 December 2008
c
c  Author:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart
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
c    Input, double precision X(*), the vector to be examined.
c
c    Input, integer INCX, the increment between successive entries of X.
c    INCX must not be negative.
c
c    Output, double precision DASUM, the sum of the absolute values of X.
c
      implicit none

      double precision dasum
      double precision dtemp
      double precision dx(*)
      integer i
      integer incx
      integer m
      integer n
      integer nincx

      dasum = 0.0D+00
      dtemp = 0.0D+00

      if( n .le. 0 ) then
        return
      end if

      if ( incx .le. 0 ) then
        return
      end if

      if ( incx .ne. 1 ) then

        nincx = n * incx
        do i = 1, nincx, incx
          dtemp = dtemp + dabs ( dx(i) )
        end do

      else

        m = mod ( n, 6 )

        do i = 1,m
          dtemp = dtemp + dabs ( dx(i) )
        end do

        do i = m + 1, n, 6
          dtemp = dtemp 
     &      + dabs ( dx(i) ) 
     &      + dabs ( dx(i+1) ) 
     &      + dabs ( dx(i+2) )
     &      + dabs ( dx(i+3) ) 
     &      + dabs ( dx(i+4) ) 
     &      + dabs ( dx(i+5) )
        end do

      end if

      dasum = dtemp

      return
      end
      subroutine daxpy ( n, da, dx, incx, dy, incy )

c*********************************************************************72
c
cc DAXPY computes constant times a vector plus a vector.
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
c
c    This routine uses unrolled loops for increments equal to one.
c
c  Modified:
c
c    18 December 2008
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
c    Input, integer N, the number of elements in DX and DY.
c
c    Input, double precision DA, the multiplier of DX.
c
c    Input, double precision DX(*), the first vector.
c
c    Input, integer INCX, the increment between successive entries of DX.
c
c    Input/output, double precision DY(*), the second vector.
c    On output, DY(*) has been replaced by DY(*) + DA * DX(*).
c
c    Input, integer INCY, the increment between successive entries of DY.
c
      implicit none

      double precision da
      double precision dx(*)
      double precision dy(*)
      integer i
      integer incx
      integer incy
      integer ix
      integer iy
      integer m
      integer n

      if ( n .le. 0 ) then
        return
      end if

      if ( da .eq. 0.0d0 ) then
        return
      end if

      if ( incx .ne. 1 .or. incy .ne. 1 ) then

        ix = 1
        iy = 1
        if ( incx .lt. 0 ) ix = (-n+1)*incx + 1
        if ( incy .lt. 0 ) iy = (-n+1)*incy + 1

        do i = 1, n
          dy(iy) = dy(iy) + da*dx(ix)
          ix = ix + incx
          iy = iy + incy
        end do

      else

        m = mod(n,4)

        do i = 1, m
          dy(i) = dy(i) + da*dx(i)
        end do

        do i = m + 1, n, 4
          dy(i) = dy(i) + da*dx(i)
          dy(i + 1) = dy(i + 1) + da*dx(i + 1)
          dy(i + 2) = dy(i + 2) + da*dx(i + 2)
          dy(i + 3) = dy(i + 3) + da*dx(i + 3)
        end do

      end if

      return
      end
      subroutine dcopy ( n, dx, incx, dy, incy )

c*********************************************************************72
c
cc DCOPY copies a vector.
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
c
c    The routine uses unrolled loops for increments equal to one.
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
c    Input, integer N, the number of elements in DX and DY.
c
c    Input, double precision DX(*), the first vector.
c
c    Input, integer INCX, the increment between successive entries of DX.
c
c    Output, double precision DY(*), the second vector.
c
c    Input, integer INCY, the increment between successive entries of DY.
c
      implicit none

      double precision dx(*)
      double precision dy(*)
      integer i
      integer incx
      integer incy
      integer ix
      integer iy
      integer m
      integer n

      if ( n .le. 0 ) then
        return
      end if

      if(incx.eq.1.and.incy.eq.1)go to 20
c
c  code for unequal increments or equal increments not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
      end do

      return
c
c  code for both increments equal to 1
c
c
c  clean-up loop
c
   20 m = mod(n,7)

      do i = 1,m
        dy(i) = dx(i)
      end do

      do i = m + 1, n, 7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
      end do

      return
      end
      function ddot ( n, dx, incx, dy, incy )

c*********************************************************************72
c
cc DDOT forms the dot product of two vectors.
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
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
c    Input, double precision DX(*), the first vector.
c
c    Input, integer INCX, the increment between successive entries in DX.
c
c    Input, double precision DY(*), the second vector.
c
c    Input, integer INCY, the increment between successive entries in DY.
c
c    Output, double precision DDOT, the sum of the product of the 
c    corresponding entries of DX and DY.
c
      implicit none

      double precision ddot
      double precision dx(*)
      double precision dy(*)
      double precision dtemp
      integer i,incx,incy,ix,iy,m,n

      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c  code for unequal increments or equal increments not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
      end do
      ddot = dtemp
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
        dtemp = dtemp + dx(i)*dy(i)
      end do
      if( n .lt. 5 ) go to 60
   40 continue
      do i = m+1, n, 5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     &   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
      end do

   60 ddot = dtemp

      return
      end
      function dnrm2 ( n, x, incx )

c*********************************************************************72
c
cc DNRM2 returns the euclidean norm of a vector. 
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
c
c      DNRM2 ( X ) = sqrt ( X' * X )
c
c  Author:
c
c    Sven Hammarling
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
c    Input, double precision X(*), the vector whose norm is to be computed.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Output, double precision DNRM2, the Euclidean norm of X.
c
      implicit none

      integer                           incx, n

      double precision dnrm2 
      double precision                  x( * )


      double precision      one         , zero
      parameter           ( one = 1.0d+0, zero = 0.0d+0 )

      integer               ix
      double precision      absxi, norm, scale, ssq

      intrinsic             abs, sqrt

      if( n.lt.1 .or. incx.lt.1 )then
         norm  = zero
      else if( n.eq.1 )then
         norm  = abs( x( 1 ) )
      else
         scale = zero
         ssq   = one
c
c  The following loop is equivalent to this call to the LAPACK
c  auxiliary routine:
c  call dlassq( n, x, incx, scale, ssq )
c
         do ix = 1, 1 + ( n - 1 )*incx, incx
            if( x( ix ).ne.zero )then
               absxi = abs( x( ix ) )
               if( scale.lt.absxi )then
                  ssq   = one   + ssq*( scale/absxi )**2
                  scale = absxi
               else
                  ssq   = ssq   +     ( absxi/scale )**2
               end if
            end if
         end do
         norm  = scale * sqrt( ssq )
      end if

      dnrm2 = norm

      return
      end
      subroutine dscal ( n, da, dx, incx )

c*********************************************************************72
c
cc DSCAL scales a vector by a constant.
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
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
c    Input, double precision SA, the multiplier.
c
c    Input/output, double precision X(*), the vector to be scaled.
c
c    Input, integer INCX, the increment between successive entries of X.
c
      implicit none

      double precision da
      double precision dx(*)
      integer i,incx,m,n,nincx

      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c  code for increment not equal to 1
c
      nincx = n*incx
      do i = 1,nincx,incx
        dx(i) = da*dx(i)
      end do
      return
c
c  code for increment equal to 1
c
c
c  clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        dx(i) = da*dx(i)
      end do
      if( n .lt. 5 ) return
   40 continue
      do i = m+1, n, 5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
      end do

      return
      end
      subroutine dswap ( n, dx, incx, dy, incy )

c*********************************************************************72
c
cc DSWAP interchanges two vectors.
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
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
c    Input/output, double precision X(*), one of the vectors to swap.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Input/output, double precision Y(*), one of the vectors to swap.
c
c    Input, integer INCY, the increment between successive elements of Y.
c
      implicit none

      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,n

      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c  code for unequal increments or equal increments not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
      end do
      return
c
c  code for both increments equal to 1
c
c
c  clean-up loop
c
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
      end do
      if( n .lt. 3 ) return
   40 continue

      do i = m+1, n, 3
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
        dtemp = dx(i + 1)
        dx(i + 1) = dy(i + 1)
        dy(i + 1) = dtemp
        dtemp = dx(i + 2)
        dx(i + 2) = dy(i + 2)
        dy(i + 2) = dtemp
      end do

      return
      end
      function idamax ( n, dx, incx )

c*********************************************************************72
c
cc IDAMAX finds the index of element having maximum absolute value.
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
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
c    Input, double precision X(*), the vector to be examined.
c
c    Input, integer INCX, the increment between successive entries of SX.
c
c    Output, integer IDAMAX, the index of the element of SX of maximum
c    absolute value.
c
      implicit none

      double precision dx(*),dmax
      integer idamax
      integer i,incx,ix,n

      idamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c  code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do  i = 2,n
         if ( dmax .lt. dabs ( dx(ix) ) ) then
           idamax = i
           dmax = dabs(dx(ix))
         end if
         ix = ix + incx
      end do
      return
c
c  code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do i = 2,n
        if( dmax .lt. dabs(dx(i)) ) then
          idamax = i
          dmax = dabs(dx(i))
        end if
      end do

      return
      end
      subroutine timer ( ttime )

c*********************************************************************72
c
cc TIMER returns an estimate for the elapsed user time.
c
c  Discussion:
c
c    TIMER calls ETIME, which may or may not be available as a system
c    library routine.
c
c  Modified:
c
c    18 December 2008
c
c  Parameters:
c
c    Output, double precision TTIME, the elapsed user time in seconds.
c
      implicit none

      real temp
      real tarray(2)
      real etime
      double precision ttime

      temp = etime ( tarray ) 

      ttime = dble ( tarray(1) )
 
      return

      end
