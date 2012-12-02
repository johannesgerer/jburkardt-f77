c-------------------------------------------------------------------------!
c                                                                         !
c        N  A  S     P A R A L L E L     B E N C H M A R K S  3.3         !
c                                                                         !
c                      S E R I A L     V E R S I O N                      !
c                                                                         !
c                                   E P                                   !
c                                                                         !
c-------------------------------------------------------------------------!
c                                                                         !
c    This benchmark is a serial version of the NPB EP code.               !
c    Refer to NAS Technical Reports 95-020 for details.                   !
c                                                                         !
c    Permission to use, copy, distribute and modify this software         !
c    for any purpose with or without fee is hereby granted.  We           !
c    request, however, that all derived work reference the NAS            !
c    Parallel Benchmarks 3.3. This software is provided "as is"           !
c    without express or implied warranty.                                 !
c                                                                         !
c    Information on NPB 3.3, including the technical report, the          !
c    original specifications, source code, results and information        !
c    on how to submit new results, is available at:                       !
c                                                                         !
c           http://www.nas.nasa.gov/Software/NPB/                         !
c                                                                         !
c    Send comments or suggestions to  npb@nas.nasa.gov                    !
c                                                                         !
c          NAS Parallel Benchmarks Group                                  !
c          NASA Ames Research Center                                      !
c          Mail Stop: T27A-1                                              !
c          Moffett Field, CA   94035-1000                                 !
c                                                                         !
c          E-mail:  npb@nas.nasa.gov                                      !
c          Fax:     (650) 604-3957                                        !
c                                                                         !
c-------------------------------------------------------------------------!

      program main

c*********************************************************************72
c
cc MAIN is the main program for EP_SERIAL.
c
c  Discussion:
c
c    This is the serial version of the APP Benchmark 1,
c    the "embarassingly parallel" benchmark.
c
c    M is the Log_2 of the number of complex pairs of uniform (0, 1) random
c    numbers.  MK is the Log_2 of the size of each batch of uniform random
c    numbers.  MK can be set for convenience on a given system, since it does
c    not affect the results.
c
c  Author:
c
c    P. O. Frederickson, D. H. Bailey, A. C. Woo
c
      implicit none

      include 'npbparams.h'

      double precision Mops, epsilon, a, s, t1, t2, t3, t4, x, x1, 
     >                 x2, q, sx, sy, tm, an, tt, gc, dum(3)
      double precision sx_verify_value, sy_verify_value, sx_err, sy_err
      integer          mk, mm, nn, nk, nq, np, 
     >                 i, ik, kk, l, k, nit,
     >                 k_offset, j, fstatus
      logical          verified, timers_enabled
      external         randlc, timer_read
      double precision randlc, timer_read
      character*15     size

      parameter (mk = 16, mm = m - mk, nn = 2 ** mm,
     >           nk = 2 ** mk, nq = 10, epsilon=1.d-8,
     >           a = 1220703125.d0, s = 271828183.d0)

      common/storage/ x(2*nk), q(0:nq-1)
      data             dum /1.d0, 1.d0, 1.d0/


      open(unit=2, file='timer.flag', status='old', iostat=fstatus)
      if (fstatus .eq. 0) then
         timers_enabled = .true.
         close(2)
      else
         timers_enabled = .false.
      endif
c
c   Because the size of the problem is too large to store in a 32-bit
c   integer for some classes, we put it into a string (for printing).
c   Have to strip off the decimal point put in there by the floating
c   point print statement (internal file)
c
      write(*, 1000)
      write(size, '(f15.0)' ) 2.d0**(m+1)
      j = 15
      if (size(j:j) .eq. '.') j = j - 1
      write (*,1001) size(1:j)
      write (*,*)

 1000 format(//,' NAS Parallel Benchmarks (NPB3.3-SER)',
     >          ' - EP Benchmark', /)
 1001 format(' Number of random numbers generated: ', a15)

      verified = .false.
c
c   Compute the number of "batches" of random number pairs generated 
c   per processor. Adjust if the number of processors does not evenly 
c   divide the total number
c
      np = nn 
c
c   Call the random number generator functions and initialize
c   the x-array to reduce the effects of paging on the timings.
c   Also, call all mathematical functions that are used. Make
c   sure these initializations cannot be eliminated as dead code.
c
      call vranlc(0, dum(1), dum(2), dum(3))
      dum(1) = randlc(dum(2), dum(3))

      do i = 1, 2*nk
         x(i) = -1.d99
      end do

      Mops = log(sqrt(abs(max(1.d0,1.d0))))

      
      call timer_clear(1)
      call timer_clear(2)
      call timer_clear(3)
      call timer_start(1)

      t1 = a
      call vranlc(0, t1, a, x)
c
c   Compute AN = A ^ (2 * NK) (mod 2^46).
c
      t1 = a

      do i = 1, mk + 1
         t2 = randlc(t1, t1)
      end do

      an = t1
      tt = s
      gc = 0.d0
      sx = 0.d0
      sy = 0.d0

      do i = 0, nq - 1
         q(i) = 0.d0
      end do
c
c   Each instance of this loop may be performed independently. We compute
c   the k offsets separately to take into account the fact that some nodes
c   have more numbers to generate than others
c
      k_offset = -1

      do k = 1, np

         kk = k_offset + k 
         t1 = s
         t2 = an
c
c  Find starting seed t1 for this kk.
c
         do i = 1, 100
            ik = kk / 2
            if (2 * ik .ne. kk) t3 = randlc(t1, t2)
            if (ik .eq. 0) goto 130
            t3 = randlc(t2, t2)
            kk = ik
         end do
c
c  Compute uniform pseudorandom numbers.
c
 130     continue

         if (timers_enabled) call timer_start(3)
         call vranlc(2 * nk, t1, a, x)
         if (timers_enabled) call timer_stop(3)
c
c  Compute Gaussian deviates by acceptance-rejection method and 
c  tally counts in concentric square annuli.  This loop is not 
c  vectorizable. 
c
         if (timers_enabled) call timer_start(2)

         do i = 1, nk
            x1 = 2.d0 * x(2*i-1) - 1.d0
            x2 = 2.d0 * x(2*i) - 1.d0
            t1 = x1 ** 2 + x2 ** 2
            if (t1 .le. 1.d0) then
               t2   = sqrt(-2.d0 * log(t1) / t1)
               t3   = (x1 * t2)
               t4   = (x2 * t2)
               l    = max(abs(t3), abs(t4))
               q(l) = q(l) + 1.d0
               sx   = sx + t3
               sy   = sy + t4
            endif
         end do

         if (timers_enabled) call timer_stop(2)

      end do

      do i = 0, nq - 1
        gc = gc + q(i)
      end do

      call timer_stop(1)
      tm  = timer_read(1)

      nit=0
      verified = .true.
      if (m.eq.24) then
         sx_verify_value = -3.247834652034740D+3
         sy_verify_value = -6.958407078382297D+3
      elseif (m.eq.25) then
         sx_verify_value = -2.863319731645753D+3
         sy_verify_value = -6.320053679109499D+3
      elseif (m.eq.28) then
         sx_verify_value = -4.295875165629892D+3
         sy_verify_value = -1.580732573678431D+4
      elseif (m.eq.30) then
         sx_verify_value =  4.033815542441498D+4
         sy_verify_value = -2.660669192809235D+4
      elseif (m.eq.32) then
         sx_verify_value =  4.764367927995374D+4
         sy_verify_value = -8.084072988043731D+4
      elseif (m.eq.36) then
         sx_verify_value =  1.982481200946593D+5
         sy_verify_value = -1.020596636361769D+5
      elseif (m.eq.40) then
         sx_verify_value = -5.319717441530D+05
         sy_verify_value = -3.688834557731D+05
      else
         verified = .false.
      endif

      if (verified) then
         sx_err = abs((sx - sx_verify_value)/sx_verify_value)
         sy_err = abs((sy - sy_verify_value)/sy_verify_value)
         verified = ((sx_err.le.epsilon) .and. (sy_err.le.epsilon))
      endif
      Mops = 2.d0**(m+1)/tm/1000000.d0

      write (6,11) tm, m, gc, sx, sy, (i, q(i), i = 0, nq - 1)
 11   format ('EP Benchmark Results:'//'CPU Time =',f10.4/'N = 2^',
     &        i5/'No. Gaussian Pairs =',f15.0/'Sums = ',1p,2d25.15/
     &        'Counts:'/(i3,0p,f15.0))

      call print_results('EP', class, m+1, 0, 0, nit,
     &                   tm, Mops, 
     &                   'Random numbers generated', 
     &                   verified, npbversion, compiletime, cs1,
     &                   cs2, cs3, cs4, cs5, cs6, cs7)


      if (timers_enabled) then
         if (tm .le. 0.d0) tm = 1.0
         tt = timer_read(1)
         print 810, 'Total time:    ', tt, tt*100./tm
         tt = timer_read(2)
         print 810, 'Gaussian pairs:', tt, tt*100./tm
         tt = timer_read(3)
         print 810, 'Random numbers:', tt, tt*100./tm
810      format(1x,a,f9.3,' (',f6.2,'%)')
      endif


      end
      subroutine print_results(name, class, n1, n2, n3, niter, 
     &               t, mops, optype, verified, npbversion, 
     &               compiletime, cs1, cs2, cs3, cs4, cs5, cs6, cs7)

c*********************************************************************72
c
cc PRINT_RESULTS prints the results.
c
      implicit none

      character name*(*)
      character class*1
      integer   n1, n2, n3, niter, j
      double precision t, mops
      character optype*24, size*15
      logical   verified
      character*(*) npbversion, compiletime, 
     &              cs1, cs2, cs3, cs4, cs5, cs6, cs7

         write (*, 2) name 
 2       format(//, ' ', A, ' Benchmark Completed.')

         write (*, 3) Class
 3       format(' Class           = ', 12x, a12)
c
c   If this is not a grid-based problem (EP, FT, CG), then
c   we only print n1, which contains some measure of the
c   problem size. In that case, n2 and n3 are both zero.
c   Otherwise, we print the grid size n1xn2xn3
c
         if ((n2 .eq. 0) .and. (n3 .eq. 0)) then
            if (name(1:2) .eq. 'EP') then
               write(size, '(f15.0)' ) 2.d0**n1
               j = 15
               if (size(j:j) .eq. '.') then
                  size(j:j) = ' '
                  j = j - 1
               endif
               write (*,42) size(1:j)
 42            format(' Size            = ',9x, a15)
            else
               write (*,44) n1
 44            format(' Size            = ',12x, i12)
            endif
         else
            write (*, 4) n1,n2,n3
 4          format(' Size            =  ',9x, i4,'x',i4,'x',i4)
         endif

         write (*, 5) niter
 5       format(' Iterations      = ', 12x, i12)
         
         write (*, 6) t
 6       format(' Time in seconds = ',12x, f12.2)
         
         write (*,9) mops
 9       format(' Mop/s total     = ',12x, f12.2)

         write(*, 11) optype
 11      format(' Operation type  = ', a24)

         if (verified) then 
            write(*,12) '  SUCCESSFUL'
         else
            write(*,12) 'UNSUCCESSFUL'
         endif
 12      format(' Verification    = ', 12x, a)

         write(*,13) npbversion
 13      format(' Version         = ', 12x, a12)

         write(*,14) compiletime
 14      format(' Compile date    = ', 12x, a12)


         write (*,121) cs1
 121     format(/, ' Compile options:', /, 
     >          '    F77          = ', A)

         write (*,122) cs2
 122     format('    FLINK        = ', A)

         write (*,123) cs3
 123     format('    F_LIB        = ', A)

         write (*,124) cs4
 124     format('    F_INC        = ', A)

         write (*,125) cs5
 125     format('    FFLAGS       = ', A)

         write (*,126) cs6
 126     format('    FLINKFLAGS   = ', A)

         write(*, 127) cs7
 127     format('    RAND         = ', A)
        
         write (*,130)
 130     format(//' Please send all errors/feedbacks to:'//
     >            ' NPB Development Team'/
     >            ' npb@nas.nasa.gov'//)

         return
         end
      double precision function randlc(x, a)

c*********************************************************************72
c
cc RANDLC returns a uniform pseudorandom double precision number.
c
c  Discussion:
c
c    The number returned is in the range (0, 1).  
c
c    The algorithm uses the linear congruential generator
c
c   x_{k+1} = a x_k  (mod 2^46)
c
c   where 0 < x_k < 2^46 and 0 < a < 2^46.  This scheme generates 2^44 numbers
c   before repeating.  The argument A is the same as 'a' in the above formula,
c   and X is the same as x_0.  A and X must be odd double precision integers
c   in the range (1, 2^46).  The returned value RANDLC is normalized to be
c   between 0 and 1, i.e. RANDLC = 2^(-46) * x_1.  X is updated to contain
c   the new seed x_1, so that subsequent calls to RANDLC using the same
c   arguments will generate a continuous sequence.
c
      implicit none

      double precision x, a
      integer*8 i246m1, Lx, La
      double precision d2m46

      parameter(d2m46=0.5d0**46)

      save i246m1
      data i246m1/X'00003FFFFFFFFFFF'/

      Lx = X
      La = A

      Lx   = iand(Lx*La,i246m1)
      randlc = d2m46*dble(Lx)
      x    = dble(Lx)
      return
      end
      subroutine vranlc (n, x, a, y)

c*********************************************************************72
c
cc VRANLC returns a vector of uniform pseudorandom double precision numbers.
c
      implicit none

      integer n, i
      double precision x, a, y(*)
      integer*8 i246m1, Lx, La
      double precision d2m46
c
c  This doesn't work, because the compiler does the calculation in 32
c  bits and overflows. No standard way (without f90 stuff) to specify
c  that the rhs should be done in 64 bit arithmetic. 
c      parameter(i246m1=2**46-1)

      parameter(d2m46=0.5d0**46)

      save i246m1
      data i246m1/X'00003FFFFFFFFFFF'/
c
c  Note that the v6 compiler on an R8000 does something stupid with
c  the above. Using the following instead (or various other things)
c  makes the calculation run almost 10 times as fast. 
c 
c      save d2m46
c      data d2m46/0.0d0/
c      if (d2m46 .eq. 0.0d0) then
c         d2m46 = 0.5d0**46
c      endif
c
      Lx = X
      La = A
      do i = 1, N
         Lx   = iand(Lx*La,i246m1)
         y(i) = d2m46*dble(Lx)
      end do
      x    = dble(Lx)

      return
      end      
      subroutine timer_clear ( n )

c*********************************************************************72
c
cc TIMER_CLEAR clears the timer.
c
      implicit none

      integer n
      
      double precision start(64), elapsed(64)
      common /tt/ start, elapsed

      elapsed(n) = 0.0
      return
      end
      subroutine timer_start ( n )

c*********************************************************************72
c
cc TIMER_START starts the timer.
c
      implicit none

      external         elapsed_time
      double precision elapsed_time
      integer n
      double precision start(64), elapsed(64)
      common /tt/ start, elapsed

      start(n) = elapsed_time()

      return
      end
      subroutine timer_stop ( n )

c*********************************************************************72
c
cc TIMER_STOP stops the timer.
c
      implicit none

      external         elapsed_time
      double precision elapsed_time
      integer n
      double precision start(64), elapsed(64)
      common /tt/ start, elapsed
      double precision t, now
      now = elapsed_time()
      t = now - start(n)
      elapsed(n) = elapsed(n) + t

      return
      end
      double precision function timer_read ( n )

c*********************************************************************72
c
cc TIMER_READ reads the timer.
c

      implicit none

      integer n
      double precision start(64), elapsed(64)
      common /tt/ start, elapsed
      
      timer_read = elapsed(n)
      return
      end
      double precision function elapsed_time ( )

c*********************************************************************72
c
cc ELAPSED_TIME measures wall clock time.
c
      implicit none

      double precision t
c
c  This function must measure wall clock time, not CPU time. 
c  Since there is no portable timer in Fortran (77)
c  we call a routine compiled in C (though the C source may have
c  to be tweaked). 
c
      call wtime(t)
c
c  The following is not ok for "official" results because it reports
c  CPU time not wall clock time. It may be useful for developing/testing
c  on timeshared Crays, though. 
c     call second(t)
c
      elapsed_time = t

      return
      end

