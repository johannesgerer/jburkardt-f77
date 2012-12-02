!-------------------------------------------------------------------------!
!                                                                         !
!        N  A  S     P A R A L L E L     B E N C H M A R K S  3.3         !
!                                                                         !
!                      S E R I A L     V E R S I O N                      !
!                                                                         !
!                                   B T                                   !
!                                                                         !
!-------------------------------------------------------------------------!
!                                                                         !
!    This benchmark is a serial version of the NPB BT code.               !
!    Refer to NAS Technical Reports 95-020 and 99-011 for details.        !
!                                                                         !
!    Permission to use, copy, distribute and modify this software         !
!    for any purpose with or without fee is hereby granted.  We           !
!    request, however, that all derived work reference the NAS            !
!    Parallel Benchmarks 3.3. This software is provided "as is"           !
!    without express or implied warranty.                                 !
!                                                                         !
!    Information on NPB 3.3, including the technical report, the          !
!    original specifications, source code, results and information        !
!    on how to submit new results, is available at:                       !
!                                                                         !
!           http://www.nas.nasa.gov/Software/NPB/                         !
!                                                                         !
!    Send comments or suggestions to  npb@nas.nasa.gov                    !
!                                                                         !
!          NAS Parallel Benchmarks Group                                  !
!          NASA Ames Research Center                                      !
!          Mail Stop: T27A-1                                              !
!          Moffett Field, CA   94035-1000                                 !
!                                                                         !
!          E-mail:  npb@nas.nasa.gov                                      !
!          Fax:     (650) 604-3957                                        !
!                                                                         !
!-------------------------------------------------------------------------!

c---------------------------------------------------------------------
c
c Authors: R. Van der Wijngaart
c          T. Harris
c          M. Yarrow
c          H. Jin
c
c---------------------------------------------------------------------

c---------------------------------------------------------------------
       program BT
c---------------------------------------------------------------------

       include  'header.h'
      
       integer i, niter, step, fstatus
       double precision navg, mflops, n3

       external timer_read
       double precision tmax, timer_read, t, trecs(t_last)
       logical verified
       character class
       character        t_names(t_last)*8

c---------------------------------------------------------------------
c      Root node reads input file (if it exists) else takes
c      defaults from parameters
c---------------------------------------------------------------------
          
       open (unit=2,file='timer.flag',status='old', iostat=fstatus)
       if (fstatus .eq. 0) then
         timeron = .true.
         t_names(t_total) = 'total'
         t_names(t_rhsx) = 'rhsx'
         t_names(t_rhsy) = 'rhsy'
         t_names(t_rhsz) = 'rhsz'
         t_names(t_rhs) = 'rhs'
         t_names(t_xsolve) = 'xsolve'
         t_names(t_ysolve) = 'ysolve'
         t_names(t_zsolve) = 'zsolve'
         t_names(t_rdis1) = 'redist1'
         t_names(t_rdis2) = 'redist2'
         t_names(t_add) = 'add'
         close(2)
       else
         timeron = .false.
       endif

       write(*, 1000)
       open (unit=2,file='inputbt.data',status='old', iostat=fstatus)

       if (fstatus .eq. 0) then
         write(*,233) 
 233     format(' Reading from input file inputbt.data')
         read (2,*) niter
         read (2,*) dt
         read (2,*) grid_points(1), grid_points(2), grid_points(3)
         close(2)
       else
         write(*,234) 
         niter = niter_default
         dt    = dt_default
         grid_points(1) = problem_size
         grid_points(2) = problem_size
         grid_points(3) = problem_size
       endif
 234   format(' No input file inputbt.data. Using compiled defaults')

       write(*, 1001) grid_points(1), grid_points(2), grid_points(3)
       write(*, 1002) niter, dt
       write(*, *)

 1000  format(//, ' NAS Parallel Benchmarks (NPB3.3-SER)',
     >            ' - BT Benchmark',/)
 1001  format(' Size: ', i4, 'x', i4, 'x', i4)
 1002  format(' Iterations: ', i4, '    dt: ', F10.6)

       if ( (grid_points(1) .gt. IMAX) .or.
     >      (grid_points(2) .gt. JMAX) .or.
     >      (grid_points(3) .gt. KMAX) ) then
             print *, (grid_points(i),i=1,3)
             print *,' Problem size too big for compiled array sizes'
             goto 999
       endif


       call set_constants

       call initialize
       do i = 1, t_last
          call timer_clear(i)
       end do

c       call lhsinit

       call exact_rhs

c---------------------------------------------------------------------
c      do one time step to touch all code, and reinitialize
c---------------------------------------------------------------------
       call adi
       call initialize

       do i = 1, t_last
          call timer_clear(i)
       end do
       call timer_start(1)

       do  step = 1, niter

          if (mod(step, 20) .eq. 0 .or. 
     >        step .eq. 1) then
             write(*, 200) step
 200         format(' Time step ', i4)
          endif

          call adi

       end do

       call timer_stop(1)
       tmax = timer_read(1)
       
       call verify(niter, class, verified)

       n3 = 1.0d0*grid_points(1)*grid_points(2)*grid_points(3)
       navg = (grid_points(1)+grid_points(2)+grid_points(3))/3.0
       if( tmax .ne. 0. ) then
          mflops = 1.0e-6*float(niter)*
     >  (3478.8*n3-17655.7*navg**2+28023.7*navg)
     >  / tmax
       else
          mflops = 0.0
       endif
       call print_results('BT', class, grid_points(1), 
     >  grid_points(2), grid_points(3), niter,
     >  tmax, mflops, '          floating point', 
     >  verified, npbversion,compiletime, cs1, cs2, cs3, cs4, cs5, 
     >  cs6, '(none)')

c---------------------------------------------------------------------
c      More timers
c---------------------------------------------------------------------
       if (.not.timeron) goto 999

       do i=1, t_last
          trecs(i) = timer_read(i)
       end do

       if (tmax .eq. 0.0) tmax = 1.0
       write(*,800)
 800   format('  SECTION   Time (secs)')
       do i=1, t_last
          write(*,810) t_names(i), trecs(i), trecs(i)*100./tmax
          if (i.eq.t_rhs) then
             t = trecs(t_rhsx) + trecs(t_rhsy) + trecs(t_rhsz)
             write(*,820) 'sub-rhs', t, t*100./tmax
             t = trecs(t_rhs) - t
             write(*,820) 'rest-rhs', t, t*100./tmax
          elseif (i.eq.t_zsolve) then
             t = trecs(t_zsolve) - trecs(t_rdis1) - trecs(t_rdis2)
             write(*,820) 'sub-zsol', t, t*100./tmax
          elseif (i.eq.t_rdis2) then
             t = trecs(t_rdis1) + trecs(t_rdis2)
             write(*,820) 'redist', t, t*100./tmax
          endif
 810      format(2x,a8,':',f9.3,'  (',f6.2,'%)')
 820      format('    --> total ',a8,':',f9.3,'  (',f6.2,'%)')
       end do

 999   continue

       end
      subroutine  add

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     addition of update to the vector u
c---------------------------------------------------------------------

      include 'header.h'

      integer i, j, k, m

      if (timeron) call timer_start(t_add)
      do     k = 1, grid_points(3)-2
         do     j = 1, grid_points(2)-2
            do     i = 1, grid_points(1)-2
               do    m = 1, 5
                  u(m,i,j,k) = u(m,i,j,k) + rhs(m,i,j,k)
               enddo
            enddo
         enddo
      enddo
      if (timeron) call timer_stop(t_add)

      return
      end
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine  adi

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      call compute_rhs

      call x_solve

      call y_solve

      call z_solve

      call add

      return
      end

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine error_norm(rms)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     this function computes the norm of the difference between the
c     computed solution and the exact solution
c---------------------------------------------------------------------

      include 'header.h'

      integer i, j, k, m, d
      double precision xi, eta, zeta, u_exact(5), rms(5), add

      do m = 1, 5 
         rms(m) = 0.0d0
      enddo

      do k = 0, grid_points(3)-1
         zeta = dble(k) * dnzm1
         do j = 0, grid_points(2)-1
            eta = dble(j) * dnym1
            do i = 0, grid_points(1)-1
               xi = dble(i) * dnxm1
               call exact_solution(xi, eta, zeta, u_exact)

               do m = 1, 5
                  add = u(m,i,j,k)-u_exact(m)
                  rms(m) = rms(m) + add*add
               enddo
            enddo
          enddo
       enddo

      do m = 1, 5
         do d = 1, 3
            rms(m) = rms(m) / dble(grid_points(d)-2)
         enddo
         rms(m) = dsqrt(rms(m))
      enddo

      return
      end


c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine rhs_norm(rms)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      include 'header.h'

      integer i, j, k, d, m
      double precision rms(5), add

      do m = 1, 5
         rms(m) = 0.0d0
      enddo 

      do k = 1, grid_points(3)-2
         do j = 1, grid_points(2)-2
            do i = 1, grid_points(1)-2
               do m = 1, 5
                  add = rhs(m,i,j,k)
                  rms(m) = rms(m) + add*add
               enddo 
            enddo 
         enddo 
      enddo 

      do m = 1, 5
         do d = 1, 3
            rms(m) = rms(m) / dble(grid_points(d)-2)
         enddo 
         rms(m) = dsqrt(rms(m))
      enddo 

      return
      end


c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine exact_rhs

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     compute the right hand side based on exact solution
c---------------------------------------------------------------------

      include 'header.h'

      double precision dtemp(5), xi, eta, zeta, dtpp
      integer m, i, j, k, ip1, im1, jp1, jm1, km1, kp1

c---------------------------------------------------------------------
c     initialize                                  
c---------------------------------------------------------------------
      do k= 0, grid_points(3)-1
         do j = 0, grid_points(2)-1
            do i = 0, grid_points(1)-1
               do m = 1, 5
                  forcing(m,i,j,k) = 0.0d0
               enddo
            enddo
         enddo
      enddo

c---------------------------------------------------------------------
c     xi-direction flux differences                      
c---------------------------------------------------------------------
      do k = 1, grid_points(3)-2
         zeta = dble(k) * dnzm1
         do j = 1, grid_points(2)-2
            eta = dble(j) * dnym1

            do i=0, grid_points(1)-1
               xi = dble(i) * dnxm1

               call exact_solution(xi, eta, zeta, dtemp)
               do m = 1, 5
                  ue(i,m) = dtemp(m)
               enddo

               dtpp = 1.0d0 / dtemp(1)

               do m = 2, 5
                  buf(i,m) = dtpp * dtemp(m)
               enddo

               cuf(i)   = buf(i,2) * buf(i,2)
               buf(i,1) = cuf(i) + buf(i,3) * buf(i,3) + 
     >                 buf(i,4) * buf(i,4) 
               q(i) = 0.5d0*(buf(i,2)*ue(i,2) + buf(i,3)*ue(i,3) +
     >                 buf(i,4)*ue(i,4))

            enddo
               
            do i = 1, grid_points(1)-2
               im1 = i-1
               ip1 = i+1

               forcing(1,i,j,k) = forcing(1,i,j,k) -
     >                 tx2*( ue(ip1,2)-ue(im1,2) )+
     >                 dx1tx1*(ue(ip1,1)-2.0d0*ue(i,1)+ue(im1,1))

               forcing(2,i,j,k) = forcing(2,i,j,k) - tx2 * (
     >                 (ue(ip1,2)*buf(ip1,2)+c2*(ue(ip1,5)-q(ip1)))-
     >                 (ue(im1,2)*buf(im1,2)+c2*(ue(im1,5)-q(im1))))+
     >                 xxcon1*(buf(ip1,2)-2.0d0*buf(i,2)+buf(im1,2))+
     >                 dx2tx1*( ue(ip1,2)-2.0d0* ue(i,2)+ue(im1,2))

               forcing(3,i,j,k) = forcing(3,i,j,k) - tx2 * (
     >                 ue(ip1,3)*buf(ip1,2)-ue(im1,3)*buf(im1,2))+
     >                 xxcon2*(buf(ip1,3)-2.0d0*buf(i,3)+buf(im1,3))+
     >                 dx3tx1*( ue(ip1,3)-2.0d0*ue(i,3) +ue(im1,3))
                  
               forcing(4,i,j,k) = forcing(4,i,j,k) - tx2*(
     >                 ue(ip1,4)*buf(ip1,2)-ue(im1,4)*buf(im1,2))+
     >                 xxcon2*(buf(ip1,4)-2.0d0*buf(i,4)+buf(im1,4))+
     >                 dx4tx1*( ue(ip1,4)-2.0d0* ue(i,4)+ ue(im1,4))

               forcing(5,i,j,k) = forcing(5,i,j,k) - tx2*(
     >                 buf(ip1,2)*(c1*ue(ip1,5)-c2*q(ip1))-
     >                 buf(im1,2)*(c1*ue(im1,5)-c2*q(im1)))+
     >                 0.5d0*xxcon3*(buf(ip1,1)-2.0d0*buf(i,1)+
     >                 buf(im1,1))+
     >                 xxcon4*(cuf(ip1)-2.0d0*cuf(i)+cuf(im1))+
     >                 xxcon5*(buf(ip1,5)-2.0d0*buf(i,5)+buf(im1,5))+
     >                 dx5tx1*( ue(ip1,5)-2.0d0* ue(i,5)+ ue(im1,5))
            enddo

c---------------------------------------------------------------------
c     Fourth-order dissipation                         
c---------------------------------------------------------------------

            do m = 1, 5
               i = 1
               forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
     >                    (5.0d0*ue(i,m) - 4.0d0*ue(i+1,m) +ue(i+2,m))
               i = 2
               forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
     >                    (-4.0d0*ue(i-1,m) + 6.0d0*ue(i,m) -
     >                    4.0d0*ue(i+1,m) +       ue(i+2,m))
            enddo

            do i = 3, grid_points(1)-4
               do m = 1, 5
                  forcing(m,i,j,k) = forcing(m,i,j,k) - dssp*
     >                    (ue(i-2,m) - 4.0d0*ue(i-1,m) +
     >                    6.0d0*ue(i,m) - 4.0d0*ue(i+1,m) + ue(i+2,m))
               enddo
            enddo

            do m = 1, 5
               i = grid_points(1)-3
               forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
     >                    (ue(i-2,m) - 4.0d0*ue(i-1,m) +
     >                    6.0d0*ue(i,m) - 4.0d0*ue(i+1,m))
               i = grid_points(1)-2
               forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
     >                    (ue(i-2,m) - 4.0d0*ue(i-1,m) + 5.0d0*ue(i,m))
            enddo

         enddo
      enddo

c---------------------------------------------------------------------
c     eta-direction flux differences             
c---------------------------------------------------------------------
      do k = 1, grid_points(3)-2          
         zeta = dble(k) * dnzm1
         do i=1, grid_points(1)-2
            xi = dble(i) * dnxm1

            do j=0, grid_points(2)-1
               eta = dble(j) * dnym1

               call exact_solution(xi, eta, zeta, dtemp)
               do m = 1, 5 
                  ue(j,m) = dtemp(m)
               enddo
                  
               dtpp = 1.0d0/dtemp(1)

               do m = 2, 5
                  buf(j,m) = dtpp * dtemp(m)
               enddo

               cuf(j)   = buf(j,3) * buf(j,3)
               buf(j,1) = cuf(j) + buf(j,2) * buf(j,2) + 
     >                 buf(j,4) * buf(j,4)
               q(j) = 0.5d0*(buf(j,2)*ue(j,2) + buf(j,3)*ue(j,3) +
     >                 buf(j,4)*ue(j,4))
            enddo

            do j = 1, grid_points(2)-2
               jm1 = j-1
               jp1 = j+1
                  
               forcing(1,i,j,k) = forcing(1,i,j,k) -
     >                 ty2*( ue(jp1,3)-ue(jm1,3) )+
     >                 dy1ty1*(ue(jp1,1)-2.0d0*ue(j,1)+ue(jm1,1))

               forcing(2,i,j,k) = forcing(2,i,j,k) - ty2*(
     >                 ue(jp1,2)*buf(jp1,3)-ue(jm1,2)*buf(jm1,3))+
     >                 yycon2*(buf(jp1,2)-2.0d0*buf(j,2)+buf(jm1,2))+
     >                 dy2ty1*( ue(jp1,2)-2.0* ue(j,2)+ ue(jm1,2))

               forcing(3,i,j,k) = forcing(3,i,j,k) - ty2*(
     >                 (ue(jp1,3)*buf(jp1,3)+c2*(ue(jp1,5)-q(jp1)))-
     >                 (ue(jm1,3)*buf(jm1,3)+c2*(ue(jm1,5)-q(jm1))))+
     >                 yycon1*(buf(jp1,3)-2.0d0*buf(j,3)+buf(jm1,3))+
     >                 dy3ty1*( ue(jp1,3)-2.0d0*ue(j,3) +ue(jm1,3))

               forcing(4,i,j,k) = forcing(4,i,j,k) - ty2*(
     >                 ue(jp1,4)*buf(jp1,3)-ue(jm1,4)*buf(jm1,3))+
     >                 yycon2*(buf(jp1,4)-2.0d0*buf(j,4)+buf(jm1,4))+
     >                 dy4ty1*( ue(jp1,4)-2.0d0*ue(j,4)+ ue(jm1,4))

               forcing(5,i,j,k) = forcing(5,i,j,k) - ty2*(
     >                 buf(jp1,3)*(c1*ue(jp1,5)-c2*q(jp1))-
     >                 buf(jm1,3)*(c1*ue(jm1,5)-c2*q(jm1)))+
     >                 0.5d0*yycon3*(buf(jp1,1)-2.0d0*buf(j,1)+
     >                 buf(jm1,1))+
     >                 yycon4*(cuf(jp1)-2.0d0*cuf(j)+cuf(jm1))+
     >                 yycon5*(buf(jp1,5)-2.0d0*buf(j,5)+buf(jm1,5))+
     >                 dy5ty1*(ue(jp1,5)-2.0d0*ue(j,5)+ue(jm1,5))
            enddo

c---------------------------------------------------------------------
c     Fourth-order dissipation                      
c---------------------------------------------------------------------
            do m = 1, 5
               j = 1
               forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
     >                    (5.0d0*ue(j,m) - 4.0d0*ue(j+1,m) +ue(j+2,m))
               j = 2
               forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
     >                    (-4.0d0*ue(j-1,m) + 6.0d0*ue(j,m) -
     >                    4.0d0*ue(j+1,m) +       ue(j+2,m))
            enddo

            do j = 3, grid_points(2)-4
               do m = 1, 5
                  forcing(m,i,j,k) = forcing(m,i,j,k) - dssp*
     >                    (ue(j-2,m) - 4.0d0*ue(j-1,m) +
     >                    6.0d0*ue(j,m) - 4.0d0*ue(j+1,m) + ue(j+2,m))
               enddo
            enddo

            do m = 1, 5
               j = grid_points(2)-3
               forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
     >                    (ue(j-2,m) - 4.0d0*ue(j-1,m) +
     >                    6.0d0*ue(j,m) - 4.0d0*ue(j+1,m))
               j = grid_points(2)-2
               forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
     >                    (ue(j-2,m) - 4.0d0*ue(j-1,m) + 5.0d0*ue(j,m))

            enddo

         enddo
      enddo

c---------------------------------------------------------------------
c     zeta-direction flux differences                      
c---------------------------------------------------------------------
      do j=1, grid_points(2)-2
         eta = dble(j) * dnym1
         do i = 1, grid_points(1)-2
            xi = dble(i) * dnxm1

            do k=0, grid_points(3)-1
               zeta = dble(k) * dnzm1

               call exact_solution(xi, eta, zeta, dtemp)
               do m = 1, 5
                  ue(k,m) = dtemp(m)
               enddo

               dtpp = 1.0d0/dtemp(1)

               do m = 2, 5
                  buf(k,m) = dtpp * dtemp(m)
               enddo

               cuf(k)   = buf(k,4) * buf(k,4)
               buf(k,1) = cuf(k) + buf(k,2) * buf(k,2) + 
     >                 buf(k,3) * buf(k,3)
               q(k) = 0.5d0*(buf(k,2)*ue(k,2) + buf(k,3)*ue(k,3) +
     >                 buf(k,4)*ue(k,4))
            enddo

            do k=1, grid_points(3)-2
               km1 = k-1
               kp1 = k+1
                  
               forcing(1,i,j,k) = forcing(1,i,j,k) -
     >                 tz2*( ue(kp1,4)-ue(km1,4) )+
     >                 dz1tz1*(ue(kp1,1)-2.0d0*ue(k,1)+ue(km1,1))

               forcing(2,i,j,k) = forcing(2,i,j,k) - tz2 * (
     >                 ue(kp1,2)*buf(kp1,4)-ue(km1,2)*buf(km1,4))+
     >                 zzcon2*(buf(kp1,2)-2.0d0*buf(k,2)+buf(km1,2))+
     >                 dz2tz1*( ue(kp1,2)-2.0d0* ue(k,2)+ ue(km1,2))

               forcing(3,i,j,k) = forcing(3,i,j,k) - tz2 * (
     >                 ue(kp1,3)*buf(kp1,4)-ue(km1,3)*buf(km1,4))+
     >                 zzcon2*(buf(kp1,3)-2.0d0*buf(k,3)+buf(km1,3))+
     >                 dz3tz1*(ue(kp1,3)-2.0d0*ue(k,3)+ue(km1,3))

               forcing(4,i,j,k) = forcing(4,i,j,k) - tz2 * (
     >                 (ue(kp1,4)*buf(kp1,4)+c2*(ue(kp1,5)-q(kp1)))-
     >                 (ue(km1,4)*buf(km1,4)+c2*(ue(km1,5)-q(km1))))+
     >                 zzcon1*(buf(kp1,4)-2.0d0*buf(k,4)+buf(km1,4))+
     >                 dz4tz1*( ue(kp1,4)-2.0d0*ue(k,4) +ue(km1,4))

               forcing(5,i,j,k) = forcing(5,i,j,k) - tz2 * (
     >                 buf(kp1,4)*(c1*ue(kp1,5)-c2*q(kp1))-
     >                 buf(km1,4)*(c1*ue(km1,5)-c2*q(km1)))+
     >                 0.5d0*zzcon3*(buf(kp1,1)-2.0d0*buf(k,1)
     >                 +buf(km1,1))+
     >                 zzcon4*(cuf(kp1)-2.0d0*cuf(k)+cuf(km1))+
     >                 zzcon5*(buf(kp1,5)-2.0d0*buf(k,5)+buf(km1,5))+
     >                 dz5tz1*( ue(kp1,5)-2.0d0*ue(k,5)+ ue(km1,5))
            enddo

c---------------------------------------------------------------------
c     Fourth-order dissipation                        
c---------------------------------------------------------------------
            do m = 1, 5
               k = 1
               forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
     >                    (5.0d0*ue(k,m) - 4.0d0*ue(k+1,m) +ue(k+2,m))
               k = 2
               forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
     >                    (-4.0d0*ue(k-1,m) + 6.0d0*ue(k,m) -
     >                    4.0d0*ue(k+1,m) +       ue(k+2,m))
            enddo

            do k = 3, grid_points(3)-4
               do m = 1, 5
                  forcing(m,i,j,k) = forcing(m,i,j,k) - dssp*
     >                    (ue(k-2,m) - 4.0d0*ue(k-1,m) +
     >                    6.0d0*ue(k,m) - 4.0d0*ue(k+1,m) + ue(k+2,m))
               enddo
            enddo

            do m = 1, 5
               k = grid_points(3)-3
               forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
     >                    (ue(k-2,m) - 4.0d0*ue(k-1,m) +
     >                    6.0d0*ue(k,m) - 4.0d0*ue(k+1,m))
               k = grid_points(3)-2
               forcing(m,i,j,k) = forcing(m,i,j,k) - dssp *
     >                    (ue(k-2,m) - 4.0d0*ue(k-1,m) + 5.0d0*ue(k,m))
            enddo

         enddo
      enddo

c---------------------------------------------------------------------
c     now change the sign of the forcing function, 
c---------------------------------------------------------------------
      do k = 1, grid_points(3)-2
         do j = 1, grid_points(2)-2
            do i = 1, grid_points(1)-2
               do m = 1, 5
                  forcing(m,i,j,k) = -1.d0 * forcing(m,i,j,k)
               enddo
            enddo
         enddo
      enddo


      return
      end
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine exact_solution(xi,eta,zeta,dtemp)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     this function returns the exact solution at point xi, eta, zeta  
c---------------------------------------------------------------------

      include 'header.h'

      double precision  xi, eta, zeta, dtemp(5)
      integer m

      do m = 1, 5
         dtemp(m) =  ce(m,1) +
     >     xi*(ce(m,2) + xi*(ce(m,5) + xi*(ce(m,8) + xi*ce(m,11)))) +
     >     eta*(ce(m,3) + eta*(ce(m,6) + eta*(ce(m,9) + eta*ce(m,12))))+
     >     zeta*(ce(m,4) + zeta*(ce(m,7) + zeta*(ce(m,10) + 
     >     zeta*ce(m,13))))
      enddo

      return
      end


c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine  initialize

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     This subroutine initializes the field variable u using 
c     tri-linear transfinite interpolation of the boundary values     
c---------------------------------------------------------------------

      include 'header.h'
      
      integer i, j, k, m, ix, iy, iz
      double precision  xi, eta, zeta, Pface(5,3,2), Pxi, Peta, 
     >     Pzeta, temp(5)

c---------------------------------------------------------------------
c  Later (in compute_rhs) we compute 1/u for every element. A few of 
c  the corner elements are not used, but it convenient (and faster) 
c  to compute the whole thing with a simple loop. Make sure those 
c  values are nonzero by initializing the whole thing here. 
c---------------------------------------------------------------------
      do k = 0, grid_points(3)-1
         do j = 0, grid_points(2)-1
            do i = 0, grid_points(1)-1
               do m = 1, 5
                  u(m,i,j,k) = 1.0
               end do
            end do
         end do
      end do
c---------------------------------------------------------------------



c---------------------------------------------------------------------
c     first store the "interpolated" values everywhere on the grid    
c---------------------------------------------------------------------

      do k = 0, grid_points(3)-1
         zeta = dble(k) * dnzm1
         do j = 0, grid_points(2)-1
            eta = dble(j) * dnym1
            do i = 0, grid_points(1)-1
               xi = dble(i) * dnxm1
                  
               do ix = 1, 2
                  call exact_solution(dble(ix-1), eta, zeta, 
     >                    Pface(1,1,ix))
               enddo

               do iy = 1, 2
                  call exact_solution(xi, dble(iy-1) , zeta, 
     >                    Pface(1,2,iy))
               enddo

               do iz = 1, 2
                  call exact_solution(xi, eta, dble(iz-1),   
     >                    Pface(1,3,iz))
               enddo

               do m = 1, 5
                  Pxi   = xi   * Pface(m,1,2) + 
     >                    (1.0d0-xi)   * Pface(m,1,1)
                  Peta  = eta  * Pface(m,2,2) + 
     >                    (1.0d0-eta)  * Pface(m,2,1)
                  Pzeta = zeta * Pface(m,3,2) + 
     >                    (1.0d0-zeta) * Pface(m,3,1)
                     
                  u(m,i,j,k) = Pxi + Peta + Pzeta - 
     >                    Pxi*Peta - Pxi*Pzeta - Peta*Pzeta + 
     >                    Pxi*Peta*Pzeta

               enddo
            enddo
         enddo
      enddo

c---------------------------------------------------------------------
c     now store the exact values on the boundaries        
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     west face                                                  
c---------------------------------------------------------------------
      i = 0
      xi = 0.0d0
      do k = 0, grid_points(3)-1
         zeta = dble(k) * dnzm1
         do j = 0, grid_points(2)-1
            eta = dble(j) * dnym1
            call exact_solution(xi, eta, zeta, temp)
            do m = 1, 5
               u(m,i,j,k) = temp(m)
            enddo
         enddo
      enddo

c---------------------------------------------------------------------
c     east face                                                      
c---------------------------------------------------------------------

      i = grid_points(1)-1
      xi = 1.0d0
      do k = 0, grid_points(3)-1
         zeta = dble(k) * dnzm1
         do j = 0, grid_points(2)-1
            eta = dble(j) * dnym1
            call exact_solution(xi, eta, zeta, temp)
            do m = 1, 5
               u(m,i,j,k) = temp(m)
            enddo
         enddo
      enddo

c---------------------------------------------------------------------
c     south face                                                 
c---------------------------------------------------------------------
      j = 0
      eta = 0.0d0
      do k = 0, grid_points(3)-1
         zeta = dble(k) * dnzm1
         do i = 0, grid_points(1)-1
            xi = dble(i) * dnxm1
            call exact_solution(xi, eta, zeta, temp)
            do m = 1, 5
               u(m,i,j,k) = temp(m)
            enddo
         enddo
      enddo


c---------------------------------------------------------------------
c     north face                                    
c---------------------------------------------------------------------
      j = grid_points(2)-1
      eta = 1.0d0
      do k = 0, grid_points(3)-1
         zeta = dble(k) * dnzm1
         do i = 0, grid_points(1)-1
            xi = dble(i) * dnxm1
            call exact_solution(xi, eta, zeta, temp)
            do m = 1, 5
               u(m,i,j,k) = temp(m)
            enddo
         enddo
      enddo

c---------------------------------------------------------------------
c     bottom face                                       
c---------------------------------------------------------------------
      k = 0
      zeta = 0.0d0
      do j = 0, grid_points(2)-1
         eta = dble(j) * dnym1
         do i =0, grid_points(1)-1
            xi = dble(i) *dnxm1
            call exact_solution(xi, eta, zeta, temp)
            do m = 1, 5
               u(m,i,j,k) = temp(m)
            enddo
         enddo
      enddo

c---------------------------------------------------------------------
c     top face     
c---------------------------------------------------------------------
      k = grid_points(3)-1
      zeta = 1.0d0
      do j = 0, grid_points(2)-1
         eta = dble(j) * dnym1
         do i =0, grid_points(1)-1
            xi = dble(i) * dnxm1
            call exact_solution(xi, eta, zeta, temp)
            do m = 1, 5
               u(m,i,j,k) = temp(m)
            enddo
         enddo
      enddo

      return
      end


c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine lhsinit(lhs, size)
      implicit none
      integer size
      double precision lhs(5,5,3,0:size)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      integer i, m, n

      i = size
c---------------------------------------------------------------------
c     zero the whole left hand side for starters
c---------------------------------------------------------------------
      do m = 1, 5
         do n = 1, 5
            lhs(m,n,1,0) = 0.0d0
            lhs(m,n,2,0) = 0.0d0
            lhs(m,n,3,0) = 0.0d0
            lhs(m,n,1,i) = 0.0d0
            lhs(m,n,2,i) = 0.0d0
            lhs(m,n,3,i) = 0.0d0
         enddo
      enddo

c---------------------------------------------------------------------
c     next, set all diagonal values to 1. This is overkill, but convenient
c---------------------------------------------------------------------
      do m = 1, 5
         lhs(m,m,2,0) = 1.0d0
         lhs(m,m,2,i) = 1.0d0
      enddo

      return
      end




      subroutine print_results(name, class, n1, n2, n3, niter, 
     >               t, mops, optype, verified, npbversion, 
     >               compiletime, cs1, cs2, cs3, cs4, cs5, cs6, cs7)
      
      implicit none
      character name*(*)
      character class*1
      integer   n1, n2, n3, niter, j
      double precision t, mops
      character optype*24, size*15
      logical   verified
      character*(*) npbversion, compiletime, 
     >              cs1, cs2, cs3, cs4, cs5, cs6, cs7

         write (*, 2) name 
 2       format(//, ' ', A, ' Benchmark Completed.')

         write (*, 3) Class
 3       format(' Class           = ', 12x, a12)

c   If this is not a grid-based problem (EP, FT, CG), then
c   we only print n1, which contains some measure of the
c   problem size. In that case, n2 and n3 are both zero.
c   Otherwise, we print the grid size n1xn2xn3

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

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine compute_rhs

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      include 'header.h'

      integer i, j, k, m
      double precision rho_inv, uijk, up1, um1, vijk, vp1, vm1,
     >     wijk, wp1, wm1


      if (timeron) call timer_start(t_rhs)
c---------------------------------------------------------------------
c     compute the reciprocal of density, and the kinetic energy, 
c     and the speed of sound.
c---------------------------------------------------------------------
      do k = 0, grid_points(3)-1
         do j = 0, grid_points(2)-1
            do i = 0, grid_points(1)-1
               rho_inv = 1.0d0/u(1,i,j,k)
               rho_i(i,j,k) = rho_inv
               us(i,j,k) = u(2,i,j,k) * rho_inv
               vs(i,j,k) = u(3,i,j,k) * rho_inv
               ws(i,j,k) = u(4,i,j,k) * rho_inv
               square(i,j,k)     = 0.5d0* (
     >                 u(2,i,j,k)*u(2,i,j,k) + 
     >                 u(3,i,j,k)*u(3,i,j,k) +
     >                 u(4,i,j,k)*u(4,i,j,k) ) * rho_inv
               qs(i,j,k) = square(i,j,k) * rho_inv
            enddo
         enddo
      enddo

c---------------------------------------------------------------------
c copy the exact forcing term to the right hand side;  because 
c this forcing term is known, we can store it on the whole grid
c including the boundary                   
c---------------------------------------------------------------------

      do k = 0, grid_points(3)-1
         do j = 0, grid_points(2)-1
            do i = 0, grid_points(1)-1
               do m = 1, 5
                  rhs(m,i,j,k) = forcing(m,i,j,k)
               enddo
            enddo
         enddo
      enddo


      if (timeron) call timer_start(t_rhsx)
c---------------------------------------------------------------------
c     compute xi-direction fluxes 
c---------------------------------------------------------------------
      do k = 1, grid_points(3)-2
         do j = 1, grid_points(2)-2
            do i = 1, grid_points(1)-2
               uijk = us(i,j,k)
               up1  = us(i+1,j,k)
               um1  = us(i-1,j,k)

               rhs(1,i,j,k) = rhs(1,i,j,k) + dx1tx1 * 
     >                 (u(1,i+1,j,k) - 2.0d0*u(1,i,j,k) + 
     >                 u(1,i-1,j,k)) -
     >                 tx2 * (u(2,i+1,j,k) - u(2,i-1,j,k))

               rhs(2,i,j,k) = rhs(2,i,j,k) + dx2tx1 * 
     >                 (u(2,i+1,j,k) - 2.0d0*u(2,i,j,k) + 
     >                 u(2,i-1,j,k)) +
     >                 xxcon2*con43 * (up1 - 2.0d0*uijk + um1) -
     >                 tx2 * (u(2,i+1,j,k)*up1 - 
     >                 u(2,i-1,j,k)*um1 +
     >                 (u(5,i+1,j,k)- square(i+1,j,k)-
     >                 u(5,i-1,j,k)+ square(i-1,j,k))*
     >                 c2)

               rhs(3,i,j,k) = rhs(3,i,j,k) + dx3tx1 * 
     >                 (u(3,i+1,j,k) - 2.0d0*u(3,i,j,k) +
     >                 u(3,i-1,j,k)) +
     >                 xxcon2 * (vs(i+1,j,k) - 2.0d0*vs(i,j,k) +
     >                 vs(i-1,j,k)) -
     >                 tx2 * (u(3,i+1,j,k)*up1 - 
     >                 u(3,i-1,j,k)*um1)

               rhs(4,i,j,k) = rhs(4,i,j,k) + dx4tx1 * 
     >                 (u(4,i+1,j,k) - 2.0d0*u(4,i,j,k) +
     >                 u(4,i-1,j,k)) +
     >                 xxcon2 * (ws(i+1,j,k) - 2.0d0*ws(i,j,k) +
     >                 ws(i-1,j,k)) -
     >                 tx2 * (u(4,i+1,j,k)*up1 - 
     >                 u(4,i-1,j,k)*um1)

               rhs(5,i,j,k) = rhs(5,i,j,k) + dx5tx1 * 
     >                 (u(5,i+1,j,k) - 2.0d0*u(5,i,j,k) +
     >                 u(5,i-1,j,k)) +
     >                 xxcon3 * (qs(i+1,j,k) - 2.0d0*qs(i,j,k) +
     >                 qs(i-1,j,k)) +
     >                 xxcon4 * (up1*up1 -       2.0d0*uijk*uijk + 
     >                 um1*um1) +
     >                 xxcon5 * (u(5,i+1,j,k)*rho_i(i+1,j,k) - 
     >                 2.0d0*u(5,i,j,k)*rho_i(i,j,k) +
     >                 u(5,i-1,j,k)*rho_i(i-1,j,k)) -
     >                 tx2 * ( (c1*u(5,i+1,j,k) - 
     >                 c2*square(i+1,j,k))*up1 -
     >                 (c1*u(5,i-1,j,k) - 
     >                 c2*square(i-1,j,k))*um1 )
            enddo
         enddo

c---------------------------------------------------------------------
c     add fourth order xi-direction dissipation               
c---------------------------------------------------------------------
         do j = 1, grid_points(2)-2
            i = 1
            do m = 1, 5
               rhs(m,i,j,k) = rhs(m,i,j,k)- dssp * 
     >                    ( 5.0d0*u(m,i,j,k) - 4.0d0*u(m,i+1,j,k) +
     >                    u(m,i+2,j,k))
            enddo

            i = 2
            do m = 1, 5
               rhs(m,i,j,k) = rhs(m,i,j,k) - dssp * 
     >                    (-4.0d0*u(m,i-1,j,k) + 6.0d0*u(m,i,j,k) -
     >                    4.0d0*u(m,i+1,j,k) + u(m,i+2,j,k))
            enddo
         enddo

         do j = 1, grid_points(2)-2
            do i = 3,grid_points(1)-4
               do m = 1, 5
                  rhs(m,i,j,k) = rhs(m,i,j,k) - dssp * 
     >                    (  u(m,i-2,j,k) - 4.0d0*u(m,i-1,j,k) + 
     >                    6.0*u(m,i,j,k) - 4.0d0*u(m,i+1,j,k) + 
     >                    u(m,i+2,j,k) )
               enddo
            enddo
         enddo
         
         do j = 1, grid_points(2)-2
            i = grid_points(1)-3
            do m = 1, 5
               rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
     >                    ( u(m,i-2,j,k) - 4.0d0*u(m,i-1,j,k) + 
     >                    6.0d0*u(m,i,j,k) - 4.0d0*u(m,i+1,j,k) )
            enddo

            i = grid_points(1)-2
            do m = 1, 5
               rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
     >                    ( u(m,i-2,j,k) - 4.d0*u(m,i-1,j,k) +
     >                    5.d0*u(m,i,j,k) )
            enddo
         enddo
      enddo
      if (timeron) call timer_stop(t_rhsx)

      if (timeron) call timer_start(t_rhsy)
c---------------------------------------------------------------------
c     compute eta-direction fluxes 
c---------------------------------------------------------------------
      do k = 1, grid_points(3)-2
         do j = 1, grid_points(2)-2
            do i = 1, grid_points(1)-2
               vijk = vs(i,j,k)
               vp1  = vs(i,j+1,k)
               vm1  = vs(i,j-1,k)
               rhs(1,i,j,k) = rhs(1,i,j,k) + dy1ty1 * 
     >                 (u(1,i,j+1,k) - 2.0d0*u(1,i,j,k) + 
     >                 u(1,i,j-1,k)) -
     >                 ty2 * (u(3,i,j+1,k) - u(3,i,j-1,k))
               rhs(2,i,j,k) = rhs(2,i,j,k) + dy2ty1 * 
     >                 (u(2,i,j+1,k) - 2.0d0*u(2,i,j,k) + 
     >                 u(2,i,j-1,k)) +
     >                 yycon2 * (us(i,j+1,k) - 2.0d0*us(i,j,k) + 
     >                 us(i,j-1,k)) -
     >                 ty2 * (u(2,i,j+1,k)*vp1 - 
     >                 u(2,i,j-1,k)*vm1)
               rhs(3,i,j,k) = rhs(3,i,j,k) + dy3ty1 * 
     >                 (u(3,i,j+1,k) - 2.0d0*u(3,i,j,k) + 
     >                 u(3,i,j-1,k)) +
     >                 yycon2*con43 * (vp1 - 2.0d0*vijk + vm1) -
     >                 ty2 * (u(3,i,j+1,k)*vp1 - 
     >                 u(3,i,j-1,k)*vm1 +
     >                 (u(5,i,j+1,k) - square(i,j+1,k) - 
     >                 u(5,i,j-1,k) + square(i,j-1,k))
     >                 *c2)
               rhs(4,i,j,k) = rhs(4,i,j,k) + dy4ty1 * 
     >                 (u(4,i,j+1,k) - 2.0d0*u(4,i,j,k) + 
     >                 u(4,i,j-1,k)) +
     >                 yycon2 * (ws(i,j+1,k) - 2.0d0*ws(i,j,k) + 
     >                 ws(i,j-1,k)) -
     >                 ty2 * (u(4,i,j+1,k)*vp1 - 
     >                 u(4,i,j-1,k)*vm1)
               rhs(5,i,j,k) = rhs(5,i,j,k) + dy5ty1 * 
     >                 (u(5,i,j+1,k) - 2.0d0*u(5,i,j,k) + 
     >                 u(5,i,j-1,k)) +
     >                 yycon3 * (qs(i,j+1,k) - 2.0d0*qs(i,j,k) + 
     >                 qs(i,j-1,k)) +
     >                 yycon4 * (vp1*vp1       - 2.0d0*vijk*vijk + 
     >                 vm1*vm1) +
     >                 yycon5 * (u(5,i,j+1,k)*rho_i(i,j+1,k) - 
     >                 2.0d0*u(5,i,j,k)*rho_i(i,j,k) +
     >                 u(5,i,j-1,k)*rho_i(i,j-1,k)) -
     >                 ty2 * ((c1*u(5,i,j+1,k) - 
     >                 c2*square(i,j+1,k)) * vp1 -
     >                 (c1*u(5,i,j-1,k) - 
     >                 c2*square(i,j-1,k)) * vm1)
            enddo
         enddo

c---------------------------------------------------------------------
c     add fourth order eta-direction dissipation         
c---------------------------------------------------------------------
         j = 1
         do i = 1, grid_points(1)-2
            do m = 1, 5
               rhs(m,i,j,k) = rhs(m,i,j,k)- dssp * 
     >                    ( 5.0d0*u(m,i,j,k) - 4.0d0*u(m,i,j+1,k) +
     >                    u(m,i,j+2,k))
            enddo
         enddo

         j = 2
         do i = 1, grid_points(1)-2
            do m = 1, 5
               rhs(m,i,j,k) = rhs(m,i,j,k) - dssp * 
     >                    (-4.0d0*u(m,i,j-1,k) + 6.0d0*u(m,i,j,k) -
     >                    4.0d0*u(m,i,j+1,k) + u(m,i,j+2,k))
            enddo
         enddo

         do j = 3, grid_points(2)-4
            do i = 1,grid_points(1)-2
               do m = 1, 5
                  rhs(m,i,j,k) = rhs(m,i,j,k) - dssp * 
     >                    (  u(m,i,j-2,k) - 4.0d0*u(m,i,j-1,k) + 
     >                    6.0*u(m,i,j,k) - 4.0d0*u(m,i,j+1,k) + 
     >                    u(m,i,j+2,k) )
               enddo
            enddo
         enddo
         
         j = grid_points(2)-3
         do i = 1, grid_points(1)-2
            do m = 1, 5
               rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
     >                    ( u(m,i,j-2,k) - 4.0d0*u(m,i,j-1,k) + 
     >                    6.0d0*u(m,i,j,k) - 4.0d0*u(m,i,j+1,k) )
            enddo
         enddo

         j = grid_points(2)-2
         do i = 1, grid_points(1)-2
            do m = 1, 5
               rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
     >                    ( u(m,i,j-2,k) - 4.d0*u(m,i,j-1,k) +
     >                    5.d0*u(m,i,j,k) )
            enddo
         enddo
      enddo
      if (timeron) call timer_stop(t_rhsy)

      if (timeron) call timer_start(t_rhsz)
c---------------------------------------------------------------------
c     compute zeta-direction fluxes 
c---------------------------------------------------------------------
      do k = 1, grid_points(3)-2
         do j = 1, grid_points(2)-2
            do i = 1, grid_points(1)-2
               wijk = ws(i,j,k)
               wp1  = ws(i,j,k+1)
               wm1  = ws(i,j,k-1)

               rhs(1,i,j,k) = rhs(1,i,j,k) + dz1tz1 * 
     >                 (u(1,i,j,k+1) - 2.0d0*u(1,i,j,k) + 
     >                 u(1,i,j,k-1)) -
     >                 tz2 * (u(4,i,j,k+1) - u(4,i,j,k-1))
               rhs(2,i,j,k) = rhs(2,i,j,k) + dz2tz1 * 
     >                 (u(2,i,j,k+1) - 2.0d0*u(2,i,j,k) + 
     >                 u(2,i,j,k-1)) +
     >                 zzcon2 * (us(i,j,k+1) - 2.0d0*us(i,j,k) + 
     >                 us(i,j,k-1)) -
     >                 tz2 * (u(2,i,j,k+1)*wp1 - 
     >                 u(2,i,j,k-1)*wm1)
               rhs(3,i,j,k) = rhs(3,i,j,k) + dz3tz1 * 
     >                 (u(3,i,j,k+1) - 2.0d0*u(3,i,j,k) + 
     >                 u(3,i,j,k-1)) +
     >                 zzcon2 * (vs(i,j,k+1) - 2.0d0*vs(i,j,k) + 
     >                 vs(i,j,k-1)) -
     >                 tz2 * (u(3,i,j,k+1)*wp1 - 
     >                 u(3,i,j,k-1)*wm1)
               rhs(4,i,j,k) = rhs(4,i,j,k) + dz4tz1 * 
     >                 (u(4,i,j,k+1) - 2.0d0*u(4,i,j,k) + 
     >                 u(4,i,j,k-1)) +
     >                 zzcon2*con43 * (wp1 - 2.0d0*wijk + wm1) -
     >                 tz2 * (u(4,i,j,k+1)*wp1 - 
     >                 u(4,i,j,k-1)*wm1 +
     >                 (u(5,i,j,k+1) - square(i,j,k+1) - 
     >                 u(5,i,j,k-1) + square(i,j,k-1))
     >                 *c2)
               rhs(5,i,j,k) = rhs(5,i,j,k) + dz5tz1 * 
     >                 (u(5,i,j,k+1) - 2.0d0*u(5,i,j,k) + 
     >                 u(5,i,j,k-1)) +
     >                 zzcon3 * (qs(i,j,k+1) - 2.0d0*qs(i,j,k) + 
     >                 qs(i,j,k-1)) +
     >                 zzcon4 * (wp1*wp1 - 2.0d0*wijk*wijk + 
     >                 wm1*wm1) +
     >                 zzcon5 * (u(5,i,j,k+1)*rho_i(i,j,k+1) - 
     >                 2.0d0*u(5,i,j,k)*rho_i(i,j,k) +
     >                 u(5,i,j,k-1)*rho_i(i,j,k-1)) -
     >                 tz2 * ( (c1*u(5,i,j,k+1) - 
     >                 c2*square(i,j,k+1))*wp1 -
     >                 (c1*u(5,i,j,k-1) - 
     >                 c2*square(i,j,k-1))*wm1)
            enddo
         enddo
      enddo

c---------------------------------------------------------------------
c     add fourth order zeta-direction dissipation                
c---------------------------------------------------------------------
      k = 1
      do j = 1, grid_points(2)-2
         do i = 1, grid_points(1)-2
            do m = 1, 5
               rhs(m,i,j,k) = rhs(m,i,j,k)- dssp * 
     >                    ( 5.0d0*u(m,i,j,k) - 4.0d0*u(m,i,j,k+1) +
     >                    u(m,i,j,k+2))
            enddo
         enddo
      enddo

      k = 2
      do j = 1, grid_points(2)-2
         do i = 1, grid_points(1)-2
            do m = 1, 5
               rhs(m,i,j,k) = rhs(m,i,j,k) - dssp * 
     >                    (-4.0d0*u(m,i,j,k-1) + 6.0d0*u(m,i,j,k) -
     >                    4.0d0*u(m,i,j,k+1) + u(m,i,j,k+2))
            enddo
         enddo
      enddo

      do k = 3, grid_points(3)-4
         do j = 1, grid_points(2)-2
            do i = 1,grid_points(1)-2
               do m = 1, 5
                  rhs(m,i,j,k) = rhs(m,i,j,k) - dssp * 
     >                    (  u(m,i,j,k-2) - 4.0d0*u(m,i,j,k-1) + 
     >                    6.0*u(m,i,j,k) - 4.0d0*u(m,i,j,k+1) + 
     >                    u(m,i,j,k+2) )
               enddo
            enddo
         enddo
      enddo
         
      k = grid_points(3)-3
      do j = 1, grid_points(2)-2
         do i = 1, grid_points(1)-2
            do m = 1, 5
               rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
     >                    ( u(m,i,j,k-2) - 4.0d0*u(m,i,j,k-1) + 
     >                    6.0d0*u(m,i,j,k) - 4.0d0*u(m,i,j,k+1) )
            enddo
         enddo
      enddo

      k = grid_points(3)-2
      do j = 1, grid_points(2)-2
         do i = 1, grid_points(1)-2
            do m = 1, 5
               rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
     >                    ( u(m,i,j,k-2) - 4.d0*u(m,i,j,k-1) +
     >                    5.d0*u(m,i,j,k) )
            enddo
         enddo
      enddo
      if (timeron) call timer_stop(t_rhsz)

      do k = 1, grid_points(3)-2
         do j = 1, grid_points(2)-2
            do i = 1, grid_points(1)-2
               do m = 1, 5
                  rhs(m,i,j,k) = rhs(m,i,j,k) * dt
               enddo
            enddo
         enddo
      enddo
      if (timeron) call timer_stop(t_rhs)

      return
      end




c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine  set_constants

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      include 'header.h'
      
      ce(1,1)  = 2.0d0
      ce(1,2)  = 0.0d0
      ce(1,3)  = 0.0d0
      ce(1,4)  = 4.0d0
      ce(1,5)  = 5.0d0
      ce(1,6)  = 3.0d0
      ce(1,7)  = 0.5d0
      ce(1,8)  = 0.02d0
      ce(1,9)  = 0.01d0
      ce(1,10) = 0.03d0
      ce(1,11) = 0.5d0
      ce(1,12) = 0.4d0
      ce(1,13) = 0.3d0
      
      ce(2,1)  = 1.0d0
      ce(2,2)  = 0.0d0
      ce(2,3)  = 0.0d0
      ce(2,4)  = 0.0d0
      ce(2,5)  = 1.0d0
      ce(2,6)  = 2.0d0
      ce(2,7)  = 3.0d0
      ce(2,8)  = 0.01d0
      ce(2,9)  = 0.03d0
      ce(2,10) = 0.02d0
      ce(2,11) = 0.4d0
      ce(2,12) = 0.3d0
      ce(2,13) = 0.5d0

      ce(3,1)  = 2.0d0
      ce(3,2)  = 2.0d0
      ce(3,3)  = 0.0d0
      ce(3,4)  = 0.0d0
      ce(3,5)  = 0.0d0
      ce(3,6)  = 2.0d0
      ce(3,7)  = 3.0d0
      ce(3,8)  = 0.04d0
      ce(3,9)  = 0.03d0
      ce(3,10) = 0.05d0
      ce(3,11) = 0.3d0
      ce(3,12) = 0.5d0
      ce(3,13) = 0.4d0

      ce(4,1)  = 2.0d0
      ce(4,2)  = 2.0d0
      ce(4,3)  = 0.0d0
      ce(4,4)  = 0.0d0
      ce(4,5)  = 0.0d0
      ce(4,6)  = 2.0d0
      ce(4,7)  = 3.0d0
      ce(4,8)  = 0.03d0
      ce(4,9)  = 0.05d0
      ce(4,10) = 0.04d0
      ce(4,11) = 0.2d0
      ce(4,12) = 0.1d0
      ce(4,13) = 0.3d0

      ce(5,1)  = 5.0d0
      ce(5,2)  = 4.0d0
      ce(5,3)  = 3.0d0
      ce(5,4)  = 2.0d0
      ce(5,5)  = 0.1d0
      ce(5,6)  = 0.4d0
      ce(5,7)  = 0.3d0
      ce(5,8)  = 0.05d0
      ce(5,9)  = 0.04d0
      ce(5,10) = 0.03d0
      ce(5,11) = 0.1d0
      ce(5,12) = 0.3d0
      ce(5,13) = 0.2d0

      c1 = 1.4d0
      c2 = 0.4d0
      c3 = 0.1d0
      c4 = 1.0d0
      c5 = 1.4d0

      dnxm1 = 1.0d0 / dble(grid_points(1)-1)
      dnym1 = 1.0d0 / dble(grid_points(2)-1)
      dnzm1 = 1.0d0 / dble(grid_points(3)-1)

      c1c2 = c1 * c2
      c1c5 = c1 * c5
      c3c4 = c3 * c4
      c1345 = c1c5 * c3c4

      conz1 = (1.0d0-c1c5)

      tx1 = 1.0d0 / (dnxm1 * dnxm1)
      tx2 = 1.0d0 / (2.0d0 * dnxm1)
      tx3 = 1.0d0 / dnxm1

      ty1 = 1.0d0 / (dnym1 * dnym1)
      ty2 = 1.0d0 / (2.0d0 * dnym1)
      ty3 = 1.0d0 / dnym1
      
      tz1 = 1.0d0 / (dnzm1 * dnzm1)
      tz2 = 1.0d0 / (2.0d0 * dnzm1)
      tz3 = 1.0d0 / dnzm1

      dx1 = 0.75d0
      dx2 = 0.75d0
      dx3 = 0.75d0
      dx4 = 0.75d0
      dx5 = 0.75d0

      dy1 = 0.75d0
      dy2 = 0.75d0
      dy3 = 0.75d0
      dy4 = 0.75d0
      dy5 = 0.75d0

      dz1 = 1.0d0
      dz2 = 1.0d0
      dz3 = 1.0d0
      dz4 = 1.0d0
      dz5 = 1.0d0

      dxmax = dmax1(dx3, dx4)
      dymax = dmax1(dy2, dy4)
      dzmax = dmax1(dz2, dz3)

      dssp = 0.25d0 * dmax1(dx1, dmax1(dy1, dz1) )

      c4dssp = 4.0d0 * dssp
      c5dssp = 5.0d0 * dssp

      dttx1 = dt*tx1
      dttx2 = dt*tx2
      dtty1 = dt*ty1
      dtty2 = dt*ty2
      dttz1 = dt*tz1
      dttz2 = dt*tz2

      c2dttx1 = 2.0d0*dttx1
      c2dtty1 = 2.0d0*dtty1
      c2dttz1 = 2.0d0*dttz1

      dtdssp = dt*dssp

      comz1  = dtdssp
      comz4  = 4.0d0*dtdssp
      comz5  = 5.0d0*dtdssp
      comz6  = 6.0d0*dtdssp

      c3c4tx3 = c3c4*tx3
      c3c4ty3 = c3c4*ty3
      c3c4tz3 = c3c4*tz3

      dx1tx1 = dx1*tx1
      dx2tx1 = dx2*tx1
      dx3tx1 = dx3*tx1
      dx4tx1 = dx4*tx1
      dx5tx1 = dx5*tx1
      
      dy1ty1 = dy1*ty1
      dy2ty1 = dy2*ty1
      dy3ty1 = dy3*ty1
      dy4ty1 = dy4*ty1
      dy5ty1 = dy5*ty1
      
      dz1tz1 = dz1*tz1
      dz2tz1 = dz2*tz1
      dz3tz1 = dz3*tz1
      dz4tz1 = dz4*tz1
      dz5tz1 = dz5*tz1

      c2iv  = 2.5d0
      con43 = 4.0d0/3.0d0
      con16 = 1.0d0/6.0d0
      
      xxcon1 = c3c4tx3*con43*tx3
      xxcon2 = c3c4tx3*tx3
      xxcon3 = c3c4tx3*conz1*tx3
      xxcon4 = c3c4tx3*con16*tx3
      xxcon5 = c3c4tx3*c1c5*tx3

      yycon1 = c3c4ty3*con43*ty3
      yycon2 = c3c4ty3*ty3
      yycon3 = c3c4ty3*conz1*ty3
      yycon4 = c3c4ty3*con16*ty3
      yycon5 = c3c4ty3*c1c5*ty3

      zzcon1 = c3c4tz3*con43*tz3
      zzcon2 = c3c4tz3*tz3
      zzcon3 = c3c4tz3*conz1*tz3
      zzcon4 = c3c4tz3*con16*tz3
      zzcon5 = c3c4tz3*c1c5*tz3

      return
      end

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine matvec_sub(ablock,avec,bvec)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     subtracts bvec=bvec - ablock*avec
c---------------------------------------------------------------------

      implicit none

      double precision ablock,avec,bvec
      dimension ablock(5,5),avec(5),bvec(5)

c---------------------------------------------------------------------
c            rhs(i,ic,jc,kc) = rhs(i,ic,jc,kc) 
c     $           - lhs(i,1,ablock,ia)*
c---------------------------------------------------------------------
         bvec(1) = bvec(1) - ablock(1,1)*avec(1)
     >                     - ablock(1,2)*avec(2)
     >                     - ablock(1,3)*avec(3)
     >                     - ablock(1,4)*avec(4)
     >                     - ablock(1,5)*avec(5)
         bvec(2) = bvec(2) - ablock(2,1)*avec(1)
     >                     - ablock(2,2)*avec(2)
     >                     - ablock(2,3)*avec(3)
     >                     - ablock(2,4)*avec(4)
     >                     - ablock(2,5)*avec(5)
         bvec(3) = bvec(3) - ablock(3,1)*avec(1)
     >                     - ablock(3,2)*avec(2)
     >                     - ablock(3,3)*avec(3)
     >                     - ablock(3,4)*avec(4)
     >                     - ablock(3,5)*avec(5)
         bvec(4) = bvec(4) - ablock(4,1)*avec(1)
     >                     - ablock(4,2)*avec(2)
     >                     - ablock(4,3)*avec(3)
     >                     - ablock(4,4)*avec(4)
     >                     - ablock(4,5)*avec(5)
         bvec(5) = bvec(5) - ablock(5,1)*avec(1)
     >                     - ablock(5,2)*avec(2)
     >                     - ablock(5,3)*avec(3)
     >                     - ablock(5,4)*avec(4)
     >                     - ablock(5,5)*avec(5)


      return
      end

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine matmul_sub(ablock, bblock, cblock)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     subtracts a(i,j,k) X b(i,j,k) from c(i,j,k)
c---------------------------------------------------------------------

      implicit none

      double precision ablock, bblock, cblock
      dimension ablock(5,5), bblock(5,5), cblock(5,5)


         cblock(1,1) = cblock(1,1) - ablock(1,1)*bblock(1,1)
     >                             - ablock(1,2)*bblock(2,1)
     >                             - ablock(1,3)*bblock(3,1)
     >                             - ablock(1,4)*bblock(4,1)
     >                             - ablock(1,5)*bblock(5,1)
         cblock(2,1) = cblock(2,1) - ablock(2,1)*bblock(1,1)
     >                             - ablock(2,2)*bblock(2,1)
     >                             - ablock(2,3)*bblock(3,1)
     >                             - ablock(2,4)*bblock(4,1)
     >                             - ablock(2,5)*bblock(5,1)
         cblock(3,1) = cblock(3,1) - ablock(3,1)*bblock(1,1)
     >                             - ablock(3,2)*bblock(2,1)
     >                             - ablock(3,3)*bblock(3,1)
     >                             - ablock(3,4)*bblock(4,1)
     >                             - ablock(3,5)*bblock(5,1)
         cblock(4,1) = cblock(4,1) - ablock(4,1)*bblock(1,1)
     >                             - ablock(4,2)*bblock(2,1)
     >                             - ablock(4,3)*bblock(3,1)
     >                             - ablock(4,4)*bblock(4,1)
     >                             - ablock(4,5)*bblock(5,1)
         cblock(5,1) = cblock(5,1) - ablock(5,1)*bblock(1,1)
     >                             - ablock(5,2)*bblock(2,1)
     >                             - ablock(5,3)*bblock(3,1)
     >                             - ablock(5,4)*bblock(4,1)
     >                             - ablock(5,5)*bblock(5,1)
         cblock(1,2) = cblock(1,2) - ablock(1,1)*bblock(1,2)
     >                             - ablock(1,2)*bblock(2,2)
     >                             - ablock(1,3)*bblock(3,2)
     >                             - ablock(1,4)*bblock(4,2)
     >                             - ablock(1,5)*bblock(5,2)
         cblock(2,2) = cblock(2,2) - ablock(2,1)*bblock(1,2)
     >                             - ablock(2,2)*bblock(2,2)
     >                             - ablock(2,3)*bblock(3,2)
     >                             - ablock(2,4)*bblock(4,2)
     >                             - ablock(2,5)*bblock(5,2)
         cblock(3,2) = cblock(3,2) - ablock(3,1)*bblock(1,2)
     >                             - ablock(3,2)*bblock(2,2)
     >                             - ablock(3,3)*bblock(3,2)
     >                             - ablock(3,4)*bblock(4,2)
     >                             - ablock(3,5)*bblock(5,2)
         cblock(4,2) = cblock(4,2) - ablock(4,1)*bblock(1,2)
     >                             - ablock(4,2)*bblock(2,2)
     >                             - ablock(4,3)*bblock(3,2)
     >                             - ablock(4,4)*bblock(4,2)
     >                             - ablock(4,5)*bblock(5,2)
         cblock(5,2) = cblock(5,2) - ablock(5,1)*bblock(1,2)
     >                             - ablock(5,2)*bblock(2,2)
     >                             - ablock(5,3)*bblock(3,2)
     >                             - ablock(5,4)*bblock(4,2)
     >                             - ablock(5,5)*bblock(5,2)
         cblock(1,3) = cblock(1,3) - ablock(1,1)*bblock(1,3)
     >                             - ablock(1,2)*bblock(2,3)
     >                             - ablock(1,3)*bblock(3,3)
     >                             - ablock(1,4)*bblock(4,3)
     >                             - ablock(1,5)*bblock(5,3)
         cblock(2,3) = cblock(2,3) - ablock(2,1)*bblock(1,3)
     >                             - ablock(2,2)*bblock(2,3)
     >                             - ablock(2,3)*bblock(3,3)
     >                             - ablock(2,4)*bblock(4,3)
     >                             - ablock(2,5)*bblock(5,3)
         cblock(3,3) = cblock(3,3) - ablock(3,1)*bblock(1,3)
     >                             - ablock(3,2)*bblock(2,3)
     >                             - ablock(3,3)*bblock(3,3)
     >                             - ablock(3,4)*bblock(4,3)
     >                             - ablock(3,5)*bblock(5,3)
         cblock(4,3) = cblock(4,3) - ablock(4,1)*bblock(1,3)
     >                             - ablock(4,2)*bblock(2,3)
     >                             - ablock(4,3)*bblock(3,3)
     >                             - ablock(4,4)*bblock(4,3)
     >                             - ablock(4,5)*bblock(5,3)
         cblock(5,3) = cblock(5,3) - ablock(5,1)*bblock(1,3)
     >                             - ablock(5,2)*bblock(2,3)
     >                             - ablock(5,3)*bblock(3,3)
     >                             - ablock(5,4)*bblock(4,3)
     >                             - ablock(5,5)*bblock(5,3)
         cblock(1,4) = cblock(1,4) - ablock(1,1)*bblock(1,4)
     >                             - ablock(1,2)*bblock(2,4)
     >                             - ablock(1,3)*bblock(3,4)
     >                             - ablock(1,4)*bblock(4,4)
     >                             - ablock(1,5)*bblock(5,4)
         cblock(2,4) = cblock(2,4) - ablock(2,1)*bblock(1,4)
     >                             - ablock(2,2)*bblock(2,4)
     >                             - ablock(2,3)*bblock(3,4)
     >                             - ablock(2,4)*bblock(4,4)
     >                             - ablock(2,5)*bblock(5,4)
         cblock(3,4) = cblock(3,4) - ablock(3,1)*bblock(1,4)
     >                             - ablock(3,2)*bblock(2,4)
     >                             - ablock(3,3)*bblock(3,4)
     >                             - ablock(3,4)*bblock(4,4)
     >                             - ablock(3,5)*bblock(5,4)
         cblock(4,4) = cblock(4,4) - ablock(4,1)*bblock(1,4)
     >                             - ablock(4,2)*bblock(2,4)
     >                             - ablock(4,3)*bblock(3,4)
     >                             - ablock(4,4)*bblock(4,4)
     >                             - ablock(4,5)*bblock(5,4)
         cblock(5,4) = cblock(5,4) - ablock(5,1)*bblock(1,4)
     >                             - ablock(5,2)*bblock(2,4)
     >                             - ablock(5,3)*bblock(3,4)
     >                             - ablock(5,4)*bblock(4,4)
     >                             - ablock(5,5)*bblock(5,4)
         cblock(1,5) = cblock(1,5) - ablock(1,1)*bblock(1,5)
     >                             - ablock(1,2)*bblock(2,5)
     >                             - ablock(1,3)*bblock(3,5)
     >                             - ablock(1,4)*bblock(4,5)
     >                             - ablock(1,5)*bblock(5,5)
         cblock(2,5) = cblock(2,5) - ablock(2,1)*bblock(1,5)
     >                             - ablock(2,2)*bblock(2,5)
     >                             - ablock(2,3)*bblock(3,5)
     >                             - ablock(2,4)*bblock(4,5)
     >                             - ablock(2,5)*bblock(5,5)
         cblock(3,5) = cblock(3,5) - ablock(3,1)*bblock(1,5)
     >                             - ablock(3,2)*bblock(2,5)
     >                             - ablock(3,3)*bblock(3,5)
     >                             - ablock(3,4)*bblock(4,5)
     >                             - ablock(3,5)*bblock(5,5)
         cblock(4,5) = cblock(4,5) - ablock(4,1)*bblock(1,5)
     >                             - ablock(4,2)*bblock(2,5)
     >                             - ablock(4,3)*bblock(3,5)
     >                             - ablock(4,4)*bblock(4,5)
     >                             - ablock(4,5)*bblock(5,5)
         cblock(5,5) = cblock(5,5) - ablock(5,1)*bblock(1,5)
     >                             - ablock(5,2)*bblock(2,5)
     >                             - ablock(5,3)*bblock(3,5)
     >                             - ablock(5,4)*bblock(4,5)
     >                             - ablock(5,5)*bblock(5,5)

              
      return
      end



c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine binvcrhs( lhs,c,r )

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     
c---------------------------------------------------------------------

      implicit none

      double precision pivot, coeff, lhs
      dimension lhs(5,5)
      double precision c(5,5), r(5)

c---------------------------------------------------------------------
c     
c---------------------------------------------------------------------

      pivot = 1.00d0/lhs(1,1)
      lhs(1,2) = lhs(1,2)*pivot
      lhs(1,3) = lhs(1,3)*pivot
      lhs(1,4) = lhs(1,4)*pivot
      lhs(1,5) = lhs(1,5)*pivot
      c(1,1) = c(1,1)*pivot
      c(1,2) = c(1,2)*pivot
      c(1,3) = c(1,3)*pivot
      c(1,4) = c(1,4)*pivot
      c(1,5) = c(1,5)*pivot
      r(1)   = r(1)  *pivot

      coeff = lhs(2,1)
      lhs(2,2)= lhs(2,2) - coeff*lhs(1,2)
      lhs(2,3)= lhs(2,3) - coeff*lhs(1,3)
      lhs(2,4)= lhs(2,4) - coeff*lhs(1,4)
      lhs(2,5)= lhs(2,5) - coeff*lhs(1,5)
      c(2,1) = c(2,1) - coeff*c(1,1)
      c(2,2) = c(2,2) - coeff*c(1,2)
      c(2,3) = c(2,3) - coeff*c(1,3)
      c(2,4) = c(2,4) - coeff*c(1,4)
      c(2,5) = c(2,5) - coeff*c(1,5)
      r(2)   = r(2)   - coeff*r(1)

      coeff = lhs(3,1)
      lhs(3,2)= lhs(3,2) - coeff*lhs(1,2)
      lhs(3,3)= lhs(3,3) - coeff*lhs(1,3)
      lhs(3,4)= lhs(3,4) - coeff*lhs(1,4)
      lhs(3,5)= lhs(3,5) - coeff*lhs(1,5)
      c(3,1) = c(3,1) - coeff*c(1,1)
      c(3,2) = c(3,2) - coeff*c(1,2)
      c(3,3) = c(3,3) - coeff*c(1,3)
      c(3,4) = c(3,4) - coeff*c(1,4)
      c(3,5) = c(3,5) - coeff*c(1,5)
      r(3)   = r(3)   - coeff*r(1)

      coeff = lhs(4,1)
      lhs(4,2)= lhs(4,2) - coeff*lhs(1,2)
      lhs(4,3)= lhs(4,3) - coeff*lhs(1,3)
      lhs(4,4)= lhs(4,4) - coeff*lhs(1,4)
      lhs(4,5)= lhs(4,5) - coeff*lhs(1,5)
      c(4,1) = c(4,1) - coeff*c(1,1)
      c(4,2) = c(4,2) - coeff*c(1,2)
      c(4,3) = c(4,3) - coeff*c(1,3)
      c(4,4) = c(4,4) - coeff*c(1,4)
      c(4,5) = c(4,5) - coeff*c(1,5)
      r(4)   = r(4)   - coeff*r(1)

      coeff = lhs(5,1)
      lhs(5,2)= lhs(5,2) - coeff*lhs(1,2)
      lhs(5,3)= lhs(5,3) - coeff*lhs(1,3)
      lhs(5,4)= lhs(5,4) - coeff*lhs(1,4)
      lhs(5,5)= lhs(5,5) - coeff*lhs(1,5)
      c(5,1) = c(5,1) - coeff*c(1,1)
      c(5,2) = c(5,2) - coeff*c(1,2)
      c(5,3) = c(5,3) - coeff*c(1,3)
      c(5,4) = c(5,4) - coeff*c(1,4)
      c(5,5) = c(5,5) - coeff*c(1,5)
      r(5)   = r(5)   - coeff*r(1)


      pivot = 1.00d0/lhs(2,2)
      lhs(2,3) = lhs(2,3)*pivot
      lhs(2,4) = lhs(2,4)*pivot
      lhs(2,5) = lhs(2,5)*pivot
      c(2,1) = c(2,1)*pivot
      c(2,2) = c(2,2)*pivot
      c(2,3) = c(2,3)*pivot
      c(2,4) = c(2,4)*pivot
      c(2,5) = c(2,5)*pivot
      r(2)   = r(2)  *pivot

      coeff = lhs(1,2)
      lhs(1,3)= lhs(1,3) - coeff*lhs(2,3)
      lhs(1,4)= lhs(1,4) - coeff*lhs(2,4)
      lhs(1,5)= lhs(1,5) - coeff*lhs(2,5)
      c(1,1) = c(1,1) - coeff*c(2,1)
      c(1,2) = c(1,2) - coeff*c(2,2)
      c(1,3) = c(1,3) - coeff*c(2,3)
      c(1,4) = c(1,4) - coeff*c(2,4)
      c(1,5) = c(1,5) - coeff*c(2,5)
      r(1)   = r(1)   - coeff*r(2)

      coeff = lhs(3,2)
      lhs(3,3)= lhs(3,3) - coeff*lhs(2,3)
      lhs(3,4)= lhs(3,4) - coeff*lhs(2,4)
      lhs(3,5)= lhs(3,5) - coeff*lhs(2,5)
      c(3,1) = c(3,1) - coeff*c(2,1)
      c(3,2) = c(3,2) - coeff*c(2,2)
      c(3,3) = c(3,3) - coeff*c(2,3)
      c(3,4) = c(3,4) - coeff*c(2,4)
      c(3,5) = c(3,5) - coeff*c(2,5)
      r(3)   = r(3)   - coeff*r(2)

      coeff = lhs(4,2)
      lhs(4,3)= lhs(4,3) - coeff*lhs(2,3)
      lhs(4,4)= lhs(4,4) - coeff*lhs(2,4)
      lhs(4,5)= lhs(4,5) - coeff*lhs(2,5)
      c(4,1) = c(4,1) - coeff*c(2,1)
      c(4,2) = c(4,2) - coeff*c(2,2)
      c(4,3) = c(4,3) - coeff*c(2,3)
      c(4,4) = c(4,4) - coeff*c(2,4)
      c(4,5) = c(4,5) - coeff*c(2,5)
      r(4)   = r(4)   - coeff*r(2)

      coeff = lhs(5,2)
      lhs(5,3)= lhs(5,3) - coeff*lhs(2,3)
      lhs(5,4)= lhs(5,4) - coeff*lhs(2,4)
      lhs(5,5)= lhs(5,5) - coeff*lhs(2,5)
      c(5,1) = c(5,1) - coeff*c(2,1)
      c(5,2) = c(5,2) - coeff*c(2,2)
      c(5,3) = c(5,3) - coeff*c(2,3)
      c(5,4) = c(5,4) - coeff*c(2,4)
      c(5,5) = c(5,5) - coeff*c(2,5)
      r(5)   = r(5)   - coeff*r(2)


      pivot = 1.00d0/lhs(3,3)
      lhs(3,4) = lhs(3,4)*pivot
      lhs(3,5) = lhs(3,5)*pivot
      c(3,1) = c(3,1)*pivot
      c(3,2) = c(3,2)*pivot
      c(3,3) = c(3,3)*pivot
      c(3,4) = c(3,4)*pivot
      c(3,5) = c(3,5)*pivot
      r(3)   = r(3)  *pivot

      coeff = lhs(1,3)
      lhs(1,4)= lhs(1,4) - coeff*lhs(3,4)
      lhs(1,5)= lhs(1,5) - coeff*lhs(3,5)
      c(1,1) = c(1,1) - coeff*c(3,1)
      c(1,2) = c(1,2) - coeff*c(3,2)
      c(1,3) = c(1,3) - coeff*c(3,3)
      c(1,4) = c(1,4) - coeff*c(3,4)
      c(1,5) = c(1,5) - coeff*c(3,5)
      r(1)   = r(1)   - coeff*r(3)

      coeff = lhs(2,3)
      lhs(2,4)= lhs(2,4) - coeff*lhs(3,4)
      lhs(2,5)= lhs(2,5) - coeff*lhs(3,5)
      c(2,1) = c(2,1) - coeff*c(3,1)
      c(2,2) = c(2,2) - coeff*c(3,2)
      c(2,3) = c(2,3) - coeff*c(3,3)
      c(2,4) = c(2,4) - coeff*c(3,4)
      c(2,5) = c(2,5) - coeff*c(3,5)
      r(2)   = r(2)   - coeff*r(3)

      coeff = lhs(4,3)
      lhs(4,4)= lhs(4,4) - coeff*lhs(3,4)
      lhs(4,5)= lhs(4,5) - coeff*lhs(3,5)
      c(4,1) = c(4,1) - coeff*c(3,1)
      c(4,2) = c(4,2) - coeff*c(3,2)
      c(4,3) = c(4,3) - coeff*c(3,3)
      c(4,4) = c(4,4) - coeff*c(3,4)
      c(4,5) = c(4,5) - coeff*c(3,5)
      r(4)   = r(4)   - coeff*r(3)

      coeff = lhs(5,3)
      lhs(5,4)= lhs(5,4) - coeff*lhs(3,4)
      lhs(5,5)= lhs(5,5) - coeff*lhs(3,5)
      c(5,1) = c(5,1) - coeff*c(3,1)
      c(5,2) = c(5,2) - coeff*c(3,2)
      c(5,3) = c(5,3) - coeff*c(3,3)
      c(5,4) = c(5,4) - coeff*c(3,4)
      c(5,5) = c(5,5) - coeff*c(3,5)
      r(5)   = r(5)   - coeff*r(3)


      pivot = 1.00d0/lhs(4,4)
      lhs(4,5) = lhs(4,5)*pivot
      c(4,1) = c(4,1)*pivot
      c(4,2) = c(4,2)*pivot
      c(4,3) = c(4,3)*pivot
      c(4,4) = c(4,4)*pivot
      c(4,5) = c(4,5)*pivot
      r(4)   = r(4)  *pivot

      coeff = lhs(1,4)
      lhs(1,5)= lhs(1,5) - coeff*lhs(4,5)
      c(1,1) = c(1,1) - coeff*c(4,1)
      c(1,2) = c(1,2) - coeff*c(4,2)
      c(1,3) = c(1,3) - coeff*c(4,3)
      c(1,4) = c(1,4) - coeff*c(4,4)
      c(1,5) = c(1,5) - coeff*c(4,5)
      r(1)   = r(1)   - coeff*r(4)

      coeff = lhs(2,4)
      lhs(2,5)= lhs(2,5) - coeff*lhs(4,5)
      c(2,1) = c(2,1) - coeff*c(4,1)
      c(2,2) = c(2,2) - coeff*c(4,2)
      c(2,3) = c(2,3) - coeff*c(4,3)
      c(2,4) = c(2,4) - coeff*c(4,4)
      c(2,5) = c(2,5) - coeff*c(4,5)
      r(2)   = r(2)   - coeff*r(4)

      coeff = lhs(3,4)
      lhs(3,5)= lhs(3,5) - coeff*lhs(4,5)
      c(3,1) = c(3,1) - coeff*c(4,1)
      c(3,2) = c(3,2) - coeff*c(4,2)
      c(3,3) = c(3,3) - coeff*c(4,3)
      c(3,4) = c(3,4) - coeff*c(4,4)
      c(3,5) = c(3,5) - coeff*c(4,5)
      r(3)   = r(3)   - coeff*r(4)

      coeff = lhs(5,4)
      lhs(5,5)= lhs(5,5) - coeff*lhs(4,5)
      c(5,1) = c(5,1) - coeff*c(4,1)
      c(5,2) = c(5,2) - coeff*c(4,2)
      c(5,3) = c(5,3) - coeff*c(4,3)
      c(5,4) = c(5,4) - coeff*c(4,4)
      c(5,5) = c(5,5) - coeff*c(4,5)
      r(5)   = r(5)   - coeff*r(4)


      pivot = 1.00d0/lhs(5,5)
      c(5,1) = c(5,1)*pivot
      c(5,2) = c(5,2)*pivot
      c(5,3) = c(5,3)*pivot
      c(5,4) = c(5,4)*pivot
      c(5,5) = c(5,5)*pivot
      r(5)   = r(5)  *pivot

      coeff = lhs(1,5)
      c(1,1) = c(1,1) - coeff*c(5,1)
      c(1,2) = c(1,2) - coeff*c(5,2)
      c(1,3) = c(1,3) - coeff*c(5,3)
      c(1,4) = c(1,4) - coeff*c(5,4)
      c(1,5) = c(1,5) - coeff*c(5,5)
      r(1)   = r(1)   - coeff*r(5)

      coeff = lhs(2,5)
      c(2,1) = c(2,1) - coeff*c(5,1)
      c(2,2) = c(2,2) - coeff*c(5,2)
      c(2,3) = c(2,3) - coeff*c(5,3)
      c(2,4) = c(2,4) - coeff*c(5,4)
      c(2,5) = c(2,5) - coeff*c(5,5)
      r(2)   = r(2)   - coeff*r(5)

      coeff = lhs(3,5)
      c(3,1) = c(3,1) - coeff*c(5,1)
      c(3,2) = c(3,2) - coeff*c(5,2)
      c(3,3) = c(3,3) - coeff*c(5,3)
      c(3,4) = c(3,4) - coeff*c(5,4)
      c(3,5) = c(3,5) - coeff*c(5,5)
      r(3)   = r(3)   - coeff*r(5)

      coeff = lhs(4,5)
      c(4,1) = c(4,1) - coeff*c(5,1)
      c(4,2) = c(4,2) - coeff*c(5,2)
      c(4,3) = c(4,3) - coeff*c(5,3)
      c(4,4) = c(4,4) - coeff*c(5,4)
      c(4,5) = c(4,5) - coeff*c(5,5)
      r(4)   = r(4)   - coeff*r(5)


      return
      end



c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine binvrhs( lhs,r )

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     
c---------------------------------------------------------------------

      implicit none

      double precision pivot, coeff, lhs
      dimension lhs(5,5)
      double precision r(5)

c---------------------------------------------------------------------
c     
c---------------------------------------------------------------------


      pivot = 1.00d0/lhs(1,1)
      lhs(1,2) = lhs(1,2)*pivot
      lhs(1,3) = lhs(1,3)*pivot
      lhs(1,4) = lhs(1,4)*pivot
      lhs(1,5) = lhs(1,5)*pivot
      r(1)   = r(1)  *pivot

      coeff = lhs(2,1)
      lhs(2,2)= lhs(2,2) - coeff*lhs(1,2)
      lhs(2,3)= lhs(2,3) - coeff*lhs(1,3)
      lhs(2,4)= lhs(2,4) - coeff*lhs(1,4)
      lhs(2,5)= lhs(2,5) - coeff*lhs(1,5)
      r(2)   = r(2)   - coeff*r(1)

      coeff = lhs(3,1)
      lhs(3,2)= lhs(3,2) - coeff*lhs(1,2)
      lhs(3,3)= lhs(3,3) - coeff*lhs(1,3)
      lhs(3,4)= lhs(3,4) - coeff*lhs(1,4)
      lhs(3,5)= lhs(3,5) - coeff*lhs(1,5)
      r(3)   = r(3)   - coeff*r(1)

      coeff = lhs(4,1)
      lhs(4,2)= lhs(4,2) - coeff*lhs(1,2)
      lhs(4,3)= lhs(4,3) - coeff*lhs(1,3)
      lhs(4,4)= lhs(4,4) - coeff*lhs(1,4)
      lhs(4,5)= lhs(4,5) - coeff*lhs(1,5)
      r(4)   = r(4)   - coeff*r(1)

      coeff = lhs(5,1)
      lhs(5,2)= lhs(5,2) - coeff*lhs(1,2)
      lhs(5,3)= lhs(5,3) - coeff*lhs(1,3)
      lhs(5,4)= lhs(5,4) - coeff*lhs(1,4)
      lhs(5,5)= lhs(5,5) - coeff*lhs(1,5)
      r(5)   = r(5)   - coeff*r(1)


      pivot = 1.00d0/lhs(2,2)
      lhs(2,3) = lhs(2,3)*pivot
      lhs(2,4) = lhs(2,4)*pivot
      lhs(2,5) = lhs(2,5)*pivot
      r(2)   = r(2)  *pivot

      coeff = lhs(1,2)
      lhs(1,3)= lhs(1,3) - coeff*lhs(2,3)
      lhs(1,4)= lhs(1,4) - coeff*lhs(2,4)
      lhs(1,5)= lhs(1,5) - coeff*lhs(2,5)
      r(1)   = r(1)   - coeff*r(2)

      coeff = lhs(3,2)
      lhs(3,3)= lhs(3,3) - coeff*lhs(2,3)
      lhs(3,4)= lhs(3,4) - coeff*lhs(2,4)
      lhs(3,5)= lhs(3,5) - coeff*lhs(2,5)
      r(3)   = r(3)   - coeff*r(2)

      coeff = lhs(4,2)
      lhs(4,3)= lhs(4,3) - coeff*lhs(2,3)
      lhs(4,4)= lhs(4,4) - coeff*lhs(2,4)
      lhs(4,5)= lhs(4,5) - coeff*lhs(2,5)
      r(4)   = r(4)   - coeff*r(2)

      coeff = lhs(5,2)
      lhs(5,3)= lhs(5,3) - coeff*lhs(2,3)
      lhs(5,4)= lhs(5,4) - coeff*lhs(2,4)
      lhs(5,5)= lhs(5,5) - coeff*lhs(2,5)
      r(5)   = r(5)   - coeff*r(2)


      pivot = 1.00d0/lhs(3,3)
      lhs(3,4) = lhs(3,4)*pivot
      lhs(3,5) = lhs(3,5)*pivot
      r(3)   = r(3)  *pivot

      coeff = lhs(1,3)
      lhs(1,4)= lhs(1,4) - coeff*lhs(3,4)
      lhs(1,5)= lhs(1,5) - coeff*lhs(3,5)
      r(1)   = r(1)   - coeff*r(3)

      coeff = lhs(2,3)
      lhs(2,4)= lhs(2,4) - coeff*lhs(3,4)
      lhs(2,5)= lhs(2,5) - coeff*lhs(3,5)
      r(2)   = r(2)   - coeff*r(3)

      coeff = lhs(4,3)
      lhs(4,4)= lhs(4,4) - coeff*lhs(3,4)
      lhs(4,5)= lhs(4,5) - coeff*lhs(3,5)
      r(4)   = r(4)   - coeff*r(3)

      coeff = lhs(5,3)
      lhs(5,4)= lhs(5,4) - coeff*lhs(3,4)
      lhs(5,5)= lhs(5,5) - coeff*lhs(3,5)
      r(5)   = r(5)   - coeff*r(3)


      pivot = 1.00d0/lhs(4,4)
      lhs(4,5) = lhs(4,5)*pivot
      r(4)   = r(4)  *pivot

      coeff = lhs(1,4)
      lhs(1,5)= lhs(1,5) - coeff*lhs(4,5)
      r(1)   = r(1)   - coeff*r(4)

      coeff = lhs(2,4)
      lhs(2,5)= lhs(2,5) - coeff*lhs(4,5)
      r(2)   = r(2)   - coeff*r(4)

      coeff = lhs(3,4)
      lhs(3,5)= lhs(3,5) - coeff*lhs(4,5)
      r(3)   = r(3)   - coeff*r(4)

      coeff = lhs(5,4)
      lhs(5,5)= lhs(5,5) - coeff*lhs(4,5)
      r(5)   = r(5)   - coeff*r(4)


      pivot = 1.00d0/lhs(5,5)
      r(5)   = r(5)  *pivot

      coeff = lhs(1,5)
      r(1)   = r(1)   - coeff*r(5)

      coeff = lhs(2,5)
      r(2)   = r(2)   - coeff*r(5)

      coeff = lhs(3,5)
      r(3)   = r(3)   - coeff*r(5)

      coeff = lhs(4,5)
      r(4)   = r(4)   - coeff*r(5)


      return
      end



c---------------------------------------------------------------------
c---------------------------------------------------------------------
      
      subroutine timer_clear(n)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      implicit none
      integer n
      
      double precision start(64), elapsed(64)
      common /tt/ start, elapsed

      elapsed(n) = 0.0
      return
      end


c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine timer_start(n)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      implicit none
      external         elapsed_time
      double precision elapsed_time
      integer n
      double precision start(64), elapsed(64)
      common /tt/ start, elapsed

      start(n) = elapsed_time()

      return
      end
      

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine timer_stop(n)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

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


c---------------------------------------------------------------------
c---------------------------------------------------------------------

      double precision function timer_read(n)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      implicit none
      integer n
      double precision start(64), elapsed(64)
      common /tt/ start, elapsed
      
      timer_read = elapsed(n)
      return
      end


c---------------------------------------------------------------------
c---------------------------------------------------------------------

      double precision function elapsed_time()

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      implicit none

      double precision t

c This function must measure wall clock time, not CPU time. 
c Since there is no portable timer in Fortran (77)
c we call a routine compiled in C (though the C source may have
c to be tweaked). 
      call wtime(t)
c The following is not ok for "official" results because it reports
c CPU time not wall clock time. It may be useful for developing/testing
c on timeshared Crays, though. 
c     call second(t)

      elapsed_time = t

      return
      end


c---------------------------------------------------------------------
c---------------------------------------------------------------------

        subroutine verify(no_time_steps, class, verified)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c  verification routine                         
c---------------------------------------------------------------------

        include 'header.h'

        double precision xcrref(5),xceref(5),xcrdif(5),xcedif(5), 
     >                   epsilon, xce(5), xcr(5), dtref
        integer m, no_time_steps
        character class
        logical verified

c---------------------------------------------------------------------
c   tolerance level
c---------------------------------------------------------------------
        epsilon = 1.0d-08


c---------------------------------------------------------------------
c   compute the error norm and the residual norm, and exit if not printing
c---------------------------------------------------------------------
        call error_norm(xce)
        call compute_rhs

        call rhs_norm(xcr)

        do m = 1, 5
           xcr(m) = xcr(m) / dt
        enddo


        class = 'U'
        verified = .true.

        do m = 1,5
           xcrref(m) = 1.0
           xceref(m) = 1.0
        end do

c---------------------------------------------------------------------
c    reference data for 12X12X12 grids after 60 time steps, with DT = 1.0d-02
c---------------------------------------------------------------------
        if ( (grid_points(1)  .eq. 12     ) .and. 
     >       (grid_points(2)  .eq. 12     ) .and.
     >       (grid_points(3)  .eq. 12     ) .and.
     >       (no_time_steps   .eq. 60    ))  then

           class = 'S'
           dtref = 1.0d-2

c---------------------------------------------------------------------
c  Reference values of RMS-norms of residual.
c---------------------------------------------------------------------
         xcrref(1) = 1.7034283709541311d-01
         xcrref(2) = 1.2975252070034097d-02
         xcrref(3) = 3.2527926989486055d-02
         xcrref(4) = 2.6436421275166801d-02
         xcrref(5) = 1.9211784131744430d-01

c---------------------------------------------------------------------
c  Reference values of RMS-norms of solution error.
c---------------------------------------------------------------------
         xceref(1) = 4.9976913345811579d-04
         xceref(2) = 4.5195666782961927d-05
         xceref(3) = 7.3973765172921357d-05
         xceref(4) = 7.3821238632439731d-05
         xceref(5) = 8.9269630987491446d-04

c---------------------------------------------------------------------
c    reference data for 24X24X24 grids after 200 time steps, with DT = 0.8d-3
c---------------------------------------------------------------------
        elseif ( (grid_points(1) .eq. 24) .and. 
     >           (grid_points(2) .eq. 24) .and.
     >           (grid_points(3) .eq. 24) .and.
     >           (no_time_steps .eq. 200) ) then

           class = 'W'
           dtref = 0.8d-3
c---------------------------------------------------------------------
c  Reference values of RMS-norms of residual.
c---------------------------------------------------------------------
           xcrref(1) = 0.1125590409344d+03
           xcrref(2) = 0.1180007595731d+02
           xcrref(3) = 0.2710329767846d+02
           xcrref(4) = 0.2469174937669d+02
           xcrref(5) = 0.2638427874317d+03

c---------------------------------------------------------------------
c  Reference values of RMS-norms of solution error.
c---------------------------------------------------------------------
           xceref(1) = 0.4419655736008d+01
           xceref(2) = 0.4638531260002d+00
           xceref(3) = 0.1011551749967d+01
           xceref(4) = 0.9235878729944d+00
           xceref(5) = 0.1018045837718d+02


c---------------------------------------------------------------------
c    reference data for 64X64X64 grids after 200 time steps, with DT = 0.8d-3
c---------------------------------------------------------------------
        elseif ( (grid_points(1) .eq. 64) .and. 
     >           (grid_points(2) .eq. 64) .and.
     >           (grid_points(3) .eq. 64) .and.
     >           (no_time_steps .eq. 200) ) then

           class = 'A'
           dtref = 0.8d-3
c---------------------------------------------------------------------
c  Reference values of RMS-norms of residual.
c---------------------------------------------------------------------
         xcrref(1) = 1.0806346714637264d+02
         xcrref(2) = 1.1319730901220813d+01
         xcrref(3) = 2.5974354511582465d+01
         xcrref(4) = 2.3665622544678910d+01
         xcrref(5) = 2.5278963211748344d+02

c---------------------------------------------------------------------
c  Reference values of RMS-norms of solution error.
c---------------------------------------------------------------------
         xceref(1) = 4.2348416040525025d+00
         xceref(2) = 4.4390282496995698d-01
         xceref(3) = 9.6692480136345650d-01
         xceref(4) = 8.8302063039765474d-01
         xceref(5) = 9.7379901770829278d+00

c---------------------------------------------------------------------
c    reference data for 102X102X102 grids after 200 time steps,
c    with DT = 3.0d-04
c---------------------------------------------------------------------
        elseif ( (grid_points(1) .eq. 102) .and. 
     >           (grid_points(2) .eq. 102) .and.
     >           (grid_points(3) .eq. 102) .and.
     >           (no_time_steps .eq. 200) ) then

           class = 'B'
           dtref = 3.0d-4

c---------------------------------------------------------------------
c  Reference values of RMS-norms of residual.
c---------------------------------------------------------------------
         xcrref(1) = 1.4233597229287254d+03
         xcrref(2) = 9.9330522590150238d+01
         xcrref(3) = 3.5646025644535285d+02
         xcrref(4) = 3.2485447959084092d+02
         xcrref(5) = 3.2707541254659363d+03

c---------------------------------------------------------------------
c  Reference values of RMS-norms of solution error.
c---------------------------------------------------------------------
         xceref(1) = 5.2969847140936856d+01
         xceref(2) = 4.4632896115670668d+00
         xceref(3) = 1.3122573342210174d+01
         xceref(4) = 1.2006925323559144d+01
         xceref(5) = 1.2459576151035986d+02

c---------------------------------------------------------------------
c    reference data for 162X162X162 grids after 200 time steps,
c    with DT = 1.0d-04
c---------------------------------------------------------------------
        elseif ( (grid_points(1) .eq. 162) .and. 
     >           (grid_points(2) .eq. 162) .and.
     >           (grid_points(3) .eq. 162) .and.
     >           (no_time_steps .eq. 200) ) then

           class = 'C'
           dtref = 1.0d-4

c---------------------------------------------------------------------
c  Reference values of RMS-norms of residual.
c---------------------------------------------------------------------
         xcrref(1) = 0.62398116551764615d+04
         xcrref(2) = 0.50793239190423964d+03
         xcrref(3) = 0.15423530093013596d+04
         xcrref(4) = 0.13302387929291190d+04
         xcrref(5) = 0.11604087428436455d+05

c---------------------------------------------------------------------
c  Reference values of RMS-norms of solution error.
c---------------------------------------------------------------------
         xceref(1) = 0.16462008369091265d+03
         xceref(2) = 0.11497107903824313d+02
         xceref(3) = 0.41207446207461508d+02
         xceref(4) = 0.37087651059694167d+02
         xceref(5) = 0.36211053051841265d+03

c---------------------------------------------------------------------
c    reference data for 408x408x408 grids after 250 time steps,
c    with DT = 0.2d-04
c---------------------------------------------------------------------
        elseif ( (grid_points(1) .eq. 408) .and. 
     >           (grid_points(2) .eq. 408) .and.
     >           (grid_points(3) .eq. 408) .and.
     >           (no_time_steps .eq. 250) ) then

           class = 'D'
           dtref = 0.2d-4

c---------------------------------------------------------------------
c  Reference values of RMS-norms of residual.
c---------------------------------------------------------------------
         xcrref(1) = 0.2533188551738d+05
         xcrref(2) = 0.2346393716980d+04
         xcrref(3) = 0.6294554366904d+04
         xcrref(4) = 0.5352565376030d+04
         xcrref(5) = 0.3905864038618d+05

c---------------------------------------------------------------------
c  Reference values of RMS-norms of solution error.
c---------------------------------------------------------------------

         xceref(1) = 0.3100009377557d+03
         xceref(2) = 0.2424086324913d+02
         xceref(3) = 0.7782212022645d+02
         xceref(4) = 0.6835623860116d+02
         xceref(5) = 0.6065737200368d+03


c---------------------------------------------------------------------
c    reference data for 1020x1020x1020 grids after 250 time steps,
c    with DT = 0.4d-05
c---------------------------------------------------------------------
        elseif ( (grid_points(1) .eq. 1020) .and. 
     >           (grid_points(2) .eq. 1020) .and.
     >           (grid_points(3) .eq. 1020) .and.
     >           (no_time_steps .eq. 250) ) then

           class = 'E'
           dtref = 0.4d-5

c---------------------------------------------------------------------
c  Reference values of RMS-norms of residual.
c---------------------------------------------------------------------
         xcrref(1) = 0.9795372484517d+05
         xcrref(2) = 0.9739814511521d+04
         xcrref(3) = 0.2467606342965d+05
         xcrref(4) = 0.2092419572860d+05
         xcrref(5) = 0.1392138856939d+06

c---------------------------------------------------------------------
c  Reference values of RMS-norms of solution error.
c---------------------------------------------------------------------

         xceref(1) = 0.4327562208414d+03
         xceref(2) = 0.3699051964887d+02
         xceref(3) = 0.1089845040954d+03
         xceref(4) = 0.9462517622043d+02
         xceref(5) = 0.7765512765309d+03


        else
           verified = .false.
        endif

c---------------------------------------------------------------------
c    verification test for residuals if gridsize is one of 
c    the defined grid sizes above (class .ne. 'U')
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c    Compute the difference of solution values and the known reference values.
c---------------------------------------------------------------------
        do m = 1, 5
           
           xcrdif(m) = dabs((xcr(m)-xcrref(m))/xcrref(m)) 
           xcedif(m) = dabs((xce(m)-xceref(m))/xceref(m))
           
        enddo

c---------------------------------------------------------------------
c    Output the comparison of computed results to known cases.
c---------------------------------------------------------------------

        if (class .ne. 'U') then
           write(*, 1990) class
 1990      format(' Verification being performed for class ', a)
           write (*,2000) epsilon
 2000      format(' accuracy setting for epsilon = ', E20.13)
           verified = (dabs(dt-dtref) .le. epsilon)
           if (.not.verified) then  
              class = 'U'
              write (*,1000) dtref
 1000         format(' DT does not match the reference value of ', 
     >                 E15.8)
           endif
        else 
           write(*, 1995)
 1995      format(' Unknown class')
        endif


        if (class .ne. 'U') then
           write (*, 2001) 
        else
           write (*, 2005)
        endif

 2001   format(' Comparison of RMS-norms of residual')
 2005   format(' RMS-norms of residual')
        do m = 1, 5
           if (class .eq. 'U') then
              write(*, 2015) m, xcr(m)
           else if (xcrdif(m) .le. epsilon) then
              write (*,2011) m,xcr(m),xcrref(m),xcrdif(m)
           else 
              verified = .false.
              write (*,2010) m,xcr(m),xcrref(m),xcrdif(m)
           endif
        enddo

        if (class .ne. 'U') then
           write (*,2002)
        else
           write (*,2006)
        endif
 2002   format(' Comparison of RMS-norms of solution error')
 2006   format(' RMS-norms of solution error')
        
        do m = 1, 5
           if (class .eq. 'U') then
              write(*, 2015) m, xce(m)
           else if (xcedif(m) .le. epsilon) then
              write (*,2011) m,xce(m),xceref(m),xcedif(m)
           else
              verified = .false.
              write (*,2010) m,xce(m),xceref(m),xcedif(m)
           endif
        enddo
        
 2010   format(' FAILURE: ', i2, E20.13, E20.13, E20.13)
 2011   format('          ', i2, E20.13, E20.13, E20.13)
 2015   format('          ', i2, E20.13)
        
        if (class .eq. 'U') then
           write(*, 2022)
           write(*, 2023)
 2022      format(' No reference values provided')
 2023      format(' No verification performed')
        else if (verified) then
           write(*, 2020)
 2020      format(' Verification Successful')
        else
           write(*, 2021)
 2021      format(' Verification failed')
        endif

        return


        end

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine x_solve

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     
c     Performs line solves in X direction by first factoring
c     the block-tridiagonal matrix into an upper triangular matrix, 
c     and then performing back substitution to solve for the unknow
c     vectors of each line.  
c     
c     Make sure we treat elements zero to cell_size in the direction
c     of the sweep.
c     
c---------------------------------------------------------------------

      include 'header.h'
      include 'work_lhs.h'

      integer i,j,k,m,n,isize

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      if (timeron) call timer_start(t_xsolve)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     This function computes the left hand side in the xi-direction
c---------------------------------------------------------------------

      isize = grid_points(1)-1

c---------------------------------------------------------------------
c     determine a (labeled f) and n jacobians
c---------------------------------------------------------------------
      do k = 1, grid_points(3)-2
         do j = 1, grid_points(2)-2
            do i = 0, isize

               tmp1 = rho_i(i,j,k)
               tmp2 = tmp1 * tmp1
               tmp3 = tmp1 * tmp2
c---------------------------------------------------------------------
c     
c---------------------------------------------------------------------
               fjac(1,1,i) = 0.0d+00
               fjac(1,2,i) = 1.0d+00
               fjac(1,3,i) = 0.0d+00
               fjac(1,4,i) = 0.0d+00
               fjac(1,5,i) = 0.0d+00

               fjac(2,1,i) = -(u(2,i,j,k) * tmp2 * 
     >              u(2,i,j,k))
     >              + c2 * qs(i,j,k)
               fjac(2,2,i) = ( 2.0d+00 - c2 )
     >              * ( u(2,i,j,k) / u(1,i,j,k) )
               fjac(2,3,i) = - c2 * ( u(3,i,j,k) * tmp1 )
               fjac(2,4,i) = - c2 * ( u(4,i,j,k) * tmp1 )
               fjac(2,5,i) = c2

               fjac(3,1,i) = - ( u(2,i,j,k)*u(3,i,j,k) ) * tmp2
               fjac(3,2,i) = u(3,i,j,k) * tmp1
               fjac(3,3,i) = u(2,i,j,k) * tmp1
               fjac(3,4,i) = 0.0d+00
               fjac(3,5,i) = 0.0d+00

               fjac(4,1,i) = - ( u(2,i,j,k)*u(4,i,j,k) ) * tmp2
               fjac(4,2,i) = u(4,i,j,k) * tmp1
               fjac(4,3,i) = 0.0d+00
               fjac(4,4,i) = u(2,i,j,k) * tmp1
               fjac(4,5,i) = 0.0d+00

               fjac(5,1,i) = ( c2 * 2.0d0 * square(i,j,k)
     >              - c1 * u(5,i,j,k) )
     >              * ( u(2,i,j,k) * tmp2 )
               fjac(5,2,i) = c1 *  u(5,i,j,k) * tmp1 
     >              - c2
     >              * ( u(2,i,j,k)*u(2,i,j,k) * tmp2
     >              + qs(i,j,k) )
               fjac(5,3,i) = - c2 * ( u(3,i,j,k)*u(2,i,j,k) )
     >              * tmp2
               fjac(5,4,i) = - c2 * ( u(4,i,j,k)*u(2,i,j,k) )
     >              * tmp2
               fjac(5,5,i) = c1 * ( u(2,i,j,k) * tmp1 )

               njac(1,1,i) = 0.0d+00
               njac(1,2,i) = 0.0d+00
               njac(1,3,i) = 0.0d+00
               njac(1,4,i) = 0.0d+00
               njac(1,5,i) = 0.0d+00

               njac(2,1,i) = - con43 * c3c4 * tmp2 * u(2,i,j,k)
               njac(2,2,i) =   con43 * c3c4 * tmp1
               njac(2,3,i) =   0.0d+00
               njac(2,4,i) =   0.0d+00
               njac(2,5,i) =   0.0d+00

               njac(3,1,i) = - c3c4 * tmp2 * u(3,i,j,k)
               njac(3,2,i) =   0.0d+00
               njac(3,3,i) =   c3c4 * tmp1
               njac(3,4,i) =   0.0d+00
               njac(3,5,i) =   0.0d+00

               njac(4,1,i) = - c3c4 * tmp2 * u(4,i,j,k)
               njac(4,2,i) =   0.0d+00 
               njac(4,3,i) =   0.0d+00
               njac(4,4,i) =   c3c4 * tmp1
               njac(4,5,i) =   0.0d+00

               njac(5,1,i) = - ( con43 * c3c4
     >              - c1345 ) * tmp3 * (u(2,i,j,k)**2)
     >              - ( c3c4 - c1345 ) * tmp3 * (u(3,i,j,k)**2)
     >              - ( c3c4 - c1345 ) * tmp3 * (u(4,i,j,k)**2)
     >              - c1345 * tmp2 * u(5,i,j,k)

               njac(5,2,i) = ( con43 * c3c4
     >              - c1345 ) * tmp2 * u(2,i,j,k)
               njac(5,3,i) = ( c3c4 - c1345 ) * tmp2 * u(3,i,j,k)
               njac(5,4,i) = ( c3c4 - c1345 ) * tmp2 * u(4,i,j,k)
               njac(5,5,i) = ( c1345 ) * tmp1

            enddo
c---------------------------------------------------------------------
c     now jacobians set, so form left hand side in x direction
c---------------------------------------------------------------------
            call lhsinit(lhs, isize)
            do i = 1, isize-1

               tmp1 = dt * tx1
               tmp2 = dt * tx2

               lhs(1,1,aa,i) = - tmp2 * fjac(1,1,i-1)
     >              - tmp1 * njac(1,1,i-1)
     >              - tmp1 * dx1 
               lhs(1,2,aa,i) = - tmp2 * fjac(1,2,i-1)
     >              - tmp1 * njac(1,2,i-1)
               lhs(1,3,aa,i) = - tmp2 * fjac(1,3,i-1)
     >              - tmp1 * njac(1,3,i-1)
               lhs(1,4,aa,i) = - tmp2 * fjac(1,4,i-1)
     >              - tmp1 * njac(1,4,i-1)
               lhs(1,5,aa,i) = - tmp2 * fjac(1,5,i-1)
     >              - tmp1 * njac(1,5,i-1)

               lhs(2,1,aa,i) = - tmp2 * fjac(2,1,i-1)
     >              - tmp1 * njac(2,1,i-1)
               lhs(2,2,aa,i) = - tmp2 * fjac(2,2,i-1)
     >              - tmp1 * njac(2,2,i-1)
     >              - tmp1 * dx2
               lhs(2,3,aa,i) = - tmp2 * fjac(2,3,i-1)
     >              - tmp1 * njac(2,3,i-1)
               lhs(2,4,aa,i) = - tmp2 * fjac(2,4,i-1)
     >              - tmp1 * njac(2,4,i-1)
               lhs(2,5,aa,i) = - tmp2 * fjac(2,5,i-1)
     >              - tmp1 * njac(2,5,i-1)

               lhs(3,1,aa,i) = - tmp2 * fjac(3,1,i-1)
     >              - tmp1 * njac(3,1,i-1)
               lhs(3,2,aa,i) = - tmp2 * fjac(3,2,i-1)
     >              - tmp1 * njac(3,2,i-1)
               lhs(3,3,aa,i) = - tmp2 * fjac(3,3,i-1)
     >              - tmp1 * njac(3,3,i-1)
     >              - tmp1 * dx3 
               lhs(3,4,aa,i) = - tmp2 * fjac(3,4,i-1)
     >              - tmp1 * njac(3,4,i-1)
               lhs(3,5,aa,i) = - tmp2 * fjac(3,5,i-1)
     >              - tmp1 * njac(3,5,i-1)

               lhs(4,1,aa,i) = - tmp2 * fjac(4,1,i-1)
     >              - tmp1 * njac(4,1,i-1)
               lhs(4,2,aa,i) = - tmp2 * fjac(4,2,i-1)
     >              - tmp1 * njac(4,2,i-1)
               lhs(4,3,aa,i) = - tmp2 * fjac(4,3,i-1)
     >              - tmp1 * njac(4,3,i-1)
               lhs(4,4,aa,i) = - tmp2 * fjac(4,4,i-1)
     >              - tmp1 * njac(4,4,i-1)
     >              - tmp1 * dx4
               lhs(4,5,aa,i) = - tmp2 * fjac(4,5,i-1)
     >              - tmp1 * njac(4,5,i-1)

               lhs(5,1,aa,i) = - tmp2 * fjac(5,1,i-1)
     >              - tmp1 * njac(5,1,i-1)
               lhs(5,2,aa,i) = - tmp2 * fjac(5,2,i-1)
     >              - tmp1 * njac(5,2,i-1)
               lhs(5,3,aa,i) = - tmp2 * fjac(5,3,i-1)
     >              - tmp1 * njac(5,3,i-1)
               lhs(5,4,aa,i) = - tmp2 * fjac(5,4,i-1)
     >              - tmp1 * njac(5,4,i-1)
               lhs(5,5,aa,i) = - tmp2 * fjac(5,5,i-1)
     >              - tmp1 * njac(5,5,i-1)
     >              - tmp1 * dx5

               lhs(1,1,bb,i) = 1.0d+00
     >              + tmp1 * 2.0d+00 * njac(1,1,i)
     >              + tmp1 * 2.0d+00 * dx1
               lhs(1,2,bb,i) = tmp1 * 2.0d+00 * njac(1,2,i)
               lhs(1,3,bb,i) = tmp1 * 2.0d+00 * njac(1,3,i)
               lhs(1,4,bb,i) = tmp1 * 2.0d+00 * njac(1,4,i)
               lhs(1,5,bb,i) = tmp1 * 2.0d+00 * njac(1,5,i)

               lhs(2,1,bb,i) = tmp1 * 2.0d+00 * njac(2,1,i)
               lhs(2,2,bb,i) = 1.0d+00
     >              + tmp1 * 2.0d+00 * njac(2,2,i)
     >              + tmp1 * 2.0d+00 * dx2
               lhs(2,3,bb,i) = tmp1 * 2.0d+00 * njac(2,3,i)
               lhs(2,4,bb,i) = tmp1 * 2.0d+00 * njac(2,4,i)
               lhs(2,5,bb,i) = tmp1 * 2.0d+00 * njac(2,5,i)

               lhs(3,1,bb,i) = tmp1 * 2.0d+00 * njac(3,1,i)
               lhs(3,2,bb,i) = tmp1 * 2.0d+00 * njac(3,2,i)
               lhs(3,3,bb,i) = 1.0d+00
     >              + tmp1 * 2.0d+00 * njac(3,3,i)
     >              + tmp1 * 2.0d+00 * dx3
               lhs(3,4,bb,i) = tmp1 * 2.0d+00 * njac(3,4,i)
               lhs(3,5,bb,i) = tmp1 * 2.0d+00 * njac(3,5,i)

               lhs(4,1,bb,i) = tmp1 * 2.0d+00 * njac(4,1,i)
               lhs(4,2,bb,i) = tmp1 * 2.0d+00 * njac(4,2,i)
               lhs(4,3,bb,i) = tmp1 * 2.0d+00 * njac(4,3,i)
               lhs(4,4,bb,i) = 1.0d+00
     >              + tmp1 * 2.0d+00 * njac(4,4,i)
     >              + tmp1 * 2.0d+00 * dx4
               lhs(4,5,bb,i) = tmp1 * 2.0d+00 * njac(4,5,i)

               lhs(5,1,bb,i) = tmp1 * 2.0d+00 * njac(5,1,i)
               lhs(5,2,bb,i) = tmp1 * 2.0d+00 * njac(5,2,i)
               lhs(5,3,bb,i) = tmp1 * 2.0d+00 * njac(5,3,i)
               lhs(5,4,bb,i) = tmp1 * 2.0d+00 * njac(5,4,i)
               lhs(5,5,bb,i) = 1.0d+00
     >              + tmp1 * 2.0d+00 * njac(5,5,i)
     >              + tmp1 * 2.0d+00 * dx5

               lhs(1,1,cc,i) =  tmp2 * fjac(1,1,i+1)
     >              - tmp1 * njac(1,1,i+1)
     >              - tmp1 * dx1
               lhs(1,2,cc,i) =  tmp2 * fjac(1,2,i+1)
     >              - tmp1 * njac(1,2,i+1)
               lhs(1,3,cc,i) =  tmp2 * fjac(1,3,i+1)
     >              - tmp1 * njac(1,3,i+1)
               lhs(1,4,cc,i) =  tmp2 * fjac(1,4,i+1)
     >              - tmp1 * njac(1,4,i+1)
               lhs(1,5,cc,i) =  tmp2 * fjac(1,5,i+1)
     >              - tmp1 * njac(1,5,i+1)

               lhs(2,1,cc,i) =  tmp2 * fjac(2,1,i+1)
     >              - tmp1 * njac(2,1,i+1)
               lhs(2,2,cc,i) =  tmp2 * fjac(2,2,i+1)
     >              - tmp1 * njac(2,2,i+1)
     >              - tmp1 * dx2
               lhs(2,3,cc,i) =  tmp2 * fjac(2,3,i+1)
     >              - tmp1 * njac(2,3,i+1)
               lhs(2,4,cc,i) =  tmp2 * fjac(2,4,i+1)
     >              - tmp1 * njac(2,4,i+1)
               lhs(2,5,cc,i) =  tmp2 * fjac(2,5,i+1)
     >              - tmp1 * njac(2,5,i+1)

               lhs(3,1,cc,i) =  tmp2 * fjac(3,1,i+1)
     >              - tmp1 * njac(3,1,i+1)
               lhs(3,2,cc,i) =  tmp2 * fjac(3,2,i+1)
     >              - tmp1 * njac(3,2,i+1)
               lhs(3,3,cc,i) =  tmp2 * fjac(3,3,i+1)
     >              - tmp1 * njac(3,3,i+1)
     >              - tmp1 * dx3
               lhs(3,4,cc,i) =  tmp2 * fjac(3,4,i+1)
     >              - tmp1 * njac(3,4,i+1)
               lhs(3,5,cc,i) =  tmp2 * fjac(3,5,i+1)
     >              - tmp1 * njac(3,5,i+1)

               lhs(4,1,cc,i) =  tmp2 * fjac(4,1,i+1)
     >              - tmp1 * njac(4,1,i+1)
               lhs(4,2,cc,i) =  tmp2 * fjac(4,2,i+1)
     >              - tmp1 * njac(4,2,i+1)
               lhs(4,3,cc,i) =  tmp2 * fjac(4,3,i+1)
     >              - tmp1 * njac(4,3,i+1)
               lhs(4,4,cc,i) =  tmp2 * fjac(4,4,i+1)
     >              - tmp1 * njac(4,4,i+1)
     >              - tmp1 * dx4
               lhs(4,5,cc,i) =  tmp2 * fjac(4,5,i+1)
     >              - tmp1 * njac(4,5,i+1)

               lhs(5,1,cc,i) =  tmp2 * fjac(5,1,i+1)
     >              - tmp1 * njac(5,1,i+1)
               lhs(5,2,cc,i) =  tmp2 * fjac(5,2,i+1)
     >              - tmp1 * njac(5,2,i+1)
               lhs(5,3,cc,i) =  tmp2 * fjac(5,3,i+1)
     >              - tmp1 * njac(5,3,i+1)
               lhs(5,4,cc,i) =  tmp2 * fjac(5,4,i+1)
     >              - tmp1 * njac(5,4,i+1)
               lhs(5,5,cc,i) =  tmp2 * fjac(5,5,i+1)
     >              - tmp1 * njac(5,5,i+1)
     >              - tmp1 * dx5

            enddo

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     performs guaussian elimination on this cell.
c     
c     assumes that unpacking routines for non-first cells 
c     preload C' and rhs' from previous cell.
c     
c     assumed send happens outside this routine, but that
c     c'(IMAX) and rhs'(IMAX) will be sent to next cell
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     outer most do loops - sweeping in i direction
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     multiply c(0,j,k) by b_inverse and copy back to c
c     multiply rhs(0) by b_inverse(0) and copy to rhs
c---------------------------------------------------------------------
            call binvcrhs( lhs(1,1,bb,0),
     >                        lhs(1,1,cc,0),
     >                        rhs(1,0,j,k) )

c---------------------------------------------------------------------
c     begin inner most do loop
c     do all the elements of the cell unless last 
c---------------------------------------------------------------------
            do i=1,isize-1

c---------------------------------------------------------------------
c     rhs(i) = rhs(i) - A*rhs(i-1)
c---------------------------------------------------------------------
               call matvec_sub(lhs(1,1,aa,i),
     >                         rhs(1,i-1,j,k),rhs(1,i,j,k))

c---------------------------------------------------------------------
c     B(i) = B(i) - C(i-1)*A(i)
c---------------------------------------------------------------------
               call matmul_sub(lhs(1,1,aa,i),
     >                         lhs(1,1,cc,i-1),
     >                         lhs(1,1,bb,i))


c---------------------------------------------------------------------
c     multiply c(i,j,k) by b_inverse and copy back to c
c     multiply rhs(1,j,k) by b_inverse(1,j,k) and copy to rhs
c---------------------------------------------------------------------
               call binvcrhs( lhs(1,1,bb,i),
     >                        lhs(1,1,cc,i),
     >                        rhs(1,i,j,k) )

            enddo

c---------------------------------------------------------------------
c     rhs(isize) = rhs(isize) - A*rhs(isize-1)
c---------------------------------------------------------------------
            call matvec_sub(lhs(1,1,aa,isize),
     >                         rhs(1,isize-1,j,k),rhs(1,isize,j,k))

c---------------------------------------------------------------------
c     B(isize) = B(isize) - C(isize-1)*A(isize)
c---------------------------------------------------------------------
            call matmul_sub(lhs(1,1,aa,isize),
     >                         lhs(1,1,cc,isize-1),
     >                         lhs(1,1,bb,isize))

c---------------------------------------------------------------------
c     multiply rhs() by b_inverse() and copy to rhs
c---------------------------------------------------------------------
            call binvrhs( lhs(1,1,bb,isize),
     >                       rhs(1,isize,j,k) )


c---------------------------------------------------------------------
c     back solve: if last cell, then generate U(isize)=rhs(isize)
c     else assume U(isize) is loaded in un pack backsub_info
c     so just use it
c     after call u(istart) will be sent to next cell
c---------------------------------------------------------------------

            do i=isize-1,0,-1
               do m=1,BLOCK_SIZE
                  do n=1,BLOCK_SIZE
                     rhs(m,i,j,k) = rhs(m,i,j,k) 
     >                    - lhs(m,n,cc,i)*rhs(n,i+1,j,k)
                  enddo
               enddo
            enddo

         enddo
      enddo
      if (timeron) call timer_stop(t_xsolve)

      return
      end
      
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine y_solve

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     Performs line solves in Y direction by first factoring
c     the block-tridiagonal matrix into an upper triangular matrix, 
c     and then performing back substitution to solve for the unknow
c     vectors of each line.  
c     
c     Make sure we treat elements zero to cell_size in the direction
c     of the sweep.
c---------------------------------------------------------------------

      include 'header.h'
      include 'work_lhs.h'

      integer i, j, k, m, n, jsize

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      if (timeron) call timer_start(t_ysolve)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     This function computes the left hand side for the three y-factors   
c---------------------------------------------------------------------

      jsize = grid_points(2)-1

c---------------------------------------------------------------------
c     Compute the indices for storing the tri-diagonal matrix;
c     determine a (labeled f) and n jacobians for cell c
c---------------------------------------------------------------------
      do k = 1, grid_points(3)-2
         do i = 1, grid_points(1)-2
            do j = 0, jsize

               tmp1 = rho_i(i,j,k)
               tmp2 = tmp1 * tmp1
               tmp3 = tmp1 * tmp2

               fjac(1,1,j) = 0.0d+00
               fjac(1,2,j) = 0.0d+00
               fjac(1,3,j) = 1.0d+00
               fjac(1,4,j) = 0.0d+00
               fjac(1,5,j) = 0.0d+00

               fjac(2,1,j) = - ( u(2,i,j,k)*u(3,i,j,k) )
     >              * tmp2
               fjac(2,2,j) = u(3,i,j,k) * tmp1
               fjac(2,3,j) = u(2,i,j,k) * tmp1
               fjac(2,4,j) = 0.0d+00
               fjac(2,5,j) = 0.0d+00

               fjac(3,1,j) = - ( u(3,i,j,k)*u(3,i,j,k)*tmp2)
     >              + c2 * qs(i,j,k)
               fjac(3,2,j) = - c2 *  u(2,i,j,k) * tmp1
               fjac(3,3,j) = ( 2.0d+00 - c2 )
     >              *  u(3,i,j,k) * tmp1 
               fjac(3,4,j) = - c2 * u(4,i,j,k) * tmp1 
               fjac(3,5,j) = c2

               fjac(4,1,j) = - ( u(3,i,j,k)*u(4,i,j,k) )
     >              * tmp2
               fjac(4,2,j) = 0.0d+00
               fjac(4,3,j) = u(4,i,j,k) * tmp1
               fjac(4,4,j) = u(3,i,j,k) * tmp1
               fjac(4,5,j) = 0.0d+00

               fjac(5,1,j) = ( c2 * 2.0d0 * square(i,j,k)
     >              - c1 * u(5,i,j,k) )
     >              * u(3,i,j,k) * tmp2
               fjac(5,2,j) = - c2 * u(2,i,j,k)*u(3,i,j,k) 
     >              * tmp2
               fjac(5,3,j) = c1 * u(5,i,j,k) * tmp1 
     >              - c2 
     >              * ( qs(i,j,k)
     >              + u(3,i,j,k)*u(3,i,j,k) * tmp2 )
               fjac(5,4,j) = - c2 * ( u(3,i,j,k)*u(4,i,j,k) )
     >              * tmp2
               fjac(5,5,j) = c1 * u(3,i,j,k) * tmp1 

               njac(1,1,j) = 0.0d+00
               njac(1,2,j) = 0.0d+00
               njac(1,3,j) = 0.0d+00
               njac(1,4,j) = 0.0d+00
               njac(1,5,j) = 0.0d+00

               njac(2,1,j) = - c3c4 * tmp2 * u(2,i,j,k)
               njac(2,2,j) =   c3c4 * tmp1
               njac(2,3,j) =   0.0d+00
               njac(2,4,j) =   0.0d+00
               njac(2,5,j) =   0.0d+00

               njac(3,1,j) = - con43 * c3c4 * tmp2 * u(3,i,j,k)
               njac(3,2,j) =   0.0d+00
               njac(3,3,j) =   con43 * c3c4 * tmp1
               njac(3,4,j) =   0.0d+00
               njac(3,5,j) =   0.0d+00

               njac(4,1,j) = - c3c4 * tmp2 * u(4,i,j,k)
               njac(4,2,j) =   0.0d+00
               njac(4,3,j) =   0.0d+00
               njac(4,4,j) =   c3c4 * tmp1
               njac(4,5,j) =   0.0d+00

               njac(5,1,j) = - (  c3c4
     >              - c1345 ) * tmp3 * (u(2,i,j,k)**2)
     >              - ( con43 * c3c4
     >              - c1345 ) * tmp3 * (u(3,i,j,k)**2)
     >              - ( c3c4 - c1345 ) * tmp3 * (u(4,i,j,k)**2)
     >              - c1345 * tmp2 * u(5,i,j,k)

               njac(5,2,j) = (  c3c4 - c1345 ) * tmp2 * u(2,i,j,k)
               njac(5,3,j) = ( con43 * c3c4
     >              - c1345 ) * tmp2 * u(3,i,j,k)
               njac(5,4,j) = ( c3c4 - c1345 ) * tmp2 * u(4,i,j,k)
               njac(5,5,j) = ( c1345 ) * tmp1

            enddo

c---------------------------------------------------------------------
c     now joacobians set, so form left hand side in y direction
c---------------------------------------------------------------------
            call lhsinit(lhs, jsize)
            do j = 1, jsize-1

               tmp1 = dt * ty1
               tmp2 = dt * ty2

               lhs(1,1,aa,j) = - tmp2 * fjac(1,1,j-1)
     >              - tmp1 * njac(1,1,j-1)
     >              - tmp1 * dy1 
               lhs(1,2,aa,j) = - tmp2 * fjac(1,2,j-1)
     >              - tmp1 * njac(1,2,j-1)
               lhs(1,3,aa,j) = - tmp2 * fjac(1,3,j-1)
     >              - tmp1 * njac(1,3,j-1)
               lhs(1,4,aa,j) = - tmp2 * fjac(1,4,j-1)
     >              - tmp1 * njac(1,4,j-1)
               lhs(1,5,aa,j) = - tmp2 * fjac(1,5,j-1)
     >              - tmp1 * njac(1,5,j-1)

               lhs(2,1,aa,j) = - tmp2 * fjac(2,1,j-1)
     >              - tmp1 * njac(2,1,j-1)
               lhs(2,2,aa,j) = - tmp2 * fjac(2,2,j-1)
     >              - tmp1 * njac(2,2,j-1)
     >              - tmp1 * dy2
               lhs(2,3,aa,j) = - tmp2 * fjac(2,3,j-1)
     >              - tmp1 * njac(2,3,j-1)
               lhs(2,4,aa,j) = - tmp2 * fjac(2,4,j-1)
     >              - tmp1 * njac(2,4,j-1)
               lhs(2,5,aa,j) = - tmp2 * fjac(2,5,j-1)
     >              - tmp1 * njac(2,5,j-1)

               lhs(3,1,aa,j) = - tmp2 * fjac(3,1,j-1)
     >              - tmp1 * njac(3,1,j-1)
               lhs(3,2,aa,j) = - tmp2 * fjac(3,2,j-1)
     >              - tmp1 * njac(3,2,j-1)
               lhs(3,3,aa,j) = - tmp2 * fjac(3,3,j-1)
     >              - tmp1 * njac(3,3,j-1)
     >              - tmp1 * dy3 
               lhs(3,4,aa,j) = - tmp2 * fjac(3,4,j-1)
     >              - tmp1 * njac(3,4,j-1)
               lhs(3,5,aa,j) = - tmp2 * fjac(3,5,j-1)
     >              - tmp1 * njac(3,5,j-1)

               lhs(4,1,aa,j) = - tmp2 * fjac(4,1,j-1)
     >              - tmp1 * njac(4,1,j-1)
               lhs(4,2,aa,j) = - tmp2 * fjac(4,2,j-1)
     >              - tmp1 * njac(4,2,j-1)
               lhs(4,3,aa,j) = - tmp2 * fjac(4,3,j-1)
     >              - tmp1 * njac(4,3,j-1)
               lhs(4,4,aa,j) = - tmp2 * fjac(4,4,j-1)
     >              - tmp1 * njac(4,4,j-1)
     >              - tmp1 * dy4
               lhs(4,5,aa,j) = - tmp2 * fjac(4,5,j-1)
     >              - tmp1 * njac(4,5,j-1)

               lhs(5,1,aa,j) = - tmp2 * fjac(5,1,j-1)
     >              - tmp1 * njac(5,1,j-1)
               lhs(5,2,aa,j) = - tmp2 * fjac(5,2,j-1)
     >              - tmp1 * njac(5,2,j-1)
               lhs(5,3,aa,j) = - tmp2 * fjac(5,3,j-1)
     >              - tmp1 * njac(5,3,j-1)
               lhs(5,4,aa,j) = - tmp2 * fjac(5,4,j-1)
     >              - tmp1 * njac(5,4,j-1)
               lhs(5,5,aa,j) = - tmp2 * fjac(5,5,j-1)
     >              - tmp1 * njac(5,5,j-1)
     >              - tmp1 * dy5

               lhs(1,1,bb,j) = 1.0d+00
     >              + tmp1 * 2.0d+00 * njac(1,1,j)
     >              + tmp1 * 2.0d+00 * dy1
               lhs(1,2,bb,j) = tmp1 * 2.0d+00 * njac(1,2,j)
               lhs(1,3,bb,j) = tmp1 * 2.0d+00 * njac(1,3,j)
               lhs(1,4,bb,j) = tmp1 * 2.0d+00 * njac(1,4,j)
               lhs(1,5,bb,j) = tmp1 * 2.0d+00 * njac(1,5,j)

               lhs(2,1,bb,j) = tmp1 * 2.0d+00 * njac(2,1,j)
               lhs(2,2,bb,j) = 1.0d+00
     >              + tmp1 * 2.0d+00 * njac(2,2,j)
     >              + tmp1 * 2.0d+00 * dy2
               lhs(2,3,bb,j) = tmp1 * 2.0d+00 * njac(2,3,j)
               lhs(2,4,bb,j) = tmp1 * 2.0d+00 * njac(2,4,j)
               lhs(2,5,bb,j) = tmp1 * 2.0d+00 * njac(2,5,j)

               lhs(3,1,bb,j) = tmp1 * 2.0d+00 * njac(3,1,j)
               lhs(3,2,bb,j) = tmp1 * 2.0d+00 * njac(3,2,j)
               lhs(3,3,bb,j) = 1.0d+00
     >              + tmp1 * 2.0d+00 * njac(3,3,j)
     >              + tmp1 * 2.0d+00 * dy3
               lhs(3,4,bb,j) = tmp1 * 2.0d+00 * njac(3,4,j)
               lhs(3,5,bb,j) = tmp1 * 2.0d+00 * njac(3,5,j)

               lhs(4,1,bb,j) = tmp1 * 2.0d+00 * njac(4,1,j)
               lhs(4,2,bb,j) = tmp1 * 2.0d+00 * njac(4,2,j)
               lhs(4,3,bb,j) = tmp1 * 2.0d+00 * njac(4,3,j)
               lhs(4,4,bb,j) = 1.0d+00
     >              + tmp1 * 2.0d+00 * njac(4,4,j)
     >              + tmp1 * 2.0d+00 * dy4
               lhs(4,5,bb,j) = tmp1 * 2.0d+00 * njac(4,5,j)

               lhs(5,1,bb,j) = tmp1 * 2.0d+00 * njac(5,1,j)
               lhs(5,2,bb,j) = tmp1 * 2.0d+00 * njac(5,2,j)
               lhs(5,3,bb,j) = tmp1 * 2.0d+00 * njac(5,3,j)
               lhs(5,4,bb,j) = tmp1 * 2.0d+00 * njac(5,4,j)
               lhs(5,5,bb,j) = 1.0d+00
     >              + tmp1 * 2.0d+00 * njac(5,5,j) 
     >              + tmp1 * 2.0d+00 * dy5

               lhs(1,1,cc,j) =  tmp2 * fjac(1,1,j+1)
     >              - tmp1 * njac(1,1,j+1)
     >              - tmp1 * dy1
               lhs(1,2,cc,j) =  tmp2 * fjac(1,2,j+1)
     >              - tmp1 * njac(1,2,j+1)
               lhs(1,3,cc,j) =  tmp2 * fjac(1,3,j+1)
     >              - tmp1 * njac(1,3,j+1)
               lhs(1,4,cc,j) =  tmp2 * fjac(1,4,j+1)
     >              - tmp1 * njac(1,4,j+1)
               lhs(1,5,cc,j) =  tmp2 * fjac(1,5,j+1)
     >              - tmp1 * njac(1,5,j+1)

               lhs(2,1,cc,j) =  tmp2 * fjac(2,1,j+1)
     >              - tmp1 * njac(2,1,j+1)
               lhs(2,2,cc,j) =  tmp2 * fjac(2,2,j+1)
     >              - tmp1 * njac(2,2,j+1)
     >              - tmp1 * dy2
               lhs(2,3,cc,j) =  tmp2 * fjac(2,3,j+1)
     >              - tmp1 * njac(2,3,j+1)
               lhs(2,4,cc,j) =  tmp2 * fjac(2,4,j+1)
     >              - tmp1 * njac(2,4,j+1)
               lhs(2,5,cc,j) =  tmp2 * fjac(2,5,j+1)
     >              - tmp1 * njac(2,5,j+1)

               lhs(3,1,cc,j) =  tmp2 * fjac(3,1,j+1)
     >              - tmp1 * njac(3,1,j+1)
               lhs(3,2,cc,j) =  tmp2 * fjac(3,2,j+1)
     >              - tmp1 * njac(3,2,j+1)
               lhs(3,3,cc,j) =  tmp2 * fjac(3,3,j+1)
     >              - tmp1 * njac(3,3,j+1)
     >              - tmp1 * dy3
               lhs(3,4,cc,j) =  tmp2 * fjac(3,4,j+1)
     >              - tmp1 * njac(3,4,j+1)
               lhs(3,5,cc,j) =  tmp2 * fjac(3,5,j+1)
     >              - tmp1 * njac(3,5,j+1)

               lhs(4,1,cc,j) =  tmp2 * fjac(4,1,j+1)
     >              - tmp1 * njac(4,1,j+1)
               lhs(4,2,cc,j) =  tmp2 * fjac(4,2,j+1)
     >              - tmp1 * njac(4,2,j+1)
               lhs(4,3,cc,j) =  tmp2 * fjac(4,3,j+1)
     >              - tmp1 * njac(4,3,j+1)
               lhs(4,4,cc,j) =  tmp2 * fjac(4,4,j+1)
     >              - tmp1 * njac(4,4,j+1)
     >              - tmp1 * dy4
               lhs(4,5,cc,j) =  tmp2 * fjac(4,5,j+1)
     >              - tmp1 * njac(4,5,j+1)

               lhs(5,1,cc,j) =  tmp2 * fjac(5,1,j+1)
     >              - tmp1 * njac(5,1,j+1)
               lhs(5,2,cc,j) =  tmp2 * fjac(5,2,j+1)
     >              - tmp1 * njac(5,2,j+1)
               lhs(5,3,cc,j) =  tmp2 * fjac(5,3,j+1)
     >              - tmp1 * njac(5,3,j+1)
               lhs(5,4,cc,j) =  tmp2 * fjac(5,4,j+1)
     >              - tmp1 * njac(5,4,j+1)
               lhs(5,5,cc,j) =  tmp2 * fjac(5,5,j+1)
     >              - tmp1 * njac(5,5,j+1)
     >              - tmp1 * dy5

            enddo

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     performs guaussian elimination on this cell.
c     
c     assumes that unpacking routines for non-first cells 
c     preload C' and rhs' from previous cell.
c     
c     assumed send happens outside this routine, but that
c     c'(JMAX) and rhs'(JMAX) will be sent to next cell
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     multiply c(i,0,k) by b_inverse and copy back to c
c     multiply rhs(0) by b_inverse(0) and copy to rhs
c---------------------------------------------------------------------
            call binvcrhs( lhs(1,1,bb,0),
     >                        lhs(1,1,cc,0),
     >                        rhs(1,i,0,k) )

c---------------------------------------------------------------------
c     begin inner most do loop
c     do all the elements of the cell unless last 
c---------------------------------------------------------------------
            do j=1,jsize-1

c---------------------------------------------------------------------
c     subtract A*lhs_vector(j-1) from lhs_vector(j)
c     
c     rhs(j) = rhs(j) - A*rhs(j-1)
c---------------------------------------------------------------------
               call matvec_sub(lhs(1,1,aa,j),
     >                         rhs(1,i,j-1,k),rhs(1,i,j,k))

c---------------------------------------------------------------------
c     B(j) = B(j) - C(j-1)*A(j)
c---------------------------------------------------------------------
               call matmul_sub(lhs(1,1,aa,j),
     >                         lhs(1,1,cc,j-1),
     >                         lhs(1,1,bb,j))

c---------------------------------------------------------------------
c     multiply c(i,j,k) by b_inverse and copy back to c
c     multiply rhs(i,1,k) by b_inverse(i,1,k) and copy to rhs
c---------------------------------------------------------------------
               call binvcrhs( lhs(1,1,bb,j),
     >                        lhs(1,1,cc,j),
     >                        rhs(1,i,j,k) )

            enddo


c---------------------------------------------------------------------
c     rhs(jsize) = rhs(jsize) - A*rhs(jsize-1)
c---------------------------------------------------------------------
            call matvec_sub(lhs(1,1,aa,jsize),
     >                         rhs(1,i,jsize-1,k),rhs(1,i,jsize,k))

c---------------------------------------------------------------------
c     B(jsize) = B(jsize) - C(jsize-1)*A(jsize)
c     call matmul_sub(aa,i,jsize,k,c,
c     $              cc,i,jsize-1,k,c,bb,i,jsize,k)
c---------------------------------------------------------------------
            call matmul_sub(lhs(1,1,aa,jsize),
     >                         lhs(1,1,cc,jsize-1),
     >                         lhs(1,1,bb,jsize))

c---------------------------------------------------------------------
c     multiply rhs(jsize) by b_inverse(jsize) and copy to rhs
c---------------------------------------------------------------------
            call binvrhs( lhs(1,1,bb,jsize),
     >                       rhs(1,i,jsize,k) )


c---------------------------------------------------------------------
c     back solve: if last cell, then generate U(jsize)=rhs(jsize)
c     else assume U(jsize) is loaded in un pack backsub_info
c     so just use it
c     after call u(jstart) will be sent to next cell
c---------------------------------------------------------------------
      
            do j=jsize-1,0,-1
               do m=1,BLOCK_SIZE
                  do n=1,BLOCK_SIZE
                     rhs(m,i,j,k) = rhs(m,i,j,k) 
     >                    - lhs(m,n,cc,j)*rhs(n,i,j+1,k)
                  enddo
               enddo
            enddo

         enddo
      enddo
      if (timeron) call timer_stop(t_ysolve)

      return
      end


c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine z_solve

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     Performs line solves in Z direction by first factoring
c     the block-tridiagonal matrix into an upper triangular matrix, 
c     and then performing back substitution to solve for the unknow
c     vectors of each line.  
c     
c     Make sure we treat elements zero to cell_size in the direction
c     of the sweep.
c---------------------------------------------------------------------

      include 'header.h'
      include 'work_lhs.h'

      integer i, j, k, m, n, ksize
      
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      if (timeron) call timer_start(t_zsolve)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     This function computes the left hand side for the three z-factors   
c---------------------------------------------------------------------

      ksize = grid_points(3)-1

c---------------------------------------------------------------------
c     Compute the indices for storing the block-diagonal matrix;
c     determine c (labeled f) and s jacobians
c---------------------------------------------------------------------
      do j = 1, grid_points(2)-2
         do i = 1, grid_points(1)-2
            do k = 0, ksize

               tmp1 = 1.0d+00 / u(1,i,j,k)
               tmp2 = tmp1 * tmp1
               tmp3 = tmp1 * tmp2

               fjac(1,1,k) = 0.0d+00
               fjac(1,2,k) = 0.0d+00
               fjac(1,3,k) = 0.0d+00
               fjac(1,4,k) = 1.0d+00
               fjac(1,5,k) = 0.0d+00

               fjac(2,1,k) = - ( u(2,i,j,k)*u(4,i,j,k) ) 
     >              * tmp2 
               fjac(2,2,k) = u(4,i,j,k) * tmp1
               fjac(2,3,k) = 0.0d+00
               fjac(2,4,k) = u(2,i,j,k) * tmp1
               fjac(2,5,k) = 0.0d+00

               fjac(3,1,k) = - ( u(3,i,j,k)*u(4,i,j,k) )
     >              * tmp2 
               fjac(3,2,k) = 0.0d+00
               fjac(3,3,k) = u(4,i,j,k) * tmp1
               fjac(3,4,k) = u(3,i,j,k) * tmp1
               fjac(3,5,k) = 0.0d+00

               fjac(4,1,k) = - (u(4,i,j,k)*u(4,i,j,k) * tmp2 ) 
     >              + c2 * qs(i,j,k)
               fjac(4,2,k) = - c2 *  u(2,i,j,k) * tmp1 
               fjac(4,3,k) = - c2 *  u(3,i,j,k) * tmp1
               fjac(4,4,k) = ( 2.0d+00 - c2 )
     >              *  u(4,i,j,k) * tmp1 
               fjac(4,5,k) = c2

               fjac(5,1,k) = ( c2 * 2.0d0 * square(i,j,k) 
     >              - c1 * u(5,i,j,k) )
     >              * u(4,i,j,k) * tmp2
               fjac(5,2,k) = - c2 * ( u(2,i,j,k)*u(4,i,j,k) )
     >              * tmp2 
               fjac(5,3,k) = - c2 * ( u(3,i,j,k)*u(4,i,j,k) )
     >              * tmp2
               fjac(5,4,k) = c1 * ( u(5,i,j,k) * tmp1 )
     >              - c2
     >              * ( qs(i,j,k)
     >              + u(4,i,j,k)*u(4,i,j,k) * tmp2 )
               fjac(5,5,k) = c1 * u(4,i,j,k) * tmp1

               njac(1,1,k) = 0.0d+00
               njac(1,2,k) = 0.0d+00
               njac(1,3,k) = 0.0d+00
               njac(1,4,k) = 0.0d+00
               njac(1,5,k) = 0.0d+00

               njac(2,1,k) = - c3c4 * tmp2 * u(2,i,j,k)
               njac(2,2,k) =   c3c4 * tmp1
               njac(2,3,k) =   0.0d+00
               njac(2,4,k) =   0.0d+00
               njac(2,5,k) =   0.0d+00

               njac(3,1,k) = - c3c4 * tmp2 * u(3,i,j,k)
               njac(3,2,k) =   0.0d+00
               njac(3,3,k) =   c3c4 * tmp1
               njac(3,4,k) =   0.0d+00
               njac(3,5,k) =   0.0d+00

               njac(4,1,k) = - con43 * c3c4 * tmp2 * u(4,i,j,k)
               njac(4,2,k) =   0.0d+00
               njac(4,3,k) =   0.0d+00
               njac(4,4,k) =   con43 * c3 * c4 * tmp1
               njac(4,5,k) =   0.0d+00

               njac(5,1,k) = - (  c3c4
     >              - c1345 ) * tmp3 * (u(2,i,j,k)**2)
     >              - ( c3c4 - c1345 ) * tmp3 * (u(3,i,j,k)**2)
     >              - ( con43 * c3c4
     >              - c1345 ) * tmp3 * (u(4,i,j,k)**2)
     >              - c1345 * tmp2 * u(5,i,j,k)

               njac(5,2,k) = (  c3c4 - c1345 ) * tmp2 * u(2,i,j,k)
               njac(5,3,k) = (  c3c4 - c1345 ) * tmp2 * u(3,i,j,k)
               njac(5,4,k) = ( con43 * c3c4
     >              - c1345 ) * tmp2 * u(4,i,j,k)
               njac(5,5,k) = ( c1345 )* tmp1


            enddo

c---------------------------------------------------------------------
c     now jacobians set, so form left hand side in z direction
c---------------------------------------------------------------------
            call lhsinit(lhs, ksize)
            do k = 1, ksize-1

               tmp1 = dt * tz1
               tmp2 = dt * tz2

               lhs(1,1,aa,k) = - tmp2 * fjac(1,1,k-1)
     >              - tmp1 * njac(1,1,k-1)
     >              - tmp1 * dz1 
               lhs(1,2,aa,k) = - tmp2 * fjac(1,2,k-1)
     >              - tmp1 * njac(1,2,k-1)
               lhs(1,3,aa,k) = - tmp2 * fjac(1,3,k-1)
     >              - tmp1 * njac(1,3,k-1)
               lhs(1,4,aa,k) = - tmp2 * fjac(1,4,k-1)
     >              - tmp1 * njac(1,4,k-1)
               lhs(1,5,aa,k) = - tmp2 * fjac(1,5,k-1)
     >              - tmp1 * njac(1,5,k-1)

               lhs(2,1,aa,k) = - tmp2 * fjac(2,1,k-1)
     >              - tmp1 * njac(2,1,k-1)
               lhs(2,2,aa,k) = - tmp2 * fjac(2,2,k-1)
     >              - tmp1 * njac(2,2,k-1)
     >              - tmp1 * dz2
               lhs(2,3,aa,k) = - tmp2 * fjac(2,3,k-1)
     >              - tmp1 * njac(2,3,k-1)
               lhs(2,4,aa,k) = - tmp2 * fjac(2,4,k-1)
     >              - tmp1 * njac(2,4,k-1)
               lhs(2,5,aa,k) = - tmp2 * fjac(2,5,k-1)
     >              - tmp1 * njac(2,5,k-1)

               lhs(3,1,aa,k) = - tmp2 * fjac(3,1,k-1)
     >              - tmp1 * njac(3,1,k-1)
               lhs(3,2,aa,k) = - tmp2 * fjac(3,2,k-1)
     >              - tmp1 * njac(3,2,k-1)
               lhs(3,3,aa,k) = - tmp2 * fjac(3,3,k-1)
     >              - tmp1 * njac(3,3,k-1)
     >              - tmp1 * dz3 
               lhs(3,4,aa,k) = - tmp2 * fjac(3,4,k-1)
     >              - tmp1 * njac(3,4,k-1)
               lhs(3,5,aa,k) = - tmp2 * fjac(3,5,k-1)
     >              - tmp1 * njac(3,5,k-1)

               lhs(4,1,aa,k) = - tmp2 * fjac(4,1,k-1)
     >              - tmp1 * njac(4,1,k-1)
               lhs(4,2,aa,k) = - tmp2 * fjac(4,2,k-1)
     >              - tmp1 * njac(4,2,k-1)
               lhs(4,3,aa,k) = - tmp2 * fjac(4,3,k-1)
     >              - tmp1 * njac(4,3,k-1)
               lhs(4,4,aa,k) = - tmp2 * fjac(4,4,k-1)
     >              - tmp1 * njac(4,4,k-1)
     >              - tmp1 * dz4
               lhs(4,5,aa,k) = - tmp2 * fjac(4,5,k-1)
     >              - tmp1 * njac(4,5,k-1)

               lhs(5,1,aa,k) = - tmp2 * fjac(5,1,k-1)
     >              - tmp1 * njac(5,1,k-1)
               lhs(5,2,aa,k) = - tmp2 * fjac(5,2,k-1)
     >              - tmp1 * njac(5,2,k-1)
               lhs(5,3,aa,k) = - tmp2 * fjac(5,3,k-1)
     >              - tmp1 * njac(5,3,k-1)
               lhs(5,4,aa,k) = - tmp2 * fjac(5,4,k-1)
     >              - tmp1 * njac(5,4,k-1)
               lhs(5,5,aa,k) = - tmp2 * fjac(5,5,k-1)
     >              - tmp1 * njac(5,5,k-1)
     >              - tmp1 * dz5

               lhs(1,1,bb,k) = 1.0d+00
     >              + tmp1 * 2.0d+00 * njac(1,1,k)
     >              + tmp1 * 2.0d+00 * dz1
               lhs(1,2,bb,k) = tmp1 * 2.0d+00 * njac(1,2,k)
               lhs(1,3,bb,k) = tmp1 * 2.0d+00 * njac(1,3,k)
               lhs(1,4,bb,k) = tmp1 * 2.0d+00 * njac(1,4,k)
               lhs(1,5,bb,k) = tmp1 * 2.0d+00 * njac(1,5,k)

               lhs(2,1,bb,k) = tmp1 * 2.0d+00 * njac(2,1,k)
               lhs(2,2,bb,k) = 1.0d+00
     >              + tmp1 * 2.0d+00 * njac(2,2,k)
     >              + tmp1 * 2.0d+00 * dz2
               lhs(2,3,bb,k) = tmp1 * 2.0d+00 * njac(2,3,k)
               lhs(2,4,bb,k) = tmp1 * 2.0d+00 * njac(2,4,k)
               lhs(2,5,bb,k) = tmp1 * 2.0d+00 * njac(2,5,k)

               lhs(3,1,bb,k) = tmp1 * 2.0d+00 * njac(3,1,k)
               lhs(3,2,bb,k) = tmp1 * 2.0d+00 * njac(3,2,k)
               lhs(3,3,bb,k) = 1.0d+00
     >              + tmp1 * 2.0d+00 * njac(3,3,k)
     >              + tmp1 * 2.0d+00 * dz3
               lhs(3,4,bb,k) = tmp1 * 2.0d+00 * njac(3,4,k)
               lhs(3,5,bb,k) = tmp1 * 2.0d+00 * njac(3,5,k)

               lhs(4,1,bb,k) = tmp1 * 2.0d+00 * njac(4,1,k)
               lhs(4,2,bb,k) = tmp1 * 2.0d+00 * njac(4,2,k)
               lhs(4,3,bb,k) = tmp1 * 2.0d+00 * njac(4,3,k)
               lhs(4,4,bb,k) = 1.0d+00
     >              + tmp1 * 2.0d+00 * njac(4,4,k)
     >              + tmp1 * 2.0d+00 * dz4
               lhs(4,5,bb,k) = tmp1 * 2.0d+00 * njac(4,5,k)

               lhs(5,1,bb,k) = tmp1 * 2.0d+00 * njac(5,1,k)
               lhs(5,2,bb,k) = tmp1 * 2.0d+00 * njac(5,2,k)
               lhs(5,3,bb,k) = tmp1 * 2.0d+00 * njac(5,3,k)
               lhs(5,4,bb,k) = tmp1 * 2.0d+00 * njac(5,4,k)
               lhs(5,5,bb,k) = 1.0d+00
     >              + tmp1 * 2.0d+00 * njac(5,5,k) 
     >              + tmp1 * 2.0d+00 * dz5

               lhs(1,1,cc,k) =  tmp2 * fjac(1,1,k+1)
     >              - tmp1 * njac(1,1,k+1)
     >              - tmp1 * dz1
               lhs(1,2,cc,k) =  tmp2 * fjac(1,2,k+1)
     >              - tmp1 * njac(1,2,k+1)
               lhs(1,3,cc,k) =  tmp2 * fjac(1,3,k+1)
     >              - tmp1 * njac(1,3,k+1)
               lhs(1,4,cc,k) =  tmp2 * fjac(1,4,k+1)
     >              - tmp1 * njac(1,4,k+1)
               lhs(1,5,cc,k) =  tmp2 * fjac(1,5,k+1)
     >              - tmp1 * njac(1,5,k+1)

               lhs(2,1,cc,k) =  tmp2 * fjac(2,1,k+1)
     >              - tmp1 * njac(2,1,k+1)
               lhs(2,2,cc,k) =  tmp2 * fjac(2,2,k+1)
     >              - tmp1 * njac(2,2,k+1)
     >              - tmp1 * dz2
               lhs(2,3,cc,k) =  tmp2 * fjac(2,3,k+1)
     >              - tmp1 * njac(2,3,k+1)
               lhs(2,4,cc,k) =  tmp2 * fjac(2,4,k+1)
     >              - tmp1 * njac(2,4,k+1)
               lhs(2,5,cc,k) =  tmp2 * fjac(2,5,k+1)
     >              - tmp1 * njac(2,5,k+1)

               lhs(3,1,cc,k) =  tmp2 * fjac(3,1,k+1)
     >              - tmp1 * njac(3,1,k+1)
               lhs(3,2,cc,k) =  tmp2 * fjac(3,2,k+1)
     >              - tmp1 * njac(3,2,k+1)
               lhs(3,3,cc,k) =  tmp2 * fjac(3,3,k+1)
     >              - tmp1 * njac(3,3,k+1)
     >              - tmp1 * dz3
               lhs(3,4,cc,k) =  tmp2 * fjac(3,4,k+1)
     >              - tmp1 * njac(3,4,k+1)
               lhs(3,5,cc,k) =  tmp2 * fjac(3,5,k+1)
     >              - tmp1 * njac(3,5,k+1)

               lhs(4,1,cc,k) =  tmp2 * fjac(4,1,k+1)
     >              - tmp1 * njac(4,1,k+1)
               lhs(4,2,cc,k) =  tmp2 * fjac(4,2,k+1)
     >              - tmp1 * njac(4,2,k+1)
               lhs(4,3,cc,k) =  tmp2 * fjac(4,3,k+1)
     >              - tmp1 * njac(4,3,k+1)
               lhs(4,4,cc,k) =  tmp2 * fjac(4,4,k+1)
     >              - tmp1 * njac(4,4,k+1)
     >              - tmp1 * dz4
               lhs(4,5,cc,k) =  tmp2 * fjac(4,5,k+1)
     >              - tmp1 * njac(4,5,k+1)

               lhs(5,1,cc,k) =  tmp2 * fjac(5,1,k+1)
     >              - tmp1 * njac(5,1,k+1)
               lhs(5,2,cc,k) =  tmp2 * fjac(5,2,k+1)
     >              - tmp1 * njac(5,2,k+1)
               lhs(5,3,cc,k) =  tmp2 * fjac(5,3,k+1)
     >              - tmp1 * njac(5,3,k+1)
               lhs(5,4,cc,k) =  tmp2 * fjac(5,4,k+1)
     >              - tmp1 * njac(5,4,k+1)
               lhs(5,5,cc,k) =  tmp2 * fjac(5,5,k+1)
     >              - tmp1 * njac(5,5,k+1)
     >              - tmp1 * dz5

            enddo

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     performs guaussian elimination on this cell.
c     
c     assumes that unpacking routines for non-first cells 
c     preload C' and rhs' from previous cell.
c     
c     assumed send happens outside this routine, but that
c     c'(KMAX) and rhs'(KMAX) will be sent to next cell.
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     outer most do loops - sweeping in i direction
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     multiply c(i,j,0) by b_inverse and copy back to c
c     multiply rhs(0) by b_inverse(0) and copy to rhs
c---------------------------------------------------------------------
            call binvcrhs( lhs(1,1,bb,0),
     >                        lhs(1,1,cc,0),
     >                        rhs(1,i,j,0) )


c---------------------------------------------------------------------
c     begin inner most do loop
c     do all the elements of the cell unless last 
c---------------------------------------------------------------------
            do k=1,ksize-1

c---------------------------------------------------------------------
c     subtract A*lhs_vector(k-1) from lhs_vector(k)
c     
c     rhs(k) = rhs(k) - A*rhs(k-1)
c---------------------------------------------------------------------
               call matvec_sub(lhs(1,1,aa,k),
     >                         rhs(1,i,j,k-1),rhs(1,i,j,k))

c---------------------------------------------------------------------
c     B(k) = B(k) - C(k-1)*A(k)
c     call matmul_sub(aa,i,j,k,c,cc,i,j,k-1,c,bb,i,j,k)
c---------------------------------------------------------------------
               call matmul_sub(lhs(1,1,aa,k),
     >                         lhs(1,1,cc,k-1),
     >                         lhs(1,1,bb,k))

c---------------------------------------------------------------------
c     multiply c(i,j,k) by b_inverse and copy back to c
c     multiply rhs(i,j,1) by b_inverse(i,j,1) and copy to rhs
c---------------------------------------------------------------------
               call binvcrhs( lhs(1,1,bb,k),
     >                        lhs(1,1,cc,k),
     >                        rhs(1,i,j,k) )

            enddo

c---------------------------------------------------------------------
c     Now finish up special cases for last cell
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     rhs(ksize) = rhs(ksize) - A*rhs(ksize-1)
c---------------------------------------------------------------------
            call matvec_sub(lhs(1,1,aa,ksize),
     >                         rhs(1,i,j,ksize-1),rhs(1,i,j,ksize))

c---------------------------------------------------------------------
c     B(ksize) = B(ksize) - C(ksize-1)*A(ksize)
c     call matmul_sub(aa,i,j,ksize,c,
c     $              cc,i,j,ksize-1,c,bb,i,j,ksize)
c---------------------------------------------------------------------
            call matmul_sub(lhs(1,1,aa,ksize),
     >                         lhs(1,1,cc,ksize-1),
     >                         lhs(1,1,bb,ksize))

c---------------------------------------------------------------------
c     multiply rhs(ksize) by b_inverse(ksize) and copy to rhs
c---------------------------------------------------------------------
            call binvrhs( lhs(1,1,bb,ksize),
     >                       rhs(1,i,j,ksize) )


c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     back solve: if last cell, then generate U(ksize)=rhs(ksize)
c     else assume U(ksize) is loaded in un pack backsub_info
c     so just use it
c     after call u(kstart) will be sent to next cell
c---------------------------------------------------------------------

            do k=ksize-1,0,-1
               do m=1,BLOCK_SIZE
                  do n=1,BLOCK_SIZE
                     rhs(m,i,j,k) = rhs(m,i,j,k) 
     >                    - lhs(m,n,cc,k)*rhs(n,i,j,k+1)
                  enddo
               enddo
            enddo

         enddo
      enddo
      if (timeron) call timer_stop(t_zsolve)

      return
      end
