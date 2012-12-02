      program main

c*********************************************************************72
c
cc MAIN is the main program for EXPOKIT_PRB.
c
c  Author:
c
c    Roger Sidje
c
c  Reference:
c
c    Roger Sidje,
c    EXPOKIT: Software Package for Computing Matrix Exponentials,
c    ACM Transactions on Mathematical Software,
c    Volume 24, Number 1, 1998, pages 130-156.
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EXPOKIT:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the EXPOKIT library.'

      call sample_b ( )
      call sample_d ( )
      call sample_g ( )
      call sample_m ( )
      call sample_p ( )
      call sample_z ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EXPOKIT:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine sample_b ( )

c*********************************************************************72
c
cc SAMPLE_B illustrates the use of DSEXPV.
c
c  Discussion:
c
c    Forward-backward problem (Example 6.4 in the Expokit report).
c
      implicit none

      external dgcoov, dgcrsv, dgccsv

      double precision tic, tac, clock
c
c  matrix data.
c  BEWARE: these values must match those in dgmatv.f
c
      integer n, nz, nmax, nzmax
      parameter( nmax = 5000, nzmax = 600000 )
      integer ia(nzmax), ja(nzmax)
      double precision a(nzmax)
      common /RMAT/ a, ia, ja, nz, n
c
c  arguments variables.
c
      integer m, mmax, lwsp, liwsp
      parameter( mmax = 50 )
      parameter( lwsp = nmax*(mmax+2)+5*(mmax+2)**2+7, liwsp = nmax )
      integer iwsp(liwsp)
      double precision t, tol, anorm, tmp
      double precision v(nmax), w(nmax), wsp(lwsp)

      integer i, j, itrace, iflag
      double precision ZERO, ONE
      parameter( ZERO=0.0d0, ONE=1.0d0 )

      intrinsic ABS
c
c  Executable statements.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SAMPLE_B'
c
c  Load the matrix.
c
      n = nmax
      nz = nzmax
      call loadhb ( 'gr3030.txt', 'crs', n,nz,ia,ja,a, iwsp )
c
c  Compute the infinity norm of A.
c
      anorm = 0.0d0
      do i = 1,n
         tmp = 0.0d0
         do j = ia(i),ia(i+1)-1
            tmp = tmp + ABS( a(j) )
         end do
         if ( anorm.lt.tmp ) anorm = tmp
      end do
      write(*,FMT='(A,E8.2)') '||A||_inf= ',anorm
c
c  The operand vector v is set to (1,..., 1)'.
c
      do i = 1,n
         v(i) = ONE
      end do
c
c  Set other input arguments.
c
      t = 1.0d0
      tol = 1.0d-10
      m = 30
      itrace = 0
c
c  Compute exp(t*A)v with CRS format.
c
      tic = clock()
      call DSEXPV( n, m, t,v,w, tol, anorm,
     &             wsp,lwsp, iwsp,liwsp, dgcrsv, itrace, iflag )
      tac = clock()

      print 9001,'----------------------------------------------------'
      print 9001,'DSEXPV (Forward) has completed:'
      print 9001,'----------------------------------------------------'
      print 9001,'w(1:10) ='
      do i = 1,10
         print*,w(i)
      end do
c
c  Display some statistics if desired.
c
      print 9001,'final report----------------------------------------'
      print 9002,'runtime   = ',tac-tic
      print 9002,'||A||_inf = ',anorm
      print 9003,'nz        =',nz
      print 9003,'n         =',n
      print 9003,'m         =',m
      print 9003,'itrace    =',itrace
      print 9003,'iflag     =',iflag
      print 9003,'ibrkflag  =',iwsp(6)
      print 9003,'mbrkdwn   =',iwsp(7)
      print 9003,'nstep     =',iwsp(4)
      print 9003,'nreject   =',iwsp(5)
      print 9003,'nmult     =',iwsp(1)
      print 9003,'nexph     =',iwsp(2)
      print 9003,'nscale    =',iwsp(3)

      print 9002,'tol       = ',tol
      print 9002,'t         = ',t
      print 9002,'tbrkdwn   = ',wsp(7)
      print 9002,'step_min  = ',wsp(1)
      print 9002,'step_max  = ',wsp(2)
      print 9002,'max_round = ',wsp(3)
      print 9002,'sum_round = ',wsp(4)
      print 9002,'max_error = ',wsp(5)
      print 9002,'sum_error = ',wsp(6)
      print 9002,'hump      = ',wsp(9)
      print 9002,'scale-norm= ',wsp(10)
c
c  Compute exp(-t*A)w with CRS format.
c
      t = -t
      call DCOPY( n, w,1, v,1 )
      tic = clock()
      call DSEXPV( n, m, t,v,w, tol, anorm,
     &             wsp,lwsp, iwsp,liwsp, dgcrsv, itrace, iflag )
      tac = clock()

      print 9001,'----------------------------------------------------'
      print 9001,'DSEXPV (Backward) has completed:'
      print 9001,'----------------------------------------------------'
      print 9001,'w(1:10) ='
      do i = 1,10
         print*,w(i)
      end do
c
c  Display some statistics if desired.
c
      print 9001,'final report----------------------------------------'
      print 9002,'runtime   = ',tac-tic
      print 9002,'||A||_inf = ',anorm
      print 9003,'nz        =',nz
      print 9003,'n         =',n
      print 9003,'m         =',m
      print 9003,'itrace    =',itrace
      print 9003,'iflag     =',iflag
      print 9003,'ibrkflag  =',iwsp(6)
      print 9003,'mbrkdwn   =',iwsp(7)
      print 9003,'nstep     =',iwsp(4)
      print 9003,'nreject   =',iwsp(5)
      print 9003,'nmult     =',iwsp(1)
      print 9003,'nexph     =',iwsp(2)
      print 9003,'nscale    =',iwsp(3)

      print 9002,'tol       = ',tol
      print 9002,'t         = ',t
      print 9002,'tbrkdwn   = ',wsp(7)
      print 9002,'step_min  = ',wsp(1)
      print 9002,'step_max  = ',wsp(2)
      print 9002,'max_round = ',wsp(3)
      print 9002,'sum_round = ',wsp(4)
      print 9002,'max_error = ',wsp(5)
      print 9002,'sum_error = ',wsp(6)
      print 9002,'hump      = ',wsp(9)
      print 9002,'scale-norm= ',wsp(10)

 9001 format(A)
 9002 format(A,E8.2)
 9003 format(A,I9)
      END
      subroutine sample_d ( )

c*********************************************************************72
c
cc SAMPLE_D demonstrates EXPOKIT on small dense matrices.
c
c  Discussion:
c
c    Sample program illustrating the computation of small matrix
c    exponentials in full with Expokit.  Refer to the Expokit
c    documentation for more details about the methods, and 
c    especially the domain of applicability of the Chebyshev scheme.
c
      implicit none

      integer m,i,j,k,ideg,mprint,lda,ldh,lwsp,liwsp,ns,iflag,iexp
      parameter ( ideg=6, lda=50, ldh=lda )
      parameter ( lwsp=4*ldh*ldh+ideg+1, liwsp=ldh )

      integer iwsp(liwsp), iseed(4)
      double precision t, A(lda,lda), H(ldh,ldh), y(ldh)
      double precision wsp(lwsp), s1, s2
      complex*16 Hc(ldh,ldh), yc(ldh), wspc(lwsp)

      double precision DLARAN
      intrinsic CMPLX, CONJG, MIN, DBLE, IMAG

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SAMPLE_D'
c
c  REAL CASE
c
c  set A = random symmetric negative define matrix.
c
      t = 2.0d0
      m = 5
      iseed(1) = 3
      iseed(2) = 7
      iseed(3) = 3
      iseed(4) = 7
      do j = 1,m
         do i = j,m
            A(i,j) = DLARAN( iseed )
            A(j,i) = A(i,j)
         end do
         A(j,j) = -2.5d0 + A(j,j)
      end do
c
c  maximum number of rows/columns to be printed.
c
      mprint = MIN(5,m)

      print*,"t ="
      print*,t
      print 9000,"REAL SYMMETRIC CASE","*******************************"
      print*,"A = "
      print 9001,( (A(i,j), j=1,mprint), i=1,mprint )

 9000 format( /,A,/,A )
 9001 format( 5(1X,D11.4) )
c
c  Pade.
c
      call DGPADM( ideg, m, t, A,lda, wsp,lwsp, iwsp,iexp,ns, iflag )
      print 9000,"With DGPADM:","exp(t*A) ="
      print 9001,( (wsp(iexp+(j-1)*m+i-1), j=1,mprint), i=1,mprint )

      call DSPADM( ideg, m, t, A,lda, wsp,lwsp, iwsp,iexp,ns, iflag )
      print 9000,"With DSPADM:","exp(t*A) ="
      print 9001,( (wsp(iexp+(j-1)*m+i-1), j=1,mprint), i=1,mprint )
c
c  Chebyshev.
c
      do i = 1,m
         y(i) = 0.0d0
      end do
      y(1) = 1.0d0
      call DGCHBV( m,t , A,lda, y, wsp, iwsp, iflag )
      print 9000,"With DGCHBV:","exp(t*A)e_1 ="
      do i = 1,mprint
         print*,y(i)
      end do

      do i = 1,m
         y(i) = 0.0d0
      end do
      y(1) = 1.0d0
      call DSCHBV( m,t , A,lda, y, wsp, iwsp, iflag )
      print 9000,"With DSCHBV:","exp(t*A)e_1 ="
      do i = 1,mprint
         print*,y(i)
      end do
c
c  set H = upper Hessenberg part of A.
c
      do j = 1,m
         do i = 1,MIN(j+1,m)
            H(i,j) = A(i,j)
         end do
         do k = i,m
            H(k,j) = 0.0d0
         end do
      end do

      print 9000,"REAL UPPER HESSENBERG CASE","************************"
      print*,"H ="
      print 9001,( (H(i,j), j=1,mprint), i=1,mprint )
c
c  Pade.
c
      call DGPADM( ideg, m, t, H,ldh, wsp,lwsp, iwsp,iexp,ns, iflag )
      print 9000,"With DGPADM:","exp(t*H) ="
      print 9001,( (wsp(iexp+(j-1)*m+i-1), j=1,mprint), i=1,mprint )
c
c  Chebyshev
c
      do i = 1,m
         y(i) = 0.0d0
      end do
      y(1) = 1.0d0
      call DNCHBV( m,t, A,lda, y, wsp, iflag )
      print 9000,"With DNCHBV:","exp(t*A)e_1 ="
      do i = 1,mprint
         print*,y(i)
      end do
c
c  COMPLEX CASE
c
c  generate the diagonal.
c
      do i = 1,m
         s1 = DLARAN( iseed ) - 2.0d0
         Hc(i,i) = s1
      end do
c
c  generate the lower part.
c
      do j = 1,m
         do i = j+1,m
            s1 = DLARAN( iseed ) - 0.5d0
            s2 = DLARAN( iseed ) - 0.5d0
            Hc(i,j) = CMPLX( s1,s2 )
         end do
      end do
c
c  include the conjugate upper part explicitly.
c
      do j = 1,m
         do i = j+1,m
            Hc(j,i) = CONJG( Hc(i,j) )
         end do
      end do

      print 9000,"COMPLEX HERMITIAN CASE","****************************"
      print 9000," ","Re(H) ="
      print 9001,( (DBLE(Hc(i,j)), j=1,mprint), i=1,mprint )
      print 9000," ","Im(H) ="
      print 9001,( (IMAG(Hc(i,j)), j=1,mprint), i=1,mprint )

      call ZGPADM( ideg, m,t, Hc,ldh, wspc,lwsp,iwsp, iexp, ns,iflag )
      print 9000,"With ZGPADM:","exp(t*H)e_1 ="
      do i = 1,mprint
         print*,wspc(iexp+i-1)
      end do

      call ZHPADM( ideg, m,t, Hc,ldh, wspc,lwsp,iwsp, iexp, ns,iflag )
      print 9000,"With ZHPADM:","exp(t*H)_e_1 ="
      do i = 1,mprint
         print*,wspc(iexp+i-1)
      end do
c
c  Chebyshev
c
      do i = 1,m
         yc(i) = 0.0d0
      end do
      yc(1) = 1.0d0
      call ZGCHBV( m,t, Hc,ldh, yc, wsp, iwsp, iflag )
      print 9000,"With ZGCHBV:","exp(t*H)e_1 ="
      do i = 1,mprint
         print*,yc(i)
      end do
c
c  Note: the Hermitian feature is not useful vis-a-vis Chebyshev
c
      END
      subroutine sample_g ( )

c*********************************************************************72
c
cc SAMPLE_G illustrates the use of DGEXPV.
c
c  Discussion:
c
c    Non-symmetric problem (Example 6.3 in the Expokit report).
c
      implicit none

      external dgcoov, dgcrsv, dgccsv

      double precision tic, tac, clock
c
c  matrix data. 
c  BEWARE: these values should match those in dgmatv.f
c
      integer n, nz, nmax, nzmax
      parameter( nmax = 5000, nzmax = 600000 )
      integer ia(nzmax), ja(nzmax)
      double precision a(nzmax)
      common /RMAT/ a, ia, ja, nz, n
c
c  arguments variables.
c
      integer m, mmax, lwsp, liwsp
      parameter( mmax = 50 )
      parameter( lwsp = nmax*(mmax+2)+5*(mmax+2)**2+7, liwsp = nmax )
      integer iwsp(liwsp)
      double precision t, tol, anorm
      double precision v(nmax), w(nmax), wsp(lwsp)

      integer i, itrace, iflag
      double precision ZERO, ONE
      parameter( ZERO=0.0d0, ONE=1.0d0 )

      intrinsic ABS
c
c  Executable statements.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SAMPLE_G'
c
c  load a Harwell-Boeing matrix.
c
      n = nmax
      nz = nzmax
      call loadhb( 'orani678.txt', 'coo', n,nz,ia,ja,a, iwsp )
c
c  compute the infinity norm of A.
c
      do i = 1,n
         wsp(i) = ZERO
      end do
      do i = 1,nz
         wsp(ia(i)) = wsp(ia(i)) + ABS( a(i) )
      end do
      anorm = wsp(1)
      do i = 2,n
         if ( anorm.lt.wsp(i) ) anorm =  wsp(i)
      end do
      write(UNIT=*,FMT='(A,E8.2)') '||A||_inf= ',anorm
c
c  the operand vector v is set to (1, ..., 1)^T.
c
      do i = 1,n
         v(i) = ONE
      end do
c
c  set other input arguments.
c
      t = 10.0d0
      tol = 0.0d0
      m = 30
      itrace = 0
c
c  compute w = exp(t*A)v with COO format.
c
      tic = clock()
      call DGEXPV( n, m, t,v,w, tol, anorm,
     &             wsp,lwsp, iwsp,liwsp, dgcoov, itrace, iflag ) 
      tac = clock()

      print 9001,'----------------------------------------------------'
      print 9001,'DGEXPV (COO) has completed:'
      print 9001,'----------------------------------------------------'
      print 9001,'w(1:10) ='
      do i = 1,10
         print*,w(i)
      end do
c
c  display some statistics if desired.
c
      print 9001,'final report----------------------------------------'
      print 9002,'runtime   = ',tac-tic
      print 9002,'||A||_inf = ',anorm
      print 9003,'nz        =',nz
      print 9003,'n         =',n
      print 9003,'m         =',m
      print 9003,'itrace    =',itrace
      print 9003,'iflag     =',iflag
      print 9003,'ibrkflag  =',iwsp(6)
      print 9003,'mbrkdwn   =',iwsp(7)
      print 9003,'nstep     =',iwsp(4)
      print 9003,'nreject   =',iwsp(5)
      print 9003,'nmult     =',iwsp(1)
      print 9003,'nexph     =',iwsp(2)
      print 9003,'nscale    =',iwsp(3)

      print 9002,'tol       = ',tol
      print 9002,'t         = ',t
      print 9002,'tbrkdwn   = ',wsp(7)
      print 9002,'step_min  = ',wsp(1)
      print 9002,'step_max  = ',wsp(2)
      print 9002,'max_round = ',wsp(3)
      print 9002,'sum_round = ',wsp(4)
      print 9002,'max_error = ',wsp(5)
      print 9002,'sum_error = ',wsp(6)
      print 9002,'hump      = ',wsp(9)
      print 9002,'scale-norm= ',wsp(10)
c
c  convert from COO to CCS.
c
      call dgcnvr( 'coo','ccs','n', n,n, nz, ia, ja, a, iwsp )
c
c  compute w = exp(t*A)v with CCS format.
c
      tic = clock()
      call DGEXPV( n, m, t,v,w, tol, anorm,
     &             wsp,lwsp, iwsp,liwsp, dgccsv, itrace, iflag ) 
      tac = clock()

      print 9001,'----------------------------------------------------'
      print 9001,'DGEXPV (CCS) has completed:'
      print 9001,'----------------------------------------------------'
      print 9001,'w(1:10) ='
      do i = 1,10
         print*,w(i)
      end do
c
c  display some statistics if desired.
c
      print 9001,'final report----------------------------------------'
      print 9002,'runtime   = ',tac-tic
      print 9002,'||A||_inf = ',anorm
      print 9003,'nz        =',nz
      print 9003,'n         =',n
      print 9003,'m         =',m
      print 9003,'itrace    =',itrace
      print 9003,'iflag     =',iflag
      print 9003,'ibrkflag  =',iwsp(6)
      print 9003,'mbrkdwn   =',iwsp(7)
      print 9003,'nstep     =',iwsp(4)
      print 9003,'nreject   =',iwsp(5)
      print 9003,'nmult     =',iwsp(1)
      print 9003,'nexph     =',iwsp(2)
      print 9003,'nscale    =',iwsp(3)

      print 9002,'tol       = ',tol
      print 9002,'t         = ',t
      print 9002,'tbrkdwn   = ',wsp(7)
      print 9002,'step_min  = ',wsp(1)
      print 9002,'step_max  = ',wsp(2)
      print 9002,'max_round = ',wsp(3)
      print 9002,'sum_round = ',wsp(4)
      print 9002,'max_error = ',wsp(5)
      print 9002,'sum_error = ',wsp(6)
      print 9002,'hump      = ',wsp(9)
      print 9002,'scale-norm= ',wsp(10)
c
c  convert from CCS to CRS.
c
      call dgcnvr( 'ccs','crs','n', n,n, nz, ia, ja, a, iwsp )
c
c  compute w = exp(t*A)v with CRS format.
c
      tic = clock()
      call DGEXPV( n, m, t,v,w, tol, anorm,
     &             wsp,lwsp, iwsp,liwsp, dgcrsv, itrace, iflag ) 
      tac = clock()

      print 9001,'----------------------------------------------------'
      print 9001,'DGEXPV (CRS) has completed:'
      print 9001,'----------------------------------------------------'
      print 9001,'w(1:10) ='
      do i = 1,10
         print*,w(i)
      end do
c
c  display some statistics if desired.
c
      print 9001,'final report----------------------------------------'
      print 9002,'runtime   = ',tac-tic
      print 9002,'||A||_inf = ',anorm
      print 9003,'nz        =',nz
      print 9003,'n         =',n
      print 9003,'m         =',m
      print 9003,'itrace    =',itrace
      print 9003,'iflag     =',iflag
      print 9003,'ibrkflag  =',iwsp(6)
      print 9003,'mbrkdwn   =',iwsp(7)
      print 9003,'nstep     =',iwsp(4)
      print 9003,'nreject   =',iwsp(5)
      print 9003,'nmult     =',iwsp(1)
      print 9003,'nexph     =',iwsp(2)
      print 9003,'nscale    =',iwsp(3)

      print 9002,'tol       = ',tol
      print 9002,'t         = ',t
      print 9002,'tbrkdwn   = ',wsp(7)
      print 9002,'step_min  = ',wsp(1)
      print 9002,'step_max  = ',wsp(2)
      print 9002,'max_round = ',wsp(3)
      print 9002,'sum_round = ',wsp(4)
      print 9002,'max_error = ',wsp(5)
      print 9002,'sum_error = ',wsp(6)
      print 9002,'hump      = ',wsp(9)
      print 9002,'scale-norm= ',wsp(10)

 9001 format(A)
 9002 format(A,E8.2)
 9003 format(A,I9)
      END
      subroutine sample_m ( )

c*********************************************************************72
c
cc SAMPLE_M illustrates the use of DMEXPV and DGEXPV.
c
c  Discussion:
c
c    Binary Markov Model (Example 6.1 in the Expokit report).
c
      implicit none

      external dgcoov, dgcrsv, dgccsv

      double precision tic, tac, clock
c
c  matrix data.
c  BEWARE: these values must match those in dgmatv.f
c
      integer n, nz, nmax, nzmax
      parameter( nmax=5000, nzmax = 600000 )
      integer ia(nzmax), ja(nzmax)
      double precision a(nzmax)
      common /RMAT/ a, ia, ja, nz, n
c
c  arguments variables.
c
      integer m, mmax, lwsp, liwsp
      parameter( mmax = 50 )
      parameter( lwsp = nmax*(mmax+2)+5*(mmax+2)**2+7, liwsp = nmax )
      integer iwsp(liwsp)
      double precision t, tol, anorm
      double precision v(nmax), w(nmax), wsp(lwsp)

      integer i, j, itrace, iflag
      double precision ZERO, ONE, tmp
      parameter( ZERO=0.0d0, ONE=1.0d0 )

      intrinsic ABS
c
c  Executable statements.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SAMPLE_M'
c
c  load the infinitesimal generator (CRS format)
c
      open( UNIT=7,STATUS='old',IOSTAT=iflag,FILE='c1024.txt')
      if ( iflag.ne.0 ) stop 'Error - matrix could not be loaded.'
      read( UNIT=7,FMT=* ) n, nz
      if ( nz.gt.nzmax ) stop 'Please increase nzmax.'
      if ( n.gt.nmax ) stop 'Please increase nmax.'
      read( UNIT=7,FMT=* ) (ia(i), i=1,n+1)
      read( UNIT=7,FMT=* ) (ja(i), a(i), i=1,nz)
      close( UNIT=7 )
c
c  make sure the infinitesimal generator is transposed,
c     this encoded check prevents from falling in the famous
c     (or rather infamous) `transpose trap' !
c
      do j = 1,2*n
         wsp(j) = 0.0d0
      end do
      do i = 1,n
         do j = ia(i),ia(i+1)-1
            wsp(i) = wsp(i) + a(j)
            wsp(n+ja(j)) = wsp(n+ja(j)) + a(j)
         end do
      end do
      wsp(1) = ABS( wsp(1) )
      wsp(n+1) = ABS( wsp(n+1) )
      do i = 2,n
         wsp(1) = wsp(1) + ABS( wsp(i) )
         wsp(n+1) = wsp(n+1) + ABS( wsp(n+i) )
      end do
      if ( wsp(n+1).gt.wsp(1) ) then
         print*,'Transposing the input matrix... '
         call tnspos( n, nz, ia, ja, a, iwsp )
      endif
c
c  compute the infinity norm of A.
c
      anorm = 0.0d0
      do i = 1,n
         tmp = 0.0d0
         do j = ia(i),ia(i+1)-1
            tmp = tmp + ABS( a(j) )
         end do
         if ( anorm.lt.tmp ) anorm = tmp
      end do
      write(UNIT=*,FMT='(A,E8.2)') '||A||_inf= ',anorm
c
c  the operand vector v is set to the first unit basis vector.
c
      v(1) = ONE
      do i = 2,n
         v(i) = ZERO
      end do
c
c  set other input arguments.
c
      t = 10.0d0
      tol = 1.0d-10
      m = 30
      itrace = 0
c
c  compute w = exp(t*A)v with DMEXPV.
c
      tic = clock()
      call DMEXPV( n, m, t,v,w, tol, anorm,
     &             wsp,lwsp, iwsp,liwsp, dgcrsv, itrace, iflag )
      tac = clock()

      print 9001,'----------------------------------------------------'
      print 9001,'DMEXPV has completed:'
      print 9001,'----------------------------------------------------'
      print 9001,'w(1:10) ='
      do i = 1,10
         print*,w(i)
      end do
c
c  display some statistics if desired.
c
      print 9001,'final report----------------------------------------'
      print 9002,'runtime   = ',tac-tic
      print 9002,'||A||_inf = ',anorm
      print 9003,'nz        =',nz
      print 9003,'n         =',n
      print 9003,'m         =',m
      print 9003,'itrace    =',itrace
      print 9003,'iflag     =',iflag
      print 9003,'ibrkflag  =',iwsp(6)
      print 9003,'mbrkdwn   =',iwsp(7)
      print 9003,'nstep     =',iwsp(4)
      print 9003,'nreject   =',iwsp(5)
      print 9003,'nmult     =',iwsp(1)
      print 9003,'nexph     =',iwsp(2)
      print 9003,'nscale    =',iwsp(3)

      print 9002,'tol       = ',tol
      print 9002,'t         = ',t
      print 9002,'tbrkdwn   = ',wsp(7)
      print 9002,'step_min  = ',wsp(1)
      print 9002,'step_max  = ',wsp(2)
      print 9002,'max_round = ',wsp(3)
      print 9002,'sum_round = ',wsp(4)
      print 9002,'max_error = ',wsp(5)
      print 9002,'sum_error = ',wsp(6)
      print 9002,'hump      = ',wsp(9)
      print 9002,'scale-norm= ',wsp(10)
c
c  compute w = exp(t*A)v with DGEXPV.
c
      tic = clock()
      call DGEXPV( n, m, t,v,w, tol, anorm,
     &             wsp,lwsp, iwsp,liwsp, dgcrsv, itrace, iflag )
      tac = clock()

      print 9001,'----------------------------------------------------'
      print 9001,'DGEXPV has completed:'
      print 9001,'----------------------------------------------------'
      print 9001,'w(1:10) ='
      do i = 1,10
         print*,w(i)
      end do
c
c  display some statistics if desired.
c
      print 9001,'final report----------------------------------------'
      print 9002,'runtime   = ',tac-tic
      print 9002,'||A||_inf = ',anorm
      print 9003,'nz        =',nz
      print 9003,'n         =',n
      print 9003,'m         =',m
      print 9003,'itrace    =',itrace
      print 9003,'iflag     =',iflag
      print 9003,'ibrkflag  =',iwsp(6)
      print 9003,'mbrkdwn   =',iwsp(7)
      print 9003,'nstep     =',iwsp(4)
      print 9003,'nreject   =',iwsp(5)
      print 9003,'nmult     =',iwsp(1)
      print 9003,'nexph     =',iwsp(2)
      print 9003,'nscale    =',iwsp(3)

      print 9002,'tol       = ',tol
      print 9002,'t         = ',t
      print 9002,'tbrkdwn   = ',wsp(7)
      print 9002,'step_min  = ',wsp(1)
      print 9002,'step_max  = ',wsp(2)
      print 9002,'max_round = ',wsp(3)
      print 9002,'sum_round = ',wsp(4)
      print 9002,'max_error = ',wsp(5)
      print 9002,'sum_error = ',wsp(6)
      print 9002,'hump      = ',wsp(9)
      print 9002,'scale-norm= ',wsp(10)

 9001 format(A)
 9002 format(A,E8.2)
 9003 format(A,I9)
      END
      subroutine sample_p ( )

c*********************************************************************72
c
cc SAMPLE_P illustrates the use of DGPHIV.
c
c  Discussion:
c
c    Nonhomogeneous problem (Example 6.5 in the Expokit report).
c
      implicit none

      external dgcoov, dgcrsv, dgccsv

      double precision tic, tac, clock
c
c  matrix data.
c  BEWARE: these values must match those in dgmatv.f
c
      integer n, nz, nmax, nzmax
      parameter( nmax = 5000, nzmax = 600000 )
      integer ia(nzmax), ja(nzmax)
      double precision a(nzmax)
      common /RMAT/ a, ia, ja, nz, n
c
c  arguments variables.
c
      integer m, mmax, lwsp, liwsp
      parameter( mmax = 50 )
      parameter( lwsp = nmax*(mmax+3)+5*(mmax+3)**2+7, liwsp = nmax )
      integer iwsp(liwsp)
      double precision t, tol, anorm, tmp
      double precision u(nmax),v(nmax),w(nmax), wsave(nmax), wsp(lwsp)

      integer i, itrace, iflag
      double precision ZERO, ONE
      parameter( ZERO=0.0d0, ONE=1.0d0 )

      intrinsic ABS
      double precision DNRM2
c
c  Executable statements.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SAMPLE_P'
c
c  load a Harwell-Boeing matrix.
c
      n = nmax
      nz = nzmax
      call loadhb( 'orani678.txt', 'coo', n,nz,ia,ja,a, iwsp )
c
c  compute the infinity norm of A.
c
      do i = 1,n
         wsp(i) = ZERO
      end do
      do i = 1,nz
         wsp(ia(i)) = wsp(ia(i)) + ABS( a(i) )
      end do
      anorm = wsp(1)
      do i = 2,n
         if ( anorm.lt.wsp(i) ) anorm =  wsp(i)
      end do
      write(UNIT=*,FMT='(A,E8.2)') '||A||_inf= ',anorm
c
c  back to CCS format.
c
      call dgcnvr( 'coo','ccs','n', n,n, nz, ia, ja, a, iwsp )
c
c  set other input arguments.
c
      t = 10.0d0
      tol = 0.0d0
      m = 30
      itrace = 0
c
c  First Run
c
c  the operand vector u is set to zero.
c
      do i = 1,n
         u(i) = ZERO
      end do
c
c  the operand vector v is set to (1, ..., 1)'.
c
      do i = 1,n
         v(i) = ONE
      end do
c
c  compute w = exp(t*A)*v with DGPHIV.
c
      tic = clock()
      call DGPHIV( n, m, t,u,v,w, tol, anorm,
     &             wsp,lwsp, iwsp,liwsp, dgccsv, itrace, iflag ) 
      tac = clock()

      print 9001,'----------------------------------------------------'
      print 9001,'DGPHIV (CCS) has completed:'
      print 9001,'----------------------------------------------------'
      print 9001,'w(1:10) ='
      do i = 1,10
         write(*,*) w(i)
      end do
      call DCOPY( n, w,1, wsave,1 )
c
c  display some statistics if desired.
c
      print 9001,'final report----------------------------------------'
      print 9002,'runtime   = ',tac-tic
      print 9002,'||A||_inf = ',anorm
      print 9003,'nz        =',nz
      print 9003,'n         =',n
      print 9003,'m         =',m
      print 9003,'itrace    =',itrace
      print 9003,'iflag     =',iflag
      print 9003,'ibrkflag  =',iwsp(6)
      print 9003,'mbrkdwn   =',iwsp(7)
      print 9003,'nstep     =',iwsp(4)
      print 9003,'nreject   =',iwsp(5)
      print 9003,'nmult     =',iwsp(1)
      print 9003,'nexph     =',iwsp(2)
      print 9003,'nscale    =',iwsp(3)

      print 9002,'tol       = ',tol
      print 9002,'t         = ',t
      print 9002,'tbrkdwn   = ',wsp(7)
      print 9002,'step_min  = ',wsp(1)
      print 9002,'step_max  = ',wsp(2)
      print 9002,'max_round = ',wsp(3)
      print 9002,'sum_round = ',wsp(4)
      print 9002,'max_error = ',wsp(5)
      print 9002,'sum_error = ',wsp(6)
c
c  Second Run
c
c  the operand vector u is set to (1,..., 1)'.
c
      do i = 1,n
         u(i) = ONE
      end do
c
c  the operand vector v is set to zero.
c
      do i = 1,n
         v(i) = ZERO
      end do
c
c  compute w = t*phi(t*A)*u with DGPHIV.
c
      tic = clock()
      call DGPHIV( n, m, t,u,v,w, tol, anorm,
     &             wsp,lwsp, iwsp,liwsp, dgccsv, itrace, iflag ) 
      tac = clock()

      print 9001,'----------------------------------------------------'
      print 9001,'DGPHIV (CCS) has completed:'
      print 9001,'----------------------------------------------------'
      print 9001,'w(1:10) ='
      do i = 1,10
         write(*,*) w(i)
      end do
c
c  display some statistics if desired.
c
      print 9001,'final report----------------------------------------'
      print 9002,'runtime   = ',tac-tic
      print 9002,'||A||_inf = ',anorm
      print 9003,'nz        =',nz
      print 9003,'n         =',n
      print 9003,'m         =',m
      print 9003,'itrace    =',itrace
      print 9003,'iflag     =',iflag
      print 9003,'ibrkflag  =',iwsp(6)
      print 9003,'mbrkdwn   =',iwsp(7)
      print 9003,'nstep     =',iwsp(4)
      print 9003,'nreject   =',iwsp(5)
      print 9003,'nmult     =',iwsp(1)
      print 9003,'nexph     =',iwsp(2)
      print 9003,'nscale    =',iwsp(3)

      print 9002,'tol       = ',tol
      print 9002,'t         = ',t
      print 9002,'tbrkdwn   = ',wsp(7)
      print 9002,'step_min  = ',wsp(1)
      print 9002,'step_max  = ',wsp(2)
      print 9002,'max_round = ',wsp(3)
      print 9002,'sum_round = ',wsp(4)
      print 9002,'max_error = ',wsp(5)
      print 9002,'sum_error = ',wsp(6)
c
c  now check ||(u + A * w) - wsave||/||wsave||. Due to the specific
c  previous settings, the answer should agree with tol.
c
      call dgccsv( w, v )
      call DAXPY( n, 1.0d0, u,1, v,1 )
      call DAXPY( n, -1.0d0, wsave,1, v,1 )
      tmp = DNRM2( n, v,1 ) / DNRM2( n, wsave, 1 )
      print*
      print*,"relative difference (phi vs. exp) =", tmp
      print*
c
c  Third Run.
c
c  the operand vector v is set to (1, ..., 1)'.
c
      do i = 1,n
         v(i) = ONE
      end do
c
c  the operand vector u is set to (1,..., 1)'.
c
      do i = 1,n
         v(i) = ONE
      end do
c
c  compute w = exp(t*A)*v + t*phi(t*A)*u with DGPHIV.
c
      tic = clock()
      call DGPHIV( n, m, t,u,v,w, tol, anorm,
     &             wsp,lwsp, iwsp,liwsp, dgccsv, itrace, iflag ) 
      tac = clock()

      print 9001,'----------------------------------------------------'
      print 9001,'DGPHIV (CCS) has completed:'
      print 9001,'----------------------------------------------------'
      print 9001,'w(1:10) ='
      do i = 1,10
         print*,w(i)
      end do
c
c  display some statistics if desired.
c
      print 9001,'final report----------------------------------------'
      print 9002,'runtime   = ',tac-tic
      print 9002,'||A||_inf = ',anorm
      print 9003,'nz        =',nz
      print 9003,'n         =',n
      print 9003,'m         =',m
      print 9003,'itrace    =',itrace
      print 9003,'iflag     =',iflag
      print 9003,'ibrkflag  =',iwsp(6)
      print 9003,'mbrkdwn   =',iwsp(7)
      print 9003,'nstep     =',iwsp(4)
      print 9003,'nreject   =',iwsp(5)
      print 9003,'nmult     =',iwsp(1)
      print 9003,'nexph     =',iwsp(2)
      print 9003,'nscale    =',iwsp(3)

      print 9002,'tol       = ',tol
      print 9002,'t         = ',t
      print 9002,'tbrkdwn   = ',wsp(7)
      print 9002,'step_min  = ',wsp(1)
      print 9002,'step_max  = ',wsp(2)
      print 9002,'max_round = ',wsp(3)
      print 9002,'sum_round = ',wsp(4)
      print 9002,'max_error = ',wsp(5)
      print 9002,'sum_error = ',wsp(6)

 9001 format(A)
 9002 format(A,E8.2)
 9003 format(A,I9)
      END
      subroutine sample_z ( )

c*********************************************************************72
c
cc SAMPLE_Z illustrates the use of ZGEXPV and ZHEXPV.
c
c  Discussion:
c
c    Hermitian problem (Example 6.2 in the Expokit report).
c
      implicit none

      external zgcoov, zgcrsv, zgccsv

      double precision tic, tac, clock
c
c  matrix data.
c  BEWARE: these values must match those in zgmatv.f
c
      integer n, nz, nmax, nzmax
      parameter( nmax = 5500, nzmax = 50000 )
      integer ia(nzmax), ja(nzmax)
      complex*16 a(nzmax)
      common /CMAT/ a, ia, ja, nz, n
c
c  arguments variables.
c
      integer m, mmax, lwsp, liwsp
      parameter( mmax = 50 )
      parameter( lwsp = nmax*(mmax+2)+5*(mmax+2)**2+7, liwsp = mmax+2 )
      integer iwsp(liwsp)
      double precision t, tol, anorm, s1, s2
      complex*16 v(nmax), w(nmax), wsp(lwsp)

      integer i, j, nnz, itrace, iflag, iseed(4)
      complex*16 ZERO, ONE
      parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )

      double precision DLARAN
      intrinsic ABS, CMPLX, CONJG, DBLE
c
c  Executable statements.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SAMPLE_Z'
c
c  Load a symmetric pattern.
c
      n = nmax
      nz = nzmax/2
      call getpat ( 'bcspwr10.txt', n, nz, ia, ja )
c
c  For the purpose of the experiments, expand to COOrdinates.
c
      nnz = nz
      do j = n,1,-1
         do i = 1,ja(j+1)-ja(j)
            ja(nnz) = j
            nnz = nnz-1
         end do
      end do
c
c  Fill-in an Hermitian matrix -- the conjugate part is included.
c
      iseed(1) = 0
      iseed(2) = 0
      iseed(3) = 0
      iseed(4) = 5
      nnz = nz
      do i = 1,nz
         if ( ia(i).ne.ja(i) ) then
            s1 = 10.0d0*DLARAN( iseed ) - 5.0d0
            s2 = 10.0d0*DLARAN( iseed ) - 5.0d0
            a(i) = CMPLX( s1,s2 )
            nnz = nnz + 1
            a(nnz) = CONJG( a(i) )
            ia(nnz) = ja(i)
            ja(nnz) = ia(i)
        else
            s1 = 10.0d0*DLARAN( iseed ) - 5.0d0
            a(i) = CMPLX( s1,0.0d0 )
         endif
      end do
      nz = nnz
c
c  Compute the infinity norm of A.
c
      do i = 1,n
         wsp(i) = ZERO
      end do
      do i = 1,nz
         wsp(ia(i)) = wsp(ia(i)) + ABS( a(i) )
      end do
      anorm = wsp(1)
      do i = 2,n
        if ( anorm.lt.DBLE(wsp(i)) ) then
          anorm =  wsp(i)
        end if
      end do
c
c  Convert from COO to CRS.
c
c      call zgcnvr( 'coo','crs','n', n,n, nz, ia, ja, a, iwsp )
c
c  Compute the infinity norm of A.
c
c      anorm = 0.0d0
c      do i = 1,n
c         s1 = 0.0d0
c         do j = ia(i),ia(i+1)-1
c            s1 = s1 + ABS( a(j) )
c         end do
c         if ( anorm.lt.tmp ) anorm = s1
c      end do

      write(UNIT=*,FMT='(A,E8.2)') '||A||_inf= ',anorm
c
c  The operand vector v is set to e_1 + e_n.
c
      do i = 1,n
         v(i) = ZERO
      end do
      v(1) = ONE
      v(n) = ONE
c
c  Set other input arguments.
c
      t = 1.0d0
      tol = 1.0d-5
      m = 30
      itrace = 0
c
c  Compute w = exp(t*A)v with ZGEXPV.
c
      tic = clock()
      call ZGEXPV( n, m, t,v,w, tol, anorm,
     &             wsp,lwsp, iwsp,liwsp, zgcoov, itrace, iflag )
      tac = clock()

      print 9001,'----------------------------------------------------'
      print 9001,'ZGEXPV has completed:'
      print 9001,'----------------------------------------------------'
      print 9001,'w(1:10) ='
      do i = 1,10
         print*,w(i)
      end do
c
c  Display some statistics if desired.
c
      print 9001,'final report----------------------------------------'
      print 9002,'runtime   = ',tac-tic
      print 9002,'||A||_inf = ',anorm
      print 9003,'nz        =',nz
      print 9003,'n         =',n
      print 9003,'m         =',m
      print 9003,'itrace    =',itrace
      print 9003,'iflag     =',iflag
      print 9003,'ibrkflag  =',iwsp(6)
      print 9003,'mbrkdwn   =',iwsp(7)
      print 9003,'nstep     =',iwsp(4)
      print 9003,'nreject   =',iwsp(5)
      print 9003,'nmult     =',iwsp(1)
      print 9003,'nexph     =',iwsp(2)
      print 9003,'nscale    =',iwsp(3)

      print 9002,'tol       = ',tol
      print 9002,'t         = ',t
      print 9002,'tbrkdwn   = ',DBLE( wsp(7) )
      print 9002,'step_min  = ',DBLE( wsp(1) )
      print 9002,'step_max  = ',DBLE( wsp(2) )
      print 9002,'max_round = ',DBLE( wsp(3) )
      print 9002,'sum_round = ',DBLE( wsp(4) )
      print 9002,'max_error = ',DBLE( wsp(5) )
      print 9002,'sum_error = ',DBLE( wsp(6) )
      print 9002,'hump      = ',DBLE( wsp(9) )
      print 9002,'scale-norm= ',DBLE( wsp(10) )
c
c  Compute w = exp(t*A)v with ZHEXPV.
c
      tic = clock()
      call ZHEXPV( n, m, t,v,w, tol, anorm,
     &             wsp,lwsp, iwsp,liwsp, zgcoov, itrace, iflag )
      tac = clock()

      print 9001,'----------------------------------------------------'
      print 9001,'ZHEXPV has completed:'
      print 9001,'----------------------------------------------------'
      print 9001,'w(1:10) ='
      do i = 1,10
         print*,w(i)
      end do
c
c  Display some statistics if desired.
c
      print 9001,'final report----------------------------------------'
      print 9002,'runtime   = ',tac-tic
      print 9002,'||A||_inf = ',anorm
      print 9003,'nz        =',nz
      print 9003,'n         =',n
      print 9003,'m         =',m
      print 9003,'itrace    =',itrace
      print 9003,'iflag     =',iflag
      print 9003,'ibrkflag  =',iwsp(6)
      print 9003,'mbrkdwn   =',iwsp(7)
      print 9003,'nstep     =',iwsp(4)
      print 9003,'nreject   =',iwsp(5)
      print 9003,'nmult     =',iwsp(1)
      print 9003,'nexph     =',iwsp(2)
      print 9003,'nscale    =',iwsp(3)

      print 9002,'tol       = ',tol
      print 9002,'t         = ',t
      print 9002,'tbrkdwn   = ',DBLE( wsp(7) )
      print 9002,'step_min  = ',DBLE( wsp(1) )
      print 9002,'step_max  = ',DBLE( wsp(2) )
      print 9002,'max_round = ',DBLE( wsp(3) )
      print 9002,'sum_round = ',DBLE( wsp(4) )
      print 9002,'max_error = ',DBLE( wsp(5) )
      print 9002,'sum_error = ',DBLE( wsp(6) )
      print 9002,'hump      = ',DBLE( wsp(9) )
      print 9002,'scale-norm= ',DBLE( wsp(10) )

 9001 format(A)
 9002 format(A,E8.2)
 9003 format(A,I9)
      END
