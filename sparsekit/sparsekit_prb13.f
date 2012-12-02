      program main

c*********************************************************************72
c
cc MAIN is the test program for ilut preconditioned gmres.
c
c  Discussion:
c
c    This program generates a sparse matrix using
c    matgen and then solves a linear system with an
c    artificial right hand side.
c
      implicit none

      integer nmax
      parameter ( nmax = 10 * 10 + 1 )
c
c  NZMAX would normally be 5*NMAX, reflecting the fact that 5 entries
c  are nonzero in each row of a matrix derived from the Laplacian.
c  However, ILUT requires that we allow a certain amount of fillin
c  during a partial factorization.
c
      integer nzmax
      parameter ( nzmax = 19 * nmax )

      double precision a(nzmax)
      double precision au(nzmax)
      double precision alpha
      double precision amax
      double precision eps
      double precision gammax
      double precision gammay
      integer i
      integer ia(nmax)
      integer ierr
      integer im
      parameter ( im = 10 )
      integer iout
      parameter ( iout = 6 )
      integer iw(nmax,3)
      integer j
      integer ja(nzmax)
      integer jau(nzmax)
      integer ju(nmax)
      integer k
      integer lfil
      integer maxits
      parameter ( maxits = 100 )
      integer meth
      integer n
      integer nw k
      integer nx
      parameter ( nx = 10 )
      integer ny
      parameter ( ny = 10 )
      integer nz 
      parameter ( nz = 1 )
      double precision random
      double precision tol
      double precision vv(nmax,20)
      double precision x(nmax)
      double precision xran(nmax)
      double precision y(nmax)

      common /func/ gammax, gammay, alpha

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPARSEKIT_PRB13'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Test the preconditioners and iterative'
      write ( *, '(a)' ) 'solvers in SPARSKIT.'
c
c  The PDE to be discretized is :
c
c -Lap u + gammax exp (xy)delx u + gammay exp (-xy) dely u +alpha u
c
c where Lap = 2-D laplacean, delx = part. der. wrt x,
c dely = part. der. wrt y.
c gammax, gammay, and alpha are passed via the commun func.
c
c data for PDE:
c
      alpha = -60.0
      gammax = 10.0
      gammay = 10.0
c
c  data for preconditioner
c
      nwk = nzmax
c
c  data for pgmres
c
      eps = 1.0E-07
c
c  same initial guess for gmres
c
      do j = 1, nmax
        xran(j) = random()
      end do
c
c  call gen57 to generate matrix in compressed sparse row format
c
      call gen57pt ( nx, ny, nz, a, ja, ia, ju, x )
c
c  define N.
c
      n = nx * ny * nz
c
c test all different methods:
c ILU0, MILU0, ILUT and with different values of tol and lfil
c ( from cheaper to more expensive preconditioners)
c The more accurate the preconditioner the fewer iterations
c are required in pgmres, in general.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Methods 1-5 '
      write ( *, '(a)' ) ' '

      do meth = 1, 5

        if ( meth == 1 ) then

          write ( *, '(a)' ) '1. Try ILU(0) Preconditioner'
          call ilu0 (n, a, ja, ia, au, jau, ju, iw, ierr)

        else if ( meth == 2 ) then

          write ( *, '(a)' ) '2. Try MILU(0) Preconditioner'
          call milu0 (n, a, ja, ia, au, jau, ju, iw, ierr)

        else if ( meth == 3 ) then

          write ( *, '(a)' ) '3. Try ILUT Preconditioner'
          write ( *, '(a)' ) 'with tol = 0.001, lfil=1.'
          tol  = 0.001
          lfil = 1

          call ilut (n,a,ja,ia,lfil,tol,au,jau,ju,nwk, 
     &      vv,vv(1,2),iw,iw(1,2),iw(1,3),ierr)

        else if ( meth == 4 ) then

          write ( *, '(a)' ) '4. Try ILUT Preconditioner'
          write ( *, '(a)' ) 'with tol = 0.001, lfil=5.'
          tol = 0.001
          lfil = 5

          call ilut (n,a,ja,ia,lfil,tol,au,jau,ju,nwk, 
     &      vv,vv(1,2),iw,iw(1,2),iw(1,3),ierr)

        else if ( meth == 5 ) then

          write ( *, '(a)' ) '5. Try ILUT Preconditioner'
          write ( *, '(a)' ) 'with tol = .0001, lfil=7.'
          tol = 0.0001
          lfil = 7

          call ilut (n,a,ja,ia,lfil,tol,au,jau,ju,nwk, 
     &      vv,vv(1,2),iw,iw(1,2),iw(1,3),ierr)

        end if
c
c  Check that return was succesful
c
        write ( *, * ) ' Precon set-up returned with ierr ', ierr

        if ( ierr /= 0 ) then
          continue
        end if
c
c  Generate right hand side = A * (1,2,3,...n)**T
c
        do k = 1, n
          x(k) = real ( k, kind = 8 )
        end do

        call ope ( n, x, y, a, ja, ia )
c
c  Generate initial guess.
c
        do j = 1, n
          x(j) = xran(j)
        end do

        call pgmres (n, im, y, x, vv, eps, maxits, iout, 
     &    a, ja, ia, au, jau, ju, ierr)

        write ( *, '(a,i6)' ) ' pgmres returned with ierr = ',ierr

        amax = 0.0
        do i = 1, n
          amax = max ( amax, abs ( x(i) - dble ( i ) ) )
        end do

        write ( *, '(a)' ) ' '
        write (*,*) 'Maximum error in solution = ', amax

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPARSEKIT_PRB13'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      function afun (x,y,z)

c*********************************************************************72
c
      implicit none

      double precision afun, x,y, z

      afun = -1.0

      return
      end
      function bfun (x,y,z)

c*********************************************************************72
c
      implicit none

      double precision bfun, x,y, z

      bfun = -1.0

      return
      end
      function cfun (x,y,z)

c*********************************************************************72
c
      implicit none

      double precision cfun, x,y, z

      cfun = -1.0

      return
      end
      function dfun (x,y,z)

c*********************************************************************72
c
      implicit none

      double precision dfun, x,y, z, gammax, gammay, alpha

      common /func/ gammax, gammay, alpha

      dfun = gammax*exp(x*y)

      return
      end
      function efun (x,y,z)

c*********************************************************************72
c
      implicit none

      double precision efun, x,y, z, gammax, gammay, alpha
      common /func/ gammax, gammay, alpha

      efun = gammay*exp(-x*y)

      return
      end
      function ffun (x,y,z)

c*********************************************************************72
c
      implicit none

      double precision ffun, x,y, z

      ffun = 0.0

      return
      end
      function gfun (x,y,z)

c*********************************************************************72
c
      implicit none

      double precision gfun, x,y, z, gammax, gammay, alpha

      common /func/ gammax, gammay, alpha
      gfun = alpha

      return
      end
      subroutine afunbl (nfree,x,y,z,coeff)

c*********************************************************************72
c
      implicit none

      double precision coeff(*)
      integer nfree
      double precision x
      double precision y
      double precision z

      return
      end
      subroutine bfunbl (nfree,x,y,z,coeff)

c*********************************************************************72
c
      implicit none

      double precision coeff(*)
      integer nfree
      double precision x
      double precision y
      double precision z

      return
      end
      subroutine cfunbl (nfree,x,y,z,coeff)

c*********************************************************************72
c
      implicit none

      double precision coeff(*)
      integer nfree
      double precision x
      double precision y
      double precision z

      return
      end
      subroutine dfunbl (nfree,x,y,z,coeff)

c*********************************************************************72
c
      implicit none

      double precision coeff(*)
      integer nfree
      double precision x
      double precision y
      double precision z

      return
      end
      subroutine efunbl (nfree,x,y,z,coeff)

c*********************************************************************72
c
      implicit none

      double precision coeff(*)
      integer nfree
      double precision x
      double precision y
      double precision z

      return
      end
      subroutine ffunbl (nfree,x,y,z,coeff)
     
c*********************************************************************72
c
      implicit none

      double precision coeff(*)
      integer nfree
      double precision x
      double precision y
      double precision z

      return
      end
      subroutine gfunbl (nfree,x,y,z,coeff)

c*********************************************************************72
c
      implicit none

      double precision coeff(*)
      integer nfree
      double precision x
      double precision y
      double precision z

      return
      end
      function random ( )
          
c*********************************************************************72
c
c  This routine was extracted from ELEFUNT.
c
      integer iy
      double precision random

      save iy

      data iy / 100001 /

      iy = iy * 125
      iy = iy - (iy/2796203) * 2796203
      random = dble ( iy ) / 2796203.0

      return
      end
