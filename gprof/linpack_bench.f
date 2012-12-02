      program main

      double precision a(1001,1000),b(1000),x(1000)
      double precision time(6),cray,ops,total,norma,normx
      double precision resid,residn,eps,epslon
      integer ipvt(1000)
      lda = 1001
c
c     this program was updated on 10/12/92 to correct a
c     problem with the random number generator. The previous
c     random number generator had a short period and produced
c     singular matrices occasionally.
c
      call timestamp ( )
      write ( *, '(a)' ) ' '

      n = 1000
      cray = .056
      write(6,1)
    1 format(' Please send the results of this run to:'//
     $       ' Jack J. Dongarra'/
     $       ' Computer Science Department'/
     $       ' University of Tennessee'/
     $       ' Knoxville, Tennessee 37996-1300'//
     $       ' Fax: 615-974-8296'//
     $       ' Internet: dongarra@cs.utk.edu'/)
      ops = (2.0d0*dfloat(n)**3)/3.0d0 + 2.0d0*dfloat(n)**2
c
         call matgen(a,lda,n,b,norma)
c
c******************************************************************
c******************************************************************
c        you should replace the call to dgefa and dgesl
c        by calls to your linear equation solver.
c******************************************************************
c******************************************************************
c
         t1 = second()
         call dgefa(a,lda,n,ipvt,info)
         time(1) = second() - t1
         t1 = second()
         call dgesl(a,lda,n,ipvt,b,0)
         time(2) = second() - t1
         total = time(1) + time(2)
c******************************************************************
c******************************************************************
c
c     compute a residual to verify results.
c
         do 10 i = 1,n
            x(i) = b(i)
   10    continue
         call matgen(a,lda,n,b,norma)
         do 20 i = 1,n
            b(i) = -b(i)
   20    continue
         call dmxpy(n,b,n,lda,x,a)
         resid = 0.0
         normx = 0.0
         do 30 i = 1,n
            resid = dmax1( resid, dabs(b(i)) )
            normx = dmax1( normx, dabs(x(i)) )
   30    continue
         eps = epslon(1.0d0)
         residn = resid/( n*norma*normx*eps )
         write(6,40)
   40    format('     norm. resid      resid           machep',
     $          '         x(1)          x(n)')
         write(6,50) residn,resid,eps,x(1),x(n)
   50    format(1p5e16.8)
c
         write(6,60) n
   60    format(//'    times are reported for matrices of order ',i5)
         write(6,70)
   70    format(6x,'factor',5x,'solve',6x,'total',5x,'mflops',7x,'unit',
     $         6x,'ratio')
c
         time(3) = total
         time(4) = ops/(1.0d6*total)
         time(5) = 2.0d0/time(4)
         time(6) = total/cray
         write(6,80) lda
   80    format(' times for array with leading dimension of',i4)
         write(6,110) (time(i),i=1,6)
  110    format(6(1pe11.3))
         write(6,*)' end of tests -- this version dated 10/12/92'
c

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine matgen(a,lda,n,b,norma)
      integer lda,n,init(4),i,j
      double precision a(lda,1),b(1),norma,random_value
c
      init(1) = 1
      init(2) = 2
      init(3) = 3
      init(4) = 1325
      norma = 0.0
      do 30 j = 1,n
         do 20 i = 1,n
            a(i,j) = random_value(init) - .5
            norma = dmax1(dabs(a(i,j)), norma)
   20    continue
   30 continue
      do 35 i = 1,n
          b(i) = 0.0
   35 continue
      do 50 j = 1,n
         do 40 i = 1,n
            b(i) = b(i) + a(i,j)
   40    continue
   50 continue
      return
      end
      subroutine dgefa(a,lda,n,ipvt,info)
      integer lda,n,ipvt(1),info
      double precision a(lda,1)
c
c     dgefa factors a double precision matrix by gaussian elimination.
c
c     dgefa is usually called by dgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
c
c     on entry
c
c        a       double precision(lda, n)
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
c                     indicate that dgesl or dgedi will divide by zero
c                     if called.  use  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c
c     internal variables
c
      double precision t
      integer idamax,j,k,kp1,l,nm1
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
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0d0) go to 40
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
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end
      subroutine dgesl(a,lda,n,ipvt,b,job)
      integer lda,n,ipvt(1),job
      double precision a(lda,1),b(1)
c
c     dgesl solves the double precision system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        b       double precision(n)
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
c        called correctly and if dgeco has set rcond .gt. 0.0
c        or dgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call dgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c
c     internal variables
c
      double precision ddot,t
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
            call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call daxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            t = ddot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
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
      subroutine daxpy(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end
      double precision function ddot(n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      ddot = 0.0d0
      dtemp = 0.0d0
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
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end
      subroutine  dscal(n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision da,dx(1)
      integer i,incx,m,mp1,n,nincx
c
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end
      integer function idamax(n,dx,incx)
c
c     finds the index of element having max. dabsolute value.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dmax
      integer i,incx,ix,n
c
      idamax = 0
      if( n .lt. 1 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end
      double precision function epslon (x)
      double precision x
c
c     estimate unit roundoff in quantities of size x.
c
      double precision a,b,c,eps
c
c     this program should function properly on all systems
c     satisfying the following two assumptions,
c        1.  the base used in representing dfloating point
c            numbers is not a power of three.
c        2.  the quantity  a  in statement 10 is represented to 
c            the accuracy used in dfloating point variables
c            that are stored in memory.
c     the statement number 10 and the go to 10 are intended to
c     force optimizing compilers to generate code satisfying 
c     assumption 2.
c     under these assumptions, it should be true that,
c            a  is not exactly equal to four-thirds,
c            b  has a zero for its last bit or digit,
c            c  is not exactly equal to one,
c            eps  measures the separation of 1.0 from
c                 the next larger dfloating point number.
c     the developers of eispack would appreciate being informed
c     about any systems where these assumptions do not hold.
c
c     *****************************************************************
c     this routine is one of the auxiliary routines used by eispack iii
c     to avoid machine dependencies.
c     *****************************************************************
c
c     this version dated 4/6/83.
c
      a = 4.0d0/3.0d0
   10 b = a - 1.0d0
      c = b + b + b
      eps = dabs(c-1.0d0)
      if (eps .eq. 0.0d0) go to 10
      epslon = eps*dabs(x)
      return
      end
      subroutine mm (a, lda, n1, n3, b, ldb, n2, c, ldc)
      double precision a(lda,*), b(ldb,*), c(ldc,*)
c
c   purpose:
c     multiply matrix b times matrix c and store the result in matrix a.
c
c   parameters:
c
c     a double precision(lda,n3), matrix of n1 rows and n3 columns
c
c     lda integer, leading dimension of array a
c
c     n1 integer, number of rows in matrices a and b
c
c     n3 integer, number of columns in matrices a and c
c
c     b double precision(ldb,n2), matrix of n1 rows and n2 columns
c
c     ldb integer, leading dimension of array b
c
c     n2 integer, number of columns in matrix b, and number of rows in
c         matrix c
c
c     c double precision(ldc,n3), matrix of n2 rows and n3 columns
c
c     ldc integer, leading dimension of array c
c
c ----------------------------------------------------------------------
c
      do 20 j = 1, n3
         do 10 i = 1, n1
            a(i,j) = 0.0
   10    continue
         call dmxpy (n2,a(1,j),n1,ldb,c(1,j),b)
   20 continue
c
      return
      end
      subroutine dmxpy (n1, y, n2, ldm, x, m)
      double precision y(*), x(*), m(ldm,*)
c
c   purpose:
c     multiply matrix m times vector x and add the result to vector y.
c
c   parameters:
c
c     n1 integer, number of elements in vector y, and number of rows in
c         matrix m
c
c     y double precision(n1), vector of length n1 to which is added 
c         the product m*x
c
c     n2 integer, number of elements in vector x, and number of columns
c         in matrix m
c
c     ldm integer, leading dimension of array m
c
c     x double precision(n2), vector of length n2
c
c     m double precision(ldm,n2), matrix of n1 rows and n2 columns
c
c ----------------------------------------------------------------------
c
c   cleanup odd vector
c
      j = mod(n2,2)
      if (j .ge. 1) then
         do 10 i = 1, n1
            y(i) = (y(i)) + x(j)*m(i,j)
   10    continue
      endif
c
c   cleanup odd group of two vectors
c
      j = mod(n2,4)
      if (j .ge. 2) then
         do 20 i = 1, n1
            y(i) = ( (y(i))
     $             + x(j-1)*m(i,j-1)) + x(j)*m(i,j)
   20    continue
      endif
c
c   cleanup odd group of four vectors
c
      j = mod(n2,8)
      if (j .ge. 4) then
         do 30 i = 1, n1
            y(i) = ((( (y(i))
     $             + x(j-3)*m(i,j-3)) + x(j-2)*m(i,j-2))
     $             + x(j-1)*m(i,j-1)) + x(j)  *m(i,j)
   30    continue
      endif
c
c   cleanup odd group of eight vectors
c
      j = mod(n2,16)
      if (j .ge. 8) then
         do 40 i = 1, n1
            y(i) = ((((((( (y(i))
     $             + x(j-7)*m(i,j-7)) + x(j-6)*m(i,j-6))
     $             + x(j-5)*m(i,j-5)) + x(j-4)*m(i,j-4))
     $             + x(j-3)*m(i,j-3)) + x(j-2)*m(i,j-2))
     $             + x(j-1)*m(i,j-1)) + x(j)  *m(i,j)
   40    continue
      endif
c
c   main loop - groups of sixteen vectors
c
      jmin = j+16
      do 60 j = jmin, n2, 16
         do 50 i = 1, n1
            y(i) = ((((((((((((((( (y(i))
     $             + x(j-15)*m(i,j-15)) + x(j-14)*m(i,j-14))
     $             + x(j-13)*m(i,j-13)) + x(j-12)*m(i,j-12))
     $             + x(j-11)*m(i,j-11)) + x(j-10)*m(i,j-10))
     $             + x(j- 9)*m(i,j- 9)) + x(j- 8)*m(i,j- 8))
     $             + x(j- 7)*m(i,j- 7)) + x(j- 6)*m(i,j- 6))
     $             + x(j- 5)*m(i,j- 5)) + x(j- 4)*m(i,j- 4))
     $             + x(j- 3)*m(i,j- 3)) + x(j- 2)*m(i,j- 2))
     $             + x(j- 1)*m(i,j- 1)) + x(j)   *m(i,j)
   50    continue
   60 continue
      return
      end
      DOUBLE PRECISION FUNCTION RANDOM_VALUE( ISEED )
c
c     modified from the LAPACK auxiliary routine 10/12/92 JD
c  -- LAPACK auxiliary routine (version 1.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c     February 29, 1992
c
c     .. Array Arguments ..
      INTEGER            ISEED( 4 )
c     ..
c
c  Purpose
c  =======
c
c  DLARAN returns a random real number from a uniform (0,1)
c  distribution.
c
c  Arguments
c  =========
c
c  ISEED   (input/output) INTEGER array, dimension (4)
c          On entry, the seed of the random number generator; the array
c          elements must be between 0 and 4095, and ISEED(4) must be
c          odd.
c          On exit, the seed is updated.
c
c  Further Details
c  ===============
c
c  This routine uses a multiplicative congruential method with modulus
c  2**48 and multiplier 33952834046453 (see G.S.Fishman,
c  'Multiplicative congruential random number generators with modulus
c  2**b: an exhaustive analysis for b = 32 and a partial analysis for
c  b = 48', Math. Comp. 189, pp 331-344, 1990).
c
c  48-bit integers are stored in 4 integer array elements with 12 bits
c  per element. Hence the routine is portable across machines with
c  integers of 32 bits or more.
c
c     .. Parameters ..
      INTEGER            M1, M2, M3, M4
      PARAMETER          ( M1 = 494, M2 = 322, M3 = 2508, M4 = 2549 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
      INTEGER            IPW2
      DOUBLE PRECISION   R
      PARAMETER          ( IPW2 = 4096, R = ONE / IPW2 )
c     ..
c     .. Local Scalars ..
      INTEGER            IT1, IT2, IT3, IT4
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MOD
c     ..
c     .. Executable Statements ..
c
c     multiply the seed by the multiplier modulo 2**48
c
      IT4 = ISEED( 4 )*M4
      IT3 = IT4 / IPW2
      IT4 = IT4 - IPW2*IT3
      IT3 = IT3 + ISEED( 3 )*M4 + ISEED( 4 )*M3
      IT2 = IT3 / IPW2
      IT3 = IT3 - IPW2*IT2
      IT2 = IT2 + ISEED( 2 )*M4 + ISEED( 3 )*M3 + ISEED( 4 )*M2
      IT1 = IT2 / IPW2
      IT2 = IT2 - IPW2*IT1
      IT1 = IT1 + ISEED( 1 )*M4 + ISEED( 2 )*M3 + ISEED( 3 )*M2 +
     $      ISEED( 4 )*M1
      IT1 = MOD( IT1, IPW2 )
c
c     return updated seed
c
      ISEED( 1 ) = IT1
      ISEED( 2 ) = IT2
      ISEED( 3 ) = IT3
      ISEED( 4 ) = IT4
c
c     convert 48-bit integer to a real number in the interval (0,1)
c
      RANDOM_VALUE = R*( DBLE( IT1 )+R*( DBLE( IT2 )+R*( DBLE( IT3 )
     &  +R*( DBLE( IT4 ) ) ) ) )
      RETURN
c
c     End of RAN
c
      END
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
