      subroutine cchdc(a,lda,p,work,jpvt,job,info)

c*********************************************************************72

      integer lda,p,jpvt(1),job,info
      complex a(lda,1),work(1)
c
c     cchdc computes the cholesky decomposition of a positive definite
c     matrix.  a pivoting option allows the user to estimate the
c     condition of a positive definite matrix or determine the rank
c     of a positive semidefinite matrix.
c
c     on entry
c
c         a      complex(lda,p).
c                a contains the matrix whose decomposition is to
c                be computed.  onlt the upper half of a need be stored.
c                the lower part of the array a is not referenced.
c
c         lda    integer.
c                lda is the leading dimension of the array a.
c
c         p      integer.
c                p is the order of the matrix.
c
c         work   complex.
c                work is a work array.
c
c         jpvt   integer(p).
c                jpvt contains integers that control the selection
c                of the pivot elements, if pivoting has been requested.
c                each diagonal element a(k,k)
c                is placed in one of three classes according to the
c                value of jpvt(k).
c
c                   if jpvt(k) .gt. 0, then x(k) is an initial
c                                      element.
c
c                   if jpvt(k) .eq. 0, then x(k) is a free element.
c
c                   if jpvt(k) .lt. 0, then x(k) is a final element.
c
c                before the decomposition is computed, initial elements
c                are moved by symmetric row and column interchanges to
c                the beginning of the array a and final
c                elements to the end.  both initial and final elements
c                are frozen in place during the computation and only
c                free elements are moved.  at the k-th stage of the
c                reduction, if a(k,k) is occupied by a free element
c                it is interchanged with the largest free element
c                a(l,l) with l .ge. k.  jpvt is not referenced if
c                job .eq. 0.
c
c        job     integer.
c                job is an integer that initiates column pivoting.
c                if job .eq. 0, no pivoting is done.
c                if job .ne. 0, pivoting is done.
c
c     on return
c
c         a      a contains in its upper half the cholesky factor
c                of the matrix a as it has been permuted by pivoting.
c
c         jpvt   jpvt(j) contains the index of the diagonal element
c                of a that was moved into the j-th position,
c                provided pivoting was requested.
c
c         info   contains the index of the last positive diagonal
c                element of the cholesky factor.
c
c     for positive definite matrices info = p is the normal return.
c     for pivoting with positive semidefinite matrices info will
c     in general be less than p.  however, info may be greater than
c     the rank of a, since rounding error can cause an otherwise zero
c     element to be positive. indefinite systems will always cause
c     info to be less than p.
c
c     linpack. this version dated 03/19/79 .
c     j.j. dongarra and g.w. stewart, argonne national laboratory and
c     university of maryland.
c
c
c     blas caxpy,cswap
c     fortran sqrt,real,conjg
c
c     internal variables
c
      integer pu,pl,plp1,i,j,jp,jt,k,kb,km1,kp1,l,maxl
      complex temp
      real maxdia
      logical swapk,negk
c
      pl = 1
      pu = 0
      info = p
      if (job .eq. 0) go to 160
c
c        pivoting has been requested. rearrange the
c        the elements according to jpvt.
c
         do 70 k = 1, p
            swapk = jpvt(k) .gt. 0
            negk = jpvt(k) .lt. 0
            jpvt(k) = k
            if (negk) jpvt(k) = -jpvt(k)
            if (.not.swapk) go to 60
               if (k .eq. pl) go to 50
                  call cswap(pl-1,a(1,k),1,a(1,pl),1)
                  temp = a(k,k)
                  a(k,k) = a(pl,pl)
                  a(pl,pl) = temp
                  a(pl,k) = conjg(a(pl,k))
                  plp1 = pl + 1
                  if (p .lt. plp1) go to 40
                  do 30 j = plp1, p
                     if (j .ge. k) go to 10
                        temp = conjg(a(pl,j))
                        a(pl,j) = conjg(a(j,k))
                        a(j,k) = temp
                     go to 20
   10                continue
                     if (j .eq. k) go to 20
                        temp = a(k,j)
                        a(k,j) = a(pl,j)
                        a(pl,j) = temp
   20                continue
   30             continue
   40             continue
                  jpvt(k) = jpvt(pl)
                  jpvt(pl) = k
   50          continue
               pl = pl + 1
   60       continue
   70    continue
         pu = p
         if (p .lt. pl) go to 150
         do 140 kb = pl, p
            k = p - kb + pl
            if (jpvt(k) .ge. 0) go to 130
               jpvt(k) = -jpvt(k)
               if (pu .eq. k) go to 120
                  call cswap(k-1,a(1,k),1,a(1,pu),1)
                  temp = a(k,k)
                  a(k,k) = a(pu,pu)
                  a(pu,pu) = temp
                  a(k,pu) = conjg(a(k,pu))
                  kp1 = k + 1
                  if (p .lt. kp1) go to 110
                  do 100 j = kp1, p
                     if (j .ge. pu) go to 80
                        temp = conjg(a(k,j))
                        a(k,j) = conjg(a(j,pu))
                        a(j,pu) = temp
                     go to 90
   80                continue
                     if (j .eq. pu) go to 90
                        temp = a(k,j)
                        a(k,j) = a(pu,j)
                        a(pu,j) = temp
   90                continue
  100             continue
  110             continue
                  jt = jpvt(k)
                  jpvt(k) = jpvt(pu)
                  jpvt(pu) = jt
  120          continue
               pu = pu - 1
  130       continue
  140    continue
  150    continue
  160 continue
      do 270 k = 1, p
c
c        reduction loop.
c
         maxdia = real(a(k,k))
         kp1 = k + 1
         maxl = k
c
c        determine the pivot element.
c
         if (k .lt. pl .or. k .ge. pu) go to 190
            do 180 l = kp1, pu
               if (real(a(l,l)) .le. maxdia) go to 170
                  maxdia = real(a(l,l))
                  maxl = l
  170          continue
  180       continue
  190    continue
c
c        quit if the pivot element is not positive.
c
         if (maxdia .gt. 0.0e0) go to 200
            info = k - 1
c     ......exit
            go to 280
  200    continue
         if (k .eq. maxl) go to 210
c
c           start the pivoting and update jpvt.
c
            km1 = k - 1
            call cswap(km1,a(1,k),1,a(1,maxl),1)
            a(maxl,maxl) = a(k,k)
            a(k,k) = cmplx(maxdia,0.0e0)
            jp = jpvt(maxl)
            jpvt(maxl) = jpvt(k)
            jpvt(k) = jp
            a(k,maxl) = conjg(a(k,maxl))
  210    continue
c
c        reduction step. pivoting is contained across the rows.
c
         work(k) = cmplx(sqrt(real(a(k,k))),0.0e0)
         a(k,k) = work(k)
         if (p .lt. kp1) go to 260
         do 250 j = kp1, p
            if (k .eq. maxl) go to 240
               if (j .ge. maxl) go to 220
                  temp = conjg(a(k,j))
                  a(k,j) = conjg(a(j,maxl))
                  a(j,maxl) = temp
               go to 230
  220          continue
               if (j .eq. maxl) go to 230
                  temp = a(k,j)
                  a(k,j) = a(maxl,j)
                  a(maxl,j) = temp
  230          continue
  240       continue
            a(k,j) = a(k,j)/work(k)
            work(j) = conjg(a(k,j))
            temp = -a(k,j)
            call caxpy(j-k,temp,work(kp1),1,a(kp1,j),1)
  250    continue
  260    continue
  270 continue
  280 continue
      return
      end
      subroutine cchdd(r,ldr,p,x,z,ldz,nz,y,rho,c,s,info)

c*********************************************************************72

      integer ldr,p,ldz,nz,info
      complex r(ldr,1),x(1),z(ldz,1),y(1),s(1)
      real rho(1),c(1)
c
c     cchdd downdates an augmented cholesky decomposition or the
c     triangular factor of an augmented qr decomposition.
c     specifically, given an upper triangular matrix r of order p,  a
c     row vector x, a column vector z, and a scalar y, cchdd
c     determineds a unitary matrix u and a scalar zeta such that
c
c                        (r   z )     (rr  zz)
c                    u * (      )  =  (      ) ,
c                        (0 zeta)     ( x   y)
c
c     where rr is upper triangular.  if r and z have been obtained
c     from the factorization of a least squares problem, then
c     rr and zz are the factors corresponding to the problem
c     with the observation (x,y) removed.  in this case, if rho
c     is the norm of the residual vector, then the norm of
c     the residual vector of the downdated problem is
c     sqrt(rho**2 - zeta**2). cchdd will simultaneously downdate
c     several triplets (z,y,rho) along with r.
c     for a less terse description of what cchdd does and how
c     it may be applied, see the linpack guide.
c
c     the matrix u is determined as the product u(1)*...*u(p)
c     where u(i) is a rotation in the (p+1,i)-plane of the
c     form
c
c                       ( c(i)  -conjg(s(i)) )
c                       (                    ) .
c                       ( s(i)       c(i)    )
c
c     the rotations are chosen so that c(i) is real.
c
c     the user is warned that a given downdating problem may
c     be impossible to accomplish or may produce
c     inaccurate results.  for example, this can happen
c     if x is near a vector whose removal will reduce the
c     rank of r.  beware.
c
c     on entry
c
c         r      complex(ldr,p), where ldr .ge. p.
c                r contains the upper triangular matrix
c                that is to be downdated.  the part of  r
c                below the diagonal is not referenced.
c
c         ldr    integer.
c                ldr is the leading dimension fo the array r.
c
c         p      integer.
c                p is the order of the matrix r.
c
c         x      complex(p).
c                x contains the row vector that is to
c                be removed from r.  x is not altered by cchdd.
c
c         z      complex(ldz,nz), where ldz .ge. p.
c                z is an array of nz p-vectors which
c                are to be downdated along with r.
c
c         ldz    integer.
c                ldz is the leading dimension of the array z.
c
c         nz     integer.
c                nz is the number of vectors to be downdated
c                nz may be zero, in which case z, y, and rho
c                are not referenced.
c
c         y      complex(nz).
c                y contains the scalars for the downdating
c                of the vectors z.  y is not altered by cchdd.
c
c         rho    real(nz).
c                rho contains the norms of the residual
c                vectors that are to be downdated.
c
c     on return
c
c         r
c         z      contain the downdated quantities.
c         rho
c
c         c      real(p).
c                c contains the cosines of the transforming
c                rotations.
c
c         s      complex(p).
c                s contains the sines of the transforming
c                rotations.
c
c         info   integer.
c                info is set as follows.
c
c                   info = 0  if the entire downdating
c                             was successful.
c
c                   info =-1  if r could not be downdated.
c                             in this case, all quantities
c                             are left unaltered.
c
c                   info = 1  if some rho could not be
c                             downdated.  the offending rhos are
c                             set to -1.
c
c     linpack. this version dated 08/14/78 .
c     g.w. stewart, university of maryland, argonne national lab.
c
c     cchdd uses the following functions and subprograms.
c
c     fortran aimag,cabs,conjg
c     blas cdotc, scnrm2
c
      integer i,ii,j
      real a,alpha,azeta,norm,scnrm2
      complex cdotc,t,zeta,b,xx
c
c     solve the system ctrans(r)*a = x, placing the result
c     in the array s.
c
      info = 0
      s(1) = conjg(x(1))/conjg(r(1,1))
      if (p .lt. 2) go to 20
      do 10 j = 2, p
         s(j) = conjg(x(j)) - cdotc(j-1,r(1,j),1,s,1)
         s(j) = s(j)/conjg(r(j,j))
   10 continue
   20 continue
      norm = scnrm2(p,s,1)
      if (norm .lt. 1.0e0) go to 30
         info = -1
      go to 120
   30 continue
         alpha = sqrt(1.0e0-norm**2)
c
c        determine the transformations.
c
         do 40 ii = 1, p
            i = p - ii + 1
            scale = alpha + cabs(s(i))
            a = alpha/scale
            b = s(i)/scale
            norm = sqrt(a**2+real(b)**2+aimag(b)**2)
            c(i) = a/norm
            s(i) = conjg(b)/norm
            alpha = scale*norm
   40    continue
c
c        apply the transformations to r.
c
         do 60 j = 1, p
            xx = (0.0e0,0.0e0)
            do 50 ii = 1, j
               i = j - ii + 1
               t = c(i)*xx + s(i)*r(i,j)
               r(i,j) = c(i)*r(i,j) - conjg(s(i))*xx
               xx = t
   50       continue
   60    continue
c
c        if required, downdate z and rho.
c
         if (nz .lt. 1) go to 110
         do 100 j = 1, nz
            zeta = y(j)
            do 70 i = 1, p
               z(i,j) = (z(i,j) - conjg(s(i))*zeta)/c(i)
               zeta = c(i)*zeta - s(i)*z(i,j)
   70       continue
            azeta = cabs(zeta)
            if (azeta .le. rho(j)) go to 80
               info = 1
               rho(j) = -1.0e0
            go to 90
   80       continue
               rho(j) = rho(j)*sqrt(1.0e0-(azeta/rho(j))**2)
   90       continue
  100    continue
  110    continue
  120 continue
      return
      end
      subroutine cchex(r,ldr,p,k,l,z,ldz,nz,c,s,job)

c*********************************************************************72

      integer ldr,p,k,l,ldz,nz,job
      complex r(ldr,1),z(ldz,1),s(1)
      real c(1)
c
c     cchex updates the cholesky factorization
c
c                   a = ctrans(r)*r
c
c     of a positive definite matrix a of order p under diagonal
c     permutations of the form
c
c                   trans(e)*a*e
c
c     where e is a permutation matrix.  specifically, given
c     an upper triangular matrix r and a permutation matrix
c     e (which is specified by k, l, and job), cchex determines
c     a unitary matrix u such that
c
c                           u*r*e = rr,
c
c     where rr is upper triangular.  at the users option, the
c     transformation u will be multiplied into the array z.
c     if a = ctrans(x)*x, so that r is the triangular part of the
c     qr factorization of x, then rr is the triangular part of the
c     qr factorization of x*e, i.e. x with its columns permuted.
c     for a less terse description of what cchex does and how
c     it may be applied, see the linpack guide.
c
c     the matrix q is determined as the product u(l-k)*...*u(1)
c     of plane rotations of the form
c
c                           (    c(i)       s(i) )
c                           (                    ) ,
c                           ( -conjg(s(i))  c(i) )
c
c     where c(i) is real, the rows these rotations operate on
c     are described below.
c
c     there are two types of permutations, which are determined
c     by the value of job.
c
c     1. right circular shift (job = 1).
c
c         the columns are rearranged in the following order.
c
c                1,...,k-1,l,k,k+1,...,l-1,l+1,...,p.
c
c         u is the product of l-k rotations u(i), where u(i)
c         acts in the (l-i,l-i+1)-plane.
c
c     2. left circular shift (job = 2).
c         the columns are rearranged in the following order
c
c                1,...,k-1,k+1,k+2,...,l,k,l+1,...,p.
c
c         u is the product of l-k rotations u(i), where u(i)
c         acts in the (k+i-1,k+i)-plane.
c
c     on entry
c
c         r      complex(ldr,p), where ldr.ge.p.
c                r contains the upper triangular factor
c                that is to be updated.  elements of r
c                below the diagonal are not referenced.
c
c         ldr    integer.
c                ldr is the leading dimension of the array r.
c
c         p      integer.
c                p is the order of the matrix r.
c
c         k      integer.
c                k is the first column to be permuted.
c
c         l      integer.
c                l is the last column to be permuted.
c                l must be strictly greater than k.
c
c         z      complex(ldz,nz), where ldz.ge.p.
c                z is an array of nz p-vectors into which the
c                transformation u is multiplied.  z is
c                not referenced if nz = 0.
c
c         ldz    integer.
c                ldz is the leading dimension of the array z.
c
c         nz     integer.
c                nz is the number of columns of the matrix z.
c
c         job    integer.
c                job determines the type of permutation.
c                       job = 1  right circular shift.
c                       job = 2  left circular shift.
c
c     on return
c
c         r      contains the updated factor.
c
c         z      contains the updated matrix z.
c
c         c      real(p).
c                c contains the cosines of the transforming rotations.
c
c         s      complex(p).
c                s contains the sines of the transforming rotations.
c
c     linpack. this version dated 08/14/78 .
c     g.w. stewart, university of maryland, argonne national lab.
c
c     cchex uses the following functions and subroutines.
c
c     extended blas crotg
c     fortran min0
c
      integer i,ii,il,iu,j,jj,km1,kp1,lmk,lm1
      complex rjp1j,t
c
c     initialize
c
      km1 = k - 1
      kp1 = k + 1
      lmk = l - k
      lm1 = l - 1
c
c     perform the appropriate task.
c
      go to (10,130), job
c
c     right circular shift.
c
   10 continue
c
c        reorder the columns.
c
         do 20 i = 1, l
            ii = l - i + 1
            s(i) = r(ii,l)
   20    continue
         do 40 jj = k, lm1
            j = lm1 - jj + k
            do 30 i = 1, j
               r(i,j+1) = r(i,j)
   30       continue
            r(j+1,j+1) = (0.0e0,0.0e0)
   40    continue
         if (k .eq. 1) go to 60
            do 50 i = 1, km1
               ii = l - i + 1
               r(i,k) = s(ii)
   50       continue
   60    continue
c
c        calculate the rotations.
c
         t = s(1)
         do 70 i = 1, lmk
            call crotg(s(i+1),t,c(i),s(i))
            t = s(i+1)
   70    continue
         r(k,k) = t
         do 90 j = kp1, p
            il = max0(1,l-j+1)
            do 80 ii = il, lmk
               i = l - ii
               t = c(ii)*r(i,j) + s(ii)*r(i+1,j)
               r(i+1,j) = c(ii)*r(i+1,j) - conjg(s(ii))*r(i,j)
               r(i,j) = t
   80       continue
   90    continue
c
c        if required, apply the transformations to z.
c
         if (nz .lt. 1) go to 120
         do 110 j = 1, nz
            do 100 ii = 1, lmk
               i = l - ii
               t = c(ii)*z(i,j) + s(ii)*z(i+1,j)
               z(i+1,j) = c(ii)*z(i+1,j) - conjg(s(ii))*z(i,j)
               z(i,j) = t
  100       continue
  110    continue
  120    continue
      go to 260
c
c     left circular shift
c
  130 continue
c
c        reorder the columns
c
         do 140 i = 1, k
            ii = lmk + i
            s(ii) = r(i,k)
  140    continue
         do 160 j = k, lm1
            do 150 i = 1, j
               r(i,j) = r(i,j+1)
  150       continue
            jj = j - km1
            s(jj) = r(j+1,j+1)
  160    continue
         do 170 i = 1, k
            ii = lmk + i
            r(i,l) = s(ii)
  170    continue
         do 180 i = kp1, l
            r(i,l) = (0.0e0,0.0e0)
  180    continue
c
c        reduction loop.
c
         do 220 j = k, p
            if (j .eq. k) go to 200
c
c              apply the rotations.
c
               iu = min0(j-1,l-1)
               do 190 i = k, iu
                  ii = i - k + 1
                  t = c(ii)*r(i,j) + s(ii)*r(i+1,j)
                  r(i+1,j) = c(ii)*r(i+1,j) - conjg(s(ii))*r(i,j)
                  r(i,j) = t
  190          continue
  200       continue
            if (j .ge. l) go to 210
               jj = j - k + 1
               t = s(jj)
               call crotg(r(j,j),t,c(jj),s(jj))
  210       continue
  220    continue
c
c        apply the rotations to z.
c
         if (nz .lt. 1) go to 250
         do 240 j = 1, nz
            do 230 i = k, lm1
               ii = i - km1
               t = c(ii)*z(i,j) + s(ii)*z(i+1,j)
               z(i+1,j) = c(ii)*z(i+1,j) - conjg(s(ii))*z(i,j)
               z(i,j) = t
  230       continue
  240    continue
  250    continue
  260 continue
      return
      end
      subroutine cchud(r,ldr,p,x,z,ldz,nz,y,rho,c,s)

c*********************************************************************72

      integer ldr,p,ldz,nz
      real rho(1),c(1)
      complex r(ldr,1),x(1),z(ldz,1),y(1),s(1)
c
c     cchud updates an augmented cholesky decomposition of the
c     triangular part of an augmented qr decomposition.  specifically,
c     given an upper triangular matrix r of order p, a row vector
c     x, a column vector z, and a scalar y, cchud determines a
c     untiary matrix u and a scalar zeta such that
c
c
c                              (r  z)     (rr   zz )
c                         u  * (    )  =  (        ) ,
c                              (x  y)     ( 0  zeta)
c
c     where rr is upper triangular.  if r and z have been
c     obtained from the factorization of a least squares
c     problem, then rr and zz are the factors corresponding to
c     the problem with the observation (x,y) appended.  in this
c     case, if rho is the norm of the residual vector, then the
c     norm of the residual vector of the updated problem is
c     sqrt(rho**2 + zeta**2).  cchud will simultaneously update
c     several triplets (z,y,rho).
c     for a less terse description of what cchud does and how
c     it may be applied see the linpack guide.
c
c     the matrix u is determined as the product u(p)*...*u(1),
c     where u(i) is a rotation in the (i,p+1) plane of the
c     form
c
c                       (     c(i)      s(i) )
c                       (                    ) .
c                       ( -conjg(s(i))  c(i) )
c
c     the rotations are chosen so that c(i) is real.
c
c     on entry
c
c         r      complex(ldr,p), where ldr .ge. p.
c                r contains the upper triangular matrix
c                that is to be updated.  the part of r
c                below the diagonal is not referenced.
c
c         ldr    integer.
c                ldr is the leading dimension of the array r.
c
c         p      integer.
c                p is the order of the matrix r.
c
c         x      complex(p).
c                x contains the row to be added to r.  x is
c                not altered by cchud.
c
c         z      complex(ldz,nz), where ldz .ge. p.
c                z is an array containing nz p-vectors to
c                be updated with r.
c
c         ldz    integer.
c                ldz is the leading dimension of the array z.
c
c         nz     integer.
c                nz is the number of vectors to be updated
c                nz may be zero, in which case z, y, and rho
c                are not referenced.
c
c         y      complex(nz).
c                y contains the scalars for updating the vectors
c                z.  y is not altered by cchud.
c
c         rho    real(nz).
c                rho contains the norms of the residual
c                vectors that are to be updated.  if rho(j)
c                is negative, it is left unaltered.
c
c     on return
c
c         rc
c         rho    contain the updated quantities.
c         z
c
c         c      real(p).
c                c contains the cosines of the transforming
c                rotations.
c
c         s      complex(p).
c                s contains the sines of the transforming
c                rotations.
c
c     linpack. this version dated 08/14/78 .
c     g.w. stewart, university of maryland, argonne national lab.
c
c     cchud uses the following functions and subroutines.
c
c     extended blas crotg
c     fortran conjg,sqrt
c
      integer i,j,jm1
      real azeta,scale
      complex t,xj,zeta
c
c     update r.
c
      do 30 j = 1, p
         xj = x(j)
c
c        apply the previous rotations.
c
         jm1 = j - 1
         if (jm1 .lt. 1) go to 20
         do 10 i = 1, jm1
            t = c(i)*r(i,j) + s(i)*xj
            xj = c(i)*xj - conjg(s(i))*r(i,j)
            r(i,j) = t
   10    continue
   20    continue
c
c        compute the next rotation.
c
         call crotg(r(j,j),xj,c(j),s(j))
   30 continue
c
c     if required, update z and rho.
c
      if (nz .lt. 1) go to 70
      do 60 j = 1, nz
         zeta = y(j)
         do 40 i = 1, p
            t = c(i)*z(i,j) + s(i)*zeta
            zeta = c(i)*zeta - conjg(s(i))*z(i,j)
            z(i,j) = t
   40    continue
         azeta = cabs(zeta)
         if (azeta .eq. 0.0e0 .or. rho(j) .lt. 0.0e0) go to 50
            scale = azeta + rho(j)
            rho(j) = scale*sqrt((azeta/scale)**2+(rho(j)/scale)**2)
   50    continue
   60 continue
   70 continue
      return
      end
      subroutine cgbco(abd,lda,n,ml,mu,ipvt,rcond,z)

c*********************************************************************72

      integer lda,n,ml,mu,ipvt(1)
      complex abd(lda,1),z(1)
      real rcond
c
c     cgbco factors a complex band matrix by gaussian
c     elimination and estimates the condition of the matrix.
c
c     if  rcond  is not needed, cgbfa is slightly faster.
c     to solve  a*x = b , follow cgbco by cgbsl.
c     to compute  inverse(a)*c , follow cgbco by cgbsl.
c     to compute  determinant(a) , follow cgbco by cgbdi.
c
c     on entry
c
c        abd     complex(lda, n)
c                contains the matrix in band storage.  the columns
c                of the matrix are stored in the columns of  abd  and
c                the diagonals of the matrix are stored in rows
c                ml+1 through 2*ml+mu+1 of  abd .
c                see the comments below for details.
c
c        lda     integer
c                the leading dimension of the array  abd .
c                lda must be .ge. 2*ml + mu + 1 .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c                0 .le. ml .lt. n .
c
c        mu      integer
c                number of diagonals above the main diagonal.
c                0 .le. mu .lt. n .
c                more efficient if  ml .le. mu .
c
c     on return
c
c        abd     an upper triangular matrix in band storage and
c                the multipliers which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        rcond   real
c                an estimate of the reciprocal condition of  a .
c                for the system  a*x = b , relative perturbations
c                in  a  and  b  of size  epsilon  may cause
c                relative perturbations in  x  of size  epsilon/rcond .
c                if  rcond  is so small that the logical expression
c                           1.0 + rcond .eq. 1.0
c                is true, then  a  may be singular to working
c                precision.  in particular,  rcond  is zero  if
c                exact singularity is detected or the estimate
c                underflows.
c
c        z       complex(n)
c                a work vector whose contents are usually unimportant.
c                if  a  is close to a singular matrix, then  z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .
c
c     band storage
c
c           if  a  is a band matrix, the following program segment
c           will set up the input.
c
c                   ml = (band width below the diagonal)
c                   mu = (band width above the diagonal)
c                   m = ml + mu + 1
c                   do 20 j = 1, n
c                      i1 = max0(1, j-mu)
c                      i2 = min0(n, j+ml)
c                      do 10 i = i1, i2
c                         k = i - j + m
c                         abd(k,j) = a(i,j)
c                10    continue
c                20 continue
c
c           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
c           in addition, the first  ml  rows in  abd  are used for
c           elements generated during the triangularization.
c           the total number of rows needed in  abd  is  2*ml+mu+1 .
c           the  ml+mu by ml+mu  upper left triangle and the
c           ml by ml  lower right triangle are not referenced.
c
c     example..  if the original matrix is
c
c           11 12 13  0  0  0
c           21 22 23 24  0  0
c            0 32 33 34 35  0
c            0  0 43 44 45 46
c            0  0  0 54 55 56
c            0  0  0  0 65 66
c
c      then  n = 6, ml = 1, mu = 2, lda .ge. 5  and abd should contain
c
c            *  *  *  +  +  +  , * = not used
c            *  * 13 24 35 46  , + = used for pivoting
c            * 12 23 34 45 56
c           11 22 33 44 55 66
c           21 32 43 54 65  *
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     linpack cgbfa
c     blas caxpy,cdotc,csscal,scasum
c     fortran abs,aimag,amax1,cmplx,conjg,max0,min0,real
c
c     internal variables
c
      complex cdotc,ek,t,wk,wkm
      real anorm,s,scasum,sm,ynorm
      integer is,info,j,ju,k,kb,kp1,l,la,lm,lz,m,mm
c
      complex zdum,zdum1,zdum2,csign1
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
      csign1(zdum1,zdum2) = cabs1(zdum1)*(zdum2/cabs1(zdum2))
c
c     compute 1-norm of a
c
      anorm = 0.0e0
      l = ml + 1
      is = l + mu
      do 10 j = 1, n
         anorm = amax1(anorm,scasum(l,abd(is,j),1))
         if (is .gt. ml + 1) is = is - 1
         if (j .le. mu) l = l + 1
         if (j .ge. n - ml) l = l - 1
   10 continue
c
c     factor
c
      call cgbfa(abd,lda,n,ml,mu,ipvt,info)
c
c     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
c     estimate = norm(z)/norm(y) where  a*z = y  and  ctrans(a)*y = e .
c     ctrans(a)  is the conjugate transpose of a .
c     the components of  e  are chosen to cause maximum local
c     growth in the elements of w  where  ctrans(u)*w = e .
c     the vectors are frequently rescaled to avoid overflow.
c
c     solve ctrans(u)*w = e
c
      ek = (1.0e0,0.0e0)
      do 20 j = 1, n
         z(j) = (0.0e0,0.0e0)
   20 continue
      m = ml + mu + 1
      ju = 0
      do 100 k = 1, n
         if (cabs1(z(k)) .ne. 0.0e0) ek = csign1(ek,-z(k))
         if (cabs1(ek-z(k)) .le. cabs1(abd(m,k))) go to 30
            s = cabs1(abd(m,k))/cabs1(ek-z(k))
            call csscal(n,s,z,1)
            ek = cmplx(s,0.0e0)*ek
   30    continue
         wk = ek - z(k)
         wkm = -ek - z(k)
         s = cabs1(wk)
         sm = cabs1(wkm)
         if (cabs1(abd(m,k)) .eq. 0.0e0) go to 40
            wk = wk/conjg(abd(m,k))
            wkm = wkm/conjg(abd(m,k))
         go to 50
   40    continue
            wk = (1.0e0,0.0e0)
            wkm = (1.0e0,0.0e0)
   50    continue
         kp1 = k + 1
         ju = min0(max0(ju,mu+ipvt(k)),n)
         mm = m
         if (kp1 .gt. ju) go to 90
            do 60 j = kp1, ju
               mm = mm - 1
               sm = sm + cabs1(z(j)+wkm*conjg(abd(mm,j)))
               z(j) = z(j) + wk*conjg(abd(mm,j))
               s = s + cabs1(z(j))
   60       continue
            if (s .ge. sm) go to 80
               t = wkm - wk
               wk = wkm
               mm = m
               do 70 j = kp1, ju
                  mm = mm - 1
                  z(j) = z(j) + t*conjg(abd(mm,j))
   70          continue
   80       continue
   90    continue
         z(k) = wk
  100 continue
      s = 1.0e0/scasum(n,z,1)
      call csscal(n,s,z,1)
c
c     solve ctrans(l)*y = w
c
      do 120 kb = 1, n
         k = n + 1 - kb
         lm = min0(ml,n-k)
         if (k .lt. n) z(k) = z(k) + cdotc(lm,abd(m+1,k),1,z(k+1),1)
         if (cabs1(z(k)) .le. 1.0e0) go to 110
            s = 1.0e0/cabs1(z(k))
            call csscal(n,s,z,1)
  110    continue
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
  120 continue
      s = 1.0e0/scasum(n,z,1)
      call csscal(n,s,z,1)
c
      ynorm = 1.0e0
c
c     solve l*v = y
c
      do 140 k = 1, n
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
         lm = min0(ml,n-k)
         if (k .lt. n) call caxpy(lm,t,abd(m+1,k),1,z(k+1),1)
         if (cabs1(z(k)) .le. 1.0e0) go to 130
            s = 1.0e0/cabs1(z(k))
            call csscal(n,s,z,1)
            ynorm = s*ynorm
  130    continue
  140 continue
      s = 1.0e0/scasum(n,z,1)
      call csscal(n,s,z,1)
      ynorm = s*ynorm
c
c     solve  u*z = w
c
      do 160 kb = 1, n
         k = n + 1 - kb
         if (cabs1(z(k)) .le. cabs1(abd(m,k))) go to 150
            s = cabs1(abd(m,k))/cabs1(z(k))
            call csscal(n,s,z,1)
            ynorm = s*ynorm
  150    continue
         if (cabs1(abd(m,k)) .ne. 0.0e0) z(k) = z(k)/abd(m,k)
         if (cabs1(abd(m,k)) .eq. 0.0e0) z(k) = (1.0e0,0.0e0)
         lm = min0(k,m) - 1
         la = m - lm
         lz = k - lm
         t = -z(k)
         call caxpy(lm,t,abd(la,k),1,z(lz),1)
  160 continue
c     make znorm = 1.0
      s = 1.0e0/scasum(n,z,1)
      call csscal(n,s,z,1)
      ynorm = s*ynorm
c
      if (anorm .ne. 0.0e0) rcond = ynorm/anorm
      if (anorm .eq. 0.0e0) rcond = 0.0e0
      return
      end
      subroutine cgbdi(abd,lda,n,ml,mu,ipvt,det)

c*********************************************************************72

      integer lda,n,ml,mu,ipvt(1)
      complex abd(lda,1),det(2)
c
c     cgbdi computes the determinant of a band matrix
c     using the factors computed by cgbco or cgbfa.
c     if the inverse is needed, use cgbsl  n  times.
c
c     on entry
c
c        abd     complex(lda, n)
c                the output from cgbco or cgbfa.
c
c        lda     integer
c                the leading dimension of the array  abd .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c
c        mu      integer
c                number of diagonals above the main diagonal.
c
c        ipvt    integer(n)
c                the pivot vector from cgbco or cgbfa.
c
c     on return
c
c        det     complex(2)
c                determinant of original matrix.
c                determinant = det(1) * 10.0**det(2)
c                with  1.0 .le. cabs1(det(1)) .lt. 10.0
c                or  det(1) = 0.0 .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     fortran abs,aimag,cmplx,real
c
c     internal variables
c
      real ten
      integer i,m
c
      complex zdum
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
c
      m = ml + mu + 1
      det(1) = (1.0e0,0.0e0)
      det(2) = (0.0e0,0.0e0)
      ten = 10.0e0
      do 50 i = 1, n
         if (ipvt(i) .ne. i) det(1) = -det(1)
         det(1) = abd(m,i)*det(1)
c     ...exit
         if (cabs1(det(1)) .eq. 0.0e0) go to 60
   10    if (cabs1(det(1)) .ge. 1.0e0) go to 20
            det(1) = cmplx(ten,0.0e0)*det(1)
            det(2) = det(2) - (1.0e0,0.0e0)
         go to 10
   20    continue
   30    if (cabs1(det(1)) .lt. ten) go to 40
            det(1) = det(1)/cmplx(ten,0.0e0)
            det(2) = det(2) + (1.0e0,0.0e0)
         go to 30
   40    continue
   50 continue
   60 continue
      return
      end
      subroutine cgbfa(abd,lda,n,ml,mu,ipvt,info)

c*********************************************************************72

      integer lda,n,ml,mu,ipvt(1),info
      complex abd(lda,1)
c
c     cgbfa factors a complex band matrix by elimination.
c
c     cgbfa is usually called by cgbco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c
c     on entry
c
c        abd     complex(lda, n)
c                contains the matrix in band storage.  the columns
c                of the matrix are stored in the columns of  abd  and
c                the diagonals of the matrix are stored in rows
c                ml+1 through 2*ml+mu+1 of  abd .
c                see the comments below for details.
c
c        lda     integer
c                the leading dimension of the array  abd .
c                lda must be .ge. 2*ml + mu + 1 .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c                0 .le. ml .lt. n .
c
c        mu      integer
c                number of diagonals above the main diagonal.
c                0 .le. mu .lt. n .
c                more efficient if  ml .le. mu .
c     on return
c
c        abd     an upper triangular matrix in band storage and
c                the multipliers which were used to obtain it.
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
c                     indicate that cgbsl will divide by zero if
c                     called.  use  rcond  in cgbco for a reliable
c                     indication of singularity.
c
c     band storage
c
c           if  a  is a band matrix, the following program segment
c           will set up the input.
c
c                   ml = (band width below the diagonal)
c                   mu = (band width above the diagonal)
c                   m = ml + mu + 1
c                   do 20 j = 1, n
c                      i1 = max0(1, j-mu)
c                      i2 = min0(n, j+ml)
c                      do 10 i = i1, i2
c                         k = i - j + m
c                         abd(k,j) = a(i,j)
c                10    continue
c                20 continue
c
c           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
c           in addition, the first  ml  rows in  abd  are used for
c           elements generated during the triangularization.
c           the total number of rows needed in  abd  is  2*ml+mu+1 .
c           the  ml+mu by ml+mu  upper left triangle and the
c           ml by ml  lower right triangle are not referenced.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas caxpy,cscal,icamax
c     fortran abs,aimag,max0,min0,real
c
c     internal variables
c
      complex t
      integer i,icamax,i0,j,ju,jz,j0,j1,k,kp1,l,lm,m,mm,nm1
c
      complex zdum
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
c
      m = ml + mu + 1
      info = 0
c
c     zero initial fill-in columns
c
      j0 = mu + 2
      j1 = min0(n,m) - 1
      if (j1 .lt. j0) go to 30
      do 20 jz = j0, j1
         i0 = m + 1 - jz
         do 10 i = i0, ml
            abd(i,jz) = (0.0e0,0.0e0)
   10    continue
   20 continue
   30 continue
      jz = j1
      ju = 0
c
c     gaussian elimination with partial pivoting
c
      nm1 = n - 1
      if (nm1 .lt. 1) go to 130
      do 120 k = 1, nm1
         kp1 = k + 1
c
c        zero next fill-in column
c
         jz = jz + 1
         if (jz .gt. n) go to 50
         if (ml .lt. 1) go to 50
            do 40 i = 1, ml
               abd(i,jz) = (0.0e0,0.0e0)
   40       continue
   50    continue
c
c        find l = pivot index
c
         lm = min0(ml,n-k)
         l = icamax(lm+1,abd(m,k),1) + m - 1
         ipvt(k) = l + k - m
c
c        zero pivot implies this column already triangularized
c
         if (cabs1(abd(l,k)) .eq. 0.0e0) go to 100
c
c           interchange if necessary
c
            if (l .eq. m) go to 60
               t = abd(l,k)
               abd(l,k) = abd(m,k)
               abd(m,k) = t
   60       continue
c
c           compute multipliers
c
            t = -(1.0e0,0.0e0)/abd(m,k)
            call cscal(lm,t,abd(m+1,k),1)
c
c           row elimination with column indexing
c
            ju = min0(max0(ju,mu+ipvt(k)),n)
            mm = m
            if (ju .lt. kp1) go to 90
            do 80 j = kp1, ju
               l = l - 1
               mm = mm - 1
               t = abd(l,j)
               if (l .eq. mm) go to 70
                  abd(l,j) = abd(mm,j)
                  abd(mm,j) = t
   70          continue
               call caxpy(lm,t,abd(m+1,k),1,abd(mm+1,j),1)
   80       continue
   90       continue
         go to 110
  100    continue
            info = k
  110    continue
  120 continue
  130 continue
      ipvt(n) = n
      if (cabs1(abd(m,n)) .eq. 0.0e0) info = n
      return
      end
      subroutine cgbsl(abd,lda,n,ml,mu,ipvt,b,job)

c*********************************************************************72

      integer lda,n,ml,mu,ipvt(1),job
      complex abd(lda,1),b(1)
c
c     cgbsl solves the complex band system
c     a * x = b  or  ctrans(a) * x = b
c     using the factors computed by cgbco or cgbfa.
c
c     on entry
c
c        abd     complex(lda, n)
c                the output from cgbco or cgbfa.
c
c        lda     integer
c                the leading dimension of the array  abd .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c
c        mu      integer
c                number of diagonals above the main diagonal.
c
c        ipvt    integer(n)
c                the pivot vector from cgbco or cgbfa.
c
c        b       complex(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  ctrans(a)*x = b , where
c                            ctrans(a)  is the conjugate transpose.
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
c        called correctly and if cgbco has set rcond .gt. 0.0
c        or cgbfa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call cgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call cgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas caxpy,cdotc
c     fortran conjg,min0
c
c     internal variables
c
      complex cdotc,t
      integer k,kb,l,la,lb,lm,m,nm1
c
      m = mu + ml + 1
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve l*y = b
c
         if (ml .eq. 0) go to 30
         if (nm1 .lt. 1) go to 30
            do 20 k = 1, nm1
               lm = min0(ml,n-k)
               l = ipvt(k)
               t = b(l)
               if (l .eq. k) go to 10
                  b(l) = b(k)
                  b(k) = t
   10          continue
               call caxpy(lm,t,abd(m+1,k),1,b(k+1),1)
   20       continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/abd(m,k)
            lm = min0(k,m) - 1
            la = m - lm
            lb = k - lm
            t = -b(k)
            call caxpy(lm,t,abd(la,k),1,b(lb),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  ctrans(a) * x = b
c        first solve  ctrans(u)*y = b
c
         do 60 k = 1, n
            lm = min0(k,m) - 1
            la = m - lm
            lb = k - lm
            t = cdotc(lm,abd(la,k),1,b(lb),1)
            b(k) = (b(k) - t)/conjg(abd(m,k))
   60    continue
c
c        now solve ctrans(l)*x = y
c
         if (ml .eq. 0) go to 90
         if (nm1 .lt. 1) go to 90
            do 80 kb = 1, nm1
               k = n - kb
               lm = min0(ml,n-k)
               b(k) = b(k) + cdotc(lm,abd(m+1,k),1,b(k+1),1)
               l = ipvt(k)
               if (l .eq. k) go to 70
                  t = b(l)
                  b(l) = b(k)
                  b(k) = t
   70          continue
   80       continue
   90    continue
  100 continue
      return
      end
      subroutine cgeco(a,lda,n,ipvt,rcond,z)

c*********************************************************************72

      integer lda,n,ipvt(1)
      complex a(lda,1),z(1)
      real rcond
c
c     cgeco factors a complex matrix by gaussian elimination
c     and estimates the condition of the matrix.
c
c     if  rcond  is not needed, cgefa is slightly faster.
c     to solve  a*x = b , follow cgeco by cgesl.
c     to compute  inverse(a)*c , follow cgeco by cgesl.
c     to compute  determinant(a) , follow cgeco by cgedi.
c     to compute  inverse(a) , follow cgeco by cgedi.
c
c     on entry
c
c        a       complex(lda, n)
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
c        rcond   real
c                an estimate of the reciprocal condition of  a .
c                for the system  a*x = b , relative perturbations
c                in  a  and  b  of size  epsilon  may cause
c                relative perturbations in  x  of size  epsilon/rcond .
c                if  rcond  is so small that the logical expression
c                           1.0 + rcond .eq. 1.0
c                is true, then  a  may be singular to working
c                precision.  in particular,  rcond  is zero  if
c                exact singularity is detected or the estimate
c                underflows.
c
c        z       complex(n)
c                a work vector whose contents are usually unimportant.
c                if  a  is close to a singular matrix, then  z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     linpack cgefa
c     blas caxpy,cdotc,csscal,scasum
c     fortran abs,aimag,amax1,cmplx,conjg,real
c
c     internal variables
c
      complex cdotc,ek,t,wk,wkm
      real anorm,s,scasum,sm,ynorm
      integer info,j,k,kb,kp1,l
c
      complex zdum,zdum1,zdum2,csign1
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
      csign1(zdum1,zdum2) = cabs1(zdum1)*(zdum2/cabs1(zdum2))
c
c     compute 1-norm of a
c
      anorm = 0.0e0
      do 10 j = 1, n
         anorm = amax1(anorm,scasum(n,a(1,j),1))
   10 continue
c
c     factor
c
      call cgefa(a,lda,n,ipvt,info)
c
c     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
c     estimate = norm(z)/norm(y) where  a*z = y  and  ctrans(a)*y = e .
c     ctrans(a)  is the conjugate transpose of a .
c     the components of  e  are chosen to cause maximum local
c     growth in the elements of w  where  ctrans(u)*w = e .
c     the vectors are frequently rescaled to avoid overflow.
c
c     solve ctrans(u)*w = e
c
      ek = (1.0e0,0.0e0)
      do 20 j = 1, n
         z(j) = (0.0e0,0.0e0)
   20 continue
      do 100 k = 1, n
         if (cabs1(z(k)) .ne. 0.0e0) ek = csign1(ek,-z(k))
         if (cabs1(ek-z(k)) .le. cabs1(a(k,k))) go to 30
            s = cabs1(a(k,k))/cabs1(ek-z(k))
            call csscal(n,s,z,1)
            ek = cmplx(s,0.0e0)*ek
   30    continue
         wk = ek - z(k)
         wkm = -ek - z(k)
         s = cabs1(wk)
         sm = cabs1(wkm)
         if (cabs1(a(k,k)) .eq. 0.0e0) go to 40
            wk = wk/conjg(a(k,k))
            wkm = wkm/conjg(a(k,k))
         go to 50
   40    continue
            wk = (1.0e0,0.0e0)
            wkm = (1.0e0,0.0e0)
   50    continue
         kp1 = k + 1
         if (kp1 .gt. n) go to 90
            do 60 j = kp1, n
               sm = sm + cabs1(z(j)+wkm*conjg(a(k,j)))
               z(j) = z(j) + wk*conjg(a(k,j))
               s = s + cabs1(z(j))
   60       continue
            if (s .ge. sm) go to 80
               t = wkm - wk
               wk = wkm
               do 70 j = kp1, n
                  z(j) = z(j) + t*conjg(a(k,j))
   70          continue
   80       continue
   90    continue
         z(k) = wk
  100 continue
      s = 1.0e0/scasum(n,z,1)
      call csscal(n,s,z,1)
c
c     solve ctrans(l)*y = w
c
      do 120 kb = 1, n
         k = n + 1 - kb
         if (k .lt. n) z(k) = z(k) + cdotc(n-k,a(k+1,k),1,z(k+1),1)
         if (cabs1(z(k)) .le. 1.0e0) go to 110
            s = 1.0e0/cabs1(z(k))
            call csscal(n,s,z,1)
  110    continue
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
  120 continue
      s = 1.0e0/scasum(n,z,1)
      call csscal(n,s,z,1)
c
      ynorm = 1.0e0
c
c     solve l*v = y
c
      do 140 k = 1, n
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
         if (k .lt. n) call caxpy(n-k,t,a(k+1,k),1,z(k+1),1)
         if (cabs1(z(k)) .le. 1.0e0) go to 130
            s = 1.0e0/cabs1(z(k))
            call csscal(n,s,z,1)
            ynorm = s*ynorm
  130    continue
  140 continue
      s = 1.0e0/scasum(n,z,1)
      call csscal(n,s,z,1)
      ynorm = s*ynorm
c
c     solve  u*z = v
c
      do 160 kb = 1, n
         k = n + 1 - kb
         if (cabs1(z(k)) .le. cabs1(a(k,k))) go to 150
            s = cabs1(a(k,k))/cabs1(z(k))
            call csscal(n,s,z,1)
            ynorm = s*ynorm
  150    continue
         if (cabs1(a(k,k)) .ne. 0.0e0) z(k) = z(k)/a(k,k)
         if (cabs1(a(k,k)) .eq. 0.0e0) z(k) = (1.0e0,0.0e0)
         t = -z(k)
         call caxpy(k-1,t,a(1,k),1,z(1),1)
  160 continue
c     make znorm = 1.0
      s = 1.0e0/scasum(n,z,1)
      call csscal(n,s,z,1)
      ynorm = s*ynorm
c
      if (anorm .ne. 0.0e0) rcond = ynorm/anorm
      if (anorm .eq. 0.0e0) rcond = 0.0e0
      return
      end
      subroutine cgedi(a,lda,n,ipvt,det,work,job)

c*********************************************************************72

      integer lda,n,ipvt(1),job
      complex a(lda,1),det(2),work(1)
c
c     cgedi computes the determinant and inverse of a matrix
c     using the factors computed by cgeco or cgefa.
c
c     on entry
c
c        a       complex(lda, n)
c                the output from cgeco or cgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from cgeco or cgefa.
c
c        work    complex(n)
c                work vector.  contents destroyed.
c
c        job     integer
c                = 11   both determinant and inverse.
c                = 01   inverse only.
c                = 10   determinant only.
c
c     on return
c
c        a       inverse of original matrix if requested.
c                otherwise unchanged.
c
c        det     complex(2)
c                determinant of original matrix if requested.
c                otherwise not referenced.
c                determinant = det(1) * 10.0**det(2)
c                with  1.0 .le. cabs1(det(1)) .lt. 10.0
c                or  det(1) .eq. 0.0 .
c
c     error condition
c
c        a division by zero will occur if the input factor contains
c        a zero on the diagonal and the inverse is requested.
c        it will not occur if the subroutines are called correctly
c        and if cgeco has set rcond .gt. 0.0 or cgefa has set
c        info .eq. 0 .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas caxpy,cscal,cswap
c     fortran abs,aimag,cmplx,mod,real
c
c     internal variables
c
      complex t
      real ten
      integer i,j,k,kb,kp1,l,nm1
c
      complex zdum
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
c
c     compute determinant
c
      if (job/10 .eq. 0) go to 70
         det(1) = (1.0e0,0.0e0)
         det(2) = (0.0e0,0.0e0)
         ten = 10.0e0
         do 50 i = 1, n
            if (ipvt(i) .ne. i) det(1) = -det(1)
            det(1) = a(i,i)*det(1)
c        ...exit
            if (cabs1(det(1)) .eq. 0.0e0) go to 60
   10       if (cabs1(det(1)) .ge. 1.0e0) go to 20
               det(1) = cmplx(ten,0.0e0)*det(1)
               det(2) = det(2) - (1.0e0,0.0e0)
            go to 10
   20       continue
   30       if (cabs1(det(1)) .lt. ten) go to 40
               det(1) = det(1)/cmplx(ten,0.0e0)
               det(2) = det(2) + (1.0e0,0.0e0)
            go to 30
   40       continue
   50    continue
   60    continue
   70 continue
c
c     compute inverse(u)
c
      if (mod(job,10) .eq. 0) go to 150
         do 100 k = 1, n
            a(k,k) = (1.0e0,0.0e0)/a(k,k)
            t = -a(k,k)
            call cscal(k-1,t,a(1,k),1)
            kp1 = k + 1
            if (n .lt. kp1) go to 90
            do 80 j = kp1, n
               t = a(k,j)
               a(k,j) = (0.0e0,0.0e0)
               call caxpy(k,t,a(1,k),1,a(1,j),1)
   80       continue
   90       continue
  100    continue
c
c        form inverse(u)*inverse(l)
c
         nm1 = n - 1
         if (nm1 .lt. 1) go to 140
         do 130 kb = 1, nm1
            k = n - kb
            kp1 = k + 1
            do 110 i = kp1, n
               work(i) = a(i,k)
               a(i,k) = (0.0e0,0.0e0)
  110       continue
            do 120 j = kp1, n
               t = work(j)
               call caxpy(n,t,a(1,j),1,a(1,k),1)
  120       continue
            l = ipvt(k)
            if (l .ne. k) call cswap(n,a(1,k),1,a(1,l),1)
  130    continue
  140    continue
  150 continue
      return
      end
      subroutine cgefa(a,lda,n,ipvt,info)

c*********************************************************************72

      integer lda,n,ipvt(1),info
      complex a(lda,1)
c
c     cgefa factors a complex matrix by gaussian elimination.
c
c     cgefa is usually called by cgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for cgeco) = (1 + 9/n)*(time for cgefa) .
c
c     on entry
c
c        a       complex(lda, n)
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
c                     indicate that cgesl or cgedi will divide by zero
c                     if called.  use  rcond  in cgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas caxpy,cscal,icamax
c     fortran abs,aimag,real
c
c     internal variables
c
      complex t
      integer icamax,j,k,kp1,l,nm1
c
      complex zdum
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
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
         l = icamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (cabs1(a(l,k)) .eq. 0.0e0) go to 40
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
            t = -(1.0e0,0.0e0)/a(k,k)
            call cscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call caxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (cabs1(a(n,n)) .eq. 0.0e0) info = n
      return
      end
      subroutine cgesl(a,lda,n,ipvt,b,job)

c*********************************************************************72

      integer lda,n,ipvt(1),job
      complex a(lda,1),b(1)
c
c     cgesl solves the complex system
c     a * x = b  or  ctrans(a) * x = b
c     using the factors computed by cgeco or cgefa.
c
c     on entry
c
c        a       complex(lda, n)
c                the output from cgeco or cgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from cgeco or cgefa.
c
c        b       complex(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  ctrans(a)*x = b  where
c                            ctrans(a)  is the conjugate transpose.
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
c        called correctly and if cgeco has set rcond .gt. 0.0
c        or cgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call cgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call cgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas caxpy,cdotc
c     fortran conjg
c
c     internal variables
c
      complex cdotc,t
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
            call caxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call caxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  ctrans(a) * x = b
c        first solve  ctrans(u)*y = b
c
         do 60 k = 1, n
            t = cdotc(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/conjg(a(k,k))
   60    continue
c
c        now solve ctrans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + cdotc(n-k,a(k+1,k),1,b(k+1),1)
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
      subroutine cgtsl(n,c,d,e,b,info)

c*********************************************************************72

      integer n,info
      complex c(1),d(1),e(1),b(1)
c
c     cgtsl given a general tridiagonal matrix and a right hand
c     side will find the solution.
c
c     on entry
c
c        n       integer
c                is the order of the tridiagonal matrix.
c
c        c       complex(n)
c                is the subdiagonal of the tridiagonal matrix.
c                c(2) through c(n) should contain the subdiagonal.
c                on output c is destroyed.
c
c        d       complex(n)
c                is the diagonal of the tridiagonal matrix.
c                on output d is destroyed.
c
c        e       complex(n)
c                is the superdiagonal of the tridiagonal matrix.
c                e(1) through e(n-1) should contain the superdiagonal.
c                on output e is destroyed.
c
c        b       complex(n)
c                is the right hand side vector.
c
c     on return
c
c        b       is the solution vector.
c
c        info    integer
c                = 0 normal value.
c                = k if the k-th element of the diagonal becomes
c                    exactly zero.  the subroutine returns when
c                    this is detected.
c
c     linpack. this version dated 08/14/78 .
c     jack dongarra, argonne national laboratory.
c
c     no externals
c     fortran abs,aimag,real
c
c     internal variables
c
      integer k,kb,kp1,nm1,nm2
      complex t
      complex zdum
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
c     begin block permitting ...exits to 100
c
         info = 0
         c(1) = d(1)
         nm1 = n - 1
         if (nm1 .lt. 1) go to 40
            d(1) = e(1)
            e(1) = (0.0e0,0.0e0)
            e(n) = (0.0e0,0.0e0)
c
            do 30 k = 1, nm1
               kp1 = k + 1
c
c              find the largest of the two rows
c
               if (cabs1(c(kp1)) .lt. cabs1(c(k))) go to 10
c
c                 interchange row
c
                  t = c(kp1)
                  c(kp1) = c(k)
                  c(k) = t
                  t = d(kp1)
                  d(kp1) = d(k)
                  d(k) = t
                  t = e(kp1)
                  e(kp1) = e(k)
                  e(k) = t
                  t = b(kp1)
                  b(kp1) = b(k)
                  b(k) = t
   10          continue
c
c              zero elements
c
               if (cabs1(c(k)) .ne. 0.0e0) go to 20
                  info = k
c     ............exit
                  go to 100
   20          continue
               t = -c(kp1)/c(k)
               c(kp1) = d(kp1) + t*d(k)
               d(kp1) = e(kp1) + t*e(k)
               e(kp1) = (0.0e0,0.0e0)
               b(kp1) = b(kp1) + t*b(k)
   30       continue
   40    continue
         if (cabs1(c(n)) .ne. 0.0e0) go to 50
            info = n
         go to 90
   50    continue
c
c           back solve
c
            nm2 = n - 2
            b(n) = b(n)/c(n)
            if (n .eq. 1) go to 80
               b(nm1) = (b(nm1) - d(nm1)*b(n))/c(nm1)
               if (nm2 .lt. 1) go to 70
               do 60 kb = 1, nm2
                  k = nm2 - kb + 1
                  b(k) = (b(k) - d(k)*b(k+1) - e(k)*b(k+2))/c(k)
   60          continue
   70          continue
   80       continue
   90    continue
  100 continue
c
      return
      end
      subroutine chico(a,lda,n,kpvt,rcond,z)

c*********************************************************************72

      integer lda,n,kpvt(1)
      complex a(lda,1),z(1)
      real rcond
c
c     chico factors a complex hermitian matrix by elimination with
c     symmetric pivoting and estimates the condition of the matrix.
c
c     if  rcond  is not needed, chifa is slightly faster.
c     to solve  a*x = b , follow chico by chisl.
c     to compute  inverse(a)*c , follow chico by chisl.
c     to compute  inverse(a) , follow chico by chidi.
c     to compute  determinant(a) , follow chico by chidi.
c     to compute  inertia(a), follow chico by chidi.
c
c     on entry
c
c        a       complex(lda, n)
c                the hermitian matrix to be factored.
c                only the diagonal and upper triangle are used.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     output
c
c        a       a block diagonal matrix and the multipliers which
c                were used to obtain it.
c                the factorization can be written  a = u*d*ctrans(u)
c                where  u  is a product of permutation and unit
c                upper triangular matrices , ctrans(u) is the
c                conjugate transpose of  u , and  d  is block diagonal
c                with 1 by 1 and 2 by 2 blocks.
c
c        kpvt    integer(n)
c                an integer vector of pivot indices.
c
c        rcond   real
c                an estimate of the reciprocal condition of  a .
c                for the system  a*x = b , relative perturbations
c                in  a  and  b  of size  epsilon  may cause
c                relative perturbations in  x  of size  epsilon/rcond .
c                if  rcond  is so small that the logical expression
c                           1.0 + rcond .eq. 1.0
c                is true, then  a  may be singular to working
c                precision.  in particular,  rcond  is zero  if
c                exact singularity is detected or the estimate
c                underflows.
c
c        z       complex(n)
c                a work vector whose contents are usually unimportant.
c                if  a  is close to a singular matrix, then  z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     linpack chifa
c     blas caxpy,cdotc,csscal,scasum
c     fortran abs,aimag,amax1,cmplx,conjg,iabs,real
c
c     internal variables
c
      complex ak,akm1,bk,bkm1,cdotc,denom,ek,t
      real anorm,s,scasum,ynorm
      integer i,info,j,jm1,k,kp,kps,ks
c
      complex zdum,zdum2,csign1
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
      csign1(zdum,zdum2) = cabs1(zdum)*(zdum2/cabs1(zdum2))
c
c     find norm of a using only upper half
c
      do 30 j = 1, n
         z(j) = cmplx(scasum(j,a(1,j),1),0.0e0)
         jm1 = j - 1
         if (jm1 .lt. 1) go to 20
         do 10 i = 1, jm1
            z(i) = cmplx(real(z(i))+cabs1(a(i,j)),0.0e0)
   10    continue
   20    continue
   30 continue
      anorm = 0.0e0
      do 40 j = 1, n
         anorm = amax1(anorm,real(z(j)))
   40 continue
c
c     factor
c
      call chifa(a,lda,n,kpvt,info)
c
c     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
c     estimate = norm(z)/norm(y) where  a*z = y  and  a*y = e .
c     the components of  e  are chosen to cause maximum local
c     growth in the elements of w  where  u*d*w = e .
c     the vectors are frequently rescaled to avoid overflow.
c
c     solve u*d*w = e
c
      ek = (1.0e0,0.0e0)
      do 50 j = 1, n
         z(j) = (0.0e0,0.0e0)
   50 continue
      k = n
   60 if (k .eq. 0) go to 120
         ks = 1
         if (kpvt(k) .lt. 0) ks = 2
         kp = iabs(kpvt(k))
         kps = k + 1 - ks
         if (kp .eq. kps) go to 70
            t = z(kps)
            z(kps) = z(kp)
            z(kp) = t
   70    continue
         if (cabs1(z(k)) .ne. 0.0e0) ek = csign1(ek,z(k))
         z(k) = z(k) + ek
         call caxpy(k-ks,z(k),a(1,k),1,z(1),1)
         if (ks .eq. 1) go to 80
            if (cabs1(z(k-1)) .ne. 0.0e0) ek = csign1(ek,z(k-1))
            z(k-1) = z(k-1) + ek
            call caxpy(k-ks,z(k-1),a(1,k-1),1,z(1),1)
   80    continue
         if (ks .eq. 2) go to 100
            if (cabs1(z(k)) .le. cabs1(a(k,k))) go to 90
               s = cabs1(a(k,k))/cabs1(z(k))
               call csscal(n,s,z,1)
               ek = cmplx(s,0.0e0)*ek
   90       continue
            if (cabs1(a(k,k)) .ne. 0.0e0) z(k) = z(k)/a(k,k)
            if (cabs1(a(k,k)) .eq. 0.0e0) z(k) = (1.0e0,0.0e0)
         go to 110
  100    continue
            ak = a(k,k)/conjg(a(k-1,k))
            akm1 = a(k-1,k-1)/a(k-1,k)
            bk = z(k)/conjg(a(k-1,k))
            bkm1 = z(k-1)/a(k-1,k)
            denom = ak*akm1 - 1.0e0
            z(k) = (akm1*bk - bkm1)/denom
            z(k-1) = (ak*bkm1 - bk)/denom
  110    continue
         k = k - ks
      go to 60
  120 continue
      s = 1.0e0/scasum(n,z,1)
      call csscal(n,s,z,1)
c
c     solve ctrans(u)*y = w
c
      k = 1
  130 if (k .gt. n) go to 160
         ks = 1
         if (kpvt(k) .lt. 0) ks = 2
         if (k .eq. 1) go to 150
            z(k) = z(k) + cdotc(k-1,a(1,k),1,z(1),1)
            if (ks .eq. 2)
     *         z(k+1) = z(k+1) + cdotc(k-1,a(1,k+1),1,z(1),1)
            kp = iabs(kpvt(k))
            if (kp .eq. k) go to 140
               t = z(k)
               z(k) = z(kp)
               z(kp) = t
  140       continue
  150    continue
         k = k + ks
      go to 130
  160 continue
      s = 1.0e0/scasum(n,z,1)
      call csscal(n,s,z,1)
c
      ynorm = 1.0e0
c
c     solve u*d*v = y
c
      k = n
  170 if (k .eq. 0) go to 230
         ks = 1
         if (kpvt(k) .lt. 0) ks = 2
         if (k .eq. ks) go to 190
            kp = iabs(kpvt(k))
            kps = k + 1 - ks
            if (kp .eq. kps) go to 180
               t = z(kps)
               z(kps) = z(kp)
               z(kp) = t
  180       continue
            call caxpy(k-ks,z(k),a(1,k),1,z(1),1)
            if (ks .eq. 2) call caxpy(k-ks,z(k-1),a(1,k-1),1,z(1),1)
  190    continue
         if (ks .eq. 2) go to 210
            if (cabs1(z(k)) .le. cabs1(a(k,k))) go to 200
               s = cabs1(a(k,k))/cabs1(z(k))
               call csscal(n,s,z,1)
               ynorm = s*ynorm
  200       continue
            if (cabs1(a(k,k)) .ne. 0.0e0) z(k) = z(k)/a(k,k)
            if (cabs1(a(k,k)) .eq. 0.0e0) z(k) = (1.0e0,0.0e0)
         go to 220
  210    continue
            ak = a(k,k)/conjg(a(k-1,k))
            akm1 = a(k-1,k-1)/a(k-1,k)
            bk = z(k)/conjg(a(k-1,k))
            bkm1 = z(k-1)/a(k-1,k)
            denom = ak*akm1 - 1.0e0
            z(k) = (akm1*bk - bkm1)/denom
            z(k-1) = (ak*bkm1 - bk)/denom
  220    continue
         k = k - ks
      go to 170
  230 continue
      s = 1.0e0/scasum(n,z,1)
      call csscal(n,s,z,1)
      ynorm = s*ynorm
c
c     solve ctrans(u)*z = v
c
      k = 1
  240 if (k .gt. n) go to 270
         ks = 1
         if (kpvt(k) .lt. 0) ks = 2
         if (k .eq. 1) go to 260
            z(k) = z(k) + cdotc(k-1,a(1,k),1,z(1),1)
            if (ks .eq. 2)
     *         z(k+1) = z(k+1) + cdotc(k-1,a(1,k+1),1,z(1),1)
            kp = iabs(kpvt(k))
            if (kp .eq. k) go to 250
               t = z(k)
               z(k) = z(kp)
               z(kp) = t
  250       continue
  260    continue
         k = k + ks
      go to 240
  270 continue
c     make znorm = 1.0
      s = 1.0e0/scasum(n,z,1)
      call csscal(n,s,z,1)
      ynorm = s*ynorm
c
      if (anorm .ne. 0.0e0) rcond = ynorm/anorm
      if (anorm .eq. 0.0e0) rcond = 0.0e0
      return
      end
      subroutine chidi(a,lda,n,kpvt,det,inert,work,job)

c*********************************************************************72

      integer lda,n,job
      complex a(lda,1),work(1)
      real det(2)
      integer kpvt(1),inert(3)
c
c     chidi computes the determinant, inertia and inverse
c     of a complex hermitian matrix using the factors from chifa.
c
c     on entry
c
c        a       complex(lda,n)
c                the output from chifa.
c
c        lda     integer
c                the leading dimension of the array a.
c
c        n       integer
c                the order of the matrix a.
c
c        kpvt    integer(n)
c                the pivot vector from chifa.
c
c        work    complex(n)
c                work vector.  contents destroyed.
c
c        job     integer
c                job has the decimal expansion  abc  where
c                   if  c .ne. 0, the inverse is computed,
c                   if  b .ne. 0, the determinant is computed,
c                   if  a .ne. 0, the inertia is computed.
c
c                for example, job = 111  gives all three.
c
c     on return
c
c        variables not requested by job are not used.
c
c        a      contains the upper triangle of the inverse of
c               the original matrix.  the strict lower triangle
c               is never referenced.
c
c        det    real(2)
c               determinant of original matrix.
c               determinant = det(1) * 10.0**det(2)
c               with 1.0 .le. abs(det(1)) .lt. 10.0
c               or det(1) = 0.0.
c
c        inert  integer(3)
c               the inertia of the original matrix.
c               inert(1)  =  number of positive eigenvalues.
c               inert(2)  =  number of negative eigenvalues.
c               inert(3)  =  number of zero eigenvalues.
c
c     error condition
c
c        a division by zero may occur if the inverse is requested
c        and  chico  has set rcond .eq. 0.0
c        or  chifa  has set  info .ne. 0 .
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab
c
c     subroutines and functions
c
c     blas caxpy,ccopy,cdotc,cswap
c     fortran abs,cabs,cmplx,conjg,iabs,mod,real
c
c     internal variables.
c
      complex akkp1,cdotc,temp
      real ten,d,t,ak,akp1
      integer j,jb,k,km1,ks,kstep
      logical noinv,nodet,noert
c
      noinv = mod(job,10) .eq. 0
      nodet = mod(job,100)/10 .eq. 0
      noert = mod(job,1000)/100 .eq. 0
c
      if (nodet .and. noert) go to 140
         if (noert) go to 10
            inert(1) = 0
            inert(2) = 0
            inert(3) = 0
   10    continue
         if (nodet) go to 20
            det(1) = 1.0e0
            det(2) = 0.0e0
            ten = 10.0e0
   20    continue
         t = 0.0e0
         do 130 k = 1, n
            d = real(a(k,k))
c
c           check if 1 by 1
c
            if (kpvt(k) .gt. 0) go to 50
c
c              2 by 2 block
c              use det (d  s)  =  (d/t * c - t) * t  ,  t = cabs(s)
c                      (s  c)
c              to avoid underflow/overflow troubles.
c              take two passes through scaling.  use  t  for flag.
c
               if (t .ne. 0.0e0) go to 30
                  t = cabs(a(k,k+1))
                  d = (d/t)*real(a(k+1,k+1)) - t
               go to 40
   30          continue
                  d = t
                  t = 0.0e0
   40          continue
   50       continue
c
            if (noert) go to 60
               if (d .gt. 0.0e0) inert(1) = inert(1) + 1
               if (d .lt. 0.0e0) inert(2) = inert(2) + 1
               if (d .eq. 0.0e0) inert(3) = inert(3) + 1
   60       continue
c
            if (nodet) go to 120
               det(1) = d*det(1)
               if (det(1) .eq. 0.0e0) go to 110
   70             if (abs(det(1)) .ge. 1.0e0) go to 80
                     det(1) = ten*det(1)
                     det(2) = det(2) - 1.0e0
                  go to 70
   80             continue
   90             if (abs(det(1)) .lt. ten) go to 100
                     det(1) = det(1)/ten
                     det(2) = det(2) + 1.0e0
                  go to 90
  100             continue
  110          continue
  120       continue
  130    continue
  140 continue
c
c     compute inverse(a)
c
      if (noinv) go to 270
         k = 1
  150    if (k .gt. n) go to 260
            km1 = k - 1
            if (kpvt(k) .lt. 0) go to 180
c
c              1 by 1
c
               a(k,k) = cmplx(1.0e0/real(a(k,k)),0.0e0)
               if (km1 .lt. 1) go to 170
                  call ccopy(km1,a(1,k),1,work,1)
                  do 160 j = 1, km1
                     a(j,k) = cdotc(j,a(1,j),1,work,1)
                     call caxpy(j-1,work(j),a(1,j),1,a(1,k),1)
  160             continue
                  a(k,k) = a(k,k)
     *                     + cmplx(real(cdotc(km1,work,1,a(1,k),1)),
     *                             0.0e0)
  170          continue
               kstep = 1
            go to 220
  180       continue
c
c              2 by 2
c
               t = cabs(a(k,k+1))
               ak = real(a(k,k))/t
               akp1 = real(a(k+1,k+1))/t
               akkp1 = a(k,k+1)/t
               d = t*(ak*akp1 - 1.0e0)
               a(k,k) = cmplx(akp1/d,0.0e0)
               a(k+1,k+1) = cmplx(ak/d,0.0e0)
               a(k,k+1) = -akkp1/d
               if (km1 .lt. 1) go to 210
                  call ccopy(km1,a(1,k+1),1,work,1)
                  do 190 j = 1, km1
                     a(j,k+1) = cdotc(j,a(1,j),1,work,1)
                     call caxpy(j-1,work(j),a(1,j),1,a(1,k+1),1)
  190             continue
                  a(k+1,k+1) = a(k+1,k+1)
     *                         + cmplx(real(cdotc(km1,work,1,a(1,k+1),
     *                                            1)),0.0e0)
                  a(k,k+1) = a(k,k+1) + cdotc(km1,a(1,k),1,a(1,k+1),1)
                  call ccopy(km1,a(1,k),1,work,1)
                  do 200 j = 1, km1
                     a(j,k) = cdotc(j,a(1,j),1,work,1)
                     call caxpy(j-1,work(j),a(1,j),1,a(1,k),1)
  200             continue
                  a(k,k) = a(k,k)
     *                     + cmplx(real(cdotc(km1,work,1,a(1,k),1)),
     *                             0.0e0)
  210          continue
               kstep = 2
  220       continue
c
c           swap
c
            ks = iabs(kpvt(k))
            if (ks .eq. k) go to 250
               call cswap(ks,a(1,ks),1,a(1,k),1)
               do 230 jb = ks, k
                  j = k + ks - jb
                  temp = conjg(a(j,k))
                  a(j,k) = conjg(a(ks,j))
                  a(ks,j) = temp
  230          continue
               if (kstep .eq. 1) go to 240
                  temp = a(ks,k+1)
                  a(ks,k+1) = a(k,k+1)
                  a(k,k+1) = temp
  240          continue
  250       continue
            k = k + kstep
         go to 150
  260    continue
  270 continue
      return
      end
      subroutine chifa(a,lda,n,kpvt,info)

c*********************************************************************72

      integer lda,n,kpvt(1),info
      complex a(lda,1)
c
c     chifa factors a complex hermitian matrix by elimination
c     with symmetric pivoting.
c
c     to solve  a*x = b , follow chifa by chisl.
c     to compute  inverse(a)*c , follow chifa by chisl.
c     to compute  determinant(a) , follow chifa by chidi.
c     to compute  inertia(a) , follow chifa by chidi.
c     to compute  inverse(a) , follow chifa by chidi.
c
c     on entry
c
c        a       complex(lda,n)
c                the hermitian matrix to be factored.
c                only the diagonal and upper triangle are used.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       a block diagonal matrix and the multipliers which
c                were used to obtain it.
c                the factorization can be written  a = u*d*ctrans(u)
c                where  u  is a product of permutation and unit
c                upper triangular matrices , ctrans(u) is the
c                conjugate transpose of  u , and  d  is block diagonal
c                with 1 by 1 and 2 by 2 blocks.
c
c        kpvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if the k-th pivot block is singular. this is
c                     not an error condition for this subroutine,
c                     but it does indicate that chisl or chidi may
c                     divide by zero if called.
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab.
c
c     subroutines and functions
c
c     blas caxpy,cswap,icamax
c     fortran abs,aimag,amax1,cmplx,conjg,real,sqrt
c
c     internal variables
c
      complex ak,akm1,bk,bkm1,denom,mulk,mulkm1,t
      real absakk,alpha,colmax,rowmax
      integer imax,imaxp1,j,jj,jmax,k,km1,km2,kstep,icamax
      logical swap
c
      complex zdum
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
c
c     initialize
c
c     alpha is used in choosing pivot block size.
      alpha = (1.0e0 + sqrt(17.0e0))/8.0e0
c
      info = 0
c
c     main loop on k, which goes from n to 1.
c
      k = n
   10 continue
c
c        leave the loop if k=0 or k=1.
c
c     ...exit
         if (k .eq. 0) go to 200
         if (k .gt. 1) go to 20
            kpvt(1) = 1
            if (cabs1(a(1,1)) .eq. 0.0e0) info = 1
c     ......exit
            go to 200
   20    continue
c
c        this section of code determines the kind of
c        elimination to be performed.  when it is completed,
c        kstep will be set to the size of the pivot block, and
c        swap will be set to .true. if an interchange is
c        required.
c
         km1 = k - 1
         absakk = cabs1(a(k,k))
c
c        determine the largest off-diagonal element in
c        column k.
c
         imax = icamax(k-1,a(1,k),1)
         colmax = cabs1(a(imax,k))
         if (absakk .lt. alpha*colmax) go to 30
            kstep = 1
            swap = .false.
         go to 90
   30    continue
c
c           determine the largest off-diagonal element in
c           row imax.
c
            rowmax = 0.0e0
            imaxp1 = imax + 1
            do 40 j = imaxp1, k
               rowmax = amax1(rowmax,cabs1(a(imax,j)))
   40       continue
            if (imax .eq. 1) go to 50
               jmax = icamax(imax-1,a(1,imax),1)
               rowmax = amax1(rowmax,cabs1(a(jmax,imax)))
   50       continue
            if (cabs1(a(imax,imax)) .lt. alpha*rowmax) go to 60
               kstep = 1
               swap = .true.
            go to 80
   60       continue
            if (absakk .lt. alpha*colmax*(colmax/rowmax)) go to 70
               kstep = 1
               swap = .false.
            go to 80
   70       continue
               kstep = 2
               swap = imax .ne. km1
   80       continue
   90    continue
         if (amax1(absakk,colmax) .ne. 0.0e0) go to 100
c
c           column k is zero.  set info and iterate the loop.
c
            kpvt(k) = k
            info = k
         go to 190
  100    continue
         if (kstep .eq. 2) go to 140
c
c           1 x 1 pivot block.
c
            if (.not.swap) go to 120
c
c              perform an interchange.
c
               call cswap(imax,a(1,imax),1,a(1,k),1)
               do 110 jj = imax, k
                  j = k + imax - jj
                  t = conjg(a(j,k))
                  a(j,k) = conjg(a(imax,j))
                  a(imax,j) = t
  110          continue
  120       continue
c
c           perform the elimination.
c
            do 130 jj = 1, km1
               j = k - jj
               mulk = -a(j,k)/a(k,k)
               t = conjg(mulk)
               call caxpy(j,t,a(1,k),1,a(1,j),1)
               a(j,j) = cmplx(real(a(j,j)),0.0e0)
               a(j,k) = mulk
  130       continue
c
c           set the pivot array.
c
            kpvt(k) = k
            if (swap) kpvt(k) = imax
         go to 190
  140    continue
c
c           2 x 2 pivot block.
c
            if (.not.swap) go to 160
c
c              perform an interchange.
c
               call cswap(imax,a(1,imax),1,a(1,k-1),1)
               do 150 jj = imax, km1
                  j = km1 + imax - jj
                  t = conjg(a(j,k-1))
                  a(j,k-1) = conjg(a(imax,j))
                  a(imax,j) = t
  150          continue
               t = a(k-1,k)
               a(k-1,k) = a(imax,k)
               a(imax,k) = t
  160       continue
c
c           perform the elimination.
c
            km2 = k - 2
            if (km2 .eq. 0) go to 180
               ak = a(k,k)/a(k-1,k)
               akm1 = a(k-1,k-1)/conjg(a(k-1,k))
               denom = 1.0e0 - ak*akm1
               do 170 jj = 1, km2
                  j = km1 - jj
                  bk = a(j,k)/a(k-1,k)
                  bkm1 = a(j,k-1)/conjg(a(k-1,k))
                  mulk = (akm1*bk - bkm1)/denom
                  mulkm1 = (ak*bkm1 - bk)/denom
                  t = conjg(mulk)
                  call caxpy(j,t,a(1,k),1,a(1,j),1)
                  t = conjg(mulkm1)
                  call caxpy(j,t,a(1,k-1),1,a(1,j),1)
                  a(j,k) = mulk
                  a(j,k-1) = mulkm1
                  a(j,j) = cmplx(real(a(j,j)),0.0e0)
  170          continue
  180       continue
c
c           set the pivot array.
c
            kpvt(k) = 1 - k
            if (swap) kpvt(k) = -imax
            kpvt(k-1) = kpvt(k)
  190    continue
         k = k - kstep
      go to 10
  200 continue
      return
      end
      subroutine chisl(a,lda,n,kpvt,b)

c*********************************************************************72

      integer lda,n,kpvt(1)
      complex a(lda,1),b(1)
c
c     chisl solves the complex hermitian system
c     a * x = b
c     using the factors computed by chifa.
c
c     on entry
c
c        a       complex(lda,n)
c                the output from chifa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        kpvt    integer(n)
c                the pivot vector from chifa.
c
c        b       complex(n)
c                the right hand side vector.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero may occur if  chico  has set rcond .eq. 0.0
c        or  chifa  has set info .ne. 0  .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call chifa(a,lda,n,kpvt,info)
c           if (info .ne. 0) go to ...
c           do 10 j = 1, p
c              call chisl(a,lda,n,kpvt,c(1,j))
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab.
c
c     subroutines and functions
c
c     blas caxpy,cdotc
c     fortran conjg,iabs
c
c     internal variables.
c
      complex ak,akm1,bk,bkm1,cdotc,denom,temp
      integer k,kp
c
c     loop backward applying the transformations and
c     d inverse to b.
c
      k = n
   10 if (k .eq. 0) go to 80
         if (kpvt(k) .lt. 0) go to 40
c
c           1 x 1 pivot block.
c
            if (k .eq. 1) go to 30
               kp = kpvt(k)
               if (kp .eq. k) go to 20
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
   20          continue
c
c              apply the transformation.
c
               call caxpy(k-1,b(k),a(1,k),1,b(1),1)
   30       continue
c
c           apply d inverse.
c
            b(k) = b(k)/a(k,k)
            k = k - 1
         go to 70
   40    continue
c
c           2 x 2 pivot block.
c
            if (k .eq. 2) go to 60
               kp = iabs(kpvt(k))
               if (kp .eq. k - 1) go to 50
c
c                 interchange.
c
                  temp = b(k-1)
                  b(k-1) = b(kp)
                  b(kp) = temp
   50          continue
c
c              apply the transformation.
c
               call caxpy(k-2,b(k),a(1,k),1,b(1),1)
               call caxpy(k-2,b(k-1),a(1,k-1),1,b(1),1)
   60       continue
c
c           apply d inverse.
c
            ak = a(k,k)/conjg(a(k-1,k))
            akm1 = a(k-1,k-1)/a(k-1,k)
            bk = b(k)/conjg(a(k-1,k))
            bkm1 = b(k-1)/a(k-1,k)
            denom = ak*akm1 - 1.0e0
            b(k) = (akm1*bk - bkm1)/denom
            b(k-1) = (ak*bkm1 - bk)/denom
            k = k - 2
   70    continue
      go to 10
   80 continue
c
c     loop forward applying the transformations.
c
      k = 1
   90 if (k .gt. n) go to 160
         if (kpvt(k) .lt. 0) go to 120
c
c           1 x 1 pivot block.
c
            if (k .eq. 1) go to 110
c
c              apply the transformation.
c
               b(k) = b(k) + cdotc(k-1,a(1,k),1,b(1),1)
               kp = kpvt(k)
               if (kp .eq. k) go to 100
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
  100          continue
  110       continue
            k = k + 1
         go to 150
  120    continue
c
c           2 x 2 pivot block.
c
            if (k .eq. 1) go to 140
c
c              apply the transformation.
c
               b(k) = b(k) + cdotc(k-1,a(1,k),1,b(1),1)
               b(k+1) = b(k+1) + cdotc(k-1,a(1,k+1),1,b(1),1)
               kp = iabs(kpvt(k))
               if (kp .eq. k) go to 130
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
  130          continue
  140       continue
            k = k + 2
  150    continue
      go to 90
  160 continue
      return
      end
      subroutine chpco(ap,n,kpvt,rcond,z)

c*********************************************************************72

      integer n,kpvt(1)
      complex ap(1),z(1)
      real rcond
c
c     chpco factors a complex hermitian matrix stored in packed
c     form by elimination with symmetric pivoting and estimates
c     the condition of the matrix.
c
c     if  rcond  is not needed, chpfa is slightly faster.
c     to solve  a*x = b , follow chpco by chpsl.
c     to compute  inverse(a)*c , follow chpco by chpsl.
c     to compute  inverse(a) , follow chpco by chpdi.
c     to compute  determinant(a) , follow chpco by chpdi.
c     to compute  inertia(a), follow chpco by chpdi.
c
c     on entry
c
c        ap      complex (n*(n+1)/2)
c                the packed form of a hermitian matrix  a .  the
c                columns of the upper triangle are stored sequentially
c                in a one-dimensional array of length  n*(n+1)/2 .
c                see comments below for details.
c
c        n       integer
c                the order of the matrix  a .
c
c     output
c
c        ap      a block diagonal matrix and the multipliers which
c                were used to obtain it stored in packed form.
c                the factorization can be written  a = u*d*ctrans(u)
c                where  u  is a product of permutation and unit
c                upper triangular matrices , ctrans(u) is the
c                conjugate transpose of  u , and  d  is block diagonal
c                with 1 by 1 and 2 by 2 blocks.
c
c        kpvt    integer(n)
c                an integer vector of pivot indices.
c
c        rcond   real
c                an estimate of the reciprocal condition of  a .
c                for the system  a*x = b , relative perturbations
c                in  a  and  b  of size  epsilon  may cause
c                relative perturbations in  x  of size  epsilon/rcond .
c                if  rcond  is so small that the logical expression
c                           1.0 + rcond .eq. 1.0
c                is true, then  a  may be singular to working
c                precision.  in particular,  rcond  is zero  if
c                exact singularity is detected or the estimate
c                underflows.
c
c        z       complex(n)
c                a work vector whose contents are usually unimportant.
c                if  a  is close to a singular matrix, then  z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .
c
c     packed storage
c
c          the following program segment will pack the upper
c          triangle of a hermitian matrix.
c
c                k = 0
c                do 20 j = 1, n
c                   do 10 i = 1, j
c                      k = k + 1
c                      ap(k) = a(i,j)
c             10    continue
c             20 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     linpack chpfa
c     blas caxpy,cdotc,csscal,scasum
c     fortran abs,aimag,amax1,cmplx,conjg,iabs,real
c
c     internal variables
c
      complex ak,akm1,bk,bkm1,cdotc,denom,ek,t
      real anorm,s,scasum,ynorm
      integer i,ij,ik,ikm1,ikp1,info,j,jm1,j1
      integer k,kk,km1k,km1km1,kp,kps,ks
c
      complex zdum,zdum2,csign1
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
      csign1(zdum,zdum2) = cabs1(zdum)*(zdum2/cabs1(zdum2))
c
c     find norm of a using only upper half
c
      j1 = 1
      do 30 j = 1, n
         z(j) = cmplx(scasum(j,ap(j1),1),0.0e0)
         ij = j1
         j1 = j1 + j
         jm1 = j - 1
         if (jm1 .lt. 1) go to 20
         do 10 i = 1, jm1
            z(i) = cmplx(real(z(i))+cabs1(ap(ij)),0.0e0)
            ij = ij + 1
   10    continue
   20    continue
   30 continue
      anorm = 0.0e0
      do 40 j = 1, n
         anorm = amax1(anorm,real(z(j)))
   40 continue
c
c     factor
c
      call chpfa(ap,n,kpvt,info)
c
c     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
c     estimate = norm(z)/norm(y) where  a*z = y  and  a*y = e .
c     the components of  e  are chosen to cause maximum local
c     growth in the elements of w  where  u*d*w = e .
c     the vectors are frequently rescaled to avoid overflow.
c
c     solve u*d*w = e
c
      ek = (1.0e0,0.0e0)
      do 50 j = 1, n
         z(j) = (0.0e0,0.0e0)
   50 continue
      k = n
      ik = (n*(n - 1))/2
   60 if (k .eq. 0) go to 120
         kk = ik + k
         ikm1 = ik - (k - 1)
         ks = 1
         if (kpvt(k) .lt. 0) ks = 2
         kp = iabs(kpvt(k))
         kps = k + 1 - ks
         if (kp .eq. kps) go to 70
            t = z(kps)
            z(kps) = z(kp)
            z(kp) = t
   70    continue
         if (cabs1(z(k)) .ne. 0.0e0) ek = csign1(ek,z(k))
         z(k) = z(k) + ek
         call caxpy(k-ks,z(k),ap(ik+1),1,z(1),1)
         if (ks .eq. 1) go to 80
            if (cabs1(z(k-1)) .ne. 0.0e0) ek = csign1(ek,z(k-1))
            z(k-1) = z(k-1) + ek
            call caxpy(k-ks,z(k-1),ap(ikm1+1),1,z(1),1)
   80    continue
         if (ks .eq. 2) go to 100
            if (cabs1(z(k)) .le. cabs1(ap(kk))) go to 90
               s = cabs1(ap(kk))/cabs1(z(k))
               call csscal(n,s,z,1)
               ek = cmplx(s,0.0e0)*ek
   90       continue
            if (cabs1(ap(kk)) .ne. 0.0e0) z(k) = z(k)/ap(kk)
            if (cabs1(ap(kk)) .eq. 0.0e0) z(k) = (1.0e0,0.0e0)
         go to 110
  100    continue
            km1k = ik + k - 1
            km1km1 = ikm1 + k - 1
            ak = ap(kk)/conjg(ap(km1k))
            akm1 = ap(km1km1)/ap(km1k)
            bk = z(k)/conjg(ap(km1k))
            bkm1 = z(k-1)/ap(km1k)
            denom = ak*akm1 - 1.0e0
            z(k) = (akm1*bk - bkm1)/denom
            z(k-1) = (ak*bkm1 - bk)/denom
  110    continue
         k = k - ks
         ik = ik - k
         if (ks .eq. 2) ik = ik - (k + 1)
      go to 60
  120 continue
      s = 1.0e0/scasum(n,z,1)
      call csscal(n,s,z,1)
c
c     solve ctrans(u)*y = w
c
      k = 1
      ik = 0
  130 if (k .gt. n) go to 160
         ks = 1
         if (kpvt(k) .lt. 0) ks = 2
         if (k .eq. 1) go to 150
            z(k) = z(k) + cdotc(k-1,ap(ik+1),1,z(1),1)
            ikp1 = ik + k
            if (ks .eq. 2)
     *         z(k+1) = z(k+1) + cdotc(k-1,ap(ikp1+1),1,z(1),1)
            kp = iabs(kpvt(k))
            if (kp .eq. k) go to 140
               t = z(k)
               z(k) = z(kp)
               z(kp) = t
  140       continue
  150    continue
         ik = ik + k
         if (ks .eq. 2) ik = ik + (k + 1)
         k = k + ks
      go to 130
  160 continue
      s = 1.0e0/scasum(n,z,1)
      call csscal(n,s,z,1)
c
      ynorm = 1.0e0
c
c     solve u*d*v = y
c
      k = n
      ik = n*(n - 1)/2
  170 if (k .eq. 0) go to 230
         kk = ik + k
         ikm1 = ik - (k - 1)
         ks = 1
         if (kpvt(k) .lt. 0) ks = 2
         if (k .eq. ks) go to 190
            kp = iabs(kpvt(k))
            kps = k + 1 - ks
            if (kp .eq. kps) go to 180
               t = z(kps)
               z(kps) = z(kp)
               z(kp) = t
  180       continue
            call caxpy(k-ks,z(k),ap(ik+1),1,z(1),1)
            if (ks .eq. 2) call caxpy(k-ks,z(k-1),ap(ikm1+1),1,z(1),1)
  190    continue
         if (ks .eq. 2) go to 210
            if (cabs1(z(k)) .le. cabs1(ap(kk))) go to 200
               s = cabs1(ap(kk))/cabs1(z(k))
               call csscal(n,s,z,1)
               ynorm = s*ynorm
  200       continue
            if (cabs1(ap(kk)) .ne. 0.0e0) z(k) = z(k)/ap(kk)
            if (cabs1(ap(kk)) .eq. 0.0e0) z(k) = (1.0e0,0.0e0)
         go to 220
  210    continue
            km1k = ik + k - 1
            km1km1 = ikm1 + k - 1
            ak = ap(kk)/conjg(ap(km1k))
            akm1 = ap(km1km1)/ap(km1k)
            bk = z(k)/conjg(ap(km1k))
            bkm1 = z(k-1)/ap(km1k)
            denom = ak*akm1 - 1.0e0
            z(k) = (akm1*bk - bkm1)/denom
            z(k-1) = (ak*bkm1 - bk)/denom
  220    continue
         k = k - ks
         ik = ik - k
         if (ks .eq. 2) ik = ik - (k + 1)
      go to 170
  230 continue
      s = 1.0e0/scasum(n,z,1)
      call csscal(n,s,z,1)
      ynorm = s*ynorm
c
c     solve ctrans(u)*z = v
c
      k = 1
      ik = 0
  240 if (k .gt. n) go to 270
         ks = 1
         if (kpvt(k) .lt. 0) ks = 2
         if (k .eq. 1) go to 260
            z(k) = z(k) + cdotc(k-1,ap(ik+1),1,z(1),1)
            ikp1 = ik + k
            if (ks .eq. 2)
     *         z(k+1) = z(k+1) + cdotc(k-1,ap(ikp1+1),1,z(1),1)
            kp = iabs(kpvt(k))
            if (kp .eq. k) go to 250
               t = z(k)
               z(k) = z(kp)
               z(kp) = t
  250       continue
  260    continue
         ik = ik + k
         if (ks .eq. 2) ik = ik + (k + 1)
         k = k + ks
      go to 240
  270 continue
c     make znorm = 1.0
      s = 1.0e0/scasum(n,z,1)
      call csscal(n,s,z,1)
      ynorm = s*ynorm
c
      if (anorm .ne. 0.0e0) rcond = ynorm/anorm
      if (anorm .eq. 0.0e0) rcond = 0.0e0
      return
      end
      subroutine chpdi(ap,n,kpvt,det,inert,work,job)

c*********************************************************************72

      integer n,job
      complex ap(1),work(1)
      real det(2)
      integer kpvt(1),inert(3)
c
c     chpdi computes the determinant, inertia and inverse
c     of a complex hermitian matrix using the factors from chpfa,
c     where the matrix is stored in packed form.
c
c     on entry
c
c        ap      complex (n*(n+1)/2)
c                the output from chpfa.
c
c        n       integer
c                the order of the matrix a.
c
c        kpvt    integer(n)
c                the pivot vector from chpfa.
c
c        work    complex(n)
c                work vector.  contents ignored.
c
c        job     integer
c                job has the decimal expansion  abc  where
c                   if  c .ne. 0, the inverse is computed,
c                   if  b .ne. 0, the determinant is computed,
c                   if  a .ne. 0, the inertia is computed.
c
c                for example, job = 111  gives all three.
c
c     on return
c
c        variables not requested by job are not used.
c
c        ap     contains the upper triangle of the inverse of
c               the original matrix, stored in packed form.
c               the columns of the upper triangle are stored
c               sequentially in a one-dimensional array.
c
c        det    real(2)
c               determinant of original matrix.
c               determinant = det(1) * 10.0**det(2)
c               with 1.0 .le. abs(det(1)) .lt. 10.0
c               or det(1) = 0.0.
c
c        inert  integer(3)
c               the inertia of the original matrix.
c               inert(1)  =  number of positive eigenvalues.
c               inert(2)  =  number of negative eigenvalues.
c               inert(3)  =  number of zero eigenvalues.
c
c     error condition
c
c        a division by zero will occur if the inverse is requested
c        and  chpco  has set rcond .eq. 0.0
c        or  chpfa  has set  info .ne. 0 .
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab.
c
c     subroutines and functions
c
c     blas caxpy,ccopy,cdotc,cswap
c     fortran abs,cabs,cmplx,conjg,iabs,mod,real
c
c     internal variables.
c
      complex akkp1,cdotc,temp
      real ten,d,t,ak,akp1
      integer ij,ik,ikp1,iks,j,jb,jk,jkp1
      integer k,kk,kkp1,km1,ks,ksj,kskp1,kstep
      logical noinv,nodet,noert
c
      noinv = mod(job,10) .eq. 0
      nodet = mod(job,100)/10 .eq. 0
      noert = mod(job,1000)/100 .eq. 0
c
      if (nodet .and. noert) go to 140
         if (noert) go to 10
            inert(1) = 0
            inert(2) = 0
            inert(3) = 0
   10    continue
         if (nodet) go to 20
            det(1) = 1.0e0
            det(2) = 0.0e0
            ten = 10.0e0
   20    continue
         t = 0.0e0
         ik = 0
         do 130 k = 1, n
            kk = ik + k
            d = real(ap(kk))
c
c           check if 1 by 1
c
            if (kpvt(k) .gt. 0) go to 50
c
c              2 by 2 block
c              use det (d  s)  =  (d/t * c - t) * t  ,  t = cabs(s)
c                      (s  c)
c              to avoid underflow/overflow troubles.
c              take two passes through scaling.  use  t  for flag.
c
               if (t .ne. 0.0e0) go to 30
                  ikp1 = ik + k
                  kkp1 = ikp1 + k
                  t = cabs(ap(kkp1))
                  d = (d/t)*real(ap(kkp1+1)) - t
               go to 40
   30          continue
                  d = t
                  t = 0.0e0
   40          continue
   50       continue
c
            if (noert) go to 60
               if (d .gt. 0.0e0) inert(1) = inert(1) + 1
               if (d .lt. 0.0e0) inert(2) = inert(2) + 1
               if (d .eq. 0.0e0) inert(3) = inert(3) + 1
   60       continue
c
            if (nodet) go to 120
               det(1) = d*det(1)
               if (det(1) .eq. 0.0e0) go to 110
   70             if (abs(det(1)) .ge. 1.0e0) go to 80
                     det(1) = ten*det(1)
                     det(2) = det(2) - 1.0e0
                  go to 70
   80             continue
   90             if (abs(det(1)) .lt. ten) go to 100
                     det(1) = det(1)/ten
                     det(2) = det(2) + 1.0e0
                  go to 90
  100             continue
  110          continue
  120       continue
            ik = ik + k
  130    continue
  140 continue
c
c     compute inverse(a)
c
      if (noinv) go to 270
         k = 1
         ik = 0
  150    if (k .gt. n) go to 260
            km1 = k - 1
            kk = ik + k
            ikp1 = ik + k
            kkp1 = ikp1 + k
            if (kpvt(k) .lt. 0) go to 180
c
c              1 by 1
c
               ap(kk) = cmplx(1.0e0/real(ap(kk)),0.0e0)
               if (km1 .lt. 1) go to 170
                  call ccopy(km1,ap(ik+1),1,work,1)
                  ij = 0
                  do 160 j = 1, km1
                     jk = ik + j
                     ap(jk) = cdotc(j,ap(ij+1),1,work,1)
                     call caxpy(j-1,work(j),ap(ij+1),1,ap(ik+1),1)
                     ij = ij + j
  160             continue
                  ap(kk) = ap(kk)
     *                     + cmplx(real(cdotc(km1,work,1,ap(ik+1),1)),
     *                             0.0e0)
  170          continue
               kstep = 1
            go to 220
  180       continue
c
c              2 by 2
c
               t = cabs(ap(kkp1))
               ak = real(ap(kk))/t
               akp1 = real(ap(kkp1+1))/t
               akkp1 = ap(kkp1)/t
               d = t*(ak*akp1 - 1.0e0)
               ap(kk) = cmplx(akp1/d,0.0e0)
               ap(kkp1+1) = cmplx(ak/d,0.0e0)
               ap(kkp1) = -akkp1/d
               if (km1 .lt. 1) go to 210
                  call ccopy(km1,ap(ikp1+1),1,work,1)
                  ij = 0
                  do 190 j = 1, km1
                     jkp1 = ikp1 + j
                     ap(jkp1) = cdotc(j,ap(ij+1),1,work,1)
                     call caxpy(j-1,work(j),ap(ij+1),1,ap(ikp1+1),1)
                     ij = ij + j
  190             continue
                  ap(kkp1+1) = ap(kkp1+1)
     *                         + cmplx(real(cdotc(km1,work,1,
     *                                            ap(ikp1+1),1)),0.0e0)
                  ap(kkp1) = ap(kkp1)
     *                       + cdotc(km1,ap(ik+1),1,ap(ikp1+1),1)
                  call ccopy(km1,ap(ik+1),1,work,1)
                  ij = 0
                  do 200 j = 1, km1
                     jk = ik + j
                     ap(jk) = cdotc(j,ap(ij+1),1,work,1)
                     call caxpy(j-1,work(j),ap(ij+1),1,ap(ik+1),1)
                     ij = ij + j
  200             continue
                  ap(kk) = ap(kk)
     *                     + cmplx(real(cdotc(km1,work,1,ap(ik+1),1)),
     *                             0.0e0)
  210          continue
               kstep = 2
  220       continue
c
c           swap
c
            ks = iabs(kpvt(k))
            if (ks .eq. k) go to 250
               iks = (ks*(ks - 1))/2
               call cswap(ks,ap(iks+1),1,ap(ik+1),1)
               ksj = ik + ks
               do 230 jb = ks, k
                  j = k + ks - jb
                  jk = ik + j
                  temp = conjg(ap(jk))
                  ap(jk) = conjg(ap(ksj))
                  ap(ksj) = temp
                  ksj = ksj - (j - 1)
  230          continue
               if (kstep .eq. 1) go to 240
                  kskp1 = ikp1 + ks
                  temp = ap(kskp1)
                  ap(kskp1) = ap(kkp1)
                  ap(kkp1) = temp
  240          continue
  250       continue
            ik = ik + k
            if (kstep .eq. 2) ik = ik + k + 1
            k = k + kstep
         go to 150
  260    continue
  270 continue
      return
      end
      subroutine chpfa(ap,n,kpvt,info)

c*********************************************************************72

      integer n,kpvt(1),info
      complex ap(1)
c
c     chpfa factors a complex hermitian matrix stored in
c     packed form by elimination with symmetric pivoting.
c
c     to solve  a*x = b , follow chpfa by chpsl.
c     to compute  inverse(a)*c , follow chpfa by chpsl.
c     to compute  determinant(a) , follow chpfa by chpdi.
c     to compute  inertia(a) , follow chpfa by chpdi.
c     to compute  inverse(a) , follow chpfa by chpdi.
c
c     on entry
c
c        ap      complex (n*(n+1)/2)
c                the packed form of a hermitian matrix  a .  the
c                columns of the upper triangle are stored sequentially
c                in a one-dimensional array of length  n*(n+1)/2 .
c                see comments below for details.
c
c        n       integer
c                the order of the matrix  a .
c
c     output
c
c        ap      a block diagonal matrix and the multipliers which
c                were used to obtain it stored in packed form.
c                the factorization can be written  a = u*d*ctrans(u)
c                where  u  is a product of permutation and unit
c                upper triangular matrices , ctrans(u) is the
c                conjugate transpose of  u , and  d  is block diagonal
c                with 1 by 1 and 2 by 2 blocks.
c
c        kpvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if the k-th pivot block is singular. this is
c                     not an error condition for this subroutine,
c                     but it does indicate that chpsl or chpdi may
c                     divide by zero if called.
c
c     packed storage
c
c          the following program segment will pack the upper
c          triangle of a hermitian matrix.
c
c                k = 0
c                do 20 j = 1, n
c                   do 10 i = 1, j
c                      k = k + 1
c                      ap(k)  = a(i,j)
c             10    continue
c             20 continue
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab.
c
c     subroutines and functions
c
c     blas caxpy,cswap,icamax
c     fortran abs,aimag,amax1,cmplx,conjg,real,sqrt
c
c     internal variables
c
      complex ak,akm1,bk,bkm1,denom,mulk,mulkm1,t
      real absakk,alpha,colmax,rowmax
      integer icamax,ij,ijj,ik,ikm1,im,imax,imaxp1,imim,imj,imk
      integer j,jj,jk,jkm1,jmax,jmim,k,kk,km1,km1k,km1km1,km2,kstep
      logical swap
c
      complex zdum
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
c
c     initialize
c
c     alpha is used in choosing pivot block size.
      alpha = (1.0e0 + sqrt(17.0e0))/8.0e0
c
      info = 0
c
c     main loop on k, which goes from n to 1.
c
      k = n
      ik = (n*(n - 1))/2
   10 continue
c
c        leave the loop if k=0 or k=1.
c
c     ...exit
         if (k .eq. 0) go to 200
         if (k .gt. 1) go to 20
            kpvt(1) = 1
            if (cabs1(ap(1)) .eq. 0.0e0) info = 1
c     ......exit
            go to 200
   20    continue
c
c        this section of code determines the kind of
c        elimination to be performed.  when it is completed,
c        kstep will be set to the size of the pivot block, and
c        swap will be set to .true. if an interchange is
c        required.
c
         km1 = k - 1
         kk = ik + k
         absakk = cabs1(ap(kk))
c
c        determine the largest off-diagonal element in
c        column k.
c
         imax = icamax(k-1,ap(ik+1),1)
         imk = ik + imax
         colmax = cabs1(ap(imk))
         if (absakk .lt. alpha*colmax) go to 30
            kstep = 1
            swap = .false.
         go to 90
   30    continue
c
c           determine the largest off-diagonal element in
c           row imax.
c
            rowmax = 0.0e0
            imaxp1 = imax + 1
            im = imax*(imax - 1)/2
            imj = im + 2*imax
            do 40 j = imaxp1, k
               rowmax = amax1(rowmax,cabs1(ap(imj)))
               imj = imj + j
   40       continue
            if (imax .eq. 1) go to 50
               jmax = icamax(imax-1,ap(im+1),1)
               jmim = jmax + im
               rowmax = amax1(rowmax,cabs1(ap(jmim)))
   50       continue
            imim = imax + im
            if (cabs1(ap(imim)) .lt. alpha*rowmax) go to 60
               kstep = 1
               swap = .true.
            go to 80
   60       continue
            if (absakk .lt. alpha*colmax*(colmax/rowmax)) go to 70
               kstep = 1
               swap = .false.
            go to 80
   70       continue
               kstep = 2
               swap = imax .ne. km1
   80       continue
   90    continue
         if (amax1(absakk,colmax) .ne. 0.0e0) go to 100
c
c           column k is zero.  set info and iterate the loop.
c
            kpvt(k) = k
            info = k
         go to 190
  100    continue
         if (kstep .eq. 2) go to 140
c
c           1 x 1 pivot block.
c
            if (.not.swap) go to 120
c
c              perform an interchange.
c
               call cswap(imax,ap(im+1),1,ap(ik+1),1)
               imj = ik + imax
               do 110 jj = imax, k
                  j = k + imax - jj
                  jk = ik + j
                  t = conjg(ap(jk))
                  ap(jk) = conjg(ap(imj))
                  ap(imj) = t
                  imj = imj - (j - 1)
  110          continue
  120       continue
c
c           perform the elimination.
c
            ij = ik - (k - 1)
            do 130 jj = 1, km1
               j = k - jj
               jk = ik + j
               mulk = -ap(jk)/ap(kk)
               t = conjg(mulk)
               call caxpy(j,t,ap(ik+1),1,ap(ij+1),1)
               ijj = ij + j
               ap(ijj) = cmplx(real(ap(ijj)),0.0e0)
               ap(jk) = mulk
               ij = ij - (j - 1)
  130       continue
c
c           set the pivot array.
c
            kpvt(k) = k
            if (swap) kpvt(k) = imax
         go to 190
  140    continue
c
c           2 x 2 pivot block.
c
            km1k = ik + k - 1
            ikm1 = ik - (k - 1)
            if (.not.swap) go to 160
c
c              perform an interchange.
c
               call cswap(imax,ap(im+1),1,ap(ikm1+1),1)
               imj = ikm1 + imax
               do 150 jj = imax, km1
                  j = km1 + imax - jj
                  jkm1 = ikm1 + j
                  t = conjg(ap(jkm1))
                  ap(jkm1) = conjg(ap(imj))
                  ap(imj) = t
                  imj = imj - (j - 1)
  150          continue
               t = ap(km1k)
               ap(km1k) = ap(imk)
               ap(imk) = t
  160       continue
c
c           perform the elimination.
c
            km2 = k - 2
            if (km2 .eq. 0) go to 180
               ak = ap(kk)/ap(km1k)
               km1km1 = ikm1 + k - 1
               akm1 = ap(km1km1)/conjg(ap(km1k))
               denom = 1.0e0 - ak*akm1
               ij = ik - (k - 1) - (k - 2)
               do 170 jj = 1, km2
                  j = km1 - jj
                  jk = ik + j
                  bk = ap(jk)/ap(km1k)
                  jkm1 = ikm1 + j
                  bkm1 = ap(jkm1)/conjg(ap(km1k))
                  mulk = (akm1*bk - bkm1)/denom
                  mulkm1 = (ak*bkm1 - bk)/denom
                  t = conjg(mulk)
                  call caxpy(j,t,ap(ik+1),1,ap(ij+1),1)
                  t = conjg(mulkm1)
                  call caxpy(j,t,ap(ikm1+1),1,ap(ij+1),1)
                  ap(jk) = mulk
                  ap(jkm1) = mulkm1
                  ijj = ij + j
                  ap(ijj) = cmplx(real(ap(ijj)),0.0e0)
                  ij = ij - (j - 1)
  170          continue
  180       continue
c
c           set the pivot array.
c
            kpvt(k) = 1 - k
            if (swap) kpvt(k) = -imax
            kpvt(k-1) = kpvt(k)
  190    continue
         ik = ik - (k - 1)
         if (kstep .eq. 2) ik = ik - (k - 2)
         k = k - kstep
      go to 10
  200 continue
      return
      end
      subroutine chpsl(ap,n,kpvt,b)

c*********************************************************************72

      integer n,kpvt(1)
      complex ap(1),b(1)
c
c     chisl solves the complex hermitian system
c     a * x = b
c     using the factors computed by chpfa.
c
c     on entry
c
c        ap      complex(n*(n+1)/2)
c                the output from chpfa.
c
c        n       integer
c                the order of the matrix  a .
c
c        kpvt    integer(n)
c                the pivot vector from chpfa.
c
c        b       complex(n)
c                the right hand side vector.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero may occur if  chpco  has set rcond .eq. 0.0
c        or  chpfa  has set info .ne. 0  .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call chpfa(ap,n,kpvt,info)
c           if (info .ne. 0) go to ...
c           do 10 j = 1, p
c              call chpsl(ap,n,kpvt,c(1,j))
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab.
c
c     subroutines and functions
c
c     blas caxpy,cdotc
c     fortran conjg,iabs
c
c     internal variables.
c
      complex ak,akm1,bk,bkm1,cdotc,denom,temp
      integer ik,ikm1,ikp1,k,kk,km1k,km1km1,kp
c
c     loop backward applying the transformations and
c     d inverse to b.
c
      k = n
      ik = (n*(n - 1))/2
   10 if (k .eq. 0) go to 80
         kk = ik + k
         if (kpvt(k) .lt. 0) go to 40
c
c           1 x 1 pivot block.
c
            if (k .eq. 1) go to 30
               kp = kpvt(k)
               if (kp .eq. k) go to 20
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
   20          continue
c
c              apply the transformation.
c
               call caxpy(k-1,b(k),ap(ik+1),1,b(1),1)
   30       continue
c
c           apply d inverse.
c
            b(k) = b(k)/ap(kk)
            k = k - 1
            ik = ik - k
         go to 70
   40    continue
c
c           2 x 2 pivot block.
c
            ikm1 = ik - (k - 1)
            if (k .eq. 2) go to 60
               kp = iabs(kpvt(k))
               if (kp .eq. k - 1) go to 50
c
c                 interchange.
c
                  temp = b(k-1)
                  b(k-1) = b(kp)
                  b(kp) = temp
   50          continue
c
c              apply the transformation.
c
               call caxpy(k-2,b(k),ap(ik+1),1,b(1),1)
               call caxpy(k-2,b(k-1),ap(ikm1+1),1,b(1),1)
   60       continue
c
c           apply d inverse.
c
            km1k = ik + k - 1
            kk = ik + k
            ak = ap(kk)/conjg(ap(km1k))
            km1km1 = ikm1 + k - 1
            akm1 = ap(km1km1)/ap(km1k)
            bk = b(k)/conjg(ap(km1k))
            bkm1 = b(k-1)/ap(km1k)
            denom = ak*akm1 - 1.0e0
            b(k) = (akm1*bk - bkm1)/denom
            b(k-1) = (ak*bkm1 - bk)/denom
            k = k - 2
            ik = ik - (k + 1) - k
   70    continue
      go to 10
   80 continue
c
c     loop forward applying the transformations.
c
      k = 1
      ik = 0
   90 if (k .gt. n) go to 160
         if (kpvt(k) .lt. 0) go to 120
c
c           1 x 1 pivot block.
c
            if (k .eq. 1) go to 110
c
c              apply the transformation.
c
               b(k) = b(k) + cdotc(k-1,ap(ik+1),1,b(1),1)
               kp = kpvt(k)
               if (kp .eq. k) go to 100
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
  100          continue
  110       continue
            ik = ik + k
            k = k + 1
         go to 150
  120    continue
c
c           2 x 2 pivot block.
c
            if (k .eq. 1) go to 140
c
c              apply the transformation.
c
               b(k) = b(k) + cdotc(k-1,ap(ik+1),1,b(1),1)
               ikp1 = ik + k
               b(k+1) = b(k+1) + cdotc(k-1,ap(ikp1+1),1,b(1),1)
               kp = iabs(kpvt(k))
               if (kp .eq. k) go to 130
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
  130          continue
  140       continue
            ik = ik + k + k + 1
            k = k + 2
  150    continue
      go to 90
  160 continue
      return
      end
      subroutine cpbco(abd,lda,n,m,rcond,z,info)

c*********************************************************************72

      integer lda,n,m,info
      complex abd(lda,1),z(1)
      real rcond
c
c     cpbco factors a complex hermitian positive definite matrix
c     stored in band form and estimates the condition of the matrix.
c
c     if  rcond  is not needed, cpbfa is slightly faster.
c     to solve  a*x = b , follow cpbco by cpbsl.
c     to compute  inverse(a)*c , follow cpbco by cpbsl.
c     to compute  determinant(a) , follow cpbco by cpbdi.
c
c     on entry
c
c        abd     complex(lda, n)
c                the matrix to be factored.  the columns of the upper
c                triangle are stored in the columns of abd and the
c                diagonals of the upper triangle are stored in the
c                rows of abd .  see the comments below for details.
c
c        lda     integer
c                the leading dimension of the array  abd .
c                lda must be .ge. m + 1 .
c
c        n       integer
c                the order of the matrix  a .
c
c        m       integer
c                the number of diagonals above the main diagonal.
c                0 .le. m .lt. n .
c
c     on return
c
c        abd     an upper triangular matrix  r , stored in band
c                form, so that  a = ctrans(r)*r .
c                if  info .ne. 0 , the factorization is not complete.
c
c        rcond   real
c                an estimate of the reciprocal condition of  a .
c                for the system  a*x = b , relative perturbations
c                in  a  and  b  of size  epsilon  may cause
c                relative perturbations in  x  of size  epsilon/rcond .
c                if  rcond  is so small that the logical expression
c                           1.0 + rcond .eq. 1.0
c                is true, then  a  may be singular to working
c                precision.  in particular,  rcond  is zero  if
c                exact singularity is detected or the estimate
c                underflows.  if info .ne. 0 , rcond is unchanged.
c
c        z       complex(n)
c                a work vector whose contents are usually unimportant.
c                if  a  is singular to working precision, then  z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .
c                if  info .ne. 0 , z  is unchanged.
c
c        info    integer
c                = 0  for normal return.
c                = k  signals an error condition.  the leading minor
c                     of order  k  is not positive definite.
c
c     band storage
c
c           if  a  is a hermitian positive definite band matrix,
c           the following program segment will set up the input.
c
c                   m = (band width above diagonal)
c                   do 20 j = 1, n
c                      i1 = max0(1, j-m)
c                      do 10 i = i1, j
c                         k = i-j+m+1
c                         abd(k,j) = a(i,j)
c                10    continue
c                20 continue
c
c           this uses  m + 1  rows of  a , except for the  m by m
c           upper left triangle, which is ignored.
c
c     example..  if the original matrix is
c
c           11 12 13  0  0  0
c           12 22 23 24  0  0
c           13 23 33 34 35  0
c            0 24 34 44 45 46
c            0  0 35 45 55 56
c            0  0  0 46 56 66
c
c     then  n = 6 , m = 2  and  abd  should contain
c
c            *  * 13 24 35 46
c            * 12 23 34 45 56
c           11 22 33 44 55 66
c
c     linpack.  this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     linpack cpbfa
c     blas caxpy,cdotc,csscal,scasum
c     fortran abs,aimag,amax1,cmplx,conjg,max0,min0,real
c
c     internal variables
c
      complex cdotc,ek,t,wk,wkm
      real anorm,s,scasum,sm,ynorm
      integer i,j,j2,k,kb,kp1,l,la,lb,lm,mu
c
      complex zdum,zdum2,csign1
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
      csign1(zdum,zdum2) = cabs1(zdum)*(zdum2/cabs1(zdum2))
c
c     find norm of a
c
      do 30 j = 1, n
         l = min0(j,m+1)
         mu = max0(m+2-j,1)
         z(j) = cmplx(scasum(l,abd(mu,j),1),0.0e0)
         k = j - l
         if (m .lt. mu) go to 20
         do 10 i = mu, m
            k = k + 1
            z(k) = cmplx(real(z(k))+cabs1(abd(i,j)),0.0e0)
   10    continue
   20    continue
   30 continue
      anorm = 0.0e0
      do 40 j = 1, n
         anorm = amax1(anorm,real(z(j)))
   40 continue
c
c     factor
c
      call cpbfa(abd,lda,n,m,info)
      if (info .ne. 0) go to 180
c
c        rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
c        estimate = norm(z)/norm(y) where  a*z = y  and  a*y = e .
c        the components of  e  are chosen to cause maximum local
c        growth in the elements of w  where  ctrans(r)*w = e .
c        the vectors are frequently rescaled to avoid overflow.
c
c        solve ctrans(r)*w = e
c
         ek = (1.0e0,0.0e0)
         do 50 j = 1, n
            z(j) = (0.0e0,0.0e0)
   50    continue
         do 110 k = 1, n
            if (cabs1(z(k)) .ne. 0.0e0) ek = csign1(ek,-z(k))
            if (cabs1(ek-z(k)) .le. real(abd(m+1,k))) go to 60
               s = real(abd(m+1,k))/cabs1(ek-z(k))
               call csscal(n,s,z,1)
               ek = cmplx(s,0.0e0)*ek
   60       continue
            wk = ek - z(k)
            wkm = -ek - z(k)
            s = cabs1(wk)
            sm = cabs1(wkm)
            wk = wk/abd(m+1,k)
            wkm = wkm/abd(m+1,k)
            kp1 = k + 1
            j2 = min0(k+m,n)
            i = m + 1
            if (kp1 .gt. j2) go to 100
               do 70 j = kp1, j2
                  i = i - 1
                  sm = sm + cabs1(z(j)+wkm*conjg(abd(i,j)))
                  z(j) = z(j) + wk*conjg(abd(i,j))
                  s = s + cabs1(z(j))
   70          continue
               if (s .ge. sm) go to 90
                  t = wkm - wk
                  wk = wkm
                  i = m + 1
                  do 80 j = kp1, j2
                     i = i - 1
                     z(j) = z(j) + t*conjg(abd(i,j))
   80             continue
   90          continue
  100       continue
            z(k) = wk
  110    continue
         s = 1.0e0/scasum(n,z,1)
         call csscal(n,s,z,1)
c
c        solve  r*y = w
c
         do 130 kb = 1, n
            k = n + 1 - kb
            if (cabs1(z(k)) .le. real(abd(m+1,k))) go to 120
               s = real(abd(m+1,k))/cabs1(z(k))
               call csscal(n,s,z,1)
  120       continue
            z(k) = z(k)/abd(m+1,k)
            lm = min0(k-1,m)
            la = m + 1 - lm
            lb = k - lm
            t = -z(k)
            call caxpy(lm,t,abd(la,k),1,z(lb),1)
  130    continue
         s = 1.0e0/scasum(n,z,1)
         call csscal(n,s,z,1)
c
         ynorm = 1.0e0
c
c        solve ctrans(r)*v = y
c
         do 150 k = 1, n
            lm = min0(k-1,m)
            la = m + 1 - lm
            lb = k - lm
            z(k) = z(k) - cdotc(lm,abd(la,k),1,z(lb),1)
            if (cabs1(z(k)) .le. real(abd(m+1,k))) go to 140
               s = real(abd(m+1,k))/cabs1(z(k))
               call csscal(n,s,z,1)
               ynorm = s*ynorm
  140       continue
            z(k) = z(k)/abd(m+1,k)
  150    continue
         s = 1.0e0/scasum(n,z,1)
         call csscal(n,s,z,1)
         ynorm = s*ynorm
c
c        solve  r*z = w
c
         do 170 kb = 1, n
            k = n + 1 - kb
            if (cabs1(z(k)) .le. real(abd(m+1,k))) go to 160
               s = real(abd(m+1,k))/cabs1(z(k))
               call csscal(n,s,z,1)
               ynorm = s*ynorm
  160       continue
            z(k) = z(k)/abd(m+1,k)
            lm = min0(k-1,m)
            la = m + 1 - lm
            lb = k - lm
            t = -z(k)
            call caxpy(lm,t,abd(la,k),1,z(lb),1)
  170    continue
c        make znorm = 1.0
         s = 1.0e0/scasum(n,z,1)
         call csscal(n,s,z,1)
         ynorm = s*ynorm
c
         if (anorm .ne. 0.0e0) rcond = ynorm/anorm
         if (anorm .eq. 0.0e0) rcond = 0.0e0
  180 continue
      return
      end
      subroutine cpbdi(abd,lda,n,m,det)

c*********************************************************************72

      integer lda,n,m
      complex abd(lda,1)
      real det(2)
c
c     cpbdi computes the determinant
c     of a complex hermitian positive definite band matrix
c     using the factors computed by cpbco or cpbfa.
c     if the inverse is needed, use cpbsl  n  times.
c
c     on entry
c
c        abd     complex(lda, n)
c                the output from cpbco or cpbfa.
c
c        lda     integer
c                the leading dimension of the array  abd .
c
c        n       integer
c                the order of the matrix  a .
c
c        m       integer
c                the number of diagonals above the main diagonal.
c
c     on return
c
c        det     real(2)
c                determinant of original matrix in the form
c                determinant = det(1) * 10.0**det(2)
c                with  1.0 .le. det(1) .lt. 10.0
c                or  det(1) .eq. 0.0 .
c
c     linpack.  this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     fortran real
c
c     internal variables
c
      real s
      integer i
c
c     compute determinant
c
      det(1) = 1.0e0
      det(2) = 0.0e0
      s = 10.0e0
      do 50 i = 1, n
         det(1) = real(abd(m+1,i))**2*det(1)
c     ...exit
         if (det(1) .eq. 0.0e0) go to 60
   10    if (det(1) .ge. 1.0e0) go to 20
            det(1) = s*det(1)
            det(2) = det(2) - 1.0e0
         go to 10
   20    continue
   30    if (det(1) .lt. s) go to 40
            det(1) = det(1)/s
            det(2) = det(2) + 1.0e0
         go to 30
   40    continue
   50 continue
   60 continue
      return
      end
      subroutine cpbfa(abd,lda,n,m,info)

c*********************************************************************72

      integer lda,n,m,info
      complex abd(lda,1)
c
c     cpbfa factors a complex hermitian positive definite matrix
c     stored in band form.
c
c     cpbfa is usually called by cpbco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c
c     on entry
c
c        abd     complex(lda, n)
c                the matrix to be factored.  the columns of the upper
c                triangle are stored in the columns of abd and the
c                diagonals of the upper triangle are stored in the
c                rows of abd .  see the comments below for details.
c
c        lda     integer
c                the leading dimension of the array  abd .
c                lda must be .ge. m + 1 .
c
c        n       integer
c                the order of the matrix  a .
c
c        m       integer
c                the number of diagonals above the main diagonal.
c                0 .le. m .lt. n .
c
c     on return
c
c        abd     an upper triangular matrix  r , stored in band
c                form, so that  a = ctrans(r)*r .
c
c        info    integer
c                = 0  for normal return.
c                = k  if the leading minor of order  k  is not
c                     positive definite.
c
c     band storage
c
c           if  a  is a hermitian positive definite band matrix,
c           the following program segment will set up the input.
c
c                   m = (band width above diagonal)
c                   do 20 j = 1, n
c                      i1 = max0(1, j-m)
c                      do 10 i = i1, j
c                         k = i-j+m+1
c                         abd(k,j) = a(i,j)
c                10    continue
c                20 continue
c
c     linpack.  this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas cdotc
c     fortran aimag,cmplx,conjg,max0,real,sqrt
c
c     internal variables
c
      complex cdotc,t
      real s
      integer ik,j,jk,k,mu
c     begin block with ...exits to 40
c
c
         do 30 j = 1, n
            info = j
            s = 0.0e0
            ik = m + 1
            jk = max0(j-m,1)
            mu = max0(m+2-j,1)
            if (m .lt. mu) go to 20
            do 10 k = mu, m
               t = abd(k,j) - cdotc(k-mu,abd(ik,jk),1,abd(mu,j),1)
               t = t/abd(m+1,jk)
               abd(k,j) = t
               s = s + real(t*conjg(t))
               ik = ik - 1
               jk = jk + 1
   10       continue
   20       continue
            s = real(abd(m+1,j)) - s
c     ......exit
            if (s .le. 0.0e0 .or. aimag(abd(m+1,j)) .ne. 0.0e0)
     *         go to 40
            abd(m+1,j) = cmplx(sqrt(s),0.0e0)
   30    continue
         info = 0
   40 continue
      return
      end
      subroutine cpbsl(abd,lda,n,m,b)

c*********************************************************************72

      integer lda,n,m
      complex abd(lda,1),b(1)
c
c     cpbsl solves the complex hermitian positive definite band
c     system  a*x = b
c     using the factors computed by cpbco or cpbfa.
c
c     on entry
c
c        abd     complex(lda, n)
c                the output from cpbco or cpbfa.
c
c        lda     integer
c                the leading dimension of the array  abd .
c
c        n       integer
c                the order of the matrix  a .
c
c        m       integer
c                the number of diagonals above the main diagonal.
c
c        b       complex(n)
c                the right hand side vector.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains
c        a zero on the diagonal.  technically this indicates
c        singularity but it is usually caused by improper subroutine
c        arguments.  it will not occur if the subroutines are called
c        correctly and  info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call cpbco(abd,lda,n,rcond,z,info)
c           if (rcond is too small .or. info .ne. 0) go to ...
c           do 10 j = 1, p
c              call cpbsl(abd,lda,n,c(1,j))
c        10 continue
c
c     linpack.  this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas caxpy,cdotc
c     fortran min0
c
c     internal variables
c
      complex cdotc,t
      integer k,kb,la,lb,lm
c
c     solve ctrans(r)*y = b
c
      do 10 k = 1, n
         lm = min0(k-1,m)
         la = m + 1 - lm
         lb = k - lm
         t = cdotc(lm,abd(la,k),1,b(lb),1)
         b(k) = (b(k) - t)/abd(m+1,k)
   10 continue
c
c     solve r*x = y
c
      do 20 kb = 1, n
         k = n + 1 - kb
         lm = min0(k-1,m)
         la = m + 1 - lm
         lb = k - lm
         b(k) = b(k)/abd(m+1,k)
         t = -b(k)
         call caxpy(lm,t,abd(la,k),1,b(lb),1)
   20 continue
      return
      end
      subroutine cpoco(a,lda,n,rcond,z,info)

c*********************************************************************72

      integer lda,n,info
      complex a(lda,1),z(1)
      real rcond
c
c     cpoco factors a complex hermitian positive definite matrix
c     and estimates the condition of the matrix.
c
c     if  rcond  is not needed, cpofa is slightly faster.
c     to solve  a*x = b , follow cpoco by cposl.
c     to compute  inverse(a)*c , follow cpoco by cposl.
c     to compute  determinant(a) , follow cpoco by cpodi.
c     to compute  inverse(a) , follow cpoco by cpodi.
c
c     on entry
c
c        a       complex(lda, n)
c                the hermitian matrix to be factored.  only the
c                diagonal and upper triangle are used.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix  r  so that  a =
c                ctrans(r)*r where  ctrans(r)  is the conjugate
c                transpose.  the strict lower triangle is unaltered.
c                if  info .ne. 0 , the factorization is not complete.
c
c        rcond   real
c                an estimate of the reciprocal condition of  a .
c                for the system  a*x = b , relative perturbations
c                in  a  and  b  of size  epsilon  may cause
c                relative perturbations in  x  of size  epsilon/rcond .
c                if  rcond  is so small that the logical expression
c                           1.0 + rcond .eq. 1.0
c                is true, then  a  may be singular to working
c                precision.  in particular,  rcond  is zero  if
c                exact singularity is detected or the estimate
c                underflows.  if info .ne. 0 , rcond is unchanged.
c
c        z       complex(n)
c                a work vector whose contents are usually unimportant.
c                if  a  is close to a singular matrix, then  z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .
c                if  info .ne. 0 , z  is unchanged.
c
c        info    integer
c                = 0  for normal return.
c                = k  signals an error condition.  the leading minor
c                     of order  k  is not positive definite.
c
c     linpack.  this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     linpack cpofa
c     blas caxpy,cdotc,csscal,scasum
c     fortran abs,aimag,amax1,cmplx,conjg,real
c
c     internal variables
c
      complex cdotc,ek,t,wk,wkm
      real anorm,s,scasum,sm,ynorm
      integer i,j,jm1,k,kb,kp1
c
      complex zdum,zdum2,csign1
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
      csign1(zdum,zdum2) = cabs1(zdum)*(zdum2/cabs1(zdum2))
c
c     find norm of a using only upper half
c
      do 30 j = 1, n
         z(j) = cmplx(scasum(j,a(1,j),1),0.0e0)
         jm1 = j - 1
         if (jm1 .lt. 1) go to 20
         do 10 i = 1, jm1
            z(i) = cmplx(real(z(i))+cabs1(a(i,j)),0.0e0)
   10    continue
   20    continue
   30 continue
      anorm = 0.0e0
      do 40 j = 1, n
         anorm = amax1(anorm,real(z(j)))
   40 continue
c
c     factor
c
      call cpofa(a,lda,n,info)
      if (info .ne. 0) go to 180
c
c        rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
c        estimate = norm(z)/norm(y) where  a*z = y  and  a*y = e .
c        the components of  e  are chosen to cause maximum local
c        growth in the elements of w  where  ctrans(r)*w = e .
c        the vectors are frequently rescaled to avoid overflow.
c
c        solve ctrans(r)*w = e
c
         ek = (1.0e0,0.0e0)
         do 50 j = 1, n
            z(j) = (0.0e0,0.0e0)
   50    continue
         do 110 k = 1, n
            if (cabs1(z(k)) .ne. 0.0e0) ek = csign1(ek,-z(k))
            if (cabs1(ek-z(k)) .le. real(a(k,k))) go to 60
               s = real(a(k,k))/cabs1(ek-z(k))
               call csscal(n,s,z,1)
               ek = cmplx(s,0.0e0)*ek
   60       continue
            wk = ek - z(k)
            wkm = -ek - z(k)
            s = cabs1(wk)
            sm = cabs1(wkm)
            wk = wk/a(k,k)
            wkm = wkm/a(k,k)
            kp1 = k + 1
            if (kp1 .gt. n) go to 100
               do 70 j = kp1, n
                  sm = sm + cabs1(z(j)+wkm*conjg(a(k,j)))
                  z(j) = z(j) + wk*conjg(a(k,j))
                  s = s + cabs1(z(j))
   70          continue
               if (s .ge. sm) go to 90
                  t = wkm - wk
                  wk = wkm
                  do 80 j = kp1, n
                     z(j) = z(j) + t*conjg(a(k,j))
   80             continue
   90          continue
  100       continue
            z(k) = wk
  110    continue
         s = 1.0e0/scasum(n,z,1)
         call csscal(n,s,z,1)
c
c        solve r*y = w
c
         do 130 kb = 1, n
            k = n + 1 - kb
            if (cabs1(z(k)) .le. real(a(k,k))) go to 120
               s = real(a(k,k))/cabs1(z(k))
               call csscal(n,s,z,1)
  120       continue
            z(k) = z(k)/a(k,k)
            t = -z(k)
            call caxpy(k-1,t,a(1,k),1,z(1),1)
  130    continue
         s = 1.0e0/scasum(n,z,1)
         call csscal(n,s,z,1)
c
         ynorm = 1.0e0
c
c        solve ctrans(r)*v = y
c
         do 150 k = 1, n
            z(k) = z(k) - cdotc(k-1,a(1,k),1,z(1),1)
            if (cabs1(z(k)) .le. real(a(k,k))) go to 140
               s = real(a(k,k))/cabs1(z(k))
               call csscal(n,s,z,1)
               ynorm = s*ynorm
  140       continue
            z(k) = z(k)/a(k,k)
  150    continue
         s = 1.0e0/scasum(n,z,1)
         call csscal(n,s,z,1)
         ynorm = s*ynorm
c
c        solve r*z = v
c
         do 170 kb = 1, n
            k = n + 1 - kb
            if (cabs1(z(k)) .le. real(a(k,k))) go to 160
               s = real(a(k,k))/cabs1(z(k))
               call csscal(n,s,z,1)
               ynorm = s*ynorm
  160       continue
            z(k) = z(k)/a(k,k)
            t = -z(k)
            call caxpy(k-1,t,a(1,k),1,z(1),1)
  170    continue
c        make znorm = 1.0
         s = 1.0e0/scasum(n,z,1)
         call csscal(n,s,z,1)
         ynorm = s*ynorm
c
         if (anorm .ne. 0.0e0) rcond = ynorm/anorm
         if (anorm .eq. 0.0e0) rcond = 0.0e0
  180 continue
      return
      end
      subroutine cpodi(a,lda,n,det,job)

c*********************************************************************72

      integer lda,n,job
      complex a(lda,1)
      real det(2)
c
c     cpodi computes the determinant and inverse of a certain
c     complex hermitian positive definite matrix (see below)
c     using the factors computed by cpoco, cpofa or cqrdc.
c
c     on entry
c
c        a       complex(lda, n)
c                the output  a  from cpoco or cpofa
c                or the output  x  from cqrdc.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        job     integer
c                = 11   both determinant and inverse.
c                = 01   inverse only.
c                = 10   determinant only.
c
c     on return
c
c        a       if cpoco or cpofa was used to factor  a  then
c                cpodi produces the upper half of inverse(a) .
c                if cqrdc was used to decompose  x  then
c                cpodi produces the upper half of inverse(ctrans(x)*x)
c                where ctrans(x) is the conjugate transpose.
c                elements of  a  below the diagonal are unchanged.
c                if the units digit of job is zero,  a  is unchanged.
c
c        det     real(2)
c                determinant of  a  or of  ctrans(x)*x  if requested.
c                otherwise not referenced.
c                determinant = det(1) * 10.0**det(2)
c                with  1.0 .le. det(1) .lt. 10.0
c                or  det(1) .eq. 0.0 .
c
c     error condition
c
c        a division by zero will occur if the input factor contains
c        a zero on the diagonal and the inverse is requested.
c        it will not occur if the subroutines are called correctly
c        and if cpoco or cpofa has set info .eq. 0 .
c
c     linpack.  this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas caxpy,cscal
c     fortran conjg,mod,real
c
c     internal variables
c
      complex t
      real s
      integer i,j,jm1,k,kp1
c
c     compute determinant
c
      if (job/10 .eq. 0) go to 70
         det(1) = 1.0e0
         det(2) = 0.0e0
         s = 10.0e0
         do 50 i = 1, n
            det(1) = real(a(i,i))**2*det(1)
c        ...exit
            if (det(1) .eq. 0.0e0) go to 60
   10       if (det(1) .ge. 1.0e0) go to 20
               det(1) = s*det(1)
               det(2) = det(2) - 1.0e0
            go to 10
   20       continue
   30       if (det(1) .lt. s) go to 40
               det(1) = det(1)/s
               det(2) = det(2) + 1.0e0
            go to 30
   40       continue
   50    continue
   60    continue
   70 continue
c
c     compute inverse(r)
c
      if (mod(job,10) .eq. 0) go to 140
         do 100 k = 1, n
            a(k,k) = (1.0e0,0.0e0)/a(k,k)
            t = -a(k,k)
            call cscal(k-1,t,a(1,k),1)
            kp1 = k + 1
            if (n .lt. kp1) go to 90
            do 80 j = kp1, n
               t = a(k,j)
               a(k,j) = (0.0e0,0.0e0)
               call caxpy(k,t,a(1,k),1,a(1,j),1)
   80       continue
   90       continue
  100    continue
c
c        form  inverse(r) * ctrans(inverse(r))
c
         do 130 j = 1, n
            jm1 = j - 1
            if (jm1 .lt. 1) go to 120
            do 110 k = 1, jm1
               t = conjg(a(k,j))
               call caxpy(k,t,a(1,j),1,a(1,k),1)
  110       continue
  120       continue
            t = conjg(a(j,j))
            call cscal(j,t,a(1,j),1)
  130    continue
  140 continue
      return
      end
      subroutine cpofa(a,lda,n,info)

c*********************************************************************72

      integer lda,n,info
      complex a(lda,1)
c
c     cpofa factors a complex hermitian positive definite matrix.
c
c     cpofa is usually called by cpoco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for cpoco) = (1 + 18/n)*(time for cpofa) .
c
c     on entry
c
c        a       complex(lda, n)
c                the hermitian matrix to be factored.  only the
c                diagonal and upper triangle are used.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix  r  so that  a =
c                ctrans(r)*r where  ctrans(r)  is the conjugate
c                transpose.  the strict lower triangle is unaltered.
c                if  info .ne. 0 , the factorization is not complete.
c
c        info    integer
c                = 0  for normal return.
c                = k  signals an error condition.  the leading minor
c                     of order  k  is not positive definite.
c
c     linpack.  this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas cdotc
c     fortran aimag,cmplx,conjg,real,sqrt
c
c     internal variables
c
      complex cdotc,t
      real s
      integer j,jm1,k
c     begin block with ...exits to 40
c
c
         do 30 j = 1, n
            info = j
            s = 0.0e0
            jm1 = j - 1
            if (jm1 .lt. 1) go to 20
            do 10 k = 1, jm1
               t = a(k,j) - cdotc(k-1,a(1,k),1,a(1,j),1)
               t = t/a(k,k)
               a(k,j) = t
               s = s + real(t*conjg(t))
   10       continue
   20       continue
            s = real(a(j,j)) - s
c     ......exit
            if (s .le. 0.0e0 .or. aimag(a(j,j)) .ne. 0.0e0) go to 40
            a(j,j) = cmplx(sqrt(s),0.0e0)
   30    continue
         info = 0
   40 continue
      return
      end
      subroutine cposl(a,lda,n,b)

c*********************************************************************72

      integer lda,n
      complex a(lda,1),b(1)
c
c     cposl solves the complex hermitian positive definite system
c     a * x = b
c     using the factors computed by cpoco or cpofa.
c
c     on entry
c
c        a       complex(lda, n)
c                the output from cpoco or cpofa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        b       complex(n)
c                the right hand side vector.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains
c        a zero on the diagonal.  technically this indicates
c        singularity but it is usually caused by improper subroutine
c        arguments.  it will not occur if the subroutines are called
c        correctly and  info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call cpoco(a,lda,n,rcond,z,info)
c           if (rcond is too small .or. info .ne. 0) go to ...
c           do 10 j = 1, p
c              call cposl(a,lda,n,c(1,j))
c        10 continue
c
c     linpack.  this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas caxpy,cdotc
c
c     internal variables
c
      complex cdotc,t
      integer k,kb
c
c     solve ctrans(r)*y = b
c
      do 10 k = 1, n
         t = cdotc(k-1,a(1,k),1,b(1),1)
         b(k) = (b(k) - t)/a(k,k)
   10 continue
c
c     solve r*x = y
c
      do 20 kb = 1, n
         k = n + 1 - kb
         b(k) = b(k)/a(k,k)
         t = -b(k)
         call caxpy(k-1,t,a(1,k),1,b(1),1)
   20 continue
      return
      end
      subroutine cppco(ap,n,rcond,z,info)

c*********************************************************************72

      integer n,info
      complex ap(1),z(1)
      real rcond
c
c     cppco factors a complex hermitian positive definite matrix
c     stored in packed form
c     and estimates the condition of the matrix.
c
c     if  rcond  is not needed, cppfa is slightly faster.
c     to solve  a*x = b , follow cppco by cppsl.
c     to compute  inverse(a)*c , follow cppco by cppsl.
c     to compute  determinant(a) , follow cppco by cppdi.
c     to compute  inverse(a) , follow cppco by cppdi.
c
c     on entry
c
c        ap      complex (n*(n+1)/2)
c                the packed form of a hermitian matrix  a .  the
c                columns of the upper triangle are stored sequentially
c                in a one-dimensional array of length  n*(n+1)/2 .
c                see comments below for details.
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        ap      an upper triangular matrix  r , stored in packed
c                form, so that  a = ctrans(r)*r .
c                if  info .ne. 0 , the factorization is not complete.
c
c        rcond   real
c                an estimate of the reciprocal condition of  a .
c                for the system  a*x = b , relative perturbations
c                in  a  and  b  of size  epsilon  may cause
c                relative perturbations in  x  of size  epsilon/rcond .
c                if  rcond  is so small that the logical expression
c                           1.0 + rcond .eq. 1.0
c                is true, then  a  may be singular to working
c                precision.  in particular,  rcond  is zero  if
c                exact singularity is detected or the estimate
c                underflows.  if info .ne. 0 , rcond is unchanged.
c
c        z       complex(n)
c                a work vector whose contents are usually unimportant.
c                if  a  is singular to working precision, then  z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .
c                if  info .ne. 0 , z  is unchanged.
c
c        info    integer
c                = 0  for normal return.
c                = k  signals an error condition.  the leading minor
c                     of order  k  is not positive definite.
c
c     packed storage
c
c          the following program segment will pack the upper
c          triangle of a hermitian matrix.
c
c                k = 0
c                do 20 j = 1, n
c                   do 10 i = 1, j
c                      k = k + 1
c                      ap(k) = a(i,j)
c             10    continue
c             20 continue
c
c     linpack.  this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     linpack cppfa
c     blas caxpy,cdotc,csscal,scasum
c     fortran abs,aimag,amax1,cmplx,conjg,real
c
c     internal variables
c
      complex cdotc,ek,t,wk,wkm
      real anorm,s,scasum,sm,ynorm
      integer i,ij,j,jm1,j1,k,kb,kj,kk,kp1
c
      complex zdum,zdum2,csign1
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
      csign1(zdum,zdum2) = cabs1(zdum)*(zdum2/cabs1(zdum2))
c
c     find norm of a
c
      j1 = 1
      do 30 j = 1, n
         z(j) = cmplx(scasum(j,ap(j1),1),0.0e0)
         ij = j1
         j1 = j1 + j
         jm1 = j - 1
         if (jm1 .lt. 1) go to 20
         do 10 i = 1, jm1
            z(i) = cmplx(real(z(i))+cabs1(ap(ij)),0.0e0)
            ij = ij + 1
   10    continue
   20    continue
   30 continue
      anorm = 0.0e0
      do 40 j = 1, n
         anorm = amax1(anorm,real(z(j)))
   40 continue
c
c     factor
c
      call cppfa(ap,n,info)
      if (info .ne. 0) go to 180
c
c        rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
c        estimate = norm(z)/norm(y) where  a*z = y  and  a*y = e .
c        the components of  e  are chosen to cause maximum local
c        growth in the elements of w  where  ctrans(r)*w = e .
c        the vectors are frequently rescaled to avoid overflow.
c
c        solve ctrans(r)*w = e
c
         ek = (1.0e0,0.0e0)
         do 50 j = 1, n
            z(j) = (0.0e0,0.0e0)
   50    continue
         kk = 0
         do 110 k = 1, n
            kk = kk + k
            if (cabs1(z(k)) .ne. 0.0e0) ek = csign1(ek,-z(k))
            if (cabs1(ek-z(k)) .le. real(ap(kk))) go to 60
               s = real(ap(kk))/cabs1(ek-z(k))
               call csscal(n,s,z,1)
               ek = cmplx(s,0.0e0)*ek
   60       continue
            wk = ek - z(k)
            wkm = -ek - z(k)
            s = cabs1(wk)
            sm = cabs1(wkm)
            wk = wk/ap(kk)
            wkm = wkm/ap(kk)
            kp1 = k + 1
            kj = kk + k
            if (kp1 .gt. n) go to 100
               do 70 j = kp1, n
                  sm = sm + cabs1(z(j)+wkm*conjg(ap(kj)))
                  z(j) = z(j) + wk*conjg(ap(kj))
                  s = s + cabs1(z(j))
                  kj = kj + j
   70          continue
               if (s .ge. sm) go to 90
                  t = wkm - wk
                  wk = wkm
                  kj = kk + k
                  do 80 j = kp1, n
                     z(j) = z(j) + t*conjg(ap(kj))
                     kj = kj + j
   80             continue
   90          continue
  100       continue
            z(k) = wk
  110    continue
         s = 1.0e0/scasum(n,z,1)
         call csscal(n,s,z,1)
c
c        solve r*y = w
c
         do 130 kb = 1, n
            k = n + 1 - kb
            if (cabs1(z(k)) .le. real(ap(kk))) go to 120
               s = real(ap(kk))/cabs1(z(k))
               call csscal(n,s,z,1)
  120       continue
            z(k) = z(k)/ap(kk)
            kk = kk - k
            t = -z(k)
            call caxpy(k-1,t,ap(kk+1),1,z(1),1)
  130    continue
         s = 1.0e0/scasum(n,z,1)
         call csscal(n,s,z,1)
c
         ynorm = 1.0e0
c
c        solve ctrans(r)*v = y
c
         do 150 k = 1, n
            z(k) = z(k) - cdotc(k-1,ap(kk+1),1,z(1),1)
            kk = kk + k
            if (cabs1(z(k)) .le. real(ap(kk))) go to 140
               s = real(ap(kk))/cabs1(z(k))
               call csscal(n,s,z,1)
               ynorm = s*ynorm
  140       continue
            z(k) = z(k)/ap(kk)
  150    continue
         s = 1.0e0/scasum(n,z,1)
         call csscal(n,s,z,1)
         ynorm = s*ynorm
c
c        solve r*z = v
c
         do 170 kb = 1, n
            k = n + 1 - kb
            if (cabs1(z(k)) .le. real(ap(kk))) go to 160
               s = real(ap(kk))/cabs1(z(k))
               call csscal(n,s,z,1)
               ynorm = s*ynorm
  160       continue
            z(k) = z(k)/ap(kk)
            kk = kk - k
            t = -z(k)
            call caxpy(k-1,t,ap(kk+1),1,z(1),1)
  170    continue
c        make znorm = 1.0
         s = 1.0e0/scasum(n,z,1)
         call csscal(n,s,z,1)
         ynorm = s*ynorm
c
         if (anorm .ne. 0.0e0) rcond = ynorm/anorm
         if (anorm .eq. 0.0e0) rcond = 0.0e0
  180 continue
      return
      end
      subroutine cppdi(ap,n,det,job)

c*********************************************************************72

      integer n,job
      complex ap(1)
      real det(2)
c
c     cppdi computes the determinant and inverse
c     of a complex hermitian positive definite matrix
c     using the factors computed by cppco or cppfa .
c
c     on entry
c
c        ap      complex (n*(n+1)/2)
c                the output from cppco or cppfa.
c
c        n       integer
c                the order of the matrix  a .
c
c        job     integer
c                = 11   both determinant and inverse.
c                = 01   inverse only.
c                = 10   determinant only.
c
c     on return
c
c        ap      the upper triangular half of the inverse .
c                the strict lower triangle is unaltered.
c
c        det     real(2)
c                determinant of original matrix if requested.
c                otherwise not referenced.
c                determinant = det(1) * 10.0**det(2)
c                with  1.0 .le. det(1) .lt. 10.0
c                or  det(1) .eq. 0.0 .
c
c     error condition
c
c        a division by zero will occur if the input factor contains
c        a zero on the diagonal and the inverse is requested.
c        it will not occur if the subroutines are called correctly
c        and if cpoco or cpofa has set info .eq. 0 .
c
c     linpack.  this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas caxpy,cscal
c     fortran conjg,mod,real
c
c     internal variables
c
      complex t
      real s
      integer i,ii,j,jj,jm1,j1,k,kj,kk,kp1,k1
c
c     compute determinant
c
      if (job/10 .eq. 0) go to 70
         det(1) = 1.0e0
         det(2) = 0.0e0
         s = 10.0e0
         ii = 0
         do 50 i = 1, n
            ii = ii + i
            det(1) = real(ap(ii))**2*det(1)
c        ...exit
            if (det(1) .eq. 0.0e0) go to 60
   10       if (det(1) .ge. 1.0e0) go to 20
               det(1) = s*det(1)
               det(2) = det(2) - 1.0e0
            go to 10
   20       continue
   30       if (det(1) .lt. s) go to 40
               det(1) = det(1)/s
               det(2) = det(2) + 1.0e0
            go to 30
   40       continue
   50    continue
   60    continue
   70 continue
c
c     compute inverse(r)
c
      if (mod(job,10) .eq. 0) go to 140
         kk = 0
         do 100 k = 1, n
            k1 = kk + 1
            kk = kk + k
            ap(kk) = (1.0e0,0.0e0)/ap(kk)
            t = -ap(kk)
            call cscal(k-1,t,ap(k1),1)
            kp1 = k + 1
            j1 = kk + 1
            kj = kk + k
            if (n .lt. kp1) go to 90
            do 80 j = kp1, n
               t = ap(kj)
               ap(kj) = (0.0e0,0.0e0)
               call caxpy(k,t,ap(k1),1,ap(j1),1)
               j1 = j1 + j
               kj = kj + j
   80       continue
   90       continue
  100    continue
c
c        form  inverse(r) * ctrans(inverse(r))
c
         jj = 0
         do 130 j = 1, n
            j1 = jj + 1
            jj = jj + j
            jm1 = j - 1
            k1 = 1
            kj = j1
            if (jm1 .lt. 1) go to 120
            do 110 k = 1, jm1
               t = conjg(ap(kj))
               call caxpy(k,t,ap(j1),1,ap(k1),1)
               k1 = k1 + k
               kj = kj + 1
  110       continue
  120       continue
            t = conjg(ap(jj))
            call cscal(j,t,ap(j1),1)
  130    continue
  140 continue
      return
      end
      subroutine cppfa(ap,n,info)

c*********************************************************************72

      integer n,info
      complex ap(1)
c
c     cppfa factors a complex hermitian positive definite matrix
c     stored in packed form.
c
c     cppfa is usually called by cppco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for cppco) = (1 + 18/n)*(time for cppfa) .
c
c     on entry
c
c        ap      complex (n*(n+1)/2)
c                the packed form of a hermitian matrix  a .  the
c                columns of the upper triangle are stored sequentially
c                in a one-dimensional array of length  n*(n+1)/2 .
c                see comments below for details.
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        ap      an upper triangular matrix  r , stored in packed
c                form, so that  a = ctrans(r)*r .
c
c        info    integer
c                = 0  for normal return.
c                = k  if the leading minor of order  k  is not
c                     positive definite.
c
c
c     packed storage
c
c          the following program segment will pack the upper
c          triangle of a hermitian matrix.
c
c                k = 0
c                do 20 j = 1, n
c                   do 10 i = 1, j
c                      k = k + 1
c                      ap(k) = a(i,j)
c             10    continue
c             20 continue
c
c     linpack.  this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas cdotc
c     fortran aimag,cmplx,conjg,real,sqrt
c
c     internal variables
c
      complex cdotc,t
      real s
      integer j,jj,jm1,k,kj,kk
c     begin block with ...exits to 40
c
c
         jj = 0
         do 30 j = 1, n
            info = j
            s = 0.0e0
            jm1 = j - 1
            kj = jj
            kk = 0
            if (jm1 .lt. 1) go to 20
            do 10 k = 1, jm1
               kj = kj + 1
               t = ap(kj) - cdotc(k-1,ap(kk+1),1,ap(jj+1),1)
               kk = kk + k
               t = t/ap(kk)
               ap(kj) = t
               s = s + real(t*conjg(t))
   10       continue
   20       continue
            jj = jj + j
            s = real(ap(jj)) - s
c     ......exit
            if (s .le. 0.0e0 .or. aimag(ap(jj)) .ne. 0.0e0) go to 40
            ap(jj) = cmplx(sqrt(s),0.0e0)
   30    continue
         info = 0
   40 continue
      return
      end
      subroutine cppsl(ap,n,b)

c*********************************************************************72

      integer n
      complex ap(1),b(1)
c
c     cppsl solves the complex hermitian positive definite system
c     a * x = b
c     using the factors computed by cppco or cppfa.
c
c     on entry
c
c        ap      complex (n*(n+1)/2)
c                the output from cppco or cppfa.
c
c        n       integer
c                the order of the matrix  a .
c
c        b       complex(n)
c                the right hand side vector.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains
c        a zero on the diagonal.  technically this indicates
c        singularity but it is usually caused by improper subroutine
c        arguments.  it will not occur if the subroutines are called
c        correctly and  info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call cppco(ap,n,rcond,z,info)
c           if (rcond is too small .or. info .ne. 0) go to ...
c           do 10 j = 1, p
c              call cppsl(ap,n,c(1,j))
c        10 continue
c
c     linpack.  this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas caxpy,cdotc
c
c     internal variables
c
      complex cdotc,t
      integer k,kb,kk
c
      kk = 0
      do 10 k = 1, n
         t = cdotc(k-1,ap(kk+1),1,b(1),1)
         kk = kk + k
         b(k) = (b(k) - t)/ap(kk)
   10 continue
      do 20 kb = 1, n
         k = n + 1 - kb
         b(k) = b(k)/ap(kk)
         kk = kk - k
         t = -b(k)
         call caxpy(k-1,t,ap(kk+1),1,b(1),1)
   20 continue
      return
      end
      subroutine cptsl(n,d,e,b)

c*********************************************************************72

      integer n
      complex d(1),e(1),b(1)
c
c     cptsl given a positive definite tridiagonal matrix and a right
c     hand side will find the solution.
c
c     on entry
c
c        n        integer
c                 is the order of the tridiagonal matrix.
c
c        d        complex(n)
c                 is the diagonal of the tridiagonal matrix.
c                 on output d is destroyed.
c
c        e        complex(n)
c                 is the offdiagonal of the tridiagonal matrix.
c                 e(1) through e(n-1) should contain the
c                 offdiagonal.
c
c        b        complex(n)
c                 is the right hand side vector.
c
c     on return
c
c        b        contains the soultion.
c
c     linpack. this version dated 08/14/78 .
c     jack dongarra, argonne national laboratory.
c
c     no externals
c     fortran conjg,mod
c
c     internal variables
c
      integer k,kbm1,ke,kf,kp1,nm1,nm1d2
      complex t1,t2
c
c     check for 1 x 1 case
c
      if (n .ne. 1) go to 10
         b(1) = b(1)/d(1)
      go to 70
   10 continue
         nm1 = n - 1
         nm1d2 = nm1/2
         if (n .eq. 2) go to 30
            kbm1 = n - 1
c
c           zero top half of subdiagonal and bottom half of
c           superdiagonal
c
            do 20 k = 1, nm1d2
               t1 = conjg(e(k))/d(k)
               d(k+1) = d(k+1) - t1*e(k)
               b(k+1) = b(k+1) - t1*b(k)
               t2 = e(kbm1)/d(kbm1+1)
               d(kbm1) = d(kbm1) - t2*conjg(e(kbm1))
               b(kbm1) = b(kbm1) - t2*b(kbm1+1)
               kbm1 = kbm1 - 1
   20       continue
   30    continue
         kp1 = nm1d2 + 1
c
c        clean up for possible 2 x 2 block at center
c
         if (mod(n,2) .ne. 0) go to 40
            t1 = conjg(e(kp1))/d(kp1)
            d(kp1+1) = d(kp1+1) - t1*e(kp1)
            b(kp1+1) = b(kp1+1) - t1*b(kp1)
            kp1 = kp1 + 1
   40    continue
c
c        back solve starting at the center, going towards the top
c        and bottom
c
         b(kp1) = b(kp1)/d(kp1)
         if (n .eq. 2) go to 60
            k = kp1 - 1
            ke = kp1 + nm1d2 - 1
            do 50 kf = kp1, ke
               b(k) = (b(k) - e(k)*b(k+1))/d(k)
               b(kf+1) = (b(kf+1) - conjg(e(kf))*b(kf))/d(kf+1)
               k = k - 1
   50       continue
   60    continue
         if (mod(n,2) .eq. 0) b(1) = (b(1) - e(1)*b(2))/d(1)
   70 continue
      return
      end
      subroutine cqrdc(x,ldx,n,p,qraux,jpvt,work,job)

c*********************************************************************72

      integer ldx,n,p,job
      integer jpvt(1)
      complex x(ldx,1),qraux(1),work(1)
c
c     cqrdc uses householder transformations to compute the qr
c     factorization of an n by p matrix x.  column pivoting
c     based on the 2-norms of the reduced columns may be
c     performed at the users option.
c
c     on entry
c
c        x       complex(ldx,p), where ldx .ge. n.
c                x contains the matrix whose decomposition is to be
c                computed.
c
c        ldx     integer.
c                ldx is the leading dimension of the array x.
c
c        n       integer.
c                n is the number of rows of the matrix x.
c
c        p       integer.
c                p is the number of columns of the matrix x.
c
c        jpvt    integer(p).
c                jpvt contains integers that control the selection
c                of the pivot columns.  the k-th column x(k) of x
c                is placed in one of three classes according to the
c                value of jpvt(k).
c
c                   if jpvt(k) .gt. 0, then x(k) is an initial
c                                      column.
c
c                   if jpvt(k) .eq. 0, then x(k) is a free column.
c
c                   if jpvt(k) .lt. 0, then x(k) is a final column.
c
c                before the decomposition is computed, initial columns
c                are moved to the beginning of the array x and final
c                columns to the end.  both initial and final columns
c                are frozen in place during the computation and only
c                free columns are moved.  at the k-th stage of the
c                reduction, if x(k) is occupied by a free column
c                it is interchanged with the free column of largest
c                reduced norm.  jpvt is not referenced if
c                job .eq. 0.
c
c        work    complex(p).
c                work is a work array.  work is not referenced if
c                job .eq. 0.
c
c        job     integer.
c                job is an integer that initiates column pivoting.
c                if job .eq. 0, no pivoting is done.
c                if job .ne. 0, pivoting is done.
c
c     on return
c
c        x       x contains in its upper triangle the upper
c                triangular matrix r of the qr factorization.
c                below its diagonal x contains information from
c                which the unitary part of the decomposition
c                can be recovered.  note that if pivoting has
c                been requested, the decomposition is not that
c                of the original matrix x but that of x
c                with its columns permuted as described by jpvt.
c
c        qraux   complex(p).
c                qraux contains further information required to recover
c                the unitary part of the decomposition.
c
c        jpvt    jpvt(k) contains the index of the column of the
c                original matrix that has been interchanged into
c                the k-th column, if pivoting was requested.
c
c     linpack. this version dated 08/14/78 .
c     g.w. stewart, university of maryland, argonne national lab.
c
c     cqrdc uses the following functions and subprograms.
c
c     blas caxpy,cdotc,cscal,cswap,scnrm2
c     fortran abs,aimag,amax1,cabs,cmplx,csqrt,min0,real
c
c     internal variables
c
      integer j,jp,l,lp1,lup,maxj,pl,pu
      real maxnrm,scnrm2,tt
      complex cdotc,nrmxl,t
      logical negj,swapj
c
      complex csign,zdum,zdum1,zdum2
      real cabs1
      csign(zdum1,zdum2) = cabs(zdum1)*(zdum2/cabs(zdum2))
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
c
      pl = 1
      pu = 0
      if (job .eq. 0) go to 60
c
c        pivoting has been requested.  rearrange the columns
c        according to jpvt.
c
         do 20 j = 1, p
            swapj = jpvt(j) .gt. 0
            negj = jpvt(j) .lt. 0
            jpvt(j) = j
            if (negj) jpvt(j) = -j
            if (.not.swapj) go to 10
               if (j .ne. pl) call cswap(n,x(1,pl),1,x(1,j),1)
               jpvt(j) = jpvt(pl)
               jpvt(pl) = j
               pl = pl + 1
   10       continue
   20    continue
         pu = p
         do 50 jj = 1, p
            j = p - jj + 1
            if (jpvt(j) .ge. 0) go to 40
               jpvt(j) = -jpvt(j)
               if (j .eq. pu) go to 30
                  call cswap(n,x(1,pu),1,x(1,j),1)
                  jp = jpvt(pu)
                  jpvt(pu) = jpvt(j)
                  jpvt(j) = jp
   30          continue
               pu = pu - 1
   40       continue
   50    continue
   60 continue
c
c     compute the norms of the free columns.
c
      if (pu .lt. pl) go to 80
      do 70 j = pl, pu
         qraux(j) = cmplx(scnrm2(n,x(1,j),1),0.0e0)
         work(j) = qraux(j)
   70 continue
   80 continue
c
c     perform the householder reduction of x.
c
      lup = min0(n,p)
      do 200 l = 1, lup
         if (l .lt. pl .or. l .ge. pu) go to 120
c
c           locate the column of largest norm and bring it
c           into the pivot position.
c
            maxnrm = 0.0e0
            maxj = l
            do 100 j = l, pu
               if (real(qraux(j)) .le. maxnrm) go to 90
                  maxnrm = real(qraux(j))
                  maxj = j
   90          continue
  100       continue
            if (maxj .eq. l) go to 110
               call cswap(n,x(1,l),1,x(1,maxj),1)
               qraux(maxj) = qraux(l)
               work(maxj) = work(l)
               jp = jpvt(maxj)
               jpvt(maxj) = jpvt(l)
               jpvt(l) = jp
  110       continue
  120    continue
         qraux(l) = (0.0e0,0.0e0)
         if (l .eq. n) go to 190
c
c           compute the householder transformation for column l.
c
            nrmxl = cmplx(scnrm2(n-l+1,x(l,l),1),0.0e0)
            if (cabs1(nrmxl) .eq. 0.0e0) go to 180
               if (cabs1(x(l,l)) .ne. 0.0e0)
     *            nrmxl = csign(nrmxl,x(l,l))
               call cscal(n-l+1,(1.0e0,0.0e0)/nrmxl,x(l,l),1)
               x(l,l) = (1.0e0,0.0e0) + x(l,l)
c
c              apply the transformation to the remaining columns,
c              updating the norms.
c
               lp1 = l + 1
               if (p .lt. lp1) go to 170
               do 160 j = lp1, p
                  t = -cdotc(n-l+1,x(l,l),1,x(l,j),1)/x(l,l)
                  call caxpy(n-l+1,t,x(l,l),1,x(l,j),1)
                  if (j .lt. pl .or. j .gt. pu) go to 150
                  if (cabs1(qraux(j)) .eq. 0.0e0) go to 150
                     tt = 1.0e0 - (cabs(x(l,j))/real(qraux(j)))**2
                     tt = amax1(tt,0.0e0)
                     t = cmplx(tt,0.0e0)
                     tt = 1.0e0
     *                    + 0.05e0*tt*(real(qraux(j))/real(work(j)))**2
                     if (tt .eq. 1.0e0) go to 130
                        qraux(j) = qraux(j)*csqrt(t)
                     go to 140
  130                continue
                        qraux(j) = cmplx(scnrm2(n-l,x(l+1,j),1),0.0e0)
                        work(j) = qraux(j)
  140                continue
  150             continue
  160          continue
  170          continue
c
c              save the transformation.
c
               qraux(l) = x(l,l)
               x(l,l) = -nrmxl
  180       continue
  190    continue
  200 continue
      return
      end
      subroutine cqrsl(x,ldx,n,k,qraux,y,qy,qty,b,rsd,xb,job,info)

c*********************************************************************72

      integer ldx,n,k,job,info
      complex x(ldx,1),qraux(1),y(1),qy(1),qty(1),b(1),rsd(1),xb(1)
c
c     cqrsl applies the output of cqrdc to compute coordinate
c     transformations, projections, and least squares solutions.
c     for k .le. min(n,p), let xk be the matrix
c
c            xk = (x(jpvt(1)),x(jpvt(2)), ... ,x(jpvt(k)))
c
c     formed from columnns jpvt(1), ... ,jpvt(k) of the original
c     n x p matrix x that was input to cqrdc (if no pivoting was
c     done, xk consists of the first k columns of x in their
c     original order).  cqrdc produces a factored unitary matrix q
c     and an upper triangular matrix r such that
c
c              xk = q * (r)
c                       (0)
c
c     this information is contained in coded form in the arrays
c     x and qraux.
c
c     on entry
c
c        x      complex(ldx,p).
c               x contains the output of cqrdc.
c
c        ldx    integer.
c               ldx is the leading dimension of the array x.
c
c        n      integer.
c               n is the number of rows of the matrix xk.  it must
c               have the same value as n in cqrdc.
c
c        k      integer.
c               k is the number of columns of the matrix xk.  k
c               must nnot be greater than min(n,p), where p is the
c               same as in the calling sequence to cqrdc.
c
c        qraux  complex(p).
c               qraux contains the auxiliary output from cqrdc.
c
c        y      complex(n)
c               y contains an n-vector that is to be manipulated
c               by cqrsl.
c
c        job    integer.
c               job specifies what is to be computed.  job has
c               the decimal expansion abcde, with the following
c               meaning.
c
c                    if a.ne.0, compute qy.
c                    if b,c,d, or e .ne. 0, compute qty.
c                    if c.ne.0, compute b.
c                    if d.ne.0, compute rsd.
c                    if e.ne.0, compute xb.
c
c               note that a request to compute b, rsd, or xb
c               automatically triggers the computation of qty, for
c               which an array must be provided in the calling
c               sequence.
c
c     on return
c
c        qy     complex(n).
c               qy conntains q*y, if its computation has been
c               requested.
c
c        qty    complex(n).
c               qty contains ctrans(q)*y, if its computation has
c               been requested.  here ctrans(q) is the conjugate
c               transpose of the matrix q.
c
c        b      complex(k)
c               b contains the solution of the least squares problem
c
c                    minimize norm2(y - xk*b),
c
c               if its computation has been requested.  (note that
c               if pivoting was requested in cqrdc, the j-th
c               component of b will be associated with column jpvt(j)
c               of the original matrix x that was input into cqrdc.)
c
c        rsd    complex(n).
c               rsd contains the least squares residual y - xk*b,
c               if its computation has been requested.  rsd is
c               also the orthogonal projection of y onto the
c               orthogonal complement of the column space of xk.
c
c        xb     complex(n).
c               xb contains the least squares approximation xk*b,
c               if its computation has been requested.  xb is also
c               the orthogonal projection of y onto the column space
c               of x.
c
c        info   integer.
c               info is zero unless the computation of b has
c               been requested and r is exactly singular.  in
c               this case, info is the index of the first zero
c               diagonal element of r and b is left unaltered.
c
c     the parameters qy, qty, b, rsd, and xb are not referenced
c     if their computation is not requested and in this case
c     can be replaced by dummy variables in the calling program.
c     to save storage, the user may in some cases use the same
c     array for different parameters in the calling sequence.  a
c     frequently occuring example is when one wishes to compute
c     any of b, rsd, or xb and does not need y or qty.  in this
c     case one may identify y, qty, and one of b, rsd, or xb, while
c     providing separate arrays for anything else that is to be
c     computed.  thus the calling sequence
c
c          call cqrsl(x,ldx,n,k,qraux,y,dum,y,b,y,dum,110,info)
c
c     will result in the computation of b and rsd, with rsd
c     overwriting y.  more generally, each item in the following
c     list contains groups of permissible identifications for
c     a single callinng sequence.
c
c          1. (y,qty,b) (rsd) (xb) (qy)
c
c          2. (y,qty,rsd) (b) (xb) (qy)
c
c          3. (y,qty,xb) (b) (rsd) (qy)
c
c          4. (y,qy) (qty,b) (rsd) (xb)
c
c          5. (y,qy) (qty,rsd) (b) (xb)
c
c          6. (y,qy) (qty,xb) (b) (rsd)
c
c     in any group the value returned in the array allocated to
c     the group corresponds to the last member of the group.
c
c     linpack. this version dated 08/14/78 .
c     g.w. stewart, university of maryland, argonne national lab.
c
c     cqrsl uses the following functions and subprograms.
c
c     blas caxpy,ccopy,cdotc
c     fortran abs,aimag,min0,mod,real
c
c     internal variables
c
      integer i,j,jj,ju,kp1
      complex cdotc,t,temp
      logical cb,cqy,cqty,cr,cxb
c
      complex zdum
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
c
c     set info flag.
c
      info = 0
c
c     determine what is to be computed.
c
      cqy = job/10000 .ne. 0
      cqty = mod(job,10000) .ne. 0
      cb = mod(job,1000)/100 .ne. 0
      cr = mod(job,100)/10 .ne. 0
      cxb = mod(job,10) .ne. 0
      ju = min0(k,n-1)
c
c     special action when n=1.
c
      if (ju .ne. 0) go to 40
         if (cqy) qy(1) = y(1)
         if (cqty) qty(1) = y(1)
         if (cxb) xb(1) = y(1)
         if (.not.cb) go to 30
            if (cabs1(x(1,1)) .ne. 0.0e0) go to 10
               info = 1
            go to 20
   10       continue
               b(1) = y(1)/x(1,1)
   20       continue
   30    continue
         if (cr) rsd(1) = (0.0e0,0.0e0)
      go to 250
   40 continue
c
c        set up to compute qy or qty.
c
         if (cqy) call ccopy(n,y,1,qy,1)
         if (cqty) call ccopy(n,y,1,qty,1)
         if (.not.cqy) go to 70
c
c           compute qy.
c
            do 60 jj = 1, ju
               j = ju - jj + 1
               if (cabs1(qraux(j)) .eq. 0.0e0) go to 50
                  temp = x(j,j)
                  x(j,j) = qraux(j)
                  t = -cdotc(n-j+1,x(j,j),1,qy(j),1)/x(j,j)
                  call caxpy(n-j+1,t,x(j,j),1,qy(j),1)
                  x(j,j) = temp
   50          continue
   60       continue
   70    continue
         if (.not.cqty) go to 100
c
c           compute ctrans(q)*y.
c
            do 90 j = 1, ju
               if (cabs1(qraux(j)) .eq. 0.0e0) go to 80
                  temp = x(j,j)
                  x(j,j) = qraux(j)
                  t = -cdotc(n-j+1,x(j,j),1,qty(j),1)/x(j,j)
                  call caxpy(n-j+1,t,x(j,j),1,qty(j),1)
                  x(j,j) = temp
   80          continue
   90       continue
  100    continue
c
c        set up to compute b, rsd, or xb.
c
         if (cb) call ccopy(k,qty,1,b,1)
         kp1 = k + 1
         if (cxb) call ccopy(k,qty,1,xb,1)
         if (cr .and. k .lt. n) call ccopy(n-k,qty(kp1),1,rsd(kp1),1)
         if (.not.cxb .or. kp1 .gt. n) go to 120
            do 110 i = kp1, n
               xb(i) = (0.0e0,0.0e0)
  110       continue
  120    continue
         if (.not.cr) go to 140
            do 130 i = 1, k
               rsd(i) = (0.0e0,0.0e0)
  130       continue
  140    continue
         if (.not.cb) go to 190
c
c           compute b.
c
            do 170 jj = 1, k
               j = k - jj + 1
               if (cabs1(x(j,j)) .ne. 0.0e0) go to 150
                  info = j
c           ......exit
                  go to 180
  150          continue
               b(j) = b(j)/x(j,j)
               if (j .eq. 1) go to 160
                  t = -b(j)
                  call caxpy(j-1,t,x(1,j),1,b,1)
  160          continue
  170       continue
  180       continue
  190    continue
         if (.not.cr .and. .not.cxb) go to 240
c
c           compute rsd or xb as required.
c
            do 230 jj = 1, ju
               j = ju - jj + 1
               if (cabs1(qraux(j)) .eq. 0.0e0) go to 220
                  temp = x(j,j)
                  x(j,j) = qraux(j)
                  if (.not.cr) go to 200
                     t = -cdotc(n-j+1,x(j,j),1,rsd(j),1)/x(j,j)
                     call caxpy(n-j+1,t,x(j,j),1,rsd(j),1)
  200             continue
                  if (.not.cxb) go to 210
                     t = -cdotc(n-j+1,x(j,j),1,xb(j),1)/x(j,j)
                     call caxpy(n-j+1,t,x(j,j),1,xb(j),1)
  210             continue
                  x(j,j) = temp
  220          continue
  230       continue
  240    continue
  250 continue
      return
      end
      subroutine csico(a,lda,n,kpvt,rcond,z)

c*********************************************************************72

      integer lda,n,kpvt(1)
      complex a(lda,1),z(1)
      real rcond
c
c     csico factors a complex symmetric matrix by elimination with
c     symmetric pivoting and estimates the condition of the matrix.
c
c     if  rcond  is not needed, csifa is slightly faster.
c     to solve  a*x = b , follow csico by csisl.
c     to compute  inverse(a)*c , follow csico by csisl.
c     to compute  inverse(a) , follow csico by csidi.
c     to compute  determinant(a) , follow csico by csidi.
c
c     on entry
c
c        a       complex(lda, n)
c                the symmetric matrix to be factored.
c                only the diagonal and upper triangle are used.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       a block diagonal matrix and the multipliers which
c                were used to obtain it.
c                the factorization can be written  a = u*d*trans(u)
c                where  u  is a product of permutation and unit
c                upper triangular matrices , trans(u) is the
c                transpose of  u , and  d  is block diagonal
c                with 1 by 1 and 2 by 2 blocks.
c
c        kpvt    integer(n)
c                an integer vector of pivot indices.
c
c        rcond   real
c                an estimate of the reciprocal condition of  a .
c                for the system  a*x = b , relative perturbations
c                in  a  and  b  of size  epsilon  may cause
c                relative perturbations in  x  of size  epsilon/rcond .
c                if  rcond  is so small that the logical expression
c                           1.0 + rcond .eq. 1.0
c                is true, then  a  may be singular to working
c                precision.  in particular,  rcond  is zero  if
c                exact singularity is detected or the estimate
c                underflows.
c
c        z       complex(n)
c                a work vector whose contents are usually unimportant.
c                if  a  is close to a singular matrix, then  z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     linpack csifa
c     blas caxpy,cdotu,csscal,scasum
c     fortran abs,aimag,amax1,cmplx,iabs,real
c
c     internal variables
c
      complex ak,akm1,bk,bkm1,cdotu,denom,ek,t
      real anorm,s,scasum,ynorm
      integer i,info,j,jm1,k,kp,kps,ks
c
      complex zdum,zdum2,csign1
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
      csign1(zdum,zdum2) = cabs1(zdum)*(zdum2/cabs1(zdum2))
c
c     find norm of a using only upper half
c
      do 30 j = 1, n
         z(j) = cmplx(scasum(j,a(1,j),1),0.0e0)
         jm1 = j - 1
         if (jm1 .lt. 1) go to 20
         do 10 i = 1, jm1
            z(i) = cmplx(real(z(i))+cabs1(a(i,j)),0.0e0)
   10    continue
   20    continue
   30 continue
      anorm = 0.0e0
      do 40 j = 1, n
         anorm = amax1(anorm,real(z(j)))
   40 continue
c
c     factor
c
      call csifa(a,lda,n,kpvt,info)
c
c     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
c     estimate = norm(z)/norm(y) where  a*z = y  and  a*y = e .
c     the components of  e  are chosen to cause maximum local
c     growth in the elements of w  where  u*d*w = e .
c     the vectors are frequently rescaled to avoid overflow.
c
c     solve u*d*w = e
c
      ek = (1.0e0,0.0e0)
      do 50 j = 1, n
         z(j) = (0.0e0,0.0e0)
   50 continue
      k = n
   60 if (k .eq. 0) go to 120
         ks = 1
         if (kpvt(k) .lt. 0) ks = 2
         kp = iabs(kpvt(k))
         kps = k + 1 - ks
         if (kp .eq. kps) go to 70
            t = z(kps)
            z(kps) = z(kp)
            z(kp) = t
   70    continue
         if (cabs1(z(k)) .ne. 0.0e0) ek = csign1(ek,z(k))
         z(k) = z(k) + ek
         call caxpy(k-ks,z(k),a(1,k),1,z(1),1)
         if (ks .eq. 1) go to 80
            if (cabs1(z(k-1)) .ne. 0.0e0) ek = csign1(ek,z(k-1))
            z(k-1) = z(k-1) + ek
            call caxpy(k-ks,z(k-1),a(1,k-1),1,z(1),1)
   80    continue
         if (ks .eq. 2) go to 100
            if (cabs1(z(k)) .le. cabs1(a(k,k))) go to 90
               s = cabs1(a(k,k))/cabs1(z(k))
               call csscal(n,s,z,1)
               ek = cmplx(s,0.0e0)*ek
   90       continue
            if (cabs1(a(k,k)) .ne. 0.0e0) z(k) = z(k)/a(k,k)
            if (cabs1(a(k,k)) .eq. 0.0e0) z(k) = (1.0e0,0.0e0)
         go to 110
  100    continue
            ak = a(k,k)/a(k-1,k)
            akm1 = a(k-1,k-1)/a(k-1,k)
            bk = z(k)/a(k-1,k)
            bkm1 = z(k-1)/a(k-1,k)
            denom = ak*akm1 - 1.0e0
            z(k) = (akm1*bk - bkm1)/denom
            z(k-1) = (ak*bkm1 - bk)/denom
  110    continue
         k = k - ks
      go to 60
  120 continue
      s = 1.0e0/scasum(n,z,1)
      call csscal(n,s,z,1)
c
c     solve trans(u)*y = w
c
      k = 1
  130 if (k .gt. n) go to 160
         ks = 1
         if (kpvt(k) .lt. 0) ks = 2
         if (k .eq. 1) go to 150
            z(k) = z(k) + cdotu(k-1,a(1,k),1,z(1),1)
            if (ks .eq. 2)
     *         z(k+1) = z(k+1) + cdotu(k-1,a(1,k+1),1,z(1),1)
            kp = iabs(kpvt(k))
            if (kp .eq. k) go to 140
               t = z(k)
               z(k) = z(kp)
               z(kp) = t
  140       continue
  150    continue
         k = k + ks
      go to 130
  160 continue
      s = 1.0e0/scasum(n,z,1)
      call csscal(n,s,z,1)
c
      ynorm = 1.0e0
c
c     solve u*d*v = y
c
      k = n
  170 if (k .eq. 0) go to 230
         ks = 1
         if (kpvt(k) .lt. 0) ks = 2
         if (k .eq. ks) go to 190
            kp = iabs(kpvt(k))
            kps = k + 1 - ks
            if (kp .eq. kps) go to 180
               t = z(kps)
               z(kps) = z(kp)
               z(kp) = t
  180       continue
            call caxpy(k-ks,z(k),a(1,k),1,z(1),1)
            if (ks .eq. 2) call caxpy(k-ks,z(k-1),a(1,k-1),1,z(1),1)
  190    continue
         if (ks .eq. 2) go to 210
            if (cabs1(z(k)) .le. cabs1(a(k,k))) go to 200
               s = cabs1(a(k,k))/cabs1(z(k))
               call csscal(n,s,z,1)
               ynorm = s*ynorm
  200       continue
            if (cabs1(a(k,k)) .ne. 0.0e0) z(k) = z(k)/a(k,k)
            if (cabs1(a(k,k)) .eq. 0.0e0) z(k) = (1.0e0,0.0e0)
         go to 220
  210    continue
            ak = a(k,k)/a(k-1,k)
            akm1 = a(k-1,k-1)/a(k-1,k)
            bk = z(k)/a(k-1,k)
            bkm1 = z(k-1)/a(k-1,k)
            denom = ak*akm1 - 1.0e0
            z(k) = (akm1*bk - bkm1)/denom
            z(k-1) = (ak*bkm1 - bk)/denom
  220    continue
         k = k - ks
      go to 170
  230 continue
      s = 1.0e0/scasum(n,z,1)
      call csscal(n,s,z,1)
      ynorm = s*ynorm
c
c     solve trans(u)*z = v
c
      k = 1
  240 if (k .gt. n) go to 270
         ks = 1
         if (kpvt(k) .lt. 0) ks = 2
         if (k .eq. 1) go to 260
            z(k) = z(k) + cdotu(k-1,a(1,k),1,z(1),1)
            if (ks .eq. 2)
     *         z(k+1) = z(k+1) + cdotu(k-1,a(1,k+1),1,z(1),1)
            kp = iabs(kpvt(k))
            if (kp .eq. k) go to 250
               t = z(k)
               z(k) = z(kp)
               z(kp) = t
  250       continue
  260    continue
         k = k + ks
      go to 240
  270 continue
c     make znorm = 1.0
      s = 1.0e0/scasum(n,z,1)
      call csscal(n,s,z,1)
      ynorm = s*ynorm
c
      if (anorm .ne. 0.0e0) rcond = ynorm/anorm
      if (anorm .eq. 0.0e0) rcond = 0.0e0
      return
      end
      subroutine csidi(a,lda,n,kpvt,det,work,job)

c*********************************************************************72

      integer lda,n,job
      complex a(lda,1),det(2),work(1)
      integer kpvt(1)
c
c     csidi computes the determinant and inverse
c     of a complex symmetric matrix using the factors from csifa.
c
c     on entry
c
c        a       complex(lda,n)
c                the output from csifa.
c
c        lda     integer
c                the leading dimension of the array a.
c
c        n       integer
c                the order of the matrix a.
c
c        kpvt    integer(n)
c                the pivot vector from csifa.
c
c        work    complex(n)
c                work vector.  contents destroyed.
c
c        job     integer
c                job has the decimal expansion  ab  where
c                   if  b .ne. 0, the inverse is computed,
c                   if  a .ne. 0, the determinant is computed,
c
c                for example, job = 11  gives both.
c
c     on return
c
c        variables not requested by job are not used.
c
c        a      contains the upper triangle of the inverse of
c               the original matrix.  the strict lower triangle
c               is never referenced.
c
c        det    complex(2)
c               determinant of original matrix.
c               determinant = det(1) * 10.0**det(2)
c               with 1.0 .le. abs(det(1)) .lt. 10.0
c               or det(1) = 0.0.
c
c     error condition
c
c        a division by zero may occur if the inverse is requested
c        and  csico  has set rcond .eq. 0.0
c        or  csifa  has set  info .ne. 0 .
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab
c
c     subroutines and functions
c
c     blas caxpy,ccopy,cdotu,cswap
c     fortran abs,cmplx,iabs,mod,real
c
c     internal variables.
c
      complex ak,akp1,akkp1,cdotu,d,t,temp
      real ten
      integer j,jb,k,km1,ks,kstep
      logical noinv,nodet
c
      complex zdum
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
c
      noinv = mod(job,10) .eq. 0
      nodet = mod(job,100)/10 .eq. 0
c
      if (nodet) go to 100
         det(1) = (1.0e0,0.0e0)
         det(2) = (0.0e0,0.0e0)
         ten = 10.0e0
         t = (0.0e0,0.0e0)
         do 90 k = 1, n
            d = a(k,k)
c
c           check if 1 by 1
c
            if (kpvt(k) .gt. 0) go to 30
c
c              2 by 2 block
c              use det (d  t)  =  (d/t * c - t) * t
c                      (t  c)
c              to avoid underflow/overflow troubles.
c              take two passes through scaling.  use  t  for flag.
c
               if (cabs1(t) .ne. 0.0e0) go to 10
                  t = a(k,k+1)
                  d = (d/t)*a(k+1,k+1) - t
               go to 20
   10          continue
                  d = t
                  t = (0.0e0,0.0e0)
   20          continue
   30       continue
c
            det(1) = d*det(1)
            if (cabs1(det(1)) .eq. 0.0e0) go to 80
   40          if (cabs1(det(1)) .ge. 1.0e0) go to 50
                  det(1) = cmplx(ten,0.0e0)*det(1)
                  det(2) = det(2) - (1.0e0,0.0e0)
               go to 40
   50          continue
   60          if (cabs1(det(1)) .lt. ten) go to 70
                  det(1) = det(1)/cmplx(ten,0.0e0)
                  det(2) = det(2) + (1.0e0,0.0e0)
               go to 60
   70          continue
   80       continue
   90    continue
  100 continue
c
c     compute inverse(a)
c
      if (noinv) go to 230
         k = 1
  110    if (k .gt. n) go to 220
            km1 = k - 1
            if (kpvt(k) .lt. 0) go to 140
c
c              1 by 1
c
               a(k,k) = (1.0e0,0.0e0)/a(k,k)
               if (km1 .lt. 1) go to 130
                  call ccopy(km1,a(1,k),1,work,1)
                  do 120 j = 1, km1
                     a(j,k) = cdotu(j,a(1,j),1,work,1)
                     call caxpy(j-1,work(j),a(1,j),1,a(1,k),1)
  120             continue
                  a(k,k) = a(k,k) + cdotu(km1,work,1,a(1,k),1)
  130          continue
               kstep = 1
            go to 180
  140       continue
c
c              2 by 2
c
               t = a(k,k+1)
               ak = a(k,k)/t
               akp1 = a(k+1,k+1)/t
               akkp1 = a(k,k+1)/t
               d = t*(ak*akp1 - (1.0e0,0.0e0))
               a(k,k) = akp1/d
               a(k+1,k+1) = ak/d
               a(k,k+1) = -akkp1/d
               if (km1 .lt. 1) go to 170
                  call ccopy(km1,a(1,k+1),1,work,1)
                  do 150 j = 1, km1
                     a(j,k+1) = cdotu(j,a(1,j),1,work,1)
                     call caxpy(j-1,work(j),a(1,j),1,a(1,k+1),1)
  150             continue
                  a(k+1,k+1) = a(k+1,k+1)
     *                         + cdotu(km1,work,1,a(1,k+1),1)
                  a(k,k+1) = a(k,k+1) + cdotu(km1,a(1,k),1,a(1,k+1),1)
                  call ccopy(km1,a(1,k),1,work,1)
                  do 160 j = 1, km1
                     a(j,k) = cdotu(j,a(1,j),1,work,1)
                     call caxpy(j-1,work(j),a(1,j),1,a(1,k),1)
  160             continue
                  a(k,k) = a(k,k) + cdotu(km1,work,1,a(1,k),1)
  170          continue
               kstep = 2
  180       continue
c
c           swap
c
            ks = iabs(kpvt(k))
            if (ks .eq. k) go to 210
               call cswap(ks,a(1,ks),1,a(1,k),1)
               do 190 jb = ks, k
                  j = k + ks - jb
                  temp = a(j,k)
                  a(j,k) = a(ks,j)
                  a(ks,j) = temp
  190          continue
               if (kstep .eq. 1) go to 200
                  temp = a(ks,k+1)
                  a(ks,k+1) = a(k,k+1)
                  a(k,k+1) = temp
  200          continue
  210       continue
            k = k + kstep
         go to 110
  220    continue
  230 continue
      return
      end
      subroutine csifa(a,lda,n,kpvt,info)

c*********************************************************************72

      integer lda,n,kpvt(1),info
      complex a(lda,1)
c
c     csifa factors a complex symmetric matrix by elimination
c     with symmetric pivoting.
c
c     to solve  a*x = b , follow csifa by csisl.
c     to compute  inverse(a)*c , follow csifa by csisl.
c     to compute  determinant(a) , follow csifa by csidi.
c     to compute  inverse(a) , follow csifa by csidi.
c
c     on entry
c
c        a       complex(lda,n)
c                the symmetric matrix to be factored.
c                only the diagonal and upper triangle are used.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       a block diagonal matrix and the multipliers which
c                were used to obtain it.
c                the factorization can be written  a = u*d*trans(u)
c                where  u  is a product of permutation and unit
c                upper triangular matrices , trans(u) is the
c                transpose of  u , and  d  is block diagonal
c                with 1 by 1 and 2 by 2 blocks.
c
c        kpvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if the k-th pivot block is singular. this is
c                     not an error condition for this subroutine,
c                     but it does indicate that csisl or csidi may
c                     divide by zero if called.
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab.
c
c     subroutines and functions
c
c     blas caxpy,cswap,icamax
c     fortran abs,aimag,amax1,real,sqrt
c
c     internal variables
c
      complex ak,akm1,bk,bkm1,denom,mulk,mulkm1,t
      real absakk,alpha,colmax,rowmax
      integer imax,imaxp1,j,jj,jmax,k,km1,km2,kstep,icamax
      logical swap
c
      complex zdum
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
c
c     initialize
c
c     alpha is used in choosing pivot block size.
      alpha = (1.0e0 + sqrt(17.0e0))/8.0e0
c
      info = 0
c
c     main loop on k, which goes from n to 1.
c
      k = n
   10 continue
c
c        leave the loop if k=0 or k=1.
c
c     ...exit
         if (k .eq. 0) go to 200
         if (k .gt. 1) go to 20
            kpvt(1) = 1
            if (cabs1(a(1,1)) .eq. 0.0e0) info = 1
c     ......exit
            go to 200
   20    continue
c
c        this section of code determines the kind of
c        elimination to be performed.  when it is completed,
c        kstep will be set to the size of the pivot block, and
c        swap will be set to .true. if an interchange is
c        required.
c
         km1 = k - 1
         absakk = cabs1(a(k,k))
c
c        determine the largest off-diagonal element in
c        column k.
c
         imax = icamax(k-1,a(1,k),1)
         colmax = cabs1(a(imax,k))
         if (absakk .lt. alpha*colmax) go to 30
            kstep = 1
            swap = .false.
         go to 90
   30    continue
c
c           determine the largest off-diagonal element in
c           row imax.
c
            rowmax = 0.0e0
            imaxp1 = imax + 1
            do 40 j = imaxp1, k
               rowmax = amax1(rowmax,cabs1(a(imax,j)))
   40       continue
            if (imax .eq. 1) go to 50
               jmax = icamax(imax-1,a(1,imax),1)
               rowmax = amax1(rowmax,cabs1(a(jmax,imax)))
   50       continue
            if (cabs1(a(imax,imax)) .lt. alpha*rowmax) go to 60
               kstep = 1
               swap = .true.
            go to 80
   60       continue
            if (absakk .lt. alpha*colmax*(colmax/rowmax)) go to 70
               kstep = 1
               swap = .false.
            go to 80
   70       continue
               kstep = 2
               swap = imax .ne. km1
   80       continue
   90    continue
         if (amax1(absakk,colmax) .ne. 0.0e0) go to 100
c
c           column k is zero.  set info and iterate the loop.
c
            kpvt(k) = k
            info = k
         go to 190
  100    continue
         if (kstep .eq. 2) go to 140
c
c           1 x 1 pivot block.
c
            if (.not.swap) go to 120
c
c              perform an interchange.
c
               call cswap(imax,a(1,imax),1,a(1,k),1)
               do 110 jj = imax, k
                  j = k + imax - jj
                  t = a(j,k)
                  a(j,k) = a(imax,j)
                  a(imax,j) = t
  110          continue
  120       continue
c
c           perform the elimination.
c
            do 130 jj = 1, km1
               j = k - jj
               mulk = -a(j,k)/a(k,k)
               t = mulk
               call caxpy(j,t,a(1,k),1,a(1,j),1)
               a(j,k) = mulk
  130       continue
c
c           set the pivot array.
c
            kpvt(k) = k
            if (swap) kpvt(k) = imax
         go to 190
  140    continue
c
c           2 x 2 pivot block.
c
            if (.not.swap) go to 160
c
c              perform an interchange.
c
               call cswap(imax,a(1,imax),1,a(1,k-1),1)
               do 150 jj = imax, km1
                  j = km1 + imax - jj
                  t = a(j,k-1)
                  a(j,k-1) = a(imax,j)
                  a(imax,j) = t
  150          continue
               t = a(k-1,k)
               a(k-1,k) = a(imax,k)
               a(imax,k) = t
  160       continue
c
c           perform the elimination.
c
            km2 = k - 2
            if (km2 .eq. 0) go to 180
               ak = a(k,k)/a(k-1,k)
               akm1 = a(k-1,k-1)/a(k-1,k)
               denom = 1.0e0 - ak*akm1
               do 170 jj = 1, km2
                  j = km1 - jj
                  bk = a(j,k)/a(k-1,k)
                  bkm1 = a(j,k-1)/a(k-1,k)
                  mulk = (akm1*bk - bkm1)/denom
                  mulkm1 = (ak*bkm1 - bk)/denom
                  t = mulk
                  call caxpy(j,t,a(1,k),1,a(1,j),1)
                  t = mulkm1
                  call caxpy(j,t,a(1,k-1),1,a(1,j),1)
                  a(j,k) = mulk
                  a(j,k-1) = mulkm1
  170          continue
  180       continue
c
c           set the pivot array.
c
            kpvt(k) = 1 - k
            if (swap) kpvt(k) = -imax
            kpvt(k-1) = kpvt(k)
  190    continue
         k = k - kstep
      go to 10
  200 continue
      return
      end
      subroutine csisl(a,lda,n,kpvt,b)

c*********************************************************************72

      integer lda,n,kpvt(1)
      complex a(lda,1),b(1)
c
c     csisl solves the complex symmetric system
c     a * x = b
c     using the factors computed by csifa.
c
c     on entry
c
c        a       complex(lda,n)
c                the output from csifa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        kpvt    integer(n)
c                the pivot vector from csifa.
c
c        b       complex(n)
c                the right hand side vector.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero may occur if  csico  has set rcond .eq. 0.0
c        or  csifa  has set info .ne. 0  .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call csifa(a,lda,n,kpvt,info)
c           if (info .ne. 0) go to ...
c           do 10 j = 1, p
c              call csisl(a,lda,n,kpvt,c(1,j))
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab.
c
c     subroutines and functions
c
c     blas caxpy,cdotu
c     fortran iabs
c
c     internal variables.
c
      complex ak,akm1,bk,bkm1,cdotu,denom,temp
      integer k,kp
c
c     loop backward applying the transformations and
c     d inverse to b.
c
      k = n
   10 if (k .eq. 0) go to 80
         if (kpvt(k) .lt. 0) go to 40
c
c           1 x 1 pivot block.
c
            if (k .eq. 1) go to 30
               kp = kpvt(k)
               if (kp .eq. k) go to 20
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
   20          continue
c
c              apply the transformation.
c
               call caxpy(k-1,b(k),a(1,k),1,b(1),1)
   30       continue
c
c           apply d inverse.
c
            b(k) = b(k)/a(k,k)
            k = k - 1
         go to 70
   40    continue
c
c           2 x 2 pivot block.
c
            if (k .eq. 2) go to 60
               kp = iabs(kpvt(k))
               if (kp .eq. k - 1) go to 50
c
c                 interchange.
c
                  temp = b(k-1)
                  b(k-1) = b(kp)
                  b(kp) = temp
   50          continue
c
c              apply the transformation.
c
               call caxpy(k-2,b(k),a(1,k),1,b(1),1)
               call caxpy(k-2,b(k-1),a(1,k-1),1,b(1),1)
   60       continue
c
c           apply d inverse.
c
            ak = a(k,k)/a(k-1,k)
            akm1 = a(k-1,k-1)/a(k-1,k)
            bk = b(k)/a(k-1,k)
            bkm1 = b(k-1)/a(k-1,k)
            denom = ak*akm1 - 1.0e0
            b(k) = (akm1*bk - bkm1)/denom
            b(k-1) = (ak*bkm1 - bk)/denom
            k = k - 2
   70    continue
      go to 10
   80 continue
c
c     loop forward applying the transformations.
c
      k = 1
   90 if (k .gt. n) go to 160
         if (kpvt(k) .lt. 0) go to 120
c
c           1 x 1 pivot block.
c
            if (k .eq. 1) go to 110
c
c              apply the transformation.
c
               b(k) = b(k) + cdotu(k-1,a(1,k),1,b(1),1)
               kp = kpvt(k)
               if (kp .eq. k) go to 100
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
  100          continue
  110       continue
            k = k + 1
         go to 150
  120    continue
c
c           2 x 2 pivot block.
c
            if (k .eq. 1) go to 140
c
c              apply the transformation.
c
               b(k) = b(k) + cdotu(k-1,a(1,k),1,b(1),1)
               b(k+1) = b(k+1) + cdotu(k-1,a(1,k+1),1,b(1),1)
               kp = iabs(kpvt(k))
               if (kp .eq. k) go to 130
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
  130          continue
  140       continue
            k = k + 2
  150    continue
      go to 90
  160 continue
      return
      end
      subroutine cspco(ap,n,kpvt,rcond,z)

c*********************************************************************72

      integer n,kpvt(1)
      complex ap(1),z(1)
      real rcond
c
c     cspco factors a complex symmetric matrix stored in packed
c     form by elimination with symmetric pivoting and estimates
c     the condition of the matrix.
c
c     if  rcond  is not needed, cspfa is slightly faster.
c     to solve  a*x = b , follow cspco by cspsl.
c     to compute  inverse(a)*c , follow cspco by cspsl.
c     to compute  inverse(a) , follow cspco by cspdi.
c     to compute  determinant(a) , follow cspco by cspdi.
c
c     on entry
c
c        ap      complex (n*(n+1)/2)
c                the packed form of a symmetric matrix  a .  the
c                columns of the upper triangle are stored sequentially
c                in a one-dimensional array of length  n*(n+1)/2 .
c                see comments below for details.
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        ap      a block diagonal matrix and the multipliers which
c                were used to obtain it stored in packed form.
c                the factorization can be written  a = u*d*trans(u)
c                where  u  is a product of permutation and unit
c                upper triangular matrices , trans(u) is the
c                transpose of  u , and  d  is block diagonal
c                with 1 by 1 and 2 by 2 blocks.
c
c        kpvt    integer(n)
c                an integer vector of pivot indices.
c
c        rcond   real
c                an estimate of the reciprocal condition of  a .
c                for the system  a*x = b , relative perturbations
c                in  a  and  b  of size  epsilon  may cause
c                relative perturbations in  x  of size  epsilon/rcond .
c                if  rcond  is so small that the logical expression
c                           1.0 + rcond .eq. 1.0
c                is true, then  a  may be singular to working
c                precision.  in particular,  rcond  is zero  if
c                exact singularity is detected or the estimate
c                underflows.
c
c        z       complex(n)
c                a work vector whose contents are usually unimportant.
c                if  a  is close to a singular matrix, then  z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .
c
c     packed storage
c
c          the following program segment will pack the upper
c          triangle of a symmetric matrix.
c
c                k = 0
c                do 20 j = 1, n
c                   do 10 i = 1, j
c                      k = k + 1
c                      ap(k) = a(i,j)
c             10    continue
c             20 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     linpack cspfa
c     blas caxpy,cdotu,csscal,scasum
c     fortran abs,aimag,amax1,cmplx,iabs,real
c
c     internal variables
c
      complex ak,akm1,bk,bkm1,cdotu,denom,ek,t
      real anorm,s,scasum,ynorm
      integer i,ij,ik,ikm1,ikp1,info,j,jm1,j1
      integer k,kk,km1k,km1km1,kp,kps,ks
c
      complex zdum,zdum2,csign1
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
      csign1(zdum,zdum2) = cabs1(zdum)*(zdum2/cabs1(zdum2))
c
c     find norm of a using only upper half
c
      j1 = 1
      do 30 j = 1, n
         z(j) = cmplx(scasum(j,ap(j1),1),0.0e0)
         ij = j1
         j1 = j1 + j
         jm1 = j - 1
         if (jm1 .lt. 1) go to 20
         do 10 i = 1, jm1
            z(i) = cmplx(real(z(i))+cabs1(ap(ij)),0.0e0)
            ij = ij + 1
   10    continue
   20    continue
   30 continue
      anorm = 0.0e0
      do 40 j = 1, n
         anorm = amax1(anorm,real(z(j)))
   40 continue
c
c     factor
c
      call cspfa(ap,n,kpvt,info)
c
c     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
c     estimate = norm(z)/norm(y) where  a*z = y  and  a*y = e .
c     the components of  e  are chosen to cause maximum local
c     growth in the elements of w  where  u*d*w = e .
c     the vectors are frequently rescaled to avoid overflow.
c
c     solve u*d*w = e
c
      ek = (1.0e0,0.0e0)
      do 50 j = 1, n
         z(j) = (0.0e0,0.0e0)
   50 continue
      k = n
      ik = (n*(n - 1))/2
   60 if (k .eq. 0) go to 120
         kk = ik + k
         ikm1 = ik - (k - 1)
         ks = 1
         if (kpvt(k) .lt. 0) ks = 2
         kp = iabs(kpvt(k))
         kps = k + 1 - ks
         if (kp .eq. kps) go to 70
            t = z(kps)
            z(kps) = z(kp)
            z(kp) = t
   70    continue
         if (cabs1(z(k)) .ne. 0.0e0) ek = csign1(ek,z(k))
         z(k) = z(k) + ek
         call caxpy(k-ks,z(k),ap(ik+1),1,z(1),1)
         if (ks .eq. 1) go to 80
            if (cabs1(z(k-1)) .ne. 0.0e0) ek = csign1(ek,z(k-1))
            z(k-1) = z(k-1) + ek
            call caxpy(k-ks,z(k-1),ap(ikm1+1),1,z(1),1)
   80    continue
         if (ks .eq. 2) go to 100
            if (cabs1(z(k)) .le. cabs1(ap(kk))) go to 90
               s = cabs1(ap(kk))/cabs1(z(k))
               call csscal(n,s,z,1)
               ek = cmplx(s,0.0e0)*ek
   90       continue
            if (cabs1(ap(kk)) .ne. 0.0e0) z(k) = z(k)/ap(kk)
            if (cabs1(ap(kk)) .eq. 0.0e0) z(k) = (1.0e0,0.0e0)
         go to 110
  100    continue
            km1k = ik + k - 1
            km1km1 = ikm1 + k - 1
            ak = ap(kk)/ap(km1k)
            akm1 = ap(km1km1)/ap(km1k)
            bk = z(k)/ap(km1k)
            bkm1 = z(k-1)/ap(km1k)
            denom = ak*akm1 - 1.0e0
            z(k) = (akm1*bk - bkm1)/denom
            z(k-1) = (ak*bkm1 - bk)/denom
  110    continue
         k = k - ks
         ik = ik - k
         if (ks .eq. 2) ik = ik - (k + 1)
      go to 60
  120 continue
      s = 1.0e0/scasum(n,z,1)
      call csscal(n,s,z,1)
c
c     solve trans(u)*y = w
c
      k = 1
      ik = 0
  130 if (k .gt. n) go to 160
         ks = 1
         if (kpvt(k) .lt. 0) ks = 2
         if (k .eq. 1) go to 150
            z(k) = z(k) + cdotu(k-1,ap(ik+1),1,z(1),1)
            ikp1 = ik + k
            if (ks .eq. 2)
     *         z(k+1) = z(k+1) + cdotu(k-1,ap(ikp1+1),1,z(1),1)
            kp = iabs(kpvt(k))
            if (kp .eq. k) go to 140
               t = z(k)
               z(k) = z(kp)
               z(kp) = t
  140       continue
  150    continue
         ik = ik + k
         if (ks .eq. 2) ik = ik + (k + 1)
         k = k + ks
      go to 130
  160 continue
      s = 1.0e0/scasum(n,z,1)
      call csscal(n,s,z,1)
c
      ynorm = 1.0e0
c
c     solve u*d*v = y
c
      k = n
      ik = n*(n - 1)/2
  170 if (k .eq. 0) go to 230
         kk = ik + k
         ikm1 = ik - (k - 1)
         ks = 1
         if (kpvt(k) .lt. 0) ks = 2
         if (k .eq. ks) go to 190
            kp = iabs(kpvt(k))
            kps = k + 1 - ks
            if (kp .eq. kps) go to 180
               t = z(kps)
               z(kps) = z(kp)
               z(kp) = t
  180       continue
            call caxpy(k-ks,z(k),ap(ik+1),1,z(1),1)
            if (ks .eq. 2) call caxpy(k-ks,z(k-1),ap(ikm1+1),1,z(1),1)
  190    continue
         if (ks .eq. 2) go to 210
            if (cabs1(z(k)) .le. cabs1(ap(kk))) go to 200
               s = cabs1(ap(kk))/cabs1(z(k))
               call csscal(n,s,z,1)
               ynorm = s*ynorm
  200       continue
            if (cabs1(ap(kk)) .ne. 0.0e0) z(k) = z(k)/ap(kk)
            if (cabs1(ap(kk)) .eq. 0.0e0) z(k) = (1.0e0,0.0e0)
         go to 220
  210    continue
            km1k = ik + k - 1
            km1km1 = ikm1 + k - 1
            ak = ap(kk)/ap(km1k)
            akm1 = ap(km1km1)/ap(km1k)
            bk = z(k)/ap(km1k)
            bkm1 = z(k-1)/ap(km1k)
            denom = ak*akm1 - 1.0e0
            z(k) = (akm1*bk - bkm1)/denom
            z(k-1) = (ak*bkm1 - bk)/denom
  220    continue
         k = k - ks
         ik = ik - k
         if (ks .eq. 2) ik = ik - (k + 1)
      go to 170
  230 continue
      s = 1.0e0/scasum(n,z,1)
      call csscal(n,s,z,1)
      ynorm = s*ynorm
c
c     solve trans(u)*z = v
c
      k = 1
      ik = 0
  240 if (k .gt. n) go to 270
         ks = 1
         if (kpvt(k) .lt. 0) ks = 2
         if (k .eq. 1) go to 260
            z(k) = z(k) + cdotu(k-1,ap(ik+1),1,z(1),1)
            ikp1 = ik + k
            if (ks .eq. 2)
     *         z(k+1) = z(k+1) + cdotu(k-1,ap(ikp1+1),1,z(1),1)
            kp = iabs(kpvt(k))
            if (kp .eq. k) go to 250
               t = z(k)
               z(k) = z(kp)
               z(kp) = t
  250       continue
  260    continue
         ik = ik + k
         if (ks .eq. 2) ik = ik + (k + 1)
         k = k + ks
      go to 240
  270 continue
c     make znorm = 1.0
      s = 1.0e0/scasum(n,z,1)
      call csscal(n,s,z,1)
      ynorm = s*ynorm
c
      if (anorm .ne. 0.0e0) rcond = ynorm/anorm
      if (anorm .eq. 0.0e0) rcond = 0.0e0
      return
      end
      subroutine cspdi(ap,n,kpvt,det,work,job)

c*********************************************************************72

      integer n,job
      complex ap(1),work(1),det(2)
      integer kpvt(1)
c
c     cspdi computes the determinant and inverse
c     of a complex symmetric matrix using the factors from cspfa,
c     where the matrix is stored in packed form.
c
c     on entry
c
c        ap      complex (n*(n+1)/2)
c                the output from cspfa.
c
c        n       integer
c                the order of the matrix a.
c
c        kpvt    integer(n)
c                the pivot vector from cspfa.
c
c        work    complex(n)
c                work vector.  contents ignored.
c
c        job     integer
c                job has the decimal expansion  ab  where
c                   if  b .ne. 0, the inverse is computed,
c                   if  a .ne. 0, the determinant is computed,
c
c                for example, job = 11  gives both.
c
c     on return
c
c        variables not requested by job are not used.
c
c        ap     contains the upper triangle of the inverse of
c               the original matrix, stored in packed form.
c               the columns of the upper triangle are stored
c               sequentially in a one-dimensional array.
c
c        det    complex(2)
c               determinant of original matrix.
c               determinant = det(1) * 10.0**det(2)
c               with 1.0 .le. abs(det(1)) .lt. 10.0
c               or det(1) = 0.0.
c
c     error condition
c
c        a division by zero will occur if the inverse is requested
c        and  cspco  has set rcond .eq. 0.0
c        or  cspfa  has set  info .ne. 0 .
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab.
c
c     subroutines and functions
c
c     blas caxpy,ccopy,cdotu,cswap
c     fortran abs,cmplx,iabs,mod,real
c
c     internal variables.
c
      complex ak,akkp1,akp1,cdotu,d,t,temp
      real ten
      integer ij,ik,ikp1,iks,j,jb,jk,jkp1
      integer k,kk,kkp1,km1,ks,ksj,kskp1,kstep
      logical noinv,nodet
c
      complex zdum
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
c
      noinv = mod(job,10) .eq. 0
      nodet = mod(job,100)/10 .eq. 0
c
      if (nodet) go to 110
         det(1) = (1.0e0,0.0e0)
         det(2) = (0.0e0,0.0e0)
         ten = 10.0e0
         t = (0.0e0,0.0e0)
         ik = 0
         do 100 k = 1, n
            kk = ik + k
            d = ap(kk)
c
c           check if 1 by 1
c
            if (kpvt(k) .gt. 0) go to 30
c
c              2 by 2 block
c              use det (d  t)  =  (d/t * c - t) * t
c                      (t  c)
c              to avoid underflow/overflow troubles.
c              take two passes through scaling.  use  t  for flag.
c
               if (cabs1(t) .ne. 0.0e0) go to 10
                  ikp1 = ik + k
                  kkp1 = ikp1 + k
                  t = ap(kkp1)
                  d = (d/t)*ap(kkp1+1) - t
               go to 20
   10          continue
                  d = t
                  t = (0.0e0,0.0e0)
   20          continue
   30       continue
c
            if (nodet) go to 90
               det(1) = d*det(1)
               if (cabs1(det(1)) .eq. 0.0e0) go to 80
   40             if (cabs1(det(1)) .ge. 1.0e0) go to 50
                     det(1) = cmplx(ten,0.0e0)*det(1)
                     det(2) = det(2) - (1.0e0,0.0e0)
                  go to 40
   50             continue
   60             if (cabs1(det(1)) .lt. ten) go to 70
                     det(1) = det(1)/cmplx(ten,0.0e0)
                     det(2) = det(2) + (1.0e0,0.0e0)
                  go to 60
   70             continue
   80          continue
   90       continue
            ik = ik + k
  100    continue
  110 continue
c
c     compute inverse(a)
c
      if (noinv) go to 240
         k = 1
         ik = 0
  120    if (k .gt. n) go to 230
            km1 = k - 1
            kk = ik + k
            ikp1 = ik + k
            if (kpvt(k) .lt. 0) go to 150
c
c              1 by 1
c
               ap(kk) = (1.0e0,0.0e0)/ap(kk)
               if (km1 .lt. 1) go to 140
                  call ccopy(km1,ap(ik+1),1,work,1)
                  ij = 0
                  do 130 j = 1, km1
                     jk = ik + j
                     ap(jk) = cdotu(j,ap(ij+1),1,work,1)
                     call caxpy(j-1,work(j),ap(ij+1),1,ap(ik+1),1)
                     ij = ij + j
  130             continue
                  ap(kk) = ap(kk) + cdotu(km1,work,1,ap(ik+1),1)
  140          continue
               kstep = 1
            go to 190
  150       continue
c
c              2 by 2
c
               kkp1 = ikp1 + k
               t = ap(kkp1)
               ak = ap(kk)/t
               akp1 = ap(kkp1+1)/t
               akkp1 = ap(kkp1)/t
               d = t*(ak*akp1 - (1.0e0,0.0e0))
               ap(kk) = akp1/d
               ap(kkp1+1) = ak/d
               ap(kkp1) = -akkp1/d
               if (km1 .lt. 1) go to 180
                  call ccopy(km1,ap(ikp1+1),1,work,1)
                  ij = 0
                  do 160 j = 1, km1
                     jkp1 = ikp1 + j
                     ap(jkp1) = cdotu(j,ap(ij+1),1,work,1)
                     call caxpy(j-1,work(j),ap(ij+1),1,ap(ikp1+1),1)
                     ij = ij + j
  160             continue
                  ap(kkp1+1) = ap(kkp1+1)
     *                         + cdotu(km1,work,1,ap(ikp1+1),1)
                  ap(kkp1) = ap(kkp1)
     *                       + cdotu(km1,ap(ik+1),1,ap(ikp1+1),1)
                  call ccopy(km1,ap(ik+1),1,work,1)
                  ij = 0
                  do 170 j = 1, km1
                     jk = ik + j
                     ap(jk) = cdotu(j,ap(ij+1),1,work,1)
                     call caxpy(j-1,work(j),ap(ij+1),1,ap(ik+1),1)
                     ij = ij + j
  170             continue
                  ap(kk) = ap(kk) + cdotu(km1,work,1,ap(ik+1),1)
  180          continue
               kstep = 2
  190       continue
c
c           swap
c
            ks = iabs(kpvt(k))
            if (ks .eq. k) go to 220
               iks = (ks*(ks - 1))/2
               call cswap(ks,ap(iks+1),1,ap(ik+1),1)
               ksj = ik + ks
               do 200 jb = ks, k
                  j = k + ks - jb
                  jk = ik + j
                  temp = ap(jk)
                  ap(jk) = ap(ksj)
                  ap(ksj) = temp
                  ksj = ksj - (j - 1)
  200          continue
               if (kstep .eq. 1) go to 210
                  kskp1 = ikp1 + ks
                  temp = ap(kskp1)
                  ap(kskp1) = ap(kkp1)
                  ap(kkp1) = temp
  210          continue
  220       continue
            ik = ik + k
            if (kstep .eq. 2) ik = ik + k + 1
            k = k + kstep
         go to 120
  230    continue
  240 continue
      return
      end
      subroutine cspfa(ap,n,kpvt,info)

c*********************************************************************72

      integer n,kpvt(1),info
      complex ap(1)
c
c     cspfa factors a complex symmetric matrix stored in
c     packed form by elimination with symmetric pivoting.
c
c     to solve  a*x = b , follow cspfa by cspsl.
c     to compute  inverse(a)*c , follow cspfa by cspsl.
c     to compute  determinant(a) , follow cspfa by cspdi.
c     to compute  inverse(a) , follow cspfa by cspdi.
c
c     on entry
c
c        ap      complex (n*(n+1)/2)
c                the packed form of a symmetric matrix  a .  the
c                columns of the upper triangle are stored sequentially
c                in a one-dimensional array of length  n*(n+1)/2 .
c                see comments below for details.
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        ap      a block diagonal matrix and the multipliers which
c                were used to obtain it stored in packed form.
c                the factorization can be written  a = u*d*trans(u)
c                where  u  is a product of permutation and unit
c                upper triangular matrices , trans(u) is the
c                transpose of  u , and  d  is block diagonal
c                with 1 by 1 and 2 by 2 blocks.
c
c        kpvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if the k-th pivot block is singular. this is
c                     not an error condition for this subroutine,
c                     but it does indicate that cspsl or cspdi may
c                     divide by zero if called.
c
c     packed storage
c
c          the following program segment will pack the upper
c          triangle of a symmetric matrix.
c
c                k = 0
c                do 20 j = 1, n
c                   do 10 i = 1, j
c                      k = k + 1
c                      ap(k)  = a(i,j)
c             10    continue
c             20 continue
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab.
c
c     subroutines and functions
c
c     blas caxpy,cswap,icamax
c     fortran abs,aimag,amax1,real,sqrt
c
c     internal variables
c
      complex ak,akm1,bk,bkm1,denom,mulk,mulkm1,t
      real absakk,alpha,colmax,rowmax
      integer icamax,ij,ijj,ik,ikm1,im,imax,imaxp1,imim,imj,imk
      integer j,jj,jk,jkm1,jmax,jmim,k,kk,km1,km1k,km1km1,km2,kstep
      logical swap
c
      complex zdum
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
c
c     initialize
c
c     alpha is used in choosing pivot block size.
      alpha = (1.0e0 + sqrt(17.0e0))/8.0e0
c
      info = 0
c
c     main loop on k, which goes from n to 1.
c
      k = n
      ik = (n*(n - 1))/2
   10 continue
c
c        leave the loop if k=0 or k=1.
c
c     ...exit
         if (k .eq. 0) go to 200
         if (k .gt. 1) go to 20
            kpvt(1) = 1
            if (cabs1(ap(1)) .eq. 0.0e0) info = 1
c     ......exit
            go to 200
   20    continue
c
c        this section of code determines the kind of
c        elimination to be performed.  when it is completed,
c        kstep will be set to the size of the pivot block, and
c        swap will be set to .true. if an interchange is
c        required.
c
         km1 = k - 1
         kk = ik + k
         absakk = cabs1(ap(kk))
c
c        determine the largest off-diagonal element in
c        column k.
c
         imax = icamax(k-1,ap(ik+1),1)
         imk = ik + imax
         colmax = cabs1(ap(imk))
         if (absakk .lt. alpha*colmax) go to 30
            kstep = 1
            swap = .false.
         go to 90
   30    continue
c
c           determine the largest off-diagonal element in
c           row imax.
c
            rowmax = 0.0e0
            imaxp1 = imax + 1
            im = imax*(imax - 1)/2
            imj = im + 2*imax
            do 40 j = imaxp1, k
               rowmax = amax1(rowmax,cabs1(ap(imj)))
               imj = imj + j
   40       continue
            if (imax .eq. 1) go to 50
               jmax = icamax(imax-1,ap(im+1),1)
               jmim = jmax + im
               rowmax = amax1(rowmax,cabs1(ap(jmim)))
   50       continue
            imim = imax + im
            if (cabs1(ap(imim)) .lt. alpha*rowmax) go to 60
               kstep = 1
               swap = .true.
            go to 80
   60       continue
            if (absakk .lt. alpha*colmax*(colmax/rowmax)) go to 70
               kstep = 1
               swap = .false.
            go to 80
   70       continue
               kstep = 2
               swap = imax .ne. km1
   80       continue
   90    continue
         if (amax1(absakk,colmax) .ne. 0.0e0) go to 100
c
c           column k is zero.  set info and iterate the loop.
c
            kpvt(k) = k
            info = k
         go to 190
  100    continue
         if (kstep .eq. 2) go to 140
c
c           1 x 1 pivot block.
c
            if (.not.swap) go to 120
c
c              perform an interchange.
c
               call cswap(imax,ap(im+1),1,ap(ik+1),1)
               imj = ik + imax
               do 110 jj = imax, k
                  j = k + imax - jj
                  jk = ik + j
                  t = ap(jk)
                  ap(jk) = ap(imj)
                  ap(imj) = t
                  imj = imj - (j - 1)
  110          continue
  120       continue
c
c           perform the elimination.
c
            ij = ik - (k - 1)
            do 130 jj = 1, km1
               j = k - jj
               jk = ik + j
               mulk = -ap(jk)/ap(kk)
               t = mulk
               call caxpy(j,t,ap(ik+1),1,ap(ij+1),1)
               ijj = ij + j
               ap(jk) = mulk
               ij = ij - (j - 1)
  130       continue
c
c           set the pivot array.
c
            kpvt(k) = k
            if (swap) kpvt(k) = imax
         go to 190
  140    continue
c
c           2 x 2 pivot block.
c
            km1k = ik + k - 1
            ikm1 = ik - (k - 1)
            if (.not.swap) go to 160
c
c              perform an interchange.
c
               call cswap(imax,ap(im+1),1,ap(ikm1+1),1)
               imj = ikm1 + imax
               do 150 jj = imax, km1
                  j = km1 + imax - jj
                  jkm1 = ikm1 + j
                  t = ap(jkm1)
                  ap(jkm1) = ap(imj)
                  ap(imj) = t
                  imj = imj - (j - 1)
  150          continue
               t = ap(km1k)
               ap(km1k) = ap(imk)
               ap(imk) = t
  160       continue
c
c           perform the elimination.
c
            km2 = k - 2
            if (km2 .eq. 0) go to 180
               ak = ap(kk)/ap(km1k)
               km1km1 = ikm1 + k - 1
               akm1 = ap(km1km1)/ap(km1k)
               denom = 1.0e0 - ak*akm1
               ij = ik - (k - 1) - (k - 2)
               do 170 jj = 1, km2
                  j = km1 - jj
                  jk = ik + j
                  bk = ap(jk)/ap(km1k)
                  jkm1 = ikm1 + j
                  bkm1 = ap(jkm1)/ap(km1k)
                  mulk = (akm1*bk - bkm1)/denom
                  mulkm1 = (ak*bkm1 - bk)/denom
                  t = mulk
                  call caxpy(j,t,ap(ik+1),1,ap(ij+1),1)
                  t = mulkm1
                  call caxpy(j,t,ap(ikm1+1),1,ap(ij+1),1)
                  ap(jk) = mulk
                  ap(jkm1) = mulkm1
                  ijj = ij + j
                  ij = ij - (j - 1)
  170          continue
  180       continue
c
c           set the pivot array.
c
            kpvt(k) = 1 - k
            if (swap) kpvt(k) = -imax
            kpvt(k-1) = kpvt(k)
  190    continue
         ik = ik - (k - 1)
         if (kstep .eq. 2) ik = ik - (k - 2)
         k = k - kstep
      go to 10
  200 continue
      return
      end
      subroutine cspsl(ap,n,kpvt,b)

c*********************************************************************72

      integer n,kpvt(1)
      complex ap(1),b(1)
c
c     csisl solves the complex symmetric system
c     a * x = b
c     using the factors computed by cspfa.
c
c     on entry
c
c        ap      complex(n*(n+1)/2)
c                the output from cspfa.
c
c        n       integer
c                the order of the matrix  a .
c
c        kpvt    integer(n)
c                the pivot vector from cspfa.
c
c        b       complex(n)
c                the right hand side vector.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero may occur if  cspco  has set rcond .eq. 0.0
c        or  cspfa  has set info .ne. 0  .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call cspfa(ap,n,kpvt,info)
c           if (info .ne. 0) go to ...
c           do 10 j = 1, p
c              call cspsl(ap,n,kpvt,c(1,j))
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab.
c
c     subroutines and functions
c
c     blas caxpy,cdotu
c     fortran iabs
c
c     internal variables.
c
      complex ak,akm1,bk,bkm1,cdotu,denom,temp
      integer ik,ikm1,ikp1,k,kk,km1k,km1km1,kp
c
c     loop backward applying the transformations and
c     d inverse to b.
c
      k = n
      ik = (n*(n - 1))/2
   10 if (k .eq. 0) go to 80
         kk = ik + k
         if (kpvt(k) .lt. 0) go to 40
c
c           1 x 1 pivot block.
c
            if (k .eq. 1) go to 30
               kp = kpvt(k)
               if (kp .eq. k) go to 20
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
   20          continue
c
c              apply the transformation.
c
               call caxpy(k-1,b(k),ap(ik+1),1,b(1),1)
   30       continue
c
c           apply d inverse.
c
            b(k) = b(k)/ap(kk)
            k = k - 1
            ik = ik - k
         go to 70
   40    continue
c
c           2 x 2 pivot block.
c
            ikm1 = ik - (k - 1)
            if (k .eq. 2) go to 60
               kp = iabs(kpvt(k))
               if (kp .eq. k - 1) go to 50
c
c                 interchange.
c
                  temp = b(k-1)
                  b(k-1) = b(kp)
                  b(kp) = temp
   50          continue
c
c              apply the transformation.
c
               call caxpy(k-2,b(k),ap(ik+1),1,b(1),1)
               call caxpy(k-2,b(k-1),ap(ikm1+1),1,b(1),1)
   60       continue
c
c           apply d inverse.
c
            km1k = ik + k - 1
            kk = ik + k
            ak = ap(kk)/ap(km1k)
            km1km1 = ikm1 + k - 1
            akm1 = ap(km1km1)/ap(km1k)
            bk = b(k)/ap(km1k)
            bkm1 = b(k-1)/ap(km1k)
            denom = ak*akm1 - 1.0e0
            b(k) = (akm1*bk - bkm1)/denom
            b(k-1) = (ak*bkm1 - bk)/denom
            k = k - 2
            ik = ik - (k + 1) - k
   70    continue
      go to 10
   80 continue
c
c     loop forward applying the transformations.
c
      k = 1
      ik = 0
   90 if (k .gt. n) go to 160
         if (kpvt(k) .lt. 0) go to 120
c
c           1 x 1 pivot block.
c
            if (k .eq. 1) go to 110
c
c              apply the transformation.
c
               b(k) = b(k) + cdotu(k-1,ap(ik+1),1,b(1),1)
               kp = kpvt(k)
               if (kp .eq. k) go to 100
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
  100          continue
  110       continue
            ik = ik + k
            k = k + 1
         go to 150
  120    continue
c
c           2 x 2 pivot block.
c
            if (k .eq. 1) go to 140
c
c              apply the transformation.
c
               b(k) = b(k) + cdotu(k-1,ap(ik+1),1,b(1),1)
               ikp1 = ik + k
               b(k+1) = b(k+1) + cdotu(k-1,ap(ikp1+1),1,b(1),1)
               kp = iabs(kpvt(k))
               if (kp .eq. k) go to 130
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
  130          continue
  140       continue
            ik = ik + k + k + 1
            k = k + 2
  150    continue
      go to 90
  160 continue
      return
      end
      subroutine csvdc(x,ldx,n,p,s,e,u,ldu,v,ldv,work,job,info)

c*********************************************************************72

      integer ldx,n,p,ldu,ldv,job,info
      complex x(ldx,1),s(1),e(1),u(ldu,1),v(ldv,1),work(1)
c
c
c     csvdc is a subroutine to reduce a complex nxp matrix x by
c     unitary transformations u and v to diagonal form.  the
c     diagonal elements s(i) are the singular values of x.  the
c     columns of u are the corresponding left singular vectors,
c     and the columns of v the right singular vectors.
c
c     on entry
c
c         x         complex(ldx,p), where ldx.ge.n.
c                   x contains the matrix whose singular value
c                   decomposition is to be computed.  x is
c                   destroyed by csvdc.
c
c         ldx       integer.
c                   ldx is the leading dimension of the array x.
c
c         n         integer.
c                   n is the number of rows of the matrix x.
c
c         p         integer.
c                   p is the number of columns of the matrix x.
c
c         ldu       integer.
c                   ldu is the leading dimension of the array u
c                   (see below).
c
c         ldv       integer.
c                   ldv is the leading dimension of the array v
c                   (see below).
c
c         work      complex(n).
c                   work is a scratch array.
c
c         job       integer.
c                   job controls the computation of the singular
c                   vectors.  it has the decimal expansion ab
c                   with the following meaning
c
c                        a.eq.0    do not compute the left singular
c                                  vectors.
c                        a.eq.1    return the n left singular vectors
c                                  in u.
c                        a.ge.2    returns the first min(n,p)
c                                  left singular vectors in u.
c                        b.eq.0    do not compute the right singular
c                                  vectors.
c                        b.eq.1    return the right singular vectors
c                                  in v.
c
c     on return
c
c         s         complex(mm), where mm=max(n+1,p).
c                   the first min(n,p) entries of s contain the
c                   singular values of x arranged in descending
c                   order of magnitude.
c
c         e         complex(mm), where mm=max(n+1,p).
c                   e ordinarily contains zeros.  however see the
c                   discussion of info for exceptions.
c
c         u         complex(ldu,k), where ldu.ge.n.  if joba.eq.1 then
c                                   k.eq.n, if joba.ge.2 then
c                                   k.eq.min(n,p).
c                   u contains the matrix of left singular vectors.
c                   u is not referenced if joba.eq.0.  if n.le.p
c                   or if joba.gt.2, then u may be identified with x
c                   in the subroutine call.
c
c         v         complex(ldv,p), where ldv.ge.p.
c                   v contains the matrix of right singular vectors.
c                   v is not referenced if jobb.eq.0.  if p.le.n,
c                   then v may be identified whth x in the
c                   subroutine call.
c
c         info      integer.
c                   the singular values (and their corresponding
c                   singular vectors) s(info+1),s(info+2),...,s(m)
c                   are correct (here m=min(n,p)).  thus if
c                   info.eq.0, all the singular values and their
c                   vectors are correct.  in any event, the matrix
c                   b = ctrans(u)*x*v is the bidiagonal matrix
c                   with the elements of s on its diagonal and the
c                   elements of e on its super-diagonal (ctrans(u)
c                   is the conjugate-transpose of u).  thus the
c                   singular values of x and b are the same.
c
c     linpack. this version dated 03/19/79 .
c              correction to shift calculation made 2/85.
c     g.w. stewart, university of maryland, argonne national lab.
c
c     csvdc uses the following functions and subprograms.
c
c     external csrot
c     blas caxpy,cdotc,cscal,cswap,scnrm2,srotg
c     fortran abs,aimag,amax1,cabs,cmplx
c     fortran conjg,max0,min0,mod,real,sqrt
c
c     internal variables
c
      integer i,iter,j,jobu,k,kase,kk,l,ll,lls,lm1,lp1,ls,lu,m,maxit,
     *        mm,mm1,mp1,nct,nctp1,ncu,nrt,nrtp1
      complex cdotc,t,r
      real b,c,cs,el,emm1,f,g,scnrm2,scale,shift,sl,sm,sn,smm1,t1,test,
     *     ztest
      logical wantu,wantv
c
      complex csign,zdum,zdum1,zdum2
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
      csign(zdum1,zdum2) = cabs(zdum1)*(zdum2/cabs(zdum2))
c
c     set the maximum number of iterations.
c
      maxit = 30
c
c     determine what is to be computed.
c
      wantu = .false.
      wantv = .false.
      jobu = mod(job,100)/10
      ncu = n
      if (jobu .gt. 1) ncu = min0(n,p)
      if (jobu .ne. 0) wantu = .true.
      if (mod(job,10) .ne. 0) wantv = .true.
c
c     reduce x to bidiagonal form, storing the diagonal elements
c     in s and the super-diagonal elements in e.
c
      info = 0
      nct = min0(n-1,p)
      nrt = max0(0,min0(p-2,n))
      lu = max0(nct,nrt)
      if (lu .lt. 1) go to 170
      do 160 l = 1, lu
         lp1 = l + 1
         if (l .gt. nct) go to 20
c
c           compute the transformation for the l-th column and
c           place the l-th diagonal in s(l).
c
            s(l) = cmplx(scnrm2(n-l+1,x(l,l),1),0.0e0)
            if (cabs1(s(l)) .eq. 0.0e0) go to 10
               if (cabs1(x(l,l)) .ne. 0.0e0) s(l) = csign(s(l),x(l,l))
               call cscal(n-l+1,1.0e0/s(l),x(l,l),1)
               x(l,l) = (1.0e0,0.0e0) + x(l,l)
   10       continue
            s(l) = -s(l)
   20    continue
         if (p .lt. lp1) go to 50
         do 40 j = lp1, p
            if (l .gt. nct) go to 30
            if (cabs1(s(l)) .eq. 0.0e0) go to 30
c
c              apply the transformation.
c
               t = -cdotc(n-l+1,x(l,l),1,x(l,j),1)/x(l,l)
               call caxpy(n-l+1,t,x(l,l),1,x(l,j),1)
   30       continue
c
c           place the l-th row of x into  e for the
c           subsequent calculation of the row transformation.
c
            e(j) = conjg(x(l,j))
   40    continue
   50    continue
         if (.not.wantu .or. l .gt. nct) go to 70
c
c           place the transformation in u for subsequent back
c           multiplication.
c
            do 60 i = l, n
               u(i,l) = x(i,l)
   60       continue
   70    continue
         if (l .gt. nrt) go to 150
c
c           compute the l-th row transformation and place the
c           l-th super-diagonal in e(l).
c
            e(l) = cmplx(scnrm2(p-l,e(lp1),1),0.0e0)
            if (cabs1(e(l)) .eq. 0.0e0) go to 80
               if (cabs1(e(lp1)) .ne. 0.0e0) e(l) = csign(e(l),e(lp1))
               call cscal(p-l,1.0e0/e(l),e(lp1),1)
               e(lp1) = (1.0e0,0.0e0) + e(lp1)
   80       continue
            e(l) = -conjg(e(l))
            if (lp1 .gt. n .or. cabs1(e(l)) .eq. 0.0e0) go to 120
c
c              apply the transformation.
c
               do 90 i = lp1, n
                  work(i) = (0.0e0,0.0e0)
   90          continue
               do 100 j = lp1, p
                  call caxpy(n-l,e(j),x(lp1,j),1,work(lp1),1)
  100          continue
               do 110 j = lp1, p
                  call caxpy(n-l,conjg(-e(j)/e(lp1)),work(lp1),1,
     *                       x(lp1,j),1)
  110          continue
  120       continue
            if (.not.wantv) go to 140
c
c              place the transformation in v for subsequent
c              back multiplication.
c
               do 130 i = lp1, p
                  v(i,l) = e(i)
  130          continue
  140       continue
  150    continue
  160 continue
  170 continue
c
c     set up the final bidiagonal matrix or order m.
c
      m = min0(p,n+1)
      nctp1 = nct + 1
      nrtp1 = nrt + 1
      if (nct .lt. p) s(nctp1) = x(nctp1,nctp1)
      if (n .lt. m) s(m) = (0.0e0,0.0e0)
      if (nrtp1 .lt. m) e(nrtp1) = x(nrtp1,m)
      e(m) = (0.0e0,0.0e0)
c
c     if required, generate u.
c
      if (.not.wantu) go to 300
         if (ncu .lt. nctp1) go to 200
         do 190 j = nctp1, ncu
            do 180 i = 1, n
               u(i,j) = (0.0e0,0.0e0)
  180       continue
            u(j,j) = (1.0e0,0.0e0)
  190    continue
  200    continue
         if (nct .lt. 1) go to 290
         do 280 ll = 1, nct
            l = nct - ll + 1
            if (cabs1(s(l)) .eq. 0.0e0) go to 250
               lp1 = l + 1
               if (ncu .lt. lp1) go to 220
               do 210 j = lp1, ncu
                  t = -cdotc(n-l+1,u(l,l),1,u(l,j),1)/u(l,l)
                  call caxpy(n-l+1,t,u(l,l),1,u(l,j),1)
  210          continue
  220          continue
               call cscal(n-l+1,(-1.0e0,0.0e0),u(l,l),1)
               u(l,l) = (1.0e0,0.0e0) + u(l,l)
               lm1 = l - 1
               if (lm1 .lt. 1) go to 240
               do 230 i = 1, lm1
                  u(i,l) = (0.0e0,0.0e0)
  230          continue
  240          continue
            go to 270
  250       continue
               do 260 i = 1, n
                  u(i,l) = (0.0e0,0.0e0)
  260          continue
               u(l,l) = (1.0e0,0.0e0)
  270       continue
  280    continue
  290    continue
  300 continue
c
c     if it is required, generate v.
c
      if (.not.wantv) go to 350
         do 340 ll = 1, p
            l = p - ll + 1
            lp1 = l + 1
            if (l .gt. nrt) go to 320
            if (cabs1(e(l)) .eq. 0.0e0) go to 320
               do 310 j = lp1, p
                  t = -cdotc(p-l,v(lp1,l),1,v(lp1,j),1)/v(lp1,l)
                  call caxpy(p-l,t,v(lp1,l),1,v(lp1,j),1)
  310          continue
  320       continue
            do 330 i = 1, p
               v(i,l) = (0.0e0,0.0e0)
  330       continue
            v(l,l) = (1.0e0,0.0e0)
  340    continue
  350 continue
c
c     transform s and e so that they are real.
c
      do 380 i = 1, m
         if (cabs1(s(i)) .eq. 0.0e0) go to 360
            t = cmplx(cabs(s(i)),0.0e0)
            r = s(i)/t
            s(i) = t
            if (i .lt. m) e(i) = e(i)/r
            if (wantu) call cscal(n,r,u(1,i),1)
  360    continue
c     ...exit
         if (i .eq. m) go to 390
         if (cabs1(e(i)) .eq. 0.0e0) go to 370
            t = cmplx(cabs(e(i)),0.0e0)
            r = t/e(i)
            e(i) = t
            s(i+1) = s(i+1)*r
            if (wantv) call cscal(p,r,v(1,i+1),1)
  370    continue
  380 continue
  390 continue
c
c     main iteration loop for the singular values.
c
      mm = m
      iter = 0
  400 continue
c
c        quit if all the singular values have been found.
c
c     ...exit
         if (m .eq. 0) go to 660
c
c        if too many iterations have been performed, set
c        flag and return.
c
         if (iter .lt. maxit) go to 410
            info = m
c     ......exit
            go to 660
  410    continue
c
c        this section of the program inspects for
c        negligible elements in the s and e arrays.  on
c        completion the variables kase and l are set as follows.
c
c           kase = 1     if s(m) and e(l-1) are negligible and l.lt.m
c           kase = 2     if s(l) is negligible and l.lt.m
c           kase = 3     if e(l-1) is negligible, l.lt.m, and
c                        s(l), ..., s(m) are not negligible (qr step).
c           kase = 4     if e(m-1) is negligible (convergence).
c
         do 430 ll = 1, m
            l = m - ll
c        ...exit
            if (l .eq. 0) go to 440
            test = cabs(s(l)) + cabs(s(l+1))
            ztest = test + cabs(e(l))
            if (ztest .ne. test) go to 420
               e(l) = (0.0e0,0.0e0)
c        ......exit
               go to 440
  420       continue
  430    continue
  440    continue
         if (l .ne. m - 1) go to 450
            kase = 4
         go to 520
  450    continue
            lp1 = l + 1
            mp1 = m + 1
            do 470 lls = lp1, mp1
               ls = m - lls + lp1
c           ...exit
               if (ls .eq. l) go to 480
               test = 0.0e0
               if (ls .ne. m) test = test + cabs(e(ls))
               if (ls .ne. l + 1) test = test + cabs(e(ls-1))
               ztest = test + cabs(s(ls))
               if (ztest .ne. test) go to 460
                  s(ls) = (0.0e0,0.0e0)
c           ......exit
                  go to 480
  460          continue
  470       continue
  480       continue
            if (ls .ne. l) go to 490
               kase = 3
            go to 510
  490       continue
            if (ls .ne. m) go to 500
               kase = 1
            go to 510
  500       continue
               kase = 2
               l = ls
  510       continue
  520    continue
         l = l + 1
c
c        perform the task indicated by kase.
c
         go to (530, 560, 580, 610), kase
c
c        deflate negligible s(m).
c
  530    continue
            mm1 = m - 1
            f = real(e(m-1))
            e(m-1) = (0.0e0,0.0e0)
            do 550 kk = l, mm1
               k = mm1 - kk + l
               t1 = real(s(k))
               call srotg(t1,f,cs,sn)
               s(k) = cmplx(t1,0.0e0)
               if (k .eq. l) go to 540
                  f = -sn*real(e(k-1))
                  e(k-1) = cs*e(k-1)
  540          continue
               if (wantv) call csrot(p,v(1,k),1,v(1,m),1,cs,sn)
  550       continue
         go to 650
c
c        split at negligible s(l).
c
  560    continue
            f = real(e(l-1))
            e(l-1) = (0.0e0,0.0e0)
            do 570 k = l, m
               t1 = real(s(k))
               call srotg(t1,f,cs,sn)
               s(k) = cmplx(t1,0.0e0)
               f = -sn*real(e(k))
               e(k) = cs*e(k)
               if (wantu) call csrot(n,u(1,k),1,u(1,l-1),1,cs,sn)
  570       continue
         go to 650
c
c        perform one qr step.
c
  580    continue
c
c           calculate the shift.
c
            scale = amax1(cabs(s(m)),cabs(s(m-1)),cabs(e(m-1)),
     *                    cabs(s(l)),cabs(e(l)))
            sm = real(s(m))/scale
            smm1 = real(s(m-1))/scale
            emm1 = real(e(m-1))/scale
            sl = real(s(l))/scale
            el = real(e(l))/scale
            b = ((smm1 + sm)*(smm1 - sm) + emm1**2)/2.0e0
            c = (sm*emm1)**2
            shift = 0.0e0
            if (b .eq. 0.0e0 .and. c .eq. 0.0e0) go to 590
               shift = sqrt(b**2+c)
               if (b .lt. 0.0e0) shift = -shift
               shift = c/(b + shift)
  590       continue
            f = (sl + sm)*(sl - sm) + shift
            g = sl*el
c
c           chase zeros.
c
            mm1 = m - 1
            do 600 k = l, mm1
               call srotg(f,g,cs,sn)
               if (k .ne. l) e(k-1) = cmplx(f,0.0e0)
               f = cs*real(s(k)) + sn*real(e(k))
               e(k) = cs*e(k) - sn*s(k)
               g = sn*real(s(k+1))
               s(k+1) = cs*s(k+1)
               if (wantv) call csrot(p,v(1,k),1,v(1,k+1),1,cs,sn)
               call srotg(f,g,cs,sn)
               s(k) = cmplx(f,0.0e0)
               f = cs*real(e(k)) + sn*real(s(k+1))
               s(k+1) = -sn*e(k) + cs*s(k+1)
               g = sn*real(e(k+1))
               e(k+1) = cs*e(k+1)
               if (wantu .and. k .lt. n)
     *            call csrot(n,u(1,k),1,u(1,k+1),1,cs,sn)
  600       continue
            e(m-1) = cmplx(f,0.0e0)
            iter = iter + 1
         go to 650
c
c        convergence.
c
  610    continue
c
c           make the singular value  positive
c
            if (real(s(l)) .ge. 0.0e0) go to 620
               s(l) = -s(l)
               if (wantv) call cscal(p,(-1.0e0,0.0e0),v(1,l),1)
  620       continue
c
c           order the singular value.
c
  630       if (l .eq. mm) go to 640
c           ...exit
               if (real(s(l)) .ge. real(s(l+1))) go to 640
               t = s(l)
               s(l) = s(l+1)
               s(l+1) = t
               if (wantv .and. l .lt. p)
     *            call cswap(p,v(1,l),1,v(1,l+1),1)
               if (wantu .and. l .lt. n)
     *            call cswap(n,u(1,l),1,u(1,l+1),1)
               l = l + 1
            go to 630
  640       continue
            iter = 0
            m = m - 1
  650    continue
      go to 400
  660 continue
      return
      end
      subroutine ctrco(t,ldt,n,rcond,z,job)

c*********************************************************************72

      integer ldt,n,job
      complex t(ldt,1),z(1)
      real rcond
c
c     ctrco estimates the condition of a complex triangular matrix.
c
c     on entry
c
c        t       complex(ldt,n)
c                t contains the triangular matrix. the zero
c                elements of the matrix are not referenced, and
c                the corresponding elements of the array can be
c                used to store other information.
c
c        ldt     integer
c                ldt is the leading dimension of the array t.
c
c        n       integer
c                n is the order of the system.
c
c        job     integer
c                = 0         t  is lower triangular.
c                = nonzero   t  is upper triangular.
c
c     on return
c
c        rcond   real
c                an estimate of the reciprocal condition of  t .
c                for the system  t*x = b , relative perturbations
c                in  t  and  b  of size  epsilon  may cause
c                relative perturbations in  x  of size  epsilon/rcond .
c                if  rcond  is so small that the logical expression
c                           1.0 + rcond .eq. 1.0
c                is true, then  t  may be singular to working
c                precision.  in particular,  rcond  is zero  if
c                exact singularity is detected or the estimate
c                underflows.
c
c        z       complex(n)
c                a work vector whose contents are usually unimportant.
c                if  t  is close to a singular matrix, then  z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas caxpy,csscal,scasum
c     fortran abs,aimag,amax1,cmplx,conjg,real
c
c     internal variables
c
      complex w,wk,wkm,ek
      real tnorm,ynorm,s,sm,scasum
      integer i1,j,j1,j2,k,kk,l
      logical lower
      complex zdum,zdum1,zdum2,csign1
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
      csign1(zdum1,zdum2) = cabs1(zdum1)*(zdum2/cabs1(zdum2))
c
      lower = job .eq. 0
c
c     compute 1-norm of t
c
      tnorm = 0.0e0
      do 10 j = 1, n
         l = j
         if (lower) l = n + 1 - j
         i1 = 1
         if (lower) i1 = j
         tnorm = amax1(tnorm,scasum(l,t(i1,j),1))
   10 continue
c
c     rcond = 1/(norm(t)*(estimate of norm(inverse(t)))) .
c     estimate = norm(z)/norm(y) where  t*z = y  and  ctrans(t)*y = e .
c     ctrans(t)  is the conjugate transpose of t .
c     the components of  e  are chosen to cause maximum local
c     growth in the elements of y .
c     the vectors are frequently rescaled to avoid overflow.
c
c     solve ctrans(t)*y = e
c
      ek = (1.0e0,0.0e0)
      do 20 j = 1, n
         z(j) = (0.0e0,0.0e0)
   20 continue
      do 100 kk = 1, n
         k = kk
         if (lower) k = n + 1 - kk
         if (cabs1(z(k)) .ne. 0.0e0) ek = csign1(ek,-z(k))
         if (cabs1(ek-z(k)) .le. cabs1(t(k,k))) go to 30
            s = cabs1(t(k,k))/cabs1(ek-z(k))
            call csscal(n,s,z,1)
            ek = cmplx(s,0.0e0)*ek
   30    continue
         wk = ek - z(k)
         wkm = -ek - z(k)
         s = cabs1(wk)
         sm = cabs1(wkm)
         if (cabs1(t(k,k)) .eq. 0.0e0) go to 40
            wk = wk/conjg(t(k,k))
            wkm = wkm/conjg(t(k,k))
         go to 50
   40    continue
            wk = (1.0e0,0.0e0)
            wkm = (1.0e0,0.0e0)
   50    continue
         if (kk .eq. n) go to 90
            j1 = k + 1
            if (lower) j1 = 1
            j2 = n
            if (lower) j2 = k - 1
            do 60 j = j1, j2
               sm = sm + cabs1(z(j)+wkm*conjg(t(k,j)))
               z(j) = z(j) + wk*conjg(t(k,j))
               s = s + cabs1(z(j))
   60       continue
            if (s .ge. sm) go to 80
               w = wkm - wk
               wk = wkm
               do 70 j = j1, j2
                  z(j) = z(j) + w*conjg(t(k,j))
   70          continue
   80       continue
   90    continue
         z(k) = wk
  100 continue
      s = 1.0e0/scasum(n,z,1)
      call csscal(n,s,z,1)
c
      ynorm = 1.0e0
c
c     solve t*z = y
c
      do 130 kk = 1, n
         k = n + 1 - kk
         if (lower) k = kk
         if (cabs1(z(k)) .le. cabs1(t(k,k))) go to 110
            s = cabs1(t(k,k))/cabs1(z(k))
            call csscal(n,s,z,1)
            ynorm = s*ynorm
  110    continue
         if (cabs1(t(k,k)) .ne. 0.0e0) z(k) = z(k)/t(k,k)
         if (cabs1(t(k,k)) .eq. 0.0e0) z(k) = (1.0e0,0.0e0)
         i1 = 1
         if (lower) i1 = k + 1
         if (kk .ge. n) go to 120
            w = -z(k)
            call caxpy(n-kk,w,t(i1,k),1,z(i1),1)
  120    continue
  130 continue
c     make znorm = 1.0
      s = 1.0e0/scasum(n,z,1)
      call csscal(n,s,z,1)
      ynorm = s*ynorm
c
      if (tnorm .ne. 0.0e0) rcond = ynorm/tnorm
      if (tnorm .eq. 0.0e0) rcond = 0.0e0
      return
      end
      subroutine ctrdi(t,ldt,n,det,job,info)

c*********************************************************************72

      integer ldt,n,job,info
      complex t(ldt,1),det(2)
c
c     ctrdi computes the determinant and inverse of a complex
c     triangular matrix.
c
c     on entry
c
c        t       complex(ldt,n)
c                t contains the triangular matrix. the zero
c                elements of the matrix are not referenced, and
c                the corresponding elements of the array can be
c                used to store other information.
c
c        ldt     integer
c                ldt is the leading dimension of the array t.
c
c        n       integer
c                n is the order of the system.
c
c        job     integer
c                = 010       no det, inverse of lower triangular.
c                = 011       no det, inverse of upper triangular.
c                = 100       det, no inverse.
c                = 110       det, inverse of lower triangular.
c                = 111       det, inverse of upper triangular.
c
c     on return
c
c        t       inverse of original matrix if requested.
c                otherwise unchanged.
c
c        det     complex(2)
c                determinant of original matrix if requested.
c                otherwise not referenced.
c                determinant = det(1) * 10.0**det(2)
c                with  1.0 .le. cabs1(det(1)) .lt. 10.0
c                or  det(1) .eq. 0.0 .
c
c        info    integer
c                info contains zero if the system is nonsingular
c                and the inverse is requested.
c                otherwise info contains the index of
c                a zero diagonal element of t.
c
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas caxpy,cscal
c     fortran abs,aimag,cmplx,mod,real
c
c     internal variables
c
      complex temp
      real ten
      integer i,j,k,kb,km1,kp1
c
      complex zdum
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
c     begin block permitting ...exits to 180
c
c        compute determinant
c
         if (job/100 .eq. 0) go to 70
            det(1) = (1.0e0,0.0e0)
            det(2) = (0.0e0,0.0e0)
            ten = 10.0e0
            do 50 i = 1, n
               det(1) = t(i,i)*det(1)
c           ...exit
               if (cabs1(det(1)) .eq. 0.0e0) go to 60
   10          if (cabs1(det(1)) .ge. 1.0e0) go to 20
                  det(1) = cmplx(ten,0.0e0)*det(1)
                  det(2) = det(2) - (1.0e0,0.0e0)
               go to 10
   20          continue
   30          if (cabs1(det(1)) .lt. ten) go to 40
                  det(1) = det(1)/cmplx(ten,0.0e0)
                  det(2) = det(2) + (1.0e0,0.0e0)
               go to 30
   40          continue
   50       continue
   60       continue
   70    continue
c
c        compute inverse of upper triangular
c
         if (mod(job/10,10) .eq. 0) go to 170
            if (mod(job,10) .eq. 0) go to 120
c              begin block permitting ...exits to 110
                  do 100 k = 1, n
                     info = k
c              ......exit
                     if (cabs1(t(k,k)) .eq. 0.0e0) go to 110
                     t(k,k) = (1.0e0,0.0e0)/t(k,k)
                     temp = -t(k,k)
                     call cscal(k-1,temp,t(1,k),1)
                     kp1 = k + 1
                     if (n .lt. kp1) go to 90
                     do 80 j = kp1, n
                        temp = t(k,j)
                        t(k,j) = (0.0e0,0.0e0)
                        call caxpy(k,temp,t(1,k),1,t(1,j),1)
   80                continue
   90                continue
  100             continue
                  info = 0
  110          continue
            go to 160
  120       continue
c
c              compute inverse of lower triangular
c
               do 150 kb = 1, n
                  k = n + 1 - kb
                  info = k
c     ............exit
                  if (cabs1(t(k,k)) .eq. 0.0e0) go to 180
                  t(k,k) = (1.0e0,0.0e0)/t(k,k)
                  temp = -t(k,k)
                  if (k .ne. n) call cscal(n-k,temp,t(k+1,k),1)
                  km1 = k - 1
                  if (km1 .lt. 1) go to 140
                  do 130 j = 1, km1
                     temp = t(k,j)
                     t(k,j) = (0.0e0,0.0e0)
                     call caxpy(n-k+1,temp,t(k,k),1,t(k,j),1)
  130             continue
  140             continue
  150          continue
               info = 0
  160       continue
  170    continue
  180 continue
      return
      end
      subroutine ctrsl(t,ldt,n,b,job,info)

c*********************************************************************72

      integer ldt,n,job,info
      complex t(ldt,1),b(1)
c
c
c     ctrsl solves systems of the form
c
c                   t * x = b
c     or
c                   ctrans(t) * x = b
c
c     where t is a triangular matrix of order n. here ctrans(t)
c     denotes the conjugate transpose of the matrix t.
c
c     on entry
c
c         t         complex(ldt,n)
c                   t contains the matrix of the system. the zero
c                   elements of the matrix are not referenced, and
c                   the corresponding elements of the array can be
c                   used to store other information.
c
c         ldt       integer
c                   ldt is the leading dimension of the array t.
c
c         n         integer
c                   n is the order of the system.
c
c         b         complex(n).
c                   b contains the right hand side of the system.
c
c         job       integer
c                   job specifies what kind of system is to be solved.
c                   if job is
c
c                        00   solve t*x=b, t lower triangular,
c                        01   solve t*x=b, t upper triangular,
c                        10   solve ctrans(t)*x=b, t lower triangular,
c                        11   solve ctrans(t)*x=b, t upper triangular.
c
c     on return
c
c         b         b contains the solution, if info .eq. 0.
c                   otherwise b is unaltered.
c
c         info      integer
c                   info contains zero if the system is nonsingular.
c                   otherwise info contains the index of
c                   the first zero diagonal element of t.
c
c     linpack. this version dated 08/14/78 .
c     g. w. stewart, university of maryland, argonne national lab.
c
c     subroutines and functions
c
c     blas caxpy,cdotc
c     fortran abs,aimag,conjg,mod,real
c
c     internal variables
c
      complex cdotc,temp
      integer case,j,jj
c
      complex zdum
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
c     begin block permitting ...exits to 150
c
c        check for zero diagonal elements.
c
         do 10 info = 1, n
c     ......exit
            if (cabs1(t(info,info)) .eq. 0.0e0) go to 150
   10    continue
         info = 0
c
c        determine the task and go to it.
c
         case = 1
         if (mod(job,10) .ne. 0) case = 2
         if (mod(job,100)/10 .ne. 0) case = case + 2
         go to (20,50,80,110), case
c
c        solve t*x=b for t lower triangular
c
   20    continue
            b(1) = b(1)/t(1,1)
            if (n .lt. 2) go to 40
            do 30 j = 2, n
               temp = -b(j-1)
               call caxpy(n-j+1,temp,t(j,j-1),1,b(j),1)
               b(j) = b(j)/t(j,j)
   30       continue
   40       continue
         go to 140
c
c        solve t*x=b for t upper triangular.
c
   50    continue
            b(n) = b(n)/t(n,n)
            if (n .lt. 2) go to 70
            do 60 jj = 2, n
               j = n - jj + 1
               temp = -b(j+1)
               call caxpy(j,temp,t(1,j+1),1,b(1),1)
               b(j) = b(j)/t(j,j)
   60       continue
   70       continue
         go to 140
c
c        solve ctrans(t)*x=b for t lower triangular.
c
   80    continue
            b(n) = b(n)/conjg(t(n,n))
            if (n .lt. 2) go to 100
            do 90 jj = 2, n
               j = n - jj + 1
               b(j) = b(j) - cdotc(jj-1,t(j+1,j),1,b(j+1),1)
               b(j) = b(j)/conjg(t(j,j))
   90       continue
  100       continue
         go to 140
c
c        solve ctrans(t)*x=b for t upper triangular.
c
  110    continue
            b(1) = b(1)/conjg(t(1,1))
            if (n .lt. 2) go to 130
            do 120 j = 2, n
               b(j) = b(j) - cdotc(j-1,t(1,j),1,b(1),1)
               b(j) = b(j)/conjg(t(j,j))
  120       continue
  130       continue
  140    continue
  150 continue
      return
      end
