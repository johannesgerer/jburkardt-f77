      subroutine banfac ( w, nroww, nrow, nbandl, nbandu, iflag )

c*********************************************************************72
c
cc BANFAC factors a banded matrix without pivoting.
c
c  from  * a practical guide to splines *  by c. de boor    
c  returns in  w  the lu-factorization (without pivoting) of the banded
c  matrix  a  of order  nrow  with  (nbandl + 1 + nbandu) bands or diag-
c  onals in the work array  w .
c
c******  i n p u t  ******
c  w.....work array of size  (nroww,nrow)  containing the interesting
c        part of a banded matrix  a , with the diagonals or bands of  a
c        stored in the rows of  w , while columns of  a  correspond to
c        columns of  w . this is the storage mode used in  linpack  and
c        results in efficient innermost loops.
c           explicitly,  a  has  nbandl  bands below the diagonal
c                            +     1     (main) diagonal
c                            +   nbandu  bands above the diagonal
c        and thus, with    middle = nbandu + 1,
c          a(i+j,j)  is in  w(i+middle,j)  for i=-nbandu,...,nbandl
c                                              j=1,...,nrow .
c        for example, the interesting entries of a (1,2)-banded matrix
c        of order  9  would appear in the first  1+1+2 = 4  rows of  w
c        as follows.
c
c
c
c
c
c        all other entries of  w  not identified in this way with an en-
c        try of  a  are never referenced .
c  nroww.....row dimension of the work array  w .
c        must be  .ge.  nbandl + 1 + nbandu  .
c  nbandl.....number of bands of  a  below the main diagonal
c  nbandu.....number of bands of  a  above the main diagonal .
c
c******  o u t p u t  ******
c  iflag.....integer indicating success( = 1) or failure ( = 2) .
c     if  iflag = 1, then
c  w.....contains the lu-factorization of  a  into a unit lower triangu-
c        lar matrix  l  and an upper triangular matrix  u (both banded)
c        and stored in customary fashion over the corresponding entries
c        of  a . this makes it possible to solve any particular linear
c        system  a*x = b  for  x  by a
c              call banslv ( w, nroww, nrow, nbandl, nbandu, b )
c        with the solution x  contained in  b  on return .
c     if  iflag = 2, then
c        one of  nrow-1, nbandl,nbandu failed to be nonnegative, or else
c        one of the potential pivots was found to be zero indicating
c        that  a  does not have an lu-factorization. this implies that
c        a  is singular in case it is totally positive .
c
c******  m e t h o d  ******
c     gauss elimination  w i t h o u t  pivoting is used. the routine is
c  intended for use with matrices  a  which do not require row inter-
c  changes during factorization, especially for the  t o t a l l y
c  p o s i t i v e  matrices which occur in spline calculations.
c     the routine should not be used for an arbitrary banded matrix.
c
      implicit none

      integer iflag,nbandl,nbandu,nrow,nroww,   i,ipk,j,jmax,k,kmax
     &                                        ,middle,midmk,nrowm1
      double precision w(nroww,nrow),   factor,pivot
c
      iflag = 1
      middle = nbandu + 1
c                         w(middle,.) contains the main diagonal of  a .
      nrowm1 = nrow - 1
      if (nrowm1)                       999,900,1
    1 if (nbandl .gt. 0)                go to 10
c                a is upper triangular. check that diagonal is nonzero .
      do 5 i=1,nrowm1
         if (w(middle,i) .eq. 0.0D+00)       go to 999
    5    continue
                                        go to 900
   10 if (nbandu .gt. 0)                go to 20
c              a is lower triangular. check that diagonal is nonzero and
c                 divide each column by its diagonal .
      do 15 i=1,nrowm1
         pivot = w(middle,i)
         if(pivot .eq. 0.0D+00)              go to 999
         jmax = min0(nbandl, nrow - i)
         do 15 j=1,jmax
   15       w(middle+j,i) = w(middle+j,i)/pivot
                                        return
c
c        a  is not just a triangular matrix. construct lu factorization
   20 do 50 i=1,nrowm1
c                                  w(middle,i)  is pivot for i-th step .
         pivot = w(middle,i)
         if (pivot .eq. 0.0D+00)             go to 999
c                 jmax  is the number of (nonzero) entries in column  i
c                     below the diagonal .
         jmax = min0(nbandl,nrow - i)
c              divide each entry in column  i  below diagonal by pivot .
         do 32 j=1,jmax
   32       w(middle+j,i) = w(middle+j,i)/pivot
c                 kmax  is the number of (nonzero) entries in row  i  to
c                     the right of the diagonal .
         kmax = min0(nbandu,nrow - i)
c                  subtract  a(i,i+k)*(i-th column) from (i+k)-th column
c                  (below row  i ) .
         do 40 k=1,kmax
            ipk = i + k
            midmk = middle - k
            factor = w(midmk,ipk)
            do 40 j=1,jmax
   40          w(midmk+j,ipk) = w(midmk+j,ipk) - w(middle+j,i)*factor
   50    continue
c                                       check the last diagonal entry .
  900 if (w(middle,nrow) .ne. 0.0D+00)       return
  999 iflag = 2
                                        return
      end
      subroutine banslv ( w, nroww, nrow, nbandl, nbandu, b )

c*********************************************************************72
c
cc BANSLV solves a banded linear system A * X = B factored by BANFAC.
c
c  from  * a practical guide to splines *  by c. de boor    
c  companion routine to  banfac . it returns the solution  x  of the
c  linear system  a*x = b  in place of  b , given the lu-factorization
c  for  a  in the workarray  w .
c
c******  i n p u t  ******
c  w, nroww,nrow,nbandl,nbandu.....describe the lu-factorization of a
c        banded matrix  a  of order  nrow  as constructed in  banfac .
c        for details, see  banfac .
c  b.....right side of the system to be solved .
c
c******  o u t p u t  ******
c  b.....contains the solution  x , of order  nrow .
c
c******  m e t h o d  ******
c     (with  a = l*u, as stored in  w,) the unit lower triangular system
c  l(u*x) = b  is solved for  y = u*x, and  y  stored in  b . then the
c  upper triangular system  u*x = y  is solved for  x  . the calcul-
c  ations are so arranged that the innermost loops stay within columns.
c
      implicit none

      integer nbandl,nbandu,nrow,nroww,   i,j,jmax,middle,nrowm1
      double precision w(nroww,nrow),b(nrow)
      middle = nbandu + 1
      if (nrow .eq. 1)                  go to 49
      nrowm1 = nrow - 1
      if (nbandl .eq. 0)                go to 30
c                                 forward pass
c            for i=1,2,...,nrow-1, subtract  right side(i)*(i-th column
c            of  l )  from right side  (below i-th row) .
      do 21 i=1,nrowm1
         jmax = min0(nbandl, nrow-i)
         do 21 j=1,jmax
   21       b(i+j) = b(i+j) - b(i)*w(middle+j,i)
c                                 backward pass
c            for i=nrow,nrow-1,...,1, divide right side(i) by i-th diag-
c            onal entry of  u, then subtract  right side(i)*(i-th column
c            of  u)  from right side  (above i-th row).
   30 if (nbandu .gt. 0)                go to 40
c                                a  is lower triangular .
      do 31 i=1,nrow
   31    b(i) = b(i)/w(1,i)
                                        return
   40 i = nrow    
   41    b(i) = b(i)/w(middle,i)
         jmax = min0(nbandu,i-1)
         do 45 j=1,jmax
   45       b(i-j) = b(i-j) - b(i)*w(middle-j,i)
         i = i - 1
         if (i .gt. 1)                  go to 41
   49 b(1) = b(1)/w(middle,1)
                                        return
      end
      subroutine bchfac ( w, nbands, nrow, diag )

c*********************************************************************72
c
cc BCHFAC constructs a Cholesky factorization of a matrix.
c
c  from  * a practical guide to splines *  by c. de boor    
constructs cholesky factorization
c                     c  =  l * d * l-transpose
c  with l unit lower triangular and d diagonal, for given matrix c of
c  order  n r o w , in case  c  is (symmetric) positive semidefinite
c  and  b a n d e d , having  n b a n d s  diagonals at and below the
c  main diagonal.
c
c******  i n p u t  ******
c  nrow.....is the order of the matrix  c .
c  nbands.....indicates its bandwidth, i.e.,
c          c(i,j) = 0 for i-j .ge. nbands .
c  w.....workarray of size (nbands,nrow)  containing the  nbands  diago-
c        nals in its rows, with the main diagonal in row  1 . precisely,
c        w(i,j)  contains  c(i+j-1,j), i=1,...,nbands, j=1,...,nrow.
c          for example, the interesting entries of a seven diagonal sym-
c        metric matrix  c  of order  9  would be stored in  w  as
c
c
c
c
c
c
c        all other entries of  w  not identified in this way with an en-
c        try of  c  are never referenced .
c  diag.....is a work array of length  nrow .
c
c******  o u t p u t  ******
c  w.....contains the cholesky factorization  c = l*d*l-transp, with
c        w(1,i) containing  1/d(i,i)
c        and  w(i,j)  containing  l(i-1+j,j), i=2,...,nbands.
c
c******  m e t h o d  ******
c   gauss elimination, adapted to the symmetry and bandedness of  c , is
c   used .
c     near zero pivots are handled in a special way. the diagonal ele-
c  ment c(n,n) = w(1,n) is saved initially in  diag(n), all n. at the n-
c  th elimination step, the current pivot element, viz.  w(1,n), is com-
c  pared with its original value, diag(n). if, as the result of prior
c  elimination steps, this element has been reduced by about a word
c  length, (i.e., if w(1,n)+diag(n) .le. diag(n)), then the pivot is de-
c  clared to be zero, and the entire n-th row is declared to be linearly
c  dependent on the preceding rows. this has the effect of producing
c   x(n) = 0  when solving  c*x = b  for  x, regardless of  b. justific-
c  ation for this is as follows. in contemplated applications of this
c  program, the given equations are the normal equations for some least-
c  squares approximation problem, diag(n) = c(n,n) gives the norm-square
c  of the n-th basis function, and, at this point,  w(1,n)  contains the
c  norm-square of the error in the least-squares approximation to the n-
c  th basis function by linear combinations of the first n-1 . having
c  w(1,n)+diag(n) .le. diag(n) signifies that the n-th function is lin-
c  early dependent to machine accuracy on the first n-1 functions, there
c  fore can safely be left out from the basis of approximating functions
c     the solution of a linear system
c                       c*x = b
c   is effected by the succession of the following  t w o  calls:
c     call bchfac ( w, nbands, nrow, diag )       , to get factorization
c     call bchslv ( w, nbands, nrow, b )          , to solve for x.
c
      implicit none

      integer nbands,nrow,   i,imax,j,jmax,n
      double precision w(nbands,nrow),diag(nrow),   ratio
      if (nrow .gt. 1)                  go to 9
      if (w(1,1) .gt. 0.0D+00) w(1,1) = 1.0D+00/w(1,1)
                                        return
c                                        store diagonal of  c  in  diag.
    9 do 10 n=1,nrow
   10    diag(n) = w(1,n)
c                                                        factorization .
      do 20 n=1,nrow
         if (w(1,n)+diag(n) .gt. diag(n)) go to 15
         do 14 j=1,nbands
   14       w(j,n) = 0.0D+00
                                        go to 20
   15    w(1,n) = 1./w(1,n)
         imax = min0(nbands-1,nrow - n)
         if (imax .lt. 1)               go to 20
         jmax = imax
         do 18 i=1,imax
            ratio = w(i+1,n)*w(1,n)
            do 17 j=1,jmax
   17          w(j,n+i) = w(j,n+i) - w(j+i,n)*ratio
            jmax = jmax - 1
   18       w(i+1,n) = ratio
   20    continue
                                        return
      end
      subroutine bchslv ( w, nbands, nrow, b )

c*********************************************************************72
c
cc BCHSLV solves a banded symmetric positive definite system.
c
c  from  * a practical guide to splines *  by c. de boor    
c  solves the linear system     c*x = b   of order  n r o w  for  x
c  provided  w  contains the cholesky factorization for the banded (sym-
c  metric) positive definite matrix  c  as constructed in the subroutine
c    b c h f a c  (quo vide).
c
c******  i n p u t  ******
c  nrow.....is the order of the matrix  c .
c  nbands.....indicates the bandwidth of  c .
c  w.....contains the cholesky factorization for  c , as output from
c        subroutine bchfac  (quo vide).
c  b.....the vector of length  n r o w  containing the right side.
c
c******  o u t p u t  ******
c  b.....the vector of length  n r o w  containing the solution.
c
c******  m e t h o d  ******
c  with the factorization  c = l*d*l-transpose  available, where  l  is
c  unit lower triangular and  d  is diagonal, the triangular system
c  l*y = b  is solved for  y (forward substitution), y is stored in  b,
c  the vector  d**(-1)*y is computed and stored in  b, then the triang-
c  ular system  l-transpose*x = d**(-1)*y is solved for  x (backsubstit-
c  ution).
c
      implicit none

      integer nbands,nrow,   j,jmax,n,nbndm1
      double precision w(nbands,nrow),b(nrow)
      if (nrow .gt. 1)                  go to 21
      b(1) = b(1)*w(1,1)
                                        return
c
c                 forward substitution. solve l*y = b for y, store in b.
   21 nbndm1 = nbands - 1
      do 30 n=1,nrow
         jmax = min0(nbndm1,nrow-n)
         if (jmax .lt. 1)               go to 30
         do 25 j=1,jmax
   25       b(j+n) = b(j+n) - w(j+1,n)*b(n)
   30    continue
c
c     backsubstitution. solve l-transp.x = d**(-1)*y  for x, store in b.
      n = nrow    
   39    b(n) = b(n)*w(1,n)   
         jmax = min0(nbndm1,nrow-n)
         if (jmax .lt. 1)               go to 40
         do 35 j=1,jmax
   35       b(n) = b(n) - w(j+1,n)*b(j+n)
   40    n = n-1  
         if (n .gt. 0)                  go to 39
                                        return
      end
      subroutine bsplpp ( t, bcoef, n, k, scrtch, break, coef, l )

c*********************************************************************72
c
cc BSPLPP converts from B-spline to piecewise polynomial form.
c
c  from  * a practical guide to splines *  by c. de Boor (7 may 92)
calls  bsplvb
c
converts the b-representation  t, bcoef, n, k  of some spline into its
c  pp-representation  break, coef, l, k .
c
c******  i n p u t  ******
c  t.....knot sequence, of length  n+k
c  bcoef.....b-spline coefficient sequence, of length  n
c  n.....length of  bcoef  and  dimension of spline space  spline(k,t)
c  k.....order of the spline
c
c  w a r n i n g  . . .  the restriction   k .le. kmax (= 20)   is impo-
c        sed by the arbitrary dimension statement for  biatx  below, but
c        is  n o w h e r e   c h e c k e d   for.
c
c******  w o r k   a r e a  ******
c  scrtch......of size  (k,k) , needed to contain bcoeffs of a piece of
c        the spline and its  k-1  derivatives
c
c******  o u t p u t  ******
c  break.....breakpoint sequence, of length  l+1, contains (in increas-
c        ing order) the distinct points in the sequence  t(k),...,t(n+1)
c  coef.....array of size (k,l), with  coef(i,j) = (i-1)st derivative of
c        spline at break(j) from the right
c  l.....number of polynomial pieces which make up the spline in the in-
c        terval  (t(k), t(n+1))
c
c******  m e t h o d  ******
c     for each breakpoint interval, the  k  relevant b-coeffs of the
c  spline are found and then differenced repeatedly to get the b-coeffs
c  of all the derivatives of the spline on that interval. the spline and
c  its first  k-1  derivatives are then evaluated at the left end point
c  of that interval, using  bsplvb  repeatedly to obtain the values of
c  all b-splines of the appropriate order at that point.
c
      implicit none

      integer k,l,n,   i,j,jp1,kmax,kmj,left,lsofar
      parameter (kmax = 20)
      double precision bcoef(n),break(l+1),coef(k,l),t(n+k), scrtch(k,k)
     &                                      ,biatx(kmax),diff,factor,sum
c
      lsofar = 0
      break(1) = t(k)
      do 50 left=k,n
c                                find the next nontrivial knot interval.
         if (t(left+1) .eq. t(left))    go to 50
         lsofar = lsofar + 1
         break(lsofar+1) = t(left+1)
         if (k .gt. 1)                  go to 9
         coef(1,lsofar) = bcoef(left)
                                        go to 50
c        store the k b-spline coeff.s relevant to current knot interval
c                             in  scrtch(.,1) .
    9    do 10 i=1,k
   10       scrtch(i,1) = bcoef(left-k+i)
c
c        for j=1,...,k-1, compute the  k-j  b-spline coeff.s relevant to
c        current knot interval for the j-th derivative by differencing
c        those for the (j-1)st derivative, and store in scrtch(.,j+1) .
         do 20 jp1=2,k
            j = jp1 - 1
            kmj = k - j
            do 20 i=1,kmj
               diff = t(left+i) - t(left+i - kmj)
               if (diff .gt. 0.0D+00)  scrtch(i,jp1) =
     &                       (scrtch(i+1,j)-scrtch(i,j))/diff
   20          continue
c
c        for  j = 0, ..., k-1, find the values at  t(left)  of the  j+1
c        b-splines of order  j+1  whose support contains the current
c        knot interval from those of order  j  (in  biatx ), then comb-
c        ine with the b-spline coeff.s (in scrtch(.,k-j) ) found earlier
c        to compute the (k-j-1)st derivative at  t(left)  of the given
c        spline.
c           note. if the repeated calls to  bsplvb  are thought to gene-
c        rate too much overhead, then replace the first call by
c           biatx(1) = 1.
c        and the subsequent call by the statement
c           j = jp1 - 1
c        followed by a direct copy of the lines
c           deltar(j) = t(left+j) - x
c                  ......
c           biatx(j+1) = saved
c        from  bsplvb . deltal(kmax)  and  deltar(kmax)  would have to
c        appear in a dimension statement, of course.
c
         call bsplvb ( t, 1, 1, t(left), left, biatx )
         coef(k,lsofar) = scrtch(1,k)
         do 30 jp1=2,k
            call bsplvb ( t, jp1, 2, t(left), left, biatx )
            kmj = k+1 - jp1
            sum = 0.0D+00
            do 28 i=1,jp1
   28          sum = biatx(i)*scrtch(i,kmj) + sum
   30       coef(kmj,lsofar) = sum
   50    continue
      l = lsofar
	  if (k .eq. 1)                     return
	  factor = 1.0D+00
	  do 60 i=2,k
		 factor = factor*dble(k+1-i)
		 do 60 j=1,lsofar
   60       coef(i,j) = coef(i,j)*factor
                                        return
      end
      subroutine bsplvb ( t, jhigh, index, x, left, biatx )

c*********************************************************************72
c
cc BSPLVB evaluates B-splines at a point X with a given knot sequence.
c
c  from  * a practical guide to splines *  by c. de boor    
calculates the value of all possibly nonzero b-splines at  x  of order
c
c               jout  =  max( jhigh , (j+1)*(index-1) )
c
c  with knot sequence  t .
c
c******  i n p u t  ******
c  t.....knot sequence, of length  left + jout  , assumed to be nonde-
c        creasing.  a s s u m p t i o n . . . .
c                       t(left)  .lt.  t(left + 1)   .
c   d i v i s i o n  b y  z e r o  will result if  t(left) = t(left+1)
c  jhigh,
c  index.....integers which determine the order  jout = max(jhigh,
c        (j+1)*(index-1))  of the b-splines whose values at  x  are to
c        be returned.  index  is used to avoid recalculations when seve-
c        ral columns of the triangular array of b-spline values are nee-
c        ded (e.g., in  bsplpp  or in  bsplvd ). precisely,
c                     if  index = 1 ,
c        the calculation starts from scratch and the entire triangular
c        array of b-spline values of orders 1,2,...,jhigh  is generated
c        order by order , i.e., column by column .
c                     if  index = 2 ,
c        only the b-spline values of order  j+1, j+2, ..., jout  are ge-
c        nerated, the assumption being that  biatx , j , deltal , deltar
c        are, on entry, as they were on exit at the previous call.
c           in particular, if  jhigh = 0, then  jout = j+1, i.e., just
c        the next column of b-spline values is generated.
c
c  w a r n i n g . . .  the restriction   jout .le. jmax (= 20)  is im-
c        posed arbitrarily by the dimension statement for  deltal  and
c        deltar  below, but is  n o w h e r e  c h e c k e d  for .
c
c  x.....the point at which the b-splines are to be evaluated.
c  left.....an integer chosen (usually) so that
c                  t(left) .le. x .le. t(left+1)  .
c
c******  o u t p u t  ******
c  biatx.....array of length  jout , with  biatx(i)  containing the val-
c        ue at  x  of the polynomial of order  jout  which agrees with
c        the b-spline  b(left-jout+i,jout,t)  on the interval (t(left),
c        t(left+1)) .
c
c******  m e t h o d  ******
c  the recurrence relation
c
c                       x - t(i)              t(i+j+1) - x
c     b(i,j+1)(x)  =  -----------b(i,j)(x) + ---------------b(i+1,j)(x)
c                     t(i+j)-t(i)            t(i+j+1)-t(i+1)
c
c  is used (repeatedly) to generate the (j+1)-vector  b(left-j,j+1)(x),
c  ...,b(left,j+1)(x)  from the j-vector  b(left-j+1,j)(x),...,
c  b(left,j)(x), storing the new values in  biatx  over the old. the
c  facts that
c            b(i,1) = 1  if  t(i) .le. x .lt. t(i+1)
c  and that
c            b(i,j)(x) = 0  unless  t(i) .le. x .lt. t(i+j)
c  are used. the particular organization of the calculations follows al-
c  gorithm  (8)  in chapter x of the text.
c
      implicit none

      integer index,jhigh,left,   i,j,jmax,jp1
      parameter (jmax = 20)
      double precision biatx(jhigh),t(1),x,deltal(jmax),
     &  deltar(jmax),saved,term
c
c     dimension biatx(jout), t(left+jout)
current fortran standard makes it impossible to specify the length of
c  t  and of  biatx  precisely without the introduction of otherwise
c  superfluous additional arguments.
      data j/1/
      save j,deltal,deltar 
c
                                        go to (10,20), index
   10 j = 1
      biatx(1) = 1.0D+00
      if (j .ge. jhigh)                 go to 99
c
   20    jp1 = j + 1
         deltar(j) = t(left+j) - x
         deltal(j) = x - t(left+1-j)
         saved = 0.0D+00
         do 26 i=1,j
            term = biatx(i)/(deltar(i) + deltal(jp1-i))
            biatx(i) = saved + deltar(i)*term
   26       saved = deltal(jp1-i)*term
         biatx(jp1) = saved
         j = jp1
         if (j .lt. jhigh)              go to 20
c
   99                                   return
      end
      subroutine bsplvd ( t, k, x, left, a, dbiatx, nderiv )

c*********************************************************************72
c
cc BSPLVD calculates the nonvanishing B-splines and derivatives at X.
c
c  from  * a practical guide to splines *  by c. de Boor (7 may 92)    
calls bsplvb
calculates value and deriv.s of all b-splines which do not vanish at x
c
c******  i n p u t  ******
c  t     the knot array, of length left+k (at least)
c  k     the order of the b-splines to be evaluated
c  x     the point at which these values are sought
c  left  an integer indicating the left endpoint of the interval of
c        interest. the  k  b-splines whose support contains the interval
c               (t(left), t(left+1))
c        are to be considered.
c  a s s u m p t i o n  - - -  it is assumed that
c               t(left) .lt. t(left+1)
c        division by zero will result otherwise (in  b s p l v b ).
c        also, the output is as advertised only if
c               t(left) .le. x .le. t(left+1) .
c  nderiv   an integer indicating that values of b-splines and their
c        derivatives up to but not including the  nderiv-th  are asked
c        for. ( nderiv  is replaced internally by the integer  m h i g h
c        in  (1,k)  closest to it.)
c
c******  w o r k   a r e a  ******
c  a     an array of order (k,k), to contain b-coeff.s of the derivat-
c        ives of a certain order of the  k  b-splines of interest.
c
c******  o u t p u t  ******
c  dbiatx   an array of order (k,nderiv). its entry  (i,m)  contains
c        value of  (m-1)st  derivative of  (left-k+i)-th  b-spline of
c        order  k  for knot sequence  t , i=1,...,k, m=1,...,nderiv.
c
c******  m e t h o d  ******
c  values at  x  of all the relevant b-splines of order k,k-1,...,
c  k+1-nderiv  are generated via  bsplvb  and stored temporarily in
c  dbiatx .  then, the b-coeffs of the required derivatives of the b-
c  splines of interest are generated by differencing, each from the pre-
c  ceding one of lower order, and combined with the values of b-splines
c  of corresponding order in  dbiatx  to produce the desired values .
c
      implicit none

      integer k,left,nderiv,   i,ideriv,il,j,jlow,jp1mid,kp1,kp1mm
     *                        ,ldummy,m,mhigh
      double precision a(k,k),dbiatx(k,nderiv),t(1),x,
     &   factor,fkp1mm,sum
      mhigh = max0(min0(nderiv,k),1)
c     mhigh is usually equal to nderiv.
      kp1 = k+1
      call bsplvb(t,kp1-mhigh,1,x,left,dbiatx)
      if (mhigh .eq. 1)                 go to 99
c     the first column of  dbiatx  always contains the b-spline values
c     for the current order. these are stored in column k+1-current
c     order  before  bsplvb  is called to put values for the next
c     higher order on top of it.
      ideriv = mhigh
      do 15 m=2,mhigh
         jp1mid = 1
         do 11 j=ideriv,k
            dbiatx(j,ideriv) = dbiatx(jp1mid,1)
   11       jp1mid = jp1mid + 1
         ideriv = ideriv - 1
         call bsplvb(t,kp1-ideriv,2,x,left,dbiatx)
   15    continue
c
c     at this point,  b(left-k+i, k+1-j)(x) is in  dbiatx(i,j) for
c     i=j,...,k and j=1,...,mhigh ('=' nderiv). in particular, the
c     first column of  dbiatx  is already in final form. to obtain cor-
c     responding derivatives of b-splines in subsequent columns, gene-
c     rate their b-repr. by differencing, then evaluate at  x.
c
      jlow = 1
      do 20 i=1,k
         do 19 j=jlow,k
   19       a(j,i) = 0.0D+00
         jlow = i
   20    a(i,i) = 1.0D+00
c     at this point, a(.,j) contains the b-coeffs for the j-th of the
c     k  b-splines of interest here.
c
      do 40 m=2,mhigh
         kp1mm = kp1 - m
         fkp1mm = dble(kp1mm)
         il = left
         i = k
c
c        for j=1,...,k, construct b-coeffs of  (m-1)st  derivative of
c        b-splines from those for preceding derivative by differencing
c        and store again in  a(.,j) . the fact that  a(i,j) = 0  for
c        i .lt. j  is used.
         do 25 ldummy=1,kp1mm
            factor = fkp1mm/(t(il+kp1mm) - t(il))
c           the assumption that t(left).lt.t(left+1) makes denominator
c           in  factor  nonzero.
            do 24 j=1,i
   24          a(i,j) = (a(i,j) - a(i-1,j))*factor
            il = il - 1
   25       i = i - 1
c
c        for i=1,...,k, combine b-coeffs a(.,i) with b-spline values
c        stored in dbiatx(.,m) to get value of  (m-1)st  derivative of
c        i-th b-spline (of interest here) at  x , and store in
c        dbiatx(i,m). storage of this value over the value of a b-spline
c        of order m there is safe since the remaining b-spline derivat-
c        ives of the same order do not use this value due to the fact
c        that  a(j,i) = 0  for j .lt. i .
         do 40 i=1,k
            sum = 0.0D+00
            jlow = max0(i,m)
            do 35 j=jlow,k
   35          sum = a(j,i)*dbiatx(j,m) + sum
   40       dbiatx(i,m) = sum
   99                                   return
      end
      subroutine bspp2d ( t, bcoef, n, k, m, scrtch, break, coef)       ********    

c*********************************************************************72
c
cc BSPP2D converts from B-spline to piecewise polynomial representation.
c
c  from  * a practical guide to splines *  by c. de Boor    
calls  bsplvb     
c  this is an extended version of  bsplpp  for use with tensor products --------    
c     
converts the b-representation  t, bcoef(.,j), n, k  of some spline into --------    
c  its pp-representation  break, coef(j,.,.), l, k , j=1, ..., m  .     --------    
c     
c******  i n p u t  ******    
c  t.....knot sequence, of length  n+k    
c  bcoef(.,j) b-spline coefficient sequence, of length  n ,j=1,...,m    --------    
c  n.....length of  bcoef  and  dimension of spline space  spline(k,t)  
c  k.....order of the spline  
c     
c  w a r n i n g  . . .  the restriction   k .le. kmax (= 20)   is impo-
c        sed by the arbitrary dimension statement for  biatx  below, but
c        is  n o w h e r e   c h e c k e d   for.     
c     
c  m     number of data sets                                            ********    
c     
c******  w o r k   a r e a  ******  
c  scrtch   of size  (k,k,m), needed to contain bcoeffs of a piece of   ********    
c        the spline and its  k-1  derivatives   for each of the m sets  --------    
c     
c******  o u t p u t  ******  
c  break.....breakpoint sequence, of length  l+1, contains (in increas- 
c        ing order) the distinct points in the sequence  t(k),...,t(n+1)
c  coef(mm,.,.)  array of size (k,n), with  coef(mm,i,j) = (i-1)st der- ********    
c        ivative of  mm-th  spline at break(j) from the right, mm=1,.,m ********    
c     
c******  m e t h o d  ******  
c     for each breakpoint interval, the  k  relevant b-coeffs of the    
c  spline are found and then differenced repeatedly to get the b-coeffs 
c  of all the derivatives of the spline on that interval. the spline and
c  its first  k-1  derivatives are then evaluated at the left end point 
c  of that interval, using  bsplvb  repeatedly to obtain the values of  
c  all b-splines of the appropriate order at that point.    
c     
      implicit none

      integer k,m,n,   i,j,jp1,kmax,kmj,left,lsofar,mm                  ********
      parameter (kmax = 20)     
      double precision bcoef(n,m),break(n+2-k),coef(m,k,n+1-k),t(n+k)               ********
     &                   ,scrtch(k,k,m),biatx(kmax),diff,fkmj,sum
c
      lsofar = 0  
      break(1) = t(k)   
      do 50 left=k,n    
c                                find the next nontrivial knot interval.
         if (t(left+1) .eq. t(left))    go to 50
         lsofar = lsofar + 1  
         break(lsofar+1) = t(left+1)
         if (k .gt. 1)                  go to 9 
         do 5 mm=1,m                                                    ********    
    5       coef(mm,1,lsofar) = bcoef(left,mm)                          ********    
                                        go to 50
c        store the k b-spline coeff.s relevant to current knot interval 
c                             in  scrtch(.,1) . 
    9    do 10 i=1,k    
            do 10 mm=1,m                                                ********    
   10          scrtch(i,1,mm) = bcoef(left-k+i,mm)                      ********    
c     
c        for j=1,...,k-1, compute the  k-j  b-spline coeff.s relevant to
c        current knot interval for the j-th derivative by differencing  
c        those for the (j-1)st derivative, and store in scrtch(.,j+1) . 
         do 20 jp1=2,k  
            j = jp1 - 1 
            kmj = k - j 
            fkmj = dble(kmj) 
            do 20 i=1,kmj     
               diff = (t(left+i) - t(left+i - kmj))/fkmj                --------    
               if (diff .le. 0.0D+00)         go to 20                       --------    
               do 15 mm=1,m                                             ********    
   15             scrtch(i,jp1,mm) =                                    ********    
     *            (scrtch(i+1,j,mm) - scrtch(i,j,mm))/diff              ********    
   20          continue 
c     
c        for  j = 0, ..., k-1, find the values at  t(left)  of the  j+1 
c        b-splines of order  j+1  whose support contains the current    
c        knot interval from those of order  j  (in  biatx ), then comb- 
c        ine with the b-spline coeff.s (in scrtch(.,k-j) ) found earlier
c        to compute the (k-j-1)st derivative at  t(left)  of the given  
c        spline.  
c           note. if the repeated calls to  bsplvb  are thought to gene-
c        rate too much overhead, then replace the first call by   
c           biatx(1) = 1.     
c        and the subsequent call by the statement     
c           j = jp1 - 1 
c        followed by a direct copy of the lines 
c           deltar(j) = t(left+j) - x     
c                  ......     
c           biatx(j+1) = saved
c        from  bsplvb . deltal(kmax)  and  deltar(kmax)  would have to  
c        appear in a dimension statement, of course.  
c     
         call bsplvb ( t, 1, 1, t(left), left, biatx )
         do 25 mm=1,m                                                   ********    
   25       coef(mm,k,lsofar) = scrtch(1,k,mm)                          ********    
         do 30 jp1=2,k  
            call bsplvb ( t, jp1, 2, t(left), left, biatx ) 
            kmj = k+1 - jp1   
            do 30 mm=1,m                                                ********    
               sum = 0.0D+00                                                --------    
               do 28 i=1,jp1                                            --------    
   28             sum = biatx(i)*scrtch(i,kmj,mm) + sum                 ********    
   30          coef(mm,kmj,lsofar) = sum                                ********    
   50    continue 
                                        return  
      end   
      double precision function bvalue ( t, bcoef, n, k, x, jderiv )

c*********************************************************************72
c
cc BVALUE evaluates a derivative of a spline from its B-spline representation.
c
c  from  * a practical guide to splines *  by c. de boor    
calls  interv
c
calculates value at  x  of  jderiv-th derivative of spline from b-repr.
c  the spline is taken to be continuous from the right, EXCEPT at the
c  rightmost knot, where it is taken to be continuous from the left.
c
c******  i n p u t ******
c  t, bcoef, n, k......forms the b-representation of the spline  f  to
c        be evaluated. specifically,
c  t.....knot sequence, of length  n+k, assumed nondecreasing.
c  bcoef.....b-coefficient sequence, of length  n .
c  n.....length of  bcoef  and dimension of spline(k,t),
c        a s s u m e d  positive .
c  k.....order of the spline .
c
c  w a r n i n g . . .   the restriction  k .le. kmax (=20)  is imposed
c        arbitrarily by the dimension statement for  aj, dl, dr  below,
c        but is  n o w h e r e  c h e c k e d  for.
c
c  x.....the point at which to evaluate .
c  jderiv.....integer giving the order of the derivative to be evaluated
c        a s s u m e d  to be zero or positive.
c
c******  o u t p u t  ******
c  bvalue.....the value of the (jderiv)-th derivative of  f  at  x .
c
c******  m e t h o d  ******
c     The nontrivial knot interval  (t(i),t(i+1))  containing  x  is lo-
c  cated with the aid of  interv . The  k  b-coeffs of  f  relevant for
c  this interval are then obtained from  bcoef (or taken to be zero if
c  not explicitly available) and are then differenced  jderiv  times to
c  obtain the b-coeffs of  (d**jderiv)f  relevant for that interval.
c  Precisely, with  j = jderiv, we have from x.(12) of the text that
c
c     (d**j)f  =  sum ( bcoef(.,j)*b(.,k-j,t) )
c
c  where
c                   / bcoef(.),                     ,  j .eq. 0
c                   /
c    bcoef(.,j)  =  / bcoef(.,j-1) - bcoef(.-1,j-1)
c                   / ----------------------------- ,  j .gt. 0
c                   /    (t(.+k-j) - t(.))/(k-j)
c
c     Then, we use repeatedly the fact that
c
c    sum ( a(.)*b(.,m,t)(x) )  =  sum ( a(.,x)*b(.,m-1,t)(x) )
c  with
c                 (x - t(.))*a(.) + (t(.+m-1) - x)*a(.-1)
c    a(.,x)  =    ---------------------------------------
c                 (x - t(.))      + (t(.+m-1) - x)
c
c  to write  (d**j)f(x)  eventually as a linear combination of b-splines
c  of order  1 , and the coefficient for  b(i,1,t)(x)  must then be the
c  desired number  (d**j)f(x). (see x.(17)-(19) of text).
c
      implicit none

      integer jderiv,k,n,   i,ilo,imk,j,jc,jcmin,jcmax,jj,kmax,kmj,km1
     &                     ,mflag,nmi,jdrvp1
      parameter (kmax = 20)
      double precision bcoef(n),t(n+k),x,   
     &  aj(kmax),dl(kmax),dr(kmax),fkmj
      bvalue = 0.0D+00
      if (jderiv .ge. k)                go to 99
c
c  *** Find  i   s.t.   1 .le. i .lt. n+k   and   t(i) .lt. t(i+1)   and
c      t(i) .le. x .lt. t(i+1) . If no such i can be found,  x  lies
c      outside the support of  the spline  f , hence  bvalue = 0.
c      (The asymmetry in this choice of  i  makes  f  rightcontinuous, except
c      at  t(n+k) where it is leftcontinuous.)
      call interv ( t, n+k, x, i, mflag )
      if (mflag .ne. 0)                 go to 99
c  *** if k = 1 (and jderiv = 0), bvalue = bcoef(i).
      km1 = k - 1
      if (km1 .gt. 0)                   go to 1
      bvalue = bcoef(i)
                                        go to 99
c
c  *** store the k b-spline coefficients relevant for the knot interval
c     (t(i),t(i+1)) in aj(1),...,aj(k) and compute dl(j) = x - t(i+1-j),
c     dr(j) = t(i+j) - x, j=1,...,k-1 . set any of the aj not obtainable
c     from input to zero. set any t.s not obtainable equal to t(1) or
c     to t(n+k) appropriately.
    1 jcmin = 1
      imk = i - k
      if (imk .ge. 0)                   go to 8
      jcmin = 1 - imk
      do 5 j=1,i
    5    dl(j) = x - t(i+1-j)
      do 6 j=i,km1
         aj(k-j) = 0.0D+00
    6    dl(j) = dl(i)
                                        go to 10
    8 do 9 j=1,km1
    9    dl(j) = x - t(i+1-j)
c
   10 jcmax = k
      nmi = n - i
      if (nmi .ge. 0)                   go to 18
      jcmax = k + nmi
      do 15 j=1,jcmax
   15    dr(j) = t(i+j) - x
      do 16 j=jcmax,km1
         aj(j+1) = 0.0D+00
   16    dr(j) = dr(jcmax)
                                        go to 20
   18 do 19 j=1,km1
   19    dr(j) = t(i+j) - x
c
   20 do 21 jc=jcmin,jcmax
   21    aj(jc) = bcoef(imk + jc)
c
c               *** difference the coefficients  jderiv  times.
      if (jderiv .eq. 0)                go to 30
      do 23 j=1,jderiv
         kmj = k-j
         fkmj = dble(kmj)
         ilo = kmj
         do 23 jj=1,kmj
            aj(jj) = ((aj(jj+1) - aj(jj))/(dl(ilo) + dr(jj)))*fkmj
   23       ilo = ilo - 1
c
c  *** compute value at  x  in (t(i),t(i+1)) of jderiv-th derivative,
c     given its relevant b-spline coeffs in aj(1),...,aj(k-jderiv).
   30 if (jderiv .eq. km1)              go to 39
      jdrvp1 = jderiv + 1     
      do 33 j=jdrvp1,km1
         kmj = k-j
         ilo = kmj
         do 33 jj=1,kmj
            aj(jj) = (aj(jj+1)*dl(ilo) + aj(jj)*dr(jj))/(dl(ilo)+dr(jj))
   33       ilo = ilo - 1
   39 bvalue = aj(1)
c
   99                                   return
      end
      subroutine chol1d ( p, v, qty, npoint, ncol, u, qu )

c*********************************************************************72
c
cc CHOL1D sets up and solves linear systems needed by SMOOTH.
c
c  from  * a practical guide to splines *  by c. de boor    
c  from  * a practical guide to splines *  by c. de boor    
c to be called in  s m o o t h
constructs the upper three diags. in v(i,j), i=2,,npoint-1, j=1,3, of
c  the matrix  6*(1-p)*q-transp.*(d**2)*q + p*r, then computes its
c  l*l-transp. decomposition and stores it also in v, then applies
c  forward and backsubstitution to the right side q-transp.*y in  qty
c  to obtain the solution in  u .
c
      implicit none

      integer ncol,npoint,   i,npm1,npm2
      double precision p,qty(npoint),qu(npoint),u(npoint),v(npoint,7),
     &   prev,ratio
     &    ,six1mp,twop
      npm1 = npoint - 1
c     construct 6*(1-p)*q-transp.*(d**2)*q  +  p*r
      six1mp = 6.0D+00*(1.0D+00-p)
      twop = 2.0D+00*p
      do 2 i=2,npm1
         v(i,1) = six1mp*v(i,5) + twop*(v(i-1,4)+v(i,4))
         v(i,2) = six1mp*v(i,6) + p*v(i,4)
    2    v(i,3) = six1mp*v(i,7)
      npm2 = npoint - 2
      if (npm2 .ge. 2)                  go to 10
      u(1) = 0.0D+00
      u(2) = qty(2)/v(2,1)
      u(3) = 0.0D+00
                                        go to 41
c  factorization
   10 do 20 i=2,npm2
         ratio = v(i,2)/v(i,1)
         v(i+1,1) = v(i+1,1) - ratio*v(i,2)
         v(i+1,2) = v(i+1,2) - ratio*v(i,3)
         v(i,2) = ratio
         ratio = v(i,3)/v(i,1)
         v(i+2,1) = v(i+2,1) - ratio*v(i,3)
   20    v(i,3) = ratio
c
c  forward substitution
      u(1) = 0.0D+00
      v(1,3) = 0.0D+00
      u(2) = qty(2)
      do 30 i=2,npm2
   30    u(i+1) = qty(i+1) - v(i,2)*u(i) - v(i-1,3)*u(i-1)
c  back substitution
      u(npoint) = 0.0D+00
      u(npm1) = u(npm1)/v(npm1,1)
      i = npm2    
   40    u(i) = u(i)/v(i,1)-u(i+1)*v(i,2)-u(i+2)*v(i,3)
         i = i - 1
         if (i .gt. 1)                  go to 40
c  construct q*u
   41 prev = 0.0D+00
      do 50 i=2,npoint
         qu(i) = (u(i) - u(i-1))/v(i-1,4)
         qu(i-1) = qu(i) - prev
   50    prev = qu(i)
      qu(npoint) = -qu(npoint)
                                        return
      end
      subroutine colloc(aleft,aright,lbegin,iorder,ntimes,addbrk,relerr)

c*********************************************************************72
c
cc COLLOC solves an ordinary differential equation by collocation.
c
c  from  * a practical guide to splines *  by c. de boor    
chapter xv, example. solution of an ode by collocation.
calls colpnt, difequ(ppvalu(interv)), knots, eqblok(putit(difequ*,
c     bsplvd(bsplvb)))), slvblk(various subprograms), bsplpp(bsplvb*),
c     newnot
c
c******  i n p u t  ******
c  aleft, aright  endpoints of interval of approximation
c  lbegin   initial number of polynomial pieces in the approximation.
c           a uniform breakpoint sequence is chosen.
c  iorder   order of polynomial pieces in the approximation
c  ntimes   number of passes through  n e w n o t  to be made
c  addbrk   the number (possibly fractional) of breaks to be added per
c           pass through newnot. e.g., if addbrk = .33334, then a break-
c           point will be added at every third pass through newnot.
c  relerr   a tolerance. newton iteration is stopped if the difference
c           between the b-coeffs of two successive iterates is no more
c           than  relerr*(absol.largest b-coefficient).
c
c******  p r i n t e d   o u t p u t  ******
c     consists of the pp-representation of the approximate solution,
c     and of the error at selected points.
c
c******  m e t h o d  ******
c  the m-th order ordinary differential equation with  m  side condit-
c  ions, to be specified in subroutine  d i f e q u , is solved approx-
c  imately by collocation.
c  the approximation  f  to the solution  g  is pp of order  k+m  with
c  l  pieces and  m-1 continuous derivatives.  f  is determined by the
c  requirement that it satisfy the d.e. at  k  points per interval (to
c  be specified in  c o l p n t ) and the  m  side conditions.
c     this usually nonlinear system of equations for  f  is solved by
c  newton's method. the resulting linear system for the b-coeffs of an
c  iterate is constructed appropriately in  e q b l o k  and then solved
c  in  s l v b l k , a program designed to solve  a l m o s t  b l o c k
c  d i a g o n a l  linear systems efficiently.
c     there is an opportunity to attempt improvement of the breakpoint
c  sequence (both in number and location) through use of  n e w n o t .
c
      implicit none

      integer npiece
      parameter (npiece=100)
      integer iorder,lbegin,ntimes,   i,iflag,ii,integs(3,npiece),iside
     &  ,iter,itermx,j,k,kmax,kpm,l,lenblk,lnew,m,n,nbloks
     &                  ,ndim,ncoef,nncoef,nt
      parameter (ndim=200,kmax=20,ncoef=npiece*kmax,lenblk=ncoef)
      double precision addbrk,aleft,aright,relerr,
     &     a(ndim),amax,asave(ndim)
     &     ,b(ndim),bloks(lenblk),break,coef,dx,err,rho,t(ndim)
     &     ,templ(lenblk),temps(ndim),xside
      equivalence (bloks,templ)
      common /approx/ break(npiece), coef(ncoef), l,kpm
      common /side/ m, iside, xside(10)
      common /other/ itermx,k,rho(kmax-1)
c
      kpm = iorder
      if (lbegin*kpm .gt. ncoef)        go to 999
c  *** set the various parameters concerning the particular dif.equ.
c     including a first approx. in case the de is to be solved by
c     iteration ( itermx .gt. 0) .
      call difequ (1, temps(1), temps )
c  *** obtain the  k  collocation points for the standard interval.
      k = kpm - m
      call colpnt ( k, rho )
c  *** the following five statements could be replaced by a read in or-
c     der to obtain a specific (nonuniform) spacing of the breakpnts.
      dx = (aright - aleft)/dble(lbegin)
      temps(1) = aleft
      do 4 i=2,lbegin
    4    temps(i) = temps(i-1) + dx
      temps(lbegin+1) = aright
c  *** generate, in knots, the required knots t(1),...,t(n+kpm).
      call knots ( temps, lbegin, kpm, t, n )
      nt = 1
c  *** generate the almost block diagonal coefficient matrix  bloks  and
c     right side  b  from collocation equations and side conditions.
c     then solve via  slvblk , obtaining the b-representation of the ap-
c     proximation in  t , a ,  n  , kpm  .
   10    call eqblok(t,n,kpm,temps,a,bloks,lenblk,integs,nbloks,b)
         call slvblk(bloks,integs,nbloks,b,temps,a,iflag)
         iter = 1
         if (itermx .le. 1)             go to 30
c  *** save b-spline coeff. of current approx. in  asave , then get new
c     approx. and compare with old. if coeff. are more than  relerr
c     apart (relatively) or if no. of iterations is less than  itermx ,
c     continue iterating.
   20       call bsplpp(t,a,n,kpm,templ,break,coef,l)
            do 25 i=1,n
   25          asave(i) = a(i)
            call eqblok(t,n,kpm,temps,a,bloks,lenblk,integs,nbloks,b)
            call slvblk(bloks,integs,nbloks,b,temps,a,iflag)
            err = 0.0D+00
            amax = 0.0D+00
            do 26 i=1,n
               amax = dmax1(amax,dabs(a(i)))
   26          err = dmax1(err,dabs(a(i)-asave(i)))
            if (err .le. relerr*amax)   go to 30
            iter = iter+1
            if (iter .lt. itermx)       go to 20
c  *** iteration (if any) completed. print out approx. based on current
c     breakpoint sequence, then try to improve the sequence.
   30    print 630,kpm,l,n,(break(i),i=2,l)
  630    format(47h approximation from a space of splines of order,i3
     &         ,4h on ,i3,11h intervals,/13h of dimension,i4
     &         ,16h.  breakpoints -/(5e20.10))
         if (itermx .gt. 0)  print 635,iter,itermx
  635    format(6h after,i3,3h of,i3,20h allowed iterations,)
         call bsplpp(t,a,n,kpm,templ,break,coef,l)
         print 637
  637    format(46h the pp representation of the approximation is)
         do 38 i=1,l
            ii = (i-1)*kpm
   38       print 638, break(i),(coef(ii+j),j=1,kpm)
  638    format(f9.3,e13.6,10e11.3)
c  *** the following call is provided here for possible further analysis
c     of the approximation specific to the problem being solved.
c     it is, of course, easily omitted.
         call difequ ( 4, temps(1), temps )
c
         if (nt .gt. ntimes)            return
c  *** from the pp-rep. of the current approx., obtain in  newnot  a new
c     (and possibly better) sequence of breakpoints, adding (on the ave-
c     rage)  a d d b r k  breakpoints per pass through newnot.
         lnew = lbegin + int(dble(nt)*addbrk)
         if (lnew*kpm .gt. ncoef)       go to 999
         call newnot(break,coef,l,kpm,temps,lnew,templ)
         call knots(temps,lnew,kpm,t,n)
         nt = nt + 1
                                        go to 10
  999 nncoef = ncoef
      print 699,nncoef
  699 format(11h **********/23h the assigned dimension,i5
     &      ,25h for  coef  is too small.)
                                        return
      end
      subroutine colpnt(k,rho)

c*********************************************************************72
c
cc COLPNT supplies collocation points.
c
c  from  * a practical guide to splines *  by c. de boor    
c  the  k collocation points for the standard interval (-1,1) are sup-
c  plied here as the zeros of the legendre polynomial of degree k ,
c  provided  k .le. 8 . otherwise, uniformly spaced points are given.
c
      implicit none

      integer k,  j
      double precision rho(k),  fkm1o2
      if (k .gt. 8)                    go to 99
                                       go to (10,20,30,40,50,60,70,80),k
   10 rho(1) = 0.0D+00
                                       return
   20 rho(2) = .577350269189626D+00
      rho(1) = - rho(2)
                                       return
   30 rho(3) = .774596669241483D+00
      rho(1) = - rho(3)
      rho(2) = 0.0D+00
                                       return
   40 rho(3) = .339981043584856D+00
      rho(2) = - rho(3)
      rho(4) = .861136311594053D+00
      rho(1) = - rho(4)
                                       return
   50 rho(4) = .538469310105683D+00
      rho(2) = - rho(4)
      rho(5) = .906179845938664D+00
      rho(1) = - rho(5)
      rho(3) = 0.0D+00
                                       return
   60 rho(4) = .238619186083197D+00
      rho(3) = - rho(4)
      rho(5) = .66120 9386466265D+00
      rho(2) = - rho(5)
      rho(6) = .93246 9514203152D+00
      rho(1) = - rho(6)
                                       return
   70 rho(5) = .405845151377397D+00
      rho(3) = - rho(5)
      rho(6) = .741531185599394D+00
      rho(2) = - rho(6)
      rho(7) = .949107912342759D+00
      rho(1) = - rho(7)
      rho(4) = 0.0D+00
                                       return
   80 rho(5) = .183434642495650D+00
      rho(4) = - rho(5)
      rho(6) = .525532409916329D+00
      rho(3) = - rho(6)
      rho(7) = .796666477413627D+00
      rho(2) = - rho(7)
      rho(8) = .960289856497536D+00
      rho(1) = - rho(8)
                                       return
c  if k .gt. 8, use equispaced points, but print warning
   99 print 699,k
  699 format(11h **********/
     &       49h equispaced collocation points are used since k =,i2,
     &       19h is greater than 8.)
      fkm1o2 = dble(k-1)/2.0D+00
      do 100 j=1,k
  100    rho(j) = dble(j-1)/fkm1o2 - 1.0D+00
                                       return
      end
      subroutine cubspl ( tau, c, n, ibcbeg, ibcend )

c*********************************************************************72
c
cc CUBSPL defines an interpolatory cubic spline.
c
c  from  * a practical guide to splines *  by c. de boor    
c     ************************  input  ***************************
c     n = number of data points. assumed to be .ge. 2.
c     (tau(i), c(1,i), i=1,...,n) = abscissae and ordinates of the
c        data points. tau is assumed to be strictly increasing.
c     ibcbeg, ibcend = boundary condition indicators, and
c     c(2,1), c(2,n) = boundary condition information. specifically,
c        ibcbeg = 0  means no boundary condition at tau(1) is given.
c           in this case, the not-a-knot condition is used, i.e. the
c           jump in the third derivative across tau(2) is forced to
c           zero, thus the first and the second cubic polynomial pieces
c           are made to coincide.)
c        ibcbeg = 1  means that the slope at tau(1) is made to equal
c           c(2,1), supplied by input.
c        ibcbeg = 2  means that the second derivative at tau(1) is
c           made to equal c(2,1), supplied by input.
c        ibcend = 0, 1, or 2 has analogous meaning concerning the
c           boundary condition at tau(n), with the additional infor-
c           mation taken from c(2,n).
c     ***********************  output  **************************
c     c(j,i), j=1,...,4; i=1,...,l (= n-1) = the polynomial coefficients
c        of the cubic interpolating spline with interior knots (or
c        joints) tau(2), ..., tau(n-1). precisely, in the interval
c        (tau(i), tau(i+1)), the spline f is given by
c           f(x) = c(1,i)+h*(c(2,i)+h*(c(3,i)+h*c(4,i)/3.)/2.)
c        where h = x - tau(i). the function program *ppvalu* may be
c        used to evaluate f or its derivatives from tau,c, l = n-1,
c        and k=4.
c
      implicit none

      integer ibcbeg,ibcend,n,   i,j,l,m
      double precision c(4,n),tau(n),   divdf1,divdf3,dtau,g
c****** a tridiagonal linear system for the unknown slopes s(i) of
c  f  at tau(i), i=1,...,n, is generated and then solved by gauss elim-
c  ination, with s(i) ending up in c(2,i), all i.
c     c(3,.) and c(4,.) are used initially for temporary storage.
      l = n-1
compute first differences of tau sequence and store in c(3,.). also,
compute first divided difference of data and store in c(4,.).
      do 10 m=2,n
         c(3,m) = tau(m) - tau(m-1)
   10    c(4,m) = (c(1,m) - c(1,m-1))/c(3,m)
construct first equation from the boundary condition, of the form
c             c(4,1)*s(1) + c(3,1)*s(2) = c(2,1)
      if (ibcbeg-1)                     11,15,16
   11 if (n .gt. 2)                     go to 12
c     no condition at left end and n = 2.
      c(4,1) = 1.0D+00
      c(3,1) = 1.0D+00
      c(2,1) = 2.0D+00*c(4,2)
                                        go to 25
c     not-a-knot condition at left end and n .gt. 2.
   12 c(4,1) = c(3,3)
      c(3,1) = c(3,2) + c(3,3)
      c(2,1) =((c(3,2)+2.0D+00*c(3,1))*c(4,2)*c(3,3)
     &  +c(3,2)**2*c(4,3))/c(3,1)
                                        go to 19
c     slope prescribed at left end.
   15 c(4,1) = 1.0D+00
      c(3,1) = 0.0D+00
                                        go to 18
c     second derivative prescribed at left end.
   16 c(4,1) = 2.0D+00
      c(3,1) = 1.0D+00
      c(2,1) = 3.0D+00*c(4,2) - c(3,2)/2.0D+00*c(2,1)
   18 if(n .eq. 2)                      go to 25
c  if there are interior knots, generate the corresp. equations and car-
c  ry out the forward pass of gauss elimination, after which the m-th
c  equation reads    c(4,m)*s(m) + c(3,m)*s(m+1) = c(2,m).
   19 do 20 m=2,l
         g = -c(3,m+1)/c(4,m-1)
         c(2,m) = g*c(2,m-1) + 3.0D+00*(c(3,m)*c(4,m+1)+c(3,m+1)*c(4,m))
   20    c(4,m) = g*c(3,m-1) + 2.0D+00*(c(3,m) + c(3,m+1))
construct last equation from the second boundary condition, of the form
c           (-g*c(4,n-1))*s(n-1) + c(4,n)*s(n) = c(2,n)
c     if slope is prescribed at right end, one can go directly to back-
c     substitution, since c array happens to be set up just right for it
c     at this point.
      if (ibcend-1)                     21,30,24
   21 if (n .eq. 3 .and. ibcbeg .eq. 0) go to 22
c     not-a-knot and n .ge. 3, and either n.gt.3 or  also not-a-knot at
c     left end point.
      g = c(3,n-1) + c(3,n)
      c(2,n) = ((c(3,n)+2.0D+00*g)*c(4,n)*c(3,n-1)
     *            + c(3,n)**2*(c(1,n-1)-c(1,n-2))/c(3,n-1))/g
      g = -g/c(4,n-1)
      c(4,n) = c(3,n-1)
                                        go to 29
c     either (n=3 and not-a-knot also at left) or (n=2 and not not-a-
c     knot at left end point).
   22 c(2,n) = 2.0D+00*c(4,n)
      c(4,n) = 1.0D+00
                                        go to 28
c     second derivative prescribed at right endpoint.
   24 c(2,n) = 3.0D+00*c(4,n) + c(3,n)/2.0D+00*c(2,n)
      c(4,n) = 2.0D+00
                                        go to 28
   25 if (ibcend-1)                     26,30,24
   26 if (ibcbeg .gt. 0)                go to 22
c     not-a-knot at right endpoint and at left endpoint and n = 2.
      c(2,n) = c(4,n)
                                        go to 30
   28 g = -1.0D+00/c(4,n-1)
complete forward pass of gauss elimination.
   29 c(4,n) = g*c(3,n-1) + c(4,n)
      c(2,n) = (g*c(2,n-1) + c(2,n))/c(4,n)
carry out back substitution
   30 j = l 
   40    c(2,j) = (c(2,j) - c(3,j)*c(2,j+1))/c(4,j)
         j = j - 1
         if (j .gt. 0)                  go to 40
c****** generate cubic coefficients in each interval, i.e., the deriv.s
c  at its left endpoint, from value and slope at its endpoints.
      do 50 i=2,n
         dtau = c(3,i)
         divdf1 = (c(1,i) - c(1,i-1))/dtau
         divdf3 = c(2,i-1) + c(2,i) - 2.0D+00*divdf1
         c(3,i-1) = 2.0D+00*(divdf1 - c(2,i-1) - divdf3)/dtau
   50    c(4,i-1) = (divdf3/dtau)*(6.0D+00/dtau)
                                        return
      end
      subroutine cwidth ( w,b,nequ,ncols,integs,nbloks, d, x,iflag )    

c*********************************************************************72
c
cc CWIDTH solves an almost block diagonal linear system.
c
c  this program is a variation of the theme in the algorithm bandet1    
c  by martin and wilkinson (numer.math. 9(1976)279-307). it solves
c  the linear system    
c                           a*x  =  b     
c  of  nequ  equations in case  a  is almost block diagonal with all    
c  blocks having  ncols  columns using no more storage than it takes to 
c  store the interesting part of  a . such systems occur in the determ- 
c  ination of the b-spline coefficients of a spline approximation.
c     
c                           parameters    
c  w     on input, a two-dimensional array of size (nequ,ncols) contain-
c        ing the interesting part of the almost block diagonal coeffici-
c        ent matrix  a (see description and example below). the array   
c        integs  describes the storage scheme.  
c        on output, w  contains the upper triangular factor  u  of the  
c        lu factorization of a possibly permuted version of  a . in par-
c        ticular, the determinant of  a  could now be found as    
c            iflag*w(1,1)*w(2,1)* ... * w(nequ,1)  .  
c  b     on input, the right side of the linear system, of length  nequ.
c        the contents of  b  are changed during execution.  
c  nequ  number of equations in system    
c  ncols block width, i.e., number of columns in each block.
c  integs  integer array, of size (2,nequ), describing the block struct-
c        ure of  a .    
c           integs(1,i) = no. of rows in block i               =  nrow  
c           integs(2,i) = no. of elimination steps in block i     
c                       = overhang over next block             =  last  
c  nbloks  number of blocks   
c  d     work array, to contain row sizes . if storage is scarce, the   
c        array  x  could be used in the calling sequence for  d . 
c  x     on output, contains computed solution (if iflag .ne. 0), of    
c        length  nequ . 
c  iflag  on output, integer  
c         = (-1)**(no.of interchanges during elimination)   
c                if  a  is invertible     
c         =  0   if  a  is singular 
c     
c        ------  block structure of  a  ------  
c     the interesting part of  a  is taken to consist of  nbloks  con-  
c  secutive blocks, with the i-th block made up of  nrowi = integs(1,i) 
c  consecutive rows and  ncols  consecutive columns of  a , and with    
c  the first  lasti = integs(2,i) columns to the left of the next block.
c  these blocks are stored consecutively in the workarray  w .    
c     for example, here is an 11th order matrix and its arrangement in  
c  the workarray  w . (the interesting entries of  a  are indicated by  
c  their row and column index modulo 10.) 
c     
c                  ---   a   ---                          ---   w   --- 
c     
c                     nrow1=3 
c          11 12 13 14                                     11 12 13 14  
c          21 22 23 24                                     21 22 23 24  
c          31 32 33 34      nrow2=2                        31 32 33 34  
c   last1=2      43 44 45 46                               43 44 45 46  
c                53 54 55 56         nrow3=3               53 54 55 56  
c         last2=3         66 67 68 69                      66 67 68 69  
c                         76 77 78 79                      76 77 78 79  
c                         86 87 88 89   nrow4=1            86 87 88 89  
c                  last3=1   97 98 99 90   nrow5=2         97 98 99 90  
c                     last4=1   08 09 00 01                08 09 00 01  
c                               18 19 10 11                18 19 10 11  
c                        last5=4    
c     
c  for this interpretation of  a  as an almost block diagonal matrix,   
c  we have  nbloks = 5 , and the integs array is
c     
c                        i=  1   2   3   4   5  
c                  k=   
c  integs(k,i) =      1      3   2   3   1   2  
c                     2      2   3   1   1   4  
c     
c       --------  method  --------  
c     gauss elimination with scaled partial pivoting is used, but mult- 
c  ipliers are  n o t  s a v e d  in order to save storage. rather, the 
c  right side is operated on during elimination.
c     the two parameters
c                  i p v t e q   and  l a s t e q     
c  are used to keep track of the action.  ipvteq is the index of the    
c  variable to be eliminated next, from equations  ipvteq+1,...,lasteq, 
c  using equation  ipvteq (possibly after an interchange) as the pivot  
c  equation. the entries in the pivot column are  a l w a y s  in column
c  1 of  w . this is accomplished by putting the entries in rows  
c  ipvteq+1,...,lasteq  revised by the elimination of the  ipvteq-th    
c  variable one to the left in  w . in this way, the columns of the     
c  equations in a given block (as stored in  w ) will be aligned with   
c  those of the next block at the moment when these next equations be-  
c  come involved in the elimination process.    
c     thus, for the above example, the first elimination steps proceed  
c  as follows.    
c     
c  *11 12 13 14    11 12 13 14    11 12 13 14    11 12 13 14
c  *21 22 23 24   *22 23 24       22 23 24       22 23 24   
c  *31 32 33 34   *32 33 34      *33 34          33 34
c   43 44 45 46    43 44 45 46   *43 44 45 46   *44 45 46        etc.   
c   53 54 55 56    53 54 55 56   *53 54 55 56   *54 55 56   
c   66 67 68 69    66 67 68 69    66 67 68 69    66 67 68 69
c        .              .              .              .     
c     
c     in all other respects, the procedure is standard, including the   
c  scaled partial pivoting.   
c     
      implicit none

      integer nbloks

      integer iflag,integs(2,nbloks),ncols,nequ,   i,ii,icount,ipvteq   
     &  ,ipvtp1,istar,j,jmax,lastcl,lasteq,lasti,nexteq,nrowad   
      double precision b(nequ),d(nequ),w(nequ,ncols),x(nequ),
     &   awi1od,colmax,ratio,rowmax,sum,temp   
      iflag = 1   
      ipvteq = 0  
      lasteq = 0  
c                                 the i-loop runs over the blocks 
      do 50 i=1,nbloks  
c     
c        the equations for the current block are added to those current-
c        ly involved in the elimination process, by increasing  lasteq  
c        by  integs(1,i) after the rowsize of these equations has been  
c        recorded in the array  d . 
c     
         nrowad = integs(1,i) 
         do 10 icount=1,nrowad
            nexteq = lasteq + icount
            rowmax = 0.0D+00 
            do 5 j=1,ncols    
    5          rowmax = dmax1(rowmax,dabs(w(nexteq,j)))
            if (rowmax .eq. 0.0D+00)         go to 999     
   10       d(nexteq) = rowmax
         lasteq = lasteq + nrowad   
c     
c        there will be  lasti = integs(2,i)  elimination steps before   
c        the equations in the next block become involved. further,
c        l a s t c l  records the number of columns involved in the cur-
c        rent elimination step. it starts equal to  ncols  when a block 
c        first becomes involved and then drops by one after each elim-  
c        ination step.  
c     
         lastcl = ncols 
         lasti = integs(2,i)  
         do 30 icount=1,lasti 
            ipvteq = ipvteq + 1     
            if (ipvteq .lt. lasteq)     go to 11
            if ( dabs(w(ipvteq,1))+d(ipvteq) .gt. d(ipvteq) )
     *                                  go to 50
                                        go to 999     
c     
c        determine the smallest  i s t a r  in  (ipvteq,lasteq)  for    
c        which  abs(w(istar,1))/d(istar)  is as large as possible, and  
c        interchange equations  ipvteq  and  istar  in case  ipvteq     
c        .lt. istar .   
c     
   11       colmax = dabs(w(ipvteq,1))/d(ipvteq) 
            istar = ipvteq    
            ipvtp1 = ipvteq + 1     
            do 13 ii=ipvtp1,lasteq  
               awi1od = dabs(w(ii,1))/d(ii)
               if (awi1od .le. colmax)  go to 13
               colmax = awi1od
               istar = ii     
   13          continue 
            if ( dabs(w(istar,1))+d(istar) .eq. d(istar) )   
     &                                  go to 999     
            if (istar .eq. ipvteq)      go to 16
            iflag = -iflag    
            temp = d(istar)   
            d(istar) = d(ipvteq)    
            d(ipvteq) = temp  
            temp = b(istar)   
            b(istar) = b(ipvteq)    
            b(ipvteq) = temp  
            do 14 j=1,lastcl  
               temp = w(istar,j)    
               w(istar,j) = w(ipvteq,j)   
   14          w(ipvteq,j) = temp   
c     
c        subtract the appropriate multiple of equation  ipvteq  from    
c        equations  ipvteq+1,...,lasteq to make the coefficient of the  
c        ipvteq-th unknown (presently in column 1 of  w ) zero, but     
c        store the new coefficients in  w  one to the left from the old.
c     
   16       do 20 ii=ipvtp1,lasteq  
               ratio = w(ii,1)/w(ipvteq,1)
               do 18 j=2,lastcl     
   18             w(ii,j-1) = w(ii,j) - ratio*w(ipvteq,j)   
               w(ii,lastcl) = 0.0D+00   
   20          b(ii) = b(ii) - ratio*b(ipvteq)  
   30       lastcl = lastcl - 1     
   50    continue 
c     
c  at this point,  w  and  b  contain an upper triangular linear system 
c  equivalent to the original one, with  w(i,j) containing entry  
c  (i, i-1+j ) of the coefficient matrix. solve this system by backsub- 
c  stitution, taking into account its block structure.
c     
c                          i-loop over the blocks, in reverse order     
      i = nbloks  
   59    lasti = integs(2,i)  
         jmax = ncols - lasti 
         do 70 icount=1,lasti 
            sum = 0.0D+00   
            if (jmax .eq. 0)            go to 61
            do 60 j=1,jmax    
   60          sum = sum + x(ipvteq+j)*w(ipvteq,j+1)  
   61       x(ipvteq) = (b(ipvteq)-sum)/w(ipvteq,1)   
            jmax = jmax + 1   
   70       ipvteq = ipvteq - 1     
         i = i - 1
         if (i .gt. 0)                  go to 59
                                        return  
  999 iflag = 0   
                                        return  
      end   
      subroutine difequ ( mode, xx, v )

c*********************************************************************72
c
cc DIFEQU returns information about a differential equation.
c
c  from  * a practical guide to splines *  by c. de boor    
calls  ppvalu(interv)
c  to be called by   c o l l o c ,  p u t i t
c  information about the differential equation is dispensed from here
c
c******  i n p u t  ******
c  mode  an integer indicating the task to be performed.
c        = 1   initialization
c        = 2   evaluate  de  at  xx
c        = 3   specify the next side condition
c        = 4   analyze the approximation
c  xx a point at which information is wanted
c
c******  o u t p u t  ******
c  v  depends on the  mode  . see comments below
c
      implicit none

      integer mode,   i,iside,itermx,k,kmax,kpm,l,m,ncoef,npiece
      parameter (npiece=100, kmax=20, ncoef=npiece*kmax)
      double precision v(kmax),xx,
     &   break,coef,eps,ep1,ep2,error,factor,ppvalu,rho,solutn
     &                ,s2ovep,un,x,xside
      common /approx/ break(npiece),coef(ncoef),l,kpm
      common /side/ m,iside,xside(10)
      common /other/ itermx,k,rho(kmax-1)
      save eps,factor,s2ovep
c
c  this sample of  difequ  is for the example in chapter xv. it is a
c  nonlinear second order two point boundary value problem.
c
                                        go to (10,20,30,40),mode
c  initialize everything
c  i.e. set the order  m  of the dif.equ., the nondecreasing sequence
c  xside(i),i=1,...,m, of points at which side cond.s are given and
c  anything else necessary.
   10 m = 2
      xside(1) = 0.0D+00
      xside(2) = 1.0D+00
c  *** print out heading
      print 499
  499 format(' carrier,s nonlinear perturb. problem')
      eps = 0.5D-02
      print 610, eps
  610 format(' eps ',e20.10)
c  *** set constants used in formula for solution below.
      factor = (dsqrt(2.0D+00) + dsqrt(3.0D+00))**2
      s2ovep = dsqrt(2.0D+00/eps)
c  *** initial guess for newton iteration. un(x) = x*x - 1.
      l = 1
      break(1) = 0.0D+00
      do 16 i=1,kpm
   16    coef(i) = 0.0D+00
      coef(1) = -1.0D+00
      coef(3) = 2.0D+00
      itermx = 10
                                        return
c
c  provide value of left side coeff.s and right side at  xx .
c  specifically, at  xx  the dif.equ. reads
c        v(m+1)d**m + v(m)d**(m-1) + ... + v(1)d**0  =  v(m+2)
c  in terms of the quantities v(i),i=1,...,m+2, to be computed here.
   20 continue
      v(3) = eps
      v(2) = 0.0D+00
      un = ppvalu(break,coef,l,kpm,xx,0)
      v(1) = 2.0D+00*un
      v(4) = un**2 + 1.0D+00
                                        return
c
c  provide the  m  side conditions. these conditions are of the form
c        v(m+1)d**m + v(m)d**(m-1) + ... + v(1)d**0  =  v(m+2)
c  in terms of the quantities v(i),i=1,...,m+2, to be specified here.
c  note that v(m+1) = 0  for customary side conditions.
   30 v(m+1) = 0.0D+00
                                        go to (31,32,39),iside
   31 v(2) = 1.0D+00
      v(1) = 0.0D+00
      v(4) = 0.0D+00
                                        go to 38
   32 v(2) = 0.0D+00
      v(1) = 1.0D+00
      v(4) = 0.0D+00
   38 iside = iside + 1
   39                                   return
c
c  calculate the error near the boundary layer at  1.
   40 continue
      print 640
  640 format(' x, g(x)  and  g(x)-f(x)  at selected points')
      x = .75
      do 41 i=1,9
         ep1 = exp(s2ovep*(1.0D+00-x))*factor
         ep2 = exp(s2ovep*(1.0D+00+x))*factor
         solutn = 12.0D+00/(1.0D+00+ep1)**2*ep1 
     &    + 12.0D+00/(1.0D+00+ep2)**2*ep2 - 1.0D+00
         error = solutn - ppvalu(break,coef,l,kpm,x,0)
         print 641,x,solutn,error
  641    format(3e20.10)
   41    x = x + 0.03125D+00
                                        return
      end
      subroutine dtblok ( bloks, integs, nbloks, ipivot, iflag,
     &                    detsgn, detlog )

c*********************************************************************72
c
cc DTBLOK gets the determinant of an almost block diagonal matrix.
c
c  computes the determinant of an almost block diagonal matrix whose
c  plu factorization has been obtained previously in fcblok.
c  *** the logarithm of the determinant is computed instead of the
c  determinant itself to avoid the danger of overflow or underflow
c  inherent in this calculation.
c
c parameters
c    bloks, integs, nbloks, ipivot, iflag  are as on return from fcblok.
c            in particular, iflag = (-1)**(number of interchanges dur-
c            ing factorization) if successful, otherwise iflag = 0.
c    detsgn  on output, contains the sign of the determinant.
c    detlog  on output, contains the natural logarithm of the determi-
c            nant if determinant is not zero. otherwise contains 0.
c
      implicit none

      integer nbloks

      integer index,integs(3,nbloks),ipivot(1),iflag, i,indexp,ip,
     &  k,last,nrow
      double precision bloks(1),detsgn,detlog
c
      detsgn = iflag
      detlog = 0.0D+00
      if (iflag .eq. 0)                 return
      index = 0
      indexp = 0
      do 2 i=1,nbloks
         nrow = integs(1,i)
         last = integs(3,i)
         do 1 k=1,last
            ip = index + nrow*(k-1) + ipivot(indexp+k)
            detlog = detlog + dlog(dabs(bloks(ip)))
    1       detsgn = detsgn*dsign(1.0D+00,bloks(ip))
         index = nrow*integs(2,i) + index
    2    indexp = indexp + nrow
                                        return
      end
      subroutine eqblok ( t, n, kpm,  work1, work2,
     &                 bloks, lenblk, integs, nbloks,  b )

c*********************************************************************72
c
cc EQBLOK is to be called in COLLOC.
c
c  from  * a practical guide to splines *  by c. de boor    
calls putit(difequ,bsplvd(bsplvb))
c  to be called in  c o l l o c
c
c******  i n p u t  ******
c  t   the knot sequence, of length n+kpm
c  n   the dimension of the approximating spline space, i.e., the order
c      of the linear system to be constructed.
c  kpm = k+m, the order of the approximating spline
c  lenblk   the maximum length of the array  bloks  as allowed by the
c           dimension statement in  colloc .
c
c******  w o r k   a r e a s  ******
c  work1    used in  putit, of size (kpm,kpm)
c  work2    used in  putit, of size (kpm,m+1)
c
c******  o u t p u t  ******
c  bloks    the coefficient matrix of the linear system, stored in al-
c           most block diagonal form, of size
c              kpm*sum(integs(1,i) , i=1,...,nbloks)
c  integs   an integer array, of size (3,nbloks), describing the block
c           structure.
c           integs(1,i)  =  number of rows in block  i
c           integs(2,i)  =  number of columns in block  i
c           integs(3,i)  =  number of elimination steps which can be
c                       carried out in block  i  before pivoting might
c                       bring in an equation from the next block.
c  nbloks   number of blocks, equals number of polynomial pieces
c  b   the right side of the linear system, stored corresponding to the
c      almost block diagonal form, of size sum(integs(1,i) , i=1,...,
c      nbloks).
c
c******  m e t h o d  ******
c  each breakpoint interval gives rise to a block in the linear system.
c  this block is determined by the  k  colloc.equations in the interval
c  with the side conditions (if any) in the interval interspersed ap-
c  propriately, and involves the  kpm  b-splines having the interval in
c  their support. correspondingly, such a block has  nrow = k + isidel
c  rows, with  isidel = number of side conditions in this and the prev-
c  ious intervals, and  ncol = kpm  columns.
c     further, because the interior knots have multiplicity  k, we can
c  carry out (in slvblk)  k  elimination steps in a block before pivot-
c  ing might involve an equation from the next block. in the last block,
c  of course, all kpm elimination steps will be carried out (in slvblk).
c
c  see the detailed comments in the solveblok package for further in-
c  formation about the almost block diagonal form used here.
c
      implicit none

      integer integs(3,1),kpm,lenblk,n,nbloks,   i,index,indexb,iside
     &                                      ,isidel,itermx,k,left,m,nrow
      double precision b(1),bloks(1),t(1),work1(1),work2(1),   rho,xside
      common /side/ m, iside, xside(10)
      common /other/ itermx,k,rho(19)
      index = 1
      indexb = 1
      i = 0
      iside = 1
      do 20 left=kpm,n,k
         i = i+1
c        determine integs(.,i)
         integs(2,i) = kpm
         if (left .lt. n)               go to 14
         integs(3,i) = kpm
         isidel = m
                                        go to 16
   14    integs(3,i) = k
c        at this point,  iside-1  gives the number of side conditions
c        incorporated so far. adding to this the side conditions in the
c        current interval gives the number  isidel .
         isidel = iside-1
   15    if (isidel .eq. m)             go to 16
         if (xside(isidel+1) .ge. t(left+1))
     &                                  go to 16
         isidel = isidel+1
                                        go to 15
   16    nrow = k + isidel
         integs(1,i) = nrow
c        the detailed equations for this block are generated and put
c        together in  p u t i t .
         if (lenblk .lt. index+nrow*kpm-1)go to 999
         call putit(t,kpm,left,work1,work2,bloks(index),nrow,b(indexb))
         index = index + nrow*kpm
   20    indexb = indexb + nrow
      nbloks = i
                                        return
  999 print 699,lenblk
  699 format(11h **********/23h the assigned dimension,i5
     &        ,38h for  bloks  in  colloc  is too small.)
                                        stop
      end
      subroutine evnnot ( break, coef, l, k, brknew, lnew, coefg )

c*********************************************************************72
c
cc EVNNOT is a version of NEWNOT returning uniform knots.
c
c  from  * a practical guide to splines *  by c. de boor    
c  this is a fake version of  n e w n o t  , of use in example 3 of     
c  chapter xiv .  
c  returns  lnew+1  knots in  brknew  which are equidistributed on (a,b)
c  = (break(1),break(l+1)) .  
c
      implicit none

      integer k,l,lnew,   i
      double precision break(1),brknew(1),coef(k,l),coefg(2,l),   step
c******  i n p u t ******
c  break, coef, l, k.....contains the pp-representation of a certain
c        function  f  of order  k . specifically,
c        d**(k-1)f(x) = coef(k,i)  for  break(i).le. x .lt.break(i+1)
c  lnew.....number of intervals into which the interval (a,b) is to be
c        sectioned by the new breakpoint sequence  brknew .
c
c******  o u t p u t  ******
c  brknew.....array of length  lnew+1  containing the new breakpoint se-
c        quence
c  coefg.....the coefficient part of the pp-repr.  break, coefg, l, 2
c        for the monotone p.linear function  g  wrto which  brknew  will
c        be equidistributed.
c
      brknew(1) = break(1)
      brknew(lnew+1) = break(l+1)
      step = (break(l+1) - break(1))/dble(lnew)
      do 93 i=2,lnew
   93    brknew(i) = break(1) + dble(i-1)*step
                                        return
      end
      subroutine factrb ( w, ipivot, d, nrow, ncol, last, iflag )

c*********************************************************************72
c
cc FACTRB constructs a partial PLU factorization.
c
c  adapted from p.132 of 'element.numer.analysis' by conte-de boor
c
c  constructs a partial plu factorization, corresponding to steps 1,...,
c   l a s t   in gauss elimination, for the matrix  w  of order
c   ( n r o w ,  n c o l ), using pivoting of scaled rows.
c
c  parameters
c    w       contains the (nrow,ncol) matrix to be partially factored
c            on input, and the partial factorization on output.
c    ipivot  an integer array of length nrow containing a record of the
c            pivoting strategy used; row ipivot(i) is used during the
c            i-th elimination step, i=1,...,last.
c    d       a work array of length nrow used to store row sizes
c            temporarily.
c    nrow    number of rows of w.
c    ncol    number of columns of w.
c    last    number of elimination steps to be carried out.
c    iflag   on output, equals iflag on input times (-1)**(number of
c            row interchanges during the factorization process), in
c            case no zero pivot was encountered.
c            otherwise, iflag = 0 on output.
c
      implicit none

      integer nrow

      integer ipivot(nrow),ncol,last,iflag, i,ipivi,ipivk,j,k,kp1
      double precision w(nrow,ncol),d(nrow),
     & awikdi,colmax,ratio,rowmax
c  initialize ipivot, d
      do 10 i=1,nrow
         ipivot(i) = i
         rowmax = 0.0D+00
         do 9 j=1,ncol
    9       rowmax = dmax1(rowmax, dabs(w(i,j)))
         if (rowmax .eq. 0.0D+00)            go to 999
   10    d(i) = rowmax
c gauss elimination with pivoting of scaled rows, loop over k=1,.,last
      k = 1
c        as pivot row for k-th step, pick among the rows not yet used,
c        i.e., from rows ipivot(k),...,ipivot(nrow), the one whose k-th
c        entry (compared to the row size) is largest. then, if this row
c        does not turn out to be row ipivot(k), redefine ipivot(k) ap-
c        propriately and record this interchange by changing the sign
c        of  i f l a g .
   11    ipivk = ipivot(k)
         if (k .eq. nrow)               go to 21
         j = k
         kp1 = k+1
         colmax = dabs(w(ipivk,k))/d(ipivk)
c              find the (relatively) largest pivot
         do 15 i=kp1,nrow
            ipivi = ipivot(i)
            awikdi = dabs(w(ipivi,k))/d(ipivi)
            if (awikdi .le. colmax)     go to 15
               colmax = awikdi
               j = i
   15       continue
         if (j .eq. k)                  go to 16
         ipivk = ipivot(j)
         ipivot(j) = ipivot(k)
         ipivot(k) = ipivk
         iflag = -iflag
   16    continue
c        if pivot element is too small in absolute value, declare
c        matrix to be noninvertible and quit.
         if (dabs(w(ipivk,k))+d(ipivk) .le. d(ipivk))
     &                                  go to 999
c        otherwise, subtract the appropriate multiple of the pivot
c        row from remaining rows, i.e., the rows ipivot(k+1),...,
c        ipivot(nrow), to make k-th entry zero. save the multiplier in
c        its place.
         do 20 i=kp1,nrow
            ipivi = ipivot(i)
            w(ipivi,k) = w(ipivi,k)/w(ipivk,k)
            ratio = -w(ipivi,k)
            do 20 j=kp1,ncol
   20          w(ipivi,j) = ratio*w(ipivk,j) + w(ipivi,j)
         k = kp1
c        check for having reached the next block.
         if (k .le. last)               go to 11
                                        return
c     if  last  .eq. nrow , check now that pivot element in last row
c     is nonzero.
   21 if( dabs(w(ipivk,nrow))+d(ipivk) .gt. d(ipivk) )
     *                                  return
c                   singularity flag set
  999 iflag = 0
                                        return
      end
      subroutine fcblok ( bloks, integs, nbloks, ipivot, scrtch, iflag )

c*********************************************************************72
c
cc FCBLOK supervises the PLU factorization of an almost block diagonal matrix.
c
calls subroutines  f a c t r b  and  s h i f t b .
c
c   f c b l o k  supervises the plu factorization with pivoting of
c  scaled rows of the almost block diagonal matrix stored in the arrays
c   b l o k s  and  i n t e g s .
c
c   factrb = subprogram which carries out steps 1,...,last of gauss
c            elimination (with pivoting) for an individual block.
c   shiftb = subprogram which shifts the remaining rows to the top of
c            the next block
c
c parameters
c    bloks   an array that initially contains the almost block diagonal
c            matrix  a  to be factored, and on return contains the com-
c            puted factorization of  a .
c    integs  an integer array describing the block structure of  a .
c    nbloks  the number of blocks in  a .
c    ipivot  an integer array of dimension  sum (integs(1,n) ; n=1,
c            ...,nbloks) which, on return, contains the pivoting stra-
c            tegy used.
c    scrtch  work area required, of length  max (integs(1,n) ; n=1,
c            ...,nbloks).
c    iflag   output parameter;
c            = 0  in case matrix was found to be singular.
c            otherwise,
c            = (-1)**(number of row interchanges during factorization)
c
      implicit none

      integer nbloks

      integer integs(3,nbloks),ipivot(1),iflag, i,index,indexb,indexn,
     &        last,ncol,nrow
      double precision bloks(1),scrtch(1)
      iflag = 1
      indexb = 1
      indexn = 1
      i = 1
c                        loop over the blocks.  i  is loop index
   10    index = indexn
         nrow = integs(1,i)
         ncol = integs(2,i)
         last = integs(3,i)
c        carry out elimination on the i-th block until next block
c        enters, i.e., for columns 1,...,last  of i-th block.
         call factrb(bloks(index),ipivot(indexb),scrtch,nrow,ncol,last,
     &            iflag)
c         check for having reached a singular block or the last block
         if (iflag .eq. 0 .or. i .eq. nbloks)
     &                                  return
         i = i+1
         indexn = nrow*ncol + index
c              put the rest of the i-th block onto the next block
         call shiftb(bloks(index),ipivot(indexb),nrow,ncol,last,
     &            bloks(indexn),integs(1,i),integs(2,i))
         indexb = indexb + nrow
                                        go to 10
      end
      subroutine interv ( xt, lxt, x, left, mflag )

c*********************************************************************72
c
cc INTERV brackets a real value in an ascending vector of values.
c
c  from  * a practical guide to splines *  by C. de Boor    
computes  left = max( i :  xt(i) .lt. xt(lxt) .and.  xt(i) .le. x )  .
c
c******  i n p u t  ******
c  xt.....a real sequence, of length  lxt , assumed to be nondecreasing
c  lxt.....number of terms in the sequence  xt .
c  x.....the point whose location with respect to the sequence  xt  is
c        to be determined.
c
c******  o u t p u t  ******
c  left, mflag.....both integers, whose value is
c
c   1     -1      if               x .lt.  xt(1)
c   i      0      if   xt(i)  .le. x .lt. xt(i+1)
c   i      0      if   xt(i)  .lt. x .eq. xt(i+1) .eq. xt(lxt)
c   i      1      if   xt(i)  .lt.        xt(i+1) .eq. xt(lxt) .lt. x
c
c        In particular,  mflag = 0  is the 'usual' case.  mflag .ne. 0
c        indicates that  x  lies outside the CLOSED interval
c        xt(1) .le. y .le. xt(lxt) . The asymmetric treatment of the
c        intervals is due to the decision to make all pp functions cont-
c        inuous from the right, but, by returning  mflag = 0  even if
C        x = xt(lxt), there is the option of having the computed pp function
c        continuous from the left at  xt(lxt) .
c
c******  m e t h o d  ******
c  The program is designed to be efficient in the common situation that
c  it is called repeatedly, with  x  taken from an increasing or decrea-
c  sing sequence. This will happen, e.g., when a pp function is to be
c  graphed. The first guess for  left  is therefore taken to be the val-
c  ue returned at the previous call and stored in the  l o c a l  varia-
c  ble  ilo . A first check ascertains that  ilo .lt. lxt (this is nec-
c  essary since the present call may have nothing to do with the previ-
c  ous call). Then, if  xt(ilo) .le. x .lt. xt(ilo+1), we set  left =
c  ilo  and are done after just three comparisons.
c     Otherwise, we repeatedly double the difference  istep = ihi - ilo
c  while also moving  ilo  and  ihi  in the direction of  x , until
c                      xt(ilo) .le. x .lt. xt(ihi) ,
c  after which we use bisection to get, in addition, ilo+1 = ihi .
c  left = ilo  is then returned.
c
      implicit none

      integer left,lxt,mflag,   ihi,ilo,istep,middle
      double precision x,xt(lxt)
      data ilo /1/
      save ilo  
      ihi = ilo + 1
      if (ihi .lt. lxt)                 go to 20
         if (x .ge. xt(lxt))            go to 110
         if (lxt .le. 1)                go to 90
         ilo = lxt - 1
         ihi = lxt
c
   20 if (x .ge. xt(ihi))               go to 40
      if (x .ge. xt(ilo))               go to 100
c
c              **** now x .lt. xt(ilo) . decrease  ilo  to capture  x .
      istep = 1
   31    ihi = ilo
         ilo = ihi - istep
         if (ilo .le. 1)                go to 35
         if (x .ge. xt(ilo))            go to 50
         istep = istep*2
                                        go to 31
   35 ilo = 1
      if (x .lt. xt(1))                 go to 90
                                        go to 50
c              **** now x .ge. xt(ihi) . increase  ihi  to capture  x .
   40 istep = 1
   41    ilo = ihi
         ihi = ilo + istep
         if (ihi .ge. lxt)              go to 45
         if (x .lt. xt(ihi))            go to 50
         istep = istep*2
                                        go to 41
   45 if (x .ge. xt(lxt))               go to 110
      ihi = lxt
c
c           **** now xt(ilo) .le. x .lt. xt(ihi) . narrow the interval.
   50 middle = (ilo + ihi)/2
      if (middle .eq. ilo)              go to 100
c     note. it is assumed that middle = ilo in case ihi = ilo+1 .
      if (x .lt. xt(middle))            go to 53
         ilo = middle
                                        go to 50
   53    ihi = middle
                                        go to 50
c**** set output and return.
   90 mflag = -1
      left = 1
                                        return
  100 mflag = 0
      left = ilo
                                        return
  110 mflag = 1
	  if (x .eq. xt(lxt)) mflag = 0
      left = lxt
  111 if (left .eq. 1)                  return
	  left = left - 1
	  if (xt(left) .lt. xt(lxt))        return
										go to 111
      end
      subroutine knots ( break, l, kpm, t, n )

c*********************************************************************72
c
cc KNOTS is to be called in COLLOC.
c
c  from  * a practical guide to splines *  by c. de boor    
c  to be called in  c o l l o c .
c  constructs from the given breakpoint sequence  b r e a k  the knot
c  sequence  t  so that
c  spline(k+m,t) = pp(k+m,break) with  m-1  continuous derivatives .
c  this means that
c  t(1),...,t(n+kpm)  =  break(1) kpm times, then break(2),...,
c        break(l) each  k  times, then, finally, break(l+1) kpm times.
c
c******  i n p u t  ******
c  break(1),...,break(l+1)  breakpoint sequence
c  l  number of intervals or pieces
c  kpm   = k + m, order of the pp function or spline
c
c******  o u t p u t  ******
c  t(1),...,t(n+kpm)  the knot sequence.
c  n     = l*k + m  =  dimension of  spline(k+m,t).
c
      implicit none

      integer l,kpm,n,   iside,j,jj,jjj,k,ll,m
      double precision break(1),t(1),   xside
      common /side/ m,iside,xside(10)
      k = kpm - m
      n = l*k + m
      jj = n + kpm
      jjj = l + 1
      do 11 ll=1,kpm
         t(jj) = break(jjj)
   11    jj = jj - 1
      do 12 j=1,l
      jjj = jjj - 1
         do 12 ll=1,k
            t(jj) = break(jjj)
   12       jj = jj - 1
      do 13 ll=1,kpm
   13    t(ll) = break(1)
                                        return
      end
      subroutine l2appr ( t, n, k, q, diag, bcoef )

c*********************************************************************72
c
cc L2APPR constructs a weighted L2 spline approximation to given data.
c
c  from  * a practical guide to splines *  by c. de boor    
c  to be called in main program  l 2 m a i n .
calls subprograms  bsplvb, bchfac/slv
c
constructs the (weighted discrete) l2-approximation by splines of order
c  k  with knot sequence  t(1), ..., t(n+k)  to given data points
c  ( tau(i), gtau(i) ), i=1,...,ntau. the b-spline coefficients
c  b c o e f   of the approximating spline are determined from the
c  normal equations using cholesky's method.
c
c******  i n p u t  ******
c  t(1), ..., t(n+k)  the knot sequence
c  n.....the dimension of the space of splines of order k with knots t.
c  k.....the order
c
c  w a r n i n g  . . .  the restriction   k .le. kmax (= 20)   is impo-
c        sed by the arbitrary dimension statement for  biatx  below, but
c        is  n o w h e r e   c h e c k e d   for.
c
c******  w o r k  a r r a y s  ******
c  q....a work array of size (at least) k*n. its first  k  rows are used
c       for the  k  lower diagonals of the gramian matrix  c .
c  diag.....a work array of length  n  used in bchfac .
c
c******  i n p u t  via  c o m m o n  /data/  ******
c  ntau.....number of data points
c  (tau(i),gtau(i)), i=1,...,ntau     are the  ntau  data points to be
c        fitted .
c  weight(i), i=1,...,ntau    are the corresponding weights .
c
c******  o u t p u t  ******
c  bcoef(1), ..., bcoef(n)  the b-spline coeffs. of the l2-appr.
c
c******  m e t h o d  ******
c  the b-spline coefficients of the l2-appr. are determined as the sol-
c  ution of the normal equations
c     sum ( (b(i),b(j))*bcoef(j) : j=1,...,n)  = (b(i),g),
c                                               i = 1, ..., n .
c  here,  b(i)  denotes the i-th b-spline,  g  denotes the function to
c  be approximated, and the  i n n e r   p r o d u c t  of two funct-
c  ions  f  and  g  is given by
c      (f,g)  :=  sum ( f(tau(i))*g(tau(i))*weight(i) : i=1,...,ntau) .
c  the arrays  t a u  and  w e i g h t  are given in common block
c   d a t a , as is the array  g t a u  containing the sequence
c  g(tau(i)), i=1,...,ntau.
c  the relevant function values of the b-splines  b(i), i=1,...,n, are
c  supplied by the subprogram  b s p l v b .
c     the coeff.matrix  c , with
c           c(i,j)  :=  (b(i), b(j)), i,j=1,...,n,
c  of the normal equations is symmetric and (2*k-1)-banded, therefore
c  can be specified by giving its k bands at or below the diagonal. for
c  i=1,...,n,  we store
c   (b(i),b(j))  =  c(i,j)  in  q(i-j+1,j), j=i,...,min0(i+k-1,n)
c  and the right side
c   (b(i), g )  in  bcoef(i) .
c  since b-spline values are most efficiently generated by finding sim-
c  ultaneously the value of  e v e r y  nonzero b-spline at one point,
c  the entries of  c  (i.e., of  q ), are generated by computing, for
c  each ll, all the terms involving  tau(ll)  simultaneously and adding
c  them to all relevant entries.
c
      implicit none

      integer k,n,   i,j,jj,kmax,left,leftmk,ll,mm,ntau,ntmax
      parameter (kmax=20,ntmax=200)
      double precision bcoef(n),diag(n),q(k,n),t(n+k),
     &  biatx(kmax),dw,gtau,tau,weight
      common / i4data / ntau
      common / r8data / tau(ntmax),gtau(ntmax),weight(ntmax)
c
      do 7 j=1,n
         bcoef(j) = 0.0D+00
         do 7 i=1,k
    7       q(i,j) = 0.0D+00
      left = k
      leftmk = 0
      do 20 ll=1,ntau
c                   locate  l e f t  s.t. tau(ll) in (t(left),t(left+1))
   10       if (left .eq. n)            go to 15
            if (tau(ll) .lt. t(left+1)) go to 15
            left = left+1
            leftmk = leftmk + 1
                                        go to 10
   15    call bsplvb ( t, k, 1, tau(ll), left, biatx )
c        biatx(mm) contains the value of b(left-k+mm) at tau(ll).
c        hence, with  dw := biatx(mm)*weight(ll), the number dw*gtau(ll)
c        is a summand in the inner product
c           (b(left-k+mm), g)  which goes into  bcoef(left-k+mm)
c        and the number biatx(jj)*dw is a summand in the inner product
c           (b(left-k+jj), b(left-k+mm)), into  q(jj-mm+1,left-k+mm)
c        since  (left-k+jj) - (left-k+mm) + 1  =  jj - mm + 1 .
         do 20 mm=1,k
            dw = biatx(mm)*weight(ll)
            j = leftmk + mm
            bcoef(j) = dw*gtau(ll) + bcoef(j)
            i = 1
            do 20 jj=mm,k
               q(i,j) = biatx(jj)*dw + q(i,j)
   20          i = i + 1
c
c             construct cholesky factorization for  c  in  q , then use
c             it to solve the normal equations
c                    c*x  =  bcoef
c             for  x , and store  x  in  bcoef .
      call bchfac ( q, k, n, diag )
      call bchslv ( q, k, n, bcoef )
                                        return
      end
      subroutine l2err ( prfun , ftau , error )

c*********************************************************************72
c
cc L2ERR computes the errors of an L2 approximation.
c
c  from  * a practical guide to splines *  by c. de boor    
c  this routine is to be called in the main program  l 2 m a i n .
calls subprogram  ppvalu(interv)
c  this subroutine computes various errors of the current l2-approxi-
c  mation , whose pp-repr. is contained in common block  approx  ,
c  to the given data contained in common block  data . it prints out
c  the average error  e r r l 1 , the l2-error  e r r l 2,  and the
c  maximum error  e r r m a x .
c
c******  i n p u t  ******
c  prfun  a logical value.  if prfun = TRUE, the routine prints out
c          the value of the approximation as well as its error at
c          every data point.
c
c******  o u t p u t  ******
c  ftau(1), ..., ftau(ntau),  with  ftau(i)  the approximation  f at
c          tau(i), all i.
c  error(1), ..., error(ntau),  with  error(i) = scale*(g - f)
c          at tau(i), all i. here,  s c a l e  equals  1. in case
c          prfun .ne. 'ON' , or the abs.error is greater than 100 some-
c          where. otherwise, s c a l e  is such that the maximum of
c          abs(error))  over all  i  lies between  10  and  100. this
c          makes the printed output more illustrative.
c
      implicit none

      logical prfun
      integer ie,k,l,ll,lpkmax,ltkmax,ntau,ntmax,on
      double precision ftau(ntau),error(ntau),
     &  break,coef,err,errmax,errl1,errl2,gtau,ppvalu,
     &  scale,tau,totalw,weight

      parameter (lpkmax=100,ntmax=200,ltkmax=2000)

      common / i4data / ntau
      common / r8data / tau(ntmax),gtau(ntmax),weight(ntmax),totalw
      common /approx/ break(lpkmax),coef(ltkmax),l,k

      errl1 = 0.0D+00
      errl2 = 0.0D+00
      errmax = 0.0D+00

      do ll=1,ntau
         ftau(ll) = ppvalu (break, coef, l, k, tau(ll), 0 )
         error(ll) = gtau(ll) - ftau(ll)
         err = dabs(error(ll))
         if (errmax .lt. err)   errmax = err
         errl1 = errl1 + err*weight(ll)
         errl2 = errl2 + err**2*weight(ll)
      end do

      errl1 = errl1/totalw
      errl2 = dsqrt(errl2/totalw)
      print 615,errl2,errl1,errmax
  615 format(///' least square error =',e20.6/
     &          ' average error      =',e20.6/
     &          ' maximum error      =',e20.6//)
c
c  Scale error curve and print.
c
      if ( prfun ) then

        ie = 0
        scale = 1.0D+00
        if (errmax .ge. 10.0D+00)              go to 18

        do ie=1,9
           scale = scale*10.0D+00
           if (errmax*scale .ge. 10.0D+00)     go to 18
        end do

   18   do 19 ll=1,ntau
   19    error(ll) = error(ll)*scale

        print 620,ie,(ll,tau(ll),ftau(ll),error(ll),ll=1,ntau)
  620 format(///14x,'approximation and scaled error curve'/7x,
     &'data point',7x,'approximation',3x,'deviation x 10**',i1/
     &(i4,f16.8,f16.8,f17.6))

      end if

      return
      end
      subroutine l2knts ( break, l, k, t, n )

c*********************************************************************72
c
cc L2KNTS converts breakpoints to knots.
c
c  from  * a practical guide to splines *  by c. de boor    
c  to be called in main program  l 2 m a i n .
converts the breakpoint sequence  b r e a k   into a corresponding knot
c  sequence  t  to allow the repr. of a pp function of order  k  with
c  k-2 continuous derivatives as a spline of order  k  with knot
c  sequence  t . this means that
c  t(1), ..., t(n+k) =  break(1) k times, then break(i), i=2,...,l, each
c                       once, then break(l+1) k times .
c  therefore,  n = k-1 + l.
c
c******  i n p u t  ******
c  k     the order
c  l     the number of polynomial pieces
c  break(1), ...,break(l+1)  the breakpoint sequence
c
c******  o u t p u t  ******
c  t(1),...,t(n+k)   the knot sequence
c  n     the dimension of the corresp. spline space of order  k .
c
      implicit none

      integer k,l,n,   i,km1
      double precision break(1),t(1)
c     dimension break(l+1),t(n+k)
      km1 = k - 1
      do 5 i=1,km1
    5    t(i) = break(1)
      do 6 i=1,l
    6    t(km1+i) = break(i)
      n = km1 + l
      do 7 i=1,k
    7    t(n+i) = break(l+1)
                                        return
      end
      subroutine newnot ( break, coef, l, k, brknew, lnew, coefg )

c*********************************************************************72
c
cc NEWNOT returns LNEW+1 knots which are equidistributed on (A,B)
c
c  from  * a practical guide to splines *  by c. de boor    
c  returns  lnew+1  knots in  brknew  which are equidistributed on (a,b)
c  = (break(1),break(l+1)) wrto a certain monotone fctn  g  related to
c  the k-th root of the k-th derivative of the pp function  f  whose pp-
c  representation is contained in  break, coef, l, k .
c
c******  i n p u t ******
c  break, coef, l, k.....contains the pp-representation of a certain
c        function  f  of order  k . Specifically,
c        d**(k-1)f(x) = coef(k,i)  for  break(i).le. x .lt.break(i+1)
c  lnew.....number of intervals into which the interval (a,b) is to be
c        sectioned by the new breakpoint sequence  brknew .
c
c******  o u t p u t  ******
c  brknew.....array of length  lnew+1  containing the new breakpoint se-
c        quence
c  coefg.....the coefficient part of the pp-repr.  break, coefg, l, 2
c        for the monotone p.linear function  g  wrto which  brknew  will
c        be equidistributed.
c
c******  optional  p r i n t e d  o u t p u t  ******
c  coefg.....the pp coeffs of  g  are printed out if  iprint  is set
c        .gt. 0  in data statement below.
c
c******  m e t h o d  ******
c     The k-th derivative of the given pp function  f  does not exist
c  (except perhaps as a linear combination of delta functions). Never-
c  theless, we construct a p.constant function  h  with breakpoint se-
c  quence  break  which is approximately proportional to abs(d**k(f)).
c  Specifically, on  (break(i), break(i+1)),
c
c     abs(jump at break(i) of pc)    abs(jump at break(i+1) of pc)
c h = --------------------------  +  ----------------------------
c       break(i+1) - break(i-1)         break(i+2) - break(i)
c
c  with  pc  the p.constant (k-1)st derivative of  f .
c      Then, the p.linear function  g  is constructed as
c
c    g(x)  =  integral of  h(y)**(1/k)  for  y  from  a  to  x
c
c  and its pp coeffs. stored in  coefg .
c     then  brknew  is determined by
c
c        brknew(i)  =  a + g**(-1)((i-1)*step) , i=1,...,lnew+1
c
c  where  step = g(b)/lnew  and  (a,b) = (break(1),break(l+1)) .
c     In the event that  pc = d**(k-1)(f) is constant in  (a,b)  and
c  therefore  h = 0 identically,  brknew  is chosen uniformly spaced.
c
      implicit none

      integer k,l,lnew,   i,iprint,j
      double precision break(1),brknew(1),coef(k,l),coefg(2,l),
     &   dif,difprv,oneovk
     &                                               ,step,stepi
c     dimension break(l+1), brknew(lnew+1)
current fortran standard makes it impossible to specify the dimension
c  of  break   and  brknew  without the introduction of additional
c  otherwise superfluous arguments.
      data iprint /0/
c
      brknew(1) = break(1)
      brknew(lnew+1) = break(l+1)
c                               if  g  is constant,  brknew  is uniform.
      if (l .le. 1)                     go to 90
c                        construct the continuous p.linear function  g .
      oneovk = 1.0D+00/dble(k)
      coefg(1,1) = 0.0D+00
      difprv = dabs(coef(k,2) - coef(k,1))/(break(3)-break(1))
      do 10 i=2,l
         dif = dabs(coef(k,i) - coef(k,i-1))/(break(i+1) - break(i-1))
         coefg(2,i-1) = (dif + difprv)**oneovk
         coefg(1,i) = coefg(1,i-1)+coefg(2,i-1)*(break(i)-break(i-1))
   10    difprv = dif
      coefg(2,l) = (2.0D+00*difprv)**oneovk
c                                                     step  =  g(b)/lnew
      step = (coefg(1,l)+coefg(2,l)*(break(l+1)-break(l)))/dble(lnew)
c
      if (iprint .gt. 0) print 600, step,(i,coefg(1,i),coefg(2,i),i=1,l)
  600 format(7h step =,e16.7/(i5,2e16.5))
c                              if  g  is constant,  brknew  is uniform .
      if (step .le. 0.0D+00)                 go to 90
c
c      For i=2,...,lnew, construct  brknew(i) = a + g**(-1)(stepi),
c      with  stepi = (i-1)*step .  this requires inversion of the p.lin-
c      ear function  g .  For this,  j  is found so that
c         g(break(j)) .le. stepi .le. g(break(j+1))
c      and then
c         brknew(i)  =  break(j) + (stepi-g(break(j)))/dg(break(j)) .
c      The midpoint is chosen if  dg(break(j)) = 0 .
      j = 1
      do 30 i=2,lnew
         stepi = dble(i-1)*step
   21       if (j .eq. l)               go to 27
            if (stepi .le. coefg(1,j+1))go to 27
            j = j + 1
                                        go to 21
   27    if (coefg(2,j) .eq. 0.0D+00)        go to 29
            brknew(i) = break(j) + (stepi - coefg(1,j))/coefg(2,j)
                                        go to 30
   29       brknew(i) = (break(j) + break(j+1))/2.
   30    continue
                                        return
c
c                              if  g  is constant,  brknew  is uniform .
   90 step = (break(l+1) - break(1))/dble(lnew)
      do 93 i=2,lnew
   93    brknew(i) = break(1) + dble(i-1)*step
                                        return
      end
      double precision function ppvalu (break, coef, l, k, x, jderiv )

c*********************************************************************72
c
cc PPVALU evaluates a piecewise polynomial function or its derivative.
c
c  from  * a practical guide to splines *  by c. de boor    
calls  interv
calculates value at  x  of  jderiv-th derivative of pp fct from pp-repr
c
c******  i n p u t  ******
c  break, coef, l, k.....forms the pp-representation of the function  f
c        to be evaluated. specifically, the j-th derivative of  f  is
c        given by
c
c     (d**j)f(x) = coef(j+1,i) + h*(coef(j+2,i) + h*( ... (coef(k-1,i) +
c                             + h*coef(k,i)/(k-j-1))/(k-j-2) ... )/2)/1
c
c        with  h = x - break(i),  and
c
c       i  =  max( 1 , max( j ,  break(j) .le. x , 1 .le. j .le. l ) ).
c
c  x.....the point at which to evaluate.
c  jderiv.....integer giving the order of the derivative to be evaluat-
c        ed.  a s s u m e d  to be zero or positive.
c
c******  o u t p u t  ******
c  ppvalu.....the value of the (jderiv)-th derivative of  f  at  x.
c
c******  m e t h o d  ******
c     the interval index  i , appropriate for  x , is found through a
c  call to  interv . the formula above for the  jderiv-th derivative
c  of  f  is then evaluated (by nested multiplication).
c
      implicit none

      integer jderiv,k,l,   i,m,ndummy
      double precision break(l+1),coef(k,l),x,   fmmjdr,h
      ppvalu = 0.0D+00
      fmmjdr = k - jderiv
c              derivatives of order  k  or higher are identically zero.
      if (fmmjdr .le. 0.0D+00)               go to 99
c
c              find index  i  of largest breakpoint to the left of  x .
      call interv ( break, l+1, x, i, ndummy )
c
c      Evaluate  jderiv-th derivative of  i-th polynomial piece at  x .
      h = x - break(i)
      m = k 
    9    ppvalu = (ppvalu/fmmjdr)*h + coef(m,i)
         m = m - 1
         fmmjdr = fmmjdr - 1.0D+00
         if (fmmjdr .gt. 0.0D+00)            go to 9 
   99                                   return
      end
      subroutine putit ( t, kpm, left, scrtch, dbiatx, q, nrow, b )

c*********************************************************************72
c
cc PUTIT puts together one block of the collocation equation system.
c
c  from  * a practical guide to splines *  by c. de Boor (7 may 92)    
calls  bsplvd(bsplvb),difequ(*)
c  to be called by  e q b l o k .
c
c  puts together one block of the collocation equation system
c
c******  i n p u t  ******
c  t     knot sequence, of size  left+kpm (at least)
c  kpm   order of spline
c  left  integer indicating interval of interest, viz the interval
c           (t(left), t(left+1))
c  nrow  number of rows in block to be put together
c
c******  w o r k  a r e a  ******
c  scrtch   used in bsplvd, of size (kpm,kpm)
c  dbiatx   used to contain derivatives of b-splines, of size (kpm,m+1)
c           with dbiatx(j,i+1) containing the i-th derivative of the
c           j-th b-spline of interest
c
c******  o u t p u t  ******
c  q  the block, of size (nrow,kpm)
c  b  the corresponding piece of the right side, of size (nrow)
c
c******  m e t h o d  ******
c  the  k  collocation equations for the interval  (t(left),t(left+1))
c  are constructed with the aid of the subroutine  d i f e q u ( 2, .,
c  . ) and interspersed (in order) with the side conditions (if any) in
c  this interval, using  d i f e q u ( 3, ., . )  for the information.
c     the block  q  has  kpm  columns, corresponding to the  kpm  b-
c  splines of order  kpm  which have the interval (t(left),t(left+1))
c  in their support. the block's diagonal is part of the diagonal of the
c  total system. the first equation in this block not overlapped by the
c  preceding block is therefore equation  l o w r o w , with lowrow = 1+
c  number of side conditions in preceding intervals (or blocks).
c
      implicit none

      integer kpm,left,nrow,   i,irow,iside,itermx,j,k,ll,lowrow,m,mode
     &                        ,mp1
      double precision b(1),dbiatx(kpm,1),q(nrow,kpm),scrtch(1),t(1),
     &   dx,rho,sum
     &                                                ,v(20),xm,xside,xx
      common /side/ m, iside, xside(10)
      common /other/ itermx,k,rho(19)
      mp1 = m+1
      do 10 j=1,kpm
         do 10 i=1,nrow
   10       q(i,j) = 0.0D+00
      xm = (t(left+1)+t(left))/2.0D+00
      dx = (t(left+1)-t(left))/2.0D+00
c
      ll = 1
      lowrow = iside
      do 30 irow=lowrow,nrow
         if (ll .gt. k)                 go to 22
         mode = 2
c        next collocation point is ...
         xx = xm + dx*rho(ll)
         ll = ll + 1
c        the corresp.collocation equation is next unless the next side
c        condition occurs at a point at, or to the left of, the next
c        collocation point.
         if (iside .gt. m)              go to 24
         if (xside(iside) .gt. xx)      go to 24
         ll = ll - 1
   22    mode = 3
         xx = xside(iside)
   24    call difequ ( mode, xx, v )
c        the next equation, a collocation equation (mode = 2) or a side
c        condition (mode = 3), reads
c   (*)   (v(m+1)*d**m + v(m)*d**(m-1) +...+ v(1)*d**0)f(xx) = v(m+2)
c        in terms of the info supplied by  difequ . the corresponding
c        equation for the b-coeffs of  f  therefore has the left side of
c        (*), evaluated at each of the  kpm  b-splines having  xx  in
c        their support, as its  kpm  possibly nonzero coefficients.
         call bsplvd ( t, kpm, xx, left, scrtch, dbiatx, mp1 )
         do 26 j=1,kpm
            sum = 0.0D+00
            do 25 i=1,mp1
   25          sum = v(i)*dbiatx(j,i) + sum
   26       q(irow,j) = sum
   30    b(irow) = v(m+2)
                                        return
      end
      double precision function round ( x )

c*********************************************************************72
c
cc ROUND is called to add some noise to data.
c
c  from  * a practical guide to splines *  by c. de boor    
called in example 1  of chapter xiii
c
      implicit none

      double precision x,   flip,size
      common /rount/ size
      data flip /-1.0D+00/
      flip = -flip
      round = x + flip*size
                                        return
      end
      subroutine sbblok ( bloks, integs, nbloks, ipivot, b, x )

c*********************************************************************72
c
cc SBBLOK solves a linear system that was factored by FCBLOK.
c
calls subroutines  s u b f o r  and  s u b b a k .
c
c  supervises the solution (by forward and backward substitution) of
c  the linear system  a*x = b  for x, with the plu factorization of  a
c  already generated in  f c b l o k .  individual blocks of equations
c  are solved via  s u b f o r  and  s u b b a k .
c
c parameters
c    bloks, integs, nbloks, ipivot    are as on return from fcblok.
c    b       the right side, stored corresponding to the storage of
c            the equations. see comments in  s l v b l k  for details.
c    x       solution vector
c
      implicit none

      integer nbloks

      integer integs(3,nbloks),ipivot(1), i,index,indexb,indexx,j,last,
     &        nbp1,ncol,nrow
      double precision bloks(1),b(1),x(1)
c
c      forward substitution pass
c
      index = 1
      indexb = 1
      indexx = 1
      do 20 i=1,nbloks
         nrow = integs(1,i)
         last = integs(3,i)
         call subfor(bloks(index),ipivot(indexb),nrow,last,b(indexb),
     &               x(indexx))
         index = nrow*integs(2,i) + index
         indexb = indexb + nrow
   20    indexx = indexx + last
c
c     back substitution pass
c
      nbp1 = nbloks + 1
      do 30 j=1,nbloks
         i = nbp1 - j
         nrow = integs(1,i)
         ncol = integs(2,i)
         last = integs(3,i)
         index = index - nrow*ncol
         indexb = indexb - nrow
         indexx = indexx - last
   30    call subbak(bloks(index),ipivot(indexb),nrow,ncol,last,
     &               x(indexx))
                                        return
      end
      subroutine setupq ( x, dx, y, npoint, v, qty )

c*********************************************************************72
c
cc SETUPQ is to be called in SMOOTH.
c
c  from  * a practical guide to splines *  by c. de boor    
c  to be called in  s m o o t h
c  put  delx = x(.+1) - x(.)  into  v(.,4),
c  put  the three bands of  q-transp*d  into  v(.,1-3), and
c  put the three bands of  (d*q)-transp*(d*q)  at and above the diagonal
c     into  v(.,5-7) .
c     here,  q is  the tridiagonal matrix of order (npoint-2,npoint)
c  with general row  1/delx(i) , -1/delx(i) - 1/delx(i+1) , 1/delx(i+1)
c  and   d  is the diagonal matrix  with general row  dx(i) .
c
      implicit none

      integer npoint,   i,npm1
      double precision dx(npoint),qty(npoint),v(npoint,7),
     &  x(npoint),y(npoint),
     &                                                       diff,prev
      npm1 = npoint - 1
      v(1,4) = x(2) - x(1)
      do 11 i=2,npm1
         v(i,4) = x(i+1) - x(i)
         v(i,1) = dx(i-1)/v(i-1,4)
         v(i,2) = - dx(i)/v(i,4) - dx(i)/v(i-1,4)
   11    v(i,3) = dx(i+1)/v(i,4)
       v(npoint,1) = 0.0D+00
      do 12 i=2,npm1
   12    v(i,5) = v(i,1)**2 + v(i,2)**2 + v(i,3)**2
      if (npm1 .lt. 3)                  go to 14
      do 13 i=3,npm1
   13    v(i-1,6) = v(i-1,2)*v(i,1) + v(i-1,3)*v(i,2)
   14 v(npm1,6) = 0.0D+00
      if (npm1 .lt. 4)                  go to 16
      do 15 i=4,npm1
   15    v(i-2,7) = v(i-2,3)*v(i,1)
   16 v(npm1-1,7) = 0.0D+00
      v(npm1,7) = 0.0D+00
construct  q-transp. * y  in  qty.
      prev = (y(2) - y(1))/v(1,4)
      do 21 i=2,npm1
         diff = (y(i+1)-y(i))/v(i,4)
         qty(i) = diff - prev
   21    prev = diff
                                        return
      end
      subroutine shiftb ( ai, ipivot, nrowi, ncoli, last,
     &                    ai1, nrowi1, ncoli1 )

c*********************************************************************72
c
cc SHIFTB shifts the rows in the current block.
c
c  shifts the rows in current block, ai, not used as pivot rows, if
c  any, i.e., rows ipivot(last+1),...,ipivot(nrowi), onto the first
c  mmax = nrow-last rows of the next block, ai1, with column last+j of
c  ai  going to column j , j=1,...,jmax=ncoli-last. the remaining col-
c  umns of these rows of ai1 are zeroed out.
c
c                             picture
c
c       original situation after         results in a new block i+1
c       last = 2 columns have been       created and ready to be
c       done in factrb (assuming no      factored by next factrb call.
c       interchanges of rows)
c                   1
c              x  x 1x  x  x           x  x  x  x  x
c                   1
c              0  x 1x  x  x           0  x  x  x  x
c  block i          1                       ---------------
c  nrowi = 4   0  0 1x  x  x           0  0 1x  x  x  0  01
c  ncoli = 5        1                       1             1
c  last = 2    0  0 1x  x  x           0  0 1x  x  x  0  01
c  -------------------------------          1             1   new
c                   1x  x  x  x  x          1x  x  x  x  x1  block
c                   1                       1             1   i+1
c  block i+1        1x  x  x  x  x          1x  x  x  x  x1
c  nrowi1= 5        1                       1             1
c  ncoli1= 5        1x  x  x  x  x          1x  x  x  x  x1
c  -------------------------------          1-------------1
c                   1
c
      implicit none

      integer ncoli
      integer ncoli1
      integer nrowi
      integer nrowi1

      integer ipivot(nrowi),last, ip,j,jmax,jmaxp1,m,mmax
      double precision ai(nrowi,ncoli),ai1(nrowi1,ncoli1)

      mmax = nrowi - last
      jmax = ncoli - last
      if (mmax .lt. 1 .or. jmax .lt. 1) return
c              put the remainder of block i into ai1
      do 10 m=1,mmax
         ip = ipivot(last+m)
         do 10 j=1,jmax
   10       ai1(m,j) = ai(ip,last+j)
      if (jmax .eq. ncoli1)             return
c              zero out the upper right corner of ai1
      jmaxp1 = jmax + 1
      do 20 j=jmaxp1,ncoli1
         do 20 m=1,mmax
   20       ai1(m,j) = 0.0D+00
                                        return
      end
      subroutine slvblk ( bloks, integs, nbloks, b, ipivot, x, iflag )

c*********************************************************************72
c
cc SLVBLK solves the almost block diagonal linear system A * x = b. 
c
c    this program solves  the  linear system  a*x = b  where a is an
c  almost block diagonal matrix.  such almost block diagonal matrices
c  arise naturally in piecewise polynomial interpolation or approx-
c  imation and in finite element methods for two-point boundary value
c  problems.  the plu factorization method is implemented here to take
c  advantage of the special structure of such systems for savings in
c  computing time and storage requirements.
c
c                  parameters
c  bloks   a one-dimenional array, of length
c                   sum( integs(1,i)*integs(2,i) ; i = 1,nbloks )
c          on input, contains the blocks of the almost block diagonal
c          matrix  a  .  the array integs (see below and the example)
c          describes the block structure.
c          on output, contains correspondingly the plu factorization
c          of  a  (if iflag .ne. 0).  certain of the entries into bloks
c          are arbitrary (where the blocks overlap).
c  integs  integer array description of the block structure of  a .
c            integs(1,i) = no. of rows of block i        =  nrow
c            integs(2,i) = no. of colums of block i      =  ncol
c            integs(3,i) = no. of elim. steps in block i =  last
c                          i  = 1,2,...,nbloks
c          the linear system is of order
c                n  =  sum ( integs(3,i) , i=1,...,nbloks ),
c          but the total number of rows in the blocks is
c              nbrows = sum( integs(1,i) ; i = 1,...,nbloks)
c  nbloks  number of blocks
c  b       right side of the linear system, array of length nbrows.
c          certain of the entries are arbitrary, corresponding to
c          rows of the blocks which overlap (see block structure and
c          the example below).
c  ipivot  on output, integer array containing the pivoting sequence
c          used. length is nbrows
c  x       on output, contains the computed solution (if iflag .ne. 0)
c          length is n.
c  iflag   on output, integer
c            = (-1)**(no. of interchanges during factorization)
c                   if  a  is invertible
c            = 0    if  a  is singular
c
c                   auxiliary programs
c  fcblok (bloks,integs,nbloks,ipivot,scrtch,iflag)  factors the matrix
c           a , and is used for this purpose in slvblk. its arguments
c          are as in slvblk, except for
c              scrtch = a work array of length max(integs(1,i)).
c
c  sbblok (bloks,integs,nbloks,ipivot,b,x)  solves the system a*x = b
c          once  a  is factored. this is done automatically by slvblk
c          for one right side b, but subsequent solutions may be
c          obtained for additional b-vectors. the arguments are all
c          as in slvblk.
c
c  dtblok (bloks,integs,nbloks,ipivot,iflag,detsgn,detlog) computes the
c          determinant of  a  once slvblk or fcblok has done the fact-
c          orization.the first five arguments are as in slvblk.
c              detsgn  = sign of the determinant
c              detlog  = natural log of the determinant
c
c             ------ block structure of  a  ------
c  the nbloks blocks are stored consecutively in the array  bloks .
c  the first block has its (1,1)-entry at bloks(1), and, if the i-th
c  block has its (1,1)-entry at bloks(index(i)), then
c         index(i+1) = index(i)  +  nrow(i)*ncol(i) .
c    the blocks are pieced together to give the interesting part of  a
c  as follows.  for i = 1,2,...,nbloks-1, the (1,1)-entry of the next
c  block (the (i+1)st block ) corresponds to the (last+1,last+1)-entry
c  of the current i-th block.  recall last = integs(3,i) and note that
c  this means that
c      a. every block starts on the diagonal of  a .
c      b. the blocks overlap (usually). the rows of the (i+1)st block
c         which are overlapped by the i-th block may be arbitrarily de-
c         fined initially. they are overwritten during elimination.
c    the right side for the equations in the i-th block are stored cor-
c  respondingly as the last entries of a piece of  b  of length  nrow
c  (= integs(1,i)) and following immediately in  b  the corresponding
c  piece for the right side of the preceding block, with the right side
c  for the first block starting at  b(1) . in this, the right side for
c  an equation need only be specified once on input, in the first block
c  in which the equation appears.
c
c             ------ example and test driver ------
c    the test driver for this package contains an example, a linear
c  system of order 11, whose nonzero entries are indicated in the fol-
c  lowing schema by their row and column index modulo 10. next to it
c  are the contents of the  integs  arrray when the matrix is taken to
c  be almost block diagonal with  nbloks = 5, and below it are the five
c  blocks.
c
c                      nrow1 = 3, ncol1 = 4
c           11 12 13 14
c           21 22 23 24   nrow2 = 3, ncol2 = 3
c           31 32 33 34
c  last1 = 2      43 44 45
c                 53 54 55            nrow3 = 3, ncol3 = 4
c        last2 = 3         66 67 68 69   nrow4 = 3, ncol4 = 4
c                          76 77 78 79      nrow5 = 4, ncol5 = 4
c                          86 87 88 89
c                 last3 = 1   97 98 99 90
c                    last4 = 1   08 09 00 01
c                                18 19 10 11
c                       last5 = 4
c
c         actual input to bloks shown by rows of blocks of  a .
c      (the ** items are arbitrary, this storage is used by slvblk)
c
c  11 12 13 14  / ** ** **  / 66 67 68 69  / ** ** ** **  / ** ** ** **
c  21 22 23 24 /  43 44 45 /  76 77 78 79 /  ** ** ** ** /  ** ** ** **
c  31 32 33 34/   53 54 55/   86 87 88 89/   97 98 99 90/   08 09 00 01
c                                                           18 19 10 11
c
c  index = 1      index = 13  index = 22     index = 34     index = 46
c
c         actual right side values with ** for arbitrary values
c  b1 b2 b3 ** b4 b5 b6 b7 b8 ** ** b9 ** ** b10 b11
c
c  (it would have been more efficient to combine block 3 with block 4)
c
      implicit none

      integer nbloks

      integer integs(3,nbloks),ipivot(1),iflag
      double precision bloks(1),b(1),x(1)
c     in the call to fcblok,  x  is used for temporary storage.
      call fcblok(bloks,integs,nbloks,ipivot,x,iflag)
      if (iflag .eq. 0)                 return
      call sbblok(bloks,integs,nbloks,ipivot,b,x)
                                        return
      end
      double precision function smooth ( x, y, dy, npoint, s, v, a )

c*********************************************************************72
c
cc SMOOTH constructs the cubic smoothing spline to given data.
c
c  from  * a practical guide to splines *  by c. de boor    
calls  setupq, chol1d
c
c  constructs the cubic smoothing spline  f  to given data  (x(i),y(i)),
c  i=1,...,npoint, which has as small a second derivative as possible
c  while
c  s(f) = sum( ((y(i)-f(x(i)))/dy(i))**2 , i=1,...,npoint ) .le. s .
c
c******  i n p u t  ******
c  x(1),...,x(npoint)   data abscissae,  a s s u m e d  to be strictly
c        increasing .
c  y(1),...,y(npoint)     corresponding data ordinates .
c  dy(1),...,dy(npoint)     estimate of uncertainty in data,  a s s u m-
c        e d  to be positive .
c  npoint.....number of data points,  a s s u m e d  .gt. 1
c  s.....upper bound on the discrete weighted mean square distance of
c        the approximation  f  from the data .
c
c******  w o r k  a r r a y s  *****
c  v.....of size (npoint,7)
c  a.....of size (npoint,4)
c
c*****  o u t p u t  *****
c  a(.,1).....contains the sequence of smoothed ordinates .
c  a(i,j) = f^(j-1)(x(i)), j=2,3,4, i=1,...,npoint-1 ,  i.e., the
c        first three derivatives of the smoothing spline  f  at the
c        left end of each of the data intervals .
c     w a r n i n g . . .   a  would have to be transposed before it
c        could be used in  ppvalu .
c
c******  m e t h o d  ******
c     The matrices  Q-transp*d  and  Q-transp*D**2*Q  are constructed in
c   s e t u p q  from  x  and  dy , as is the vector  qty = Q-transp*y .
c  Then, for given  p , the vector  u  is determined in  c h o l 1 d  as
c  the solution of the linear system
c               (6(1-p)Q-transp*D**2*Q + p*Q)u  = qty  .
c  From  u , the smoothing spline  f  (for this choice of smoothing par-
c  ameter  p ) is obtained in the sense that
c                        f(x(.))  =  y - 6(1-p)D**2*Q*u        and
c                      f''(x(.))  =  6*p*u                      .
c     The smoothing parameter  p  is found (if possible) so that
c                sf(p)  =  s ,
c  with  sf(p) = s(f) , where  f  is the smoothing spline as it depends
c  on  p .  if  s = 0, then p = 1 . if  sf(0) .le. s , then p = 0 .
c  Otherwise, the secant method is used to locate an appropriate  p  in
c  the open interval  (0,1) . However, straightforward application of
c  the secant method, as done in the original version of this program,
c  can be very slow and is influenced by the units in which  x  and  y 
c  are measured, as C. Reinsch has pointed out. Instead, on recommend-
c  ation from C. Reinsch, the secant method is applied to the function
c           g:q |--> 1/sqrt{sfq(q)} - 1/sqrt{s} ,
c  with  sfq(q) := sf(q/(1+q)), since  1/sqrt{sfq}  is monotone increasing
c  and close to linear for larger  q . One starts at  q = 0  with a
c  Newton step, i.e., 
c                q_0 = 0,  q_1 = -g(0)/g'(0)
c  with  g'(0) = -(1/2) sfq(0)^{-3/2} dsfq, where dsfq = -12*u-transp*r*u ,
c  and  u  as obtained for  p = 0 . Iteration terminates as soon as 
c   abs(sf - s) .le. .01*s .
c
      implicit none
c
c     logical test 
c     parameter (test = .true.)
c     integer itercnt
      integer npoint,   i,npm1
      double precision a(npoint,4),dy(npoint),s,v(npoint,7),
     &  x(npoint),y(npoint)
     &     ,change,ooss,oosf,p,prevsf,prevq,q,sfq,sixp,six1mp,utru
      call setupq(x,dy,y,npoint,v,a(1,4))
      if (s .gt. 0.0D+00)                    go to 20
   10 p = 1.0D+00
      call chol1d(p,v,a(1,4),npoint,1,a(1,3),a(1,1))
      sfq = 0.0D+00
                                        go to 60
   20 p = 0.0D+00
      call chol1d(p,v,a(1,4),npoint,1,a(1,3),a(1,1))
      sfq = 0.0D+00
      do 21 i=1,npoint
   21    sfq = sfq + (a(i,1)*dy(i))**2
      sfq = sfq*36.0D+00
      if (sfq .le. s)                   go to 60
      utru = 0.0D+00
      do 25 i=2,npoint
   25    utru = utru + v(i-1,4)*(a(i-1,3)*(a(i-1,3)+a(i,3))+a(i,3)**2)
      ooss = 1.0D+00/dsqrt(s)
      oosf = 1.0D+00/dsqrt(sfq)
      q = -(oosf-ooss)*sfq/(6.0D+00*utru*oosf)
c  secant iteration for the determination of p starts here.
c     itercnt = 0
      prevq = 0.0D+00
      prevsf = oosf
   30 call chol1d(q/(1.0D+00+q),v,a(1,4),npoint,1,a(1,3),a(1,1))
      sfq = 0.0D+00
      do 35 i=1,npoint
   35    sfq = sfq + (a(i,1)*dy(i))**2
      sfq = sfq*36.0D+00/(1.0D+00+q)**2
      if (dabs(sfq-s) .le. 0.01D+00*s)        go to 59
      oosf = 1.0D+00/dsqrt(sfq)
      change = (q-prevq)/(oosf-prevsf)*(oosf-ooss)
      prevq = q
      q = q - change
      prevsf = oosf
c     itercnt = itercnt + 1
                                        go to 30
   59 p = q/(1.0D+00+q)
correct value of p has been found.
compute pol.coefficients from  Q*u (in a(.,1)).
   60 smooth = sfq
c     if (test) then
c        print *, 'number of iterations = ', itercnt
c     end if
c     six1mp = 6./(1.0D+00+q)
      do 61 i=1,npoint
c  61    a(i,1) = y(i) - six1mp*dy(i)**2*a(i,1)
   61    a(i,1) = y(i) - 6.0D+00 * ( 1.0D+00 - p ) *dy(i)**2*a(i,1)
      sixp = 6.0D+00*p
      do 62 i=1,npoint
   62    a(i,3) = a(i,3)*sixp
      npm1 = npoint - 1
      do 63 i=1,npm1
         a(i,4) = (a(i+1,3)-a(i,3))/v(i,4)
   63    a(i,2) = (a(i+1,1)-a(i,1))/v(i,4)
     &    - (a(i,3)+a(i,4)/3.0D+00*v(i,4))/2.0D+00*v(i,4)
                                        return
      end
      subroutine spli2d ( tau, gtau, t, n, k, m, work, q, bcoef, iflag )********    

c*********************************************************************72
c
cc SPLI2D produces a interpolatory tensor product spline.
c
c  from  * a practical guide to splines *  by c. de boor    
calls bsplvb, banfac/slv
c  this is an extended version of  splint , for the use in tensor prod- --------    
c  uct interpolation.                                                   --------    
c     
c   spli2d  produces the b-spline coeff.s  bcoef(j,.)  of the spline of --------    
c   order  k  with knots  t (i), i=1,..., n + k , which takes on the    --------    
c   value  gtau (i,j)  at  tau (i), i=1,..., n , j=1,..., m .           --------    
c     
c******  i n p u t  ******    
c  tau.....array of length  n , containing data point abscissae.  
c    a s s u m p t i o n . . .  tau  is strictly increasing 
c  gtau(.,j)..corresponding array of length  n , containing data point  --------    
c        ordinates, j=1,...,m                                           --------    
c  t.....knot sequence, of length  n+k    
c  n.....number of data points and dimension of spline space  s(k,t)    
c  k.....order of spline
c  m.....number of data sets                                            ********    
c     
c******  w o r k   a r e a  ******                                      ********    
c  work  a vector of length  n                                          ********    
c     
c******  o u t p u t  ******  
c  q.....array of size  (2*k-1)*n , containing the triangular factoriz- 
c        ation of the coefficient matrix of the linear system for the b-
c        coefficients of the spline interpolant.
c           the b-coeffs for the interpolant of an additional data set  
c        (tau(i),htau(i)), i=1,...,n  with the same data abscissae can  
c        be obtained without going through all the calculations in this 
c        routine, simply by loading  htau  into  bcoef  and then execut-
c        ing the    call banslv ( q, 2*k-1, n, k-1, k-1, bcoef )  
c  bcoef.....the b-coefficients of the interpolant, of length  n  
c  iflag.....an integer indicating success (= 1)  or failure (= 2)
c        the linear system to be solved is (theoretically) invertible if
c        and only if    
c              t(i) .lt. tau(i) .lt. tau(i+k),    all i.    
c        violation of this condition is certain to lead to  iflag = 2 . 
c     
c******  m e t h o d  ******  
c     the i-th equation of the linear system  a*bcoef = b  for the b-co-
c  effs of the interpolant enforces interpolation at  tau(i), i=1,...,n.
c  hence,  b(i) = gtau(i), all i, and  a  is a band matrix with  2k-1   
c   bands (if it is invertible).    
c     the matrix  a  is generated row by row and stored, diagonal by di-
c  agonal, in the  r o w s  of the array  q , with the main diagonal go-
c  ing into row  k .  see comments in the program below.    
c     the banded system is then solved by a call to  banfac (which con- 
c  structs the triangular factorization for  a  and stores it again in  
c   q ), followed by a call to  banslv (which then obtains the solution 
c   bcoef  by substitution).  
c     banfac  does no pivoting, since the total positivity of the matrix
c  a  makes this unnecessary. 
c     
      implicit none

      integer iflag,k,m,n,   i,ilp1mx,j,jj,km1,kpkm2,left,lenq,np1      ********    
      double precision bcoef(m,n),gtau(n,m),q(1),
     &  t(1),tau(n),work(n),   taui
c     dimension q(2*k-1,n), t(n+k)  
current fortran standard makes it impossible to specify precisely the   
c  dimension of  q  and  t  without the introduction of otherwise super-
c  fluous additional arguments.     
      np1 = n + 1 
      km1 = k - 1 
      kpkm2 = 2*km1     
      left = k    
c                zero out all entries of q
      lenq = n*(k+km1)  
      do 5 i=1,lenq     
    5    q(i) = 0.0D+00
c     
c  ***   loop over i to construct the  n  interpolation equations 
      do 30 i=1,n 
         taui = tau(i)  
         ilp1mx = min0(i+k,np1)     
c        *** find  left  in the closed interval (i,i+k-1) such that     
c                t(left) .le. tau(i) .lt. t(left+1)   
c        matrix is singular if this is not possible   
         left = max0(left,i)  
         if (taui .lt. t(left))         go to 998     
   15       if (taui .lt. t(left+1))    go to 16
            left = left + 1   
            if (left .lt. ilp1mx)       go to 15
         left = left - 1
         if (taui .gt. t(left+1))       go to 998     
c        *** the i-th equation enforces interpolation at taui, hence    
c        a(i,j) = b(j,k,t)(taui), all j. only the  k  entries with  j = 
c        left-k+1,...,left actually might be nonzero. these  k  numbers 
c        are returned, in  work  (used for temp.storage here), by the   --------    
c        following
   16    call bsplvb ( t, k, 1, taui, left, work )                      --------    
c        we therefore want  work(j) = b(left -k+j)(taui) to go into     --------    
c        a(i,left-k+j), i.e., into  q(i-(left+j)+2*k,(left+j)-k) since  
c        a(i+j,j)  is to go into  q(i+k,j), all i,j,  if we consider  q 
c        as a two-dim. array , with  2*k-1  rows (see comments in 
c        banfac). in the present program, we treat  q  as an equivalent 
c        one-dimensional array (because of fortran restrictions on
c        dimension statements) . we therefore want  work(j) to go into  --------    
c        entry    
c            i -(left+j) + 2*k + ((left+j) - k-1)*(2*k-1)   
c                   =  i-left+1 + (left -k)*(2*k-1) + (2*k-2)*j   
c        of  q .  
         jj = i-left+1 + (left-k)*(k+km1) 
         do 30 j=1,k    
            jj = jj+kpkm2     
   30       q(jj) = work(j)                                             --------    
c     
c     ***obtain factorization of  a  , stored again in  q.  
      call banfac ( q, k+km1, n, km1, km1, iflag )    
                                        go to (40,999), iflag     
c     *** solve  a*bcoef = gtau  by backsubstitution  
   40 do 50 j=1,m                                                       ********    
         do 41 i=1,n                                                    --------    
   41       work(i) = gtau(i,j)                                         ********    
         call banslv ( q, k+km1, n, km1, km1, work )                    --------    
         do 50 i=1,n                                                    ********    
   50       bcoef(j,i) = work(i)                                        ********    
                                        return  
  998 iflag = 2   
  999 print 699   
  699 format(41h linear system in  splint  not invertible)  
                                        return  
      end   
      subroutine splint ( tau, gtau, t, n, k, q, bcoef, iflag )

c*********************************************************************72
c
cc SPLINT produces the B-spline coefficients BCOEF of an interpolating spline.
c
c  from  * a practical guide to splines *  by c. de boor    
calls bsplvb, banfac/slv
c
c   splint  produces the b-spline coeff.s  bcoef  of the spline of order
c   k  with knots  t(i), i=1,..., n + k , which takes on the value
c   gtau(i) at  tau(i), i=1,..., n .
c
c******  i n p u t  ******
c  tau.....array of length  n , containing data point abscissae.
c    a s s u m p t i o n . . .  tau  is strictly increasing
c  gtau.....corresponding array of length  n , containing data point or-
c        dinates
c  t.....knot sequence, of length  n+k
c  n.....number of data points and dimension of spline space  s(k,t)
c  k.....order of spline
c
c******  o u t p u t  ******
c  q.....array of size  (2*k-1)*n , containing the triangular factoriz-
c        ation of the coefficient matrix of the linear system for the b-
c        coefficients of the spline interpolant.
c           the b-coeffs for the interpolant of an additional data set
c        (tau(i),htau(i)), i=1,...,n  with the same data abscissae can
c        be obtained without going through all the calculations in this
c        routine, simply by loading  htau  into  bcoef  and then execut-
c        ing the    call banslv ( q, 2*k-1, n, k-1, k-1, bcoef )
c  bcoef.....the b-coefficients of the interpolant, of length  n
c  iflag.....an integer indicating success (= 1)  or failure (= 2)
c        the linear system to be solved is (theoretically) invertible if
c        and only if
c              t(i) .lt. tau(i) .lt. t(i+k),    all i.
c        violation of this condition is certain to lead to  iflag = 2 .
c
c******  m e t h o d  ******
c     the i-th equation of the linear system  a*bcoef = b  for the b-co-
c  effs of the interpolant enforces interpolation at  tau(i), i=1,...,n.
c  hence,  b(i) = gtau(i), all i, and  a  is a band matrix with  2k-1
c   bands (if it is invertible).
c     the matrix  a  is generated row by row and stored, diagonal by di-
c  agonal, in the  r o w s  of the array  q , with the main diagonal go-
c  ing into row  k .  see comments in the program below.
c     the banded system is then solved by a call to  banfac (which con-
c  structs the triangular factorization for  a  and stores it again in
c   q ), followed by a call to  banslv (which then obtains the solution
c   bcoef  by substitution).
c     banfac  does no pivoting, since the total positivity of the matrix
c  a  makes this unnecessary.
c
      implicit none

      integer iflag,k,n,   i,ilp1mx,j,jj,km1,kpkm2,left,lenq,np1
      double precision bcoef(n),gtau(n),q(1),t(1),tau(n),   taui
c     dimension q(2*k-1,n), t(n+k)
current fortran standard makes it impossible to specify precisely the
c  dimension of  q  and  t  without the introduction of otherwise super-
c  fluous additional arguments.
      np1 = n + 1
      km1 = k - 1
      kpkm2 = 2*km1
      left = k
c                zero out all entries of q
      lenq = n*(k+km1)
      do 5 i=1,lenq
    5    q(i) = 0.0D+00
c
c  ***   loop over i to construct the  n  interpolation equations
      do 30 i=1,n
         taui = tau(i)
         ilp1mx = min0(i+k,np1)
c        *** find  left  in the closed interval (i,i+k-1) such that
c                t(left) .le. tau(i) .lt. t(left+1)
c        matrix is singular if this is not possible
         left = max0(left,i)
         if (taui .lt. t(left))         go to 998
   15       if (taui .lt. t(left+1))    go to 16
            left = left + 1
            if (left .lt. ilp1mx)       go to 15
         left = left - 1
         if (taui .gt. t(left+1))       go to 998
c        *** the i-th equation enforces interpolation at taui, hence
c        a(i,j) = b(j,k,t)(taui), all j. only the  k  entries with  j =
c        left-k+1,...,left actually might be nonzero. these  k  numbers
c        are returned, in  bcoef (used for temp.storage here), by the
c        following
   16    call bsplvb ( t, k, 1, taui, left, bcoef )
c        we therefore want  bcoef(j) = b(left-k+j)(taui) to go into
c        a(i,left-k+j), i.e., into  q(i-(left+j)+2*k,(left+j)-k) since
c        a(i+j,j)  is to go into  q(i+k,j), all i,j,  if we consider  q
c        as a two-dim. array , with  2*k-1  rows (see comments in
c        banfac). in the present program, we treat  q  as an equivalent
c        one-dimensional array (because of fortran restrictions on
c        dimension statements) . we therefore want  bcoef(j) to go into
c        entry
c            i -(left+j) + 2*k + ((left+j) - k-1)*(2*k-1)
c                   =  i-left+1 + (left -k)*(2*k-1) + (2*k-2)*j
c        of  q .
         jj = i-left+1 + (left-k)*(k+km1)
         do 30 j=1,k
            jj = jj+kpkm2
   30       q(jj) = bcoef(j)
c
c     ***obtain factorization of  a  , stored again in  q.
      call banfac ( q, k+km1, n, km1, km1, iflag )
                                        go to (40,999), iflag
c     *** solve  a*bcoef = gtau  by backsubstitution
   40 do 41 i=1,n
   41    bcoef(i) = gtau(i)
      call banslv ( q, k+km1, n, km1, km1, bcoef )
                                        return
  998 iflag = 2
  999 print 699
  699 format(41h linear system in  splint  not invertible)
                                        return
      end
      subroutine splopt ( tau, n, k, scrtch, t, iflag )

c*********************************************************************72
c
cc SPLOPT computes the knots for an optimal recovery scheme. 
c
c  from  * a practical guide to splines *  by c. de boor    
calls bsplvb, banfac/slv
computes the knots  t  for the optimal recovery scheme of order  k
c  for data at  tau(i), i=1,...,n .
c
c******  i n p u t  ******
c  tau.....array of length  n , containing the interpolation points.
c     a s s u m e d  to be nondecreasing, with tau(i).lt.tau(i+k),all i.
c  n.....number of data points .
c  k.....order of the optimal recovery scheme to be used .
c
c******  w o r k  a r r a y  *****
c  scrtch.....array of length  (n-k)(2k+3) + 5k + 3 . the various
c        contents are specified in the text below .
c
c******  o u t p u t  ******
c  iflag.....integer indicating success (=1) or failure (=2) .
c     if iflag = 1, then
c  t.....array of length  n+k  containing the optimal knots ready for
c        use in optimal recovery. specifically,  t(1) = ... = t(k) =
c        tau(1)  and  t(n+1) = ... = t(n+k) = tau(n) , while the  n-k
c        interior knots  t(k+1), ..., t(n)  are calculated as described
c        below under  *method* .
c     if iflag = 2, then
c        k .lt. 3, or n .lt. k, or a certain linear system was found to
c        be singular.
c
c******  p r i n t e d  o u t p u t  ******
c  a comment will be printed in case  iflag = 2  or newton iterations
c  failed to converge in  n e w t m x  iterations .
c
c******  m e t h o d  ******
c     the (interior) knots  t(k+1), ..., t(n)  are determined by newtons
c  method in such a way that the signum function which changes sign at
c   t(k+1), ..., t(n)  and nowhere else in  (tau(1),tau(n)) is orthogon-
c  al to the spline space  spline( k , tau )  on that interval .
c     let  xi(j)  be the current guess for  t(k+j), j=1,...,n-k. then
c  the next newton iterate is of the form
c              xi(j)  +  (-)**(n-k-j)*x(j)  ,  j=1,...,n-k,
c  with  x  the solution of the linear system
c                        c*x  =  d  .
c  here,  c(i,j) = b(i)(xi(j)), all j, with  b(i)  the i-th b-spline of
c  order  k  for the knot sequence  tau , all i, and  d  is the vector
c  given by  d(i) = sum( -a(j) , j=i,...,n )*(tau(i+k)-tau(i))/k, all i,
c  with  a(i) = sum ( (-)**(n-k-j)*b(i,k+1,tau)(xi(j)) , j=1,...,n-k )
c  for i=1,...,n-1, and  a(n) = -.5 .
c     (see chapter  xiii  of text and references there for a derivation)
c     the first guess for  t(k+j)  is  (tau(j+1)+...+tau(j+k-1))/(k-1) .
c     iteration terminates if  max(abs(x(j))) .lt. t o l  , with
c                 t o l  =  t o l r t e *(tau(n)-tau(1))/(n-k) ,
c  or else after  n e w t m x  iterations , currently,
c                 newtmx, tolrte / 10, .000001
c
      implicit none

      integer iflag,k,n,   i,id,index,j,km1,kpk,kpkm1,kpn,kp1,l,left
     &,leftmk,lenw,ll,llmax,llmin,na,nb,nc,nd,newtmx,newton,nmk,nmkm1,nx
      double precision scrtch(1),t(1),tau(n),
     &   del,delmax,floatk,sign,signst,sum
     &                             ,tol,tolrte,xij
c     dimension scrtch((n-k)*(2*k+3)+5*k+3), t(n+k)
current fortran standard makes it impossible to specify the precise dim-
c  ensions of  scrtch  and  t  without the introduction of otherwise
c  superfluous additional arguments .
      data newtmx,tolrte / 10,0.000001D+00/
      nmk = n-k
      if (nmk)                          1,56,2
    1 print 601,n,k
  601 format(13h argument n =,i4,29h in  splopt  is less than k =,i3)
                                        go to 999
    2 if (k .gt. 2)                     go to 3
      print 602,k
  602 format(13h argument k =,i3,27h in  splopt  is less than 3)
                                        go to 999
   3  nmkm1 = nmk - 1
      floatk = k
      kpk = k+k
      kp1 = k+1
      km1 = k-1
      kpkm1 = kpk-1
      kpn = k+n
      signst = -1.0D+00
      if (nmk .gt. (nmk/2)*2)  signst = 1.0D+00
c  scrtch(i) = tau-extended(i), i=1,...,n+k+k
      nx = n+kpk+1
c  scrtch(i+nx) = xi(i),i=0,...,n-k+1
      na = nx + nmk + 1
c  scrtch(i+na) = -a(i), i=1,...,n
      nd = na + n
c  scrtch(i+nd) = x(i) or d(i), i=1,...,n-k
      nb = nd + nmk
c  scrtch(i+nb) = biatx(i),i=1,...,k+1
      nc = nb + kp1
c  scrtch(i+(j-1)*(2k-1) + nc) = w(i,j) = c(i-k+j,j), i=j-k,...,j+k,
c                                                     j=1,...,n-k.
      lenw = kpkm1*nmk
c  extend  tau  to a knot sequence and store in scrtch.
      do 5 j=1,k
         scrtch(j) = tau(1)
    5    scrtch(kpn+j) = tau(n)
      do 6 j=1,n
    6    scrtch(k+j) = tau(j)
c  first guess for  scrtch (.+nx)  =  xi .
      scrtch(nx) = tau(1)
      scrtch(nmk+1+nx) = tau(n)
      do 10 j=1,nmk
         sum = 0.0D+00
         do 9 l=1,km1
    9       sum = sum + tau(j+l)
   10    scrtch(j+nx) = sum/dble(km1)
c  last entry of  scrtch (.+na)  =  - a  is always ...
      scrtch(n+na) = .5
c  start newton iteration.
      newton = 1
      tol = tolrte*(tau(n) - tau(1))/dble(nmk)
c  start newton step
compute the 2k-1 bands of the matrix c and store in scrtch(.+nc),
c  and compute the vector  scrtch(.+na) = -a.
   20 do 21 i=1,lenw
   21    scrtch(i+nc) = 0.0D+00
      do 22 i=2,n
   22    scrtch(i-1+na) = 0.0D+00
      sign = signst
      left = kp1
      do 28 j=1,nmk
         xij = scrtch(j+nx)
   23       if (xij .lt. scrtch(left+1))go to 25
            left = left + 1
            if (left .lt. kpn)          go to 23
            left = left - 1
   25    call bsplvb(scrtch,k,1,xij,left,scrtch(1+nb))
c        the tau sequence in scrtch is preceded by  k  additional knots
c        therefore,  scrtch(ll+nb)  now contains  b(left-2k+ll)(xij)
c        which is destined for  c(left-2k+ll,j), and therefore for
c            w(left-k-j+ll,j)= scrtch(left-k-j+ll + (j-1)*kpkm1 + nc)
c        since we store the 2k-1 bands of  c  in the 2k-1  r o w s  of
c        the work array w, and  w  in turn is stored in  s c r t c h ,
c        with  w(1,1) = scrtch(1 + nc) .
c            also, c  being of order  n-k, we would want  1 .le.
c        left-2k+ll .le. n-k  or
c           llmin = 2k-left  .le.  ll  .le.  n-left+k = llmax .
         leftmk = left-k
         index = leftmk-j + (j-1)*kpkm1 + nc
         llmin = max0(1,k-leftmk)
         llmax = min0(k,n-leftmk)
         do 26 ll=llmin,llmax
   26       scrtch(ll+index) = scrtch(ll+nb)
         call bsplvb(scrtch,kp1,2,xij,left,scrtch(1+nb))
         id = max0(0,leftmk-kp1)
         llmin = 1 - min0(0,leftmk-kp1)
         do 27 ll=llmin,kp1
            id = id + 1
   27        scrtch(id+na) = scrtch(id+na) - sign*scrtch(ll+nb)
   28    sign = -sign
      call banfac(scrtch(1+nc),kpkm1,nmk,km1,km1,iflag)
                                        go to (45,44),iflag
   44 print 644
  644 format(32h c in  splopt  is not invertible)
                                        return
compute  scrtch (.+nd) =  d  from  scrtch (.+na) = - a .
   45 i=n   
   46    scrtch(i-1+na) = scrtch(i-1+na) + scrtch(i+na)
         i = i-1  
         if (i .gt. 1)                  go to 46
      do 49 i=1,nmk
   49    scrtch(i+nd) = scrtch(i+na)*(tau(i+k)-tau(i))/floatk
compute  scrtch (.+nd) =  x .
      call banslv(scrtch(1+nc),kpkm1,nmk,km1,km1,scrtch(1+nd))
compute  scrtch (.+nd) = change in  xi . modify, if necessary, to
c  prevent new  xi  from moving more than 1/3 of the way to its
c  neighbors. then add to  xi  to obtain new  xi  in scrtch(.+nx).
      delmax = 0.0D+00
      sign = signst
      do 53 i=1,nmk
         del = sign*scrtch(i+nd)
         delmax = dmax1(delmax,dabs(del))
         if (del .gt. 0.0D+00)               go to 51
         del = dmax1(del,(scrtch(i-1+nx)-scrtch(i+nx))/3.0D+00)
                                        go to 52
   51    del = dmin1(del,(scrtch(i+1+nx)-scrtch(i+nx))/3.0D+00)
   52    sign = -sign
   53    scrtch(i+nx) = scrtch(i+nx) + del
call it a day in case change in  xi  was small enough or too many
c  steps were taken.
      if (delmax .lt. tol)              go to 54
      newton = newton + 1
      if (newton .le. newtmx)           go to 20
      print 653,newtmx
  653 format(33h no convergence in  splopt  after,i3,14h newton steps.)
   54 do 55 i=1,nmk
   55    t(k+i) = scrtch(i+nx)
   56 do 57 i=1,k
         t(i) = tau(1)
   57    t(n+i) = tau(n)
                                        return
  999 iflag = 2
                                        return
      end
      subroutine subbak ( w, ipivot, nrow, ncol, last, x )

c*********************************************************************72
c
cc SUBBAK carries out back substitution for the current block.
c
c  carries out backsubstitution for current block.
c
c parameters
c    w, ipivot, nrow, ncol, last  are as on return from factrb.
c    x(1),...,x(ncol)  contains, on input, the right side for the
c            equations in this block after backsubstitution has been
c            carried up to but not including equation ipivot(last).
c            means that x(j) contains the right side of equation ipi-
c            vot(j) as modified during elimination, j=1,...,last, while
c            for j .gt. last, x(j) is already a component of the solut-
c            ion vector.
c    x(1),...,x(ncol) contains, on output, the components of the solut-
c            ion corresponding to the present block.
c
      implicit none

      integer ncol
      integer nrow

      integer ipivot(nrow),last,  ip,j,k,kp1
      double precision w(nrow,ncol),x(ncol), sum
      k = last
      ip = ipivot(k)
      sum = 0.0D+00
      if (k .eq. ncol)                  go to 4
      kp1 = k+1
    2    do 3 j=kp1,ncol
    3       sum = w(ip,j)*x(j) + sum
    4    x(k) = (x(k) - sum)/w(ip,k)
         if (k .eq. 1)                  return
         kp1 = k
         k = k-1
         ip = ipivot(k)
         sum = 0.0D+00
                                        go to 2
      end
      subroutine subfor ( w, ipivot, nrow, last, b, x )

c*********************************************************************72
c
cc SUBFOR carries out the forward pass of substitution for the current block.
c
c  carries out the forward pass of substitution for the current block,
c  i.e., the action on the right side corresponding to the elimination
c  carried out in  f a c t r b  for this block.
c     at the end, x(j) contains the right side of the transformed
c  ipivot(j)-th equation in this block, j=1,...,nrow. then, since
c  for i=1,...,nrow-last, b(nrow+i) is going to be used as the right
c  side of equation  i  in the next block (shifted over there from
c  this block during factorization), it is set equal to x(last+i) here.
c
c parameters
c    w, ipivot, nrow, last  are as on return from factrb.
c    b(j)   is expected to contain, on input, the right side of j-th
c           equation for this block, j=1,...,nrow.
c    b(nrow+j)   contains, on output, the appropriately modified right
c           side for equation j in next block, j=1,...,nrow-last.
c    x(j)   contains, on output, the appropriately modified right
c           side of equation ipivot(j) in this block, j=1,...,last (and
c           even for j=last+1,...,nrow).
c
      implicit none

      integer last
      integer nrow

      integer ipivot(nrow), ip,j,jmax,k,lastp1,nrowml
c     dimension b(nrow + nrow-last)
      double precision w(nrow,last),b(1),sum,x(nrow)
      ip = ipivot(1)
      x(1) = b(ip)
      if (nrow .eq. 1)                  go to 99
      do 15 k=2,nrow
         ip = ipivot(k)
         jmax = amin0(k-1,last)
         sum = 0.0D+00
         do 14 j=1,jmax
   14       sum = w(ip,j)*x(j) + sum
   15    x(k) = b(ip) - sum
c
c     transfer modified right sides of equations ipivot(last+1),...,
c     ipivot(nrow) to next block.
      nrowml = nrow - last
      if (nrowml .eq. 0)                go to 99
      lastp1 = last+1
      do 25 k=lastp1,nrow
   25    b(nrowml+k) = x(k)
   99                                   return
      end
      subroutine tautsp ( tau, gtau, ntau, gamma, s,
     &                    break, coef, l, k, iflag )

c*********************************************************************72
c
cc TAUTSP constructs a cubic spline interpolant to given data.
c
c  from  * a practical guide to splines *  by c. de boor    
constructs cubic spline interpolant to given data
c         tau(i), gtau(i), i=1,...,ntau.
c  if  gamma .gt. 0., additional knots are introduced where needed to
c  make the interpolant more flexible locally. this avoids extraneous
c  inflection points typical of cubic spline interpolation at knots to
c  rapidly changing data.
c
c  parameters
c            input
c  tau      sequence of data points. must be strictly increasing.
c  gtau     corresponding sequence of function values.
c  ntau     number of data points. must be at least  4 .
c  gamma    indicates whether additional flexibility is desired.
c          = 0., no additional knots
c          in (0.,3.), under certain conditions on the given data at
c                points i-1, i, i+1, and i+2, a knot is added in the
c                i-th interval, i=2,...,ntau-2. see description of meth-
c                od below. the interpolant gets rounded with increasing
c                gamma. a value of  2.5  for gamma is typical.
c          in (3.,6.), same , except that knots might also be added in
c                intervals in which an inflection point would be permit-
c                ted.  a value of  5.5  for gamma is typical.
c            output
c  break, coef, l, k  give the pp-representation of the interpolant.
c          specifically, for break(i) .le. x .le. break(i+1), the
c        interpolant has the form
c  f(x) = coef(1,i) +dx(coef(2,i) +(dx/2)(coef(3,i) +(dx/3)coef(4,i)))
c        with  dx = x - break(i) and i=1,...,l .
c  iflag   = 1, ok
c          = 2, input was incorrect. a printout specifying the mistake
c            was made.
c            workspace
c  s     is required, of size (ntau,6). the individual columns of this
c        array contain the following quantities mentioned in the write-
c        up and below.
c     s(.,1) = dtau = tau(.+1) - tau
c     s(.,2) = diag = diagonal in linear system
c     s(.,3) = u = upper diagonal in linear system
c     s(.,4) = r = right side for linear system (initially)
c            = fsecnd = solution of linear system , namely the second
c                       derivatives of interpolant at  tau
c     s(.,5) = z = indicator of additional knots
c     s(.,6) = 1/hsecnd(1,x) with x = z or = 1-z. see below.
c
c  ------  m e t h o d  ------
c  on the i-th interval, (tau(i), tau(i+1)), the interpolant is of the
c  form
c  (*)  f(u(x)) = a + b*u + c*h(u,z) + d*h(1-u,1-z) ,
c  with  u = u(x) = (x - tau(i))/dtau(i). here,
c       z = z(i) = addg(i+1)/(addg(i) + addg(i+1))
c  (= .5, in case the denominator vanishes). with
c       addg(j) = abs(ddg(j)), ddg(j) = dg(j+1) - dg(j),
c       dg(j) = divdif(j) = (gtau(j+1) - gtau(j))/dtau(j)
c  and
c       h(u,z) = alpha*u**3 + (1 - alpha)*(max(((u-zeta)/(1-zeta)),0)**3
c  with
c       alpha(z) = (1-gamma/3)/zeta
c       zeta(z) = 1 - gamma*min((1 - z), 1/3)
c  thus, for 1/3 .le. z .le. 2/3,  f  is just a cubic polynomial on
c  the interval i. otherwise, it has one additional knot, at
c         tau(i) + zeta*dtau(i) .
c  as  z  approaches  1, h(.,z) has an increasingly sharp bend  near 1,
c  thus allowing  f  to turn rapidly near the additional knot.
c     in terms of f(j) = gtau(j) and
c       fsecnd(j) = 2.derivative of  f  at  tau(j),
c  the coefficients for (*) are given as
c       a = f(i) - d
c       b = (f(i+1) - f(i)) - (c - d)
c       c = fsecnd(i+1)*dtau(i)**2/hsecnd(1,z)
c       d = fsecnd(i)*dtau(i)**2/hsecnd(1,1-z)
c  hence can be computed once fsecnd(i),i=1,...,ntau, is fixed.
c   f  is automatically continuous and has a continuous second derivat-
c  ive (except when z = 0 or 1 for some i). we determine fscnd(.) from
c  the requirement that also the first derivative of  f  be contin-
c  uous. in addition, we require that the third derivative be continuous
c  across  tau(2) and across  tau(ntau-1) . this leads to a strictly
c  diagonally dominant tridiagonal linear system for the fsecnd(i)
c  which we solve by gauss elimination without pivoting.
c
      implicit none

      integer iflag,k,l,ntau,   i,method,ntaum1
      double precision alph,break(1),coef(4,1),gamma,
     &  gtau(ntau),s(ntau,6),tau(ntau)
     &    ,alpha,c,d,del,denom,divdif,entry,entry3,factor,factr2,gam
     &    ,onemg3,onemzt,ratio,sixth,temp,x,z,zeta,zt2
      alph(x) = dmin1(1.0D+00,onemg3/x)
c
c  there must be at least  4  interpolation points.
      if (ntau .ge. 4)                  go to 5
      print 600,ntau
  600 format(8h ntau = ,i4,20h. should be .ge. 4 .)
                                        go to 999
c
construct delta tau and first and second (divided) differences of data
c
    5 ntaum1 = ntau - 1
      do 6 i=1,ntaum1
         s(i,1) = tau(i+1) - tau(i)
         if (s(i,1) .gt. 0.0D+00)            go to 6
         print 610,i,tau(i),tau(i+1)
  610    format(7h point ,i3,13h and the next,2e15.6,15h are disordered)
                                        go to 999
    6    s(i+1,4) = (gtau(i+1)-gtau(i))/s(i,1)
      do 7 i=2,ntaum1
    7    s(i,4) = s(i+1,4) - s(i,4)
c
construct system of equations for second derivatives at  tau. at each
c  interior data point, there is one continuity equation, at the first
c  and the last interior data point there is an additional one for a
c  total of  ntau  equations in  ntau  unknowns.
c
      i = 2
      s(2,2) = s(1,1)/3.0D+00
      sixth = 1.0D+00/6.0D+00
      method = 2
      gam = gamma
      if (gam .le. 0.0D+00)   method = 1
      if ( gam .le. 3.0D+00)                 go to 9
      method = 3
      gam = gam - 3.0D+00
    9 onemg3 = 1.0D+00 - gam/3.0D+00
c                 ------ loop over i ------
   10 continue
c          construct z(i) and zeta(i)
      z = 0.5D+00
                                        go to (19,11,12),method
   11 if (s(i,4)*s(i+1,4) .lt. 0.0D+00)      go to 19
   12 temp = dabs(s(i+1,4))
      denom = dabs(s(i,4)) + temp
      if (denom .eq. 0.0D+00)                go to 19
      z = temp/denom
      if (dabs(z - 0.5D+00) .le. sixth)  z = 0.5D+00
   19 s(i,5) = z
c   ******set up part of the i-th equation which depends on
c         the i-th interval
      if (z - 0.5D+00)                       21,22,23
   21 zeta = gam*z
      onemzt = 1.0D+00 - zeta
      zt2 = zeta**2
      alpha = alph(onemzt)
      factor = zeta/(alpha*(zt2-1.0D+00) + 1.0D+00)
      s(i,6) = zeta*factor/6.0D+00
      s(i,2) = s(i,2) + s(i,1)
     &  *((1.0D+00-alpha*onemzt)*factor/2.0D+00 - s(i,6))
c     if z = 0 and the previous z = 1, then d(i) = 0. since then
c     also u(i-1) = l(i+1) = 0, its value does not matter. reset
c     d(i) = 1 to insure nonzero pivot in elimination.
      if (s(i,2) .le. 0.0D+00) s(i,2) = 1.0D+00
      s(i,3) = s(i,1)/6.0D+00
                                        go to 25
   22 s(i,2) = s(i,2) + s(i,1)/3.0D+00
      s(i,3) = s(i,1)/6.0D+00
                                        go to 25
   23 onemzt = gam*(1.0D+00 - z)
      zeta = 1.0D+00 - onemzt
      alpha = alph(zeta)
      factor = onemzt/(1.0D+00 - alpha*zeta*(1.0D+00+onemzt))
      s(i,6) = onemzt*factor/6.0D+00
      s(i,2) = s(i,2) + s(i,1)/3.
      s(i,3) = s(i,6)*s(i,1)
   25 if (i .gt. 2)                     go to 30
      s(1,5) = 0.5D+00
c  ******the first two equations enforce continuity of the first and of
c        the third derivative across tau(2).
      s(1,2) = s(1,1)/6.0D+00
      s(1,3) = s(2,2)
      entry3 = s(2,3)
      if (z - 0.5D+00)                       26,27,28
   26 factr2 = zeta*(alpha*(zt2-1.0D+00) 
     &  + 1.0D+00)/(alpha*(zeta*zt2-1.0D+00)+1.0D+00)
      ratio = factr2*s(2,1)/s(1,2)
      s(2,2) = factr2*s(2,1) + s(1,1)
      s(2,3) = -factr2*s(1,1)
                                        go to 29
   27 ratio = s(2,1)/s(1,2)
      s(2,2) = s(2,1) + s(1,1)
      s(2,3) = -s(1,1)
                                        go to 29
   28 ratio = s(2,1)/s(1,2)
      s(2,2) = s(2,1) + s(1,1)
      s(2,3) = -s(1,1)*6.0D+00*alpha*s(2,6)
c       at this point, the first two equations read
c              diag(1)*x1 + u(1)*x2 + entry3*x3 = r(2)
c       -ratio*diag(1)*x1 + diag(2)*x2 + u(2)*x3 = 0.0D+00
c       eliminate first unknown from second equation
   29 s(2,2) = ratio*s(1,3) + s(2,2)
      s(2,3) = ratio*entry3 + s(2,3)
      s(1,4) = s(2,4)
      s(2,4) = ratio*s(1,4)
                                        go to 35
   30 continue
c  ******the i-th equation enforces continuity of the first derivative
c        across tau(i). it has been set up in statements 35 up to 40
c        and 21 up to 25 and reads now
c         -ratio*diag(i-1)*xi-1 + diag(i)*xi + u(i)*xi+1 = r(i) .
c        eliminate (i-1)st unknown from this equation
      s(i,2) = ratio*s(i-1,3) + s(i,2)
      s(i,4) = ratio*s(i-1,4) + s(i,4)
c
c  ******set up the part of the next equation which depends on the
c        i-th interval.
   35 if (z - 0.5D+00)                       36,37,38
   36 ratio = -s(i,6)*s(i,1)/s(i,2)
      s(i+1,2) = s(i,1)/3.0D+00
                                        go to 40
   37 ratio = -(s(i,1)/6.0D+00)/s(i,2)
      s(i+1,2) = s(i,1)/3.0D+00
                                        go to 40
   38 ratio = -(s(i,1)/6.0D+00)/s(i,2)
      s(i+1,2) = s(i,1)*((1.0D+00-zeta*alpha)*factor/2. - s(i,6))
c         ------  end of i loop ------
   40 i = i+1
      if (i .lt. ntaum1)                go to 10
      s(i,5) = 0.5D+00
c
c        ------  last two equations  ------
c  the last two equations enforce continuity of third derivative and
c  of first derivative across  tau(ntau-1).
      entry = ratio*s(i-1,3) + s(i,2) + s(i,1)/3.0D+00
      s(i+1,2) = s(i,1)/6.0D+00
      s(i+1,4) = ratio*s(i-1,4) + s(i,4)
      if (z - 0.5D+00)                       41,42,43
   41 ratio = s(i,1)*6.0D+00*s(i-1,6)*alpha/s(i-1,2)
      s(i,2) = ratio*s(i-1,3) + s(i,1) + s(i-1,1)
      s(i,3) = -s(i-1,1)
                                        go to 45
   42 ratio = s(i,1)/s(i-1,2)
      s(i,2) = ratio*s(i-1,3) + s(i,1) + s(i-1,1)
      s(i,3) = -s(i-1,1)
                                        go to 45
   43 factr2 = onemzt*(alpha*(onemzt**2-1.0D+00)+1.0D+00)/
     *               (alpha*(onemzt**3-1.0D+00)+1.0D+00)
      ratio = factr2*s(i,1)/s(i-1,2)
      s(i,2) = ratio*s(i-1,3) + factr2*s(i-1,1) + s(i,1)
      s(i,3) = -factr2*s(i-1,1)
c     at this point, the last two equations read
c             diag(i)*xi + u(i)*xi+1 = r(i)
c      -ratio*diag(i)*xi + diag(i+1)*xi+1 = r(i+1)
c     eliminate xi from last equation
   45 s(i,4) = ratio*s(i-1,4)
      ratio = -entry/s(i,2)
      s(i+1,2) = ratio*s(i,3) + s(i+1,2)
      s(i+1,4) = ratio*s(i,4) + s(i+1,4)
c
c        ------ back substitution ------
c
      s(ntau,4) = s(ntau,4)/s(ntau,2)
   50    s(i,4) = (s(i,4) - s(i,3)*s(i+1,4))/s(i,2)
         i = i - 1
         if (i .gt. 1)                  go to 50
      s(1,4) = (s(1,4)-s(1,3)*s(2,4)-entry3*s(3,4))/s(1,2)
c
c        ------ construct polynomial pieces ------
c
      break(1) = tau(1)
      l = 1
      do 70 i=1,ntaum1
         coef(1,l) = gtau(i)
         coef(3,l) = s(i,4)
         divdif = (gtau(i+1)-gtau(i))/s(i,1)
         z = s(i,5)
         if (z - 0.5D+00)                    61,62,63
   61    if (z .eq. 0.0D+00)                 go to 65
         zeta = gam*z
         onemzt = 1.0D+00 - zeta
         c = s(i+1,4)/6.0D+00
         d = s(i,4)*s(i,6)
         l = l + 1
         del = zeta*s(i,1)
         break(l) = tau(i) + del
         zt2 = zeta**2
         alpha = alph(onemzt)
         factor = onemzt**2*alpha
         coef(1,l) = gtau(i) + divdif*del
     &             + s(i,1)**2*(d*onemzt*(factor-1.0D+00)
     &    +c*zeta*(zt2-1.0D+00))
         coef(2,l) = divdif + s(i,1)*(d*(1.0D+00-3.0D+00*factor)
     &  +c*(3.0D+00*zt2-1.0D+00))
         coef(3,l) = 6.0D+00*(d*alpha*onemzt + c*zeta)
         coef(4,l) = 6.0D+00*(c - d*alpha)/s(i,1)
         coef(4,l-1) = coef(4,l) - 6.0D+00*d*(1.0D+00-alpha)/(del*zt2)
         coef(2,l-1) = coef(2,l) - del*(coef(3,l) 
     &     -(del/2.0D+00)*coef(4,l-1))
                                        go to 68
   62    coef(2,l) = divdif - s(i,1)*(2.0D+00*s(i,4) + s(i+1,4))/6.0D+00
         coef(4,l) = (s(i+1,4)-s(i,4))/s(i,1)
                                        go to 68
   63    onemzt = gam*(1.0D+00 - z)
         if (onemzt .eq. 0.0D+00)            go to 65
         zeta = 1.0D+00 - onemzt
         alpha = alph(zeta)
         c = s(i+1,4)*s(i,6)
         d = s(i,4)/6.0D+00
         del = zeta*s(i,1)
         break(l+1) = tau(i) + del
         coef(2,l) = divdif - s(i,1)*(2.0D+00*d + c)
         coef(4,l) = 6.0D+00*(c*alpha - d)/s(i,1)
         l = l + 1
         coef(4,l) = coef(4,l-1) 
     &     + 6.0D+00*(1.0D+00-alpha)*c/(s(i,1)*onemzt**3)
         coef(3,l) = coef(3,l-1) + del*coef(4,l-1)
         coef(2,l) = coef(2,l-1)
     &     +del*(coef(3,l-1)+(del/2.0D+00)*coef(4,l-1))
         coef(1,l) = coef(1,l-1)
     &     +del*(coef(2,l-1)+(del/2.0D+00)*(coef(3,l-1)
     &                  +(del/3.0D+00)*coef(4,l-1)))
                                        go to 68
   65    coef(2,l) = divdif
         coef(3,l) = 0.0D+00
         coef(4,l) = 0.0D+00
   68    l = l + 1
   70    break(l) = tau(i+1)
      l = l - 1
      k = 4
      iflag = 1
                                        return
  999 iflag = 2
                                        return
      end
      subroutine titand ( tau, gtau, n )

c*********************************************************************72
c
cc TITAND represents a temperature dependent property of titanium.
c
c  from  * a practical guide to splines *  by c. de boor    
c  these data represent a property of titanium as a function of
c  temperature. they have been used extensively as an example in spline
c  approximation with variable knots.
c
      implicit none

      integer n,   i
      double precision gtau(49),tau(49),gtitan(49)
      data gtitan /
     &  0.644D+00, 0.622D+00, 0.638D+00, 0.649D+00, 0.652D+00,
     &  0.639D+00, 0.646D+00, 0.657D+00, 0.652D+00, 0.655D+00,
     &  0.644D+00, 0.663D+00, 0.663D+00, 0.668D+00, 0.676D+00,
     &  0.676D+00, 0.686D+00, 0.679D+00, 0.678D+00, 0.683D+00,
     &  0.694D+00, 0.699D+00, 0.710D+00, 0.730D+00, 0.763D+00,
     &  0.812D+00, 0.907D+00, 1.044D+00, 1.336D+00, 1.881D+00,
     &  2.169D+00, 2.075D+00, 1.598D+00, 1.211D+00, 0.916D+00,
     &  0.746D+00, 0.672D+00, 0.627D+00, 0.615D+00, 0.607D+00,
     &  0.606D+00, 0.609D+00, 0.603D+00, 0.601D+00, 0.603D+00,
     &  0.601D+00, 0.611D+00, 0.601D+00, 0.608D+00 /
      n = 49
      do 10 i=1,n
         tau(i) = 585.0D+00 + 10.0D+00*dble(i)
   10    gtau(i) = gtitan(i)
                                        return
      end
