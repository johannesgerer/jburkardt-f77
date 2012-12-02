      program main

c*********************************************************************72
c
cc MAIN is the main program for TOMS672_PRB.
c
c  Discussion:
c
c    The following templates demonstrate the use of the procedure EXTEND
c    to generate sequences of extended quadrature rules for various weight
c    functions given the definitions of the orthogonal polynomial 3-term
c    recurrence relations.
c
c    ICASE selects the required initial sequence as follows:
c    * 1, 3 point Gauss-Legendre in [-1,1]
c    * 2, 2 point Gauss-Lobatto in [-1,1]
c    * 3, 6 point Radau in [-1,1]
c    * 4, 2 point Gauss-Laguerre in [0,+oo)
c    * 5, 3 point Gauss-Hermite in (-oo,+oo)
c    * 6, 3 point Gauss-Jacobi in [0,1]
c
c  Modified:
c
c    04 March 2009
c
c  Author:
c
c    Original FORTRAN77 version by Thomas Patterson.
c
      implicit none

      integer icase
      integer nseq

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS672_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS672 library.'
c
c  Get data from terminal.
c
c  Select demonstration.
c
10    continue

      write ( *, '(a)' ) ' Enter ICASE, between 1 and 6, or -1 to stop:'
      read ( *, * ) icase

      if ( icase .lt. 1 .or. 6 .lt. icase ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TOMS672_PRB:'
        write ( *, '(a)' ) '  ICASE out of bounds.'
        write ( *, '(a)' ) '  Normal end of execution.'
        stop
      end if
c
c  Select number of iterative extensions to be performed.
c
      write ( *, '(a)' )
     &  '  Enter NSEQ, the number of nested rules to compute:'

      read ( *, * ) nseq

      if ( nseq .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TOMS672_PRB:'
        write ( *, '(a)' ) '  NSEQ out of bounds.'
        write ( *, '(a)' ) '  Normal end of execution.'
        stop
      end if

      if ( icase .eq. 1 ) then
        call test01 ( nseq )
      else if ( icase .eq. 2 ) then
        call test02 ( nseq )
      else if ( icase .eq. 3 ) then
        call test03 ( nseq )
      else if ( icase .eq. 4 ) then
        call test04 ( nseq )
      else if ( icase .eq. 5 ) then
        call test05 ( nseq )
      else if ( icase .eq. 6 ) then
        call test06 ( nseq )
      end if

      go to 10

      end
      subroutine test01 ( nseq )

c*********************************************************************72
c
cc TEST01: extension of 3 point Gauss-Legendre rule.
c
c  Discussion:
c
c    Generate a 3 point Gauss initially from 0 point rule,
c    No pre-assigned nodes.  Symmetry exploited.
c
c  Modified:
c
c    03 March 2009
c
c  Author:
c
c    Original FORTRAN77 version by Thomas Patterson.
c
c  Parameters:
c
c    Input, integer NSEQ, the number of nested sequences to compute.
c
      implicit none

      integer lda
      parameter ( lda = 1023 )
      integer ldb
      parameter ( ldb = 2 * lda + 1 )
      integer ntest
      parameter ( ntest = 4 )

      double precision err(lda)
      double precision ext(0:lda)
      double precision h0
      integer i
      integer ideg
      integer idigit
      parameter ( idigit = 8 )
      integer ierr
      integer iflag
      integer iwork(lda)
      integer k
      integer m
      integer m0
      integer n
      integer nexp
      parameter ( nexp = 38 )
      integer nodes
      integer nseq
      double precision pnodes(lda)
      double precision qi(lda)
      double precision qr(lda)
      external recura
      logical start
      logical symmet
      double precision t(0:lda)
      double precision test(0:ntest)
      double precision worka(lda,lda)
      double precision workb(ldb,3)
      double precision wt(lda)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Extension of a 3 point Gauss-Legendre rule.'

      n = 0
      m = 3
      m0 = 0
      t(0) = 1.0D+00
      symmet = .true.
      start = .false.

      do i = 1, nseq
c
c  Calculate the extension.
c
        h0 = 2.0D+00
        ideg = n + 2 * m - 1

        call extend ( n, m, m0, t, recura, symmet, start, pnodes, h0,
     &    nexp, idigit, wt, nodes, qr, qi, err, ext, iwork, worka, lda,
     &    workb, ldb, iflag )
c
c  Tests.
c
        if ( iflag .eq. 0 .or. iflag .eq. 6 ) then
          do k = 0, min ( ntest, ideg / 2 )
            call check ( n, pnodes, wt, k, h0, recura, test(k), ierr )
          end do
        end if
c
c  Print results.
c
        call results ( err, ext, i, ideg, iflag, iwork, lda, m, m0,
     &    n, nodes, ntest, pnodes, qi, qr, symmet, t, test, wt )
c
c  On next iteration, add N+1 nodes using the pre-assigned nodes defined by
c  the T polynomial generated in the previous cycle.
c
        m = n + 1
        start = .false.

      end do

      return
      end
      subroutine test02 ( nseq )

c*********************************************************************72
c
cc TEST02: extension of a 2 point Lobatto rule.
c
c  Discussion:
c
c    Add one node, pre-assign -1.0 and 1.0.  Symmetry exploited.
c
c  Modified:
c
c    03 March 2009
c
c  Author:
c
c    Original FORTRAN77 version by Thomas Patterson.
c
c  Parameters:
c
c    Input, integer NSEQ, the number of nested sequences to compute.
c
      implicit none

      integer lda
      parameter (lda = 257)
      integer ldb
      parameter ( ldb = 2 * lda + 1 )
      integer ntest
      parameter ( ntest = 4 )

      double precision err(lda)
      double precision ext(0:lda)
      double precision h0
      integer i
      integer ideg
      integer idigit
      parameter ( idigit = 8 )
      integer ierr
      integer iflag
      integer iwork(lda)
      integer k
      integer m
      integer m0
      integer n
      integer nexp
      parameter ( nexp = 38 )
      integer nodes
      integer nseq
      double precision pnodes(lda)
      double precision qi(lda)
      double precision qr(lda)
      external recura
      logical start
      logical symmet
      double precision t(0:lda)
      double precision test(0:ntest)
      double precision worka(lda,lda)
      double precision workb(ldb,3)
      double precision wt(lda)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  Extension of a Lobatto 2 point rule.'

      n = 2
      m = 1
      pnodes(1) = 1.0D+00
      pnodes(2) = - 1.0D+00
      symmet = .true.
      start = .true.

      do i = 1, nseq

        h0 = 2.0D+00
        ideg = n + 2 * m - 1

        call extend ( n, m, m0, t, recura, symmet, start, pnodes, h0,
     &    nexp, idigit, wt, nodes, qr, qi, err, ext, iwork, worka, lda,
     &    workb, ldb, iflag )
c
c  Tests.
c
        if ( iflag .eq. 0 .or. iflag .eq. 6 ) then
          do k  =  0, min ( ntest, ideg/2 )
            call check ( n, pnodes, wt, k, h0, recura, test(k), ierr )
          end do
        end if
c
c  Print results.
c
        call results ( err, ext, i, ideg, iflag, iwork, lda, m, m0,
     &    n, nodes, ntest, pnodes, qi, qr, symmet, t, test, wt )
c
c  On next iteration, add N-1 nodes using the pre-assigned nodes defined
c  by the T polynomial generated in the previous cycle.
c
        m = n - 1
        start = .false.

      end do

      return
      end
      subroutine test03 ( nseq )

c*********************************************************************72
c
cc TEST03: extension of 6-point Radau rule.
c
c  Discussion:
c
c    Add five nodes.  Pre-assign -1.0.  No symmetry.
c
c  Modified:
c
c    03 March 2009
c
c  Author:
c
c    Original FORTRAN77 version by Thomas Patterson.
c
c  Parameters:
c
c    Input, integer NSEQ, the number of nested sequences to compute.
c
      implicit none

      integer lda
      parameter (lda = 257)
      integer ldb
      parameter ( ldb = 2 * lda + 1 )
      integer ntest
      parameter ( ntest = 4 )

      double precision err(lda)
      double precision ext(0:lda)
      double precision h0
      integer i
      integer ideg
      integer idigit
      parameter ( idigit = 8 )
      integer ierr
      integer iflag
      integer iwork(lda)
      integer k
      integer m
      integer m0
      integer n
      integer nexp
      parameter ( nexp = 38 )
      integer nodes
      integer nseq
      double precision pnodes(lda)
      double precision qi(lda)
      double precision qr(lda)
      external recura
      logical start
      logical symmet
      double precision t(0:lda)
      double precision test(0:ntest)
      double precision worka(lda,lda)
      double precision workb(ldb,3)
      double precision wt(lda)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  Extension of a 6 point Radau rule.'

      n = 1
      m = 5
      pnodes(1) = - 1.0D+00
      symmet = .false.
      start = .true.

      do i = 1, nseq

        h0 = 2.0D+00
        ideg = n + 2 * m - 1

        call extend ( n, m, m0, t, recura, symmet, start, pnodes, h0,
     &    nexp, idigit, wt, nodes, qr, qi, err, ext, iwork, worka, lda,
     &    workb, ldb, iflag )
c
c  Tests.
c
        if ( iflag .eq. 0 .or. iflag .eq. 6 ) then
          do k = 0, min ( ntest, ideg / 2 )
            call check ( n, pnodes, wt, k, h0, recura, test(k), ierr )
          end do
        end if
c
c  Print results.
c
        call results ( err, ext, i, ideg, iflag, iwork, lda, m, m0,
     &    n, nodes, ntest, pnodes, qi, qr, symmet, t, test, wt )
c
c  Add N+1 nodes using the pre-assigned nodes defined by the T polynomial
c  generated in the previous cycle.
c
        m = n + 1
        start = .false.

      end do

      return
      end
      subroutine test04 ( nseq )

c*********************************************************************72
c
cc TEST04: extension of 2 point Gauss-Laguerre rule.
c
c  Discussion:
c
c    Generate a 2 point rule initially from a 0 point rule.
c    No pre-assigned nodes are used.
c
c    There is no symmetry.
c
c  Modified:
c
c    03 March 2009
c
c  Author:
c
c    Original FORTRAN77 version by Thomas Patterson.
c
c  Parameters:
c
c    Input, integer NSEQ, the number of nested sequences to compute.
c
      implicit none

      integer lda
      parameter (lda = 257)
      integer ldb
      parameter ( ldb = 2 * lda + 1 )
      integer ntest
      parameter ( ntest = 4 )

      double precision err(lda)
      double precision ext(0:lda)
      double precision h0
      integer i
      integer ideg
      integer idigit
      parameter ( idigit = 8 )
      integer ierr
      integer iflag
      integer iwork(lda)
      integer k
      integer m
      integer m0
      integer n
      integer nexp
      parameter ( nexp = 38 )
      integer nodes
      integer nseq
      double precision pnodes(lda)
      double precision qi(lda)
      double precision qr(lda)
      external recurb
      logical start
      logical symmet
      double precision t(0:lda)
      double precision test(0:ntest)
      double precision worka(lda,lda)
      double precision workb(ldb,3)
      double precision wt(lda)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) '  Extension of 2 point Gauss-Laguerre rule.'

      n = 0
      m = 2
      m0 = 0
      t(0) = 1.0D+00
      symmet = .false.
      start = .false.

      do i = 1, nseq
c
c  Calculate extension.
c
        h0 = 1.0D+00
        ideg = n + 2 * m - 1

        call extend ( n, m, m0, t, recurb, symmet, start, pnodes, h0,
     &    nexp, idigit, wt, nodes, qr, qi, err, ext, iwork, worka, lda,
     &    workb, ldb, iflag )
c
c  Tests.
c
        if ( iflag .eq. 0 .or. iflag .eq. 6 ) then
          do k = 0, min ( ntest, ideg / 2 )
            call check ( n, pnodes, wt, k, h0, recurb, test(k), ierr )
          end do
        end if
c
c  Print results.
c
        call results ( err, ext, i, ideg, iflag, iwork, lda, m, m0,
     &    n, nodes, ntest, pnodes, qi, qr, symmet, t, test, wt )
c
c  On the next iteration, add N+1 nodes using the pre-assigned nodes defined
c  by the T polynomial generated in the previous cycle.
c
        m = n + 1
        start = .false.

      end do

      return
      end
      subroutine test05 ( nseq )

c*********************************************************************72
c
cc TEST05: extension of 3-point Gauss-Hermite rule.
c
c  Discussion:
c
c    Generate a 3-point rule initially from a zero point rule.
c    That is, there are no pre-assigned nodes.
c
c    Symmetry is exploited.
c
c  Modified:
c
c    03 March 2009
c
c  Author:
c
c    Original FORTRAN77 version by Thomas Patterson.
c
c  Parameters:
c
c    Input, integer NSEQ, the number of nested sequences to compute.
c
      implicit none

      integer lda
      parameter (lda = 257)
      integer ldb
      parameter ( ldb = 2 * lda + 1 )
      integer ntest
      parameter ( ntest = 4 )

      double precision err(lda)
      double precision ext(0:lda)
      double precision h0
      integer i
      integer ideg
      integer idigit
      parameter ( idigit = 8 )
      integer ierr
      integer iflag
      integer iwork(lda)
      integer k
      integer m
      integer m0
      integer n
      integer nexp
      parameter ( nexp = 38 )
      integer nodes
      integer nseq
      double precision pnodes(lda)
      double precision qi(lda)
      double precision qr(lda)
      external recurc
      logical start
      logical symmet
      double precision t(0:lda)
      double precision test(0:ntest)
      double precision worka(lda,lda)
      double precision workb(ldb,3)
      double precision wt(lda)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) '  Extension of a 3 point Gauss-Hermite rule.'

      m = 3
      n = 0
      m0 = 0
      t(0) = 1.0D+00
      symmet = .true.
      start = .false.

      do i = 1, nseq
c
c  Calculate extension.
c  Zero moment integral = sqrt(pi)
c
        h0 = 2.0D+00 * sqrt ( atan ( 1.0D+00 ) )
        ideg = n + 2 * m - 1

        call extend ( n, m, m0, t, recurc, symmet, start, pnodes, h0,
     &    nexp, idigit, wt, nodes, qr, qi, err, ext, iwork, worka, lda,
     &    workb, ldb, iflag )
c
c  Tests.
c
        if ( iflag .eq. 0 .or. iflag .eq. 6 ) then
          do k = 0, min ( ntest, ideg / 2 )
            call check ( n, pnodes, wt, k, h0, recurc, test(k), ierr )
          end do
        end if
c
c  Print results.
c
        call results ( err, ext, i, ideg, iflag, iwork, lda, m, m0,
     &    n, nodes, ntest, pnodes, qi, qr, symmet, t, test, wt )
c
c  On the next iteration, add N+1 nodes using the pre-assigned nodes defined
c  by the T polynomial enerated in the previous cycle.
c
        m = n + 1
        start = .false.

      end do

      return
      end
      subroutine test06 ( nseq )

c*********************************************************************72
c
cc TEST06: extension of 3-point Gauss-Jacobi rule for weight sqrt(x) in [0,1].
c
c  Discussion:
c
c    Generate 3-point rule initially from zero point rule.
c
c    No pre-assigned nodes, no symmetry.
c
c  Modified:
c
c    03 March 2009
c
c  Author:
c
c    Original FORTRAN77 version by Thomas Patterson.
c
c  Parameters:
c
c    Input, integer NSEQ, the number of nested sequences to compute.
c
      implicit none

      integer lda
      parameter (lda = 257)
      integer ldb
      parameter ( ldb = 2 * lda + 1 )
      integer ntest
      parameter ( ntest = 4 )

      double precision err(lda)
      double precision ext(0:lda)
      double precision h0
      integer i
      integer ideg
      integer idigit
      parameter ( idigit = 8 )
      integer ierr
      integer iflag
      integer iwork(lda)
      integer k
      integer m
      integer m0
      integer n
      integer nexp
      parameter ( nexp = 38 )
      integer nodes
      integer nseq
      double precision pnodes(lda)
      double precision qi(lda)
      double precision qr(lda)
      external recurd
      logical start
      logical symmet
      double precision t(0:lda)
      double precision test(0:ntest)
      double precision worka(lda,lda)
      double precision workb(ldb,3)
      double precision wt(lda)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' ) '  Extension of a 3 point Gauss-Jacobi rule.'

      m = 3
      n = 0
      m0 = 0
      t(0) = 1.0D+00
      symmet = .false.
      start = .false.

      do i = 1, nseq
c
c  Calculate extension.
c
        h0 = 2.0D+00 / 3.0D+00
        ideg = n + 2 * m - 1

        call extend ( n, m, m0, t, recurd, symmet, start, pnodes, h0,
     &    nexp, idigit, wt, nodes, qr, qi, err, ext, iwork, worka, lda,
     &    workb, ldb, iflag )
c
c  Tests.
c
        if ( iflag .eq. 0 .or. iflag .eq. 6 ) then
          do k = 0, min ( ntest, ideg / 2 )
           call check ( n, pnodes, wt, k, h0, recurd, test(k), ierr )
          end do
        end if
c
c  Print results.
c
        call results ( err, ext, i, ideg, iflag, iwork, lda, m, m0,
     &    n, nodes, ntest, pnodes, qi, qr, symmet, t, test, wt )
c
c  On the next iteration, add N+1 nodes using the pre-assigned nodes defined
c  by the T polynomial generated in the previous cycle.
c
        m = n + 1
        start = .false.

      end do

      return
      end
      subroutine recura ( k, c, d, e )

c*********************************************************************72
c
cc RECURA is the recurrence used for Gauss, Lobatto and Radau.
c
c  Discussion:
c
c    This is an example of a user supplied subroutine to define the
c    orthogonal polynomials.
c
c    RECURA ( K, C, D, E ) gives the coefficients C, D and E such that:
c
c      P(K+1,X) = ( C * X + D ) * P(K,X) + E * P(K-1,X)
c
c  Modified:
c
c    02 March 2009
c
c  Author:
c
c    Original FORTRAN77 version by Thomas Patterson.
c
c  Parameters:
c
c    Input, integer K, the index.
c
c    Output, double precision C, D, E, the recurrence relation parameters.
c
      implicit none

      double precision c
      double precision d
      double precision e
      double precision f
      integer k

      f  =  dble ( k + 1 )
      c  =  dble ( 2 * k + 1 ) / f
      d  =  0.0D+00
      e  =  - dble ( k ) / f

      return
      end
      subroutine recurb ( k, c, d, e )

c*********************************************************************72
c
cc RECURB is the recurrence used for Laguerre rules.
c
c  Discussion:
c
c    This is an example of a user supplied subroutine to define the
c    orthogonal polynomials.
c
c    RECURB ( K, C, D, E ) gives the coefficients C, D and E such that:
c
c      P(K+1,X) = ( C * X + D ) * P(K,X) + E * P(K-1,X)
c
c  Modified:
c
c    02 March 2009
c
c  Author:
c
c    Original FORTRAN77 version by Thomas Patterson.
c
c  Parameters:
c
c    Input, integer K, the index.
c
c    Output, double precision C, D, E, the recurrence relation parameters.
c
      implicit none

      double precision c
      double precision d
      double precision e
      double precision f
      integer k

      f = dble ( k + 1 )
      c = - 1.0D+00 / f
      d = dble ( 2 * k + 1 ) / f
      e = - dble ( k ) / f

      return
      end
      subroutine recurc ( k, c, d, e )

c*********************************************************************72
c
cc RECURC is the recurrence used for Hermite.
c
c  Discussion:
c
c    This is an example of a user supplied subroutine to define the
c    orthogonal polynomials.
c
c    RECURC ( K, C, D, E ) gives the coefficients C, D and E such that:
c
c      P(K+1,X) = ( C * X + D ) * P(K,X) + E * P(K-1,X)
c
c  Modified:
c
c    02 March 2009
c
c  Author:
c
c    Original FORTRAN77 version by Thomas Patterson.
c
c  Parameters:
c
c    Input, integer K, the index.
c
c    Output, double precision C, D, E, the recurrence relation parameters.
c
      implicit none

      double precision c
      double precision d
      double precision e
      integer k

      c = 1.0D+00
      d = 0.0D+00
      e = - dble ( k ) / 2.0D+00

      return
      end
      subroutine recurd ( k, c, d, e )

c*********************************************************************72
c
cc RECURD is the recurrence used for Jacobi.
c
c  Discussion:
c
c    Jacobi polynomials in [0,1].
c
c    The weight function is (1-x)^(p-q) * x^(q-1).
c
c    This case for weight sqrt(x), p = 3/2 and q = p.
c
c    This is an example of a user supplied subroutine to define the
c    orthogonal polynomials.
c
c    RECURD ( K, C, D, E ) gives the coefficients C, D and E such that:
c
c      P(K+1,X) = ( C * X + D ) * P(K,X) + E * P(K-1,X)
c
c  Modified:
c
c    02 March 2009
c
c  Author:
c
c    Original FORTRAN77 version by Thomas Patterson.
c
c  Parameters:
c
c    Input, integer K, the index.
c
c    Output, double precision C, D, E, the recurrence relation parameters.
c
      implicit none

      double precision a1
      double precision a2
      double precision a3
      double precision a4
      double precision b
      double precision bp1
      double precision c
      double precision d
      double precision e
      double precision f2k
      double precision fk
      integer k
      double precision p
      double precision q
      double precision x3

      p = 1.5D+00
      q = p
      fk = dble ( k )
      f2k = 2.0D+00 * fk
      b = f2k + p
      bp1 = b + 1.0D+00
      x3 = ( b - 2.0D+00 ) * ( b - 1.0D+00 ) * b
      a1 = ( b - 1.0D+00 ) * bp1 * x3
      a2 = - x3 * ( f2k * ( fk + p ) + q * ( p - 1.0D+00 ) )
      a3 = x3 * bp1 * ( b - 1.0D+00 )
      a4 = fk * ( fk + q - 1.0D+00 ) * ( fk + p - 1.0D+00 )
     &  * ( fk + p - q ) * bp1
      c = a3 / a1
      d = a2 / a1
      e = - a4 / a1

      return
      end
      subroutine results ( err, ext, i, ideg, iflag, iwork, lda, m, m0,
     &  n, nodes, ntest, pnodes, qi, qr, symmet, t, test, wt )

c*********************************************************************72
c
cc RESULTS prints the results.
c
c  Modified:
c
c    04 March 2009
c
c  Author:
c
c    Original FORTRAN77 version by Thomas Patterson.
c
c  Parameters:
c
c    Input, double precision ERR(M), a measure of the relative error in the
c    nodes.  This may be inspected if the convergence error flag has been raised
c    (IFLAG = 3) to decide if the nodes in question are acceptable.  (ERR(*)
c    actually gives the mean last correction to the quadratic factor in the
c    generalized Bairstow root finder (see BAIR).
c
c    Input, double precision EXT(0:M), the coefficients of the polynomial whose
c    roots are the  extended nodes (QRNODES(*),QINODES(*)) and expressed as:
c      EXT = SUM (I = 0 to M) EXT(I)*P(I,X).
c
c    Input, integer I, the current stage of the set of nested quadrature rule.
c
c    Input, integer IDEG, the expected polynomial accuracy of the quadrature rule.
c
c    Input, integer IFLAG, error flag.
c    * 0, No error detected;
c    * 1, The linear system of equations defining the polynomial
c      whose roots are the extended nodes became singular or
c      very  ill-conditioned.   (FATAL).
c    * 2, The linear system of equations used to generate the
c      polynomial T when START is TRUE became singular
c      or very ill-conditioned. (FATAL).
c    * 3, Poor convergence has been detected in the calculation
c      of the roots of EXT corresponding to the new
c      nodes or all nodes have not been found (M not equal
c      to NODES). See also ERR(*).
c    * 4, Possible imaginary nodes detected.
c    * 5, Value of N and M incompatible for SYMMET = TRUE.
c      Both cannot be odd. (FATAL)
c    * 6, Test of new quadrature rule has failed.
c
c    Input, integer IWORK(max(M,N)), convergence flags.  Elements 1 to NODES
c    give information on the convergence of the roots of the polynomial EXT
c    corresponding to each extended node.
c    * 0, Convergence of I th root satisfactory;
c    * 1, Convergence of I th root unsatisfactory.
c
c    Input, integer LDA, the leading dimension of WORKA in the calling program.
c
c    Input, integer M, the number of nodes to be optimally added.
c
c    Input, integer M0, the lower limit to the expansion of T.
c
c    Input, integer N, the number of nodes in the current quadrature rule.
c
c    Input, integer NODES, the number of extended nodes found.  Normally this
c    equals M, but NODES will be less than M in cases where the computation
c    was terminated because an error condition was encountered.
c
c    Input, integer NTEST, the number of tests made of the quadrature rule.
c
c    Input, double precision PNODES(N), the nodes of the extended quadrature rule.
c
c    Input, double precision QINODE(M), the imaginary parts of the extended
c    nodes.
c
c    Input, double precision QRNODE(M), the real parts of the extended nodes.
c
c    Input, logical SYMMET,
c    * FALSE, if no advantage is to be taken of symmetry, if any,
c      about x = 0 in the interval of integration and the
c      orthogonality  weight function. Note that if symmetry in
c      fact does exist setting this parameter to zero will still
c      produce correct results - only efficiency is effected.
c    * TRUE, if the extended rule computations should
c      exploit symmetry about x = 0 in the interval of
c      integration and the orthogonality  weight function.
c      This reduces the size of the system of linear equations
c      determining EXT by a factor of about 2 (see WORKA). If
c      symmetry does not in fact exist erroneous results will be
c      produced.
c
c    Input, double precision T(0:M); the coefficients TI of the orthogonal
c    expansion whose roots are the nodes of the extended quadrature rule
c    (that is, the pre-assigned nodes plus the extended nodes) and expressed as:
c      SUM (I = M to N+M) (TI/HI)*P(I,X)
c    T(I-M) holds the value of TI.
c
c    Input, double precision TEST(NTEST), the results of tests on the
c    quadrature rule.
c
c    Input, double precision WT(N), the quadrature weights for
c    the extended rule associated with the nodes in PNODES.
c
      implicit none

      integer lda
      integer ntest

      double precision err(lda)
      double precision ext(0:lda)
      integer i
      integer ideg
      integer iflag
      integer iwork(lda)
      integer j
      integer k
      integer m
      integer m0
      integer n
      integer nodes
      integer num
      double precision pnodes(lda)
      double precision qi(lda)
      double precision qr(lda)
      logical symmet
      double precision t(0:lda)
      double precision test(0:ntest)
      double precision wt(lda)
c
c  Display results.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a,i3)' ) 'Iteration', i
      write ( *, '(a)' )
     &  'Coefficients of expansion whose roots are the new nodes:'
      write ( *, '(d25.16,a,i3,a)' )
     &  ( ext(j), '*P(', j, ',X)', j = 0, m )

      write ( *, '(a)' ) 'New nodes'
      write ( *, '(a,a)' )
     &  '                     Real                Imaginary',
     &  ' Flag       Err'
      write ( *, '(2d25.16,i5,d10.1)')
     &  ( qr(k), qi(k), iwork(k), err(k), k = 1, nodes )

      write ( *, '(a)' ) 'New full extended expansion'
      write ( *, '(d25.16,a,i3,a)')
     &  ( t(j-m0), '*P(',j, ',X)/HI', j = m0, n )

      write ( *, '(a,i2,a,i3,a,i1,a,i3)' )
     &  '  Complete extended rule: STEP = ', i,
     &  '  POINTS = ', n,
     &  '  IFLAG = ', iflag,
     &  '  Nodes added = ', nodes

      if ( iflag .ne. 0 .and. iflag .ne. 6 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Terminated prematurely.  See IFLAG.'
        return
      end if
c
c  Print rule (positive nodes only if symmetry present).
c
      if ( symmet ) then
        num = ( n + 1 ) / 2
      else
        num = n
      end if

      write ( *, '(a)' )
     &  '  No.                     Node                   Weight'
      write ( *, '(i5,2d25.16)' ) ( j, pnodes(j), wt(j), j = 1, num )
c
c  Display test results.
c
      do k = 0, min ( ntest, ideg / 2 )
        write ( *, '(a,i2,a,d25.16)' ) '  Test(', k, ') = ', test(k)
      end do

      if ( iflag .eq. 6 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Rule test is unsatisfactory.'
      end if

      return
      end
