      subroutine enum ( r, c, n, m, eval, ifault )

c*********************************************************************72
c
cc ENUM generates contingency tables with given shape and row and column sums.
c
c  Discussion:
c
c   The routine enumerates all M by N contingency tables with given row
c   and column totals, and calculates the hypergeometric probability of
c   each table.
c
c   For tables having two or more row sums repeated, equivalent
c   tables differing only by a row permutation are not separately
c   enumerated.  A representative of each equivalence class is enumerated
c   and the multiplicity of each class calculated.
c
c   For each table enumerated, subroutine EVAL is called to carry out
c   calculations on the table.
c
c   Note that the entries in the column sum and row sum vectors will
c   be (implicitly) sorted into ascending order, and results will be 
c   returned as though these orderings were being used!
c
c  Modified:
c
c    29 November 2006
c
c  Author:
c
c    Ian Saunders
c    Modifications by John Burkardt
c
c  Reference:
c
c    Ian Saunders,
c    Algorithm AS 205:
c    Enumeration of R x C Tables with Repeated Row Totals,
c    Applied Statistics,
c    Volume 33, Number 3, 1984, pages 340-352.
c
c  Parameters:
c
c    Input, integer R, the number of rows.
c
c    Input, integer C, the number of columns.
c
c    Input, integer N(R), the row sums.
c
c    Input, integer M(C), the column sums.
c
c    Input, external EVAL, the name of the user supplied routine
c    that will be called each time a new table is determined.
c    The routine has the form:
c
c      subroutine eval ( iflag, table, r, c, n, m, prob, mult )
c
c    See the dummy version of EVAL for an explanation of the meaning
c    of the arguments.
c
c    Output, integer IFAULT, error flag.
c    0, no error occurred.
c
c  Local parameters:
c
c    table(I,j)   - (I,J)-th entry of current table
c
c    NTOTAL - total of table entries
c
c    bound(I,J) - current upper bound on table(i,j) to satisfy
c    row and column totals
c
c    rept(I) - LOGICAL = true if row totals N(I), N(I-1) are equal
c    = false otherwise
c
c    reps(I) - number of previous rows equal to row I
c
c    MAX  - maximum dimension of TABLE
c
c    NMAX - maximum number of observations + 1
c
c    mult2(I) - maximum number of equivalent tables given first i rows
c
c    Z(J) - lower bound on sum of entries used by algorithm C
c
c    PROB2(I,J) - partial sum of terms in log(P)
c
      implicit none

      integer c
      integer mmax
      parameter ( mmax = 10 )
      integer nmax
      parameter ( nmax = 201 )
      integer r

      integer bound(mmax,mmax)
      integer c2
      integer cm
      external eval
      double precision factlm(nmax)
      double precision flogm(nmax)
      integer i
      logical ieqim
      integer ifault
      integer iflag
      integer ii
      integer iim
      integer iip
      integer ij
      integer j
      integer jj
      integer jjm
      integer jkeep
      integer jm
      integer jnext
      integer k
      integer keep
      integer left
      integer m(c)
      integer m2(mmax)
      integer maxrc
      integer mult
      integer mult2(mmax)
      integer multc
      integer multr
      integer mult0
      integer n(r)
      integer n2(mmax)
      integer ntot
      integer ntotal
      double precision prob
      double precision prob0
      double precision prob2(mmax,mmax)
      integer r2
      integer reps(mmax)
      integer repsc
      integer repsr
      logical rept(mmax)
      logical reptc(mmax)
      integer rm
      integer rowbnd
      integer rowsum
      integer table(mmax,mmax)
      integer table2(mmax,mmax)
      integer z(mmax)
c
c  Check the input values.
c
      ifault = 0

      if ( mmax .lt. r ) then
        ifault = 1
        return
      end if

      if ( mmax .lt. c ) then
        ifault = 1
        return
      end if

      if ( r .le. 0 ) then
        ifault = 1
        return
      end if

      if ( c .le. 0 ) then
        ifault = 1
        return
      end if

      ntotal = 0

      do i = 1, r
        if ( n(i) .le. 0 ) then
          ifault = 3
          return
        end if
        ntotal = ntotal + n(i)
      end do

      ntot = 0
      do j = 1, c
        if ( m(j) .le. 0 ) then
          ifault = 3
          return
        end if
        ntot = ntot + m(j)
      end do

      if ( ntot .ne. ntotal ) then
        ifault = 2
        return
      end if

      if ( nmax .le. ntotal ) then
        ifault = 4
        return
      end if
c
c  We need to copy the input parameters.  
c  One reason is that we may end up "transposing" the matrix,
c  so we need to use arrays that are large enough to work with
c  either the row or column dimension.
c
      c2 = c
      r2 = r

      do i = 1, r2
        n2(i) = n(i)
      end do

      do j = 1, c2
        m2(j) = m(j)
      end do
c
c  The first call to EVAL is simply to allow EVAL to initialize.
c
      iflag = 1
      prob = 0.0D+00
      mult = 0
      do j = 1, c
        do i = 1, r
          table(i,j) = 0
        end do
      end do

      call eval ( iflag, table, r, c, n, m, prob, mult ) 
c
c  Initialize FLOGM(K) = LOG(K-1), FACTLM(K) = LOG(K-1 FACTORIAL)
c
      flogm(1) = 0.0D+00
      factlm(1) = 0.0D+00
      do k = 1, ntotal
        flogm(k+1) = dlog ( dble ( k ) )
        factlm(k+1) = factlm(k) + flogm(k+1)
      end do
c
c  Constants.
c
      rm = r2 - 1
      cm = c2 - 1
c
c  Sort rows and columns into ascending order.
c
      do i = 1, r2 - 1
        do ii = i+1, r2
          if ( n2(ii) .lt. n2(i) ) then
            keep = n2(i)
            n2(i) = n2(ii)
            n2(ii) = keep
          end if
        end do
      end do

      do j = 1, c2 - 1
        do jj = j+1, c2
          if ( m2(jj) .lt. m2(j) ) then
            keep = m2(j)
            m2(j) = m2(jj)
            m2(jj) = keep
          end if
        end do
      end do
c
c  Calculate multiplicities of rows and columns.
c
c  REPTC(J) = TRUE if columns J and J-1 have the same total.
c  REPT(I)  = TRUE if rows I and I-1 have the same total.
c
      multc = 1
      repsc = 1
      reptc(1) = .false.

      do j = 2, c2

        reptc(j) = ( m2(j) .eq. m2(j-1) )

        if ( reptc(j) ) then
          repsc = repsc + 1
          multc = multc * repsc
        else
          repsc = 1
        end if

      end do

      multr = 1
      repsr = 1
      rept(1) = .false.

      do i = 2, r2

        rept(i) = ( n2(i) .eq. n2(i-1) )

        if ( rept(i) ) then
          repsr = repsr + 1
          multr = multr * repsr
        else
          repsr = 1
        end if

      end do
c
c  If column multiplicity exceeds row multiplicity, transpose the table.
c
      if ( multr .lt. multc ) then

        maxrc = max ( r2, c2 )

        do ij = 1, maxrc
          keep = n2(ij)
          n2(ij) = m2(ij)
          m2(ij) = keep
        end do

        keep = r2
        r2 = c2
        c2 = keep
        rm = r2 - 1
        cm = c2 - 1

        do i = 1, r2
          rept(i) = reptc(i)
        end do

        multr = multc

      end if
c
c  Set up the initial table.
c
c  Maximum multiplicity.
c
      mult2(1) = multr
      reps(1) = 1
c
c  Constant term in probability.
c
      prob0 = - factlm(ntotal+1)

      do i = 1, r2
        ii = n2(i)
        prob0 = prob0 + factlm(ii+1)
      end do

      do j = 1, c2
        jj = m2(j)
        prob0 = prob0 + factlm(jj+1)
      end do
c
c  Calculate bounds on row 1.
c
      do j = 1, c2
        bound(1,j) = m2(j)
      end do
c
c  For each I, find the greatest I-TH row satisfying the bounds.
c
      do i = 1, r2

        if ( i .ne. 1 ) then
          prob0 = prob2(i-1,c2)
        end if

        left = n2(i)
c
c  Elements of row I.
c
        ieqim = rept(i)

        do j = 1, cm

          ij = min ( left, bound(i,j) )
          table2(i,j) = ij

          if ( j .eq. 1 ) then
            prob2(i,j) = prob0 - factlm(ij+1)
          else
            prob2(i,j) = prob2(i,j-1) - factlm(ij+1)
          end if

          left = left - table2(i,j)

          if ( i .lt. r2 ) then
            bound(i+1,j) = bound(i,j) - table2(i,j)
          end if

          if ( left .eq. 0 ) then
            do jj = j+1, c2
              table2(i,jj) = 0
              prob2(i,jj) = prob2(i,jj-1)
              bound(i+1,jj) = bound(i,jj)
            end do
            go to 123
          end if

          if ( ieqim ) then
            ieqim = table2(i,j) .eq. table2(i-1,j)
          end if

        end do

        table2(i,c2) = left
        prob2(i,c2) = prob2(i,cm) - factlm(left+1)

        if ( i .lt. r2 ) then
          bound(i+1,c2) = bound(i,c2) - left
        end if

123     continue

        if ( i .ne. 1 ) then

          mult2(i) = mult2(i-1)
          reps(i) = 1

          if ( ieqim ) then
            reps(i) = reps(i-1) + 1
            mult2(i) = mult2(i) / reps(i)
          end if

        end if

      end do
c
c  Call EVAL for the first table.
c
      iflag = 2
      prob = prob2(r2,c2)
      mult = mult2(r2)

      if ( r .eq. r2 .and. c .eq. c2 ) then

        do j = 1, c
          do i = 1, r
            table(i,j) = table2(i,j)
          end do
        end do

      else

        do j = 1, c
          do i = 1, r
            table(i,j) = table2(j,i)
          end do
        end do

      end if

      call eval ( iflag, table, r, c, n, m, prob, mult )
c
c  Commence enumeration of remaining tables.
c
200   continue

      i = r2

210   continue

      i = i - 1
c
c  If I = 0, no more tables are possible.
c
      if ( i .eq. 0 ) then

        iflag = 3
        prob = 0.0D+00
        mult = 0
        do j = 1, c
          do i = 1, r
            table(i,j) = 0
          end do
        end do

        call eval ( iflag, table, r, c, n, m, prob, mult ) 

        return

      end if

      j = cm
      left = table2(i,c2)
      rowbnd = bound(i,c2)
c
c  Try to decrease element (I,J).
c
220   continue

      if ( 0 .lt. table2(i,j) .and. left .lt. rowbnd ) then
        go to 230
      end if
c
c  Element (I,J) cannot be decreased.  Try (I,J-1).
c
      if ( j .eq. 1 ) then
        go to 210
      end if

      left = left + table2(i,j)
      rowbnd = rowbnd + bound(i,j)
      j = j - 1
      go to 220
c
c  Decrease element (I,J).
c
230   continue

      ij = table2(i,j)
      prob2(i,j) = prob2(i,j) + flogm(ij+1)
      table2(i,j) = table2(i,j) - 1
      bound(I+1,j) = bound(i+1,j) + 1
c
c  If row I was the same as row I-1, it is no longer.
c
      if ( reps(i) .ne. 1 ) then
        reps(i) = 1
        mult2(i) = mult2(i-1)
      end if
c
c  Complete row I with the largest possible values.
c
      ii = i
      iip = ii + 1
      iim = ii - 1
      jnext = j + 1
      left = left + 1
      go to 380
c
c  Fill up remaining rows.
c
300   continue

      ii = ii + 1
c
c  The last row is treated separately.
c
      if ( ii .eq. r2 ) then
        go to 400
      end if

      iip = ii + 1
      iim = ii - 1
c
c  Row total N(II) is not a repeat.  Make row II as large as possible.
c
      if ( .not. rept(ii) ) then
        left = n2(ii)
        jnext = 1
        go to 380
      end if
c
c  Repeated row totals.
c
c  (I) If row II-1 satisfies the bounds on row II, repeat it.
c
      do j = 1, c2

        if ( bound(ii,j) .lt. table2(iim,j) ) then
          go to 330
        end if

        ij = table2(iim,j)
        table2(ii,j) = ij
        bound(iip,j) = bound(ii,j) - table2(ii,j)

        if ( j .eq. 1 ) then
          prob2(ii,j) = prob2(iim,c2) - factlm(ij+1)
        else
          prob2(ii,j) = prob2(ii,j-1) - factlm(ij+1)
        end if

      end do
c
c  Row II is a repeat of row II-1.
c
      reps(ii) = reps(iim) + 1
      mult2(ii) = mult2(iim) / reps(ii)
      go to 300
c
c  Element J of row II-1 was too big.
c
c  Construct the sequence Z(J) of lower bounds.
c
330   continue
c
c  If J = 1, the bounds are satisfied automatically.
c
      if ( j .eq. 1 ) then

        ij = bound(ii,1)
        table2(ii,1) = ij
        prob2(ii,1) = prob2(iim,c2) - factlm(ij+1)
        jnext = 2
        left = n2(ii) - table2(ii,1)
        bound(iip,1) = 0
        go to 380

      end if

      z(j) = n2(ii)
      jm = j - 1

      do jj = j+1, c2
        z(j) = z(j) - bound(ii,jj)
      end do

      do jjm = 1, jm
        jj = j - jjm
        z(jj) = z(jj+1) - bound(ii,jj+1)
      end do
c
c  (II) If the cumulative totals of row II-1 all exceed the bounds Z(J),
c  make element (II,J) equal to its bound.
c
      rowsum = 0
      jkeep = 0

      do jj = 1, jm

        rowsum = rowsum + table2(iim,jj)

        if ( rowsum .lt. z(jj) ) then
          go to 360
        end if

        if ( z(jj) .lt. rowsum .and. 0 .lt. table2(iim,jj) ) then
          jkeep = jj
        end if

      end do

      table2(ii,j) = bound(ii,j)
      bound(iip,j) = 0
      ij = table2(ii,j)
      prob2(ii,j) = prob2(ii,jm) - factlm(ij+1)
      reps(ii) = 1
      mult2(ii) = mult2(iim)
c
c  Complete row II with the largest possible elements.
c
      jnext = j + 1
      left = n2(ii)

      do jj = 1, j
        left = left - table2(ii,jj)
      end do

      go to 380
c
c  (III) The cumulative sums violate the bounds.
c  If no element of row II-1 can be changed to satisfy the bounds,
c  no suitable row II is possible.
c  In that case, go back and try decreasing row II-1.
c
360   continue

      if ( jkeep .eq. 0 ) then
        i = ii
        go to 210
      end if
c
c  Element (II,JKEEP) can be decreased.
c
370   continue

      bound(iip,jkeep) = bound(iip,jkeep) + 1
      ij = table2(ii,jkeep)
      prob2(ii,jkeep) = prob2(ii,jkeep) + flogm(ij+1)
      table2(ii,jkeep) = table2(ii,jkeep) - 1
c
c  Complete the row.
c
      jnext = jkeep + 1
      left = n2(ii)
      do jj = 1, jkeep
        left = left - table2(ii,jj)
      end do
c
c  Row II is complete up to element JNEXT-1.
c  Make the remaining elements as large as possible.
c  (This section of code is used for every row, repeated or not.)
c
380   continue

      do j = jnext, cm

        table2(ii,j) = min ( left, bound(ii,j) )
        left = left - table2(ii,j)
        bound(iip,j) = bound(ii,j) - table2(ii,j)
        ij = table2(ii,j)

        if ( j .eq. 1 ) then
          prob2(ii,j) = prob2(iim,c2) - factlm(ij+1)
        else
          prob2(ii,j) = prob2(ii,j-1) - factlm(ij+1)
        end if

        if ( left .eq. 0 ) then
          do jj = j+1, c2
            table2(ii,jj) = 0
            prob2(ii,jj) = prob2(ii,jj-1)
            bound(iip,jj) = bound(ii,jj)
          end do
          go to 393
        end if

      end do

      table2(ii,c2) = left
      prob2(ii,c2) = prob2(ii,cm) - factlm(left+1)
      bound(iip,c2) = bound(ii,c2) - left

393   continue

      reps(ii) = 1

      if ( 1 .lt. ii ) then
        mult2(ii) = mult2(iim)
      end if

      go to 300
c
c  The final row.
c
400   continue

      if ( .not. rept(r2) ) then
c
c  Not a repeat.  Set row R equal to its bounds.
c
        ij = bound(r2,1)
        table2(r2,1) = ij
        prob2(r2,1) = prob2(rm,c2) - factlm(ij+1)

        do j = 2, c2
          ij = bound(r2,j)
          table2(r2,j) = ij
          prob2(r2,j) = prob2(r2,j-1) - factlm(ij+1)
        end do

        mult2(r2) = mult2(rm)
        go to 500

      end if
c
c  Row total R is a repeat.  Ensure that it is less than row R-1.
c
      do j = 1, c2
c
c  If row R would be bigger than row R-1 go back and try
c  decreasing row R-2.
c
        if ( table2(rm,j) .lt. bound(r2,j) ) then
          i = rm
          go to 210
        end if

        ij = bound(r2,j)
        table2(r2,j) = ij

        if ( j .eq. 1 ) then
          prob2(r2,j) = prob2(rm,c2) - factlm(ij+1)
        else
          prob2(r2,j) = prob2(r2,j-1) - factlm(ij+1)
        end if
c
c  Row R is already less than row R-1, so no more checks are needed.
c
        if ( table2(r2,j) .ne. table2(rm,j) ) then

          do jj = j+1, c2
            ij = bound(r2,jj)
            table2(r2,jj) = ij
            prob2(r2,jj) = prob2(r2,jj-1) - factlm(ij+1)
          end do

          mult2(r2) = mult2(rm)
          go to 500
        end if

      end do
c
c  Row R is a repeat of row R-1.
c
      reps(r2) = reps(rm) + 1
      mult2(r2) = mult2(rm) / reps(r2)
c
c  The table is complete.
c
500   continue

      iflag = 2
      prob = prob2(r2,c2)
      mult = mult2(r2)

      if ( r .eq. r2 .and. c .eq. c2 ) then

        do j = 1, c
          do i = 1, r
            table(i,j) = table2(i,j)
          end do
        end do

      else

        do j = 1, c
          do i = 1, r
            table(i,j) = table2(j,i)
          end do
        end do

      end if

      call eval ( iflag, table, r, c, n, m, prob, mult )

      go to 200

      end
      subroutine eval ( iflag, table, m, n, rowsum, colsum, prob, mult )

c*********************************************************************72
c
cc EVAL is called by ENUM every time a new contingency table is determined.
c
c  Discussion:
c
c    This is a dummy version of the routine.
c
c    The user might wish to print out each contingency table, or collect
c    some statistics.
c
c  Modified:
c
c    27 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Ian Saunders,
c    Algorithm AS 205,
c    Enumeration of R x C Tables with Repeated Row Totals,
c    Applied Statistics,
c    Volume 33, Number 3, pages 340-352, 1984.
c
c  Parameters:
c
c    Input, integer IFLAG, input flag.
c    1, this is the first call.  No table is input.
c    2, this is a call with a new table.
c    3, this is the last call.  No table is input.
c
c    Input, integer TABLE(M,N), the current contingency table.
c
c    Input, integer M, the number of rows.
c
c    Input, integer N, the number of columns.
c
c    Input, integer ROWSUM(M), the row sums.
c
c    Input, integer COLSUM(N), the column sums.
c
c    Input, double precision PROB, the logarithm of the hypergeometric probability
c    of this table.
c
c    Input, integer MULT, the multiplicity of this table, that is,
c    the number of different tables that still have the same set of
c    entries, but differ by a permutation of some rows and columns.
c
      implicit none

      integer mmax
      parameter ( mmax = 10 )

      integer m
      integer n

      integer colsum(n)
      integer count1
      integer count2
      integer iflag
      integer mult
      double precision prob
      double precision psum
      integer rowsum(m)
      integer table(mmax,mmax)

      save count1
      save count2
      save psum

      data count1 / 0 /
      data count2 / 0 /
      data psum / 0.0D+00 /
c
c  First call, no table, initialize.
c
      if ( iflag .eq. 1 ) then

        count1 = 0
        count2 = 0
        psum = 0.0D+00
c
c  Call with a new table.
c
      else if ( iflag .eq. 2 ) then

        count1 = count1 + 1
        count2 = count2 + mult
        psum = psum + dble ( mult ) * dexp ( prob )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EVAL:'
        write ( *, '(i3,i3,g14.6)' ) count1, mult, prob

        call i4mat_print ( m, n, table, '  Table' )
c
c  Last call, no table.
c
      else if ( iflag .eq. 3 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EVAL summary'
        write ( *, '(a,i8)' ) 
     &  '  Number of cases (ignoring multiplicity):', count1
        write ( *, '(a,i8)' ) 
     &  '  Number of cases (allowing multiplicity):', count2
        write ( *, '(a,g14.6)' ) '  Probability sum = ', psum

      end if

      return
      end
      subroutine i4_swap ( i, j )

c*********************************************************************72
c
cc I4_SWAP switches two I4's.
c
c  Modified:
c
c    06 January 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer I, J.  On output, the values of I and
c    J have been interchanged.
c
      implicit none

      integer i
      integer j
      integer k

      k = i
      i = j
      j = k

      return
      end
      subroutine i4mat_print ( m, n, a, title )

c*********************************************************************72
c
cc I4MAT_PRINT prints an I4MAT.
c
c  Discussion:
c
c    An I4MAT is a rectangular array of integer values.
c
c  Modified:
c
c    30 June 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows in A.
c
c    Input, integer N, the number of columns in A.
c
c    Input, integer A(M,N), the matrix to be printed.
c
c    Input, character*(*) TITLE, a title to be printed first.
c    TITLE may be blank.
c
      implicit none

      integer mmax
      parameter ( mmax = 10 )

      integer m
      integer n

      integer a(mmax,mmax)
      integer ihi
      integer ilo
      integer jhi
      integer jlo
      character*(*) title

      ilo = 1
      ihi = m
      jlo = 1
      jhi = n

      call i4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

      return
      end
      subroutine i4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

c*********************************************************************72
c
cc I4MAT_PRINT_SOME prints some of an I4MAT.
c
c  Discussion:
c
c    An I4MAT is a rectangular array of integer values.
c
c  Modified:
c
c    04 November 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, integer A(M,N), an M by N matrix to be printed.
c
c    Input, integer ILO, JLO, the first row and column to print.
c
c    Input, integer IHI, JHI, the last row and column to print.
c
c    Input, character*(*) TITLE, an optional title.
c
      implicit none

      integer mmax
      parameter ( mmax = 10 )

      integer incx
      parameter ( incx = 10 )
      integer m
      integer n

      integer a(mmax,mmax)
      character*(8) ctemp(incx)
      integer i
      integer i2hi
      integer i2lo
      integer ihi
      integer ilo
      integer inc
      integer j
      integer j2
      integer j2hi
      integer j2lo
      integer jhi
      integer jlo
      integer s_len_trim
      character*(*) title

      if ( 0 .lt. s_len_trim ( title ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) title
      end if

      do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

        j2hi = j2lo + incx - 1
        j2hi = min ( j2hi, n )
        j2hi = min ( j2hi, jhi )

        inc = j2hi + 1 - j2lo

        write ( *, '(a)' ) ' '

        do j = j2lo, j2hi
          j2 = j + 1 - j2lo
          write ( ctemp(j2), '(i8)' ) j
        end do

        write ( *, '(''  Col '',10a8)' ) ( ctemp(j), j = 1, inc )
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) ' '

        i2lo = max ( ilo, 1 )
        i2hi = min ( ihi, m )

        do i = i2lo, i2hi

          do j2 = 1, inc

            j = j2lo - 1 + j2

            write ( ctemp(j2), '(i8)' ) a(i,j)

          end do

          write ( *, '(i5,1x,10a8)' ) i, ( ctemp(j), j = 1, inc )

        end do

      end do

      write ( *, '(a)' ) ' '

      return
      end
      subroutine i4vec_print ( n, a, title )

c*********************************************************************72
c
cc I4VEC_PRINT prints an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of integer values.
c
c  Modified:
c
c    27 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of components of the vector.
c
c    Input, integer A(N), the vector to be printed.
c
c    Input, character*(*) TITLE, a title to be printed first.
c    TITLE may be blank.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer s_len_trim
      character*(*) title
      integer title_length

      title_length = s_len_trim ( title )
      if ( 0 .lt. title_length ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) title(1:title_length)
      end if

      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,2x,i12)' ) i, a(i)
      end do

      return
      end
      function s_len_trim ( s )

c*********************************************************************72
c
cc S_LEN_TRIM returns the length of a string to the last nonblank.
c
c  Modified:
c
c    05 March 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*(*) S, a string.
c
c    Output, integer S_LEN_TRIM, the length of the string to the last nonblank.
c
      implicit none

      integer i
      character*(*) s
      integer s_len_trim

      do i = len ( s ), 1, -1

        if ( s(i:i) .ne. ' ' ) then
          s_len_trim = i
          return
        end if

      end do

      s_len_trim = 0

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

