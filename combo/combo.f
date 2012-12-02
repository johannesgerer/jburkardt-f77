      subroutine backtrack ( l, iarray, indx, k, nstack, stack, 
     &  maxstack )

c*********************************************************************72
c
cc BACKTRACK supervises a backtrack search.
c
c  Discussion:
c
c    The routine builds a vector, one element at a time, which is
c    required to satisfy some condition.
c
c    At any time, the partial vector may be discovered to be
c    unsatisfactory, but the routine records information about where the
c    last arbitrary choice was made, so that the search can be
c    carried out efficiently, rather than starting out all over again.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 April 1999
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input/output, integer L, the length of the completed
c    candidate vector.
c
c    Input/output, integer IARRAY(L), the candidate vector.
c
c    Input/output, integer INDX.
c    On input, set INDX = 0 to start a search.
c    On output:
c    1, a complete output vector has been determined.
c    2, candidates are needed.
c    3, no more possible vectors exist.
c
c    Input/output, integer K, the current length of the candidate
c    vector.
c
c    Input/output, integer NSTACK, the current length of the stack.
c
c    Input/output, integer STACK(MAXSTACK), a list of candidates
c    for positions 1 through K.
c
c    Input, integer MAXSTACK, the maximum length of the stack.
c
      implicit none

      integer l
      integer maxstack

      integer iarray(l)
      integer indx
      integer k
      integer nstack
      integer stack(maxstack)
c
c  If this is the first call, request a candidate for position 1.
c
      if ( indx .eq. 0 ) then
        k = 1
        nstack = 0
        indx = 2
        return
      end if
c
c  Examine the stack.
c
10    continue

        nstack = nstack - 1
c
c  If there are candidates for position K, take the first available
c  one off the stack, and increment K.
c
c  This may cause K to reach the desired value of L, in which case
c  we need to signal the user that a complete set of candidates
c  is being returned.
c
        if ( stack(nstack+1) .ne. 0 ) then

          iarray(k) = stack(nstack)
          stack(nstack) = stack(nstack+1) - 1

          if ( k .ne. l ) then
            k = k + 1
            indx = 2
          else
            indx = 1
          end if

          go to 20
c
c  If there are no candidates for position K, then decrement K.
c  If K is still positive, repeat the examination of the stack.
c
        else

          k = k - 1

          if ( k .le. 0 ) then
            indx = 3
            go to 20
          end if

        end if

      go to 10

20    continue

      return
      end
      subroutine bal_seq_check ( n, t, ierror )

c*********************************************************************72
c
cc BAL_SEQ_CHECK checks a balanced sequence.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of 0's (and 1's) in the sequence.
c    N must be positive.
c
c    Input, integer T(2*N), a balanced sequence.
c
c    Output, integer IERROR, error flag.
c    0, no error.
c    -1, N is not positive.
c    I, the I-th entry of T is illegal.
c    2*N+1, there are not the same number of 1's as 0's.
c
      implicit none

      integer n

      integer i
      integer ierror
      integer one_count
      integer t(2*n)
      integer zero_count

      ierror = 0

      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BAL_SEQ_CHECK - Fatal error!'
        write ( *, '(a)' ) '  N .lt. 1.'
        stop
      end if

      one_count = 0
      zero_count = 0

      do i = 1, 2 * n

        if ( t(i) .eq. 0 ) then
          zero_count = zero_count + 1
        else if ( t(i) .eq. 1 ) then
          one_count = one_count + 1
        else
          ierror = i
          return
        end if

        if ( zero_count .lt. one_count ) then
          ierror = 1
          return
        end if

      end do

      if ( one_count .ne. zero_count ) then
        ierror = 2 * n + 1
      end if

      return
      end
      subroutine bal_seq_enum ( n, nseq )

c*********************************************************************72
c
cc BAL_SEQ_ENUM enumerates the balanced sequences.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of 0's (and 1's) in the sequence.
c    N must be nonnegative.
c
c    Output, integer NSEQ, the number of balanced sequences.
c
      implicit none

      integer binomial
      integer n
      integer nseq

      nseq = binomial ( 2 * n, n ) / ( n + 1 )

      return
      end
      subroutine bal_seq_rank ( n, t, rank )

c*********************************************************************72
c
cc BAL_SEQ_RANK ranks a balanced sequence.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of 0's (and 1's) in the sequence.
c    N must be positive.
c
c    Input, integer T(2*N), a balanced sequence.
c
c    Output, integer RANK, the rank of the balanced sequence.
c
      implicit none

      integer n

      integer ierror
      integer mxy
      integer rank
      integer t(2*n)
      integer x
      integer y
c
c  Check.
c
      call bal_seq_check ( n, t, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BAL_SEQ_RANK - Fatal error!'
        write ( *, '(a)' ) '  The input array is illegal.'
        write ( *, '(a,i8)' ) '  IERROR = ', ierror
        stop
      end if

      y = 0
      rank = 0

      do x = 1, 2 * n - 1

        if ( t(x) .eq. 0 ) then
          y = y + 1
        else
          call mountain ( n, x, y + 1, mxy )
          rank = rank + mxy
          y = y - 1
        end if

      end do

      return
      end
      subroutine bal_seq_successor ( n, t, rank )

c*********************************************************************72
c
cc BAL_SEQ_SUCCESSOR computes the lexical balanced sequence successor.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of 0's (and 1's) in the sequence.
c    N must be positive.
c
c    Input/output, integer T(2*N), on input, a balanced sequence,
c    and on output, its lexical successor.
c
c    Input/output, integer RANK, the rank.
c    If RANK = -1 on input, then the routine understands that this is
c    the first call, and that the user wishes the routine to supply
c    the first element in the ordering, which has RANK = 0.
c    In general, the input value of RANK is increased by 1 for output,
c    unless the very last element of the ordering was input, in which
c    case the output value of RANK is 0.
c
      implicit none

      integer n

      integer i
      integer ierror
      integer j
      integer open
      integer open_index
      integer rank
      integer slot
      integer slot_index
      integer slot_ones
      integer t(2*n)
c
c  Return the first element.
c
      if ( rank .eq. -1 ) then
        do i = 1, n
          t(i) = 0
        end do
        do i = n + 1, 2 * n
          t(i) = 1
        end do
        rank = 0
        return
      end if
c
c  Check.
c
      call bal_seq_check ( n, t, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BAL_SEQ_SUCCESSOR - Fatal error!'
        write ( *, '(a)' ) '  The input array is illegal.'
        write ( *, '(a,i8)' ) '  IERROR = ', ierror
        stop
      end if
c
c  After the I-th 0 there is a 'slot' with the capacity to
c  hold between 0 and I ones.
c
c  The first element of the sequence has all the 1's cowering
c  behind the N-th 0.
c
c  We seek to move a 1 to the left, and to do it lexically,
c  we will move a 1 to the rightmost slot that is under capacity.
c
c  Find the slot.
c
      slot = 0
      slot_index = 0
      slot_ones = 0

      open = 0
      open_index = 0

      do i = 1, 2 * n

        if ( t(i) .eq. 0 ) then

          if ( 0 .lt. slot ) then
            if ( slot_ones .lt. slot ) then
              open = slot
              open_index = slot_index
            end if
          end if

          slot = slot + 1
          slot_index = i

c     slot_ones = 0

        else
          slot_ones = slot_ones + 1
        end if

      end do
c
c  If OPEN is not 0, then preserve the string up to the OPEN-th 0,
c  preserve the 1's that follow, but then write a 1, then
c  all the remaining 0's and all the remaining 1's.
c
      if ( open .ne. 0 ) then

        j = open_index + 1

10      continue

        if ( t(j) .eq. 1 ) then
          j = j + 1
          go to 10
        end if

        t(j) = 1

        do i = open + 1, n
          j = j + 1
          t(j) = 0
        end do

        do i = j + 1, 2 * n
          t(i) = 1
        end do
c
c  If OPEN is 0, the last element was input.
c  Return the first one.
c
      else

        do i = 1, n
          t(i) = 0
        end do
        do i = n + 1, 2 * n
          t(i) = 1
        end do
        rank = 0
        return

      end if

      rank = rank + 1

      return
      end
      subroutine bal_seq_to_tableau ( n, t, tab )

c*********************************************************************72
c
cc BAL_SEQ_TO_TABLEAU converts a balanced sequence to a 2 by N tableau.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of 0's (and 1's) in the sequence.
c    N must be positive.
c
c    Input, integer T(2*N), a balanced sequence.
c
c    Output, integer TAB(2,N), a 2 by N tableau.
c
      implicit none

      integer n

      integer c(2)
      integer i
      integer ierror
      integer r
      integer t(2*n)
      integer tab(2,n)
c
c  Check.
c
      call bal_seq_check ( n, t, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BAL_SEQ_TO_TABLEAU - Fatal error!'
        write ( *, '(a)' ) '  The input array is illegal.'
        write ( *, '(a,i8)' ) '  IERROR = ', ierror
        stop
      end if

      c(1) = 0
      c(2) = 0
      do i = 1, 2 * n
        r = t(i) + 1
        c(r) = c(r) + 1
        tab(r,c(r)) = i
      end do

      return
      end
      subroutine bal_seq_unrank ( rank, n, t )

c*********************************************************************72
c
cc BAL_SEQ_UNRANK unranks a balanced sequence.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer RANK, the rank of the balanced sequence.
c
c    Input, integer N, the number of 0's (and 1's) in the sequence.
c    N must be positive.
c
c    Output, integer T(2*N), a balanced sequence.
c
      implicit none

      integer n

      integer low
      integer m
      integer nseq
      integer rank
      integer t(2*n)
      integer x
      integer y
c
c  Check.
c
      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BAL_SEQ_UNRANK - Fatal error!'
        write ( *, '(a)' ) '  Input N is illegal.'
        stop
      end if

      call bal_seq_enum ( n, nseq )

      if ( rank .lt. 0 .or. nseq .lt. rank ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BAL_SEQ_UNRANK - Fatal error!'
        write ( *, '(a)' ) '  The input rank is illegal.'
        stop
      end if

      y = 0
      low = 0

      do x = 1, 2 * n

        call mountain ( n, x, y + 1, m )

        if ( rank .le. low + m - 1 ) then
          y = y + 1
          t(x) = 0
        else
          low = low + m
          y = y - 1
          t(x) = 1
        end if

      end do

      return
      end
      subroutine bell_numbers ( m, b )

c*********************************************************************72
c
cc BELL_NUMBERS computes the Bell numbers.
c
c  Discussion:
c
c    There are B(M) restricted growth functions of length M.
c
c    There are B(M) partitions of a set of M objects.
c
c    B(M) is the sum of the Stirling numbers of the second kind,
c    S(M,N), for N = 0 to M.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer M, indicates how many Bell numbers are to
c    compute.  M must be nonnegative.
c
c    Output, integer B(0:M), the first M+1 Bell numbers.
c
      implicit none

      integer m

      integer b(0:m)
      integer binomial
      integer i
      integer j

      b(0) = 1
      do j = 1, m
        b(j) = 0
        do i = 0, j - 1
          b(j) = b(j) + binomial ( j - 1, i ) * b(i)
        end do
      end do

      return
      end
      subroutine bell_values ( n_data, n, c )

c*********************************************************************72
c
cc BELL_VALUES returns some values of the Bell numbers for testing.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Milton Abramowitz, Irene Stegun,
c    Handbook of Mathematical Functions,
c    National Bureau of Standards, 1964,
c    ISBN: 0-486-61272-4,
c    LC: QA47.A34.
c
c  Parameters:
c
c    Input/output, integer N_DATA.
c    On input, if N_DATA is 0, the first test data is returned, and N_DATA
c    is set to 1.  On each subsequent call, the input value of N_DATA is
c    incremented and that test data item is returned, if available.  When
c    there is no more test data, N_DATA is set to 0.
c
c    Output, integer N, the order of the Bell number.
c
c    Output, integer C, the value of the Bell number.
c
      implicit none

      integer n_max
      parameter ( n_max = 11 )

      integer c
      integer c_vec(n_max)
      integer n
      integer n_data
      integer n_vec(n_max)

      save c_vec
      save n_vec

      data c_vec /
     &  1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147, 115975 /
      data n_vec /
     &   0,  1,  2,  3,  4, 5,  6,  7,  8,  9,  10 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        n = 0
        c = 0
      else
        n = n_vec(n_data)
        c = c_vec(n_data)
      end if

      return
      end
      function binomial ( n, k )

c*********************************************************************72
c
cc BINOMIAL computes the binomial coefficient C(N,K).
c
c  Discussion:
c
c    BINOMIAL(N,K) = C(N,K) = N! / ( K! * (N-K)! )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    ML Wolfson, HV Wright,
c    Algorithm 160:
c    Combinatorial of M Things Taken N at a Time,
c    Communications of the ACM,
c    Volume 6, Number 4, April 1963, page 161.
c
c  Parameters:
c
c    Input, integer N, K, are the values of N and K.
c
c    Output, integer BINOMIAL, the number of combinations of N
c    things taken K at a time.
c
      implicit none

      integer binomial
      integer i
      integer icnk
      integer k
      integer mn
      integer mx
      integer n

      mn = min ( k, n - k )

      if ( mn .lt. 0 ) then

        icnk = 0

      else if ( mn .eq. 0 ) then

        icnk = 1

      else

        mx = max ( k, n - k )
        icnk = mx + 1

        do i = 2, mn
          icnk = ( icnk * ( mx + i ) ) / i
        end do

      end if

      binomial = icnk

      return
      end
      subroutine combin ( n, k, cnk )

c*********************************************************************72
c
cc COMBIN computes the combinatorial coefficient C(N,K).
c
c  Discussion:
c
c    Real arithmetic is used, and C(N,K) is computed directly, via
c    Gamma functions, rather than recursively.
c
c    C(N,K) is the number of distinct combinations of K objects
c    chosen from a set of N distinct objects.  A combination is
c    like a set, in that order does not matter.
c
c    C(N,K) = N! / ( (N-K)! * K! )
c
c  Example:
c
c    The number of combinations of 2 things chosen from 5 is 10.
c
c    C(5,2) = ( 5 * 4 * 3 * 2 * 1 ) / ( ( 3 * 2 * 1 ) * ( 2 * 1 ) ) = 10.
c
c    The actual combinations may be represented as:
c
c      (1,2), (1,3), (1,4), (1,5), (2,3),
c      (2,4), (2,5), (3,4), (3,5), (4,5).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 June 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the value of N.
c
c    Input, integer K, the value of K.
c
c    Output, double precision CNK, the value of C(N,K)
c
      implicit none

      double precision arg
      double precision cnk
      double precision fack
      double precision facn
      double precision facnmk
      integer k
      integer n
      double precision r8_gamma_log

      if ( n .lt. 0 ) then

        cnk = 0.0D+00

      else if ( k .eq. 0 ) then

        cnk = 1.0D+00

      else if ( k .eq. 1 ) then

        cnk = dble ( n )

      else if ( 1 .lt. k .and. k .lt. n - 1 ) then

        arg = dble ( n + 1 )
        facn = r8_gamma_log ( arg )

        arg = dble ( k + 1 )
        fack = r8_gamma_log ( arg )

        arg = dble ( n - k + 1 )
        facnmk = r8_gamma_log ( arg )

        cnk = dble ( nint ( exp ( facn - fack - facnmk ) ) )

      else if ( k .eq. n - 1 ) then

        cnk = dble ( n )

      else if ( k .eq. n ) then

        cnk = 1.0D+00

      else

        cnk = 0.0D+00

      end if

      return
      end
      subroutine cycle_check ( n, ncycle, t, index, ierror )

c*********************************************************************72
c
cc CYCLE_CHECK checks a permutation in cycle form.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 July 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of items permuted.
c    N must be positive.
c
c    Input, integer NCYCLE, the number of cycles.
c    1 .le. NCYCLE .le. N.
c
c    Input, integer T(N), INDEX(NCYCLE), describes the permutation
c    as a collection of NCYCLE cycles.  The first cycle is
c    T(1) -> T(2) -> ... -> T(INDEX(1)) -> T(1).
c
c    Output, integer IERROR, error flag.
c    0, no error.
c    -1, N is less than 1.
c    -2, NCYCLE is less than 1 or greater than N.
c    -3, an entry of INDEX is illegal.
c    -4, the entries of INDEX do not add up to N.
c    I, entry I of T is illegal.
c
      implicit none

      integer n
      integer ncycle

      integer i
      integer i4vec_sum
      integer ierror
      integer ifind
      integer index(ncycle)
      integer iseek
      integer t(n)

      ierror = 0
c
c  N must be at least 1.
c
      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CYCLE_CHECK - Fatal error!'
        write ( *, '(a)' ) '  N .lt. 1.'
        stop
      end if
c
c  1 .le. NCYCLE .le. N.
c
      if ( ncycle .lt. 1 .or. n .lt. ncycle ) then
        ierror = -2
        return
      end if
c
c  1 .le. INDEX(I) .le. N.
c
      do i = 1, ncycle
        if ( index(i) .lt. 1 .or. n .lt. index(i) ) then
          ierror = -3
          return
        end if
      end do
c
c  The INDEX(I)'s sum to N.
c
      if ( i4vec_sum ( ncycle, index ) .ne. n ) then
        ierror = -4
        return
      end if
c
c  1 .le. T(I) .le. N.
c
      do i = 1, n
        if ( t(i) .lt. 1 .or. n .lt. t(i) ) then
          ierror = i
          return
        end if
      end do
c
c  Verify that every value from 1 to N occurs in T.
c
      do iseek = 1, n

        ifind = 0

        do i = 1, n
          if ( t(i) .eq. iseek ) then
            ifind = i
            go to 10
          end if
        end do

10      continue

        if ( ifind .eq. 0 ) then
          ierror = iseek
          return
        end if

      end do

      return
      end
      subroutine cycle_to_perm ( n, ncycle, t, index, p )

c*********************************************************************72
c
cc CYCLE_TO_PERM converts a permutation from cycle to array form.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of items permuted.
c    N must be positive.
c
c    Input, integer NCYCLE, the number of cycles.
c    1 .le. NCYCLE .le. N.
c
c    Input, integer T(N), INDEX(NCYCLE), describes the permutation
c    as a collection of NCYCLE cycles.  The first cycle is
c    T(1) -> T(2) -> ... -> T(INDEX(1)) -> T(1).
c
c    Output, integer P(N), describes the permutation using a
c    single array.  For each index I, I -> P(I).
c
      implicit none

      integer n
      integer ncycle

      integer i
      integer ierror
      integer index(ncycle)
      integer j
      integer jhi
      integer jlo
      integer p(n)
      integer t(n)
c
c  Check.
c
      call cycle_check ( n, ncycle, t, index, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CYCLE_TO_PERM - Fatal error!'
        write ( *, '(a)' ) '  The input array is illegal.'
        write ( *, '(a,i8)' ) '  IERROR = ', ierror
        stop
      end if

      jhi = 0

      do i = 1, ncycle

        jlo = jhi + 1
        jhi = jhi + index(i)

        do j = jlo, jhi

          if ( j .lt. jhi ) then
            p(t(j)) = t(j+1)
          else
            p(t(j)) = t(jlo)
          end if

        end do

      end do

      return
      end
      subroutine dist_enum ( k, m, num_dist )

c*********************************************************************72
c
cc DIST_ENUM returns the number of distributions of indistinguishable objects.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 June 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer K, the number of distinguishable "slots".
c
c    Input, integer M, the number of indistinguishable objects.
c
c    Output, integer NUM_DIST, the number of distributions of M
c    indistinguishable objects about K distinguishable slots.
c
      implicit none

      double precision cnk
      integer k
      integer m
      integer num_dist

      call combin ( m + k - 1, m, cnk )

      num_dist = nint ( cnk )

      return
      end
      subroutine dist_next ( k, m, q, more )

c*********************************************************************72
c
cc DIST_NEXT returns the next distribution of indistinguishable objects.
c
c  Discussion:
c
c    A distribution of M objects into K parts is an ordered sequence
c    of K nonnegative integers which sum to M.  This is similar to
c    a partition of a set into K subsets, except that here the order
c    matters.  That is, (1,1,2) and (1,2,1) are considered to be
c    different distributions.
c
c    On the first call to this routine, the user should set MORE = FALSE,
c    to signal that this is a startup for the given computation.  The routine
c    will return the first distribution, and set MORE = TRUE.
c
c    If the user calls again, with MORE = TRUE, the next distribution
c    is being requested.  If the routine returns with MORE = TRUE, then
c    that distribution was found and returned.  However, if the routine
c    returns with MORE = FALSE, then no more distributions were found;
c    the enumeration of distributions has terminated.
c
c    A "distribution of M indistinguishable objects into K slots" is
c    sometimes called a "composition of the integer M into K parts".
c
c  Example:
c
c    K = 3, M = 5
c
c    0           0           5
c    0           1           4
c    0           2           3
c    0           3           2
c    0           4           1
c    0           5           0
c    1           0           4
c    1           1           3
c    1           2           2
c    1           3           1
c    1           4           0
c    2           0           3
c    2           1           2
c    2           2           1
c    2           3           0
c    3           0           2
c    3           1           1
c    3           2           0
c    4           0           1
c    4           1           0
c    5           0           0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 November 2005
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Fenichel,
c    Algorithm 329:
c    Distribution of Indistinguishable Objects into
c    Distinguishable Slots,
c    Communications of the ACM,
c    Volume 11, Number 6, June 1968, page 430.
c
c  Parameters:
c
c    Input, integer K, the number of distinguishable "slots".
c
c    Input, integer M, the number of indistinguishable objects.
c
c    Input/output, integer Q(K), the number of objects in each
c    slot.
c
c    Input/output, logical MORE, used by the user to start the computation,
c    and by the routine to stop the computation.
c
      implicit none

      integer k

      integer i
      integer leftmost
      integer m
      logical more
      integer q(k)

      save leftmost

      data leftmost / 1 /
c
c  The startup call.
c
      if ( .not. more ) then

        more = .true.

        do i = 1, k - 1
          q(i) = 0
        end do
        q(k) = m

        leftmost = k + 1
c
c  There are no more distributions.
c  Reset Q to the first distribution in the sequence.
c
      else if ( q(1) .eq. m ) then

        more = .false.

        do i = 1, k - 1
          q(i) = 0
        end do
        q(k) = m

        leftmost = k + 1

      else if ( leftmost .lt. k + 1 ) then

        leftmost = leftmost - 1
        q(k) = q(leftmost) - 1
        q(leftmost) = 0
        q(leftmost-1) = q(leftmost-1) + 1
        if ( q(k) .ne. 0 ) then
          leftmost = k + 1
        end if

      else

        if ( q(k) .eq. 1 ) then
          leftmost = k
        end if

        q(k) = q(k) - 1
        q(k-1) = q(k-1) + 1

      end if

      return
      end
      subroutine edge_check ( n_node, n_edge, t, ierror )

c*********************************************************************72
c
cc EDGE_CHECK checks a graph stored by edges.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N_NODE, the number of nodes in the graph.
c    N_NODE must be positive.
c
c    Input, integer N_EDGE, the number of edges in the graph.
c    N_EDGE must be positive.
c
c    Input, integer T(2,N_EDGE), describes the edges of the tree
c    as pairs of nodes.
c
c    Output, integer IERROR, error flag.
c    -1, N_NODE is not positive.
c    -2, N_EDGE is not positive.
c    0, no error.
c    I, edge T(1,I), T(2,I) is illegal.
c
      implicit none

      integer n_edge
      integer n_node

      integer i
      integer ierror
      integer j
      integer j2
      integer t(2,n_edge)

      ierror = 0

      if ( n_node .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EDGE_CHECK - Fatal error!'
        write ( *, '(a)' ) '  N_NODE .lt. 1.'
        stop
      end if

      if ( n_edge .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EDGE_CHECK - Fatal error!'
        write ( *, '(a)' ) '  N_EDGE .lt. 1.'
        stop
      end if
c
c  Every edge must join two legal nodes.
c
      do i = 1, 2
        do j = 1, n_edge
          if ( t(i,j) .lt. 1 .or. n_node .lt. t(i,j) ) then
            ierror = i
            return
          end if
        end do
      end do
c
c  Every edge must join distinct nodes.
c
      do j = 1, n_edge
        if ( t(1,j) .eq. t(2,j) ) then
          ierror = i
          return
        end if
      end do
c
c  Every edge must be distinct.
c
      do j = 1, n_edge - 1
        do j2 = j + 1, n_edge
          if ( t(1,j) .eq. t(1,j2) .and. t(2,j) .eq. t(2,j2) ) then
            ierror = j2
            return
          else if ( t(1,j) .eq. t(2,j2) .and. 
     &              t(2,j) .eq. t(1,j2) ) then
            ierror = j2
            return
          end if
        end do
      end do

      return
      end
      subroutine edge_degree ( n_node, n_edge, t, d )

c*********************************************************************72
c
cc EDGE_DEGREE returns the degree of the nodes of a graph stored by edges.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N_NODE, the number of nodes in the graph.
c    N_NODE must be positive.
c
c    Input, integer N_EDGE, the number of edges in the graph.
c    N_EDGE must be positive.
c
c    Input, integer T(2,N_EDGE), describes the edges of the tree
c    as pairs of nodes.
c
c    Output, integer D(N_NODE), the degree of each node.
c
      implicit none

      integer n_edge
      integer n_node

      integer d(n_node)
      integer i
      integer ierror
      integer j
      integer t(2,n_edge)
c
c  Check.
c
      call edge_check ( n_node, n_edge, t, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EDGE_DEGREE - Fatal error!'
        write ( *, '(a)' ) '  The input array is illegal.'
        write ( *, '(a,i8)' ) '  IERROR = ', ierror
        stop
      end if
c
c  Compute the degree of each node.
c
      do i = 1, n_node
        d(i) = 0
      end do

      do j = 1, n_edge
        d(t(1,j)) = d(t(1,j)) + 1
        d(t(2,j)) = d(t(2,j)) + 1
      end do

      return
      end
      subroutine edge_enum ( n_node, nedge )

c*********************************************************************72
c
cc EDGE_ENUM enumerates the maximum number of edges in a graph on N_NODE nodes.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N_NODE, the number of nodes in the graph.
c    N_NODE must be positive.
c
c    Output, integer NEDGE, the maximum number of edges in a graph
c    on N_NODE nodes.
c
      implicit none

      integer n_node
      integer nedge

      nedge = ( n_node * ( n_node - 1 ) ) / 2

      return
      end
      function fall ( x, n )

c*********************************************************************72
c
cc FALL computes the falling factorial function [X]_N.
c
c  Discussion:
c
c    The number of "injections" or 1-to-1 mappings from
c    a set of N elements to a set of M elements is [M]_N.
c
c    The number of permutations of N objects out of M is [M}_N.
c
c    The Stirling numbers of the first kind can be used
c    to convert a falling factorial into a polynomial, as follows:
c
c      [X]_N = S^0_N + S^1_N * X + S^2_N * X^2 + ... + S^N_N X^N.
c
c  Formula:
c
c    [X]_N = X * ( X - 1 ) * ( X - 2 ) * ... * ( X - N + 1 ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer X, the argument of the falling factorial
c    function.
c
c    Input, integer N, the order of the falling factorial function.
c    If N = 0, FALL = 1, if N = 1, FALL = X.  Note that if N is
c    negative, a "rising" factorial will be computed.
c
c    Output, integer FALL, the falling factorial function.
c
      implicit none

      integer arg
      integer fall
      integer i
      integer n
      integer x

      fall = 1

      arg = x

      if ( 0 .lt. n ) then

        do i = 1, n
          fall = fall * arg
          arg = arg - 1
        end do

      else if ( n .lt. 0 ) then

        do i = n, -1
          fall = fall * arg
          arg = arg + 1
        end do

      end if

      return
      end
      subroutine gray_code_check ( n, t, ierror )

c*********************************************************************72
c
cc GRAY_CODE_CHECK checks a Gray code element.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of digits in each element.
c    N must be positive.
c
c    Input, integer T(N), an element of the Gray code.
c    Each entry T(I) is either 0 or 1.
c
c    Output, integer IERROR, error flag.
c    0, no error, T represents a Gray code element.
c    -1, N is not positive.
c    I, error, T(I) is an illegal value for a Gray code element.
c
      implicit none

      integer n

      integer i
      integer ierror
      integer t(n)

      ierror = 0

      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRAY_CODE_CHECK - Fatal error!'
        write ( *, '(a)' ) '  N .lt. 1.'
        stop
      end if

      do i = 1, n

        if ( t(i) .ne. 0 .and. t(i) .ne. 1 ) then
          ierror = i
          return
        end if

      end do

      return
      end
      subroutine gray_code_enum ( n, ngray )

c*********************************************************************72
c
cc GRAY_CODE_ENUM enumerates the Gray codes on N digits.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of digits in each element.
c    N must be nonnegative.
c
c    Output, integer NGRAY, the number of distinct elements.
c
      implicit none

      integer n
      integer ngray

      ngray = 2**n

      return
      end
      subroutine gray_code_rank ( n, t, rank )

c*********************************************************************72
c
cc GRAY_CODE_RANK computes the rank of a Gray code element.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of digits in each element.
c    N must be positive.
c
c    Input, integer T(N), an element of the Gray code.
c    Each entry T(I) is either 0 or 1.
c
c    Output, integer RANK, the rank of the element.
c
      implicit none

      integer n

      integer b
      integer i
      integer ierror
      integer rank
      integer t(n)
c
c  Check.
c
      call gray_code_check ( n, t, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRAY_CODE_RANK - Fatal error!'
        write ( *, '(a)' ) '  The input array is illegal.'
        write ( *, '(a,i8)' ) '  IERROR = ', ierror
        stop
      end if

      rank = 0
      b = 0

      do i = n - 1, 0, -1

        if ( t(n-i) .ne. 0 ) then
          b = 1 - b
        end if

        if ( b .eq. 1 ) then
          rank = rank + 2 ** i
        end if

      end do

      return
      end
      subroutine gray_code_successor ( n, t, rank )

c*********************************************************************72
c
cc GRAY_CODE_SUCCESSOR computes the binary reflected Gray code successor.
c
c  Example:
c
c    000, 001, 011, 010, 110, 111, 101, 100,
c    after which the sequence repeats.
c
c  Discussion:
c
c    In the original code, the successor of the element that has an
c    initial 1 followed by N-1 zeroes is undefined.  In this version,
c    the successor is the element with N zeroes.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of digits in each element.
c    N must be positive.
c
c    Input/output, integer T(N).
c    On input, T contains an element of the Gray code, that is,
c    each entry T(I) is either 0 or 1.
c    On output, T contains the successor to the input value; this
c    is an element of the Gray code, which differs from the input
c    value in a single position.
c
c    Input/output, integer RANK, the rank.
c    If RANK = -1 on input, then the routine understands that this is
c    the first call, and that the user wishes the routine to supply
c    the first element in the ordering, which has RANK = 0.
c    In general, the input value of RANK is increased by 1 for output,
c    unless the very last element of the ordering was input, in which
c    case the output value of RANK is 0.
c
      implicit none

      integer n

      integer i
      integer i4vec_sum
      integer ierror
      integer rank
      integer t(n)
      integer weight
c
c  Return the first element.
c
      if ( rank .eq. -1 ) then

        do i = 1, n
          t(i) = 0
        end do
        rank = 0
        return

      end if
c
c  Check.
c
      call gray_code_check ( n, t, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRAY_CODE_SUCCESSOR - Fatal error!'
        write ( *, '(a)' ) '  The input array is illegal.'
        write ( *, '(a,i8)' ) '  IERROR = ', ierror
        stop
      end if

      weight = i4vec_sum ( n, t )

      if ( mod ( weight, 2 ) .eq. 0 ) then

        if ( t(n) .eq. 0 ) then
          t(n) = 1
        else
          t(n) = 0
        end if

        rank = rank + 1
        return

      else

        do i = n, 2, -1
          if ( t(i) .eq. 1 ) then
            if ( t(i-1) .eq. 0 ) then
              t(i-1) = 1
            else
              t(i-1) = 0
            end if
            rank = rank + 1
            return
          end if
        end do
c
c  The final element was input.
c  Return the first element.
c
        do i = 1, n
          t(i) = 0
        end do
        rank = 0

      end if

      return
      end
      subroutine gray_code_unrank ( rank, n, t )

c*********************************************************************72
c
cc GRAY_CODE_UNRANK computes the Gray code element of given rank.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer RANK, the rank of the element.
c    0 .le. RANK .le. 2**N.
c
c    Input, integer N, the number of digits in each element.
c    N must be positive.
c
c    Output, integer T(N), the element of the Gray code which has
c    the given rank.
c
      implicit none

      integer n

      integer b
      integer bprime
      integer i
      integer ngray
      integer rank
      integer rank_copy
      integer t(n)
c
c  Check.
c
      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRAY_CODE_UNRANK - Fatal error!'
        write ( *, '(a)' ) '  Input N is illegal.'
        stop
      end if

      call gray_code_enum ( n, ngray )

      if ( rank .lt. 0 .or. ngray .lt. rank ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRAY_CODE_UNRANK - Fatal error!'
        write ( *, '(a)' ) '  The input rank is illegal.'
        stop
      end if

      rank_copy = rank
      do i = 1, n
        t(i) = 0
      end do
      bprime = 0

      do i = n - 1, 0, -1

        b = rank_copy / 2**i

        if ( b .ne. bprime ) then
          t(n-i) = 1
        end if

        bprime = b
        rank_copy = rank_copy - b * 2 ** i

      end do

      return
      end
      function i4_factorial ( n )

c*********************************************************************72
c
cc I4_FACTORIAL computes the factorial of N.
c
c  Discussion:
c
c    factorial ( N ) = product ( 1 <= I <= N ) I
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 June 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the argument of the factorial function.
c    If N is less than 1, the function value is returned as 1.
c    0 <= N <= 13 is required.
c
c    Output, integer I4_FACTORIAL, the factorial of N.
c
      implicit none

      integer i
      integer i4_factorial
      integer n

      i4_factorial = 1

      if ( 13 .lt. n ) then
        i4_factorial = - 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_FACTORIAL - Fatal error!'
        write ( *, '(a)' )
     &  '  I4_FACTORIAL(N) cannot be computed as an integer'
        write ( *, '(a)' ) '  for 13 < N.'
        write ( *, '(a,i8)' ) '  Input value N = ', n
        stop
      end if

      do i = 1, n
        i4_factorial = i4_factorial * i
      end do

      return
      end
      subroutine i4_factorial_values ( n_data, n, fn )

c*********************************************************************72
c
cc I4_FACTORIAL_VALUES returns values of the factorial function.
c
c  Discussion:
c
c    0! = 1
c    I! = Product ( 1 <= J <= I ) I
c
c    In Mathematica, the function can be evaluated by:
c
c      n!
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 March 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Milton Abramowitz, Irene Stegun,
c    Handbook of Mathematical Functions,
c    National Bureau of Standards, 1964,
c    ISBN: 0-486-61272-4,
c    LC: QA47.A34.
c
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Cambridge University Press, 1999,
c    ISBN: 0-521-64314-7,
c    LC: QA76.95.W65.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, integer N, the argument of the function.
c
c    Output, integer FN, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 13 )

      integer fn_vec(n_max)
      integer fn
      integer n
      integer n_data
      integer n_vec(n_max)

      save fn_vec
      save n_vec

      data fn_vec /
     &          1,
     &          1,
     &          2,
     &          6,
     &         24,
     &        120,
     &        720,
     &       5040,
     &      40320,
     &     362880,
     &    3628800,
     &   39916800,
     &  479001600 /
      data n_vec /
     &   0,  1,  2,  3,
     &   4,  5,  6,  7,
     &   8,  9, 10, 11,
     &  12 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        n = 0
        fn = 0
      else
        n = n_vec(n_data)
        fn = fn_vec(n_data)
      end if

      return
      end
      function i4_huge ( )

c*********************************************************************72
c
cc I4_HUGE returns a "huge" I4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer I4_HUGE, a huge number.
c
      implicit none

      integer i4_huge

      i4_huge = 2147483647

      return
      end
      function i4_uniform ( a, b, seed )

c*********************************************************************72
c
cc I4_UNIFORM returns a scaled pseudorandom I4.
c
c  Discussion:
c
c    An I4 is an integer value.
c
c    The pseudorandom number should be uniformly distributed
c    between A and B.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Springer Verlag, pages 201-202, 1983.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley Interscience, page 95, 1998.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, pages 362-376, 1986.
c
c    Peter Lewis, Allen Goodman, James Miller
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, pages 136-143, 1969.
c
c  Parameters:
c
c    Input, integer A, B, the limits of the interval.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, integer I4_UNIFORM, a number between A and B.
c
      implicit none

      integer a
      integer b
      integer i4_uniform
      integer k
      real r
      integer seed
      integer value

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_UNIFORM - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + 2147483647
      end if

      r = real ( seed ) * 4.656612875E-10
c
c  Scale R to lie between A-0.5 and B+0.5.
c
      r = ( 1.0E+00 - r ) * ( real ( min ( a, b ) ) - 0.5E+00 )
     &  +             r   * ( real ( max ( a, b ) ) + 0.5E+00 )
c
c  Use rounding to convert R to an integer between A and B.
c
      value = nint ( r )

      value = max ( value, min ( a, b ) )
      value = min ( value, max ( a, b ) )

      i4_uniform = value

      return
      end
      subroutine i4vec_backtrack ( n, x, indx, k, nstack, stack, 
     &  maxstack, ncan )

c*********************************************************************72
c
cc I4VEC_BACKTRACK supervises a backtrack search for an I4VEC.
c
c  Discussion:
c
c    The routine tries to construct an integer vector one index at a time,
c    using possible candidates as supplied by the user.
c
c    At any time, the partially constructed vector may be discovered to be
c    unsatisfactory, but the routine records information about where the
c    last arbitrary choice was made, so that the search can be
c    carried out efficiently, rather than starting out all over again.
c
c    First, call the routine with INDX = 0 so it can initialize itself.
c
c    Now, on each return from the routine, if INDX is:
c      1, you've just been handed a complete candidate vector;
c         Admire it, analyze it, do what you like.
c      2, please determine suitable candidates for position X(K).
c         Return the number of candidates in NCAN(K), adding each
c         candidate to the end of STACK, and increasing NSTACK.
c      3, you're done.  Stop calling the routine;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 July 2000
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of positions to be filled in
c    the vector.
c
c    Input/output, integer X(N), the partial or complete candidate
c    vector.
c
c    Input/output, integer INDX, a communication flag.
c    On input,
c      0 to start a search.
c    On output:
c      1, a complete output vector has been determined and returned in X(1:N);
c      2, candidates are needed for position X(K);
c      3, no more possible vectors exist.
c
c    Output, integer K, if INDX=2, the current vector index
c    being considered.
c
c    Input/output, integer NSTACK, the current length of the stack.
c
c    Input, integer STACK(MAXSTACK), a list of all current
c    candidates for all positions 1 through K.
c
c    Input, integer MAXSTACK, the maximum length of the stack.
c
c    Input/output, integer NCAN(N), lists the current number of
c    candidates for positions 1 through K.
c
      implicit none

      integer n
      integer maxstack

      integer indx
      integer k
      integer ncan(n)
      integer nstack
      integer stack(maxstack)
      integer x(n)
c
c  If this is the first call, request a candidate for position 1.
c
      if ( indx .eq. 0 ) then
        k = 1
        nstack = 0
        indx = 2
        return
      end if
c
c  Examine the stack.
c
10    continue
c
c  If there are candidates for position K, take the first available
c  one off the stack, and increment K.
c
c  This may cause K to reach the desired value of N, in which case
c  we need to signal the user that a complete set of candidates
c  is being returned.
c
        if ( 0 .lt. ncan(k) ) then

          x(k) = stack(nstack)
          nstack = nstack - 1

          ncan(k) = ncan(k) - 1

          if ( k .ne. n ) then
            k = k + 1
            indx = 2
          else
            indx = 1
          end if

          go to 20
c
c  If there are no candidates for position K, then decrement K.
c  If K is still positive, repeat the examination of the stack.
c
        else

          k = k - 1

          if ( k .le. 0 ) then
            indx = 3
            go to 20
          end if

        end if

      go to 10

20    continue

      return
      end
      subroutine i4vec_indicator ( n, a )

c*********************************************************************72
c
cc I4VEC_INDICATOR sets an I4VEC to the indicator vector.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of elements of A.
c
c    Output, integer A(N), the array to be initialized.
c
      implicit none

      integer n

      integer a(n)
      integer i

      do i = 1, n
        a(i) = i
      end do

      return
      end
      subroutine i4vec_max ( n, a, amax )

c*********************************************************************72
c
cc I4VEC_MAX computes the maximum element of an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input, integer A(N), the array.
c
c    Output, integer AMAX, the value of the largest entry.
c
      implicit none

      integer n

      integer a(n)
      integer amax
      integer i

      amax = a(1)
      do i = 2, n
        amax = max ( amax, a(i) )
      end do

      return
      end
      subroutine i4vec_part1 ( n, npart, x )

c*********************************************************************72
c
cc I4VEC_PART1 partitions an integer N into NPART parts.
c
c  Example:
c
c    Input:
c
c      N = 17, NPART = 5
c
c    Output:
c
c      X = ( 13, 1, 1, 1, 1 ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the integer to be partitioned.  N
c    may be positive, zero, or negative.
c
c    Input, integer NPART, the number of entries in the array.
c    1 .le. NPART .le. N.
c
c    Output, integer X(NPART), the partition of N.  The entries of
c    X add up to N.  X(1) = N + 1 - NPART, and all other entries
c    are equal to 1.
c
      implicit none

      integer npart

      integer i
      integer n
      integer x(npart)

      if ( npart .lt. 1 .or. n .lt. npart ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4VEC_PART1 - Fatal error!'
        write ( *, '(a)' ) '  The input value of NPART is illegal.'
        stop
      end if

      x(1) = n + 1 - npart
      do i = 2, npart
        x(i) = 1
      end do

      return
      end
      subroutine i4vec_part2 ( n, npart, x )

c*********************************************************************72
c
cc I4VEC_PART2 partitions an integer N into NPART nearly equal parts.
c
c  Example:
c
c    Input:
c
c      N = 17, NPART = 5
c
c    Output:
c
c      X = ( 4, 4, 3, 3, 3 ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the integer to be partitioned.  N
c    may be positive, zero, or negative.
c
c    Input, integer NPART, the number of entries in the array.
c    1 .le. NPART
c
c    Output, integer X(NPART), the partition of N.  The entries of
c    X add up to N.  The entries of X are either all equal, or
c    differ by at most 1.  The entries of X all have the same sign
c    as N, and the "largest" entries occur first.
c
      implicit none

      integer npart

      integer i
      integer j
      integer n
      integer x(npart)

      if ( npart .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4VEC_PART2 - Fatal error!'
        write ( *, '(a)' ) '  The input value of NPART is illegal.'
        stop
      end if

      do i = 1, npart
        x(i) = 0
      end do

      if ( 0 .lt. n ) then

        j = 1
        do i = 1, n
          x(j) = x(j) + 1
          j = j + 1
          if ( npart .lt. j ) then
            j = 1
          end if
        end do

      else if ( n .lt. 0 ) then

        j = 1
        do i = n, - 1
          x(j) = x(j) - 1
          j = j + 1
          if ( npart .lt. j ) then
            j = 1
          end if
        end do

      end if

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
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
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
c    Input, character*(*) TITLE, a title.
c
      implicit none

      integer n

      integer a(n)
      integer i
      character*(*) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,a,1x,i12)' ) i, ':', a(i)
      end do

      return
      end
      subroutine i4vec_reverse ( n, a )

c*********************************************************************72
c
cc I4VEC_REVERSE reverses the elements of an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    In FORTRAN90, call I4VEC_REVERSE is equivalent to:
c
c      A(1:N) = A(N:1:-1)
c
c  Example:
c
c    Input:
c
c      N = 5,
c      A = ( 11, 12, 13, 14, 15 ).
c
c    Output:
c
c      A = ( 15, 14, 13, 12, 11 ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input/output, integer A(N), the array to be reversed.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer t

      do i = 1, n / 2
        t        = a(i)
        a(i)     = a(n+1-i)
        a(n+1-i) = t
      end do

      return
      end
      subroutine i4vec_search_binary_a ( n, a, b, indx )

c*********************************************************************72
c
cc I4VEC_SEARCH_BINARY_A searches an ascending sorted I4VEC for a value.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    Binary search is used.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Algorithm 1.9,
c    Combinatorial Algorithms,
c    CRC Press, 1998, page 26.
c
c  Parameters:
c
c    Input, integer N, the number of elements in the vector.
c
c    Input, integer A(N), the array to be searched.  A must
c    be sorted in ascending order.
c
c    Input, integer B, the value to be searched for.
c
c    Output, integer INDX, the result of the search.
c    -1, B does not occur in A.
c    I, A(I) = B.
c
      implicit none

      integer n

      integer a(n)
      integer b
      integer high
      integer indx
      integer low
      integer mid

      indx = - 1

      low = 1
      high = n

10    continue

      if ( low .le. high ) then

        mid = ( low + high ) / 2

        if ( a(mid) .eq. b ) then
          indx = mid
          go to 20
        else if ( a(mid) .lt. b ) then
          low = mid + 1
        else if ( b .lt. a(mid) ) then
          high = mid - 1
        end if

        go to 10

      end if

20    continue

      return
      end
      subroutine i4vec_search_binary_d ( n, a, b, indx )

c*********************************************************************72
c
cc I4VEC_SEARCH_BINARY_D searches a descending sorted I4VEC for a value.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    Binary search is used.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Algorithm 1.9,
c    Combinatorial Algorithms,
c    CRC Press, 1998, page 26.
c
c  Parameters:
c
c    Input, integer N, the number of elements in the vector.
c
c    Input, integer A(N), the array to be searched.  A must
c    be sorted in descending order.
c
c    Input, integer B, the value to be searched for.
c
c    Output, integer INDX, the result of the search.
c    -1, B does not occur in A.
c    I, A(I) = B.
c
      implicit none

      integer n

      integer a(n)
      integer b
      integer high
      integer indx
      integer low
      integer mid

      indx = - 1

      low = 1
      high = n

10    continue

      if ( low .le. high ) then

        mid = ( low + high ) / 2

        if ( a(mid) .eq. b ) then
          indx = mid
          go to 20
        else if ( b .lt. a(mid) ) then
          low = mid + 1
        else if ( a(mid) .lt. b ) then
          high = mid - 1
        end if

        go to 10

      end if

20    continue

      return
      end
      subroutine i4vec_sort_insert_a ( n, a )

c*********************************************************************72
c
cc I4VEC_SORT_INSERT_A uses an ascending insertion sort on an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Algorithm 1.1,
c    Combinatorial Algorithms,
c    CRC Press, 1998, page 11.
c
c  Parameters:
c
c    Input, integer N, the number of items in the vector.
c    N must be positive.
c
c    Input/output, integer A(N).
c    On input, A contains data to be sorted.
c    On output, the entries of A have been sorted in ascending order.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer j
      integer x

      do i = 2, n

        x = a(i)

        j = i - 1

10      continue

        if ( 1 .le. j ) then

          if ( a(j) .le. x ) then
            go to 20
          end if

          a(j+1) = a(j)
          j = j - 1

          go to 10

        end if

20      continue

        a(j+1) = x

      end do

      return
      end
      subroutine i4vec_sort_insert_d ( n, a )

c*********************************************************************72
c
cc I4VEC_SORT_INSERT_D uses a descending insertion sort on an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Algorithm 1.1,
c    Combinatorial Algorithms,
c    CRC Press, 1998, page 11.
c
c  Parameters:
c
c    Input, integer N, the number of items in the vector.
c    N must be positive.
c
c    Input/output, integer A(N).
c    On input, A contains data to be sorted.
c    On output, the entries of A have been sorted in descending order.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer j
      integer x

      do i = 2, n

        x = a(i)

        j = i - 1

10      continue

        if ( 1 .le. j ) then

          if ( x .le. a(j) ) then
            go to 20
          end if

          a(j+1) = a(j)
          j = j - 1
          go to 10

        end if

20      continue

        a(j+1) = x

      end do

      return
      end
      function i4vec_sum ( n, a )

c*********************************************************************72
c
cc I4VEC_SUM returns the sum of the entries of an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    In FORTRAN90, this facility is offered by the built in
c    SUM function:
c
c      I4VEC_SUM ( N, A ) = SUM ( A(1:N) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input, integer A(N), the array.
c
c    Output, integer I4VEC_SUM, the sum of the entries.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer i4vec_sum

      i4vec_sum = 0

      do i = 1, n
        i4vec_sum = i4vec_sum + a(i)
      end do

      return
      end
      subroutine knapsack_01 ( n, mass_limit, p, w, x, mass, profit )

c*********************************************************************72
c
cc KNAPSACK_01 solves the 0/1 knapsack problem.
c
c  Discussion:
c
c    The 0/1 knapsack problem is as follows:
c
c      Given:
c        a set of N objects,
c        a profit P(I) and weight W(I) associated with each object,
c        and a weight limit MASS_LIMIT,
c      Determine:
c        a set of choices X(I) which are 0 or 1, that maximizes the profit
c          P = Sum ( 1 .le. I .le. N ) P(I) * X(I)
c        subject to the constraint
c          Sum ( 1 .le. I .le. N ) W(I) * X(I) .le. MASS_LIMIT.
c
c    This routine assumes that the objects have already been sorted
c    in order of decreasing "profit density", P(I)/W(I).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 July 2000
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of objects.
c
c    Input, double precision MASS_LIMIT, the weight limit of the
c    chosen objects.
c
c    Input/output, double precision P(N), the "profit" or value of each object.
c    P is assumed to be nonnegative.
c
c    Input/output, double precision W(N), the "weight" or cost of each object.
c    W is assumed to be  nonnegative.
c
c    Output, double precision X(N), the choice function for the objects.
c    0, the object was not taken.
c    1, the object was taken.
c
c    Output, double precision MASS, the total mass of the objects taken.
c
c    Output, double precision PROFIT, the total profit of the objects taken.
c
      implicit none

      integer maxstack
      parameter ( maxstack = 100 )

      integer n

      integer i
      integer indx
      integer k
      double precision mass
      double precision mass_1
      double precision mass_2
      double precision mass_best
      double precision mass_limit
      double precision mass_remaining
      integer ncan(n)
      integer nstack
      double precision p(n)
      double precision profit
      double precision profit_1
      double precision profit_2
      double precision profit_best
      double precision r8vec_dot_product
      double precision stack(maxstack)
      double precision w(n)
      double precision x(n)
      double precision x_best(n)

      nstack = 0
c
c  Initialize the "best so far" data.
c
      do i = 1, n
        x_best(i) = 0.0D+00
      end do
      profit_best = 0.0D+00
      mass_best = 0
c
c  Begin the backtracking solution.
c
      indx = 0

10    continue

        call r8vec_backtrack ( n, x, indx, k, nstack, stack, 
     &    maxstack, ncan )
c
c  Got a new candidate.  Compare it to the best so far.
c
        if ( indx .eq. 1 ) then

          profit = r8vec_dot_product ( n, p, x )
          mass = r8vec_dot_product ( n, w, x )

          if ( profit_best .lt. profit .or. 
     &       ( profit .eq. profit_best .and. 
     &         mass .lt. mass_best ) ) then
            profit_best = profit
            mass_best = mass
            do i = 1, n
              x_best(i) = x(i)
            end do
          end if
c
c  Need candidates for X(K).
c
c  X(K) = 1 is possible if:
c
c    * adding W(K) to our mass doesn't put us over our mass limit;
c    * and adding P(K) to our current profit, and taking the best we
c      could get using rational X for the remainder would put us over
c      our current best.
c
c  X(K) = 0 is always possible.
c
        else if ( indx .eq. 2 ) then

          ncan(k) = 0

          mass_1 = r8vec_dot_product ( k - 1, w, x ) + w(k)

          if ( mass_1 .le. mass_limit ) then

            mass_remaining = mass_limit - mass_1

            profit_1 = r8vec_dot_product ( k - 1, p, x ) + p(k)

            if ( k .lt. n ) then
              call knapsack_rational ( n - k, mass_remaining, p(k+1), 
     &          w(k+1), x(k+1), mass_2, profit_2 )
            else
              profit_2 = 0.0D+00
            end if

            if ( profit_best .lt. profit_1 + profit_2 ) then
              if ( maxstack .le. nstack ) then
                write ( *, '(a)' ) ' '
                write ( *, '(a)' ) 'KNAPSACK_01 - Fatal error!'
                write ( *, '(a)' ) '  Exceeded stack space.'
                return
              end if
              ncan(k) = ncan(k) + 1
              nstack = nstack + 1
              stack(nstack) = 1.0D+00
            end if

          end if

          if ( maxstack .le. nstack ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'KNAPSACK_01 - Fatal error!'
            write ( *, '(a)' ) '  Exceeded stack space.'
            return
          end if

          ncan(k) = ncan(k) + 1
          nstack = nstack + 1
          stack(nstack) = 0.0D+00
c
c  Done.  Return the best solution.
c
        else

          profit = profit_best
          mass = mass_best
          do i = 1, n
            x(i) = x_best(i)
          end do
          go to 20

        end if

      go to 10

20    continue

      return
      end

      subroutine knapsack_rational ( n, mass_limit, p, w, x, mass, 
     &  profit )

c*********************************************************************72
c
cc KNAPSACK_RATIONAL solves the rational knapsack problem.
c
c  Discussion:
c
c    The rational knapsack problem is a generalization of the 0/1 knapsack
c    problem.  It is mainly used to derive a bounding function for the
c    0/1 knapsack problem.
c
c    The 0/1 knapsack problem is as follows:
c
c      Given:
c        a set of N objects,
c        a profit P(I) and weight W(I) associated with each object,
c        and a weight limit MASS_LIMIT,
c      Determine:
c        a set of choices X(I) which are 0 or 1, that maximizes the profit
c          P = Sum ( 1 .le. I .le. N ) P(I) * X(I)
c        subject to the constraint
c          Sum ( 1 .le. I .le. N ) W(I) * X(I) .le. MASS_LIMIT.
c
c    By contrast, the rational knapsack problem allows the values X(I)
c    to be any value between 0 and 1.  A solution for the rational knapsack
c    problem is known.  Arrange the objects in order of their "profit density"
c    ratios P(I)/W(I), and then take in order as many of these as you can.
c    If you still have "room" in the weight constraint, then you should
c    take the maximal fraction of the very next object, which will complete
c    your weight limit, and maximize your profit.
c
c    If should be obvious that, given the same data, a solution for
c    the rational knapsack problem will always have a profit that is
c    at least as high as for the 0/1 problem.  Since the rational knapsack
c    maximum profit is easily computed, this makes it a useful bounding
c    function.
c
c    Note that this routine assumes that the objects have already been
c    arranged in order of the "profit density".
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 July 2000
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of objects.
c
c    Input, double precision MASS_LIMIT, the weight limit of the
c    chosen objects.
c
c    Input, double precision P(N), the "profit" or value of each object.
c    The entries of P are assumed to be nonnegative.
c
c    Input, double precision W(N), the "weight" or cost of each object.
c    The entries of W are assumed to be nonnegative.
c
c    Output, double precision X(N), the choice function for the objects.
c    0.0, the object was not taken.
c    1.0, the object was taken.
c    R, where 0 .lt. R .lt. 1, a fractional amount of the object was taken.
c
c    Output, double precision MASS, the total mass of the objects taken.
c
c    Output, double precision PROFIT, the total profit of the objects taken.
c
      implicit none

      integer n

      integer i
      double precision mass
      double precision mass_limit
      double precision p(n)
      double precision profit
      double precision w(n)
      double precision x(n)

      mass = 0.0D+00
      profit = 0.0D+00

      do i = 1, n

        if ( mass_limit .le. mass ) then
          x(i) = 0.0D+00
        else if ( mass + w(i) .le. mass_limit ) then
          x(i) = 1.0D+00
          mass = mass + w(i)
          profit = profit + p(i)
        else
          x(i) = ( mass_limit - mass ) / w(i)
          mass = mass_limit
          profit = profit + p(i) * x(i)
        end if

      end do

      return
      end
      subroutine knapsack_reorder ( n, p, w )

c*********************************************************************72
c
cc KNAPSACK_REORDER reorders the knapsack data by "profit density".
c
c  Discussion:
c
c    This routine must be called to rearrange the data before calling
c    routines that handle a knapsack problem.
c
c    The "profit density" for object I is defined as P(I)/W(I).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 July 2000
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of objects.
c
c    Input/output, double precision P(N), the "profit" or value of each object.
c
c    Input/output, double precision W(N), the "weight" or cost of each object.
c
      implicit none

      integer n

      integer i
      integer j
      double precision p(n)
      double precision t
      double precision w(n)
c
c  Rearrange the objects in order of "profit density".
c
      do i = 1, n
        do j = i + 1, n
          if ( p(i) * w(j) .lt. p(j) * w(i) ) then

            t    = p(i)
            p(i) = p(j)
            p(j) = t

            t    = w(i)
            w(i) = w(j)
            w(j) = t

          end if
        end do
      end do

      return
      end
      subroutine ksubset_colex_check ( k, n, t, ierror )

c*********************************************************************72
c
cc KSUBSET_COLEX_CHECK checks a K subset in colex form.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer K, the number of elements each K subset must
c    have. 1 .le. K .le. N.
c
c    Input, integer N, the number of elements in the master set.
c    N must be positive.
c
c    Input, integer T(K), describes a K subset.  T(I) is the I-th
c    element of the K subset.  The elements must be listed in
c    DESCENDING order.
c
c    Output, integer IERROR, error flag.
c    0, no error.
c    -1, N is not positive.
c    -2, K is not positive.
c    I, entry I is illegal.
c
      implicit none

      integer k

      integer i
      integer ierror
      integer n
      integer t(k)
      integer tmax

      ierror = 0

      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'KSUBSET_COLEX_CHECK - Fatal error!'
        write ( *, '(a)' ) '  N .lt. 1.'
        stop
      end if

      if ( k .lt. 1 .or. n .lt. k ) then
        ierror = -2
        return
      end if

      tmax = n + 1

      do i = 1, k

        if ( t(i) .le. 0 .or. tmax .le. t(i) ) then
          ierror = i
          return
        end if

        tmax = t(i)

      end do

      return
      end
      subroutine ksubset_colex_rank ( k, n, t, rank )

c*********************************************************************72
c
cc KSUBSET_COLEX_RANK computes the colex rank of a K subset.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer K, the number of elements each K subset must
c    have.  1 .le. K .le. N.
c
c    Input, integer N, the number of elements in the master set.
c    N must be positive.
c
c    Input, integer T(K), describes a K subset.  T(I) is the I-th
c    element of the K subset.  The elements must be listed in DESCENDING order.
c
c    Output, integer RANK, the rank of the subset.
c
      implicit none

      integer k

      integer binomial
      integer i
      integer ierror
      integer n
      integer rank
      integer t(k)
c
c  Check.
c
      call ksubset_colex_check ( k, n, t, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'KSUBSET_COLEX_CHECK - Fatal error!'
        write ( *, '(a)' ) '  The input array is illegal.'
        write ( *, '(a,i8)' ) '  IERROR = ', ierror
        stop
      end if

      rank = 0

      do i = 1, k
        rank = rank + binomial ( t(i) - 1, k + 1 - i )
      end do

      return
      end
      subroutine ksubset_colex_successor ( k, n, t, rank )

c*********************************************************************72
c
cc KSUBSET_COLEX_SUCCESSOR computes the K subset colex successor.
c
c  Discussion:
c
c    In the original code, there is a last element with no successor.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer K, the number of elements each K subset must
c    have.  1 .le. K .le. N.
c
c    Input, integer N, the number of elements in the master set.
c    N must be positive.
c
c    Input/output, integer T(K), describes a K subset.  T(I) is the
c    I-th element.  The elements must be listed in DESCENDING order.
c    On input, T describes a K subset.
c    On output, T describes the next K subset in the ordering.
c    If the input T was the last in the ordering, then the output T
c    will be the first.
c
c    Input/output, integer RANK, the rank.
c    If RANK = -1 on input, then the routine understands that this is
c    the first call, and that the user wishes the routine to supply
c    the first element in the ordering, which has RANK = 0.
c    In general, the input value of RANK is increased by 1 for output,
c    unless the very last element of the ordering was input, in which
c    case the output value of RANK is 0.
c
      implicit none

      integer k

      integer i
      integer ierror
      integer n
      integer rank
      integer t(k)
c
c  Return the first element.
c
      if ( rank .eq. -1 ) then
        do i = 1, k
          t(i) = k + 1 - i
        end do
        rank = 0
        return
      end if
c
c  Check.
c
      call ksubset_colex_check ( k, n, t, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'KSUBSET_COLEX_SUCCESSOR - Fatal error!'
        write ( *, '(a)' ) '  The input array is illegal.'
        write ( *, '(a,i8)' ) '  IERROR = ', ierror
        stop
      end if

      do i = k - 1, 1, -1
        if ( t(k+1-i) + 1 .lt. t(k-i) ) then
          t(k+1-i) = t(k+1-i) + 1
          rank = rank + 1
          return
        end if
      end do

      if ( t(1) .lt. n ) then
        t(1) = t(1) + 1
        do i = 1, k - 1
          t(k+1-i) = i
        end do
        rank = rank + 1
        return
      end if
c
c  The last K subset was input.
c  Return the first one.
c
      do i = 1, k
        t(i) = k + 1 - i
      end do

      rank = 0

      return
      end
      subroutine ksubset_colex_unrank ( rank, k, n, t )

c*********************************************************************72
c
cc KSUBSET_COLEX_UNRANK computes the K subset of given colex rank.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer RANK, the rank of the K subset.
c
c    Input, integer K, the number of elements each K subset must
c    have.  1 .le. K .le. N.
c
c    Input, integer N, the number of elements in the master set.
c    N must be positive.
c
c    Output, integer T(K), describes the K subset of the given
c    rank.  T(I) is the I-th element.  The elements must be listed in
c    DESCENDING order.
c
      implicit none

      integer k

      integer binomial
      integer i
      integer n
      integer nksub
      integer rank
      integer rank_copy
      integer t(k)
      integer x
c
c  Check.
c
      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'KSUBSET_COLEX_UNRANK - Fatal error!'
        write ( *, '(a)' ) '  Input N is illegal.'
        stop
      end if

      if ( k .lt. 1 .or. n .lt. k ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'KSUBSET_COLEX_UNRANK - Fatal error!'
        write ( *, '(a)' ) '  Input K is illegal.'
        stop
      end if

      call ksubset_enum ( k, n, nksub )

      if ( rank .lt. 0 .or. nksub .lt. rank ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'KSUBSET_COLEX_UNRANK - Fatal error!'
        write ( *, '(a)' ) '  The input rank is illegal.'
        stop
      end if

      rank_copy = rank

      x = n

      do i = 1, k

10      continue

        if ( rank_copy .lt. binomial ( x, k + 1 - i ) ) then
          x = x - 1
          go to 10
        end if

        t(i) = x + 1
        rank_copy = rank_copy - binomial ( x, k + 1 - i )

      end do

      return
      end
      subroutine ksubset_enum ( k, n, nksub )

c*********************************************************************72
c
cc KSUBSET_ENUM enumerates the K element subsets of an N set.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer K, the number of elements each K subset must
c    have. 0 .le. K .le. N.
c
c    Input, integer N, the number of elements in the master set.
c    0 .le. N.
c
c    Output, integer NKSUB, the number of distinct elements.
c
      implicit none

      integer binomial
      integer k
      integer n
      integer nksub

      nksub = binomial ( n, k )

      return
      end
      subroutine ksubset_lex_check ( k, n, t, ierror )

c*********************************************************************72
c
cc KSUBSET_LEX_CHECK checks a K subset in lex form.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer K, the number of elements each K subset must
c    have. 1 .le. K .le. N.
c
c    Input, integer N, the number of elements in the master set.
c    N must be positive.
c
c    Input, integer T(K), describes a K subset.  T(I) is the I-th
c    element of the K subset.  The elements must be listed in
c    DESCENDING order.
c
c    Output, integer IERROR, error flag.
c    0, no error.
c    -1, N is illegal.
c    -2, K is illegal.
c    I, entry I is illegal.
c
      implicit none

      integer k

      integer i
      integer ierror
      integer n
      integer t(k)
      integer tmin

      ierror = 0

      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'KSUBSET_LEX_CHECK - Fatal error!'
        write ( *, '(a)' ) '  N .lt. 1.'
        stop
      end if

      if ( k .lt. 1 .or. n .lt. k ) then
        ierror = -2
        return
      end if

      tmin = 0

      do i = 1, k

        if ( t(i) .le. tmin .or. n .lt. t(i) ) then
          ierror = i
          return
        end if

        tmin = t(i)

      end do

      return
      end
      subroutine ksubset_lex_rank ( k, n, t, rank )

c*********************************************************************72
c
cc KSUBSET_LEX_RANK computes the lexicographic rank of a K subset.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer K, the number of elements each K subset must
c    have.  1 .le. K .le. N.
c
c    Input, integer N, the number of elements in the master set.
c    N must be positive.
c
c    Input, integer T(K), describes a K subset.  T(I) is the I-th
c    element.  The elements must be listed in ascending order.
c
c    Output, integer RANK, the rank of the K subset.
c
      implicit none

      integer k

      integer binomial
      integer i
      integer ierror
      integer j
      integer n
      integer rank
      integer t(k)
      integer tim1
c
c  Check.
c
      call ksubset_lex_check ( k, n, t, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'KSUBSET_LEX_RANK - Fatal error!'
        write ( *, '(a)' ) '  The input array is illegal.'
        stop
      end if

      rank = 0

      do i = 1, k

        if ( i .eq. 1 ) then
          tim1 = 0
        else
          tim1 = t(i-1)
        end if

        if ( tim1 + 1 .le. t(i) - 1 ) then
          do j = tim1 + 1, t(i) - 1
            rank = rank + binomial ( n - j, k - i )
          end do
        end if

      end do

      return
      end
      subroutine ksubset_lex_successor ( k, n, t, rank )

c*********************************************************************72
c
cc KSUBSET_LEX_SUCCESSOR computes the K subset lexicographic successor.
c
c  Discussion:
c
c    In the original code, there is a last element with no successor.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 February 2001
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer K, the number of elements each K subset must
c    have. 1 .le. K .le. N.
c
c    Input, integer N, the number of elements in the master set.
c    N must be positive.
c
c    Input/output, integer T(K), describes a K subset.  T(I) is
c    the I-th element.  The elements must be listed in ascending order.
c    On input, T describes a K subset.
c    On output, T describes the next K subset in the ordering.
c    If the input T was the last in the ordering, then the output T
c    will be the first.
c
c    Input/output, integer RANK, the rank.
c    If RANK = -1 on input, then the routine understands that this is
c    the first call, and that the user wishes the routine to supply
c    the first element in the ordering, which has RANK = 0.
c    In general, the input value of RANK is increased by 1 for output,
c    unless the very last element of the ordering was input, in which
c    case the output value of RANK is 0.
c
      implicit none

      integer k

      integer i
      integer ierror
      integer isave
      integer j
      integer n
      integer rank
      integer t(k)
c
c  Return the first element.
c
      if ( rank .eq. -1 ) then
        call i4vec_indicator ( k, t )
        rank = 0
        return
      end if
c
c  Check.
c
      call ksubset_lex_check ( k, n, t, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'KSUBSET_LEX_SUCCESSOR - Fatal error!'
        write ( *, '(a)' ) '  The input array is illegal.'
        write ( *, '(a,i8)' ) '  IERROR = ', ierror
        stop
      end if

      isave = 0

      do i = k, 1, -1
        if ( t(i) .ne. n - k + i ) then
          isave = i
          go to 10
        end if
      end do

10    continue
c
c  The last K subset was input.
c  Return the first one.
c
      if ( isave .eq. 0 ) then
        call i4vec_indicator ( k, t )
        rank = 0
      else

        do j = k, isave, -1
          t(j) = t(isave) + 1 + j - isave
        end do

        rank = rank + 1

      end if

      return
      end
      subroutine ksubset_lex_unrank ( rank, k, n, t )

c*********************************************************************72
c
cc KSUBSET_LEX_UNRANK computes the K subset of given lexicographic rank.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer RANK, the rank of the K subset.
c
c    Input, integer K, the number of elements each K subset must
c    have.  1 .le. K .le. N.
c
c    Input, integer N, the number of elements in the master set.
c    N must be positive.
c
c    Output, integer T(K), describes the K subset of the given
c    rank.  T(I) is the I-th element.  The elements must be listed in
c    ascending order.
c
      implicit none

      integer k

      integer binomial
      integer i
      integer n
      integer nksub
      integer rank
      integer rank_copy
      integer t(k)
      integer x
c
c  Check.
c
      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'KSUBSET_LEX_UNRANK - Fatal error!'
        write ( *, '(a)' ) '  Input N is illegal.'
        stop
      end if

      if ( k .lt. 1 .or. n .lt. k ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'KSUBSET_LEX_UNRANK - Fatal error!'
        write ( *, '(a)' ) '  Input K is illegal.'
        stop
      end if

      call ksubset_enum ( k, n, nksub )

      if ( rank .lt. 0 .or. nksub .lt. rank ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'KSUBSET_LEX_UNRANK - Fatal error!'
        write ( *, '(a)' ) '  Input rank is illegal.'
        stop
      end if

      rank_copy = rank

      x = 1

      do i = 1, k

10      continue

        if ( binomial ( n - x, k - i ) .le. rank_copy ) then
          rank_copy = rank_copy - binomial ( n - x, k - i )
          x = x + 1
          go to 10
        end if

        t(i) = x
        x = x + 1

      end do

      return
      end
      subroutine ksubset_revdoor_rank ( k, n, t, rank )

c*********************************************************************72
c
cc KSUBSET_REVDOOR_RANK computes the revolving door rank of a K subset.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer K, the number of elements each K subset must
c    have.  1 .le. K .le. N.
c
c    Input, integer N, the number of elements in the master set.
c    N must be positive.
c
c    Input, integer T(K), describes a K subset.  T(I) is the I-th
c    element.  The elements must be listed in ascending order.
c
c    Output, integer RANK, the rank of the K subset.
c
      implicit none

      integer k

      integer binomial
      integer i
      integer ierror
      integer n
      integer rank
      integer s
      integer t(k)
c
c  Check.
c
      call ksubset_lex_check ( k, n, t, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'KSUBSET_REVDOOR_RANK - Fatal error!'
        write ( *, '(a)' ) '  The input array is illegal.'
        write ( *, '(a,i8)' ) '  IERROR = ', ierror
        stop
      end if

      if ( mod ( k, 2 ) .eq. 0 ) then

        rank = 0

      else

        rank = - 1

      end if

      s = 1

      do i = k, 1, -1
        rank = rank + s * binomial ( t(i), i )
        s = - s
      end do

      return
      end
      subroutine ksubset_revdoor_successor ( k, n, t, rank )

c*********************************************************************72
c
cc KSUBSET_REVDOOR_SUCCESSOR computes the K subset revolving door successor.
c
c  Discussion:
c
c    After numerous attempts to implement the algorithm published in
c    Kreher and Stinson, the Nijenhuis and Wilf version was implemented
c    instead.  The K and S algorithm is supposedly based on the N and W one.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 March 2001
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer K, the number of elements each K subset must
c    have.  1 .le. K .le. N.
c
c    Input, integer N, the number of elements in the master set.
c    N must be positive.
c
c    Input/output, integer T(K), describes a K subset.  T(I) is the
c    I-th element.  The elements must be listed in ascending order.
c    On input, T describes a K subset.
c    On output, T describes the next K subset in the ordering.
c    If the input T was the last in the ordering, then the output T
c    will be the first.
c
c    Input/output, integer RANK, the rank.
c    If RANK = -1 on input, then the routine understands that this is
c    the first call, and that the user wishes the routine to supply
c    the first element in the ordering, which has RANK = 0.
c    In general, the input value of RANK is increased by 1 for output,
c    unless the very last element of the ordering was input, in which
c    case the output value of RANK is 0.
c
      implicit none

      integer k

      integer ierror
      integer j
      integer n
      integer rank
      integer t(k)
c
c  Return the first element.
c
      if ( rank .eq. - 1 ) then
        call i4vec_indicator ( k, t )
        rank = 0
        return
      end if
c
c  Check.
c
      call ksubset_lex_check ( k, n, t, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'KSUBSET_REVDOOR_SUCCESSOR - Fatal error!'
        write ( *, '(a)' ) '  The input array is illegal.'
        write ( *, '(a,i8)' ) '  IERROR = ', ierror
        stop
      end if

      j = 0

10    continue

        if ( 0 .lt. j .or. mod ( k, 2 ) .eq. 0 ) then

          j = j + 1

          if ( k .lt. j ) then
            t(k) = k
            rank = 0
            return
          end if

          if ( t(j) .ne. j ) then

            t(j) = t(j) - 1

            if ( j .ne. 1 ) then
              t(j-1) = j - 1
            end if

            rank = rank + 1
            return

          end if

        end if

        j = j + 1

        if ( j .lt. k ) then
          if ( t(j) .ne. t(j+1) - 1 ) then
            go to 20
          end if
        else
          if ( t(j) .ne. n ) then
            go to 20
          end if
        end if

      go to 10

20    continue

      t(j) = t(j) + 1

      if ( j .ne. 1 ) then
        t(j-1) = t(j) - 1
      end if

      rank = rank + 1

      return
      end
      subroutine ksubset_revdoor_unrank ( rank, k, n, t )

c*********************************************************************72
c
cc KSUBSET_REVDOOR_UNRANK computes the K subset of given revolving door rank.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer RANK, the rank of the K subset.
c
c    Input, integer K, the number of elements each K subset must
c    have.  1 .le. K .le. N.
c
c    Input, integer N, the number of elements in the master set.
c    N must be positive.
c
c    Output, integer T(K), describes the K subset of the given
c    rank.  T(I) is the I-th element.  The elements must be listed in
c    ascending order.
c
      implicit none

      integer k

      integer binomial
      integer i
      integer n
      integer nksub
      integer rank
      integer rank_copy
      integer t(k)
      integer x
c
c  Check.
c
      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'KSUBSET_REVDOOR_UNRANK - Fatal error!'
        write ( *, '(a)' ) '  Input N is illegal.'
        stop
      end if

      if ( k .lt. 1 .or. n .lt. k ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'KSUBSET_REVDOOR_UNRANK - Fatal error!'
        write ( *, '(a)' ) '  Input K is illegal.'
        stop
      end if

      call ksubset_enum ( k, n, nksub )

      if ( rank .lt. 0 .or. nksub .lt. rank ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'KSUBSET_REVDOOR_UNRANK - Fatal error!'
        write ( *, '(a)' ) '  The input rank is illegal.'
        stop
      end if

      rank_copy = rank

      x = n

      do i = k, 1, -1

10      continue

        if ( rank_copy .lt. binomial ( x, i ) ) then
          x = x - 1
          go to 10
        end if

        t(i) = x + 1
        rank_copy = binomial ( x + 1, i ) - rank_copy - 1

      end do

      return
      end
      subroutine marriage ( n, prefer, rank, fiancee, next )

c*********************************************************************72
c
cc MARRIAGE finds a stable set of marriages for given preferences.
c
c  Discussion:
c
c    Given a set of N men and N women who must be married in pairs,
c    and information defining the relative rankings that each person
c    assigns to the candidates of the opposite sex, this routine finds
c    a stable set of marriages for them.
c
c    A stable set of marriages is a pairing of the men and women with
c    the stability property: if M1 marries W1 and M2 marries W2, then
c    it is never the case that M1 and W2 would both prefer to be married
c    to each other.
c
c    An important application of stable marriage algorithms occurs in
c    the annual matching of medical residents to hospitals.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 February 2001
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Sedgewick,
c    Algorithms in C,
c    Addison-Wesley, 1990,
c    ISBN: 0-201-51425-7,
c    LC: QA76.73.C15S43.
c
c  Parameters:
c
c    Input, integer N, the number of pairs of men and women.
c
c    Input, integer PREFER(N,N); for man I, the value of
c    PREFER(I,J) represents his J-th preference for a wife.
c
c    Input, integer RANK(N,N); for woman I, the value of RANK(I,J)
c    represents her ranking of man number J.  A value of 1 for RANK(I,J)
c    means woman I ranks man J most preferable, while a value of N
c    would mean she ranked him least preferable.
c
c    Output, integer FIANCEE(N); for woman I, FIANCEE(I) is the
c    man to whom she is now engaged.
c
c    Output, integer NEXT(N); for man I, NEXT(I) is his preference
c    ranking for the woman to whom he is now engaged.  A value of 1 represents
c    his first choice, a value of N his last.
c
      implicit none

      integer n

      integer fiancee(n)
      integer i
      integer m
      integer next(n)
      integer prefer(n,n)
      integer rank(n,n)
      integer temp
      integer w
c
c  For man I, NEXT(I) is the woman I has most recently proposed to,
c  and hence NEXT(I)+1 is the next one to try.
c
      do i = 1, n
        next(i) = 0
      end do
c
c  For woman I, FIANCEE(I) is the man she has agree to marry,
c  or 0 if she has not agreed to any man yet.
c
      do i = 1, n
        fiancee(i) = 0
      end do
c
c  Start with an unengaged man, and end with an engaged woman.
c
      do i = 1, n

        m = i

10      continue

          next(m) = next(m) + 1

          w = prefer(m,next(m))

          if ( fiancee(w) .eq. 0 ) then
            fiancee(w) = m
            go to 20
          end if

          if ( rank (w,m) .lt. rank(w,fiancee(w)) ) then
            temp       = fiancee(w)
            fiancee(w) = m
            m          = temp
          end if

        go to 10

20      continue

      end do

      return
      end
      subroutine mountain ( n, x, y, mxy )

c*********************************************************************72
c
cc MOUNTAIN enumerates the mountains.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, ...
c    N must be positive.
c
c    Input, integer X, Y, ...
c    0 .le. X .le. 2 * N,
c
c    Output, integer MXY, the value of the "mountain function"
c    M ( N, X, Y ), which is the number of all mountain ranges from
c    (X,Y) to (2*N,0) which do not drop below sea level.
c
      implicit none

      integer a
      integer b
      integer binomial
      integer c
      integer mxy
      integer n
      integer x
      integer y
c
c  Check.
c
      if ( n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MOUNTAIN - Fatal error!'
        write ( *, '(a)' ) '  N .le. 0.'
        write ( *, '(a,i8)' ) '  N = ', n
        stop
      else if ( x .lt. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MOUNTAIN - Fatal error!'
        write ( *, '(a)' ) '  X .lt. 0.'
        write ( *, '(a,i8)' ) '  X = ', x
        stop
      else if ( 2 * n .lt. x ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MOUNTAIN - Fatal error!'
        write ( *, '(a)' ) '  2 * N .lt. X.'
        write ( *, '(a,i8)' ) '  X = ', x
        write ( *, '(a,i8)' ) '  N = ', n
        stop
      end if
c
c  Special cases.
c
      if ( y .lt. 0 ) then
        mxy = 0
        return
      end if

      if ( 2 * n .lt. x + y ) then
        mxy = 0
        return
      end if

      if ( mod ( x + y, 2 ) .eq. 1 ) then
        mxy = 0
        return
      end if

      a = 2 * n - x
      b = n - ( x + y ) / 2
      c = n - 1 - ( x + y ) / 2

      mxy = binomial ( a, b ) - binomial ( a, c )

      return
      end
      subroutine npart_enum ( n, npart, npartitions )

c*********************************************************************72
c
cc NPART_ENUM enumerates the number of partitions of N with NPART parts.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the integer to be partitioned.
c    Normally N must be positive, but for this routine any
c    N is allowed.
c
c    Input, integer NPART, the number of parts of the partition.
c    Normally, 1 .le. NPART .le. N is required,
c    but for this routine any value of NPART is allowed.
c
c    Output, integer NPARTITIONS is the number of partitions of N
c    with NPART parts.
c
      implicit none

      integer n

      integer npart
      integer npartitions
      integer p(0:n,0:n)

      if ( n .le. 0 ) then

        npartitions = 0

      else if ( npart .le. 0 .or. n .lt. npart ) then

        npartitions = 0

      else

        call npart_table ( n, npart, n, p )

        npartitions = p(n,npart)

      end if

      return
      end
      subroutine npart_rsf_lex_random ( n, npart, seed, a )

c*********************************************************************72
c
cc NPART_RSF_LEX_RANDOM returns a random RSF NPART partition.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 May 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the integer to be partitioned.
c    N must be positive.
c
c    Input, integer NPART, the number of parts of the partition.
c    1 .le. NPART .le. N.
c
c    Input/output, integer SEED, a seed for the random number
c    generator.
c
c    Output, integer A(NPART), contains the partition.
c    A(1) through A(NPART) contain the nonzero integers which
c    sum to N.
c
      implicit none

      integer n
      integer npart

      integer a(npart)
      integer i4_uniform
      integer npartitions
      integer rank
      integer seed

      call npart_enum ( n, npart, npartitions )

      rank = i4_uniform ( 1, npartitions, seed )

      call npart_rsf_lex_unrank ( rank, n, npart, a )

      return
      end
      subroutine npart_rsf_lex_rank ( n, npart, a, rank )

c*********************************************************************72
c
cc NPART_RSF_LEX_RANK computes the lex rank of an RSF NPART partition.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the integer to be partitioned.
c    N must be positive.
c
c    Input, integer NPART, the number of parts of the partition.
c    1 .le. NPART .le. N.
c
c    Input, integer A(NPART), contains the partition.
c    A(1) through A(NPART) contain the nonzero integers which
c    sum to N.
c
c    Output, integer RANK, the rank of the partition.
c
      implicit none

      integer n
      integer npart

      integer a(npart)
      integer b(npart)
      integer i
      integer ierror
      integer ncopy
      integer npartcopy
      integer p(0:n,0:npart)
      integer rank
c
c  Check.
c
      call part_rsf_check ( n, npart, a, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'NPART_RSF_LEX_RANK - Fatal error!'
        write ( *, '(a)' ) '  The input array is illegal.'
        write ( *, '(a,i8)' ) '  IERROR = ', ierror
        stop
      end if
c
c  Get the table of partitions of N with NPART parts.
c
      call npart_table ( n, npart, n, p )
c
c  Copy the partition "backwards".
c
      do i = 1, npart
        b(i) = a(npart+1-i)
      end do

      rank = 0
      ncopy = n
      npartcopy = npart

10    continue

      if ( 0 .lt. ncopy .and. 0 .lt. npartcopy ) then

        if ( b(npartcopy) .eq. 1 ) then

          ncopy = ncopy - 1
          npartcopy = npartcopy - 1

        else

          do i = 1, npartcopy
            b(i) = b(i) - 1
          end do
          rank = rank + p(ncopy-1,npartcopy-1)
          ncopy = ncopy - npartcopy

        end if

        go to 10

      end if

      return
      end
      subroutine npart_rsf_lex_successor ( n, npart, a, rank )

c*********************************************************************72
c
cc NPART_RSF_LEX_SUCCESSOR computes the RSF lex successor for NPART partitions.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 February 2001
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the integer to be partitioned.
c    N must be at least 1.
c
c    Input, integer NPART, the number of parts of the partition.
c    1 .le. NPART .le. N.
c
c    Input/output, integer A(NPART), contains the partition.
c    A(1) through A(NPART) contain the nonzero integers which
c    sum to N.
c
c    Input/output, integer RANK, the rank.
c    If RANK = -1 on input, then the routine understands that this is
c    the first call, and that the user wishes the routine to supply
c    the first element in the ordering, which has RANK = 0.
c    In general, the input value of RANK is increased by 1 for output,
c    unless the very last element of the ordering was input, in which
c    case the output value of RANK is 0.
c
      implicit none

      integer npart

      integer a(npart)
      integer d
      integer i
      integer ierror
      integer j
      integer n
      integer rank
c
c  Return the first element.
c
      if ( rank .eq. -1 ) then

        if ( npart .lt. 1 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'NPART_RSF_LEX_SUCCESSOR - Fatal error!'
          write ( *, '(a)' ) '  NPART .lt. 1.'
          stop
        end if

        do i = 1, npart - 1
          a(i) = 1
        end do
        a(npart) = n - ( npart - 1 )

        rank = 0
        return

      end if
c
c  Check.
c
      call part_rsf_check ( n, npart, a, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'NPART_RSF_LEX_SUCCESSOR - Fatal error!'
        write ( *, '(a)' ) '  The input array is illegal.'
        write ( *, '(a,i8)' ) '  IERROR = ', ierror
        stop
      end if
c
c  Find the first index I for which A(NPART+1-I) + 1 .lt. A(NPART).
c
      i = 2

10    continue

        if ( npart .lt. i ) then
          go to 20
        end if

        if ( a(npart+1-i) + 1 .lt. a(npart) ) then
          go to 20
        end if

        i = i + 1

      go to 10

20    continue
c
c  If no such index, we've reached the end of the line.
c
      if ( i .eq. npart + 1 ) then

        do j = 1, npart - 1
          a(j) = 1
        end do
        a(npart) = n - ( npart - 1 )

        rank = 0
        return
c
c  Otherwise, increment A(NPART+1-I), and adjust other entries.
c
      else

        a(npart+1-i) = a(npart+1-i) + 1
        d = - 1

        do j = i - 1, 2, -1
          d = d + a(npart+1-j) - a(npart+1-i)
          a(npart+1-j) = a(npart+1-i)
        end do

        a(npart) = a(npart) + d

      end if

      rank = rank + 1

      return
      end
      subroutine npart_rsf_lex_unrank ( rank, n, npart, a )

c*********************************************************************72
c
cc NPART_RSF_LEX_UNRANK unranks an RSF NPART partition in the lex ordering.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 April 2001
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer RANK, the rank of the partition.
c
c    Input, integer N, the integer to be partitioned.
c    N must be positive.
c
c    Input, integer NPART, the number of parts of the partition.
c    1 .le. NPART .le. N.
c
c    Output, integer A(NPART), contains the partition.
c    A(1) through A(NPART) contain the nonzero integers which
c    sum to N.
c
      implicit none

      integer n
      integer npart

      integer a(npart)
      integer i
      integer ncopy
      integer npartcopy
      integer npartitions
      integer p(0:n,0:npart)
      integer rank
      integer rank_copy
c
c  Check.
c
      if ( n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'NPART_RSF_LEX_UNRANK - Fatal error!'
        write ( *, '(a)' ) '  The input N is illegal.'
        stop
      end if

      if ( npart .lt. 1 .or. n .lt. npart ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'NPART_RSF_LEX_UNRANK - Fatal error!'
        write ( *, '(a)' ) '  The input NPART is illegal.'
        stop
      end if

      call npart_enum ( n, npart, npartitions )

      if ( rank .lt. 0 .or. npartitions .lt. rank ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'NPART_RSF_LEX_UNRANK - Fatal error!'
        write ( *, '(a)' ) '  The input rank is illegal.'
        stop
      end if
c
c  Get the table of partitions of N with NPART parts.
c
      call npart_table ( n, npart, n, p )

      do i = 1, npart
        a(i) = 0
      end do

      rank_copy = rank
      ncopy = n
      npartcopy = npart

10    continue

      if ( 0 .lt. ncopy ) then

        if ( rank_copy .lt. p(ncopy-1,npartcopy-1) ) then
          a(npart+1-npartcopy) = a(npart+1-npartcopy) + 1
          ncopy = ncopy - 1
          npartcopy = npartcopy - 1
        else
          do i = 1, npartcopy
            a(npart+1-i) = a(npart+1-i) + 1
          end do
          rank_copy = rank_copy - p(ncopy-1,npartcopy-1)
          ncopy = ncopy - npartcopy
        end if

        go to 10

      end if

      return
      end
      subroutine npart_sf_lex_successor ( n, npart, a, rank )

c*********************************************************************72
c
cc NPART_SF_LEX_SUCCESSOR computes SF NPART partition.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 February 2001
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the integer to be partitioned.
c    N must be positive.
c
c    Input, integer NPART, the number of parts of the partition.
c    1 .le. NPART .le. N.
c
c    Input/output, integer A(NPART), contains the partition.
c    A(1) through A(NPART) contain the nonzero integers which
c    sum to N.  The values in A must be in DESCENDING order.
c
c    Input/output, integer RANK, the rank.
c    If RANK = -1 on input, then the routine understands that this is
c    the first call, and that the user wishes the routine to supply
c    the first element in the ordering, which has RANK = 0.
c    In general, the input value of RANK is increased by 1 for output,
c    unless the very last element of the ordering was input, in which
c    case the output value of RANK is 0.
c
      implicit none

      integer npart

      integer a(npart)
      integer i
      integer ierror
      integer indx
      integer n
      integer rank
      integer temp
c
c  Return the first element.
c
      if ( rank .eq. -1 ) then
        call i4vec_part2 ( n, npart, a )
        rank = 0
        return
      end if
c
c  Check.
c
      call part_sf_check ( n, npart, a, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'NPART_SF_LEX_SUCCESSOR - Fatal error!'
        write ( *, '(a)' ) '  The input array is illegal.'
        write ( *, '(a,i8)' ) '  IERROR = ', ierror
        stop
      end if
c
c  Find the last entry that is 2 or more.
c
      do i = npart, 1, - 1
        if ( 1 .lt. a(i) ) then
          indx = i
          go to 10
        end if
      end do

10    continue
c
c  As long as the last nonunit occurs after the first position,
c  have it donate 1 to the left.
c
      if ( 1 .lt. indx ) then

        a(indx) = a(indx) - 1
        a(indx-1) = a(indx-1) + 1
        indx = indx - 1

20      continue

          if ( indx .le. 1 ) then
            go to 30
          end if

          if ( a(indx) .le. a(indx-1) ) then
            go to 30
          end if

          temp      = a(indx)
          a(indx)   = a(indx-1)
          a(indx-1) = temp

          indx = indx - 1

        go to 20

30      continue
c
c  Sum the tail.
c
        temp = 0
        do i = indx + 1, npart
          temp = temp + a(i)
        end do
c
c  Partition the tail sum equally over the tail.
c
        call i4vec_part2 ( temp, npart - indx, a(indx+1) )

        rank = rank + 1
c
c  If A(2) through A(NPART) are 1, then this is the last element.
c  Return the first one.
c
      else

        call i4vec_part2 ( n, npart, a )
        rank = 0

      end if

      return
      end
      subroutine npart_table ( n, npart, nmax, p )

c*********************************************************************72
c
cc NPART_TABLE tabulates the number of partitions of N having NPART parts.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the integer to be partitioned.
c    N must be positive.
c
c    Input, integer NPART, the number of parts of the partition.
c    1 .le. NPART .le. N.
c
c    Input, integer NMAX, the leading dimension of P.
c
c    Output, integer P(0:NMAX,0:NPART), P(I,J) is the number of
c    partitions of I having J parts.
c
      implicit none

      integer nmax
      integer npart

      integer i
      integer j
      integer n
      integer p(0:nmax,0:npart)

      p(0,0) = 1
      p(1:n,0) = 0

      do i = 1, n
        do j = 1, npart
          if ( i .lt. j ) then
            p(i,j) = 0
          else if ( i .lt. 2 * j ) then
            p(i,j) = p(i-1,j-1)
          else
            p(i,j) = p(i-1,j-1) + p(i-j,j)
          end if
        end do
      end do

      return
      end
      subroutine part_enum ( n, npartitions )

c*********************************************************************72
c
cc PART_ENUM enumerates the number of partitions of N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the integer to be partitioned.
c    Normally N must be positive, but for this routine any
c    N is allowed.
c
c    Output, integer NPARTITIONS is the number of partitions of N.
c
      implicit none

      integer n

      integer npartitions
      integer p(0:n)

      if ( n .lt. 0 ) then

        npartitions = 0

      else

        call part_table ( n, p )

        npartitions = p(n)

      end if

      return
      end
      subroutine part_rsf_check ( n, npart, a, ierror )

c*********************************************************************72
c
cc PART_RSF_CHECK checks a reverse standard form partition of an integer.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the integer to be partitioned.
c    N must be positive.
c
c    Input, integer NPART, the number of parts of the partition.
c    1 .le. NPART .le. N.
c
c    Input, integer A(NPART), contains the partition.
c    A(1) through A(NPART) contain the nonzero integers which
c    sum to N.  The entries must be in ASCENDING order.
c
c    Output, integer IERROR, error flag.
c    0, no error.
c    -1, N is illegal.
c    -2, NPART is illegal.
c    -3, the entries do not add up to N.
c    I, the I-th entry of A is illegal.
c
      implicit none

      integer npart

      integer a(npart)
      integer i
      integer i4vec_sum
      integer ierror
      integer n

      ierror = 0

      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PART_RSF_CHECK - Fatal error!'
        write ( *, '(a)' ) '  N .lt. 1.'
        stop
      end if

      if ( npart .lt. 1 .or. n .lt. npart ) then
        ierror = -2
        return
      end if
c
c  Every entry must lie between 1 and N.
c
      do i = 1, npart
        if ( a(i) .lt. 1 .or. n .lt. a(i) ) then
          ierror = i
          return
        end if
      end do
c
c  The entries must be in ascending order.
c
      do i = 2, npart
        if ( a(i) .lt. a(i-1) ) then
          ierror = i
          return
        end if
      end do
c
c  The entries must add up to N.
c
      if ( i4vec_sum ( npart, a ) .ne. n ) then
        ierror = -3
      end if

      return
      end
      subroutine part_sf_check ( n, npart, a, ierror )

c*********************************************************************72
c
cc PART_SF_CHECK checks a standard form partition of an integer.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 February 2001
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the integer to be partitioned.
c    N must be positive.
c
c    Input, integer NPART, the number of parts of the partition.
c    1 .le. NPART .le. N.
c
c    Input, integer A(NPART), contains the partition.
c    A(1) through A(NPART) contain the nonzero integers which
c    sum to N.  The entries must be in DESCENDING order.
c
c    Output, integer IERROR, error flag.
c    0, no error.
c    -1, N is illegal.
c    -2, NPART is illegal.
c    -3, the entries do not add up to N.
c    I, the I-th entry of A is illegal.
c
      implicit none

      integer npart

      integer a(npart)
      integer i
      integer i4vec_sum
      integer ierror
      integer n

      ierror = 0

      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PART_SF_CHECK - Fatal error!'
        write ( *, '(a)' ) '  N .lt. 1.'
        stop
      end if

      if ( npart .lt. 1 .or. n .lt. npart ) then
        ierror = -2
        return
      end if
c
c  Every entry must lie between 1 and N.
c
      do i = 1, npart
        if ( a(i) .lt. 1 .or. n .lt. a(i) ) then
          ierror = i
          return
        end if
      end do
c
c  The entries must be in descending order.
c
      do i = 2, npart
        if ( a(i-1) .lt. a(i) ) then
          ierror = i
          return
        end if
      end do
c
c  The entries must add up to N.
c
      if ( i4vec_sum ( npart, a ) .ne. n ) then
        ierror = -3
      end if

      return
      end
      subroutine part_sf_conjugate ( n, npart, a, npart2, b )

c*********************************************************************72
c
cc PART_SF_CONJUGATE computes the conjugate of a partition.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the integer to be partitioned.
c    N must be positive.
c
c    Input, integer NPART, the number of parts of the partition.
c    1 .le. NPART .le. N.
c
c    Input, integer A(N), contains the partition.
c    A(1) through A(NPART) contain the nonzero integers which
c    sum to N.
c
c    Output, integer NPART2, the number of parts of the conjugate
c    partition.
c
c    Output, integer B(N), contains the conjugate partition.
c
      implicit none

      integer n

      integer a(n)
      integer b(n)
      integer i
      integer ierror
      integer j
      integer npart
      integer npart2
c
c  Check.
c
      call part_sf_check ( n, npart, a, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PART_SF_CONJUGATE - Fatal error!'
        write ( *, '(a)' ) '  The input array is illegal.'
        write ( *, '(a,i8)' ) '  IERROR = ', ierror
        stop
      end if

      npart2 = a(1)
      do i = 1, npart2
        b(i) = 0
      end do

      do i = 1, npart
        do j = 1, a(i)
          b(j) = b(j) + 1
        end do
      end do

      return
      end
      subroutine part_sf_majorize ( n, nparta, a, npartb, b, result )

c*********************************************************************72
c
cc PART_SF_MAJORIZE determines if partition A majorizes partition B.
c
c  Discussion:
c
c    The partitions must be in standard form.
c
c    If A, with NPARTA parts, and B, with NPARTB parts, are both partitions
c    of the same positive integer N, then we say that A majorizes B if,
c    for every index K from 1 to N, it is true that
c
c      sum ( 1 .le. I .le. K ) B(I) .le. sum ( 1 .le. I .le. K ) A(I)
c
c    where entries of A beyond index NPARTA, and of B beyond BPARTB
c    are assumed to be 0.  We say that A strictly majorizes B if
c    A majorizes B, and for at least one index K the inequality is strict.
c
c    For any two partitions of N, it is possible that A majorizes B,
c    B majorizes A, both partitions majorize each other (in which case
c    they are equal), or that neither majorizes the other.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Jack vanLint, Richard Wilson,
c    A Course in Combinatorics,
c    Cambridge, 1992,
c    ISBN: 0-521-42260-4,
c    LC: QA164.L56.
c
c  Parameters:
c
c    Input, integer N, the integer to be partitioned.
c    N must be positive.
c
c    Input, integer NPARTA, the number of parts in partition A.
c    1 .le. NPARTA .le. N.
c
c    Input, integer A(NPARTA), contains partition A in standard
c    form.  A(1) through A(NPARTA) contain nonzero integers which sum to N.
c
c    Input, integer NPARTB, the number of parts in partition B.
c    1 .le. NPARTB .le. N.
c
c    Input, integer B(NPARTB), contains partition B in standard
c    form.  B(1) through B(NPARTB) contain nonzero integers which sum to N.
c
c    Output, integer RESULT, the result of the comparison.
c    -2, A and B are incomparable, but would have been +1.
c    -1, A .lt. B, (A is strictly majorized by B),
c     0, A = B, (A and B are identical),
c    +1, A > B, (A strictly majorizes B),
c    +2, A and B are incomparable, but would have been +1.
c
      implicit none

      integer nparta
      integer npartb

      integer a(nparta)
      integer b(npartb)
      integer i
      integer ierror
      integer n
      integer result
      integer suma
      integer sumb
c
c  Check.
c
      call part_sf_check ( n, nparta, a, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PART_SF_MAJORIZE - Fatal error!'
        write ( *, '(a)' ) '  The input array A is illegal.'
        stop
      end if

      call part_sf_check ( n, npartb, b, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PART_SF_MAJORIZE - Fatal error!'
        write ( *, '(a)' ) '  The input array B is illegal.'
        stop
      end if

      result = 0
      suma = 0
      sumb = 0

      do i = 1, min ( nparta, npartb )

        if ( i .le. nparta ) then
          suma = suma + a(i)
        end if

        if ( i .le. npartb ) then
          sumb = sumb + b(i)
        end if

        if ( result .eq. -1 ) then

          if ( sumb .lt. suma ) then
            result = -2
            return
          end if

        else if ( result .eq. 0 ) then

          if ( suma .lt. sumb ) then
            result = -1
          else if ( sumb .lt. suma ) then
            result = +1
          end if

        else if ( result .eq. + 1 ) then

          if ( suma .lt. sumb ) then
            result = +2
            return
          end if

        end if

      end do

      return
      end
      subroutine part_successor ( n, npart, a, rank )

c*********************************************************************72
c
cc PART_SUCCESSOR computes the lexicographic partition successor.
c
c  Discussion:
c
c    PART_SUCCESSOR is "inspired by" the GenPartitions algorithm,
c    but instead of relying on recursion, generates the partitions
c    one at a time.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the integer to be partitioned.
c    N must be positive.
c
c    Input/output, integer NPART, the number of parts of the
c    partition.  1 .le. NPART .le. N.
c
c    Input/output, integer A(N), contains the partition.
c    A(1) through A(NPART) contain the nonzero integers which
c    sum to N.
c
c    Input/output, integer RANK, the rank.
c    If RANK = -1 on input, then the routine understands that this is
c    the first call, and that the user wishes the routine to supply
c    the first element in the ordering, which has RANK = 0.
c    In general, the input value of RANK is increased by 1 for output,
c    unless the very last element of the ordering was input, in which
c    case the output value of RANK is 0.
c
      implicit none

      integer n

      integer a(n)
      integer asum
      integer i
      integer ierror
      integer ihi
      integer j
      integer npart
      integer rank
c
c  Return the first element.
c
      if ( rank .eq. -1 ) then
        do i = 1, n
          a(i) = 1
        end do
        npart = n
        rank = 0
        return
      end if
c
c  Check.
c
      call part_sf_check ( n, npart, a, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PART_SUCCESSOR - Fatal error!'
        write ( *, '(a)' ) '  The input array is illegal.'
        stop
      end if
c
c  If possible, increment the first intermediate position that
c  is less than its left hand neighbor, and has at least one
c  right hand neighbor.
c
      ihi = npart - 1

      do i = ihi, 2, -1

        if ( a(i) .lt. a(i-1) ) then
          asum = - 1
          do j = i + 1, npart
            asum = asum + a(j)
          end do
          a(i) = a(i) + 1
          do j = i + 1, npart
            a(j) = 0
          end do
          npart = i + asum
          do j = i + 1, npart
            a(j) = 1
          end do
          rank = rank + 1
          return
        end if

      end do
c
c  A) there are two or more parts
c  Increment the first, replace the rest by 1's.
c
      if ( 2 .le. npart ) then
        a(1) = a(1) + 1
        do i = 2, npart
          a(i) = 0
        end do
        npart = n - a(1) + 1
        do i = 2, npart
          a(i) = 1
        end do
        rank = rank + 1
c
c  B) there's only one part.
c  We've reached the last item.
c  Return the first one.
c
      else if ( npart .eq. 1 ) then
        do i = 1, n
          a(i) = 1
        end do
        npart = n
        rank = 0
      end if

      return
      end
      subroutine part_table ( n, p )

c*********************************************************************72
c
cc PART_TABLE tabulates the number of partitions of N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the integer to be partitioned.
c    N must be positive.
c
c    Output, integer P(0:N), P(I) is the number of partitions of I.
c
      implicit none

      integer n

      integer i
      integer j
      integer p(0:n)
      integer psum
      integer sign
      integer w
      integer wprime

      p(0) = 1
      p(1) = 1

      do i = 2, n

        sign = 1
        psum = 0
        w = 1
        j = 1
        wprime = w + j

10      continue

        if ( w .lt. n ) then

          if ( 0 .le. i - w ) then
            if ( sign .eq. 1 ) then
              psum = psum + p(i-w)
            else
              psum = psum - p(i-w)
            end if
          end if

          if ( wprime .le. i ) then

            if ( sign .eq. 1 ) then
              psum = psum + p(i-wprime)
            else
              psum = psum - p(i-wprime)
            end if

          end if

          w = w + 3 * j + 1
          j = j + 1
          wprime = w + j
          sign = - sign

          go to 10

        end if

        p(i) = psum

      end do

      return
      end
      subroutine partition_greedy ( n, a, indx )

c*********************************************************************72
c
cc PARTITION_GREEDY attacks the partition problem with a greedy algorithm.
c
c  Discussion:
c
c    Given a collection of N not necessarily distinct positive integers A(I),
c    it is desired to partition the values into two groups, whose sums are
c    as close as possible.
c
c  Algorithm:
c
c    Begin with sets 1 and 2 empty.
c
c    Process the data in descending order of magnitude.
c
c    The next item A(I) is added to set 1 or set 2, whichever has the
c    smallest current sum.
c
c    Stop as soon as all items have been allocated.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Brian Hayes,
c    The Easiest Hard Problem,
c    American Scientist,
c    Volume 90, Number 2, March-April 2002, pages 113-117.
c
c  Parameters:
c
c    Input, integer N, the number of values.  N must be positive.
c
c    Input/output, integer A(N), a collection of positive values.
c    On output, A has been sorted into descending order.
c
c    Output, integer INDX(N); INDX(I) is 1 if A(I) is part of
c    set 1, and 2 if it is assigned to set 2.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer indx(n)
      integer j
      integer sums(2)

      sums(1) = 0
      sums(2) = 0

      call i4vec_sort_insert_d ( n, a )

      do i = 1, n

        if ( sums(1) .lt. sums(2) ) then
          j = 1
        else
          j = 2
        end if

        indx(i) = j
        sums(j) = sums(j) + a(i)

      end do

      return
      end
      subroutine partn_enum ( n, nmax, npartitions )

c*********************************************************************72
c
cc PARTN_ENUM enumerates the partitions of N with maximum element NMAX.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the integer to be partitioned.
c    Normally N must be positive, but for this routine any
c    N is allowed.
c
c    Input, integer NMAX, the maximum element in the partition.
c    Normally, 1 .le. NMAX .le. N is required,
c    but for this routine any value of NMAX is allowed.
c
c    Output, integer NPARTITIONS is the number of partitions of N
c    with maximum element NMAX.
c
      implicit none

      integer nbig
      parameter ( nbig = 25 )

      integer n
      integer nmax
      integer npartitions
      integer p(0:nbig,0:nbig)

      if ( n .le. 0 ) then

        npartitions = 0

      else if ( nmax .le. 0 .or. n .lt. nmax ) then

        npartitions = 0

      else

        call npart_table ( n, nmax, nbig, p )

        npartitions = p(n,nmax)

      end if

      return
      end
      subroutine partn_sf_check ( n, nmax, npart, a, ierror )

c*********************************************************************72
c
cc PARTN_SF_CHECK checks an SF partition of an integer with largest entry NMAX.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the integer to be partitioned.
c    N must be positive.
c
c    Input, integer NMAX, the value of the largest entry.
c    1 .le. NMAX .le. N.
c
c    Input, integer NPART, the number of parts of the partition.
c    1 .le. NPART .le. N.
c
c    Input, integer A(NPART), contains the partition.
c    A(1) through A(NPART) contain the nonzero integers which
c    sum to N.  The entries must be in DESCENDING order.
c
c    Output, integer IERROR, error flag.
c    0, no error.
c    -1, N is illegal.
c    -2, NMAX is illegal.
c    -3, NPART is illegal.
c    -3, the entries do not add up to N.
c    I, the I-th entry of A is illegal.
c
      implicit none

      integer npart

      integer a(npart)
      integer asum
      integer i
      integer ierror
      integer n
      integer nmax

      ierror = 0

      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PARTN_SF_CHECK - Fatal error!'
        write ( *, '(a)' ) '  N .lt. 1.'
        stop
      end if

      if ( nmax .lt. 1 .or. n .lt. nmax ) then
        ierror = -2
        return
      end if

      if ( npart .lt. 1 .or. n .lt. npart ) then
        ierror = -3
        return
      end if
c
c  Entry 1 must be NMAX.
c
      if ( a(1) .ne. nmax ) then
        ierror = 1
        return
      end if
c
c  Every entry must lie between 1 and N.
c
      do i = 1, npart
        if ( a(i) .lt. 1 .or. n .lt. a(i) ) then
          ierror = i
          return
        end if
      end do
c
c  The entries must be in descending order.
c
      do i = 2, npart
        if ( a(i-1) .lt. a(i) ) then
          ierror = i
          return
        end if
      end do
c
c  The entries must add up to N.
c
      asum = 0
      do i = 1, npart
        asum = asum + a(i)
        if ( n .lt. asum ) then
          ierror = i
          return
        end if
      end do

      if ( asum .ne. n ) then
        ierror = -3
      end if

      return
      end
      subroutine partn_successor ( n, nmax, npart, a, rank )

c*********************************************************************72
c
cc PARTN_SUCCESSOR computes partitions whose largest part is NMAX.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the integer to be partitioned.
c    N must be positive.
c
c    Input, integer NMAX, the maximum size of any part of the
c    partition.  1 .le. NMAX .le. N.
c
c    Input/output, integer NPART, the number of parts of the
c    partition.  1 .le. NPART .le. N.
c
c    Input/output, integer A(N), contains the partition.
c    A(1) through A(NPART) contain the nonzero integers which
c    sum to N.
c
c    Input/output, integer RANK, the rank.
c    If RANK = -1 on input, then the routine understands that this is
c    the first call, and that the user wishes the routine to supply
c    the first element in the ordering, which has RANK = 0.
c    In general, the input value of RANK is increased by 1 for output,
c    unless the very last element of the ordering was input, in which
c    case the output value of RANK is 0.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer ierror
      integer index
      integer nmax
      integer npart
      integer rank
      integer temp
c
c  Return the first element.
c
      if ( rank .eq. -1 ) then
        a(1) = nmax
        npart = n + 1 - nmax
        do i = 2, npart
          a(i) = 1
        end do
        rank = 0
        return
      end if
c
c  Check.
c
      call partn_sf_check ( n, nmax, npart, a, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PARTN_SUCCESSOR - Fatal error!'
        write ( *, '(a)' ) '  The input array is illegal.'
        write ( *, '(a,i8)' ) '  IERROR = ', ierror
        stop
      end if
c
c  If there are at least two parts, and the next to last is not NMAX,
c  then rob the last part and pay the next to the last part.
c  Then, if the next to last part is too big, swap it leftwards.
c
      if ( 1 .lt. npart ) then

        if ( a(npart-1) .lt. nmax ) then

          a(npart) = a(npart) - 1
          a(npart-1) = a(npart-1) + 1
          index = npart - 1

10        continue

            if ( index .le. 1 ) then
              go to 20
            end if

            if ( a(index) .le. a(index-1) ) then
              go to 20
            end if

            temp       = a(index-1)
            a(index-1) = a(index)
            a(index)   = temp

            index = index - 1

          go to 10

20        continue
c
c  Sum the tail.
c
          temp = 0
          do i = index + 1, npart
            temp = temp + a(i)
          end do
c
c  Spread the sum as 1's.
c
          npart = index + temp
          do i = index + 1, npart
            a(i) = 1
          end do
          rank = rank + 1
          return

        end if
c
c  Otherwise, we've reached the last item.
c  Return the first one.
c
      else

        npart = n + 1 - nmax
        a(1) = nmax
        do i = 2, npart
          a(i) = 1
        end do
        rank = 0
        return

      end if

      return
      end
      subroutine perm_check ( n, p )

c*********************************************************************72
c
cc PERM_CHECK checks a representation of a permutation.
c
c  Discussion:
c
c    The routine is given N and P, a vector of length N.
c    P is a legal represention of a permutation of the integers from
c    1 to N if and only if every integer from 1 to N occurs
c    as a value of P(I) for some I between 1 and N.
c
c    If an error is observed in the input, this routine prints a message
c    and stops.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 July 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of values being permuted.
c    N must be positive.
c
c    Input, integer P(N), the array to check.
c
      implicit none

      integer n

      integer i
      integer ifind
      integer iseek
      integer p(n)

      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PERM_CHECK - Fatal error!'
        write ( *, '(a)' ) '  N .lt. 1.'
        write ( *, '(a,i12)' ) '  N = ', n
        stop
      end if

      do i = 1, n
        if ( p(i) .lt. 1 .or. n .lt. p(i) ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'PERM_CHECK - Fatal error!'
          write ( *, '(a)' ) '  P(I) .lt. 1 or N .lt. P(I).'
          write ( *, '(a,i12)' ) '  N = ', n
          write ( *, '(a,i12)' ) '  I = ', i
          write ( *, '(a,i12)' ) '  P(I) = ', p(i)
          stop
        end if
      end do

      do iseek = 1, n

        ifind = -1

        do i = 1, n
          if ( p(i) .eq. iseek ) then
            ifind = i
            go to 10
          end if
        end do

10      continue

        if ( ifind .eq. -1 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'PERM_CHECK - Fatal error!'
          write ( *, '(a)' ) '  Every I from 1 to N must occur in P.'
          write ( *, '(a)' ) '  Could not find occurrence of ', iseek
          stop
        end if

      end do

      return
      end
      subroutine perm_enum ( n, nperm )

c*********************************************************************72
c
cc PERM_ENUM enumerates the permutations on N digits.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of values being permuted.
c    N must be nonnegative.
c
c    Output, integer NPERM, the number of distinct elements.
c
      implicit none

      integer i4_factorial
      integer n
      integer nperm

      nperm = i4_factorial ( n )

      return
      end
      subroutine perm_inv ( n, p, pinv  )

c*********************************************************************72
c
cc PERM_INV computes the inverse of a permutation.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of values being permuted.
c    N must be positive.
c
c    Input, integer P(N), describes the permutation.
c    P(I) is the item which is permuted into the I-th place
c    by the permutation.
c
c    Output, integer PINV(N), the inverse permutation.
c
      implicit none

      integer n

      integer i
      integer p(n)
      integer pinv(n)
c
c  Check.
c
      call perm_check ( n, p )

      do i = 1, n
        pinv(p(i)) = i
      end do

      return
      end
      subroutine perm_lex_rank ( n, p, rank )

c*********************************************************************72
c
cc PERM_LEX_RANK computes the lexicographic rank of a permutation.
c
c  Discussion:
c
c    The original code altered the input permutation.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of values being permuted.
c    N must be positive.
c
c    Input, integer P(N), describes the permutation.
c    P(I) is the item which is permuted into the I-th place
c    by the permutation.
c
c    Output, integer RANK, the rank of the permutation.
c
      implicit none

      integer n

      integer i
      integer i4_factorial
      integer j
      integer p(n)
      integer pcopy(n)
      integer rank
c
c  Check.
c
      call perm_check ( n, p )

      rank = 0
      do i = 1, n
        pcopy(i) = p(i)
      end do

      do j = 1, n

        rank = rank + ( pcopy(j) - 1 ) * i4_factorial ( n - j )

        do i = j + 1, n
          if ( pcopy(j) .lt. pcopy(i) ) then
            pcopy(i) = pcopy(i) - 1
          end if
        end do

      end do

      return
      end
      subroutine perm_lex_successor ( n, p, rank )

c*********************************************************************72
c
cc PERM_LEX_SUCCESSOR computes the lexicographic permutation successor.
c
c  Example:
c
c    RANK  Permutation
c
c       0  1 2 3 4
c       1  1 2 4 3
c       2  1 3 2 4
c       3  1 3 4 2
c       4  1 4 2 3
c       5  1 4 3 2
c       6  2 1 3 4
c       ...
c      23  4 3 2 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of values being permuted.
c    N must be positive.
c
c    Input/output, integer P(N), describes the permutation.
c    P(I) is the item which is permuted into the I-th place
c    by the permutation.
c
c    Input/output, integer RANK, the rank.
c    If RANK = -1 on input, then the routine understands that this is
c    the first call, and that the user wishes the routine to supply
c    the first element in the ordering, which has RANK = 0.
c    In general, the input value of RANK is increased by 1 for output,
c    unless the very last element of the ordering was input, in which
c    case the output value of RANK is 0.
c
      implicit none

      integer n

      integer i
      integer j
      integer p(n)
      integer rank
      integer temp
c
c  Return the first element.
c
      if ( rank .eq. -1 ) then
        call i4vec_indicator ( n, p )
        rank = 0
        return
      end if
c
c  Check.
c
      call perm_check ( n, p )
c
c  Seek I, the highest index for which the next element is bigger.
c
      i = n - 1

10    continue

        if ( i .le. 0 ) then
          go to 20
        end if

        if ( p(i) .le. p(i+1) ) then
          go to 20
        end if

        i = i - 1

      go to 10

20    continue
c
c  If no I could be found, then we have reach the final permutation,
c  N, N-1, ..., 2, 1.  Time to start over again.
c
      if ( i .eq. 0 ) then
        call i4vec_indicator ( n, p )
        rank = 0
      else
c
c  Otherwise, look for the the highest index after I whose element
c  is bigger than I's.  We know that I+1 is one such value, so the
c  loop will never fail.
c
        j = n

30      continue

        if ( p(j) .lt. p(i) ) then
          j = j - 1
          go to 30
        end if
c
c  Interchange elements I and J.
c
        temp = p(i)
        p(i) = p(j)
        p(j) = temp
c
c  Reverse the elements from I+1 to N.
c
        call i4vec_reverse ( n - i, p(i+1:n) )
        rank = rank + 1

      end if

      return
      end
      subroutine perm_lex_unrank ( rank, n, p )

c*********************************************************************72
c
cc PERM_LEX_UNRANK computes the permutation of given lexicographic rank.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer RANK, the rank of the permutation.
c
c    Input, integer N, the number of values being permuted.
c    N must be positive.
c
c    Output, integer P(N), describes the permutation.
c
      implicit none

      integer n

      integer d
      integer i
      integer i4_factorial
      integer j
      integer nperm
      integer p(n)
      integer rank
      integer rank_copy
c
c  Check.
c
      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PERM_LEX_UNRANK - Fatal error!'
        write ( *, '(a)' ) '  Input N is illegal.'
        stop
      end if

      call perm_enum ( n, nperm )

      if ( rank .lt. 0 .or. nperm .lt. rank ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PERM_LEX_UNRANK - Fatal error!'
        write ( *, '(a)' ) '  The input rank is illegal.'
        stop
      end if

      rank_copy = rank

      p(n) = 1

      do j = 1, n - 1

        d = mod ( rank_copy, i4_factorial ( j + 1 ) ) 
     &    / i4_factorial ( j )
        rank_copy = rank_copy - d * i4_factorial ( j )
        p(n-j) = d + 1

        do i = n - j + 1, n

          if ( d .lt. p(i) ) then
            p(i) = p(i) + 1
          end if

        end do

      end do

      return
      end
      subroutine perm_mul ( n, p, q, r  )

c*********************************************************************72
c
cc PERM_MUL computes the product of two permutations.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,inson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of values being permuted.
c    N must be positive.
c
c    Input, integer P(N), Q(N), describes the permutation factors.
c
c    Output, integer R(N), the product permutation P * Q.
c    R(I) = P(Q(I)).
c
      implicit none

      integer n

      integer i
      integer p(n)
      integer q(n)
      integer r(n)
      integer s(n)
c
c  Check.
c
      call perm_check ( n, p )

      call perm_check ( n, q )
c
c  Use a temporary vector for the result, to avoid problems if
c  some arguments are actually identified.
c
      do i = 1, n
        s(i) = p(q(i))
      end do

      do i = 1, n
        r(i) = s(i)
      end do

      return
      end
      subroutine perm_parity ( n, p, parity )

c*********************************************************************72
c
cc PERM_PARITY computes the parity of a permutation.
c
c  Discussion:
c
c    The routine requires the use of a temporary array.
c
c    A permutation is called "even" or "odd", depending on whether
c    it is equivalent to an even or odd number of pairwise
c    transpositions.  This is known as the "parity" of the
c    permutation.
c
c    The "sign" of a permutation is +1 if it has even parity,
c    and -1 if it has odd parity.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of values being permuted.
c    N must be positive.
c
c    Input, integer P(N), describes the permutation.
c    P(I) is the item which is permuted into the I-th place
c    by the permutation.
c
c    Output, integer PARITY, the parity of the permutation.
c    0, the permutation has even parity.
c    1, the permutation has odd parity.
c
      implicit none

      integer n

      integer a(n)
      integer c
      integer i
      integer j
      integer parity
      integer p(n)
c
c  Check.
c
      call perm_check ( n, p )

      do i = 1, n
        a(i) = 0
      end do

      c = 0

      do j = 1, n

        if ( a(j) .eq. 0 ) then

          c = c + 1
          a(j) = 1
          i = j

10        continue

          if ( p(i) .ne. j ) then
            i = p(i)
            a(i) = 1
            go to 10
          end if

        end if

      end do

      parity = mod ( n - c, 2 )

      return
      end
      subroutine perm_print ( n, p, title )

c*********************************************************************72
c
cc PERM_PRINT prints a permutation.
c
c  Example:
c
c    Input:
c
c      P = 7 2 4 1 5 3 6
c
c    Printed output:
c
c      "This is the permutation:"
c
c      1 2 3 4 5 6 7
c      7 2 4 1 5 3 6
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 July 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of objects permuted.
c
c    Input, integer P(N), the permutation, in standard index form.
c
c    Input, character * ( * ) TITLE, a title.
c    If no title is supplied, then only the permutation is printed.
c
      implicit none

      integer n

      integer i
      integer ihi
      integer ilo
      integer inc
      parameter ( inc = 20 )
      integer p(n)
      character * ( * ) title

      if ( len_trim ( title ) .ne. 0 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) trim ( title )

        do ilo = 1, n, inc
          ihi = min ( n, ilo + inc - 1 )
          write ( *, '(a)' ) ' '
          write ( *, '(20i4)' ) ( i, i = ilo, ihi )
          write ( *, '(20i4)' ) ( p(i), i = ilo, ihi )
        end do

      else

        do ilo = 1, n, inc
          ihi = min ( n, ilo + inc - 1 )
          write ( *, '(20i4)' ) ( p(i), i = ilo, ihi )
        end do

      end if

      return
      end
      subroutine perm_tj_rank ( n, p, rank )

c*********************************************************************72
c
cc PERM_TJ_RANK computes the Trotter-Johnson rank of a permutation.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of values being permuted.
c    N must be positive.
c
c    Input, integer P(N), describes the permutation.
c    P(I) is the item which is permuted into the I-th place
c    by the permutation.
c
c    Output, integer RANK, the rank of the permutation.
c
      implicit none

      integer n

      integer i
      integer j
      integer k
      integer p(n)
      integer rank
c
c  Check.
c
      call perm_check ( n, p )

      rank = 0

      do j = 2, n

        k = 1
        i = 1

10      continue

        if ( p(i) .ne. j ) then
          if ( p(i) .lt. j ) then
            k = k + 1
          end if
          i = i + 1
          go to 10
        end if

        if ( mod ( rank, 2 ) .eq. 0 ) then
          rank = j * rank + j - k
        else
          rank = j * rank + k - 1
        end if

      end do

      return
      end
      subroutine perm_tj_successor ( n, p, rank )

c*********************************************************************72
c
cc PERM_TJ_SUCCESSOR computes the Trotter-Johnson permutation successor.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of values being permuted.
c    N must be positive.
c
c    Input/output, integer P(N), describes the permutation.
c    P(I) is the item which is permuted into the I-th place
c    by the permutation.
c
c    Input/output, integer RANK, the rank.
c    If RANK = -1 on input, then the routine understands that this is
c    the first call, and that the user wishes the routine to supply
c    the first element in the ordering, which has RANK = 0.
c    In general, the input value of RANK is increased by 1 for output,
c    unless the very last element of the ordering was input, in which
c    case the output value of RANK is 0.
c
      implicit none

      integer n

      integer d
      logical done
      integer i
      integer m
      integer p(n)
      integer par
      integer q(n)
      integer rank
      integer st
      integer temp
c
c  Return the first element.
c
      if ( rank .eq. -1 ) then
        call i4vec_indicator ( n, p )
        rank = 0
        return
      end if
c
c  Check.
c
      call perm_check ( n, p )

      st = 0
      do i = 1, n
        q(i) = p(i)
      end do
      done = .false.
      m = n

10    continue

      if ( 1 .lt. m .and. .not. done ) then

        d = 1

20      continue

        if ( q(d) .ne. m ) then
          d = d + 1
          go to 20
        end if

        do i = d, m - 1
          q(i) = q(i+1)
        end do

        call perm_parity ( m - 1, q, par )

        if ( par .eq. 1 ) then

          if ( d .eq. m ) then
            m = m - 1
          else
            temp      = p(st+d)
            p(st+d)   = p(st+d+1)
            p(st+d+1) = temp
            done = .true.
          end if

        else

          if ( d .eq. 1 ) then
            m = m - 1
            st = st + 1
          else
            temp      = p(st+d)
            p(st+d)   = p(st+d-1)
            p(st+d-1) = temp
            done = .true.
          end if

        end if

        go to 10

      end if
c
c  Last element was input.  Return first one.
c
      if ( m .eq. 1 ) then
        call i4vec_indicator ( n, p )
        rank = 0
        return
      end if

      rank = rank + 1

      return
      end
      subroutine perm_tj_unrank ( rank, n, p )

c*********************************************************************72
c
cc PERM_TJ_UNRANK computes the permutation of given Trotter-Johnson rank.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer RANK, the rank of the permutation.
c
c    Input, integer N, the number of values being permuted.
c    N must be positive.
c
c    Output, integer P(N), describes the permutation.
c
      implicit none

      integer n

      integer i
      integer i4_factorial
      integer j
      integer k
      integer jhi
      integer nperm
      integer p(n)
      integer rank
      integer r1
      integer r2
c
c  Check.
c
      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PERM_TJ_UNRANK - Fatal error!'
        write ( *, '(a)' ) '  Input N is illegal.'
        stop
      end if

      call perm_enum ( n, nperm )

      if ( rank .lt. 0 .or. nperm .lt. rank ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PERM_TJ_UNRANK - Fatal error!'
        write ( *, '(a)' ) '  The input rank is illegal.'
        stop
      end if

      p(1) = 1
      r2 = 0

      do j = 2, n
c
c  Replace this ratio of factorials!
c
        r1 = ( rank * i4_factorial ( j ) ) / i4_factorial ( n )
        k = r1 - j * r2

        if ( mod ( r2, 2 ) .eq. 0 ) then
          jhi = j - k
        else
          jhi = k + 1
        end if

        do i = j - 1, jhi, -1
          p(i+1) = p(i)
        end do
        p(jhi) = j

        r2 = r1

      end do

      return
      end
      subroutine perm_to_cycle ( n, p, ncycle, t, index )

c*********************************************************************72
c
cc PERM_TO_CYCLE converts a permutation from array to cycle form.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of values being permuted.
c    N must be positive.
c
c    Input, integer P(N), describes the permutation using a
c    single array.  For each index I, I -> P(I).
c
c    Output, integer NCYCLE, the number of cycles.
c    1 .le. NCYCLE .le. N.
c
c    Output, integer T(N), INDEX(N), describes the permutation
c    as a collection of NCYCLE cycles.  The first cycle is
c    T(1) -> T(2) -> ... -> T(INDEX(1)) -> T(1).
c
      implicit none

      integer n
      integer ncycle

      integer i
      integer index(n)
      integer j
      integer nset
      integer p(n)
      integer t(n)
c
c  Check.
c
      call perm_check ( n, p )
c
c  Initialize.
c
      ncycle = 0
      do i = 1, n
        index(i) = 0
      end do
      do i = 1, n
        t(i) = 0
      end do
      nset = 0
c
c  Find the next unused entry.
c
      do i = 1, n

        if ( 0 .lt. p(i) ) then

          ncycle = ncycle + 1
          index(ncycle) = 1

          nset = nset + 1
          t(nset) = p(i)
          p(i) = - p(i)

10        continue

            j = t(nset)

            if ( p(j) .lt. 0 ) then
              go to 20
            end if

            index(ncycle) = index(ncycle) + 1

            nset = nset + 1
            t(nset) = p(j)
            p(j) = - p(j)

          go to 10

20        continue

        end if

      end do
c
c  If no unused entries remain, we are done.
c  Restore the sign of the permutation and return.
c
      do i = 1, n
        p(i) = - p(i)
      end do

      return
      end
      subroutine pruefer_check ( n, p, ierror )

c*********************************************************************72
c
cc PRUEFER_CHECK checks a Pruefer code.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of nodes in the tree.
c    N must be at least 3.
c
c    Input, integer P(N-2), the Pruefer code for the tree.
c
c    Output, integer IERROR, error flag.
c    0, no error.
c    -1, N is less than 3.
c    J, the element P(J) is illegal.
c
      implicit none

      integer n

      integer i
      integer ierror
      integer p(n-2)

      ierror = 0

      if ( n .lt. 3 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PRUEFER_CHECK - Fatal error!'
        write ( *, '(a)' ) '  N .lt. 3.'
        stop
      end if

      do i = 1, n - 2
        if ( p(i) .lt. 1 .or. n .lt. p(i) ) then
          ierror = i
          return
        end if
      end do

      return
      end
      subroutine pruefer_enum ( n, ncode )

c*********************************************************************72
c
cc PRUEFER_ENUM enumerates the Pruefer codes on N-2 digits.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of digits in the code, plus 2.
c    N must be at least 3.
c
c    Output, integer NCODE, the number of distinct elements.
c
      implicit none

      integer n
      integer ncode

      ncode = n ** ( n - 2 )

      return
      end
      subroutine pruefer_rank ( n, p, rank )

c*********************************************************************72
c
cc PRUEFER_RANK ranks a Pruefer code.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of nodes in the tree.
c    N must be at least 3.
c
c    Input, integer P(N-2), the Pruefer code for the tree.
c
c    Output, integer RANK, the rank of the Pruefer code.
c
      implicit none

      integer n

      integer i
      integer ierror
      integer k
      integer p(n-2)
      integer rank
c
c  Check.
c
      call pruefer_check ( n, p, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PRUEFER_RANK - Fatal error!'
        write ( *, '(a)' ) '  Input array is illegal.'
        write ( *, '(a,i8)' ) '  Error code = ', ierror
        stop
      end if

      rank = 0
      k = 1
      do i = n - 2, 1, -1
        rank = rank + k * ( p(i) - 1 )
        k = k * n
      end do

      return
      end
      subroutine pruefer_successor ( n, p, rank )

c*********************************************************************72
c
cc PRUEFER_SUCCESSOR computes the lexical Pruefer sequence successor.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 February 2001
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of nodes in the tree.
c    N must be at least 3.
c
c    Input/output, integer P(N-2), on input, the Pruefer code
c    for a tree, and on output, its lexical successor.
c
c    Input/output, integer RANK, the rank.
c    If RANK = -1 on input, then the routine understands that this is
c    the first call, and that the user wishes the routine to supply
c    the first element in the ordering, which has RANK = 0.
c    In general, the input value of RANK is increased by 1 for output,
c    unless the very last element of the ordering was input, in which
c    case the output value of RANK is 0.
c
      implicit none

      integer n

      integer i
      integer ierror
      integer j
      integer rank
      integer p(n-2)
c
c  Return the first element.
c
      if ( rank .eq. -1 ) then
        do i = 1, n - 2
          p(i) = 1
        end do
        rank = 0
        return
      end if
c
c  Check.
c
      call pruefer_check ( n, p, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PRUEFER_SUCCESSOR - Fatal error!'
        write ( *, '(a)' ) '  The input array is illegal.'
        write ( *, '(a,i8)' ) '  IERROR = ', ierror
        stop
      end if

      j = n - 2

10    continue

        if ( p(j) .ne. n ) then
          go to 20
        end if

        j = j - 1

        if ( j .le. 0 ) then
          go to 20
        end if

      go to 10

20    continue

      if ( j .ne. 0 ) then
        p(j) = p(j) + 1
        do i = j + 1, n - 2
          p(i) = 1
        end do
        rank = rank + 1
      else
        do i = 1, n - 2
          p(i) = 1
        end do
        rank = 0
      end if

      return
      end
      subroutine pruefer_to_tree ( n, p, t )

c*********************************************************************72
c
cc PRUEFER_TO_TREE converts a Pruefer code to a tree.
c
c  Discussion:
c
c    The original code attempts to tack on an extra entry to P.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of nodes in the tree.
c    N must be at least 3.
c
c    Input, integer P(N-2), the Pruefer code for the tree.
c
c    Output, integer T(2,N-1), describes the edges of the tree
c    as pairs of nodes.
c
      implicit none

      integer n

      integer d(n)
      integer i
      integer ierror
      integer j
      integer p(n-2)
      integer t(2,n-1)
      integer x
      integer y
c
c  Check.
c
      call pruefer_check ( n, p, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PRUEFER_TO_TREE - Fatal error!'
        write ( *, '(a)' ) '  The input array is illegal!'
        write ( *, '(a,i8)' ) '  IERROR = ', ierror
        stop
      end if
c
c  Initialize the tree to 0.
c
      do j = 1, n - 1
        do i = 1, 2
          t(i,j) = 0
        end do
      end do

      do i = 1, n
        d(i) = 1
      end do

      do i = 1, n - 2
        d(p(i)) = d(p(i)) + 1
      end do

      do i = 1, n - 1

        x = n

10      continue

        if ( d(x) .ne. 1 ) then
          x = x - 1
          go to 10
        end if

        if ( i .eq. n - 1 ) then
          y = 1
        else
          y = p(i)
        end if

        d(x) = d(x) - 1
        d(y) = d(y) - 1

        t(1,i) = x
        t(2,i) = y

      end do

      return
      end
      subroutine pruefer_unrank ( rank, n, p )

c*********************************************************************72
c
cc PRUEFER_UNRANK unranks a Pruefer code.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer RANK, the rank of the Pruefer code.
c
c    Input, integer N, the number of nodes in the tree.
c    N must be at least 3.
c
c    Output, integer P(N-2), the Pruefer code for the tree.
c
      implicit none

      integer n

      integer i
      integer ncode
      integer p(n-2)
      integer rank
      integer rank_copy
c
c  Check.
c
      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PRUEFER_UNRANK - Fatal error!'
        write ( *, '(a)' ) '  Input N is illegal.'
        stop
      end if

      call pruefer_enum ( n, ncode )

      if ( rank .lt. 0 .or. ncode .lt. rank ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PRUEFER_UNRANK - Fatal error!'
        write ( *, '(a)' ) '  The input rank is illegal.'
        stop
      end if

      rank_copy = rank

      do i = n - 2, 1, -1
        p(i) = mod ( rank_copy, n ) + 1
        rank_copy = ( rank_copy - p(i) + 1 ) / n
      end do

      return
      end
      subroutine queens ( n, iarray, k, nstack, istack, maxstack )

c*********************************************************************72
c
cc QUEENS finds possible positions for the K-th nonattacking queen.
c
c  Discussion:
c
c    The chessboard is N by N, and is being filled one column at a time,
c    with a tentative solution to the nonattacking queen problem.  So
c    far, K-1 rows have been chosen, and we now need to provide a list
c    of all possible rows that might be used in column K.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 July 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the total number of queens to place, and
c    the length of a side of the chessboard.
c
c    Input, integer IARRAY(N).  The first K-1 entries of IARRAY
c    record the rows into which queens have already been placed.
c
c    Input, integer K, the column for which we need possible
c    row positions for the next queen.
c
c    Input/output, integer NSTACK, the current length of stack.
c    On output, this has been updated.
c
c    Input/output, integer ISTACK(MAXSTACK).  On output, we
c    have added the candidates, and the number of candidates, to the end
c    of the stack.
c
c    Input, integer MAXSTACK, maximum dimension of ISTACK.
c
      implicit none

      integer n
      integer maxstack

      logical diag
      integer iarray(n)
      integer irow
      integer istack(maxstack)
      integer jcol
      integer k
      integer ncan
      integer nstack
      logical row

      ncan = 0

      do irow = 1, n
c
c  If row IROW has already been used, that's it.
c
        row = .false.

        do jcol = 1, k - 1
          if ( iarray(jcol) .eq. irow ) then
            row = .true.
          end if
        end do

        if ( .not. row ) then

          diag = .false.

          do jcol = 1, k - 1

            if ( irow .eq. iarray(jcol) + k - jcol .or. 
     &           irow .eq. iarray(jcol) - ( k - jcol ) ) then

              diag = .true.

            end if

          end do

          if ( .not. diag ) then
            ncan = ncan + 1
            nstack = nstack + 1
            istack(nstack) = irow
          end if

        end if

      end do

      nstack = nstack + 1
      istack(nstack) = ncan

      return
      end
      function r8_gamma_log ( x )

c*********************************************************************72
c
cc R8_GAMMA_LOG evaluates log ( Gamma ( X ) ) for a real argument.
c
c  Discussion:
c
c    This routine calculates the LOG(GAMMA) function for a positive real
c    argument X.  Computation is based on an algorithm outlined in
c    references 1 and 2.  The program uses rational functions that
c    theoretically approximate LOG(GAMMA) to at least 18 significant
c    decimal digits.  The approximation for X > 12 is from reference
c    3, while approximations for X < 12.0 are similar to those in
c    reference 1, but are unpublished.
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody, Laura Stoltz
c
c  Reference:
c
c    William Cody, Kenneth Hillstrom,
c    Chebyshev Approximations for the Natural Logarithm of the 
c    Gamma Function,
c    Mathematics of Computation,
c    Volume 21, Number 98, April 1967, pages 198-203.
c
c    Kenneth Hillstrom,
c    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
c    May 1969.
c
c    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, 
c    Charles Mesztenyi, John Rice, Henry Thatcher, 
c    Christoph Witzgall,
c    Computer Approximations,
c    Wiley, 1968,
c    LC: QA297.C64.
c
c  Parameters:
c
c    Input, double precision X, the argument of the function.
c
c    Output, double precision R8_GAMMA_LOG, the value of the function.
c
      implicit none

      double precision c(7)
      double precision corr
      double precision d1
      double precision d2
      double precision d4
      double precision eps
      double precision frtbig
      double precision four
      double precision half
      integer i
      double precision one
      double precision pnt68
      double precision p1(8)
      double precision p2(8)
      double precision p4(8)
      double precision q1(8)
      double precision q2(8)
      double precision q4(8)
      double precision r8_gamma_log
      double precision res
      double precision sqrtpi
      double precision thrhal
      double precision twelve
      double precision two
      double precision x
      double precision xbig
      double precision xden
      double precision xinf
      double precision xm1
      double precision xm2
      double precision xm4
      double precision xnum
      double precision y
      double precision ysq
      double precision zero
c
c  Mathematical constants
c
      data one /1.0D+00/
      data half /0.5D+00/
      data twelve /12.0D+00/
      data zero /0.0D+00/
      data four /4.0D+00/
      data thrhal /1.5D+00/
      data two /2.0D+00/
      data pnt68 /0.6796875D+00/
      data sqrtpi /0.9189385332046727417803297D+00/
c
c  Machine dependent parameters
c
      data xbig /2.55D+305/
      data xinf /1.79D+308/
      data eps /2.22D-16/
      data frtbig /2.25D+76/
c
c  Numerator and denominator coefficients for rational minimax
c  approximation over (0.5,1.5).
c
      data d1/-5.772156649015328605195174D-01/
      data p1/
     &   4.945235359296727046734888D+00,
     &   2.018112620856775083915565D+02,
     &   2.290838373831346393026739D+03,
     &   1.131967205903380828685045D+04,
     &   2.855724635671635335736389D+04,
     &   3.848496228443793359990269D+04,
     &   2.637748787624195437963534D+04,
     &   7.225813979700288197698961D+03/
      data q1/
     &   6.748212550303777196073036D+01,
     &   1.113332393857199323513008D+03,
     &   7.738757056935398733233834D+03,
     &   2.763987074403340708898585D+04,
     &   5.499310206226157329794414D+04,
     &   6.161122180066002127833352D+04,
     &   3.635127591501940507276287D+04,
     &   8.785536302431013170870835D+03/
c
c  Numerator and denominator coefficients for rational minimax
c  Approximation over (1.5,4.0).
c
      data d2/4.227843350984671393993777D-01/
      data p2/
     &   4.974607845568932035012064D+00,
     &   5.424138599891070494101986D+02,
     &   1.550693864978364947665077D+04,
     &   1.847932904445632425417223D+05,
     &   1.088204769468828767498470D+06,
     &   3.338152967987029735917223D+06,
     &   5.106661678927352456275255D+06,
     &   3.074109054850539556250927D+06/
      data q2/
     &   1.830328399370592604055942D+02,
     &   7.765049321445005871323047D+03,
     &   1.331903827966074194402448D+05,
     &   1.136705821321969608938755D+06,
     &   5.267964117437946917577538D+06,
     &   1.346701454311101692290052D+07,
     &   1.782736530353274213975932D+07,
     &   9.533095591844353613395747D+06/
c
c  Numerator and denominator coefficients for rational minimax
c  Approximation over (4.0,12.0).
c
      data d4/1.791759469228055000094023D+00/
      data p4/
     &   1.474502166059939948905062D+04,
     &   2.426813369486704502836312D+06,
     &   1.214755574045093227939592D+08,
     &   2.663432449630976949898078D+09,
     &   2.940378956634553899906876D+10,
     &   1.702665737765398868392998D+11,
     &   4.926125793377430887588120D+11,
     &   5.606251856223951465078242D+11/
      data q4/
     &   2.690530175870899333379843D+03,
     &   6.393885654300092398984238D+05,
     &   4.135599930241388052042842D+07,
     &   1.120872109616147941376570D+09,
     &   1.488613728678813811542398D+10,
     &   1.016803586272438228077304D+11,
     &   3.417476345507377132798597D+11,
     &   4.463158187419713286462081D+11/
c
c  Coefficients for minimax approximation over (12, INF).
c
      data c/
     &  -1.910444077728D-03,
     &   8.4171387781295D-04,
     &  -5.952379913043012D-04,
     &   7.93650793500350248D-04,
     &  -2.777777777777681622553D-03,
     &   8.333333333333333331554247D-02,
     &   5.7083835261D-03/

      y = x

      if ( zero .lt. y .and. y .le. xbig ) then

        if ( y .le. eps ) then

          res = - dlog ( y )
c
c  EPS < X <= 1.5.
c
        else if ( y .le. thrhal ) then

          if ( y .lt. pnt68 ) then
            corr = - dlog ( y )
            xm1 = y
          else
            corr = zero
            xm1 = ( y - half ) - half
          end if

          if ( y .le. half .or. pnt68 .le. y ) then

            xden = one
            xnum = zero
            do i = 1, 8
              xnum = xnum * xm1 + p1(i)
              xden = xden * xm1 + q1(i)
            end do

            res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) )

          else

            xm2 = ( y - half ) - half
            xden = one
            xnum = zero
            do i = 1, 8
              xnum = xnum * xm2 + p2(i)
              xden = xden * xm2 + q2(i)
            end do
            res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) )

          end if
c
c  1.5 < X <= 4.0.
c
        else if ( y .le. four ) then

          xm2 = y - two
          xden = one
          xnum = zero
          do i = 1, 8
            xnum = xnum * xm2 + p2(i)
            xden = xden * xm2 + q2(i)
          end do

          res = xm2 * ( d2 + xm2 * ( xnum / xden ) )
c
c  4.0 < X <= 12.0.
c
        else if ( y .le. twelve ) then

          xm4 = y - four
          xden = -one
          xnum = zero
          do i = 1, 8
            xnum = xnum * xm4 + p4(i)
            xden = xden * xm4 + q4(i)
          end do

          res = d4 + xm4 * ( xnum / xden )
c
c  Evaluate for 12 <= argument.
c
        else

          res = zero

          if ( y .le. frtbig ) then

            res = c(7)
            ysq = y * y

            do i = 1, 6
              res = res / ysq + c(i)
            end do

          end if

          res = res / y
          corr = dlog ( y )
          res = res + sqrtpi - half * corr
          res = res + y * ( corr - one )

        end if
c
c  Return for bad arguments.
c
      else

        res = xinf

      end if
c
c  Final adjustments and return.
c
      r8_gamma_log = res

      return
      end
      subroutine r8vec_backtrack ( n, x, indx, k, nstack, stack, 
     &  maxstack, ncan )

c*********************************************************************72
c
cc R8VEC_BACKTRACK supervises a backtrack search for an R8VEC.
c
c  Discussion:
c
c    The routine tries to construct a real vector one index at a time,
c    using possible candidates as supplied by the user.
c
c    At any time, the partially constructed vector may be discovered to be
c    unsatisfactory, but the routine records information about where the
c    last arbitrary choice was made, so that the search can be
c    carried out efficiently, rather than starting out all over again.
c
c    First, call the routine with INDX = 0 so it can initialize itself.
c
c    Now, on each return from the routine, if INDX is:
c      1, you've just been handed a complete candidate vector;
c         Admire it, analyze it, do what you like.
c      2, please determine suitable candidates for position X(K).
c         Return the number of candidates in NCAN(K), adding each
c         candidate to the end of STACK, and increasing NSTACK.
c      3, you're done.  Stop calling the routine;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 July 2000
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of positions to be filled in
c    the vector.
c
c    Input/output, double precision X(N), the partial or complete
c    candidate vector.
c
c    Input/output, integer INDX, a communication flag.
c    On input,
c      0 to start a search.
c    On output:
c      1, a complete output vector has been determined and returned in X(1:N);
c      2, candidates are needed for position X(K);
c      3, no more possible vectors exist.
c
c    Input/output, integer K, if INDX=2, the current vector index
c    being considered.
c
c    Input/output, integer NSTACK, the current length of the stack.
c
c    Input/output, double precision STACK(MAXSTACK), a list of all current
c    candidates for all positions 1 through K.
c
c    Input, integer MAXSTACK, the maximum length of the stack.
c
c    Input/output, integer NCAN(N), lists the current number
c    of candidates for positions 1 through K.
c
      implicit none

      integer n
      integer maxstack

      integer indx
      integer k
      integer ncan(n)
      integer nstack
      double precision stack(maxstack)
      double precision x(n)
c
c  If this is the first call, request a candidate for position 1.
c
      if ( indx .eq. 0 ) then
        k = 1
        nstack = 0
        indx = 2
        return
      end if
c
c  Examine the stack.
c
10    continue
c
c  If there are candidates for position K, take the first available
c  one off the stack, and increment K.
c
c  This may cause K to reach the desired value of N, in which case
c  we need to signal the user that a complete set of candidates
c  is being returned.
c
        if ( 0 .lt. ncan(k) ) then

          x(k) = stack(nstack)
          nstack = nstack - 1

          ncan(k) = ncan(k) - 1

          if ( k .ne. n ) then
            k = k + 1
            indx = 2
          else
            indx = 1
          end if

          go to 20
c
c  If there are no candidates for position K, then decrement K.
c  If K is still positive, repeat the examination of the stack.
c
        else

          k = k - 1

          if ( k .le. 0 ) then
            indx = 3
            go to 20
          end if

        end if

      go to 10

20    continue

      return
      end
      function r8vec_dot_product ( n, v1, v2 )

c*********************************************************************72
c
cc R8VEC_DOT_PRODUCT finds the dot product of a pair of R8VEC's.
c
c  Discussion:
c
c    An R8VEC is a vector of R8 values.
c
c    In FORTRAN90, the system routine DOT_PRODUCT should be called
c    directly.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of the vectors.
c
c    Input, double precision V1(N), V2(N), the vectors.
c
c    Output, double precision R8VEC_DOT_PRODUCT, the dot product.
c
      implicit none

      integer n

      integer i
      double precision r8vec_dot_product
      double precision v1(n)
      double precision v2(n)
      double precision value

      value = 0.0D+00
      do i = 1, n
        value = value + v1(i) * v2(i)
      end do

      r8vec_dot_product = value

      return
      end
      subroutine rgf_check ( m, f, ierror )

c*********************************************************************72
c
cc RGF_CHECK checks a restricted growth function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer M, the domain of the RGF is the integers
c    from 1 to M.  M must be positive.
c
c    Input, integer F(M), the restricted growth function.
c
c    Output, integer IERROR, error flag.
c    0, no error.
c    -1, M is illegal.
c    I, entry I of the restricted growth function is illegal.
c
      implicit none

      integer m

      integer f(m)
      integer fmax
      integer i
      integer ierror

      ierror = 0

      if ( m .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RGF_CHECK - Fatal error!'
        write ( *, '(a)' ) '  M .lt. 1.'
        stop
      end if

      fmax = 0
      do i = 1, m
        if ( f(i) .le. 0 .or. fmax + 1 .lt. f(i) ) then
          ierror = i
          return
        end if
        fmax = max ( fmax, f(i) )
      end do

      return
      end
      subroutine rgf_enum ( m, nrgf )

c*********************************************************************72
c
cc RGF_ENUM enumerates the restricted growth functions on M.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 April 2001
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer M, the domain of the RGF is the integers
c    from 1 to M.  M must be positive.  However, for the enumeration routine
c    only, it is legal to call with any value of M.
c
c    Output, integer NRGF, the number of restricted growth
c    functions.
c
      implicit none

      integer m

      integer b(0:m)
      integer binomial
      integer i
      integer j
      integer nrgf

      if ( m .lt. 0 ) then

        nrgf = 0

      else if ( m .eq. 0 ) then

        nrgf = 1

      else

        b(0) = 1
        do j = 1, m
          b(j) = 0
          do i = 0, j - 1
            b(j) = b(j) + binomial ( j - 1, i ) * b(i)
          end do
        end do

        nrgf = b(m)

      end if

      return
      end
      subroutine rgf_g_table ( m, mmax, d )

c*********************************************************************72
c
cc RGF_G_TABLE tabulates the generalized restricted growth functions.
c
c  Example:
c
c    M = 6
c
c    D =  1    1    1    1    1    1    1
c         1    2    3    4    5    6    0
c         2    5   10   17   26    0    0
c         5   15   37   77    0    0    0
c        15   52  151    0    0    0    0
c        52  203    0    0    0    0    0
c       203    0    0    0    0    0    0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer M, indicates how many rows and columns are to
c    be computed.  M must be nonnegative.
c
c    Input, integer MMAX, the value used to allocate space for the
c    D array.  MMAX must be at least M.
c
c    Output, integer D(0:MMAX,0:MMAX), the first M+1 rows and
c    M+1 columns of the table of the number of generalized restricted growth
c    functions.  D(I,J) is the number of GRGF's of length I with restriction
c    parameter J.
c
      implicit none

      integer mmax

      integer d(0:mmax,0:mmax)
      integer i
      integer j
      integer m

      do j = 0, m
        d(0,j) = 1
      end do

      do i = 1, m
        do j = 0, m
          if ( j .le. m - i ) then
            d(i,j) = j * d(i-1,j) + d(i-1,j+1)
          else
            d(i,j) = 0
          end if
        end do
      end do

      return
      end
      subroutine rgf_rank ( m, f, rank )

c*********************************************************************72
c
cc RGF_RANK ranks a restricted growth function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer M, the domain of the RGF is the integers
c    from 1 to M.  M must be positive.
c
c    Input, integer F(M), the restricted growth function.
c
c    Output, integer RANK, the rank of the restricted growth
c    function.
c
      implicit none

      integer m

      integer d(0:m,0:m)
      integer f(m)
      integer i
      integer ierror
      integer j
      integer rank
c
c  Check.
c
      call rgf_check ( m, f, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RGF_RANK - Fatal error!'
        write ( *, '(a)' ) '  The input array is illegal!'
        write ( *, '(a,i8)' ) '  IERROR = ', ierror
        stop
      end if
c
c  Get the generalized restricted growth function table.
c
      call rgf_g_table ( m, m, d )

      rank = 0
      j = 1
      do i = 2, m
        rank = rank + ( f(i) - 1 ) * d(m-i,j)
        j = max ( j, f(i) )
      end do

      return
      end
      subroutine rgf_successor ( m, f, rank )

c*********************************************************************72
c
cc RGF_SUCCESSOR generates the next restricted growth function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer M, the domain of the RGF is the integers
c    from 1 to M.  M must be positive.
c
c    Input/output, integer F(M), the restricted growth function.
c
c    Input/output, integer RANK, the rank.
c    If RANK = -1 on input, then the routine understands that this is
c    the first call, and that the user wishes the routine to supply
c    the first element in the ordering, which has RANK = 0.
c    In general, the input value of RANK is increased by 1 for output,
c    unless the very last element of the ordering was input, in which
c    case the output value of RANK is 0.
c
      implicit none

      integer m

      integer f(m)
      integer fmax
      integer i
      integer ierror
      integer j
      integer rank
c
c  Return the first element.
c
      if ( rank .eq. -1 ) then
        do i = 1, m
          f(i) = 1
        end do
        rank = 0
        return
      end if
c
c  Check.
c
      call rgf_check ( m, f, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RGF_SUCCESSOR - Fatal error!'
        write ( *, '(a)' ) '  The input array is illegal!'
        write ( *, '(a,i8)' ) '  IERROR = ', ierror
        stop
      end if
c
c  Find the first position from the right which can be incremented.
c
      do i = m, 2, -1

        fmax = 1
        do j = 2, i - 1
          fmax = max ( fmax, f(j) )
        end do
c
c  Increment the function at this position, and set later entries to 1.
c
        if ( f(i) .ne. fmax + 1 ) then
          f(i) = f(i) + 1
          do j = i + 1, m
            f(j) = 1
          end do
          rank = rank + 1
          return
        end if

      end do
c
c  The final element was input.
c  Return the first element.
c
      do j = 1, m
        f(j) = 1
      end do
      rank = 0

      return
      end
      subroutine rgf_to_setpart ( m, f, nsub, s, index )

c*********************************************************************72
c
cc RGF_TO_SETPART converts a restricted growth function to a set partition.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 December 2001
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer M, the domain of the RGF is the integers
c    from 1 to M.  M must be positive.
c
c    Input, integer F(M), the restricted growth function.
c
c    Output, integer NSUB, the number of nonempty subsets into
c    which the set is partitioned.
c
c    Output, integer S(M), describes the partition of a set of
c    M objects into NSUB nonempty subsets.  If element I of the
c    superset belongs to subset J, then S(I) = J.
c
c    Output, integer INDEX(M), lists the location in S of the last
c    element of each subset.  Thus, the elements of subset 1
c    are S(1) through S(INDEX(1)), the elements of subset 2
c    are S(INDEX(1)+1) through S(INDEX(2)) and so on.
c
      implicit none

      integer m

      integer f(m)
      integer i
      integer ierror
      integer index(m)
      integer j
      integer k
      integer nsub
      integer s(m)
c
c  Check.
c
      call rgf_check ( m, f, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RGF_TO_SETPART - Fatal error!'
        write ( *, '(a)' ) '  The input array is illegal!'
        write ( *, '(a,i8)' ) '  IERROR = ', ierror
        stop
      end if
c
c  Determine the number of subsets.
c
      call i4vec_max ( m, f, nsub )
c
c  Initialize.
c
      do i = 1, m
        s(i) = 0
      end do

      do i = 1, nsub
        index(i) = 0
      end do
c
c  For each subset I, collect the indices of F which have value I.
c  These are the elements of the I-th subset.
c
      k = 0
      do i = 1, nsub
        do j = 1, m
          if ( f(j) .eq. i ) then
            k = k + 1
            s(k) = j
          end if
        end do
        index(i) = k
      end do

      return
      end
      subroutine rgf_unrank ( rank, m, f )

c*********************************************************************72
c
cc RGF_UNRANK returns the restricted growth function of a given rank.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 April 2001
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer RANK, the rank of the restricted growth
c    function.
c
c    Input, integer M, the domain of the RGF is the integers
c    from 1 to M.  M must be positive.
c
c    Output, integer F(M), the restricted growth function.
c
      implicit none

      integer m

      integer d(0:m,0:m)
      integer f(m)
      integer i
      integer j
      integer nrgf
      integer rank
      integer rank_copy
c
c  Check.
c
      if ( m .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RGF_UNRANK - Fatal error!'
        write ( *, '(a)' ) '  Input M is illegal.'
        stop
      end if

      call rgf_enum ( m, nrgf )

      if ( rank .lt. 0 .or. nrgf .lt. rank ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RGF_UNRANK - Fatal error!'
        write ( *, '(a)' ) '  The input rank is illegal.'
        stop
      end if
c
c  Get the generalized restricted growth function table.
c
      call rgf_g_table ( m, m, d )

      rank_copy = rank
      j = 1
      f(1) = 1

      do i = 2, m

        if ( j * d(m-i,j) .le. rank_copy ) then
          f(i) = j + 1
          rank_copy = rank_copy - j * d(m-i,j)
          j = j + 1
        else
          f(i) = 1 + ( rank_copy / d(m-i,j) )
          rank_copy = mod ( rank_copy, d(m-i,j) )
        end if

      end do

      return
      end
      subroutine setpart_check ( m, nsub, s, index, ierror )

c*********************************************************************72
c
cc SETPART_CHECK checks a set partition.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer M, the number of elements of the set.
c    M must be positive.
c
c    Input, integer NSUB, the number of nonempty subsets into
c    which the set is partitioned.  1 .le. NSUB .le. M.
c
c    Input, integer INDEX(NSUB), lists the location in S of the
c    last element of each subset.  Thus, the elements of subset 1
c    are S(1) through S(INDEX(1)), the elements of subset 2
c    are S(INDEX(1)+1) through S(INDEX(2)) and so on.
c
c    Input, integer S(M), contains the integers from 1 to M,
c    grouped into subsets as described by INDEX.
c
c    Output, integer IERROR, error flag.
c    0, no error.
c    -I, the I-th element of INDEX is illegal.
c    +I, the I-th element of S is illegal.
c
      implicit none

      integer m
      integer nsub

      integer i
      integer ierror
      integer imin
      integer index(nsub)
      integer j
      integer s(m)

      ierror = 0
c
c  Check M.
c
      if ( m .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SETPART_CHECK - Fatal error!'
        write ( *, '(a)' ) '  M .lt. 1.'
        stop
      end if
c
c  Check INDEX.
c
      imin = 0
      do i = 1, nsub
        if ( index(i) .le. imin .or. m .lt. index(i) ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SETPART_CHECK - Fatal error!'
          write ( *, '(a)' ) '  INDEX(I) <= IMIN or'
          write ( *, '(a)' ) '  M < INDEX(I).'
          write ( *, '(a,i4)' ) '  IMIN = ', imin
          write ( *, '(a,i4)' ) '  M = ', m
          write ( *, '(a,i4,a,i4)' ) '  INDEX(', i, ') = ', index(i)
          ierror = -i
          return
        end if
        imin = index(i)
      end do
c
c  Check the elements of S.
c
      do i = 1, nsub

        if ( s(i) .le. 0 .or. m .lt. s(i) ) then
          ierror = i
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SETPART_CHECK - Fatal error!'
          write ( *, '(a)' ) '  S(I) <= 0 .or. M < S(I)'
          write ( *, '(a,i4)' ) '  M = ', m
          write ( *, '(a,i4,a,i4)' ) '  S(', i, ') = ', s(i)
          return
        end if

        do j = 1, i - 1
          if ( s(j) .eq. s(i) ) then
            ierror = i
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'SETPART_CHECK - Fatal error!'
            write ( *, '(a,i4)' ) '  S(J) == S(I)'
            write ( *, '(a,i4,a,i4)' ) '  S(', i, ') = ', s(i)
            write ( *, '(a,i4,a,i4)' ) '  S(', j, ') = ', s(j)
            return
          end if
        end do

      end do

      return
      end
      subroutine setpart_enum ( m, npart )

c*********************************************************************72
c
cc SETPART_ENUM enumerates the partitions of a set of M elements.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 April 2001
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer M, the number of elements in the set.
c    M must be positive.  However, for the enumeration routine only,
c    it is legal to call with any value of M.
c
c    Output, integer NPART, the number of partitions of the set.
c
      implicit none

      integer m

      integer b(0:m)
      integer binomial
      integer i
      integer j
      integer npart

      if ( m .lt. 0 ) then

        npart = 0

      else if ( m .eq. 0 ) then

        npart = 1

      else

        b(0) = 1
        do j = 1, m
          b(j) = 0
          do i = 0, j - 1
            b(j) = b(j) + binomial ( j - 1, i ) * b(i)
          end do
        end do

        npart = b(m)

      end if

      return
      end
      subroutine setpart_to_rgf ( m, nsub, s, index, f )

c*********************************************************************72
c
cc SETPART_TO_RGF converts a set partition to a restricted growth function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer M, the number of elements of the set.
c    M must be positive.
c
c    Input, integer NSUB, the number of nonempty subsets into
c    which the set is partitioned.  1 .le. NSUB .le. M.
c
c    Input, integer INDEX(NSUB), lists the location in S of the
c    last element of each subset.  Thus, the elements of subset 1
c    are S(1) through S(INDEX(1)), the elements of subset 2
c    are S(INDEX(1)+1) through S(INDEX(2)) and so on.
c
c    Input, integer S(M), contains the integers from 1 to M,
c    grouped into subsets as described by INDEX.
c
c    Output, integer F(M), the restricted growth function from
c    M to NSUB.
c
      implicit none

      integer m
      integer nsub

      integer f(m)
      integer i
      integer ierror
      integer index(nsub)
      integer k
      integer khi
      integer klo
      integer s(m)
c
c  Check.
c
      call setpart_check ( m, nsub, s, index, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SETPART_TO_RGF - Fatal error!'
        write ( *, '(a)' ) '  The input array is illegal.'
        stop
      end if

      khi = 0
      do i = 1, nsub
        klo = khi + 1
        khi = index(i)
        do k = klo, khi
          f(s(k)) = i
        end do
      end do

      return
      end
      subroutine stirling_numbers1 ( m, n, s )

c*********************************************************************72
c
cc STIRLING_NUMBERS1 computes Stirling numbers of the first kind.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer M, the maximum row to compute.
c    M must be nonnegative.
c
c    Input, integer N, the maximum column to compute.
c    N must be nonnegative.
c
c    Output, integer S(0:M,0:N), the first M+1 rows and N+1 columns
c    of the table of Stirling numbers of the first kind.
c
      implicit none

      integer m
      integer n

      integer i
      integer j
      integer s(0:m,0:n)

      do j = 0, n
        do i = 0, m
          s(i,j) = 0
        end do
      end do

      s(0,0) = 1

      do i = 1, m
        do j = 1, n
          if ( j .le. i ) then
            s(i,j) = s(i-1,j-1) - ( i - 1 ) * s(i-1,j)
          end if
        end do
      end do

      return
      end
      subroutine stirling_numbers2 ( m, n, s )

c*********************************************************************72
c
cc STIRLING_NUMBERS2 computes Stirling numbers of the second kind.
c
c  Discussion:
c
c    The reference has a typographical error, referring to
c    S(I-J,J-1) instead of S(I-1,J-1).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer M, the maximum row to compute.
c    M must be nonnegative.
c
c    Input, integer N, the maximum column to compute.
c    N must be nonnegative.
c
c    Output, integer S(0:M,0:N), the first M+1 rows and N+1 columns
c    of the table of Stirling numbers of the second kind.
c
      implicit none

      integer m
      integer n

      integer i
      integer j
      integer s(0:m,0:n)

      do j = 0, n
        do i = 0, m
          s(i,j) = 0
        end do
      end do

      s(0,0) = 1

      do i = 1, m
        do j = 1, n
          if ( j .le. i ) then
            s(i,j) = j * s(i-1,j) + s(i-1,j-1)
          end if
        end do
      end do

      return
      end
      subroutine subset_check ( n, t )

c*********************************************************************72
c
cc SUBSET_CHECK checks a subset.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of elements in the master set.
c    N must be positive.
c
c    Input, integer T(N), the subset.  If T(I) = 0, item I is
c    not in the subset; if T(I) = 1, item I is in the subset.
c
      implicit none

      integer n

      integer i
      integer t(n)

      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SUBSET_CHECK - Fatal error!'
        write ( *, '(a)' ) '  N .lt. 1.'
        stop
      end if

      do i = 1, n

        if ( t(i) .ne. 0 .and. t(i) .ne. 1 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SUBSET_CHECK - Fatal error!'
          write ( *, '(a)' ) '  T(I) is not 0 and not 1.'
          stop
        end if

      end do

      return
      end
      subroutine subset_colex_rank ( n, t, rank )

c*********************************************************************72
c
cc SUBSET_COLEX_RANK computes the colexicographic rank of a subset.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 August 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of items in the master set.
c    N must be positive.
c
c    Input, integer T(N), the subset.  If T(I) = 0, item I is
c    not in the subset; if T(I) = 1, item I is in the subset.
c
c    Output, integer RANK, the rank of the subset.
c
      implicit none

      integer n

      integer i
      integer rank
      integer t(n)
c
c  Check.
c
      call subset_check ( n, t )

      rank = 0

      do i = 1, n

        if ( t(i) .eq. 1 ) then
          rank = rank + 2 ** ( i - 1 )
        end if

      end do

      return
      end
      subroutine subset_colex_successor ( n, t, rank )

c*********************************************************************72
c
cc SUBSET_COLEX_SUCCESSOR computes the subset colexicographic successor.
c
c  Discussion:
c
c    In the original code, there is a last element with no successor.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of elements in the master set.
c    N must be positive.
c
c    Input/output, integer T(N), describes a subset.  T(I) is 0 if
c    the I-th element of the master set is not in the subset, and is
c    1 if the I-th element is part of the subset.
c    On input, T describes a subset.
c    On output, T describes the next subset in the ordering.
c    If the input T was the last in the ordering, then the output T
c    will be the first.
c
c    Input/output, integer RANK, the rank.
c    If RANK = -1 on input, then the routine understands that this is
c    the first call, and that the user wishes the routine to supply
c    the first element in the ordering, which has RANK = 0.
c    In general, the input value of RANK is increased by 1 for output,
c    unless the very last element of the ordering was input, in which
c    case the output value of RANK is 0.
c
      implicit none

      integer n

      integer i
      integer rank
      integer t(n)
c
c  Return the first element.
c
      if ( rank .eq. -1 ) then
        do i = 1, n
          t(i) = 0
        end do
        rank = 0
        return
      end if
c
c  Check.
c
      call subset_check ( n, t )

      do i = 1, n

        if ( t(i) .eq. 0 ) then
          t(i) = 1
          rank = rank + 1
          return
        else
          t(i) = 0
        end if

      end do

      rank = 0

      return
      end
      subroutine subset_colex_unrank ( rank, n, t )

c*********************************************************************72
c
cc SUBSET_COLEX_UNRANK computes the subset of given colexicographic rank.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer RANK, the rank of the subset.
c
c    Input, integer N, the number of items in the master set.
c    N must be positive.
c
c    Output, integer T(N), the subsetof the given rank.
c    If T(I) = 0, item I is not in the subset; if T(I) = 1, item I is
c    in the subset.
c
      implicit none

      integer n

      integer i
      integer nsub
      integer rank
      integer rank_copy
      integer t(n)
c
c  Check.
c
      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SUBSET_COLEX_UNRANK - Fatal error!'
        write ( *, '(a)' ) '  Input N is illegal.'
        stop
      end if

      call subset_enum ( n, nsub )

      if ( rank .lt. 0 .or. nsub .lt. rank ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SUBSET_COLEX_UNRANK - Fatal error!'
        write ( *, '(a)' ) '  The input rank is illegal.'
        stop
      end if

      rank_copy = rank

      do i = 1, n
        if ( mod ( rank_copy, 2 ) .eq. 1 ) then
          t(i) = 1
        else
          t(i) = 0
        end if

        rank_copy = rank_copy / 2

      end do

      return
      end
      subroutine subset_complement ( n, a, b )

c*********************************************************************72
c
cc SUBSET_COMPLEMENT computes the complement of a set.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the order of the master set, of which A is
c    a subset.  N must be positive.
c
c    Input, integer A(N), a subset of the master set.
c    A(I) = 0 if the I-th element is in the subset A, and is
c    1 otherwise.
c
c    Output, integer B(N), the complement of A.
c
      implicit none

      integer n

      integer a(n)
      integer b(n)
      integer i
c
c  Check.
c
      call subset_check ( n, a )

      do i = 1, n
        b(i) = 1 - a(i)
      end do

      return
      end
      subroutine subset_distance ( n, t1, t2, dist )

c*********************************************************************72
c
cc SUBSET_DISTANCE computes the Hamming distance between two sets.
c
c  Discussion:
c
c    The sets T1 and T2 are assumed to be subsets of a set of N elements.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the order of the master set, of which T1 and
c    T2 are subsets.  N must be positive.
c
c    Input, integer T1(N), T2(N), two subsets of the master set.
c    T1(I) = 0 if the I-th element is in the subset T1, and is
c    1 otherwise; T2 is defined similarly.
c
c    Output, integer DIST, the Hamming distance between T1 and T2,
c    defined as the number of elements of the master set which are
c    in either T1 or T2 but not both.
c
      implicit none

      integer n

      integer dist
      integer i
      integer t1(n)
      integer t2(n)
c
c  Check.
c
      call subset_check ( n, t1 )

      call subset_check ( n, t2 )

      dist = 0

      do i = 1, n

        if ( ( t1(i) .eq. 0 .and. t2(i) .ne. 0 ) .or. 
     &       ( t1(i) .ne. 0 .and. t2(i) .eq. 0 ) ) then
          dist = dist + 1
        end if

      end do

      return
      end
      subroutine subset_enum ( n, nsub )

c*********************************************************************72
c
cc SUBSET_ENUM enumerates the subsets of a set with N elements.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of elements in the set.
c    N must be at least 0.
c
c    Output, integer NSUB, the number of distinct elements.
c
      implicit none

      integer n
      integer nsub

      nsub = 2**n

      return
      end
      subroutine subset_intersect ( n, a, b, c )

c*********************************************************************72
c
cc SUBSET_INTERSECT computes the intersection of two sets.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 July 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the order of the master set, of which A and
c    B are subsets.  N must be positive.
c
c    Input, integer A(N), B(N), two subsets of the master set.
c    A(I) = 0 if the I-th element is in the subset A, and is
c    1 otherwise; B is defined similarly.
c
c    Output, integer C(N), the intersection of A and B.
c
      implicit none

      integer n

      integer a(n)
      integer b(n)
      integer c(n)
      integer i
c
c  Check.
c
      call subset_check ( n, a )

      call subset_check ( n, b )

      do i = 1, n
        c(i) = min ( a(i), b(i) )
      end do

      return
      end
      subroutine subset_lex_rank ( n, t, rank )

c*********************************************************************72
c
cc SUBSET_LEX_RANK computes the lexicographic rank of a subset.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of items in the master set.
c    N must be positive.
c
c    Input, integer T(N), the subset.  If T(I) = 0, item I is
c    not in the subset; if T(I) = 1, item I is in the subset.
c
c    Output, integer RANK, the rank of the subset.
c
      implicit none

      integer n

      integer i
      integer rank
      integer t(n)
c
c  Check.
c
      call subset_check ( n, t )

      rank = 0

      do i = 1, n

        if ( t(i) .eq. 1 ) then
          rank = rank + 2**( n - i )
        end if

      end do

      return
      end
      subroutine subset_lex_successor ( n, t, rank )

c*********************************************************************72
c
cc SUBSET_LEX_SUCCESSOR computes the subset lexicographic successor.
c
c  Discussion:
c
c    In the original code, there is a last element with no successor.
c
c  Example:
c
c    On initial call, N is 5 and the input value of RANK is -1.
c    Then here are the successive outputs from the program:
c
c   Rank   T1   T2   T3   T4   T5
c   ----   --   --   --   --   --
c      0    0    0    0    0    0
c      1    0    0    0    0    1
c      2    0    0    0    1    0
c      3    0    0    0    1    1
c     ..   ..   ..   ..   ..   ..
c     30    1    1    1    1    0
c     31    1    1    1    1    1
c      0    0    0    0    0    0  <-- Cycle restarts with first element.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of elements in the master set.
c    N must be positive.
c
c    Input/output, integer T(N), describes a subset.  T(I) is 0 if
c    the I-th element of the master set is not in the subset, and is
c    1 if the I-th element is part of the subset.
c    On input, T describes a subset.
c    On output, T describes the next subset in the ordering.
c    If the input T was the last in the ordering, then the output T
c    will be the first.
c
c    Input/output, integer RANK, the rank.
c    If RANK = -1 on input, then the routine understands that this is
c    the first call, and that the user wishes the routine to supply
c    the first element in the ordering, which has RANK = 0.
c    In general, the input value of RANK is increased by 1 for output,
c    unless the very last element of the ordering was input, in which
c    case the output value of RANK is 0.
c
      implicit none

      integer n

      integer i
      integer rank
      integer t(n)
c
c  Return the first element.
c
      if ( rank .eq. -1 ) then
        do i = 1, n
          t(i) = 0
        end do
        rank = 0
        return
      end if
c
c  Check.
c
      call subset_check ( n, t )

      do i = n, 1, -1

        if ( t(i) .eq. 0 ) then
          t(i) = 1
          rank = rank + 1
          return
        else
          t(i) = 0
        end if

      end do

      rank = 0

      return
      end
      subroutine subset_lex_unrank ( rank, n, t )

c*********************************************************************72
c
cc SUBSET_LEX_UNRANK computes the subset of given lexicographic rank.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer RANK, the rank of the subset.
c
c    Input, integer N, the number of items in the master set.
c    N must be positive.
c
c    Output, integer T(N), the subset of the given rank.
c    If T(I) = 0, item I is not in the subset; if T(I) = 1, item I is in
c    the subset.
c
      implicit none

      integer n

      integer i
      integer nsub
      integer rank
      integer rank_copy
      integer t(n)
c
c  Check.
c
      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SUBSET_LEX_UNRANK - Fatal error!'
        write ( *, '(a)' ) '  Input N is illegal.'
        stop
      end if

      call subset_enum ( n, nsub )

      if ( rank .lt. 0 .or. nsub .lt. rank ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SUBSET_LEX_UNRANK - Fatal error!'
        write ( *, '(a)' ) '  The input rank is illegal.'
        stop
      end if

      rank_copy = rank

      do i = n, 1, -1

        if ( mod ( rank_copy, 2 ) .eq. 1 ) then
          t(i) = 1
        else
          t(i) = 0
        end if

        rank_copy = rank_copy / 2

      end do

      return
      end
      subroutine subset_union ( n, a, b, c )

c*********************************************************************72
c
cc SUBSET_UNION computes the union of two sets.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the order of the master set, of which A and
c    B are subsets.  N must be positive.
c
c    Input, integer A(N), B(N), two subsets of the master set.
c    A(I) = 0 if the I-th element is in the subset A, and is
c    1 otherwise; B is defined similarly.
c
c    Output, integer C(N), the union of A and B.
c
      implicit none

      integer n

      integer a(n)
      integer b(n)
      integer c(n)
      integer i
c
c  Check.
c
      call subset_check ( n, a )

      call subset_check ( n, b )

      do i = 1, n
        c(i) = max ( a(i), b(i) )
      end do

      return
      end
      subroutine subset_weight ( n, t, weight )

c*********************************************************************72
c
cc SUBSET_WEIGHT computes the Hamming weight of a set.
c
c  Discussion:
c
c    The Hamming weight is simply the number of elements in the set.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the order of the master set, of which T
c    is a subset.  N must be positive.
c
c    Input, integer T(N), defines the subset T.
c    T(I) is 1 if I is an element of T, and 0 otherwise.
c
c    Output, integer WEIGHT, the Hamming weight of the subset T.
c
      implicit none

      integer n

      integer i4vec_sum
      integer t(n)
      integer weight
c
c  Check.
c
      call subset_check ( n, t )

      weight = i4vec_sum ( n, t )

      return
      end
      subroutine subset_xor ( n, a, b, c )

c*********************************************************************72
c
cc SUBSET_XOR computes the symmetric difference of two sets.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the order of the master set, of which A and
c    B are subsets.  N must be positive.
c
c    Input, integer A(N), B(N), two subsets of the master set.
c    A(I) = 0 if the I-th element is in the subset A, and is
c    1 otherwise; B is defined similarly.
c
c    Output, integer C(N), the symmetric difference of A and B.
c
      implicit none

      integer n

      integer a(n)
      integer b(n)
      integer c(n)
      integer i
c
c  Check.
c
      call subset_check ( n, a )

      call subset_check ( n, b )

      do i = 1, n
        c(i) = max ( a(i), b(i) ) - min ( a(i), b(i) )
      end do

      return
      end
      subroutine subsetsum_swap ( n, a, sum_desired, index, 
     &  sum_achieved )

c*********************************************************************72
c
cc SUBSETSUM_SWAP seeks a solution of the subset sum problem by swapping.
c
c  Discussion:
c
c    Given a collection of N not necessarily distinct positive integers A(I),
c    and a positive integer SUM_DESIRED, select a subset of the values so that
c    their sum is as close as possible to SUM_DESIRED without exceeding it.
c
c  Algorithm:
c
c    Start with no values selected, and SUM_ACHIEVED = 0.
c
c    Consider each element A(I):
c
c      If A(I) is not selected and SUM_ACHIEVED + A(I) .le. SUM_DESIRED,
c        select A(I).
c
c      If A(I) is still not selected, and there is a selected A(J)
c      such that SUM_GOT .lt. SUM_ACHIEVED + A(I) - A(J),
c        select A(I) and deselect A(J).
c
c      If no items were selected on this sweep,
c        exit.
c      Otherwise,
c        repeat the search.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 March 2001
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of values.  N must be positive.
c
c    Input/output, integer A(N), a collection of positive values.
c    On output, A has been sorted into descending order.
c
c    Input, integer SUM_DESIRED, the desired sum.
c
c    Output, integer INDEX(N); INDEX(I) is 1 if A(I) is part of the
c    sum, and 0 otherwise.
c
c    Output, integer SUM_ACHIEVED, the sum of the selected
c    elements.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer index(n)
      integer j
      integer nmove
      integer sum_achieved
      integer sum_desired
c
c  Initialize.
c
      sum_achieved = 0

      do i = 1, n
        index(i) = 0
      end do
c
c  Sort into descending order.
c
      call i4vec_sort_insert_d ( n, a )

10    continue

        nmove = 0

        do i = 1, n

          if ( index(i) .eq. 0 ) then

            if ( sum_achieved + a(i) .le. sum_desired ) then
              index(i) = 1
              sum_achieved = sum_achieved + a(i)
              nmove = nmove + 1
              go to 30
            end if

          end if

          if ( index(i) .eq. 0 ) then

            do j = 1, n

              if ( index(j) .eq. 1 ) then

                if ( sum_achieved .lt. sum_achieved + a(i) - a(j) .and.
     &            sum_achieved + a(i) - a(j) .le. sum_desired ) then
                  index(j) = 0
                  index(i) = 1
                  nmove = nmove + 2
                  sum_achieved = sum_achieved + a(i) - a(j)
                  go to 20
                end if

              end if

            end do

20          continue

          end if

30        continue

        end do

        if ( nmove .le. 0 ) then
          go to 40
        end if

      go to 10

40    continue

      return
      end
      subroutine tableau_check ( n, tab )

c*********************************************************************72
c
cc TABLEAU_CHECK checks a 2 by N tableau.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of columns in the tableau.
c    N must be positive.
c
c    Input, integer TAB(2,N), a 2 by N tableau.
c
      implicit none

      integer n

      integer i
      integer j
      integer tab(2,n)

      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TABLEAU_CHECK - Fatal error!'
        write ( *, '(a)' ) '  N .lt. 1.'
        stop
        return
      end if
c
c  The entries must be between 0 and 2*N.
c
      do i = 1, 2
        do j = 1, n
          if ( tab(i,j) .lt. 1 .or. 2 * n .lt. tab(i,j) ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'TABLEAU_CHECK - Fatal error!'
            write ( *, '(a)' ) '  TAB(I,J) .lt. 1 or N .lt. TAB(I,J).'
            stop
          end if
        end do
      end do
c
c  The entries must be increasing to the right.
c
      do i = 1, 2
        do j = 2, n
          if ( tab(i,j) .le. tab(i,j-1) ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'TABLEAU_CHECK - Fatal error!'
            write ( *, '(a)' ) '  TAB(I,J) .le. TAB(I,J-1).'
            stop
          end if
        end do
      end do
c
c  The entries must be increasing down.
c
      i = 2
      do j = 1, n
        if ( tab(i,j) .le. tab(i-1,j) ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'TABLEAU_CHECK - Fatal error!'
          write ( *, '(a)' ) '  TAB(I,J) .le. TAB(I-1,J).'
          stop
        end if
      end do

      return
      end
      subroutine tableau_enum ( n, ntab )

c*********************************************************************72
c
cc TABLEAU_ENUM enumerates the 2 by N standard tableaus.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of columns in the tableau.
c    N must be nonnegative.
c
c    Output, integer NTAB, the number of 2 by N standard tableaus.
c
      implicit none

      integer binomial
      integer n
      integer ntab

      ntab = binomial ( 2 * n, n ) / ( n + 1 )

      return
      end
      subroutine tableau_to_bal_seq ( n, tab, t )

c*********************************************************************72
c
cc TABLEAU_TO_BAL_SEQ converts a 2 by N tableau to a balanced sequence.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of 0's (and 1's) in the sequence.
c    N must be positive.
c
c    Input, integer TAB(2,N), a 2 by N tableau.
c
c    Output, integer T(2*N), a balanced sequence.
c
      implicit none

      integer n

      integer i
      integer j
      integer t(2*n)
      integer tab(2,n)
c
c  Check.
c
      call tableau_check ( n, tab )

      do i = 1, 2
        do j = 1, n
          t(tab(i,j)) = i - 1
        end do
      end do

      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
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
      subroutine tree_check ( n, t, ierror )

c*********************************************************************72
c
cc TREE_CHECK checks a tree.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 February 2001
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of nodes in the tree.
c    N must be positive.
c
c    Input, integer T(2,N-1), describes the edges of the tree
c    as pairs of nodes.
c
c    Output, integer IERROR, error flag.
c    0, no error.
c    -1, N was illegal.
c    J, the edge T(1,J) to T(2,J) is illegal.
c
      implicit none

      integer n

      integer d(n)
      integer i
      integer ierror
      integer j
      integer k
      integer t(2,n-1)
      integer x
      integer y

      ierror = 0

      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TREE_CHECK - Fatal error!'
        write ( *, '(a)' ) '  N .lt. 1.'
        stop
      end if

      do i = 1, 2
        do j = 1, n - 1
          if ( t(i,j) .lt. 1 .or. n .lt. t(i,j) ) then
            ierror = j
            return
          end if
        end do
      end do
c
c  Compute the degree of each node.
c
      call edge_degree ( n, n - 1, t, d )
c
c  Delete a node of degree 1, N-1 times.
c
      do k = 1, n - 1

        x = 1

10      continue

        if ( d(x) .ne. 1 ) then
          x = x + 1
          if ( n .lt. x ) then
            ierror = -1
            return
          end if
          go to 10
        end if
c
c  Find its neighbor.
c
        j = 1

20      continue

          if ( t(1,j) .eq. x ) then
            y = t(2,j)
            go to 30
          end if

          if ( t(2,j) .eq. x ) then
            y = t(1,j)
            go to 30
          end if

          j = j + 1

          if ( n .lt. j ) then
            ierror = -1
            return
          end if

        go to 20

30      continue
c
c  Delete the edge.
c
        t(1,j) = - t(1,j)
        t(2,j) = - t(2,j)

        d(x) = d(x) - 1
        d(y) = d(y) - 1

      end do

      do j = 1, n - 1
        do i = 1, 2
          t(i,j) = - t(i,j)
        end do
      end do

      return
      end
      subroutine tree_enum ( n, ntree )

c*********************************************************************72
c
cc TREE_ENUM enumerates the trees on N nodes.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of nodes in each tree.
c    N must normally be at least 3, but for this routine,
c    any value of N is allowed.
c
c    Output, integer NTREE, the number of distinct elements.
c
      implicit none

      integer n
      integer ntree

      if ( n .lt. 1 ) then
        ntree = 0
      else if ( n .eq. 1 ) then
        ntree = 1
      else if ( n .eq. 2 ) then
        ntree = 1
      else
        ntree = n**( n - 2 )
      end if

      return
      end
      subroutine tree_rank ( n, t, rank )

c*********************************************************************72
c
cc TREE_RANK ranks a tree.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of nodes in the tree.
c    N must be at least 3.
c
c    Input, integer T(2,N-1), describes the edges of the tree
c    as pairs of nodes.
c
c    Output, integer RANK, the rank of the tree.
c
      implicit none

      integer n

      integer i
      integer ierror
      integer k
      integer p(n-2)
      integer rank
      integer t(2,n-1)
c
c  Check the tree.
c
      call tree_check ( n, t, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TREE_RANK - Fatal error!'
        write ( *, '(a)' ) '  Input tree is illegal.'
        write ( *, '(a,i8)' ) '  Error code = ', ierror
        stop
      end if
c
c  Convert the tree to a Pruefer code.
c
      call tree_to_pruefer ( n, t, p )
c
c  Find the rank of the Pruefer code.
c
      call pruefer_rank ( n, p, rank )

      return
      end
      subroutine tree_successor ( n, t, rank )

c*********************************************************************72
c
cc TREE_SUCCESSOR returns the successor of a tree.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of nodes in the tree.
c    N must be at least 3.
c
c    Input/output, integer T(2,N-1), describes the edges of the
c    tree as pairs of nodes.  On output, the input tree has been replaced
c    by its successor.
c
c    Input/output, integer RANK, the rank of the tree.
c
      implicit none

      integer n

      integer i
      integer ierror
      integer k
      integer p(n-2)
      integer rank
      integer t(2,n-1)
c
c  Return the first element.
c
      if ( rank .eq. -1 ) then
        do i = 1, n - 2
          p(i) = 1
        end do
        call pruefer_to_tree ( n, p, t )
        rank = 0
        return
      end if
c
c  Check the tree.
c
      call tree_check ( n, t, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TREE_SUCCESSOR - Fatal error!'
        write ( *, '(a)' ) '  Input tree is illegal.'
        write ( *, '(a,i8)' ) '  Error code = ', ierror
        stop
      end if
c
c  Convert the tree to a Pruefer code.
c
      call tree_to_pruefer ( n, t, p )
c
c  Find the successor of the Pruefer code.
c
      call pruefer_successor ( n, p, rank )
c
c  Convert the Pruefer code to the tree.
c
      call pruefer_to_tree ( n, p, t )

      return
      end
      subroutine tree_to_pruefer ( n, t, p )

c*********************************************************************72
c
cc TREE_TO_PRUEFER converts a tree to a Pruefer code.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 February 2001
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of nodes in the tree.
c    N must be positive.
c
c    Input, integer T(2,N-1), describes the edges of the tree
c    as pairs of nodes.
c
c    Output, integer P(N-2), the Pruefer code for the tree.
c
      implicit none

      integer n

      integer d(n)
      integer i
      integer ierror
      integer j
      integer k
      integer p(n-2)
      integer t(2,n-1)
      integer x
      integer y
c
c  Check.
c
      call tree_check ( n, t, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TREE_TO_PRUEFER - Fatal error!'
        write ( *, '(a)' ) '  The input array is illegal!'
        write ( *, '(a,i8)' ) '  IERROR = ', ierror
        stop
      end if
c
c  Compute the degree of each node.
c
      call edge_degree ( n, n - 1, t, d )

      do j = 1, n - 2
c
c  Find a node of degree 1.
c
        x = n

10      continue

        if ( d(x) .ne. 1 ) then
          x = x - 1
          go to 10
        end if
c
c  Find its neighbor.
c
        k = 1

20      continue

          if ( t(1,k) .eq. x ) then
            y = t(2,k)
            go to 30
          end if

          if ( t(2,k) .eq. x ) then
            y = t(1,k)
            go to 30
          end if

          k = k + 1

        go to 20

30      continue
c
c  Store the neighbor.
c
        p(j) = y
c
c  Delete the edge from the tree.
c
        d(x) = d(x) - 1
        d(y) = d(y) - 1

        t(1,k) = - t(1,k)
        t(2,k) = - t(2,k)

      end do
c
c  Remove the negative signs from the first N-2 columns of the tree.
c
      do j = 1, n - 2
        do i = 1, 2
          t(i,j) = - t(i,j)
        end do
      end do

      return
      end
      subroutine tree_unrank ( rank, n, t )

c*********************************************************************72
c
cc TREE_UNRANK unranks a tree.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer RANK, the rank of the tree.
c
c    Input, integer N, the number of nodes in the tree.
c    N must be at least 3.
c
c    Output, integer T(2,N-1), describes the edges of the tree
c    as pairs of nodes.
c
      implicit none

      integer n

      integer i
      integer p(n-2)
      integer rank
      integer rank_copy
      integer t(2,n-1)
      integer tree_num
c
c  Check.
c
      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TREE_UNRANK - Fatal error!'
        write ( *, '(a)' ) '  Input N is illegal.'
        stop
      end if

      call tree_enum ( n, tree_num )

      if ( rank .lt. 0 .or. tree_num .lt. rank ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TREE_UNRANK - Fatal error!'
        write ( *, '(a)' ) '  The input rank is illegal.'
        stop
      end if
c
c  Unrank the Pruefer code.
c
      call pruefer_unrank ( rank, n, p )
c
c  Convert the Pruefer code to a tree.
c
      call pruefer_to_tree ( n, p, t )

      return
      end
