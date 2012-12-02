      subroutine i4mat_print ( m, n, a, title )

c*********************************************************************72
c
cc I4MAT_PRINT prints an I4MAT.
c
c  Discussion:
c
c    An I4MAT is a rectangular array of integer values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
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

      integer m
      integer n

      integer a(m,n)
      integer ihi
      integer ilo
      integer jhi
      integer jlo
      character * ( * ) title

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
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
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

      integer incx
      parameter ( incx = 10 )
      integer m
      integer n

      integer a(m,n)
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
      function r8_uniform_01 ( seed )

c*********************************************************************72
c
cc R8_UNIFORM_01 returns a unit pseudorandom R8.
c
c  Discussion:
c
c    This routine implements the recursion
c
c      seed = 16807 * seed mod ( 2**31 - 1 )
c      r8_uniform_01 = seed / ( 2**31 - 1 )
c
c    The integer arithmetic never requires more than 32 bits,
c    including a sign bit.
c
c    If the initial seed is 12345, then the first three computations are
c
c      Input     Output      R8_UNIFORM_01
c      SEED      SEED
c
c         12345   207482415  0.096616
c     207482415  1790989824  0.833995
c    1790989824  2035175616  0.947702
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 August 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Second Edition,
c    Springer, 1987,
c    ISBN: 0387964673,
c    LC: QA76.9.C65.B73.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, December 1986, pages 362-376.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley, 1998,
c    ISBN: 0471134031,
c    LC: T57.62.H37.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, Number 2, 1969, pages 136-143.
c
c  Parameters:
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R8_UNIFORM_01, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer k
      double precision r8_uniform_01
      integer seed

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + i4_huge
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r8_uniform_01 = dble ( seed ) * 4.656612875D-10

      return
      end
      subroutine rcont2 ( nrow, ncol, nrowt, ncolt, jwork, key, seed, 
     &  matrix, ifault )

c*********************************************************************72
c
cc RCONT2 generates a random two-way table with given marginal totals.
c
c  Discussion:
c
c    See the journal article for the significance of the lines
c    which start with c*.
c
c  Modified:
c
c    10 March 2009
c
c  Author:
c
c    WM Patefield
c
c  Reference:
c
c    WM Patefield,
c    Algorithm AS 159:
c    An Efficient Method of Generating RXC Tables with
c    Given Row and Column Totals,
c    Applied Statistics,
c    Volume 30, Number 1, 1981, pages 91-97.
c
c  Parameters:
c
c    Input, integer NROW, NCOL, the number of rows and columns in the table.
c    NROW and NCOL must each be at least 2.
c
c    Input, integer NROWT(NROW), NCOLT(NCOL), the row and column sums.
c    Each entry must be positive.
c
c    Workspace, integer JWORK(NCOL).
c
c    Input/output, logical KEY, a flag that indicates whether data has
c    been initialized for this problem.  Set KEY = .FALSE. before the first
c    call.
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, integer MATRIX(NROW,NCOL), the matrix.
c
c    Output, integer IFAULT, an error flag, which is returned as 0
c    if no error occurred.
c
      implicit none

      integer ncol
      integer nrow

      double precision dummy
      double precision fact(5001)
      integer i
      integer ia
      integer iap
      integer ib
      integer ic
      integer id
      integer idp
      integer ie
      integer ifault
      integer igp
      integer ihp
      integer ii
      integer iip
      integer j
      integer jc
      integer jwork(ncol)
      logical key
      integer l
      logical lsm
      logical lsp
      integer m
      integer matrix(nrow,ncol)
      integer maxtot
      parameter ( maxtot = 5000 )
      integer ncolm
      integer ncolt(ncol)
      integer nll
      integer nlm
      integer nlmp
      integer nrowm
      integer nrowt(nrow)
      integer nrowtl
      integer ntotal
      double precision r8_uniform_01
      integer seed
      double precision sumprb
      double precision x
      double precision y

      common / b / ntotal, nrowm, ncolm, fact
c
c*      common /tempry/ hop
c
      ifault = 0
      if (key) go to 103
c
c  Set KEY for subsequent calls.
c
      key = .true.
c
c  Check for faults and prepare for future calls.
c
      if ( nrow .le. 1 ) then
        ifault = 1
        return
      end if

      if (ncol .le. 1) go to 213
      nrowm = nrow - 1
      ncolm = ncol - 1
      do i = 1, nrow
        if (nrowt(i) .le. 0) go to 214
      end do

      ntotal = 0
      do j = 1, ncol
        if (ncolt(j) .le. 0) go to 215
        ntotal = ntotal + ncolt(j)
      end do

      if (ntotal .gt. maxtot) go to 216
c
c  Calculate log-factorials
c
      x = 0.0D+00
      fact(1) = 0.0D+00
      do i = 1, ntotal
        x = x + log ( dble ( i ) )
        fact(i+1) = x
      end do
c
c  Construct random matrix
c
  103 continue

      do j = 1, ncolm
        jwork(j) = ncolt(j)
      end do
      jc = ntotal
c
c*      hop = 1.0D+00
c
      do 190 l = 1, nrowm
        nrowtl = nrowt(l)
        ia = nrowtl
        ic = jc
        jc = jc - nrowtl
        do 180 m = 1, ncolm
          id = jwork(m)
          ie = ic
          ic = ic - id
          ib = ie - ia
          ii = ib - id
c
c  Test for zero entries in matrix
c
          if (ie .ne. 0) go to 130

          do j = m, ncol
            matrix(l,j) = 0
          end do

          go to 190
c
c     Generate pseudo-random number
c
  130      dummy = r8_uniform_01 ( seed )
c
c     compute conditional expected value of matrix(l, m)
c
  131     nlm = int ( dble ( ia * id ) / dble ( ie ) + 0.5D+00 )
          iap = ia + 1
          idp = id + 1
          igp = idp - nlm
          ihp = iap - nlm
      nlmp = nlm + 1
      iip = ii + nlmp
          x = exp(fact(iap) + fact(ib+1) + fact(ic+1) + fact(idp) -
     *      fact(ie+1) - fact(nlmp) - fact(igp) - fact(ihp) - fact(iip))
          if (x .ge. dummy) go to 160
      sumprb = x
      y = x
      nll = nlm
      lsp = .false.
      lsm = .false.
c
c  Increment entry in row l, column m.
c
  140 continue

      j = (id - nlm) * (ia - nlm)
      if (j .eq. 0) go to 156
      nlm = nlm + 1
      x = x * dble ( j ) / dble ( nlm * (ii + nlm))
      sumprb = sumprb + x
      if (sumprb .ge. dummy) go to 160

  150 continue

      if (lsm) go to 155
c
c     decrement entry in row l, column m
c
      j = nll * (ii + nll)
      if (j .eq. 0) go to 154
      nll = nll - 1
      y = y * dble ( j ) / dble ( ( id - nll ) * ( ia - nll ) )
      sumprb = sumprb + y
      if (sumprb .ge. dummy) go to 159
      if (.not. lsp) go to 140
          go to 150
  154     lsm = .true.
  155     if (.not. lsp) go to 140
          dummy = sumprb * r8_uniform_01 ( seed )
          go to 131
  156     lsp = .true.
          go to 150
  159     nlm = nll
c
c*          hop = hop * y
c*          go to 161
c*160       hop = hop * x
c*161       matrix(l,m) = nlm
c
  160     matrix(l,m) = nlm
          ia = ia - nlm
          jwork(m) = jwork(m) - nlm
  180   continue
        matrix(l,ncol) = ia
  190 continue
c
c  compute entries in last row of matrix
c
      do m = 1, ncolm
        matrix(nrow,m) = jwork(m)
      end do

      matrix(nrow,ncol) = ib - matrix(nrow,ncolm)
      return
c
c  set faults
c
  213 ifault = 2
      return
  214 ifault = 3
      return
  215 ifault = 4
      return
  216 ifault = 5
      return
      end
      function s_len_trim ( s )

c*********************************************************************72
c
cc S_LEN_TRIM returns the length of a string to the last nonblank.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
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
