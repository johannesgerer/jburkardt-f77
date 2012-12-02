      subroutine rcont ( nrow, ncol, nrowt, ncolt, nsubt, matrix, key,
     &  ifault )

c*********************************************************************72
c
cc RCONT generates a random two-way table with given marginal totals.
c
c  Discussion:
c
c    Each time the program is called, another table will be randomly
c    generated.
c
c    Note that it should be the case that the sum of the row totals
c    is equal to the sum of the column totals.  However, this program
c    does not check for that condition.
c
c  Modified:
c
c    06 April 2006
c
c  Author:
c
c    Original FORTRAN77 version by James Boyett.
c    Modifications by John Burkardt.
c
c  Reference:
c
c    James Boyett,
c    Algorithm AS 144: 
c    Random R x C Tables with Given Row and Column Totals,
c    Applied Statistics,
c    Volume 28, Number 3, pages 329-332, 1979.
c
c  Parameters:
c
c    Input, integer NROW, the number of rows in the observed matrix.
c
c    Input, integer NCOL, the number of columns in the observed matrix.
c
c    Input, integer NROWT(NROW), the row totals of the observed matrix.
c
c    Input, integer NCOLT(NCOL), the column totals of the observed matrix.
c
c    Input/output, integer NSUBT(NCOL), used by RCONT for partial column sums.
c    Must not be changed by the calling program.
c
c    Output, integer MATRIX(NROW,NCOL), the randomly generated matrix.
c
c    Input/output, logical KEY, should be set to FALSE by the user before
c    the initial call.  RCONT will reset it to TRUE, and it should be left
c    at that value for subsequent calls in which the same values of NROW,
c    NCOL, NROWT and NCOLT are being used.
c
c    Output, integer IFAULT, fault indicator.
c    0, no error occured.
c    1, NROW <= 0.
c    2, NCOL <= 1.
c    3, some entry of NROWT is less than 0.
c    4, some entry of NCOLT is less than 0.
c    5, the sample size (sum of the column totals) is too large.
c 
      implicit none

      integer ncol
      integer nrow

      integer nvec_max
      parameter ( nvec_max = 200 )

      integer i
      integer ifault
      integer ii
      integer j
      integer k
      logical key
      integer limit
      integer matrix(nrow,ncol)
      integer ncolt(ncol)
      integer nnvect(nvec_max)
      integer noct
      integer nrowt(nrow)
      integer nsubt(ncol)
      integer ntemp
      integer ntotal
      integer nvect(nvec_max)
      double precision r8_uniform_01
      integer seed

      save ntotal
      save nvect
      save seed

      ifault = 0

      if ( .not. key ) then
c
c  Set KEY for subsequent calls.
c
        key = .true.
        seed = 123456789
c
c  Check for faults and prepare for future calls.
c
        if ( nrow .le. 0 ) then
          ifault = 1
          return
        end if

        if ( ncol .le. 1 ) then
          ifault = 2
          return
        end if

        do i = 1, nrow
          if ( nrowt(i) .le. 0 ) then
            ifault = 3
            return
          end if
        end do

        if ( ncolt(1) .le. 0 ) then
          ifault = 4
          return
        end if

        nsubt(1) = ncolt(1)

        do j = 2, ncol

          if ( ncolt(j) .le. 0 ) then
            ifault = 4
            return
          end if

          nsubt(j) = nsubt(j-1) + ncolt(j)

        end do

        ntotal = nsubt(ncol)

        if ( nvec_max .lt. ntotal ) then
          ifault = 5
          return
        end if
c
c  Initialize vector to be permuted.
c
        do i = 1, ntotal
          nvect(i) = i
        end do

      end if
c
c  Initialize vector to be permuted.
c
      do i = 1, ntotal
        nnvect(i) = nvect(i)
      end do
c
c  Permute vector.
c
      ntemp = ntotal
      do i = 1, ntotal
        noct = int ( r8_uniform_01 ( seed ) * dble ( ntemp ) + 1.0D+00 )
        nvect(i) = nnvect(noct)
        nnvect(noct) = nnvect(ntemp)
        ntemp = ntemp - 1
      end do
c
c  Construct random matrix.
c
      do i = 1, nrow
        do j = 1, ncol
          matrix(i,j) = 0
        end do
      end do

      ii = 1

      do i = 1, nrow

        limit = nrowt(i)

        do k = 1, limit

          do j = 1, ncol
            if ( nvect(ii) .le. nsubt(j) ) then
              ii = ii + 1
              matrix(i,j) = matrix(i,j) + 1
              go to 208
            end if
          end do

208       continue

        end do

      end do

      return
      end
      subroutine i4mat_print ( m, n, a, title )

c*********************************************************************72
c
cc I4MAT_PRINT prints an I4MAT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 May 2000
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
c    Input, character * ( * ) TITLE, a title to be printed first.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer i
      integer j
      integer jhi
      integer jlo
      character * ( * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) title

      do jlo = 1, n, 10
        jhi = min ( jlo + 9, n )
        write ( *, '(a)' ) ' '
        write ( *, '(6x,10(i7))' ) ( j, j = jlo, jhi )
        write ( *, * ) ' '
        do i = 1, m
          write ( *, '(i6,10i7)' ) i, ( a(i,j), j = jlo, jhi )
        end do
      end do

      return
      end
      subroutine i4vec_print ( n, a, title )

c*********************************************************************72
c
cc I4VEC_PRINT prints an I4VEC.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 November 2000
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
c    Input, character * ( * ) TITLE, a title to be printed first.
c    TITLE may be blank.
c
      implicit none

      integer n

      integer a(n)
      integer big
      integer i
      character * ( * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) title

      big = 0.0
      do i = 1, n
        big = max ( big, abs ( a(i) ) )
      end do

      write ( *, '(a)' ) ' '
      if ( big < 1000 ) then
        do i = 1, n
          write ( *, '(i6,1x,i4)' ) i, a(i)
        end do
      else if ( big < 1000000 ) then
        do i = 1, n
          write ( *, '(i6,1x,i7)' ) i, a(i)
        end do
      else
        do i = 1, n
          write ( *, '(i6,i11)' ) i, a(i)
        end do
      end if

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
