      program main

c*********************************************************************72
c
cc MAIN is the main program for ASA205_PRB.
c
c  Discussion:
c
c    ASA205_PRB tests ASA205.
c
c  Modified:
c
c    29 November 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA205_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the ASA205 library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA205_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 examines a simple case with no repeated sum values.
c
c  Modified:
c
c    29 November 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 3 )
      integer n
      parameter ( n = 4 )

      external eval01
      integer ifault
      integer colsum(n)
      integer rowsum(m)

      data colsum / 2, 3, 1, 4 /
      data rowsum / 5, 3, 2 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' )
     &  '  The tables will not have any multiplicities.'

      call i4vec_print ( m, rowsum, '  The row sums:' )

      call i4vec_print ( n, colsum, '  The column sums:' )

      call enum ( m, n, rowsum, colsum, eval01, ifault )

      if ( ifault .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' )
     &  '  ENUM returned error flag IFAULT = ', ifault
      end if

      return
      end
      subroutine eval01 ( iflag, table, m, n, rowsum, colsum, prob,
     &  mult )

c*********************************************************************72
c
cc EVAL01 is called by ENUM every time a new contingency table is determined.
c
c  Modified:
c
c    29 November 2006
c
c  Author:
c
c    John Burkardt
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

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EVAL01'
        write ( *, '(a)' ) '  Only first ten cases will be printed.'
        write ( *, '(a)' ) ' '
c
c  Call with a new table.
c
      else if ( iflag .eq. 2 ) then

        count1 = count1 + 1
        count2 = count2 + mult
        psum = psum + mult * dexp ( prob )

        if ( count1 .le. 10 ) then

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'EVAL01:'
          write ( *, '(i3,i3,g14.6)' ) count1, mult, prob

          call i4mat_print ( m, n, table, '  Table' )

        end if
c
c  Last call, no table.
c
      else if ( iflag .eq. 3 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EVAL01 summary'
        write ( *, '(a,i8)' )
     &  '  Number of cases (ignoring multiplicity):', count1
        write ( *, '(a,i8)' )
     &  '  Number of cases (allowing multiplicity):', count2
        write ( *, '(a,g14.6)' ) '  Probability sum = ', psum

      end if

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 examines a case where a sum value is repeated.
c
c  Modified:
c
c    29 November 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 2 )
      integer n
      parameter ( n = 3 )

      external eval02
      integer ifault
      integer colsum(n)
      integer rowsum(m)

      data colsum / 1, 2, 1 /
      data rowsum / 3, 1 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  The data will have multiplicities.'

      call i4vec_print ( m, rowsum, '  The row sums:' )

      call i4vec_print ( n, colsum, '  The column sums:' )

      call enum ( m, n, rowsum, colsum, eval02, ifault )

      if ( ifault .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' )
     &  '  ENUM returned error flag IFAULT = ', ifault
      end if

      return
      end
      subroutine eval02 ( iflag, table, m, n, rowsum, colsum, prob,
     &  mult )

c*********************************************************************72
c
cc EVAL02 is called by ENUM every time a new contingency table is determined.
c
c  Modified:
c
c    29 November 2006
c
c  Author:
c
c    John Burkardt
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
c    Input, double precision PROB, the logarithm of the hypergeometric
c    probability of this table.
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
        psum = psum + mult * dexp ( prob )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EVAL02:'
        write ( *, '(i3,i3,g14.6)' ) count1, mult, prob

        call i4mat_print ( m, n, table, '  Table' )
c
c  Last call, no table.
c
      else if ( iflag .eq. 3 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EVAL02 summary'
        write ( *, '(a,i8)' )
     &  '  Number of cases (ignoring multiplicity):', count1
        write ( *, '(a,i8)' )
     &  '  Number of cases (allowing multiplicity):', count2
        write ( *, '(a,g14.6)' ) '  Probability sum = ', psum

      end if

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 examines a test case from the paper, with 118489 tables.
c
c  Modified:
c
c    29 November 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 4 )
      integer n
      parameter ( n = 3 )

      external eval03
      integer ifault
      integer colsum(n)
      integer rowsum(m)

      data colsum / 4, 57, 59 /
      data rowsum / 3, 38, 39, 40 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  Big problem test from the paper.'

      call i4vec_print ( m, rowsum, '  The row sums:' )

      call i4vec_print ( n, colsum, '  The column sums:' )

      call enum ( m, n, rowsum, colsum, eval03, ifault )

      if ( ifault .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' )
     &  '  ENUM returned error flag IFAULT = ', ifault
      end if

      return
      end
      subroutine eval03 ( iflag, table, m, n, rowsum, colsum, prob,
     &  mult )

c*********************************************************************72
c
cc EVAL03 is called by ENUM every time a new contingency table is determined.
c
c  Modified:
c
c    29 November 2006
c
c  Author:
c
c    John Burkardt
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
        psum = psum + mult * dexp ( prob )
c
c  Last call, no table.
c
      else if ( iflag .eq. 3 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'EVAL03 summary'
        write ( *, '(a,i8)' )
     &  '  Number of cases (ignoring multiplicity): ', count1
        write ( *, '(a,i8)' )
     &  '  Number of cases (allowing multiplicity): ', count2
        write ( *, '(a,g14.6)' ) '  Probability sum = ', psum
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Result from paper:'
        write ( *, '(a,i8)' )
     &  '  Number of cases (ignoring multiplicity): ', 118489

      end if

      return
      end
