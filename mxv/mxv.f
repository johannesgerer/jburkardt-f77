      program main

c*********************************************************************72
c
cc MAIN is the main program for MXV.
c
c  Discussion:
c
c    MXV computes a matrix-vector product in a number of ways, and reports
c    the elapsed CPU time.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Usage:
c
c    mxv m n
c
c  Parameters:
c
c    Command line argument, integer M, the number of rows in the matrix.
c
c    Command line argument, integer N, the number of columns in the matrix.
c
      implicit none

      integer maxm
      parameter ( maxm = 100000 )
      integer maxn
      parameter ( maxn = 100000 )
      integer maxmn
      parameter ( maxmn = 1000000 )

      double precision a(maxmn)
      integer arg_num
      double precision cpu_seconds
      integer flop_count
      integer iarg
      integer iargc
      integer ierror
      integer last
      integer m
      double precision mflops
      integer n
      character * ( 80 ) string
      double precision x(maxn)
      double precision y(maxm)

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MXV:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Compute matrix vector products y = A*x.'
c
c  Get the number of command line arguments.
c
      arg_num = iargc ( )
c
c  Get the number of rows, M.
c
      if ( 1 .le. arg_num ) then
      
        iarg = 1
        call getarg ( iarg, string )
        call s_to_i4 ( string, m, ierror, last )
        
      else
      
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Enter the number of rows, M'
        read ( *, * ) m
        
      end if
c
c  Get the number of columns, N.
c
      if ( 2 .le. arg_num ) then
      
        iarg = 2
        call getarg ( iarg, string )
        call s_to_i4 ( string, n, ierror, last )
        
      else
      
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Enter the number of rows, N'
        read ( *, * ) n
        
      end if
c
c  Reject dimensions that are too large.
c
      if ( maxm .lt. m ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MXV - Fatal error!'
        write ( *, '(a)' ) '  MaxM < M'
        write ( *, '(a,i8)' ) '  MaxM = ', maxm
        write ( *, '(a,i8)' ) '  M =    ', m
        stop
      end if

      if ( maxn .lt. n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MXV - Fatal error!'
        write ( *, '(a)' ) '  MaxN < N'
        write ( *, '(a,i8)' ) '  MaxN = ', maxn
        write ( *, '(a,i8)' ) '  N =    ', n
        stop
      end if

      if ( maxmn .lt. m * n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MXV - Fatal error!'
        write ( *, '(a)' ) '  MaxMN < M * N'
        write ( *, '(a,i8)' ) '  MaxMN = ', maxmn
        write ( *, '(a,i8)' ) '  M * N = ', m * n
        stop
      end if
c
c  Record the amount of work.
c  Each of the M entries of Y requires N multiplies and N adds.
c
      flop_count = 2 * m * n

      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)'  ) 
     &  '  Number of matrix rows M =             ', m
      write ( *, '(a,i12)'  ) 
     &  '  Number of matrix columns N =          ', n
      write ( *, '(a,i12)'  ) 
     &  '  Number of floating point operations = ', flop_count
c
c  Set A and X.
c
      call matgen ( m, n, a, x )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Method     Cpu Seconds       MegaFlopS'
      write ( *, '(a)' ) '  ------  --------------  --------------'
c
c  DOIDOJ
c
      call mxv_doidoj ( m, n, a, x, y, cpu_seconds )

      if ( 0.0D+00 .lt. cpu_seconds ) then
        mflops = dble ( flop_count ) / cpu_seconds / 1000000.0D+00
      else
        mflops = -1.0D+00
      end if

      write ( *, '(2x,a6,2x,g14.6,2x,g14.6)' ) 
     &  'DOIDOJ', cpu_seconds, mflops
c
c  DOJDOI
c
      call mxv_dojdoi ( m, n, a, x, y, cpu_seconds )

      if ( 0.0D+00 .lt. cpu_seconds ) then
        mflops = dble ( flop_count ) / cpu_seconds / 1000000.0D+00
      else
        mflops = -1.0D+00
      end if

      write ( *, '(2x,a6,2x,g14.6,2x,g14.6)' ) 
     &  'DOJDOI', cpu_seconds, mflops
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MXV:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine matgen ( m, n, a, x )

c*********************************************************************72
c
cc MATGEN generates a random matrix A and vector X.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns 
c    of the matrix.
c
c    Output, double precision A(M,N), the matrix.
c
c    Output, double precision X(N), the vector.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      integer i
      integer j
      integer seed
      double precision x(n)

      seed = 1325
c
c  Set the matrix A.
c
      do j = 1, n
        do i = 1, m
          seed = mod ( ( 3125 * seed ), 65536 )
          a(i,j) = ( seed - 32768.0 ) / 16384.0
        end do
      end do
c
c  Set X.
c
      do i = 1, n
        x(i) = i
      end do

      return
      end
      subroutine mxv_doidoj ( m, n, a, x, y, cpu_seconds )

c*********************************************************************72
c
cc MXV_DOIDOJ computes y = A * x, using DO I, DO J loops.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns
c    of the matrix.
c
c    Input, double precision A(M,N), the matrix.
c
c    Input, double precision X(N), the vector to be multiplied.
c
c    Output, double precision Y(M), the product vector.
c
c    Output, double precision CPU_SECONDS, the elapsed CPU time.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision cpu_seconds
      integer i
      integer j
      double precision time1
      double precision time2
      double precision x(n)
      double precision y(m)

      call cpu_time ( time1 )

      do i = 1, m
        y(i) = 0.0
        do j = 1, n
          y(i) = y(i) + a(i,j) * x(j)
        end do
      end do

      call cpu_time ( time2 )

      cpu_seconds = time2 - time1

      return
      end
      subroutine mxv_dojdoi ( m, n, a, x, y, cpu_seconds )

c*********************************************************************72
c
cc MXV_DOJDOI computes y = A * x, using DO J, DO I loops.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns
c    of the matrix.
c
c    Input, double precision A(M,N), the matrix.
c
c    Input, double precision X(N), the vector to be multiplied.
c
c    Output, double precision Y(M), the product vector.
c
c    Output, double precision CPU_SECONDS, the elapsed CPU time.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision cpu_seconds
      integer i
      integer j
      double precision time1
      double precision time2
      double precision x(n)
      double precision y(m)

      call cpu_time ( time1 )

      do i = 1, m
        y(i) = 0.0
      end do

      do j = 1, n
        do i = 1, m
          y(i) = y(i) + a(i,j) * x(j)
        end do
      end do

      call cpu_time ( time2 )

      cpu_seconds = time2 - time1

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
      subroutine s_to_i4 ( s, ival, ierror, length )

c*********************************************************************72
c
cc S_TO_I4 reads an I4 from a string.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) S, a string to be examined.
c
c    Output, integer IVAL, the integer value read from the string.
c    If the string is blank, then IVAL will be returned 0.
c
c    Output, integer IERROR, an error flag.
c    0, no error.
c    1, an error occurred.
c
c    Output, integer LENGTH, the number of characters of S
c    used to make IVAL.
c
      implicit none

      character c
      integer i
      integer ierror
      integer isgn
      integer istate
      integer ival
      integer length
      character * ( * ) s
      integer s_len_trim

      ierror = 0
      istate = 0
      isgn = 1
      ival = 0

      do i = 1, s_len_trim ( s )

        c = s(i:i)
c
c  Haven't read anything.
c
        if ( istate .eq. 0 ) then

          if ( c .eq. ' ' ) then

          else if ( c .eq. '-' ) then
            istate = 1
            isgn = -1
          else if ( c .eq. '+' ) then
            istate = 1
            isgn = + 1
          else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
            istate = 2
            ival = ichar ( c ) - ichar ( '0' )
          else
            ierror = 1
            return
          end if
c
c  Have read the sign, expecting digits.
c
        else if ( istate .eq. 1 ) then

          if ( c .eq. ' ' ) then

          else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
            istate = 2
            ival = ichar ( c ) - ichar ( '0' )
          else
            ierror = 1
            return
          end if
c
c  Have read at least one digit, expecting more.
c
        else if ( istate .eq. 2 ) then

          if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
            ival = 10 * ival + ichar ( c ) - ichar ( '0' )
          else
            ival = isgn * ival
            length = i - 1
            return
          end if

        end if

      end do
c
c  If we read all the characters in the string, see if we're OK.
c
      if ( istate .eq. 2 ) then
        ival = isgn * ival
        length = s_len_trim ( s )
      else
        ierror = 1
        length = 0
      end if

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
