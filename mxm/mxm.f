      program main

c*********************************************************************72
c
cc MAIN is the main program for MXM.
c
c  Discussion:
c
c    MXV computes a matrix-matrix product in a number of ways, and reports
c    the elapsed CPU time.
c
c    The multiplication carried out is
c
c      A(1:N1,1:N3) = B(1:N1,1:N2) * C(1:N2,1:N3)
c
c    where B and C are real double precision matrices whose entries
c    are assigned randomly.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 October 2010
c
c  Author:
c
c    John Burkardt
c
c  Usage:
c
c    mxm n1 n2 n3
c
c  Parameters:
c
c    Command line argument, integer N1, N2, N3, defines the number of
c    rows and columns in the two matrices.
c
      implicit none

      integer arg_num
      integer iarg
      integer iargc
      integer ierror
      integer last
      integer n1
      integer n2
      integer n3
      character * ( 80 ) string

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MXM:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Compute matrix-matrix product A = B * C'
c
c  Get the number of command line arguments.
c
      arg_num = iargc ( )
c
c  Get N1.
c
      if ( 1 .le. arg_num ) then
      
        iarg = 1
        call getarg ( iarg, string )
        call s_to_i4 ( string, n1, ierror, last )
        
      else
      
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Enter N1, the number of rows in B.'
        read ( *, * ) n1
        
      end if
c
c  Get N2.
c
      if ( 2 .le. arg_num ) then
      
        iarg = 2
        call getarg ( iarg, string )
        call s_to_i4 ( string, n2, ierror, last )
        
      else
      
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &    '  Enter N2, the number of columns in B and rows in C.'
        read ( *, * ) n2
        
      end if
c
c  Get N3.
c
      if ( 3 .le. arg_num ) then
      
        iarg = 3
        call getarg ( iarg, string )
        call s_to_i4 ( string, n3, ierror, last )
        
      else
      
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Enter N3, the number of columns in C.'
        read ( *, * ) n3
        
      end if
c
c  Try to get around FORTRAN77's inability to allocate memory dynamically
c  by calling a subroutine to create the arrays.
c
      call mxm_sub ( n1, n2, n3 )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MXM:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine mxm_sub ( n1, n2, n3 )

c*********************************************************************72
c
cc MXM_SUB carries out the computations.
c
c  Discussion:
c
c    We would like to use matrices whose dimensions are not specified
c    until run time.  FORTRAN77 cannot allocate memory in this way,
c    but at least some compilers will allow you to call a subroutine
c    inside of which are "automatic" arrays whose memory is allocated
c    on call, based on an input argument.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 October 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, N3, define the orders of the
c    matrices.
c
      double precision a(n1,n3)
      double precision b(n1,n2)
      double precision c(n2,n3)
      double precision cpu_seconds
      integer flop_count
      double precision mflops
      integer n1
      integer n2
      integer n3
      integer seed
      double precision time_estimate
c
c  Record the amount of work.
c  Each of the N1 * N3 entries of A requires N2 multiplies and N2 adds.
c
      flop_count = 2 * n1 * n2 * n3

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6,a,i6)' ) '  Matrix B is ', n1, ' by ', n2
      write ( *, '(a,i6,a,i6)' ) '  Matrix C is ', n2, ' by ', n3
      write ( *, '(a,i6,a,i6)' ) '  Matrix A will be ', n1, ' by ', n3
      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) 
     &  '  Number of floating point operations = ', flop_count
      time_estimate = dble ( flop_count ) / 2.6E+09
      write ( *, '(a,g10.2,a)' ) 
     &  '  Estimated CPU time is ', time_estimate, ' seconds.'
c
c  Set B and C.
c
      seed = 1325
      call matgen ( n1, n2, seed, b )
      call matgen ( n2, n3, seed, c )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Method     Cpu Seconds       MegaFlopS'
      write ( *, '(a)' ) '  ------  --------------  --------------'
c
c  IJK
c
      a(1:n1,1:n3) = 0.0D+00
      call mxm_ijk ( n1, n2, n3, a, b, c, cpu_seconds )

      if ( 0.0D+00 .lt. cpu_seconds ) then
        mflops = dble ( flop_count ) / cpu_seconds / 1000000.0D+00
      else
        mflops = -1.0D+00
      end if

      write ( *, '(2x,a6,2x,g14.6,2x,g14.6)' ) 'IJK   ', cpu_seconds, mflops
c
c  IKJ
c
      a(1:n1,1:n3) = 0.0D+00
      call mxm_ikj ( n1, n2, n3, a, b, c, cpu_seconds )

      if ( 0.0D+00 .lt. cpu_seconds ) then
        mflops = dble ( flop_count ) / cpu_seconds / 1000000.0D+00
      else
        mflops = -1.0D+00
      end if

      write ( *, '(2x,a6,2x,g14.6,2x,g14.6)' ) 'IKJ   ', cpu_seconds, mflops
c
c  JIK
c
      a(1:n1,1:n3) = 0.0D+00
      call mxm_jik ( n1, n2, n3, a, b, c, cpu_seconds )

      if ( 0.0D+00 .lt. cpu_seconds ) then
        mflops = dble ( flop_count ) / cpu_seconds / 1000000.0D+00
      else
        mflops = -1.0D+00
      end if

      write ( *, '(2x,a6,2x,g14.6,2x,g14.6)' ) 'JIK   ', cpu_seconds, mflops
c
c  JKI
c
      a(1:n1,1:n3) = 0.0D+00
      call mxm_jki ( n1, n2, n3, a, b, c, cpu_seconds )

      if ( 0.0D+00 .lt. cpu_seconds ) then
        mflops = dble ( flop_count ) / cpu_seconds / 1000000.0D+00
      else
        mflops = -1.0D+00
      end if

      write ( *, '(2x,a6,2x,g14.6,2x,g14.6)' ) 'JKI   ', cpu_seconds, mflops
c
c  KIJ
c
      a(1:n1,1:n3) = 0.0D+00
      call mxm_kij ( n1, n2, n3, a, b, c, cpu_seconds )

      if ( 0.0D+00 .lt. cpu_seconds ) then
        mflops = dble ( flop_count ) / cpu_seconds / 1000000.0D+00
      else
        mflops = -1.0D+00
      end if

      write ( *, '(2x,a6,2x,g14.6,2x,g14.6)' ) 'KIJ   ', cpu_seconds, mflops
c
c  KJI
c
      a(1:n1,1:n3) = 0.0D+00
      call mxm_kji ( n1, n2, n3, a, b, c, cpu_seconds )

      if ( 0.0D+00 .lt. cpu_seconds ) then
        mflops = dble ( flop_count ) / cpu_seconds / 1000000.0D+00
      else
        mflops = -1.0D+00
      end if

      write ( *, '(2x,a6,2x,g14.6,2x,g14.6)' ) 'KJI   ', cpu_seconds, mflops

      return
      end
      subroutine matgen ( m, n, seed, a )

c*********************************************************************72
c
cc MATGEN generates a random matrix.
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
c    Input, integer SEED, a seed for the random number
c    generator.
c
c    Output, double precision A(M,N), the matrix.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      integer i
      integer j
      integer seed
c
c  Set the matrix A.
c
      do j = 1, n
        do i = 1, m
          seed = mod ( ( 3125 * seed ), 65536 )
          a(i,j) = ( seed - 32768.0 ) / 16384.0
        end do
      end do

      return
      end
      subroutine mxm_ijk ( n1, n2, n3, a, b, c, cpu_seconds )

c*********************************************************************72
c
cc MXM_IJK computes A = B * C using DO I, DO J, DO K loops.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 October 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, N3, define the orders of the
c    matrices.
c
c    Output, double precision A(N1,N3), the result matrix.
c
c    Input, double precision B(N1,N2), C(N2,N3), the factor matrices.
c
c    Output, double precision CPU_SECONDS, the elapsed CPU time.
c
      implicit none

      integer n1
      integer n2
      integer n3

      double precision a(n1,n3)
      double precision b(n1,n2)
      double precision c(n2,n3)
      double precision cpu_seconds
      integer i
      integer j
      integer k
      double precision time1
      double precision time2

      call cpu_time ( time1 )

      do i = 1, n1
        do j = 1, n3
          do k = 1, n2
            a(i,j) = a(i,j) + b(i,k) * c(k,j)
          end do
        end do
      end do

      call cpu_time ( time2 )

      cpu_seconds = time2 - time1

      return
      end
      subroutine mxm_ikj ( n1, n2, n3, a, b, c, cpu_seconds )

c*********************************************************************72
c
cc MXM_IKJ computes A = B * C using DO I, DO K, DO J loops.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 October 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, N3, define the orders of the
c    matrices.
c
c    Output, double precision A(N1,N3), the result matrix.
c
c    Input, double precision B(N1,N2), C(N2,N3), the factor matrices.
c
c    Output, double precision CPU_SECONDS, the elapsed CPU time.
c
      implicit none

      integer n1
      integer n2
      integer n3

      double precision a(n1,n3)
      double precision b(n1,n2)
      double precision c(n2,n3)
      double precision cpu_seconds
      integer i
      integer j
      integer k
      double precision time1
      double precision time2

      call cpu_time ( time1 )

      do i = 1, n1
        do k = 1, n2
          do j = 1, n3
            a(i,j) = a(i,j) + b(i,k) * c(k,j)
          end do
        end do
      end do

      call cpu_time ( time2 )

      cpu_seconds = time2 - time1

      return
      end
      subroutine mxm_jik ( n1, n2, n3, a, b, c, cpu_seconds )

c*********************************************************************72
c
cc MXM_JIK computes A = B * C using DO J, DO I, DO K loops.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 October 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, N3, define the orders of the
c    matrices.
c
c    Output, double precision A(N1,N3), the result matrix.
c
c    Input, double precision B(N1,N2), C(N2,N3), the factor matrices.
c
c    Output, double precision CPU_SECONDS, the elapsed CPU time.
c
      implicit none

      integer n1
      integer n2
      integer n3

      double precision a(n1,n3)
      double precision b(n1,n2)
      double precision c(n2,n3)
      double precision cpu_seconds
      integer i
      integer j
      integer k
      double precision time1
      double precision time2

      call cpu_time ( time1 )

      do j = 1, n3
        do i = 1, n1
          do k = 1, n2
            a(i,j) = a(i,j) + b(i,k) * c(k,j)
          end do
        end do
      end do

      call cpu_time ( time2 )

      cpu_seconds = time2 - time1

      return
      end
      subroutine mxm_jki ( n1, n2, n3, a, b, c, cpu_seconds )

c*********************************************************************72
c
cc MXM_JKI computes A = B * C using DO J, DO K, DO I loops.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 October 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, N3, define the orders of the
c    matrices.
c
c    Output, double precision A(N1,N3), the result matrix.
c
c    Input, double precision B(N1,N2), C(N2,N3), the factor matrices.
c
c    Output, double precision CPU_SECONDS, the elapsed CPU time.
c
      implicit none

      integer n1
      integer n2
      integer n3

      double precision a(n1,n3)
      double precision b(n1,n2)
      double precision c(n2,n3)
      double precision cpu_seconds
      integer i
      integer j
      integer k
      double precision time1
      double precision time2

      call cpu_time ( time1 )

      do j = 1, n3
        do k = 1, n2
          do i = 1, n1
            a(i,j) = a(i,j) + b(i,k) * c(k,j)
          end do
        end do
      end do

      call cpu_time ( time2 )

      cpu_seconds = time2 - time1

      return
      end
      subroutine mxm_kij ( n1, n2, n3, a, b, c, cpu_seconds )

c*********************************************************************72
c
cc MXM_KIJ computes A = B * C using DO K, DO I, DO J loops.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 October 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, N3, define the orders of the
c    matrices.
c
c    Output, double precision A(N1,N3), the result matrix.
c
c    Input, double precision B(N1,N2), C(N2,N3), the factor matrices.
c
c    Output, double precision CPU_SECONDS, the elapsed CPU time.
c
      implicit none

      integer n1
      integer n2
      integer n3

      double precision a(n1,n3)
      double precision b(n1,n2)
      double precision c(n2,n3)
      double precision cpu_seconds
      integer i
      integer j
      integer k
      double precision time1
      double precision time2

      call cpu_time ( time1 )

      do k = 1, n2
        do i = 1, n1
          do j = 1, n3
            a(i,j) = a(i,j) + b(i,k) * c(k,j)
          end do
        end do
      end do

      call cpu_time ( time2 )

      cpu_seconds = time2 - time1

      return
      end
      subroutine mxm_kji ( n1, n2, n3, a, b, c, cpu_seconds )

c*********************************************************************72
c
cc MXM_KJI computes A = B * C using DO K, DO J, DO I loops.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 October 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, N3, define the orders of the
c    matrices.
c
c    Output, double precision A(N1,N3), the result matrix.
c
c    Input, double precision B(N1,N2), C(N2,N3), the factor matrices.
c
c    Output, double precision CPU_SECONDS, the elapsed CPU time.
c
      implicit none

      integer n1
      integer n2
      integer n3

      double precision a(n1,n3)
      double precision b(n1,n2)
      double precision c(n2,n3)
      double precision cpu_seconds
      integer i
      integer j
      integer k
      double precision time1
      double precision time2

      call cpu_time ( time1 )

      do k = 1, n2
        do j = 1, n3
          do i = 1, n1
            a(i,j) = a(i,j) + b(i,k) * c(k,j)
          end do
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
