      subroutine simplex_lattice_point_next ( n, t, more, x )

c*********************************************************************72
c
cc SIMPLEX_LATTICE_POINT_NEXT generates lattice points in a simplex.
c
c  Discussion:
c
c    The simplex is defined by N-dimensional points X such that:
c
c        0 <= X(1:N)
c
c    and
c
c      sum ( X(1:N) ) <= T
c
c    where T is an integer.
c
c    Lattice points are points X which satisfy the simplex conditions and
c    for which all the components are integers.
c
c    This routine generates all the lattice points in a given simplex, one at 
c    a time, in a reverse lexicographic order.
c
c    To use the routine, initialize by setting N and T to appropriate values, 
c    and MORE to FALSE.  The initial value of X is not important.
c
c    Call the routine. On return, X will contain the first lattice point in 
c    the simplex.  If MORE is TRUE, then the routine may be called again to 
c    get the next point.  In fact, as long as the output value of MORE is 
c    TRUE, there is at least one more lattice point that can be found by 
c    making another call.  When MORE is returned as FALSE, then there are no 
c    more lattice points; the value of X returned at that time is the 
c    "last" such point.
c
c    During the computation of a sequence of lattice points, the user should 
c    not change the values of N, T, MORE or X.  
c
c    The output for N = 3, T = 4 would be:
c
c       1    4  0  0
c       2    3  1  0
c       3    3  0  1
c       4    2  2  0
c       5    2  1  1
c       6    2  0  2
c       7    1  3  0
c       8    1  2  1
c       9    1  1  2
c      10    1  0  3
c      11    0  4  0
c      12    0  3  1
c      13    0  2  2
c      14    0  1  3
c      15    0  0  4
c
c    The number of lattice points is (T+N-1)! / ( T! * (N-1)! )
c
c    Thus, for N = 3, T = 4, we have
c
c      number = 6! / ( 4! * 2! ) = 720 / ( 24 * 2 ) = 15
c
c  Modified:
c
c    10 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Scott Chasalow, Richard Brand,
c    Algorithm AS 299:
c    Generation of Simplex Lattice Points,
c    Applied Statistics,
c    Volume 44, Number 4, 1995, pages 534-545.
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
c    Input, integer N, the spatial dimension.
c    N must be positive.
c
c    Input, integer T, the characteristic of the simplex.
c    T must be nonnegative.
c
c    Input/output, logical MORE, initialized to FALSE by the user to
c    begin a sequence of calculations, returned by the routine as TRUE,
c    if there are more values of X that can be calculated, or FALSE
c    if the accompanying value of X is the last one for this sequence.
c
c    Input/output, integer X(N), not initialized by the user, but not
c    changed by the user on subsequent calls.  The routine returns
c    a new point on each call, and on subsequent calls uses the input
c    value (old point) to compute the output value (next point).
c
      implicit none

      integer n

      integer i
      integer j
      integer k
      logical more
      integer t
      integer x(n)

      if ( .not. more ) then

        if ( n .lt. 1 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SIMPLEX_LATTICE_POINT_NEXT - Fatal error!'
          write ( *, '(a)' ) '  N < 1.'
          stop
        end if

        if ( t .lt. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SIMPLEX_LATTICE_POINT_NEXT - Fatal error!'
          write ( *, '(a)' ) '  T < 0.'
          stop
        end if

        more = .true.
        j = 1

        x(1) = t
        do i = 2, n
          x(i) = 0
        end do
c
c  The first point can actually also be the last!
c
        if ( n .eq. 1 ) then
          more = .false.
        end if

      else
c
c  Search X(N-1 down to 1) for the first nonzero element.
c  If none, then terminate.  (This should not happen.)
c  Otherwise, set J to this index.
c  Decrement X(J) by 1.
c  Set X(J+1:N) to (T-X(1:J),0,0,...0).
c
        j = n

        do i = n-1, 1, -1

          if ( 0 .lt. x(i) ) then

            j = i
            x(j) = x(j) - 1

            x(j+1) = t
            do k = 1, j
              x(j+1) = x(j+1) - x(k)
            end do

            do k = j+2, n
              x(k) = 0
            end do

            if ( x(n) .eq. t ) then
              more = .false.
            end if

            return

          end if

        end do

        if ( j .eq. n ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SIMPLEX_LATTICE_POINT_NEXT - Fatal error!'
          write ( *, '(a)' ) '  The input X vector is nonpositive in'
          write ( *, '(a)' ) '  all entries except possibly the last'
          write ( *, '(a)' ) '  one.'
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Perhaps the user has miscalled the'
          write ( *, '(a)' ) '  routine or altered data between calls.'
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'ABNORMAL TERMINATION.'
          stop
        end if

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

