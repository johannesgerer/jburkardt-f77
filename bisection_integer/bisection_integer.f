      subroutine bisection_integer ( f, a, b, c, fc )

c*********************************************************************72
c
cc BISECTION_INTEGER seeks an integer root using bisection.
c
c  Discussion:
c
c    A function F(X) confined to integer arguments is given, with an
c    interval [A,B] over which F changes sign.  An integer C is sought
c    such that A <= C <= B and F(C) = 0.
c
c    Because we are restricted to integer arguments, it may the case that
c    there is no such C.
c
c    This routine proceeds by a form of bisection, in which the enclosing
c    interval is restricted to be defined by integer values.
c
c    If the user has given a true change of sign interval [A,B], and if,
c    in the interval, there is a single integer value C for which F(C) = 0,
c    with the additional restrictions that F(C-1) and F(C+1) are of opposite
c    signs, then this procedure should locate and return C.
c
c    In particular, if the function F is monotone, and there is an integer
c    solution C in the interval, then this procedure will find it.
c
c    However, in general, even if there is an integer C in the interval,
c    such that F(C) = 0, this procedure may be unable to find it, particularly
c    if there are also nonintegral solutions within the same interval.
c
c    While any integer function can be used with this program, the bisection
c    approach is most useful if the integer function is monotone, or
c    varies slowly, or can be regarded as the restriction to integer arguments
c    of a continuous (and smoothly varying) function of a real argument.
c    In such cases, knowing that F is negative at A and positive at B
c    suggests that F generally increases from A to B, and might attain 
c    the value 0 at some intermediate argument C.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, external integer F, the name of a user-supplied 
c    procedure that evaluates the function, of the form
c      function f ( c )
c      integer c, f
c
c   Input, integer A, B, two arguments that define a change of
c   sign interval for F.  In other words, F(A) and F(B) must be of opposite
c   sign.
c
c   Output, integer C, FC, the candidate for the root, as 
c   determined by the program, and its function value.  If FC is not zero,
c   then the procedure did not find a root in the interval, and C is only
c   an "approximate" root.
c
      integer a
      integer b
      integer c
      integer f
      external f
      integer fa
      integer fb
      integer fc
      integer t
c
c  Ensure that F(A) < 0 < F(B).
c
      fa = f(a)
      fb = f(b)

      if ( fa .eq. 0 ) then
        c = a
        fc = fa
      else if ( fb .eq. 0 ) then
        c = b
        fc = fb
      else if ( fa .lt. 0 .and. 0 .lt. fb ) then

      else if ( fb .lt. 0 .and. 0 .lt. fa ) then
        t = a
        a = b
        b = t
        t = fa
        fa = fb
        fb = t
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BISECTION_INTEGER - Fatal error!'
        write ( *, '(a)' ) '  No change of sign interval supplied.'
        write ( *, '(a,i8,a,i8)' ) '  F(', a, ') = ', fa
        write ( *, '(a,i8,a,i8)' ) '  F(', b, ') = ', fb
        stop
      end if
c
c  Bisection.
c
10    continue

      if ( 1 .lt. abs ( b - a ) ) then

        c = ( a + b ) / 2
        fc = f(c)

        if ( fc .eq. 0 ) then
          return
        else if ( fc .lt. 0 ) then
          a = c
          fa = fc
        else if ( 0 .lt. fc ) then
          b = c
          fb = fc
        end if

        go to 10

      end if
c
c  Interval is empty, with FA < 0 and 0 < FB.
c  Bisection did not produce an integer solution.
c  Return the argument with smallest function norm.
c
      if ( - fa .lt. fb ) then
        c = a
        fc = fa
      else
        c = b
        fc = fb
      end if

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
