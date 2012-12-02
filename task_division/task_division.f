      function i4_div_rounded ( a, b )

c*********************************************************************72
c
cc I4_DIV_ROUNDED computes the rounded result of I4 division.
c
c  Discussion:
c
c    This routine computes C = A / B, where A, B and C are integers
c    and C is the closest integer value to the exact real result.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer A, B, the number to be divided,
c    and the divisor.
c
c    Output, integer I4_DIV_ROUNDED, the rounded result
c    of the division.
c
      implicit none

      integer a
      integer a_abs
      integer b
      integer b_abs
      integer i4_div_rounded
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer value

      if ( a .eq. 0 .and. b .eq. 0 ) then

        value = i4_huge
 
      else if ( a .eq. 0 ) then

        value = 0

      else if ( b .eq. 0 ) then

        if ( a .lt. 0 ) then
          value = - i4_huge
        else
          value = + i4_huge
        end if

      else

        a_abs = abs ( a )
        b_abs = abs ( b )

        value = a_abs / b_abs
c
c  Round the value.
c
        if ( ( 2 * value + 1 ) * b_abs .lt. 2 * a_abs ) then
          value = value + 1
        end if
c
c  Set the sign.
c
        if ( ( a .lt. 0 .and. 0 .lt. b ) .or. 
     &       ( 0 .lt. a .and. b .lt. 0 ) ) then
          value = - value
        end if

      end if

      i4_div_rounded = value

      return
      end
      subroutine task_division ( task_number, proc_first, proc_last )

c*********************************************************************72
c
cc TASK_DIVISION divides tasks among processors.
c
c  Discussion:
c
c    This routine assigns each of T tasks to P processors, assuming that 
c    the assignment is to be beforehand.
c
c    In that case, we just want to make sure that we assign each task
c    to a processor, that we assign about the same number of tasks
c    to each processor, and that we assign each processor a contiguous
c    range of tasks, say tasks I_LO to I_HI.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer TASK_NUMBER, the number of tasks.
c
c    Input, integer PROC_FIRST, PROC_LAST, the first and last processors.
c
      implicit none

      integer i_hi
      integer i_lo
      integer i4_div_rounded
      integer proc
      integer proc_first
      integer proc_last
      integer proc_number
      integer proc_remain
      integer task_number
      integer task_proc
      integer task_remain

      proc_number = proc_last + 1 - proc_first

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TASK_DIVISION'
      write ( *, '(a)' ) '  Divide T tasks among P processors.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of tasks T = ', task_number
      write ( *, '(a,i8)' ) '  Number of processors P = ', proc_number
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  P_FIRST = ', proc_first
      write ( *, '(a,i8)' ) '  P_LAST =  ', proc_last
      write ( *, '(a)' ) ' '
      write ( *, '(a)') '              Number of   First     Last'
      write ( *, '(a)' ) ' Processor     Tasks     Last      Task'
      write ( *, '(a)' ) ' '

      i_hi = 0

      task_remain = task_number
      proc_remain = proc_number

      do proc = proc_first, proc_last

        task_proc = i4_div_rounded ( task_remain, proc_remain )

        proc_remain = proc_remain - 1
        task_remain = task_remain - task_proc

        i_lo = i_hi + 1
        i_hi = i_hi + task_proc

        write ( *, '(2x,i8,2x,i8,2x,i8,2x,i8)' ) 
     &    proc, task_proc, i_lo, i_hi

      end do

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
