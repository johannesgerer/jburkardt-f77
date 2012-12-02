      program main

c*********************************************************************72
c
cc MAIN is the main program for SATISFY.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 March 2009
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Michael Quinn,
c    Parallel Programming in C with MPI and OpenMP,
c    McGraw-Hill, 2004,
c    ISBN13: 978-0071232654,
c    LC: QA76.73.C15.Q55.
c
      implicit none

      integer n
      parameter ( n = 23 )

      integer bvec(n)
      integer circuit_value
      integer i
      integer ihi
      integer j
      integer solution_num
      integer value

      write ( *, '(a)' ) ' '
      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SATISFY'
      write ( *, '(a)' ) '  FORTRAN77 version' 
      write ( *, '(a)' ) 
     &  '  We have a logical function of N logical arguments.'
      write ( *, '(a)' ) 
     &  '  We do an exhaustive search of all 2^N possibilities,'
      write ( *, '(a)' ) 
     &  '  seeking those inputs that make the function TRUE.'
c
c  Compute the number of binary vectors to check.
c
      ihi = 2**n

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 
     &  '  The number of logical variables is N = ',  n
      write ( *, '(a,i8)' ) 
     &  '  The number of input vectors to check is ', ihi
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '   #       Index    ' // 
     &  '---------Input Values------------------------'
      write ( *, '(a)' ) ' '
c
c  Check every possible input vector.
c
      solution_num = 0

      do i = 0, ihi - 1

        call i4_to_bvec ( i, n, bvec )

        value = circuit_value ( n, bvec )

        if ( value == 1 ) then
          solution_num = solution_num + 1

          write ( *, '(2x,i2,2x,i10,3x,23i2)' )  
     &      solution_num, i, ( bvec(j), j = 1, n )
        end if

      end do
c
c  Report.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 
     &  '  Number of solutions found was ', solution_num
c
c  Shut down.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SATISFY'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      function circuit_value ( n, bvec )

c*********************************************************************72
c
cc CIRCUIT_VALUE returns the value of a circuit for a given input set.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 March 2009
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Michael Quinn,
c    Parallel Programming in C with MPI and OpenMP,
c    McGraw-Hill, 2004,
c    ISBN13: 978-0071232654,
c    LC: QA76.73.C15.Q55.
c
c  Parameters:
c
c    Input, integer N, the length of the input vector.
c
c    Input, integer BVEC(N), the binary inputs.
c
c    Output, integer CIRCUIT_VALUE, the output of the circuit.
c
      implicit none

      integer n

      integer bvec(n)
      integer circuit_value
      logical value

      value = ( bvec(1)  == 1 .or. bvec(2)  == 1 ) 
     &  .and. ( bvec(2)  == 0 .or. bvec(4)  == 0 ) 
     &  .and. ( bvec(3)  == 1 .or. bvec(4)  == 1 ) 
     &  .and. ( bvec(4)  == 0 .or. bvec(5)  == 0 ) 
     &  .and. ( bvec(5)  == 1 .or. bvec(6)  == 0 ) 
     &  .and. ( bvec(6)  == 1 .or. bvec(7)  == 0 ) 
     &  .and. ( bvec(6)  == 1 .or. bvec(7)  == 1 ) 
     &  .and. ( bvec(7)  == 1 .or. bvec(16) == 0 ) 
     &  .and. ( bvec(8)  == 1 .or. bvec(9)  == 0 ) 
     &  .and. ( bvec(8)  == 0 .or. bvec(14) == 0 ) 
     &  .and. ( bvec(9)  == 1 .or. bvec(10) == 1 ) 
     &  .and. ( bvec(9)  == 1 .or. bvec(10) == 0 ) 
     &  .and. ( bvec(10) == 0 .or. bvec(11) == 0 ) 
     &  .and. ( bvec(10) == 1 .or. bvec(12) == 1 ) 
     &  .and. ( bvec(11) == 1 .or. bvec(12) == 1 ) 
     &  .and. ( bvec(13) == 1 .or. bvec(14) == 1 ) 
     &  .and. ( bvec(14) == 1 .or. bvec(15) == 0 ) 
     &  .and. ( bvec(15) == 1 .or. bvec(16) == 1 ) 
     &  .and. ( bvec(15) == 1 .or. bvec(17) == 1 ) 
     &  .and. ( bvec(18) == 1 .or. bvec(2)  == 1 ) 
     &  .and. ( bvec(19) == 1 .or. bvec(1)  == 0 ) 
     &  .and. ( bvec(20) == 1 .or. bvec(2)  == 1 ) 
     &  .and. ( bvec(20) == 1 .or. bvec(19) == 0 ) 
     &  .and. ( bvec(20) == 0 .or. bvec(10) == 0 ) 
     &  .and. ( bvec(1)  == 1 .or. bvec(18) == 1 ) 
     &  .and. ( bvec(2)  == 0 .or. bvec(21) == 1 ) 
     &  .and. ( bvec(22) == 0 .or. bvec(21) == 1 ) 
     &  .and. ( bvec(23) == 0 .or. bvec(21) == 1 ) 
     &  .and. ( bvec(22) == 0 .or. bvec(21) == 0 ) 
     &  .and. ( bvec(23) == 1 .or. bvec(21) == 0 )

      if ( value ) then
        circuit_value = 1
      else
        circuit_value = 0
      end if

      return
      end
      subroutine i4_to_bvec ( i4, n, bvec )

c*********************************************************************72
c
cc I4_TO_BVEC converts an integer into a binary vector.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 March 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I4, the integer.
c
c    Input, integer N, the dimension of the vector.
c
c    Output, integer BVEC(N), the vector of binary remainders.
c
      implicit none

      integer n

      integer bvec(n)
      integer i
      integer i4
      integer i4_copy

      i4_copy = i4

      do i = n, 1, -1
        bvec(i) = mod ( i4_copy, 2 )
        i4_copy = i4_copy / 2
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
