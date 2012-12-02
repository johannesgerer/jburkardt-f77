      program main

c*********************************************************************72
c
cc POISSON_SIMULATION_TEST tests POISSON_SIMULATION.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'POISSON_SIMULATION_TEST'
      write ( *, '(a)' ) '  FORTRAN77 version.'
      write ( *, '(a)' ) '  Test the POISSON_SIMULATION library.'

      call test01 ( )
      call test02 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'POISSON_SIMULATION_TEST'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ''
      call timestamp ( )

      return
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 simulates waiting for a given number of events.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer bin_num
      parameter ( bin_num = 30 )
      integer event_num
      parameter ( event_num = 1000 )

      character * ( 80 ) command_filename
      integer command_unit
      character * ( 80 ) data_filename
      integer data_unit
      integer f_bin(bin_num)
      integer i
      integer j
      double precision lambda
      integer seed
      double precision t(0:event_num)
      double precision w(0:event_num)
      double precision w_ave
      double precision w_bin(bin_num)
      double precision w_max
      double precision w_min
      double precision width

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' ) 
     &  '  POISSON_FIXED_EVENTS simulates a Poisson process'
      write ( *, '(a)' ) 
     &  '  until a given number of events have occurred.'
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 
     &  '  Simulate a Poisson process, for which, on average,'
      write ( *, '(a)' ) '  LAMBDA events occur per unit time.'
      write ( *, '(a)' ) 
     &  '  Run until you have observed EVENT_NUM events.'
     
      lambda = 0.5D+00
      seed = 123456789

      write ( *, '(a)' ) ''
      write ( *, '(a,g14.6)' ) '  LAMBDA = ', lambda
      write ( *, '(a,i6)' ) '  EVENT_NUM = ', event_num

      call poisson_fixed_events ( lambda, event_num, seed, t, w )

      call r8vec_min ( event_num, w, w_max )
      call r8vec_mean ( event_num, w, w_ave )
      call r8vec_max ( event_num, w, w_min )

      write ( *, '(a)' ) ''
      write ( *, '(a,g14.6)' ) '  Minimum wait = ', w_min
      write ( *, '(a,g14.6)' ) '  Average wait = ', w_ave
      write ( *, '(a,g14.6)' ) '  Maximum wait = ', w_max

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) ' Count            Time            Wait'
      write ( *, '(a)' ) ''
      do i = 0, 5
        write ( *, '(2x,i4,2x,g14.6,2x,g14.6)' ) i, t(i), w(i)
      end do
      write ( *, '(a)' ) '  ....  ..............  ..............'
      do i = event_num - 5, event_num 
        write ( *, '(2x,i4,2x,g14.6,2x,g14.6)' ) i, t(i), w(i)
      end do
c
c  Create the data file.
c
      call get_unit ( data_unit )

      data_filename = 'poisson_timeline_data.txt'

      open ( unit = data_unit, file = data_filename, 
     &  status = 'replace' )

      do i = 0, event_num
        write ( data_unit, '(2x,g14.6,2x,i8)' ) t(i), i
      end do
      close ( unit = data_unit )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Data stored in "' 
     &  // trim ( data_filename ) // '".'
c
c  Create the command file.
c
      call get_unit ( command_unit )

      command_filename = 'poisson_timeline_commands.txt'

      open ( unit = command_unit, file = command_filename, 
     &  status = 'replace' )

      write ( command_unit, '(a)' ) 
     &  '# poisson_timeline_commands.txt'
      write ( command_unit, '(a)' ) '#'
      write ( command_unit, '(a)' ) '# Usage:'
      write ( command_unit, '(a)' ) 
     &  '#  gnuplot < poisson_timeline_commands.txt'
      write ( command_unit, '(a)' ) '#'
      write ( command_unit, '(a)' ) 'set term png'
      write ( command_unit, '(a)' ) 'set output "poisson_timeline.png"'
      write ( command_unit, '(a)' ) 'set style data lines'
      write ( command_unit, '(a)' ) 'set xlabel "Time"'
      write ( command_unit, '(a)' ) 'set ylabel "Number of Events"'
      write ( command_unit, '(a)' ) 
     &  'set title "Observation of Fixed Number of Poisson Events"'
      write ( command_unit, '(a)' ) 'set grid'
      write ( command_unit, '(a,f8.2,a)' ) 
     &  'plot "poisson_timeline_data.txt" using 1:2 lw 2'
      write ( command_unit, '(a)' ) 'quit'

      close ( unit = command_unit )

      write ( *, '(a)' ) '  Plot commands stored in "' 
     &  // trim ( command_filename ) // '".'
c
c  Determine bin information.
c
      call r8vec_midspace ( bin_num, w_min, w_max, w_bin )

      do i = 1, bin_num
        f_bin(i) = 0
      end do

      do i = 0, event_num
        j = 1 + int ( dble ( bin_num ) * ( w(i) - w_min ) 
     *    / ( w_max - w_min ) )
        j = min ( j, bin_num )
        f_bin(j) = f_bin(j) + 1
      end do
c
c  Create the data file.
c
      call get_unit ( data_unit )

      data_filename = 'poisson_times_data.txt'

      open ( unit = data_unit, file = data_filename, 
     &  status = 'replace' )

      do i = 1, bin_num
        write ( data_unit, '(2x,g14.6,2x,i6)' ) w_bin(i), f_bin(i)
      end do
      close ( unit = data_unit )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Data stored in "' 
     &  // trim ( data_filename ) // '".'
c
c  Create the command file.
c
      call get_unit ( command_unit )

      command_filename = 'poisson_times_commands.txt'

      open ( unit = command_unit, file = command_filename, 
     &  status = 'replace' )

      write ( command_unit, '(a)' ) '# poisson_times_commands.txt'
      write ( command_unit, '(a)' ) '#'
      write ( command_unit, '(a)' ) '# Usage:'
      write ( command_unit, '(a)' ) 
     &  '#  gnuplot < poisson_times_commands.txt'
      write ( command_unit, '(a)' ) '#'
      write ( command_unit, '(a)' ) 'set term png'
      write ( command_unit, '(a)' ) 'set output "poisson_times.png"'
      write ( command_unit, '(a)' ) 'set xlabel "Waiting Time"'
      write ( command_unit, '(a)' ) 'set ylabel "Frequency"'
      write ( command_unit, '(a)' ) 
     &  'set title "Waiting Times Observed Over Fixed Time"'
      write ( command_unit, '(a)' ) 'set grid'
      write ( command_unit, '(a)' ) 'set style fill solid'
      width = 0.85D+00 * ( w_max - w_min ) / dble ( bin_num )
      write ( command_unit, '(a,f8.2,a)' ) 
     &  'plot "poisson_times_data.txt" using 1:2:(', 
     &  width, ') with boxes'
      write ( command_unit, '(a)' ) 'quit'

      close ( unit = command_unit )

      write ( *, '(a)' ) '  Plot commands stored in "' 
     &  // trim ( command_filename ) // '".'

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 simulates waiting for a given length of time.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer bin_num
      parameter ( bin_num = 30 )
      integer test_num
      parameter ( test_num = 20000 )

      character * ( 80 ) command_filename
      integer command_unit
      character * ( 80 ) data_filename
      integer data_unit
      double precision f_bin(bin_num)
      integer i
      double precision lambda
      integer n(test_num)
      double precision n_bin(bin_num)
      integer n_max
      double precision n_max_r8
      double precision n_mean
      integer n_min
      double precision n_min_r8
      double precision n_var
      integer seed
      double precision t
      integer test
      double precision w

      lambda = 0.5D+00
      t = 1000.0D+00
      seed = 123456789

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST02:'
      write ( *, '(a)' ) 
     &  '  POISSON_FIXED_EVENTS simulates a Poisson process'
      write ( *, '(a)' ) 
     &  '  counting the number of events that occur during'
      write ( *, '(a)' ) '  a given time.'
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 
     &  '  Simulate a Poisson process, for which, on average,'
      write ( *, '(a)' ) '  LAMBDA events occur per unit time.'
      write ( *, '(a,g14.6,a)' ) 
     &  '  Run for a total of ', t, ' time units.'
      write ( *, '(a,g14.6)' ) '  LAMBDA = ', lambda

      do test = 1, test_num
        call poisson_fixed_time ( lambda, t, seed, n(test) )
      end do

      call i4vec_mean ( test_num, n, n_mean )
      call i4vec_variance ( test_num, n, n_var )
      write ( *, '(a)' ) ''
      write ( *, '(a,g14.6)' ) '  Mean number of events = ', n_mean
      write ( *, '(a,g14.6)' ) '  Variance = ', n_var  
      write ( *, '(a,g14.6)' ) '  STD = ', sqrt ( n_var )

      call i4vec_min ( test_num, n, n_min )
      call i4vec_max ( test_num, n, n_max )
      n_min_r8 = dble ( n_min )
      n_max_r8 = dble ( n_max )

      call r8vec_midspace ( bin_num, n_min_r8, n_max_r8, n_bin )

      do i = 1, bin_num
        f_bin(i) = 0
      end do

      do test = 1, test_num
        i = 1 + int ( dble ( bin_num * ( n(test) - n_min ) ) 
     &    / dble ( n_max - n_min ) )
        i = min ( i, bin_num )
        f_bin(i) = f_bin(i) + 1
      end do
c
c  Create the data file.
c
      call get_unit ( data_unit )

      data_filename = 'poisson_events_data.txt'

      open ( unit = data_unit, file = data_filename, 
     &  status = 'replace' )

      do i = 1, bin_num
        write ( data_unit, '(2(2x,g14.6))' ) n_bin(i), f_bin(i)
      end do
      close ( unit = data_unit )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Data stored in "' 
     &  // trim ( data_filename ) // '".'
c
c  Create the command file.
c
      call get_unit ( command_unit )

      command_filename = 'poisson_events_commands.txt'

      open ( unit = command_unit, file = command_filename, 
     &  status = 'replace' )

      write ( command_unit, '(a)' ) 
     &  '# poisson_events_commands.txt'
      write ( command_unit, '(a)' ) '#'
      write ( command_unit, '(a)' ) '# Usage:'
      write ( command_unit, '(a)' ) 
     &  '#  gnuplot < poisson_events_commands.txt'
      write ( command_unit, '(a)' ) '#'
      write ( command_unit, '(a)' ) 'set term png'
      write ( command_unit, '(a)' ) 'set output "poisson_events.png"'
      write ( command_unit, '(a)' ) 'set xlabel "Number of Events"'
      write ( command_unit, '(a)' ) 'set ylabel "Frequency"'
      write ( command_unit, '(a)' ) 
     &  'set title "Number of Poisson Events Over Fixed Time"'
      write ( command_unit, '(a)' ) 'set grid'
      write ( command_unit, '(a)' ) 'set style fill solid'
      w = 0.85D+00 * ( n_max - n_min ) / dble ( bin_num )
      write ( command_unit, '(a,f8.2,a)' ) 
     &  'plot "poisson_events_data.txt" using 1:2:(', w, ') with boxes'
      write ( command_unit, '(a)' ) 'quit'

      close ( unit = command_unit )

      write ( *, '(a)' ) '  Plot commands stored in "' 
     &  // trim ( command_filename ) // '".'
      
      return
      end
