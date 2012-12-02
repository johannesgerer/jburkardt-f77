      program main

c*********************************************************************72
c
cc TASK_DIVISION_PRB tests TASK_DIVISION.
c
c  Discussion:
c
c    This program simply demonstrates how one might automate the
c    assignment of T tasks to P processors, assuming that the assignment
c    is to be beforehand.
c
c    In that case, we just want to make sure that we assign each task
c    to a processor, that we assign about the same number of tasks
c    to each processor, and that we assign each processor a contiguous
c    range of tasks, say tasks I_LO to I_HI.
c
c    The routine that is called simulates this process.
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
      implicit none

      integer proc_first
      integer proc_last
      integer task_number

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TASK_DIVISION_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Demonstrate how to automate the division of'
      write ( *, '(a)' ) '  of T tasks among a range of P processors'
      write ( *, '(a)' ) '  indexed from PROC_FIRST to PROC_LAST.'

      task_number = 23
      proc_first = 0
      proc_last = 3
      call task_division ( task_number, proc_first, proc_last )

      task_number = 17
      proc_first = 1
      proc_last = 6
      call task_division ( task_number, proc_first, proc_last )

      task_number = 17
      proc_first = 4
      proc_last = 6
      call task_division ( task_number, proc_first, proc_last )

      task_number = 5
      proc_first = -2
      proc_last = 6
      call task_division ( task_number, proc_first, proc_last )

      task_number = 5
      proc_first = 0
      proc_last = 4
      call task_division ( task_number, proc_first, proc_last )

      task_number = 5
      proc_first = 0
      proc_last = 0
      call task_division ( task_number, proc_first, proc_last )

      task_number = 1000
      proc_first = 1
      proc_last = 17
      call task_division ( task_number, proc_first, proc_last )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TASK_DIVISION_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
