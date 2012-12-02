      program main

c*********************************************************************72
c
cc MAIN is the main program for SELECT_PRB.
c
cc  Discussion:
c
c    SELECT_PRB tests the routines in SELECT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer family

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SELECT_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the SELECT library.'

      do family = 1, 7

        write ( *, '(a)' ) ' '
        write ( *, '(a,i1)' ) '  FAMILY = ', family

        if ( family .eq. 1 ) then
          write ( *, * ) '  K subsets of an N set.'
        else if ( family .eq. 2 ) then
          write ( *, * ) '  Partitions of N objects into K classes.'
        else if ( family .eq. 3 ) then
          write ( *, * ) '  Permutations of N objects with K cycles.'
        else if ( family .eq. 4 ) then
          write ( *, * ) '  Vector subspaces of dimension K '
          write ( *, * ) '  over N-dimensional space over the field '
          write ( *, * ) '  of order Q (where Q = 2 ).'
        else if ( family .eq. 5 ) then
          write ( *, * ) '  Permutations of N letters with K runs.'
        else if ( family .eq. 6 ) then
          write ( *, * ) '  Partitions of N whose largest part is K.'
        else if ( family .eq. 7 ) then
          write ( *, * ) '  Compositions of N into K parts.'
        end if

        call test01 ( family )

      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SELECT_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( family )

c*********************************************************************72
c
cc TEST01 carries out all the TASKS for a particular FAMILY.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer family
      integer task
c
c  Tasks 2, 3 and 4 not ready yet!
c
      do task = 1, 1

        write ( *, '(a)' ) ' '
        write ( *, '(a,i1)' ) '  TASK = ', task

        if ( task .eq. 1 ) then
          write ( *, * ) '  Present each object of the family.'
        else if ( task .eq. 2 ) then
          write ( *, * ) '  Rank a given object of the family.'
        else if ( task .eq. 3 ) then
          write ( *, * ) '  Produce an object of given rank.'
        else if ( task .eq. 4 ) then
          write ( *, * ) '  Select an object at random.'
        end if
        write ( *, '(a)' ) ' '

        call test02 ( family, task )

      end do

      return
      end
      subroutine test02 ( family, task )

c*********************************************************************72
c
cc TEST02 carries out a particular TASK for a particular FAMILY.
c
c  Discussion:
c
c    The SELECT program seems to write data into entries (N+1) through (N+K)
c    of EDGE, MU and NU, beyond the original formal dimension of N.  This extra data
c    is of no interest to the user, but it must be allocated in the
c    calling program to avoid memory overwriting!
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer k
      integer n

      parameter ( n = 5 )
      parameter ( k = 4 )

      integer b(10,10)
      integer edge(n+k)
      integer family
      integer i
      integer j
      integer m
      integer mu(n+k)
      logical newone
      integer nu(n+k)
      integer rank
      integer task

      if ( task .eq. 1 ) then

        newone = .false.

        do i = 1, n
          edge(i) = 0
          mu(i) = 0
          nu(i) = 0
        end do

        do

          call select ( family, task, n, k, mu, nu, edge, m,
     &        newone, rank, b )

          if ( .not. newone ) then
            exit
          end if

          write ( *, '(i3,2x,5i1,2x,6i1,2x,6i1)' )
     &        rank, ( edge(j), j = 1,n),
     &        n, ( mu(j), j = 1,n),
     &        k, ( nu(j), j = 1,n)

        end do

      else if ( task .eq. 2 ) then
      else if ( task .eq. 3 ) then
      else if ( task .eq. 4 ) then

      end if

      return
      end

