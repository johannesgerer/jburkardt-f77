      program main

c*********************************************************************72
c
cc TOMS456_PRB tests ROUTNG.
c
c  Modified:
c
c    12 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      integer n

      parameter ( m = 15 )
      parameter ( n = 15 )

      integer d(m,m)
      integer en
      integer i
      integer ip1
      integer j
      integer l
      integer n2
      integer p(n)
      integer r
      integer sn
      integer total

      save d

      data d /
     &   0, 29, 82, 46, 68, 52, 72, 42, 51, 55, 29, 74, 23, 72, 46,
     &  29,  0, 55, 46, 42, 43, 43, 23, 23, 31, 41, 51, 11, 52, 21,
     &  82, 55,  0, 68, 46, 55, 23, 43, 41, 29, 79, 21, 64, 31, 51,
     &  46, 46, 68,  0, 82, 15, 72, 31, 62, 42, 21, 51, 51, 43, 64,
     &  68, 42, 46, 82,  0, 74, 23, 52, 21, 46, 82, 58, 46, 65, 23,
     &  52, 43, 55, 15, 74,  0, 61, 23, 55, 31, 33, 37, 51, 29, 59,
     &  72, 43, 23, 72, 23, 61,  0, 42, 23, 31, 77, 37, 51, 46, 33,
     &  42, 23, 43, 31, 52, 23, 42,  0, 33, 15, 37, 33, 33, 31, 37,
     &  51, 23, 41, 62, 21, 55, 23, 33,  0, 29, 62, 46, 29, 51, 11,
     &  55, 31, 29, 42, 46, 31, 31, 15, 29,  0, 51, 21, 41, 23, 37,
     &  29, 41, 79, 21, 82, 33, 77, 37, 62, 51,  0, 65, 42, 59, 61,
     &  74, 51, 21, 51, 58, 37, 37, 33, 46, 21, 65,  0, 61, 11, 55,
     &  23, 11, 64, 51, 46, 51, 51, 33, 29, 41, 42, 61,  0, 62, 23,
     &  72, 52, 31, 43, 65, 29, 46, 31, 51, 23, 59, 11, 62,  0, 59,
     &  46, 21, 51, 64, 23, 59, 33, 37, 11, 37, 61, 55, 23, 59,  0 /

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS456_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS 456 library.'

      sn = 1
      en = 1
      r = 5
      do i = 1, n
        p(i) = i
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Start node SN =    ', sn
      write ( *, '(a,i8)' ) '  End node EN =      ', en
      write ( *, '(a,i8)' ) '  Number of trials = ', r

      call routng ( n, p, sn, en, m, d, l, r )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' )
     &  '  The length of the optimal connection is ', l

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Connection sequence:'
      write ( *, '(a)' ) ' '

      total = 0.0
      do i = 1, n
        write ( *, '(2x,i8,2x,i8)' ) p(i), total
        if ( i .lt. n ) then
          total = total + d(p(i),p(i+1))
        end if
      end do

      sn = 1
      en = 13
      r = 5
      do i = 1, n
        p(i) = i
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Start node SN =    ', sn
      write ( *, '(a,i8)' ) '  End node EN =      ', en
      write ( *, '(a,i8)' ) '  Number of trials = ', r

      call routng ( n, p, sn, en, m, d, l, r )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' )
     &  '  The length of the optimal connection is ', l

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Connection sequence:'
      write ( *, '(a)' ) ' '

      total = 0.0
      do i = 1, n
        write ( *, '(2x,i8,2x,i8)' ) p(i), total
        if ( i .lt. n ) then
          total = total + d(p(i),p(i+1))
        end if
      end do

      sn = 1
      en = 0
      r = 5
      do i = 1, n
        p(i) = i
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Start node SN =    ', sn
      write ( *, '(a,i8)' ) '  End node EN =      ', en
      write ( *, '(a,i8)' ) '  Number of trials = ', r

      call routng ( n, p, sn, en, m, d, l, r )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' )
     &  '  The length of the optimal connection is ', l

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Connection sequence:'
      write ( *, '(a)' ) ' '

      total = 0.0
      do i = 1, n
        write ( *, '(2x,i8,2x,i8)' ) p(i), total
        if ( i .lt. n ) then
          total = total + d(p(i),p(i+1))
        end if
      end do

      n2 = 5
      sn = 1
      en = 5
      r = 5
      do i = 1, n2
        p(i) = i
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Start node SN =    ', sn
      write ( *, '(a,i8)' ) '  End node EN =      ', en
      write ( *, '(a,i8)' ) '  Number of trials = ', r

      call routng ( n2, p, sn, en, m, d, l, r )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' )
     &  '  The length of the optimal connection is ', l

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Connection sequence:'
      write ( *, '(a)' ) ' '

      total = 0.0
      do i = 1, n2
        write ( *, '(2x,i8,2x,i8)' ) p(i), total
        if ( i .lt. n2 ) then
          total = total + d(p(i),p(i+1))
        end if
      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS456_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
