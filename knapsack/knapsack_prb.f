      program main

c*********************************************************************72
c
cc MAIN is the main program for KNAPSACK_PRB.
c
c  Discussion:
c
c    KNAPSACK_PRB calls sample problems for the KNAPSACK library.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 December 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'KNAPSACK_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the KNAPSACK library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
      call test05 ( )
      call test06 ( )
      call test07 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'KNAPSACK_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests MT1 for the 01 knapsack problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 December 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
c     parameter ( n = 10 )
c     parameter ( n = 5 )
c     parameter ( n = 6 )
c     parameter ( n = 7 )
      parameter ( n = 8 )
      integer jdim
      parameter ( jdim = n + 1 )

      integer c
      integer i
      integer jck
      integer mass
      integer min(jdim)
      integer p(jdim)      
      integer profit
      integer psign(jdim)
      integer w(jdim)
      integer wsign(jdim)
      integer x(jdim)
      integer xx(jdim)
      integer z
      integer zsign(jdim)

c     data p(1:n) /  92,  57,  49,  68,  60,  43,  67,  84,  87,  72 /
c     data p(1:n) /  24,  13,  23,  15,  16 /
c     data p(1:n) /  50,  50,  64,  46,  50,   5 /
c     data p(1:n) /  70,  20,  39,  37,   7,   5,  10 /
      data p(1:n) / 350, 400, 450,  20,  70,   8,   5,   5 /

c     data w(1:n) /  23,  31,  29,  44,  53,  38,  63,  85,  89,  82 /
c     data w(1:n) /  12,   7,  11,   8,   9 /
c     data w(1:n) /  56,  59,  80,  64,  75,  17 /
c     data w(1:n) /  31,  10,  20,  19,   4,   3,   6/
      data w(1:n) /  25,  35,  45,   5,  25,   3,   2,   2 /

c     c = 165
c     c = 26
c     c = 190
c     c = 50
      c = 104

      jck = 1

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  MT1 solves the 0/1 knapsack problem.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Knapsack capacity is ', c

      call knapsack_reorder ( n, p, w )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Object  Profit    Mass'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2x,i6,2x,i6)' ) 
     &  i, p(i), w(i)
      end do

      call mt1 ( n, p, w, c, z, x, jdim, jck, xx, min, psign, wsign, 
     &  zsign )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Contents of Knapsack'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Object  Profit    Mass'
      write ( *, '(a)' ) ' '
      mass = 0
      do i = 1, n
        if ( x(i) .eq. 1 ) then
          mass = mass + w(i)
          write ( *, '(2x,i6,2x,i6,2x,i6)' ) i, p(i), w(i)
        end if
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(8x,2x,i6,2x,i6)' ) z, mass

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests MT2 for the 01 knapsack problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 December 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
c     parameter ( n = 10 )
c     parameter ( n = 5 )
c     parameter ( n = 6 )
c     parameter ( n = 7 )
      parameter ( n = 8 )
      integer jdim
      parameter ( jdim = n + 3 )

      integer c
      integer i
      integer ia1(jdim)
      integer ia2(jdim)
      integer ia3(jdim)
      integer ia4(jdim)
      integer jck
      integer jfo
      integer jfs
      integer jub
      integer mass
      integer p(jdim)      
      integer profit
      real    ra(jdim)
      integer w(jdim)
      integer x(jdim)
      integer z

c     data p(1:n) /  92,  57,  49,  68,  60,  43,  67,  84,  87,  72 /
c     data p(1:n) /  24,  13,  23,  15,  16 /
c     data p(1:n) /  50,  50,  64,  46,  50,   5 /
c     data p(1:n) /  70,  20,  39,  37,   7,   5,  10 /
      data p(1:n) / 350, 400, 450,  20,  70,   8,   5,   5 /

c     data w(1:n) /  23,  31,  29,  44,  53,  38,  63,  85,  89,  82 /
c     data w(1:n) /  12,   7,  11,   8,   9 /
c     data w(1:n) /  56,  59,  80,  64,  75,  17 /
c     data w(1:n) /  31,  10,  20,  19,   4,   3,   6/
      data w(1:n) /  25,  35,  45,   5,  25,   3,   2,   2 /

c     c = 165
c     c = 26
c     c = 190
c     c = 50
      c = 104

      jfo = 1
      jfs = 1
      jck = 1

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  MT2 solves the 0/1 knapsack problem.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Knapsack capacity is ', c

      call knapsack_reorder ( n, p, w )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Object  Profit    Mass'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2x,i6,2x,i6)' ) 
     &  i, p(i), w(i)
      end do

      call mt2 ( n, p, w, c, z, x, jdim, jfo, jfs, jck, jub,
     &  ia1, ia2, ia3, ia4, ra )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Contents of Knapsack'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Object  Profit    Mass'
      write ( *, '(a)' ) ' '
      mass = 0
      do i = 1, n
        if ( x(i) .eq. 1 ) then
          mass = mass + w(i)
          write ( *, '(2x,i6,2x,i6,2x,i6)' ) i, p(i), w(i)
        end if
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(8x,2x,i6,2x,i6)' ) z, mass

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests MTC2 for the change making problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 December 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )
c
c  JDL is suggested to be set to max ( W(:) ) - 1.
c
      integer jdl
      parameter ( jdl = 10 )
      integer jdn
      parameter ( jdn = n + 1 )

      integer c
      integer i
      integer jck
      integer jfo
      integer l(jdl)
      integer m(jdl)
      integer mass
      integer pr(jdn)
      integer w(jdn)
      integer wr(jdn)
      integer x(jdn)
      integer xx(jdn)
      integer z

      data w(1:n) / 1, 4, 5, 8, 11 /

      c = 29

      jfo = 1
      jck = 1

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  MTC2 solves the change making problem.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Amount to make change for ', c

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    Coin   Value'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2x,i6)' ) i, w(i)
      end do

      call mtc2 ( n, w, c, z, x, jdn, jdl, jfo, jck, xx,wr, pr, 
     &  m, l )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    Coin   Value  Number'
      write ( *, '(a)' ) ' '
      mass = 0
      do i = 1, n
        mass = mass + x(i) * w(i)
        write ( *, '(2x,i6,2x,i6,2x,i6)' ) i, w(i), x(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(8x,2x,i6,2x,i6)' ) mass, z

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests MTG for the generalized assignment problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 December 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 2 )
      integer n
      parameter ( n = 7 )

      integer back
      integer c(10)
      integer i
      integer j
      integer jb
      integer jck
      integer mass
      integer minmax
      integer p(10,100)
      integer profit
      integer w(10,100)
      integer xstar(100)
      integer z

      data c(1:m) / 11, 22 /

      p(1,1) = 6
      p(1,2) = 9
      p(1,3) = 4
      p(1,4) = 2
      p(1,5) = 10
      p(1,6) = 3
      p(1,7) = 6
      p(2,1) = 4
      p(2,2) = 8
      p(2,3) = 9
      p(2,4) = 1
      p(2,5) = 7
      p(2,6) = 5
      p(2,7) = 4

      w(1,1) = 4
      w(1,2) = 1
      w(1,3) = 2
      w(1,4) = 1
      w(1,5) = 4
      w(1,6) = 3
      w(1,7) = 8
      w(2,1) = 9
      w(2,2) = 9
      w(2,3) = 8
      w(2,4) = 1
      w(2,5) = 3
      w(2,6) = 8
      w(2,7) = 7

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) 
     &  '  MTG solves the generalized assignment problem.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Knapsack capacities:'
      write ( *, '(a)' ) ' '
      do j = 1, m
        write ( *, '(2x,i6,2x,i6)' ) j, c(j)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Weight Matrix'
      write ( *, '(a)' ) ' '
      do i = 1, m
        write ( *, '(7(2x,i6))' ) ( w(i,j), j = 1, n )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Profit Matrix'
      write ( *, '(a)' ) ' '
      do i = 1, m
        write ( *, '(7(2x,i6))' ) ( p(i,j), j = 1, n )
      end do

      minmax = 2
      back = -1
      jck = 1

      call mtg ( n, m, p, w, c, minmax, z, xstar, back, jck, jb )

      profit = 0
      do i = 1, m
        mass = 0
        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' ) '  Knapsack ', i
        write ( *, '(a)' ) '  Object  Profit    Mass'
        do j = 1, n
          if ( xstar(j) == i ) then
            mass = mass + w(i,j)
            profit = profit + p(i,j)
            write ( *, '(2x,i6,2x,i6,2x,i6)' ) j, p(i,j), w(i,j)
          end if
        end do
        write ( *, '(2x,a,2x,6x,2x,i6)' ) 'Total ', mass
        write ( *, '(2x,a,2x,6x,2x,i6)' ) 'Limit ', c(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Total profit = ', profit

      return
      end

      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 tests MTM for the multiple knapsack problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 December 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 2 )
      integer n
      parameter ( n = 10 )

      integer back
      integer c(m)
      integer i
      integer j
      integer jck
      integer jub
      integer mass
      integer p(n)      
      integer profit
      integer w(n)
      integer x(n)
      integer z

      data c(1:m) / 70, 127 /

      data p(1:n) / 92, 57, 49, 68, 60, 43, 67, 84, 87, 72 /

      data w(1:n) / 23, 31, 29, 44, 53, 38, 63, 85, 89, 82 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) 
     &  '  MTM solves the 0/1 multiple knapsack problem.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Knapsack capacities:'
      write ( *, '(a)' ) ' '
      do j = 1, m
        write ( *, '(2x,i6,2x,i6)' ) j, c(j)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Object  Profit    Mass'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2x,i6,2x,i6)' ) 
     &  i, p(i), w(i)
      end do

      back = -1
      jck = 1

      call mtm ( n, m, p, w, c, z, x, back, jck, jub )

      profit = 0

      do j = 1, m
        mass = 0
        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' ) '  Knapsack ', j
        write ( *, '(a)' ) '  Object  Profit    Mass'
        do i = 1, n
          if ( x(i) == j ) then
            mass = mass + w(i)
            profit = profit + p(i)
            write ( *, '(2x,i6,2x,i6,2x,i6)' ) i, p(i), w(i)
          end if
        end do
        write ( *, '(2x,a,2x,6x,2x,i6)' ) 'Total ', mass
        write ( *, '(2x,a,2x,6x,2x,i6)' ) 'Limit ', c(j)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Total profit = ', profit

      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 tests MTP for the bin packing problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 December 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 9 )
      integer jdim
      parameter ( jdim = n )

      integer back
      integer bin_num
      integer c
      integer dum(jdim)
      integer fixit(jdim)
      integer i
      integer j
      integer jck
      integer kfix(jdim)
      integer lb
      integer ls(jdim)
      integer lsb(jdim)
      integer mass
      integer r(jdim)
      integer rel(jdim)
      integer res(jdim)
      integer w(n)
      integer wa(jdim)
      integer wb(jdim)
      integer wr(jdim)
      integer x(jdim)
      integer xheu(jdim)
      integer xred(jdim)
      integer xstar(n)
      integer xstarr(jdim)
      integer z

      data w(1:n) / 70, 60, 50, 33, 33, 33, 11,  7,  3 /

      c = 100

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' ) '  MTP solves the bin packing problem.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Bin capacity is ', c

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Object  Mass'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2x,i6,2x,i6)' ) 
     &  i, w(i)
      end do

      back = -1
      jck= 1

      call mtp ( n, w, c, z, xstar, jdim, back, jck, lb, wr,
     &  xstarr, dum, res, rel, x, r, wa, wb, kfix, fixit, xred, ls, 
     &  lsb, xheu )

      bin_num = 0
      do i = 1, n
        bin_num = max ( bin_num, xstar(i) )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Number of bins required is ', bin_num


      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Assignment of objects to bins:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Object     Bin'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2x,i6)' ) i, xstar(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     Bin    Mass'
      write ( *, '(a)') ' '

      do j = 1, bin_num
        mass = 0
        do i = 1, n
          if ( xstar(i) .eq. j ) then
            mass = mass + w(i)
          end if
        end do
        write ( *, '(2x,i6,2x,i6)' ) j, mass
      end do

      return
      end
      subroutine test07 ( )

c*********************************************************************72
c
cc TEST07 tests MTSL for the subset sum problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 December 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )
      integer jdd
      parameter ( jdd = 5000 )
      integer jdn
      parameter ( jdn = n + 1)
      integer itmm
      parameter ( itmm = jdn )

      integer c
      integer i
      integer ind(itmm)
      integer jck
      integer sum(itmm)
      integer td1(jdd,2)
      integer td2(jdd,2)
      integer td3(jdd,2)
      integer w(n)
      integer wo(itmm)
      integer ws(itmm)
      integer x(n)
      integer xx(itmm)
      integer z
      integer zs(itmm)

      data w(1:n) /  41, 34, 21, 20,  8,  7,  7,  4,  3,  3 /

      c = 50
      jck = 1

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07'
      write ( *, '(a)' ) '  MTSL solves the subset sum problem.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Desired subset sum is ', c

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Object    Mass'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2x,i6)' ) 
     &  i, w(i)
      end do

      call mtsl ( n, w, c, z, x, jdn, jdd, itmm, jck, wo, ind, xx, ws, 
     &  zs, sum, td1, td2, td3 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Selected subset:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Object    Mass'
      write ( *, '(a)' ) ' '
      do i = 1, n
        if ( x(i) .eq. 1 ) then
          write ( *, '(2x,i6,2x,i6)' ) i, w(i)
        end if
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(2x,a6,2x,i6)' ) ' Total', z

      return
      end
