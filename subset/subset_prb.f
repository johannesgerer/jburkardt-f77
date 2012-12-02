      program main

c*********************************************************************72
c
cc MAIN is the main program for SUBSET_PRB.
c
c  Discussion:
c
c    SUBSET_PRB calls the SUBSET test routines.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 June 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SUBSET_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the SUBSET library.'

      call test001 ( )
      call test002 ( )
      call test003 ( )
      call test9001 ( )
      call test9002 ( )
      call test9003 ( )
      call test005 ( )
      call test006 ( )
      call test007 ( )
      call test008 ( )
      call test009 ( )

      call test010 ( )
      call test011 ( )
      call test012 ( )
      call test013 ( )
      call test014 ( )
      call test015 ( )
      call test016 ( )
      call test017 ( )
      call test0174 ( )
      call test0175 ( )
      call test018 ( )
      call test019 ( )

      call test020 ( )
      call test021 ( )
      call test022 ( )
      call test023 ( )
      call test024 ( )
      call test025 ( )
      call test026 ( )
      call test027 ( )
      call test028 ( )
      call test029 ( )

      call test021a ( )
      call test022a ( )
      call test023a ( )
      call test024a ( )
      call test025a ( )
      call test026a ( )
      call test026b ( )
      call test026c ( )
      call test026d ( )
      call test027a ( )
      call test028a ( )
      call test029a ( )
      call test0295 ( )

      call test030 ( )
      call test0304 ( )
      call test0305 ( )
      call test031 ( )
      call test032 ( )
      call test0321 ( )
      call test0322 ( )
      call test03225 ( )
      call test0323 ( )
      call test0324 ( )
      call test0325 ( )
      call test0327 ( )

      call test058 ( )
      call test059 ( )

      call test060 ( )
      call test061 ( )
      call test0615 ( )
      call test062 ( )
      call test06225 ( )

      call test033 ( )
      call test034 ( )
      call test0625 ( )
      call test035 ( )
      call test0364 ( )
      call test0627 ( )
      call test036 ( )
      call test037 ( )
      call test038 ( )
      call test039 ( )

      call test040 ( )
      call test041 ( )
      call test042 ( )
      call test047 ( )
      call test043 ( )
      call test044 ( )
      call test045 ( )
      call test046 ( )
      call test048 ( )
      call test049 ( )

      call test050 ( )
      call test051 ( )
      call test052 ( )
      call test053 ( )
      call test054 ( )
      call test055 ( )
      call test056 ( )
      call test057 ( )

      call test063 ( )
      call test064 ( )
      call test065 ( )
      call test066 ( )
      call test067 ( )
      call test0675 ( )
      call test068 ( )
      call test0683 ( )
      call test0685 ( )
      call test0686 ( )
      call test0687 ( )
      call test0688 ( )
      call test0689 ( )
      call test06895 ( )
      call test069 ( )

      call test070 ( )
      call test071 ( )
      call test072 ( )
      call test073 ( )
      call test074 ( )
      call test075 ( )
      call test076 ( )
      call test077 ( )
      call test0771 ( )
      call test07715 ( )
      call test0772 ( )
      call test0773 ( )
      call test078 ( )
      call test079 ( )

      call test080 ( )
      call test081 ( )
      call test082 ( )
      call test083 ( )
      call test0835 ( )
      call test084 ( )
      call test085 ( )
      call test093 ( )
      call test086 ( )
      call test087 ( )
      call test088 ( )
      call test089 ( )
 
      call test090 ( )
      call test091 ( )
      call test092 ( )
      call test094 ( )
      call test095 ( )
      call test0955 ( )
      call test096 ( )
      call test097 ( )
      call test098 ( )
      call test099 ( )

      call test100 ( )
      call test101 ( )
      call test102 ( )
      call test103 ( )
      call test104 ( )
      call test105 ( )
      call test106 ( )
      call test107 ( )
      call test108 ( )
      call test1085 ( )
      call test109 ( )

      call test110 ( )
      call test111 ( )
      call test112 ( )
      call test113 ( )
      call test114 ( )
      call test115 ( )
      call test116 ( )
      call test1163 ( )
      call test1165 ( )
      call test117 ( )
      call test118 ( )
      call test119 ( )

      call test120 ( )
      call test121 ( )
      call test1215 ( )
      call test122 ( )
      call test123 ( )
      call test124 ( )
      call test1245 ( )
      call test125 ( )
      call test126 ( )
      call test127 ( )
      call test128 ( )
      call test129 ( )

      call test130 ( )
      call test131 ( )
      call test132 ( )
      call test133 ( )
      call test134 ( )
      call test135 ( )
      call test136 ( )
      call test137 ( )
      call test138 ( )
      call test139 ( )
      call test1395 ( )

      call test140 ( )
      call test141 ( )
      call test142 ( )
      call test143 ( )
      call test1435 ( )
      call test144 ( )
      call test145 ( )
      call test146 ( )
      call test147 ( )
      call test1475 ( )
      call test1476 ( )
      call test1477 ( )
      call test148 ( )
      call test149 ( )

      call test150 ( )
      call test151 ( )
      call test152 ( )
      call test153 ( )
      call test1531 ( )
      call test0626 ( )
      call test154 ( )
      call test155 ( )
      call test156 ( )
      call test1565 ( )
      call test1566 ( )
      call test1567 ( )
      call test1568 ( )
      call test1569 ( )
      call test15695 ( )
      call test15696 ( )
      call test15698 ( )
      call test157 ( )
      call test158 ( )
      call test159 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SUBSET_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test001 ( )

c*********************************************************************72
c
cc TEST001 tests ASM_ENUM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 7 )

      integer asm_num
      integer n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST001'
      write ( *, '(a)' ) 
     &  '  ASM_ENUM returns the number of alternating sign'
      write ( *, '(a)' ) '  matrices of a given order.'

      write ( *, '(a)' ) ' '
      do n = 0, n_max
        call asm_enum ( n, asm_num )
        write ( *, '(2x,i2,2x,i8)' ) n, asm_num
      end do

      return
      end
      subroutine test002 ( )

c*********************************************************************72
c
cc TEST002 tests ASM_TRIANGLE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 7 )

      integer a(n_max+1)
      integer i
      integer n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST002'
      write ( *, '(a)' ) '  ASM_TRIANGLE returns a row of the'
      write ( *, '(a)' ) '  alternating sign matrix triangle.'
      write ( *, '(a)' ) ' '

      do n = 0, n_max
        call asm_triangle ( n, a )
        write ( *, '(2x,i2,2x,8i8)' ) n, ( a(i), i = 1, n+1 )
      end do

      return
      end
      subroutine test003 ( )

c*********************************************************************72
c
cc TEST003 tests BELL and BELL_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer c
      integer c2(0:10)
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST003'
      write ( *, '(a)' ) '  BELL computes Bell numbers.'
      write ( *, '(a)' ) '  BELL_VALUES returns some exact values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  N  exact C(I)  computed C(I)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bell_values ( n_data, n, c )

        if ( n_data .eq. 0 ) then
          go to 20
        end if
    
        call bell ( n, c2 )

        write ( *, '(2x,i4,2i10)' ) n, c, c2(n)

      go to 10

20    continue

      return
      end
      subroutine test9001 ( )

c*********************************************************************72
c
cc TEST9001 tests BVEC_ADD and BVEC_SUB;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer bvec1(n)
      integer bvec2(n)
      integer bvec3(n)
      integer bvec4(n)
      integer i
      integer i4_uniform
      integer j
      integer k
      integer l
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST9001'
      write ( *, '(a)' ) 
     &  '  BVEC_ADD adds binary vectors representing integers;'
      write ( *, '(a)' ) 
     &  '  BVEC_SUB subtracts binary vectors representing integers;'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '        I        J        K = I + J    L = I - J'
      write ( *, '(a)' ) ' '

      do test = 1, test_num
        
        i = i4_uniform ( -100, 100, seed )
        j = i4_uniform ( -100, 100, seed )

        write ( *, '(a)' ) ' '

        write ( *, '(2x,i8,2x,i8)' ) i, j

        k = i + j
        l = i - j

        write ( *, '(a20,2x,i8,2x,i8)' ) '  Directly:         ', k, l

        call i4_to_bvec ( i, n, bvec1 )
        call i4_to_bvec ( j, n, bvec2 )

        call bvec_add ( n, bvec1, bvec2, bvec3 )
        call bvec_to_i4 ( n, bvec3, k )

        call bvec_sub ( n, bvec1, bvec2, bvec4 )
        call bvec_to_i4 ( n, bvec4, l )

        write ( *, '(a20,2x,i8,2x,i8)' ) '  BVEC_ADD, BVEC_SUB', k, l

      end do

      return
      end
      subroutine test9002 ( )

c*********************************************************************72
c
cc TEST9002 tests BVEC_COMPLEMENT2;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer bvec1(n)
      integer bvec2(n)
      integer i
      integer j
      integer i4_uniform
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 5 )

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST9002'
      write ( *, '(a)' ) 
     &  '  BVEC_COMPLEMENT2 returns the two''s complement'
      write ( *, '(a)' ) '  of a (signed) binary vector;'
      write ( *, '(a)' ) ' '

      do test = 1, test_num
    
        i = i4_uniform ( -100, 100, seed )

        call i4_to_bvec ( i, n, bvec1 )

        call bvec_complement2 ( n, bvec1, bvec2 )

        call bvec_to_i4 ( n, bvec2, j )

        write ( *, '(a)' ) ' '
        write ( *, '(a,2x,i8)' ) '  I = ', i
        write ( *, '(a,2x,i8)' ) '  J = ', j
        call bvec_print ( n, bvec1, ' ' )
        call bvec_print ( n, bvec2, ' ' )

      end do

      return
      end
      subroutine test9003 ( )

c*********************************************************************72
c
cc TEST9003 tests BVEC_MUL;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 15 )

      integer bvec1(n)
      integer bvec2(n)
      integer bvec3(n)
      integer i
      integer i4_uniform
      integer j
      integer k
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST9003'
      write ( *, '(a)' ) '  BVEC_MUL multiplies binary vectors '
      write ( *, '(a)' ) '  representing integers;'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        I        J        K = I * J'
      write ( *, '(a)' ) ' '

      do test = 1, test_num
    
        i = i4_uniform ( -100, 100, seed )
        j = i4_uniform ( -100, 100, seed )

        write ( *, '(a)' ) ' '

        write ( *, '(2x,i8,2x,i8)' ) i, j

        k = i * j

        write ( *, '(a20,2x,i8)' ) '  Directly:         ', k

        call i4_to_bvec ( i, n, bvec1 )
        call i4_to_bvec ( j, n, bvec2 )
        call bvec_mul ( n, bvec1, bvec2, bvec3 )
        call bvec_to_i4 ( n, bvec3, k )

        write ( *, '(a20,2x,i8)' ) '  BVEC_MUL          ', k

      end do
    
      return
      end
      subroutine test005 ( )

c*********************************************************************72
c
cc TEST005 tests CATALAN and CATALAN_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer c
      integer c2(11)
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST005'
      write ( *, '(a)' ) '  CATALAN computes Catalan numbers.'
      write ( *, '(a)' ) '  CATALAN_VALUES returns some exact values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  N  exact C(I)  computed C(I)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call catalan_values ( n_data, n, c )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        call catalan ( n, c2 )

        write ( *, '(2x,i4,2i8)' ) n, c, c2(n+1)

      go to 10

20    continue

      return
      end
      subroutine test006 ( )

c*********************************************************************72
c
cc TEST006 tests CATALAN_ROW_NEXT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer c(0:n)
      integer i
      integer ido
      integer j

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST006'
      write ( *, '(a)' ) 
     &  '  CATALAN_ROW_NEXT computes a row of Catalan''s triangle.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  First, compute row 7:'

      ido = 0
      i = 7
      call catalan_row_next ( ido, i, c )
      write ( *, '(2x,i2,2x,11i6)' ) i,  ( c(j), j = 0, i )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Now compute rows one at a time:'
      write ( *, '(a)' ) ' '

      ido = 0

      do i = 0, n
        call catalan_row_next ( ido, i, c )
        ido = 1
        write ( *, '(2x,i2,2x,11i6)' ) i, ( c(j), j = 0, i )
      end do

      return
      end
      subroutine test007 ( )

c*********************************************************************72
c
cc TEST007 tests CFRAC_TO_RAT and RAT_TO_CFRAC.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Discussion:
c
c    Compute the continued fraction form of 4096/15625.
c
c  Modified:
c
c    09 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 10 )

      integer a(m)
      integer bot
      integer i
      integer ierror
      integer n
      integer p(m)
      integer q(m)
      integer top

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST007'
      write ( *, '(a)' ) 
     &  '  RAT_TO_CFRAC fraction => continued fraction,'
      write ( *, '(a)' ) 
     &  '  CFRAC_TO_RAT continued fraction => fraction.'
      write ( *, '(a)' ) ' '
      top = 4096
      bot = 15625
      write ( *, '(a,i8,a,i8)' ) 
     &  '  Regular fraction is ', top, ' / ', bot
 
      call rat_to_cfrac ( top, bot, m, n, a, ierror )
 
      call i4vec_print ( n, a, '  Continued fraction coefficients:' )

      call cfrac_to_rat ( n, a, p, q )
     
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The continued fraction convergents.'
      write ( *, '(a)' ) 
     &  '  The last row contains the value of the continued'
      write ( *, '(a)' ) '  fraction, written as a common fraction.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  I, P(I), Q(I), P(I)/Q(I)'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,i3,2i8,g14.6)' ) i, p(i), q(i), 
     &    dble ( p(i) ) / dble ( q(i) )
      end do
 
      return
      end
      subroutine test008 ( )

c*********************************************************************72
c
cc TEST008 tests CFRAC_TO_RFRAC and RFRAC_TO_CFRAC.
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
      implicit none

      integer maxm
      parameter ( maxm = 10 )

      double precision g(2*maxm)
      double precision h(2*maxm)
      integer i
      integer ierror
      integer m
      double precision p(maxm)
      double precision q(maxm+1)

      m = 3

      p(1) = 1.0D+00
      p(2) = 1.0D+00
      p(3) = 2.0D+00

      q(1) = 1.0D+00
      q(2) = 3.0D+00
      q(3) = 1.0D+00
      q(4) = 1.0D+00
    
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST008'
      write ( *, '(a)' ) 
     &  '  CFRAC_TO_RFRAC: continued fraction to ratio;'
      write ( *, '(a)' ) 
     &  '  RFRAC_TO_CFRAC: ratio to continued fration.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Rational polynomial fraction coefficients:'
      write ( *, '(a)' ) ' '
      write ( *, '(a,5f12.4)' ) '  P:  ', ( p(i), i = 1, m )
      write ( *, '(a,5f12.4)' ) '  Q:  ', ( q(i), i = 1, m+1 )
 
      call rfrac_to_cfrac ( m, p, q, h, ierror )
 
      call r8vec_print ( 2*m, h, '  Continued fraction coefficients:' )

      do i = 1, 2*m
        g(i) = 1.0D+00
      end do

      call cfrac_to_rfrac ( 2*m, g, h, p, q )
 
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Recovered rational polynomial:'
      write ( *, '(a)' ) ' '
      write ( *, '(a,5f12.4)' ) '  P:  ', ( p(i), i = 1, m )
      write ( *, '(a,5f12.4)' ) '  Q:  ', ( q(i), i = 1, m+1 )
 
      return
      end
      subroutine test009 ( )

c*********************************************************************72
c
cc TEST009 tests CHANGE_GREEDY.
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
      implicit none

      integer coin_num
      parameter ( coin_num = 6 )

      integer change(100)
      integer change_num
      integer coin_value(coin_num)
      integer i
      integer j
      integer total
      integer total2

      save coin_value

      data coin_value / 1, 5, 10, 25, 50, 100 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST009'
      write ( *, '(a)' ) '  CHANGE_GREEDY makes change using the'
      write ( *, '(a)' ) '  biggest coins first.'

      total = 73

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 
     &  '  The total for which change is to be made: ', total
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The available coins are:'
      write ( *, '(a)' ) ' '
      do i = 1, coin_num
        write ( *, '(2x,i8)' ) coin_value(i)
      end do

      call change_greedy ( total, coin_num, coin_value, change_num,
     &  change )

      write ( *, '(a)' ) ' '
      write ( *, '(2x,i8,4x,(20i3))' ) 
     &  change_num, ( change(j), j = 1, change_num )
      total2 = 0
      do j = 1, change_num
        total2 = total2 + coin_value(change(j))
      end do
      write ( *, '(2x,i8,4x,(20i3))' ) 
     &  total2, ( coin_value(change(j)), j = 1, change_num )

      return
      end
      subroutine test010 ( )

c*********************************************************************72
c
cc TEST010 tests CHANGE_NEXT.
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
      implicit none

      integer coin_num
      parameter ( coin_num = 6 )

      integer change(100)
      integer change_num
      integer coin_value(coin_num)
      logical done
      integer i
      integer j
      integer total

      save coin_value

      data coin_value / 1, 5, 10, 25, 50, 100 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST010'
      write ( *, '(a)' ) '  CHANGE_NEXT displays the next possible'
      write ( *, '(a)' ) '  way to make change for a given total'

      total = 50

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 
     &  '  The total for which change is to be made: ', total
      write ( *, '(a)' ) ' '

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The available coins are:'
      write ( *, '(a)' ) ' '
      do i = 1, coin_num
        write ( *, '(2x,i8)' ) coin_value(i)
      end do

      done = .true.
      i = 0
    
10    continue

        call change_next ( total, coin_num, coin_value, change_num, 
     &    change, done )

        if ( done .or. 9 .lt. i ) then
          go to 20
        end if

        i = i + 1
        write ( *, '(a)' ) ' '
        write ( *, '(i3, ":")' ) i
        write ( *, '(2x,25i3)' ) 
     &  ( coin_value(change(j)), j = 1, change_num )

      go to 10

20    continue

      return
      end
      subroutine test011 ( )

c*********************************************************************72
c
cc TEST011 tests COMB_NEXT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 April 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      integer a(n)
      logical done
      integer i
      integer k

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST011'
      write ( *, '(a)' ) '  COMB_NEXT produces combinations.'
      write ( *, '(a)' ) ' '
 
      do k = 1, n

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Combinations of size K = ', k
        write ( *, '(a)' ) ' '

        done = .true.

10      continue

          call comb_next ( n, k, a, done )
 
          if ( done ) then
            go to 20
          end if

          write ( *, '(2x,5i3)' ) ( a(i), i = 1, k )

        go to 10

20      continue

      end do

      return
      end
      subroutine test012 ( )

c*********************************************************************72
c
cc TEST012 tests COMB_ROW.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    16 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer c(0:n)
      integer i
      integer ido
      integer j

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST012'
      write ( *, '(a)' ) 
     &  '  COMB_ROW computes a row of Pascal''s triangle.'
      write ( *, '(a)' ) ' '
 
      ido = 0
 
      do i = 0, n
        call comb_row ( ido, i, c )
        ido = 1
        write ( *, '(2x,i2,2x,11i5)' ) i, ( c(j), j = 0, i )
      end do
 
      return
      end
      subroutine test013 ( )

c*********************************************************************72
c
cc TEST013 tests COMB_UNRANK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    16 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      integer a(n)
      integer cnk
      integer i
      integer m
      parameter ( m = 10 )
      integer rank

      call combin2 ( m, n, cnk )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST013'
      write ( *, '(a)' ) 
     &  '  COMB_UNRANK returns a combination of N things'
      write ( *, '(a)' ) '  out of M, given the lexicographic rank.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The total set size is M = ', m
      write ( *, '(a,i8)' ) '  The subset size is N =    ', n
      write ( *, '(a,i8)' ) 
     &  '  The number of combinations of N out of M is ', cnk
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   Rank	  Combination'
      write ( *, '(a)' ) ' '
 
      do rank = 1, 3
        call comb_unrank ( m, n, rank, a )
        write ( *, '(2x,i3,3x,5i4)' ) rank, ( a(i), i = 1, n )
      end do
 
      do rank = 6, 8
        call comb_unrank ( m, n, rank, a )
        write ( *, '(2x,i3,3x,5i4)' ) rank, ( a(i), i = 1, n )
      end do
 
      do rank = 250, 252
        call comb_unrank ( m, n, rank, a )
        write ( *, '(2x,i3,3x,5i4)' ) rank, ( a(i), i = 1, n )
      end do
 
      return
      end
      subroutine test014 ( )

c*********************************************************************72
c
cc TEST014 tests COMBIN.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    16 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision cnk
      integer k
      integer n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST014'
      write ( *, '(a)' ) '  COMBIN evaluates C(N,K).'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N         K    CNK'
      write ( *, '(a)' ) ' '

      do n = 0, 4
        do k = 0, n
          call combin ( n, k, cnk )
          write ( *, '(2x,i8,2x,i8,g14.6)' ) n, k, cnk
        end do
      end do

      return
      end
      subroutine test015 ( )

c*********************************************************************72
c
cc TEST015 tests COMBIN2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer cnk
      integer k
      integer n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST015'
      write ( *, '(a)' ) '  COMBIN2 evaluates C(N,K).'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N       K      CNK'
      write ( *, '(a)' ) ' '

      do n = 0, 4
        do k = 0, n
          call combin2 ( n, k, cnk )
          write ( *, '(2x,i8,i8,i8)' ) n, k, cnk
        end do
      end do

      return
      end
      subroutine test016 ( )

c*********************************************************************72
c
cc TEST016 tests COMP_NEXT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 July 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer k
      parameter ( k = 3 )

      integer a(k)
      integer h
      integer i
      logical more
      integer n
      integer t

      n = 6
      more = .false.

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST016'
      write ( *, '(a)' ) '  COMP_NEXT generates compositions.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Seeking all compositions of N = ', n
      write ( *, '(a,i8,a)' ) '  using ', k, ' parts.'
      write ( *, '(a)' ) ' '

10    continue

      call comp_next ( n, k, a, more, h, t )

      write ( *, '(2x,8i4)' ) ( a(i), i = 1, k )

      if ( more )  then
        go to 10
      end if
 
      return
      end
      subroutine test017 ( )

c*********************************************************************72
c
cc TEST017 tests COMP_RANDOM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer k
      parameter ( k = 5 )

      integer a(k)
      integer i
      integer j
      integer n
      integer seed
      integer test_num
      parameter ( test_num = 5 )

      n = 10
      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST017'
      write ( *, '(a)' ) '  COMP_RANDOM generates random compositions.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Seeking random compositions of N = ', n
      write ( *, '(a,i8,a)' ) '  using ', k, ' parts.'
      write ( *, '(a)' ) ' '

      do i = 1, test_num
        call comp_random ( n, k, seed, a )
        write ( *, '(2x,8i4)' ) ( a(j), j = 1, k )
      end do
 
      return
      end
      subroutine test0174 ( )

c*********************************************************************72
c
cc TEST0174 tests COMPNZ_NEXT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer k
      parameter ( k = 3 )

      integer a(k)
      integer i
      integer j
      logical more
      integer n
      integer number

      more = .false.
      n = 6

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0174'
      write ( *, '(a)' ) 
     &  '  COMPNZ_NEXT generates compositions with nonzero parts.'
      write ( *, '(a)' ) '  COMPNZ_ENUM enumerates them.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Seeking all compositions of N = ', n
      write ( *, '(a,i8,a)' ) '  using ', k, ' nonzero parts.'
    
      call compnz_enum ( n, k, number )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 
     &  '  The number of these compositions is ', number
      write ( *, '(a)' ) ' '

      i = 0

10    continue

        call compnz_next ( n, k, a, more )

        i = i + 1

        write ( *, '(2x,i4,2x,8i4)' ) i, ( a(j), j = 1, k )

        if ( .not. more )  then
          go to 20
        end if

      go to 10

20    continue
 
      return
      end
      subroutine test0175 ( )

c*********************************************************************72
c
cc TEST0175 tests COMPNZ_RANDOM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer k
      parameter ( k = 5 )

      integer a(k)
      integer i
      integer j
      integer n
      integer seed
      integer test_num
      parameter ( test_num = 5 )

      n = 10
      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0175'
      write ( *, '(a)' ) '  COMPNZ_RANDOM generates random'
      write ( *, '(a)' ) '  compositions with nonzero parts.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 
     &  '  Seeking random compositions of N = ', n
      write ( *, '(a,i8,a)' ) '  using ', k, ' nonzero parts.'
      write ( *, '(a)' ) ' '

      do i = 1, test_num
        call compnz_random ( n, k, seed, a )
        write ( *, '(2x,8i4)' ) ( a(j), j = 1, k )
      end do
 
      return
      end
      subroutine test018 ( )

c*********************************************************************72
c
cc TEST018 tests CONGRUENCE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer test_num
      parameter ( test_num = 21 )

      integer a
      integer a_test(test_num)
      integer b
      integer b_test(test_num)
      integer c
      integer c_test(test_num)
      integer i4_modp
      integer ierror
      integer result
      integer test
      integer x

      save a_test
      save b_test
      save c_test

      data a_test /
     &   1027,  1027,  1027,  1027, -1027, 
     &  -1027, -1027, -1027,     6,     0, 
     &      0,     0,     1,     1,     1, 
     &   1024,     0,     0,     5,     2, 
     &      7 /
      data b_test /
     &    712,   712,  -712,  -712,   712, 
     &    712,  -712,  -712,     8,     0, 
     &      1,     1,     0,     0,     1, 
     & -15625,     0,     3,     0,     4, 
     &     19 /
      data c_test /
     &      7,    -7,     7,    -7,     7, 
     &     -7,     7,    -7,    50,     0, 
     &      0,     1,     0,     1,     0, 
     &  11529,     1,    11,    19,     7, 
     &      1 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST018'
      write ( *, '(a)' ) '  CONGRUENCE solves a congruence equation:'
      write ( *, '(a)' ) '    A * X = C mod ( B )'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) 
     &  '   I        A         B         C',
     &  '           X     Mod ( A*X-C,B)'
      write ( *, '(a)' ) ' '

      do test = 1, test_num

        a = a_test(test)
        b = b_test(test)
        c = c_test(test)

        call congruence ( a, b, c, ierror, x )

        if ( ierror .ne. 0 ) then
          write ( *, '(2x,i2,2x,3i10,a,i10)' ) test, a, b, c, 
     &      ' Error code = ', ierror
        else
          if ( b .ne. 0 ) then
            result = i4_modp ( a * x - c, b )
          else
            result = 0
          end if
          write ( *, '(2x,i2,2x,5i10)' ) test, a, b, c, x, result
        end if

      end do

      return
      end
      subroutine test019 ( )

c*********************************************************************72
c
cc TEST019 tests COUNT_POSE_RANDOM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer blocks(6)
      integer goal
      integer i
      integer j 
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST019'
      write ( *, '(a)' ) '  COUNT_POSE_RANDOM poses a random problem'
      write ( *, '(a)' ) '  for the game "The Count is Good".'

      seed = 123456789

      do i = 1, 5

        call count_pose_random ( seed, blocks, goal )

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Problem #', i
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '    The goal number:'
        write ( *, '(a)' ) ' '
        write ( *, '(6x,i8)' ) goal
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '    The available numbers are '
        write ( *, '(a)' ) ' '
        write ( *, '(6x,6i4)' ) ( blocks(j), j = 1, 6 )

      end do

      return
      end
      subroutine test020 ( )

c*********************************************************************72
c
cc TEST020 tests R8_TO_CFRAC.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer a(0:n)
      double precision error
      integer i
      integer p(-1:n)
      integer q(-1:n)
      double precision r
      double precision r8_pi
      double precision temp

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST020'
      write ( *, '(a)' ) '  R8_TO_CFRAC converts a real number to'
      write ( *, '(a)' ) '  a sequence of continued fraction '
      write ( *, '(a)' ) '  convergents.'

      r = 2.0D+00 * r8_pi ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Use the real number R = ', r

      call r8_to_cfrac ( r, n, a, p, q )

      write ( *, '(a)' ) ' '

      do i = 0, n
        temp = dble ( p(i) ) / dble ( q(i) )
        error = r - temp
        write ( *, '(2x,3i12,2g14.6)' ) a(i), p(i), q(i), temp, error
      end do

      return
      end
      subroutine test021 ( )

c*********************************************************************72
c
cc TEST021 tests DEBRUIJN.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer test_num
      parameter ( test_num = 3 )

      integer i
      integer m
      integer m_test(test_num)
      integer n
      integer n_test(test_num)
      integer string(27)
      integer test
    
      save m_test
      save n_test

      data m_test / 2, 3, 2 /
      data n_test / 3, 3, 4 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST021'
      write ( *, '(a)' ) '  DEBRUIJN computes a de Bruijn string.'

      do test = 1, test_num

        m = m_test(test)
        n = n_test(test)

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  The alphabet size is M = ', m
        write ( *, '(a,i8)' ) '  The string length is N = ', n

        call debruijn ( m, n, string )

        write ( *, '(a)' ) ' '
        write ( *, '(4x,80i1)' ) ( string(i), i = 1, m**n )

      end do

      return
      end
      subroutine test022 ( )

c*********************************************************************72
c
cc TEST022 tests DEC_ADD.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer abot
      integer atop
      integer bbot
      integer btop
      integer cbot
      integer ctop
      integer dec_digit
      integer ierror
      character * ( 15 ) string
      integer string_length

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST022'
      write ( *, '(a)' ) '  DEC_ADD adds two decimals.'

      dec_digit = 3

      atop = 128
      abot = - 1
      btop = 438
      bbot = - 2

      call dec_add ( atop, abot, btop, bbot, dec_digit, ctop, 
     &  cbot, ierror )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of decimal places is ', dec_digit
      write ( *, '(a)' ) ' '

      call dec_to_s ( atop, abot, string )
      string_length = len_trim ( string )
      write ( *, '(a,a)' ) '  A = ', string(1:string_length)

      call dec_to_s ( btop, bbot, string )
      string_length = len_trim ( string )
      write ( *, '(a,a)' ) '  B = ', string(1:string_length)

      call dec_to_s ( ctop, cbot, string )
      string_length = len_trim ( string )
      write ( *, '(a,a)' ) '  C = ', string(1:string_length)
 
      return
      end
      subroutine test023 ( )

c*********************************************************************72
c
cc TEST023 tests DEC_DIV.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer abot
      integer atop
      integer bbot
      integer btop
      integer cbot
      integer ctop
      integer dec_digit
      integer ierror
      character * ( 15 ) string
      integer string_length

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST023'
      write ( *, '(a)' ) '  DEC_DIV divides two decimals.'

      dec_digit = 3

      atop = 523
      abot = -1
      btop = 134
      bbot = 2

      call dec_div ( atop, abot, btop, bbot, dec_digit, ctop, cbot, 
     &  ierror )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of decimal places is ', dec_digit
      write ( *, '(a)' ) ' '

      call dec_to_s ( atop, abot, string )
      string_length = len_trim ( string )
      write ( *, '(a,a)' ) '  A = ', string(1:string_length)

      call dec_to_s ( btop, bbot, string )
      string_length = len_trim ( string )
      write ( *, '(a,a)' ) '  B = ', string(1:string_length)

      call dec_to_s ( ctop, cbot, string )
      string_length = len_trim ( string )
      write ( *, '(a,a)' ) '  C = ', string(1:string_length)
     
      return
      end
      subroutine test024 ( )

c*********************************************************************72
c
cc TEST024 tests DEC_MUL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer abot
      integer atop
      integer bbot
      integer btop
      integer cbot
      integer ctop
      integer dec_digit
      character * ( 15 ) string
      integer string_length

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST024'
      write ( *, '(a)' ) '  DEC_MUL multiplies two decimals.'

      dec_digit = 2

      atop = 14
      abot = - 4
      btop = 16
      bbot = 2

      call dec_mul ( atop, abot, btop, bbot, dec_digit, ctop, cbot )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 
     &  '  Number of decimal places is ', dec_digit
      write ( *, '(a)' ) ' '

      call dec_to_s ( atop, abot, string )
      string_length = len_trim ( string )
      write ( *, '(a,a)' ) '  A = ', string(1:string_length)

      call dec_to_s ( btop, bbot, string )
      string_length = len_trim ( string )
      write ( *, '(a,a)' ) '  B = ', string(1:string_length)

      call dec_to_s ( ctop, cbot, string )
      string_length = len_trim ( string )
      write ( *, '(a,a)' ) '  C = ', string(1:string_length)
     
      return
      end
      subroutine test025 ( )

c*********************************************************************72
c
cc TEST025 tests DEC_ROUND.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer test_num
      parameter ( test_num = 7 )

      integer d_test(test_num)
      integer dec_digit
      integer exponent
      integer exponent_test(test_num)
      integer mantissa
      integer mantissa_test(test_num)
      integer test

      save d_test
      save exponent_test
      save mantissa_test

      data d_test / 1, 2, 3, 4, 2, 3, 4 /
      data exponent_test / -1,  -1, -1, -1, 2, 2, 2 /
      data mantissa_test / 523, 523, 523, 523, 6340, 6340, 6340 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST025'
      write ( *, '(a)' ) 
     &  '  DEC_ROUND "rounds" a decimal to a number of digits.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '           -----Before-------  -----After--------'
      write ( *, '(a)' ) 
     &  '  Digits   Mantissa  Exponent  Mantissa  Exponent'
      write ( *, '(a)' ) ' '

      do test = 1, test_num

        dec_digit = d_test(test)

        mantissa = mantissa_test(test)
        exponent = exponent_test(test)

        call dec_round ( mantissa, exponent, dec_digit, mantissa, 
     &    exponent )

        write ( *, '(2x,3i8,4x,2i8)' ) 
     &    dec_digit, mantissa_test(test), exponent_test(test), 
     &    mantissa, exponent

      end do

      return
      end
      subroutine test026 ( )

c*********************************************************************72
c
cc TEST026 tests DEC_TO_RAT and RAT_TO_DEC.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer exponent
      integer i
      integer i4_uniform
      integer mantissa
      double precision r1
      double precision r2
      double precision r3
      integer rat_bot
      integer rat_bot2
      integer rat_top
      integer rat_top2
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST026'
      write ( *, '(a)' ) '  RAT_TO_DEC fraction => decimal,'
      write ( *, '(a)' ) '  DEC_TO_RAT decimal => fraction.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  In this test, choose the top and bottom'
      write ( *, '(a)' ) '  of a rational at random, and compute the'
      write ( *, '(a)' ) '  equivalent real number.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Then convert to decimal, and the equivalent real.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Then convert back to rational and the equivalent real.'
  
      seed = 123456789

      do i = 1, 10

        rat_top = i4_uniform ( -1000, 1000, seed )

        rat_bot = i4_uniform (     1, 1000, seed )

        r1 = dble ( rat_top ) / dble ( rat_bot )

        call rat_to_dec ( rat_top, rat_bot, mantissa, exponent )

        r2 = dble ( mantissa ) * 10.0D+00**( exponent )
 
        call dec_to_rat ( mantissa, exponent, rat_top2, rat_bot2 )
        r3 = dble ( rat_top2 ) / dble ( rat_bot2 )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f10.6,a,i12,a,i12)' ) 
     &    r1, '=', rat_top, '/', rat_bot
        write ( *, '(2x,f10.6,a,i12,a,i12)' ) 
     &    r2, '=', mantissa, '*10^', exponent
        write ( *, '(2x,f10.6,a,i12,a,i12)' ) 
     &    r3, '=', rat_top2, '/', rat_bot2

      end do
 
      return
      end
      subroutine test027 ( )

c*********************************************************************72
c
cc TEST027 tests DEC_TO_S.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 December 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer exponent
      integer mantissa
      character * ( 100 ) s
      integer s_length

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST027'
      write ( *, '(a)' ) '  DEC_TO_S prints a decimal value.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Mantissa  Exponent  String'
      write ( *, '(a)' ) ' '

      mantissa = 523
      exponent = -1
      call dec_to_s ( mantissa, exponent, s )
      s_length = len_trim ( s )
      write ( *, '(2x,i8,2x,i8,2x,a)' ) 
     &  mantissa, exponent, s(1:s_length)

      mantissa = 134
      exponent = 2
      call dec_to_s ( mantissa, exponent, s )
      s_length = len_trim ( s )
      write ( *, '(2x,i8,2x,i8,2x,a)' ) 
     &  mantissa, exponent, s(1:s_length)

      mantissa = -134
      exponent = 2
      call dec_to_s ( mantissa, exponent, s )
      s_length = len_trim ( s )
      write ( *, '(2x,i8,2x,i8,2x,a)' ) 
     &  mantissa, exponent, s(1:s_length)

      mantissa = 0
      exponent = 10
      call dec_to_s ( mantissa, exponent, s )
      s_length = len_trim ( s )
      write ( *, '(2x,i8,2x,i8,2x,a)' ) 
     &  mantissa, exponent, s(1:s_length)

      do exponent = -8, 3
        mantissa = 123456
        call dec_to_s ( mantissa, exponent, s )
        s_length = len_trim ( s )
        write ( *, '(2x,i8,2x,i8,2x,a)' ) 
     &    mantissa, exponent, s(1:s_length)
      end do

      return
      end
      subroutine test028 ( )

c*********************************************************************72
c
cc TEST028 tests DEC_WIDTH.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dec_width
      integer exponent
      integer i
      integer mantissa

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST028'
      write ( *, '(a)' ) 
     &  '  DEC_WIDTH determines the "width" of a decimal.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Mantissa  Exponent  Width'
      write ( *, '(a)' ) ' '

      mantissa = 523
      exponent = -1
      i = dec_width ( mantissa, exponent )
      write ( *, '(2x,i8,2x,i8,2x,i8)' ) mantissa, exponent, i

      mantissa = 134
      exponent = 2
      i = dec_width ( mantissa, exponent )
      write ( *, '(2x,i8,2x,i8,2x,i8)' ) mantissa, exponent, i

      mantissa = -134
      exponent = 2
      i = dec_width ( mantissa, exponent )
      write ( *, '(2x,i8,2x,i8,2x,i8)' ) mantissa, exponent, i

      mantissa = 0
      exponent = 10
      i = dec_width ( mantissa, exponent )
      write ( *, '(2x,i8,2x,i8,2x,i8)' ) mantissa, exponent, i

      do exponent = -8, 3
        mantissa = 123456
        i = dec_width ( mantissa, exponent )
        write ( *, '(2x,i8,2x,i8,2x,i8)' ) mantissa, exponent, i
      end do

      return
      end
      subroutine test029 ( )

c*********************************************************************72
c
cc TEST029 tests DECMAT_DET.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n3
      parameter ( n3 = 3 )
      integer n4
      parameter ( n4 = 4 )

      integer dec_digit
      integer i
      integer a3(n3,n3)
      integer a4(n4,n4)
      integer b3(n3,n3)
      integer b4(n4,n4)
      integer dbot
      integer dtop
      integer ierror
      integer j
      integer k
      double precision r

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST029'
      write ( *, '(a)' ) 
     &  '  DECMAT_DET: determinant of a decimal matrix.'
      write ( *, '(a)' ) ' '
 
      dec_digit = 7

      k = 0
      do i = 1, n3
        do j = 1, n3
          k = k + 1
          a3(i,j) = k
        end do
      end do

      do j = 1, n3
        do i = 1, n3
          b3(i,j) = 0
        end do
      end do
 
      call decmat_print ( n3, n3, a3, b3, '  The 123/456/789 matrix:' )

      call decmat_det ( n3, a3, b3, dec_digit, dtop, dbot, ierror )
 
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Determinant of the 123/456/789 matrix'
      write ( *, '(2x,i8,a,i8)' ) dtop, ' * 10** ', dbot
 
      do i = 1, n4
        do j = 1, n4
          r = 1.0D+00 / dble ( i + j )
          call r8_to_dec ( r, dec_digit, a4(i,j), b4(i,j) )
        end do
      end do
 
      call decmat_print ( n4, n4, a4, b4, '  The Hilbert matrix:' )

      call decmat_det ( n4, a4, b4, dec_digit, dtop, dbot, ierror )
 
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Determinant of the Hilbert matrix:'
      write ( *, '(2x,i8,a,i8)' ) dtop, ' * 10 ** ', dbot 

      do i = 1, n3
        do j = 1, n3
          if ( i .eq. j ) then
            a3(i,j) = 2
          else if ( i .eq. j+1 .or. i .eq. j-1 ) then
            a3(i,j) = -1
          else
            a3(i,j) = 0
          end if
          b3(i,j) = 0
        end do
      end do
 
      call decmat_print ( n3, n3, a3, b3, '  The -1,2,-1 matrix:' )

      call decmat_det ( n3, a3, b3, dec_digit, dtop, dbot, ierror )
 
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Determinant of the -1,2,-1 matrix:'
      write ( *, '(2x,i8,a,i8)' ) dtop, ' * 10 ** ', dbot
 
      return
      end
      subroutine test021a ( )

c*********************************************************************72
c
cc TEST021a tests DERANGE_ENUM, DERANGE_ENUM2 and DERANGE_ENUM3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer d(0:n)
      integer derange_enum
      integer derange_enum3
      integer i

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST021a'
      write ( *, '(a)' ) '  DERANGE_ENUM counts derangements;'
      write ( *, '(a)' ) '  DERANGE_ENUM2 counts derangements.'
      write ( *, '(a)' ) '  DERANGE_ENUM3 counts derangements.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  N    # of derangements'
      write ( *, '(a)' ) ' '

      call derange_enum2 ( n, d )

      do i = 0, n
        write ( *, '(2x,i8,2x,3i10)' ) 
     &    i, derange_enum(i), d(i), derange_enum3(i)
      end do

      return
      end
      subroutine test022a ( )

c*********************************************************************72
c
cc TEST022a tests DERANGE_BACK_NEXT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      integer a(n)
      integer i
      logical more
      integer number

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST022a'
      write ( *, '(a)' ) '  DERANGE_BACK_NEXT generates derangements'
      write ( *, '(a)' ) '  using backtracking.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 
     &  '  Here, we seek all derangments of order N = ', n
      write ( *, '(a)' ) ' '

      more = .false.
      number = 0

10    continue

        call derange_back_next ( n, a, more )

        if ( .not. more ) then
          go to 20
        end if

        number = number + 1
        write ( *, '(2x,i4,4x,8i4)' ) number, ( a(i), i = 1, n )

      go to 10

20    continue

      return
      end
      subroutine test023a ( )

c*********************************************************************72
c
cc TEST023a tests DERANGE_WEED_NEXT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      integer a(n)
      integer i
      logical more
      integer number

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST023a'
      write ( *, '(a)' ) '  DERANGE_WEED_NEXT generates derangements'
      write ( *, '(a)' ) 
     &  '  by generating ALL permutations, and "weeding out"'
      write ( *, '(a)' ) '  the ones that are not derangements.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 
     &  '  Here, we seek all derangments of order N = ', n
      write ( *, '(a)' ) ' '
    
      more = .false.
      number = 0
 
10    continue

        call derange_weed_next ( n, a, more )

        number = number + 1
        write ( *, '(2x,i4,4x,8i4)' ) number, ( a(i), i = 1, n )

        if ( .not. more ) then
          go to 20
        end if
 
      go to 10

20    continue

      return
      end
      subroutine test024a ( )

c*********************************************************************72
c
cc TEST024a calls DIGRAPH_ARC_EULER.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nedge
      parameter ( nedge = 7 )
      integer nnode
      parameter ( nnode = 5 )

      integer i
      integer in
      integer inode(nedge)
      integer j
      integer jnode(nedge)
      integer jp1
      logical success
      integer trail(nedge)

      data inode / 2, 1, 2, 1, 3, 5, 4 /
      data jnode / 5, 4, 3, 2, 1, 1, 2 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST024a'
      write ( *, '(a)' ) '  DIGRAPH_ARC_EULER finds an Euler '
      write ( *, '(a)' ) '  circuit of a digraph.'

      call digraph_arc_print ( nedge, inode, jnode, 
     &  '  The arc list of the digraph:' )

      call digraph_arc_euler ( nnode, nedge, inode, jnode, 
     &  success, trail )

      if ( success ) then

        call i4vec_print ( nedge, trail, 
     &    '  The edge list of the Euler circuit:' )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  The node list of the Euler circuit:'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '    I  Edge  Node'
        write ( *, '(a)' ) ' '

        do i = 1, nedge

          j = trail(i)

          if ( i .eq. nedge ) then
            jp1 = trail(1)
          else
            jp1 = trail(i+1)
          end if

          if ( jnode(j) .eq. inode(jp1) ) then
            in = jnode(j)
          else
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'The circuit has failedc'
            exit
          end if

          write ( *, '(2x,3i8)' ) i, j, in

        end do

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  The digraph is not eulerian.'
        write ( *, '(a)' ) ' '

      end if

      return
      end
      subroutine test025a ( )

c*********************************************************************72
c
cc TEST025a tests DIOPHANTINE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer test_num
      parameter ( test_num = 20 )

      integer a
      integer a_test(test_num)
      integer b
      integer b_test(test_num)
      integer c
      integer c_test(test_num)
      integer error
      integer ierror
      integer test
      integer x
      integer y

      data a_test /
     &   1027,  1027,  1027, 1027, -1027, 
     &  -1027, -1027, -1027,    6,     0, 
     &      0,     0,     1,    1,     1, 
     &   1024,     0,     0,    5,     2 /
      data b_test /
     &     712,   712, -712, -712, 712, 
     &     712,  -712, -712,    8,   0, 
     &       1,     1,    0,    0,   1, 
     &  -15625,     0,    3,    0,   4 /
      data c_test /
     &       7,    -7,    7,   -7,   7, 
     &      -7,     7,   -7,   50,   0, 
     &       0,     1,    0,    1,   0, 
     &   11529,     1,   11,   19,   7 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST025a'
      write ( *, '(a)' ) 
     &  '  DIOPHANTINE solves a Diophantine equation:'
      write ( *, '(a)' ) '    A * X + B * Y = C'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '        A         B         C         X     Y     Error'
      write ( *, '(a)' ) ' '

      do test = 1, test_num

        a = a_test(test)
        b = b_test(test)
        c = c_test(test)

        call diophantine ( a, b, c, ierror, x, y )

        if ( ierror .ne. 0 ) then
          write ( *, '(2x,3i10,a,i10)' ) 
     &      a, b, c, ' Error code = ', ierror
        else
          error = a * x + b * y - c
          write ( *, '(2x,6i10)' ) a, b, c, x, y, error
        end if

      end do

      return
      end
      subroutine test026a ( )

c*********************************************************************72
c
cc TEST026a tests DIOPHANTINE_SOLUTION_MINIMIZE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer a
      integer b
      integer c
      integer r
      integer x
      integer y

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST026a'
      write ( *, '(a)' ) 
     &  '  DIOPHANTINE_SOLUTION_MINIMIZE computes a minimal'
      write ( *, '(a)' ) 
     &  '  Euclidean norm solution of a Diophantine equation:'
      write ( *, '(a)' ) '    A * X + B * Y = C'

      a = 4096
      b = - 15625
      c = 46116
      x = 665499996
      y = 174456828

      r = a * x + b * y - c

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Coefficients:'
      write ( *, '(a,i12)' ) '    A = ', a
      write ( *, '(a,i12)' ) '    B = ', b
      write ( *, '(a,i12)' ) '    C = ', c
      write ( *, '(a)' ) '  Solution:'
      write ( *, '(a,i12)' ) '    X = ', x
      write ( *, '(a,i12)' ) '    Y = ', y
      write ( *, '(a)' ) '  Residual R = A * X + B * Y - C:'
      write ( *, '(a,i12)' ) '    R = ', r

      call diophantine_solution_minimize ( a, b, x, y )

      r = a * x + b * y - c

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The minimized solution:'
      write ( *, '(a,i12)' ) '    X = ', x
      write ( *, '(a,i12)' ) '    Y = ', y
      write ( *, '(a)' ) '  Residual R = A * X + B * Y - C:'
      write ( *, '(a,i12)' ) '    R = ', r

      x = 15621
      y = 4092

      r = a * x + b * y - c

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The minimal positive solution:'
      write ( *, '(a,i12)' ) '    X = ', x
      write ( *, '(a,i12)' ) '    Y = ', y
      write ( *, '(a)' ) '  Residual R = A * X + B * Y - C:'
      write ( *, '(a,i12)' ) '    R = ', r

      return
      end
      subroutine test026b ( )

c*********************************************************************72
c
cc TEST026B tests DVEC_ADD and DVEC_SUB;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer dvec1(n)
      integer dvec2(n)
      integer dvec3(n)
      integer dvec4(n)
      integer i
      integer i4_uniform
      integer j
      integer k
      integer l
      integer seed
      integer test
      integer test_num

      seed = 123456789
      test_num = 10

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST026B'
      write ( *, '(a)' ) '  DVEC_ADD adds decimal vectors '
      write ( *, '(a)' ) '  representing integers;'
      write ( *, '(a)' ) '  DVEC_SUB subtracts decimal vectors '
      write ( *, '(a)' ) '  representing integers;'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '        I        J        K = I + J    L = I - J'
      write ( *, '(a)' ) ' '

      do test = 1, test_num
        
        i = i4_uniform ( -100, 100, seed )
        j = i4_uniform ( -100, 100, seed )

        write ( *, '(a)' ) ' '

        write ( *, '(2x,i8,2x,i8)' ) i, j

        k = i + j
        l = i - j

        write ( *, '(a20,2x,i8,2x,i8)' ) 
     &    '  Directly:         ', k, l

        call i4_to_dvec ( i, n, dvec1 )
        call i4_to_dvec ( j, n, dvec2 )

        call dvec_add ( n, dvec1, dvec2, dvec3 )
        call dvec_to_i4 ( n, dvec3, k )

        call dvec_sub ( n, dvec1, dvec2, dvec4 )
        call dvec_to_i4 ( n, dvec4, l )

        write ( *, '(a20,2x,i8,2x,i8)' ) '  DVEC_ADD, DVEC_SUB', k, l

      end do

      return
      end
      subroutine test026c ( )

c*********************************************************************72
c
cc TEST026C tests DVEC_COMPLEMENTX;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer dvec1(n)
      integer dvec2(n)
      integer i
      integer j
      integer i4_uniform
      integer seed
      integer test
      integer test_num

      seed = 123456789
      test_num = 5

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST026C'
      write ( *, '(a)' ) '  DVEC_COMPLEMENTX returns the ten''s'
      write ( *, '(a)' ) '  complement of a (signed) decimal vector;'
      write ( *, '(a)' ) ' '

      do test = 1, test_num
        
        i = i4_uniform ( -100, 100, seed )

        call i4_to_dvec ( i, n, dvec1 )

        call dvec_complementx ( n, dvec1, dvec2 )

        call dvec_to_i4 ( n, dvec2, j )

        write ( *, '(a)' ) ' '
        write ( *, '(a,2x,i8)' ) '  I = ', i
        write ( *, '(a,2x,i8)' ) '  J = ', j
        call dvec_print ( n, dvec1, ' ' )
        call dvec_print ( n, dvec2, ' ' )

      end do

      return
      end
      subroutine test026d ( )

c*********************************************************************72
c
cc TEST026D tests DVEC_MUL;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer dvec1(n)
      integer dvec2(n)
      integer dvec3(n)
      integer i
      integer i4_uniform
      integer j
      integer k
      integer n2
      integer seed
      integer test
      integer test2
      integer test_num
      integer test2_num

      seed = 123456789
      test_num = 10
      test2_num = 2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST026D'
      write ( *, '(a)' ) '  DVEC_MUL multiplies decimal vectors '
      write ( *, '(a)' ) '  representing integers;'

      do test2 = 1, test2_num

        if ( test2 .eq. 1 ) then

          n2 = n

        else if ( test2 .eq. 2 ) then

          n2 = 6

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  NOW REPEAT THE TEST...'
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  but use too few digits to represent '
          write ( *, '(a)' ) '  big products.'
          write ( *, '(a)' ) '  This corresponds to an "overflow".'
          write ( *, '(a)' ) '  The result here should get the final'
          write ( *, '(a)' ) '  decimal digits correctly, though.'
        end if

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '        I        J        K = I * J'
        write ( *, '(a)' ) ' '

        do test = 1, test_num
          
          i = i4_uniform ( -1000, 1000, seed )
          j = i4_uniform ( -1000, 1000, seed )

          write ( *, '(a)' ) ' '

          write ( *, '(2x,i8,2x,i8)' ) i, j

          k = i * j

          write ( *, '(a20,2x,i8)' ) '  Directly:         ', k

          call i4_to_dvec ( i, n2, dvec1 )
          call i4_to_dvec ( j, n2, dvec2 )

          call dvec_mul ( n2, dvec1, dvec2, dvec3 )
          call dvec_to_i4 ( n2, dvec3, k )

          write ( *, '(a20,2x,i8)' ) '  DVEC_MUL          ', k

        end do

      end do

      return
      end
      subroutine test027a ( )

c*********************************************************************72
c
cc TEST027a tests EQUIV_NEXT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      integer a(n)
      integer i
      integer jarray(n)
      logical more
      integer nc
      integer rank

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST027a'
      write ( *, '(a)' ) '  EQUIV_NEXT generates all partitions '
      write ( *, '(a)' ) '  of a set.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,6i4)' ) '  Rank/element:', ( i, i = 1, n )
      write ( *, '(a)' ) ' '
     
      rank = 0
      more = .false.

10    continue     
     
        call equiv_next ( n, nc, jarray, a, more )
     
        rank = rank + 1
        write ( *, '(2x,15i4)' ) rank, ( a(i), i = 1, n )
     
        if ( .not. more ) then
          go to 20
        end if

      go to 10

20    continue

      return
      end
      subroutine test028a ( )

c*********************************************************************72
c
cc TEST028a tests EQUIV_NEXT2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      integer a(n)
      logical done
      integer i
      integer rank

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST028a'
      write ( *, '(a)' ) '  EQUIV_NEXT2 generates all partitions '
      write ( *, '(a)' ) '  of a set.'
      write ( *, '(a,i8)' ) '  Here, N = ', n
      rank = 0
      done = .true.
      write ( *, '(a)' ) ' '
      write ( *, '(a,6i4)' ) '  Rank/element:', ( i, i = 1, n )
      write ( *, '(a)' ) ' '
     
10    continue

        call equiv_next2 ( n, a, done )

        if ( done ) then
          go to 20
        end if

        rank = rank + 1
        write ( *, '(2x,i4,10x,6i4)' ) rank, ( a(i), i = 1, n )

      go to 10

20    continue

      return 
      end
      subroutine test029a ( )

c*********************************************************************72
c
cc TEST029a tests EQUIV_RANDOM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      integer a(n)
      double precision b(n)
      integer i
      integer npart
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST029a'
      write ( *, '(a)' ) '  EQUIV_RANDOM selects a random set '
      write ( *, '(a)' ) '  partition.'
     
      seed = 123456789

      do i = 1, 5
     
        call equiv_random ( n, seed, npart, a, b )

        call equiv_print ( n, a, '  The partition:' )
     
      end do
     
      return
      end
      subroutine test0295 ( )

c*********************************************************************72
c
cc TEST0295 tests EULER.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nmax
      parameter ( nmax = 9 )

      integer i
      integer ieuler(0:nmax)
      integer n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0295'
      write ( *, '(a)' ) '  EULER gets rows of Euler''s triangle.'
      write ( *, '(a)' ) ' '

      do n = 0, nmax
        call euler ( n, ieuler )
        write ( *, '(2x,10i7)' ) ( ieuler(i), i = 0, n )
      end do
     
      return
      end
      subroutine test030 ( )

c*********************************************************************72
c
cc TEST030 tests R8_FALL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 December 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      double precision r8_fall
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST030'
      write ( *, '(a)' ) '  R8_FALL computes the falling factorial'
      write ( *, '(a)' ) '  function.'
      write ( *, '(a)' ) 
     &  '  [X]_N = X * (X-1) * (X-2) * ... * ( X-N+1).'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    X          N  R8_FALL(X,N)'
      write ( *, '(a)' ) ' '

      x = 4.0D+00

      do n = -2, 5
     
        write ( *, '(2x,f10.6,i8,f10.6)' ) x, n, r8_fall ( x, n )
     
      end do
     
      return
      end
      subroutine test0304

c*********************************************************************72
c
cc TEST0304 tests FROBENIUS_NUMBER_ORDER2 and FROBENIUS_NUMBER_ORDER2_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer c1
      integer c2
      integer f1
      integer f2
      integer frobenius_number_order2
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0304'
      write ( *, '(a,a)' ) 
     &  '  FROBENIUS_NUMBER_ORDER2 computes Frobenius numbers',
     &  ' of order 2.'
      write ( *, '(a,a)' ) 
     &  '  FROBENIUS_NUMBER_ORDER2_VALUES',
     &  ' returns some exact values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        C1        C1   exact F  comput F'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call frobenius_number_order2_values ( n_data, c1, c2, f1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        f2 = frobenius_number_order2 ( c1, c2 )

        write ( *, '(2x,i8,2x,i8,2x,i8,2x,i8)' ) c1, c2, f1, f2

      go to 10

20    continue

      return
      end
      subroutine test0305

c*********************************************************************72
c
cc TEST0305 tests GAMMA_LOG and GAMMA_LOG_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      double precision fx2
      double precision gamma_log
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0305:'
      write ( *, '(a)' ) '  GAMMA_LOG evaluates the logarithm of the '
      write ( *, '(a)' ) '    Gamma function.'
      write ( *, '(a)' ) '  GAMMA_LOG_VALUES returns exact values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     X       Exact F       GAMMA_LOG(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call gamma_log_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = gamma_log ( x )

        write ( *, '(2x,f8.4,2g18.10)' ) x, fx, fx2

      go to 10

20    continue

      return
      end
      subroutine test031

c*********************************************************************72
c
cc TEST031 tests GRAY_NEXT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      integer a(n)
      integer change
      integer i
      integer k

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST031'
      write ( *, '(a)' ) '  GRAY_NEXT returns the index of the single '
      write ( *, '(a)' ) '  item to be changed in order to get the '
      write ( *, '(a)' ) '  next Gray code.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   K  Change  Gray Code'
      write ( *, '(a)' ) ' '

      change = -n
      k = 0

10    continue

        call gray_next ( n, change, a )

        if ( change .eq. -n ) then
          go to 20
        end if

        write ( *, '(2x,i2,2x,i8,2x,4i1)' ) 
     &    k, change, ( a(i), i = 1, n )

        k = k + 1

      go to 10

20    continue

      return
      end
      subroutine test032

c*********************************************************************72
c
cc TEST032 tests GRAY_RANK and GRAY_UNRANK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer gray
      integer rank
      integer rank2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST032'
      write ( *, '(a)' ) '  GRAY_RANK ranks a Gray code;'
      write ( *, '(a)' ) '  GRAY_UNRANK unranks a Gray code.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    R  =                         RANK'
      write ( *, '(a)' ) '    G  =            GRAY_UNRANK2(RANK)'
      write ( *, '(a)' ) '    R2 = GRAY_RANK2(GRAY_UNRANK2(RANK))'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         R         G         R2'
      write ( *, '(a)' ) ' '
     
      do rank = 0, 24
        call gray_unrank ( rank, gray )
        call gray_rank ( gray, rank2 )
        write ( *, '(2x,i8,2x,i8,2x,i8)' ) rank, gray, rank2
      end do

      return
      end
      subroutine test0321

c*********************************************************************72
c
cc TEST0321 tests GRAY_RANK2 and GRAY_UNRANK2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer gray
      integer rank
      integer rank2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0321'
      write ( *, '(a)' ) '  GRAY_RANK2 ranks a Gray code;'
      write ( *, '(a)' ) '  GRAY_UNRANK2 unranks a Gray code.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    R  =                         RANK'
      write ( *, '(a)' ) '    G  =            GRAY_UNRANK2(RANK)'
      write ( *, '(a)' ) '    R2 = GRAY_RANK2(GRAY_UNRANK2(RANK))'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         R         G         R2'
      write ( *, '(a)' ) ' '
     
      do rank = 0, 24
        call gray_unrank2 ( rank, gray )
        call gray_rank2 ( gray, rank2 )
        write ( *, '(2x,i8,2x,i8,2x,i8)' ) rank, gray, rank2
      end do

      return
      end
      subroutine test0322

c*********************************************************************72
c
cc TEST0322 tests I4_BCLR.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer test_num
      parameter ( test_num = 2 )

      integer i4
      integer i4_test(test_num)
      integer i4_bclr
      integer ivec(0:31)
      integer j1
      integer pos
      integer test

      save i4_test

      data i4_test / 101, -31 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0322'
      write ( *, '(a)' ) '  I4_BCLR sets a given bit to 0.'

      do test = 1, test_num

        i4 = i4_test(test)

        call i4_to_bvec ( i4, 32, ivec )

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Working on I4 = ', i4
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '       Pos     Digit       I4_BCLR'
        write ( *, '(a)' ) ' '

        do pos = 0, 31
  
          j1 = i4_bclr ( i4, pos )

          write ( *, '(2x,i8,2x,i8,2x,i12)' ) pos, ivec(pos), j1

        end do

      end do

      return
      end
      subroutine test03225

c*********************************************************************72
c
cc TEST03225 tests I4_BSET.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer test_num
      parameter ( test_num = 2 )

      integer i4
      integer i4_test(test_num)
      integer i4_bset
      integer ivec(0:31)
      integer j1
      integer pos
      integer test

      save i4_test

      data i4_test / 101, -31 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03225'
      write ( *, '(a)' ) '  I4_BSET sets a given bit to 0.'

      do test = 1, test_num

        i4 = i4_test(test)

        call i4_to_bvec ( i4, 32, ivec )

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Working on I4 = ', i4
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '       Pos     Digit       I4_BSET'
        write ( *, '(a)' ) ' '

        do pos = 0, 31
  
          j1 = i4_bset ( i4, pos )

          write ( *, '(2x,i8,2x,i8,2x,i12)' ) pos, ivec(pos), j1

        end do

      end do

      return
      end
      subroutine test0323

c*********************************************************************72
c
cc TEST0323 tests I4_BTEST.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer test_num
      parameter ( test_num = 2 )

      integer i4
      integer i4_test(test_num)
      logical i4_btest
      integer ivec(0:31)
      logical j1
      integer pos
      integer test

      save i4_test

      data i4_test / 101, -31 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0323'
      write ( *, '(a)' ) '  I4_BTEST reports if a given bit is 0 or 1.'

      do test = 1, test_num

        i4 = i4_test(test)

        call i4_to_bvec ( i4, 32, ivec )
    
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Analyze the integer I4 = ', i4
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '       Pos     Digit  I4_BTEST'
        write ( *, '(a)' ) ' '

        do pos = 0, 31
  
          j1 = i4_btest ( i4, pos )

          write ( *, '(2x,i8,2x,i8,2x,7x,l1,2x)' ) pos, ivec(pos), j1

        end do

      end do

      return
      end
      subroutine test0324

c*********************************************************************72
c
cc TEST0324 tests I4_FACTOR.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer factor_max
      parameter ( factor_max = 10 )

      integer factor(factor_max)
      integer i
      integer n
      integer nfactor
      integer nleft
      integer power(factor_max)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0324'
      write ( *, '(a)' ) '  I4_FACTOR factors an integer,'

      n = 2**2 * 17 * 37

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The integer is ', n

      call i4_factor ( n, factor_max, nfactor, factor, power, nleft )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Prime representation:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         I  FACTOR(I)  POWER(I)'
      write ( *, '(a)' ) ' '
      if ( abs ( nleft ) .ne. 1 ) then
        write ( *, '(2x,i8,i8,a)' ) 0, nleft, ' (UNFACTORED PORTION)'
      end if

      do i = 1, nfactor
        write ( *, '(2x,3i8)' ) i, factor(i), power(i)
      end do

      return
      end
      subroutine test0325

c*********************************************************************72
c
cc TEST0325 tests I4_GCD.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      integer i4_gcd
      integer i4_uniform
      integer j
      integer k
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0325'
      write ( *, '(a)' ) '  I4_GCD computes the greatest common'
      write ( *, '(a)' ) '  divisor of two integers.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         I       J    I4_GCD(I,J)'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do k = 1, 15
        i = i4_uniform ( -5, 15, seed )
        j = i4_uniform ( 1, 15, seed )
        write ( *, '(2x,3i8)' ) i, j, i4_gcd ( i, j )
      end do

      return
      end
      subroutine test0327

c*********************************************************************72
c
cc TEST0327 tests I4_LOG_10.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 21 )

      integer i
      integer i4_log_10
      integer x(n)

      data x /
     &    0,    1,    2,  3,  9, 10, 11,  99, 100, 101, 
     &  999, 1000, 1001, -1, -2, -3, -9, -10, -11, -99, 
     & -101 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0327'
      write ( *, '(a)' ) '  I4_LOG_10: whole part of log base 10,'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         X        I4_LOG_10'
      write ( *, '(a)' ) ' '

      do i = 1, n

        write ( *, '(2x,i8,i12)' ) x(i), i4_log_10 ( x(i) )

      end do

      return
      end
      subroutine test058

c*********************************************************************72
c
cc TEST058 tests I4_PARTITION_CONJ.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 14 )
      integer npart1
      parameter ( npart1 = 4 )

      integer a1(npart1)
      integer a2(n)
      integer mult1(npart1)
      integer mult2(n)
      integer npart2

      data a1 / 2, 5, 1, 4 /
      data mult1 / 1, 1, 3, 1 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST058'
      write ( *, '(a)' ) '  I4_PARTITION_CONJ conjugates an integer'
      write ( *, '(a)' ) '  partition.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Original partition:'
      write ( *, '(a)' ) ' '

      call i4_partition_print ( n, npart1, a1, mult1 )

      call i4_partition_conj ( n, a1, mult1, npart1, a2, mult2, 
     &  npart2 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Conjugate partition:'
      write ( *, '(a)' ) ' '

      call i4_partition_print ( n, npart2, a2, mult2 )

      return
      end
      subroutine test059

c*********************************************************************72
c
cc TEST059 tests I4_PARTITION_NEXT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 7 )

      integer a(n)
      logical done
      integer mult(n)
      integer npart
      integer rank

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST059'
      write ( *, '(a)' ) '  I4_PARTITION_NEXT generates partitions '
      write ( *, '(a)' ) '  of an integer.'
      write ( *, '(a,i8)' ) '  Here N = ', n
      write ( *, '(a)' ) ' '

      rank = 0
      done = .true.
     
10    continue

        call i4_partition_next ( n, npart, a, mult, done )
     
        if ( done ) then
          go to 20 
        end if

        rank = rank + 1

        call i4_partition_print ( n, npart, a, mult )

      go to 10

20    continue
     
      return
      end
      subroutine test060

c*********************************************************************72
c
cc TEST060 tests I4_PARTITION_NEXT2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 7 )

      integer a(n)
      logical more
      integer mult(n)
      integer npart

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST060'
      write ( *, '(a)' ) '  I4_PARTITION_NEXT2 produces partitions'
      write ( *, '(a)' ) '  of an integer.'
      write ( *, '(a)' ) ' '

      more = .false.

10    continue

        call i4_partition_next2 ( n, npart, a, mult, more )

        call i4_partition_print ( n, npart, a, mult )

        if ( .not. more ) then
          go to 20
        end if

      go to 10

20    continue
      
      return
      end
      subroutine test061

c*********************************************************************72
c
cc TEST061 tests I4_PARTITION_COUNT and I4_PARTITION_COUNT_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 20 )

      integer n
      integer n_data
      integer p
      integer p2(0:n_max)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST061'
      write ( *, '(a)' ) '  I4_PARTITION_COUNT counts partitions'
      write ( *, '(a)' ) '  of an integer.'
      write ( *, '(a)' ) '  I4_PARTITION_COUNT_VALUES returns some'
      write ( *, '(a)' ) '  exact values.'

      n_data = 0

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   N     Exact     Count'
      write ( *, '(a)' ) ' '

10    continue

        call i4_partition_count_values ( n_data, n, p )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        call i4_partition_count ( n, p2 )
     
        write ( *, '(2x,i4,2i10)' ) n, p, p2(n)

      go to 10

20    continue

      return
      end
      subroutine test0615

c*********************************************************************72
c
cc TEST0615 tests I4_PARTITION_COUNT2 and I4_PARTITION_COUNT_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 20 )

      integer n
      integer n_data
      integer p
      integer p2(0:n_max)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0615'
      write ( *, '(a)' ) '  I4_PARTITION_COUNT2 counts partitions'
      write ( *, '(a)' ) '  of an integer.'
      write ( *, '(a)' ) '  I4_PARTITION_COUNT_VALUES returns '
      write ( *, '(a)' ) '  some exact values.'

      n_data = 0

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   N     Exact     Count'
      write ( *, '(a)' ) ' '

10    continue

        call i4_partition_count_values ( n_data, n, p )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        call i4_partition_count2 ( n, p2 )
     
        write ( *, '(2x,i4,2i10)' ) n, p, p2(n)

      go to 10

20    continue

      return
      end
      subroutine test062

c*********************************************************************72
c
cc TEST062 tests I4_PARTITION_RANDOM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 8 )

      integer a(n)
      integer i
      integer j
      integer mult(n)
      integer npart
      integer seed
      integer table(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST062'
      write ( *, '(a)' ) '  I4_PARTITION_RANDOM generates a '
      write ( *, '(a)' ) '  random partition.'
      write ( *, '(a)' ) '  I4_PARTITION_COUNT2 sets up a partition '
      write ( *, '(a)' ) '  count table needed by I4_PARTITION_RANDOM.'
      write ( *, '(a)' ) ' '

      seed = 123456789
c
c  Get the partition table.
c
      call i4_partition_count2 ( n, table )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The number of partitions of N'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N    Number of partitions'
      write ( *, '(a)' ) ' '

      do j = 1, n
        write ( *, '(2x,i8,4x,i8)' ) j, table(j)
      end do

      write ( *, '(a)' ) ' '

      do i = 1, 5

        call i4_partition_random ( n, table, seed, a, mult, npart )

        call i4_partition_print ( n, npart, a, mult )

      end do
     
      return
      end
      subroutine test06225 ( )

c*********************************************************************72
c
cc TEST06225 tests I4_PARTITIONS_NEXT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 August 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer s
      parameter ( s = 3 )

      integer i
      integer j
      integer m(s)
      integer msum

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06225'
      write ( *, '(a)' ) '  I4_PARTITIONS_NEXT produces the next'
      write ( *, '(a)' ) '  nondecreasing partitions of an integer, and'
      write ( *, '(a)' ) 
     &  '  if necessary, increments the integer to keep on going.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   I Sum    Partition'
      write ( *, '(a)' ) ' '
      i = 0

      m(1) = 0
      m(2) = 0
      m(3) = 0

      call i4vec_sum ( s, m, msum )

      write ( *, '(2x,i2,2x,i2,4x,10i2)' ) i, msum, ( m(j), j = 1, s )

      do i = 1, 15

        call i4_partitions_next ( s, m )

        call i4vec_sum ( s, m, msum )

        write ( *, '(2x,i2,2x,i2,4x,10i2)' ) i, msum, ( m(j), j = 1, s )

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  You can start from any legal partition.'
      write ( *, '(a)' ) '  Here, we restart at ( 2, 1, 0 ).'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   I Sum    Partition'
      write ( *, '(a)' ) ' '

      i = 0
      m(1) = 2
      m(2) = 1
      m(3) = 0

      call i4vec_sum ( s, m, msum )

      write ( *, '(2x,i2,2x,i2,4x,10i2)' ) i, msum, ( m(j), j = 1, s )

      do i = 1, 15

        call i4_partitions_next ( s, m )

        call i4vec_sum ( s, m, msum )

        write ( *, '(2x,i2,2x,i2,4x,10i2)' ) i, msum, ( m(j), j = 1, s )

      end do

      return
      end
      subroutine test033

c*********************************************************************72
c
cc TEST033 tests I4_SQRT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer q
      integer r

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST033'
      write ( *, '(a)' ) '  I4_SQRT computes the square root'
      write ( *, '(a)' ) '  of an integer.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N  Sqrt(N) Remainder'
      write ( *, '(a)' ) ' '

      do n = -5, 20

        call i4_sqrt ( n, q, r )
        write ( *, '(2x,3i9)' ) n, q, r

      end do

      return
      end
      subroutine test034

c*********************************************************************72
c
cc TEST034 tests I4_SQRT_CF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer max_term
      parameter ( max_term = 100 )

      integer b(0:max_term)
      integer i
      integer n
      integer n_term

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST034'
      write ( *, '(a)' ) '  I4_SQRT_CF computes the continued'
      write ( *, '(a)' ) '  fraction form of the square root of'
      write ( *, '(a)' ) '  an integer.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   N  Period  Whole  Repeating Part'
      write ( *, '(a)' ) ' '
      do n = 1, 20

        call i4_sqrt_cf ( n, max_term, n_term, b )
        write ( *, '(2x,i5,3x,i5,2x,i5,10i5)' ) 
     &    n, n_term, b(0), ( b(i), i = 1, n_term )

      end do

      return
      end
      subroutine test0625

c*********************************************************************72
c
cc TEST0625 tests I4_TO_BVEC and BVEC_TO_I4;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer bvec(n)
      integer i
      integer i2
      integer j

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0625'
      write ( *, '(a)' ) '  I4_TO_BVEC converts an integer to a '
      write ( *, '(a)' ) '    signed binary vector;'
      write ( *, '(a)' ) '  BVEC_TO_I4 converts a signed binary vector'
      write ( *, '(a)' ) '    to an integer;'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  I --> BVEC  -->  I'
      write ( *, '(a)' ) ' '
      do i = -3, 10
        call i4_to_bvec ( i, n, bvec )
        call bvec_to_i4 ( n, bvec, i2 )
        write ( *, '(2x,i3,2x,10i1,2x,i3)' ) 
     &    i, ( bvec(j), j = 1, n ), i2
      end do

      return
      end
      subroutine test035

c*********************************************************************72
c
cc TEST035 tests I4_TO_CHINESE and CHINESE_TO_I4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      integer j
      integer j2
      integer m(n)
      integer r(n)

      save m

      data m / 3, 4, 5, 7 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST035'
      write ( *, '(a)' ) '  I4_TO_CHINESE computes the Chinese'
      write ( *, '(a)' ) '  Remainder representation of an integer.'
      write ( *, '(a)' ) '  CHINESE_TO_I4 computes an integer with '
      write ( *, '(a)' ) '  the given Chinese Remainder '
      write ( *, '(a)' ) '  representation.'

      call i4vec_print ( n, m, '  The moduli:' )

      j = 37

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The number being analyzed is ', j

      call i4_to_chinese ( j, n, m, r )

      call i4vec_print ( n, r, '  The remainders:' )

      call chinese_to_i4 ( n, m, r, j2 )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The reconstructed number is ', j2

      call i4_to_chinese ( j2, n, m, r )

      call i4vec_print ( n, r, 
     &  '  The remainders of the reconstructed number are:' )

      return
      end
      subroutine test0364

c*********************************************************************72
c
cc TEST0364 tests I4_TO_I4POLY and I4_TO_VAN_DER_CORPUT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer degree_max
      parameter ( degree_max = 10 )

      integer base
      integer degree
      double precision h1
      double precision h2
      integer i
      integer i4poly(0:degree_max)
      integer j
      integer n_test
      integer prime
      integer value

      n_test = 10

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0364'
      write ( *, '(a)' ) '  I4_TO_VAN_DER_CORPUT computes elements '
      write ( *, '(a)' ) '    of a van der Corput sequence.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  I4_TO_I4POLY converts an integer to an'
      write ( *, '(a)' ) '    integer polynomial in some base, and'
      write ( *, '(a)' ) '    can be used to mimic the van der Corput'
      write ( *, '(a)' ) '    calculation.'
      write ( *, '(a)' ) ' '

      do j = 1, 3

        base = prime(j)
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  BASE = ', base
        write ( *, '(a)' ) ' '

        do i = 1, n_test

          call i4_to_van_der_corput ( i, base, h1 )

          call i4_to_i4poly ( i, base, degree_max, degree, i4poly )

          call i4vec_reverse ( degree+1, i4poly )

          call i4poly_to_i4 ( degree, i4poly, base, value )

          h2 = dble ( value ) / base**( degree + 1 )

          write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) i, h1,  h2

        end do

      end do

      return
      end
      subroutine test0627

c*********************************************************************72
c
cc TEST0627 tests I4_TO_I4POLY and I4POLY_TO_I4;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer degree_max
      parameter ( degree_max = 5 )
      integer test_num
      parameter ( test_num = 9 )

      integer a(0:degree_max)
      integer base
      integer base_test(test_num)
      integer degree
      integer i
      integer intval
      integer intval2
      integer intval_test(test_num)
      integer test

      data base_test / 2, 2, 2, 3, 4, 5, 6, 23, 24 /
      data intval_test / 1, 6, 23, 23, 23, 23, 23, 23, 23 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0627'
      write ( *, '(a)' ) '  I4_TO_I4POLY converts an integer to a'
      write ( *, '(a)' ) '  polynomial in a given base;'
      write ( *, '(a)' ) '  I4POLY_TO_I4 evaluates an integer '
      write ( *, '(a)' ) '  polynomial.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       I    BASE  DEGREE  Coefficients'
      write ( *, '(a)' ) ' '
      do test = 1, test_num
        intval = intval_test(test)
        base = base_test(test)
        call i4_to_i4poly ( intval, base, degree_max, degree, a )
        write ( *, '(2x,i4,2x,i4,2x,i4,6x,10i8)' ) 
     &    intval, base, degree, ( a(i), i = degree, 0, -1 )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Let I4_TO_I4POLY convert I to a polynomial,'
      write ( *, '(a)' ) '  use I4POLY_TO_I4 to evaluate, and compare.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       I    I2'
      write ( *, '(a)' ) ' '
      do test = 1, test_num
        intval = intval_test(test)
        base = base_test(test)
        call i4_to_i4poly ( intval, base, degree_max, degree, a )
        call i4poly_to_i4 ( degree, a, base, intval2 )
        write ( *, '(2x,i8,2x,i8)' ) intval, intval2
      end do

      return
      end
      subroutine test036

c*********************************************************************72
c
cc TEST036 tests I4_TO_VAN_DER_CORPUT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_prime
      parameter ( n_prime = 5 )

      double precision h(n_prime)
      integer i
      integer j
      integer n_test
      parameter ( n_test = 10 )
      integer p
      integer prime

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST036'
      write ( *, '(a)' ) '  I4_TO_VAN_DER_CORPUT computes the elements '
      write ( *, '(a)' ) '  of a van der Corput sequence.'
      write ( *, '(a)' ) '  The sequence depends on the prime numbers'
      write ( *, '(a)' ) '  used as a base.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Bases:'
      write ( *, '(a)' ) ' '
      write ( *, '(2x,5(8x,i8))' ) ( prime(j), j = 1, n_prime )
      write ( *, '(a)' ) ' '

      do i = 1, n_test
        do j = 1, n_prime
          p = prime(j)
          call i4_to_van_der_corput ( i, p, h(j) )
        end do
        write ( *, '(2x,i3,5g14.6)' ) i, ( h(j), j = 1, n_prime )
      end do

      return
      end
      subroutine test037

c*********************************************************************72
c
cc TEST037 tests I4MAT_01_ROWCOLSUM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 5 )
      integer n
      parameter ( n = 5 )

      integer a(m,n)
      integer c(n)
      integer ierror
      integer r(m)

      save c
      save r

      data c / 2, 2, 2, 2, 1 /
      data r / 3, 2, 2, 1, 1 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST037'
      write ( *, '(a)' ) '  I4MAT_01_ROWCOLSUM constructs a 01 matrix'
      write ( *, '(a)' ) '  with given row and column sums.'
      
      call i4vec_print ( m, r, '  The rowsum vector:' )
      call i4vec_print ( n, c, '  The columnsum vector: ' )

      call i4mat_01_rowcolsum ( m, n, r, c, a, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) 
     &    '  I4MAT_01_ROWCOLSUM returned IERROR = ', ierror
      else
        call i4mat_print ( m, n, a, '  The rowcolsum matrix:' )
      end if

      return
      end
      subroutine test038

c*********************************************************************72
c
cc TEST038 tests I4MAT_01_ROWCOLSUM2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 5 )
      integer n
      parameter ( n = 5 )

      integer a(m,n)
      integer c(n)
      integer ierror
      integer r(m)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST038'
      write ( *, '(a)' ) '  I4MAT_01_ROWCOLSUM2 constructs a 01 matrix'
      write ( *, '(a)' ) '  with given row and column sums.'

      c(1) = 2
      c(2) = 1
      c(3) = 2
      c(4) = 2
      c(5) = 2

      r(1) = 2
      r(2) = 1
      r(3) = 3
      r(4) = 1
      r(5) = 2

      call i4vec_print ( m, r, '  The rowsum vector:' )
      call i4vec_print ( n, c, '  The columnsum vector: ' )

      call i4mat_01_rowcolsum2 ( m, n, r, c, a, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) 
     &    '  I4MAT_01_ROWCOLSUM2 returned error flag IERROR = ', ierror
        write ( *, '(a)' ) 
     &  '  The matrix returned is not an exact solution.'
      end if

      call i4mat_print ( m, n, a, '  The rowcolsum matrix:' )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Now repeat, with data for which there is not'
      write ( *, '(a)' ) 
     &  '  a solution.  The program will try its best anyway.'

      c(1) = 1
      c(2) = 4
      c(3) = 1
      c(4) = 5
      c(5) = 1

      r(1) = 1
      r(2) = 3
      r(3) = 4
      r(4) = 1
      r(5) = 3

      call i4vec_print ( m, r, '  The rowsum vector:' )
      call i4vec_print ( n, c, '  The columnsum vector: ' )

      call i4mat_01_rowcolsum2 ( m, n, r, c, a, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) 
     &    '  I4MAT_01_ROWCOLSUM2 returned IERROR = ', ierror
        write ( *, '(a)' ) '  The matrix is not an exact solution.'
      end if

      call i4mat_print ( m, n, a, '  The rowcolsum matrix:' )

      return
      end
      subroutine test039

c*********************************************************************72
c
cc TEST039 tests I4MAT_01_ROWCOLSUM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 13 )
      integer n
      parameter ( n = 17 )

      integer a(m,n)
      integer c(n)
      integer c2(n)
      integer c3(n)
      integer cperm(n)
      integer i
      integer ierror
      integer j
      integer r(m)
      integer r2(m)
      integer r3(m)
      integer rperm(m)
      integer total

      save c
      save r

      data c /  4,  4, 11, 10, 10,  
     &          8,  9, 10,  8,  9,  
     &          3, 10,  4,  7,  9,  
     &          3,  3 /
      data r / 14, 13, 14, 10, 12,  
     &          2, 10,  1, 10, 11,  
     &          6,  2, 17 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST039'
      write ( *, '(a)' ) '  I4MAT_01_ROWCOLSUM constructs a 01 matrix'
      write ( *, '(a)' ) '  with given row and column sums.'

      call i4vec_print ( m, r, '  The rowsum vector R:' )
      call i4vec_sort_heap_index_d ( m, r, rperm )

      do i = 1, m
        r2(i) = r(rperm(i))
      end do

      call i4vec_print ( n, c, '  The columnsum vector C: ' )
      call i4vec_sort_heap_index_d ( n, c, cperm )

      do i = 1, n
        c2(i) = c(cperm(i))
      end do

      call i4mat_01_rowcolsum ( m, n, r2, c2, a, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &    '  I4MAT_01_ROWCOLSUM returned error flag IERROR = ', 
     &    ierror
        return
      end if
c
c  Invert the R and C permutations.
c
      call i4mat_perm2 ( m, n, a, rperm, cperm )

      call i4mat_print ( m, n, a, '  The rowcolsum matrix:' )

      do i = 1, m
        total = 0
        do j = 1, n
          total = total + a(i,j)
        end do
        r3(i) = total
      end do

      call i4vec_print ( m, r3, '  Computed row sums' )

      do j = 1, n
        total = 0
        do i = 1, m
          total = total + a(i,j)
        end do
        c3(j) = total
      end do

      call i4vec_print ( n, c3, '  Computed column sums:' )

      return
      end
      subroutine test040

c*********************************************************************72
c
cc TEST040 tests I4MAT_U1_INVERSE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 6 )

      integer a(n,n)
      integer b(n,n)

      save a

      data a /
     &  1, 0, 0, 0, 0, 0, 
     &  1, 1, 0, 0, 0, 0, 
     &  0, 0, 1, 0, 0, 0, 
     &  0, 0, 1, 1, 0, 0, 
     &  0, 0, 0, 0, 1, 0, 
     & 75, 0, 0, 0, 1, 1 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST040'
      write ( *, '(a)' ) '  I4MAT_U1_INVERSE inverts a '
      write ( *, '(a)' ) '  unit upper triangular matrix.'

      call i4mat_print ( n, n, a, '  The input matrix:' )
     
      call i4mat_u1_inverse ( n, a, b )
     
      call i4mat_print ( n, n, b, '  The inverse:' )
     
      return
      end
      subroutine test041

c*********************************************************************72
c
cc TEST041 test I4MAT_PERM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 9 )

      integer a(n,n)
      integer i
      integer j
      integer p(n)

      save p

      data p / 2,3,9,6,7,8,5,4,1 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST041'
      write ( *, '(a)' ) '  I4MAT_PERM reorders an integer matrix '
      write ( *, '(a)' ) '  in place.  The rows and columns use the '
      write ( *, '(a)' ) '  same permutation.'

      do i = 1, n
        do j = 1, n
          a(i,j) = i * 10 + j
        end do
      end do

      call i4mat_print ( n, n, a, '  The input matrix:' )
     
      call perm_print ( n, p, '  The row and column permutation:' )
     
      call i4mat_perm ( n, a, p )
     
      call i4mat_print ( n, n, a, '  The permuted matrix:' )
     
      return
      end
      subroutine test042

c*********************************************************************72
c
cc TEST042 test I4MAT_PERM2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 9 )
      integer n
      parameter ( n = 7 )

      integer a(m,n)
      integer i
      integer p(m)
      integer q(n)
      integer j

      save p
      save q

      data p / 2, 3, 9, 6, 7, 8, 5, 4, 1 /
      data q / 3, 4, 5, 6, 7, 1, 2 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST042'
      write ( *, '(a)' ) '  I4MAT_PERM2 reorders a matrix in place.'
      write ( *, '(a)' ) 
     &  '  Rows and columns use different permutations.'

      do i = 1, m
        do j = 1, n
          a(i,j) = i * 10 + j
        end do
      end do
     
      call i4mat_print ( m, n, a, '  The input matrix:' )
     
      call perm_print ( m, p, '  The row permutation:' )

      call perm_print ( n, q, '  The column permutation:' )
     
      call i4mat_perm2 ( m, n, a, p, q )
     
      call i4mat_print ( m, n, a, '  The permuted matrix:' )

      return
      end
      subroutine test047

c*********************************************************************72
c
cc TEST047 tests I4MAT_PERM2 and TRIANG.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer a(n,n)
      integer p(n)

      save a

      data a /
     &  1,0,1,0,1,0,1,0,0,1, 
     &  0,1,0,0,1,0,0,0,0,0, 
     &  0,0,1,0,1,0,1,0,0,1, 
     &  0,1,1,1,1,1,1,1,0,1, 
     &  0,0,0,0,1,0,0,0,0,0, 
     &  0,1,0,0,1,1,1,0,0,0, 
     &  0,0,0,0,1,0,1,0,0,0, 
     &  0,1,0,0,1,1,1,1,0,1, 
     &  0,0,0,0,0,0,0,0,0,0, 
     &  0,0,0,0,1,0,1,0,0,1 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST047'
      write ( *, '(a)' ) '  TRIANG relabels elements for a '
      write ( *, '(a)' ) '  partial ordering,'
      write ( *, '(a)' ) '  I4MAT_PERM2 reorders an integer matrix'
      write ( *, '(a)' ) '  in place.'

      call i4mat_print ( n, n, a, '  The input matrix:' )
     
      call triang ( n, a, p )
     
      call perm_print ( n, p, '  The new ordering:' )

      call i4mat_perm2 ( n, n, a, p, p )
     
      call i4mat_print ( n, n, a, '  The reordered matrix:' )
     
      return
      end
      subroutine test043

c*********************************************************************72
c
cc TEST043 tests INDEX_BOX_NEXT_2D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      integer j
      logical more
      integer n1
      parameter ( n1 = 5 )
      integer n2
      parameter ( n2 = 3 )
      integer n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST043'
      write ( *, '(a)' ) '  INDEX_BOX_NEXT_2D produces IJ indices that'
      write ( *, '(a)' ) '  lie on the surface of a box in 2D.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The box has logical dimensions:'
      write ( *, '(3x,3i3)' ) n1, n2
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   #    I   J'
      write ( *, '(a)' ) ' '

      more = .false.
      n = 0

10    continue

        call index_box_next_2d ( n1, n2, i, j, more )

        if ( .not. more ) then
          go to 20
        end if

        n = n + 1
        write ( *, '(2x,4i3)' ) n, i, j

      go to 10

20    continue

      return
      end
      subroutine test044

c*********************************************************************72
c
cc TEST044 tests INDEX_BOX_NEXT_3D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      integer j
      integer k
      logical more
      integer n1
      parameter ( n1 = 5 )
      integer n2
      parameter ( n2 = 3 )
      integer n3
      parameter ( n3 = 4 )
      integer n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST044'
      write ( *, '(a)' ) '  INDEX_BOX_NEXT_3D produces IJK indices that'
      write ( *, '(a)' ) '  lie on the surface of a box.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The box has logical dimensions:'
      write ( *, '(3x,3i3)' ) n1, n2, n3
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     #    I   J   K'
      write ( *, '(a)' ) ' '

      more = .false.
      n = 0

10    continue

        call index_box_next_3d ( n1, n2, n3, i, j, k, more )

        if ( .not. more ) then
          go to 20
        end if

        n = n + 1
        write ( *, '(2x,4i3)' ) n, i, j, k

      go to 10

20    continue

      return
      end
      subroutine test045

c*********************************************************************72
c
cc TEST045 tests INDEX_BOX2_NEXT_2D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      integer ic
      parameter ( ic = 10 )
      integer j
      integer jc
      parameter ( jc = 20 )
      logical more
      integer n1
      parameter ( n1 = 4 )
      integer n2
      parameter ( n2 = 3 )
      integer n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST045'
      write ( *, '(a)' ) '  INDEX_BOX2_NEXT_2D produces IJ indices that'
      write ( *, '(a)' ) '  lie on the surface of a box2 in 2D.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The box has half-widths:'
      write ( *, '(3x,3i3)' ) n1, n2
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  and has center cell:'
      write ( *, '(3x,2i3)' ) ic, jc
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     #    I   J'
      write ( *, '(a)' ) ' '

      more = .false.
      n = 0

10    continue

        call index_box2_next_2d ( n1, n2, ic, jc, i, j, more )

        if ( .not. more ) then
          go to 20
        end if

        n = n + 1
        write ( *, '(2x,4i3)' ) n, i, j

      go to 10

20    continue

      return
      end
      subroutine test046

c*********************************************************************72
c
cc TEST046 tests INDEX_BOX2_NEXT_3D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      integer ic
      parameter ( ic = 10 )
      integer j
      integer jc
      parameter ( jc = 20 )
      integer k
      integer kc
      parameter ( kc = 30 )
      logical more
      integer n1
      parameter ( n1 = 5 )
      integer n2
      parameter ( n2 = 3 )
      integer n3
      parameter ( n3 = 4 )
      integer n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST046'
      write ( *, '(a)' ) '  INDEX_BOX2_NEXT_3D produces IJK indices'
      write ( *, '(a)' ) '  that lie on the surface of a box.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The box has half widths:'
      write ( *, '(3x,3i3)' ) n1, n2, n3
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  and central cell:'
      write ( *, '(3x,3i3)' ) ic, jc, kc
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  We will only print a PORTION of the datac'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   #    I   J   K'
      write ( *, '(a)' ) ' '

      more = .false.
      n = 0

10    continue

        call index_box2_next_3d ( n1, n2, n3, ic, jc, kc, i, j, 
     &    k, more )

        if ( .not. more ) then
          go to 20
        end if

        n = n + 1

        if ( n .le. 10 .or. 370 .le. n ) then
          write ( *, '(2x,4i3)' ) n, i, j, k
        end if

      go to 10

20    continue

      return
      end
      subroutine test048

c*********************************************************************72
c
cc TEST048 tests INDEX_NEXT0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      integer a(n)
      integer hi
      parameter ( hi = 3 )
      integer i
      logical more

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST048'
      write ( *, '(a)' ) '  INDEX_NEXT0 generates all indices of an'
      write ( *, '(a)' ) '  array of given shape, with'
      write ( *, '(a)' ) '  lower limit 1 and given upper limit.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of index entries = ', n
      write ( *, '(a,i8)' ) '  Coordinate maximum HI =   ', hi
     
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Index arrays:'
      write ( *, '(a)' ) ' '

      more = .false.

10    continue

        call index_next0 ( n, hi, a, more )

        write ( *, '(2x,3i4)' ) ( a(i), i = 1, n )

        if ( .not. more ) then
          go to 20
        end if

      go to 10

20    continue

      return
      end
      subroutine test049

c*********************************************************************72
c
cc TEST049 tests INDEX_NEXT1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      integer a(n)
      integer hi(n)
      integer i
      logical more

      save hi

      data hi / 4, 2, 3 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST049'
      write ( *, '(a)' ) '  INDEX_NEXT1 generates all indices of an'
      write ( *, '(a)' ) '  array of given shape, with'
      write ( *, '(a)' ) '  lower limit 1 and given upper limits.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of index entries = ', n

      call i4vec_print ( n, hi, '  Coordinate maximum indices:' )
     
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Index arrays:'
      write ( *, '(a)' ) ' '

      more = .false.

10    continue

        call index_next1 ( n, hi, a, more )

        write ( *, '(2x,3i4)' ) ( a(i), i = 1, n )

        if ( .not. more ) then
          go to 20
        end if

      go to 10

20    continue

      return
      end
      subroutine test050

c*********************************************************************72
c
cc TEST050 tests INDEX_NEXT2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      integer a(n)
      integer hi(n)
      integer i
      integer lo(n)
      logical more

      save hi
      save lo

      data hi / 11, -3, 1 /
      data lo / 10, -5, 0 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST050'
      write ( *, '(a)' ) '  INDEX_NEXT2 generates all indices of an'
      write ( *, '(a)' ) '  array of given shape with given'
      write ( *, '(a)' ) '  lower and upper limits.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of index entries = ', n
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Coordinate, Maximum Index'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(3i10)' ) i, lo(i), hi(i)
      end do
     
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Index arrays:'
      write ( *, '(a)' ) ' '

      more = .false.

10    continue

        call index_next2 ( n, lo, hi, a, more )

        write ( *, '(2x,3i4)' ) ( a(i), i = 1, n )

        if ( .not. more ) then
          go to 20
        end if

      go to 10

20    continue

      return
      end
      subroutine test051

c*********************************************************************72
c
cc TEST051 tests INDEX_RANK0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      integer a(n)
      integer hi
      parameter ( hi = 3 )
      integer rank

      save a

      data a / 3, 1, 2 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST051'
      write ( *, '(a)' ) '  INDEX_RANK0 ranks an index with'
      write ( *, '(a)' ) '  lower limit 1 and given upper limit.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of index entries = ', n
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Coordinate maximum Index = ', hi
      write ( *, '(a)' ) ' '

      call i4vec_print ( n, a, '  The index array:' )

      call index_rank0 ( n, hi, a, rank )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The rank of this object is ', rank

      return
      end
      subroutine test052

c*********************************************************************72
c
cc TEST052 tests INDEX_RANK1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      integer a(n)
      integer hi(n)
      integer i
      integer rank

      save a
      save hi

      data a / 4, 1, 2 /
      data hi / 4, 2, 3 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST052'
      write ( *, '(a)' ) '  INDEX_RANK1 ranks an index with'
      write ( *, '(a)' ) '  lower limit 1 and given upper limits.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of index entries = ', n
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Coordinate, Maximum Index'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,2i10)' ) i, hi(i)
      end do
     
      call i4vec_print ( n, a, '  The index array:' )

      call index_rank1 ( n, hi, a, rank )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The rank of this object is ', rank

      return
      end
      subroutine test053

c*********************************************************************72
c
cc TEST053 tests INDEX_RANK2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      integer a(n)
      integer hi(n)
      integer i
      integer lo(n)
      integer rank

      save a
      save hi
      save lo

      data a  / 1, 11, 5 /
      data hi / 2, 11, 6 /
      data lo / 1, 10, 4 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST053'
      write ( *, '(a)' ) '  INDEX_RANK2 ranks an index with given'
      write ( *, '(a)' ) '  lower and upper limits.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of index entries = ', n
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Coordinate, Minimum index, Maximum Index'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,3i10)' ) i, lo, hi(i)
      end do
     
      call i4vec_print ( n, a, '  The index array:' )

      call index_rank2 ( n, lo, hi, a, rank )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The rank of this object is ', rank

      return
      end
      subroutine test054

c*********************************************************************72
c
cc TEST054 tests INDEX_UNRANK0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      integer a(n)
      integer hi
      parameter ( hi = 3 )
      integer i
      integer maxrank
      integer rank

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST054'
      write ( *, '(a)' ) '  INDEX_UNRANK0 unranks a multi-index.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The multi-index has dimension ', n
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The upper limit is HI = ', hi
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Rank, Multi-Index:'
      write ( *, '(a)' ) ' '
     
      maxrank = hi**n

      do rank = 1, maxrank
        call index_unrank0 ( n, hi, rank, a )
        write ( *, '(2x,i3,3i8)' ) rank, ( a(i), i = 1, n )
      end do
     
      return
      end
      subroutine test055

c*********************************************************************72
c
cc TEST055 tests INDEX_UNRANK1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      integer a(n)
      integer hi(n)
      integer i
      integer i4vec_product
      integer maxrank
      integer rank

      save hi

      data hi / 4, 2, 3 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST055'
      write ( *, '(a)' ) '  INDEX_UNRANK1 unranks a multi-index.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The multi-index has dimension ', n
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The upper limits are:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,2i10)' ) i, hi(i)
      end do
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Rank, Multi-Index:'
      write ( *, '(a)' ) ' '
     
      maxrank = i4vec_product ( n, hi )

      do rank = 1, maxrank
        call index_unrank1 ( n, hi, rank, a )
        write ( *, '(2x,i3,3i8)' ) rank, ( a(i), i = 1, n )
      end do
     
      return
      end
      subroutine test056

c*********************************************************************72
c
cc TEST056 tests INDEX_UNRANK2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      integer a(n)
      integer hi(n)
      integer i
      integer lo(n)
      integer rank

      save hi
      save lo

      data hi / 2, 11, 6 /
      data lo / 1, 10, 4 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST056'
      write ( *, '(a)' ) '  INDEX_UNRANK2 unranks a multi-index.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The multi-index has dimension ', n
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The lower and upper limits are:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,3i10)' ) i, lo(i), hi(i)
      end do
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Rank, Multi-Index:'
      write ( *, '(a)' ) ' '
     
      rank = 7

      call index_unrank2 ( n, lo, hi, rank, a )
      write ( *, '(2x,i3,3i8)' ) rank, ( a(i), i = 1, n )
     
      return
      end
      subroutine test057

c*********************************************************************72
c
cc TEST057 tests INS_PERM and PERM_INS.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      integer i
      integer ins(n)
      integer perm(n)
      integer perm2(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST057'
      write ( *, '(a)' ) '  PERM_INS computes the inversion sequence.'
      write ( *, '(a)' ) '  INS_PERM recovers the permutation.'
      write ( *, '(a)' ) ' '

      perm(1) = 3
      perm(2) = 5
      perm(3) = 1
      perm(4) = 4
      perm(5) = 2

      call perm_ins ( n, perm, ins )

      call ins_perm ( n, ins, perm2 )

      write ( *, '(2x,6i3)' ) ( i, i = 1, n )
      write ( *, '(2x,6i3)' ) ( perm(i), i = 1, n )
      write ( *, '(2x,6i3)' ) ( ins(i), i = 1, n )
      write ( *, '(2x,6i3)' ) ( perm2(i), i = 1, n )
     
      return
      end
      subroutine test063

c*********************************************************************72
c
cc TEST063 tests INVOLUTE_ENUM;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer i
      integer s(0:n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST063'
      write ( *, '(a)' ) '  INVOLUTE_ENUM counts involutions;'
      write ( *, '(a)' ) ' '

      call involute_enum ( n, s )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  N    # of involutions'
      write ( *, '(a)' ) ' '

      do i = 0, n
        write ( *, '(2x,2i10)' ) i, s(i)
      end do

      return
      end
      subroutine test064

c*********************************************************************72
c
cc TEST064 test I4POLY.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 6 )
      integer test_num
      parameter ( test_num = 6 )

      integer a(n)
      integer i
      integer iopt
      integer test
      integer val
      integer x0

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST064'
      write ( *, '(a)' ) '  I4POLY converts between power sum'
      write ( *, '(a)' ) '  factorial, and Taylor forms, and can '
      write ( *, '(a)' ) '  evaluate a polynomial'
      write ( *, '(a)' ) ' '
     
      do test = 1, test_num

        if ( test .eq. 1 ) then
          iopt = -3
        else if ( test .eq. 2 ) then
          iopt = -2
        else if ( test .eq. 3 ) then
          iopt = -1
          x0 = 2
        else if ( test .eq. 4 ) then
          iopt = 0
          x0 = 2
        else if ( test .eq. 5 ) then
          iopt = 6
          x0 = 2
        else if ( test .eq. 6 ) then
          iopt = 6
          x0 = -2
        end if

        a(1) = 0
        a(2) = 0
        a(3) = 0
        a(4) = 0
        a(5) = 0
        a(6) = 1

        if ( test .eq. 1 ) then
          write ( *, '(a)' ) '  All calls have input A as follows'
          write ( *, '(10i4)' ) ( a(i), i = 1, n )
        end if
     
        call i4poly ( n, a, x0, iopt, val )
     
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Option IOPT = ', iopt
        if ( -1 .le. iopt ) then
          write ( *, '(a,i8)' ) '  X0 = ', x0
        end if

        if ( iopt .eq. -3 .or. iopt .eq. -2 .or. 0 .lt. iopt ) then
          write ( *, '(a)' ) '  Output array = '
          write ( *, '(2x,10i4)' ) ( a(i), i = 1, n )
        end if

        if ( iopt .eq. -1 .or. iopt .eq. 0 ) then
          write ( *, '(a,i8)' ) '  Value = ', val
        end if
     
      end do

      return
      end
      subroutine test065

c*********************************************************************72
c
cc TEST065 tests I4POLY_CYCLO.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer max_n
      parameter ( max_n = 10 )

      integer phi(0:max_n)
      integer n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST065'
      write ( *, '(a)' ) 
     &  '  I4POLY_CYCLO computes cyclotomic polynomials.'

      do n = 0, max_n

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  N = ', n
        write ( *, '(a)' ) ' '
        call i4poly_cyclo ( n, phi )

        call i4poly_print ( n, phi, '  The cyclotomic polynomial:' )

      end do

      return
      end
      subroutine test066

c*********************************************************************72
c
cc TEST066 tests I4POLY_DIV.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer a(0:10)
      integer b(0:10)
      integer na
      integer nb
      integer nq
      integer nr
      integer test_num
      integer q(0:10)
      integer r(0:10)
      integer test

      test_num = 2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST066'
      write ( *, '(a)' ) '  I4POLY_DIV computes the quotient and'
      write ( *, '(a)' ) '  remainder for polynomial division.'
      write ( *, '(a)' ) ' '
c
c  1: Divide X**3 + 2*X**2 - 5*X - 6  by X-2.  
c     Quotient is 3+4*X+X**2, remainder is 0.
c
c  2: Divide X**4 + 3*X**3 + 2*X**2 - 2  by  X**2 + X - 3.
c     Quotient is X**2 + 2*X + 3, remainder 8*X + 7.
c
      do test = 1, test_num

        if ( test .eq. 1 ) then
          na = 3
          a(0) = -6
          a(1) = -5
          a(2) =  2
          a(3) =  1
          nb = 1
          b(0) = -2
          b(1) =  1
        else if ( test .eq. 2 ) then
          na = 4
          a(0) = -2
          a(1) =  5
          a(2) =  2
          a(3) =  3
          a(4) =  1
          nb = 2
          b(0) = -3
          b(1) =  1
          b(2) =  1
        end if

        call i4poly_print ( na, a, 
     &    '  The polynomial to be divided, A:' )
        call i4poly_print ( nb, b, '  The divisor polynomial, B:' )

        call i4poly_div ( na, a, nb, b, nq, q, nr, r )
     
        call i4poly_print ( nq, q, '  The quotient polynomial, Q:' )
        call i4poly_print ( nr, r, '  The remainder polynomial, R:' )

      end do

      return
      end
      subroutine test067

c*********************************************************************72
c
cc TEST067 tests I4POLY_MUL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer maxn
      parameter ( maxn = 5 )

      integer a(0:maxn)
      integer b(0:maxn)
      integer c(0:maxn)
      integer test
      integer na
      integer nb
      integer test_num

      test_num = 2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST067'
      write ( *, '(a)' ) '  I4POLY_MUL multiplies two polynomials.'
      write ( *, '(a)' ) ' '
c
c  1: Multiply (1+X) times (1-X).  Answer is 1-X**2.
c  2: Multiply (1+2*X+3*X**2) by (1-2*X). Answer is 1 + 0*X - X**2 - 6*X**3
c
      do test = 1, test_num

        if ( test .eq. 1 ) then
          na = 1
          a(0) = 1
          a(1) = 1
          nb = 1
          b(0) =  1
          b(1) = -1
        else if ( test .eq. 2 ) then
          na = 2
          a(0) = 1
          a(1) = 2
          a(2) = 3
          nb = 1
          b(0) =  1
          b(1) = -2
        end if

        call i4poly_mul ( na, a, nb, b, c )

        call i4poly_print ( na, a, '  The factor A:' )

        call i4poly_print ( nb, b, '  The factor B:' )

        call i4poly_print ( na+nb, c, '  The product C = A*B:' )

      end do

      return
      end
      subroutine test0675

c*********************************************************************72
c
cc TEST0675 tests I4VEC_DESCENDS;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c

      implicit none

      integer n
      parameter ( n = 4 )

      integer a(n)
      logical descends
      integer i
      logical i4vec_descends
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0675'
      write ( *, '(a)' ) '  I4VEC_DESCENDS is true if an integer '
      write ( *, '(a)' ) '  vector decreases.'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do i = 1, 5

        call i4vec_uniform ( n, 1, n, seed, a )

        call i4vec_print ( n, a, '  The integer array to search:' )
     
        descends = i4vec_descends ( n, a )

        if ( descends ) then
          write ( *, '(a)' ) '  The vector is descending.'
        else
          write ( *, '(a)' ) '  The vector is not descending.'
        end if

      end do

      return
      end
      subroutine test068

c*********************************************************************72
c
cc TEST068 tests I4VEC_FRAC;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c

      implicit none

      integer n
      parameter ( n = 10 )

      integer a(n)
      integer afrac
      integer k
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST068'
      write ( *, '(a)' ) '  I4VEC_FRAC: K-th smallest entry.'
      write ( *, '(a)' ) ' '

      seed = 123456789

      call i4vec_uniform ( n, 1, 2*n, seed, a )

      call i4vec_print ( n, a, '  The integer array to search:' )
     
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       K       K-th smallest'
      write ( *, '(a)' ) ' '

      do k = 1, n

        call i4vec_frac ( n, a, k, afrac )
        write ( *, '(2x,i8,2x,i8)' ) k, afrac

      end do

      return
      end
      subroutine test0683

c*********************************************************************72
c
cc TEST0683 tests I4VEC_INDEX;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n 
      parameter ( n = 20 )

      integer a(n)
      integer aval
      integer first
      integer i4vec_index
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0683'
      write ( *, '(a)' ) '  I4VEC_INDEX returns the index of the '
      write ( *, '(a)' ) '  first occurrence of a given value in an '
      write ( *, '(a)' ) '  integer vector.'
      write ( *, '(a)' ) ' '

      seed = 123456789

      call i4vec_uniform ( n, 1, n/2, seed, a )

      aval = a(n/2)

      call i4vec_print ( n, a, '  The integer array to search:' )
     
      first = i4vec_index ( n, a, aval )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The value searched for is ', aval
      write ( *, '(a,i8)' ) 
     &  '  The index of first occurrence is ', first

      return
      end
      subroutine test0685

c*********************************************************************72
c
cc TEST0685 tests I4VEC_MAXLOC_LAST;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 20 )

      integer a(n)
      integer i4vec_maxloc_last
      integer last
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0685'
      write ( *, '(a)' ) '  I4VEC_MAXLOC_LAST: index of the last'
      write ( *, '(a)' ) '  maximal entry in an integer vector.'
      write ( *, '(a)' ) ' '

      seed = 123456789

      call i4vec_uniform ( n, 1, n/4, seed, a )

      call i4vec_print ( n, a, '  The integer array to search:' )
     
      last = i4vec_maxloc_last ( n, a )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Index of last maximal entry is ', last

      return
      end
      subroutine test0686

c*********************************************************************72
c
cc TEST0686 tests I4VEC_PAIRWISE_PRIME;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      integer a(n)
      integer i
      logical i4vec_pairwise_prime
      logical pairwise_prime
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0686'
      write ( *, '(a)' ) '  I4VEC_PAIRWISE_PRIME is true if an integer'
      write ( *, '(a)' ) '  vector is pairwise prime.'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do i = 1, 5

        call i4vec_uniform ( n, 1, n, seed, a )

        call i4vec_print ( n, a, '  The integer array to check:' )
     
        pairwise_prime = i4vec_pairwise_prime ( n, a )

        if ( pairwise_prime ) then
          write ( *, '(a)' ) '  The vector is pairwise prime.'
        else
          write ( *, '(a)' ) '  The vector is not pairwise prime.'
        end if

      end do

      return
      end
      subroutine test0687

c*********************************************************************72
c
cc TEST0687 tests I4VEC_REVERSE;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      integer a(n)
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0687'
      write ( *, '(a)' ) '  I4VEC_REVERSE reverses an integer vector.'
      write ( *, '(a)' ) ' '

      seed = 123456789

      call i4vec_uniform ( n, 1, n, seed, a )

      call i4vec_print ( n, a, '  The integer array:' )
     
      call i4vec_reverse ( n, a )

      call i4vec_print ( n, a, '  The reversed integer array:' )

      return
      end
      subroutine test0688

c*********************************************************************72
c
cc TEST0688 tests I4VEC_SORT_BUBBLE_A.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 20 )

      integer a(n)
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0688'
      write ( *, '(a)' ) '  I4VEC_SORT_BUBBLE_A ascending sorts an'
      write ( *, '(a)' ) '  integer vector using bubble sort.'
      write ( *, '(a)' ) ' '

      seed = 123456789

      call i4vec_uniform ( n, 0, 3*n, seed, a )

      call i4vec_print ( n, a, '  Unsorted array:' )

      call i4vec_sort_bubble_a ( n, a )

      call i4vec_print ( n, a, '  Sorted array:' )

      return
      end
      subroutine test0689

c*********************************************************************72
c
cc TEST0689 tests I4VEC_SORT_HEAP_INDEX_D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 20 )

      integer a(n)
      integer i
      integer indx(n)
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0689'
      write ( *, '(a)' ) '  I4VEC_SORT_HEAP_INDEX_D descending'
      write ( *, '(a)' ) '  index-sorts an integer vector using'
      write ( *, '(a)' ) '  heap sort.'
      write ( *, '(a)' ) ' '

      seed = 123456789

      call i4vec_uniform ( n, 0, 3*n, seed, a )

      call i4vec_print ( n, a, '  Unsorted array:' )

      call i4vec_sort_heap_index_d ( n, a, indx )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       I    INDX    A(INDX)'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,3i8)' ) i, indx(i), a(indx(i))
      end do

      return
      end
      subroutine test06895

c*********************************************************************72
c
cc TEST06895 tests INVERSE_MOD_N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 November 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer b
      integer n
      integer y
      integer z

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06895'
      write ( *, '(a)' ) 
     &  '  INVERSE_MOD_N seeks Y, the inverse of B mod N,'
      write ( *, '(a)' ) '  so that mod ( B * Y, N ) = 1, but returns 0'
      write ( *, '(a)' ) '  if the inverse does not exist.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     B     N     Y     Z = mod ( B * Y, N )'
      write ( *, '(a)' ) ' '

      do n = 1, 10
        do b = 1, n - 1
          call inverse_mod_n ( b, n, y )
          z = mod ( b * y, n )
          write ( *, '(2x,i4,2x,i4,2x,i4,2x,i4)' ) b, n, y, z
        end do
      end do

      return
      end
      subroutine test069 ( )

c*********************************************************************72
c
cc TEST069 tests JFRAC_TO_RFRAC and RFRAC_TO_JFRAC.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 October 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer maxm
      parameter ( maxm = 10 )

      integer i
      integer m
      double precision p(maxm)
      double precision q(maxm)
      double precision r(maxm)
      double precision s(maxm)
      integer seed
c
c  Generate the data, but force Q(M+1) to be 1.  
c  That will make it easier to see that the two operations are inverses
c  of each other.  JFRAC_TO_RFRAC is free to scale its output, and chooses
c  a scaling in which Q(M+1) is 1.
c
      seed = 123456789

      m = 6
      call r8vec_uniform_01 ( m, seed, p )
      call r8vec_uniform_01 ( m + 1, seed, q )

      do i = 1, m
        q(i) = q(i) / q(m+1)
      end do
      q(m+1) = 1.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST069'
      write ( *, '(a)' ) '  RFRAC_TO_JFRAC converts a rational'
      write ( *, '(a)' ) '  polynomial fraction to a J fraction.'
      write ( *, '(a)' ) '  JFRAC_TO_RFRAC converts a J fraction'
      write ( *, '(a)' ) '  to a rational polynomial fraction.'
      write ( *, '(a)' ) ' '

      write ( *, '(a)' ) '  The original rational coefficients:'
      write ( *, '(a)' ) ' '
      write ( *, '(2x,5g14.6)' ) ( p(i), i = 1, m )
      write ( *, '(2x,5g14.6)' ) ( q(i), i = 1, m+1 )
     
      call rfrac_to_jfrac ( m, p, q, r, s )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The J fraction coefficients:'
      write ( *, '(a)' ) ' '
      write ( *, '(2x,5g14.6)' ) ( r(i), i = 1, m )
      write ( *, '(2x,5g14.6)' ) ( s(i), i = 1, m )
     
      call jfrac_to_rfrac ( m, r, s, p, q )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The recovered rational polynomial:'
      write ( *, '(a)' ) ' '
      write ( *, '(2x,5g14.6)' ) ( p(i), i = 1, m )
      write ( *, '(2x,5g14.6)' ) ( q(i), i = 1, m+1 )

      return
      end
      subroutine test070

c*********************************************************************72
c
cc TEST070 tests JOSEPHUS.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer k
      integer m
      integer n
      integer x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST070'
      write ( *, '(a)' ) '  JOSEPHUS solves Josephus problems.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    N    M    K	 X'
      write ( *, '(a)' ) ' '

      m = 3
      n = 41
      k = 41
      call josephus ( n, m, k, x )
      write ( *, '(2x,4i5)' ) n, m, k, x

      m = -38
      n = 41
      k = 41
      call josephus ( n, m, k, x )

      write ( *, '(2x,4i5)' ) n, m, k, x

      m = 3
      n = 41
      k = 40
      call josephus ( n, m, k, x )

      write ( *, '(2x,4i5)' ) n, m, k, x

      m = 2
      n = 64
      k = 64
      call josephus ( n, m, k, x )

      write ( *, '(2x,4i5)' ) n, m, k, x

      m = 2
      n = 1000
      k = 1000
      call josephus ( n, m, k, x )

      write ( *, '(2x,4i5)' ) n, m, k, x

      return
      end
      subroutine test071

c*********************************************************************72
c
cc TEST071 tests KSUB_NEXT.
c

c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c

      implicit none

      integer k
      parameter ( k = 3 )

      integer a(k)
      integer i
      logical more
      integer n
      integer rank

      n = 5
      do i = 1, k
        a(i) = 0
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST071'
      write ( *, '(a)' ) '  KSUB_NEXT generates all K subsets of an'
      write ( *, '(a)' ) '  N set in lexicographic order.'
      write ( *, '(a)' ) ' '

      more = .false.
      rank = 0
     
10    continue

        call ksub_next ( n, k, a, more )

        rank = rank + 1
        write ( *, '(2x,i4,4x,8i4)' ) rank, ( a(i), i = 1, k )

        if ( .not. more ) then
          go to 20
        end if

      go to 10

20    continue
     
      return
      end
      subroutine test072

c*********************************************************************72
c
cc TEST072 tests KSUB_NEXT2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer k
      parameter ( k = 3 )

      integer a(k)
      integer i
      integer in
      integer iout
      logical more
      integer n
      integer rank

      n = 5

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST072'
      write ( *, '(a)' ) '  KSUB_NEXT2 generates the next K subset of'
      write ( *, '(a)' ) '  an N set by the revolving door method.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Rank  Subset  Added  Removed'
      write ( *, '(a)' ) ' '
c
c  KSUB_NEXT2 doesn't have a good way of stopping.  
c  We will save the starting subset, and stop when the
c  new subset is the same as the starting one.
c
      in = 0
      iout = 0
      rank = 0
     
      call i4vec_indicator ( k, a )
     
10    continue
     
        rank = rank + 1
        write ( *, '(2x,i4,2x,3i2,3x,i2,7x,i2)' ) 
     &    rank, ( a(i), i = 1, k ), in, iout
     
        call ksub_next2 ( n, k, a, in, iout )
     
        more = .false.

        do i = 1, k
          if ( a(i) .ne. i ) then
            more = .true.
          end if
        end do

        if ( .not. more ) then
          go to 20
        end if

      go to 10

20    continue
     
      return
      end
      subroutine test073

c*********************************************************************72
c
cc TEST073 tests KSUB_NEXT3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer k
      parameter ( k = 3 )

      integer a(k)
      integer i
      integer in
      integer iout
      logical more
      integer n
      integer rank

      n = 5

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST073'
      write ( *, '(a)' ) '  KSUB_NEXT3 generates all K subsets of an'
      write ( *, '(a)' ) '  N set using the revolving door method.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Rank    Subset  Added Removed'
      write ( *, '(a)' ) ' '

      rank = 0
      more = .false.
     
10    continue

        call ksub_next3 ( n, k, a, more, in, iout )

        rank = rank + 1
        write ( *, '(2x,i4,2x,3i2,3x,i2,7x,i2)' ) 
     &    rank, ( a(i), i = 1, k ), in, iout

        if ( .not. more ) then
          go to 20
        end if

      go to 10

20    continue

      return
      end
      subroutine test074

c*********************************************************************72
c
cc TEST074 tests KSUB_NEXT4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer k
      parameter ( k = 3 )

      integer a(k)
      logical done
      integer i
      integer n
      integer rank

      n = 5

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST074'
      write ( *, '(a)' ) '  KSUB_NEXT4 generates K subsets of an N set.'
      write ( *, '(a,i8)' ) '  N = ', n
      write ( *, '(a,i8)' ) '  K=  ', k
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Rank    Subset'
      write ( *, '(a)' ) ' '

      done = .true.
      rank = 0
     
10    continue
     
        call ksub_next4 ( n, k, a, done )
     
        if ( done ) then
          go to 20
        end if

        rank = rank + 1
        write ( *, '(2x,i4,4x,3i4)' ) rank, ( a(i), i = 1, k )

      go to 10

20    continue
     
      return
      end
      subroutine test075

c*********************************************************************72
c
cc TEST075 tests KSUB_RANDOM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 July 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer k
      parameter ( k = 3 )

      integer a(k)
      integer i
      integer j
      integer n
      integer seed

      n = 5

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST075'
      write ( *, '(a)' ) '  KSUB_RANDOM generates a random K subset'
      write ( *, '(a)' ) '  of an N set.'
      write ( *, '(a,i8)' ) '  Set size is N =    ', n
      write ( *, '(a,i8)' ) '  Subset size is K = ', k
      write ( *, '(a)' ) ' '

      seed = 123456789

      do i = 1, 5

        call ksub_random ( n, k, seed, a )

        write ( *, '(2x,8i3)' ) ( a(j), j = 1, k )

      end do
     
      return
      end
      subroutine test076

c*********************************************************************72
c
cc TEST076 tests KSUB_RANDOM2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 July 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer k
      parameter ( k = 3 )

      integer a(k)
      integer i
      integer j
      integer n
      integer seed

      n = 5

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST076'
      write ( *, '(a)' ) '  KSUB_RANDOM2 generates a random K subset '
      write ( *, '(a)' ) '  of an N set.'
      write ( *, '(a,i8)' ) '  Set size is N =    ', n
      write ( *, '(a,i8)' ) '  Subset size is K = ', k
      write ( *, '(a)' ) ' '
     
      seed = 123456789

      do i = 1, 5
        call ksub_random2 ( n, k, seed, a )
        write ( *, '(2x,8i3)' ) ( a(j), j = 1, k )
      end do
     
      return
      end
      subroutine test077

c*********************************************************************72
c
cc TEST077 tests KSUB_RANDOM3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      integer a(n)
      integer i
      integer j
      integer k
      integer seed

      k = 3

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST077'
      write ( *, '(a)' ) '  KSUB_RANDOM3 generates a random K-subset '
      write ( *, '(a)' ) '  of an N-set.'
      write ( *, '(a,i8)' ) '  Set size is N =    ', n
      write ( *, '(a,i8)' ) '  Subset size is K = ', k
      write ( *, '(a)' ) ' '

      seed = 123456789

      do i = 1, 10
        call ksub_random3 ( n, k, seed, a )
        write ( *, '(2x,15i3)' ) ( a(j), j = 1, n )
      end do
     
      return
      end
      subroutine test0771

c*********************************************************************72
c
cc TEST0771 tests KSUB_RANDOM4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 July 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer k 
      parameter ( k = 3 )

      integer a(k)
      integer i
      integer j
      integer n
      integer seed

      n = 5

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0771'
      write ( *, '(a)' ) '  KSUB_RANDOM4 generates a random K subset'
      write ( *, '(a)' ) '  of an N set.'
      write ( *, '(a,i8)' ) '  Set size is N =    ', n
      write ( *, '(a,i8)' ) '  Subset size is K = ', k
      write ( *, '(a)' ) ' '

      seed = 123456789

      do i = 1, 5

        call ksub_random4 ( n, k, seed, a )

        write ( *, '(2x,8i3)' ) ( a(j), j = 1, k )

      end do
     
      return
      end
      subroutine test07715

c*********************************************************************72
c
cc TEST07715 tests KSUB_RANDOM5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 June 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer k 
      parameter ( k = 5 )

      integer a(k)
      integer i
      integer j
      integer n
      integer seed

      n = 52

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07715'
      write ( *, '(a)' ) '  KSUB_RANDOM5 generates a random K subset'
      write ( *, '(a)' ) '  of an N set.'
      write ( *, '(a,i8)' ) '  Set size is N =    ', n
      write ( *, '(a,i8)' ) '  Subset size is K = ', k
      write ( *, '(a)' ) ' '

      seed = 123456789

      do i = 1, 5

        call ksub_random5 ( n, k, seed, a )

        write ( *, '(2x,8i3)' ) ( a(j), j = 1, k )

      end do
     
      return
      end
      subroutine test0772

c*********************************************************************72
c
cc TEST0772 tests KSUB_RANK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )
      integer k
      parameter ( k = 3 )

      integer a(k)
      integer i
      integer rank

      save a

      data a / 1, 3, 5 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0772'
      write ( *, '(a)' ) '  KSUB_RANK: rank of a K subset of an N set.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  For N = ', n
      write ( *, '(a,i8)' ) '  and K = ', k
      write ( *, '(a)' ) '  the subset is:'
      write ( *, '(5i4)' ) ( a(i), i = 1, k )
     
      call ksub_rank ( k, a, rank )
     
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The rank is ', rank
     
      return
      end
      subroutine test0773

c*********************************************************************72
c
cc TEST0773 tests KSUB_UNRANK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer k
      parameter ( k = 3 )

      integer a(k)
      integer i
      integer n
      integer rank

      n = 5
      rank = 8

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0773'
      write ( *, '(a)' ) '  KSUB_UNRANK: find the K-subset of an N set'
      write ( *, '(a)' ) '  of a given rank.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  N is ', n
      write ( *, '(a,i8)' ) '  K is ', k
      write ( *, '(a,i8)' ) '  and the desired rank is ', rank
     
      call ksub_unrank ( k, rank, a )
     
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The subset of the given rank is:'
      write ( *, '(2x,5i4)' ) ( a(i), i = 1, k )
     
      return
      end
      subroutine test078

c*********************************************************************72
c
cc TEST078 tests MATRIX_PRODUCT_OPT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 6 )

      integer cost
      integer i
      integer order(n-1)
      integer rank(n+1)

      save rank

      data rank / 4, 2, 3, 1, 2, 2, 3 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST078'
      write ( *, '(a)' ) '  MATRIX_PRODUCT_OPT seeks the optimal'
      write ( *, '(a)' ) '  order for a chain of matrix products.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Matrix ranks:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   I    R    C'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,7i5)' ) i, rank(i), rank(i+1)
      end do

      call matrix_product_opt ( n, rank, cost, order )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Optimal cost is ', cost

      call i4vec_print ( n-1, order, '  Ordering:' )

      return
      end
      subroutine test079

c*********************************************************************72
c
cc TEST079 tests MOEBIUS_MATRIX.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 11 )

      integer ih(n,n)
      integer matrix(n,n)

      save ih

      data ih /
     &  0,0,1,1,0,0,0,0,0,0,0, 
     &  0,0,0,0,0,0,0,1,0,0,0, 
     &  0,1,0,0,0,0,0,0,0,0,0, 
     &  0,1,0,0,0,0,0,0,0,0,0, 
     &  0,0,0,1,0,0,0,0,0,0,0, 
     &  1,0,0,0,1,0,0,0,1,0,0, 
     &  0,0,0,0,0,1,0,0,0,1,1, 
     &  0,0,0,0,0,0,0,0,0,0,0, 
     &  0,0,1,1,0,0,0,0,0,0,0, 
     &  1,0,0,0,0,0,0,0,1,0,0, 
     &  0,0,0,0,0,0,0,0,1,0,0 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST079'
      write ( *, '(a)' ) '  MOEBIUS_MATRIX computes the Moebius matrix.'
     
      call i4mat_print ( n, n, ih, '  The input matrix:' )

      call moebius_matrix ( n, ih, matrix )
     
      call i4mat_print ( n, n, matrix, '  The Moebius matrix:' )
     
      return
      end
      subroutine test080

c*********************************************************************72
c
cc TEST080 tests MORSE_THUE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 100 )

      integer i
      integer ihi
      integer ilo
      integer s(0:n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST080'
      write ( *, '(a)' ) '  MORSE_THUE computes the Morse-Thue numbers.'
      write ( *, '(a)' ) ' '

      do i = 0, n
        call morse_thue ( i, s(i) )
      end do

      do ilo = 0, n, 10
        ihi = min ( ilo + 9, n )
        write ( *, '(4x,40i1)' ) ( s(i), i = ilo, ihi )
      end do

      return
      end
      subroutine test081

c*********************************************************************72
c
cc TEST081 tests MULTINOMIAL_COEF1 and MULTINOMIAL_COEF2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer maxfactor
      parameter ( maxfactor = 5 )

      integer factor(maxfactor)
      integer i
      integer j
      integer n
      integer ncomb1
      integer ncomb2
      integer nfactor

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST081'
      write ( *, '(a)' ) '  MULTINOMIAL_COEF1 computes multinomial'
      write ( *, '(a)' ) '    coefficients using the Gamma function;'
      write ( *, '(a)' ) '  MULTINOMIAL_COEF2 computes multinomial'
      write ( *, '(a)' ) '    coefficients directly.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Line 10 of the BINOMIAL table:'
      write ( *, '(a)' ) ' '

      n = 10
      nfactor = 2

      do i = 0, n

        factor(1) = i
        factor(2) = n - i

        call multinomial_coef1 ( nfactor, factor, ncomb1 )

        call multinomial_coef2 ( nfactor, factor, ncomb2 )

        write ( *, '(2x,i4,i4,3x,i5,i5)' ) 
     &    factor(1), factor(2), ncomb1, ncomb2

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Level 5 of the TRINOMIAL coefficients:'

      n = 5
      nfactor = 3

      do i = 0, n

        factor(1) = i

        write ( *, '(a)' ) ' '

        do j = 0, n - factor(1)

          factor(2) = j
          factor(3) = n - factor(1) - factor(2)

          call multinomial_coef1 ( nfactor, factor, ncomb1 )

          call multinomial_coef2 ( nfactor, factor, ncomb2 )

          write ( *, '(2x,i4,i4,i4,3x,i5,i5)' ) 
     &      factor(1), factor(2), factor(3), ncomb1, ncomb2

        end do

      end do

      return
      end
      subroutine test082

c*********************************************************************72
c
cc TEST082 tests NETWORK_FLOW_MAX.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nnode
      parameter ( nnode = 6 )
      integer nedge
      parameter ( nedge = 20 )

      integer i
      integer icut(nnode)
      integer icpflo(2,nedge)
      integer iendpt(2,nedge)
      integer node_flow(nnode)
      integer sink
      integer source

      save icpflo
      save iendpt

      data icpflo /
     &  3,0,7,0,2,0,5,0,4,0,1,0,4,0,2,0,8,0,3,0, 
     &  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /
      data iendpt /
     &  1,2, 1,3, 2,3, 2,4, 2,5, 3,4, 3,5, 4,5, 4,6, 5,6, 
     &  2,1, 3,1, 3,2, 4,2, 5,2, 4,3, 5,3, 5,4, 6,4, 6,5 /

      sink = 6
      source = 1

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST082'
      write ( *, '(a)' ) '  NETWORK_FLOW_MAX finds the maximum '
      write ( *, '(a)' ) '  flow on a network.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The source is node ', source
      write ( *, '(a,i8)' ) '  The sink is node   ', sink
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Endpoint array:'
      write ( *, '(a)' ) ' '
      write ( *, '(2x,20i3)' ) ( iendpt(1,i), i = 1, nedge )
      write ( *, '(2x,20i3)' ) ( iendpt(2,i), i = 1, nedge )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Input edge capacity array:'
      write ( *, '(a)' ) ' '
      write ( *, '(2x,20i3)' ) ( icpflo(1,i), i = 1, nedge )
     
      call network_flow_max ( nnode, nedge, iendpt, icpflo, source, 
     &  sink, icut, node_flow )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Reordered endpoint array:'
      write ( *, '(a)' ) ' '
      write ( *, '(2x,20i3)' ) ( iendpt(1,i), i = 1, nedge )
      write ( *, '(2x,20i3)' ) ( iendpt(2,i), i = 1, nedge )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Output edge capacity/flow array:'
      write ( *, '(a)' ) ' '
      write ( *, '(2x,20i3)' ) ( icpflo(1,i), i = 1, nedge )
      write ( *, '(2x,20i3)' ) ( icpflo(2,i), i = 1, nedge )

      call i4vec_print ( nnode, icut, '  Minimal node cut vector:' )

      call i4vec_print ( nnode, node_flow, '  Nodal flow vector:' )

      return
      end
      subroutine test083

c*********************************************************************72
c
cc TEST083 tests NIM_SUM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 32 )

      integer i
      integer i4_uniform
      integer i1
      integer i1vec(n)
      integer i2
      integer i2vec(n)
      integer i3
      integer i3vec(n)
      integer ihi
      integer ilo
      integer test_num
      integer seed

      ihi = 1000
      ilo = 0
      test_num = 5

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST083'
      write ( *, '(a)' ) 
     &  '  NIM_SUM computes the Nim sum of two integers.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    I    J    Nim(I+J)'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do i = 1, test_num

        i1 = i4_uniform ( ilo, ihi, seed )
        call ui4_to_ubvec ( i1, n, i1vec )

        i2 = i4_uniform ( ilo, ihi, seed )
        call ui4_to_ubvec ( i2, n, i2vec )

        call nim_sum ( i1, i2, i3 )
        call ui4_to_ubvec ( i3, n, i3vec )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  I1, I2, I3 in decimal:'
        write ( *, '(a)' ) ' '
        write ( *, '(i5)' ) i1
        write ( *, '(i5)' ) i2
        write ( *, '(i5)' ) i3
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  I1, I2, I3 in binary:'
        write ( *, '(a)' ) ' '
        call bvec_print ( n, i1vec, ' ' )
        call bvec_print ( n, i2vec, ' ' )
        call bvec_print ( n, i3vec, ' ' )

      end do

      return
      end
      subroutine test0835

c*********************************************************************72
c
cc TEST0835 tests PADOVAN.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 15 )

      integer i
      integer p(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0835'
      write ( *, '(a)' ) '  PADOVAN computes the Padovan numbers.'
      write ( *, '(a)' ) ' '

      call padovan ( n, p )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N        P(N)'
      write ( *, '(a)' ) ' '

      do i = 0, n-1
        write ( *, '(2x,i8,2x,i10)' ) i, p(i+1)
      end do

      return
      end
      subroutine test084

c*********************************************************************72
c
cc TEST084 tests PELL_BASIC and PELL_NEXT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer d
      integer q
      integer r
      integer x0
      integer x1
      integer y0
      integer y1

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST084'
      write ( *, '(a)' ) '  PELL_BASIC solves the basic Pell equation.'
      write ( *, '(a)' ) '  PELL_NEXT computes the "next" solution.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       D       X        Y         R'
      write ( *, '(a)' ) ' '

      do d = 2, 20

        call i4_sqrt ( d, q, r )

        if ( r /= 0 ) then

          call pell_basic ( d, x0, y0 )

          r = x0**2 - d * y0**2

          write ( *, '(2x,4i9)' ) d, x0, y0, r

          call pell_next ( d, x0, y0, x0, y0, x1, y1 )

          r = x1**2 - d * y1**2

          write ( *, '(2x,9x,3i9)' ) x1, y1, r

        end if

      end do

      return
      end
      subroutine test085

c*********************************************************************72
c
cc TEST085 tests PENT_ENUM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer i
      integer pi

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST085'
      write ( *, '(a)' ) '  PENT_ENUM counts points in pentagons.'
      write ( *, '(a)' ) ' '

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  N    Pent(N)'
      write ( *, '(a)' ) ' '

      do i = 0, n
        call pent_enum ( i, pi )
        write ( *, '(2x,2i10)' ) i, pi
      end do

      return
      end
      subroutine test093

c*********************************************************************72
c
cc TEST093 tests PERM_ASCEND.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 9 )

      integer length
      integer p(n)
      integer subseq(n)

      save p

      data p / 2,3,9,6,7,8,5,4,1 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST093'
      write ( *, '(a)' ) '  PERM_ASCEND determines the length of '
      write ( *, '(a)' ) '  the longest increasing subsequence in '
      write ( *, '(a)' ) '  a permutation.'

      call perm_print ( n, p, '  The permutation:' )

      call perm_ascend ( n, p, length, subseq )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 
     &  '  The length of the longest increasing subsequence is ', 
     &  length

      call i4vec_print ( length, subseq, 
     &  '  A longest increasing subsequence:' )

      return
      end
      subroutine test086

c*********************************************************************72
c
cc TEST086 tests PERM_BREAK_COUNT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 6 )

      integer break_count
      integer p(n)

      save p

      data p / 4, 5, 2, 1, 6, 3 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST086'
      write ( *, '(a)' ) '  PERM_BREAK_COUNT counts the breaks'
      write ( *, '(a)' ) '  in a permutation.'
     
      call perm_print ( n, p, '  The permutation:' )
     
      call perm_break_count ( n, p, break_count )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The number of breaks is ', break_count

      return
      end
      subroutine test087

c*********************************************************************72
c
cc TEST087 tests PERM_CANON_TO_CYCLE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 6 )

      integer p1(n)
      integer p2(n)

      save p1

      data p1 / 4, 5, 2, 1, 6, 3 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST087'
      write ( *, '(a)' ) '  PERM_CANON_TO_CYCLE converts a permutation'
      write ( *, '(a)' ) '  from canonical to cycle form.'
     
      call perm_print ( n, p1, '  The permutation in canonical form:' )
     
      call perm_canon_to_cycle ( n, p1, p2 )

      call perm_print ( n, p2, '  The permutation in cycle form:' )
     
      return
      end
      subroutine test088

c*********************************************************************72
c
cc TEST088 tests PERM_CYCLE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 9 )

      integer iopt
      integer isgn
      integer ncycle
      integer p(n)

      save p

      data p / 2, 3, 9, 6, 7, 8, 5, 4, 1 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST088'
      write ( *, '(a)' ) '  PERM_CYCLE analyzes a permutation.'
     
      call perm_print ( n, p, '  The permutation:' )
     
      iopt = 1
      call perm_cycle ( n, iopt, p, isgn, ncycle )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  NCYCLE = ', ncycle
      write ( *, '(a,i8)' ) '  ISGN =   ', isgn

      call perm_print ( n, p, '  The permutation in cycle form:' )
     
      return
      end
      subroutine test089

c*********************************************************************72
c
cc TEST089 tests PERM_CYCLE_TO_CANON.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 6 )

      integer p1(n)
      integer p2(n)

      save p1

      data p1 / -6, 3, 1, -5, 4, -2 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST089'
      write ( *, '(a)' ) '  PERM_CYCLE_TO_CANON converts a permutation'
      write ( *, '(a)' ) '  from cycle to canonical form.'
     
      call perm_print ( n, p1, '  The permutation in cycle form:' )
     
      call perm_cycle_to_canon ( n, p1, p2 )

      call perm_print ( n, p2, '  The permutation in canonical form:' )
     
      return
      end
      subroutine test090

c*********************************************************************72
c
cc TEST090 tests PERM_CYCLE_TO_INDEX and PERM_INDEX_TO_CYCLE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 9 )

      integer p1(n)
      integer p2(n)
      integer p3(n)

      save p1

      data p1 / 2, 3, 9, 6, 7, 8, 5, 4, 1 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST090'
      write ( *, '(a)' ) '  PERM_CYCLE_TO_INDEX converts a permutation'
      write ( *, '(a)' ) '    from cycle to standard index form.'
      write ( *, '(a)' ) '  PERM_INDEX_TO_CYCLE converts a permutation'
      write ( *, '(a)' ) '    from standard index to cycle form.'
     
      call perm_print ( n, p1, 
     &  '  The standard index form permutation:' )
     
      call perm_index_to_cycle ( n, p1, p2 )

      call perm_print ( n, p2, '  The permutation in cycle form:' )

      call perm_cycle_to_index ( n, p2, p3 )
     
      call perm_print ( n, p3, 
     &  '  The standard index form permutation:' )
     
      return
      end
      subroutine test091

c*********************************************************************72
c
cc TEST091 tests PERM_DISTANCE
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer k11
      integer k12
      integer k13
      integer k21
      integer k23
      integer p1(n)
      integer p2(n)
      integer p3(n)
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST091'
      write ( *, '(a)' ) '  PERM_DISTANCE computes the Ulam metric'
      write ( *, '(a)' ) '  distance between two permutations.'

      seed = 123456789

      call perm_random3 ( n, seed, p1 )
      call perm_print ( n, p1, '  Permutation P1' )
      call perm_random3 ( n, seed, p2 )
      call perm_print ( n, p2, '  Permutation P2' )
      call perm_random3 ( n, seed, p3 )
      call perm_print ( n, p3, '  Permutation P3' )

      call perm_distance ( n, p1, p1, k11 )
      call perm_distance ( n, p1, p2, k12 )
      call perm_distance ( n, p2, p1, k21 )
      call perm_distance ( n, p1, p3, k13 )
      call perm_distance ( n, p2, p3, k23 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  K(P1,P1) should be 0.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  K(P1,P1) = ', k11
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  K(P1,P2) should equal K(P2,P1).'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  K(P1,P2) = ', k12
      write ( *, '(a,i8)' ) '  K(P2,P1) = ', k21
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  K(P1,P3) <= K(P1,P2) + K(P2,P3).'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  K(P1,P3) = ', k13
      write ( *, '(a,i8)' ) '  K(P1,P2) = ', k12
      write ( *, '(a,i8)' ) '  K(P2,P3) = ', k23
      write ( *, '(a,i8)' ) '  K(P1,P2) + K(P2,P3) = ', k12 + k23

      return
      end
      subroutine test092

c*********************************************************************72
c
cc TEST092 tests PERM_FIXED_ENUM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer fnm
      integer m

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST092'
      write ( *, '(a)' ) '  PERM_FIXED_ENUM enumerates the permutations'
      write ( *, '(a)' ) '  of N objects that leave M unchanged.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  For this test, N = ', n
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  M    F(N,M)'
      write ( *, '(a)' ) ' '

      do m = 0, n

        call perm_fixed_enum ( n, m, fnm )
        write ( *, '(2x,i3,2x,i8)' ) m, fnm

      end do

      return
      end
      subroutine test094

c*********************************************************************72
c
cc TEST094 tests PERM_INVERSE;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 7 )

      integer p(n)

      save p

      data p / 4, 3, 5, 1, 7, 6, 2 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST094'
      write ( *, '(a)' ) 
     &  '  PERM_INVERSE inverts a permutation in place;'
      write ( *, '(a)' ) ' '

      call perm_print ( n, p, '  The original permutation:' )
     
      call perm_inverse ( n, p )
     
      call perm_print ( n, p, '  The inverted permutation:' )
     
      return
      end
      subroutine test095

c*********************************************************************72
c
cc TEST095 tests PERM_INVERSE2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 7 )

      integer p(n)

      save p

      data p / 4, 3, 5, 1, 7, 6, 2 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST095'
      write ( *, '(a)' ) 
     &  '  PERM_INVERSE2 inverts a permutation in place.'

      call perm_print ( n, p, '  The original permutation:' )
     
      call perm_inverse2 ( n, p )
     
      call perm_print ( n, p, '  The inverted permutation:' )
     
      return
      end
      subroutine test0955

c*********************************************************************72
c
cc TEST0955 tests PERM_INVERSE3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 7 )

      integer p(n)
      integer p_inv(n)

      save p

      data p / 4, 3, 5, 1, 7, 6, 2 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0955'
      write ( *, '(a)' ) '  PERM_INVERSE3 inverts a permutation.'

      call perm_print ( n, p, '  The original permutation:' )
     
      call perm_inverse3 ( n, p, p_inv )
     
      call perm_print ( n, p_inv, '  The inverted permutation:' )
     
      return
      end
      subroutine test096

c*********************************************************************72
c
cc TEST096 tests PERM_LEX_NEXT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      logical more
      integer p(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST096'
      write ( *, '(a)' ) 
     &  '  PERM_LEX_NEXT generates permutations in order.'
      write ( *, '(a)' ) ' '
      more = .false.
     
10    continue

        call perm_lex_next ( n, p, more )

        if ( .not. more ) then
          go to 20
        end if

        call perm_print ( n, p, ' ' )

      go to 10

20    continue
     
      return
      end
      subroutine test097

c*********************************************************************72
c
cc TEST097 tests PERM_LEX_NEXT and PERM_SIGN.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      integer i
      logical more
      integer p(n)
      integer rank
      integer p_sign

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST097'
      write ( *, '(a)' ) 
     &  '  PERM_LEX_NEXT generates permutations in order.'
      write ( *, '(a)' ) 
     &  '  PERM_SIGN computes the sign of a permutation.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  RANK  SIGN  Permutation'
      write ( *, '(a)' ) ' '

      more = .false.
      rank = 0 

10    continue

        call perm_lex_next ( n, p, more )
        call perm_sign ( n, p, p_sign )

        if ( .not. more ) then
          go to 20
        end if

        write ( *, '(2x,i4,2x,i4,2x,10i4)' ) 
     &    rank, p_sign, ( p(i), i = 1, n )

        rank = rank + 1

      go to 10

20    continue
     
      return
      end
      subroutine test098

c*********************************************************************72
c
cc TEST098 tests PERM_MUL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      integer p1(n)
      integer p2(n)
      integer p3(n)
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST098'
      write ( *, '(a)' ) '  PERM_MUL multiplies two permutations.'
      write ( *, '(a)' ) ' '

      seed = 123456789

      call perm_random ( n, seed, p1 )
      call perm_random ( n, seed, p2 )

      call perm_print ( n, p1, '  Permutation P1:' )

      call perm_print ( n, p2, '  Permutation P2:' )

      call perm_mul ( n, p1, p2, p3 )

      call perm_print ( n, p3, '  Product permutation:' )

      return
      end
      subroutine test099

c*********************************************************************72
c
cc TEST099 tests PERM_NEXT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      logical even
      logical more
      integer p(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST099'
      write ( *, '(a)' ) '  PERM_NEXT generates permutations.'
      write ( *, '(a)' ) ' '
      more = .false.
     
10    continue

        call perm_next ( n, p, more, even )

        call perm_print ( n, p, ' ' )

        if ( .not. more ) then
          go to 20
        end if
     
      go to 10

20    continue

      return
      end
      subroutine test100

c*********************************************************************72
c
cc TEST100 tests PERM_NEXT2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      logical done
      integer iactiv(n)
      integer idir(n)
      integer invers(n)
      integer p(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST100'
      write ( *, '(a)' ) '  PERM_NEXT2 generates permutations in order.'
      write ( *, '(a)' ) ' '
      done = .true.
     
10    continue

        call perm_next2 ( n, p, done, iactiv, idir, invers )
     
        if ( done ) then
          go to 20
        end if

        call perm_print ( n, p, ' ' )

      go to 10

20    continue
     
      return
      end
      subroutine test101

c*********************************************************************72
c
cc TEST101 tests PERM_NEXT3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      logical more
      integer p(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST101'
      write ( *, '(a)' ) '  PERM_NEXT3 generates permutations in order.'
      write ( *, '(a)' ) ' '
      more = .false.
     
10    continue

        call perm_next3 ( n, p, more )

        call perm_print ( n, p, ' ' )

        if ( .not. more ) then
          go to 20
        end if

      go to 10
     
20    continue

      return
      end
      subroutine test102

c*********************************************************************72
c
cc TEST102 tests PERM_RANDOM;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      integer i
      integer p(n)
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST102'
      write ( *, '(a)' ) '  PERM_RANDOM produces a random permutation;'
      write ( *, '(a,i8)' ) '  For this test, N = ', n
      write ( *, '(a)' ) ' '

      seed = 123456789

      do i = 1, 5
        call perm_random ( n, seed, p )
        call perm_print ( n, p, ' ' )
      end do
     
      return
      end
      subroutine test103

c*********************************************************************72
c
cc TEST103 tests PERM_RANDOM2;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      integer i
      integer p(n)
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST103'
      write ( *, '(a)' ) '  PERM_RANDOM2 produces a random permutation '
      write ( *, '(a)' ) '  of labels;'
      write ( *, '(a,i8)' ) '  For this test, N = ', n
      write ( *, '(a)' ) ' '

      seed = 123456789

      do i = 1, 5

        p(1) = 101
        p(2) = 202
        p(3) = 303
        p(4) = 404

        call perm_random2 ( n, seed, p )
        call perm_print ( n, p, ' ' )
      end do
     
      return
      end
      subroutine test104

c*********************************************************************72
c
cc TEST104 tests PERM_RANDOM3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      integer i
      integer p(n)
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST104'
      write ( *, '(a)' ) '  PERM_RANDOM3 produces a random permutation.'
      write ( *, '(a,i8)' ) '  For this test, N = ', n
      write ( *, '(a)' ) ' '
     
      seed = 123456789

      do i = 1, 5
        call perm_random3 ( n, seed, p )
        call perm_print ( n, p, ' ' )
      end do
     
      return
      end
      subroutine test105

c*********************************************************************72
c
cc TEST105 tests PERM_RANK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      integer p(n)
      integer rank

      save p

      data p / 1, 4, 2, 3 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST105'
      write ( *, '(a)' ) '  PERM_RANK ranks a permutation.'

      call perm_print ( n, p, '  The permutation:' )
     
      call perm_rank ( n, p, rank )
     
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The rank is:', rank
     
      return
      end
      subroutine test106

c*********************************************************************72
c
cc TEST106 tests PERM_TO_EQUIV.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 9 )

      integer a(n)
      integer jarray(n)
      integer npart
      integer p(n)

      save p

      data p / 2,3,9,6,7,8,5,4,1 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST106'
      write ( *, '(a)' ) '  PERM_TO_EQUIV returns the set partition'
      write ( *, '(a)' ) '  or equivalence classes determined by a'
      write ( *, '(a)' ) '  permutation.'

      call perm_print ( n, p, '  The input permutation:' )
     
      call perm_to_equiv ( n, p, npart, jarray, a )

      call equiv_print ( n, a, '  The partition:' )
     
      return
      end
      subroutine test107

c*********************************************************************72
c
cc TEST107 tests PERM_TO_YTB.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 7 )

      integer a(n)
      integer lambda(n)
      integer p(n)

      save p

      data p / 7, 2, 4, 1, 5, 3, 6 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST107'
      write ( *, '(a)' ) '  PERM_TO_YTB converts a permutation to a'
      write ( *, '(a)' ) '  Young table.'

      call perm_print ( n, p, '  The permutation:' )
     
      call perm_to_ytb ( n, p, lambda, a )

      call ytb_print ( n, a, '  The Young table:' )
     
      return
      end
      subroutine test108

c*********************************************************************72
c
cc TEST108 tests PERM_UNRANK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      integer p(n)
      integer rank

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST108'
      write ( *, '(a)' ) '  PERM_UNRANK, given a rank, computes the'
      write ( *, '(a)' ) '  corresponding permutation.'
      write ( *, '(a)' ) ' '
      rank = 6
      write ( *, '(a,i8)' ) '  The requested rank is ', rank
     
      call perm_unrank ( n, rank, p )
     
      call perm_print ( n, p, '  The permutation:' )
     
      return
      end
      subroutine test1085

c*********************************************************************72
c
cc TEST1085 tests PERRIN.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 15 )

      integer i
      integer p(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1085'
      write ( *, '(a)' ) '  PERRIN computes the Perrin numbers.'
      write ( *, '(a)' ) ' '

      call perrin ( n, p )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     N        P(N)'
      write ( *, '(a)' ) ' '

      do i = 0, n-1
        write ( *, '(2x,i8,i10)' ) i, p(i+1)
      end do

      return
      end
      subroutine test109

c*********************************************************************72
c
cc TEST109 tests POWER_MOD.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer a
      integer m
      integer n
      integer x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST109'
      write ( *, '(a)' ) '  POWER_MOD computes the remainder of a power'
      write ( *, '(a)' ) '  of an integer modulo another integer.'

      a = 7
      n = 50
      m = 11

      call power_mod ( a, n, m, x )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  A = ', a
      write ( *, '(a,i8)' ) '  N = ', n
      write ( *, '(a,i8)' ) '  M = ', m
      write ( *, '(a,i8)' ) '  mod ( A**N, M ) = ', x

      a = 3
      n = 118
      m = 119

      call power_mod ( a, n, m, x )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  A = ', a
      write ( *, '(a,i8)' ) '  N = ', n
      write ( *, '(a,i8)' ) '  M = ', m
      write ( *, '(a,i8)' ) '  mod ( A**N, M ) = ', x

      return
      end
      subroutine test110

c*********************************************************************72
c
cc TEST110 tests POWER_SERIES1;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      double precision a(n)
      double precision alpha
      double precision b(n)
      integer i

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST110'
      write ( *, '(a)' ) '  POWER_SERIES1 composes a power series;'

      alpha = 7.0D+00
     
      a(1) = 1.0D+00
      do i = 2, n
        a(i) = 0.0D+00
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Power series of G(x) = (1+F(x))**alpha'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  N = ', n
      write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Series for F(x):'
      write ( *, '(a)' ) ' '
      write ( *, '(2x,5g14.6)' ) ( a(i), i = 1, n )
     
      call power_series1 ( n, alpha, a, b )
     
      write ( *, '(a)' ) '  Series for G(x):'
      write ( *, '(a)' ) ' '
      write ( *, '(2x,5g14.6)' ) ( b(i), i = 1, n )
     
      return
      end
      subroutine test111

c*********************************************************************72
c
cc TEST111 tests POWER_SERIES2;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      double precision a(n)
      double precision b(n)
      integer i

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST111'
      write ( *, '(a)' ) '  POWER_SERIES2 composes a power series;'
      write ( *, '(a)' ) '  Here we compute the power series of '
      write ( *, '(a)' ) '  G(x) = exp(F(x))-1'
      write ( *, '(a,i8)' ) '  The number of terms is N = ', n

      a(1) = -4.0D+00
      do i = 2, n
        a(i) = 0.0D+00
      end do
     
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Series for F(x):'
      write ( *, '(a)' ) ' '
      write ( *, '(2x,5g14.6)' ) ( a(i), i = 1, n )
     
      call power_series2 ( n, a, b )
     
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Series for G(x):'
      write ( *, '(a)' ) ' '
      write ( *, '(2x,5g14.6)' ) ( b(i), i = 1, n )
     
      return
      end
      subroutine test112

c*********************************************************************72
c
cc TEST112 tests POWER_SERIES3;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      double precision a(n)
      double precision b(n)
      double precision c(n)
      integer i

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST112'
      write ( *, '(a)' ) '  POWER_SERIES3 composes a power series;'
     
      a(1) = 1.0D+00
      a(2) = 1.0D+00
      do i = 3, n
        a(i) = 0.0D+00
      end do
     
      b(1) = 1.0D+00
      b(2) = 1.0D+00
      do i = 3, n
        b(i) = 0.0D+00
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Power series of H(x) = G(F(x))'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of terms, N = ', n
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Series for F(x):'
      write ( *, '(a)' ) ' '
      write ( *, '(2x,5g14.6)' ) ( a(i), i = 1, n )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Series for G(x):'
      write ( *, '(a)' ) ' '
      write ( *, '(2x,5g14.6)' ) ( b(i), i = 1, n )
     
      call power_series3 ( n, a, b, c )
     
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Series for H(x):'
      write ( *, '(a)' ) ' '
      write ( *, '(2x,5g14.6)' ) ( c(i), i = 1, n )
     
      return
      end
      subroutine test113

c*********************************************************************72
c
cc TEST113 tests POWER_SERIES4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      double precision a(n)
      double precision b(n)
      double precision c(n)
      integer i

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST113'
      write ( *, '(a)' ) '  POWER_SERIES4 composes a power series;'

      do i = 1, n
        a(i) = 1.0D+00 / dble ( i )
      end do

      b(1) = 1.0D+00
      do i = 2, n
        b(i) = 0.0D+00
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Power series of H(x) = G(1/F(x))'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  N = ', n
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Series for F(x):'
      write ( *, '(a)' ) ' '
      write ( *, '(2x,5g14.6)' ) ( a(i), i = 1, n )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Series for G(x):'
      write ( *, '(a)' ) ' '
      write ( *, '(2x,5g14.6)' ) ( b(i), i = 1, n )
     
      call power_series4 ( n, a, b, c )
     
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Series for H(x):'
      write ( *, '(a)' ) ' '
      write ( *, '(2x,5g14.6)' ) ( c(i), i = 1, n )
     
      return
      end
      subroutine test114

c*********************************************************************72
c
cc TEST114 tests PYTHAG_TRIPLE_NEXT
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer a
      integer b
      integer c
      integer d
      integer e
      integer i
      integer j
      integer k

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST114'
      write ( *, '(a)' ) '  PYTHAG_TRIPLE_NEXT computes the "next"'
      write ( *, '(a)' ) '    Pythagorean triple.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   I   J   A   B   C  A^2+B^2   C^2'
      write ( *, '(a)' ) ' '

      i = 0
      j = 0

      do k = 0, 20
        call pythag_triple_next ( i, j, a, b, c )
        d = a**2 + b**2
        e = c**2
        write ( *, '(2x,5i4,2i8)' ) i, j, a, b, c, d, e
      end do

      return
      end
      subroutine test115

c*********************************************************************72
c
cc TEST115 tests R8_AGM
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      integer i4_uniform
      integer j
      double precision r8_agm
      integer seed
      double precision x
      double precision y
      double precision z

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST115'
      write ( *, '(a)' ) '  R8_AGM computes the arithmetic-geometric'
      write ( *, '(a)' ) '  mean (AGM) of two nonnegative real numbers.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    X        Y    R8_AGM(X,Y)'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do i = 1, 10
        j = i4_uniform ( 1, 10, seed )
        x = dble ( j )
        j = i4_uniform ( 1, 10, seed )
        y = dble ( j )
        z = r8_agm ( x, y )
        write ( *, '(2x,f8.4,2x,f8.4,2x,f8.4)' ) x, y, z
      end do

      return
      end
      subroutine test116

c*********************************************************************72
c
cc TEST116 tests R8_TO_CFRAC.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 7 )

      integer a(0:n)
      double precision error
      integer i
      integer p(-1:n)
      integer q(-1:n)
      double precision r
      double precision r8_pi
      double precision temp

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST116'
      write ( *, '(a)' ) '  R8_TO_CFRAC converts a real number to a'
      write ( *, '(a)' ) '  sequence of continued fraction convergents.'

      r = 2.0D+00 * r8_pi ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Use the real number R = ', r

      call r8_to_cfrac ( r, n, a, p, q )

      write ( *, '(a)' ) ' '

      do i = 0, n
        temp = dble ( p(i) ) / dble ( q(i) )
        error = r - temp
        write ( *, '(2x,3i8,2g14.6)' ) a(i), p(i), q(i), temp, error
      end do

      return
      end
      subroutine test1163

c*********************************************************************72
c
cc TEST1163 tests R8_TO_DEC and DEC_TO_R8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer a
      integer b
      integer dec_digit
      integer i
      double precision r
      double precision r2
      double precision r8_uniform_01
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1163'
      write ( *, '(a)' ) '  R8_TO_DEC converts a real to a decimal;'
      write ( *, '(a)' ) '  DEC_TO_R8 converts a decimal to a real.'

      dec_digit = 5

      write ( *, '(a)' ) ' '
      write ( *, '(a,i3)' ) 
     &  '  The number of decimal digits is ', dec_digit

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     R   =>  A * 10^B  =>  R2'
      write ( *, '(a)' ) ' '

      do i = 1, 10
        r = r8_uniform_01 ( seed )
        r = 10.0D+00 * ( r - 0.25D+00 )
        call r8_to_dec ( r, dec_digit, a, b )
        call dec_to_r8 ( a, b, r2 )
        write ( *, '(2x,f10.6,2x,i8,2x,i8,2x,f10.6)' ) r, a, b, r2
      end do

      return
      end
      subroutine test1165

c*********************************************************************72
c
cc TEST1165 tests R8_TO_RAT and RAT_TO_R8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer a
      integer b
      integer i
      integer ndig
      double precision r
      double precision r2
      double precision r8_uniform_01
      integer seed

      ndig = 4

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1165'
      write ( *, '(a)' ) '  R8_TO_RAT converts a real to a rational;'
      write ( *, '(a)' ) '  RAT_TO_R8 converts a rational to a real.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i3)' ) 
     &  '  The maximum number of digits allowed is ', ndig

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     R   =>  A / B  =>  R2'
      write ( *, '(a)' ) ' '

      do i = 1, 10
        r = r8_uniform_01 ( seed )
        r = 10.0D+00 * ( r - 0.25D+00 )
        call r8_to_rat ( r, ndig, a, b )
        call rat_to_r8 ( a, b, r2 )
        write ( *, '(2x,f10.6,i8,2x,i8,f10.6)' ) r, a, b, r2
      end do

      return
      end
      subroutine test117

c*********************************************************************72
c
cc TEST117 tests RAT_ADD.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer abot
      integer atop
      integer bbot
      integer btop
      integer cbot
      integer ctop
      integer ierror
      character * ( 22 ) string
      integer string_length

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST117'
      write ( *, '(a)' ) '  RAT_ADD adds two rationals.'

      atop = 3
      abot = 4
      btop = 10
      bbot = 7

      call rat_add ( atop, abot, btop, bbot, ctop, cbot, ierror )

      write ( *, '(a)' ) ' '
      call rat_to_s_left ( atop, abot, string )
      string_length = len_trim ( string )
      write ( *, '(a,a)' ) '  A = ', string(1:string_length)

      call rat_to_s_left ( btop, bbot, string )
      string_length = len_trim ( string )
      write ( *, '(a,a)' ) '  B = ', string(1:string_length)

      call rat_to_s_left ( ctop, cbot, string )
      string_length = len_trim ( string )
      write ( *, '(a,a)' ) '  C = A + B = ', string(1:string_length)
     
      return
      end
      subroutine test118

c*********************************************************************72
c
cc TEST118 tests RAT_DIV.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer abot
      integer atop
      integer bbot
      integer btop
      integer cbot
      integer ctop
      integer ierror
      character * ( 22 ) string
      integer string_length

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST118'
      write ( *, '(a)' ) '  RAT_DIV divides two rationals.'

      atop = 3
      abot = 4
      btop = 10
      bbot = 7

      call rat_div ( atop, abot, btop, bbot, ctop, cbot, ierror )

      write ( *, '(a)' ) ' '
      call rat_to_s_left ( atop, abot, string )
      string_length = len_trim ( string )
      write ( *, '(a,a)' ) '  A = ', string(1:string_length)

      call rat_to_s_left ( btop, bbot, string )
      string_length = len_trim ( string )
      write ( *, '(a,a)' ) '  B = ', string(1:string_length)

      call rat_to_s_left ( ctop, cbot, string )
      string_length = len_trim ( string )
      write ( *, '(a,a)' ) '  C = A / B = ', string(1:string_length)
     
      return
      end
      subroutine test119

c*********************************************************************72
c
cc TEST119 tests RAT_FAREY.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer max_frac
      parameter ( max_frac = 20 )

      integer a(max_frac)
      integer b(max_frac)
      integer i
      integer ihi
      integer ilo
      integer n
      integer num_frac

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST119'
      write ( *, '(a)' ) '  RAT_FAREY computes a row of '
      write ( *, '(a)' ) '  the Farey fraction table.'

      do n = 1, 7

        call rat_farey ( n, max_frac, num_frac, a, b )
     
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Row ', n
        write ( *, '(a,i8)' ) '  Number of fractions: ', num_frac

        do ilo = 1, num_frac, 20
          ihi = min ( ilo+20-1, num_frac )
          write ( *, '(a)' ) ' '
          write ( *, '(2x,20i3)' ) ( a(i), i = ilo, ihi )
          write ( *, '(2x,20i3)' ) ( b(i), i = ilo, ihi )
        end do

      end do

      return
      end
      subroutine test120

c*********************************************************************72
c
cc TEST120 tests RAT_FAREY2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit  none

      integer max_n
      parameter ( max_n = 4 )

      integer a(2**max_n+1)
      integer b(2**max_n+1)
      integer i
      integer ihi
      integer ilo
      integer n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST120'
      write ( *, '(a)' ) '  RAT_FAREY2 computes a row of'
      write ( *, '(a)' ) '  the Farey fraction table.'

      do n = 0, max_n

        call rat_farey2 ( n, a, b )
     
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Row ', n+1

        do ilo = 1, 2**n+1, 20
          ihi = min ( ilo+20-1, 2**n+1 )
          write ( *, '(a)' ) ' '
          write ( *, '(2x,20i3)' ) ( a(i), i = ilo, ihi )
          write ( *, '(2x,20i3)' ) ( b(i), i = ilo, ihi )
        end do

      end do

      return
      end
      subroutine test121

c*********************************************************************72
c
cc TEST121 tests RAT_MUL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer abot
      integer atop
      integer bbot
      integer btop
      integer cbot
      integer ctop
      integer ierror
      character * ( 22 ) string
      integer string_length

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST121'
      write ( *, '(a)' ) '  RAT_MUL multiplies two rationals.'

      atop = 3
      abot = 4
      btop = 10
      bbot = 7

      call rat_mul ( atop, abot, btop, bbot, ctop, cbot, ierror )

      write ( *, '(a)' ) ' '
      call rat_to_s_left ( atop, abot, string )
      string_length = len_trim ( string )
      write ( *, '(a,a)' ) '  A = ', string(1:string_length)

      call rat_to_s_left ( btop, bbot, string )
      string_length = len_trim ( string )
      write ( *, '(a,a)' ) '  B = ', string(1:string_length)

      call rat_to_s_left ( ctop, cbot, string )
      string_length = len_trim ( string )
      write ( *, '(a,a)' ) '  C = A * B = ', string(1:string_length)
     
      return
      end
      subroutine test1215

c*********************************************************************72
c
cc TEST1215 tests RAT_WIDTH.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_test
      parameter ( n_test = 17 )

      integer a
      integer a_test(n_test)
      integer b
      integer b_test(n_test)
      integer i
      integer rat_width
      integer width

      data a_test /
     &  1000, 1000, 1000, 1000, 1000, 1, -1, -10, -100, -1000, 
     &  1, 10, 100, 1000, 10000, 17, 4000000 /
      data b_test /
     &  3, 40, 500, 6000, 70000, 1, 200, 200, 200, 200, 
     & -200, -200, -200, -200, -200, 3000, 4000000 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1215'
      write ( *, '(a)' ) '  RAT_WIDTH determines the "width"'
      write ( *, '(a)' ) '  of a rational.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     Top    Bottom  Width'
      write ( *, '(a)' ) ' '

      do i = 1, n_test
        a = a_test(i)
        b = b_test(i)
        width = rat_width ( a, b )
        write ( *, '(2x,3i8)' ) a, b, width
      end do

      return
      end
      subroutine test122

c*********************************************************************72
c
cc TEST122 tests RATMAT_DET.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 6 )

      integer a(0:n,n+1)
      integer b(0:n,n+1)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST122'
      write ( *, '(a)' ) '  RAT_SUM_FORMULA computes the coefficients'
      write ( *, '(a)' ) '  for the formulas for the sums of powers '
      write ( *, '(a)' ) '  of integers.'
      
      call rat_sum_formula ( n, a, b )

      call ratmat_print ( n+1, n+1, a, b, '  Power Sum Coefficients:' )

      return
      end
      subroutine test123

c*********************************************************************72
c
cc TEST123 tests RATMAT_DET.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n3
      parameter ( n3 = 3 )

      integer a3(n3,n3)
      integer b3(n3,n3)
      integer i
      integer idbot
      integer idtop
      integer ierror
      integer j
      integer k

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST123'
      write ( *, '(a)' ) '  RATMAT_DET:'
      write ( *, '(a)' ) '    determinant of a rational matrix.'
      write ( *, '(a)' ) ' '
     
      k = 0
      do i = 1, n3
        do j = 1, n3
          k = k + 1
          a3(i,j) = k
        end do
      end do

      do i = 1, n3
        do j = 1, n3
          b3(i,j) = 1
        end do
      end do
     
      call ratmat_print ( n3, n3, a3, b3, '  The 123/456/789 matrix:' )

      call ratmat_det ( n3, a3, b3, idtop, idbot, ierror )
     
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Determinant of the 123/456/789 matrix:'
      write ( *, '(2x,i8,a,i8)' ) idtop, ' / ', idbot
     
      do i = 1, n3
        do j = 1, n3
          a3(i,j) = 1
          b3(i,j) = i + j
        end do
      end do
     
      call ratmat_print ( n3, n3, a3, b3, '  The Hilbert matrix:' )

      call ratmat_det ( n3, a3, b3, idtop, idbot, ierror )
     
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Determinant of the Hilbert matrix:'
      write ( *, '(2x,i8,a,i8)' ) idtop, ' / ', idbot
     
      do i = 1, n3
        do j = 1, n3
          if ( i .eq. j ) then
            a3(i,j) = 2
          else if ( i .eq. j+1 .or. i .eq. j-1 ) then
            a3(i,j) = -1
          else
            a3(i,j) = 0
          end if
          b3(i,j) = 1
        end do
      end do
     
      call ratmat_print ( n3, n3, a3, b3, '  The -1 2 -1 matrix:' )

      call ratmat_det ( n3, a3, b3, idtop, idbot, ierror )
     
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Determinant of the -1,2,-1 matrix:'
      write ( *, '(2x,i8,a,i8)' ) idtop, ' / ', idbot
     
      return
      end
      subroutine test124

c*********************************************************************72
c
cc TEST124 tests REGRO_NEXT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      logical done
      integer i
      integer rank
      integer v(n)
      integer vmax(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST124'
      write ( *, '(a)' ) '  REGRO_NEXT generates all restricted growth '
      write ( *, '(a)' ) '  functions.'
      write ( *, '(a)' ) ' '

      rank = 0

      done = .true.
     
10    continue

        call regro_next ( n, v, vmax, done )

        if ( done ) then
          go to 20
        end if

        rank = rank + 1
        write ( *, '(2x,5i3)' ) rank, ( v(i), i = 1, n )

      go to 10

20    continue
     
      return
      end
      subroutine test1245

c*********************************************************************72
c
cc TEST1245 tests R8_RISE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      double precision r8_rise
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1245'
      write ( *, '(a)' ) 
     &  '  R8_RISE computes the rising factorial function.'
      write ( *, '(a)' ) 
     &  '  [X]_N = X * (X+1) * (X+2) * ... * ( X+N-1).'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    X          N  R8_RISE(X,N)'
      write ( *, '(a)' ) ' '

      x = 4.0D+00

      do n = -2, 5
     
        write ( *, '(2x,f10.6,i8,f12.6)' ) x, n, r8_rise ( x, n )
     
      end do
     
      return
      end
      subroutine test125

c*********************************************************************72
c
cc TEST125 tests R8MAT_DET.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n3
      parameter ( n3 = 3 )
      integer n4
      parameter ( n4 = 4 )

      double precision a3(n3,n3)
      double precision a4(n4,n4)
      double precision det
      integer i
      integer j
      integer k

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST125'
      write ( *, '(a)' ) '  R8MAT_DET: determinant of a real matrix.'
      write ( *, '(a)' ) ' '
     
      k = 0
      do i = 1, n3
        do j = 1, n3
          k = k+1
          a3(i,j) = dble ( k )
        end do
      end do
     
      call r8mat_print ( n3, n3, a3, '  The 123/456/789 matrix:' )

      call r8mat_det ( n3, a3, det )
     
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) 
     &  '  Determinant of the 123/456/789 matrix is ', det
     
      do i = 1, n4
        do j = 1, n4
          a4(i,j) = 1.0D+00 / dble ( i + j )
        end do
      end do
     
      call r8mat_print ( n4, n4, a4, '  The Hilbert matrix:' )

      call r8mat_det ( n4, a4, det )
     
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) 
     &  '  Determinant of the Hilbert matrix is ', det
     
      do i = 1, n3
        do j = 1, n3
          if ( i .eq. j ) then
            a3(i,j) = 2.0D+00
          else if ( i .eq. j+1 .or. i .eq. j-1 ) then
            a3(i,j) = -1.0D+00
          else
            a3(i,j) = 0.0D+00
          end if
        end do
      end do
     
      call r8mat_print ( n3, n3, a3, '  The -1,2,-1 matrix:' )

      call r8mat_det ( n3, a3, det )
     
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) 
     &  '  Determinant of the -1,2,-1 matrix is ', det
     
      return
      end
      subroutine test126

c*********************************************************************72
c
cc TEST126 tests R8MAT_PERM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 9 )

      double precision a(n,n)
      integer i
      integer p(n)
      integer j

      save p

      data p / 2, 3, 9, 6, 7, 8, 5, 4, 1 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST126'
      write ( *, '(a)' ) '  R8MAT_PERM reorders a real matrix in place.'
      write ( *, '(a)' ) 
     &  '  The rows and columns use the same permutation.'
     
      do i = 1, n
        do j = 1, n
          a(i,j) = dble ( i * 10 + j )
        end do
      end do
     
      call r8mat_print ( n, n, a, '  The original matrix' )
     
      call perm_print ( n, p, '  The row and column permutation:' )
     
      call r8mat_perm ( n, a, p )
     
      call r8mat_print ( n, n, a, '  The permuted matrix' )
     
      return
      end
      subroutine test127

c*********************************************************************72
c
cc TEST127 tests R8MAT_PERM2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 9 )
      integer n
      parameter ( n = 7 )

      double precision a(m,n)
      integer i
      integer j
      integer p(m)
      integer q(n)

      save p
      save q

      data p / 2, 3, 9, 6, 7, 8, 5, 4, 1 /
      data q / 3, 4, 5, 6, 7, 1, 2 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST127'
      write ( *, '(a)' ) 
     &  '  R8MAT_PERM2 reorders a real matrix in place.'
      write ( *, '(a)' ) 
     &  '  Rows and columns use different permutations.'
     
      do i = 1, m
        do j = 1, n
          a(i,j) = dble ( i * 10 + j )
        end do
      end do
     
      call r8mat_print ( m, n, a, '  The original matrix' )
     
      call perm_print ( m, p, '  The row permutation:' )
     
      call perm_print ( n, q, '  The column permutation:' )

      call r8mat_perm2 ( m, n, a, p, q )
     
      call r8mat_print ( m, n, a, '  The permuted matrix' )
     
      return
      end
      subroutine test128

c*********************************************************************72
c
cc TEST128 tests R8MAT_PERMANENT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 12 )

      double precision a(n_max*n_max)
      integer i
      integer j
      integer n
      double precision perm

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST128'
      write ( *, '(a)' ) 
     &  '  R8MAT_PERMANENT: the matrix permanent function.'
      write ( *, '(a)' ) '  We will analyze matrices with 0 diagonal'
      write ( *, '(a)' ) '  and 1 on all offdiagonals.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Order	    Permanent.'
      write ( *, '(a)' ) ' '
     
      do n = 2, n_max
     
        do j = 1, n
          do i = 1, n
            if ( i .eq. j ) then
              a(i+(j-1)*n) = 0.0D+00
            else
              a(i+(j-1)*n) = 1.0D+00
            end if
          end do
        end do
     
        call r8mat_permanent ( n, a, perm )
     
        write ( *, '(7x,i2,8x,g18.10)' ) n, perm
     
      end do
     
      return
      end
      subroutine test129

c*********************************************************************72
c
cc TEST129 test R8POLY.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 6 )

      double precision a(n)
      integer i
      integer iopt
      integer test
      double precision val
      double precision x0

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST129'
      write ( *, '(a)' ) '  R8POLY converts between power sum,'
      write ( *, '(a)' ) '  factorial, and Taylor forms, '
      write ( *, '(a)' ) '  and can evaluate a polynomial'
      write ( *, '(a)' ) ' '
     
      do test = 1, 6

        if ( test .eq. 1 ) then
          iopt = -3
        else if ( test .eq. 2 ) then
          iopt = -2
        else if ( test .eq. 3 ) then
          iopt = -1
          x0 = 2.0D+00
        else if ( test .eq. 4 ) then
          iopt = 0
          x0 = 2.0D+00
        else if ( test .eq. 5 ) then
          iopt = 6
          x0 = 2.0D+00
        else if ( test .eq. 6 ) then
          iopt = 6
          x0 = -2.0D+00
        end if

        do i = 1, n-1
          a(i) = 0.0D+00
        end do
        a(n) = 1.0D+00

        if ( test .eq. 1 ) then
          write ( *, '(a)' ) '  All calls have input A as follows'
          write ( *, '(2x,6f7.2)' ) ( a(i), i = 1, n )
        end if
     
        call r8poly ( n, a, x0, iopt, val )
     
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Option IOPT = ', iopt
        if ( -1 .le. iopt ) then
          write ( *, '(a,g14.6)' ) '  X0 = ', x0
        end if

        if ( iopt .eq. -3 .or. iopt .eq. -2 .or. 0 .lt. iopt ) then
          write ( *, '(a)' ) '  Output array = '
          write ( *, '(2x,6f7.2)' ) ( a(i), i = 1, n )
        end if

        if ( iopt .eq. -1 .or. iopt .eq. 0 ) then
          write ( *, '(a,g14.6)' ) '  Value = ', val
        end if
     
      end do

      return
      end
      subroutine test130

c*********************************************************************72
c
cc TEST130 tests R8POLY_DIV.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a(0:10)
      double precision b(0:10)
      integer na
      integer nb
      integer nq
      integer nr
      double precision q(0:10)
      double precision r(0:10)
      integer test
      integer test_num

      test_num = 2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST130'
      write ( *, '(a)' ) '  R8POLY_DIV computes the quotient and'
      write ( *, '(a)' ) '  remainder for polynomial division.'
      write ( *, '(a)' ) ' '
c
c  1: Divide X**3 + 2*X**2 - 5*X - 6  by X-2.  
c     Quotient is 3+4*X+X**2, remainder is 0.
c
c  2: Divide X**4 + 3*X**3 + 2*X**2 - 2  by  X**2 + X - 3.
c     Quotient is X**2 + 2*X + 3, remainder 8*X + 7.
c
      do test = 1, test_num

        if ( test .eq. 1 ) then
          na = 3
          a(0) = -6.0D+00
          a(1) = -5.0D+00
          a(2) =  2.0D+00
          a(3) =  1.0D+00
          nb = 1
          b(0) = -2.0D+00
          b(1) =  1.0D+00
        else if ( test .eq. 2 ) then
          na = 4
          a(0) = -2.0D+00
          a(1) =  5.0D+00
          a(2) =  2.0D+00
          a(3) =  3.0D+00
          a(3) =  1.0D+00
          nb = 2
          b(0) = -3.0D+00
          b(1) =  1.0D+00
          b(2) =  1.0D+00
        end if

        call r8poly_print ( na, a,
     &    '  The polynomial to be divided, A:' )
        call r8poly_print ( nb, b, '  The divisor polynomial, B:' )

        call r8poly_div ( na, a, nb, b, nq, q, nr, r )
     
        call r8poly_print ( nq, q, '  The quotient polynomial, Q:' )
        call r8poly_print ( nr, r, '  The remainder polynomial, R:' )

      end do

      return
      end
      subroutine test131

c*********************************************************************72
c
cc TEST131 tests R8POLY_F2P and R8POLY_P2F.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      double precision a(n)

      call r8vec_indicator ( n, a )
     
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST131'
      write ( *, '(a)' ) '  R8POLY_P2F: power sum => factorial;'
      write ( *, '(a)' ) '  R8POLY_F2P: factorial => power sum.'

      call r8poly_print ( n-1, a, '  The power sum polynomial:' )
     
      call r8poly_p2f ( n, a )
     
      call r8vec_print ( n, a, 
     &  '  The factorial polynomial coefficients:' )
     
      call r8poly_f2p ( n, a )
     
      call r8poly_print ( n-1, a, 
     &  '  The recovered power sum polynomial:' )
     
      return
      end
      subroutine test132

c*********************************************************************72
c
cc TEST132 tests R8POLY_FVAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      double precision a(n)
      double precision val
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST132'
      write ( *, '(a)' ) '  R8POLY_FVAL evaluates a polynomial'
      write ( *, '(a)' ) '  in factorial form.'

      call r8vec_indicator ( n, a )
     
      call r8vec_print ( n, a, 
     &  '  The factorial polynomial coefficients:' )

      x = 2.0D+00

      call r8poly_fval ( n, a, x, val )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,g14.6)' ) '  R8POLY (', x, ' ) = ', val
      write ( *, '(a)' ) '  The correct value is 11.'
     
      return
      end
      subroutine test133

c*********************************************************************72
c
cc TEST133 tests R8POLY_MUL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer maxn
      parameter ( maxn = 5 )

      double precision a(0:maxn)
      double precision b(0:maxn)
      double precision c(0:maxn)
      integer na
      integer nb
      integer test
      integer test_num

      test_num = 2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST133'
      write ( *, '(a)' ) '  R8POLY_MUL multiplies two polynomials.'
      write ( *, '(a)' ) ' '
c
c  1: Multiply (1+X) times (1-X).  Answer is 1-X**2.
c  2: Multiply (1+2*X+3*X**2) by (1-2*X). Answer is 1 + 0*X - X**2 - 6*X**3
c
      do test = 1, test_num

        if ( test .eq. 1 ) then
          na = 1
          a(0) = 1.0D+00
          a(1) = 1.0D+00
          nb = 1
          b(0) =  1.0D+00
          b(1) = -1.0D+00
        else if ( test .eq. 2 ) then
          na = 2
          a(0) = 1.0D+00
          a(1) = 2.0D+00
          a(2) = 3.0D+00
          nb = 1
          b(0) =  1.0D+00
          b(1) = -2.0D+00
        end if

        call r8poly_mul ( na, a, nb, b, c )

        call r8poly_print ( na, a, '  The factor A:' )

        call r8poly_print ( nb, b, '  The factor B:' )

        call r8poly_print ( na+nb, c, '  The product C = A*B:' )

      end do

      return
      end
      subroutine test134

c*********************************************************************72
c
cc TEST134 tests R8POLY_N2P and R8POLY_P2N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      double precision a(n)
      double precision a2(n)
      integer i

      call r8vec_indicator ( n, a )

      do i = 1, n
        a2(i) = 2.0D+00 * a(i)
      end do
     
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST134'
      write ( *, '(a)' ) '  R8POLY_N2P: Newton => power sum;'
      write ( *, '(a)' ) '  R8POLY_P2N: Power sum => Newton.'

      call r8poly_print ( n-1, a, '  The power sum polynomial:' )
     
      call r8poly_p2n ( n, a, a2 )
     
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Newton polynomial coefficients:'
      write ( *, '(a)' ) ' '
      write ( *, '(6f12.4)' ) ( a(i), i = 1, n )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Newton polynomial abscissas:'
      write ( *, '(a)' ) ' '
      write ( *, '(6f12.4)' ) ( a2(i), i = 1, n )
     
      call r8poly_n2p ( n, a, a2 )
     
      call r8poly_print ( n-1, a, 
     &  '  The recovered power sum polynomial:' )

      return
      end
      subroutine test135

c*********************************************************************72
c
cc TEST135 tests R8POLY_NVAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      double precision a(n)
      double precision a2(n-1)
      integer i
      double precision val
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST135'
      write ( *, '(a)' ) '  R8POLY_NVAL evaluates a Newton polynomial.'

      call r8vec_indicator ( n, a )

      do i = 1, n-1
        a2(i) = a(i) - 1.0D+00
      end do
     
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Newton polynomial coefficients:'
      write ( *, '(a)' ) ' '
      write ( *, '(2x,6f12.4)' ) ( a(i), i = 1, n )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Newton polynomial abscissas:'
      write ( *, '(a)' ) ' '
      write ( *, '(2x,6f12.4)' ) ( a2(i), i = 1, n )
     
      x = 2.0D+00
     
      call r8poly_nval ( n, a, a2, x, val )
     
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,g14.6)' ) '  R8POLY ( ', x,' ) = ', val
      write ( *, '(a)' ) '  The correct value is 11.'
     
      return
      end
      subroutine test136

c*********************************************************************72
c
cc TEST136 tests R8POLY_P2T and R8POLY_T2P.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      double precision a(n)
      integer i
      double precision x

      call r8vec_indicator ( n, a )

      x = 2.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST136'
      write ( *, '(a)' ) '  R8POLY_T2P: Taylor => Power sum;'
      write ( *, '(a)' ) '  R8POLY_P2T: Power sum => Taylor.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Taylor expansion point is X = ', x
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The Taylor coefficients:'
      write ( *, '(a)' ) ' '
      write ( *, '(2x,6f12.4)' ) ( a(i), i = 1, n )

      call r8poly_t2p ( n, a, x )

      call r8poly_print ( n-1, a, '  The power sum polynomial:' )

      call r8poly_p2t ( n, a, x )
     
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The recovered Taylor coefficients:'
      write ( *, '(a)' ) ' '
      write ( *, '(2x,6f12.4)' ) ( a(i), i = 1, n )
     
      return
      end 
      subroutine test137

c*********************************************************************72
c
cc TEST137 tests R8POLY_POWER.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer lmax
      parameter ( lmax = 10 )

      double precision a(0:lmax)
      double precision b(0:10)
      integer na
      integer p

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST137'
      write ( *, '(a)' ) 
     &  '  R8POLY_POWER takes a polynomial to a power.'
c
c  Cube (2-X).  Answer is 8-12*X+6*X**2-X**3.
c
      na = 1
      a(0) =  2.0D+00
      a(1) = -1.0D+00
      p = 3

      call r8poly_print ( na, a, '  The polynomial A:' )
     
      call r8poly_power ( na, a, p, b )
     
      call r8poly_print ( p*na, b, '  Raised to the power 3:' )
c
c  Square X+X**2
c
      na = 2

      a(0) = 0.0D+00
      a(1) = 1.0D+00
      a(2) = 1.0D+00

      p = 2

      call r8poly_print ( na, a, '  The polynomial A:' )
     
      call r8poly_power ( na, a, p, b )
     
      call r8poly_print ( p*na, b, '  Raised to the power 2:' )
     
      return
      end
      subroutine test138

c*********************************************************************72
c
cc TEST138 tests R8POLY_PVAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      integer i
      double precision a(0:n)
      double precision val
      double precision x

      do i = 0, n
        a(i) = dble ( i + 1 )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST138'
      write ( *, '(a)' ) '  R8POLY_PVAL evaluates a polynomial'
      write ( *, '(a)' ) '  in power sum form.'

      call r8poly_print ( n, a, '  The polynomial to be evaluated:' )

      x = 2.0D+00
     
      call r8poly_pval ( n, a, x, val )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  At X = ', x
      write ( *, '(a,g14.6)' ) '  Computed polynomial value is ', val
      write ( *, '(a)' ) '  Correct value is 129.'
     
      return
      end
      subroutine test139

c*********************************************************************72
c
cc TEST139 tests R8VEC_FRAC.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      double precision a(n)
      double precision ahi
      double precision alo
      double precision afrac
      integer k
      integer seed

      ahi = 10.0D+00
      alo = 0.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST139'
      write ( *, '(a)' ) 
     &  '  R8VEC_FRAC: K-th smallest real vector entry;'

      seed = 123456789

      call r8vec_uniform ( n, alo, ahi, seed, a )

      call r8vec_print ( n, a, '  The real array to search: ' )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Frac   R8VEC_FRAC'
      write ( *, '(a)' ) ' '

      do k = 1, n

        call r8vec_frac ( n, a, k, afrac )
        write ( *, '(2x,i4,2x,g14.6)' ) k, afrac

      end do

      return
      end
      subroutine test1395

c*********************************************************************72
c
cc TEST1395 tests R8VEC_MIRROR_NEXT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      double precision a(n)
      logical done

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1395'
      write ( *, '(a)' ) '  R8VEC_MIRROR_NEXT generates all sign'
      write ( *, '(a)' ) '  variations of a real vector.'

      a(1) = 1.0D+00
      a(2) = 2.0D+00
      a(3) = 3.0D+00

10    continue

        call r8vec_print ( n, a, '  Next vector:' )

        call r8vec_mirror_next ( n, a, done )

        if ( done ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Done.'
          go to 20
        end if

      go to 10

20    continue

      a(1) = 1.0D+00
      a(2) = 0.0D+00
      a(3) = 3.0D+00

30    continue

        call r8vec_print ( n, a, '  Next vector:' )

        call r8vec_mirror_next ( n, a, done )

        if ( done ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Done.'
          go to 40
        end if

      go to 30

40    continue

      return
      end
      subroutine test140

c*********************************************************************72
c
cc TEST140 tests SCHROEDER.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer i
      integer s(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST140'
      write ( *, '(a)' ) '  SCHROEDER computes the Schroeder numbers.'
      write ( *, '(a)' ) ' '

      call schroeder ( n, s )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N        S(N)'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,i8,i10)' ) i, s(i)
      end do

      return
      end
      subroutine test141

c*********************************************************************72
c
cc TEST141 tests SORT_HEAP_EXTERNAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 20 )

      integer a(n)
      integer i
      integer indx
      integer isgn
      integer j
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST141'
      write ( *, '(a)' ) 
     &  '  SORT_HEAP_EXTERNAL sorts objects externally.'
      write ( *, '(a)' ) ' '

      indx = 0
      i = 0
      j = 0
      isgn = 0
      seed = 123456789

      call i4vec_uniform ( n, 1, n, seed, a )
     
      call i4vec_print ( n, a, '  Unsorted array:' )
     
10    continue

        call sort_heap_external ( n, indx, i, j, isgn )
     
        if ( indx .lt. 0 ) then
          isgn = 1
          if ( a(i) .le. a(j) ) then
            isgn = -1
          end if
        else if ( 0 .lt. indx ) then
          call i4_swap ( a(i), a(j) )
        else
          go to 20
        end if

      go to 10

20    continue

      call i4vec_print ( n, a, '  Sorted array:' )
     
      return
      end
      subroutine test142

c*********************************************************************72
c
cc TEST142 tests SUB_BY_SIZE_NEXT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      integer a(n)
      integer i
      logical more
      integer rank
      integer size

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST142'
      write ( *, '(a)' ) '  SUB_BY_SIZE_NEXT generates all subsets'
      write ( *, '(a)' ) '  of an N set.'
      write ( *, '(a)' ) ' '

      more = .false.
      rank = 0

10    continue

        call sub_by_size_next ( n, a, size, more )

        rank = rank + 1

        if ( 0 .lt. size ) then
          write ( *, '(2x,i4,4x,5i2)' ) rank, ( a(i), i = 1, size )
        else
          write ( *, '(2x,i4,4x,a)' ) rank, 'The empty set'
        end if

        if ( .not. more ) then
          go to 20
        end if

      go to 10

20    continue

      return
      end
      subroutine test143

c*********************************************************************72
c
cc TEST143 tests SUB_LEX_NEXT without size restrictions.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      integer a(n)
      integer i
      integer k
      logical ltest
      integer ndim
      integer rank

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST143'
      write ( *, '(a)' ) '  SUB_LEX_NEXT generates all subsets'
      write ( *, '(a)' ) '  of an N set.'
      write ( *, '(a)' ) ' '

      ndim = n
      k = 0
      rank = 0
      ltest = .false.
     
10    continue
     
        call sub_lex_next ( n, ltest, ndim, k, a )
     
        rank = rank + 1

        if ( 0 .lt. k ) then
          write ( *, '(2x,i4,4x,6i2)' ) rank, ( a(i), i = 1, k )
        else
          write ( *, '(2x,i4,4x,a)' ) rank, 'The empty set.'
        end if
     
        if ( k .eq. 0 ) then
          go to 20
        end if

      go to 10

20    continue
     
      return
      end
      subroutine test1435

c*********************************************************************72
c
cc TEST1435 tests SUB_LEX_NEXT with size restrictions.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer ndim
      parameter ( ndim = 3 )

      integer a(ndim)
      integer i
      integer k
      logical ltest
      integer n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1435'
      write ( *, '(a)' ) '  SUB_LEX_NEXT generates all subsets'
      write ( *, '(a)' ) '  of an N set.'
      write ( *, '(a)' ) '  The user can impose a restriction on'
      write ( *, '(a)' ) '  the maximum size of the subsets.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Here, we require the subsets to be no'
      write ( *, '(a,i8)' ) '  larger than ', ndim

      n = 5
      k = 0
     
10    continue
     
        ltest = ( k .eq. ndim )

        call sub_lex_next ( n, ltest, ndim, k, a )
     
        if ( 0 .lt. k ) then
          write ( *, '(2x,6i2)' ) ( a(i), i = 1, k )
        else
          write ( *, '(a)' ) '  The empty set.'
        end if
     
        if ( k .eq. 0 ) then
          go to 20
        end if

      go to 10

20    continue
     
      return
      end
      subroutine test144

c*********************************************************************72
c
cc TEST144 tests SUB_GRAY_NEXT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      integer a(n)
      integer i
      integer iadd
      logical more
      integer ncard
      integer rank

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST144'
      write ( *, '(a)' ) '  SUB_GRAY_NEXT generates all subsets of an'
      write ( *, '(a)' ) '  N set using the Gray code ordering:'
      write ( *, '(a)' ) '  0 0 1 0 1 means the subset has 3 and 5.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Gray code'
      write ( *, '(a)' ) ' '
     
      rank = 0
      more = .false.
     
10    continue
     
        call sub_gray_next ( n, a, more, ncard, iadd )

        rank = rank + 1 
        write ( *, '(2x,i4,4x,5i2)' ) rank, ( a(i), i = 1, n )

        if ( .not. more ) then
          go to 20
        end if

      go to 10

20    continue

      return
      end
      subroutine test145

c*********************************************************************72
c
cc TEST145 tests SUB_RANDOM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n 
      parameter ( n = 5 )

      integer a(n)
      integer i
      integer j
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST145'
      write ( *, '(a)' ) '  SUB_RANDOM picks a subset at random.'
      write ( *, '(a,i8)' ) '  The number of elements available is ', n
      write ( *, '(a)' ) ' '

      seed = 123456789

      do j = 1, 5
        call sub_random ( n, seed, a )
        write ( *, '(2x,40i2)' ) ( a(i), i = 1, n )
      end do
     
      return
      end
      subroutine test146

c*********************************************************************72
c
cc TEST146 tests SUB_GRAY_RANK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      integer a(n)
      integer i
      integer rank

      save a

      data a / 1, 0, 1, 1, 0 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST146'
      write ( *, '(a)' ) '  SUB_GRAY_RANK returns rank of a subset of'
      write ( *, '(a)' ) '  an N set using the Gray code ordering.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  For N = ', n
      write ( *, '(a)' ) '  the subset is:'
      write ( *, '(2x,5i2)' ) ( a(i), i = 1, n )
     
      call sub_gray_rank ( n, a, rank )
     
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The rank is ', rank
     
      return
      end
      subroutine test147

c*********************************************************************72
c
cc TEST147 tests SUB_GRAY_UNRANK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      integer a(n)
      integer i
      integer rank

      rank = 8
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST147'
      write ( *, '(a)' ) 
     &  '  SUB_GRAY_UNRANK finds the subset of an N set'
      write ( *, '(a)' ) 
     &  '  of a given rank under the Gray code ordering.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  N is ', n
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Rank   Subset'
      write ( *, '(a)' ) ' '

      do rank = 1, 10
     
        call sub_gray_unrank ( rank, n, a )

        write ( *, '(2x,i4,4x,5i2)' ) rank, ( a(i), i = 1, n )

      end do
     
      return
      end
      subroutine test1475

c*********************************************************************72
c
cc TEST1475 tests SUBCOMP_NEXT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 July 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer k
      parameter ( k = 3 )

      integer a(k)
      integer count
      integer h
      integer i
      logical more
      integer n
      integer t
      integer total

      n = 6

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1475'
      write ( *, '(a)' ) '  SUBCOMP_NEXT generates subcompositions.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Seek all subcompositions of N = ', n
      write ( *, '(a,i8,a)' ) '  using K = ', k, ' parts.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     #   Sum'
      write ( *, '(a)' ) ' '

      more = .false.
      count = 0

10    continue

        call subcomp_next ( n, k, a, more, h, t )

        count = count + 1
        call i4vec_sum ( k, a, total )
        write ( *, '(2x,i4,2x,i4,2x,8i4)' ) 
     &    count, total, ( a(i), i = 1, k)

        if ( .not. more )  then
          go to 20
        end if

      go to 10

20    continue
     
      return
      end
      subroutine test1476

c*********************************************************************72
c
cc TEST1476 tests SUBCOMPNZ_NEXT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer k
      parameter ( k = 3 )

      integer a(k)
      integer count
      integer i
      logical more
      integer n
      integer total

      n = 6

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1476'
      write ( *, '(a)' ) '  SUBCOMPNZ_NEXT generates subcompositions'
      write ( *, '(a)' ) '  using nonzero parts.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Seek all subcompositions of N = ', n
      write ( *, '(a,i8,a)' ) '  using K = ', k, ' nonzero parts.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     #   Sum'
      write ( *, '(a)' ) ' '

      more = .false.
      count = 0

10    continue

        call subcompnz_next ( n, k, a, more )

        count = count + 1
        call i4vec_sum ( k, a, total )
        write ( *, '(2x,i4,2x,i4,2x,8i4)' ) 
     &    count, total, ( a(i), i = 1, k )

        if ( .not. more )  then
          go to 20
        end if

      go to 10

20    continue
     
      return
      end
      subroutine test1477

c*********************************************************************72
c
cc TEST1477 tests SUBCOMPNZ2_NEXT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer k
      parameter ( k = 3 )

      integer a(k)
      integer count
      integer i
      logical more
      integer n
      integer n_hi
      integer n_lo

      n_hi = 7
      n_lo = 5

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1477'
      write ( *, '(a)' ) '  SUBCOMPNZ2_NEXT generates subcompositions'
      write ( *, '(a)' ) '  using nonzero parts.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Seek all subcompositions of N'
      write ( *, '(a,i8,a)' ) '  using K = ', k, ' nonzero parts.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8,a,i8)' ) '  N ranges from ', n_lo, ' to ', n_hi
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     #     N'
      write ( *, '(a)' ) ' '

      more = .false.
      count = 0

10    continue

        call subcompnz2_next ( n_lo, n_hi, k, a, more )

        count = count + 1
        call i4vec_sum ( k, a, n )
        write ( *, '(2x,i4,2x,i4,2x,8i4)' ) 
     &    count, n, ( a(i), i = 1, k )

        if ( .not. more )  then
          go to 20
        end if

      go to 10

20    continue
     
      return
      end
      subroutine test148

c*********************************************************************72
c
cc TEST148 tests THUE_BINARY_NEXT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 100 )

      integer i
      integer j
      integer n
      integer thue(n_max)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST148'
      write ( *, '(a)' ) '  THUE_BINARY_NEXT returns the next '
      write ( *, '(a)' ) '  Thue binary sequence.'
      write ( *, '(a)' ) ' '

      n = 1
      thue(1) = 0
      write ( *, '(2x,i4,4x,80i1)' ) n, ( thue(i), i = 1, n )

      do j = 1, 6
        call thue_binary_next ( n, thue )
        write ( *, '(2x,i4,4x,80i1)' ) n, ( thue(i), i = 1, n )
      end do

      return
      end
      subroutine test149

c*********************************************************************72
c
cc TEST149 tests THUE_TERNARY_NEXT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 100 )

      integer i
      integer j
      integer n
      integer thue(n_max)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST149'
      write ( *, '(a)' ) '  THUE_TERNARY_NEXT returns the next '
      write ( *, '(a)' ) '  Thue ternary sequence.'
      write ( *, '(a)' ) ' '

      n = 1
      thue(1) = 1
      write ( *, '(2x,i4,4x,80i1)' ) n, ( thue(i), i = 1, n )

      do j = 1, 5
        call thue_ternary_next ( n, thue )
        write ( *, '(2x,i4,4x,80i1)' ) n, ( thue(i), i = 1, n )
      end do

      return
      end
      subroutine test150

c*********************************************************************72
c
cc TEST150 tests TUPLE_NEXT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 2 )

      integer i
      integer m1
      integer m2
      integer rank
      integer x(n)

      m1 = 2
      m2 = 4

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST150'
      write ( *, '(a)' ) 
     &  '  TUPLE_NEXT returns the next "tuple", that is,'
      write ( *, '(a)' ) 
     &  '  a vector of N integers, each between M1 and M2.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  M1 = ', m1
      write ( *, '(a,i8)' ) '  M2 = ', m2
      write ( *, '(a,i8)' ) '  N =  ', n
      write ( *, '(a)' ) ' '

      rank = 0

10    continue

        call tuple_next ( m1, m2, n, rank, x )

        if ( rank .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i4,2x,10i3)' ) rank, ( x(i), i = 1, n )

      go to 10

20    continue

      return
      end
      subroutine test151

c*********************************************************************72
c
cc TEST151 tests TUPLE_NEXT_FAST.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 2 )

      integer i
      integer m
      integer rank
      integer x(n)

      m = 3
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST151'
      write ( *, '(a)' ) 
     &  '  TUPLE_NEXT_FAST returns the next "tuple", that is,'
      write ( *, '(a)' ) 
     &  '  a vector of N integers, each between 1 and M.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  M = ', m
      write ( *, '(a,i8)' ) '  N = ', n
      write ( *, '(a)' ) ' '
c
c  Initialize.
c
      rank = -1
      call tuple_next_fast ( m, n, rank, x )

      do rank = 0, (m**n)-1

        call tuple_next_fast ( m, n, rank, x )

        write ( *, '(2x,i4,2x,10i3)' ) rank, ( x(i), i = 1, n )

      end do

      return
      end
      subroutine test152

c*********************************************************************72
c
cc TEST152 tests TUPLE_NEXT_GE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      integer i
      integer m
      integer rank
      integer x(n)

      m = 3

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST152'
      write ( *, '(a)' ) 
     &  '  TUPLE_NEXT_GE returns the next "tuple", that is,'
      write ( *, '(a)' ) 
     &  '  a vector of N integers, each between 1 and M,'
      write ( *, '(a)' ) '  with the constraint that the entries be'
      write ( *, '(a)' ) '  nondecreasing.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  M = ', m
      write ( *, '(a,i8)' ) '  N = ', n
      write ( *, '(a)' ) ' '

      rank = 0

10    continue

        call tuple_next_ge ( m, n, rank, x )

        if ( rank .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i4,2x,10i3)' ) rank, ( x(i), i = 1, n )

      go to 10

20    continue

      return
      end
      subroutine test153

c*********************************************************************72
c
cc TEST153 tests TUPLE_NEXT2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      integer i
      integer rank
      integer x(n)
      integer xmax(n)
      integer xmin(n)

      save xmax
      save xmin

      data xmax / 4, 3, 5 /
      data xmin / 2, 3, 8 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST153'
      write ( *, '(a)' ) 
     &  '  TUPLE_NEXT2 returns the next "tuple", that is,'
      write ( *, '(a)' ) '  a vector of N integers.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  N = ', n
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The minimum tuple is '
      write ( *, '(2x,5i8)' ) ( xmin(i), i = 1, n )
      write ( *, '(a)' ) '  The maximum tuple is '
      write ( *, '(2x,5i8)' ) ( xmax(i), i = 1, n )
      write ( *, '(a)' ) ' '

      rank = 0

10    continue

        call tuple_next2 ( n, xmin, xmax, rank, x )

        if ( rank .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i4,2x,10i3)' ) rank, ( x(i), i = 1, n )

      go to 10

20    continue

      return
      end
      subroutine test1531

c*********************************************************************72
c
cc TEST1531 tests UBVEC_ADD;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer bvec1(n)
      integer bvec2(n)
      integer bvec3(n)
      integer i
      integer i4_uniform
      integer j
      integer k
      integer seed
      integer test
      integer test_num

      test_num = 10
      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1531'
      write ( *, '(a)' ) '  UBVEC_ADD adds unsigned binary vectors'
      write ( *, '(a)' ) '  representing unsigned integers;'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        I        J        K = I + J'
      write ( *, '(a)' ) ' '

      do test = 1, test_num
        
        i = i4_uniform ( 0, 100, seed )
        j = i4_uniform ( 0, 100, seed )

        write ( *, '(a)' ) ' '

        write ( *, '(2x,i8,2x,i8)' ) i, j

        k = i + j

        write ( *, '(a20,2x,i8)' ) '  Directly:         ', k

        call ui4_to_ubvec ( i, n, bvec1 )
        call ui4_to_ubvec ( j, n, bvec2 )

        call ubvec_add ( n, bvec1, bvec2, bvec3 )
        call ubvec_to_ui4 ( n, bvec3, k )

        write ( *, '(a20,2x,i8)' ) '  UBVEC_ADD         ', k

      end do

      return
      end
      subroutine test0626

c*********************************************************************72
c
cc TEST0626 tests UI4_TO_UBVEC and UBVEC_TO_UI4;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n 
      parameter ( n = 10 )

      integer bvec(n)
      integer i
      integer i4
      integer j4

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0626'
      write ( *, '(a)' ) '  UI4_TO_UBVEC converts an unsigned integer'
      write ( *, '(a)' ) '    to an unsigned binary vector;'
      write ( *, '(a)' ) '  UBVEC_TO_UI4 converts an unsigned binary'
      write ( *, '(a)' ) '    vector to an unsigned integer;'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  I --> BVEC  -->  I'
      write ( *, '(a)' ) ' '
      do i4 = 0, 10
        call ui4_to_ubvec ( i4, n, bvec )
        call ubvec_to_ui4 ( n, bvec, j4 )
        write ( *, '(2x,i3,2x,10i1,2x,i3)' ) 
     &    i4, ( bvec(i), i = 1, n ), j4
      end do

      return
      end
      subroutine test154

c*********************************************************************72
c
cc TEST154 tests VEC_NEXT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      integer a(n)
      integer base
      integer i
      logical more

      base = 3
      more = .false.
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST154'
      write ( *, '(a)' ) '  VEC_NEXT generates all N-vectors in a'
      write ( *, '(a,i8)' ) 
     &  '  given base.  Here we use base ', base
      write ( *, '(a)' ) ' '
     
10    continue

        call vec_next ( n, base, a, more )

        write ( *, '(2x,3i4)' ) ( a(i), i = 1, n )

        if ( .not. more ) then
          go to 20
        end if

      go to 10

20    continue

      return
      end
      subroutine test155

c*********************************************************************72
c
cc TEST155 tests VEC_GRAY_NEXT, VEC_RANK and VEC_UNRANK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      integer a(n)
      integer base(n)
      integer change
      logical done
      integer i
      integer i4vec_product
      integer prod
      integer rank

      save base

      data base / 2, 2, 1, 4 /

      prod = i4vec_product ( n, base )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST155'
      write ( *, '(a)' ) 
     &  '  VEC_GRAY_NEXT generates product space elements.'
      write ( *, '(a)' ) '  VEC_RANK ranks them.'
      write ( *, '(a)' ) '  VEC_UNRANK unranks them.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The number of components is ', n
      write ( *, '(a,i8)' ) '  The number of elements is ', prod
      write ( *, '(a)' ) 
     &  '  Each component has its own number of degrees of'
      write ( *, '(a)' ) '  freedom, which, for this example, are:'
      write ( *, '(a)' ) ' '
      write ( *, '(a,6i4)' ) 
     &  '  Rank Change     ', ( base(i), i = 1, n )
      write ( *, '(a)' ) ' '
      rank = 0
      done = .true.
     
10    continue
     
        rank = rank + 1
     
        call vec_gray_next ( n, base, a, done, change )
     
        if ( done ) then
          go to 20
        end if

        write ( *, '(2x,i4,2x,i4,2x,4x,6i4)' ) 
     &    rank, change, ( a(i), i = 1, n )

      go to 10

20    continue
     
      do i = 1, n
        a(i) = base(i) / 2
      end do

      call vec_rank ( n, base, a, rank )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  VEC_RANK reports the element '
      write ( *, '(a)' ) ' '
      write ( *, '(4x,3x,6i4)' ) ( a(i), i = 1, n )
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  has rank ', rank

      rank = 7
      call vec_unrank ( n, base, rank, a )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 
     &  '  VEC_UNRANK reports the element of rank ', rank
      write ( *, '(a)' ) '  is:'
      write ( *, '(a)' ) ' '
      write ( *, '(4x,3x,6i4)' ) ( a(i), i = 1, n )
      write ( *, '(a)' ) ' '

      return
      end
      subroutine test156

c*********************************************************************72
c
cc TEST156 tests VEC_RANDOM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      integer a(n)
      integer base
      integer i
      integer j
      integer seed

      base = 3
      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST156'
      write ( *, '(a)' ) '  VEC_RANDOM generates a random N-vector'
      write ( *, '(a)' ) '  in a given base.'
      write ( *, '(a,i8)' ) '  Here, we use base ', base
      write ( *, '(a)' ) ' '

      do j = 1, 5
        call vec_random ( n, base, seed, a )
        write ( *, '(2x,3i4)' ) ( a(i), i = 1, n )
      end do
     
      return
      end
      subroutine test1565

c*********************************************************************72
c
cc TEST1565 tests VECTOR_CONSTRAINED_NEXT
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      integer constraint
      integer i
      integer i4vec_product
      logical more
      integer rank
      integer x(n)
      integer x_max(n)
      integer x_min(n)
      integer x_prod

      save x_max
      save x_min

      data x_max / 4, 5, 3 /
      data x_min / 2, 2, 1 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1565'
      write ( *, '(a)' ) '  VECTOR_CONSTRAINED_NEXT:'
      write ( *, '(a)' ) '  Consider vectors:'
      write ( *, '(a)' ) '    X_MIN(1:N) <= X(1:N) <= X_MAX(1:N),'
      write ( *, '(a)' ) '  Set'
      write ( *, '(a)' ) '    P = Product X_MAX(1:N)'
      write ( *, '(a)' ) '  Accept only vectors for which:'
      write ( *, '(a)' ) '    sum ( (X(1:N)-1) * P / X_MAX(1:N) ) <= P'

      more = .false.

      write ( *, '(a)' ) ' '
      write ( *, '(a,3i4)' ) '  X_MIN:', ( x_min(i), i = 1, n )
      write ( *, '(a,3i4)' ) '  X_MAX:', ( x_max(i), i = 1, n )

      rank = 0
      x_prod = i4vec_product ( n, x_max )

      write ( *, '(a)' ) ' '

      write ( *, '(a,i12)' ) 
     &  '  Maximum allowed CONSTRAINT = P = ', x_prod
      write ( *, '(a)' ) ' '

10    continue

        call vector_constrained_next ( n, x_min, x_max, x, 
     &    constraint, more )

        if ( .not. more ) then
          go to 20
        end if

        rank = rank + 1
        write ( *, '(2x,i8,2x,i12,2x,i8,2x,i8,2x,i8)' ) 
     &    rank, constraint, ( x(i), i = 1, n )

      go to 10

20    continue

      return
      end
      subroutine test1566

c*********************************************************************72
c
cc TEST1566 tests VECTOR_CONSTRAINED_NEXT2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 3 )

      integer constraint
      integer i
      integer i4vec_product
      logical more
      integer n
      integer rank
      integer x(n_max)
      integer x_max(n_max)
      integer x_min(n_max)
      integer x_prod

      save x_max
      save x_min

      data x_max / 5, 6, 4 /
      data x_min / 1, 1, 1 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1566'
      write ( *, '(a)' ) '  VECTOR_CONSTRAINED_NEXT2:'
      write ( *, '(a)' ) '  Consider vectors:'
      write ( *, '(a)' ) '    X_MIN(1:N) <= X(1:N) <= X_MAX(1:N),'
      write ( *, '(a)' ) '  Set'
      write ( *, '(a)' ) '    P = Product X_MAX(1:N)'
      write ( *, '(a)' ) '  Accept only vectors for which:'
      write ( *, '(a)' ) '    sum ( X(1:N) * P / X_MAX(1:N) ) <= P'

      do n = 2, n_max

        more = .false.

        write ( *, '(a)' ) ' '
        write ( *, '(a,3i4)' ) '  X_MIN:', ( x_min(i), i = 1, n )
        write ( *, '(a,3i4)' ) '  X_MAX:', ( x_max(i), i = 1, n )

        rank = 0
        x_prod = i4vec_product ( n, x_max )

        write ( *, '(a)' ) ' '
        write ( *, '(a,i12)' ) 
     &    '  Maximum allowed CONSTRAINT = P = ', x_prod
        write ( *, '(a)' ) ' '

10      continue

          call vector_constrained_next2 ( n, x_min, x_max, x, 
     &      constraint, more )

          if ( .not. more ) then
            go to 20
          end if

          rank = rank + 1
          write ( *, '(2x,i8,2x,i12,2x,i8,2x,i8,2x,i8)' ) 
     &      rank, constraint, ( x(i), i = 1, n )

        go to 10

20      continue

      end do

      return
      end
      subroutine test1567

c*********************************************************************72
c
cc TEST1567 tests VECTOR_CONSTRAINED_NEXT3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 3 )

      double precision constraint
      integer i
      logical more
      integer n
      integer rank
      integer x(n_max)
      integer x_max(n_max)
      integer x_min(n_max)
      integer x_prod

      save x_max
      save x_min

      data x_max / 5, 6, 4 /
      data x_min / 1, 1, 1 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1567'
      write ( *, '(a)' ) '  VECTOR_CONSTRAINED_NEXT3:'
      write ( *, '(a)' ) '  Consider vectors:'
      write ( *, '(a)' ) '    X_MIN(1:N) <= X(1:N) <= X_MAX(1:N),'
      write ( *, '(a)' ) '  Set'
      write ( *, '(a)' ) '    CONSTRAINT = sum ( X(1:N) / X_MAX(1:N) )'
      write ( *, '(a)' ) '  Accept only vectors for which:'
      write ( *, '(a)' ) '    CONSTRAINT <= 1'

      do n = 2, n_max

        more = .false.

        write ( *, '(a)' ) ' '
        write ( *, '(a,3i4)' ) '  X_MIN:', ( x_min(i), i = 1, n )
        write ( *, '(a,3i4)' ) '  X_MAX:', ( x_max(i), i = 1, n )
        write ( *, '(a)' ) ' '

        rank = 0

10      continue

          call vector_constrained_next3 ( n, x_min, x_max, x, 
     &      constraint, more )

          if ( .not. more ) then
            go to 20
          end if

          rank = rank + 1
          write ( *, '(2x,i8,2x,g14.6,2x,i8,2x,i8,2x,i8)' ) 
     &      rank, constraint, ( x(i), i = 1, n )

        go to 10

20      continue

      end do

      return
      end
      subroutine test1568

c*********************************************************************72
c
cc TEST1568 tests VECTOR_CONSTRAINED_NEXT4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 3 )

      double precision alpha(n_max)
      integer i
      logical more
      integer n
      double precision q
      integer rank
      double precision total
      integer x(n_max)
      integer x_max(n_max)
      integer x_min(n_max)

      save alpha
      save x_max
      save x_min

      data alpha / 4.0D+00, 3.0D+00, 5.0D+00 /
      data x_max / 2, 6, 4 /
      data x_min / 1, 0, 1 /

      q = 20.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1568'
      write ( *, '(a)' ) '  VECTOR_CONSTRAINED_NEXT4:'
      write ( *, '(a)' ) '  Consider vectors:'
      write ( *, '(a)' ) '    X_MIN(1:N) <= X(1:N) <= X_MAX(1:N),'
      write ( *, '(a)' ) '  Set'
      write ( *, '(a)' ) '    TOTAL = sum ( ALPHA(1:N) * X(1:N) )'
      write ( *, '(a)' ) '  Accept only vectors for which:'
      write ( *, '(a)' ) '    TOTAL <= Q'

      do n = 2, n_max

        more = .false.

        write ( *, '(a)' ) ' '
        write ( *, '(a,3g14.6)' ) 
     &    '  ALPHA:', ( alpha(i), i = 1, n )
        write ( *, '(a, g14.6)' ) '  Q:    ', q
        write ( *, '(a,3i4)'    ) 
     &    '  X_MIN:', ( x_min(i), i = 1, n )
        write ( *, '(a,3i4)'    ) 
     &    '  X_MAX:', ( x_max(i), i = 1, n )
        write ( *, '(a)' ) ' '

        rank = 0

10      continue

          call vector_constrained_next4 ( n, alpha, x_min, x_max, 
     &      x, q, more )

          if ( .not. more ) then
            go to 20
          end if

          total = 0.0D+00
          do i = 1, n
            total = total + alpha(i) * dble ( x(i) )
          end do

          rank = rank + 1
          write ( *, '(2x,i8,2x,g14.6,2x,i8,2x,i8,2x,i8)' ) 
     &      rank, total, ( x(i), i = 1, n )

        go to 10

20      continue

      end do

      return
      end
      subroutine test1569

c*********************************************************************72
c
cc TEST1569 tests VECTOR_CONSTRAINED_NEXT5
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      integer i
      logical more
      integer rank
      integer sum_max
      integer sum_min
      integer x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1569'
      write ( *, '(a)' ) '  VECTOR_CONSTRAINED_NEXT5:'
      write ( *, '(a)' ) '  Generate integer vectors X such that:'
      write ( *, '(a)' ) '    SUM_MIN <= sum ( X(1:N) ) <= SUM_MAX,'
      write ( *, '(a)' ) '  We require every X(I) to be at least 1.'

      more = .false.
      sum_min = 5
      sum_max = 7

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  N =       ', n
      write ( *, '(a,i8)' ) '  SUM_MIN = ', sum_min
      write ( *, '(a,i8)' ) '  SUM_MAX = ', sum_max
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '         #        X(1)      X(2)      X(3)'
      write ( *, '(a)' ) ' '

      rank = 0

10    continue

        call vector_constrained_next5 ( n, x, sum_min, sum_max, 
     &    more )

        if ( .not. more ) then
          go to 20
        end if

        rank = rank + 1
        write ( *, '(2x,i8,2x,i8,2x,i8,2x,i8)' ) 
     &    rank, ( x(i), i = 1, n )

      go to 10

20    continue

      return
      end
      subroutine test15695

c*********************************************************************72
c
cc TEST15695 tests VECTOR_CONSTRAINED_NEXT4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 February 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 3 )

      double precision alpha(n_max)
      integer i
      integer indx
      logical more
      integer n
      double precision q_max
      double precision q_min
      double precision total
      integer x(n_max)
      integer x_max(n_max)
      integer x_min(n_max)

      save alpha
      save x_max
      save x_min

      data alpha / 4.0D+00, 3.0D+00, 5.0D+00 /
      data x_max / 2, 6, 4 /
      data x_min / 1, 0, 1 /

      q_min = 16.0D+00
      q_max = 20.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST15695'
      write ( *, '(a)' ) '  VECTOR_CONSTRAINED_NEXT6:'
      write ( *, '(a)' ) '  Consider vectors:'
      write ( *, '(a)' ) '    X_MIN(1:N) <= X(1:N) <= X_MAX(1:N),'
      write ( *, '(a)' ) '  Set'
      write ( *, '(a)' ) '    TOTAL = sum ( ALPHA(1:N) * X(1:N) )'
      write ( *, '(a)' ) '  Accept only vectors for which:'
      write ( *, '(a)' ) '    Q_MIN <= TOTAL <= Q_MAX'

      do n = 2, n_max

        more = .false.

        write ( *, '(a)' ) ' '
        write ( *, '(a,3g14.6)' ) '  ALPHA:', 
     &  ( alpha(i), i = 1, n )
        write ( *, '(a, g14.6)' ) '  Q_MIN:', q_min
        write ( *, '(a, g14.6)' ) '  Q_MAX:', q_max
        write ( *, '(a,3i4)'    ) '  X_MIN:', 
     &  ( x_min(i), i = 1, n )
        write ( *, '(a,3i4)'    ) '  X_MAX:', 
     &  ( x_max(i), i = 1, n )
        write ( *, '(a)' ) ' '

        indx = 0

10      continue

          call vector_constrained_next6 ( n, alpha, x_min, x_max, 
     &      x, q_min, q_max, more )

          if ( .not. more ) then
            go to 20
          end if

          total = 0.0D+00
          do i = 1, n
            total = total + alpha(i) * dble ( x(i) )
          end do

          indx = indx + 1
          write ( *, '(2x,i8,2x,g14.6,2x,i8,2x,i8,2x,i8)' ) 
     &    indx, total, ( x(i), i = 1, n )

        go to 10

20      continue

      end do

      return
      end
      subroutine test15696 ( )

c*********************************************************************72
c
cc TEST15696 tests VECTOR_CONSTRAINED_NEXT7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 3 )

      double precision alpha(n_max)
      integer i
      integer j
      logical more
      integer n
      double precision q_max
      double precision q_min
      double precision total
      integer x(n_max)
      integer x_max(n_max)

      save alpha
      save x_max

      data alpha / 4.0D+00, 3.0D+00, 5.0D+00 /
      data x_max / 2, 6, 4 /

      q_min = 16.0D+00
      q_max = 20.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST15696'
      write ( *, '(a)' ) '  VECTOR_CONSTRAINED_NEXT7:'
      write ( *, '(a)' ) '  Consider vectors:'
      write ( *, '(a)' ) '    0 <= X(1:N) <= X_MAX(1:N),'
      write ( *, '(a)' ) '  Set'
      write ( *, '(a)' ) '    TOTAL = sum ( ALPHA(1:N) * X(1:N) )'
      write ( *, '(a)' ) '  Accept only vectors for which:'
      write ( *, '(a)' ) '    Q_MIN <= TOTAL <= Q_MAX'

      do n = 2, n_max

        more = .false.

        write ( *, '(a)' ) ' '
        write ( *, '(a,3g14.6)' ) '  ALPHA:', ( alpha(j), j = 1, n )
        write ( *, '(a, g14.6)' ) '  Q_MIN:', q_min
        write ( *, '(a, g14.6)' ) '  Q_MAX:', q_max
        write ( *, '(a,3i4)'    ) '  X_MAX:', ( x_max(j), j = 1, n )
        write ( *, '(a)' ) ' '

        i = 0

10      continue

          call vector_constrained_next7 ( n, alpha, x_max, x, q_min, 
     &      q_max, more )

          if ( .not. more ) then
            go to 20
          end if

          total = 0.0D+00
          do j = 1, n
            total = total + alpha(j) * dble ( x(j) )
          end do
          i = i + 1
          write ( *, '(2x,i8,2x,g14.6,2x,i8,2x,i8,2x,i8)' )
     &      i, total, ( x(j), j = 1, n )

        go to 10

20      continue

      end do

      return
      end
      subroutine test15698 ( )

!*********************************************************************72
!
!! TEST15698 tests VECTOR_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 July 2010
!
!  Author:
!
!    John Burkardt
!
      implicit none

      integer n_max 
      parameter ( n_max = 3 )

      integer i
      integer j
      logical more
      integer n
      integer x(n_max)
      integer x_max(n_max)
      integer x_min(n_max)

      save x_max
      save x_min

      data x_max / 2, 6, 4 /
      data x_min / 1, 4, 3 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST15698'
      write ( *, '(a)' ) '  VECTOR_NEXT:'
      write ( *, '(a)' ) '  Generate all vectors X such that:'
      write ( *, '(a)' ) '    X_MIN(1:N) <= X(1:N) <= X_MAX(1:N),'

      do n = 2, n_max

        more = .false.

        write ( *, '(a)' ) ' '
        write ( *, '(2x,a4,4x,i8,2x,i8,2x,i8)' ) 
     &    'XMIN', ( x_min(j), j = 1, n )

        i = 0

10      continue

          call vector_next ( n, x_min, x_max, x, more )

          if ( .not. more ) then
            go to 20
          end if

          i = i + 1
          write ( *, '(2x,i4,4x,i8,2x,i8,2x,i8)' ) 
     &      i, ( x(j), j = 1, n )

        go to 10

20      continue

        write ( *, '(2x,a4,4x,i8,2x,i8,2x,i8)' ) 
     &    'XMAX', ( x_max(j), j = 1, n )

      end do

      return
      end
      subroutine test157

c*********************************************************************72
c
cc TEST157 tests YTB_ENUM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer i
      integer pi

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST157'
      write ( *, '(a)' ) '  YTB_ENUM counts Young table.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  N    YTB(N)'
      write ( *, '(a)' ) ' '

      do i = 0, n
        call ytb_enum ( i, pi )
        write ( *, '(2x,2i10)' ) i, pi
      end do

      return
      end
      subroutine test158

c*********************************************************************72
c
cc TEST158 tests YTB_NEXT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 6 )

      integer a(n)
      integer i
      integer lambda(n)
      logical more

      save lambda

      data lambda / 3, 2, 1, 0, 0, 0 /

      do i = 1, n
        a(i) = 0
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST158'
      write ( *, '(a)' ) '  YTB_NEXT generates Young tables.'
      write ( *, '(a)' ) ' '
      more = .false.
     
10    continue
     
        call ytb_next ( n, lambda, a, more )
     
        call ytb_print ( n, a, ' ' )
     
        if ( .not. more ) then
          go to 20
        end if

      go to 10

20    continue

      return
      end
      subroutine test159

c*********************************************************************72
c
cc TEST159 tests YTB_RANDOM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 6 )

      integer a(n)
      integer i
      integer lambda(n)
      integer seed

      save lambda

      data lambda / 3, 2, 1, 0, 0, 0 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST159'
      write ( *, '(a)' ) '  YTB_RANDOM generates a random Young table'

      seed = 123456789

      do i = 1, 5
     
        call ytb_random ( n, lambda, seed, a )

        call ytb_print ( n, a, ' ' )
     
      end do
     
      return
      end
