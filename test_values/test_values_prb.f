      program main

c*********************************************************************72
c
cc MAIN is the main program for TEST_VALUES_PRB.
c
c  Discussion:
c
c    TEST_VALUES_PRB calls the TEST_VALUE routines.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 March 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_VALUES_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TEST_VALUES library.'

      call test001 ( )
      call test002 ( )
      call test003 ( )
      call test0035 ( )
      call test004 ( )
      call test005 ( )
      call test006 ( )
      call test007 ( )
      call test008 ( )
      call test009 ( )

      call test010 ( )
      call test011 ( )
      call test012 ( )
      call test0127 ( )
      call test0128 ( )
      call test013 ( )
      call test0134 ( )
      call test0135 ( )
      call test014 ( )
      call test015 ( )
      call test016 ( )
      call test017 ( )
      call test018 ( )
      call test0185 ( )
      call test019 ( )
      call test0195 ( )

      call test020 ( )
      call test0205 ( )
      call test021 ( )
      call test022 ( )
      call test023 ( )
      call test024 ( )
      call test025 ( )
      call test026 ( )
      call test0265 ( )
      call test027 ( )
      call test028 ( )
      call test029 ( )

      call test030 ( )
      call test0305 ( )
      call test031 ( )
      call test032 ( )
      call test033 ( )
      call test034 ( )
      call test035 ( )
      call test036 ( )
      call test0365 ( )
      call test037 ( )
      call test038 ( )
      call test039 ( )
      call test0395 ( )

      call test040 ( )
      call test041 ( )
      call test042 ( )
      call test0425 ( )
      call test043 ( )
      call test044 ( )
      call test045 ( )
      call test046 ( )
      call test047 ( )
      call test048 ( )
      call test049 ( )

      call test050 ( )
      call test051 ( )
      call test05125 ( )
      call test0515 ( )
      call test052 ( )
      call test053 ( )
      call test054 ( )
      call test055 ( )
      call test056 ( )
      call test057 ( )
      call test058 ( )
      call test059 ( )

      call test060 ( )
      call test061 ( )
      call test062 ( )
      call test063 ( )
      call test064 ( )
      call test065 ( )
      call test066 ( )
      call test0665 ( )
      call test067 ( )
      call test068 ( )
      call test0685 ( )
      call test069 ( )

      call test070 ( )
      call test071 ( )
      call test072 ( )
      call test073 ( )
      call test074 ( )
      call test075 ( )
      call test0755 ( )
      call test076 ( )
      call test077 ( )
      call test078 ( )
      call test079 ( )

      call test080 ( )
      call test081 ( )
      call test082 ( )
      call test083 ( )
      call test0835 ( )
      call test084 ( )
      call test0843 ( )
      call test0845 ( )
      call test085 ( )
      call test0855 ( )
      call test086 ( )
      call test087 ( )
      call test088 ( )
      call test089 ( )

      call test090 ( )
      call test091 ( )
      call test092 ( )
      call test093 ( )
      call test094 ( )
      call test0945 ( )
      call test095 ( )
      call test096 ( )
      call test097 ( )
      call test0972 ( )
      call test0973 ( )
      call test0974 ( )
      call test0975 ( )
      call test098 ( )
      call test099 ( )
      call test0995 ( )

      call test100 ( )
      call test101 ( )
      call test1015 ( )
      call test1016 ( )
      call test102 ( )
      call test103 ( )
      call test104 ( )
      call test1037 ( )
      call test105 ( )
      call test106 ( )
      call test107 ( )
      call test108 ( )
      call test10825 ( )
      call test10850 ( )
      call test109 ( )

      call test110 ( )
      call test1105 ( )
      call test111 ( )
      call test112 ( )
      call test113 ( )
      call test1135 ( )
      call test114 ( )
      call test115 ( )
      call test116 ( )
      call test117 ( )
      call test118 ( )
      call test1185 ( )
      call test119 ( )

      call test120 ( )
      call test121 ( )
      call test122 ( )
      call test123 ( )
      call test124 ( )
      call test125 ( )
      call test1255 ( )
      call test126 ( )
      call test127 ( )
      call test128 ( )
      call test1285 ( )
      call test129 ( )

      call test131 ( )
      call test132 ( )
      call test1325 ( )
      call test130 ( )
      call test133 ( )
      call test134 ( )
      call test135 ( )
      call test136 ( )
      call test137 ( )
      call test138 ( )
      call test139 ( )

      call test140 ( )
      call test141 ( )
      call test1415 ( )
      call test142 ( )
      call test143 ( )
      call test144 ( )
      call test145 ( )
      call test146 ( )
      call test1465 ( )
      call test147 ( )
      call test148 ( )
      call test149 ( )

      call test150 ( )
      call test151 ( )
      call test152 ( )
      call test153 ( )
      call test154 ( )
      call test1545 ( )
      call test155 ( )
      call test156 ( )
      call test157 ( )
      call test1575 ( )
      call test158 ( )
      call test159 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_VALUES_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test001 ( )

c*********************************************************************72
c
cc TEST001 demonstrates the use of ABRAM0_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 February 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST001:'
      write ( *, '(a)' ) '  ABRAM0_VALUES returns values of '
      write ( *, '(a)' ) '  the Abramowitz function of order 0'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X           Abram0'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call abram0_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test002 ( )

c*********************************************************************72
c
cc TEST002 demonstrates the use of ABRAM1_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 February 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST002:'
      write ( *, '(a)' ) '  ABRAM1_VALUES returns values of '
      write ( *, '(a)' ) '  the Abramowitz function of order 1'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X           Abram1'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call abram1_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test003 ( )

c*********************************************************************72
c
cc TEST003 demonstrates the use of ABRAM2_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 February 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST003:'
      write ( *, '(a)' ) '  ABRAM2_VALUES returns values of '
      write ( *, '(a)' ) '  the Abramowitz function of order 2'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X           Abram2'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call abram2_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test0035 ( )

c*********************************************************************72
c
cc TEST0035 demonstrates the use of AGM_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision b
      double precision fx
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0035:'
      write ( *, '(a)' ) '  AGM_VALUES returns'
      write ( *, '(a)' ) '  values of the arithmetic geometric '
      write ( *, '(a)' ) '  mean function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '          A               B            AGM(A,B)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call agm_values ( n_data, a, b, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,f14.6,2x,g24.16)' ) a, b, fx

      go to 10

20    continue

      return
      end
      subroutine test004 ( )

c*********************************************************************72
c
cc TEST004 demonstrates the use of AIRY_AI_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 February 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST004:'
      write ( *, '(a)' ) '  AIRY_AI_VALUES returns values of '
      write ( *, '(a)' ) '  the Airy function Ai(X)'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X           Ai'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call airy_ai_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test005 ( )

c*********************************************************************72
c
cc TEST005 demonstrates the use of AIRY_AI_INT_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST005:'
      write ( *, '(a)' ) '  AIRY_AI_INT_VALUES returns values of '
      write ( *, '(a)' ) '  the integral of the Airy function Ai(X)'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X           Ai_Int'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call airy_ai_int_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test006 ( )

c*********************************************************************72
c
cc TEST006 demonstrates the use of AIRY_AI_PRIME_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST006:'
      write ( *, '(a)' ) '  AIRY_AI_PRIME_VALUES returns values of '
      write ( *, '(a)' ) '  the derivative of the Airy function Ai(X)'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X           AiP'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call airy_ai_prime_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test007 ( )

c*********************************************************************72
c
cc TEST007 demonstrates the use of AIRY_BI_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST007:'
      write ( *, '(a)' ) '  AIRY_BI_VALUES returns values of '
      write ( *, '(a)' ) '  the Airy function Bi(X)'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X           Bi'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call airy_bi_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test008 ( )

c*********************************************************************72
c
cc TEST008 demonstrates the use of AIRY_BI_INT_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST008:'
      write ( *, '(a)' ) '  AIRY_BI_INT_VALUES returns values of '
      write ( *, '(a)' ) '  the integral of the Airy function Bi(X)'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X           Bi_Int'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call airy_bi_int_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test009 ( )

c*********************************************************************72
c
cc TEST009 demonstrates the use of AIRY_BI_PRIME_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST009:'
      write ( *, '(a)' ) '  AIRY_BI_PRIME_VALUES returns values of '
      write ( *, '(a)' ) '  the derivative of the Airy function Bi(X)'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X           BiP'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call airy_bi_prime_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test010 ( )

c*********************************************************************72
c
cc TEST010 demonstrates the use of AIRY_GI_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST010:'
      write ( *, '(a)' ) '  AIRY_GI_VALUES returns values of '
      write ( *, '(a)' ) '  the modified Airy function Gi(X)'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X           Gi'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call airy_gi_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test011 ( )

c*********************************************************************72
c
cc TEST011 demonstrates the use of AIRY_HI_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST011:'
      write ( *, '(a)' ) '  AIRY_HI_VALUES returns values of '
      write ( *, '(a)' ) '  the modified Airy function Hi(X)'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X           Hi'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call airy_hi_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test012 ( )

c*********************************************************************72
c
cc TEST012 demonstrates the use of ARCTAN_INT_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST012:'
      write ( *, '(a)' ) '  ARCTAN_INT_VALUES returns values of '
      write ( *, '(a)' ) '  the arctangent integral.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X           F(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call arctan_int_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test0127 ( )

c*********************************************************************72
c
cc TEST0127 demonstrates the use of BEI0_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0127:'
      write ( *, '(a)' ) '  BEI0_VALUES returns values of '
      write ( *, '(a)' ) '  the Kelvin BEI function of order 0.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X                     FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bei0_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f24.16,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test0128 ( )

c*********************************************************************72
c
cc TEST0128 demonstrates the use of BEI1_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0128'
      write ( *, '(a)' ) '  BEI1 VALUES returns values of '
      write ( *, '(a)' ) '  the Kelvin BEI function of order 1.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X                     FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bei1_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f24.16,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test013 ( )

c*********************************************************************72
c
cc TEST013 demonstrates the use of BELL_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer c
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST013:'
      write ( *, '(a)' ) '  BELL_VALUES returns values of '
      write ( *, '(a)' ) '  the Bell numbers.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         N        BELL(N)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bell_values ( n_data, n, c )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,i12)' ) n, c

      go to 10

20    continue

      return
      end
      subroutine test0134 ( )

c*********************************************************************72
c
cc TEST0134 demonstrates the use of BER0_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0134:'
      write ( *, '(a)' ) '  BER0_VALUES returns values of '
      write ( *, '(a)' ) '  the Kelvin BER function of order 0.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X                     FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call ber0_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f24.16,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test0135 ( )

c*********************************************************************72
c
cc TEST0135 demonstrates the use of BER1_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0135:'
      write ( *, '(a)' ) '  BER1_VALUES returns values of '
      write ( *, '(a)' ) '  the Kelvin BER function of order 1.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X                     FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bei1_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f24.16,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test014 ( )

c*********************************************************************72
c
cc TEST014 demonstrates the use of BERNOULLI_NUMBER_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision c
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST014:'
      write ( *, '(a)' ) '  BERNOULLI_NUMBER_VALUES returns values'
      write ( *, '(a)' ) '  of the Bernoulli numbers.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         N        BERNOULLI_NUMBER(N)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bernoulli_number_values ( n_data, n, c )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,6x,g24.16)' ) n, c

      go to 10

20    continue

      return
      end
      subroutine test015 ( )

c*********************************************************************72
c
cc TEST015 demonstrates the use of BERNOULLI_POLY_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST015:'
      write ( *, '(a)' ) '  BERNOULLI_POLY_VALUES returns values of '
      write ( *, '(a)' ) '  the Bernoulli polynomials.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '         N          X            BERNOULLI_POLY(N,X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bernoulli_poly_values ( n_data, n, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,f14.6,2x,g24.16)' ) n, x, fx

      go to 10

20    continue

      return
      end
      subroutine test016 ( )

c*********************************************************************72
c
cc TEST016 demonstrates the use of BERNSTEIN_POLY_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision b
      integer k
      integer n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST016:'
      write ( *, '(a)' ) '  BERNSTEIN_POLY_VALUES returns values of '
      write ( *, '(a)' ) '  the Bernstein polynomials.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '         N         K          X           BERNSTEIN(N,K)(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bernstein_poly_values ( n_data, n, k, x, b )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,i8,2x,f14.6,2x,g24.16)' ) n, k, x, b

      go to 10

20    continue

      return
      end
      subroutine test017 ( )

c*********************************************************************72
c
cc TEST017 demonstrates the use of BESSEL_I0_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST017:'
      write ( *, '(a)' ) '  BESSEL_I0_VALUES returns values of '
      write ( *, '(a)' ) '  the Bessel I0 function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X                     FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_i0_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f24.16,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test018 ( )

c*********************************************************************72
c
cc TEST018 demonstrates the use of BESSEL_I0_INT_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST018:'
      write ( *, '(a)' ) '  BESSEL_I0_INT_VALUES returns values of '
      write ( *, '(a)' ) '  the integral of the Bessel I0 function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X                     FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_i0_int_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f24.16,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test0185 ( )

c*********************************************************************72
c
cc TEST0185 demonstrates the use of BESSEL_I0_SPHERICAL_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0185:'
      write ( *, '(a)' ) '  BESSEL_I0_SPHERICAL_VALUES returns values'
      write ( *, '(a)' ) '  of the spherical Bessel i0 function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X                     FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_i0_spherical_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f24.16,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test019 ( )

c*********************************************************************72
c
cc TEST019 demonstrates the use of BESSEL_I1_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST019:'
      write ( *, '(a)' ) '  BESSEL_I1_VALUES returns values of '
      write ( *, '(a)' ) '  the Bessel I1 function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X                     FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_i1_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f24.16,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test0195 ( )

c*********************************************************************72
c
cc TEST0195 demonstrates the use of BESSEL_I1_SPHERICAL_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0195:'
      write ( *, '(a)' ) '  BESSEL_I1_SPHERICAL_VALUES returns values'
      write ( *, '(a)' ) '  of the spherical Bessel i1 function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X                     FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_i1_spherical_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f24.16,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test020 ( )

c*********************************************************************72
c
cc TEST020 demonstrates the use of BESSEL_IN_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST020:'
      write ( *, '(a)' ) '  BESSEL_IN_VALUES returns values of '
      write ( *, '(a)' ) '  the Bessel In function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N          X                     FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_in_values ( n_data, n, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i6,2x,f24.16,2x,g24.16)' ) n, x, fx

      go to 10

20    continue

      return
      end
      subroutine test0205 ( )

c*********************************************************************72
c
cc TEST0205 demonstrates the use of BESSEL_IX_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      double precision n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0205:'
      write ( *, '(a)' ) '  BESSEL_IX_VALUES returns values of '
      write ( *, '(a)' ) '  the Bessel In function with NONINTEGER'
      write ( *, '(a)' ) '  order.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '          N             X                     FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_ix_values ( n_data, n, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f12.6,2x,f24.16,2x,g24.16)' ) n, x, fx

      go to 10

20    continue

      return
      end
      subroutine test021 ( )

c*********************************************************************72
c
cc TEST021 demonstrates the use of BESSEL_J0_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST021:'
      write ( *, '(a)' ) '  BESSEL_J0_VALUES returns values of '
      write ( *, '(a)' ) '  the Bessel J0 function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X                     FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_j0_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f24.16,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test022 ( )

c*********************************************************************72
c
cc TEST022 demonstrates the use of BESSEL_J0_INT_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST022:'
      write ( *, '(a)' ) '  BESSEL_J0_INT_VALUES returns values of '
      write ( *, '(a)' ) '  the integral of the Bessel J0 function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X                     FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_j0_int_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f24.16,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test023 ( )

c*********************************************************************72
c
cc TEST023 demonstrates the use of BESSEL_J0_SPHERICAL_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST023:'
      write ( *, '(a)' ) '  BESSEL_J0_SPHERICAL_VALUES returns values'
      write ( *, '(a)' ) '  of the spherical Bessel j0 function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X                     FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_j0_spherical_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f24.16,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test024 ( )

c*********************************************************************72
c
cc TEST024 demonstrates the use of BESSEL_J1_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST024:'
      write ( *, '(a)' ) '  BESSEL_J1_VALUES returns values of '
      write ( *, '(a)' ) '  the Bessel J10 function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X                     FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_j1_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f24.16,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test025 ( )

c*********************************************************************72
c
cc TEST025 demonstrates the use of BESSEL_J1_SPHERICAL_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST025:'
      write ( *, '(a)' ) '  BESSEL_J1_SPHERICAL_VALUES returns values'
      write ( *, '(a)' ) '  of the spherical Bessel j1 function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X                     FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_j1_spherical_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f24.16,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test026 ( )

c*********************************************************************72
c
cc TEST026 demonstrates the use of BESSEL_JN_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST026:'
      write ( *, '(a)' ) '  BESSEL_JN_VALUES returns values of '
      write ( *, '(a)' ) '  the Bessel Jn function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N          X                     FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_jn_values ( n_data, n, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i6,2x,f24.16,2x,g24.16)' ) n, x, fx

      go to 10

20    continue

      return
      end
      subroutine test0265 ( )

c*********************************************************************72
c
cc TEST0265 demonstrates the use of BESSEL_JX_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      double precision n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0265:'
      write ( *, '(a)' ) '  BESSEL_JX_VALUES returns values of '
      write ( *, '(a)' ) '  the Bessel Jn function with NONINTEGER'
      write ( *, '(a)' ) '  order.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '          N             X                     FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_jx_values ( n_data, n, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f12.6,2x,f24.16,2x,g24.16)' ) n, x, fx

      go to 10

20    continue

      return
      end
      subroutine test027 ( )

c*********************************************************************72
c
cc TEST027 demonstrates the use of BESSEL_K0_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST027:'
      write ( *, '(a)' ) '  BESSEL_K0_VALUES returns values of '
      write ( *, '(a)' ) '  the Bessel K0 function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X                     FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_k0_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f24.16,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test028 ( )

c*********************************************************************72
c
cc TEST028 demonstrates the use of BESSEL_K0_INT_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST028:'
      write ( *, '(a)' ) '  BESSEL_K0_INT_VALUES returns values of '
      write ( *, '(a)' ) '  the integral of the Bessel K0 function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X                     FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_k0_int_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f24.16,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test029 ( )

c*********************************************************************72
c
cc TEST029 demonstrates the use of BESSEL_K1_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST029:'
      write ( *, '(a)' ) '  BESSEL_K1_VALUES returns values of '
      write ( *, '(a)' ) '  the Bessel K1 function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X                     FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_k1_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f24.16,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test030 ( )

c*********************************************************************72
c
cc TEST030 demonstrates the use of BESSEL_KN_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST030:'
      write ( *, '(a)' ) '  BESSEL_KN_VALUES returns values of '
      write ( *, '(a)' ) '  the Bessel Kn function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N          X                     FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_kn_values ( n_data, n, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i6,2x,f24.16,2x,g24.16)' ) n, x, fx

      go to 10

20    continue

      return
      end
      subroutine test0305 ( )

c*********************************************************************72
c
cc TEST0305 demonstrates the use of BESSEL_KX_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      double precision n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0305:'
      write ( *, '(a)' ) '  BESSEL_KX_VALUES returns values of '
      write ( *, '(a)' ) '  the Bessel Kn function with NONINTEGER'
      write ( *, '(a)' ) '  order.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '          N             X                     FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_kx_values ( n_data, n, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f12.6,2x,f24.16,2x,g24.16)' ) n, x, fx

      go to 10

20    continue

      return
      end
      subroutine test031 ( )

c*********************************************************************72
c
cc TEST031 demonstrates the use of BESSEL_Y0_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST031:'
      write ( *, '(a)' ) '  BESSEL_Y0_VALUES returns values of '
      write ( *, '(a)' ) '  the Bessel Y0 function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X                     FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_y0_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f24.16,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test032 ( )

c*********************************************************************72
c
cc TEST032 demonstrates the use of BESSEL_Y0_INT_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST032:'
      write ( *, '(a)' ) '  BESSEL_Y0_INT_VALUES returns values of '
      write ( *, '(a)' ) '  the integral of the Bessel Y0 function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X                     FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_y0_int_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f24.16,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test033 ( )

c*********************************************************************72
c
cc TEST033 demonstrates the use of BESSEL_Y0_SPHERICAL_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST033:'
      write ( *, '(a)' ) '  BESSEL_Y0_SPHERICAL_VALUES returns values'
      write ( *, '(a)' ) '  of the spherical Bessel y0 function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X                     FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_y0_spherical_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f24.16,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test034 ( )

c*********************************************************************72
c
cc TEST034 demonstrates the use of BESSEL_Y1_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST034:'
      write ( *, '(a)' ) '  BESSEL_Y1_VALUES returns values of '
      write ( *, '(a)' ) '  the Bessel Y1 function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X                     FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_y1_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f24.16,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test035 ( )

c*********************************************************************72
c
cc TEST035 demonstrates the use of BESSEL_Y1_SPHERICAL_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST035:'
      write ( *, '(a)' ) '  BESSEL_Y1_SPHERICAL_VALUES returns values'
      write ( *, '(a)' ) '  of the spherical Bessel y1 function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X                     FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_y1_spherical_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f24.16,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test036 ( )

c*********************************************************************72
c
cc TEST036 demonstrates the use of BESSEL_YN_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST036:'
      write ( *, '(a)' ) '  BESSEL_YN_VALUES returns values of '
      write ( *, '(a)' ) '  the Bessel Yn function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N          X                     FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_yn_values ( n_data, n, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i6,2x,f24.16,2x,g24.16)' ) n, x, fx

      go to 10

20    continue

      return
      end
      subroutine test0365 ( )

c*********************************************************************72
c
cc TEST0365 demonstrates the use of BESSEL_YX_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      double precision n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0365:'
      write ( *, '(a)' ) '  BESSEL_YX_VALUES returns values of '
      write ( *, '(a)' ) '  the Bessel Yn function with NONINTEGER'
      write ( *, '(a)' ) '  order.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '          N             X                     FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_yx_values ( n_data, n, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f12.6,2x,f24.16,2x,g24.16)' ) n, x, fx

      go to 10

20    continue

      return
      end
      subroutine test037 ( )

c*********************************************************************72
c
cc TEST037 demonstrates the use of BETA_CDF_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision b
      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST037:'
      write ( *, '(a)' ) '  BETA_CDF_VALUES returns values of '
      write ( *, '(a)' ) '  the Beta CDF.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) 
     &  '      A             B                 X',
     &  '                     CDF'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call beta_cdf_values ( n_data, a, b, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f12.8,2x,f12.8,2x,f24.16,2x,g24.16)' ) 
     &    a, b, x, fx

      go to 10

20    continue

      return
      end
      subroutine test038 ( )

c*********************************************************************72
c
cc TEST038 demonstrates the use of BETA_INC_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision b
      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST038:'
      write ( *, '(a)' ) '  BETA_INC_VALUES returns values of '
      write ( *, '(a)' ) '  the incomplete Beta function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '          A               B               X           CDF'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call beta_inc_values ( n_data, a, b, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f12.8,2x,f12.8,2x,f24.16,2x,g24.16)' ) 
     &    a, b, x, fx

      go to 10

20    continue

      return
      end
      subroutine test039 ( )

c*********************************************************************72
c
cc TEST039 demonstrates the use of BETA_LOG_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fxy
      integer n_data
      double precision x
      double precision y

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST039:'
      write ( *, '(a)' ) '  BETA_LOG_VALUES returns values of '
      write ( *, '(a)' ) '  the logarithm of the Beta function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '          X               Y           Log(BETA(X,Y))'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call beta_log_values ( n_data, x, y, fxy )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f12.8,2x,f12.8,2x,g24.16)' ) 
     &    x, y, fxy

      go to 10

20    continue

      return
      end
      subroutine test0395 ( )

c*********************************************************************72
c
cc TEST0395 demonstrates the use of BETA_NONCENTRAL_CDF_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 January 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision b
      double precision fx
      double precision lambda
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0395:'
      write ( *, '(a)' ) '  BETA_NONCENTRAL_CDF_VALUES returns values'
      write ( *, '(a)' ) '  of the noncentral Beta CDF.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) 
     &  '      A             B           LAMBDA              ',
     &  'X                     CDF'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call beta_noncentral_cdf_values ( n_data, a, b, lambda, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f12.8,2x,f12.8,2x,f12.8,2x,f24.16,2x,g24.16)' )
     &    a, b, lambda, x, fx

      go to 10

20    continue

      return
      end
      subroutine test040

c*********************************************************************72
c
cc TEST040 demonstrates the use of BETA_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fxy
      integer n_data
      double precision x
      double precision y

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST040:'
      write ( *, '(a)' ) '  BETA_VALUES returns values of '
      write ( *, '(a)' ) '  the Beta function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '          X               Y            BETA(X,Y)'
      write ( *, '(a)' ) ' '

      n_data = 0
 
10    continue

        call beta_values ( n_data, x, y, fxy )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,f14.6,2x,g24.16)' ) x, y, fxy

      go to 10

20    continue

      return
      end
      subroutine test041

c*********************************************************************72
c
cc TEST041 demonstrates the use of BINOMIAL_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer a
      integer b
      integer c
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST041:'
      write ( *, '(a)' ) '  BINOMIAL_VALUES returns values of '
      write ( *, '(a)' ) '  the binomial numbers.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         A         B        C(A,B)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call binomial_values ( n_data, a, b, c )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,i8,2x,i12)' ) a, b, c

      go to 10

20    continue

      return
      end
      subroutine test042

c*********************************************************************72
c
cc TEST042 demonstrates the use of BINOMIAL_CDF_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer a
      double precision b
      double precision fx
      integer n_data
      integer x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST042:'
      write ( *, '(a)' ) '  BINOMIAL_CDF_VALUES returns values of '
      write ( *, '(a)' ) 
     &  '  the Binomial Cumulative Density Function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         A        B            X       CDF'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call binomial_cdf_values ( n_data, a, b, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,f10.4,2x,i8,2x,g24.16)' ) a, b, x, fx

      go to 10

20    continue

      return
      end
      subroutine test0425

c*********************************************************************72
c
cc TEST0425 demonstrates the use of BIVARIATE_NORMAL_CDF_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 May 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fxy
      integer n_data
      double precision r
      double precision x
      double precision y

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0425:'
      write ( *, '(a)' ) '  BIVARIATE_NORMAL_CDF_VALUES returns values'
      write ( *, '(a)' ) '  of the bivariate normal CDF.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     & '          X               Y               R           F(R)(X,Y)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bivariate_normal_cdf_values ( n_data, x, y, r, fxy )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,f14.6,2x,f14.6,2x,g14.6)' )
     &    x, y, r, fxy

      go to 10

20    continue

      return
      end
      subroutine test043

c*********************************************************************72
c
cc TEST043 demonstrates the use of CATALAN_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer c
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST043:'
      write ( *, '(a)' ) '  CATALAN_VALUES returns values of '
      write ( *, '(a)' ) '  the Catalan numbers.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         N          C(N)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call catalan_values ( n_data, n, c )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,i12)' ) n, c

      go to 10

20    continue

      return
      end
      subroutine test044

c*********************************************************************72
c
cc TEST044 demonstrates the use of CAUCHY_CDF_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      double precision mu
      integer n_data
      double precision sigma
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST044:'
      write ( *, '(a)' ) '  CAUCHY_CDF_VALUES returns values of '
      write ( *, '(a)' ) '  the Cauchy Cumulative Density Function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      MU     SIGMA    X                  CDF'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call cauchy_cdf_values ( n_data, mu, sigma, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g24.16)' ) 
     &    mu, sigma, x, fx

      go to 10

20    continue

      return
      end
      subroutine test045

c*********************************************************************72
c
cc TEST045 demonstrates the use of CHEBY_T_POLY_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST045:'
      write ( *, '(a)' ) '  CHEBY_T_POLY_VALUES returns values of '
      write ( *, '(a)' ) '  the Chebyshev T polynomials.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N       X      T(N)(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call cheby_t_poly_values ( n_data, n, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,f10.4,2x,g24.16)' ) n, x, fx

      go to 10

20    continue

      return
      end
      subroutine test046

c*********************************************************************72
c
cc TEST046 demonstrates the use of CHEBY_U_POLY_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST046:'
      write ( *, '(a)' ) '  CHEBY_U_POLY_VALUES returns values of '
      write ( *, '(a)' ) '  the Chebyshev U polynomials.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N       X      U(N)(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call cheby_u_poly_values ( n_data, n, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,f10.4,2x,g24.16)' ) n, x, fx

      go to 10

20    continue

      return
      end
      subroutine test047

c*********************************************************************72
c
cc TEST047 demonstrates the use of CHI_SQUARE_CDF_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer a
      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST047:'
      write ( *, '(a)' ) '  CHI_SQUARE_CDF_VALUES returns values of the'
      write ( *, '(a)' ) '  Chi-Squared Cumulative Density Function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N       X    CDF'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call chi_square_cdf_values ( n_data, a, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,f10.4,2x,g24.16)' ) a, x, fx

      go to 10

20    continue

      return
      end
      subroutine test048

c*********************************************************************72
c
cc TEST048 demonstrates the use of CHI_SQUARE_NONCENTRAL_CDF_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 December 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer df
      double precision fx
      double precision lambda
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST048:'
      write ( *, '(a)' ) '  CHI_SQUARE_NONCENTRAL_CDF_VALUES returns'
      write ( *, '(a)' ) '  values of the noncentral Chi-Squared '
      write ( *, '(a)' ) '  Cumulative Density Function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     DF      LAMBDA       X     CDF'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call chi_square_noncentral_cdf_values ( n_data, df, lambda, 
     &    x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,f10.4,2x,g24.16,2x,g24.16)' ) 
     &    df, lambda, x, fx

      go to 10

20    continue

      return
      end
      subroutine test049

c*********************************************************************72
c
cc TEST049 demonstrates the use of CI_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST049:'
      write ( *, '(a)' ) '  CI_VALUES returns values of '
      write ( *, '(a)' ) '  the Cosine Integral function CI(X).'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X            CI(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call ci_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test050

c*********************************************************************72
c
cc TEST050 demonstrates the use of CIN_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST050:'
      write ( *, '(a)' ) '  CIN_VALUES returns values of '
      write ( *, '(a)' ) '  the Cosine Integral function CIN(X).'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X            CIN(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call cin_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test051

c*********************************************************************72
c
cc TEST051 demonstrates the use of CLAUSEN_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST051:'
      write ( *, '(a)' ) '  CLAUSEN_VALUES returns values of '
      write ( *, '(a)' ) '  Clausen''s integral.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X            FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call clausen_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test05125

c*********************************************************************72
c
cc TEST05125 demonstrates the use of CLEBSCH_GORDAN_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 February 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      double precision j1
      double precision j2
      double precision j3
      double precision m1
      double precision m2
      double precision m3
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05125:'
      write ( *, '(a)' ) '  CLEBSCH_GORDAN_VALUES returns values of '
      write ( *, '(a)' ) '  the Clebsch-Gordan coefficient.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) 
     &  '      J1      J2      J3      ',
     &  'M1      M2      M3        CG'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call clebsch_gordan_values ( n_data, j1, j2, j3, m1, m2, 
     &  m3, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, 
     &  '(2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,g24.16)' ) 
     &    j1, j2, j3, m1, m2, m3, fx

      go to 10

20    continue

      return
      end
      subroutine test0515

c*********************************************************************72
c
cc TEST0515 demonstrates the use of COLLATZ_COUNT_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 March 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer count
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0515:'
      write ( *, '(a)' ) '  COLLATZ_COUNT_VALUES returns values of '
      write ( *, '(a)' ) '  the length of the Collatz sequence that'
      write ( *, '(a)' ) '  starts at N.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         N      COLLATZ_COUNT(N)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call collatz_count_values ( n_data, n, count )
    
        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,i12)' ) n, count

      go to 10

20    continue

      return
      end
      subroutine test052

c*********************************************************************72
c
cc TEST052 demonstrates the use of CP_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision cp
      integer n_data
      double precision p
      double precision tc

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST052:'
      write ( *, '(a)' ) '  CP_VALUES returns values of '
      write ( *, '(a)' ) '  the specific heat CP '
      write ( *, '(a)' ) '  as a function of temperature and pressure.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      T            P            CP(T,P)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call cp_values ( n_data, tc, p, cp )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,f14.6,2x,g14.6)' ) tc, p, cp

      go to 10

20    continue

      return
      end
      subroutine test053

c*********************************************************************72
c
cc TEST053 demonstrates the use of DAWSON_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST053:'
      write ( *, '(a)' ) '  DAWSON_VALUES returns values of '
      write ( *, '(a)' ) '  Dawson''s Integral function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X          DAWSON(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call dawson_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test054

c*********************************************************************72
c
cc TEST054 demonstrates the use of DEBYE1_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST054:'
      write ( *, '(a)' ) '  DEBYE1_VALUES returns values of '
      write ( *, '(a)' ) '  Debye''s function of order 1.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X          DEBYE1(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call debye1_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test055

c*********************************************************************72
c
cc TEST055 demonstrates the use of DEBYE2_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST055:'
      write ( *, '(a)' ) '  DEBYE2_VALUES returns values of '
      write ( *, '(a)' ) '  Debye''s function of order 2.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X          DEBYE2(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call debye2_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test056

c*********************************************************************72
c
cc TEST056 demonstrates the use of DEBYE3_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST056:'
      write ( *, '(a)' ) '  DEBYE3_VALUES returns values of '
      write ( *, '(a)' ) '  Debye''s function of order 3.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X          DEBYE3(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call debye3_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test057

c*********************************************************************72
c
cc TEST057 demonstrates the use of DEBYE4_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST057:'
      write ( *, '(a)' ) '  DEBYE4_VALUES returns values of '
      write ( *, '(a)' ) '  Debye''s function of order 4.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X          DEBYE4(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call debye4_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test058

c*********************************************************************72
c
cc TEST058 demonstrates the use of DIELECTRIC_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision eps
      integer n_data
      double precision p
      double precision tc

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST058:'
      write ( *, '(a)' ) '  DIELECTRIC_VALUES returns values of '
      write ( *, '(a)' ) '  the dielectric function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      T           P            EPS(T,P)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call dielectric_values ( n_data, tc, p, eps )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,f14.6,2x,g14.6)' ) tc, p, eps

      go to 10

20    continue

      return
      end
      subroutine test059

c*********************************************************************72
c
cc TEST059 demonstrates the use of DILOGARITHM_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST059:'
      write ( *, '(a)' ) '  DILOGARITHM_VALUES returns values of'
      write ( *, '(a)' ) '  the dilogarithm function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X          DILOGARITHM(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call dilogarithm_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test060

c*********************************************************************72
c
cc TEST060 demonstrates the use of E1_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST060:'
      write ( *, '(a)' ) '  E1_VALUES returns values of'
      write ( *, '(a)' ) '  the exponential integral function E1(X).'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X          E1(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call e1_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test061

c*********************************************************************72
c
cc TEST061 demonstrates the use of EI_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST061:'
      write ( *, '(a)' ) '  EI_VALUES returns values of'
      write ( *, '(a)' ) '  the exponential integral function EI(X).'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X          EI(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call ei_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test062

c*********************************************************************72
c
cc TEST062 demonstrates the use of ELLIPTIC_EA_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST062:'
      write ( *, '(a)' ) '  ELLIPTIC_EA_VALUES returns values of '
      write ( *, '(a)' ) 
     &  '  the complete elliptic integral of the second'
      write ( *, '(a)' ) 
     &  '  kind, with parameter angle ALPHA in degrees.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    ALPHA        EA(ALPHA)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call elliptic_ea_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test063

c*********************************************************************72
c
cc TEST063 demonstrates the use of ELLIPTIC_EM_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST063:'
      write ( *, '(a)' ) '  ELLIPTIC_EM_VALUES returns values of '
      write ( *, '(a)' ) 
     &  '  the complete elliptic integral of the second'
      write ( *, '(a)' ) '  kind, with parameter modulus M.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      M            EM(M)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call elliptic_em_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test064

c*********************************************************************72
c
cc TEST064 demonstrates the use of ELLIPTIC_KA_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST064:'
      write ( *, '(a)' ) '  ELLIPTIC_KA_VALUES returns values of '
      write ( *, '(a)' ) 
     &  '  the complete elliptic integral of the first'
      write ( *, '(a)' ) 
     &  '  kind, with parameter angle ALPHA in degrees.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    ALPHA        KA(ALPHA)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call elliptic_ka_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test065

c*********************************************************************72
c
cc TEST065 demonstrates the use of ELLIPTIC_KM_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST065:'
      write ( *, '(a)' ) '  ELLIPTIC_KM_VALUES returns values of '
      write ( *, '(a)' ) '  the complete elliptic integral of the first'
      write ( *, '(a)' ) '  kind, with parameter modulus M.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      M            KM(M)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call elliptic_km_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test066

c*********************************************************************72
c
cc TEST066 demonstrates the use of ERF_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST066:'
      write ( *, '(a)' ) '  ERF_VALUES returns values of '
      write ( *, '(a)' ) '  the error function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X               ERF(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call erf_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f12.8,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test0665

c*********************************************************************72
c
cc TEST0665 demonstrates the use of ERFC_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 May 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0665:'
      write ( *, '(a)' ) '  ERFC_VALUES returns values of '
      write ( *, '(a)' ) '  the complementary error function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X               ERFC(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call erfc_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f12.8,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test067

c*********************************************************************72
c
cc TEST067 demonstrates the use of EULER_NUMBER_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer c
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST067:'
      write ( *, '(a)' ) '  EULER_NUMBER_VALUES returns values of '
      write ( *, '(a)' ) '  the Euler numbers.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N        EULER_NUMBER(N)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call euler_number_values ( n_data, n, c )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,i12)' ) n, c

      go to 10

20    continue

      return
      end
      subroutine test068

c*********************************************************************72
c
cc TEST068 demonstrates the use of EULER_POLY_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST068:'
      write ( *, '(a)' ) '  EULER_POLY_VALUES returns values of '
      write ( *, '(a)' ) '  the Euler polynomials.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N     X             EULER_POLY(N,X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call euler_poly_values ( n_data, n, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,f14.6,2x,g24.16)' ) n, x, fx

      go to 10

20    continue

      return
      end
      subroutine test0685

c*********************************************************************72
c
cc TEST0685 demonstrates the use of EXP_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 March 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0685:'
      write ( *, '(a)' ) '  EXP_VALUES returns values of '
      write ( *, '(a)' ) '  the exponential function'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X           exp(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call exp_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test069

c*********************************************************************72
c
cc TEST069 demonstrates the use of EXP3_INT_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST069:'
      write ( *, '(a)' ) '  EXP3_INT_VALUES returns values of '
      write ( *, '(a)' ) '  the EXP3 Integral function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X          EXP3_INT(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call exp3_int_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test070

c*********************************************************************72
c
cc TEST070 demonstrates the use of EXPONENTIAL_CDF_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      double precision lambda
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST070:'
      write ( *, '(a)' ) '  EXPONENTIAL_CDF_VALUES returns values of'
      write ( *, '(a)' ) '  the Exponential CDF.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      LAMBDA       X          CDF'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call exponential_cdf_values ( n_data, lambda, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,f14.6,2x,g24.16)' ) lambda, x, fx

      go to 10

20    continue

      return
      end
      subroutine test071

c*********************************************************************72
c
cc TEST071 demonstrates the use of EXTREME_VALUES_CDF_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision alpha
      double precision beta
      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST071:'
      write ( *, '(a)' ) '  EXTREME_VALUES_CDF_VALUES returns values'
      write ( *, '(a)' ) 
     &  '  of the Extreme Values Cumulative Density Function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      ALPHA  BETA     X                  CDF'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call extreme_values_cdf_values ( n_data, alpha, beta, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g24.16)' ) 
     &    alpha, beta, x, fx

      go to 10

20    continue

      return
      end
      subroutine test072

c*********************************************************************72
c
cc TEST072 demonstrates the use of F_CDF_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer a
      integer b
      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST072:'
      write ( *, '(a)' ) '  F_CDF_VALUES returns values of'
      write ( *, '(a)' ) '  the F cumulative density function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       A       B      X            CDF'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call f_cdf_values ( n_data, a, b, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,i8,2x,f14.6,2x,g24.16)' ) a, b, x, fx

      go to 10

20    continue

      return
      end
      subroutine test073

c*********************************************************************72
c
cc TEST073 demonstrates the use of F_NONCENTRAL_CDF_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer a
      integer b
      double precision fx
      double precision lambda
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST073:'
      write ( *, '(a)' ) '  F_NONCENTRAL_CDF_VALUES returns values of'
      write ( *, '(a)' ) '  the F cumulative density function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '       A       B      LAMBDA    X            CDF'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call f_noncentral_cdf_values ( n_data, a, b, lambda, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,i8,2x,f10.6,2x,g14.6,2x,g14.6)' ) 
     &    a, b, lambda, x, fx

      go to 10

20    continue

      return
      end
      subroutine test074

c*********************************************************************72
c
cc TEST074 demonstrates the use of FRESNEL_COS_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST074:'
      write ( *, '(a)' ) '  FRESNEL_COS_VALUES returns values of'
      write ( *, '(a)' ) '  the Fresnel cosine integral C(X).'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X           C(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call fresnel_cos_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test075

c*********************************************************************72
c
cc TEST075 demonstrates the use of FRESNEL_SIN_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST075:'
      write ( *, '(a)' ) '  FRESNEL_SIN_VALUES returns values of'
      write ( *, '(a)' ) '  the Fresnel sine integral S(X).'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X           S(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call fresnel_sin_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test0755

c***********************************************************************72
c
cc TEST0755 demonstrates the use of FROBENIUS_NUMBER_ORDER2_VALUES.
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
      integer f
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0755:'
      write ( *, '(a)' ) '  FROBENIUS_NUMBER_ORDER2_VALUES returns'
      write ( *, '(a)' ) '  values of the Frobenius number of order 2.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        C1        C2     F(C1,C2)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call frobenius_number_order2_values ( n_data, c1, c2, f )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

       write ( *, '(2x,i8,2x,i8,2x,i8)' ) c1, c2, f

      go to 10

20    continue

      return
      end
      subroutine test076

c*********************************************************************72
c
cc TEST076 demonstrates the use of GAMMA_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST076:'
      write ( *, '(a)' ) 
     &  '  GAMMA_VALUES returns values of the Gamma function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X            GAMMA(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call gamma_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test077

c*********************************************************************72
c
cc TEST077 demonstrates the use of GAMMA_CDF_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      double precision mu
      integer n_data
      double precision sigma
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST077:'
      write ( *, '(a)' ) '  GAMMA_CDF_VALUES returns values of '
      write ( *, '(a)' ) '  the GAMMA Cumulative Density Function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      MU     SIGMA    X                  CDF'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call gamma_cdf_values ( n_data, mu, sigma, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g24.16)' ) 
     &    mu, sigma, x, fx

      go to 10

20    continue

      return
      end
      subroutine test078

c*********************************************************************72
c
cc TEST078 demonstrates the use of GAMMA_INC_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST078:'
      write ( *, '(a)' ) '  GAMMA_INC_VALUES returns values of '
      write ( *, '(a)' ) '  the incomplete Gamma function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '          A               X           CDF'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call gamma_inc_values ( n_data, a, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f12.8,2x,f24.16,2x,g24.16)' ) 
     &    a, x, fx

      go to 10

20    continue

      return
      end
      subroutine test079

c*********************************************************************72
c
cc TEST079 demonstrates the use of GAMMA_LOG_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST079:'
      write ( *, '(a)' ) '  GAMMA_LOG_VALUES returns values of '
      write ( *, '(a)' ) '  the logarithm of the Gamma function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X           FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call gamma_log_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f24.16,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test080

c*********************************************************************72
c
cc TEST080 demonstrates the use of GEGENBAUER_POLY_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision fx
      integer n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST080:'
      write ( *, '(a)' ) '  GEGENBAUER_POLY_VALUES returns values of '
      write ( *, '(a)' ) '  the Gegenbauer polynomials.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         N       A            X     G(N,A)(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call gegenbauer_poly_values ( n_data, n, a, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,f10.4,2x,f10.4,2x,g14.6)' ) n, a, x, fx

      go to 10

20    continue

      return
      end
      subroutine test081

c*********************************************************************72
c
cc TEST081 demonstrates the use of GEOMETRIC_CDF_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision cdf
      integer n_data
      double precision p
      integer x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST081:'
      write ( *, '(a)' ) '  GEOMETRIC_CDF_VALUES returns values of '
      write ( *, '(a)' ) '  the Geometric Probability '
      write ( *, '(a)' ) '  Cumulative Density Function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       X      P       CDF'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call geometric_cdf_values ( n_data, x, p, cdf )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,f10.4,2x,g24.16)' ) x, p, cdf

      go to 10

20    continue

      return
      end
      subroutine test082

c*********************************************************************72
c
cc TEST082 demonstrates the use of GOODWIN_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST082:'
      write ( *, '(a)' ) '  GOODWIN_VALUES returns values of'
      write ( *, '(a)' ) '  the Goodwin and Staton function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X          F(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call goodwin_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test083

c*********************************************************************72
c
cc TEST083 demonstrates the use of GUD_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST083:'
      write ( *, '(a)' ) '  GUD_VALUES returns values of '
      write ( *, '(a)' ) '  the Gudermannian function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X            GUD(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call gud_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test0835

c*********************************************************************72
c
cc TEST0835 demonstrates the use of HERMITE_FUNCTION_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0835:'
      write ( *, '(a)' ) '  HERMITE_FUNCTION_VALUES returns values of '
      write ( *, '(a)' ) '  the Hermite function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N      X            Hf(N,X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call hermite_function_values ( n_data, n, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,g14.6,2x,g24.16)' ) n, x, fx

      go to 10

20    continue

      return
      end
      subroutine test084

c*********************************************************************72
c
cc TEST084 demonstrates the use of HERMITE_POLY_PHYS_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST084:'
      write ( *, '(a)' ) '  HERMITE_POLY_PHYS_VALUES returns values of '
      write ( *, '(a)' ) '  the physicist''s Hermite polynomials.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N      X            H(N,X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call hermite_poly_phys_values ( n_data, n, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,g14.6,2x,g24.16)' ) n, x, fx

      go to 10

20    continue

      return
      end
      subroutine test0843

c*********************************************************************72
c
cc TEST0843 demonstrates the use of HERMITE_POLY_PROB_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0843:'
      write ( *, '(a)' ) '  HERMITE_POLY_PROB_VALUES returns values of '
      write ( *, '(a)' ) '  the probabilist''s Hermite polynomials.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N      X            He(N,X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call hermite_poly_prob_values ( n_data, n, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,g14.6,2x,g24.16)' ) n, x, fx

      go to 10

20    continue

      return
      end
      subroutine test0845

c*********************************************************************72
c
cc TEST0845 tests HYPER_2F1_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 September 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision b
      double precision c
      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0845:'
      write ( *, '(a)' ) '  HYPER_2F1_VALUES returns values of '
      write ( *, '(a)' ) '  the hypergeometric 2F1 function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) 
     &  '         A           B           C            ',
     &  'X       Hyper_2F1(A,B,C,X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call hyper_2f1_values ( n_data, a, b, c, x, fx )

        if ( n_data == 0 ) then
          go to 20
        end if

        write ( *, '(2x,f10.6,2x,f10.6,2x,f10.6,2x,f10.6,2x,g24.16)' ) 
     &    a, b, c, x, fx

      go to 10

20    continue

      return
      end
      subroutine test085

c*********************************************************************72
c
cc TEST085 demonstrates the use of HYPERGEOMETRIC_CDF_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      integer pop
      integer sam
      integer suc
      integer x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST085:'
      write ( *, '(a)' ) '  HYPERGEOMETRIC_CDF_VALUES returns values'
      write ( *, '(a)' ) '  of the Hypergeometric '
      write ( *, '(a)' ) '  Cumulative Density Function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '       SAM      SUC     POP       X   HyperCDF(S,S,P)(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call hypergeometric_cdf_values ( n_data, sam, suc, pop, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,i8,2x,i8,2x,i8,2x,g24.16)' ) 
     &  sam, suc, pop, x, fx

      go to 10

20    continue

      return
      end
      subroutine test0855

c*********************************************************************72
c
cc TEST0855 demonstrates the use of HYPERGEOMETRIC_PDF_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 January 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      integer pop
      integer sam
      integer suc
      integer x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0855:'
      write ( *, '(a)' ) '  HYPERGEOMETRIC_PDF_VALUES returns values'
      write ( *, '(a)' ) '  of the Hypergeometric '
      write ( *, '(a)' ) '  Cumulative Density Function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '       SAM      SUC     POP       X   HyperPDF(S,S,P)(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call hypergeometric_pdf_values ( n_data, sam, suc, pop, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,i8,2x,i8,2x,i8,2x,g24.16)' ) 
     &  sam, suc, pop, x, fx

      go to 10

20    continue

      return
      end
      subroutine test086

c*********************************************************************72
c
cc TEST086 demonstrates the use of FACTORIAL_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer fn
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST086:'
      write ( *, '(a)' ) '  FACTORIAL_VALUES returns values of '
      write ( *, '(a)' ) '  the factorial function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N         Nc'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call factorial_values ( n_data, n, fn )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,i12)' ) n, fn

      go to 10

20    continue

      return
      end
      subroutine test087

c*********************************************************************72
c
cc TEST087 demonstrates the use of FACTORIAL2_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer fn
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST087:'
      write ( *, '(a)' ) '  FACTORIAL2_VALUES returns values of '
      write ( *, '(a)' ) '  the double factorial function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N         Ncc'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call factorial2_values ( n_data, n, fn )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,i12)' ) n, fn

      go to 10

20    continue

      return
      end
      subroutine test088

c*********************************************************************72
c
cc TEST088 demonstrates the use of FACTORIAL_RISING_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer fmn
      integer m
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST088:'
      write ( *, '(a)' ) '  FACTORIAL_RISING_VALUES returns some exact'
      write ( *, '(a)' ) '  values of the rising factorial function:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       M       N      Fctorial_rising(M,N)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call factorial_rising_values ( n_data, m, n, fmn )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,i8,2x,i12)' ) m, n, fmn

      go to 10

20    continue

      return
      end
      subroutine test089

c*********************************************************************72
c
cc TEST089 demonstrates the use of I0ML0_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST089:'
      write ( *, '(a)' ) '  I0ML0_VALUES returns values of'
      write ( *, '(a)' ) '  the I0ML0 function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X          F(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call i0ml0_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test090

c*********************************************************************72
c
cc TEST090 demonstrates the use of I1ML1_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST090:'
      write ( *, '(a)' ) '  I1ML1_VALUES returns values of'
      write ( *, '(a)' ) '  the I1ML1 function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X          F(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call i1ml1_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test091

c*********************************************************************72
c
cc TEST091 demonstrates the use of JACOBI_CN_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST091:'
      write ( *, '(a)' ) '  JACOBI_CN_VALUES returns values of '
      write ( *, '(a)' ) '  the Jacobi elliptic CN function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      A         X       CN(A,X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call jacobi_cn_values ( n_data, a, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f10.4,2x,f10.4,2x,g24.16)' ) a, x, fx

      go to 10

20    continue

      return
      end
      subroutine test092

c*********************************************************************72
c
cc TEST092 demonstrates the use of JACOBI_DN_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST092:'
      write ( *, '(a)' ) '  JACOBI_DN_VALUES returns values of '
      write ( *, '(a)' ) '  the Jacobi elliptic DN function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      A         X       DN(A,X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call jacobi_dn_values ( n_data, a, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f10.4,2x,f10.4,2x,g24.16)' ) a, x, fx

      go to 10

20    continue

      return
      end
      subroutine test093 ( )

c*********************************************************************72
c
cc TEST093 demonstrates the use of JACOBI_POLY_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 April 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision b
      double precision fx
      integer n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST093:'
      write ( *, '(a)' ) '  JACOBI_POLY_VALUES returns values of '
      write ( *, '(a)' ) '  the Jacobi polynomial.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '       N       A       B      X       J(N,A,B)(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call jacobi_poly_values ( n_data, n, a, b, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i6,2x,f8.4,2x,f8.4,2x,f24.16,2x,g24.16)' ) 
     &    n, a, b, x, fx

      go to 10

20    continue

      return
      end
      subroutine test094 ( )

c*********************************************************************72
c
cc TEST094 demonstrates the use of JACOBI_SN_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST094:'
      write ( *, '(a)' ) '  JACOBI_SN_VALUES returns values of '
      write ( *, '(a)' ) '  the Jacobi elliptic SN function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      A         X       SN(A,X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call jacobi_sn_values ( n_data, a, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f10.4,2x,f10.4,2x,g14.6)' ) a, x, fx

      go to 10

20    continue

      return
      end
      subroutine test0945 ( )

c*********************************************************************72
c
cc TEST0945 demonstrates the use of JED_CE_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 May 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer d
      double precision f
      double precision jed
      integer n_data
      integer m
      integer y

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0945:'
      write ( *, '(a)' ) '  JED_CE_VALUES returns:'
      write ( *, '(a)' ) '  JED, a Julian Ephemeris Date, and'
      write ( *, '(a)' ) 
     &  '  YMDF, the corresponding year, month, day, fraction.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        JED          Y   M   D    F'
      write ( *, '(a)' ) ' '

      n_data = 0
    
10    continue

        call jed_ce_values ( n_data, jed, y, m, d, f )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f12.2,2x,i6,2x,i2,2x,i2,2x,f6.4)' ) 
     &    jed, y, m, d, f

      go to 10

20    continue

      return
      end
      subroutine test095 ( )

c*********************************************************************72
c
cc TEST095 demonstrates the use of JED_MJD_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision jed
      integer n_data
      double precision mjd

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST095:'
      write ( *, '(a)' ) '  JED_MJD_VALUES returns:'
      write ( *, '(a)' ) '  JED, a Julian Ephemeris Date, and'
      write ( *, '(a)' ) 
     &  '  MJD, the corresponding Modified Julian Day count.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        JED           MJD'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call jed_mjd_values ( n_data, jed, mjd )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f12.2,2x,f12.2)' ) jed, mjd

      go to 10

20    continue

      return
      end
      subroutine test096 ( )

c*********************************************************************72
c
cc TEST096 demonstrates the use of JED_RD_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision jed
      integer n_data
      double precision rd

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST096:'
      write ( *, '(a)' ) '  JED_RD_VALUES returns:'
      write ( *, '(a)' ) '  JED, a Julian Ephemeris Date, and'
      write ( *, '(a)' ) '  RD, the corresponding '
      write ( *, '(a)' ) '  Reingold Dershowitz Day count.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        JED            RD'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call jed_rd_values ( n_data, jed, rd )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f12.2,2x,f12.2)' ) jed, rd

      go to 10

20    continue

      return
      end
      subroutine test097 ( )

c*********************************************************************72
c
cc TEST097 demonstrates the use of JED_WEEKDAY_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision jed
      integer n_data
      integer weekday
      character*9 weekday_name(7)

      save weekday_name

      data weekday_name /
     &  'Sunday   ', 'Monday   ', 
     &  'Tuesday  ', 'Wednesday', 
     &  'Thursday ', 
     &  'Friday   ', 'Saturday ' /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST097:'
      write ( *, '(a)' ) '  JED_WEEKDAY_VALUES returns '
      write ( *, '(a)' ) '  Julian Ephemeris Dates '
      write ( *, '(a)' ) '  (JED) and the corresponding weekday'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        JED     #  Weekday'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call jed_weekday_values ( n_data, jed, weekday )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f12.2,2x,i1,2x,a9)' ) 
     &    jed, weekday, weekday_name(weekday)

      go to 10

20    continue

      return
      end
      subroutine test0972 ( )

c*********************************************************************72
c
cc TEST0972 demonstrates the use of KEI0_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0972:'
      write ( *, '(a)' ) '  KEI0_VALUES returns values of '
      write ( *, '(a)' ) '  the Kelvin KEI function of order 0.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X           FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call kei0_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f24.16,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test0973 ( )

c*********************************************************************72
c
cc TEST0973 demonstrates the use of KEI1_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0973:'
      write ( *, '(a)' ) '  KEI1_VALUES returns values of '
      write ( *, '(a)' ) '  the Kelvin KEI function of order 1.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X           FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call kei1_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f24.16,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test0974 ( )

c*********************************************************************72
c
cc TEST0974 demonstrates the use of KER0_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0974:'
      write ( *, '(a)' ) '  KER0_VALUES returns values of '
      write ( *, '(a)' ) '  the Kelvin KER function of order 0.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X           FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call ker0_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f24.16,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test0975 ( )

c*********************************************************************72
c
cc TEST0975 demonstrates the use of KER1_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0975:'
      write ( *, '(a)' ) '  KER1_VALUES returns values of '
      write ( *, '(a)' ) '  the Kelvin KER function of order 1.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X           FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call ker1_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f24.16,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test098 ( )

c*********************************************************************72
c
cc TEST098 demonstrates the use of LAGUERRE_ASSOCIATED_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer m
      integer n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST098:'
      write ( *, '(a)' ) '  LAGUERRE_ASSOCIATED_VALUES returns values'
      write ( *, '(a)' ) '  of the associated Laguerre polynomials.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N       M      X            L(N,M)(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call laguerre_associated_values ( n_data, n, m, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,i8,2x,f14.6,2x,g24.16)' ) n, m, x, fx

      go to 10

20    continue

      return
      end
      subroutine test099 ( )

c*********************************************************************72
c
cc TEST099 demonstrates the use of LAGUERRE_POLYNOMIAL_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST099:'
      write ( *, '(a)' ) '  LAGUERRE_POLYNOMIAL_VALUES returns values'
      write ( *, '(a)' ) '  of the Laguerre polynomials.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N     X            L(N)(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call laguerre_polynomial_values ( n_data, n, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,f14.6,2x,g24.16)' ) n, x, fx

      go to 10

20    continue

      return
      end
      subroutine test0995 ( )

c*********************************************************************72
c
cc TEST0995 demonstrates the use of LAMBERT_W_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0995:'
      write ( *, '(a)' ) '  LAMBERT_W_VALUES returns values of '
      write ( *, '(a)' ) '  the Lambert W function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '          X               W(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call lambert_w_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f12.8,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test100 ( )

c*********************************************************************72
c
cc TEST100 demonstrates the use of LAPLACE_CDF_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision beta
      double precision fx
      double precision mu
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST100:'
      write ( *, '(a)' ) '  LAPLACE_CDF_VALUES returns values of '
      write ( *, '(a)' ) '  the Laplace CDF.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     MU     BETA      X            CDF'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call laplace_cdf_values ( n_data, mu, beta, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,f14.6,2x,f14.6,2x,g24.16)' ) 
     &    mu, beta, x, fx

      go to 10

20    continue

      return
      end
      subroutine test101 ( )

c*********************************************************************72
c
cc TEST101 demonstrates the use of LEGENDRE_ASSOCIATED_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer m
      integer n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST101:'
      write ( *, '(a)' ) '  LEGENDRE_ASSOCIATED_VALUES returns values'
      write ( *, '(a)' ) '  of the associated Legendre polynomials.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N       M    X             P(N,M)(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call legendre_associated_values ( n_data, n, m, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,i8,2x,f14.6,2x,g24.16)' ) n, m, x, fx

      go to 10

20    continue

      return
      end
      subroutine test1015 ( )

c*********************************************************************72
c
cc TEST1015 demonstrates the use of LEGENDRE_ASSOCIATED_NORMALIZED_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 September 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer m
      integer n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1015:'
      write ( *, '(a)' ) '  LEGENDRE_ASSOCIATED_NORMALIZED_VALUES'
      write ( *, '(a)' ) '  returns values of the normalized associated'
      write ( *, '(a)' ) '  Legendre polynomials.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N       M    X             P(N,M)(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call legendre_associated_normalized_values ( n_data, n, m, x, 
     &    fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,i8,2x,f14.6,2x,g24.16)' ) n, m, x, fx

      go to 10

20    continue

      return
      end
      subroutine test1016 ( )

c*********************************************************************72
c
cc TEST1016 demonstrates the use of LEGENDRE_ASSOCIATED_NORMALIZED_SPHERE_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 March 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer m
      integer n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1016:'
      write ( *, '(a)' ) 
     &  '  LEGENDRE_ASSOCIATED_NORMALIZED_SPHERE_VALUES'
      write ( *, '(a)' ) '  returns values of the associated Legendre'
      write ( *, '(a)' ) '  polynomials, normalized for a sphere'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N       M    X             P(N,M)(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call legendre_associated_normalized_sphere_values ( n_data, n, 
     &    m, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,i8,2x,f14.6,2x,g24.16)' ) n, m, x, fx

      go to 10

20    continue

      return
      end
      subroutine test102 ( )

c*********************************************************************72
c
cc TEST102 demonstrates the use of LEGENDRE_POLY_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST102:'
      write ( *, '(a)' ) '  LEGENDRE_POLY_VALUES returns values of '
      write ( *, '(a)' ) '  the Legendre PN polynomials.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N    X             P(N)(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call legendre_poly_values ( n_data, n, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,f14.6,2x,g24.16)' ) n, x, fx

      go to 10

20    continue

      return
      end
      subroutine test103 ( )

c*********************************************************************72
c
cc TEST103 demonstrates the use of LEGENDRE_FUNCTION_Q_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST103:'
      write ( *, '(a)' ) '  LEGENDRE_FUNCTION_Q_VALUES returns values'
      write ( *, '(a)' ) '  of the Legendre QN polynomials.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N    X             Q(N)(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call legendre_function_q_values ( n_data, n, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,f14.6,2x,g14.6)' ) n, x, fx

      go to 10

20    continue

      return
      end
      subroutine test1035 ( )

c*********************************************************************72
c
cc TEST1035 demonstrates the use of LERCH_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision fx
      integer n_data
      integer s
      double precision z

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1035:'
      write ( *, '(a)' ) '  LERCH_VALUES returns values of '
      write ( *, '(a)' ) '  the Lerch transcendent function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     Z        S        A       CDF'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call lerch_values ( n_data, z, s, a, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f10.4,2x,i8,2x,f10.4,2x,g24.16)' ) 
     &    z, s, a, fx

      go to 10

20    continue

      return
      end
      subroutine test104 ( )

c*********************************************************************72
c
cc TEST104 demonstrates the use of LOBACHEVSKY_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST104:'
      write ( *, '(a)' ) '  LOBACHEVSKY_VALUES returns values of'
      write ( *, '(a)' ) '  the Lobachevsky function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X          F(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call lobachevsky_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test1037 ( )

c*********************************************************************72
c
cc TEST1037 demonstrates the use of LOG_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1037:'
      write ( *, '(a)' ) '  LOG_VALUES returns values of'
      write ( *, '(a)' ) '  the natural logarithm function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X          F(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call log_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test105 ( )

c*********************************************************************72
c
cc TEST105 demonstrates the use of LOG_NORMAL_CDF_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      double precision mu
      integer n_data
      double precision sigma
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST105:'
      write ( *, '(a)' ) '  LOG_NORMAL_CDF_VALUES returns values of '
      write ( *, '(a)' ) '  the Log Normal Cumulative Density Function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      MU     SIGMA    X                  CDF'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call log_normal_cdf_values ( n_data, mu, sigma, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g24.16)' ) 
     &    mu, sigma, x, fx

      go to 10

20    continue

      return
      end
      subroutine test106 ( )

c*********************************************************************72
c
cc TEST106 demonstrates the use of LOG_SERIES_CDF_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n
      integer n_data
      double precision t

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST106:'
      write ( *, '(a)' ) '  LOG_SERIES_CDF_VALUES returns values of '
      write ( *, '(a)' ) '  the Log Series Cumulative Density Function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     T          N   CDF'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call log_series_cdf_values ( n_data, t, n, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f10.4,2x,i8,2x,g24.16)' ) t, n, fx

      go to 10

20    continue

      return
      end
      subroutine test107 ( )

c*********************************************************************72
c
cc TEST107 demonstrates the use of LOGARITHMIC_INTEGRAL_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST107:'
      write ( *, '(a)' ) '  LOGARITHMIC_INTEGAL_VALUES returns values'
      write ( *, '(a)' ) '  of the logarithmic integral function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X            LI(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call logarithmic_integral_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test108 ( )

c*********************************************************************72
c
cc TEST108 demonstrates the use of LOGISTIC_CDF_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision beta
      double precision fx
      double precision mu
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST108:'
      write ( *, '(a)' ) '  LOGISTIC_CDF_VALUES returns values of '
      write ( *, '(a)' ) '  the Logistic Cumulative Density Function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      MU     BETA     X                  CDF'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call logistic_cdf_values ( n_data, mu, beta, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g24.16)' ) 
     &    mu, beta, x, fx

      go to 10

20    continue

      return
      end
      subroutine test10825 ( )

c*********************************************************************72
c
cc TEST10825 demonstrates the use of MATHIEU_EVEN_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      integer n_data
      integer q
      integer r

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST10825:'
      write ( *, '(a)' ) '  MATHIEU_EVEN_VALUES returns values of the'
      write ( *, '(a)' ) '  eigenvalues of Mathieu''s differential'
      write ( *, '(a)' ) '  equation associated with even periodic'
      write ( *, '(a)' ) '  solutions.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       R       Q            A(R,Q)'
      write ( *, '(a)' ) ' '

      n_data = 0

 10   continue

      call mathieu_even_values ( n_data, r, q, a )

      if ( n_data .eq. 0 ) then
        go to 20
      end if

      write ( *, '(2x,i6,2x,i6,2x,g24.16)' ) r, q, a

      go to 10

20    continue

      return
      end
      subroutine test10850 ( )

c*********************************************************************72
c
cc TEST10850 demonstrates the use of MATHIEU_ODD_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision b
      integer n_data
      integer q
      integer r

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST10850:'
      write ( *, '(a)' ) '  MATHIEU_ODD_VALUES returns values of the'
      write ( *, '(a)' ) '  eigenvalues of Mathieu''s differential'
      write ( *, '(a)' ) '  equation associated with odd periodic'
      write ( *, '(a)' ) '  solutions.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       R       Q            B(R,Q)'
      write ( *, '(a)' ) ' '

      n_data = 0

 10   continue

      call mathieu_odd_values ( n_data, r, q, b )

      if ( n_data .eq. 0 ) then
        go to 20
      end if

      write ( *, '(2x,i6,2x,i6,2x,g24.16)' ) r, q, b

      go to 10

20    continue

      return
      end
      subroutine test109 ( )

c*********************************************************************72
c
cc TEST109 demonstrates the use of MOEBIUS_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer fn
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST109:'
      write ( *, '(a)' ) '  MOEBIUS_VALUES returns values of '
      write ( *, '(a)' ) '  the Moebius function.'
      write ( *, '(a)' ) ' ' 
      write ( *, '(a)' ) '       N         MU(N)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call moebius_values ( n_data, n, fn )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,i12)' ) n, fn

      go to 10

20    continue

      return
      end
      subroutine test110 ( )

c*********************************************************************72
c
cc TEST110 demonstrates the use of NEGATIVE_BINOMIAL_CDF_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision cdf
      integer f
      integer n_data
      double precision p
      integer s

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST110:'
      write ( *, '(a)' ) '  NEGATIVE_BINOMIAL_CDF_VALUES returns'
      write ( *, '(a)' ) '  values of the Negative Binomial '
      write ( *, '(a)' ) '  Cumulative Density Function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       F       S         P         CDF'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call negative_binomial_cdf_values ( n_data, f, s, p, cdf )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g14.6)' ) f, s, p, cdf

      go to 10

20    continue

      return
      end
      subroutine test1105 ( )

c*********************************************************************72
c
cc TEST1105 demonstrates the use of NINE_J_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 February 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      double precision j1
      double precision j2
      double precision j3
      double precision j4
      double precision j5
      double precision j6
      double precision j7
      double precision j8
      double precision j9
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1105:'
      write ( *, '(a)' ) '  NINE_J_VALUES returns values of '
      write ( *, '(a)' ) '  the Wigner 9J coefficient.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) 
     &  '      J1      J2      J3      J4      J5      J6',
     &  '      J7      J8      J9        NINE_J'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call nine_j_values ( n_data, j1, j2, j3, j4, j5, j6, 
     &    j7, j8, j9, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, 
     &  '(2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,
     &  f6.2,2x,f6.2,2x,f6.2,2x,g24.16)' ) 
     &  j1, j2, j3, j4, j5, j6, j7, j8, j9, fx

      go to 10

20    continue

      return
      end
      subroutine test111 ( )

c*********************************************************************72
c
cc TEST111 demonstrates the use of NORMAL_CDF_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      double precision mu
      integer n_data
      double precision sigma
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST111:'
      write ( *, '(a)' ) '  NORMAL_CDF_VALUES returns values of '
      write ( *, '(a)' ) '  the Normal Cumulative Density Function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      MU     SIGMA    X                  CDF'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call normal_cdf_values ( n_data, mu, sigma, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g24.16)' ) 
     &    mu, sigma, x, fx

      go to 10

20    continue

      return
      end
      subroutine test112 ( )

c*********************************************************************72
c
cc TEST112 demonstrates the use of NORMAL_01_CDF_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST112:'
      write ( *, '(a)' ) '  NORMAL_01_CDF_VALUES returns values of '
      write ( *, '(a)' ) '  the Normal Cumulative Density Function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X                  CDF'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call normal_01_cdf_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test113 ( )

c*********************************************************************72
c
cc TEST113 demonstrates the use of OMEGA_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer fn
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST113:'
      write ( *, '(a)' ) '  OMEGA_VALUES returns values of '
      write ( *, '(a)' ) '  the Omega function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     N           OMEGA(N)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call omega_values ( n_data, n, fn )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i12,2x,i12)' ) n, fn

      go to 10

20    continue

      return
      end
      subroutine test1135 ( )

c*********************************************************************72
c
cc TEST1135 demonstrates the use of OWEN_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision h
      integer n_data
      double precision t

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1135:'
      write ( *, '(a)' ) '  OWEN_VALUES returns values of '
      write ( *, '(a)' ) '  the Owen T function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      H       A       T'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call owen_values ( n_data, h, a, t )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,f14.6,2x,g24.16)' ) h, a, t

      go to 10

20    continue

      return
      end
      subroutine test114 ( )

c*********************************************************************72
c
cc TEST114 demonstrates the use of PARTITION_COUNT_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer fn
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST114:'
      write ( *, '(a)' ) '  PARTITION_COUNT_VALUES returns values of '
      write ( *, '(a)' ) '  the integer partition count function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N         P(N)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call partition_count_values ( n_data, n, fn )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,i12)' ) n, fn

      go to 10

20    continue

      return
      end
      subroutine test115 ( )

c*********************************************************************72
c
cc TEST115 demonstrates the use of PARTITION_DISTINCT_COUNT_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer fn
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST115:'
      write ( *, '(a)' ) '  PARTITION_DISTINCT_COUNT_VALUES returns '
      write ( *, '(a)' ) '  values of the integer distinct partition'
      write ( *, '(a)' ) '  count function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N         Q(N)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call partition_distinct_count_values ( n_data, n, fn )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,i12)' ) n, fn

      go to 10

20    continue

      return
      end
      subroutine test116 ( )

c*********************************************************************72
c
cc TEST116 demonstrates the use of PHI_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer fn
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST116:'
      write ( *, '(a)' ) '  PHI_VALUES returns values of '
      write ( *, '(a)' ) '  the PHI function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N         PHI(N)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call phi_values ( n_data, n, fn )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,i12)' ) n, fn

      go to 10

20    continue

      return
      end
      subroutine test117 ( )

c*********************************************************************72
c
cc TEST117 demonstrates the use of PI_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer fn
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST117:'
      write ( *, '(a)' ) '  PI_VALUES returns values of '
      write ( *, '(a)' ) '  the PI function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             N         PI(N)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call pi_values ( n_data, n, fn )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i12,2x,i12)' ) n, fn

      go to 10

20    continue

      return
      end
      subroutine test118 ( )

c*********************************************************************72
c
cc TEST118 demonstrates the use of POISSON_CDF_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision fx
      integer n_data
      integer x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST118:'
      write ( *, '(a)' ) '  POISSON_CDF_VALUES returns values of '
      write ( *, '(a)' ) '  the Poisson Cumulative Density Function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      A          X    CDF'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call poisson_cdf_values ( n_data, a, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f10.4,2x,i8,2x,g24.16)' ) a, x, fx

      go to 10

20    continue

      return
      end
      subroutine test1185 ( )

c*********************************************************************72
c
cc TEST1185 demonstrates the use of POLYLOGARITHM_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n
      integer n_data
      double precision z

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1185:'
      write ( *, '(a)' ) '  POLYLOGARITHM_VALUES returns values of '
      write ( *, '(a)' ) '  the polylogarithm.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N      Z        FX'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call polylogarithm_values ( n_data, n, z, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,f10.4,2x,g24.16)' ) n, z, fx

      go to 10

20    continue

      return
      end
      subroutine test119 ( )

c*********************************************************************72
c
cc TEST119 demonstrates the use of PRANDTL_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_data
      double precision p
      double precision pr
      double precision tc

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST119:'
      write ( *, '(a)' ) '  PRANDTL_VALUES returns values of '
      write ( *, '(a)' ) '  the Prandtl number of water '
      write ( *, '(a)' ) '  as a function of temperature and pressure.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      T            P            Pr(T,P)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call prandtl_values ( n_data, tc, p, pr )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,f14.6,2x,g14.6)' ) tc, p, pr

      go to 10

20    continue

      return
      end
      subroutine test120 ( )

c*********************************************************************72
c
cc TEST120 demonstrates the use of PRIME_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer n_data
      integer p

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST120:'
      write ( *, '(a)' ) '  PRIME_VALUES returns values of '
      write ( *, '(a)' ) '  the prime function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '           N          P[N]'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call prime_values ( n_data, n, p )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i12,2x,i12)' ) n, p

      go to 10

20    continue

      return
      end
      subroutine test121 ( )

c*********************************************************************72
c
cc TEST121 demonstrates the use of PSAT_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_data
      double precision psat
      double precision tc

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST121:'
      write ( *, '(a)' ) '  PSAT_VALUES returns values of '
      write ( *, '(a)' ) '  the saturation pressure of water '
      write ( *, '(a)' ) '  as a function of temperature.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      T            PSAT(T)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call psat_values ( n_data, tc, psat )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g14.6)' ) tc, psat

      go to 10

20    continue

      return
      end
      subroutine test122 ( )

c*********************************************************************72
c
cc TEST122 demonstrates the use of PSI_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST122:'
      write ( *, '(a)' ) '  PSI_VALUES returns values of '
      write ( *, '(a)' ) '  the Psi function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '          X               PSI(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call psi_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f12.8,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test123 ( )

c*********************************************************************72
c
cc TEST123 demonstrates the use of R8_FACTORIAL_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fn
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST123:'
      write ( *, '(a)' ) '  R8_FACTORIAL_VALUES returns values of '
      write ( *, '(a)' ) '  the factorial function '
      write ( *, '(a)' ) '  (using real arithmetic).'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        N       Factorial(N)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call r8_factorial_values ( n_data, n, fn )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,g24.16)' ) n, fn

      go to 10

20    continue

      return
      end
      subroutine test124 ( )

c*********************************************************************72
c
cc TEST124 demonstrates the use of R8_FACTORIAL_LOG_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fn
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST124:'
      write ( *, '(a)' ) '  R8_FACTORIAL_LOG_VALUES returns values of '
      write ( *, '(a)' ) '  the logarithm of the factorial function '
      write ( *, '(a)' ) '  (using real arithmetic).'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        N       Log(Factorial(N))'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call r8_factorial_log_values ( n_data, n, fn )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,g24.16)' ) n, fn

      go to 10

20    continue

      return
      end
      subroutine test125 ( )

c*********************************************************************72
c
cc TEST125 demonstrates the use of SECVIR_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_data
      double precision tc
      double precision vir

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST125:'
      write ( *, '(a)' ) '  SECVIR_VALUES returns values of '
      write ( *, '(a)' ) '  the second virial coefficient of water '
      write ( *, '(a)' ) '  as a function of temperature.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      T            VIR(T)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call secvir_values ( n_data, tc, vir )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g14.6)' ) tc, vir

      go to 10

20    continue

      return
      end
      subroutine test1255 ( )

c*********************************************************************72
c
cc TEST1255 demonstrates the use of SHI_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 January 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1255:'
      write ( *, '(a)' ) '  SHI_VALUES returns values of '
      write ( *, '(a)' ) '  the hyperbolic sine integral function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X            SHI(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call shi_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test126 ( )

c*********************************************************************72
c
cc TEST126 demonstrates the use of SI_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST126:'
      write ( *, '(a)' ) '  SI_VALUES returns values of '
      write ( *, '(a)' ) '  the sine integral function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X            SI(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call si_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test127 ( )

c*********************************************************************72
c
cc TEST127 demonstrates the use of SIGMA_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer fn
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST127:'
      write ( *, '(a)' ) '  SIGMA_VALUES returns values of '
      write ( *, '(a)' ) '  the SIGMA function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N         SIGMA(N)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call sigma_values ( n_data, n, fn )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,i12)' ) n, fn

      go to 10

20    continue

      return
      end
      subroutine test128 ( )

c*********************************************************************72
c
cc TEST128 demonstrates the use of SIN_POWER_INT_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision b
      double precision fx
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST128:'
      write ( *, '(a)' ) '  SIN_POWER_INT returns values of '
      write ( *, '(a)' ) '  the integral of SIN(X)^N from A to B.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) 
     &  '          A               B              N',
     &  '    SIN_POWER_INT(A,B,N)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call sin_power_int_values ( n_data, a, b, n, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,f14.6,2x,i8,2x,g14.6)' ) a, b, n, fx

      go to 10

20    continue

      return
      end
      subroutine test1285 ( )

c*********************************************************************72
c
cc TEST1285 demonstrates the use of SIX_J_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 February 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      double precision j1
      double precision j2
      double precision j3
      double precision j4
      double precision j5
      double precision j6
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1285:'
      write ( *, '(a)' ) '  SIX_J_VALUES returns values of '
      write ( *, '(a)' ) '  the Wigner 6J coefficient.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '      J1      J2      J3      J4      J5      J6        SIX_J'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call six_j_values ( n_data, j1, j2, j3, j4, j5, j6, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, 
     &  '(2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,g24.16)' ) 
     &  j1, j2, j3, j4, j5, j6, fx

      go to 10

20    continue

      return
      end
      subroutine test129 ( )

c*********************************************************************72
c
cc TEST129 demonstrates the use of SOUND_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision c
      integer n_data
      double precision p
      double precision tc

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST129:'
      write ( *, '(a)' ) '  SOUND_VALUES returns values of '
      write ( *, '(a)' ) '  the spead of sound in water '
      write ( *, '(a)' ) '  as a function of temperature and pressure.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      T            P            C(T,P)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call sound_values ( n_data, tc, p, c )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,f14.6,2x,g14.6)' ) tc, p, c

      go to 10

20    continue

      return
      end
      subroutine test131

c*********************************************************************72
c
cc TEST131 demonstrates the use of SPHERE_UNIT_AREA_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 December 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      integer n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST131:'
      write ( *, '(a)' ) '  SPHERE_UNIT_AREA_VALUES returns values of '
      write ( *, '(a)' ) '  the area of the unit sphere'
      write ( *, '(a)' ) '  in various dimensions.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      N            Area'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call sphere_unit_area_values ( n_data, n, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,g24.16)' ) n, fx

      go to 10

20    continue

      return
      end
      subroutine test132

c*********************************************************************72
c
cc TEST132 demonstrates the use of SPHERE_UNIT_VOLUME_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 December 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      integer n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST132:'
      write ( *, '(a)' ) '  SPHERE_UNIT_VOLUME_VALUES returns values'
      write ( *, '(a)' ) '  of the volume of the unit sphere '
      write ( *, '(a)' ) '  in various dimensions.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      N            Volume'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call sphere_unit_volume_values ( n_data, n, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i4,2x,g24.16)' ) n, fx

      go to 10

20    continue

      return
      end
      subroutine test1325

c*********************************************************************72
c
cc TEST1325 demonstrates the use of SPHERICAL_HARMONIC_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer l
      integer m
      integer n_data
      double precision phi
      double precision theta
      double precision yi
      double precision yr

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1325:'
      write ( *, '(a)' ) '  SPHERICAL_HARMONIC_VALUES returns values'
      write ( *, '(a)' ) '  of the spherical harmonic functions.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) 
     &  '   L   M    THETA       PHI       Yr',
     &  '                           Yi'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call spherical_harmonic_values ( n_data, l, m, theta, 
     &    phi, yr, yi )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, 
     &    '(2x,i2,2x,i2,2x,f8.4,2x,f8.4,2x,g24.16,2x,g24.16)' ) 
     &    l, m, theta, phi, yr, yi

      go to 10

20    continue

      return
      end
      subroutine test130

c*********************************************************************72
c
cc TEST130 demonstrates the use of SQRT_VALUES.
c
c  Discussion:
c
c    In this example, we suggest how the tabulated values can be
c    used to look for large discrepancies.  (In fact, in this case,
c    we suppose there will be nonec).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx
      double precision fx2
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST130:'
      write ( *, '(a)' ) '  SQRT evaluates the square root function.'
      write ( *, '(a)' ) '  SQRT_VALUES returns some exact values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '         X      Exact F       SQRT(X)        Diff'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call sqrt_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = sqrt ( x )

        diff = abs ( fx - fx2 )

        write ( *, '(2x,f14.4,2x,g14.6,2x,g14.6,2x,g14.6)' )
     &    x, fx, fx2, diff

      go to 10

20    continue

      return
      end
      subroutine test133

c*********************************************************************72
c
cc TEST133 demonstrates the use of STIRLING1_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer n_data
      integer m
      integer s1

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST133:'
      write ( *, '(a)' ) '  STIRLING1_VALUES returns values of '
      write ( *, '(a)' ) '  the Stirling numbers of the first kind.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         N         M        S1(N,M)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call stirling1_values ( n_data, n, m, s1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,i8,2x,i12)' ) n, m, s1

      go to 10

20    continue

      return
      end
      subroutine test134

c*********************************************************************72
c
cc TEST134 demonstrates the use of STIRLING2_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer n_data
      integer m
      integer s2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST134:'
      write ( *, '(a)' ) '  STIRLING2_VALUES returns values of '
      write ( *, '(a)' ) '  the Stirling numbers of the first kind.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N       M        S2(N,M)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call stirling2_values ( n_data, n, m, s2 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,i8,2x,i12)' ) n, m, s2

      go to 10

20    continue

      return
      end
      subroutine test135

c*********************************************************************72
c
cc TEST135 demonstrates the use of STROMGEN_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST135:'
      write ( *, '(a)' ) '  STROMGEN_VALUES returns values of'
      write ( *, '(a)' ) '  the Stromgen function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X          F(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call stromgen_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test136

c*********************************************************************72
c
cc TEST136 demonstrates the use of STRUVE_H0_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST136:'
      write ( *, '(a)' ) '  STRUVE_H0_VALUES returns values of '
      write ( *, '(a)' ) '  the Struve H0 function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X            H0(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call struve_h0_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g14.6)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test137

c*********************************************************************72
c
cc TEST137 demonstrates the use of STRUVE_H1_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST137:'
      write ( *, '(a)' ) '  STRUVE_H1_VALUES returns values of '
      write ( *, '(a)' ) '  the Struve H1 function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X            H1(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call struve_h1_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g14.6)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test138

c*********************************************************************72
c
cc TEST138 demonstrates the use of STRUVE_L0_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST138:'
      write ( *, '(a)' ) '  STRUVE_L0_VALUES returns values of '
      write ( *, '(a)' ) '  the Struve L0 function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X            L0(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call struve_l0_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g14.6)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test139

c*********************************************************************72
c
cc TEST139 demonstrates the use of STRUVE_L1_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST139:'
      write ( *, '(a)' ) '  STRUVE_L1_VALUES returns values of '
      write ( *, '(a)' ) '  the Struve L1 function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X            L1(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call struve_l1_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g14.6)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test140

c*********************************************************************72
c
cc TEST140 demonstrates the use of STUDENT_CDF_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 January 2006
c
      implicit none

      double precision c
      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST140:'
      write ( *, '(a)' ) '  STUDENT_CDF_VALUES returns values of '
      write ( *, '(a)' ) '  the Student T Cumulative Density Function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      C     X       CDF'
      write ( *, '(a)' ) ' '

      n_data = 0

 10   continue

      call student_cdf_values ( n_data, c, x, fx )

      if ( n_data .eq. 0 ) then
        go to 20
      end if

      write ( *, '(2x,f10.4,2x,f10.4,2x,g24.16)' ) c, x, fx

      go to 10

20    continue

      return
      end
      subroutine test141

c*********************************************************************72
c
cc TEST141 demonstrates the use of STUDENT_NONCENTRAL_CDF_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer df
      double precision fx
      double precision lambda
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST141:'
      write ( *, '(a)' ) '  STUDENT_NONCENTRAL_CDF_VALUES returns'
      write ( *, '(a)' ) '  values of the noncentral Student T '
      write ( *, '(a)' ) '  Cumulative Density Function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X      LAMBDA       DF     CDF'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call student_noncentral_cdf_values ( n_data, df, lambda, 
     &    x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f10.4,2x,f10.4,2x,i8,2x,g14.6)' ) 
     &    x, lambda, df, fx

      go to 10

20    continue

      return
      end
      subroutine test1415

c*********************************************************************72
c
cc TEST1415 demonstrates the use of SUBFACTORIAL_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer fn
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1415:'
      write ( *, '(a)' ) '  SUBFACTORIAL_VALUES returns values of '
      write ( *, '(a)' ) '  the subfactorial function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N     Subfactorial[N]'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call subfactorial_values ( n_data, n, fn )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,i12)' ) n, fn

      go to 10

20    continue

      return
      end
      subroutine test142

c*********************************************************************72
c
cc TEST142 demonstrates the use of SURTEN_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_data
      double precision sigma
      double precision tc

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST142:'
      write ( *, '(a)' ) '  SURTEN_VALUES returns values of '
      write ( *, '(a)' ) '  the surface tension of water '
      write ( *, '(a)' ) '  as a function of temperature.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      T            SIGMA(T)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call surten_values ( n_data, tc, sigma )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g14.6)' ) tc, sigma

      go to 10

20    continue

      return
      end
      subroutine test143

c*********************************************************************72
c
cc TEST143 demonstrates the use of SYNCH1_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST143:'
      write ( *, '(a)' ) '  SYNCH1_VALUES returns values of '
      write ( *, '(a)' ) '  the synchrotron radiation function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X            S1(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call synch1_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g14.6)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test144

c*********************************************************************72
c
cc TEST144 demonstrates the use of SYNCH2_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST144:'
      write ( *, '(a)' ) '  SYNCH2_VALUES returns values of '
      write ( *, '(a)' ) '  the synchrotron radiation function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X            S2(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call synch2_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g14.6)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test145

c*********************************************************************72
c
cc TEST145 demonstrates the use of TAU_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer fn
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST145:'
      write ( *, '(a)' ) '  TAU_VALUES returns values of '
      write ( *, '(a)' ) '  the TAU function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N         TAU(N)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call tau_values ( n_data, n, fn )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,i12)' ) n, fn

      go to 10

20    continue

      return
      end
      subroutine test146

c*********************************************************************72
c
cc TEST146 demonstrates the use of THERCON_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision lambda
      integer n_data
      double precision p
      double precision tc

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST146:'
      write ( *, '(a)' ) '  THERCON_VALUES returns values of '
      write ( *, '(a)' ) '  the thermal conductivity of water '
      write ( *, '(a)' ) '  as a function of temperature and pressure.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      T            P            LAMBDA(T,P)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call thercon_values ( n_data, tc, p, lambda )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,f14.6,2x,g14.6)' ) tc, p, lambda

      go to 10

20    continue

      return
      end
      subroutine test1465

c*********************************************************************72
c
cc TEST1465 demonstrates the use of THREE_J_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 February 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      double precision j1
      double precision j2
      double precision j3
      double precision m1
      double precision m2
      double precision m3
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1465:'
      write ( *, '(a)' ) '  THREE_J_VALUES returns values of '
      write ( *, '(a)' ) '  the Wigner 3J coefficient.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) 
     &  '      J1      J2      J3      ',
     &  'M1      M2      M3        THREE_J'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call three_j_values ( n_data, j1, j2, j3, m1, m2, m3, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, 
     &  '(2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,g24.16)' ) 
     &    j1, j2, j3, m1, m2, m3, fx

      go to 10

20    continue

      return
      end
      subroutine test147

c*********************************************************************72
c
cc TEST147 demonstrates the use of TRAN02_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST147:'
      write ( *, '(a)' ) '  TRAN02_VALUES returns values of '
      write ( *, '(a)' ) '  the order 2 transportation function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X            T2(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call tran02_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g14.6)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test148

c*********************************************************************72
c
cc TEST148 demonstrates the use of TRAN03_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST148:'
      write ( *, '(a)' ) '  TRAN03_VALUES returns values of '
      write ( *, '(a)' ) '  the order 3 transportation function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X            T3(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call tran03_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g14.6)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test149

c*********************************************************************72
c
cc TEST149 demonstrates the use of TRAN04_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST149:'
      write ( *, '(a)' ) '  TRAN04_VALUES returns values of '
      write ( *, '(a)' ) '  the order 4 transportation function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X            T4(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call tran04_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g14.6)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test150

c*********************************************************************72
c
cc TEST150 demonstrates the use of TRAN05_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST150:'
      write ( *, '(a)' ) '  TRAN05_VALUES returns values of '
      write ( *, '(a)' ) '  the order 5 transportation function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X            T5(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call tran05_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g14.6)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test151

c*********************************************************************72
c
cc TEST151 demonstrates the use of TRAN06_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST151:'
      write ( *, '(a)' ) '  TRAN06_VALUES returns values of '
      write ( *, '(a)' ) '  the order 6 transportation function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X            T6(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call tran06_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g14.6)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test152

c*********************************************************************72
c
cc TEST152 demonstrates the use of TRAN07_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST152:'
      write ( *, '(a)' ) '  TRAN07_VALUES returns values of '
      write ( *, '(a)' ) '  the order 7 transportation function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X            T7(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call tran07_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g14.6)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test153

c*********************************************************************72
c
cc TEST153 demonstrates the use of TRAN08_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST153:'
      write ( *, '(a)' ) '  TRAN08_VALUES returns values of '
      write ( *, '(a)' ) '  the order 8 transportation function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X            T8(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call tran08_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g14.6)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test154

c*********************************************************************72
c
cc TEST154 demonstrates the use of TRAN09_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST154:'
      write ( *, '(a)' ) '  TRAN09_VALUES returns values of '
      write ( *, '(a)' ) '  the order 9 transportation function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X            T9(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call tran09_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g14.6)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test1545

c*********************************************************************72
c
cc TEST1545 demonstrates the use of TRIGAMMA_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1545:'
      write ( *, '(a)' ) '  TRIGAMMA_VALUES returns values of '
      write ( *, '(a)' ) '  the TriGamma function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X            F(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call trigamma_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

      go to 10

20    continue

      return
      end
      subroutine test155

c*********************************************************************72
c
cc TEST155 demonstrates the use of TSAT_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_data
      double precision p
      double precision tc

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST155:'
      write ( *, '(a)' ) '  TSAT_VALUES returns values of '
      write ( *, '(a)' ) '  the saturation temperature '
      write ( *, '(a)' ) '  as a function of pressure.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      P           Tsat(P)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call tsat_values ( n_data, p, tc )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,g14.6)' ) p, tc

      go to 10

20    continue

      return
      end
      subroutine test156

c*********************************************************************72
c
cc TEST156 demonstrates the use of VAN_DER_CORPUT_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer base
      integer n_data
      integer seed
      double precision value

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST156:'
      write ( *, '(a)' ) '  VAN_DER_CORPUT_VALUES returns values of '
      write ( *, '(a)' ) '  the van der Corput sequence '
      write ( *, '(a)' ) '  in a given base.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      BASE      SEED    VDC(BASE,SEED)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call van_der_corput_values ( n_data, base, seed, value )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,i8,2x,g16.8)' ) base, seed, value

      go to 10

20    continue

      return
      end
      subroutine test157

c*********************************************************************72
c
cc TEST157 demonstrates the use of VISCOSITY_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision eta
      integer n_data
      double precision p
      double precision tc

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST157:'
      write ( *, '(a)' ) '  VISCOSITY_VALUES returns values of '
      write ( *, '(a)' ) '  the viscosity of water '
      write ( *, '(a)' ) '  as a function of temperature and pressure.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      T            P            ETA(T,P)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call viscosity_values ( n_data, tc, p, eta )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,f14.6,2x,g14.6)' ) tc, p, eta

      go to 10

20    continue

      return
      end
      subroutine test1575

c*********************************************************************72
c
cc TEST1575 demonstrates the use of VON_MISES_CDF_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision b
      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1575:'
      write ( *, '(a)' ) '  VON_MISES_CDF_VALUES returns values of '
      write ( *, '(a)' ) '  the von Mises CDF.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '      A            B            X            CDF'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call von_mises_cdf_values ( n_data, a, b, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,f14.6,2x,f14.6,2x,f14.6,2x,g14.6)' ) 
     &    a, b, x, fx

      go to 10

20    continue

      return
      end
      subroutine test158

c*********************************************************************72
c
cc TEST158 demonstrates the use of WEIBULL_CDF_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision alpha
      double precision beta
      double precision fx
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST158:'
      write ( *, '(a)' ) '  WEIBULL_CDF_VALUES returns values of '
      write ( *, '(a)' ) '  the Weibull Cumulative Density Function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '      ALPHA   BETA    X                  CDF'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call weibull_cdf_values ( n_data, alpha, beta, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g24.16)' ) 
     &    alpha, beta, x, fx

      go to 10

20    continue

      return
      end
      subroutine test159

c*********************************************************************72
c
cc TEST159 demonstrates the use of ZETA_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer n_data
      double precision zeta

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST159:'
      write ( *, '(a)' ) '  ZETA_VALUES returns values of '
      write ( *, '(a)' ) '  the Riemann Zeta function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         N        ZETA(N)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call zeta_values ( n_data, n, zeta )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,i8,2x,g24.16)' ) n, zeta

      go to 10

20    continue

      return
      end

