      program main

c*********************************************************************72
c
cc MAIN is the main program for G77_INTRINSICS.
c
c  Discussion:
c
c    G77_INTRINSICS calls some of the G77 intrinsic routines.
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

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'G77_INTRINSICS'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the G77 intrinsic library.'

      call test_besj0
      call test_besj1
      call test_besjn
      call test_besy0
      call test_besy1
      call test_besyn
      call test_dbesj0
      call test_dbesj1
      call test_dbesjn
      call test_dbesy0
      call test_dbesy1
      call test_dbesyn
      call test_derf
      call test_derfc
      call test_erf
      call test_erfc
      call test_rand
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'G77_INTRINSICS'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test_besj0

c*********************************************************************72
c
cc TEST_BESJ0 checks BESJ0 against BESSEL_J0_VALUES.
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
      real fx2
      integer n_data
      double precision x
      real x2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_BESJ0:'
      write ( *, '(a)' ) '  BESJ0 computes the Bessel J0 function.'
      write ( *, '(a)' ) '  BESSEL_J0_VALUES returns selected values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '      X             FX',
     &  '                        FX2'
      write ( *, '(a,a)' )
     &  '                   (table)',
     &  '                   (BESJ0)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_j0_values ( n_data, x, fx )

        if ( n_data == 0 ) then
          go to 20
        end if

        x2 = real ( x )
        fx2 = besj0 ( x2 )

        write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' )
     &    x, fx, fx2

      go to 10

20    continue

      return
      end
      subroutine test_besj1

c*********************************************************************72
c
cc TEST_BESJ1 checks BESJ1 against BESSEL_J1_VALUES.
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
      real fx2
      integer n_data
      double precision x
      real x2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_BESJ1:'
      write ( *, '(a)' ) '  BESJ1 computes the Bessel J1 function.'
      write ( *, '(a)' ) '  BESSEL_J1_VALUES returns selected values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '      X             FX',
     &  '                        FX2'
      write ( *, '(a,a)' )
     &  '                   (table)',
     &  '                   (BESJ1)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_j1_values ( n_data, x, fx )

        if ( n_data == 0 ) then
          go to 20
        end if

        x2 = real ( x )
        fx2 = besj1 ( x2 )

        write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' )
     &    x, fx, fx2

      go to 10

20    continue

      return
      end
      subroutine test_besjn

c*********************************************************************72
c
cc TEST_BESJN checks BESJN against BESSEL_JN_VALUES.
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
      real fx2
      integer n
      integer n_data
      double precision x
      real x2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_BESJN:'
      write ( *, '(a)' ) '  BESJN computes the Bessel Jn function.'
      write ( *, '(a)' ) '  BESSEL_JN_VALUES returns selected values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '      N     X             FX',
     &  '                        FX2'
      write ( *, '(a,a)' )
     &  '                         (table)',
     &  '                   (BESJN)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_jn_values ( n_data, n, x, fx )

        if ( n_data == 0 ) then
          go to 20
        end if

        x2 = real ( x )
        fx2 = besjn ( n, x2 )

        write ( *, '(2x,i4,2x,f12.8,2x,g24.16,2x,g24.16)' )
     &    n, x, fx, fx2

      go to 10

20    continue

      return
      end
      subroutine test_besy0

c*********************************************************************72
c
cc TEST_BESY0 checks BESY0 against BESSEL_Y0_VALUES.
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
      real fx2
      integer n_data
      double precision x
      real x2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_BESY0:'
      write ( *, '(a)' ) '  BESY0 computes the Bessel Y0 function.'
      write ( *, '(a)' ) '  BESSEL_Y0_VALUES returns selected values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '      X             FX',
     &  '                        FX2'
      write ( *, '(a,a)' )
     &  '                   (table)',
     &  '                   (BESY0)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_y0_values ( n_data, x, fx )

        if ( n_data == 0 ) then
          go to 20
        end if

        x2 = real ( x )
        fx2 = besy0 ( x2 )

        write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' )
     &    x, fx, fx2

      go to 10

20    continue

      return
      end
      subroutine test_besy1

c*********************************************************************72
c
cc TEST_BESY1 checks BESY1 against BESSEL_Y1_VALUES.
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
      real fx2
      integer n_data
      double precision x
      real x2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_BESY1:'
      write ( *, '(a)' ) '  BESY1 computes the Bessel Y1 function.'
      write ( *, '(a)' ) '  BESSEL_Y1_VALUES returns selected values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '      X             FX',
     &  '                        FX2'
      write ( *, '(a,a)' )
     &  '                   (table)',
     &  '                   (BESY1)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_y1_values ( n_data, x, fx )

        if ( n_data == 0 ) then
          go to 20
        end if

        x2 = real ( x )
        fx2 = besy1 ( x2 )

        write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' )
     &    x, fx, fx2

      go to 10

20    continue

      return
      end
      subroutine test_besyn

c*********************************************************************72
c
cc TEST_BESYN checks BESYN against BESSEL_YN_VALUES.
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
      real fx2
      integer n
      integer n_data
      double precision x
      real x2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_BESYN:'
      write ( *, '(a)' ) '  BESYN computes the Bessel Yn function.'
      write ( *, '(a)' ) '  BESSEL_YN_VALUES returns selected values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '      N     X             FX',
     &  '                        FX2'
      write ( *, '(a,a)' )
     &  '                         (table)',
     &  '                   (BESYN)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_yn_values ( n_data, n, x, fx )

        if ( n_data == 0 ) then
          go to 20
        end if

        x2 = real ( x )
        fx2 = besyn ( n, x2 )

        write ( *, '(2x,i4,2x,f12.8,2x,g24.16,2x,g24.16)' )
     &    n, x, fx, fx2

      go to 10

20    continue

      return
      end
      subroutine test_dbesj0

c*********************************************************************72
c
cc TEST_DBESJ0 checks DBESJ0 against BESSEL_J0_VALUES.
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
      double precision fx2
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_DBESJ0:'
      write ( *, '(a)' ) '  DBESJ0 computes the Bessel J0 function.'
      write ( *, '(a)' ) '  BESSEL_J0_VALUES returns selected values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '      X             FX',
     &  '                        FX2'
      write ( *, '(a,a)' )
     &  '                   (table)',
     &  '                   (DBESJ0)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_j0_values ( n_data, x, fx )

        if ( n_data == 0 ) then
          go to 20
        end if

        fx2 = dbesj0 ( x )

        write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' )
     &    x, fx, fx2

      go to 10

20    continue

      return
      end
      subroutine test_dbesj1

c*********************************************************************72
c
cc TEST_DBESJ1 checks DBESJ1 against BESSEL_J1_VALUES.
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
      double precision fx2
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_DBESJ1:'
      write ( *, '(a)' ) '  DBESJ1 computes the Bessel J1 function.'
      write ( *, '(a)' ) '  BESSEL_J1_VALUES returns selected values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '      X             FX',
     &  '                        FX2'
      write ( *, '(a,a)' )
     &  '                   (table)',
     &  '                   (DBESJ1)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_j1_values ( n_data, x, fx )

        if ( n_data == 0 ) then
          go to 20
        end if

        fx2 = dbesj1 ( x )

        write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' )
     &    x, fx, fx2

      go to 10

20    continue

      return
      end
      subroutine test_dbesjn

c*********************************************************************72
c
cc TEST_DBESJN checks DBESJN against BESSEL_JN_VALUES.
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
      double precision fx2
      integer n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_DBESJN:'
      write ( *, '(a)' ) '  DBESJN computes the Bessel Jn function.'
      write ( *, '(a)' ) '  BESSEL_JN_VALUES returns selected values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '      N     X             FX',
     &  '                        FX2'
      write ( *, '(a,a)' )
     &  '                         (table)',
     &  '                   (DBESJN)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_jn_values ( n_data, n, x, fx )

        if ( n_data == 0 ) then
          go to 20
        end if

        fx2 = dbesjn ( n, x )

        write ( *, '(2x,i4,2x,f12.8,2x,g24.16,2x,g24.16)' )
     &    n, x, fx, fx2

      go to 10

20    continue

      return
      end
      subroutine test_dbesy0

c*********************************************************************72
c
cc TEST_DBESY0 checks DBESY0 against BESSEL_Y0_VALUES.
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
      double precision fx2
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_DBESY0:'
      write ( *, '(a)' ) '  DBESY0 computes the Bessel Y0 function.'
      write ( *, '(a)' ) '  BESSEL_Y0_VALUES returns selected values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '      X             FX',
     &  '                        FX2'
      write ( *, '(a,a)' )
     &  '                   (table)',
     &  '                   (DBESY0)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_y0_values ( n_data, x, fx )

        if ( n_data == 0 ) then
          go to 20
        end if

        fx2 = dbesy0 ( x )

        write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' )
     &    x, fx, fx2

      go to 10

20    continue

      return
      end
      subroutine test_dbesy1

c*********************************************************************72
c
cc TEST_DBESY1 checks DBESY1 against BESSEL_Y1_VALUES.
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
      double precision fx2
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_DBESY1:'
      write ( *, '(a)' ) '  DBESY1 computes the Bessel Y1 function.'
      write ( *, '(a)' ) '  BESSEL_Y1_VALUES returns selected values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '      X             FX',
     &  '                        FX2'
      write ( *, '(a,a)' )
     &  '                   (table)',
     &  '                   (DBESY1)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_y1_values ( n_data, x, fx )

        if ( n_data == 0 ) then
          go to 20
        end if

        fx2 = dbesy1 ( x )

        write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' )
     &    x, fx, fx2

      go to 10

20    continue

      return
      end
      subroutine test_dbesyn

c*********************************************************************72
c
cc TEST_DBESYN checks DBESYN against BESSEL_YN_VALUES.
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
      double precision fx2
      integer n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_DBESYN:'
      write ( *, '(a)' ) '  DBESYN computes the Bessel Yn function.'
      write ( *, '(a)' ) '  BESSEL_YN_VALUES returns selected values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '      N     X             FX',
     &  '                        FX2'
      write ( *, '(a,a)' )
     &  '                         (table)',
     &  '                   (DBESYN)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bessel_yn_values ( n_data, n, x, fx )

        if ( n_data == 0 ) then
          go to 20
        end if

        fx2 = dbesyn ( n, x )

        write ( *, '(2x,i4,2x,f12.8,2x,g24.16,2x,g24.16)' )
     &    n, x, fx, fx2

      go to 10

20    continue

      return
      end
      subroutine test_derf

c*********************************************************************72
c
cc TEST_DERF checks DERF against ERF_VALUES.
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
      double precision fx2
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_DERF:'
      write ( *, '(a)' ) '  DERF computes the error function.'
      write ( *, '(a)' ) '  ERF_VALUES returns selected values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '      X             FX',
     &  '                        FX2'
      write ( *, '(a,a)' )
     &  '                   (table)',
     &  '                   (DERF)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call erf_values ( n_data, x, fx )

        if ( n_data == 0 ) then
          go to 20
        end if

        fx2 = derf ( x )

        write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' )
     &    x, fx, fx2

      go to 10

20    continue

      return
      end
      subroutine test_derfc

c*********************************************************************72
c
cc TEST_DERFC checks DERFC against ERFC_VALUES.
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
      double precision fx2
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_DERFC:'
      write ( *, '(a)' )
     &  '  DERFC computes the complementary error function.'
      write ( *, '(a)' ) '  ERFC_VALUES returns selected values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '      X             FX',
     &  '                        FX2'
      write ( *, '(a,a)' )
     &  '                   (table)',
     &  '                   (DERFC)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call erfc_values ( n_data, x, fx )

        if ( n_data == 0 ) then
          go to 20
        end if

        fx2 = derfc ( x )

        write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' )
     &    x, fx, fx2

      go to 10

20    continue

      return
      end
      subroutine test_erf

c*********************************************************************72
c
cc TEST_ERF checks ERF against ERF_VALUES.
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
      real fx2
      integer n_data
      double precision x
      real x2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_ERF:'
      write ( *, '(a)' ) '  ERF computes the error function.'
      write ( *, '(a)' ) '  ERF_VALUES returns selected values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '      X             FX',
     &  '                        FX2'
      write ( *, '(a,a)' )
     &  '                   (table)',
     &  '                   (ERF)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call erf_values ( n_data, x, fx )

        if ( n_data == 0 ) then
          go to 20
        end if

        x2 = real ( x )
        fx2 = erf ( x2 )

        write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' )
     &    x, fx, fx2

      go to 10

20    continue

      return
      end
      subroutine test_erfc

c*********************************************************************72
c
cc TEST_ERFC checks ERFC against ERFC_VALUES.
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
      real fx2
      integer n_data
      double precision x
      real x2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_ERFC:'
      write ( *, '(a)' )
     &  '  ERFC computes the complementary error function.'
      write ( *, '(a)' ) '  ERFC_VALUES returns selected values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '      X             FX',
     &  '                        FX2'
      write ( *, '(a,a)' )
     &  '                   (table)',
     &  '                   (ERFC)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call erfc_values ( n_data, x, fx )

        if ( n_data == 0 ) then
          go to 20
        end if

        x2 = real ( x )
        fx2 = erfc ( x2 )

        write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' )
     &    x, fx, fx2

      go to 10

20    continue

      return
      end
      subroutine test_rand

c*********************************************************************72
c
cc TEST_RAND tests RAND.
c
c  Discussion:
c
c    The seed is optional.
c
c    If it is 0, it is ignored and the next value in the random sequence
c    is returned.
c
c    If it is 1, the generator is restarted by calling SRAND(0).
c
c    Otherwise, it is used as input to SRAND to reseed the sequence.
c
c    Normal usage of RAND would be to call once with the seed nonzero,
c    and then to call repeatedly with the seed 0 or not specified.
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

      integer i
      real r
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_RAND'
      write ( *, '(a)' ) '  RAND returns a real random value.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The sequence is "seeded" by:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    R = RAND ( 1 ): reseeds by SRAND ( 0 )'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    or'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    R = RAND ( S ): reseeds by SRAND ( S )'
      write ( *, '(a)' ) '                    assuming S not 0, not 1.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The sequence is used by:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    R = RAND ( ): returns next value in'
      write ( *, '(a)' ) '                  current random sequence.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    or'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    R = RAND ( 0 ): same as R = RAND ( );'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    or'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    R = RAND ( S ): reseeds sequence, and'
      write ( *, '(a)' ) '                    returns first value in'
      write ( *, '(a)' ) '                    new sequence,'
      write ( *, '(a)' ) '                    assuming S is not 0.'

      seed = 123456789
      call srand ( seed )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Call SRAND(123456789)'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  do i = 1, 10'
      write ( *, '(a)' ) '    R = RAND ( )'
      write ( *, '(a)' ) '  end do'
      write ( *, '(a)' ) ' '

      do i = 1, 10
        r = rand ( )
        write ( *, '(2x,i8,2x,g14.6)' ) i, r
      end do

      seed = 987654321
      call srand ( seed )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  CHANGING THE SEED CHANGES THE SEQUENCE:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Call SRAND(987654321)'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  do i = 1, 10'
      write ( *, '(a)' ) '    R = RAND ( )'
      write ( *, '(a)' ) '  end do'
      write ( *, '(a)' ) ' '

      do i = 1, 10
        r = rand ( )
        write ( *, '(2x,i8,2x,g14.6)' ) i, r
      end do

      seed = 123456789
      call srand ( seed )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  RESTORING THE OLD SEED RESTARTS THE OLD SEQUENCE:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Call SRAND(123456789)'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  do i = 1, 10'
      write ( *, '(a)' ) '    R = RAND ( )'
      write ( *, '(a)' ) '  end do'
      write ( *, '(a)' ) ' '

      do i = 1, 10
        r = rand ( )
        write ( *, '(2x,i8,2x,g14.6)' ) i, r
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  WE CAN PASS THE SEED IN THROUGH RAND:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  R = RAND(987654321)'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  do i = 2, 10'
      write ( *, '(a)' ) '    R = RAND ( )'
      write ( *, '(a)' ) '  end do'
      write ( *, '(a)' ) ' '

      r = rand ( 987654321 )

      i = 1
      write ( *, '(2x,i8,2x,g14.6)' ) i, r

      do i = 2, 10
        r = rand ( )
        write ( *, '(2x,i8,2x,g14.6)' ) i, r
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  WE CAN GET A "RANDOM" SEED:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  R = RAND(1)'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  do i = 2, 10'
      write ( *, '(a)' ) '    R = RAND ( )'
      write ( *, '(a)' ) '  end do'
      write ( *, '(a)' ) ' '

      r = rand ( 1 )

      i = 1
      write ( *, '(2x,i8,2x,g14.6)' ) i, r

      do i = 2, 10
        r = rand ( )
        write ( *, '(2x,i8,2x,g14.6)' ) i, r
      end do

      seed = 123456789
      call srand ( seed )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  CALLING WITH RAND(0) IS THE SAME AS RAND():'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Call SRAND(123456789)'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  do i = 1, 10'
      write ( *, '(a)' ) '    R = RAND ( )'
      write ( *, '(a)' ) '  end do'
      write ( *, '(a)' ) ' '

      do i = 1, 10
        r = rand ( 0 )
        write ( *, '(2x,i8,2x,g14.6)' ) i, r
      end do

      return
      end
      subroutine bessel_j0_values ( n_data, x, fx )

c*********************************************************************72
c
cc BESSEL_J0_VALUES returns some values of the J0 Bessel function.
c
c  Discussion:
c
c    In Mathematica, the function can be evaluated by:
c
c      BesselJ[0,x]
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
c  Reference:
c
c    Milton Abramowitz, Irene Stegun,
c    Handbook of Mathematical Functions,
c    National Bureau of Standards, 1964,
c    ISBN: 0-486-61272-4,
c    LC: QA47.A34.
c
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Cambridge University Press, 1999,
c    ISBN: 0-521-64314-7,
c    LC: QA76.95.W65.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, double precision X, the argument of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 21 )

      double precision fx
      double precision fx_vec(n_max)
      integer n_data
      double precision x
      double precision x_vec(n_max)

      save fx_vec
      save x_vec

      data fx_vec /
     &  -0.1775967713143383D+00,
     &  -0.3971498098638474D+00,
     &  -0.2600519549019334D+00,
     &   0.2238907791412357D+00,
     &   0.7651976865579666D+00,
     &   0.1000000000000000D+01,
     &   0.7651976865579666D+00,
     &   0.2238907791412357D+00,
     &  -0.2600519549019334D+00,
     &  -0.3971498098638474D+00,
     &  -0.1775967713143383D+00,
     &   0.1506452572509969D+00,
     &   0.3000792705195556D+00,
     &   0.1716508071375539D+00,
     &  -0.9033361118287613D-01,
     &  -0.2459357644513483D+00,
     &  -0.1711903004071961D+00,
     &   0.4768931079683354D-01,
     &   0.2069261023770678D+00,
     &   0.1710734761104587D+00,
     &  -0.1422447282678077D-01 /
      data x_vec /
     &  -5.0D+00,
     &  -4.0D+00,
     &  -3.0D+00,
     &  -2.0D+00,
     &  -1.0D+00,
     &   0.0D+00,
     &   1.0D+00,
     &   2.0D+00,
     &   3.0D+00,
     &   4.0D+00,
     &   5.0D+00,
     &   6.0D+00,
     &   7.0D+00,
     &   8.0D+00,
     &   9.0D+00,
     &  10.0D+00,
     &  11.0D+00,
     &  12.0D+00,
     &  13.0D+00,
     &  14.0D+00,
     &  15.0D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        x = 0.0D+00
        fx = 0.0D+00
      else
        x = x_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      subroutine bessel_j1_values ( n_data, x, fx )

c*********************************************************************72
c
cc BESSEL_J1_VALUES returns some values of the J1 Bessel function.
c
c  Discussion:
c
c    In Mathematica, the function can be evaluated by:
c
c      BesselJ[1,x]
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
c  Reference:
c
c    Milton Abramowitz, Irene Stegun,
c    Handbook of Mathematical Functions,
c    National Bureau of Standards, 1964,
c    ISBN: 0-486-61272-4,
c    LC: QA47.A34.
c
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Cambridge University Press, 1999,
c    ISBN: 0-521-64314-7,
c    LC: QA76.95.W65.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, double precision X, the argument of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 21 )

      double precision fx
      double precision fx_vec(n_max)
      integer n_data
      double precision x
      double precision x_vec(n_max)

      save fx_vec
      save x_vec

      data fx_vec /
     &   0.3275791375914652D+00,
     &   0.6604332802354914D-01,
     &  -0.3390589585259365D+00,
     &  -0.5767248077568734D+00,
     &  -0.4400505857449335D+00,
     &   0.0000000000000000D+00,
     &   0.4400505857449335D+00,
     &   0.5767248077568734D+00,
     &   0.3390589585259365D+00,
     &  -0.6604332802354914D-01,
     &  -0.3275791375914652D+00,
     &  -0.2766838581275656D+00,
     &  -0.4682823482345833D-02,
     &   0.2346363468539146D+00,
     &   0.2453117865733253D+00,
     &   0.4347274616886144D-01,
     &  -0.1767852989567215D+00,
     &  -0.2234471044906276D+00,
     &  -0.7031805212177837D-01,
     &   0.1333751546987933D+00,
     &   0.2051040386135228D+00 /
      data x_vec /
     &  -5.0D+00,
     &  -4.0D+00,
     &  -3.0D+00,
     &  -2.0D+00,
     &  -1.0D+00,
     &   0.0D+00,
     &   1.0D+00,
     &   2.0D+00,
     &   3.0D+00,
     &   4.0D+00,
     &   5.0D+00,
     &   6.0D+00,
     &   7.0D+00,
     &   8.0D+00,
     &   9.0D+00,
     &  10.0D+00,
     &  11.0D+00,
     &  12.0D+00,
     &  13.0D+00,
     &  14.0D+00,
     &   15.0D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        x = 0.0D+00
        fx = 0.0D+00
      else
        x = x_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      subroutine bessel_jn_values ( n_data, nu, x, fx )

c*********************************************************************72
c
cc BESSEL_JN_VALUES returns some values of the Jn Bessel function.
c
c  Discussion:
c
c    In Mathematica, the function can be evaluated by:
c
c      BesselJ[n,x]
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
c  Reference:
c
c    Milton Abramowitz, Irene Stegun,
c    Handbook of Mathematical Functions,
c    National Bureau of Standards, 1964,
c    ISBN: 0-486-61272-4,
c    LC: QA47.A34.
c
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Cambridge University Press, 1999,
c    ISBN: 0-521-64314-7,
c    LC: QA76.95.W65.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, integer NU, the order of the function.
c
c    Output, double precision X, the argument of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 20 )

      double precision fx
      double precision fx_vec(n_max)
      integer n_data
      integer nu
      integer nu_vec(n_max)
      double precision x
      double precision x_vec(n_max)

      save fx_vec
      save nu_vec
      save x_vec

      data fx_vec /
     &   0.1149034849319005D+00,
     &   0.3528340286156377D+00,
     &   0.4656511627775222D-01,
     &   0.2546303136851206D+00,
     &  -0.5971280079425882D-01,
     &   0.2497577302112344D-03,
     &   0.7039629755871685D-02,
     &   0.2611405461201701D+00,
     &  -0.2340615281867936D+00,
     &  -0.8140024769656964D-01,
     &   0.2630615123687453D-09,
     &   0.2515386282716737D-06,
     &   0.1467802647310474D-02,
     &   0.2074861066333589D+00,
     &  -0.1138478491494694D+00,
     &   0.3873503008524658D-24,
     &   0.3918972805090754D-18,
     &   0.2770330052128942D-10,
     &   0.1151336924781340D-04,
     &  -0.1167043527595797D+00 /
      data nu_vec /
     &   2,  2,  2,  2,
     &   2,  5,  5,  5,
     &   5,  5, 10, 10,
     &  10, 10, 10, 20,
     &  20, 20, 20, 20 /
      data x_vec /
     &   1.0D+00,
     &   2.0D+00,
     &   5.0D+00,
     &  10.0D+00,
     &  50.0D+00,
     &   1.0D+00,
     &   2.0D+00,
     &   5.0D+00,
     &  10.0D+00,
     &  50.0D+00,
     &   1.0D+00,
     &   2.0D+00,
     &   5.0D+00,
     &  10.0D+00,
     &  50.0D+00,
     &   1.0D+00,
     &   2.0D+00,
     &   5.0D+00,
     &  10.0D+00,
     &  50.0D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        nu = 0
        x = 0.0D+00
        fx = 0.0D+00
      else
        nu = nu_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      subroutine bessel_y0_values ( n_data, x, fx )

c*********************************************************************72
c
cc BESSEL_Y0_VALUES returns some values of the Y0 Bessel function.
c
c  Discussion:
c
c    In Mathematica, the function can be evaluated by:
c
c      BesselY[0,x]
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
c  Reference:
c
c    Milton Abramowitz, Irene Stegun,
c    Handbook of Mathematical Functions,
c    National Bureau of Standards, 1964,
c    ISBN: 0-486-61272-4,
c    LC: QA47.A34.
c
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Cambridge University Press, 1999,
c    ISBN: 0-521-64314-7,
c    LC: QA76.95.W65.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, double precision X, the argument of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 16 )

      double precision fx
      double precision fx_vec(n_max)
      integer n_data
      double precision x
      double precision x_vec(n_max)

      save fx_vec
      save x_vec

      data fx_vec /
     &  -0.1534238651350367D+01,
     &   0.8825696421567696D-01,
     &   0.5103756726497451D+00,
     &   0.3768500100127904D+00,
     &  -0.1694073932506499D-01,
     &  -0.3085176252490338D+00,
     &  -0.2881946839815792D+00,
     &  -0.2594974396720926D-01,
     &   0.2235214893875662D+00,
     &   0.2499366982850247D+00,
     &   0.5567116728359939D-01,
     &  -0.1688473238920795D+00,
     &  -0.2252373126343614D+00,
     &  -0.7820786452787591D-01,
     &   0.1271925685821837D+00,
     &   0.2054642960389183D+00 /
      data x_vec /
     &   0.1D+00,
     &   1.0D+00,
     &   2.0D+00,
     &   3.0D+00,
     &   4.0D+00,
     &   5.0D+00,
     &   6.0D+00,
     &   7.0D+00,
     &   8.0D+00,
     &   9.0D+00,
     &  10.0D+00,
     &  11.0D+00,
     &  12.0D+00,
     &  13.0D+00,
     &  14.0D+00,
     &  15.0D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        x = 0.0D+00
        fx = 0.0D+00
      else
        x = x_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      subroutine bessel_y1_values ( n_data, x, fx )

c*********************************************************************72
c
cc BESSEL_Y1_VALUES returns some values of the Y1 Bessel function.
c
c  Discussion:
c
c    In Mathematica, the function can be evaluated by:
c
c      BesselY[1,x]
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
c  Reference:
c
c    Milton Abramowitz, Irene Stegun,
c    Handbook of Mathematical Functions,
c    National Bureau of Standards, 1964,
c    ISBN: 0-486-61272-4,
c    LC: QA47.A34.
c
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Cambridge University Press, 1999,
c    ISBN: 0-521-64314-7,
c    LC: QA76.95.W65.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, double precision X, the argument of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 16 )

      double precision fx
      double precision fx_vec(n_max)
      integer n_data
      double precision x
      double precision x_vec(n_max)

      save fx_vec
      save x_vec

      data fx_vec /
     &  -0.6458951094702027D+01,
     &  -0.7812128213002887D+00,
     &  -0.1070324315409375D+00,
     &   0.3246744247918000D+00,
     &   0.3979257105571000D+00,
     &   0.1478631433912268D+00,
     &  -0.1750103443003983D+00,
     &  -0.3026672370241849D+00,
     &  -0.1580604617312475D+00,
     &   0.1043145751967159D+00,
     &   0.2490154242069539D+00,
     &   0.1637055374149429D+00,
     &  -0.5709921826089652D-01,
     &  -0.2100814084206935D+00,
     &  -0.1666448418561723D+00,
     &   0.2107362803687351D-01 /
      data x_vec /
     &   0.1D+00,
     &   1.0D+00,
     &   2.0D+00,
     &   3.0D+00,
     &   4.0D+00,
     &   5.0D+00,
     &   6.0D+00,
     &   7.0D+00,
     &   8.0D+00,
     &   9.0D+00,
     &  10.0D+00,
     &  11.0D+00,
     &  12.0D+00,
     &  13.0D+00,
     &  14.0D+00,
     &  15.0D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        x = 0.0D+00
        fx = 0.0D+00
      else
        x = x_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      subroutine bessel_yn_values ( n_data, nu, x, fx )

c*********************************************************************72
c
cc BESSEL_YN_VALUES returns some values of the Yn Bessel function.
c
c  Discussion:
c
c    In Mathematica, the function can be evaluated by:
c
c      BesselY[n,x]
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
c  Reference:
c
c    Milton Abramowitz, Irene Stegun,
c    Handbook of Mathematical Functions,
c    National Bureau of Standards, 1964,
c    ISBN: 0-486-61272-4,
c    LC: QA47.A34.
c
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Cambridge University Press, 1999,
c    ISBN: 0-521-64314-7,
c    LC: QA76.95.W65.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, integer NU, the order of the function.
c
c    Output, double precision X, the argument of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 20 )

      double precision fx
      double precision fx_vec(n_max)
      integer n_data
      integer nu
      integer nu_vec(n_max)
      double precision x
      double precision x_vec(n_max)

      save fx_vec
      save nu_vec
      save x_vec

      data fx_vec /
     &  -0.1650682606816254D+01,
     &  -0.6174081041906827D+00,
     &   0.3676628826055245D+00,
     &  -0.5868082442208615D-02,
     &   0.9579316872759649D-01,
     &  -0.2604058666258122D+03,
     &  -0.9935989128481975D+01,
     &  -0.4536948224911019D+00,
     &   0.1354030476893623D+00,
     &  -0.7854841391308165D-01,
     &  -0.1216180142786892D+09,
     &  -0.1291845422080393D+06,
     &  -0.2512911009561010D+02,
     &  -0.3598141521834027D+00,
     &   0.5723897182053514D-02,
     &  -0.4113970314835505D+23,
     &  -0.4081651388998367D+17,
     &  -0.5933965296914321D+09,
     &  -0.1597483848269626D+04,
     &   0.1644263394811578D-01 /
      data nu_vec /
     &   2,  2,  2,  2,
     &   2,  5,  5,  5,
     &   5,  5, 10, 10,
     &  10, 10, 10, 20,
     &  20, 20, 20, 20 /
      data x_vec /
     &   1.0D+00,
     &   2.0D+00,
     &   5.0D+00,
     &  10.0D+00,
     &  50.0D+00,
     &   1.0D+00,
     &   2.0D+00,
     &   5.0D+00,
     &  10.0D+00,
     &  50.0D+00,
     &   1.0D+00,
     &   2.0D+00,
     &   5.0D+00,
     &  10.0D+00,
     &  50.0D+00,
     &   1.0D+00,
     &   2.0D+00,
     &   5.0D+00,
     &  10.0D+00,
     &  50.0D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        nu = 0
        x = 0.0D+00
        fx = 0.0D+00
      else
        nu = nu_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      subroutine erf_values ( n_data, x, fx )

c*********************************************************************72
c
cc ERF_VALUES returns some values of the ERF or "error" function for testing.
c
c  Discussion:
c
c    The error function is defined by:
c
c      ERF(X) = ( 2 / sqrt ( PI ) * integral ( 0 <= T <= X ) exp ( - T^2 ) dT
c
c    In Mathematica, the function can be evaluated by:
c
c      Erf[x]
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 March 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Milton Abramowitz, Irene Stegun,
c    Handbook of Mathematical Functions,
c    National Bureau of Standards, 1964,
c    ISBN: 0-486-61272-4,
c    LC: QA47.A34.
c
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Cambridge University Press, 1999,
c    ISBN: 0-521-64314-7,
c    LC: QA76.95.W65.
c
c  Parameters:
c
c    Input/output, integer N_DATA.
c    On input, if N_DATA is 0, the first test data is returned, and
c    N_DATA is set to the index of the test data.  On each subsequent
c    call, N_DATA is incremented and that test data is returned.  When
c    there is no more test data, N_DATA is set to 0.
c
c    Output, double precision X, the argument of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 21 )

      double precision bvec ( n_max )
      double precision fx
      integer n_data
      double precision x
      double precision xvec ( n_max )

      data bvec /
     & 0.0000000000000000D+00,
     & 0.1124629160182849D+00,
     & 0.2227025892104785D+00,
     & 0.3286267594591274D+00,
     & 0.4283923550466685D+00,
     & 0.5204998778130465D+00,
     & 0.6038560908479259D+00,
     & 0.6778011938374185D+00,
     & 0.7421009647076605D+00,
     & 0.7969082124228321D+00,
     & 0.8427007929497149D+00,
     & 0.8802050695740817D+00,
     & 0.9103139782296354D+00,
     & 0.9340079449406524D+00,
     & 0.9522851197626488D+00,
     & 0.9661051464753107D+00,
     & 0.9763483833446440D+00,
     & 0.9837904585907746D+00,
     & 0.9890905016357307D+00,
     & 0.9927904292352575D+00,
     & 0.9953222650189527D+00 /
      data xvec /
     &  0.0D+00,
     &  0.1D+00,
     &  0.2D+00,
     &  0.3D+00,
     &  0.4D+00,
     &  0.5D+00,
     &  0.6D+00,
     &  0.7D+00,
     &  0.8D+00,
     &  0.9D+00,
     &  1.0D+00,
     &  1.1D+00,
     &  1.2D+00,
     &  1.3D+00,
     &  1.4D+00,
     &  1.5D+00,
     &  1.6D+00,
     &  1.7D+00,
     &  1.8D+00,
     &  1.9D+00,
     &  2.0D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        x = 0.0D+00
        fx = 0.0D+00
      else
        x = xvec(n_data)
        fx = bvec(n_data)
      end if

      return
      end
      subroutine erfc_values ( n_data, x, fx )

c*********************************************************************72
c
cc ERFC_VALUES returns some values of the ERFC or "complementary error" function for testing.
c
c  Discussion:
c
c    The complementary error function is defined by:
c
c      ERFC(X) = 1 - ( 2 / sqrt ( PI ) * integral ( 0 <= T <= X ) exp ( - T^2 ) dT
c
c    In Mathematica, the function can be evaluated by:
c
c      Erfc[x]
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
c  Reference:
c
c    Milton Abramowitz, Irene Stegun,
c    Handbook of Mathematical Functions,
c    National Bureau of Standards, 1964,
c    ISBN: 0-486-61272-4,
c    LC: QA47.A34.
c
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Cambridge University Press, 1999,
c    ISBN: 0-521-64314-7,
c    LC: QA76.95.W65.
c
c  Parameters:
c
c    Input/output, integer N_DATA.
c    On input, if N_DATA is 0, the first test data is returned, and
c    N_DATA is set to the index of the test data.  On each subsequent
c    call, N_DATA is incremented and that test data is returned.  When
c    there is no more test data, N_DATA is set to 0.
c
c    Output, double precision X, the argument of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 21 )

      double precision bvec ( n_max )
      double precision fx
      integer n_data
      double precision x
      double precision xvec ( n_max )

      data bvec /
     &  1.000000000000000D+00,
     &  0.7772974107895215D+00,
     &  0.5716076449533315D+00,
     &  0.3961439091520741D+00,
     &  0.2578990352923395D+00,
     &  0.1572992070502851D+00,
     &  0.08968602177036462D+00,
     &  0.04771488023735119D+00,
     &  0.02365161665535599D+00,
     &  0.01090949836426929D+00,
     &  0.004677734981047266D+00,
     &  0.001862846297981891D+00,
     &  0.0006885138966450786D+00,
     &  0.0002360344165293492D+00,
     &  0.00007501319466545902D+00,
     &  0.00002209049699858544D+00,
     &  6.025761151762095D-06,
     &  1.521993362862285D-06,
     &  3.558629930076853D-07,
     &  7.700392745696413D-08,
     &  1.541725790028002D-08 /
      data xvec /
     &  0.0D+00,
     &  0.2D+00,
     &  0.4D+00,
     &  0.6D+00,
     &  0.8D+00,
     &  1.0D+00,
     &  1.2D+00,
     &  1.4D+00,
     &  1.6D+00,
     &  1.8D+00,
     &  2.0D+00,
     &  2.2D+00,
     &  2.4D+00,
     &  2.6D+00,
     &  2.8D+00,
     &  3.0D+00,
     &  3.2D+00,
     &  3.4D+00,
     &  3.6D+00,
     &  3.8D+00,
     &  4.0D+00  /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        x = 0.0D+00
        fx = 0.0D+00
      else
        x = xvec(n_data)
        fx = bvec(n_data)
      end if

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
