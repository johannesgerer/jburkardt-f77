      program main

c*********************************************************************72
c
cc MAIN is the main program for FN_PRB.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FN_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( * ,'(a)' ) '  Test the FN library.'

      call acos_test ( )
      call acosh_test ( )
      call ai_test ( )
      call aid_test ( )
      call asin_test ( )
      call asinh_test ( )
      call atan_test ( )
      call atan2_test ( )
      call atanh_test ( )
      call besi0_test ( )
      call besi1_test ( )
      call besj0_test ( )
      call besj1_test ( )
      call besk_test ( )
      call besk0_test ( )
      call besk1_test ( )
      call besy0_test ( )
      call besy1_test ( )
      call beta_test ( )
      call betai_test ( )
      call bi_test ( )
      call bid_test ( )
      call binom_test ( )
      call cbrt_test ( )
      call chi_test ( )
      call chu_test ( )
      call ci_test ( )
      call cin_test ( )
      call cinh_test ( )
      call cos_test ( )
      call cos_deg_test ( )
      call cosh_test ( )
      call cot_test ( )
      call dawson_test ( )
      call e1_test ( )
      call ei_test ( )
      call erf_test ( )
      call erfc_test ( )
      call exp_test ( )
      call fac_test ( )
      call gamma_test ( )
      call gamma_inc_test ( )
      call gamma_inc_tricomi_test ( )
      call int_test ( )
      call lbeta_test ( )
      call li_test ( )
      call lngam_test ( )
      call log_test ( )
      call log10_test ( )
      call poch_test ( )
      call psi_test ( )
      call rand_test ( )
      call shi_test ( )
      call si_test ( )
      call sin_test ( )
      call sin_deg_test ( )
      call sinh_test ( )
      call spence_test ( )
      call sqrt_test ( )
      call tan_test ( )
      call tanh_test ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FN_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine acos_test ( )

c*********************************************************************72
c
cc ACOS_TEST tests R4_ACOS and R8_ACOS.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_acos
      double precision r8_acos
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ACOS_TEST:'
      write ( *, '(a)' ) '  Test ARCCOS_VALUES, R4_ACOS, R8_ACOS.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X      ARCCOS(X)'
      write ( *, '(a)' ) '                   R4_ACOS(X)         Diff'
      write ( *, '(a)' ) '                   R8_ACOS(X)         Diff'

      n_data = 0

10    continue

        call arccos_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_acos ( real ( x ) )
        fx3 = r8_acos ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine acosh_test ( )

c*********************************************************************72
c
cc ACOSH_TEST tests R4_ACOSH and R8_ACOSH.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_acosh
      double precision r8_acosh
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ACOSH_TEST:'
      write ( *, '(a)' ) '  Test ARCCOSH_VALUES, R4_ACOSH, R8_ACOSH'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X      ARCCOSH(X)'
      write ( *, '(a)' ) '                   R4_ACOSH(X)        Diff'
      write ( *, '(a)' ) '                   R8_ACOSH(X)        Diff'

      n_data = 0

10    continue

        call arccosh_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_acosh ( real ( x ) )
        fx3 = r8_acosh ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine ai_test ( )

c*********************************************************************72
c
cc AI_TEST tests R4_AI and R8_AI.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_ai
      double precision r8_ai
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'AI_TEST:'
      write ( *, '(a)' ) '  Test AIRY_AI_VALUES, R4_AI, R8_AI'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X      AIRY_AI(X)'
      write ( *, '(a)' ) '                      R4_AI(X)        Diff'
      write ( *, '(a)' ) '                      R8_AI(X)        Diff'

      n_data = 0

10    continue

        call airy_ai_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_ai ( real ( x ) )
        fx3 = r8_ai ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine aid_test ( )

c*********************************************************************72
c
cc AID_TEST tests R4_AID and R8_AID.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_aid
      double precision r8_aid
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'AID_TEST:'
      write ( *, '(a)' ) '  Test AIRY_AI_PRIME_VALUES, R4_AID, R8_AID'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X     AIRY_AID(X)'
      write ( *, '(a)' ) '                     R4_AID(X)        Diff'
      write ( *, '(a)' ) '                     R8_AID(X)        Diff'

      n_data = 0

10    continue

        call airy_ai_prime_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_aid ( real ( x ) )
        fx3 = r8_aid ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine asin_test ( )

c*********************************************************************72
c
cc ASIN_TEST tests R4_ASIN and R8_ASIN.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_asin
      double precision r8_asin
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASIN_TEST:'
      write ( *, '(a)' ) '  Test ARCSIN_VALUES, R4_ASIN, R8_ASIN'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X      ARCSIN(X)'
      write ( *, '(a)' ) '                   R4_ASIN(X)         Diff'
      write ( *, '(a)' ) '                   R8_ASIN(X)         Diff'

      n_data = 0

10    continue

        call arcsin_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_asin ( real ( x ) )
        fx3 = r8_asin ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine asinh_test ( )

c*********************************************************************72
c
cc ASINH_TEST tests R4_ASINH and R8_ASINH.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_asinh
      double precision r8_asinh
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASINH_TEST:'
      write ( *, '(a)' ) '  Test ARCSINH_VALUES, R4_ASINH, R8_ASINH'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X     ARCSINH(X)'
      write ( *, '(a)' ) '                  R4_ASINH(X)         Diff'
      write ( *, '(a)' ) '                  R8_ASINH(X)         Diff'

      n_data = 0

10    continue

        call arcsinh_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_asinh ( real ( x ) )
        fx3 = r8_asinh ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine atan_test ( )

c*********************************************************************72
c
cc ATAN_TEST tests R4_ATAN and R8_ATAN.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_atan
      double precision r8_atan
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ATAN_TEST:'
      write ( *, '(a)' ) '  Test ARCTAN_VALUES, R4_ATAN, R8_ATAN'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X      ARCTAN(X)'
      write ( *, '(a)' ) '                   R4_ATAN(X)         Diff'
      write ( *, '(a)' ) '                   R8_ATAN(X)         Diff'

      n_data = 0

10    continue

        call arctan_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_atan ( real ( x ) )
        fx3 = r8_atan ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine atan2_test ( )

c*********************************************************************72
c
cc ATAN2_TEST tests R4_ATAN2 and R8_ATAN2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_atan2
      double precision r8_atan2
      double precision x
      double precision y

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ATAN2_TEST:'
      write ( *, '(a)' ) '  Test ARCTAN2_VALUES, R4_ATAN2, R8_ATAN2'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X               Y   ARCTAN2(Y,X)'
      write ( *, '(a)' )
     &  '                                R4_ATAN2(Y,X)         Diff'
      write ( *, '(a)' )
     &  '                                R8_ATAN2(Y,X)         Diff'

      n_data = 0

10    continue

        call arctan2_values ( n_data, x, y, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_atan2 ( real ( y ), real ( x ) )
        fx3 = r8_atan2 ( y, x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,f14.4,2x,g14.6)' )     x, y, fx1
        write ( *, '(32x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(32x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine atanh_test ( )

c*********************************************************************72
c
cc ATANH_TEST tests R4_ATANH and R8_ATANH.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_atanh
      double precision r8_atanh
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ATANH_TEST:'
      write ( *, '(a)' ) '  Test ARCTANH_VALUES, R4_ATANH, R8_ATANH'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X     ARCTANH(X)'
      write ( *, '(a)' ) '                  R4_ATANH(X)         Diff'
      write ( *, '(a)' ) '                  R8_ATANH(X)         Diff'

      n_data = 0

10    continue

        call arctanh_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_atanh ( real ( x ) )
        fx3 = r8_atanh ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine besi0_test ( )

c*********************************************************************72
c
cc BESI0_TEST tests R4_BESI0 and R8_BESI0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_besi0
      double precision r8_besi0
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BESI0_TEST:'
      write ( *, '(a)' ) '  Test BESI0_VALUES, R4_BESI0, R8_BESI0'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X       BESI0(X)'
      write ( *, '(a)' ) '                  R4_BESI0(X)         Diff'
      write ( *, '(a)' ) '                  R8_BESI0(X)         Diff'

      n_data = 0

10    continue

        call bessel_i0_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_besi0 ( real ( x ) )
        fx3 = r8_besi0 ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine besi1_test ( )

c*********************************************************************72
c
cc BESI1_TEST tests R4_BESI1 and R8_BESI1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_besi1
      double precision r8_besi1
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BESI1_TEST:'
      write ( *, '(a)' ) '  Test BESI1_VALUES, R4_BESI1, R8_BESI1'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X       BESI1(X)'
      write ( *, '(a)' ) '                  R4_BESI1(X)         Diff'
      write ( *, '(a)' ) '                  R8_BESI1(X)         Diff'

      n_data = 0

10    continue

        call bessel_i1_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_besi1 ( real ( x ) )
        fx3 = r8_besi1 ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine besj0_test ( )

c*********************************************************************72
c
cc BESJ0_TEST tests R4_BESJ0 and R8_BESJ0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_besj0
      double precision r8_besj0
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BESJ0_TEST:'
      write ( *, '(a)' ) '  Test BESJ0_VALUES, R4_BESJ0, R8_BESJ0'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X       BESJ0(X)'
      write ( *, '(a)' ) '                  R4_BESJ0(X)         Diff'
      write ( *, '(a)' ) '                  R8_BESJ0(X)         Diff'

      n_data = 0

10    continue

        call bessel_j0_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_besj0 ( real ( x ) )
        fx3 = r8_besj0 ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine besj1_test ( )

c*********************************************************************72
c
cc BESJ1_TEST tests R4_BESJ1 and R8_BESJ1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_besj1
      double precision r8_besj1
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BESJ1_TEST:'
      write ( *, '(a)' ) '  Test BESJ1_VALUES, R4_BESJ1, R8_BESJ1'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X       BESJ1(X)'
      write ( *, '(a)' ) '                  R4_BESJ1(X)         Diff'
      write ( *, '(a)' ) '                  R8_BESJ1(X)         Diff'

      n_data = 0

10    continue

        call bessel_j1_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_besj1 ( real ( x ) )
        fx3 = r8_besj1 ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine besk_test ( )

c*********************************************************************72
c
cc BESK_TEST tests R4_BESK and R8_BESK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 November 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      double precision nu
      real r4_besk
      double precision r8_besk
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BESK_TEST:'
      write ( *, '(a)' ) '  Test BESK_KX_VALUES, R4_BESK, R8_BESK'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '              NU             X       BESK(X)'
      write ( *, '(a)' ) 
     &  '                                  R4_BESK(X)         Diff'
      write ( *, '(a)' ) 
     &  '                                  R8_BESK(X)         Diff'

      n_data = 0

10    continue

        call bessel_kx_values ( n_data, nu, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_besk ( real ( nu ), real ( x ) )
        fx3 = r8_besk ( nu, x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,f14.4,2x,g14.6)' ) 
     &    nu, x, fx1
        write ( *, '(16x,16x,2x,g14.6,2x,g14.6)' )
     &    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,16x,2x,g14.6,2x,g14.6)' )
     &    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine besk0_test ( )

c*********************************************************************72
c
cc BESK0_TEST tests R4_BESK0 and R8_BESK0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_besk0
      double precision r8_besk0
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BESK0_TEST:'
      write ( *, '(a)' ) '  Test BESK0_VALUES, R4_BESK0, R8_BESK0'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X       BESK0(X)'
      write ( *, '(a)' ) '                  R4_BESK0(X)         Diff'
      write ( *, '(a)' ) '                  R8_BESK0(X)         Diff'

      n_data = 0

10    continue

        call bessel_k0_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_besk0 ( real ( x ) )
        fx3 = r8_besk0 ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine besk1_test ( )

c*********************************************************************72
c
cc BESK1_TEST tests R4_BESK1 and R8_BESK1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_besk1
      double precision r8_besk1
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BESK1_TEST:'
      write ( *, '(a)' ) '  Test BESK1_VALUES, R4_BESK1, R8_BESK1'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X       BESK1(X)'
      write ( *, '(a)' ) '                  R4_BESK1(X)         Diff'
      write ( *, '(a)' ) '                  R8_BESK1(X)         Diff'

      n_data = 0

10    continue

        call bessel_k1_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_besk1 ( real ( x ) )
        fx3 = r8_besk1 ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine besy0_test ( )

c*********************************************************************72
c
cc BESY0_TEST tests R4_BESY0 and R8_BESY0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_besy0
      double precision r8_besy0
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BESY0_TEST:'
      write ( *, '(a)' ) '  Test BESY0_VALUES, R4_BESY0, R8_BESY0'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X       BESY0(X)'
      write ( *, '(a)' ) '                  R4_BESY0(X)         Diff'
      write ( *, '(a)' ) '                  R8_BESY0(X)         Diff'

      n_data = 0

10    continue

        call bessel_y0_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_besy0 ( real ( x ) )
        fx3 = r8_besy0 ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine besy1_test ( )

c*********************************************************************72
c
cc BESY1_TEST tests R4_BESY1 and R8_BESY1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_besy1
      double precision r8_besy1
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BESY1_TEST:'
      write ( *, '(a)' ) '  Test BESY1_VALUES, R4_BESY1, R8_BESY1'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X       BESY1(X)'
      write ( *, '(a)' ) '                  R4_BESY1(X)         Diff'
      write ( *, '(a)' ) '                  R8_BESY1(X)         Diff'

      n_data = 0

10    continue

        call bessel_y1_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_besy1 ( real ( x ) )
        fx3 = r8_besy1 ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine beta_test ( )

c*********************************************************************72
c
cc BETA_TEST tests R4_BETA and R8_BETA.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision b
      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_beta
      double precision r8_beta

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BETA_TEST:'
      write ( *, '(a)' ) '  Test BETA_VALUES, R4_BETA, R8_BETA.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             A               B      BETA(A,B)'
      write ( *, '(a)' )
     &  '                                 R4_BETA(A,B)         Diff'
      write ( *, '(a)' )
     &  '                                 R8_BETA(A,B)         Diff'

      n_data = 0

10    continue

        call beta_values ( n_data, a, b, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_beta ( real ( a ), real ( b ) )
        fx3 = r8_beta ( a, b )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.4,2x,g14.6)' )  a, b, fx1
        write ( *, '(2x,14x,2x,14x,2x,g14.6,2x,g14.6)' )
     &    fx2, abs ( fx1 - fx2 )
        write ( *, '(2x,14x,2x,14x,2x,g14.6,2x,g14.6)' )
     &    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine betai_test ( )

c*********************************************************************72
c
cc BETAI_TEST tests R4_BETAI and R8_BETAI.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision b
      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_betai
      double precision r8_betai
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BETAI_TEST:'
      write ( *, '(a)' ) '  Test BETA_INC_VALUES, R4_BETAI, R8_BETAI.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '             X               A               B     BETAI(A,B)'
      write ( *, '(a)' )
     &  '                                                R4_BETAI(A,B)'
     &  // '         Diff'
      write ( *, '(a)' )
     &  '                                                R8_BETAI(A,B)'
     &  // '         Diff'

      n_data = 0

10    continue

        call beta_inc_values ( n_data, a, b, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_betai ( real ( x ), real ( a ), real ( b ) )
        fx3 = r8_betai ( x, a, b )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,f14.4,2x,g14.4,2x,g14.6)' )
     &    x, a, b, fx1
        write ( *, '(2x,14x,2x,14x,2x,14x,2x,g14.6,2x,g14.6)' )
     &    fx2, abs ( fx1 - fx2 )
        write ( *, '(2x,14x,2x,14x,2x,14x,2x,g14.6,2x,g14.6)' )
     &    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine bi_test ( )

c*********************************************************************72
c
cc BI_TEST tests R4_BI and R8_BI.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_bi
      double precision r8_bi
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BI_TEST:'
      write ( *, '(a)' ) '  Test AIRY_BI_VALUES, R4_BI, R8_BI'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X      AIRY_BI(X)'
      write ( *, '(a)' ) '                      R4_BI(X)        Diff'
      write ( *, '(a)' ) '                      R8_BI(X)        Diff'

      n_data = 0

10    continue

        call airy_bi_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_bi ( real ( x ) )
        fx3 = r8_bi ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine bid_test ( )

c*********************************************************************72
c
cc BID_TEST tests R4_BID and R8_BID.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_bid
      double precision r8_bid
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BID_TEST:'
      write ( *, '(a)' ) '  Test AIRY_BI_PRIME_VALUES, R4_BID, R8_BID'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X     AIRY_BID(X)'
      write ( *, '(a)' ) '                     R4_BID(X)        Diff'
      write ( *, '(a)' ) '                     R8_BID(X)        Diff'

      n_data = 0

10    continue

        call airy_bi_prime_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_bid ( real ( x ) )
        fx3 = r8_bid ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine binom_test ( )

c*********************************************************************72
c
cc BINOM_TEST tests R4_BINOM and R8_BINOM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer a
      integer b
      double precision diff
      integer fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_binom
      double precision r8_binom

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BINOM_TEST:'
      write ( *, '(a)' ) '  Test BINOM_VALUES, R4_BINOM, R8_BINOM.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             A               B     BINOM(A,B)'
      write ( *, '(a)' )
     &  '                                R4_BINOM(A,B)         Diff'
      write ( *, '(a)' )
     &  '                                R8_BINOM(A,B)         Diff'

      n_data = 0

10    continue

        call binomial_values ( n_data, a, b, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_binom ( a, b )
        fx3 = r8_binom ( a, b )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,i14,2x,i14,2x,i14)' )  a, b, fx1
        write ( *, '(2x,14x,2x,14x,2x,g14.6,2x,g14.6)' )
     &    fx2, abs ( real ( fx1 ) - fx2 )
        write ( *, '(2x,14x,2x,14x,2x,g14.6,2x,g14.6)' )
     &    fx3, abs ( dble ( fx1 ) - fx3 )

      go to 10

20    continue

      return
      end
      subroutine cbrt_test ( )

c*********************************************************************72
c
cc CBRT_TEST tests R4_CBRT and R8_CBRT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_cbrt
      double precision r8_cbrt
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CBRT_TEST:'
      write ( *, '(a)' ) '  Test CBRT_VALUES, R4_CBRT, R8_CBRT.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X        CBRT(X)'
      write ( *, '(a)' ) '                   R4_CBRT(X)         Diff'
      write ( *, '(a)' ) '                   R8_CBRT(X)         Diff'

      n_data = 0

10    continue

        call cbrt_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_cbrt ( real ( x ) )
        fx3 = r8_cbrt ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,g14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine chu_test ( )

c*********************************************************************72
c
cc CHU_TEST tests R4_CHU and R8_CHU.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 October 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision b
      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_chu
      double precision r8_chu
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CHU_TEST:'
      write ( *, '(a)' ) 
     &  '  Test HYPERGEOMETRIC_U_VALUES, R4_CHU, R8_CHU.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '             A               B               X     CHU(A,B,X)'
      write ( *, '(a)' )
     &  '                                                R4_CHU(A,B,X)'
     &  // '         Diff'
      write ( *, '(a)' )
     &  '                                                R8_CHU(A,B,X)'
     &  // '         Diff'

      n_data = 0

10    continue

        call hypergeometric_u_values ( n_data, a, b, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_chu ( real ( a ), real ( b ), real ( x ) )
        fx3 = r8_chu ( a, b, x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,f14.4,2x,g14.4,2x,g14.6)' )
     &    a, b, x, fx1
        write ( *, '(2x,14x,2x,14x,2x,14x,2x,g14.6,2x,g14.6)' )
     &    fx2, abs ( fx1 - fx2 )
        write ( *, '(2x,14x,2x,14x,2x,14x,2x,g14.6,2x,g14.6)' )
     &    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine chi_test ( )

c*********************************************************************72
c
cc CHI_TEST tests R4_CHI and R8_CHI.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_chi
      double precision r8_chi
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CHI_TEST:'
      write ( *, '(a)' )
     &  '  Test CHI_VALUES, R4_CHI, R8_CHI.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X          CHI(X)'
      write ( *, '(a)' ) '                     R4_CHI(X)        Diff'
      write ( *, '(a)' ) '                     R8_CHI(X)        Diff'

      n_data = 0

10    continue

        call chi_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_chi ( real ( x ) )
        fx3 = r8_chi ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine ci_test ( )

c*********************************************************************72
c
cc CI_TEST tests R4_CI and R8_CI.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_ci
      double precision r8_ci
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CI_TEST:'
      write ( *, '(a)' )
     &  '  Test CI_VALUES, R4_CI, R8_CI.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X           CI(X)'
      write ( *, '(a)' ) '                      R4_CI(X)        Diff'
      write ( *, '(a)' ) '                      R8_CI(X)        Diff'

      n_data = 0

10    continue

        call ci_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_ci ( real ( x ) )
        fx3 = r8_ci ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine cin_test ( )

c*********************************************************************72
c
cc CIN_TEST tests R4_CIN and R8_CIN.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_cin
      double precision r8_cin
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CIN_TEST:'
      write ( *, '(a)' )
     &  '  Test CIN_VALUES, R4_CIN, R8_CIN.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X          CIN(X)'
      write ( *, '(a)' ) '                     R4_CIN(X)        Diff'
      write ( *, '(a)' ) '                     R8_CIN(X)        Diff'

      n_data = 0

10    continue

        call cin_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_cin ( real ( x ) )
        fx3 = r8_cin ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine cinh_test ( )

c*********************************************************************72
c
cc CINH_TEST tests R4_CINH and R8_CINH.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_cinh
      double precision r8_cinh
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CINH_TEST:'
      write ( *, '(a)' )
     &  '  Test CINH_VALUES, R4_CINH, R8_CINH.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X         CINH(X)'
      write ( *, '(a)' ) '                    R4_CINH(X)        Diff'
      write ( *, '(a)' ) '                    R8_CINH(X)        Diff'

      n_data = 0

10    continue

        call cinh_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_cinh ( real ( x ) )
        fx3 = r8_cinh ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine cos_test ( )

c*********************************************************************72
c
cc COS_TEST tests R4_COS and R8_COS.
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

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_cos
      double precision r8_cos
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'COS_TEST:'
      write ( *, '(a)' ) '  Test COS_VALUES, R4_COS, R8_COS.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X         COS(X)'
      write ( *, '(a)' ) '                    R4_COS(X)         Diff'
      write ( *, '(a)' ) '                    R8_COS(X)         Diff'

      n_data = 0

10    continue

        call cos_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_cos ( real ( x ) )
        fx3 = r8_cos ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine cos_deg_test ( )

c*********************************************************************72
c
cc COS_DEG_TEST tests R4_COS_DEG and R8_COS_DEG.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_cos_deg
      double precision r8_cos_deg
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'COS_DEG_TEST:'
      write ( *, '(a)' )
     &  '  Test COS_DEGREE_VALUES, R4_COS_DEG, R8_COS_DEG.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X     COS_DEG(X)'
      write ( *, '(a)' ) '                R4_COS_DEG(X)         Diff'
      write ( *, '(a)' ) '                R8_COS_DEG(X)         Diff'

      n_data = 0

10    continue

        call cos_degree_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_cos_deg ( real ( x ) )
        fx3 = r8_cos_deg ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine cosh_test ( )

c*********************************************************************72
c
cc COSH_TEST tests R4_COSH and R8_COSH.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_cosh
      double precision r8_cosh
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'COSH_TEST:'
      write ( *, '(a)' ) '  Test COSH_VALUES, R4_COSH, R8_COSH.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X         COSH(X)'
      write ( *, '(a)' ) '                    R4_COSH(X)        Diff'
      write ( *, '(a)' ) '                    R8_COSH(X)        Diff'

      n_data = 0

10    continue

        call cosh_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_cosh ( real ( x ) )
        fx3 = r8_cosh ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine cot_test ( )

c*********************************************************************72
c
cc COT_TEST tests R4_COT and R8_COT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_cot
      double precision r8_cot
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'COT_TEST:'
      write ( *, '(a)' ) '  Test COT_VALUES, R4_COT, R8_COT.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X         COT(X)'
      write ( *, '(a)' ) '                    R4_COT(X)         Diff'
      write ( *, '(a)' ) '                    R8_COT(X)         Diff'

      n_data = 0

10    continue

        call cot_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_cot ( real ( x ) )
        fx3 = r8_cot ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine dawson_test ( )

c*********************************************************************72
c
cc DAWSON_TEST tests R4_DAWSON and R8_DAWSON.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_dawson
      double precision r8_dawson
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DAWSON_TEST:'
      write ( *, '(a)' ) '  Test DAWSON_VALUES, R4_DAWSON, R8_DAWSON.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X      DAWSON(X)'
      write ( *, '(a)' ) '                 R4_DAWSON(X)         Diff'
      write ( *, '(a)' ) '                 R8_DAWSON(X)         Diff'

      n_data = 0

10    continue

        call dawson_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_dawson ( real ( x ) )
        fx3 = r8_dawson ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine e1_test ( )

c*********************************************************************72
c
cc E1_TEST tests R4_E1 and R8_E1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_e1
      double precision r8_e1
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'E1_TEST:'
      write ( *, '(a)' )
     &  '  Test E1_VALUES, R4_E1, R8_E1.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X           E1(X)'
      write ( *, '(a)' ) '                      R4_E1(X)        Diff'
      write ( *, '(a)' ) '                      R8_E1(X)        Diff'

      n_data = 0

10    continue

        call e1_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_e1 ( real ( x ) )
        fx3 = r8_e1 ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine ei_test ( )

c*********************************************************************72
c
cc EI_TEST tests R4_EI and R8_EI.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_ei
      double precision r8_ei
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EI_TEST:'
      write ( *, '(a)' )
     &  '  Test EI_VALUES, R4_EI, R8_EI.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X           EI(X)'
      write ( *, '(a)' ) '                      R4_EI(X)        Diff'
      write ( *, '(a)' ) '                      R8_EI(X)        Diff'

      n_data = 0

10    continue

        call ei_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_ei ( real ( x ) )
        fx3 = r8_ei ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine erf_test ( )

c*********************************************************************72
c
cc ERF_TEST tests R4_ERF and R8_ERF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_erf
      double precision r8_erf
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ERF_TEST:'
      write ( *, '(a)' ) '  Test ERF_VALUES, R4_ERF R8_ERF'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X          ERF(X)'
      write ( *, '(a)' ) '                     R4_ERF(X)        Diff'
      write ( *, '(a)' ) '                     R8_ERF(X)        Diff'

      n_data = 0

10    continue

        call erf_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_erf ( real ( x ) )
        fx3 = r8_erf ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine erfc_test ( )

c*********************************************************************72
c
cc ERFC_TEST tests R4_ERFC and R8_ERFC.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_erfc
      double precision r8_erfc
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ERFC_TEST:'
      write ( *, '(a)' ) '  Test ERFD_VALUES, R4_ERFC R8_ERFC'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X         ERFC(X)'
      write ( *, '(a)' ) '                    R4_ERFC(X)        Diff'
      write ( *, '(a)' ) '                    R8_ERFC(X)        Diff'

      n_data = 0

10    continue

        call erfc_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_erfc ( real ( x ) )
        fx3 = r8_erfc ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine exp_test ( )

c*********************************************************************72
c
cc EXP_TEST tests R4_EXP and R8_EXP.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_exp
      double precision r8_exp
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EXP_TEST:'
      write ( *, '(a)' ) '  Test EXP_VALUES, R4_EXP, R8_EXP'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X          EXP(X)'
      write ( *, '(a)' ) '                     R4_EXP(X)        Diff'
      write ( *, '(a)' ) '                     R8_EXP(X)        Diff'

      n_data = 0

10    continue

        call exp_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_exp ( real ( x ) )
        fx3 = r8_exp ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine fac_test ( )

c*********************************************************************72
c
cc FAC_TEST tests R4_FAC and R8_FAC.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      integer fx1
      real fx2
      double precision fx3
      integer n
      integer n_data
      real r4_fac
      double precision r8_fac

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FAC_TEST:'
      write ( *, '(a)' ) '  Test FACTORIAL_VALUES, R4_FAC, R8_FAC.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             N          FAC(N)'
      write ( *, '(a)' ) '                     R4_FAC(N)        Diff'
      write ( *, '(a)' ) '                     R8_FAC(N)        Diff'

      n_data = 0

10    continue

        call factorial_values ( n_data, n, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_fac ( n )
        fx3 = r8_fac ( n )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,i14,2x,i14)' )       n, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine gamma_test ( )

c*********************************************************************72
c
cc GAMMA_TEST tests R4_GAMMA and R8_GAMMA.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_gamma
      double precision r8_gamma
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GAMMA_TEST:'
      write ( *, '(a)' ) '  Test GAMMA_VALUES, R4_GAMMA, R8_GAMMA.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X        GAMMA(X)'
      write ( *, '(a)' ) '                   R4_GAMMA(X)        Diff'
      write ( *, '(a)' ) '                   R8_GAMMA(X)        Diff'

      n_data = 0

10    continue

        call gamma_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_gamma ( real ( x ) )
        fx3 = r8_gamma ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine gamma_inc_test ( )

c*********************************************************************72
c
cc GAMMA_INC_TEST tests R4_GAMIC and R8_GAMIC.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_gamic
      double precision r8_gamic
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GAMMA_INC_TEST:'
      write ( *, '(a)' ) '  Test GAMMA_INC_VALUES, R4_GAMIC, R8_GAMIC.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             A               X     GAMIC(A,X)'
      write ( *, '(a)' )
     &  '                                R4_GAMIC(A,X)         Diff'
      write ( *, '(a)' )
     &  '                                R8_GAMIC(A,X)         Diff'

      n_data = 0

10    continue

        call gamma_inc_values ( n_data, a, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_gamic ( real ( a ), real ( x ) )
        fx3 = r8_gamic ( a, x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.4,2x,g14.6)' )  a, x, fx1
        write ( *, '(2x,14x,2x,14x,2x,g14.6,2x,g14.6)' )
     &    fx2, abs ( fx1 - fx2 )
        write ( *, '(2x,14x,2x,14x,2x,g14.6,2x,g14.6)' )
     &    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine gamma_inc_tricomi_test ( )

c*********************************************************************72
c
cc GAMMA_INC_TRICOMI_TEST tests R4_GAMIT and R8_GAMIT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_gamit
      double precision r8_gamit
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GAMMA_INC_TRICOMI_TEST:'
      write ( *, '(a)' )
     &  '  Test GAMMA_INC_TRICOMI_VALUES, R4_GAMIT, R8_GAMIT.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             A               X     GAMIT(A,X)'
      write ( *, '(a)' )
     &  '                                R4_GAMIT(A,X)         Diff'
      write ( *, '(a)' )
     &  '                                R8_GAMIT(A,X)         Diff'

      n_data = 0

10    continue

        call gamma_inc_tricomi_values ( n_data, a, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_gamit ( real ( a ), real ( x ) )
        fx3 = r8_gamit ( a, x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.4,2x,g14.6)' )  a, x, fx1
        write ( *, '(2x,14x,2x,14x,2x,g14.6,2x,g14.6)' )
     &    fx2, abs ( fx1 - fx2 )
        write ( *, '(2x,14x,2x,14x,2x,g14.6,2x,g14.6)' )
     &    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine int_test ( )

c*********************************************************************72
c
cc INT_TEST tests R4_INT and R8_INT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 October 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_int
      double precision r8_int
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'INT_TEST:'
      write ( *, '(a)' ) '  Test INT_VALUES, R4_INT, R8_INT.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X         INT(X)'
      write ( *, '(a)' ) '                    R4_INT(X)         Diff'
      write ( *, '(a)' ) '                    R8_INT(X)         Diff'

      n_data = 0

10    continue

        call int_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_int ( real ( x ) )
        fx3 = r8_int ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine lbeta_test ( )

c*********************************************************************72
c
cc LBETA_TEST tests R4_LBETA and R8_LBETA.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision b
      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_lbeta
      double precision r8_lbeta

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LBETA_TEST:'
      write ( *, '(a)' ) '  Test BETA_LOG_VALUES, R4_LBETA, R8_LBETA.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             A               B     LBETA(A,B)'
      write ( *, '(a)' )
     &  '                                R4_LBETA(A,B)         Diff'
      write ( *, '(a)' )
     &  '                                R8_LBETA(A,B)         Diff'

      n_data = 0

10    continue

        call beta_log_values ( n_data, a, b, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_lbeta ( real ( a ), real ( b ) )
        fx3 = r8_lbeta ( a, b )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.4,2x,g14.6)' )  a, b, fx1
        write ( *, '(2x,14x,2x,14x,2x,g14.6,2x,g14.6)' )
     &    fx2, abs ( fx1 - fx2 )
        write ( *, '(2x,14x,2x,14x,2x,g14.6,2x,g14.6)' )
     &    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine li_test ( )

c*********************************************************************72
c
cc LI_TEST tests R4_LI and R8_LI.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_li
      double precision r8_li
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LI_TEST:'
      write ( *, '(a)' )
     &  '  Test LOGARITHMIC_INTEGRAL_VALUES, R4_LI, R8_LI.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X           LI(X)'
      write ( *, '(a)' ) '                      R4_LI(X)        Diff'
      write ( *, '(a)' ) '                      R8_LI(X)        Diff'

      n_data = 0

10    continue

        call logarithmic_integral_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_li ( real ( x ) )
        fx3 = r8_li ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine lngam_test ( )

c*********************************************************************72
c
cc LNGAM_TEST tests R4_LNGAM and R8_LNGAM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_lngam
      double precision r8_lngam
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LNGAM_TEST:'
      write ( *, '(a)' ) '  Test GAMMA_LOG_VALUES, R4_LNGAM, R8_LNGAM.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X        LNGAM(X)'
      write ( *, '(a)' ) '                   R4_LNGAM(X)        Diff'
      write ( *, '(a)' ) '                   R8_LNGAM(X)        Diff'

      n_data = 0

10    continue

        call gamma_log_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_lngam ( real ( x ) )
        fx3 = r8_lngam ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine log_test ( )

c*********************************************************************72
c
cc LOG_TEST tests R4_LOG and R8_LOG.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_log
      double precision r8_log
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LOG_TEST:'
      write ( *, '(a)' ) '  Test LOG_VALUES, R4_LOG, R8_LOG.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X          LOG(X)'
      write ( *, '(a)' ) '                     R4_LOG(X)        Diff'
      write ( *, '(a)' ) '                     R8_LOG(X)        Diff'

      n_data = 0

10    continue

        call log_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_log ( real ( x ) )
        fx3 = r8_log ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine log10_test ( )

c*********************************************************************72
c
cc LOG10_TEST tests R4_LOG10 and R8_LOG10.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_log10
      double precision r8_log10
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LOG10_TEST:'
      write ( *, '(a)' ) '  Test LOG10_VALUES, R4_LOG10, R8_LOG10.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X        LOG10(X)'
      write ( *, '(a)' ) '                   R4_LOG10(X)        Diff'
      write ( *, '(a)' ) '                   R8_LOG10(X)        Diff'

      n_data = 0

10    continue

        call log10_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_log10 ( real ( x ) )
        fx3 = r8_log10 ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine poch_test ( )

c*********************************************************************72
c
cc POCH_TEST tests R4_POCH and R8_POCH.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_poch
      double precision r8_poch
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POCH_TEST:'
      write ( *, '(a)' ) '  Test POCHHAMMER_VALUES, R4_POCH, R8_POCH.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             A               X      POCH(A,X)'
      write ( *, '(a)' )
     &  '                                 R4_POCH(A,X)         Diff'
      write ( *, '(a)' )
     &  '                                 R8_POCH(A,X)         Diff'

      n_data = 0

10    continue

        call pochhammer_values ( n_data, a, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_poch ( real ( a ), real ( x ) )
        fx3 = r8_poch ( a, x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.4,2x,g14.6)' )  a, x, fx1
        write ( *, '(2x,14x,2x,14x,2x,g14.6,2x,g14.6)' )
     &    fx2, abs ( fx1 - fx2 )
        write ( *, '(2x,14x,2x,14x,2x,g14.6,2x,g14.6)' )
     &    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine psi_test ( )

c*********************************************************************72
c
cc PSI_TEST tests R4_PSI and R8_PSI.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_psi
      double precision r8_psi
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PSI_TEST:'
      write ( *, '(a)' ) '  Test PSI_VALUES, R4_PSI, R8_PSI.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X         PSI(X)'
      write ( *, '(a)' ) '                    R4_PSI(X)         Diff'
      write ( *, '(a)' ) '                    R8_PSI(X)         Diff'

      n_data = 0

10    continue

        call psi_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        fx2 = r4_psi ( real ( x ) )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        fx3 = r8_psi ( x )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine rand_test ( )

c*********************************************************************72
c
cc RAND_TEST tests R4_RAND.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 October 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      real average
      integer i
      integer i_value(7)
      integer k
      real r
      real r_value(7)
      real r4_rand
      real variance

      save i_value
      save r_value

      data i_value / 1, 2, 3, 4, 10, 100, 1000 /
      data r_value / 
     &  0.0004127026E+00, 
     &  0.6750836372E+00, 
     &  0.1614754200E+00, 
     &  0.9086198807E+00,
     &  0.5527787209E+00,
     &  0.3600893021E+00,  
     &  0.2176990509E+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RAND_TEST:'
      write ( *, '(a)' ) '  Test R4_RAND.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '             I         R4_RAND        Expected'
      write ( *, '(a)' ) ' '

      k = 1

      do i = 1, 1000

        r = r4_rand ( 0.0E+00 )

        if ( i == i_value(k) ) then
          write ( *, '(2x,i14,2x,g14.6,2x,g14.6)' ) i, r, r_value(k)
          k = k + 1
        end if

      end do

      average = 0.0E+00
      do i = 1, 1000000
        r = r4_rand ( 0.0E+00 )
        average = average + r
      end do
      average = average / 1000000.0E+00
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,2x,g14.6)' ) 
     &  '       Average =  ', average, 0.5

      variance = 0.0E+00
      do i = 1, 1000000
        r = r4_rand ( 0.0 )
        variance = variance + ( r - average )**2
      end do
      variance = variance / 1000000.0E+00
      write ( *, '(a,g14.6,2x,g14.6)' ) 
     &  '       Variance = ', variance, 1.0E+00 / 12.0E+00

      return
      end
      subroutine shi_test ( )

c*********************************************************************72
c
cc SHI_TEST tests R4_SHI and R8_SHI.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_shi
      double precision r8_shi
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SHI_TEST:'
      write ( *, '(a)' )
     &  '  Test SHI_VALUES, R4_SHI, R8_SHI.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X          SHI(X)'
      write ( *, '(a)' ) '                     R4_SHI(X)        Diff'
      write ( *, '(a)' ) '                     R8_SHI(X)        Diff'

      n_data = 0

10    continue

        call shi_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_shi ( real ( x ) )
        fx3 = r8_shi ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine si_test ( )

c*********************************************************************72
c
cc SI_TEST tests R4_SI and R8_SI.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_si
      double precision r8_si
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SI_TEST:'
      write ( *, '(a)' )
     &  '  Test SI_VALUES, R4_SI, R8_SI.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X           SI(X)'
      write ( *, '(a)' ) '                      R4_SI(X)        Diff'
      write ( *, '(a)' ) '                      R8_SI(X)        Diff'

      n_data = 0

10    continue

        call si_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_si ( real ( x ) )
        fx3 = r8_si ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine sin_test ( )

c*********************************************************************72
c
cc SIN_TEST tests R4_SIN and R8_SIN.
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

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_sin
      double precision r8_sin
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SIN_TEST:'
      write ( *, '(a)' ) '  Test SIN_VALUES, R4_SIN, R8_SIN.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X         SIN(X)'
      write ( *, '(a)' ) '                    R4_SIN(X)         Diff'
      write ( *, '(a)' ) '                    R8_SIN(X)         Diff'

      n_data = 0

10    continue

        call sin_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_sin ( real ( x ) )
        fx3 = r8_sin ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine sin_deg_test ( )

c*********************************************************************72
c
cc SIN_DEG_TEST tests R4_SIN_DEG and R8_SIN_DEG.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_sin_deg
      double precision r8_sin_deg
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SIN_DEG_TEST:'
      write ( *, '(a)' )
     &  '  Test SIN_DEGREE_VALUES, R4_SIN_DEG, R8_SIN_DEG.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X     SIN_DEG(X)'
      write ( *, '(a)' ) '                R4_SIN_DEG(X)         Diff'
      write ( *, '(a)' ) '                R8_SIN_DEG(X)         Diff'

      n_data = 0

10    continue

        call sin_degree_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_sin_deg ( real ( x ) )
        fx3 = r8_sin_deg ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine sinh_test ( )

c*********************************************************************72
c
cc SINH_TEST tests R4_SINH and R8_SINH.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_sinh
      double precision r8_sinh
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SINH_TEST:'
      write ( *, '(a)' ) '  Test SINH_VALUES, R4_SINH, R8_SINH.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X         SINH(X)'
      write ( *, '(a)' ) '                    R4_SINH(X)        Diff'
      write ( *, '(a)' ) '                    R8_SINH(X)        Diff'

      n_data = 0

10    continue

        call sinh_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_sinh ( real ( x ) )
        fx3 = r8_sinh ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine spence_test ( )

c*********************************************************************72
c
cc SPENCE_TEST tests R4_SPENCE and R8_SPENCE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_spence
      double precision r8_spence
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPENCE_TEST:'
      write ( *, '(a)' )
     &  '  Test DILOGARITHM_VALUES, R4_SPENCE, R8_SPENCE.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X       SPENCE(X)'
      write ( *, '(a)' ) '                  R4_SPENCE(X)        Diff'
      write ( *, '(a)' ) '                  R8_SPENCE(X)        Diff'

      n_data = 0

10    continue

        call dilogarithm_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_spence ( real ( x ) )
        fx3 = r8_spence ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine sqrt_test ( )

c*********************************************************************72
c
cc SQRT_TEST tests R4_SQRT and R8_SQRT.
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

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_sqrt
      double precision r8_sqrt
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SQRT_TEST:'
      write ( *, '(a)' ) '  Test SQRT_VALUES, R4_SQRT, R8_SQRT.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X         SQRT(X)'
      write ( *, '(a)' ) '                    R4_SQRT(X)        Diff'
      write ( *, '(a)' ) '                    R8_SQRT(X)        Diff'

      n_data = 0

10    continue

        call sqrt_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_sqrt ( real ( x ) )
        fx3 = r8_sqrt ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine tan_test ( )

c*********************************************************************72
c
cc TAN_TEST tests R4_TAN and R8_TAN.
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

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_tan
      double precision r8_tan
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TAN_TEST:'
      write ( *, '(a)' ) '  Test TAN_VALUES, R4_TAN, R8_TAN.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X         TAN(X)'
      write ( *, '(a)' ) '                    R4_TAN(X)         Diff'
      write ( *, '(a)' ) '                    R8_TAN(X)         Diff'

      n_data = 0

10    continue

        call tan_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_tan ( real ( x ) )
        fx3 = r8_tan ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
      subroutine tanh_test ( )

c*********************************************************************72
c
cc TANH_TEST tests R4_TANH and R8_TANH.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision diff
      double precision fx1
      real fx2
      double precision fx3
      integer n_data
      real r4_tanh
      double precision r8_tanh
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TANH_TEST:'
      write ( *, '(a)' ) '  Test TANH_VALUES, R4_TANH, R8_TANH.'
      write ( *, '(a)' ) '  TANH_VALUES returns some exact values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             X         TANH(X)'
      write ( *, '(a)' ) '                    R4_TANH(X)         Diff'
      write ( *, '(a)' ) '                    R8_TANH(X)         Diff'

      n_data = 0

10    continue

        call tanh_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r4_tanh ( real ( x ) )
        fx3 = r8_tanh ( x )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
        write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

      go to 10

20    continue

      return
      end
