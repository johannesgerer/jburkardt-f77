      function d1mach ( i )

c*********************************************************************72
c
cc D1MACH returns double precision machine-dependent constants.
c
c  Discussion:
c
c    D1MACH can be used to obtain machine-dependent parameters
c    for the local machine environment.  It is a function
c    with one input argument, and can be called as follows:
c
c      D = D1MACH ( I )
c
c    where I=1,...,5.  The output value of D above is
c    determined by the input value of I:.
c
c    D1MACH ( 1) = B**(EMIN-1), the smallest positive magnitude.
c    D1MACH ( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
c    D1MACH ( 3) = B**(-T), the smallest relative spacing.
c    D1MACH ( 4) = B**(1-T), the largest relative spacing.
c    D1MACH ( 5) = LOG10(B)
c
c  Modified:
c
c    06 December 2006
c
c  Author:
c
c    Phyllis Fox, Andrew Hall, Norman Schryer
c
c  Reference:
c
c    Phyllis Fox, Andrew Hall, Norman Schryer,
c    Algorithm 528:
c    Framework for a Portable Library,
c    ACM Transactions on Mathematical Software,
c    Volume 4, Number 2, June 1978, page 176-188.
c
c  Parameters:
c
c    Input, integer I, the index of the desired constant.
c
c    Output, double precision D1MACH, the value of the constant.
c
      implicit none

      double precision d1mach
      integer diver(4)
      double precision dmach(5)
      integer i
      integer large(4)
      integer log10(4)
      integer right(4)
      integer small(4)

      equivalence ( dmach(1), small(1) )
      equivalence ( dmach(2), large(1) )
      equivalence ( dmach(3), right(1) )
      equivalence ( dmach(4), diver(1) )
      equivalence ( dmach(5), log10(1) )
c
c     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
c     3B SERIES AND MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
c     PC 7300), IN WHICH THE MOST SIGNIFICANT BYTE IS STORED FIRST.
c
c === MACHINE = IEEE.MOST-SIG-BYTE-FIRST
c === MACHINE = SUN
c === MACHINE = 68000
c === MACHINE = ATT.3B
c === MACHINE = ATT.7300
c     DATA SMALL(1),SMALL(2) /    1048576,          0 /
c     DATA LARGE(1),LARGE(2) / 2146435071,         -1 /
c     DATA RIGHT(1),RIGHT(2) / 1017118720,          0 /
c     DATA DIVER(1),DIVER(2) / 1018167296,          0 /
c     DATA LOG10(1),LOG10(2) / 1070810131, 1352628735 /
c
c     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES AND 8087-BASED
c     MICROS, SUCH AS THE IBM PC AND AT&T 6300, IN WHICH THE LEAST
c     SIGNIFICANT BYTE IS STORED FIRST.
c
c === MACHINE = IEEE.LEAST-SIG-BYTE-FIRST
c === MACHINE = 8087
c === MACHINE = IBM.PC
c === MACHINE = ATT.6300
c
       data small(1),small(2) /          0,    1048576 /
       data large(1),large(2) /         -1, 2146435071 /
       data right(1),right(2) /          0, 1017118720 /
       data diver(1),diver(2) /          0, 1018167296 /
       data log10(1),log10(2) / 1352628735, 1070810131 /
c
c     MACHINE CONSTANTS FOR AMDAHL MACHINES.
c
c === MACHINE = AMDAHL
c      DATA SMALL(1),SMALL(2) /    1048576,          0 /
c      DATA LARGE(1),LARGE(2) / 2147483647,         -1 /
c      DATA RIGHT(1),RIGHT(2) /  856686592,          0 /
c      DATA DIVER(1),DIVER(2) /  873463808,          0 /
c      DATA LOG10(1),LOG10(2) / 1091781651, 1352628735 /
c
c     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
c
c === MACHINE = BURROUGHS.1700
c      DATA SMALL(1) / ZC00800000 /
c      DATA SMALL(2) / Z000000000 /
c      DATA LARGE(1) / ZDFFFFFFFF /
c      DATA LARGE(2) / ZFFFFFFFFF /
c      DATA RIGHT(1) / ZCC5800000 /
c      DATA RIGHT(2) / Z000000000 /
c      DATA DIVER(1) / ZCC6800000 /
c      DATA DIVER(2) / Z000000000 /
c      DATA LOG10(1) / ZD00E730E7 /
c      DATA LOG10(2) / ZC77800DC0 /
c
c     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.
c
c === MACHINE = BURROUGHS.5700
c      DATA SMALL(1) / O1771000000000000 /
c      DATA SMALL(2) / O0000000000000000 /
c      DATA LARGE(1) / O0777777777777777 /
c      DATA LARGE(2) / O0007777777777777 /
c      DATA RIGHT(1) / O1461000000000000 /
c      DATA RIGHT(2) / O0000000000000000 /
c      DATA DIVER(1) / O1451000000000000 /
c      DATA DIVER(2) / O0000000000000000 /
c      DATA LOG10(1) / O1157163034761674 /
c      DATA LOG10(2) / O0006677466732724 /
c
c     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.
c
c === MACHINE = BURROUGHS.6700
c === MACHINE = BURROUGHS.7700
c      DATA SMALL(1) / O1771000000000000 /
c      DATA SMALL(2) / O7770000000000000 /
c      DATA LARGE(1) / O0777777777777777 /
c      DATA LARGE(2) / O7777777777777777 /
c      DATA RIGHT(1) / O1461000000000000 /
c      DATA RIGHT(2) / O0000000000000000 /
c      DATA DIVER(1) / O1451000000000000 /
c      DATA DIVER(2) / O0000000000000000 /
c      DATA LOG10(1) / O1157163034761674 /
c      DATA LOG10(2) / O0006677466732724 /
c
c     MACHINE CONSTANTS FOR THE CONVEX C-120 (NATIVE MODE)
c     WITH OR WITHOUT -R8 OPTION
c
c === MACHINE = CONVEX.C1
c === MACHINE = CONVEX.C1.R8
c      DATA DMACH(1) / 5.562684646268007D-309 /
c      DATA DMACH(2) / 8.988465674311577D+307 /
c      DATA DMACH(3) / 1.110223024625157D-016 /
c      DATA DMACH(4) / 2.220446049250313D-016 /
c      DATA DMACH(5) / 3.010299956639812D-001 /
c
c     MACHINE CONSTANTS FOR THE CONVEX C-120 (IEEE MODE)
c     WITH OR WITHOUT -R8 OPTION
c
c === MACHINE = CONVEX.C1.IEEE
c === MACHINE = CONVEX.C1.IEEE.R8
c      DATA DMACH(1) / 2.225073858507202D-308 /
c      DATA DMACH(2) / 1.797693134862315D+308 /
c      DATA DMACH(3) / 1.110223024625157D-016 /
c      DATA DMACH(4) / 2.220446049250313D-016 /
c      DATA DMACH(5) / 3.010299956639812D-001 /
c
c     MACHINE CONSTANTS FOR THE CYBER 170/180 SERIES USING NOS (FTN5).
c
c === MACHINE = CYBER.170.NOS
c === MACHINE = CYBER.180.NOS
c      DATA SMALL(1) / O"00604000000000000000" /
c      DATA SMALL(2) / O"00000000000000000000" /
c      DATA LARGE(1) / O"37767777777777777777" /
c      DATA LARGE(2) / O"37167777777777777777" /
c      DATA RIGHT(1) / O"15604000000000000000" /
c      DATA RIGHT(2) / O"15000000000000000000" /
c      DATA DIVER(1) / O"15614000000000000000" /
c      DATA DIVER(2) / O"15010000000000000000" /
c      DATA LOG10(1) / O"17164642023241175717" /
c      DATA LOG10(2) / O"16367571421742254654" /
c
c     MACHINE CONSTANTS FOR THE CDC 180 SERIES USING NOS/VE
c
c === MACHINE = CYBER.180.NOS/VE
c      DATA SMALL(1) / Z"3001800000000000" /
c      DATA SMALL(2) / Z"3001000000000000" /
c      DATA LARGE(1) / Z"4FFEFFFFFFFFFFFE" /
c      DATA LARGE(2) / Z"4FFE000000000000" /
c      DATA RIGHT(1) / Z"3FD2800000000000" /
c      DATA RIGHT(2) / Z"3FD2000000000000" /
c      DATA DIVER(1) / Z"3FD3800000000000" /
c      DATA DIVER(2) / Z"3FD3000000000000" /
c      DATA LOG10(1) / Z"3FFF9A209A84FBCF" /
c      DATA LOG10(2) / Z"3FFFF7988F8959AC" /
c
c     MACHINE CONSTANTS FOR THE CYBER 205
c
c === MACHINE = CYBER.205
c      DATA SMALL(1) / X'9000400000000000' /
c      DATA SMALL(2) / X'8FD1000000000000' /
c      DATA LARGE(1) / X'6FFF7FFFFFFFFFFF' /
c      DATA LARGE(2) / X'6FD07FFFFFFFFFFF' /
c      DATA RIGHT(1) / X'FF74400000000000' /
c      DATA RIGHT(2) / X'FF45000000000000' /
c      DATA DIVER(1) / X'FF75400000000000' /
c      DATA DIVER(2) / X'FF46000000000000' /
c      DATA LOG10(1) / X'FFD04D104D427DE7' /
c      DATA LOG10(2) / X'FFA17DE623E2566A' /
c
c     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.
c
c === MACHINE = CDC.6000
c === MACHINE = CDC.7000
c      DATA SMALL(1) / 00604000000000000000B /
c      DATA SMALL(2) / 00000000000000000000B /
c      DATA LARGE(1) / 37767777777777777777B /
c      DATA LARGE(2) / 37167777777777777777B /
c      DATA RIGHT(1) / 15604000000000000000B /
c      DATA RIGHT(2) / 15000000000000000000B /
c      DATA DIVER(1) / 15614000000000000000B /
c      DATA DIVER(2) / 15010000000000000000B /
c      DATA LOG10(1) / 17164642023241175717B /
c      DATA LOG10(2) / 16367571421742254654B /
c
c     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
c
c === MACHINE = CRAY
c      DATA SMALL(1) / 201354000000000000000B /
c      DATA SMALL(2) / 000000000000000000000B /
c      DATA LARGE(1) / 577767777777777777777B /
c      DATA LARGE(2) / 000007777777777777776B /
c      DATA RIGHT(1) / 376434000000000000000B /
c      DATA RIGHT(2) / 000000000000000000000B /
c      DATA DIVER(1) / 376444000000000000000B /
c      DATA DIVER(2) / 000000000000000000000B /
c      DATA LOG10(1) / 377774642023241175717B /
c      DATA LOG10(2) / 000007571421742254654B /
c
c     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
c
c     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING LINE -
c     STATIC DMACH(5)
c
c === MACHINE = DATA_GENERAL.ECLIPSE.S/200
c      DATA SMALL/20K,3*0/,LARGE/77777K,3*177777K/
c      DATA RIGHT/31420K,3*0/,DIVER/32020K,3*0/
c      DATA LOG10/40423K,42023K,50237K,74776K/
c
c     ELXSI 6400
c
c === MACHINE = ELSXI.6400
c      DATA SMALL(1), SMALL(2) / '00100000'X,'00000000'X /
c      DATA LARGE(1), LARGE(2) / '7FEFFFFF'X,'FFFFFFFF'X /
c      DATA RIGHT(1), RIGHT(2) / '3CB00000'X,'00000000'X /
c      DATA DIVER(1), DIVER(2) / '3CC00000'X,'00000000'X /
c      DATA LOG10(1), DIVER(2) / '3FD34413'X,'509F79FF'X /
c
c     MACHINE CONSTANTS FOR THE HARRIS 220
c     MACHINE CONSTANTS FOR THE HARRIS SLASH 6 AND SLASH 7
c
c === MACHINE = HARRIS.220
c === MACHINE = HARRIS.SLASH6
c === MACHINE = HARRIS.SLASH7
c      DATA SMALL(1),SMALL(2) / '20000000, '00000201 /
c      DATA LARGE(1),LARGE(2) / '37777777, '37777577 /
c      DATA RIGHT(1),RIGHT(2) / '20000000, '00000333 /
c      DATA DIVER(1),DIVER(2) / '20000000, '00000334 /
c      DATA LOG10(1),LOG10(2) / '23210115, '10237777 /
c
c     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES.
c     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
c
c === MACHINE = HONEYWELL.600/6000
c === MACHINE = HONEYWELL.DPS.8/70
c      DATA SMALL(1),SMALL(2) / O402400000000, O000000000000 /
c      DATA LARGE(1),LARGE(2) / O376777777777, O777777777777 /
c      DATA RIGHT(1),RIGHT(2) / O604400000000, O000000000000 /
c      DATA DIVER(1),DIVER(2) / O606400000000, O000000000000 /
c      DATA LOG10(1),LOG10(2) / O776464202324, O117571775714 /
c
c      MACHINE CONSTANTS FOR THE HP 2100
c      3 WORD DOUBLE PRECISION OPTION WITH FTN4
c
c === MACHINE = HP.2100.3_WORD_DP
c      DATA SMALL(1), SMALL(2), SMALL(3) / 40000B,       0,       1 /
c      DATA LARGE(1), LARGE(2), LARGE(3) / 77777B, 177777B, 177776B /
c      DATA RIGHT(1), RIGHT(2), RIGHT(3) / 40000B,       0,    265B /
c      DATA DIVER(1), DIVER(2), DIVER(3) / 40000B,       0,    276B /
c      DATA LOG10(1), LOG10(2), LOG10(3) / 46420B,  46502B,  77777B /
c
c      MACHINE CONSTANTS FOR THE HP 2100
c      4 WORD DOUBLE PRECISION OPTION WITH FTN4
c
c === MACHINE = HP.2100.4_WORD_DP
c      DATA SMALL(1), SMALL(2) /  40000B,       0 /
c      DATA SMALL(3), SMALL(4) /       0,       1 /
c      DATA LARGE(1), LARGE(2) /  77777B, 177777B /
c      DATA LARGE(3), LARGE(4) / 177777B, 177776B /
c      DATA RIGHT(1), RIGHT(2) /  40000B,       0 /
c      DATA RIGHT(3), RIGHT(4) /       0,    225B /
c      DATA DIVER(1), DIVER(2) /  40000B,       0 /
c      DATA DIVER(3), DIVER(4) /       0,    227B /
c      DATA LOG10(1), LOG10(2) /  46420B,  46502B /
c      DATA LOG10(3), LOG10(4) /  76747B, 176377B /
c
c     HP 9000
c
c      D1MACH(1) = 2.8480954D-306
c      D1MACH(2) = 1.40444776D+306
c      D1MACH(3) = 2.22044605D-16
c      D1MACH(4) = 4.44089210D-16
c      D1MACH(5) = 3.01029996D-1
c
c === MACHINE = HP.9000
c      DATA SMALL(1), SMALL(2) / 00040000000B, 00000000000B /
c      DATA LARGE(1), LARGE(2) / 17737777777B, 37777777777B /
c      DATA RIGHT(1), RIGHT(2) / 07454000000B, 00000000000B /
c      DATA DIVER(1), DIVER(2) / 07460000000B, 00000000000B /
c      DATA LOG10(1), LOG10(2) / 07764642023B, 12047674777B /
c
c     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
c     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
c     THE INTERDATA 3230 AND INTERDATA 7/32.
c
c === MACHINE = IBM.360
c === MACHINE = IBM.370
c === MACHINE = XEROX.SIGMA.5
c === MACHINE = XEROX.SIGMA.7
c === MACHINE = XEROX.SIGMA.9
c === MACHINE = SEL.85
c === MACHINE = SEL.86
c === MACHINE = INTERDATA.3230
c === MACHINE = INTERDATA.7/32
c      DATA SMALL(1),SMALL(2) / Z00100000, Z00000000 /
c      DATA LARGE(1),LARGE(2) / Z7FFFFFFF, ZFFFFFFFF /
c      DATA RIGHT(1),RIGHT(2) / Z33100000, Z00000000 /
c      DATA DIVER(1),DIVER(2) / Z34100000, Z00000000 /
c      DATA LOG10(1),LOG10(2) / Z41134413, Z509F79FF /
c
c     MACHINE CONSTANTS FOR THE INTERDATA 8/32
c     WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.
c
c     FOR THE INTERDATA FORTRAN VII COMPILER REPLACE
c     THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S.
c
c === MACHINE = INTERDATA.8/32.UNIX
c      DATA SMALL(1),SMALL(2) / Z'00100000', Z'00000000' /
c      DATA LARGE(1),LARGE(2) / Z'7EFFFFFF', Z'FFFFFFFF' /
c      DATA RIGHT(1),RIGHT(2) / Z'33100000', Z'00000000' /
c      DATA DIVER(1),DIVER(2) / Z'34100000', Z'00000000' /
c      DATA LOG10(1),LOG10(2) / Z'41134413', Z'509F79FF' /
c
c     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
c
c === MACHINE = PDP-10.KA
c      DATA SMALL(1),SMALL(2) / "033400000000, "000000000000 /
c      DATA LARGE(1),LARGE(2) / "377777777777, "344777777777 /
c      DATA RIGHT(1),RIGHT(2) / "113400000000, "000000000000 /
c      DATA DIVER(1),DIVER(2) / "114400000000, "000000000000 /
c      DATA LOG10(1),LOG10(2) / "177464202324, "144117571776 /
c
c     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
c
c === MACHINE = PDP-10.KI
c      DATA SMALL(1),SMALL(2) / "000400000000, "000000000000 /
c      DATA LARGE(1),LARGE(2) / "377777777777, "377777777777 /
c      DATA RIGHT(1),RIGHT(2) / "103400000000, "000000000000 /
c      DATA DIVER(1),DIVER(2) / "104400000000, "000000000000 /
c      DATA LOG10(1),LOG10(2) / "177464202324, "047674776746 /
c
c     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
c     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
c
c === MACHINE = PDP-11.32-BIT
c      DATA SMALL(1),SMALL(2) /    8388608,           0 /
c      DATA LARGE(1),LARGE(2) / 2147483647,          -1 /
c      DATA RIGHT(1),RIGHT(2) /  612368384,           0 /
c      DATA DIVER(1),DIVER(2) /  620756992,           0 /
c      DATA LOG10(1),LOG10(2) / 1067065498, -2063872008 /
c
c      DATA SMALL(1),SMALL(2) / O00040000000, O00000000000 /
c      DATA LARGE(1),LARGE(2) / O17777777777, O37777777777 /
c      DATA RIGHT(1),RIGHT(2) / O04440000000, O00000000000 /
c      DATA DIVER(1),DIVER(2) / O04500000000, O00000000000 /
c      DATA LOG10(1),LOG10(2) / O07746420232, O20476747770 /
c
c     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
c     16-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
c
c === MACHINE = PDP-11.16-BIT
c      DATA SMALL(1),SMALL(2) /    128,      0 /
c      DATA SMALL(3),SMALL(4) /      0,      0 /
c      DATA LARGE(1),LARGE(2) /  32767,     -1 /
c      DATA LARGE(3),LARGE(4) /     -1,     -1 /
c      DATA RIGHT(1),RIGHT(2) /   9344,      0 /
c      DATA RIGHT(3),RIGHT(4) /      0,      0 /
c      DATA DIVER(1),DIVER(2) /   9472,      0 /
c      DATA DIVER(3),DIVER(4) /      0,      0 /
c      DATA LOG10(1),LOG10(2) /  16282,   8346 /
c      DATA LOG10(3),LOG10(4) / -31493, -12296 /
c
c      DATA SMALL(1),SMALL(2) / O000200, O000000 /
c      DATA SMALL(3),SMALL(4) / O000000, O000000 /
c      DATA LARGE(1),LARGE(2) / O077777, O177777 /
c      DATA LARGE(3),LARGE(4) / O177777, O177777 /
c      DATA RIGHT(1),RIGHT(2) / O022200, O000000 /
c      DATA RIGHT(3),RIGHT(4) / O000000, O000000 /
c      DATA DIVER(1),DIVER(2) / O022400, O000000 /
c      DATA DIVER(3),DIVER(4) / O000000, O000000 /
c      DATA LOG10(1),LOG10(2) / O037632, O020232 /
c      DATA LOG10(3),LOG10(4) / O102373, O147770 /
c
c     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000
c
c === MACHINE = SEQUENT.BALANCE.8000
c      DATA SMALL(1),SMALL(2) / $00000000,  $00100000 /
c      DATA LARGE(1),LARGE(2) / $FFFFFFFF,  $7FEFFFFF /
c      DATA RIGHT(1),RIGHT(2) / $00000000,  $3CA00000 /
c      DATA DIVER(1),DIVER(2) / $00000000,  $3CB00000 /
c      DATA LOG10(1),LOG10(2) / $509F79FF,  $3FD34413 /
c
c     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES. FTN COMPILER
c
c === MACHINE = UNIVAC.1100
c      DATA SMALL(1),SMALL(2) / O000040000000, O000000000000 /
c      DATA LARGE(1),LARGE(2) / O377777777777, O777777777777 /
c      DATA RIGHT(1),RIGHT(2) / O170540000000, O000000000000 /
c      DATA DIVER(1),DIVER(2) / O170640000000, O000000000000 /
c      DATA LOG10(1),LOG10(2) / O177746420232, O411757177572 /
c
c     MACHINE CONSTANTS FOR VAX 11/780
c     (EXPRESSED IN INTEGER AND HEXADECIMAL)
c    *** THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS***
c
c === MACHINE = VAX.11/780
c      DATA SMALL(1), SMALL(2) /        128,           0 /
c      DATA LARGE(1), LARGE(2) /     -32769,          -1 /
c      DATA RIGHT(1), RIGHT(2) /       9344,           0 /
c      DATA DIVER(1), DIVER(2) /       9472,           0 /
c      DATA LOG10(1), LOG10(2) /  546979738,  -805796613 /
c
c    ***THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSYEMS***
c      DATA SMALL(1), SMALL(2) / Z00000080, Z00000000 /
c      DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
c      DATA RIGHT(1), RIGHT(2) / Z00002480, Z00000000 /
c      DATA DIVER(1), DIVER(2) / Z00002500, Z00000000 /
c      DATA LOG10(1), LOG10(2) / Z209A3F9A, ZCFF884FB /
c
c   MACHINE CONSTANTS FOR VAX 11/780 (G-FLOATING)
c     (EXPRESSED IN INTEGER AND HEXADECIMAL)
c    *** THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS***
c
c      DATA SMALL(1), SMALL(2) /         16,           0 /
c      DATA LARGE(1), LARGE(2) /     -32769,          -1 /
c      DATA RIGHT(1), RIGHT(2) /      15552,           0 /
c      DATA DIVER(1), DIVER(2) /      15568,           0 /
c      DATA LOG10(1), LOG10(2) /  1142112243, 2046775455 /
c
c    ***THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSYEMS***
c      DATA SMALL(1), SMALL(2) / Z00000010, Z00000000 /
c      DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
c      DATA RIGHT(1), RIGHT(2) / Z00003CC0, Z00000000 /
c      DATA DIVER(1), DIVER(2) / Z00003CD0, Z00000000 /
c      DATA LOG10(1), LOG10(2) / Z44133FF3, Z79FF509F /
c
      if ( i .lt. 1  .or.  5 .lt. i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'D1MACH - Fatal error!'
        write ( *, '(a)' ) '  I OUT OF BOUNDS'
        stop
      end if

      d1mach = dmach(i)

      return
      end
      subroutine triex ( f, ver, epsabs, epsrel, limev, result, abserr,
     & ier )

c*********************************************************************72
c
cc TRIEX approximates the integral of a function over a triangle.
c
c  Discussion:
c
c    The routine calculates an approximation  RESULT  to a
c    double integral   I = integral of  F  over the triangle
c    with vertices (VER(1,K),VER(2,K)) for K=1,2,3,
c    hopefully satisfying following claim for accuracy
c    ABS(I-RESULT) .LE. MAX(EPSABS,EPSREL*ABS(I)).
c
c  Modified:
c
c    06 December 2006
c
c  Author:
c
c    Elise deDoncker, Ian Robinson
c
c  Reference:
c
c    Elise deDoncker, Ian Robinson,
c    Algorithm 612:
c    Integration over a Triangle Using Nonlinear Extrapolation,
c    ACM Transactions on Mathematical Software,
c    Volume 10, Number 1, March 1984, pages 17-22.
c
c  Parameters:
c
c            F      - DOUBLE PRECISION
c                     FUNCTION DEFINING THE INTEGRAND F(X,Y)
c                     THE ACTUAL NAME FOR F NEEDS TO BE DECLARED
c                     E X T E R N A L IN THE DRIVER PROGRAM.
c
c            VER    - DOUBLE PRECISION
c                     ARRAY OF DIMENSION (2,3), WHERE VER(1,K) AND
c                     VER(2,K) REPRESENT THE ABSCISSAE AND THE
c                     ORDINATES RESPECTIVELY OF THE K-TH VERTEX
c
c            EPSABS - DOUBLE PRECISION
c                     ABSOLUTE ACCURACY REQUESTED
c
c            EPSREL - DOUBLE PRECISION
c                     RELATIVE ACCURACY REQUESTED
c                     IF EPSABS .LE. 0 AND EPSREL .LE. 0,
c                     THE ROUTINE WILL END WITH IER .NE. 0.
c
c            LIMEV  - LIMIT ON THE NUMBER OF INTEGRAND EVALUATIONS
c                     IF MORE THAN XXX EVALUATIONS ARE TO BE ALLOWED,
c                     DIMENSION STATEMENTS AND  DATA N  IN TRIEX
c                     (WHERE N IS THE MAXIMUM NUMBER OF TRIANGLES
c                     ALLOWED) NEED TO BE ADJUSTED ACCORDING TO
c                     XXX .LE. (N-1)*178/3 + 46,
c                     I.E. PRESENT LIMIT N = 1000 IMPLIES
c                     LIMEV .LE. 59320.
c
c            RESULT - DOUBLE PRECISION
c                     APPROXIMATION TO THE INTEGRAL I
c
c            ABSERR - DOUBLE PRECISION
c                     ESTIMATE OF THE ABSOLUTE ERROR ABS(I-RESULT)
c
c            IER    - IER  = 0 NORMAL AND RELIABLE TERMINATION OF
c                              THE ROUTINE. IT IS ASSUMED THAT THE
c                              REQUESTED ACCURACY HAS BEEN OBTAINED.
c                   - IER .NE. 0 ABNORMAL TERMINATION OF THE ROUTINE
c                              IT IS ASSUMED THAT THE REQUESTED
c                              ACCURACY HAS NOT BEEN OBTAINED.
c                          = 1 MAXIMUM NUMBER OF FUNCTION EVALUATIONS
c                              ALLOWED HAS BEEN REACHED.
c                          = 2 THE PRESENCE OF ROUNDOFF ERROR HAS BEEN
c                              DETECTED, WHICH PREVENTS THE REQUESTED
c                              ACCURACY FROM BEING ACHIEVED.
c                          = 3 IT IS PRESUMED THAT THE TOLERANCE CANNOT
c                              BE ACHIEVED WITHIN THE SPECIFIED NUMBER
c                              OF INTEGRAND EVALUATIONS, AND THAT THE
c                              RESULT RETURNED IS THE BEST WHICH CAN BE
c                              OBTAINED WITHIN THAT LIMIT.
c
      implicit none

      double precision abserr, abser1, area, area4, bndec, dovflo,
     * drelpr, dres, dresc, dres0, dunflo, d1mach, elist, emach, eprn,
     * epsabs, epsrel, erlarg, erlast, errbnd, errmax, erro4, errsum,
     * f, ofrn, q, rcopy, resex, result, resul1, resla, rex, rlast,
     * rlist, stor, ufrn, ver, x, y
      double precision dabs, dmax1, dmin1
      integer i, icnt, id, ier, iord, j, j1, j2, j3, jpoint, jroff,
     * jslast, jstage, jvrtx, k, kbiger, last, limev, limit, l4,
     * maxerr, maxer0, n, nrmax, ntr, numrl2, nvrtx, n2
      integer max0, min0
      logical extrap
      dimension elist(1000), iord(1000), jpoint(4), jstage(1000),
     * jvrtx(1000,3), rcopy(52), rex(1000), rlist(1000), ver(2,3),
     * x(3000), y(3000)

      common /mach/ emach, eprn, ufrn, ofrn
c
c            DIMENSIONS
c            ----------
c            N  IS THE MAXIMUM NUMBER OF TRIANGLES ALLOWED.
c           INCREASING N REQUIRES INCREASING DIMENSIONS IN
c           TRIEX (AND IN SUBROUTINE DIVIN CORRESPONDINGLY)
c           DIMENSIONS OF ELIST(N),IORD(N),JSTAGE(N),JVRTX(N,3),
c           REX(N),RLIST(N),X(3*N),Y(3*N)
c           SHOULD BE SUCH THAT
c           LIMEV .LE. (N-1)*178/3 + 46
      external f
      data n /1000/
c
c            RCOPY  IS OF DIMENSION AT LEAST (LIMEXP+2)
c           WITH  DATA LIMEXP  GIVEN IN SUBROUTINE EPSALG
c
c
c            MACHINE DEPENDENT CONSTANTS
c            ---------------------------
c            DRELPR  IS THE LARGEST RELATIVE SPACING
c            DUNFLO  IS THE SMALLEST POSITIVE MAGNITUDE
c            DOVFLO  IS THE LARGEST MAGNITUDE
c
      drelpr = d1mach(4)
      dunflo = d1mach(1)
      dovflo = d1mach(2)
c
c            LIST OF MAJOR VARIABLES
c            -----------------------
c
c           X         - LIST OF ABSCISSAE OF ALL VERTICES
c           Y         - LIST OF ORDINATES OF ALL VERTICES
c           JVRTX     - ARRAY OF NUMBERS OF ALL VERTICES
c           RLIST     - LIST OF DEGREE-11 INTEGRAL APPROXIMATIONS FOR
c                       ALL SUBTRIANGLES
c           REX       - LIST OF DEGREE-9 INTEGRAL APPROXIMATIONS FOR
c                       ALL SUBTRIANGLES
c           ELIST     - LIST OF ERROR ESTIMATES FOR ALL SUBTRIANGLES
c           JSTAGE    - LIST OF SUBDIVISION LEVELS FOR ALL
c                       SUBTRIANGLES
c           AREA      - CURRENT SUM OF DEGREE-11 RESULTS
c           ERRSUM    - CURRENT SUM OF ERROR ESTIMATES
c           RESEX     - APPROXIMATION TO THE TOTAL INTEGRAL COMPUTED
c                       AS A SUM OF DEGREE-9 RESULTS OVER THE SMALLEST
c                       SUBTRIANGLES (HIGHEST LEVEL OF SUBDIVISION) AND
c                       OF DEGREE-11 RESULTS OVER THE LARGER TRIANGLES
c           ERRBND    - ESTIMATE OF THE GLOBAL ABSOLUTE TOLERANCE
c           BNDEC     - GLOBAL TOLERANCE OVER THE SET OF LARGER
c                       TRIANGLES
c           ERLARG    - SUM OF ERROR ESTIMATES OVER THE CURRENT SET OF
c                       LARGE TRIANGLES
c           NTR       - CURRENT NUMBER OF TRIANGLES IN THE PARTITION
c                       OF THE ORIGINAL INTEGRATION REGION
c           MAXERR    - POINTER TO SUBTRIANGLE WITH LARGEST ERROR
c                       ESTIMATE
c           RCOPY     - LIST OF THE ENTRIES TO THE EXTRAPOLATION
c                       PROCEDURE
c           NUMRL2    - NUMBER OF ELEMENTS IN RCOPY
c           EXTRAP    - LOGICAL VARIABLE DENOTING THAT THE ROUTINE
c                       IS ATTEMPTING TO DECREASE THE VALUE OF ERLARG
c                       THROUGH SUBDIVISION OF LARGE (LOW-LEVEL)
c                       TRIANGLES, BEFORE A NEW EXTRAPOLATION IS
c                       CARRIED OUT
c
c  First approximation to the integral
c
      ier = 0
      i = 0
      k = 1
      eprn = 0.5d+01*drelpr
      ufrn = dunflo/eprn
      stor = 0.0d+00
      do j=1,3
        x(j) = ver(1,j)
        y(j) = ver(2,j)
        jvrtx(1,j) = j
      end do
      call quatri(f, x(1), y(1), x(2), y(2), x(3), y(3), resex, result,
     & abserr, dresc, i, k, stor)
      dres = dabs(resex)
      errbnd = dmax1(epsabs,epsrel*dres,eprn*dres)
      if (abserr.le.errbnd) go to 170
c
c  Initialization
c
      rlist(1) = result
      rcopy(1) = resex
      elist(1) = abserr
      ofrn = dovflo
      jstage(1) = 1
      area = result
      errsum = abserr
      errmax = abserr
      emach = drelpr
      maxerr = 1
      ntr = 1
      nrmax = 1
      nvrtx = 3
      numrl2 = 2
      n2 = 2
      jroff = 0
      extrap = .false.
c
c  LIMIT = number of triangles allowed in final partition of
c  original triangle
c
      limit = max0(4,min0(n,1+((limev-46)/178)*3))
      l4 = (limit-8)/4
c
c  Main DO loop
c
      do 130 last=2,limit,3
c
c  Subdivide the triangle with NRMAX-th largest error estimate
c  and integrate over its subtriangles
c
        call divin(maxerr, x, y, jvrtx, ntr, nvrtx, jpoint)
        maxer0 = maxerr
        rlast = rlist(maxerr)
        jslast = jstage(maxerr)
        erlast = errmax
        area4 = 0.0d+00
        erro4 = 0.0d+00
        k = 1
        icnt = 0
        stor = 0.0d+00

        do i=1,4
          j = jpoint(i)
          j1 = jvrtx(j,1)
          j2 = jvrtx(j,2)
          j3 = jvrtx(j,3)
          call quatri(f, x(j1), y(j1), x(j2), y(j2), x(j3), y(j3),
     &     rex(j), rlist(j), elist(j), dresc, i, k, stor)
          area4 = area4 + rlist(j)
          erro4 = erro4 + elist(j)
          jstage(j) = jslast + 1
          if (elist(j).eq.dresc) then
            icnt = 1
          end if
        end do
c
c  Improve previous integral approximation and error
c  estimate and test for accuracy
c
        area = area + area4 - rlast
        errsum = errsum + erro4 - errmax
        if (errsum.le.errbnd) go to 150
c
c  Test for roundoff error and eventually set error flag
c
        if (icnt.eq.1) go to 30
        q = 0.0d+00
        if (area4.ne.0.0d+00) q = rlast/area4
        if (0.9999d+00.lt.q .and. q.lt.1.0001d+00 .and.
     &   errmax.le.erro4) jroff = jroff + 1
        if (jroff.ge.10) ier = 2
c
c  Set error flag in the case that maximum number of
c  subtriangles is reached
c
   30   if (last.ge.limit-2) ier = 1
        if (ier.ne.0) go to 140
c
c  Call subroutine ORDER to maintain the descending
c  ordering in the list of error estimates, and to select
c  the subtriangle with NRMAX-th largest error estimate
c  (to be subdivided next)
c
        call order(limit, last, maxerr, errmax, elist, iord, nrmax, ntr)
        if (last.eq.2) go to 120
c
c  Check whether the triangles obtained in the foregoing
c  subdivision are large, and adjust ERLARG if so
c
        erlarg = erlarg - erlast
        if (jstage(maxer0).le.numrl2) erlarg = erlarg + erro4
        if (extrap) go to 40
c
c  Check whether the following subdivision would provide
c  large triangles, and proceed immediately to this
c  subdivision if so
c
        if (jstage(maxerr).le.numrl2) go to 130
c
c  Decrease the integration error over the collection of
c  large triangles
c
        extrap = .true.
        nrmax = 2
   40   bndec = dmax1(eprn*dres,dmin1(0.1d+00*errbnd,0.1d-02*dres))
        if (erlarg.le.bndec) go to 60
        id = nrmax
        do j=id,ntr
          maxerr = iord(nrmax)
          errmax = elist(maxerr)
          if (jstage(maxerr).le.numrl2) go to 130
          nrmax = nrmax + 1
        end do
c
c  Compute next entry for the extrapolation procedure
c
   60   resex = 0.0d+00
        do j=1,ntr
          if (jstage(j).gt.numrl2) resex = resex + rex(j)
          if (jstage(j).le.numrl2) resex = resex + rlist(j)
        end do
        dres = dabs(resex)
        errbnd = dmax1(epsabs,epsrel*dres,eprn*dres)
c
c  Extrapolate and test for accuracy
c
        numrl2 = numrl2 + 1
        n2 = n2 + 1
        rcopy(numrl2) = resex
        call epsalg(numrl2, rcopy, resul1, abser1, resla)
        abser1 = abser1 + erlarg
        if (numrl2.eq.n2 .or. n2.eq.3) go to 80
        if (abser1.gt.abserr) go to 90
   80   result = resul1
        abserr = abser1
        dres0 = dabs(result)
        if (abserr.le.dmax1(epsabs,epsrel*dres0,eprn*dres0) .and.
     &   n2.gt.3) go to 170
c
c  In the case that the extrapolation error estimate is
c  considerably smaller than the error sum over all
c  subtriangles, check if enough steps are left over
c  in order to complete another extrapolation
c  stage. Set error flag for abnormal termination
c  if necessary
c
   90   if (abserr*dabs(resex).ge.0.5d-01*errsum*dres0 .or. last.le.l4)
     &   go to 110
        kbiger = 0
        do j=1,ntr
          if (elist(j).gt.bndec) kbiger = kbiger + 1
        end do
        if (kbiger.ge.(limit-last-2)/3) ier = 3
  110   if (ier.ne.0) go to 140
c
c  Start new extrapolation stage, by integrating over
c  subtriangle with largest error estimate
c
        maxerr = iord(1)
        errmax = elist(maxerr)
        nrmax = 1
        erlarg = errsum
        go to 130
c
c  Second extrapolation entry
c
  120   erlarg = errsum
        resex = rex(1) + rex(2) + rex(3) + rex(4)
        rcopy(2) = resex
        dres = dabs(resex)
        resla = resex
  130 continue
c
c  In the case that IER .NE. 0, select the best approximation
c
c  (RESULT (extrapolated), or AREA (sum of integrals over
c  subtriangles) )
c
  140 dres = dabs(area)
      if (dres*abserr.le.dres0*errsum) go to 170
      if (dabs(area-resex)*dres0.gt.dabs(result-resex)*dres) go to 170
  150 result = 0.0d+00
      do j=1,ntr
        result = result + rlist(j)
      end do
      abserr = errsum
  170 continue

      return
      end
      subroutine divin ( max, x, y, ivrtx, ntr, nvrtx, ipoint )

c*********************************************************************72
c
cc DIVIN subdivides a triangle into four subtriangles.
c
c  Discussion:
c
c    The routine subdivides the given triangle with index MAX
c    into 4 similar subtriangles, by halving the sides of the
c    given triangle. it calculates the coordinates of the new
c    vertices, and defines the subtriangles by defining pointers
c    to these triangles and to their vertices.
c
c  Modified:
c
c    06 December 2006
c
c  Author:
c
c    Elise deDoncker, Ian Robinson
c
c  Reference:
c
c    Elise deDoncker, Ian Robinson,
c    Algorithm 612:
c    Integration over a Triangle Using Nonlinear Extrapolation,
c    ACM Transactions on Mathematical Software,
c    Volume 10, Number 1, March 1984, pages 17-22.
c
c  Parameters:
c
c           MAX     - INTEGER
c                     INDEX OF THE TRIANGLE WHICH MUST BE SUBDIVIDED
c
c           X       - DOUBLE PRECISION
c                     ABSCISSAE OF THE VERTICES OF ALL TRIANGLES IN THE
c                     PARTITION OF THE INTEGRATION DOMAIN
c
c           Y       - DOUBLE PRECISION
c                     ORDINATES OF THE VERTICES OF ALL TRIANGLES IN THE
c                     PARTITION OF THE INTEGRATION DOMAIN
c
c           IVRTX   - INTEGER
c                     IVRTX(J,K) POINTS TO THE K-TH VERTEX (K=1,2,3)
c                     OF THE TRIANGLE WITH INDEX J
c
c           NTR     - INTEGER
c                     TOTAL NUMBER OF TRIANGLES IN THE PARTITION
c                     NTR IS INCREASED BY 3 DURING THE COURSE OF THIS
c                     SUBROUTINE
c
c           NVRTX   - INTEGER
c                     TOTAL NUMBER OF VERTICES IN THE PARTITION
c                     NVRTX IS INCREASED BY 3 DURING THIS CALL.
c
c           IPOINT  - INTEGER
c                     POINTERS TO THE 4 NEW SUBTRIANGLES, SUCH THAT
c                     IPOINT(1)=MAX, IPOINT(2)=NTR+1, IPOINT(3)=NTR+2,
c                     IPOINT(4)=NTR+3
c
      implicit none

      double precision x, y
      integer i, i1, i2, i3, i4, i5, i6, ip, ipoint, ivrtx, max, ntr,
     * nvrtx

      dimension ipoint(4), ivrtx(1000,3), x(3000), y(3000)
c
c           LIST OF MAJOR VARIABLES
c           -----------------------
c           I1,...,I6 - INDICES OF VERTICES OF THE 4 NEW SUBTRIANGLES
c
c
c  Set pointers to the 4 new subtriangles
c
      do i=1,3
        ip = ntr + i
        ipoint(i+1) = ip
      end do
      ipoint(1) = max
c
c  Define vertices of the subtriangles
c
      i1 = ivrtx(max,1)
      i2 = ivrtx(max,2)
      i3 = ivrtx(max,3)
      i4 = nvrtx + 1
      i5 = nvrtx + 2
      i6 = nvrtx + 3
      x(i4) = (x(i2)+x(i3))*0.5d+00
      y(i4) = (y(i2)+y(i3))*0.5d+00
      x(i5) = (x(i3)+x(i1))*0.5d+00
      y(i5) = (y(i3)+y(i1))*0.5d+00
      x(i6) = (x(i1)+x(i2))*0.5d+00
      y(i6) = (y(i1)+y(i2))*0.5d+00
c
c  Adjust NVRTX
c
      NVRTX = I6
c
c  The subdivision is organized as is indicated in the
c  pictures below. the pointers to the vertices of the
c  subtriangles are set in such a way that
c  subdivision of the triangle with index MAX produces the
c  triangles with indices MAX, NTR+1, NTR+2, NTR+3
c  (determined by (I1,I6,I5), (I6,I2,I4), (I5,I4,I3) and
c  (I4,I6,I5) respectively).
c
c                       I3
c
c
c                       I5       I4
c
c
c                       I1       I6       I2
c
      ivrtx(max,1) = i1
      ivrtx(max,2) = i6
      ivrtx(max,3) = i5
      ntr = ntr + 1
      ivrtx(ntr,1) = i6
      ivrtx(ntr,2) = i2
      ivrtx(ntr,3) = i4
      ntr = ntr + 1
      ivrtx(ntr,1) = i5
      ivrtx(ntr,2) = i4
      ivrtx(ntr,3) = i3
      ntr = ntr + 1
      ivrtx(ntr,1) = i4
      ivrtx(ntr,2) = i6
      ivrtx(ntr,3) = i5

      return
      end
      subroutine quatri ( f, u1, v1, u2, v2, u3, v3, res9, res11, est,
     * dresc, i, k, stor )

c*********************************************************************72
c
cc QUATRI returns a quadrature rule for a triangle.
c
c  Discussion:
c
c    To compute - I(F) = integral of F over the triangle
c    with vertices (U1,V1),(U2,V2),(U3,V3), and error
c    estimate, - integral of abs(F) over this triangle
c
c  Modified:
c
c    06 December 2006
c
c  Author:
c
c    Elise deDoncker, Ian Robinson
c
c  Reference:
c
c    Elise deDoncker, Ian Robinson,
c    Algorithm 612:
c    Integration over a Triangle Using Nonlinear Extrapolation,
c    ACM Transactions on Mathematical Software,
c    Volume 10, Number 1, March 1984, pages 17-22.
c
c  PARAMETERS:
c
c           F       - DOUBLE PRECISION
c                     FUNCTION SUBPROGRAM DEFINING THE INTEGRAND
c                     F(X,Y)  THE ACTUAL NAME FOR F NEEDS TO BE
c                     DECLARED E X T E R N A L IN THE CALLING
c                     PROGRAM.
c
c           U1,U2,U3- DOUBLE PRECISION
c                     ABSCISSAE OF VERTICES
c
c           V1,V2,V3- DOUBLE PRECISION
c                     ORDINATES OF VERTICES
c
c           RES9    - DOUBLE PRECISION
c                     APPROXIMATION TO THE INTEGRAL IF, OBTAINED BY THE
c                     LYNESS AND JESPERSEN RULE OF DEGREE 9, USING
c                     19 POINTS
c
c           RES11   - DOUBLE PRECISION
c                     APPROXIMATION TO THE INTEGRAL IF, OBTAINED BY THE
c                     LYNESS AND JESPERSEN RULE OF DEGREE 11,
c                     USING 28 POINTS
c
c           EST     - DOUBLE PRECISION
c                     ESTIMATE OF THE ABSOLUTE ERROR
c
c           DRESC   - DOUBLE PRECISION
c                     APPROXIMATION TO THE INTEGRAL OF ABS(F- IF/DJ),
c                     OBTAINED BY THE RULE OF DEGREE 9, AND USED FOR
c                     THE COMPUTATION OF EST
c
c           I       - INTEGER
c                     NUMBER OF TRIANGLES IN CURRENT SUBDIVISION
c                     (1 .LE. I .LE. 4)
c
c           K       - INTEGER
c                     INDEX OF POINTS TO BE SHARED WITH 4TH TRIANGLE
c                     OF SUBDIVISION
c
c           STOR    - DOUBLE PRECISION
c                     SUM OF FUNCTION VALUES AT THESE POINTS
c
      implicit none

      double precision dj, df0, dresc, emach, eprn, est, f, fstor, fv,
     & f0, ofrn, resab9, res9, res11, stor, u1, u2, u3, ufrn, v1, v2,
     & v3, w, w90, w110, x, y, zeta1, zeta2, z1, z2, z3
      double precision dmax1, dmin1
      integer i, j, jend, k, kount, l

      dimension fstor(3), fv(19), w(15), x(3), y(3), zeta1(15),
     & zeta2(15)

      common /mach/ emach, eprn, ufrn, ofrn
c
c  First homogeneous coordinates of points in degree-9
c  and degree-11 formula, taken with multiplicity 3
c
      data zeta1(1), zeta1(2), zeta1(3), zeta1(4), zeta1(5), zeta1(6),
     & zeta1(7), zeta1(8), zeta1(9), zeta1(10), zeta1(11), zeta1(12),
     & zeta1(13), zeta1(14), zeta1(15) /0.2063496160252593d-01,
     & 0.1258208170141290d+00,0.6235929287619356d+00,
     & 0.9105409732110941d+00,0.3683841205473626d-01,
     & 0.7411985987844980d+00,0.9480217181434233d+00,
     & 0.8114249947041546d+00,0.1072644996557060d-01,
     & 0.5853132347709715d+00,0.1221843885990187d+00,
     & 0.4484167758913055d-01,0.6779376548825902d+00,0.0d+00,
     & 0.8588702812826364d+00/
c
c  Second homogeneous coordinates of points in degree-9
c  and degree-11 formula, taken with multiplicity 3
c
      data zeta2(1), zeta2(2), zeta2(3), zeta2(4), zeta2(5), zeta2(6),
     & zeta2(7), zeta2(8), zeta2(9), zeta2(10), zeta2(11), zeta2(12),
     & zeta2(13), zeta2(14), zeta2(15) /0.4896825191987370d+00,
     & 0.4370895914929355d+00,0.1882035356190322d+00,
     & 0.4472951339445297d-01,0.7411985987844980d+00,
     & 0.3683841205473626d-01,0.2598914092828833d-01,
     & 0.9428750264792270d-01,0.4946367750172147d+00,
     & 0.2073433826145142d+00,0.4389078057004907d+00,
     & 0.6779376548825902d+00,0.4484167758913055d-01,
     & 0.8588702812826364d+00,0.0d+00/
c
c  Weights of midpoint of triangle in degree-9
c  resp. degree-11 formulae
c
      data w90 /0.9713579628279610d-01/
      data w110 /0.8797730116222190d-01/
c
c  Weights in degree-9 and degree-11 rule
c
      data w(1), w(2), w(3), w(4), w(5), w(6), w(7), w(8), w(9), w(10),
     & w(11), w(12), w(13), w(14), w(15) /0.3133470022713983d-01,
     & 0.7782754100477543d-01,0.7964773892720910d-01,
     & 0.2557767565869810d-01,0.4328353937728940d-01,
     & 0.4328353937728940d-01,0.8744311553736190d-02,
     & 0.3808157199393533d-01,0.1885544805613125d-01,
     & 0.7215969754474100d-01,0.6932913870553720d-01,
     & 0.4105631542928860d-01,0.4105631542928860d-01,
     & 0.7362383783300573d-02,0.7362383783300573d-02/
c
c  LIST OF MAJOR VARIABLES
c
c           DJ      - AREA OF THE TRIANGLE
c           DRESC   - APPROXIMATION TO INTEGRAL OF
c                     ABS(F- IF/DJ)  OVER THE TRIANGLE
c           RESAB9  - APPROXIMATION TO INTEGRAL OF
c                     ABS(F) OVER THE TRIANGLE
c           X       - CARTESIAN ABSCISSAE OF THE INTEGRATION
c                     POINTS
c           Y       - CARTESIAN ORDINATES OF THE INTEGRATION
c                     POINTS
c           FV      - FUNCTION VALUES
c           FSTOR
c
c  Compute degree-9 and degree-11 results for IF/DJ and
c  degree-9 approximation for ABS(F)
c
      dj = dabs(u1*v2-u2*v1-u1*v3+v1*u3+u2*v3-v2*u3)*0.5d+00
      f0 = f((u1+u2+u3)/3.0d+00,(v1+v2+v3)/3.0d+00)
      res9 = f0*w90
      resab9 = dabs(f0)*w90
      fv(1) = f0
      kount = 1
      res11 = f0*w110
      if (i.eq.4) res11 = res11 + stor*w(15)
      jend = 15
      if (i.eq.4) jend = 13

      do 50 j=1,jend

        z1 = zeta1(j)
        z2 = zeta2(j)
        z3 = 1.0d+00 - z1 - z2
        x(1) = z1*u1 + z2*u2 + z3*u3
        y(1) = z1*v1 + z2*v2 + z3*v3
        x(2) = z2*u1 + z3*u2 + z1*u3
        y(2) = z2*v1 + z3*v2 + z1*v3
        x(3) = z3*u1 + z1*u2 + z2*u3
        y(3) = z3*v1 + z1*v2 + z2*v3
        if (j.gt.6) go to 20

        f0 = 0.0d+00
        df0 = 0.0d+00
        do l=1,3
          kount = kount + 1
          fv(kount) = f(x(l),y(l))
          f0 = f0 + fv(kount)
          df0 = df0 + dabs(fv(kount))
        end do

        res9 = res9 + f0*w(j)
        resab9 = resab9 + df0*w(j)
        go to 50

   20   f0 = 0.0d+00
        do l=1,3
          fstor(l) = f(x(l),y(l))
          f0 = f0 + fstor(l)
        end do

        if (j.lt.14) go to 40
        stor = stor + fstor(k)
        k = k + 1
        if (k.gt.3) k = 1

   40   continue

        res11 = res11 + f0*w(j)

   50 continue
c
c  Compute degree-9 approximation for the integral of
c  ABS(F-IF/DJ)
c
      dresc = dabs(fv(1)-res9)*w90
      kount = 2
      do j=1,6
        dresc = dresc + (dabs(fv(kount)-res9)+dabs(fv(kount+1)-res9)
     &   +dabs(fv(kount+2)-res9))*w(j)
        kount = kount + 3
      end do
c
c  Compute degree-9 and degree-11 approximations for if,
c  and error estimate
c
      res9 = res9*dj
      res11 = res11*dj
      resab9 = resab9*dj
      dresc = dresc*dj
      est = dabs(res9-res11)
      if (dresc.ne.0.0d+00) est = dmax1(est,dresc*dmin1(1.0d+00,
     & (20.0d+00*est/dresc)**1.5d+00))
      if (resab9.gt.ufrn) est = dmax1(eprn*resab9,est)

      return
      end
      subroutine epsalg ( n, epstab, result, abserr, resla )

c*********************************************************************72
c
cc EPSALG carries out the epsilon algorithm.
c
c  Discussion:
c
c    The routine transforms a given sequence using the
c    epsilon algorithm of P. Wynn.
c    An estimate of the absolute error is also given.
c    the condensed epsilon table is computed. only those
c    elements needed for the computation of the next diagonal
c    are preserved.
c
c  Modified:
c
c    06 December 2006
c
c  Author:
c
c    Elise deDoncker, Ian Robinson
c
c  Reference:
c
c    Elise deDoncker, Ian Robinson,
c    Algorithm 612:
c    Integration over a Triangle Using Nonlinear Extrapolation,
c    ACM Transactions on Mathematical Software,
c    Volume 10, Number 1, March 1984, pages 17-22.
c
c  Parameters:
c
c              N      - INTEGER
c                       EPSTAB(N) CONTAINS THE NEW ELEMENT IN THE
c                       FIRST COLUMN OF THE EPSILON TABLE.
c
c              EPSTAB - DOUBLE PRECISION
c                       ONE DIMENSIONAL ARRAY CONTAINING THE
c                       ELEMENTS OF THE LAST TWO DIAGONALS OF
c                       THE TRIANGULAR EPSILON TABLE
c                       THE ELEMENTS ARE NUMBERED STARTING AT THE
c                       RIGHT-HAND CORNER OF THE TRIANGLE.
c                       THE DIMENSION SHOULD BE AT LEAST (LIMEXP+2)
c                       (SEE DATA LIMEXP).
c
c              RESULT - DOUBLE PRECISION
c                       RESULTING APPROXIMATION TO THE INTEGRAL
c
c              ABSERR - DOUBLE PRECISION
c                       ESTIMATE OF THE ABSOLUTE ERROR
c
c              RESLA  - DOUBLE PRECISION
c                       PREVIOUS RESULT
c
      implicit none

      double precision ABSERR, DELTA1, DELTA2, DELTA3, EMACH, EPRN,
     & EPSINF, EPSTAB, ERROR, ERR1, ERR2, ERR3, E0, E1, E2, E3, E1ABS,
     & OFRN, RES, RESULT, RESLA, SS, TOL1, TOL2, TOL3, UFRN
      double precision DABS, DMAX1
      integer I, IB, IB2, IE, IN, K1, K2, K3, LIMEXP, N, NEWELM, NUM
      DIMENSION EPSTAB(52)

      COMMON /MACH/ EMACH, EPRN, UFRN, OFRN
c
c            LIMEXP  IS THE MAXIMUM NUMBER OF ELEMENTS THE EPSILON
c           TABLE CAN CONTAIN. IF THIS NUMBER IS REACHED, THE UPPER
c           DIAGONAL OF THE EPSILON TABLE IS DELETED.
c           ( EPSTAB  IS OF DIMENSION (LIMEXP+2) AT LEAST.)
c
      data limexp /50/
c
c           LIST OF MAJOR VARIABLES
c           -----------------------
c           E0     - THE 4 ELEMENTS ON WHICH THE
c           E1       COMPUTATION OF A NEW ELEMENT IN
c           E2       THE EPSILON TABLE IS BASED
c           E3                 E0
c                        E3    E1    NEW
c                              E2
c           NEWELM - NUMBER OF ELEMENTS TO BE COMPUTED IN THE NEW
c                    DIAGONAL
c           ERROR  -    ERROR = ABS(E0-E1)+ABS(E1-E2)+ABS(E2-NEW)
c           RESULT - THE ELEMENT IN THE NEW DIAGONAL WITH LEAST
c                       ERROR
c
      epstab(n+2) = epstab(n)
      newelm = (n-1)/2
      epstab(n) = ofrn
      abserr = ofrn
      num = n
      k1 = n

      do 40 i=1,newelm

        k2 = k1 - 1
        k3 = k1 - 2
        res = epstab(k1+2)
        e0 = epstab(k3)
        e1 = epstab(k2)
        e2 = res
        e1abs = dabs(e1)
        delta2 = e2 - e1
        err2 = dabs(delta2)
        tol2 = dmax1(dabs(e2),e1abs)*emach
        delta3 = e1 - e0
        err3 = dabs(delta3)
        tol3 = dmax1(e1abs,dabs(e0))*emach
c
c  IF E0, E1 AND E2 ARE EQUAL TO WITHIN MACHINE
c  ACCURACY, CONVERGENCE IS ASSUMED.
c  RESULT = E2
c  ABSERR = ABS(E1-E0)+ABS(E2-E1)
c
        if (err2.gt.tol2 .or. err3.gt.tol3) go to 10
        result = res
        abserr = err2 + err3
        go to 90

   10   e3 = epstab(k1)
        epstab(k1) = e1
        delta1 = e1 - e3
        err1 = dabs(delta1)
        tol1 = dmax1(e1abs,dabs(e3))*emach
c
c  IF TWO ELEMENTS ARE VERY CLOSE TO EACH OTHER, OMIT
c  A PART OF THE TABLE BY ADJUSTING THE VALUE OF N.
c
        if (err1.le.tol1 .or. err2.le.tol2 .or. err3.le.tol3) go to 20
        ss = 0.1d+01/delta1 + 0.1d+01/delta2 - 0.1d+01/delta3
        epsinf = dabs(ss*e1)
c
c  TEST TO DETECT IRREGULAR BEHAVIOUR IN THE TABLE, AND
c  EVENTUALLY OMIT A PART OF THE TABLE ADJUSTING THE VALUE
c  OF N.
c
        if ( epsinf .gt. 0.1d-03 ) go to 30
   20   n = i + i - 1
        go to 50
c
c  COMPUTE A NEW ELEMENT AND EVENTUALLY ADJUST
c  THE VALUE OF RESULT
c
   30   res = e1 + 0.1d+01 / ss
        epstab(k1) = res
        k1 = k1 - 2
        error = err2 + dabs(res-e2) + err3
        if (error.gt.abserr) go to 40
        abserr = error
        result = res

   40 continue
c
c  SHIFT THE TABLE
c
   50 if (n.eq.limexp) n = 2*(limexp/2) - 1
      ib = 1
      if ((num/2)*2.eq.num) ib = 2
      ie = newelm + 1

      do i=1,ie
        ib2 = ib + 2
        epstab(ib) = epstab(ib2)
        ib = ib2
      end do

      if (num.eq.n) go to 80
      in = num - n + 1
      do i=1,n
        epstab(i) = epstab(in)
        in = in + 1
      end do
c
c  COMPUTE ERROR ESTIMATE
c
   80 abserr = dabs(result-resla)
      resla = result
   90 abserr = dmax1(abserr,eprn*dabs(result))
      return
      end
      subroutine order ( limit, last, maxerr, ermax, elist, iord, nrmax,
     * navail )

c*********************************************************************72
c
cc ORDER carries out a sorting algorithm.
c
c  Discussion:
c
c    This routine maintains the descending ordering
c    in the list of the error estimates resulting from the
c    subdivision process. At each call 4 error estimates
c    are inserted using the sequential search method
c    top-down for the largest estimate, bottom up
c    for the other ones.
c
c  Modified:
c
c    06 December 2006
c
c  Author:
c
c    Elise deDoncker, Ian Robinson
c
c  Reference:
c
c    Elise deDoncker, Ian Robinson,
c    Algorithm 612:
c    Integration over a Triangle Using Nonlinear Extrapolation,
c    ACM Transactions on Mathematical Software,
c    Volume 10, Number 1, March 1984, pages 17-22.
c
c  Parameters:
c
c           LIMIT     - INTEGER
c                       MAXIMUM NUMBER OF ERROR ESTIMATES
c                       THE LIST CAN CONTAIN
c
c           LAST      - INTEGER
c                       NUMBER OF ERROR ESTIMATES CURRENTLY IN
c                       THE LIST  ELIST(LAST) CONTAINS THE SMALLEST
c                       ERROR ESTIMATE
c
c           MAXERR    - INTEGER
c                       POINTS TO THE NRMAX-TH LARGEST ERROR
c                       ESTIMATE CURRENTLY IN THE LIST
c
c           ERMAX     - DOUBLE PRECISION
c                       NRMAX-TH LARGEST ERROR ESTIMATE
c                       ERMAX = ELIST(MAXERR)
c
c           ELIST     - DOUBLE PRECISION
c                       ARRAY OF DIMENSION LIMIT CONTAINING
c                       THE ERROR ESTIMATES
c
c           IORD      - INTEGER
c                       ARRAY CONTAINING POINTERS TO ELIST SO
c                       THAT IORD(1) POINTS TO THE LARGEST
c                       ERROR ESTIMATE,..., IORD(LAST) TO THE
c                       SMALLEST ERROR ESTIMATE
c
c           NRMAX     - INTEGER
c                            MAXERR = IORD(NRMAX)
c
c           NAVAIL    - INTEGER
c                            NAVAIL = LAST+2
c
      implicit none

      double precision elist
      double precision ermax, errmax, errmin
      integer i, ibeg, ii, incr, ind, indi, indj, iord, iord4, isucc,
     & j, jbn, k, last, limit, maxerr, minerr, navail, nrmax, ntest, num
      dimension elist(limit), iord(limit), iord4(4)
c
c  SORT THE 4 NEW ELIST-ENTRIES INTO DESCENDING ORDERING
c
      iord4(1) = maxerr
      iord4(2) = last
      iord4(3) = last + 1
      iord4(4) = last + 2
      do 20 k=1,3
        i = k
        do ii=1,3
          j = i + 1
          indi = iord4(i)
          indj = iord4(j)
          if (elist(indi).ge.elist(indj)) go to 20
          iord4(i) = indj
          iord4(j) = indi
          i = i - 1
          if (i.eq.0) go to 20
        end do
   20 continue
c
c  IF THE LIST CONTAINS ONLY 4 ELEMENTS, RETURN
c
      if (last.gt.2) go to 40
      do i=1,4
        iord(i) = iord4(i)
      end do
      go to 170
c
c  THIS PART OF THE ROUTINE IS ONLY EXECUTED IF, DUE TO A
c  DIFFICULT INTEGRAND, SUBDIVISION INCREASED THE ERROR
c  ESTIMATE. IN THE NORMAL CASE THE INSERT PROCEDURE
c  SHOULD START AFTER THE NRMAX-TH LARGEST ERROR ESTIMATE.
c
   40 maxerr = iord4(1)
      errmax = elist(maxerr)
      if (nrmax.eq.1) go to 60
      ii = nrmax - 1
      do i=1,ii
        isucc = iord(nrmax-1)
        if (errmax.le.elist(isucc)) go to 60
        iord(nrmax) = isucc
        nrmax = nrmax - 1
      end do
c
c  INSERT ERRMAX BY TRAVERSING THE LIST TOP-DOWN
c  STARTING COMPARISON FROM THE ELEMENT ELIST (IORD(NRMAX+1))
c
   60 navail = last + 2
      jbn = navail - 3
      ibeg = nrmax + 1
      if (ibeg.gt.jbn) go to 80

      do i=ibeg,jbn
        isucc = iord(i)
        if (errmax.ge.elist(isucc)) go to 100
        iord(i-1) = isucc
      end do

   80 ind = jbn

      do i=1,4
        iord(ind) = iord4(i)
        ind = ind + 1
      end do

      go to 170
c
c  INSERT THE OTHER ESTIMATES BY TRAVERSING THE LIST BOTTOM-UP
c
  100 continue

      iord(i-1) = maxerr
      k = jbn
      num = jbn - i + 1
      incr = 3

      do ii=1,3

        ntest = num
        minerr = iord4(incr+1)
        errmin = elist(minerr)

        do j=1,ntest
          isucc = iord(k)
          if (errmin.lt.elist(isucc)) go to 150
          ind = k + incr
          iord(ind) = isucc
          k = k - 1
          num = num - 1
        end do

        go to (140, 130, 120), incr
  120   iord(i+2) = iord4(4)
  130   iord(i+1) = iord4(3)
  140   iord(i) = iord4(2)
        go to 170
  150   ind = k + incr
        iord(ind) = minerr
        incr = incr - 1

      end do
c
c  SET MAXERR AND ERMAX
c
  170 maxerr = iord(nrmax)
      ermax = elist(maxerr)

      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
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
