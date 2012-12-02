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
        write ( *, '(a)' ) '  I out of bounds.'
        stop
      end if

      d1mach = dmach(i)

      return
      end
      subroutine dchtri(numfun,mdiv,ver,numtri,minpts,maxpts,epsabs,
     &                  epsrel,lenver,nw,restar,maxsub,minsub,ifail)

c*********************************************************************72
c
cc DCHTRI checks the validity of the input parameters to DCUTRI.
c
c  Discussion:
c
c    DCHTRI computes MAXSUB, MINSUB and IFAIL as
c    functions of the input parameters to DCUTRI,
c    and checks the validity of the input parameters to DCUTRI.
c
c  Modified:
c
c    09 December 2006
c
c  Author:
c
c    Jarle Berntsen, Terje Espelid
c
c  Reference:
c
c    Jarle Berntsen, Terje Espelid,
c    Algorithm 706:
c    DCUTRI, an Algorithm for Adaptive Cubature over a Collection
c    of Triangles,
c    ACM Transactions on Mathematical Software,
c    Volume 18, Number 3, September 1992, pages 329-342.
c
c  Parameters:
c
c    Input, NUMFUN Integer.
c            Number of components of the integral.
c
c    Input, MDIV   Integer.
c            MDIV is the number of triangles that are divided in
c            each subdivision step in DADTRI.
c            MDIV is chosen default to 1.
c            For efficient execution on parallel computers
c            with NPROC processors MDIV should be set equal to
c            the smallest integer such that MOD(4*MDIV,NPROC) = 0.
c
c    Input, VER    Real array of dimension (2,3,LENVER).
c            VER(1,K,L) and VER(2,K,L) are the x and y coordinates
c            respectively of vertex K in triangle L.
c            On entry VER(*,*,L) must contain the vertices of the
c            NUMTRI user specified triangles for L = 1,2,...,NUMTRI.
c
c    Input, NUMTRI Integer.
c            The number of triangles.
c
c    Input, MINPTS Integer.
c            Minimum number of function evaluations.
c
c    Input, MAXPTS Integer.
c            Maximum number of function evaluations.
c            The number of function evaluations over each subregion
c            is 37.
c            MAXPTS >= 37*NUMTRI and MAXPTS >= MINPTS
c
c    Input, EPSABS Real.
c            Requested absolute error.
c
c    Input, EPSREL Real.
c            Requested relative error.
c
c    Input, LENVER Integer.
c            Defines the length of the array VER.
c
c            We let
c            MAXSUB denote the maximum allowed number of subregions
c            for the given values of MAXPTS.
c            MAXSUB = 3*((MAXPTS-37*NUMTRI)/(4*37)) + NUMTRI
c            LENVER should then be greater or equal to MAXSUB.
c
c    Input, NW     Integer.
c            Defines the length of the working array WORK.
c
c            We let
c            MAXSUB denote the maximum allowed number of subregions
c            for the given values of MAXPTS.
c            MAXSUB = 3*((MAXPTS-37*NUMTRI)/(4*37)) + NUMTRI
c            NW should then be greater or equal to
c            MAXSUB*(2*NUMFUN+1) + MAX(32*MDIV,8*NUMTRI)*NUMFUN + 1
c
c    Input, RESTAR Integer.
c            If RESTAR = 0, this is the first attempt to compute
c            the integral.
c            If RESTAR = 1,
c            then we restart a previous attempt.
c            In this case the only parameters for DCUTRI that may
c            be changed (with respect to the previous call of DCUTRI)
c            are MINPTS, MAXPTS, EPSABS, EPSREL and RESTAR.
c
c    Output, MAXSUB Integer.
c            The maximum allowed number of subregions for the
c            given values of MAXPTS.
c
c    Output, MINSUB Integer.
c            The minimum allowed number of subregions for the given
c            values of MINPTS.
c
c    Output, IFAIL  Integer.
c            IFAIL = 0 for normal exit.
c            IFAIL = 2 if NUMFUN is less than 1.
c            IFAIL = 3 if area of one of the initially given
c                      triangles is zero.
c            IFAIL = 4 if MAXPTS is less than 37*NUMTRI.
c            IFAIL = 5 if MAXPTS is less than MINPTS.
c            IFAIL = 6 if EPSABS <= 0 and EPSREL <= 0.
c            IFAIL = 7 if LENVER is less than MAXSUB.
c            IFAIL = 8 if NW is less than MAXSUB*(2*NUMFUN+1) +
c                      NUMFUN*MAX(32*MDIV,8*NUMTRI) + 1.
c            IFAIL = 9 if illegal RESTAR.
c
      implicit none
c
c   Global variables.
c
      integer numfun,mdiv,numtri,minpts,maxpts,nw,maxsub,minsub,restar,
     &        lenver,ifail
      double precision ver(2,3,numtri),epsabs,epsrel
c
c   Local variables.
c
      integer limit,j
      double precision area

      ifail = 0
c
c   We compute MAXSUB and MINSUB.
c
      maxsub = 3* ((maxpts-37*numtri)/ (4*37)) + numtri
      minsub = 3* ((minpts-37*numtri)/ (4*37)) + numtri
      if (mod((minpts-37*numtri),4*37).gt.0) then
          minsub = minsub + 3
      end if
      minsub = max(numtri,minsub)
c
c   Check on positive NUMFUN.
c
      if (numfun.lt.1) then
          ifail = 2
          return
      end if
c
c   Check on legal area of the region of integration.
c
      do j = 1,numtri

          area = abs(ver(1,1,j)*ver(2,2,j)-ver(1,2,j)*ver(2,1,j)-
     &           ver(1,1,j)*ver(2,3,j)+ver(1,3,j)*ver(2,1,j)+
     &           ver(1,2,j)*ver(2,3,j)-ver(1,3,j)*ver(2,2,j))*0.5

          if (area.eq.0) then
              ifail = 3
              return
          end if

      end do
c
c   Check on MAXPTS >= 37*NUMTRI.
c
      if (maxpts.lt.37*numtri) then
          ifail = 4
          return
      end if
c
c   Check on MAXPTS >= MINPTS.
c
      if (maxpts.lt.minpts) then
          ifail = 5
          return
      end if
c
c   Check on legal accuracy requests.
c
      if (epsabs.le.0 .and. epsrel.le.0) then
          ifail = 6
          return
      end if
c
c   Check on big enough LENVER.
c
      if (lenver.lt.maxsub) then
          ifail = 7
          return
      end if
c
c   Check on big enough NW.
c
      limit = maxsub* (2*numfun+1) + max(32*mdiv,8*numtri)*numfun + 1
      if (nw.lt.limit) then
          ifail = 8
          return
      end if
c
c    Check on legal RESTAR.
c
      if (restar.ne.0 .and. restar.ne.1) then
          ifail = 9
          return
      end if

      return
      end
      subroutine dcutri(funsub,numfun,ver,numtri,minpts,maxpts,epsabs,
     &                 epsrel,lenver,nw,restar,result,abserr,neval,
     &                 ifail,work,iwork)

c*********************************************************************72
c
cc DCUTRI integrates a function over a triangulated region.
c
c  Discussion:
c
c            The routine calculates an approximation to a given
c            vector of definite integrals
c
c            I  I (F ,F ,...,F )     DX(2)DX(1),
c                   1  2      NUMFUN
c
c            where the region of integration is a selection of
c            NUMTRI triangles and
c            where F = F (X(1),X(2)), J = 1,2,...,NUMFUN.
c                   J   J
c
c
c            A globally adaptive strategy is applied in order to
c            compute approximations RESULT(K)
c            hopefully satisfying for each component of I the following
c            claim for accuracy:
c            ABS(I(K)-RESULT(K)).LE.MAX(EPSABS,EPSREL*ABS(I(K)))
c
c            DCUTRI is a driver for the integration routine
c            DADTRI.
c
c            DADTRI repeatedly
c            subdivides the triangles with greatest estimated  errors
c            and estimates the integrals and the
c            errors over the new subtriangles
c            until either the error request
c            is met or MAXPTS function evaluations have been used.
c
c            A 37 point integration rule
c            with all evaluation points inside the triangle
c            is applied. The rule has polynomial degree 13.
c
c            If the values of the input parameters EPSABS
c            or EPSREL are selected great enough,
c            an integration rule is applied over each triangle and
c            the results are added up to give the approximations
c            RESULT(K). No further subdivision
c            of the triangles will then be applied.
c
c            When DCUTRI computes estimates to a vector of
c            integrals, all components of the vector are given
c            the same treatment. That is, I(F ) and I(F ) for
c                                            J         K
c            J not equal to K, are estimated with the same
c            subdivision of the region of integration.
c            For integrals with enough similarity, we may save
c            time by applying DCUTRI to all integrands in one call.
c            For integrals that varies continuously as functions of
c            some parameter, the estimates produced by DCUTRI will
c            also vary continuously when the same subdivision is
c            applied to all components. This will generally not be
c            the case when the different components are given
c            separate treatment.
c
c            On the other hand this feature should be used with
c            caution when the different components of the integrals
c            require clearly different subdivisions.
c
c  Modified:
c
c    09 December 2006
c
c  Author:
c
c    Jarle Berntsen, Terje Espelid
c
c  Reference:
c
c    Jarle Berntsen, Terje Espelid,
c    Algorithm 706:
c    DCUTRI, an Algorithm for Adaptive Cubature over a Collection
c    of Triangles,
c    ACM Transactions on Mathematical Software,
c    Volume 18, Number 3, September 1992, pages 329-342.
c
c  Parameters:
c
c   ON ENTRY
c
c     FUNSUB Externally declared subroutine for computing
c            all components of the integrand at the given
c            evaluation point.
c            It must have parameters (X,NUMFUN,FUNVLS)
c            Input parameters:
c              X(1)   The x-coordinate of the evaluation point.
c              X(2)   The y-coordinate of the evaluation point.
c              NUMFUN Integer that defines the number of
c                     components of I.
c            Output parameter:
c              FUNVLS Real array of dimension NUMFUN
c                     that defines NUMFUN components of the integrand.
c
c     NUMFUN Integer.
c            Number of components of the integral.
c     VER    Real array of dimension (2,3,LENVER).
c            VER(1,K,L) and VER(2,K,L) are the x and y coordinates
c            respectively of vertex K in triangle L.
c            On entry VER(*,*,L) must contain the vertices of the
c            NUMTRI user specified triangles for L = 1,2,...,NUMTRI.
c            VER may be changed on exit.
c
c     NUMTRI Integer.
c            The number of triangles.
c     MINPTS Integer.
c            Minimum number of function evaluations.
c     MAXPTS Integer.
c            Maximum number of function evaluations.
c            The number of function evaluations over each subregion
c            is 37.
c
c            MAXPTS >= 37*NUMTRI and MAXPTS >= MINPTS
c
c     EPSABS Real.
c            Requested absolute error.
c     EPSREL Real.
c            Requested relative error.
c     LENVER Integer.
c            Defines the length of the array VER.
c
c            We let
c            MAXSUB denote the maximum allowed number of subregions
c            for the given value of MAXPTS.
c            MAXSUB = 3*((MAXPTS-37*NUMTRI)/(4*37)) + NUMTRI
c            LENVER should be greater or equal to MAXSUB.
c
c     NW     Integer.
c            Defines the length of the working array WORK.
c
c            If LENVER is chosen correctly as a function of MAXPTS
c            and NUMTRI and NW is selected to be greater or equal to
c            LENVER*(2*NUMFUN+1) + MAX(32*MDIV,8*NUMTRI)*NUMFUN + 1 ,
c            then the size of the working storage will be great enough.
c            MDIV is the number of triangles that are divided in
c            each subdivision step.
c            MDIV is default set internally in DCUTRI equal to 1.
c            For efficient execution on parallel computers
c            with NPROC processors MDIV should be set equal to
c            the smallest integer such that MOD(4*MDIV,NPROC) = 0.
c
c     RESTAR Integer.
c            If RESTAR = 0, this is the first attempt to compute
c            the integral.
c            If RESTAR = 1,
c            then we restart a previous attempt.
c            In this case the only parameters for DCUTRI that may
c            be changed (with respect to the previous call of DCUTRI)
c            are MINPTS, MAXPTS, EPSABS, EPSREL and RESTAR.
c            Note: If MAXPTS was too small in the previous call,
c            then we get a second chance to continue the computations
c            with MAXPTS increased.
c            LENVER can not be changed, therefore the new value of MAXPTS
c            should not be chosen greater than the value of LENVER
c            allows us to do.
c     WORK   Real array of dimension NW.
c            Used as working storage.
c     IWORK  Integer array of dimension LENVER + MDIV.
c            Used as working storage.
c
c   ON RETURN
c
c     RESULT Real array of dimension NUMFUN.
c            Approximations to all components of the integral.
c     ABSERR Real array of dimension NUMFUN.
c            Estimates of absolute errors.
c     NEVAL  Integer.
c            Number of function evaluations used by DCUTRI.
c     IFAIL  Integer.
c            IFAIL = 0 for normal exit.
c
c              ABSERR(K) <=  EPSABS or
c              ABSERR(K) <=  ABS(RESULT(K))*EPSREL with MAXPTS or less
c              function evaluations for all values of K,
c              1 <= K <= NUMFUN .
c
c            IFAIL = 1 if MAXPTS was too small for DCUTRI
c              to obtain the required accuracy. In this case DCUTRI
c              returns values of RESULT with estimated absolute
c              errors ABSERR.
c            IFAIL = 2 if NUMFUN is less than 1.
c            IFAIL = 3 if area of one of the initially given
c                      triangles is zero.
c            IFAIL = 4 if MAXPTS is less than 37*NUMTRI.
c            IFAIL = 5 if MAXPTS is less than MINPTS.
c            IFAIL = 6 if EPSABS <= 0 and EPSREL <= 0.
c            IFAIL = 7 if LENVER is less than MAXSUB.
c            IFAIL = 8 if NW is less than MAXSUB*(2*NUMFUN+1) +
c                      NUMFUN*MAX(32*MDIV,8*NUMTRI) + 1.
c            IFAIL = 9 if unlegal RESTAR.
c     VER    Real array of dimension (2,3,LENVER).
c            On exit VER
c            contains the vertices of the triangles used to produce
c            the approximations to the integrals.
c     WORK   Real array of dimension NW.
c            Used as working storage.
c            WORK(NW) = NSUB, the number of subregions in the data
c            structure.
c            Let WRKSUB=(NW-1-NUMFUN*MAX(32*MDIV,8*NUMTRI))/
c                        (2*NUMFUN+1).
c            WORK(1),...,WORK(NUMFUN*WRKSUB) contain
c              the estimated components of the integrals over the
c              subregions.
c            WORK(NUMFUN*WRKSUB+1),...,WORK(2*NUMFUN*WRKSUB) contain
c              the estimated errors over the subregions.
c            WORK(2*NUMFUN*WRKSUB+1),...,WORK(2*NUMFUN*
c              WRKSUB+WRKSUB) contain the greatest errors
c              in each subregion.
c            WORK((2*NUMFUN+1)*WRKSUB),...,WORK(NW - 1) is used as
c              temporary storage in DADTRI.
c     IWORK  Integer array of dimension LENVER + MDIV.
c            Used as working storage.
c
c           SAMPLE PROGRAM
c   The following program will approximate the integral of
c                     exp(x*x+y*y)
c   over the triangle (0.,0.),(2.,0.),(0.,2.) and print out the
c   values of the estimated integral, the estimated error, the number
c   of function evaluations, and IFAIL.
c     PROGRAM DTEST1
c     INTEGER NUMFUN,NW,MDIV,MINPTS,LENVER,NUMTRI
c     PARAMETER (NUMFUN=1,MDIV=1,NUMTRI=1,MINPTS=37)
c     PARAMETER ( LENVER = 100 )
c
c   We need to choose LENVER  such that:
c   LENVER >= 3*((MAXPTS-37*NUMTRI)/(4*37)) + NUMTRI
c   This simply means that we have enough space to achieve MAXPTS
c   function evaluations. By choosing LENVER bigger than
c   this lower bound we may later change MAXPTS and RESTAR and restart
c   the computations from the point where the last run stopped.
c   The reason for stopping may have been that MAXPTS was too small
c   to achieve the requested error.
c          With our choice LENVER = 100 any value of MAXPTS
c   between 37 and 5068 will be accepted by the code. Choosing
c   MAXPTS = 2000 allows us to restart with a greater value
c   if necessary.
c
c     PARAMETER (NW = LENVER*(2*NUMFUN+1) +
c    +                MAX(32*MDIV,8*NUMTRI)*NUMFUN + 1)
c     DOUBLE PRECISION VER(2,3,LENVER),RESULT(NUMFUN),
c    +                 ABSERR(NUMFUN),WORK(NW),EPSABS,EPSREL
c     INTEGER RESTAR,NEVAL,IFAIL,MDIV,IWORK(LENVER+MDIV),MAXPTS
c     EXTERNAL F
c     VER(1,1,1) = 0.D0
c     VER(1,2,1) = 2.D0
c     VER(1,3,1) = 0.D0
c     VER(2,1,1) = 0.D0
c     VER(2,2,1) = 0.D0
c     VER(2,3,1) = 2.D0
c     EPSABS = 0.D0
c     EPSREL = 1.D-5
c     RESTAR = 0
c     MAXPTS = 2000
c     CALL DCUTRI(F,NUMFUN,VER,NUMTRI,MINPTS,MAXPTS,EPSABS,
c    +            EPSREL,LENVER,NW,RESTAR,RESULT,ABSERR,NEVAL,
c    +            IFAIL,WORK,IWORK)
c     WRITE(*,*)'RESULT=',RESULT(1)
c     WRITE(*,*)'ERROR =',ABSERR(1)
c     WRITE(*,*)'NEVAL =',NEVAL
c     WRITE(*,*)'IFAIL =',IFAIL
c     END
c     SUBROUTINE F(X,NUMFUN,FUNVLS)
c     INTEGER NUMFUN
c     DOUBLE PRECISION X(2),FUNVLS(NUMFUN)
c     FUNVLS(1) = EXP(X(1)*X(1)+X(2)*X(2))
c     RETURN
c     END
c
c   Output produced on a SUN SPARC station 1.
c
c       RESULT=   11.181284417019
c       ERROR =    4.7009048669038D-06
c       NEVAL =  185
c       IFAIL =  0
c
c
c***LONG DESCRIPTION
c
c   The information for each triangle is contained in the
c   data structures VER, WORK and IWORK.
c   VER contains the coordinates of the triangles.
c   When passed on to DADTRI, WORK is split into four
c   arrays VALUES, ERRORS, GREATE and WORK2.
c   VALUES contains the estimated values of the integrals.
c   ERRORS contains the estimated errors.
c   GREATE contains the greatest estimated error for each triangle.
c   The data structures for the triangles are in DADTRI organized
c   as a heap, and the size of GREATE(I) defines the position of
c   triangle I in the heap. The heap is maintained by the program
c   DTRTRI and we use a partially ordered list of pointers, saved
c   in IWORK. When passed on to DADTRI, IWORK is split into two
c   arrays LIST and VACANT. LIST is a partially ordered list
c   such that GREATE(LIST(1)) is the maximum error estimate for
c   all sub-triangles in our datastructure. VACANT is a work space
c   needed in the updating process of the list.
c
c   The subroutine DADTRI is written for efficient execution on shared
c   memory parallel computer. On a computer with NPROC processors we will
c   in each subdivision step divide MDIV triangles, where MDIV is
c   chosen such that MOD(4*MDIV,NPROC) = 0, in totally 4*MDIV new triangles.
c   Each processor will then compute estimates of the integrals and errors
c   over 4*MDIV/NPROC triangles in each subdivision step.
c   The subroutine for estimating the integral and the error over
c   each subregion, DRLTRI, uses WORK2 as a work array.
c   We must make sure that each processor writes its results to
c   separate parts of the memory, and therefore the sizes of WORK and
c   WORK2 are functions of MDIV.
c   In order to achieve parallel processing of triangles, compiler
c   directives should be placed in front of the DO 20 and the DO 200
c   loops in DADTRI on machines like Alliant and CRAY.
c
c  Local Parameters:
c
c   MDIV   Integer.
c          MDIV is the number of triangles that are divided in
c          each subdivision step in DADTRI.
c          MDIV is chosen default to 1.
c          For efficient execution on parallel computers
c          with NPROC processors MDIV should be set equal to
c          the smallest integer such that MOD(4*MDIV,NPROC) = 0.
c
c   MAXSUB Integer.
c          The maximum allowed number of subregions
c          for the given values of MAXPTS.
c
c   MINSUB Integer.
c          The minimum allowed number of subregions for the given
c          values of MINPTS.
c
c   WRKSUB Integer.
c          Determines the length of the main work arrays.
c
      implicit none

      integer mdiv
      parameter (mdiv=1)

      external funsub
      integer numfun,numtri,minpts,maxpts,lenver,nw,restar,neval,ifail,
     &        iwork(lenver+mdiv)
      double precision ver(2,3,lenver),epsabs,epsrel,result(numfun),
     &                 abserr(numfun),work(nw)

      integer maxsub,minsub,nsub,lenw
      integer wrksub,i1,i2,i3,i4
c
c   Compute MAXSUB and MINSUB,
c   and check the input parameters.
c
      call dchtri(numfun,mdiv,ver,numtri,minpts,maxpts,epsabs,epsrel,
     &            lenver,nw,restar,maxsub,minsub,ifail)
      wrksub = (nw-1-numfun*max(32*mdiv,8*numtri))/ (2*numfun+1)
      if (ifail.ne.0) then
          return
      end if
c
c   Split up the work space.
c
      i1 = 1
      i2 = i1 + wrksub*numfun
      i3 = i2 + wrksub*numfun
      i4 = i3 + wrksub
c
c   On restart runs the number of subregions from the
c   previous call is assigned to NSUB.
c
      if (restar.eq.1) then
          nsub = work(nw)
      end if
c
c   Compute the size of the temporary work space needed in DADTRI.
c
      lenw = max(32*mdiv,8*numtri)*numfun
      call dadtri(numfun,mdiv,ver,numtri,minsub,maxsub,funsub,epsabs,
     &            epsrel,lenver,restar,lenw,result,abserr,neval,nsub,
     &            ifail,work(i1),work(i2),work(i3),work(i4),
     &            iwork(1),iwork(1+lenver))
      work(nw) = nsub

      return
      end
      subroutine drltri(ver,numfun,funsub,null,basval,rgnerr,greate)

c*********************************************************************72
c
cc DRLTRI computes basic integration rules and error estimates.
c
c  Discussion:
c
c    DRLTRI computes basic integration rule values
c    for a vector of integrands over a triangular region.
c    DRLTRI also computes estimates for the errors by
c    using several null rule approximations.
c
c  Modified:
c
c    09 December 2006
c
c  Author:
c
c    Jarle Berntsen, Terje Espelid
c
c  Reference:
c
c    Jarle Berntsen, Terje Espelid,
c    Algorithm 706:
c    DCUTRI, an Algorithm for Adaptive Cubature over a Collection
c    of Triangles,
c    ACM Transactions on Mathematical Software,
c    Volume 18, Number 3, September 1992, pages 329-342.
c
c    Jarle Berntsen, Terje Espelid,
c    Degree 13 Symmetric Quadrature Rules for the Triangle, 
c    Reports in Informatics 44, 
c    Department of Informatics, University of Bergen, 1990.
c
c  Parameters:
c
c   ON ENTRY
c
c   VER    Real array of dimension (2,3).
c          The coordinates of the vertices of the triangle.
c   NUMFUN Integer.
c          Number of components of the vector integrand.
c   FUNSUB Externally declared subroutine for computing
c            all components of the integrand at the given
c            evaluation point.
c            It must have parameters (X,NUMFUN,FUNVLS)
c            Input parameters:
c              X(1)      The x-coordinate of the evaluation point.
c              X(2)      The y-coordinate of the evaluation point.
c              NUMFUN Integer that defines the number of
c                     components of I.
c            Output parameter:
c              FUNVLS Real array of dimension NUMFUN
c                     that defines NUMFUN components of the integrand.
c
c   NULL   Real array of dimension (NUMFUN,8).
c          A work array.
c
c   ON RETURN
c
c   BASVAL Real array of dimension NUMFUN.
c          The values for the basic rule for each component
c          of the integrand.
c   RGNERR Real array of dimension NUMFUN.
c          The error estimates for each component of the integrand.
c   GREATE Real.
c          The greatest error component of the integrand.
c
c  Local parameters:
c
c   WTLENG The number of weights of the integration rule.
c
c   G      Real array of dimension (2,WTLENG).
c          The homogeneous coordinates for the generators of
c          the evaluation points.
c          The integration rule is using symmetric evaluation
c          points and has the structure (1,6,3). That is,
c          1 point of multiplicity 1,
c          6 sets of points of multiplicity 3 and
c          3 sets of points of multiplicity 6.
c          This gives totally 37 evaluation points.
c          In order to reduce the number of loops in DRLTRI,
c          the 3 loops for the sets of multiplicity 6 are split
c          into 6 loops and added to the loops for the sets of
c          multiplicity 3.
c          The number of weights we have to give with
c          this splitting is WTLENG = 13.
c
c   W      Real array of dimension (9,WTLENG).
c          The weights of the basic rule and the null rules.
c          W(1,1),...,W(1,WTLENG) are weights for the basic rule.
c          W(I,1),...,W(I,WTLENG) for I>1 are null rule weights.
c
      implicit none

      external funsub
      integer numfun
      double precision ver(2,3),basval(numfun),rgnerr(numfun),
     &                 null(numfun,8),greate,noise,d1mach,tres

      integer wtleng
      parameter (wtleng=13)
      double precision area,x(2,3),z1,z2,z3,r1,r2,r3,r,deg7,deg5,deg3,
     &                 deg1,g(2,13),w(9,13)
      integer i,j,k,l
c
c  The abscissas are given in homogeneous coordinates.
c
      data (g(1,i),i=1,13)/0.333333333333333333333333333333d0,
     &     0.950275662924105565450352089520d0,
     &     0.171614914923835347556304795551d0,
     &     0.539412243677190440263092985511d0,
     &     0.772160036676532561750285570113d0,
     &     0.009085399949835353883572964740d0,
     &     0.062277290305886993497083640527d0,
     &     0.022076289653624405142446876931d0,
     &     0.018620522802520968955913511549d0,
     &     0.096506481292159228736516560903d0,
     &     0.851306504174348550389457672223d0,
     &     0.689441970728591295496647976487d0,
     &     0.635867859433872768286976979827d0/
      data (g(2,i),i=1,13)/0.333333333333333333333333333333d0,
     &     0.024862168537947217274823955239d0,
     &     0.414192542538082326221847602214d0,
     &     0.230293878161404779868453507244d0,
     &     0.113919981661733719124857214943d0,
     &     0.495457300025082323058213517632d0,
     &     0.468861354847056503251458179727d0,
     &     0.851306504174348550389457672223d0,
     &     0.689441970728591295496647976487d0,
     &     0.635867859433872768286976979827d0,
     &     0.022076289653624405142446876931d0,
     &     0.018620522802520968955913511549d0,
     &     0.096506481292159228736516560903d0/
c
c   Weights of the degree 13 quadrature rule.
c
      data (w(1,i),i=1,13)/0.051739766065744133555179145422d0,
     &     0.008007799555564801597804123460d0,
     &     0.046868898981821644823226732071d0,
     &     0.046590940183976487960361770070d0,
     &     0.031016943313796381407646220131d0,
     &     0.010791612736631273623178240136d0,
     &     0.032195534242431618819414482205d0,
     &     0.015445834210701583817692900053d0,
     &     0.017822989923178661888748319485d0,
     &     0.037038683681384627918546472190d0,
     &     0.015445834210701583817692900053d0,
     &     0.017822989923178661888748319485d0,
     &     0.037038683681384627918546472190d0/
c
c   Weights of the first null rule of degree 7.
c
      data (w(2,i),i=1,13)/-0.077738051051462052051304462750d0,
     &     0.001640389740236881582083124927d0,
     &     0.078124083459915167386776552733d0,
     &     -0.030706528522391137165581298102d0,
     &     0.010246307817678312345028512621d0,
     &     0.012586300774453821540476193059d0,
     &     -0.043630506151410607808929481439d0,
     &     -0.004567055157220063810223671248d0,
     &     0.003393373439889186878847613140d0,
     &     0.000000000000000000000000000000d0,
     &     -0.004567055157220063810223671248d0,
     &     0.003393373439889186878847613140d0,
     &     0.000000000000000000000000000000d0/
c
c   Weights of the second null rule of degree 7.
c
      data (w(3,i),i=1,13)/-0.064293709240668260928898888457d0,
     &     0.003134665264639380635175608661d0,
     &     0.007822550509742830478456728602d0,
     &     0.048653051907689492781049400973d0,
     &     0.032883327334384971735434067029d0,
     &     -0.017019508374229390108580829589d0,
     &     0.025973557893399824586684707198d0,
     &     -0.010716753326806275930657622320d0,
     &     0.018315629578968063765722278290d0,
     &     -0.047607080313197299401024682666d0,
     &     -0.010716753326806275930657622320d0,
     &     0.018315629578968063765722278290d0,
     &     -0.047607080313197299401024682666d0/
c
c   Weights of the first degree 5 null rule.
c
c
      data (w(4,i),i=1,13)/0.021363205584741860993131879186d0,
     &     0.022716410154120323440432428315d0,
     &     -0.026366191271182090678117381002d0,
     &     0.029627021479068212693155637482d0,
     &     0.004782834546596399307634111034d0,
     &     0.004178667433984132052378990240d0,
     &     -0.065398996748953861618846710897d0,
     &     -0.033589813176131630980793760168d0,
     &     0.033018320112481615757912576257d0,
     &     0.012241086002709814125707333127d0,
     &     -0.033589813176131630980793760168d0,
     &     0.033018320112481615757912576257d0,
     &     0.012241086002709814125707333127d0/
c
c   Weights of the second degree 5 null rule.
c
      data (w(5,i),i=1,13)/-0.046058756832790538620830792345d0,
     &     0.005284159186732627192774759959d0,
     &     0.009325799301158899112648198129d0,
     &     -0.006101110360950124560783393745d0,
     &     -0.056223328794664871336486737231d0,
     &     -0.062516479198185693171971930698d0,
     &     0.022428226812039547178810743269d0,
     &     -0.000026014926110604563130107142d0,
     &     0.032882099937471182365626663487d0,
     &     0.018721740987705986426812755881d0,
     &     -0.000026014926110604563130107142d0,
     &     0.032882099937471182365626663487d0,
     &     0.018721740987705986426812755881d0/
c
c   Weights of first degree 3 null rule.
c
      data (w(6,i),i=1,13)/0.080867117677405246540283712799d0,
     &     -0.033915806661511608094988607349d0,
     &     0.014813362053697845461526433401d0,
     &     0.001442315416337389214102507204d0,
     &     -0.024309696484708683486455879210d0,
     &     -0.005135085639122398522835391664d0,
     &     -0.034649417896235909885490654650d0,
     &     0.035748423431577326597742956780d0,
     &     0.024548155266816447583155562333d0,
     &     -0.032897267038856299280541675107d0,
     &     0.035748423431577326597742956780d0,
     &     0.024548155266816447583155562333d0,
     &     -0.032897267038856299280541675107d0/
c
c   Weights of second degree 3 null rule.
c
      data (w(7,i),i=1,13)/-0.038457863913548248582247346193d0,
     &     -0.055143631258696406147982448269d0,
     &     -0.021536994314510083845999131455d0,
     &     0.001547467894857008228010564582d0,
     &     0.057409361764652373776043522086d0,
     &     -0.040636938884669694118908764512d0,
     &     -0.020801144746964801777584428369d0,
     &     0.019490770404993674256256421103d0,
     &     0.002606109985826399625043764771d0,
     &     0.023893703367437102825618048130d0,
     &     0.019490770404993674256256421103d0,
     &     0.002606109985826399625043764771d0,
     &     0.023893703367437102825618048130d0/
c
c   Weights of first degree 1 null rule.
c
      data (w(8,i),i=1,13)/0.074839568911184074117081012527d0,
     &     -0.004270103034833742737299816615d0,
     &     0.049352639555084484177095781183d0,
     &     0.048832124609719176627453278550d0,
     &     0.001042698696559292759051590242d0,
     &     -0.044445273029113458906055765365d0,
     &     -0.004670751812662861209726508477d0,
     &     -0.015613390485814379318605247424d0,
     &     -0.030581651696100000521074498679d0,
     &     0.010801113204340588798240297593d0,
     &     -0.015613390485814379318605247424d0,
     &     -0.030581651696100000521074498679d0,
     &     0.010801113204340588798240297593d0/
c
c   Weights of second degree 1 null rule.
c
      data (w(9,i),i=1,13)/0.009373028261842556370231264134d0,
     &     -0.074249368848508554545399978725d0,
     &     0.014709707700258308001897299938d0,
     &     0.009538502545163567494354463302d0,
     &     -0.014268362488069444905870465047d0,
     &     0.040126396495352694403045023109d0,
     &     0.028737181842214741174950928350d0,
     &     -0.031618075834734607275229608099d0,
     &     0.016879961075872039084307382161d0,
     &     0.010878914758683152984395046434d0,
     &     -0.031618075834734607275229608099d0,
     &     0.016879961075872039084307382161d0,
     &     0.010878914758683152984395046434d0/

      tres = 50*d1mach(4)
c
c  Compute area of triangle.
c
      area = abs(ver(1,1)*ver(2,2)-ver(1,2)*ver(2,1)-ver(1,1)*ver(2,3)+
     &       ver(1,3)*ver(2,1)+ver(1,2)*ver(2,3)-ver(1,3)*ver(2,2))/2
c
c   Compute contributions from the center of the triangle.
c
      x(1,1) = (ver(1,1)+ver(1,2)+ver(1,3))/3
      x(2,1) = (ver(2,1)+ver(2,2)+ver(2,3))/3

      call funsub(x,numfun,rgnerr)

      do j = 1,numfun
          basval(j) = w(1,1)*rgnerr(j)
          do k = 1,8
              null(j,k) = w(k+1,1)*rgnerr(j)
          end do
      end do
c
c   Compute the contributions from the points with
c   multiplicity 3.
c
      do i = 2,wtleng
          z1 = g(1,i)
          z2 = g(2,i)
          z3 = 1 - z1 - z2
          x(1,1) = z1*ver(1,1) + z2*ver(1,2) + z3*ver(1,3)
          x(2,1) = z1*ver(2,1) + z2*ver(2,2) + z3*ver(2,3)
          x(1,2) = z2*ver(1,1) + z3*ver(1,2) + z1*ver(1,3)
          x(2,2) = z2*ver(2,1) + z3*ver(2,2) + z1*ver(2,3)
          x(1,3) = z3*ver(1,1) + z1*ver(1,2) + z2*ver(1,3)
          x(2,3) = z3*ver(2,1) + z1*ver(2,2) + z2*ver(2,3)
          do l = 1,3
              call funsub(x(1,l),numfun,rgnerr)
              do j = 1,numfun
                  basval(j) = basval(j) + w(1,i)*rgnerr(j)
                  do k = 1,8
                      null(j,k) = null(j,k) + w(k+1,i)*rgnerr(j)
                  end do
              end do
          end do
      end do
c
c    Compute the errors.
c
      greate = 0

      do j = 1,numfun

          deg7 = area*sqrt(null(j,1)**2+null(j,2)**2)
          deg5 = area*sqrt(null(j,3)**2+null(j,4)**2)
          deg3 = area*sqrt(null(j,5)**2+null(j,6)**2)
          deg1 = area*sqrt(null(j,7)**2+null(j,8)**2)

          if (deg5.ne.0) then
              r1 = deg7/deg5
          else
              r1 = 1
          end if

          if (deg3.ne.0) then
              r2 = deg5/deg3
          else
              r2 = 1
          end if

          if (deg1.ne.0) then
              r3 = deg3/deg1
          else
              r3 = 1
          end if

          r = max(r1,r2,r3)

          if (r.ge.1) then
              rgnerr(j) = 10*max(deg1,deg3,deg5,deg7)
          else if (r.ge.0.5d0) then
              rgnerr(j) = 10*r*deg7
          else
              rgnerr(j) = 40* (r**3)*deg7
          end if

          basval(j) = area*basval(j)
          noise = abs(basval(j))*tres
c
c  The following 6 statements added as the result of authors remark
c
c    The following statement is included to set the error to the noise
c    level if the two error estimates, assumed to be the best, are both
c    below or on the noise level
c
          if ((deg7.le.noise).and.(deg5.le.noise)) rgnerr(j)=noise
          rgnerr(j) = max(noise,rgnerr(j))

          if (rgnerr(j).gt.greate) then
              greate = rgnerr(j)
          end if

      end do

      return
      end
      subroutine dtrtri ( dvflag, sbrgns, greate, list, new )

c*********************************************************************72
c
cc DTRTRI maintains a heap of subregions.
c
c  Discussion:
c
c    DTRTRI maintains a heap of subregions.  The subregions are stored
c    in a partially sorted binary tree, ordered according to the size
c    of the greatest error estimates of each subregion(GREATE).
c    The subregion with greatest error estimate is in the
c    first position of the heap.
c
c  Modified:
c
c    09 December 2006
c
c  Author:
c
c    Jarle Berntsen, Terje Espelid
c
c  Reference:
c
c    Jarle Berntsen, Terje Espelid,
c    Algorithm 706:
c    DCUTRI, an Algorithm for Adaptive Cubature over a Collection
c    of Triangles,
c    ACM Transactions on Mathematical Software,
c    Volume 18, Number 3, September 1992, pages 329-342.
c
c  Parameters:
c
c     DVFLAG Integer.
c            If DVFLAG = 1, we remove the subregion with
c            greatest error from the heap.
c            If DVFLAG = 2, we insert a new subregion in the heap.
c
c     SBRGNS Integer.
c            Number of subregions in the heap.
c
c     GREATE Real array of dimension SBRGNS.
c            Used to store the greatest estimated errors in
c            all subregions.
c
c     LIST   Integer array of dimension SBRGNS.
c            Used as a partially ordered list of pointers to the
c            different subregions. This list is a heap where the
c            element on top of the list is the subregion with the
c            greatest error estimate.
c
c     NEW    Integer.
c            Index to the new region to be inserted in the heap.
c
c  Local parameters:
c
c   GREAT  is used as intermediate storage for the greatest error of a
c          subregion.
c
c   SUBRGN Position of child/parent subregion in the heap.
c
c   SUBTMP Position of parent/child subregion in the heap.
c
      implicit none

      integer dvflag,new,sbrgns,list(*)
      double precision greate(*)

      integer subrgn,subtmp
      double precision great
c
c    If DVFLAG = 1, we will reduce the partial ordered list by the
c    element with greatest estimated error. Thus the element in
c    in the heap with index LIST(1) is vacant and can be used later.
c    Reducing the heap by one element implies that the last element
c    should be re-positioned.
c
      if (dvflag.eq.1) then
          great = greate(list(sbrgns))
          sbrgns = sbrgns - 1
          subrgn = 1
   20     subtmp = 2*subrgn
          if (subtmp.le.sbrgns) then
              if (subtmp.ne.sbrgns) then
c
c   Find max. of left and right child.
c
                  if (greate(list(subtmp)).lt.
     &                greate(list(subtmp+1))) then
                      subtmp = subtmp + 1
                  end if
              end if
c
c   Compare max.child with parent.
c   If parent is max., then done.
c
              if (great.lt.greate(list(subtmp))) then
c
c   Move the pointer at position SUBTMP up the heap.
c
                  list(subrgn) = list(subtmp)
                  subrgn = subtmp
                  go to 20
              end if
          end if
c
c   Update the pointer.
c
          if (sbrgns.gt.0) then
              list(subrgn) = list(sbrgns+1)
          end if
      else if (dvflag.eq.2) then
c
c   If DVFLAG = 2, find the position for the NEW region in the heap.
c
          great = greate(new)
          subrgn = sbrgns
   40     subtmp = subrgn/2
          if (subtmp.ge.1) then
c
c   Compare max.child with parent.
c   If parent is max, then done.
c
              if (great.gt.greate(list(subtmp))) then
c
c   Move the pointer at position SUBTMP down the heap.
c
                  list(subrgn) = list(subtmp)
                  subrgn = subtmp
                  go to 40
              end if
          end if
c
c    Set the pointer to the new region in the heap.
c
          list(subrgn) = new
      end if
      return
      end
      subroutine dadtri(numfun,mdiv,ver,numtri,minsub,maxsub,funsub,
     &                  epsabs,epsrel,lenver,restar,lenw,result,abserr,
     &                  neval,nsub,ifail,values,errors,greate,work2,
     &                  list,vacant)

c*********************************************************************72
c
cc DADTRI computes integrals over triangular regions.
c
c  Discussion:
c
c    DADTRI repeatedly subdivides the triangles with greatest estimated
c    errors and estimates the integrals and the errors over the new
c    subtriangles until either the error request is met or MAXPTS
c    function evaluations have been used.
c
c    A 37 point integration rule with all evaluation points inside
c    the triangle is applied.  The rule has polynomial degree 13.
c
c  Modified:
c
c    09 December 2006
c
c  Author:
c
c    Jarle Berntsen, Terje Espelid
c
c  Reference:
c
c    Jarle Berntsen, Terje Espelid,
c    Algorithm 706:
c    DCUTRI, an Algorithm for Adaptive Cubature over a Collection
c    of Triangles,
c    ACM Transactions on Mathematical Software,
c    Volume 18, Number 3, September 1992, pages 329-342.
c
c  Parameters:
c
c   ON ENTRY
c
c     NUMFUN Integer.
c            Number of components of the integral.
c
c     MDIV   Integer.
c            MDIV is the number of triangles that are divided in
c            each subdivision step in DADTRI.
c            MDIV is chosen by default to be 1.
c            For efficient execution on parallel computers
c            with NPROC processors MDIV should be set equal to
c            the smallest integer such that MOD(4*MDIV,NPROC) = 0.
c
c     VER    Real array of dimension (2,3,LENVER).
c            VER(1,K,L) and VER(2,K,L) are the x and y coordinates
c            respectively of vertex K in triangle L.
c            On entry VER(*,*,L) must contain the vertices of the
c            NUMTRI user specified triangles for L = 1,2,...,NUMTRI.
c
c     NUMTRI Integer.
c            The number of triangles.
c
c     MINSUB Integer.
c            The minimum allowed number of subregions.
c
c     MAXSUB Integer.
c            The maximum allowed number of subregions.
c
c     FUNSUB Externally declared subroutine for computing
c            all components of the integrand at the given
c            evaluation point.
c            It must have parameters (X,NUMFUN,FUNVLS)
c            Input parameters:
c              X(1)      The x-coordinate of the evaluation point.
c              X(2)      The y-coordinate of the evaluation point.
c              NUMFUN Integer that defines the number of
c                     components of I.
c            Output parameter:
c              FUNVLS Real array of dimension NUMFUN
c                     that defines NUMFUN components of the integrand.
c
c     EPSABS Real.
c            Requested absolute error.
c
c     EPSREL Real.
c            Requested relative error.
c
c     LENVER Integer.
c            Defines the length of the array VER.
c
c            We let
c            MAXSUB denote the maximum allowed number of subregions
c            for the given values of MAXPTS.
c            MAXSUB = 3*((MAXPTS-37*NUMTRI)/(4*37)) + NUMTRI
c            LENVER should then be greater or equal to MAXSUB.
c
c     RESTAR Integer.
c            If RESTAR = 0, this is the first attempt to compute
c            the integral.
c            If RESTAR = 1,
c            then we restart a previous attempt.
c            In this case the only parameters for DCUTRI that may
c            be changed (with respect to the previous call of DCUTRI)
c            are MINPTS, MAXPTS, EPSABS, EPSREL and RESTAR.
c
c     LENW   Integer.
c            Length of the workspace WORK2.
c
c   ON RETURN
c
c     RESULT Real array of dimension NUMFUN.
c            Approximations to all components of the integral.
c
c     ABSERR Real array of dimension NUMFUN.
c            Estimates of absolute errors.
c
c     NEVAL  Integer.
c            Number of function evaluations used by DCUTRI.
c
c     NSUB   Integer.
c            The number of triangles in the data structure.
c
c     IFAIL  Integer.
c            IFAIL = 0 for normal exit.
c
c              ABSERR(K) <=  EPSABS or
c              ABSERR(K) <=  ABS(RESULT(K))*EPSREL with MAXPTS or less
c              function evaluations for all values of K,
c              1 <= K <= NUMFUN .
c
c            IFAIL = 1 if MAXSUB was too small for DADTRI
c              to obtain the required accuracy. In this case DADTRI
c              returns values of RESULT with estimated absolute
c              errors ABSERR.
c
c     VALUES Real array of dimension (NUMFUN,MAXSUB).
c            The estimated components of the integrals over the
c            subregions.
c
c     ERRORS Real array of dimension (NUMFUN,MAXSUB).
c            The estimated errors over the subregions.
c
c     GREATE Real array of dimension MAXSUB.
c            The greatest errors in each subregion.
c
c     WORK2   Real array of dimension LENW.
c            Work array used in DTRTRI and DRLTRI.
c
c     LIST   Integer array used in DTRTRI of dimension LENVER.
c            Is a partially sorted list, where LIST(1) is the top
c            element in a heap of subregions.
c
c     VACANT Integer array of dimension MDIV.
c            Used as intermediate storage of pointers.
c
c  Local parameters:
c
c   SBRGNS is the number of stored subregions.
c
c   NDIV   The number of subregions to be divided in each main step.
c
c   POINTR Pointer to the position in the data structure where
c          the new subregions are to be stored.
c
c   TOP    is a pointer to the top element in the heap of subregions.
c
      implicit none

      external funsub
      integer numfun,mdiv,numtri,minsub,maxsub,lenw,restar,lenver,neval,
     &        nsub,ifail,list(lenver),vacant(mdiv)
      double precision ver(2,3,lenver),epsabs,epsrel,result(numfun),
     &                 abserr(numfun),values(numfun,maxsub),
     &                 errors(numfun,maxsub),greate(maxsub),work2(lenw)

      integer i,j,k,l
      integer sbrgns,i1,i2,i3,i4,size
      integer l1
      integer ndiv,pointr,index,top
      double precision verold(2,3)

      if (restar.eq.1) then
          sbrgns = nsub
          go to 110
      end if
c
c   Initiate SBRGNS, RESULT, ABSERR and NEVAL.
c
      sbrgns = 0

      do j = 1,numfun
          result(j) = 0
          abserr(j) = 0
      end do

      neval = 0
c
c   Apply DRLTRI over the NUMTRI triangles.
c   This loop may be run in parallel.
c
      do i = 1,numtri
          l1 = 1 + (i-1)*8*numfun
          call drltri(ver(1,1,i),numfun,funsub,work2(l1),values(1,i),
     &                errors(1,i),greate(i))
          neval = neval + 37
          sbrgns = sbrgns + 1
      end do
c
c   Add the computed values to RESULT and ABSERR.
c
      do i = 1,numtri
          do j = 1,numfun
              result(j) = result(j) + values(j,i)
              abserr(j) = abserr(j) + errors(j,i)
          end do
      end do
c
c   Store results in heap.
c
      do i = 1,numtri
          index = i
          call dtrtri(2,index,greate,list,i)
      end do
c
c   We check for termination.
c
      if (sbrgns.lt.minsub) then
          go to 110
      end if

      do j = 1,numfun
          if (abserr(j).gt.epsrel*abs(result(j)) .and.
     &        abserr(j).gt.epsabs) then
              go to 110
          end if
      end do

      ifail = 0
      go to 499
c
c  End initiation.
c
c  Begin loop while the error is too great,
c  and SBRGNS+3 is less than MAXSUB.
c
  110 continue

      if (sbrgns+3.le.maxsub) then
c
c   If we are allowed to divide further,
c   prepare to apply basic rule over  the triangles produced
c   by dividing the
c   NDIV subtriangles with greatest errors.
c   If MAXSUB and SBRGNS are great enough, NDIV = MDIV.
c
          ndiv = maxsub - sbrgns
          ndiv = min(ndiv,mdiv,sbrgns)
c
c   Divide the NDIV subtriangles in four new subtriangles, and compute
c   integral and error over each.
c   When we pick a triangle to divide in four, then one of the new
c   sub-triangles is stored in the original triangle's position in the
c   datastructure. Thus we get 3*NDIV new elements in the datastructure
c   after the subdivision. The size of the datastructure before the
c   subdivision is stored in the variable SIZE, while SBRGNS is the
c   size of the heap at any time.
c   Hence.
c
          size = sbrgns

          do i = 1,ndiv

              pointr = size + 3* (ndiv+1-i)
c
c   Adjust RESULT and ABSERR. TOP is a pointer to the top of the heap.
c
              top = list(1)
              vacant(i) = top
              do j = 1,numfun
                  result(j) = result(j) - values(j,top)
                  abserr(j) = abserr(j) - errors(j,top)
              end do
c
c   Save the vertices.
c
              do l = 1,2
                  do k = 1,3
                      verold(l,k) = ver(l,k,top)
                  end do
              end do
c
c   Adjust the heap.
c
              call dtrtri(1,sbrgns,greate,list,k)
c
c   Compute the four new triangles.
c
              i1 = top
              i2 = pointr - 2
              i3 = pointr - 1
              i4 = pointr
              ver(1,1,i1) = verold(1,1)
              ver(2,1,i1) = verold(2,1)
              ver(1,2,i1) = (verold(1,1)+verold(1,2))/2
              ver(2,2,i1) = (verold(2,1)+verold(2,2))/2
              ver(1,3,i1) = (verold(1,1)+verold(1,3))/2
              ver(2,3,i1) = (verold(2,1)+verold(2,3))/2
              ver(1,1,i2) = ver(1,2,i1)
              ver(2,1,i2) = ver(2,2,i1)
              ver(1,2,i2) = verold(1,2)
              ver(2,2,i2) = verold(2,2)
              ver(1,3,i2) = (verold(1,2)+verold(1,3))/2
              ver(2,3,i2) = (verold(2,2)+verold(2,3))/2
              ver(1,1,i3) = ver(1,3,i1)
              ver(2,1,i3) = ver(2,3,i1)
              ver(1,2,i3) = ver(1,3,i2)
              ver(2,2,i3) = ver(2,3,i2)
              ver(1,3,i3) = verold(1,3)
              ver(2,3,i3) = verold(2,3)
              ver(1,1,i4) = ver(1,3,i2)
              ver(2,1,i4) = ver(2,3,i2)
              ver(1,2,i4) = ver(1,3,i1)
              ver(2,2,i4) = ver(2,3,i1)
              ver(1,3,i4) = ver(1,2,i1)
              ver(2,3,i4) = ver(2,2,i1)

          end do
c
c   Apply basic rule over 4*NDIV triangles.
c   This loop may be run in parallel.
c
          do i = 1,4*ndiv
              if (i.le.ndiv) then
                  index = vacant(i)
              else
                  index = sbrgns + i
              end if
              l1 = 1 + (i-1)*8*numfun
              call drltri(ver(1,1,index),numfun,funsub,work2(l1),
     &                    values(1,index),errors(1,index),greate(index))
          end do

          neval = neval + 4*ndiv*37
c
c   Add new contributions to RESULT and ABSERR.
c
          do i = 1,4*ndiv
              if (i.le.ndiv) then
                  index = vacant(i)
              else
                  index = sbrgns + i
              end if

              do j = 1,numfun
                  result(j) = result(j) + values(j,index)
                  abserr(j) = abserr(j) + errors(j,index)
              end do

          end do
c
c   Store results in heap.
c
          do i = 1,4*ndiv
              if (i.le.ndiv) then
                  index = vacant(i)
              else
                  index = sbrgns + i
              end if
              j = sbrgns + i
              call dtrtri(2,j,greate,list,index)
          end do

          sbrgns = sbrgns + 4*ndiv
c
c   Check for termination.
c
          if (sbrgns.lt.minsub) then
              go to 110
          end if

          do j = 1,numfun
              if (abserr(j).gt.epsrel*abs(result(j)) .and.
     &            abserr(j).gt.epsabs) then
                  go to 110
              end if
          end do

          ifail = 0
c
c   Else we did not succeed with the
c   given value of MAXSUB.
c
      else
          ifail = 1
      end if
c
c   Compute more accurate values of RESULT and ABSERR.
c
 499  continue

      do j = 1,numfun
          result(j) = 0
          abserr(j) = 0
      end do

      do i = 1,sbrgns
          do j = 1,numfun
              result(j) = result(j) + values(j,i)
              abserr(j) = abserr(j) + errors(j,i)
          end do
      end do

      nsub = sbrgns
      return
      end
