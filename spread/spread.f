      program main

c*********************************************************************72
c
cc MAIN is the main program for SPREAD.
c
c  Discussion:
c
c    SPREAD analyzes all of the sequences in a sequence data
c    file to identify the tyrosine residues that are most likely to
c    undergo tryosine sulfation as a posttranslational modification.
c
c    Subroutine DBSite reads the sequences one at a time from the data
c    file and finds all of the tyrosine residues in the sequences.
c    The subroutine returns a each potential site as a subsequence of
c    up to 15 residues as a character string in array ChrSeq and as
c    a numerical list of residue indentities in array IntSeq.
c
c  Modified:
c
c    02 January 2009
c
      implicit none

      integer bprop
      integer eprop
      integer ksitep
      integer msitep
      integer naap

      parameter ( bprop = -7 )
      parameter ( eprop = 7 )
      parameter ( ksitep = 73 )
      parameter ( msitep = 1000 )
      parameter ( naap = 24 )

      character*1 aachar(24)
      integer aanumb( 0:127 )
      integer bwin
      character*15 chrseq(msitep)
      real cuts
      integer dbfile
      integer ewin
      logical fertig
      real grade(msitep)
      integer hitknt(ksitep)
      integer i
      integer intseq(bprop:eprop,msitep)
      integer ispos
      integer k
      integer ls
      integer ns
      integer nsites
      character*100 pir_file
      integer pirfl
      character*6 posid(0:ksitep)
      integer posn15
      character*15 postbl(0:ksitep)
      real profil(bprop:eprop,naap)
      integer result
      character*100 score_file
      character*100 seq_file
      character*6 seqac
      character*10 seqid
      character*100 site_file
      integer tyrnum(msitep)

      data  AAChar  / 'a', 'r', 'n', 'd', 'c', 'q', 'e', 'g', 'h',
     &                'i', 'l', 'k', 'm', 'f', 'p', 's', 't', 'w',
     &                'y', 'v', 'b', 'z', 'x', '-'  /

      data AANumb / 32 * -1,100,-1,-1,19,3*-1,3*0,200,-1,0, 24,  0,  0,
     &             -1, 300,302, 10 * -1, 0, -1, -1, -1,  1, 21,  5,  4,
     &              7,  14,  8,  9,  10, 0, 12, 11, 13,  3,  0, 15,  6,
     &              2,  16,  17,  0,  20, 18, 23, 19, 22, 6 * -1,    1,
     &              21,  5,  4,   7,  14,  8,  9, 10, 0, 12, 11, 13, 3,
     &              0,  15,  6,   2,  16, 17,  0, 20, 18, 23, 19, 22,
     &              5 * -1   /

      data  PosID / 'Null00', 'FICAN2', 'FICAN3', 'FIHOR3', 'FIRAB4',
     &              'FIPIG4', 'FILAM4', 'FIELE4', 'FITAP4', 'FICAP5',
     &              'LSLEU5', 'FIBIS6', 'LSLEU6', 'FIANA6', 'FISyN6',
     &              'FIBOV6', 'FIODO6', 'FIRAN6', 'FICER6', 'FIMUN6',
     &              'FIANT6', 'COB200', 'CM1416', 'SGH341', 'SGM348',
     &              'CM1413', 'F8H737', 'FAB200', 'G1H292', 'F8H738',
     &              'G1H295', 'CHIK47', 'CH1419', 'HEHU79', 'FIH444',
     &              'ITHI64', 'CAX174', 'ITHI61', 'ITHI63', 'G1H294',
     &              'F8H742', 'CCR111', 'CCH111', 'CCP110', 'AMR965',
     &              'CCFR41', 'CCRA63', 'CCHU97', 'HEHU92', 'FH1699',
     &              'GAFE12', 'GAHU87', 'GACA29', 'GARA87', 'CICIO2',
     &              'CM1417', 'FAH365', 'SGH151', 'SGR153', 'VID172',
     &              'CICIO3', 'CCTR46', 'CAX156', 'CCR113',
     &              'CCH113', 'CCP112', 'GAR103', 'FG1683', 'CH1417',
     &              'GABO12', 'FIPE13', 'A2H484', 'GACH28', 'GAOP28' /

      data (PosTbl(i),i=0,37)/ '               ', '------hyyddtdee',
     &      '-----hyyddtdeee', '-----ldydheeedg', '----addyddevlpd',
     &      '----aidydededgr', '----atdydeeeddr', '----atdyedeefpg',
     &      '----lsdydeeeder', '---gyldydevddnr', '---qsddyghmrf--',
     &      '--efptdydegeddr', '--eqfedyghmrf--', '--qastdyddedest',
     &      '--qfptdydegeddr', '--qfptdydegqddr', '--qhladydevdddr',
     &      '--qhladydeveddr', '--qhstdydeeeedr', '--qhstdydeveddr',
     &      '--qpsydydeeeddr', 'aaraeleyglvaeae', 'anedyedyydmpaad',
     &      'aseeepeygeeikgy', 'aseeepeygeesrsy', 'awdanedyedyydmp',
     &      'cdkntgdyyedsyed', 'daselehydpadlsp', 'degdtdlydyypeed',
     &      'dkntgdyyedsyedi', 'dtdlydyypeedteg', 'ealhdhfypdwmdf-',
     &      'eanedyeydelpakd', 'egeedddyldlekif', 'ehpaeteydslyped',
     &      'epipedayde-----', 'fadgqqdytgwmdfg', 'fdpipeeyls-----' /

      data (PosTbl(i),i=38,73) /
     &      'feeipeeylq-----', 'gdtdlydyypeedte', 'gdyyedsyedisayl',
     &      'grrsaedyeyps---', 'grrsaeeyeyps---', 'grrsaeeyeyts---',
     &      'gteseeeysaplpkp', 'hpmrdrdyagwmdf-', 'hrindrdymgwmdf-',
     &      'hrisdrdymgwmdfg', 'ifsedddyidivdsl', 'kkedfdiydedenqs',
     &      'leeeeaaygwmdfgr', 'leeeeeaygwmdfgr', 'meeeeaaygwmdf--',
     &      'meeeeeaygwmdfgr', 'mqrmdrnyygwmdfg', 'nedyedyydmpaadd',
     &      'nneeaedydddltds', 'pmdmsddyetqqwpe', 'pvdtpddyetqqwpe',
     &      'qpyettdysneeqsq', 'qrmdrnyygwmdfgk', 'rplhdhdypgwmdf-',
     &      'rrdgqqdytgwmdfg', 'rsaedyeyps-----',
     &      'rsaeeyeyps-----', 'rsaeeyeyts-----', 'saeeedqyn------',
     &      'sdqeeidyddtisve', 'tmeanedyeydelpa', 'veeeeaaygwmdf--',
     &      'vgqpendydtgddbt', 'vppmeedypqfgspk', 'waeeeaaygwmdf--',
     &      'wleeeeaygwmdf--'  /
c
c  Log-Odds position specific scoring matrix for tyrosine sulfation
c  sites.  The matrix is in units of (whole) bits.  The matrix was
c  computed using 44 known tyrosine sulfation sites in a variety of
c  peptides and proteins.  The background counts were derived from
c  sites around tyrosines in the same proteins and peptides which
c  were not found to accept tyrosine sulfation.  Pseudocounts to
c  fill zero cells in the table were derived from the amino acid
c  composition totaled over the positions from -7 to +7 amino acids
c  from the sulfated tyrosines.  A total of 13 pseudocounts were
c  added to each position (column) in the counts table.
c
c  The most effective range to use for identifying new sites with
c  this matrix is a window of from -5 (5 amino acids N-terminal
c  to the sulfated tyrosine) to +5 (5 amino acids C-terminal).
c
c  Alanine
c
      data (Profil(i,1), i = -7,7,1)/  0.266,  0.839, -0.662,  0.044,
     &                -0.378, -0.956, -0.609,  0.000, -1.199, -0.882,
     &                -0.880, -0.442,  0.665,  0.262, -0.818 /
c
c  Arginine
c
      data (Profil(i,2), i = -7,7,1)/ -0.439, -0.005,  0.057, -0.568,
     &                -1.106, -0.444, -3.063,  0.000, -3.067, -2.790,
     &                -3.107,  0.291, -2.485, -2.466,  0.318 /
c
c  Asparagine
c
      data (Profil(i,3), i = -7,7,1)/  0.936,  0.456,  0.676, -0.298,
     &                -0.044, -0.234, -0.617,  0.000, -1.307, -0.596,
     &                -1.559, -2.502, -0.718, -0.141, -2.580 /
c
c  Aspartic Acid
c
      data (Profil(i,4), i = -7,7,1)/  1.124,  0.915,  0.965,  1.469,
     &                 2.111,  2.373,  3.652,  0.000,  2.444,  2.369,
     &                 1.130,  1.532,  2.973,  1.493,  1.667 /
c
c  Cysteine
c
      data (Profil(i,5), i = -7,7,1)/ -0.577, -5.486, -5.010, -5.344,
     &                -5.064, -5.361, -4.690,  0.000, -4.694, -4.668,
     &                -5.285, -4.978, -4.559, -5.893, -5.402 /
c
c  Glutamine
c
      data (Profil(i,6), i = -7,7,1)/  0.196, -0.692,  0.408, -1.330,
     &                -1.281,  0.223, -0.846,  0.000, -1.847,  0.068,
     &                -0.634, -0.846, -1.088, -1.069, -0.710 /
c
c  Glutamic Acid
c
      data (Profil(i,7), i = -7,7,1)/  1.249,  1.716,  1.013,  1.925,
     &                 2.619,  2.233,  1.842,  0.000,  2.380,  2.416,
     &                 1.128,  1.473,  0.498,  1.036,  1.060 /
c
c  Glycine
c
      data (Profil(i,8), i = -7,7,1)/  0.945, -0.262,  0.008, -0.136,
     &                -0.903, -1.118, -2.518,  0.000,  0.802,  0.675,
     &                -0.950, -1.236, -1.766, -0.081,  0.568 /
c
c  Histidine
c
      data (Profil(i,9), i = -7,7,1)/  0.697, -0.708, -2.151,  0.958,
     &                -1.663,  0.573, -0.015,  0.000, -1.684,  0.656,
     &                -2.952, -3.386, -2.646, -2.464, -3.049 /
c
c  Isoleucine
c
      data (Profil(i,10),i = -7,7,1)/ -0.678, -2.243,  0.622, -0.895,
     &                -3.245, -1.014, -0.938,  0.000, -1.316, -2.930,
     &                 0.100, -0.300, -2.821,  0.057, -0.194 /
c
c  Leucine
c
      data (Profil(i,11),i = -7,7,1)/ -3.472, -3.224, -1.043, -2.124,
     &                -0.607, -0.482, -2.127,  0.000, -1.015, -2.185,
     &                -0.060, -0.329, -1.292, -2.800, -0.874 /
c
c  Lysine
c
      data (Profil(i,12),i = -7,7,1)/ -1.483, -0.026, -1.817, -3.789,
     &                -3.732, -4.226, -3.659,  0.000, -3.468, -4.036,
     &                -3.730, -3.006, -0.245, -0.553,  0.059 /
c
c  Methionine
c
      data (Profil(i,13),i = -7,7,1)/  1.046,  0.866,  0.967,  2.185,
     &                -0.260, -1.391, -1.205,  0.000, -0.028, -1.290,
     &                 1.846,  3.517, -1.782,  0.126, -0.038 /
c
c  Phenlyalanine
c
      data (Profil(i,14),i = -7,7,1)/ -0.341, -0.750, -2.652, -1.285,
     &                -0.473, -1.222, -1.000,  0.000, -3.323, -2.144,
     &                -0.198, -2.422,  0.505,  1.871, -0.062 /
c
c  Proline
c
      data (Profil(i,15),i = -7,7,1)/ -0.893,  0.441,  0.006, -0.021,
     &                 0.149, -0.264, -1.625,  0.000,  0.814, -0.132,
     &                -0.214,  1.181, -0.333,  0.700,  1.046 /
c
c  Serine
c
      data (Profil(i,16),i = -7,7,1)/ -0.550, -0.343, -0.427,  0.339,
     &                -1.505, -3.229, -0.964,  0.000, -0.502, -0.402,
     &                -1.299, -0.734, -0.293, -0.080, -0.258 /
c
c  Threonine
c
      data (Profil(i,17),i = -7,7,1)/ -0.772, -0.367, -1.024, -1.051,
     &                -0.023,  0.978, -2.630,  0.000, -0.199,  0.132,
     &                -0.344, -1.285, -0.153, -1.303, -0.487 /
c
c  Tryptophan
c
      data (Profil(i,18),i = -7,7,1)/  0.188,  1.750,  0.109, -1.323,
     &                -2.516, -2.038, -2.317,  0.000, -2.891,  0.585,
     &                 3.502, -1.264,  0.491, -1.217, -2.931 /
c
c  Tyrosine
c
      data (Profil(i,19),i = -7,7,1)/ -0.391, -1.327,  0.702,  1.019,
     &                 1.105,  0.405,  0.725,  0.000,  1.000,  0.400,
     &                 0.735,  1.743, -0.054,  0.110,  1.287 /
c
c  Valine
c
      data (Profil(i,20),i = -7,7,1)/  0.027, -3.830, -3.805, -4.387,
     &                -4.162, -4.968, -4.430,  0.000, -4.648, -4.483,
     &                 0.137, -2.035, -3.853, -1.702, -4.266 /
c
c  Asartic Acid or Asparagine  --  ambiguous (Asx)
c
      data (Profil(i,21),i = -7,7,1)/  1.070,  0.765,  0.886,  0.803,
     &                 1.042,  1.036,  1.443,  0.000,  0.416,  1.068,
     &                 0.487, -0.033,  1.072,  0.817, -0.245 /
c
c  Glutamic Acid or Glutamine  --  ambiguous (Glx)
c
      data (Profil(i,22),i = -7,7,1)/  0.827,  0.427,  0.853,  0.306,
     &                 0.381,  1.641,  0.599,  0.000,  0.535,  1.204,
     &                 0.443,  0.489, -0.193,  0.096,  0.370 /
c
c  Unknown or completely ambiguous amino acid  --  X
c
      data (Profil(i,23),i = -7,7,1)/ -0.420, -0.559, -0.468, -0.853,
     &                -0.929, -1.249, -1.415,  0.000, -1.178, -1.107,
     &                -0.716, -0.822, -0.750, -0.593, -0.666 /
c
c  Gap (beyond the N-terminal or C-terminal end of the sequence).
c
      Data (Profil(i,24),i = -7,7,1) /   15 * 0.000   /

    5 format( A1, 1X, A10, 1X, A6, I6, F9.3:, 2X, A6)
    6 format( ' ', A6, I6 )
    7 format(' The number of sites in the database that were found',
     &  /' to be identical to the specified known positive site.', / )
    8 format(//' Please enter a cutoff score - sequences for all',
     &         ' sites with scores',
     &        /' greater than or equal to this will be reported.'/)
   10 format( A1, 1X, A10, 1X, A6, I6, F9.3, 2X, A6, 2X, A15 )
   11 format( A1, 1X, A10, 1X, A6, I6, F9.3, 10X, A15 )
   12 format( '>F1;',A10, /'  ', A1, '  SeqID: ', A10,'  Ac Num: ',A6,
     &                     '  Tyr at:',I6, '  Score =' F9.3,
     &                    / A15 '*' )

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPREAD'

      dbfile = 10
      bwin = -5
      ewin = 5
      result = 2
      pirfl = 3

      seq_file = 'fibb.seq'
      score_file = 'fibb.score'
      pir_file = 'fibb.pir'
      site_file = 'fibb_identical.sites'
      
      do i = 1, ksitep
        hitknt(i) = 0
      end do

      write ( *, 8 )
      read ( *, * ) cuts

      open ( 
     &  unit = dbfile, 
     &  file = seq_file, 
     &  status = 'old' )

      call file_delete ( score_file )

      open ( 
     &  unit = result, 
     &  file = score_file, 
     &  status = 'new',
     &  access = 'sequential', 
     &  form = 'formatted' )

      call file_delete ( pir_file )

      open ( 
     &  unit = pirfl, 
     &  file = pir_file, 
     &  status = 'new',
     &  access = 'sequential', 
     &  form = 'formatted' )

  150 continue

      nsites = 0
      seqid = '          '
      seqac = '      '

      call dbsite ( dbfile, ls, bprop, eprop, nsites,
     &             msitep, intseq, chrseq, tyrnum, aanumb,
     &             aachar, seqid, seqac, fertig )

      if ( fertig ) then
        go to 1000
      end if
c
c  Process the sites
c
      call score ( bprop, eprop, naap, msitep, nsites, profil, intseq,
     &  grade, bwin, ewin )

      do ns = 1, nsites

        ispos = posn15 ( ksitep, chrseq(ns), postbl )

        if ( ispos .gt. 0 )    then

          write(result,10)  '+', seqid, seqac, tyrnum(ns),
     &                             grade(ns), posid(ispos),
     &                             chrseq(ns)
          write( pirfl, 12 )  seqid, '+', seqid, seqac, tyrnum(ns),
     &                          grade(ns), chrseq(ns)
          hitknt(ispos) = hitknt(ispos) + 1

        else

          if ( grade(ns) .ge. cuts )    then

            write( result, 11 )  '-', seqid, seqac, tyrnum(ns),
     &                                   grade(ns), chrseq(ns)
            write( pirfl, 12 )  seqid, '-', seqid,seqac,tyrnum(ns),
     &                                  grade(ns), chrseq(ns)

          else

            write(result,5)  '-', seqid, seqac,tyrnum(ns),grade(ns)

          end if

        end if

      end do

      go to 150

 1000 continue

      close ( unit = dbfile, status = 'keep' )
      close ( unit = result, status = 'keep' )
      close ( unit = pirfl,  status = 'keep' )

      call file_delete ( site_file )

      open ( 
     &  unit = result, 
     &  file = site_file,
     &  status = 'new',
     &  access = 'sequential', 
     &  form = 'formatted' )

      write ( result, 7 )

      do k = 1, ksitep
        write ( result, 6 )  posid(k), hitknt(k)
      end do

      close ( unit = result, status = 'keep' )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPREAD:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop ' swiss-prot sequences scored.'

      end
      subroutine dbsite ( dbfile, ls, bprop, eprop, nsites, msitep, 
     &  intseq, chrseq, tyrnum, aanumb, aachar, seqid, seqac, fertig )

c*********************************************************************72
c
cc DBSITE finds tyrosine residues in a sequence read from a file.
c
c  Modified:
c
c    02 January 2009
c
      implicit none

      integer bprop
      integer eprop
      integer msitep

      character*1 aachar(24)
      integer aanumb(0:127)
      character*15 chrseq( msitep )
      integer dbfile
      logical fertig
      integer i
      integer intseq( bprop:eprop, msitep )
      integer j
      integer k
      integer l
      integer lb
      integer le
      character*80 line
      integer ls
      character*1 lseq( 25000 )
      integer nl
      integer nlines
      integer nsites
      integer seq( 25000 )
      character*6 seqac
      character*10 seqid
      integer tyrnum( msitep )
      integer wbgn
      integer wend
      integer wp

    3 format(/' Accession line is missing after the following line:',
     &       /' :', A60 )
    4 format(/, 9X, 5( 1X, 10A1: ) )
c
c  Read file until an "ID" (identifier) line is detected and get
c  the 10 character identifier and the length
c
  100 continue

      read ( dbfile, '(a)', end = 1000 )  line
      fertig = .false.

      if ( line(1:2) .eq. 'ID' )    then
        seqid = line(6:15)
        read ( line(40:45), '(i6)' )  ls
      else
        go to 100
      end if
c
c  The next line should be the "AC" (accession number) line.  Save
c  the accession number since it is the only guarenteed means of
c  locating and identifying the sequence.
c
      read ( dbfile, '(a)' ) line

      if ( line(1:2) .eq. 'AC' ) then
        seqac = line( 6:12 )
      else
        write ( *, 3 ) line(1:60)
        write ( *, * ) ' bad accession number line in database.'
        stop
      end if
c
c  Skip lines (should be only one) until we get to the line marking
c  the start of the sequence data.
c
  110 continue

      read ( dbfile, '(a)' )  line

      if ( 
     &  index ( line, '..' ) .le. 0 .or.
     &  index ( line, 'Type: P' ) .le. 0 ) then
        go to 110
      end if
c
c  Now read the sequence
c
      nlines = ( ls + 49 ) / 50
      lb = -49

      do nl = 1, nlines
        lb = lb + 50
        le = lb + 49
        le = min ( le, ls )
        read ( dbfile, 4 )   ( lseq(l), l = lb, le, 1 )
      end do
c
c  Convert the sequence characters to a numerical equivalent.
c
      do l = 1, ls
        seq(l) = aanumb( ichar( lseq( l ) ) )
      end do
c
c  Find the potential tyrosine sulfation sites in the sequence.
c
      do l = 1, ls

        if ( seq(l) .eq. 19 )    then

          nsites = nsites + 1
          tyrnum( nsites ) = l
          chrseq( nsites ) = '---------------'
          wbgn = l + bprop
          j = bprop

          if ( wbgn .lt. 1 )    then
            j = bprop -( wbgn - 1 )
            wbgn = 1
          end if

          k = j + 7
          wend = l + eprop
          if( wend .gt. ls )    wend = ls

          do i = bprop, eprop, 1
            intseq( i, nsites ) = 24
          end do

          do wp = wbgn, wend, 1
            k = k + 1
            chrseq( nsites )( k:k ) = aachar( seq( wp ) )
            intseq( j, nsites ) = seq( wp )
            j = j + 1
          end do

        end if

      end do

      return

 1000 continue

      fertig = .true.

      return
      end
      subroutine file_delete ( file_name )

c*********************************************************************72
c
cc FILE_DELETE deletes a file if it exists.
c
c  Discussion:
c
c    You might want to call this routine to get rid of any old copy
c    of a file, before trying to open a new copy with the OPEN argument:
c
c      status = 'new'.
c
c    It's not always safe to open a file with " STATUS = 'UNKNOWN' ".
c
c    For instance, on the SGI, the most recent version of the FORTRAN
c    compiler seems to go crazy when I open an unformatted direct
c    access file this way.  It creates an enormous file (of somewhat
c    random size).  The problem goes away if I delete any old copy
c    using this routine, and then open a fresh copy with
c    " STATUS = 'NEW' ".  It's a scary world.
c
c  Modified:
c
c    25 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*(*) FILE_NAME, the name of the file to be deleted.
c
      implicit none

      character*80 ctemp
      character*(*) file_name
      logical lfile
      integer s_len_trim
      integer unit
c
c  Does the file exist?
c
      inquire (
     &  file = file_name,
     &  exist = lfile )

      if ( .not. lfile ) then
        return
      end if
c
c  Can we get a FORTRAN unit number?
c
      call get_unit ( unit )

      if ( unit .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FILE_DELETE - Error!'
        write ( *, '(a)' ) '  A free FORTRAN unit could not be found.'
        return
      end if
c
c  Can we open the file?
c  
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FILE_DELETE:'
      write ( *, '(a)' ) '  Open "' // trim ( file_name ) // '".'

      open (
     &  unit = unit,
     &  file = file_name,
     &  status = 'old',
     &  err = 10 )
c
c  Can we close the file with "Delete" status?
c
      write ( *, '(a)' ) '  Delete "' // trim ( file_name ) // '".'

      close (
     &  unit = unit,
     &  status = 'delete',
     &  err = 20 )

      return

10    continue

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FILE_DELETE - Error!'
      write ( *, '(a)' ) 
     &  '  Could not open "' // trim ( file_name ) // '".'
      return

20    continue

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FILE_DELETE - Error!'
      write ( *, '(a)' ) 
     &  '  Could not delete "' // trim ( file_name ) // '".'
      return

      end
      subroutine get_unit ( unit )

c*********************************************************************72
c
cc GET_UNIT returns a free FORTRAN unit number.
c
c  Discussion:
c
c    A "free" FORTRAN unit number is an integer between 1 and 99 which
c    is not currently associated with an I/O device.  A free FORTRAN unit
c    number is needed in order to open a file with the OPEN command.
c
c    If UNIT = 0, then no free FORTRAN unit could be found, although
c    all 99 units were checked (except for units 5, 6 and 9, which
c    are commonly reserved for console I/O).
c
c    Otherwise, UNIT is an integer between 1 and 99, representing a
c    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
c    are special, and will never return those values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 March 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer UNIT, the free unit number.
c
      implicit none

      integer i
      integer unit

      unit = 0

      do i = 1, 99

        if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

          open ( unit = i, err = 10 )
          close ( unit = i )

          unit = i

          return
        end if

10      continue

      end do

      return
      end
      function posn15 ( size, string, table )

c*********************************************************************72
c
cc POSN15 compares a data file site to a table of known tyrosine sulfation sites.
c
c  Discussion:
c
c    This routine compares the site from the data file to a
c    table of sites that are experimentally known to be tyrosine
c    sulfation sites.  
c
c    This result, in variable IsPos, is used to label the output for ROC 
c    plot analysis.  
c
c    Sites that are identical to one of the sites in array PosTbl (Positives Table) 
c    are labeled with a plus sign, "+", and those that do not match any of the
c    known sites are lableled with a minus sign, "-".  This output is
c    used as input for the ROC plot program.
c
c    POSN15 is a binary search of the table STDPOS.
c
c    search algorithm is an implementation of Knuth's algorithm U
c    (vol 3., ch. 6.2.1, p 411; see vol 1., ch. 1.2.4, p 37 for
c    details of the notation).
c
c    INDEX = The current position in the search table.
c
c    RATE  = The change in INDEX if a match is not found on the current
c    cycle.  If RATE becomes zero the value passed to the
c    function is not in the search table.
c
c  Modified:
c
c    02 January 2009
c
      implicit none

      integer size

      integer indx
      integer posn15
      integer rate
      character*15 string
      character*15 table(0:size)

      rate = size  / 2
      indx = ( size + 1 ) / 2

10    continue

      if ( llt ( string, table(indx) ) )    then
        indx = indx - ( ( rate + 1 ) / 2 )
      else if ( lgt ( string, table(indx) ) )    then
        indx = indx + ( ( rate + 1 ) / 2 )
      else if ( string .eq. table(indx) )    then
        posn15 = indx
        return
      end if

      if ( rate .eq. 0 )    then
        posn15 = 0
        return
      end if

      rate = rate / 2

      go to 10

      end
      function s_len_trim ( s )

c*********************************************************************72
c
cc S_LEN_TRIM returns the length of a string to the last nonblank.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 March 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*(*) S, a string.
c
c    Output, integer S_LEN_TRIM, the length of the string to the last nonblank.
c
      implicit none

      integer i
      character*(*) s
      integer s_len_trim

      do i = len ( s ), 1, -1

        if ( s(i:i) .ne. ' ' ) then
          s_len_trim = i
          return
        end if

      end do

      s_len_trim = 0

      return
      end
      subroutine score ( bprop, eprop, naap, maxsqp, nsites, profil,
     &  seq, grade, bwin, ewin )

c*********************************************************************72
c
cc SCORE scores each potential site in the array.
c
c  Discussion:
c
c    SCORE scores each potential site in array IntSeq using
c    the position specific score stored in array Profil and returns
c    a score in array Grade that reflects how similar the site is to
c    a set of sites that are experimentally known to be tyrosine
c    sulfation sites.
c
c  Modified:
c
c    02 January 2009
c
      implicit none

      integer bprop
      integer eprop
      integer maxsqp

      integer bwin
      integer ewin
      real grade( maxsqp )
      integer naap
      integer ns
      integer nsites
      integer pos
      real profil( bprop:eprop, naap )
      integer seq( bprop:eprop, maxsqp )

      do ns = 1, nsites
        grade(ns) = 0.0
        do pos = bwin, ewin
          grade(ns) = grade(ns) + profil( pos, seq(pos,ns) )
        end do
      end do

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
