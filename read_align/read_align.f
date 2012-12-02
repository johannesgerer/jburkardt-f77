C---------------------------------------------------------------------
C Des Higgins, EMBL Data Library.      APRIL 1992                    I
C---------------------------------------------------------------------
C Program to TRY to read multiple alignments and write them out in
C PIR format.
C Optionally read in pairs of numbers in free format from a file to
C specify fragments.
        program readal3
        implicit integer(a-z)
        parameter (MAXL=5000, MAXN=500, LINELEN=200,NAMELEN=16)
        character seqs(maxn,maxl), names(MAXN)*(NAMELEN)
        character fname*100,line*(linelen),ch,linbuff*(linelen)
        logical blankline,segments

        lin = 44
        lout= 45
        lseg= 46
        print*
        print*
1       fname = ' '
        write(*,'(a$)') ' Input file ? (alignment    ) >> '
        read(*,'(a)') fname
        open(unit=lin,file=fname,status='old',err=1)

2       fname = ' '
        print*
        write(*,'(a$)') ' Output file ? (PIR format) >> '
        read(*,'(a)') fname
        open(unit=lout,file=fname,status='new',err=2)

        print*
4       fname = ' '
        write(*,'(a)')  ' Press <RETURN> to ignore this:'
        write(*,'(a)')  ' Use segments ? '
        write(*,'(a$)') ' File of fragment numbers     >> '
        read(*,'(a)') fname
        if(fname.ne.' ') then
                open(unit=lseg,file=fname,status='old',err=4)
                segments = .true.
        end if

C skip the first few lines of text ..... read until completely blank line
3       line = ' '
        read(lin,'(a)',end=50) line
        if(.not.blankline(line)) go to 3

        seqlen = 0
C skip any more blank lines
5       line = ' '
        read(lin,'(a)',end=50) line
        if(blankline(line)) go to 5

        backspace(lin)
        nseqs  = 0
C now read a block of alignment
        line = ' '
7       read(lin,'(a)',end=50) line
        if(blankline(line)) then
          seqlen = seqlen + found
          go to 5
        end if
        found = 0
        NSEQS = NSEQS + 1

        do 9 i = linelen,1,-1
        if(line(i:i).ne.' ') then
          blockend = i
          go to 13
        end if
9       continue
13      continue

c find first non blank character (start of name)
        do 20 i = 1,linelen
        if(line(i:i).ne.' ') then
                namestart = i
                go to 21
        end if
20      continue
21      continue

c find first blank character after the name
        do 22 i = namestart,linelen
        if(line(i:i).eq.' ') then
                namend = i
                go to 23
        end if
22      continue
23      continue

        do 10 i = blockend,namend,-1
        ch = line(i:i)
        if(ch.eq.' ') go to 10
        found = found + 1
        linbuff(found:found) = ch
10      continue
11      continue

        do 12 i = 1,found
12      seqs(nseqs,seqlen+i) = linbuff(found-i+1:found-i+1)
        if(seqlen.eq.0) names(nseqs) = line(1:min(namelen,namend))
        go to 7


50      print*
        print*,' No. of seqs. read = ',nseqs
        print*,' Alignment length  = ',seqlen

        do 60 i = 1,nseqs
        write(lout,'(a,a)') '>P1;',names(i)
        write(lout,'(a)') ' '
        if(.not.segments) then
                write(lout,'(1x,60a1)') (seqs(i,j),j=1,seqlen)
        else
                rewind(lseg)
51              read(lseg,*,end=52) segfrom, segto
                write(lout,'(1x,60a1)') (seqs(i,j),j=segfrom,segto)
                go to 51
        end if
52      continue
        write(lout,'(a)') '*'
60      continue

        stop
        end

C-----------------------------------------------------------------------
        logical function blankline(line)
C-----------------------------------------------------------------------
C Looks for a blank line i.e. a line that contains only spaces, dots,
C asterisks or digits.  ONLY WORKS WITH ASCII!!!!!!!
        implicit integer(a-z)
        character line*(*),ch

        blankline = .false.
        linlen = len(line)

        do 1 i = 1,linlen
        ch = line(i:i)
        if((ch.eq.' ').or.(ch.eq.'.').or.(ch.eq.'*')) go to 1
C ASCII only!! Look for digits.
        ich = ichar(ch)
        if((ich.ge.48).and.(ich.le.57)) go to 1
        return
1       continue

        blankline = .true.
        return
        end




