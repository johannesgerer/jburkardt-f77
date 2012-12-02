      program main

c*********************************************************************72
c
cc MAIN is the main program for LAWSON_PRB4,
c
C  DEMONSTRATE SINGULAR VALUE ANALYSIS. 
c
c  Modified:
c
c    19 October 2008
c
c  Author:
c
c    Charles Lawson, Richard Hanson
c
      integer MDA, MX
      parameter(MDA = 15, MX = 5)
      integer I, J, KPVEC(4), M, N
      double precision A(MDA,MX), B(MDA), D(MX), SING(MX), WORK(2*MX)
      character NAMES(MX)*8
      data KPVEC / 1, 111111, -1, 76 /
      data NAMES/'  Fire',' Water',' Earth','   Air','Cosmos'/

      M = MDA
      N = MX

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LAWSON_PRB4'

      write (*,'(a)')
     * ' DEMONSTRATE SINGULAR VALUE ANALYSIS',
     * ' Read data from file: lawson_prb4_input.txt'
      open(unit= 10, file= 'lawson_prb4_input.txt', status= 'old')
      read (10,'(6F12.0)') ((A(I,J),J=1,N),B(I),I=1,M)  
      write (*,'(/a)')
     * ' LISTING OF INPUT MATRIX, A, AND VECTOR, B, FOLLOWS..'
      write (*,'(/(1x,f11.7,4F12.7,F13.4))') ((A(I,J),J=1,N),B(I),I=1,M)
      write (*,'(/)')  

      call SVA (A,MDA,M,N,MDA,B,SING,KPVEC, NAMES,1,D, WORK)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LAWSON_PRB4:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      end  
