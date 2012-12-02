      program hello2

c*********************************************************************72

      integer i
      integer imax

      i=0
      imax=12

10    continue

      i=i+1

      write(*,*)'Hello, person number ',i
      if (i.lt.imax) go to 10

      write(*,*)'Goodbye!'

      stop
      end
