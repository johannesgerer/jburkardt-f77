      program mix

c*********************************************************************72

      integer i
      real x
      real z

      i = 2
      j = 3.5

      x = 2
      y = 3.5

      write (*,*) 'I=', i, ' and X=', x

      i = 2.3
      x = 2.3

      write (*,*) 'I=', i, ' and X=', x

      i = y
      x = j

      write (*,*) 'I=', i, ' and X=', x

      i = j
      x = y

      write (*,*) 'I=', i, ' and X=', x

      i = 2.3 * j * y
      x = 2.3 * j * y

      write (*,*) 'I=', i, ' and X=', x

      stop
      end
