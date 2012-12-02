      program sumup

c*********************************************************************72

      real average
      real sum

      sum = 0
      sum = sum + 2
      sum = sum + 16
      sum = sum + 3

      average = sum / 3

      write (*,*) 'The sum is ', sum
      write (*,*) 'The average value is ', average

      stop
      end
