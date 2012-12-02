      Program Topcat

c*********************************************************************72

      Write (*,*) 'This is Topcat speaking!  Is Garfield home?'
      Call Garfield

      Write (*,*) 'This is Topcat again.  Is Underdog around?'
      Call Underdog

      Write (*,*) 'This is Topcat again.  I''m ready for a nap!'

      Stop
      End

      Subroutine Garfield

c*********************************************************************72

      Write (*,*) 'Garfield speaking.  I''m in charge now!'

      Call Underdog

      Write (*,*) 'This is Garfield, and I''m stepping down now!'

      Return
      End

      Subroutine Underdog

c*********************************************************************72

      Write (*,*) 'Underdog here!  Leave me alone!'

      Return
      End
