c CLASS = W
c  
c  
c  This file is generated automatically by the setparams utility.
c  It sets the number of processors and the class of the NPB
c  in this directory. Do not modify it by hand.
c  
        integer nx_default, ny_default, nz_default
        parameter (nx_default=128, ny_default=128, nz_default=128)
        integer nit_default, lm, lt_default
        parameter (nit_default=4, lm = 7, lt_default=7)
        integer debug_default
        parameter (debug_default=0)
        integer ndim1, ndim2, ndim3
        parameter (ndim1 = 7, ndim2 = 7, ndim3 = 7)
        integer one, nv, nr, ir
        parameter (one=1)
        logical  convertdouble
        parameter (convertdouble = .false.)
        character compiletime*11
        parameter (compiletime='28 May 2012')
        character npbversion*3
        parameter (npbversion='3.3')
        character cs1*8
        parameter (cs1='gfortran')
        character cs2*6
        parameter (cs2='$(F77)')
        character cs3*6
        parameter (cs3='(none)')
        character cs4*6
        parameter (cs4='(none)')
        character cs5*2
        parameter (cs5='-O')
        character cs6*2
        parameter (cs6='-O')
        character cs7*6
        parameter (cs7='randi8')
