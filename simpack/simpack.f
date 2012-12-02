      subroutine smpint ( nd, nf, mincls, maxcls, funsub, epsabs, 
     &  epsrel, key, sbrgns, wrklen, vrtwrk, restar, value, error, 
     &  funcls, inform )

c*********************************************************************72
c
cc SMPINT integrates a vector function over a collection of simplexes.
c
c***BEGIN PROLOGUE SMPINT
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Alan Genz
c    Department of Mathematics
c    Washington State University
c    Pullman, WA 99164-3113, USA
c
c***KEYWORDS automatic multidimensional integrator,
c            n-dimensional simplex, general purpose, global adaptive
c***PURPOSE  To calculate an approximation to a vector of integrals
c
c               I = I (F ,F ,...,F   ) DS
c                    S  1  2      NF
c
c            where S is a collection ND-dimensional simplices,
c            and  F = F (X ,X ,...,X  ), K = 1, 2, ..., NF,
c                  K   K  1  2      ND
c            and try to satisfy for each component I(K) of I
c              ABS( I(K) - VALUE(K) ) < MAX( EPSABS, EPSREL*ABS(I(K)) )
c
c***DESCRIPTION Computation of integrals over simplical regions.
c            SMPINT is a driver for the integration routine SMPSAD,
c            which repeatedly subdivides the region of integration and
c            estimates the integrals and the errors over the subregions
c            with greatest estimated errors until the error request
c            is met or MAXCLS function evaluations have been used.
c
c   ON ENTRY
c
c     ND     Integer, number of variables. 1 < ND
c     NF     Integer, number of components of the integral.
c     MINCLS Integer, minimum number of FUNSUB calls.
c     MAXCLS Integer, maximum number of FUNSUB calls.
c            RULCLS is number FUNSUB calls for each subregion (see WRKLEN),
c            DIFCLS = 1 + 2*ND*( ND + 1 ).
c            If RESTAR = 0, MAXCLS must be >= MAX(SBRGNS*RULCLS,MINCLS).
c            If RESTAR = 1, MAXCLS must be >= MAX(4*RULCLS+DIFCLS,MINCLS).
c     FUNSUB Externally declared subroutine for computing components of
c            the integrand at the given evaluation point.
c            It must have parameters (ND,X,NF,FUNVLS)
c            Input parameters:
c              ND   Integer that gives the dimension of I
c              X      Real array of dimension ND that contains the
c                     evaluation point.
c              NF Integer that gives the number of components of I.
c            Output parameter:
c              FUNVLS Real array of dimension NF that contains the
c                     components of the integrand.
c     EPSABS Real.
c            Requested absolute accuracy.
c     EPSREL Real requested relative accuracy.
c     KEY    Integer, key to selected local integration rule.
c            KEY = 3 gives the user a (default) degree 7 integration rule.
c            KEY = 1 gives the user a degree 3 integration rule.
c            KEY = 2 gives the user a degree 5 integration rule.
c            KEY = 3 gives the user a degree 7 integration rule.
c            KEY = 4 gives the user a degree 9 integration rule.
c     WRKLEN Integer, length of the working array VRTWRK.
c             WRKLEN should be >= WRKSBS*( ND*(ND+1) + 2*NF + 3 )
c                                + (ND+1)*(ND+2) + 7*NF,      where
c             WRKSBS = SBRGNS + 3*( MAXCLS/RULCLS - SBRGNS*(1-RESTAR) )/4.
c            If
c              KEY = 0, RULCLS = (ND+4)*(ND+3)*(ND+2)/6 + (ND+2)*(ND+1);
c              KEY = 1, RULCLS = 2*ND+3;
c              KEY = 2, RULCLS = (ND+3)*(ND+2)/2 + 2*(ND+1);
c              KEY = 3, RULCLS = (ND+4)*(ND+3)*(ND+2)/6 + (ND+2)*(ND+1);
c              KEY = 4, RULCLS = (ND+5)*(ND+4)*(ND+3)*(ND+2)/24
c                                 + 5*(ND+2)*(ND+1)/2 .
c     SBRGNS Integer, initial number of simplices.
c     VRTWRK Real array of dimension WRKLEN.
c            Work should contain the simplex vertices for SBRGNS
c            simplices; the coordinates of vertex J for simplex K
c            must be in VRTWRK(I+J*ND+(K-1)*ND*(ND+1)),
c            for I = 1, ..., ND; J = 0, ..., ND; K = 1, ..., SBRGNS.
c            The rest of VRTWRK is used as working storage; see below.
c     RESTAR Integer.
c            If RESTAR = 0, this is the first attempt to compute
c             the integrals over SBRGNS simplices.
c            If RESTAR = 1, then we restart a previous attempt.
c             In this case the only parameters for SMPINT that may
c             be changed (with respect to the previous call of SMPINT)
c             are MINCLS, MAXCLS, EPSABS, EPSREL, KEY and RESTAR.
c   ON RETURN
c
c     SBRGNS Integer.
c            SBRGNS contains the current number of simplices. They
c            were obtained by subdividing the input simplicies.
c     VRTWRK   Real array of dimension WRKLEN.
c            Used as working storage.
c            Let MAXSUB = (WRKLEN-(ND+1)*(ND+2)-7*NF)/(ND*(ND+1)+2*NF+3).
c            VRTWRK(1), ..., VRTWRK(ND*(ND+1)), ...,
c             VRTWRK(ND*(ND+1)*SBRGNS) contain subregion vertices.
c            VRTWRK(ND*(ND+1)*MAXSUB+1), ...,
c             VRTWRK(ND*(ND+1)*MAXSUB+NF*SBRGNS) contain
c             estimated components of the integrals over the subregions.
c            VRTWRK(ND*(ND+1)*MAXSUB+NF*MAXSUB+1), ...,
c             VRTWRK(ND*(ND+1)*MAXSUB+NF*(MAXSUB+SBRGNS))
c             contain estimated errors for the subregions.
c            VRTWRK(ND*(ND+1)*MAXSUB+2*NF*MAXSUB+1), ...,
c             VRTWRK(ND*(ND+1)*MAXSUB+2*NF*MAXSUB+SBRGNS))
c             contain volumes of the subregions (scaled by ND!).
c            VRTWRK(ND*(ND+1)*MAXSUB+(2*NF+1)*MAXSUB+1), ...,
c             VRTWRK(ND*(ND+1)*MAXSUB+(2*NF+1)*MAXSUB+SBRGNS)
c             contain greatest errors in each subregion.
c            VRTWRK(ND*(ND+1)*MAXSUB+(2*NF+2)*MAXSUB+1), ...,
c             VRTWRK(ND*(ND+1)*MAXSUB+(2*NF+2)*MAXSUB+SBRGNS)
c             contain pointers for the subregion heap.
c            The rest of VRTWRK is used as temporary storage in SMPSAD.
c     VALUE  Real array of dimension NF of integral approximations.
c     ERROR  Real array of dimension NF, of absolute accuracy estimates.
c     FUNCLS Integer, number of FUNSUB calls used by SMPINT.
c     INFORM Integer.
c            INFORM = 0 for normal exit, when ERROR(K) <=  EPSABS or
c              ERROR(K) <=  ABS(VALUE(K))*EPSREL with MAXCLS or less
c              function evaluations for all values of K, 1 <= K <= NF.
c            INFORM = 1 if MAXCLS was too small for SMPINT to obtain
c              the required accuracy. In this case SMPINT returns
c              values VALUE with estimated absolute accuracies ERROR.
c            INFORM = 2 if KEY < 0 or KEY > 4,
c            INFORM = 3 if ND < 2,
c            INFORM = 4 if NF < 1,
c            INFORM = 5 if EPSABS < 0 and EPSREL < 0,
c            INFORM = 6 if WRKLEN is too small,
c            INFORM = 7 if RESTAR < 0 or RESTAR > 1,
c            INFORM = 8 if SBRGNS <= 0.
c
c***ROUTINES CALLED SMPCHC,SMPSAD
c***END PROLOGUE SMPINT
c
c   Global variables.
c
      external funsub
      integer nd, nf, mincls, maxcls, sbrgns
      integer key, wrklen, restar, funcls, inform
      double precision epsabs, epsrel
      double precision value(nf), error(nf), vrtwrk(wrklen)
c
c   Local variables.
c
c   MAXSUB Integer, maximum allowed number of subdivisions
c          for the given values of KEY, ND and NF.
c
      integer maxsub, rulcls, i1, i2, i3, i4, i5, i6
c
c***FIRST PROCESSING STATEMENT SMPINT
c
c  Compute MAXSUB and RULCLS, and check the input parameters.
c
c
      call smpchc ( nd, nf, mincls, maxcls, epsabs, epsrel, sbrgns,
     &  key, wrklen, restar, rulcls, maxsub, inform )

      if ( inform .eq. 0 ) then
c
c  Split up the work space and call SMPSAD.
c
         i1 =  1 + maxsub*nd*(nd+1)
         i2 = i1 + maxsub*nf
         i3 = i2 + maxsub*nf
         i4 = i3 + maxsub
         i5 = i4 + maxsub
         i6 = i5 + maxsub
         call smpsad( nd, nf, funsub, mincls, maxcls, epsabs, epsrel,
     &        restar, key, rulcls, maxsub, sbrgns, vrtwrk,
     &        vrtwrk(i1), vrtwrk(i2), vrtwrk(i3), vrtwrk(i4),
     &        vrtwrk(i5), vrtwrk(i6), value, error, funcls, inform )
      else
         funcls = 0
      end if

      return
      end
      subroutine smpchc ( nd, nf, mincls, maxcls, epsabs, epsrel, 
     &  sbrgns, key, wrklen, restar, rulcls, maxsub, inform )

c*********************************************************************72
c
cc SMPCHC checks the input parameters to SMPINT.
c
c***BEGIN PROLOGUE SMPCHC
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Alan Genz
c    Department of Mathematics
c    Washington State University
c    Pullman, WA 99164-3113, USA
c
c***PURPOSE  SMPCHC checks validity of input parameters for SMPINT.
c***DESCRIPTION
c            SMPCHC computes MAXSUB, RULCLS and INFORM as functions of
c             input parameters for SMPINT, and checks the validity of
c             input parameters for SMPINT.
c
c   ON ENTRY
c
c     ND   Integer, number of variables,  ND > 1.
c     NF Integer, number of components of the integral.
c     MINCLS Integer, minimum number of new FUNSUB calls.
c     MAXCLS Integer, maximum number of new FUNSUB calls.
c            The number of function values for each subregion is RULCLS.
c            If
c             KEY = 0, RULCLS = (ND+4)*(ND+3)*(ND+2)/6 + (ND+2)*(ND+1);
c             KEY = 1, RULCLS = 2*ND+3;
c             KEY = 2, RULCLS = (ND+3)*(ND+2)/2 + 2*(ND+1);
c             KEY = 3, RULCLS = (ND+4)*(ND+3)*(ND+2)/6 + (ND+2)*(ND+1);
c             KEY = 4, RULCLS = (ND+5)*(ND+4)*(ND+3)*(ND+2)/24
c                               + 5*(ND+2)*(ND+1)/2 .
c            DIFCLS = 1 + 2*ND*( ND + 1 ).
c            If RESTAR = 0, MAXCLS must be >= MAX(SBRGNS*RULCLS,MINCLS).
c            If RESTAR = 1, MAXCLS must be >= MAX(4*RULCLS+DIFCLS,MINCLS).
c     EPSABS Real, requested absolute accuracy.
c     EPSREL Real, requested relative accuracy.
c     SBRGNS Integer, initial number of simplices.
c     KEY    Integer, key to selected local integration rule.
c            KEY = 0 gives the user a (default)degree 7 integration rule.
c            KEY = 1 gives the user a degree 3 integration rule.
c            KEY = 2 gives the user a degree 5 integration rule.
c            KEY = 3 gives the user a degree 7 integration rule.
c            KEY = 4 gives the user a degree 9 integration rule.
c     WRKLEN Integer, length of the working array WORK.
c             WRKLEN should be >= WRKSBS*( ND*(ND+1) + 2*NF + 3 )
c                                + (ND+1)*(ND+2) + 7*NF, where
c             WRKSBS = SBRGNS + 3*( MAXCLS/RULCLS - SBRGNS*(1-RESTAR) )/4.
c     RESTAR Integer.
c            If RESTAR = 0, this is the first attempt to compute
c             the integral over the SBRGNS input simplices.
c            If RESTAR = 1, then we restart a previous attempt.
c
c   ON RETURN
c
c     RULCLS Integer, number of function values for each subregion.
c            If
c             KEY = 0, RULCLS = (ND+4)*(ND+3)*(ND+2)/6 + (ND+2)*(ND+1);
c             KEY = 1, RULCLS = 2*ND+3;
c             KEY = 2, RULCLS = (ND+3)*(ND+2)/2 + 2*(ND+1);
c             KEY = 3, RULCLS = (ND+4)*(ND+3)*(ND+2)/6 + (ND+2)*(ND+1);
c             KEY = 4, RULCLS = (ND+5)*(ND+4)*(ND+3)*(ND+2)/24
c                               + 5*(ND+2)*(ND+1)/2 .
c     MAXSUB Integer, maximum allowed number of subregions for the
c            given values of MAXCLS, WRKLEN, KEY and ND.
c     INFORM Integer.
c            INFORM = 0 for normal exit,
c            INFORM = 2 if KEY < 0 or KEY > 4,
c            INFORM = 3 if ND < 2,
c            INFORM = 4 if NF < 1,
c            INFORM = 5 if EPSABS < 0 and EPSREL < 0.,
c            INFORM = 6 if WRKLEN is too small,
c            INFORM = 7 if RESTAR < 0 or RESTAR > 1,
c            INFORM = 8 if SBRGNS <= 0.
c
c***END PROLOGUE SMPCHC
c
c   Global variables.
c
      integer nd, nf, mincls, maxcls, key, maxsub
      integer wrklen, inform, restar, rulcls, sbrgns
      double precision epsabs, epsrel
c
c     Local variables.
c
      integer wrkdif, difcls
c
c***FIRST PROCESSING STATEMENT SMPCHC
c
      inform = 0
c
c  Check valid KEY.
c
      if ( key .lt. 0 .or. key .gt. 4 ) inform = 2
c
c  Check valid ND.
c
      if ( nd .lt. 2 ) inform = 3
c
c  Check positive NF.
c
      if ( nf .lt. 1 ) inform = 4
c
c  Check valid accuracy requests.
c
      if ( epsabs .lt. 0 .and. epsrel .lt. 0 ) inform = 5
c
c  Check workspace.
c
      wrkdif = (nd+1)*(nd+2) + 7*nf
      maxsub = ( wrklen - wrkdif )/( (nd+1)*nd + 2*nf + 3 )
      if ( maxsub .le. sbrgns ) inform = 6
c
c  Check valid RESTAR.
c
      if ( restar .ne. 0 .and. restar .ne. 1 ) inform = 7
c
c  Check valid SBRGNS.
c
      if ( sbrgns .le. 0 ) inform = 8
c
c  Compute RULCLS as a function of KEY and ND and check MAXCLS.
c
      if ( inform .eq. 0 ) then
         difcls = 1 + 2*nd*( nd + 1 )
         if (key .eq. 0) rulcls = (nd+4)*(nd+3)*(nd+2)/6 + (nd+2)*(nd+1)
         if (key .eq. 1) rulcls = 2*nd + 3
         if (key .eq. 2) rulcls = (nd+3)*(nd+2)/2 + 2*(nd+1)
         if (key .eq. 3) rulcls = (nd+4)*(nd+3)*(nd+2)/6 + (nd+2)*(nd+1)
         if (key .eq. 4) rulcls = (nd+5)*(nd+4)*(nd+3)*(nd+2)/24
     &                           + 5*(nd+2)*(nd+1)/2
         if ( restar.eq.0 .and. maxcls.lt.max(sbrgns*rulcls,mincls) .or.
     &        restar.eq.1 .and. maxcls.lt.max(4*rulcls+difcls,mincls) )
     &        inform = 1
      end if

      return
      end
      subroutine smpsad ( nd, nf, funsub, mincls, maxcls, epsabs, 
     &  epsrel, restar, key, rulcls, maxsub, sbrgns, vertcs, values, 
     &  errors, volums, greats, pontrs, work, value, error, funcls, 
     &  inform )

c*********************************************************************72
c
cc SMPSAD integrates a vector function over a hyperrectangle.
c
c***BEGIN PROLOGUE SMPSAD
c***KEYWORDS automatic multidimensional integrator,
c            n-dimensional simplex,
c            general purpose, global adaptive
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Alan Genz
c    Department of Mathematics
c    Washington State University
c    Pullman, WA 99164-3113, USA
c
c***PURPOSE  The routine calculates an approximation to a given
c            vector of definite integrals, I, over a hyper-rectangular
c            region hopefully satisfying for each component of I the
c            following claim for accuracy:
c            ABS( I(K) - VALUE(K) ) .LE. MAX( EPSABS, EPSREL*ABS(I(K) ) )
c***DESCRIPTION Computation of integrals over hyper-rectangular regions.
c            SMPSAD repeatedly subdivides the regions of integration
c            and estimates the integrals and the errors over the
c            subregions with  greatest estimated errors until the error
c            request is met or MAXSUB subregions are stored. The regions
c            are divided into three or four equally sized parts along
c            the direction(s) with greatest absolute fourth difference.
c
c   ON ENTRY
c
c     ND     Integer, number of variables, ND > 1.
c     NF     Integer, number of components of the integral.
c     FUNSUB Externally declared subroutine for computing components of
c            the integrand at the given evaluation point.
c            It must have parameters (ND,X,NF,FUNVLS)
c            Input parameters:
c              ND Integer that gives the dimension of I
c              X  Real array of dimension ND that contains the
c                     evaluation point.
c              NF Integer that gives the number of components of I.
c            Output parameter:
c              FUNVLS Real array of dimension NF that contains the
c                     components of the integrand.
c     MINCLS Integer.
c            The computations proceed until there are at least
c            MINCLS FUNSUB calls.
c     MAXCLS Integer.
c            The computations proceed until further subdivision would
c            require more than MAXCLS FUNSUB calls. When RESTAR = 1,
c            this is the number of new FUNSUB calls.
c     EPSABS Real, requested absolute accuracy.
c     EPSREL Real, requested relative accuracy.
c     RESTAR Integer.
c            If RESTAR = 0, this is the first attempt to compute
c             the integral.
c            If RESTAR = 1, then we restart a previous attempt.
c             (In this case the output parameters from SMPSAD
c             must not be changed since the last exit from SMPSAD.)
c     KEY    Integer, key to selected local integration rule.
c            KEY = 1 gives the user a degree 3 integration rule.
c            KEY = 2 gives the user a degree 5 integration rule.
c            KEY = 3 gives the user a degree 7 integration rule.
c            KEY = 4 gives the user a degree 9 integration rule.
c     RULCLS Integer, number of FUNSUB calls needed for each subregion.
c     MAXSUB Integer; computations proceed until there are at most
c            MAXSUB subregions in the data structure.
c     SBRGNS Integer.
c            If RESTAR = 0, then SBRGNS must specify the number
c            of subregions stored in a previous call to SMPSAD.
c     VERTCS Real array of dimension (ND,0:ND,*).
c            Simplex vertices for each subregion; for subregion K vertex
c            J must have components VERTEX(I,J,K), I = 1, 2, ..., ND.
c     VALUES Real array of dimension (NF,*), for estimated values of the
c             integrals over the subregions.
c     ERRORS Real array of dimension (NF,*).
c            Used to store the corresponding estimated errors.
c            Used to store the half widths of the stored subregions.
c     GREATS Real array of dimension (*).
c            Used to store the greatest estimated errors in subregions.
c     PONTRS Real array of dimension (*), for the pointers from the
c             subregion heap to the actual subregions.
c     WORK   Real array, used in SMPVOL, SMPDFS, and SMPRUL.
c
c   ON RETURN
c
c     SBRGNS Integer, number of stored subregions.
c     VALUE  Real array of dimension NF.
c            Approximations to all components of the integral.
c     ERROR  Real array of dimension NF.
c            Estimates of absolute accuracies.
c     FUNCLS Integer, number of new FUNSUB calls used by SMPSAD.
c     INFORM Integer.
c            INFORM = 0 for normal exit, when ERROR(K) <=  EPSABS or
c              ERROR(K) <=  ABS(VALUE(K))*EPSREL, 1 <= K <= NF,
c              with MAXSUB or fewer subregions processed.
c            INFORM = 1 if MAXSUB was too small for SMPSAD
c              to obtain the required accuracy. In this case SMPSAD
c              returns values of VALUE with estimated absolute
c              accuracies ERROR.
c
c***REFERENCES
c***ROUTINES CALLED SMPSTR, SMPVOL, SMPRUL
c***END PROLOGUE SMPSAD
c
c   Global variables.
c
      external funsub
      integer nd, nf, rulcls, mincls, maxcls, maxsub, key, restar
      integer funcls, sbrgns, inform
      double precision epsabs, epsrel, value(nf), error(nf)
      double precision values(nf,*), errors(nf,*), vertcs(nd,0:nd,*)
      double precision volums(*), greats(*), pontrs(*), work(*)
c
c   Local variables.
c
c
c   MXNWSB is the maxiumum number of new subregions per subdivision.
c
      integer i, index, j, top, mxnwsb, newsbs, dfcost, rgncls
      parameter ( mxnwsb = 4 )
      double precision smpvol, tune
      parameter( tune = 1 )
c
c***FIRST PROCESSING STATEMENT SMPSAD
c
c
c  Initialize for rule parameters.
c
      funcls = 0
      dfcost = 1 + 2*nd*( nd + 1 )
c
c  If RESTAR = 0, initialize for first call.
c
      if ( restar .eq. 0 ) then
c
c  Initialize FUNCLS, and VALUE and ERROR arrays.
c
         do j = 1, nf
            value(j) = 0
            error(j) = 0
         end do
         do index = 1, sbrgns
c
c  Call SMPVOL to compute the simplex volume(s).
c
            volums(index) = smpvol( nd, vertcs(1,0,index), work )
c
c  Apply basic rule over each simplex.
c
            call smprul( tune, nd, vertcs(1,0,index), volums(index),
     &           nf, funsub, key, values(1,index), errors(1,index),
     &           greats(index), work, work(2*nd+2) )
c
c  Add new contributions to VALUE and ERROR.
c  Store results in heap.
c
            do j = 1, nf
               value(j) = value(j) + values(j,index)
               error(j) = error(j) + errors(j,index)
            end do
            call smpstr( index, index, pontrs, greats )
            funcls = funcls + rulcls
         end do
      end if
      inform = max( 0, min( mincls - funcls, 1 ) )
      do j = 1, nf
         if( error(j) .gt. max(epsabs,epsrel*abs(value(j))) ) inform = 1
      end do
c
c  End initialisation.
c
      do while ( inform .gt. 0 .and. sbrgns + mxnwsb - 1 .le. maxsub
     &           .and. funcls + dfcost + mxnwsb*rulcls .le. maxcls )
c
c  Begin loop while error is too large, and FUNCLS and SBRGNS
c  are not too large.
c
c  Adjust VALUE and ERROR.
c
         top = pontrs(1)
         do j = 1, nf
            value(j) = value(j) - values(j,top)
            error(j) = error(j) - errors(j,top)
         end do
c
c  Determine NEWSBS new subregions.
c
         call smpdfs( nd, nf, funsub, top, sbrgns, vertcs,
     &                volums, work, work(nd+1), work(2*nd+1),
     &                work(3*nd+1), work(3*nd+1+5*nf), newsbs )
c
c  Apply basic rule, store results in heap and
c  add new contributions to VALUE and ERROR.
c
         index = top
         do i = 1, newsbs
            call smprul( tune, nd, vertcs(1,0,index), volums(index), nf,
     &                   funsub, key, values(1,index), errors(1,index),
     &                   greats(index), work, work(2*nd+2) )
            call smpstr( index, sbrgns+i-1, pontrs, greats )
            do j = 1, nf
               value(j) = value(j) + values(j,index)
               error(j) = error(j) + errors(j,index)
            end do
            index = sbrgns + i
         end do
         funcls = funcls + dfcost + newsbs*rulcls
         sbrgns = sbrgns + newsbs - 1
c
c  Check for error termination.
c
         inform = max( 0, min( mincls - funcls, 1 ) )
         do j = 1, nf
            if( error(j) .gt. max(epsabs,epsrel*abs(value(j))) )
     &           inform = 1
         end do
      end do
c
c  Compute more accurate values of VALUE and ERROR.
c
      do i = 1, nf
         value(i) = 0
         error(i) = 0
         do j = 1, sbrgns
            value(i) = value(i) + values(i,j)
            error(i) = error(i) + errors(i,j)
         end do
      end do

      return
      end
      function smpvol ( nd, vertex, work )

c*********************************************************************72
c
cc SMPVOL computes the scaled volume of a simplex.
c
c  Discussion:
c
c    This routine computes the volume of an ND-simplex scaled by 
c    using Gauss elimination to compute a determinant.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Alan Genz
c    Department of Mathematics
c    Washington State University
c    Pullman, WA 99164-3113, USA
c
c  Parameters:
c
c    Input, integer ND, the number of variables.
c
c    Input, double precision VERTEX(ND,0:ND), the simplex vertices;
c    The simplex vertex J must have components VERTEX(I,J), I = 1, 2, ..., ND.
c
c    Input, double precision WORK(ND,ND).
c
c    Output, double precision SMPVOL, the volume.
c
      implicit none

      integer nd

      integer i
      integer j
      integer k
      double precision mult
      integer pivpos
      double precision smpvol
      double precision vertex(nd,0:nd)
      double precision vol
      double precision work(nd,nd)
      double precision wtemp
c
c  Copy vertex differences to WORK array.
c
      do j = 1, nd
        do i = 1, nd
          work(i,j) = vertex(i,j) - vertex(i,0)
        end do
      end do
c
c  Use Gauss elimination with partial pivoting.
c
      vol = 1

      do k = 1, nd

        pivpos = k
        do j = k+1, nd
          if ( abs ( work(k,pivpos) ) .lt. abs ( work(k,j) ) ) then
            pivpos = j
          end if
        end do

        do i = k, nd
          wtemp = work(i,k)
          work(i,k) = work(i,pivpos)
          work(i,pivpos) = wtemp
        end do

        vol = vol * work(k,k) / k

        do j = k+1, nd
          mult = work(k,j) / work(k,k)
          do i = k+1, nd
            work(i,j) = work(i,j) - mult * work(i,k)
          end do
        end do

      end do

      smpvol = abs ( vol )

      return
      end
      subroutine smprul ( tune, nd, vertex, volume, nf, intgnd,
     &  inkey, basval, rgnerr, great, gt, rule )

c*********************************************************************72
c
cc SMPRUL computes the integration rules.
c
c***BEGIN PROLOGUE SMPRUL
c***KEYWORDS basic numerical integration rule
c***PURPOSE  To compute basic integration rule values.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Alan Genz
c    Department of Mathematics
c    Washington State University
c    Pullman, WA 99164-3113, USA
c
c***DESCRIPTION SMPRUL computes basic integration rule values for a
c            vector of integrands over a hyper-rectangular region.
c            These are estimates for the integrals. SMPRUL also computes
c            estimates for the errors.
c
c   ON ENTRY
c
c     TUNE   Real, tuning parameter, with 0 <= TUNE <= 1, with
c            TUNE = 1 for the most conservative error estimates.
c            If TUNE < 0, only the rule parameters are computed.
c     ND    Integer, number of variables.
c     VERTEX Real array of dimension (ND,0:ND).
c            The simplex vertices; vertex J must have components
c            VERTEX(I,J), I = 1, 2, ..., ND.
c     NF     Integer, number of components for the vector integrand.
c     INTGND Subroutine for computing components of the integrand at Z.
c            It must have parameters (ND,X,NF,FUNVLS)
c            Input parameters:
c              ND    Integer that gives the dimension.
c              X      Real array of dimension ND that contains the
c                     evaluation point.
c              NF     Integer that gives the number of components of I.
c            Output parameter:
c              FUNVLS Real array of dimension NF that contains the
c                     components of the integrand.
c     INKEY  Integer rule parameter.
c            If INKEY .GT. 0 and INKEY .LT. 5 then a rule of degree
c            2*INKEY + 1; otherwise default degree 7 rule is used.
c     GT     Real work array of length 2*ND+1.
c     RULE   Real work array of dimension (NF,7).
c
c   ON RETURN
c
c     BASVAL Real array of length NF, values for the basic rule for
c            each component of the integrand.
c     RGNERR Real array of length NF, error estimates for BASVAL.
c     GREAT  Real, maximum component of RGNERR.
c
c
c***ROUTINES CALLED: SMPRMS, SYMRUL
c
c***END PROLOGUE SMPRUL
c
c   Global variables.
c
      external intgnd
      integer nf, nd, inkey
      double precision vertex(nd,0:nd), basval(nf), rgnerr(nf)
      double precision volume, tune, great
c
c   Local variables.
c
c   WTS    Integer number of weights in the integration rules.
c   W      Real array of dimension (WTS,RLS).
c          The weights for the basic and null rules.
c          W(1,1),...,W(WTS,1) are weights for the basic rule.
c          W(1,I),...,W(WTS,I), for I > 1 are null rule weights.
c   G      Real array of dimension (0:4, WTS).
c          The fully symmetric sum generators for the rules.
c
      integer key, numnul, rls, wts, mxw, mxrls, mxg
      parameter( mxw = 21, mxrls = 7, mxg = 4  )
      double precision w( mxw, mxrls ), g( 0:mxg, mxw ), wtsum
      double precision gt( 0:2*nd ), rule( nf, mxrls )
      double precision normcf, normnl, normcp, alpha(mxrls)
      double precision ratio, errcof, ratmin, small, smprod
      parameter( ratmin = 1d-1, small = 1d-12 )
      integer i, j, k, oldkey, oldn, pts(mxw)
      save oldkey, oldn, key, pts, w, g, rls, wts
      data oldkey, oldn/ -1, 0 /
c
c***FIRST PROCESSING STATEMENT SMPRUL
c
      if ( oldkey .ne. inkey .or. oldn .ne. nd ) then
         oldn = nd
         oldkey = inkey
         if ( inkey .gt. 0 .and. inkey .lt. 5 ) then
            key = inkey
         else
            key = 3
         end if
c
c  Compute WTS, RLS, weights, generators, ERRCOF and PTS.
c
         call smprms( nd, key, mxw, w, mxg, g, wts, rls, pts )
c
c  Orthogonalize and normalize null rules.
c
         normcf = smprod( wts, pts, w(1,1), w(1,1) )
         do k = 2, rls
            do j = 2, k-1
               alpha(j) = -smprod( wts, pts, w(1,j), w(1,k) )
            end do
            do i = 1, wts
               wtsum = 0
               do j = 2, k-1
                  wtsum = wtsum + w(i,j)*alpha(j)
               end do
               w(i,k) = w(i,k) + wtsum/normcf
            end do
            normnl = smprod( wts, pts, w(1,k), w(1,k) )
            do i = 1, wts
               w(i,k) = w(i,k)*sqrt( normcf/normnl )
            end do
         end do
      end if
      if ( tune .ge. 0 ) then
c
c  Compute the rule values.
c
         do i = 1, nf
            do j = 1, rls
               rule(i,j) = 0
            end do
         end do
         do k = 1, wts
            if ( pts(k) .gt. 0 ) then
               do i = 0, min(nd,mxg-1)
                  gt(i) = g(i,k)
               end do
               if ( nd .ge. mxg ) call smpcpy( mxg, nd, gt, g(mxg,k) )
               call smpsms( nd, vertex, nf, intgnd, gt, basval,
     &                                              gt(nd+1), rgnerr )
               do j = 1, rls
                  do i = 1, nf
                     rule(i,j) = rule(i,j) + w(k,j)*basval(i)
                  end do
               end do
            end if
         end do
c
c  Scale integral values and compute the error estimates.
c
         errcof = ( 8*tune + ( 1 - tune ) )
         great = 0
         do i = 1, nf
            basval(i) = rule(i,1)
            normcf = abs( basval(i) )
            rgnerr(i) = 0
            ratio = ratmin
            do k = rls, 3, -2
               normnl = max( abs( rule(i,k) ), abs( rule(i,k-1) ) )
               if ( normnl .gt. small*normcf .and. k .lt. rls )
     &              ratio = max( normnl/normcp, ratio )
               rgnerr(i) = max( normnl, rgnerr(i) )
               normcp = normnl
            end do
            if( ratio .ge. 1 ) then
               rgnerr(i) = tune*rgnerr(i) + ( 1 - tune )*normcp
            else if ( key .gt. 1 ) then
               rgnerr(i) = ratio*normcp
            end if
            rgnerr(i) = volume*max( errcof*rgnerr(i), small*normcf )
            basval(i) = volume*basval(i)
            great = max( great, rgnerr(i) )
         end do
      end if

      return
      end
      function smprod ( n, w, x, y )

c*********************************************************************72
c
cc SMPROD computes a weighted dot product.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Alan Genz
c    Department of Mathematics
c    Washington State University
c    Pullman, WA 99164-3113, USA
c
      integer n, i, w(*)
      double precision smprod
      double precision x(*), y(*), sum

      sum = 0.0D+00
      do i = 1, n
         sum = sum + w(i)*x(i)*y(i)
      end do
      smprod = sum

      return
      end
      subroutine smprms ( nd, key, mxw, w, mxg, g, wts, rls, pts )

c*********************************************************************72
c
cc SMPRMS initializes a basic rule and some null rules.
c
c***BEGIN PROLOGUE SMPRMS
c***KEYWORDS basic integration rule, degree 2*KEY+1
c***PURPOSE  To initialize a degree 2*KEY+1 basic rule and null rules.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Alan Genz
c    Department of Mathematics
c    Washington State University
c    Pullman, WA 99164-3113, USA
c
c***DESCRIPTION  SMPRMS initializes a degree 2*KEY+1 rule, and
c                and max(2*KEY,2) lower degree null rules.
c
c   ON ENTRY
c
c   ND    Integer, number of variables.
c   KEY    Integer, < 5 and >= 0, rule parameter.
c          If KEY > 0 a degree 2*KEY+1 rule is initialized.
c          If KEY = 0 a degree 7 rule is initialized.
c
c   ON RETURN
c   RLS    Integer, total number of rules.
c   WTS    Integer, total number of weights in each of the rules.
c   W      Real array of dimension (MXW,*).
c          The weights for the basic and null rules.
c          W(1,1),...,W(WTS,1) are weights for the basic rule.
c          W(I,1),...,W(WTS,I) for I .GT. 1 are null rule weights.
c   G      Real array of dimension (0:MXG,MXW).
c          The fully symmetric sum generators for the rules.
c          G(0,J), ..., G(MXG,J) are the generators for the
c          points associated with the Jth weights.
c   PTS    Integer array of length (MXW). PTS(J) is number of integrand
c          values needed for generator J.
c
c***REFERENCES
c
c    Axel Grundmann, Michael Moeller,
c    Invariant Integration Formulas for the N-Simplex
c    by Combinatorial Methods,
c    SIAM Journal on Numerical Analysis,
c    Volume 15, Number 2, April 1978, pages 282-290.
c
c    Arthur Stroud,
c    A Fifth Degree Integration Formula for the n-Simplex,
c    SIAM Journal on Numerical Analysis,
c    Volume 6, Number 1, March 1969, pages 90-98.
c
c    Ivan Mysovskikh,
c    On a cubature formula for the simplex,
c    Voprosy Vycislitelnoj i Prikl. Mat.,
c    Volume 51, 1978, pages 74-90.
c
c***END PROLOGUE SMPRMS
c
c   Global variables
c
      integer nd, key, wts, mxw, rls, mxg
      integer pts(mxw)
      double precision w(mxw,*), g(0:mxg,*)
c
c   Local Variables
c
      double precision one, ffteen
      parameter( one = 1, ffteen = 15 )
      double precision dr, dr2, dr4, dr6, dr8
      double precision r1, s1, r2, s2, u1, v1, u2, v2, l1, l2, d1, d2
      double precision a1, a2, a3, p0, p1, p2, p3, u5, u6, u7, sg
      double precision r, a, p, q, th, tp
      integer iw, gms, i, j
c
c***FIRST PROCESSING STATEMENT SMPRMS
c
c
c     Initialize RLS and GMS.
c
      if ( key .eq. 1 ) then
         rls = 3
         gms = 2
         wts = 3
      else if ( key .eq. 2 ) then
         rls = 5
         gms = 4
         wts = 6
      else if ( key .eq. 3 .or. key .eq. 0 ) then
         rls = 7
         gms = 7
         wts = 11
      else if ( key .eq. 4 ) then
         rls = 7
         if ( nd .eq. 2 ) then
            gms = 11
            wts = 20
         else
            gms = 12
            wts = 21
         end if
      end if
c
c  Initialize generators, weights and PTS.
c
      do i = 1, wts
         do j = 1, rls
            w(i,j) = 0
         end do
         pts(i) = 0
      end do
c
c  Compute generator, PTS and weight values for all rules.
c
      dr = nd
      dr2 =     ( dr + 1 )*( dr + 2 )
      dr4 = dr2*( dr + 3 )*( dr + 4 )
      dr6 = dr4*( dr + 5 )*( dr + 6 )
      dr8 = dr6*( dr + 7 )*( dr + 8 )
      call smpcpy( 0, mxg, g(0,1), 1/( dr + 1 ) )
      pts(1) = 1
      r1 = ( dr + 4 - sqrt(ffteen) )/( dr*dr + 8*dr + 1 )
      s1 = 1 - dr*r1
      l1 = s1 - r1
      g(0   ,gms+1) = s1
      call smpcpy( 1, mxg, g(0,gms+1), r1 )
      do i = 1, mxg
         g(i,gms+1) = r1
      end do
      pts(gms+1) = dr + 1
      iw = rls
      if ( key .lt. 4 )  then
c
c  Compute weights for special degree 1 rule.
c
         w(1,iw) = 1
         iw = iw - 1
         w(gms+1,iw) = 1/( dr + 1 )
         iw = iw - 1
      end if
c
c  Compute weights, generators and PTS for degree 3 rule.
c
      g(0,2) = 3/( dr + 3 )
      call smpcpy( 1, mxg, g(0,2), 1/( dr + 3 ) )
      pts(2) = dr + 1
        w(2,iw) = ( dr + 3 )**3/( 4*dr2*( dr + 3 ) )
      if ( key .gt. 1 ) then
         iw = iw - 1
c
c  Compute weights, generators and PTS for degree 3 and 5 rules.
c
         if ( nd .eq. 2 ) then
c
c  Special degree 3 rule.
c
            l2 = .62054648267200632589046034361711d0
            l1 = -sqrt( one/2 - l2**2 )
            r1 = ( 1 - l1 )/3
            s1 = 1 - 2*r1
            g(0,gms+1) = s1
            call smpcpy( 1, mxg, g(0,gms+1), r1 )
            pts(gms+1) = 3
              w(gms+1,iw) = one/6
            r2 = ( 1 - l2 )/3
            s2 = 1 - 2*r2
            g(0,gms+2) = s2
            call smpcpy( 1, mxg, g(0,gms+2), r2 )
            pts(gms+2) = 3
              w(gms+2,iw) = one/6
         else
c
c  Degree 3 rule using Stroud points.
c
            r2 = ( dr + 4 + sqrt(ffteen) )/( dr*dr + 8*dr + 1 )
            s2 = 1 - dr*r2
            l2 = s2 - r2
            g(0,gms+2) = s2
            call smpcpy( 1, mxg, g(0,gms+2), r2 )
            pts(gms+2) = dr + 1
              w(gms+2,iw) = ( 2/(dr+3) - l1 )/( dr2*(l2-l1)*l2**2 )
              w(gms+1,iw) = ( 2/(dr+3) - l2 )/( dr2*(l1-l2)*l1**2 )
         end if
         iw = iw - 1
c
c        Grundmann-Moller degree 5 rule.
c
         g(0,3) = 5/( dr + 5 )
         call smpcpy( 1, mxg, g(0,3), 1/( dr + 5 ) )
         pts(3) = dr + 1
         g(0,4) = 3/( dr + 5 )
         g(1,4) = 3/( dr + 5 )
         call smpcpy( 2, mxg, g(0,4), 1/( dr + 5 ) )
         pts(4) = ( ( dr + 1 )*dr )/2
           w(2,iw) = -( dr + 3 )**5/( 16*dr4 )
           w(3,iw) =  ( dr + 5 )**5/( 16*dr4*( dr + 5 ) )
           w(4,iw) =  ( dr + 5 )**5/( 16*dr4*( dr + 5 ) )
      end if
      if ( key .gt. 2 )  then
         iw = iw - 1
c
c  Compute weights, generators and PTS for degree 5 and 7 rules.
c
c
c  Stroud degree 5 rule.
c
         u1 = ( dr + 7 + 2*sqrt(ffteen) )/( dr*dr + 14*dr - 11 )
         v1 = ( 1 - ( dr - 1 )*u1 )/2
         d1 = v1 - u1
         g(0,gms+3) = v1
         g(1,gms+3) = v1
         call smpcpy( 2, mxg, g(0,gms+3), u1 )
         pts(gms+3) = ( ( dr + 1 )*dr )/2
         u2 = ( dr + 7 - 2*sqrt(ffteen) )/( dr*dr + 14*dr - 11 )
         v2 = ( 1 - ( dr - 1 )*u2 )/2
         d2 = v2 - u2
         g(0,gms+4) = v2
         g(1,gms+4) = v2
         call smpcpy( 2, mxg, g(0,gms+4), u2 )
         pts(gms+4) = ( ( dr + 1 )*dr )/2
         if ( nd .eq. 2 ) then
            w(gms+3,iw) = ( 155 - sqrt(ffteen) )/1200
            w(gms+4,iw) = ( 155 + sqrt(ffteen) )/1200
            w(1,    iw) = 1 - 3*( w(gms+3,iw) + w(gms+4,iw) )
         else if ( nd .eq. 3 ) then
            w(gms+1,iw) = ( 2665 + 14*sqrt(ffteen) )/37800
            w(gms+2,iw) = ( 2665 - 14*sqrt(ffteen) )/37800
            w(gms+3,iw) = 2*ffteen/567
            pts(gms+4) = 0
         else
            w(gms+1,iw) = ( 2*(27-dr)/(dr+5)-l2*(13-dr) )
     &                       /( l1**4*(l1-l2)*dr4 )
            w(gms+2,iw) = ( 2*(27-dr)/(dr+5)-l1*(13-dr) )
     &                       /( l2**4*(l2-l1)*dr4 )
            w(gms+3,iw)=( 2/( dr + 5 ) - d2 )/( dr4*( d1 - d2 )*d1**4 )
            w(gms+4,iw)=( 2/( dr + 5 ) - d1 )/( dr4*( d2 - d1 )*d2**4 )
         end if
         iw = iw - 1
c
c  Grundmann-Moller degree 7 rule.
c
         g(0,5) = 7/( dr + 7 )
         call smpcpy( 1, mxg, g(0,5), 1/( dr + 7 ) )
         pts(5) = dr + 1
         g(0,6) = 5/( dr + 7 )
         g(1,6) = 3/( dr + 7 )
         call smpcpy( 2, mxg, g(0,6), 1/( dr + 7 ) )
         pts(6) = ( dr + 1 )*dr
         g(0,7) = 3/( dr + 7 )
         g(1,7) = 3/( dr + 7 )
         g(2,7) = 3/( dr + 7 )
         call smpcpy( 3, mxg, g(0,7), 1/( dr + 7 ) )
         pts(7) = ( ( dr + 1 )*dr*( dr - 1 ) )/6
         w(2,iw) =  ( dr + 3 )**7/( 2*64*dr4*( dr + 5 ) )
         w(3,iw) = -( dr + 5 )**7/(   64*dr6 )
         w(4,iw) = -( dr + 5 )**7/(   64*dr6 )
         w(5,iw) =  ( dr + 7 )**7/(   64*dr6*( dr + 7 ) )
         w(6,iw) =  ( dr + 7 )**7/(   64*dr6*( dr + 7 ) )
         w(7,iw) =  ( dr + 7 )**7/(   64*dr6*( dr + 7 ) )
      end if
      if ( key .eq. 4 )  then
         iw = iw - 1
c
c  Compute weights, generators and PTS for degree 7, 9 rules.
c
c  Mysovskikh degree 7 rule.
c
         sg = 1/( 23328*dr6 )
         u5 = -6**3*sg*( 52212 - dr*( 6353 + dr*( 1934 - dr*27 ) ) )
         u6 =  6**4*sg*(  7884 - dr*( 1541 - dr*9 ) )
         u7 = -6**5*sg*(  8292 - dr*( 1139 - dr*3 ) )/( dr + 7 )
         p0 = -144*( 142528 + dr*( 23073 - dr*115 ) )
         p1 = -12*( 6690556 + dr*( 2641189 + dr*( 245378 - dr*1495 ) ) )
         p2 = -16*(6503401 + dr*(4020794+dr*(787281+dr*(47323-dr*385))))
         p3 = -( 6386660 + dr*(4411997+dr*(951821+dr*(61659-dr*665))) )
     &        *( dr + 7 )
         a = p2/( 3*p3 )
         p = a*( p1/p2 - a )
         q = a*( 2*a*a - p1/p3 ) + p0/p3
         r = sqrt( -p**3 )
         th = acos( -q/( 2*r ) )/3
         r = 2*r**( one/3 )
         tp = 2*acos(-one)/3
         a1 = -a + r*cos( th )
         a2 = -a + r*cos( th + tp + tp )
         a3 = -a + r*cos( th + tp )
         g(0,gms+5) = ( 1 - dr*a1 )/( dr + 1 )
         call smpcpy( 1, mxg, g(0,gms+5), ( 1 + a1 )/( dr + 1 ) )
         pts(gms+5) = dr + 1
         g(0,gms+6) = ( 1 - dr*a2 )/( dr + 1 )
         call smpcpy( 1, mxg, g(0,gms+6), ( 1 + a2 )/( dr + 1 ) )
         pts(gms+6) = dr + 1
         g(0,gms+7) = ( 1 - dr*a3 )/( dr + 1 )
         call smpcpy( 1, mxg, g(0,gms+7), ( 1 + a3 )/( dr + 1 ) )
         pts(gms+7) = dr + 1
           w(gms+5,iw) = ( u7-(a2+a3)*u6+a2*a3*u5 )
     &                  /( a1**2-(a2+a3)*a1+a2*a3 )/a1**5
           w(gms+6,iw) = ( u7-(a1+a3)*u6+a1*a3*u5 )
     &                  /( a2**2-(a1+a3)*a2+a1*a3 )/a2**5
           w(gms+7,iw) = ( u7-(a2+a1)*u6+a2*a1*u5 )
     &                  /( a3**2-(a2+a1)*a3+a2*a1 )/a3**5
         g(0,gms+8) = 4/( dr + 7 )
         g(1,gms+8) = 4/( dr + 7 )
         call smpcpy( 2, mxg, g(0,gms+8), 1/( dr + 7 ) )
         pts(gms+8) = ( ( dr + 1 )*dr )/2
           w(gms+8,iw) = 10*(dr+7)**6/( 729*dr6 )
         g(0,gms+9) = 11/( dr + 7 )/2
         g(1,gms+9) =  5/( dr + 7 )/2
         call smpcpy( 2, mxg, g(0,gms+9), 1/( dr + 7 ) )
         pts(gms+9) = ( ( dr + 1 )*dr )
           w(gms+9,iw) = 64*(dr+7)**6/( 6561*dr6 )
           w(    4,iw) = w(4,iw+1)
           w(    7,iw) = w(7,iw+1)
         iw = iw - 1
c
c  Grundmann-Moller degree 9 rule.
c
         g(0,8) = 9/( dr + 9 )
         call smpcpy( 1, mxg, g(0, 8), 1/( dr + 9 ) )
         pts(8) = dr + 1
         g(0,9) = 7/( dr + 9 )
         g(1,9) = 3/( dr + 9 )
         call smpcpy( 2, mxg, g(0, 9), 1/( dr + 9 ) )
         pts(9) = ( dr + 1 )*dr
         g(0,10) = 5/( dr + 9 )
         g(1,10) = 5/( dr + 9 )
         call smpcpy( 2, mxg, g(0,10), 1/( dr + 9 ) )
         pts(10) = ( ( dr + 1 )*dr )/2
         g(0,11) = 5/( dr + 9 )
         g(1,11) = 3/( dr + 9 )
         g(2,11) = 3/( dr + 9 )
         call smpcpy( 3, mxg, g(0,11), 1/( dr + 9 ) )
         pts(11) = ( ( dr + 1 )*dr*( dr - 1 ) )/2
           w(2 ,iw) = -( dr + 3 )**9/( 6*256*dr6 )
           w(3 ,iw) =  ( dr + 5 )**9/( 2*256*dr6*( dr + 7 ) )
           w(4 ,iw) =  ( dr + 5 )**9/( 2*256*dr6*( dr + 7 ) )
           w(5 ,iw) = -( dr + 7 )**9/(   256*dr8 )
           w(6 ,iw) = -( dr + 7 )**9/(   256*dr8 )
           w(7 ,iw) = -( dr + 7 )**9/(   256*dr8 )
           w(8 ,iw) =  ( dr + 9 )**9/(   256*dr8*( dr + 9 ) )
           w(9 ,iw) =  ( dr + 9 )**9/(   256*dr8*( dr + 9 ) )
           w(10,iw) =  ( dr + 9 )**9/(   256*dr8*( dr + 9 ) )
           w(11,iw) =  ( dr + 9 )**9/(   256*dr8*( dr + 9 ) )
         if ( nd .gt. 2 ) then
            g(0,12) = 3/( dr + 9 )
            g(1,12) = 3/( dr + 9 )
            g(2,12) = 3/( dr + 9 )
            g(3,12) = 3/( dr + 9 )
            call smpcpy( 4, mxg, g(0,12), 1/( dr + 9 ) )
            pts(12) = ( ( dr + 1 )*dr*( dr - 1 )*( dr - 2 ) )/24
              w(12,iw) = w(8,iw)
         end if
      end if
c
c  Compute constant weight values.
c
      do j = 1, rls
         w(1,j) = 1
         do i = 2, wts
            w(1,j) = w(1,j) - pts(i)*w(i,j)
         end do
      end do
c
c  Compute final weight values; null rule weights are computed as
c  differences between weights from highest degree and lower degree rules.
c
      do j = 2, rls
         do i = 1, wts
            w(i,j) = w(i,j) - w(i,1)
         end do
      end do

      return
      end
      subroutine smpcpy ( start, end, param, value )

c*********************************************************************72
c
cc SMPCPY sets a vector to a constant.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Alan Genz
c    Department of Mathematics
c    Washington State University
c    Pullman, WA 99164-3113, USA
c
      double precision value, param(0:*)
      integer start, end, i
      do i = start, end
         param(i) = value
      end do

      return
      end
      subroutine smpsms ( n, vertex, nf, f, g, symsms, x, funvls )

c*********************************************************************72
c
cc SMPSMS computes a symmetric sum over a simplex.
c
c***BEGIN PROLOGUE SMPSMS
c***KEYWORDS fully symmetric sum
c***PURPOSE  To compute fully symmetric basic rule sums
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Alan Genz
c    Department of Mathematics
c    Washington State University
c    Pullman, WA 99164-3113, USA
c
c***DESCRIPTION SMPSMS computes a fully symmetric sum for a vector
c            of integrand values over a simplex. The sum is taken over
c            all permutations of the generators for the sum.
c
c   ON ENTRY
c
c   N       Integer, number of variables.
c   VERTEX  Real array of dimension (N,0:N)
c           The vertices of the simplex, one vertex per column.
c   NF      Integer, number of components for the vector integrand.
c   F       Subroutine for computing components of the integrand at X.
c            It must have parameters ( N, X, NF, FUNVLS );
c            Input parameters:
c              N      Integer dimension of integral.
c              X      Real array of length N, the evaluation point.
c              NF     Integer number of components of the integrand.
c            Output parameter:
c             FUNVLS  Real array of length NF, the integrand values at X.
c   G       Real Array of dimension (0:N).
c           The generators for the fully symmetric sum.
c
c   ON RETURN
c
c   SYMSMS  Real array of length NF, the values for the fully symmetric
c            sums for each component of the integrand.
c
c***ROUTINES CALLED: Integrand
c
c***END PROLOGUE SMPSMS
c
c   Global variables.
c
      integer n, nf
      double precision vertex(n,0:n),g(0:n), symsms(nf),funvls(nf), x(n)
c
c   Local variables.
c
      integer ix, lx, i, j, k, l
      double precision gl, gi
c
c***FIRST PROCESSING STATEMENT SymSum
c
      do i = 1, nf
         symsms(i) = 0
      end do
c
c  Sort generators if necessary
c
      k = 0
      do i = 1, n
         if ( g(i) .gt. g(i-1) ) k = 1
      end do
      if ( k .gt. 0 ) then
         do i = 1, n
            k = i - 1
            do j = i, n
               if ( g(j) .gt. g(k) ) k = j
            end do
            if ( k .ge. i ) then
               gi = g(i-1)
               g(i-1) = g(k)
               g(k) = gi
            end if
         end do
      end if
c
c  Compute integrand value for permutations of G
c
 10   continue

      do i = 1, n
         x(i) = vertex(i,0)*g(0)
         do j = 1, n
            x(i) = x(i) + vertex(i,j)*g(j)
         end do
      end do
      call f( n, x, nf, funvls )
      do j = 1, nf
         symsms(j) = symsms(j) + funvls(j)
      end do
c
c  Find next distinct permuation of G and loop back for value.
c  Permutations are generated in reverse lexicographic order.
c
      do i = 1, n
         if ( g(i-1) .gt. g(i) ) then
            gi = g(i)
            ix = i - 1
            do l = 0, i/2-1
               gl = g(l)
               g(l) = g(i-l-1)
               g(i-l-1) = gl
               if (  gl .le. gi ) ix = ix - 1
               if ( g(l) .gt. gi ) lx = l
            end do
            if ( g(ix) .le. gi ) ix = lx
            g(i) = g(ix)
            g(ix) = gi
            go to 10
         end if
      end do

      return
      end
      subroutine smpstr ( pointr, sbrgns, pontrs, rgners )

c*********************************************************************72
c
cc SMPSTR maintains a heap to keep track of subregions.
c
c***BEGIN PROLOGUE SMPSTR
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Alan Genz
c    Department of Mathematics
c    Washington State University
c    Pullman, WA 99164-3113, USA
c
c***PURPOSE SMPSTR maintains a heap for subregions.
c***DESCRIPTION SMPSTR maintains a heap for subregions.
c            The subregions are ordered according to the size of the
c            greatest error estimates of each subregion (RGNERS).
c
c   PARAMETERS
c
c     POINTR Integer.
c            The index for the subregion to be inserted in the heap.
c     SBRGNS Integer.
c            Number of subregions in the heap.
c     PONTRS Real array of dimension SBRGNS.
c            Used to store the indices for the greatest estimated errors
c            for each subregion.
c     RGNERS Real array of dimension SBRGNS.
c            Used to store the greatest estimated errors for each
c            subregion.
c
c***ROUTINES CALLED NONE
c***END PROLOGUE SMPSTR
c
c   Global variables.
c
      integer pointr, sbrgns
      double precision pontrs(*), rgners(*)
c
c   Local variables.
c
c   RGNERR Intermediate storage for the greatest error of a subregion.
c   SUBRGN Position of child/parent subregion in the heap.
c   SUBTMP Position of parent/child subregion in the heap.
c
      integer subrgn, subtmp, pt, ptp
      double precision rgnerr
c
c***FIRST PROCESSING STATEMENT SMPSTR
c
      rgnerr = rgners(pointr)
      if ( pointr .eq. pontrs(1) ) then
c
c  Move the new subregion inserted at the top of the heap
c  to its correct position in the heap.
c
         subrgn = 1

 10      continue

         subtmp = 2*subrgn

         if ( subtmp .le. sbrgns ) then
            if ( subtmp .ne. sbrgns ) then
c
c  Find maximum of left and right child.
c
               pt = pontrs(subtmp)
               ptp = pontrs(subtmp+1)
               if ( rgners(pt) .lt. rgners(ptp) ) subtmp = subtmp + 1
            end if
c
c  Compare maximum child with parent.
c  If parent is maximum, then done.
c
            pt = pontrs(subtmp)
            if ( rgnerr .lt. rgners(pt) ) then
c
c  Move the pointer at position subtmp up the heap.
c
               pontrs(subrgn) = pt
               subrgn = subtmp
               go to 10
            end if
         end if
      else
c
c  Insert new subregion in the heap.
c
         subrgn = sbrgns

 20      continue

         subtmp = subrgn/2
         if ( subtmp .ge. 1 ) then
c
c  Compare child with parent. If parent is maximum, then done.
c
            pt = pontrs(subtmp)
            if ( rgnerr .gt. rgners(pt) ) then
c
c  Move the pointer at position subtmp down the heap.
c
               pontrs(subrgn) = pt
               subrgn = subtmp
               go to 20
            end if
         end if
      end if
      pontrs(subrgn) = pointr

      return
      end
      subroutine smpdfs ( nd, nf, funsub, top, sbrgns, vertcs, volums,
     &  x, h, center, work, frthdf, newsbs )

c*********************************************************************72
c
cc SMPDFS determines how to subdivide the subregion.
c
c***BEGIN PROLOGUE SMPDFS
c***PURPOSE  To compute new subregions
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Alan Genz
c    Department of Mathematics
c    Washington State University
c    Pullman, WA 99164-3113, USA
c
c***DESCRIPTION SMPDFS computes fourth differences along each edge
c            direction. It uses these differences to determine a
c            subdivision of the orginal subregion into either three or
c            four new subregions.
c
c   ON ENTRY
c
c   ND   Integer, number of variables.
c   NF   Integer, number of components for the vector integrand.
c   FUNSUB Externally declared subroutine.
c          For computing the components of the integrand at a point X.
c          It must have parameters (ND, X, NF, FUNVLS).
c           Input Parameters:
c            X  Real array of dimension ND, the evaluation point.
c            ND Integer, number of variables for the integrand.
c            NF Integer, number of components for the vector integrand.
c           Output Parameters:
c            FUNVLS Real array of dimension NF.
c                   The components of the integrand at the point X.
c   TOP    Integer, location in VERTCS array for original subregion.
c   SBRGNS Integer, number of subregions in VERTCS BEFORE subdivision.
c   VERTCS Real array of dimension (ND,0:ND,*), vertices of orginal
c          subregion must be in VERTCS(1:ND,0:ND,TOP).
c   VOLUMS Real array of dimension (*) of volumes for subregions.
c   X      Real work array of dimension ND.
c   H      Real work array of dimension ND.
c   CENTER Real work array of dimension (0:ND).
c   WORK   Real work array of dimension 5*NF.
c   FRTHDF Real work array of dimension (0:ND-1,ND).
c
c   ON RETURN
c
c   NEWSBS Integer, number of new subregions (3 or 4).
c   FUNCLS Integer, number of FUNSUB calls used by SMPDFS.
c   VERTCS Real array of dimension (ND,0:ND,*).
c          The vertices of the of new subegions will be at locations
c          TOP, SBRGNS+1, ..., SBRGNS+NEWSBS-1.
c   VOLUMS Real Array of dimension (*).
c          VOLUMS has been updated for new subregions.
c
c***ROUTINES CALLED: FUNSUB
c
c***END PROLOGUE SMPDFS
c
      external funsub
      integer nd, nf, top, sbrgns, newsbs
      double precision vertcs(nd,0:nd,*), volums(*), work(nf,*)
      double precision x(nd), h(nd), center(nd), frthdf(0:nd-1,nd)
      double precision differ, difmax, difmid, difnxt, ewidth, edgmax
      double precision cuttf, cuttb, difil, diflj, dfsmax, vti, vtj, vtl
      parameter ( cuttf = 2, cuttb = 8 )
      integer i, j, k, l, ie, je,  is, js, ls, it, jt, lt
      double precision smpvol
c
c***FIRST PROCESSING STATEMENT SMPDFS
c
c
c       Compute the differences.
c
      is = 0
      js = 1
      difmax = 0
      edgmax = 0
      do k = 1, nd
         center(k) = vertcs(k,0,top)
         do l = 1, nd
            center(k) = center(k) + vertcs(k,l,top)
         end do
         center(k) = center(k)/( nd + 1 )
      end do
      call funsub(nd, center, nf, work(1,3))
      do i = 0, nd-1
         do j = i+1, nd
            ewidth = 0
            do k = 1, nd
               h(k) = 2*( vertcs(k,i,top)-vertcs(k,j,top) )/( 5*(nd+1) )
               ewidth = ewidth + abs( h(k) )
               x(k) = center(k) - 3*h(k)
            end do
            do l = 1, 5
               do k = 1, nd
                  x(k) = x(k) + h(k)
               end do
               if ( l. ne. 3 ) call funsub(nd, x, nf, work(1,l))
            end do
            if ( ewidth .ge. edgmax ) then
               ie = i
               je = j
               edgmax = ewidth
            end if
            differ = 0
            difmid = 0
            do k = 1, nf
               difmid = difmid + abs( work(k,3) )
               differ = differ + abs( work(k,1) + work(k,5)+ 6*work(k,3)
     &                                - 4*( work(k,2) + work(k,4) ) )
            end do
            if ( difmid + differ/8 .eq. difmid ) differ = 0
            differ = differ*ewidth
            frthdf(i,j) = differ
            if ( differ .ge. difmax ) then
               it = is
               jt = js
               difnxt = difmax
               is = i
               js = j
               difmax = differ
            else if ( differ .ge. difnxt ) then
               it = i
               jt = j
               difnxt = differ
            end if
         end do
      end do
c
c  Determine whether to compute three or four new subregions.
c
      if ( difnxt .gt. difmax/cuttf ) then
         newsbs = 4
      else
         newsbs = 3
         if ( difmax .eq. 0 ) then
            is = ie
            js = je
         else
            dfsmax = 0
            do l = 0, nd
               if ( l .ne. is .and. l .ne. js ) then
                  it = min( l, is, js )
                  jt = max( l, is, js )
                  lt = is + js + l - it - jt
                  differ =  frthdf(it,lt) + frthdf(lt,jt)
                  if ( differ .ge. dfsmax ) then
                     dfsmax = differ
                     ls = l
                  end if
               end if
            end do
            difil = frthdf( min(is,ls), max(is,ls) )
            diflj = frthdf( min(js,ls), max(js,ls) )
            difnxt = difil + diflj - min( difil,diflj )
            if ( difmax/cuttb .lt. difnxt .and. difil .gt. diflj ) then
               it = is
               is = js
               js = it
            end if
         end if
      end if
c
c  Copy vertices and volume for TOP to new subregions
c
      volums(top) = volums(top)/newsbs
      do l = sbrgns + 1, sbrgns + newsbs - 1
         volums(l) = volums(top)
         do j = 0, nd
            do k = 1, nd
               vertcs(k,j,l) = vertcs(k,j,top)
            end do
         end do
      end do
      do k = 1, nd
         vti = vertcs(k,is,top)
         vtj = vertcs(k,js,top)
         if ( newsbs .eq. 4 ) then
c
c  Compute four new subregions.
c
            vertcs(k,js,top)      = ( vti + vtj )/2
            vertcs(k,is,sbrgns+1) = vti
            vertcs(k,js,sbrgns+1) = ( vti + vtj )/2
            vertcs(k,is,sbrgns+2) = ( vti + vtj )/2
            vertcs(k,js,sbrgns+2) = vtj
            vertcs(k,is,sbrgns+3) = ( vti + vtj )/2
            vertcs(k,js,sbrgns+3) = vtj
            vti = vertcs(k,it,top)
            vtj = vertcs(k,jt,top)
            vertcs(k,jt,top)      = ( vti + vtj )/2
            vertcs(k,it,sbrgns+1) = ( vti + vtj )/2
            vertcs(k,jt,sbrgns+1) = vtj
            vti = vertcs(k,it,sbrgns+2)
            vtj = vertcs(k,jt,sbrgns+2)
            vertcs(k,jt,sbrgns+2) = ( vti + vtj )/2
            vertcs(k,it,sbrgns+3) = ( vti + vtj )/2
            vertcs(k,jt,sbrgns+3) = vtj
         else
c
c  Compute three new subregions.
c
            vertcs(k,js,top)      = ( 2*vti + vtj )/3
            vertcs(k,is,sbrgns+1) = ( 2*vti + vtj )/3
            if ( difmax/cuttf .lt. difnxt ) then
               vertcs(k,js,sbrgns+1) = vtj
               vertcs(k,is,sbrgns+2) = ( 2*vti + vtj )/3
               vertcs(k,js,sbrgns+2) = vtj
               vtj = vertcs(k,js,sbrgns+1)
               vtl = vertcs(k,ls,sbrgns+1)
               vertcs(k,ls,sbrgns+1) = ( vtj + vtl )/2
               vertcs(k,js,sbrgns+2) = ( vtj + vtl )/2
               vertcs(k,ls,sbrgns+2) = vtl
            else
               vertcs(k,js,sbrgns+1) = ( vti + 2*vtj )/3
               vertcs(k,is,sbrgns+2) = ( vti + 2*vtj )/3
               vertcs(k,js,sbrgns+2) = vtj
            end if
         end if
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
