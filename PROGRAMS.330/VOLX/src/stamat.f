c-----
c       CHANGES
c
c-----
        subroutine stamat(id,nfilt,inarg,wref,invcsl,nf1,invdep,
     1      lstinv,sigtp)
c-----
c       id - command
c       nfilt
c       inarg
c       wref        - reference frequency for causality
c       invcsl
c       nf1
c       invdep      - 1  invert for velocity change
c                     0  invert for layer thickness change
c       lstinv
c       sigtp
c-----

c-----
c     a general purpose reformatting file
c-----
c
c     list partial derivatives or computed travel time values
c-----
      parameter(LER=0,LIN=5,LOT=6)
        integer NL, NL2, NLAY
        parameter(NL=200,NLAY=200,NL2=NL+NL)

c-----
c       LIN - unit for FORTRAN read from terminal
c       LOT - unit for FORTRAN write to terminal
c       LER - unit for FORTRAN error output to terminal
c       NL  - number of layers in model
c       NL2 - number of columns in model (first NL2/2 are
c           - velocity parameters, second NL2/2 are Q values)
c-----
        common/ctrl/numa,d(NL),a(NL),b(NL),r(NL),rat(NL),dd(NL2),x(NL2),
     $      h(NL2),u(NL2),c(NL2),v(NL2,NL2),qbinv(NL),qainv(NL),
     2      wc(NL2)
        logical wc

        real gaussalp
        integer nlow,nhgh,mm,mm2
        real sum2o, sum2p, sum2r,redv

        common/isomod/dl(NLAY),va(NLAY),vb(NLAY),rho(NLAY),
     1      qa(NLAY),qb(NLAY),etap(NLAY),etas(NLAY), 
     2      frefp(NLAY), frefs(NLAY)
        common/depref/refdep
        common/modtit/title
        character*80 title

        character kstnm*8
        integer nzyear, nzjday, nzhour, nzmin, nzmon, nzday

        real invwgt
        real tdd(NL2)
c-----
c       COMMON FOR WEIGHTING
c-----
        common/datvar/numt, sumrr, sumrd, rnorm,
     1                nums, sumss, sumsd, snorm

c-----
c       get earth and q model to set up causal partials
c---- 
        call getmod(2,'tmpsrfi.17',mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,.false.)
        m = mmax
c-----
c       open prediction/partial file
c-----
        if(id.eq.2.or.id.eq.4)then
            call setupmat(id,nfilt,inarg,wref,invcsl,nf1,invdep,
     1      lstinv,sigtp)
        else if(id.eq.11. or. id.eq.12)then
            call lstttime(id,nfilt,inarg,wref,invcsl,nf1,invdep,
     1      lstinv,sigtp)
        endif
        return
        end

       subroutine lstttime(id,nfilt,inarg,wref,invcsl,nf1,invdep,
     1      lstinv,sigtp)
c-----
c      list travel time information
c-----
        parameter(LER=0,LIN=5,LOT=6)
        integer NL, NL2, NLAY
        parameter(NL=200,NLAY=200,NL2=NL+NL)

        real tdd(NL2)
        character*2 ostr

c-----
c       get inversion control
c-----
        call ddcur(m,m2)
            open(4,file='tmpsrfi.21',form='unformatted',
     1          access='sequential',status='unknown')
            rewind 4
        if(id.eq.11)then
              if(invdep.eq.1) then
                   WRITE(LOT,*)'Partial with respect to velocity'
              else if(invdep.eq.1) then
                   WRITE(LOT,*)'Partial with respect to layer thickness'
              endif
        else if(id.eq.12)then
              WRITE(LOT,*)
     1' Ph  Dist(km)  Tobs(s)   Tpre(s)   StdObs(s) Evdp(km)  Stdp(km)'
        endif
 1000   continue
        read(4,end=9000) iphase,dist, (tdd(i),i=1,m),obs,pre,err,
     1        evdp, stdp,dy
        if(iphase.eq.0)then
              ostr = 'P '
        else if(iphase.eq.3)then
              ostr = 'S '
        else if(iphase.eq.1)then
              ostr = 'SV'
        else if(iphase.eq.2)then
              ostr = 'SH'
        endif
        if(id.eq.11)then
            WRITE(LOT,'(1x,a2,2f10.4)')
     1         ostr,dist,pre
            WRITE(LOT,*)(tdd(i),i=1,m),dy
        else if(id.eq.12)then
           WRITE(LOT,'(1x,a2,6f10.4)')
     1         ostr,dist, obs,pre,err,evdp,stdp
      
        endif
        go to 1000
 9000   continue
         close(4)
      return
      end

       subroutine setupmat(id,nfilt,inarg,wref,invcsl,nf1,invdep,
     1      lstinv,sigtp)
c-----
c       id - command
c       nfilt
c       inarg
c       wref        - reference frequency for causality
c       invcsl
c       nf1
c       invdep      - 1  invert for velocity change
c                     0  invert for layer thickness change
c       lstinv
c       sigtp
c-----

c-----
c     set up travel time inversion equations
c-----
      parameter(LER=0,LIN=5,LOT=6)
        integer NL, NL2, NLAY
        parameter(NL=200,NLAY=200,NL2=NL+NL)

c-----
c       LIN - unit for FORTRAN read from terminal
c       LOT - unit for FORTRAN write to terminal
c       LER - unit for FORTRAN error output to terminal
c       NL  - number of layers in model
c       NL2 - number of columns in model (first NL2/2 are
c           - velocity parameters, second NL2/2 are Q values)
c-----
        common/ctrl/numa,d(NL),a(NL),b(NL),r(NL),rat(NL),dd(NL2),x(NL2),
     $      h(NL2),u(NL2),c(NL2),v(NL2,NL2),qbinv(NL),qainv(NL),
     2      wc(NL2)
        logical wc

        real gaussalp
        integer nlow,nhgh,mm,mm2
        real sum2o, sum2p, sum2r,redv

        common/isomod/dl(NLAY),va(NLAY),vb(NLAY),rho(NLAY),
     1      qa(NLAY),qb(NLAY),etap(NLAY),etas(NLAY), 
     2      frefp(NLAY), frefs(NLAY)
        common/depref/refdep
        common/modtit/title
        character*80 title

        character kstnm*8
        integer nzyear, nzjday, nzhour, nzmin, nzmon, nzday

        real invwgt
        real tdd(NL2)
c-----
c       COMMON FOR WEIGHTING
c-----
        common/datvar/numt, sumrr, sumrd, rnorm,
     1                nums, sumss, sumsd, snorm

c-----
c       get inversion control
c-----
        call ddcur(m,m2)
c-----
c           tmpsrfi.01 is for temporary storage 
c-----
            open(2,file='tmpsrfi.01',form='unformatted',
     1          access='sequential',status='unknown')
            rewind 2
c-----
c       set up file for inversion
c           tmpmrgt.9 is the input file to surfinv
c-----
            open(3,file='tmpmrgt.9',form='unformatted',
     1          access='sequential')
            rewind 3
            open(4,file='tmpsrfi.21',form='unformatted',
     1          access='sequential',status='unknown')
            rewind 4
            write(3)m2,nfilt
            sumv0 = 0.0
            sumv1 = 0.0
            sumv2 = 0.0
        numt = 0
        numPt = 0
        numSt = 0

        sumwt = 0.0
        sumwr = 0.0
        sumwr2 = 0.0
        sumrr = 0.0
        sumrd = 0.0
        sumrw2 = 0.0
        sumwat1 = 0.0

        Psumwt = 0.0
        Psumwr = 0.0
        Psumwr2 = 0.0
        Psumrr = 0.0
        Psumrd = 0.0
        Psumrw2 = 0.0
        Psumwat1 = 0.0

        Ssumwt = 0.0
        Ssumwr = 0.0
        Ssumwr2 = 0.0
        Ssumrr = 0.0
        Ssumrd = 0.0
        Ssumrw2 = 0.0
        Ssumwat1 = 0.0

        rnorm = 0.0
 1000   continue
        read(4,end=9000) iphase,dist, (tdd(i),i=1,m),obs,pre,err,
     1        evdp, stdp,dy
C        if(err.gt.0.0)then
C                 WT = 1./(err)
C        else
                 WT = 1.0
C        endif
CRBH        WRITE(6,*) iphase,dist,obs,pre,err,evdp, stdp,dy
CRBH        WRITE(6,*) (tdd(i),i=1,m)
                do 1101 ii=1,m
                    tdd(ii) = tdd(ii)*wt
 1101           continue
                call zero(dd,1,m2)
                do i=1,m
                      if(invdep.eq.1)then
c-----
c                       invert for velocity and
c                       set up partials with respect to S and P
c-----
                        if(iphase.eq.0)then
c-----
c                       P is second half of the row
c-----
                            dd(i+m) = tdd(i)
                            Psumwt = Psumwt + wt
                            Psumwr  = Psumwr  + dy*wt
                            Psumwr2 = Psumwr2 + dy*dy*wt
                            Psumrr = Psumrr + ( dy*wt)**2
                            Psumrd = Psumrd + (obs*wt)**2
                            Psumrw2 = Psumrw2 + (wt)**2
                            Psumwat1 = Psumwat1 + abs(dy*wt)
                            numPt = numPt + 1
                        else
c-----
c                       S is the first half of the row
c-----
                            dd(i) = tdd(i)
                            Ssumwt = Ssumwt + wt
                            Ssumwr  = Ssumwr  + dy*wt
                            Ssumwr2 = Ssumwr2 + dy*dy*wt
                            Ssumrr = Ssumrr + ( dy*wt)**2
                            Ssumrd = Ssumrd + (obs*wt)**2
                            Ssumrw2 = Ssumrw2 + (wt)**2
                            Ssumwat1 = Ssumwat1 + abs(dy*wt)
                            numSt = numSt + 1
                        endif
                      else if(invdep.eq.0)then
c-----
c                       set up partials with respect to layer thickness
c                       which occupy first m positions of array
c-----
                            dd(i) = tdd(i)
                        if(iphase.eq.0)then
c-----
c                       P is second half of the row
c-----
                            Psumwt = Psumwt + wt
                            Psumwr  = Psumwr  + dy*wt
                            Psumwr2 = Psumwr2 + dy*dy*wt
                            Psumrr = Psumrr + ( dy*wt)**2
                            Psumrd = Psumrd + (obs*wt)**2
                            Psumrw2 = Psumrw2 + (wt)**2
                            Psumwat1 = Psumwat1 + abs(dy*wt)
                            numPt = numPt + 1
                        else
c-----
c                       S is the first half of the row
c-----
                            Ssumwt = Ssumwt + wt
                            Ssumwr  = Ssumwr  + dy*wt
                            Ssumwr2 = Ssumwr2 + dy*dy*wt
                            Ssumrr = Ssumrr + ( dy*wt)**2
                            Ssumrd = Ssumrd + (obs*wt)**2
                            Ssumrw2 = Ssumrw2 + (wt)**2
                            Ssumwat1 = Ssumwat1 + abs(dy*wt)
                            numSt = numSt + 1
                        endif
                      endif
                enddo
                write(2)(dd(j),j=1,m2),dy*wt
                sumwt = sumwt + wt
                sumwr  = sumwr  + dy*wt
                sumwr2 = sumwr2 + dy*dy*wt
                sumrr = sumrr + ( dy*wt)**2
                sumrd = sumrd + (obs*wt)**2
                sumrw2 = sumrw2 + (wt)**2
                sumwat1 = sumwat1 + abs(dy*wt)
                numt = numt + 1
        go to 1000
 9000   continue
c-----
c       estimate standard deviation
c-----
c       if unweighted
c
c            N SUM x^2 - ( SUM X ) ^2
c       sd = ------------------------
c                     N^2
c       
c       for weighted
c
c            SUM w ( SUM w x^2 ) - ( SUM w  X ) ^2
c       sd = ------------------------
c                     (SUM w ) ^ 2
c       note this would give the proper units in terms of the
c       orignial receiver function  
c-----
c       adjust everything so that I get correct variance
c       this is done by scaling sumwt up to numt
c-----
        tbar = sumwr/sumwt
        sdr = ( sumwt * sumwr2 - sumwr * sumwr * tbar )/( sumwt*sumwt)
        sdr = abs(sdr)
        sumat1 = sumwat1 / sumwt
        if(numPt .gt. 0)then
             Ptbar = Psumwr/Psumwt
             Psdr = ( Psumwt * Psumwr2 - Psumwr * Psumwr * Ptbar )/
     1           ( Psumwt*Psumwt)
             Psdr = abs(Psdr)
             Psumat1 = Psumwat1 / Psumwt
            write(LOT,'(a)')
     1          '--------------'
             write(LOT,'(a,f10.4,a)')
     1      ' P Travel time fit             std err  :',
     1      Psdr, ' (sec)'
             write(LOT,'(a,f10.4,a)')
     1      ' P Travel time fit        mean residual :',
     1      Ptbar   ,' (sec)'
             write(LOT,'(a,f10.4,a)')
     1      ' P Travel time fit        avg |residual|:',
     1      Psumat1 ,' (sec)'
        endif
        if(numSt .gt. 0)then
             Stbar = Ssumwr/Ssumwt
             Ssdr = ( Ssumwt * Ssumwr2 - Ssumwr * Ssumwr * Stbar )/
     1           ( Ssumwt*Ssumwt)
             Ssdr = abs(Ssdr)
             Ssumat1 = Ssumwat1 / Ssumwt
            write(LOT,'(a)')
     1          '--------------'
             write(LOT,'(a,f10.4,a)')
     1      ' S Travel time fit             std err  :',
     1      Ssdr, ' (sec)'
             write(LOT,'(a,f10.4,a)')
     1      ' S Travel time fit        mean residual :',
     1      Stbar   ,' (sec)'
             write(LOT,'(a,f10.4,a)')
     1      ' S Travel time fit        avg |residual|:',
     1      Ssumat1 ,' (sec)'
        endif


c-----
c       now weight by expected std err - not this could eventually
c       be done by time point
c-----
c-----
c       now apply the priors sigtp 
c-----
            if(sigtp .gt. sdr)then
                sdr = sigtp
            endif
c-----
c       in applying the correction to get a unit variance, we
c       actually need to use 1, w dy and ( w dy )^2 for sums
c       since we have already corrected for things but mess to use the
c       floor
c
c       or,
c       the computed sdr is from the entire data set with mixed alphas.
c       I have tried to correct for the alpha effect on anplitude by
c       dividing the partials and residuals by alpha before computing
c       sdr - Now when I compare to the user sigtp the question is
c       what does sigtp mean - It is the expected sigma on the
c       gauss alpha = 1.0 traces.  Note this is OK since
c       by using many traces the sigma should be on the distribution
c       of the scatter between traces and not the scatter on the mean,
c       If only one stacked function is fit then the sigma should be
c       from the variability of the stack
c-----
        fac = ( numt * sumrr - sumwr*sumwr )/(numt*numt)
        fac = sqrt(abs(fac))
        sdr = sdr/fac
        rewind (2)
 2000   continue
                read(2,end=9100)(dd(j),j=1,m2),dy
                dy = dy / sdr
                do 2100 i=1,m2
                    dd(i) = dd(i) / sdr
 2100           continue
                write(3)(dd(j),j=1,m2),dy
CRBH                WRITE(6,*)(dd(j),j=1,m2),dy
c-----
c       get the row norm
c-----
        rowr = 0.0
        do j=1,m2
           rowr = rowr + abs(dd(j))
        enddo
        if(rowr .gt. rnorm) rnorm = rowr
        go to 2000
 9100   continue
c-----
c       close output devices
c-----
C        close(1)
        close(2)
        close(3)
        close(4)
        return
        end

