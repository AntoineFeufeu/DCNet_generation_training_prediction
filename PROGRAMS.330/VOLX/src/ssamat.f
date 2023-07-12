c-----
c       CHANGES
c
c       02 10 2000 - added a sumav1 and a sumag1 
c           for goodness of fit information
c           This will be the sum of | residual |
c           this will be an indicator of goodness of fit
c
c       20 MAY 2002  -  16 output correct mode number
c           which in internal_mode - 1
c
c       we will compute average SUM | res |
c
c       19 JUN 2009 - changed order of code near line 45 so that the
c            data statement appears after the type declarations
c-----
      subroutine ssamat(id,nfilt,inarg,wref,invcsl,nf1,invdep,
     1      lstinv,sigv,sigts)
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
c       sigv
c       sigts
c-----

c-----
c     a general purpose reformatting file
c-----
c
c     list partial derivatives or computed dispersion values
c
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
     $      h(NL2),u(NL2),ct(NL2),v(NL2,NL2),qbinv(NL),qainv(NL),
     2      wc(NL2)
        logical wc
        character*50 names
        character ostr*12
        save iwant, ihave


        common/isomod/dl(NLAY),va(NLAY),vb(NLAY),rho(NLAY),
     1      qa(NLAY),qb(NLAY),etap(NLAY),etas(NLAY), 
     2      frefp(NLAY), frefs(NLAY)
        common/depref/refdep
        common/modtit/title
        character*80 title

c-----
c       COMMON FOR WEIGHTING
c-----
        common/datvar/numr, sumrr, sumrd, rnorm,
     1                nums, sumss, sumsd, snorm

        character*1 lorr(2), porg(3)
        data lorr/'L','R'/,porg/'C','U','G'/
c-----
c       get earth and q model to set up causal partials
c---- 
        call getmod(2,'tmpsrfi.17',mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,.false.)
        m = mmax
c-----
c       save current values
c-----
        do 39 i=1,mmax
            if(qb(i).gt.1.0)then
                qbinv(i) = 1.0/qb(i)
            else
                qbinv(i) =     qb(i)
            endif
            if(qa(i).gt.1.0)then
                qainv(i) = 1.0/qa(i)
            else
                qainv(i) =     qa(i)
            endif
            b(i) = vb(i)
            a(i) = va(i)
            d(i) = dl(i)
            r(i) = rho(i)
   39   continue
c-----
c       get inversion control
c-----
        open(1,file='tmpsrfi.04',form='unformatted',access='sequential')
        rewind 1
        read(1) nd,m
        m2 = m + m
        read(1)(dd(i),i=1,m2)
        read(1)(wc(i),i=1,m2)
c-----
c       open dispersion file
c-----
        open(2,file='tmpsrfi.08',form='unformatted',access='sequential')
        rewind 2
        if(id.eq.2 .or. id.eq.5)then
c-----
c       set up file for inversion
c           tmpmrgs.9 is the input file to surfinv
c           tmpsrfi.1 is a temporary file used here only
c-----
            open(3,file='tmpmrgs.9',form='unformatted',
     1          access='sequential')
            rewind 3
c-----
c       tmpsrfi.04 is a array of observed dispersion
c-----
            open(4,file='tmpsrfi.01',form='unformatted',
     1          access='sequential',status='unknown')
            rewind 4
            write(3)m2,nfilt
            sumwav1 = 0.0
            sumwv0 = 0.0
            sumwv1 = 0.0
            sumwv2 = 0.0
            snorm = 0.0
        elseif(id.eq.22 .or. id.eq.27)then
c-----
c       set up file for dispersion listing
c-----
            if(numa.ne.0) then
                inarg=inarg+1
                call getarg(inarg,names)
            else
                write(LOT,*)'Enter File Name for Dispersion Data'
                read(LIN,'(a)') names
            endif
            open(4,file=names,form='formatted',status='unknown',
     1          access='sequential')
            rewind 4
            ifrper=0
        elseif(id.eq.12 .or. id.eq.17)then
c-----
c       list dispersion on the terminal
c-----
            write(LOT,'(a,a)')'      Mode        Period     ',
     1          '  Observed      Predicted Observed Sigma'
        elseif(id.eq.16 )then
            if(invdep.eq.1)then
            write(LOT,*)'Period, Partial with beta, Residual:'
            else if(invdep.eq.0)then
            write(LOT,*)'Period, Partial with depth, Residual:'
            endif
            write(LOT,*)m,'  Layers'
        endif
c-----
c       read in output of disper, derivl, derivr (unit 1) tmpsrfi.04
c       read in observed data (unit 2) tmpsrfi.08
c
c       Note there is no one-to-one correspondence in the between
c       observed and theoretical. There can be more theoretical than
c       observed because of missing periods at various modes.
c       In addition, observed data can have multiple observations
c       at each period for a given mode
c
c       The mitigating factor is that both are arranged in order
c       of increasing period, as
c           LOVE
c               PHASE
c                   MODES
c               GROUP
c                   MODES
c           RAYLEIGH
c               PHASE
c                   MODES
c               GROUP
c                   MODES
c
c-----
c
c       initialize the search
c-----
        tp1 = 0.0
        md1 = 0
        k1 = 0
        ihave = -1
   10   continue
c-----
c       get observed data
c-----
        read(1,end=40)ifn,k,md,tp,cobs,sd
        if(cobs.eq.0.0)go to 10
        iwant = ifn*1000 +mod(k,2)*100 + md
        call gttheo(ifn,k,md,tp,itst,k1,md1,tp1,cpred,iret,
     1      m,wref,invcsl,invdep,iwant,ihave)
        if(ihave.gt.iwant)then
            rewind 2
            go to 10
        endif
        if(iret.lt.0)then
            rewind 2
            tp1 = 0.0
            md1 = 0
            k1 = 0
            ihave = -1
            go to 10
        endif

c-----
c       if model does not support a higher mode drop from listing
c-----
        if(itst.eq.0)go to 10
c-----
c       now output matched data
c-----
        if(id.eq.2 .or. id.eq.5)then
c-----
c       calculate statistics on prediction
c-----
            dy = cobs - cpred
            sd = abs(sd)
            if(nf1.eq.1)then
                sd = 1.0
            endif

            if(k.le.2)then
c-----
c               Dispersion observation
c-----
                sumv0 = sumv0 + 1.0
                if(nf1.eq.0)then
                    sumwv0 = sumwv0 + 1.0/sd
                    sumwv1 = sumwv1 + dy/sd
                    sumwv2 = sumwv2 + dy*dy/sd
                else
                    sumwv0 = sumwv0 + 1.0
                    sumwv1 = sumwv1 + dy
                    sumwv2 = sumwv2 + dy*dy
                endif
                sumwav1 = sumwav1 + abs(dy)/sd
            endif
            if(nf1.eq.0)then
                do 15 j=1,m2
                    dd(j) = dd(j)/sd
   15           continue
            endif
CRBH HERE BE CAREFUL IN FILLING THE DD
CRBH gttheo fills the dd array - make sure that it is done correctly
            write(4)(dd(j),j=1,m2),dy,cobs/sd,cpred/sd,
     1             ifn,k,md,tp
        elseif(id.eq.17 .and. k.le.2)then
            write(LOT,'(1x,a1,1x,a1,i5,4f15.7)') lorr(ifn),porg(k)
     1          ,md1 -1 ,tp1,cobs,cpred,sd
        elseif(id.eq.16.and.k.lt.3 )then
            write(LOT,*)ifn,k,md1-1,tp1,cpred
            dy = cobs - cpred
            write(LOT,*)(dd(j),j=1,m+m),dy
        elseif(id.eq.27 .and. k.le.2)then
            ostr(1:7) = 'SURF96 '
            ostr(8:8) = lorr(ifn)
            ostr(9:9) = ' '
            ostr(10:10) = porg(k)
            ostr(11:12) = ' T'
            obserr = 0.0
            write(4,'(a12,1x,i3,1x,3f11.4)') ostr,md1-1,tp1,cpred,obserr
        elseif(id.eq.22 .and. k.eq.3)then
            ostr(1:7) = 'SURF96 '
            ostr(8:8) = lorr(ifn)
            ostr(9:9) = ' '
            ostr(10:10) = porg(k)
            ostr(11:12) = ' T'
            obserr = 0.0
            write(4,'(a12,1x,i3,1x,3f11.4)') ostr,md1-1,tp1,cpred,obserr
        endif
        go to 10
      write(LOT,*)'mismatch between observed and computed values'
      write(LOT,*)'in main program - job aborted'
      stop
  40  continue
      close(1,status='keep')
      close(2,status='keep')
      if(id.eq.2  .or. id.eq.5)then
c-----
c       rescale rows of the A matrix in Ax = b so that
c       uncoupled simultaneous inversion yields essentially the
c       same covariances as the individual inversions
c-----
            call stdev(sumwv0,sumwv1,sumwv2,vbar,sdv,0)
            if(sumv0.gt.0.0)then
                sumav1 = sumwav1 / sumwv0
            endif
        if(id.eq.2 .or. id.eq.5)then
            write(LOT,'(a)')
     1          '--------------'
            write(LOT,'(a,f10.4,a)')
     1          ' Dispersion  fit (vel)         std err  :',
     1          sdv    ,' (km/s)'
            write(LOT,'(a,f10.4,a)')
     1          ' Dispersion  fit (vel)    mean residual :',
     1          vbar   ,' (km/s)'
            write(LOT,'(a,f10.4,a)')
     1          ' Dispersion  fit (vel)    avg |residual|:',
     1          sumav1 ,' (km/s)'
        endif
c-----
c       now apply the priors sigv  for the case 
c          nf1 = 1 which means
c       that the residuals are used to estimate the variance
c-----
C        WRITE(6,*)'nf1:',nf1
        if(nf1.eq.1)then
            if(sigv .gt. sdv)then
                sdv = sigv
            endif
        endif
        call weight(sdv,wv)
c-----
c       The following is s two stage process. We attempt to
c       construct a data set which has unit variance. The
c       weights computed above are only approximate, since the
c       new data set will have a new mean, and hence possibly a
c       different variance in the case of mixed velocity and gamma
c       observations.
c
c       The rescaling is a two stage process. A rescaling is done, and
c       the variance of the rescaled data set is determined. Using this,
c       the scaling factors are modified so that the data variance
c       is in fact unity
c-----
            sum = 0.0
            sumx = 0.0
            sumxx = 0.0
            rewind 4
 1010       continue 
            read(4,end=1011)(dd(j),j=1,m2),dy,k,cobs,cpred
            if(nf1.eq.0)then
                wv = 1.0
            endif
            if(k.le.2)then
                call rescl(dd,m2,dy,wv)
            endif
            sum = sum + 1.0
            sumx = sumx + dy
            sumxx = sumxx + dy**2
            go to 1010
 1011       continue
            if(id.eq.2 .or. id.eq.5)then
                call stdev(sum,sumx,sumxx,xbar,sd,0)
            endif
            if(sd.gt.0.0 .and. nf1.eq.1)then
                wv=wv/sd
            endif
c-----
c       apply weighting and write output file for inversion
c-----
            sum = 0.0
            sumx = 0.0
            sumxx = 0.0
            nums  = 0
            sumss = 0.0
            sumsd = 0.0
            sumsw2 = 0.0
            rewind 4
 1020       continue 
            read(4,end=1021)(dd(j),j=1,m2),dy,cobs,cpred,
     1            ifn,k,md,tp

            if(k.le.2 .and. wv.gt.0.0)then
                call rescl(dd,m2,dy,wv)
                nums = nums + 1
                sumss = sumss + dy*dy
                sumsd = sumsd +(cobs*wv)**2
                sumsw2 = sumsw2 + (wv)**2
                write(3)(dd(j),j=1,m2),dy
CRBH            WRITE(6,*)'SSAMAT:',ifn,k,md,tp,(dd(j),j=1,m2),dy
c-----
c          determine the snorm
c-----
                rows = 0.0
                do j=1,m2
                    rows = rows + abs(dd(j))
                enddo
                if(rows.gt.snorm)snorm = rows
            endif
            sum = sum + 1.0
            sumx = sumx + dy
            sumxx = sumxx + dy**2
            go to 1020
 1021   continue
            close(4,status='delete')
            close(3,status='keep')
        endif
        if(id.eq.22 .or. id.eq.27) close(4,status='keep')
        return
        end

        subroutine gttheo(ifn,k,md,tp,itst,k1,md1,tp1,c1,iret,
     1  m,wref,invcsl,invdep,iwant,ihave)
        parameter (NL=200,NL2=NL+NL)
        common/ctrl/numa,d(NL),a(NL),b(NL),r(NL),rat(NL),dd(NL2),x(NL2),
     $      h(NL2),u(NL2),ct(NL2),v(NL2,NL2),qbinv(NL),qainv(NL),
     2      wc(NL2)
        logical wc
        real*4 dcdb(NL2),dcda(NL2),dgdq(NL2),dudb(NL2),duda(NL2),
     1      dudq(NL2),dcdq(NL2),dgdv(NL2),dcdh(NL2), dudh(NL2)
        save dcdb,dcda,dgdq,dudb,duda,dudq,
     1      dcdq,dgdv,dcdh,dudh
        save cvel,gvel,gam
c-----
c       search through theoretical file to find a match
c       to a particular observation
c
c       OBSERVED
c           ifn : 1 Love 2 Rayl
c           k   : 1 Phase 2 Group
c           md  : Mode 
c           tp  : Period
c       THEORETICAL
c           itst    : Whether dispersion point found
c               : 0 mode does not exist, 1 = Love, 2 = Rayleigh
c           k1  : 1 Phase 2 Group
c           md1 : Mode
c           tp1 : Period
c           c1  : Velocity
c           dd  : Partial Derivatives
c-----
c
c       check for possibility of a double observation
c
c-----
c-----
c       initialize the dd array
c-----
c-----
c       test to see if request is for information already in hand
c       NOTE: the test mod(k,2) == mod(k1,2) arises since
c       we only need dcdb, phase velocity and perehaps dcda to get
c       the anelastic attenutation coefficient gamma
c-----
        iret = 0
        dif = abs( (tp-tp1)/tp )
        jret = -1
        if(ifn.eq.itst .and. md.eq.md1 .and. dif.lt.1.0e-5)then
            if(mod(k,2).eq.mod(k1,2))then
c-----
c       should be .eq. 1 but to permit velocity or gamma
c       dispersion immediately after dispersion computation
c       use .ne. 0
c-----
                if(invdep.ne.0)then
                    if(k.eq.1)then
                        call zero(dd,1,m+m)
                        call cpy(dcdb,dd,0,m)
                        if(itst.eq.2)then
CRBH        WRITE(6,*)'C RAYLEIGH DCDA'
                           call cpy(dcda,dd,m,m)
                        endif
                        c1 = cvel
                    elseif(k.eq.2)then
                        call zero(dd,1,m+m)
                        call cpy(dudb,dd,0,m)
                        if(itst.eq.2)then
                           call cpy(duda,dd,m,m)
CRBH        WRITE(6,*)'U RAYLEIGH DUDA'
                        endif
                        c1 = gvel
                    endif
                else if(invdep.eq.0)then
                    if(k.eq.1)then
                        call zero(dd,1,m+m)
                        call cpy(dcdh,dd,0,m)
                        c1 = cvel
                    elseif(k.eq.2)then
                        call zero(dd,1,m+m)
                        call cpy(dudh,dd,0,m)
                        c1 = gvel
                    endif
                endif
                jret=1
            endif
        endif
c-----
c       test to see if information is in hand for a mode
c       not computed
c       E.g., 1 st mode not found, but we are asked for
c       phase and gamma, so read in another observed data
c       value rather than reading another theoretical value
c       and getting the two lists out of synchronization
c-----
        if(itst.eq.0 .and. md.eq.md1 .and. dif.lt.1.0e-5
     1      .and. mod(k,2).eq.mod(k1,2) )then
            jret = 1
        endif
        if(jret.eq.1)return
c-----
c       since we will not reuse previous value
c       read in new theoretical observation
c-----
   10   continue
        iret = 0
        read(2,end=30) itst,k1,md1,tp1
c-----
c      itst   - 1, Love, 2 Rayleigh
c      k1     - 1 phase velocity, 2 group velocity sequence
c      md1    - mode
c      tp1    - period
c-----
CRBH        WRITE(6,*)'SW  :', itst,k1,md1,tp1
        ihave = 1000*itst + 100*mod(k1,2) + md1
        dif = abs ( (tp-tp1)/tp)
        if(k1.eq.2)then
            read(2,end=30)gvel,(dudb(j),j=1,m)
            read(2,end=30) (dudh(j),j=1,m)
            read(2,end=30) cvel,(dcdb(j),j=1,m)
            read(2,end=30) (dcdh(j),j=1,m)
CRBH        WRITE(6,*)'   U:',gvel 
CRBH        WRITE(6,*)'DUDB:',(dudb(j),j=1,m)
CRBH        WRITE(6,*)'DUDH:',(dudh(j),j=1,m)
CRBH        WRITE(6,*)'   C:',cvel 
CRBH        WRITE(6,*)'DCDB:',(dcdb(j),j=1,m)
CRBH        WRITE(6,*)'DCDH:',(dcdh(j),j=1,m)
            if(itst.eq.2)then
                 read(2,end=30)(dcda(j),j=1,m)
                 read(2,end=30)(duda(j),j=1,m)
CRBH        WRITE(6,*)'DCDA:',(dcda(j),j=1,m)
CRBH        WRITE(6,*)'DUDA:',(duda(j),j=1,m)
            endif
        else
            read(2,end=30)cvel,(dcdb(j),j=1,m)
            read(2,end=30) (dcdh(j),j=1,m)
CRBH        WRITE(6,*)'   C:',cvel 
CRBH        WRITE(6,*)'DCDB:',(dcdb(j),j=1,m)
CRBH        WRITE(6,*)'DCDH:',(dcdh(j),j=1,m)
            if(itst.eq.2)then
                 read(2,end=30)(dcda(j),j=1,m)
CRBH        WRITE(6,*)'DCDA:',(dcda(j),j=1,m)
            endif
        endif
        if(mod(k,2).ne.mod(k1,2).or.md.ne.md1.or.dif.gt.1.0e-5)then
            go to 10
        endif
c-----
c       define correct parameter
c-----
        call getgam(dcdb,b,qbinv,dcda,a,qainv,m,itst,cvel,tp,
     1      gam,dgdq)
        call zero(dcdq,1,m)
        call zero(dudq,1,m)
        call zero(dgdv,1,m)
        if(invcsl.ge.1)call cslmod(gvel,cvel,gam,wref,tp,m,
     1      qbinv,qainv,dcdb,dudb,dudb,itst,dcdq,dudq,dgdv,dcda
     2      ,a,b,k1)
c-----
c       setup causal, uncoupled if invcsl = 1
c-----
        if(invcsl.eq.1 .or. invdep.eq.0)then
            call zero(dudq,1,m)
            call zero(dcdq,1,m)
            call zero(dgdv,1,m)
        endif
c-----
c       initialize the partial array
c-----
        call zero(dd,1,m+m)
        if(invdep.ne.0)then
c-----
c           invert for velocity - here first half of dd array
c           is S-velocity partials and the second is P partial
c           we must be careful with respec to the type
c=====
            if(k1.eq.2)then
                call zero(dd,1,m+m)
                call cpy(dudb,dd,0,m)
                if(itst.eq.2)then
                    call cpy(duda,dd,m,m)
CRBH        WRITE(6,*)'U RAYLEIGH DUDA'
CRBH        WRITE(6,*)(dd(i),i=1,m+m)
                endif
                c1 = gvel
            elseif(k1.eq.1)then
                if(k.eq.1)then
                    call zero(dd,1,m+m)
                    c1 = cvel
                    call cpy(dcdb,dd,0,m)
                    if(itst.eq.2)then
CRBH        WRITE(6,*)'C RAYLEIGH DCDA'
                        call cpy(dcda,dd,m,m)
CRBH        WRITE(6,*)(dd(i),i=1,m+m)
                    endif
                endif
            endif
        else if(invdep.eq.0)then
c-----
c           invert for layer thickness
c-----
            if(k1.eq.2)then
                call zero(dd,1,m+m)
                call cpy(dudh,dd,0,m)
                c1 = gvel
            elseif(k1.eq.1)then
                if(k.eq.1)then
                    call zero(dd,1,m+m)
                    c1 = cvel
                    call cpy(dcdh,dd,0,m)
                endif
            endif
        endif
        if(ihave.lt.iwant)go to 10
        return
   30   continue
            iret = -1
        return
        end

        subroutine cslmod(u,c,gam,wref,tp,m,
     1      qbinv,qainv,dcdb,dudb,duda,itst,dcdq,dudq,dgdv,dcda
     2      ,a,b,k1)
        real*4 dcdb(*),dudb(*),duda(*),qbinv(*),dcdq(*),dudq(*)
     1      ,dgdv(*),dcda(*),a(*),b(*),qainv(*)
c-----
c       modify phase, group velocities and partials for
c       constant Q causality
c
c       itst 1 = L, 2 = R
c-----
        c0 = c
        u0 = u
        f  = 1./tp
        pi = 3.1415927
        omega = 6.2831853*f
        faclog = alog(f/wref)/pi
        facg = 2.*gam*u0/(pi*omega)
c-----
c       get correct phase velocity
c-----
        c = c0 + 2.*gam*c0*c0*faclog/omega
c-----
c       get correct group velocity
c-----
        cmc0c0 = (c-c0)/c0
        if(k1.eq.2)then
            u0c0 = u0/c0
            u = u0 *( 1. + (2. - u0c0)*cmc0c0 + facg ) 
            uu0 = u/u0
        endif
c-----
c       get correct partials
c-----
        do 100 i=1,m
            if(b(i).eq.0.0)then
                qaiqbi = 1.0e-6
            else 
                if(qbinv(i).gt.0.0)then
                    qaiqbi = qainv(i)/qbinv(i)
                else
                    qaiqbi = 1.0/(0.75*(a(i)/b(i))**2)
                endif
            endif
            dcdp = 0
            if(itst.eq.2)dcdp=dcda(i)
            dcds = dcdb(i)
            dcdb(i) = dcds*(1. + faclog*qbinv(i) )
            if(k1.eq.2)then
            dudv = dudb(i)
            dudb(i) = dudv*(uu0 -u0c0*cmc0c0 + facg )
     1       + dcds*u0c0*(-2.*facg +u0c0*qbinv(i)/pi +
     2          (2.-u0c0)*(faclog*qbinv(i) -cmc0c0) +
     3          u0c0 *cmc0c0 )
            dudv = duda(i)
            duda(i) = dudv*(uu0 -u0c0*cmc0c0 + facg )
     1       + dcds*u0c0*(-2.*facg +u0c0*qbinv(i)/pi +
     2          (2.-u0c0)*(faclog*qbinv(i) -cmc0c0) +
     3          u0c0 *cmc0c0 )
            dudq(i) = u0c0*(2.-u0c0)*dcdq(i) +
     1          u0c0*u0c0*(dcds*b(i)+dcdp*a(i)*qaiqbi)/pi
     2          - 2*gam*dcds/c0
            endif
            dcdq(i) = faclog*(dcds*b(i) +dcdp*a(i)*qaiqbi)
            dgdv(i) = 0.5*omega*dcds*qbinv(i)/(c0*c0)
  100   continue
        return
        end

        subroutine getgam(dcdb,b,qbinv,dcda,a,qainv,m,ilvry,cc,tp,
     1          gam,dgdq)
        real*4 dcdb(*),b(*),dcda(*),a(*),dgdq(*),qbinv(*), qainv(*)
c-----
c       determine spatial attenuation at constant frequency
c       and partial derivatives
c-----
        omega = 6.2831853/tp
        factr = 0.5*omega/(cc*cc)
        sum = 0.0
        do 100 i=1,m
            sum = sum + dcdb(i)*b(i)*qbinv(i)
            dgdq(i) = factr * dcdb(i)*b(i)
            if(ilvry.eq.2)then
                if(b(i).eq.0.0)then
                    qaiqbi = 1.0e-6
                else 
                    if(qbinv(i).gt.0.0)then
                        qaiqbi = qainv(i)/qbinv(i)
                    else
                        qaiqbi = 0.75*(a(i)/b(i))**2
                    endif
                endif
                sum = sum + dcda(i)*a(i)*qainv(i)
                dgdq(i) = dgdq(i) + factr*dcda(i)*a(i)*qaiqbi
            endif
  100   continue
        gam = factr * sum
        return
        end

        subroutine weight(sdv,wv)
c-----
c       given standard error estimates of the velocity and
c       gamma residuals, determine weights for velocity and
c       gamma data so that the covariance matrix of the
c       uncoupled simultaneous solution will be identical to
c       that for the velocity alone and q alone inversions
c       for the velocity and q model covariances, respectively
c
c       This may be viewed as a guarantee that the same 
c          physical units are
c       used is estimating the standard deviation of fit
c
c       The weights are inversely proportional to the standard error
c       estimates so that the normalized Variance is close to unity
c-----
c-----
c       check for infinite standard errors
c-----
        if(sdv.le.0.0 )then
            wv = 0.0
        else
            wv = 1.0/sdv
        endif
        return
        end

        subroutine rescl(x,n,y,wt)
c-----
c       rescale x vector and y by sd if sd > 0
c-----
        real*4 x(*)
        if(wt.gt.0.0)then
            do 100 i=1,n
                x(i)=x(i)*wt
  100       continue
            y = y*wt
        endif
        return
        end

        subroutine zero(x,m,n)
        real*4 x(*)
c-----
c       set x(m) = 0, x(m+1)=0, ..., x(n)=0
c-----
        do 100 i=m,n
            x(i) = 0.0
  100   continue
        return
        end

        subroutine cpy(x,y,m,n)
        real*4 x(*),y(*)
c-----
c       copy contents of array x into array y with offset m
c-----
        do 100 i=1,n
            y(i+m) = x(i)
  100   continue
        return
        end

        subroutine stdev(sum0,sumx,sumxx,sbar,sd,m)
c-----
c       estimate standard error - note since
c       this must also work with weighted data
c       we will use an N in the denominator instead of N-1
c-----
c       sum = number of observations
c       sumx= sum of x(i)
c       sumxx= sum of x(i)**2
c-----
        if(sum0.le.0.0)then
            if(sum0 .le. 0.0)then
                sbar = 0.0
            else
                sbar = sumx/sum0
            endif
            sd = -1.0
        else
            sbar = sumx/sum0
            sd = sqrt(abs( sumxx - sum0*sbar*sbar)/sum0)
        endif
        return
        end

