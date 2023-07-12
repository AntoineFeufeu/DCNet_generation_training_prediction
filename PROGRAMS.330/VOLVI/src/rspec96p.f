        program hspec96p
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME VI                                                      c
c                                                                     c
c      PROGRAM: RSPEC96P                                              c
c                                                                     c
c      COPYRIGHT 2019                                                 c
c      R. B. Herrmann                                                 c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63013                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
c       CHANGES
c       11 SEP 2000 - build in P, SV and SH first arrival times
c       14 JAN 2001 - put in A, C, F, L and N constants into
c           trace header
c       08 JUN 2001 - error in mapping of ksrc - fixed to
c           ksrc(i) = jsrc(lsrc(i))
c       18 JUL 2001 - fixed swpecification for KSRC when source and
c           receiver are in fluids. 
c       09 APR 2002 - another hack at fixing the 
c       JSRC when fluids are involved.
c           commented out lines in gethsp which was not correct, 
c           modified output loop
c       05 MAY 2002 - corected putting A C F L and N into header
c       29 OCT 2002 - Fixed evalg and evalh for fluid layers
c       29 JAN 2003 - corrected determination of Source/Receiver 
c               velocities by using srclay and not srclyr
c
c       THIS IS BUILT ON hspec96p BUT USES THE GENERALIZED REFLECTION
c       TRANSMISSION MATRICES INSTEAD OF PROPAGATOR MATRICES. FOR TH#E
c       P-SV PROBLEM THIS MEANS THAT IT CAN HANDLE ARBITRARY SEQUENCES
c       OF FLUID AND SOLID LAYERS. THE PART WILL GIVE RESULTS IF THERE
c       ARE NO INTERVENING FLUID LAYERS BETWEEN THE SOURCE AND RECEIVER
c       21 AUG 2019 - implemented. Note that there is some redundancy
c                     in the coding so that the net parts can be
c                     easily placed into hspec96 to make mspec96
c
c                     Since the result is a ptau time series, there is not
c                     first arrival time
c    04 JAN 2021 - change all do XXX yyy ... XXX a=b 
c                  to   do yy ... a=b .. enddo
c           to be compatible with gfortran 9.3
c-----
        parameter(LER=0, LIN=5, LOT=6)
        parameter(NL=200)
        parameter(NSOURCE=100,NRECEIVER=100,NSR=100)
        common/source/depths(NSOURCE),lmaxs(NSOURCE),mdpths
        common/receiv/depthr(NRECEIVER),lmaxr(NRECEIVER),mdpthr
        real*4 depths, depthr
        integer lmaxs, lmaxr, mdpths, mdpthr
        common/bound/vmin,vamin,vamax,vbmin,vbmax
        common/damp/alpha,ieqex
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        common/modelt/dt(NL),at(NL),bt(NL),rhot(NL),
     1      mmaxt,qat(NL),qbt(NL),etapt(NL),etast(NL),
     1      frefpt(NL),frefst(NL)
        common/jout/jsrc(21) , jbdrys, jbdryh
c-----
c       matrix components in layers and boundaries saved
c-----
        common/water/iwater(NL),iwats(2),iwatb(2)
        common/updnsm/equalu(NL), equald(NL)
        logical equalu, equald
        common/lyrctl/lyrins
        logical lyrins
        common/rlimit/rlim
        real*4 rlim
c-----c
        real*4 ffreq(8)
        parameter(NDST=500)
        real*4 r(NDST), tshift(NDST), vred(NDST)
        complex smo(21)
        complex ztmp, zdata
        common/c/pmin,pmax,dp,pcntrl
        common/frlim/fl,fu,df,fwhich
        common/cntrl/ishank,hnkarg,dstcor,dokjar
c-----
c       ishank  L   - .true. use Hankel function and not Bessel
c                    not that asymptotic tricks are not used
c       hnkarg  R*4 - (default 6.0) For kr > hnkarg use the Hankel
c                   function, otherwise the Bessel
c       isfref  L   - .true. reference frequency uis not the 
c                   default of 1.0 Hz
c       fref    R*4 - reference frequency for Causal Q
c       dstcor  R*4 - 0 first time sample is tshift + r/vred
c                 1 first time sample is tshift + z/vred where
c                   z = abs (source depth - receiver depth)
c                 2 first time sample is tshift + 
c                   sqrt(z*z + r*r) /vred
c       dokjar  L   - .false. use Futterman causal Q
c                 .true.  use Kjartansson Causal Q
c-----
        integer*4 dstcor
        real*4 hnkarg
        logical ishank,dokjar
        common/ctlrud/rpup,rpdn,rsup,rsdn,dorud
        logical rpup, rpdn, rsup, rsdn, dorud
        common/ctlsud/spup,spdn,ssup,ssdn,dosud
        logical spup, spdn, ssup, ssdn, dosud

        logical ext
        character mname*80, title*80

        common/depref/refdep
        common/srcinf/vsa, vsb, vsr

C       COMMON/DEBUG/VERBY
C       LOGICAL VERBY
c-----
c       ksrc = temporary array for jsrc for output
c       here jsrc != 0 reflects source information
c       if receiver is not in a fluid DO NOT output pressure field
c-----
        integer ksrc(21)
c-----
c       lsrc maps jsrc to output Green s functions. e.g., if
c       jsrc(8) = radial explosion term, but in final output it
c       occupies position 10, or jsrc(lsrc(10)) = computed
c-----
        integer*4 lsrc(21)
        data lsrc/1,2,3,4,13,5,6,14,7,8,9,10,11,12,15,16,17,18,19,20,21/
        fwhich = -1.0
c-----
c       call machine dependent initialization
c-----
        call mchdep()
c-----
c       parse command line arguments
c-----
        call gcmdln()
c-----
c       set up tolerances
c           rlim is the distance which is effectively zero
c-----
            rlim = 1.0e-07
c-----
c       See if the file hspec96p.dat exists, if it does
c       open it for all control information
c-----
        inquire(file='hspec96p.dat',exist=ext)
        if(.not. ext)then
                write(LER,*)'Control file ',
     1          'hspec96p.dat does not exist'
                go to 9999
        endif
        open(3,file='hspec96p.dat',access='sequential',
     1      form='formatted',status='unknown')
        rewind 3
        call gethsp(3,mname,fl,fu,delt,n,n1,n2,xleng,xfac,
     1      ndist,r,tshift,vred)
        close (3)
c-----
c       process
c-----
        df = 1./(n*delt)
        nyq = n/2 + 1
        nyq2 = 2*nyq
        write(LOT,2)  fl,fu,df,n1,n2,n,
     2      vmin,vamin,vamax,vbmin,vbmax
c-----
c       UNIX output - no carriage control
c-----
    2   format('fl =',f10.5,5x,'fu =',f10.5,5x,'df =',f10.5,/
     1      4x,'n1 =',i4,5x,'n2 =',i4,5x, ' n =',i5/
     2  'vmin  =',f10.5,' vamin =',f10.5,' vamax =',f10.5/
     3  '       ',10x  ,' vbmin =',f10.5,' vbmax =',f10.5)
 2021   format('depths =',f20.6)
 2031   format('depthr =',f14.6)
 2041   format('     r =',f14.6,' tshift=',f10.2,' vred=',f10.2)
    3   format(8f10.5)
    4   format('frequencies for which response computed     ')
    5   format('alpha =',f10.5,5x,'dt =',f10.3)
c-----
c       MSDOS MICROSOFT FORTRAN OUTPUT
c-----
c    2   format(1x,'fl =',f10.5,5x,'fu =',f10.5,5x,'df =',f10.5,/
c     1          1x,4x,'n1 =',i4,5x,'n2 =',i4,5x, ' n =',i5/
c     2  1x,'vmin  =',f10.5,' vamin =',f10.5,' vamax =',f10.5/
c     3  1x,'       ',10x  ,' vbmin =',f10.5,' vbmax =',f10.5)
c 2021   format(1x,'depths =',f20.6)
c 2031   format(1x,'depthr =',f14.6)
c 2041   format(1x,'     r =',f14.6,' tshift=',f10.2,' vred=',f10.2)
c    3   format(1x,8f10.5)
c    4   format(1x,'frequencies for which response computed     ')
c    5   format(1x,'alpha =',f10.5,5x,'dt =',f10.3)
c-----
        write(LOT,*)'SOURCE DEPTH (',mdpths,')'
        do 2020 i=1,mdpths
            write(LOT,2021)depths(i)
 2020   continue
        write(LOT,*)'RECEIVER DEPTH (',mdpthr,')'
        do 2030 i=1,mdpthr
            write(LOT,2031)depthr(i)
 2030   continue
        write(LOT,*)'RECEIVER DISTANCES (',ndist,')'
        do 2040 i=1,ndist
            write(LOT,2041)r(i), tshift(i), vred(i)
 2040   continue
        write(LOT,5)alpha,delt
        if(dokjar)then
        write(LOT,*)'Kjartansson Constant Q operator used'
        endif
        if(dosud)then
            write(LOT,*)'Wavefield separation at Source'
            write(LOT,*)'    Upward P-wave  ',spup
            write(LOT,*)'    Upward S-wave  ',ssup
            write(LOT,*)'    Downward P-wave',spdn
            write(LOT,*)'    Downward S-wave',ssdn
        endif
        if(dorud)then
            write(LOT,*)'Wavefield separation at Receiver'
            write(LOT,*)'    Upward P-wave  ',rpup
            write(LOT,*)'    Upward S-wave  ',rsup
            write(LOT,*)'    Downward P-wave',rpdn
            write(LOT,*)'    Downward S-wave',rsdn
        endif
        write(LOT,4)
c-----
c     open output file for hspec96
c-----
      open(unit=2,status='scratch',form=
     1            'unformatted',access='sequential')
      rewind 2
c-----
c       process the frequencies
c-----
        do  i=1,8
            ffreq(i)=-1.0
        enddo
        ndist = 1
        call bufini(1,ierr)
        do 100 ii = n1,n2
            freq=(ii-1)*df
            if(freq.lt.df) freq = 0.01*df
C       IF(II.EQ.N1.OR.II.EQ.N2)THEN
C           VERBY=.TRUE.
C           WRITE(6,*)'ii=',ii,freq
C       ELSE
C           VERBY=.FALSE.
C       ENDIF
            call excit(freq,xleng,xfac,dk,nk,omega)
            index=mod(ii,8)
            if(index.eq.0)index=8
            ffreq(index)=freq
            if (index.eq.8) then
                    write(LOT,3)ffreq
                    do 102 ij=1,8
                        ffreq(ij)=-1.
  102               continue
            endif
  100   continue
        call buflsh()
c-----
c       output the final spectrum as a function of distance
c-----
        open(unit=4,file='hspec96.grn',status='unknown',
     1      form='unformatted',access='sequential')
        rewind 4
c-----
c       iprog   1 hspec96
c           2 hspec96p
c           3 hwhole96
c-----
        iprog = 2
        write(4)iprog
        write(4) alpha,fl,fu,delt,n,n1,n2,df,nyq2
        write(4)mname
c-----
c       now output the spectrum for each distance
c
c       The order in the temporary file 'hspec91.tmp' is
c           FREQ
c               RAY PARAMETER
c                   SOURCE_DEPTH
c                       RECEIVER_DEPTH
c       This must be rearranged to form
c           RAY PARAMETER
c               SOURCE_DEPTH
c                   RECEIVER_DEPTH
c                       FREQ
c
c-----
c       because first arrival may be non-causal, shift time series
c       by 40*DT seconds
c-----
        tshft = -40.0*delt
c-----
c       TP TSV and TSH have no meaning in ordeinary sense
c-----
        TP = -12345.
        TSV = -12345.
        TSH = -12345.
c-----
c       output by ray parameter
c-----
        do 5000 jd=1,nk
        p  = (jd-1)*dp + pmin
        k = 0
        do 5005 js=1,mdpths
        do 5010 jr=1,mdpthr
            k = k + 1
            rewind 2
            call bufini(0,ierr)
c-----
c       safety do an insert
c----
            call srclay(depths(js),lmax,dph)
            VSB = b(lmax)
            VSA = a(lmax)
            VSR = rho(lmax)
c-----
c       define TI constants
c-----
            SR = VSR
            SL = VSR*VSB*VSB
            SN = VSR*VSB*VSB
            SA = VSR*VSA*VSA
            SC = VSR*VSA*VSA
            SF = SA - 2.0*SN
            
            write(4)p,tshft,depths(js)-refdep,
     1          depthr(jr)-refdep,
     2          TP, TSV, TSH, 
     3          SA, SC, SF, SL, SN, SR
            call srclay(depthr(jr),lmax,dph)
            VRB = b(lmax)
            VRA = a(lmax)
c-----
c       Dislocations and forces must act in a solid source region
c       if the receiver is in a fluid, 
c       then permit pressure field output
c-----
            do 5011 i=1,21
                ksrc(i) = jsrc(lsrc(i))
                if(ksrc(i).gt.0)then
                if(i.ge.1.and.i.le.8)then
                    if(VSB .le. 0.0001*VSA)then
                        ksrc(i) = 0
                    endif
                else if(i.ge.11.and.i.le.15)then
                    if(VSB .le. 0.0001*VSA)then
                        ksrc(i) = 0
                    endif
                elseif(i.eq.16)then
                    if(VRB .lt. 0.0001*VRA)then
                        ksrc(i) = 1
                    else
                        ksrc(i) = 0
                    endif
                else if(i.ge.17)then
                    if(VRB.lt.0.0001*VRA.and.VSB.gt.0.0)then
                        ksrc(i) = 1
                    else
                        ksrc(i) = 0
                    endif
                endif
                endif
 5011       continue
            write(4)ksrc
            do 5200 i=n1,n2
                freq=(i-1)*df
                fac = 6.2831853*freq*tshft
                ztmp = cmplx(cos(fac), sin(fac) )
                do 5300 jjd=1,nk
                kk = 0
                do 5301 jjs=1,mdpths
                do 5302 jjr=1,mdpthr
                    kk = kk + 1
                    do 5401 jj=1,21
                        if(jsrc(jj).eq.1)then
                            call bufrd(xr,ierr)
                            call bufrd(xi,ierr)
                            smo(jj)=cmplx(xr,xi)
                        endif
 5401               continue
                    do 5400 jj=1,21
                        if(jsrc(lsrc(jj)).eq.1)then
                        if(jjd.eq.jd .and. k.eq.kk.and.
     1                      ksrc(jj).ne.0)then
                            zdata=ztmp*smo(lsrc(jj))
                            datar= real(zdata)
                            datai=aimag(zdata)
                            write(4)datar,datai
                        endif
                        endif
 5400               continue
 5302           continue
 5301           continue
 5300           continue
 5200       continue
 5010   continue
 5005   continue
 5000   continue
        rr = -1.0
        tt0 = 0.0
        write(4)rr,tt0
        close (4)
        close(2)
 9999   continue
        end

        subroutine gethsp(lun,mname,fl,fu,delt,n,n1,n2,xleng,xfac,
     1      ndist,r,tshift,vred)
c-----
c       read in data file of hspec8 commands
c-----
        parameter(LER=0, LIN=5, LOT=6)
        parameter(NL=200)
        parameter(NSOURCE=100,NRECEIVER=100,NSR=100)
        common/source/depths(NSOURCE),lmaxs(NSOURCE),mdpths
        common/receiv/depthr(NRECEIVER),lmaxr(NRECEIVER),mdpthr
        real*4 depths, depthr
        integer lmaxs, lmaxr, mdpths, mdpthr
        common/bound/vmin,vamin,vamax,vbmin,vbmax
        common/damp/alpha,ieqex
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        common/jout/jsrc(21) , jbdrys, jbdryh
        parameter(NDST=500)
        real*4 r(NDST), tshift(NDST), vred(NDST)
        common/c/pmin,pmax,dp,pcntrl
        character ostr*80
        character mname*80, title*80
        common/lyrctl/lyrins
        logical lyrins
        logical ext

        common/cntrl/ishank,hnkarg,dstcor,dokjar
        integer*4 dstcor
        real*4 hnkarg
        logical ishank,dokjar

        common/depref/refdep
c-----
c       UNIX FORTRAN - NO CARRIAGE CONTROL
c-----
   21   format(11i5)
   22   format(10e10.3)
   24   format('XLENG=',e15.7,' XFAC=',e15.7)
   30   format(2e15.7/1x,3i10)
c-----
c       MSDOS MICROSOFT FORTRAN OUTPUT - CARRIAGE CONTROL
c-----
c   21  format(1x,11i5)
c   22  format(1x,10e10.3)
c   24  format(1x,'XLENG=',e15.7,' XFAC=',e15.7)
c   30  format(1x,2e15.7/1x,3i10)
c-----
c       read in distance correction
c-----
        read(lun,*)idcor
        if(idcor.ge.0 .and. idcor.le.2)then
            dstcor = idcor
        else
            dstcor = 0
        endif
c-----
c       read in time domain damping, sampling interval
c-----
        read(lun,*)alphat,delt
c-----
c       read in number of time samples, frequency limits
c-----
        read(lun,*)n,n1,n2  
        alpha = alphat/(n*delt)
        write(LOT,30) alpha,delt,n,n1,n2
        df = 1.0/(n*delt)
        fl = (n1-1)*df
        fu = (n2-1)*df
c-----
c       Specify desired output Green's functions
c-----
        read(lun,*)ieqex
        if(ieqex.gt.6 .or. ieqex.lt.0)ieqex = 2
        if(ieqex .eq. 0)then
            write(LOT,*)'ieqex= ',ieqex,' EARTHQUAKE + EXPLOSION'
        else if(ieqex.eq.1)then
            write(LOT,*)'ieqex= ',ieqex,' POINT FORCES + EXPLOSION'
        else if(ieqex.eq.2)then
            write(LOT,*)'ieqex= ',ieqex,' ALL GREEN'
        else if(ieqex.eq.3)then
            write(LOT,*)'ieqex= ',ieqex,' EXPLOSION ONLY'
        else if(ieqex.eq.4)then
            write(LOT,*)'ieqex= ',ieqex,' EARTHQUAKE ONLY'
        else if(ieqex.eq.5)then
            write(LOT,*)'ieqex= ',ieqex,' POINT FORCES ONLY'
        else if(ieqex.eq.6)then
            write(LOT,*)'ieqex= ',ieqex,' SH ONLY'
        endif
c-----
c       provide names for output Green's functions in order of output
c-----
        if(ieqex.eq.0)then
c-----
c           EARTHQUAKE + EXPLOSION
c-----
            do 1234 i=1,21
                if(i.le.8)then
                    jsrc(i) = 1
                else
                    jsrc(i) = 0
                endif
 1234       continue
            jsrc(13) = 1
            jsrc(14) = 1
            jsrc(16) = 1
            jsrc(17) = 1
            jsrc(18) = 1
            jsrc(19) = 1
            jsrc(20) = 1
            jsrc(21) = 1
        else if(ieqex.eq.1)then
c-----
c           POINT FORCES + EXPLOSION
c-----
            do 1235 i=1,21
                if(i.ge.7)then
                    jsrc(i) = 1
                else
                    jsrc(i) = 0
                endif
 1235       continue
            jsrc(13) = 0
            jsrc(14) = 0
            jsrc(17) = 0
            jsrc(18) = 0
            jsrc(19) = 0
        else if(ieqex.eq.2)then
c-----
c           EXPLOSION ONLY
c-----
            do 1236 i=1,21
                jsrc(i) = 1
 1236       continue
        else if(ieqex.eq.3)then
            do 1237 i=1,21
                if(i.eq.7 .or. i.eq.8 .or. i.eq.16)then
                    jsrc(i) = 1
                else
                    jsrc(i) = 0
                endif
 1237       continue
        else if(ieqex.eq.4)then
c-----
c           EARTHQUAKE ONLY
c-----
            do 1238 i=1,21
                if(i.le.6 .or. i.eq.13 .or. i.eq.14)then
                    jsrc(i) = 1
                else
                    jsrc(i) = 0
                endif
 1238       continue
        else if(ieqex.eq.5)then
c-----
c           POINT FORCES ONLY
c-----
            do 1239 i=1,21
                if(i.ge.9)then
                    jsrc(i) = 1
                else
                    jsrc(i) = 0
                endif
 1239       continue
            jsrc(13) = 0
            jsrc(14) = 0
            jsrc(16) = 0
        else if(ieqex.eq.6)then
c-----
c           SH ONLY
c-----
            do 1240 i=1,21
                if(i.eq.13 .or. i.eq.14 .or. i.eq.15)then
                    jsrc(i) = 1
                else
                    jsrc(i) = 0
                endif
 1240       continue
        endif
c-----
c       jsrc - array giving Green's functions to be evaluated
c           this controls computations, gives far field terms
c           or course true solution for radial may involve transverse
c       jsrc(lsrc) maps into Green,s functions, e.g.,
c           For REP=10, lsrc(10) = 8, and if jsrc(8) = 1
c           P-SV contribution to explosion radial 
c       time history is computed
c
c       ieqex = 0 
c       Earthquake + Explosion
c       1-ZDD   2-RDD   3-ZDS   4-RDS   5-TDS   6-ZSS
c       7-RSS   8-TSS   9-ZEP   10-REP
c
c       ieqex = 1 
c       Point Forces + Explosion
c       11-ZVF  12-RVF  13-ZHF  14-RHF  15-THF  16-PEP
c       9-ZEP   10-REP (others have no meaning)
c
c       ieqex = 2
c       All Green's functions
c       1-ZDD   2-RDD   3-ZDS   4-RDS   5-TDS   6-ZSS
c       7-RSS   8-TSS   9-ZEP   10-REP
c       11-ZVF  12-RVF  13-ZHF  14-RHF  15-THF  16-PEP
c
c       ieqex = 3
c       Explosion Only
c       9-ZEP   10-REP
c
c       ieqex = 4 
c       Earthquake
c       1-ZDD   2-RDD   3-ZDS   4-RDS   5-TDS   6-ZSS
c       7-RSS   8-TSS
c
c       ieqex = 5 
c       Point Forces Only
c       11-ZVF  12-RVF  13-ZHF  14-RHF  15-THF  
c
c       ieqex = 6
c       SH Solution Only
c       5-TDS 8-TSS 15-THF
c
c       If fluid layer for receiver, 16 is forced to be fluid 
c           stress due to explosion
c-----
c       input jbdry = 10*surface + halfspace
c       surface   = 0 - elastic halfspace   = 0 - elastic
c                   1 - free              1 - free 
c                   2 - rigid             2 - rigid
c-----
        read(lun,*)jbdry
        write(LOT,21) jbdry
        if(jbdry.lt.0)jbdry = 0
        ibdrys = jbdry / 10
        if(ibdrys.eq.0)then
            jbdrys = 0
        else if(ibdrys.eq.1)then
            jbdrys = 1
        else if(ibdrys.eq.2)then
            jbdrys = -1
        else
            jbdrys = 0
        endif
        ibdryh = mod(jbdry,10)
        if(ibdryh.eq.0)then
            jbdryh = 0
        else if(ibdryh.eq.1)then
            jbdryh = 1
        else if(ibdryh.eq.2)then
            jbdryh = -1
        else
            jbdryh = 0
        endif
c-----
c       jbdrys  =  surface boundary condition
c           = -1 top surface is rigid
c           =  0 really a halfspace with parameters of top layer    
c           =  1 free surface
c       jbdryh  = halfspace boundary condition
c           = -1  RIGID
c           =  0  ELASTIC
c           = +1  FREE SURFACE
c-----
        if(jbdrys.eq.1)then
            write(LOT,*)' TOP  OF MODEL IS FREE SURFACE  '
        else if(jbdrys.eq.0)then
            write(LOT,*)' TOP  OF MODEL IS HALFSPACE WITH',
     1          ' PROPERTIES OF FIRST LAYER'
        else if(jbdrys.eq.-1)then
            write(LOT,*)' TOP  OF MODEL IS RIGID'
        endif
        if(jbdryh.eq.0)then
            write(LOT,*)' BASE OF MODEL IS HALFSPACE WITH',
     1          ' PROPERTIES OF BOTTOM LAYER'
        else if(jbdryh.eq.-1)then
            write(LOT,*)' BASE OF MODEL IS RIGID'
        else if(jbdryh.eq.1)then
            write(LOT,*)' BASE OF MODEL IS FREE'
        endif


c-----
c       read in the earth model name
c-----
        read(lun,'(a)')mname
        lmnm = lgstr(mname)
        write(LOT,*)mname(1:lmnm)
c-----
c       read in the earth model in model96 format
c-----
        inquire(file=mname,exist=ext)
        if(.not. ext)then
            write(LER,*)'Model file not located'
            stop
        endif
        call getmod(1,mname,mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,.true.)
c-----
c       check the appropriateness of the model file
c-----
c               - -1 file does not exist
c               - -2 file is not a model file
c                 -3 error in the model file
        if(ierr.eq. -1)then
            call usage('Model file does not exist')
        else if(ierr.eq. -2)then
            call usage('Model file given is not a model96 file')
        else if(ierr.eq. -3)then
            call usage('Error in model file')
        endif
c-----
c       error checking
c-----
        if(idimen.ne.1)then
            call usage('1-D velocity model required')
        endif
        if(icnvel.ne.0)then
            call usage('Constant velocity model required')
        endif
        if(iiso.ne.0)then
            call usage('Isotropic velocity model required')
        endif
        if(iflsph.ne.0)then
            call usage('Flat earth velocity model required')
        endif
c-----
c       do not permit Q < 1 If qa or qb is entered > 1 
c       invert to form q inverse
c-----
        do 3007 i=1,mmax
            if(qa(i) .gt. 1.0)qa(i) = 1.0/qa(i)
            if(qb(i) .gt. 1.0)qb(i) = 1.0/qb(i)
            if(frefp(i) .le. 0.0)frefp(i) = 1.0
            if(frefs(i) .le. 0.0)frefs(i) = 1.0
 3007   continue
        call modcpy(.true.)
        call velbnd()
c-----
c       check model for inconsistencies
c-----
        call chkmod()
c-----
c       read in controls for wavenumber integration
c-----
        read(lun,*) xleng, xfac
        write(LOT,24)xleng, xfac
c-----
c       read in the  source depth
c-----
        read(lun,*)mdpths
        mtmp = NSOURCE
        if(mdpths .gt. NSOURCE)then
            write(ostr,3008)mdpths,mtmp
 3008   format(' NUMBER OF SOURCE DEPTHS',i5,' EXCEEDS DIMENSION ',i5)
            call werror(ostr)
        endif
c-----
c       get maximum source or receiver depth
c-----
        do 2008 i=1,mdpths
            read(lun,*)depths(i)
            depths(i) = depths(i) + refdep
 2008   continue
c-----
c       read in the receiver depth
c-----
        read(lun,*)mdpthr
        mtmp = NRECEIVER
        if(mdpthr .gt. NRECEIVER)then
            write(ostr,3009)mdpthr,mtmp
 3009   format(' NUMBER OF RECEIVER DEPTHS',i5,' EXCEEDS DIMENSION ',i5)
            call werror(ostr)
        endif
        do 2009 i=1,mdpthr
            read(lun,*)depthr(i)
            depthr(i) = depthr(i) + refdep
 2009   continue
c-----
c       check for filling the final depth array
c-----
        mtmp = mdpths * mdpthr
        ntmp = NSR
        if(mtmp .gt. ntmp)then
            write(ostr,3011)mtmp,ntmp
 3011   format(' NUMBER SOURCE-RECEIVER COMB',
     1      i5,' EXCEEDS DIMENSION ',i5)
            call werror(ostr)
        endif
c-----
c       read in the distances
c-----
        read(lun,*)ndist
        mtmp = NDST
        if(mdpthr .gt. NDST)then
            write(ostr,3010)ndist,mtmp
 3010   format(' NUMBER OF DISTANCES',i5,' EXCEEDS DIMENSION ',i5)
            call werror(ostr)
        endif
        do 2010 i=1,ndist
            read(lun,*)rr,tshf,vr
            r(i)=rr
            tshift(i) = tshf
            vred(i) = vr
 2010   continue
c-----
c       get phase velocity limits
c-----
        read(lun,*)pmin,pmax,dp,pcntrl
        write(LOT,*)'[pmin,pmax,dp,pcntrl]=[' ,pmin, ',' ,pmax, ',' ,dp,
     1      pcntrl, ']'
c-----
c       here pmin = mimimum ray parameter
c       pmax = maximum ray parameter
c       dp = ray parameter increment
c       pcntrl <= 0.0  modified for time series
c           > 0.0 true p-tau response
c-----
        if(pcntrl .le. 0.0)then
            write(LOT,*)' Modifed p-tau, for causal P arrivals'
        else
            write(LOT,*)' True p-tau'
        endif
        if(dp.le.0.0)dp = 1.0
c-----
c       For reasons of efficiency, decide whether to
c       add all layers at once, to the model
c       or to evaluate each layer source-receiver
c       combination separately. 
c
c       Roughly if the number of unique source and receiver depths
c       are mdpths+mdpthr if we insert layers, then we
c       end up with roughly mmax+mdpths+mdpthr layers, and
c       hence layer multiplication of this many matrices
c       for each source-receiver combination. Of course, for
c       equally spaced depth points, some economy arises
c       in avoiding matrix recomputation.
c
c       So if mdpths+mdpthr > 2*mmax we do not make a big model
c       other wise we do
c-----
c       adjust the model so that additional layers are added
c       to permit source and receiver at top of a give layer
c-----
c       lyrins = .true.
c       if(mdpths+mdpthr .gt. 2.0*mmax .and. lyrins.eq. .false.)then
        if(mdpths+mdpthr .gt. 2.0*mmax )then
            lyrins = .false.
            write(LOT,*)' LAYER INSERTION NOT DONE'
        else
            lyrins = .true. 
            write(LOT,*)' LAYER INSERTION DONE'
            do 2108 i=1,mdpths
                call insert(depths(i))
 2108       continue
            do 2109 i=1,mdpthr
                call insert(depthr(i))
 2109       continue
            call dezero()
c-----
c           check whether neighboring layers are identical
c           to avoid redundant evaluation
c----
            call equlck()
        endif
c-----
c       verify the new model parameters
c-----
        write(LOT,*)'mmax=',mmax
        write(LOT,22)(d(i),a(i),b(i),rho(i),
     1      qa(i),qb(i),etap(i),etas(i),frefp(i),frefs(i),i=1,mmax)
C-----
C removed 01 APR 2002 since messes up a mixed fluid/solid modium
Cc-----
Cc      Guarantee that no time wasted if any source in in the water
Cc      since there can only be a center of expansion source
Cc-----
C removed 01 APR 2002 since messes up a mixed fluid/solid modium
C       do 2019 i=1,mdpths
C           call srclay(depths(i), lmaxs(i), dphs)
C           if(b(lmaxs(i)).le.1.0e-04)then
C               do 2091 ii=1,21
C                   if(ii.ne.7 .and. ii.ne.8 
C     1                 .and. ii.ne.16)then
C                       jsrc(ii) = 0
C                   endif
C 2091          continue
C           endif
C 2019  continue
c-----
c       determine position of source and receive in new layers
c-----
        if(lyrins)then
            do 2020 i=1,mdpths
                call srclyr(depths(i), lmaxs(i), dphs)
 2020       continue
            do 2021 i=1,mdpthr
                call srclyr(depthr(i), lmaxr(i), dphs)
 2021       continue
        endif
        return
        end

        subroutine excit(freq,xleng,xfac,dk,nk,omega)
c-----
c     sample response for all wavenumbers at a given frequency
c     using Bouchon equal wavenumber sampling = dk
c     with offset of 0.218dk
c-----
        parameter(LER=0,LIN=5,LOT=6)
        parameter(NL=200)
        parameter(NSOURCE=100,NRECEIVER=100,NSR=100)
        common/source/depths(NSOURCE),lmaxs(NSOURCE),mdpths
        common/receiv/depthr(NRECEIVER),lmaxr(NRECEIVER),mdpthr
        real*4 depths, depthr
        integer lmaxs, lmaxr, mdpths, mdpthr
        common/bound/vmin,vamin,vamax,vbmin,vbmax
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        common/jout/jsrc(21) , jbdrys, jbdryh
        common/damp/alpha,ieqex
c-----
c       set up common blocks for wavenumber sampling at
c       suitable depths. This is necessary since the matrix
c       evaluation is done here for all source-receiver pairs
c       The source-receiver distance is important for the
c       wavenumber sampling at low frequencies
c-----
        common/kint1/gasymp
            logical gasymp(NSR)
        common/kint2/mkup
            integer mkup(NSR)
        common/kint3/wave
            real*4 wave(NSR,2)
        common/kint4/aa,bb,cc
            complex aa(NSR,21),bb(NSR,21),cc(NSR,21)

      common/result/gg
      complex gg(21)

        complex wvn,om, wvn2, om2
        complex zeye
c-----
c       matrix components in layers and boundaries saved
c-----
        common/water/iwater(NL),iwats(2),iwatb(2)
        common/lyrctl/lyrins
        logical lyrins
        common/c/pmin,pmax,dp,pcntrl

      common/srcrec/isrc,irec
      integer isrc,irec

      common/bcpsv/alppsv,betpsv
      complex alppsv(4,4), betpsv(4,4)

      common/bcsh/alpsh,betsh
      complex alpsh(2,2), betsh(2,2)

        common/jbdry/jtop,jbot
        integer jtop, jbot

        common/layertype/isfluid(NL),allsolid
        logical isfluid,allsolid

        common/shcontrol/ltop,lbot
        integer ltop, lbot

        logical proceed

        zeye = cmplx(0.0,1.0)
        omega=6.2831853*freq
        om=cmplx((omega),-(alpha))
        om2 = om * om
c-----
c       output by ray parameter
c-----
        nk = (pmax - pmin)/dp
        if(nk.lt.0)nk = 0
        nk = nk + 1
        do 3998 ii=1,nk
            p  = (ii-1)*dp + pmin
c-----
c       safety
c-----
            if(p.eq.0.0)p = 0.01*dp
            wvn=cmplx((p),0.0)*om
            wvn2 = wvn*wvn
            if(jbdrys.eq.0)then
                jtop = 3
            else if(jbdrys.eq.1)then
                jtop = 1
            else if(jbdrys.eq.-1)then
                jtop = 2
            endif
            if(jbdryh.eq.0)then
                jbot = 3
            else if(jbdryh.eq.1)then
                jbot = 1
            else if(jbdryh.eq.-1)then
                jbot = 2
            endif
c-----
c           
c           call dobotpsv(alppsv,jbot)
c           call dotoppsv(betpsv,jtop)
c           call botuppsv(betpsv,jtop)
c           call topdownpsv(betpsv,jtop)
c           call dobotpsh(alppsh,jbot)
c           call dotoppsh(betpsh,jtop)
c           call botuppsh(betpsh,jtop)
c           call topdownpsh(betpsh,jtop)
c-----
c       evaluate matrices first
c-----
C           if(lyrins)then
C               call evlmat(om,wvn,jbdrys,jbdryh,wvn2,om2)
C           endif
c-----
c       now evaluate for a specific source, receiver position
c-----
            k = 0
            do 4000 js=1,mdpths
                do 4010 jr=1,mdpthr
                    k = k + 1
                    if(.not.lyrins)then
c-----
c                   evaluate matrices first
c                   for currently defined layering
c-----
                       call modcpy(.false.)
                       call insert(depths(js))
                       call insert(depthr(jr))
                       call srclyr(depths(js), 
     1                     lmaxs(js), dphs)
                       call srclyr(depthr(jr), 
     1                     lmaxr(jr), dphr)
                       call dezero()
                       call equlck()
                    endif
c-----
c               initialize all GG to 0 + 0i
c-----
                   do i=1,21
                      gg(i) = cmplx(0.0,0.0)
                   enddo
c-----
c       determine if the layer is a fluid
c-----
                   do i=1,mmax
                      if(b(i).eq.0.0)then
                         isfluid(i) = .true.
                      else
                         isfluid(i) = .false.
                      endif
                   enddo
c-----
c       for SH there must be no fluid between the source
c       and receiver
c-----
                   allsolid = .true.
                   do i=min(isrc,irec),max(isrc,irec)
                      if(isfluid(i))then
                         allsolid = .false.
                      endif
                   enddo

                isrc = lmaxs(js)
                irec = lmaxr(jr)
c-----
c                        for SH processing all layers between source and
c                        receiver must be solid. We must also save the boundaries of 
c                        this solid region.   These will be know ltop, lbot
c
c                        first look between the source and receiver
c-----
                         allsolid = .true.
                         do i=min(isrc,irec),max(isrc,irec)
                           if(isfluid(i))then
                                allsolid = .false.
                           endif
                         enddo
C                        write(6,*)' Model can be processsed for SHd:',allsolid
                         if(allsolid)then
c-----
c                           get the bounds after setting defaults
c-----
                            ltop = 1
                            lbot = mmax 
                            proceed = .true.
                            do i=min(isrc,irec),1,-1
                               if(isfluid(i) .and. proceed)then
                                   ltop = i + 1
                                   proceed = .false.
                               endif
                            enddo
                            proceed = .true.
                            do i=max(isrc,irec),mmax
                               if(isfluid(i) .and. proceed)then
                                   lbot = i - 1
                                   proceed = .false.
                               endif
                            enddo
C                           write(6,*)'Effective model is in layers ',ltop,'-',lbot
                         endif
                call dopsv(OM,WVN,gg)
c               GG 01 to 12
                call dosh(om,wvn,gg)
c               GG 13 14 15
c-----
c       To make radial look pulse like for small ray parameters
c       Also take time derivative of the point force and 
c           pressure fields
c       The multiplication by 'i' accomplishes this for 
c           responses with a linear term in wavenumber
c       Note that these are the integrands and not the 
c           output Green s functions. Thus
c       the TDS contribution (output=5) is integrand (13)
c-----
                    if(pcntrl .le.0.0)then
                        GG(2)  = - GG(2) * zeye
                        GG(3)  =   GG(3) * zeye
                        GG(5)  = - GG(5)
                        GG(6)  =   GG(6) * zeye
                        GG(14) =   GG(14) * zeye
                        GG(8)  = - GG(8) * zeye
                        GG(9)  =   GG(9) * zeye * om
                        GG(10) = - GG(10) * zeye * om * zeye
                        GG(11) =   GG(11) * zeye * om * zeye
                        GG(12) =   GG(12) * zeye * om 
                        GG(15) =   GG(15) * zeye * om 
                    endif
c-----
c       output
c-----
                    do 3999 j=1,21
                        if(jsrc(j).eq.1)then
                            call bufwr(real(gg(j)))
                            call bufwr(aimag(gg(j)))
                        endif
 3999               continue
 4010           continue
 4000       continue
 3998   continue
        return
        end

        subroutine gcmdln()
c-----
c       parse command line arguments
c       requires subroutine mgtarg() and function mnmarg()
c-----
        common/cntrl/ishank,hnkarg,dstcor,dokjar
        integer*4 dstcor
        real*4 hnkarg
        logical ishank,dokjar
        common/ctlrud/rpup,rpdn,rsup,rsdn,dorud
        logical rpup, rpdn, rsup, rsdn, dorud
        common/ctlsud/spup,spdn,ssup,ssdn,dosud
        logical spup, spdn, ssup, ssdn, dosud

        integer*4 mnmarg
        character*50 name
        dokjar = .false.
        dosud = .false.
        spup = .false.
        spdn = .false.
        ssup = .false.
        ssdn = .false.
        dorud = .false.
        rpup = .false.
        rpdn = .false.
        rsup = .false.
        rsdn = .false.
        nmarg=mnmarg()
        i = 0
   11   i = i + 1
            if(i.gt.nmarg)goto 13
            call mgtarg(i,name)
            if(name(1:3).eq.'-SU')then
                spup = .true.
                ssup = .true.
                ssdn = .false.
                spdn = .false.
                dosud = .true.
            else if(name(1:3).eq.'-SD')then
                spup = .false.
                ssup = .false.
                ssdn = .true.
                spdn = .true.
                dosud = .true.
            else if(name(1:5).eq.'-SPUP')then
                spup = .true.
                dosud = .true.
            else if(name(1:5).eq.'-SSUP')then
                ssup = .true.
                dosud = .true.
            else if(name(1:5).eq.'-SPDN')then
                spdn = .true.
                dosud = .true.
            else if(name(1:5).eq.'-SSDN')then
                ssdn = .true.
                dosud = .true.
            else if(name(1:5).eq.'-RPUP')then
                rpup = .true.
                dorud = .true.
            else if(name(1:5).eq.'-RSUP')then
                rsup = .true.
                dorud = .true.
            else if(name(1:5).eq.'-RPDN')then
                rpdn = .true.
                dorud = .true.
            else if(name(1:5).eq.'-RSDN')then
                rsdn = .true.
                dorud = .true.
            else if(name(1:3).eq.'-RD')then
                rpup = .false.
                rsup = .false.
                rsdn = .true.
                rpdn = .true.
                dorud = .true.
            else if(name(1:3).eq.'-RU')then
                rpup = .true.
                rsup = .true.
                rsdn = .false.
                rpdn = .false.
                dorud = .true.
            else if(name(1:2).eq.'-K')then
                dokjar = .true.
            else if(name(1:2).eq.'-?')then
                call usage(' ')
            else if(name(1:2).eq.'-h')then
                call usage(' ')
            endif
        goto 11
   13   continue
c-----
c       safety check
c-----
        return
        end

        subroutine usage(ostr)
c------
c       write out program syntax
c-----
        character ostr*(*)
        integer LER,LOT,LIN
        parameter(LER=0, LIN=5, LOT=6)
        if(ostr.ne. ' ' )then
            lostr = lgstr(ostr)
            write(LER,*)ostr(1:lostr)
        endif
        write(LER,*)'USAGE: ',
     1  'hspec96p [-K] ',
     2      '[-SU] [-SD] [-SPUP] [-SSUP] [-SPDN] [-SSDN] ',
     3      '[-RU] [-RD] [-RPUP] [-RSUP] [-RPDN] [-RSDN] [-?] [-h]'
        write(LER,*)
     1  '-K      (default Futterman) use Kjartansson Causal Q'
        write(LER,*)
     1  'The following govern wavefield at source. The default is',
     2  ' the entire wavefield'
        write(LER,*)
     1  '-SU      (default whole wavefield) Compute only upgoing ',
     2  '               wavefield from the source'
        write(LER,*)
     1  '-SD      (default whole wavefield) Compute only downgoing ',
     2  '               wavefield from the source'
        write(LER,*)
     1  ' -SPUP  Include upward P at source'
        write(LER,*)
     1  ' -SSUP  Include upward S at source'
        write(LER,*)
     1  ' -SPDN  Include downward P at source'
        write(LER,*)
     1  ' -SSDN  Include downward S at source'
        write(LER,*)
     1  'The following govern wavefield at receiver. The default is',
     2  ' the entire wavefield'
        write(LER,*)
     1  ' -RD    Only downgoing waves at receiver'
        write(LER,*)
     1  ' -RU    Only upgoing waves at receiver'
        write(LER,*)
     1  ' -RPUP  Include upward P at receiver'
        write(LER,*)
     1  ' -RSUP  Include upward S at receiver'
        write(LER,*)
     1  ' -RPDN  Include downward P at receiver'
        write(LER,*)
     1  ' -RSDN  Include downward S at receiver'
        write(LER,*)
     1  '-?                   Display this usage message'
        write(LER,*)
     1  '-h                   Display this usage message'
        stop 
        end

        subroutine equlck()
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        common/updnsm/equalu(NL), equald(NL)
        logical equalu, equald
c-----
c       To avoid repeated computation, check to see if 
c       neighboring layers
c       are identical, once for going up and another for going down
c-----
c       First check top down
c-----
        do 100 m=1,mmax
            if(m.eq.1)then
                equald(m) = .false.
            else if(m.gt.1
     1          .and. a(m).eq.a(m-1) 
     2          .and. b(m).eq.b(m-1)
     3          .and. d(m).eq.d(m-1) 
     4          .and. rho(m).eq.rho(m-1)
     5          .and. qa(m).eq.qa(m-1)
     6          .and. qb(m).eq.qb(m-1) )then
                equald(m) = .true.
            else
                equald(m) = .false.
            endif
  100   continue
c-----
c       check bottom up
c-----
        do 200 m=1,mmax
            if(m.eq.mmax)then
                equalu(m) = .false.
            else if(m.lt.mmax
     1          .and. a(m).eq.a(m+1) 
     2          .and. b(m).eq.b(m+1)
     3          .and. d(m).eq.d(m+1) 
     4          .and. rho(m).eq.rho(m+1)
     5          .and. qa(m).eq.qa(m+1)
     6          .and. qb(m).eq.qb(m+1) )then
                equalu(m) = .true.
            else
                equalu(m) = .false.
            endif
  200   continue
        return
        end

        subroutine insert(dph)
        parameter (LER=0, LIN=5, LOT=6)
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
c-----
c       Insert a depth point into the model by splitting a layer
c       so that the point appears at the top boundary of the layer
c       dph = depth of interest
c-----
c       determine the layer in which the depth dph lies.
c       if necessary, adjust  layer thickness at the base
c-----
c       Here determine layer corresponding to specific depth dph
c       If the bottom layer is not thick enough, extend it
c
c       dep - depth to bottom of layer
c       dphh    - height of specific depth above bottom of the layer
c-----
        if(dph.le.0)then
            d(1) = d(1) - dph
            return
        else if(dph.ge.0)then
            dep = 0.0 
            dp = 0.0 
            dphh = -1.0
            do 100 m = 1,mmax 
                dp = dp + d(m) 
                dphh = dp - dph 
                if(m.eq.mmax)then
                    if(d(mmax).le.0.0 .or. dphh.lt.0.0)then
                        d(mmax) = (dph - dp)
                    endif
                endif
                dep = dep + d(m) 
                dphh = dep - dph 
                ls = m 
                if(dphh.ge.0.0) go to 101 
  100       continue 
  101       continue 
        endif
c-----
c       In the current model, the depth point is in the ls layer
c       with a distance dphh to the bottom of the layer
c
c       Do not create unnecessary layers, e.g., 
c       at surface and internally
c       However do put in a zero thickness layer 
c       at the base if necessary
c-----
        if(dphh .eq. 0.0 .and. ls.ne.mmax)then
            return
        else
c-----
c           adjust layering
c-----
             do 102 m = mmax,ls,-1
                d(m+1) = d(m)
                a(m+1) = a(m)
                b(m+1) = b(m)
                rho(m+1) = rho(m)
                qa(m+1) = qa(m)
                qb(m+1) = qb(m)
                etap(m+1) = etap(m)
                etas(m+1) = etas(m)
                frefp(m+1) = frefp(m)
                frefs(m+1) = frefs(m)
  102       continue
            hsave = d(ls)
            d(ls) = hsave - dphh
            d(ls+1) = dphh
            ls = ls + 1
            mmax = mmax + 1
            if(d(mmax).lt.0.0)d(mmax)=0.0
        endif
        return
        end

        subroutine dezero()
        parameter (LER=0, LIN=5, LOT=6)
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
c-----
c       ultimately get rid of zero thickness layers - this
c       will require readjusting the model from top down, and
c       also readjusting the source and receiver indices.
c----
c       Here just guarantee that the halfspace is not of zero thickness
c-----
        dmin = 1.0e+30
        do 100 i=1,mmax-1
            if(d(i) .lt. dmin .and. d(i).gt.0.0)dmin = d(i)
  100   continue
c       if(d(mmax).le.0.0)then
c           d(mmax) = 0.1*dmin
c       endif
        return
        end

        subroutine srclyr(depth,lmax,dph)
        parameter (LER=0, LIN=5, LOT=6)
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
c-----
c       Find source/receiver boundary. It is assumed that
c       it will lie upon a boundary
c
c       lmax = source layer 
c            = 0 is the free surface 
c       depth = source depth 
c       dph = height of  source above lmax + 1 interface 
c       lmax = 0 is the free surface 
c-----
        if(depth.le.0.0)then
            lmax = 1
            dph = 0.0
        else
            dep = 0.0 
            do 100 m = 2,mmax
                dep = dep + d(m-1) 
                dph = dep - depth 
                lmax = m 
                if(abs(dph).lt. 0.0001*d(m-1) .or.
     1              abs(dph).lt.1.0e-6)go to 101
  100       continue 
  101   continue 
        endif
        return 
        end 

        subroutine srclay(depth,lmax,dph)
        parameter (LER=0, LIN=5, LOT=6)
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        common/lyrctl/lyrins
        logical lyrins
        if(.not.lyrins)then
            call modcpy(.false.)
            call insert(depth)
        endif
        call srclyr(depth,lmax,dph)
        return
        end

        subroutine velbnd() 
        parameter (LER=0, LIN=5, LOT=6)
c-----
c       get bounds on earth model 
c-----
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        common/bound/vmin,vamin,vamax,vbmin,vbmax
c-----
c       initialize bound search
c-----
        vamin = 1.0e+38
        vbmin = 1.0e+38
        vmin  = 1.0e+38
        vamax = 0.0
        vbmax = 0.0
        write(LOT,2) 
    2   format(' ',7x,'d',9x,'a',9x,'b',9x,'rho',6x,'1/qa',6x,'1/qb')
    3   format(' ',4f10.3,2f10.6) 
        do 20 i = 1,mmax 
            if(a(i).gt.vamax)vamax=a(i)
            if(b(i).gt.vbmax)vbmax=b(i)
            if(a(i).lt.vamin)vamin=a(i)
            if(b(i).lt.vbmin .and. b(i).gt.0.0)vbmin=b(i)
            if(b(i).gt.0.1)then
                if(b(i).lt.vmin)vmin=b(i)
            else
                if(a(i).lt.vmin)vmin=a(i)
            endif
            if(i.lt.mmax)then
            write(LOT,3)d(i),a(i),b(i),rho(i),qa(i),qb(i)
            endif
   20   continue 
    5   format(' ',10x,3f10.3,2f10.6/' ') 
        write(LOT,5)a(mmax),b(mmax),rho(mmax),qa(mmax),qb(mmax) 
c-----
c     obtain extreme velocity limits
c-----
      return 
      end 

        subroutine bufini(irdwr,ierr)
c-----
c       initialize buffer pointer
c       irdwr = 0 read initialize
c       irdwr = 1 write initialize
c-----
        integer BUFMAX
        parameter(BUFMAX=2000)
        common/buf/iptr,max,buffer(BUFMAX)
        iptr = 1
        if(irdwr.eq.0)call getbuf(ierr)
        return
        end

        subroutine buflsh
c-----
c       flush output buffer
c-----
        integer BUFMAX
        parameter(BUFMAX=2000)
        common/buf/iptr,max,buffer(BUFMAX)
        ipt = iptr -1
        if(ipt.gt.0)write(2)ipt,(buffer(i),i=1,ipt)
        iptr = 1
        return
        end

        subroutine bufwr(x)
c-----
c       fill buffer with floating point variable x,
c       flush buffer as necessary
c-----
        integer BUFMAX
        parameter(BUFMAX=2000)
        common/buf/iptr,max,buffer(BUFMAX)
        buffer(iptr) = x
        iptr = iptr + 1
        if(iptr.gt.BUFMAX)call buflsh
        return
        end

        subroutine getbuf(ierr)
c-----
c       read in file contents into buffer, taking care not to
c       read beyond the contents of the file
c-----
        integer BUFMAX
        parameter(BUFMAX=2000)
        common/buf/iptr,max,buffer(BUFMAX)
c-----
c       ierr = 0 successful read
c            = 1 read error
c            = 2 end of file
c-----
        read(2,err=1000,end=2000)max,(buffer(i),i=1,max)
        iptr = 1
        ierr = 0
        return
 1000   ierr = 1
        return
 2000   ierr = 2
        return
        end

        subroutine bufrd(x,ierr)
c-----
c       retrieve a value from buffer array, red in new array
c       as necessary
c       iptr is here the next array element to be read
c       it is always >= 1. We do not worry the upper limit
c       since the calling program must worry about this
c       because read always follows a complete write
c-----
        integer BUFMAX
        parameter(BUFMAX=2000)
        common/buf/iptr,max,buffer(BUFMAX)
c-----
c       only yank in new data if actually required
c-----
        if(iptr.gt.max)call getbuf(ierr)
        x = buffer(iptr)
        iptr = iptr + 1
        return
        end

        subroutine werror(ostr)
c-----
c       output error message and terminate program
c-----
        parameter(LER=0, LIN=5, LOT=6)
        character ostr*(*)
        write(LER,*)'PROGRAM TERMINATION'
        write(LER,*)ostr
        stop
        end

        subroutine modcpy(totmp) 
        logical totmp
c-----
c       copy model to temporary array
c-----
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        common/modelt/dt(NL),at(NL),bt(NL),rhot(NL),
     1      mmaxt,qat(NL),qbt(NL),etapt(NL),etast(NL),
     1      frefpt(NL),frefst(NL)
c-----
c       copy to temporary array
c-----
        if(totmp)then
            do 20 i = 1,mmax 
                dt(i) = d(i)
                at(i) = a(i)
                bt(i) = b(i)
                rhot(i) = rho(i)
                qat(i) = qa(i)
                qbt(i) = qb(i)
                etapt(i) = etap(i)
                etast(i) = etas(i)
                frefpt(i) = frefp(i)
                frefst(i) = frefs(i)
   20       continue 
            mmaxt = mmax
        else
            do 30 i = 1,mmaxt 
                d(i) = dt(i)
                a(i) = at(i)
                b(i) = bt(i)
                rho(i) = rhot(i)
                qa(i) = qat(i)
                qb(i) = qbt(i)
                etap(i) = etapt(i)
                etas(i) = etast(i)
                frefp(i) = frefpt(i)
                frefs(i) = frefst(i)
   30       continue 
            mmax = mmaxt
        endif
        return 
        end 

        subroutine chkmod()
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        common/jout/jsrc(21) , jbdrys, jbdryh
        common/lwater/lfluid
        logical lfluid
c-----
c       check model for inconsistencies
c-----
c-----
c       Model cannot consist entirely of water layers
c       Also determine first solid layer from top
c-----
        iw = 0  
        isoldd = 0
        do 100 i=1,mmax
            if(b(i).eq.0.0)then
                iw = iw + 1
            else
                if(isoldd .eq.0)isoldd=i
            endif
  100   continue
        if(iw .eq. mmax)then
            lfluid = .true.
C           call werror('MODEL CONSISTS ONLY OF LIQUID LAYERS')
        else
            lfluid = .false.
        endif
c-----
c       Determine first solid layer from bottom
c-----
        iw = 0  
        isoldu = 0
        do 101 i=mmax,1,-1
            if(b(i).eq.0.0)then
                iw = iw + 1
            else
                if(isoldu .eq.0)isoldu=i
            endif
  101   continue
c-----
c       Check for interior water layer
c-----
C       if(iw.gt.0 .and. .not. lfluid)then
C           do 102 i=isoldd,isoldu
C               if(b(i).eq.0.0)then
C               call werror('MODEL HAS INTERIOR  FLUID LAYERS')
C               endif
C 102       continue
C       endif
c-----
c       If boundary condition is rigid, and the adjacent layer is
c       fluid, terminate 
c-----
C       if(b(1).le.1.0e-04 .and. jbdrys.eq.-1 .and. lfluid)then
C           call werror('TOP LAYER IS FLUID AND RIGID BOUNDARY')
C       endif
C       if(b(mmax).le.1.0e-04 .and. jbdryh.eq.-1 .and. lfluid)then
C           call werror('BOTTOM LAYER IS FLUID AND RIGID BOUNDARY')
C       endif
        return
        end
