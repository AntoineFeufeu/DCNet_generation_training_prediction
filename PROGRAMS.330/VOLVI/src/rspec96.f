        program rspec96
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME VI                                                      c
c                                                                     c
c      PROGRAM: rSPEC96                                               c
c                                                                     c
c      COPYRIGHT                                                  c
c      R. B. Herrmann                                                 c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63013                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
c       CHANGES TO HSPEC96
c       11 SEP 2000 - build in P, SV and SH first arrival times
c       14 JAN 2001 - put in A, C, F, L and N constants into
c           trace header
c       18 MAR 2001 - new formalism from Herrmann (2001) 
c           Seismic Waves in Layered Media
c       18 JUL 2001 - fixed specification for KSRC when source and
c           receiver are in fluids. 
c       09 APR 2002 - another hack at fixing the 
c           JSRC when fluids are involved.
c           commented out lines in gethsp which was not correct, 
c           modified output loop 
c       29 OCT 2002 - Fixed evalg and evalh for fluid layers
c       18 APR 2005 - put in simple minded earth flattening 
c           for J. Ritsema
c       16 JUN 2005 - travel time for P and S work for spherical
c       04 AUG 2006 - corrected error in first arrival pick that
c               falsely gave the refraction time instead of the
c               direct time because the refraction arrival was
c               unphysical
c       15 DEC 2007 - introduced an earth flattening correction
c               that uses a different density mapping for the SH and
c               P-SV problems 
c               The sphericity correction invoked for a spherical model
c               is as follows
c               
c               Spherical              Flat
c                   r                  x=ln(a/r)   where a is radius of sphere
c               velocity               velocity*(a/r)
c               density                density*(r/a)**POWER
c                                        POWER=5 for SH
c                                        POWER=2.275 for P-SV
c       14 JAN 2008 - put Radius of Earth into common/earth/radius for
c               generality
c               correct layer insertion for spherical earth
c               NOTE refdep handled correctly for a spherical earth
c       22 JAN 2008 - rearranged routines, put in new common
c                     isotropic travel time routines from time96
c       25 JAN 2008 - put Radius of Earth into common/earth/radius for
c                      generality
c                   -  define a separate common block for the
c                      SH velocity and density
c                   -  have sphericity correction
c                      work on common blocks instead of procedure call
c                   -  create a default adomod  to fill the SH for a flat model
c                      note the separation of SH is important for wavenumber
c                      integration code
c                      NOTE - the separation into qsh bsh is incomplete 
c                        since the current mapping requires only a 
c                        different density mapping for
c                        SH - we would have to create a varsh 
c                        routine to be consistent,
c                        but that would increase computational time slightly
c       08 FEB 2008 -  subtle change in fstarr for source receiver in same layer -
c                      spherical mapping was not done
c       18 FEB 2009 -  took care to add common/depref/refdep in frstar
c       24 MAR 2012 -  corrected code to output pressure field in a fluid
c       01 JUN 2013 -  Modified subroutine to prevent NaN for direct ray by replacing
c                      pnew = 0.999d+00 * pupper to
c                      pnew = 0.99d+00 * pupper
c                      also modified subroutine frstar to have the dogeom argument.
c                      We do not want to compute teleseismic geometrical spreading
c       09 OCT 2013 -  Corrected code at line 452 - this C syntax works with
c                      gfortran but not with Solaris f77 - sorry
c       16 OCT 2018 -  used double precision to specify source/receiver depth in
c                      flattened model
c
c       THIS IS BUILT ON hspec96 BUT USES THE GENERALIZED REFLECTION
c       TRANSMISSION MATRICES INSTEAD OF PROPAGATOR MATRICES. FOR TH#E
c       P-SV PROBLEM THIS MEANS THAT IT CAN HANDLE ARBITRARY SEQUENCES
c       OF FLUID AND SOLID LAYERS. THE SH PART WILL GIVE RESULTS IF THERE
c       ARE NO INTERVENING FLUID LAYERS BETWEEN THE SOURCE AND RECEIVER
c       21 AUG 2019 - implemented. Note that there is some redundancy
c                     in the coding so that to permit easy merging
c                     into the old code
c       26 DEC 2019 - implemented the stress field PDD PDS PDD PEX PVF PHF
c                     when the receiver is in a fluid. The PDD PDS PSS PVF and PHF
c                     only exist if the source is in a solid
c       25 JAN 2022 -  if source is in fluid do not compute S time
c       27 JAN 2022 -  Changed screen output to handle high frequencies
c                      and short distances
c-----
        parameter(LER=0, LIN=5, LOT=6)
        parameter(NL=200)
        parameter(NSOURCE=100,NRECEIVER=100,NSR=100)
        common/source/depths(NSOURCE),lmaxs(NSOURCE),mdpths
        common/receiv/depthr(NRECEIVER),lmaxr(NRECEIVER),mdpthr
        common/sourcesv/depthssv(NSOURCE)
        common/receivsv/depthrsv(NRECEIVER)
        real*4 depths, depthr
        integer lmaxs, lmaxr, mdpths, mdpthr
        common/bound/vmin,vamin,vamax,vbmin,vbmax
        common/damp/alpha,ieqex
        common/modctl/iflsph
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
        real bsh, rhosh, qbsh
        common/modlly/mmax
        common/modelt/dt(NL),at(NL),bt(NL),rhot(NL),
     1      mmaxt,qat(NL),qbt(NL),etapt(NL),etast(NL),
     1      frefpt(NL),frefst(NL),bsht(NL),rhosht(NL),qbsht(NL)
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
        character*3 istat3
        logical ixst3
        real*4 ffreq(8)
        parameter(NDST=500)
        real*4 r(NDST), tshift(NDST), vred(NDST)
        complex smm(NSR,21)
        complex ztmp, zdata
        common/c/cmax,c1,c2,cmin
        common/frlim/fl,fu,df,fwhich
        common/cntrl/ishank,hnkarg,dstcor,dokjar,docausal

c-----
c       ishank  L   - .true. use Hankel function and not Bessel
c                    note that asymptotic tricks are not used
c       hnkarg  R*4 - (default 6.0) For kr > hnkarg use the Hankel
c                   function, otherwise the Bessel
c       isfref  L   - .true. reference frequency uis not the 
c                   default of 1.0 Hz
c       fref    R*4 - reference frequency for Causal Q
c       dstcor  I*4 - 0 first time sample is tshift + r/vred
c                 1 first time sample is tshift + z/vred where
c                   z = abs (source depth - receiver depth)
c                 2 first time sample is tshift + 
c                   sqrt(z*z + r*r) /vred
c       dokjar  L   - .false. use Futterman causal Q
c       docausal L  - .true. Use causal Q
c                 .true.  use Kjartansson Causal Q
c-----
        integer*4 dstcor
        real*4 hnkarg
        logical ishank,dokjar,docausal
        common/ctlrud/rpup,rpdn,rsup,rsdn,dorud
        logical rpup, rpdn, rsup, rsdn, dorud
        common/ctlsud/spup,spdn,ssup,ssdn,dosud
        logical spup, spdn, ssup, ssdn, dosud

        logical ext
        character mname*80, title*80

        common/lwater/lfluid
        logical lfluid

        common/depref/refdep
        real refdep

        common/earth/radius
        real radius

        real tp, tsv, tsh, vra, vrb, vsa, vsb, vsr, rayp, geom, tstar
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
c       initialize
c-----
        radius = 6371.
c-----
c       parse command line arguments
c-----
        call gcmdln()
        if(ishank)write(LOT,2042)hnkarg
c-----
c       set up tolerances
c           rlim is the distance which is effectively zero
c-----
            rlim = 1.0e-07
c-----
c       See if the file hspec96.dat exists, if it does
c       open it for all control information
c-----
        inquire(file='hspec96.dat',exist=ext)
        if(.not. ext)then
                write(LER,*)'Control file ',
     1          'hspec96.dat does not exist'
                go to 9999
        endif
        open(3,file='hspec96.dat',access='sequential',
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
     2      vmin,vamin,vamax
        if(.not.lfluid)write(LOT,21)vbmin,vbmax
c-----
c       UNIX output - no carriage control
c-----
    2   format('fl =',1pg12.5,5x,'fu =',1pg12.5,5x,'df =',1pg12.5,/
     1      4x,'n1 =',i4,5x,'n2 =',i4,5x, ' n =',i5/
     2  'vmin  =',1pg12.5,' vamin =',1pg12.5,' vamax =',1pg12.5)
   21   format(
     1  '       ',10x  ,' vbmin =',1pg12.5,' vbmax =',1pg12.5)
 2021   format('depths =',f20.6,' (',f20.6,')')
 2031   format('depthr =',f14.6,' (',f20.6,')')
 2041   format('     r =',f14.6,' tshift=',f10.2,' vred=',f10.2)
    3   format(8(1pg12.5))
    4   format('frequencies for which response computed     ')
    5   format('alpha =',1pg15.7,5x,'dt =',1pg15.3)
 2042   format('Hankel function used for kr >=',f10.2)

        write(LOT,*)'SOURCE DEPTH IN WORKING AND ORIGINAL MODEL (',
     1       mdpths,')'
        do 2020 i=1,mdpths
            write(LOT,2021)depths(i), depthssv(i) - refdep
 2020   continue
        write(LOT,*)'RECEIVER DEPTH IN WORKING AND ORIGINAL MODEL ('
     1       ,mdpthr,')'
        do 2030 i=1,mdpthr
            write(LOT,2031)depthr(i), depthrsv(i) - refdep
 2030   continue
        write(LOT,*)'RECEIVER DISTANCES (',ndist,')'
        do 2040 i=1,ndist
            write(LOT,2041)r(i), tshift(i), vred(i)
 2040   continue
        write(LOT,5)alpha,delt
        if(docausal)then
            if(dokjar)then
            write(LOT,*)'Kjartansson Constant Q operator used'
            endif
        else
            write(LOT,*)'Non-causal Q operator used'
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
      open(unit=2,file='hspec96.tmp',status='unknown',form=
     1            'unformatted',access='sequential')
      rewind 2
c-----
c     open temporary distance multiplexed file which
c     will be examined in case it exists already for
c     work already done
c-----
      inquire(file='hspec96.bin',exist=ixst3)
      if(ixst3)then
            istat3 = 'old'
      else
            istat3 = 'new'
      endif
      open(unit=3,file='hspec96.bin',status=istat3,form=
     1            'unformatted',access='sequential')
      rewind 3
      ifreq = 0
      if(ixst3)then
            call reset3(ifreq,jsrc,lsrc,ndist,n1,n2)
            write(LER,*)'Partial processing - starting at',
     1          ifreq, ' in range [',n1,'-' ,n2,']'
      endif
c-----
c       process the frequencies
c-----
        do 101 i=1,8
            ffreq(i)=-1.0
  101   continue
        n11 = n1
        if(ifreq.gt.n1)n11 = ifreq + 1
        do 100 ii = n11,n2
            freq=(ii-1)*df
            if(freq.lt.df) freq = 0.01*df
            call excit(freq,xleng,xfac,dk,nk)
            do 200 jd=1,ndist
                call setup(r(jd))
                call wvint(r(jd),smm,dk,nk)
                k = 0
                do 303 is=1,mdpths
                    do 304 ir=1,mdpthr
                        call gett0(t0,r(jd),depthssv(is),
     1                      depthrsv(ir),tshift(jd),
     2                      vred(jd),dstcor)
                        fac = 6.2831853*freq*t0
                        ztmp = cmplx(cos(fac), sin(fac) )
                        k = k + 1
                        do 301 jj=1,21
                            zdata = ztmp*smm(k,jj)
                            if(jsrc(lsrc(jj)).eq.1)then
                            datar= real(zdata)
                            datai=aimag(zdata)
                            write(3)datar,datai
                            endif
  301                   continue
  304               continue
  303           continue
  200       continue
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
        write(LOT,3)ffreq
        write(LOT,*)'Transition to non-asymptotic at',fwhich
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
        iprog = 1
        write(4)iprog
        write(4) alpha,fl,fu,delt,n,n1,n2,df,nyq2
        write(4)mname
c-----
c       now output the spectrum for each distance
c
c       The order in the temporary file 'hspec96.bin' is
c           FREQ
c               DIST
c                   SOURCE_DEPTH
c                       RECEIVER_DEPTH
c       This must be rearranged to form
c           DIST
c               SOURCE_DEPTH
c                   RECEIVER_DEPTH
c                       FREQ
c
c-----
        write(LOT,'(a,a)')
     1  '     r(km)   t0(sec)    hs(km)    hr(km)    Tp(sec)   Tsv(sec)'
     2  ,'   Tsh(sec)'
        do 5000 jd=1,ndist
        if (iflsph.eq.0)then
                gfac = 1.0
        else
                if(r(jd).eq.0.0)then
                        gfac = 1.0
                else
                        gfac = (r(jd)/radius)/sin(r(jd)/radius)
                        gfac = sqrt(gfac)
                endif
        endif
        k = 0
        do 5005 js=1,mdpths
        do 5010 jr=1,mdpthr
            k = k + 1
            rewind 3
            call gett0(t0,r(jd),depthssv(js),
     1          depthrsv(jr),tshift(jd),vred(jd),dstcor)
c-----
c       to get first arrivals, we must account for the fact that
c       there may be a a difference reference depth which is
c       indicated by the consistent use of negative layer thicknesses
c       in the upper part of the model. The getmod() subroutine
c       notes this and converts to a model that has a surface reference
c       depth. This program adds this reference depth to all
c       source ande receiver depths input. This in calling FRSTAR,
c       which rereads the model file, we use the adjusted source
c       depths. When we output the results, 
c          we convert back to the actual
c       physical depths
c-----
            CALL FRSTAR(R(JD),DEPTHSSV(JS)-refdep,DEPTHRSV(JR)-refdep,
     1          MNAME,1,TP ,VRA, VRB, den,
     1          VSA, VSB, VSR, rayp, geom, tstar, .false. ,.false.)
            CALL FRSTAR(R(JD),DEPTHSSV(JS)-refdep,DEPTHRSV(JR)-refdep,
     1          MNAME,2,TSV,VRA, VRB, den,
     1          VSA, VSB, VSR, rayp, geom, tstar, .false. ,.false.)
            CALL FRSTAR(R(JD),DEPTHSSV(JS)-refdep,DEPTHRSV(JR)-refdep,
     1          MNAME,3,TSH,VRA, VRB, den,
     1          VSA, VSB, VSR, rayp, geom, tstar, .false. ,.false.)
c-----
c      Fix 24 MARCH 2012 - when there is a fluid, the TSH will be
c      a large number. Use the Sac not defined as the value
c-----
            if(TSH .gt. 1.0e+29)TSH=-12345.
            if(TSV .gt. 1.0e+29)TSV=-12345.
            if(VSB .eq. 0.0    )TSH=-12345.
            if(VSB .eq. 0.0    )TSV=-12345.
c-----
c       define TI constants
c-----
            SR = VSR
            SL = VSR*VSB*VSB
            SN = VSR*VSB*VSB
            SA = VSR*VSA*VSA
            SC = VSR*VSA*VSA
            SF = SA - 2.0*SN
            write(4)r(jd),t0,depthssv(js)-refdep,
     1              depthrsv(jr)-refdep,
     2              TP,TSV,TSH, 
     3              SA, SC, SF, SL, SN, SR
c-----
c       Dislocations and forces must act in a solid source region
c       if the receiver is in a fluid, then permit pressure field output
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
                    if(VRB .lt. 0.0001*VRA.and.VSB.gt.0.0)then
                        ksrc(i) = 1
                    else
                        ksrc(i) = 0
                    endif
                endif
                endif
 5011       continue
            write(4)ksrc
            write(LOT,'(4(1pg11.3),3(1pg11.3))')
     1          r(jd),t0,depthssv(js)-refdep,depthrsv(jr)-refdep,
     2          TP,TSV,TSH
c-----
c       search through the Green s functions and output
c       only is permitted by ksrc. 
c           Note pressure only output for fluid receiver      
c-----
            do 5200 i=n1,n2
                do 5300 jjd=1,ndist
                kk = 0
                do 5301 jjs=1,mdpths
                do 5302 jjr=1,mdpthr
                    kk = kk + 1
                    do 5400 jj=1,21
                        if(jsrc(lsrc(jj)).eq.1)then
                        read(3)datar,datai
                        if(jjd.eq.jd .and. k.eq.kk.and.
     1                      ksrc(jj).ne.0)then
                            write(4)datar*gfac,datai*gfac
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
        close(2,status='delete')
        close(3,status='delete')
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
        common/sourcesv/depthssv(NSOURCE)
        common/receivsv/depthrsv(NRECEIVER)
        real*4 depthssv, depthrsv
        integer lmaxs, lmaxr, mdpths, mdpthr
        common/bound/vmin,vamin,vamax,vbmin,vbmax
        common/damp/alpha,ieqex
        common/modctl/iflsph
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
        real bsh, rhosh, qbsh
        common/modlly/mmax
        common/jout/jsrc(21) , jbdrys, jbdryh
        parameter(NDST=500)
        real*4 r(NDST), tshift(NDST), vred(NDST)
        common/c/cmax,c1,c2,cmin
        character ostr*80
        character mname*80, title*80
        common/lyrctl/lyrins
        logical lyrins
        logical ext

        common/cntrl/ishank,hnkarg,dstcor,dokjar,docausal
        integer*4 dstcor
        real*4 hnkarg
        logical ishank,dokjar,docausal

        common/depref/refdep
        real refdep
        common/earth/radius
        real radius
c-----
c       UNIX FORTRAN - NO CARRIAGE CONTROL
c-----
   21   format(11i5)
   22   format(6e11.4)
   24   format('XLENG=',e15.7,' XFAC=',e15.7)
   30   format(2(1pf15.7)/1x,3i10)
 4012       format('WAVENUMBER FILTERING bounded in phase velocites'/
     1          '[cmax,c1,c2,cmin]=','[', f10.3, ',', 
     2          f10.3, ',', f10.3, ',', f10.3,']' /
     3          '( -1.0 means 0.0 for cmin and infinity for cmax)')
 4013   format('WAVENUMBER FILTERING NOT DONE')
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
            jsrc(20) = 0
            jsrc(21) = 0
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
c           ALL
c-----
            do 1236 i=1,21
                jsrc(i) = 1
 1236       continue
        else if(ieqex.eq.3)then
c-----
c           EXPLOSION ONLY
c-----
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
            jsrc(17) = 1
            jsrc(18) = 1
            jsrc(10) = 1

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
c           P-SV contribution to explosion radial time 
c           history is computed
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
        call adomod()
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
c-----
c       transform the spherical model to a flat model
c-----
        if(iflsph.ne.0)then
            call adosph()
        endif
c-----
c       do not permit Q < 1 If qa or qb is entered > 1 
c       invert to form q inverse
c-----
        do 3007 i=1,mmax
            if(qa(i).lt.0.0)qa(i) = 0.0
            if(qb(i).lt.0.0)qb(i) = 0.0
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
c       transform all spherical source depths to the equivalent depths
c       in flat model - but do not change flat model depths
c-----
        do 1191 i=1,mdpths
            depthssv(i) = depths(i)
            if(iflsph.ne.0)then
C               depths(i) = radius*alog(radius/(radius-depths(i)))
                depths(i) = sngl(
     1               dble(radius)
     2              *dlog(dble(radius)/dble((radius-depths(i)))))
            endif
 1191       continue
            do 1192 i=1,mdpthr
            depthrsv(i) = depthr(i)
            if(iflsph.ne.0)then
C               depthr(i) = radius*alog(radius/(radius-depthr(i)))
                depthr(i) = sngl(
     1               dble(radius)
     2              *dlog(dble(radius)/dble((radius-depthr(i)))))
            endif
 1192       continue
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
c       read in phase velocity limits for wavenumber filtering
c       Wavenumber filtering will consist of following
c
c       |*ZERO*|-COSINE TAPER-|*ALL PASS*|-COSINE TAPER-|*ZERO
c       |      |              |          |              |
c            omega          omega      omega          omega
c k =   0    -----          -----      -----          -----   infinity
c            cmax            c1         c2             cmin
c-----
c       If c2 or cmin <= 0, then upper wavenumber limit is infinite
c       If c1 or cmax <= 0, then lower wavenumber limit is zero
c-----
        read(lun,*,end=4010,err=4010)cmax,c1,c2,cmin
            if(c1.le.0.0)cmax = -1.0
            if(c2.le.0.0)cmin = -1.0
            write(LOT,4012)cmax,c1,c2,cmin
        goto 4011
 4010   continue
            cmax = -1.0
            c1 = -1.0
            c2 = -1.0
            cmin = -1.0
        write(LOT,4013)
 4011   continue
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
        if(mdpths+mdpthr .gt. 2*mmax )then
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
        write(LOT,22)(d(i),a(i),b(i),rho(i),qa(i),qb(i),i=1,mmax)
C-----
C removed 01 APR 2002 since messes up a mixed fluid/solid modium       
Cc-----
Cc      Guarantee that no time wasted if any source is in the water
Cc      since there can only be a center of expansion source
Cc-----
C       do 2019 i=1,mdpths
C           call srclay(depths(i), lmaxs(i), dphs)
C           if(b(lmaxs(i)).le.1.0e-04)then
C               do 2091 ii=1,21
C                   if(ii.ne.7 .and. ii.ne.8 .and. ii.ne.16)then
C                       jsrc(ii) = 0
C                   endif
C 2091          continue
C           endif
C 2019  continue
c-----
c       determine position of source and receiver in new layers
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

        subroutine excit(freq,xleng,xfac,dk,nk)
c-----
c     sample response for all wavenumbers at a given frequency
c     using Bouchon equal wavenumber sampling = dk
c     with offset of 0.218dk
c-----
        parameter(LER=0,LIN=5,LOT=6)
        parameter(NL=200)
        parameter(NSOURCE=100,NRECEIVER=100,NSR=100)
        common/source/depths(NSOURCE),lmaxs(NSOURCE),mdpths
        real depths
        integer lmaxs, mdpths
        common/receiv/depthr(NRECEIVER),lmaxr(NRECEIVER),mdpthr
        real depthr
        integer lmaxr, mdpthr
        common/bound/vmin,vamin,vamax,vbmin,vbmax
        real vmin, vamin, vamax, vbmin, vbmax

        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        real d,a,b,rho,qa,qb,etap,etas,frefp,frefs

C       common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
C       real bsh, rhosh, qbsh

        common/modlly/mmax
        integer mmax
        common/jout/jsrc(21) , jbdrys, jbdryh
        integer jsrc, jbdrys, jdyrh
        common/damp/alpha,ieqex
        real alpha
        integer ieqex
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

        common/layertype/isfluid(NL),allsolid
        logical isfluid,allsolid

        common/shcontrol/ltop,lbot
        integer ltop, lbot

        complex wvn,om, wvn2, om2
        complex gg(21)
c-----
c       matrix components in layers and boundaries saved
c-----
        common/water/iwater(NL),iwats(2),iwatb(2)
        common/lyrctl/lyrins
        logical lyrins
        common/wvflt/wvmn,wvc1,wvc2,wvmx
        real wvmn, wvc1,wvc2, wvmx
        common/c/cmax,c1,c2,cmin
        real cmin, c1, c2, cmax

      common/srcrec/isrc,irec
      integer isrc,irec

        common/jbdry/jtop,jbot
        integer jtop, jbot

        real wv
        integer ii
        logical proceed

c-----


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
c       define angular frequency
c-----
        omega=6.2831853*freq
        om=cmplx((omega),-(alpha))
        om2 = om * om
c-----
c       define wavenumber tapering limits
c-----
        wvmn=omega/cmax
        wvc1=omega/c1
        wvc2=omega/c2
        wvmx=omega/cmin

c-----
c       evaluate wavenumber integration limits
c       and asymptotic coefficients
c-----
        call wvlimit(nk,omega,dk,xfac,xleng)
        rewind 2
c-----
c       output wavenumber in reverse order
c-----
        call bufini(1,ierr)
        do 3998 ii=nk,1,-1
            wv = (ii-1)*dk + 0.218*dk
            call wvfilt(wvmn,wvc1,wvc2,wvmx,wv,fact0)
            wvn=cmplx((wv),0.0d+00)
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
c       evaluate matrices first
c-----
C           if(lyrins)then
C               call evlmat(om,wvn,jbdrys,jbdryh,wvn2,om2)
C           endif
c-----
c       now evaluate for a specific source, receiver position
c-----
            call bufwr(wv)
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
     1                  lmaxs(js), dphs)
                    call srclyr(depthr(jr), 
     1                  lmaxr(jr), dphr)
                    call dezero()
                    call equlck()
C                   call evlmat(om,wvn,jbdrys,jbdryh,wvn2,om2)
                    endif
c               initialize all GG to 0 + 0i
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
                do i=1,21
                   gg(i) = cmplx(0.0,0.0)
                enddo
                    if(ii.le.mkup(k))then
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
                         ltop = 1
                         lbot = mmax 
                         if(allsolid)then
c-----
c                           get the bounds after setting defaults
c-----
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
                         call dosh(om,wvn,gg)
c-----
c                   adjust 12 and 15 to get avoid asymptotic problems
c-----
                         gg(12) = gg(12) * real(wvn)
                         gg(15) = gg(15) * real(wvn)
c-----
c       if asymptotic is invoked, modify the integrand here
c-----
                        if(gasymp(k))then
                             depth = abs(depths(js)-depthr(jr))
                             call gasym(gg,k,wv,depth,jsrc)
                        endif
                    else
                        do 3997 j=1,21
                            gg(j) = cmplx(0.0,0.0)
 3997                   continue
                    endif
                    do 3999 j=1,21
                        if(jsrc(j).eq.1)then
                            gg(j) = gg(j)*fact0
                            call bufwr(real(gg(j)))
                            call bufwr(aimag(gg(j)))
                        endif
 3999               continue
 4010           continue
 4000       continue
 3998   continue
        call buflsh()
        return
        end

        subroutine gcmdln()
c-----
c       parse command line arguments
c       requires subroutine mgtarg() and function mnmarg()
c-----
c       ishank  L   - .true. use Hankel function and not Bessel
c                    note that asymptotic tricks are not used
c       hnkarg  R*4 - (default 6.0) For kr > hnkarg use the Hankel
c                   function, otherwise the Bessel
c       dokjar  L   - .false. use Futterman causal Q
c                 .true.  use Kjartansson Causal Q
c                 .true.  use Kjartansson Causal Q
c       docausal L  - .true. Use causal Q
c-----
        common/cntrl/ishank,hnkarg,dstcor,dokjar,docausal
        integer*4 dstcor
        real*4 hnkarg
        logical ishank,dokjar,docausal
        common/ctlrud/rpup,rpdn,rsup,rsdn,dorud
        logical rpup, rpdn, rsup, rsdn, dorud
        common/ctlsud/spup,spdn,ssup,ssdn,dosud
        logical spup, spdn, ssup, ssdn, dosud

        integer*4 mnmarg
        character*50 name


        ishank = .false.
        hnkarg = 6.0
        dokjar = .false.
        docausal = .true.
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
            if(name(1:2).eq.'-H')then
                ishank = .true.
            else if(name(1:2).eq.'-A')then
                i=i+1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')hnkarg
            else if(name(1:2).eq.'-R' .and. name(1:3).ne.'-RD'
     1              .and. name(1:3).ne.'-RU'
     2              .and. name(1:3).ne.'-RP'
     3              .and. name(1:3).ne.'-RS')then
                dstcor = 2
            else if(name(1:3).eq.'-SU')then
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
            else if(name(1:2).eq.'-N')then
                docausal = .false.
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
        if(hnkarg.lt.3.0)hnkarg=3.0
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
     1  'hspec96 [-H] [-A arg] [-K] [-N]',
     2      '[-SU] [-SD] [-SPUP] [-SSUP] [-SPDN] [-SSDN] ',
     3      '[-RU] [-RD] [-RPUP] [-RSUP] [-RPDN] [-RSDN] [-?] [-h]'
        write(LER,*)
     1  '-H (default false)   Use Hankel function not Bessel'
        write(LER,*)
     1  '-A arg (default arg=3.0) value of kr where Hn(kr) replaces'
        write(LER,*)
     1  '            Jn(kr) in integration - only used when -H is used'
        write(LER,*)
     1  '-K      (default Futterman) use Kjartansson Causal Q'
        write(LER,*)
     1  '-N      (default causal) use non-causal Q'
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

        subroutine reset3(ifreq,jsrc,lsrc,ndist,n1,n2)
c-----
c       routine to reset temporary output file on unit 3
c       if it already exists due to an aborted run
c
c       the file is read until an error is found, which
c       indicates the total number of correct frequencies on the
c       output file. The file is rewound and the correct frequencies
c       are read in to reposition the output file. The
c       parameter ifreq indicates the number of complete frequency sets
c-----
        integer*4 lsrc(*), jsrc(21)
        parameter(NSOURCE=100,NRECEIVER=100,NSR=100)
        common/source/depths(NSOURCE),lmaxs(NSOURCE),mdpths
        common/receiv/depthr(NRECEIVER),lmaxr(NRECEIVER),mdpthr
        real*4 depths, depthr
        integer lmaxs, lmaxr, mdpths, mdpthr
c-----
c       find the correct number by duplicating the 
c       writes of the main program
c-----
        ifreq = 0
        rewind 3
        do 2000 i = n1, n2
            do 2100 jd=1,ndist
            k = 0
            do 1000 js=1,mdpths
            do 1010 jr=1,mdpthr
                k = k + 1
                do 2200 jj=1,21
                if(jsrc(lsrc(jj)).eq.1)then
                    read(3,err=9998,end=9999)xx,yy
                endif
 2200           continue
 1010       continue
 1000       continue
 2100       continue
            ifreq = ifreq + 1
 2000   continue
 9998   continue
 9999   continue
c-----
c       there are now ifreq known data sets out there
c       position write pointer on the output file
c-----
        rewind 3
        do 5000 i = 1,ifreq
            do 5100 jd=1,ndist
            k = 0
            do 4000 js=1,mdpths
            do 4010 jr=1,mdpthr
                k = k + 1
                do 5200 jj=1,21
                if(jsrc(lsrc(jj)).eq.1)then
                    read(3,err=9998,end=9999)xx,yy
                endif
 5200           continue
 4010       continue
 4000       continue
 5100       continue
 5000   continue
        ifreq = ifreq -1 + n1
        return
        end

        subroutine wvint(r,smm,dk,nk)
c-----
c       to work with potentially large disk files, we cannot read in
c       all wavenumbers at once. We only work with neighboring
c       points at any time. The first two are for the DC correction,
c       followed by wavenumbers in decreasing order
c-----
        parameter(NSOURCE=100,NRECEIVER=100,NSR=100)
        common/source/depths(NSOURCE),lmaxs(NSOURCE),mdpths
        common/receiv/depthr(NRECEIVER),lmaxr(NRECEIVER),mdpthr
        real*4 depths, depthr
        integer lmaxs, lmaxr, mdpths, mdpthr
        common/jout/jsrc(21) , jbdrys, jbdryh
        common/asym/j0k0(NSR),j0k1(NSR),j0k2(NSR),j0k3(NSR),
     1      j1k0(NSR),j1k1(NSR),j1k2(NSR),j1k3(NSR),
     2      j2k0(NSR),j2k1(NSR),j2k2(NSR),j2k3(NSR),
     3      j1k0r(NSR),j1k1r(NSR),j1k2r(NSR),j1k3r(NSR),
     4      j1km1r(NSR)
        real*4 j0k0,j0k1,j0k2,j0k3,
     1      j1k0,j1k1,j1k2,j1k3,
     2      j2k0,j2k1,j2k2,j2k3,
     3      j1k0r,j1k1r,j1k2r,j1k3r,
     4      j1km1r
        common/kint4/aa,bb,cc
            complex aa(NSR,21),bb(NSR,21),cc(NSR,21)
        complex gg1(NSR),sumc(NSR)
        complex g(NSR,21)
        complex smm(NSR,21)
        complex sumd(NSR)
        real*4 wvn
        complex h0, h1

        common/cntrl/ishank,hnkarg,dstcor,dokjar,docausal
        integer*4 dstcor
        real*4 hnkarg
        logical ishank,dokjar,docausal
c-----
c       process
c-----
        rewind 2
        call bufini(0,ierr)
c-----
c       initialize integral
c-----
        call intini(smm,r)
c-----
c       we now can procede with integration
c-----
c       in the variables below the t0,j0,j1,sum refer to upper limit
c       of integration and t01,j01,j11 and sum1 refer 
c       to the lower limit
c-----
        do 200 ik = nk,1,-1
            call getgk(g,jsrc,wvn)
            t01 = wvn * r
            dkk = dk
            call hank(t01,h0,h1)
            if(ishank .and. t01 .le. hnkarg)then
                h0 = cmplx(0.0,0.0)
                h1 = cmplx(0.0,0.0)
            endif
c-----
c       perform windowing in wavenumber domain to pass
c       certain ranges of phase velocity
c-----
            if(jsrc(1).eq.1)then
c-----
c       ZDD
c-----
                call fmake(1,0,wvn,1,g,gg1)
                call wint(sumd,gg1,h0,h1,t01,0,dkk,r,wvn)
                call smmadd(1,smm,sumd)
            endif
            if(jsrc(2).eq.1)then
c-----
c       RDD
c-----
                call fmake(2,0,wvn,1,g,gg1)
                call wint(sumd,gg1,h0,h1,t01,1,dkk,r,wvn)
                call smmadd(2,smm,sumd)
            endif
            if(jsrc(3).eq.1)then
c-----
c       ZDS
c-----
                call fmake(3,0,wvn,1,g,gg1)
                call wint(sumd,gg1,h0,h1,t01,1,dkk,r,wvn)
                call smmadd(3,smm,sumd)
            endif
            call smmzer(sumc)
c-----
c           only include near field term if both SH and P-SV
c           computed
c-----
            if(jsrc(4).eq.1.and.jsrc(13).eq.1)then
                call fmake(4,13,wvn,0,g,gg1)
                call wint(sumc,gg1,h0,h1,t01,10,dkk,r,wvn)
            endif
            if(jsrc(4).eq.1)then
c-----
c       RDS
c-----
                call fmake(4,0,wvn,1,g,gg1)
                call wint(sumd,gg1,h0,h1,t01,0,dkk,r,wvn)
                call smmadd(4,smm,sumd)
                call smmadd(4,smm,sumc)
            endif
            if(jsrc(13).eq.1)then
c-----
c       TDS
c-----
                call fmake(13,0,wvn,1,g,gg1)
                call wint(sumd,gg1,h0,h1,t01,0,dkk,r,wvn)
                call smmadd(5,smm,sumd)
                call smmadd(5,smm,sumc)
            endif
            if(jsrc(5).eq.1)then
c-----
c       ZSS
c-----
                call fmake(5,0,wvn,1,g,gg1)
                call wint(sumd,gg1,h0,h1,t01,2,dkk,r,wvn)
                call smmadd(6,smm,sumd)
            endif
            call smmzer(sumc)
c-----
c           only include near field term if both SH and P-SV
c           computed
c-----
            if(jsrc(6).eq.1.and.jsrc(14).eq.1)then
                call fmake(6,14,wvn,0,g,gg1)
                call wint(sumc,gg1,h0,h1,t01,20,dkk,r,wvn)
            endif
            if(jsrc(6).eq.1)then
c-----
c       RSS
c-----
                call fmake(6,0,wvn,1,g,gg1)
                call wint(sumd,gg1,h0,h1,t01,1,dkk,r,wvn)
                call smmadd(7,smm,sumd)
                call smmadd(7,smm,sumc)
            endif
            if(jsrc(14).eq.1)then
c-----
c       TSS
c-----
                call fmake(14,0,wvn,1,g,gg1)
                call wint(sumd,gg1,h0,h1,t01,1,dkk,r,wvn)
                call smmadd(8,smm,sumd)
                call smmadd(8,smm,sumc)
            endif
            if(jsrc(7).eq.1)then
c-----
c       ZEX
c-----
                call fmake(7,0,wvn,1,g,gg1)
                call wint(sumd,gg1,h0,h1,t01,0,dkk,r,wvn)
                call smmadd(9,smm,sumd)
            endif
            if(jsrc(8).eq.1)then
c-----
c       REX
c-----
                call fmake(8,0,wvn,1,g,gg1)
                call wint(sumd,gg1,h0,h1,t01,1,dkk,r,wvn)
                call smmadd(10,smm,sumd)
            endif
            if(jsrc(9).eq.1)then
c-----
c       ZVF
c-----
                call fmake( 9,0,wvn,1,g,gg1)
                call wint(sumd,gg1,h0,h1,t01,0,dkk,r,wvn)
                call smmadd(11,smm,sumd)
            endif
            if(jsrc(10).eq.1)then
c-----
c       RVF
c-----
                call fmake(10,0,wvn,1,g,gg1)
                call wint(sumd,gg1,h0,h1,t01,1,dkk,r,wvn)
                call smmadd(12,smm,sumd)
            endif
            if(jsrc(11).eq.1)then
c-----
c       ZHF
c-----
                call fmake(11,0,wvn,1,g,gg1)
                call wint(sumd,gg1,h0,h1,t01,1,dkk,r,wvn)
                call smmadd(13,smm,sumd)
            endif
            call smmzer(sumc)
c-----
c           include near field term only if both SH and P-SV
c           computed
c-----
            if(jsrc(12).eq.1.and.jsrc(15).eq.1)then
                call fmake(12,15,wvn,0,g,gg1)
                call wint(sumc,gg1,h0,h1,t01,-10,dkk,r,wvn)
            endif
            if(jsrc(12).eq.1)then
c-----
c       RHF
c-----
                call fmake(12,0,wvn,0,g,gg1)
                call wint(sumd,gg1,h0,h1,t01,0,dkk,r,wvn)
                call smmadd(14,smm,sumd)
                call smmadd(14,smm,sumc)
            endif
            if(jsrc(15).eq.1)then
c-----
c       THF
c-----
                call fmake(15,0,wvn,0,g,gg1)
                call wint(sumd,gg1,h0,h1,t01,0,dkk,r,wvn)
                call smmadd(15,smm,sumd)
                call smmadd(15,smm,sumc)
            endif
            if(jsrc(16).eq.1)then
c-----
c       PEX
c-----
                call fmake(16,0,wvn,1,g,gg1)
                call wint(sumd,gg1,h0,h1,t01,0,dkk,r,wvn)
                call smmadd(16,smm,sumd)
            endif
            if(jsrc(17).eq.1)then
c-----
c       PDD
c-----
                call fmake(17,0,wvn,1,g,gg1)
                call wint(sumd,gg1,h0,h1,t01,0,dkk,r,wvn)
                call smmadd(17,smm,sumd)
            endif
            if(jsrc(18).eq.1)then
c-----
c       PDS
c-----
                call fmake(18,0,wvn,1,g,gg1)
                call wint(sumd,gg1,h0,h1,t01,1,dkk,r,wvn)
                call smmadd(18,smm,sumd)
            endif
            if(jsrc(19).eq.1)then
c-----
c       PSS
c-----
                call fmake(19,0,wvn,1,g,gg1)
                call wint(sumd,gg1,h0,h1,t01,2,dkk,r,wvn)
                call smmadd(19,smm,sumd)
            endif
            if(jsrc(20).eq.1)then
c-----
c       PVF
c-----
                call fmake(20,0,wvn,1,g,gg1)
                call wint(sumd,gg1,h0,h1,t01,0,dkk,r,wvn)
                call smmadd(20,smm,sumd)
            endif
            if(jsrc(21).eq.1)then
c-----
c       PHF
c-----
                call fmake(21,0,wvn,1,g,gg1)
                call wint(sumd,gg1,h0,h1,t01,1,dkk,r,wvn)
                call smmadd(21,smm,sumd)
            endif
  200   continue
c-----
c       sign change due to k j(-1)
c-----
        call smmflp(smm)
        return
        end

        subroutine wint(smm,g1,h0,h1,t01,n,dkk,r,wvn)
c-----
c       perform the wavenumber integration for the particular
c       integrand
c-----
        parameter(NSOURCE=100,NRECEIVER=100,NSR=100)
        common/source/depths(NSOURCE),lmaxs(NSOURCE),mdpths
        common/receiv/depthr(NRECEIVER),lmaxr(NRECEIVER),mdpthr
        real*4 depths, depthr
        integer lmaxs, lmaxr, mdpths, mdpthr
        complex smm(NSR),g1(NSR)
        complex h0, h1
        k = 0
        do 1000 js=1,mdpths
        do 1010 jr=1,mdpthr
            k = k + 1
            call wwint(smm(k),g1(k),h0,h1,t01,n,dkk,r,wvn)
 1010   continue
 1000   continue
        return
        end

        subroutine wwint(smm,g1,h0,h1,t01,n,dkk,r,wvn)
        complex smm,g1
        complex h0, h1
        complex h2
        real j21
        common/rlimit/rlim
        real*4 rlim
        common/cntrl/ishank,hnkarg,dstcor,dokjar,docausal
        integer*4 dstcor
        real*4 hnkarg
        logical ishank,dokjar,docausal
c-----
c       test for Hankel Function
c-----
        if(ishank)then
c-----
c           rectangular rule
c-----
            if(n.eq.0)then
c-----
c           integral (c + d z) * j0(z) dz
c-----
                smm =  (g1 * (h0 * dkk))
            elseif(n.eq.1)then
c-----
c           integral (c + d z) j1(z) dz
c-----
                smm =  (g1 * (h1 * dkk))
c-----
c           integral (c + d z) j2(z) dz
c-----
            elseif(n.eq.2)then
                if(t01.eq.0.0)then
                    h2 = (0.0,0.0)
                else
                    h2 = 2.*h1/t01 - h0
                endif
                smm =  (g1 * (h2 * dkk))
c-----
c            - integral (c + d z) j1(kr) dk / r
c-----
            elseif(n.eq.10)then
                if(r.le.rlim)then
                    smm =  (g1 * (0.5 * dkk)) * wvn
                else
                    smm =  (g1 * (h1 * dkk))/r
                endif
                smm = - smm
        SMM =  CMPLX(0.0,0.0)
c-----
c            - integral (c + d z) j1(kr) dk / kr
c-----
            elseif(n.eq.-10)then
                if(r.le.rlim)then
                    smm =  (g1 * (0.5 * dkk)) 
                else
                    smm =  (g1 * (h1 * dkk))/(r*wvn)
                endif
                smm = - smm
        SMM =  CMPLX(0.0,0.0)
c-----
c            - 2 integral (c + d z) j2(kr) dk / r
c-----
            elseif(n.eq.20)then
                if(t01.eq.0.0)then
                    h2 = (0.0,0.0)
                else
                    h2 = 2.*h1/t01 - h0
                endif
                if(r.le.rlim)then
                    smm = (0.0,0.0)
                else
                    smm =  (g1 * (h2 * dkk)) / r
                endif
                smm = -2.0*smm
        SMM =  CMPLX(0.0,0.0)
            endif
        else
c-----
c           rectangular rule
c-----
            if(n.eq.0)then
c-----
c           integral (c + d z) * j0(z) dz
c-----
                smm =  (g1 * (real(h0) * dkk))
            elseif(n.eq.1)then
c-----
c           integral (c + d z) j1(z) dz
c-----
                smm =  (g1 * (real(h1) * dkk))
c-----
c           integral (c + d z) j2(z) dz
c-----
            elseif(n.eq.2)then
                if(t01.eq.0.0)then
                    j21 = (0.0,0.0)
                else
                    j21 = 2.*real(h1)/t01 - real(h0)
                endif
                smm =  (g1 * (j21 * dkk))
c-----
c            - integral (c + d z) j1(kr) dk / r
c-----
            elseif(n.eq.10)then
                if(r.le.rlim)then
                    smm =  (g1 * (0.5 * dkk)) * wvn
                else
                    smm =  (g1 * (real(h1) * dkk))/r
                endif
                smm = - smm
c-----
c            - integral (c + d z) j1(kr) dk / kr
c-----
            elseif(n.eq.-10)then
                if(r.le.rlim)then
                    smm =  (g1 * (0.5 * dkk)) 
                else
                    smm =  (g1 * (real(h1) * dkk))/(r*wvn)
                endif
                smm = - smm
c-----
c            - 2 integral (c + d z) j2(kr) dk / r
c-----
            elseif(n.eq.20)then
                if(t01.eq.0.0)then
                    j21 = (0.0,0.0)
                else
                    j21 = 2.*real(h1)/t01 - real(h0)
                endif
                if(r.le.rlim)then
                    smm = (0.0,0.0)
                else
                    smm =  (g1 * (j21 * dkk)) / r
                endif
                smm = -2.0*smm
            endif
        endif
        return
        end

        subroutine hank(z,h0,h1)
c-----
c       evaluate Bessel functions using Abromiwitz and Stegun
c-----
        complex h0, h1
        real z
        real j0, j1
        real j1z
        common/cntrl/ishank,hnkarg,dstcor,dokjar,docausal
        integer*4 dstcor
        real*4 hnkarg
        logical ishank,dokjar,docausal
c-----
        h0 = cmplx(0.0,0.0)
        h1 = cmplx(0.0,0.0)
        if(z.eq.0.0)then
            h0 = cmplx(1.0,0.0)
        elseif(z.gt.0.0 .and. z.le.3.0)then
            x = (z/3.)*(z/3.)
            j0 = 1.-x*(2.2499997-x*(1.2656208-x*(.3163866-x*(
     1            .0444479-x*(.0039444-x*(.0002100))))))
            j1z = 0.5-x*(.56249985-x*(.21093573-x*(.03954289-x*(
     1      .00443319-x*(.00031761-x*(.00001109))))))
            j1 = z * j1z
            h0 = cmplx(j0,0.0)
            h1 = cmplx(j1,0.0)
        else
            x = 3./z
            fac = 1./sqrt(z)
            f0 = .79788456+x*(-.00000077 + x*(-.00552740 + x*(
     1      -.00009512+x*(.00137237+x*(-.00072805+x*(.00014476))))
     2            ))
            t0 = z - .78539816+x*(-.04166397+x*(-.00003954+x*(
     1      .00262573+x*(-.00054125+x*(-.00029333+x*(.00013558))))
     2            ))
            f1 = .79788456+x*(.00000156+x*(.01659667+x*(.00017105+
     1            x*(-.00249511+x*(.00113653+x*(-.00020033))))))
            t1 = z-2.35619449+x*(.12499612+x*(.00005650+x*(
     1       -.00637879+x*(.00074348+x*(.00079824+x*(-.00029166)))
     2            )))
            if(ishank .and. z .ge. hnkarg)then
                h0 = 0.5 * fac * f0 * cmplx(cos(t0), -sin(t0))
                h1 = 0.5 * fac * f1 * cmplx(cos(t1), -sin(t1))
            else
                h0 =       fac * f0 * cmplx(cos(t0),0.0)
                h1 =       fac * f1 * cmplx(cos(t1),0.0)
            endif
        endif
        return
        end

        subroutine solu(y1,y2,x1,x2,h,j,a,b,c)
c-----
c       Using two data points, determine the coefficients of
c       a function
c       y(k) = [ a + bk + ck^2 ] exp [ -kh]
c
c       Given only two data points, we will constrain one of the a,b,c
c       to be zero, depending on the nature of the 
c       elastic wave integrand
c-----
c       we do not solve for a,b,c together, only two at most
c       thus we only need two values of wavenumber, x1 and x2
c
c
c       Since the program permits consideration of more than a
c       single depth, we must check for overflow here and 
c       underflow in gasym
c-----
        complex y1,y2,a,b,c
        integer imap(21)
C       data imap/3,2,3,2,3,2,3,2,3,2,3,3,2,2,3,3,3,3,3,3,3/
        data imap/2,2,3,2,3,3,2,3,2,2,2,3,2,2,3,2,2,2,2,2,2/
        a=cmplx(0.0,0.0)
        b=cmplx(0.0,0.0)
        c=cmplx(0.0,0.0)
        wfac = x1*h
c       if(wfac.gt.10.0)return
        ii = imap(j)
        if(ii.eq.1)then
c-----
c       a exp(-kh)
c-----
            b=cmplx(0.0,0.0)
            a=y1*exp(x1*h)
        else if(ii.eq.2)then
c-----
c       [ a + b k  ]exp(-kh)
c-----
            u1=x1*h
            u2=x2*h
            det=x2-x1
            a= x2*(y1*exp(u1))-x1*(y2*exp(u2))
            a=a/det
            b= y2*exp(u2) - y1*exp(u1)
            b=b/det
        else if(ii.eq.3)then
c-----
c       [ b + c k  ]  k exp(-kh)
c-----
            u1=x1*h
            u2=x2*h
            det=x2-x1
            a = cmplx(0.0,0.0)
            b = x2*(y1*exp(u1))/x1 - x1*(y2*exp(u2))/x2
            b = b/det
            c= y2*exp(u2)/x2 - y1*exp(u1)/x1
            c = c/det
        else if(ii.eq.4)then
c-----
c       [ c k*k ] exp(-kh)
c-----
            a = cmplx(0.0,0.0)
            b = cmplx(0.0,0.0)
            c = y1 * exp(x1*h)/(x1)**2
        else if(ii.eq.5)then
c-----
c       [ b k ] exp(-kh)
c-----
            a = cmplx(0.0,0.0)
            b = y1 * exp(x1*h)/ x1
        endif
        return
        end

        subroutine setup(rr) 
c---------------------------------------------------------- 
c 
c       jnkm =  integral exp(-kh) krsup m j sub n (kr) dk 
c
c       This is used in the fit of low frequency information
c
c       integral f(k) Jn(kr) dk = 
c           integral [ f(k) - (a+bk+ck^2)e^{-kh} ] Jn(kr) dk
c           +integral [  (a+bk+ck^2)e^{-kh} ] Jn(kr) dk
c
c       The last integral is obtained analytically
c
c       Special care is taken when r=0, especially for a near field
c       TDS, RDS term
c       j1kmr = lim r->0 integral exp(-kh) k rsup m j sub 1 (kr) dk
c
c       Herrmann, R. B., and C. Y. Wang (1985).
c       A comparison of synthetic seismograms,
c       Bull. Seism. Soc. Am. 75, 41-56.
c 
c---------------------------------------------------------- 
        implicit double precision (a-h,o-z)
        parameter(NSOURCE=100,NRECEIVER=100,NSR=100)
        common/source/depths(NSOURCE),lmaxs(NSOURCE),mdpths
        common/receiv/depthr(NRECEIVER),lmaxr(NRECEIVER),mdpthr
        real*4 depths, depthr
        integer lmaxs, lmaxr, mdpths, mdpthr
        real*4 rr,zz
        common/rlimit/rlim
        real*4 rlim
        common/asym/j0k0(NSR),j0k1(NSR),j0k2(NSR),j0k3(NSR),
     1      j1k0(NSR),j1k1(NSR),j1k2(NSR),j1k3(NSR),
     2      j2k0(NSR),j2k1(NSR),j2k2(NSR),j2k3(NSR),
     3      j1k0r(NSR),j1k1r(NSR),j1k2r(NSR),j1k3r(NSR),
     4      j1km1r(NSR)
        real*4 j0k0,j0k1,j0k2,j0k3,
     1      j1k0,j1k1,j1k2,j1k3,
     2      j2k0,j2k1,j2k2,j2k3,
     3      j1k0r,j1k1r,j1k2r,j1k3r,
     4      j1km1r
        k = 0
        do 2500 i=1,mdpths
            do 2600 j=1,mdpthr
                zz = abs(depths(i) - depthr(j))
                k = k + 1
                    r = dble(rr)
                    z = dble(zz)
                    dist=dsqrt(r*r + z*z) 
c-----
c       if distance == 0 , force small answers
c-----
                if(dist.le.rlim)dist=rlim
                    dist3=dist**3 
                    dist5=dist**5 
                    dist7=dist**7 
                zdist = z + dist
                    rz=r*z 
                    z2=z*z 
                    r2=r*r 
                    r3=r*r2 
                    z3=z*z2 
                    rz2=r*z2 
                    rz3=r*z3 
                    zor = z/dist
                    zor2= zor*zor
                    zor3= zor*zor2
                    j0k0(k) = sngl( 1.0/dist   )
                    j0k1(k) = sngl( z/dist3   )
                    j0k2(k) = sngl( (2.0*z2 - r2)/dist5   )
                    j0k3(k) = sngl( (6.0*z3 - 9.0*z*r2)/dist7   )
                if(rr.le.rlim)then
                    j1k0(k) = 0.0
                    j1k1(k) = 0.0
                    j1k2(k) = 0.0
                    j1k3(k) = 0.0
                    if(zz .le.rlim)then
                        j1km1r(k)=sngl(1.0d+00/zdist)
                        j1k0r(k) = sngl(0.5/(dist*dist))
                    else
                        j1km1r(k) = sngl(0.5/z)
                        j1k0r(k) = sngl(0.5/z2)
                    endif
                    j1k1r(k) = sngl(1.0/dist3)
                        j1k2r(k) = sngl( 3.0*z/dist5   )
                        j1k3r(k) = sngl( 12.0*z2 /dist7   )
                    
                else
c                       j1k0(k) = sngl(( 1.0 -z/dist   )/r       )
                        j1k0(k) = sngl((r/(dist*zdist) )      )
                        j1k1(k) = sngl( r/dist3                  )
                        j1k2(k) = sngl( 3.0*z*r/dist5            )
                        j1k3(k) = sngl( 3.0*r*(4.0*z2 - r2)/dist7)
                    j1km1r(k)=sngl(1.0d+00/(zdist)          )
c                   j1k0r(k) = j1k0(k)/rr
                        j1k0r(k) = sngl(1.0d+00/(dist*zdist))
                    j1k1r(k) = j1k1(k)/rr
                    j1k2r(k) = j1k2(k)/rr
                    j1k3r(k) = j1k3(k)/rr
                endif
                if(rr.le.rlim)then
                        j2k0(k) = 0.0
                        j2k1(k) = 0.0
                else
                        j2k0(k)=sngl((1.0 -zor)*(1.0-zor)*(dist/r2))
c                       j2k0(k)=sngl(((r/zdist)**2)/dist)
                        j2k1(k)=sngl((1.0-zor)*(1.0-zor)*(2.0+zor)/r2)
c                       j2k1(k)=j2k0(k)*sngl( (zdist+dist)/dist3)
                endif
                    j2k2(k) = sngl( 3.0*r2/dist5   )
                    j2k3(k) = sngl( 15.0*z*r2/dist7   )
 2600       continue
 2500   continue
        return 
        end 

        subroutine fmake(j,k,wvn,l,g,gg1)
c-----
c       make the proper integrand
c-----
        parameter(NSOURCE=100,NRECEIVER=100,NSR=100)
        common/source/depths(NSOURCE),lmaxs(NSOURCE),mdpths
        common/receiv/depthr(NRECEIVER),lmaxr(NRECEIVER),mdpthr
        real*4 depths, depthr
        integer lmaxs, lmaxr, mdpths, mdpthr
        real*4 wvn
        complex g(NSR,21)
        complex gg1(NSR)
        kk = 0
        do 1000 js=1,mdpths
        do 1010 jr=1,mdpthr
            kk = kk + 1
            gg1(kk)= g(kk,j)
            if(k.gt.0)then
                    gg1(kk)=gg1(kk) + g(kk,k)
            endif
            if(l.gt.0)then
                    gg1(kk)=gg1(kk) * wvn
            endif
 1010   continue
 1000   continue
        return
        end

        subroutine getgk(g,jsrc,wvno)
c-----
c       read input to obtain elements of g(21,j) array
c-----
        parameter(NSOURCE=100,NRECEIVER=100,NSRC=100)
        common/source/depths(NSOURCE),lmaxs(NSOURCE),mdpths
        common/receiv/depthr(NRECEIVER),lmaxr(NRECEIVER),mdpthr
        real*4 depths, depthr
        integer lmaxs, lmaxr, mdpths, mdpthr
        complex g(NSRC,21)
        dimension jsrc(21)
        call bufrd(wvno,ierr)
        k = 0
        do 1000 js=1,mdpths
        do 1010 jr=1,mdpthr
            k = k + 1
            do 101 i=1,21
                if(jsrc(i).eq.1)then
                    call bufrd(xr,ierr)
                    call bufrd(xi,ierr)
                    g(k,i)=cmplx(xr,xi)
                endif
  101       continue
 1010   continue
 1000   continue
        return
        end

        subroutine intini(smm,r)
        parameter(NSOURCE=100,NRECEIVER=100,NSR=100)
        common/source/depths(NSOURCE),lmaxs(NSOURCE),mdpths
        common/receiv/depthr(NRECEIVER),lmaxr(NRECEIVER),mdpthr
        real*4 depths, depthr
        integer lmaxs, lmaxr, mdpths, mdpthr
        common/jout/jsrc(21) , jbdrys, jbdryh
        common/asym/j0k0(NSR),j0k1(NSR),j0k2(NSR),j0k3(NSR),
     1      j1k0(NSR),j1k1(NSR),j1k2(NSR),j1k3(NSR),
     2      j2k0(NSR),j2k1(NSR),j2k2(NSR),j2k3(NSR),
     3      j1k0r(NSR),j1k1r(NSR),j1k2r(NSR),j1k3r(NSR),
     4      j1km1r(NSR)
        real*4 j0k0,j0k1,j0k2,j0k3,
     1      j1k0,j1k1,j1k2,j1k3,
     2      j2k0,j2k1,j2k2,j2k3,
     3      j1k0r,j1k1r,j1k2r,j1k3r,
     4      j1km1r
        common/kint4/aa,bb,cc
            complex aa(NSR,21),bb(NSR,21),cc(NSR,21)
        complex smm(NSR,21),sumd
        common/kint1/gasymp
            logical gasymp(NSR)
        common/rlimit/rlim
        real*4 rlim
        k = 0
        do 1000 js=1,mdpths
        do 1010 jr=1,mdpthr
            k = k + 1
            if(gasymp(k))then
c-----
c       set up sum arrays, but put in asymptotic value now
c       of setting to zero and then resetting
c-----
                smm(k,1)=         aa(k,1)*j0k1(k) 
     1                  + bb(k,1)*j0k2(k) 
     2                  + cc(k,1)*j0k3(k)
                smm(k,2)=         aa(k,2)*j1k1(k) 
     1                  + bb(k,2)*j1k2(k) 
     2                  + cc(k,2)*j1k3(k)
                smm(k,3)=         aa(k,3)*j1k1(k) 
     1                  + bb(k,3)*j1k2(k) 
     2                  + cc(k,3)*j1k3(k)
                if(jsrc(4).eq.1 .and. jsrc(13).eq.1)then
                    sumd  = (aa(k,4)+aa(k,13))*j1k0r(k) 
     1                  + (bb(k,4)+bb(k,13))*j1k1r(k) 
     2                  + (cc(k,4)+cc(k,13))*j1k2r(k)
                    sumd = - sumd
                else
                    sumd = cmplx(0.0,0.0)
                endif
                smm(k,4)= sumd  + aa(k,4) *j0k1(k) 
     1                  + bb(k,4) *j0k2(k) 
     2                  + cc(k,4) *j0k3(k)
                smm(k,5)= sumd  + aa(k,13)*j0k1(k) 
     1                  + bb(k,13)*j0k2(k) 
     2                  + cc(k,13)*j0k3(k)
                smm(k,6)=     aa(k,5)*j2k1(k) 
     1                  + bb(k,5)*j2k2(k) 
     2                  + cc(k,5)*j2k3(k)
                if(jsrc(6).eq.1 .and. jsrc(14).eq.1 .and.
     1              r.gt.rlim)then
                    sumd= (aa(k,6)+aa(k,14))*j2k0(k) 
     1                  + (bb(k,6)+bb(k,14))*j2k1(k) 
     2                  + (cc(k,6)+cc(k,14))*j2k2(k)
                    sumd = -2.*sumd/r
                else
                    sumd = cmplx(0.0,0.0)
                endif
                smm(k,7)= sumd  + aa(k,6) *j1k1(k) 
     1                  + bb(k,6) *j1k2(k) 
     2                  + cc(k,6) *j1k3(k)
                smm(k,8)= sumd  + aa(k,14)*j1k1(k)
     1                  + bb(k,14)*j1k2(k)
     2                  + cc(k,14)*j1k3(k)
                smm(k,9)=         aa(k,7) *j0k1(k)        
     1                  + bb(k,7) *j0k2(k) 
     2                  + cc(k,7) *j0k3(k)
                smm(k,10)=        aa(k,8) *j1k1(k) 
     1                  + bb(k,8) *j1k2(k) 
     2                  + cc(k,8) *j1k3(k)
                smm(k,11)=        aa(k,9) *j0k1(k) 
     1                  + bb(k,9) *j0k2(k) 
     2                  + cc(k,9) *j0k3(k)
                smm(k,12)=        aa(k,10)*j1k1(k) 
     1                  + bb(k,10)*j1k2(k) 
     2                  + cc(k,10)*j1k3(k)
                smm(k,13)=        aa(k,11)*j1k1(k) 
     1                  + bb(k,11)*j1k2(k) 
     2                  + cc(k,11)*j1k3(k)
                if(jsrc(12).eq.1 .and. jsrc(15).eq.1)then
                    sumd  = (aa(k,12)+aa(k,15))*j1km1r(k) 
     1                  + (bb(k,12)+bb(k,15))*j1k0r(k) 
     2                  + (cc(k,12)+cc(k,15))*j1k1r(k)
                    sumd = - sumd
                else
                    sumd = cmplx(0.0,0.0)
                endif
                smm(k,14)= sumd + aa(k,12)*j0k0(k) 
     1                  + bb(k,12)*j0k1(k) 
     2                  + cc(k,12)*j0k2(k)
                smm(k,15)= sumd + aa(k,15)*j0k0(k) 
     1                  + bb(k,15)*j0k1(k) 
     2                  + cc(k,15)*j0k2(k)
                smm(k,16)=        aa(k,16)*j0k1(k) 
     1                  + bb(k,16)*j0k2(k) 
     2                  + cc(k,16)*j0k3(k)
                smm(k,17)=        aa(k,17)*j0k1(k) 
     1                  + bb(k,17)*j0k2(k) 
     2                  + cc(k,17)*j0k3(k)
                smm(k,18)=        aa(k,18)*j1k1(k) 
     1                  + bb(k,18)*j1k2(k) 
     2                  + cc(k,18)*j1k3(k)
                smm(k,19)=        aa(k,19)*j2k1(k) 
     1                  + bb(k,19)*j2k2(k) 
     2                  + cc(k,19)*j2k3(k)
                smm(k,20)=        aa(k,20)*j0k1(k) 
     1                  + bb(k,20)*j0k2(k) 
     2                  + cc(k,20)*j0k3(k)
                smm(k,21)=        aa(k,21)*j1k1(k) 
     1                  + bb(k,21)*j1k2(k) 
     2                  + cc(k,21)*j1k3(k)
            else
                do 100 i=1,21
                    smm(k,i)=cmplx(0.0,0.0)
  100           continue
            endif
 1010   continue
 1000   continue
        return
        end

        subroutine wvfilt(wvmn,wvc1,wvc2,wvmx,wvno,fact)
        common/c/cmax,c1,c2,cmin
c-----
c       apply a cosine taper wavenumber filter
c-----
        pi = 3.1415927
c-----
c       test if no filtering is to be done
c-----
        if(cmin .le. 0.0 .and. cmax.le.0.0)then
            fact = 1.0
c-----
c       test if high pass in wavenumber (low phase velocity)
c-----
        else if(cmin .le. 0.0)then
            if(wvno.ge.wvc1)then
                fact=1.0
            else if(wvno.ge.wvmn.and.wvno.lt.wvc1)then
                fact=(1.-cos(pi*(wvno-wvmn)/ (wvc1-wvmn)))/2.
            else if(wvno.lt.wvmn)then
                fact = 0.0
            endif
c-----
c       test if low pass in wavenumber (high phase velocity)
c-----
        else if(cmax .le. 0.0)then
            if(wvno.le.wvc2)then
                fact = 1.0 
            else if(wvno.gt.wvc2.and.wvno.le.wvmx)then
                fact=(1.-cos(pi*(wvno-wvmx)/ (wvc2-wvmx)))/2.
            else
                fact = 0.0
            endif
c-----
c       test if band pass in wavenumber
c-----
        else
            if(wvno.ge.wvc1.and.wvno.le.wvc2)then
                fact=1.0
            elseif(wvno.ge.wvmn.and.wvno.lt.wvc1)then
                fact=(1.-cos(pi*(wvno-wvmn)/ (wvc1-wvmn)))/2.
            elseif(wvno.gt.wvc2.and.wvno.le.wvmx)then
                fact=(1.-cos(pi*(wvno-wvmx)/ (wvc2-wvmx)))/2.
            else
                fact = 0.0
            endif
        endif
        return
        end

        subroutine smmadd(i,smm,sumd)
        parameter(NSOURCE=100,NRECEIVER=100,NSR=100)
        common/source/depths(NSOURCE),lmaxs(NSOURCE),mdpths
        common/receiv/depthr(NRECEIVER),lmaxr(NRECEIVER),mdpthr
        real*4 depths, depthr
        integer lmaxs, lmaxr, mdpths, mdpthr
        complex smm(NSR,21), sumd(NSR)
c-----
c       add sumd vector to the particular entry of smm
c-----
        k = 0
        do 1000 js=1,mdpths
            do 1010 jr=1,mdpthr
                k = k + 1
                smm(k,i) = smm(k,i) +  sumd(k)
 1010       continue
 1000   continue
        return
        end

        subroutine smmzer(sumc)
c-----
c       zero an array for all depths
c-----
        parameter(NSOURCE=100,NRECEIVER=100,NSR=100)
        common/source/depths(NSOURCE),lmaxs(NSOURCE),mdpths
        common/receiv/depthr(NRECEIVER),lmaxr(NRECEIVER),mdpthr
        real*4 depths, depthr
        integer lmaxs, lmaxr, mdpths, mdpthr
        complex sumc(NSR)
c-----
c       zero a vector over all depths
c-----
        k = 0
        do 1000 js=1,mdpths
            do 1010 jr=1,mdpthr
                k = k + 1
                sumc(k) = cmplx(0.0,0.0)
 1010       continue
 1000   continue
        return
        end

        subroutine smmscl(sumc,scl)
c-----
c       scale an array for all depths
c-----
        parameter(NSOURCE=100,NRECEIVER=100,NSR=100)
        common/source/depths(NSOURCE),lmaxs(NSOURCE),mdpths
        common/receiv/depthr(NRECEIVER),lmaxr(NRECEIVER),mdpthr
        real*4 depths, depthr
        integer lmaxs, lmaxr, mdpths, mdpthr
        complex sumc(NSR)
        real*4 scl
c-----
c       add sumd vector to the particular entry of smm
c-----
        k = 0
        do 1000 js=1,mdpths
            do 1010 jr=1,mdpthr
                k = k + 1
                sumc(k) = sumc(k) * scl
 1010       continue
 1000   continue
        return
        end

        subroutine smmflp(smm)
c-----
c       flip the -1 Bessel function values
c-----
        parameter(NSOURCE=100,NRECEIVER=100,NSR=100)
        common/source/depths(NSOURCE),lmaxs(NSOURCE),mdpths
        common/receiv/depthr(NRECEIVER),lmaxr(NRECEIVER),mdpthr
        real*4 depths, depthr
        integer lmaxs, lmaxr, mdpths, mdpthr
        complex smm(NSR,21)
        k = 0
        do 1000 js=1,mdpths
            do 1010 jr=1,mdpthr
                k = k + 1
                smm(k,2)  = -smm(k,2)
                smm(k,10) = -smm(k,10)
                smm(k,12) = -smm(k,12)
 1010       continue
 1000   continue
        return
        end

        subroutine gett0(t0,r,depths,depthr,tshift,vred,dstcor)
        real*4 t0, depths, depthr, tshift, vred
        integer*4 dstcor
c-----
c       compute time of first sample of the time series
c-----
            if(dstcor.eq.0)then
                rr = r
            else if(dstcor.eq.1)then
                rr = abs(depthr - depths)
            else if(dstcor.eq.2)then
                rr = sqrt(r*r + (depthr-depths)*(depthr-depths))
            endif
            if(vred.eq.0.0)then
                t0 = tshift 
            else
                t0 = tshift + rr/vred
            endif
        return
        end

        subroutine wvlimit(nk,omega,dk,xfac,xleng)
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
        real bsh, rhosh, qbsh
        common/modlly/mmax
        common/jout/jsrc(21) , jbdrys, jbdryh
        parameter(NSOURCE=100,NRECEIVER=100,NSR=100)
        common/source/depths(NSOURCE),lmaxs(NSOURCE),mdpths
        common/receiv/depthr(NRECEIVER),lmaxr(NRECEIVER),mdpthr
        common/frlim/fl,fu,df,fwhich
        real*4 depths, depthr
        integer lmaxs, lmaxr, mdpths, mdpthr
        common/bound/vmin,vamin,vamax,vbmin,vbmax
        common/damp/alpha,ieqex
        common/lyrctl/lyrins
        logical lyrins
        common/wvflt/wvmn,wvc1,wvc2,wvmx
        common/c/cmax,c1,c2,cmin
c-----
c       set up common blocks for wavenumber sampling at
c       suitable depths. This is necessary since the matrix
c       evaluation is done here for all source-receiver pairs
c       The source-receiver distance is important for the
c       wavenumber sampling at low frequencies
c-----
c       If wavenumber filtering is used and if 
c          the upper wavenumber limit
c       is not defined, then asymptotic may be permissible
c-----
        common/kint1/gasymp
            logical gasymp(NSR)
        common/kint2/mkup
            integer mkup(NSR)
        common/kint3/wave
            real*4 wave(NSR,2)
        common/kint4/aa,bb,cc
            complex aa(NSR,21),bb(NSR,21),cc(NSR,21)
        complex wvn,om, wvn2, om2
        complex g1(21), g2(21)
c-----
c       command line arguments
c-----
        common/cntrl/ishank,hnkarg,dstcor,dokjar,docausal
        integer*4 dstcor
        real*4 hnkarg
        logical ishank,dokjar,docausal

      common/srcrec/isrc,irec
      integer isrc,irec

        common/layertype/isfluid(NL),allsolid
        logical isfluid,allsolid

        common/shcontrol/ltop,lbot
        integer ltop, lbot

        logical proceed

c-----
c       get average  layer thickness
c-----
        deep = 0.0
        do 10 i=mmax,1,-1
            deep = deep + d(i)
   10   continue
        deep = deep/mmax
        k = 0
        nk = 0
        wvfu = 6.2831853*fu/vmin
        wvbm = omega/vmin
        dk = 6.2831853/xleng
        om=cmplx((omega),-(alpha))
        do 1000 js=1,mdpths
            do 1010 jr=1,mdpthr
                k = k + 1
                depth = abs(depths(js) - depthr(jr))
                if(depth.lt. deep)then
                    dpth = deep
                else
                    dpth = depth
                endif
c-----
c       check for wavenumber filtering. If done, do not use asymptotic
c-----
                if(cmin .le. 0.0)then
c-----
c       new logic for asymptotic as of 083192
c       If kh never gets large enough, the asymptotic
c       must be forced. Also the integration must at least go beyond
c       4 wvbm. Note that below we actually integrate in this case
c       over a square region in the omega-k plane rather than
c       the triangular when kh gets large 
c-----
                    wvmm = (5.0/dpth) + xfac*wvbm
                    wvzmx = wvmm * depth
                    if(wvzmx.gt.5.0)then
                        wave(k,1) = 6.0/depth
                        wave(k,2) = 2.5/depth
                        mk = wvmm / dk
                    else
                        wave(k,1)=20.0*wvfu
                        wave(k,2)=5.0*wvfu
                        if(wave(k,1)*depth .gt. 5.0)then
                            wave(k,1) = 6.0/depth
                            wave(k,2) = 
     1                      (5.0/dpth)+4.0*wvbm
                        endif
                        mk = wave(k,2) / dk
                    endif
                    if(mk .lt.1)mk = 1
                    wv = (mk-1)*dk + 0.218*dk
                    if(wv.gt.wave(k,1))then
                        gasymp(k) = .false.
                        if(fwhich.lt.0.0)then
                        fwhich=omega/6.2831853
                        endif
                    else
                        gasymp(k) = .true.
                    endif
                else
                    mk = wvmx /dk
                    gasymp(k) = .false.
                    if(fwhich.lt.0.0)then
                        fwhich=omega/6.2831853
                    endif
                endif
                if(ishank)gasymp(k) = .false.
c-----
c       define upper index of wavenumber
c-----
                mkup(k) = mk
                if(mk.gt.nk)nk = mk
c-----
c       now evaluate asymptotic coefficients, if appropriate
c-----
                if(gasymp(k))then
                    if(.not. lyrins)then
                        call modcpy(.false.)
                        call insert(depths(js))
                        call insert(depthr(jr))
                        call srclyr(depths(js),
     1                      lmaxs(js), dphs)
                        call srclyr(depthr(jr),
     1                      lmaxr(jr), dphr)
                        call dezero()
                        call equlck()
                    endif
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

                    om2 = om*om
                    wvn=cmplx((wave(k,1)),0.0)
                    wvn2 = wvn*wvn

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
                         ltop = 1
                         lbot = mmax 
                         if(allsolid)then
c-----
c                           get the bounds after setting defaults
c-----
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
                    call dopsv(OM,WVN,g1)
                    call dosh(om,wvn,g1)

C                   call evlmat(om,wvn,jbdrys,jbdryh,wvn2,om2)
C                   call rshof(g1,om,wvn,lmaxs(js),lmaxr(jr),wvn2,om2)
c-----
c                   adjust 12 and 15 to get avoid asymptotic problems
c-----
                    g1(12) = g1(12) * real(wvn)
                    g1(15) = g1(15) * real(wvn)

                    wvn=cmplx((wave(k,2)),0.0)
                    wvn2 = wvn*wvn

                    isrc = lmaxs(js)
                    irec = lmaxr(jr)
                    call dopsv(OM,WVN,g2)
                    call dosh(om,wvn,g2)

C                   call evlmat(om,wvn,jbdrys,jbdryh,wvn2,om2)
C                   call rshof(g2,om,wvn,lmaxs(js),lmaxr(jr),wvn2,om2)
c-----
c                   adjust 12 and 15 to get avoid asymptotic problems
c-----
                    g2(12) = g2(12) * real(wvn)
                    g2(15) = g2(15) * real(wvn)
c-----
c                    find asymptotic coefficients using these values
c-----
                    do 102 j=1,21
                            call solu(g1(j),g2(j),wave(k,1),
     1                  wave(k,2),depth,j,
     2                      aa(k,j),bb(k,j),cc(k,j) )
  102               continue
c-----
                else
                    do 1101 i=1,21
                        aa(k,i) = cmplx(0.0,0.0)
                        bb(k,i) = cmplx(0.0,0.0)
                        cc(k,i) = cmplx(0.0,0.0)
 1101               continue
                endif
 1010       continue
 1000   continue
        return
        end

        subroutine gasym(g,k,wvno,depth,jsrc)
c-----
c       remove asymptotic trend from integrands
c-----
        parameter(NSOURCE=100,NRECEIVER=100,NSR=100)
        common/kint4/aa,bb,cc
            complex aa(NSR,21),bb(NSR,21),cc(NSR,21)
        complex g(21)
        dimension jsrc(21)
        ex = exp(-wvno*depth)
        do 1000 j=1,21
            if(jsrc(j).eq.1)then
c           if(j.eq.12 .or. j.eq.15)then
c               g(j)=g(j) - ex*(aa(k,j)/wvno +bb(k,j))
c           else
                g(j)=g(j) - ex*(aa(k,j)+wvno*(bb(k,j)+
     1              wvno*(cc(k,j))))
c           endif
            endif
 1000   continue
c----
c       the fix for j=12, the P-SV horizontal force was made 8/26/92
c       this is OK since we never use k=0. Also, note that the
c       integrands involve J1 or a kJ0
c-----
        return
        end

        subroutine equlck()
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
        real bsh, rhosh, qbsh
        common/modlly/mmax
        common/updnsm/equalu(NL), equald(NL)
        logical equalu, equald
c-----
c       To avoid repeated computation, 
c       check to see if neighboring layers
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
     6          .and. qb(m).eq.qb(m-1) 
     2          .and. bsh(m).eq.bsh(m-1)
     4          .and. rhosh(m).eq.rhosh(m-1)
     5          .and. qbsh(m).eq.qbsh(m-1))then
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
     6          .and. qb(m).eq.qb(m+1) 
     2          .and. bsh(m).eq.bsh(m+1)
     4          .and. rhosh(m).eq.rhosh(m+1)
     5          .and. qbsh(m).eq.qbsh(m+1) )then
                equalu(m) = .true.
            else
                equalu(m) = .false.
            endif
  200   continue
        return
        end

        subroutine dezero()
        parameter (LER=0, LIN=5, LOT=6)
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
        real bsh, rhosh, qbsh
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

        subroutine srclay(depth,lmax,dph)
        parameter (LER=0, LIN=5, LOT=6)
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
        real bsh, rhosh, qbsh
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
        common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
        real bsh, rhosh, qbsh
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
    2   format(' Working model'/
     1      ' ',7x,'d',9x,'a',9x,'b',9x,'rho',6x,'1/qa',6x,'1/qb',
     2      7x,'bsh',4x,'1/qbsh')
    3   format(' ',4f10.3,2f10.6,2f10.3) 
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
            write(LOT,3)d(i),a(i),b(i),rho(i),qa(i),qb(i),
     1          bsh(i),rhosh(i)
            endif
   20   continue 
    5   format(' ',10x,3f10.3,2f10.6,2f10.3/' ') 
        write(LOT,5)a(mmax),b(mmax),rho(mmax),qa(mmax),qb(mmax) ,
     1          bsh(mmax),rhosh(mmax)
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
        implicit none
        integer irdwr, ierr
        integer BUFMAX
        parameter(BUFMAX=8192*64)
        common/buf/iptr,mxbuf,buffer(BUFMAX)
        integer iptr, mxbuf
        real buffer
        iptr = 1
        if(irdwr.eq.0)call getbuf(ierr)
        return
        end

        subroutine buflsh()
c-----
c       flush output buffer
c-----
        implicit none
        integer BUFMAX
        parameter(BUFMAX=8192*64)
        common/buf/iptr,mxbuf,buffer(BUFMAX)
        integer iptr, mxbuf
        real buffer
        integer ipt, i
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
        implicit none
        real x
        integer BUFMAX
        parameter(BUFMAX=8192*64)
        common/buf/iptr,mxbuf,buffer(BUFMAX)
        integer iptr, mxbuf
        real buffer
        buffer(iptr) = x
        iptr = iptr + 1
        if(iptr.gt.BUFMAX)call buflsh()
        return
        end

        subroutine getbuf(ierr)
c-----
c       read in file contents into buffer, taking care not to
c       read beyond the contents of the file
c-----
        implicit none
        integer ierr
        integer BUFMAX
        parameter(BUFMAX=8192*64)
        common/buf/iptr,mxbuf,buffer(BUFMAX)
        integer iptr, mxbuf
        real buffer
        integer i
c-----
c       ierr = 0 successful read
c            = 1 read error
c            = 2 end of file
c-----
        read(2,err=1000,end=2000)mxbuf,(buffer(i),i=1,mxbuf)
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
        implicit none
        real x
        integer ierr
        integer BUFMAX
        parameter(BUFMAX=8192*64)
        common/buf/iptr,mxbuf,buffer(BUFMAX)
        integer iptr, mxbuf
        real buffer
c-----
c       only yank in new data if actually required
c-----
        if(iptr.gt.mxbuf)call getbuf(ierr)
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
        common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
        real bsh, rhosh, qbsh
        common/modlly/mmax
        common/modelt/dt(NL),at(NL),bt(NL),rhot(NL),
     1      mmaxt,qat(NL),qbt(NL),etapt(NL),etast(NL),
     1      frefpt(NL),frefst(NL),bsht(NL),rhosht(NL),qbsht(NL)
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
                bsht(i) = bsh(i)
                rhosht(i) = rhosh(i)
                qbsht(i) = qbsh(i)
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
                bsh(i) = bsht(i)
                rhosh(i) = rhosht(i)
                qbsh(i) = qbsht(i)
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
        common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
        real bsh, rhosh, qbsh
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

        subroutine insert(dph)
        implicit none
        real dph
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
        real bsh, rhosh, qbsh

        real dp, dphh, hsave, dep
        integer m, ls
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
c           at surface and internally
c       However do put in a zero thickness layer 
c           at the base if necessary
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
                bsh(m+1) = bsh(m)
                qbsh(m+1) = qbsh(m)
                rhosh(m+1) = rhosh(m)
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

        subroutine frstar(r,hs,hr,mname,ipsvsh,time,pvel,svel,den,
     1      vsa, vsb, vsr, rayp, geom, tstar, dolock, dogeom)
c-----
c       r   R   Epicentral distance
c       hs  R   Source depth
c       hr  R   Receiver depth
c       mname   Ch*(*)  Name of model file
c       ipsvsh  I*4 1 - get P time
c               2 - get SV time
c               3 - get SH time
c               4 - get pP time
c               5 - get sP time
c       time    R   First arrival time
c       pvel    R   Velocity of P wave at receiver
c       svel    R   Velocity of S wave at receiver
c       den     R   Density at receiver
c       vsa R   P-wave velocity at source
c       vsb R   S-wave velocity at source
c       vsr R   Density at source
c       rayp R   Ray parameter in sec/km
c       geom R   geometrical spreading factor
c       tstar R   geometrical spreading factor
c       dolock L .true. apply locked mode which means to ignore the
c                bottom layer
c       dogeom L .true. compute geometrical spreading. Only do this for
c                teleseisms
c-----
        real r, hs, hr, time, pvel, svel, vsa, vsb, vsr
        real rayp, geom, tstar
        logical dolock, dogeom
        character mname*(*)
        integer ipsvsh
c-----
c-----
c       internal variables
c-----
        integer lmaxs, lmaxr, lmaxref
        real tm, rp(2),ts
        real dpdx, deginrad, sinis, sinir, cosir, cosis, sindeg
        real vs, vr, rs, rr

        
        common/earth/radius
        real radius

        common/depref/refdep
        real refdep
c-----
c       compute the travel time
c-----
        call fstarr(r,time,lmaxs, lmaxr, lmaxref,
     1      hs, hr, ipsvsh,iflsph, rayp,
     2      tstar, dolock,
     3      mname, varec, vbrec, rhorec,
     4      vasrc, vbsrc, rhosrc)

        svel = vbrec
        pvel = varec
        den  = rhorec
        vsa  = vasrc
        vsb  = vbsrc
        vsr  = rhosrc
c-----
C       compute the geometrical spreading
C
C       since a flat earth is always used, we use the flat layered
C       geometrical spreading from Officer
C       The geometrical spreading is dimensionless and gives the decrease in amplitude from
C       a distance of 1 km from the source to the receiver
C                        2                            2
C            ( rhos vs  a sin Is Vs                  d T      )
C       sqrt |  -------------------------------     -----     |
C            |                                         2      |
C            ( rhor vr sin DEG  cos Ir Rs Cos Is     dx       )
C
C       where p (sec/km) = dT/dx
C       a = radius of sphere about source - we use 1.0 km
C       Is= incident angle at source
C       Ir= incident angle at receiver
C       Rs= distance from center of Earth to source
C       Rr= distance from center of Earth to receiver
C       DEG=epicental distance in degrees
C       rhos and vs are the density and wave velocity at the source depth
c       rhor and vr are the density and wave velocity at the receiver depth
c
c       To get the dp/dx, we determine p at different distances, and then form
c       an interpolating polynomial
c
        if(dogeom)then
              call fstarr(r-500,tm,lmaxs, lmaxr, lmaxref,
     1            hs, hr, ipsvsh,iflsph, rp(1),
     2            ts, dolock,
     3            mname, varec, vbrec, rhorec,
     4            vasrc, vbsrc, rhosrc)

              call fstarr(r+500,tm,lmaxs, lmaxr, lmaxref,
     1            hs, hr, ipsvsh,iflsph, rp(2),
     2            ts, dolock,
     3            mname, varec, vbrec, rhorec,
     4            vasrc, vbsrc, rhosrc)


              dpdx = abs(rp(1) - rp(2))/(500.0 - ( - 500.0))
              deginrad = r/radius
              sindeg = sin(deginrad)
  
              if(ipsvsh.eq.1)then
                  vs = vsa
                  vr = pvel
                  rs = vsr
                  rr = den
              else if(ipsvsh.eq.2)then
                  vs = vsb
                  vr = svel
                  rs = vsr
                  rr = den
              else if(ipsvsh.eq.3)then
                  vs = vsb
                  vr = svel
                  rs = vsr
                  rr = den
              else if(ipsvsh.eq.4)then
                  vs = vsa
                  vr = pvel
                  rs = vsr
                  rr = den
              else if(ipsvsh.eq.5)then
                  vs = vsb
                  vr = pvel
                  rs = vsr
                  rr = den
              endif
              sinis = rayp*vs
              cosis = sqrt(1.0-sinis*sinis)
              sinir = rayp*vr
              cosir = sqrt(1.0-sinir*sinir)

              fac = (rs*vs*sinis*vs*dpdx)/
     1           (rr*vr*sindeg*cosir*(radius-hs)*cosis)
              geom = sqrt(abs(fac))
      else
c-----
c     default to have a defined return value
c-----
              geom = 1.0
      endif
              

c-----
        return
        end

        subroutine fstarr(dist,tfirst,lmaxs,lmaxr,lmaxref,
     1      hs,hr,ipsvsh,iflsph, rayp,
     2      tstar, dolock,
     3      mname, varec, vbrec, rhorec,
     4      vasrc, vbsrc, rhosrc)
c-----
c       given a distance, the source depth, receiver depth,
c       get time of first arrival of P
c-----
c       dist    R   - distance
c       tfirst  R   - first arrival time
c       mmax    I*4 - number of layers in model
c       lmaxs   I*4 - layer index for source
c       lmaxr   I*4 - layer index for receiver
c       lmaxref I*4 - layer index for reference depth,
c                     used only for pP and sS
c       hs      R   - depth of source
c       hs      R   - depth of receiver
c       ipsvsh  I*4 1 - get P time
c               2 - get SV time
c               3 - get SH time
c               4 - get pP time
c               5 - get sP time
c       iflsph  I*4 - 0 flat earth model
c                 1 spherical earth model
c       rayp    R   - ray parameter in sec/km
c       tstar   R   T* operator
c       dolock L .true. apply locked mode which means to ignore the
c                bottom layer
c       mname   Ch* - name of the velocity model
c       varec   R - P velocity at receiver (untransformed)
c       vbrec   R - S velocity at receiver (untransformed)
c       rhorec  R - Density at receiver (untransformed)
c       vasrc   R - P velocity at source (untransformed)
c       vbsrc   R - S velocity at source (untransformed)
c       rhosrc  R - Density at source (untransformed)
c-----
c       since this routine is to be used for omega-k,
c       we will approximate the direct arrival
c
c       18 JAN 2008 - everything is straightforward. The addition of
c          the request for pP and sP changes the logic in that
c          the direct arrival is ignored, and that the upgoing refraction 
c          from the source is ignored. We handle this by just setting
c          a very large tfirst before trying to do the modified 
c          downward path refraction to avoid another level of
c          if/then/else/endif
c       24 MAR 2008 - modified to place the model read into this
c          routine instead of in frstar
c-----
        real dist, tfirst, dphs, dphr, hs, hr, depths, depthr
        real rayp, tstar
        real varec, vbrec, rhorec
        real vasrc, vbsrc, rhosrc
        logical dolock
        integer lmaxs, lmaxr, lmaxref
        character mname*(*)

        integer NL
        parameter (NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmmax
        integer mmmax
        common/depref/refdep
        real refdep
        common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
        real bsh, rhosh, qbsh

        integer mmax

        integer iunit, iiso, iflsph, idimen, icnvel, ierr
        character title*80

        real v(NL), h(NL), qi(NL)

        real*8 c, s, t, x, p, tint, dxdp, vel, pnew, pupper
        real*8 ts
        real*8 sumx, sumt
        logical ext

        real tds, tdr
        common/earth/radius
        real radius

c-----
c       first read in the model and determine the medium parameters at the
c       source and receiver depths
c-----
c       get the earth model
c-----
        inquire(file=mname,exist=ext)
        if(.not. ext) call usage('Model file does not exist')
        l = lgstr(mname)

                call getmod(1,mname,mmax,title,iunit,iiso,iflsph,
     1          idimen,icnvel,ierr,.false.)
        mmmax = mmax
        if(ierr .lt. 0)return
                call adomod()
c-----
c       insert the source and receiver depths into the model
c       placing the source and receiver on a layer boundary
c-----
        call insert(hs+refdep)
        call insert(hr+refdep)
        call insert(   refdep)

c-----
c       get the layer in which the source lies
c-----
        call srclyr(hs+refdep, lmaxs, dphs)
        call srclyr(hr+refdep, lmaxr, dphr)
        call srclyr(   refdep, lmaxref, dphref)
c-----
c       get the medium parameters at the source and reciever depths
c-----
        varec = a(lmaxr)
        vbrec = b(lmaxr)
        rhorec = rho(lmaxr)
        vasrc = a(lmaxs)
        vbsrc = b(lmaxs)
        rhosrc = rho(lmaxs)
c-----
c       prepare for the computations
c-----
        depths = hs + refdep
        depthr = hr + refdep

c-----
c       set up default
c-----
        tfirst = 1.0e+30
c-----
c       special case for locked mode
c-----
        if(dolock)then
            mmax = mmmax -1
        else
            mmax = mmmax
        endif

c-----
c       get specifics about upward and downward distances
c       with a layer. We need his to define ray paths
c       We will also use the fact that the source/receiver are
c       on layer boundaries
c
c       lmn = layer number of shallowest of source/receiver
c       lmx = layer number of deepest    of source/receiver
c-----
        lmn = min(lmaxs,lmaxr)
        lmx = max(lmaxs,lmaxr)

c-----
c       perform spherical -> flat earth model transformation
c-----
        if(iflsph.ne.0)then
            call adosph()
            tds = radius*alog(radius/(radius-hs))
            tdr = radius*alog(radius/(radius-hr))
        else
            tds = depths
            tdr = depthr
        endif
c-----
c       now fill in velocity array according to desired first arrival
c       for SH there can be no water layer
c       for SV can be a water layer
c       Also define the Q for the T* analysis. Note we define
c        eventually q = 1/Q based on whether the given Q > or < 1
c-----
        do 100 i=1,mmax
            if(ipsvsh.eq.1)then
                v(i) = a(i)
                qi(i) = qa(i)
            else if(ipsvsh.eq.2)then
                v(i) = b(i)
                qi(i) = qb(i)
                if(b(i).le.0.001)then
                    v(i) = a(i)
                    qi(i) = qa(i)
                endif
            else if(ipsvsh.eq.3)then
                v(i) = bsh(i)
                qi(i) = qbsh(i)
            else if(ipsvsh.eq.4)then
                v(i) = a(i)
                qi(i) = qa(i)
            else if(ipsvsh.eq.5)then
                v(i) = a(i)
                qi(i) = qa(i)
            endif
            if(qi(i) .gt. 1.0)then
                qi(i) = 1.0 / qi(i)
            endif
            h(i) = d(i)
 100    continue
c-----
c       For the computations we look at four cases
c       1) direct path between source and receiver 
c          a) source and receiver in the same layer
c          b) source and receiver in different layers
c       2) refracted arrivals       
c          a) path is downward from source and then up to
c             receiver
c          b) path is upward from the source and then down to
c             receiver
c          This recognized the possibility that velocity does
c          not increase uniformly with depth
c-----
                    
c-----
c       direct arrival source/receiver at same layer
c-----
        if(v(lmaxs).eq.0.0)return
        if(v(lmaxr).eq.0.0)return
        if(lmaxs .eq. lmaxr)then
            tfirst = sqrt(dist**2 + abs(tds - tdr)**2)/
     1          v(lmaxs)
            rayp = (dist/sqrt(dist**2 + abs(tds - tdr)**2))/
     1          v(lmaxs)
            tstar = tfirst*qi(lmaxs)
        else
c-----
c       direct arrival source/receiver in different layers
c-----
c       Newton Iteration for direct arrival source/receiver at
c           different depths
c           
c           x = SUM h  tan theta
c                    i          i
c
c           t = SUM h  / v  cos theta
c                    i    i          i
c                                                          2 2
c       where sin theta  = p V  , cos theta  = sqrt ( 1 - p V )
c                      i      i                              i
c       and p is the ray parameter bounded by [0, 1/V  ] where V
c                                                    sr         sr
c       is the wave velocity at the starting point of the ray. 
c       Since the ray must also reach the receiver, we consider
c       that possibility too. The upper bound is MIN ( 1/Vs, 1/Vr)
c       Also we test for a real ray path, between source and receiver
c
c       Because source/receiver at top of a layer boundary, we have
c
c           -----------X----------
c           h(lmn)      \
c           ----------------------
c                      ....
c           ----------------------
c           h(lmx-1)        \
c                            \
c           ------------------X---
c            
c-----
            ps = 1.0/v(lmaxs)
            pr = 1.0/v(lmaxr)
            if(ps.lt.pr)then
                pupper = ps
            else
                pupper = pr
            endif
            do 1000 l=lmn,lmx
                if(v(l).eq.0.0)return
                p = 1.0/v(l)
                if(p.lt.pupper)pupper = p
 1000       continue
            p = 0.5d+00  * pupper
            do 1111 iter=1,20
                x = 0.0d+00
                t = 0.0d+00
                ts = 0.0d+00
                tint = 0.0d+00
                dxdp = 0.0d+00
                do 1500 l=lmn,lmx - 1
                    vel = dble(v(l))
                    s = p*vel
                    c = dsqrt(1.0d+00 - s*s)
                    t = t + dble(h(l)) /(vel*c)
                    x = x + dble(h(l)) * s / c
                    dxdp  = dxdp + dble(h(l)) *
     1                  vel/(c*c*c)
                    tint = tint + dble(h(l)) * c / vel
                    ts = ts + qi(l) * dble(h(l))/(c*vel)
 1500           continue
                pnew = p - (x-dble(dist))/dxdp
c-----
c       safety - we must have a real ray, with upper bound
c       of  min[ 1/v(src), 1/v(rec)]
c-----  
                if(pnew .gt. pupper)then
                    if(iter.lt.20)then
                        pnew = 0.99d+00 * pupper
                    else
c-----
c       this is propably working like a refraction, so stop iterations
c-----  
                        t = tint + p * (dist)
                        go to 1112
                    endif
                endif
                p = pnew
 1111       continue
 1112       continue
            tfirst = t
            rayp = p
            tstar = ts
        endif
c-----
c       now proceed through the possible refracted arrivals
c       considering first upward rays from the source
c-----  
        if(lmn.gt.1)then
        do 3020 m=1,lmn-1
c-----
c       m is the refracting layer
c
c       get velocity of refractor
c-----
            vel = v(m)
            if(v(m).eq.0.0)return
            p = 1.0/vel
c-----
c
c           --------------------------------
c           h(1)
c           --------------------------------
c                      ....
c           --------------------------------
c           h(m)
c           ----------------...-------------
c           h(m+1)         /   \
c           --------------------------------
c                         /     \
c                      ....
c           --------------------------------
c           h(lmn-1)              \
c           -----------------------X--------
c               
c           h(lmn)     /    
c           --------------------------------
c                      ....
c           --------------------------------
c           h(lmx-1) /
c           --------X-----------------------
c
c       safety check, velocity at source or receiver must be less than
c       refraction velocity
c-----
        if(v(lmn).ge.vel)go to 3020
        if(v(lmx).ge.vel)go to 3020
c-----
c       single leg
c-----
        sumt = 0.0
        sumx = 0.0
        ts = 0.0
            do 3021 l=1,lmx-1,lmn
                if(v(l).gt.vel)go to 3020
                cs = sqrt(abs(1.0 - p*p*v(l)*v(l)))
                sumt = sumt + h(l)*cs/v(l)
                sumx = sumx + h(l)*p/cs
                ts = ts + qi(l)*h(l)/(cs * v(l))
 3021       continue
            do 3022 l=m+1,lmn-1
                if(v(l).gt.vel)go to 3020
                cs = sqrt(abs(1.0 - p*p*v(l)*v(l)))
                sumt = sumt + 2.0*h(l)*cs/v(l)
                sumx = sumx + 2.0*h(l)*p/cs
                ts = ts + 2.0*qi(l)*h(l)/(cs * v(l))
 3022       continue
            tint = sumt
            tt = tint + dist / vel
            ts = ts + qi(m)*(dist-sumx)/v(m)
            if(tt .lt. tfirst .and. dist.ge.sumx)then
                  tfirst = tt
                  rayp = p
                 tstar = ts
            endif
 3020       continue
        endif
c-----
c       For the special case of the depth phases, ignore previous
c       first arrival times
c-----
        if(ipsvsh.eq.4 .or. ipsvsh.eq.5)then
             tfirst = 1.0e+30
        endif
c-----
c       now proceed through the possible refracted arrivals
c       considering first downward rays from the source
c
c       We start considering the deepest point since we place
c       a source/receiver position just below a layer boundary
c       and thus should consider a horizontal ray
c
c       The refraction is accepted only if the desired distance >
c       first refraction from the source - this puts physics in the problem
c           
c           x = SUM h  tan theta
c                    i          i
c
c           t = SUM h  cos theta / V
c                    i          i   i
c                                                          2 2
c       where sin theta  = p V  , cos theta  = sqrt ( 1 - p V )
c                      i      i                              i
c       For the T* computation we need to follow the path, e.g.,
c       SUM h qi / ( cos theta  / V ) + qi (dist -  SUM h tan theta / V )/V
c            i  i             i    i      i              i         i   i   r
c-----  
        do 2020 m=lmx+1, mmax
c-----
c       m is the refracting layer
c
c       get velocity of refractor
c-----
            vel = v(m)
            if(v(m).eq.0.0)return
            p = 1.0/vel
c-----
c
c           -----------X--------------------
c           h(lmn)      \
c           --------------------------------
c                      ....
c           --------------------------------
c           h(lmx-1)        \             
c                            \           
c           ------------------X--------X----
c           h(lmx)             \       /
c           --------------------\-----/-----
c                      ....      \   /
c           ----------------------...-------
c           h(m)
c
c-----
c       safety check, velocity at source or receiver must be less than
c       refraction velocity
c-----
        if(v(lmn).ge.vel)go to 2020
        if(v(lmx).ge.vel)go to 2020
c-----
c       single leg
c-----
        sumx = 0.0
        sumt = 0.0
        ts = 0.0
c-----
c       special case for depth phases
c-----
            if(ipsvsh.eq.4)then
c-----
c               pP
c-----
                  do  l=lmaxref,lmaxs - 1
                      if(a(l).gt.vel)go to 2020
                      cs = sqrt(abs(1.0 - p*p*a(l)*a(l)))
                      sumt = sumt + 2.*h(l)*cs/a(l)
                      sumx = sumx + 2.*h(l)*p*a(l)/cs
                      if(qa(l).gt.1.0)qa(l) = 1.0/qa(l)
                      ts = ts + 2.*qa(l)*h(l)/(cs * a(l))
                  enddo
            else if(ipsvsh.eq.5)then
c-----
c               sP
c-----
                  do  l=lmaxref,lmaxs - 1
                      if(a(l).gt.vel)go to 2020
                      if(b(l).gt.vel)go to 2020
                      csa = sqrt(abs(1.0 - p*p*a(l)*a(l)))
                      csb = sqrt(abs(1.0 - p*p*b(l)*b(l)))
                      sumt = sumt + h(l)*csa/a(l)
     1                        +h(l)*csb/b(l)
                      sumx = sumx + 2.*h(l)*p*a(l)/csa
                      if(qa(l).gt.1.0)qa(l) = 1.0/qa(l)
                      if(qb(l).gt.1.0)qb(l) = 1.0/qb(l)
                      ts = ts + qa(l)*h(l)/(csa * a(l))
     1                        + qb(l)*h(l)/(csb * b(l))
                  enddo
            endif
c-----
c       continue
c-----
            do 2021 l=lmn,lmx - 1
                if(v(l).gt.vel)go to 2020
                if(v(l).eq.0.0)return
                cs = sqrt(abs(1.0 - p*p*v(l)*v(l)))
                sumt = sumt + h(l)*cs/v(l)
                sumx = sumx + h(l)*p*v(l)/cs
                ts = ts + qi(l)*h(l)/(cs * v(l))
 2021       continue
c-----
c       double leg
c-----
            do 2022 l=lmx,m-1
                if(v(l).gt.vel)go to 2020
                if(v(l).eq.0.0)return
                cs = sqrt(abs(1.0 - p*p*v(l)*v(l)))
                sumt = sumt + 2.0*h(l)*cs/v(l)
                sumx = sumx + 2.0*h(l)*p*v(l)/cs
                ts = ts + 2.*qi(l)*h(l)/(cs * v(l))
 2022       continue
            tint = sumt
            tt = tint + dist / vel
            ts = ts + qi(m)*(dist-sumx)/vel
            if(tt .lt. tfirst .and. dist.ge.sumx)then
                 tfirst = tt
                 rayp = p
                 tstar = ts
            endif
 2020       continue
             if(tfirst .eq. 1.0e+30)then
                tfirst = -12345.
                rayp   = -12345.
                tstar  = -12345.
             endif
        return
        end

        subroutine adosph()
c-----
c       Transform spherical earth to flat earth
c
c       Schwab, F. A., and L. Knopoff (1972). 
c           Fast surface wave and free
c       mode computations, in  
c           Methods in Computational Physics, Volume 11,
c       Seismology: Surface Waves and Earth Oscillations,  
c           B. A. Bolt (ed),
c       Academic Press, New York
c
c       Love Wave Equations  44, 45 , 41 pp 112-113
c       Rayleigh Wave Equations 102, 108, 109 pp 142, 144
c
c
c-----
c       mmax    I*4 number of layers
c       ipsvsh  I*4     1 - get P time
c                       2 - get SV time
c                       3 - get SH time
c-----
        implicit none
        integer NL
        parameter (NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/depref/refdep
        real refdep
        common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
        real bsh, rhosh, qbsh

        double precision z0,z1,r0,r1,dr,ar,tmp
        real btmp, qbtmp, rhotmp
        integer i

        common/earth/radius
        real radius

        ar=radius
        dr=0.0d0
        r0=ar + refdep
        d(mmax)=1.0
        do 10 i=1,mmax
            r1=r0-dble(d(i))
            z0=ar*dlog(ar/r0)
            z1=ar*dlog(ar/r1)
            d(i)=z1-z0
c-----
c        attempt 7 15 2007 - use standard rule but at mid layer depth as per DGH
c-----
            TMP=(ar+ar)/(r0+r1)

            a(i)=a(i)*tmp
            btmp = b(i)
            b(i)=btmp*tmp
            bsh(i)=btmp*tmp
            qbtmp = qb(i)
            qbsh(i)=qbtmp
            rhotmp=rho(i)
            rhosh(i) = rhotmp * tmp **(-5.0)
            rho(i) = rhotmp * tmp **(-2.275)
            r0 = r1
   10   continue
        d(mmax)=0.0
        return
        end

        subroutine adomod()
c-----
c       just fill the rhosh, bsh and qbsh arrays 
c-----
        implicit none
        integer NL
        parameter (NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/depref/refdep
        real refdep
        common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
        real bsh, rhosh, qbsh

        integer i

        do  i=1,mmax
            bsh(i)=b(i)
            qbsh(i)=qb(i)
            rhosh(i) = rho(i) 
        enddo
        return
        end

        subroutine srclyr(depth,lmax,dph)
        implicit none
        real depth, dph
        integer lmax
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
        real bsh, rhosh, qbsh

        integer m
        real dep
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

