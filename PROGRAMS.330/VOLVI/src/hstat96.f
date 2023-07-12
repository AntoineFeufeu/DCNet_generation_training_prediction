        program hstat96
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME VI                                                      c
c                                                                     c
c      PROGRAM: HSTAT96                                               c
c                                                                     c
c      COPYRIGHT 2011                                                 c
c      R. B. Herrmann                                                 c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63013                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
c       CHANGES
c       02 JUN 2011 - created this program based on hspec96
c         This program computes permanent deformation for
c         a step change in moment or force.  Perhaps this
c         will be made a part of hspec96. Looking forward to,
c         the real double precision will be maintained
c       01 AUG 2013 - add -V flag to force  verbose text output
c       21 OCT 2018 - correct output flags for file96 for the case
c         that -EQEX in hspec96.  The -ALL worked correctly
c       11 APR 2022 - updated the setlim subroutine 
c
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
        real*8 har(NL,4,4), dar(NL,5,5), hsr(2,5), gbr(2,5), 
     1      hal(NL,2,2), hsl(2,2), gbl(2,2)
        real*8 hex(NL), lex(NL), dex(NL), hexw(NL)
        common/hamat/har
        common/damat/dar
        common/hsrfr/hsr
        common/gbrfr/gbr
        common/hlmat/hal
        common/hsrfl/hsl
        common/gbrfl/gbl
        common/hexex/hex
        common/hexexw/hexw
        common/dexex/dex
        common/lexex/lex 
        common/water/iwater(NL),iwats(2),iwatb(2)
        common/updnsm/equalu(NL), equald(NL)
        logical equalu, equald
        common/lyrctl/lyrins
        logical lyrins
        common/rlimit/rlim
        real*4 rlim
c-----c
        character*7 istat3
        logical ixst3
        parameter(NDST=500)
        real*4 r(NDST), tshift(NDST), vred(NDST)
        real smm(NSR,21)
        real stmp, sdata
        common/c/cmax,c1,c2,cmin
        common/frlim/fl,fu,df,fwhich
        common/cntrl/ishank,hnkarg,dstcor,dokjar,docausal
        common/kint1/doasym
        logical doasym
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

         common/num/numer(21)
         real numer


        common/ihdr96/ ftype, iftype, iobsyn, itmfrq, iunit, idecon,
     1      keyear, kemon, keday, kehour, kemin,
     2      ksyear, ksmon, ksday, kshour, ksmin,
     3      jjsrc, junit
        integer*4 ftype, iftype, iobsyn, itmfrq, iunit, idecon
        integer*4 keyear, kemon, keday, kehour, kemin
        integer*4 ksyear, ksmon, ksday, kshour, ksmin
        integer*4 jjsrc(21), junit

        common/rhdr96/ esec, distkm, distdg, evstaz, stevaz,
     1      evlat, evlon, evdep, stlat, stlon, stelev,
     2      TP, TSV, TSH, SA, SC, SF, SL, SN, SR
        real*4 esec, distkm, distdg, evstaz, stevaz
        real*4 evlat, evlon, evdep 
        real*4 stlat, stlon, stelev 
        REAL*4 TP, TSV, TSH
        REAL*4 SA, SC, SF, SL, SN, SR

        common/chdr96/ stname, cfilt, cpulse, ccomnt
        character*8 stname*8, cfilt*80, cpulse*80, ccomnt*80    

        character stcomp*8
        real*4 cmpaz, cmpinc, cmpdt
        real*4 ssec

        logical outtxt
        common/verby/verbose
        logical verbose

c-----
c       ksrc = temporary array for jsrc for output
c       here jsrc != 0 reflects source information
c       if receiver is not in a fluid DO NOT output pressure field
c-----
        integer ksrc(21)
c-----
c       lsrc maps jsrc to output Green s functions. e.g., if
c       jsrc(8) = radial explosion term, but in final output it
c       ccupies position 10, or jsrc(lsrc(10)) = computed
c-----
        integer*4 lsrc(21)
        integer*4 iszrt(21)
        character*8 ost(21)
        data iszrt/1,4, 1,4,5, 1,4,5, 1,4, 1,4, 1,4,5,1,1,1,1,1,1/
        data lsrc/1,2,3,4,13,5,6,14,7,8,9,10,11,12,15,16,17,18,19,20,21/
        data ost/'ZDD     ', 'RDD     ', 'ZDS     ', 'RDS     ',
     1       'TDS     ', 'ZSS     ', 'RSS     ', 'TSS     ',
     2       'ZEX     ', 'REX     ', 'ZVF     ', 'RVF     ',
     3       'ZHF     ', 'RHF     ', 'THF     ', 'PEX     ',
     4       'PDD     ', 'PDS     ', 'PSS     ', 'PVF     ',
     5       'PHF     '/

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
        call gcmdln(outtxt)
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
c-----
c       UNIX output - no carriage control
c-----
 2021   format('depths =',f20.6,' (',f20.6,')')
 2031   format('depthr =',f14.6,' (',f20.6,')')
 2041   format('     r =',f14.6,' tshift=',f10.2,' vred=',f10.2)
    3   format(8f10.5)
    4   format('frequencies for which response computed     ')
    5   format('alpha =',f10.5,5x,'dt =',f10.3)

        if(verbose)then
           write(LOT,*)'SOURCE DEPTH IN WORKING AND ORIGINAL MODEL (',
     1          mdpths,')'
           do i=1,mdpths
               write(LOT,2021)depths(i), depthssv(i)
           enddo
           write(LOT,*)'RECEIVER DEPTH IN WORKING AND ORIGINAL MODEL ('
     1          ,mdpthr,')'
           do i=1,mdpthr
               write(LOT,2031)depthr(i), depthrsv(i)
           enddo
           write(LOT,*)'RECEIVER DISTANCES (',ndist,')'
           do i=1,ndist
               write(LOT,2041)r(i), tshift(i), vred(i)
           enddo
        endif
c-----
c     open the output file file96
c-----
      open(8,file='file96',status='unknown',form='formatted',
     1 access='sequential')
      rewind 8
      do i=1,16
         if(jsrc(lsrc(i)).eq.1)then
             jjsrc(i) = 1
         else
             jjsrc(i) = 0
         endif
      enddo
      jjsrc(16) = 0
c-----
c     process everything
c-----

      do jd=1,ndist
         k = 0
         do is=1,mdpths
c-----
c           for this source depth we need the velocity and
c           density to compute the SA, SC< SF, SL, SN, SR
c-----  
            call srcpar(is)
            do ir=1,mdpthr
               cpulse='hstat96'
               iunit = 2
               iftype = 16
               iobsyn = 2
               itmfrq = 1
               cfilt = 'hstat96'
               keyear = 0
               kemon = 0
               keday = 0
               kehour = 0
               kemin = 0
               esec = 0.0
               evlat = -12345
               evlon = -12345
               evdep = depthssv(1)
               ksyear = 0
               ksmon = 0
               ksday = 0
               kshour = 0
               ksmin = 0
               ssec = 0.0
               tp = -12345.
               tsv = -12345.
               tsh = -12345.
               stname = 'GRN21'
               stlat  = -12345
               stlon = -12345
               stelev = depthrsv(1)
               distkm = r(jd)
               distdg = r(jd)/111.195
               evstaz = 0.0
               stevaz = 180.0
               lmnm = lgstr(mname)
               ccomnt = mname(1:lmnm)
               cmpdt = 1.0
               call wrhd96(8,nerr)
               k = k + 1
               call numerical(r(jd),is,ir)
               do  jj=1,15
c-----
c                 vertical
c-----
                  if(iszrt(jj).eq.1)then
                     cmpinc = -90.0
                     cmpaz  =   0.0
c-----
c                 radial
c-----
                  else if(iszrt(jj).eq.4)then
                     cmpinc = 0.0
                     cmpaz  =   0.0
c-----
c                    transverse
c-----
                  else if(iszrt(jj).eq.5)then
                     cmpinc = 0.0
                     cmpaz  =  90.0
                  endif
                  if(jsrc(lsrc(jj)).eq.1)then
                     if(verbose)then
                        write(LOT,*)'HSTATIC:',
     1                     jj,ost(jj)(1:3),lsrc(jj),numer(jj)
                     endif
                     npts = 1
                     call wrtr96(8,ost(jj)(1:3),cmpinc,cmpaz,cmpdt,npts,
     1                  ksyear, ksmon, ksday, kshour, ksmin, ssec,
     2                  numer(jj), nerr, 1)
                   endif
              enddo
            enddo
         enddo
      enddo
      close(8)
 9999 continue
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
        common/verby/verbose
        logical verbose
c-----
c       UNIX FORTRAN - NO CARRIAGE CONTROL
c-----
   21   format(11i5)
   22   format(6e11.4)
   24   format('XLENG=',e15.7,' XFAC=',e15.7)
   30   format(2e15.7/1x,3i10)
 4012       format('WAVENUMBER FILTERING bounded in phase velocites'/
     1          '[cmax,c1,c2,cmin]=','[', f10.3, ',', 
     2          f10.3, ',', f10.3, ',', f10.3,']' /
     3          '( -1.0 means 0.0 for cmin and infinity for cmax)')
 4013   format('WAVENUMBER FILTERING NOT DONE')
c-----
c       MSDOS MICROSOFT FORTRAN OUTPUT - CARRIAGE CONTROL
c-----
c   21  format(1x,11i5)
c   22  format(1x,6e11.4)
c   24  format(1x,'XLENG=',e15.7,' XFAC=',e15.7)
c   30  format(1x,2e15.7/1x,3i10)
c 4012      format(1x,'WAVENUMBER FILTERING bounded in phase velocites'/
c     1         1x,'[cmax,c1,c2,cmin]=',
c     2         1x,'[', f10.3, ',', f10.3, ',', f10.3, ',', f10.3,']' /
c     3         1x,'( -1.0 means 0.0 for cmin and infinity for cmax)')
c 4013  format(1x,'WAVENUMBER FILTERING NOT DONE')
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
        if(verbose)then
             write(LOT,30) alpha,delt,n,n1,n2
        endif
        df = 1.0/(n*delt)
        fl = (n1-1)*df
        fu = (n2-1)*df
c-----
c       Specify desired output Green s functions
c-----
        read(lun,*)ieqex
        if(ieqex.gt.6 .or. ieqex.lt.0)ieqex = 2
        if(verbose)then
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
        endif
c-----
c       provide names for output Green s functions in order of output
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
        if(verbose)then
           write(LOT,21) jbdry
        endif
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
        if(verbose)then
           if(jbdrys.eq.1)then
               write(LOT,*)' TOP  OF MODEL IS FREE SURFACE  '
           else if(jbdrys.eq.0)then
               write(LOT,*)' TOP  OF MODEL IS HALFSPACE WITH',
     1             ' PROPERTIES OF FIRST LAYER'
           else if(jbdrys.eq.-1)then
               write(LOT,*)' TOP  OF MODEL IS RIGID'
           endif
           if(jbdryh.eq.0)then
               write(LOT,*)' BASE OF MODEL IS HALFSPACE WITH',
     1             ' PROPERTIES OF BOTTOM LAYER'
           else if(jbdryh.eq.-1)then
               write(LOT,*)' BASE OF MODEL IS RIGID'
           else if(jbdryh.eq.1)then
               write(LOT,*)' BASE OF MODEL IS FREE'
           endif
        endif


c-----
c       read in the earth model name
c-----
        read(lun,'(a)')mname
        lmnm = lgstr(mname)
        if(verbose)then
           write(LOT,*)mname(1:lmnm)
        endif
c-----
c       read in the earth model in model96 format
c-----
        inquire(file=mname,exist=ext)
        if(.not. ext)then
            write(LER,*)'Model file not located'
            stop
        endif
        call getmod(1,mname,mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,verbose)
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
        if(verbose)then
           write(LOT,24)xleng, xfac
        endif
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
                depths(i) = radius*alog(radius/(radius-depths(i)))
            endif
 1191       continue
            do 1192 i=1,mdpthr
            depthrsv(i) = depthr(i)
            if(iflsph.ne.0)then
                depthr(i) = radius*alog(radius/(radius-depthr(i)))
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
           if(verbose)then
               write(LOT,4012)cmax,c1,c2,cmin
           endif
        goto 4011
 4010   continue
            cmax = -1.0
            c1 = -1.0
            c2 = -1.0
            cmin = -1.0
        if(verbose)then
           write(LOT,4013)
        endif
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
CRBH        if(mdpths+mdpthr .gt. 2*mmax )then
            lyrins = .false.
            if(verbose)write(LOT,*)' LAYER INSERTION NOT DONE'
CRBH        else
CRBH            lyrins = .true. 
CRBH            if(verbose)write(LOT,*)' LAYER INSERTION DONE'
CRBH            do 2108 i=1,mdpths
CRBH                call insert(depths(i))
CRBH 2108       continue
CRBH            do 2109 i=1,mdpthr
CRBH                call insert(depthr(i))
CRBH 2109       continue
CRBH            call dezero()
CRBHc-----
CRBHc           check whether neighboring layers are identical
CRBHc           to avoid redundant evaluation
CRBHc----
CRBH            call equlck()
CRBH        endif
c-----
c       verify the new model parameters
c-----
        if(verbose)then
           write(LOT,*)'mmax=',mmax
           write(LOT,22)(d(i),a(i),b(i),rho(i),qa(i),qb(i),i=1,mmax)
        endif
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

      subroutine svar(k,z,ex,exa,c,s,skz,ckz)
      implicit none
      real*8 k, z, ex, exa, c, s, skz, ckz
      common/sovrflw/a0,c2,s2,sc,kz,k2z2
      real*8 a0,c2,s2,sc,kz,k2z2
c-----
c     here compute the functions required for shska and sdnka
c     for the shska we factor out an exp(kz) which is saved as
c     the parameter ex, e.g.,
c
c     ex == kz
c     cosh(kz) = exp(kz) * [ 1 + exp (-2kz) ]/2 == exp(kz) c 
c     sinh(kz) = exp(kz) * [ 1 - exp (-2kz) ]/2 == exp(kz) s 
c     skz = s * kz
c     ckz = c * kz
c     
c     Note that for shska the product kz does only appears in combination 
c     with cosh or sinh
c
c     For the compund matrix we must worry about the combination of double
c     exponential together with non exponentials. For example,
c     
c                            2      2                            2 2 2
c     A33 = 1 * (mu/lam+2 mu)   sinh kz - (lam + mu)/(lam + s mu) k z
c     
c     For stability reasons, we write this as
c   
c                                      2                             2  
c     A33 = exp(2kz)[a0 * (mu/lam+2 mu)   s2 - (lam + mu)/(lam + s mu) k2z2]
c     
c     where a0 = exp(-2kz)
c     s2 = s*s          2 2
c     k2z2 = exp(-2kz)*k z
c     exa = 2*ex
c-----
c     internal variables
c     dfac used to avoid underflow in exp(-2kz)
c-----
      real*8 dfac

      ex = k*z
      exa = ex + ex

      if(ex .gt. 40.0d+00)then
c-----
c          exp(-2ka) < 1e-13
c-----
           dfac = 0.0d+00
      else
           dfac = dexp(-ex-ex)
      endif
c-----
c     terms used by shska
c-----
      c = ( 1. + dfac)/2.
      s = ( 1. - dfac)/2.

C      ex = 0
C      exa = 0
C      dfac = 1
C      c = (exp(k*z) + exp(-k*z))/2.
C      s = (exp(k*z) - exp(-k*z))/2.
      skz = k*z*s
      ckz = k*z*c
c-----
c     terms used by sdnka
c-----
      a0 = dfac
      c2 = c*c
      s2 = s*s
      sc = s*c
      kz = dfac*k*z
      k2z2 = dfac*k*z*k*z

      return
      end

      subroutine shska(saa,k,z,c,s,skz,ckz,lam,mu)
      implicit none
      real*8 saa(4,4),k,z,c,s,skz,ckz,lam,mu
c-----
c     The propagator matrix is actually
c       A = exp(kz) AA
c     because of the factorization within lvar. Both A and AA
c     have a symmetry, e.g.,
c
c     
c     | a11   a12   a13   a14 |
c     | a21   a22   a23  -a13 |
c     | a31   a32   a33  -a12 |
c     | a14  -a31  -a21   a11 |
c-----
c     internal variables
c-----
      real*8 lm, l2m, l3m, mk
      lm = lam + mu
      l2m = lam + mu + mu
      l3m = lam + mu + mu + mu
      mk = mu*k
      saa(1,1) = c + (lm/l2m) * skz

      saa(1,2) = - ( lm *ckz + mu*s)/l2m
      saa(1,3) = - (lm/l2m)* s*z/(2.*mu)
      saa(1,4) = ( l3m *s + lm*ckz)/(l2m*2*mk)

      saa(2,1) = ( -mu * s + lm*ckz)/l2m
      saa(2,2) = c - (lm/l2m)*skz
      saa(2,3) =  ( lm/l2m)* ( s*l3m/lm - ckz)/(2.*mk)
      saa(2,4) = - saa(1,3)

      saa(3,1) = (lm/l2m)*2*mk*skz
      saa(3,2) = (lm/l2m)*2*mk*(s - ckz)
      saa(3,3) =   saa(2,2)
      saa(3,4) = - saa(1,2)

      saa(4,1) = (lm/l2m)*2*mk*(s + ckz)
      saa(4,2) = - saa(3,1)
      saa(4,3) = - saa(2,1)
      saa(4,4) =   saa(1,1)

      return
      end

      subroutine shskl(hl,k,z,c,s,mu)
      implicit none
      real*8 hl(2,2),k,z,c,s,mu
c-----
c     The propagator matrix is actually
c       A = exp(kz) AA
c     because of the factorization within lvar. Both A and AA
c     have a symmetry, e.g.,
c
c-----
c     internal variables
c-----
      real*8  mk
      mk = mu*k

      hl(1,1) = c
      if(k.eq.0.0d+00)then
           hl(1,2) = z/mu
      else
           hl(1,2) = s/mk
      endif
      hl(2,1) = mk * s
      hl(2,2) = c

      return
      end

      
      subroutine sdnka(ca,z,lam,mu,k)
      implicit none
      real*8 ca(5,5),k,z,lam,mu
      common/sovrflw/a0,c2,s2,sc,kz,k2z2
      real*8 a0,c2,s2,sc,kz,k2z2
c-----
c     The propagator matrix is actually
c       Comp(A) = exp(2kz) CA
c     because of the factorization within lvar. 
c     
c     In addition the properties of the compound (6,6) are
c     used to reduce it to a 5x5 to reduce the number of multiplications
c
c     The 6x6 is of the form
c   
c-----
c        A11     A12     A13    -A13     A15     A16
c        A21     A22     A23    -A23     A25     A15
c        A31     A32     A33    1-A33   -A23    -A13
c       -A31    -A32    1-A33    A33     A23     A13
c        A51     A52    -A32     A32     A22     A12
c        A61     A51    -A31     A31     A21     A11
c-----
c       this will be multipled on the left by the G matrix
c
c       [ G11   G12 G13 -G13    G15 G16 ]
c
c-----
c       or on the right by
c
c       [ H11   H21 H31 -H31    H51 H61  ] ^T
c-----
c       the number of multiplications can be reduced from 36 to 25 
c       if we define a new matrices
c       related to the original matrices by
c-----
c         A11     A12     A13         A15     A16
c         A21     A22     A23         A25     A15
c        2 A31   2 A32   2 A33 -1   -2 A23  -2 A13
c         A51     A52    -A32         A22     A12
c         A61     A51    -A31         A21     A11
c-----
c
c       [ G11   G12  G13    G15 G16  ]
c       [ H11   H21 2 H31   H51 H61  ] ^T
c
c-----
c       this means that some of the original definitions of the 
c       Aij elements must be changed for the
c       definition of the modified 5x5 compount A matrix
c
c       old 6x6                 new 5x5
c
c       old 6x6                 new 5x5
c       A11 = 1 - A33 + A22     1 - (1/2)(new A33 + 1) + new A2)
c       A53 = -A32              A43 = - (1/2) new A32
c       A63 = -A31              A53 = - (1/2) new A31
c-----
c       To recover the needed elements, 
c          we note that the old G14 = -old G14 = new G13
c-----
c     internal variables
c-----
      real*8 lm, l2m, l3m, mk
      lm = lam + mu
      l2m = lam + mu + mu
      l3m = lam + mu + mu + mu
      mk = mu*k

c-----
c     the full 6x6 is given here commented out since the
c     reduced form will actually be used
c-----
C         cd(1,1) = c2 + (lm/l2m)*(lm/l2m)*k2z2 - (mu/l2m)*(mu/l2m)*s2
C         cd(1,2) = (lm/l2m)*(1./(2.*mk))*((l3m/lm)*sc - kz)
C         cd(1,3) = 1./(l2m*l2m*2.*mk)*(l3m*mu*s2 - lm*lm*k2z2)
C         cd(1,4) = -cd(1,3)
C         cd(1,5) = - 1./(l2m*2.*mk)*( l3m*SC + lm*kz)
C         cd(1,6) = ( -l3m*l3m*S2 +lm*lm*k2z2)/(l2m*l2m*2*mk*2*mk)
C 
C         cd(2,1) = 2*mk*(lm/l2m)*(sc - kz)
C         cd(2,2) = c2
C         cd(2,3) = ( lm*kz + mu*sc)/l2m
C         cd(2,4) = - cd(2,3)
C         cd(2,5) = - s2
C         cd(2,6) =   cd(1,5)
C         cd(2,6) = - ( l3m*sc +  lm*kz)/(l2m * 2*mk)
C 
C         cd(3,1) = (lm*2*mk/l2m)*(mu*s2/l2m + lm*k2z2/l2m)
C         cd(3,2) = ( mu*sc - lm*kz)/l2m
C         cd(3,3) =  a0 + mu*mu*s2/(l2m*l2m)
C     1     - (lm/l2m)*(lm/l2m)*k2z2 
C         cd(3,4) =   a0 - cd(3,3)
C         cd(3,5) = - cd(2,3)
C         cd(3,6) = - cd(1,3)
C         cd(4,1) = - cd(3,1)
C         cd(4,2) = - cd(3,2)
C         cd(4,3) =   a0 - cd(3,3)
C         cd(4,4) =   cd(3,3)
C         cd(4,5) =   cd(2,3)
C         cd(4,6) =   cd(1,3)
C         cd(5,1) = - lm*(2*mk)*(sc + kz)/l2m
C         cd(5,2) = - s2
C         cd(5,3) = - cd(3,2)
C         cd(5,4) =   cd(3,2)
C         cd(5,5) =   cd(2,2)
C         cd(5,6) =   cd(1,2)
C         cd(6,1) = - (lm*2*mk/l2m)*(lm*2*mk/l2m)*(s2 -k2z2)
C         cd(6,2) =   cd(5,1)
C         cd(6,3) = - cd(3,1)
C         cd(6,4) =   cd(3,1)
C         cd(6,5) =   cd(2,1)
C         cd(6,6) =   cd(1,1)
c-----
c        now apply the relations to make a 5x5 matrix
c-----
         ca(1,1) = c2 + (lm/l2m)*(lm/l2m)*k2z2 - (mu/l2m)*(mu/l2m)*s2
         ca(1,2) = (lm/l2m)*(1./(2.*mk))*((l3m/lm)*sc - kz)
         ca(1,3) = 1./(l2m*l2m*2.*mk)*(l3m*mu*s2 - lm*lm*k2z2)
         ca(1,4) = - 1./(l2m*2.*mk)*( l3m*SC + lm*kz)
         ca(1,5) = ( -l3m*l3m*S2 +lm*lm*k2z2)/(l2m*l2m*2*mk*2*mk)

         ca(2,1) =  2*mk*(lm/l2m)*(sc - kz)
         ca(2,2) =   c2
         ca(2,3) =  (lm*kz + mu*sc)/l2m
         ca(2,4) = - s2
         ca(2,5) =   ca(1,4)

         ca(3,1) = 2.* (lm*2*mk/l2m)*(mu*s2/l2m + lm*k2z2/l2m)
         ca(3,2) = 2.* ( mu*sc - lm*kz)/l2m
         ca(3,3) = a0 + 2.0d+00*(mu*mu*s2/(l2m*l2m)
     1     - (lm/l2m)*(lm/l2m)*k2z2 )
         ca(3,4) = -2.0d+00* ca(2,3)
         ca(3,5) = -2.0d+00* ca(1,3)

         ca(4,1) = - lm*(2*mk)*(sc + kz)/l2m
         ca(4,2) =  -s2
         ca(4,3) =  - ca(3,2)/2.0d+00
         ca(4,4) =  ca(2,2)
         ca(4,5) =  ca(1,2)

         ca(5,1) =  - (lm*2*mk/l2m)*(lm*2*mk/l2m)*(s2 -k2z2)
         ca(5,2) =  ca(4,1)
         ca(5,3) = -ca(3,1)/2.0d+00
         ca(5,4) =  ca(2,1)
         ca(5,5) =  ca(1,1)


      return
      end


      subroutine  sevalg(jbdry,gr,gl,lam,mu,mush,k,im)
c-----
c     The original G matrix of 6 elements is reduced to
c     5 elements by the relation
c       [ G11   G12  G13    G15 G16  ]
c-----
      implicit none
      real*8 gr(2,5),gl(2,2), lam,mu,mush, k
      integer im,jbdry
c-----
c     internal variables
c-----
      integer i
      real*8 dfac

        if(jbdry.lt.0)then
c-----
c               ELASTIC ABOVE - RIGID
c-----
             gr(im,1) = 1.0d+00
             gr(im,2) = 0.0d+00
             gr(im,3) = 0.0d+00
             gr(im,4) = 0.0d+00
             gr(im,5) = 0.0d+00
             gl(im,1) = 1.0d+00
             gl(im,2) = 0.0d+00
        else if(jbdry.eq.0)then
             dfac = (lam + mu)/(2.*k*(lam + 2*mu))
             dfac = dfac*dfac

c-----
C            G(im,1) = -1
C            G(im,2) = -(lam+2*mu)/((lam + mu)*2*k*mu)
C            G(im,3) = - 1/((lam + mu)*2*k)
C            G(im,4) = 1./(2*k*(lam + mu ))
C            G(im,5) = (lam + 2*mu)/(2*mu*k*(lam+mu))
C            G(im,6) = (lam + 3*mu)/(4.*mu*mu*k*k*(lam+mu))
c-----
             GR(im,1) = -1
             GR(im,2) = -(lam+2*mu)/((lam + mu)*2*k*mu)
             GR(im,3) = - 1/((lam + mu)*2*k)
             GR(im,4) = (lam + 2*mu)/(2*mu*k*(lam+mu))
             GR(im,5) = (lam + 3*mu)/(4.*mu*mu*k*k*(lam+mu))
             do i=1,5
                GR(im,i) = GR(im,i) * dfac
             enddo
             GL(im,1) = mush*k
             GL(im,2) = 1
        else if(jbdry.eq.1)then
c-----
c       FREE - check properties of layer above
c-----
             gr(im,1) = 0.0d+00
             gr(im,2) = 0.0d+00
             gr(im,3) = 0.0d+00
             gr(im,4) = 0.0d+00
             gr(im,5) = 1.0d+00
             gl(im,1) = 0.0d+00
             gl(im,2) = 1.0d+00
        endif
      return
      end

      subroutine  sevalh(jbdry,hr,hl,lam,mu,mush,k,im)
c-----
c     The original H matrix of 6 elements is reduced to
c     5 elements by the relation
c       [ H11   H21 2 H31   H51 H61  ] ^T
c-----
      real*8 hr(2,5),hl(2,2), lam,mu,mush, k
      integer im, jbdry
c-----
c       do top surface condition
c
c           = -1  RIGID
c           =  0  ELASTIC
c           = +1  FREE SURFACE
c-----
       if(jbdry.eq.0)then
c-----
c           ELASTIC ABOVE SOLID
c-----
c-----
C              H(im,1) = -k*k*(lam+3.*mu)/(lam + mu)
C              H(im,2) = -2.*mu*k*k*k*(lam +  2.*mu)/(lam + mu)
C              H(im,3) = - 2.*mu*k*k*k*mu/(lam + mu)
C              H(im,4) = 2.*k*k*k*mu*mu/(lam + mu)
C              H(im,5) = 2.*k*k*k*mu*(lam + 2.*mu)/(lam + mu)
C              H(im,6) = 4.*mu*mu*k*k*k*k
c-----
               HR(im,1) = -k*k*(lam+3.*mu)/(lam + mu)
               HR(im,2) = -2.*mu*k*k*k*(lam +  2.*mu)/(lam + mu)
               HR(im,3) = 2.*( - 2.*mu*k*k*k*mu/(lam + mu) )
               HR(im,4) = 2.*k*k*k*mu*(lam + 2.*mu)/(lam + mu)
               HR(im,5) = 4.*mu*mu*k*k*k*k
      
               HL(im,1) = 0.5*k
               HL(im,2) = 0.5*mush*k*k
            else if(jbdry.eq.-1)then
c-----
c           RIGID ABOVE SOLID
c-----
                hr(im,1) = 0.0d+00
                hr(im,2) = 0.0d+00
                hr(im,3) = 0.0d+00
                hr(im,4) = 0.0d+00
                hr(im,5) = 1.0d+00
                hl(im,1) = 0.0d+00
                hl(im,2) = 1.0d+00
            else if(jbdry.eq.1)then
c-----
c           FREE ABOVE SOLID
c-----
                hr(im,1) = 1.0d+00
                hr(im,2) = 0.0d+00
                hr(im,3) = 0.0d+00
                hr(im,4) = 0.0d+00
                hr(im,5) = 0.0d+00
                hl(im,1) = 1.0d+00
                hl(im,2) = 0.0d+00
            endif

      return
      end

        subroutine scopy5(sca,sdar,m,itofrm,dex,exa)
        implicit none
        integer NL
        parameter (NL=200)
        real*8 sdar(NL,5,5)
        real*8 sca(5,5)
        integer itofrm,m
        real*8 dex(NL)
        real*8 exa
        integer i,j
c-----
c       copy from ca to dar
c-----
        if(itofrm.eq.0)then
            do 100 j=1,5
                do 110 i=1,5
                    sdar(m,i,j) = sca(i,j)
  110           continue
                dex(m) = exa
  100       continue
c-----
c       copy from dar to ca
c-----
        else
            do 200 j=1,5
                do 210 i=1,5
                    sca(i,j) = sdar(m,i,j)
  210           continue
                exa = dex(m)
  200       continue
        endif
        return
        end

        subroutine scopy4(saa,shar,m,itofrm,hex,ex)
c-----
c       copy between aa and har arrays for m'th layer
c
c       aa  R   - 4x4 Haskell matrix array
c       har R   - NLx4x4 storage array
c       m   I   - layer index
c       itofrm  I   - 0 copy aa to har, ex to hex
c                 !=0 copy from har to aa, hex to ex
c       hex R   - NL array - storage for exponent
c       ex  R   - exponent value
c-----
        integer NL
        parameter (NL=200)
        real*8 shar(NL,4,4)
        real*8 saa(4,4)
        integer itofrm,m
        real*8 hex(NL)
        real*8 ex
        integer i,j
c-----
c       copy from aa to har
c-----
        if(itofrm.eq.0)then
            do 100 j=1,4
                do 110 i=1,4
                    shar(m,i,j) = saa(i,j)
  110           continue
                hex(m) = ex
  100       continue
c-----
c       copy from har to aa
c-----
        else
            do 200 j=1,4
                do 210 i=1,4
                    saa(i,j) = shar(m,i,j)
  210           continue
                ex = hex(m)
  200       continue
        endif
        return
        end

        subroutine scopy2(shl,shal,m,itofrm,lex,exb)
c-----
c       itofrm     0 place simple 2x2 array into storage of all layers
c                  1 read from the storage of all layers
c-----
        implicit none
        integer NL
        parameter (NL=200)
        real*8 shal(NL,2,2)
        real*8 shl(2,2)
        integer itofrm,m
        real*8 lex(NL)
        real*8 exb

        integer j,i
c-----
c       copy from hl to hal
c-----
        if(itofrm.eq.0)then
            do 100 j=1,2
                do 110 i=1,2
                    shal(m,i,j) = shl(i,j)
  110           continue
                lex(m) = exb
  100       continue
c-----
c       copy from hal to hl
c-----
        else
            do 200 j=1,2
                do 210 i=1,2
                    shl(i,j) = shal(m,i,j)
  210           continue
                exb = lex(m)
  200       continue
        endif
        return
        end
        subroutine scmult(se,sca,exa,exe)
        implicit none
c-----
c       FORM EC where se(1x5) sc(5x5)
        real*8 sca(5,5)
        real*8 exa,exe,eval
        real*8 xnorm
        real*8 se(5)
        real*8 sc, see(5)

        integer i,j,IUP

            IUP = 5
        do 1350 i=1,IUP
            sc = 0.0d+00
            do 1349 j=1,IUP
                sc=sc+ se(j) * sca(j,i)
 1349       continue
            see(i)=sc
 1350   continue
        exe = exe + exa
        call snormc(see,eval,xnorm)
        do 1351 i=1,IUP
            se(i) = see(i)*xnorm
 1351   continue
        exe = exe + eval
        return
        end

        subroutine srcmult(sy,sc,exa,exe)
c-----
c       FORM YC where y(5x5) c(5x5) RETURN Y
c-----
        real*8 sc(5,5)
        real*8 exa,exe,eval
        real*8 sxnorm
        real*8 sy(5,5)
        real*8 stmp, see(5,5)

        integer i,j,IUP,k
            IUP = 5
        do 1350 i=1,IUP
            do 1351 j=1,IUP
                stmp = 0.0d+00
                do 1349 k=1,IUP
                    stmp=stmp+ sy(i,k) * sc(k,j)
 1349           continue
                see(i,j)=stmp
 1351       continue
 1350   continue
        exe = exe + exa
        call rnormc(see,eval,sxnorm)
        do 1353 j=1,IUP
            do 1352 i=1,IUP
                sy(i,j) = see(i,j)*sxnorm
 1352       continue
 1353   continue
        exe = exe + eval
        return
        end

        subroutine dmult(sda,saa)
        implicit none
c-----
c       propagator up
c       FORM D = DA
c-----
        real*8 saa(4,4)
        real*8 sumd,sea(4,4),sda(4,4)

        integer i,j,jj
        do 1360 i=1,4
            do 1361 j=1,4
                sumd = 0.0d+00
                do 1362 jj=1,4
                    sumd=sumd+sda(i,jj) * saa(jj,j)
 1362           continue
                sea(i,j)=sumd
 1361       continue
 1360   continue
        do 1363 j=1,4
            do 1364 i=1,4
                sda(i,j)=sea(i,j)
 1364       continue
 1363   continue
        return
        end

        subroutine strans4(sa)
        implicit none
        real*8 sa(4,4)
c-----
c       from a transpose 
c-----
        real*8 sdum
        integer i,j
        do 100 i=1,4
            do 101 j=i,4
                sdum = sa(i,j)
                sa(i,j) = sa(j,i)
                sa(j,i) = sdum
  101       continue
  100   continue
        return
        end

        subroutine lmult(d11,d12,d21,d22,hl,exel,exb,icomp)
        implicit none
c-----
c       multiply SH matrix by a row vector on left
c-----
        real*8 d11,d12,d21,d22,hl(2,2),e1,e2
        real*8 exel, exb
        logical icomp
c-----
c       elastic layer
c-----
            e1=d11
            e2=d12
c-----
c           a11 = cosql
c           a12 = yl
c           a21 = zl
c           a22 = cosql
c-----
            d11=e1*hl(1,1) + e2*hl(2,1)
            d12=e1*hl(1,2) + e2*hl(2,2)
            exel = exel + exb
            if(icomp)then
                e1=d21
                e2=d22
                d21=e1*hl(1,1) + e2*hl(2,1)
                d22=e1*hl(1,2) + e2*hl(2,2)
            endif
        return
        end

        subroutine snormc(se,ex,xnorm)
        implicit none
        real*8 ex
        real*8 testt,xnorm
        real*8 se(*)
        integer i, iup
        testt = 0.0D+00
            IUP = 5
        do i = 1,IUP
            if(dabs(se(i)).gt.testt) testt =dabs(se(i))
        enddo
        if(testt.lt.1.0e-30)testt=1.0
        xnorm = 1./testt
        ex = - dlog(xnorm)
        return
        end

        subroutine rnormc(e,ex,xnorm)
        implicit none
        real*8 ex, xnorm
        real *8 testt
        real*8 e(5,5)
        integer i,j,iup
        testt = 0.0D+00
            IUP = 5
        do  j=1,IUP
            do  i = 1,IUP
            if(dabs(e(i,j)).gt.testt) testt =dabs(e(i,j))
            enddo
        enddo
        if(testt.lt.1.0e-30)testt=1.0
        xnorm = 1./testt
        ex = -dlog(xnorm)
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
        if(iw.gt.0 .and. .not. lfluid)then
            do 102 i=isoldd,isoldu
                if(b(i).eq.0.0)then
                call werror('MODEL HAS INTERIOR  FLUID LAYERS')
                endif
  102       continue
        endif
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
        integer NL
        parameter (NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        integer mmax
        common/depref/refdep
        real refdep
        common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
        real bsh, rhosh, qbsh

        double precision z0,z1,r0,r1,dr,ar,tmp

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
        integer NL
        parameter (NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        integer mmax
        common/depref/refdep
        real refdep
        common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
        real bsh, rhosh, qbsh

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
        common/verby/verbose
        logical verbose
c-----
c       initialize bound search
c-----
        vamin = 1.0e+38
        vbmin = 1.0e+38
        vmin  = 1.0e+38
        vamax = 0.0
        vbmax = 0.0
        if(verbose)write(LOT,2) 
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
            if(verbose)then
            if(i.lt.mmax)then
               write(LOT,3)d(i),a(i),b(i),rho(i),qa(i),qb(i),
     1             bsh(i),rhosh(i)
            endif
            endif
   20   continue 
    5   format(' ',10x,3f10.3,2f10.6,2f10.3/' ') 
        if(verbose)then
           write(LOT,5)a(mmax),b(mmax),rho(mmax),qa(mmax),qb(mmax) ,
     1          bsh(mmax),rhosh(mmax)
        endif
c-----
c     obtain extreme velocity limits
c-----
      return 
      end 

        subroutine scoef(cd,da,fr,exe,exl,exwu,wvno,
     1      fl,d11,d12,exel,exll,llmaxs,llmaxr)
c-----
c       llmaxs  I*4 - source layer index in original model
c       llmaxr  I*4 - receiver layer index in original model
c-----
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
        real bsh, rhosh, qbsh
        common/modlly/mmax
        common/damp/alpha,ieqex
        real*8  da(4,4),ca(5,5)
        real*8   cd(5),e(5),fr
        real*8 d11,d12,e1,e2,e21, e22, fl
        real*8 exe,exl,exel,exll,ex,exa,exb,exwu
        real*8 dzero
        real*8 wvno
        real*8 zdum
        real*8 aa(4,4)
        real*8 cy(5,5)
        real*8 zero, zone
        real*8 y11, y12, y21, y22, sd11, sd21
        real*8 hl(2,2)
c-----
c       matrix components in layers and boundaries saved
c-----
        real*8 har(NL,4,4), dar(NL,5,5), hsr(2,5), gbr(2,5), 
     1      hal(NL,2,2), hsl(2,2), gbl(2,2)
        real*8 hex(NL), lex(NL), dex(NL), hexw(NL)
        common/hamat/har
        common/damat/dar
        common/hsrfr/hsr
        common/gbrfr/gbr
        common/hlmat/hal
        common/hsrfl/hsl
        common/gbrfl/gbl
        common/hexex/hex
        common/hexexw/hexw
        common/dexex/dex
        common/lexex/lex 
        common/water/iwater(NL),iwats(2),iwatb(2)
        common/updnsm/equalu(NL), equald(NL)
        logical equalu, equald
        logical retrieve
c-----
c       check for decomposition at wavefield at receiver
c-----
        common/ctlrud/rpup,rpdn,rsup,rsdn,dorud
        logical rpup, rpdn, rsup, rsdn, dorud
        common/ctlsud/spup,spdn,ssup,ssdn,dosud
        logical spup, spdn, ssup, ssdn, dosud

        REAL*8 DDEX, DXNORM

c-----
c       this routine computes the layer response. 
c       To simplify the mathematics
c       of the case of receiver above or beneath the source, the
c       layer is internally flipped.
c
c       llmaxs  I*4 - source layer index in original model
c       llmaxr  I*4 - receiver layer index in original model
c       mm  I*4 - pointer to layer in original model
c
c       lmaxs   I*4 - source layer index in current model
c       lmaxr   I*4 - receiver layer index in current model
c       mm  I*4 - pointer to layer in current model
c
c       in  I*4 - 1 use current model (source beneath receiver)
c               - 2 use inverted model (receiver beneath source)
c-----
c       initialize matrices
c-----
        zero = 0.0d+00
        zone  = 1.0d+00
        dzero = 0.0d+00
        exe=0.0
        exl=0.0
        exwu = 0.0
        do 2 j = 1,4
            do 3 i = 1,4
                da(i,j)=zero
    3       continue
            da(j,j) = zone
    2   continue
        do 12 j=1,5
            do 13 i=1,5
                cy(i,j) = zero
   13       continue
            cy(j,j) = zone
   12   continue
        y11 = zone
        y12 = zero
        y21 = zero
        y22 = zone
        exel = 0.0
        exll = 0.0
c-----
c     set up halfspace conditions
c-----
        if(llmaxs .ge. llmaxr)then
            im = 1
        else
            im = 2
        endif
        do 100 i=1,5
            e(i) = gbr(im,i)
  100   continue
c-DEBUG
        call snormc(e,DDex,Dxnorm)
        do i=1,5
          e(i) = e(i)*Dxnorm
        enddo
c-DEND DEBUG
        e1 = gbl(im,1)
        e2 = gbl(im,2)
        
        do 11 i=1,5
            cd(i)=e(i)
   11   continue
        d11=e1
        d12=e2
c-----
c       set up limits on the layer stacking
c-----
        if(llmaxs .ge. llmaxr)then
            lmaxs = llmaxs
            lmaxr = llmaxr
        else
            lmaxs = mmax - llmaxs + 2
            lmaxr = mmax - llmaxr + 2
        endif
c-----
c       matrix multiplication from bottom layer upward
c-----
        do 1340 mm = mmax,1,-1
            if(llmaxs .ge. llmaxr)then
                m = mm
                if(equalu(m))then
                    retrieve = .false.
                else
                    retrieve = .true.
                endif
            else
                m = mmax + 1 - mm
                if(equald(m))then
                    retrieve = .false.
                else
                    retrieve = .true.
                endif
            endif
            iwat = iwater(m)
            if(retrieve)then
                call scopy5(ca,dar,m,1,dex,exa)
                call scopy2(hl,hal,m,1,lex,exb)
                call scopy4(aa,har,m,1,hex,ex)
            endif
            call scmult(e,ca,exa,exe)
            call lmult(e1,e2,e21,e22,hl,exel,exb,.false.)
            if(mm.lt.lmaxr)then
                call srcmult(cy,ca,exa,exl)
                call lmult(y11,y12,y21,y22,hl,
     1              exll,exb,.true.)
            else if(mm.ge.lmaxr .and. mm.lt.lmaxs) then
                call dmult(da,aa)
                    exl = exl + ex
c-----
c       save values at top of source layer
c-----
            else if(mm.eq.lmaxs) then
                    do 1352 i=1,5
                    cd(i)=e(i)
 1352           continue
                exl=exe
                exll = exel
                d11=e1
                d12=e2
            endif
            if(mm.eq.1)then
                do 200 i=1,5
                    ca(i,1) = hsr(im,i)
  200           continue
                sd11 = hsl(im,1)
                sd21 = hsl(im,2)
                zdum = e1
                e1 = zdum*sd11 + e2*sd21
c               e2 = zdum*sd11 - e2*sd21
                zdum = y11
                y11 = zdum*sd11 + y12*sd21
                zdum = y21
                y21 = zdum*sd11 + y22*sd21
                zdum = 0.0d+00
                do 1402 i=1,5
                    zdum = zdum + e(i)*ca(i,1)
 1402           continue
                e(1) = zdum
                call srcmult(cy,ca,dzero,exl)
            endif
 1340 continue
c-----
c       get final matrices
c-----
c-SH
        fl=e1
c-P-SV
c       form x(l,m)y(ij|12)
c-----take care of x(i,j) y(1j|12) and replace the da
            aa(1,1) =   zero
            aa(2,1) =   cy(1,1)
            aa(3,1) =   cy(2,1)
            aa(4,1) =   cy(3,1) /(2.0d+00 )
c change sign 0430 1200
            aa(1,2) = - aa(2,1)
            aa(2,2) =   zero
            aa(3,2) = - aa(4,1)
            aa(4,2) =   cy(4,1)
            aa(1,3) = - aa(3,1)
            aa(2,3) = - aa(3,2)
            aa(3,3) =   zero
            aa(4,3) =   cy(5,1)
            aa(1,4) = - aa(4,1)
            aa(2,4) = - aa(4,2)
            aa(3,4) = - aa(4,3)
            aa(4,4) =   zero
            call dmult(da,aa)
            fr=e(1)
            d11 = y11*d11
            d12 = y11*d12
        return
        end

        subroutine evlmat(wvno,jbdrys,jbdryh)
        implicit none
c-----
c       routine carugments
c-----
        real*8 wvno
        integer jbdrys,jbdryh
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
        real bsh, rhosh, qbsh
        common/modlly/mmax
        integer mmax
        integer NSOURCE, NRECEIVER, NSR
        parameter(NSOURCE=100,NRECEIVER=100,NSR=100)
        common/source/depths(NSOURCE),lmaxs(NSOURCE),mdpths
        real depths
        integer lmaxs, mdpths
        common/receiv/depthr(NRECEIVER),lmaxr(NRECEIVER),mdpthr
        real depthr
        integer lmaxr, mdpthr
        real*8  ca(5,5), hl(2,2)
        real*8 ex,exa
        real*8 aa(4,4)
c-----
c       matrix components in layers and boundaries saved
c-----
        common/hamat/har
        real*8 har(NL,4,4)
        common/damat/dar
        real*8 dar(NL,5,5)
        common/hsrfr/hsr
        real*8 hsr(2,5)
        common/gbrfr/gbr
        real*8 gbr(2,5)
        common/hlmat/hal
        real*8 hal(NL,2,2)
        common/hsrfl/hsl
        real*8 hsl(2,2)
        common/gbrfl/gbl
        real*8 gbl(2,2)
        common/hexex/hex
        real*8 hex(NL)
        common/hexexw/hexw
        real*8 hexw(NL)
        common/dexex/dex
        real*8 dex(NL)
        common/lexex/lex 
        real*8 lex(NL)
        common/water/iwater(NL),iwats(2),iwatb(2)
        integer iwater, iwats, iwatb
        common/updnsm/equalu(NL), equald(NL)
        logical equalu, equald

        logical compute

        real*8 mu, lam
        real*8 mush
        real*8 c,s,skz,ckz

        integer m
c-----
c       evaluate the G matrix 
c           gbr(1,x) is for normal stack
c           gbr(2,x) is for inverted stack
c-----
        mu = rho(mmax)*b(mmax)*b(mmax)
        lam = rho(mmax)*a(mmax)*a(mmax) -2.*mu
        mush  = rhosh(mmax)*bsh(mmax)*bsh(mmax)
        call sevalg(jbdryh,gbr,gbl,lam,mu,mush,wvno,1)

        mu  = rho(1)*b(1)*b(1)
        lam = rho(1)*a(1)*a(1) -2.*mu
        mush  = rhosh(1)*bsh(1)*bsh(1)
        call sevalg(jbdrys,gbr,gbl,lam,mu,mush,wvno,2)
c-----
c       evaluate the H matrix
c           hsr(1,x) is for normal stack
c           hsr(2,x) is for inverted stack
c-----
        mu  = rho(1)*b(1)*b(1)
        lam = rho(1)*a(1)*a(1) -2.*mu
        mush  = rhosh(1)*bsh(1)*bsh(1)
        call sevalh(jbdrys,hsr,hsl,lam,mu,mush,wvno,1)

        mu  = rho(mmax)*b(mmax)*b(mmax)
        lam = rho(mmax)*a(mmax)*a(mmax) -2.*mu
        mush  = rhosh(mmax)*bsh(mmax)*bsh(mmax)
        call sevalh(jbdryh,hsr,hsl,lam,mu,mush,wvno,2)
c-----

c-----
c       matrix multiplication from bottom layer upward
c-----
        do 1340 m = 1,mmax,1
c-----
c       first check to see if computations already done
c-----
            if(equald(m))then
                compute = .false.
            else
                compute = .true.
            endif
            if(compute)then
                mu = rho(m)*b(m)*b(m)
                lam = rho(m)*a(m)*a(m) -2.*mu
                call svar (wvno,dble(d(m)),ex,exa,c,s,skz,ckz)
                call shska(aa,wvno,dble(d(m)),c,s,skz,ckz,lam,mu)
                call sdnka(ca,dble(d(m)),lam,mu,wvno)
                mush = rhosh(m)*bsh(m)*bsh(m)
                call shskl(hl,wvno,dble(d(m)),c,s,mush)

            endif
            call scopy5(ca,dar,m,0,dex,exa)
            call scopy4(aa,har,m,0,hex,ex)
            call scopy2(hl,hal,m,0,lex,ex)
 1340   continue
        return
        end


        subroutine srshof(gg,wvno, lmaxs, lmaxr) 
        implicit none
        real gg(21)
        real*8 wvno
        integer lmaxs, lmaxr

        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        real d,a,b,rho,qa,qb,etap,etas,frefp,frefs
        common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
        real bsh, rhosh, qbsh
        common/modlly/mmax
        integer mmax
        common/jout/jsrc(21) , jbdrys, jbdryh
        integer jsrc, jbdrys, jbdryh
c-----
c       gus - surface displacements or potentials or top of layer
c-----
        real*8 gus(21)
        real*8 cd(5),da(4,4),fr,y(4,4)
        real*8 d11,d12,fl 
        real*8 s21,s32,s14,s34,s32e,s34e 
        real*8 s24,s33
        real*8 wv4pi
        real*8 zero
        real*8 fact,exe,exl,exel,exll,elj
        real *8 exwu
C        real*8 sdd(4), sds(4), sss(4), sep(4), svf(4), shf(4)
        real*8 shds(2), shss(2), shhf(2)
        real*8 fourpi

        real*8 lam, mu, l2m
        integer i, j, k, jj


c-----
c       Initialization
c-----
        fourpi=12.5663706d+00
        zero  = 0.0d+00
c-----
c       do not evaluate for wvno = 0.0
c-----
        if(wvno.eq.0.0d+00) then
            do 102 i=1,21
                gg(i) = 0.0
  102       continue
        else
c-----
c       process for this wavenumber and frequency
c-----
            do 101 i = 1,21
                gus(i) = zero
  101       continue
            wv4pi = 2.0d+00 * wvno / fourpi
            mu = rho(lmaxs)*b(lmaxs)*b(lmaxs)
            l2m = rho(lmaxs)*a(lmaxs)*a(lmaxs)
            lam = l2m - mu - mu
c-----
c           this follows Wang and Herrmann (1980), 
c               A numerical study of P-, SV-, and SH-wave 
c               generation in a plane layerd medium,
c               Bull. Seism. Soc. Am. 70, 1015-1036.
c           For p-SV we use equation (8). Note below in "do 60" 
c               that the Z component (odd index) is reversed 
c               from (8) so that +z is up
c
c       using the trick of reducing the 6x6 compound matrices to 5x5, 
c           the
c       following correspondence between what is given here 
c       and that in Wang
c       and Herrmann is in effect
c-----
c       Compound Matrix (6x6)   Here (5x5)
c
c       X|12/12 ==      cd(1)
c       X|12/13 ==      cd(2)
c       X|12/14 ==      cd(3)
c       X|12/23 ==      -cd(3)
c       X|12/24 ==      cd(4)
c       X|12/34 ==      cd(5)
c       X|12/21 = -X|12/12
c       X|12/ij = -X|12/ji
c
c       where
c
c        |12
c       X|  = X|12/ij
c        |ij
c-----
            call scoef(cd,da,fr,exe,exl,exwu,wvno,
     1          fl,d11,d12,exel,exll,lmaxs,lmaxr)
C        WRITE(6,*)'G( WVNO:',wvno
C        WRITE(6,*)'G( CD:',cd
C        WRITE(6,*)'G( FR:',fr,exe,exl
C        WRITE(6,*)'G( DA:',da
C
C        WRITE(6,*)'wvno,fl,d11,d12,exel,exll:',wvno,fl,d11,d12,exel,exll
c-----
c       Form X|12/ij x Zj2 in Equation 8 of Wang and Herrmann
c----
c KLUDGE to CHANGE ORDER AND ALSO GET UR CORRECT FOR WATER LAYER
            do 50 k=1,2
c-----
c       k =1 UR     k=2 UZ
c       jj=1 = UZ, jj=2 = UR
c-----
                if(k.eq.1)then
                    jj = 2
                else if(k.eq.2)then
                    jj = 1
                endif
                    j = k
                y(1,jj)= cd(1)*da(2,j) + cd(2)*da(3,j) 
     1              + cd(3)*da(4,j)
                y(2,jj)= -cd(1)*da(1,j) - cd(3)*da(3,j) 
     1              + cd(4)*da(4,j)
                y(3,jj)=-cd(2)*da(1,j) + cd(3)*da(2,j) 
     1              + cd(5)*da(4,j)
                y(4,jj)= -cd(3)*da(1,j) - cd(4)*da(2,j) 
     1              - cd(5)*da(3,j)
   50       continue
c-----
c           evaluate different Green's functions
c           apply source terms
c----- 
c-----
c           START OF P-SV
c-----
c           First compute the DELTA displacement-stress source terms
c           for inverted model, the UZ, TR elements change, These will
c           be only those required for dipoles and forces
c           
c           Stress-displacement discontinuities for Green's functions
c           Grn dUr dUz dTz dTr
c           DD      s32     s34
c           DS  s21
c           SS              s14
c           EX      s32e    s34e
c           VF          s33
c           HF              s24
c-----
                s32  = 4./(fourpi * l2m)
                s34  = - wv4pi*(3.*lam+2.*mu)/l2m

                s21  = 2.0d+00/(fourpi * mu)

                s14  = -wv4pi

                s32e = 2.0d+00/(fourpi * l2m) 
                s34e = 2.0d+00*wv4pi*(mu/l2m)

                s33  = -2.0d+00/fourpi

                s24  = -2.0d+00/fourpi
c-----
c           receiver beneath the source
c-----
            if(lmaxr .gt. lmaxs)then
                s14  = - s14
                s24  = - s24
                s32  = - s32
                s32e = - s32e
                s34  = - s34
                s34e = - s34e
            endif
c-----
c           For complete wavefield computation do not
c           waste cycles computing W matrix elements
c-----
                do 61 j=1,2
c       DD
                    gus(j   )=s32 *y(2,j)+ s34*y(4,j)
c       DS
                    gus(j+ 2)=s21 *y(1,j)             
c       SS
                    gus(j+ 4)=             s14*y(4,j)
c       EX
                    gus(j+ 6)=s32e*y(2,j)+s34e*y(4,j)
c       VF
                    gus(j+ 8)=s33 *y(3,j)
c       HF
                    gus(j+10)=s24 *y(4,j)
   61           continue
c-----
c       invert the vertical
c-----
            do 62 j=1,12,1
                gus(j) = -gus(j)
   62       continue
c-----
c           if receiver beneath the source unflip radial
c-----
            if(lmaxr .gt. lmaxs)then
                do 63 j=2,12,2
                    gus(j) = - gus(j)
   63           continue
            endif
c-----
c           END OF P-SV
c-----
c           START OF SH
c-----
                shds(1) = - 2.0d+00/(rhosh(lmaxs)*
     1              12.5663706d+00*
     1              bsh(lmaxs)*bsh(lmaxs))
                shds(2) = zero
                shss(1) = zero
                shss(2) =  2.0d+00*wvno/12.5663706d+00
                shhf(1) = zero
                shhf(2) =  2.0d+00/12.5663706d+00
                    gus(13) = - ( d11*shds(1)           )
                    gus(14) = - (            d12*shss(2))
                    gus(15) = - (            d12*shhf(2))
                if(lmaxr .gt. lmaxs)then
                    gus(13) = - gus(13)
                endif
c-----
c           END OF SH
c-----
c-----
c           do final scaling for exponential
c-----
c-----
c           SV
c-----
c           fix for radial derived from vertical for fluid
c-----
            do 71 k=1,2
                elj = -exe + exl 
                fact = 0.0D+00
                if(elj.gt.-55.) fact=dexp(elj)
                do 72 i=0,10,2
                    j = i + k
                    gg(j) = ( gus(j) * fact/fr)
c-----
c           flip UZ to make vertical positive up
c-----
                    if(k.eq.1)then
                        gg(j) = -gg(j)
                    endif
   72           continue
c----
c           do pressure field
c----
                if(k.eq.1)then
                    gg(16) = - (gus(16)*fact/fr)
                    gg(17) = - (gus(17)*fact/fr)
                    gg(18) = - (gus(18)*fact/fr)
                    gg(19) = - (gus(19)*fact/fr)
                    gg(20) = - (gus(20)*fact/fr)
                    gg(21) = - (gus(21)*fact/fr)
                endif
   71       continue
c-----
c           SH
c-----
            elj=-exel+exll
                if(elj.gt.-55.) then
                    fact = dexp(elj)
                    gg(13)=(gus(13)*fact)/(fl)
                    gg(14)=(gus(14)*fact)/(fl)
                    gg(15)=(gus(15)*fact)/(fl)
                else
                    gg(13) = 0.0
                    gg(14) = 0.0
                    gg(15) = 0.0
                endif
        endif
        return
        end

        subroutine gasym(g,wvno,depth)
        common/kint4/aa,bb,cc
            real aa(21),bb(21),cc(21)
        real g(21)
        ex = exp(-wvno*depth)
        do  j=1,21
                g(j)=g(j) - ex*(aa(j)+wvno*(bb(j)+
     1              wvno*(cc(j))))
        enddo
        return
        end

        subroutine intini(smm,r)
        common/asym/j0k0,j0k1,j0k2,j0k3,
     1      j1k0,j1k1,j1k2,j1k3,
     2      j2k0,j2k1,j2k2,j2k3,
     3      j1k0r,j1k1r,j1k2r,j1k3r,
     4      j1km1r,j2k0r,j2k1r,j2k2r
        real*4 j0k0,j0k1,j0k2,j0k3,
     1      j1k0,j1k1,j1k2,j1k3,
     2      j2k0,j2k1,j2k2,j2k3,
     3      j1k0r,j1k1r,j1k2r,j1k3r,
     4      j1km1r,j2k0r,j2k1r,j2k2r
        common/kint4/aa,bb,cc
            real aa(21),bb(21),cc(21)
        real smm(21),sumd
        common/kint1/gasymp
            logical gasymp
        common/rlimit/rlim
        real*4 rlim
            if(gasymp)then
c-----
c       set up sum arrays, but put in asymptotic value now
c       of setting to zero and then resetting
c-----
                smm(1)=         aa(1)*j0k1 
     1                  + bb(1)*j0k2 
     2                  + cc(1)*j0k3
                smm(2)=         aa(2)*j1k1 
     1                  + bb(2)*j1k2 
     2                  + cc(2)*j1k3
                smm(3)=         aa(3)*j1k1 
     1                  + bb(3)*j1k2 
     2                  + cc(3)*j1k3
                    sumd  = (aa(4)+aa(13))*j1k0r 
     1                  + (bb(4)+bb(13))*j1k1r 
     2                  + (cc(4)+cc(13))*j1k2r
                    sumd = - sumd
                smm(4)= sumd  + aa(4) *j0k1 
     1                  + bb(4) *j0k2 
     2                  + cc(4) *j0k3
                smm(5)= sumd  + aa(13)*j0k1 
     1                  + bb(13)*j0k2 
     2                  + cc(13)*j0k3
                smm(6)=     aa(5)*j2k1 
     1                  + bb(5)*j2k2 
     2                  + cc(5)*j2k3
                     sumd= (aa(6)+aa(14))*j2k0r
     1                  + (bb(6)+bb(14))*j2k1r
     2                  + (cc(6)+cc(14))*j2k2r
                    sumd = -2.*sumd
                smm(7)= sumd  + aa(6) *j1k1 
     1                  + bb(6) *j1k2 
     2                  + cc(6) *j1k3
                smm(8)= sumd  + aa(14)*j1k1
     1                  + bb(14)*j1k2
     2                  + cc(14)*j1k3
                smm(9)=         aa(7) *j0k1        
     1                  + bb(7) *j0k2 
     2                  + cc(7) *j0k3
                smm(10)=        aa(8) *j1k1 
     1                  + bb(8) *j1k2 
     2                  + cc(8) *j1k3
                smm(11)=        aa(9) *j0k1 
     1                  + bb(9) *j0k2 
     2                  + cc(9) *j0k3
                smm(12)=        aa(10)*j1k1 
     1                  + bb(10)*j1k2 
     2                  + cc(10)*j1k3
                smm(13)=        aa(11)*j1k1 
     1                  + bb(11)*j1k2 
     2                  + cc(11)*j1k3
c----- 
c   attempt to look at functional form A/k + b + Ck
C                    sumd  = (aa(12)+aa(15))*j1km1r 
C     1                  + (bb(12)+bb(15))*j1k0r 
C     2                  + (cc(12)+cc(15))*j1k1r
C                    sumd = - sumd
C                smm(14)= sumd + aa(12)*j0k0 
C     1                  + bb(12)*j0k1 
C     2                  + cc(12)*j0k2
C                smm(15)= sumd + aa(15)*j0k0 
C     1                  + bb(15)*j0k1 
C     2                  + cc(15)*j0k2
c-----
                    sumd  = (aa(12)+aa(15))*j1k0r 
     1                  + (bb(12)+bb(15))*j1k1r 
     2                  + (cc(12)+cc(15))*j1k2r
                    sumd = - sumd
                smm(14)= sumd + aa(12)*j0k1 
     1                  + bb(12)*j0k2 
     2                  + cc(12)*j0k3
                smm(15)= sumd + aa(15)*j0k1 
     1                  + bb(15)*j0k2 
     2                  + cc(15)*j0k3
                smm(16)=        aa(16)*j0k1 
     1                  + bb(16)*j0k2 
     2                  + cc(16)*j0k3
                smm(17)=        aa(17)*j0k1 
     1                  + bb(17)*j0k2 
     2                  + cc(17)*j0k3
                smm(18)=        aa(18)*j1k1 
     1                  + bb(18)*j1k2 
     2                  + cc(18)*j1k3
                smm(19)=        aa(19)*j2k1 
     1                  + bb(19)*j2k2 
     2                  + cc(19)*j2k3
                smm(20)=        aa(20)*j0k1 
     1                  + bb(20)*j0k2 
     2                  + cc(20)*j0k3
                smm(21)=        aa(21)*j1k1 
     1                  + bb(21)*j1k2 
     2                  + cc(21)*j1k3
            else
                do 100 i=1,21
                    smm(i)=0.0
  100           continue
            endif
 1010   continue
 1000   continue
        return
        end

        subroutine setup(rr,zz,rlim) 
c---------------------------------------------------------- 
c 
c       jnkm =  integral exp(-kh) krsup m j sub n (kr) dk 
c
c       This is used in the fit of low frequency information
c
c       integral f Jn(kr) dk = 
c           integral [ f - (a+bk+ck^2)e^{-kh} ] Jn(kr) dk
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
        real*4 rr,zz
        real*4 rlim
        common/asym/j0k0,j0k1,j0k2,j0k3,
     1      j1k0,j1k1,j1k2,j1k3,
     2      j2k0,j2k1,j2k2,j2k3,
     3      j1k0r,j1k1r,j1k2r,j1k3r,
     4      j1km1r,j2k0r,j2k1r,j2k2r
        real*4 j0k0,j0k1,j0k2,j0k3,
     1      j1k0,j1k1,j1k2,j1k3,
     2      j2k0,j2k1,j2k2,j2k3,
     3      j1k0r,j1k1r,j1k2r,j1k3r,
     4      j1km1r,j2k0r,j2k1r,j2k2r
                    r = dble(rr)
                    z = dble(zz)
                    dist=sqrt(r*r + z*z) 
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
                    j0k0 = ( 1.0/dist   )
                    j0k1 = ( z/dist3   )
                    j0k2 = ( (2.0*z2 - r2)/dist5   )
                    j0k3 = ( (6.0*z3 - 9.0*z*r2)/dist7   )
                if(rr.le.rlim)then
                    j1k0 = 0.0
                    j1k1 = 0.0
                    j1k2 = 0.0
                    j1k3 = 0.0
                    if(zz .le.rlim)then
                        j1km1r=(1.0d+00/zdist)
                        j1k0r = (0.5/(dist*dist))
                    else
                        j1km1r = (0.5/z)
                        j1k0r = (0.5/z2)
                    endif
                    j1k1r = (1.0/dist3)
                        j1k2r = ( 3.0*z/dist5   )
                        j1k3r = ( 12.0*z2 /dist7   )
                    
                else
c                       j1k0 = (( 1.0 -z/dist   )/r       )
                        j1k0 = ((r/(dist*zdist) )      )
                        j1k1 = ( r/dist3                  )
                        j1k2 = ( 3.0*z*r/dist5            )
                        j1k3 = ( 3.0*r*(4.0*z2 - r2)/dist7)
                    j1km1r=(1.0d+00/(zdist)          )
c                   j1k0r = j1k0/rr
                        j1k0r = (1.0d+00/(dist*zdist))
                    j1k1r = j1k1/rr
                    j1k2r = j1k2/rr
                    j1k3r = j1k3/rr
                endif
                if(rr.le.rlim)then
                        j2k0 = 0.0
                        j2k1 = 0.0
                        j2k2 = 0.0
                        j2k0r = 0.0
                        j2k1r = 0.0
                        j2k2r = 0.0
                else
                        j2k0=((1.0 -zor)*(1.0-zor)*(dist/r2))
c                       j2k0=(((r/zdist)**2)/dist)
                        j2k1=((1.0-zor)*(1.0-zor)*(2.0+zor)/r2)
c                       j2k1=j2k0*( (zdist+dist)/dist3)
                        j2k2 = ( 3.0*r2/dist5   )
                        j2k3 = ( 15.0*z*r2/dist7   )
                        j2k0r = j2k0/r
                        j2k1r = j2k1/r
                        j2k2r = j2k2/r
                endif
        return 
        end 

        subroutine solu(y1,y2,y3,x1,x2,x3,h,j,a,b,c)
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
        real y1,y2,y3
        real a,b,c
        real x1,x2,x3
        integer imap(21)
        complex double precision da, db, dc
         double precision eu1, eu2, eu3
         double precision a11, a12, a13, a22, a23, a33, ddet
COK        data imap/2,2,2,2,2,2,2,2,3,3,3,2,2,2,2,2,2,2,2,2,2/
C        data imap/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
        data imap/9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9/
        a=0.0
        b=0.0
        c=0.0
        wfac = x1*h
c       if(wfac.gt.10.0)return
        ii = imap(j)
        if(ii.eq.1)then
c-----
c       a exp(-kh)
c-----
            b=0.0
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
            a = 0.0
            b = x2*(y1*exp(u1))/x1 - x1*(y2*exp(u2))/x2
            b = b/det
            c= y2*exp(u2)/x2 - y1*exp(u1)/x1
            c = c/det
        else if(ii.eq.4)then
c-----
c       [ c k*k ] exp(-kh)
c-----
            a = 0.0
            b = 0.0
            c = y1 * exp(x1*h)/(x1)**2
        else if(ii.eq.5)then
c-----
c       [ b k ] exp(-kh)
c-----
            a = 0.0
            b = y1 * exp(x1*h)/ x1
        else if(ii.eq.9)then
c-----
c       [ a + b k + c k^2  ]exp(-kh)
c-----
c        solve
c            | 1  x1   x1^2 | |a|   |y1|
c            | 1  x2   x2^2 | |b| = |y2|
c            | 1  x3   x3^2 | |c|   |y3|
c
c        subtract row 1 from row2 amd row 1 from row 2 
c-----
             eu1 = dexp(dble(x1*h))
             eu2 = dexp(dble(x2*h))
             eu3 = dexp(dble(x3*h))
             a11 = 1
             a12 = x1
             a13 = x1*x1
             a22 = x2 - x1
             a23 = x2*x2 - x1*x1
             a32 = x3 - x1
             a33 = x3*x3 - x1*x1
             b1 = eu1*y1
             b2 = eu2*y2 - eu1*y1
             b3 = eu3*y3 - eu1*y1
             ddet = a22*a33 - a23*a32
             db = (  a33*b2 - a23*b3)/ddet
             dc = ( -a32*b2 + a22*b3)/ddet
             da = b1 - x1*(db + x1*dc) 
             a = da
             b = db
             c = dc
        endif
        return
        end

        subroutine setlim(r,z,nk,dk,kmax,doasym,wvno1,wvno2,wvno3)
        real r, z, dk, kmax, wvno1, wvno2,wvno3
        integer nk
        logical doasym
c-----
c       provide the logic for the integration
c-----
        doasym = .true.
        if(abs(z).gt.r)then
             delk = 6.2831853/(32.*abs(z))
             kmax = 16.*3.1415927/abs(z)
        else 
             delk = 6.2831853/(32.*abs(r))
             kmax = 16.*3.1415927/abs(r)
        endif
        nk = max0(100, int(kmax/delk) )
        dk = kmax/nk
        wvno1 = 0.3*kmax
        wvno2 = 0.6*kmax
        wvno3 = 0.9*kmax
        WRITE(6,*)'DK=',dk,'  NK=',nk,' Kmax=',kmax, 'doasym=',doasym
        WRITE(6,*)'wvno1=',wvno1,' wvno2=',wvno2,' wvno3=',wvno3
        return
        end

        subroutine numerical(r,is,ir)
        real r
        integer is,ir
        parameter(NSOURCE=100,NRECEIVER=100,NSR=100)
        common/sourcesv/depthssv(NSOURCE)
        common/receivsv/depthrsv(NRECEIVER)

        common/num/numer(21)
        real numer
        real kmax
        common /whohlf/liswho, ldosurf
        logical liswho, ldosurf

        real gg(21)
        real g1(21), g2(21), g3(21)
        real dfac

        common/kint1/doasym
        logical doasym

        common/kint4/aa,bb,cc
            real aa(21),bb(21),cc(21)
        z = abs(depthrsv(ir)-depthssv(is))
        call setlim(r,z,nk,dk,kmax,doasym,wvno1,wvno2,wvno3)

        call setup(r,abs(z),1.0e-5)
c-----
c       if asymptotic then
c       evaluate the g functions twice and call solu
c-----
        if(doasym)then
           call getspc(g1,ir,is,wvno1)
           call getspc(g2,ir,is,wvno2)
           call getspc(g3,ir,is,wvno3)
           do i=1,16
             call solu(g1(i),g2(i),g3(i),wvno1,wvno2,wvno3,
     1         abs(z),i,aa(i),bb(i),cc(i))
           enddo
c-----
c          initialize the integrals
c          applying asymptotic if required
c-----
           call intini(numer,r)
        else
           do i=1,16
              numer(i) = 0.0
           enddo
        endif
c-----
c       now do the integration from large to small k
c-----
        do i=nk,1,-1
           wvno = (i)*dk - 0.5*dk
           call getspc(gg,ir,is,wvno)
           call hank(wvno*r,h0,h1,h2,h1kr,h2kr)
           if(doasym)then
               call gasym(gg,wvno,abs(z))
           endif
c ZDD
           numer(1)  = numer(1)  + dk*gg(1)*h0*wvno
c - RDD
           numer(2)  = numer(2)  + dk*gg(2)*h1*wvno
c ZDS
           numer(3)  = numer(3)  + dk*gg(3)*h1*wvno
           dfac = dk*(gg(4)+gg(13))*h1kr*wvno
c RDS
           numer(4)  = numer(4)  + dk*gg(4)*h0*wvno
     1           -dfac
c TDS
           numer(5)  = numer(5)  + dk*gg(13)*h0*wvno
     1           -dfac
c ZSS
           numer(6)  = numer(6)  + dk*gg(5)*h2*wvno
           dfac = 2.*dk*(gg(6)+gg(14))*h2kr*wvno
c RSS
           numer(7)  = numer(7)  + dk*gg(6)*h1*wvno
     1           -dfac
c TSS
           numer(8)  = numer(8)  + dk*gg(14)*h1*wvno
     1           -dfac
c ZEX
           numer(9)  = numer(9)  + dk*gg(7)*h0*wvno
c REX
           numer(10) = numer(10) + dk*gg(8)*h1*wvno
c ZVF
           numer(11) = numer(11) + dk*gg(9)*h0*wvno
c RVF
           numer(12) = numer(12) + dk*gg(10)*h1*wvno
c ZHF
           numer(13) = numer(13) + dk*gg(11)*h1*wvno
           dfac = dk*(gg(12)+gg(15))*h1kr*wvno
c RHF
           numer(14) = numer(14) + dk*gg(12)*h0*wvno
     1           -dfac
c THF
           numer(15) = numer(15) + dk*gg(15)*h0*wvno
     1           -dfac
c som
           numer(16) = numer(16) + dk*gg(16)*h0*wvno

        enddo
c-----
c       flip RDD REX RVF because of the J-1 = - J1
c-----
        numer(2) = - numer(2)
        numer(10) = - numer(10)
        numer(12) = - numer(12)
        
        return
        end
 
        subroutine  getspc(g,jr,js,wvno)
        real g(21)
        real wvno
        integer ir,is
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
        common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
        real bsh, rhosh, qbsh
        common/jout/jsrc(21) , jbdrys, jbdryh
        integer jsrc, jbdrys, jbdryh
c-----
        common/lyrctl/lyrins
        logical lyrins

        if(.not.lyrins)then
c-----
c              evaluate matrices first
c              for currently defined layering
c-----
               call modcpy(.false.)
               call insert(depths(js))
               call insert(depthr(jr))
               call srclyr(depths(js),
     1             lmaxs(js), dphs)
               call srclyr(depthr(jr),
     1             lmaxr(jr), dphr)
               call dezero()
               call equlck()
               call evlmat(dble(wvno),jbdrys,jbdryh)
        endif
        call srshof(g,dble(wvno), lmaxs(js), lmaxr(jr)) 
        return
        end

        subroutine hank(z,h0,h1,h2,h1z,h2z)
c-----
c       evaluate Bessel functions using Abromiwitz and Stegun
c       this returns J0, J1, J2 J1/r and J2/r
c       we use h0 etc to avoid linking problems with
c       the math library
c-----
        real h0, h1, h2, h1z, h2z
        real z
        real j0, j1
        real j1z
c-----
        h0 = 0.0
        h1 = 0.0
        if(z.eq.0.0)then
            h0 = 1.0
            h1 = 0.0
            h1z = 0.5
            h2 = 0.0
            h2z = 0.0
        elseif(z.gt.0.0 .and. z.le.3.0)then
            x = (z/3.)*(z/3.)
            j0 = 1.-x*(2.2499997-x*(1.2656208-x*(.3163866-x*(
     1            .0444479-x*(.0039444-x*(.0002100))))))
            j1z = 0.5-x*(.56249985-x*(.21093573-x*(.03954289-x*(
     1      .00443319-x*(.00031761-x*(.00001109))))))
            j1 = z * j1z
            h0 = j0
            h1 = j1
            h1z = j1z
            h2 = 2.*h1/z - h0
            h2z = h2/z
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
            h0 =       fac * f0 * cos(t0)
            h1 =       fac * f1 * cos(t1)
            h1z = h1/z
            h2 = 2.*h1/z - h0
            h2z = h2/z
        endif
        return
        end

        subroutine gcmdln(outtxt)
c-----
c       parse command line arguments
c       requires subroutine mgtarg() and function mnmarg()
c-----
c       outtxt  L   - .true. create the file
c                     hstat96.txt in addition to the file96c
c-----
        logical outtxt
        common/verby/verbose
        logical verbose

        integer*4 mnmarg
        character*50 name

        outtxt = .false.

        verbose = .false.
        nmarg=mnmarg()
        i = 0
   11   i = i + 1
            if(i.gt.nmarg)goto 13
            call mgtarg(i,name)
            if(name(1:2).eq.'-?')then
                call usage(' ')
            else if(name(1:2).eq.'-h')then
                call usage(' ')
            else if(name(1:2).eq.'-V')then
                verbose = .true.
C            else if(name(1:4).eq.'-TXT')then
C            outtxt = .true.
            endif
        goto 11
   13   continue
        return
        end

        subroutine usage(str)
c------
c       write out program syntax
c-----
        integer LER,LOT,LIN
        parameter(LER=0, LIN=5, LOT=6)
        character str*(*)
        lstr = lgstr (str)
        write(LER,*)str(1:lstr)
        write(LER,*)'USAGE: ',
     1  'hstat96   [-?] [-h] '
        write(LER,*)
     1  '-V     (default .false.) Force  verbose output'
        write(LER,*)
     1  '-?                   Display this usage message'
        write(LER,*)
     1  '-h                   Display this usage message'
        stop 
        end
      
        subroutine srcpar(is)
c-----
c       given the source depth, get the parameters for
c       the TI specification needed by common rhdr96
c-----
        implicit none

        integer is

        integer NL
        integer NSOURCE, NRECEIVER, NSR
        parameter(NL=200)
        parameter(NSOURCE=100,NRECEIVER=100,NSR=100)
        common/source/depths(NSOURCE),lmaxs(NSOURCE),mdpths
        common/receiv/depthr(NRECEIVER),lmaxr(NRECEIVER),mdpthr
        real*4 depths, depthr
        integer lmaxs, lmaxr, mdpths, mdpthr

        common/rhdr96/ esec, distkm, distdg, evstaz, stevaz,
     1      evlat, evlon, evdep, stlat, stlon, stelev,
     2      TP, TSV, TSH, SA, SC, SF, SL, SN, SR
        real*4 esec, distkm, distdg, evstaz, stevaz
        real*4 evlat, evlon, evdep 
        real*4 stlat, stlon, stelev 
        REAL*4 TP, TSV, TSH
        REAL*4 SA, SC, SF, SL, SN, SR

        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        real d,a,b, rho,qa,qb,etap,etas,frefp,frefs
        common/depref/refdep
        real refdep
        common/modlly/mmax
        integer mmax

        real dphs

        real VSA, VSB, VSR

        integer i

        call modcpy(.false.)
        call insert(depths(is))
        call srclyr(depths(is),
     1     lmaxs(is), dphs)
        
        VSA = a(lmaxs(is))
        VSB = b(lmaxs(is))
        VSR = rho(lmaxs(is))

            SR = VSR
            SL = VSR*VSB*VSB
            SN = VSR*VSB*VSB
            SA = VSR*VSA*VSA
            SC = VSR*VSA*VSA
            SF = SA - 2.0*SN
        

        return 
        end
