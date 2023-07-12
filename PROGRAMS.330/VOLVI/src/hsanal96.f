
        program hsanal96
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME VI                                                      c
c                                                                     c
c      PROGRAM: HSANAL96                                               c
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
        character*3 istat3
        logical ixst3
        real*4 ffreq(8)
        parameter(NDST=500)
        real*4 r(NDST), tshift(NDST), vred(NDST)
        real smm(NSR,21)
        real stmp, sdata
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

        common/anal/analy(21)
        real analy
        common/cart/carte(21)
        real carte

        common /whohlf/liswho, ldosurf
        logical liswho, ldosurf

        logical ext
        character mname*80, title*80

        common/lwater/lfluid
        logical lfluid

        common/depref/refdep
        real refdep

        common/earth/radius
        real radius


        real*8 gg(21)

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
c       occupies position 10, or jsrc(lsrc(10)) = computed
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
c       define whether wholespace or halfspace
c-----
        if(jbdrys .eq.1)then
           ldosurf = .true.
           liswho = .false.
        else
           ldosurf = .false.
           liswho = .true.
        endif
        
c-----
c       process
c-----
c-----
c       UNIX output - no carriage control
c-----
    2   format('fl =',f10.5,5x,'fu =',f10.5,5x,'df =',f10.5,/
     1      4x,'n1 =',i4,5x,'n2 =',i4,5x, ' n =',i5/
     2  'vmin  =',f10.5,' vamin =',f10.5,' vamax =',f10.5)
   21   format(
     1  '       ',10x  ,' vbmin =',f10.5,' vbmax =',f10.5)
 2021   format('depths =',f20.6,' (',f20.6,')')
 2031   format('depthr =',f14.6,' (',f20.6,')')
 2041   format('     r =',f14.6,' tshift=',f10.2,' vred=',f10.2)
    3   format(8f10.5)
    4   format('frequencies for which response computed     ')
    5   format('alpha =',f10.5,5x,'dt =',f10.3)
 2042   format('Hankel function used for kr >=',f10.2)
c-----
c       MSDOS MICROSOFT FORTRAN OUTPUT
c-----
c    2   format(1x,'fl =',f10.5,5x,'fu =',f10.5,5x,'df =',f10.5,/
c     1          1x,4x,'n1 =',i4,5x,'n2 =',i4,5x, ' n =',i5/
c     2  1x,'vmin  =',f10.5,' vamin =',f10.5,' vamax =',f10.5)
c   21  format(
c     1 1x,'       ',10x  ,' vbmin =',f10.5,' vbmax =',f10.5)
c 2021   format(1x,'depths =',f20.6)
c 2031   format(1x,'depthr =',f14.6)
c 2041   format(1x,'     r =',f14.6,' tshift=',f10.2,' vred=',f10.2)
c    3   format(1x,8f10.5)
c    4   format(1x,'frequencies for which response computed     ')
c    5   format(1x,'alpha =',f10.5,5x,'dt =',f10.3)
c 2042  format(1x, 'Hankel function used for kr >=',f10.2)
c-----
        if(verbose)then

            write(LOT,*)'SOURCE DEPTH IN WORKING AND ORIGINAL MODEL (',
     1           mdpths,')'
            do i=1,mdpths
                write(LOT,2021)depths(i), depthssv(i)
            enddo
            write(LOT,*)'RECEIVER DEPTH IN WORKING AND ORIGINAL MODEL ('
     1           ,mdpthr,')'
            do i=1,mdpthr
                write(LOT,2031)depthr(i), depthrsv(i)
            enddo
            write(LOT,*)'RECEIVER DISTANCES (',ndist,')'
            do i=1,ndist
                write(LOT,2041)r(i), tshift(i), vred(i)
            enddo
            write(LOT,4)
            WRITE(LOT,*)'==================='
        endif
c-----
c     open the output file file96
c-----
      open(8,file='file96',status='unknown',form='formatted',
     1 access='sequential')
      rewind 8
      do i=1,16
         jjsrc(i) = 1
      enddo
      jjsrc(16) = 0

      do jd=1,ndist
         k = 0
         do is=1,mdpths
c-----
c     for this program the source and receiver
c     density and velocity are for laeyr 1
c-----
            VSA = a(1)
            VSB = b(1)
            VSR = rho(1)

            SR = VSR
            SL = VSR*VSB*VSB
            SN = VSR*VSB*VSB
            SA = VSR*VSA*VSA
            SC = VSR*VSA*VSA
            SF = SA - 2.0*SN

            do ir=1,mdpthr
               cpulse='hstat96'
               iunit = 2
               iftype = 16
               iobsyn = 2
               itmfrq = 1
               cfilt = 'hsanal96'
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
               if(verbose)then
                  WRITE(LOT,*)'r=',r(jd),' HS=',depths(is), 
     1               'HR=',depthr(ir)
               endif
               call hlfwho(r(jd),depths(is),depthr(ir),ldosurf)
               do jj=1,15
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
                  if(verbose)then
                     if(ldosurf)then
                        WRITE(LOT,*)'TSTATIC:',
     1                     jj,ost(jj)(1:3),lsrc(jj),analy(jj)
                      else
                        WRITE(LOT,*)'TSTATIC:',jj,ost(jj)(1:3),
     1                     lsrc(jj),analy(jj),carte(jj)
                      endif
                   endif
                   npts = 1
                   call wrtr96(8,ost(jj)(1:3),cmpinc,cmpaz,cmpdt,npts,
     1                ksyear, ksmon, ksday, kshour, ksmin, ssec,
     2                analy(jj), nerr, 1)

               enddo
            enddo
         enddo
      enddo
      close(8)
 9999 continue
      end

       subroutine hlfwho(r,depths,depthr,ldosurf)
       real r, deptths, depthr
       logical ldosurf
c-----
c      use downward positive z axis
c      depth negative means receiver is about source
c-----
       depth = - (depths - depthr)
       if(ldosurf)then
          call analytic(r,depth)
       else
          call cartesian(r,depth)
          call analytic(r,depth)
       endif
       return
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
        if(mdpths+mdpthr .gt. 2*mmax )then
            lyrins = .false.
            if(verbose)write(LOT,*)' LAYER INSERTION NOT DONE'
        else
            lyrins = .true. 
            if(verbose)write(LOT,*)' LAYER INSERTION DONE'
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
        integer ifunc

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
     1             bsh(mmax),rhosh(mmax)
        endif
c-----
c     obtain extreme velocity limits
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
C            if(name(1:4).eq.'-TXT')then
C                outtxt = .true.
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
     1  'hsanal96  [-?] [-h] '
        write(LER,*)
     1  '-V     (default .false.) Force verbose output'
        write(LER,*)
     1  '-?                   Display this usage message'
        write(LER,*)
     1  '-h                   Display this usage message'
        stop 
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

        subroutine analytic(r,z)
        implicit none
        real r, z
        common/anal/analy(21)
        real analy
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        real d,a,b,rho,qa,qb,etap,etas,frefp,frefs
        common/modlly/mmax
        integer mmax
        common /whohlf/liswho, ldosurf
        logical liswho, ldosurf

        real j0k0, j0k1, j0k2,j1k0,j1k1,j1k2,j2k0,j2k1,j2k2
        real j1k0r, j1k1r, j1km1r, j2k0r, j2k1r

        real*8 rs1,rs2, rs3, rs5, rs7,  absz, z2, z3, rsph
        real*8 pi4r
        real*8 ba2, ba2m1, h, a2, b2
        real*8 fac
        integer i
        real rlim
        common/verby/verbose
        logical verbose

        if(verbose)then
           WRITE(6,*)'analytic: r=',r,' z=',z
           WRITE(6,*)'a=',a(1),' b=',b(1),' rho=',rho(1)
        endif
        rlim = 1.0 e-07
        rsph = sqrt(r*r + z*z)
c------
c       never permit a zero distance
c------
        if(rsph .lt. rlim)then
           rs1=rlim
        else
           rs1 = rsph
        endif
        rs2 = rsph**2
        rs3 = rsph**3
        rs5 = rsph**5
        rs7 = rsph**7
        absz = abs(z)
        z2 = z*z
        z3 = z*z2
        pi4r = 1./(4.*3.1415927*rho(1))
        do i=1,21
          analy(i) = 0.0
        enddo
c-----
c     WHOLESPACE
c-----
        if(liswho)then
        analy( 1) =  pi4r*(  2.*z/(a(1)**2*RS3) - 
     1     (3./2.)*(1./a(1)**2 -1./b(1)**2)
     2     *(3.*z*z*z/RS5 - z/RS3))
        analy( 2) =  pi4r*(  2.*r/(a(1)**2*RS3) - 
     1     (3./2.)*(1./a(1)**2 -1./b(1)**2)
     2     *(3.*r*z*z/RS5 ) - 1.5*(1./a(1)**2 + 1./b(1)**2)*r/RS3 )
        analy( 3) = -pi4r * ( 
     1       (1./a(1)**2 - 1./b(1)**2) * 3*r*z*z/RS5  
     2       - ( 1 / a(1)**2 + 1 / b(1)**2)* r/RS3
     3       + 1./b(1)**2 * r/RS3
     4        )
        analy( 4) = pi4r *( 
     1      (1./a(1)**2 - 1./b(1)**2)*(3.*z*z*z/RS5 - 2.*z/RS3) 
     2      + 1 / b(1)**2 * z/RS3
     3      )
        analy( 5) =  - pi4r * (1./a(1)**2) * z/RS3
        analy( 6) = - pi4r* 0.5* (1./a(1)**2 - 1./b(1)**2)
     1     * 3*r*r*z/RS5
        analy( 7) = pi4r * 0.5 * (
     1   (1./a(1)**2 - 1./b(1)**2) 
     2     *( 3*z*z/RS2 - 2)
     3  + (1./a(1)**2 + 1./b(1)**2)
     1  ) * r / RS3
        analy( 8) = -pi4r * (
     1    (1./a(1)**2) *r/RS3 
     1    )
        analy( 9) =  pi4r * (1./a(1)**2) * z/RS3
        analy(10) =  pi4r * (1./a(1)**2) * r/RS3
        analy(11) =  pi4r * 0.5 *(1./a(1)**2 + 1./b(1)**2)/RS1
     1    - pi4r * 0.5 * (1./a(1)**2 - 1./b(1)**2) * z*z/RS3
        analy(12) = - pi4r*0.5* (1./a(1)**2 - 1./b(1)**2) * r * z/RS3
        analy(13) = - pi4r*0.5*(1./a(1)**2 - 1./b(1)**2) * r * z/RS3
        analy(14) = pi4r * 0.5 *(
     1    - (1./a(1)**2 - 1./b(1)**2) *r*r/RS3 
     1    + (1./a(1)**2 + 1./b(1)**2)/RS1 
     1 )
        analy(15) = - pi4r * 0.5 *(1 / a(1)**2 + 1./b(1)**2) / RS1
        else
c-----
c       HALFSPACE SURFACE ANALYTIC
c-----
        pi4r = (4.*3.1415927*rho(1))
        h = abs(z)
        absz = h
        b2 = b(1)*b(1)
        a2 = a(1)*a(1)
        ba2 = b2/a2
        ba2m1 = ba2 -1
c------
c        special care for  for short epicentral distances
c-----
        if(r.lt.rlim)then
              j0k0 = 1./RS1
              j0k1 = h/RS3
              j0k2 = (2.*h*h -r*r)/RS5
              j1km1r = 1./(RS1 + absz)
              j1k0 = 0.0
              j1k0r = 0.5/(absz*absz)
              j1k1 = 0.0
              j1k1r = 1/RS3
              j1k2 = 0.0
              j2k0 = 0.0
              j2k0r = 0.0
              j2k1 = 0.0
              j2k1r = 0.0
              j2k2 = 0.0

        else
              j0k0 = 1./RS1
              j0k1 = h/RS3
              j0k2 = (2.*h*h -r*r)/RS5
              j1km1r = 1./(RS1 + absz)
              j1k0 = r/(RS1*(RS1 + absz))
              j1k0r = 1/(RS1*(RS1 + absz))
              j1k1 = r/RS3
              j1k1r = 1/RS3
              j1k2 = 3.*r*h/RS5
              j2k0 = (1.-2.*(h/RS1) + (h/RS1)**2) *RS1/(r*r)
              j2k0r = (1.-2.*(h/RS1) + (h/RS1)**2) *RS1/(r*r*r)
              j2k1 = (2.-3.*(h/RS1) + (h/RS1)**3)     /(r*r)
              j2k1r = (2.-3.*(h/RS1) + (h/RS1)**3)    /(r*r*r)
              j2k2 = 3.*r*r/RS5
        endif
c-DD
        analy( 1) = (ba2)*j0k1/(pi4r*b2*ba2m1) -
     1       3.* absz*ba2m1*j0k2/(pi4r*b2*ba2m1)
        analy( 2) = (3.-4.*ba2)*j1k1/(pi4r*b2*ba2m1)
     1     +3.*absz*ba2m1*j1k2/(pi4r*b2*ba2m1)
c-DS
        analy( 3) = 2.*absz*j1k2/(pi4r*b2)
        analy( 4) = -2.*(j0k1 - absz*j0k2)/(pi4r*b2)
     1    - (1.)*(-2.*j1k0r +2.*absz*j1k1r)/(pi4r*b2)
     1    - (1.)*(2.*j1k0r)/(pi4r*b2)
C     1    - (1./r)*(-2.*j1k0 +2.*absz*j1k1)/(pi4r*b2)
C     1    - (1./r)*(2.*j1k0)/(pi4r*b2)
C analytically TDS = 0 at the surface
        analy( 5) = 2.*(j0k1) /(pi4r*b2)
     1    - (1.)*(-2.*j1k0r +2.*absz*j1k1r)/(pi4r*b2)
     1    - (1.)*(2.*j1k0r)/(pi4r*b2)
C     1    - (1./r)*(-2.*j1k0 +2.*absz*j1k1)/(pi4r*b2)
C     1    - (1./r)*(2.*j1k0)/(pi4r*b2)
c-SS
        analy( 6) = -ba2*j2k1/(pi4r*b2*ba2m1)
     1      -absz*ba2m1*j2k2/(pi4r*b2*ba2m1)
        analy( 7) = -j1k1/(pi4r*b2*ba2m1) - absz*j1k2/(pi4r*b2)
     1    -(2.)*(- j2k0r - ba2m1*absz*j2k1r)/(pi4r*b2*ba2m1)
     2    -(2.)*(-2.*j2k0r)/(pi4r*b2)
C     1    -(2./r)*(- j2k0 - ba2m1*absz*j2k1)/(pi4r*b2*ba2m1)
C     2    -(2./r)*(-2.*j2k0)/(pi4r*b2)
        analy( 8) = -2.*j1k1/(pi4r*b2)
     1    -(2.)*(- j2k0r - ba2m1*absz*j2k1r)/(pi4r*b2*ba2m1)
     2    -(2.)*(-2.*j2k0r)/(pi4r*b2)
C     1    -(2./r)*(- j2k0 - ba2m1*absz*j2k1)/(pi4r*b2*ba2m1)
C     2    -(2./r)*(-2.*j2k0)/(pi4r*b2)
c-EX
        analy( 9) = 2.*j0k1/(pi4r*(b2-a2))
        analy(10) = -2.*j1k1/(pi4r*(b2-a2))
c-VF
        analy(11) = - j0k0/(pi4r*b2*ba2m1)
     1      +absz*j0k1/(pi4r*b2)
        analy(12) =  ba2*j1k0/(pi4r*b2*ba2m1)
     1      -absz*j1k1/(pi4r*b2)
c-HF
        analy(13) = - ba2*j1k0/(pi4r*b2*ba2m1)
     1      -absz*ba2m1*j1k1/(pi4r*b2*ba2m1)
C        fac = -(1./r)*(1./(pi4r*b2))*(-2.*j1km1
C     1      - (1./ba2m1)*j1km1 - absz*j1k0)
        fac = -(1./(pi4r*b2))*(-2.*j1km1r
     1      - (1./ba2m1)*j1km1r - absz*j1k0r)
        analy(14) = - (j0k0 + absz*ba2m1*j0k1)/(pi4r*b2*ba2m1)
     1     + fac
        analy(15) = -2.*j0k0/(pi4r*b2)
     1     + fac
        
        continue
       
        endif
c-----
c       invert all Z s to make compatible with CPS usage 
c       that z-positive is up
c-----
        analy( 1) = -analy(1)
        analy( 3) = -analy(3)
        analy( 6) = -analy(6)
        analy( 9) = -analy(9)
        analy(11) = -analy(11)
        analy(13) = -analy(13)


        return
        end
  
        
        subroutine cartesian(r,z)
        implicit none
        real r, z
        common/cart/carte(21)
        real carte
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        real d,a,b,rho,qa,qb,etap,etas,frefp,frefs
        common/modlly/mmax
        integer mmax
        common /whohlf/liswho, ldosurf
        logical liswho, ldosurf


        integer i
        real rlim
        real rsph
        real gij, gijk
c-----
c       direction cosine for the observation point
c-----
        real gamm(3)
c-----
c       azimuthal angle to observation point, taken here as
c       in x1 direction
c-----
        real cost, sint

        cost = 1.0
        sint = 0.0

        rsph = sqrt(r*r + z*z)
c------
c       never permit a zero distance
c------
           gamm(1) = r/rsph
           gamm(2) = 0./rsph
           gamm(3) = z/rsph
        do i=1,21
          carte(i) = 0.0
        enddo
c-----
c     WHOLESPACE
c-----
        if(liswho)then
c-----
c       analytic solution in cartesian coordinates
c       note that in this case that the gij,k are not the full moment 
c       tensor. For the off-diagonal I must include the effect of
c       the force symmetry of the Mij. For a rece3iver azimuth in the x1 direction,
c       e.g., cost = 1, this affectsa ZDS, RDS, TSS, TDS
c-----
c-----
c   ZDD
c-----
           carte( 1) = -gijk(3,1,1,a(1),b(1),rho(1),rsph,gamm) 
     1            -gijk(3,2,2,a(1),b(1),rho(1),rsph,gamm) +
     2            2.*gijk(3,3,3,a(1),b(1),rho(1),rsph,gamm) 
c-----
c   RDD
c-----
           carte( 2) = -gijk(1,1,1,a(1),b(1),rho(1),rsph,gamm) 
     1            -gijk(1,2,2,a(1),b(1),rho(1),rsph,gamm) +
     2            2.*gijk(1,3,3,a(1),b(1),rho(1),rsph,gamm) 
c-----
c   ZDS
c-----
           carte( 3) =  gijk(3,3,1,a(1),b(1),rho(1),rsph,gamm) 
     1                + gijk(3,1,3,a(1),b(1),rho(1),rsph,gamm)
c-----
c   RDS
c-----
           carte( 4) =  gijk(1,1,3,a(1),b(1),rho(1),rsph,gamm) 
     1                + gijk(1,3,1,a(1),b(1),rho(1),rsph,gamm)
c-----
c   TDS
c-----
           carte( 5) = -(gijk(2,2,3,a(1),b(1),rho(1),rsph,gamm) 
     1                     + gijk(2,3,2,a(1),b(1),rho(1),rsph,gamm))
c-----
c   ZSS
c-----
           carte( 6) =  gijk(3,1,1,a(1),b(1),rho(1),rsph,gamm) 
     1                 -gijk(3,2,2,a(1),b(1),rho(1),rsph,gamm)
c-----
c   RSS
c-----
           carte( 7) =  gijk(1,1,1,a(1),b(1),rho(1),rsph,gamm) 
     1                 -gijk(1,2,2,a(1),b(1),rho(1),rsph,gamm)
c-----
c   TSS
c-----
           carte( 8) = -(gijk(2,1,2,a(1),b(1),rho(1),rsph,gamm) 
     1                    + gijk(2,2,1,a(1),b(1),rho(1),rsph,gamm))
c-----
c   ZEX
c-----
           carte( 9) = gijk(3,1,1,a(1),b(1),rho(1),rsph,gamm) +
     1            gijk(3,2,2,a(1),b(1),rho(1),rsph,gamm) +
     2            gijk(3,3,3,a(1),b(1),rho(1),rsph,gamm) 
c-----
c   REX
c-----
           carte(10) = gijk(1,1,1,a(1),b(1),rho(1),rsph,gamm) +
     1            gijk(1,2,2,a(1),b(1),rho(1),rsph,gamm) +
     2            gijk(1,3,3,a(1),b(1),rho(1),rsph,gamm) 
c-----
c   ZVF
c-----
           carte(11) = gij(3,3,a(1),b(1),rho(1),rsph,gamm)
c-----
c   RVF
c-----
           carte(12) = gij(1,3,a(1),b(1),rho(1),rsph,gamm)
c-----
c   ZHF
c-----
           carte(13) = gij(3,1,a(1),b(1),rho(1),rsph,gamm)
c-----
c   RHF =   g11*cost + g21*sint
c-----
           carte(14) = gij(1,1,a(1),b(1),rho(1),rsph,gamm)
c-----
c   THF  = -g11*sint + g21*cost
c-----
           carte(15) = -gij(2,2,a(1),b(1),rho(1),rsph,gamm)
        else
            continue
        endif
c-----
c       invert all Z s to make compatible with CPS usage that z-pos is up
c-----
        carte( 1) = -carte(1)
        carte( 3) = -carte(3)
        carte( 6) = -carte(6)
        carte( 9) = -carte(9)
        carte(11) = -carte(11)
        carte(13) = -carte(13)


        return
        end
  
        real function  gij(i,j,a,b,rho,r,g)
        integer i, j
        real a, b, rho, r
        real g(3)
        real delta
        gij = (0.5*(g(i)*g(j) - delta(i,j))*(1./(b*b) - 1./(a*a))
     1         + delta(i,j)/(b*b))/(4.*3.1415927*rho*r)
        return
        end
  
        real function  gijk(i,j,k,a,b,rho,r,g)
        integer i, j, k
        real a, b, rho, r
        real g(3)
        real delta
        gijk = (( g(i)*delta(j,k) + g(j)*delta(i,k)+g(k)*delta(i,j)
     1          - 3.*g(i)*g(j)*g(k))/(2.*a*a)
     2        - (g(i)*delta(j,k)+g(j)*delta(i,k)-g(k)*delta(i,j)
     1          - 3.*g(i)*g(j)*g(k))/(2.*b*b) 
     1        ) /(4.*3.1415927*rho*r*r)
        return
        end
 
 
        function delta(i,j)
        real delta
             if(i.eq.j)then
                 delta = 1.0
             else 
                 delta = 0.0
             endif
        return
        end
