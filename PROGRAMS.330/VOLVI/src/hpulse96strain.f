        program hpulse96strain
c---------------------------------------------------------------------c
c                                                                     c
c        COMPUTER PROGRAMS IN SEISMOLOGY                              c
c        VOLUME V                                                     c
c                                                                     c
c        PROGRAM: HPULSE96STRAIN                                      c
c                                                                     c
c        COPYRIGHT 1996 R. B. Herrmann                                c
c                                                                     c
c        Department of Earth and Atmospheric Sciences                 c
c        Saint Louis University                                       c
c        221 North Grand Boulevard                                    c
c        St. Louis, Missouri 63103                                    c
c        U. S. A.                                                     c
c                                                                     c
c---------------------------------------------------------------------c
c       CHANGES to hspec96
c       11 SEP 2000 - build in P, SV and SH first arrival times
c       14 JAN 2001 - put in A, C, F, L and N constants into
c           trace header
c       05 FEB 2004 - modified to be slightly more tolerant about DT for
c           user supplied pulse - now issues WARNING and not termination
c           mlarocca@ov.ingv.it
c       07 FEB 2005 - add a -Z flag to indicate that the 
c           internal parabolic 
c           or triangular pulses are to be zero phase 
c       12 DEC 2007 - set header values evlat, evlon, stlat, stlon to -12345
c           for compatibility of resulting SAC traces
c       23 MAR 2008 - define iprog = 4 for hudson96.f
c       23 JUL 2009 - ensured consistent call to FFT routine with
c               complex arguments. We now call zfour
c       05 JUL 2012 - changed the IO for the user defined pulse to be
c                   list directed IO instead of a specific format
c
c       CHANGES to hpulse96strain.f
c       03 MAR 2021    modified hpulse96 to create sac files of strain 
c             for a given source. This works together with hspec96strain
c             
c             One problem with computing strain is that we must be 
c             very careful about the units. Since we are creating the
c             synthetics and the strain from surface-wave modal superposition,
c             we will need terms such as dU/dr and dU/dz at the receiver.
c             
c             
c             The UZ UR UT  depend  on the units of the model. Typically the model 
c             is given as KM KM/S GM/CM^3 in which case the these are
c             CM, CM/S or CM/S/S.  To compute strain we must worry about
C             CM/KM 
c
c             Because of this we will have to be careful about the strains, e.g.,
c             The units will be DIsplacement (which is another question)/KM
c
c             So for a moment of 1.0e+20 the output will be cm so the strain will be
c             multiplied 0.01 / 1000  = 1.0e+6
c
c            Note the original spulse96 was defined to have Uz positive up
c            we must be careful of this sign change for the comptuation of strain
c            which uses a coordinate system of  z positive down
c
c       24 JAN 2022    corrected flags -ZZ was interpreted at -Z which
c            gave the wrong moment and forced zero phase
c                 
c            Also note a problem with the predicted P time for a fluid
c            This was correctly defined in hspec96strain. The problem
c            occured since I needed to compute material parameters
c            at the receiver to get correct stresses. Thus the
c            calls to outtest1 outstrain outstress outrotate outgreen
c            now use the TP TSV and TSH from hspec96strain and not from
c            the local call to FRSTAR
c
c-----

c-----
c       command line arguments
c-----
        character stcomp*8, stname*8
        integer*4 ntau, ipt, idva
        real*4 alp
        character*80 rfile
        logical dozero, dotest1
        logical dostep, dosdr

        common/grnfmt/outfmt
        integer outfmt

        common/stresstrain/dostrain,dostress,dorotate,dogreen,
     1     az, baz, xmt(3,3), xmom,fx, fy,fz,
     1     caz, saz, c2az, s2az, stk, dip, rake
        logical dostrain, dostress,dorotate,dogreen
        real az, baz, xmt, xmom,fx, fy,fz
        real caz, saz, c2az, s2az, stk, dip, rake

c-----
c       variables local to the program
c-----
        integer LER, LIN, LOT
        parameter (LER=6, LIN=5, LOT=6)
        integer NGRN
        parameter (NGRN=45)
        integer NSAMP, NFREQ
        parameter (NSAMP=16384,NFREQ=8193)
        common/xx/x
        real*4 x(NSAMP)
        common/zzx/zdata
        complex zdata(NSAMP)
        common/zzy/datas
        complex datas(NFREQ)
        integer*4 ksrc(NGRN)
        common/srctim/ src(NSAMP)
        real src
        complex xx(NGRN)
        real ar, ai
        character*80 mname
        logical ext
        complex spec(NSAMP,NGRN)
        integer*4 iszrt(45)
        character*8 ost(45)
        complex ctmp
        data iszrt/1,4, 1,4,5, 1,4,5, 1,4, 1,4, 1,4,5,
     1             1,4, 1,4,5, 1,4,5, 1,4, 1,4, 1,4,5,
     1             1,4, 1,4,5, 1,4,5, 1,4, 1,4, 1,4,5/
        data ost(1:15)/'ZDD     ', 'RDD     ', 'ZDS     ', 'RDS     ',
     1  'TDS     ', 'ZSS     ', 'RSS     ', 'TSS     ','ZEX     ', 
     1  'REX     ', 'ZVF     ', 'RVF     ', 'ZHF     ', 'RHF     ', 
     1  'THF     ' /
        data ost(16:30)/'ZDDz    ', 'RDDz    ', 'ZDSz    ', 'RDSz    ',
     1  'TDSz    ', 'ZSSz    ', 'RSSz    ', 'TSSz    ','ZEXz    ', 
     1  'REXz    ', 'ZVFz    ', 'RVFz    ', 'ZHFz    ', 'RHFz    ', 
     1  'THFz    ' /
        data ost(31:45)/'ZDDr    ', 'RDDr    ', 'ZDSr    ', 'RDSr    ',
     1  'TDSr    ', 'ZSSr    ', 'RSSr    ', 'TSSr    ','ZEXr    ', 
     1  'REXr    ', 'ZVFr    ', 'RVFr    ', 'ZHFr    ', 'RHFr    ', 
     1  'THFr    ' /

c-----
c       call machine dependent initialization
c-----
        call mchdep()
c-----
c       parse command line arguments
c-----
        call gcmdln(ntau,ipt,alp,idva,
     1      rfile,
     2      dozero,
     1      dostrain,dostress,dorotate,dogreen,
     1      az,baz,stk,dip,rake,xmom,xmt,fx,fy,fz,
     1      dotest1,dostep,outfmt,dosdr)
        WRITE(LER,*)'ntau              :',ntau
        WRITE(LER,*)'ipt               :',ipt 
        WRITE(LER,*)'alp               :',alp 
        WRITE(LER,*)'idva              :',idva
        WRITE(LER,*)'rfile             :',ntau
        WRITE(LER,*)'dozero            :',dozero
        WRITE(LER,*)'dostrain          :',dostrain
        WRITE(LER,*)'dostress          :',dostress
        WRITE(LER,*)'dorotate          :',dorotate
        WRITE(LER,*)'dogreen           :',dogreen 
        WRITE(LER,*)'dotest1           :',dotest1
        WRITE(LER,*)'az                :',az
        WRITE(LER,*)'baz               :',baz
        WRITE(LER,*)'xmom              :',xmom
        if(dosdr)then
        WRITE(LER,*)'strike            :',stk
        WRITE(LER,*)'dip               :',dip   
        WRITE(LER,*)'rake              :',rake   
        else
           WRITE(LER,*)'xmt               :',xmt(1,1),xmt(1,2),xmt(1,3)
           WRITE(LER,*)'                  :',xmt(2,1),xmt(2,2),xmt(2,3)
           WRITE(LER,*)'                  :',xmt(3,1),xmt(3,2),xmt(3,3)
        endif
        WRITE(LER,*)'forcex            :',fx
        WRITE(LER,*)'forcey            :',fy
        WRITE(LER,*)'forcez            :',fz
        WRITE(LER,*)'dostep            :',dostep
        WRITE(LER,*)'outfmt            :',outfmt
c-----
c       get the azimuths for the source
c-----
        degrad = 3.1415927/180.0
        caz = cos(degrad*az)
        saz = sin(degrad*az)
        c2az = cos(2*degrad*az)
        s2az = sin(2*degrad*az)
c-----
c       ensure the existence of the Green s functions generated by
c       hspec96(VI)
c-----
        inquire(file='hspec96ss.grn',exist=ext)
        if(.not. ext)then
                write(LER,*)'Green s function file ',
     1          'hspec96ss.grn does not exist'
                go to 9000
        endif
        open(4,file='hspec96ss.grn',access='sequential',
     1      form='unformatted', status='unknown')
        rewind 4
c-----
c       open scratch temporary file
c-----
        open(unit=10,status='scratch',form='unformatted')
        rewind 10
c-----
c       get header information
c-----
        read(4,end=9999,err=9999)iprog
        read(4,end=9999,err=9999) alpha,fl,fu,dt,n,n1,n2,df,nyq2
        read(4,end=9999,err=9999)mname
C       WRITE(0,*)iprog
C       WRITE(0,*)alpha,fl,fu,dt,n,n1,n2,df,nyq2 
C       WRITE(0,*)MNAME

c-----
c       define the source pulse whose spectra is stroed in datas
c-----
        call srcpul(ntau,ipt,alp,dt,n,rfile,duration,alpha)
        if(.not. dozero)then
            duration = 0.0
        endif

c-----
c       start reading Green s functions. However, only define the pulses
c       with the first read. The data are in the following order:
c           DIST
c               SOURCE_DEPTH
c                   RECEIVER_DEPTH
c                       FREQ
c                           GRN01 .. GRN21 .. GRN45
c-----
 1000   continue
            read(4,err=9999,end=9999)r,t0,depths,depthr,
     1          TP,TSV,TSH, 
     2          SA, SC, SF, SL, SN, SR
            WRITE(0,*)'r,t0,depths,depthr,TP,TSV,TSH:',
     1          r,t0,depths,depthr,TP,TSV,TSH
            if(r.lt.0.0) go to 9999 
            read(4,end=9999,err=9999)ksrc

c-----
c           get Green s functions for this distance, 
c           apply the source time function
c           which will have to be demultiplexed later
c-----
            rewind 10
            do i=n1,n2
                do j=1,NGRN
                    xx(j)=cmplx(0.0,0.0)
                    if(ksrc(j).ne.0)then
                        read(4,err=9999,end=9999)ar,ai
                            xx(j) = cmplx(ar,ai)*datas(i)
                    endif
                enddo
c-----
c        The output of hspec96strain  is GRN, d/dz GRN, d/dr GRN
c        The output of spulse96strain is GRN, d/dr GRN, d/dz GRN
c        To use the same subroutines here for both hpulse96strain and spulse96strain
c        reorder to make the order GRN d/dr GRN and d/dz GRN
c-----
                do k=1,15
                     ctmp = xx(k+15)
                     xx(k+15) = xx(k+30)
                     xx(k+30) = ctmp
                enddo
                write(10)(xx(k),k=1,NGRN)
            enddo
c-----
c     the input order from hspec96strain is
c      1  ZDD   16  d ZDD/dz   31 d ZDD/dr
c      2  RDS
c      3  ZDS
c      4  RDS
c      5  ZSS
c      6  RSS
c      7  ZEX
c      8  REX
c      9  ZVF
c     10  RVF
c     11  ZHF
c     12  RHF
c     13  TDS
c     14  TSS
c     15  THF   30  d THF/dz   45 d THF/dr
c-----
c           output the functions 
c-----
            if(dotest1)then
c-----
c             this is for testing and outputs the Green's functions
c             and the d/dr and d/dz for the Green's functions. However
c             the Uz is positive up
c-----
                 call outtest1(ntau,idva,r,dt,n,
     1               t0-0.5*duration,depthr,depths,
     2               mname,dostep,alpha,TP,TSV,TSH)
            else
                 WRITE(0,*)'DURATION:',DURATION
                 if(dosdr)then
                 call  gtensor(mname,depths,depthr,r)
                 endif
                 if(dostrain)call outstrain(ntau,idva,r,dt,n,
     1              t0-0.5*duration,depthr,depths,az,
     2              mname,dostep,alpha,dosdr,TP,TSV,TSH)
                 if(dostress)then
                    call outstress(ntau,idva,r,dt,n,
     1              t0-0.5*duration,depthr,depths,az,
     2              mname,dostep,alpha,dosdr,TP,TSV,TSH)
                 endif
                 if(dorotate)then
                    call outrotate(ntau,idva,r,dt,n,
     1              t0-0.5*duration,depthr,depths,az,
     2              mname,dostep,alpha,dosdr,TP,TSV,TSH)
                 endif
                 if(dogreen)then
                 call outgreen(ntau,idva,r,dt,n,
     1               t0-0.5*duration,depthr,depths,
     2               mname,dostep,alpha,TP,TSV,TSH)
                 endif
            endif
        go to 1000
 9999   continue
        close (4)
        close (10)
 9000   continue
        end

        subroutine chtofp(str,fout)
c------
c       routine to convert string to floating point
c       The E format is accepted as well as free form
c       input
c
c       This assumes that the string is 20 characters long
c-----
        character*(*) str
        real*4 fout
        integer*4 lgstr
        logical hase
        l = lgstr(str)
c------
c       If the string str contains an E or e, then
c       we parse the result using an E format.
c------
        hase = .false.
        do 100 i=1,l
            if(str(i:i) .eq.'e' .or. str(i:i).eq.'E')then
                hase = .true.
            endif
  100   continue
c-----
c       read floating point number
c-----
        if(hase)then
            read(str,'(bn,e20.0)')fout
        else
            read(str,'(bn,f20.0)') fout
        endif
        return
        end

        subroutine gtensor(mname,depths,depthr,dist)
c-----
c       calculate moment tensor for a double couple mechanism
c       or for an explostion
c-----
c-----
c       Changes
c       13 APR 2021 - set up to work with TI media
c          Use TI coefficients A, C, F, L, N
c          normalize MIJ matrix then rescale 
c          with xmom defined as
c               |            2  | (1/2)
c          M  = | (1/2) SUM M   |
c           o   |            ij |
c          Silver, P. G. and T. H. Jordan (1982). Optimal estimation of scalar
c             seismi moemtn, Geophys. J. Roy. Astr.  Soc. 70, 755-787.
c-----
        implicit none
        integer LOT
        parameter (LOT=6)
        character mname*(*)
        real depths, depthr,dist

        common/stresstrain/dostrain,dostress,dorotate,dogreen,
     1     az, baz, xmt(3,3), xmom,fx, fy,fz,
     1     caz, saz, c2az, s2az, stk, dip, rake
        logical dostrain, dostress,dorotate,dogreen
        real az, baz, xmt, xmom,fx, fy,fz
        real caz, saz, c2az, s2az, stk, dip, rake

        real degrad, dp, st, rk
        real sind,cosd,sinr,cosr, sins, coss, sin2d, cos2d,sin2s,cos2s
        real s1, s2, s3, n1,n2,n3
        real tol, xmax, thresh
        real mijnorm
        real ttp
        real SA, SC, SF, SL, SN, SR
        real RA, RC, RF, RL, RN, RR
        real rayp, geom, tstar
        integer iiso

        integer i,j
c-----
c       get the medium parameters by a call to travel time
c-----
            CALL FRSTAR(dist,depths,depthr,
     1          MNAME,1,tTP ,
     1          SA, SC, SF, SL, SN, SR,
     1          RA, RC, RF, RL, RN, RR,
     1          rayp, geom, tstar, iiso, .false.)
        


            degrad=0.0174532925
            tol = 1.0e-7
            dp = degrad*dip
            st = degrad*stk
            rk = degrad*rake
            sind=sin(dp)
            cosd=cos(dp)
            sinr=sin(rk)
            cosr=cos(rk)
            sins=sin(st)
            coss=cos(st)
            sin2d=sin(2.*dp)
            cos2d=cos(2.*dp)
            sin2s=sin(2.*st)
            cos2s=cos(2.*st)
            s1 = cosr*coss + sinr*cosd*sins
            s2 = cosr*sins - sinr*cosd*coss
            s3 = - sinr*sind
            n1 = -sins*sind
            n2 =  coss*sind
            n3 = -cosd
            xmt(1,1)=SA*s1*n1+(SA-2*SN)*s2*n2 + SF*s3*n3
            xmt(2,2)=SA*s2*n2+(SA-2*SN)*s1*n1 + SF*s3*n3
            xmt(3,3)=SF*s1*n1 +SF*s2*n2 + SC*s3*n3
            xmt(1,2)=SN*(s1*n2+s2*n1)
            xmt(1,3)=SL*(s1*n3+s3*n1)
            xmt(2,3)=SL*(s2*n3+s3*n2)
            xmt(2,1) = xmt(1,2)
            xmt(3,1) = xmt(1,3)
            xmt(3,2) = xmt(2,3)
c-----
c           get the norm of the matrix and then adjust to
c           the desired moment
c----- 
            mijnorm = 0.0 
            do i=1,3
               do j=1,3
                    mijnorm = mijnorm + xmt(i,j)**2
              enddo
            enddo
            mijnorm = sqrt(mijnorm/2.0)
c-----
c           Silver and jordan
c-----
            do i=1,3
               do j=1,3
                    xmt(i,j) = xmt(i,j) * xmom / mijnorm
               enddo
            enddo
c-----
c           clean up small values
c-----
      
            xmax=-1.0e+37
            do  i=1,3
                do  j=1,3
                    if(abs(xmt(i,j)).gt.xmax)xmax = abs(xmt(i,j))
                enddo
            enddo
            thresh = tol * xmax
            do  i=1,3
                do  j=1,3
                    if(abs(xmt(i,j)).lt.thresh) xmt(i,j) = 0.0
                enddo
            enddo

c-----
c           write out the information
c-----
           WRITE(LOT,*)' Strike           :',stk
           WRITE(LOT,*)' Dip              :',dip
           WRITE(LOT,*)' Rake             :',rake
           WRITE(LOT,*)' Moment           :',xmom
           WRITE(LOT,*)' xmt              :',xmt(1,1),xmt(1,2),xmt(1,3)
           WRITE(LOT,*)'                  :',xmt(2,1),xmt(2,2),xmt(2,3)
           WRITE(LOT,*)'                  :',xmt(3,1),xmt(3,2),xmt(3,3)
           WRITE(LOT,*)' Trace            :',xmt(1,1)+xmt(2,2)+xmt(3,3)

      
        return
        end

        subroutine gcmdln(ntau,ipt,alp,idva,
     1      rfile,
     2      dozero,
     1      dostrain,dostress,dorotate,dogreen,
     1      az,baz,stk,dip,rake,xmom,xmt,fx,fy,fz,
     1      dotest1,dostep,outfmt,dosdr)
c-----
c       parse command line arguments
c       requires subroutine mgtarg() and function mnmarg()
c-----
c
c       ntau    I*4 - pulse duration factor for parabolic and triangular
c               pulses
c       ipt I*4 - pulse type
c               0 - triangular
c               1 - parabolic
c               2 - Ohnaka
c               3 - Dirac delta function
c               4 = user pulse in file rfile
c       alp R*4 - Ohnaka pulse shape factor, fc= alp / 2 pi
c       idva    I*4 - time history type
c               0 - displacement
c               1 - velocity
c               2 - acceleration
c       rfile   C*80- name of user provided pulse 
c               note dt must be the same as specified in dfile
c       dozero  L   .true. for triangular and parabolic pulses, center
c               at zero lag
c       dostrain  - L .true. output strain
c       dostress  - L .true. output stress
c       dorotate  - L .true. output rotation
c       dogreen   - L .true. output Green as in f96tosac
c       az        - R  aximuth
c       baz       - R back azimuth not used
c       stk       - R fault strike
c       dip       - R fault dip   
c       rake      - R fault rake  
c       xmom      - R moment in dyne-cm
c       xmt       - moment tensor (3,3)
c       fx, fy, fz - for specification north, east, down (dynes)
c       dotest    - L
c       dostep    - L   .true.  source time function is step-like
c                         whose derivative is the triangular,parabolic
c                         Ohnaka or impulse
c                   .false. source time function is impulsive.
c                         The impulse is approximated by the triangular,
c                         parabolic, Ohnaka or impulse function
c                         When -STEP -D are given this is the Green's 
c                         function.
c       outfmt    -I specification of output file name
c                    -FMT 1      DDDDDd_HHHh_ZZZz.cmp
c                                e.g. 005001_1234_0045.Uz
c                    -FMT 2      DDDDDddd_HHHhhh_ZZZzzz.cmp
c                                e.g. 00500123_123456_004578.Erf
c                    -FMT 3      DDDDDdHHHh.grn(default)
c                                e.g. 0050010041.ZVF
c                    -FMT 4      DDDDdHHHh.grn
c                                e.g. 050010045.Srz
c                    -FMT 5      DDDdddHhhh.grn
c                                e.g. 5001234578.Err
c-----
        character*(*) rfile
        integer*4 ntau, ipt, idva
        real*4 alp
        logical  dozero, dotest1,dostep
        logical dostrain, dostress,dorotate,dogreen,dosdr
        real az, baz,stk,dip,rake, xmt(3,3), xmom,forcex, forcey,forcez
        integer outfmt

        integer LER, LIN, LOT
        parameter (LER=6, LIN=5, LOT=6)

        integer*4 mnmarg
        character*80 name
        integer i,j
        logical lmij, lex, lsdr, lforce
        integer lgstr
c-----
c       initialization
c-----
        ntau = -1
        ipt = -1
        alp = -1.0
        idva = 1
        rfile = ' '
        dozero = .false.
        xmom = 1.0
        stk    = 0
        dip = 0
        rake = 0
        xmom = 1.0
        lsdr   = .false.
        lmij   = .false.
        lex    = .false.
        lforce = .false.
        dostrain =.false.
        dostress = .false.
        dorotate = .false.
        dogreen  = .false.
        dotest1 = .false.
        dostep = .true.
        dosdr = .false. 
        outfmt = 1
        isds = -1
        fx = 0.0
        fy = 0.0
        fz = 0.0
        az = 0
        baz = 0
        do i=1,3
           do j=1,3
              xmt(i,j) = 0.0
           enddo
        enddo
c-----
c     process the command line arguments
c----
        nmarg=mnmarg()
        i = 0
   11   i = i + 1
            if(i.gt.nmarg)go to 13
            call mgtarg(i,name)
            if(name(1:2).eq.'-l')then
                i=i+1
                call mgtarg(i,name)
                read(name,'(bn,i10)')ntau
            else if(name(1:2).eq.'-a')then
                i=i+1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')alp
            else if(name(1:2).eq.'-t')then
                ipt = 0
            else if(name(1:2).eq.'-p')then
                ipt = 1
            else if(name(1:2).eq.'-o')then
                ipt = 2
            else if(name(1:2).eq.'-i')then
                ipt = 3
            else if(name(1:2).eq.'-Z'.and.
     1          name(1:3).ne.'-ZZ')then
                dozero = .true.
            else if(name(1:2).eq.'-F' .and. name(1:3).ne.'-FU'
     1         .and. name(1:3).ne.'-FX'
     1         .and. name(1:3).ne.'-FY'
     1         .and. name(1:3).ne.'-FZ'
     1         .and. name(1:3).ne.'-FM'
     1          )then
                ipt = 4
                        i=i+1
                        call mgtarg(i,rfile)
            else if(name(1:4).eq.'-FMT')then
                i=i+1
                call mgtarg(i,name)
                read(name,'(bn,i10.0)')outfmt
                if(outfmt.lt.1 .or. outfmt.gt.5)then
                   outfmt = 4
                endif
            else if(name(1:2).eq.'-D' .and.
     1          name(1:3).ne.'-DI')then
                idva = 0
            else if(name(1:2).eq.'-V')then
                idva = 1
            else if(name(1:2).eq.'-A' .and. name(1:4).ne.'-ALL'
     1          .and. name(1:3).ne.'-AZ')then
                idva = 2
            else if(name(1:5).eq.'-STEP')then
                dostep = .true.
            else if(name(1:4).eq.'-IMP')then
                dostep = .false.
            else if(name(1:4).eq.'-DIP')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,dip)
                lsdr = .true.
                dosdr = .true.
                isds = 0
            else if(name(1:4).eq.'-STK')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,stk   )
                lsdr = .true.
                dosdr = .true.
                isds = 0
            else if(name(1:5).eq.'-RAKE')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,rake)
                lsdr = .true.
                dosdr = .true.
                isds = 0
            else if(name(1:3).eq.'-M0' .or. name(1:3).eq.'-MO')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,xmom)
            else if(name(1:3).eq.'-MW' .or. name(1:3).eq.'-Mw')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,xmom)
                xmom = 10.**(1.5*xmom + 16.10)
            else if(name(1:3).eq.'-AZ')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,az)
            else if(name(1:4).eq.'-BAZ')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,baz)
            else if(name(1:3).eq.'-EX')then
                isds = 2
            else if(name(1:3).eq.'-xx' .or. name(1:3).eq.'-XX')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,xmt(1,1))
                lmij = .true.
                isds = 1
            else if(name(1:3).eq.'-yy' .or. name(1:3).eq.'-YY')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,xmt(2,2))
                lmij = .true.
                isds = 1
            else if(name(1:3).eq.'-zz' .or. name(1:3).eq.'-ZZ')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,xmt(3,3))
                lmij = .true.
                isds = 1
            else if(name(1:3).eq.'-xy' .or. name(1:3).eq.'-XY')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,xmt(1,2))
                xmt(2,1) = xmt(1,2)
                lmij = .true.
                isds = 1
            else if(name(1:3).eq.'-xz' .or. name(1:3).eq.'-XZ')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,xmt(1,3))
                xmt(3,1) = xmt(1,3)
                lmij = .true.
                isds = 1
            else if(name(1:3).eq.'-yz' .or. name(1:3).eq.'-YZ')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,xmt(2,3))
                xmt(3,2) = xmt(2,3)
                lmij = .true.
                isds = 1
            else if(name(1:3).eq.'-fx' .or. name(1:3).eq.'-FX')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,fx)
                lforce = .true.
            else if(name(1:3).eq.'-fy' .or. name(1:3).eq.'-FY')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,fy)
                lforce = .true.
            else if(name(1:3).eq.'-fz' .or. name(1:3).eq.'-FZ')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,fz)
                lforce = .true.
            else if(name(1:7).eq.'-STRAIN')then
                dostrain = .true.
            else if(name(1:7).eq.'-STRESS')then
                dostress =.true.
            else if(name(1:4).eq.'-ROT')then
                dorotate = .true.
            else if(name(1:4).eq.'-GRN')then
                dogreen  = .true.
            else if(name(1:6).eq.'-TEST1')then
                dotest1 = .true.
            else if(name(1:2).eq.'-?')then
                call usage('Help')
            else if(name(1:2).eq.'-h')then
                call usage('Help')
            endif
        go to 11
   13   continue
c-----
c     do some error checking, e.g., we cannot generate velocity
c     for triangular pulse
c-----
        if(ipt.ge.0 .and. ipt .le.1 .and . ntau.le.0)ntau = 1
        if(ipt.eq.2 .and. alp .le.0.0)alp = 1.0
        if(ipt.eq.2 .and. alp.lt.0.0)
     1      call usage('No alpha for Ohnaka pulse')
        if(ipt.lt.0)
     1      call usage('No pulse shape defined')
c-----
c       Priorties for the source specification
c       Mij > Strike, Dip, Slip >  Explosion > Force
c       This means that if Mij  is specified, the others are ignored
c       If Mij not given, but strike,dip,slip are this overrides the -EX
c
c       not if anything if Mij turn off the force
c-----
        if(isds.eq.2)then
           xmt(1,1) = xmom
           xmt(2,2) = xmom
           xmt(3,3) = xmom
                WRITE(LER,*)' Mij not given, Explosion, no force'
                fx = 0.0
                fy = 0.0
                fz = 0.0
        endif
        if(.not.lmij)then
                if(lforce)then
                   WRITE(6,*)'Force only'
                endif
        else
c-----
c          turn off parameters that are not related to MT source
c-----     
           WRITE(LER,*)' Mij given no force'
           fx = 0.0
           fy = 0.0
           fz = 0.0
        endif
        WRITE(6,*)'STK,DIP,RAKE:',STK,DIP,RAKE
        WRITE(6,*)'XMOM,AZ:',XMOM,AZ
        return
        end

        subroutine usage(str)
        parameter (LER=0, LIN=5, LOT=6)
        character str*(*)
        write(LER,*)'hpulse96strain:',str
        write(LER,'(a/a/a/a/a/a/a/a)')'USAGE: ',
     1  'hpulse96strain -d Distance_File [ -t -o -p -i ] [-a alpha]',
     2  '    -l L [ -D|-V |A]  [-F rfile ] [ -m mult] [-STEP|-IMP]',
     3  '    [-STRESS  -STRAIN -ROTATE -GRN] [-FUND] [-HIGH] [-Z] ',
     4  '    [-LAT] [-2] [ -M mode ] [-LOCK] -FMT ifmt ',
     1  '    [-M0 moment ] [-MW mw] [-STK stk -DIP dip -RAKE rake]',
     2  '    [-FX fx -FY fy -FZ fz] ',
     3  '    [-XX Mxx ... -ZZ Mzz] [-?] [-h]'
        write(LER,*)'TIME FUNCTION SPECIFICATION'
        write(LER,'(a/a/a/a/a/a/a/a/a)')
     1  '  -t           Triangular pulse of base 2 L dt',
     1  '  -p           Parabolic Pulse of base  4 L dt',
     1  '                  -p -l 1 recommended',
     1  '  -l L         (default 1 )duration control parameter',
     1  '  -o           Ohnaka pulse with parameter alpha',
     1  '  -i           Dirac Delta function',
     1  '  -a alpha     Shape parameter for Ohnaka pulse',
     1  '  -F rfile     User supplied pulse',
     1  '  -m mult      Multiplier (default 1.0)'
        write(LER,'(a/a/a/a,a)')
     1  '  -STEP        (default)',
     1  '  -IMP         ',
     1  '              By default the source time function is ',
     1  '              steplike. -IMP forces impulse like. -D -IMP ',
     1  'is Green s function'
        write(LER,*)'OUTPUT FILE NAME'
        write(LER,'(a/a/a/a/a/a/a/a/a)')
     1  '  The format for the name of the binary output attempts to',
     1  '  give information on epicentral distance (km), ',
     1  '  source depth (km), and receiver depth(km). The options are',
     1  '  -FMT 1      DDDDDd_HHHh_ZZZz.cmp',
     1  '              e.g. 005001_1234_0045.Uz',
     1  '  -FMT 2      DDDDDddd_HHHhhh_ZZZzzz.cmp',
     1  '              e.g. 00500123_123456_004578.Erf',
     1  '  -FMT 3      DDDDDdHHHh.grn(default)',
     1  '              e.g. 0050010041.ZVF'
        write(LER,'(a/a/a/a/a/a/a/a/a)')
     1  '  -FMT 4      DDDDdHHHh.grn',
     1  '              e.g. 050010045.Srz',
     1  '  -FMT 5      DDDdddHhhh.grn',
     1  '              e.g. 5001234578.Err',
     1  '  where D is for epicentral distance, H source depth, and',
     1  '  Z receiver depth. The lower case indicates the digits ',
     1  '  to the right of the decimal place.  The examples above',
     1  '  are for an epicentral distance is 500.123 km, source',
     1  '  depth 123.456 km and receiver depth 4.578 km.'

        WRITE(LER,*)'OUTPUT TIMESERIES FOR SOURCE as Ur, Ut, Uz',
     1  ' components with strain, stress optional'
        write(LER,'(a/a/a)')
     1  '  -D           Output is ground displacement        (m)',
     1  '  -V           Output is ground velocity (default) (m/s)',
     1  '  -A           Output is ground acceleration       (m/s^2)'
        write(LER,'(a/a/a/a/a/a)')
     1  '  -STRESS (default .false. ) output stress for mechanism',
     1  '    units are Pa, with suffix Srr, Srf, Srz, Stt, Sfz, Szz',
     1  '  -STRAIN (default .false. ) output strain for mechanism',
     1  '       with suffix, Err, Erf, Erz, Eff, Efz, Ezz',
     1  '    NOTE the Ur, Ut, Uz components are created with -STRESS',
     1  '       or -STRAIN flags. The Uz is positive down.',
     1  '  -ROTATE (default .false. ) output rotation for mechanism',
     1  '       with suffix, Wfz, Wrz, Wrf'
        write(LER,'(a/a/a/a/a)')
     1  '  -GRN    (default false) Output Green;s functions',
     1  '    hpulse96strain -STEP -V -p -l 1 -GRN -FMT 4  is same as',
     1  '     hpulse96 -V -p -l 1 | f96tosac -G . For KM,KM/S,GM/CM^3',
     1  '     model, output will be CM/S for moment of 1.0e+20 dyne-cm',
     1  '     of force of 1.0e+15 dyne',
     1  '    NOTE the Z component of ZDD ... ZHF is positive up'
        write(LER,'(a/a/a/a/a)')
        write(LER,'(a/a/a/a/a)')
     1  '  -TEST1  (default .false.) output CPS Green functions ,e.g.,',
     1  '       ZDS RDS ... RHF THF for use with moment tensor codes',
     1  '       and gsac MT command. This is equivalent to ',
     1  '       hpulse96 -V -p -l 1 | f96tosac -G if -FMT 4 is used ',
     1  '       with hpulse96strain'

        write(LER,*)'COMPUTATIONS'
        write(LER,'(a)')
     1  '  -Z         (default false) zero phase '
c------
c      strain and stress and rotation
c------
        write(LER,*)'SOURCE MECHANISM SPECIFICATION'
        write(LER,'(a/a/a/a/a/a/a,a)')
     1  '  -DIP dip               dip of fault plane',
     1  '  -STK Strike            strike of fault plane',
     1  '  -RAKE Rake              slip angle on fault plane',
     1  '  -M0 Moment (def=1.0) Seismic moment in units of dyne-cm',
     1  '  -MW mw            Moment Magnitude  ',
     1  '            moment (dyne-cm) from log10 Mom = 16.10 + 1.5 Mw',
     1  '            For strike,dip,rake source mw or Moment must',
     1  ' be specified'
        write(LER,'(a/a/a/a,a/a/a/a,a)')
     1  '  -EX                  Explosion',
     1  '  -AZ Az                Source to Station Azimuth',
     1  '  -BAZ Baz               Station to Source azimuth',
     1  '  -fx FX -fy Fy -fZ fz  Point force amplitudes ',
     2  ' (N,E,down) in  dynes',
     1  '  -XX Mxx -YY Myy -ZZ Mzz  Moment tensor elements in units of',
     1  '  -XY Mxy -XZ Mxz -YZ Myz    dyne-cm',
     1  ' The moment tensor coordinates are typically X = north',
     2  ' Y = east and Z = down'
        write(LER,'(a/a)')
     1  ' If by accident more than one source specification is used,',
     1  ' the hierarchy is Mij > Strike,dip,rake > Explosion > Force'
        write(LER,*)
     1  '--------------------------------------------------------------'
        write(LER,*)
     1  'NOTE: The output units are related tot he model specification.'
        write(LER,*)
     1  'To have the desired units the model must be in KM, KM/S',
     1  '  and GM/CM^3'
        write(LER,*)
     1  '--------------------------------------------------------------'

        write(LER,*)
     1  ' -?           Write this help message'
        write(LER,*)
     1  ' -h           Write this help message'
        stop
        end

        subroutine zfour(zarr,nn,isign,dt,df) 
c-----
c     THE input is a complex array
c     which has numbers stored in memory as
c     R1, I1, R2, I2, ..., Rnn, Inn
c     where nn must be a power of 2 R and I are the real and imaginary
c     parts of the complex number
c
c     For isign -1 this is a complex time series
c     For isign +1 this is a complex frequency series with
c        index 1 (in fortran corresponding to f=0
c              2                              f=df
c            nn/2 + 1                         f = 1/2dt = Nyquist
c            nn - 1                           f = -2df
c             nn                              f = -df

c-----
c     the cooley-tookey fast fourier transform in usasi basic fortran
c     transform(j) = sum(datc(i)*w**((i-1)(j-1)), where i and j run
c     from 1 to nn and w = exp(isign*2*pi*sqrt(-1)/nn).  datc is a one-
c     dimensional complex array (i.e., the real and imaginary parts of
c     datc are located immediately adjacent in storage, such as fortran
c     places them) whose length nn is a power of two.  isign
c     is +1 or -1, giving the sign of the transform.  transform values
c     are returned in array datc, replacing the input datc.  the time is
c     proportional to n*log2(n), rather than the usual n**2
c     rms resolution error being bounded by 6*sqrt(i)*log2(nn)*2**(-b),
c     b is the number of bits in the floating point fraction.
c
c     the program computes df from dt, dt from df and checks to see
c     if they are consistent. In addition, the transforms are multiplied
c     by dt or df to make the results dimensionally correct
c
c     This is a slightly modified version of the original Brenner routine
c     The changes were to add physical dimensions to the transform
c     and to make it all complex
c-----
        complex zarr(*) 
        integer nn, isign
        real dt, df

        complex ztemp
c-----
c       ensure that the dt and df are defined and
c       consistent
c-----
        if(dt.eq.0.0) dt = 1./(nn*df) 
        if(df.eq.0.0) df = 1./(nn*dt) 
        if(dt.ne.(nn*df)) df = 1./(nn*dt) 
c-----
c       now begin the transform
c-----
        jj = 1
        do 5 ii=1,nn 
        if(ii .lt. jj) then
              ztemp = zarr(jj)
              zarr(jj) = zarr(ii)
              zarr(ii) = ztemp
        endif
        mm = nn/2
    3   continue
        if(jj.le.mm) then
            go to 55
        else 
              jj = jj - mm
              mm = mm /2
              if(mm.lt.1)then
                  go to 55
              else
                  go to 3
              endif
        endif
   55   continue
        jj = jj + mm
    5   continue
        mmax = 1 
c-----
    6 continue
        if(mmax .lt. nn)then
            go to 7
        else if(mmax .ge. nn)then
            go to 10
        endif
    7   continue
        istep= 2*mmax 
        theta = 6.283185307/(isign*2.0*mmax) 
        sinth=sin(theta/2.) 
        wstpr=-2.*sinth*sinth 
        wstpi=sin(theta) 
        wr=1.0 
        wi=0.0 
        do 9 mm=1,mmax
              do 8 ii=mm,nn,istep
                    jj=ii+mmax
                    ztemp=cmplx(wr,wi)*zarr(jj)
                    zarr(jj) = zarr(ii) - ztemp
                    zarr(ii) = zarr(ii) + ztemp
    8         continue
c-----
c       use trig relations to compute the next sin/cos
c       without actually calling sin() or cos()
c-----
              tempr = wr 
              wr = wr*wstpr-wi*wstpi + wr 
              wi = wi*wstpr+tempr*wstpi + wi 
    9   continue
        mmax = istep 
        go to 6 
c-----
c       transform is done
c-----
   10   continue 
c-----
c     give the arrays the proper physical dimensions
c-----
        if(isign.lt.0)then
c-----
c             time to frequency domain
c-----
              do  ii = 1,nn
                    zarr(ii) = zarr(ii) * dt
              enddo
        else
c-----
c             frequency to time domain
c-----
              do ii = 1,nn
                    zarr(ii) = zarr(ii) * df
              enddo
        endif
        return
        end

        subroutine outstrain(ntau,idva,dist,dt,nptin,
     1      tfirst,depthr,depths,az,
     2      mname,dostep,alpha,dosdr,TP,TSV,TSH)
c-----
c       write the Green s functions for this distance
c-----
c       ntau    I*4 - pulse duration parameter
c       idva    I*4 - displacmeent, velocity or acceleration
c       dist    R*4 - distance
c       dt  R*4 - sample rate
c       nptin   I*4 - number of points in time series
c       tfirst  R*4 - time of first sample
c       depthr  R*4 - receiver depth
c       depths  R*4 - source depth
c       az      R   - source to station azimuth
c       mname   C*(*)   - name of model file
c-----
        real stk, dip, rake,xmom
        logical dosdr

        character stcomp*8, stname*8
        real*4 cmpaz, cmpinc, cmpdt
        integer*4 npts
        real*4 ssec

        integer NL
        parameter (NL=200)
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL),
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax

        common/receiver/l2mu,mu
        real l2mu, mu
        common/grnfmt/outfmt
        integer outfmt

c-----
c       internal variables
c-----
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
        integer NGRN
        parameter (NGRN=45)
        integer*4 ntau, idva,   nptin
        real*4 dt, tfirst, depths, depthr

        integer kerr, system

        integer NSAMP, NFREQ
        parameter (NSAMP=16384,NFREQ=8193)
        common/dd/datc
        complex datc(NSAMP)
        common/xxx/x
        real*4 x(NSAMP)
        complex xx(NGRN)
        complex z
         integer iszrt(45)
        character*8 ostnm(2)
        character*8 ocmpnm(10)

        character mname*80
        logical dostep
        common/zzy/datas
        complex datas(NFREQ)

        real svel, pvel, vsa, vsb, vsr, geom, rayp, tstar, den

        real caz, saz, c2az, s2az

        complex uz,ur,ut,duzdz,durdz,dufdz,
     1      duzdr,durdr,dufdr,duzdf,durdf,dufdf

        complex err,efr,erz,eff,efz,ezz,del
        complex xout(10)

        character ofile*80
        integer ls
        

        data ostnm/'MT      ', 'F0      '/
        data ocmpnm/'Uz      ','Ur      ','Ut      ',
     1       'Err     ', 'Erf     ',
     1       'Erz     ', 'Eff     ', 'Efz     ', 'Ezz     ',
     1       'Del     '/


        tau = ntau * dt
c-----
c       get first arrival time - NOTE for implementation of
c-----
            tTP  = -12345.
            tTSV = -12345.
            tTSH = -12345.
            CALL FRSTAR(dist,depths,depthr,
     1          MNAME,1,tTP ,
     1          SSA, SSC, SSF, SSL, SSN, SSR,
     1          RRA, RRC, RRF, RRL, RRN, RRR,
     1          rayp, geom, tstar, iiso, .false.)
            CALL FRSTAR(dist,depths,depthr,
     1          MNAME,2,tTSV,
     1          SSA, SSC, SSF, SSL, SSN, SSR,
     1          RRA, RRC, RRF, RRL, RRN, RRR,
     1          rayp, geom, tstar, iiso, .false.)
            CALL FRSTAR(dist,depths,depthr,
     1          MNAME,3,tTSH,
     1          SSA, SSC, SSF, SSL, SSN, SSR,
     1          RRA, RRC, RRF, RRL, RRN, RRR,
     1          rayp, geom, tstar, iiso, .false.)
c-----
c       Now that all spectra for a given distance are in
c       make time history and output in desired order
c-----
        np2 = nptin/2 + 1
        df = 1.0/(nptin*dt)
c-----
c       perform inverse FFT
c       If the number of data points is too large, then use supfft
c-----
c       here we will form the synthetics
c       there will be a total of 9, e.g.,
c       Uz  Ur Ut Err Ert Erz Ett Etz Ezz
c-----
c-----
c       jout   1  Uz
c              2  Ur
c              3  Ut
c              4  Err
c              5  Ert
c              6  Erz
c              7  Ett
c              8  Etz
c              9  Ezz
c----
c       This is a two pass process
c       Starting with  ur uz tz durdz and other derivatives on LUN=10
c       compute ur uz ut and strains an then output them in LUN=9
c       
c-----
        open(9,file='file9.tmp',access='sequential', 
     1   form='unformatted', status='unknown')
        rewind 9
        rewind 10
        do  k=1,np2
              read(10) (xx(i),i=1,NGRN)
              call domom(xx,uz,ur,ut,duzdz,durdz,dufdz,
     1            duzdr,durdr,dufdr,duzdf,durdf,dufdf)
c-----
c     dist is in km 
c     the following is the strain formula
c     be careful of units especially the eff
c     Some derivatives are dU/dr or dU/dz obviously
c     we must convert from km to m
c-----
              r = dist
              err = durdr
              eff = (dufdf + ur)/r
              ezz = duzdz
              efr = 0.5*(durdf/r  + dufdr - ut/r)
              erz = 0.5*(duzdr + durdz)
              efz = 0.5*(dufdz + duzdf/r)
              del = err + eff + ezz
C         ezz = duzdf
C         err = durdf
C         eff = dufdf
C         efr = ur
c-----
c     CONVERT FROM CM TO M
c-----
             xout(1) = uz/1.0e+02
             xout(2) = ur/1.0e+02
             xout(3) = ut/1.0e+02
c-----
c     CONVERT FROM KM TO M for distance and CM to M for displacement
c-----
             xout(4) = err/1.0e+05
             xout(5) = efr/1.0e+05
             xout(6) = erz/1.0e+05
             xout(7) = eff/1.0e+05
             xout(8) = efz/1.0e+05
             xout(9) = ezz/1.0e+05
             xout(10) = del/1.0e+05
              write(9)(xout(i),i=1,10)
         enddo
c-----
c     now make 9 passes to make the synthetics
c     for ur uz ut err ert erz ett etz ezz
c-----
         jkup = 10
         rewind 9
         do 1300 jk=1,jkup
                rewind 9
                do k=1,np2
                    read(9) (xout(i),i=1,10)
                    datc(k)=xout(jk)
                    if(k.gt.1)then
                        datc(nptin+2-k)=conjg(datc(k))
                    endif
                enddo
c-----
c       perform inverse FFT
c-----
c       also take the opportunity to perform the correct
c       spherical spreading for an earth with radius of 6371.0 km
c-----
                    call zfour(datc,nptin,+1,dt,df)
                    fac = exp(alpha*tfirst)
                    dfac = exp(alpha*dt)
c-----
c     undamp 
c-----
                    do 1303 k=1,nptin
                        x(k) = real(datc(k))*fac
                        fac = fac*dfac
 1303               continue
                    npts = nptin
c-----
c           convert velocity to displacement or acceleration
c           if required. Since we are using surface waves, the
c           DC component must be zero!
c
c           note the special case for the pressure field
c-----
                IF(JK.LE.3)THEN
                    if(dostep)then
                       if(idva.eq.0)then
                           call integ(x,npts,dt,.false.)
                       else if(idva.eq.2)then
                           call deriv(x,npts,dt)
                       endif
                    else
                       if(idva.eq.1)then
                           call deriv(x,npts,dt)
                       else if(idva.eq.2)then
                           call deriv(x,npts,dt)
                           call deriv(x,npts,dt)
                       endif
                    endif
                else
                    if(dostep)then
                        call integ(x,npts,dt,.false.)
                    endif
                endif
c-----
c     component orientations SEED convention. Later 90 is added to 
c     make Sac convention. Difference is that Sac measures from upward
c     vertical and SEED from downward dip 
c            Component SEED     Sac
c     CMPAZ     Z        0       0
c     CMPINC           -90       0
c
c     CMPAZ     E       90      90
c     CMPINC             0      90
c
c     CMPAZ     N        0       0
c     CMPINC             0      90
c
c     These are Sac convention
c     Here the Uz is positive downward, e.g., SAC CMPINC 180
c-----
c           vertical
c-----
                if(jk.eq.1)then
                    cmpinc = 180.0 
                    cmpaz  =   0.0
c-----
c           radial
c-----
                else if(jk.eq.2)then
                    cmpinc =  90.0 
                    cmpaz  =   0.0 + az
c-----
c           transverse
c-----
                else if(jk.eq.3)then
                    cmpinc =  90.0 
                    cmpaz  =  amod(90.0 + az, 360.0)
                else
                    cmpinc = -12345.
                    cmpaz  = -12345.
                endif
                stcomp = 'SYN     '
    
                cmpdt  = dt
                ksyear = 1970
                ksdoy  = 1
                ksmon  = 1
                ksday  = j
                kshour = 0
                ksmin  = 0
                ssec =  0.0

                o = 0
                distdg = -12345.
c-----
c
c-----
                stcomp = ocmpnm(jk)(1:3)
c-----
c     create the output file name
c-----
                 if(outfmt.eq.1)then
c                   DDDDDd_HHHh_ZZZz.cmp
                    ievdep = 10.*depths
                    istel  = 10.*depthr
                    idist  = 10.0*dist
                    write(ofile,11)idist,ievdep,istel,stcomp
   11               format(i6.6,'_',i4.4,'_',i4.4,'.',a3)
                 else if (outfmt.eq.2)then
c                   DDDDDddd_HHHhhh_ZZZzzz.cmp
                    ievdep = 1000.*depths
                    istel  = 1000.*depthr
                    idist  = 1000.*dist
                    write(ofile,12)idist,ievdep,istel,stcomp
   12               format(i7.7,'_',i6.6,'_',i6.6,'.',a3)
                 else if (outfmt.eq.3)then
c                   DDDDDdHHHh.grn
                    ievdep = 10.*depths
                    idist  = 10.0*dist
                    write(ofile,13)idist,ievdep,stcomp
   13               format(i6.6,i4.4,'.',a3)
                 else if (outfmt.eq.4)then
c                   DDDDdHHHh.grn
                    ievdep = 10.*depths
                    idist  = 10.0*dist
                    write(ofile,14)idist,ievdep,stcomp
   14               format(i5.5,i4.4,'.',a3)
                 else if (outfmt.eq.5)then
c                   DDDdddHhhh.grn
                    ievdep = 1000.*depths
                    idist  = 1000.*dist
                    write(ofile,15)idist,ievdep,stcomp
   15               format(i6.6,i4.4,'.',a3)
                 endif
                 ls = lgstr(ofile)
                call scmxmn(x,npts,depmax,depmin,depmen,indmax,indmin)
                call newhdr()
                call setfhv('DELTA   ',cmpdt,nerr)
                    bb = tfirst
                    ee = bb + (npts -1 ) * cmpdt
                    call setfhv('A       ',tp    ,nerr)
                    call setkhv('KA      ','P       ',nerr)
                    call setfhv('T0      ',tsv   ,nerr)
                    call setfhv('T1      ',tsh   ,nerr)
                    if(iiso.eq.0)then
                       call setkhv('KT0     ','S       ',nerr)
                       call setkhv('KT1     ','S       ',nerr)
                    else
                       call setkhv('KT0     ','SV      ',nerr)
                       call setkhv('KT1     ','SH      ',nerr)
                    endif
                call setfhv('B       ',bb    ,nerr)
                call setfhv('E       ',ee    ,nerr)
                call setfhv('O       ',o     ,nerr)
                call setkhv('KO      ','O       ',nerr)
                call setfhv('STLA    ',stlat ,nerr)
                call setfhv('STLO    ',stlon ,nerr)
                call setfhv('TIMMAX  ',bb+ indmax*cmpdt,nerr)
                call setfhv('TIMMIN  ',bb+ indmin*cmpdt,nerr)
                call setfhv('DEPMIN  ',depmin,nerr)
                call setfhv('DEPMAX  ',depmax,nerr)
c-----
c               our elevation is KM, SAC is set to KILOMETERS
c-----
                call setfhv('STEL    ',depthr,nerr)
                call setfhv('EVLA    ',evlat ,nerr)
                call setfhv('EVLO    ',evlon ,nerr)
c-----
c               our depth is KM, SAC is set to KILOMETERS
c-----
                call setfhv('EVDP    ',depths ,nerr)
                call setfhv('DIST    ',dist,nerr)
                evstaz =   az
                stevaz = amod(180.0+az, 360.0)
                call setfhv('AZ      ',evstaz,nerr)
                call setfhv('BAZ     ',stevaz,nerr)
                call setfhv('GCARC   ',distdg,nerr)
                call setfhv('CMPAZ   ',cmpaz ,nerr)
                call setfhv('CMPINC  ',cmpinc,nerr)
                call setfhv('DEPMEN  ',depmen,nerr)
                call setfhv('STLA    ',-12345.,nerr)
                call setfhv('STLO    ',-12345.,nerr)
                call setfhv('EVLA    ',-12345.,nerr)
                call setfhv('EVLO    ',-12345.,nerr)
c-----
c       set integer header value
c-----
                call setnhv('NZYEAR  ',ksyear,nerr)
                call setnhv('NZJDAY  ',ksdoy ,nerr)
                call setnhv('NZHOUR  ',kshour,nerr)
                call setnhv('NZMIN   ',ksmin ,nerr)
                isec = ssec
                smsec = (ssec - isec) * 1000.0
                call setnhv('NZSEC   ',isec  ,nerr)
                jsmsec = smsec
                call setnhv('NZMSEC  ',jsmsec,nerr)
                call setnhv('NPTS    ',npts, nerr)
c-----
c       set logical header value
c-----
                call setlhv('LEVEN   ',.true.,nerr)
                call setlhv('LPSPOL  ',.true.,nerr)
                call setlhv('LOVROK  ',.true.,nerr)
                call setlhv('LCALDA  ',.false.,nerr)
c-----
c       set character header value
c-----
                call setkhv('KSTNM    ','SYN     ' ,nerr)
                if(iobsyn.eq.2)then
                    call setkhv('KEVNM   ','SYNTHETI',nerr)
                    call setkhv('KEVNMC  ','C       ',nerr)
                else
                    call setkhv('KEVNM   ','OBSERVED',nerr)
                    call setkhv('KEVNMC  ','        ',nerr)
                endif
                call setkhv('KCMPNM   ',stcomp ,nerr)
c-----
c       set enumerated header value
c-----
                call setihv('IFTYPE   ','ITIME   ',nerr)
                call setihv('IZTYPE   ','IB      ',nerr)
                    call bwsac(12,npts,ofile(1:ls),x)
 1300   continue
        close(9,status='delete')
        return
        end

        subroutine outtest1(ntau,idva,dist,dt,nptin,
     1      tfirst,depthr,depths,
     2      mname,dostep,alpha,TP,TSV,TSH)
c-----
c       write the Green s functions for this distance
c-----
c       ntau    I*4 - pulse duration parameter
c       idva    I*4 - displacmeent, velocity or acceleration
c       dist    R*4 - distance
c       dt  R*4 - sample rate
c       nptin   I*4 - number of points in time series
c       tfirst  R*4 - time of first sample
c       depthr  R*4 - receiver depth
c       depths  R*4 - source depth
c       mname   C*(*)   - name of model file
c-----
        character stcomp*8, stname*8
        real*4 cmpaz, cmpinc, cmpdt
        integer*4 npts
        real*4 ssec

        integer NL
        parameter (NL=200)
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL),
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax

        common/receiver/l2mu,mu
        real l2mu, mu

c-----
c       internal variables
c-----
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
        integer NGRN
        parameter (NGRN=45)
        integer*4 ntau, idva,   nptin
        real*4 dt, tfirst, depths, depthr

        common/grnfmt/outfmt
        integer outfmt

        integer kerr, system

        integer NSAMP, NFREQ
        parameter (NSAMP=16384,NFREQ=8193)
        common/dd/datc
        complex datc(NSAMP)
        common/xxx/x
        real*4 x(NSAMP)
        complex xx(NGRN)
        complex z
         integer iszrt(45)
        character*8 ostnm(2)
        character*8 ocmpnm1(15)
        character*3 ocmpnm2(3)

        character mname*80
        logical dostep
        common/zzy/datas
        complex datas(NFREQ)

        real svel, pvel, vsa, vsb, vsr, geom, rayp, tstar, den

        real caz, saz, c2az, s2az

        complex uz,ur,ut,duzdz,durdz,dufdz,
     1      duzdr,durdr,dutdr,duzdt,durdt,dufdt

        complex err,ert,erz,ett,etz,ezz
        complex xout(9)

        character ofile*80
        character oofile*80
        integer ls, lgstr
        

        data ostnm/'MT      ', 'F0      '/
        data ocmpnm1/'ZDD     ','RDD     ',
     1       'ZDS     ','RDS     ','TDS     ',
     1       'ZSS     ','RSS     ','TSS     ',
     1       'ZEX     ','REX     ',
     1       'ZVF     ','RVF     ',
     1       'ZHF     ','RHF     ','THF     '/
        data ocmpnm2/'   ','_dr','_dz'/
        data iszrt/1,4, 1,4,5, 1,4,5, 1,4, 1,4, 1,4,5,
     1             1,4, 1,4,5, 1,4,5, 1,4, 1,4, 1,4,5,
     1             1,4, 1,4,5, 1,4,5, 1,4, 1,4, 1,4,5/



        tau = ntau * dt
c-----
c       get first arrival time - NOTE for implementation of
c-----
            tTP  = -12345.
            tTSV = -12345.
            tTSH = -12345.
            CALL FRSTAR(dist,depths,depthr,
     1          MNAME,1,tTP ,
     1          SSA, SSC, SSF, SSL, SSN, SSR,
     1          RRA, RRC, RRF, RRL, RRN, RRR,
     1          rayp, geom, tstar, iiso, .false.)
            CALL FRSTAR(dist,depths,depthr,
     1          MNAME,2,tTSV,
     1          SSA, SSC, SSF, SSL, SSN, SSR,
     1          RRA, RRC, RRF, RRL, RRN, RRR,
     1          rayp, geom, tstar, iiso, .false.)
            CALL FRSTAR(dist,depths,depthr,
     1          MNAME,3,tTSH,
     1          SSA, SSC, SSF, SSL, SSN, SSR,
     1          RRA, RRC, RRF, RRL, RRN, RRR,
     1          rayp, geom, tstar, iiso, .false.)
c-----
c       Now that all spectra for a given distance are in
c       make time history and output in desired order
c-----
        np2 = nptin/2 + 1
        df = 1.0/(nptin*dt)
c-----
c       perform inverse FFT
c       If the number of data points is too large, then use supfft
c-----
c       here we will form the synthetics
c       there will be a total of 9, e.g.,
c       Uz  Ur Ut Err Ert Erz Ett Etz Ezz
c-----
c-----
c       jout   1  Uz
c              2  Ur
c              3  Ut
c              4  Err
c              5  Ert
c              6  Erz
c              7  Ett
c              8  Etz
c              9  Ezz
c----
c-----
c     now make 9 passes to make the synthetics
c     for ur uz ut err ert erz ett etz ezz
c-----
         do 1300 jk1=1,15
              do 1301 jk2=1,3
                jk = (jk2-1)*15 + jk1
                rewind 10
                do k=1,np2
                    read(10) (xx(i),i=1,NGRN)
                    datc(k)=xx(jk)
                    if(k.gt.1)then
                        datc(nptin+2-k)=conjg(datc(k))
                    endif
                enddo
c-----
c       perform inverse FFT
c-----
c       also take the opportunity to perform the correct
c       spherical spreading for an earth with radius of 6371.0 km
c-----
                    call zfour(datc,nptin,+1,dt,df)
                    fac = exp(alpha*tfirst)
                    dfac = exp(alpha*dt)
c-----
c     undamp 
c-----
                    do k=1,nptin
                        x(k) = real(datc(k))*fac
                        fac = fac*dfac
                    enddo
                    npts = nptin
c-----
c           convert velocity to displacement or acceleration
c           if required. Since we are using surface waves, the
c           DC component must be zero!
c
c           note the special case for the pressure field
c-----
C               if(jk.le.3)then
                IF(JK2.LE.1)THEN
                    if(idva.eq.0)then
                        call integ(x,npts,dt,.false.)
                    else if(idva.eq.2)then
                        call deriv(x,npts,dt)
                    endif
                else
                    call integ(x,npts,dt,.false.)
                endif
c-----
c     These orientation are in the Sac convention, e.g., Z is positive down here
c     Even though these are ZRT, the AZ is set to 0
c-----
                    if(iszrt(jk).eq.1)then
c-----
c               vertical
c-----
                        cmpinc = 180.0
                        cmpaz  =   0.0
c                     make vertical positive up
c             
                       do i=1,npts
                          x(i) = - x(i)
                       enddo
         
                    else if(iszrt(jk).eq.4)then
c-----
c               radial
c-----
                        cmpinc = 90.0
                        cmpaz  =   0.0
                    else if(iszrt(jk).eq.5)then
c-----
c               transverse
c-----
                        cmpinc = 90.0
                        cmpaz  =  90.0
                    endif
    
                cmpdt  = dt
                ksyear = 1970
                ksdoy  = 1
                ksmon  = 1
                ksday  = j
                kshour = 0
                ksmin  = 0
                ssec =  0.0

                o = 0
                distdg = -12345.
c
c-----
c-----
c     create the output file name
c-----
                 if(outfmt.eq.1)then
c                   DDDDDd_HHHh_ZZZz.cmp
                    ievdep = 10.*depths
                    istel  = 10.*depthr
                    idist  = 10.0*dist
                    write(ofile,11)idist,ievdep,istel
   11               format(i6.6,'_',i4.4,'_',i4.4,'.')
                 else if (outfmt.eq.2)then
c                   DDDDDddd_HHHhhh_ZZZzzz.cmp
                    ievdep = 1000.*depths
                    istel  = 1000.*depthr
                    idist  = 1000.*dist
                    write(ofile,12)idist,ievdep,istel
   12               format(i7.7,'_',i6.6,'_',i6.6,'.')
                 else if (outfmt.eq.3)then
c                   DDDDDdHHHh.grn
                    ievdep = 10.*depths
                    idist  = 10.0*dist
                    write(ofile,13)idist,ievdep
   13               format(i6.6,'_',i4.4,'.')
                 else if (outfmt.eq.4)then
c                   DDDDdHHHh.grn
                    ievdep = 10.*depths
                    idist  = 10.0*dist
                    write(ofile,14)idist,ievdep
   14               format(i5.5,'_',i4.4,'.')
                 else if (outfmt.eq.5)then
c                   DDDdddHhhh.grn
                    ievdep = 1000.*depths
                    idist  = 1000.*dist
                    write(ofile,15)idist,ievdep
   15               format(i6.6,'_',i4.4,'.')
                 endif
                ls = lgstr(ofile)
                oofile=ofile(1:ls)//ocmpnm1(jk1)(1:3)//ocmpnm2(jk2)(1:3)
                
                call scmxmn(x,npts,depmax,depmin,depmen,indmax,indmin)
                call newhdr()
                call setfhv('DELTA   ',cmpdt,nerr)
                    bb = tfirst
                    ee = bb + (npts -1 ) * cmpdt
                    call setfhv('A       ',tp    ,nerr)
                    call setkhv('KA      ','P       ',nerr)
                    call setfhv('T0      ',tsv   ,nerr)
                    call setfhv('T1      ',tsh   ,nerr)
                    if(iiso.eq.0)then
                       call setkhv('KT0     ','S       ',nerr)
                       call setkhv('KT1     ','S       ',nerr)
                    else
                       call setkhv('KT0     ','SV      ',nerr)
                       call setkhv('KT1     ','SH      ',nerr)
                    endif
                call setfhv('B       ',bb    ,nerr)
                call setfhv('E       ',ee    ,nerr)
                call setfhv('O       ',o     ,nerr)
                call setkhv('KO      ','O       ',nerr)
                call setfhv('STLA    ',stlat ,nerr)
                call setfhv('STLO    ',stlon ,nerr)
                call setfhv('TIMMAX  ',bb+ indmax*cmpdt,nerr)
                call setfhv('TIMMIN  ',bb+ indmin*cmpdt,nerr)
                call setfhv('DEPMIN  ',depmin,nerr)
                call setfhv('DEPMAX  ',depmax,nerr)
c-----
c               our elevation is KM, SAC is set to KILOMETERS
c-----
                call setfhv('STEL    ',depthr,nerr)
                call setfhv('EVLA    ',evlat ,nerr)
                call setfhv('EVLO    ',evlon ,nerr)
c-----
c               our depth is KM, SAC is set to KILOMETERS
c-----
                call setfhv('EVDP    ',depths ,nerr)
                call setfhv('DIST    ',dist,nerr)
                evstaz =   0.0
                stevaz = 180.0
                call setfhv('AZ      ',evstaz,nerr)
                call setfhv('BAZ     ',stevaz,nerr)
                call setfhv('GCARC   ',distdg,nerr)
                call setfhv('CMPAZ   ',cmpaz ,nerr)
                call setfhv('CMPINC  ',cmpinc,nerr)
                call setfhv('DEPMEN  ',depmen,nerr)
                call setfhv('STLA    ',-12345.,nerr)
                call setfhv('STLO    ',-12345.,nerr)
                call setfhv('EVLA    ',-12345.,nerr)
                call setfhv('EVLO    ',-12345.,nerr)
c-----
c       set integer header value
c-----
                call setnhv('NZYEAR  ',ksyear,nerr)
                call setnhv('NZJDAY  ',ksdoy ,nerr)
                call setnhv('NZHOUR  ',kshour,nerr)
                call setnhv('NZMIN   ',ksmin ,nerr)
                isec = ssec
                smsec = (ssec - isec) * 1000.0
                call setnhv('NZSEC   ',isec  ,nerr)
                jsmsec = smsec
                call setnhv('NZMSEC  ',jsmsec,nerr)
                call setnhv('NPTS    ',npts, nerr)
c-----
c       set logical header value
c-----
                call setlhv('LEVEN   ',.true.,nerr)
                call setlhv('LPSPOL  ',.true.,nerr)
                call setlhv('LOVROK  ',.true.,nerr)
                call setlhv('LCALDA  ',.false.,nerr)
c-----
c       set character header value
c-----
                call setkhv('KSTNM    ','GRN     ',nerr)
                if(iobsyn.eq.2)then
                    call setkhv('KEVNM   ','SYNTHETI',nerr)
                    call setkhv('KEVNMC  ','C       ',nerr)
                else
                    call setkhv('KEVNM   ','OBSERVED',nerr)
                    call setkhv('KEVNMC  ','        ',nerr)
                endif
                call setkhv('KCMPNM  ',
     1            ocmpnm1(jk1)(1:3)//ocmpnm2(jk2)(1:3)//'  ' ,nerr)
c-----
c       set enumerated header value
c-----
                call setihv('IFTYPE   ','ITIME   ',nerr)
                call setihv('IZTYPE   ','IB      ',nerr)
                ls = lgstr(oofile)
                    call bwsac(12,npts,oofile(1:ls),x)
 1301   continue
 1300   continue
        return
        end


        subroutine outstress(ntau,idva,dist,dt,nptin,
     1      tfirst,depthr,depths,az,
     2      mname,dostep,alpha,dosdr,TP,TSV,TSH)
c-----
c       write the Green s functions for this distance
c-----
c       ntau    I*4 - pulse duration parameter
c       idva    I*4 - displacmeent, velocity or acceleration
c       dist    R*4 - distance
c       dt  R*4 - sample rate
c       nptin   I*4 - number of points in time series
c       tfirst  R*4 - time of first sample
c       depthr  R*4 - receiver depth
c       depths  R*4 - source depth
c       az      R   - source to station azimuth
c       mname   C*(*)   - name of model file
c-----
        real stk, dip, rake,xmom
        logical dosdr

        character stcomp*8, stname*8
        real*4 cmpaz, cmpinc, cmpdt
        integer*4 npts
        real*4 ssec

        integer NL
        parameter (NL=200)
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL),
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax

        common/receiver/l2mu,mu
        real l2mu, mu
        common/grnfmt/outfmt
        integer outfmt

c-----
c       internal variables
c-----
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
        integer NGRN
        parameter (NGRN=45)
        integer*4 ntau, idva,   nptin
        real*4 dt, tfirst, depths, depthr

        integer kerr, system

        integer NSAMP, NFREQ
        parameter (NSAMP=16384,NFREQ=8193)
        common/dd/datc
        complex datc(NSAMP)
        common/xxx/x
        real*4 x(NSAMP)
        complex xx(NGRN)
        complex z
         integer iszrt(45)
        character*8 ostnm(2)
        character*8 ocmpnm(9)

        character mname*80
        logical dostep
        common/zzy/datas
        complex datas(NFREQ)

        real svel, pvel, vsa, vsb, vsr, geom, rayp, tstar, den

        real caz, saz, c2az, s2az

        complex uz,ur,ut,duzdz,durdz,dufdz,
     1      duzdr,durdr,dufdr,duzdf,durdf,dufdf

        complex srr,srf,srz,sff,sfz,szz
        complex err,efr,erz,eff,efz,ezz
        complex del
        complex xout(9)

        character ofile*80
        integer ls, lgstr
        

        data ostnm/'MT      ', 'F0      '/
        data ocmpnm/'Uz      ','Ur      ','Ut      ',
     1       'Srr     ', 'Srf     ',
     1       'Srz     ', 'Sff     ', 'Sfz     ', 'Szz     '/


        tau = ntau * dt
c-----
c       get first arrival time - NOTE for implementation of
c-----
            tTP  = -12345.
            tTSV = -12345.
            tTSH = -12345.
            CALL FRSTAR(dist,depths,depthr,
     1          MNAME,1,tTP ,
     1          SSA, SSC, SSF, SSL, SSN, SSR,
     1          RRA, RRC, RRF, RRL, RRN, RRR,
     1          rayp, geom, tstar, iiso, .false.)
            CALL FRSTAR(dist,depths,depthr,
     1          MNAME,2,tTSV,
     1          SSA, SSC, SSF, SSL, SSN, SSR,
     1          RRA, RRC, RRF, RRL, RRN, RRR,
     1          rayp, geom, tstar, iiso, .false.)
            CALL FRSTAR(dist,depths,depthr,
     1          MNAME,3,ttTSH,
     1          SSA, SSC, SSF, SSL, SSN, SSR,
     1          RRA, RRC, RRF, RRL, RRN, RRR,
     1          rayp, geom, tstar, iiso, .false.)
c-----
c       Now that all spectra for a given distance are in
c       make time history and output in desired order
c-----
        np2 = nptin/2 + 1
        df = 1.0/(nptin*dt)
c-----
c       perform inverse FFT
c       If the number of data points is too large, then use supfft
c-----
c       here we will form the synthetics
c       there will be a total of 9, e.g.,
c       Uz  Ur Ut Err Ert Erz Ett Etz Ezz
c-----
c-----
c       jout   1  Uz
c              2  Ur
c              3  Ut
c              4  Err
c              5  Ert
c              6  Erz
c              7  Ett
c              8  Etz
c              9  Ezz
c----
c       This is a two pass process
c       Starting with  ur uz tz durdz and other derivatives on LUN=10
c       compute ur uz ut and strains an then output them in LUN=9
c       
c-----
        open(9,file='file9.tmp',access='sequential', 
     1   form='unformatted', status='unknown')
        rewind 9
        rewind 10
        do  k=1,np2
              read(10) (xx(i),i=1,NGRN)
              call domom(xx,uz,ur,ut,duzdz,durdz,dufdz,
     1            duzdr,durdr,dufdr,duzdf,durdf,dufdf)
c-----
c     dist is in km 
c     the following is the strain formula
c     be careful of units especially the eff
c     Some derivatives are dU/dr or dU/dz obviously
c     we must convert from km to m
c-----
              r = dist
              err = durdr
              eff = (dufdf + ur)/r
              ezz = duzdz
              efr = 0.5*(durdf/r  + dufdr - ut/r)
              erz = 0.5*(duzdr + durdz)
              efz = 0.5*(dufdz + duzdf/r)

              del = err + eff + ezz
              srz = 2*RRL*erz
              sfz = 2*RRL*efz
              srf = 2*RRN*efr
              srr = (RRA -RRN -RRN)*del + 2*RRN*err
     1            + (RRF-RRA+2*RRN)*ezz
              sff = (RRA -RRN -RRN)*del + 2*RRN*eff
     1            + (RRF-RRA+2*RRN)*ezz
              szz = RRF*del + (RRC-RRF)*ezz
c-----
c     CONVERT FROM CM TO M
c-----
             xout(1) = uz/1.0e+02
             xout(2) = ur/1.0e+02
             xout(3) = ut/1.0e+02
c-----
c     CONVERT FROM GM/CM_^3 and KM/S to KG/M^3 and M/S
c     err = err/1.0E+05
c     srr = srr * 1.0e+09
c-----
             xout(4) = srr*1.0e+04
             xout(5) = srf*1.0e+04
             xout(6) = srz*1.0e+04
             xout(7) = sff*1.0e+04
             xout(8) = sfz*1.0e+04
             xout(9) = szz*1.0e+04
              write(9)(xout(i),i=1,9)
         enddo
c-----
c     now make 9 passes to make the synthetics
c     for ur uz ut err ert erz ett etz ezz
c-----
         rewind 9
         do 1300 jk=1,9
                rewind 9
                do k=1,np2
                    read(9) (xout(i),i=1,9)
                    datc(k)=xout(jk)
                    if(k.gt.1)then
                        datc(nptin+2-k)=conjg(datc(k))
                    endif
                enddo
c-----
c       perform inverse FFT
c-----
c       also take the opportunity to perform the correct
c       spherical spreading for an earth with radius of 6371.0 km
c-----
                    call zfour(datc,nptin,+1,dt,df)
                    fac = exp(alpha*tfirst)
                    dfac = exp(alpha*dt)
c-----
c     undamp 
c-----
                    do 1303 k=1,nptin
                        x(k) = real(datc(k))*fac
                        fac = fac*dfac
 1303               continue
                    npts = nptin
c-----
c           convert velocity to displacement or acceleration
c           if required. Since we are using surface waves, the
c           DC component must be zero!
c
c           note the special case for the pressure field
c-----
                IF(JK.LE.3)THEN
                    if(dostep)then
                       if(idva.eq.0)then
                           call integ(x,npts,dt,.false.)
                       else if(idva.eq.2)then
                           call deriv(x,npts,dt)
                       endif
                    else
                       if(idva.eq.1)then
                           call deriv(x,npts,dt)
                       else if(idva.eq.2)then
                           call deriv(x,npts,dt)
                           call deriv(x,npts,dt)
                       endif
                    endif
                else
                    if(dostep)then
                        call integ(x,npts,dt,.false.)
                    endif
                endif
c-----
c     component orientations SEED convention. Later 90 is added to 
c     make Sac convention. Difference is that Sac measures from upward
c     vertical and SEED from downward dip 
c            Component SEED     Sac
c     CMPAZ     Z        0       0
c     CMPINC           -90       0
c
c     CMPAZ     E       90      90
c     CMPINC             0      90
c
c     CMPAZ     N        0       0
c     CMPINC             0      90
c
c     These are Sac convention
c     Here the Uz is positive downward, e.g., SAC CMPINC 180
c-----
c           vertical
c-----
                if(jk.eq.1)then
                    cmpinc = 180.0 
                    cmpaz  =   0.0
c-----
c           radial
c-----
                else if(jk.eq.2)then
                    cmpinc =  90.0 
                    cmpaz  =   0.0 + az
c-----
c           transverse
c-----
                else if(jk.eq.3)then
                    cmpinc =  90.0 
                    cmpaz  =  amod(90.0 + az, 360.0)
                else
                    cmpinc = -12345.
                    cmpaz  = -12345.
                endif
                stcomp = 'SYN     '
    
                cmpdt  = dt
                ksyear = 1970
                ksdoy  = 1
                ksmon  = 1
                ksday  = j
                kshour = 0
                ksmin  = 0
                ssec =  0.0

                o = 0
                distdg = -12345.
c-----
c
c-----
                stcomp = ocmpnm(jk)(1:3)
c-----
c     create the output file name
c-----
                 if(outfmt.eq.1)then
c                   DDDDDd_HHHh_ZZZz.cmp
                    ievdep = 10.*depths
                    istel  = 10.*depthr
                    idist  = 10.0*dist
                    write(ofile,11)idist,ievdep,istel,stcomp
   11               format(i6.6,'_',i4.4,'_',i4.4,'.',a3)
                 else if (outfmt.eq.2)then
c                   DDDDDddd_HHHhhh_ZZZzzz.cmp
                    ievdep = 1000.*depths
                    istel  = 1000.*depthr
                    idist  = 1000.*dist
                    write(ofile,12)idist,ievdep,istel,stcomp
   12               format(i7.7,'_',i6.6,'_',i6.6,'.',a3)
                 else if (outfmt.eq.3)then
c                   DDDDDdHHHh.grn
                    ievdep = 10.*depths
                    idist  = 10.0*dist
                    write(ofile,13)idist,ievdep,stcomp
   13               format(i6.6,'_',i4.4,'.',a3)
                 else if (outfmt.eq.4)then
c                   DDDDdHHHh.grn
                    ievdep = 10.*depths
                    idist  = 10.0*dist
                    write(ofile,14)idist,ievdep,stcomp
   14               format(i5.5,'_',i4.4,'.',a3)
                 else if (outfmt.eq.5)then
c                   DDDdddHhhh.grn
                    ievdep = 1000.*depths
                    idist  = 1000.*dist
                    write(ofile,15)idist,ievdep,stcomp
   15               format(i6.6,'_',i4.4,'.',a3)
                 endif
                ls = lgstr(ofile)
                call scmxmn(x,npts,depmax,depmin,depmen,indmax,indmin)
                call newhdr()
                call setfhv('DELTA   ',cmpdt,nerr)
                    bb = tfirst
                    ee = bb + (npts -1 ) * cmpdt
                    call setfhv('A       ',tp    ,nerr)
                    call setkhv('KA      ','P       ',nerr)
                    call setfhv('T0      ',tsv   ,nerr)
                    call setfhv('T1      ',tsh   ,nerr)
                    if(iiso.eq.0)then
                       call setkhv('KT0     ','S       ',nerr)
                       call setkhv('KT1     ','S       ',nerr)
                    else
                       call setkhv('KT0     ','SV      ',nerr)
                       call setkhv('KT1     ','SH      ',nerr)
                    endif
                call setfhv('B       ',bb    ,nerr)
                call setfhv('E       ',ee    ,nerr)
                call setfhv('O       ',o     ,nerr)
                call setkhv('KO      ','O       ',nerr)
                call setfhv('STLA    ',stlat ,nerr)
                call setfhv('STLO    ',stlon ,nerr)
                call setfhv('TIMMAX  ',bb+ indmax*cmpdt,nerr)
                call setfhv('TIMMIN  ',bb+ indmin*cmpdt,nerr)
                call setfhv('DEPMIN  ',depmin,nerr)
                call setfhv('DEPMAX  ',depmax,nerr)
c-----
c               our elevation is KM, SAC is set to KILOMETERS
c-----
                call setfhv('STEL    ',depthr,nerr)
                call setfhv('EVLA    ',evlat ,nerr)
                call setfhv('EVLO    ',evlon ,nerr)
c-----
c               our depth is KM, SAC is set to KILOMETERS
c-----
                call setfhv('EVDP    ',depths ,nerr)
                call setfhv('DIST    ',dist,nerr)
                evstaz =   az
                stevaz = amod(180.0+az, 360.0)
                call setfhv('AZ      ',evstaz,nerr)
                call setfhv('BAZ     ',stevaz,nerr)
                call setfhv('GCARC   ',distdg,nerr)
                call setfhv('CMPAZ   ',cmpaz ,nerr)
                call setfhv('CMPINC  ',cmpinc,nerr)
                call setfhv('DEPMEN  ',depmen,nerr)
                call setfhv('STLA    ',-12345.,nerr)
                call setfhv('STLO    ',-12345.,nerr)
                call setfhv('EVLA    ',-12345.,nerr)
                call setfhv('EVLO    ',-12345.,nerr)
c-----
c       set integer header value
c-----
                call setnhv('NZYEAR  ',ksyear,nerr)
                call setnhv('NZJDAY  ',ksdoy ,nerr)
                call setnhv('NZHOUR  ',kshour,nerr)
                call setnhv('NZMIN   ',ksmin ,nerr)
                isec = ssec
                smsec = (ssec - isec) * 1000.0
                call setnhv('NZSEC   ',isec  ,nerr)
                jsmsec = smsec
                call setnhv('NZMSEC  ',jsmsec,nerr)
                call setnhv('NPTS    ',npts, nerr)
c-----
c       set logical header value
c-----
                call setlhv('LEVEN   ',.true.,nerr)
                call setlhv('LPSPOL  ',.true.,nerr)
                call setlhv('LOVROK  ',.true.,nerr)
                call setlhv('LCALDA  ',.false.,nerr)
c-----
c       set character header value
c-----
                call setkhv('KSTNM    ','SYN     ' ,nerr)
                if(iobsyn.eq.2)then
                    call setkhv('KEVNM   ','SYNTHETI',nerr)
                    call setkhv('KEVNMC  ','C       ',nerr)
                else
                    call setkhv('KEVNM   ','OBSERVED',nerr)
                    call setkhv('KEVNMC  ','        ',nerr)
                endif
                call setkhv('KCMPNM   ',stcomp ,nerr)
c-----
c       set enumerated header value
c-----
                call setihv('IFTYPE   ','ITIME   ',nerr)
                call setihv('IZTYPE   ','IB      ',nerr)
                    call bwsac(12,npts,ofile(1:ls),x)
 1300   continue
        close(9,status='delete')
        return
        end

        subroutine outrotate(ntau,idva,dist,dt,nptin,
     1      tfirst,depthr,depths,az,
     2      mname,dostep,alpha,dosdr,TP,TSV,TSH)
c-----
c       write the Green s functions for this distance
c-----
c       ntau    I*4 - pulse duration parameter
c       idva    I*4 - displacmeent, velocity or acceleration
c       dist    R*4 - distance
c       dt  R*4 - sample rate
c       nptin   I*4 - number of points in time series
c       tfirst  R*4 - time of first sample
c       depthr  R*4 - receiver depth
c       depths  R*4 - source depth
c       az      R   - source to station azimuth
c       mname   C*(*)   - name of model file
c-----
        real stk, dip, rake,xmom
        logical dosdr

        character stcomp*8, stname*8
        real*4 cmpaz, cmpinc, cmpdt
        integer*4 npts
        real*4 ssec

        integer NL
        parameter (NL=200)
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL),
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax

        common/receiver/l2mu,mu
        real l2mu, mu
        common/grnfmt/outfmt
        integer outfmt

c-----
c       internal variables
c-----
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
        integer NGRN
        parameter (NGRN=45)
        integer*4 ntau, idva,   nptin
        real*4 dt, tfirst, depths, depthr

        integer kerr, system

        integer NSAMP, NFREQ
        parameter (NSAMP=16384,NFREQ=8193)
        common/dd/datc
        complex datc(NSAMP)
        common/xxx/x
        real*4 x(NSAMP)
        complex xx(NGRN)
        complex z
         integer iszrt(45)
        character*8 ostnm(2)
        character*8 ocmpnm(3)

        character mname*80
        logical dostep, dostrain
        common/zzy/datas
        complex datas(NFREQ)

        real svel, pvel, vsa, vsb, vsr, geom, rayp, tstar, den

        real caz, saz, c2az, s2az

        complex uz,ur,ut,duzdz,durdz,dufdz,
     1      duzdr,durdr,dufdr,duzdf,durdf,dufdf

        complex wfz, wrz, wrf
        complex xout(9)

        character ofile*80
        integer ls
        

        data ostnm/'MT      ', 'F0      '/
        data ocmpnm/ 'Wfz     ','Wrz     ','Wrf     '/


        tau = ntau * dt
        if(idva.eq.0)then
            iunit = 2
        else if(idva.eq.1)then
            iunit = 3
        else if(idva.eq.2)then
            iunit = 4
        endif
c-----
c       get first arrival time - NOTE for implementation of
c-----
            tTP  = -12345.
            tTSV = -12345.
            tTSH = -12345.
            CALL FRSTAR(dist,depths,depthr,
     1          MNAME,1,tTP ,
     1          SSA, SSC, SSF, SSL, SSN, SSR,
     1          RRA, RRC, RRF, RRL, RRN, RRR,
     1          rayp, geom, tstar, iiso, .false.)
            CALL FRSTAR(dist,depths,depthr,
     1          MNAME,2,tTSV,
     1          SSA, SSC, SSF, SSL, SSN, SSR,
     1          RRA, RRC, RRF, RRL, RRN, RRR,
     1          rayp, geom, tstar, iiso, .false.)
            CALL FRSTAR(dist,depths,depthr,
     1          MNAME,3,tTSH,
     1          SSA, SSC, SSF, SSL, SSN, SSR,
     1          RRA, RRC, RRF, RRL, RRN, RRR,
     1          rayp, geom, tstar, iiso, .false.)
c-----
c       Now that all spectra for a given distance are in
c       make time history and output in desired order
c-----
        np2 = nptin/2 + 1
        df = 1.0/(nptin*dt)
c-----
c       perform inverse FFT
c       If the number of data points is too large, then use supfft
c-----
c       here we will form the synthetics
c-----
c-----
c       jout   1  Wfz
c              2  Wrz
c              3  Wrf
c----
c       This is a two pass process
c       Starting with  ur uz tz durdz and other derivatives on LUN=10
c       compute ur uz ut and strains an then output them in LUN=9
c       
c-----
        open(9,file='file9.tmp',access='sequential', 
     1   form='unformatted', status='unknown')
        rewind 9
        rewind 10
        do  k=1,np2
              read(10) (xx(i),i=1,NGRN)
              call domom(xx,uz,ur,ut,duzdz,durdz,dufdz,
     1            duzdr,durdr,dufdr,duzdf,durdf,dufdf)
c-----
c     dist is in km 
c     the following is the strain formula
c     be careful of units especially the eff
c     Some derivatives are dU/dr or dU/dz obviously
c     we must convert from km to m
c-----
              r = dist
              wfz = 0.5*(dufdz - duzdf/r)
              wrz = 0.5*(durdz - duzdr  )
              wrf = 0.5*(durdf/r - dufdr - ut/r)
c-----
c     CONVERT FROM KM TO M for distance and CM to M for displacement
c-----
             xout(1) = wfz/1.0e+05
             xout(2) = wrz/1.0e+05
             xout(3) = wrf/1.0e+05
              write(9)(xout(i),i=1,3)
         enddo
c-----
c     now make 3 passes to make the synthetics
c     for wrf wrz wfz
c-----
         rewind 9
         do 1300 jk=1,3
                rewind 9
                do k=1,np2
                    read(9) (xout(i),i=1,3)
                    datc(k)=xout(jk)
                    if(k.gt.1)then
                        datc(nptin+2-k)=conjg(datc(k))
                    endif
                enddo
c-----
c       perform inverse FFT
c-----
c       also take the opportunity to perform the correct
c       spherical spreading for an earth with radius of 6371.0 km
c-----
                    call zfour(datc,nptin,+1,dt,df)
                    fac = exp(alpha*tfirst)
                    dfac = exp(alpha*dt)
c-----
c     undamp 
c-----
                    do 1303 k=1,nptin
                        x(k) = real(datc(k))*fac
                        fac = fac*dfac
 1303               continue
                    npts = nptin
c-----
c           convert velocity to displacement or acceleration
c           if required. Since we are using surface waves, the
c           DC component must be zero!
c
c           note the special case for the pressure field
c-----
                    if(dostep)then
                        call integ(x,npts,dt,.false.)
                    endif
                    cmpinc = -12345.
                    cmpaz  = -12345.
                stcomp = 'SYN     '
    
                cmpdt  = dt
                ksyear = 1970
                ksdoy  = 1
                ksmon  = 1
                ksday  = j
                kshour = 0
                ksmin  = 0
                ssec =  0.0

                o = 0
                distdg = -12345.
c-----
c
c-----
                stcomp = ocmpnm(jk)(1:3)
c-----
c     create the output file name
c-----
                 if(outfmt.eq.1)then
c                   DDDDDd_HHHh_ZZZz.cmp
                    ievdep = 10.*depths
                    istel  = 10.*depthr
                    idist  = 10.0*dist
                    write(ofile,11)idist,ievdep,istel,stcomp
   11               format(i6.6,'_',i4.4,'_',i4.4,'.',a3)
                 else if (outfmt.eq.2)then
c                   DDDDDddd_HHHhhh_ZZZzzz.cmp
                    ievdep = 1000.*depths
                    istel  = 1000.*depthr
                    idist  = 1000.*dist
                    write(ofile,12)idist,ievdep,istel,stcomp
   12               format(i7.7,'_',i6.6,'_',i6.6,'.',a3)
                 else if (outfmt.eq.3)then
c                   DDDDDdHHHh.grn
                    ievdep = 10.*depths
                    idist  = 10.0*dist
                    write(ofile,13)idist,ievdep,stcomp
   13               format(i6.6,i4.4,'.',a3)
                 else if (outfmt.eq.4)then
c                   DDDDdHHHh.grn
                    ievdep = 10.*depths
                    idist  = 10.0*dist
                    write(ofile,14)idist,ievdep,stcomp
   14               format(i5.5,i4.4,'.',a3)
                 else if (outfmt.eq.5)then
c                   DDDdddHhhh.grn
                    ievdep = 1000.*depths
                    idist  = 1000.*dist
                    write(ofile,15)idist,ievdep,stcomp
   15               format(i6.6,i4.4,'.',a3)
                 endif
                 ls = lgstr(ofile)
                call scmxmn(x,npts,depmax,depmin,depmen,indmax,indmin)
                call newhdr()
                call setfhv('DELTA   ',cmpdt,nerr)
                    bb = tfirst
                    ee = bb + (npts -1 ) * cmpdt
                    call setfhv('A       ',tp    ,nerr)
                    call setkhv('KA      ','P       ',nerr)
                    call setfhv('T0      ',tsv   ,nerr)
                    call setfhv('T1      ',tsh   ,nerr)
                    if(iiso.eq.0)then
                       call setkhv('KT0     ','S       ',nerr)
                       call setkhv('KT1     ','S       ',nerr)
                    else
                       call setkhv('KT0     ','SV      ',nerr)
                       call setkhv('KT1     ','SH      ',nerr)
                    endif
                call setfhv('B       ',bb    ,nerr)
                call setfhv('E       ',ee    ,nerr)
                call setfhv('O       ',o     ,nerr)
                call setkhv('KO      ','O       ',nerr)
                call setfhv('STLA    ',stlat ,nerr)
                call setfhv('STLO    ',stlon ,nerr)
                call setfhv('TIMMAX  ',bb+ indmax*cmpdt,nerr)
                call setfhv('TIMMIN  ',bb+ indmin*cmpdt,nerr)
                call setfhv('DEPMIN  ',depmin,nerr)
                call setfhv('DEPMAX  ',depmax,nerr)
c-----
c               our elevation is KM, SAC is set to KILOMETERS
c-----
                call setfhv('STEL    ',depthr,nerr)
                call setfhv('EVLA    ',evlat ,nerr)
                call setfhv('EVLO    ',evlon ,nerr)
c-----
c               our depth is KM, SAC is set to KILOMETERS
c-----
                call setfhv('EVDP    ',depths ,nerr)
                call setfhv('DIST    ',dist,nerr)
                evstaz =   az
                stevaz = amod(180.0+az, 360.0)
                call setfhv('AZ      ',evstaz,nerr)
                call setfhv('BAZ     ',stevaz,nerr)
                call setfhv('GCARC   ',distdg,nerr)
                call setfhv('CMPAZ   ',cmpaz ,nerr)
                call setfhv('CMPINC  ',cmpinc,nerr)
                call setfhv('DEPMEN  ',depmen,nerr)
                call setfhv('STLA    ',-12345.,nerr)
                call setfhv('STLO    ',-12345.,nerr)
                call setfhv('EVLA    ',-12345.,nerr)
                call setfhv('EVLO    ',-12345.,nerr)
c-----
c       set integer header value
c-----
                call setnhv('NZYEAR  ',ksyear,nerr)
                call setnhv('NZJDAY  ',ksdoy ,nerr)
                call setnhv('NZHOUR  ',kshour,nerr)
                call setnhv('NZMIN   ',ksmin ,nerr)
                isec = ssec
                smsec = (ssec - isec) * 1000.0
                call setnhv('NZSEC   ',isec  ,nerr)
                jsmsec = smsec
                call setnhv('NZMSEC  ',jsmsec,nerr)
                call setnhv('NPTS    ',npts, nerr)
c-----
c       set logical header value
c-----
                call setlhv('LEVEN   ',.true.,nerr)
                call setlhv('LPSPOL  ',.true.,nerr)
                call setlhv('LOVROK  ',.true.,nerr)
                call setlhv('LCALDA  ',.false.,nerr)
c-----
c       set character header value
c-----
                call setkhv('KSTNM    ','SYN     ' ,nerr)
                if(iobsyn.eq.2)then
                    call setkhv('KEVNM   ','SYNTHETI',nerr)
                    call setkhv('KEVNMC  ','C       ',nerr)
                else
                    call setkhv('KEVNM   ','OBSERVED',nerr)
                    call setkhv('KEVNMC  ','        ',nerr)
                endif
                call setkhv('KCMPNM   ',stcomp ,nerr)
c-----
c       set enumerated header value
c-----
                call setihv('IFTYPE   ','ITIME   ',nerr)
                call setihv('IZTYPE   ','IB      ',nerr)
                    call bwsac(12,npts,ofile(1:ls),x)
 1300   continue
        close(9,status='delete')
        return
        end

        subroutine integ(x,nt,dt,dodc)
        parameter (LER=0, LIN=5, LOT=6)
c-----
c       Form the integral x(t) dt numerically
c       x   R*4 - array to be integrated
c       nt  I*4 - number of data points in time series
c       dt  R*4 - sampling interval
c       dodc    L*4 - .true. remve x(1) before summing
c-----
        
        real x(nt)
        integer*4 nt
        real*4 dt
        logical dodc
        sum=0.0
        if(dodc)then
            do 100 k = 1,nt
                sum = sum + dt*(x(k) - x(1))
                x(k)=sum
  100       continue
        else
            do 200 k = 1,nt
                sum = sum + dt*x(k)
                x(k)=sum
  200       continue
        endif
        return
        end

        subroutine srcpul(ntau,ipt,alp,dt,npts,rfile,duration,
     1     alpha)
c-----
c       define the source time function
c-----
c       ntau    I*4 - pulse duration factor for parabolic and triangular
c               pulses
c       ipt I*4 - pulse type
c               0 - triangular
c               1 - parabolic
c               2 - Ohnaka
c               3 - Dirac delta function
c               4 = user pulse in file rfile
c       alp R*4 - Ohnaka pulse shape factor, fc= alp / 2 pi
c       dt  R*4 - sample rate for time series
c       npts    I*4 - number of samples in time series
c       rfile   C*80- name of user provided pulse 
c               note dt must be the same as specified in dfile
c       duration R*4 - duration of triangular/parabolic source pulse
c       datas   Complex - source spectra
c       alp     R - time domain damping factor
c-----
        parameter (LIN=5, LOT=6, LER=0)
        integer*4 ntau, ipt, npts
        real alp, dt, alpha
      
        character*80 rfile
        integer NSAMP, NFREQ
        parameter (NSAMP=16384,NFREQ=8193)
        common/srctim/ src(NSAMP)
        common/zzy/datas
        complex datas(NFREQ)
        do 50 i=1,NSAMP
            src(i)=0.0
   50   continue
        if(ipt.eq.0)then
c-----
c           triangluar
c-----
            tau = ntau * dt
            nt = 2*ntau + 1
            call pultd(tau,dt,nt,l)
            duration = 2.*ntau*dt
        else if(ipt.eq.1)then
c-----
c           parabolic
c-----
            tau = ntau * dt
            nt = 4*ntau + 1
            call pulpd(tau,dt,nt,l)
            duration = 4.*ntau * dt
        else if(ipt.eq.2)then
c-----
c           Ohnaka
c-----
            if(npts.gt.NSAMP)then
                nt=NSAMP
            else
                nt = npts
            endif
            call  pulod(alp,dt,nt,l)
            duration = 0.0
        else if(ipt.eq.3)then
c-----
c           Dirac
c-----
            if(npts.gt.NSAMP)then
                nt=NSAMP
            else
                nt = npts
            endif
            call puldd(dt,nt,l)
            duration = 0.0
        else if(ipt.eq.4)then
c-----
c           User specified
c-----
            if(npts.gt.NSAMP)then
                nt=NSAMP
            else
                nt = npts
            endif
            call pulud(rfile,nt,tau,ddt)
                if(ddt.ne.dt)then
                    write(LER,*)'Warning: RFILE dt',ddt,
     1              ' is not same as dt',dt,
     1                  ' specified for synthetic'
                endif

            ntau = tau/dt
            duration = 0.0
        endif
c-----
c       now dampen
c-----
        fac = 1.0
        dfac = exp(-alpha*dt)
        do i=1,npts
             if(i.le.nt)then
             datas(i) = cmplx(fac*src(i),0.0)
             else
             datas(i) = cmplx(0.0,0.0)
             endif
             fac = fac * dfac
        enddo
        call zfour(datas,npts,-1,dt,df)
        return
        end

        subroutine puldd(dt,n,l)
        integer NSAMP
        parameter (NSAMP=16384)
        common/srctim/d(NSAMP)
c-----
c       Dirac Delta Pulse
c-----
        do 100 i=1,n
            d(i) = 0.0
  100   continue
        d(1) = 1.0/dt
        l = 2
        return
        end

        subroutine pulod(alp,dt,nt,l)
c-----
c       unit area far field displacement Ohnaka pulse
c           Harkrider (1976) Geophys J. 47, p 97.
c-----
c       alp R*4 - shape parameter, corner frequency
c               fc = alp / 2 pi
c       dt  R*4 - sample rate
c       nt  I*4 - number of points in time series
c       l   I*4 - pulse d(i) = 0 for i >=l
c-----
        real*4 alp, dt
        integer*4 nt, l
        integer NSAMP
        parameter (NSAMP=16384)
        common/srctim/d(NSAMP)
        ltest = 0
        al2=alp*alp
        do 100 i=1,nt
            t=(i-1)*dt
            d(i)=0.0
            arg= alp*t
            if(arg.le.25.0)then
                d(i)= al2*t*exp(-arg)
            else
                ltest = ltest +1
                if(ltest.eq.1)l = i
            endif
  100   continue
        return
        end

        subroutine pulpd(tau,dt,nt,l)
c-----
c       unit area far field displacement parabolic pulse
c-----
c       tau R*4 - duration parameter, total duration = 4 tau
c       dt  R*4 - sample rate
c       nt  I*4 - number of points for time series
c       l   I*4 - pulse d(i) = 0 for i >=l
c-----
        real*4 tau, dt
        integer*4 nt, l
        integer NSAMP
        parameter (NSAMP=16384)
        common /srctim/ d(NSAMP)
        ltest = 0
        tl = tau
        t1 = 0.01*dt
        t2 = t1 + tau
        t3 = t2 + tau
        t4 = t3 + tau
        t5 = t4 + tau
        do 100 i = 1,nt
            y=(i-1)*dt
            z = y - t1
            d(i) = 0.0
            if(y.ge.t1 .and. y.lt.t2)then
                d(i) = 0.5*(z/tl)**2
            else if(y.ge.t2 .and. y.le.t3)then
                d(i)= -0.5*(z/tl)**2 + 2.*(z/tl) - 1
            else if(y.ge.t3 .and. y.le.t4)then
                d(i)= -0.5*(z/tl)**2 + 2.*(z/tl) - 1.
            else if(y.ge.t4 .and. y.le.t5)then
                d(i)= 0.5*(z/tl)**2 - 4.*(z/tl) + 8.
            else
                ltest = ltest + 1
                if(ltest.eq.1) l = i
            endif
  100   continue
c-----
c       pulse normalized so first integral has area of unity
c-----
        do 200 i = 1,nt
            d(i) = d(i)/(2.*tl)
  200   continue
        return
        end

        subroutine pultd(tau,dt,nt,l)
c-----
c       unit area far field displacement triangular pulse
c-----
c       tau R*4 - duration parameter, total duration = 2 tau
c       dt  R*4 - sample rate
c       nt  I*4 - number of points for time series
c       l   I*4 - pulse d(i) = 0 for i >=l
c-----
        real*4 tau, dt
        integer*4 nt, l
        integer NSAMP
        parameter (NSAMP=16384)
        common /srctim/ d(NSAMP)
        ltest = 0
        fac = 1./tau
        t1 = tau
        t2 = tau + tau
        do 100 i=1,nt
            t = (i-1)*dt
            d(i)=0.0
            if(t.le.t1)then
                z = t - 0.0
                d(i) = z*fac
            elseif(t.gt.t1.and.t.le.t2)then
                z = t - t1
                d(i)= fac - z*fac
            elseif(t.gt.t2)then
                d(i)=0.0
                ltest = ltest + 1
                if(ltest.eq.1)l = i
            endif
  100   continue
        return
        end

        subroutine pulud(rfile,n,tau,dtt)
        character rfile*(*)
        integer NSAMP
        parameter (NSAMP=16384)
        common/srctim/d(NSAMP)
        do 100 i=1,n
            d(i) = 0.0
  100   continue
        open(4,file=rfile,status='unknown',form='formatted',
     1      access='sequential')
        rewind 4
        read(4,*)np,dtt
        read(4,*)(d(i),i=1,np)
        close(4)
        tau = np*dtt
        return
        end

        subroutine deriv(y,n,dt)
c-----
c       Use centered difference to take derivative of time series
c-----
c       y   R*4 - time series to be differentiated and returned as y
c       n   I*4 - length of time series
c       dt  R*4 - sample interval
c-----
        real*4 y(1)
        b0=1.0/dt
        xm1=y(1)
        xm2=0.
        y1=0.
        do 11 i=1,n
            y2=b0*(y(i)-xm1) 
            xm2=xm1
            xm1=y(i)
            y(i)=y2
   11   continue
        return
        end

      subroutine domom(xx,uz,ur,ut,duzdz,durdz,dutdz,
     1      duzdr,durdr,dutdr,duzdt,durdt,dutdt)
c-----
c     apply the moment tensor/force to this set of Green's functions
c     at one frequency
c-----
        implicit none
        integer NGRN
        parameter (NGRN=45)
        complex xx(NGRN)
        complex uz,ur,ut,duzdz,durdz,dutdz,
     1      duzdr,durdr,dutdr,duzdt,durdt,dutdt

        common/stresstrain/dostrain,dostress,dorotate,dogreen,
     1     az, baz, oxmt(3,3), xmom,fx, fy,fz,
     1     caz, saz, c2az, s2az, stk, dip, rake
        logical dostrain, dostress,dorotate,dogreen
        real az, baz, oxmt, xmom,fx, fy,fz
        real caz, saz, c2az, s2az, stk, dip, rake

        real xmt(3,3)
        real forcex, forcey, forcez
        integer i,j
c-----
c     the input order from hspec96strain is
c      1  ZDD   16  d ZDD/dz   31 d ZDD/dr
c      2  RDS
c      3  ZDS
c      4  RDS
c      5  ZSS
c      6  RSS
c      7  ZEX
c      8  REX
c      9  ZVF
c     10  RVF
c     11  ZHF
c     12  RHF
c     13  TDS
c     14  TSS
c     15  THF   30  d THF/dz   45 d THF/dr
c-----
c       assume that the model in KM KM/S GM/CM^3
c       and that the force is in dyne and the moment in dyne-cm
c
c-----
        do i=1,3
           do j=1,3
              xmt(i,j) = oxmt(i,j)/1.0e+20
           enddo
        enddo
        forcex = fx/1.0e+15
        forcey = fy/1.0e+15
        forcez = fz/1.0e+15

        uz = forcex*caz*xx(13) + forcey*saz*xx(13) + forcez*xx(11)
     1     + xmt(1,1)*( c2az*xx(6)/2. - xx(1)/6. + xx(9)/3.)
     1     + xmt(2,2)*(-c2az*xx(6)/2. - xx(1)/6. + xx(9)/3.)
     1     + xmt(3,3)*( xx(1)/3. + xx(9)/3.)
     1     + xmt(1,2)*s2az*xx(6)
     1     + xmt(1,3)* caz*xx(3)
     1     + xmt(2,3)* saz*xx(3)
        ur = forcex*caz*xx(14) + forcey*saz*xx(14) + forcez*xx(12)
     1     + xmt(1,1)*( c2az*xx(7)/2. - xx(2)/6. + xx(10)/3.)
     1     + xmt(2,2)*(-c2az*xx(7)/2. - xx(2)/6. + xx(10)/3.)
     1     + xmt(3,3)*( xx(2)/3. + xx(10)/3.)
     1     + xmt(1,2)*s2az*xx(7)
     1     + xmt(1,3)* caz*xx(4)
     1     + xmt(2,3)* saz*xx(4)
c check the force since this is old?
c-----
c       note this does not follow text because of the way that THF
c       is defined. The difference is a minus sign
c-----
        ut = ( forcex*saz-forcey*caz)*xx(15)
     1     + xmt(1,1)*s2az*xx(8)/2.
     1     - xmt(2,2)*s2az*xx(8)/2.
     1     - xmt(1,2)*c2az*xx(8)
     1     + xmt(1,3)*saz*xx(5)
     1     - xmt(2,3)*caz*xx(5)

        duzdr = forcex*caz*xx(13+15) + forcey*saz*xx(13+15) 
     1     + forcez*xx(11+15)
     1     + xmt(1,1)*( c2az*xx(6+15)/2. - xx(1+15)/6. + xx(9+15)/3.)
     1     + xmt(2,2)*(-c2az*xx(6+15)/2. - xx(1+15)/6. + xx(9+15)/3.)
     1     + xmt(3,3)*( xx(1+15)/3. + xx(9+15)/3.)
     1     + xmt(1,2)*s2az*xx(6+15)
     1     + xmt(1,3)* caz*xx(3+15)
     1     + xmt(2,3)* saz*xx(3+15)
        durdr = forcex*caz*xx(14+15) + forcey*saz*xx(14+15) 
     1     + forcez*xx(12+15)
     1     + xmt(1,1)*( c2az*xx(7+15)/2. - xx(2+15)/6. + xx(10+15)/3.)
     1     + xmt(2,2)*(-c2az*xx(7+15)/2. - xx(2+15)/6. + xx(10+15)/3.)
     1     + xmt(3,3)*( xx(2+15)/3. + xx(10+15)/3.)
     1     + xmt(1,2)*s2az*xx(7+15)
     1     + xmt(1,3)*xx(4+15)*caz
     1     + xmt(2,3)*xx(4+15)*saz
        dutdr = ( forcex*saz-forcey*caz)*xx(15+15)
     1     + xmt(1,1)*s2az*xx(8+15)/2.
     1     - xmt(2,2)*s2az*xx(8+15)/2.
     1     - xmt(1,2)*c2az*xx(8+15)
     1     + xmt(1,3)*saz*xx(5+15)
     1     - xmt(2,3)*caz*xx(5+15)

        duzdz = forcex*caz*xx(13+30) + forcey*saz*xx(13+30) 
     1     + forcez*xx(11+30)
     1     + xmt(1,1)*(c2az*xx(6+30)/2. - xx(1+30)/6. + xx(9+30)/3.)
     1     + xmt(2,2)*(-c2az*xx(6+30)/2. - xx(1+30)/6. + xx(9+30)/3.)
     1     + xmt(3,3)*( xx(1+30)/3. + xx(9+30)/3.)
     1     + xmt(1,2)*s2az*xx(6+30)
     1     + xmt(1,3)*xx(3+30)*caz
     1     + xmt(2,3)*xx(3+30)*saz
        durdz = forcex*caz*xx(14+30) + forcey*saz*xx(14+30) 
     1     + forcez*xx(12+30)
     1     + xmt(1,1)*(c2az*xx(7+30)/2. - xx(2+30)/6. + xx(10+30)/3.)
     1     + xmt(2,2)*(-c2az*xx(7+30)/2. - xx(2+30)/6. + xx(10+30)/3.)
     1     + xmt(3,3)*( xx(2+30)/3. + xx(10+30)/3.)
     1     + xmt(1,2)*s2az*xx(7+30)
     1     + xmt(1,3)*xx(4+30)*caz
     1     + xmt(2,3)*xx(4+30)*saz
        dutdz = ( forcex*saz-forcey*caz)*xx(15+30)
     1     + xmt(1,1)*s2az*xx(8+30)/2.
     1     - xmt(2,2)*s2az*xx(8+30)/2.
     1     - xmt(1,2)*c2az*xx(8+30)
     1     + xmt(1,3)*saz*xx(5+30)
     1     - xmt(2,3)*caz*xx(5+30)

        duzdt = -forcex*saz*xx(13) + forcey*caz*xx(13) 
     1     + xmt(1,1)*(-2*s2az*xx(6)/2. )
     1     + xmt(2,2)*( 2*s2az*xx(6)/2. )
     1     + xmt(1,2)*2*c2az*xx(6)
     1     - xmt(1,3)*xx(3)*saz
     1     + xmt(2,3)*xx(3)*caz
        durdt = -forcex*saz*xx(14) + forcey*caz*xx(14) 
     1     + xmt(1,1)*(-2*s2az*xx(7)/2. )
     1     + xmt(2,2)*( 2*s2az*xx(7)/2. )
     1     + xmt(1,2)*2*c2az*xx(7)
     1     - xmt(1,3)*xx(4)*saz
     1     + xmt(2,3)*xx(4)*caz
        dutdt = ( forcex*caz+forcey*saz)*xx(15)
     1     + xmt(1,1)*2*c2az*xx(8)/2.
     1     - xmt(2,2)*2*c2az*xx(8)/2.
     1     + xmt(1,2)*2*s2az*xx(8)
     1     + xmt(1,3)*caz*xx(5)
     1     + xmt(2,3)*saz*xx(5)




        return
        end


        subroutine outgreen(ntau,idva,dist,dt,nptin,
     1               tfirst,depthr,depths,
     2               mname,dostep,alpha,TP,TSV,TSH)
c-----
c       write the Green s functions for this distance
c-----
c       ntau    I*4 - pulse duration parameter
c       idva    I*4 - displacmeent, velocity or acceleration
c       dist    R*4 - distance
c       dt  R*4 - sample rate
c       nptin   I*4 - number of points in time series
c       tfirst  R*4 - time of first sample
c       depthr  R*4 - receiver depth
c       depths  R*4 - source depth
c       mname   C*(*)   - name of model file
c       logical outstep
c       alpha   R   - time domain damping factor
c-----
        character stcomp*8, staname*8
        real*4 cmpaz, cmpinc, cmpdt
        integer*4 npts
        real*4 ssec

        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax

        common/receiver/l2mu,mu
        real l2mu, mu
        common/grnfmt/outfmt
        integer outfmt

c-----
c       internal variables
c-----
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
        integer NGRN
        parameter (NGRN=45)
        integer*4 ntau, idva,   nptin
        real*4 dt, tfirst, depths, depthr

        integer kerr, system

        integer NSAMP2
        parameter(NSAMP2=16384)
        common/dd/datc
        complex datc(NSAMP2)
        common/xxx/x
        real*4 x(NSAMP2)
        complex xx(NGRN)
        complex z
         integer iszrt(15)
        character*8 ostnm(1)
        character*8 ocmpnm(15)

        character mname*80
        logical dostep

        real svel, pvel, vsa, vsb, vsr, geom, rayp, tstar, den

        complex xout(NGRN)

        character ofile*80
        integer ls

        data iszrt/1,4, 1,4,5, 1,4,5, 1,4, 1,4, 1,4,5/
        

        data ostnm/'GRN     '/
        data ocmpnm/'ZDD     ','RDD     ',
     1       'ZDS     ','RDS     ', 'TDS     ',
     1       'ZSS     ','RSS     ', 'TSS     ',
     1       'ZEX     ','REX     ', 
     1       'ZVF     ','RVF     ', 
     1       'ZHF     ','RHF     ', 'THF     '/


        tau = ntau * dt
c-----
c       get first arrival time - NOTE for implementation of
c-----
            tTP  = -12345.
            tTSV = -12345.
            tTSH = -12345.
            CALL FRSTAR(dist,depths,depthr,
     1          MNAME,1,tTP ,
     1          SSA, SSC, SSF, SSL, SSN, SSR,
     1          RRA, RRC, RRF, RRL, RRN, RRR,
     1          rayp, geom, tstar, iiso, .false.)
            CALL FRSTAR(dist,depths,depthr,
     1          MNAME,2,tTSV,
     1          SSA, SSC, SSF, SSL, SSN, SSR,
     1          RRA, RRC, RRF, RRL, RRN, RRR,
     1          rayp, geom, tstar, iiso, .false.)
            CALL FRSTAR(dist,depths,depthr,
     1          MNAME,3,tTSH,
     1          SSA, SSC, SSF, SSL, SSN, SSR,
     1          RRA, RRC, RRF, RRL, RRN, RRR,
     1          rayp, geom, tstar, iiso, .false.)
c-----
c       Now that all spectra for a given distance are in
c       make time history and output in desired order
c-----
        np2 = nptin/2 + 1
        df = 1.0/(nptin*dt)
c-----
c       perform inverse FFT
c       If the number of data points is too large, then use supfft
c-----
c----
c       This is a two pass process
c       Starting with  ur uz tz durdz and other derivatives on LUN=10
c       compute ur uz ut and strains an then output them in LUN=9
c       
c-----
         do 1300 jk=1,15
                rewind 10
                do k=1,np2
                    read(10) (xout(i),i=1,NGRN)
                    datc(k)=xout(jk)
                    if(k.gt.1)then
                        datc(nptin+2-k)=conjg(datc(k))
                    endif
                enddo
c-----
c       perform inverse FFT
c-----
c       also take the opportunity to perform the correct
c       spherical spreading for an earth with radius of 6371.0 km
c-----
                    call zfour(datc,nptin,+1,dt,df)
                    fac = exp(alpha*tfirst)
                    dfac = exp(alpha*dt)
c-----
c     undamp 
c-----
                    do 1303 k=1,nptin
                        x(k) = real(datc(k))*fac
                        fac = fac*dfac
 1303               continue

                    npts = nptin
c-----
c           convert velocity to displacement or acceleration
c           if required. Since we are using surface waves, the
c           DC component must be zero!
c
c           note the special case for the pressure field
c-----
                    if(dostep)then
                       if(idva.eq.0)then
                           call integ(x,npts,dt,.false.)
                       else if(idva.eq.2)then
                           call deriv(x,npts,dt)
                       endif
                    else
                       if(idva.eq.1)then
                           call deriv(x,npts,dt)
                       else if(idva.eq.2)then
                           call deriv(x,npts,dt)
                           call deriv(x,npts,dt)
                       endif
                    endif
c-----
c     These orientation are in the Sac convention, e.g., Z is positive up
c     Even though these are ZRT, the AZ is set to 0
c-----
                    if(iszrt(jk).eq.1)then
c-----
c               vertical
c-----
                        cmpinc = 0.0
                        cmpaz  =   0.0
c                     make vertical positive up
c             
                       do i=1,npts
                          x(i) = - x(i)
                       enddo
         
                    else if(iszrt(jk).eq.4)then
c-----
c               radial
c-----
                        cmpinc = 90.0
                        cmpaz  =   0.0
                    else if(iszrt(jk).eq.5)then
c-----
c               transverse
c-----
                        cmpinc = 90.0
                        cmpaz  =  90.0
                    endif
                stcomp = 'SYN     '
    
                cmpdt  = dt
                ksyear = 1970
                ksdoy  = 1
                ksmon  = 1
                ksday  = j
                kshour = 0
                ksmin  = 0
                ssec =  0.0

                o = 0
                distdg = -12345.
c-----
c
c-----
                stcomp = ocmpnm(jk)(1:3)
c-----
c     create the output file name
c-----
                 if(outfmt.eq.1)then
c                   DDDDDd_HHHh_ZZZz.cmp
                    ievdep = 10.*depths
                    istel  = 10.*depthr
                    idist  = 10.0*dist
                    write(ofile,11)idist,ievdep,istel,stcomp
   11               format(i6.6,'_',i4.4,'_',i4.4,'.',a3)
                 else if (outfmt.eq.2)then
c                   DDDDDddd_HHHhhh_ZZZzzz.cmp
                    ievdep = 1000.*depths
                    istel  = 1000.*depthr
                    idist  = 1000.*dist
                    write(ofile,12)idist,ievdep,istel,stcomp
   12               format(i7.7,'_',i6.6,'_',i6.6,'.',a3)
                 else if (outfmt.eq.3)then
c                   DDDDDdHHHh.grn
                    ievdep = 10.*depths
                    idist  = 10.0*dist
                    write(ofile,13)idist,ievdep,stcomp
   13               format(i6.6,i4.4,'.',a3)
                 else if (outfmt.eq.4)then
c                   DDDDdHHHh.grn
                    ievdep = 10.*depths
                    idist  = 10.0*dist
                    write(ofile,14)idist,ievdep,stcomp
   14               format(i5.5,i4.4,'.',a3)
                 else if (outfmt.eq.5)then
c                   DDDdddHhhh.grn
                    ievdep = 1000.*depths
                    idist  = 1000.*dist
                    write(ofile,15)idist,ievdep,stcomp
   15               format(i6.6,i4.4,'.',a3)
                 endif
                 ls = lgstr(ofile)
                call scmxmn(x,npts,depmax,depmin,depmen,indmax,indmin)
                call newhdr()
                call setfhv('DELTA   ',cmpdt,nerr)
                    bb = tfirst
                    ee = bb + (npts -1 ) * cmpdt
                    call setfhv('A       ',tp    ,nerr)
                    call setkhv('KA      ','P       ',nerr)
                    call setfhv('T0      ',tsv   ,nerr)
                    call setfhv('T1      ',tsh   ,nerr)
                    if(iiso.eq.0)then
                       call setkhv('KT0     ','S       ',nerr)
                       call setkhv('KT1     ','S       ',nerr)
                    else
                       call setkhv('KT0     ','SV      ',nerr)
                       call setkhv('KT1     ','SH      ',nerr)
                    endif
                call setfhv('B       ',bb    ,nerr)
                call setfhv('E       ',ee    ,nerr)
                call setfhv('O       ',o     ,nerr)
                call setkhv('KO      ','O       ',nerr)
                call setfhv('STLA    ',stlat ,nerr)
                call setfhv('STLO    ',stlon ,nerr)
                call setfhv('TIMMAX  ',bb+ indmax*cmpdt,nerr)
                call setfhv('TIMMIN  ',bb+ indmin*cmpdt,nerr)
                call setfhv('DEPMIN  ',depmin,nerr)
                call setfhv('DEPMAX  ',depmax,nerr)
c-----
c               our elevation is KM, SAC is set to KILOMETERS
c-----
                call setfhv('STEL    ',depthr,nerr)
                call setfhv('EVLA    ',evlat ,nerr)
                call setfhv('EVLO    ',evlon ,nerr)
c-----
c               our depth is KM, SAC is set to KILOMETERS
c-----
                call setfhv('EVDP    ',depths ,nerr)
                call setfhv('DIST    ',dist,nerr)
                evstaz =   0.0
                stevaz = 180.0
                call setfhv('AZ      ',evstaz,nerr)
                call setfhv('BAZ     ',stevaz,nerr)
                call setfhv('GCARC   ',distdg,nerr)
                call setfhv('CMPAZ   ',cmpaz ,nerr)
                call setfhv('CMPINC  ',cmpinc,nerr)
                call setfhv('DEPMEN  ',depmen,nerr)
                call setfhv('STLA    ',-12345.,nerr)
                call setfhv('STLO    ',-12345.,nerr)
                call setfhv('EVLA    ',-12345.,nerr)
                call setfhv('EVLO    ',-12345.,nerr)
c-----
c       set integer header value
c-----
                call setnhv('NZYEAR  ',ksyear,nerr)
                call setnhv('NZJDAY  ',ksdoy ,nerr)
                call setnhv('NZHOUR  ',kshour,nerr)
                call setnhv('NZMIN   ',ksmin ,nerr)
                isec = ssec
                smsec = (ssec - isec) * 1000.0
                call setnhv('NZSEC   ',isec  ,nerr)
                jsmsec = smsec
                call setnhv('NZMSEC  ',jsmsec,nerr)
                call setnhv('NPTS    ',npts, nerr)
c-----
c       set logical header value
c-----
                call setlhv('LEVEN   ',.true.,nerr)
                call setlhv('LPSPOL  ',.true.,nerr)
                call setlhv('LOVROK  ',.true.,nerr)
                call setlhv('LCALDA  ',.false.,nerr)
c-----
c       set character header value
c-----
                call setkhv('KSTNM    ','SYN     ' ,nerr)
                if(iobsyn.eq.2)then
                    call setkhv('KEVNM   ','SYNTHETI',nerr)
                    call setkhv('KEVNMC  ','C       ',nerr)
                else
                    call setkhv('KEVNM   ','OBSERVED',nerr)
                    call setkhv('KEVNMC  ','        ',nerr)
                endif
                call setkhv('KCMPNM   ',stcomp ,nerr)
c-----
c       set enumerated header value
c-----
                call setihv('IFTYPE   ','ITIME   ',nerr)
                call setihv('IZTYPE   ','IB      ',nerr)
                    call bwsac(12,npts,ofile(1:ls),x)
 1300   continue
        close(9,status='delete')
        return
        end

        subroutine FRSTAR(r,hs,hr,mname,ipsvsh,time,
     1          SSA, SSC, SSF, SSL, SSN, SSR,
     1          RRA, RRC, RRF, RRL, RRN, RRR,
     1          rayp, geom, tstar, iiso, dolock)
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
c       SSA     R   A at the source
c       SSC     R   C at the source
c       SSF     R   F at the source
c       SSL     R   L at the source
c       SSN     R   N at the source
c       SSR     R - density at the source
c       RRA     R   A at the receiver
c       RRC     R   C at the receiver
c       RRF     R   F at the receiver
c       RRL     R   L at the receiver
c       RRN     R   N at the receiver
c       RRR     R - density at the receiver
c       rayp R   Ray parameter in sec/km
c       geom R   geometrical spreading factor
c       tstar R   geometrical spreading factor
c       dolock L .true. apply locked mode which means to ignore the
c                bottom layer
c-----
        implicit none
        real r, hs, hr, time
        real SSA, SSC, SSF, SSL, SSN, SSR
        real RRA, RRC, RRF, RRL, RRN, RRR
        real rayp, geom, tstar
        logical dolock
        character mname*(*)
        integer ipsvsh
        logical ext
c-----
c-----
c       internal variables
c-----
        real depths, depthr
        real dphs, dphr, dphref
        integer lmaxs, lmaxr, lmaxref

        integer NL
        parameter (NL=200)
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL),
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax

        common/depref/refdep
        real refdep

        integer l, lgstr
        integer iunit, iiso, iflsph, idimen, icnvel, ierr
        character title*80 
        
c-----
c       get the earth model
c-----
        inquire(file=mname,exist=ext)
        if(.not. ext)call usage('Model file does not exist')
        l = lgstr(mname)

                call getmod(1,mname,mmax,title,iunit,iiso,iflsph,
     1          idimen,icnvel,ierr,.false.)
        if(ierr .lt. 0)return      
                call tdomod()
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

        RRA = TA(lmaxr)
        RRC = TC(lmaxr)
        RRF = TF(lmaxr)
        RRL = TL(lmaxr)
        RRN = TN(lmaxr)
        RRR = TRho(lmaxr)
        SSA = TA(lmaxs)
        SSC = TC(lmaxs)
        SSF = TF(lmaxs)
        SSL = TL(lmaxs)
        SSN = TN(lmaxs)
        SSR = TRho(lmaxs)

c-----
c       compute the travel time
c-----
        call fstarr(r,time,lmaxs, lmaxr, lmaxref,
     1      hs, hr, ipsvsh,iflsph, rayp,
     2      tstar, dolock)
        return
        end

        subroutine fstarr(dist,tfirst,lmaxs,lmaxr,lmaxref,
     1      depths,depthr,ipsvsh,iflsph, rayp,
     2      tstar, dolock)
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
c       depths  R   - depth of source
c       depthr  R   - depth of receiver
c       ipsvsh  I*4 1 - get P time
c               2 - get SV time
c               3 - get SH time
c               4 - get pP time
c               5 - get sP time
c       iflsph  I*4 - 0 flat earth model
c                 1 spherical earth model
c       rayp    R   - ray parameter in sec/km
c       geom R   geometrical spreading factor
c       dolock L .true. apply locked mode which means to ignore the
c                bottom layer
c-----
c
c       18 JAN 2008 - everything is straightforward. The addition of
c          the request for pP and sP changes the logic in that
c          the direct arrival is ignored, and that the upgoing refraction 
c          from the source is ignored. We handle this by just setting
c          a very large tfirst before trying to do the modified 
c          downward path refraction to avoid another level of
c          if/then/else/endif
c-----
        real dist, tfirst, depths, depthr
        real rayp
        integer lmaxs, lmaxr, lmaxref, ipsvsh, iflsph
        logical dolock

        integer NL
        parameter (NL=200)
        common/timodel/TD(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL),
     3      frefp(NL), frefs(NL)
        real TD, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmmax
        integer mmmax
        common/depref/refdep
        real refdep
        common/shwave/TLsh(NL), TNsh(NL), TRhosh(NL),qbsh(NL)
        real TLsh, TNsh, TRhosh,qbsh

        integer mmax

        real*4  h(NL)

        real*8   pupper
        complex*16 p
        integer lmx, lmn
        integer i, l
        real sumx, sumt, tt
        real time

        COMPLEX*16 A,B,C,D,E,F
        COMPLEX*16 wvno2, omega2, wvn, omg
        COMPLEX*16 rsh, rp, rsv
c-----
c       get the modified eigen vectors x11 and x31 
c           have lambda1 (rp ) factored out
c               modified eigen vectors x22 and x42 
c           have lambda2 (rsv) factored out
c-----
        COMPLEX*16 X11, X21, X31, X41
        COMPLEX*16 X12, X22, X32, X42
        COMPLEX*16 NP, NSV
        real getvel

        COMPLEX*16 dtdp
        complex*16 pold, pcur, dp
        complex*16 dtdpold, dtdpcur
        integer ilast

        logical baseisp, layerisp

c-----
c       initialize
c-----
        omg = dcmplx(1.0d+00, 0.0d+00)
        omega2 = omg *omg
        tstar = 0.0


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
            call tdosph()
        endif
c-----
c       now fill in velocity array according to desired first arrival
c       for SH there can be no water layer
c       for SV can be a water layer
c       Also define the Q for the T* analysis. Note we define
c        eventually q = 1/Q based on whether the given Q > or < 1
c-----
        do i=1,mmax
            if(qa(i) .gt. 1.0)then
                qa(i) = 1.0 / qa(i)
            endif
            if(qb(i) .gt. 1.0)then
                qb(i) = 1.0 / qb(i)
            endif
            h(i) = td(i)
        enddo

c-----
c       For the computations we look at four cases
c       1) direct path between source and receiver 
c       2) refracted arrivals       
c          a) path is downward from source and then up to
c             receiver
c          b) path is upward from the source and then down to
c             receiver
c          This recognized the possibility that velocity does
c          not increase uniformly with depth
c-----
                    
c-----
c       direct arrival 
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
c          reflection occurs when dt/dp = 0, so search for the p value 
c          numerically. The travel time is just
c               t = p r + Sum eta h
c-----
            ps = 1.0/getvel(TA,TL,TN,TRho,lmaxs,ipsvsh)
            pr = 1.0/getvel(TA,TL,TN,TRho,lmaxr,ipsvsh)
            if(ps.lt.pr)then
                pupper = ps
            else
                pupper = pr
            endif
            do 1000 l=lmn,lmx
                vl = getvel(TA,TL,TN,TRho,l,ipsvsh)
                if(vl.eq.0.0)return
                p = dcmplx(1.0/vl, 0.0d+00)
                if(dreal(p).lt.pupper)pupper = p
 1000       continue
            pold = dcmplx(0.0d+00, 0.0d+00)
            dp =  dcmplx(pupper/100.0, 0.0d+00)
            do  i=0,100
                ilast = i
                p = i*dp
                wvn = p * omg
                wvno2 = wvn * wvn
            call gttdtdp(omega2,wvno2,omg,wvn,p,ipsvsh,
     1          dist,lmn,lmx-1,time,dtdp,tstar)
            if(dreal(dtdp).lt. 0.0d+00)then
c----
c                       refine the root
c-----
                        pold = p - dp
                        pcur = p
                        dtdpcur = dtdp
                ilast = i -1
                        go to 2000
                endif
                dtdpold = dtdp
        enddo

c-----
c       assume we always get here
c-----
 2000   continue
c-----
c       use interval halving until I can compute the d2t/dp2!
c       also as a fallback, do not do this if the maximum index about was 100
c-----
        if(ilast.ne.100)then
        do 3000 i=1,20
            p = 0.5*(pold + pcur)
            wvn = p*omg
            wvno2 = wvn*wvn
            call gttdtdp(omega2,wvno2,omg,wvn,p,ipsvsh,
     1          dist,lmn,lmx-1,time,dtdp,tstar)
            if(dsign(1.0d+00,dreal(dtdpcur)).eq.
     1          dsign(1.0d+00,dreal(dtdp)))then
                pcur = p
                dtdpcur = dtdp
            else
                pold = pcur
                dtdpold = dtdp
            endif

 3000       continue
        endif
            tfirst = time
            rayp = dreal(p)
c-----
c       now proceed through the possible refracted arrivals
c       considering first upward rays from the source
c-----  
        if(lmn.gt.1)then
        do 3020 m=1,lmn-1
c-----
c       m is the refracting layer
c
c       get velocity of refractor testing if fluid
c-----
            if(ipsvsh.eq.1 .or. ipsvsh.eq.4 .or.ipsvsh.eq.5)then
                 vel = getvel(TA,TL,TN,TRho,m,1)
                 if(vel.eq.0.0)go to 3040
                 baseisp = .true.
            else  
                 vel = getvel(TA,TL,TN,TRho,m,ipsvsh)
                 if(vel.eq.0.0.and.ipsvsh.eq.2)then
c-----
c                     permit S-P conversion, e.g., SKS
c-----
                      vel = getvel(TA,TL,TN,TRho,m,1)
                      baseisp = .true.
                 else
                      baseisp = .false.
                 endif
            endif
            if(vel.eq.0.0)goto 3040
            p = 1.0/vel
            wvn = p * omg
            wvno2 = wvn * wvn
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
        if(ipsvsh.le.3)then
             vlmn = getvel(TA,TL,TN,TRho,lmn,ipsvsh)
             vlmx = getvel(TA,TL,TN,TRho,lmx,ipsvsh)
        else
c-----
c       this is just a subterfuge since we will not use the results
c-----
             vlmn = getvel(TA,TL,TN,TRho,lmn,1)
             vlmx = getvel(TA,TL,TN,TRho,lmx,1)
        endif
        if(vlmx.ge.vel)go to 3020
        if(vlmn.ge.vel)go to 3020
c-----
c       single leg
c-----
            sumx = 0.0
            sumt = 0.0
            ts = 0.0
            do 3021 l=lmn,lmx-1
                if(ipsvsh.eq.1 .or. ipsvsh.eq.4 .or.ipsvsh.eq.5)then
                     vl = getvel(TA,TL,TN,TRho,l,1)
                     layerisp = .true.
                else  
                     vl = getvel(TA,TL,TN,TRho,l,ipsvsh)
                     layerisp = .false.
                     if(vl.eq.0.0.and.ipsvsh.eq.2)then
c-----
c                         permit S-P conversion, e.g., SKS
c-----
                          vl = getvel(TA,TL,TN,TRho,l,1)
                          layerisp = .true.
                     endif
                endif

                if(vl.gt.vel)go to 3020
                if(vl.eq.0.0)go to 3040
                call tgetegn(a,b,c,d,e,f,omega2,wvno2,rp, rsv, rsh,
     1              x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,l,omg,wvn)
                if(ipsvsh.eq.1.or.ipsvsh.eq.4.or.ipsvsh.eq.5)then
C                    if(dimag(rp).lt. 0.0d+00)then
C                         rp = - rp
C                    endif
                    sumt = sumt + h(l) * rp /dcmplx(0.0d+00,1.0d+00)
                    sumx = sumx + h(l) * 
     1                  p/(rp/dcmplx(0.0d+00, 1.0d+00))
                    ts = ts + qa(l)*h(l)*(p*p - rp*rp)
     1                   /(rp/dcmplx(0.0d+00, 1.0d+00))
                else if(ipsvsh.eq.2)then
C                    if(dimag(rsv).lt. 0.0d+00)then
C                         rsv = - rsv
C                    endif
                     if(layerisp)then
                    sumt = sumt + h(l) * rp /dcmplx(0.0d+00,1.0d+00)
                    sumx = sumx + h(l) * 
     1                  p/(rp/dcmplx(0.0d+00, 1.0d+00))
                    ts = ts + qa(l)*h(l)*(p*p - rp*rp)
     1                   /(rp/dcmplx(0.0d+00, 1.0d+00))
                      else
                    sumt = sumt + h(l) * rsv/dcmplx(0.0d+00,1.0d+00)
                    sumx = sumx + h(l) * 
     1                  p/(rsv/dcmplx(0.0d+00, 1.0d+00))
                    ts = ts + qb(l)*h(l)*(p*p - rsv*rsv)
     1                   /(rsv/dcmplx(0.0d+00, 1.0d+00))
                      endif
                else if(ipsvsh.eq.3)then
C                    if(dimag(rsh).lt. 0.0d+00)then
C                         rsh = - rsh
C                    endif
                    sumt = sumt + h(l) * rsh/dcmplx(0.0d+00,1.0d+00)
                    sumx = sumx + h(l) * 
     1                  p/(rsh/dcmplx(0.0d+00, 1.0d+00))
                    ts = ts + qb(l)*h(l)*(p*p - rsh*rsh)
     1                   /(rsh/dcmplx(0.0d+00, 1.0d+00))
                endif
 3021       continue
            do 3022 l=m+1,lmn-1
                if(ipsvsh.eq.1 .or. ipsvsh.eq.4 .or.ipsvsh.eq.5)then
                     vl = getvel(TA,TL,TN,TRho,l,1)
                     layerisp = .true.
                else  
                     vl = getvel(TA,TL,TN,TRho,l,ipsvsh)
                     layerisp = .false.
                     if(vl.eq.0.0.and.ipsvsh.eq.2)then
c-----
c                         permit S-P conversion, e.g., SKS
c-----
                          vl = getvel(TA,TL,TN,TRho,l,1)
                          layerisp = .true.
                     endif
                endif
                if(vl.gt.vel)go to 3020
                if(vl.eq.0.0)go to 3040
                call tgetegn(a,b,c,d,e,f,omega2,wvno2,rp, rsv, rsh,
     1              x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,l,omg,wvn)
                if(ipsvsh.eq.1.or.ipsvsh.eq.4.or.ipsvsh.eq.5)then
C                    if(dimag(rp).lt. 0.0d+00)then
C                         rp = - rp
C                    endif
                    sumt = sumt + 2.0*h(l) * rp /dcmplx(0.0d+00,1.0d+00)
                    sumx = sumx + 2.0*h(l) * 
     1                  p/(rp/dcmplx(0.0d+00, 1.0d+00))
                    ts = ts + 2.*qa(l)*h(l)*(p*p - rp*rp)
     1                   /(rp/dcmplx(0.0d+00, 1.0d+00))
                else if(ipsvsh.eq.2)then
C                    if(dimag(rsv).lt. 0.0d+00)then
C                         rsv = - rsv
C                    endif
                     if(layerisp)then
                    sumt = sumt + 2.0*h(l) * rp /dcmplx(0.0d+00,1.0d+00)
                    sumx = sumx + 2.0*h(l) * 
     1                  p/(rp/dcmplx(0.0d+00, 1.0d+00))
                    ts = ts + 2.*qa(l)*h(l)*(p*p - rp*rp)
     1                   /(rp/dcmplx(0.0d+00, 1.0d+00))
                      else
                    sumt = sumt + 2.0*h(l) * rsv/dcmplx(0.0d+00,1.0d+00)
                    sumx = sumx + 2.0*h(l) * 
     1                  p/(rsv/dcmplx(0.0d+00, 1.0d+00))
                    ts = ts + 2.*qb(l)*h(l)*(p*p - rsv*rsv)
     1                   /(rsv/dcmplx(0.0d+00, 1.0d+00))
                      endif
                else if(ipsvsh.eq.3)then
C                    if(dimag(rsh).lt. 0.0d+00)then
C                         rsh = - rsh
C                    endif
                    sumt = sumt + 2.0*h(l) * rsh/dcmplx(0.0d+00,1.0d+00)
                    sumx = sumx + 2.0*h(l) * 
     1                  p/(rsh/dcmplx(0.0d+00, 1.0d+00))
                    ts = ts + 2.*qb(l)*h(l)*(p*p - rsh*rsh)
     1                   /(rsh/dcmplx(0.0d+00, 1.0d+00))
                endif
 3022       continue
            tint = sumt
            tt = tint + dist / vel
            if(baseisp)then
                 ts = ts + qa(m)*(dist-sumx)/vel
            else
                 ts = ts + qb(m)*(dist-sumx)/vel
            endif
            if(tt .lt. tfirst .and. dist.ge.sumx)then
                  tfirst = tt
                  rayp = dreal(p)
                 tstar = ts
            endif
 3020       continue
 3040       continue
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
c       get velocity of refractor testing if fluid
c-----
            if(ipsvsh.eq.1 .or. ipsvsh.eq.4 .or.ipsvsh.eq.5)then
                 vel = getvel(TA,TL,TN,TRho,m,1)
                 if(vel.eq.0.0)go to 2040
                 baseisp = .true.
            else  
                 vel = getvel(TA,TL,TN,TRho,m,ipsvsh)
                 if(vel.eq.0.0.and.ipsvsh.eq.2)then
c-----
c                     permit S-P conversion, e.g., SKS
c-----
                      vel = getvel(TA,TL,TN,TRho,m,1)
                      baseisp = .true.
                 else
                      baseisp = .false.
                 endif
            endif
            if(vel.eq.0.0)go to 2040
            p = 1.0/vel
            wvn = p * omg
            wvno2 = wvn * wvn
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
c       refraction velocity otherwise there will be no real ray
c-----
        if(ipsvsh.le.3)then
             vlmn = getvel(TA,TL,TN,TRho,lmn,ipsvsh)
             vlmx = getvel(TA,TL,TN,TRho,lmx,ipsvsh)
        else
c-----
c            this is just a subterfuge since 
c            we will not use the results for pP sP
c-----
             vlmn = getvel(TA,TL,TN,TRho,lmn,1)
             vlmx = getvel(TA,TL,TN,TRho,lmx,1)
        endif
        if(vlmx.ge.vel)go to 2020
        if(vlmn.ge.vel)go to 2020
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
                      vp = getvel(TA,TL,TN,TRho,l,1)
                      if(vp.gt.vel)go to 2020
                      if(vp.eq.0.0)go to 2040
                call tgetegn(a,b,c,d,e,f,omega2,wvno2,rp, rsv, rsh,
     1              x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,l,omg,wvn)
                    sumt = sumt + 2.*h(l) * rp /dcmplx(0.0d+00, 1.0d+00)
                    sumx = sumx + 2.*h(l) *
     1                  p/(rp/dcmplx(0.0d+00, 1.0d+00))
                    ts = ts + qa(l)*2.*h(l)*(p*p - rp*rp)
     1                   /(rp/dcmplx(0.0d+00, 1.0d+00))
                  enddo
            else if(ipsvsh.eq.5)then
c-----
c               sP
c-----
                  do  l=lmaxref,lmaxs - 1
                      vp = getvel(TA,TL,TN,TRho,l,1)
                      vs = getvel(TA,TL,TN,TRho,l,2)
                      if(vp.gt.vel)go to 2020
                      if(vp.eq.0.0)go to 2040
                      if(vs.gt.vel)go to 2020
                      if(vs.eq.0.0)go to 2040
                call tgetegn(a,b,c,d,e,f,omega2,wvno2,rp, rsv, rsh,
     1              x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,l,omg,wvn)
                    sumt = sumt + h(l) * rp /dcmplx(0.0d+00, 1.0d+00)
     1                          + h(l) * rsv/dcmplx(0.0d+00, 1.0d+00)
                    sumx = sumx 
     1                  + h(l)*p/(rp/dcmplx(0.0d+00, 1.0d+00))
     1                  + h(l)*p/(rsv/dcmplx(0.0d+00, 1.0d+00))
                    ts = ts + qa(l)*h(l)*(p*p - rp*rp)
     1                   /(rp/dcmplx(0.0d+00, 1.0d+00))
     1                      + qb(l)*h(l)*(p*p - rsv*rsv)
     1                   /(rsv/dcmplx(0.0d+00, 1.0d+00))
                  enddo
            endif
            do 2021 l=lmn,lmx - 1
                if(ipsvsh.eq.1 .or. ipsvsh.eq.4 .or.ipsvsh.eq.5)then
                     vl = getvel(TA,TL,TN,TRho,l,1)
                     layerisp = .true.
                else  
                     vl = getvel(TA,TL,TN,TRho,l,ipsvsh)
                     layerisp = .false.
                     if(vl.eq.0.0.and.ipsvsh.eq.2)then
c-----
c                         permit S-P conversion, e.g., SKS
c-----
                          vl = getvel(TA,TL,TN,TRho,l,1)
                          layerisp = .true.
                     endif
                endif

                if(vl.gt.vel)go to 2020
                if(vl.eq.0.0)go to 2040
                call tgetegn(a,b,c,d,e,f,omega2,wvno2,rp, rsv, rsh,
     1              x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,l,omg,wvn)
c-----
c      KLUDGE - to fix the case when the imaginary part of the wavenumber is
c            negative - it must be positive
c-----
                if(ipsvsh.eq.1.or.ipsvsh.eq.4.or.ipsvsh.eq.5)then
C                    if(dimag(rp).lt. 0.0d+00)then
C                         rp = - rp
C                    endif
                    sumt = sumt + h(l) * rp /dcmplx(0.0d+00, 1.0d+00)
                    sumx = sumx + h(l) * 
     1                  p/(rp/dcmplx(0.0d+00, 1.0d+00))
                    ts = ts + qa(l)*h(l)*(p*p - rp*rp)
     1                   /(rp/dcmplx(0.0d+00, 1.0d+00))
                else if(ipsvsh.eq.2)then
C                    if(dimag(rsv).lt. 0.0d+00)then
C                         rsv = - rsv
C                    endif
                     if(layerisp)then
                    sumt = sumt + h(l) * rp /dcmplx(0.0d+00,1.0d+00)
                    sumx = sumx + h(l) * 
     1                  p/(rp/dcmplx(0.0d+00, 1.0d+00))
                    ts = ts + qa(l)*h(l)*(p*p - rp*rp)
     1                   /(rp/dcmplx(0.0d+00, 1.0d+00))
                      else
                    sumt = sumt + h(l) * rsv/dcmplx(0.0d+00,1.0d+00)
                    sumx = sumx + h(l) * 
     1                  p/(rsv/dcmplx(0.0d+00, 1.0d+00))
                    ts = ts + qb(l)*h(l)*(p*p - rsv*rsv)
     1                   /(rsv/dcmplx(0.0d+00, 1.0d+00))
                      endif
                else if(ipsvsh.eq.3)then
C                    if(dimag(rsh).lt. 0.0d+00)then
C                         rsh = - rsh
C                    endif
                    sumt = sumt + h(l) * rsh/dcmplx(0.0d+00, 1.0d+00)
                    sumx = sumx + h(l) * 
     1                  p/(rsh/dcmplx(0.0d+00, 1.0d+00))
                    ts = ts + qb(l)*h(l)*(p*p - rsh*rsh)
     1                   /(rsh/dcmplx(0.0d+00, 1.0d+00))
                endif
 2021       continue
c-----
c       double leg
c-----

            do 2022 l=lmx,m-1
                if(ipsvsh.eq.1 .or. ipsvsh.eq.4 .or.ipsvsh.eq.5)then
                     vl = getvel(TA,TL,TN,TRho,l,1)
                     layerisp = .true.
                else  
                     vl = getvel(TA,TL,TN,TRho,l,ipsvsh)
                     layerisp = .false.
                     if(vl.eq.0.0.and.ipsvsh.eq.2)then
c-----
c                         permit S-P conversion, e.g., SKS
c-----
                          vl = getvel(TA,TL,TN,TRho,l,1)
                          layerisp = .true.
                     endif
                endif

                if(vl.gt.vel)go to 2020
                if(vl.eq.0.0)go to 2040
                call tgetegn(a,b,c,d,e,f,omega2,wvno2,rp, rsv, rsh,
     1              x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,l,omg,wvn)
                if(ipsvsh.eq.1.or.ipsvsh.eq.4.or.ipsvsh.eq.5)then
C                    if(dimag(rp).lt. 0.0d+00)then
C                         rp = - rp
C                    endif
                    sumt = sumt + 2.0*h(l) * rp /dcmplx(0.0d+00,1.0d+00)
                    sumx = sumx + 2.0*h(l) * 
     1                  p/(rp/dcmplx(0.0d+00, 1.0d+00))
                    ts = ts + 2.*qa(l)*h(l)*(p*p - rp*rp)
     1                   /(rp/dcmplx(0.0d+00, 1.0d+00))
                else if(ipsvsh.eq.2)then
C                    if(dimag(rsv).lt. 0.0d+00)then
C                         rsv = - rsv
C                    endif
                     if(layerisp)then
                    sumt = sumt + 2.0*h(l) * rp /dcmplx(0.0d+00,1.0d+00)
                    sumx = sumx + 2.0*h(l) * 
     1                  p/(rp/dcmplx(0.0d+00, 1.0d+00))
                    ts = ts + 2.*qa(l)*h(l)*(p*p - rp*rp)
     1                   /(rp/dcmplx(0.0d+00, 1.0d+00))
                      else
                    sumt = sumt + 2.0*h(l) * rsv/dcmplx(0.0d+00,1.0d+00)
                    sumx = sumx + 2.0*h(l) * 
     1                  p/(rsv/dcmplx(0.0d+00, 1.0d+00))
                    ts = ts + 2.*qb(l)*h(l)*(p*p - rsv*rsv)
     1                   /(rsv/dcmplx(0.0d+00, 1.0d+00))
                      endif
                else if(ipsvsh.eq.3)then
C                    if(dimag(rsh).lt. 0.0d+00)then
C                         rsh = - rsh
C                    endif
                    sumt = sumt + 2.0*h(l) * rsh/dcmplx(0.0d+00,1.0d+00)
                    sumx = sumx + 2.0*h(l) * 
     1                  p/(rsh/dcmplx(0.0d+00, 1.0d+00))
                    ts = ts + 2.*qb(l)*h(l)*(p*p - rsh*rsh)
     1                   /(rsh/dcmplx(0.0d+00, 1.0d+00))
                endif
 2022       continue
            tint = sumt
            tt = tint + dist / vel
            if(baseisp)then
                 ts = ts + qa(m)*(dist-sumx)/vel
            else
                 ts = ts + qb(m)*(dist-sumx)/vel
            endif
            if(tt .lt. tfirst .and. dist.ge.sumx)then
                 tfirst = tt
                 rayp = dreal(p)
                 tstar = ts
            endif

                 vp = getvel(TA,TL,TN,TRho,m,1)
                 vsv = getvel(TA,TL,TN,TRho,m,2)
                 vsh = getvel(TA,TL,TN,TRho,m,3)
 2020       continue
 2040       continue
             if(tfirst .eq. 1.0e+30)then
                tfirst = -12345.
                tstar  = -12345.
                rayp   = -12345.
             endif
        return
        end

        subroutine tdosph()
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
c       We will treat all as P-SV for the heck of it
c       This requires more work
c-----
c       mmax    I*4 number of layers
c       TA     R   A 
c       TC     R   C 
c       TF     R   F 
c       TL     R   L 
c       TN     R   N 
c                  note  density not required
c       TD     R   layer thickness
c       v() R   array of velocities
c       h() R   array of layer thicknesses
c       ipsvsh  I       1 - get P time
c                       2 - get SV time
c                       3 - get SH time
c       refdep R   Reference depth for the model specification
c
c       Note we need the constants here.  Since the velocities
c       must increase with depth, e.g., vf = vs (a/r)
c       and that density  varies
c       as rhof = rhos (a/r)^-P, [not the TI surface wave code has not yet
c        been written], then using the model that m = rho beta^2, we have
c
c       TA = rho VA^2,
c       TAf = rhof * VAf^2 = rhos (a/r)^-P VAs^2 (a/r)^2
c           = (a/r)^2-P TAs
c-----
        integer NL
        parameter (NL=200)
        common/timodel/TD(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL),
     3      frefp(NL), frefs(NL)
        real TD, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/depref/refdep
        real refdep
        common/shwave/TLsh(NL), TNsh(NL), TRhosh(NL), qbsh(NL)
        real TLsh, TNsh, TRhosh, qbsh


        double precision z0,z1,r0,r1,ar,tmp

        common/earth/radius
        real radius

        ar=radius
        r0=ar + refdep
        td(mmax)=1.0
        do 10 i=1,mmax
            r1=r0-dble(td(i))
            z0=ar*dlog(ar/r0)
            z1=ar*dlog(ar/r1)
            td(i)=z1-z0
c-----
c        attempt 7 15 2007 - use standard rule but at mid layer depth as per DGH
c-----
            TMP=(ar+ar)/(r0+r1)
c-----
c                SV
c-----
                 rhosph    = trho(i)
                 trho(i)   = rhosph * tmp**(-2.275)
                 trhosh(i) = rhosph * tmp**(-5)

                 ta(i)=ta(i)*tmp**(-0.2750)
                 tc(i)=tc(i)*tmp**(-0.2750)
                 tf(i)=tf(i)*tmp**(-0.2750)

                 elsph = tl(i)
                 tl(i)  =elsph*tmp**(-0.2750)
                 tlsh(i)=elsph*tmp**(-3.0)
                 ensph = tn(i)

                 tn(i)=ensph*tmp**(-0.2750)
                 tnsh(i)=ensph*tmp**(-3.0)
            r0 = r1
   10   continue
        td(mmax)=0.0
        return
        end

        subroutine tdomod()
c-----
c       just fill in the TRhosh, TLsh, TNsh and qbsh arrays
c-----
        integer NL
        parameter (NL=200)
        common/timodel/TD(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL),
     3      frefp(NL), frefs(NL)
        real TD, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/depref/refdep
        real refdep
        common/shwave/TLsh(NL), TNsh(NL), TRhosh(NL), qbsh(NL)
        real TLsh, TNsh, TRhosh, qbsh

        do i=1,mmax
           TLsh(i) = TL(i)
           TNsh(i) = TN(i)
           TRhosh(i) = TRho(i)
           qbsh(i) = qbsh(i)
        enddo
        return
        end

        subroutine srclyr(depth,lmax,dph)
        parameter (LER=0, LIN=5, LOT=6)
        parameter(NL=200)
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
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

        subroutine getabc(m,omg,wvn,a,b,c,d,e,f)
        implicit none
        integer m
        COMPLEX*16 omg,wvn
        COMPLEX*16 a, b, c, d, e, f
        integer NL
        parameter (NL=200)
        common/timodel/H(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real H, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        a = wvn * TF(m) / TC(m)
        b = 1.0/TC(m)
        c = - TRho(m)*omg*omg + wvn*wvn *(TA(m) -TF(m)*TF(m)/TC(m))
        d = - wvn
        e = 1.0/TL(m)
        f = - TRho(m)*omg*omg
        return
        end                                               

        subroutine tgetegn(a,b,c,d,e,f,omega2,wvno2,rp, rsv, rsh,
     1      x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,m,omg,wvn)
        implicit none
        COMPLEX*16 A,B,C,D,E,F
        COMPLEX*16 wvno2, omega2
        COMPLEX*16 rsh, rp, rsv
c-----
c       get the modified eigen vectors x11 and x31 have 
c           lambda1 (rp ) factored out
c               modified eigen vectors x22 and x42 have 
c           lambda2 (rsv) factored out
c-----
        COMPLEX*16 X11, X21, X31, X41
        COMPLEX*16 X12, X22, X32, X42
c-----
c       norms
c-----
        COMPLEX*16 NP, NSV
        integer m
        COMPLEX*16 omg, wvn
        COMPLEX*16 xka2, xkb2

        integer NL
        parameter (NL=200)
        common/timodel/H(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real H, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/shwave/TLsh(NL), TNsh(NL), TRhosh(NL), qbsh(NL)
        real TLsh, TNsh, TRhosh, qbsh
c-----
c       internal variables
c-----
        COMPLEX*16 L2(2)
        COMPLEX*16 bb, cc
        COMPLEX*16 CDSQRT 

        COMPLEX*16 ZFAC
c-----
c       first test to see if a fluid layer - if it is fluid, the
c       eigenfunctions are specially computed and we need only the
c       rp
c-----
        if(TL(m).eq.0.0 .or. TN(m).eq.0.0)then
            rp = cdsqrt(wvno2 -omega2*TRho(m)/TA(m))
            rsv = dcmplx(0.0d+000, 0.0d+00)
            rsh = dcmplx(0.0d+000, 0.0d+00)
            return
        endif
        
        call getabc(m,omg,wvn,a,b,c,d,e,f)
c-----
c       Do the SH
c-----
        rsh = CDSQRT(TNsh(m)*wvno2/TLsh(m) - Trhosh(m)*omega2/TLsh(m)) 
        if( dimag(rsh) .lt. 0.0)then
                rsh = - rsh
        endif
c-----
c       Do the P and SV
c-----
c-----
c       The characteristic equation to be solved is
c
c       L^4 + L^2[ -2 ad -ec -fb ] + [ (d^2+ef)(a^2+bc)] = 0
c-----
        bb = -2.0d+00 * a*d - e*c -f*b
        cc = ( d*d + e*f)*(a*a + b*c)
        L2(1) = ( - bb + CDSQRT(bb*bb - 4.000*cc))/2.0d+00
        L2(2) = ( - bb - CDSQRT(bb*bb - 4.000*cc))/2.0d+00

        L2(1) = cc/L2(2)
c-----
c       Use the Lambda^2 values to form
c       xka^2 == k^2 - L(1)^2
c       xkb^2 == k^2 - L(2)^2
c       Associate the smallest xka, xkb with the P!
c-----
        xka2 = wvno2 - L2(1)
        xkb2 = wvno2 - L2(2)
        if(cdabs(xkb2) .lt. cdabs(xka2))THEN
                ZFAC = L2(1)
                L2(1) = L2(2)
                L2(2) = ZFAC
        endif
        rp  = CDSQRT(L2(1))
        rsv = CDSQRT(L2(2))
        if( dimag(rp) .lt. 0.0)then
                rp = - rp
        endif
        if( dimag(rsv) .lt. 0.0)then
                rsv = - rsv
        endif
c-----
c       get the norms - note that the true norm will be 
c           2  NP amd 2 L(2) NSV
c       The factorization permits us to use the sin nz/n or n sin nz
c-----
        NP  = (  L2(1)*(-2*a*b*d + 2*a*a*e + b*c*e - b*b*f)
     1      + (a*a+b*c)*(2*b*d*d - 2*a*d*e + b*e*f - c*e*e) )
        NSV = (- L2(2)*(2*b*d*d - 2*a*d*e - c*e*e + b*e*f)
     1      + (d*d+e*f)*(2*a*b*d - 2*a*a*e + b*b*f - b*c*e) )
c-----
c       get the modified eigen vectors x11 and x31 have 
c           lambda1 (rp ) factored out
c               modified eigen vectors x22 and x42 have 
c           lambda2 (rsv) factored out
c-----
        x11 =              (b*d - a*e)
        x21 =  b*L2(1) - e*(b*c + a*a)
        x31 =    L2(1) -   (a*d + c*e)
        x41 = -a*L2(1) + d*(b*c + a*a)

        x12 = -e*L2(2) + b*(d*d + e*f)
        x22 = ( b*d - a*e)
        x32 = d*L2(2) - a*(d*d + e*f)
        x42 = - ( L2(2) -  a*d - b*f)
c-----
c       TEST
c       Force the eigenfunctions to be as given in 5.4.4
c-----
        zfac = rp / x21
        x11  = x11 *zfac
        x21  = x21 *zfac
        x31  = x31 *zfac
        x41  = x41 *zfac

        zfac = rsv / x12
        x12  = rsv
        x22  = x22 * zfac
        x32  = x32 * zfac
        x42  = x42 * zfac
        
        np   = x11*x41 - x21*x31
        nsv  = x12*x42 - x22*x32

        return
        end

        subroutine gttdtdp(omega2,wvno2,omg,wvn,p,ipsvsh,r,
     1      llow,lhgh,time,dtdp,tstar)
        integer ipsvsh, llow, lhgh
        real r, time, tstar
        COMPLEX*16 p
        COMPLEX*16 A,B,C,D,E,F
        COMPLEX*16 wvno2, omega2, wvn, omg
        COMPLEX*16 rsh, rp, rsv
c-----
c       omega2  C - angular frequency squared
c       wvno2   C - wavenumber squared
c       omg     C - angular frequency
c       wvn     C - wavenumber
c       p       C - ray parameter
c       ipsvsh  I - 1 P, 2 SV, 3 SH, 4 pP, 5 sP
c                  since this is for the direct arrival pP and sP not considered
c       r       C - distance
c       llow    I - layer interface indices
c       lhgh    I - layer interface indices
c       time    R - travel time
c       dtdp    C - This must be zero for the direct arrival
c       tstar   R - attenuation operator
c-----
c       get the modified eigen vectors x11 and x31 have 
c           lambda1 (rp ) factored out
c               modified eigen vectors x22 and x42 have 
c           lambda2 (rsv) factored out
c-----
        COMPLEX*16 X11, X21, X31, X41
        COMPLEX*16 X12, X22, X32, X42
        COMPLEX*16 NP, NSV

        COMPLEX*16 detadp
        COMPLEX*16 dtdp

        integer NL
        parameter (NL=200)
        common/timodel/TD(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real TD, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax

        dtdp = dcmplx(dble(r), 0.0d+00)
        time = p*r
        ts = 0.0
        do 1000 l=llow,lhgh

            call tgetegn(a,b,c,d,e,f,omega2,wvno2,rp, rsv, rsh,
     1          x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,l,omg,wvn)     
        if(ipsvsh.eq.1 .or. ipsvsh.eq.4 .or. ipsvsh.eq.5)then
C               if(dimag(rp).lt. 0.0d+00)then
C                    rp = - rp
C               endif
        dtdp  = dtdp + 
     1      TD(l)*detadp(p,rp *x11,x21,rp *x31,x41,NP ,l,wvn,omg,rp )
        time = time + rp *TD(l)/dcmplx(0.0d+00, 1.0d+00)
C       write(6,*)'l,deta:',l,TD(l),detadp(p,rp *x11,x21,rp *x31,x41,NP ,l,wvn,omg,rp ),dtdp
                    ts = ts + qa(l)*TD(l)*(p*p - rp*rp)
     1                   /(rp/dcmplx(0.0d+00, 1.0d+00))

        else if(ipsvsh.eq.2)then
C               if(dimag(rsv).lt. 0.0d+00)then
C                    rsv = - rsv
C               endif
        dtdp = dtdp + TD(l)*detadp(p,x12,rsv*x22,x32,
     1      rsv*x42,NSV,l,wvn,omg,rsv)
        time = time + rsv*TD(l)/dcmplx(0.0d+00, 1.0d+00)
                    ts = ts + qb(l)*TD(l)*(p*p - rsv*rsv)
     1                   /(rsv/dcmplx(0.0d+00, 1.0d+00))
        else if(ipsvsh.eq.3)then
C               if(dimag(rsh).lt. 0.0d+00)then
C                    rsh = - rsh
C               endif
        dtdp = dtdp + TD(l)*((omg*wvn*TN(l)/TL(l))/rsh)/
     1      cmplx(0.0d+00, 1.0d+00)
        time = time + rsh*TD(l)/dcmplx(0.0d+00, 1.0d+00)
                    ts = ts + qb(l)*TD(l)*(p*p - rsh*rsh)
     1                   /(rsh/dcmplx(0.0d+00, 1.0d+00))
        endif
 1000   continue
        tstar = ts
        return
        end

        function detadp(p,x1,x2,x3,x4,NORM,m,wvn,omg,nu)
c-----
c       van der Hijden  6.109
c-----
        implicit none
        complex*16 p,x1,x2,x3,x4,NORM,nu
        complex*16 detadp
        integer m
        complex*16 wvn, omg
        integer NL
        parameter(NL=200)
        common/timodel/H(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real H, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        complex*16 da, db, dc, dd, de, df
        da =   omg *TF(m)/TC(m)
        db =   dcmplx(0.0d+00, 0.0d+00)
        dc =   2.0d+00*wvn*omg*(TA(m) - TF(m)*TF(m)/TC(m)) 
        dd = - omg
        de =   dcmplx(0.0d+00, 0.0d+00)
        df =   dcmplx(0.0d+00, 0.0d+00) 
        detadp =  x4 * (       dd*x2         + de*x4  )
     1      - x3 * (da*x1        + db*x3          )
     1      - x2 * (       df*x2         - dd*x4  )
     1      + x1 * (dc*x1        - da*x3          )
        detadp = detadp /(2.0d+00 * nu * NORM )
        detadp = detadp/dcmplx(00d+00, 1.0d+00)
        return
        end

        function getvel(TA,TL,TN,TRho,m,ipsvsh)
c-----
c     this determines the horizontally propagating velocity
c     This is useful for the refraction and for determining the
c     limits on a reflected arrival
c-----
            real TA(*), TL(*), TN(*), TRho(*)
            integer m, ipsvsh
            real getvel
        
            if(ipsvsh.eq.1)then
                getvel = sqrt(TA(m)/TRho(m))
            else if(ipsvsh.eq.2)then
                getvel = sqrt(TL(m)/TRho(m))
            else if(ipsvsh.eq.3)then
                getvel = sqrt(TN(m)/TRho(m))
            else
                getvel = 1.0
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
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/shwave/TLsh(NL), TNsh(NL), TRhosh(NL), qbsh(NL)
        real TLsh, TNsh, TRhosh, qbsh

        integer m
        real dep, dp, dphh, hsave
        integer ls
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
c       Do not create unnecessary layers, e.g., at 
c            surface and internally
c       However do put in a zero thickness layer at the 
c            base if necessary
c-----
        if(dphh .eq. 0.0 .and. ls.ne.mmax)then
            return
        else
c-----
c           adjust layering
c-----
             do 102 m = mmax,ls,-1
                d(m+1) = d(m)
                TA(m+1) = TA(m)
                TC(m+1) = TC(m)
                TF(m+1) = TF(m)
                TL(m+1) = TL(m)
                TN(m+1) = TN(m)
                TRho(m+1) = TRho(m)
                TLsh(m+1) = TLsh(m)
                TNsh(m+1) = TNsh(m)
                TRhosh(m+1) = TRhosh(m)
                qa(m+1) = qa(m)
                qb(m+1) = qb(m)
                qbsh(m+1) = qbsh(m)
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
