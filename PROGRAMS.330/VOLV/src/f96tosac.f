        program f96tosac
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME V                                                       c
c                                                                     c
c      PROGRAM: f96tosac                                              c
c                                                                     c
c      COPYRIGHT 1996 R. B. Herrmann                                  c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c------------------------------------- -------------------------------c
c-----
c       11 SEP 2000 add TP, TSV TSH to FILE96 format
c       14 JAN 2001 add A, C, F, L, N and density information
c       30 JAN 2001 error in setlhv('LOVROK  - was 'LOVERK
c       27 MAR 2002 DEFINE SAC DEPTH IN KILOMETERS NOT METERS
c       07 AUG 2002 NEW option for naming the output sac files
c               HHHHhDDDDd.grn where where the digits
c               represent the source depths of HHHH.h km
c               and epicentral distance of DDDD.d km and
c               grn is ZDS ZEX etc
c       14 JUN 2005 option of naming sac files HHHHHhDDDDDd.grn
c       18 MAR 2007 corrected error in setting LPSPOL - previously
c              spelling error in call - I used LSPOL
c       04 FEB 2013 changed file naming with -E flag to
C              permit file name to indicate nearest meter
c       11 FEB 2021 introduce more file name types
c-----
c       Convert file96(V) format to SAC
c           If the file is SYNTHETIC seismogram
c           the SAC file names will be of the form
c               XXXyyyP.sac
c           where
c               XXX will be 001 002 ..., e.g.,
c                    the position in the file
c               yyy will be ZDD ... PEX, a Green s
c                   function name
c               P wil be Z R or T
c
c           Else, the name will be of the form
c
c               XXXyyyP.sac
c           where
c               XXX will be 001 002 ..., e.g.,
c                    the position in the file
c               P wil be Z R T N or E
c               yyy will be up to three characters
c                   if the Station Name
c-----
c-----
c       variables in common blocks
c
c       iftype  I*4 File type
c               1 - single trace
c               3 - three component
c               16 - Green's function
c               21 - Green's function
c
c       iobsyn  I*4 1 - observed
c               2 - synthetic
c
c       itmfrq  I*4 1 - time series
c               2 - Fourier spectra (not implemented)
c
c       iunit   I*4 1   - counts 
c               2   - cm
c               3   - cm/sec
c               4   - cm/sec/sec
c               5   - m
c               6   - m/sec
c               7   - m/sec/sec
c               8   - microns
c               9   - microns/sec
c               10  - microns/sec/sec
c       junit   I*4
c               11  - Pa  (nt/m^2)
c               12  - MPa  (mega Pascal)
c
c       cfilt   C*80    comment on filtering operations
c
c       keyear  I*4 event year
c       kemon   I*4 event mon
c       keday   I*4 event day
c       kehour  I*4 event hour
c       kemin   I*4 event minute
c       esec    R*4 event second
c       evlat   R*4 event latitude
c       evlon   R*4 event longitude
c       evdep   R*4 event depth
c       
c       stname  C*8 station name
c               NOTE: AH(6), SEED(5), CSS(6), SAC(8)
c       
c       stlat   R*4 station latitude
c       stlon   R*4 station longitude
c       stelev  R*4 station elevation
c
c
c       distkm  R*4 epicentral distance in km
c       distdg  R*4 epicentral distance in degrees (111.195 km/degree)
c       evstaz  R*4 event -> station azimuth, degrees east of north
c       stevaz  R*4 station -> event azimuth, degrees east of north
c       
c       cpulse  C*80    pulse description
c
c       ccomnt  C*80    comment
c
c       jsrc    I*4 Array of indices to indicate if traces are 
c                   present
c               iftype =  1  jsrc(1) indicates if trace is 
c                   present
c               iftype =  3  jsrc(i),i=1,5 indicates if trace 
c                   is present
c                   i = Z, 2=N, 3=E, 4=R, 5=T
c                   (we could carry all five, but the real 
c                   purpose  is to permit no more than 
c                   3 traces, but 2 may all that are 
c                   present in a real data set)
c
c-----

        common/ihdr96/ ftype, iftype, iobsyn, itmfrq, iunit, idecon,
     1      keyear, kemon, keday, kehour, kemin,
     2      ksyear, ksmon, ksday, kshour, ksmin,
     3      jsrc, junit
        integer*4 ftype, iftype, iobsyn, itmfrq, iunit, idecon
        integer*4 keyear, kemon, keday, kehour, kemin
        integer*4 ksyear, ksmon, ksday, kshour, ksmin
        integer*4 jsrc(21), junit

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
        integer*4 npts
        real*4 ssec
c-----
c       command line arguments
c-----
        character fname*80
        logical verby
        logical sacbin
        integer outfmt
c-----
c       internal program variables
c-----
        parameter (LIN=5, LOT=6, LER=0)
        parameter (NSAMP=16384)
        real*4 x(NSAMP)
        logical ext
        real*8 eepoch, sepoch
        character*50 str
        integer*4 doy, date
        character ofile*80
        integer idist, ievdep, istel
        
c-----
c       machine dependent initialization
c-----
        call mchdep()
c-----
c       parse command line arguments
c-----
        call gcmdln(fname,verby,sacbin,outfmt)
c-----
c       determine whether the file comes from the standard input
c       or from the specified file
c-----
        if(fname.ne. ' ' )then
            inquire(file=fname,exist=ext)
            lnm = lgstr(fname)
            if(.not. ext)then
                write(LER,*)'Input file: ',fname(1:lnm),
     1              '  does not exist'
                call usage()
            endif
            lunin = 1
            open(lunin,file=fname,status='unknown',form=
     1          'formatted',access='sequential')
            rewind lunin
        else
            lunin = LIN
        endif
c-----
c       process the input
c-----
        ntrace = 0
 1000   continue
            call rdhd96(lunin,nerr)
            if(nerr.lt.0)go to 9000
            ntrace = ntrace + 1
            ndone = 0
            do 1100 jtr = 1, 21
                if(jsrc(jtr).ne.0)then
                    call rdtr96(lunin,stcomp,cmpinc,cmpaz,
     1                  cmpdt, npts, ksyear, ksmon, 
     2                  ksday, kshour, ksmin, ssec, 
     3                  x,nerr,NSAMP)
                    if(nerr .lt. 0 )go to 9999
                    ndone = ndone + 1
                    if(ndone.gt.iftype)go to 201
c-----
c       create a SAC file for this trace
c-----
c-----
c       set time
c-----
        if(ksyear.eq.0)ksyear = 1970
        if(keyear.eq.0)keyear = 1970
        if(ksmon.eq.0)ksmon = 1
        if(kemon.eq.0)kemon = 1
        if(ksday.eq.0)ksday = 1
        if(keday.eq.0)keday = 1
        if(iobsyn.eq.2)then
            tsec = ssec
            ssec = 0.0
        endif
        call htoe(ksyear,ksmon,ksday,kshour,ksmin,ssec,sepoch)
        call htoe(keyear,kemon,keday,kehour,kemin,esec,eepoch)
        call etoh(sepoch,date,str,doy,
     1      ksyear,ksmon,ksday,kshour,ksmin,ssec)
        o = eepoch - sepoch
c-----
c       set floating point header value
c-----
                call scmxmn(x,npts,depmax,depmin,depmen,indmax,indmin)
                call newhdr()
                call setfhv('DELTA   ',cmpdt,nerr)
                if(iobsyn.eq.2)then
                    b = tsec
                    e = b + (npts -1 ) * cmpdt
                    call setfhv('A       ',tp    ,nerr)
                    call setkhv('KA      ','P       ',nerr)
                    if(stcomp(1:1).eq.'Z' .or.
     1                  stcomp(1:1).eq.'R')then
                    call setfhv('T0      ',tsv   ,nerr)
                    call setkhv('KT0     ','SV      ',nerr)
                    else if(stcomp(1:1).eq.'T')then
                    call setkhv('KT1     ','SH      ',nerr)
                    call setfhv('T1      ',tsh   ,nerr)
                    endif
                else
                    b = 0.0
                    e = b + (npts -1 ) * cmpdt
                    call setfhv('A       ',tp    ,nerr)
                    call setfhv('T0      ',tsv   ,nerr)
                    call setfhv('T1      ',tsh   ,nerr)
                    call setkhv('KA      ','P       ',nerr)
                    call setkhv('KT0     ','SV      ',nerr)
                    call setkhv('KT1     ','SH      ',nerr)
                endif
                call setfhv('B       ',b     ,nerr)
                call setfhv('E       ',e     ,nerr)
                call setfhv('O       ',o     ,nerr)
                call setkhv('KO      ','O       ',nerr)
                call setfhv('STLA    ',stlat ,nerr)
                call setfhv('STLO    ',stlon ,nerr)
                call setfhv('TIMMAX  ',b + indmax*cmpdt,nerr)
                call setfhv('TIMMIN  ',b + indmin*cmpdt,nerr)
                call setfhv('DEPMIN  ',depmin,nerr)
                call setfhv('DEPMAX  ',depmax,nerr)
c-----
c               our elevation is KM, SAC is set to KILOMETERS
c-----
                call setfhv('STEL    ',stelev,nerr)
                call setfhv('EVLA    ',evlat ,nerr)
                call setfhv('EVLO    ',evlon ,nerr)
c-----
c               our depth is KM, SAC is set to KILOMETERS
c-----
                call setfhv('EVDP    ',evdep ,nerr)
                call setfhv('DIST    ',distkm,nerr)
                call setfhv('AZ      ',evstaz,nerr)
                call setfhv('BAZ     ',stevaz,nerr)
                call setfhv('GCARC   ',distdg,nerr)
                call setfhv('CMPAZ   ',cmpaz ,nerr)
c----
c       convert from SEED to SAC component incidence
c-----
                cmpinc = cmpinc + 90.0
                call setfhv('CMPINC  ',cmpinc,nerr)
                call setfhv('DEPMEN  ',depmen,nerr)
c-----
c       set integer header value
c-----
                call setnhv('NZYEAR  ',ksyear,nerr)
                call setnhv('NZJDAY  ',doy   ,nerr)
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
                call setkhv('KSTNM    ',stname ,nerr)
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
c-----
c       create a unique file name
c-----
c       outfmt format   string
c       0       10   
c       1       11   
c       2       12   
c       3       13   
c       4       14   
c       5       15   
c-----
c       AIIIjjjGRN.sac
c       BIIIjjjGRN.sac
c-----
    9   format(a1,i3,i2,a3,'.asc')
   10   format(a1,i3,i2,a3,'.sac')

c-----
c       DDDDDddd_HHHhhh_ZZZzzz.cmp
c-----
   11   format(i6.6,'_',i4.4,'_',i4.4,'.',a3)

c-----
c       DDDDDd_HHHh_ZZZz.cmp
c-----
   12   format(i7.7,'_',i6.6,'_',i6.6,'.',a3)

c-----
c       DDDDDdHHHh.grn
c-----
   13   format(i6.6,i4.4,'.',a3)

c-----
c       DDDDdHHHh.grn
c-----
   14   format(i5.5,i4.4,'.',a3)

c-----
c       DDDdddHhhh.grn
c-----
   15   format(i6.6,i4.4,'.',a3)

                if(iobsyn.eq.2)then
                    if(sacbin)then
                      if(outfmt.eq.0)then
                       write(ofile,10)'B',
     1                   ntrace,jtr,stcomp(1:3)
                      else if(outfmt.eq.1)then
c                         DDDDDd_HHHh_ZZZz.cmp
                          ievdep = 10.*evdep
                          istel = 10.*stelev 
                          ievdep = 10.0*evdep
                          idist  = 10.0*distkm
                          write(ofile,11)idist,ievdep,
     1                       istel,stcomp(1:3)
                      else if(outfmt.eq.2)then
c                         DDDDDddd_HHHhhh_ZZZzzz.cmp
                          ievdep = 1000.0*evdep
                          istel  = 1000.0*stelev 
                          idist  = 1000.0*distkm
                          write(ofile,12)idist,ievdep,
     1                       istel,stcomp(1:3)
                      else if(outfmt.eq.3)then
c                         DDDDDdHHHh.grn
                          ievdep = 10.0*evdep
                          idist  = 10.0*distkm
                          write(ofile,13)idist,ievdep,
     1                       stcomp(1:3)
                      else if(outfmt.eq.4)then
c                        DDDdddHhhh.grn
                         ievdep = 10.0*evdep
                         idist = 10.0*distkm
                         write(ofile,14)idist,ievdep,
     1                       stcomp(1:3)
                      else if(outfmt.eq.5)then
c                        DDDdddHhhh.grn
                         ievdep = 1000.*evdep
                         idist = 1000.0*distkm
                         write(ofile,15)idist,ievdep,
     1                       stcomp(1:3)
                      endif
                    else
                      write(ofile,9)'A',
     1                  ntrace,jtr,stcomp(1:3)
                    endif
                else
                    if(sacbin)then
                    write(ofile,10)'B',
     1                  stname(1:4),stcomp(1:3)
                    else
                    write(ofile,9)'A',
     1                  stname(1:4),stcomp(1:3)
                    endif
                endif
c-----
c       zero blanks
c-----
                do 1101 i=1,12
                    if(ofile(i:i).eq.' ')ofile(i:i)='0'
 1101           continue
                    if(sacbin)then
                        call bwsac(3,npts,ofile,x)
                    else
                        call awsac(3,npts,ofile,x)
                    endif
                endif
  201           continue
 1100       continue
        go to 1000
c-----
c       close open files
c-----
 9000   continue
 9999   continue
        if(lunin.ne.LIN)close (lunin)
        end

        subroutine gcmdln(fname,verby,sacbin,outfmt)
c-----
c       parse command line arguments
c
c       requires subroutine mgtarg() and function mnmarg()
c
c-----
c       fname   C*80 - file name for conversion
c       verby   L    - verbose output on standard error
c       sacbin  L    - .true. if SAC binary output, else ascii
c       outfmt  I    - 
c        -A          SAC alphanumeric file, else binary
c        -B          SAC binary (default)
c        -G          Output names in form DDDDdHHHh.grn  
c        -T          Output names in form DDDDDdHHHh.grn 
c        -E          Output names in form DDDdddHhhh.grn 
c         The format for the name of the binary output attempts to
c         give information on epicentral distance (km), 
c         source depth (km), and receiver depth(km). The options are
c        -FMT 1      DDDDDd_HHHh_ZZZz.cmp
c                    e.g. 005001_1234_0045.Uz
c        -FMT 2      DDDDDddd_HHHhhh_ZZZzzz.cmp
c                   e.g. 00500123_123456_004578.Erf
c        -FMT 3      DDDDDdHHHh.grn  (same as -T)
c                    e.g. 0050010041.ZVF
c        -FMT 4      DDDDdHHHh.grn   (same as -G)
c                    e.g. 050010045.Srz
c        -FMT 5      DDDdddHhhh.grn  (same as -E)
c                    e.g. 5001234578.Err
c         where D is for epicentral distance, H source depth, and
c         Z receiver depth. The lower case indicates the digits 
c         to the right of the decimal place.  The examples above
c         are for an epicentral distance is 500.123 km, source
c         depth 123.456 km and receiver depth 4.578 km.
c         The grn are ZDD,RDD,ZDS,RDS,TDS,ZSS,RSS,TSS,ZEX,REX,
c         ZVF,RVF,ZHF,RHF,THF,PEX,PDD,PDS,PSS,PVF,PHF

c           where DDDDd=12345 -> distance of 1234.5 km
c           where  HHHh=1234  -> source depth of 123.4 km
c           where GRN is ZDD RDD ZDS RDS TDS ZSS RSS TSS ZEX REX
c               ZVF RVF ZHF RHF THF PEX 
c                   - if 2 DISTANCE FIELD is one digit larger to permit
c               > 9999.9 epicentral distance
c                   - if 3 Use exploration format of form DDDdddHhhh
c               which gives 1 merter resolution
c                    - if 0 use file BRECRgrGRN.sac
c-----
        character*25 name
        character fname*(*)
        logical verby
        logical sacbin
        integer outfmt

        verby = .false.
        fname = ' '
        sacbin = .true.
        outfmt = 4
        nmarg=mnmarg()
        i = 0
   11   i = i + 1
        if(i.gt.nmarg)go to 13
            call mgtarg(i,name)
            if(name(1:2).eq.'-?' .or. name(1:2).eq.'-h')then
                call usage()
C           else if(name(1:2).eq.'-v')then
C               verby = .true.
            else if(name(1:2).eq.'-A')then
                sacbin = .false.
            else if(name(1:2).eq.'-B')then
                sacbin = .true.
                outfmt = 0
            else if(name(1:2).eq.'-G')then
                outfmt = 4
                sacbin = .true.
            else if(name(1:2).eq.'-T')then
                outfmt = 3
                sacbin = .true.
            else if(name(1:2).eq.'-E')then
                outfmt = 5
                sacbin = .true.
            else if(name(1:4).eq.'-FMT')then
                i=i+1
                call mgtarg(i,name)
                read(name,'(bn,i10.0)')outfmt
                if(outfmt.lt.0 .or. outfmt.gt.5)then
                   outfmt = 4
                endif
            else
                fname = name
            endif
        goto 11
   13   continue
        return
        end
        
        subroutine usage()
        parameter (LER=0, LIN=5, LOT=6)
        write(LER,*)'f96tosac  [-A | -B | -G | -T  | -E | -FMT ifmt ]',
     1    ' [ -? | -h ] file_name'
        write(LER,*)
     1  ' file_name   Name of file96 file to be converted'
        write(LER,*)
     1  '             If not given, input is stdin'
        write(LER,*)
     1  'OUTPUT FILE NAME'
        write(LER,*)
     1  ' -A           SAC alphanumeric file, else binary'
        write(LER,*)
     1  ' -B           SAC binary (default)'
        write(LER,*)
     1  ' -G (default) Output names in form DDDDdHHHh.grn ' 
     2  , '(binary)'
        write(LER,*)
     1  ' -T           Output names in form DDDDDdHHHh.grn' 
     2  , '(binary)'
        write(LER,*)
     1  ' -E           Output names in form DDDdddHhhh.grn' 
     2  , '(binary)'
        write(LER,*)
     1  ' The format for the name of the binary output attempts to'
        write(LER,*)
     1  ' give information on epicentral distance (km), '
        write(LER,*)
     1  ' source depth (km), and receiver depth(km). The options are'
        write(LER,*)
     1  ' -FMT 1      DDDDDd_HHHh_ZZZz.grn'
        write(LER,*)
     1  '             e.g. 005001_1234_0045.Uz'
        write(LER,*)
     1  ' -FMT 2      DDDDDddd_HHHhhh_ZZZzzz.grn'
        write(LER,*)
     1  '             e.g. 00500123_123456_004578.ZVF'
        write(LER,*)
     1  ' -FMT 3      DDDDDdHHHh.grn  (same as -T)'
        write(LER,*)
     1  '             e.g. 0050010041.ZVF'
        write(LER,*)
     1  ' -FMT 4      DDDDdHHHh.grn   (same as -G)'
        write(LER,*)
     1  '             e.g. 050010045.ZVF'
        write(LER,*)
     1  ' -FMT 5      DDDdddHhhh.grn  (same as -E)'
        write(LER,*)
     1  '             e.g. 5001234578.ZVF'
        write(LER,*)
     1  ' where D is for epicentral distance, H source depth, and'
        write(LER,*)
     1  ' Z receiver depth. The lower case indicates the digits '
        write(LER,*)
     1  ' to the right of the decimal place.  The examples above'
        write(LER,*)
     1  ' are for an epicentral distance is 500.123 km, source'
        write(LER,*)
     1  ' depth 123.456 km and receiver depth 4.578 km.'
        write(LER,*)
     1  ' The -FMT x naming convention is compatible with the'
        write(LER,*)
     1  ' mt command of gsac. Howver the model must be in units of'
        write(LER,*)
     1  ' km, km/s and g/cm^3'
        write(LER,*)
     1  ' -?          This online help'
        write(LER,*)
     1  ' -h          This online help'
        stop
        end

        

