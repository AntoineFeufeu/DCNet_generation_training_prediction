        program tpulse96
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME III                                                     c
c                                                                     c
c      PROGRAM: TPULSE96                                              c
c                                                                     c
c      COPYRIGHT 1996 ,2010                                           c
c      R. B. Herrmann                                                 c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
c       CHANGES
c       06 SEP 2000 - build in P, SV and SH first arrival times
c       14 JAN 2001 - use new FILE96.2 with sa, sc, sf, sl, sn, sr
c       19 FEB 2001 - fixed obscure error having that messed up
c               spectra in fundamental mode only at long period
c               by forcing kkf = 4 for all runs. Also
c               redefined mode =0 for fundamental in -M mode
c       30 JAN 2002 - caught error in excitr for ZDD had 
c            d1=w1*v1*fact1 should
c               be d1=d1*v1*fact1
c       17 OCT 2002 - Added description of dfile format to usage routine
c       10 SEP 2003 - error in pressure field in fluid due to a
c               point horizontal force
c           BAD  cvp(6) = cvp(6) + dkp(6)*xmom*cmplx(st2,-ct2)
c           GOOD cvp(6) = cvp(6) + dkp(6)*xmom*cmplx(st1,-ct1)
c           also changes notation slightly dk, dkk, dkp in excitr
c               to dkz, dkr , dkp
c       05 FEB 2004 - modified to be slightly more tolerant about DT for
c           user supplied pulse - now issues WARNING and not termination
c           mlarocca@ov.ingv.it
c       07 FEB 2005 - add a -Z flag to indicate that the 
c            internal parabolic 
c           or triangular pulses are to be zero phase 
c       04 AUG 2006 - corrected error in first arrival pick that
c               falsely gave the refraction time instead of the
c               direct time because the refraction arrival was
c               unphysical
c       12 DEC 2006 - set header values evlat, evlon, stlat, stlon to -12345
c           for compatibility of resulting SAC traces
c       26 MAY 2007 - careful reworking to make the slat2d work
c          FIXED exceed array bound problem with -2 flat. Everything
c          now defaults to NSAMP as the maximum of anything
c       21 JAN 2008 - upgraded travel time computation so that all code
c          is the same for CPS isotropic media
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
c       08 FEB 2008 -  subtle change in fstarr for source receiver in same layer -
c                      spherical mapping was not done
c       18 FEB 2009 -  caught egregious error in frstar where refdep was not set
c       23 JUL 2009 - ensured consistent call to FFT routine with
c               complex arguments. We now call zfour
c       17 MAR 2022  because of problems in working with a model consisting
c         of CGS velocities, densities and thicknesses, change instances of lines like\
c                 t1 = cdabs(cvt(j))
c                 if(t1 .le. 1.0d-20)cvt(j) = 0.0d+00
c         in subroutinex excitr and excitl to
c                 if(t1 .le. 1.0d-40)cvt(j) = 0.0d+00
c         This did not prolems with CGS units but make a pure MKS model bette
c-----
c       This program takes the output of tregn96 and tlegn96
c       and makes surface wave time series in the file96(V) format
c
c       PROGRAM CONTROL
c           INPUT  - from command line
c           OUTPUT - standard output
c-----
c       spulse96 -d Distance_File -v  -t -o -p -i -a alpha \
c           -l dur  -D -V -A -F rfile -m mult \
c       
c        -d Distance_File Distance control information
c        -v           Verbose output 
c        -t           Triangular pulse of base 2 dur dt 
c        -p           Parabolic Pulse of base  4 dur dt 
c        -o           Ohnaka pulse with parameter alpha 
c        -i           Dirac Delta function 
c        -a alpha     Shape parameter for Ohnaka pulse 
c        -D           Output is ground displacment 
c        -V           Output is ground velocity (default) 
c        -A           Output is ground acceleration 
c        -F rfile     User supplied pulse 
c        -m mult      Multiplier (default 1.0) 
c        -?           Write this help message 
c        -h           Write this help message 
c        -EX          Explosion and point force green s functions
c        -EQ          Earthquake and double couple green s functions
c        -ALL         Earthquake, Explosion and Point Force 
c            green s functions
c               
c-----
        parameter (LER=0, LIN=5, LOT=6)
        character*80 dfile, ryfile, lvfile, rfile
        integer*4 ntau, ipt, idva, iodva, ieqex
        real*4 xmult, alp
        logical dolat, dodble, dolock, dozero
        character ostr*80

        character mname*80
        integer ipar(20)
        real*4 fpar(20)
c-----
c       call machine dependent initialization
c-----
        call mchdep()
c-----
c       parse command line arguments
c-----
        call gcmdln(ntau,ipt,alp,dfile,idva,iodva,
     1      ieqex,xmult,rfile,ostr,dolat,dodble,nmode,dolock,
     1      dozero)
        if(dolat)then
            lvfile = 'tlatl96.egn'
            ryfile = 'tlatr96.egn'
        else
            lvfile = 'tlegn96.egn'
            ryfile = 'tregn96.egn'
        endif
c-----
c       process to make time series
c-----
        call process(ntau,ipt,alp,dfile,idva,iodva,
     1      ieqex,xmult,rfile,ryfile,lvfile,ostr,
     2      mname,ipar,fpar,dodble,nmode,dolock,dozero,dolat)
        end

        subroutine process(ntau,ipt,alp,dfile,idva,iodva,
     1      ieqex,xmult,rfile,ryfile,lvfile,ostr,
     2      mname,ipar,fpar,dodble,nmode,dolock,dozero,dolat)
c-----
c       make Green s functions for a given source pulse and
c       distance
c-----
        parameter (LER=0, LIN=5, LOT=6)
        character*80 dfile, ryfile, lvfile, rfile
        integer*4 ntau, ipt, idva, iodva, ieqex, nmode
        real*4 xmult, alp
        logical dolock, dozero,dolat

        real*4 dist, dt, tfirst
        integer npts
        logical ext
        character ostr*(*)

        character mname*80
        integer ipar(20)
        real*4 fpar(20)
        logical dodble
c-----
c       open file containing distance information
c-----
        inquire(file=dfile, exist=ext)
        if(.not. ext)then
            call usage('Distance Information File does not exist')
        endif
        open(3,file=dfile, status='old',form='formatted',
     1      access='sequential')
        rewind 3
c-----
c       read the data from the distance file, and process each distance
c       separately
c-----
 1000   continue
            read(3,*,err=9999,end=9999)dist,dt,npts,t0,vred
c-----
c       dist    R*4 - desired distance from the source
c               cannot be zero
c       dt  R*4 - desired sampling interval
c       npts    R*4 - desired number of points in the time series
c       t0,vred R*4 - control for the time of the first sample
c               if vred > 0, then the
c                   tfirst = t0 + dist/vred
c               else
c                   tfirst = t0
c-----
            dist = abs(dist)
            if(vred.gt.0.0)then
                tfirst = t0 + dist/vred
            else
                tfirst = t0
            endif
            call maksyn(ntau,ipt,alp,idva,iodva,ieqex,
     1          xmult,rfile,ryfile,lvfile,dist,dt,npts,
     2          tfirst,ostr,
     3          mname,ipar,fpar,dodble,nmode,dolock,dozero,dolat)
        go to 1000
 9999   continue
        close(3)
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
        read(4,'(i5,f10.3)')np,dtt
        read(4,10)(d(i),i=1,np)
   10   format(4e15.7)
        close(4)
        tau = np*dtt
        return
        end

        subroutine gcmdln(ntau,ipt,alp,dfile,idva,iodva,
     1      ieqex,xmult,rfile,ostr,dolat,dodble,nmode,
     2      dolock,dozero)
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
c       dfile   C*80- name of distance file
c       idva    I*4 - time history type
c               0 - displacement
c               1 - velocity
c               2 - acceleration
c       iodva   I*4 - force output time history type
c               Useful if internal pulse is not to reprsent a
c               step in fault displacement, explosion pressure,
c               or force
c               0 - displacement
c               1 - velocity
c               2 - acceleration
c       ieqex   I*4 -   Green s function choice
c               0 = earthquake and explosion Green s functions
c               1 = explosion and point force Green s functions
c               2 = earthquake, explosion and point force Green s func
c       xmult   R*4 - moment scaling factor
c       rfile   C*80- name of user provided pulse 
c               note dt must be the same as specified in dfile
c       dolat   L   .true. do laterally varying
c       nmode   I*4 - mode number to plot if >= 0
c               -  0 do fundamental only
c               - -1 do all
c               - -2 do fundamental only
c               - -3 do all higher
c       dolock  L   .true. locked mode used so do not use
c                   bottom layer for first arrival computation
c               .false.  (default)
c       dozero  L   .true. for triangular and parabolic pulses, center
c               at zero lag
c-----
        character*(*) dfile, rfile
        integer*4 ntau, ipt, idva, iodva, ieqex 
        real*4 xmult, alp
        logical  dolat, dodble, dolock, dozero
        character ostr*(*)

        integer*4 mnmarg
        character*80 name
c-----
c       initialization
c-----
        ntau = -1
        ipt = -1
        alp = -1.0
        idva = 1
        iodva = -1
        ieqex = 2
        rfile = ' '
        xmult = 1.0
        dfile = ' '
        ostr = ' '
        dolat = .false.
        dodble = .false.
        nmode = -1
        dolock = .false.
        dozero = .false.
        nmarg=mnmarg()
        i = 0
   11   i = i + 1
            if(i.gt.nmarg)go to 13
            call mgtarg(i,name)
            if(name(1:2).eq.'-l')then
                i=i+1
                call mgtarg(i,name)
                read(name,'(bn,i10)')ntau
            elseif(name(1:2).eq.'-d')then
                i=i+1
                call mgtarg(i,dfile)
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
            else if(name(1:2).eq.'-2')then
                dodble = .true.
            else if(name(1:2).eq.'-Z')then
                dozero = .true.
            else if(name(1:2).eq.'-F' .and. name(1:3).ne.'-FU')then
                ipt = 4
                        i=i+1
                        call mgtarg(i,rfile)
            else if(name(1:2).eq.'-D')then
                idva = 0
            else if(name(1:2).eq.'-V')then
                idva = 1
            else if(name(1:2).eq.'-A' .and. name(1:4).ne.'-ALL')then
                idva = 2
            else if(name(1:3).eq.'-OD')then
                iodva = 0
            else if(name(1:3).eq.'-OV')then
                iodva = 1
            else if(name(1:3).eq.'-OA')then
                iodva = 2
            else if(name(1:3).eq.'-EQ')then
                ieqex = 0
            else if(name(1:3).eq.'-EX')then
                ieqex = 1
            else if(name(1:4).eq.'-ALL')then
                ieqex = 2
            else if(name(1:2).eq.'-m')then
                i=i+1
                call mgtarg(i,name)
                read(name,'(bn,f20.0)')xmult
            else if(name(1:4).eq.'-LAT')then
                dolat = .true.
            else if(name(1:4).eq.'-LOC')then
                dolock = .true.
            else if(name(1:2).eq.'-M')then
                i=i+1
                call mgtarg(i,name)
                read(name,'(bn,i10)')nmode
                if(nmode.ge.0)nmode = nmode + 1
            else if(name(1:3).eq.'-FU')then
                nmode = -2 
            else if(name(1:3).eq.'-HI')then
                nmode = -3 
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
        if(dfile .eq. ' ')call usage('No distance control data file')
        if(ipt.eq.2 .and. alp.lt.0.0)
     1      call usage('No alpha for Ohnaka pulse')
        if(ipt.lt.0)
     1      call usage('No pulse shape defined')
        if(iodva.lt.0)iodva = idva
        lostr = len(ostr)
        ostr = 'spulse96'
        do 14 i=1,nmarg
            call getarg(i,name)
            L = lgstr(name)
            lo = lgstr(ostr)
            lomax = lo + 1 + L
            if(lomax .lt. lostr)then
                ostr = ostr(1:lo)//' '//name(1:L)
            endif
   14   continue
        return
        end

        subroutine usage(str)
        parameter (LER=0, LIN=5, LOT=6)
        character str*(*)
        write(LER,*)'spulse96:',str
        write(LER,*)'USAGE: ',
     1  'tpulse96 -d Distance_File [ -v  ] [ -t -o -p -i ] -a alpha',
     2  ' -l L [ -D -V -A]  [-F rfile ] [ -m mult] ',
     3  ' [ -OD -OV -OA ] [-FUND] [-HIGH] [-Z] ',
     4  ' [-EQEX -EXF -ALL] [-LAT] [-2] [ -M mode ] [-LOCK] [-?] [-h]'
        write(LER,*)
     1  ' -d Distance_File      Distance control file'
        write(LER,*)
     1  '   This contains one of more lines with following entries'
        write(LER,*)
     1  '       DIST(km) DT(sec) NPTS T0(sec) VRED(km/s)',
     2  '           first time point is T0 + DIST/VRED',
     3  '           VRED=0 means infinite velocity though'
        write(LER,*)
     1  ' -v           Verbose output'
        write(LER,*)
     1  ' -t           Triangular pulse of base 2 L dt'
        write(LER,*)
     1  ' -p           Parabolic Pulse of base  4 L dt'
        write(LER,*)
     1  ' -l L          (default 1 )duration control parameter'
        write(LER,*)
     1  ' -o           Ohnaka pulse with parameter alpha'
        write(LER,*)
     1  ' -i           Dirac Delta function'
        write(LER,*)
     1  ' -a alpha     Shape parameter for Ohnaka pulse'
        write(LER,*)
     1  ' -D           Output is ground displacment'
        write(LER,*)
     1  ' -V           Output is ground velocity (default)'
        write(LER,*)
     1  ' -A           Output is ground acceleration'
        write(LER,*)
     1  ' -F rfile     User supplied pulse'
        write(LER,*)
     1  ' -m mult      Multiplier (default 1.0)'
        write(LER,*)
     1  ' -OD           Output is ground displacement'
        write(LER,*)
     1  ' -OV           Output is ground velocity'
        write(LER,*)
     1  ' -OA           Output is ground acceleration'
        write(LER,*)
     1  ' -EXF         Explosion and point force green s functions'
        write(LER,*)
     1  ' -EQEX        Earthquake and double couple green s functions'
        write(LER,*)
     1  ' -ALL         Earthquake, Explosion and Point Force '
        write(LER,*)
     1  ' -LAT         (default false) Laterally varying eigenfunctions'
        write(LER,*)
     1  ' -2           (default false) Use double length  internally'
        write(LER,*)
     1  ' -M [ nmode ] (default all) mode number to ',
     2  'compute[0=fund,1=1st]'
        write(LER,*)
     1  ' -Z           (default false) zero phase ',
     2  'triangular/parabolic pulse'
        write(LER,*)
     1  ' -FUND        (default all) fundamental modes only  '
        write(LER,*)
     1  ' -HIGH        (default all) all higher modes only  '
        write(LER,*)
     1  ' -LOCK        (default false) locked mode used  '
        write(LER,*)
     1  ' -?           Write this help message'
        write(LER,*)
     1  ' -h           Write this help message'
        stop
        end

        subroutine maksyn(ntau,ipt,alp,idva,iodva,ieqex,
     1      xmult,rfile,ryfile,lvfile,dist,dt,npts,
     2      tfirst,ostr,
     2      mname,ipar,fpar,dodble,nmode,dolock,dozero,dolat)
        
c-----
c                                                                       
c     This program generates sixteen basic types of synthetic          
c     spectrum after combining the eigenfunctions from                 
c     sregn96 and/or slegn96  with the source spectrum.              
c                                                                    
c     The output is a time series in file96 format
c                                                                 
c     When number of points <= NSAMP, do FFT (call sub four) here,
c     otherwise, it will need 'supfft96' (through system call). 
c                                                              
c     The maximum dimension for time histories is LIMITED  by dimension
c     The maximum mode number at any period is NOT LIMITED.  
c                                                           
c     The Q is independent of frequency                   
c                                                        
c     The program uses a large number of files in order to function 
c     The file Logical Unit Numbers and their use are as follow:   
c                                                                 
c     1       Love Eigenfunction Source Depth                    
c     2       Rayleigh Eigenfunction Source Depth               
c     3       Distance control file                            
c     4       Temporary source pulse file used in pulud       
c     5       Standard Input                                 
c     6       Standard Output                               
c     7       Temporary Love Mode File                     
c     8       Temporary Rayleigh Mode File                
c     9       Internal Spectra File and then Output time series File 
c                                                                   
c     The program may seem complicated, but the computational structure
c     is designed for efficiency and relies upon the order of the output
c      of sdisp96, sregn96 and slegn96
c                                                                     
c     For each frequency, the 16 Green s functions are computed
c     at each distance.
c-----
        parameter (LER=0, LIN=5, LOT=6)
        character*80 ryfile, lvfile, rfile
        integer*4 ntau, ipt, idva,iodva, ieqex, npts
        real*4 xmult, alp, dist, dt, tfirst
        character ostr*(*)
        logical dodble, dolock, dozero, dolat
C-----
c       MSDOS, and SUN F77
c-----
        integer kerr, system
C-----
c       END MSDOS
c-----

        integer NSAMP2
        parameter(NSAMP2=16384)

        logical ext
        parameter (NGRN=21)
        integer*4 kkr,kkl,kkf,nper
        integer*4 np0
        integer*4 isign
        integer NSAMP
        parameter (NSAMP=16384)
          common/srctim/ src(NSAMP)
        parameter (NMD=100)
        common/srcspc/ssrc(NSAMP)
        complex ssrc
        common/rayl/
     1      ur(NMD,2),dur(NMD,2),uz(NMD,2),duz(NMD,2),
     2      ur0(NMD,2),uz0(NMD,2),tr0(NMD,2),tz0(NMD,2),
     3      wvmr(NMD,2),ares(NMD,2),arer(NMD,2),gamr(NMD,2),
     4      wvmrr(NMD,2), wvmrs(NMD,2)
        common/love/
     1      ut(NMD,2),dut(NMD,2),ut0(NMD,2),dut0(NMD,2),
     2      wvml(NMD,2),ales(NMD,2),aler(NMD,2),gaml(NMD,2),
     3      wvmlr(NMD,2), wvmls(NMD,2)
        common/ctrl/   sr,si,xmom,zeta,ms,np0,kkr,kkl,
     1                 kkf,itrigl,itrigr,mods,modt,modi
        common/resp/   df,bmax,modes1(2),modes2(2),pers1(2),pers2(2),
     1                 ipart1(2),ipart2(2)

        character mname*80
        integer ipar(20)
        real*4 fpar(20)

        data iopen0/0/
        itrigl = 0
        itrigr = 0
c-----
c       check for the existence of the Rayleigh and Love wave
c       eigenfunction files
c-----
        inquire(file=lvfile,exist=ext)
        if(ext .and. lvfile.ne.' ')then
            kkl=1
            open(1,file=lvfile,status='old',form='unformatted',
     1          access='sequential')
        else
            kkl=0
        endif
        inquire(file=ryfile,exist=ext)
        if(ext .and. ryfile.ne.' ')then
            kkr=1
            open(2,file=ryfile,status='old',form='unformatted',
     1          access='sequential')
        else
            kkr=0
        endif
c-----
c       enter source and other parameters:
c-----
c       Mode selection
c       kkf =   1   all modes
c           2   fundamental only
c           3   higher modes only
c           4   range
c                   mods = first
c                   modt = last
c                   modi = increment
c-----
        if(nmode.gt.0)then
            kkf = 4
            mods = nmode
            modt = nmode
            modi = 1
        else if(nmode.eq.-1)then
            kkf = 4
            mods = 1
            modt = 10000
            modi = 1
        else if(nmode.eq.-2)then
            kkf = 4
            mods = 1
            modt = 1
            modi = 1
        else if(nmode.eq.-3)then
            kkf = 4
            mods = 2
            modt = 10000
            modi = 1
        else
            kkf = 4
            mods = 1
            modt = 10000
            modi = 1
        endif
        xmom = xmult
c-----
c       define the source spectrum
c-----
        call srcpul(ntau,ipt,alp,dt,npts,rfile,duration)
        if(.not. dozero)then
            duration = 0.0
        endif
        zeta = 0.0
c-----
c       check for power of 2
c-----
        call npow2(npts)
        if(dodble)npts = npts * 2
        if(npts.gt.NSAMP)npts=NSAMP
c-----
c       read the eigenfunction headers
c-----
        call gtemod(depths,depthr,nper,mname,ipar,fpar)
c-----
c       process the spectra for this distance
c-----
        np0=npts
c-----
c       output the file16 header for this distance
c-----
        df=1./(float(np0)*dt)
c------
c       generate the source spectrum.
c       Due to limitation of available space, when np0>NSAMP
c       use 'supfft' to do FFT.
c------
            do 330 i=1,NSAMP
                ssrc(i) = cmplx(0.0,0.0)
  330       continue
            do 340 i=1,NSAMP
                ssrc(i)=cmplx(src(i),0.0)
  340       continue
            call zfour(ssrc,np0,-1,dt,df)
c------
c       open temporary output file for spectra
c-----
        open(10,file='keep.pl0',status='unknown',form='unformatted')
        iopen0=1
        rewind 10

c------
c       the main call to process the spectra for this distance.
c------
        call gtegn(dist,tfirst,ipar)
c------
c       output.
c------
        if(itrigl.eq.1) close(7,status='delete')
        if(itrigr.eq.1) close(8,status='delete')
        if(kkl.eq.1) close(1)
        if(kkr.eq.1) close(2)
        call output(ntau,idva,iodva,ieqex,dist,dt,npts,
     1      tfirst-0.5*duration,depthr,depths,kkl,kkr,ostr,
     2      mname,ipar,fpar,dodble,dolock,dolat)
        if(iopen0.eq.1) close(10,status='delete')
        return
        end

        subroutine npow2(npts)
c-----
c       Given npts, determine the N=2**m such that N >= npts
c       return the new ntps
c-----
        integer*4 nsmp, npts
        nsmp = npts
        npts = 1
 1000   continue
        npts = 2*npts
        if(npts.lt.nsmp)goto 1000
        return
        end

        subroutine output(ntau,idva,iodva,ieqex,dist,dt,nptin,
     1      tfirst,depthr,depths,kkl,kkr,ostr,
     2      mname,ipar,fpar,dodble,dolock,dolat)
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
c       kkl I*4 - = 1 indicates Love wave eigenfunctions
c       kkr I*4 - = 1 indicates Rayleigh wave eigenfunctions
c       ostr    C*(*)   - command line string
c       mname   C*(*)   - name of model file
c       ipar    I*4 - arrau of integer command parameters from sdisp96
c       fpar    R*4 - arrau of float command parameters from sdisp96
c       dodble  L   - if .true. double length of time series for FFT
c       dolock  L   - if .true. do not use bottom layer for 
c                   first arrival computation
c-----
c-----
c       variables in common blocks
c
c       iftype  I*4 File type
c               1 - single trace
c               3 - three component
c               16 - Green s function
c               21 - Green s function
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
     2      tp, tsv, tsh,
     3      sa, sc, sf, sl, sn, sr
        real*4 esec, distkm, distdg, evstaz, stevaz
        real*4 evlat, evlon, evdep 
        real*4 stlat, stlon, stelev
        real*4 tp, tsv, tsh
        real*4 sa, sc, sf, sl, sn, sr

        common/chdr96/ stname, cfilt, cpulse, ccomnt
        character*8 stname*8, cfilt*80, cpulse*80, ccomnt*80    

        character stcomp*8
        real*4 cmpaz, cmpinc, cmpdt
        integer*4 npts
        real*4 ssec

        parameter(NL=200)
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL),
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax

c-----
c       internal variables
c-----
        parameter (LER=0, LIN=5, LOT=6)
        parameter (NGRN=21)
        integer*4 ntau, idva, iodva, ieqex, nptin
        real*4 dt, tfirst, depths, depthr
        character ostr*(*)

        integer kerr, system

        integer NSAMP2
        parameter(NSAMP2=16384)
        common/dd/datc
        complex datc(NSAMP2)
        common/xx/x
        real*4 x(NSAMP2)
        complex xx(NGRN)
        complex z
        integer iszrt(NGRN)
        character*8 ost(NGRN)

        character mname*80
        integer ipar(20)
        real*4 fpar(20)
        logical dodble, dolock, dolat

        real svel, pvel, vsa, vsb, vsr, geom, rayp, tstar, den

        data iszrt/1,4, 1,4,5, 1,4,5, 1,4, 1,4, 1,4,5,1,1,1,1,1,1/
        data ost/'ZDD     ', 'RDD     ', 'ZDS     ', 'RDS     ',
     1       'TDS     ', 'ZSS     ', 'RSS     ', 'TSS     ',
     2       'ZEX     ', 'REX     ', 'ZVF     ', 'RVF     ',
     3       'ZHF     ', 'RHF     ', 'THF     ', 'PEX     ',
     4       'PDD     ', 'PDS     ', 'PSS     ', 'PVF     ', 
     5       'PHF     '/ 

c-----
c       establish the jsrc array for the file16(V) header
c-----
        do 100 i=1,21
            jsrc(i) = iszrt(i)
            if(ieqex.eq.0 .and. i.ge.11)then
                jsrc(i) = 0
            else if(ieqex.eq.1 .and. i.le.8)then
                jsrc(i) = 0
            endif
  100   continue
c-----
c       receiver in fluid
c-----
        if(ipar(3).eq.1)then
            jsrc(16) = 1
            jsrc(17) = 1
            jsrc(18) = 1
            jsrc(19) = 1
            jsrc(20) = 1
            jsrc(21) = 1
        else
            jsrc(16) = 0
            jsrc(17) = 0
            jsrc(18) = 0
            jsrc(19) = 0
            jsrc(20) = 0
            jsrc(21) = 0
        endif
c-----
c       no Love waves
c-----
        if(kkl.eq.0)then
            jsrc(5) = 0
            jsrc(8) = 0
            jsrc(15) = 0
        endif
c-----
c       no Rayleigh waves
c-----
        if(kkr.eq.0)then
            jsrc( 1) = 0
            jsrc( 2) = 0
            jsrc( 3) = 0
            jsrc( 4) = 0
            jsrc( 6) = 0
            jsrc( 7) = 0
            jsrc( 9) = 0
            jsrc(10) = 0
            jsrc(11) = 0
            jsrc(12) = 0
            jsrc(13) = 0
            jsrc(14) = 0
            jsrc(16) = 0
            jsrc(17) = 0
            jsrc(18) = 0
            jsrc(19) = 0
            jsrc(20) = 0
            jsrc(21) = 0
        endif
c-----
c       source is in fluid
c-----
        if(ipar(2).eq.1)then
            jsrc(1) = 0
            jsrc(2) = 0
            jsrc(3) = 0
            jsrc(4) = 0
            jsrc(5) = 0
            jsrc(6) = 0
            jsrc(7) = 0
            jsrc(8) = 0
            jsrc(9) = 1
            jsrc(10) = 1
            jsrc(11) = 0
            jsrc(12) = 0
            jsrc(13) = 0
            jsrc(14) = 0
            jsrc(15) = 0
            jsrc(16) = 1
            jsrc(17) = 0
            jsrc(18) = 0
            jsrc(19) = 0
            jsrc(20) = 0
            jsrc(21) = 0
        endif
        tau = ntau * dt
        los = lgstr(ostr)
        cpulse = ostr(1:los)
        if(idva.eq.0)then
            iunit = 2
        else if(idva.eq.1)then
            iunit = 3
        else if(idva.eq.2)then
            iunit = 4
        endif
        if(iodva.eq.0)then
            iunit = 2
        else if(iodva.eq.1)then
            iunit = 3
        else if(iodva.eq.2)then
            iunit = 4
        endif
        junit = 11
            iftype = NGRN
            iobsyn = 2
            itmfrq = 1
            cfilt = 'None'
            keyear = 0
            kemon = 0
            keday = 0
            kehour = 0
            kemin = 0
            esec = 0.0
            evlat = -12345
            evlon = -12345
            evdep = depths

            stname = 'GRN21'
            stlat  = -12345
            stlon = -12345
            stelev = depthr
            distkm = dist
            distdg = dist/111.195
            evstaz = 0.0
            stevaz = 180.0
            lmnm = lgstr(mname)
            ccomnt = mname(1:lmnm)
c-----
c       get first arrival time - NOTE for implementation of
c-----
            CALL FRSTAR(dist,depths,depthr,
     1          MNAME,1,TP ,
     1          SSA, SSC, SSF, SSL, SSN, SSR,
     1          RRA, RRC, RRF, RRL, RRN, RRR,
     1          rayp, geom, tstar, dolock)
            CALL FRSTAR(dist,depths,depthr,
     1          MNAME,2,TSV,
     1          SSA, SSC, SSF, SSL, SSN, SSR,
     1          RRA, RRC, RRF, RRL, RRN, RRR,
     1          rayp, geom, tstar, dolock)
            CALL FRSTAR(dist,depths,depthr,
     1          MNAME,3,TSH,
     1          SSA, SSC, SSF, SSL, SSN, SSR,
     1          RRA, RRC, RRF, RRL, RRN, RRR,
     1          rayp, geom, tstar, dolock)
           

            sa = SSA
            sc = SSC
            sl = SSL
            sn = SSN
            sf = SSF
            sr = SSR
c-----
c       lat2d96 we cannot easily predict the first arrival time
c       of a 2-D model. So we set this value to the default for
c       sac, e.g., -12345
c       JUST IGNORE PREVIOUSE - however we must preserve the sa sc cl sn sf sr
c-----
            if(dolat)then
                tp   = -12345.
                tsv   = -12345.
                tsh   = -12345.
            endif
C       write(0,*)dist,depths,depthr,tp,tsv,tsh
c-----
c           Output header
c-----
            call wrhd96(LOT,nerr)
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
        do 1300 jk=1,NGRN
            if(jsrc(jk).ne.0)then
                rewind 10
                do 1301 k=np2,1,-1
                    read(10) (xx(i),i=1,NGRN)
                        datc(k)=xx(jk)
                        if(k.gt.1)then
                        datc(nptin+2-k)=conjg(datc(k))
                        endif
 1301           continue
c-----
c       perform inverse FFT
c-----
c       also take the opportunity to perform the correct
c       spherical spreading for an earth with radius of 6371.0 km
c-----
                if(ipar(1).gt.0)then
                    rad = distdg*3.1415927/180.0
                    fac = sqrt(dist/(6371.0*sin(rad)) )
                else
                    fac = 1.0
                endif
                    call zfour(datc,nptin,+1,dt,df)
                    do 1303 k=1,nptin
                        x(k) = real(datc(k))*fac
 1303               continue
                if(dodble)then
                    npts = nptin / 2
                else
                    npts = nptin
                endif
c-----
c           convert velocity to displacement or acceleration
c           if required. Since we are using surface waves, the
c           DC component must be zero!
c
c           note the special case for the pressure field
c-----
                if(jk.le.15)then
                if(idva.eq.0)then
                    call integ(x,npts,dt,.false.)
                else if(idva.eq.2)then
                    call deriv(x,npts,dt)
                endif
                else if(jk.ge.16)then
                    call integ(x,npts,dt,.false.)
                endif
c-----
c           vertical
c-----
                if(iszrt(jk).eq.1)then
                    cmpinc = -90.0
                    cmpaz  =   0.0
c-----
c           radial
c-----
                else if(iszrt(jk).eq.4)then
                    cmpinc = 0.0
                    cmpaz  =   0.0
c-----
c           transverse
c-----
                else if(iszrt(jk).eq.5)then
                    cmpinc = 0.0
                    cmpaz  =  90.0
                endif
                stcomp = ost(jk)
                cmpdt = dt
                ksyear = 0
                ksmon = 0
                ksday = 0
                kshour = 0
                ksmin = 0
                ssec = tfirst
                call wrtr96(LOT,stcomp,cmpinc,cmpaz,
     1              cmpdt, npts, ksyear, ksmon, 
     2              ksday, kshour, ksmin, ssec, 
     3              x,nerr,NSAMP2)
            endif
 1300   continue
        return
        end


        subroutine gtemod(depths,depthr,nper,mname,ipar,fpar)
c------
c       read in the earth model. consistency check between files.
c-----
        parameter (LIN=5, LOT=6, LER=0)
        parameter (NGRN=21, NLAY=200)
        real*4 d(NLAY),ta(NLAY),tc(NLAY),tf(NLAY),tl(NLAY),tn(NLAY),
     1       rho(NLAY),qa1(NLAY),qb1(NLAY)
        integer*4 kkr,kkl,kkf,mmax,nperl,mperl,nperr,mperr,nper,np0
        common/ctrl/   sr,si,xmom,zeta,ms,np0,kkr,kkl,
     *                 kkf,itrigl,itrigr,mods,modt,modi
        common/resp/   df,bmax,modes1(2),modes2(2),pers1(2),pers2(2),
     *                 ipart1(2),ipart2(2)

        character mname*80
        integer ipar(20)
        real*4 fpar(20)
        common/water/rrho
c-----
        dphsr=0.0
        dphre=0.0
        mperr = 0
        mperl = 0
        ijk=1
        if(kkl.eq.1) then
            rewind 1
C            call getmdl(1,mmax,d,a,b,rho,qa1,qb1,nperl,dphsrl,dphrel,
C     2              mname,ipar,fpar)
            call getmdt(1,mmax,d,ta,tc,tf,
     1           tl,tn,rho,qa1,qb1,nperl,dphsrl,dpthrel,
     1           mname,ipar,fpar)

            depths=dphsrl
            depthr=dphrel
            nper=nperl
        endif
        if(kkr.eq.1) then
            rewind 2
            call getmdt(2,mmax,d,ta,tc,tf,
     1           tl,tn,rho,qa1,qb1,nperr,dphsrr,dpthrer,
     1           mname,ipar,fpar)
C            call getmdl(2,mmax,d,a,b,rho,qa1,qb1,nperr,dphsrr,dphrer,
C     2              mname,ipar,fpar)
            depths=dphsrr
            depthr=dphrer
            nper=nperr
            deplw = 0.0
            depup = 0.0
            do 1000 i=1,mmax
                depup = deplw + d(i)
                if(depthr.ge.deplw .and. depthr.lt.depup)then
                    rrho = rho(i)
                endif
                deplw = depup
 1000       continue
        endif
c       write(LER,8)
c       write(LER,10)
C       write(LER,20) (i,d(i),a(i),b(i),rho(i),qa1(i),qb1(i),i=1,mmax)
C 07/29/2010 get rid of bmax - never used
        bmax=0.0
c       write(LER,30) dphsr,dphre,nperl,mperl,nperr,mperr
        return
        end

        subroutine gtegn(rr,tshft,ipar)
c-----
c       set the end values for interpolation.
c-----
        integer ipar(20)
        parameter (LIN=5, LOT=6, LER=0)
        parameter (NGRN=21)
        parameter (NMD=100)
        integer*4 itrig(2)
        integer*4 kkr,kkl,kkf,ifunc(2),mode2,np0
        integer*4 npx,isign
        integer NSAMP
        parameter (NSAMP=16384)
        common/srcspc/ssrc(NSAMP)
        complex ssrc
        common/rayl/
     1      ur(NMD,2),dur(NMD,2),uz(NMD,2),duz(NMD,2),
     2      ur0(NMD,2),uz0(NMD,2),tr0(NMD,2),tz0(NMD,2),
     3      wvmr(NMD,2),ares(NMD,2),arer(NMD,2),gamr(NMD,2),
     4      wvmrr(NMD,2), wvmrs(NMD,2)
        common/love/
     1      ut(NMD,2),dut(NMD,2),ut0(NMD,2),dut0(NMD,2),
     2      wvml(NMD,2),ales(NMD,2),aler(NMD,2),gaml(NMD,2),
     3      wvmlr(NMD,2), wvmls(NMD,2)
        common/ctrl/   sr,si,xmom,zeta,ms,np0,kkr,kkl,
     *      kkf,itrigl,itrigr,mods,modt,modi
        common/resp/   df,bmax,modes1(2),modes2(2),pers1(2),pers2(2),
     *      ipart1(2),ipart2(2)
        complex xx(NGRN)
        np2=np0/2+1
c------
c     connect with the input data files.
c------
        call readin(0,jj,ifunc,mode2,ipart,per)
        do 150 jj=1,2
            j=jj
            if(jj.eq.1.and.kkl.ne.1) go to 150
            if(jj.eq.2.and.kkr.ne.1) go to 150
            itrig(j)=0
            igo=j
            call readin(igo,j,ifunc,mode2,ipar2,per2)
            modes2(j)=mode2
            ipart2(j)=ipar2
            pers2(j)=per2
  150   continue
c------
c     start the job.
c------
        do 600 ii=np2,1,-1
            do 200 i=1,NGRN
                xx(i)=cmplx(0.0,0.0)
  200       continue
            per=-1.0
            if(ms.ne.1) then
                    sr = real(ssrc(ii))
                    si = aimag(ssrc(ii))
            endif
            if(ii.eq.1) then
                per=1.e+5
                if(ii.eq.1) go to 320
            endif
            per=1./((ii-1)*df)
            if(ii.eq.2) go to 320
            dper=1./((ii-2)*df)-per
  320       continue
            per0=per-0.005*dper
            if(ii.eq.np2) per0=per+0.005*dper
c
            do 500 jj=1,2
                j=jj
                if(jj.eq.1.and.kkl.ne.1) go to 500
                if(jj.eq.2.and.kkr.ne.1) go to 500
c------
c     find the data at two end points. 
c            Interpolation is done between them.
c------
  330           if(per0.lt.pers2(j)) go to 400
                if(itrig(j).eq.2) go to 500
                pers1(j)=pers2(j)
                modes1(j)=modes2(j)
                ipart1(j)=ipart2(j)
                igo=j
                call readin(igo,j,ifunc,mode2,ipar2,per2)
                modes2(j)=mode2
                ipart2(j)=ipar2
                pers2(j)=per2
                itrig(j)=1
                if(ifunc(igo).le.0) itrig(j)=2
                go to 330
  400           if(itrig(j).eq.0) go to 500
                mode0=modes2(j)
                if(mode0.gt.modes1(j)) mode0=modes1(j)
                if(jj.eq.1) then
                     call excitl(per,mode0,rr,tshft,xx,ipar)
                else if(jj.eq.2) then
                     call excitr(per,mode0,rr,tshft,xx,ipar)
                endif
  500       continue
c------
c     output
c     16 complex values correspond to ten Green s functions:
c       01: ZDD     02: RDD     17: PDD
c       03: ZDS     04: RDS     18: PDS
c       05: TDS     06: ZSS     19: PSS
c       07: RSS     08: TSS     20: PVF
c       09: ZEX     10: REX     21: PHF
c       11: ZVF     12: RVF
c       13: ZHF     14: RHF
c       15: THF     16: PEX
c------
            write(10) (xx(i),i=1,NGRN)
  600   continue
        return
        write(LOT,*) 'Period not agree between source and receiver'
        write(LOT,*) 'files at period=',per2,perr2
        stop
        end


        subroutine srcpul(ntau,ipt,alp,dt,npts,rfile,duration)
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
c-----
        parameter (LIN=5, LOT=6, LER=0)
        integer*4 ntau, ipt, npts
        real*4 alp, dt
        character*80 rfile
        integer*4 kkr,kkl,kkf,np0
        common/ctrl/   sr,si,xmom,zeta,ms,np0,kkr,kkl,
     1      kkf,itrigl,itrigr,mods,modt,modi
        integer NSAMP
        parameter (NSAMP=16384)
        common/srctim/ src(NSAMP)
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
        return
        end
 
        subroutine readin(igo,jj,ifunc,mode,ipart,per)
c-----
c       read in the eigenfunction files generated by
c       reigen85 and leigen85
c
c       use the definitions of the eigenfunction and stress to get
c       the necessary vertical derivatives
c
c       Instead of reading only a single mode at a time, groups of
c       NMD are read in - the purpose of this is to 
c            significantly reduce
c       IO time. If not for this the code would be much simpler
c-----
c       igo I*4 - 0
c                 1
c                 2
c       jj  I*4 -
c       ifunc   I*4 - ifunc for L and R
c       mode    I*4 -
c       ipart   I*4 -
c       per R*4 -
c-----
        parameter (NGRN=21)
        parameter (NMD=100)
        integer*4 kkr,kkl,kkf,ifunc(2),mode,np0
        integer NSAMP
        parameter (NSAMP=16384)
        common/srcspc/ssrc(NSAMP)
        complex ssrc
        common/rayl/
     1      ur(NMD,2),dur(NMD,2),uz(NMD,2),duz(NMD,2),
     2      ur0(NMD,2),uz0(NMD,2),tr0(NMD,2),tz0(NMD,2),
     3      wvmr(NMD,2),ares(NMD,2),arer(NMD,2),gamr(NMD,2),
     4      wvmrr(NMD,2), wvmrs(NMD,2)
        common/love/
     1      ut(NMD,2),dut(NMD,2),ut0(NMD,2),dut0(NMD,2),
     2      wvml(NMD,2),ales(NMD,2),aler(NMD,2),gaml(NMD,2),
     3      wvmlr(NMD,2), wvmls(NMD,2)
        common/ctrl/   sr,si,xmom,zeta,ms,np0,kkr,kkl,
     *                 kkf,itrigl,itrigr,mods,modt,modi
        common/resp/   df,bmax,modes1(2),modes2(2),pers1(2),pers2(2),
     *                 ipart1(2),ipart2(2)
        common/where/  pp(2,2)
        parameter (NLAY=200)
        real*4 d(NLAY),ta(NLAY),tc(NLAY),tf(NLAY), tl(NLAY), tn(NLAY),
     1     rho(NLAY),qa1(NLAY),qb1(NLAY)

        character mname*80
        integer ipar(20)
        real*4 fpar(20)
c-----
c       set fpar(1) = refdep
c       set ipar(1) = 1 if medium is spherical
c       set ipar(2) = 1 if source is in fluid
c       set ipar(3) = 1 if receiver is in fluid
c       set ipar(4) = 1 if eigenfunctions are output with -DER flag
c       set ipar(5) = 1 if dc/dh are output with -DER flag
c       set ipar(6) = 1 if dc/da are output with -DER flag
c       set ipar(7) = 1 if dc/db are output with -DER flag
c       set ipar(8) = 1 if dc/dr are output with -DER flag
c-----
c       initialize
c-----
        if(igo.eq.0)then
            do 90 k=1,2
                ipart1(k)=0
                ipart2(k)=0
   90       continue
c------
c       skip the first several records, which contain the model.
c       This is required only for proper input positioning
c------
            if(kkl.eq.1) then
                rewind 1
            call getmdt(1,mmax,d,ta,tc,tf,
     1           tl,tn,rho,qa1,qb1,nper,depths,depthr,
     1           mname,ipar,fpar)
C            call getmdl(1,mmax,d,a,b,rho,qa1,qb1,nperl,dphsrl,dphrel,
C     2              mname,ipar,fpar)
            endif
            if(kkr.eq.1) then
                rewind 2
            call getmdt(2,mmax,d,ta,tc,tf,
     1           tl,tn,rho,qa1,qb1,nper,depths,depthr,
     1           mname,ipar,fpar)
C            call getmdl(2,mmax,d,a,b,rho,qa1,qb1,nperr,dphsrl,dphrel,
C     2              mname,ipar,fpar)
            endif
        else
c-----
c       start reading the eigenfunction file 
c       The structure of the file, per read, is
c       earth model
c       number_periods, source_depth
c            | ifunc,mode,per [ifunc=1 for Lv, 
c            2 for RY, -1 for end of file]
c           p|m|        [ mode - number of modes ]
c           e|o|        [ per period of this entry, < 0 
c            end of per set ]
c           r|d| eigenfunctions for this mode and period
c-----
            icode=igo
            call gethed(icode,ifunc(icode),mode,per,ierr)
            if(ifunc(icode).le.0) return
c------
c     Use a buffer with NMD data length to handle large mode number.
c     IF f77 allows bigger dimension, the change in the parameter NMD
c     number will make this program faster.
c
c       Current Love IO has  7 variables, need recl= 7*4*NMD
c       Current Rayl IO has 11 variables, need recl=11*4*NMD
c------
            ipart=(mode-1)/NMD+1
            if(mode.le.0) return
            if(ipart.gt.1) then
                if(itrigl.eq.0.and.jj.eq.1) then
                    locl=0
                    itrigl=1
                    open(7,file='keep.pl7',
     1                  status='unknown',
     2                  form='unformatted',recl=4096
     3                  ,access='direct')
                endif
                if(itrigr.eq.0.and.jj.eq.2) then
                    locr=0
                    itrigr=1
                    open(8,file='keep.pl8',
     1                  status='unknown',
     2                  form='unformatted',recl=6144
     3                  ,access='direct')
                endif
            endif
c------
c     use double buffer. This is a difficult control.
c------
            if(itrigl.eq.1.and.jj.eq.1.and.ipart.gt.1.and.igo.le.2)
     1          locl=locl+1
            if(itrigr.eq.1.and.jj.eq.2.and.ipart.gt.1.and.igo.le.2)
     1          locr=locr+1
c
            do 300 i=1,ipart
                j2=NMD
                if(i.eq.ipart) j2=mod(mode,NMD)
                if(j2.eq.0) j2=NMD
c
                do 200 j=1,j2
                    if(igo.eq.1)then
c-----
c                       READ LOVE
c-----
                        wvml(j,1) = wvml(j,2) 
                        wvmls(j,1) = wvmls(j,2) 
                        wvmlr(j,1) = wvmlr(j,2) 
                        ales(j,1)  = ales(j,2) 
                        aler(j,1)  = aler(j,2) 
                        gaml(j,1)  = gaml(j,2) 
                        ut(j,1)    = ut(j,2) 
                        ut0(j,1)   = ut0(j,2) 
                        dut(j,1)   = dut(j,2) 
                        dut0(j,1)  = dut0(j,2) 
                        call getegn(icode,ifunc(icode),1,wvno,u,gamma,
     1                      sur,sdur,suz,sduz,sare,wvnsrc,sur0,
     2                      rur,rdur,ruz,rduz,rare,wvnrec,rur0,
     3                      sumkr,sumgr,sumgv,ierr)
                        wvml(j,2) = wvno
c-----
c                       Get Love - 1 = normal, 3 = lat2d
c                       actually this is all the same
c-----
                        wvmlr(j,2) = wvnrec
                        wvmls(j,2) = wvnsrc
                        ales(j,2) = sare
                        aler(j,2) = rare

                        gaml(j,2)  = gamma
                        ut(j,2)    = sur
                        dut(j,2)   = sdur
                        ut0(j,2)   = rur
                        dut0(j,2)  = rdur
                    else if(igo.eq.2)then
c-----
c                       READ RAYLEIGH
c-----
                        wvmr(j,1)  = wvmr(j,2)
                        wvmrs(j,1) = wvmrs(j,2)
                        wvmrr(j,1) = wvmrr(j,2)
                        ares(j,1)  = ares(j,2)
                        arer(j,1)  = arer(j,2)
                        gamr(j,1)  = gamr(j,2)
                        ur(j,1)    = ur(j,2)
                        dur(j,1)   = dur(j,2)
                        uz(j,1)    = uz(j,2)
                        duz(j,1)   = duz(j,2)
                        ur0(j,1)   = ur0(j,2)
                        tr0(j,1)  = tr0(j,2)
                        uz0(j,1)   = uz0(j,2)
                        tz0(j,1)  = tz0(j,2)
                        call getegn(icode,ifunc(icode),1,wvno,u,gamma,
     1                      sur,sdur,suz,sduz,sare,wvnsrc,sur0,
     2                      rur,rdur,ruz,rduz,rare,wvnrec,rur0,
     3                      sumkr,sumgr,sumgv,ierr)
                        wvmr(j,2) = wvno
c-----
c                       Get Rayleigh - 2 = normal, 4 = lat2d
c                       actually this is all the same
c-----
                        wvmrr(j,2) = wvnrec
                        wvmrs(j,2) = wvnsrc
                        arer(j,2) = rare
                        ares(j,2) = sare
                        
                        gamr(j,2) = gamma
                        ur(j,2)   = sur
                        dur(j,2)  = sdur
                        uz(j,2)   = suz
                        duz(j,2)  = sduz

                        ur0(j,2)  = rur
                        tr0(j,2) = rdur
                        uz0(j,2)  = ruz
                        tz0(j,2) = rduz
                    endif
  200           continue
                if(ipart.gt.1) then
                    if(jj.eq.1) then
                        irec=i
                        if(locl.ge.3) locl=1
                        if(locl.eq.2) irec=irec+40
                        pp(1,locl)=per
                        write(7,rec=irec) (wvml(k,2),
     1              ales(k,2),aler(k,2),gaml(k,2),ut(k,2),
     2                      dut(k,2),ut0(k,2),
     3                      dut0(k,2),wvmlr(k,2),
     4                      wvmls(k,2),
     5                      k=1,NMD)
                    endif
                    if(jj.eq.2) then
                        irec=i
                        if(locr.ge.3) locr=1
                        if(locr.eq.2) irec=irec+40
                        pp(2,locr)=per
                        write(8,rec=irec) (wvmr(k,2),
     1              ares(k,2),arer(k,2),gamr(k,2),ur(k,2),
     2                  dur(k,2),uz(k,2),duz(k,2),
     3                  ur0(k,2),tr0(k,2),uz0(k,2),
     4                  tz0(k,2),
     5                  wvmrs(k,2),wvmrr(k,2),
     6                      k=1,NMD)
                    endif
                endif
c-----
c               DONE WITH WAVE TYPE
c-----
  300       continue
        endif
        return
        end
 
        subroutine finds(m,ii,j1,j2,j11,j22,j33,n1,n2,n3,mode)
c-----
c       In order to have no limit on the maximum number of modes
c       the modal output is written in groups of 100 modes
c       find the part where the modes>100 are stored.
c-----
        parameter (NGRN=21)
        parameter (NMD=100)
        integer*4 kkr,kkl,kkf,np0
        integer NSAMP
        parameter (NSAMP=16384)
        common/srcspc/ssrc(NSAMP)
        complex ssrc
        common/rayl/
     1      ur(NMD,2),dur(NMD,2),uz(NMD,2),duz(NMD,2),
     2      ur0(NMD,2),uz0(NMD,2),tr0(NMD,2),tz0(NMD,2),
     3      wvmr(NMD,2),ares(NMD,2),arer(NMD,2),gamr(NMD,2),
     4      wvmrr(NMD,2), wvmrs(NMD,2)
        common/love/
     1      ut(NMD,2),dut(NMD,2),ut0(NMD,2),dut0(NMD,2),
     2      wvml(NMD,2),ales(NMD,2),aler(NMD,2),gaml(NMD,2),
     3      wvmlr(NMD,2), wvmls(NMD,2)
        common/ctrl/   sr,si,xmom,zeta,ms,np0,kkr,kkl,
     *      kkf,itrigl,itrigr,mods,modt,modi
        common/resp/   df,bmax,modes1(2),modes2(2),pers1(2),pers2(2),
     *      ipart1(2),ipart2(2)
        common/where/  pp(2,2)
        data epi/1.0e-8/
c-----
c
c-----
        if(ipart1(m).gt.1) then
            n1=n1+1
            krec0=0
c------
c     If want max mode < 5000, change krec0=40 to krec0=100
c       if(pp(m,2).eq.pers1(m)) krec0=40
c------
            if(pers1(m).ge.(pp(m,2)-epi)
     1          .and.pers1(m).le.(pp(m,2)+epi) )then
                    krec0=40
            endif
            irec=ii+krec0
            if(m.eq.1) then
                read(7,rec=irec) (wvml(k,1),
     1          ales(k,1),aler(k,1),gaml(k,1),ut(k,1),
     2              dut(k,1),ut0(k,1),
     3              dut0(k,1),wvmlr(k,1),wvmls(k,1),
     4              k=1,NMD)
            endif
            if(m.eq.2) then
                read(8,rec=irec) (wvmr(k,1),
     1          ares(k,1),arer(k,1),gamr(k,1),ur(k,1),
     2              dur(k,1),uz(k,1),duz(k,1),
     3              ur0(k,1),tr0(k,1),uz0(k,1),
     4              tz0(k,1),
     5              wvmrs(k,1),wvmrr(k,1),
     6                  k=1,NMD)
            endif
        endif
        if(ipart2(m).gt.1) then
            n2=n2+1
            krec0=0
c------
c     If want max mode < 5000, change krec0=40 to krec0=100
c       if(pp(m,2).eq.pers2(m)) krec0=40
c------
            if(pers2(m).ge.(pp(m,2)-epi)
     1          .and.pers2(m).le.(pp(m,2)+epi) )then
                krec0=40
            endif
            irec=ii+krec0
            if(m.eq.1) then
                read(7,rec=irec) (wvml(k,2),
     1          ales(k,2),aler(k,2),gaml(k,2),ut(k,2),
     2              dut(k,2),ut0(k,2),
     3              dut0(k,2),wvmlr(k,2),wvmls(k,2),
     4              k=1,NMD)
            endif
            if(m.eq.2) then
                read(8,rec=irec) (wvmr(k,2),
     1          ares(k,2),arer(k,2),gamr(k,2),ur(k,2),
     2              dur(k,2),uz(k,2),duz(k,2),
     3              ur0(k,2),tr0(k,2),uz0(k,2),
     4              tz0(k,2),
     5              wvmrs(k,2),wvmrr(k,2),
     6                  k=1,NMD)
            endif
        endif
        if(n1*n2.eq.0) then
            j11=j1
            j22=j2
            j33=j2
            n3=1
        else
            j11=1
c------
c     find the range of mode to be treated.
c------
            j22=NMD
            if(ipart1(m).le.ipart2(m)) then
                if(n1.eq.ipart1(m)) j22=mod(mode,NMD)
                n3=n1
            else
                if(n2.eq.ipart2(m)) j22=mod(mode,NMD)
                n3=n2
            endif
            if(j22.eq.0) j22=NMD
            j33=(n3-1)*NMD+j22
        endif
        return
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

        subroutine excitr(per,mode,rr,tshft,xx,ipar)
c-----
c       ieqex   I*4 - 0 Earthquake + Explosion
c               - 1 Explosion + Point Force
c               - 2 Earthquake, explosion, point force
c-----
c
c     This generates the Z and R components of
c     seismogram for
c      1: 45-deg dip-slip source  2: strike-slip source
c      3: dip-slip source         4: explosion source.
c
c     number of mode is unlimited.
c
        integer ipar(20)
        parameter (NGRN=21)
        parameter (NMD=100)
        real*8 rat,atn1,atn2,fact1,fact2,v1,v2,w1,w2,u1,u2
        real*8 d1, d2
        real*8 p1,p2
        real*8 t1,t2,ct1,ct2,st1,st2,sr,si,xmom,rx,tx
        real*8 wvno,omega,pi,pi4,pi34
        real*8 dkz(6),dkr(6),dkp(6)
        complex*16 cvz(6), cvr(6), cvp(6)
        integer*4 kkr,kkl,kkf,np0,np
        integer NSAMP
        parameter (NSAMP=16384)
        common/srcspc/ssrc(NSAMP)
        complex ssrc
        common/rayl/
     1      ur(NMD,2),dur(NMD,2),uz(NMD,2),duz(NMD,2),
     2      ur0(NMD,2),uz0(NMD,2),tr0(NMD,2),tz0(NMD,2),
     3      wvmr(NMD,2),ares(NMD,2),arer(NMD,2),gamr(NMD,2),
     4      wvmrr(NMD,2), wvmrs(NMD,2)
        common/ctrl/   sr0,si0,xmom0,zeta,ms,np0,kkr,kkl,
     *                 kkf,itrigl,itrigr,mods,modt,modi
        common/resp/   df,bmax,modes1(2),modes2(2),pers1(2),pers2(2),
     *                 ipart1(2),ipart2(2)

        common/water/rrho

        complex xx(NGRN)
        data pi/3.141592653589793d+00/
c------
        xmom=dble(xmom0)/dsqrt(2.d+00*pi)
        pi4=pi/4.d+00
        pi34=3.d+00*pi/4.d+00
        do 91 j=1,6
            cvz(j)  = dcmplx(0.0d+00,0.0d+00)
            cvr(j)  = dcmplx(0.0d+00,0.0d+00)
            cvp(j)  = dcmplx(0.0d+00,0.0d+00)
   91   continue

        j1=mods
        j2=modt
        if(j2.gt.mode) j2=mode
        j3=modi

        if(j1.gt.j2) return
        omega=2.d+00*pi/dble(per)
        per1=pers1(2)
        per2=pers2(2)
        n1=0
        n2=0
        ii=0
  100   continue
        ii=ii+1
c------
c       find where the data locate.
c------
        call finds(2,ii,j1,j2,j11,j22,j33,n1,n2,n3,mode)
c------
c       put different elements into the final solution 
c            after interpolation.
c------
        icnt=0
        do 400 j=j11,j22
        if(ares(j,1).lt.1.e-30.and.ares(j,2).lt.1.0e-30 .and.
     1  arer(j,1).lt.1.e-30.and.arer(j,2).lt.1.0e-30) then
          icnt=icnt+1
        else
          icnt=0
        endif
        if(icnt.gt.3) go to 410
        if(ares(j,1).lt.1.e-30.or.ares(j,2).lt.1.0e-30 .and.
     1  arer(j,1).lt.1.e-30.or.arer(j,2).lt.1.0e-30) go to 400
        k=(n3-1)*NMD+j
        if(k.ge.j1.and.k.le.j2) then
          if(mod(k-j1,j3).eq.0) go to 150
        endif
        go to 400
  150   continue
        rat=0.0d+00
        rat=(dble(per)-dble(per1))/(dble(per2)-dble(per1))
        wvno=(dble(wvmr(j,1)*per1) + 
     1      dble(wvmr(j,2)*per2-wvmr(j,1)*per1)*rat)/per
c------
c       For different distance using same npt.
c------
        rx=dble(rr)
        tx=dble(tshft)
c
        if(zeta.lt.0.0) then
          fact1=1.d+00
          fact2=1.d+00
          if(fact2.eq.1.d+00) go to 160
        endif
        atn1=1.0d+00
        if(zeta.gt.2.0) then
            atn1=zeta/100.0
        else
            if(per.lt.10.0) then
                atn1=dble(per)**dble(zeta)
            else
                atn1=dble(10.0)**dble(zeta)
            endif
        endif
        atn2=atn1*dble(gamr(j,2))*rx
        atn1=atn1*dble(gamr(j,1))*rx
        fact1=0.0d+00
        if(atn1.lt.80.0) fact1=1./dexp(atn1)
        fact2=0.0d+00
        if(atn2.lt.80.0) fact2=1./dexp(atn2)
  160   continue
c
        v1 = dsqrt(dble(ares(j,1))*dble(arer(j,1)))
     1      /dsqrt(dble(wvmr(j,1))*rx)
        v2 = dsqrt(dble(ares(j,2))*dble(arer(j,2)))
     1      /dsqrt(dble(wvmr(j,2))*rx)
c------
c       DD component.
c------
        d1    = 2.0*dble(duz(j,1))+dble(wvmrs(j,1))*dble(ur(j,1))
        d1    = d1*v1*fact1
        u1    = d1*dble(ur0(j,1))
        w1    = d1*dble(uz0(j,1))
        p1    = -d1*dble(tz0(j,1))
        d2    = 2.0*dble(duz(j,2))+dble(wvmrs(j,2))*dble(ur(j,2))
        d2    = d2*v2*fact2
        u2    = d2*dble(ur0(j,2))
        w2    = d2*dble(uz0(j,2))
        p2    = -d2*dble(tz0(j,2))
        dkz(1) = w1+(w2-w1)*rat
        dkr(1) = u1+(u2-u1)*rat
        dkp(1) = p1+(p2-p1)*rat
c------
c       SS component.
c------
        d1    = dble(wvmrs(j,1))*dble(ur(j,1))
        d1    = d1*v1*fact1
        u1    = d1*dble(ur0(j,1))
        w1    = d1*dble(uz0(j,1))
        p1    = -d1*dble(tz0(j,1))
        d2    = dble(wvmrs(j,2))*dble(ur(j,2))
        d2    = d2*v2*fact2
        u2    = d2*dble(ur0(j,2))
        w2    = d2*dble(uz0(j,2))
        p2    = -d2*dble(tz0(j,2))
        dkz(2) = w1+(w2-w1)*rat
        dkr(2) = u1+(u2-u1)*rat
        dkp(2) = p1+(p2-p1)*rat
c------
c       DS component.
c------
        d1    = dble(wvmrs(j,1))*dble(uz(j,1)) + dble(dur(j,1))
        d1    = d1*v1*fact1
        u1    = d1*dble(ur0(j,1))
        w1    = d1*dble(uz0(j,1))
        p1    = -d1*dble(tz0(j,1))
        d2    = dble(wvmrs(j,2))*dble(uz(j,2)) + dble(dur(j,2))
        d2    = d2*v2*fact2
        u2    = d2*dble(ur0(j,2))
        w2    = d2*dble(uz0(j,2))
        p2    = -d2*dble(tz0(j,2))
        dkz(3) = w1+(w2-w1)*rat
        dkr(3) = u1+(u2-u1)*rat
        dkp(3) = p1+(p2-p1)*rat
c------
c       EXPLOSION.
c------
        d1    = dble(duz(j,1))-dble(wvmrs(j,1))*dble(ur(j,1))
        d1    = d1*v1*fact1
        u1    = d1*dble(ur0(j,1))
        w1    = d1*dble(uz0(j,1))
        p1    = -d1*dble(tz0(j,1))
        d2    = dble(duz(j,2))-dble(wvmrs(j,2))*dble(ur(j,2))
        d2    = d2*v2*fact2
        u2    = d2*dble(ur0(j,2))
        w2    = d2*dble(uz0(j,2))
        p2    = -d2*dble(tz0(j,2))
c-----
        dkz(4) = w1+(w2-w1)*rat
        dkr(4) = u1+(u2-u1)*rat
        dkp(4) = p1+(p2-p1)*rat
c------
c       VF component.
c------
        d1    = dble(uz(j,1))
        d1    = d1*v1*fact1
        u1    = d1*dble(ur0(j,1))
        w1    = d1*dble(uz0(j,1))
        p1    = -d1*dble(tz0(j,1))
        d2    = dble(uz(j,2))
        d2    = d2*v2*fact2
        u2    = d2*dble(ur0(j,2))
        w2    = d2*dble(uz0(j,2))
        p2    = -d2*dble(tz0(j,2))
        dkz(5) = w1+(w2-w1)*rat
        dkr(5) = u1+(u2-u1)*rat
        dkp(5) = p1+(p2-p1)*rat
c------
c       HF component.
c------
        d1    = dble(ur(j,1))
        d1    = d1*v1*fact1
        u1    = d1*dble(ur0(j,1))
        w1    = d1*dble(uz0(j,1))
        p1    = -d1*dble(tz0(j,1))
        d2    = dble(ur(j,2))
        d2    = d2*v2*fact2
        u2    = d2*dble(ur0(j,2))
        w2    = d2*dble(uz0(j,2))
        p2    = -d2*dble(tz0(j,2))
        dkz(6) = w1+(w2-w1)*rat
        dkr(6) = u1+(u2-u1)*rat
        dkp(6) = p1+(p2-p1)*rat
c------
c       1/sqrt(2 pi).  exp(-i pi/4).  exp(-i 3pi/4).
c       A far-field approximation of Hankel function.
c------
c-----
c       introduce time shift and also phase shift from eigenfunction
c-----
        t1 =omega*tx-wvno*rx-pi4
        t2 =omega*tx-wvno*rx-pi34
        ct1=dcos(t1)
        st1=dsin(t1)
        ct2=dcos(t2)
        st2=dsin(t2)
c-----
c       DD - requires real eigenfunction
c-----
        cvz(1) = cvz(1) - dkz(1)*xmom*cmplx(ct1,st1)
        cvr(1) = cvr(1) + dkr(1)*xmom*cmplx(ct2,st2)
        cvp(1) = cvp(1) - dkp(1)*xmom*cmplx(ct1,st1)
c-----
c       DS - requires -i  eigenfunction
c-----
        cvz(3) = cvz(3) + dkz(3)*xmom*cmplx(st1,-ct1)
        cvr(3) = cvr(3) - dkr(3)*xmom*cmplx(st2,-ct2)
        cvp(3) = cvp(3) + dkp(3)*xmom*cmplx(st1,-ct1)
c-----
c       SS - requires real eigenfunction
c-----
        cvz(2) = cvz(2) - dkz(2)*xmom*cmplx(ct1,st1)
        cvr(2) = cvr(2) + dkr(2)*xmom*cmplx(ct2,st2)
        cvp(2) = cvp(2) - dkp(2)*xmom*cmplx(ct1,st1)
c-----
c       EX - requires real eigenfunction
c-----
        cvz(4) = cvz(4) - dkz(4)*xmom*cmplx(ct1,st1)
        cvr(4) = cvr(4) + dkr(4)*xmom*cmplx(ct2,st2)
        cvp(4) = cvp(4) - dkp(4)*xmom*cmplx(ct1,st1)
c-----
c       VF - requires real eigenfunction
c-----
        cvz(5) = cvz(5) - dkz(5)*xmom*cmplx(ct1,st1)
        cvr(5) = cvr(5) + dkr(5)*xmom*cmplx(ct2,st2)
        cvp(5) = cvp(5) - dkp(5)*xmom*cmplx(ct1,st1)
c-----
c       HF - requires -i eigenfunction
c-----
        cvz(6) = cvz(6) + dkz(6)*xmom*cmplx(st1,-ct1)
        cvr(6) = cvr(6) - dkr(6)*xmom*cmplx(st2,-ct2)
        cvp(6) = cvp(6) + dkp(6)*xmom*cmplx(st1,-ct1)
  400   continue
        if(j2.gt.j33) go to 100
c------
c       The standard output is impulse response type. i.e.,
c       ground velocity, no matter what source type is specified.
c------
410   continue
c-----
c       include source spectrum
c-----
        if(ms.eq.1)then
            sr=1.0d+00
            si=0.0d+00
        else
            sr=dble(sr0)
            si=dble(si0)
        endif
        do 500 j=1,6
            cvz(j) = cvz(j) * dcmplx(sr,si)
            cvr(j) = cvr(j) * dcmplx(sr,si)
            cvp(j) = cvp(j) * dcmplx(sr,si)
  500   continue
        do 601 j=1,6
            t1 = cdabs(cvz(j))
            if(t1.le.1.0d-40)cvz(j) = dcmplx(0.0d+00,0.0d+00)
            t2 = cdabs(cvr(j))
            if(t2.le.1.0d-40)cvr(j) = dcmplx(0.0d+00,0.0d+00)
            t2 = cdabs(cvp(j))
            if(t2.le.1.0d-40)cvp(j) = dcmplx(0.0d+00,0.0d+00)
  601   continue
        if(ipar(2).eq.0)then
c-----
c       source in solid
c-----
            xx(1)  =     cvz(1)
            xx(2)  =     cvr(1)
            xx(3)  =     cvz(3)
            xx(4)  =     cvr(3)
            xx(6)  =    -cvz(2)
            xx(7)  =    -cvr(2)
            xx(9)  =     cvz(4)
            xx(10) =     cvr(4)
            xx(11) =     cvz(5)
            xx(12) =     cvr(5)
            xx(13) =     cvz(6)
            xx(14) =     cvr(6)
            if(ipar(3).eq.1)then
c-----
c       receiver in fluid
c-----
                xx(16) =     cvp(4)
                xx(17) =     cvp(1)
                xx(18) =     cvp(3)
                xx(19) =    -cvp(2)
                xx(20) =     cvp(5)
                xx(21) =     cvp(6)
            else
c-----
c       receiver in solid
c-----
                xx(16) = 0.0
                xx(17) = 0.0
                xx(18) = 0.0
                xx(19) = 0.0
                xx(20) = 0.0
                xx(21) = 0.0
            endif
        else 
c-----
c       source in fluid
c-----
            xx(1) = 0.0
            xx(2) = 0.0
            xx(3) = 0.0
            xx(4) = 0.0
            xx(6) = 0.0
            xx(7) = 0.0
            xx(9) =     cvz(4)
            xx(10)=     cvr(4)
            xx(11)= 0.0
            xx(12)= 0.0
            xx(13)= 0.0
            xx(14)= 0.0
            if(ipar(3).eq.1)then
c-----
c       receiver in fluid
c-----
                xx(16) =    cvp(4)
                xx(17) =    0.0
                xx(18) =    0.0
                xx(19) =    0.0
                xx(20) =    0.0
                xx(21) =    0.0
            else
c-----
c       receiver in solid
c-----
                xx(16) = 0.0
                xx(17) = 0.0
                xx(18) = 0.0
                xx(19) = 0.0
                xx(20) = 0.0
                xx(21) = 0.0
            endif
        endif
            
      return
      end
c
c  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
        subroutine excitl(per,mode,rr,tshft,xx,ipar)
c-----
c       ieqex   I*4 - 0 Earthquake + Explosion
c               - 1 Explosion + Point Force
c               - 2 Earthquake, explosion, point force
c-----
c     This generates the T component of seismograms
c     for  1: dip-slip source     2: strike-slip source.
c
c     vt(i,j)   i=1,2  real & imag part    j=1,2  source
c---
        integer ipar(20)
        parameter (NGRN=21)
        complex xx(NGRN)
        parameter (NMD=100)
        common/love/
     1      ut(NMD,2),dut(NMD,2),ut0(NMD,2),dut0(NMD,2),
     2      wvml(NMD,2),ales(NMD,2),aler(NMD,2),gaml(NMD,2),
     3      wvmlr(NMD,2), wvmls(NMD,2)
        common/ctrl/   sr0,si0,xmom0,zeta,ms,np0,kkr,kkl,
     *                 kkf,itrigl,itrigr,mods,modt,modi
        common/resp/   df,bmax,modes1(2),modes2(2),pers1(2),pers2(2),
     *                 ipart1(2),ipart2(2)
        real*8 rx,tx,rat,wvno,atn1,atn2,fact1,fact2
        real*8 v1,v2,w1,w2,t1,ct1,st1,sr,si,xmom,pi,pi4
        real*8 dk(3),omega
        complex*16 cvt(3)
        integer*4 kkr,kkl,kkf,np0
        data pi/3.141592653589793d+00/
c
        xx(5) = cmplx(0.0,0.0)
        xx(8) = cmplx(0.0,0.0)
        xx(15) = cmplx(0.0,0.0)
c-----
c       if source or receiver is in fluid, there is no 
c            far-field SH motion
c-----
        if(ipar(2).eq.1 .or. ipar(3).eq.1)return
        xmom=dble(xmom0)/dsqrt(2.d+00*pi)
        pi4=pi/4.d+00
        do 90 j=1,3
            cvt(j)=dcmplx(0.0d+00,0.0d+00)
   90   continue
c-----
c       initialize processing
c-----

            j1=mods
            j2=modt
            if(j2.gt.mode) j2=mode
            j3=modi

        if(j1.gt.j2) return
        omega=2.d+00*pi/dble(per)
        per1=pers1(1)
        per2=pers2(1)
        n1=0
        n2=0
        ii=0
  100   continue
        ii=ii+1
        call finds(1,ii,j1,j2,j11,j22,j33,n1,n2,n3,mode)
        icnt=0
        do 300 j=j11,j22
            if(ales(j,1).lt.1.e-30.and.ales(j,2).lt.1.e-30 .and.
     1      aler(j,1).lt.1.e-30.and.aler(j,2).lt.1.e-30) then
                icnt=icnt+1
            else
                icnt=0
            endif
            if(icnt.gt.3) go to 310
            if(ales(j,1).lt.1.e-30.or.ales(j,2).lt.1.e-30 .and.
     1      aler(j,1).lt.1.e-30.or.aler(j,2).lt.1.e-30) go to 300
            k=(n3-1)*NMD+j
            if(k.ge.j1.and.k.le.j2) then
                if(mod(k-j1,j3).eq.0) go to 150
            endif
            go to 300
  150       continue
            rat=0.0d+00
            rat=(dble(per)-dble(per1))/(dble(per2)-dble(per1))
            wvno=(dble(wvml(j,1)*per1) + 
     1          dble(wvml(j,2)*per2-wvml(j,1)*per1)*rat)/per
c-----
c       process this distance
c-----
            rx=dble(rr)
            tx=dble(tshft)
            if(zeta.lt.0.0) then
                fact1=1.d+00
                fact2=1.d+00
                if(fact2.eq.1.d+00) go to 160
            endif
            atn1=1.0d+00
            if(zeta.gt.2.0) then
                atn1=zeta/100.0
            else
                if(per.lt.10.0) then
                    atn1=dble(per)**dble(zeta)
                else
                    atn1=dble(10.0)**dble(zeta)
                endif
            endif
            atn2=atn1*dble(gaml(j,2))*rx
            atn1=atn1*dble(gaml(j,1))*rx
            fact1=0.0d+00
            if(atn1.lt.80.0) fact1=1./dexp(atn1)
            fact2=0.0d+00
            if(atn2.lt.80.0) fact2=1./dexp(atn2)
  160       continue
c------
c       DS component.
c------
            
            v1 = dsqrt(dble(ales(j,1))*dble(aler(j,1)))
     1          /dsqrt(dble(wvml(j,1))*rx)
            v2 = dsqrt(dble(ales(j,2))*dble(aler(j,2)))
     1          /dsqrt(dble(wvml(j,2))*rx)
            w1 = dble(dut(j,1))*v1*fact1
            w2 = dble(dut(j,2))*v2*fact2
            w1= dble(ut0(j,1))*w1
            w2= dble(ut0(j,2))*w2
            dk(1) = w1+(w2-w1)*rat
c------
c       SS component.
c------
            w1 = dble(wvmls(j,1))*dble(ut(j,1))
            w1 = w1*v1*fact1
            w2 = dble(wvmls(j,2))*dble(ut(j,2))
            w2 = w2*v2*fact2
            w1= dble(ut0(j,1))*w1
            w2= dble(ut0(j,2))*w2
            dk(2) = w1+(w2-w1)*rat
c-----
c       HF Component
c-----
            w1 = dble(ut(j,1))*v1*fact1
            w2 = dble(ut(j,2))*v2*fact2
            w1= dble(ut0(j,1))*w1
            w2= dble(ut0(j,2))*w2
            dk(3) = w1+(w2-w1)*rat
c-----
c       Introduce time shift and also phase shift from eigenfunctions
c-----
            t1=omega*tx-wvno*rx+pi4
            ct1=dcos(t1)
            st1=dsin(t1)
c-----
c       TDS - requires -i in front of eigenfunction
c-----
            cvt(1) = cvt(1) + dk(1)*xmom*cmplx(st1,-ct1)
c-----
c       TSS - eigenfunction excitation is pure real
c-----
            cvt(2) = cvt(2) + dk(2)*xmom*cmplx(ct1,st1)
c-----
c       THF - requires -i in front of eigenfunction
c-----
            cvt(3) = cvt(3) + dk(3)*xmom*cmplx(st1,-ct1)
  300   continue
        if(j2.gt.j33) go to 100
  310   continue
c-----
c       introduce source time function
c-----
        if(ms.eq.1)then
            sr=1.0d+00
            si=0.0d+00
        else
            sr=dble(sr0)
            si=dble(si0)
        endif
        do 400 j=1,3
            cvt(j) = cvt(j)*dcmplx(sr,si)
  400   continue
        do 501 j=1,3
            t1 = cdabs(cvt(j))
            if(t1 .le. 1.0d-40)cvt(j) = 0.0d+00
  501   continue
c-----
c       5 - TDS
c       8 - TSS
c       15 - THF
c-----
        xx(5) = -cvt(1)
        xx(8) = -cvt(2)
        xx(15)= -cvt(3)
        return
        end


        subroutine frstar(r,hs,hr,mname,ipsvsh,time,
     1          SSA, SSC, SSF, SSL, SSN, SSR,
     1          RRA, RRC, RRF, RRL, RRN, RRR,
     1          rayp, geom, tstar, dolock)
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

        common/earth/radius
        real radius

        common/depref/refdep
        real refdep

        integer l, lgstr
        integer iunit, iiso, iflsph, idimen, icnvel, ierr
        character title*80 
        
        radius = 6371.
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
     1      hs+refdep, hr+refdep, ipsvsh,iflsph, rayp,
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
        do 3000 i=1,10
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
        integer NL
        parameter(NL=200)
            real TA(NL), TL(NL), TN(NL), TRho(NL)
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
