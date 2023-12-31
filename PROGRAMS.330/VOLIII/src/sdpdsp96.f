c----------------------------------------------------------------------c
c                                                                      c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                 c
c      VOLUME III                                                      c
c                                                                      c
c      PROGRAM: SDPDSP96                                               c
c                                                                      c
c      COPYRIGHT 1996                                                  c
c      R. B. Herrmann, Chien-Ying Wang                                 c
c      Department of Earth and Atmospheric Sciences                    c
c      Saint Louis University                                          c
c      221 North Grand Boulevard                                       c
c      St. Louis, Missouri 63103                                       c
c      U. S. A.                                                        c
c                                                                      c
c----------------------------------------------------------------------c
c Revision history:
c       06 AMR 2001 - new program 
c           purpose plot observed and theoretical dispersion curves
c           from a surf96 file[s]
c       03 MAY 2004 - changed maximum number of observations to 
c           5000 from 1000
c       07 JUL 2004 - color flag is now -K kolor - also changed the
c           override on the maximum range for log period/frequency axis
c       07 JUL 2005 - added return to subroutine outplt
c       08 NOV 2006 - error in logic for xmin/xmax/ymin/ymax
c       09 FEB 2017 - change label from sec to s
c       04 JUN 2019 - when plotting error bars and shaded circles, plot
c                     circular outline last for solid symbol
c       09 NOV 2020 - When plotting a continuous curvedi1, The -S solsiz 
c                     sets the curve width
c----------------------------------------------------------------------c
        integer LIN, LER, LOT
        parameter (LIN=5,LOT=6,LER=0)

        common/ctrl/ tmin,tmax,cmin,cmax

c-----
c       command line arguments
c-----
        logical dofreq, dobox,  doxlog, doylog
        integer*4 ilorr, iporg
        real*4 xmin, xmax, ymin, ymax, x0, y0, xlen, ylen
        integer*4 kolor

        parameter (NDAT=100)
        character*80 datfil(NDAT)
        logical doerr(NDAT), dosolid(NDAT)
        logical doblack
        logical verby
        integer ls

c-----
c       machine dependent initialization
c-----
        call mchdep()
c-----
c       get control parameters from command line
c-----
        call gcmdln(ilorr, dofreq, dobox, xmin, xmax, ymin, ymax, 
     1      x0, y0, xlen, ylen, solsiz,doblack,verby,
     2      kolor, iporg,doxlog,datfil,doerr,dosolid,ndatfl)
        if(ilorr.lt.1)call usage('specify -L or -R')
        if(iporg.lt.1)call usage('specify -C or -U or -G')
        if(iporg.eq.3)then
            doylog = .true.
        else
            doylog = .false.
        endif
        if(verby)then
            WRITE(LER,*)'ilorr      :',ilorr
            WRITE(LER,*)'dofreq     :',dofreq
            WRITE(LER,*)'dobox      :',dobox
            WRITE(LER,*)'xmin       :',xmin
            WRITE(LER,*)'ymin       :',ymin
            WRITE(LER,*)'xmax       :',xmax
            WRITE(LER,*)'ymax       :',ymax
            WRITE(LER,*)'x0         :',x0
            WRITE(LER,*)'y0         :',y0
            WRITE(LER,*)'xlen       :',xlen
            WRITE(LER,*)'ylen       :',ylen
            WRITE(LER,*)'solsiz     :',solsiz
            WRITE(LER,*)'doblack    :',doblack
            WRITE(LER,*)'kolor      :',kolor
            WRITE(LER,*)'iporg      :',iporg
            WRITE(LER,*)'doxlog     :',doxlog
            WRITE(LER,*)'ndatfl     :',ndatfl
            WRITE(LER,*)'doerr dosolid dispersion_file'
            do 1000 i=1,ndatfl
                ls = lgstr(datfil(i))
                WRITE(LER,'(l5,l8,1x,a)')
     1              doerr(i),dosolid(i),datfil(i)(1:ls)
 1000       continue
        endif
c------
c       get control parameters and open input file.
c------
        call limt(dofreq,iporg,datfil,ndatfl)
        if(xmin .lt.0 .or. xmax .lt.0)then
            xmin = tmin
            xmax = tmax
        endif
        if(ymin .lt.0 .or. ymax .lt.0)then
            ymin = cmin
            ymax = cmax
        endif
        if(verby)then
            write(LER,*)'After examining dispersion data sets'
            WRITE(LER,*)'xmin       :',xmin
            WRITE(LER,*)'ymin       :',ymin
            WRITE(LER,*)'xmax       :',xmax
            WRITE(LER,*)'ymax       :',ymax
        endif
c-----


        call outplt(ilorr, dofreq, dobox, xmin, xmax, ymin, ymax, 
     1      x0, y0, xlen, ylen, solsiz,doblack,
     2      kolor,iporg,doxlog,doylog,datfil,doerr,dosolid,ndatfl)

        close (1)
        end

        subroutine outplt(ilorr, dofreq, dobox, xmin, xmax, ymin, ymax, 
     1      x0, y0, xlen, ylen, solsiz,doblack,
     2      kolor,iporg,doxlog,doylog,datfil,doerr,dosolid,ndatfl)
        implicit none
        integer LIN, LER, LOT
        parameter (LIN=5,LOT=6,LER=0)
c-----
c       plot the dispersion curves
c
c       ilorr   I*4 1 = Love, 2 = Rayleigh
c       dofreq  L   .true. horizontal axis is frequency
c       dobox   L   .true. put frame around plot
c       xmin    R*4 minimum value of independent variable
c       xmax    R*4 maximum value of independent variable
c       ymin    R*4 minimum value of   dependent variable
c       ymax    R*4 maximum value of   dependent variable
c       x0  R*4 
c       y0  R*4 lower left coordinates of plot box
c       xlen    R*4 length of X axis in inches
c       ylen    R*4 length of Y axis in inches
c       kolor   I*4 CALPLOT pen value
c       iporg   I*4 1 phase velocity
c               2 group velocity
c               3 gamma
c       datfil  Ch*(*)  User provided dispersion points
c       doblack L   .true. put black edge around solid symbol
c-----
        integer*4 ilorr
        logical dofreq, dobox, doxlog, doylog
        real*4 xmin, xmax, ymin, ymax
        real*4 x0, y0, xlen, ylen, solsiz
        integer*4 kolor, iporg, ndatfl
        integer NDAT
        parameter (NDAT=100)
        character*80 datfil(NDAT)
        logical doerr(NDAT), dosolid(NDAT)
        logical doblack

        real permin, permax

c-----
c       user provided data
c-----
        integer ierr, ndsp,  NOBS
        parameter (NOBS=160000)
        integer jlorr(NOBS), jobs(NOBS), jobsyn(NOBS), jmode(NOBS) 
        real*4  fper(NOBS)   , fobs(NOBS)   , fobserr(NOBS)

        integer nd, i, NPLT, ipen
        real v,yy,yyp,t,yym,xx,dv


        real tmin,tmax,cmin,cmax
        common/ctrl/ tmin,tmax,cmin,cmax

        logical inside

c------
c       plot the dispersion values.
c------
c------
c       plot the frame and find the mode range to be plotted.
c------
c       To ensure useful logarithmic plots, ensure that there
c       is at least one complete cycle
c       except for x-axis = do not want to break sacmft96 and
c       sacpom96 - so commented out 22 JUL 2004
c-----
C       if(doylog)then
C           if( ymin .gt. 0.1*ymax)ymin = 0.099*ymax
C       endif
C       if(doxlog)then
C           if( xmin .gt. 0.2*xmax)xmin = 0.099*xmax
C       endif

        call start(dobox,iporg,ilorr,x0,y0,xlen,ylen,
     1      xmin,xmax,ymin,ymax,
     1      doxlog,doylog,dofreq)
c-----
c       plot the dispersion curves
c-----
        do 1000 nd=1,ndatfl
c-----
c       plot observation points
c-----
            if(datfil(nd) .ne. ' ')then
c-----
c       get dispersion values
c-----
            call rddisp(datfil(nd),3,ierr,ndsp,NOBS,
     1          jlorr,jobs,jobsyn,jmode,
     2          fper,fobs,fobserr)
            if(dosolid(nd))then
                call newpen(1)
            else
                call newpen(kolor)
            endif
c-----
c       plot the observed dispersion values
c-----
            ipen = 3
            NPLT = 0
            permin =  1.0e+38
            permax = -1.0e+38

            do 2000 i=1,ndsp
                if(dofreq)then
                    t = 1.0/fper(i)
                else
                    t = fper(i)
                endif
                if(iporg .eq. jobs(i) .and. ilorr.eq.jlorr(i))then
                    v = fobs(i)
                    dv = fobserr(i)
                    if(.not.doerr(nd))dv = 0.0
                    if(inside(t,v,xmin,xmax,ymin,ymax))then
                    call dotran(doxlog,xmin,xmax,t,
     1                  x0,xlen,xx,
     1                  doylog,ymin,ymax,v,dv,
     1                  y0,ylen,yy,yyp,yym)
                    if(dosolid(nd))then
                        call gwidth(solsiz)
                        call plot(xx,yym,ipen)
                        ipen = 2
                    else
                        call newpen(kolor)
                        call gsolid(xx,yy,solsiz,jmode(i)+1)
                        if(doerr(nd))then
                            call plot(xx,yym,3)
                            call plot(xx,yyp,2)
                        endif
                           call newpen(1)
                        if(doblack)then
                           call gpoly(xx,yy,solsiz,jmode(i)+1) 
                        endif
                    endif
                    NPLT = NPLT + 1
                    if(fper(i).gt.permax)permax = fper(i)
                    if(fper(i).lt.permin)permin = fper(i)
                    else
c-----
c                   if clipped should continue line to
c                   clip edge but too much work now
c                   requires clipping
c-----
                    ipen = 3
                    endif
                endif
                call gwidth(0.0)
 2000       continue
            endif
 1000   continue
        call newpen(1)
        call pend()
        write(LOT,*)'NDSP:',ndsp,' NPLT:',NPLT, ' PERMIN-PERMAX:',
     1      permin,'-',permax
        return
        end
c
c-----------------------------------------------------------------
c
        subroutine start(dobox,iporg,ilorr,x0,y0,xlen,ylen,
     1      xmn,xmx,ymn,ymx,
     2      doxlog,doylog, dofreq)
        parameter (LIN=5,LOT=6,LER=0)

        logical dobox, doxlog, doylog, dofreq
        character ylabel*12, xlabel*14
c------
c       open the plotter and read in the plot control parameter.
c------
        if(dofreq)then
            xlabel = 'Frequency (Hz)'
        else
            xlabel = 'Period (s)'
        endif
        if(ilorr.eq.1)then
            if(iporg.eq.1)then
                call pinitf('SLDSPC.PLT')
                ylabel = 'C (km/s)'
            else if(iporg.eq.2)then
                call pinitf('SLDSPU.PLT')
                ylabel = 'U (km/s)'
            else if(iporg.eq.3)then
                call pinitf('SLDSPG.PLT')
                ylabel = 'Gamma (1/km)'
            endif
        else if(ilorr.eq.2)then
            if(iporg.eq.1)then
                call pinitf('SRDSPC.PLT')
                ylabel = 'C (km/s)'
            else if(iporg.eq.2)then
                call pinitf('SRDSPU.PLT')
                ylabel = 'U (km/s)'
            else if(iporg.eq.3)then
                call pinitf('SRDSPG.PLT')
                ylabel = 'Gamma (1/km)'
            endif
        endif
c-----
c       ensure that the plot limits are in increasing order
c-----
        if(xmn.gt.xmx)then
            xmin = xmx
            xmax = xmn
        else
            xmin = xmn
            xmax = xmx
        endif
        if(ymn.gt.ymx)then
            ymin = ymx
            ymax = ymn
        else
            ymin = ymn
            ymax = ymx
        endif
c-----
c       Plot the frame.
c-----
        lx = lgstr(xlabel)
        ly = lgstr(ylabel)
        if(dobox)then
            call gbox(x0+0.0,y0+0.0,x0+xlen,y0+ylen)
            if(doylog)then
                call dology(x0+0.0 ,y0+0.0,ylen,ymax,ymin,
     1              0.14,.false.,.true.,.true.,ly,
     2              ylabel)
                call dology(x0+xlen,y0+0.0,ylen,ymax,ymin,
     1              0.14,.true.,.false.,.false.,ly,
     2              ylabel)
            else
                call doliny(x0+0.0 ,y0+0.0,ylen,ymax,ymin,
     1              0.14,.false.,.true.,.true.,ly,
     2              ylabel)
                call doliny(x0+xlen,y0+0.0,ylen,ymax,ymin,
     1              0.14,.true.,.false.,.false.,ly,
     2              ylabel)
            endif
            if(doxlog)then
                call dologx(x0+0.0,y0+0.0,xlen,xmax,xmin,
     1              0.14,.true.,.false.,.true.,lx,
     2              xlabel)
                call dologx(x0+0.0,y0+ylen,xlen,xmax,xmin,
     1              0.14,.false.,.false.,.false.,lx,
     2              xlabel)
            else
                call dolinx(x0+0.0 ,y0+0.0,xlen,xmax,xmin,
     1              0.14,.true.,.false.,.true.,lx,
     2              xlabel)
                call dolinx(x0+0.0,y0+ylen,xlen,xmax,xmin,
     1              0.14,.false.,.true.,.false.,lx,
     2              xlabel)
            endif
        else
c-----
c       put in corner markers
c-----
            call plot(x0+0.14,y0+0.00,3)
            call plot(x0+0.00,y0+0.00,2)
            call plot(x0+0.00,y0+0.14,2)

            call plot(x0+xlen-0.14,y0+0.00,3)
            call plot(x0+xlen-0.00,y0+0.00,2)
            call plot(x0+xlen-0.00,y0+0.14,2)

            call plot(x0+xlen-0.14,y0+ylen-0.00,3)
            call plot(x0+xlen-0.00,y0+ylen-0.00,2)
            call plot(x0+xlen-0.00,y0+ylen-0.14,2)

            call plot(x0+0.14,y0+ylen-0.00,3)
            call plot(x0+0.00,y0+ylen-0.00,2)
            call plot(x0+0.00,y0+ylen-0.14,2)

        endif
        return
        end

        subroutine limt(dofreq,iporg,datfil,ndatfl)
c-----
c       find the range of mode numbers to be plotted.
c       and also plotting limits
c-----
        implicit none
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)
        logical dofreq
        integer iporg, ndatfl

        integer NDAT
        parameter(NDAT=100)
        character*80 datfil(NDAT)


        real tp,vp

        real tmin,tmax,cmin,cmax
        common/ctrl/ tmin,tmax,cmin,cmax

c-----
c       user provided data
c-----
        integer NOBS
        parameter (NOBS=160000)
        integer jlorr(NOBS), jobs(NOBS), jobsyn(NOBS), jmode(NOBS)  
        real*4  fper(NOBS)   , fobs(NOBS), fobserr(NOBS)

        integer i, nd, ierr, ndsp
c------
c       find tmin and tmax in which data will be plotted.
c------
        tmin=1.0e+10
        tmax=-1.0e+10
        cmin = 1.0e+10
        cmax = -1.0e+10
        do 1000 nd=1,ndatfl
            call rddisp(datfil(nd),3,ierr,ndsp,NOBS,
     1          jlorr,jobs,jobsyn,jmode,
     2          fper,fobs,fobserr)
c-----
c           now get the limits
c-----
            do 2000 i=1,ndsp
                if(iporg .eq. jobs(i))then
                    if(dofreq)then
                        tp=1.0/fper(i)
                    else
                        tp=fper(i)
                    endif
                    if(tp.gt.tmax)tmax = tp
                    if(tp.lt.tmin)tmin = tp
                    vp = fobs(i)
                    if(vp.gt.cmax)cmax = vp
                    if(vp.lt.cmin)cmin = vp
                endif
 2000       continue
        write(0,*)'tmin,tmax,cmin,cmax:',tmin,tmax,cmin,cmax
 1000   continue
        return
        end

        subroutine gcmdln(ilorr, dofreq, dobox, xmin, xmax, 
     1      cmin, cmax, 
     1      x0, y0, xlen, ylen, solsiz,doblack,verby,
     2      kolor, iporg,doxlog,datfil,doerr,dosolid,ndatfl)
        integer*4 ilorr, iporg
        logical dofreq, dobox,  doxlog
        real*4 xmin, xmax, cmin, cmax
        real*4 x0, y0, xlen, ylen
        integer*4 kolor
        integer ndatfl 
        logical doblack
        logical verby

        parameter (NDAT=100)
        character*80 datfil(NDAT)
        logical doerr(NDAT) ,dosolid(NDAT)

        character names*80
        real*4 solsiz
c-----
c       defaults
c-----
        ilorr = 0
        iporg = 0
        dofreq = .true.
        dobox = .true.
        xmin = -1.0
        xmax = -1.0
        cmin = -1.0
        cmax = -1.0
        x0 = 2.0
        y0 = 1.0
        xlen = 6.0
        ylen = 6.0
        kolor = 1
        doxlog = .false.
        solsiz = 0.03
        doblack = .true.
        verby = .false.
        
        ndatfl = 0

        nmarg = mnmarg()
        i = 0
 1000   continue
            i = i + 1
            if(i.gt.nmarg)go to 2000
            call mgtarg(i,names)
            if(names(1:2).eq.'-L' )then
                ilorr = 1
            else if(names(1:2).eq.'-R')then
                ilorr = 2
            else if(names(1:2).eq.'-C')then
                iporg = 1
            else if(names(1:2).eq.'-U')then
                iporg = 2
            else if(names(1:2).eq.'-G')then
                iporg = 3
            else if(names(1:5).eq.'-FREQ')then
                dofreq = .true. 
            else if(names(1:4).eq.'-PER')then
                dofreq = .false.    
            else if(names(1:5).eq.'-XMIN')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f20.0)')xmin
            else if(names(1:5).eq.'-XMAX')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f20.0)')xmax
            else if(names(1:5).eq.'-YMIN')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f20.0)')cmin
            else if(names(1:5).eq.'-YMAX')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f20.0)')cmax
            else if(names(1:3).eq.'-X0')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f20.0)')x0
            else if(names(1:3).eq.'-Y0')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f20.0)')y0
            else if(names(1:5).eq.'-XLEN')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f20.0)')xlen
            else if(names(1:5).eq.'-YLEN')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f20.0)')ylen
            else if(names(1:2).eq.'-K')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,i20)')kolor
            else if(names(1:2).eq.'-S')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f20.0)')solsiz
            else if(names(1:6).eq.'-NOBOX')then
                dobox = .false.
            else if(names(1:5).eq.'-XLOG')then
                doxlog = .true.
            else if(names(1:3).eq.'-DE')then
                i = i + 1
                call mgtarg(i,names)
                ndatfl = ndatfl + 1
                datfil(ndatfl) = names
                doerr(ndatfl) = .true.
                dosolid(ndatfl) = .false.
            else if(names(1:3).eq.'-DC')then
                i = i + 1
                call mgtarg(i,names)
                ndatfl = ndatfl + 1
                datfil(ndatfl) = names
                doerr(ndatfl) = .false.
                dosolid(ndatfl) = .true.
            else if(names(1:2).eq.'-D'.and.names(1:3).ne.'-DE'
     1          .and.names(1:3).ne.'-DC')then
                i = i + 1
                call mgtarg(i,names)
                ndatfl = ndatfl + 1
                datfil(ndatfl) = names
                doerr(ndatfl) = .false.
                dosolid(ndatfl) = .false.
            else if(names(1:5).eq.'-NOBL')then
                doblack = .false.
            else if(names(1:2).eq.'-V')then
                verby = .true.
            else if(names(1:2).eq.'-?')then
                call usage(' ')
            else if(names(1:2).eq.'-h')then
                call usage(' ')
            endif
        go to 1000
 2000   continue
        if(xmin .lt. 0.0 .or. xmax .lt. 0.0 .or. cmin .lt. 0.0
     1      .or. cmax .lt. 0.0)then
            call usage('XMIN/XMAX YMIN/YMAX limits required')
        endif
        return
        end

        subroutine usage(str)
        character str*(*)
        parameter (LER=0,LIN=5,LOT=6)
        write(LER,*)'sdpsrf96: ',str
        write(LER,*)'USAGE:',
     1  'sdpdsp96 ',
     1  '[-L -R] [-C -U -G] [-FREQ -PER]',
     1  '-XMIN xmin -XMAX xmax -YMIN ymin -YMAX ymax ',
     1  '-X0 x0 -Y0 y0 -K kolor -NOBOX -TXT -ASC -XLOG ',
     1  '-D dispfile -DE dispfile -DC dispfile -NOBLACK -V -? -h'

        write(LER,*)
     1  '-L                        Love waves'
        write(LER,*)
     1  '-R                        Rayleigh waves'
        write(LER,*)
     1  '     Note one of -L or -R is required'
        write(LER,*)
     1  '-C                        Phase velocity'
        write(LER,*)
     1  '-U                Group velocity'
        write(LER,*)
     1  '-G        Anelastic attenuation coefficient'
        write(LER,*)
     1  '     Note one of -C  -U or -G is required'
        write(LER,*)
     1  '-FREQ      (default true) X-Axis is frequency'
        write(LER,*)
     1  '-PER       (default false)X-Axis is period'
        write(LER,*)
     1  '-XMIN xmin (default 0.0)  minimum value of X-Axis'
        write(LER,*)
     1  '-XMAX xmax (default    )  maximum value of X-Axis'
        write(LER,*)
     1  '-YMIN ymin (default 0.0)  minimum value of Y-Axis'
        write(LER,*)
     1  '-YMAX ymax (default 0.0)  maximum value of Y-Axis'
        write(LER,*)
     1  '-X0 x0     (default 2.0)  lower left corner of plot'
        write(LER,*)
     1  '-Y0 y0     (default 1.0)  bottom left corner of plot'
        write(LER,*)
     1  '-XLEN xlen (default 6.0)  length of X-Axis'
        write(LER,*)
     1  '-YLEN ylen (default 6.0)  length of Y-Axis'
        write(LER,*)
     1  '-S solsiz  (default 0.03) size of observed symbol'
        write(LER,*)
     1  '-K kolor   (default 1  )  color for curves'
        write(LER,*)
     1  '-NOBOX     (default false) do not plot axes'
        write(LER,*)
     1  '-XLOG      (default linear) X axis is logarithmic'
        write(LER,*)
     1  '-D dispfile (default ignore) plot SURF96 dispersion values'
        write(LER,*)
     1  '-DE dispfile (default ignore) plot SURF96 dispersion/error'
        write(LER,*)
     1  '-DC dispfile (default ignore) plot SURF96',
     1        ' continuous dispersion'
        write(LER,*)
     1  '-NOBLACK    (default black) do not but black around symbol'
        write(LER,*)
     1  '-V          (default false) verby output'
        write(LER,*)
     1  '-?          (default false) online help'
        write(LER,*)
     1  '-h          (default false) online help'
        stop 
        end


        subroutine dotran(doxlog,xmin,xmax,xc1,x0,xlen,xx,
     1      doylog,ymin,ymax,yc1,ydv,y0,ylen,yy,yyp,yym)
c------
c       transform from user space to plot space
c----
c       doxlog  L   .true. X-axis is logarithmic
c       xmin    R*4 minimum value of X-axis
c       xmax    R*4 maximum value of X-axis
c       xc1 R*4 X-coordinate
c       x0  R*4 Physical position of xmin
c       xlen    R*4 Length of X-axis
c       xx  R*4 Physical position of point xc1
c       doylog  L   .true. Y-axis is logarithmic
c       ymin    R*4 minimum value of Y-axis
c       ymax    R*4 maximum value of Y-axis
c       yc1 R*4 Y-coordinate
c       ydv R*4 error in y-coordinate
c       y0  R*4 Physical position of xmin
c       ylen    R*4 Length of Y-axis
c       yy  R*4 Physical position of point yc1
c       yyp R*4 Physical position of point yc1 + ydv
c       yym R*4 Physical position of point yc1 - ydv
c-----
        logical doxlog, doylog
        real*4 xmin, xmax, xc1, yc1, x0, y0, xlen, ylen, xx, yy
                if(doylog)then
                     yy = y0 + ylen*alog10(yc1/ymin)
     1                  /alog10(ymax/ymin)
                     yyp = y0 + ylen*alog10((yc1+ydv)/ymin)
     1                  /alog10(ymax/ymin)
                     yym = y0 + ylen*alog10((yc1-ydv)/ymin)
     1                  /alog10(ymax/ymin)
                else
                     yy = y0 + ylen*(yc1 - ymin)/(ymax - ymin)
                     yyp = y0 + ylen*(yc1+ydv - ymin)/(ymax - ymin)
                     yym = y0 + ylen*(yc1-ydv - ymin)/(ymax - ymin)
                endif
                if(doxlog)then
                     xx = x0 + xlen*alog10(xc1/xmin)
     1                  /alog10(xmax/xmin)
                else
                     xx = x0 + xlen*(xc1 - xmin)/(xmax - xmin)
                endif
        return
        end

        function inside(x,y,xmin,xmax,ymin,ymax)
        logical inside
        real*4 x, y, xmin, xmax, ymin, ymax
        if(x.ge.xmin .and. x.le.xmax .and. 
     1          y.ge.ymin .and. y.le.ymax)then
            inside = .true.
        else
            inside = .false.
        endif
        return
        end

