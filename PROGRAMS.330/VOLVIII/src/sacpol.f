        program sacpol
c--------------------------------------------------------------------c 
c                                                                    c 
c     COMPUTER PROGRAMS IN SEISMOLOGY                                c 
c     VOLUME VIII                                                    c 
c                                                                    c
c     PROGRAM: SACPOL                                                c 
c                                                                    c
c     COPYRIGHT 2013                                                 c
c     R. B. Herrmann                                                 c
c     Department of Earth and Atmospheric Sciences                   c
c     Saint Louis University                                         c
c     221 North Grand Boulevard                                      c
c     St. Louis, Missouri 63103                                      c
c     U. S. A.                                                       c
c                                                                    c
c--------------------------------------------------------------------c
c      CHANGES
c      07 JUN 2013 - created this program to plot particle motions
c                  from sac files
c      02 MAY 2014 - corrected errors in call to getsac - declaration
c                  for the second trace were missing Larry Baker, USGS 
c      14 NOV 2016 - -ARROW for indicating increasing time
c      24 JUN 2017 - related arrow length to the length of the axes
c      22 FEB 2018 - correct labeling error for -X1RT
c                    draw 4 arrows and plot at end of segment
c      06 DEC 2021 - corrected Usage hint
c-----
        integer LER
        parameter (LER=0)
c-----
c       command line parameters
c-----
        real x0, y0, xlen, ylen
        integer kolor, x1axisdir, x2axisdir
        logical  doarrow
        character*80 lab1, lab2
        integer MXFILE
        parameter (MXFILE=2)
        character*80 fname(MXFILE)

c-----
c       trace variables
c-----
        integer MXPTS
        parameter (MXPTS = 64000)
        real seis1(MXPTS), seis2(MXPTS)
        character sta1*8, comp1*8, cdate1*12
        integer npts1, npts2
        integer nzyear1, nzjday1, nzhour1, nzmin1, mzsec1
        real dist1, deg1, az1, baz1, t01, dt1a,
     1     tp1, ts1, user91, user51,beg1
        character sta2*8, comp2*8, cdate2*12
        integer nzyear2, nzjday2, nzhour2, nzmin2, mzsec2
        real dist2, deg2, az2, baz2, t02, dt2a,
     1     tp2, ts2, user92, user52,beg2
c-----
c       plot variables
c-----
        real x(MXPTS), y(MXPTS)



c-----
c      get plot controls
c-----
        call gcmdln(x0,y0,xlen,ylen,
     1      kolor,fname,MXFILE,x1axisdir,x2axisdir,
     2      lab1, lab2,doarrow)
        
        
        call pinitf('SACPOL.PLT')
        call newpen(1)
c-----
c      open the sacfiles
c-----
       WRITE(6,*)'1: ',fname(1)
       WRITE(6,*)'2: ',fname(2)
            call getsac(fname(1),npts1,dist1,deg1,az1,baz1,t01,dt1,sta1,
     1              comp1,cdate1,seis1,tp1,ts1,user91,user51,
     2              nzyear1,nzjday1,nzhour1,nzmin1,nzsec1,beg1)
            call getsac(fname(2),npts2,dist2,deg2,az2,baz2,t02,dt2,sta2,
     2              comp2,cdate2,seis2,tp2,ts2,user92,user52,
     2              nzyear2,nzjday2,nzhour2,nzmin2,nzsec2,beg2)
            if(npts1.ne.npts2)then
                 WRITE(LER,*)'Unequal number samples'
                 call usage(' ')
            endif
            if(dt1.ne.dt1)then
                 WRITE(LER,*)'Unequal sample intervals'
                 call usage(' ')
            endif
c-----
c      now plot the traces
c-----
        call plotxy(seis1,seis2,npts1,x0,y0,xlen,ylen,
     1      kolor, x1axisdir,x2axisdir,
     2      lab1, lab2,x,y,doarrow)
        call pend()
        end

        subroutine plotxy(seis1,seis2,npts,x0,y0,xlen,ylen,
     1      kolor, x1axisdir,x2axisdir,
     2      lab1, lab2,x,y,doarrow)
c-----
c      seis1     R(*)  array to plot on x1 axis
c      seis2     R(*)  array to plot on x2 axis
c      npts      I     number of trace points to plot
c      x0        R   - lower left corner of plot frame
c      y0        R   - lower left corner of plot frame
c      xlen      R   - width  of plot frame on page
c      ylen      R   - height of plot frame on page
c      kolor     I   - color of trace
c      x1axisdir I     direction of positive x1 axis
c                      0 down, 1 up, 2 right, 3 left
c      x2axisdir I     direction of positive x2 axis
c                      0 down, 1 up, 2 right, 3 left
c      lab1      C*80  label for x1 axis
c      lab2      C*80  label for x2 axis
c      x         R(*)  array for final x-axis plot
c      y         R(*)  array for final y-axis plot
c-----
        real seis1(*), seis2(*)
        integer npts
        real x0, y0, xlen, ylen
        integer kolor
        integer x1axisdir,x2axisdir
        logical doarrow
        character lab1*(*), lab2*(*)
        real x(*), y(*)

        real vmax
        real xx, yy
        real xc, yc

c-----
c       the bounding box is [x0,y0], [x0+xlen,y0+ylen]
c       define the center
c-----
        xc = x0 + 0.5*xlen 
        yc = y0 + 0.5*ylen 
c-----
c       put in reference lines
c-----
        call plot (xc-0.5*xlen,yc      ,3)
        call plotd(xc+0.5*xlen,yc,7.0,0.01*xlen)
        call plot (xc,yc-0.5*ylen      ,3)
        call plotd(xc,yc+0.5*ylen,7.0,0.01*ylen)
c-----
c       put in axis labels
c-----
        call putaxlab(x0,y0,xc,yc,xlen,ylen,
     1    x1axisdir,lab1)
        call putaxlab(x0,y0,xc,yc,xlen,ylen,
     1    x2axisdir,lab2)
c-----
c       get the extreme values
c       eventually to will be overridden by a command line
c       value
c-----
        vmax = 0.0
c-----
c       put arrows at npts/4 npts/2 npts/4
        do i=1,npts
           if(abs(seis1(i)) .gt.vmax)vmax = abs(seis1(i))
           if(abs(seis2(i)) .gt.vmax)vmax = abs(seis2(i))
        enddo
c-----
c       safety
c-----
        if(vmax.eq.0.0)vmax = 1.0
c-----
c       fill and scale the x and y arrays 
c       x and y are scaled to be [-1, 1]
c-----
c-----
c      x1axisdir I     direction of positive x1 axis
c                      0 down, 1 up, 2 right, 3 left
c      x2axisdir I     direction of positive x2 axis
c                      0 down, 1 up, 2 right, 3 left
c-----
        do i=1,npts
           if(x1axisdir.eq.0)then
              y(i) = - seis1(i)/vmax
           else if(x1axisdir.eq.1)then
              y(i) =   seis1(i)/vmax
           else if(x1axisdir.eq.2)then
              x(i) =   seis1(i)/vmax
           else if(x1axisdir.eq.3)then
              x(i) = - seis1(i)/vmax
           endif
           if(x2axisdir.eq.0)then
              y(i) = - seis2(i)/vmax
           else if(x2axisdir.eq.1)then
              y(i) =   seis2(i)/vmax
           else if(x2axisdir.eq.2)then
              x(i) =   seis2(i)/vmax
           else if(x2axisdir.eq.3)then
              x(i) = - seis2(i)/vmax
           endif
        enddo
c-----
c       now plot these but do not use the entire
c       frame window since I need space for the
c       labels and arrows
c-----
        ipen = 3
        call newpen(kolor)
        do i=1,npts
           xx = xc + 0.4*xlen*x(i)
           yy = yc + 0.4*ylen*y(i)
           call plot(xx,yy,ipen)
           ipen = 2
        enddo
        call plot(xx,yy,3)
        call newpen(1)
c-----
c       plot 3 arrows at (1/4) (2/4) (3/4) and (4/4) of length of trace
c       the arrow goes from (xx1,yy1) to (xx2,yy2)
c       the the arrow size is keyed to the xlen ylen
c-----
        if(doarrow .and. npts.ge.8)then
           xx0 = xc + 0.4*xlen*x(npts/4 -1)
           yy0 = yc + 0.4*ylen*y(npts/4 -1)
           xx1 = xc + 0.4*xlen*x(npts/4)
           yy1 = yc + 0.4*ylen*y(npts/4)
           call makearrow(xx0,yy0,xx1,yy1,.false.,xlen,ylen)
           xx0 = xc + 0.4*xlen*x(2*npts/4 -1)
           yy0 = yc + 0.4*ylen*y(2*npts/4 -1)
           xx1 = xc + 0.4*xlen*x(2*npts/4)
           yy1 = yc + 0.4*ylen*y(2*npts/4)
           call makearrow(xx0,yy0,xx1,yy1,.false.,xlen,ylen)
           xx0 = xc + 0.4*xlen*x(3*npts/4 -1)
           yy0 = yc + 0.4*ylen*y(3*npts/4 -1)
           xx1 = xc + 0.4*xlen*x(3*npts/4)
           yy1 = yc + 0.4*ylen*y(3*npts/4)
           call makearrow(xx0,yy0,xx1,yy1,.false.,xlen,ylen)
           xx0 = xc + 0.4*xlen*x(4*npts/4 -1)
           yy0 = yc + 0.4*ylen*y(4*npts/4 -1)
           xx1 = xc + 0.4*xlen*x(4*npts/4)
           yy1 = yc + 0.4*ylen*y(4*npts/4)
           call makearrow(xx0,yy0,xx1,yy1,.false.,xlen,ylen)
        endif
        call newpen(1)
        return
        end
 
        subroutine makearrow(x1,y1,x2,y2,dobold,xlen,ylen)
c-----
c       arrow starts at (x1,y1) in direction of (x2,y2)
c       when calling subroutine arrow, start at (x1,y1) but
c       make length proportional to xlen,ylen
c       arrow is plotted in 1-2 direction but at (x2,y2)
c-----
        real x1,y1,x2,y2,xlen,ylen
        logical dobold
        r = sqrt( (x1-x2)**2 + (y1-y2)**2)
        if(r.lt.1.0e-6)return
        c = (x2-x1)/r
        s = (y2-y1)/r
        ht = 0.05 * xlen
        xx1 = x2 - ht*c
        yy1 = y2 - ht*s
        xx2 = x2 
        yy2 = y2 
        call arrow(xx1,yy1,xx2,yy2,dobold)
        return
        end


        subroutine getsac(name,npts,dist,deg,az,baz,t0,dt,sta,
     1                 comp,cdate,seis,tp,ts,user9,user5,
     2          nzyear,nzjday,nzhour,nzmin,nzsec,beg)
c-----
c
c      name    - file name to write
c      n   - number of points in FFT must be power of 2
c      n21 - number of frequencies = n/2 + 1
c      npts    - number of points in original time series
c          - which may have been zero filled to make power of 2
c      dist    - epicentral distance in km
c      deg - epicentral distance in degrees
c      az  - source - receiver azimuth in degrees
c      baz - receiver-source back azimuth
c      t0  - time of first sample after origin
c      dt  - sampling interval
c      sta - C*4 station name string
c      comp    - C*4 component name string
c      cdate   - C*12 date string
c      z   - COMPLEX array of spectra
c          nzyear,nzjday,nzhour,nzmin,nzsec)
c-----
        integer MXPTS
        parameter (MXPTS = 64000)
        integer  npts
        real dist, deg, az, baz, t0, dt
        character sta*8, comp*8, cdate*12
        character kstnm*8, kcmpnm*8
        character name*(*)
        real seis(MXPTS)
        integer*4  nzyear,nzjday,nzhour,nzmin,nzsec
*
        real evla, evlo, evdp, stla, stlo, stel, beg, origtime

        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*4 chdr(48)
*
        
        call brsac (3,MXPTS,name,seis,ierr)
*
        call getfhv('AZ      ',az,nerr)
        call getfhv('BAZ     ',baz,nerr)
        call getfhv('DIST    ',dist,nerr)
        call getfhv('GCARC   ',deg,nerr)
        call getfhv('DELTA   ', dt, nerr)
        call getnhv('NPTS    ',npts,nerr)
        call getfhv('B       ', beg, nerr)
        call getfhv('O       ',origtime,nerr)
        call getfhv('A       ',tp,nerr)
        call getfhv('T0      ',ts,nerr)
        call getfhv('USER5   ',user5,nerr)
        call getfhv('USER9   ',user9,nerr)
        call getfhv('EVLA    ',evla,nerr)
        call getfhv('EVLO    ',evlo,nerr)
        call getfhv('STLA    ',stla,nerr)
        call getfhv('STLO    ',stlo,nerr)
        call getkhv('KSTNM   ',kstnm,nerr)
        call getkhv('KCMPNM  ',kcmpnm,nerr)
        call getnhv('NZYEAR', nzyear,nerr)
        call getnhv('NZJDAY', nzjday,nerr)
        call getnhv('NZHOUR', nzhour,nerr)
        call getnhv('NZMIN', nzmin,nerr)
        call getnhv('NZSEC', nzsec,nerr)
C       write(6,*)'az,baz,dist,deg,dt,npts,b,o,nerr,tp,ts',
C     1     az,baz,dist,deg,dt,npts,beg,origtime,nerr,tp,ts
C       write(6,*)'evla,evlo,stla,stlo,kstnm,kcmpnm',
C     1     evla,evlo,stla,stlo,kstnm,kcmpnm
*
        if(nerr .eq. 0 .and. origtime .ne. -12345)then
            t0 = beg - origtime
        else
            t0 = beg
        end if
        ttp = tp - origtime
*
C       write(6,*)'ttp=',ttp
C       write(6,*)'dist=',dist
        sta = kstnm(1:8)
        comp = kcmpnm(1:8)
C       write(6,*)'name,npts,dist,deg,az,baz,t0,dt,sta',
C     1     name,npts,dist,deg,az,baz,t0,dt,sta
        cdate = ' '
*
*
        return
        end

        subroutine dmxmn(seis,n1,n2,ampmx,ampmn)
        real seis(*)
        integer n1, n2
        ampmx = -1.0e+38
        ampmn =  1.0e+38
        do 1000 i=n1,n2
            if(seis(i).gt.ampmx)ampmx=seis(i)
            if(seis(i).lt.ampmn)ampmn=seis(i)
 1000   continue
        return 
        end

        subroutine gcmdln(x0,y0,xlen,ylen,
     1      kolor,fname,MXFILE,x1axisdir,x2axisdir,
     2      lab1, lab2,doarrow)
c-----
c      parse the command line arguments
c-----
c      x0        R   - lower left corner of plot frame
c      y0        R   - lower left corner of plot frame
c      xlen      R   - width  of plot frame on page
c      ylen      R   - height of plot frame on page
c      kolor     I   - color of trace
c      fname     C*80    - character array of SAC file names
c      MXFILE    I   - maximum number of SAC files permitted
c      x1axisdir I   - direction of positive x1 axis
c                      0 down, 1 up, 2 right, 3 left
c      x2axisdir I   - direction of positive x2 axis
c                      0 down, 1 up, 2 right, 3 left
c      lab1      C*80  label for x1 axis
c      lab2      C*80  label for x2 axis
c      doarrow   L   - 0 none 1 plot arrows
c-----
        real x0, y0, xlen, ylen
        integer kolor, x1axisdir, x2axisdir
        logical  doarrow
        character lab1*(*), lab2*(*)
        integer MXFILE
        character*80 fname(MXFILE)

        integer mnmarg
        character name*80
        integer i

        x0 = 2.0
        y0 = 2.0
        xlen = 2.0
        ylen = 2.0
        kolor = 1
        x1axisdir = 2
        x2axisdir = 1
        lab1 = ' '
        lab2 = ' '
        fname(1) = ' '
        fname(2) = ' '
        doarrow = .false.
        
        nmarg = mnmarg()
        i = 0
 1000   continue
            i = i + 1
            if(i.gt.nmarg)go to 2000
            call mgtarg(i,name)
            if(name(1:2).eq.'-K')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(i10)')kolor
            else if(name(1:3).eq.'-X0')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')x0
            else if(name(1:3).eq.'-Y0')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')y0
            else if(name(1:5).eq.'-XLEN')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')xlen
            else if(name(1:5).eq.'-YLEN')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')ylen
            else if(name(1:3).eq.'-F1')then
                i = i + 1
                call mgtarg(i,fname(1))
            else if(name(1:3).eq.'-F2')then
                i = i + 1
                call mgtarg(i,fname(2))
            else if(name(1:5).eq.'-LAB1')then
                i = i + 1
                call mgtarg(i,lab1)
            else if(name(1:5).eq.'-LAB2')then
                i = i + 1
                call mgtarg(i,lab2)
            else if(name(1:5).eq.'-X1DN')then
                x1axisdir = 0
            else if(name(1:5).eq.'-X1UP')then
                x1axisdir = 1
            else if(name(1:5).eq.'-X1RT')then
                x1axisdir = 2
            else if(name(1:5).eq.'-X1LT')then
                x1axisdir = 3
            else if(name(1:5).eq.'-X2DN')then
                x2axisdir = 0
            else if(name(1:5).eq.'-X2UP')then
                x2axisdir = 1
            else if(name(1:5).eq.'-X2RT')then
                x2axisdir = 2
            else if(name(1:5).eq.'-X2LT')then
                x2axisdir = 3
            else if(name(1:6).eq.'-ARROW')then
                doarrow = .true.
            else if(name(1:2) .eq. '-?')then
                call usage(' ')
            else if(name(1:2) .eq. '-h')then
                call usage(' ')
            endif
        go to 1000
 2000   continue
c-----
c      safety
c----
        if(fname(1).eq.' ' .or. fname(2).eq.' ')then
            call usage(' ')
        endif

        return
        end

        subroutine usage(emsg)
        implicit none
        character emsg*(*)
        integer LER
        parameter (LER=0)
        integer lgstr
        integer ls
        ls = lgstr(emsg)
        if(emsg .ne. ' ')write(LER,*)emsg(1:ls)
        write(LER,*)'Usage: sacpol -XLEN xlen -YLEN ylen',
     1      ' -X0 x0 -Y0 y0 ',
     2      ' -F1 sacfile1 ',
     2      ' -F2 sacfile2 ',
     2      ' -X1DN -X1UP -X1RT -X1LT ',
     2      ' -X2DN -X2UP -X2RT -X2LT ',
     1      ' -LAB1 lab1 -LAB2 lab2 ',
     1      ' -ARROW ',
     3      ' -K IPEN -DOAMP  -h -?'
        write(LER,*)'Plot particle motions '
        write(LER,*)
     1  '-XLEN xlen (default 6.0  ) Length X-axis'
        write(LER,*)
     1  '-YLEN ylen (default 6.0  ) Length Y-axis'
        write(LER,*)
     1  '-X0 x0     (default  2.0 )  x-position of lower left corner'
        write(LER,*)
     1  '-Y0 y0     (default  7.0 )  y-position of lower left corner'
        write(LER,*)
     1  '-K PEN     (default  1)  Use color for '
        write(LER,*)
     1  '           if kolor < 0 use red->blue progression'
        write(LER,*)
     1  '-F1 sacfile1 (mandatory) plotted on x1 axis'
        write(LER,*)
     1  '-F2 sacfile2 (mandatory) plotted on x2 axis'
        write(LER,*)
     1  '    Controls to define axis directions'
        write(LER,*)
     1  '-X1DN  (default false) x1 positive down'
        write(LER,*)
     1  '-X1UP  (default false) x1 positive up'
        write(LER,*)
     1  '-X1RT  (default true ) x1 positive right'
        write(LER,*)
     1  '-X1LT  (default false) x1 positive left'
        write(LER,*)
     1  '-X2DN  (default false) x2 positive down'
        write(LER,*)
     1  '-X2UP  (default true ) x2 positive up'
        write(LER,*)
     1  '-X2RT  (default false) x2 positive right'
        write(LER,*)
     1  '-X2LT  (default false) x2 positive left'
        write(LER,*)
     1  '-LAB1 lab1 (default blank) quote delimited label for 1 axis'
        write(LER,*)
     1  '-LAB2 lab2 (default blank) quote delimited label for 2 axis'
        write(LER,*)
     1  '-ARROW  (default false) plot arrow indicating positime time'
        write(LER,*)
     1  '-?        (default none )  this help message '
        write(LER,*)
     1  '-h        (default none )  this help message '
        stop
        end
        
        subroutine putaxlab(x0,y0,xc,yc,xlen,ylen,
     1    axisdir,lab)
        real x0,y0,xc,yc,xlen,ylen
        integer axisdir
        character lab*(*)


        integer lgstr, ll
        real ht, xpos, ypos

        ht = 0.05 * xlen
c-----
c          since the axis is labeled, put in an arrow
c-----
        if(lab.ne. ' ')then
           ll = lgstr(lab)
           if(axisdir.eq.0)then
               xpos = xc - 0.5*ll*ht
               ypos = y0 -1.5*ht
               call arrow(xc,y0+1.*ht,xc,y0,.false.)
           else if(axisdir.eq.1)then
               xpos = xc - 0.5*ll*ht
               ypos = y0 + ylen + 0.5*ht
               call arrow(xc,y0+ylen -1.*ht,xc,y0+ylen,.false.)
           else if(axisdir.eq.2)then
               xpos = x0 + xlen +0.5*ht
               ypos = yc -0.5*ht
               call arrow(x0+xlen-1.*ht,yc,x0+xlen,yc,.false.)
           else if(axisdir.eq.3)then
               xpos = x0  -0.5*ht - ll*ht
               ypos = yc -0.5*ht
               call arrow(x0+1.*ht,yc,x0,yc,.false.)
           endif
c-----
c          special case
c-----
           if(lab.eq.'X1' .or. lab.eq.'X2' .or. lab.eq.'X3' .or.
     1      lab.eq.'UZ' .or. lab.eq.'UR' .or.  lab.eq.'UT' .or.
     1      lab.eq.'UH' .or. lab.eq.'UE' )then
              call gsubsc(xpos,ypos,ht,lab(1:1),1,lab(2:2),1)
           else
              call gleft(xpos,ypos,ht,lab(1:ll),0.0)

           endif
        endif
        return
        end

c-----
c     arrow routine from PROGRAMS.330/CALPLOT/Utility/calplt.f
c-----

        subroutine arrow(xx,yy,xxx,yyy,bold)
c-----
c       note height of arrow is sqrt[xx-xxx)^2 + (yy-yyy)^2]
        logical bold
        real x(3), y(3)
        if(bold)then
            call gwidth(0.06)
        else
            call gwidth(0.02)
        endif
        call plot(xx,yy,3)
        call plot(xxx,yyy,2)
        call plot(xx,yy,2)
c-----
c       get direction cosines for arrow tip
c-----
        dx = xxx - xx
        dy = yyy - yy
        r = sqrt(dx**2 + dy**2)
        if(r .le. 1.0e-5)return
        c = dx/r
        s = dy/r
        r = 0.5*r
C       dx = -0.2*c - 0.1*s
C       dy = -0.2*s + 0.1*c
        dx = r*(-2*c - 1*s)
        dy = r*(-2*s + 1*c)
        x(1) = xxx + dx
        y(1) = yyy + dy
        call plot(xxx+dx,yyy+dy,3)
C       dx = -0.2*c + 0.1*s
C       dy = -0.2*s - 0.1*c
        dx = r*(-2*c + 1*s)
        dy = r*(-2*s - 1*c)
        x(2) = xxx + dx
        y(2) = yyy + dy
C        call plot(xxx, yyy, 2)
C        call plot(xxx+dx,yyy+dy,3)
C        call plot(xxx, yyy, 2)
C        call plot(xxx, yyy, 3)
        x(3) = xxx
        y(3) = yyy
        call shadep(3,x,y)
        if(bold)call gwidth(0.02)

        return
        end
