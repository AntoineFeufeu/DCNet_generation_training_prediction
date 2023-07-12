c----------------------------------------------------------------------c
c                                                                    c
c    COMPUTER PROGRAMS IN SEISMOLOGY                                 c
c    VOLUME X                                                        c
c                                                                    c
c    PROGRAM: TTINVPV96                                               c
c                                                                    c
c    COPYRIGHT 2011                                                  c
c    R. B. Herrmann                                                  c
c    Department of Earth and Atmospheric Sciences                    c
c    Saint Louis University                                          c
c    221 North Grand Boulevard                                       c
c    St. Louis, Missouri 63103                                       c
c    U. S. A.                                                        c
c                                                                    c
c----------------------------------------------------------------------c
        program ttinvpv96
c-----
c     CHANGES
c     27 JAN 2011 - created a program to plot the travel times
c        from ttinv96 and also ttsrf96
c        The command line will say whether to plot P or S
c     15 MAR 2011 - forced the origin to be plotted - note
c        cannot do this for dovred
c-----
        implicit none
c-----
c     LIN - unit for FORTRAN read from terminal
c     LOT - unit for FORTRAN write to terminal
c     LER - unit for FORTRAN error output to terminal
c-----
        integer LER, LIN, LOT
        parameter(LER=0,LIN=5,LOT=6)
c-----
c     NL  - number of layers in model
c     NL2 - number of columns in model (first NL2/2 are
c         - velocity parameters, second NL2/2 are Q values)
c-----
        integer NL, NL2
        parameter(NL=200,NL2=NL+NL)
        integer nf10(NL)

        integer nf
        logical dop
        real vred, wid
        character kvreds*20

        integer iprog, itot, nf1,nf2,nf34,nf5,nf67,nfilt,
     1      nup, invcsl, lstinv
        real dlam,qaqb,wref
        integer iter, nurftn, idtwo,idum2,idum3,idum4,idum5
        real twnmin,twnmax,pval, sigv, sigtp, sigts
        real rdum1,rdum2,rdum3,rdum4,rdum5

        integer nid
c-----
c     machine dependent initialization
c-----
        call mchdep()
c-----
c     get the command line 
c-----
       call gcmdln(dop,vred,kvreds,wid)
c-----
c-----
c     iprog is a binary OR: 2= rftn, 1=surf
c-----
        call gttmp0(iprog,itot,nf1,nf2,nf34,nf5,nf67,nf10,nfilt,
     1      nup, dlam,qaqb,wref,invcsl,
     2      lstinv,twnmin,twnmax,iter,nurftn,pval,
     3      sigv,sigtp,sigts,
     4      idtwo,idum2,idum3,idum4,idum5,
     5      rdum1,rdum2,rdum3,rdum4,rdum5)
c-----
c     iprog   = inversion type. Logical or of  rftn 2, surf 1 - if does
c             not match, terminate this run
c     itot    = total number of iterations [37]
c     nf1     = 1 estimated stdev computed from residuals
c               0 no scaling by residuals
c     nf2     = TOTAL number of Love Wave Gamma Modes
c               0 DO NOT PROCESS Love Wave Gamma Data for Q
c     nf34    = TOTAL number of Love Wave modes to process 
c          (from C and U)
c     nf5     = TOTAL number of Rayleigh Wave Gamma Modes
c               0 DO NOT PROCESS Rayleigh Wave Gamma Data for Q
c     nf67    = TOTAL number of Rayleigh Wave modes to process 
c          (from C and U)
c     nf10    = Input Format (from model file)
c               0 - Inversion a, rho fixed
c               1 - Inversion Poisson Ratio Fixed, Rho computed from Vp
c     nfilt   = smoothing parameter 
c               0  No model Weight  No smoothing
c               1  Model Weight     No smoothing
c               2  No model weight  Smoothing
c               3  Model weight     Smoothing
c     nup     = state
c             =1 partials computed 
c             =0 surfinv run =1 before
c             =2 on update of invcsl>=1 or 3
c     dlam    = damping [32]
c     qaqb    = 2.25 default  [34]
c     wref    = reference frequency 1.0 is default [33]
c     invcsl  = 0 acausal, 1 uncoupled causel, 2 coupled causal [35]
c     invdep  = 0 last inversion was for depth
c             = 1 last inversion was for velocity and Q inverse
c     lstinv  = 2,3,4,5 depending on the last inversion
c             invdep = 1 for 2,3,4 and 0 for 5
c     twnmin  These give the receiver function window
c     twnmax
c     iter    Current iteration
c     nurftn  Number of receiver functions to be read
c-----
CRBH        if((iprog/2).ne.1)then
CRBH            WRITE(LOT,*)'rftnpv96 requires a receiver ',
CRBH     1         'function inversion'
CRBH            STOP
CRBH        endif
            

c-----
c-----
c     get solution
c-----
      if(dop)then
            call pinitf('TTINVPVP.PLT')
      else
            call pinitf('TTINVPVS.PLT')
      endif
c-----
c     now plot the velocity model(s) for reference
c     nid = 0 forces velocity model plot
c-----
        nid = 0
        call pltmod(nid,dop)
c-----
c     Now systematically plot the travel times
c     by opening the partial derivative file, selecting and
c     then sorting
c-----
        call plttcv(dop,vred,kvreds,wid)
        call newpen(1)
        call pend()
        end

        subroutine plttcv(dop,vred,kvreds,wid)
c-----
c       dop    L - .true. plot P velocity
c                  .false. plot S velocity
c       vred   R - if > 0 use reducted travel time
c       kvreds C*20 - string giving velocity for plot
c-----
        logical dop
        real vred, wid
        character kvreds*(*)
c-----
c       open the file containing the partial derivatives
c       obervations and predictions
c-----
c-----
c     NL  - number of layers in model
c     NL2 - number of columns in model (first NL2/2 are
c         - velocity parameters, second NL2/2 are Q values)
c-----
        integer NL, NL2
        parameter(NL=200,NL2=NL+NL)
        real dtd(NL)
c-----
c       define a large array for arrival times and distances
c-----
        integer NOBSMX
        parameter (NOBSMX=4000)
        real tmp(NOBSMX), r(NOBSMX), tobs(NOBSMX), tpred(NOBSMX),
     1    errobs(NOBSMX)
        integer itmp(NOBSMX)

        integer nobs
c-----
c       since the model is read in first
c       get the number of layers in the model
c-----
        integer mmax
        common/modlly/mmax

        logical doplot

        real evdp, stdp

c-----
c     assume that tmpsrfi.21 exists
c-----
c-----
c           open binary file for the partial derivatives
c           open tmpsrfi.21
c           output for observation
c               TYPE DTDV OBS PRED ERR EVDP STDP RES
c           or
c               TYPE DTDH OBS PRED ERR EVDP STDP RES
c-----
            open(4,file='tmpsrfi.21',status='unknown',
     1          form='unformatted',access='sequential')
            rewind 4
            xmin = 1.0e+38
            xmax = -xmin
            tmin = 1.0e+38
            tmax = -tmin
            nobs = 0
 1000       continue
c-----
c           initial scan to get extreme values for the plot
c-----
            read(4,end=1001)iphase,dist,(DTD(i),i=1,mmax),
     1          time,TT,err,evdp,stdp,res
CRBH            WRITE(6,*)iphase,dist,(DTD(i),i=1,mmax),
CRBH     1          time,TT,err,evdp,stdp,res
            if(dop .and. iphase .eq.0)then
                  if(vred.gt.0.0)then
                     time = time - dist/vred
                     TT   = TT   - dist/vred
                  endif
                  nobs = nobs + 1
                  r(nobs) = dist
                  tobs(nobs) = time
                  tpred(nobs) = TT
                  errobs(nobs) = err     
                  if(dist.gt.xmax)xmax=dist
                  if(dist.lt.xmin)xmin=dist
                  if(time.gt.tmax)tmax=time
                  if(time.lt.tmin)tmin=time
                  if(TT  .gt.tmax)tmax=TT  
                  if(TT  .lt.tmin)tmin=TT  
                  
            endif
            if(.not. dop .and. iphase .gt.0)then
                  if(vred.gt.0.0)then
                     time = time - dist/vred
                     TT   = TT   - dist/vred
                  endif
                  nobs = nobs + 1
                  r(nobs) = dist
                  tobs(nobs) = time
                  tpred(nobs) = TT
                  errobs(nobs) = err     
                  if(dist.gt.xmax)xmax=dist
                  if(dist.lt.xmin)xmin=dist
                  if(time.gt.tmax)tmax=time
                  if(time.lt.tmin)tmin=time
                  if(TT  .gt.tmax)tmax=TT  
                  if(TT  .lt.tmin)tmin=TT  
            endif

        go to 1000
 1001   continue
c-----
        doplot = .false.
        if(xmin .ge. 0.0 .and. xmax .ge.0)then
              doplot = .true.
        endif
c-----
c       for a nice looking plot adjust the tmin/tmax 
c-----
        if(vred.gt.0.0)then
           tmin = tmin - 0.1*abs(tmin)
           tmax = tmax + 0.1*abs(tmax)
        else
           tmax = tmax + 0.1*abs(tmax)
           TMIN = 0.0
           XMIN = 0.0
        endif
c-----
c       now plot
c-----
        if (doplot)then
           x0 = 3.0
           y0 = 2.0
           xlen = 6.0
           ylen = 4.0
           call newpen(1)
           call doframe(x0,y0,xlen,ylen,tmin,tmax,xmin,xmax,
     1      vred,kvreds)
c-----
c          plot the predicted first arrival curves
c          since there is not certainty that the values are ordered with
c          increasing distance, we must sort, or rather set up some
c          pointers to the correct order
c-----
           do i=1,nobs
              tmp(i)  = r(i)
           enddo
           call sort(tmp,itmp,nobs)
           call pltsort(x0,y0,xlen,ylen,xmin,xmax,tmin,tmax,wid,
     1          nobs,r,tpred,itmp)
c-----
c          plot the observed data
c-----
           call dopltobs(x0,y0,xlen,ylen,xmin,xmax,tmin,tmax,
     1      r,tobs,errobs,nobs)
        endif
        close(4)

c-----
c       annotate the plot
c-----
        if(dop)then
           call gcent(x0+0.5*xlen, y0+ylen + 0.15,0.1,
     1   'P-wave first arrival time',0.0)
        else
           call gcent(x0+0.5*xlen, y0+ylen + 0.15,0.1,
     1   'S-wave first arrival time',0.0)
        endif
        call newpen(1)
        call plot(x0+0.25, y0+ylen - 0.30,3)
        call plot(x0+0.75, y0+ylen - 0.30,2)
        call gleft(x0+0.85,y0+ylen - 0.30 -0.5*0.10, 0.10,
     1   'Predicted', 0.0)
        call newpen(2)
        call fillit('CI',0.05,x0+2.0, y0+ylen - 0.30)
        call newpen(1)
        call curvit('CI',0.05,x0+2.0, y0+ylen - 0.30)
        call gleft(x0+2.20,y0+ylen - 0.30 -0.5*0.10, 0.10,
     1   'Observed', 0.0)
        call newpen(1)
        return
        end

       subroutine sort(x,key,n)
c-----
c    Starting with x(1) ,,, x(n)
c    return   the xarray sorted in increasing order
c    also return the pointers key to the initial array. 
c    For example given x = [ 3, 1, 2 ]
c    the returned values are
c                    x = [ 1, 2, 3 ]        
c                  key = [ 2, 3, 1 ]
c-----
c    Reference: http://en.wikipedia.org/wiki/Bubble_sort
c-----
       integer n
       real x(n)
       integer key(n)

       do i=1,n
           key(i) = i
       enddo
       do i = n, 1, -1
           do j = 1 , i -1
               if(x(j) .gt. x(j+1))then
                   tmp = x(j)
                   x(j) = x(j+1)
                   x(j+1) = tmp
                   ktmp = key(j)
                   key(j) = key(j+1)
                   key(j+1) = ktmp
               endif
           enddo
       enddo
       return
       end


        subroutine dopltobs(x0,y0,xlen,ylen,xmin,xmax,tmin,tmax,
     1      r,tobs,errobs,nobs)
        real x0,y0,xlen,ylen,xmin,xmax,tmin,tmax
        integer nobs
        real r(nobs), tobs(nobs), errobs(nobs)

        real dx, dy, xx, yy, yp, ym
        integer i
        dx =  xlen/(xmax-xmin)
        dy =  ylen/(tmax-tmin)
        do i = 1,nobs
                 xx = x0 + ( r(i) - xmin)*dx
                 yy = y0 + ( tobs(i) - tmin)*dy
                 yp = y0 + ( tobs(i) + errobs(i) - tmin)*dy
                 ym = y0 + ( tobs(i) - errobs(i) - tmin)*dy
                 call plot(xx,yp,3)
                 call plot(xx,ym,2)
                 call plot(xx,ym,3)
                 call newpen(2)
                 call fillit('CI',0.05,xx,yy)
                 call newpen(1)
                 call curvit('CI',0.05,xx,yy)
        enddo
        return 
        end

        subroutine pltsort(x0,y0,xlen,ylen,xmin,xmax,tmin,tmax,wid,
     1      npt,x,y,key)
        real x0, y0, xlen, ylen, tmin, tmax, xmin, xmax,wid
        integer npt
        real x(npt), y(npt)
        integer key(npt)

        real xx, yy, dx, dy
        integer i,j

        dx =  xlen/(xmax-xmin)
        dy =  ylen/(tmax-tmin)
        call gwidth(wid)
        do 1000 i=1,npt
            j = key(i)
            xx = x0 + ( x(j) - xmin)*dx
            yy = y0 + ( y(j) - tmin)*dy
            if(i.eq.1)then
                call plot(xx,yy,3)
            else
                call plot(xx,yy,2)
            endif
 1000   continue
        call gwidth(0.0)
        return
        end

        subroutine doframe(x0,y0,xlen,ylen,tmin,tmax,xmin,xmax,
     1      vred,kvreds)
        implicit none
        real*4 x0,y0,xlen,ylen,tmin,tmax,xmin,xmax,vred
        character kvreds*(*)
        integer lgstr, ls, lx, ly
        real size
        character*80 ystr, xstr
        ls = lgstr(kvreds)
        if(vred .gt. 0.0)then
            ystr = 'T - X /'
            ystr(7+1:7+ls) =kvreds(1:ls)
            xstr = 'X (km)'
        else
            ystr = 'Time (sec)'
            xstr = 'Distance (km)'
        endif
        lx = lgstr(xstr)
        ly = lgstr(ystr)
c-----
c       compute the symbol size so that it all fits nicely
c-----
        size = 0.6* min(xlen/lx, ylen/ly, 0.05*xlen, 0.04*ylen)
            call gbox(x0,y0,x0+xlen,y0+ylen)
            call dolinx(x0,y0     ,xlen,xmin,xmax,size,
     1          .true.,.false.,.true.,
     2          lx,xstr)
            call dolinx(x0,y0+ylen,xlen,xmin,xmax,size,
     1          .false.,.true.,.false. ,
     2          0,' ')
            call doliny(x0     ,y0,ylen,tmin,tmax,size,
     1          .false.,.true. ,.true. ,ly,ystr)
            call doliny(x0+xlen,y0,ylen,tmin,tmax,size,
     1          .true. ,.true. ,.false.,0,' ')
        return
        end


      subroutine gttmp0(iprog,itot,nf1,nf2,nf34,nf5,nf67,nf10,nfilt,
     1      nup,dlam,qaqb,wref,invcsl,
     2      lstinv,twnmin,twnmax,iter,nurftn,pval,
     3      sigv,sigtp,sigts,
     3      idtwo,idum2,idum3,idum4,idum5,
     4      rdum1,rdum2,rdum3,rdum4,rdum5)
        integer NL
        parameter (NL=200)
        integer nf10(NL)
c-----
c     read control file
c-----
      open(1,file='tmpsrfi.00',form='unformatted',access='sequential')
      rewind 1
      read(1) iprog,itot,nf1,nf2,nf34,nf5,nf67,nf10,nfilt,nup,dlam
     1      ,qaqb,wref,invcsl,lstinv,twnmin,twnmax,iter,nurftn,pval,
     3      sigv,sigtp,sigts,
     3      idtwo,idum2,idum3,idum4,idum5,
     4      rdum1,rdum2,rdum3,rdum4,rdum5
      close(1,status='keep')
        return
        end

        subroutine pltmod(nid,dop)
c-----
c     plot the P or S velocity 
c     nid   I - 0 for velocity plot
c     dop   L - .true. for P
c               .false. for S
c-----
       integer nid
       logical dop
c-----
c     common for igetmod
c-----
        integer NL, NL2, NLAY
        parameter(NL=200,NLAY=200,NL2=NL+NL)

        common/isomod/d(NLAY),va(NLAY),vb(NLAY),rho(NLAY),
     1      qa(NLAY),qb(NLAY),etap(NLAY),etas(NLAY), 
     2      frefp(NLAY), frefs(NLAY)
        common/depref/refdep
        integer mmax, iunit, iiso, iflsph, idimen, icnvel
        common/modtit/title
        character title*80
c-----
c       place to save the number of layers
c-----
        integer mlyr
        common/modlly/mlyr

c-----
c     LIN - unit for FORTRAN read from terminal
c     LOT - unit for FORTRAN write to terminal
c     LER - unit for FORTRAN error output to terminal
c     NL  - number of layers in model
c-----
        real varr(4), qbarr(4), darr(4)
        real vel, vmax, vmin
c-----
c     set the plot
c-----  
        x0 = 0.5
        y0 = 2.0
        call plot(x0,y0,-3)
        xlen = 1.5
        ylen = 4.0
c-----
c     get the velocity model twice, once for the
c     current model and then for the initial model
c     however we use both to get the plot extremes
c-----
c-----
c-----
c     get bounds for the plot
c-----
        vmin = 1.0e+38
        vmax = 0.0
        qbmin = 1.0e+38
        qbmax = 0.0
        zmin = 0.0
        zmax = 0.0
        do 2000 jmd=1,2
c-----
c         get velocity model and q(beta) inverse model
c-----
            if(jmd.eq.1)then
        call getmod(2,'tmpsrfi.17',mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,.false.)
            else
        call getmod(2,'tmpmod96.000',mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,.false.)
            endif
        mlyr = mmax
            depth = 0.0
c-----
c       get the bounds for the velocity
c-----
        do 1000 i=1,mlyr
            if(nid.eq.0)then
                if(dop)then
                     vel = va(i)
                else 
                     vel = vb(i)
                endif
                if(vel.gt.vmax)vmax = vel
                if(vel.lt.vmin)vmin = vel
            else if(nid.eq.1)then
                if(qb(i) .gt.1.0)qb(i) = 1.0/qb(i)
                if(qb(i).gt.qbmax)qbmax = qb(i)
                if(qb(i).lt.qbmin)qbmin = qb(i)
            endif
            if(i.lt.mlyr)depth = depth + d(i)
 1000   continue
        if(depth+d(mlyr-1).gt.zmax)zmax = depth+d(mlyr-1)
 2000   continue
c-----
c     safety test - force a minimum separation in km/sec units
c-----
            vmax = vmax + 0.125
            vmin = vmin - 0.125
        
        darr(1) = 0.0
        darr(2) = zmax
c-----
c     now scale the velocity/Q inverse axis
c-----
        varr(1) = vmin
        varr(2) = vmax
        call dvscale(vmin,vmax)
        varr(3) = vmin
        varr(4) = (vmax - vmin)/xlen
        qbarr(1) = 0.95*qbmin
        qbarr(2) = 1.05*qbmax
        call gscale(qbarr,xlen,2,1)
        qbmin = qbarr(3)
        qbmax = qbarr(3) + xlen * qbarr(4)

        yylen = 10.0
        call gscale(darr,yylen,2,1)
        dmin = 0.0
C       dmax = darr(3) + yylen*darr(4)
        dmax = darr(2)
c-----
c     create bounding box
c-----
        call box(0.0,0.0,xlen,ylen)
c-----
c     label the depth axis
c-----
            call labz(dmin,dmax,10,0.0,ylen,ylen,0,dd)
c-----
c     plot the velocity/q inverse axis
c-----
        if(nid.eq.0)then
            if(dop)then
              call dolinx(0.0,ylen,xlen,vmax,vmin,
     1          0.05,.true.,.true.,.true.,9,'VP (KM/S)')
            else
              call dolinx(0.0,ylen,xlen,vmax,vmin,
     1          0.05,.true.,.true.,.true.,9,'VS (KM/S)')
            endif
        else
            call myaxis(0.0,ylen,'QBINV'    , 5,xlen,0.0,
     1          qbarr(3),qbarr(4))
        endif
c-----
c     plot the velocity and Qb models
c-----
c-----
c     now read the model file two times, once for original
c     and then for current
c-----
        do 3000 jmd=1,2
            if(jmd.eq.1)then
        call getmod(2,'tmpmod96.000',mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,.false.)
            call newpen(4)
            else
        call getmod(2,'tmpsrfi.17',mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,.false.)
            call newpen(2)
            endif
        mlyr = mmax
        if(nid .eq. 0)then
            dx = 1.0/varr(4)
c-----
c     plot velocity model
c-----
            if(dop)then
            call vplt(0.0,ylen,mlyr,va,d,vmin,vmax,dx,
     1          dmin,dmax,dd,jmd)
            else
            call vplt(0.0,ylen,mlyr,vb,d,vmin,vmax,dx,
     1          dmin,dmax,dd,jmd)
            endif
c-----
c     plot Qb^-1 model
c-----
        else if(nid.eq.1)then
            dx = 1.0/qbarr(4)
            call vplt(0.0,ylen,mlyr,qb,d,qbmin,qbmax,dx,
     1          dmin,dmax,dd,jmd)
        endif
 3000   continue
c-----
c     annotate
c-----
        call newpen(2)
        call plot(0.0,-0.25,3)
        call plot(0.25*xlen, -0.25, 2)
        call symbol(0.30*xlen, -0.25, 0.10, ' Current',0.0,8)
        call newpen(4)
        call plot(0.0,-0.50,3)
        call plotd(0.25*xlen,-0.50,21,0.05)
        call symbol(0.30*xlen, -0.50, 0.10, ' Initial',0.0,8)
c-----
c     reset the plot
c-----
        call plot(-x0,-y0,-3)
        return
        end

        subroutine box(xl,yl,xu,yu)
        call plot(xl,yl,3)
        call plot(xu,yl,2)
        call plot(xu,yu,2)
        call plot(xl,yu,2)
        call plot(xl,yl,2)
        return
        end
c---------------------------------------------------------------------c
c                                                                   c
c    COMPUTER PROGRAMS IN SEISMOLOGY                                c
c    VOLUME I                                                       c
c                                                                   c
c    PROGRAM: AXIS                                                  c
c                                                                   c
c    COPYRIGHT (C)  1986 R. B. Herrmann                             c
c                                                                   c
c    Department of Earth and Atmospheric Sciences                   c
c    Saint Louis University                                         c
c    221 North Grand Boulevard                                      c
c    St. Louis, Missouri 63103                                      c
c    U. S. A.                                                       c
c                                                                   c
c---------------------------------------------------------------------c
        subroutine myaxis(xpage,ypage,ititl,nchar,axlen,angle,
     1          firstv,deltav)
        real HTNUM,SYMHT,TICHT
       parameter (HTNUM=0.07,SYMHT=0.10,TICHT=0.07)
c-----
c     xpage,ypage  coordinates of starting point of axis in
c                  inches relative to current origin
c     ititl   axis title
c     nchar   number of characters in title
c             >0 for tic marks, numbers and title on
c                counterclockwise side
c             <0 for tic marks, numbers and title on
c                clockwise side
c     axlen   floating point axis length in inches
c     angle   angle (degrees) that x axis makes with
c             horizontal x-direction
c     firstv  scale value at (xpage,ypage)
c     deltav  change in slace between tic marks per unit
c             length
c             <0 value at xpage,ypage is a max, and values
c             decrease along the axis
c-----
c-----
c     we do not worry about the vagaries of storage of characters
c     since the address is passed on down to subroutine symbol
c-----
        character*(*) ititl
        a=1.0
        kn=nchar
        if(kn.lt.0)then
                a= -a
                kn= -kn
        endif
c-----
c     if deltav is too large invoke scientific notation
c-----
        ex = 0.0
        if(deltav.ne.0.0)then
                yex= alog10(abs(deltav))
                if(yex.lt.-2.0)then
                ex = aint(yex+0.01) - 1.0
            else if(yex.lt.-1.0)then
                ex = aint(yex+0.01) 
            endif
                if(yex.ge.2.0) ex = aint(yex + 0.01)
        endif
        xval = firstv*10.0**(-ex)
        xdel = deltav*10.0**(-ex)
        ct = cos(angle*0.01745329)
        st = sin(angle*0.01745329)
        ntic = axlen + 1.0
c-----
c     first put in numbers and title
c     adjust offset for numbers
c-----
        dx = - HTNUM
        dy = 1.5*a*HTNUM - 0.5*HTNUM
c-----
c     find initial position given rotation
c-----
        xn = xpage + dx*ct - dy*st
        yn = ypage + dy*ct + dx*st
        mtic = ntic/2
        do 100 i=1,ntic
            if(i.eq.1)then
                if(xval.lt.0.0)then
                    dll = 2.0*HTNUM
                else
                    dll = 1.0*HTNUM
                endif
                if(abs(xval).ge.10.0)dll = dll + HTNUM
                dxx = dll*ct
                dyy = dll*st
                    call number(xn+dxx,yn+dyy,HTNUM,xval,angle,2)
            else if(i.eq.ntic)then
                if(xval.lt.0.0)then
                    dll = 4.0*HTNUM
                else
                    dll = 3.0*HTNUM
                endif
                if(abs(xval).ge.10.0)dll = dll + HTNUM
                dxx = -dll*ct
                dyy = -dll*st
                    call number(xn+dxx,yn+dyy,HTNUM,xval,angle,2)
            else
                    call number(xn,yn,HTNUM,xval,angle,2)
            endif
                xval=xval+xdel
                xn=xn+ct
                yn=yn+st
c-----
c     halfway down axis put in title
c-----
                if(i.eq.mtic)then
                        z=kn
                        if(ex.ne.0.0)z=z+7.0
                        dx = -0.5*SYMHT*z + axlen*0.5
                        dy = (2.5*a-0.5)*SYMHT
                        xt=xpage +dx*ct-dy*st
                        yt=ypage +dy*ct+dx*st
                        call symbol(xt,yt,SYMHT,ititl,angle,kn)
                        if(ex.ne.0.0)then
                                z=kn+2
                                xt=xt+z*ct*SYMHT
                                yt=yt+z*st*SYMHT
                                call symbol(xt,yt,SYMHT,'*10',angle,3)
                                xt=xt+(3.0*ct-0.8*st)*SYMHT
                                yt=yt+(3.0*st+0.8*ct)*SYMHT
                                call number(xt,yt,0.7*SYMHT,ex,angle,
     1                                  -1)
                        endif
                endif
  100   continue
c-----
c     now put in tic marks
c-----
        call plot(xpage+axlen*ct,ypage+axlen*st,3)
        dx = - TICHT*st*a
        dy =   TICHT*ct*a
        a = ntic -1
        xn = xpage + a*ct
        yn = ypage + a*st
        do 200 i=1,ntic
                call plot(xn,yn,2)
                call plot(xn+dx,yn+dy,2)
                call plot(xn,yn,2)
                xn=xn-ct
                yn=yn-st
  200   continue
        return
        end

        subroutine labz(dmin,dmax,nd,x0,y0,ylen,id,ddd)
c-----
c     label the depth axis - however 
c-----
c     dmin    R*4 - minimum depth
c     dmax    R*4 - maximum depth
c     nd  I*4 - number of segments
c     x0  R*4
c     Y0  R*4 - upper left corner
c     ylen    R*4 - length of Y-axis
c     id  I*4 - 0 single
c               1 double
c-----
        character ostr*10, fmt*10
        dd = (dmax - dmin) / real (nd)
        if(id.eq.1)then 
            ndo = 2
        else
            ndo = 1
        endif
        if(ndo.eq.1)then
            dz = ylen / real(nd)
        else
            dz = 0.5 * ylen / real(nd)
        endif
        ddd = dz / dd
        if(dmax.gt.10 .and. dmax.lt.100)then
            ndec = 2
            fmt = '(f10.2)'
        else if(dmax.ge.100)then
            ndec = 0
            fmt = '(f10.0)'
        else
            ndec = 3
            fmt = '(f10.3)'
        endif
        do 100 n=1,ndo
            do 110 m=1,nd-1
                ypos = ylen - real(m)*dz - (n-1)*0.5*ylen
                yval = dmin + m*dd
                call plot(x0,ypos,3)
                call plot(x0+0.10,ypos,2)
                write(ostr,fmt)yval
                ht = 0.05   
                call symbol(x0-11.0*ht,ypos-0.5*ht,ht,ostr,0.0,10)
  110       continue
  100   continue
        hht = 0.07
            nchar = 10
            call symbol(x0-7.0*ht,ylen/2.0 - nchar*0.5*hht,
     1          hht,'DEPTH (KM)',90.0,nchar)
        return
        end
                
        subroutine gcmdln(dop,vred,kvreds,wid)
        implicit none
c-----
c     parse the command line parameters
c-----
c     dop       L - .true. plot P times
c                   .false. plot S times
c     vred      R - reduction velocity for plot
c                   if negative do not do reduced travel time
c     kvreds    C*20 the vred string for graph annotation
c     wid       R - width of the predicted travel time curve
c                   default is 0.0
c-----
        logical dop
        real vred, wid
        character kvreds*20
        integer mnmarg, nmarg, i
        character names*20
        nmarg = mnmarg()
c-----
c     initialize defaults
c-----
        dop = .false.
        vred = -1.0
        kvreds = ' '
        wid = 0.0
c-----
c     loop through command line arguments
c-----

        i = 0
 1000   continue
            i = i + 1
            if(i.gt.nmarg)go to 2000
            call mgtarg(i,names)
            if(names(1:2).eq.'-P' .or. names(1:2).eq.'-p')then
                dop = .true.
            else if (names(1:2).eq.'-S' .or. names(1:2).eq.'-s')then
                dop = .false.
            else if (names(1:2).eq.'-V' .or. names(1:2).eq.'-v')then
                i = i + 1
                call mgtarg(i,names)
                kvreds = names
                read(names,'(bn,f10.0)')vred
            else if (names(1:2).eq.'-W' .or. names(1:2).eq.'-w')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f10.0)')wid
            else if (names(1:2) .eq. '-h')then
                call usage()
            else if (names(1:2) .eq. '-?')then
                call usage()
            endif
        go to 1000
 2000   continue
c-----
c       safety checks
c-----
        if(wid.lt.0.0)wid = 0.0
        return
        end

        subroutine usage()
        implicit none
c------
c       write out program syntax
c-----
        integer LER,LOT,LIN
        parameter(LER=0, LIN=5, LOT=6)
        write(LER,*)'USAGE: ',
     1  'ttinvpv96 [-P] [-S] [-W width] [-VRED vred] ',
     1      ' [-h] [-?] '
        write(LER,*)
     1  '-P           (default true )  plot P time'
        write(LER,*)
     1  '-S           (default false)  plot S time'
        write(LER,*)
     1  '-W   width (default 0.00) Width of line (inch) for model plot'
        write(LER,*)
     1  '                    0.05 is large enough '
        write(LER,*)
     1  '-VRED vred (default not used) reduction velocity'
        write(LER,*)
     1  '-?                   Display this usage message'
        write(LER,*)
     1  '-h                   Display this usage message'
        stop
        end


                
        subroutine vplt(xpos,ypos,n, x, y,xmin,xmax,dx,
     1          ymin,ymax,dy,iskp)
c-----
c     y   R*4 - array of parameters to be plotted along y-axis
c     x   R*4 - array of parameters to be plotted along x-axis
c     iskp    I*4 - 2 solid line
c               1 dashed line
c-----
c     plotting velocity model is odd since have to force the steps
c-----
        parameter(NL=200,NL2=NL+NL)
        real*4 x(*), y(*)
        dphl = 0.0
        dphu = y(1)
        do 100 i=2,n
            xx1 = xpos + (x(i-1) - xmin)*dx
            xx2 = xpos + (x(i  ) - xmin)*dx
            yy1 = ypos - (dphl - ymin)*dy
            yy2 = ypos - (dphu - ymin)*dy
            call plot(xx1,yy1,3)
            if(iskp.eq.2)then
                call plot(xx1,yy2,2)
                call plot(xx2,yy2,2)
            else
                call plotd(xx1,yy2,21,0.05)
                call plotd(xx2,yy2,21,0.05)
            endif
            dphl = dphu
            dphu = dphl + y(i)
  100   continue
        dphu = 1.10*dphl
        yy2 = ypos - (dphu - ymin)*dy
        if(iskp.eq.2)then
            call plot(xx2,yy2,2)
        else
            call plotd(xx2,yy2,21,0.05)
        endif
        return
        end

        subroutine dvscale(vmin,vmax)
c-----
c     rescale xval so that a plot will look nice
c     this is currently for velocity which can never be negative
c-----
c     vmin    R   - value to be scaled
c     vmax    R   - value to be scaled
c-----
c-----
c     first get size
c-----
        size = max(vmin,vmax)
        if(size.gt.1)then
            lpow = alog10(size)
            tmp = 10.0**lpow
            if(size.lt.tmp)then
                tmp = 0.1 * tmp
            endif
        else
            lpow = alog10(size) + 1
            tmp = 10.0**lpow
            if(size.lt.tmp)then
                tmp = 0.1 * tmp
            endif
        endif
        vmn = vmin
        do 1000 i=0,20
            xx = real(i) * 0.5 *tmp
            if(xx .lt. vmin)then
                vmn = xx 
            endif
 1000   continue
        vmx = vmax
        do 2000 i=20,0,-1
            xx = real(i) * 0.5 *tmp
            if(xx .gt. vmax)then
                vmx = xx 
            endif
 2000   continue
C       write(6,*)vmin,vmax,vmn,vmx,size,lpow,tmp
        if(vmx.gt.vmn)then
            vmin = vmn
            vmax = vmx
        endif
C       write(6,*)vmin,vmax,vmn,vmx,size,lpow,tmp

        return
        end
