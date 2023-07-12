        program refmod96
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME V                                                       c
c                                                                     c
c      PROGRAM: REFMOD96                                              c
c                                                                     c
c      COPYRIGHT 2005 R. B. Herrmann                                  c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
c       This program plot model96 models files
c
c       refmod96 -XLEN xlen -YLEN ylen -X0 x0 -Y0 y0 
c           -VMIN vmin -VMAX vmax -ZMIN zmin -ZMAX zmax
c           -KF k1 -KR k2 -W width
c           model96_file(s)
c-----
c       CHANGES
c       12 AUG 2005 - created building on timmod96
c       04 AUG 2006 - corrected error in first arrival pick that
c               falsely gave the refraction time instead of the
c               direct time because the refraction arrival was
c               unphysical
c       22 JAN 2008 - rearranged routines, put in new common
c                     isotropic travel time routines from time96
c       25 JAN 2008 - put Radius of Earth into common/earth/radius for
c                      generality
c                   -  define a separate common block for the
c                      SH velocity and density
c                   -  have sphericity correction
c                      work on common blocks instead of procedure call
c                   -  create a default adoflt  to fill the SH for a flat model
c                      note the separation of SH is important for wavenumber
c                      integration code
c       25 JAN 2008 -  permit xmin xmax to be negative so that a split spread can
c                      be plotted
c       08 FEB 2008 -  subtle change in fstarr for source receiver in same layer
c                      spherical mapping was not done
c       04 APR 2008 -  permitted source depth - note however that the 
c                      computation for multiples assumes that the source and
c                      receiver are both at the surface
c       18 FEB 2009 -  took car to add common/depref/refdep in frstar
c       01 JUN 2013 -  Modified subroutine to prevent NaN for direct ray by replacing
c                      pnew = 0.999d+00 * pupper to
c                      pnew = 0.99d+00 * pupper
c                      also modified subroutine frstar to have the dogeom argument.
c                      We do not want to compute teleseismic geometrical spreading
c       21 JAN 2014    Corrected plotting and TXT output for reflection
c                      multiples. In subroutine dorefl
c                      -TXT flag creates REFMOD96.TXT tabulation
c       01 JUN 2016    Added clipping in the doplttim routine
c       14 AUG 2017    correctly implemented -HS source depth for reflection multiples
c                      Thoughts on multiples reflections in layer 
c
c  
c-----
        implicit none
c-----
c       LIN I*4 - logical unit for standard input
c       LOT I*4 - logical unit for standard output
c-----
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)
c-----
c       command line arguments
c-----
        real x0, y0, xlen, ylen, tmin, tmax, xmin, xmax
        real vred, hs, hr, wid
        integer k1,k2
        character kvreds*80
        character*80  mname
        logical doleg, dobox
        integer kpar
        logical dop, dosv, dosh, dotxt
        real dx, dist
        integer i,m,ipen, nmult
        real pvel, svel, vsa, vsb, vsr, rayp, geom, tstar, den
        real xleg, yleg, dy, ht
        integer lgstr, ls

        integer NDIST
        parameter (NDIST=201)
        real*4 t(NDIST), r(NDIST)
        logical doit

        logical ext
        
        common/earth/radius
        real radius
        real kmdeg
c-----
c       initialize
c-----
        radius = 6371.
        kmdeg=radius*3.1415927/180.0
c-----
c       parse the command line
c-----
        call gcmdln(x0,y0,xlen,ylen,tmin,tmax,xmin,xmax,
     1      vred,mname,k1,k2,hs,hr,wid,dop,dosv,dosh,
     2      doleg,dobox,kvreds,nmult,dotxt)
c-----
c       open REFMOD96.TXT if -TXT in command laine
c-----
        if(dotxt)then
           open(4,file='REFMOD96.TXT',access='sequential',
     1       form='formatted',status='unknown')
           rewind 4
        endif
c-----
c       open plot file
c-----
        call pinitf('REFMOD96.PLT')
        call newpen(1)
        if(dobox)then
            call doframe(x0,y0,xlen,ylen,tmin,tmax,xmin,xmax,
     1          vred,kvreds)
            else
c-----
c               alignment tics
c-----
                call plot(x0,y0+0.1,3)
                call plot(x0,y0,2)
                call plot(x0+0.1,y0,2)

                call plot(x0+xlen,y0+0.1,3)
                call plot(x0+xlen,y0,2)
                call plot(x0+xlen-0.1,y0,2)

                call plot(x0,y0+ylen-0.1,3)
                call plot(x0,y0+ylen,2)
                call plot(x0+0.1,y0+ylen,2)

                call plot(x0+xlen,y0+ylen-0.1,3)
                call plot(x0+xlen,y0+ylen,2)
                call plot(x0+xlen-0.1,y0+ylen,2)


        endif
c-----
c       compute the first arrival times for the models
c-----
        dx = (xmax - xmin)/(NDIST -1)
        yleg = y0 + ylen
        xleg = x0 + xlen
c-----
c       check to see that the model exists
c-----
        inquire(file=mname,exist=ext)
        if(.not.ext)then
            call usage('Error: model does not exist')
        endif
c-----
c       beginning of big loop for P SV and SH
c-----
        do 1234 kpar = 1,3
            doit = .false.
            if(kpar.eq.1 .and. dop )doit = .true.
            if(kpar.eq.2 .and. dosv)doit = .true.
            if(kpar.eq.3 .and. dosh)doit = .true.
            if(doit)then
            do 100 i=1,NDIST
                dist = xmin + ( i-1)*dx
                r(i) = dist
                if(kpar.eq.1)then
                 call frstar(abs(dist),hs,hr,mname,1,t(i),pvel,svel,
     1              den,vsa, vsb, vsr, rayp, geom,tstar,.false.,.false.)

                    if(vred.gt.0)then
                        t(i)  = t(i)  - abs(dist)/vred
                    endif
                else if(kpar.eq.2)then
                 call frstar(abs(dist),hs,hr,mname,2,t(i),pvel,svel,
     1              den,vsa, vsb, vsr, rayp, geom,tstar,.false.,.false.)
                    if(vred.gt.0)then
                        t(i) = t(i) - abs(dist)/vred
                    endif
                else if(kpar.eq.3)then
                 call frstar(abs(dist),hs,hr,mname,3,t(i),pvel,svel,
     1              den,vsa, vsb, vsr, rayp, geom,tstar,.false.,.false.)
                    if(vred.gt.0)then
                        t(i) = t(i) - abs(dist)/vred
                    endif
                endif
                if(dotxt)then
                   if(kpar.eq.1)then
                      write(4,*)'P  ','FirstArrival ',r(i),t(i)
                   else if(kpar.eq.2)then
                      write(4,*)'SV ','FirstArrival ',r(i),t(i)
                   else if(kpar.eq.3)then
                      write(4,*)'SH ','FirstArrival ',r(i),t(i)
                   endif
                endif
 100            continue
c-----
c               PLOT FIRST ARRIVAL
c-----
                if(k1.gt.1)then
                    call newpen(k1)
                else
                    call newpen(2)
                endif
                if(kpar.eq.1)then
                call doplttim(x0,y0,xlen,ylen,xmin,xmax,tmin,tmax,
     1              r,t,NDIST,wid)
                else if(kpar.eq.2)then
                call doplttim(x0,y0,xlen,ylen,xmin,xmax,tmin,tmax,
     1              r,t,NDIST,wid)
                else if(kpar.eq.3)then
                call doplttim(x0,y0,xlen,ylen,xmin,xmax,tmin,tmax,
     1              r,t,NDIST,wid)
                endif
C               if(doleg)then
C                   call newpen(ipen)
C                   call plot(xleg+0.1,yleg,3)
C                   call gwidth(wid)
C                   call plot(xleg+0.3,yleg,2)
C                   call newpen(1)
C                   dy = 0.14
C                   ht = 0.71 * dy
C                   ls = lgstr(mname)
C                   call gwidth(0.001)
C                   call symbol(xleg+0.4,yleg-0.5*ht,ht,
C     1                 mname,0.0,ls)
C                   yleg = yleg - dy
C                   call gwidth(0.001)
C               endif
                call newpen(1)
c-----
c               do reflections
c               reflections will be plot using dashed lines
c               eventually add multiples
c-----      
                if(k2.eq.0)then
                    call newpen(4)
                else
                    call newpen(k2)
                endif
                call dorefl(nmult,kpar,x0,y0,xlen,ylen,vred,
     1              tmin,tmax,xmin,xmax,hs,wid,dotxt)
c-----
c               close open files
c-----
                call newpen(1)
            endif
 1234   continue
        if(dotxt)then
            close(4)
        endif

            call pend()
C           write(LER,*)'vred   :',vred
C           write(LER,*)'xmin   :',xmin
C           write(LER,*)'xmax   :',xmax
C           write(LER,*)'tmin   :',tmin
C           write(LER,*)'tmax   :',tmax
        end

        subroutine dorefl(nmult,kpar,x0,y0,xlen,ylen,vred,
     1          tmin,tmax,xmin,xmax,hs,wid,dotxt)
c-----
c       nmult   I   - number of multiple reflections
c       kpar    I   - 1 P, 2 SV, 3 SH
c               note:   no coversions
c                   surface source, receiver not VSP
c                   clipping invoked
c                   check for velocity contrast
c       vred    R   - reduction velocity
c       hs      R   - source depth
c       wid     R   - pen width for plot
c       dotxt   L   - .true. put times in REFMOD96.TXT
c-----
        implicit none 

        real x0, y0, xlen, ylen
        real vred
        real tmin, tmax, xmin,xmax, hs, wid
        integer nmult, kpar, iret
        logical dotxt

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
        
        integer NRAY
        parameter (NRAY=201)

        real pmin, pmax, dp, p
        real rcoef
        integer i, jmult, jupdn, ipen, m

        real time, dist
        real xx, yy, dx, dy
c-----
c       define scaling
c-----
        dx = xlen/(xmax - xmin)
        dy = ylen/(tmax - tmin)

        
c-----
c       define the clipping region
c-----
        call gclip('on' ,x0,y0,x0+xlen,y0+ylen)
c-----
c       determine upper bound of ray parameter which will be
c       taken to be 0.9999 of the surface ray parameter
c-----
        pmin = 0.0
        if(kpar.eq.1)then
            pmax = 0.9999/a(1)
        else if(kpar.eq.2)then
            pmax = 0.9999/b(1)
        else if(kpar.eq.3)then
            pmax = 0.9999/b(1)
        endif
        if(xmin.lt.0)then
            pmin = -pmax
        endif
        dp = (pmax - pmin)/(NRAY-1)
        do m=1,mmax-1
            do jmult=0,nmult-1
               do jupdn = 0,1
c----
c this only works for a surface source and receiver change this
c-----
                ipen = 3
c-----
c     call this for vertical incidence to get reflection coefficient
c-----
                call ptau(0.0,jmult,
     1                       jupdn,kpar,m,iret,time,dist,hs,rcoef)
                call gwidth(wid*abs(rcoef))
                do i=1,NRAY
                    p = pmin + ( i-1)*dp
                    call ptau(p,jmult,
     1                       jupdn,kpar,m,iret,time,dist,hs,rcoef)

                    if(iret.gt.0)then
                        time = time +p*dist 
                        if(dotxt)then
            write(4,*)'Layer,multiple,(r,t):',m,jmult,jupdn,dist,time
                        endif
                        if(vred.gt.0.0)then
                            time = time - dist/vred
                        endif
                        xx = x0 + (dist-xmin)*dx
                        yy = y0 + (time-tmin)*dy
                        call plot(xx,yy,ipen)
                        ipen = 2
                    endif
                enddo
                call gwidth(0.0)
               enddo
            enddo
        enddo

        call gclip('off',x0,y0,x0+xlen,y0+ylen)
        return
        end

        subroutine ptau(p,jmult,jupdn,kpar,mm,iret,time,dist,hs,rcoef)
c-----
c       p       R   - ray parameter
c       jmult   I   - multiple in layer
c                     for source layers
c       jupdn   I   - 1 go down from hs, 0 go up from hs
c (jmult,jupdn)= (0,0) (0,1)  (1,0)   (1,1)      (2,0)      (2,1)
c                 /       /   /\  /    /\  /  /\  /\  /    /\  /\  /       
c                       \/      \/   \/  \/     \/  \/   \/  \/  \/ 
c       kpar    I   - 1 P, 2 SV, 3 SH
c       mm  I   - index of deepest layer but not halfspace
c       iret    I   - .gt.0 success
c       time    R
c       dist    R
c       hs      R - source depth
c       rcoef   R - acoustic reflection coefficient 
c-----
c       note for the deepest layer the velocities between mm mm+1 must 
c       differ
c-----
        implicit none
        real p, time, dist, hs, rcoef
        integer kpar, mm, iret

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

        real v, v1, v2, pv
        real r1, r2
        double precision sr,srf
        integer i, jmult, jupdn

        real zboundary
        
        iret = 0
c-----
c       test for a reflection
c-----
        if(kpar.eq.1)then
            v1 = a(mm)
            v2 = a(mm+1)
            r1 = rho(mm)
            r2 = rho(mm+1)
        else if(kpar.eq.2)then
            v1 = b(mm)
            v2 = b(mm+1)
            r1 = rho(mm)
            r2 = rho(mm+1)
        else if(kpar.eq.3)then
            v1 = b(mm)
            v2 = b(mm+1)
            r1 = rho(mm)
            r2 = rho(mm+1)
        endif
        rcoef = (r2*v2 - r1*v1)/(r2*v2 + r1*v1)
        if(v1.eq.v2)then
            return
        endif
c------
c       OK let s compute P-tau
c-----
        time = 0.0
        dist = 0.0
        zboundary = 0.0
        do i=1,mm
            
            if(kpar.eq.1)then
                v = a(i)
            else if(kpar.eq.2)then
                v = b(i)
            else if(kpar.eq.3)then
                v = b(i)
            endif
            pv = p*v
            srf = 1.0d+00 - dble(pv)*dble(pv)
            if(srf.le.0.0d+00)return
            sr = dsqrt(srf)
            if(jupdn.eq.0)then
c (jmult,jupdn)= (0,0) (0,1)  (1,0)   (1,1)      (2,0)      (2,1)
c                 /       /   /\  /    /\  /  /\  /\  /    /\  /\  /       
c                       \/      \/   \/  \/     \/  \/   \/  \/  \/ 
                 if(hs .gt. zboundary)then
                      dist = dist + (1 +2*(jmult))*d(i)*pv/sr
                      time = time + (1 +2*(jmult))*d(i)*sr/v
                 else
                      dist = dist + (   2*(jmult))*d(i)*pv/sr
                      time = time + (   2*(jmult))*d(i)*sr/v
                 endif
            else
                 if(hs .gt. zboundary)then
                      dist = dist + (1 +2*(jmult))*d(i)*pv/sr
                      time = time + (1 +2*(jmult))*d(i)*sr/v
                 else
                      dist = dist + (2 +2*(jmult))*d(i)*pv/sr
                      time = time + (2 +2*(jmult))*d(i)*sr/v
                 endif
            endif
            zboundary = zboundary + d(i)
        enddo
        iret = 1
            
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
        size = 0.8* min(xlen/lx, ylen/ly, 0.05*xlen, 0.04*ylen)
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


        subroutine doplttim(x0,y0,xlen,ylen,xmin,xmax,tmin,tmax,
     1      r,tp,NDIST,wid)
        implicit none
        real*4 x0, y0, xlen, ylen, tmin, tmax, xmin, xmax
        integer NDIST
        real r(ndist), tp(ndist)
        real wid
        real dx, dy, xx, yy
        integer i
c-----
        dx =  xlen/(xmax-xmin)
        dy =  ylen/(tmax-tmin)
        
        call gclip('on',x0,y0,x0+xlen,y0+ylen)
        call gwidth(wid)
        do 100 i=1,NDIST
            xx = x0 + ( r (i) - xmin)*dx
            yy = y0 + ( tp(i) - tmin)*dy
            if(i.eq.1)then
                call plot(xx,yy,3)
            else
                call plot(xx,yy,2)
            endif
 100    continue
        call gwidth(0.0)
        call gclip('off',x0,y0,x0+xlen,y0+ylen)
        return
        end
        

        subroutine gcmdln(x0,y0,xlen,ylen,tmin,tmax,xmin,xmax,vred,
     1      mname,k1,k2,hs,hr,wid, dop,dosv,dosh,
     2      doleg,dobox,kvreds,nmult,dotxt)
c-----
c       parse the command line arguments
c-----
c       x0  R*4 - lower left corner of plot frame
c       y0  R*4 - lower left corner of plot frame
c       xlen    R*4 - width  of plot frame on page
c       ylen    R*4 - height of plot frame on page
c       k1  I   - pen color for first arrival
c       k2  I   - pen color for reflection
c       tmin    R*4 - minimum time for plot
c       tmax    R*4 - maximum time for plot
c       xmin    R*4 - minimum distance for plot
c       xmax    R*4 - mazimum distance for plot
c       vred    R*4 - reduction velocity (-10000 == do not apply)
c       mname   Ch*80   - array of model files
c       hs  R*4 - source depth
c       hr  R*4 - receiver depth
c       wid R   - line width for travel time curve
c       dop L   - .true. plot P arrival
c       dosv    L   - .true. plot SV arrival
c       dosh    L   - .true. plot SH arrival
c       kvreds  Ch*80   - Vred string for plot annotation
c       doleg   L   - .true. put in file legend to right of plot
c       dobox   L   - .true. put in frame
c       nmult   I   - number of multiples - default = 1
c       dotxt   L   - .true. create REFMOD96.TXT file
c-----
        implicit none
        real*4 x0, y0, xlen, ylen, xmin, xmax, tmin, tmax, vred
        real hs, hr, wid
        integer k1,k2, nmult
        character*80 mname
        character*80 kvreds
        logical doleg, dobox
        logical dop, dosv, dosh, dotxt

        integer mnmarg, i, nmarg
        character name*80
        real tmp

        hs = 0.0
        hr = 0.0
        x0 = 2.0
        y0 = 1.0
        xlen = 6.0
        ylen = 6.0
        k1 = 0
        k2 = 0
        tmin = 0.0
        tmax = 100.0
        xmin = 0.0
        xmax = 300.0
        vred = -10000.0
        wid = 0.001
        dop = .false.
        dosv = .false.
        dosh = .false.
        dotxt = .false.
        kvreds = ' '
        doleg = .false.
        dobox = .true.
        mname = ' '
        nmult = 1
        
        nmarg = mnmarg()
        i = 0
 1000   continue
            i = i + 1
            if(i.gt.nmarg)go to 2000
            call mgtarg(i,name)
            if(name(1:3).eq.'-KF')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,i10)')k1
            else if(name(1:3).eq.'-KR')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,i10)')k2
            else if(name(1:5).eq.'-TMIN')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')tmin
            else if(name(1:5).eq.'-TMAX')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')tmax
            else if(name(1:5).eq.'-XMIN')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')xmin
            else if(name(1:5).eq.'-XMAX')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')xmax
            else if(name(1:5).eq.'-VRED')then
                i = i + 1
                call mgtarg(i,name)
                kvreds = name
                read(name,'(bn,f10.0)')vred
            else if(name(1:3).eq.'-X0')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(f10.0)')x0
            else if(name(1:3).eq.'-Y0')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(f10.0)')y0
            else if(name(1:5).eq.'-XLEN')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')xlen
            else if(name(1:5).eq.'-YLEN')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')ylen
            else if(name(1:2).eq.'-M')then
                i = i + 1
                call mgtarg(i,mname)
            else if(name(1:2).eq.'-W')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')wid
                if(wid .le. 0.0)wid = 0.001
            else if(name(1:2) .eq. '-P')then
                dop  = .true.
            else if(name(1:3) .eq. '-SV')then
                dosv = .true.
            else if(name(1:3) .eq. '-SH')then
                dosh = .true.
            else if(name(1:4).eq.'-LEG')then
                doleg = .true.
            else if(name(1:3).eq.'-NO')then
                dobox = .false.
            else if(name(1:3).eq.'-HS')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(f10.0)')hs
            else if(name(1:3).eq.'-HR')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(f10.0)')hr
            else if(name(1:3).eq.'-NM')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,i10)')nmult
            else if(name(1:4).eq.'-TXT')then
                dotxt = .true.
            else if(name(1:2) .eq. '-?')then
                call usage(' ')
            else if(name(1:2) .eq. '-h')then
                call usage(' ')
            endif
        go to 1000
 2000   continue
c-----
c       error checking
c-----
        if(mname .eq. ' ')then
            call usage('Error: no model given')
        endif
        if(.not.dop .and. .not. dosv .and. .not. dosh)then 
            call usage('Error: no phase given,eg., -P,-SH')
        endif
        if(tmax.lt.tmin)then
              tmp = tmax
              tmax = tmin
              tmin = tmp
        endif
        if(xmax.lt.xmin)then
              tmp = xmax
              xmax = xmin
              xmin = tmp
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
        write(LER,*)'Usage: refmod96 -XLEN xlen -YLEN ylen',
     1      ' -X0 x0 -Y0 y0 [-P -SV -SH ]',
     2      '-TMIN tmin -TMAX tmax -XMIN xmin -XMAX xmax -W width',
     3      '-KR kr -KF kf -VRED vred -M  model96_file -NOBOX ',
     4      '-NMULT nmult'
        write(LER,*)'Plot P-wave first arrival times'
        write(LER,*)
     1  '-XLEN xlen (default 6.0  ) Length X-axis'
        write(LER,*)
     1  '-YLEN ylen (default 6.0  ) Length Y-axis'
        write(LER,*)
     1  '-VRED vred (default not used) reduction velocity'
        write(LER,*)
     1  '-TMIN tmin (default 0.0  )  Minimum value of time'
        write(LER,*)
     1  '-TMAX tmax (default 100.0)  Maximum value of time'
        write(LER,*)
     1  '-XMIN xmin (default  0.0 )  Minimum value of distance'
        write(LER,*)
     1  '-XMAX xmax (default 300.0)  Maximum value of distance'
        write(LER,*)
     1  '-X0 x0     (default  2.0 )  x-position of lower left corner'
        write(LER,*)
     1  '-Y0 y0     (default  1.0 )  y-position of lower left corner'
        write(LER,*)
     1  '-W   width (default 0.001) Width of line (inch) for model plot'
        write(LER,*)
     1  '-KR kr   (default  4   ) Color of reflection'
        write(LER,*)
     1  '-KF kf   (default  2   ) Color of first arrival'
        write(LER,*)
     1  '-M model   (required    ) model96 file name'
        write(LER,*)
     1  '-P        (default      )  Plot P times'
        write(LER,*)
     1  '-SV       (default P    )  Plot SV times'
        write(LER,*)
     1  '-SH       (default P    )  Plot SH times'
        write(LER,*)
     1  '-NOBOX    (default none) do not plot bounding frame'
        write(LER,*)
     1  '-NMULT nmult (default 1) number of multiples'
        write(LER,*)
     1  '-LEG       (default none) Put in file legend'
        write(LER,*)
     1  '-HS       (default 0.0  )  source depth in km - in model'
        write(LER,*)
     1  '-?        (default none )  this help message '
        write(LER,*)
     1  '-h        (default none )  this help message '
        stop
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
                bsh(m+1) = b(m)
                qbsh(m+1) = qb(m)
                rhosh(m+1) = rho(m)
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
        real depths, depthr
        real dphs, dphr, dphref
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
     1            hs+refdep, hr+refdep, ipsvsh,iflsph, rp(1),
     2            ts, dolock,
     3            mname, varec, vbrec, rhorec,
     4            vasrc, vbsrc, rhosrc)

              call fstarr(r+500,tm,lmaxs, lmaxr, lmaxref,
     1            hs+refdep, hr+refdep, ipsvsh,iflsph, rp(2),
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
c       geom R   geometrical spreading factor
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
        real rayp
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
            do 1111 iter=1,10
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
                    if(iter.lt.10)then
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
            b(i)=b(i)*tmp
            bsh(i)=b(i)
            qbsh(i)=qb(i)
            rhosph=rho(i)
            rhosh(i) = rhosph * tmp **(-5.0)
            rho(i) = rhosph * tmp **(-2.275)
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

