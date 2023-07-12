        program fmplot
c-----
c       Program to plot focal mechanism
c       Original version in 
C       Pujol, J., and R. B. Herrmann (1990).
C       A student s guide to point sources in homogeneous media,
C       Seism. Res. Letters 61, 209-224.
C-----
C       CHANGES
C       31 AUG 2003 - 
c           in the effort to remove any Numerical Recipe routines
C           The eigenvectors for plotting P and T axes was not
c           properly used
c       09 JAN 2004 - added -TP -tP options  
c       19 AUG 2004 - error in tensor, replaced xmax=a0 to xmax=-1.0e+37
c       27 JAN 2005 - changed default - if strike,dip,rake,force or
c           moment tensor is NOT given, then only first motions are 
c           plotted
c     13 FEB 2008 - In an effort to both comply with and promulgate
c          international standards, the NEIC Policy and Procedures
c          Working Group has made the decision to adopt the
c          international standard for converting scalar moment
c          to moment magnitude proposed by the IASPEI Commission on
c          Seismological Observation and Interpretation (CoSOI)
c          Working Group on Magnitude --
c
c         old     Mw = 2/3 log Mo - 10.7    Mo in dyne-cm
c         new     Mw = 2/3(log Mo - 16.1)   Mo in dyne-cm
c                 Mw = 2/3(log Mo -  9.1)   Mo in N-m
c     19 NOV 2011 - added a -V flag to see the program controls,
c         rewrote the gcmdln parsing of the commands starting
c         with '-f' or '-F'
c     07 FEB 2014 - corrected the test for inconsistent amplitude
c         for P-wave first motion. The "Inconsistent" was not 
c         properly defined for an pper hemispehre projection.
c         This only affected the output and
c         did not affect the plots
c     22 AUG 2018 - added -NOSYMMETRIC to permit a non-symmetric 
c         moment tensor
c     28 AUG 2020 - add an argument to -K. before this turned on colors for
c         positive and negative contours. Now the -KC  does this if
c         and -K kval uses kval for the fill or nodal planes
c     04 JAN 2021 - change all do XXX yyy ... XXX a=b 
c                  to   do yy ... a=b .. enddo
c           to be compatible with gfortran 9.3
c     15 SEP 2021 - also output LUNE.PLT for discrimination
c           thansk to C. J. ammon of Penn State
c         
c-----

        integer LOT
        parameter (LOT=6)
        integer NTR
        parameter (NTR=90 )
        real*4 val(NTR,2), tval(NTR)
        real*4 mom(3,3), gamma(3), phi(3), theta(3), f(3)
        common/pltcon/x0,y0,rad,eqarea
        logical eqarea
        common/pltmap/tr0,dtr,pl0,dpl

        parameter (mp=3,np=3)
        real*8 xmt(np,np),ev(np),ev1(np),z(np,np)

        character dfile*80
        character title*80
        character subtit*80
        character project*10

        integer kval
        logical ifm
        logical fmamp
        integer fmfill
        logical fmtype
        logical verbos
        logical nullit
        logical hastit
        logical hassub
        logical pltmch
        logical kolor
        logical hemis
        logical doptax, doptlabel
        logical domom
        logical dofmtp
        logical donodal
        logical verby
        logical donode
        
        real*4 titsiz

        character sym(5)*3
        data sym/'  P',' SV',' SH','POL',' S '/
        data mom/9*0.0/, f/3*0.0/


        call gcmdln(eqarea,mom,itype,rad,x0,y0,dfile,fmfill,
     1      fmamp,fmtype,verbos,nullit,title,hastit,subtit,
     2      pltmch,hassub,titsiz,kolor,hemis,doptax,domom,f,dofmtp,
     3      donodal,doptlabel,verby,donode,kval)
       if(verby)then
            WRITE(6,*)'eqarea  :',eqarea
            WRITE(6,*)'mom     :',mom
            WRITE(6,*)'rad     :',rad
            WRITE(6,*)'itype   :',itype
            WRITE(6,*)'x0      :',x0
            WRITE(6,*)'y0      :',y0
            WRITE(6,*)'dfile   :',dfile
            WRITE(6,*)'fmfil   :',fmfill
            WRITE(6,*)'fmamp   :',fmamp
            WRITE(6,*)'fmtype  :',fmtype
            WRITE(6,*)'verbos  :',verbos
            WRITE(6,*)'nullit  :',nullit
            WRITE(6,*)'title   :',title
            WRITE(6,*)'hastit  :',hastit
            WRITE(6,*)'subtit  :',subtit
            WRITE(6,*)'pltmch  :',pltmch
            WRITE(6,*)'hassub  :',hassub
            WRITE(6,*)'titsiz  :',titsiz
            WRITE(6,*)'kolor   :',kolor 
            WRITE(6,*)'kval    :',kval  
            WRITE(6,*)'hemis   :',hemis 
            WRITE(6,*)'doptax  :',doptax 
            WRITE(6,*)'domom   :',domom 
            WRITE(6,*)'dofmtp  :',dofmtp 
            WRITE(6,*)'donodal :',donodal 
            WRITE(6,*)'f       :',f 
            WRITE(6,*)'dofmtp  :',dofmtp 
            WRITE(6,*)'donode  :',donode 
       endif

c-----
c       check for first motion data
c-----
        if(dfile(1:1).ne.' ')then
            ifm = .true.
            open(2,file=dfile,status='unknown',form='formatted',
     1          access='sequential')
            rewind 2
        else
            ifm = .false.
        endif
c-----
c       check for title length
c-----
        if(hastit)then
            ltit = lgstr(title)
        endif
        if(hassub)then
            lsub = lgstr(subtit)
        endif
        if(domom)then
            call gtev(xmt,mom,np,ev,ev1,ilg,ism,z)
        endif

c-----
c       initialize graphics
c-----

        call pinitf('FMPLOT.PLT')
c-----
c       sample the radiation pattern on focal sphere, 
c       and save value of amplitude on the focal sphere
c-----
        wid = 0.02*rad
        if(wid .lt. 0.005)wid = 0.005
        if(fmfill .gt. 0)then
            npl = 3*(rad/wid )
        else
            npl = 20
        endif
        dtr = 360.0 / (NTR - 1)
        dpl = 90.0 / (npl - 1)
        tr0 = 0.0
        pl0 = 0.0
        open(1,status='scratch',form='unformatted',access='sequential')
        rewind 1
        vmax = -1.0e+38
        vmin =  1.0e+38
        if(pltmch)then
        do 1000 ip = 1, npl
            if(fmfill .gt. 0)then
                r = ip*wid/3.0
                if(eqarea)then
                    pl = 2.0*asin(0.707*r/rad)
                else
                    pl = 2.0*atan(r/rad)
                endif
            else
                pl = pl0 + (ip-1)*dpl
                pl = pl * 3.1415927/180.0
            endif
            do 1100 it = 1, NTR
                tr = tr0 + (it-1)*dtr
c-----
c           convert to spherical coordinates on lower unit hemisphere
c-----
                tr = tr * 3.1415927/180.0
                if(hemis)then
                call getmtn(tr,pl,gamma,phi,theta,mom,sump,sumsv,sumsh,
     1              domom,f)
                else
                call getmtn(tr,-pl,gamma,phi,theta,mom,sump,sumsv,sumsh,
     1              domom,f)
                endif
c-----
c       get  P-wave pattern for itype = 1 
c       get SV-wave pattern for itype = 2 
c       get SH-wave pattern for itype = 3 
c       get S  polarization for itype = 4
c       get S  wave pattern for itype = 5
c-----
            if(itype.eq.1)then
                sum = sump
            else if(itype.eq.2)then
                sum = sumsv
            else if(itype.eq.3)then
                sum = sumsh
            else if(itype.eq.4)then
                sum = sump
                if(mod(it,5).eq.1 .and. mod(ip,3).eq.0)then
                if(sumsh.ne.0.0 .or. sumsv.ne.0.0)then
                    if(hemis)then
                        pol = atan2(sumsh,sumsv)
                    else
                        pol = -atan2(sumsh,sumsv)
                    endif
                        call pltpol(ip,it,pol,tr)
                endif
                endif
            else if(itype.eq.5)then
                sum = sqrt(sumsh**2 + sumsv**2)
            endif

            tval (it) = sum
            if(sum.gt.vmax)vmax = sum
            if(sum.lt.vmin)vmin = sum
 1100       continue
            write(1)(tval(it),it=1,NTR)
 1000   continue
            call circle(x0,y0,rad,.false.)
        endif
c-----
c       define symbol height
c-----
        ht = 0.14*rad
        hht = titsiz
c-----
c       put on the North symbol if verbos
c-----
        if(verbos)then
            call gfont (3)
            call symbol(x0-0.28*ht,y0+1.1*rad,0.7*ht,'N',0.0,1)
            call gfont (3)
            call plot(x0,y0+0.95*rad,3)
            call plot(x0,y0+1.05*rad,2)
            call plot(x0,y0+1.05*rad,3)
        endif
        call gwidth(0.5*wid)
c-----
c       put on the focal mechanism type
c-----
        call gfont(7)
        ht = 0.14*rad
        if(fmtype)then
            call symbol(x0+0.7*rad,y0+0.8*rad,ht,sym(itype),0.0,3)
        endif
c-----
c       annotate with titles
c-----
        yl = 1.15*rad
c-----
c       center string, but note that there is always a blank space
c       between characters, and that one character is 0.6 Ht wide and
c       space is 0.4 ht
c-----
        xpos = x0-0.5*ltit*hht +0.2*hht
        if(nullit)then
            call newpen(0)
            ipatx = 0
            ipaty = 0
            xlen = 0.01
            ylen = 0.01
            ymin = y0+yl+hht - 0.8*hht
            ymax = y0+yl+hht + 1.2*hht
            xmin = x0-0.5*ltit*hht -0.5*hht
            xmax = x0+0.5*ltit*hht +0.5*hht
            if(hastit)then
            call shader(xmin,ymin,xmax,ymax,
     1              ipatx,ipaty,xlen,ylen)
            endif
            ymin = y0+yl-hht - 0.6*hht
            ymax = y0+yl-hht + 1.2*hht
            if(hassub)then
            call shader(xmin,ymin,xmax,ymax,
     1              ipatx,ipaty,xlen,ylen)
            endif
            call newpen(1)
        endif
c-----
c       put up the titles
c-----
        if(hastit)then
            xpos = x0-0.5*ltit*hht +0.2*hht
            call symbol(xpos,y0+yl+hht,hht,title ,0.0,ltit)
        endif
        if(hassub)then
            xpos = x0-0.5*lsub*hht +0.2*hht
            call symbol(xpos,y0-yl-hht,hht,subtit,0.0,lsub)
        endif
c-----
c       put in P and T symbols
c-----

        if(dofmtp)then
        if((itype.eq.4 .and. pltmch) .or. doptax)then
            if(doptlabel)then
            call plpole(ev(ilg),z(1,ilg),z(2,ilg),z(3,ilg),'T',hemis)
            call plpole(ev(ism),z(1,ism),z(2,ism),z(3,ism),'P',hemis)
            else
            call plpole(ev(ilg),z(1,ilg),z(2,ilg),z(3,ilg),' ',hemis)
            call plpole(ev(ism),z(1,ism),z(2,ism),z(3,ism),' ',hemis)
            endif
        endif
        endif
c-----
c       plot focal circle here, so that it overwrites the title
c       string
c-----
        call gwidth(wid)
        call circle(x0,y0,rad,nullit)

        if(pltmch .and.donodal)then
c-----
c       go through focal mechanism options - 
c           they are not necessarily all exclusive
c-----
c-----
c       plot in the regions of positive functional value
c-----
        cval = 0.0
        if(itype.ge.1 .and. (itype.le.3 .or. itype.eq.5)
     1       .and. .not. ifm)then
            if(abs(vmin).gt.vmax)vmax = abs(vmin)
            if(fmfill.gt.0)then
                call conshd(val,tval,1,NTR,npl,cval,NTR,kval)
            else if(fmfill.lt.0)then
                call consh1(val,tval,1,NTR,npl,vmax,NTR)
            endif
        endif
c-----
c       put in contours of equal relative amplitude on the focal sphere,
c       from 0.75 to -0.75 of the maximum in increments of 0.25
c       
c       Otherwise just plot the zero contour
c-----
        if(fmamp .and. .not. ifm)then
            do 3000 ival=-3,3
                cval = 0.25*(ival)
                if(kolor)then
                    if(ival.lt.0)then
                        call newpen(4)
                    else if(ival.gt.0)then
                        call newpen(2)
                    endif
                else
                    call newpen(kval)
                endif
                call contur(val,tval,2,NTR,npl,cval,NTR)
                if(kolor)call newpen(1)
 3000       continue
        else
           if(donode)then
            cval = 0.0
            call contur(val,tval,2,NTR,npl,cval,NTR)
           endif
        endif
        endif
c-----
c       plot first motion data that is contained in the data file
c-----
        if(ifm)then
            call gwidth(wid)
            call circle(x0,y0,rad,.false.)
            ht = 0.14*rad
            if(itype.eq.1 )then
                call gfont (3)
                call symbol(x0-0.28*ht,y0+1.1*rad,0.7*ht,'N',0.0,1)
                call gfont (3)
            endif
            call gwidth(wid*0.5)
            if(pltmch)call contur(val,tval,2,NTR,npl,0.0,NTR)
            rewind 2
            isump = 0
            isumd = 0
            isumc = 0
            isumi = 0
            WRITE(LOT,*)'P-WAVE FIRST MOTION DATA'
            WRITE(LOT,'(a)')'  TR   PL   AZ   Io   Pol    Amp'
 1234       continue
c-----
c       P-wave first motion information
c       TREND   TAKEOFF_ANGLE   Polarity
c
c       TREND   - angle measured from north (azimuth)
c       TAKEOFF_ANGLE - measured upward from downward vertical.
c             0 degrees is a ray going straight down
c            90 degrees is a ray going horizontal
c           180 degrees is a ray going upward
c       POLARITY
c           +1 = compression  -> octogon (circle)`
c           -1 = dilatation   -> triangle
c           other positive compression( e.g., +2) -> + (plus)
c           other negative dilatation ( e.g., -2) -> - (minus)
c           0 = X
c-----
                read(2,*,end=1235,err=1235)tr,pi0,ip
                trsv = tr
                pisv = pi0
                if(.not.hemis)then
                    pi0 = 180.0 - pi0
                endif
                if(pi0.gt.90.0)then
                    tr0 = tr + 180.0
                    pl0 =  90.0 - pi0
                    pl0 = - pl0
                else
                    tr0 = tr
                    pl0 = 90.0 - pi0
                endif
C        write(6,*)trsv,pisv,tr,pi0,ip,tr0,pl0
                call pplot(1.0,1.0,xx,yy)
                pl = pl0 * 3.1415927/180.0
                tr = tr0 * 3.1415927/180.0
C                call getmtn(tr,pl,gamma,phi,theta,mom,
C     1              sump,sumsv,sumsh,
C     1              domom,f)
                if(hemis)then
                call getmtn(tr,pl,gamma,phi,theta,mom,sump,sumsv,sumsh,
     1              domom,f)
                else
                call getmtn(tr,-pl,gamma,phi,theta,mom,sump,sumsv,sumsh,
     1              domom,f)
                endif
                if(ip.gt.0)then
                    isump = isump + 1
                    if(ip.eq.1)then
                        isymb = 1
                    else
                        isymb = 3
                    endif
                    if(sump.lt.0.0)then
                        isumi = isumi + 1
                    else 
                        isumc = isumc + 1
                    endif
                    call symbol(xx,yy,0.07*rad,char(isymb),0.0,-1)
                else if(ip.lt.0)then
                    isumd = isumd +1
                    if(ip.eq.-1)then
                        isymb = 2
                    else
                        isymb = 29
                    endif
                    if(sump.gt.0.0)then
                        isumi = isumi + 1
                    else 
                        isumc = isumc + 1
                    endif
                    call symbol(xx,yy,0.07*rad,char(isymb),0.0,-1)
                else if(ip.eq.0)then
                    isymb = 4
                    call symbol(xx,yy,0.07*rad,char(isymb),0.0,-1)
                endif
                sgn = ip*sump
                if(sgn .gt. 0)then
                    write(LOT,'(4f5.0,i5,1x,g10.3)')
     1              mod(tr0,360.0),pl0,
     1              trsv,pisv,ip,sump
                else
                    write(LOT,'(4f5.0,i5,1x,g10.3,a)')
     1              mod(tr0,360.0),pl0,
     1              trsv,pisv,ip,sump,' Inconsistent'
                endif
            goto 1234
 1235       continue
            if(verbos)then
            call symbol(x0-0.5*rad,y0+1.8*rad,0.8*ht,'+=',0.0,2)
            call number(999.0,y0+1.8*rad,0.8*ht,real(isump),0.0,-1)
            call symbol(x0-0.5*rad,y0+1.6*rad,0.8*ht,'-=',0.0,2)
            call number(999.0,y0+1.6*rad,0.8*ht,real(isumd),0.0,-1)
            call symbol(x0-0.5*rad,y0+1.4*rad,0.8*ht,'C=',0.0,2)
            call number(999.0,y0+1.4*rad,0.8*ht,real(isumc),0.0,-1)
            call symbol(x0-0.5*rad,y0+1.2*rad,0.8*ht,'I=',0.0,2)
            call number(999.0,y0+1.2*rad,0.8*ht,real(isumi),0.0,-1)
            endif
            close(2)
        endif
c-----
        close(1)
        call pend()
        call pinitf('LUNE.PLT')
        if(eqarea)then
              project = 'equal'
        else
              project = 'stereo'
        endif
        call lunegrid(x0,y0,rad,project)
        call dopt(x0,y0,rad,sngl(ev(1)),sngl(ev(2)),sngl(ev(3)),' ', 
     1       'red',project, 'LB')
c-----
c     annotate certain points
c-----
        call dopt(x0,y0,rad,3.0,1.0,1.0,'+Crack', 
     1       'black',project, 'RC')
        call dopt(x0,y0,rad,-3.0,-1.0,-1.0,'-Crack', 
     1       'black',project, 'LC')
        call dopt(x0,y0,rad,1.0,1.0,1.0,'+ISO', 
     1       'black',project, 'MB')
        call dopt(x0,y0,rad,-1.0,-1.0,-1.0,'-ISO', 
     1       'black',project, 'MT')
        call dopt(x0,y0,rad,2.0,-1.0,-1.0,'+CLVD', 
     1       'black',project, 'RC')
        call dopt(x0,y0,rad,-2.0,1.0,1.0,'-CLVD', 
     1       'black',project, 'LC')
        call pend()
        end

        subroutine getmtn(tr,pl,gamma,phi,theta,mom,sump,sumsv,sumsh,
     1              domom,f)
        real gamma(3), theta(3), phi(3), mom(3,3)
        real f(3)
        logical domom
                coss = cos(tr)
                sins = sin(tr)
                cosd = cos(pl)
                sind = sin(pl)
                gamma(1) = coss*cosd
                gamma(2) = sins*cosd
                gamma(3) = sind
                theta(1) = sind*coss
                theta(2) = sind*sins
                theta(3) = -cosd
                phi(1) = -sins
                phi(2) = coss
                phi(3) = 0.0
                sump = 0.0
                sumsv = 0.0
                sumsh = 0.0
                if(domom)then
                do 1200 i=1,3
                    do 1201 j=1,3
                        sump=sump+gamma(j)*gamma(i)*
     1                      mom(i,j)
                        sumsv=sumsv+theta(i)*gamma(j)*
     1                      mom(i,j)
                        sumsh=sumsh+phi(i)*gamma(j)*
     1                      mom(i,j)
 1201               continue
 1200           continue
                else
                    do 1301 j=1,3
                        sump=sump+gamma(j)*f(j)
                        sumsv=sumsv+theta(j)*f(j)
                        sumsh=sumsh+phi(j)*f(j)
 1301               continue
                endif
        return
        end

        subroutine contur(f,g,i1,imax,jmax,val,npts)
        integer LIN, LOT, LER
        parameter (LIN=5, LER=0, LOT=6)
        integer  NTR
        parameter (NTR=90)
        real*4 f(NTR,2),g(NTR)
c-----
c       contouring program
c           plot a contour of level val from array f(imax,jmax)
c-----
        do  jj = 2,jmax
            call getf(f,g,jj,i1,imax,npts)
             do i = i1,imax
c-------
c
c            F(i-1,j-1)      F(i-1,j)
c
c            F(i,j-1)    F(i,j)
c
c-----
c-----
c            to save memory and to permit this to
c            work in a non-virtual machine, the  function f(time,freq)
c            are stored sequentially on UNIT 01. At any one time we
c            need only access 2 columns, so we use j = 2
c-----
                 j = 2
c-----
c            get center value, assuming Laplacian = 0
c-----
                 cval = (F(i,j)+F(i-1,j)+F(i-1,j-1)+F(i,j-1))/4.0
c-----
c            now sequentially process the four triangles
c-----
                 call lookt(i,jj,i  ,j  ,i-1,j  ,val,cval,F)
                 call lookt(i,jj,i-1,j  ,i-1,j-1,val,cval,F)
                 call lookt(i,jj,i-1,j-1,i  ,j-1,val,cval,F)
                 call lookt(i,jj,i  ,j-1,i  ,j  ,val,cval,F)
             enddo
        enddo
        return
        end

        subroutine lookt(ii,jj,i1 ,j1 ,i2,j2 ,val,cval,F)
        integer  NTR
        parameter (NTR=90)
        real*4 f(NTR,2)
        real*4 xc(3), yc(3)
        logical inarea
        logical interp
        nc = 0
c-----
c       sequentially go through corners of triangle
c-----
        inarea = interp(F(i1,j1),F(i2,j2),val,pos)
        if(inarea)then
            nc = nc + 1
            xc(nc) =  (i1 + (i2-i1)*pos)
            yc(nc) =  (j1 + (j2-j1)*pos) + (jj -2)
        endif
        inarea = interp(cval,F(i2,j2),val,pos)
        if(inarea)then
            nc = nc + 1
            xc(nc) =  (ii-0.5 + (i2-(ii-0.5))*pos)
            yc(nc) =  (1.5 + (j2-1.5)*pos) + (jj -2)
        endif
        inarea = interp(F(i1,j1),cval,val,pos)
        if(inarea)then
            nc = nc + 1
            xc(nc) =  (i1 + ((ii-0.5)-i1)*pos)
            yc(nc) =  (j1 + (1.5-j1)*pos) + (jj -2)
        endif
        ipen = 3
        iret = 1
        if(nc.gt.1)then
            do 100 i=1,nc
            call pplot(xc(i),yc(i),x,y)
            if(iret.gt.0)then
                call plot(x,y,ipen)
                ipen = 2
            endif
  100       continue
        endif
        return
        end

        function interp(y0,y1,val,x)
        logical interp
        real y1,y0,val,x
        denom = y1-y0
        if(denom.eq.0.0 .and. val.ne.y0)then
            interp = .false.
        elseif(denom.eq.0.0 .and. val.eq.y0)then
            interp = .true.
            x = 0.0
        else
            x = (val-y0)/denom
            if(x.lt.0.0 .or. x.gt.1.0)then
                interp = .false.
            else
                interp = .true.
            endif
        endif
        return
        end

        subroutine getf(f,g,jj,n1,n,npts)
      integer LIN, LOT, LER
      parameter (LIN=5, LER=0, LOT=6)
      integer  NTR
      parameter (NTR=90)
        real*4 f(NTR,2),g(NTR)
c-----
c       get the traces envelopes from UNIT 01
c       then update the array f(4096,2) so that column 1 refers
c       always to the leftmost column, even after updating
c-----
c-----
c       use special care for the first entry, jj = 2
c-----
        if(jj.eq.2)then
            rewind 1
            read(1)(g(j),j=1,npts)
            DO 5009 KK = 1,npts
                f(KK,2) = g(KK)
 5009       continue
        endif
c-----
c       normal processing, no check for EOF on read
c
c       1. read in new column
c       2. determine maximum amplitude of new column
c       3. normalize new column, put in f(KK,2)
c          move old f(KK,2) to f(KK,1)
c-----
            read(1)(g(j),j=1,npts)
            DO 5011 KK = 1,npts
                temp = f(KK,2)
                f(KK,2) = g(KK)
                f(KK,1) = temp
 5011       continue
        return
        end

        subroutine pplot(xi,yi,xx,yy)
        real*4 xx,yy
        common/pltcon/x0,y0,rad,eqarea
        logical eqarea
        common/pltmap/tr0,dtr,pl0,dpl
c-----
c       This routine converts an xi, yi coordinate to equivalent value
c       of trend and plunge, which are then projected onto a plane
c       using an equal area or stereographic projection
c-----
c       xi  R*4 -
c       yi  R*4 -
c       tr0 R*4 - value of trend corresponding to index xi = 1.0
c       dtr R*4 - trend increment
c       pl0 R*4 - value of plunge corresponding to index yi = 1.0
c       dpl R*4 - plunge increment
c       xx  R*4 - returned plot coordinate
c       yy  R*4 - returned plot coordinate
c       rad R*4 - radius of sphere
c       x0  R*4 - plot origin offset in x
c       y0  R*4 - plot origin offset in y
c       eqarea  L   - .true. equal area projection
c                   - .false. stereographic projection
c-----
        tr = tr0 + (xi - 1.0)*dtr
        pl = pl0 + (yi - 1.0)*dpl
        tr = tr * 3.1415927/180.0
        pl = pl * 3.1415927/180.0
c-----
c       convert to spherical coordinates on lower unit hemisphere
c-----
        x = cos(tr)*cos(pl)
        y = sin(tr)*cos(pl)
        z = sin(pl)
c-----
c       obtain angle between downward Z axis and vector
c       This is really just pi/2 - pl, but we will compute here
c       because there may be times when just (x,y,z) are given
c-----
        theta = atan2(sqrt(x**2 + y**2), abs(z))
c-----
c       get projection of the point onto the z=0 plane, determining
c       its distance from the origin
c-----
        if(eqarea)then
            r = sqrt(2.0)* sin(theta/2.0)
        else
            r = sin(theta/2.0)/cos(theta/2.0)
        endif
        if(x .ne. 0.0 .and. y .ne. 0.0)then
            alpha=atan2(y,x)
        else
            alpha = 0.0
        endif
c-----
c       get plotting point, noting that plotting system (xx,yy) is with
c       positive xx to right and yy up
c       and the fact that the trend =0 is equal to yy direction
c-----
        yy = y0 + rad * r * cos(alpha)
        xx = x0 + rad * r * sin(alpha)
        return
        end

        subroutine consh1(f,g,i1,imax,jmax,vmax,npts)
c-----
c       plot polarity information as + or - signs, with the amplitude
c       keyed to the amplitude on the focal sphere
c-----
        integer LIN, LOT, LER
        parameter (LIN=5, LER=0, LOT=6)
        integer  NTR
        parameter (NTR=90)
        real*4 f(NTR,2),g(NTR)
        common/pltcon/x0,y0,rad,eqarea
        logical eqarea

        do 1000 jj = 2,jmax
            call getf(f,g,jj,i1,imax,npts)
c-----
c       search for regions of f greater than cval
c       If this condition is met, interpolate and plot
c       We must be careful since the f(,) is a 2-D array, and
c       the first and last j indices do not wrap around, thus
c       we must handle jj=2 case carefully
c-----
            if(mod(jj,2).eq.0)then
            do 1100 ii=1,imax,3
            xi = ii
            xj = jj
            call pplot(xi,xj,xx,yy)
            cal = f(ii,1)
            ssize = 0.05 * rad * abs(f(ii,1))/vmax
            dx = ssize/2.0
            dy = ssize/2.0
c-----
c       Make a plus sign
c_____
            if(cal.gt.0.0)then
                call plot(xx+dx,yy,3)
                call plot(xx-dx,yy,2)
                call plot(xx,yy+dy,3)
                call plot(xx,yy-dy,2)
                call plot(xx,yy,3)
c-----
c       Make a minus sign
c-----
            elseif(cal.lt.0.0)then
                call plot(xx+dx,yy,3)
                call plot(xx-dx,yy,2)
            endif

 1100   continue
        endif
 1000   continue
        return
        end

        subroutine conshd(f,g,i1,imax,jmax,cval,npts,kval)
        integer LIN, LOT, LER
        parameter (LIN=5, LER=0, LOT=6)
        integer  NTR
        parameter (NTR=90)
        real*4 f(NTR,2),g(NTR)
        call newpen(kval)
        do 1000 jj = 2,jmax
            call getf(f,g,jj,i1,imax,npts)
c-----
c       search for regions of f greater than cval
c       If this condition is met, interpolate and plot
c       We must be careful since the f(,) is a 2-D array, and
c       the first and last j indices do not wrap around, thus
c       we must handle jj=2 case carefully
c-----
            if(jj.eq.2)then
                call pltpos(f,1,i1,jj-1,imax,cval)
            endif
            call pltpos(f,2,i1,jj,imax,cval)
 1000   continue
        call newpen(1)
        return
        end

        subroutine pltpos(f,j,i1,jj,imax,cval)
        integer  NTR
        parameter (NTR=90)
        real*4 f(NTR,2)
        logical interp
        yi = jj
        xi = i1
        ipen = 3
        call pplot(xi,yi,xold,yold)
        call plot(xold,yold,ipen)
        vold = f(i1,j)
        do 1100 i = i1+1,imax
            vnew = f(i,j)
            xi = i
            if(vold.ge.cval .and. vnew.ge.cval)then
                call pplot(xi,yi,xx,yy)
                call plot(xold,yold,3)
                call plot(xx,yy,2)
            else if(vold.lt. cval .and. vnew.le.cval)then
                call pplot(xi,yi,xx,yy)
                call plot(xx,yy,3)
            else if(interp(vold,vnew,cval,x))then
                call pplot(xi+x-1.0,yi,xx,yy)
                if(vold.ge.cval)then
                    call plot(xold,yold,3)
                    call plot(xx,yy,2)
                    call plot(xx,yy,3)
                else
                    call plot(xx,yy,3)
                    call pplot(xi,yi,xx,yy)
                    call plot(xx,yy,2)
                endif
            endif
            xold = xx
            yold = yy
            vold = vnew
 1100   continue
        return
        end

        subroutine circle(x0,y0,rad,nullit)
        real*4 x0, y0, rad
        logical nullit
        degrad = 3.1415927/180.0
c-----
c       if necessary, reset the region behind the focal circle to the
c       background color
c-----
        if(nullit)then
            call newpen(0)
            ipatx=0
            ipaty=0
            xlen=0.01
            ylen=0.01
            x1 = x0
            y1 = y0
            x2 = x0 + rad
            y2 = y0
            do 100 i=0,360,5
                arg = i*degrad
                x3 = x0 + 1.15*rad*cos(arg)
                y3 = y0 + 1.15*rad*sin(arg)
                call shadet(x1,y1,x2,y2,x3,y3,
     1              ipatx,ipaty,xlen,ylen)
                x2 = x3
                y2 = y3
  100       continue
            call newpen(1)
        endif
        ipen = 3
        do 200 i=0,360
            arg = i*degrad
            x = x0 + rad*cos(arg)
            y = y0 + rad*sin(arg)
            call plot(x,y,ipen)
            ipen = 2
  200   continue
        call plot(x,y,3)
        return
        end

        subroutine plpole(ev,xmt1,xmt2,xmt3,sym,hemis)
        real*8 ev, xmt1, xmt2, xmt3
        character sym*1
        logical hemis
        
        common/pltmap/tr0,dtr,pl0,dpl
        common/pltcon/x0,y0,rad,eqarea
        logical eqarea
        real*8 r
c-----
c       plot pole on focal sphere
c-----
c       convert to trend and plunge, lower hemisphere
c-----
        degrad = 3.1415927/180.0
        if(xmt3 .lt. 0.0d+00)then
            xmt1 = - xmt1
            xmt2 = - xmt2
            xmt3 = - xmt3
        endif
        r = xmt1**2 + xmt2**2 + xmt3**2
        r = dsqrt(r)
        xmt1 = xmt1/r
        xmt2 = xmt2/r
        xmt3 = xmt3/r
        pl = 90.0 - sngl(dacos(xmt3))/degrad
        r = xmt1**2 + xmt2**2 
        r = dsqrt(r)
        if(r.le.1.0d-06)then
            tr = 0.0
        else
            tr = sngl(datan2(xmt2,xmt1))/degrad
        endif
        if(.not.hemis)then
            pl =  - pl
            tr = 180 + tr
        endif
        xi = (tr - tr0)/dtr + 1.0
        yi = (pl - pl0)/dpl + 1.0
        call pplot(xi,yi,xx,yy)
        call symfil(xx,yy,0.05*rad)
        ht = 0.12*rad
        if(sym.ne.' ')then
        call symbol(xx+0.6*ht,yy+ht,ht,sym,0.0,1)
        endif
        return
        end

        subroutine symfil(xx,yy,ht)
        real*4 x(9),y(9)
        data x/0.9238,0.3827,-0.3827,-0.9238,-0.9238,-0.3827,
     1      0.3827,0.9238,0.9238/
        data y/0.3827,0.9238,0.9238,0.3827,-0.3827,-0.9238,
     1      -0.9238,-0.3827,0.3827/
c-----
c       put in a filled octagon
c-----
        do 100 i=1,8
            xxx = xx + ht*x(i)
            yyy = yy + ht*y(i)
            xxxx= xx + ht*x(i+1)
            yyyy= yy + ht*y(i+1)
            call newpen(1)
            call shadet(xx,yy,xxx,yyy,xxxx,yyyy,0,0,0.01,0.01)
            call newpen(1000)
            call plot(xxx,yyy,3)
            call plot(xxxx,yyyy,2)
  100   continue
        call newpen(1)
        return
        end

        subroutine gcmdln(eqarea,m,itype,rad,x0,y0,dfile,fmfill,
     1      fmamp,fmtype,verbos,nullit,title,hastit,subtit,
     2      pltmch,hassub,titsiz,kolor,hemis,doptax,domom,f,dofmtp,
     3      donodal,doptlabel,verby,donode,kval)
c-----
c       eqarea  L   - .true. equal area projection
c                 .false. stereographic projection
c       m(3,3)  R*4 - moment tensor
c       itype   I*4 - type of plot
c                   1 = P, 2= SV, 3=SH,. 4=POL, 5=S ampl
c       rad R*4 - radius of projection circle
c       x0  R*4 - x - coordinate of center of projection circle
c       y0  R*4 - y - coordinate of center of projection circle
c       dfile   C*80    - name of data file having P-wave 
c                  first motion data
c                  entries are
c                   az, i0 (0=down, 180 = up), pol ( +1,0,-1),'sta'
c       fmfill  I*4 - > 0 solid fill positive area
c                 = 0 no fill
c                 < 0 fill with +- signs
c       fmamp   L   - .true. plot in contours of 
c                   focal mechanism amplitude
c                   from -0.75 to +0.75 of the maximum 
c                   in increments of 0.25
c       fmtype  L   - .true. annotate with symbol 
c                   `indicating type of plot,
c                   e.g., P, SV, SH, POL or S
c       verbos  L   - .true. annotate mechanism plot containing observed
c                   P -wave data with consistencies, inconsistencies
c       nullit  L   - .true. use shadet() to fill focal circle 
c                   with background. This is useful if this output
c                   is to overlay other output
c       title   C*80    - character string giving title at top of plot
c       hastit  L   - .true. if title is present
c       subtit  C*80    - character string giving title at base of plot
c       hassub  L   - .true. if  subtitle is present
c       pltmch  L   -  .true. plot mechanism, else only first motion
c       titsiz  R*4 - height of title and subtitle
c                   default is 0.28 * radius
c       kolor   L   - .true. use red to indicate 0 line of amplitude
c       hemis   L   - .true. lower hemisphere
c                 .false. upper hemisphere
c       doptax  L   - .true. plot P and T symbols on focal sphere
c       domom   L   - .true. Do moment tensor plot, else do force
c       f(3)    R*4 -  point force components
c       dofmtp    L
c       donodal L   - .true. plot nodal planes -TP sets false
c       doptlabel L - .true. put P and T on P T axes
c       verby   L   - .true. output more debug information
c       donode  L   - .true. node node curve
c-----
        integer LOT
        parameter (LOT=6)
        
        character names*80, dfile*80, title*80, subtit*80
        real*4 m(3,3), f(3)
        integer mnmarg
        logical eqarea
        integer fmfill
        logical fmamp
        logical fmtype
        logical verbos
        logical nullit
        logical hastit
        logical hassub
        logical pltmch
        logical kolor
        logical hemis
        logical doptax
        logical domom
        logical dofmtp
        logical donodal
        logical doptlabel
        logical verby
        logical issymmetric
        logical donode
        integer kval

        eqarea = .true.
        fmfill = 0
        fmamp = .false.
        fmtype = .false.
        verbos = .false.
        nullit = .false.
        hastit = .false.
        hassub = .false.
        pltmch = .false.
        kolor = .true.
        hemis = .true.
        doptax = .false.
        domom = .true.
        dofmtp = .false.
        donodal = .true.
        doptlabel = .true.
        verby = .false.
        donode = .true.

        dfile = ' '
        title = ' '
        subtit = ' '
        itype = 1
        rad = 2.0
        x0 = 4.0
        y0 = 4.0
        imom = 0
        stk = 0
        dip = 0
        rake = 0
        xmom = 1.0
        titsiz = -1.0
        kval = 1
        nmarg = mnmarg()
        f(1) = 0.0
        f(2) = 0.0
        f(3) = 0.0
        m(1,1) = 0.0
        m(1,2) = 0.0
        m(1,3) = 0.0
        m(2,1) = 0.0
        m(2,2) = 0.0
        m(2,3) = 0.0
        m(3,1) = 0.0
        m(3,2) = 0.0
        m(3,3) = 0.0
        issymmetric = .true.
        isds = -1
        i = 0
 1000   continue
        i = i + 1
        if(i.gt.nmarg)go to 9999
        call mgtarg(i,names)
c-----
c       source specification
c-----
            if(names(1:2).eq.'-D')then
                i=i+1
                imom = 1
                call mgtarg(i,names)
                call chtofp(names,dip)
                isds = 0
            else if(names(1:2).eq.'-S' .and. names(1:3).ne.'-SH'
     1          .and. names(1:3).ne.'-SV')then
                i=i+1
                imom = 1
                call mgtarg(i,names)
                call chtofp(names,stk)
                isds = 0
            else if(names(1:2).eq.'-R'.and.names(1:3).ne.'-RA')then
                i=i+1
                imom = 1
                call mgtarg(i,names)
                call chtofp(names,rake)
                isds = 0
c-----
c       moment for fault source
c-----
            else if(names(1:3).eq.'-M0' .or. names(1:3).eq.'-MO')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,xmom)
            else if(names(1:3).eq.'-MW' .or. names(1:3).eq.'-Mw')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,xmom)
                xmom = 10.**(1.5*xmom + 16.10)
c-----
c       moment tensor
c-----
            else if(names(1:3).eq.'-xx' .or. names(1:3).eq.'-XX')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,m(1,1))
                isds = 1
            else if(names(1:3).eq.'-yy' .or. names(1:3).eq.'-YY')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,m(2,2))
                isds = 1
            else if(names(1:3).eq.'-zz' .or. names(1:3).eq.'-ZZ')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,m(3,3))
                isds = 1
            else if(names(1:3).eq.'-xy' .or. names(1:3).eq.'-XY')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,m(1,2))
                isds = 1
            else if(names(1:3).eq.'-xz' .or. names(1:3).eq.'-XZ')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,m(1,3))
                isds = 1
            else if(names(1:3).eq.'-yz' .or. names(1:3).eq.'-YZ')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,m(2,3))
                isds = 1
            else if(names(1:3).eq.'-yx' .or. names(1:3).eq.'-YX')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,m(2,1))
                isds = 1
            else if(names(1:3).eq.'-zx' .or. names(1:3).eq.'-ZX')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,m(3,1))
                isds = 1
            else if(names(1:3).eq.'-zy' .or. names(1:3).eq.'-ZY')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,m(3,2))
                isds = 1
            else if(names(1:2) .eq. '-F' .or. names(1:2).eq.'-f')then
c-----
c               point force
c-----
                if(names(1:3).eq.'-fx'.or.names(1:3).eq.'-FX')then
                    i=i+1
                    call mgtarg(i,names)
                    call chtofp(names,f(1))
                    isds = 2
                else if(names(1:3).eq.'-fy'.or.names(1:3).eq.'-FY')then
                    i=i+1
                    call mgtarg(i,names)
                    call chtofp(names,f(2))
                    isds = 2
                else if(names(1:3).eq.'-fz'.or.names(1:3).eq.'-FZ')then
                    i=i+1
                    call mgtarg(i,names)
                    call chtofp(names,f(3))
                    isds = 2
                else if(names(1:3).eq.'-F1'.or.names(1:3).eq.'-f1')then
                    domom = .false.
                    i = i + 1
                    call mgtarg(i,names)
                    read(names,'(bn,f20.0)')f(1)
                    isds = 2
                else if(names(1:3).eq.'-F2'.or.names(1:3).eq.'-f2')then
                    domom = .false.
                    i = i + 1
                    call mgtarg(i,names)
                    read(names,'(bn,f20.0)')f(2)
                    isds = 2
                else if(names(1:3).eq.'-F3'.or.names(1:3).eq.'-f3')then
                    domom = .false.
                    i = i + 1
                    call mgtarg(i,names)
                    read(names,'(bn,f20.0)')f(3)
                    isds = 2
c-----
c               specific plot type, e.g., plus/minus, solid
c-----
                else if(names(1:4).eq.'-FMF' .or. 
     1                      names(1:4).eq.'-fmf')then
                    fmfill = 1
                else if(names(1:4).eq.'-FMP' .or. 
     1                      names(1:4).eq.'-fmp')then
                    fmfill = -1
                else if(names(1:4).eq.'-FMA' .or. 
     1                      names(1:4).eq.'-fma')then
                    fmamp = .true.
                else
c-----
c                    first motion data file
c-----
                    i = i + 1
                    call mgtarg(i,dfile)
                endif
c-----
c       plot type
c-----
            else if(names(1:2).eq.'-P')then
                itype = 1
            else if(names(1:3).eq.'-SV')then
                itype = 2
            else if(names(1:3).eq.'-SH')then
                itype = 3
            else if(names(1:4).eq.'-pol')then
                itype = 4
c-----
c       plot projection
c-----
            else if(names(1:3).eq.'-st')then
                    eqarea = .false.               
            else if(names(1:3).eq.'-eq')then
                    eqarea = .true.               
c-----
c       plot positioning
c-----
            else if(names(1:3).eq.'-RA')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f20.0)')rad
            else if(names(1:3).eq.'-X0')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f20.0)')x0
            else if(names(1:3).eq.'-Y0')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f20.0)')y0
c-----
c       plot annotation
c-----
            else if(names(1:4).eq.'-ANN' .or. names(1:4).eq.'-ann')then
                fmtype = .true.
            else if(names(1:2).eq.'-Z')then
                nullit = .true.
            else if(names(1:3).eq.'-TT')then
                hastit = .true.
                i = i + 1
                call mgtarg(i,title)
            else if(names(1:3).eq.'-TB')then
                hassub = .true.
                i = i + 1
                call mgtarg(i,subtit)
            else if(names(1:3).eq.'-TS')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f20.0)')titsiz
            else if(names(1:4).eq.'-NOS')then
                issymmetric = .false.
            else if(names(1:4).eq.'-NON')then
                donode = .false.
            else if(names(1:3).eq.'-NM')then
                pltmch = .false.
            else if(names(1:7).eq.'-NOFMTP')then
                dofmtp = .false.
            else if(names(1:2).eq.'-K'.and.names(1:3).ne.'-KC')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,i10)')kval
            else if(names(1:3).eq.'-KC')then
                kolor = .true.
            else if(names(1:2).eq.'-U' .or. names(1:2).eq.'-u')then
                hemis = .false.
            else if(names(1:3).eq.'-tp')then
                dofmtp = .true.
                doptax = .true.
                doptlabel = .true.
            else if(names(1:3).eq.'-TP')then
                dofmtp = .true.
                doptax = .true.
                donodal = .false.
            else if(names(1:3).eq.'-tP')then
                doptax = .true.
                dofmtp = .true.
                donodal = .false.
                doptlabel = .false.
            else if(names(1:2).eq.'-V')then
                verby = .true.
                verbos =.true.
            else if(names(1:2).eq.'-?')then
                call usage()
            else if(names(1:2).eq.'-h')then
                call usage()
            endif

        go to 1000
 9999   continue
c-----
c       if a stk, dip, rake force or moment tensor 
c       specified plot nodal planes
c-----
        if(isds.ge.0)then
            pltmch = .true.
        endif
        if(domom)then
c-----
c           if strike, dip, rake specified, calculate
c           double couple moment tensor
c-----
            if(imom.gt.0)then
                call tensor(stk,dip,rake,xmom,m)
            endif
c-----
c           force symmetry of moment tensor - this is valid for
c           double couples. 
c           If the field of a single couple is desired, do not force
c           this symmetry
c-----
            if(issymmetric)then
            m(3,1) = m(1,3)
            m(3,2) = m(2,3)
            m(2,1) = m(1,2)
            endif
            write(LOT,*)'  M(1,1)    =',m(1,1)
            write(LOT,*)'  M(1,2)    =',m(1,2)
            write(LOT,*)'  M(1,3)    =',m(1,3)
            write(LOT,*)'  M(2,1)    =',m(2,1)
            write(LOT,*)'  M(2,2)    =',m(2,2)
            write(LOT,*)'  M(2,3)    =',m(2,3)
            write(LOT,*)'  M(3,1)    =',m(3,1)
            write(LOT,*)'  M(3,2)    =',m(3,2)
            write(LOT,*)'  M(3,3)    =',m(3,3)
        else
c-----
c           normalize the force
c-----
            
            write(LOT,*)' Force Specification ' 
            write(LOT,*)'  F(1)      =',f(1)
            write(LOT,*)'  F(2)      =',f(2)
            write(LOT,*)'  F(3)      =',f(3)
            fnorm = f(1)*f(1) + f(2)*f(2) + f(3)*f(3)
            fnorm = sqrt(fnorm)
            if(fnorm.gt.0.0)then
                f(1) = f(1)/fnorm
                f(2) = f(2)/fnorm
                f(3) = f(3)/fnorm
            endif
            
        endif
        if(titsiz.le.0.0)titsiz = 0.28*rad
        return
        end

        subroutine tensor(stk,dip,rake,xmom,m)
c-----
c       calculate moment tensor for a double couple mechanism
c-----
c       stk - R*4   - strike, measured clockwise from north
c                 when looking down, e.g., east = 90
c       dip - R*4   - dip, measured from horizontal when looking
c                 in direction of strike, fault dips down to right
c                 0 <= dip <= 90
c       rake    - R*4   - rake. Measured from horizontal with respect
c                 to strike. If -180 < rake <= 180 is taken as
c                 the convention, then a negative rake, implies 
c                 that the P-wave first motion in the center of the
c                 focal sphere is negative (e.g., dilatation)
c       xmom    - R*4   - moment value
c       m(3,3)  - R*4   - moment tensor
c-----
        integer LOT
        parameter (LOT=6)
        real*4 stk, dip, rake, xmom, m(3,3)
            degrad=0.0174532925
            tol = 1.0e-7
            dp = degrad*dip
            st = degrad*stk
            sl = degrad*rake
            sind=sin(dp)
            cosd=cos(dp)
            sinr=sin(sl)
            cosr=cos(sl)
            sins=sin(st)
            coss=cos(st)
            sin2d=sin(2.*dp)
            cos2d=cos(2.*dp)
            sin2s=sin(2.*st)
            cos2s=cos(2.*st)
            m(1,1)=(-sind*cosr*sin2s-sin2d*sinr*sins*sins)*xmom
            m(2,2)=(sind*cosr*sin2s-sin2d*sinr*coss*coss)*xmom
            m(3,3)=(sin2d*sinr)*xmom
            m(1,2)=(sind*cosr*cos2s+0.5*sin2d*sinr*sin2s)*xmom
            m(1,3)=(-cosd*cosr*coss-cos2d*sinr*sins)*xmom
            m(2,3)=(-cosd*cosr*sins+cos2d*sinr*coss)*xmom
            m(2,1) = m(1,2)
            m(3,1) = m(1,3)
            m(3,2) = m(2,3)
c-----
c       clean up small values
c-----
        xmax=-1.0e+37
        do 4 i=1,3
            do 5 j=1,3
                if(abs(m(i,j)).gt.xmax)xmax = abs(m(i,j))
    5       continue
    4   continue
        thresh = tol * xmax
        do 6 i=1,3
            do 7 j=1,3
                if(abs(m(i,j)).lt.thresh) m(i,j) = 0.0
    7       continue
    6   continue

c-----
c       write out the information
c-----
        write(LOT,*)' STRIKE    =',stk
        write(LOT,*)' DIP       =',dip
        write(LOT,*)' RAKE      =',rake
        write(LOT,*)' MOMENT    =',xmom
        write(LOT,*)' Equivalent Moment Tensor' 
        return
        end

        subroutine pltpol(ip,it,pol,tr)
c-----
c       plot polarization angle
c-----
        common/pltcon/x0,y0,rad,eqarea
        logical eqarea
            xi = it
            yi = ip
            call pplot(xi,yi,xx,yy)
            ang = tr + pol 
            c = cos(ang)
            s = sin(ang)
            xxx = 0.05*rad*s
            yyy = 0.05*rad*c
            call plot(xx+xxx,yy+yyy,3)
            call plot(xx-xxx,yy-yyy,2)
            call plot(xx-xxx,yy-yyy,3)
        return
        end
        function lgstr(str)
c-----
c      function to find the length of a string
c      this will only be used with file system path names
c      thus the first blank 
c      indicates the end of the string
c-----
        character*(*) str
        integer lgstr
        n = len(str)
        lgstr = 1
        do 1000 i=n,1,-1
            lgstr = i
            if(str(i:i).ne.' ')goto 100
 1000   continue
  100   continue
        return
        end
      
        subroutine usage()
        parameter (LOT=6, LER=0)
        write(LER,*)
     1  'fmplot -eq -st -XX Mxx -XY Mxy -XZ Mxz -YY Myy'
        write(LER,*)
     1  '       -YZ Myz -ZZ Mzz -P -SV -SH -S -pol -tp -TP -tP '
        write(LER,*)
     1  '       -F1 f1 -F2 f2 -F3 f3'
        write(LER,*)
     1  '       -S stk -D dip -R rake -MOM mom -MW Mw'
        write(LER,*)
     1  '       -RAD rad -X0 x0 -Y0 y0 '
        write(LER,*)
     1  '       -FMFILL -FMAMP -FMPLMN -F fmfile -ANN -V -Z '
        write(LER,*)
     1  '       -TT title -TB subtitle -TS titlesize -NM -UP -?'
        write(LER,*)
     1  ' -NOSYMMETRIC (default false) single couple'
        write(LER,*)
     1  '       -NONODAL (default false) Do not draw nodal plane'
        write(LER,*)
     1  ' -eq          Equal area projection (default)'
        write(LER,*)
     1  ' -st          Stereographic projection'
        write(LER,*)
     1  ' -XX  Mxx    (1,1) component of moment tensor'
        write(LER,*)
     1  ' -XY  Mxy    (1,2) component of moment tensor'
        write(LER,*)
     1  ' -XZ  Mxz    (1,3) component of moment tensor'
        write(LER,*)
     1  ' -YY  Myy    (2,2) component of moment tensor'
        write(LER,*)
     1  ' -YZ  Myz    (2,3) component of moment tensor'
        write(LER,*)
     1  ' -ZZ  Mzz    (3,3) component of moment tensor'
        write(LER,*)
     1  ' -P          P-wave display'
        write(LER,*)
     1  ' -SV         SV-wave display'
        write(LER,*)
     1  ' -SH         SH-wave display'
        write(LER,*)
     1  ' -pol        S-wave polarization angle'
        write(LER,*)
     1  ' -tp         put T- and P- axes on plot'
        write(LER,*)
     1  ' -TP         only plot TP, no nodal planes'
        write(LER,*)
     1  ' -tP         only plot TP, no nodal planes, no TP label'
        write(LER,*)
     1  ' -NOFMTP  (default false) '
        write(LER,*)
     1  ' -RAD rad    Radius of circle (default 2.0 in)'
        write(LER,*)
     1  ' -X0 x0      x-coordinate of center of circle (default 4.0 in)'
        write(LER,*)
     1  ' -Y0 y0      y-coordinate of center of circle (default 4.0 in)'
        write(LER,*)
     1  ' -S          Strike of fault plane'
        write(LER,*)
     1  ' -D          Dip of fault plane'
        write(LER,*)
     1  ' -R          Rake or rake angle on plane'
        write(LER,*)
     1  ' -MOM Mom    Seismic moment in dyne-cm )default 1.0)'
        write(LER,*)
     1  ' -MW Mw      Moment Magnitude'
        write(LER,*)
     1  ' -FMFILL     Solid Fill region of positive amplitude',
     2      '    (default = .false.)'
        write(LER,*)
     1  ' -FMPLMN     Fill region with +- signs related to',
     2  '    amplitude  (default = false)'
        write(LER,*)
     1  ' -FMAMP      Display amplitude contour ',
     2          '    (default = .false.)'
        write(LER,*)
     1  ' -F file    file contains P-wave first motion data'
        write(LER,*)
     1  '     Trend Takeoff-angle ID '
        write(LER,*)
     1  '        where ID =+-1 > Circle/triangle,',
     2  '   +-2 -> + or - sign'
        write(LER,*)
     1  ' -ANN or -ann  Annotate plot with type: P, SV, SH, S or pol'
        write(LER,*)
     1  ' -Z    Clears background - useful for 0verlays (default=false)'
        write(LER,*)
     1  ' -TT title   Title above plot'
        write(LER,*)
     1  ' -TB subtitle Title below plot'
        write(LER,*)
     1  ' -TS titlesize (inches. Default=0.28*rad)'
        write(LER,*)
     1  ' -NM         No mechanism only circle and first motion data'
        write(LER,*)
     1  ' -KC         For -FMAMP amplitude plots use red for zero line'
        write(LER,*)
     1  ' -K kval     default black) Use color kval for fill and ',
     1     'nodal planes.'
        write(LER,*)
     1  ' -UP         Upper hemisphere projection (default lower)'
        write(LER,*)
     1  ' -F1 f1 -F2 f2 -F3 f3 (default 0.0) point force'
        write(LER,*)
     1  ' -NOSYMMETRIC (default is symmetric moment tensor'
        write(LER,*)
     1  ' -?          Usage query, but no execution'
        write(LER,*)
     1  ' -h          Usage query, but no execution'

        write(LER,*)' OUTPUT IS IN FILE FMPLOT.PLT'
        stop
        end

        subroutine gtev(xmt,mom,np,ev,ev1,ilg,ism,z)
c-----
c       xmt   - R* symmetric moment tensor 
c       mom   - R symmetric moment tensor R
c       np    - I dimension of mom and xmt
c       ev    - R*8
c       ev1   - R*8
c       ilg   - pointer to largest eigenvalue
c       ism   - pointer to smallest eigenvalue
c       z     - R*8  eigenvectors
c-----
        real*8 xmt(np,np),ev(np),ev1(np)
        real*8 z(3,3)
        real*4 mom(3,3)
        integer ierr
        integer LOT
        parameter (LOT=6)
c-----
c       compute eigenvalues and eigenvectors of moment tensor matrix
c       Get index of largest and smallest eigenvalues
c-----
        do 100 i=1,3
            do 110 j=1,3
                xmt(i,j) = dble(mom(i,j))
  110       continue
  100   continue
        call tred2(np,np,xmt,ev,ev1,z)
        call tql2(np,np,ev,ev1,z,ierr)
        elg = -1.0e+38
        esm =  1.0e+38
        ilg = 1
        ism = 1
        write(LOT,*)' EIGENVALUE, AND EIGENVECTOR OF M(i,j)'
    2   format(' ',e11.4,'  (',e11.4,',',e11.4,',',e11.4,')')
        do 120 j=1,3
            write(LOT,2)ev(j),(z(i,j),i=1,3)
            if(ev(j).gt.elg)then
                elg = ev(j)
                ilg = j
            endif
            if(ev(j).lt.esm)then
                esm = ev(j)
                ism = j
            endif
  120   continue
        return
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
        integer lgstr
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

      subroutine dopt(x0,y0,rad,ol1,ol2,ol3,str,ccolor,project,pos)
      implicit none
      real x0, y0, rad
      real ol1, ol2, ol3
      character str*(*)
      character ccolor*(*),project*(*), pos*(*)
c-----
c     x0, y0 center of projection
c     rad    radius of projection
c     l1, l2, l3 eigenvalues - do not have to be ordered
c     str    string to be plotted
c     ccolor - C*(*) - pen color  'red' 'black' or 'blue'
c     project 'stereo' or 'equal'
c     pos    position of annotation string
c     LT - left top    MT - center top    RT right top
c     LC - left center MC - center center RC right center
c     LB - left bottom MB - center bottom RB right bottom
c       note for use here LC MC RC would overwrite the symbol 
c       and would be ugly
c-----

      real x, y, ht
      real xoff, yoff
      real l1, l2 , l3
      real snorm

       real bet, gam, mnorm
       real size
       real degrad
       DEGRAD=0.01745329251994

      call order(l1,l2,l3,ol1,ol2,ol3)
      

c-----c
c     Tape, W. and Tape, C. (2012). A geometric comparison of source-type plots for 
c     moment tensors, Geophys. J. Int. 190, 499-510
c-----    Equation 4
c   before doing the square root, normalize. This will not affect the
c   determination of angles
C   otherwise overflow in sqrt
c-----
       snorm=amax1(abs(l1),abs(l2),abs(l3))
       l1 = l1 / snorm
       l2 = l2 / snorm
       l3 = l3 / snorm
       mnorm = sqrt(l1*l1 + l2*l2 + l3*l3)
       gam = atan2(-l1 + 2*l2 - l3,sqrt(3.0)*(l1-l3))
       bet = acos( (l1+l2+l3)/(sqrt(3.0)*mnorm))
       gam = gam/DEGRAD
       bet = bet/DEGRAD

       call luneproj(bet,gam,x0,y0,rad,x,y,project)
      if(ccolor.eq.'red')then
         call newpen(2)
         call fillit('CI',0.07,x,y)
         call newpen(1)
         call curvit('CI',0.07,x,y)
      else if(ccolor.eq.'blue')then
         call newpen(4)
         call fillit('CI',0.05,x,y)
         call newpen(1)
         call curvit('CI',0.05,x,y)
      else if(ccolor.eq.'black')then
         call newpen(1)
         call fillit('CI',0.05,x,y)
         call newpen(1)
         call curvit('CI',0.05,x,y)
      endif
       ht = 0.07
       if(pos(1:2).eq.'LT')then
            xoff =  ht + ht
            yoff = -ht -ht
       else if(pos(1:2).eq.'LC')then
            xoff =  ht + ht
            yoff = -0.5*ht
       else if(pos(1:2).eq.'LB')then
            xoff =  ht + ht
            yoff =  ht + ht
       else if(pos(1:2).eq.'MT')then
            xoff = 0.0
            yoff = -ht -ht -0.5*ht
       else if(pos(1:2).eq.'MC')then
            xoff = 0.0
            yoff = 0.5*ht
       else if(pos(1:2).eq.'MB')then
            xoff = 0.0
            yoff =  ht + ht
       else if(pos(1:2).eq.'RT')then
            xoff = -ht -ht
            yoff = -ht -ht
       else if(pos(1:2).eq.'RC')then
            xoff = -ht -ht
            yoff = -0.5*ht
       else if(pos(1:2).eq.'RB')then
            xoff = -ht -ht
            yoff =  ht + ht
       endif
       if(pos(1:1) .eq. 'M')then
            call gcent (x+xoff,y+yoff,ht,str,0.0)
       else if(pos(1:1) .eq. 'R')then
            call gright(x+xoff,y+yoff,ht,str,0.0)
       else if(pos(1:1) .eq. 'L')then
            call gleft( x+xoff,y+yoff,ht,str,0.0)
       endif
       return
       end

      subroutine lunegrid(x0,y0,rad,ctype)
c-----
c     plot the line with [ - pi/6] to [pi/6] in longitude
c      x0, y0   - R  coordinates of center of project
c      rad      - R  radius of project
c      ctype    - C*(*) 'stereo' or 'equal'
c-----
      implicit none
      real x0, y0,rad
      character ctype*(*)

      integer NPTS
      parameter (NPTS=181)
      real x(NPTS), y(NPTS)
      integer igam, ibet
      real bet, gam
      real xx, yy
      real xd, yd, zd


      do igam = -30,30,10
         gam = igam
         do ibet = 1,NPTS
            bet = ibet-1
            call luneproj(bet,gam,x0,y0,rad,xx,yy,ctype)
            x(ibet) = xx
            y(ibet) = yy
          enddo
c-----
c         plot the curve
c-----
             call doplotline(x,y,NPTS,'black','solid',rad)
      enddo
c-----
c     plot the lines of constant beta
c-----
      do ibet=10,170,10
         bet = ibet
         call luneproj(bet,30.0,x0,y0,rad,xx,yy,ctype)
         x(1) = xx
         y(1) = yy
         call luneproj(bet,-30.0,x0,y0,rad,xx,yy,ctype)
         x(2) = xx
         y(2) = yy
         if(ibet.eq.90)then
              call doplotline(x,y,2,'red','solid',rad)
         else
              call doplotline(x,y,2,'black','solid',rad)
         endif
      enddo
      return
      end

      subroutine doplotline(x,y,n,ccolor,ctype,rad)
c-----
c     x   - R array of x values
c     y   - R array of y values
c     n   - I number of points to plot
c     ccolor - C*(*) - pen color  'red' 'black' or 'blue'
c     ctype  - C*(*) - line type  'solid' or 'dashed' or 'dotted'
c     rad - R radius of the lune project - this is used to set
c     the dashed length
c-----
      implicit none
      integer n
      real x(n), y(n), rad
      character ccolor*(*), ctype*(*)

      real xx, yy

      integer ipat
      real xlen
      integer i

      xlen = rad/60.


      if(ccolor.eq.'red')then
         call newpen(2)
      else if(ccolor.eq.'blue')then
         call newpen(4)
      else if(ccolor.eq.'black')then
         call newpen(1)
      endif

      xx = x(1)
      yy = y(1)
c-----
c     lift pen to first point
c-----
      call plot(xx,yy,3)
      do i=1,n
         xx = x(i)
         yy = y(i)
         if(ctype.eq.'solid')then
              call plot(xx,yy,2)
         else if(ctype.eq.'dotted')then
              ipat = 9
              call plotd(xx,yy,ipat,xlen)
         else if(ctype.eq.'dashed')then
              ipat = 7
              call plotd(xx,yy,ipat,xlen)
         endif
      enddo
c-----
c     lift pen at last point
c-----
      call plot(xx,yy,3)
c-----
c     reset pen color
c-----
      call newpen(1)
      return
      end

            

      subroutine luneproj(degcolat,deglon,x0,y0,rad,x,y,project)
c-----
c      Input values
c      degcolat - R  co-latitude, e.g, 0 = top (N) and 180 = bottom(S)
c                    = beta
c      deglon   - R  longitude 
c                    = gamma
c      x0, y0   - R  coordinates of center of project
c      rad      - R  radius of project
c      project  - C*(*) - 'stereo' for equal antlke 
c                          'equal' for equal area
c      Return values
c      x,y      - R  output x,y coordinates for the plot
c-----
       implicit none
       real degcolat, deglon, x0, y0, rad, x, y
       character project*(*)

       real radcolat, radlon
       real degrad
       real k

       real deglat
       DEGRAD=0.01745329251994
     
       radcolat = degcolat*degrad
       radlon   = deglon*degrad

c-----
c     Ammon Mathematica code 
c-----
       deglat = 90 - degcolat
       if(project.eq.'stereo')then
           k = (2.)/(1 + cos(degrad*deglat)*cos(degrad*deglon))
           k = k /2.
           x = x0 + rad*k*cos(degrad*deglat)*sin(degrad*deglon)
           y = y0 + rad*k*sin(degrad*deglat)
       else if(project.eq.'equal')then
           k = (2.)/(1 + cos(degrad*deglat)*cos(degrad*deglon))
           k = k /2.
           k = sqrt(k)
           x = x0 + rad*k*cos(degrad*deglat)*sin(degrad*deglon)
           y = y0 + rad*k*sin(degrad*deglat)
       endif
      
       return
       end

       subroutine order(l1,l2,l3,ol1,ol2,ol3)
       implicit none
       real l1,l2,l3,ol1,ol2,ol3

       real arr(3)
       real tmp
       integer i,j
c-----
c      order the original values ol1,ol2,ol3
c      so that they are l1 >= l2 >=l3
c      Since there are only thee values, do not use sort routine
c-----
       arr(1) = ol1
       arr(2) = ol2
       arr(3) = ol3
C      l1 = amax1(ol1,ol2,ol3)
C      l3 = amin1(ol1,ol2,ol3)
       do i=1,3
          do j=1,3
            if(arr(j).lt.arr(i))then
               tmp = arr(i)
               arr(i) = arr(j)
               arr(j) = tmp
            endif
          enddo
       enddo
     
       l1 = arr(1)
       l2 = arr(2)
       l3 = arr(3)
       return
       end
               
