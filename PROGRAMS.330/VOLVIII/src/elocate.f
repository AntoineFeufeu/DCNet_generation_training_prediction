        program elocate
c-----
c    COMPUTER PROGRAMS IN SEISMOLOGY
c
c    PROGRAM ELOCATE
c
c    COPYRIGHT (C) 1990
c    R. B. Herrmann
c    Department of Earth and Atmospheric Sciences
c    3507 Laclede Avenue
c    St. Louis, MO 63103
c    USA
c    FAX (314) 658-3117
c    TEL (314) 658-3131
c-----
c    Changes
c        08 DEC 2004
c        As per comment of zhixiany@cea-igp.ac.cn 
c            who wanted to fix the location,
c        I changed the initial lat/lon - 
c            if -F and if -LAT and -LON are specified
c        then I DO NOT make tlat = tlat + 0.1 and tlon = tlon + 0.1
c    09 JAN 2005 - emares not defined in subroutine 
c            sumup - now set to -9.99
c    06 DEC 2005 - output GAP as per hypo71 
c    30 JUL 2009 - modified Event Time output to permit event time in 
c              OCAL format for quick cut-and-paste in the gsac command
c              ch ocal 2009 07 29 10 00 37 123
c    17 FEB 2010 - initial listing of phases also shows lat/lon
c              also add line to facility ch evla event_lat evlo event_lon      
c    06 JUN 2010 - added the syntax for ch o gmt YEAR DOY HR MN SEC MSEC
c              for use with sac2000 since sac2000 does not support ch o cal
c              The appears after the output summary
c    30 AUG 2011 - modified elocate.sum for more precision in the depth
c    20 NOV 2014 - added references for confidence regions in errelp()
c    26 NOV 2014 - added -C flag to get covariance matrix output
c                  documented and corrected the determination of the
c                  rotation angle for the confidence ellipse.
c                  Note: the confidence ellipse is 1 sigma and the results
c                  must be multiplied by 
c                  sqrt(2 var F(w,2,N-2))
c    03 FEB 2016 - corrected error in etoh for years before 1970
c                  if(doy .gt. 0 ) go to 2000 replaced by
c                  if(doy .ge. 0 ) go to 2000 replaced by
c    22 MAR 2016 - corrected error ellipse to agree with Flinn.
c                  also checked by real data set. Note that the 
c                  program does not apply the F-statistic to get the
c                  actual confidence limits.
c-----
        character tratbl*8
        real*4 tlat, tlon, z
        real*4 l1, l2
        integer*4 nr
        real*8 t0
        logical query
c-----
c    gcmdln arguments
c-----
        integer model
        logical fixed, verbose, inter,veryverbose
        integer fmask(4), mxiter
        character auth*15

c-----
c    parse command line argument
c-----
        call gcmdln(tlat,tlon,z,t0,model,fixed,verbose,fmask,
     1      auth,mxiter,inter,veryverbose)
        imodel = model

  333   continue
            call modlin(tratbl,model,verbose,inter)
            call setup(nr,verbose)
   97       continue 
                call setinl(nr,t0,tlat,tlon,tratbl,
     1              verbose,fmask,inter)
                if(tratbl.eq.'BEAM')then
                    call beamit(t0,tlat,tlon,nr,z,
     1                  tratbl,dlat,dlon,dz,dt,
     2                  verbose,
     3                  fmask,rms,inter)
                else
                    call locate(t0,tlat,tlon,nr,z,
     1                  tratbl,dlat,dlon,dz,dt,
     2                  l1,l2,theta,fixed,verbose,veryverbose,
     3                  fmask,mxiter,rms,inter)
                endif
                if(verbose .and. inter .and. imodel.lt.0)then
                    if(query('FINISHED WITH DATA SET? (y/n)'))
     1                  go to 93
                    if(query('USE DIFFERENT MODEL? (y/n)'))
     1                  go to 333
                    go to 97
                endif
   93   continue
c-----
c    generate output files
c-----
        call sumup(t0,tlat,tlon,nr,tratbl,z,dlat,dlon,dz,dt,
     1  l1,l2,theta,fmask,fixed,rms)
        end 

        subroutine  errelp(sum2,a11,a12,a22,theta,l1,l2,
     1      verbose) 
c-----
c    given elements of variance-covariance matrix,
c    find the orientation of the errelpe in the
c    deltaz=0 plane. This is not the projection, but rather the
c    intersection. The projection incorporating the third
c    dimension will give a larger ellipse
c
c    Flinn, E. A. (1965). Confidence regions and error determinations for
c    seismic event location, Reviews of Geophysics 3, 157-185
c    [Note: Flinn has the order of the solution vector (latitude, longitude,
c           depth, time)
c
c    Iterative earthquake location entails solving, assuming a stagle inverse
c                T  -1  T
c    DELTA p = (A A)   A R   (9) these are equation numbers in Flinn (1965)
c     where DELTA p would be the change in
c     hypocenter and origin time, and R is the vector of residuals
c                                          2   T  -1
c   The variance-covariance matrix is S = s  (A A)   where    (12)
c      2           -1      2
c     s   = (N - K)   SUM R  = Q/(N-K)                        (11)
c                                                        1/2
c   Estimate of error on change in parameter is s = (S  )     (13)
c                                                k    kk
c                          i
c   Confidence region is
c             T   -1            2
c     (x - x0')  S  (x-x0') = cK         (14) 
c   where the primed  are the least squares solution and the
c   unprimed define the confidence region
c        2
c      cK  = K Q F(w;K,N-K)/(N-K)   (15) and F is F-statistic with N-K
c   degrees of freedom at the 100 w % confidence level
c
c   To estimate the confidence region for L < K, e.g., in the jm plane, consider
c             T  -1           2
c     (x - x0)  E  (x-x0) = c2         (16) 
c   where
c          | Sjj  Sjm |
c      E = |          |                (17)
c          | Smj  Smm |
c   
c  to compute the angle of rotation start with (16). Now define a coordinate system
c  rotation such that 
c  | x |   | cos theta   -sin theta | | x' |     | x' |
c  |   | = |                        | |    | = T |    |
c  | y |   | sin theta    cos theta | | y' |     | y' |
c  Here x = N, y = E and theta in measuered from N to E
c  insert into (16) and focus on the T E^-1 T^T matrix
c  If theta is chosen correctly, then the off diagonal terms are zero
c  and (16)  becomes
c  
c             T             2
c     (x' - x0) diag(A', B')   (x'-x0) = c2         
c  where theta makes
c     2 E12^-1 cos 2 theta + (E22^-1 - E11^-1) sin 2 thata = 0
c  and
c    A'=E11^-1 cos^2 theta + 2 E12^-1 cos theta sin theta + E22^-1 sin^2 theta
c    B'=E22^-1 cos^2 theta - 2 E12^-1 cos theta sin theta + E11^-1 sin^2 theta
c    The lengths of the semi-axes are c2/sqrt(A') and c2/sqrt(B')
c----
        implicit none
        real sum2,a11,a12,a22,theta,l1,l2 
        logical verbose
        real deter, d11i, d22i, d12i,c,s,a,b
            deter=a11*a22-a12*a12 
            d11i=a22/deter 
            d22i=a11/deter 
            d12i=-a12/deter 
C        WRITE(6,*)d11i,d12i
C        WRITE(6,*)d12i,d22i
            theta=0.5*atan2(2.*d12i,d11i-d22i) 
            c=cos(theta) 
            s=sin(theta) 
            a=c*c*d11i + 2.*c*s*d12i + s*s*d22i 
            b=s*s*d11i - 2.*c*s*d12i + c*c*d22i 
            l1=sqrt(sum2/a) 
            l2=sqrt(sum2/b) 
            theta=theta*57.2957795 
            theta=amod(theta,360.)
            if(theta.lt.0.0)theta = theta + 360.
c-----
c    theta is the angle from north to the x-axis of the errelpe 
c    The length of this axis is givne by X
c-----
            if(verbose)write(6,1) l1,l2,theta 
    1   format(/' ','Error Ellipse  X=',f9.4,' km  Y=',
     1      f7.4,' km  Theta =',f9.4,' deg (azimuth of X axis)'/) 
c-----
c    theta is the angle from north to the x-axis of the errelpe 
c    The length of this axis is givne by X
c-----
        return 
        end 

        subroutine etoh(epoch,date,str,doy)
c-----
c    convert from epoch time to human time
c
c    epoch   - R*8 time in seconds relative to 0000 1 Jan 1970
c    date    - I*4 Julian date
c    str - C*  printable string
c-----
        real*8 epoch
        integer*4 date
        character str*(*)
        integer*4 diy, doy, hour, minute, year, month
        integer*4 day
        real*4 second
        real*8 seclft

        str=' '
        doy = int(epoch/86400.0d+00)
        seclft = dmod(epoch, 86400.0d+00)
        hour = 0
        minute = 0
        second = 0.00
c-----
c    compute hours minutes seconds
c-----
        if(seclft .ne. 0.00d+00)then
c-----
c                before 1970 subtract and add a day
c-----
            if(seclft .lt. 0.0d+00)then
                doy = doy - 1
                seclft = seclft + 86400.00d+00
            endif
            hour = int(seclft/3600.00d+00)
            seclft = dmod(seclft,3600.0d+00)
            minute = int(seclft/60.0d+00)
            second = dmod(seclft,60.0d+00)
        endif

        if(doy .ge. 0)then
            year = 1970
 1000       continue
                diy =  leapdy(year)
                if(doy .lt. diy)go to 2000
                doy = doy - diy
                year = year + 1
            go to 1000
        else
            year = 1969
 1100       continue
                diy =  leapdy(year)
                doy = doy + diy
                if( doy .ge. 0 ) go to 2000
                year = year - 1
            go to 1100
        endif
 2000   continue
        doy = doy + 1
        date = year*1000 + doy
        call mnthdy(year,doy,month,day)
        write(str,110) year,month,day,hour,minute,second
  110   format(i4,i2,i2,i2,i2,f6.3)
c-----
c    guarantee that there are no blanks in the string str
c-----
        do 2100 i=1,17
            if(str(i:i).eq.' ')str(i:i)='0'
 2100   continue
        return
        end

        function leapdy(yr)
        integer*4 yr
        logical t1, t2, t3
        t1 = mod(yr,4).ne.0
        t2 = mod(yr,100).ne.0
        t3 = mod(yr,400).ne.0
        if( .not.t1 .and. t2)then
            isleap = 1
            leapdy = 366
        elseif( .not.t3)then
            isleap = 1
            leapdy = 366
        else
            isleap = 0
            leapdy = 365
        endif
        return
        end

        subroutine mnthdy(year,doy,month,day)
        integer*4 year, doy, month, day
        integer*4 i, dim, leap
        integer*4 dmnth(12)
        data dmnth/31,28,31,30,31,30,31,31,30,31,30,31/
        if(leapdy(year).eq.366)then
            leap = 1
        else
            leap = 0
        endif
        day = doy
        do 100 i=1,12
            month = i
            dim = dmnth(i)
            if(leap.eq.1 .and. i.eq.2)dim = dim + 1
            if(day .le.dim)goto 1000
            day = day - dim 
  100   continue
 1000   continue
        return
        end


        subroutine htoe(year,month,day,hour,minute,second,epoch)
c-----
c    convert calendar date to epoch time since January 1, 1970
c-----
c    year    - I*4   year
c    month   - I*4   month
c    day - I*4   day
c    hour    - I*4   hour
c    minute  - I*4   minute c    second  - I*4   second
c    second  - R*4   seconds
c    epoch   - R*8   time in seconds relative to 00:00 01 Jan 1970
c-----
        integer*4 year, month, day, hour, minute, date, diy
        real*4 second
        real*8 epoch, dtoepo
        integer*4 daymon(12)
        data daymon/0, 31, 59, 90, 120, 151, 181, 212, 243, 273,
     1      304, 334/
        diy = daymon(month) + day
        if(leapdy(year).eq.366 .and. month.gt.2)diy=diy+1
        date = 1000*year + diy
c    write(6,*)'date=',date
        epoch = dtoepo(date) + hour * 3600.0d+00 + 
     1      minute * 60.0d+00 +dble(second)
        return
        end

c-----
c    convert julian date to epoch time
c-----
        function dtoepo(date)
        real*8 dtoepo
        integer*4 date, diy, cnt, days

        cnt = date / 1000
        days = 0
        if (cnt .gt. 1970)then
 1000       continue
            cnt = cnt -1
            if(cnt.lt.1970)go to 2000
                days = days + leapdy(cnt)
            go to 1000
        else if (cnt .lt. 1970)then
 1100       continue
            if(cnt.ge.1970)goto 2000
                days = days - leapdy(cnt)
                cnt = cnt + 1
            go to 1100
        endif
 2000   continue
        diy = (date -1) / 1000
        diy = (date -1 ) -  1000*diy
c    write(6,*)'days=',days,' diy=',diy
        dtoepo = (days + diy) * 86400.0d+00
        return
        end

        function iyesno() 
        logical iyesno
        character*1 ik 
c-----
c    this subroutine gives .true.  for a y(yes) reply 
c                          .false. for a n(no) reply 
c-----
        ikj = +1
  200   continue 
            if(ikj.eq.0)write(6,*)' RESPOND (y/n) '
            read(5,'(a)') ik 
            if(ik.eq.'n'.or.ik.eq.'N')then
                ikj = -1
                iyesno = .false.
            else if(ik.eq.'y'.or.ik.eq.'Y')then
                ikj = +1
                iyesno = .true.
            else
                ikj = 0
            endif
        if(ikj.eq.0)goto 200
        return 
        end 

        function query(str)
        logical query
        character str*(*)
        logical iyesno
        integer LOT
        parameter (LOT=6)
            write(LOT,*)str
            query = iyesno()
        return
        end

        subroutine jbtime(del,caz,saz,time,dtdlat,dtdlon,dtdel,
     1      dtdz,h,phase,anin)
c-----
c    del R*4 - distance in degrees
c    caz R*4 - cos(azimuth)
c    saz R*4 - sin(azimuth)  
c    time    R*4 - travel time in seconds (returned)
c    dtdlat  R*4 - change in travel time for one km NS position change 
c    dtdlon  R*4 - change in travel time for one km EW position change
c    dtdel   R*4 - change in time for one degree in distance
c    h   R*4 - focal depth in km
c    phase   C*8 - Phase identification, e.g., 'P      '
c    anin    R*8 - sin of takeoff angle
c-----
        real*4 del, caz, saz, time, dtdlat, dtdlon, dtdel, h
        character phase*8

      common /tabl/ tp,ts,bofp,bofs,cp,cs,pv,sv
      real*4 tp(14,21),ts(14,24),bofp(17),bofs(20),cp(21),cs(24)
      real*4 pv(14),sv(14)
c-----
c    get velocity at source depth
c-----
      if(h.le.33.) then
        j=1
        vp=pv(j)+(pv(j+1)-pv(j))*h/33.
        vs=sv(j)+(sv(j+1)-sv(j))*h/33.
      else
        do 310 kl=1,12
           ak=kl
           q=33.+ak*63.38
C           if(h-q)315,315,310
           if( h .le. q)then
               go to 315
           endif
  310      continue
  315      j=kl+1
           r=h-q
           r = abs(r)
           vp=pv(j+1)+(pv(j)-pv(j+1))*r/63.38
           vs=sv(j+1)+(sv(j)-sv(j+1))*r/63.38
      endif
c------------------------------------------------------------------
c    calculate the p and s arrival times 
c------------------------------------------------------------------
        call pstime(h,del,0,pt,st)
        if(h.lt.60.0)then
            call pstime(h     ,del,0,ptm,stm)
            call pstime(h+30.0,del,0,ptp,stp)
            dh = 30.0
        else if(h.ge.60.0 .and. h.le.730.0)then
            call pstime(h-63.38,del,0,ptm,stm)
            call pstime(h+63.38,del,0,ptp,stp)
            dh = 2.0*63.38
        else if(h.gt.793.0)then
            call pstime(h     ,del,0,ptp,stp)
            call pstime(h-30.0,del,0,ptm,stm)
            dh = 30.0
        endif
            
        if(phase.eq.'P')then
            time = pt
            dtdz = (ptp - ptm)/(dh)
        else if(phase.eq.'S')then
            time = st
            dtdz = (stp - stm)/(dh)
        endif
c------------------------------------------------------------------
c    The application of the J-B table is limited. If the depth is
c    too large, it is applied to smaller distance.
c------------------------------------------------------------------
c    S does not exist
c-----
        if(h.gt.97..and.del.gt.107.) go to 147 
        if(h.gt.287..and.del.gt.106.) go to 147 
        if(h.gt.540..and.del.gt.105.) go to 147 
        if(h.gt.667..and.del.gt.104.) go to 147 
        if(h.gt.800..and.del.gt.103.) go to 147 
        goto 148
  147   continue
            st = -1.0
  148   continue
c-----
c    P does not exist
c-----

        if(h.gt.33..and.del.gt.105.) go to 145 
        if(h.gt.287..and.del.gt.104.) go to 145 
        if(h.gt.540..and.del.gt.103.) go to 145 
        if(h.gt.800..and.del.gt.102.) go to 145 
        go to 146
  145   continue
            pt = -1.0
  146   continue
      continue
      call psang(del,san,pan,h,vp,vs,dtddp,dtdds)
c-----
c    get information on d2d/dDelta2
c    If positive ray goes up and anin > 0
c    If negative ray goes down and anin < 0
c-----
      call psang(del+1.0,san1,pan1,h,vp,vs,dtddp1,dtdds1)
        if(phase.eq.'P')then
            dtdlat = -dtddp * caz / 111.195
            dtdlon = -dtddp * saz / 111.195
            dtdel = dtddp 
            anin = sin(3.1415927*pan/180.0)
            if(dtddp1 .lt. dtddp)anin = - anin
        else if(phase.eq.'S')then
            dtdlat = -dtdds * caz / 111.195
            dtdlon = -dtdds * saz / 111.195
            dtdel = dtdds 
            anin = sin(3.1415927*san/180.0)
            if(dtdds1 .lt. dtdds)anin = - anin
        endif
        return
      end

      subroutine pstime(h,angi,nfd,pt,st)
c-----------------------------------------------------------------------
c    this subroutine calculate the p & s travel time 
c-----------------------------------------------------------------------
c   Cubic spline interpolation is employed to determine the travel time 
c-----------------------------------------------------------------------
            call bfit(h,angi,nfd,pt,st)
      return
      end

      subroutine psang(dis,san,pan,h,vp,vs,dtddp,dtdds)
c----------------------------------------------------------------------
c    employ bfit subroutine to calculate the 1st derivative
c----------------------------------------------------------------------
      call bfit(h,dis,1,yp,ys)
      degrad = 57.29577951
        dtdds = ys
      u = ys*vs*degrad/(6371.-h)
      u = abs(u)
      if(u.lt.1.) then
         san = asin(u)
      else
         san = 1.5708
      endif
      san=san*57.2958

        dtddp = yp
      c=yp*vp*degrad/(6371.-h)
      c = abs(c)
      if(c.lt.1.) then
         pan = asin(c)
      else
         pan = 1.5708
      endif
      pan=pan*57.2958
      return
      end

      subroutine bfit(h,dis,nfd,yp,ys)
      dimension ut(4),ss(4)
      dimension tp(14,21),ts(14,24),bofp(17),bofs(20),cp(21),cs(24)
      dimension pv(14),sv(14)
      common /tabl/ tp,ts,bofp,bofs,cp,cs,pv,sv
            p1(t)=.25*t**3
            p2(t)=-(1.-t)**2*(1.+t)*.75+1.
            cp1(t)=t**4/16.
            cp2(t)=t*(1./4.+t*(3./8.+t*(1./4.-3./16.*t)))
            dp1(t)=.75*t**2
            dp2(t)=.75 + t*(1.5-2.25*t)
            ddp1(t)=1.5*t
            ddp2(t)=-4.5*t+1.5
c-----------------------------------------------------------------------
c    calculate the coefficients at particular depth
c    if the depth is not exactly at 0.01R, 0.02R ...
c    we use the linear interporation for those coefficients
c-----------------------------------------------------------------------
         if (h .le. 0.) then
              do 11 j=1,21
              cp(j)=tp(1,j)/10.
  11         continue
              do 12 j=1,24
              cs(j)=ts(1,j)/10.
  12         continue
         elseif ( h .gt. 0. .and. h .le. 33.) then
             do 13 j=1,21
             ptemp=(h - 0.)*(tp(2,j)-tp(1,j))/33.
             cp(j) = (tp(1,j) + ptemp)/10.
  13         continue
             do 14 j=1,24
             stemp=(h-0.)*(ts(2,j)-ts(1,j))/33.
             cs(j)=(ts(1,j)+stemp)/10.
  14         continue
         elseif (h.gt.33..and.h.le.797.52) then
            m=0
            r1=63.71
  111       r2=33.+m*r1
            r3=33.+(m+1)*r1
             if (h.gt.r2.and. h.le.r3) then
                do 20 j=1,21
                ptemp=(h-r2)*(tp(3+m,j)-tp(2+m,j))/r1
                cp(j)=(tp(2+m,j)+ptemp)/10.
   20           continue
                do 21 j=1,24
                stemp=(h-r2)*(ts(3+m,j)-ts(2+m,j))/r1
                cs(j)=(ts(2+m,j)+stemp)/10.
   21           continue
             else 
                m=m+1
                if (m .lt. 12) go to 111
             endif
         elseif (h .gt. 797.52) then
              do 25 j=1,21
               cp(j)=(tp(14,j))/10.
   25         continue
              do 26 j=1,24
               cs(j)=(ts(14,j))/10.
   26         continue
         endif
c-----------------------------------------------------------------------
c    calculate the break points
c-----------------------------------------------------------------------
          h1p=(35.0-0.0)/9.
          bofp(1)=0.
          do 30 n=1,9
          bofp(n+1)=bofp(n)+h1p
  30     continue
          h2p=(107.0-35.0)/6.
          bofp(11)=35.0
          do 31 n=11,16
          bofp(n+1)=bofp(n)+h2p
  31    continue
          h1s=(70.0-0.0)/14.
          bofs(1)=0.
          do 32 n=1,14
          bofs(n+1)=bofs(n)+h1s
  32     continue
          h2s=(109.0-70.0)/4.
          bofs(16)=70.0
          do 33 n=16,19
          bofs(n+1)=bofs(n)+h2s
  33    continue
c-----------------------------------------------------------------------
c    determine the distance between particular break points
c-----------------------------------------------------------------------
        if ( dis .le. 35.0 ) then
              jn=1
 100          if( dis .le. bofp(jn+1)) go to 10
               jn=jn+1
               go to 100
 10            u=(dis-bofp(jn))/h1p
        else
              jn=11
 200          if( dis .le. bofp(jn+1)) go to 210
               jn=jn+1
               go to 200
 210           u=(dis-bofp(jn))/h2p
        endif
        if ( dis .le. 70.0 ) then
              jns=1
 300          if( dis .le. bofs(jns+1)) go to 310
               jns=jns+1
               go to 300
 310            us=(dis-bofs(jns))/h1s
        else
              jns=16
 400          if( dis.le.bofs(jns+1)) go to 410
               jns=jns+1
               go to 400
 410           us=(dis-bofs(jns))/h2s
        endif
c-----------------------------------------------------------------------
c  calculate 
c    i. nfd=0 -- interpolation-- travel time (sec)
c   ii. nfd=1 -- 1st derivative --slowness (sec/degree)
c    iii. nfd=2 -- 2nd derivatives --(sec/degree)**2
c-----------------------------------------------------------------------
        v=1.-u
        vv=1.-us
          if(nfd.eq.0) then
             ut(1)=p1(v)
             ut(2)=p2(v)
             ut(3)=p2(u)
             ut(4)=p1(u)
             ss(1)=p1(vv)
             ss(2)=p2(vv)
             ss(3)=p2(us)
             ss(4)=p1(us)
          elseif (nfd.eq.1) then
            if(dis .le. 35.0) then
              ut(1)=-dp1(v)/h1p
              ut(2)=-dp2(v)/h1p
              ut(3)=dp2(u)/h1p
              ut(4)=dp1(u)/h1p
            elseif(dis .gt. 35.0) then
              ut(1)=-dp1(v)/h2p
              ut(2)=-dp2(v)/h2p
              ut(3)=dp2(u)/h2p
              ut(4)=dp1(u)/h2p
            endif
            if(dis .le. 70.0) then
              ss(1)=-dp1(vv)/h1s
              ss(2)=-dp2(vv)/h1s
              ss(3)=dp2(us)/h1s
              ss(4)=dp1(us)/h1s
            elseif(dis .gt. 70.0) then
              ss(1)=-dp1(vv)/h2s
              ss(2)=-dp2(vv)/h2s
              ss(3)=dp2(us)/h2s
              ss(4)=dp1(us)/h2s
            endif
          elseif (nfd.eq.2) then
            if (dis .le. 35.0) then
             ut(1)=ddp1(v)/h1p**2
             ut(2)=ddp2(v)/h1p**2
             ut(3)=ddp2(u)/h1p**2
             ut(4)=ddp1(u)/h1p**2
            else 
             ut(1)=ddp1(v)/h2p**2
             ut(2)=ddp2(v)/h2p**2
             ut(3)=ddp2(u)/h2p**2
             ut(4)=ddp1(u)/h2p**2
            endif
            if (dis .le. 70.0) then
             ss(1)=ddp1(vv)/h1s**2
             ss(2)=ddp2(vv)/h1s**2
             ss(3)=ddp2(us)/h1s**2
             ss(4)=ddp1(us)/h1s**2
            else 
             ss(1)=ddp1(vv)/h2s**2
             ss(2)=ddp2(vv)/h2s**2
             ss(3)=ddp2(us)/h2s**2
             ss(4)=ddp1(us)/h2s**2
            endif
          elseif (nfd.eq.3) then
            if (dis .le. 35.0) then
             ut(1)=-cp1(v)*h1p
             ut(2)=-cp2(v)*h1p
             ut(3)=cp2(u)*h1p
             ut(4)=cp1(u)*h1p
            else 
             ut(1)=-cp1(v)*h2p
             ut(2)=-cp2(v)*h2p
             ut(3)=cp2(u)*h2p
             ut(4)=cp1(u)*h2p
            endif
            if (dis .le. 70.0) then
             ss(1)=-cp1(vv)*h1s
             ss(2)=-cp2(vv)*h1s
             ss(3)=cp2(us)*h1s
             ss(4)=cp1(us)*h1s
            else 
             ss(1)=-cp1(vv)*h2s
             ss(2)=-cp2(vv)*h2s
             ss(3)=cp2(us)*h2s
             ss(4)=cp1(us)*h2s
            endif
          endif       
        if (dis .gt. 35.0) jn=jn+2
        if (dis .gt. 70.0) jns=jns+2
      yp=cp(jn)*ut(1)+cp(jn+1)*ut(2)+cp(jn+2)*ut(3)+cp(jn+3)*ut(4)
      ys=cs(jns)*ss(1)+cs(jns+1)*ss(2)+cs(jns+2)*ss(3)+cs(jns+3)*ss(4)
        return
        end


      blockdata
      common /tabl/ tp,ts,bofp,bofs,cp,cs,pv,sv
      real*4 tp(14,21),ts(14,24),bofp(17),bofs(20),cp(21),cs(24)
      real*4 pv(14),sv(14)
        data ((tp(i,j),j=1,21),i=1,3)/
     1  -617.,   46.,  417.,  781., 1142., 1481., 1816., 2072., 2319.,
     2  2548., 2775., 2991., 2085., 2779., 3442., 4038., 4554., 5002.,
     3  5373., 5732., 6071.,
     4  -246.,   11.,  397.,  756., 1118., 1456., 1788., 2040., 2289.,
     5  2516., 2746., 2950., 2052., 2747., 3410., 4005., 4519., 4968.,
     6  5337., 5698., 6033.,
     7    65.,   10.,  399.,  745., 1100., 1435., 1754., 1999., 2249.,
     8  2474., 2703., 2906., 2017., 2703., 3366., 3959., 4471., 4918.,
     9  5286., 5647., 5978./
        data ((tp(i,j),j=1,21),i=4,6)/
     1   246.,   47.,  400.,  741., 1083., 1418., 1718., 1963., 2209.,
     2  2436., 2661., 2867., 1969., 2664., 3322., 3913., 4425., 4869.,
     3  5237., 5597., 5929.,
     4   342.,  100.,  407.,  741., 1072., 1402., 1681., 1929., 2170.,
     5  2397., 2621., 2833., 1930., 2626., 3281., 3867., 4379., 4821.,
     6  5188., 5550., 5877.,
     7   394.,  160.,  421.,  742., 1064., 1385., 1646., 1898., 2132.,
     8  2362., 2581., 2800., 1894., 2587., 3241., 3823., 4333., 4775.,
     9  5141., 5504., 5829./
        data ((tp(i,j),j=1,21),i=7,9)/
     1   445.,  216.,  443.,  745., 1062., 1360., 1617., 1867., 2098.,
     2  2327., 2544., 2768., 1869., 2549., 3205., 3782., 4290., 4731.,
     3  5096., 5458., 5779.,
     4   478.,  274.,  467.,  753., 1059., 1333., 1591., 1836., 2066.,
     5  2293., 2512., 2730., 1840., 2516., 3170., 3745., 4249., 4687.,
     6  5052., 5413., 5737.,
     7   508.,  329.,  492.,  763., 1045., 1313., 1569., 1808., 2039.,
     8  2265., 2481., 2701., 1806., 2487., 3136., 3708., 4211., 4645.,
     9  5010., 5370., 5699./
        data ((tp(i,j),j=1,21),i=10,12)/
     1   529.,  382.,  520.,  771., 1038., 1299., 1549., 1785., 2015.,
     2  2238., 2455., 2669., 1790., 2459., 3106., 3674., 4175., 4606.,
     3  4971., 5331., 5659.,
     4   569.,  428.,  551.,  780., 1035., 1290., 1532., 1767., 1994.,
     5  2215., 2431., 2644., 1750., 2439., 3077., 3643., 4139., 4573.,
     6  4928., 5310., 5459.,
     7   596.,  474.,  580.,  794., 1035., 1283., 1518., 1752., 1975.,
     8  2195., 2411., 2619., 1732., 2419., 3051., 3613., 4107., 4538.,
     9  4893., 5276., 5419./
        data ((tp(i,j),j=1,21),i=13,14)/
     1   632.,  517.,  613.,  809., 1038., 1277., 1510., 1738., 1960.,
     2  2178., 2392., 2600., 1723., 2399., 3028., 3584., 4076., 4503.,
     3  4858., 5242., 5381.,
     4   662.,  562.,  644.,  826., 1043., 1274., 1503., 1728., 1948.,
     5  2164., 2376., 2582., 1714., 2383., 3005., 3558., 4046., 4470.,
     6  4823., 5209., 5345./
        data ((ts(i,j),j=1,24),i=1,3)/
     1 -1328.,  100.,  909., 1758., 2550., 3344., 3936., 4469., 4990.,
     2  5500., 5988., 6460., 6914., 7353., 7773., 8173., 8555., 7388.,
     3  8182., 8906., 9548.,10094.,10637.,11164.,
     4  -656.,   20.,  886., 1714., 2514., 3297., 3884., 4415., 4936.,
     5  5446., 5934., 6404., 6859., 7297., 7716., 8115., 8501., 7324.,
     6  8126., 8847., 9490.,10034.,10578.,11104.,
     7    -2.,   -6.,  891., 1684., 2487., 3238., 3807., 4344., 4860.,
     8  5368., 5855., 6324., 6777., 7215., 7633., 8032., 8412., 7263.,
     9  8038., 8761., 9399., 9943.,10487.,11023./
        data ((ts(i,j),j=1,24),i=4,6)/
     1   370.,   45.,  893., 1666., 2460., 3181., 3734., 4275., 4787.,
     2  5295., 5778., 6247., 6699., 7135., 7551., 7949., 8324., 7176.,
     3  7956., 8675., 9311., 9855.,10398.,10927.,
     4   630.,  122.,  901., 1655., 2437., 3119., 3667., 4205., 4718.,
     5  5224., 5705., 6174., 6622., 7059., 7472., 7870., 8243., 7101.,
     6  7876., 8593., 9225., 9769.,10312.,10841.,
     7   795.,  218.,  916., 1651., 2414., 3056., 3607., 4137., 4651.,
     8  5155., 5635., 6101., 6549., 6984., 7396., 7793., 8165., 7026.,
     9  7799., 8513., 9141., 9685.,10228.,10755./
        data ((ts(i,j),j=1,24),i=7,9)/
     1   877.,  327.,  935., 1655., 2386., 2997., 3550., 4073., 4589.,
     2  5088., 5567., 6032., 6479., 6911., 7323., 7717., 8091., 6954.,
     3  7723., 8436., 9060., 9603.,10147.,10673.,
     4   928.,  437.,  960., 1665., 2349., 2947., 3492., 4014., 4528.,
     5  5025., 5503., 5964., 6411., 6840., 7251., 7643., 8016., 6888.,
     6  7649., 8360., 8984., 9524.,10068.,10593.,
     7  1001.,  535., 1000., 1662., 2314., 2901., 3441., 3961., 4475.,
     8  4969., 5443., 5904., 6348., 6776., 7184., 7575., 7943., 6820.,
     9  7580., 8289., 8908., 9449., 9993.,10517./
        data ((ts(i,j),j=1,24),i=10,12)/
     1  1055.,  631., 1039., 1662., 2287., 2862., 3395., 3914., 4427.,
     2  4916., 5389., 5847., 6290., 6716., 7122., 7512., 7877., 6765.,
     3  7515., 8222., 8837., 9377., 9921.,10446.,
     4  1097.,  725., 1080., 1669., 2265., 2829., 3356., 3874., 4382.,
     5  4870., 5340., 5797., 6237., 6661., 7066., 7451., 7817., 6673.,
     6  7463., 8155., 8776., 9297., 9891.,10019.,
     7  1131.,  815., 1124., 1675., 2250., 2800., 3324., 3841., 4344.,
     8  4828., 5296., 5751., 6187., 6610., 7011., 7394., 7763., 6619.,
     9  7408., 8095., 8712., 9232., 9827., 9951./
        data ((ts(i,j),j=1,24),i=13,14)/
     1  1183.,  897., 1173., 1685., 2240., 2778., 3299., 3813., 4310.,
     2  4792., 5257., 5708., 6142., 6560., 6960., 7342., 7702., 6562.,
     3  7357., 8035., 8656., 9158., 9805., 9394.,
     4  1233.,  976., 1224., 1698., 2236., 2761., 3280., 3789., 4282.,
     5  4759., 5221., 5669., 6101., 6515., 6913., 7291., 7649., 1399.,
     6   890., 1399., 2186., 2943., 3692., 4402./
      data (pv(i),sv(i),i=1,14)/
     1   6.75,  3.90,  7.75,  4.35,  7.94,  4.44,  8.13,  4.54,
     2   8.33,  4.64,  8.54,  4.74,  8.75,  4.85,  8.97,  4.96,
     3   9.50,  5.23,  9.91,  5.46, 10.26,  5.67, 10.55,  5.85,
     4  10.77,  6.00, 10.99,  6.12/
      end


        subroutine locate(t0,tlat,tlon,nr,depth,tratbl,
     1      dlat,dlon,dz,dt,l1,l2,theta,fixed,
     2      verbose,veryverbose,fmask,mxiter,rms,inter)
c-----
c    subroutine arguments
c-----
        real*8 t0
        real*4 tlat, tlon, depth, dlat, dlon, dz,dt, l1, l2, theta
        integer nr, fmask(4), mxiter
        character tratbl*8
        logical fixed, verbose, inter, veryverbose
c-----
        parameter (NLAYR=100)
        common/event/iyr,imon,idy,ihr,imin,cmnt
        common/invert/sigma(4),c(4),co(4,4),a(4,4),b(4,4),y(4),sdobs,
     1      nass,ndef
        parameter (NPHA=3000)
        common/phases/sta(NPHA),ipwt(NPHA),tp(NPHA),xlat(NPHA),
     1  xlon(NPHA),pfm(NPHA),pph(NPHA),pqual(NPHA),pres(NPHA),
     2  pchan(NPHA),pwt(NPHA),ips(NPHA,2),jain(NPHA),jaz(NPHA),
     3  jbaz(NPHA),asid(NPHA),arid(NPHA),mapd(NPHA),elev(NPHA), 
     4  ondate(NPHA),offdat(NPHA)
        integer*4 asid, arid, mapd
        character*6 sta
        integer*4 ipwt
        real*8 tp
        real*4 xlat,xlon
        character*2 pfm
        character*8 pph
        character*1 pqual
        real *4 pres,pwt
        character*2 pchan
        integer*4 ips
        real*4 jain,jaz,jbaz
        real*4 elev
        integer*4 ondate, offdat
c-----
c    VARIABLES IN COMMON BLOCK phases
c        sta - station name
c        ipwt    - integer phase weight 0-4
c        tp  - phase arrival time in seconds past midnight
c        xlat    - latitude  degrees north
c        xlon    - longitude degrees east
c        pfm - phase first motion x,X,c,C,d,D,+, or -
c        pph - phase classification, P PKP Lg
c        pqual   - phase arrival quality e,E,i,I
c        pres    - arrival time residual
c        pchan   - phase channel on station
c        pwt - final phase weight
c        ips(i,j)- (i,1) = table 1 for i th phase
c            - (i,2) = table 2 to compute phase(ips(i,1) - ips(i,2)
c        jain    - Phase takeoff angle at source in degrees
c        jaz - ray aziumth from source to receiver
c        elev    - elevation in kilometers
c        ondate  - 
c        offdat  -
c-----
        integer*4 date
        character str*20
        character*1 quald,quals
        common/a6/nl,nthpha
        real*4 x(4,NPHA)
        common/a10/anin(NPHA)
        common/a2/delta(NPHA),dx(NPHA),dy(NPHA),tt(NPHA)
        real*4 tt
        dimension temp1(NPHA),key(NPHA),temp(NPHA)
c-----
c    internal
c-----
        real*8 dsum
        integer*4 doy
c-----
c    set up trial solution slightly away from nearest
c    station
c-----
        if(tratbl.eq.'TELE')then
            del = 1.0
        else if(tratbl.eq.'BEAM')then
            del = 1.0
        else
            del = 0.1
        endif
        if(fixed)then
            if(fmask(1).eq.0 ) then
                tlat = tlat + del
            endif
            if(fmask(2).eq.0 ) then
                tlon = tlon + del
            endif
        endif
c-----
c    damping for least squares matrix 
c-----
        sigma(1)=0.0001 
        sigma(2)=0.0001 
        sigma(3)=0.0001 
        sigma(4)=0.0001 
        idamp = 0
        oldrms = 99999.
c-----
c    get depth
c-----
        if(fmask(3).eq.0 .or. depth.le. -1000.0)then
            if(verbose )then
            write(6,*)' enter depth, depth < 0 is fixed at abs(depth)'
            endif
            if(inter)then
                read(5,*)z
            else
                z = depth
            endif
        else
            z = depth
        endif
c-----
c    safety check for non-negative depth
c-----
        if(z.lt.0.0)then
            idphld = 1
            z = -z
        else
            idphld=0 
        endif
        iter = 0 
   98   continue 
c-----
c    iteration variable damping for the depth element, to
c    improve the epicenter location before determining
c    the depth
c-----
        if(iter.le.2)then
            sigma(3)=1.0
        else if(iter.gt.2 .and. iter.le.3)then
            sigma(3) = 0.10
        elseif(iter.gt.3.and.iter.le.6)then
            sigma(3)=0.01
        else
            sigma(3)=0.0001
        endif
        do 200 k=1,nr 
            i = mapd(k)
            call trtim(tlat,tlon,z,xlat(i),xlon(i),anin(i),
     1          delta(i),jaz(i),jbaz(i),ips(i,1),ips(i,2),tratbl,
     2          y,tt(i),pph(i),elev(i),deldeg)
            temp1(i) = jaz(i)
            temp(i) = delta(i)
            do 201 j=1,4
                x(j,i) = y(j)
  201       continue
  200   continue 
c-----
c    get azimuthal gaps for azimuthal weighting
c-----
        call sort(temp1,key,nr)
        gap=temp1(1)+360.-temp1(nr)
        do 220 i=2,nr
            dtemp=temp1(i)-temp1(i-1)
            if(dtemp.gt.gap)gap=dtemp
  220   continue
        igap=int(gap)
        gap= float(igap)
c-----
c    get minimum distance
c-----
        call sort(temp,key,nr)
        dmin=temp(1)
c-----
c
c-----
        avwt=0.0
        sum2=0.0
c-----
c    for stability the weighted sum of residuals should
c    equal zero, if not, adjust the origin time estimate
c
c    also set up sort on distance
c-----
        do 250 k=1,nr
            i = mapd(k)
            temp(k) = delta(i)
            call getres(tres,tp(i),tt(i),t0,ips(i,1),ips(i,2))
            if(ips(i,2).eq.0)then
                wtt = wt(ipwt(i),delta(i),tres,iter,0,tratbl)
                avwt=avwt + wtt
                sum2=sum2 + wtt*tres
            endif
  250   continue
        if(avwt.gt.0.0)then
            sum2=sum2/avwt
            t0=t0 + sum2
        endif
c-----
c    sort in order of increasing distance
c-----
        call sort(temp,key,nr)
        avwt=0.0
        sum2=0.0
        do 298 i=1,4
            y(i)=0.0
            do 299 j=1,4
                b(i,j)=0.0
  299       continue
  298   continue
        do 300 kk=1,nr
            i = mapd(kk)
            if(tt(i).ge.0.0)then
            call getres(tres,tp(i),tt(i),t0,ips(i,1),ips(i,2))
            wtt=wt(ipwt(i),delta(i),tres,iter,0,tratbl)
            pwt(i) = wtt
            avwt=avwt+wtt
            sum2=sum2+tres*tres*wtt
            do 301 j=1,4
                y(j)=y(j)+x(j,i)*tres*wtt
  301       continue
            do 400 j=1,4
                do 401 k=1,4
                    b(j,k)=b(j,k)+x(j,i)*x(k,i)*wtt
  401           continue
  400       continue
            endif
  300   continue
        avwt=avwt/nr
        sum2=sum2/avwt
        do 302 i=1,4
            y(i)=y(i)/avwt
            do 303 j=1,4
                b(i,j)=b(i,j)/avwt
  303       continue
  302   continue
        call etoh(t0,date,str,doy)
        if(verbose)write(6,6)tlat,tlon,z,str,sum2
c-----
c    check if last change was too extreme and caused
c    sum of errors to increase
c-----
c-----
        ijump = 0
        if(sum2.gt.1.2*oldrms .and. iter.gt.3.and. idamp.lt.3)then
            t0 = t0 -c(4)/4.
            z = z - c(3)/4.
            tlat = tlat - c(2)/4.
            tlon = tlon - c(1)/4.
            ijump = 1
            idamp=idamp +1
        endif
        if(idamp.eq.4)then
            ijump = 0
            fmask(3) = 1
        endif
        if(ijump.eq.1)goto 98
        oldrms = sum2
        idamp = 0
c-----
c    form a^T a + sigma(i) 
c-----
        do 600 i=1,4 
            do 601 j=1,4 
                a(i,j)=b(i,j) 
  601       continue
            a(i,i)=b(i,i)+sigma(i) 
  600   continue
c-----
c    modify matrices for fixed source parameters
c-----
        do 1600 j=1,4
            if(fmask(j).eq.1 .and. fixed)then
                call fixij(a, y,j)
            endif
 1600   continue
        if(idphld.eq.1)call fixij(a,y,3)
c-----
c    invert
c-----
        call matinv(a,4) 
c-----
c    find solution vector
c-----
        do 800 i=1,4 
            dsum=0.0d+00
            do 801 j=1,4 
                dsum=dsum+a(i,j)*y(j) 
  801       continue
            c(i)=real(dsum )
  800   continue 
c-----
c    print solution vector 
c-----
c    determine km per degree latitude and longitude
c-----
        call delaz(tlat-0.5,tlon+0.5,tlat+0.5,tlon-0.5,del,az,baz,
     1      delkm,saz,caz)
        dxx = delkm * saz
        dyy = delkm * caz
        dxx=abs(dxx) 
        dyy=abs(dyy) 
        if(tratbl.eq.'TELE')then
            if(abs(c(1)) .gt. 1000.0 .or. abs(c(1)).gt. 1000.0)then
                c(1) = sign(1000.0,c(1))
                c(2) = sign(1000.0,c(2))
                c(3) = 0.0
                c(4) = sign(50.0,c(4))
            endif
        else if(tratbl.eq.'BEAM')then
            if(abs(c(1)) .gt. 1000.0 .or. abs(c(1)).gt. 1000.0)then
                c(1) = sign(1000.0,c(1))
                c(2) = sign(1000.0,c(2))
                c(3) = 0.0
                c(4) = sign(50.0,c(4))
            endif
        endif
c-----
c    now convert change in EW and NS km  to change in EW and NS degrees
c-----
        c(1) = c(1)/dxx 
        c(2) = c(2)/dyy 
        tlat = tlat + (c(2) )
        tlon = tlon + (c(1) )
c-----
c    enforce a largest change in depth
c    use a relative measure together with absolute
c----
        cz = abs(c(3))
        if(cz.gt.5.0 .and. z.le.40.0)then
            c(3) = sign(5.0,c(3))
        else if(cz .gt. 0.25*z .and. z.gt.40.0)then
            c(3) = sign(0.25*z, c(3))
        endif
c    if(abs(c(3)).gt.5.) c(3)=sign(5.0,c(3)) 
        z=z+c(3) 
        if(z.le.0.0) z = 0.0001 
        t0 = t0 + c(4) 
c-----
c    test for convergence
c
c    convergence criteria -- terminate is latest
c        correction on lat and lon is less than 0.0001 deg,
c      on depth is less than 0.01 km , on origin time
c      is less than 0.01 sec and if at least 5 iterations
c      have been performed
c-----
        if(abs(c(1)).lt.1.e-4.and.abs(c(2)).lt.1.e-4.and. 
     1      abs(c(3)).lt.1.e-2.and.abs(c(4)).lt.1.e-2.and.iter.gt.10) 
     2           goto 920 
      iter=iter + 1
c-----
c    one needa a go to somewhere
c-----
      if(iter.lt.mxiter) go to 98 
  920 continue 
c-----
c    obtain true depth variable inverse
c-----
        if(fixed .or. idphld.eq.1)then
            do 940 i=1,4 
                do 941 j=1,4 
                    a(i,j)=b(i,j) 
  941           continue
                a(i,i)=b(i,i) + sigma(i) 
  940       continue
            call matinv(a,4) 
        endif
c-----
c    determine covariance matrix
c    Crosson, Robert S. (1976). Crstal structure modeling of
c      earthquake data 1. Simultaneous least squares estimation of hypocenter
c      and velocity parameters, J. Geophys. Res 81, 3036-3046
c            T     2  -1  T     T     2  -1
c      C = (A A + d I)   A A  (A A + d I)     d is damping (18) and (19)
c                        T     2  -1        T      T
c    below  a(i,j) is  (A A + d I)     and A A is X X
c-----
    6   format(1x ,2f10.4,f10.2,1x,a20,f10.2) 
        do 949 i=1,4 
            do 950 m=1,4 
                sum=0.0 
                do 951 kk = 1,nr
                    l = mapd(kk)
                    if(tt(l).ge.0.0)then
                
                    call getres(tres,tp(l),tt(l),t0,
     1                  ips(l,1),ips(l,2))
                    wtt=wt(ipwt(l),delta(l),tres,iter,0,tratbl)
                    do 952 j=1,4 
                        do 953 k=1,4 
                            sum = sum + 
     1                  a(i,j)*x(j,l)*x(k,l)*a(m,k)*wtt/avwt 
  953                   continue
  952               continue
                    continue
                    endif
  951           continue 
                co(i,m)=sum 
  950       continue 
  949   continue 
        if(verbose)then
            write(6,*)nr,' phases used '
        endif
        nass = nr
        sum =0.0 
        if(nr.gt.4)then
            ndef = nr -4
            if(fmask(1).eq.1)ndef = ndef + 1
            if(fmask(2).eq.1)ndef = ndef + 1
            if(fmask(3).eq.1)ndef = ndef + 1
            if(fmask(4).eq.1)ndef = ndef + 1
            sum3=sum2/ndef
        else
            sum3=sum2/nr
            ndef = nr
        endif
        dln = sqrt(sum3*abs(co(1,1))) 
        dlt = sqrt(sum3*abs(co(2,2))) 
        dlon = dln 
        dlat = dlt 
c-----
c-----
c    dlat and dlon are the location errors in degrees
c    dlt  and dln  are the location errors in km
c
c    Teleseism travel times are with respect to degrees
c    Local models are a function of km. thus the
c    covariance matrix has different units for each
c-----
c    if(tratbl.eq.'TELE')then
c        dlat = dlt
c        dlon = dln
c        dlt = dlt*dyy
c        dln = dln*dxx 
c    else if(tratbl.eq.'BEAM')then
c        dlat = dlt
c        dlon = dln
c        dlt = dlt*dyy
c        dln = dln*dxx 
c    else
            dlat = dlt/dyy 
            dlon = dln/dxx 
c    endif
        dz=sqrt(sum3*abs(co(3,3)) )
        dt=sqrt(sum3*abs(co(4,4)) )
        sum3sr=sqrt(sum3) 
        sdobs = sum3sr
        erh=sqrt(dlt**2 + dln**2)
        call qual(quals,quald,sum3sr,erh,dz,nr,gap,dmin,z)
        if(verbose)then
            write(6,*) ' ' 
            if(tratbl.eq.'TELE')then
                write(6,10)
            else if(tratbl.eq.'BEAM')then
                write(6,10)
            else
                write(6,11) 
            endif
        endif
        sumf=0.0
        numf=0
        suma=0.0
        numa=0
        sumf2 = 0.0
        suma2 = 0.0
c-----
c    readjust origin time so that sum of residuals
c    is zero
c-----
        avwt = 0.0
        sum2 = 0.0
        do 980 k=1,nr
            i = mapd(k)
            if(tt(i).gt.0.0)then
            call getres(tres,tp(i),tt(i),t0,ips(i,1),ips(i,2))
            wtt=wt(ipwt(i),delta(i),tres,iter,0,tratbl)
            avwt = avwt + wtt
            sum2 = sum2 + wtt*tres
            endif
  980   continue
        sum2 = sum2/avwt
        t0 = t0 + sum2
        do 990 k=1,nr 
            l=key(k)
            i=mapd(l)
            if(mod(k,22).eq.0 .and. verbose .and. inter)then
                write(6,*)'RETURN FOR NEXT SCREEN'
                read(5,'(1x)')
            endif
            dtp = 0.0 
            call getres(dtp,tp(i),tt(i),t0,ips(i,1),ips(i,2))
            iain=int(asin(anin(i))*180./3.1415927  )
            if(iain.lt.0) iain=iain+180 
            iain=180-iain 
            jain(i)=iain
            pres(i)=dtp
            if(tt(i).lt.0.0)pres(i)= -999.0
            if(ips(i,2).eq.0 .and. verbose )then
                call etoh(tp(i),date,str,doy)
                write(6,12)sta(i),pchan(i),delta(i),jaz(i),
     1              jain(i),str  ,pres(i),ipwt(i),
     2              pqual(i),pfm(i),pph(i),pwt(i)
            else if(ips(i,2).gt.0 .and. verbose)then
                write(6,22)sta(i),pchan(i),delta(i),jaz(i),
     1              jain(i),tp(i),pres(i),ipwt(i),
     2              pqual(i),pfm(i),pph(i),pwt(i)
            endif
  990   continue 
            if(verbose .and. inter)then
                write(6,*)'RETURN FOR NEXT SCREEN'
                read(5,'(1x)')
            endif
   10 format(1x,'STA  COMP DIS(D) AZM  AIN    ARR TIME',
     1  10x,'   RES(SEC) WT QFM PHASE    WGT')
   11 format(1x,'STA  COMP DIS(K) AZM  AIN    ARR TIME',
     1  10x,'   RES(SEC) WT QFM PHASE    WGT')
   12 format(1x,a6,1x,a2,f7.2,2f5.0,1x,a20  ,1x,f10.2,2x,i1,1x, 
     1  a1,a2,a8,f5.2)
   22 format(1x,a6,1x,a2,f7.2,2f5.0,1x,f20.3,1x,f10.2,2x,i1,1x, 
     1  a1,a2,a8,f5.2)
        if(abs(pres(i)).gt.9.99)pres(i)=sign(9.99,pres(i))
        call etoh(t0,date,str,doy)
c-----
c    invoke errelp once more to get printed outpur
c-----
        if(veryverbose)then
           WRITE(6,*)'NDEF  :',ndef
           WRITE(6,*)'STDERR:',sqrt(sum3)
           WRITE(6,*)'VARIANCE-COVARIANCE :'
           WRITE(6,*)(sum3*co(1,j),j=1,4)
           WRITE(6,*)(sum3*co(2,j),j=1,4)
           WRITE(6,*)(sum3*co(3,j),j=1,4)
           WRITE(6,*)(sum3*co(4,j),j=1,4)
        endif
c-----
c       we call this in Flinn order, e.g., lat, lon, depth, time
        call  errelp(sum3,co(2,2),co(1,2),co(1,1),theta,l1,l2,
     1      verbose)
        if(verbose)then
            write(6,7)sqrt(sum3),tratbl,tlat,dlat,dlt,tlon,dlon,dln,
     1          z,dz,t0,dt,str(1:18),dt,
     2          str(1:4),str(5:6),str(7:8),str(9:10),str(11:12),
     3          str(13:14),str(16:18),quals, quald, igap 
            write(6,8)tlat,tlon,
     1          z,
     2          str(1:4),str(5:6),str(7:8),str(9:10),str(11:12),
     3          str(13:14),str(16:18)
            write(6,9)
     2          str(1:4),doy,str(9:10),str(11:12),
     3          str(13:14),str(16:18)
        endif
    7 format(
     1  ' RMS Error        :',10x,f10.3,13x,' sec'/
     1  ' Travel_Time_Table:',10x,a/
     1  ' Latitude         :',10x,f10.4,' +-',f10.4,' N',5x,f10.4,' km'/
     2  ' Longitude        :',10x,f10.4,' +-',f10.4,' E',5x,f10.4,' km'/
     3  ' Depth            :',10x,f10.2,' +-',f10.2,' km'/
     4  ' Epoch Time       :', 5x,f15.3,' +-',f10.2,' sec'/
     6  ' Event Time       :', 2x,a18  ,' +-',f10.2,' sec'/
     6  ' Event (OCAL)     :',2x,a4,1x,a2,1x,a2,1x,a2,1x,a2,1x,a2,1x,a3/
     5  ' HYPO71 Quality   :',18x,a1,a1/
     7  ' Gap              :',16x, i4  , 3x  ,10x  ,' deg'/
     8  ' ') 
c-----
c      added 17 FEB 2010 for use with gsac and/or sac/sac2000
c-----
    8 format(
     1  '    ch evla ',f10.4,' evlo ',f10.4, ' evdp ', f10.2/
     2  '    ch ocal ',2x,a4,1x,a2,1x,a2,1x,a2,1x,a2,1x,a2,1x,a3/
     3  ' ') 
    9 format(
     2  '    ch o gmt ',2x,a4,1x,i3,1x,1x,a2,1x,a2,1x,a2,1x,a3/
     2  ' ') 
        depth = z
        rms = sum3
        return
        end

        function sd(n,x,x2)
c-----
c    standard deviation, n = number of observations
c    n   integer number of observations
c    x   real*4  mean
c    x2  real*4  sum of x(i)**2
c    sd  real*4  standard error
c-----
        sd = sqrt((x2-n*x*x)/float(n-1))
        return
        end

        subroutine getres(tres,tp,tt,t0,ips1,ips2)
c-----
c    determine travel time residual according to type
c    of phase
c-----
c    tres    - residual
c    tp  - phase arrival time
c    tt  - phase travel time
c    t0  - origin time
c    ips - phase id for travel time table
c------
        real*4 tres, tt
        real*8 tp, t0
        integer ips1, ips2
            if(ips2.eq.0)then
                tres=real(tp-tt-t0)
            elseif(ips2.gt.0 )then
                tres=real(tp-tt)
            endif
        return
        end

        subroutine matinv(a,n) 
        real*4 a(4,4),v(4) 
c-----
c    this algorithm determines the inverse of a symmetric n x n matrix 
c    the inverse is determined and returned in the same location as a 
c-----
        nm1 = n - 1 
        do 35 k = 1,n 
            pivot = 1.0/a(1,1) 
            do 20 i = 2,n 
                v(i-1) = a(1,i) 
   20       continue
            do 30 i = 1,nm1 
                yy = -v(i) * pivot 
                a(i,n) = yy 
                do 31 j = 1,nm1 
                    a(i,j) = a(i+1,j+1) + v(j) * yy 
   31           continue
   30       continue
            a(n,n) = - pivot 
   35   continue
        do 40 i = 1,n 
            do 41 j = 1,n 
                a(i,j) = - a(i,j) 
   41       continue
   40   continue
        do 50 i = 2,n 
            ii = i - 1 
            do 51 j = 1,ii 
                a(i,j) = a(j,i) 
   51       continue
   50   continue
        return 
        end 

        subroutine modlin(tratbl,model,verbose,inter)
c-----
c    subroutine to input earth model
c        tratbl  -> string to identify model
c-----
        character tratbl*8
        integer*4 model
        logical verbose, inter
c-----
        parameter (NPHA=3000)
            parameter (NLAYR=100,NMDPHA=6)
            common/a21/phathe(NMDPHA)
            character phathe*8
        common/a6/nl,nthpha

        logical ext
c-----
c    use model defined in command line, but only if it is within
c    the range of possible models
c-----
c-----
c    read user specified crustal model
c-----
        inquire(file='VEL.MOD',exist=ext)
        if(.not. ext)then
            call usage('velocity model file VEL.MOD missing')
        endif
        open(2,file='VEL.MOD',status='unknown',access='sequential', 
     1  form='formatted') 
        rewind 2 
c-----
c    The file VEL.MOD contains all models for use. There is no limit
c    to the number of possible models, however the fewer there are,
c    the simpler the search. The structure of the model file is as
c    follows:
c
c    MODEL NAME      character*4
c    NUMBER OF INTERFACES nl
c    P-velocity     Depth to top of layer
c        ( a total of nl of these)
c        ........................
c    MODEL NAME
c    NUMBER_OF_INTERFACES NUMBER_OF_PHASES
c    'Phase 1' 'Phase 2' 
c    Depth_to_top_of_layer   Velocity_Phase_1    Velocity_Phase_2
c        ( more layers)
c        ........................
c    sequence repeated until an END OF FILE
c-----
        nmodel = 0
        ione = 1
        itwo = 2
        nthpha = 2
        phathe(1) = 'P       '
        phathe(2) = 'S       '
        if(verbose .and.inter)then
            write(6,*)' CHOOSE VELOCITY MODEL'
            write(6,105)'TELE    ',ione, (phathe(i),i=1,nthpha)
            write(6,105)'BEAM    ',itwo, (phathe(i),i=1,1)
        endif
 1000   continue
            call getmod(nnl,nthpha,tratbl)
            if(nnl.le.0)goto 1001
            nmodel = nmodel + 1
            if(verbose .and.inter)then
                write(6,105)tratbl,nmodel+2,(phathe(i),i=1,nthpha)
            endif
  105   format(' ',a8,' = Model(',i2,') ',6a)
        goto 1000
 1001   continue
        if(verbose .and.inter)then
            write(6,*)' WHICH MODEL ? 1 - ',nmodel+2
        endif
        if(model.le.0 .and. inter)then
            read(5,*)jmodel
        else
            jmodel = model
        endif
        if(jmodel.le.0.or.jmodel.gt.(nmodel+2))goto 1001
c-----
c    process the appropriate model, checking for teleseism
c-----
        if(jmodel.eq.1)then
            tratbl = 'TELE    '
            nthpha = 2
            phathe(1) = 'P       '
            phathe(2) = 'S       '
        else if(jmodel.eq.2)then
            tratbl = 'BEAM    '
            nthpha = 1
            phathe(1) = 'P       '
        else
            call inimod(jmodel-2,tratbl)
        endif
        close(2)
        return
        end

        subroutine inimod(jmodel,tratbl)
        character tratbl*8
        parameter (NPHA=3000)
c--------this subroutine establishes the earth model to be used--------
            parameter (NLAYR=100,NMDPHA=6)
            common/a20/v(NLAYR,NMDPHA),d(NLAYR),
     1          vsq(NLAYR,NMDPHA),thk(NLAYR),
     2          tid(NLAYR,NLAYR,NMDPHA),did(NLAYR,NLAYR,NMDPHA)
            common/a21/phathe(NMDPHA)
            character phathe*8
      common/a6/nl,nthpha
      common/a22/f(NLAYR,NLAYR) 
        rewind 2
        do 1002 i=1,jmodel
            call getmod(nl,nthpha,tratbl)
 1002   continue
        do 1004 j=1,nthpha
        do 1003 i=1,nl
            vsq(i,j)   = v(i,j)*v(i,j)
 1003   continue
 1004   continue
        close(2)
        n1=nl-1 
c-----layer thickness, f and g matrices--------------------------------
        do 145 l=1,n1 
            thk(l)=d(l+1)-d(l) 
  145   continue
c-----
c    the following matrices are computed once and saved.
c    this is useful when more than one event is located,
c        with the same earth model
c
c    f(j,l)  Matrix to indicate number of ray paths in each layer
c        It looks like
c
c        ............ l
c        .2 2 2 2 2 2
c        .1 2 2 2 2 2
c        .1 1 2 2 2 2
c        .1 1 1 2 2 2
c        .1 1 1 1 2 2
c       j.1 1 1 1 1 2
c
c        j is the source layer index, so that, for example
c        if the source is in layer 3, there is only one
c        upgoing ray in layers 1,2 and both up- and down-going
c        rays in the layers 3 - bottom
c
c    tid(j,m) Intercept time for a refraction with velocity v(m)
c        for a source in layer j. 
c        Enhancement RBH - if tid(j,m) = negative, no refraction
c        possible, this would arise if there is a low velocity
c        layer beneath
c
c    did(j,m) Distance at which a refraction arrives, for source in
c        layer j and a refraction velocity v(m)
c        Enhancement RBH - if did(j,m) = negative, no refraction
c        possible, e.g., low velocity layer beneath
c-----
c-----compute tid and did----------------------------------------------
        do 150 j=1,nl 
            do 151 l=1,nl 
                if(l.lt.j)then
                    f(l,j) = 1.
                else
                    f(l,j) = 2.
                endif
  151       continue 
  150   continue 
c-----
c    initialize
c-----
        do 173 jj=1,nthpha
         do 165 j=1,nl 
            do 166 m=1,nl
                tid(j,m,jj) = 0.0
                did(j,m,jj) = 0.0
  166       continue
  165   continue
c-----
c    loop over source layer index j
c-----
        do 170 j=1,nl
c-----
c    loop over refractor index
c-----
            do 171 m=j+1,nl
c-----
c        get maximum velocity of layers above refractor
c-----
                vmax = 0.0
                do 172 k=1,m-1
                    if(v(k,jj).gt.vmax)vmax=v(k,jj)
  172           continue
c-----
c    test for valid refraction
c    which requires refractor velocity to be greater than
c    velocity of all layers above refractor
c-----
c        if(v(m).le.vmax .or. v(m).le. v(j))then
            if(v(m,jj) .le. vmax)then
                tid(j,m,jj) = -1.0
                did(j,m,jj) = -1.0
            else
                m1 = m-1
                do 160 l=1,m1
                   sqt=sqrt(abs(v(m,jj)-v(l,jj))*(v(m,jj)+v(l,jj)) )
                   tim=thk(l)*sqt/(v(l,jj)*v(m,jj)) 
                   dim=thk(l)*v(l,jj)/sqt 
                   tid(j,m,jj)=tid(j,m,jj)+f(l,j)*tim 
                   did(j,m,jj)=did(j,m,jj)+f(l,j)*dim 
  160           continue
            endif
  171       continue
  170   continue 
  173   continue
      return 
      end 

        subroutine getmod(nl,nthpha,tratbl)
            character tratbl*8
            parameter (NLAYR=100,NMDPHA=6)
            common/a20/v(NLAYR,NMDPHA),d(NLAYR),
     1          vsq(NLAYR,NMDPHA),thk(NLAYR),
     2          tid(NLAYR,NLAYR,NMDPHA),did(NLAYR,NLAYR,NMDPHA)
            common/a21/phathe(NMDPHA)
            character phathe*8
c-----
c    READ VEL.MOD FILE
c-----
        nl = -1
        read(2,'(a)',end=99,err=99)tratbl
        read(2,*,end=99,err=99)nl,nthpha
        read(2,*)(phathe(j),j=1,nthpha)
        do 100 i=1,nl
            read(2,*)d(i),(v(i,j),j=1,nthpha)
  100   continue
   99   continue
        return
        end


        subroutine fixij(a,y,j)
c-----
c    fix inversion parameters by forcing changes to be zero
c-----
c    a(4,4)  R*4 - symmetric least squares matrix
c    y(4)    R*4 - vector
c    j   I   - row
c    i   I   - column
c-----
        real*4 a(4,4), y(4)
            do 699 i=1,4 
                a(j,i)=0.0 
                a(i,j)=0.0 
  699       continue
            a(j,j)=1.0 
            y(j)=0.0 
        return
        end

        subroutine qual(quals,quald,rms,erh,erz,nr,gap,dmin,depth)
        character*1 quals,quald
            igap = int(gap)
            quals='D'
            if(rms.lt.0.5.and.erh.le.5.0)quals='C'
            if(rms.lt.0.3.and.erh.le.2.5.and.erz.le.5.)quals='B'
            if(rms.lt.0.15.and.erh.le.1.0.and.erz.le.2.)quals='A'
            quald='D'
            if(nr.ge.6.and.igap.le.180.0.and.dmin.le.50.)quald='C'
            d1=depth
            if(d1.lt.5.0)d1=5.0
            d2=2.*d1
            if(nr.ge.6.and.igap.le.135.0.and.dmin.le.d2)quald='B'
            if(nr.ge.6.and.igap.le.90.0.and.dmin.le.d1)quald='A'
        return
        end
    
        subroutine setinl(nr,t0,tlat,tlon,tratbl,
     1      verbose,fmask,inter)
c-----
c    subroutine arguments
c-----
        integer nr, fmask(4)
        real*8 t0
        real*4 tlat, tlon
        character tratbl*8
        logical verbose, inter
c-----
c    ESTIMATE TRIAL COORDINATES OF SULUTION
c    USING TIME AND STATION NEAREST THE SOURCE
c-----
        parameter (NPHA=3000)
        common/phases/sta(NPHA),ipwt(NPHA),tp(NPHA),xlat(NPHA),
     1  xlon(NPHA),pfm(NPHA),pph(NPHA),pqual(NPHA),pres(NPHA),
     2  pchan(NPHA),pwt(NPHA),ips(NPHA,2),jain(NPHA),jaz(NPHA),
     3  jbaz(NPHA),asid(NPHA),arid(NPHA),mapd(NPHA),elev(NPHA), 
     4  ondate(NPHA),offdat(NPHA)
        integer*4 asid, arid, mapd
        character*6 sta
        integer*4 ipwt
        real*8 tp
        real*4 xlat,xlon
        character*2 pfm
        character*8 pph
        character*1 pqual
        real *4 pres,pwt
        character*2 pchan
        integer*4 ips
        real*4 jain,jaz,jbaz
        real*4 elev
        integer*4 ondate, offdat
c-----
c    VARIABLES IN COMMON BLOCK phases
c        sta - station name
c        ipwt    - integer phase weight 0-4
c        tp  - phase arrival time in seconds past midnight
c        xlat    - latitude  degrees north
c        xlon    - longitude degrees east
c        pfm - phase first motion x,X,c,C,d,D,+, or -
c        pph - phase classification, P PKP Lg
c        pqual   - phase arrival quality e,E,i,I
c        pres    - arrival time residual
c        pchan   - phase channel on station
c        pwt - final phase weight
c        ips(i,j)- (i,1) = table 1 for i th phase
c            - (i,2) = table 2 to compute phase(ips(i,1) - ips(i,2)
c        jain    - Phase takeoff angle at source in degrees
c        jaz - ray aziumth from source to receiver
c        elev    - elevation in kilometers
c        ondate  - 
c        offdat  -
c-----
        real*8 tt0
        real*4 ttlat, ttlon
        tt0 = 1.0d+37
        do 100 k=1,nr
            i = mapd(k)
            if(tt0.gt.tp(i) .and. ips(i,2).le.0)then
                tt0 = tp(i)
                ttlat = xlat(i)
                ttlon = xlon(i)
            endif
  100   continue
        if(tratbl.eq.'TELE')then
            if(inter)then
                if(verbose)then
            write(6,*)'Enter Teleseism Test Latitude and Longitude'
                endif
                if(fmask(1).eq.0 .or. fmask(2).eq.0)then
                    read(5,*)ttlat,ttlon
                endif
            endif
        endif
c-----
c    now apply the fixed criteria, as defined in gcmdln
c-----
        if(fmask(1).eq.0)tlat = ttlat
        if(fmask(2).eq.0)tlon = ttlon
        if(fmask(4).eq.0)t0 = tt0
        end

        subroutine setup(nr,verbose)
c-----
c    get phase data, associate with particular permitted phases
c----
c    nr  I*4 - number of phase observations to be used in location
c    verbose L   - output info
c-----
        integer nr
        logical verbose 

        parameter (NPHA=3000)
        common/phases/sta(NPHA),ipwt(NPHA),tp(NPHA),xlat(NPHA),
     1  xlon(NPHA),pfm(NPHA),pph(NPHA),pqual(NPHA),pres(NPHA),
     2  pchan(NPHA),pwt(NPHA),ips(NPHA,2),jain(NPHA),jaz(NPHA),
     3  jbaz(NPHA),asid(NPHA),arid(NPHA),mapd(NPHA),elev(NPHA), 
     4  ondate(NPHA),offdat(NPHA)
        integer*4 asid, arid, mapd
        character*6 sta
        integer*4 ipwt
        real*8 tp
        real*4 xlat,xlon
        character*2 pfm
        character*8 pph
        character*1 pqual
        real *4 pres,pwt
        character*2 pchan
        integer*4 ips
        real*4 jain,jaz,jbaz
        real*4 elev
        integer*4 ondate, offdat
c-----
c    VARIABLES IN COMMON BLOCK phases
c        sta - station name
c        ipwt    - integer phase weight 0-4
c        tp  - phase arrival time in seconds past midnight
c        xlat    - latitude  degrees north
c        xlon    - longitude degrees east
c        pfm - phase first motion x,X,c,C,d,D,+, or -
c        pph - phase classification, P PKP Lg
c        pqual   - phase arrival quality e,E,i,I
c        pres    - arrival time residual
c        pchan   - phase channel on station
c        pwt - final phase weight
c        ips(i,j)- (i,1) = table 1 for i th phase
c            - (i,2) = table 2 to compute phase(ips(i,1) - ips(i,2)
c        jain    - Phase takeoff angle at source in degrees
c        jaz - ray aziumth from source to receiver
c        elev    - elevation in kilometers
c        ondate  - 
c        offdat  -
c-----
        character str*20
        character ostr*20

        logical ext

        integer*4 date
        integer*4 doy, hour, minute, year, month
        integer*4 day
        real*4 second
        
c-----
c    FORMAT STATEMENTS
c-----
c    DO NOT WORRY ABOUT CARRIAGE CONTROL ON READ
c-----
    1   format(a6,1x,a2,1x,
     1  i4.4,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2,1x,f6.3,1x,
     1  i8,1x,a1,1x,a8,1x,i2,1x,a2,1x,i8,
     1  1x,f9.4,1x,f9.4,1x,f9.4,1x,i8,1x,i8)
c-----
c    UNIX FORTRAN WRITE - DO NOT USE CARRIAGE CONTROL
c-----
    4 format(a6,1x,i2,1x,a20,  1x,2i2,1x,a,1x,a,1x,a,1x,a,f9.4,
     1f10.4,f10.2,5x,a)
   44 format(a6,1x,i2,1x,f20.3,1x,2i2,1x,a,1x,a,1x,a,1x,a,f9.4,
     1f10.4,f10.2,5x,a)
    5   format('STA   ','IWT',7x,'ARRIVAL TIME ','  PhID','QL',
     1 ' PHASE  ','FM',' CHAN','    LAT   ',
     2 '    LON   ','   ELEV   ')
c-----
c    MICROSOFT FORTRAN WRITE - USE CARRIAGE CONTROL
c-----
C   4 format(1x.a6,1x,i2,1x,a20,  1x,2i2,1x,a,1x,a,1x,a,1x,a,f9.4,
C    1f10.4,f10.2,5x,a)
C  44 format(1x.a6,1x,i2,1x,f20.3,1x,2i2,1x,a,1x,a,1x,a,1x,a,f9.4,
C    1f10.4,f10.2,5x,a)
C   5   format(1x.'STA   ','IWT',7x,'ARRIVAL TIME ','  PhID','QL',
C    1 ' PHASE  ','FM',' CHAN','    LAT   ',
C    2 '    LON   ','   ELEV   ')
c-----
c    open data file of phase times
c-----
        inquire(file='elocate.dat',exist=ext)
        if(.not. ext)then
            call usage('phase file elocate.dat missing')
        endif
        open(1,file='elocate.dat',access='sequential',status='old',
     1          form='formatted')
            rewind 1
            nrr = 1
 1000   continue    
        ips(nrr,1) = -1
        ips(nrr,2) = -1
c-----
c    note the order in the IF statements is critical
c    for proper parsing
c-----
        read(1,1,end=2000,err=2000)sta(nrr), pchan(nrr), 
     1  year, month, day,hour,minute,second,
     1  arid(nrr), pqual(nrr), 
     2  pph(nrr), ipwt(nrr), pfm(nrr),asid(nrr),xlat(nrr), xlon(nrr), 
     3  elev(nrr), ondate(nrr), offdat(nrr)
        call htoe(year,month,day,hour,minute,second,tp(nrr))

            ipwt(nrr)=mod(ipwt(nrr),5)
c-----
c    parse the phase to get correct ips
c-----
            call pparse(pph(nrr),ips(nrr,1),ips(nrr,2))
                    
            nrr = nrr + 1
            goto 1000
 2000   continue
        nrr = nrr -1
        close(1)
            if(verbose)write(6,5)
            nr = 0
            do 37 i=1,nrr
                if(ips(i,2).eq.0)then
                    if(ips(i,1).gt.0)then
                        ostr = ' '
                        nr = nr + 1
                        mapd(nr) = i
                    else
                        ostr = 'NOT USED'
                    endif
                    call etoh(tp(i),date,str,doy)
                    if(verbose)write(6,4)sta(i),ipwt(i),str,
     1              ips(i,1),ips(i,2),
     2              pqual(i),pph(i),pfm(i),pchan(i),xlat(i),xlon(i),
     3              elev(i),ostr
                else if(ips(i,2).ne.0)then
                    if(ips(i,2).gt.0)then
                        ostr = ' '
                        nr = nr + 1
                        mapd(nr) = i
                    else
                        ostr = 'NOT USED'
                    endif
                    if(verbose)write(6,44)sta(i),ipwt(i),tp(i),
     1              ips(i,1),ips(i,2),
     2              pqual(i),pph(i),pfm(i),pchan(i),xlat(i),xlon(i),
     3              elev(i),ostr
                endif
   37       continue
        return
        end

        function lgstr(str)
c-----
c    find the length of a string, counting from end
c    toward beginning, stopping at the first non-blank
c    character, which has position str(lgstr:lgstr)
c----
        character str*(*)
        l = len(str)
        lgstr = 1
        do 100 i=l,1,-1
            if(str(i:i).ne.' ')then
                lgstr = i
                return
            endif
  100   continue
        return
        end

        subroutine pparse(pph,ips1,ips2)
c-----
c    parse the string pph to define the corresponding
c    travel time model to use, e.g., S-P requires S and P
c-----
c    pph C*8 - phase to be parsed
c    ips1    I*4 - model to use for first phase
c    ips2    I*4 - model to use for second phase
c        If one phase is found, e.g., P this pair is (1,0)
c        If one phase is not found, e.g., Ps this pair is (0,0)
c        If subtraction found, e.g., S-P this pair is (2,1)
c        If subtraction, but phase not found, e.g., S-Ps pair is (2,-1)
c        Here P=1, S=2
c-----
        character pph*8
        integer ips1, ips2
        integer lgstr
        integer gtthmd
        character tph*8

        l = lgstr(pph)
c-----
c    find minus sign if present
c-----
        iminus = 0
        do 100 i=1,l
            if(pph(i:i).eq.'-')then
                iminus = i
            endif
  100   continue
        if(iminus.eq.0)then
            tph = pph
            ips1 = gtthmd(tph,0)
            ips2 = 0
        else if(iminus.gt.1)then
            tph = pph(1:iminus-1)
            ips1 = gtthmd(tph,-1)
            tph = pph(iminus+1:l)
            ips2 = gtthmd(tph,-1)
        else
            ips1 = 0
            ips2 = 0
        endif
        if(ips1.lt.0 .or. ips2.lt.0)then
            ips1 = -1
            ips2 = -1
        endif
        return
        end

        function gtthmd(tph,ival)
        integer gtthmd
        character tph*8
        parameter (NPHA=3000)
        parameter (NMDPHA=6)
        common/a21/phathe(NMDPHA)
        character phathe*8
        common/a6/nl,nthpha
        gtthmd = ival
        do 100 i=1,nthpha
            if(tph .eq. phathe(i))then
                gtthmd = i
                return
            endif
  100   continue
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

        function wt(jj,dist,res,iter,ips,tratbl)
c-----
c    jj  -> hypo71 weight for arrival quality 0-4
c    dist    -> epicentral distance in km
c    res -> arrival time residual in sec
c    iter    -> number of iteration
c    ips -> 0 for P-wave, 1 for S-wave, 2 for Lg, 3 for Lg-P, 4 for S-P
c-----
        character tratbl*8
        data xnear,xfar/50.0,500.0/
        data tres/1.0/
c-----
c    simple weighting function
c-----
        j=mod(jj,5)
        wt = 1.0
c-----
c    use HYPO71 phase weight
c-----
        if(j.ge.0.and.j.le.4)then
            wt=wt*float(4-j)/4.
        else
            wt=0
        endif
c-----
c    reduce weight of S arrivals by vpvs = 1.732
c    to reduce sensitivity of results to s-waves
c-----
c    if(ips.gt.0)wt = wt/1.732
c-----
c do not do other weighting unless we are past
c the third iteration
c-----
        if(iter.le.3)return
c-----
c    perform distance weighting
c    distance weighting is very simple
c-----
        if(dist.gt.xnear.and.tratbl.ne.'TELE'.and.tratbl.ne.'BEAM')then
c-----
c    this weight =1.0 at xnear and 0.1 at xfar
c-----
            wt=wt*(xfar-xnear)/(9.*dist+xfar-10.*xnear)
        endif
c-----
c    perform a dumb residual weighting
c-----
        wt=wt* (tres/(tres+abs(res)))**2
        return
        end

        subroutine sumup(to,tlat,tlon,nr,tratbl,z,dlat,dlon,sdepth,
     1      stime,smajax,sminax,strike,fmask,fixed,rms)
c-----
c    subroutine to output results to data base
c
c    to  - origin epoch time 
c    tlat    - latitude degrees east
c    tlon    - longitude degrees west
c    nr  - number of phases used
c    tratbl  - travel time table used
c    z   - focal depth in km
c    fmask   - flag indicating what was fixed
c    fixed   - flag to indicate that part of the solution is fixed
c    rms - RMS arrival time error
c-----
        integer fmask(4)
        logical fixed

        parameter (NPHA=3000)
        common/phases/sta(NPHA),ipwt(NPHA),tp(NPHA),xlat(NPHA),
     1  xlon(NPHA),pfm(NPHA),pph(NPHA),pqual(NPHA),pres(NPHA),
     2  pchan(NPHA),pwt(NPHA),ips(NPHA,2),jain(NPHA),jaz(NPHA),
     3  jbaz(NPHA),asid(NPHA),arid(NPHA),mapd(NPHA),elev(NPHA), 
     4  ondate(NPHA),offdat(NPHA)
        integer*4 asid, arid, mapd
        character*6 sta
        integer*4 ipwt
        real*8 tp
        real*4 xlat,xlon
        character*2 pfm
        character*8 pph
        character*1 pqual
        real *4 pres,pwt
        character*2 pchan
        integer*4 ips
        real*4 jain,jaz,jbaz
        real*4 elev
        integer*4 ondate, offdat
c-----
c    VARIABLES IN COMMON BLOCK phases
c        sta - station name
c        ipwt    - integer phase weight 0-4
c        tp  - phase arrival time in seconds past midnight
c        xlat    - latitude  degrees north
c        xlon    - longitude degrees east
c        pfm - phase first motion x,X,c,C,d,D,+, or -
c        pph - phase classification, P PKP Lg
c        pqual   - phase arrival quality e,E,i,I
c        pres    - arrival time residual
c        pchan   - phase channel on station
c        pwt - final phase weight
c        ips(i,j)- (i,1) = table 1 for i th phase
c            - (i,2) = table 2 to compute phase(ips(i,1) - ips(i,2)
c        jain    - Phase takeoff angle at source in degrees
c        jaz - ray aziumth from source to receiver
c        elev    - elevation in kilometers
c        ondate  - 
c        offdat  -
c-----
      common/a2/delta(NPHA),dx(NPHA),dy(NPHA),t(NPHA) 
c-----
c    COMMON BLOCK a2
c        delta - epicentral distance in kilometers
c        dx    - delta latitude
c        dy    - delta longitude
c        t     - travel time in sec
c-----
        common/invert/sigma(4),c(4),co(4,4),a(4,4),b(4,4),y(4),sum2,
     1      nass,ndef
c-----
c    COMMON BLOCK invert
c        sigma(4) -
c        c(4)     -
c        co(4,4)  - variance-covariance matrix
c        a(4,4)   -
c        b(4,4)   -
c        y(4)     -
c        sum2     - standard deviation of observations
c-----
        real*8 to
        character tratbl*8
c    character atype*1
        character str*20
c-----
c    PLACE THE FOLLOWING INTO assoc
c
c    pres(i) ->  resid
c    jaz(i)  ->  esaz
c    jbaz(i) ->  seaz
c    tratbl  ->  tratbl
c    pwt ->  wgt
c    delta(i)/111.195
c        ->  delta   convert km to degrees
c-----
c-----
c    declarations needed for origin error
c-----
        character newln*1
        integer*4 iflag, orid, commid
        real*4 sdobs, sxx,syy, szz, stt, sxy, sxz, syz, stx, sty, stz
        real*4 conf
        character lddate*17
c-----
c    declarations needed for origin
c-----
        integer*4 date, nass, ndef, evid, grn, srn, doy
        integer*4 ndp
        integer*4 mbid, msid, mlid
        real*8 time
        real*4 depdp, lat, lon, depth, mb, ms, ml
        character algorithm*15, auth*15, dtype*1, etype*7
        character ltype*4
c-----
c    declarations needed for association
c-----
        character azdef*1, slodef*1, timedef*1, vmodel*8
        real*4 azres, slores
c-----
c    declarations needed for output
c-----
        character rq*2, lq*2
        character fm1*1, fm2*1
        newln(1:1) = char(10)
        lq = ' '''
        rq = ''' '
c-----
c    FORMAT STATEMENTS
c-----
c    DO NOT WORRY ABOUT CARRIAGE CONTROL ON READ
c-----
c-----
c    MS FORTRAN does not use carriage control on files
c-----
  100   format(i8,1x,10(f15.4,1x),3(f9.4,1x),f6.2,1x,f9.4,1x,f8.2,1x,
     1        f5.3,1x,i8,1x,a17,1x,i5,1x,a1)
  101   format(3(f9.4,1x),a18,1x,3(i8,1x),3(i4,1x),2(i8,1x),
     1  a7,1x,f9.4,1x,a1,1x,3(f7.2,1x,i8,1x),2(a15,1x),
     2  i8,1x,a17,1x,i5,4(f15.4,1x),a1)
  102   format(i8,1x,i8,1x,a6,1x,a8,1x,f4.2,1x,f8.3,1x,f7.2,1x,
     1  f7.2,1x,f8.3,1x,a1,1x,f7.1,1x,a1,1x,f7.2,1x,a1,1x,f7.1,
     2  1x,f6.3,1x,a15,1x,i8,1x,a17,1x,i5,1x,a1)
c-----
c    open and write to summary file for solution
c    this file basically states how the program worked
c    and would be used for a script based choice of the best solution
c-----
        open(1,file='elocate.sum',status='unknown',form='formatted'
     1      ,access='sequential')
        rewind 1
        call etoh(to,date,str,doy)
        write(1,'(f10.3,1x,2f12.4,f10.2,1x,a18,1x,f20.2)')
     1      rms, tlat,tlon,z,str(1:18),to
        close (1)
c-----
c    open and fill origin temporary file
c-----
        open(1,file='orgerr.tmp',status='unknown',form='formatted'
     1      ,access='sequential')
        rewind 1
        orid = -1
        sxx = co(1,1)
        syy = co(2,2)
        szz = co(3,3)
        stt = co(4,4)
        sxy = co(1,2)
        sxz = co(1,3)
        syz = co(2,3)
        stx = co(1,4)
        sty = co(2,4)
        stz = co(3,4)
        sdobs = sum2
        conf =  1.0
        commid = -1
        lddate = '-'
        iflag = 0
        write(1,100)orid,sxx,syy,szz,stt,sxy,sxz,syz,stx,sty,stz,sdobs,
     +  smajax,sminax,strike,sdepth,stime,conf,commid,lddate,
     +  iflag,newln 
        close(1)
c    close(1,status='save')
c-----
c    PLACE THE FOLLOWING INTO origin
c    
c    date    ->  date (get from DB here) ==> want lddate !!!
c    to  ->  time
c    tlat    ->  lat  
c    tlon    ->  lon
c    z   ->  depth
c    nr  ->  ndef
c    tratbl->    algorithm
c-----
        open(1,file='origin.tmp',status='unknown',form='formatted'
     1      ,access='sequential')
        rewind 1
        time = to
        call etoh(time,date,str,doy)
        lat =tlat
        lon = tlon
        depth = z
        algorithm = ' '
        algorithm(1:8) = tratbl
        mb = -9.99
        mbid = -1
        ms = -9.99
        msid = -1
        ml = -9.99
        mlid = -1
        ndp = 0
        depdp = -1.0
        orid = -1
        evid = -1
        grn = -1
        srn = -1
        ltype = 'LOCA'
c-----
c    if any parameters fixed so note
c-----
        if(fixed)then
            dtype = 'G'
        else
            dtype = 'R'
        endif
        etype= 'E'
        auth = '-'
        write(1,101)lat,lon,depth,str(1:18),orid,evid,date,nass,ndef,
     1  ndp,grn,srn,etype,depdp,dtype,mb,mbid,ms,msid,ml,mlid,
     2  algorithm,auth,commid,lddate,iflag,sxx,sxy,syy,szz,newln
        close(1)
c    close(1,status='save')
c-----
c    PLACE THE FOLLOWING INTO assoc
c
c    pres(i) ->  resid
c    jaz(i)  ->  esaz
c    jbaz    ->  seaz
c    tratbl  ->  tratbl
c    pwt ->  wgt
c    delta(i)/111.195
c        ->  delta   convert km to degrees
c-----
        open(1,file='assoc.tmp',status='unknown',form='formatted'
     1      ,access='sequential')
        rewind 1
        timedef = ' '
        azres = -9.99
        azdef = ' '
        slores = -9.99
        slodef = ' '
        emares = -9.99
        vmodel(1:8) = tratbl
        belief = 1
        do 200 k=1,nr
            i = mapd(k)
            esaz = jaz(i)
            seaz = jbaz(i)
            del = delta(i)/111.195
            write (1,102) arid(i),orid,sta(i),pph(i),belief,del,
     1      seaz,esaz,pres(i),timedef,azres,azdef,slores,
     2      slodef,emares,pwt(i),vmodel,commid,lddate,
     3      iflag,newln
  200   continue
c    close(1,status='save')
        close(1)
c-----
c    PLACE THE FOLLOWING INTO fmplot if the phase is P or p
c
c    jaz(i)  ->  
c    jain(i) ->  
c    pfm(i)  ->  iplmn
c    sta(i)  ->  
c-----
        open(1,file='fmplot.tmp',status='unknown',form='formatted'
     1      ,access='sequential')
        rewind 1
        do 300 k=1,nr
            i = mapd(k)
            if(pph(i)(1:1).eq.'P' .or. pph(i)(1:1).eq.'p')then
                fm1 = pfm(i)(1:1)
                fm2 = pfm(i)(2:2)
                if(fm1.eq.'C' .or. fm2.eq.'C')then
                    iplmn = 1
                else if(fm1.eq.'c' .or. fm2.eq.'c')then
                    iplmn = 1
                else if(fm1.eq.'+' .or. fm2.eq.'+')then
                    iplmn = 2
                else if(fm1.eq.'D' .or. fm2.eq.'D')then
                    iplmn = -1
                else if(fm1.eq.'d' .or. fm2.eq.'d')then
                    iplmn = -1
                else if(fm1.eq.'-' .or. fm2.eq.'-')then
                    iplmn = -2
                else
                    iplmn = 0
                endif
                iiaz = int( jaz(i) )
                iain = int(jain(i) )
C               if(iplmn.ne.0)then
                write (1,103) iiaz,iain,iplmn,lq//sta(i)//rq
C               endif
            endif
  300   continue
  103   format(i5,1x,i5,1x,i5,a)
        close(1)
        return
        end

        subroutine trtim(tlat,tlon,z,xlat,xlon,anin,
     1      delta,az,baz,ips1,ips2,tratbl,y,time,pph,elev,deldeg)
c-----
c    this subroutine calls the appropriate travel time routine
c    and sets up the travel times and partial derivatives
c-----
        integer*4 ips1, ips2
        character tratbl*8, pph*8
        real*4 y(4)
        real*4 time, t
        real*4 x(4)
c-----
c    compute source - receiver distance
c-----
        call delaz(tlat,tlon,xlat,xlon,deldeg,az,baz,delkm,saz,caz)
c-----
c    now given the travel times and partials
c    use the phase parser
c-----
        call vzero(y,4)
        time = 0.0d+00
        if(ips1.gt.0)then
            call trmat(tratbl,deldeg,delkm,saz,caz,t,x,ips1,z,pph,anin)
            call vsum(y,x,1,4)
            time = time + t
        endif
        if(ips2.gt.0)then
            call trmat(tratbl,deldeg,delkm,saz,caz,t,x,ips2,z,pph,anin)
            call vsum(y,x,-1,4)
            time = time - t
        endif
        if(tratbl.eq.'TELE')then
            delta = deldeg
        else if(tratbl.eq.'BEAM')then
            delta = deldeg
        else
            delta = delkm
        endif
        return 
        end 
        
        subroutine vzero(x,n)
        real*4 x(*)
        do 100 i=1,n
            x(i) = 0.0
  100   continue
        return
        end
        
        subroutine vsum(x,y,is,n)
        real*4 x(*), y(*)
        do 100 i=1,n
            if(is.gt.0)then
                x(i) = x(i) + y(i)
            else
                x(i) = x(i) - y(i)
            endif
  100   continue
        return
        end

        subroutine trmat(tratbl,deldeg,delkm,saz,caz,
     1      time,x,ips,z,pph,anin)
c-----
c    choose the earth model and compute the travel time and
c    partial derivatives
c-----
        character tratbl*8, pph*8
        real*4 deldeg, delkm, saz, caz, x(4)
        integer*4  ips
        real*4 time
        x(4) = 1.0
        if(tratbl.ne.'TELE' .and.tratbl.ne.'BEAM')then
            call trvdrv(delkm,z,time,dtdz,anin,caz,saz,x,ips) 
            continue
        else if(tratbl.eq.'TELE')then
            call jbtime(deldeg,caz,saz,time,dtdlat,dtdlon,dtdel
     1          ,dtdz,z,pph,anin)
            x(1) = dtdlon
            x(2) = dtdlat
            x(3) = dtdz
        else if(tratbl.eq.'BEAM')then
            call jbtime(deldeg,caz,saz,time,dtdlat,dtdlon,dtdel
     1          ,dtdz,z,pph,anin)
            x(1) = dtdlon
            x(2) = dtdlat
            x(3) = 0.0
        endif
        return
        end

       subroutine trvdrv(delta,z,t,dtdh,anin,caz,saz,x,ips) 
       implicit none
       real delta,z,t,dtdh,anin,caz,saz,x(4)
       integer ips

       integer NLAYR, NMDPHA
       parameter (NLAYR=100,NMDPHA=6)
       common/a20/v(NLAYR,NMDPHA),d(NLAYR),
     1       vsq(NLAYR,NMDPHA),thk(NLAYR),
     2       tid(NLAYR,NLAYR,NMDPHA),did(NLAYR,NLAYR,NMDPHA)
       real v, d, vsq, thk, tid, did
      
       common/a6/nl,nthpha
       integer nl, nthpha
       common/a22/f(NLAYR,NLAYR) 
       real f
       double precision q,qmax,dmdp0,temp,term1,sdel,sder,sum,p
       real tinj(NLAYR),didj(NLAYR),tr(NLAYR) 

c-----
c      internal variables
c-----
       real tdir, u, vmax, xovmax, tkj, zsq, tkjsq, sqt, srr, dtdd
       real tdj1, tmin
       real*8 deps
       integer l, j1, jl, jj, j, k, itno, m
c----------------------------------------------------------------------
c     The objective is to determine the first arrival time at a distance
c     delta for a source with depth z in the model. Onece this is 
c     determined also compute the dtdx. dtdz and angle of indidence
c-----
c     Logically this is as follows
c     1.  Determine source layer
c     2.  Do preliminary work to define the refraction arrivals
c         For the refractions we can them loop through all
c         possible refractions to get the one corresponding to the
c         shorted travel time at the desired distance
c     3   Consider the possible direct arrivals
c         If source in the first layer, this is a simple computation
c         For a source in a deeper layer, then this search must
c         be done numerically. 
c         Normally this search would be over ray parameters
c         [0,1/Vsource]. However the last would correspond to an
C         infinite distance.  Geometrically this last case 
c         corresponds to a ray leaving the source horizontally, but such
c         an arrival would be slower than the refraction along the
c         layer underneath. So in that case the direct ray would not be
c         the first arrival. 
c
      deps = 0.0002d00 
      zsq=z*z 
c-----initialization for fixed layer model ----------------------------
c-----
c    determine source layer
c-----
      do  1 l=1,nl 
            if (d(l) .gt. z) go to  2 
    1 continue 
      jl=nl 
      go to  3 
    2 jj=l 
      jl=l-1 
    3 continue
c-----
c    jl = source layer
c-----
        if(z.eq.d(jl)) then
            jj=jj-1
            jl=jl-1
        endif
         if(z.eq.0) then
            jj=2
            jl=1
        endif
      tkj=z-d(jl) 
      tkjsq=tkj**2+0.000001 
c------
c    jj =  layer beneath source
c    jl = source layer
c    if there is a refraction, adjust intercept times
c        and distance term for position of source in
c        the source layer. 
c
c        tinj(l) = intercept time for source layer jl
c            for refractor in layer l
c        didj(l) = distance of first critical arrival
c            for source in layer jl and refractor l
c            this is defined as the distance at which the
c            reflection from the top of refraction layer has
c            the ray parameter = 1./refractor velocity
c-----
      if (jl .eq. nl) go to  5 
c-----
c    adjust tinj, didj for the source position in the
c    source layer for valid refractions only
c-----
        do  l=jj,nl
c-----
c    determine maximum velocity in layer above refractor
c-----
            vmax = 0.0
            do j=1,l-1
                if(vmax.lt.v(j,ips))vmax = v(j,ips)
            enddo
            if(v(l,ips) .le. vmax)then
                tinj(l)= -1.0
                didj(l)= -1.0
            else
                sqt=sqrt(vsq(l,ips)-vsq(jl,ips)) 
                tinj(l)=tid(jl,l,ips)-tkj*sqt/(v(l,ips)*v(jl,ips)) 
                didj(l)=did(jl,l,ips)-tkj*v(jl,ips)/sqt 
            endif
      enddo
c-----
c    determine first possible refraction downward from the source
c    if monotonically increasing velocity, this just
c    invloves the source layer velocities and the refractor
c    velocity immediately beneath
c-----
        xovmax = -1.0
        do 404 l=jj,nl
            if(tinj(l).ge.0.0 .and.v(l,ips).ne.v(jl,ips))then
                xovmax=v(l,ips)*v(jl,ips)*(tinj(l)-tid(jl,jl,ips))/
     1              (v(l,ips)-v(jl,ips))
            endif
            if(xovmax.ge.0.0)goto 405
  404   continue
  405   if(xovmax.lt.0.0)xovmax = 1.0e+30
    5   continue
        tmin=999999.99
      if (jl .eq. nl) go to 100 
      do 60 m=jj,nl 
        if(tinj(m).ge.0.0)then
            tr(m)=tinj(m)+delta/v(m,ips)
        else
            tr(m)=1.0e+30
        endif
   60   continue
        vmax=0.0
        do 761 j=1,jl
            vmax=amax1(v(j,ips),vmax)
  761   continue
      do 70 m=jj,nl 
          k=m
          if (tr(m) .gt. tmin) go to 70 
          if(didj(m).lt.0)goto 70
          if (didj(m) .gt. delta) go to 70 
          k=m 
          tmin=tr(m) 
   70 continue 
      if (delta .lt. xovmax) go to 90 
c-----travel time + derivatives for refracted wave---------------------
   80 t=tr(k) 
      dtdd=1.0/v(k,ips) 
      dtdh=-sqrt(vsq(k,ips)-vsq(jl,ips))/(v(k,ips)*v(jl,ips)) 
      anin=-v(jl,ips)/v(k,ips) 
      go to 260 
c-----calculation for direct wave -------------------------------------
   90 if (jl .ne. 1) go to 100 
      sqt=sqrt(zsq+delta**2) 
      tdj1=sqt/v(1,ips) 
      if (tdj1 .ge. tmin) then
            k=1
            go to 80 
      endif
c-----travel time + derivatives for direct wave in first layer---------
      t=tdj1 
      dtdd=delta/(v(1,ips)*sqt) 
      dtdh=z/(v(1,ips)*sqt) 
      anin=delta/sqt 
      go to 260 
c-----find a direct wave that will emerge at the station---------------
  100 continue
c-----
c    for a ray to leave the source, it is bounded in
c    phase velocity between the fastest velocity in all
c    layers above and including the source layer on the
c    low end and by infinity on the high end
c-----
      j1=jl-1
      vmax=v(1,ips)
      do 175 j=2,jl
      if(v(j,ips).gt.vmax)vmax = v(j,ips)
  175 continue
c -- test to find maximum value for p
      qmax = 1.d0/vmax
      q = 0.5d0*qmax
  155 continue
      q = (q+qmax)/2.d0
      sdel = v(jl,ips)*tkj/dsqrt((1.d0-q*v(jl,ips))*(1.d0+q*v(jl,ips)))
      do 160 j=1,j1
      sdel =v(j,ips)*thk(j)/dsqrt((1.d0-q*v(j,ips))*(1.d0+q*v(j,ips)))
     1 +sdel
  160 continue
         dmdp0=delta-q*sdel
C      if (dmdp0) 165,161,155
       if(dmdp0 .eq. 0.0d+00)then
           go to 161
       else if(dmdp0 .gt. 0.0d+00)then
           go to 155
       endif
c -- now perform newton convergence from top down
         itno=0
  166 continue
         itno=itno+1
         if(itno.gt.30) stop
       sdel=v(jl,ips)*tkj/dsqrt((1.d0-q*v(jl,ips))*(1.d0+q*v(jl,ips)))
      sder = sdel/((1.d0-q*v(jl,ips))*(1.d0+q*v(jl,ips)))
      do 162 j=1,j1
      temp = dsqrt((1.d0-q*v(j,ips))*(1.d0+q*v(j,ips)))
      term1 = v(j,ips)*thk(j)/temp
      sdel = term1+sdel
      sder = term1/temp**2+sder
  162 continue
      dmdp0 = delta-q*sdel
      if(dabs(dmdp0).lt.deps) go to 161
      q = dmdp0/sder+q
      if(q.ge.qmax)q = qmax - 1.0d-05
      go to 166
  161 continue
c -- p has been determined to sufficient accuracy
c -- calculate direct wave travel time by summation
         p=q
      sum = tkj/(v(jl,ips)*dsqrt((1.d0-p*v(jl,ips))*(1.d0+p*v(jl,ips))))
        do 180 j=1,j1
            sum = thk(j)/(v(j,ips)*
     1      dsqrt((1.d0-p*v(j,ips))*(1.d0+p*v(j,ips)))) +sum
  180   continue
        tdir=real(sum)
        if(tdir.ge.tmin) go to 80
c-----travel time + derivatives for direct wave below first layer-------
      t=tdir 
      u=real(p)*v(jl,ips)
        srr=sqrt((1.-u)*(1.+u))
      dtdd=real(p)
c   dtdh=(1.0-v(jl)*u*dtdd)/(v(jl)*srr) 
c    I changed to
        dtdh=srr/v(jl,ips)
      anin=u 
c-----set up partial derivatives for regression analysis --------------
  260 x(1)=-dtdd*saz
      x(2)=-dtdd*caz
      x(3)=dtdh 
      return 
      end 

        subroutine gcmdln(xlat,xlon,depth,t0,model,fixed,
     1      verbose,fmask,auth,mxiter,inter,veryverbose)
c-----
c    Parse command line arguments
c    requires subroutine mgtarg() and function mnmarg()
c
c    If command line arguments are used, then processing will
c    be in a batch mode and not in an interactive mode
c
c    elocate -LAT fixed_latitude -LON fixed_longitude 
c        -DEPTH fixed_depth -T0 fixed_origin time
c        -M model_number -TELE -BEAM -F
c        -Q -TS time_string (199107112345.235)
c        -A authority -N maxiter -BATCH
c-----
c    xlat    R*4 Initial latitude    
c            (< -1000 means no given value)
c    xlon    R*4 Initial longitude   
c            (< -1000 means no given value)
c    depth   R*4 Initial depth
c            (< -1000 means no given value)
c    t0  R*4 Initial origin time
c    model   I*4 model
c    fixed   L   Initial values are fixed if .true.
c            (default false)
c    verbose L   Verbose mode (default is true)
c            Quiet mode turned on by -Q flag
c    fmask(4)I   Mask to define what parameters are fixed
c            fmask(1) == 1 latitude fixed
c            fmask(2) == 1 longitude fixed
c            fmask(3) == 1 depth fixed
c            fmask(4) == 1 origin time fixed
c    auth    C*15    character string for location authority
c    mxiter  I   maximum number of iterations (20=default)
c    inter   L   .true. interactive (default)
c            .false. non-interactive
c    veryverbose L .true. output variance covariance
c-----
c    subroutine arguments
c-----
        real*4 xlat, xlon, depth
        real*8 t0
        integer model, fmask(4), mxiter
        logical fixed, verbose, inter,veryverbose
        character auth*15
c-----
c    internal subroutine variables
c-----
        integer*4 year, month, day, hour, minute
        real*4 second
        integer mnmarg
        character*50 name
c-----
c    default initialization
c-----
        fixed = .false.
        verbose = .true.
        veryverbose = .false.
        inter = .true.
c-----
        fmask(1) = 0
        fmask(2) = 0
        fmask(3) = 0
        fmask(4) = 0
        auth = ' '
        model = -1
        xlat = -1000.0
        xlon = -1000.0
        depth = -1000.0
        mxiter = 20
        nmarg=mnmarg()
        i = 0
   11   i = i + 1
            if(i.gt.nmarg)goto 13
            call mgtarg(i,name)
            if(name(1:2).eq.'-M' .or. name(1:2).eq.'m')then
                i=i+1
                call mgtarg(i,name)
                read(name,'(bn,i10)')model
            else if(name(1:2).eq.'-N' .or. name(1:2).eq.'n')then
                i=i+1
                call mgtarg(i,name)
                read(name,'(bn,i10)')mxiter
            else if(name(1:4).eq.'-LAT' .or.
     1              name(1:4).eq.'-lat')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f20.0)')xlat
                fmask(1) = 1
            else if(name(1:4).eq.'-LON' .or.
     1              name(1:4).eq.'-lon')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f20.0)')xlon
                fmask(2) = 1
            else if(name(1:2).eq.'-D' .or. 
     1              name(1:2).eq.'-d')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f20.0)')depth
                fmask(3) = 1
            else if(name(1:3).eq.'-T0' .or. 
     1              name(1:3).eq.'-t0')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f20.0)')t0
                fmask(4) = 1
            else if(name(1:3).eq.'-TS' .or. 
     1              name(1:3).eq.'-ts')then
                i = i + 1
                call mgtarg(i,name)
c-----
c    TODO put in a string length check
c-----
                read(name,'(bn,i4,i2,i2,i2,i2,f6.3)')year,
     1              month,day,hour,minute,second
                call htoe(year,month,day,hour,minute,second,t0)
            else if(name(1:2).eq.'-F' .or. name(1:2).eq.'-f')then
                fixed = .true.
            else if(name(1:2).eq.'-Q' .or. name(1:2).eq.'-q')then
                verbose = .false.
            else if(name(1:3).eq.'-BA' .or. name(1:3).eq.'-ba')then
                inter = .false.
            else if(name(1:2) .eq. '-A' .or. name(1:2).eq.'-a')then
                i = i + 1
                call mgtarg(i,auth)
            else if(name(1:3).eq.'-TE' .or. 
     1              name(1:3).eq.'-te')then
                model = 1
            else if(name(1:3).eq.'-BE' .or. 
     1              name(1:3).eq.'-be')then
                model = 2
            else if(name(1:2).eq.'-?' .or. name(1:2).eq.'-h')then
                call usage(' ')
            else if(name(1:2).eq.'-V' .or. name(1:2).eq.'-v')then
                call musage()
            else if(name(1:2).eq.'-C' .or. name(1:2).eq.'-c')then
                veryverbose = .true.
            endif
        goto 11
   13   continue
        return
        end

        subroutine usage(str)
        character str*(*)
        integer LER, LOT, LIN
        parameter (LER=0, LOT=6, LIN=5)
        if(str .ne. ' ')write(LER,*)'elocate:',str
        write(LER,*)'Usage: elocate -M model -N itmax ',
     1      '-LAT elat -LON elon -D depth -T t0 ',
     2      '-TS YYYYMMDDHHMMSS.SSS -F -Q -BATCH -A auth ',
     3      '-TELE -BEAM -VELMOD -C'
        write(LER,*)
     1  ' -M model  (none) define model for interactive'
        write(LER,*)
     1  ' -N itmax (default 20 ) maximum iterations'
        write(LER,*)
     1  ' -LAT elat (starting latitude for teleseism or fixed lat)'
        write(LER,*)
     1  ' -LON elon (starting longitude for teleseism or fixed depth)'
        write(LER,*)
     1  ' -D depth  (depth for teleseism beam else fixed depth)'
        write(LER,*)
     1  ' -T t0     (fixed origin time in CSS format)'
        write(LER,*)
     1  ' -TS YYYYMMDDHHMMSS.SS   (fixed origin time in human format)'
        write(LER,*)
     1  ' -F        (default false) fix source depth'
        write(LER,*)
     1  ' -Q        (default false) force quiet mode'
        write(LER,*)
     1  ' -BATCH    (default false) not interactive'
        write(LER,*)
     1  ' -A auth                   comment for location'
        write(LER,*)
     1  ' -VELMOD   (default false) output VEL.MOD prototype'
        write(LER,*)
     1  '     Run this first to get the model file VEL.MOD in ',
     2  'the current directory'
        write(LER,*)
     1  ' -C        (default false) output covariance matrix'
        stop
        end

        subroutine musage()
        integer LER, LOT, LIN
        parameter (LER=0, LOT=6, LIN=5)
c----
c    output a prototype of VEL.MOD
c-----
        open(1,file='VEL.MOD',access='sequential',form='formatted',
     1      status='unknown')
        rewind 1
        write(1,'(a)')
     1      'HALF'
        write(1,'(a)')
     1      '1       2'
        write(1,'(a)')
     1      '     ''P''  ''S'''
        write(1,'(a)')
     1      '00.0 6.00 3.46'
        write(1,'(a)')
     1      'CUS'
        write(1,'(a)')
     1      '5       3'
        write(1,'(a)')
     1      '     ''P''  ''S''   ''Lg'''
        write(1,'(a)')
     1      '00.0 5.00 2.89 3.55'
        write(1,'(a)')
     1      '01.0 6.10 3.52 3.55'
        write(1,'(a)')
     1      '10.0 6.40 3.70 3.55'
        write(1,'(a)')
     1      '20.0 6.70 3.87 3.55'
        write(1,'(a)')
     1      '40.0 8.15 4.70 3.55'
        write(1,'(a)')
     1      'UPL'
        write(1,'(a)')
     1      '5       3'
        write(1,'(a)')
     1      '     ''P''  ''S''  ''Lg'''
        write(1,'(a)')
     1      '00.0 5.60 3.23 3.55'
        write(1,'(a)')
     1      '02.0 6.15 3.55 3.55'
        write(1,'(a)')
     1      '20.0 6.70 3.87 3.55'
        write(1,'(a)')
     1      '40.0 8.18 4.72 3.55'
        write(1,'(a)')
     1      '97.0 8.37 4.83 3.55'
        write(1,'(a)')
     1      'EMBN'
        write(1,'(a)')
     1      '6       4'
        write(1,'(a)')
     1      '     ''P''  ''Ps''  ''Sp''  ''S'''
        write(1,'(a)')
     1      ' 0.0 1.80 0.40 1.80 0.40'
        write(1,'(a)')
     1      ' 0.6 5.10 5.10 2.94 2.94'
        write(1,'(a)')
     1      ' 2.0 6.15 6.15 3.55 3.55'
        write(1,'(a)')
     1      '20.0 6.70 6.70 3.87 3.87'
        write(1,'(a)')
     1      '40.0 8.18 8.18 4.72 4.72'
        write(1,'(a)')
     1      '97.0 8.37 8.37 4.83 4.83'
        write(1,'(a)')
     1      'WUS'
        write(1,'(a)')
     1      '5       3'
        write(1,'(a)')
     1      '     ''P''  ''S''  ''Lg'''
        write(1,'(a)')
     1      ' 0.0  3.41 2.01 3.55'
        write(1,'(a)')
     1      ' 1.9  5.54 3.29 3.55'
        write(1,'(a)')
     1      ' 8.0  6.27 3.74 3.55'
        write(1,'(a)')
     1      '21.0  6.41 3.77 3.55'
        write(1,'(a)')
     1      '40.0  7.90 4.62 3.55'
        write(1,'(a)')
     1      'nnCIA '
        write(1,'(a)')
     1      '8 2 '
        write(1,'(a)')
     1      '    ''P''  ''S'' '
        write(1,'(a)')
     1      '0.0     3.7497     2.1436      '
        write(1,'(a)')
     1      '1.5     4.9399     2.8210      '
        write(1,'(a)')
     1      '4.5     6.0129     3.4336      '
        write(1,'(a)')
     1      '7.5     5.5516     3.1475      '
        write(1,'(a)')
     1      '14.5    5.8805     3.3583      '
        write(1,'(a)')
     1      '29.5    7.1059     4.0081      '
        write(1,'(a)')
     1      '35.5    7.1000     3.9864      '
        write(1,'(a)')
     1      '43.5    7.9000     4.4036     '
        write(1,'(a)')
     1      'KOREA '
        write(1,'(a)')
     1      '15    3 '
        write(1,'(a)')
     1      '            ''P''       ''S''      ''Lg'' '
        write(1,'(a)')
     1      '  0.00      5.38      3.00      3.50 '
        write(1,'(a)')
     1      '  1.00      5.81      3.24      3.50 '
        write(1,'(a)')
     1      '  2.00      6.17      3.44      3.50 '
        write(1,'(a)')
     1      '  3.00      6.29      3.51      3.50 '
        write(1,'(a)')
     1      '  6.00      6.32      3.53      3.50 '
        write(1,'(a)')
     1      ' 11.00      6.42      3.58      3.50 '
        write(1,'(a)')
     1      ' 16.00      6.56      3.66      3.50 '
        write(1,'(a)')
     1      ' 20.00      6.64      3.70      3.50 '
        write(1,'(a)')
     1      ' 25.00      6.65      3.71      3.50 '
        write(1,'(a)')
     1      ' 27.50      7.10      3.96      3.50 '
        write(1,'(a)')
     1      ' 30.00      7.92      4.41      3.50 '
        write(1,'(a)')
     1      ' 32.50      7.89      4.40      3.50 '
        write(1,'(a)')
     1      ' 35.00      7.87      4.39      3.50 '
        write(1,'(a)')
     1      ' 40.00      7.57      4.22      3.50 '
        write(1,'(a)')
     1      ' 45.00      7.76      4.33      3.50 '

        
        close(1)
        stop
        end
  


        subroutine beamit(t0,tlat,tlon,nr,dpth,tratbl,dlat,
     1      dlon,dz,dt,verbose,fmask,rms,inter)
c-----
c    This routine forms a beam for teleseismic P, calculates the
c    phase velocity and angle of approach across the array,
c    form which epicentral distance and location are determined
c
c    REFERENCE:
c    
c    Herrmann, R. B. (1982). Digital processing of 
c        regional network data,
c        Bull. Seism. Soc. Am. 72, S377-S392.
c-----
c    Assume a plane wave crossing the array (small network)
c    such that the arrival time i given by
c
c    F[ (p/c)x + (q/c)y - t ] where x is latitude (km), 
c        y is longitude (km)
c    and c is the phase velocity
c-----
c    subroutine arguments
c-----
        real*8 t0
        real*4 tlat, tlon, dpth, dlat, dlon, dz, dt
        integer nr, fmask(4)
        logical verbose, inter
        character tratbl*8
c-----
c    common variables
c-----
        parameter (NLAYR=100)
        common/event/iyr,imon,idy,ihr,imin,cmnt
        common/invert/sigma(4),c(4),co(4,4),a(4,4),b(4,4),y(4),sdobs,
     1      nass,ndef
        parameter (NPHA=3000)
        common/phases/sta(NPHA),ipwt(NPHA),tp(NPHA),xlat(NPHA),
     1  xlon(NPHA),pfm(NPHA),pph(NPHA),pqual(NPHA),pres(NPHA),
     2  pchan(NPHA),pwt(NPHA),ips(NPHA,2),jain(NPHA),jaz(NPHA),
     3  jbaz(NPHA),asid(NPHA),arid(NPHA),mapd(NPHA),elev(NPHA), 
     4  ondate(NPHA),offdat(NPHA)
        integer*4 asid, arid, mapd
        character*6 sta
        integer*4 ipwt
        real*8 tp
        real*4 xlat,xlon
        character*2 pfm
        character*8 pph
        character*1 pqual
        real *4 pres,pwt
        character*2 pchan
        integer*4 ips
        real*4 jain,jaz,jbaz
        real*4 elev
        integer*4 ondate, offdat
c-----
c    VARIABLES IN COMMON BLOCK phases
c        sta - stane setition name
c        ipwt    - integer phase weight 0-4
c        tp  - phase arrival time in seconds past midnight
c        xlat    - latitude  degrees north
c        xlon    - longitude degrees east
c        pfm - phase first motion x,X,c,C,d,D,+, or -
c        pph - phase classification, P PKP Lg
c        pqual   - phase arrival quality e,E,i,I
c        pres    - arrival time residual
c        pchan   - phase channel on station
c        pwt - final phase weight
c        ips(i,j)- (i,1) = table 1 for i th phase
c            - (i,2) = table 2 to compute phase(ips(i,1) - ips(i,2)
c        jain    - Phase takeoff angle at source in degrees
c        jaz - ray aziumth from source to receiver
c        elev    - elevation in kilometers
c        ondate  - 
c        offdat  -
c-----
        integer*4 date, doy
        character str*20
        common/a6/nl,nthpha
        common/a10/anin(NPHA)
        common/a2/delta(NPHA),dx(NPHA),dy(NPHA),tt(NPHA)
        real*4 tt
c-----
c    routine variables
c-----
        real*8 tref
c-----
c    initialization for possible error return
c-----
        rms = 100000.0
c-----
c    initialize the matrices
c-----
        do 100 i=1,4
            y(i) = 0.0
            do 101 j=1,4
                a(j,i) = 0.0
  101       continue
  100   continue
c-----
c    get fixed depth
c-----
        if(fmask(3).eq.0)then
            if(verbose)then
                write(6,*)' enter depth'
            endif
            read(5,*)dpth
            if(dpth.lt.0.0)then
                depth = abs(dpth)
            else
                depth = dpth
            endif
        else
            depth = abs(dpth)
        endif
c-----
c    first get working center of array of stations
c-----
        tlat = 0.0
        tlon = 0.0
        tref = 0.0d+00
        nrr = 0
        do 200 k=1,nr
            i=mapd(k)
            tlat = tlat + xlat(i)
            tlon = tlon + xlon(i)
            tref = tref + tp(i)
            nrr = nrr + 1
  200   continue
c-----
c    safety check
c-----
        if(nrr.eq.0)return
        tlat = tlat / nrr
        tlon = tlon / nrr
        tref = tref / nrr
c-----
c    for a fixed array, let xref and yref be the 
c    latitude and longitude of the
c    reference point of the array, so that all distances 
c        are with respect to
c    this point
c-----
c    here we use tlat and tlon
c-----
        xref = tlat
        yref = tlon
        call delaz(tlat,tlon,xref,yref,deldeg,az,tbaz,delkm,saz,caz)
        dxj = delkm * caz
        dyj = delkm * saz
c-----
c    now process to get phase velocity and azimuth of 
c        approach across array
c-----
        do 210 k=1,nr
            i=mapd(k)
            call delaz(tlat,tlon,xlat(i),xlon(i),deldeg,
     1          az,tbaz,delkm,saz,caz)
            dx(i) = delkm * caz
            dy(i) = delkm * saz
  210   continue
c-----
c    process to get  matrix
c-----
        do 300 k=1,nr
            i=mapd(k)
            w=wt(ipwt(i),deldeg,0.0,1,0,tratbl)
            pwt(i) = w
            xval = dx(i) - dxj
            yval = dy(i) - dyj
            tval = real( tp(i) - tref )
            a(1,1) = a(1,1) + (w * xval*xval)
            a(2,2) = a(2,2) + (w * yval*yval)
            a(1,2) = a(1,2) + (w * xval*yval)
            a(1,3) = a(1,3) + (w * xval)
            a(2,3) = a(2,3) + (w * yval)
            a(3,3) = a(3,3) + w
            y(1) = y(1) + (w * tval * xval)
            y(2) = y(2) + (w * tval * yval)
            y(3) = y(3) + (w * tval)
  300   continue
c-----
c    damp diagonal and make matrix symmetric
c-----
        a(1,1) = a(1,1) + (10.**(-4))
        a(2,2) = a(2,2) + (10.**(-4))
        a(3,3) = a(3,3) + (10.**(-4))
        a(4,4) = 1.0
        a(2,1) = a(1,2)
        a(3,1) = a(1,3)
        a(3,2) = a(2,3)
        a(1,4) = 0.0
        a(2,4) = 0.0
        a(3,4) = 0.0
        a(4,1) = 0.0
        a(4,2) = 0.0
        a(4,3) = 0.0
c-----
c    determine inverse of matrix a
c-----
        call matinv(a,4)
c-----
c    get the solution vector 
c-----
        do 310 i=1,3
            c(i) = (y(1) * a(i,1) + y(2) * a(i,2) + y(3) * a(i,3))
  310   continue
c-----
c    get reference time at the reference latitude and longitude
c-----
        tref = tref + c(3)
c-----
c    calculate dT/dX, c, p, q, baz
c    dT/dX = sec/km  c = km/sec
c    iain = angle of incidence at free surface (P vel = 6.0 km/s)
c-----
        pi = 3.1415927
        dtdd = sqrt((c(1)**2) + (c(2)**2))
        cvel = 1.0/dtdd
        pvel = 6.0
        rsin = dtdd * pvel
        iain = int(( 180.0/pi) * asin(rsin))
c-----
c    get backazimuth
c-----
        p = c(1) * cvel
        q = c(2) * cvel
        baz = atan2(-q, -p) 
        if(baz .lt.0.0)baz = baz + 2.0*pi
        bazdeg=baz*180./pi
c-----
c    convert dtddeg to sec/deg
c-----
        dtddeg = dtdd * 111.195
c-----
c    get epicentral distance from dtddeg and fixed depth
c-----
        call getdel(dtddeg,del,depth,iret)
        if(iret.lt.0)then
            return
        endif
c-----
c    now convert the back azimuth and distance to a 
c    location by using spherical trig
c-----
        call pos(qlat,qlon,xref,yref,bazdeg,del)
        call pos(qlat1,qlon1,xref,yref,bazdeg+0.05,del)
        call pos(qlat2,qlon2,xref,yref,bazdeg,del+0.05)
c-----
c    with this location, get true travel time from JB
c    to get true origin time
c-----
        call jbtime(del,1.0,0.0,time,dtdlat,dtdlon,dtdel2,
     1          dtdz,depth,'P       ',bnin)
        t0 = tref - dble(time)
        sse = 0.0
c-----
c    calculate residuals, by still assuming plane wave across array
c    one could in fact calculate the trune distance and then 
c        the times, but for
c    small arrays this might lead to a large nimerical error
c-----
        do 400 k=1,nr
            i=mapd(k)
            pres(i) = c(1) * (dx(i)-dxj) + c(2) * (dy(i)-dyj) 
     1           - real(tp(i) - tref)
            w=wt(ipwt(i),deldeg,0.0,1,0,tratbl)
            sse = sse + (w*pres(i)**2)
  400   continue
        var = sse / (nrr - 3)
        rms = sqrt(var)
        s1 = a(1,1)*var
        s2 = a(2,2)*var
        s3 = a(3,3)*var
        sigc = cvel*cvel*sqrt( p*p*s1 + q*q*s2)
        sigdtd = sigc/(cvel*cvel)
        call getdel(dtddeg+0.1,del1,depth,iret)
        sdel = abs(10.0*(del1-del))*sigdtd*111.195
        sigp=sqrt( ((p/cvel)*sigc)**2 + cvel*cvel*s1 )
        sigq=sqrt( ((q/cvel)*sigc)**2 + cvel*cvel*s2 )
        sbaz = sqrt((q*sigp)**2 + (p*sigq)**2)*180./pi
c-----
c    dlatdd,dlondd,dlatdb,dlondb are dimensionless
c-----
        dlatdd=20.*(qlat2-qlat)
        dlondd=20.*(qlon2-qlon)
c-----
c    baz is in radians but lat and lon are degrees
c    make dimensionless
c-----
        dlatdb=20.*(qlat1-qlat)*pi/180.
        dlondb=20.*(qlon1-qlon)*pi/180.
        dlat=sqrt((dlatdd*sdel)**2 + (dlatdb*sbaz)**2)
        dlon=sqrt((dlondd*sdel)**2 + (dlondb*sbaz)**2)
        if(qlon.gt.180.0)qlon=qlon-360.
        if(qlon.lt.-180.)qlon=qlon+360.
        tlat = qlat
        tlon = qlon
c-----
c    output the location
c-----
        dz = 0.0
        dt = sqrt(s3)
c-----
c    NOTE that some of this stuff could go into 
c        CSS DB for array locations
c-----
        write(6,10)
        do 990 k=1,nr 
            i=mapd(k)
            if(mod(k,20).eq.0 .and. verbose .and. inter)then
                write(6,*)'RETURN FOR NEXT SCREEN'
                read(5,'(1x)')
            endif
            call delaz(tlat,tlon,xlat(i),xlon(i),delta(i),
     1          az,baz,delkm,saz,caz)
            jaz(i) = az
            jain(i)=iain
            call etoh(tp(i),date,str,doy)
            if(verbose)then
                write(6,12)sta(i),pchan(i),delta(i),jaz(i),
     1              jain(i),str  ,pres(i),ipwt(i),
     2              pqual(i),pfm(i),pph(i),pwt(i)
            endif
  990   continue 
        if(verbose .and. inter)then
            write(6,*)'RETURN FOR NEXT SCREEN'
            read(5,'(1x)')
        endif
   10 format(1x,'STA  COMP DIS(D) AZM  AIN    ARR TIME',
     1  10x,'   RES(SEC) WT QFM PHASE    WGT')
   12 format(1x,a6,1x,a2,f7.2,2f5.0,1x,a20  ,1x,f10.2,2x,i1,1x, 
     1  a1,a2,a8,f5.2)
        if(abs(pres(i)).gt.9.99)pres(i)=sign(9.99,pres(i))
        call etoh(t0,date,str,doy)
        if(verbose)then
            call delaz(tlat-0.5,tlon+0.5,tlat+0.5,tlon-0.5,
     1          tdel,az,baz,delkm,saz,caz)
            dxx = delkm * saz
            dyy = delkm * caz
            dxx=abs(dxx) 
            dyy=abs(dyy) 
            dlt = dyy * dlat
            dln = dxx * dlon
            write(6,7)rms,tlat,dlat,dlt,tlon,dlon,dln,
     1          depth,dz,t0,dt,str(1:18),dt, 
     2          bazdeg,sbaz,dtddeg,sigdtd,del,sdel
        endif
    7 format(/' ',10x,'rms=',f10.3,' sec'/
     1  ' ',10x,f10.4,'+-',f10.4,'N',5x,f10.4,' km'/ 
     2  ' ',10x,f10.4,'+-',f10.4,'E',5x,f10.4,' km'/ 
     3  ' ',10x,f10.2,'+-',f10.2,' km'/
     4  ' ', 5x,f15.3,'+-',f10.2,' sec'/
     5  ' ', 2x,a18  ,'+-',f10.2,' sec'/
     6  ' ',10x,f10.4,'+-',f10.2,' Back Azimuth (deg)'/
     7  ' ',10x,f10.4 ,'+-',f10.3,' DTDD (sec/deg)'/
     8  ' ',10x,f10.4 ,'+-',f10.2,' Distance (deg)'/
     9  ' ') 
        return
        end

        subroutine getdel(dtddeg,del,depth,iret)
c-----
c    do a distance search and then a refinement 
c    here we jump through the 
c        trtim -> trmat -> jbtime sequence directly to jbtime
c
c    Note we ignore the ambiguity of having the array at a 
c        short distance from
c    a deep focus earthquake. In this case there will be 
c        two distances which
c    will have the same dt/dD. This is also implicit in the search
c-----
c    DEBUG since the jbtime does not do PKP, we limit out search to P
c-----
        icross = 0
        call jbtime(10.0,1.0,0.0,time,dtdlat,dtdlon,dtdel1,
     1          dtdz,depth,'P       ',anin)
        best = -1
        do 100 ideg = 11,105,1
            del = ideg
            call jbtime(del,1.0,0.0,time,dtdlat,dtdlon,dtdel2,
     1          dtdz,depth,'P       ',anin)
            if( dtddeg.le.dtdel1 .and. dtddeg.ge.dtdel2)then
                diff = (dtdel1 - dtdel2)/1.0
                if(diff .ne.0.0)then
                    best = (dtddeg - dtdel1)/diff + (del - 1.0)
                    ttime = time
                else
                    best = del
                endif
            endif
            dtdel1 = dtdel2
  100   continue
        if(best .lt.0)then
            iret = -1
        else
            del = best
            iret = 1
        endif
        return
        end

        subroutine pos(qlat,qlon,xref,yref,bazdeg,del)
c-----
c    given reference latitude and longitude, distance and back azimuth
c    use spherical trigonometry to get location
c-----
        pi = 3.1415927
c-----
c    conversion factors: condeg - convert to degrees
c    conrad - convert to radians
c-----
        condeg = 180. / pi
        conrad = pi / 180.
        baz = bazdeg * conrad
c-----
c    Using geometry, earthquake location is calculated
c    use eccentricity factor
c-----
        eccen = (1.-(1./297.))**2
        arg = eccen * tan(xref*conrad)
        argtmp = atan(arg)
        phipr = (90. * conrad) - argtmp
        cosb = cos(phipr)
        cosc = cos(del*conrad)
        sinb = sin(phipr)
        sinc = sin(del*conrad)
        cosA = cos(baz)
**********
        cosa = cosb * cosc + (sinb * sinc * cosA)
        a = acos(cosa)
        cosC = (cosc - cosa * cosb) / (sin(a) * sinb)
        CC =  acos(cosC)
*****************************
      qlat = 90. - (a * condeg)
******** convert from phi  
        arg1 = tan(qlat*conrad)/eccen
        qlat = atan(arg1)*condeg
**********
      if (baz .le. pi) then
             qlon = yref + (CC * condeg)
      else if (baz .gt. pi) then
             qlon = yref - (CC * condeg)
      end if
      if (baz .eq. 0.0 .and. (xref+del .eq. 90.)) then
          qlat = 90.
          qlon = 0.
      else if (baz .eq. 180. .and. (del-xref .eq. 90.)) then
          qlat = -90.
          qlon = 0.
      end if
        return
        end

        subroutine delaz(evla,evlo,stla,stlo,del,az,baz,delkm,saz,caz)
c-----
c    evla    R*4 Epicenter latitude (degrees)
c    evlo    R*4 Epicenter longitude (degrees)
c    stla    R*4 Station latitude (degrees)
c    stlo    R*4 Station longitude (degrees)
c    del R*4 Epicentral Distance (degrees)
c    az  R*4 Epicenter to Station Azimuth
c    baz R*4 Station to Epicenter Backazimuth
c    delkm   R*4 Epicentral Distance (km)
c    saz R*4 Sine of Azimuth
c    caz R*4 Cosine of Azimuth
c-----
        implicit none
        real evla,evlo,stla,stlo
        real del,az,baz,delkm,saz,caz
        double precision sevla,cevla,sevlo,cevlo,ea,eb,ec
        double precision sstla,cstla,sstlo,cstlo,sa,sb,sc

        double precision cdel, fac, tmp
        double precision cbz, sbz

        double precision degrad 
        double precision ellip
        double precision a, b, f, e, e2

        double precision xp1, xp2, xp3
        double precision yp1, yp2, yp3
        double precision zp1, zp2, zp3
        double precision xpe, zpe
        double precision xps, zps
        double precision re, rs 
        double precision bp, ap
        double precision te, ts
        double precision e2p
        double precision PI


c-----
c double a = 6378.137                    /* equatorial radius */
c double b = 6356.752                    /* polar radius = a(1-f) */
c double f = 0.003352810664747481        /* f = 1 - b/a =
c                                                1/298.257223563 */
c double e = 0.081819191                         /* eccentricity (A.38)
c                                                e^2 = (a^2 - b^2)/a^2
c                                                e^2 = 2f - f^2 */
c double e2 = 0.006694379990141317       /* e^2 */
c double eps = 0.006739496742276435      /* second ellipticity (A.42)
c                                                eps = (a^2 - b^2)/b^2
c                                                = e^2/(1-e^2) */
c static float re=6371.003
c-----
        data degrad/0.0174532925199433/
        data a/6378.137/,b/6356.752/,f/0.003352810664747481/
        data e/0.081819191/, e2/0.006694379990141317/
        data PI/3.141592653589793238512808959406186204433/

        

c-----
c      major stability hack
c-----
       if(evla .eq.  90.0)evla =  89.9999 
       if(evla .eq. -90.0)evla = -89.9999 
       if(stla .eq.  90.0)stla =  89.9999 
       if(stla .eq. -90.0)stla = -89.9999 
       if(evla .eq. 0.0 .and. stla .eq. 0.0)then
              stla =  0.0001 
              evla = -0.0001 
       endif
       if(( stla .eq. - evla ) .and. (abs( evla - stlo) .eq. 180.0))then
              stla = stla +  0.001 
              stlo = stla + 0.001 
       endif
       call dircos(evla,evlo,sevla,cevla,sevlo,cevlo,ea,eb,ec)
       call dircos(stla,stlo,sstla,cstla,sstlo,cstlo,sa,sb,sc)
c-----
c    compute distance
c    Choose correct formula for short and large distances
c-----
        cdel = (ea*sa + eb*sb + ec*sc)
c-----
c    if DEL = [0,20)
c-----
        if(cdel .gt. 0.9396)then
            fac = (ea-sa)**2 + (eb-sb)**2 + (ec-sc)**2
            fac = sqrt(fac)/2.0
            del = 2.0*asin(fac)
c-----
c    if DEL = [20,160]
c-----
        else if(cdel. le. 0.9396 .and. cdel. ge. -0.9396)then
            del = acos(cdel)
c-----
c    if DEL = (160,180]
c-----
        else
            fac = (ea+sa)**2 + (eb+sb)**2 + (ec+sc)**2
            fac = sqrt(fac)/2.0
            del = 2.0*acos(fac)
        endif
        
c-----
C     check for station or epicenter at pole
c-----
        if( del .eq. 0.0 )then
                az = 0.0 
                baz = 0.0 
                delkm = 0.0
                return
        endif
        if(evla .eq. 90.0 )then
                az = 360.0 - stlo 
                baz = 0.0
                saz = dsin(degrad*(360.0 - stlo))
                caz = dcos(degrad*(360.0 - stlo))
        else if(stla .eq. 90.0 )then
                az = 0.0
                baz = 360.0 - evlo
                saz =  0.0
                caz =  1.0
        else if(evla .eq. -90.0 )then
                az =  stlo
                baz = 180.0
                saz = dsin(degrad*(stlo)) 
                caz = dcos(degrad*(stlo)) 
        else if(stla .eq. -90.0 )then
                az = 180.0 
                baz = evlo 
                saz =  0.0 
                caz = -1.0 
        else 
                saz = cevla*(cstla * sin(degrad*(stlo - evlo)))
                caz = (sstla - cdel*sevla)
                fac = sqrt((saz)*(saz) + (caz)*(caz))
                if(fac .gt. 0.0)then
                        saz = saz / fac
                        caz = caz / fac
                        az = atan2(saz,caz)
                
                        sbz = -cstla*(cevla * sin(degrad*(stlo - evlo)))
                        cbz = (sevla - cdel*sstla)
                        baz = atan2(sbz,cbz)
                else 
                        az = 0.0
                        caz = 1.0
                        saz = 0.0
                        baz = 180.0
                endif
                az  = real(az / degrad )
                baz = real(baz / degrad )
        endif
        del = real(del / degrad )

c-----
c    put az and baz in the range [0,360)
c-----
        if( az  .lt.   0.0)az  = az  + 360.0 
        if( az  .ge. 360.0)az  = az  - 360.0 
        if( baz .lt.   0.0)baz = baz + 360.0 
        if( baz .ge. 360.0)baz = baz - 360.0 

c-----
c       to compute the distance in kilometers we
c              follow the idea of Rudoe (Bomford, 1980) but do it
c              algebraically as much as possible.
c              1. E and S define a plane which goes through the origin
c              2. The intersection of this plane with the spheroid is
c                an ellipse whose major axis is 'a' but whose minor axis is
c                b <=  minor_axis <= a
c              3. Now we just define  latitudes on the ellipse 
c              4. Get the distance by an elliptic integral 
c                     or approximate formula
c       
c        define the plane:  alpha x + beta y + gamma z = 0 since
c              plane goes through the origin. The vector normal to the
c              plane is ( alpha, beta, gamma) = E x S
c-----        
       re = a*sqrt(1.0 -e2*sevla*sevla)
       rs = a*sqrt(1.0 -e2*sstla*sstla)
       yp1 =  eb * sc - sb * ec 
       yp2 =  ec * sa - ea * sc 
       yp3 =  ea * sb - eb * sa 
c-----
c        safety
c-----
       fac = sqrt ( yp1*yp1 + yp2*yp2 + yp3*yp3)
       if(abs(yp1) .lt.  1.0e-6 * fac) then
              yp1 =  1.0e-6 
       endif
       if(abs(yp2) .lt.  1.0e-6 * fac)then
              yp2 = -1.0e-6 
       endif
       fac = sqrt ( yp1*yp1 + yp2*yp2 + yp3*yp3)
       yp1 = yp1 / fac 
       yp2 = yp2 / fac 
       yp3 = yp3 / fac
c-----
c     since this defines the normal, now define other unit vectors
c-----
       xp1 =  yp2
       xp2 = -yp1 
       xp3 = 0.0 
       fac = sqrt ( xp1*xp1 + xp2*xp2 + xp3*xp3)
       xp1 = xp1 / fac 
       xp2 = xp2 / fac 
       xp3 = xp3 / fac

       zp1 = - yp1*yp3 
       zp2 = - yp2*yp3 
       zp3 = yp1*yp1 + yp2*yp2 
       fac = sqrt ( zp1*zp1 + zp2*zp2 + zp3*zp3)
       zp1 = zp1 / fac 
       zp2 = zp2 / fac 
       zp3 = zp3 / fac

c-----
c by construction the xp axis is horizontal - this has
c    an ellipse width of 2*ap = 2*a
c    the xp axis is the director of the minor axis and
c    has length of 2*bp
c    we solve the ellipse equation
c-----
       ap = a 
       bp = 1./sqrt((zp1/a)*(zp1/a) + (zp2/a)*(zp2/a) + (zp3/b)*(zp3/b))
       e2p = (1.0 - bp/ap)*(1.0 + bp/ap)

c-----
c   now get the points on the new ellipse E(3d) -> E'(2d)
c        S(3d) -> S'(2d) 
c-----

       xpe = ea*xp1 + eb*xp2 + ec*xp3 
       zpe = ea*zp1 + eb*zp2 + ec*zp3 
       te = atan2(zpe/bp , xpe/ap)
       xps = sa*xp1 + sb*xp2 + sc*xp3 
       zps = sa*zp1 + sb*zp2 + sc*zp3 

       ts = atan2(zps/bp , xps/ap)
c-----
c       compute the linear distance by integrating along the path.
c              However be careful to use the minor arc and to use
c              increasing angle 
c-----  
       if ( ts .gt. te )then
              tmp = ts 
              ts = te 
              te = tmp 
       endif
c-----
c now check minor arc 
c-----
       if( (te -ts ) .le. PI )then
c                minor arc 
              delkm = real(a) * ellip(e2p, (ts)/degrad,(te)/degrad) 
       else
              delkm = real(a) * ellip(e2p, (te)/degrad,180.0d+00) 
     1                  + real(a) * ellip(e2p, -180.0d+00,(ts)/degrad) 
       endif
       return
       end

        function ellip(k2, t1, t2)
        implicit none
        double precision ellip
        double precision k2, t1, t2
        integer N

        double precision degrad 
        double precision dt, t, ct, st, dct, dst
        integer i
        double precision sum
        double precision st2
        double precision ct2
        double precision tmp
        data degrad/ 0.0174532925199433/

        data N/500/
c-----
c      compute incomplete elliptic integral 
c           t2          2    2    1/2
c        INT   ( 1.0 - k cos t  )  dt
c           t1
c        NOTE I USE cos AND NOT sin BECAUSE OF MY COORDINATES
c-----

        sum = 0.00
        dt = (t2-t1)*degrad/N 
        st = dsin(t1*degrad)
        ct = dcos(t1*degrad)
        dct = dcos(dt)
        dst = dsin(dt)
        t = t1 
c-----
c                Use rectangular rule
c-----
        do i=0,N
                st2 = st*st 
                ct2 = ct*ct 
                if(i .eq. 0)then
                        sum = sum +  0.5*dt*sqrt(1.0-ct2*k2)
                else if(i .eq. N)then
                        sum = sum +  0.5*dt*sqrt(1.0-ct2*k2)
                else
                        sum = sum +  dt*sqrt(1.0-ct2*k2) 
                endif
                tmp = ct*dct - st*dst
                st  = st*dct + ct*dst
                ct = tmp
                t= t + dt
        enddo
        ellip = sum
        return
        end
        
        subroutine dircos(lat,lon,slat,clat,slon,clon,aa,bb,cc)
        implicit none
        real lat,lon
        double precision slat,slon,clat,clon, aa, bb, cc
        double precision degrad, e2
        double precision c, s, fac, e4
        data degrad/0.0174532925199433/,e2/0.006694379990141317/
c-----
c    convert geographic latitude to geocentric
c    Use flattening of Chovitz (1981) f= 1/298.257 adopted by IUGG 1980
c    
c    The revlaion between geocentric and geographic latitude is
c    tan phi c = ( 1 - f)^2 tan phi g
c
c    To avoid problems at poles, define sin phi c and cos phi c
c    so that the relation holds and also that s^2 + c^2 = 1
c-----
        c = dcos(degrad*lat)
        s = dsin(degrad*lat)
        e4 = (1.-e2)*(1.-e2)
        fac = dsqrt(e4 + (1.0-e4)*c*c)
        slat = (1.-e2) * s /fac
        clat =      c /fac
        slon = dsin(degrad*lon)
        clon = dcos(degrad*lon)

        aa = clat * clon
        bb = clat * slon
        cc = slat
        
        return
        end

