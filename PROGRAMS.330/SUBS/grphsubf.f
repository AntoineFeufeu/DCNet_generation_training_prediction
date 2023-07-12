c-----
c       changes
c
c       01 JAN 2004 adjusted position of Y-axis title in doliny dolny 
c           for case when numbers are negative
c       27 JUL 2020  - permit inverted triangle 'IT'
c       31 JAN 2020  - modify the choice of increments on the linear
c          axes to avoid entries like 70 140 210 instead of 50 150 250
c
c       limitations:
c       For logarithmic axes that do not extend over one order of
c          magnitude, there may not be a number label, e.g.,
c          from 2 to 8
c       For linear axes, there may ne be a good plot if the
c           limits are too close. This is related to the maximum
c           number of significal figures. This 4.0 to 4.05 should work
c           but 4.00000 to  4.0001 will not. The solution here is to
c           X - 4   instead of X
c----
        subroutine gbox(xl,yl,xh,yh)
c-----
c       draw a box bounded by the lower left and upper right corners
c
c       xl   R - x-coordinate of lower left
c       yl   R - y-coordinate of lower left
c       xh   R - x-coordinate of upper right
c       yh   R - y-coordinate of upper right
c-----
        implicit none
        real xl,yl,xh,yh
        call plot(xl,yl,3)
        call plot(xh,yl,2)
        call plot(xh,yh,2)
        call plot(xl,yh,2)
        call plot(xl,yl,2)
        return
        end

        subroutine gcent(xx,yy,ht,string,angle)
c-----
c       center a text string
c
c       xx       R - x-coordinate for center
c       yy       R - y-coordinate for center
c       ht       R - height of string. Really from
c                   baseline to top
c       string   C - string to be plotted
c       angle    R - angle in degrees with respect to 
c               horizontal z-axis
c-----
        implicit none
        real xx,yy,ht,angle
        character string*(*)

        integer lgstr

        integer il
        real t, ct,st,rl, xxx, yyy

            il = lgstr(string)
            t = angle*3.1415927/180.0
            ct = cos(t)
            st = sin(t)
            rl = 0.5*il*ht
            
            xxx = xx - rl*ct
            yyy = yy - rl*st
            call symbol(xxx,yyy,ht,string,angle,il)
        return
        end

        subroutine gright(xx,yy,ht,string,angle)
c-----
c       right justify a text string
c
c       xx       R - x-coordinate for center
c       yy       R - y-coordinate for center
c       ht       R - height of string. Really from
c                   baseline to top
c       string   C - string to be plotted
c       angle    R - angle in degrees with respect to 
c               horizontal z-axis
c-----
        implicit none
        real xx,yy,ht,angle
        character string*(*)

        integer lgstr

        integer il
        real t, ct,st, rl, xxx, yyy

            il = lgstr(string)
            t = angle*3.1415927/180.0
            ct = cos(t)
            st = sin(t)
            rl = 1.0*il*ht
            
            xxx = xx - rl*ct
            yyy = yy - rl*st
            call symbol(xxx,yyy,ht,string,angle,il)
        return
        end

        subroutine gleft(xx,yy,ht,string,angle)
c-----
c       left justify a text string
c
c       xx       R - x-coordinate for center
c       yy       R - y-coordinate for center
c       ht       R - height of string. Really from
c                   baseline to top
c       string   C - string to be plotted
c       angle    R - angle in degrees with respect to 
c               horizontal z-axis
c-----
        implicit none
        real xx,yy,ht,angle
        character string*(*)

        integer lgstr

        integer il
        real t, ct,st, rl, xxx, yyy

            il = lgstr(string)
            t = angle*3.1415927/180.0
            ct = cos(t)
            st = sin(t)
            rl = 0.0*il*ht
            
            xxx = xx - rl*ct
            yyy = yy - rl*st
            call symbol(xxx,yyy,ht,string,angle,il)
        return
        end

        subroutine dology(x0,y0,yleng,symax,symin,
     1      sizey,ticlft,lablft,dopow,ly,titley)
c-----
c       plot a logarithmic y-axis
c
c       x0      R - position of bottom side of axis
c       y0      R - position of bottom side of axis
c       yleng   R - length of Y axis
c       symax   R - maximum value of number corresponding to
c                   top
c       symin   R - minimum value of number corresponding to
c                   bottom
c       sizey   R - size of numbers, e.g., 10
c       ticlft  L   - .true. tics go left
c                 .false. tics go right
c       lablft  L   - .true. 10**N goes to left of axis
c                 .false. goes to right
c       dopow   L   - .true. put in 10**N labels
c       ly      I length of y-axis title string
c       titley  C y-axis title string
c-----
        real x0, y0, yleng, ymax
        integer nocy
        logical ticlft, lablft, dopow
        integer ly
        character titley*(*)
        character title*80

        real xnum(9)
        data xnum/1.0 ,2.0 ,3.0 ,4.0 ,5.0 ,6.0 ,7.0 ,8.0 ,9.0/

        ymin = symin
        ymax = symax
c-----
c       put in the axis
c-----
        ycen = y0 + 0.5*yleng
        call plot(x0,y0,3)
        call plot(x0,y0+yleng,2)
c-----
c       set up positioning for 10**N labels
c-----
c-----
c       compute limits for scale
c-----
        ymxlog = alog10(ymax)
        ymmin = alog10(ymin)
        nocy = INT(ymxlog - ymmin + 1)
        iy = INT(ymmin)
        if(ymmin .lt. 0)iy = iy - 1

        ymlog = alog10(ymax)
        iiy = INT(ymmin)
        if(ymlog .lt. 0)iiy = iiy - 1

        ja = max(abs(iy),abs(iiy))
        if(ja.lt.10)then
            jpow = 1
        else if(ja.ge.10 .and. ja.lt.100)then
            jpow = 2
        else if(ja.ge.100 .and. ja.lt.1000)then
            jpow = 3
        endif
        ii = min(iy,iiy)
        if(ii.lt.0)jpow = jpow + 1
        jpow = jpow + 2
        xshft = (2 + jpow*0.7)*sizey

        if(lablft)then
            if(ticlft)then
                xpos = x0 - 0.50*sizey - xshft
                xpos2 = x0 - 0.50*sizey - xshft -1.2*sizey
                xpos1 = xpos + 0.7*sizey
            else
                xpos = x0 + 0.50*sizey - xshft
                xpos2 = x0 + 0.50*sizey - xshft -1.2*sizey
                xpos1 = xpos + 0.7*sizey
            endif
        else
            if(ticlft)then
                xpos = x0 + 1.00*sizey
                xpos2 = x0 + 0.00*sizey + 0.2*sizey + xshft
                xpos1 = xpos + 0.7*sizey
            else
                xpos = x0 + 2.0*sizey
                xpos2 = x0 + 1.0*sizey + 0.2*sizey + xshft
                xpos1 = xpos + 0.7*sizey
            endif
        endif

        ifirst = 1
        ylow = y0
        yscal = yleng/alog10(ymax/ymin)
        do ii=iy,iy+nocy+2
            tenpow = 10.0**ii
            do jj=1,9
                yval = xnum(jj)*tenpow
c-----
c               put in tics where appropriate
c-----
                if(yval .ge. ymin .and. yval.le.ymax)then
                    yv = alog10(yval/ymin)
                    yy = y0 + yv*yscal
                    if(jj.eq.1)then
                        ticlen = 1.5*sizey
c-----
c                       put in 10**N
c-----
                        if(ifirst.eq.1)then
                            ypos = yy
                        else
                            ypos = yy - 0.5*sizey
                    if(yval.eq.ymax)ypos=yy-1.4*sizey
                            if(ypos.lt.y0)ypos = y0
                        endif
                        if(dopow)then
                        ypos1 = ypos + 0.7*sizey
            call symbol(xpos,ypos,sizey,'10',0.0,2)
            call number(999.,ypos1,0.7*sizey,real(ii),0.0,-1)
                        endif
                    else
                        ticlen = sizey
                    endif
                        ifirst = ifirst + 1
                    call plot(x0,yy,3)
                    if(ticlft)then  
                        call plot(x0-ticlen,yy,2)
                    else
                        call plot(x0+ticlen,yy,2)
                    endif
                endif
            enddo
            call plot(x0,yy,3)
        enddo
c-----
c       put in the title if we put in the numbers
c-----
        if(dopow .and. ly.gt.0)then
            title = ' '
            title(1:ly) = titley(1:ly)
            if(lablft)then
                call gcent(xpos2,ycen,1.2*sizey,title, 90.0)
            else
                call gcent(xpos2,ycen,1.2*sizey,title,-90.0)
            endif
        endif
        return
        end

        subroutine dologx(x0,y0,xleng,sxmax,sxmin,
     1      sizex,ticup,labtop,dopow,lx,titlex)
c-----
c       put in tic marks along the x-axis 
c
c       x0  R - position of left side of axis
c       y0  R - position of left side of axis
c       xleng   R - length of X axis
c       xmax    R - maximum value of number corresponding to
c                   far right
c       xmin    R - maximum value of number corresponding to
c                   far left
c       sizex   R - size of numbers, e.g., 10
c       ticup   L   - .true. tics go up
c                 .false. tics go down
c       labtop  L   - .true. 10**N goes above axis
c                 .false. goes below
c       dopow   L   - .true. put in 10**N labels
c       lx  I length of x-axis title string
c       titlex  Ch  x-axis title string
c-----
        real x0, y0, xleng, xmax
        integer nocx
        logical ticup, labtop, dopow
        integer lx
        character titlex*(*)
        character title*80

        real xnum(9)
        data xnum/1.0 ,2.0 ,3.0 ,4.0 ,5.0 ,6.0 ,7.0 ,8.0 ,9.0/

        xmin = sxmin
        xmax = sxmax
c-----
c       put in the axis
c-----
        xcen = x0 + 0.5*xleng
        call plot(x0,y0,3)
        call plot(x0+xleng,y0,2)
c-----
c       set up positioning for 10**N labels
c-----
        if(labtop)then
            if(ticup)then
                ypos  = y0   + 2.0 * sizex
                ypos2 = y0   + 3.7 * sizex
                ypos1 = ypos + 0.7*sizex
            else
                ypos  = y0   + sizex
                ypos2 = y0   + 2.7*sizex
                ypos1 = ypos + 0.7*sizex
            endif
        else
            if(ticup)then
                ypos  = y0   - 2.0*sizex
                ypos2 = y0   - 3.7*sizex
                ypos1 = ypos + 0.7*sizex
            else
                ypos  = y0   - 3.0*sizex
                ypos2 = y0   - 4.7*sizex
                ypos1 = ypos + 0.7*sizex
            endif
        endif
c-----
c       compute limits for scale
c-----
C       xmin = xmax / 10.0**nocx
        xmxlog = alog10(xmax)
        xmmin = alog10(xmin)
        nocx = INT(xmxlog - xmmin + 1)
        ix = INT(xmmin)
        if(xmmin .lt. 0)ix = ix - 1
        ifirst = 1
        do ii=ix,ix+nocx+2
            tenpow = 10.0**ii
            ja = abs(ii)
            if(ja.lt.10)then
                jpow = 1
            else if(ja.ge.10 .and. ja.lt.100)then
                jpow = 2
            else if(ja.ge.100 .and. ja.lt.1000)then
                jpow = 3
            endif
            if(ii.lt.0)jpow = jpow + 1
            xscal = xleng/alog10(xmax/xmin)
            do jj=1,9
                xval = xnum(jj)*tenpow
c-----
c               put in tics where appropriate
c-----
                if(xval .ge. xmin .and. xval.le.xmax)then
                    xv = alog10(xval/xmin)
                    xx = x0 + xv*xscal
                    if(jj.eq.1)then
                        ticlen = 1.6*sizex
c-----
c                       put in 10**N
c-----
                    if(ifirst.eq.1)then
                            xpos = xx
                    else
                        xshft = (2 + jpow*0.7)*sizex
                        if(xx+xshft .gt. x0+xleng)then
                            xpos = x0+xleng - xshft
                        else
                            xpos = xx - xshft/2.0
                        endif
                        if(xpos.lt.x0)xpos = x0
                    endif
                    if(dopow)then
                    call symbol(xpos,ypos,sizex,'10',0.0,2)
            call number(999.,ypos1,0.7*sizex,real(ii),0.0,-1)
                        endif
                    else if(jj.eq.5)then
                        ticlen = 1.3*sizex
                    else
                        ticlen = sizex
                    endif
                        ifirst = ifirst + 1
                    call plot(xx,y0,3)
                    if(ticup)then   
                        call plot(xx,y0+ticlen,2)
                    else
                        call plot(xx,y0-ticlen,2)
                    endif
                endif
            enddo
            call plot(xx,y0,3)
        enddo
        if(dopow .and. lx.gt.0)then
            title = ' '
            title(1:lx) = titlex(1:lx)
            call gcent(xcen,ypos2,1.2*sizex,title,0.0)
        endif
        return
        end

        subroutine dnlinx(x0,y0,xleng,xmax,xmin,
     1      sizex,ticup,labtop,dopow,llx,titlex)
c-----
c       put in tic marks along the x-axis but with
c       maximum to the left and the minimum to the right
c
c       x0  R - position of left side of axis
c       y0  R - position of left side of axis
c       xleng   R - length of X axis
c       xmax    R - maximum value of number corresponding to
c                   far right
c       xmin    R - maximum value of number corresponding to
c                   far left
c       sizex   R - size of numbers, e.g., 10
c       ticup   L   - .true. tics go up
c                 .false. tics go down
c       labtop  L   - .true. number goes above axis
c                 .false. goes below
c       dopow   L   - .true. put in number labels
c       llx I length of x-axis title string
c       titlex  Ch  x-axis title string
c-----
        real x0, y0, xleng, xmax, xmin, sizex
        logical ticup, labtop, dopow
        integer llx
        character titlex*(*)
        call dolnx(x0,y0,xleng,xmax,xmin,
     1      sizex,ticup,labtop,dopow,llx,titlex,.false.)
        return
        end

        subroutine dolinx(x0,y0,xleng,xmax,xmin,
     1      sizex,ticup,labtop,dopow,llx,titlex)
c-----
c       put in tic marks along the x-axis but with
c       minimum to the left and the maximum to the right
c
c       x0  R - position of left side of axis
c       y0  R - position of left side of axis
c       xleng   R - length of X axis
c       xmax    R - maximum value of number corresponding to
c                   far right
c       xmin    R - maximum value of number corresponding to
c                   far left
c       sizex   R - size of numbers, e.g., 10
c       ticup   L   - .true. tics go up
c                 .false. tics go down
c       labtop  L   - .true. number goes above axis
c                 .false. goes below
c       dopow   L   - .true. put in number labels
c       llx I length of x-axis title string
c       titlex  Ch  x-axis title string
c-----
        real x0, y0, xleng, xmax, xmin, sizex
        logical ticup, labtop, dopow
        integer llx
        character titlex*(*)
        call dolnx(x0,y0,xleng,xmax,xmin,
     1      sizex,ticup,labtop,dopow,llx,titlex,.true.)
        return
        end

        subroutine dolnx(x0,y0,xleng,xmax,xmin,
     1      sizex,ticup,labtop,dopow,llx,titlex,doasis)
c-----
c       put in tic marks along the x-axis 
c       this routine does all of the work
c
c       x0  R - position of left side of axis
c       y0  R - position of left side of axis
c       xleng   R - length of X axis
c       xmax    R - maximum value of number corresponding to
c                   far right
c       xmin    R - maximum value of number corresponding to
c                   far left
c       sizex   R - size of numbers, e.g., 10
c       ticup   L   - .true. tics go up
c                 .false. tics go down
c       labtop  L   - .true. number goes above axis
c                 .false. goes below
c       dopow   L   - .true. put in number labels
c       llx I length of x-axis title string
c       titlex  Ch  x-axis title string
c       doasis  L   .true. plot values
c                   .false. annotate with negative (useful for depth plot)
c-----
        implicit none
        real x0, y0, xleng, xmax, xmin, sizex
        logical ticup, labtop, dopow
        integer llx
        character titlex*(*)
        logical doasis
        character title*80

        logical dosci
        real xl, xcen, tlen1, tlen2, ypos, ypos1, ypos2, xmn, xmx, xxv
        real xval, xnorm, dxx, dxxx, xe, xb, xc, xx, xv, xxx, ticlen
        real sizexx
        integer nxnorm, nscalx, n, ifac, ilow, iup, ndxxx, i, j, nshft
        integer ndec, lx
        real oneneg


        if(xmin.lt.0.0 .or. xmax.lt.0.0)then
          oneneg = 3.5*sizex
        else
          oneneg = 2.5*sizex
        endif
        lx = llx
        if(lx.lt.0)lx = 0
c-----
c       put in the axis
c-----
        xl = x0 + xleng
        xcen = x0 + 0.5*xleng
        call plot(x0,y0,3)
        call plot(xl,y0,2)
c-----
c       set up positioning for tics
c-----
        tlen1 = 0.6*sizex
        tlen2 = 1.2*sizex
        if(labtop)then
            if(ticup)then
                ypos = y0 + 2.0 * sizex
                ypos1 = ypos + 0.7*sizex
                ypos2 = y0 + 4.0 * sizex
            else
                ypos = y0 + sizex
                ypos1 = ypos + 0.7*sizex
                ypos2 = y0 + 3.0 * sizex
            endif
        else
            if(ticup)then
                ypos = y0 - 2.0*sizex
                ypos1 = ypos + 0.7*sizex
                ypos2 = y0 - 4.5*sizex
            else
                ypos = y0 - 3.0*sizex
                ypos1 = ypos + 0.7*sizex
                ypos2 = y0 - 5.5*sizex
            endif
        endif
c-----
c       compute limits for scale
c-----
        xmn = min(xmin,xmax)
        xmx = max(xmin,xmax)
        xxv = xmn
c-----
c       safety
c-----
        if(xmn.eq.xmx)xmx = xmn + 1.0
c-----
c       get some idea of the size of the number to be plotted
c-----
        xval = max(abs(xmn),abs(xmx))
c-----
c       we set the maximum number of decimal digits as 3, 
c       if the number is greater than 1000 or less than 0.01, use
c       scientific notation
c       get the maximum number of digits for the number
c-----
        if(xval.le.0.1 .and.xval.gt.0.0 )then
            dosci = .true.
            xnorm = alog10(1.0/xval)
            nxnorm = INT(xnorm+1)
            xnorm = 10.0**( nxnorm)
            nscalx = -nxnorm
        else if( xval.ge.10000.0)then
            dosci = .true.
            xnorm = alog10(xval)
            nxnorm = INT(xnorm)
            nscalx =  nxnorm
            xnorm = 10.0**(-nxnorm)
        else
            dosci = .false.
            xnorm = 1.0
            nxnorm = 0.0
            nscalx = 0
        endif
c-----
c       choose the maximum number of tics somehow on
c       xleng, and size of numbers and necessary decimals
c       dxx is the increment of the major tics while dxxx
c           is the increment of the minor tics
c           because this switches to scientific notation for
c           large and small numbers the range of dxx is
c           like 0.007 0.07 0.7 7 70 700  and for scientific
c           notation is 0.7. This to provide some order, we will
c           override the estimate to start with the larger of
c           
c-----
        dxx = xnorm*(XMX - XMN)/5.0
        if(dxx.eq.0.0)dxx = 1.0
        call redodvv(dxx)
c-----
c       special case for 0 to 360
c-----
        if(XMN.eq.0.0 .and. XMX.eq.360.0)then
              dxx = xnorm*(XMX - XMN)/4.0
        endif
c-----
c       get start limits for search - we will place at most 10 tics
c       here
c-----
        n = INT(alog10(dxx))
        if(dxx.lt.1.0)n = n -1
        ifac = INT(dxx*10.0**(-n))
        
c----- 
c       now get a new minimum value that is some multiple 
c       of ifac 10**n 
c----- 
        dxx = ifac * 10.0**(n)
        ilow = INT(xnorm*xmn/dxx -1)
        iup  = INT(xnorm*xmx/dxx +1)
        if(ifac.eq.1)then
            ndxxx = 5
        else if(ifac.eq.2)then
            ndxxx = 5
        else
            ndxxx = ifac
        endif
        dxxx = dxx/ndxxx

c-----
c       loop for labels
c-----
c       xe is the end of the previous label
c-----
        xc = xmn
        call plot(x0,y0,3)
c-----
c       here
c       xx is the position of the plot corresponding to a value of xxv
c       xb is the position of the beginning of the number ot be plotted
c       xe is the position of the end of the last number 
c           plotted + one space
c       xl is the end of the line
c-----
c-----
c       consider numbers between ilow*dxx and iup*dxx
c-----
        xe = x0
        do  i=ilow,iup
            do j=0,ndxxx-1
            
            xxv = i*dxx + j*dxxx
 
            xxv = xxv/xnorm
            if(xxv.ge.xmn .and. xxv .le.xmx)then
c-----
c               We can plot this number/tic
c-----
                xx = x0 + (xxv - xmn)*xleng/(xmx - xmn)
                call plot(xx,y0,3)
                if(j.eq.0 )then
c-----
c       estimate the number of digits for proper labeling
c       only do this once for the plot
c           dosci   = .false.    +0.000
c               = .true.     +0.000 10^+00
c       also Never overlay plot labels
c----
c       get width of number string
c-----
                    if(doasis)then
                        xv = xxv*xnorm
                    else
                        xv = - xxv*xnorm
                    endif
                    nshft = 1
                    if(xv .lt. 0.0)nshft=nshft+1
                    if(dosci)then
                        nshft = nshft + 1
                        ndec = 2
                        nshft = nshft +  ndec
                    else 
                        if(abs(dxx) .gt. 1.0)then
                            ndec = -1
                        else
                            if(ABS(dxx)+0.0001.lt.0.1)then
                                ndec = 3
                            else
                                ndec = 2
                            endif
                            nshft = nshft + ndec
                        endif
                    if(abs(xv) .ge. 10.0) nshft=nshft+1
                    if(abs(xv) .ge. 100.0)nshft=nshft+1
                    if(abs(xv) .ge. 1000.0)nshft=nshft+1
                    endif
                    xxx = 0.5*(nshft -0.5)*sizex
                    xb = xx - xxx
                    if((xb - xe) .gt. oneneg)then
                       if(xb.ge.x0.and.(xx+xxx).le.xl)then
                          if(dopow)then
                             call number(xb,ypos,sizex,xv,0.0,ndec)
                          endif
                          xe = xb + xxx + sizex
                          endif
                         ticlen = tlen2
                    else
                         ticlen = 0.8*tlen2
                    endif
                else
                    ticlen = tlen1
                endif
                call plot(xx,y0,3)
                if(ticup)then   
                    call plot(xx,y0+ticlen,2)
                else
                    call plot(xx,y0-ticlen,2)
                endif
                call plot(xx,y0,3)
            endif
            enddo
        enddo
c-----
c       put in the title if we put in the numbers
c-----
        if(dopow )then
            sizexx = 1.2*sizex
            title = ' '
            if(lx.gt.0)then
                title(1:lx) = titlex(1:lx)
            else
                lx = 0
            endif
            if(dosci)then
                title(lx +1:lx+5)=' *10 '
                call gcent(xcen,ypos2,sizexx,title,0.0)
                xe = xcen + 0.5*sizexx*(lx+4)
                call number(xe,ypos2+0.7*sizexx,0.7*sizexx,
     1              real(nscalx),0.0,-1)
            else
                call gcent(xcen,ypos2,sizexx,title,0.0)
            endif
        endif
        return
        end

        subroutine dnliny(x0,y0,yleng,ymax,ymin,
     1      sizey,ticlft,lablft,dopow,lly,titley)
        implicit none
c-----
c       put in tic marks along the y-axis but with
c       maximum at the bottom and the minimum to the top
c       this is useful for plotting depths
c
c       x0  R - position of bottom side of axis
c       y0  R - position of bottom side of axis
c       yleng   R - length of Y axis
c       ymax    R - maximum value of number corresponding to
c                   top
c       ymax    R - maximum value of number corresponding to
c                   bottom
c       sizey   R - size of numbers, e.g., 10
c       ticlft  L   - .true. tics go left
c                 .false. tics go right
c       lablft  L   - .true. numbers goes to left of axis
c                 .false. goes to right
c       dopow   L   - .true. put in numbers labels
c       ly  I length of y-axis title string
c       titley  Ch  y-axis title string
c-----
        real x0, y0, yleng, ymax, ymin, sizey
        logical ticlft, lablft, dopow
        integer lly
        character titley*(*)
        call dolny(x0,y0,yleng,ymax,ymin,
     1      sizey,ticlft,lablft,dopow,lly,titley,.false.)
        return
        end

        subroutine doliny(x0,y0,yleng,ymax,ymin,
     1      sizey,ticlft,lablft,dopow,lly,titley)
        implicit none
c-----
c       put in tic marks along the y-axis but with
c       maximum at the top and the minimum at the bottom
c       this is the usual type of plot
c
c       x0  R - position of bottom side of axis
c       y0  R - position of bottom side of axis
c       yleng   R - length of Y axis
c       ymax    R - maximum value of number corresponding to
c                   top
c       ymax    R - maximum value of number corresponding to
c                   bottom
c       sizey   R - size of numbers, e.g., 10
c       ticlft  L   - .true. tics go left
c                 .false. tics go right
c       lablft  L   - .true. numbers goes to left of axis
c                 .false. goes to right
c       dopow   L   - .true. put in numbers labels
c       ly  I length of y-axis title string
c       titley  Ch  y-axis title string
c-----
        real x0, y0, yleng, ymax, ymin, sizey
        logical ticlft, lablft, dopow
        integer lly
        character titley*(*)
        call dolny(x0,y0,yleng,ymax,ymin,
     1      sizey,ticlft,lablft,dopow,lly,titley,.true.)
        return
        end

        subroutine dolny(x0,y0,yleng,ymax,ymin,
     1      sizey,ticlft,lablft,dopow,lly,titley,doasis)
c-----
c       put in tic marks along the y-axis 
c       this routine does all of the work
c
c       x0  R - position of bottom side of axis
c       y0  R - position of bottom side of axis
c       yleng   R - length of Y axis
c       ymax    R - maximum value of number corresponding to
c                   top
c       ymax    R - maximum value of number corresponding to
c                   bottom
c       sizey   R - size of numbers, e.g., 10
c       ticlft  L   - .true. tics go left
c                 .false. tics go right
c       lablft  L   - .true. numbers goes to left of axis
c                 .false. goes to right
c       dopow   L   - .true. put in numbers labels
c       ly  I length of y-axis title string
c       titley  Ch  y-axis title string
c       doasis  L   .true. plot values
c               .false. annotate with negative (useful for depth plot)
c-----
        implicit none
        real x0, y0, yleng, ymax, ymin, sizey
        logical ticlft, lablft, dopow
        integer lly
        character titley*(*)
        logical doasis


        real ymn, ymx, ynorm, dyy, yv
        real dyyy, ye, yyv, yp, yb, ticlen, xpos2, sizeyy, yy, yval
        real yl, ycen, tlen1, tlen2, xpos
        integer ly
        integer ifac, ndyyy, ilow, iup, i, j, n
        integer  nynorm, nscaly, maxdig, ndec
        character cmd*10
        character title*80
        logical dosci

        ly = lly
        if(ly.lt.0)ly = 0
c-----
c       put in the axis
c-----
        yl = y0 + yleng
        ycen = y0 + 0.5*yleng
        call plot(x0,y0,3) 
        call plot(x0,yl,2) 
c-----
c       set up positioning for tics
c-----
        tlen1 = 0.6*sizey
        tlen2 = 1.2*sizey
        if(lablft)then
            if(ticlft)then
                xpos   = x0 - 2*sizey
            else
                xpos   = x0 - sizey
            endif
            cmd = 'RIGHT'
        else
            if(ticlft)then
                xpos   = x0 + 2.00*sizey
            else
                xpos   = x0 + 3.0*sizey
            endif
            cmd = 'RIGHT'
        endif
c-----
c       compute limits for scale
c-----
        ymn = min(ymin,ymax)
        ymx = max(ymin,ymax)
c-----
c       safety
c-----
        if(ymn.eq.ymx)ymx = ymn + 1.0
c-----
c       get some idea of the size of the number to be plotted
c-----
        yval = max(abs(ymn),abs(ymx))
c-----
c     get some idea of size of the number
c     maxdig is the number of digits on axis
c        this will include the decimal point
c        later if the one of ymin, ymax
c        is negative, then maxdig++ for the
c        negative sign
c
c        dosci .true. or .false. to indicate if
c            scientic notation is to be used
c            in which case nscaly is the power of 10
c            and then the number is multiplied by ymorm
c            to make something between 1.00 and 9.999
c        the maxdig is very important is the number labels are
c        to the right of the axis since a left justificaiton would not
c        look good.  This would occur when dopow=.true. and lablft=.false.
c-----
        if(yval.le.0.1 .and.yval.gt.0.0 )then
            dosci = .true.
            ynorm = alog10(1.0/yval)
            nynorm = INT(ynorm+1)
            ynorm = 10.0**( nynorm)
            nscaly = -nynorm
            maxdig = 4
            ndec = 2
        else if( yval.ge.10000.0)then
            dosci = .true.
            ynorm = alog10(yval)
            nynorm = INT(ynorm)
            ynorm = 10.0**(-nynorm)
            nscaly =  nynorm
            maxdig = 4
            ndec = 2
        else
c-----     
c           make it look good for 
c           yval= 0.1000x to yval=9999.0
c----
            dosci = .false.
            ynorm = 1.0
            nynorm = 0
            nscaly = nynorm
            if(yval.ge.0.1 .and. yval.lt.1.0)then
                 maxdig=4
                 ndec = 2
            else if(yval.ge.1.0 .and. yval.lt.10.)then
                 maxdig=4
                 ndec = 2
            else if(yval.ge.10.0 .and. yval.lt.100.)then
                 maxdig=5
                 ndec = 2
            else if(yval.ge.100.0 .and. yval.lt.1000.)then
                 maxdig=5
                 ndec = 1
            else if(yval.ge.1000.0 .and. yval.lt.10000.)then
                 maxdig=4
                 ndec = -1
            endif
        endif
c-----
c     make room for the negative sign
c-----
        if(ymn.lt.0.0 .or. ymx.lt.0.0)then
            maxdig = maxdig + 1
        endif
        if(.not.lablft)then
            xpos = xpos + maxdig*sizey
        endif
c-----
c       choose the maximum number of tics somehow on
c       yleng, and size of numbers and necessary decimals
c-----
        dyy = ynorm*(ymx - ymn)/5.0
        if(dyy.eq.0.0)dyy = 1.0
       
        call redodvv(dyy)
c-----
c       get start limits for search - we will place at most 10 tics
c       between major intervals
c-----
        n = INT(alog10(dyy))
        if(dyy.lt.1.0)n = n -1
        ifac = INT(dyy*10.0**(-n))
        
c----- 
c       now get a new minimum value that is some multiple 
c       of ifac 10**n 
c
c       dyy  = major tic increment
c       dyyy = minor tic incement
c----- 
        dyy = ifac * 10.0**(n)
c-----
c       set major tic increment
c-----
        ilow = INT(ynorm*ymn/dyy -1)
        iup  = INT(ynorm*ymx/dyy +1)
        if(ifac.eq.1)then
            ndyyy = 5
        else if(ifac.eq.2)then
            ndyyy = 5
        else
            ndyyy = ifac
        endif
        dyyy = dyy/ndyyy
c-----
c       loop for ticmarks and numbers
c-----
c       ye is the end of the previous label
c-----
        call plot(x0,y0,3)
c-----
c       here
c       yy is the position of the plot corresponding to a value of yyv
c       yend is the position of the end of the last number plotted
c-----
        ye = y0 - 2.*sizey
c-----
c       consider numbers between ilow*dyy and iup*dyy
c-----
        do i=ilow,iup
           do j=0,ndyyy-1
c-----
c          yyv is the value to be plotted
c-----
              yyv = i*dyy + j*dyyy
c-----
c          now scale this to relate to original number
c----
           
              yyv = yyv/ynorm
              if(yyv.ge.ymn .and. yyv .le.ymx)then
                 yy = y0 + (yyv - ymn)*yleng/(ymx - ymn)
                 call plot(x0,yy,3)
                 if(j.eq.0)then
                    if(doasis )then
                       yv = yyv * ynorm
                    else
                       yv = - yyv * ynorm
                    endif
                    yp = yy - 0.5*sizey
                    yb = yy
c-----
c     plot the label as lon as it does not overlap with label above or
c     below
c     also see if  label fits within the y0,yl range
c
c     currently if a numebr is to be at the top o bottom, plot the
c     number but make a slight vertical adjustment so that the number
c     does not extend ouside [y0,yl]
c   
c     if you do not want this behavior, then uncommend the three CRBH lines
c     and commend the following
c                   if(yy-min(y0,yl).lt.0.5*sizey)yp = yy
c                   if(max(y0,yl)-yy.lt.0.5*sizey)yp = yy-sizey
c-----
                    if((yb - ye ) .gt. 2.0*sizey)then
CRBH                   if(yy.gt.min(y0,yl)+0.5*sizey 
CRBH 1                    .and. yy.lt.max(y0,yl)-0.5*sizey)then
                    if(yy-min(y0,yl).lt.0.5*sizey)yp = yy
                    if(max(y0,yl)-yy.lt.0.5*sizey)yp = yy-sizey
                          if(dopow)then
                             call mynum(xpos,yp,sizey,yv,0.0,ndec,cmd)
                          endif
                          ye = yb 
CRBH                   endif
                       ticlen = tlen2
                    else
                       ticlen = 0.8*tlen2
                    endif
c-----
c                       large tic mark length set
c-----
                 else 
c-----
c                       small tic mark
c-----
                        ticlen = tlen1
                 endif
                 call plot(x0,yy,3)
                 if(ticlft)then  
                    call plot(x0-ticlen,yy,2)
                 else
                    call plot(x0+ticlen,yy,2)
                 endif
              endif
           enddo
        enddo
c-----
c       put in the title if we put in the numbers
c       be careful about positioning when the labels are
c       to the right
c-----
        if(lablft)then
            if(ticlft)then
                xpos2  = x0 - maxdig* sizey -3.9*sizey
            else
                xpos2  = x0 - maxdig* sizey -2.9*sizey
            endif
        else
            if(ticlft)then
                xpos2  = x0 + maxdig* sizey +3.0*sizey
            else
                xpos2  = x0 + maxdig* sizey +4.0*sizey
            endif
        endif
        if(dopow .and. ly .gt.0 )then
            sizeyy = 1.2*sizey
            title = ' '
            if(ly.gt.0)then
                title = titley(1:ly)
            else
                ly = 0
            endif
            if(dosci)then
                title(ly +1:ly+5)=' *10 '
                if(lablft)then
                call gcent(xpos2,ycen,sizeyy,title, 90.0)
                    ye = ycen + 0.5*sizeyy*(ly+4)
                call number(xpos2-0.7*sizey,ye,0.7*sizeyy,
     1              real(nscaly), 90.0,-1)
                else
                call gcent(xpos2,ycen,sizeyy,title,-90.0)
                    ye = ycen - 0.5*sizeyy*(ly+4)
                call number(xpos2+0.7*sizeyy,ye,0.7*sizeyy,
     1              real(nscaly),-90.0,-1)
                endif
            else
                if(lablft)then
                call gcent(xpos2,ycen,sizeyy,title, 90.0)
                else
                call gcent(xpos2,ycen,sizeyy,title,-90.0)
                endif
            endif
        endif
        return
        end


        subroutine fillit(cmd,rad,x0,y0)
c-----
c       plot a filled symbol at the current color
c
c       cmd    C - string indicating the symbol
c                  'SQ' square
c                  'TR' triangle
c                  'IT' inverted triangle
c                  'HX' hexagon
c                  'DI' diamond
c                  'CI' circle
c       red    R - bounding circle radius
c       x0     R - coordinates of center of symbol
c       y0     R -
c-----
        implicit none
        character cmd*(*)
        real rad, x0, y0

        real xval(370), yval(370)
        integer ipatx, ipaty, jj, i
        real xlen, ylen, r2
        real x1,y1,x2,y2,x3,y3, ang,sa,ca,xnew,ynew
c-----
c       fill in a solid symbol
c-----
        ipatx = 0
        ipaty = 0
        xlen = 0.01
        ylen = 0.01
        r2 = rad 
        if(cmd(1:2).eq.'SQ')then
            jj = 0
            do i=45,405,90
                ang = 3.1415927*i/180.0
                sa = sin(ang)
                ca = cos(ang)
                xnew = x0 + r2*ca
                ynew = y0 + r2*sa
                jj = jj + 1
                xval(jj) = xnew
                yval(jj) = ynew
            enddo
            call shadep(4,xval,yval)
        else if(cmd(1:2).eq.'TR')then
            jj = 0
            do i=90,450,120
                ang = 3.1415927*i/180.0
                sa = sin(ang)
                ca = cos(ang)
                xnew = x0 + r2*ca
                ynew = y0 + r2*sa
                jj = jj + 1
                xval(jj) = xnew
                yval(jj) = ynew
            enddo
            call shadep(3,xval,yval)
        else if(cmd(1:2).eq.'IT')then
            jj = 0
            do  i=270,630,120
                ang = 3.1415927*i/180.0
                sa = sin(ang)
                ca = cos(ang)
                xnew = x0 + r2*ca
                ynew = y0 + r2*sa
                jj = jj + 1
                xval(jj) = xnew
                yval(jj) = ynew
            enddo
            call shadep(3,xval,yval)
        else if(cmd(1:2).eq.'HX')then
            x1 = x0 - 0.866*r2
            x2 = x0 + 0.866*r2
            x3 = x0
            y1 = y0 - 0.500*r2
            y2 = y0 - 0.500*r2
            y3 = y0 + r2
            call shadet(x1,y1,x2,y2,x3,y3,
     1          ipatx,ipaty,xlen,ylen)
            y1 = y0 + 0.500*r2
            y2 = y0 + 0.500*r2
            y3 = y0 - r2
            call shadet(x1,y1,x2,y2,x3,y3,
     1          ipatx,ipaty,xlen,ylen)
        else if(cmd(1:2).eq.'DI')then
            jj = 0
            do i=0,360,90
                ang = 3.1415927*i/180.0
                sa = sin(ang)
                ca = cos(ang)
                xnew = x0 + r2*ca
                ynew = y0 + r2*sa
                jj = jj + 1
                xval(jj) = xnew
                yval(jj) = ynew
            enddo
            call shadep(4,xval,yval)
        else if(cmd(1:2).eq.'CI')then
            jj = 0
            do i=0,360,10
                ang = 3.1415927*i/180.0
                sa = sin(ang)
                ca = cos(ang)
                xnew = x0 + r2*ca
                ynew = y0 + r2*sa
                jj = jj + 1
                xval(jj) = xnew
                yval(jj) = ynew
            enddo
            call shadep(jj,xval,yval)
        endif
        return
        end

        subroutine curvit(cmd,rad,x0,y0)
c-----
c       plot a symbol outline at the current color
c
c       cmd    C - string indicating the symbol
c                  'SQ' square
c                  'TR' triangle
c                  'IT' inverted triangle
c                  'HX' hexagon    commented out
c                  'DI' diamond
c                  'CI' circle
c       red    R - bounding circle radius
c       x0     R - coordinates of center of symbol
c       y0     R -
c-----
        implicit none
        character cmd*(*)
        real rad, x0, y0

        real xval(370), yval(370)
        integer ipatx, ipaty, jj, i
        real xlen, ylen, r2, ang, ca, sa, xnew, ynew
c-----
c       plot an outline
c-----
        ipatx = 0
        ipaty = 0
        xlen = 0.01
        ylen = 0.01
        r2 = rad 
        if(cmd(1:2).eq.'SQ')then
            jj = 0
            do  i=45,405,90
                ang = 3.1415927*i/180.0
                sa = sin(ang)
                ca = cos(ang)
                xnew = x0 + r2*ca
                ynew = y0 + r2*sa
                jj = jj + 1
                xval(jj) = xnew
                yval(jj) = ynew
            enddo
            call drawcv(jj,xval,yval)
        else if(cmd(1:2).eq.'TR')then
            jj = 0
            do  i=90,450,120
                ang = 3.1415927*i/180.0
                sa = sin(ang)
                ca = cos(ang)
                xnew = x0 + r2*ca
                ynew = y0 + r2*sa
                jj = jj + 1
                xval(jj) = xnew
                yval(jj) = ynew
            enddo
            call drawcv(jj,xval,yval)
        else if(cmd(1:2).eq.'IT')then
            jj = 0
            do  i=270,630,120
                ang = 3.1415927*i/180.0
                sa = sin(ang)
                ca = cos(ang)
                xnew = x0 + r2*ca
                ynew = y0 + r2*sa
                jj = jj + 1
                xval(jj) = xnew
                yval(jj) = ynew
            enddo
            call drawcv(jj,xval,yval)
C       else if(cmd(1:2).eq.'HX')then
C           x1 = x0 - 0.866*r2
C           x2 = x0 + 0.866*r2
C           x3 = x0
C           y1 = y0 - 0.500*r2
C           y2 = y0 - 0.500*r2
C           y3 = y0 + r2
C           call shadet(x1,y1,x2,y2,x3,y3,
C     1         ipatx,ipaty,xlen,ylen)
C           y1 = y0 + 0.500*r2
C           y2 = y0 + 0.500*r2
C           y3 = y0 - r2
C           call shadet(x1,y1,x2,y2,x3,y3,
C     1         ipatx,ipaty,xlen,ylen)
        else if(cmd(1:2).eq.'DI')then
            jj = 0
            do  i=0,360,90
                ang = 3.1415927*i/180.0
                sa = sin(ang)
                ca = cos(ang)
                xnew = x0 + r2*ca
                ynew = y0 + r2*sa
                jj = jj + 1
                xval(jj) = xnew
                yval(jj) = ynew
            enddo
            call drawcv(jj,xval,yval)
        else if(cmd(1:2).eq.'CI')then
            jj = 0
            do  i=0,360,10
                ang = 3.1415927*i/180.0
                sa = sin(ang)
                ca = cos(ang)
                xnew = x0 + r2*ca
                ynew = y0 + r2*sa
                jj = jj + 1
                xval(jj) = xnew
                yval(jj) = ynew
            enddo
            call drawcv(jj,xval,yval)
        endif
        return
        end

        subroutine drawcv(jj,xval,yval)
c-----
c       draw a curve as lienar segments
c       
c       jj   I - number of points
c       xval R - x-coordinate of point
c       yval R - y-coordinate of point
c-----
        implicit none
        integer jj
        real xval(jj),yval(jj)

        integer i

        call plot(xval(1), yval(1), 3)
        do  i=1,jj
            call plot(xval(i), yval(i), 2)
        enddo
        call plot(xval(jj), yval(jj), 3)
        return
        end
        

        subroutine rclip(xmin,xmax,ymin,ymax,xc1,yc1,xc2,yc2,
     1          x0,yy0,x1,yy1,iplt)
        implicit none
c-----
c       perform a line clipping for the plot region
c       given desire to connect (x0,yy0) -> (x1,yy1)
c       the clipped coordinates are (xc1,yc1) -> (xc2,yc2)
c
c       note the low level CALPLOT has a clipping routine
c-----
c       xmin    R - minimum value of X
c       xmax    R - maximum value of X
c       ymin    R - minimum value of Y
c       ymax    R - maximum value of Y
c       xc1     R - coordinate #1 of plotted line segment
c       yc1     R
c       xc2     R - coordinate #2 of plotted line segment
c       yc2     R
c       x0      R - first coordinate
c       yy0     R
c       x1      R - second coordinate
c       yy1     R
c       iplt    I - > 0 plot the segment, otherwise do not
c-----
        real xmin, xmax, ymin,ymax
        real xc1, yc1, xc2, yc2
        real x0, yy0, x1, yy1
        integer iplt

c-----
c       Fortran implementation of Cohen and Sutherland Line 
c           Clipping Routine
c-----
        real xx0, yz0, xx1, yz1, x, y, slope, slopeinv
        integer c0, c1, c

        iplt = 0
        xx0 = x0
        yz0 = yy0
        xx1 = x1
        yz1 = yy1
        call linecode(xx0,yz0,ymax,xmin,ymin,xmax,c0)
        call linecode(xx1,yz1,ymax,xmin,ymin,xmax,c1)
        if( xx0 .ne. xx1)slope = (yz1-yz0)/(xx1-xx0)
        if( yz0 .ne. yz1) slopeinv = (xx1-xx0)/(yz1-yz0)
 1000   continue
        if(c0 .eq. 0 .and. c1 .eq.0)go to 1001
            if( (mod(c0,2).eq.1 .and. mod(c1,2).eq.1).or.
     1      (mod(c0/2,2).eq.1 .and. mod(c1/2,2).eq.1).or.
     2      (mod(c0/4,2).eq.1 .and. mod(c1/4,2).eq.1).or.
     3      (mod(c0/8,2).eq.1 .and. mod(c1/8,2).eq.1))return

                if(c0 .eq. c1) return 
                if(c0 .eq. 0 )then
                        c = c1
                else 
                        c = c0
            endif
                if(mod(c,2).eq.1)then
                        y = yz0 + slope*(xmin-xx0)
                        x = xmin
                endif
                if(mod(c/2,2).eq.1)then 
                y = yz0 + slope*(xmax-xx0)
                        x = xmax
                endif
                if(mod(c/8,2).eq.1)then
                        x = slopeinv*(ymax-yz0) + xx0
                        y = ymax
                endif
                if(mod(c/4,2).eq.1)then
                        x = slopeinv*(ymin-yz0) + xx0
                        y = ymin
                endif
                if(c .eq. c0)then
                        xx0 = x
                        yz0 = y
                call linecode(xx0,yz0,ymax,xmin,ymin,xmax,c0)
                else 
                        xx1 = x
                        yz1 = y
                call linecode(xx1,yz1,ymax,xmin,ymin,xmax,c1)
                endif
        go to 1000
 1001   continue
        xc1 = xx0
        xc2 = xx1
        yc1 = yz0
        yc2 = yz1
        iplt = 1
        return
        end
c       
        subroutine linecode(x,y,ymax,xmin,ymin,xmax,c)
c-----
c       clipping algorithm
c
c       x      R - coordinate of point
c       y      R -
c       xmin   R - bounding box for clip region
c       ymin   R -
c       xmax   R -
c       ymax   R -
c       c      I - clipping code
c-----
        implicit none
        real x, y, ymax, xmin, ymin, xmax
        integer c
        c = 0

        if(x .lt. xmin) then
            c = 1
        else if(x .gt. xmax)then
            c = 2
        endif
        if(y .lt. ymin) then
            c = c + 4
        else if(y.gt.ymax)then
            c = c + 8
        endif
        return
        end
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: NUMBER                                                c
c                                                                     c
c      COPYRIGHT (C)  1986, 1994 R. B. Herrmann                       c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
      subroutine mynum (xpage,ypage,height,fpn,angle,ndec,cmd)
      implicit none
c-----
c.....     xpage,ypage coordinates of lower left corner of number.
c.....     height   height of plotted number.
c.....     fpn      floating point number to be plotted.
c.....     angle    angle at which number is plotted, in degrees.
c.....     ndec     number of decimal places to be drawn.
c.....     this version of number requires the symbol version with
c.....     999. x, y feature, and  nc = 0 feature.
c.....
c.....     ndec .ge. 1000 .le. 1009 invokes an exponential format
c.....     of the Fortran Format Ew.d, where d = ndec - 1000
c.....     and w = 7 + d. The output will look like
c.....     sd.dddddEsdd, where s is the sign and d is a digit
c.....     In all cases, no more than 9 significant figures are
c.....     permitted to the right of the decimal place.
c.....
c.....     ndec .ge. 2000 .le. 2009 invokes exponential but
c.....     plotted in scientific notation, e.g., 1.0 10 ^ exp
c-----
        real xpage, ypage, height, fpn, angle
        integer ndec
        character cmd*(*)

        common/Scplot/x0,y0
        real x0, y0
        common/Xcplot/xold,yold,xcur,ycur
        real xold,yold,xcur,ycur

        character num*20
        character cexp*4
        real fp, ang, ca, sa, xl, xx0, xxx0, yl, yy0, yyy0, ht, aex
        integer n, ii, mdec, nzflag, jex, iex, nexp, nman

        if(xpage .lt. 999.0) then
            x0 = xpage
        endif
        if(ypage .lt. 999.0) then
            y0 = ypage
        endif
c-----
c    to avoid dregs in memory null first character
c-----
      num(1:1)= ' '
      ii=0
c-----
      n = 0
        fp = fpn
        if(fpn.lt.0.0)then
                n=n+1
                num(n:n)='-'
                fp=abs(fp)
        endif
      if(ndec.lt.1000)then
                nzflag = 0
                mdec = -ndec
c               if(ndec.lt.0)then
c                       mdec = mdec - 1
c               endif
                if(ndec.eq.-1)then
                        nzflag = 1
                        mdec = 1
                endif
                call ffpack(num,n,fp,mdec,19,nzflag)
        else if(ndec.ge.1000)then
                       mdec = 1
                if(ndec .ge.1000 .and. ndec .le. 1009)then
                       mdec = ndec -1000
                else if(ndec .ge.2000 .and. ndec .le. 2009)then
                       mdec = ndec -2000
                else 
                       mdec = 9
                endif
                mdec = - mdec
                if(fp.eq.0.0)then
                        iex=0
                else
                        aex = alog10(fp)
                        iex = INT(aex)
c----------careful check of integer exponent
                        if(iex.lt.0.and.float(iex).gt.aex)iex=iex-1
                endif
                fp = fp * 10.** (- iex)
                nzflag = 0
                call ffpack(num,n,fp,mdec,14,nzflag)
c---put in exponent assuming iex < 100
                nman=n
                n=n+1
                num(n:n)='E'
                if(iex.lt.0)then
                        n=n+1
                        num(n:n)='-'
                        iex = -iex
                else
                        n=n+1
                        num(n:n)='+'
                endif
                nexp = n
                n=n+1
                jex = iex
                jex = mod(jex/10,10)
                num(n:n) = char(ichar('0' )+ jex)
                jex = mod(iex,10)
                n = n + 1
                num(n:n) = char(ichar('0' )+ jex)
        endif
        if(ndec .le. 1009)then
        if(cmd.eq.'LEFT')then
                call symbol(x0,y0,height,num,angle,n)
        else if(cmd.eq.'CENTER')then
                call symbol(x0-0.5*n*height,y0,height,num,
     1               angle,n)
        else if(cmd.eq.'RIGHT')then
                call symbol(x0-n*height,y0,height,num,angle,n)
        endif
        else if(ndec .ge. 2000 .and. ndec.le.2009)then
c-----
c       save the exponent, replace the E-EX with 'x10'
c-----
                cexp(1:3) = num(nexp:n)
                num(nman+1:nman+3)='x10'
            xx0 = x0
            yy0 = y0
                call symbol(x0,y0,height,num(1:nman+3),
     1               angle,nman+3)
c-----
c       get the proper position because of rotation
c-----
                ang = angle*3.1415927/180.0
                ca = cos(ang)
                sa = sin(ang)
                xl = height*(nman+3)
                yl = 0.6*height
            xx0 = xx0 + xl*ca
            yy0 = yy0 + xl*sa
                xxx0 = xx0  - yl*sa
                yyy0 = yy0  + yl*ca
                ht = 0.7*height
                call symbol(xxx0,yyy0,ht,cexp(1:3),angle,3)
c-----
c           now position at end of proper string
c-----
            xl =  3.0*ht
            xcur = xx0 + xl*ca
            ycur = yy0 + xl*sa
            x0 = xcur
            y0 = ycur
            
c-----
        endif
        return
        end

        subroutine ffpack(num,n,fpv,mdec,mwid,nzflag)
        implicit none
        character num*(*)
        integer n, mdec, mwid, nzflag
        real fpv

        real fp
        integer nzflg, m, maxm, iex, iup, ilw, k, ipos, mm


        fp = fpv
        nzflg = nzflag
c-----since we have a maximum field width of mwid and since
c-----the . takes a position, we can only have mwid -1 figures
      m = mdec 
        maxm = 9
        if(m.gt.maxm)m=maxm
        if(m.lt. -maxm)m= -maxm
        mm = m
        if(m.gt.0) mm = mm - 1
c----do a simple rounding up
        fp = fp + 0.5 * 10.** (mm)
c---- get position of first significant figure
c
c     5 4 3 2 1 0 -1 -2 -3
c
c     d d d d d .  d  d  d
c
c----

        iex = INT(alog10(fp) )
        if(fp.ge.1.0)then
                iex = iex + 1
        else
                iex = iex - 1
        endif
c----
c----   procede from most significant down
c----   but never less than 10**-9
c----
        
        if(fp.le.1.0e-9)then
                fp = 0.0
                iex = -9
        endif
        ilw = mdec
        if(fp.ge.1.0)then
                iup=iex
        else
                iup = 1
        endif
        if(m.le.0)then
                ilw= m+1
        else
                ilw = m 
        endif
        if(iex.gt.9)then
                ilw = 1
                nzflg = 1
        endif
        do ipos=iup,ilw,-1
                k = INT(fp * 10.**( -ipos +1 ))
                if(n.lt.mwid)then
                        n = n + 1
                        num(n:n) = char(ichar('0')+k)
                endif
                fp = fp - (float(k) * 10.**(ipos-1))
        if(ipos.eq.1.and.nzflg.eq.0.and.n.lt.mwid)then
                        n = n + 1
                        num(n:n) = '.'
                endif
        enddo
        return
        end

      subroutine msupsc (xpage,ypage,height,str,npow,angle)
      implicit none
c-----
c        routine to plot subscript or superscript
c-----
c          xpage   R - cooreinate os lower left corner of number
c          ypage   R - 
c          height  R - height of plotted number.
c          str     C - string
c          npow    I - power
c          angle   R - angle at which number is plotted, in degrees.
c-----
        real xpage, ypage, height, angle
        integer npow
        character str*(*)
        common/Scplot/x0,y0
        real x0, y0
        common/Xcplot/xold,yold,xcur,ycur
        real xold, yold, xcur, ycur

        integer lgstr

        real xx0, yy0, ang, sa, ca, xl, yl, xxx0, yyy0, ht
        integer il

        if(xpage .lt. 999.0) then
            x0 = xpage
        endif
        if(ypage .lt. 999.0) then
            y0 = ypage
        endif
        xx0 = x0
        yy0 = y0
        il = lgstr(str)
        call symbol(x0,y0,height,str(1:il), angle,il)
c-----
c       get the proper position because of rotation
c-----
                ang = angle*3.1415927/180.0
                ca = cos(ang)
                sa = sin(ang)
                xl = height*(il)
                yl = 0.6*height
            xx0 = xx0 + xl*ca
            yy0 = yy0 + xl*sa
                xxx0 = xx0  - yl*sa
                yyy0 = yy0  + yl*ca
                ht = 0.7*height
                call number(xxx0,yyy0,ht,real(npow),angle,-1)
c-----
c           now position at end of proper string
c-----
            xl =  3.0*ht
            xcur = xx0 + xl*ca
            ycur = yy0 + xl*sa
            x0 = xcur
            y0 = ycur
            
c-----
        return
        end

        subroutine gsubsc(x,y,ht,s1,n1,s2,n2)
c-----
c       plot a string with a subscript
c-----
c       x     R -
c       y     R -
c       ht    R -
c       s1    C - string
c       n1    I - length of string
c       s2    C - subscript
c       n2    I - length of subscript
c                 the subscript is dropped by 0.4 ht and
c                 has its height of 0.7 ht
c       if n1 or n2 are negative, the respective string 
c       is plotted in a Greek font
c-----
        implicit none
        real x, y, ht
        integer n1, n2
        character s1*(*), s2*(*)
        common/savfon/infont
        integer infont

        integer n, nn

        if(n1.lt.0)then
            call gfont(4)
            nn = -n1
        else
            nn = n1
            call gfont(infont)
        endif
        if(nn.gt.0)call symbol(x,y,ht,s1,0.0,nn)
        if(n2.lt.0)then
            call gfont(4)
            n = -n2
        else
            n = n2
            call gfont(infont)
        endif
        if(n.gt.0)call symbol(x+nn*ht+0.4*ht,y-0.4*ht,0.7*ht,s2,0.0,n)
        if(n1.lt.0 .or. n2.lt.0)call gfont(infont)
        return
        end
        
        subroutine gsupsc(x,y,ht,s1,n1,s2,n2)
c-----
c       plot a string with a superscript
c-----
c       x     R -
c       y     R -
c       ht    R -
c       s1    C - string
c       n1    I - length of string
c       s2    C - superscript
c       n2    I - length of superscript
c                 the superscript is dropped by 0.4 ht and
c                 has its height of 0.7 ht
c       if n1 or n2 are negative, the respective string 
c       is plotted in a Greek font
c-----
        implicit none
        real x, y, ht
        integer n1, n2
        character s1*(*), s2*(*)
        common/savfon/infont
        integer infont

        integer n, nn
        if(n1.lt.0)then
            call gfont(4)
            nn = -n1
        else
            nn = n1
            call gfont(infont)
        endif
        if(nn.gt.0)call symbol(x,y,ht,s1,0.0,nn)
        if(n2.lt.0)then
            call gfont(4)
            n = -n2
        else
            n = n2
            call gfont(infont)
        endif
        if(n.gt.0)call symbol(x+nn*ht+0.4*ht,y+0.4*ht,0.7*ht,s2,0.0,n)
        if(n1.lt.0 .or. n2.lt.0)call gfont(infont)
        return
        end

        subroutine gsubsup(x,y,ht,s1,n1,s2,n2,s3,n3)
c-----
c       plot a string with a subscript and superscript
c-----
c       x     R -
c       y     R -
c       ht    R -
c       s1    C - string
c       n1    I - length of string
c       s2    C - subscript
c       n2    I - length of subscript
c                 the subscript is dropped by 0.4 ht and
c                 has its height of 0.7 ht
c       s3    C - superscript
c       n3    I - length of superscript
c                 the superscript is dropped by 0.4 ht and
c                 has its height of 0.7 ht
c       if n1, n2 or n3 are negative, the respective string 
c       is plotted in a Greek font
c-----
        implicit none
        real x, y, ht
        integer n1,n2,n3
        character s1*(*), s2*(*), s3*(*)
        common/savfon/infont
        integer infont

        integer n, nn
        if(n1.lt.0)then
            call gfont(4)
            nn = -n1
        else
            nn = n1
            call gfont(infont)
        endif
        if(nn.gt.0)call symbol(x,y,ht,s1,0.0,nn)
        if(n2.lt.0)then
            call gfont(4)
            n = -n2
        else
            n = n2
            call gfont(infont)
        endif
        if(n.gt.0)call symbol(x+nn*ht+0.4*ht,y-0.4*ht,0.7*ht,s2,0.0,n)
        if(n3.lt.0)then
            call gfont(4)
            n = -n3
        else
            n = n3
            call gfont(infont)
        endif
        if(n.gt.0)call symbol(x+nn*ht+0.4*ht,y+0.4*ht,0.7*ht,s3,0.0,n)
        if(n1.lt.0 .or. n2.lt.0 .or. n3.lt.0)call gfont(infont)
        return
        end

        subroutine redodvv(dvv)
c-----
c       select a nice increment for linear plots
c       this is used by the modified dolnx and dolny
c       
c       dvv    R - input increment
c                  which is modified to give a
c                  nicer increment
c-----
        implicit none
        real dvv
        real odvv
        real oodvv
        integer n
c------
c     the dxx is based on the actual limits
c     here override that value to make nicer major tic intervales
c
c     one entry here assume that the dvv has a value between
c     0.001 and 1000
c------
      oodvv = dvv
c-----
c     set a scale to make the choice simpler
c-----
      if(dvv.ge.1000.0)then
         n = -3 
      else if(dvv.ge.100.0   .and.dvv.lt.1000.0   )then
         n = -2 
      else if(dvv.ge. 10.0   .and.dvv.lt. 100.0   )then
         n = -1 
      else if(dvv.ge.  1.0   .and.dvv.lt.  10.0   )then
         n =  0 
      else if(dvv.ge.  0.10  .and.dvv.lt.   1.0   )then
         n =  1 
      else if(dvv.ge.  0.010 .and.dvv.lt.   0.10  )then
         n =  2
      else if(dvv.ge.  0.001 .and.dvv.lt.   0.01  )then
         n =  3 
      else if(dvv.ge.  0.0001.and.dvv.lt.   0.001 )then
         n =  4
      endif
      odvv = dvv * 10.0**(n)
     
      if(odvv.lt.2.0)then
            dvv = 1.0
      else if(odvv.lt.3.5)then
           dvv = 2.5
      else if(odvv.lt.7.5)then
           dvv = 5.00
      else if(odvv.lt.10.00)then
           dvv = 10.00
      endif
      dvv = dvv * 10.0**(-n)
      return
      end

      subroutine dolnxgrid(x0,y0, y1,xleng,xmax,xmin,
     1   sizex, doasis, color, solid, dominor)
      implicit none
      real x0, y0, y1, xleng, xmax, xmin, sizex
      logical doasis, dominor
      integer color, solid
c-----
c     put in tic marks along the y-axis
c
c     x0      float - position of bottom side of axis
c     y0      float - position of bottom side of axis
c     y1      float - position of top side of axis
c     xleng   float - length of Y axis
c     xmax    float - maximum value of number corresponding to
c                     top
c     xmax    float - maximum value of number corresponding to
c                     bottom
c     sizex   float - size of numbers, e.g., 10
c                   NOT USED
c     doasis  int   .true. call call plot values
c                   .false.  annotate with negative (useful for depth call call plot)
c                   NOT USED
c     color   int   - color number for the grid lines
c     solid   int   - 1 solid
c                     2 dashed
c     dominor int  .true. connect minor tics
c-----


        logical dosci
        real xl, xval, xnorm
        integer nxnorm
        real xmx, xmn, dxx, dxxx, xe
        real xx, xxv
        integer n, ifac, ilow, iup, i, j, ndxxx

        if(color .gt. 0)then
           call newpen(color)
        endif

c-----
c     put in the axis 
c-----
      xl = x0 + xleng
      call plot(x0,y0,3)
      call plot(xl,y0,2)
c-----
c     set up positioning for tics 
c     compute limits for scale 
c-----
      xmn = MIN(xmin,xmax)
      xmx = MAX(xmin,xmax)
      xxv = xmn
c-----
c     safety 
c-----
      if(xmn .eq.xmx)xmx = xmn + 1.0
c-----
c     get some idea of the size of 
c     the number to be call plotted 
c-----
        xval = MAX(ABS(xmn),ABS(xmx))
c-----
c       we set the maximum number of decimal digits as ndig, 
c       if the number is greater than 100 or less than 0.01, use
c       scientific notation
c       get the maximum number of digits for the number
c-----
        if(xval .le. 0.1 .and.xval.gt.0.0 )then
            dosci = .true.
            xnorm = alog10(1.0/xval)
            nxnorm = INT(xnorm+1)
            xnorm = 10.0** ( nxnorm)
        else if( xval .ge. 10000.0)then
            dosci = .true.
            xnorm = alog10(xval)
            nxnorm = INT(xnorm)
            xnorm = (10.0**( -nxnorm))
        else 
            dosci = .false.
            xnorm = 1.0
            nxnorm = 0.0
        endif
c-----
c       choose the maximum number of tics somehow on
c       xleng, and size of numbers and necessary decimals
c-----
        dxx = xnorm*(xmx - xmn)/5.0
        if(dxx  .eq. 0.0)dxx = 1.0
        call redodvv(dxx)
c-----
c       special case for 0 to 360
c-----
        if(XMN.eq.0.0 .and. XMX.eq.360.0)then
              dxx = xnorm*(XMX - XMN)/4.0
        endif
c-----
c       get start limits for search - we will place at most 10 tics
c       here
c-----
        n = INT(alog10(dxx))
        if(dxx.lt.1.0)n = n -1
        ifac = INT(dxx*10.0**(-n))

c-----
c       now get a new minimum value that is some multiple 
c       of ifac 10**n 
c-----
        dxx = ifac * (10.** n)
        ilow = INT(xnorm*xmn/dxx -1)
        iup  = INT(xnorm*xmx/dxx +1)
c-----
c      define a resonable range for the numbers 
c-----
        if(ifac  .eq. 1)then
            ndxxx = 5
        else if(ifac  .eq. 2)then
            ndxxx = 5
        else
            ndxxx = ifac
        endif
        dxxx = dxx/ndxxx
c-----
c       loop for labels
c       xe is the end of the previous label
c-----
        call plot(x0,y0,3)
c-----
        xe = x0
        do i=ilow,iup
            do j=0,ndxxx-1 
                xxv = i*dxx + j*dxxx 
                xxv = xxv/xnorm
                if(xxv .ge. xmn .and. xxv  .le. xmx)then
                    xx = x0 + (xxv - xmn)*xleng/(xmx - xmn)
                    call plot(xx,y0,3)
                    if(j .eq.0 )then
c-----
c       estimate the number of digits for proper labeling
c       only do this once for the call plot
c           dosci   = .false.    +0.000
c               = .true.     +0.000 10^+00
c-----
                       call plot(xx,y0,3)
                       if(solid  .eq. 2)then
                           call plotd(xx,y1,9,0.02)
                       else
                           call plot(xx,y1,2)
                       endif
                        
                    else
                        call plot(xx,y0,3)
                        if(dominor)then
                            if(solid  .eq. 2)then
                                call plotd(xx,y1,9,0.02)
                            else
                                call plot(xx,y1,2)
                            endif
                        endif
                    endif
                endif
            enddo
        enddo
        call newpen(1)
      return
      end

      subroutine dolnygrid( x0, x1,y0,yleng,ymax,ymin,
     1  sizey, doasis, color, solid, dominor)
      implicit none
      real x0, x1, y0, yleng, ymax, ymin, sizey
      logical doasis
      integer color, solid
      logical dominor
c-----
c     put in tic marks along the y-axis
c
c     x0      R     - position of bottom left side of axis
c     x1      R     - position of bottom right side of axis
c     y0      R     - position of bottom side of axis
c     yleng   R     - length of Y axis
c     ymax    R     - maximum value of number corresponding to
c                     top
c     ymax    R     - maximum value of number corresponding to
c                     bottom
c     sizey   R     - size of numbers, e.g., 10
c                   NOT USED
c     doasis  L     .true. call plot values
c                   .false. annotate with negative (useful for depth call plot)
c                   NOT USED
c     color   I     - color number for the grid lines
c     solid   I     - 1 solid
c                     2 dashed
c     dominor L     .true. connect minor tics
c-----

      real ymn, ymx, yyv, yval, ynorm
      real dyy
      real yl, yy, dyyy
      integer n,j ,i,nynorm
      integer ilow, iup, ndyyy, ifac
      logical dosci


      if(color .gt. 0)call newpen(color)
c-----
c     put in the axis 
c-----
        yl = y0 + yleng
        call plot(x0,y0,3) 
        call plot(x0,yl,2) 
c-----
c     set up positioning for tics 
c     compute limits for scale 
c-----
        ymn = MIN(ymin,ymax)
        ymx = MAX(ymin,ymax)
        yyv = ymn
c-----
c     safety */
c-----
        if(ymn.eq.ymx)ymx = ymn + 1.0
c-----
c     get some idea of the size of the number to be call plotted 
c-----
      yval = MAX(ABS(ymn),ABS(ymx))
c-----
c
c        dosci .true. or .false. to indicate if
c            scientic notation is to be used
c            and then the number is multiplied by ymorm
c            to make something between 1.00 and 9.999
c-----
        if(yval.le.0.1 .and.yval.gt.0.0 )then
            dosci = .true.
            ynorm = alog10(1.0/yval)
            nynorm = INT(ynorm+1)
            ynorm = 10.0**( nynorm)
        else if( yval.ge.10000.0)then
            dosci = .true.
            ynorm = alog10(yval)
            nynorm = INT(ynorm)
            ynorm = 10.0**(-nynorm)
        else
c-----     
c           make it look good for 
c           yval= 0.1000x to yval=9999.0
c----
            dosci = .false.
            ynorm = 1.0
            nynorm = 0
        endif

c-----
c       choose the maximum number of tics somehow on
c       yleng, and size of numbers and necessary decimals
c-----
        dyy = ynorm*(ymx - ymn)/5.0
        if(dyy .eq. 0.0)dyy = 1.0

        call redodvv(dyy)
c-----
c       get start limits for search - we will place at most 10 tics
c       here
c-----
        n = INT(alog10(dyy))
        if(dyy.lt.1.0)n = n -1
        ifac = INT(dyy*10.0**(-n))
c-----
c       now get a new minimum value that is some multiple 
c       of ifac 10**n 
c-----
        dyy = ifac * 10.0**(n)
c-----
c       set major tic increment
c-----
        ilow = INT(ynorm*ymn/dyy -1)
        iup  = INT(ynorm*ymx/dyy +1)
        if(ifac .eq. 1)then
            ndyyy = 5
        else if(ifac .eq. 2)then
            ndyyy = 5
        else
            ndyyy = ifac
        endif
        dyyy = dyy/ndyyy
c-----
c       loop for ticmarks
c-----
        call plot(x0,y0,3)
c-----
c       here
c       yy is the position of the call plot corresponding to a 
c       value of yyv
c-----
        do i=ilow,iup
            do j=0,ndyyy-1
            
c-----
c          yyv is the value to be plotted
c-----
              yyv = i*dyy + j*dyyy
c-----
c          now scale this to relate to original number
c----
            yyv = yyv/ynorm
            if(yyv .ge. ymn .and. yyv  .le. ymx)then
            yy = y0 + (yyv - ymn)*yleng/(ymx - ymn)
                call plot(x0,yy,3)
                if(j.eq.0)then
                    if(solid .eq. 2)then
                        call plotd(x1,yy,9,0.02)
                    else
                        call plot(x1,yy,2)
                    endif
                else 
                    if(dominor)then
                        if(solid .eq. 2) then
                            call plotd(x1,yy,9,0.02)
                        else
                            call plot(x1,yy,2)
                        endif
                    endif
                endif
            endif
            enddo
       enddo
       call newpen(1)
      return
      end
      subroutine dologygrid(x0,x1,y0,yleng,ymax,ymin,
     1  sizey,color, solid, dominor)

      implicit none
      real x0, x1, y0, yleng, ymax, ymin, sizey
      integer color, solid
      logical dominor
c-----
c      put in tic marks along the y-axis 
c
c       x0      R    - position of bottom side of axis
c       y0      R    - position of bottom side of axis
c       yleng   R    - length of Y axis
c       ymax    R    - maximum value of number corresponding to
c                      top
c       ymin    R    - minimum value of number corresponding to
c                      bottom
c       sizey   R    - size of numbers, e.g., 10
c       color   I    - color of grid
c       solid   I    - 1 solid line
c                      2 dashed line
c       dominor L     - .true.  do minor tics
c-----
        integer nocy 
        integer iy,iiy
        real ymxlog, ymmin,ymlog
        integer ja, jpow,ii,jj
        real yscal
        real yy
        real tenpow, yval, yv

        real num10(9)
        data num10/1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0/


        if(color .gt. 0)call newpen(color)
c-----
c       put in the axis 
c-----
        call plot(x0,y0,3)
        call plot(x0,y0+yleng,2)
c-----
c       set up positioning for 10**N labels 
c       compute limits for scale 
c-----
        ymxlog = alog10(ymax)
        ymmin = alog10(ymin)
        nocy = INT(ymxlog - ymmin + 1)
        iy = INT(ymmin)
        if(ymmin .lt. 0)iy = iy - 1

        ymlog = alog10(ymax)
        iiy = INT(ymmin)
        if(ymlog .lt. 0)iiy = iiy - 1

        ja = MAX(ABS(iy),ABS(iiy))
        if(ja.lt.10)then
            jpow = 1
        else if(ja .ge. 10 .and. ja.lt.100)then
            jpow = 2
        else if(ja .ge. 100 .and. ja.lt.1000)then
            jpow = 3
        endif
        ii = MIN(iy,iiy)
        if(ii.lt.0)jpow = jpow + 1
        jpow = jpow + 2


        yscal = yleng/alog10(ymax/ymin)
        do ii=iy,iy+nocy+2
            tenpow = 10.0**ii
            do jj=1,9
                yval = num10(jj)*tenpow
c-----
c       put in tics where appropriate 
c-----
                if(yval  .ge.  ymin .and. yval .le. ymax)then
                    yv = alog10(yval/ymin)
                    yy = y0 + yv*yscal
                    if(jj.eq.1)then
                        if(solid .eq. 1)then
                            call plot(x0,yy,3)
                            call plot(x1,yy,2)
                        else if(solid .eq. 2)then
                            call plot(x0,yy,3)
                            call plotd(x1,yy,9,0.02)
                        endif
                    else 
                        if(dominor )then
                            if(solid .eq. 1)then
                                call plot(x0,yy,3)
                                call plot(x1,yy,2)
                            else if(solid .eq. 2)then
                                call plot(x0,yy,3)
                                call plotd(x1,yy,9,0.02)
                            endif
                        endif
                    endif
                    call plot(x0,yy,3)
                endif
            enddo
        enddo
        call plot(x0,y0,3)
        call newpen(1)
      return
      end

      subroutine  dologxgrid(x0,y0,y1,xleng,sxmax,sxmin,
     1  sizex,color, solid, dominor)
      implicit none
      real x0, y0, y1,xleng,sxmax,sxmin,sizex
      integer color, solid
      logical dominor
c-----
c       put in tic marks along the x-axis 
c
c       x0  R*4 - position of left side of axis
c       y0  R*4 - position of left side of axis
c       xleng   R*4 - length of X axis
c       xmax    R*4 - maximum value of number corresponding to
c                   far right
c       xmin    R*4 - maximum value of number corresponding to
c                   far left
c       sizex   R*4 - size of numbers, e.g., 10
c       color   I       - color of grid
c       solid   I       - 1 solid line
c                         2 dashed line
c       dominor L       - .gt. 0 do minor tics
c-----
        integer nocx 
        integer ix
        real xmxlog, xmmin
        integer ja, jpow,ii,jj
        real xscal
        real xx
        real tenpow, xval, xv
        real xmin, xmax

        real num10(9)
        data num10/1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0/

        if(color .gt. 0)call newpen(color)

        xmin = sxmin
        xmax = sxmax
c-----
c       put in the axis 
c-----
        call plot(x0,y0,3)
        call plot(x0+xleng,y0,2)
c-----
c       set up positioning for 10**N labels 
c       compute limits for scale 
c-----
        xmxlog = alog10(xmax)    
        xmmin = alog10(xmin) 
        nocx = INT(xmxlog - xmmin + 1)
        ix = INT(xmmin)
        if(xmmin .lt. 0)ix = ix - 1
        do ii=ix,ix+nocx+2
            tenpow = 10.0**ii
            ja = ABS(ii)
            if(ja.lt.10)then
                jpow = 1
            else if(ja .ge. 10 .and. ja.lt.100)then
                jpow = 2
            else if(ja .ge. 100 .and. ja.lt.1000)then
                jpow = 3
            endif
            if(ii.lt.0)jpow = jpow + 1
            xscal = xleng/alog10(xmax/xmin)
            do jj=1,9
                xval = num10(jj)*tenpow
c-----
c               put in tics where appropriate 
c-----
                if(xval  .ge.  xmin .and. xval .le. xmax)then
                    xv = alog10(xval/xmin)
                    xx = x0 + xv*xscal
                    if(jj.eq.1)then
                        if(solid .eq. 1)then
                              call plot(xx,y0,3) 
                              call plot(xx,y1,2)
                        else if(solid .eq. 2)then
                              call plot(xx,y0,3)
                              call plotd(xx,y1,9,0.02)
                        endif
                    else 
                          if(dominor )then
                              if(solid .eq. 1)then
                                 call plot(xx,y0,3)
                                 call plot(xx,y1,2)
                              else if(solid .eq. 2)then
                                 call plot(xx,y0,3)
                                 call plotd(xx,y1,9,0.02)
                              endif
                          endif
                    endif
                    call plot(xx,y0,3)
                endif
            enddo
        enddo
        call plot(xx,y0,3)
        call newpen(1)
      return
      end
