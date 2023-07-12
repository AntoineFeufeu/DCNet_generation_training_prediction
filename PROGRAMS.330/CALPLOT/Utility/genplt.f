      program genplt
c-----
c     CHANGES
c     13 OCT 2010 created
c     08 FEB 2012 merged tgenplt with genplt
c        the extra control of tgenplt is now given by the
c        -A acmdfil   while the older control is given by 
c        -C cmdfil
c        the difference is that -A acmdfil permits setting the 
c        symbol size while the -C cmdfil uses a default symbol size
c     08 DEC 2012 added option -L cmdfil which is an extension to
c        -A cmdfil to permot annotation
c     07 JUN 2013 actually implemented -NOBOX
c     14 OCT 2015 fixed error in doliny which has tilft instead of ticlft
c     17 NOV 2016 make the lines slightly longer for the NO and DA legends
c     01 OCT 2018 add DO for dotted
c     13 JAN 2020 added polygon
c     01 OCT 2020 added error bar
c-----

      integer LIN, LOT
      parameter (LIN=5,LOT=6)

      real x0, y0, xlen, ylen
      real xmin, xmax, ymin, ymax

      logical dobox, doxdown
      logical doxlin, doylin
      logical doticx, doticy
      character titlex*80
      character titley*80


      integer MXFILE
      parameter(MXFILE=500)

      integer nfile
      character*100 fname(MXFILE)
      real width(MXFILE)
      integer kolor(MXFILE)
      character*2 psymb(MXFILE)
      character*1 ftype(MXFILE)
      real symsiz(MXFILE)
      character*100 symleg(MXFILE)
      character lpos*4
      integer nleg, ncharleg
      
      real sizex, sizey
      logical ticup, labtop, dopow
      logical ticlft, lablft

      integer llx, lly
      integer lgstr


      real lx0, ly0, lht, lxx, lyy, ldy

      call gcmdln(x0,y0,xlen,ylen,titlex,titley,doxlin,doylin,
     1  xmin,xmax,ymin,ymax,nfile,fname,ftype,width,kolor,psymb,symsiz,
     1  lpos,symleg,nleg,ncharleg,dobox,doxdown,doticx,doticy)
c-----
c     define legend position according to lpos string
c     adjust sixe to fit
c-----
      if(doxdown)then
           xaxlen = ylen
           yaxlen = xlen
      else
           xaxlen = xlen
           yaxlen = ylen
      endif
      write(0,*)nleg
      if(nleg.gt.0)then
          ldy = 0.5*ylen/nleg
      else
          ldy = 0.14
      endif
      if(ldy.gt.0.14)ldy = 0.14
      if(lpos(1:1).eq.'T')then
            ly0 = y0 + yaxlen - 0.35 
      else if(lpos(1:1).eq.'M')then
            ly0 = y0 + 0.5*yaxlen +0.5*nleg*ldy
      else if(lpos(1:1).eq.'B')then
            ly0 = y0 +   0.25 + nleg*ldy
      endif
      if(ldy .gt. 0.10)then
            lht = 0.10
      else
            lht = 0.7 * ldy
      endif
      if(lpos(2:2).eq.'L')then
            lx0 = x0 + 0.40
      else if(lpos(2:2).eq.'C')then
            lx0 = x0 + 0.50*xaxlen
      else if(lpos(2:2).eq.'R')then
            lx0 = x0 + xaxlen - 0.25 -(ncharleg+3)*lht
      endif
C     WRITE(0,*)'lx0,ly0,lht,ldy:',lx0,ly0,lht,ldy

      call pinitf('GENPLT.PLT')

      llx = lgstr(titlex)
      lly = lgstr(titley)
C     WRITE(6,*)'llx=',llx
      sizex = 0.12
      sizey = 0.12
      call doframe(x0,y0,xlen,xmax,xmin,
     1     ylen,ymax,ymin,
     1     sizex, sizey,
     1     llx,titlex, lly,titley,
     1     doxlin,doylin,dobox,doxdown,doticx,doticy)

      call gclip('ON',x0,y0,x0+xaxlen,y0+yaxlen)
      do nf = 1, nfile
      WRITE(6,*)nf,nfile,fname(nf),kolor(nf),width(nf),psymb(nf),
     1     symsiz(nf),ftype(nf)
        open(1,file=fname(nf),access='sequential',
     1   form = 'formatted',status='unknown')
        rewind 1
C        if(psymb(nf).eq.'NO')then
           call gwidth(width(nf))
C        else 
C           call gwidth(0.10*symsiz(nf))
C        endif

        if(ftype(nf).eq.'P')then
          call dopoly(nf,kolor(nf),x0,xlen,xmin,xmax,
     1      y0,ylen,ymin,ymax,doxlin,doylin,doxdown,
     2      psymb(nf),symsiz(nf),width(nf),ftype(nf),
     3      symleg(nf) ,lx0,ly0,lht,ldy)
        else
          call dosym(nf,kolor(nf),x0,xlen,xmin,xmax,
     1      y0,ylen,ymin,ymax,doxlin,doylin,doxdown,
     2      psymb(nf),symsiz(nf),width(nf),ftype(nf),
     3      symleg(nf) ,lx0,ly0,lht,ldy)
        endif
      enddo
      call gclip('OFF',x0,y0,x0+xaxlen,y0+yaxlen)
c-----
c     reset the state
c-----
      call gwidth(0.0)
      call newpen(1)
      call pend()
      end

        subroutine dopoly(nf,kolor,x0,xlen,xmin,xmax,
     1      y0,ylen,ymin,ymax,doxlin,doylin,doxdown,
     2      psymb,symsiz,width,ftype,symleg,lx0,ly0,
     3      lht,ldy)
        implicit none
        integer MXPTS
        parameter (MXPTS=200)
        real xc(MXPTS), yc(MXPTS)

        integer nf, kolor
        real x0,xlen, xmin, xmax, y0,ylen, ymin,ymax
        logical doxlin, doylin, doxdown
        character psymb*(*), ftype*(*), symleg*(*)
        real symsiz,width,lx0,ly0
        real lht , ldy

        integer jpen, ipen
        integer lgstr, i
        real x, y, xx, yy
        integer ls
        real lxx, lyy

        integer npts

        ipen = 3
c-----
c     first determine the number of points
c     abort if exceed MXPTS
c-----
      rewind 1
      npts = 0
 1000 continue
      read(1,*,end=2000)x,y
      npts = npts +1
      if(npts.gt.MXPTS)then
          close(1)
          return
      endif
        if(doxdown)then
                yy = y0+xlen - xlen*(x-xmin)/(xmax-xmin)
            if(doylin)then
                xx = x0 + ylen*(y-ymin)/(ymax-ymin)
            else
                xx = x0 + ylen*alog10(y/ymin)/alog10(ymax/ymin)
            endif
        else
            if(doxlin)then
                xx = x0 + xlen*(x-xmin)/(xmax-xmin)
            else
                xx = x0 + xlen*alog10(x/xmin)/alog10(xmax/xmin)
           endif
            if(doylin)then
                yy = y0 + ylen*(y-ymin)/(ymax-ymin)
            else
                yy = y0 + ylen*alog10(y/ymin)/alog10(ymax/ymin)
            endif
        endif
      xc(npts) = xx
      yc(npts) = yy
      go to 1000
 2000 continue
      close(1)
      call newpen(kolor)
      if(psymb .eq. 'Y')then
           call shadep(npts,xc,yc)
      else
c-----
c          draw the polygon
c----
           call plot(xc(1),yc(2),3)
           do i=1,npts
               call plot(xc(i),yc(i),2)
           enddo
           call plot(xc(1),yc(2),2)
           call plot(xc(1),yc(2),3)
      endif
      call newpen(1)
      return
      end

        subroutine dosym(nf,kolor,x0,xlen,xmin,xmax,
     1      y0,ylen,ymin,ymax,doxlin,doylin,doxdown,
     2      psymb,symsiz,width,ftype,symleg,lx0,ly0,
     3      lht,ldy)
        implicit none
        integer nf, kolor
        real x0,xlen, xmin, xmax, y0,ylen, ymin,ymax
        logical doxlin, doylin, doxdown
        character psymb*(*), ftype*(*), symleg*(*)
        real symsiz,width,lx0,ly0
        real lht , ldy

        integer jpen, ipen, ipenold
        integer lgstr
        real x, y, xx, yy, dx, dy, xxp,xxm,yyp,yym
        integer ls
        real lxx, lyy

        ipen = 3
 1000   continue
        call newpen(kolor)
        if(ftype.eq.'3')then
             read(1,*,end=9999)x,y,jpen
             dx = 0
             dy = 0
        else if(ftype.eq.'E')then
             read(1,*,end=9999)x,y,dx,dy
             jpen = -1
        else
             read(1,*,end=9999)x,y
             dx = 0
             dy = 0
             jpen = -1
        endif
C        WRITE(6,*)x,y
        if(doxdown)then
C           if(doxlin)then
                yy = y0+xlen - xlen*(x-xmin)/(xmax-xmin)
C           else
C               yy = y0 + xlen*alog10(x/xmin)/alog10(xmax/xmin)
C           endif
            if(doylin)then
                xx = x0 + ylen*(y-ymin)/(ymax-ymin)
            else
                xx = x0 + ylen*alog10(y/ymin)/alog10(ymax/ymin)
            endif
C     WRITE(6,*)x,xmin,xmax,yy,y,ymin,ymax,xx,xlen,ylen,x0,y0
        else
            if(doxlin)then
                xx = x0 + xlen*(x-xmin)/(xmax-xmin)
                xxp= x0 + xlen*(x+dx-xmin)/(xmax-xmin)
                xxm= x0 + xlen*(x-dx-xmin)/(xmax-xmin)
            else
                xx = x0 + xlen*alog10(x/xmin)/alog10(xmax/xmin)
                xxp= x0 + xlen*alog10((x+dx)/xmin)/alog10(xmax/xmin)
                xxm= x0 + xlen*alog10((x-dx)/xmin)/alog10(xmax/xmin)
            endif
            if(doylin)then
                yy = y0 + ylen*(y-ymin)/(ymax-ymin)
                yyp= y0 + ylen*(y+dy-ymin)/(ymax-ymin)
                yym= y0 + ylen*(y-dy-ymin)/(ymax-ymin)
            else
                yy = y0 + ylen*alog10(y/ymin)/alog10(ymax/ymin)
                yyp= y0 + ylen*alog10((y+dy)/ymin)/alog10(ymax/ymin)
                yym= y0 + ylen*alog10((y-dy)/ymin)/alog10(ymax/ymin)
            endif
        endif
                if(psymb.eq.'NO')then
                    call plot(xx,yy,ipen)
                    ipen = 2
                else if(psymb.eq.'DA')then
                          if(ipen.eq.3)then
                             call plot(xx,yy,ipen)
                          else
                             call plotd(xx,yy,15,0.05)
                          endif
                          ipen = 2
                else if(psymb.eq.'DO')then
                          if(ipen.eq.3)then
                             call plot(xx,yy,ipen)
                          else
                             call plotd(xx,yy,18,0.05)
                          endif
                          ipen = 2
                else
                    if(ftype.eq.'E')then
                          call newpen(1)
                          call errbar(xx,yy,xxm,xxp,yym,yyp,3)
                    endif
                    if(jpen.ge.0 )then
                         call newpen(jpen)
                    else
                         call newpen(kolor)
                    endif
                    call fillit(psymb,symsiz,xx,yy)
                    if(width.gt.0.0)then
                    call newpen(1)
                    call curvit(psymb,symsiz,xx,yy)
                    endif
                endif
        go to 1000
 9999   continue
c-----
c       clean up by lifting the pen
c-----
        ipen = 3
               if(psymb.eq.'NO')then
                    call plot(xx,yy,ipen)
               else if(psymb.eq.'DA')then
                    call plot(xx,yy,ipen)
               else if(psymb.eq.'DO')then
                    call plot(xx,yy,ipen)
               endif
        close (1)
c-----
c       plot legend
c-----
        ls = lgstr(symleg)
        WRITE(0,*)symleg(1:ls),lx0,ly0,lht,ldy
        if(symleg .ne. ' ' )then
            lxx = lx0
            lyy = ly0 
                if(psymb.eq.'NO')then
                    call plot(lxx-lht,lyy+0.5*lht,3)
                    call plot(lxx+2.*lht,lyy+0.5*lht,2)
                else if(psymb.eq.'D1')then
                    call plot(lxx-lht,lyy+0.5*lht,3)
                    call plotd(lxx+2.*lht,lyy+0.5*lht,15,0.05)
                    call plot(lxx-lht,lyy+0.5*lht,3)
                else if(psymb.eq.'DO')then
                    call plot(lxx-lht,lyy+0.5*lht,3)
                    call plotd(lxx+2.*lht,lyy+0.5*lht,18,0.05)
                    call plot(lxx-lht,lyy+0.5*lht,3)
                else
                    call newpen(kolor)
                    call fillit(psymb,symsiz,lxx,lyy+0.5*lht)
                    if(width.gt.0.0)then
                    call newpen(1)
                    call curvit(psymb,symsiz,lxx,lyy+0.5*lht)
                    endif
                endif
            call newpen(1)
            call gwidth(0.0)
            call gleft(lxx+3.*lht,lyy,lht,symleg(1:ls),0.0)
            ly0 = ly0 - ldy
     
         endif
      return
      end

      subroutine errbar(xx,yy,xxm,xxp,yym,yyp,ipen)
                   if(dx.gt.0.0)then
                      call plot(xxm,yy,3)
                      call plot(xxp,yy,2)
                      call plot(xx,yy,ipen)
                   endif
                   if(dy.gt.0.0)then
                      call plot(xx,yyp,3)
                      call plot(xx,yym,2)
                      call plot(xx,yy,ipen)
                   endif
      return
      end

      subroutine doframe(x0,y0,xlen,xmax,xmin,
     1     ylen,ymax,ymin,
     1     sizex, sizey,
     1     llx,titlex, lly,titley,
     1     doxlin,doylin,dobox,doxdown,doticx,doticy)

      implicit none
      real x0,y0
      real xlen,xmax,xmin,sizex
      real ylen,ymax,ymin,sizey
      logical doxlin,doylin,dobox, doxdown,doticx,doticy
      integer llx, lly
      character titlex*(*), titley*(*)

c------
c     test if the axes are to be plotted
c-----
      if( .not. dobox)then
           return
      else
           if(doxdown)then
           call gbox(x0,y0,x0+ylen,y0+xlen)
           else
           call gbox(x0,y0,x0+xlen,y0+ylen)
           endif
      endif
c-----
c     select the method of plot, e.g.,
c     x - increasing to the right, y increasing upward, or
c     x - increasing downward, y increasing to the right
c-----


      if(doxdown)then
         call dofrmxdown(x0,y0,xlen,xmax,xmin,
     1     ylen,ymax,ymin,
     1     sizex, sizey,
     1     llx,titlex, lly,titley,
     1     doxlin,doylin,doticx,doticy)
      else
         call dofrmxrght(x0,y0,xlen,xmax,xmin,
     1     ylen,ymax,ymin,
     1     sizex, sizey,
     1     llx,titlex, lly,titley,
     1     doxlin,doylin,doticx,doticy)
      endif
      return
      end

      subroutine dofrmxrght(x0,y0,xlen,xmax,xmin,
     1     ylen,ymax,ymin,
     1     sizex, sizey,
     1     llx,titlex, lly,titley,
     1     doxlin,doylin,doticx,doticy)

      implicit none
      real x0,y0
      real xlen,xmax,xmin,sizex
      real ylen,ymax,ymin,sizey
      logical doxlin,doylin,doticx,doticy
      integer llx, lly
      character titlex*(*), titley*(*)


      integer lgstr

      logical dopow
      logical ticup, ticlft
      logical labtop, lablft

c-----
c          do x-axis -- horizontal
c-----
      if(doticx)then
           if(doxlin)then
             ticup = .true.
             labtop = .false.
             dopow = .true.
             call dolinx(x0,y0,xlen,xmax,xmin,
     1           sizex,ticup,labtop,dopow,llx,titlex)
             ticup = .false.
             labtop = .false.
             dopow = .false.
             call dolinx(x0,y0+ylen,xlen,xmax,xmin,
     1           sizex,ticup,labtop,dopow,llx,titlex)
           else
             ticup = .true.
             labtop = .false.
             dopow = .true.
             call dologx(x0,y0,xlen,xmax,xmin,
     1           sizex,ticup,labtop,dopow,llx,titlex)
             ticup = .false.
             labtop = .false.
             dopow = .false.
             call dologx(x0,y0+ylen,xlen,xmax,xmin,
     1           sizex,ticup,labtop,dopow,llx,titlex)
           endif
       endif
           sizey = 0.15
c-----
c          do y-axis -- vertical
c-----
       if(doticy)then
           if(doylin)then
             ticlft = .false.
             lablft = .true.
             dopow = .true.
             call doliny(x0,y0,ylen,ymax,ymin,
     1           sizey,ticlft,lablft,dopow,lly,titley)
             ticlft = .true.
             lablft = .false.
             dopow = .false.
             call doliny(x0+xlen,y0,ylen,ymax,ymin,
     1           sizey,ticlft,lablft,dopow,lly,titley)
           else
             ticlft = .false.
             lablft = .true.
             dopow = .true.
             call dology(x0,y0,ylen,ymax,ymin,
     1           sizey,ticlft,lablft,dopow,lly,titley)
             ticlft = .true.
             lablft = .false.
             dopow = .false.
             call dology(x0+xlen,y0,ylen,ymax,ymin,
     1           sizey,ticlft,lablft,dopow,lly,titley)
           endif
       endif
      return
      end

      subroutine dofrmxdown(x0,y0,xlen,xmax,xmin,
     1     ylen,ymax,ymin,
     1     sizex, sizey,
     1     llx,titlex, lly,titley,
     1     doxlin,doylin,doticx,doticy)
c-----
c     NOTE this will ALWAYS use linear scale for x
c-----

      implicit none
      real x0,y0
      real xlen,xmax,xmin,sizex
      real ylen,ymax,ymin,sizey
      logical doxlin,doylin,doticx,doticy
      integer llx, lly
      character titlex*(*), titley*(*)


      integer lgstr

      logical dopow
      logical ticup, ticlft
      logical labtop, lablft

c-----
c          do x-axis -- horizontal
c-----
      if(doticy)then
           if(doylin)then
             ticup = .true.
             labtop = .false.
             dopow = .false.
             call dolinx(x0,y0,ylen,ymax,ymin,
     1           sizey,ticup,labtop,dopow,lly,titley)
             ticup = .false.
             labtop = .true.
             dopow = .true.
             call dolinx(x0,y0+xlen,ylen,ymax,ymin,
     1           sizey,ticup,labtop,dopow,lly,titley)
           else
             ticup = .true.
             labtop = .false.
             dopow = .false.
             call dologx(x0,y0,ylen,ymax,ymin,
     1           sizey,ticup,labtop,dopow,lly,titley)
             ticup = .false.
             labtop = .true.
             dopow = .true.
             call dologx(x0,y0+xlen,ylen,ymax,ymin,
     1           sizey,ticup,labtop,dopow,lly,titley)
           endif
       endif
       if(doticx)then
c-----
c          do y-axis -- vertical - always linear
c-----
           if(llx.gt.0)then
               call dnliny(x0     ,y0,xlen,-xmax,-xmin,sizey,
     1             .false.,.true. ,.true. ,llx,titlex)
               call dnliny(x0+ylen,y0,xlen,-xmax,-xmin,sizey,
     1             .true. ,.true. ,.false.,1,' ')
           else
               call dnliny(x0     ,y0,xlen,-xmax,-xmin,sizey,
     1             .false.,.true. ,.false.,1 ,' ')
               call dnliny(x0+ylen,y0,xlen,-xmax,-xmin,sizey,
     1             .true. ,.true. ,.false.,1,' ')
           endif
       endif
      return
      end

      subroutine gcmdln(x0,y0,xlen,ylen,titlex,titley,doxlin,doylin,
     1  xmin,xmax,ymin,ymax,nfile,fname,ftype,width,kolor,psymb,symsiz,
     1  lpos,symleg,nleg,ncharleg,dobox,doxdown,doticx,doticy)
      real x0,y0,xlen,ylen,xmin,xmax,ymin,ymax
      logical doxlin, doylin,doticx,doticy
      character titlex*80
      character titley*80

      integer MXFILE
      parameter(MXFILE=500)

      integer nfile
      character*100 fname(MXFILE)
      real width(MXFILE)
      integer kolor(MXFILE)
      character*2 psymb(MXFILE)
      character*1 ftype(MXFILE)
      real symsiz(MXFILE)
      character lpos*4
      character*100 symleg(MXFILE)
      integer nleg,ncharleg
      logical dobox, doxdown
      

      
 
      character*80 names
      integer mnmarg
      integer i
      integer ls, ln
      character cmdfil*80
      character acmdfil*80
      character ecmdfil*80
      character lcmdfil*80
      character l3cmdfil*80
      character pcmdfil*80
      character tfile*100, tsymb*2, tleg*100
      real ifill
      logical ext

c-----
c     defaults
c-----

      x0 = 2
      y0 = 1
      xlen = 5
      ylen = 5
      doxlin = .true.
      doylin = .true.
      xmin = 1.0
      xmax = 10.0
      titlex = 'X-axis'
      titley = 'Y-axis'
      nfile = 0
      cmdfil = ' '
      acmdfil = ' '
      ecmdfil = ' '
      lcmdfil = ' '
      l3cmdfil = ' '
      pcmdfil = ' '
      nleg = 0
      lpos = ''
      ncharleg = 0
      dobox = .true.
      doxdown = .false.
      doticx = .true.
      doticy = .true.

        nmarg = mnmarg()
        i = 0
 1000   continue
            i = i + 1
            if(i.gt.nmarg)go to 2000
         names=''
            call mgtarg(i,names)
         ln = lgstr(names)
            if(names(1:5).eq.'-XMIN')then
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
                read(names,'(bn,f20.0)')ymin
            else if(names(1:5).eq.'-YMAX')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f20.0)')ymax
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
            else if(names(1:5).eq.'-XLOG')then
                doxlin = .false.
            else if(names(1:5).eq.'-XLIN')then
                doxlin = .true.
            else if(names(1:5).eq.'-YLOG')then
                doylin = .false.
            else if(names(1:5).eq.'-YLIN')then
                doylin = .true.
            else if(names(1:4).eq.'-NOB')then
                dobox = .false.
            else if(names(1:7).eq.'-NOTICX')then
                doticx = .false.
            else if(names(1:7).eq.'-NOTICY')then
                doticy = .false.
            else if(names(1:6).eq.'-XDOWN')then
                doxdown = .true.
            else if(names(1:2).eq.'-C')then
                i = i + 1
                names = ' '
                call mgtarg(i,names)
                ls = lgstr(names)
                cmdfil = names(1:ls)
            else if(names(1:2).eq.'-A')then
                i = i + 1
                names = ' '
                call mgtarg(i,names)
                ls = lgstr(names)
                acmdfil = names(1:ls)
            else if(names(1:2).eq.'-E')then
                i = i + 1
                names = ' '
                call mgtarg(i,names)
                ls = lgstr(names)
                ecmdfil = names(1:ls)
            else if(names(1:2).eq.'-P')then
                i = i + 1
                names = ' '
                call mgtarg(i,names)
                ls = lgstr(names)
                pcmdfil = names(1:ls)
            else if(names(1:2).eq.'-3')then
                      i = i + 1
                      names = ' '
                      call mgtarg(i,names)
                      ls = lgstr(names)
                      l3cmdfil = names(1:ls)
            else if(names(1:2).eq.'-L')then
                if(ln.eq.2)then
                      i = i + 1
                      names = ' '
                      call mgtarg(i,names)
                      ls = lgstr(names)
                      lcmdfil = names(1:ls)
                 else if(names(1:5).eq.'-LPOS')then
                      i = i + 1
                      call mgtarg(i,names)
                      ln = lgstr(names)
                      if(ln.eq.2)then
                         lpos = names(1:2)
                      endif
                 endif
            else if(names(1:3).eq.'-TX')then
                i=i+1
                call mgtarg(i,titlex)
            else if(names(1:3).eq.'-TY')then
                i=i+1
                call mgtarg(i,titley)
            else if(names(1:2).eq.'-?')then
                call usage(' ')
            else if(names(1:2).eq.'-h')then
                call usage(' ')
            endif
        go to 1000
 2000   continue
c-----
c       this is ugly since each type fo input must be read in
c-----
c-----
c       now open the cmdfil and read the entries
c-----
        if(cmdfil .eq. ' ' .and. acmdfil .eq. ' '
     1     .and. lcmdfil .eq. ' ' .and. l3cmdfil .eq. ' ' 
     2     .and. pcmdfil .eq. ' ' .and. ecmdfil .eq. ' ')then
            call usage('No cmdfil or acmdfil')
        endif
        nfile = 0
c-----
c       now open the lcmdfil and read the entries
c-----
        if(pcmdfil .ne. ' ' )then
            inquire(file=pcmdfil,exist=ext)
            if(ext)then
         
            
                open(1,file=pcmdfil,status='unknown',form='formatted',
     1            access='sequential')
                rewind 1
c-----
c               this is now followed by the following items
c               FILE KOLOR WIDTH PSYMB
c               FILE file name of x-y pairs to be plotted
c               KOLOR 1 - BLACK 1000=red 1050=green 1100=blue 0 = white
c-----
 3004          continue
                    read(1,*,end=4004)tfile,kkol,fill
                    nfile = nfile + 1
                    fname(nfile) = tfile
                    ftype(nfile) = 'P'
                    kolor(nfile) = kkol
                    width(nfile) = -1
                    if(fill .ne.0.0)then
                        psymb(nfile)='Y'
                    else
                        psymb(nfile)='N'
                    endif
                    symsiz(nfile) = -1
                    symleg(nfile) = ' '
                go to 3004
 4004           continue
                close(1)
            endif
        endif
        if(cmdfil .ne. ' ' )then
            inquire(file=cmdfil,exist=ext)
            if(ext)then
         
            
                open(1,file=cmdfil,status='unknown',form='formatted',
     1            access='sequential')
                rewind 1
c-----
c               this is now followed by the following items
c               FILE KOLOR WIDTH PSYMB
c               FILE file name of x-y pairs to be plotted
c               KOLOR 1 - BLACK 1000=red 1050=green 1100=blue 0 = white
c               WIDTH width of line in inches
c               PSYMB - a 2 character entry with the following meaning
c                   SQ - square
c                   TR - triangle
c                   HX - heaxgon
c                   DI - diamond
c                   CI - circle
c                   DA - dashed
c                   D0 - dotted
c                   NO - no symbol - plot a line
c-----
 3000          continue
                    read(1,*,end=4000)tfile,kkol,twid,tsymb
                    nfile = nfile + 1
                    fname(nfile) = tfile
                    ftype(nfile) = 'C'
                    kolor(nfile) = kkol
                    width(nfile) = twid
                    psymb(nfile)=tsymb
                        symsiz(nfile) = 0.1
                    symleg(nfile)=' '
                go to 3000
 4000           continue
                close(1)
            endif
        endif
c-----
c       now open the acmdfil and read the entries
c-----
        if(acmdfil .ne. ' ' )then
            inquire(file=acmdfil,exist=ext)
            if(ext)then
         
            
                open(1,file=acmdfil,status='unknown',form='formatted',
     1            access='sequential')
                rewind 1
c-----
c               this is now followed by the following items
c               FILE KOLOR WIDTH PSYMB
c               FILE file name of x-y pairs to be plotted
c               KOLOR 1 - BLACK 1000=red 1050=green 1100=blue 0 = white
c               WIDTH width of line in inches
c               PSYMB - a 2 character entry with the following meaning
c                   SQ - square
c                   TR - triangle
c                   HX - heaxgon
c                   DI - diamond
c                   CI - circle
c                   NO - no symbol - plot a line
c-----
 3001          continue
                    read(1,*,end=4001)tfile,kkol,twid,tsymb,tsiz
                    nfile = nfile + 1
                    fname(nfile) = tfile
                    ftype(nfile) = 'A'
                    kolor(nfile) = kkol
                    width(nfile) = twid
                    psymb(nfile)=tsymb
                    if(tsiz.le.0.0)tsiz = 0.05
                    symsiz(nfile) = tsiz
                    symleg(nfile)=' '
                go to 3001
 4001           continue
                close(1)
            endif
        endif
c-----
c       now open the ecmdfil and read the entries
c-----
        if(ecmdfil .ne. ' ' )then
            inquire(file=ecmdfil,exist=ext)
            if(ext)then
         
            
                open(1,file=ecmdfil,status='unknown',form='formatted',
     1            access='sequential')
                rewind 1
c-----
c               this is now followed by the following items
c               FILE KOLOR WIDTH PSYMB
c               FILE file name of x-y pairs to be plotted
c               KOLOR 1 - BLACK 1000=red 1050=green 1100=blue 0 = white
c               WIDTH width of line in inches
c               PSYMB - a 2 character entry with the following meaning
c                   SQ - square
c                   TR - triangle
c                   HX - heaxgon
c                   DI - diamond
c                   CI - circle
c                   NO - no symbol - plot a line
c-----
 3008          continue
                    read(1,*,end=4008)tfile,kkol,twid,tsymb,tsiz,tleg
                    nfile = nfile + 1
                    fname(nfile) = tfile
                    ftype(nfile) = 'E'
                    kolor(nfile) = kkol
                    width(nfile) = twid
                    psymb(nfile)=tsymb
                    if(tsiz.le.0.0)tsiz = 0.05
                    symsiz(nfile) = tsiz
                    ls = lgstr(tleg)
                    symleg(nfile) = tleg(1:ls)
                    nleg = nleg + 1
                    if(ls.gt.ncharleg)ncharleg=ls
                go to 3008
 4008           continue
                close(1)
            endif
        endif
c-----
c       now open the lcmdfil and read the entries
c-----
        if(lcmdfil .ne. ' ' )then
            inquire(file=lcmdfil,exist=ext)
            if(ext)then
         
            
                open(1,file=lcmdfil,status='unknown',form='formatted',
     1            access='sequential')
                rewind 1
c-----
c               this is now followed by the following items
c               FILE KOLOR WIDTH PSYMB
c               FILE file name of x-y pairs to be plotted
c               KOLOR 1 - BLACK 1000=red 1050=green 1100=blue 0 = white
c               WIDTH width of line in inches
c               PSYMB - a 2 character entry with the following meaning
c                   SQ - square
c                   TR - triangle
c                   HX - heaxgon
c                   DI - diamond
c                   CI - circle
c                   NO - no symbol - plot a line
c-----
 3002          continue
                    read(1,*,end=4002)tfile,kkol,twid,tsymb,tsiz,tleg
                    nfile = nfile + 1
                    fname(nfile) = tfile
                    ftype(nfile) = 'L'
                    kolor(nfile) = kkol
                    width(nfile) = twid
                    psymb(nfile)=tsymb
                    if(tsiz.le.0.0)tsiz = 0.05
                    symsiz(nfile) = tsiz
                    ls = lgstr(tleg)
                    symleg(nfile) = tleg(1:ls)
                    nleg = nleg + 1
                    if(ls.gt.ncharleg)ncharleg=ls
                go to 3002
 4002           continue
                close(1)
            endif
        endif
c-----
c       now open the l3cmdfil and read the entries
c-----
        if(l3cmdfil .ne. ' ' )then
            inquire(file=l3cmdfil,exist=ext)
            if(ext)then
         
            
                open(1,file=l3cmdfil,status='unknown',form='formatted',
     1            access='sequential')
                rewind 1
c-----
c               this is now followed by the following items
c               FILE KOLOR WIDTH PSYMB
c               FILE file name of x-y pairs to be plotted
c               KOLOR 1 - BLACK 1000=red 1050=green 1100=blue 0 = white
c               WIDTH width of line in inches
c               PSYMB - a 2 character entry with the following meaning
c                   SQ - square
c                   TR - triangle
c                   HX - heaxgon
c                   DI - diamond
c                   CI - circle
c                   NO - no symbol - plot a line
c-----
 3003          continue
                    read(1,*,end=4003)tfile,kkol,twid,tsymb,tsiz,tleg
                    nfile = nfile + 1
                    fname(nfile) = tfile
                    ftype(nfile) = '3'
                    kolor(nfile) = kkol
                    width(nfile) = twid
                    psymb(nfile)=tsymb
                    if(tsiz.le.0.0)tsiz = 0.05
                    symsiz(nfile) = tsiz
                    ls = lgstr(tleg)
                    symleg(nfile) = tleg(1:ls)
                    nleg = nleg + 1
                    if(ls.gt.ncharleg)ncharleg=ls
                go to 3003
 4003           continue
                close(1)
            endif
        endif
        return
        end

        subroutine usage(str)
        character str*(*)
        parameter (LER=0,LIN=5,LOT=6)
        write(LER,*)'genplt: ',str
        write(LER,*)'USAGE:',
     1  'genplt ',
     1  '-XMIN xmin -XMAX xmax -YMIN ymin -YMAX ymax ',
     1  '-X0 x0 -Y0 y0 -NOBOX ',
     1  '-XLIN _XLOG _YLIN YLOG ',
     1  '-XLOG -XLIN -YLOG -YLIN  ',
     1  '-TX x-title -TY y-title ',
     1  '-C cmdfil -A acmdfil -L lcmdfil -E ecmdfil -LPOS pos',
     1  '-XDOWN -NOTICX -NOTICY',
     1  '-? -h'

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
     1  '-XLEN xlen (default 5.0)  length of X-Axis'
        write(LER,*)
     1  '-YLEN ylen (default 5.0)  length of Y-Axis'
        write(LER,*)
     1  '-NOBOX     (default false) do not plot axes'
        write(LER,*)
     1  '-NOTICX    (default false) do not plot X tics and numbers'
        write(LER,*)
     1  '-NOTICY    (default false) do not plot Y tics and numbers'
        write(LER,*)
     1  '-XLIN      (default linear) X axis is linear'
        write(LER,*)
     1  '-XLOG      (default linear) X axis is logarithmic'
        write(LER,*)
     1  '-YLIN      (default linear) Y axis is linear'
        write(LER,*)
     1  '-YLOG      (default linear) Y axis is logarithmic'
        write(LER,*)
     1  '-C cmdfil  ',
     1  '   where cmdfil consists of one xy-pair file per line as'
        WRITE(LER,*)
     1  '      File Kolor Width Psymb     '
        write(LER,*)
     1  '-A acmdfil  ',
     1  '  where acmdfil consists of one xy-pair file per line as'
        WRITE(LER,*)
     1  '      File Kolor Width Psymb Size     '
        write(LER,*)
     1  '-L lcmdfil  ',
     1  '   where lcmdfil consists of one xy-pair file per line as'
        WRITE(LER,*)
     1  '      File Kolor Width Psymb Size Legend    '
        write(LER,*)
     1  '-E ecmdfil  ',
     1  '   where ecmdfil consists of one xyerr  file per line as'
        WRITE(LER,*)
     1  '      File Kolor Width Psymb Size Legend    '
        WRITE(LER,*)
     1  '      where each line of xyerr file  is x y errx erry '
        WRITE(LER,*)
     1  '-3 3cmdfil  ',
     1  '   where 3cmdfil consists of one xypen file per line as'
        WRITE(LER,*)
     1  '      File3 Kolor Width Psymb Size Legend    '
        WRITE(LER,*)
        write(LER,*)
     1  '    File: file name of x-y pairs to be plotted'
        WRITE(LER,*)
     1  '        with the File and Psymb enclosed in single quotes'
        write(LER,*)
     1  '    Kolor: (integer)1=black,1000=red,1050=green,1100=blue'
        write(LER,*)
     2  '        0=white or any CALPLOT color'
        write(LER,*)
     1  '    Width: width of line/or symbol border in inches'
        write(LER,*)
     1  '    Psymb: - a quoted 2 character entry for the following '
        write(LER,*)
     1  '           SQ - square'
        write(LER,*)
     1  '           TR - triangle'
        write(LER,*)
     1  '           HX - heaxgon'
        write(LER,*)
     1  '           DI - diamond'
        write(LER,*)
     1  '           CI - circle'
        write(LER,*)
     1  '           DA - dashed line'
        write(LER,*)
     1  '           DO - dotted line'
        write(LER,*)
     1  '           NO - no symbol - plot a line'
        write(LER,*)
     1  '    Size: - Symbol size '
        write(LER,*)
     1  '    Legend: quoted string plotted with symbol'
        WRITE(LER,*)
     1  '        as legend in position defined by -LPOS flag'
        write(LER,*)
     1  '    File3: file name of x-y-pen_color triplets to be plotted'
        WRITE(LER,*)
     1  '        with the File and Psymb enclosed in single quotes'
        write(LER,*)
     1  '-P pcmdfil  ',
     1  '   where pcmdfil consists of one xy-pair polygon  per line as'
        WRITE(LER,*)
     1  '      Polygon_file  Kolor Fill     '
        WRITE(LER,*)
     1  '      Fill = 1 is shaded with the color'
        WRITE(LER,*)
     1  '      Fill = 0 is outline only  color'

        write(LER,*)
     1  '-LPOS pos   (default none) Must be enclosed in single quotes'
        write(LER,*)
     1  '     where the legend position, pos, is TL, ML, BL, TR, MR, BR'
        write(LER,*)
     1  '-TX title-x (default none) Must be enclosed in single quotes'
        write(LER,*)
     1  '-TY title-y (default none) Must be enclosed in single quotes'
        write(LER,*)
     1  '-XDOWN      (default false) x-axis is positive down'
        write(LER,*)
     1  '-?         (default false) online help'
        write(LER,*)
     1  '-h         (default false) online help'
        stop 
        end
