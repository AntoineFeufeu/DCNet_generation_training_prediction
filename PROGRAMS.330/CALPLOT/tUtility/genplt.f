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
c-----

      integer LIN, LOT
      parameter (LIN=5,LOT=6)

      real x0, y0, xlen, ylen
      real xmin, xmax, ymin, ymax

      logical doxlin, doylin
      character titlex*80
      character titley*80


      integer MXFILE
      parameter(MXFILE=100)

      integer nfile
      character*100 fname(MXFILE)
      real width(MXFILE)
      integer kolor(MXFILE)
      character*2 psymb(MXFILE)
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
     1  xmin,xmax,ymin,ymax,nfile,fname,width,kolor,psymb,symsiz,
     1  lpos,symleg,nleg,ncharleg)
c-----
c     define legend position according to lops
c     adjust sixe to fit
c-----
      write(0,*)nleg
      if(nleg.gt.0)then
          ldy = 0.5*ylen/nleg
      else
          ldy = 0.14
      endif
      if(ldy.gt.0.14)ldy = 0.14
      if(lpos(1:1).eq.'T')then
            ly0 = y0 + ylen - 0.25 
      else if(lpos(1:1).eq.'M')then
	    ly0 = y0 + 0.5*ylen +0.5*nleg*ldy
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
            lx0 = x0 + 0.50*xlen
      else if(lpos(2:2).eq.'R')then
            lx0 = x0 + xlen - 0.25 -(ncharleg+3)*lht
      endif
      WRITE(0,*)lx0,ly0,lht,ldy

      call pinitf('GENPLT.PLT')

      llx = lgstr(titlex)
      WRITE(6,*)'llx=',llx
      sizex = 0.15
      if(doxlin)then
        ticup = .true.
        labtop = .false.
        dopow = .true.
        call dolinx(x0,y0,xlen,xmax,xmin,
     1      sizex,ticup,labtop,dopow,llx,titlex)
        ticup = .false.
        labtop = .false.
        dopow = .false.
        call dolinx(x0,y0+ylen,xlen,xmax,xmin,
     1      sizex,ticup,labtop,dopow,llx,titlex)
      else
        ticup = .true.
        labtop = .false.
        dopow = .true.
        call dologx(x0,y0,xlen,xmax,xmin,
     1      sizex,ticup,labtop,dopow,llx,titlex)
        ticup = .false.
        labtop = .false.
        dopow = .false.
        call dologx(x0,y0+ylen,xlen,xmax,xmin,
     1      sizex,ticup,labtop,dopow,llx,titlex)
      endif
      lly = lgstr(titley)
      WRITE(6,*)'lly=',lly
      sizey = 0.15
      if(doylin)then
        ticlft = .false.
        lablft = .true.
        dopow = .true.
        call doliny(x0,y0,ylen,ymax,ymin,
     1      sizey,tilft,lablft,dopow,lly,titley)
        ticlft = .true.
        lablft = .false.
        dopow = .false.
        call doliny(x0+xlen,y0,ylen,ymax,ymin,
     1      sizey,ticlft,lablft,dopow,lly,titley)
      else
        ticlft = .false.
        lablft = .true.
        dopow = .true.
        call dology(x0,y0,ylen,ymax,ymin,
     1      sizey,tilft,lablft,dopow,lly,titley)
        ticlft = .true.
        lablft = .false.
        dopow = .false.
        call dology(x0+xlen,y0,ylen,ymax,ymin,
     1      sizey,ticlft,lablft,dopow,lly,titley)
      endif

      call gclip('ON',x0,y0,x0+xlen,y0+ylen)
      do nf = 1, nfile
      WRITE(6,*)nf,nfile,fname(nf),kolor(nf),width(nf),psymb(nf),
     1     symsiz(nf)
        open(1,file=fname(nf),access='sequential',
     1   form = 'formatted',status='unknown')
        rewind 1
C        if(psymb(nf).eq.'NO')then
           call gwidth(width(nf))
C        else 
C           call gwidth(0.10*symsiz(nf))
C        endif

        ipen = 3
 1000   continue
        call newpen(kolor(nf))
        read(1,*,end=9999)x,y
C        WRITE(6,*)x,y
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
                if(psymb(nf).eq.'NO')then
                    call plot(xx,yy,ipen)
                    ipen = 2
                else
                    call newpen(kolor(nf))
                    call fillit(psymb(nf),symsiz(nf),xx,yy)
                    if(width(nf).gt.0.0)then
                    call newpen(1)
                    call curvit(psymb(nf),symsiz(nf),xx,yy)
                    endif
                endif
        go to 1000
 9999   continue
        ipen = 3
               if(psymb(nf).eq.'NO')then
                    call plot(xx,yy,ipen)
               endif
        close (1)
c-----
c       plot legend
c-----
        ls = lgstr(symleg(nf))
        WRITE(0,*)symleg(nf)(1:ls),lx0,ly0,lht,ldy
        if(symleg(nf) .ne. ' ' )then
            lxx = lx0
            lyy = ly0 
                if(psymb(nf).eq.'NO')then
                    call plot(lxx,lyy+0.5*lht,3)
                    call plot(lxx+lht,lyy+0.5*lht,2)
                else
                    call newpen(kolor(nf))
                    call fillit(psymb(nf),symsiz(nf),lxx,lyy+0.5*lht)
                    if(width(nf).gt.0.0)then
                    call newpen(1)
                    call curvit(psymb(nf),symsiz(nf),lxx,lyy+0.5*lht)
                    endif
                endif
            call newpen(1)
            call gleft(lxx+2.*lht,lyy,lht,symleg(nf)(1:ls),0.0)
            ly0 = ly0 - ldy
     
         endif
      enddo
      call gclip('OFF',x0,y0,x0+xlen,y0+ylen)
c-----
c     reset the state
c-----
      call gwidth(0.0)
      call newpen(1)
      call pend()
      end

      subroutine gcmdln(x0,y0,xlen,ylen,titlex,titley,doxlin,doylin,
     1  xmin,xmax,ymin,ymax,nfile,fname,width,kolor,psymb,symsiz,
     1  lpos,symleg,nleg,ncharleg)
      real x0,y0,xlen,ylen,xmin,xmax,ymin,ymax
      logical doxlin, doylin
      character titlex*80
      character titley*80

      integer MXFILE
      parameter(MXFILE=10)

      integer nfile
      character*100 fname(MXFILE)
      real width(MXFILE)
      integer kolor(MXFILE)
      character*2 psymb(MXFILE)
      real symsiz(MXFILE)
      character lpos*4
      character*100 symleg(MXFILE)
      integer nleg,ncharleg

      
 
      character*80 names
      integer mnmarg
      integer i
      integer ls, ln
      character cmdfil*80
      character acmdfil*80
      character lcmdfil*80
      character tfile*80, tsymb*2, tleg*80
      logical ext

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
      lcmdfil = ' '
      nleg = 0
      lpos = ''
      ncharleg = 0

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
c       now open the cmdfil and read the entries
c-----
        if(cmdfil .eq. ' ' .and. acmdfil .eq. ' '
     1     .and. lcmdfil .eq. ' ' )then
            call usage('No cmdfil or acmdfil')
        endif
        nfile = 0
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
c                   NO - no symbol - plot a line
c-----
 3000          continue
                    read(1,*,end=4000)tfile,kkol,twid,tsymb
                    nfile = nfile + 1
                    fname(nfile) = tfile
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
     1  '-C cmdfil -A acmdfil -L lcmdfil -LPOS pos',
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
     1  '-XLEN xlen (default 6.0)  length of X-Axis'
        write(LER,*)
     1  '-YLEN ylen (default 6.0)  length of Y-Axis'
        write(LER,*)
     1  '-NOBOX     (default false) do not plot axes'
        write(LER,*)
     1  '-XLIN      (default linear) X axis is linear'
        write(LER,*)
     1  '-XLOG      (default linear) X axis is logarithmic'
        write(LER,*)
     1  '-YLIN      (default linear) Y axis is linear'
        write(LER,*)
     1  '-YLOG      (default linear) Y axis is logarithmic'
        write(LER,*)
     1  '-C cmdfil  '
        WRITE(LER,*)
     1  '   where cmdfil consists of one xy-pair file per line as'
        WRITE(LER,*)
     1  '      File Kolor Width Psymb     '
        write(LER,*)
     1  '-A acmdfil  '
        WRITE(LER,*)
     1  '  where acmdfil consists of one xy-pair file per line as'
        WRITE(LER,*)
     1  '      File Kolor Width Psymb Size     '
        write(LER,*)
     1  '-L lcmdfil  '
        WRITE(LER,*)
     1  '   where lcmdfil consists of one xy-pair file per line as'
        WRITE(LER,*)
     1  '      File Kolor Width Psymb Size Legend    '
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
     1  '           NO - no symbol - plot a line'
        write(LER,*)
     1  '    Size: - Symbol size '
        write(LER,*)
     1  '    Legend: quoted string plotted with symbol'
        WRITE(LER,*)
     1  '        as legend in position defined by -LPOS flag'
        write(LER,*)
     1  '-LPOS pos   (default none) Must be enclosed in single quotes'
        write(LER,*)
     1  '     where the legend position, pos, is TL, ML, BL, TR, MR, BR'
        write(LER,*)
     1  '-TX title-x (default none) Must be enclosed in single quotes'
        write(LER,*)
     1  '-TY title-y (default none) Must be enclosed in single quotes'
        write(LER,*)
     1  '-?         (default false) online help'
        write(LER,*)
     1  '-h         (default false) online help'
        stop 
        end
