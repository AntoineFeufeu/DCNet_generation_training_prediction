        subroutine outphvctrl(xlin,ylin,sacfile,lnlg)
        logical xlin, ylin
        character*256 sacfile
        character *3 lnlg(2)

        common/pltinf/gxl(3),gyl(3),gxh(3),gyh(3),
     1      axl(3),ayl(3),axh(3),ayh(3),ixlnlg(3),iylnlg(3)
            real*4 gxl,gyl,gxh,gyh,
     1              axl,ayl,axh,ayh
            integer*4 ixlnlg,iylnlg
        common/pltcnt/laxper
            logical laxper
        common/discnt/kount, lorr, pcount, dophvel
            integer *4 kount, lorr, pcount, dophvel
        integer MAXFRQ
        parameter (MAXFRQ=100)
        common /fval/ wn, nper, wnin, mper
            real*4 wn(MAXFRQ), wnin(MAXFRQ)
            integer*4 nper, mper
        character outstr*256
        common/tit/npts,az,dt,dist,n,a0,a1
            integer*4 npts
            real*4 az,dt,dist,a0,a1

        integer ls
c-----
c       open file containing plot control information
c-----
c 249
cXAXIS-PERIOD
c 5.100  1.500  9.100  5.500   14.0  2.00      75.0  5.00    lin lin PHV96.PLT
c 32
c   75.00000       70.00000       65.00000       60.00000       55.00000
c   50.00000       48.00000       46.00000       44.00000       42.00000
c   40.00000       38.00000       36.00000       34.00000       32.00000
c   30.00000       29.00000       28.00000       27.00000       26.00000
c   25.00000       24.00000       23.00000       22.00000       21.00000
c   20.00000       19.00000       18.00000       17.00000       16.00000
c   15.00000       14.00000                                                   
c-----
c Value     Name    Meaning
c 249       kount   Number of dispersion points in mft96.dsp 
c XAXIS-PERIOD  strong  Key word for plot or XAXIS-FREQUENCY
c one lines to describe the  image on the file PHV96.PLT
c   each line has the physical coordinates of the plot (XL,YL) (XH,YH)
c   the physical values at those corners  
c   (XV,YV) at lower left and (XV,YV) at upper right
c   two strongs to describe the plotted axes - for the first figure, 
c   the spectra plot
c       the plot is log-log while for the dispersion curve (center) 
c   and trace plot(thrid) the
c       axes are linear-linear. Finally the graphics file that 
c   contains these images is
c       PHV96.PLT for all three plots
c  32   nper    number of period values
c this is now followed by 32 periods for the plot
c       
c-----
        open(3,file='mft96.ctl',access='sequential',
     1      form='formatted',status='unknown')
        rewind 3
        write(3,*)kount
        if(laxper)then
            write(3,14)
        else
            write(3,13)
        endif
        do 11 i=1,3
        write(3,10) gxl(i),gyl(i),gxh(i),gyh(i),
     1      axl(i),ayl(i),axh(i),ayh(i),
     2      lnlg(ixlnlg(i)),lnlg(iylnlg(i))
   11   continue
   10   format(4f7.3,4g11.3,' ',a3,' ',a3,' PHV96.PLT')
   13   format('XAXIS-FREQUENCY')
   14   format('XAXIS-PERIOD')
        write(3,*)nper
        do i=1,nper
            wnin(i) = 1.0/wn(i)
        enddo
        write(3,'(5g15.7)')(wnin(i),i=1,nper)
c-----
c       close the plot parameter file
c-----
        close (3)
        if(lorr.gt.0)then
c-----
c       open file for invoking  DPEGN96 control information
c-----
C-----
C       DOS
C-----
C       open(3,file='PHV96CMP.BAT',access='sequential',
C     1     form='formatted',status='unknown')
C       rewind 3
C-----
c-----
c       UNIX
c-----
        open(3,file='PHV96CMP',access='sequential',
     1      form='formatted',status='unknown')
        rewind 3
        write(3,'(a)')'#!/bin/sh'
        write(3,'(a)')' '
        lssf = lgstr(sacfile)
        write(3,17)sacfile(1:lssf), a0
        write(3,'(a)')' '
        outstr = ' '
   17 format('# Data file = ',a/'# alpha = ',f10.3)
   15   format('sdpegn96 -X0 ',f5.2,' -Y0 ',f5.2,' -XLEN ',
     1   f5.2, ' -YLEN ',f5.2,' -XMIN ',g10.3,
     1   ' -XMAX ',g10.3, ' -YMIN ',f5.2,' -YMAX ',f5.2)
        write(outstr,15)gxl(2),gyl(2),gxh(2)-gxl(2),gyh(2)-gyl(2),
     1      axl(2), axh(2), ayl(2), ayh(2)
        ls = lgstr(outstr)
        if(laxper)then
            outstr(ls+1:ls+7) = ' -PER  '
        else
            outstr(ls+1:ls+7) = ' -FREQ '
        endif
        ls = ls + 7
        if(lorr.eq.1)then
            outstr(ls+1:ls+13) = '-L -C -NOBOX '
        else if(lorr.eq.2)then
            outstr(ls+1:ls+13) = '-R -C -NOBOX '
        endif
        ls = ls + 13
        if(xlin)then
            outstr(ls+1:ls+6)='-XLIN '
        else
            outstr(ls+1:ls+6)='-XLOG '
        endif
        ls = ls + 6
        if(ylin)then
            outstr(ls+1:ls+6)='-YLIN '
        else
            outstr(ls+1:ls+6)='-YLOG '
        endif
        ls = ls + 6
        write(3,'(a)')outstr(1:ls)

c-----
c       close the SDPEGN96 invocation file
c-----
        close (3)
        endif
        return
        end
