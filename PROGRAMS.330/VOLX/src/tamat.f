c-----
c       create one file with the system of linear equations for the
c       joint inversion.
c       Read the 
c       CHANGES
c
c       02 10 2000 - added a sumav1 and a sumag1 
c           for goodness of fit information
c           This will be the sum of | residual |
c           this will be an indicator of goodness of fit
c
c       we will compute average SUM | res |
c-----
        subroutine tamat(id,nfilt,inarg,wref,invcsl,nf1,invdep,
     1      lstinv,pval,sigv,sigtp,sigts,
     1              nurftn)
c-----
c       a general purpose reformatting file
c-----
c
c       list partial derivatives or computed dispersion values
c-----
        parameter(LER=0,LIN=5,LOT=6)
        integer NL, NL2, NLAY
        parameter(NL=200,NLAY=200,NL2=NL+NL)
c-----
c       LIN - unit for FORTRAN read from terminal
c       LOT - unit for FORTRAN write to terminal
c       LER - unit for FORTRAN error output to terminal
c       NL  - number of layers in model
c       NL2 - number of columns in model (first NL2/2 are
c           - velocity parameters, second NL2/2 are Q values)
c-----
        common/ctrl/numa,d(NL),a(NL),b(NL),r(NL),rat(NL),dd(NL2),x(NL2),
     $      h(NL2),u(NL2),c(NL2),v(NL2,NL2),qbinv(NL),qainv(NL),
     2      wc(NL2)
        logical wc
        character*50 names
        character ostr*12

        integer nlow,nhgh,mm,mm2,iid
        real sum2o, sum2p, sum2r,redv

        common/isomod/dl(NLAY),va(NLAY),vb(NLAY),rho(NLAY),
     1      qa(NLAY),qb(NLAY),etap(NLAY),etas(NLAY), 
     2      frefp(NLAY), frefs(NLAY)
        common/depref/refdep
        common/modtit/title
        character*80 title

        character kstnm*8
        integer nzyear, nzjday, nzhour, nzmin, nzmon, nzday
c-----
c       COMMON FOR WEIGHTING
c-----
        common/datvar/numr, sumrr, sumrd, rnorm,
     1                nums, sumss, sumsd, snorm
c-----
c       dummy array for reading partial derivatives
c-----
        real ddd(NL2)
c-----
c       The reduction of variance parameter is defined
c       as RedVar = SUM (obs - pre)^2 / SUM obs ^2
c-----
c       Now estimate the Standard deviation
c-----
        write(6,'(a,f10.5,a,i5,a,1x,a)')
     1      'Percent of Signal Power Fit (Time)     :',
     1      (1.0 - sumrr/sumrd)*100.0, '% for',numr, ' Travel time',
     2      ' OBS'
C        if(id.ne.3)then
C        write(6,'(a,f10.5,a,i5,a)')
C     1      'Percent of Signal Power Fit (Disp)     :',
C     1      (1.0 - sumss/sumsd)*100.0, '% for',nums  , ' SW Obs'
C        endif
c-----
c       To create the weighting for the joint functional,
c       note that we have already divided the residual by the
c       std err of observation
c-----
C        p = pval
c-----
c    to make the p-wave affect the inversion more, especially when
c    the dispersion and travel time rows are of different magnitude,
c    we  adjust the P  according to the magnitudes
c-----
C        facs = sqrt(    p  ) / sqrt( real(nums))
C        fact = sqrt(1.0 - p) / sqrt( real(numr))
C         p = pval * ( rnorm/(snorm+rnorm))
C         facs = sqrt(p)
C         fact = sqrt(1.0 - p )
C        WRITE( 6,*)'p,facs,fact:',p,facs,fact,nums,numr
c-----
c       open the output file for srfinv96
c-----
           open(3,file='tmpsrfi.09',form='unformatted',
     1          access='sequential')
           rewind 3
C       sumn = 0
C       sumx = 0.0
C       sumxx = 0.0
c-----
c       get the current value of dd and save it for the inversion
c-----
        call ddcur(m,m2)
c-----
c       read the surface wave information
c-----
C       open(2,file='tmpmrgs.9',form='unformatted',access='sequential')
c-----
c	two pass get maximum fow sum for the partials
c-----
C        rowsum = 0
C        rewind 2
C            read(2)m2,nfilt
C 1001   continue
C            read(2,end=1998) (ddd(j),j=1,m2),dy
C            rows = 0.0
C            do 1101 j=1,m2
C                rows = rows + abs(ddd(j))
C 1101       continue
C            if(rows.gt.rowsum)rowsum = rows
C        go to 1001
C 1998   continue
C        write(0,*)'DISP: facs=',facs,' rowsum=',rowsum, ' snorm:',snorm
C        IF(ROWSUM.EQ.0.0)ROWSUM = 1.0
C        FACS = FACS/ROWSUM
C        WRITE(0,*)'    FACS=',FACS
c-----
c       now set up the inversion
c-----
C       rewind 2
C           read(2)m2,nfilt
C           write(3)m2,nfilt
C           write(3)(dd(i),i=1,m2)
C           write(3)(wc(i),i=1,m2)
C1000   continue
C           read(2,end=1999) (ddd(j),j=1,m2),dy
C           do 1100 j=1,m2
C               ddd(j) = ddd(j) * facs
C1100       continue
C           dy = dy * facs
C           write(3)(ddd(j),j=1,m2),dy
C        WRITE(6,*)'SHAMAT: S:',(ddd(j),j=1,m2),dy

C       go to 1000
C1999   continue
c-----
c       read the travel time information
c-----
        if(id.eq.2 .or. id.eq.4)then
        open(2,file='tmpmrgt.9',form='unformatted',access='sequential')
c-----
c	two pass get maximum fow sum for the partials
c-----
C        rowsum = 0
C        rewind 2
C            read(2)m2,nfilt
C 2001   continue
C            read(2,end=2998) (ddd(j),j=1,m2),dy
CC        WRITE(6,*)'SHAMAT: T:',(ddd(j),j=1,m2),dy
C            rows = 0.0
C            do 2101 j=1,m2
C                rows = rows + abs(ddd(j))
C 2101       continue
C            if(rows.gt.rowsum)rowsum = rows
C        go to 2001
C 2998   continue
C        write(0,*)'TIME: fact=',fact,' rowsum=',rowsum, ' rnorm:',rnorm
C        IF(ROWSUM.EQ.0.0)ROWSUM = 1.0
C        FACT = FACT/ROWSUM
C        WRITE(0,*)'    FACT=',FACT
c-----
c-----
c       now set up the inversion
c-----
        rewind 2
c-----
c           do not have to write this since we have just done it
c-----
            read(2)m2,nfilt
c-----
c       now write out the partials and residuals for the inversion
c-----
        FACT=1.0
        write(3)m2,nfilt
        write(3)(dd(i),i=1,m2)
        write(3)(wc(i),i=1,m2)

 2000   continue
            read(2,end=2999) (ddd(j),j=1,m2),dy
            do 2100 j=1,m2
                ddd(j) = ddd(j) *fact
 2100       continue
            dy = dy * fact
            write(3)(ddd(j),j=1,m2),dy
CRBH        WRITE(6,*)'SHAMAT: T:',(ddd(j),j=1,m2),dy
        go to 2000
 2999   continue
        close(2)
        endif
        close(3)
        return
        end

