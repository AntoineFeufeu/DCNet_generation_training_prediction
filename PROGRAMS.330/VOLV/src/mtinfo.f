      program mtinfo
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME V                                                       c
c                                                                     c
c      PROGRAM: MTDEC                                                 c
c                                                                     c
c      COPYRIGHT 2016 R. B. Herrmann and M. L.Jost
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
c     Changes:
c      06 OCT 2016 - in subroutine mteig, the array
c         real*8 z(3,3) was not defined. This cause a memory
c         problem
c
c     31 OCT 2016 - error in subroutine tpdss when computing rake
c          BAD if(rake0.lt.-180.)rake0=rake0+180
c      CORRECT if(rake0.lt.-180.)rake0=rake0+360.
c     22 DEC 2020 - error in output of printer plot in file mt.msg
c          The M12 element was not printed the corrected code is
c          write(FID,'(a,1pe9.2)')'       Mxy   ', Mij(1,2)
c          write(FID,'(a,1pe9.2)')'       Myy   ', Mij(2,2)

c    04 JAN 2021 - change all do XXX yyy ... XXX a=b 
c                  to   do yy ... a=b .. enddo
c           to be compatible with gfortran 9.3
c    15 SEP 2021 - add lune plot to the right of the beachball 
c           thanks to C. J. Ammon, Penn State.
c    06 JUN 2022 - cleaned output to be less redundant
c-----
      implicit real*8 (a-h,o-z)
      integer LER, LOT, LIN
      parameter (LER=0, LIN=5, LOT=6)
      integer index(3)
      real*8 z(3,3),ev(3),evd(3),evp(3)
      real*8 a1a1(3,3),a2a2(3,3),a3a3(3,3)
      real*8 m(3,3)
      real*8 xmom
      real bet, gam
      integer ndec(6)
      tol = 1.0e-12
      np=3
      ndata = 3
c-----
c
c     mtinfo -mc -cc -dc -vd -clvd -a -xx Mxx -yy Myy -zz Mzz -xy Mxy -xz Mxz -yz Myz
c
c     This program first converts the moment tensor
c     into eigenvalues and eigenvector.
c
c     Then this information is used to perform a 
c     general moment tensor decomposition.
c
c     First, the isotropic part is taken out of the moment tensor.
c     The null-, T-, and P-axis are determined from
c     the eigenvalues and eigenvectors of a seismic moment tensor.
c     Consequently, the output of MTEIG (VIII) is used as input.
c     The eigenvector to the largest positive eigenvalue is the T
c     axis. The eigenvector to the largest negative eigenvalue is
c     the P-axis. 
c
c    -mc The smallest eigenvalue in the
c        absolute sense is set to zero defining the major double 
c        couple. The corresponding eigenvector defines the null-axis.
c        From these axes, the focal mechanism i.e. strike, dip, and rake
c        as well as the trend and plunge of the X-, Y-, null-, T-, and P-
c        axes is determined following Herrmann (1975, see XYZTP (II)).
c        The minor couple is then constructed by letting the eigenvector
c        of the largest eigenvalue (in the absolute sense) be the null axis
c        (Jost and Herrmann, 1989).
c    -cc The deviatoric moment tensor is decomposed into a double couple and a
c        CLVD (Knopoff and Randall, 1970).
c    -dc The moment tensor is decomposed into three double couples following
c        Ben-Menahem and Singh (1981, eq.4.57). For each double couple,
c        the corresponding focal mechanism is determined (XYZTP(II)).
c    -vd The moment tensor is decomposed into three vector dipoles following
c        Ben-Menahem and Singh (1981, eq.4.55).
c  -clvd The moment tensor is decomposed into three compensated linear vector 
c        dipoles following Ben-Menahem and Singh (1981, eq.4.56).
c    -crack
c    -a  all decompositions are considered.
c-----

c-----
c     call machine dependent initialization
c-----
      call mchdep()
c-----
c     Get the command line arguments
c-----
      call gcmdln(ndec,m)
c-----
c     get the eigenvalectors of the z(3,3) moment tensor
c-----
      call showinput(ndata,xmom,ev,z,m,bet,gam)
      call showmij('Input:',ev,z,m)
         call isotr(ndata,index,tol,thresh,xmom,zmom,ev,evp,evd,z,
     1       bet,gam)
         call dyad(a1a1,a2a2,a3a3,z)
         if(ndec(1).eq.1) then
           call mjrcpl(ndata,np,index,tol,thresh,zmom,evp,evd,z)
         endif
         if(ndec(2).eq.1)then
           call dblcpl(ndata,np,tol,thresh,zmom,ev,a1a1,a2a2,a3a3)
         endif
         if(ndec(3).eq.1) then
           call vctdpl(ndata,tol,thresh,zmom,evp,evd,z)
         endif
         if(ndec(4).eq.1) then
c-----
c      3 CLVD's
c-----
           call clvdpl(ndata,index,tol,thresh,zmom,ev,z)
         endif
         if(ndec(5).eq.1) then
c-----
c      Double couple and CLVD
c-----
           call dcclvd(ndata,np,index,tol,thresh,xmom,zmom,ev,evp,evd,z)
         endif
         if(ndec(6).eq.1) then
c-----
c      Crack
c-----
           call crack(ndata,np,index,tol,thresh,zmom,ev,evp,evd,z)
         endif
      end

      subroutine showinput(ndata,xmom,ev,z,m,bet,gam)
c-----
c     ndate    I - dimension of moment tensor = 3
c     xmom     R - scale value of seismic moment
c     ev       R - array of eigenvalues
c     z        R - array of eigenvectors
c     m        R - input 3x3 moment tensor
c-----
      implicit real*8 (a-h,o-z)
      parameter (LER=0, LIN=5, LOT=6)
      real*8 z(3,3),ev(3), m(3,3)
      real*8 xmom

      real*8 plnt,stkt,plnp,stkp,plnn,stkn
      real*8 zdum(3,3)
      real*8 ev1(3)
      real bet, gam
      write(LOT,*) ' '
      write(LOT,'(a)') '  INPUT MOMENT TENSOR'
      do  i=1,ndata
           write(LOT,5)(m(i,j),j=1,3)
      enddo
    5 format(5e15.7)
c-----
c      SEISMIC MOMENT OF THE GIVEN MOMENT TENSOR
c-----
      xmom = 0.0
      do  i=1,ndata
           do  j=1,ndata
                xmom=xmom+m(i,j)*m(i,j)
                z(i,j) = m(i,j)
           enddo
      enddo
      xmom = dsqrt(xmom/2.0)
c-----
c     get eigenvalues and eigenvectors from the
c     moment tensor matric
c-----
      call gtev(zdum,m,ndata,ev,ev1,ilg,ism,z)
c-----
c     order the vectors according to ioncreasing eigenvalue
c-----
      write(LOT,'(a)') '   EIGENVALUES               EIGENVECTORS'
      do  i=1,ndata
           write(LOT,7)i,ev(i),(z(j,i),j=1,3)
      enddo
    7   format(' ',i5,1x,e11.4,'  (',e11.4,',',e11.4,',',e11.4,')')
      write(LOT,*) ' '
c-----
c     display the solution as a msg format
c-----
      call gettrpl(stkp,plnp,z(1,1),z(2,1),z(3,1))
      call gettrpl(stkn,plnn,z(1,2),z(2,2),z(3,2))
      call gettrpl(stkt,plnt,z(1,3),z(2,3),z(3,3))
      call  PrintPradiation('mt.msg',m,xmom,plnt,stkt,plnp,stkp,
     1       plnn,stkn,ev,bet,gam)
      return
      end

      subroutine isotr(ndata,index,tol,thresh,xmom,zmom,ev,evp,evd,z,
     1       bet,gam)
      implicit real*8 (a-h,o-z)
      parameter (LER=0, LIN=5, LOT=6)
      real*8 mdev(3,3)
      dimension z(3,3),ev(3),evp(3),evd(3),index(3)
      real bet, gam
      write(LOT,'(a,e15.7)') 
     $'  SEISMIC MOMENT OF THE INPUT MOMENT TENSOR: ',xmom
      write(LOT,*) ' ' 
      expl=0.
      do 1 i=1,ndata
      expl=expl+ev(i)
    1 continue
      write(LOT,'(a,e15.7)') 
     $'  ISOTROPIC COMPONENT (TRACE OF MOMENT TENSOR): ',expl
      write(LOT,*) ' '
      zmom = expl/3.0
      wmax=0.
      do 2 i=1,ndata
      if(abs(ev(i)).gt.wmax) wmax=abs(ev(i))
    2 continue
      thresh=tol*wmax
      if(expl.gt.thresh) then
      write(LOT,'(a)')
     $'  INPUT MOMENT TENSOR IS PARTIALLY DUE TO AN EXPLOSION'
      else if(expl.lt.-thresh) then
      write(LOT,'(a)')
     $'  INPUT MOMENT TENSOR IS PARTIALLY DUE TO AN IMPLOSION'
      else
      write(LOT,'(a)') 
     $'  INPUT MOMENT TENSOR IS PURELY DEVIATORIC'
      expl = 0.0
      endif
      do 3 i=1,ndata
      evd(i)=ev(i)-expl/float(ndata)
    3 continue
          call showmij('Deviatoric:',evd,z,mdev)
      write(LOT,'(a)') 
     $'  EIGENVALUES OF THE PURELY DEVIATORIC MOMENT TENSOR'
      write(LOT,6) (evd(i),i=1,ndata)
      write(LOT,*) ' '
c---- sorting the eigenvalues (principal stresses) into ascending order
      do 4 i=1,ndata
      index(i)=i
      evp(i)=abs(evd(i))
    4 continue
      do 8 j=2,ndata
      a=evp(j)
      nin=index(j)
      do 9 i=j-1,1,-1
      if(evp(i).le.a) goto 10
      evp(i+1)=evp(i)
      index(i+1)=index(i)
    9 continue
      i=0
   10 evp(i+1)=a
      index(i+1)=nin
    8 continue
      write(LOT,'(a)') '  SORTED EIGENVALUES'
      write(LOT,6) (evp(i),i=1,ndata)
c---- Dziewonski, Chou, Woodhouse,  JGR 1981 2825-2852
c---- eps=0 double couple
c---- eps=0.5 vector dipole
      eps=evp(1)/evp(ndata)
      eps1=eps*200.0
      eps2=100.0-eps1
      write(LOT,'(a)') '  EPSILON    % OF CLVD     % OF DC'
      write(LOT,11) eps, eps1, eps2
    6 format(5e15.7)
   11 format(f10.6,2f12.6)
      write(LOT,'(a)') '  LUNE PARAMETERS'
      WRITE(LOT,12)bet,gam
   12 format('  beta: ',f7.2,' gamma:',f7.2)
      return 
      end

      subroutine mjrcpl(ndata,np,index,tol,thresh,zmom,evp,evd,z)
      implicit real*8 (a-h,o-z)
      parameter (LER=0, LIN=5, LOT=6)
      dimension index(3),ev(3),evd(3),evp(3),z(3,3)
      jj=0
      do 1 i=1,ndata
      if(abs(evd(i)).lt.thresh) jj=jj+1
    1 continue
      write(LOT,*) ' '
      write(LOT,'(a)') 
     $'  --------------------------------------------------------------'
      write(LOT,*) ' '
      write(LOT,'(a)') 
     $'  DECOMPOSITION INTO MAJOR AND MINOR DOUBLE COUPLE'
      if(jj.eq.1) then
      write(LOT,*) ' '
      write(LOT,'(a)') 
     $'  DEVIATORIC MOMENT TENSOR IS DUE TO A DOUBLE COUPLE'
      call xyztp2(evd,z,np,ndata,tol,dip,stk,slp)
      else if(jj.eq.2) then
      write(LOT,*) ' '
      write(LOT,'(a)') 
     $'  DEVIATORIC MOMENT TENSOR IS DUE TO A VECTOR DIPOLE'
      else if(jj.eq.0) then
      write(LOT,*) ' '
      write(LOT,'(a)') 
     $'  SEISMIC MOMENT OF ISOTROPIC PART, MAJOR, AND MINOR COUPLE'
      write(LOT,5) zmom, evp(3),evp(1)
c---- MAJOR COUPLE
c---- Dziewonski, Chou, Woodhouse, 1981 p2842
      ev(index(1)) = 0.0
      ev(index(2)) = -evd(index(3))
      ev(index(3)) = evd(index(3))
      write(LOT,*) ' '
      write(LOT,'(a)') 
     $'  EIGENVALUES AND EIGENVECTORS OF THE MAJOR DC'
      do 3 i=1,ndata
      write(LOT,4)i, ev(i),(z(j,i),j=1,ndata)
    3 continue
    4   format(' ',i5,1x,e11.4,'  (',e11.4,',',e11.4,',',e11.4,')')
    5 format(5e15.7)
      call xyztp2(ev,z,np,ndata,tol,dip,stk,slp)
c---- note : the smallest eigenvalue (abs sense) = 0
c---- Minor double couple
      ev(index(3)) = 0.0
      ev(index(2)) = -evd(index(1))
      ev(index(1)) = evd(index(1))
      write(LOT,*) ' '
      write(LOT,'(a)') 
     $'  EIGENVALUES AND EIGENVECTORS OF THE MINOR DC'
      do 7 i=1,ndata
      write(LOT,4)i, ev(i),(z(j,i),j=1,ndata)
    7 continue
      call xyztp2(ev,z,np,ndata,tol,dip,stk,slp)
      endif
      return
      end

      subroutine dblcpl(ndata,np,tol,thresh,xmom,ev,a1a1,a2a2,a3a3)
      implicit real*8 (a-h,o-z)
      integer LER, LIN, LOT
      parameter (LER=0, LIN=5, LOT=6)
      real*8 ev(3),evd(3),evp(3)
      real*8 a1a1(3,3),a2a2(3,3),a3a3(3,3)
      real*8 emt1(3,3),emt2(3,3),emt3(3,3)
 
      integer iset3, j


      evp(1) = (ev(1)-ev(2))/3.0
      evp(2) = (ev(2)-ev(3))/3.0
      evp(3) = (ev(3)-ev(1))/3.0
      iset1 = 0
      iset2 = 0
      iset3 = 0
      do  i=1,ndata
         do  j=1,ndata
           if(evp(1).gt.0.0) then
                emt1(i,j) = a1a1(i,j) - a2a2(i,j)
                else
                emt1(i,j) = -(a1a1(i,j) - a2a2(i,j))
                iset1 = 1
           endif
           if(evp(2).gt.0.0) then
                emt2(i,j) = a2a2(i,j) - a3a3(i,j)
                else
                emt2(i,j) = -(a2a2(i,j) - a3a3(i,j))
                iset2 = 1
           endif
           if(evp(3).gt.0.0) then
                emt3(i,j) = a3a3(i,j) - a1a1(i,j)
                else
                emt3(i,j) = -(a3a3(i,j) - a1a1(i,j))
                iset3 = 1
           endif
         enddo
      enddo
      if(iset1.eq.1) evp(1) = -evp(1)
      if(iset2.eq.1) evp(2) = -evp(2)
      if(iset3.eq.1) evp(3) = -evp(3)
      write(LOT,*) ' '
      write(LOT,'(a)') 
     $'  --------------------------------------------------------------'
      write(LOT,*) ' '
      write(LOT,'(a)') '  DECOMPOSITION INTO THREE DOUBLE COUPLES'
      write(LOT,*) ' '
      write(LOT,'(a)') 
     $'  SEISMIC MOMENTS OF ISOTROPIC PART, EMT1, EMT2, EMT3'
      write(LOT,2) xmom,(evp(i),i=1,3)
    2 format(5e15.7)
      write(LOT,*) ' '
      write(LOT,'(a)') '  EMT1'
        call mteig(ndata,np,evd,emt1)
        call xyztp2(evd,emt1,np,ndata,tol,dip,stk,slp)
      write(LOT,*) ' '
      write(LOT,'(a)') '  EMT2'
        call mteig(ndata,np,evd,emt2)
        call xyztp2(evd,emt2,np,ndata,tol,dip,stk,slp)
      write(LOT,*) ' '
      write(LOT,'(a)') '  EMT3'
        call mteig(ndata,np,evd,emt3)
        call xyztp2(evd,emt3,np,ndata,tol,dip,stk,slp)
      return
      end

      subroutine vctdpl(ndata,tol,thresh,zmom,evp,evd,z)
      implicit real*8 (a-h,o-z)
      integer LER, LOT, LIN
      parameter (LER=0, LIN=5, LOT=6)
      dimension evp(3),evd(3),z(3,3)
      write(LOT,*) ' '
      write(LOT,'(a)') 
     $'  --------------------------------------------------------------'
      write(LOT,*) ' '
      write(LOT,'(a)') '  DECOMPOSITION INTO THREE VECTOR DIPOLES'
      write(LOT,*) ' '
      write(LOT,'(a)') 
     $'  SEISMIC MOMENTS OF ISOTROPIC PART, EMT1, EMT2, EMT3'
      write(LOT,2) zmom,(evp(i),i=1,3)
      write(LOT,*) ' '
      write(LOT,'(a)') '            DIPOLE AXES'
      do  i=1,ndata
         do  k=1,ndata
            test = abs(abs(evd(k))-evp(i))
            if(test.lt.thresh) then
               write(LOT,7) (z(j,k),j=1,3)
            endif
         enddo
      enddo
    2 format(5e15.7)
    7 format(5f11.7)
      return
      end

      subroutine clvdpl(ndata,index,tol,thresh,zmom,ev,z)
      implicit real*8 (a-h,o-z)
      integer LER, LOT, LIN
      parameter (LER=0, LIN=5, LOT=6)
      dimension ev(3),evp(3),z(3,3)
      dimension index(3),iev(3)
      evp(1) = abs(ev(1)/3.0) * sqrt(3.0)
      evp(2) = abs(ev(2)/3.0) * sqrt(3.0)
      evp(3) = abs(ev(3)/3.0) * sqrt(3.0)
      iev(1) = -1
      iev(2) = -1
      iev(3) = -1
      write(LOT,*) ' '
      write(LOT,'(a)') 
     $'  --------------------------------------------------------------'
      write(LOT,*) ' '
      write(LOT,'(a)') 
     $'  DECOMPOSITION INTO THREE COMPENSATED LINEAR VECTOR DIPOLES'
      write(LOT,*) ' '
      write(LOT,'(a)') 
     $'  SEISMIC MOMENTS OF ISOTROPIC PART, EMT1, EMT2, EMT3'
      write(LOT,2) zmom,(evp(index(i)),i=1,3)
      do 1 k=1,ndata
      iev(index(k)) = 2
      write(LOT,*) ' '
      write(LOT,'(a)') '  EIGENVALUE          AXES OF CLVD'
      do 3 i=1,ndata
      write(LOT,4) iev(i),(z(j,i),j=1,3)
    3 continue
      iev(index(k)) = -1
    1 continue
    2 format(5e15.7)
    4 format(i9,2x,5f11.7)
      return
      end

      subroutine dcclvd
     $(ndata,np,index,tol,thresh,xmom,zmom,ev,evp,evd,z)
      parameter (LER=0, LIN=5, LOT=6)
      implicit real*8 (a-h,o-z)
      dimension z(3,3),ev(3),evd(3),evp(3)
      dimension index(3)
      real*8 iev(3)
      real*8 eviso(3), ziso(3,3)
      real*8 m(3,3), miso(3,3),mdc(3,3),mclvd(3,3)
      jj=0
      do 1 i=1,ndata
      if(abs(evd(i)).lt.thresh) jj=jj+1
    1 continue
      write(LOT,*) ' '
      write(LOT,'(a)') 
     $'  --------------------------------------------------------------'
      write(LOT,*) ' '
      write(LOT,'(a)') 
     $'  DECOMPOSITION INTO A DOUBLE COUPLE AND A CLVD'
      if(jj.eq.1) then
          write(LOT,*) ' '
          write(LOT,'(a)') 
     $    '  DEVIATORIC MOMENT TENSOR IS DUE TO A DOUBLE COUPLE'
          call xyztp2(evd,z,np,ndata,tol,dip,stk,slp)
      else if(jj.eq.2) then
          write(LOT,*) ' '
          write(LOT,'(a)') 
     $    '  DEVIATORIC MOMENT TENSOR IS DUE TO A VECTOR DIPOLE'
      else if(jj.eq.0) then
          write(LOT,*) ' '
          write(LOT,'(a)') 
     $'  SEISMIC MOMENTS OF ISOTROPIC PART, DC, CLVD'
          eps=evp(1)/evp(3)
          amj = evp(3) * (1.0 - 2.0 * eps)
          bmj = evp(3) * eps
          write(LOT,5) zmom,amj,bmj
          do i=1,3
             do j=1,3
                ziso(i,j) = 0.0
             enddo
             ziso(i,i) = 1.0
             eviso(i) = zmom
          enddo
          call showmij('Isotropic:',eviso,ziso,miso)
c---- MAJOR COUPLE
          ev(index(1)) = 0.0
          ev(index(2)) = -evd(index(3))
          ev(index(3)) = evd(index(3))
          call showmij('Double couple:',ev,z,mdc)
          write(LOT,*) ' '
          write(LOT,'(a)') '  EIGENVALUES AND EIGENVECTORS OF THE DC'
          do i=1,ndata
              write(LOT,4) ev(i),(z(j,i),j=1,ndata)
          enddo
          call xyztp2(ev,z,np,ndata,tol,dip,stk,slp)
          iev(1) = -bmj
          iev(2) = -bmj
          iev(3) = -bmj
          iev(index(3)) = 2*bmj
          write(LOT,*) ' '
          call showmij('CLVD:',iev,z,mclvd)
          write(LOT,'(a)') '  EIGENVALUE          AXES OF CLVD'
          do i=1,ndata
              write(LOT,7) iev(i),(z(j,i),j=1,3)
          enddo
c-----
c         add the components to recover the initial moment tensor
c-----
          do i=1,3
            do j=1,3
              m(i,j) = miso(i,j) + mdc(i,j) + mclvd(i,j)
            enddo
          enddo
          write(6,*)'Original moment tensor'
          do  i=1,3
             write(LOT,9) (m(i,j),j=1,3)
          enddo
    4     format(e15.7,2x,5f11.7)
    5     format(5e15.7)
    7     format(e15.7,2x,5f11.7)
    9 format(10x,3e12.5)
      endif
      return
      end

      subroutine showmij(str,e,ev,mij)
c-----
c      | a1 b1 c1 | | l1       | a1 a2 a3 | 
c      | a2 b2 c2 | |    l2    | b1 b2 b3 |  =
c      | a3 b3 c3 | |       l3 | c1 c2 c3 |
c
c      | a1 b1 c1 | | l1*a1 l1*a2 l1*a3 |
c      | a2 b2 c2 | | l2*b1 l2*b2 l2*b3 |  =
c      | a3 b3 c3 | | l3*c1 l3*c2 l3*c3 |
c
c      | a1 b1 c1 | | l1*a1 l1*a2 l1*a3 |
c      | a2 b2 c2 | | l2*b1 l2*b2 l2*b3 |  =
c      | a3 b3 c3 | | l3*c1 l3*c2 l3*c3 |
c
      parameter (LER=0, LIN=5, LOT=6)
      character str*(*)
      real*8 e(3), ev(3,3)
      real*8 mij(3,3)
      real*8 xmom
      real xl(3), tmp, tmp1
c     A X = X L
c     A   = X L X^T
      do i=1,3
       do j = 1,3
         mij(i,j) = 0.0
         do k=1,3
            mij(i,j) = mij(i,j) + ev(i,k)*e(k)*ev(j,k)
         enddo
       enddo
      enddo
      xmom = 0.0
      do i=1,3
        do j=1,3
           xmom = xmom + mij(i,j)**2
        enddo
      enddo
      do i=1,3
        xl(i) = abs(e(i))
      enddo
c-----
c     get sum of two largest eigenvalues
c----
      tmp = amin1(xl(1), xl(2), xl(3) )

      tmp1 = xl(1) + xl(2) + xl(3) - tmp
      tmp1 = tmp1 / 2.
      
       
     
      xmom = sqrt(xmom/2.)
      write(6,*)str, ' moment=sqrt(sum mij**2/2)=',xmom
      write(6,*)'     (|e1| + |e2|)/2=', tmp1
      do  i=1,3
         write(LOT,7) (mij(i,j),j=1,3)
      enddo
    7 format(10x,3e12.5)
      return
      end
      

      subroutine crack
     1(ndata,np,index,tol,thresh,zmom,ev,evp,evd,z)
      parameter (LER=0, LIN=5, LOT=6)
      implicit real*8 (a-h,o-z)
      dimension z(3,3),ev(3),evd(3),evp(3)
      dimension index(3),iev(3)
      jj=0
      do 1 i=1,ndata
      if(abs(evd(i)).lt.thresh) jj=jj+1
    1 continue
      write(LOT,*) ' '
      write(LOT,'(a)') 
     1'  --------------------------------------------------------------'
      WRITE(6,*)'EVD:',evd
      WRITE(6,*)'EVP:',evp
      write(LOT,*) ' '
      write(LOT,'(a)') 
     $'  DECOMPOSITION INTO A DOUBLE COUPLE AND A CRACK'
      if(jj.eq.1) then
      write(LOT,*) ' '
      write(LOT,'(a)') 
     $'  DEVIATORIC MOMENT TENSOR IS DUE TO A DOUBLE COUPLE'
      call xyztp2(evd,z,np,ndata,tol,dip,stk,slp)
      else if(jj.eq.2) then
      write(LOT,*) ' '
      write(LOT,'(a)') 
     $'  DEVIATORIC MOMENT TENSOR IS DUE TO A VECTOR DIPOLE'
      else if(jj.eq.0) then
      write(LOT,*) ' '
      write(LOT,'(a)') 
     $'  SEISMIC MOMENTS OF ISOTROPIC PART, DC, CRACK'
      eps=evp(1)/evp(3)
      amj = evp(3) * (1.0 - 2.0 * eps)
      bmj = evp(3) * eps 
      write(LOT,5) zmom- (5./2.)*bmj,amj,bmj*3./2.
c---- MAJOR COUPLE
      ev(index(1)) = 0.0
      ev(index(2)) = -evd(index(3))
      ev(index(3)) = evd(index(3))
      write(LOT,*) ' '
      write(LOT,'(a)') '  EIGENVALUES AND EIGENVECTORS OF THE DC'
      do 3 i=1,ndata
      write(LOT,4) ev(i),(z(j,i),j=1,ndata)
    3 continue
    4 format(e15.7,2x,5f11.7)
    5 format(5e15.7)
      call xyztp2(ev,z,np,ndata,tol,dip,stk,slp)
      iev(1) = 1
      iev(2) = 1
      iev(3) = 1
      iev(index(3)) = 3
      write(LOT,*) ' '
      write(LOT,'(a)') '  EIGENVALUE          AXES OF CRACK'
      do 6 i=1,ndata
      write(LOT,7) iev(i),(z(j,i),j=1,3)
    6 continue
    7 format(i9,2x,5f11.7)
      endif
      return
      end

      subroutine dyad(a1a1,a2a2,a3a3,z)
      implicit real*8 (a-h,o-z)
      dimension z(3,3),a1a1(3,3),a2a2(3,3),a3a3(3,3)
      a1a1(1,1) = z(1,1) * z(1,1)
      a1a1(2,2) = z(2,1) * z(2,1)
      a1a1(3,3) = z(3,1) * z(3,1)
      a1a1(1,2) = z(1,1) * z(2,1)
      a1a1(1,3) = z(1,1) * z(3,1)
      a1a1(2,3) = z(2,1) * z(3,1)
      a1a1(2,1) = a1a1(1,2)
      a1a1(3,1) = a1a1(1,3)
      a1a1(3,2) = a1a1(2,3)
      a2a2(1,1) = z(1,2) * z(1,2)
      a2a2(2,2) = z(2,2) * z(2,2)
      a2a2(3,3) = z(3,2) * z(3,2)
      a2a2(1,2) = z(1,2) * z(2,2)
      a2a2(1,3) = z(1,2) * z(3,2)
      a2a2(2,3) = z(2,2) * z(3,2)
      a2a2(2,1) = a2a2(1,2)
      a2a2(3,1) = a2a2(1,3)
      a2a2(3,2) = a2a2(2,3)
      a3a3(1,1) = z(1,3) * z(1,3)
      a3a3(2,2) = z(2,3) * z(2,3)
      a3a3(3,3) = z(3,3) * z(3,3)
      a3a3(1,2) = z(1,3) * z(2,3)
      a3a3(1,3) = z(1,3) * z(3,3)
      a3a3(2,3) = z(2,3) * z(3,3)
      a3a3(2,1) = a3a3(1,2)
      a3a3(3,1) = a3a3(1,3)
      a3a3(3,2) = a3a3(2,3)
      return
      end

      subroutine xyztp2(ev,xmt,np,nmt,tol,dip,stk,slp)
c---- This program gives the nodal planes of the earthquake
c---- rupture together with the orientation of the pressure
c---- and tension axes, as well as the null vector.
c---- The eigenvectors of the moment tensor are converted
c---- to the null-vector and the tension and pressure axes.
      implicit real*8 (a-h,o-z)
      dimension ev(np),xmt(np,np)
      call trans1(stkt,plnt,stkp,plnp,ev,xmt,np,tol,nmt)
      call tpdss(stkt,plnt,stkp,plnp,dip,stk,slp)
      return
      end
 
      subroutine trans1(stkt,plnt,stkp,plnp,ev,xmt,np,tol,nmt)
      parameter (LER=0, LIN=5, LOT=6)
      implicit real*8 (a-h,o-z)
      dimension ev(np),xmt(np,np)
      wmax=0.
      do 1 i=1,nmt
      if(abs(ev(i)).gt.wmax) wmax=abs(ev(i))
    1 continue
      thresh=tol*wmax
      do 2 i=1,nmt
      if(abs(ev(i)).lt.thresh) then
      a31=xmt(1,i)
      a32=xmt(2,i)
      a33=xmt(3,i)
      endif
      if(ev(i).lt.-thresh) then
      p1=xmt(1,i)
      p2=xmt(2,i)
      p3=xmt(3,i)
      endif
      if(ev(i).gt.thresh) then
      t1=xmt(1,i)
      t2=xmt(2,i)
      t3=xmt(3,i)
      endif
    2 continue
      a11=0.5*(t1+p1)
      a12=0.5*(t2+p2)
      a13=0.5*(t3+p3)
      a21=t1-a11
      a22=t2-a12
      a23=t3-a13
      tnorm=1.0
      if(t3.le.0.0) tnorm=-tnorm
      t1=t1/tnorm
      t2=t2/tnorm
      t3=t3/tnorm
      pnorm=1.0
      if(p3.le.0.0) pnorm=-pnorm
      p1=p1/pnorm
      p2=p2/pnorm
      p3=p3/pnorm
      xnorm=sqrt(2.)
      if(a13.lt.0.0) xnorm=-xnorm
      ynorm=sqrt(2.)
      if(a23.lt.0.0) ynorm=-xnorm
      znorm=1.0
      if(a33.lt.0.0) znorm=-1.0
      a11=a11/xnorm
      a12=a12/xnorm
      a13=a13/xnorm
      a21=a21/ynorm
      a22=a22/ynorm
      a23=a23/ynorm
      a31=a31/znorm
      a32=a32/znorm
      a33=a33/znorm
c------
c     now all vectors point into the lower hemisphere.
c     To get the correct orientations, we require
c     that if the center of the focal sphere is a compression
c     that the X,Y, Z axes form a right handed coordinate system in the
c     lower hemisphere, otherwise it will form a left handed
c     coordinate system
c-----
c     obtain P-wave polarity at the center of the focal
c     sphere
c-----
      xy = a13*a23
c-----
c     determine if right handed or left handed coordinate system
c-----
      z3 = a11*a22 - a12*a21
c-----
c     make right handed coordinate system
c-----
      if(z3.lt.0.0)then
            tmp1=a11
            tmp2=a12
            tmp3=a13
            a11=a21
            a12=a22
            a13=a23
            a21=tmp1
            a22=tmp2
            a23=tmp3
      endif
      x1234 = 1.0
      if(sign(x1234,xy).ne.sign(x1234,rake))then
            tmp1=a11
            tmp2=a12
            tmp3=a13
            a11=a21
            a12=a22
            a13=a23
            a21=tmp1
            a22=tmp2
            a23=tmp3
      endif
      call phase(a11,a12,a13,delx,betx)
      call phase(a21,a22,a23,dely,bety)
      call phase(a31,a32,a33,delz,betz)
      call phase(t1,t2,t3,delt,bett)
      call phase(p1,p2,p3,delp,betp)
      plx = 90.-betx
      ply = 90.-bety
      plz = 90.-betz
      plt=90.-bett
      plp=90.-betp
      stkt=delt
      plnt=90.0-bett
      stkp=delp
      plnp=90.0-betp
      return
      end

      subroutine phase(x,y,z,del,bet)
c-----
c     This expresses a three component vector in terms
c     of two angles in a spherical coordinate system
c-----
c     x     x-coordinate of vector
c     y     y-coordinate of vector
c     z     z-coordinate of vector
c-----
c     del   direction of projection of the vector on a
c           horizontal plane, measured positively, clockwise
c           from north - the trend of the vector.
c     bet   angle that the 3-D vector makes with the
c           downward vertical. This is just 90 - plunge
c-----
      implicit real*8 (a-h,o-z)
      degrad=0.017452329
c     due to numerical inaccuracy we need:
      if(abs(z).le.1.0) then
          bet=acos(z)
      else if(z.ge.1.0) then
          bet = 0.0
c       write(6,*) ' ACOS(Z) = 0; Z = ',z
      else if(z.le.-1.0) then
          bet = 3.14159265
c       write(6,*) ' ACOS(Z) = 0; Z = ',z
      endif
      bet=bet/degrad
C     if(x)21,20,22
C  20 if(y)23,24,25
C  23 del=270.
C     return
C  24 del=0.0
C     return
C  25 del=90.
C     return
C  21 del=atan(y/x)/degrad+180.
C     go to 28
C  22 if(y)26,27,27
C  27 del=atan(y/x)/degrad
C     go to 28
C  26 del=atan(y/x)/degrad+360.
C  28 if(del.gt.360.) del=del-360.
c     if(del.lt.0.0) del=del+360.
      del = atan2(y,x)/degrad
      if(del.lt.0.0)then
          del = del+ 360.
      endif
      if(del.gt.360.0)then
          del = del- 360.
      endif
      return
      end

      subroutine tpdss(dstkt,dplnt,dstkp,dplnp,dip0,stk0,rake0)
      parameter (LER=0, LIN=5, LOT=6)
c-----
c     This takes the orientations of the tension and pressure
c     axes, and determines the two nodal planes
c-----
c
c     convert (P,T) to (dip,rake,stk) and (X,Y,Z).
c-----
      implicit real*8 (a-h,o-z)
      deg=3.141592653/180.0
c-----
c     strike angle is measured from the north clockwise.
c     plunge angle is measured from the horizon downward.
c-----
      if(dstkt.lt.-999.) go to 200
      stkp=dstkp*deg
      plnp=dplnp*deg
      stkt=dstkt*deg
      plnt=dplnt*deg
      t1=cos(plnt)*cos(stkt)
      t2=cos(plnt)*sin(stkt)
      t3=sin(plnt)
      p1=cos(plnp)*cos(stkp)
      p2=cos(plnp)*sin(stkp)
      p3=sin(plnp)
      f1=t1+p1
      f2=t2+p2
      f3=t3+p3
      yy=sqrt(f1*f1+f2*f2+f3*f3)
      f1=f1/yy
      f2=f2/yy
      f3=f3/yy
      xn1=t1-p1
      xn2=t2-p2
      xn3=t3-p3
      yy=sqrt(xn1*xn1+xn2*xn2+xn3*xn3)
      xn1=xn1/yy
      xn2=xn2/yy
      xn3=xn3/yy
      if(xn3.gt.0.0)then
            xn1=-xn1
            xn2=-xn2
            xn3=-xn3
            f1=-f1
            f2=-f2
            f3=-f3
      endif
c     due to numerical inaccuracy we need:
      if(abs(xn3).le.1.0) then
          dip=acos(-xn3)
      else if(xn3.ge.1.0) then
          dip = 3.14159265
c       write(6,*) ' DIP(xn3) = 0; xn3 = ',xn3
      else if(xn3.le.-1.0) then
          dip = 0.0
c       write(6,*) ' DIP(xn3) = 0; xn3 = ',xn3
      endif
      stk=atan2(-xn1,xn2)
      xx=f1*xn2-f2*xn1
      rake=atan2(-f3,xx)
      dip0=dip/deg
      rake0=rake/deg
c  bad
C     if(rake0.lt.-180.)rake0=rake0+180.
      if(rake0.lt.-180.)rake0=rake0+360.
      if(rake0.gt.180.)rake0=rake0-360.
      stk0=stk/deg
      if(stk0.lt.0.0) stk0=360.0+stk0
      if(stk0.lt.0.001) stk0=0.
      if(abs(dip0).lt.0.01) dip0=0.
      if(abs(rake0).lt.0.01) rake0=0.
      write(LOT,14)
   14 format(/' ',' NODAL PLANES ')
      write(LOT,'(a,f10.5)') '   STK= ',stk0
      write(LOT,'(a,f10.5)') '   DIP= ',dip0
      write(LOT,'(a,f10.5)') '  SLIP= ',rake0
      if(f3.gt.0.0)then
            xn1=-xn1
            xn2=-xn2
            xn3=-xn3
            f1=-f1
            f2=-f2
            f3=-f3
      endif
c     due to numerical inaccuracy we need:
      if(abs(f3).le.1.0) then
            dipp=acos(-f3)
      else if(f3.ge.1.0) then
            dipp = 3.14159265
c       write(6,*) ' DIPP(f3) = 0; f3 = ',f3
      else if(f3.le.-1.0) then
            dipp = 0.0
c       write(6,*) ' DIPP(f3) = 0; f3 = ',f3
      endif
      stkk=atan2(-f1,f2)
      xx=f2*xn1-f1*xn2
      rakep=atan2(-xn3,xx)
      dip0=dipp/deg
      rake0=rakep/deg
      stk0=stkk/deg
      if(rake0.lt.-180.)rake0=rake0+360.
      if(rake0.gt.180.)rake0=rake0-360.
      if(stk0.lt.0.0) stk0=360.0+stk0
      if(stk0.lt.0.001) stk0=0.
      if(abs(dip0).lt.0.01) dip0=0.
      if(abs(rake0).lt.0.01) rake0=0.
      write(LOT,*) '  OR'
      write(LOT,'(a,f10.5)') '   STK= ',stk0
      write(LOT,'(a,f10.5)') '   DIP= ',dip0
      write(LOT,'(a,f10.5)') '  SLIP= ',rake0
      call trans(dip0,stk0,rake0,stkt,plnt,stkp,plnp,1)
200   continue
      return
      end

      subroutine trans(dip,stk,rake,stkt,plnt,stkp,plnp,iverby)
      parameter (LER=0, LIN=5, LOT=6)
      implicit real*8 (a-h,o-z)
c-----
c     This takes strike, dip, and rake angles and defines a coordinate
c     transformation in which the rake vector, fault normal vector
c     and the normal vector of the auxiliary plane in terms of the
c     local NS, EW, and vertical cartesian coordinate system.
c-----
c     dip   angle of dip measured from horizontal plane.
c           0 < dip < 90
c     stk   direction of nodal plane strike. When looking
c           along stirke, fault plane dips downward to the
c           right
c           0 < stk < 360
c     rake  direction of movement along nodal plane.
c           If -180 < rake < 180, then the P-wave
c           vertically downward will a first motion polarity
c           equal to the sign of rake
c-----
c
c     REFERENCE:
c
c     Herrmann, R. B. (1975). A student's guide to the use of
c     P and S wave data, Earthquake Notes 46, 29-39
c
c-----
c     The X, Y and Z axes form a right handed coordinate
c     system, such that the compressional quadrant is
c     occurs whenever the xy > 0. The Z axis is the
c     null axis. The Y-axis is normal to the fault
c     plane and the X-axis is normal to the auxiliary plane
c     Note that for the input convention on dip, the
c     Y-axis will always point in the negative z-direction
c
c-----
      degrad=0.017452329
      sins=sin(stk*degrad)
      coss=cos(stk*degrad)
      sind=sin(dip*degrad)
      cosd=cos(dip*degrad)
      sinr=sin(rake*degrad)
      cosr=cos(rake*degrad)
c-----
c     X-axis
c-----
      a11=cosr*coss+sinr*cosd*sins
      a12=cosr*sins-sinr*cosd*coss
      a13=-sinr*sind
c-----
c     Y-axis
c-----
      a21=-sins*sind
      a22=coss*sind
      a23=-cosd
c-----
c     Z-axis
c-----
      a31=coss*sinr-cosd*cosr*sins
      a32=cosd*cosr*coss+sinr*sins
      a33=cosr*sind
      t1=a11+a21
      t2=a12+a22
      t3=a13+a23
      tnorm=sqrt(2.)
      if(t3.le.0.0) tnorm=-tnorm
      t1=t1/tnorm
      t2=t2/tnorm
      t3=t3/tnorm
      p1=a11-a21
      p2=a12-a22
      p3=a13-a23
      pnorm=sqrt(2.)
      if(p3.le.0.0) pnorm=-pnorm
      p1=p1/pnorm
      p2=p2/pnorm
      p3=p3/pnorm
      xnorm=1.0
      if(a13.lt.0.0) xnorm=-1.0
      ynorm=1.0
      if(a23.lt.0.0) ynorm=-1.0
      znorm=1.0
      if(a33.lt.0.0) znorm=-1.0
      a11=a11/xnorm
      a12=a12/xnorm
      a13=a13/xnorm
      a21=a21/ynorm
      a22=a22/ynorm
      a23=a23/ynorm
      a31=a31/znorm
      a32=a32/znorm
      a33=a33/znorm
c------
c     now all vectors point into the lower hemisphere.
c     To get the correct orientations, we require
c     that if the center of the focal sphere is a compression
c     that the X,Y, Z axes form a right handed coordinate system in the
c     lower hemisphere, otherwise it will form a left handed
c     coordinate system
c-----
c     obtain P-wave polarity at the center of the focal
c     sphere
c-----
      xy = a13*a23
c-----
c     determine if right handed or left handed coordinate system
c-----
      z3 = a11*a22 - a12*a21
c-----
c     make right handed coordinate system
c-----
      if(z3.lt.0.0)then
            tmp1=a11
            tmp2=a12
            tmp3=a13
            a11=a21
            a12=a22
            a13=a23
            a21=tmp1
            a22=tmp2
            a23=tmp3
      endif
      x1234 = 1.0
      if(sign(x1234,xy).ne.sign(x1234,rake))then
            tmp1=a11
            tmp2=a12
            tmp3=a13
            a11=a21
            a12=a22
            a13=a23
            a21=tmp1
            a22=tmp2
            a23=tmp3
      endif
      call phase(a11,a12,a13,delx,betx)
      call phase(a21,a22,a23,dely,bety)
      call phase(a31,a32,a33,delz,betz)
      call phase(t1,t2,t3,delt,bett)
      call phase(p1,p2,p3,delp,betp)
      plx = 90.-betx
      ply = 90.-bety
      plz = 90.-betz
      plt=90.-bett
      plp=90.-betp
      stkt=delt
      plnt=90.0-bett
      stkp=delp
      plnp=90.0-betp
      if(iverby.eq.1)then
      if(abs(a11).lt.0.01) a11=0.
      if(abs(a12).lt.0.01) a12=0.
      if(abs(a13).lt.0.01) a13=0.
      if(abs(a21).lt.0.01) a21=0.
      if(abs(a22).lt.0.01) a22=0.
      if(abs(a23).lt.0.01) a23=0.
      if(abs(a31).lt.0.01) a31=0.
      if(abs(a32).lt.0.01) a32=0.
      if(abs(a33).lt.0.01) a33=0.
      if(abs(t1).lt.0.01) t1=0.
      if(abs(t2).lt.0.01) t2=0.
      if(abs(t3).lt.0.01) t3=0.
      if(abs(p1).lt.0.01) p1=0.
      if(abs(p2).lt.0.01) p2=0.
      if(abs(p3).lt.0.01) p3=0.
            write(LOT,*) ' '
      write(LOT,*) '           X-DIR          Y-DIR          Z-DIR'
            write(LOT,'(a,3f15.7)') '  X: ',a11,a12,a13
            write(LOT,'(a,3f15.7)') '  Y: ',a21,a22,a23
            write(LOT,'(a,3f15.7)') '  Z: ',a31,a32,a33
            write(LOT,'(a,3f15.7)') '  T: ',t1,t2,t3
            write(LOT,'(a,3f15.7)') '  P: ',p1,p2,p3
            write(LOT,*) ' '
            write(LOT,*) '       TREND       PLUNGE'
      if(abs(delx).lt.0.01) delx=0.
      if(abs(plx).lt.0.01) plx=0.
      if(abs(dely).lt.0.01) dely=0.
      if(abs(ply).lt.0.01) ply=0.
      if(abs(delz).lt.0.01) delz=0.
      if(abs(plz).lt.0.01) plz=0.
      if(abs(stkt).lt.0.01) stkt=0.
      if(abs(plnt).lt.0.01) plnt=0.
      if(abs(stkp).lt.0.01) stkp=0.
      if(abs(plnp).lt.0.01) plnp=0.
            write(LOT,'(a,3f12.5)') '  X ',delx,plx
            write(LOT,'(a,3f12.5)') '  Y ',dely,ply
            write(LOT,'(a,3f12.5)') '  Z ',delz,plz
            write(LOT,'(a,3f12.5)') '  T ',stkt,plnt
            write(LOT,'(a,3f12.5)') '  P ',stkp,plnp
      endif
      return
      end

      subroutine mteig(ndata,np,ev,xmt)
c-----
c     ndata   I  - dimension of matrix =3
c     np      I  - dimention of matrix =3
c     ev      R  - array of eigenvectors  
c     xmt     R  - 3x3 array of moment tensor
c-----
      implicit real*8 (a-h,o-z)
      integer LER, LIN, LOT
      parameter (LER=0, LIN=5, LOT=6)
      integer ndata, np
      real*8 xmt(3,3),ev(3),ev1(3)
      real*8 z(3,3)
c
c     This program determines the eigenvalues and eigenvectors of
c     a seismic moment tensor. The program uses the Householder 
c     transformation with further QL decomposition (Press et.al. 1987).
c
      do 3 i=1,ndata
      write(LOT,4) (xmt(i,j),j=1,ndata)
    3 continue
    4 format(5e15.7)
      write(LOT,*) ' '
        call tred2(np,np,xmt,ev,ev1,z)
        call tql2(np,np,ev,ev1,z,ierr)
      call out3(ev,xmt,np,ndata)
      return
      end
  
      subroutine out3(ev,z,np,ndata)
      parameter (LER=0, LIN=5, LOT=6)
      implicit real*8 (a-h,o-z)
      dimension z(np,np),ev(np)
      write(LOT,*) ' '
      write(LOT,'(a)') '  EIGENVALUES                EIGENVECTORS'
      do 1 j=1,ndata
      write(LOT,2) ev(j),(z(i,j),i=1,ndata)
    1 continue
    2 format(e15.7,2x,5f11.7)
      return
      end
      
      subroutine gcmdln(ndec,xmt)
      parameter (LER=0, LIN=5, LOT=6)
      integer ndec(6)
      real*8 xmt(3,3)
c-----
c     parse command line arguments and return control
c     parameters
c
c     requires subroutine getarg(i,name) to return
c           the i'th argument in the string name
c
c     and the function numarg() to return the number
c           of arguments excluding program name
c           The first argument is i = 1
c
c-----
      character*20 name
      integer*4 mnmarg
      nmarg = mnmarg()
      if(nmarg.eq.0)goto 9000
      ndec(1) = 0
      ndec(2) = 0
      ndec(3) = 0
      ndec(4) = 0
      ndec(5) = 0
      ndec(6) = 0
        do i=1,3
            do j=1,3
                xmt(j,i) = 0.0
            enddo
        enddo
      i = 0
 1000 continue
      i = i + 1
      if(i.gt.nmarg)goto 9000
            call getarg(i,name)
            if(name(1:3).eq.'-xx' .or. name(1:3).eq.'-XX')then
                  i=i+1
                  call getarg(i,name)
                  call chtofp(name,xmt(1,1))
            else if(name(1:3).eq.'-yy' .or. name(1:3).eq.'-YY')then
                  i=i+1
                  call getarg(i,name)
                  call chtofp(name,xmt(2,2))
            else if(name(1:3).eq.'-zz' .or. name(1:3).eq.'-ZZ')then
                  i=i+1
                  call getarg(i,name)
                  call chtofp(name,xmt(3,3))
            else if(name(1:3).eq.'-xy' .or. name(1:3).eq.'-XY')then
                  i=i+1
                  call getarg(i,name)
                  call chtofp(name,xmt(1,2))
            else if(name(1:3).eq.'-yx' .or. name(1:3).eq.'-YX')then
                  i=i+1
                  call getarg(i,name)
                  call chtofp(name,xmt(1,2))
            else if(name(1:3).eq.'-xz' .or. name(1:3).eq.'-XZ')then
                  i=i+1
                  call getarg(i,name)
                  call chtofp(name,xmt(1,3))
            else if(name(1:3).eq.'-zx' .or. name(1:3).eq.'-ZX')then
                  i=i+1
                  call getarg(i,name)
                  call chtofp(name,xmt(1,3))
            else if(name(1:3).eq.'-yz' .or. name(1:3).eq.'-YZ')then
                  i=i+1
                  call getarg(i,name)
                  call chtofp(name,xmt(2,3))
            else if(name(1:3).eq.'-zy' .or. name(1:3).eq.'-ZY')then
                  i=i+1
                  call getarg(i,name)
                  call chtofp(name,xmt(2,3))
            else if(name(1:3).eq.'-mc')then
                  ndec(1) = 1
            else if(name(1:3).eq.'-dc')then
                  ndec(2) = 1
            else if(name(1:3).eq.'-vd')then
                  ndec(3) = 1
            else if(name(1:5).eq.'-clvd')then
                  ndec(4) = 1
            else if(name(1:3).eq.'-cc')then
                  ndec(5) = 1
            else if(name(1:6).eq.'-crack')then
                  ndec(6) = 1
            else if(name(1:2).eq.'-a')then
                  ndec(1) = 1
                  ndec(2) = 1
                  ndec(3) = 1
                  ndec(4) = 1
                  ndec(5) = 1
                  ndec(6) = 1
                  goto 9000
            else if(name(1:2).eq.'-h')then
                  call usage()
            else if(name(1:2).eq.'-?')then
                  call usage()
            endif
            go to 1000
 9000 continue
c-----
c     make the moment tensor symmetric
c-----
      xmt(3,2) = xmt(2,3)
      xmt(3,1) = xmt(1,3)
      xmt(2,1) = xmt(1,2)
c-----
c     if the output form is not specified, give everything
c-----
      if(ndec(1).eq.0 .and. ndec(2).eq.0 .and. ndec(3).eq.0.0
     1  .and. ndec(4).eq.0.0 .and. ndec(5).eq.0.0 
     2  .and. ndec(6).eq.0)then
                  ndec(1) = 1
                  ndec(2) = 1
                  ndec(3) = 1
                  ndec(4) = 1
                  ndec(5) = 1
                  ndec(6) = 1
      endif
      return
      end

      subroutine usage()
        parameter (LER=0, LIN=5, LOT=6)
        write(LER,*)'Usage: mtinfo '
        write(LER,*)' -XX Mxx -YY Myy -ZZ Mzz -XY -Mxy -XZ Mxz -YZ Myz'
        write(LER,*)' [-dc|-mc|-cc|-vd|-clvd|-a|-crack]'
        write(LER,*) ' '
        write(LER,*)
     1  ' -XX Mxx -YY Myy -ZZ Mzz  Moment tensor elements in units of'
        write(LER,*)
     1  ' -XY Mxy -XZ Mxz -YZ Myz    dyne-cm'
        write(LER,*)
     1  '               Note -xy Mxy works as well as -yx Mxy'
      write(LER,*) 
     1  ' -dc               3 double couples'
      write(LER,*) 
     1  ' -mc               major couple  '
      write(LER,*)
     1  ' -cc               double couble and  CLVD'
      write(LER,*) 
     1  ' -vd               3 vector dipoles   '
      write(LER,*)
     1  ' -clvd             3 CLVDs  '
      write(LER,*)
     1  ' -crack            opening crack'
      write(LER,*)
     1  ' -a                all of the above'
      write(LER,*)
     1  '     Note: if none are set, -a is the default'
        write(LER,*)
     1  ' -?                   This online help'
        write(LER,*)
     1  ' -h                   This online help'
        stop
        end

        subroutine chtofp(str,fout)
c------
c     routine to convert string to floating point
c     The E format is accepted as well as free form
c     input
c
c     This assumes that the string is 20 characters long
c-----
        character*(*) str
        real*8 fout
        integer*4 lgstr
        logical hase
        l = lgstr(str)
c------
c     If the string str contains an E or e, then
c     we parse the result using an E format.
c------
        hase = .false.
        do 100 i=1,l
            if(str(i:i) .eq.'e' .or. str(i:i).eq.'E')then
                hase = .true.
            endif
  100   continue
c-----
c     read floating point number
c-----
        if(hase)then
            read(str,'(bn,e20.0)')fout
        else
            read(str,'(bn,f20.0)') fout
        endif
        return
        end

        subroutine gtev(xmt,mom,np,ev,ev1,ilg,ism,z)
c-----
c       Input:
c           mom(3,3)   R - moment tensor
c           np         I - size of mom matrix, e.g., here np=3
c       Output:
c           xmt(3,3)   R - array of eigenvectors
c           ev         R - array of eigenvalues
c           ev1
c           ilg        I - index of largest eigenvalue
c           ism        I - index of smallest eigenvalue
c-----
      implicit real*8 (a-h,o-z)
        real*8 xmt(np,np),ev(np),ev1(np)
        real*8 z(3,3)
        real*8 mom(3,3)
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
        write(LOT,*)' EIGENVALUES AND EIGENVECTOR OF M(i,j)'
    2   format(' ',i5,1x,e11.4,'  (',e11.4,',',e11.4,',',e11.4,')')
        do 120 j=1,3
            write(LOT,2)j,ev(j),(z(i,j),i=1,3)
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
       
c-----
c      transliterated from Chuck Ammon s  pradiation_lpr.c
c      NOTE the outptu of Chucks ios HRV not AR
c      I input AR
c
c      BE CAREFUL HERE THE IDEA OF P N and T axes only works of some of the
c      eigenvalues are positive and the others are negative
c-----

       subroutine PrintPradiation(ofile,Mij,xmom,plnt,stkt,plnp,stkp,
     1       plnn,stkn,ev,bet,gam)
c-----
c      Print the P-wave radiation pattern
c      ofile   C - name of output file
c      Mij     R - 3x3 moment tensor
c      plnt, stkt  R - plunge and strike of T axis
c      plnp, stkp  R - plunge and strike of P axis
c      plnn, stkn  R - plunge and strike of N axis
c      ev          R - array of 3 eigenvalues
c      evec        R - 3x3 array of eigenvectors, ordered by
c                      column
c-----
      implicit none
      character ofile*(*)
      real*8 Mij(3,3)
      real*8 plnt,stkt,plnp,stkp,plnn,stkn
      real*8 xmom
      real*8 ev(3), evec(3,3)

      integer FID
      integer nrows, ncols, kMAXchar
      real*8 mw

      real*8 HMij(3,3)
      integer NR, NC
      parameter (NR=132, NC=132)
      character thePattern(NR)*132
      real x0, y0, rad
      real l1, l2 , l3
      real bet, gam

      integer nlines
c-----
c     initialize
c----
      nlines= NR
      kMAXchar = NC
      nrows=29

   
      if(ofile.eq.'stdout')then
           FID = 6
      else
           FID = 1
           open(FID,file=ofile,access='sequential',status='unknown',
     1           form='formatted')
           rewind FID
      endif
      call ar_to__hrv_mij(Mij,HMij)
C RBH c-----
C RBH c     annotate
C RBH c-----
C RBH       write(FID,'(a,2x,1pe9.2,2x,a)')' Moment (dyne-cm)',xmom,' dyne-cm'
C RBH c-----
C RBH c     IASPEI definition of Mw
C RBH c-----
C RBH       mw = 2.*(dlog10(xmom) - 16.1)/3.0
C RBH       write(FID,'(a,2x,f6.2)')' Magnitude (Mw)',mw
C RBH       write(FID,*)' '
C RBH       write(FID,'(a)')' Principal Axes:'
C RBH       write(FID,'(a)')'   Axis    Value   Plunge  Azimuth'
C RBH       write(FID,'(a,1x,1pe9.2,2(0pf8.0))')'    T ',ev(3),plnt,stkt
C RBH       write(FID,'(a,1x,1pe9.2,2(0pf8.0))')'    N ',ev(2),plnn,stkn
C RBH       write(FID,'(a,1x,1pe9.2,2(0pf8.0))')'    P ',ev(1),plnp,stkp
C RBH 
C RBH       write(FID,'(a)')' Moment Tensor: (dyne-cm) Aki-Richards'
C RBH       write(FID,'(a)')'    Component   Value'
C RBH       write(FID,'(a,1pe9.2)')'       Mxx   ', Mij(1,1)
C RBH       write(FID,'(a,1pe9.2)')'       Mxy   ', Mij(1,2)
C RBH       write(FID,'(a,1pe9.2)')'       Mxz   ', Mij(1,3)
C RBH       write(FID,'(a,1pe9.2)')'       Myy   ', Mij(2,2)
C RBH       write(FID,'(a,1pe9.2)')'       Myz   ', Mij(2,3)
C RBH       write(FID,'(a,1pe9.2)')'       Mzz   ', Mij(3,3)
c-----

c----
c     compute the radiation pattern from the moment tensor
c----
      call GetlprPradiation(thePattern,Mij,nrows,ncols,
     1     stkt,plnt,stkp,plnp,x0,y0,rad)
c----
c     now plot the Tape Lune but offset the x-axis
c-----
      x0 = x0 + 4.0*rad
      call lunegrid(x0,y0,rad,'equal',thePattern,nrows,ncols)
c-----
c     plot the value for this moment tensor
c-----
      l1 = sngl(ev(1))
      l2 = sngl(ev(2))
      l3 = sngl(ev(3))
      call Prlune(x0,y0,rad,l1,l2,l3,thePattern,nrows,ncols,bet,gam)
c-----
c     annotate
c-----
      write(FID,'(a,2x,1pe9.2,2x,a)')' Moment (dyne-cm)',xmom,' dyne-cm'
c-----
c     IASPEI definition of Mw
c-----
      mw = 2.*(dlog10(xmom) - 16.1)/3.0
      write(FID,'(a,2x,f6.2)')' Magnitude (Mw)',mw
      write(FID,*)' '
      write(FID,'(a)')' Principal Axes:'
      write(FID,'(a)')'   Axis    Value   Plunge  Azimuth'
      write(FID,'(a,1x,1pe9.2,2(0pf8.0))')'    T ',ev(3),plnt,stkt
      write(FID,'(a,1x,1pe9.2,2(0pf8.0))')'    N ',ev(2),plnn,stkn
      write(FID,'(a,1x,1pe9.2,2(0pf8.0))')'    P ',ev(1),plnp,stkp

      write(FID,'(a38x,14x,a)')' Moment Tensor: (dyne-cm) Aki-Richards',
     1    'Lune parameters'
      write(FID,'(a)')'    Component   Value'
c-----
c     note use of p scale factor to shift the exponential
c     to that 1.0e=22 appears instead of 0.100e+23
c     hwoever this has to be turned off for otuer output
c-----
      write(FID,2)           '       Mxx   ', Mij(1,1), bet
    2 format(a,1pe9.2,32x,'beta: ',0pf7.2)
    3 format(a,1pe9.2,32x,'gamma:',0pf7.2)
      write(FID,3)           '       Mxy   ', Mij(1,2), gam
      write(FID,'(a,1pe9.2)')'       Mxz   ', Mij(1,3)
      write(FID,'(a,1pe9.2)')'       Myy   ', Mij(2,2)
      write(FID,'(a,1pe9.2)')'       Myz   ', Mij(2,3)
      write(FID,'(a,1pe9.2)')'       Mzz   ', Mij(3,3)
c-----
c     output the GCMT convention
c----
       write(FID,'(a,a)')' Global CMT Convention Moment Tensor:',
     1   ' (dyne-cm)'
       write(FID,'(a)')'         R         T         F'
       write(FID,'(a,3(1pe9.2,1x))')'  R ',
     1       HMij(1,1),HMij(1,2),HMij(1,3)
       write(FID,'(a,3(1pe9.2,1x))')'  T ',
     1       HMij(2,1),HMij(2,2),HMij(2,3)
       write(FID,'(a,3(1pe9.2,1x))')'  F ',
     1       HMij(3,1),HMij(3,2),HMij(3,3)
c-----
c----
c     output the radiation pattern
c----
      call PrintRadiation(FID,thePattern,nrows,ncols)

c----
c     clean up
c----
      if(FID.ne.6)then
          close(FID)
      endif
      return
      end
          
      subroutine PrintRadiation(FID,thePattern,nrows,ncols)
      implicit none
      integer NR, NC
      parameter (NR=132, NC=132)

      integer FID
      character thePattern(NR)*132
      integer nrows, ncols, i, ls
      integer lgstr
           do i=1,nrows
                ls = lgstr(thePattern(i) )
                write(FID,'(a)')thePattern(i)(1:ls)
           enddo
      return
      return
      end

      
      subroutine GetlprPradiation(thePattern,Mij,nrows,ncols,
     1     tstkt,plnt,tstkp,plnp,x0,y0,r0)
      implicit none
      integer NR, NC
      parameter (NR=132, NC=132)

      
      character thePattern(NR)*132
      real*8 Mij(3,3)
      integer nrows,ncols
      real*8 tstkt, plnt, tstkp, plnp

      real x_to_y_ratio, sqrt2
      real x0, y0, r0
      real xx, yy, r
      integer ip, jp, it, jt
      real deg_to_rad
      real rad_to_deg
      real stkt,stkp
      real misp, mist
      real phi, theta
      real pw, svw, shw
      integer i,j

c-----
c     initialize
c-----
      deg_to_rad=0.017453292
      rad_to_deg=57.29577951
      sqrt2 = sqrt(2.0)
      x_to_y_ratio = 1.85
      misp = 1.0e+38
      mist = 1.0e+38
      if(tstkp .gt. 180.0)then
          stkp = tstkp - 360.
      else
          stkp = tstkp
      endif
      if(tstkt .gt. 180.0)then
          stkt = tstkt - 360.
      else
          stkt = tstkt
      endif
c-----
c     initialize the character image
c-----
      do i=1,NR
         thePattern(i) = ' '
      enddo
c-----
c     set up the scale for the radiation pattern image
c     The aspect ratio is meant to assume that the results
c     look good when prointed
c     and make sure the number of columns is odd
c-----
      ncols = int( float(nrows)*x_to_y_ratio)
      ncols = 2*(ncols/2) + 1
C       /* set up the origin of the focal sphere */
      x0 = ((real (ncols) * 0.5 + 1) / x_to_y_ratio)
      y0 =  (real (nrows) * 0.5)
C       /* set the radius of the focal sphere */
      r0 = sqrt(x0*x0 + y0*y0) / 1.85
c-----
c     file in the negative polarity amplitudes
c-----
      do i=1,nrows
           yy = y0 -  float(i)
           do j=1,ncols
                xx = (float (j) / x_to_y_ratio) - x0;
                r = sqrt(xx*xx + yy*yy) / r0;
                if(r .lt. 1.0)then
c-----
c                  phi is the azimuth, theta is the takeoff angle
c-----
                   phi = atan2(xx,yy) * rad_to_deg

c-----
c                  for a stereographic projection 
c                  theta = 2*atan(r) * rad_to_deg; 
c                  for an equal area projection 
c---- -
                   theta = 2*asin(r/sqrt2) * rad_to_deg
                   if(ABS(phi-stkt)+ABS(90.-theta - plnt).lt.mist)then
                        mist = ABS(phi - stkt)+ABS(90-theta - plnt)
                        it = i
                        jt = j
                   endif
                   if(ABS(phi-stkp)+ABS(90.-theta - plnp).lt.misp)then
                        misp = ABS(phi - stkp)+ABS(90-theta - plnp)
                        ip = i
                        jp = j
                   endif

c-----
c                  get the body wave radiation amplitudes */
c-----
                   call GetAmplitudes(Mij,phi,theta,pw,svw,shw);
c-----
C                  assign a character based on the P-wave polarity 
c-----
                   if(pw .gt. 0)then
                                      thePattern(i)(j:j) = '#';
                   else
                                      thePattern(i)(j:j) = '-';
                   endif
                endif
            enddo
       enddo
             do i=ip-1,ip+1
                do j=jp-1,jp+1
                   thePattern(i)(j:j) = ' '
                enddo
             enddo
             thePattern(ip)(jp:jp) = 'P'
             do i=it-1,it+1
                do j=jt-1,jt+1
                   thePattern(i)(j:j) = ' '
                enddo
             enddo
             thePattern(it)(jt:jt) = 'T'
      return
      end
    
      subroutine GetAmplitudes(Mij,phi,theta,pw,svw,shw);
c-----
c     get amplitdue from moment tensor
c-----
      real*8 Mij(3,3)
      real  phi, theta,pw,svw,shw

      real deg_to_rad
      real sth, cth, cphi, sphi, p_dir(3),sv_dir(3),sh_dir(3)
      integer i,j
      integer N,E,D
      parameter (N=1,E=2,D=3)

      deg_to_rad=0.017453292
c        /* first compute some direction cosines */

        sth = sin(theta*deg_to_rad)
        cth = cos(theta*deg_to_rad)
        cphi = cos(phi*deg_to_rad)
        sphi = sin(phi*deg_to_rad)

c       /* now construct the ray vectors */
        p_dir(N) = sth*cphi
        p_dir(E) = sth*sphi
        p_dir(D) = cth

        sv_dir(N) = cth*cphi
        sv_dir(E) = cth*sphi
        sv_dir(D) = -sth

        sh_dir(N) = -sphi
        sh_dir(E) = cphi
        sh_dir(D) = 0

c       /* now compute the radiated amplitudes */

        pw = 0
        svw = 0
        shw = 0

        do i=1,3
           do j=1,3
              pw =  pw + Mij(i,j)*p_dir(i)* p_dir(j)
             svw = svw + Mij(i,j)*p_dir(i)*sv_dir(j)
             shw = shw + Mij(i,j)*p_dir(i)*sh_dir(j)
           enddo
        enddo
      return
      end


      subroutine ar_to__hrv_mij(Mij,HMij)
c-----
c     convert from Aki Richards convention, e.g.,
C     x = north, y = east, z = down
c     to Harvard
c     R = -z, theta = -S, phi - E
c-----
      real*8 Mij(3,3), HMij(3,3)
      integer X, Y, Z, R, T, F
      parameter(X=1,Y=2,Z=3,R=1,T=2,F=3)

      HMij(R,R)=   Mij(Z,Z)
      HMij(R,T)=   Mij(Z,X)
      HMij(R,F)= - Mij(Z,Y)
      HMij(T,T)=   Mij(X,X)
      HMij(T,F)= - Mij(X,Y)
      HMij(F,F)=   Mij(Y,Y)
      HMij(T,R)=  HMij(R,T)
      HMij(F,R)=  HMij(R,F)
      HMij(F,T)=  HMij(T,F)
      

      return
      end

      subroutine gettrpl(trn,plu,z1,z2,z3)
      real*8 trn, plu,z1,z2,z3
      real*8 v1,v2,v3

      real*8 delp, betp
c-----
c     make sure vector points to lower hemispher, e.g., with z3 positive
c-----
      if(z3.lt.0.0)then
          v1 = -z1
          v2 = -z2
          v3 = -z3
      else
          v1 =  z1
          v2 =  z2
          v3 =  z3
      endif
      call phase(v1,v2,v3,delp,betp)
      plu = 90. - betp
      trn = delp
      

      return
      end
c-----
c     plot the lune as a pager plot
c-----
      subroutine Prlune(x0,y0,rad,ol1,ol2,ol3,thePattern,nrows,ncols,
     1     bet,gam)
      implicit none
      integer FID
      real ol1, ol2 , ol3

      integer NR, NC
      parameter (NR=132, NC=132)
      character thePattern(NR)*132
      integer nrows, ncols
      real bet, gam, mnorm
      real degrad
      real l1, l2 , l3
      real x0, y0, rad,x,y
      integer i,j
      real snorm

      DEGRAD=0.01745329251994
c-----
c     plot the lune grid and set up output array
c     this also defines the center of the plot space and radius
c-----
c     order eigenvalues from largest to smallest
c-----
      call order(l1,l2,l3,ol1,ol2,ol3)
      

c-----c
c     plot the particular value on the lune
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
       call luneproj(bet,gam,x0,y0,rad,x,y,'equal')
            j = nint(x)
            i = nint(y)
            i = nrows + 1 - i
            thePattern(i)(j:j) = '#'
      return
      end

      subroutine lunegrid(x0,y0,rad,project,thePattern,nrows,ncols)
c-----
c     plot the line with [ - pi/6] to [pi/6] in longitude
c      x0, y0   - R  coordinates of center of project - returned
c      rad      - R  radius of project - returned
c      project    - C*(*) 'stereo' or 'equal'
c      the Pattern the printer plot pattern
c      nrows, ncols the image dimensions
c-----
      implicit none
      integer NR, NC
      parameter (NR=132, NC=132)

      real x0, y0,rad
      character project*(*)
      character thePattern(NR)*(*)
      integer nrows, ncols



      integer NPTS
      parameter (NPTS=181)
      real x(NPTS), y(NPTS)
      integer igam, ibet
      real bet, gam
      real xx, yy
      real xd, yd, zd
      real sqrt2, x_to_y_ratio
      integer i,j
      integer il,ih,jl,jh
c-----
c     initialize 
c-----
      sqrt2 = sqrt(2.0)
      x_to_y_ratio = 1.80

c-----
c     initialize the character image
c-----
C     do i=1,NR
C        thePattern(i) = ' '
C     enddo
c-----
c     set up the scale for the radiation pattern image
c     The aspect ratio is meant to assume that the results
c     look good when prointed
c     and make sure the number of columns is odd
c-----
C     ncols = int( float(nrows)*x_to_y_ratio)
C     ncols = 2*(ncols/2) + 1
C       /* set up the origin of the focal sphere */
C     x0 = ((real (ncols) * 0.5 + 1) / x_to_y_ratio)
C     y0 =  (real (nrows) * 0.5)
C       /* set the radius of the focal sphere */
C     rad = sqrt(x0*x0 + y0*y0) / 1.85

      do igam = -30,30,15
         gam = igam
         do ibet = 1,NPTS,3
            bet = ibet-1
            call luneproj(bet,gam,x0,y0,rad,xx,yy,project)
            j = nint(xx)
            i = nint(yy)
            i = nrows + 1 - i
            if(igam.eq.-30 .or. igam.eq. 30)then
                 thePattern(i)(j:j) = ':'
            else
                 thePattern(i)(j:j) = '.'
            endif
          enddo
c-----
c         plot the curve
c-----
      enddo
c-----
c     plot the lines of constant beta
c-----
      do ibet=10,170,20
         bet = ibet
         call luneproj(bet,30.0,x0,y0,rad,xx,yy,project)
         i =  nint(yy)
         jh = nint(xx)
         call luneproj(bet,-30.0,x0,y0,rad,xx,yy,project)
         i =  nint(yy)
         jl = nint(xx)
         i = nrows + 1 - i
         do j = jl+1,jh-1
            if(ibet.eq.90)then
                 thePattern(i)(j:j) = '='
            else
                 thePattern(i)(j:j) = '-'
            endif
         enddo
      enddo
      return
      end
          
      subroutine PrintLune(FID,thePattern,nrows,ncols)
      implicit none
      integer NR, NC
      parameter (NR=132, NC=132)
      integer FID
      character thePattern(NR)*(*)
      integer nrows, ncols, i, ls
      integer lgstr
           do i=1,nrows
                ls = lgstr(thePattern(i) )
                write(FID,'(a)')thePattern(i)(1:ls)
           enddo
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
               

      subroutine luneproj(degcolat,deglon,x0,y0,rad,x,y,project)
c-----
c      Input values
c      degcolat - R  co-latitude, e.g, 0 = top (N) and 180 = bottom(S)
c                    = beta
c      deglon   - R  longitude 
c                    = gamma
c      x0, y0   - R  coordinates of center of project
c      rad      - R  radius of project
c      project     - C*(*) - 'stereo' for equal antlke 
c                          'equal' for equal area
c      Return values
c      x,y      - R  output x,y coordinates for the plot
c-----
       implicit none
       real degcolat, deglon, x0, y0, rad, x, y
       character project*(*)

       real radcolat, radlon
       real degrad, sqrt2
       real k

       real deglat
       DEGRAD=0.01745329251994
       sqrt2 = 1.414213562373095
     
       radcolat = degcolat*degrad
       radlon   = deglon*degrad



c-----
c      here x - y are plot coordinates
c      degcolat =   0 is y = y0 + rad
c                 180 is y = y0 -rad
c                        y = y0 + rad * (1 - degcolat/180))
c      deglon   
c                        x = x0 + rad*(sin(deglon)/(sin(30)) *sin(degcolat)
c-----
         y = y0 + rad * (1.0 - 2*degcolat/180.)
c        x = x0 + rad*sin(radlon)*sin(radcolat)/(3.*sin(DEGRAD*30))
         x = x0 + rad*sqrt2*sin(radlon/2.)*sin(radcolat)
c-----
c     Ammon Mathematica code 
c-----
       deglat = 90 - degcolat
       if(project.eq.'stereo')then
           k = (2.)/(1 + cos(degrad*deglat)*cos(degrad*deglon))
           k = k /2.
           x = x0 + 2*rad*k*cos(degrad*deglat)*sin(degrad*deglon)
           y = y0 + rad*k*sin(degrad*deglat)
       else if(project.eq.'equal')then
           k = (2.)/(1 + cos(degrad*deglat)*cos(degrad*deglon))
           k = k /2.
           k = sqrt(k)
           x = x0 + 2*rad*k*cos(degrad*deglat)*sin(degrad*deglon)
           y = y0 + rad*k*sin(degrad*deglat)
       endif
      
       return
       end
