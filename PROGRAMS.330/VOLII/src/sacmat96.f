c----------------------------------------------------------------------c
c                                                                      c
c        COMPUTER PROGRAMS IN SEISMOLOGY                               c
c        VOLUME II                                                     c
c                                                                      c
c        PROGRAM: SACMAT96                                             c
c                                                                      c
c        COPYRIGHT 1999                                                c
c        D. R. Russell, R. B. Herrmann                                 c
c        Department of Earth and Atmospheric Sciences                  c
c        Saint Louis University                                        c
c        221 North Grand Boulevard                                     c
c        St. Louis, Missouri 63103                                     c
c        U. S. A.                                                      c
c                                                                      c
c----------------------------------------------------------------------c
        program sacmat96
c-----
c       13 OCT 00
c       put in period bounds for output is isolated spectrum
c       20 MAR 2004 - output format change for write(LOT)
c       09 NOV 2005 - caught a big error if the input 
c           file is of mixed C U G
c           in subroutine vlstrt
c       29 MAY 2009 - fixed roundoff problem in the if (p.ge.pmin etc
c       23 JUL 2009 - ensured consistent call to FFT routine with
c               complex arguments. We now call zfour
c       11 AUG 2011 - changed the size of fdisp and fsac 
c                 to be character fsac*256, fdisp*256 to be compatible with sacmft96
c                 also changed namev, names, namer to
c                 character*256 namev,names,namer  instead of 40
c                 and names in subroutine gcmdln to be
c                 character*256 names instead of 80
c                 Meijian An, Institute of Geomechanics, 
c                      Chinese Academy of Geological Sciences, Beijing
c       04 JAN 2021 - change all do XXX yyy ... XXX a=b 
c                  to   do yy ... a=b .. enddo
c           to be compatible with gfortran 9.3
c       7 APR 2021 - the reference for the cubic spline fit used in spline() 
c       and bfit() is
c                 Lawson, C. L. and R. J. Hanson (1995). 
c                    Solving Least Squares Problems, SIAM
c                    on pages 222-225.
c-----
c
c       program match calculates corrected phase velocity and
c         group velocity dispersion curves, given an initial estimate
c         of the group velocity and the spectrum of an event.
c
c-----
        implicit none
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)
        integer NP, NP2
        parameter(NP=131072,NP2=65540)
        common/tdum/z0,x,y,yf,yg,c,zsav,df,n,n2,t0,dist,ns
        complex z0(NP)
        real x(NP2),y(NP2),yf(NP2),yg(NP2),c(NP2)
        real df, t0, dist
        integer n, n2, ns
        complex zsav(NP)
        common /pfplot/ paf(NP),waf(NP),npf,nb,pdt,pdf,nfct,ymmax
        real paf, waf, pdt, pdf, ymmax
        integer npf, nb, nfct
        common /tvar/ zc(NP),dt,ntv,nti1,nti2,ntj1,ntj2
        complex zc
        real dt
        integer ntv,nti1,nti2,ntj1,ntj2

        character fsac*256, fdisp*256

        common/auto/doauto,iter,itermx,lshift,phasecor
        logical doauto
        integer iter, itermx
        logical lshift
        real phasecor

        common/per/permin,permax
        real permin, permax

        integer ierr, jret

        iter = 0
c-----
c       machine dependent initialization
c-----
        call mchdep()
c-----
        call gcmdln(fsac,fdisp,doauto,itermx,lshift,phasecor)
        write(LOT,*)'fsac                   :',fsac
        write(LOT,*)'fdisp                  :',fdisp
        write(LOT,*)'itermx                 :',itermx
        write(LOT,*)'doauto                 :',doauto
        write(LOT,*)'output compressed trace:',lshift
        write(LOT,*)'phasecor               :',phasecor
c-----
c       initialize plot output
c-----
        if(.not.doauto)then
            call ginitf('INTER','sacmat96')
        endif
        call getspc(fsac,ierr)
        if(ierr.lt.0)go to 9999
c-----
c       call subroutine to readin initial group velocity estimates
c       and then to compute first estimate of phase velocity curves
c-----
C       write(LOT,*)'Calling vlstrt fdisp:',fdisp
        call vlstrt(fdisp,permin,permax)
        call pmatch()
c-----
c       edit group velocity curve with a linear interpolator
c-----
        if(.not. doauto)then
            call gedit()
        endif
c-----
c       since no longer needed convert original spectra and filtered
c       spectra to time series
c-----
        call zfour(zsav,n,+1,dt,df)
        call zfour(  z0,n,+1,dt,df)
        jret = 0
 1000   continue
            call mtv()
            call timser()
            call savspc(jret,permin,permax)
        if(jret.ne.0)goto 1000
c-----
c       terminate plotter output
c-----
 9999   continue
        if(.not. doauto)then
            call pend()
        endif
        stop
        end
 
        subroutine band(z,n,i1,i2,nb,nt)
        implicit none
        integer NP, NP2
        parameter(NP=131072,NP2=65540)
        complex z(NP)
        integer n, i1, i2, nb, nt
        common/iband/ib1,ib2,ib3,ib4
        integer ib1,ib2,ib3,ib4
        integer m, kl, i
        real fparz
        real fac
c-----
c       apply a parzen window to the complex function z
c
c       z - complex array to be windowed
c       n - number of points in the z array
c       i1,i2 - passband of windowing function in units
c               of array index
c             the original array is unchanged between i1 and i2
c       nb - taper width for parzen window about i1, i2
c          the array is zeroed for indices less than i1 - nb + 1
c                    and for indices greater than    i2 + nb - 1
c       nt = 0 assume spectra input, and only window
c              the positive array frequencies, e.g.,
c              the first half of the array. This is
c              in fact not a true parzen window, but
c              we use parzen half-windows to taper the
c              band edges
c          = 1 true parzen window about the center of the array
c              used here for the pseudo-autocorrelation functions
c-----
c-----
c       define center point, corners, and check for valid 
c       choices of corners
c-----
        if(nt.eq.0)then
            m = n / 2 + 1
        else
            m = n
        endif
        if(nb.lt.2)then
            kl = 2
        else
            kl = nb
        endif
        ib1 = i1 - kl + 1
        if(ib1.lt.1)ib1=1
        ib2 = i1
        ib3 = i2
        ib4 = i2 + kl -1
        if(ib4.gt.m)ib4=m
c-----
        do 100 i = 1 , m
            if(i.lt.ib1)then
                z(i) = cmplx(0.0,0.0)
            elseif(i.ge.ib1 .and. i.lt.ib2)then
                fac = fparz(ib2-i,ib2-ib1)
                z(i) = fac * z(i)
            elseif(i.ge.ib2 .and. i.lt.ib3)then
                z(i) = z(i)
            elseif(i.ge.ib3 .and. i.lt.ib4)then
                fac = fparz(i-ib3,ib4-ib3)
                z(i) = fac * z(i)
            else
                z(i) = cmplx(0.0,0.0)
            endif
  100   continue
        return
        end

        function fparz(i,iw)
        implicit none
        integer i, iw
        real fparz
c-----
c       parzen windowing function
c
c         iw = window halfwidth
c         i  = index within window
c              i = 0 corresponds to the passband = 1
c-----
        real rat, fi, fiw

            fi = i
            fiw = iw
            rat = abs(fi/fiw)
            if(2*i .lt. iw)then
                fparz = 1.0 -6.*rat*rat*(1.-rat)
            else
                fparz = 2.*(1.0-rat)**3
            endif
        return
        end

c
c       subroutine to unwrap phase from complex spectrum zc().
c
        subroutine unwrap(y,dw,i2,j1,j2,m)
        implicit none
        integer NP, NP2
        parameter(NP=131072,NP2=65540)
        real y(NP2), dw
        integer i2, j1, j2, m
        common /tvar/ zc,dt,ntv,nti1,nti2,ntj1,ntj2
        complex zc(NP)
        real dt
        integer ntv,nti1,nti2,ntj1,ntj2
        common /pfplot/ paf(NP),waf(NP),npf,nb,pdt,pdf,nfct,ymmax
        real paf, waf, ymmax, pdt, pdf, xt, yt
        integer npf,nb,nfct
        
        integer i, j, jj
        real cz0, cz1, cz2, u1, u2, u3, v1, v2, v3, x

        y(1)=0.0
        y(2)=atan2(aimag(zc(2)),real(zc(2)))
        y(3)=atan2(aimag(zc(3)),real(zc(3)))
        do 30 i=4,i2
            j=i+1
            if(i.eq.i2)j=i
            cz0=cabs(zc(j))
            cz1=cabs(zc(j-1))
            cz2=cabs(zc(j-2))
            u1=real(zc(j))/cz0
            u2=real(zc(j-1))/cz1
            u3=real(zc(j-2))/cz2
            v1=aimag(zc(j))/cz0
            v2=aimag(zc(j-1))/cz1
            v3=aimag(zc(j-2))/cz2
            x=u2*(v1-v3)-v2*(u1-u3)
      if(i.ne.i2) go to 25
        x=u1*(3.*v1-4.*v2+v3)-v1*(3.*u1-4.*u2+u3)
  25  y(i)=x/(2.*dw)
  30    continue
        jj=0
        do  i=j1,j2
           jj=jj+1
           paf(jj)=y(i)
        enddo
        npf=j2-j1+1
c
c ******** at this point, the cursor plotter program can be *********
c            called to plot group delay errors.  an editor to
c            remove spikes should be incorporated (a linear
c            interpolator between the spike endpoints is fine).
c            the group delay file to plot is in array paf(NPf).
c            the dimension of paf is npf.
c
c********************************************************************
c
        jj=0
        do  i=j1,j2
           jj=jj+1
           y(i)=paf(jj)
        enddo
        yt=y(3)
        xt=0.
        do  i=4,i2
           yt=(y(i)+xt)*dw/2.+yt
           xt=y(i)
           y(i)=yt
        enddo
        zc(1)=cmplx(0.0,0.0)
        return
        end

        subroutine bfit(x,y,c,m,nfd)
c-----
c       use spline coefficients to define
c       fit to data
c
c       x   - array of abscissae
c       y   - array to be evaluated at x values
c       c   - array of spline coefficients
c       m   - number of points in x,y arrays
c       nfd - 0 evaluate spline
c           - 1 estimate derivative at point
c           - 2 estimate integral at point
c-----
        implicit none
        integer NP2
        parameter(NP2=65540)
        real x(NP2),y(NP2),c(NP2)
        integer m, nfd
        common/bsp/b(NP2),q(4),h
        real b, q, h
        real t
        real p1, p2, cp1, cp2, dp1, dp2
        real sum, suma, fa, fb, u, yfit
        integer ifl, jt, i, l, ic
c-----
c       FORTRAN FUNCTION DEFINITIONS
c-----
        p1(t)=.25*t**2*t
        p2(t)=-(1.-t)**2*(1.+t)*.75+1.
        cp1(t)=t**4/16.
        cp2(t)=t*(1./4.+t*(3./8.+t*(1./4.-3./16.*t)))
        dp1(t)=.75*t**2
        dp2(t)=.75+t*(1.5-2.25*t)
        sum=0.
        suma=0.
        ifl=0
        jt=1
        do 110 i=1,m
  80  if(x(i).le.b(jt+1)) go to 90
        if(nfd.lt.2) go to 85
        ifl=0
        fa=-(c(jt)+11.*c(jt+1))*h/16.
        fb=(11.*c(jt+2)+c(jt+3))*h/16.
        sum=sum+fb-fa
  85  jt=jt+1
        go to 80
  90  u=(x(i)-b(jt))/h
        if(nfd.ne.1 .and. nfd.ne.2)then
            q(1)=p1(1.-u)
            q(2)=p2(1.-u)
            q(3)=p2(u)
            q(4)=p1(u)
        elseif(nfd .eq. 1)then
            q(1)=-dp1(1.-u)/h
            q(2)=-dp2(1.-u)/h
            q(3)=dp2(u)/h
            q(4)=dp1(u)/h
        elseif(nfd .eq. 2)then
            if(ifl.ne.1) then
                ifl=1
                fa=-(c(jt)+11.*c(jt+1))*h/16.
                suma=sum-fa
            endif
            q(1)=-cp1(1.-u)*h
            q(2)=-cp2(1.-u)*h
            q(3)=cp2(u)*h
            q(4)=cp1(u)*h
        endif
        yfit=0.0
        do 100 l=1,4
            ic=jt-1+l
            yfit=yfit+c(ic)*q(l)
  100   continue
        y(i)=yfit+suma
  110   continue
        return
        end

        subroutine spline(x,y,c,m,nbp)
        implicit none
        integer NP, NP2
        parameter(NP=131072,NP2=65540)
        real x(NP2),y(NP2),c(NP2)
        integer m, nbp

        common/bsp/b(NP2),q(4),h
        real b, q, h
        real g(NP2,5)
        real p1, p2, zero, u, t, rnorm
        integer mdg, nband, nc, ir, ip, i, jt, mt, ig
c-----
c       FORTRAN FUNCTION DEFINITION
c-----
        p1(t)=.25*t**2*t
        p2(t)=-(1.-t)**2*(1.+t)*.75+1.
C       write(0,*)'spline:m,nbp',m,nbp
        zero=0.
        mdg=NP2
        nband=4
        if(nbp.gt.m-2) nbp=m-4
        if(nbp.lt.2) nbp=2
        nc=nbp+2
        b(1)=x(1)
        b(nbp)=x(m)
        h=(b(nbp)-b(1))/float(nbp-1)
        if(nbp.le.2) go to 30
        do  i=3,nbp
            b(i-1)=b(i-2)+h
        enddo
  30    continue
        ir=1
        ip=1
        i=1
        jt=1
  40  mt=0
  50    continue
        if(x(i).gt.b(jt+1)) go to 60
        u=(x(i)-b(jt))/h
        ig=ir+mt
        g(ig,1)=p1(1.-u)
        g(ig,2)=p2(1.-u)
        g(ig,3)=p2(u)
        g(ig,4)=p1(u)
        g(ig,5)=y(i)
        mt=mt+1
        if(i.eq.m) go to 60
        i=i+1
        go to 50
  60    continue
        call bndacc(g,mdg,nband,ip,ir,mt,jt)
        if(i.eq.m) go to 70
        jt=jt+1
        go to 40
  70    continue
        call bndsol(1,g,mdg,nband,ip,ir,c,nc,rnorm)
        return
        end

        subroutine bndacc(g,mdg,nb,ip,ir,mt,jt)
c-----
c       c.l. lawson and r.j. hanson, jet propulsion laboratory
c         sequential algorithm for banded least squares problem.
c       accumulation phase.   for solution phase use bndsol.
c
c       the calling program must set ir=1 and ip=1 before the
c         first call to bndacc for a new case.
c
c       the second subscript of g() must be dimensioned at least
c         nb+1 in the calling program.
c-----
        parameter(NP=131072,NP2=65540)
        dimension g(NP2,5)
        zero=0.
        nbp1=nb+1
        if(mt.le.0) return
        if(jt.eq.ip) go to 70
        if(jt.le.ir) go to 30
        do  i=1,mt
           ig1=jt+mt-i
           ig2=ir+mt-i
           do  j=1,nbp1
               g(ig1,j)=g(ig2,j)
           enddo
        enddo
        ie=jt-ir
        do  i=1,ie
           ig=ir+i-1
           do  j=1,nbp1
              g(ig,j)=zero
           enddo
        enddo
        ir=jt
  30  mu=min0(nb-1,ir-ip-1)
        if(mu.eq.0) go to 60
        do  l=1,mu
           k=min0(l,jt-ip)
           lp1=l+1
           ig=ip+l
           do  i=lp1,nb
              jg=i-k
              g(ig,jg)=g(ig,i)
           enddo
           do  i=1,k
              jg=nbp1-i
              g(ig,jg)=zero
           enddo
        enddo
  60  ip=jt
  70  mh=ir+mt-ip
        kh=min0(nbp1,mh)
        do  i=1,kh
        call h12(1,i,max0(i+1,ir-ip+1),mh,g(ip,i),
     1           1,rho,g(ip,i+1),1,mdg,nbp1-i)
        enddo
        ir=ip+kh
        if(kh.lt.nbp1) go to 100
        do  i=1,nb
           g(ir-1,i)=zero
        enddo
 100    continue
        return
        end

        subroutine bndsol(mode,g,mdg,nb,ip,ir,x,n,rnorm)
c-----
c       c.l. lawson and r.j. hanson, jet propulsion laboratory
c       sequential solution of a banded least squares problem.
c       solution phase.  for the accumulation phase us bndacc.
c
c       the second subscript of g() must be dimensioned at
c       least nb+1 in the calling program
c-----
        parameter (LIN=5,LOT=6,LER=0)
        parameter(NP=131072,NP2=65540)
        dimension g(NP2,5),x(NP2)
        zero=0.
        rnorm=zero
        go to (10,90,50), mode
  10  do  j=1,n
         x(j)=g(j,nb+1)
      enddo
        rsq=zero
        np1=n+1
        irm1=ir-1
        if(NP1.gt.irm1) go to 40
        do  j=np1,irm1
           rsq=rsq+g(j,nb+1)**2
        enddo
        rnorm=sqrt(rsq)
  40    continue
  50  do  ii=1,n
        i=n+1-ii
        s=zero
        l=max0(0,i-ip)
        if(i.eq.n) go to 70
        ie=min0(n+1-i,nb)
        do  j=2,ie
           jg=j+l
           ix=i-1+j
           s=s+g(i,jg)*x(ix)
        enddo
C  70  if(g(i,l+1)) 80,130,80
   70   continue
      if(g(i,l+1) .ne. 0.0) then
            go to 80
        else
            go to 130
        endif
  80  x(i)=(x(i)-s)/g(i,l+1)
        enddo
        return
  90  do  j=1,n
        s=zero
        if(j.eq.1) go to 110
        i1=max0(1,j-nb+1)
        i2=j-1
        do  i=i1,i2
           l=j-i+1+max0(0,i-ip)
           s=s+x(i)*g(i,l)
        enddo
 110  l=max0(0,j-ip)
C       if(g(j,l+1)) 120,130,120
        if(g(j,l+1) .ne. 0.0)then
            go to 120
        else
            go to 130
        endif
 120  x(j)=(x(j)-s)/g(j,l+1)
        enddo
        return
 130    continue
c-----
c       cleanup
c-----
c-----
c       terminate program
c-----
        call pend()
        write(LER,*)'zero diagonal term in bndsol'
        stop
        end
        subroutine h12 (mode,lpivot,l1,m,u,iue,up,c,ice,icv,ncv)
c-----
c       construction and/or application of a single
c       householder transformation..     q = i+u*(u**t)/b
c
c       mode    = 1 or 2   to select algorithm h1 or h2.
c       lpivot is the index of the pivot element.
c       l1,m   if l1 .le. m   the transformation will be constructed to
c              zero elements indexed from l1 through m.   if l1 gt. m
c              the subroutine does an identity transformation.
c       u(),iuw,up    on entry to h1 u() contains the pivot vector.
c                   iue is the storage increment between elements.
c                                       on exit from h1 u() and up
c                   contain quantities defining the vector u of the
c                   householder transformation.   on entry to h2 u()
c                   and up should contain quantities previously computed
c                   by h1.  these will not be modified by h2.
c       c()    on entry to h1 or h2 c() contains a matrix which will be
c              regarded as a set of vectors to which the householder
c              transformation is to be applied.  on exit c() c
c               contains the set of transformed vectors.
c       ice    storage increment between elements of vectors in c().
c       icv    storage increment between vectors in c().
c       ncv    number of vectors in c() to be transformed. if ncv .le. 0
c              no operations will be done on c().
c-----
        parameter(NP=131072,NP2=65540)
        dimension u(iue,m),c(5*NP2)
        double precision sm,b
        one=1.
c
        if (0.ge.lpivot.or.lpivot.ge.l1.or.l1.gt.m) return
        cl=abs(u(1,lpivot))
        if (mode.eq.2) go to 60
c              ****** construct the transformation. ******
            do 10 j=l1,m
           cl=amax1(abs(u(1,j)),cl)
  10      continue
C       if (cl) 130,130,20
        if( cl .le. 0.0)then
            go to 130
        else
            go to 20
        endif
  20  clinv=one/cl
        sm=(dble(u(1,lpivot))*clinv)**2
            do 30 j=l1,m
           sm=sm+(dble(u(1,j))*clinv)**2
  30      continue
c               convert dble. prec. sm to sngl. prec. sm1
        sm1=sm
        cl=cl*sqrt(sm1)
C       if (u(1,lpivot)) 50,50,40
        if (u(1,lpivot) .le. 0.0) then
            go to 50
        else
            go to 40
        endif
  40  cl=-cl
  50  up=u(1,lpivot)-cl
        u(1,lpivot)=cl
        go to 70
c       ****** apply the transformation 1+u*(u**t)/b to c. ******
c
C  60  if (cl) 130,130,70
   60   continue
        if(cl .le. 0.0 )then
            go to 130
        else
            go to 70
        endif
  70  if (ncv.le.0) return
        b=dble(up)*u(1,lpivot)
c            b  must be nonpositive here.  if b=0., return.
c
C       if (b) 80,130,130
        if ( b .lt. 0.0 )then
            go to 80
        else
            go to 130
        endif
  80  b=one/b
        i2=1-icv+ice*(lpivot-1)
        incr=ice*(l1-lpivot)
            do 120 j=1,ncv
            i2=i2+icv
            i3=i2+incr
            i4=i3
            sm=c(i2)*dble(up)
                do  i=l1,m
                sm=sm+c(i3)*dble(u(1,i))
                i3=i3+ice
                enddo
C          if (sm) 100,120,100
            if(sm .ne. 0.0)then
                go to 100
            else
                go to 120
            endif
 100      sm=sm*b
            c(i2)=c(i2)+sm*dble(up)
                do 110 i=l1,m
                c(i4)=c(i4)+sm*dble(u(1,i))
                i4=i4+ice
 110        continue
 120      continue
 130      return
            end

        subroutine inter(sx,sy,x,y,w1,dw,num)
c-----
c       sx  array of data angular frequencies
c       sy  array of data group velocities
c
c       x   output array of interpolated natural frequencies
c       y   output array of interpolated group velocities
c-----
        parameter(NP=131072,NP2=65540)
        dimension sx(NP2),sy(NP2),x(NP2),y(NP2),xx(4),yy(4)
        x(1)=w1
        i=1
  113   continue
            if (sx(2).lt.x(i)) go to 114
            s1=(x(i)-sx(1))/(sx(2)-sx(1))
            y(i)=s1*(sy(2)-sy(1))+sy(1)
            i=i+1
            x(i)=x(i-1)+dw
        go to 113
  114   continue
        n1=num-2
        do 120 j=2,n1
            j1=j-1
            do 118 k=1,4
                xx(k)=sx(j1+k-1)
                yy(k)=sy(j1+k-1)
  118       continue
  119       if (sx(j+1).lt.x(i)) go to 120
            call pccp(xx,yy,x(i),y(i))
            i=i+1
            x(i)=x(i-1)+dw
            go to 119
  120   continue
  121   continue
            if (sx(num).lt.x(i)) go to 122
            s1=(x(i)-sx(num-1))/(sx(num)-sx(num-1))
            y(i)=s1*(sy(num)-sy(num-1))+sy(num-1)
            i=i+1
            x(i)=x(i-1)+dw
        go to 121
  122   num=i-1
        return
        end

        subroutine pccp(x,y,xd,yd)
        dimension x(4),y(4)
c  **********************************************
c  subroutine pccp performs piecewise continuous cubic
c  polynomial interpolation following wiggins(1976)
c  and akima(1970). the method requires two x and y
c  points on each side of position xd where the value of yd is
c  determined. weighted averages of the slopes are determined at the
c  knots(positions x(2) and x(3)), and the slopes and the values x(2),
c  y(2),x(3),y(3) are used to determine yd which corresponds to
c  x position xd which falls between x(2) and x(3).
c  references
c  wiggins,r.a.,bull. seism. soc. am.,v.66,p2077-2081,1976.
c  akima,h.,j.assoc.comp.mach.,v.17,p 589-602,1970.
c
c  note the x-values must be distinct!!
c  *********************************************************
        eps=0.001
c  determine slopes at x(2) and x(3)
        sx2=(y(2)-y(1))/(x(2)-x(1))
        sx3=(y(3)-y(2))/(x(3)-x(2))
        sx4=(y(4)-y(3))/(x(4)-x(3))
c  weight slopes.
        w1=1./amax1(abs(sx2),eps)
        w2=1./amax1(abs(sx3),eps)
        w3=1./amax1(abs(sx4),eps)
        s2=(w2*sx3+w3*sx4)/(w2+w3)
        s1=(w1*sx2+w2*sx3)/(w1+w2)
c  evaluate polynomial at x=xd
c  with slopes given at x(2) and x(3)
        p0=y(2)
        p1=s1
        p2=(3.*(y(3)-y(2))/(x(3)-x(2))-2.*s1-s2)/(x(3)-x(2))
        p3=(s1+s2-2.*(y(3)-y(2))/(x(3)-x(2)))/((x(3)-x(2))**2)
        xmxd=xd-x(2)
        yd=p0+p1*xmxd+p2*xmxd**2+p3*xmxd**3
        return
        end

        subroutine vlstrt(fdisp,permin,permax)
        parameter (LIN=5,LOT=6,LER=0)
        parameter(NP=131072,NP2=65540)
        implicit complex (z)
        dimension z0(NP),x(NP2),y(NP2),yf(NP2),yg(NP2),c(NP2)
        dimension zsav(NP)
        common/tdum/z0,x,y,yf,yg,c,zsav,df,n,n2,t0,dist,ns
        character*256 namev,namer,names
        common/fname/namev,names,namer
        character*2 info
        common /dpplot/ xgo(NP2),ygo(NP2),xgn(NP2),ygn(NP2),
     $                xpo(NP2),ypo(NP2),xpn(NP2),ypn(NP2),
     $                ddf,moo,mpn,ml,mh,in1,in2,nmod,m
        common /pfplot/ paf(NP),waf(NP),npf,nb,pdt,pdf,nfct,ymmax
        common /tvar/ zc(NP),dt,ntv,nti1,nti2,ntj1,ntj2
        common/vwind/i1,i2,xmx,ii1,ii2
        common/iband/ib1,ib2,ib3,ib4
        common/srfdat/iunit,ifrpr,ilr,imd
        character txt1*80

        common/auto/doauto,iter,itermx,lshift,phasecor
        logical doauto
        integer iter, itermx
        logical lshift
        real phasecor

        common/break/nbpts
        integer nbpts

        character fdisp*(*)
c-----
c       specifications for dispersion file
c-----
        integer lun, ierr, ndsp,  NOBS
        parameter(NOBS=4098)
        integer*4 jlorr(NOBS), jobs(NOBS), jobsyn(NOBS),jmode(NOBS)  
        real*4  fper(NOBS)   , fobs(NOBS)   , fobserr(NOBS)

        tupi = 2.* 3.141592654
        dw = tupi*df
        ddf = df
        ls=lgstr(fdisp)
        if(.not. doauto)then
            call gwrtxt(1.0,6.5,'Group velocity file :',0)
            lnm = 40
            call gwrtxt(4.0,6.5,fdisp(1:ls),0)
        endif
        namev=fdisp(1:ls)//'v'
   18   continue
c------
c       input initial group velocity measurements using surf(IV)
c       format. Only group velocity data will be used. However,
c       the file can only contain one wave type and one mode
c-----
c
c       x()  - periods 
c       y()  - group velocities
c       m    - total number of group velocities
c-----
        ifrpr = 0
        nmod = 1
        m = 0
        ilrold = 0
        imdold = 0
        lun = 2
C       write(6,*)'Trying to read fdisp:',fdisp
        call rddisp(fdisp,lun,ierr,ndsp,NOBS,
     1          jlorr,jobs,jobsyn,jmode,
     2          fper,fobs,fobserr)
c-----
c       convert period to frequency
c       set up group velocity array
c-----
        m = 0
        do 25 i=1,ndsp
            if(jobs(i).eq.2)then
                m = m + 1
                x(m) = fper(i)
                y(m) = fobs(i)
                imd = jmode(i)
                ilr = jlorr(i)
                if(m.gt.1)then
                if(imdold.ne.imd )then
                    call gwrtxt(1.0,6.1,
     1              'MULTIMODE DATA NOT PERMITTED',0)
                    ierr = -10
                elseif(ilrold.ne.ilr)then
                    call gwrtxt(1.0,5.7,
     1              'MIXED LOVE-RAYL NOT PERMITTED',0)
                    ierr = -10
                else
                    ierr = 0
                endif
                else
                    ilrold = ilr
                    imdold = imd
                endif
                if(ierr .lt. 0)stop
            endif
   25   continue
c-----
c       now sort the period-group velocity in order of
c       increasing period
c-----
        call sort(x,y,m)
        xmx=x(m)
        permin = x(1)
        permax = x(m)
c-----
c       eliminate any identical periods
c-----
        call uniq(x,y,m)
c-----
c reexamine the velocity file to see if it contains phase velocity data
c in the range of permin to permax, and if so output these values
c so for information for use with the next query
c-----
        tmin = permax
        tmax = permin
        iper = 0
        if(.not. doauto)then
            write(txt1,10)permin
   10       format('Group Velocity Minimum Period = ',f11.4)
   11       format('Group Velocity Maximum Period = ',f11.4)
            call gwrtxt(1.0,6.0,txt1,0)
            write(txt1,11)permax
            call gwrtxt(1.0,5.6,txt1,0)
C           if(iper.gt.0)then
C              call gwrtxt(1.0,5.2,'FOR INFORMATION. Dispersion File',0)
C              call gwrtxt(1.0,4.8,'has phase velocity information',0)
C              write(txt1,12)cmin,tmin
C   12      format('Phase velocity=',f11.4,' at ',f11.4,' seconds')
C              call gwrtxt(1.0,4.4,txt1,0)
C              write(txt1,12)cmax,tmax
C              call gwrtxt(1.0,4.0,txt1,0)
            call gwrtxt(1.0,3.0,'Enter reference period:',0)
            call grdtxt(txt1,80)
            read(txt1,'(f20.0)')f0
            call gwrtxt(1.0,2.6,'Enter reference phase velocity:',0)
            call grdtxt(txt1,80)
            read(txt1,'(f20.0)')c0
            f0=1./f0
        else
            f0 = 1.0
            c0 = 1.0E+5
        endif
        nb=int(1.5*xmx/dt) + 1
        moo=m
c-----
c       xgo(moo) - old period array
c       ygo(moo) - old group velocity array
c-----
        do 50 i=1,m
            xgo(i)=x(i)
            ygo(i)=y(i)
   50   continue
 1050 continue
        if(.not. doauto)then
            call gwrtxt(1.0,2.0,
     1      'Enter # breakpoints for group velocity spline:',0)
            call grdtxt(txt1,80)
            read(txt1,'(i5)',err=1050)nbp
        else
            nbp = nbpts
        endif
c
c       fit wiggin's spline interpolator to original group velocity
c       data in order to calculate group velocity in terms of frequency
c       sampled increments
c
        nfd=0
        WRITE(6,*)'M:',M
        do 170 i=1,m
            j=m-i+1
            yg(j)=tupi/x(i)
            yf(j)=y(i)
  170   continue
        ntw=0
        i1=int((yg(1)/dw)*1.00001)+1
        w1=float(i1-1)*dw
        call inter(yg,yf,x,y,w1,dw,m)
        i2=i1+m-1
c
c       output of wiggin's spline in x,y.  find least squares cubic
c       spline coefficients to x,y (spline) and store in c.  fit spline
c       to y (bfit) and store in yf.
c
c       nbp - number of breakpoints
c       nfd - set to zero for function fit
c
        if(nbp.lt.0) nbp=m-4
C       write(6,*)'call spline'
        call spline(x,y,c,m,nbp)
C       write(6,*)'call bfit'
        call bfit(x,yf,c,m,nfd)
c
c       enter band edges for spectrum
c
c       ns - number of array points for band edge (default = 10)
c       nt - set to zero for positive spectrum band filter
c
        nt=0
        n2=n/2+1
        ns=10
        do 1000 i=1,n
            z0(i) = zsav(i)
 1000   continue
C       write(6,*)'call band:n,i1,i2,ns,nt',n,i1,i2,ns,nt
        call band(z0,n,i1,i2,ns,nt)
c-----
c      extend group velocity data by a linear extension of the
c      data endpoints (this is done to minimize edge effects).
c
c      ii1 - lower passband limit (in terms of frequency array location)
c      ii2 - upper passband limit
c      i1 - extended lower passband limit
c      i2 - extended upper passband limit
c      m  - number of group velocity points in extended passband
c      ml - lower passband limit (in terms of period array location)
c      mh - upper passband limit
c-----
        j1=i1-ns+2
        j2=i2+ns-2
C       write(6,*)'j2,i2,ns,n2,n:',j2,i2,ns,n2,n
        if(j1.lt.2) j1=2
        if(j2.gt.n2) j2=n2-1
        ml=i1-j1
        mh=j2-i2
C       write(6,*)'mh,j2,i2,n2:',mh,j2,i2,n2
        sl=yf(2)-yf(1)
        sh=yf(m)-yf(m-1)
        yl=yf(1)
        yh=yf(m)
        do 200 i=1,m
              j=m-i+1
              x(j+ml)=x(j)
              yf(j+ml)=yf(j)
 200  continue
        do 210 i=1,ml
              dx=float(i-1-ml)
              x(i)=x(ml+1)+dx*dw
              yf(i)=sl*dx+yl
 210  continue
        do 220 i=1,mh
              dx=float(i)
              x(m+ml+i)=x(m+ml)+dx*dw
              yf(m+ml+i)=sh*dx+yh
 220  continue
        mm=m
C       write(6,*)'m,ml,mh:',m,ml,mh
        m=m+ml+mh
        ml=ml+1
        mh=m-mh
        mll=ml
        ml=m-mh+1
        mh=m-mll+1
        ii1=i1
        i1=j1
        ii2=i2
        i2=j2
        mpn=m
        ntv=n
        nti1=i1
        nti2=i2
        ntj1=ii1
        ntj2=ii2
c
c       end of band extension
c
        do 230 i=1,m
            y(i)=dist/yf(i)-t0
  230   continue
c
c       refit group velocity spline with extended data
c
C       write(6,*)'call spline:m,nbp',m,nbp
        call spline(x,y,c,m,nbp)
C       write(6,*)'call bfit:m,nfd',m,nfd
        call bfit(x,yf,c,m,nfd)
c
c       xgn - group velocity periods
c       ygn - new group velocities
c
        do 240 i=1,m
        j=m-i+1
        xgn(i)=tupi/x(j)
        ygn(i)=dist/(yf(j)+t0)
  240 continue
        do 260 i=1,m
            y(i1+i-1)=yf(i)
  260   continue
c
c       integrate spline fit of group delay for phase estimate
c
c       nfd=2 for spline integral
c
        nfd=2
        call bfit(x,yf,c,m,nfd)
c
c       calculate integration constant (ctg) based on reference
c       period and phase velocity (c0,f0)
c       note that this does not affect the time shift but rather just
c       defines a constant phase adjustment
        if(.not. doauto)then
            w=tupi*f0
            ic=int(w/dw)+1
            if(ic.gt.i2) ic=i2
            if(ic.lt.i1) ic=i1
            ic=ic-i1+1
            ctg=x(ic)*(dist/c0-t0)-yf(ic)
            do 270 i=1,m
                yf(i)=yf(i)+ctg
  270       continue
        endif
c-----
cThe purpose of CTG is to use the correct phase velocity to give
czerophase at the source, the actual shifting is done by the group delay
c-----
        
c
c       xpn - phase velocity periods
c       ypn - new phase velocity measurements
c
        do 280 i=1,m
        j=m-i+1
        xpn(i)=tupi/x(j)
        ypn(i)=dist*x(j)/(yf(j)+x(j)*t0)
  280 continue
        nmod=3
c
c *********at this point call automatic scale x-log, y-linear*********
c            plotting routine. plot:
c                xgo,ygo - dashed curve    (dimension moo)
c                xgn,ygn - solid curve     (dimension mpn)
c                xpn,ypn - solid curve     (dimension mpn)
c            on this particular plot, observe the following convention:
c            plot xgo,ygo from 1 to moo (observed data).  plot xgn,ygn
c            and xpn,ypn between the array limits ml and mh to insure
c            that band edges are not plotted.
c
c ********************************************************************
        if(.not.doauto)then
            call gpplot(info)
c
c       if spline fit is undesirable, repeat fitting procedure
c
c
            if(info(1:1).eq.'y')call gframe(1)
            if(info(1:1).eq.'y') go to 18
        endif
        j=1
        do 300 i=i1,i2
            x(i)=yf(j)
            yg(j)=float(i-1)*dw
            j=j+1
  300   continue
c-----
c       fill in negative frequency component of filtered spectrum
c-----
        do 325 i=2,n2
            j = n + 2 - i
            z0(j) = conjg(z0(i))
  325   continue
        return
        end

        subroutine sort(x,y,n)
c-----
c       This provides a sort of items
c       x(i) is the array to be sorted
c       which has associated with it an y(i) array
c       no is the number of points to be sorted
c       After the sort the x(i),y(i) pairs will be ordered
c       from the least to the largest
c-----
c-----
c      Reference: http://en.wikipedia.org/wiki/Bubble_sort
c-----
       integer n
       real x(n), y(n)

       do i = n, 1, -1
           do j = 1 , i -1
               if(x(j) .gt. x(j+1))then
                   tmp = x(j)
                   x(j) = x(j+1)
                   x(j+1) = tmp
                   tmp = y(j)
                   y(j) = y(j+1)
                   y(j+1) = tmp
                endif
           enddo
       enddo
       return
       end

        subroutine mtv()
        parameter (LIN=5,LOT=6,LER=0)
        parameter (NP=131072,NP2=65540)
        complex  zc(NP),shift,zp(NP)
        dimension rlngth(NP2),ftr(NP)
        double precision sntp,csw,csdw,snw,sndw,csr,csdr,snr,sndr
        common /tvar/ zc,dt,ntv,nti1,nti2,ntj1,ntj2
        real*4 cur(NP2)
        common /pfplot/ paf(NP),waf(NP),npf,nb,pdt,pdf,nfct,ymmax
        common/vwind/id1,id2,xmx,ii1,ii2
        character*2 info
        character txt1*80

        common/auto/doauto,iter,itermx,lshift,phasecor
        logical doauto
        integer iter, itermx
        logical lshift
        real phasecor
c-----
        nbx=nb
        nb = int(2.0*xmx/dt) + 1
        itvar = 0
c
c       read in control values
c
c       iflm = 0 for matched filter
c              1 for time variable filter (tvf)
c       alph   = number of period multiples for tvf (default=2.5)
c       ipwr   = power of cosine window (default = 1)
c       err    = maximum relative amplitude error for
c                spectrum (default=.05)
c
 1001 continue
        if(.not. doauto)then
            call gwrtxt(1.0,5.0,'enter -1 for FFT matched filter',0)
            call gwrtxt(1.0,4.6,'enter  0 for time variable filter,',0)
            call gwrtxt(1.0,4.2,'enter  1 for time variable filter ',0)
            call gwrtxt(1.0,3.8,'         curvature correction:',0)

            call gwrtxt(1.0,3.4,'RESPONSE :',0)
            call grdtxt(txt1,80)
            read(txt1,'(i5)',err=1001)iflm
            if(iflm.lt.-1 .or. iflm.gt.1)goto 1001
            if(iflm.eq.1)then
 1002           continue
                call gwrtxt(1.0,3.0,
     1          'Input Maximum Relative Error (0.1):',0)
                call grdtxt(txt1,80)
                read(txt1,'(f20.0)',err=1002)err
            elseif(iflm.eq.0)then
                      err = 0.00
            endif
        else
            iflm = -1
        endif
        alph=2.5
c-----
c       reset meaning of iflm to force time variable filter rather
c       than Discrete Fourier Transform if iflm = 0
c-----
        if(iflm.eq.0)iflm = 1
        ipwr=1
c
c  initialization
c
        do 5 i=1,NP2
            rlngth(i)=0
    5   continue
        tupi=2*3.141592654
        n=ntv
        i1=nti1
        i2=nti2
        j1=ntj1
        j2=ntj2
        npf=ntv
        nd2=n/2+1
        df=pdf
        dt=pdt
        if(iflm.ne.1) go to 15
        xpar=float(nbx-1)*dt
        cor=(6.)/xpar**2
        do 67 i=1,n
            zp(i) = cmplx(waf(i),0.0)
   67   continue
        nt =1
        k1 = nd2
        k2 = nd2
        call band(zp,n,k1,k2,nbx,nt)

c
c  calculate second derivative of spectrum
c
c-----
c       save pseudo-autocorrelation in array paf
c-----
        do 8 i=1,n
            paf(i)=real(zp(i))
    8   continue
        ifr = -1
        call zfour(zp,n,ifr,dt,df)
c-----
c calculate spectra of second moment of pseudo-autocorrelation function
c-----
        do 10 i=1,n
            t=dt*float(i-nd2)
            tmp=-paf(i)*t**2
            zc(i)=cmplx(tmp,0.0)
   10   continue
        ifr = -1
        call zfour(zc,n,ifr,dt,df)
        cmp=-1.0
        do 11 i=j1,j2
        ctmp=cabs(zp(i)-cor*zc(i))
        cmp=amax1(ctmp,cmp)
  11  continue
        do 12 i=1,nd2
            cur(i)=0.0
            if(cmp.gt.0.0)then 
                cur(i)=cabs(zc(i))/cmp
            endif
  12  continue
        cmax=-1.
        do 13 i=i1,i2
            if(cmax.ge.cur(i)) go to 13
            cmax=cur(i)
            tmax=1./(df*float(i-1))
  13  continue
        if(itvar .eq. 1)goto 14
        itvar = 1
        nb=int(alph*float(n)/float(j1-1))+1
  14  continue
        nbt=0
        if(err.gt.0.0)then
            nbt=int(tupi*sqrt(cmax/(32.*err))/dt)+1
            if(nbt.gt.nb)nb=nbt
        endif
  15  continue
c-----
c       we now know something about time variable filter windows
c       now get spectra of pseudo-autocorrelation trace
c-----
        do 16 i=1,n
            zc(i)=cmplx(waf(i),0.0)
  16    continue
        ifr = -1
        call zfour(zc,n,ifr,dt,df)
c-----
c   save only frequencies between i1 and i2 which are the extended bands
c   from the phase match filter
c-----
        do 17 i=1,nd2
            if(i.lt.i1 .or. i.gt.i2)zc(i) = cmplx(0.0,0.0)
            if(i.gt.1)zc(n+2-i) = conjg(zc(i))
   17   continue
c-----
c  output the unfiltered time trace for display
c-----
        do 77 i=1,n
            paf(i)=waf(i)
   77   continue
c
c*********at this point, the time plotter program can be ***********
c           called.  the original pseudo-autocorrelation
c           function (paf) is plotted.  the cursor can be used
c           to window the paf, and the value (nb) should be
c           passed back.  this will be the one sided width of
c           the window in array units.
c
c           the paf function is in array paf(NPf).
c
c******************************************************************
c
        iflag=0
        if(.not. doauto)then
            call swplot(iflag,info)
        endif
c-----
c  calculate the width of the window for each frequency
c  calculate discrete inverse fft
c-----
        wmax = (nb-1)*dt
c-----
c       cosine window pseudo-autocorrelation function
c       in time domain with window of fixed length
        if(iflm.eq.-1)then
            call zfour(zc,n,+1,dt,df)
            ll = nd2 - nb
            lr = nd2 + nb
            if(ll.le.0)ll=1
            if(lr.gt.n)lr=n
            do 151 i=1,n
                if(i.lt.ll .or. i.gt.lr)then
                    fac = 0.0
                else
                    rw=float(i-nd2)/float(nb)
                    rw = rw*1.57079633
                    fac=cos(rw)**ipwr
                endif
            ftr(i)=real(zc(i))*fac
  151          continue
        else
            do 105 i=i1,i2
                if(i.lt.j1)then
                    fh = df*(j1-1)
                else
                    fh = df*( i-1)
                endif
                aln=alph/fh
                if(iflm.eq.1) then 
                    if(err.gt.0.0) then
                        bln=tupi*sqrt(cur(i)/(32.*err))
                        rlngth(i)=amax1(aln,bln)
                    else
                        rlngth(i)=aln
                    endif
                    if(rlngth(i).gt.wmax)rlngth(i)=wmax
                else
                    rlngth(i)=wmax
                endif
  105       continue
c-----
c           fourier synthesis of windowed components in the passband.
c-----
            if(.not. doauto)then
                call gwrtxt(1.0,7.0,'In fourier synthesis',0)
            endif
c-----
c           perform DFT
c           For each frequency compute the sinousoidal time series
c           then window each and combine to make a windowed
c           time series
c-----
            do 7 i=1,NP
                ftr(i)=0.
    7       continue
            pi2=tupi/4.
            do 110 i=i1,i2
                fq=df*(i-1)
                mm=int(rlngth(i)/dt)+1
                ll=nd2-mm
                lr=nd2+mm
                rmm=float(mm)
                ww=pi2/rmm
                wdt=tupi*fq*dt
                csw=0.0d0
                csdw=dcos(dble(ww))
                snw=-1.0d0
                sndw=dsin(dble(ww))
                csr=dcos(dble(float(ll-1)*wdt))
                csdr=dcos(dble(wdt))
                snr=dsin(dble(float(ll-1)*wdt))
                sndr=dsin(dble(wdt))
                if(ll.le.0) ll=1
                if(lr.gt.n) lr=n
                do 111 j=ll,lr
                    ttt=real(zc(i))*csr-aimag(zc(i))*snr
                    ttt=2.*df*ttt
                    ftr(j)=ftr(j)+ttt*csw**ipwr
                    sntp=snw
                    snw=snw*csdw+sndw*csw
                    csw=csw*csdw-sntp*sndw
                    sntp=snr
                    snr=snr*csdr+sndr*csr
                    csr=csr*csdr-sntp*sndr
  111           continue
  110       continue
        endif
        tt=0.
        do 115 i=1,npf
            paf(i)=ftr(i)
  115   continue
c
c********** at this point, the windowed paf function can ***********
c             be plotted.  the function is stored in array
c                              paf(NPf)
c*******************************************************************
        if(.not. doauto)then
            iflag=1
            call swplot(iflag,info)
            if(info(1:1).eq.'n') go to 150
            go to 15
        endif
 150  continue
c
c  compute spectrum of the isolated seismogram
c
        do 155 i=1,n
            ttt=ftr(i)
            zc(i)=cmplx(ttt,0.)
  155   continue
        call zfour(zc,n,-1,dt,df)
c
c  shift the peak of the paf back to the original position
c
        do 22 i=2,nd2
            dd=-3.1415926*(i-1)
            shift=cmplx(0.,dd)
            shift=cexp(shift)
            j=n-i+2
            zc(i)=zc(i)/shift
            zc(j)=conjg(zc(i))
   22   continue
        nb=nbx
        return
        end

        subroutine gpplot(info)
c
c
c       This program plots group and phase velocity dispersion
c       curves, with automatic scaling of the axes.  Intended
c         to be used in conjunction with program 'match'.
c
        parameter (LIN=5,LOT=6,LER=0)
        parameter(NP2=65540)
        common /dpplot/ xgo(NP2),ygo(NP2),xgn(NP2),ygn(NP2),
     $                xpo(NP2),ypo(NP2),xpn(NP2),ypn(NP2),
     $                ddf,moo,mpn,ml,mh,in1,in2,nmod,m
        common/vwind/i1,i2,xmx,ii1,ii2
        character*2 info
        character*1 c

        common/auto/doauto,iter,itermx,lshift,phasecor
        logical doauto
        integer iter, itermx
        logical lshift
        real phasecor
c-----
c       set up axis dimensions for plot and position of plot
c-----
        xor=2.
        yor=2.0
        xm=5.
        ym=4.
        xm2=xm/2.0
        ym2=ym/2.0
        dym=ym/4.0
        mtx=12
c-----
c       determine the extremes in velocity and period values
c-----
        permin=9999.
        permax=-9999.
        velmin=9999.
        velmax=-9999.
        n1 = ml 
        n2 = mh 
        imod = iabs(nmod)
        kmod=imod
        if(imod.lt.3)imod=3
        if(nmod.eq.3)then
            call mxval(xgo,1,moo,permin,permax)
            call mxval(ygo,1,moo,velmin,velmax)
        else
            call mxval(xgo,n1,n2,permin,permax)
            call mxval(ygo,n1,n2,velmin,velmax)
        endif
            call mxval(xgn,n1,n2,permin,permax)
            call mxval(ygn,n1,n2,velmin,velmax)
            if(kmod.lt.3) go to 5
            call mxval(xpn,n1,n2,permin,permax)
            call mxval(ypn,n1,n2,velmin,velmax)
    5   continue
        if(imod.ne.3)then
            call mxval(xpo,n1,n2,permin,permax)
            call mxval(ypo,n1,n2,velmin,velmax)
        endif
c-----
c       set up scaling for axes
c-----
c-----
c       set up scaling for axes
c-----
c       x -axis is logarithmic
c           choose values that are reasonable
c-----
        if(permin .lt. 1.0)then
            fct =  10.0**int(alog10(permin)-1.0)
        else
            fct =  10.0**int(alog10(permin))
        endif
        permin = float(int(permin/fct))*fct
        if(permax .lt. 1.0)then
            fct = 10.0**int(alog10(permax)-1.0)
        else
            fct = 10.0**int(alog10(permax))
        endif
        permax = float(int(permax/fct))*fct + fct
        xf=permin
        xl=permax
          permin=alog10(permin)
          permax=alog10(permax)
c-----
c       y -axis is linear
c-----
        velmin=float(int(velmin))
        velmax=float(int(velmax))+1.
        scx=xm/(permax-permin)
        scy=ym/(velmax-velmin)
        call gframe(1)
        call plot(xor,yor,-3)
        call factor(1.0)
        call plot(xm,0.0,2)
        call plot(xm,ym,2)
        call plot(0.0,ym,2)
        call plot(0.0,0.0,2)
        call symbol(+0.5,ym+0.2,0.10,'NEW',0.0,+3)
        call plot(1.0,ym+0.2,3)
        call plot(1.5,ym+0.2,2)
        call symbol(2.5,ym+0.2,0.10,'OLD',0.0,+3)
        call plot(3.0,ym+0.2,3)
        call plotd(3.5,ym+0.2,21,0.10)
c-----
c       plot data
c-----
        ioldnw = 0
        if(kmod.eq.1)ioldnw=1
        if(nmod.eq.3)then
            call pltvel(xgo,ygo,1,moo,velmin,velmax,permin,permax
     1          ,scx,scy,ioldnw)
        else
            call pltvel(xgo,ygo,n1,n2,velmin,velmax,permin,permax
     1          ,scx,scy,ioldnw)
        endif
        if(kmod.eq.1) go to 10
            ioldnw = 1
            call pltvel(xgn,ygn,n1,n2,velmin,velmax,permin,permax
     1          ,scx,scy,ioldnw)
        if(kmod.eq.2) go to 10
            call pltvel(xpn,ypn,n1,n2,velmin,velmax,permin,permax
     1          ,scx,scy,ioldnw)
        if(imod.ne.3)then
            ioldnw = 0
            call pltvel(xpo,ypo,n1,n2,velmin,velmax,permin,permax
     1          ,scx,scy,ioldnw)
        endif
   10   continue
        dgp=(velmax-velmin)/4.0
        gvmx = abs(velmax)
        if(abs(velmin).gt.gvmx)gvmx = abs(velmax)
        ng = 1
        if(gvmx.gt.0)then
            ng = alog10(gvmx)+1
        endif
c-----
c       put up tics along y-axis
c-----
        do 105 i=1,5
            yy=velmin+(i-1)*dgp
            yinch=(i-1)*dym
            call symbol(-0.050,yinch,0.10,char(13),90.0,-1)
            call number(-0.5-ng*0.1,yinch-0.05,0.10,yy,0.0,2)
  105   continue
c-----
c       put up tics along x-axis
c-----
          call xaxlog(xf,xl,xm)
        call plot(-xor,-yor,-3)
        if(kmod.eq.2) return
        if(kmod.gt.2) go to 900
c-----
c           use cursors for group velocity spline edit
c-----
   51   continue
        do 50 i=1,2
        call gwrtxt(1.0,1.9-i*0.3,'enter edit point (any key)',0)
        call plot(xor,yor,-3)
        call currxy(xx,yy,c)
        xpg=10**(permin+xx/scx)
        fpg=1./xpg
        fmx=1./xgo(1)
        dfp=(fmx-fpg)/ddf
        ifp=nint(dfp)+1
        if(i.eq.1) then
            in1=ifp
        else
            in2=ifp
        endif
        call plot(-xor,-yor,-3)
  50    continue
        if(in1 .eq. in2)then
            call gwrtxt(1.0,1.0,'Choose TWO DIFFERENT POINTS',0)
        endif
        if(in1.eq.in2)goto 51
        if(in1.gt.in2) then
            itp=in1
            in1=in2
            in2=itp
        endif
        return
 900    continue
        if(nmod.gt.0)then
                if(.not.doauto)then
                    call gwrtxt(1.0,1.0,
     1              'recalculate spline? (y/n):',0)
                    call iyesno(info)
                else
                    info = 'n'
                endif
            else
                if(.not.doauto)then
                    call gwrtxt(1.0,1.0,
     1          'Enter 1,2,3 Top=1, Middle=2, Bottom=3:',0)
                    call get123(info)
                else
                    info = '2'
                endif
            endif
        return
        end
        
        subroutine get123(info)
        parameter (LIN=5,LOT=6,LER=0)
        character*2 info

        common/auto/doauto,iter,itermx,lshift,phasecor
        logical doauto
        integer iter, itermx
        logical lshift
        real phasecor

 1000   continue
            if(.not. doauto)then
                call grdtxt(info,2)
                if(info(1:1).eq.'1' .or. info(1:1).eq.'2'
     1              .or. info(1:1).eq.'3')return
            else
                info(1:1)='2'
                return
            endif
        goto 1000
        end

        subroutine mxval(x,i1,i2,xmin,xmax)
c-----
c       subroutine to determine maximum and minimum values
c       in array, within initial given limits xmin, xmax
c-----
        real*4 x(*)
        do 100 i=i1,i2
            if(x(i).gt.xmax)xmax = x(i)
            if(x(i).lt.xmin)xmin = x(i)
  100   continue
        return
        end

        subroutine pltvel(x,y,n1,n2,ymin,ymax,xmin,xmax
     1          ,scx,scy,ioldnw)
c-----
c       plot x,y array using dashed/solid line
c-----
c       x   - array of abscissa
c       y   - array of ordinates
c       n1  - first point to plot
c       n2  - last point to plot
c       xmin    - value of minimum x on plot
c       xmax    - value of maximum x on plot
c       ymin    - value of minimum y on plot
c       ymax    - value of maximum y on plot
c       ioldnw  - 0 dashed line
c           - 1 solid line
c-----
        real*4 x(*),y(*)
        do 100 i=n1,n2
            xx=(alog10(x(i))-xmin)*scx
            yy=(y(i)-ymin)*scy
            if(i.eq.n1)call plot(xx,yy,3)
            if(ioldnw.eq.1)then
                call plot(xx,yy,2)
            else
                call plotd(xx,yy,21,0.10)
            endif
  100   continue
        return
        end
        
        subroutine xaxlog (xf,xl,xln)
c-----
c       plot x-axis on a logarithmic scale
c
c       xf  - value at left end of axis
c       xl  - value at right end of axis
c       xln - length of axis in inches
c
c       this subroutine calls a subroutine nxttic to tell the position
c       of the next tic value to be plotted
c-----
        sizn = 0.20
        prf = alog10(xf)
        prl = alog10(xl)
        scx = xln/(prl - prf)
c-----
c       scan through possible axis tics, and plot then
c-----
        call plot(0.0,0.0,3)
        xcur = -1.0
 1000   continue
            call nxttic(xf,xl,xcur,ex,iret)
            if(iret.lt.0)goto 1001
            x = (alog10(xcur)-prf)*scx
            call plot(x,0.0,2)
            call plot(x,-sizn/2.0,2)
            call plot(x,0.0,2)
            if(iret.eq.0)then
c-----
c               put in power of 10
c-----
                call symbol(x-sizn,-3.0*sizn,sizn,'10',0.0,2)
                if(ex.lt.0.0 .or.ex.ge.10.0)then
                    xp = x +0.4*sizn
                else
                    xp = 999.0
                endif
                call number(xp,-1.75*sizn,0.75*sizn,ex,0.0,-1)
                call plot(x,0.0,3)
            endif
        goto 1000
 1001   continue
        call plot(xln,0.0,2)
        call plot(0.0,0.0,3)
        return
        end

        subroutine nxttic(xf,xl,xcur,ex,iret)
c-----
c       this subroutine returns the next x-value at which
c       to place a tic on a logarithmic scale
c
c       xf  - lower limit of x-values
c       xl  - upper limit of x-values
c       xcur- present value for tic mark location
c       ex  - exponent of present scale
c       iret- negative -> end of plotting
c           - zero     -> put in exponent of 10
c           - positive -> just an ordinary tic-mark
c-----
        save i1,dx,xmax
c-----
c       initialize parameters if necessary
c-----
        iret = 1
        if(xcur .lt. 0.0)then
            xcur = xf
            xtmp = alog10(xf)
            i1 = int(xtmp)
            xtmp1 = float(i1)
c-----
c           TEST EXPONENT, IF IT IS AN EXACT INTEGER, OR
c           IF IT IS POSITIVE, SET ex = i1, ELSE
c           ex = i1 -1
c-----
            if(xtmp.ne.xtmp1 .and. xtmp.lt.0.0)then
                i1 = i1 - 1
            endif
            if(xtmp.eq.xtmp1)iret=0
            dx = 10.0**i1
            ex = i1
            xmax = 10.0 * dx
        else
            xcur = xcur + dx
        endif
        if(xcur.ge.xmax)then
            i1 = i1 + 1
            ex = i1
            iret = 0
            dx = 10.0**i1
            xmax = 10.0 * dx
        endif
        if(xcur.gt.xl)then
            iret = -1
        endif
        return
        end 

        subroutine swplot(iflag,info)
        parameter (LIN=5,LOT=6,LER=0)
        parameter(NP=131072)
        character*2 info
        character*2 c
        common /pfplot/ paf(NP),waf(NP),npf,nb,pdt,pdf,nfct,ymmax
        character txt1*80

        common/auto/doauto,iter,itermx,lshift,phasecor
        logical doauto
        integer iter, itermx
        logical lshift
        real phasecor

        n=npf
        dt=pdt
        xscl=7.
        yscl=2.0
        xn=float(n)
        x0=.5
        y0=5.5
        if(iflag.gt.0) y0=2.5
        ymax=0.0
        do 10 i=1,n
            if(abs(paf(i)).gt.ymax)ymax=abs(paf(i))
   10   continue
        if(iflag.le.1) ymmax=ymax
        if(iflag.gt.1)ymax=ymmax
        if(iflag.eq.0) call gframe(1)
c-----
c       plot the trace
c-----
        call newpen(2)
        call plot(x,y,3)
        do 20 i=1,n,nfct
            x=x0+float(i-1)*xscl/xn
            y=y0+yscl*paf(i)/ymax
            if(i.eq.1)then
                call plot(x,y,3)
            else
                call plot(x,y,2)
            endif
   20   continue
c-----
c       display various prompts according to the value of iflag
c-----
        call newpen(1)
        if(iflag.eq.0)then
            call plot(x0,0.8,3)
        elseif(iflag.eq.2)then
            call plot(x0,1.0,3)
        elseif(iflag.eq.3)then
            call plot(x0,1.5,3)
        else
c-----
c                       TIC MARK AT ZERO LAG
c-----
            xmid = x0 + 0.5*xscl
            call plot(xmid,0.2,3)
            call plot(xmid,0.6,2)
            call plot(x0,1.0,3)
c-----
c                       TIC MARK AT +- WINDOW
c-----
            dx = nb*xscl/xn
            xx1 = x0 + 0.5*xscl
            call plot(xx1 + dx,y0+0.5,3)
            call plot(xx1 + dx,y0+0.2,2)
            call plot(xx1 - dx,y0+0.5,3)
            call plot(xx1 - dx,y0+0.2,2)
        endif
        if(iflag.eq.0) go to 9999
        if(iflag.lt.2) go to 23
        if(iflag.eq.3) go to 22
        if(.not. doauto)then
        call gwrtxt(1.0,1.0,'plot residual seismogram? (y/n):',0)
        call iyesno(info)
        else
            info(1:1) = 'n'
        endif
        go to 999
  22  continue
        if(.not. doauto)then
            call gwrtxt(1.0,1.0,'replot seismograms? (y/n):',0)
            call iyesno(info)
        else
            info(1:1) = 'n'
        endif
        go to 999
  23  continue
        xinc=(nb-1)*dt*2.
        call plot(x0,y0-1.0,3)
        write(txt1,25)xinc
  25  format(' width=',f6.1,' sec')
        call gwrtxt(1.0,y0-1.0,txt1,0)
        if(.not. doauto)then
            call gwrtxt(1.0,y0-1.3,'recalculate width? (y/n):',0)
            call iyesno(info)
        else
            info(1:1) = 'n'
        endif
        if(info(1:1).ne.'n') then
            call curixy(ix,iy,c)
            xx = ix/1000.0
            yy = iy/1000.0
            dx = abs(x0 + 0.5*xscl - xx)
            nb = int(xn*dx/xscl)
C       write(6,*)'swplot: ix,iy,xx,yy,dx,xscl,xn,nb',
C     1     ix,iy,xx,yy,dx,xscl,xn,nb
        endif
 999  call gframe(1)
9999  continue
        return
        end

        subroutine iyesno(info)
        parameter (LIN=5,LOT=6,LER=0)
        character*2 info
 1000   continue
            call grdtxt(info,2)
            if(info(1:1).eq.'Y')info(1:1)='y'
            if(info(1:1).eq.'N')info(1:1)='n'
            if(info(1:1).eq.'y' .or. info(1:1).eq.'n')return
        goto 1000
        end

        subroutine getspc(fsac,ierr)
        parameter (LER=0,LIN=5,LOT=6)
        parameter(NP=131072,NP2=65540)
        implicit complex (z)
        dimension z0(NP),x(NP2),y(NP2),yf(NP2),yg(NP2),c(NP2)
        dimension zsav(NP)
        common/tdum/z0,x,y,yf,yg,c,zsav,df,n,n2,tt0,ddist,ns
        character*256 namev,names,namer
        common/fname/namev,names,namer
        character*80 name
        character sta*8, comp*8, cdate*12
        common/spcsav/nn,npts,sta,comp,dist,deg,az,baz,t0,dt,cdate
        common /dpplot/ xgo(NP2),ygo(NP2),xgn(NP2),ygn(NP2),
     $                xpo(NP2),ypo(NP2),xpn(NP2),ypn(NP2),
     $                ddf,moo,mpn,ml,mh,in1,in2,nmod,m
        common /pfplot/ paf(NP),waf(NP),npf,nb,pdt,pdf,nfct,ymmax
        common /tvar/ zc(NP),ddt,ntv,nti1,nti2,ntj1,ntj2
        real tarr(NP)

        common/auto/doauto,iter,itermx,lshift,phasecor
        logical doauto
        integer iter, itermx
        logical lshift
        real phasecor

        character fsac*(*)
c
c       input event spectrum
c
c-----
c       silent vector move to get prompt at upper left screen
c-----
        name = fsac
            ls=lgstr(name)
        if(.not. doauto)then
            call gframe(1)
            call gwrtxt(1.0,7.0,'Spectrum file :',0)
            lnm = 40
            call gwrtxt(4.0,7.0,name(1:ls),0)
        endif
        namer=name(1:ls)//'r'
        names=name(1:ls)//'s'

        call gettrc(name(1:ls),n,n21,npts,dist,deg,az,baz,t0,dt,sta,
     1      comp,cdate,z0,tarr,NP,ierr)
        write(6,*)'gettrc:    n',   n
        write(6,*)'gettrc:  n21', n21
        write(6,*)'gettrc: npts',npts
        write(6,*)'gettrc: dist',dist
        write(6,*)'gettrc:  deg', deg
        write(6,*)'gettrc:   az',  az
        write(6,*)'gettrc:  baz', baz
        write(6,*)'gettrc:   t0',  t0
        write(6,*)'gettrc:   dt',  dt
        write(6,*)'gettrc:  sta', sta
        write(6,*)'gettrc: comp',comp
        write(6,*)'gettrc:cdate',cdate
        if(n.gt.NP)then
            call pend()
            write(LER,*)'Time series points exceed array size'
            write(LER,*)'Time series points = ',n
            write(LER,*)'Array size         = ',NP
            write(LER,*)'Run xspdec to reduce time series points'
            ierr = -1
        endif
        ddt = dt
        tt0 = t0
        ddist = dist
        if(ierr.lt.0)then
c----BAD USE GWRTXT
            write(LER,*)'ierr:',ierr
            write(LER,*)'Error on read: subroutine gettrc'
            if(.not. doauto)then
            call gwrtxt(1.0,6.5,'Error on reading sac file :',0)
            call gwrtxt(1.0,6.0,'      click to end        :',0)
            endif
            return
        endif
c-----
c       end of spectrum read
c
c       z0 - observed spectrum
c       n  - number of samples
c       dt - time sampling increment
c       df - frequency sampling increment
c-----
        do 16 i=1,n
            zsav(i) = z0(i)
   16   continue
        df = 1./(n*dt)
        pdt=dt
        pdf=df
        return
        end

        subroutine savspc(iret,permin,permax)
        parameter (LIN=5,LOT=6,LER=0)
        parameter(NP=131072,NP2=65540)
        implicit complex (z)
        dimension z0(NP),x(NP2),y(NP2),yf(NP2),yg(NP2),c(NP2)
        dimension zsav(NP)
        common/tdum/z0,x,y,yf,yg,c,zsav,df,n,n2,tt0,ddist,ns
        character*256 namev,names,namer
        common/fname/namev,names,namer
        character sta*8, comp*8, cdate*12
        common/spcsav/nn,npts,sta,comp,dist,deg,az,baz,t0,dt,cdate
        character*2 info, info1
        common /dpplot/ xgo(NP2),ygo(NP2),xgn(NP2),ygn(NP2),
     $                xpo(NP2),ypo(NP2),xpn(NP2),ypn(NP2),
     $                ddf,moo,mpn,ml,mh,in1,in2,nmod,m
        common /pfplot/ paf(NP),waf(NP),npf,nb,pdt,pdf,nfct,ymmax
        common /tvar/ zc(NP),ddt,ntv,nti1,nti2,ntj1,ntj2
        common/vwind/i1,i2,xmx,ii1,ii2
        common/iband/ib1,ib2,ib3,ib4
        common/srfdat/iunit,ifrpr,ilr,imd
c-----
c       arrays to store the period, phase velocity, group velocity and
c               spectral amplitude from the Fourier spectrum
c-----
        real*4 xp(NP2),cp(NP2),up(NP2),ap(NP2)
c-----
c       array to store the desired periods between [permin,permax]
c-----
        real*4 gp(NP2)

        character txt1*80
        real tarr(NP)
        real evla, evlo, evdp, stla, stlo, stel, origtime
        integer ntimes(6)
        common /shdr/ evla,evlo,evdp,stla,stlo,stel,
     &               beg,ntimes,origtime

        common/auto/doauto,iter,itermx,lshift,phasecor
        logical doauto
        integer iter, itermx
        logical lshift
        real phasecor
c---- 
        tupi=2*3.141592654
c
c********** at this point, the time plotter can be called **********
c             to plot the original seismogram and the isolated
c             seismogram. the files are stored in: paf(NPf)
c
c*******************************************************************
c
  770   continue
            f1 = ib1*df
            f2 = ib2*df
            f3 = ib3*df
            f4 = ib4*df
            if(.not.doauto)then
            call gwrtxt(1.0,6.5,
     1      'Are original signal and residual to be filtered ',0)
            call gwrtxt(1.0,6.0,' for display in the window',0)
            write(txt1,10)f1,f2,f3,f4
   10       format('(f1',e10.3,' f2 ',e10.3,' f3 ',e10.3,' f4 ',
     1            e10.3,' )y/n?')
            call gwrtxt(1.0,5.5,txt1,0)
            call gwrtxt(1.0,5.0,
     1      ' Note Residual Spectra Saved is NOT FILTERED',0)
            call gwrtxt(1.0,4.5,'RESPONSE:',0)
            call iyesno(info1)
            else
                info1 = 'n'
            endif
            do 780 i=1,n
                if(info1(1:1) .eq. 'n')then
                    paf(i) = real(zsav(i))
                else
                    paf(i) = real(z0(i))
                endif
  780       continue
            iflag=0
c-----
c           plot seismogram
c-----
            if(.not. doauto)call swplot(iflag,info)
            do 735 i=1,n
                paf(i)=real(zc(i))
  735       continue
            iflag=2
c-----
c           plot residual seismogram
c-----
            if(.not. doauto)then
                call swplot(iflag,info)
                do 800 i=1,n
                    if(info1(1:1) .eq. 'n')then
                        paf(i) = real(zsav(i))
                    else
                        paf(i) = real(z0(i))
                    endif
  800           continue
                iflag=0
                call swplot(iflag,info)
                do 810 i=1,n
                    paf(i)=paf(i)-real(zc(i))
  810           continue
                iflag=3
                call swplot(iflag,info)
                if(info(1:1).eq.'n') go to 790 
 1003           continue
                call gwrtxt(1.0,7.0,
     1          'enter decimation factor for seismograms :',0)
                call grdtxt(txt1,80)
                read(txt1,'(i5)',err=1003) nfct
                go to 770
            endif
 790        continue
c-----
c           check whether mtv should be reinvoked
c-----
            if(.not. doauto)then
                call gwrtxt(1.0,7.0,
     1          'recompute amplitude correction? (y/n):',0)
                call iyesno(info)
            else
                info = 'n'
            endif
            iret=0
            if(info(1:1).eq.'y')iret=1
            if(iret.eq.1)return
            if(.not. doauto)then
                call gwrtxt(1.0,6.5,
     1          'Writing Velocity and Spectra Results',0)
            endif
c-----
c           output results
c-----
c-----
c       output the isolated spectrum in fort.1 format
c-----
        do 840 i=1,n
            spd=real(zsav(i))-real(zc(i))
            z0(i)=cmplx(spd,0.0)
  840   continue
        ifr=-1
        call zfour(z0,n,ifr,dt,df)
        call zfour(zc,n,ifr,dt,df)
        n21 = n/2 + 1
c-----
c           exit program
c
c           updated group and phase velocity files are stored
c           in unit 3 - original group velocity file name with
c           a 'v' appended to it.
c           internal name in namev
c
c           isolated spectrum file is stored in unit 4 - original
c           spectrum file with an 's' appended to it
c           internal name in names
c
c           residual spectrum file is stored in unit 2 - original
c           spectrum file with an 'r' appended to it
c           internal name in namer
c-----
            open(3,file=namev)
            rewind 3
            open(4,file='disp.out',status='unknown',
     1          form='formatted',access='sequential')
            rewind 4
            nmod=2
c-----
c           output dispersion in MFT96 format
c-----
            n1 = ml 
            n2 = mh 
            iph = 1
            igr = 2
c-----
c       frequency if (j-1)*df
c       k is index into array for group delay
c-----
c       define output arrays
c-----
            ndat = 0
            do 650 i=n2,n1,-1
                j=i2-i+1
                k=m-i+1
                if(ifrpr.eq.0)then
                    frpr = tupi/yg(k)
                else
                    frpr = yg(k)/tupi
                endif
                uv = dist/(y(j) + t0)
                du = 1.0
                cv = dist*yg(k)/(x(j)+yg(k)*t0 + phasecor)
                dc = 1.0
            amp = cabs(zc(j))
            jpeak = 1
            ndat = ndat + 1
            xp(ndat) = frpr
            cp(ndat) = cv
            up(ndat) = uv
            ap(ndat) = amp
            if(ilr.eq.1)then
C        write(3,11)imd,frpr,cv,dc,dist,amp,evla,evlo,stla,stlo,jpeak
        write(3,12)imd,frpr,uv,du,dist,amp,evla,evlo,stla,stlo,jpeak
            else if(ilr.eq.2)then
C        write(3,21)imd,frpr,cv,dc,dist,amp,evla,evlo,stla,stlo,jpeak
        write(3,22)imd,frpr,uv,du,dist,amp,evla,evlo,stla,stlo,jpeak
            endif
  650   continue
c-----
c    now interpolate on the desired periods in the range [permin,permax]
c    after getting the array of periods for interpolation
c    output in order of decreasing period
c-----
        call getper(gp,NP2,permin,permax,nper)
c-----
c       do a gruesome linear search
c-----
        do 652 i=nper,1,-1
            per = gp(i)
c-----
c       find it
c-----
            do 653 j=1,ndat-1
                if(xp(j).ge.per .and.xp(j+1).lt.per)then
                    p = (per - xp(j))/(xp(j+1)-xp(j))
                    uv = p*up(j+1) + ( 1.0 - p)*up(j)
                    cv = p*cp(j+1) + ( 1.0 - p)*cp(j)
                    amp = p*ap(j+1) + ( 1.0 - p)*ap(j)
            if(ilr.eq.1)then
        write(3,12)imd,per,uv,du,dist,amp,evla,evlo,stla,stlo,jpeak
        write(4,13)imd,per,uv,du
        write(4,14)imd,per,cv,dc
            else if(ilr.eq.2)then
        write(3,22)imd,per,uv,du,dist,amp,evla,evlo,stla,stlo,jpeak
        write(4,23)imd,per,uv,du
        write(4,24)imd,per,cv,dc
            endif
                endif
  653       continue
  652   continue

   11   format('MFT96 L C T ',i2,' ',4g10.5,e11.4,4f12.6, 
     1          ' 0',i3,' COMMENT')
   12   format('MFT96 L U T ',i2,' ',4g10.5,e11.4,4f12.6, 
     1          ' 0',i3,' COMMENT')
   21   format('MFT96 R C T ',i2,' ',4g10.5,e11.4,4f12.6, 
     1          ' 0',i3,' COMMENT')
   22   format('MFT96 R U T ',i2,' ',4g10.5,e11.4,4f12.6, 
     1          ' 0',i3,' COMMENT')
   13   format('SURF96 L U T ',i2,1x,4g10.5)
   23   format('SURF96 R U T ',i2,1x,4g10.5)
   14   format('SURF96 L C T ',i2,1x,4g10.5)
   24   format('SURF96 R C T ',i2,1x,4g10.5)
        n21 = n/2 + 1
c-----
c       save isolated spectrum
c-----
        call puttrc(names,n,n21,npts,dist,deg,az,baz,t0,dt,sta,
     1      comp,cdate,zc,tarr,NP,ierr,permin,permax)
c-----
c       save residual spectrum
c       Note - we do not want to carry the period bounds here
c       so that we could look at higher modes
c-----

        call puttrc(namer,n,n21,npts,dist,deg,az,baz,t0,dt,sta,
     1      comp,cdate,z0,tarr,NP,ierr,-12345.0,-12345.0)
c-----
c       clean up file i/o
c-----
        close(3)
        close(4)
        if(.not. doauto)then
            call gwrtxt(1.0,6.0,'program ended',0)
        endif
        return
        end

        subroutine pmatch()
        parameter (LIN=5,LOT=6,LER=0)
        parameter(NP=131072,NP2=65540)
        implicit complex (z)
        dimension z0(NP),x(NP2),y(NP2),yf(NP2),yg(NP2),c(NP2)
        common/tdum/z0,x,y,yf,yg,c,zsav,df,n,n2,t0,dist,ns
        dimension zsav(NP)
        character*2 info
        common /dpplot/ xgo(NP2),ygo(NP2),xgn(NP2),ygn(NP2),
     $                xpo(NP2),ypo(NP2),xpn(NP2),ypn(NP2),
     $                ddf,moo,mpn,ml,mh,in1,in2,nmod,m
        common /pfplot/ paf(NP),waf(NP),npf,nb,pdt,pdf,nfct,ymmax
        common /tvar/ zc(NP),dt,ntv,nti1,nti2,ntj1,ntj2
        common/vwind/i1,i2,xmx,ii1,ii2
        character txt1*80

        common/auto/doauto,iter,itermx,lshift,phasecor
        logical doauto
        integer iter, itermx
        logical lshift
        real phasecor

        character ssta*8, ccomp*8, ccdate*12
        common/spcsav/nnn,nnpts,ssta,ccomp,ddist,ddeg,aaz,
     1      bbaz,tt0,ddt,ccdate

        common/per/ppermin,ppermax

        common/break/nbpts
        integer nbpts

c
c       calculate pseudo-autocorrelation function (zc) by removing
c       estimate of phase (x) from observed spectrum (z0)
c
        tupi = 2.* 3.141592654
        dw = tupi*df
        ts=float(n)*dt/2.
  340   continue
            zc(1)=cmplx(0.,0.)
            do 345 i=2,n2
                zr=cexp(cmplx(0.,x(i)))
                zc(i)=z0(i)*zr
                j=n-i+2
                zc(j)=conjg(zc(i))
  345       continue
c-----
c           shift zero lag to center of array
c-----
            do 400 i=2,n2
                tp=ts*float(i-1)*dw
                zp=cexp(cmplx(0.,-tp))
                zc(i)=zc(i)*zp
                j=n-i+2
                zp=conjg(zp)
                zc(j)=zc(j)*zp
  400       continue
            zc(n2)=cmplx(real(zc(n2)),0.)
c-----
c       fourier transform zc() to time domain (ifr=1)
c-----
        ifr=1
        call zfour(zc,n,ifr,dt,df)
        nfct = 1
  405   continue
        do 410 i=1,n
            paf(i)=real(zc(i))
  410   continue
        npf=n
c
c ******* at this point use plotter routine to plot time domain ********
c         pseudo-autocorrelation function... be sure plotter has
c         cursor capability to return the width of the time domain
c         window function of the correlation function.  this will
c         be symmetric about zero-lag (center of array), and in
c         terms of the one-sided distance (nb) from zero-lag.
c
c         the values for the pseudo-autocorrelation function are stored
c         in array paf(NPf).
c
c***********************************************************************
        if(.not. doauto)then
            iflag=0
            call swplot(iflag,info)
        endif
        j1=n2
        j2=n2
        nt=1
c-----
c       time window the pseudo-autocorrelation function.
c       this step is not performed when exiting the program
c       to time variable filter.
c-----
        do 424 i=1,n
            waf(i)=real(zc(i))
  424   continue
c-----
c       SAVE the shifted function
c-----
        if(lshift)then
            call putsac('shift.sac',n,ddist,ddeg,aaz,bbaz,-ts,dt,
     1      ssta,'Shifted ',ccdate,waf,ppermin,ppermax)
        endif
        call band(zc,n,j1,j2,nb,nt)
        do 425 i=1,n
            paf(i)=real(zc(i))
  425   continue
c-----
c       SAVE the windowed function
c-----
        
        if(lshift)then
            call putsac('wshift.sac',n,ddist,ddeg,aaz,bbaz,-ts,dt,
     1      ssta,'WinShift',ccdate,paf,ppermin,ppermax)
        endif
c
c ******** at this point, the windowed pseudo-autocorrelation *********
c            function can be plotted if desired. the windowed function
c            is in array:        paf(NPf)
c
c**********************************************************************
        if(.not. doauto)then
            iflag=1
            call swplot(iflag,info)
            if (info(1:1).eq.'n') go to 427
            do 426 i=1,n
                zc(i)=cmplx(waf(i),0.0)
  426       continue
            go to 405
        endif
  427  continue
c-----
c       after final iteration, exit the program
c-----
        if(.not.doauto)then
            call gwrtxt(1.0,7.0,'satisfied? (y/n):',0)
            call iyesno(info)
        else
            if(iter.lt.itermx)then
                info = 'n'
            else
                info = 'y'
            endif
        endif
        if(info(1:1).eq.'n') go to 430
        do 428 i=1,n
            zc(i)=cmplx(waf(i),0.0)
  428   continue
c
c       fourier transform back to frequency domain and reshift
c       zero lag back to its correct position.
c
  430   ifr=-1
        call zfour(zc,n,ifr,dt,df)
        zc(1)=cmplx(0.,0.)
        do 450 i=2,n2
            tp=ts*float(i-1)*dw
            zp=cexp(cmplx(0.,tp))
            zc(i)=zc(i)*zp
            j=n-i+2
            zp=conjg(zp)
            zc(j)=zc(j)*zp
  450   continue
        zc(n2) = cmplx(real(zc(n2)),0.0)
c-----
c       exit program here
c-----
        if(info(1:1).eq.'y') go to 635
c-----
c       unwrap the phase spectrum by numerical differentiation and
c       integration.  check the subroutine for the location of
c       plotting the group delay error and interpolating it for
c       getting rid of "spikes" due to spectral nulls.
c
c       unwrapped phase is returned in (yf).
c-----
  510   call unwrap(yf,dw,i2,j1,j2,m)
        j=1
        do 500 i=i1,i2
            yf(j)=yf(i)
            j=j+1
  500   continue
c-----
c       enter number of breakpoints for least squares spline fit
c       to the unwrapped phase.
c-----
        nfd=0
 1002   continue
        if(.not. doauto)then
            call gwrtxt(1.0,0.5,
     1          'Enter # breakpoints for correction spline:',0)
            call grdtxt(txt1,80)
            read(txt1,'(i5)',err=1002)nbp
        else
            if(iter.eq.0)then
                nbp = 2
            else
                nbp = nbpts
            endif
            iter = iter + 1
        endif
        if(nbp.lt.0) nbp=m-4
        call spline(yg,yf,c,m,nbp)
        call bfit(yg,yf,c,m,nfd)
c-----
c       calculate corrected phase velocity and store in xpn,ypn
c       (period, velocity).
c-----
        do 520 i=1,m
            j=i2-i+1
            k=m-i+1
            w=float(j-1)*dw
            xpo(i)=xpn(i)
            ypo(i)=ypn(i)
            xpn(i)=tupi/w
            ypn(i)=dist*w/(x(j)-yf(k)+w*t0)
  520   continue
c-----
c       differentiate the phase to calculate corrected group delay.
c       store in xgn,ygn.
c-----
        nfd=1
        call bfit(yg,yf,c,m,nfd)
        do 530 i=1,m
            j=i2-i+1
            k=m-i+1
            w=float(j-1)*dw
            xgo(i)=xgn(i)
            ygo(i)=ygn(i)
            xgn(i)=tupi/w
            ygn(i)=dist/(y(j)-yf(k)+t0)
  530   continue
        nmod=4
c ********at this point, the corrections to phase and group velocity ***
c           can be plotted. use automatic scale x-log, y-linear plotter.
c           plot:
c                      xgo,ygo - dashed curve
c                      xpo,ypo - dashed curve
c                      xgn,ygn - solid curve
c                      xpn,ypn - solid curve
c
c           plot all curves between array limits ml and mh.
c
c **********************************************************************
        if(.not. doauto)then
            call gpplot(info)
        endif
c-----
c       if the spline correction is undesirable, recalculate the fit
c-----
        if(info(1:1).eq.'n') go to 540
        do 535 i=1,m
            xpn(i)=xpo(i)
            ypn(i)=ypo(i)
            xgn(i)=xgo(i)
            ygn(i)=ygo(i)
  535   continue
        go to 510
  540  continue
c-----
c       option to specify the arbitrary 2 pi radian phase difference
c-----
        if(.not. doauto)then
            call gwrtxt(1.0,0.5,
     1  'Plot Phase Velocity at +- 2 pi radian change (y/n)?:',0)
            call iyesno(info)
        else
            info = 'n'
        endif
        iphas = 0
        if(info(1:1).eq.'y')then
            iphas = 1
            do 541 i=1,m
                j = i2 -i + 1
                k = m-i+1
                w = float(j-1)*dw
c-----
c       save some old values
c-----
                yf(i) = ygo(i)
c-----
c       generate new plot values
c-----
                xgn(i) = xpn(i)
                xgo(i) = xpn(i)
                phase = dist*w/ypn(i)
                ygo(i) = ypn(i)
                ygn(i) = dist*w/(phase + 6.2831853)
                ypn(i) = dist*w/(phase - 6.2831853)
  541       continue
            nmod = -3
            call gpplot(info)
c-----
c       reconstruct old values of phase and groupg velocity curves
c       correct phase if necessary
c-----
            do 542 i=1,m
                ygn(i) = yf(i)
                ypn(i) = ypo(i)
  542       continue
            do 543 j=i1,i2
                if(info(1:1).eq.'1')then
                    x(j) = x(j) - 6.2831853
                elseif(info(1:1).eq.'3')then
                    x(j) = x(j) + 6.2831853
                endif
  543       continue
        endif
        if(iphas.eq.1)goto 510
c-----
c       permanently update phase (x) and group delay (y) file
c-----
        nfd=0
        call bfit(yg,yf,c,m,nfd)
        do 550 j=i1,i2
            x(j)=x(j)-yf(j-i1+1)
  550   continue
        nfd=1
        call bfit(yg,yf,c,m,nfd)
        do 600 j=i1,i2
            y(j)=y(j)-yf(j-i1+1)
  600   continue
c-----
c       start next iteration
c-----
        go to 340
  635   continue
        return
        end

        subroutine timser
        implicit none
        integer NP, NP2
        parameter(NP=131072,NP2=65540)
        common/tdum/z0,x,y,yf,yg,c,zsav,df,n,n2,t0,dist,ns
        complex z0(NP)
        real x(NP2),y(NP2),yf(NP2),yg(NP2),c(NP2)
        complex zsav(NP)
        real df, t0, dist
        integer n, n2, ns
        common /tvar/ zc(NP),dt,ntv,nti1,nti2,ntj1,ntj2
        complex zc
        real dt
        integer ntv,nti1,nti2,ntj1,ntj2
        common/vwind/i1,i2,xmx,ii1,ii2
        real xmx
        integer i1, i2, ii1, ii2
        integer nt, i, j, ifr
        complex zr

        nt=0
        call band(zc,n,ii1,ii2,ns,nt)
        do 775 i=i1,i2
            zr=cexp(cmplx(0.,-x(i)))
            zc(i)=zc(i)*zr
  775   continue
        do 750 i=2,n2
            j=n-i+2
            zc(j)=conjg(zc(i))
  750   continue
        ifr=1
        call zfour(  zc,n,ifr,dt,df)
        return
        end

        subroutine gedit()
        implicit none
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)
        integer NP, NP2
        parameter (NP=131072,NP2=65540)
        character*2 info
        common /dpplot/ xgo(NP2),ygo(NP2),xgn(NP2),ygn(NP2),
     1                xpo(NP2),ypo(NP2),xpn(NP2),ypn(NP2),
     2                ddf,moo,mpn,ml,mh,in1,in2,nmod,m
        real xgo, ygo, xgn, ygn, xpo, ypo, xpn, ypn, ddf
        integer moo, mpn, ml, mh, in1, in2, nmod, m
        common/tdum/z0,x,y,yf,yg,c,zsav,df,n,n2,t0,dist,ns
        real x(NP2),y(NP2),yf(NP2),yg(NP2),c(NP2)
        complex z0(NP), zsav(NP)
        real  df, t0, dist
        integer n, n2, ns
        common/vwind/i1,i2,xmx,ii1,ii2
        integer i1, i2, ii1, ii2
        real xmx
        integer i, j
        real tupi, fr, slp

        tupi=3.141592654
        call gwrtxt(1.0,6.5,'edit group velocity? (y/n):',0)
        call iyesno(info)
        if(info(1:1).eq.'n') return
   5    do 10 i=1,m
            j=i2-i+1
            fr=float(j-1)*df
            xgo(i)=1./fr
            xgn(i)=xgo(i)
            ygo(i)=dist/(y(j)+t0)
            ygn(i)=ygo(i)
  10    continue
  15    nmod=1
        call gpplot(info)
        slp=(ygo(in2)-ygo(in1))/float(in2-in1)
        do 20 i=in1,in2
            ygn(i)=slp*float(i-in1)+ygo(in1)
  20    continue
        nmod=2
        call gpplot(info)
        do 25 i=1,m
            ygo(i)=ygn(i)
  25    continue
        call gwrtxt(1.0,0.5,'are more edits desired? (y/n):',0)
        call iyesno(info)
        if(info(1:1).eq.'y') go to 15
        call gframe(1)
        call gwrtxt(1.0,7.0,'redo the edit? (y/n):',0)
        call iyesno(info)
        if(info(1:1).eq.'y') go to 5
        call gwrtxt(1.0,6.5,'save the changes? (y/n):',0)
        call iyesno(info)
        if(info(1:1).eq.'y') then
            do 30 i=1,m
                j=i2-i+1
                y(j)=dist/ygn(i)-t0
  30        continue
        endif
        call gwrtxt(0.0,6.0,'group velocity edit finished',0)
        return
        end

        subroutine gettrc(xfile,n,n21,npts,dist,deg,az,baz,t0,dt,sta,
     1      comp,cdate,z,tarr,nptser,ierr)
c-----
c       get parameters required for processing
c-----
c       xfile   Ch* - name of data file
c       n   I*4 - array size, power of 2
c       n21 I*4 - n/2 + 1
c       npts    I*4 - original number of points (n >= npts)
c       dist    R*4 - epicentral distance (km)
c       deg R*4 - epicentral distance (deg)
c       az  R*4 - src -> rec azimuth
c       baz R*4 - rec -> src azimuth
c       t0  R*4 - time of first sample after origin time
c       dt  R*4 - sample interval
c       sta Ch*8    - station name
c       comp    Ch*8    - station component
c       cdate   Ch*12   - date string
c       z   C*4 - complex Fourier Transform array of samples
c       tarr    R*4 - temporary array  for trace
c       nptser  I*4 - length of z() and tarr() arrays
c       ierr    I*4 - error code
c-----
        implicit none
        character xfile*(*)
        integer n, n21, npts, nptser, ierr
        real dist, deg, az, baz, t0, dt
        character sta*8, comp*8, cdate*12
        complex z(nptser) 
        real tarr(nptser)
        integer ls, i
        real df
        integer lgstr

        call getsac(xfile,npts,dist,deg,az,baz,t0,dt,sta,
     1                 comp,cdate,tarr,ierr)
        ls = lgstr(xfile)
        if(npts.gt.nptser)npts = nptser
          call npow2(npts,n,n21)
        do 100 i=1,n
            if(i.le.npts)then
                z(i) = cmplx(tarr(i),0.0)
            else
                z(i) = cmplx(0.0,0.0)
            endif
  100   continue
        call zfour(z,n,-1,dt,df)
        return
        end
c
        subroutine getsac(name,npts,dist,deg,az,baz,t0,dt,sta,
     1                 comp,cdate,seis,ierr)
c-----
c
c       name    - file name to write
c       n   - number of points in FFT must be power of 2
c       n21 - number of frequencies = n/2 + 1
c       npts    - number of points in original time series
c           - which may have been zero filled to make power of 2
c       dist    - epicentral distance in km
c       deg - epicentral distance in degrees
c       az  - source - receiver azimuth in degrees
c       baz - receiver-source back azimuth
c       t0  - time of first sample after origin
c       dt  - sampling interval
c       sta - C*4 station name string
c       comp    - C*4 component name string
c       cdate   - C*12 date string
c       z   - COMPLEX array of spectra
c       ierr    I   error condition
c-----
        implicit none
        integer NP
        parameter (NP = 131072)
        integer  npts
        real dist, deg, az, baz, t0, dt
        character sta*8, comp*8, cdate*12
        character kstnm*8, kcmpnm*8
        character name*(*)
        real seis(NP)
        integer ierr
*
        common /shdr/ evla,evlo,evdp,stla,stlo,stel,
     &               beg,ntimes,origtime
        real evla, evlo, evdp, stla, stlo, stel, beg, origtime
        integer ntimes(6)

        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*4 chdr(48)

        integer LIN, LOT, LER
        parameter (LIN=5, LOT=6, LER=0)
        integer nerr

        real tp, ts
        integer lgstr
        integer ls


        ls = lgstr(name)
        write(6,*)'FILE: ',name(1:ls)
        call brsac (1,NP,name,seis,ierr)
        if(ierr.lt.0)return
        call getfhv('AZ      ',az,nerr)
        call getfhv('BAZ     ',baz,nerr)
        call getfhv('DIST    ',dist,nerr)
        call getfhv('GCARC   ',deg,nerr)
        call getfhv('DELTA   ',dt, nerr)
        call getnhv('NPTS    ',npts,nerr)
        call getfhv('B       ',beg, nerr)
        call getfhv('O       ',origtime,nerr)
        call getfhv('EVLA    ',evla,nerr)
        call getfhv('EVLO    ',evlo,nerr)
        call getfhv('STLA    ',stla,nerr)
        call getfhv('STLO    ',stlo,nerr)
        call getkhv('KSTNM   ',kstnm,nerr)
        call getkhv('KCMPNM  ',kcmpnm,nerr)
        write(6,'(a)')'  AZ   BAZ     DIST      DEG      DT      NPTS'
        write(6,'(f5.1,f6.1,f13.5,f7.2,f10.5,i7)')
     1      az,baz,dist,deg,dt,npts
        write(6,'(a)')
     1  '   EVLA      EVLO     STLA      STLO      KSTNM   KCMPNM'
        write(6,'(f9.4,f10.4,f9.4,f10.4,1x,a8,1x,a8)')
     1      evla,evlo,stla,stlo,kstnm,kcmpnm
        if(nerr .eq. 0 .and. origtime .ne. -12345)then
            t0 = beg - origtime
        else
            t0 = beg
        end if
        tp = tp - origtime
        ts = ts - origtime
        sta = kstnm(1:8)
        comp = kcmpnm(1:8)
        cdate = ' '
        return
        end

        subroutine puttrc(xfile,n,n21,npts,dist,deg,az,baz,t0,dt,sta,
     1      comp,cdate,z,tarr,nptser,ierr,permin,permax)
c-----
c       get parameters required for processing
c-----
c       xfile   Ch* - name of data file
c       n   I*4 - array size, power of 2
c       n21 I*4 - n/2 + 1
c       npts    I*4 - original number of points (n >= npts)
c       dist    R*4 - epicentral distance (km)
c       deg R*4 - epicentral distance (deg)
c       az  R*4 - src -> rec azimuth
c       baz R*4 - rec -> src azimuth
c       t0  R*4 - time of first sample after origin time
c       dt  R*4 - sample interval
c       sta Ch*8    - station name
c       comp    Ch*8    - station component
c       cdate   Ch*12   - date string
c       z   C*4 - complex Fourier Transform array of samples
c       tarr    R*4 - temporary array  for trace
c       nptser  I*4 - length of z() and tarr() arrays
c       ierr    I*4 - error code
c       permin  R*4 - minimum period for processing of dispersion
c       permax  R*4 - maximum period for processing of dispersion
c-----
        implicit none
        character xfile*(*)
        integer n, n21, npts, nptser, ierr
        real dist, deg, az, baz, t0, dt
        character sta*8, comp*8, cdate*12
        complex z(nptser) 
        real tarr(nptser)
        real permin, permax, df
        integer i

        call zfour(z,n,+1,dt,df)
        do 100 i=1,npts
            tarr(i) = real(z(i))
  100   continue
        call putsac(xfile,npts,dist,deg,az,baz,t0,dt,sta,
     1                 comp,cdate,tarr,permin,permax)
        return
        end
c
        subroutine putsac(name,npts,dist,deg,az,baz,t0,dt,sta,
     1                 comp,cdate,seis,permin,permax)
c-----
c
c       name    - file name to write
c       n   - number of points in FFT must be power of 2
c       n21 - number of frequencies = n/2 + 1
c       npts    - number of points in original time series
c           - which may have been zero filled to make power of 2
c       dist    - epicentral distance in km
c       deg - epicentral distance in degrees
c       az  - source - receiver azimuth in degrees
c       baz - receiver-source back azimuth
c       t0  - time of first sample after origin
c       dt  - sampling interval
c       sta - C*4 station name string
c       comp    - C*4 component name string
c       cdate   - C*12 date string
c       z   - COMPLEX array of spectra
c       permin  R*4 - minimum period for processing of dispersion
c       permax  R*4 - maximum period for processing of dispersion
c-----
        implicit none
        integer MXPTS
        parameter (MXPTS = 65540)
        integer  npts
        real dist, deg, az, baz, t0, dt
        character sta*8, comp*8, cdate*12
        character name*(*)
        real seis(MXPTS)
        real permin, permax
*
        real evla, evlo, evdp, stla, stlo, stel, origtime, beg
        integer ntimes(6)
        common /shdr/ evla,evlo,evdp,stla,stlo,stel,
     &               beg,ntimes,origtime

        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*4 chdr(48)
        integer ierr
        real depmax, depmin, depmen, btime
        integer indmax, indmin

        call scmxmn(seis,npts,depmax,depmin,depmen,indmax,indmin)
        call getfhv('B    ',btime,ierr)
        call setnhv('NPTS', npts, ierr)
        call setfhv('USER1', permin, ierr)
        call setfhv('USER2', permax, ierr)
        call setkhv('KUSER1', 'PER_MIN  ', ierr)
        call setkhv('KUSER2', 'PER_MAX  ', ierr)
        call setfhv('DELTA   ',dt,ierr)
        call setfhv('DEPMIN  ',depmin,ierr)
        call setfhv('DEPMAX  ',depmax,ierr)
        call setfhv('DEPMEN  ',depmen,ierr)
        call setfhv('TIMMAX  ',btime + indmax*dt,ierr)
        call setfhv('TIMMIN  ',btime + indmin*dt,ierr)
        call bwsac (1,MXPTS,name,seis)
        return
        end

        subroutine domxmn(x,npts,depmax,depmin)
        implicit none
c-----
c       get extremal values of the time series
c-----
        real*4 x(*)
        real*4 depmax,depmin
        integer*4 npts
        integer i
        real sum

        depmax = -1.0e+38
        depmin =  1.0e+38
        sum = 0.0
        do 1000 i=1, npts
            if( x(i) .gt. depmax) depmax = x(i)
            if( x(i) .lt. depmin) depmin = x(i)
            sum = sum + x(i)
 1000   continue
        return
        end

        subroutine npow2(nsamp,npts,npts21)
        implicit none
c-----
c       Given nsamp, find npts >= nsamp such that npts is a power of 2
c-----  
        integer*4 nsamp, npts, npts21
        npts = 1
 1000   continue
            npts = 2*npts
            if(npts.lt.nsamp)go to 1000
        npts21 = npts/2 + 1
        return
        end

        subroutine gcmdln(fsac,fdisp,doauto,itermx,lshift,phasecor)
        implicit none
        character fsac*(*), fdisp*(*)
        logical doauto
        logical lshift
        integer itermx
        real phasecor

        common/break/nbpts
        integer nbpts

        character names*256
        integer nmarg, i
        integer mnmarg 

        nmarg = mnmarg()
        i = 0
        fsac = ' '
        fdisp = ' '
        doauto = .false.
        lshift = .false.
        itermx = 3
        phasecor = 0.0
        nbpts=2
 1000   continue
            i = i + 1
            if(i.gt.nmarg)go to 9000
            call mgtarg(i,names)
            if(names(1:2).eq.'-F')then
                i = i + 1
                call mgtarg(i,names)
                fsac = names
            else if(names(1:2).eq.'-D')then
                i = i + 1
                call mgtarg(i,names)
                fdisp = names
            else if(names(1:2).eq.'-I')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,i10)')itermx
            else if(names(1:3).eq.'-NB')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,i10)')nbpts
            else if(names(1:2).eq.'-S')then
                lshift = .true.
            else if(names(1:4).eq.'-P4M')then
                phasecor = - 3.1415927/4.
            else if(names(1:4).eq.'-P4P')then
                phasecor =   3.1415927/4.
            else if(names(1:5).eq.'-AUTO')then
                doauto = .true.
            else if(names(1:2).eq.'-?')then
                call usage()
            else if(names(1:2).eq.'-h')then
                call usage()
            endif
        go to 1000
 9000   continue
        return
        end

        subroutine usage()
        implicit none
        integer LIN, LOT, LER
        parameter (LIN=5, LER=0, LOT=6)
           write(LER,*)
     :'Usage: sacmat96 -F sacfile_binary',
     :' -D Disp_file -I max_iter -AUTO -S -? -h'
        write(LER,*)
     :' -F sacfile_binary (default= required) SAC file'
        write(LER,*)
     :' -D Disp_file (default= required) SURF96 dispersion file'
        write(LER,*)
     :' -I  max_iter  (default 3)      maximum iterations'
        write(LER,*)
     :' -P4M         (default = false) -pi/4 correction'
        write(LER,*)
     :' -P4P         (default = false) +pi/4 correction'
        write(LER,*)
     :' -AUTO        (default = false) auto processing - no plot'
        write(LER,*)
     1' -S           (default false)   output compressed traces'
        write(LOT,*)
     1' -?           (default false)   usage'
        write(LOT,*)
     1' -h           (default false)   usage'
        stop
        end

        subroutine getper(wn,NX,pmin,pmax,nper) 
c-----
c       automatically create periods from the [pmin,pmax] limits
c       wn()    R*4 array of periods
c       NX  I   dimension of array
c       pmin    R*4 minimum period
c       pmax    R*4 maximum period
c       nper    I*4 number of periods generated
c----
        implicit none
        integer MAXFRQ
        parameter (MAXFRQ=100)
        integer NX, nper
        real wn(NX)
        real pmin, pmax

        real ymxlog, ymmin, tenpow, p
        integer iy, ii, nocy, mper, jj
        integer NP
        parameter (NP=40)
        real pfac(NP)
        data pfac/1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 
     1      2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,
     2      3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8,
     3      5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5/

        nper = 0
c-----
c       determine starting power
c-----
        ymxlog = alog10(pmax)   
        ymmin  = alog10(pmin)   
        nocy = ymxlog - ymmin + 1 
        iy = ymmin
        if(ymmin .lt. 0)iy = iy - 1
        do 100 ii=iy,iy+nocy+2
            tenpow = 10.0**ii
            do 200 jj=1,NP
                p = pfac(jj)*tenpow
                if(p .ge. 0.99*pmin .and. p .le. 1.01*pmax)then
                    nper = nper + 1
                    mper = nper
                    wn(nper) = p
                    if(nper.eq.NX)go to 1000
                endif
 200        continue
 100    continue
 1000   continue
        return
        end

        subroutine uniq(x,y,m)
c-----
c       The x,y arrays is in order of increasing period in the x array
c       go through the arrays to eliminate any duplicate periods
c       so that the interpolation does not blow up
c-----
        implicit none
        integer m
        real x(m), y(m)
        integer i, mout, j
c-----
c       do a gruesome search with array replacement
c-----
        mout = m
        i = 1
 1000   continue
            i = i + 1
            if(i.gt.mout)go to 2000

            if(x(i) .eq. x(i-1))then
c-----
c           shift left
c-----
                mout = mout - 1
                do 3000 j=i,mout
                    x(j) = x(j+1)
                    y(j) = y(j+1)
 3000           continue
                i = i - 1
            endif
        go to 1000
 2000   continue
        m = mout
        return
        end

        subroutine zfour(zarr,nn,isign,dt,df) 
c-----
c     THE input is a complex array
c     which has numbers stored in memory as
c     R1, I1, R2, I2, ..., Rnn, Inn
c     where nn must be a power of 2 R and I are the real and imaginary
c     parts of the complex number
c
c     For isign -1 this is a complex time series
c     For isign +1 this is a complex frequency series with
c        index 1 (in fortran corresponding to f=0
c              2                              f=df
c            nn/2 + 1                         f = 1/2dt = Nyquist
c            nn - 1                           f = -2df
c             nn                              f = -df

c-----
c     the cooley-tookey fast fourier transform in usasi basic fortran
c     transform(j) = sum(datc(i)*w**((i-1)(j-1)), where i and j run
c     from 1 to nn and w = exp(isign*2*pi*sqrt(-1)/nn).  datc is a one-
c     dimensional complex array (i.e., the real and imaginary parts of
c     datc are located immediately adjacent in storage, such as fortran
c     places them) whose length nn is a power of two.  isign
c     is +1 or -1, giving the sign of the transform.  transform values
c     are returned in array datc, replacing the input datc.  the time is
c     proportional to n*log2(n), rather than the usual n**2
c     rms resolution error being bounded by 6*sqrt(i)*log2(nn)*2**(-b),
c     b is the number of bits in the floating point fraction.
c
c     the program computes df from dt, dt from df and checks to see
c     if they are consistent. In addition, the transforms are multiplied
c     by dt or df to make the results dimensionally correct
c
c     This is a slightly modified version of the original Brenner routine
c     The changes were to add physical dimensions to the transform
c     and to make it all complex
c-----
        complex zarr(*) 
        integer nn, isign
        real dt, df

        complex ztemp
c-----
c       ensure that the dt and df are defined and
c       consistent
c-----
        if(dt.eq.0.0) dt = 1./(nn*df) 
        if(df.eq.0.0) df = 1./(nn*dt) 
        if(dt.ne.(nn*df)) df = 1./(nn*dt) 
c-----
c       now begin the transform
c-----
        jj = 1
        do 5 ii=1,nn 
        if(ii .lt. jj) then
              ztemp = zarr(jj)
              zarr(jj) = zarr(ii)
              zarr(ii) = ztemp
        endif
        mm = nn/2
    3   continue
        if(jj.le.mm) then
            go to 55
        else 
              jj = jj - mm
              mm = mm /2
              if(mm.lt.1)then
                  go to 55
              else
                  go to 3
              endif
        endif
   55   continue
        jj = jj + mm
    5   continue
        mmax = 1 
c-----
    6 continue
        if(mmax .lt. nn)then
            go to 7
        else if(mmax .ge. nn)then
            go to 10
        endif
    7   continue
        istep= 2*mmax 
        theta = 6.283185307/(isign*2.0*mmax) 
        sinth=sin(theta/2.) 
        wstpr=-2.*sinth*sinth 
        wstpi=sin(theta) 
        wr=1.0 
        wi=0.0 
        do 9 mm=1,mmax
              do 8 ii=mm,nn,istep
                    jj=ii+mmax
                    ztemp=cmplx(wr,wi)*zarr(jj)
                    zarr(jj) = zarr(ii) - ztemp
                    zarr(ii) = zarr(ii) + ztemp
    8         continue
c-----
c       use trig relations to compute the next sin/cos
c       without actually calling sin() or cos()
c-----
              tempr = wr 
              wr = wr*wstpr-wi*wstpi + wr 
              wi = wi*wstpr+tempr*wstpi + wi 
    9   continue
        mmax = istep 
        go to 6 
c-----
c       transform is done
c-----
   10   continue 
c-----
c     give the arrays the proper physical dimensions
c-----
        if(isign.lt.0)then
c-----
c             time to frequency domain
c-----
              do  ii = 1,nn
                    zarr(ii) = zarr(ii) * dt
              enddo
        else
c-----
c             frequency to time domain
c-----
              do ii = 1,nn
                    zarr(ii) = zarr(ii) * df
              enddo
        endif
        return
        end
