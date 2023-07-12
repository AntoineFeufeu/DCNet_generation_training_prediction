      subroutine dosh(com,cwvno,gg)
c-----
c     changes
c-----
c     27 DEC 2019 - the PDD PDS PSS PEX PVF PHF are computed for
c        receivers in a fluid
c-----      subroutine dosh(com,cwvno,gg)
      implicit none
      complex com
      complex cwvno
      complex gg(21)

        common/modlly/mmax
        integer mmax
        common/jbdry/jtop,jbot
        integer jtop, jbot

      common/bcsh/alpsh,betsh
      complex*16 alpsh(2,2), betsh(2,2)

        integer NL
        parameter (NL=200)
      common /elmat/el,elinv,es
      complex*16 el(NL,2,2)
      complex*16 elinv(NL,2,2)
      complex*16 es(NL)

        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d,a,b, rho,qa,qb,etap,etas,frefp,frefs
        common/depref/refdep
        real refdep

        common/verby/verbose, dout
        logical verbose, dout

      common/srcrec/isrc,irec
      integer isrc,irec

      common/downsh/RD,TD
      complex*16 RD(NL), TD(NL)

      common/upsh/RU,TU
      complex*16 RU(NL), TU(NL)

      common/srctype/lsrc
      character lsrc*2

      common/coeffinalsh/CDN,CUP
      complex*16 CUP(NL), CDN(NL)

      common/debug/mlam,mmu,anu,bnu
      complex*16 mlam(NL), mmu(NL)
      complex*16 anu(NL),bnu(NL)

        common/layertype/isfluid(NL),allsolid
        logical isfluid,allsolid

      complex*16 ut, tt
      complex*16 om, wvno
      real*8 rom, iom, rwvn, iwvn
      rom = real(com)
      iom = aimag(com)
      rwvn = real(cwvno)
      iwvn = aimag(cwvno)
      om = dcmplx(rom,iom)
      wvno = dcmplx(rwvn, iwvn)

c-----
c     for SH neither the source nor the receiver can be
c     in a fluid. All layers between the source and the
c     receiver must be solid
c-----
      if(.not. allsolid)then
          gg(13) = cmplx(0.0,0.0)
          gg(14) = cmplx(0.0,0.0)
          gg(15) = cmplx(0.0,0.0)
      else
c-----
c         evaluate the E Einv matrices, exponentials
c-----
          call doel(om,wvno)
c-----
c         define the boundary matrices
c-----
          call dobotsh(alpsh, jbot)
          call dotopsh(betsh, jtop)
          if(verbose) then
                WRITE(6,*)'MMAX:',mmax
                WRITE(6,*)'ALP:'
                WRITE(6,*) 1,1, alpsh(1,1)
                WRITE(6,*) 1,2, alpsh(1,2)
                WRITE(6,*) 2,1, alpsh(2,1)
                WRITE(6,*) 2,2, alpsh(2,2)
                WRITE(6,*)'BET:'
                WRITE(6,*) 1,1, betsh(1,1)
                WRITE(6,*) 1,2, betsh(1,2)
                WRITE(6,*) 2,1, betsh(2,1)
                WRITE(6,*) 2,2, betsh(2,2)
          endif
          call botupsh()
          call topdownsh()
c-----
c         note this would be faster if juse one call
c         to dosrcrec would return all gg s
c-----
          if(isrc.gt.0 .and. irec.gt.0)then
              lsrc='DS'
              call dosrcrecsh(wvno,ut,tt)
              gg(13) = (ut)
              lsrc='SS'
              call dosrcrecsh(wvno,ut,tt)
              gg(14) = (ut)
              lsrc='HF'
              call dosrcrecsh(wvno,ut,tt)
              gg(15) = (ut)
          endif
      endif
      return
      end

      subroutine dosrcrecsh(wvno,ut,tt)
c-----
c     test the source excitation
c-----
      implicit none
      complex*16 wvno

        common/modlly/mmax
        integer mmax

        integer NL
        parameter (NL=200)

        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d,a,b, rho,qa,qb,etap,etas,frefp,frefs
        common/depref/refdep
        real refdep

      common /elmat/el,elinv,es
      complex*16 el(NL,2,2)
      complex*16 elinv(NL,2,2)
      complex*16 es(NL)

      common/debug/mlam,mmu,anu,bnu
      complex*16 mlam(NL), mmu(NL)
      complex*16 anu(NL),bnu(NL)

      common/srcrec/isrc,irec
      integer isrc,irec

      common/result/gg
      complex*16 gg(21)

      common/downsh/RD,TD
      complex*16 RD(NL), TD(NL)

      common/upsh/RU,TU
      complex*16 RU(NL), TU(NL)

      common/srctype/lsrc
      character lsrc*2

        common/ctlrud/rpup,rpdn,rsup,rsdn,dorud
        logical rpup, rpdn, rsup, rsdn, dorud
        common/ctlsud/spup,spdn,ssup,ssdn,dosud
        logical spup, spdn, ssup, ssdn, dosud

        common/shcontrol/ltop,lbot
        integer ltop, lbot



      complex*16 sigmau, sigmad
      complex*16 CU, CD

      complex*16 IRuRdInv
      complex*16 IRdRuInv

      complex*16 ut, tt

      integer i
      real twopi

      complex*16 cdzero, cdone
      complex*16 cfac

      cdzero = dcmplx(0.0d+00,0.0d+00)
      cdone  = dcmplx(1.0d+00,0.0d+00)
      twopi = 6.2831853
c-----
c    safety
c-----
      if(isrc.lt.0 .or. isrc.gt.mmax)return
      if(irec.lt.0 .or. irec.gt.mmax)return
C      WRITE(6,*)'---------------------------------------'
C      WRITE(6,*)'i,TD(i),RD(i),TU(i),RU(i)'
C      do i=1,mmax
C        WRITE(6,*)i,TD(i),RD(i),TU(i),RU(i)
C      enddo


      if(lsrc.eq.'DS')then
c-----
c     DS
c-----
           sigmau =  elinv(isrc,1,1)*(-1/(twopi*mmu(isrc)))
           sigmad =  elinv(isrc,2,1)*(-1/(twopi*mmu(isrc)))
      else if(lsrc.eq.'DD')then
c-----
c     SS
c-----
           sigmau =  elinv(isrc,1,2)*( wvno/(twopi))
           sigmad =  elinv(isrc,2,2)*( wvno/(twopi))
      else if(lsrc.eq.'SS')then
c-----
c     SS
c-----
           sigmau =  elinv(isrc,1,2)*( wvno/(twopi))
           sigmad =  elinv(isrc,2,2)*( wvno/(twopi))
      else if(lsrc.eq.'HF')then
c-----
c     HF
c-----
           sigmau =  elinv(isrc,1,2)*( 1/(twopi))
           sigmad =  elinv(isrc,2,2)*( 1/(twopi))
      endif

c-----
c   if dosud is true then apply the up down filter
c   This is done by zeroing out terms
c-----
       if(dosud)then
           if(.not.ssdn)then
                 sigmad = cdzero
           endif
           if(.not.ssup)then
                 sigmau = cdzero
           endif
       endif

      if(irec.lt.isrc)then
c-----
c         get CD and CU for layer isrc-1
c         this is ugly if the receiver is at top of layer isrc-1
c-----
C         IRdRuInv = cdone/
C    1        (cdone-es(isrc)*RD(isrc)*es(isrc-1)*RU(isrc-1))
          cfac = cdone-es(isrc)*RD(isrc)*es(isrc-1)*RU(isrc-1)
          if(cdabs(cfac) .lt.1.0d-09)then
               IRdRuInv = cdzero
          else
               IRdRuInv = cdone/cfac
          endif
          IRdRuInv = cdone/
     1        (cdone-es(isrc)*RD(isrc)*es(isrc-1)*RU(isrc-1))
          CU= IRdRuInv *(es(isrc)*RD(isrc) * sigmad - sigmau)
c-----
c         not this CD is not used
c-----
          CD = RU(isrc-1)*CU
          do i=isrc-2, irec, -1
c-----
c                get the CU and CD now for layer i-1
                 CU = TU(i+1)*CU
                 CD = RU(i)*CU
          enddo
      else if(irec.ge.isrc)then
C         IRuRdInv = cdone/
C    1          (cdone-es(isrc-1)*RU(isrc-1)*es(isrc)*RD(isrc))
          cfac = cdone-es(isrc-1)*RU(isrc-1)*es(isrc)*RD(isrc)
          if(cdabs(cfac) .lt.1.0d-07)then
               IRuRdInv = cdzero
          else
               IRuRdInv = cdone/cfac
          endif
          IRuRdInv = cdone/
     1          (cdone-es(isrc-1)*RU(isrc-1)*es(isrc)*RD(isrc))
          CD = IRuRdInv * (sigmad - es(isrc-1)*RU(isrc-1)*sigmau)
          CU = RD(isrc)*CD
c-----
c         not this CU is not used
c-----
          do i=isrc+1,irec
             CD = TD(i-1)*CD
             CU = RD(i)*CD
          enddo
      endif
c-----
c   if dorud is true then apply the up down filter
c   This is done by zeroing out terms
c-----
       if(dorud)then
           if(.not.rsdn)then
                 CD = cdzero
           endif
C          if(.not.rpdn)then
C                CD(1) = cdzero
C          endif
           if(.not.rsup)then
                 CU = cdzero
           endif
C          if(.not.rpup)then
C                CU(1) = cdzero
C          endif
       endif

c-----
c     get motion at TOP of the layer
c----- 
c     ut = el(irec,1,1)*cu + el(irec,1,2)*es(irec)*cd
c     tt = el(irec,2,1)*cu + el(irec,2,2)*es(irec)*cd
c     WRITE(6,*)'BOT: isrc,irec,ut,tt:',isrc,irec,ut,tt
      ut = el(irec,1,1)*es(irec)*cu + el(irec,1,2)*1*cd
      tt = el(irec,2,1)*es(irec)*cu + el(irec,2,2)*1*cd
c     WRITE(6,*)'TOP: isrc,irec,ut,tt:',isrc,irec,ut,tt
      return
      end

      subroutine dobotsh(alpsh, jbot)
c-----
c     alp = G BN-1
c-----
      implicit none
      complex*16 alpsh(2,2)
      integer jbot

        integer NL
        parameter (NL=200)
      common /elmat/el,elinv,es
      complex*16 el(NL,2,2)
      complex*16 elinv(NL,2,2)
      complex*16 es(NL)

        common/modlly/mmax
        integer mmax
        common/shcontrol/ltop,lbot
        integer ltop, lbot


 
      complex*16 cdzero, cdone

      cdzero = dcmplx(0.0d+00,0.0d+00)
      cdone  = dcmplx(1.0d+00,0.0d+00)

    
c-----
c     to correcly handle the rigid free
c     assume that mmax is the layer just
c     above the boundary
c       alp = G E
c-----

      if(jbot.eq.1) then
c-----
c        bottom free
c-----
         alpsh(1,1) = el(mmax,2,1)
         alpsh(1,2) = el(mmax,2,2)
         alpsh(2,1) = el(mmax,1,1)
         alpsh(2,2) = el(mmax,1,2)
      else if(jbot.eq.2)then
c-----
c        bottom rigid
c-----
         alpsh(1,1) = el(mmax,1,1)
         alpsh(1,2) = el(mmax,1,2)
         alpsh(2,1) = el(mmax,2,1)
         alpsh(2,2) = el(mmax,2,2)
      else if(jbot.eq.3)then
c-----
c        bottom halfspace3
c        there is no bottom layer mmax+1 - so assume that it has the same 
c        properties as the mmax layer
c        so alp = G E = I
c-----
         alpsh(1,1) = cdone
         alpsh(1,2) = cdzero
         alpsh(2,1) = cdzero
         alpsh(2,2) = cdone
      endif
c-----
c     special case if this is a solid bounded by fluid layers
c     below. This will then force a free surface condition there
c-----
      if(lbot.lt.mmax)then
         alpsh(1,1) = el(lbot,2,1)
         alpsh(1,2) = el(lbot,2,2)
         alpsh(2,1) = el(lbot,1,1)
         alpsh(2,2) = el(lbot,1,2)
      endif
      return
      end

      subroutine dotopsh(betsh, jtop)
c-----
c     bet = H(inv) E1
c-----
      implicit none
      complex*16 betsh(2,2)
      integer jtop

        integer NL
        parameter (NL=200)
      common /elmat/el,elinv,es
      complex*16 el(NL,2,2)
      complex*16 elinv(NL,2,2)
      complex*16 es(NL)

        common/modlly/mmax
        integer mmax

        common/shcontrol/ltop,lbot
        integer ltop, lbot


      complex*16 cdzero, cdone

      cdzero = dcmplx(0.0d+00,0.0d+00)
      cdone  = dcmplx(1.0d+00,0.0d+00)

      if(jtop.eq.1) then
c-----
c        top free
c-----
         betsh(1,1) = el(1,1,1)
         betsh(1,2) = el(1,1,2)
         betsh(2,1) = el(1,2,1)
         betsh(2,2) = el(1,2,2)
      else if(jtop.eq.2)then
c-----
c        top rigid
c-----
         betsh(1,1) = el(1,2,1)
         betsh(1,2) = el(1,2,2)
         betsh(2,1) = el(1,1,1)
         betsh(2,2) = el(1,1,2)
      else if(jtop.eq.3)then
c-----
c        top halfspace - inconsistency here since there is no 0 layer
c        parameter just assume that properties are same
c       bet = H inv E1   but H = E, thus bet =I
c-----
         betsh(1,1) = cdone
         betsh(1,2) = cdzero
         betsh(2,1) = cdzero
         betsh(2,2) = cdone
      endif
c-----
c     special case if the top of the solid structure is a fluid, then
c     this will then force a free surface condition there
c-----
c     there is a free surface there
c-----
      if(jtop.gt.1)then
         betsh(1,1) = el(jtop,1,1)
         betsh(1,2) = el(jtop,1,2)
         betsh(2,1) = el(jtop,2,1)
         betsh(2,2) = el(jtop,2,2)
      endif
      return
      end

      subroutine doel(om,wvno)
      implicit none
      complex*16 om
      complex*16 wvno

        common/modlly/mmax
        integer mmax
        integer NL
        parameter (NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d,a,b, rho,qa,qb,etap,etas,frefp,frefs
        common/depref/refdep
        real refdep

      common /elmat/el,elinv,es
      complex*16 el(NL,2,2)
      complex*16 elinv(NL,2,2)
      complex*16 es(NL)


      common/debug/mlam,mmu,anu,bnu
      complex*16 mlam(NL), mmu(NL)
      complex*16 anu(NL),bnu(NL)

        common/damp/alpha,ieqex
        real alpha
        integer ieqex

      integer i
      complex*16 ka, kb, atna, atnb
      integer iwats

        common/layertype/isfluid(NL),allsolid
        logical isfluid,allsolid

      complex*16 cdzero, cdone

      cdzero = dcmplx(0.0d+00,0.0d+00)
      cdone  = dcmplx(1.0d+00,0.0d+00)


      do i=1,mmax
           if(isfluid(i))then
              el(i,1,1)    = cdzero
              el(i,1,2)    = cdzero
              el(i,2,1)    = cdzero
              el(i,2,2)    = cdzero
              elinv(i,1,1) = cdzero
              elinv(i,2,1) = cdzero
              elinv(i,1,2) = cdzero
              elinv(i,2,2) = cdzero
              es(i)        = cdzero
           else
            call aten(om,qa(i),qb(i),ka ,kb ,alpha,
     1          a(i),b(i),atna,atnb,iwats,
     2          frefp(i),frefs(i)) 

              bnu(i) = cdsqrt(wvno*wvno -kb*kb)
              mmu(i) =  rho(i)*b(i)*b(i)*atnb*atnb
              es(i) = cdexp(-bnu(i)*d(i))
              el(i,1,1)    = wvno
              el(i,1,2)    = wvno
              el(i,2,1)    =   mmu(i)*wvno*bnu(i)
              el(i,2,2)    = - mmu(i)*wvno*bnu(i)
              elinv(i,1,1) = mmu(i)*wvno*bnu(i) 
     1                        /( 2.*mmu(i)*wvno*wvno*bnu(i))
              elinv(i,2,1) = mmu(i)*wvno*bnu(i) 
     1                        /( 2.*mmu(i)*wvno*wvno*bnu(i))
              elinv(i,1,2) =   wvno / ( 2.*mmu(i)*wvno*wvno*bnu(i))
              elinv(i,2,2) = - wvno / ( 2.*mmu(i)*wvno*wvno*bnu(i))
           endif
          
      enddo
      return
      end

      subroutine botupsh()
      implicit none
        common/modlly/mmax
        integer mmax
      common/bcsh/alpsh,betsh
      complex*16 alpsh(2,2), betsh(2,2)
        integer NL
        parameter (NL=200)
      common /elmat/el,elinv,es
      complex*16 el(NL,2,2)
      complex*16 elinv(NL,2,2)
      complex*16 es(NL)

      common/downsh/RD,TD
      complex*16 RD(NL), TD(NL)

        common/verby/verbose, dout
        logical verbose, dout

        common/layertype/isfluid(NL),allsolid
        logical isfluid,allsolid

        common/shcontrol/ltop,lbot
        integer ltop, lbot

      complex*16 ee(2,2)
      integer i
      complex*16 cdzero, cdone
      complex*16 cfac

      cdzero = dcmplx(0.0d+00,0.0d+00)
      cdone  = dcmplx(1.0d+00,0.0d+00)
c-----
c     for a completely solid system define the
c     RD(mmax) amd then proceed upward recursively to
c     define RD(1) and TD(1)
c     
c     When there are mixed fluid and solid layers
c     we only compute non-zero values for the solid region
c     in which both the source and receiver lie. Thus
c     for irec > isrc, the completely solid problem would 
c     define the RD and TD at the postions
c        1    ... isrc ... irec ... mmax    for irec > isrc
c     and for the solid bounded by a fluid
c        ltop ... isrc ... irec ... lbot
c     When irec < isrc we would have
c        1    ... irec ... isrc ... mmax    for irec > isrc
c     and for the solid bounded by a fluid
c        ltop ... irec ... isrc ... lbot
c     Note that in the computations here isrc and irec 
c     are not used. We worry about the ltop and lbot
c-----
c     get the R at the bottom from
c      R = - a11(inv) a12 e(- nub d N-1)
c-----
        if(lbot.eq.mmax)then
           RD(mmax) = - alpsh(1,2) * es(mmax  ) /alpsh(1,1)
           TD(mmax) = cdzero
        else
           do i=mmax,lbot+1,-1
             RD(i) =  cdzero
             TD(i) = cdzero
           enddo
           RD(lbot) = - alpsh(1,2) * es(lbot  ) /alpsh(1,1)
           TD(lbot) = cdzero
        endif

c-----
c     march up through the layers
c     normally this is traight forward. However
c     we must work if the layer is a fluid. We
c    
c-----
      do i = min(mmax-1,lbot-1),ltop,-1
        ee(1,1) = elinv(i,1,1)*el(i+1,1,1) + elinv(i,1,2)*el(i+1,2,1)
        ee(1,2) = elinv(i,1,1)*el(i+1,1,2) + elinv(i,1,2)*el(i+1,2,2)
        ee(2,1) = elinv(i,2,1)*el(i+1,1,1) + elinv(i,2,2)*el(i+1,2,1)
        ee(2,2) = elinv(i,2,1)*el(i+1,1,2) + elinv(i,2,2)*el(i+1,2,2)
        cfac = (ee(2,1)*es(i+1)*RD(i+1) + ee(2,2) ) 
        if(cdabs(cfac).lt.1.0d-09)then
           TD(i) = cdzero
        else
           TD(i) = (cdone/cfac)*es(i)
        endif
        TD(i) = ( 1/(ee(2,1)*es(i+1)*RD(i+1) + ee(2,2) ) ) * es(i)
        RD(i) = (ee(1,1)*RD(i+1) * es(i+1) + ee(1,2))*TD(i)
        
C        if(verbose)then
C              WRITE(6,*)'EE11:',ee(1,1)
C              WRITE(6,*)'EE12:',ee(1,2)
C              WRITE(6,*)'EE21:',ee(2,1)
C              WRITE(6,*)'EE22:',ee(2,2)
C              WRITE(6,*)'RD(',i,')=',RD(i),' TD(',i,')=',TD(i)
C              WRITE(6,*)'es(',i,')=',es(i)
C        endif
      enddo
      if(ltop.gt.1)then
          do i=ltop-1,1
             TD(i) = cdzero
             RD(i) = cdzero
          enddo
      endif
      return
      end

      subroutine topdownsh()
      implicit none
        common/modlly/mmax
        integer mmax
      common/bcsh/alpsh,betsh
      complex*16 alpsh(2,2), betsh(2,2)
        integer NL
        parameter (NL=200)
      common /elmat/el,elinv,es
      complex*16 el(NL,2,2)
      complex*16 elinv(NL,2,2)
      complex*16 es(NL)

      common/upsh/RU,TU
      complex*16 RU(NL), TU(NL)

        common/verby/verbose, dout
        logical verbose, dout

        common/layertype/isfluid(NL),allsolid
        logical isfluid,allsolid
 
        common/shcontrol/ltop,lbot
        integer ltop, lbot

      complex*16 ff(2,2)
      integer i
      complex*16 cdzero, cdone
      complex*16 cfac

      cdzero = dcmplx(0.0d+00,0.0d+00)
      cdone  = dcmplx(1.0d+00,0.0d+00)
c-----
c     for a completely solid system define the
c     RU(1) amd then proceed downward recursively to
c     define RU(mmax) and TU(mmax)
c     
c     When there are mixed fluid and solid layers
c     we only comp;ute non-zero values for the solid region
c     in which both the source and receiver lie. Thus
c     for irec > isrc, the completely solid problem would 
c     define the RU and TU at the postions
c        1    ... isrc ... irec ... mmax    for irec > isrc
c     and for the solid bounded by a fluid
c        ltop ... isrc ... irec ... lbot
c     When irec < isrc we would have
c        1    ... irec ... isrc ... mmax    for irec > isrc
c     and for the solid bounded by a fluid
c        ltop ... irec ... isrc ... lbot
c     Note that in the computations here isrc and irec 
c     are not used. We worry about the ltop and lbot
c-----
c     get the R at the top from
c      R = - b22(inv) b21 e(- nub d N-1)
c-----
        if(ltop.eq.1)then
             RU(1) = - betsh(2,1) * es(1  ) /betsh(2,2)
             TU(1) = cdzero
        else
             do i=1,ltop-1
                  RU(i) = cdzero
                  TU(i) = cdzero
             enddo
             RU(ltop) = - betsh(2,1) * es(ltop  ) /betsh(2,2)
             TU(ltop) = cdzero
        endif

c-----
c     march up through the layers
c-----
      do i = ltop,min(mmax-1,lbot-1)
        ff(1,1) = elinv(i+1,1,1)*el(i,1,1) + elinv(i+1,1,2)*el(i,2,1)
        ff(1,2) = elinv(i+1,1,1)*el(i,1,2) + elinv(i+1,1,2)*el(i,2,2)
        ff(2,1) = elinv(i+1,2,1)*el(i,1,1) + elinv(i+1,2,2)*el(i,2,1)
        ff(2,2) = elinv(i+1,2,1)*el(i,1,2) + elinv(i+1,2,2)*el(i,2,2)
        cfac = ff(1,1) + ff(1,2)*es(i)*RU(i)
        if(cdabs(cfac).lt.1.0d-07)then
          TU(i+1) = cdzero
        else
          TU(i+1) = (cdone/cfac)*es(i+1)
        endif
        TU(i+1) = ( 1/(ff(1,1) + ff(1,2)*es(i)*RU(i) ) ) * es(i+1)
        RU(i+1) =     (ff(2,1) + ff(2,2)*es(i)*RU(i) )*TU(i+1)
      enddo
      if(lbot.ne.mmax)then
         do i=lbot,mmax-1
            RU(i+1) = cdzero
            RU(i+1) = cdzero
         enddo
       endif
      return
      end

      subroutine dopsv(com,cwvno,gg)
      implicit none
      complex com
      complex cwvno
      complex gg(21)

        common/modlly/mmax
        integer mmax

        common/jbdry/jtop,jbot
        integer jtop, jbot

      common/bcpsv/alppsv,betpsv
      complex*16 alppsv(4,4), betpsv(4,4)

        integer NL
        parameter (NL=200)
      common /ermat/er,erinv,ep,es,elog
      complex*16 er(NL,4,4)
      complex*16 erinv(NL,4,4)
      complex*16 ep(NL),es(NL)
      real elog(NL)

      common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
      real d,a,b, rho,qa,qb,etap,etas,frefp,frefs
      common/depref/refdep
      real refdep
        common/layertype/isfluid(NL),allsolid
        logical isfluid,allsolid

      common/verby/verbose, dout
      logical verbose, dout

      common/srcrec/isrc,irec
      integer isrc,irec

      common/downpsv/RD, TD
      complex*16 RD(NL,2,2), TD(NL,2,2)

      common/uppsv/RU, TU
      complex*16 RU(NL,2,2), TU(NL,2,2)

      common/srctype/lsrc
      character lsrc*2

      common/coeffinalpsv/CDN,CUP
      complex*16 CDN(NL,2), CUP(NL,2)

      common/debug/mlam,mmu,anu,bnu
      complex*16 mlam(NL), mmu(NL)
      complex*16 anu(NL),bnu(NL)


      complex*16 ur,uz,tz,tr
      complex*16 om, wvno
      real*8 rom, iom, rwvn, iwvn
      integer iwatr
      integer j

      rom = real(com)
      iom = aimag(com)
      rwvn = real(cwvno)
      iwvn = aimag(cwvno)
      om = dcmplx(rom,iom)
      wvno = dcmplx(rwvn, iwvn)
      if(b(irec).eq.0.0)then
          iwatr = 1
      else
          iwatr = 0
      endif

c     WRITE(6,*)'RBH dopsv om, wvno',om,wvno
c-----
c     evaluate the E Einv matrices, exponentials
c-----
      call doer(om,wvno)
c-----
c     define the boundary matrices
c-----
      call dobotpsv( alppsv, jbot)
      call dotoppsv( betpsv, jtop)
      call botuppsv()
      call topdownpsv()
c-----
c     note this would be faster if juse one call
c     to dosrcrec would return all gg s
c-----
      if(isrc.gt.0 .and. irec.gt.0)then
            lsrc='DD'
            call dosrcrecpsv(wvno,ur,uz,tz,tr)
            gg(1)  = (uz)
            gg(2)  = (ur)
            if(iwatr.eq.1)then
                gg(17) = (tz)/(om*om)
            endif
            lsrc='DS'
            call dosrcrecpsv(wvno,ur,uz,tz,tr)
            gg(3)  = (uz)
            gg(4)  = (ur)
            if(iwatr.eq.1)then
                gg(18) = (tz)/(om*om)
            endif
            lsrc='SS'
            call dosrcrecpsv(wvno,ur,uz,tz,tr)
            gg(5)  = (uz)
            gg(6)  = (ur)
            if(iwatr.eq.1)then
                gg(19) = (tz)/(om*om)
            endif
            lsrc='EX'
            call dosrcrecpsv(wvno,ur,uz,tz,tr)
            gg(7)  = (uz)
            gg(8)  = (ur)
            if(iwatr.eq.1)then
                gg(16) = (tz)/(om*om)
            endif
            lsrc='VF'
            call dosrcrecpsv(wvno,ur,uz,tz,tr)
            gg(9)  = (uz)
            gg(10) = (ur)
            if(iwatr.eq.1)then
                gg(20) = (tz)/(om*om)
            endif
            lsrc='HF'
            call dosrcrecpsv(wvno,ur,uz,tz,tr)
            gg(11) = (uz)
            gg(12) = (ur)
            if(iwatr.eq.1)then
                gg(21) = (tz)/(om*om)
            endif
      endif
c-----
c     everything up to here agrees with theory
c     that Uz is positive downward. To agree with other
c     CPS codes, e.g., the Uz is positive up, we flig the sign
c     of the vertical components of gg
c-----
c-----
c       invert the vertical
c-----
       do j=1,12,2
           gg(j) = -gg(j)
       enddo

      return
      end

      subroutine dosrcrecpsv(wvno,ur,uz,tz,tr)
c-----
c     test the source excitation
c-----
      implicit none
      complex*16 wvno
      complex*16 ur, uz, tz, tr

        common/modlly/mmax
        integer mmax

        common/verby/verbose, dout
        logical verbose, dout

        integer NL
        parameter (NL=200)

        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d,a,b, rho,qa,qb,etap,etas,frefp,frefs
        common/depref/refdep
        real refdep
        common/layertype/isfluid(NL),allsolid
        logical isfluid,allsolid

      common /ermat/er,erinv,ep,es,elog
      complex*16 er(NL,4,4)
      complex*16 erinv(NL,4,4)
      complex*16 ep(NL),es(NL)
      real elog(NL)

      common/debug/mlam,mmu,anu,bnu
      complex*16 mlam(NL), mmu(NL)
      complex*16 anu(NL),bnu(NL)

      common/srcrec/isrc,irec
      integer isrc,irec

      common/result/gg
      complex*16 gg(21)

      common/downpsv/RD, TD
      complex*16 RD(NL,2,2), TD(NL,2,2)

      common/uppsv/RU, TU
      complex*16 RU(NL,2,2), TU(NL,2,2)


      common/srctype/lsrc
      character lsrc*2

        common/ctlrud/rpup,rpdn,rsup,rsdn,dorud
        logical rpup, rpdn, rsup, rsdn, dorud
        common/ctlsud/spup,spdn,ssup,ssdn,dosud
        logical spup, spdn, ssup, ssdn, dosud


      complex*16 sigmau(2), sigmad(2)
      complex*16 CU(2), CD(2)

      complex*16 IRuRdInv(2,2)
      complex*16 IRdRuInv(2,2)
      complex*16 cmat(2,2)

      complex*16 s(4)
      complex*16 e(2)

      complex*16 ctmp1, ctmp2, cdet

      integer i
      real twopiinv
      real l2m

      complex*16 cdzero, cdone
      complex*16 cfac

      cdzero = dcmplx(0.0d+00,0.0d+00)
      cdone  = dcmplx(1.0d+00,0.0d+00)
      twopiinv = 1./(2.*3.1415927)
c     WRITE(6,*)'RBH srcrecpsv wvno isrc irec:',wvno,isrc,irec
c-----
c    safety
c-----
       if(isrc.lt.0 .or. isrc.gt.mmax)return
       if(irec.lt.0 .or. irec.gt.mmax)return
c      WRITE(6,*)'---------------------------------------'
c      WRITE(6,*)'k,TD(),RD( ),TU( ),RU( )'
c      do k=1,mmax
c        do i=1,2
c        do j=1,2
c        WRITE(6,*)k,i,j,TD(k,i,j),RD(k,i,j),TU(k,i,j),RU(k,i,j)
c        enddo
c        enddo
c      enddo


       l2m = rho(isrc)*a(isrc)*a(isrc)
       if(lsrc.eq.'DD')then
c-----
c     DD
c-----
           s(1) = cdzero
           s(2) = 2.*twopiinv/l2m
           s(3) = cdzero
           s(4) = wvno*twopiinv*((2.*b(isrc)/a(isrc))**2 -3.)
       else if(lsrc.eq.'DS')then
c-----
c     DS
c-----
           s(1) = twopiinv/mmu(isrc)
           s(2) = cdzero
           s(3) = cdzero
           s(4) = cdzero
       else if(lsrc.eq.'SS')then
c-----
c     SS
c-----
           s(1) = cdzero
           s(2) = cdzero
           s(3) = cdzero
           s(4) = - wvno*twopiinv
       else if(lsrc.eq.'EX')then
c-----
c     EX
c-----
           s(1) = cdzero
           s(2) = twopiinv/l2m
           s(3) = cdzero
           s(4) = 2.*wvno*twopiinv*(b(isrc)/a(isrc))**2
       else if(lsrc.eq.'VF')then
c-----
c     VF
c-----
           s(1) = cdzero
           s(2) = cdzero
           s(3) = - twopiinv
           s(4) = cdzero
       else if(lsrc.eq.'HF')then
c-----
c     HF
c-----
           s(1) = cdzero
           s(2) = cdzero
           s(3) = cdzero
           s(4) = - twopiinv
       endif
       sigmau(1) = erinv(isrc,1,1)*s(1) + erinv(isrc,1,2)*s(2) 
     2           + erinv(isrc,1,3)*s(3) + erinv(isrc,1,4)*s(4)
       sigmau(2) = erinv(isrc,2,1)*s(1) + erinv(isrc,2,2)*s(2) 
     2           + erinv(isrc,2,3)*s(3) + erinv(isrc,2,4)*s(4)
       sigmad(1) = erinv(isrc,3,1)*s(1) + erinv(isrc,3,2)*s(2) 
     2           + erinv(isrc,3,3)*s(3) + erinv(isrc,3,4)*s(4)
       sigmad(2) = erinv(isrc,4,1)*s(1) + erinv(isrc,4,2)*s(2) 
     2           + erinv(isrc,4,3)*s(3) + erinv(isrc,4,4)*s(4)

c-----
c   if dosud is true then apply the up down filter
c   This is done by zeroing out terms
c-----
       if(dosud)then
           if(.not.ssdn)then
                 sigmad(2) = cdzero
           endif
           if(.not.spdn)then
                 sigmad(1) = cdzero
           endif
           if(.not.ssup)then
                 sigmau(2) = cdzero
           endif
           if(.not.spup)then
                 sigmau(1) = cdzero
           endif
       endif



       if(irec.lt.isrc)then
          cmat(1,1) = cdone
     1       - ep(isrc)*RD(isrc,1,1)*ep(isrc-1)*RU(isrc-1,1,1)
     2       - ep(isrc)*RD(isrc,1,2)*es(isrc-1)*RU(isrc-1,2,1)
          cmat(1,2) = cdzero
     1       - ep(isrc)*RD(isrc,1,1)*ep(isrc-1)*RU(isrc-1,1,2)
     2       - ep(isrc)*RD(isrc,1,2)*es(isrc-1)*RU(isrc-1,2,2)
          cmat(2,1) = cdzero
     1       - es(isrc)*RD(isrc,2,1)*ep(isrc-1)*RU(isrc-1,1,1)
     2       - es(isrc)*RD(isrc,2,2)*es(isrc-1)*RU(isrc-1,2,1)
          cmat(2,2) = cdone
     1       - es(isrc)*RD(isrc,2,1)*ep(isrc-1)*RU(isrc-1,1,2)
     2       - es(isrc)*RD(isrc,2,2)*es(isrc-1)*RU(isrc-1,2,2)
c         cdet = cmat(1,1)*cmat(2,2) - cmat(1,2)*cmat(2,1)
c         IRdRuInv(1,1) =   cmat(2,2)/cdet
c         IRdRuInv(1,2) = - cmat(1,2)/cdet
c         IRdRuInv(2,1) = - cmat(2,1)/cdet
c         IRdRuInv(2,2) =   cmat(1,1)/cdet
          cdet = cmat(1,1)*cmat(2,2) - cmat(1,2)*cmat(2,1)
          if(cdabs(cdet).lt.1.0d-9)then
              cfac = cdzero
          else
              cfac = cdone/cdet
          endif
          IRdRuInv(1,1) =   cmat(2,2)*cfac
          IRdRuInv(1,2) = - cmat(1,2)*cfac
          IRdRuInv(2,1) = - cmat(2,1)*cfac
          IRdRuInv(2,2) =   cmat(1,1)*cfac
C          IRdRuInv = cdone/
C     1        (cdone-es(isrc)*RD(isrc)*es(isrc-1)*RU(isrc-1))
C          CU= IRdRuInv *(es(isrc)*RD(isrc) * sigmad - sigmau)
          e(1) = ep(isrc)*RD(isrc,1,1)*sigmad(1)
     1         + ep(isrc)*RD(isrc,1,2)*sigmad(2)
     2         - sigmau(1)
          e(2) = es(isrc)*RD(isrc,2,1)*sigmad(1)
     1         + es(isrc)*RD(isrc,2,2)*sigmad(2)
     2         - sigmau(2)
          CU(1) = IRdRuInv(1,1)*e(1) + IRdRuInv(1,2)*e(2)
          CU(2) = IRdRuInv(2,1)*e(1) + IRdRuInv(2,2)*e(2)
         
          CD(1) = RU(isrc-1,1,1)*CU(1) + RU(isrc-1,1,2)*CU(2)
          CD(2) = RU(isrc-1,2,1)*CU(1) + RU(isrc-1,2,2)*CU(2)
          do i=isrc-2, irec, -1
c-----
c                get the CU and CD now for layer i-1
c-----
             ctmp1 = TU(i+1,1,1)*CU(1) + TU(i+1,1,2)*CU(2)
             ctmp2 = TU(i+1,2,1)*CU(1) + TU(i+1,2,2)*CU(2)
             CU(1) = ctmp1
             CU(2) = ctmp2
             CD(1) = RU(i,1,1)*CU(1) + RU(i,1,2)*CU(2)
             CD(2) = RU(i,2,1)*CU(1) + RU(i,2,2)*CU(2)
          enddo
      else if(irec.ge.isrc)then
          cmat(1,1) = cdone
     1      - ep(isrc-1)*RU(isrc-1,1,1)*ep(isrc)*RD(isrc,1,1)
     2      - ep(isrc-1)*RU(isrc-1,1,2)*es(isrc)*RD(isrc,2,1)
          cmat(1,2) = cdzero
     1      - ep(isrc-1)*RU(isrc-1,1,1)*ep(isrc)*RD(isrc,1,2)
     2      - ep(isrc-1)*RU(isrc-1,1,2)*es(isrc)*RD(isrc,2,2)
          cmat(2,1) = cdzero
     1      - es(isrc-1)*RU(isrc-1,2,1)*ep(isrc)*RD(isrc,1,1)
     2      - es(isrc-1)*RU(isrc-1,2,2)*es(isrc)*RD(isrc,2,1)
          cmat(2,2) = cdone
     1      - es(isrc-1)*RU(isrc-1,2,1)*ep(isrc)*RD(isrc,1,2)
     2      - es(isrc-1)*RU(isrc-1,2,2)*es(isrc)*RD(isrc,2,2)
C         cdet = cmat(1,1)*cmat(2,2) - cmat(1,2)*cmat(2,1)
C         IRuRdInv(1,1) =   cmat(2,2)/cdet
C         IRuRdInv(1,2) = - cmat(1,2)/cdet
C         IRuRdInv(2,1) = - cmat(2,1)/cdet
C         IRuRdInv(2,2) =   cmat(1,1)/cdet
          cdet = cmat(1,1)*cmat(2,2) - cmat(1,2)*cmat(2,1)
          if(cdabs(cdet).lt.1.0d-9)then
              cfac = cdzero
          else
              cfac = cdone/cdet
          endif
          IRuRdInv(1,1) =   cmat(2,2)*cfac
          IRuRdInv(1,2) = - cmat(1,2)*cfac
          IRuRdInv(2,1) = - cmat(2,1)*cfac
          IRuRdInv(2,2) =   cmat(1,1)*cfac
C          IRuRdInv = cdone/
C     1		(cdone-es(isrc-1)*RU(isrc-1)*es(isrc)*RD(isrc))
C          CD = IRuRdInv * (sigmad - es(isrc-1)*RU(isrc-1)*sigmau)
         e(1) = - ep(isrc-1)*RU(isrc-1,1,1)*sigmau(1)
     1          - ep(isrc-1)*RU(isrc-1,1,2)*sigmau(2)
     2          + sigmad(1)
         e(2) = - es(isrc-1)*RU(isrc-1,2,1)*sigmau(1)
     1          - es(isrc-1)*RU(isrc-1,2,2)*sigmau(2)
     2          + sigmad(2)
         CD(1) = IRuRdInv(1,1)*e(1) + IRuRdInv(1,2)*e(2)
         CD(2) = IRuRdInv(2,1)*e(1) + IRuRdInv(2,2)*e(2)
         CU(1) = RD(isrc,1,1)*CD(1) + RD(isrc,1,2)*CD(2)
         CU(2) = RD(isrc,2,1)*CD(1) + RD(isrc,2,2)*CD(2)
         do i=isrc+1,irec
             ctmp1 = TD(i-1,1,1)*CD(1) + TD(i-1,1,2)*CD(2)
             ctmp2 = TD(i-1,2,1)*CD(1) + TD(i-1,2,2)*CD(2)
             CD(1) = ctmp1
             CD(2) = ctmp2
             CU(1) = RD(i,1,1)*CD(1) + RD(i,1,2)*CD(2)
             CU(2) = RD(i,2,1)*CD(1) + RD(i,2,2)*CD(2)
         enddo
      endif
c-----
c   if dorud is true then apply the up down filter
c   This is done by zeroing out terms
c-----
       if(dorud)then
           if(.not.rsdn)then
                 CD(2) = cdzero
           endif
           if(.not.rpdn)then
                 CD(1) = cdzero
           endif
           if(.not.rsup)then
                 CU(2) = cdzero
           endif
           if(.not.rpup)then
                 CU(1) = cdzero
           endif
       endif
       if(verbose)then
         WRITE(6,*)'SIGMAU:',sigmau
         WRITE(6,*)'SIGMAD:',sigmad
         WRITE(6,*)'CU    :',CU
         WRITE(6,*)'CD    :',CD
         WRITE(6,*)'S:     ',S
       endif
c-----
c     get motion
c----- 
          ur = er(irec,1,1)*ep(irec)*cu(1)+ er(irec,1,2)*es(irec)*cu(2)
     1       + er(irec,1,3)*cd(1)      + er(irec,1,4)*cd(2)
          uz = er(irec,2,1)*ep(irec)*cu(1)+ er(irec,2,2)*es(irec)*cu(2)
     1       + er(irec,2,3)*cd(1)      + er(irec,2,4)*cd(2)
          tz = er(irec,3,1)*ep(irec)*cu(1)+ er(irec,3,2)*es(irec)*cu(2)
     1       + er(irec,3,3)*cd(1)      + er(irec,3,4)*cd(2)
          tr = er(irec,4,1)*ep(irec)*cu(1)+ er(irec,4,2)*es(irec)*cu(2)
     1       + er(irec,4,3)*cd(1)      + er(irec,4,4)*cd(2)
C     WRITE(6,*)'TOP: isrc,irec,ur,uz,tz,tr:',isrc,irec,ur,uz,tz,tr
C         ur = er(irec,1,1)      *cu(1)+ er(irec,1,2)      *cu(2)
C    1       + er(irec,1,3)*ep(irec)*cd(1)+ er(irec,1,4)*es(irec)*cd(2)
C         uz = er(irec,2,1)      *cu(1)+ er(irec,2,2)      *cu(2)
C    1       + er(irec,2,3)*ep(irec)*cd(1)+ er(irec,2,4)*es(irec)*cd(2)
C         tz = er(irec,3,1)      *cu(1)+ er(irec,3,2)      *cu(2)
C    1       + er(irec,3,3)*ep(irec)*cd(1)+ er(irec,3,4)*es(irec)*cd(2)
C         tr = er(irec,4,1)      *cu(1)+ er(irec,4,2)      *cu(2)
C    1       + er(irec,4,3)*ep(irec)*cd(1)+ er(irec,4,4)*es(irec)*cd(2)
C     WRITE(6,*)'BOT: isrc,irec,ur,uz,tz,tr:',isrc,irec,ur,uz,tz,tr

      return
      end

      subroutine dobotpsv( alppsv, jbot)
c-----
c     alp = G E_N-1
c-----
      implicit none
      complex*16 alppsv(4,4)
      integer jbot

        integer NL
        parameter (NL=200)
      common /ermat/er,erinv,ep,es,elog
      complex*16 er(NL,4,4)
      complex*16 erinv(NL,4,4)
      complex*16 ep(NL),es(NL)
      real elog(NL)

        common/modlly/mmax
        integer mmax
        common/layertype/isfluid(NL),allsolid
        logical isfluid,allsolid

      common/downpsv/RD, TD
      complex*16 RD(NL,2,2), TD(NL,2,2)

      common/verby/verbose, dout
      logical verbose, dout

      integer ii,jj
      complex*16 cdet
      complex*16 a11inva12(2,2)
      complex*16 cdzero, cdone

      cdzero = dcmplx(0.0d+00,0.0d+00)
      cdone  = dcmplx(1.0d+00,0.0d+00)
c-----
c     to correcly handle the rigid free
c     assume that mmax is the layer just
c     above the boundary
c       alp = G E
c-----
      if(isfluid(mmax))then
            if(jbot.eq.1) then
c-----
c           bottom free
c-----
               alppsv(1,1) = er(mmax,3,1)
               alppsv(1,2) = er(mmax,3,3)
               alppsv(2,1) = er(mmax,2,1)
               alppsv(2,2) = er(mmax,2,3)
            else if(jbot.eq.2)then
c-----
c           bottom rigid
c-----
               alppsv(1,1) = er(mmax,2,1)
               alppsv(1,2) = er(mmax,2,3)
               alppsv(2,1) = er(mmax,3,1)
               alppsv(2,2) = er(mmax,3,3)
            else if(jbot.eq.3)then
c-----
c           bottom halfspace3
c           there is no bottom layer mmax+1 - so assume that it has the
c           same properties as the mmax layer
c           so alp = G E = I
c-----
               alppsv(1,1) = cdone
               alppsv(1,2) = cdzero
               alppsv(2,1) = cdzero
               alppsv(2,2) = cdone
            endif
c-----
c           get the R at the bottom from
c           R = - a11(inv) a12 e(- nua d N-1)
c-----
            RD(mmax,1,1) = - alppsv(1,2) * ep(mmax  ) /alppsv(1,1)
            RD(mmax,1,2) = cdzero
            RD(mmax,2,1) = cdzero
            RD(mmax,2,2) = cdzero
            TD(mmax,1,1) = cdzero
            TD(mmax,1,2) = cdzero
            TD(mmax,2,1) = cdzero
            TD(mmax,2,2) = cdzero

      else
            if(jbot.eq.1) then
c-----
c           bottom free
c-----
               alppsv(1,1) = er(mmax,3,1)
               alppsv(1,2) = er(mmax,3,2)
               alppsv(1,3) = er(mmax,3,3)
               alppsv(1,4) = er(mmax,3,4)
               alppsv(2,1) = er(mmax,4,1)
               alppsv(2,2) = er(mmax,4,2)
               alppsv(2,3) = er(mmax,4,3)
               alppsv(2,4) = er(mmax,4,4)
               alppsv(3,1) = er(mmax,1,1)
               alppsv(3,2) = er(mmax,1,2)
               alppsv(3,3) = er(mmax,1,3)
               alppsv(3,4) = er(mmax,1,4)
               alppsv(4,1) = er(mmax,2,1)
               alppsv(4,2) = er(mmax,2,2)
               alppsv(4,3) = er(mmax,2,3)
               alppsv(4,4) = er(mmax,2,4)
            else if(jbot.eq.2)then
c-----
c           bottom rigid
c-----
               alppsv(1,1) = er(mmax,1,1)
               alppsv(1,2) = er(mmax,1,2)
               alppsv(1,3) = er(mmax,1,3)
               alppsv(1,4) = er(mmax,1,4)
               alppsv(2,1) = er(mmax,2,1)
               alppsv(2,2) = er(mmax,2,2)
               alppsv(2,3) = er(mmax,2,3)
               alppsv(2,4) = er(mmax,2,4)
               alppsv(3,1) = er(mmax,3,1)
               alppsv(3,2) = er(mmax,3,2)
               alppsv(3,3) = er(mmax,3,3)
               alppsv(3,4) = er(mmax,3,4)
               alppsv(4,1) = er(mmax,4,1)
               alppsv(4,2) = er(mmax,4,2)
               alppsv(4,3) = er(mmax,4,3)
               alppsv(4,4) = er(mmax,4,4)
            else if(jbot.eq.3)then
c-----
c           bottom halfspace3
c           there is no bottom layer mmax+1 - so assume that it has the same 
c           properties as the mmax layer
c           so alp = G E = I
c-----
               alppsv(1,1) = cdone
               alppsv(1,2) = cdzero
               alppsv(1,3) = cdzero
               alppsv(1,4) = cdzero
               alppsv(2,1) = cdzero
               alppsv(2,2) = cdone
               alppsv(2,3) = cdzero
               alppsv(2,4) = cdzero
               alppsv(3,1) = cdzero
               alppsv(3,2) = cdzero
               alppsv(3,3) = cdone
               alppsv(3,4) = cdzero
               alppsv(4,1) = cdzero
               alppsv(4,2) = cdzero
               alppsv(4,3) = cdzero
               alppsv(4,4) = cdone
            endif
c-----
c           get the R at the bottom from
c            R = - a11(inv) a12 e(- nub d N-1)
c-----
            cdet = alppsv(1,1)*alppsv(2,2) 
     1            - alppsv(1,2)*alppsv(2,1)
            a11inva12(1,1)=( alppsv(2,2)*alppsv(1,3) 
     1            - alppsv(1,2)*alppsv(2,3))/cdet
            a11inva12(1,2)=( alppsv(2,2)*alppsv(1,4) 
     1            - alppsv(1,2)*alppsv(2,4))/cdet
            a11inva12(2,1)=(-alppsv(2,1)*alppsv(1,3) 
     1            + alppsv(1,1)*alppsv(2,3))/cdet
            a11inva12(2,2)=(-alppsv(2,1)*alppsv(1,4) 
     1            + alppsv(1,1)*alppsv(2,4))/cdet
            RD(mmax,1,1) = - a11inva12(1,1) * ep(mmax)
            RD(mmax,1,2) = - a11inva12(1,2) * es(mmax)
            RD(mmax,2,1) = - a11inva12(2,1) * ep(mmax)
            RD(mmax,2,2) = - a11inva12(2,2) * es(mmax)
      
c-----
c          not used but defined for neatness
c-----
            TD(mmax,1,1) = cdzero
            TD(mmax,1,2) = cdzero
            TD(mmax,2,1) = cdzero
            TD(mmax,2,2) = cdzero
         endif
            if(verbose) then
                  WRITE(6,*)'MMAX:',mmax
                  WRITE(6,*)'ALP:'
                  do ii=1,4
                     do jj=1,4
                        WRITE(6,*) ii, jj, alppsv(ii,jj)
                     enddo
                  enddo
            endif
         if(verbose)then
               WRITE(6,*)'RD11(',mmax,')=',RD(mmax,1,1)
               WRITE(6,*)'RD12(',mmax,')=',RD(mmax,1,2)
               WRITE(6,*)'RD21(',mmax,')=',RD(mmax,2,1)
               WRITE(6,*)'RD22(',mmax,')=',RD(mmax,2,2)
               WRITE(6,*)'TD11(',mmax,')=',TD(mmax,1,1)
               WRITE(6,*)'TD12(',mmax,')=',TD(mmax,1,2)
               WRITE(6,*)'TD21(',mmax,')=',TD(mmax,2,1)
               WRITE(6,*)'TD22(',mmax,')=',TD(mmax,2,2)
         endif
      return
      end

      subroutine dotoppsv( betpsv, jtop)
c-----
c     bet = H(inv) E1
c-----
      implicit none
      complex*16 betpsv(4,4)
      integer jtop

        integer NL
        parameter (NL=200)
      common /ermat/er,erinv,ep,es,elog
      complex*16 er(NL,4,4)
      complex*16 erinv(NL,4,4)
      complex*16 ep(NL),es(NL)
      real elog(NL)

        common/modlly/mmax
        integer mmax
        common/layertype/isfluid(NL),allsolid
        logical isfluid,allsolid

      common/uppsv/RU, TU
      complex*16 RU(NL,2,2), TU(NL,2,2)

      common/verby/verbose, dout
      logical verbose, dout

      integer ii,jj
      complex*16 cdet
      complex*16 b22invb21(2,2)
      complex*16 cdzero, cdone

      cdzero = dcmplx(0.0d+00,0.0d+00)
      cdone  = dcmplx(1.0d+00,0.0d+00)

c-----
c              -1
c      bet = H  E
c                 1
c-----

      if(isfluid(1))then
c-----
c         be careful in that the  FLUID E differs from the
c         ER storage
c         E11 = ER21
c         E12 = ER23
c         E21 = ER31
c         E22 = ER33
c---- 
            if(jtop.eq.1) then
c-----
c              top free
c-----
               betpsv(1,1) = er(1,2,1)
               betpsv(1,2) = er(1,2,3)
               betpsv(2,1) = er(1,3,1)
               betpsv(2,2) = er(1,3,3)
            else if(jtop.eq.2)then
c-----
c              top rigid
c-----
               betpsv(1,1) = er(1,3,1)
               betpsv(1,2) = er(1,3,3)
               betpsv(2,1) = er(1,2,1)
               betpsv(2,2) = er(1,2,3)
            else if(jtop.eq.3)then
c-----
c              top halfspace - inconsistency here since there is no 0 layer
c              parameter just assume that properties are same
c             bet = H inv E1   buz H = E, thus bet =I
c-----
               betpsv(1,1) = cdone
               betpsv(1,2) = cdzero
               betpsv(2,1) = cdzero
               betpsv(2,2) = cdone
            endif
c-----
c           get the R at the top from
c           R = - b22(inv) b21 e(- nua d N-1)
c-----
            RU(1,1,1) = - betpsv(2,1) * ep(1  ) /betpsv(2,2)
            RU(1,2,1) = cdzero
            RU(1,1,2) = cdzero
            RU(1,2,2) = cdzero
            TU(1,1,1) = cdzero
            TU(1,1,2) = cdzero
            TU(1,2,1) = cdzero
            TU(1,2,2) = cdzero

      else

            if(jtop.eq.1) then
c-----
c              top free
c-----
               betpsv(1,1) = er(1,1,1)
               betpsv(1,2) = er(1,1,2)
               betpsv(1,3) = er(1,1,3)
               betpsv(1,4) = er(1,1,4)
               betpsv(2,1) = er(1,2,1)
               betpsv(2,2) = er(1,2,2)
               betpsv(2,3) = er(1,2,3)
               betpsv(2,4) = er(1,2,4)
               betpsv(3,1) = er(1,3,1)
               betpsv(3,2) = er(1,3,2)
               betpsv(3,3) = er(1,3,3)
               betpsv(3,4) = er(1,3,4)
               betpsv(4,1) = er(1,4,1)
               betpsv(4,2) = er(1,4,2)
               betpsv(4,3) = er(1,4,3)
               betpsv(4,4) = er(1,4,4)
      
            else if(jtop.eq.2)then
c-----
c        top rigid
c-----
               betpsv(1,1) = er(1,3,1)
               betpsv(1,2) = er(1,3,2)
               betpsv(1,3) = er(1,3,3)
               betpsv(1,4) = er(1,3,4)
               betpsv(2,1) = er(1,4,1)
               betpsv(2,2) = er(1,4,2)
               betpsv(2,3) = er(1,4,3)
               betpsv(2,4) = er(1,4,4)
               betpsv(3,1) = er(1,1,1)
               betpsv(3,2) = er(1,1,2)
               betpsv(3,3) = er(1,1,3)
               betpsv(3,4) = er(1,1,4)
               betpsv(4,1) = er(1,2,1)
               betpsv(4,2) = er(1,2,2)
               betpsv(4,3) = er(1,2,3)
               betpsv(4,4) = er(1,2,4)
      
            else if(jtop.eq.3)then
c-----
c              top halfspace - inconsistency here since there is no 0 layer
c              parameter just assume that properties are same
c              bet = H inv E1   but H = E, thus bet =I
c-----
               betpsv(1,1) = cdone
               betpsv(1,2) = cdzero
               betpsv(1,3) = cdzero
               betpsv(1,4) = cdzero
               betpsv(2,1) = cdzero
               betpsv(2,2) = cdone
               betpsv(2,3) = cdzero
               betpsv(2,4) = cdzero
               betpsv(3,1) = cdzero
               betpsv(3,2) = cdzero
               betpsv(3,3) = cdone
               betpsv(3,4) = cdzero
               betpsv(4,1) = cdzero
               betpsv(4,2) = cdzero
               betpsv(4,3) = cdzero
               betpsv(4,4) = cdone
      
            endif
c-----
c           Now define the RU and TU
c-----
c-----
c           get the R at the top from
c            R = - b22(inv) b21 e(- nub d N-1)
c-----
            cdet = betpsv(3,3)*betpsv(4,4) 
     1               - betpsv(3,4)*betpsv(4,3)
            b22invb21(1,1)=( betpsv(4,4)*betpsv(3,1) 
     1               - betpsv(3,4)*betpsv(4,1))/cdet
            b22invb21(1,2)=( betpsv(4,4)*betpsv(3,2) 
     1               - betpsv(3,4)*betpsv(4,2))/cdet
            b22invb21(2,1)=(-betpsv(4,3)*betpsv(3,1) 
     1               + betpsv(3,3)*betpsv(4,1))/cdet
            b22invb21(2,2)=(-betpsv(4,3)*betpsv(3,2) 
     1               + betpsv(3,3)*betpsv(4,2))/cdet
            RU(1,1,1) = - b22invb21(1,1) * ep(1)
            RU(1,1,2) = - b22invb21(1,2) * es(1)
            RU(1,2,1) = - b22invb21(2,1) * ep(1)
            RU(1,2,2) = - b22invb21(2,2) * es(1)

c--
c    not used but defined for neatness
c-----
            TU(1,1,1) = cdzero
            TU(1,1,2) = cdzero
            TU(1,2,1) = cdzero
            TU(1,2,2) = cdzero
      endif
            if(verbose) then
                  WRITE(6,*)'MMAX:',mmax
                  WRITE(6,*)'BET:'
                  do ii=1,4
                     do jj=1,4
                        WRITE(6,*) ii, jj, betpsv(ii,jj)
                     enddo
                  enddo
            endif
            if(verbose)then
               WRITE(6,*)'RU11(',1,')=',RU(1,1,1)
               WRITE(6,*)'RU12(',1,')=',RU(1,1,2)
               WRITE(6,*)'RU21(',1,')=',RU(1,2,1)
               WRITE(6,*)'RU22(',1,')=',RU(1,2,2)
               WRITE(6,*)'TU11(',1,')=',TU(1,1,1)
               WRITE(6,*)'TU12(',1,')=',TU(1,1,2)
               WRITE(6,*)'TU21(',1,')=',TU(1,2,1)
               WRITE(6,*)'TU22(',1,')=',TU(1,2,2)
            endif
      return
      end

      subroutine doer(om,wvno)
c-----
c     form the general matrix exp(-nu d)
c     for a solid this is diagonal diag[exp(-nu_alpha d) exp)-nu_beta d)]
c     for a fluid this is just diag[exp(-nu_alpha d), 0 ]
c     This also sets up the E and Einvert matrices
c-----
      implicit none
      complex*16 om
      complex*16 wvno

        common/modlly/mmax
        integer mmax
        integer NL
        parameter (NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d,a,b, rho,qa,qb,etap,etas,frefp,frefs
        common/depref/refdep
        real refdep
        common/layertype/isfluid(NL),allsolid
        logical isfluid,allsolid
      common /ermat/er,erinv,ep,es,elog
      complex*16 er(NL,4,4)
      complex*16 erinv(NL,4,4)
      complex*16 ep(NL), es(NL)
      real elog(NL)

      complex*16 gam, gamm1, rom2

      common/debug/mlam,mmu,anu,bnu
      complex*16 mlam(NL), mmu(NL)
      complex*16 anu(NL),bnu(NL)

        common/damp/alpha,ieqex
        real alpha
        integer ieqex

      complex*16 cfac

      integer i
      complex*16 ka, kb, atna, atnb
      integer iwats
      complex*16 cdzero, cdone

      cdzero = dcmplx(0.0d+00,0.0d+00)
      cdone  = dcmplx(1.0d+00,0.0d+00)


      do i=1,mmax
            call aten(om,qa(i),qb(i),ka ,kb ,alpha,
     1          a(i),b(i),atna,atnb,iwats,
     2          frefp(i),frefs(i)) 

           anu(i) = cdsqrt(wvno*wvno - ka*ka)
           mmu(i) =  rho(i)*b(i)*b(i)*atnb*atnb
           mlam(i) =  rho(i)*(a(i)*a(i)*atna*atna 
     1        - 2.*b(i)*b(i)*atnb*atnb)
           anu(i) = cdsqrt(wvno*wvno-ka*ka)
           rom2 = rho(i)*om*om
           if(isfluid(i))then
c-----
c                  nu_beta is not used for a fluid but define it
c                  so that it is not NaN
c-----
                   bnu(i) = cdzero
                   es(i) = cdzero
           else
                   bnu(i) = cdsqrt(wvno*wvno-kb*kb)
                   es(i) = cdexp(-bnu(i)*d(i))
           endif
                   ep(i) = cdexp(-anu(i)*d(i))
           if(isfluid(i))then
c-----
c                elements to get Uz
c-----
                 er(i,2,1)    =   anu(i)
                 er(i,2,2)    =   cdzero
                 er(i,2,3)    = - anu(i)
                 er(i,2,4)    =   cdzero
c-----
c                elements to get Tz
c-----
                 cfac = -rho(i)*om*om
                 er(i,3,1) = cfac
                 er(i,3,2) = cdzero
                 er(i,3,3) = cfac
                 er(i,3,4) = cdzero
c-----
c                elements to get Ur = - rho om^2 Tz
c-----
                 er(i,1,1) = wvno
                 er(i,1,2) = cdzero
                 er(i,1,3) = wvno
                 er(i,1,4) = cdzero
c-----
c                elements to get Tr
c-----
                 er(i,4,1)    =   cdzero
                 er(i,4,2)    =   cdzero
                 er(i,4,3)    =   cdzero
                 er(i,4,4)    =   cdzero
c-----
c                for the invert we do not use the fake null
c                elements for S or for the Ur
c                This is a pseudo inverse that works for a fluid
c                er erinv  gives
c
c                | 0  0  X  0  |
c                | 0  1  0  0  | where X is never used but would be - k/ rho om^2
c                | 0  0  1  0  |
c                | 0  0  0  0  |
c-----
                 cfac = 2.*rho(i)*anu(i)*om*om
                 erinv(i,1,1) = cdzero
                 erinv(i,2,1) = cdzero
                 erinv(i,3,1) = cdzero
                 erinv(i,4,1) = cdzero

                 erinv(i,1,2) =   rho(i)*om*om/cfac
                 erinv(i,2,2) = cdzero
                 erinv(i,3,2) = - rho(i)*om*om/cfac
                 erinv(i,4,2) = cdzero

                 erinv(i,1,3) = -anu(i)/cfac
                 erinv(i,2,3) = cdzero
                 erinv(i,3,3) = -anu(i)/cfac
                 erinv(i,4,3) = cdzero

                 erinv(i,1,4) = cdzero
                 erinv(i,2,4) = cdzero
                 erinv(i,3,4) = cdzero
                 erinv(i,4,4) = cdzero
c           WRITE(6,*)'I=',I
c           WRITE(6,*)'E =',er(i,2,1), er(i,2,3),er(i,3,1),er(i,3,3)
c           WRITE(6,*)'EI=',erinv(i,1,2), erinv(i,1,3),
c    1      erinv(i,3,2),erinv(i,3,3)

           else
           gam = 2.*wvno*wvno/(kb*kb)
           gamm1 = gam - 1.

                 er(i,1,1)    =   wvno
                 er(i,1,2)    =   bnu(i)
                 er(i,1,3)    =   wvno
                 er(i,1,4)    = - bnu(i)
                 er(i,2,1)    =   anu(i)
                 er(i,2,2)    =   wvno
                 er(i,2,3)    = - anu(i)
                 er(i,2,4)    =   wvno
                 er(i,3,1)    =   rom2*gamm1
                 er(i,3,2)    =   rom2*gam*bnu(i)/wvno
                 er(i,3,3)    =   rom2*gamm1
                 er(i,3,4)    = - rom2*gam*bnu(i)/wvno
                 er(i,4,1)    =   rom2*gam*anu(i)/wvno
                 er(i,4,2)    =   rom2*gamm1
                 er(i,4,3)    = - rom2*gam*anu(i)/wvno
                 er(i,4,4)    =   rom2*gamm1
c-----
                 erinv(i,1,1) =   0.5*gam/wvno
                 erinv(i,1,2) = - 0.5*gamm1/anu(i)
                 erinv(i,1,3) = - 0.5/rom2
                 erinv(i,1,4) =   0.5*wvno/(rom2*anu(i))
                 erinv(i,2,1) = - 0.5*gamm1/bnu(i)
                 erinv(i,2,2) =   0.5*gam/wvno
                 erinv(i,2,3) =   0.5*wvno/(bnu(i)*rom2)
                 erinv(i,2,4) = - 0.5/rom2
                 erinv(i,3,1) =   0.5*gam/wvno
                 erinv(i,3,2) =   0.5*gamm1/anu(i)
                 erinv(i,3,3) = - 0.5/rom2
                 erinv(i,3,4) = - 0.5*wvno/(rom2*anu(i))
                 erinv(i,4,1) =   0.5*gamm1/bnu(i)
                 erinv(i,4,2) =   0.5*gam/wvno
                 erinv(i,4,3) = - 0.5*wvno/(rom2*bnu(i))
                 erinv(i,4,4) = - 0.5/rom2
           endif
      enddo

      return
      end

      subroutine botuppsv()
      implicit none
      common/modlly/mmax
      integer mmax
      common/bcpsv/alppsv,betpsv
      complex*16 alppsv(4,4), betpsv(4,4)

      integer NL
      parameter (NL=200)
      common /ermat/er,erinv,ep,es,elog
      complex*16 er(NL,4,4)
      complex*16 erinv(NL,4,4)
      complex*16 ep(NL),es(NL)
      real elog(NL)

      common/downpsv/RD, TD
      complex*16 RD(NL,2,2), TD(NL,2,2)

        common/layertype/isfluid(NL),allsolid
        logical isfluid,allsolid

      common/verby/verbose, dout
      logical verbose, dout

      integer i

c-----
c     march up through the layers
c-----
      do i = mmax-1,1,-1
C     WRITE(6,*)'BOTUP LAYER j=',i,isfluid(i),' above ',
C    1       isfluid(i+1)
         if(isfluid(i))then
           if(isfluid(i+1))then
C            WRITE(6,*)'call do105fluid(i)'
             call do105fluid(i)
           else
C            WRITE(6,*)'call do111(i)'
             call do111(i)
           endif
         else
           if(isfluid(i+1))then
C            WRITE(6,*)'call do117(i)'
             call do117(i)
           else
C            WRITE(6,*)'call do105solid(i)'
             call do105solid(i)
           endif
         endif
      enddo

      return
      end

      subroutine topdownpsv()
      implicit none
      common/modlly/mmax
      integer mmax
      common/bcpsv/alppsv,betpsv
      complex*16 alppsv(4,4), betpsv(4,4)

      integer NL
      parameter (NL=200)
      common /ermat/er,erinv,ep,es,elog
      complex*16 er(NL,4,4)
      complex*16 erinv(NL,4,4)
      complex*16 ep(NL),es(NL)
      real elog(NL)

      common/uppsv/RU, TU
      complex*16 RU(NL,2,2), TU(NL,2,2)

        common/layertype/isfluid(NL),allsolid
        logical isfluid,allsolid

      common/verby/verbose, dout
      logical verbose, dout

      integer i


c-----
c     march up through the layers
c-----
      do i = 1,mmax-1
C     WRITE(6,*)'TOPDOWN LAYER j=',i,isfluid(i),' above ',
C    1       isfluid(i+1)
         if(isfluid(i))then
           if(isfluid(i+1))then
C            WRITE(6,*)'call do107fluid(i)'
             call do107fluid(i)
           else
C            WRITE(6,*)'call do113(i)'
             call do113(i)
           endif
         else
           if(isfluid(i+1))then
C            WRITE(6,*)'call do119(i)'
             call do119(i)
           else
C            WRITE(6,*)'call do107solid(i)'
             call do107solid(i)
           endif
         endif
      enddo

      return
      end

      complex*16 function RDRU(m,i,j)
      integer NL
      parameter (NL=200)
      common/downpsv/RD, TD
      complex*16 RD(NL,2,2), TD(NL,2,2)

      common/uppsv/RU, TU
      complex*16 RU(NL,2,2), TU(NL,2,2)

      complex*16 eye(2,2)
      complex*16 csum
      complex*16 cdzero, cdone

      cdzero = dcmplx(0.0d+00,0.0d+00)
      cdone  = dcmplx(1.0d+00,0.0d+00)
      eye(1,1) = cdone
      eye(1,2) = cdzero
      eye(2,1) = cdzero
      eye(2,2) = cdone
 
C      WRITE(6,*)'RDRU:',m,i,j
      csum = cdzero
      do k=1,2
         csum = csum + RD(m,i,k) * RU(m,k,j)
C         WRITE(6,*)'k,RD,RU:',csum,RD(m,i,k), RU(m,k,j)
      enddo
      RDRU = eye(i,j) - csum
      return
      end

      complex*16 function RURD(m,i,j)
      integer NL
      parameter (NL=200)
      common/downpsv/RD, TD
      complex*16 RD(NL,2,2), TD(NL,2,2)

      common/uppsv/RU, TU
      complex*16 RU(NL,2,2), TU(NL,2,2)

      complex*16 eye(2,2)
      complex*16 csum
      complex*16 cdzero, cdone

      cdzero = dcmplx(0.0d+00,0.0d+00)
      cdone  = dcmplx(1.0d+00,0.0d+00)
      eye(1,1) = cdone
      eye(1,2) = cdzero
      eye(2,1) = cdzero
      eye(2,2) = cdone
 
C      WRITE(6,*)'RDRU:',m,i,j
      csum = cdzero
      do k=1,2
         csum = csum + RU(m,i,k) * RD(m,k,j)
C         WRITE(6,*)'k,RD,RU:',csum,RD(m,i,k), RU(m,k,j)
      enddo
      RURD = eye(i,j) - csum
      return
      end

      subroutine invert3x3(a)
      implicit none
c-----
c     invert a 3x3 matrix and return in the same
c     array
c-----
      complex*16 a(3,3)
c-----
c     matrix of cofactors
c-----
      complex*16 co(3,3), det

c-----
c     compute the co-factors
c-----
      co(1,1) = a(2,2)*a(3,3) - a(3,2)*a(2,3)
      co(1,2) = a(2,1)*a(3,3) - a(3,1)*a(2,3)
      co(1,3) = a(2,1)*a(3,2) - a(3,1)*a(2,2)

      co(2,1) = a(1,2)*a(3,3) - a(3,2)*a(1,3)
      co(2,2) = a(1,1)*a(3,3) - a(3,1)*a(1,3)
      co(2,3) = a(1,1)*a(3,2) - a(3,1)*a(1,2)

      co(3,1) = a(1,2)*a(2,3) - a(2,2)*a(1,3)
      co(3,2) = a(1,1)*a(2,3) - a(2,1)*a(1,3)
      co(3,3) = a(1,1)*a(2,2) - a(2,1)*a(1,2)
c-----
c     get the determinant
c-----
C     det=a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2)
C    1   -a(3,1)*a(2,2)*a(1,3)-a(3,2)*a(2,3)*a(1,1)-a(3,3)*a(2,1)*a(1,2)
      det=a(1,1)*co(1,1) - a(1,2)*co(1,2) + a(1,3)*co(1,3)

c-----
c     now form the invert
c-----
      a(1,1) =   co(1,1) / det

      a(1,2) = - co(2,1) / det
      a(1,3) =   co(3,1) / det
      a(2,1) = - co(1,2) / det
      a(2,2) =   co(2,2) / det
      a(2,3) = - co(3,2) / det
      a(3,1) =   co(1,3) / det
      a(3,2) = - co(2,3) / det
      a(3,3) =   co(3,3) / det
      return
      end

      subroutine invert2x2(a)
      implicit none
      complex*16 a(2,2)

      complex*16 cdet, b(2,2)
      integer i,j
      do i=1,2
         do j=1,2
            b(i,j) = a(i,j)
         enddo
      enddo
      cdet = b(1,1)*b(2,2) - b(1,2)*b(2,1)
      if (cdabs(cdet).lt.1.0e-7)then
           a(1,1) =   dcmplx(0.0d+00,0.0d+00)/cdet
           a(2,2) =   dcmplx(0.0d+00,0.0d+00)/cdet
           a(1,2) =   dcmplx(0.0d+00,0.0d+00)/cdet
           a(2,1) =   dcmplx(0.0d+00,0.0d+00)/cdet
      else
           a(1,1) =   b(2,2)/cdet
           a(2,2) =   b(1,1)/cdet
           a(1,2) = - b(1,2)/cdet
           a(2,1) = - b(2,1)/cdet
      endif
      return
      end

      subroutine do111(j)
      implicit none
c-----
c     layer j fluid
c     layer j+1 solid
c     Get rd td
c-----
      integer j, jp1
      integer NL
      parameter (NL=200)
      common/downpsv/RD, TD
      complex*16 RD(NL,2,2), TD(NL,2,2)

      common/uppsv/RU, TU
      complex*16 RU(NL,2,2), TU(NL,2,2)

      common /ermat/er,erinv,ep,es,elog
      complex*16 er(NL,4,4)
      complex*16 erinv(NL,4,4)
      complex*16 ep(NL), es(NL)
      real elog(NL)

      complex*16 a(3,3), y(3,2),x(3,2)
      complex*16 e11, e12, e21, e22
      integer i,k
      complex*16 cdzero, cdone

      cdzero = dcmplx(0.0d+00,0.0d+00)
      cdone  = dcmplx(1.0d+00,0.0d+00)
c-----
c     equation 1.11
c-----
      jp1 = j + 1
      e11 = er(j,2,1)
      e12 = er(j,2,3)
      e21 = er(j,3,1)
      e22 = er(j,3,3)
      y(1,1) =  e12*ep(j)
      y(2,1) =  e22*ep(j)
      y(3,1) =  cdzero

      a(1,1) = er(jp1,2,1)*ep(jp1)*RD(jp1,1,1) 
     1       + er(jp1,2,2)*es(jp1)*RD(jp1,2,1) 
     1       + er(jp1,2,3)
      a(2,1) = er(jp1,3,1)*ep(jp1)*RD(jp1,1,1) 
     1       + er(jp1,3,2)*es(jp1)*RD(jp1,2,1) 
     1       + er(jp1,3,3)
      a(3,1) = er(jp1,4,1)*ep(jp1)*RD(jp1,1,1) 
     1       + er(jp1,4,2)*es(jp1)*RD(jp1,2,1) 
     1       + er(jp1,4,3)

      a(1,2) = er(jp1,2,1)*ep(jp1)*RD(jp1,1,2) 
     1       + er(jp1,2,2)*es(jp1)*RD(jp1,2,2) 
     1       + er(jp1,2,4)
      a(2,2) = er(jp1,3,1)*ep(jp1)*RD(jp1,1,2) 
     1       + er(jp1,3,2)*es(jp1)*RD(jp1,2,2) 
     1       + er(jp1,3,4)
      a(3,2) = er(jp1,4,1)*ep(jp1)*RD(jp1,1,2) 
     1       + er(jp1,4,2)*es(jp1)*RD(jp1,2,2) 
     1       + er(jp1,4,4)

      a(1,3) = - e11
      a(2,3) = - e21
      a(3,3) = cdzero

      call invert3x3(a)
c-----
c     solve
c-----
      do i=1,3
         x(i,1) = cmplx(0.0,0.0)
         do k=1,3
            x(i,1) = x(i,1) + a(i,k)*y(k,1)
         enddo
      enddo
      do i=1,2
         do k=1,2
            RD(j,i,k) = cdzero
            TD(j,i,k) = cdzero
         enddo
      enddo
      TD(j,1,1) = x(1,1)
      TD(j,2,1) = x(2,1)
      RD(j,1,1) = x(3,1)
      return
      end
     
      subroutine do113(j)
      implicit none
c-----
c     layer j fluid
c     layer j+1 solid
c     Get ru tu
c-----
      integer j, jp1
      integer NL
      parameter (NL=200)
      common/downpsv/RD, TD
      complex*16 RD(NL,2,2), TD(NL,2,2)

      common/uppsv/RU, TU
      complex*16 RU(NL,2,2), TU(NL,2,2)

      common /ermat/er,erinv,ep,es,elog
      complex*16 er(NL,4,4)
      complex*16 erinv(NL,4,4)
      complex*16 ep(NL), es(NL)
      real elog(NL)

      complex*16 a(3,3), y(3,2),x(3,2)
      complex*16 e11, e12, e21, e22
      integer i,k,l
      complex*16 cdzero, cdone

      cdzero = dcmplx(0.0d+00,0.0d+00)
      cdone  = dcmplx(1.0d+00,0.0d+00)
c-----
c     equation 1.13
c-----
      jp1 = j + 1
      e11 = er(j,2,1)
      e12 = er(j,2,3)
      e21 = er(j,3,1)
      e22 = er(j,3,3)
      y(1,1) =  -er(jp1,2,1)*ep(jp1)
      y(2,1) =  -er(jp1,3,1)*ep(jp1)
      y(3,1) =  -er(jp1,4,1)*ep(jp1)
      y(1,2) =  -er(jp1,2,2)*es(jp1)
      y(2,2) =  -er(jp1,3,2)*es(jp1)
      y(3,2) =  -er(jp1,4,2)*es(jp1)

      a(1,1) = er(jp1,2,3)
      a(2,1) = er(jp1,3,3)
      a(3,1) = er(jp1,4,3)

      a(1,2) = er(jp1,2,4)
      a(2,2) = er(jp1,3,4)
      a(3,2) = er(jp1,4,4)

      a(1,3) = -(e11 + e12*ep(j)*RU(j,1,1))
      a(2,3) = -(e21 + e22*ep(j)*RU(j,1,1))
      a(3,3) = cdzero

      call invert3x3(a)
c-----
c     solve
c-----
      do l=1,2
          do i=1,3
             x(i,l) = cdzero
             do k=1,3
                x(i,l) = x(i,l) + a(i,k)*y(k,l)
             enddo
          enddo
      enddo
      do i=1,2
         do k=1,2
            RU(jp1,i,k) = cdzero
            TU(jp1,i,k) = cdzero
         enddo
      enddo
      RU(jp1,1,1) = x(1,1)
      RU(jp1,1,2) = x(1,2)
      RU(jp1,2,1) = x(2,1)
      RU(jp1,2,2) = x(2,2)
      TU(jp1,1,1) = x(3,1)
      TU(jp1,1,2) = x(3,2)
      return
      end

      subroutine do117(j)
      implicit none
c-----
c     layer j fluid
c     layer j+1 solid
c     Get ru tu
c-----
      integer j, jp1
      integer NL
      parameter (NL=200)
      common/downpsv/RD, TD
      complex*16 RD(NL,2,2), TD(NL,2,2)

      common/uppsv/RU, TU
      complex*16 RU(NL,2,2), TU(NL,2,2)

      common /ermat/er,erinv,ep,es,elog
      complex*16 er(NL,4,4)
      complex*16 erinv(NL,4,4)
      complex*16 ep(NL), es(NL)
      real elog(NL)

      complex*16 a(3,3), y(3,2),x(3,2)
      complex*16 e11, e12, e21, e22
      integer i,k,l
      complex*16 cdzero, cdone

      cdzero = dcmplx(0.0d+00,0.0d+00)
      cdone  = dcmplx(1.0d+00,0.0d+00)
c-----
c     equation 1.17
c-----
      jp1 = j + 1
      e11 = er(jp1,2,1)
      e12 = er(jp1,2,3)
      e21 = er(jp1,3,1)
      e22 = er(jp1,3,3)
      y(1,1) =  -er(j,2,3)*ep(j)
      y(2,1) =  -er(j,3,3)*ep(j)
      y(3,1) =  -er(j,4,3)*ep(j)
      y(1,2) =  -er(j,2,4)*es(j)
      y(2,2) =  -er(j,3,4)*es(j)
      y(3,2) =  -er(j,4,4)*es(j)

      a(1,1) = er(j,2,1)
      a(2,1) = er(j,3,1)
      a(3,1) = er(j,4,1)

      a(1,2) = er(j,2,2)
      a(2,2) = er(j,3,2)
      a(3,2) = er(j,4,2)

      a(1,3) = -( e11*ep(jp1)*RD(jp1,1,1) + e12) 
      a(2,3) = -( e21*ep(jp1)*RD(jp1,1,1) + e22) 
      a(3,3) = cmplx(0.0,0.0)

      call invert3x3(a)
c-----
c     solve
c-----
      do l=1,2
          do i=1,3
             x(i,l) = cdzero
             do k=1,3
                x(i,l) = x(i,l) + a(i,k)*y(k,l)
             enddo
          enddo
      enddo
      do i=1,2
         do k=1,2
            RD(j,i,k) = cdzero
            TD(j,i,k) = cdzero
         enddo
      enddo
      RD(j,1,1) = x(1,1)
      RD(j,1,2) = x(1,2)
      RD(j,2,1) = x(2,1)
      RD(j,2,2) = x(2,2)
      TD(j,1,1) = x(3,1)
      TD(j,1,2) = x(3,2)
      return
      end
     

      subroutine do119(j)
      implicit none
c-----
c     layer j edfluid
c     layer j+1 solid
c     Get ru tu
c-----
      integer j, jp1
      integer NL
      parameter (NL=200)
      common/downpsv/RD, TD
      complex*16 RD(NL,2,2), TD(NL,2,2)

      common/uppsv/RU, TU
      complex*16 RU(NL,2,2), TU(NL,2,2)

      common /ermat/er,erinv,ep,es,elog
      complex*16 er(NL,4,4)
      complex*16 erinv(NL,4,4)
      complex*16 ep(NL), es(NL)
      real elog(NL)

      complex*16 a(3,3), y(3,2),x(3,2)
      complex*16 e11, e12, e21, e22
      integer i,k
      complex*16 cdzero, cdone

      cdzero = dcmplx(0.0d+00,0.0d+00)
      cdone  = dcmplx(1.0d+00,0.0d+00)

c-----
c     equation 1.19
c-----
      
      jp1 = j + 1
      e11 = er(jp1,2,1)
      e12 = er(jp1,2,3)
      e21 = er(jp1,3,1)
      e22 = er(jp1,3,3)
      y(1,1) =  e11*ep(jp1)
      y(2,1) =  e21*ep(jp1)
      y(3,1) =  cmplx(0.0,0.0)

      a(1,1) = er(j,2,1) 
     1   + er(j,2,3)*ep(j)*RU(j,1,1) 
     1   + er(j,2,4)*es(j)*RU(j,2,1)
      a(2,1) = er(j,3,1) 
     1   + er(j,3,3)*ep(j)*RU(j,1,1) 
     1   + er(j,3,4)*es(j)*RU(j,2,1)
      a(3,1) = er(j,4,1) 
     1   + er(j,4,3)*ep(j)*RU(j,1,1) 
     1   + er(j,4,4)*es(j)*RU(j,2,1)

      a(1,2) = er(j,2,2) 
     1   + er(j,2,3)*ep(j)*RU(j,1,2) 
     1   + er(j,2,4)*es(j)*RU(j,2,2)
      a(2,2) = er(j,3,2) 
     1   + er(j,3,3)*ep(j)*RU(j,1,2) 
     1   + er(j,3,4)*es(j)*RU(j,2,2)
      a(3,2) = er(j,4,2) 
     1   + er(j,4,3)*ep(j)*RU(j,1,2) 
     1   + er(j,4,4)*es(j)*RU(j,2,2)

      a(1,3) = - e12
      a(2,3) = - e22
      a(3,3) =   cdzero

      call invert3x3(a)
c-----
c     solve
c-----
      do i=1,3
         x(i,1) = cdzero
         do k=1,3
            x(i,1) = x(i,1) + a(i,k)*y(k,1)
         enddo
      enddo
      do i=1,2
         do k=1,2
            RU(jp1,i,k)= cdzero
            TU(jp1,i,k)= cdzero
         enddo
      enddo
      TU(jp1,1,1) = x(1,1)
      TU(jp1,2,1) = x(2,1)
      RU(jp1,1,1) = x(3,1)
      return
      end
     
      subroutine do107fluid(j)
      implicit none
c-----
c     layer j fluid
c     layer j+1 solid
c     Get rd td
c-----
      integer j, jp1
      integer NL
      parameter (NL=200)

      common/uppsv/RU, TU
      complex*16 RU(NL,2,2), TU(NL,2,2)

      common /ermat/er,erinv,ep,es,elog
      complex*16 er(NL,4,4)
      complex*16 erinv(NL,4,4)
      complex*16 ep(NL), es(NL)
      real elog(NL)

        common/verby/verbose, dout
        logical verbose, dout

      complex*16 ff(4,4)
      integer l,m
      complex*16 f11, f12, f21, f22
      complex*16 fi11, fi12, fi21, fi22
      complex*16 cdzero, cdone

      cdzero = dcmplx(0.0d+00,0.0d+00)
      cdone  = dcmplx(1.0d+00,0.0d+00)
c-----
c     equation 1.07
c-----
      jp1 = j + 1
        f11 = er(j,2,1)
        f12 = er(j,2,3)
        f21 = er(j,3,1)
        f22 = er(j,3,3)

        fi11 = erinv(jp1,1,2)
        fi12 = erinv(jp1,1,3)
        fi21 = erinv(jp1,3,2)
        fi22 = erinv(jp1,3,3)

c       FF = Finv(j) F(jp1)

        ff(1,1) = fi11*f11 + fi12*f21
        ff(1,2) = fi11*f12 + fi12*f22
        ff(2,1) = fi21*f11 + fi22*f21
        ff(2,2) = fi21*f12 + fi22*f22

         if(verbose)then
         WRITE(6,*)'-------------------'
        WRITE(6,*)'j=',j
         WRITE(6,*)'E :',f11,f12,f21,f22
         WRITE(6,*)'EI:',fi11,fi12,fi21,fi22
               WRITE(6,*)'FF11:',ff(1,1)
               WRITE(6,*)'FF12:',ff(1,2)
               WRITE(6,*)'FF21:',ff(2,1)
               WRITE(6,*)'FF22:',ff(2,2)
               WRITE(6,*)'ep(',j,')=',ep(j)
         endif

      do l=1,2
         do m=1,2
            TU(jp1,l,m) = cdzero
            RU(jp1,l,m) = cdzero
         enddo
      enddo
      TU(jp1,1,1) = ep(jp1)/(ff(1,1) + ff(1,2)*ep(j)*RU(j,1,1))
      RU(jp1,1,1) = (ff(2,1)+ff(2,2)*ep(j)*RU(j,1,1) )*TU(jp1,1,1)
      return
      end
     
      subroutine do105fluid(j)
      implicit none
c-----
c     layer j fluid
c     layer j+1 solid
c     Get rd td
c-----
      integer j, jp1
      integer NL
      parameter (NL=200)
      common/downpsv/RD, TD
      complex*16 RD(NL,2,2), TD(NL,2,2)

      common /ermat/er,erinv,ep,es,elog
      complex*16 er(NL,4,4)
      complex*16 erinv(NL,4,4)
      complex*16 ep(NL), es(NL)
      real elog(NL)

        common/verby/verbose, dout
        logical verbose, dout

      complex*16 ee(4,4)
      integer l,m

      complex*16 e11, e12, e21, e22
      complex*16 ei11, ei12, ei21, ei22
      complex*16 cdzero, cdone

      cdzero = dcmplx(0.0d+00,0.0d+00)
      cdone  = dcmplx(1.0d+00,0.0d+00)
c-----
c     equation 1.05
c     note here we have 2x2
c
c     complications - not the E matrix was constructed to be 4x4 for 
c     a fluid and the einverse was make a pseudo inverse for the
c     sigmau sigmad source terms
c     here we form the e and ie inverse manually
c-----
      jp1 = j + 1
        e11 = er(jp1,2,1)
        e12 = er(jp1,2,3)
        e21 = er(jp1,3,1)
        e22 = er(jp1,3,3)

        ei11 = erinv(j,1,2)
        ei12 = erinv(j,1,3)
        ei21 = erinv(j,3,2)
        ei22 = erinv(j,3,3)

c       Ee = Einv(j) E(jp1)

        ee(1,1) = ei11*e11 + ei12*e21
        ee(1,2) = ei11*e12 + ei12*e22
        ee(2,1) = ei21*e11 + ei22*e21
        ee(2,2) = ei21*e12 + ei22*e22

         if(verbose)then
         WRITE(6,*)'-------------------'
         WRITE(6,*)'j=',j
         WRITE(6,*)'E :',e11,e12,e21,e22
         WRITE(6,*)'EI:',ei11,ei12,ei21,ei22
               WRITE(6,*)'EE11:',ee(1,1)
               WRITE(6,*)'EE12:',ee(1,2)
               WRITE(6,*)'EE21:',ee(2,1)
               WRITE(6,*)'EE22:',ee(2,2)
               WRITE(6,*)'ep(',j,')=',ep(j)
         endif
              do l=1,2
         do m=1,2
            TD(j,l,m) = cdzero
            RD(j,l,m) = cdzero
         enddo
      enddo
      TD(j,1,1) = ep(j)/(ee(2,1)*ep(jp1)*RD(jp1,1,1) + ee(2,2))
      RD(j,1,1) = (ee(1,1)*ep(jp1)*RD(jp1,1,1) + ee(1,2))*td(j,1,1)
      return
      end
     
      subroutine do105solid(j)
      implicit none
c-----
c     layer j fluid
c     layer j+1 solid
c     Get rd td
c-----
      integer j, jp1
      integer NL
      parameter (NL=200)
      common/downpsv/RD, TD
      complex*16 RD(NL,2,2), TD(NL,2,2)

      common /ermat/er,erinv,ep,es,elog
      complex*16 er(NL,4,4)
      complex*16 erinv(NL,4,4)
      complex*16 ep(NL), es(NL)
      real elog(NL)

      complex*16 ee(4,4)
      integer k,l,m

      complex*16 cmat(2,2)
      complex*16 EE11(2,2), EE12(2,2), EE21(2,2), EE22(2,2)
      complex*16 cdzero, cdone

      cdzero = dcmplx(0.0d+00,0.0d+00)
      cdone  = dcmplx(1.0d+00,0.0d+00)
c-----
c     equation 1.05
c     note here we have 2x2
c-----
      jp1 = j + 1
      do k=1,4
         do l=1,4
            ee(k,l) = cdzero
            do m=1,4
               ee(k,l) = ee(k,l) + erinv(j,k,m)*er(jp1,m,l)
            enddo 
         enddo 
       enddo 
c-----
c       break the 4x4 into 4 2x2 partitions
c-----
        EE11(1,1) = ee(1,1)
        EE11(1,2) = ee(1,2)
        EE11(2,1) = ee(2,1)
        EE11(2,2) = ee(2,2)

        EE12(1,1) = ee(1,3)
        EE12(1,2) = ee(1,4)
        EE12(2,1) = ee(2,3)
        EE12(2,2) = ee(2,4)

        EE21(1,1) = ee(3,1)
        EE21(1,2) = ee(3,2)
        EE21(2,1) = ee(4,1)
        EE21(2,2) = ee(4,2)

        EE22(1,1) = ee(3,3)
        EE22(1,2) = ee(3,4)
        EE22(2,1) = ee(4,3)
        EE22(2,2) = ee(4,4)

c-----
c       now define cmat = (EE21 e(ijp)RD(ijp) + EE22)
c-----
        cmat(1,1) = EE21(1,1)*ep(jp1)*RD(jp1,1,1)
     1            + EE21(1,2)*es(jp1)*RD(jp1,2,1) + EE22(1,1)
        cmat(1,2) = EE21(1,1)*ep(jp1)*RD(jp1,1,2)
     1            + EE21(1,2)*es(jp1)*RD(jp1,2,2) + EE22(1,2)
        cmat(2,1) = EE21(2,1)*ep(jp1)*RD(jp1,1,1)
     1            + EE21(2,2)*es(jp1)*RD(jp1,2,1) + EE22(2,1)
        cmat(2,2) = EE21(2,1)*ep(jp1)*RD(jp1,1,2)
     1            + EE21(2,2)*es(jp1)*RD(jp1,2,2) + EE22(2,2)
c-----
c       get the invert of cmat, use simpler variables for
c       cleaner code
c-----
        call invert2x2(cmat)

        TD(j,1,1) = cmat(1,1)*ep(j)
        TD(j,1,2) = cmat(1,2)*es(j)
        TD(j,2,1) = cmat(2,1)*ep(j)
        TD(j,2,2) = cmat(2,2)*es(j)
c-----
c       now define cmat = (EE11 e(jp1)RD(jp1) + EE12)
c-----
        cmat(1,1) = EE11(1,1)*ep(jp1)*RD(jp1,1,1)
     1            + EE11(1,2)*es(jp1)*RD(jp1,2,1) + EE12(1,1)
        cmat(1,2) = EE11(1,1)*ep(jp1)*RD(jp1,1,2)
     1            + EE11(1,2)*es(jp1)*RD(jp1,2,2) + EE12(1,2)
        cmat(2,1) = EE11(2,1)*ep(jp1)*RD(jp1,1,1)
     1            + EE11(2,2)*es(jp1)*RD(jp1,2,1) + EE12(2,1)
        cmat(2,2) = EE11(2,1)*ep(jp1)*RD(jp1,1,2)
     1            + EE11(2,2)*es(jp1)*RD(jp1,2,2) + EE12(2,2)

        RD(j,1,1) = cmat(1,1)*TD(j,1,1) + cmat(1,2)*TD(j,2,1)
        RD(j,1,2) = cmat(1,1)*TD(j,1,2) + cmat(1,2)*TD(j,2,2)
        RD(j,2,1) = cmat(2,1)*TD(j,1,1) + cmat(2,2)*TD(j,2,1)
        RD(j,2,2) = cmat(2,1)*TD(j,1,2) + cmat(2,2)*TD(j,2,2)
      return
      end

      subroutine do107solid(j)
      implicit none
c-----
c     layer j fluid
c     layer j+1 solid
c     Get rd td
c-----
      integer j, jp1
      integer NL
      parameter (NL=200)

      common/uppsv/RU, TU
      complex*16 RU(NL,2,2), TU(NL,2,2)

      common /ermat/er,erinv,ep,es,elog
      complex*16 er(NL,4,4)
      complex*16 erinv(NL,4,4)
      complex*16 ep(NL), es(NL)
      real elog(NL)

      complex*16 ff(4,4)
      complex*16 FF11(2,2), FF12(2,2), FF21(2,2), FF22(2,2)
      complex*16 cmat(2,2), csum

      integer ii, jj , kk
      complex*16 cdzero, cdone

      cdzero = dcmplx(0.0d+00,0.0d+00)
      cdone  = dcmplx(1.0d+00,0.0d+00)
c-----
c     equation 1.07
c-----
      jp1 = j + 1
c-----
c             -1
c       form E    E
c             j+1  j
c-----
        do ii=1,4
          do jj=1,4
             csum = cdzero
             do kk=1,4
                csum = csum + erinv(jp1,ii,kk)*er(j,kk,jj)
             enddo
             ff(ii,jj) = csum
          enddo
        enddo
c-----
c       break the 4x4 into 4 2x2 partitions
c-----
        FF11(1,1) = ff(1,1)
        FF11(1,2) = ff(1,2)
        FF11(2,1) = ff(2,1)
        FF11(2,2) = ff(2,2)

        FF12(1,1) = ff(1,3)
        FF12(1,2) = ff(1,4)
        FF12(2,1) = ff(2,3)
        FF12(2,2) = ff(2,4)

        FF21(1,1) = ff(3,1)
        FF21(1,2) = ff(3,2)
        FF21(2,1) = ff(4,1)
        FF21(2,2) = ff(4,2)

        FF22(1,1) = ff(3,3)
        FF22(1,2) = ff(3,4)
        FF22(2,1) = ff(4,3)
        FF22(2,2) = ff(4,4)
c-----
c       now define cmat = (FF11 + FF12 e(i)RU(i) )
c-----
        cmat(1,1) = FF11(1,1) + FF12(1,1)*ep(j)*RU(j,1,1) 
     1                        + FF12(1,2)*es(j)*RU(j,2,1) 
        cmat(1,2) = FF11(1,2) + FF12(1,1)*ep(j)*RU(j,1,2) 
     1                        + FF12(1,2)*es(j)*RU(j,2,2) 
        cmat(2,1) = FF11(2,1) + FF12(2,1)*ep(j)*RU(j,1,1) 
     1                        + FF12(2,2)*es(j)*RU(j,2,1) 
        cmat(2,2) = FF11(2,2) + FF12(2,1)*ep(j)*RU(j,1,2) 
     1                        + FF12(2,2)*es(j)*RU(j,2,2) 
c-----
c       get the invert of cmat, use simpler variables for
c       cleaner code
c-----
        call invert2x2(cmat)

        TU(jp1,1,1) = cmat(1,1)*ep(jp1)
        TU(jp1,1,2) = cmat(1,2)*es(jp1)
        TU(jp1,2,1) = cmat(2,1)*ep(jp1)
        TU(jp1,2,2) = cmat(2,2)*es(jp1)
c-----
c       now define cmat = (FF21 + FF22 e(i)RU(i) ) TU(jp1)
c-----
        cmat(1,1) = FF21(1,1) + FF22(1,1)*ep(j)*RU(j,1,1) 
     1                        + FF22(1,2)*es(j)*RU(j,2,1) 
        cmat(1,2) = FF21(1,2) + FF22(1,1)*ep(j)*RU(j,1,2) 
     1                        + FF22(1,2)*es(j)*RU(j,2,2) 
        cmat(2,1) = FF21(2,1) + FF22(2,1)*ep(j)*RU(j,1,1) 
     1                        + FF22(2,2)*es(j)*RU(j,2,1) 
        cmat(2,2) = FF21(2,2) + FF22(2,1)*ep(j)*RU(j,1,2) 
     1                        + FF22(2,2)*es(j)*RU(j,2,2) 

        RU(jp1,1,1) = cmat(1,1)*TU(jp1,1,1) + cmat(1,2)*TU(jp1,2,1)
        RU(jp1,1,2) = cmat(1,1)*TU(jp1,1,2) + cmat(1,2)*TU(jp1,2,2)
        RU(jp1,2,1) = cmat(2,1)*TU(jp1,1,1) + cmat(2,2)*TU(jp1,2,1)
        RU(jp1,2,2) = cmat(2,1)*TU(jp1,1,2) + cmat(2,2)*TU(jp1,2,2)

      return
      end

        subroutine aten(om,qa,qb,xka,xkb,alpha,a,b,atna,atnb,iwat,
     2      frefp,frefs)
        implicit none
c-----
c       make velocities complex*16, using Futterman causality operator
c-----
        real qa,qb,alpha,a,b,frefp,frefs
        complex*16 om,atna,atnb,xka,xkb
        integer iwat
        complex*16 at

        common/cntrl/ishank,hnkarg,dstcor,dokjar
        integer dstcor
        real hnkarg
        logical ishank,dokjar

        real pi, om1p, om1s, oml, fac, pi2
        real gama, gamb, rfac
      complex*16 cdzero, cdone

      cdzero = dcmplx(0.0d+00,0.0d+00)
      cdone  = dcmplx(1.0d+00,0.0d+00)
c-----
c       reference frequency is fref hz
c-----
        om1p=6.2831853*frefp
        om1s=6.2831853*frefs
        pi2 = 1.5707963
        pi=3.1415927
        if(dokjar)then
c-----
c       Kjartansson Constant Q, causal Q operator
c       Kjartansson, E. (1979). 
c          Constant Q-wave propagation and attenuation,
c       J. Geophys. Res. 84, 4737-4748.
c-----
            gama = atan(qa)/pi
            gamb = atan(qb)/pi
            if(gama.le.0.0)then
                atna = cdone
            else
                fac = pi2*gama
                rfac = sin(fac)/cos(fac)
                atna = cdone/
     1              (( (om/om1p)**(-gama) ) *
     2               dcmplx(1.00d+00,-dble(rfac)))
            endif
            if(b.gt.1.0e-04*a)then
                if(gamb.le.0.0)then
                    atnb = cdone
                else
                    fac = pi2*gamb
                    rfac = sin(fac)/cos(fac)
                    atnb = cmplx(1.00,0.00)/
     1              (( (om/om1s)**(-gamb) ) *
     2               dcmplx(1.0d+00,-dble(rfac)))
                endif
            endif
        else
c-----
c       Futterman Causal Q
c-----
c           low frequency cutoff is 0.01 hz
c-----
            oml=0.062831853
            atna=cdone
            atnb=cdone
            if(qa.gt.0.0)then
              at=cdzero
              if(cdabs(om).gt.oml) at=cdlog(om/om1p)/pi
              if(cdabs(om).le.oml) then
                fac=sqrt(oml*oml + (alpha*alpha))/oml
              at=cdlog(dcmplx(dble(oml),-dble(alpha))/(om1p*fac))/pi
              endif
              atna=(1.+(qa)*at+dcmplx(0.00d+00,dble(qa/2.)))
            endif
            if(qb.gt.0.0)then
              at=cdzero
              if(cdabs(om).gt.oml) at=cdlog(om/om1s)/pi
              if(cdabs(om).le.oml) then
                fac=sqrt(oml*oml + (alpha*alpha))/oml
              at=cdlog(dcmplx(dble(oml),-dble(alpha))/(om1s*fac))/pi
              endif
               atnb=(1.+(qb)*at+dcmplx(0.00d+00,dble(qb/2.)))
            endif
        endif
        xka=om/((a)*atna)
        if(b.le.1.0e-04*a)then
            iwat = 1
            xkb = cdzero
        else
            iwat = 0
            xkb=om/((b)*atnb)
        endif
        return
        end

