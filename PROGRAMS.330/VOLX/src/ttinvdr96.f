        program ttinvdr96
c----------------------------------------------------------------------c
c                                                                      c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                 c
c      VOLUME X                                                        c
c                                                                      c
c      PROGRAM: TTINVDR96                                               c
c                                                                      c
c      COPYRIGHT 2011                                                  c
c      R. B. Herrmann                                   c
c      Department of Earth and Atmospheric Sciences                    c
c      Saint Louis University                                          c
c      221 North Grand Boulevard                                       c
c      St. Louis, Missouri 63103                                       c
c      U. S. A.                                                        c
c                                                                      c
c----------------------------------------------------------------------c
c       This program checks the input file
c       tobs.d to process a list of phase, distance and travel times 
c       to output partial derivatives
c
c       NOTE: the program ttinvpr96 has already checked the files for
c       their validity
c-----
c       CHANGES
c       26 JAN 2011 - created
c
c       Temporary files
c       tmpsrfi.20  has the observations
c       tmpsrfi.21 has the partials
c       
c       tmpsrti.17 is the current model
c-----
        implicit none

        integer LOT
        parameter (LOT=6)
        character fname*80

        logical dop
        character mname*80

        integer mnmarg

        integer kerr, ssytem
        integer ls
        integer lgstr

c-----
c       observation file tmpsrfi.20
c-----
        integer iphase
        real dist, time, err, evdp, stdp

c-----
c       control file eventually get from gttmp0 tmpsrfi.00 and tmpsrfi.16
c-----
        integer invdep
        real twnmin,twnmax
        integer idtwo
        
c-----
c       call machine dependent initialization
c-----
        call mchdep()
c-----
c       parse command line arguments
c-----
c       If there are any command line arguments, program is run
c       only under command line control, else the program is
c       run under the command file control as part of the inversion
c-----
        if(mnmarg() .gt.0)then
c-----
c           set up default earth model file
c           can be overridden by the command line
c-----
            mname = 'tmpsrfi.17'
            call gcmdln(dop,mname)
c-----
c           force incident P 
c-----
            dop = .true.
c-----
c           compute partials with respect to layer thickness
c-----
C            call doit(' ', n, dt, rayp, gaussalp, delay,1.0,0,
C     1          twnmin,twnmax,dop,dotwo,mname,
C     2              iout,outbin,kstnm, nzyear, nzjday, 
C     3              nzhour, nzmin, nzmon, nzday)
c-----
c           compute partials with respect to layer velocities
c-----
C            call doit(' ', n, dt, rayp, gaussalp, delay,1.0,1,
C     1          twnmin,twnmax,dop,dotwo,mname,
C     2              iout,outbin,kstnm, nzyear, nzjday, 
C     3              nzhour, nzmin, nzmon, nzday)
        else
c-----
c       read in control file information concerning the variance? and
c       number of data files and the actual data window to be used
c       and whether we invert for velocity or for layer thickness
c-----
            open(1,file='tmpsrfi.16',status='unknown',
     1          form='formatted',access='sequential')
            rewind 1
            read(1,'(i5,2f20.10,i5)')invdep,twnmin,twnmax,idtwo
            close (1)
        if(invdep.eq.1)then
            write(LOT,1)
        else if(invdep.eq.0)then
            write(LOT,2)
        endif
    1   format('Processing First Arrival Time Partials for', 
     1   ' Layer Velocity ')
    2   format('Processing First Arrival Time Partials for', 
     1   ' Layer Thickness')

c-----
c           open binary file for the partial derivatives
c           open tmpsrfi.21
c           output for observation
c               TYPE DTDV OBS PRED ERR RES
c           or
c               TYPE DTDH OBS PRED ERR RES
c-----
            open(4,file='tmpsrfi.21',status='unknown',
     1          form='unformatted',access='sequential')
c-----
c           force incident P rftn, RFTN, double length FFT
c-----
            mname = 'tmpsrfi.17'
            
c-----
c           open observation file
c-----

            open(1,file='tmpsrfi.20',access='sequential',
     1          form='formatted',status='unknown')
            rewind 1

 1000       continue
            read(1,*,end=9000)iphase, dist, time, err, evdp, stdp
            call TTVH(mname,dist,iphase,time,err,invdep,evdp,stdp,4)
        go to 1000
 9000       continue
            close (1)
            rewind 4
            close(4)
c-----
c       end of COMMAND LINE/COMMAND FILE PROCESSING
c-----
        endif
        end


        subroutine gcmdln(dop,mname)
c-----
c       parse command line arguments
c       requires subroutine mgtarg() and function mnmarg()
c-----
c       dop L   - .true. P-wave incident
c                 .false. SV-wave incident
c       mname   Ch* - earth model name - required
c-----
        logical dop
        character mname*(*)

        integer*4 mnmarg
        character*50 name

        dop = .true.
C-----
c       no default for the model name - we set this externally
c----
C       mname = ' '
        nmarg=mnmarg()
        i = 0
   11   i = i + 1
            if(i.gt.nmarg)goto 13
            call mgtarg(i,name)
            if(name(1:2).eq.'-P' )then
                dop = .true.
            else if(name(1:2).eq.'-S' )then
                dop = .false.
            else if(name(1:2).eq.'-M')then
                i = i + 1
                call mgtarg(i,mname)
            else if(name(1:2).eq.'-?')then
                call usage(' ')
            else if(name(1:2).eq.'-h')then
                call usage(' ')
            endif
        goto 11
   13   continue
        return
        end

        subroutine usage(ostr)
c------
c       write out program syntax
c-----
        character ostr*(*)
        integer LER,LOT,LIN
        parameter(LER=0, LIN=5, LOT=6)
        if(ostr.ne. ' ' )then
            lostr = lgstr(ostr)
            write(LER,*)ostr(1:lostr)
        endif
        write(LER,*)'USAGE: ',
     1  'ttinvdr96 [-P] [-S] ',
     1      ' -M model '
        write(LER,*)
     1  '-P           (default true )    Incident P wave'
        write(LER,*)
     1  '-S           (default false)    Incident S wave'
        write(LER,*)
     1  '-M   model   (default none )    Earth model name'
        write(LER,*)
     1  '-?                   Display this usage message'
        write(LER,*)
     1  '-h                   Display this usage message'
        stop 
        end



       subroutine TTVH(mname,dist,iphase,time,err,invdep,evdp,stdp,lun)
c------
c      iphase, (integer) 1 is P wave, 2 is S wave
c      invdep, (integer) 1 is derivative with velocity, 
c                        0 is derivative with thickness
c      time, is observed first arrival
c      err, is observed erro.
c      lun,(integer) for print out
c------
       implicit none
       integer LER, LOT, LIN
       parameter(LER=0, LIN=5, LOT=6)
c-----
c      compute the partial derivatives and output onto unit
c      corresponding to tmpsrfi.21
c-----
c-----
c      subroutine arguments
c-----
       character mname*(*)
       real R
c-----
c       earth model information
c-----
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL) 
        real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        integer mmax
        common/modlly/mmax
        integer iunit, iiso, iflsph, idimen, icnvel, ierr
        character title*80
c-----
c      internal variables
c-----
       integer  i,j,T_F

       real S(NL), T_R(NL), DTDV(NL), DTDH(NL)
       real TT

       real tmp, sum

       integer lgstr
       integer lm, lf
        
       logical ext
       integer invdep,lun
c-----
c       observation file tmpsrfi.20
c-----
       integer iphase
       real dist, time, err, evdp, stdp

c-----
c       get the current earth model
c-----
       lm = lgstr(mname)
       inquire(file=mname,exist=ext)
       if(.not. ext)then 
           write(LER,*)'Model file does not exist:',mname(1:lm)
           return
       endif
       call getmod(2,mname,mmax,title,iunit,iiso,iflsph,
     1         idimen,icnvel,ierr,.false.)
       if(ierr .lt. 0)return
c-----
c       make sure that we use 1/Q
c-----
       do 3007 i=1,mmax
            if(qa(i).lt.0.0)qa(i) = 0.0
            if(qb(i).lt.0.0)qb(i) = 0.0
            if(qa(i) .gt. 1.0)qa(i) = 1.0/qa(i)
            if(qb(i) .gt. 1.0)qb(i) = 1.0/qb(i)
            if(frefp(i) .le. 0.0)frefp(i) = 1.0
            if(frefs(i) .le. 0.0)frefs(i) = 1.0
 3007  continue
 
      if(iphase .eq. 0) then
c-----get the slowness of P wave
      do i=1,mmax
         DTDV(i)=0.
         T_R(i)=0.
         DTDH(i)=0.
         S(i)=1.0/a(i)
      enddo              
      call DVDT(mmax,S,d,TT,DTDV,DTDH,dist)
      else
c-----get the slowness of S wave
      do i=1,mmax
         DTDV(i)=0.
         T_R(i)=0.
         DTDH(i)=0.
         S(i)=1.0/b(i)
      enddo

      call DVDT(mmax,S,d,TT,DTDV,DTDH,dist)
      endif

c------invdep = 0 , out put partial derivative with velocity
      if(invdep .eq. 1) then
          write(lun) iphase,dist, (DTDV(i),i=1,mmax),time,TT,err,
     1        evdp, stdp,time-TT
c------invdep = 1 , out put partial derivative with thickness
      elseif(invdep .eq. 0) then
          write(lun) iphase,dist, (DTDH(i),i=1,mmax),time,TT,err,
     1        evdp, stdp,time-TT
      endif

      end


c-----
c     mmax is the number of layer
c     S is the slowness of Vp or Vs
c     d is the thickness of layer
c     TT is the trval time of head wave
c     DTDV is the partial derivative with velocity
c     dist is the epicenter distance
c-----
      subroutine DVDT(mmax,S,d,TT,DTDV,DTDH,dist)

c-----
c     internal variables
c     T_F index the head wave is from T_F's layer
c     T_R is the refraction travel time from the i's layer
c------
       integer  i,j,T_F,mmax

       integer NL
       parameter(NL=200)
       real S(NL), T_R(NL), DTDV(NL), DTDH(NL),d(NL), TT, dist,tmp

c-----compute the direct arrival time for reference
       T_R(1)=dist*S(1)        
       TT=T_R(1)

c-----computer the refraction travel time, if there it is
        do i=2,mmax
            tmp=0.0

            if (S(i) .LT. S(i-1)) then
              do  j=1,i-1
                    if(S(j) .LT. S(i)) then
                     T_R(i) = 0
c                     break
c                    goto 100
                    else
                      tmp=tmp + 2 * d(j)*sqrt(S(j)*S(j) - S(i)*S(i))
                    endif
              enddo
                     T_R(i)= tmp + dist*S(i)
            else
100           T_R(i) = 0
              continue
            endif
        enddo
c-----Need consideration: there is a warning when using goto 100
c       

c-----set the head wave index, indicate it is from which layer
        T_F=1

c-----compare the refraction and direct arrive time to get the head wave
        do  i=2,mmax
            if(T_R(i) .GT. TT) then
                continue
            elseif(T_R(i) .EQ. 0) then
                continue
            else
                TT=T_R(i)
                T_F=i
            endif
        enddo
 
c-----computer the derivative of travel time 

        if(T_F .EQ. 1) then
c-----if the head wave is direct arrive
            DTDV(1)=-1*dist*S(1)*S(1)
            DTDH(1)=0
        else
c-----if the head wave is refraction 
            tmp=0.0;
            do j=1,T_F-1           
                 DTDV(j)=2*d(j)*S(j)/sqrt(S(j)*S(j)-S(T_F)*S(T_F))
                 tmp=tmp-2*d(j)*S(T_F)/sqrt(S(j)*S(j)-S(T_F)*S(T_F))
                 DTDH(j)=2*sqrt(S(j)*S(j)-S(T_F)*S(T_F))
               DTDV(j) = DTDV(j) * (-S(j)*S(j))
            enddo
            DTDV(T_F)= dist + tmp
            DTDV(T_F)=DTDV(T_F) * (-S(T_F)*S(T_F))
            DTDH(T_F)=0
        endif

       end
