      implicit double precision (a-h,o-z)
      parameter (maxlay=50)
      parameter (twopi=6.2831853072d0)
      real*8 xtoplay(maxlay)
      real*8 xbotlay(maxlay)
      integer ifluid(maxlay)
      integer ianiso(maxlay)
      integer numlev(maxlay)
      parameter (maxpts=1000)
      real*8 runradius(maxpts)
      real*8 rundelta(maxpts)
      real*8 runtime(maxpts)
      character*80 emodel
      character*80 ttfile
c
c      parameter (maxp=20500) 
c      parameter (maxp=14000)
        parameter (maxp=16000)
      real*4 p4(maxp)
      real*4 d4(maxp)
      real*4 t4(maxp)
      real*4 dddp4(maxp)
      real*4 tstar4(maxp)
      real*4 xturn4(maxp)
      real*4 top4
      real*4 bot4
      dimension ibyteray(maxp)
      dimension nraypnts(maxp)
c
      logical isotropic
      logical skip
c
      common /plevel/iprtlv
      iprtlv=0
c
c---- read in an earth model
c
      write(6,"('type name of model to read')")
      read(5,"(a)") emodel
      open(1,file=emodel)
      call readem(1,maxlay,numlay,xbotlay,xtoplay,numlev,ifluid,ianiso,ierr)
      close(1)
      do i=1,numlay
        write(6,"(i3,f10.1,f10.1,3i4)") i,xbotlay(i),xtoplay(i),
     #                                  ifluid(i),ianiso(i),numlev(i)
      enddo
c
c---- open an output file
c
      ttfile=emodel(1:lnblnk(emodel))//'.RAYS'
      open(1,file=ttfile)
      close(1)
      call openfl(1,ttfile,4,0,0,ierr,-1)
c
c---- write out the earth model
c
      call writeemfl(1,ierr)
c
c---- set up ray parameter array
c
      np=maxp
      nbytes=np*4
      call bffo(1,1,np,4,ierr,0)
      ip=1
      p4(ip)=0.
      do while (ip.lt.np)
        ip=ip+1
	if(ip.gt.maxp) then
	  write(6,"('ERROR: ip.gt.maxp')")
	  stop 
	endif
        deltap=0.0026d0
        if(p4(ip-1).gt.4.4.and.p4(ip-1).lt.4.9) then
          deltap=0.0005d0
        endif
        if(p4(ip-1).gt.8.3.and.p4(ip-1).lt.9.0) then
          deltap=0.0005d0
        endif
        if(p4(ip-1).gt.15.2.and.p4(ip-1).lt.15.7) then
          deltap=0.0005d0
        endif
        if(p4(ip-1).gt.9.6.and.p4(ip-1).lt.10.1) then
          deltap=0.0005d0
        endif
           if(p4(ip-1).gt.18.0.and.p4(ip-1).lt.18.5) then
          deltap=0.0005d0
        endif
        p4(ip)=p4(ip-1)+sngl(deltap)
      enddo
      p4(1)=0.001d0
      call bffo(1,1,p4,nbytes,ierr,0)
      write(6,"('wrote out ray parameters',g15.5)") p4(np)
c
c---- loop through calculation twice, first time to find out ending position
c
      call getflpos(1,ipospred,ierr)
      itabpos=ipospred+5
      iraypos=0
      nlaytabs=0
c
      do iloop=1,2
c
        if(iloop.eq.2) then
          call bffo(1,1,nlaytabs,4,ierr,0)
          iraypos=ipospred+5
        endif
        nonzero=0
        ntotpts=0
c
c---- loop on ray types: iani=0,1; iraytype=1,2,3 (SH,P,SV)
c
        do iani=0,1
          do iraytype=1,3
c
c---- loop on model layers and source sublayers
c
            do ilay=numlay,1,-1
              skip=.false.
              isotropic=.false.
              if(iani.eq.0) then
                isotropic=.true.
                if(iraytype.eq.3) then
                  skip=.true.
                endif
              else if(iani.eq.1) then
                if(ianiso(ilay).eq.1) then
                else 
                  skip=.true.
                endif
              endif
              nsub=1+idint((xtoplay(ilay)-xbotlay(ilay)-1.d0)/10.d0)
              do isub=1,nsub
                xbot=xbotlay(ilay)
                xtop=xtoplay(ilay)-dfloat(isub-1)*10.d0
                if(xtop.lt.(6371.d0-800.d0).and.isub.gt.1) then
                  skip=.true.
                endif
                if(.not.skip) then
                  nlev=numlev(ilay)
c
c---- loop on ray parameter
c
                  npout=0
                  do ip=1,np
                    pray=dble(p4(ip))*360.0d0/twopi
                    call wgray(pray,ilay,xbot,xtop,nlev,isotropic,iraytype,
     #                       delta,ttime,dddp,tstar,xturn,
     #                       maxpts,npts,runradius,rundelta,runtime)
                    t4(ip)=sngl(ttime)
                    d4(ip)=sngl(delta)
                    dddp4(ip)=sngl(dddp)
                    tstar4(ip)=sngl(tstar)
                    xturn4(ip)=sngl(xturn)
                    if(ttime.gt.0.001d0) then
                      npout=npout+1
                    endif
                    ibyteray(ip)=0
                    nraypnts(ip)=0
                    if(npts.gt.0) then
                      nonzero=nonzero+1
                      ntotpts=ntotpts+npts
                      if(iloop.eq.2) then
                        ibyteray(ip)=iraypos
                        nraypnts(ip)=npts
                        nrpts=8*npts
c                        call bffo(1,1,npts,4,ierr,iraypos)
c                        call bffo(1,1,runradius,nrpts,ierr,0)
c                        call bffo(1,1,rundelta,nrpts,ierr,0)
                        iraypos=iraypos+4+nrpts+nrpts
                      endif
                    endif
                  enddo
c
c---- write out the results
c
                  top4=xtop
                  bot4=xbot
                  if(iloop.eq.2) then
                    call bffo(1,1,ilay,4,ierr,itabpos)
                    call bffo(1,1,iani,4,ierr,0)
                    call bffo(1,1,iraytype,4,ierr,0)
                    call bffo(1,1,isub,4,ierr,0)
                    call bffo(1,1,top4,4,ierr,0)
                    call bffo(1,1,bot4,4,ierr,0)
                    call bffo(1,1,npout,4,ierr,0)
c
                    call bffo(1,1,t4,nbytes,ierr,0)
                    call bffo(1,1,d4,nbytes,ierr,0)
                    call bffo(1,1,dddp4,nbytes,ierr,0)
                    call bffo(1,1,tstar4,nbytes,ierr,0)
                    call bffo(1,1,xturn4,nbytes,ierr,0)
                    call bffo(1,1,nraypnts,nbytes,ierr,0)
                    itabpos=itabpos+7*4+6*nbytes
                  else
                    ipospred=ipospred+7*4+6*nbytes
                    nlaytabs=nlaytabs+1
                  endif
                  write(6,"(4i4,4f12.5,i4,i10)") iani,ilay,isub,iraytype,top4,bot4,t4(1),d4(1),nraypnts(1),ibyteray(1)
                else
                endif
              enddo
            enddo
          enddo
        enddo
        write(6,"(i7,' rays, ',i9,' points')") nonzero,ntotpts
        call getflpos(1,iposition,ierr)
        write(6,"('file position:',2i10)") iposition,ipospred
        write(6,"('itabpos,iraypos:',2i11)") itabpos,iraypos
      enddo
      call closfl(1,ierr)
      end
