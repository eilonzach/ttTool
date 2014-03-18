c---- this program test the subroutine cagcrays
c
      character*80 filename
      character*16 phase
c      
      parameter (maxarr=20)
      dimension ttime(maxarr)
      dimension slowness(maxarr)
      dimension dtddepth(maxarr)
      dimension ddeltadp(maxarr)
      dimension tstar(maxarr)
      dimension ellcorr(maxarr)
c
      parameter (maxraypts=5000)
      dimension nraypts(maxarr)
      dimension rayrad(maxraypts,maxarr)
      dimension raydelta(maxraypts,maxarr)
      dimension raylat(maxraypts,maxarr)
      dimension raylong(maxraypts,maxarr)
      dimension raytime(maxraypts,maxarr)
      dimension iraytp(maxraypts,maxarr)
      dimension iraylay(maxraypts,maxarr)
      logical exists
c
      common/plevel/iprtlv
c
c-------------------------------------------------------------------------------
c
	write(6,"('------------------------')")
	call system("date")

      lu=1


	write(6,"('type name of the .RAYS file: ')")
	read(5,"(a)") filename
	inquire(file=filename,exist=exists)
	if(.not.exists) then
	  write(6,"(a,' does not exists')")filename(1:lnblnk(filename))
	  stop 
	endif

      lfilename=lnblnk(filename)
      itypcalc=3
      write(6,"('type iprtlv')") 
      read(5,"(i1)") iprtlv
c    1 continue
      write(6,"('type phase')")
      read(5,"(a)") phase
      lphase=lnblnk(phase)

      write(6,"('type ishflag,ianisoflag')")
      read(5,*) ishflag,ianisoflag

      write(6,"('type epla,eplo,dep')")
      read(5,*) epla,eplo,depth


      write(6,"('delmin,delmax,delstep')")
      read(5,*) delmin,delmax,delstep



	stla=0.0
	stlo=delmin

    3 continue      
c      write(6,"('type stla,stlo')")
c      read(5,*) stla,stlo

	stlo=stlo+delstep


c
	geoco=0.993277
c      epla=atand(geoco*tand(epla))
c      stla=atand(geoco*tand(stla))


      call cagcrays(lu,filename,epla,eplo,depth,stla,stlo,
     #     phase,ishflag,ianisoflag,itypcalc,delta,azep,azst,
     #     maxarr,narr,ttime,slowness,dtddepth,ddeltadp,tstar,ellcorr,
     #     maxraypts,nraypts,iraytp,iraylay,
     #     rayrad,raydelta,raylat,raylong,raytime)
c
c      if(itypcalc.ge.0) then
c        write(6,"('delta,azep,azst:',3f10.4)") delta,azep,azst
c      endif
c
c      if(itypcalc.ge.1) then
c        write(6,"(a16,2x,i3,' arrivals at',f9.3,' degrees')") phase,narr,delta
c        do ia=1,narr
c          write(6,"('p,T,dDdp,tstar,ellcorr:',f9.5,f9.3,f9.3,f7.3,f7.3)") 
c     #          slowness(ia),
c     #          ttime(ia),ddeltadp(ia),tstar(ia),ellcorr(ia)
c        enddo
c      endif
c
      if(itypcalc.ge.2) then
        open(7,file='path')
        write(7,"('#PATH:',5f9.3)") epla,eplo,depth,stla,stlo
        write(7,"('#MODEL:',a)") filename(1:lfilename)
        write(7,"('#IFANI:',i2)") ianisoflag
        write(7,"('#IFSH:',i2)") ishflag
        write(7,"('#PHASE:',a)") phase(1:lphase)
        do iarr=1,narr
          write(7,"('#TTPARA:',5f12.5)") delta,slowness(iarr),ttime(iarr),
     #            ddeltadp(iarr),tstar(iarr)
          write(7,"('#START#')") 
          do i=1,nraypts(iarr)
            write(7,"(i5,i3,i4,5f10.3)") 
     #      i,iraytp(i,iarr),iraylay(i,iarr),rayrad(i,iarr),
     #      raydelta(i,iarr),raylat(i,iarr),raylong(i,iarr),raytime(i,iarr)


	write(98,*) raydelta(i,iarr),rayrad(i,iarr)-6371.
	write(99,*) rayrad(i,iarr)*sind(raydelta(i,iarr)-delta/2.0),
     #		    rayrad(i,iarr)*cosd(raydelta(i,iarr)-delta/2.0)
          enddo
          write(7,"('#END#')") 

	  write(100,*) delta,ttime(iarr)-raytime(nraypts(iarr),iarr)
	  write(101,*) slowness(iarr),ttime(iarr)-raytime(nraypts(iarr),iarr)
	  write(102,*) 1.d0/slowness(iarr)

        enddo
        close(7)
      endif

c      go to 3
	if(stlo.lt.delmax) goto3



	call system("date")
	write(6,"('------------------------')")

      end
