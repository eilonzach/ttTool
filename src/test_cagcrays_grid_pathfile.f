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
c
      character*80 systemstring
      
      parameter(nstations=5000)
      character*80 station_file
      character*8 a,station_name
      character*8 station(nstations)
      dimension slat(nstations)
      dimension slon(nstations)
      character*30 pathfilename
c
c      common/plevel/iprtlv
c
c-------------------------------------------------------------------------------
c
      lu=1

c      write(6,"('type name of the .RAYS file: ')")
c      read(5,"(a)") filename

      filename='NOH2O.RAYS'
      station_file='us_grid'
c
c-----read in station file
c      
      open(1,file=station_file)
      ios=0
      ns=0
      do while (ios.eq.0)
      	read(1,"(1f8.4,3x,1f9.4)",iostat=ios) b,c
	
	if (ios.eq.0) then
	ns=ns+1
c	station(ns)=a
	slat(ns)=b
	slon(ns)=c
	
	write(6,"(1f8.4,1x,1f9.4)") slat(ns),slon(ns)
	endif      
      enddo
      close(1)
      write(6,*) ns
      numstation=0

      lfilename=lnblnk(filename)
      itypcalc=3
      
c      systemstring='mkdir PATH_FILES'
c      ierror=system(systemstring)
      
c      write(6,"('type iprtlv')") 
c      read(5,"(i1)") iprtlv
      iprtlv=0
    1 continue
      write(6,"('type phase')")
      read(5,"(a)") phase
      lphase=lnblnk(phase)
      write(6,"('type epla,eplo,dep')")
      read(5,*) epla,eplo,depth
      open(7,file='path_grid')
      write(7,"('#NSTA: ',i5)") ns
    3 continue
      numstation=numstation+1
          
c      write(6,"('type stla,stlo')")
c      read(5,*) stla,stlo
      stla=slat(numstation)
      stlo=slon(numstation)
      write(6,*) stla,stlo
c      write(6,"('type ishflag,ianisoflag')")
c      read(5,*) ishflag,ianisoflag
      ishflag=0
      ianisoflag=0
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
      if(itypcalc.ge.0) then
        write(6,"('delta,azep,azst:',3f10.4)") delta,azep,azst
      endif
c
      if(itypcalc.ge.1) then
        write(6,"(a16,2x,i3,' arrivals at',f9.3,' degrees')") phase,narr,delta
        do ia=1,narr
          write(6,"('p,T,dDdp,tstar,ellcorr:',f9.5,f9.3,f9.3,f7.3,f7.3)") 
     #          slowness(ia),
     #          ttime(ia),ddeltadp(ia),tstar(ia),ellcorr(ia)
        enddo  
      endif
c
c      station_name=station(numstation)
c      station_name=trim(station_name)
c      pathfilename='PATH_FILES/path_'//station_name
c      write(6,*) station_name
c      write(6,*) pathfilename
c
      if(itypcalc.ge.2) then
        write(7,"('#PATH:',5f9.3)") epla,eplo,depth,stla,stlo
        write(7,"('#MODEL:',a)") filename(1:lfilename)
        write(7,"('#IFANI:',i2)") ianisoflag
        write(7,"('#IFSH:',i2)") ishflag
        write(7,"('#PHASE:',a)") phase(1:lphase)
	
	if (narr.eq.0) then
	write(7,"('#TTPARA: 0 0 0 0 0 0')")
	else
        do iarr=1,narr
          write(7,"('#TTPARA:',1i5,2x,5f12.5)") nraypts(iarr),delta,
     #			slowness(iarr),ttime(iarr),
     #            	ddeltadp(iarr),tstar(iarr)
c          write(7,"('#START#')") 
          do i=1,nraypts(iarr)
            write(7,"(i5,i3,i4,5f10.3)") 
     #      i,iraytp(i,iarr),iraylay(i,iarr),rayrad(i,iarr),
     #      raydelta(i,iarr),raylat(i,iarr),raylong(i,iarr),raytime(i,iarr)


	write(99,*) raydelta(i,iarr),rayrad(i,iarr)-6371.
c	write(99,*) rayrad(i,iarr)*sind(raydelta(i,iarr)-delta/2.0),
c     #		    rayrad(i,iarr)*cosd(raydelta(i,iarr)-delta/2.0)

          enddo
c          write(7,"('#END#')") 
        enddo
	endif
      endif
      if (numstation.lt.ns) then
      go to 3
      endif
      close(7)
      end
