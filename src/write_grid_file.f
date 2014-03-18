      character*80 outfile
      
      outfile='us_grid'
      
      latmin=25.00
      latmax=50.00
      lonmin=-125.00
      lonmax=-70.00
      
      nlat=latmax-latmin+1.0
      nlon=lonmax-lonmin+1.0
      
      write(6,*) nlat
      write(6,*) nlon
      
      open(1,file=outfile)
      do ilat=1,nlat
      xlat=latmin+(ilat-1)
      	do ilon=1,nlon	
	xlon=lonmin+(ilon-1)
	write(1,"(1f8.4,3x,1f9.4)") xlat,xlon
	enddo
      enddo
      
      end
