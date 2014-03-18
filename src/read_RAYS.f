      real*4 rr
      real*8 dd
      character*80 c80
      integer*4 ii
      integer*4 ipos
c
      open(1,file='NOH2O.RAYS',access='stream',form='unformatted')
c
      ipos=1
      read(1,pos=ipos) ii
      iversion=ii
      ipos=ipos+4
      write(6,"('version:',i6)") iversion
c
      read(1,pos=ipos) c80
      ipos=ipos+80
      write(6,"('title:',a80)") c80
c
      read(1,pos=ipos) ii
      numemlev=ii
      ipos=ipos+4
      write(6,"('number of earth model levels:',i6)") numemlev
c
c---- skip 14 variables describing the earth model at each level
c---- (double precision => 8 bytes per number)
      iskip=numemlev*8*14
      ipos=ipos+iskip
c
c---- skip 3*13 variable describing the cubic polynomials
c---- representing the earth model at any depth
c---- (double precision => 8 bytes per number)
      iskip=numemlev*8*13*3
      ipos=ipos+iskip
c
      read(1,pos=ipos) ii
      numemreg=ii
      ipos=ipos+4
      write(6,"('number of earth model regions:',i5)") numemreg
c
c---- skip 4 variables (2 double precision, 2 integer*4) describing
c---- each earth model region
c
      iskip=8*numemreg*2+4*numemreg*2
c
      ipos=ipos+iskip
c
      read(1,pos=ipos) ii
      np=ii
      ipos=ipos+4
      write(6,"('number of p values in table: ',i6)") np
c
c---- skip the array of p values
c
      iskip=np*4
      ipos=ipos+iskip
c
      read(1,pos=ipos) ii
      nlaytabs=ii
      ipos=ipos+4
      write(6,"('number of layers in table:',i4)") nlaytabs
c
      do ilay=1,nlaytabs
	read(1,pos=ipos) ii
	ilayer=ii
	ipos=ipos+4
	read(1,pos=ipos) ii
	iani=ii
	ipos=ipos+4
	read(1,pos=ipos) ii
	iraytype=ii
	ipos=ipos+4
	read(1,pos=ipos) ii
	isub=ii
	ipos=ipos+4
	read(1,pos=ipos) rr
	top4=rr
	ipos=ipos+4
	read(1,pos=ipos) rr
	bot4=rr
	ipos=ipos+4
	read(1,pos=ipos) ii
	npout=ii
	ipos=ipos+4
	write(6,100) ilayer,iani,iraytype,isub,top4,bot4,npout,np
  100   format(i4,i4,i4,i4,f10.3,f10.3,i6,i6)
c
c---- skip 6 arrays of np elements describing the variables in each 
c---- layer for each p: 
c---- time, delta, ddelta-dp, tstar, turning radius, unused
c
	iskip=6*np*4
	ipos=ipos+iskip
      enddo
      close(1)
      end
