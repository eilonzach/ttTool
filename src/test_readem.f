      implicit double precision (a-h,o-z)
      character*80 emodel
      parameter (maxlay=1000)
      double precision xb(maxlay)
      double precision xt(maxlay)
      dimension nlev(maxlay)
      dimension iflu(maxlay)
      dimension iani(maxlay)
      parameter (twopi=6.2831853072d0)
      logical isotropic
      common /plevel/iprtlv
c
      iprtlv=2
      write(6,"('type name of model to read')")
      read(5,"(a)") emodel
      open(1,file=emodel)
      call readem(1,maxlay,nlay,xb,xt,nlev,iflu,iani,ierr)
      close(1)
      do i=1,nlay
        do j=1,1+iani(i)
        isotropic=.false.
        if(j.eq.2) isotropic=.true.
        call evem(i,xb(i),isotropic,rho,vpv,vph,vsv,vsh,eta,qmu,qkappa,ierr)
        r8=xb(i)*twopi/360.d0
        write(6,"(2i3,f10.3,4f10.5)") i,iani(i),xb(i),
     #        r8/vpv,r8/vph,r8/vsv,r8/vsh
        call evem(i,xt(i),isotropic,rho,vpv,vph,vsv,vsh,eta,qmu,qkappa,ierr)
        r8=xt(i)*twopi/360.d0
        write(6,"(2i3,f10.3,4f10.5)") i,iani(i),xt(i),
     #        r8/vpv,r8/vph,r8/vsv,r8/vsh
        enddo
      enddo
      end
