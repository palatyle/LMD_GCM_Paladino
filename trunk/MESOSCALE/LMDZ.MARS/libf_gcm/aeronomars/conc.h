c**********************************************************************

c	conc.h
	
c**********************************************************************

      integer ncomptot

      parameter (ncomptot=ncomp+1)

      common/conc/mmean,Akknew,cpnew,rnew
 
      real mmean(ngridmx,nlayermx)
      real Akknew(ngridmx,nlayermx)
      real cpnew(ngridmx,nlayermx)
      real rnew(ngridmx,nlayermx)


