       SUBROUTINE conduction(ptimestep,pplay,pplev,pt,pdt,
     $                   tsurf,zzlev,zzlay,zdtconduc)

      IMPLICIT NONE

c=======================================================================
c
c   Molecular thermal conduction
c   
c   N. Descamp, F. Forget 05/1999
c
c=======================================================================

c-----------------------------------------------------------------------
c   declarations:
c-----------------------------------------------------------------------

#include "dimensions.h"
#include "dimphys.h"
#include "comcstfi.h"
#include "surfdat.h"
#include "chimiedata.h"
#include "conc.h"

c   arguments:
c   ----------

      REAL ptimestep
      REAL pplay(ngridmx,nlayermx)
      real pplev(ngridmx,nlayermx+1)
      REAL zzlay(ngridmx,nlayermx)
      real zzlev(ngridmx,nlayermx+1)
      REAL pt(ngridmx,nlayermx)
      real pdt(ngridmx,nlayermx)
      real tsurf(ngridmx)

      real zdtconduc(ngridmx,nlayermx)

c   local:
c   ------

      INTEGER i,ig,l, ngrid, nz
      real Akk,skk,phitop,m,tmean
      REAL alpha(nlayermx)
      real zt(nlayermx)
      REAL lambda(nlayermx)
      real muvol(nlayermx)
      REAL C(nlayermx)
      real D(nlayermx)
      real den(nlayermx)
      REAL pdtc(nlayermx)
      real zlay(nlayermx)
      real zlev(nlayermx+1)

c   constants used locally
c    ---------------------
c     The atmospheric conductivity is a function of temperature T :
c      conductivity = Akk* T**skk
      PARAMETER (skk=0.69)
      
      logical firstcall
      save firstcall
      data firstcall /.true./
c-----------------------------------------------------------------------
c   calcul des coefficients alpha et lambda
c-----------------------------------------------------------------------

      IF (firstcall) THEN
        write (*,*)'conduction: coeff to compute molecular',
     &             ' conductivity Akk,skk'
        write(*,*) Akk,skk
        firstcall = .false.
      END IF

! Initialize phitop
      phitop=0.0
      
      ngrid=ngridmx
      nz=nlayermx

      do ig=1,ngrid

        zt(1)=pt(ig,1)+pdt(ig,1)*ptimestep
c        zlay(1)=-log(pplay(ig,1)/pplev(ig,1))*Rnew(ig,1)*zt(1)/g
c        zlev(1)=0.0
        zlay(1)=zzlay(ig,1)
        zlev(1)=zzlev(ig,1)
      
        do i=2,nz

          zt(i)=pt(ig,i)+pdt(ig,i)*ptimestep 
c          tmean=zt(i)
c          if(zt(i).ne.zt(i-1))
c     &    tmean=(zt(i)-zt(i-1))/log(zt(i)/zt(i-1))
c          zlay(i)= zlay(i-1)
c     &          -log(pplay(ig,i)/pplay(ig,i-1))*Rnew(ig,i-1)*tmean/g
c          zlev(i)= zlev(i-1)
c     &         -log(pplev(ig,i)/pplev(ig,i-1))*Rnew(ig,i-1)*tmean/g
        zlay(i)=zzlay(ig,i)
        zlev(i)=zzlev(ig,i)
        enddo
        
c        zlev(nz+1)= zlev(nz)
c     &         -log(max(pplev(ig,nz+1),1.e-30)/pplev(ig,nz))
c     &           *Rnew(ig,nz)*tmean/g
c        if(pplev(ig,nz+1).eq.0.) 
c     &     zlev(nz+1)=zlev(nz)+(zlay(nz)-zlay(nz-1))
      
        zlev(nz+1)= zlev(nz)+10000.

        Akk=Akknew(ig,1) 
        lambda(1) = Akk*tsurf(ig)**skk/zlay(1)   

        DO i = 2 , nz
          Akk=Akknew(ig,i) 
          lambda(i)=Akk*zt(i)**skk/(zlay(i)-zlay(i-1)) 
        ENDDO
        DO i=1,nz-1
          muvol(i)=pplay(ig,i)/(rnew(ig,i)*zt(i)) 
          alpha(i)=cpnew(ig,i)*(muvol(i)/ptimestep)
     $                        *(zlev(i+1)-zlev(i))
        ENDDO

        muvol(nz)=pplay(ig,nz)/(rnew(ig,nz)*zt(nz)) 
        alpha(nz)=cpnew(ig,i)*(muvol(nz)/ptimestep)
     $                       *(zlev(nz+1)-zlev(nz))

c--------------------------------------------------------------------
c
c     calcul des coefficients C et D
c
c-------------------------------------------------------------------

      den(1)=alpha(1)+lambda(2)+lambda(1)
      C(1)=lambda(1)*(tsurf(ig)-zt(1))+lambda(2)*(zt(2)-zt(1))
      C(1)=C(1)/den(1)	     
      D(1)=lambda(2)/den(1)           
   
      DO i = 2,nz-1
        den(i)=alpha(i)+lambda(i+1)
        den(i)=den(i)+lambda(i)*(1-D(i-1))
         
        C(i) =lambda(i+1)*(zt(i+1)-zt(i)) 
     $       +lambda(i)*(zt(i-1)-zt(i)+C(i-1))    
        C(i) =C(i)/den(i)           

        D(i) =lambda(i+1) / den(i)
      ENDDO 

      den(nz)=alpha(nz) + lambda(nz) * (1-D(nz-1))
      C(nz)=C(nz-1)+zt(nz-1)-zt(nz) 
      C(nz)=(C(nz)*lambda(nz)+phitop) / den(nz) 
       	 	
c----------------------------------------------------------------------
c
c      calcul de la nouvelle temperature ptconduc
c
c----------------------------------------------------------------------

      DO i=1,nz
        pdtc(i)=0.
      ENDDO
      pdtc(nz)=C(nz)
      DO i=nz-1,1,-1
        pdtc(i)=C(i)+D(i)*pdtc(i+1)
      ENDDO 
c-----------------------------------------------------------------------
c
c     calcul de la tendance zdtconduc
c
c-----------------------------------------------------------------------
    
      DO i=1,nz
        zdtconduc(ig,i)=pdtc(i)/ptimestep
      ENDDO

      ENDDO             ! boucle sur ngrid
      RETURN
      END
