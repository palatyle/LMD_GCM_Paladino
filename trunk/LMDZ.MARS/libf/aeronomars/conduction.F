       SUBROUTINE conduction(ngrid,nlayer,ptimestep,pplay,pplev,pt,pdt,
     $                   tsurf,zzlev,zzlay,zdtconduc)

      use conc_mod, only: Akknew, rnew, cpnew
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

c   arguments:
c   ----------

      integer,intent(in) :: ngrid ! number of atmospheric columns
      integer,intent(in) :: nlayer ! number of atmospheric layers
      real,intent(in) :: ptimestep
      REAL,intent(in) :: pplay(ngrid,nlayer)
      real,intent(in) :: pplev(ngrid,nlayer+1)
      REAL,intent(in) :: zzlay(ngrid,nlayer)
      real,intent(in) :: zzlev(ngrid,nlayer+1)
      REAL,intent(in) :: pt(ngrid,nlayer)
      real,intent(in) :: pdt(ngrid,nlayer)
      real,intent(in) :: tsurf(ngrid)

      real,intent(out) :: zdtconduc(ngrid,nlayer)

c   local:
c   ------

      INTEGER i,ig,l
      real Akk
      real,save :: phitop
      real m,tmean
      REAL alpha(nlayer)
      real zt(nlayer)
      REAL lambda(nlayer)
      real muvol(nlayer)
      REAL C(nlayer)
      real D(nlayer)
      real den(nlayer)
      REAL pdtc(nlayer)
      real zlay(nlayer)
      real zlev(nlayer+1)

c   constants used locally
c    ---------------------
c     The atmospheric conductivity is a function of temperature T :
c      conductivity = Akk* T**skk
      REAL,PARAMETER :: skk=0.69
      
      logical,save :: firstcall=.true.

c-----------------------------------------------------------------------
c   calcul des coefficients alpha et lambda
c-----------------------------------------------------------------------

      IF (firstcall) THEN
!        write (*,*)'conduction: coeff to compute molecular',
!     &             ' conductivity Akk,skk'
!        write(*,*) Akk,skk
! NB: Akk is undefined at this stage
        write (*,*)'conduction: coeff to compute molecular',
     &             ' conductivity skk = ', skk

! Initialize phitop
        phitop=0.0
        
        firstcall = .false.
      ENDIF ! of IF (firstcall)

      do ig=1,ngrid

        zt(1)=pt(ig,1)+pdt(ig,1)*ptimestep
c        zlay(1)=-log(pplay(ig,1)/pplev(ig,1))*Rnew(ig,1)*zt(1)/g
c        zlev(1)=0.0
        zlay(1)=zzlay(ig,1)
        zlev(1)=zzlev(ig,1)
      
        do i=2,nlayer

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
        
c        zlev(nlayer+1)= zlev(nlayer)
c     &         -log(max(pplev(ig,nlayer+1),1.e-30)/pplev(ig,nlayer))
c     &           *Rnew(ig,nlayer)*tmean/g
c        if(pplev(ig,nlayer+1).eq.0.) 
c     &     zlev(nlayer+1)=zlev(nlayer)+(zlay(nlayer)-zlay(nlayer-1))
      
        zlev(nlayer+1)= zlev(nlayer)+10000.

        Akk=Akknew(ig,1) 
        lambda(1) = Akk*tsurf(ig)**skk/zlay(1)   

        DO i = 2 , nlayer
          Akk=Akknew(ig,i) 
          lambda(i)=Akk*zt(i)**skk/(zlay(i)-zlay(i-1)) 
        ENDDO
        DO i=1,nlayer-1
          muvol(i)=pplay(ig,i)/(rnew(ig,i)*zt(i)) 
          alpha(i)=cpnew(ig,i)*(muvol(i)/ptimestep)
     $                        *(zlev(i+1)-zlev(i))
        ENDDO

        muvol(nlayer)=pplay(ig,nlayer)/(rnew(ig,nlayer)*zt(nlayer)) 
        alpha(nlayer)=cpnew(ig,i)*(muvol(nlayer)/ptimestep)
     $                       *(zlev(nlayer+1)-zlev(nlayer))

c--------------------------------------------------------------------
c
c     calcul des coefficients C et D
c
c-------------------------------------------------------------------

        den(1)=alpha(1)+lambda(2)+lambda(1)
        C(1)=lambda(1)*(tsurf(ig)-zt(1))+lambda(2)*(zt(2)-zt(1))
        C(1)=C(1)/den(1)	     
        D(1)=lambda(2)/den(1)           
   
        DO i = 2,nlayer-1
          den(i)=alpha(i)+lambda(i+1)
          den(i)=den(i)+lambda(i)*(1-D(i-1))
           
          C(i) =lambda(i+1)*(zt(i+1)-zt(i)) 
     $         +lambda(i)*(zt(i-1)-zt(i)+C(i-1))    
          C(i) =C(i)/den(i)           

          D(i) =lambda(i+1) / den(i)
        ENDDO 

        den(nlayer)=alpha(nlayer) + lambda(nlayer) * (1-D(nlayer-1))
        C(nlayer)=C(nlayer-1)+zt(nlayer-1)-zt(nlayer) 
        C(nlayer)=(C(nlayer)*lambda(nlayer)+phitop) / den(nlayer) 
       	 	
c----------------------------------------------------------------------
c
c      calcul de la nouvelle temperature ptconduc
c
c----------------------------------------------------------------------

        DO i=1,nlayer
          pdtc(i)=0.
        ENDDO
        pdtc(nlayer)=C(nlayer)
        DO i=nlayer-1,1,-1
          pdtc(i)=C(i)+D(i)*pdtc(i+1)
        ENDDO 
c-----------------------------------------------------------------------
c
c     calcul de la tendance zdtconduc
c
c-----------------------------------------------------------------------
    
        DO i=1,nlayer
          zdtconduc(ig,i)=pdtc(i)/ptimestep
        ENDDO

      enddo ! of do ig=1,ngrid

      RETURN
      END
