       SUBROUTINE conduction(nlon, nlev,ptimestep,pplay,pplev,pt,
     $                   tsurf,zzlev,zzlay,d_t_conduc)
    
      use dimphy
      use conc, only: akknew, rnew, cpnew
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

      integer,intent(in) :: nlon ! number of atmospheric columns
      integer,intent(in) :: nlev ! number of atmospheric layers
      real,intent(in) :: ptimestep
      REAL,intent(in) :: pplay(nlon,nlev)  ! pressure at middle of layers (Pa)
      real,intent(in) :: pplev(nlon,nlev+1)
      REAL,intent(in) :: zzlay(nlon,nlev)   ! (m)
      real,intent(in) :: zzlev(nlon,nlev+1)
      REAL,intent(in) :: pt(nlon,nlev)
      real,intent(in) :: tsurf(nlon)

      real,intent(out) :: d_t_conduc(nlon,nlev)

c   local:
c   ------

      INTEGER i,ig,l
      real Akk
      real,save :: phitop
      real m,tmean
      REAL alpha(nlev)
      real zt(nlev)
      REAL lambda(nlev)
      real muvol(nlev)   ! kg m-3
      REAL C(nlev)
      real D(nlev)
      real den(nlev)
      REAL pdtc(nlev)
      real zlay(nlev)
      real zlev(nlev+1)

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

      do ig=1,nlon

c        zt(1)=pt(ig,1)+pdt(ig,1)*ptimestep
         zt(1)=pt(ig,1)

c        zlay(1)=-log(pplay(ig,1)/pplev(ig,1))*Rnew(ig,1)*zt(1)/g
c        zlev(1)=0.0
         zlay(1)=zzlay(ig,1)
         zlev(1)=zzlev(ig,1)
      
        do i=2,nlev

           zt(i)= pt(ig,i)
c           print*, zt(i)

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
        zlev(nlev+1)= zzlev(ig,nlev+1)
      
        Akk=akknew(ig,1) 
        lambda(1) = Akk*tsurf(ig)**skk/zlay(1)   

        DO i = 2 , nlev
          Akk=akknew(ig,i) 
          lambda(i)=Akk*zt(i)**skk/(zlay(i)-zlay(i-1)) 
        ENDDO
        DO i=1,nlev-1
c          print*, rnew(1,i)
           muvol(i)=pplay(ig,i)/(rnew(ig,i)*zt(i)) 
           alpha(i)=cpnew(ig,i)*(muvol(i)/ptimestep)
     $                        *(zlev(i+1)-zlev(i))
       ENDDO


c           if (ig .eq. 2) then
c              print*, '---conduction---'
c              print*, i, cpnew(ig,i), zt(i)
c           endif
C        stop

        muvol(nlev)=pplay(ig,nlev)/(rnew(ig,nlev)*zt(nlev)) 
        alpha(nlev)=cpnew(ig,i)*(muvol(nlev)/ptimestep)
     $                       *(zlev(nlev+1)-zlev(nlev))

c--------------------------------------------------------------------
c
c     calcul des coefficients C et D
c
c-------------------------------------------------------------------

        den(1)=alpha(1)+lambda(2)+lambda(1)
        C(1)=lambda(1)*(tsurf(ig)-zt(1))+lambda(2)*(zt(2)-zt(1))
        C(1)=C(1)/den(1)	     
        D(1)=lambda(2)/den(1)           
   
        DO i = 2,nlev-1
          den(i)=alpha(i)+lambda(i+1)
          den(i)=den(i)+lambda(i)*(1-D(i-1))
           
          C(i) =lambda(i+1)*(zt(i+1)-zt(i)) 
     $         +lambda(i)*(zt(i-1)-zt(i)+C(i-1))    
          C(i) =C(i)/den(i)           

          D(i) =lambda(i+1) / den(i)
        ENDDO 

        den(nlev)=alpha(nlev) + lambda(nlev) * (1-D(nlev-1))
        C(nlev)=C(nlev-1)+zt(nlev-1)-zt(nlev) 
        C(nlev)=(C(nlev)*lambda(nlev)+phitop) / den(nlev) 
       	 	
c----------------------------------------------------------------------
c
c      calcul de la nouvelle temperature ptconduc
c
c----------------------------------------------------------------------

        DO i=1,nlev
          pdtc(i)=0.
        ENDDO
        pdtc(nlev)=C(nlev)
        DO i=nlev-1,1,-1
          pdtc(i)=C(i)+D(i)*pdtc(i+1)
        ENDDO 
c-----------------------------------------------------------------------
c
c     calcul de la tendance zdtconduc
c
c-----------------------------------------------------------------------
    
        DO i=1,nlev
          d_t_conduc(ig,i)=pdtc(i)/ptimestep
c          print*, i, zdtconduc(0, i)

        ENDDO
c
      enddo ! of do ig=1,nlon

      RETURN
      END
