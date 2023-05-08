       SUBROUTINE molvis(ngrid,nlayer,ptimestep,
     &            pplay,pplev,pt,pdteuv,pdtconduc
     $           ,pvel,tsurf,zzlev,zzlay,zdvelmolvis)
      
      use conc_mod, only: cpnew, Akknew, rnew
      IMPLICIT NONE

c=======================================================================
c
c   Molecular Viscosity Diffusion
c  
c   Based on conduction.F  by N. Descamp, F. Forget 05/1999
c
c   modified by M. Angelats i Coll
c
c=======================================================================

c-----------------------------------------------------------------------
c   declarations:
c-----------------------------------------------------------------------

c   arguments:
c   ----------

      integer,intent(in) :: ngrid ! number of atmospheric columns
      integer,intent(in) :: nlayer ! number of atmospheric layers
      REAL ptimestep
      REAL pplay(ngrid,nlayer)
      REAL pplev(ngrid,nlayer+1)
      REAL zzlay(ngrid,nlayer)
      REAL zzlev(ngrid,nlayer+1)
      real pt(ngrid,nlayer)
      real tsurf(ngrid)
      REAL pvel(ngrid,nlayer)
      REAL pdvel(ngrid,nlayer)
      real pdteuv(ngrid,nlayer)
      real pdtconduc(ngrid,nlayer)

      real zdvelmolvis(ngrid,nlayer)

c   local:
c   ------

      INTEGER l,ig, nz
      real Akk,phitop,fac, m, tmean
      REAL zvel(nlayer)
      real zt(nlayer)
      REAL alpha(nlayer)
      REAL lambda(nlayer)
      real muvol(nlayer)
      REAL C(nlayer)
      real D(nlayer)
      real den(nlayer)
      REAL pdvelm(nlayer)
      REAL zlay(nlayer)
      real zlev(nlayer+1)

c   constants used locally
c    ---------------------
c     The atmospheric conductivity is a function of temperature T :
c      conductivity = Akk* T**skk
c     Molecular viscosity is related to thermal conductivity by:
c       conduc = 0.25*(9*gamma - 5)* Cv * molvis
c     where gamma = Cp/Cv.  For dry air.


      REAL,PARAMETER :: skk=0.69

      REAL,PARAMETER :: velsurf =0.0
      
      logical,save :: firstcall=.true.

c-----------------------------------------------------------------------
c   calcul des coefficients alpha et lambda
c-----------------------------------------------------------------------

      IF (firstcall) THEN
!        write(*,*)'molvis: coeff of molecular viscosity Akk,skk,factor'
!        write(*,*) Akk,skk,fac
! NB: Akk and fac are undefined at firstcall
        write(*,*)'molvis: coeff of molecular viscosity skk ', skk
        
        firstcall = .false.
      END IF

! Initialize phitop
      phitop=0.0

      nz=nlayer

      do ig=1,ngrid

        zt(1)=pt(ig,1)+(pdteuv(ig,1)+pdtconduc(ig,1))*ptimestep
        zvel(1)=pvel(ig,1)
c        zlay(1)=-log(pplay(ig,1)/pplev(ig,1))*Rnew(ig,1)*zt(1)/g
c        zlev(1)=0.0

        zlay(1)=zzlay(ig,1)
        zlev(1)=zzlev(ig,1)

        do l=2,nz
          zt(l)=pt(ig,l)+(pdteuv(ig,l)+pdtconduc(ig,l))*ptimestep
          zvel(l)=pvel(ig,l)
c          tmean=zt(l)
c          if(zt(l).ne.zt(l-1)) tmean=(zt(l)-zt(l-1))/log(zt(l)/zt(l-1))
c          zlay(l)= zlay(l-1)
c     &          -log(pplay(ig,l)/pplay(ig,l-1))*Rnew(ig,l-1)*tmean/g
c          zlev(l)= zlev(l-1)
c     &         -log(pplev(ig,l)/pplev(ig,l-1))*Rnew(ig,l-1)*tmean/g
        zlay(l)=zzlay(ig,l)
        zlev(l)=zzlev(ig,l)
        enddo

c          zlev(nz+1)= zlev(nz)
c     &         -log(max(pplev(ig,nz+1),1.e-30)/pplev(ig,nz))
c     &          *Rnew(ig,nz)*tmean/g
c          if(pplev(ig,nz+1).eq.0.)
c     &       zlev(nz+1)=zlev(nz)+(zlay(nz)-zlay(nz-1))

          zlev(nz+1)= zlev(nz)+10000.

        fac=0.25*(9.*cpnew(ig,1)-5.*(cpnew(ig,1)-rnew(ig,1))) 
        Akk=Akknew(ig,1)
        lambda(1)=Akk*tsurf(ig)**skk/zlay(1)/fac   
c        write(*,*) 'rnew(ig,nz)  ',ig , rnew(ig,nz)

        DO l=2,nz
          fac=(9.*cpnew(ig,l)-5.*(cpnew(ig,l)-rnew(ig,l)))/4. 
          Akk=Akknew(ig,l)
          lambda(l)=Akk/fac*zt(l)**skk/(zlay(l)-zlay(l-1)) 
        ENDDO
    
        DO l=1,nz-1
          muvol(l)=pplay(ig,l)/(rnew(ig,l)*zt(l)) 
          alpha(l)=(muvol(l)/ptimestep)*(zlev(l+1)-zlev(l))
        ENDDO
        muvol(nz)=pplay(ig,nz)/(rnew(ig,nz)*zt(nz)) 
        alpha(nz)=(muvol(nz)/ptimestep)*(zlev(nz+1)-zlev(nz))

c--------------------------------------------------------------------
c
c     calcul des coefficients C et D
c
c-------------------------------------------------------------------

      den(1)=alpha(1)+lambda(2)+lambda(1)
      C(1)=lambda(1)*(velsurf-zvel(1))+lambda(2)*(zvel(2)-zvel(1))
      C(1)=C(1)/den(1)	     
      D(1)=lambda(2)/den(1)           
   
      DO l = 2,nz-1
         den(l)=alpha(l)+lambda(l+1)
         den(l)=den(l)+lambda(l)*(1-D(l-1))
         
         C(l) =lambda(l+1)*(zvel(l+1)-zvel(l)) 
     $        +lambda(l)*(zvel(l-1)-zvel(l)+C(l-1))    
         C(l) =C(l)/den(l)           

         D(l) =lambda(l+1) / den(l)
      ENDDO 

      den(nz)=alpha(nz) + lambda(nz) * (1-D(nz-1))
      C(nz)=C(nz-1)+zvel(nz-1)-zvel(nz) 
      C(nz)=(C(nz)*lambda(nz)+phitop) / den(nz) 
       	 	
c----------------------------------------------------------------------
c
c      calcul de la nouvelle pdvelm 
c
c----------------------------------------------------------------------

      DO l=1,nz
         pdvelm(l)=0.
      ENDDO
         pdvelm(nz)=C(nz)
      DO l=nz-1,1,-1
         pdvelm(l)=C(l)+D(l)*pdvelm(l+1)
      ENDDO 
c-----------------------------------------------------------------------
c
c     calcul de la tendance zdvelmolvis
c
c-----------------------------------------------------------------------
    
      DO l=1,nz
        zdvelmolvis(ig,l)=pdvelm(l)/ptimestep
      ENDDO

      ENDDO             ! boucle sur ngrid

      RETURN
      END
