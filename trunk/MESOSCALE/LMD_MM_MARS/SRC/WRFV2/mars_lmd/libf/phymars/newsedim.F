      SUBROUTINE newsedim(ngrid,nlay,nsize,ptimestep,
     &  pplev,masse,epaisseur,pt,rd,rho,pqi,wq)
      IMPLICIT NONE

c=======================================================================
c
c      Compute sedimentation of 1 tracer 
c      of radius rd (m) and density rho (kg.m-3) 
c
c=======================================================================

c-----------------------------------------------------------------------
c   declarations:
c   -------------

#include "dimensions.h"
#include "dimphys.h"
#include "comcstfi.h"
c
c   arguments:
c   ----------

      INTEGER ngrid,nlay,nsize
      REAL ptimestep            ! pas de temps physique (s)
      REAL pplev(ngrid,nlay+1)   ! pression aux inter-couches (Pa)
      REAL pt(ngrid,nlay)    ! temperature au centre des couches (K)
      real masse (ngrid,nlay) ! masse d'une couche (kg)
      real epaisseur (ngrid,nlay)  ! epaisseur d'une couche (m)
      real rd(nsize)             ! particle radius (m)
      real rho             ! particle density (kg.m-3)


c    Traceurs :
      real pqi(ngrid,nlay)  ! traceur   (e.g. ?/kg)
      real wq(ngridmx,nlay+1)  ! flux de traceur durant timestep (?/m-2)
      
c   local:
c   ------

      INTEGER l,ig, k, i
      REAL rfall

      LOGICAL firstcall
      SAVE firstcall

c    Traceurs :
c    ~~~~~~~~ 
      real traversee (ngridmx,nlayermx)
      real vstokes(ngridmx,nlayermx)
      real w(ngridmx,nlayermx)
      real ptop, dztop, Ep, Stra

      DATA firstcall/.true./

c    Physical constant
c    ~~~~~~~~~~~~~~~~~
      REAL visc, molrad
c     Gas molecular viscosity (N.s.m-2)
      data visc/1.e-5/       ! CO2
c     Effective gas molecular radius (m)
      data molrad/2.2e-10/   ! CO2

c     local and saved variable
      real a,b
      save a,b



c    ** un petit test de coherence
c       --------------------------

      IF (firstcall) THEN
         IF(ngrid.NE.ngridmx) THEN
            PRINT*,'STOP dans newsedim'
            PRINT*,'probleme de dimensions :'
            PRINT*,'ngrid  =',ngrid
            PRINT*,'ngridmx  =',ngridmx
            STOP
         ENDIF
         firstcall=.false.


c       Preliminary calculations for sedimenation velocity :
c       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c       - Constant to compute stokes speed simple formulae
c        (Vstokes =  b * rho* r**2   avec   b= (2/9) * rho * g / visc
         b = 2./9. * g / visc

c       - Constant  to compute gas mean free path
c        l= (T/P)*a, avec a = (  0.707*8.31/(4*pi*molrad**2 * avogadro))
         a = 0.707*8.31/(4*3.1416* molrad**2  * 6.023e23)

c       - Correction to account for non-spherical shape (Murphy et al.  1990)
c   (correction = 0.85 for irregular particles, 0.5 for disk shaped particles)
c        a = a    * 0.85
      ENDIF





c-----------------------------------------------------------------------
c    1. initialisation
c    -----------------

c     Sedimentation velocity (m/s)
c     ~~~~~~~~~~~~~~~~~~~~~~
c     (stokes law corrected for low pressure by the Cunningham
c     slip-flow correction  according to Rossow (Icarus 36, 1-50, 1978)

        do  l=1,nlay
          do ig=1, ngrid
            if (nsize.eq.1) then 
              rfall=rd(1)
            else
              i=ngrid*(l-1)+ig
              rfall=rd(i)
            endif  
            vstokes(ig,l) = b * rho * rfall**2 *
     &      (1 + 1.333* ( a*pt(ig,l)/pplev(ig,l) )/rfall)

c           Layer crossing time (s) :
            traversee(ig,l)= epaisseur(ig,l)/vstokes(ig,l)
          end do
        end do


c     Calcul de la masse d'atmosphere correspondant a q transferee
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     (e.g. on recherche le niveau  en dessous de laquelle le traceur
c      va traverser le niveau intercouche l : "dztop" est sa hauteur
c      au dessus de l (m), "ptop" est sa pression (Pa))

      do  l=1,nlay
        do ig=1, ngrid
             
             dztop = vstokes(ig,l)*  ptimestep 
             Ep=0
             k=0

c **************************************************************
c            Simple Method
             w(ig,l) =
     &       (1- exp(-dztop*g/(r*pt(ig,l))))*pplev(ig,l) / g
cc           write(*,*) 'OK simple method l,w =', l, w(ig,l)
cc           write(*,*) 'OK simple method dztop =', dztop
c **************************************************************
cccc         Complex method :
            if (dztop.gt.epaisseur(ig,l)) then
cccc            Cas ou on "epuise" la couche l : On calcule le flux
cccc            Venant de dessus en tenant compte de la variation de Vstokes

               Ep= epaisseur(ig,l)
               Stra= traversee(ig,l)
               do while(dztop.gt.Ep.and.l+k+1.le.nlay)
                 k=k+1
                 dztop= Ep + vstokes(ig,l+k)*(ptimestep -Stra)
                 Ep = Ep + epaisseur(ig,l+k)
                 Stra = Stra + traversee(ig,l+k)
               enddo 
               Ep = Ep - epaisseur(ig,l+k)
             end if
             ptop=pplev(ig,l+k)*exp(-(dztop-Ep)*g/(r*pt(ig,l+k)))
             w(ig,l) = (pplev(ig,l) -Ptop)/g
c
cc           write(*,*) 'OK new    method l,w =', l, w(ig,l)
cc           write(*,*) 'OK new    method dztop =', dztop
cc       if(l.eq.7)write(*,*)'l=7,k,pplev,Ptop',pplev(ig,l),Ptop
cc       if(l.eq.7)write(*,*)'l=7,dztop,Ep',dztop,Ep
cc            if(l.eq.6)write(*,*)'l=6,k, w',k, w(1,l)
cc            if(l.eq.7)write(*,*)'l=7,k, w',k, w(1,l)
cc            if(l.eq.8)write(*,*)'l=8,k, w',k, w(1,l)
c **************************************************************
        end do
      end do

      call vlz_fi(ngrid,pqi,2.,masse,w,wq)
c         write(*,*) ' newsed: wq(6), wq(7), q(6)',
c    &                wq(1,6),wq(1,7),pqi(1,6)


      RETURN
      END

