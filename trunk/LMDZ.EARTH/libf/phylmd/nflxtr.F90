!
! $Id $
!
SUBROUTINE nflxtr(pdtime,pmfu,pmfd,pen_u,pde_u,pen_d,pde_d,pplay,paprs,x,dx) 
  USE dimphy
  IMPLICIT NONE 
!=====================================================================
! Objet : Melange convectif de traceurs a partir des flux de masse 
! Date : 13/12/1996 -- 13/01/97
! Auteur: O. Boucher (LOA) sur inspiration de Z. X. Li (LMD),
!         Brinkop et Sausen (1996) et Boucher et al. (1996).
! ATTENTION : meme si cette routine se veut la plus generale possible, 
!             elle a herite de certaines notations et conventions du 
!             schema de Tiedtke (1993).
! 1. En particulier, les couches sont numerotees de haut en bas !!!
!    Ceci est valable pour les flux
!    mais pas pour les entrees x, pplay, paprs !!!!
! 2. pmfu est positif, pmfd est negatif 
! 3. Tous les flux d'entrainements et de detrainements sont positifs 
!    contrairement au schema de Tiedtke d'ou les changements de signe!!!! 
!=====================================================================
!
  include "YOMCST.h"
  include "YOECUMF.h" 

  REAL,INTENT(IN) :: pdtime  ! pdtphys
!
! les flux sont definis au 1/2 niveaux 
! => pmfu(klev+1) et pmfd(klev+1) sont implicitement nuls
!
  REAL,DIMENSION(klon,klev),INTENT(IN) :: pmfu  ! flux de masse dans le panache montant 
  REAL,DIMENSION(klon,klev),INTENT(IN) :: pmfd  ! flux de masse dans le panache descendant
  REAL,DIMENSION(klon,klev),INTENT(IN) :: pen_u ! flux entraine dans le panache montant
  REAL,DIMENSION(klon,klev),INTENT(IN) :: pde_u ! flux detraine dans le panache montant
  REAL,DIMENSION(klon,klev),INTENT(IN) :: pen_d ! flux entraine dans le panache descendant
  REAL,DIMENSION(klon,klev),INTENT(IN) :: pde_d ! flux detraine dans le panache descendant

  REAL,DIMENSION(klon,klev),INTENT(IN)   :: pplay ! pression aux couches (bas en haut)
  REAL,DIMENSION(klon,klev+1),INTENT(IN) :: paprs ! pression aux 1/2 couches (bas en haut)
  REAL,DIMENSION(klon,klev),INTENT(IN)   :: x     ! q de traceur (bas en haut) 
  REAL,DIMENSION(klon,klev),INTENT(INOUT) :: dx   ! tendance de traceur  (bas en haut)

! flux convectifs mais en variables locales
  REAL,DIMENSION(klon,klev+1) :: zmfu  ! copie de pmfu avec klev+1 = 0
  REAL,DIMENSION(klon,klev+1) :: zmfd  ! copie de pmfd avec klev+1 = 0
  REAL,DIMENSION(klon,klev)   :: zen_u
  REAL,DIMENSION(klon,klev)   :: zde_u
  REAL,DIMENSION(klon,klev)   :: zen_d
  REAL,DIMENSION(klon,klev)   :: zde_d
  REAL                        :: zmfe

! variables locales      
! les flux de x sont definis aux 1/2 niveaux 
! xu et xd sont definis aux niveaux complets
  REAL,DIMENSION(klon,klev)   :: xu      ! q de traceurs dans le panache montant
  REAL,DIMENSION(klon,klev)   :: xd      ! q de traceurs dans le panache descendant
  REAL,DIMENSION(klon,klev+1) :: zmfux   ! flux de x dans le panache montant
  REAL,DIMENSION(klon,klev+1) :: zmfdx   ! flux de x dans le panache descendant
  REAL,DIMENSION(klon,klev+1) :: zmfex   ! flux de x dans l'environnement 
  INTEGER                     :: i, k 
  REAL,PARAMETER              :: zmfmin=1.E-10

! ==============================================
! Extension des flux UP et DN sur klev+1 niveaux
! ==============================================
  DO k=1,klev
     DO i=1,klon
        zmfu(i,k)=pmfu(i,k)
        zmfd(i,k)=pmfd(i,k)
     ENDDO
  ENDDO
  DO i=1,klon
     zmfu(i,klev+1)=0.
     zmfd(i,klev+1)=0.
  ENDDO
! ==========================================
! modif pour diagnostiquer les detrainements
! ==========================================
!   on privilegie l'ajustement de l'entrainement dans l'ascendance.

  DO k=1, klev
     DO i=1, klon
        zen_d(i,k)=pen_d(i,k)
        zde_u(i,k)=pde_u(i,k)
        zde_d(i,k) =-zmfd(i,k+1)+zmfd(i,k)+zen_d(i,k)
        zen_u(i,k) = zmfu(i,k+1)-zmfu(i,k)+zde_u(i,k)
     ENDDO
  ENDDO
! =========================================
! calcul des flux dans le panache montant
! =========================================
!
! Dans la premiere couche, on prend q comme valeur de qu

  DO i=1, klon
     zmfux(i,1)=0.0 
  ENDDO

! Autres couches
  DO k=1,klev
     DO i=1, klon
        IF ((zmfu(i,k+1)+zde_u(i,k)).lt.zmfmin) THEN
           xu(i,k)=x(i,k)
        ELSE
           xu(i,k)=(zmfux(i,k)+zen_u(i,k)*x(i,k))/(zmfu(i,k+1)+zde_u(i,k))
        ENDIF
        zmfux(i,k+1)=zmfu(i,k+1)*xu(i,k)
     ENDDO
  ENDDO
! ==========================================
! calcul des flux dans le panache descendant
! ==========================================
   
  DO i=1, klon
     zmfdx(i,klev+1)=0.0 
  ENDDO

  DO k=klev,1,-1
     DO i=1, klon
        IF ((zde_d(i,k)-zmfd(i,k)).lt.zmfmin) THEN
           xd(i,k)=x(i,k)
        ELSE
           xd(i,k)=(zmfdx(i,k+1)-zen_d(i,k)*x(i,k))/(zmfd(i,k)-zde_d(i,k))
        ENDIF
        zmfdx(i,k)=zmfd(i,k)*xd(i,k)
     ENDDO
  ENDDO
! ===================================================
! introduction du flux de retour dans l'environnement
! ===================================================

  DO k=2, klev
     DO i=1, klon
        zmfe=-zmfu(i,k)-zmfd(i,k)
        IF (zmfe.le.0.) then
           zmfex(i,k)= zmfe*x(i,k)
        ELSE
           zmfex(i,k)= zmfe*x(i,k-1)
        ENDIF
     ENDDO
  ENDDO

  DO i=1, klon
     zmfex(i,1)=0.
     zmfex(i,klev+1)=0.
  ENDDO
! ==========================
! calcul final des tendances
! ==========================
  DO k=1, klev
     DO i=1, klon
        dx(i,k)=RG/(paprs(i,k)-paprs(i,k+1))*pdtime*  &
             ( zmfux(i,k) - zmfux(i,k+1) +            &
             zmfdx(i,k) - zmfdx(i,k+1) +              &
             zmfex(i,k) - zmfex(i,k+1) )
     ENDDO
  ENDDO
  
END SUBROUTINE nflxtr
