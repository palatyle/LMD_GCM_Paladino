!
! $Id $
!
SUBROUTINE cvltr(pdtime,da, phi, mp,paprs,pplay,x,upd,dnd,dx)
  USE dimphy
  IMPLICIT NONE 
!=====================================================================
! Objet : convection des traceurs / KE
! Auteurs: M-A Filiberti and J-Y Grandpeix
!=====================================================================

  include "YOMCST.h"
  include "YOECUMF.h" 

! Entree
  REAL,INTENT(IN)                           :: pdtime
  REAL,DIMENSION(klon,klev),INTENT(IN)      :: da
  REAL,DIMENSION(klon,klev,klev),INTENT(IN) :: phi
  REAL,DIMENSION(klon,klev),INTENT(IN)      :: mp
  REAL,DIMENSION(klon,klev+1),INTENT(IN)    :: paprs ! pression aux 1/2 couches (bas en haut)
  REAL,DIMENSION(klon,klev),INTENT(IN)      :: pplay ! pression pour le milieu de chaque couche
  REAL,DIMENSION(klon,klev),INTENT(IN)      :: x     ! q de traceur (bas en haut) 
  REAL,DIMENSION(klon,klev),INTENT(IN)      :: upd   ! saturated updraft mass flux
  REAL,DIMENSION(klon,klev),INTENT(IN)      :: dnd   ! saturated downdraft mass flux

! Sortie
  REAL,DIMENSION(klon,klev),INTENT(OUT) :: dx ! tendance de traceur  (bas en haut)

! Variables locales     
! REAL,DIMENSION(klon,klev)       :: zed
  REAL,DIMENSION(klon,klev,klev)  :: zmd
  REAL,DIMENSION(klon,klev,klev)  :: za
  REAL,DIMENSION(klon,klev)       :: zmfd,zmfa
  REAL,DIMENSION(klon,klev)       :: zmfp,zmfu
  INTEGER                         :: i,k,j 
  REAL                            :: pdtimeRG

! =========================================
! calcul des tendances liees au downdraft
! =========================================
!cdir collapse
  DO j=1,klev
  DO i=1,klon
!   zed(i,j)=0.
    zmfd(i,j)=0.
    zmfa(i,j)=0.
    zmfu(i,j)=0.
    zmfp(i,j)=0.
  END DO
  END DO
!cdir collapse
  DO k=1,klev
  DO j=1,klev
  DO i=1,klon
    zmd(i,j,k)=0.
    za (i,j,k)=0.
  END DO
  END DO
  END DO
! entrainement
! DO k=1,klev-1
!    DO i=1,klon
!       zed(i,k)=max(0.,mp(i,k)-mp(i,k+1))
!    END DO
! END DO

! calcul de la matrice d echange
! matrice de distribution de la masse entrainee en k

  DO k=1,klev-1
     DO i=1,klon
        zmd(i,k,k)=max(0.,mp(i,k)-mp(i,k+1))
     END DO
  END DO
  DO k=2,klev
     DO j=k-1,1,-1
        DO i=1,klon
           if(mp(i,j+1).ne.0) then
              zmd(i,j,k)=zmd(i,j+1,k)*min(1.,mp(i,j)/mp(i,j+1))
           ENDif
        END DO
     END DO
  END DO
  DO k=1,klev
     DO j=1,klev-1
        DO i=1,klon
           za(i,j,k)=max(0.,zmd(i,j+1,k)-zmd(i,j,k))
        END DO
     END DO
  END DO
!
! rajout du terme lie a l ascendance induite
!
  DO j=2,klev
     DO i=1,klon
        za(i,j,j-1)=za(i,j,j-1)+mp(i,j)
     END DO
  END DO
!
! tendances
!            
  DO k=1,klev
     DO j=1,klev
        DO i=1,klon
           zmfd(i,j)=zmfd(i,j)+za(i,j,k)*(x(i,k)-x(i,j))
        END DO
     END DO
  END DO
!
! =========================================
! calcul des tendances liees aux flux satures
! =========================================
  DO j=1,klev
     DO i=1,klon
        zmfa(i,j)=da(i,j)*(x(i,1)-x(i,j))
     END DO
  END DO
  DO k=1,klev
     DO j=1,klev
        DO i=1,klon
           zmfp(i,j)=zmfp(i,j)+phi(i,j,k)*(x(i,k)-x(i,j))
        END DO
     END DO
  END DO
  DO j=1,klev-1
     DO i=1,klon
        zmfu(i,j)=max(0.,upd(i,j+1)+dnd(i,j+1))*(x(i,j+1)-x(i,j))
     END DO
  END DO
  DO j=2,klev
     DO i=1,klon
        zmfu(i,j)=zmfu(i,j)+min(0.,upd(i,j)+dnd(i,j))*(x(i,j)-x(i,j-1))
     END DO
  END DO

! =========================================
! calcul final des tendances
! =========================================
  DO k=1, klev
     DO i=1, klon
        dx(i,k)=paprs(i,k)-paprs(i,k+1)
     ENDDO
  ENDDO
  pdtimeRG=pdtime*RG
!cdir collapse
  DO k=1, klev
     DO i=1, klon
        dx(i,k)=(zmfd(i,k)+zmfu(i,k)       &
                +zmfa(i,k)+zmfp(i,k))*pdtimeRG/dx(i,k)
        !          print*,'dx',k,dx(i,k)
     ENDDO
  ENDDO

! test de conservation du traceur
!      conserv=0.
!      DO k=1, klev
!        DO i=1, klon
!         conserv=conserv+dx(i,k)*   &
!        (paprs(i,k)-paprs(i,k+1))/RG
!        ENDDO
!      ENDDO
!      print *,'conserv',conserv
     
END SUBROUTINE cvltr
