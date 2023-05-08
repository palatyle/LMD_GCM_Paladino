










! $Id: pres2lev.F 1179 2009-06-11 14:18:47Z jghattas $
!
MODULE pres2lev_mod

CONTAINS 

!******************************************************
SUBROUTINE pres2lev(varo,varn,lmo,lmn,po,pn,ni,nj,ok_invertp)
!
! interpolation lineaire pour passer
! a une nouvelle discretisation verticale pour
! les variables de GCM
! Francois Forget (01/1995)
! MOdif remy roca 12/97 pour passer de pres2sig
! Modif F.Codron 07/08 po en 3D
!**********************************************************

  IMPLICIT NONE

!   Declarations:
! ==============
!
!  ARGUMENTS
!  """""""""
  LOGICAL, INTENT(IN) :: ok_invertp
  INTEGER, INTENT(IN) :: lmo ! dimensions ancienne couches
  INTEGER, INTENT(IN) :: lmn ! dimensions nouvelle couches
  
  INTEGER, INTENT(IN) :: ni,nj ! nombre de point horizontal
  REAL, INTENT(IN) :: po(ni*nj,lmo) ! niveau de pression ancienne grille
  REAL, INTENT(IN) :: pn(ni*nj,lmn) ! niveau de pression nouvelle grille

  REAL, INTENT(IN)  :: varo(ni*nj,lmo) ! var dans l'ancienne grille
  REAL, INTENT(OUT) :: varn(ni*nj,lmn) ! var dans la nouvelle grille

  REAL :: zvaro(ni*nj,lmo),zpo(ni*nj,lmo)

! Autres variables
! """"""""""""""""
  INTEGER ::  ln ,lo, k
  REAL    :: coef


! Inversion de l'ordre des niveaux verticaux
  IF (ok_invertp) THEN
    DO lo=1,lmo
      DO k=1,ni*nj
        zpo(k,lo)=po(k,lmo+1-lo)
        zvaro(k,lo)=varo(k,lmo+1-lo)
      ENDDO
    ENDDO
  ELSE
    DO lo=1,lmo
      DO k=1,ni*nj
        zpo(k,lo)=po(k,lo)
        zvaro(k,lo)=varo(k,lo)
      ENDDO
    ENDDO
  ENDIF 

  DO ln=1,lmn
    DO lo=1,lmo-1
      DO k=1,ni*nj
        IF (pn(k,ln) >= zpo(k,1) ) THEN
          varn(k,ln) = zvaro(k,1)
        ELSE IF (pn(k,ln) <= zpo(k,lmo)) THEN
          varn(k,ln) = zvaro(k,lmo)
        ELSE IF ( pn(k,ln) <= zpo(k,lo) .AND. pn(k,ln) > zpo(k,lo+1) ) THEN
          coef = (pn(k,ln)-zpo(k,lo)) / (zpo(k,lo+1)-zpo(k,lo))
          varn(k,ln) = zvaro(k,lo) + coef*(zvaro(k,lo+1)-zvaro(k,lo))
        ENDIF
         
      ENDDO  
    ENDDO
  ENDDO                

END SUBROUTINE pres2lev    

END MODULE pres2lev_mod
