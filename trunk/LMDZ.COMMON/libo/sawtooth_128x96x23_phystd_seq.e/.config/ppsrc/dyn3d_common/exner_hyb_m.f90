










module exner_hyb_m

  USE comvert_mod, ONLY: preff
  USE comconst_mod, ONLY: jmp1,kappa,cpp,r,pi

  IMPLICIT NONE

contains

  SUBROUTINE  exner_hyb ( ngrid, ps, p, pks, pk, pkf )

    !     Auteurs :  P.Le Van  , Fr. Hourdin  .
    !    ..........
    !
    !    ....  ngrid, ps,p             sont des argum.d'entree  au sous-prog ...
    !    ....  pks,pk,pkf   sont des argum.de sortie au sous-prog ...
    !
    !   ************************************************************************
    !    Calcule la fonction d'Exner pk = Cp * (p/preff) ** kappa , aux milieux des 
    !    couches .   Pk(l) sera calcule aux milieux  des couches l ,entre les
    !    pressions p(l) et p(l+1) ,definis aux interfaces des llm couches .
    !   ************************************************************************
    !  .. N.B : Au sommet de l'atmosphere,  p(llm+1) = 0. , et ps et pks sont
    !    la pression et la fonction d'Exner  au  sol  .
    !
    !                                 -------- z
    !    A partir des relations  ( 1 ) p*dz(pk) = kappa *pk*dz(p)      et
    !                            ( 2 ) pk(l) = alpha(l)+ beta(l)*pk(l-1)
    !    ( voir note de Fr.Hourdin )  ,
    !
    !    on determine successivement , du haut vers le bas des couches, les 
    !    coef. alpha(llm),beta(llm) .,.,alpha(l),beta(l),,,alpha(2),beta(2), 
    !    puis pk(ij,1). Ensuite ,on calcule,du bas vers le haut des couches,  
    !     pk(ij,l)  donne  par la relation (2),  pour l = 2 a l = llm .
    !
    !
    !
    include "dimensions.h"
    include "paramet.h"
    include "comgeom.h"

    INTEGER  ngrid
    REAL p(ngrid,llmp1),pk(ngrid,llm)
    real, optional:: pkf(ngrid,llm)
    REAL ps(ngrid),pks(ngrid), alpha(ngrid,llm),beta(ngrid,llm)

    !    .... variables locales   ...

    INTEGER l, ij
    REAL unpl2k,dellta

    logical,save :: firstcall=.true.
    character(len=*),parameter :: modname="exner_hyb"

    ! Sanity check
    if (firstcall) then
       ! sanity checks for Shallow Water case (1 vertical layer)
       if (llm.eq.1) then
          if (kappa.ne.1) then
             call abort_gcm(modname, &
                  "kappa!=1 , but running in Shallow Water mode!!",42)
          endif
          if (cpp.ne.r) then
             call abort_gcm(modname, &
                  "cpp!=r , but running in Shallow Water mode!!",42)
          endif
       endif ! of if (llm.eq.1)

       firstcall=.false.
    endif ! of if (firstcall)

    ! Specific behaviour for Shallow Water (1 vertical layer) case:
    if (llm.eq.1) then

       ! Compute pks(:),pk(:),pkf(:)

       DO   ij  = 1, ngrid
          pks(ij) = (cpp/preff) * ps(ij)
          pk(ij,1) = .5*pks(ij)
       ENDDO

       if (present(pkf)) then
          pkf = pk
          CALL filtreg ( pkf, jmp1, llm, 2, 1, .TRUE., 1 ) 
       end if

       ! our work is done, exit routine
       return
    endif ! of if (llm.eq.1)

    ! General case:

    unpl2k    = 1.+ 2.* kappa

    !     -------------
    !     Calcul de pks
    !     -------------

    DO   ij  = 1, ngrid
       pks(ij) = cpp * ( ps(ij)/preff ) ** kappa
    ENDDO

    !    .... Calcul des coeff. alpha et beta  pour la couche l = llm ..
    !
    DO     ij      = 1, ngrid
       alpha(ij,llm) = 0.
       beta (ij,llm) = 1./ unpl2k
    ENDDO
    !
    !     ... Calcul des coeff. alpha et beta  pour l = llm-1  a l = 2 ...
    !
    DO l = llm -1 , 2 , -1
       !
       DO ij = 1, ngrid
          dellta = p(ij,l)* unpl2k + p(ij,l+1)* ( beta(ij,l+1)-unpl2k )
          alpha(ij,l)  = - p(ij,l+1) / dellta * alpha(ij,l+1)
          beta (ij,l)  =   p(ij,l  ) / dellta   
       ENDDO
    ENDDO

    !  ***********************************************************************
    !     .....  Calcul de pk pour la couche 1 , pres du sol  ....
    !
    DO   ij   = 1, ngrid
       pk(ij,1) = ( p(ij,1)*pks(ij) - 0.5*alpha(ij,2)*p(ij,2) )  / &
            (  p(ij,1)* (1.+kappa) + 0.5*( beta(ij,2)-unpl2k )* p(ij,2) )
    ENDDO
    !
    !    ..... Calcul de pk(ij,l) , pour l = 2 a l = llm  ........
    !
    DO l = 2, llm
       DO   ij   = 1, ngrid
          pk(ij,l) = alpha(ij,l) + beta(ij,l) * pk(ij,l-1)
       ENDDO
    ENDDO

    if (present(pkf)) then
       !    calcul de pkf
       pkf = pk
       CALL filtreg ( pkf, jmp1, llm, 2, 1, .TRUE., 1 )
    end if

  END SUBROUTINE exner_hyb

end module exner_hyb_m

