module exner_hyb_p_m

  IMPLICIT NONE

contains

  SUBROUTINE  exner_hyb_p ( ngrid, ps, p, pks, pk, pkf )

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
    USE parallel_lmdz
    USE comvert_mod, ONLY: preff
    USE comconst_mod, ONLY: jmp1,kappa,cpp,r
    !
    include "dimensions.h"
    include "paramet.h"
    include "comgeom.h"

    INTEGER  ngrid
    REAL p(ngrid,llmp1),pk(ngrid,llm)
    REAL, optional:: pkf(ngrid,llm)
    REAL ps(ngrid),pks(ngrid)
    REAL alpha(ngrid,llm),beta(ngrid,llm)

    !    .... variables locales   ...

    INTEGER l, ij
    REAL unpl2k,dellta

    INTEGER ije,ijb,jje,jjb
    logical,save :: firstcall=.true.
    !$OMP THREADPRIVATE(firstcall) 
    character(len=*),parameter :: modname="exner_hyb_p"

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

    !$OMP BARRIER

    ! Specific behaviour for Shallow Water (1 vertical layer) case:
    if (llm.eq.1) then

       ! Compute pks(:),pk(:),pkf(:)
       ijb=ij_begin
       ije=ij_end
       !$OMP DO SCHEDULE(STATIC)
       DO ij=ijb, ije
          pks(ij) = (cpp/preff) * ps(ij)
          pk(ij,1) = .5*pks(ij)
          if (present(pkf)) pkf(ij,1)=pk(ij,1)
       ENDDO
       !$OMP ENDDO

       !$OMP BARRIER
       if (present(pkf)) then
          jjb=jj_begin
          jje=jj_end
          CALL filtreg_p ( pkf,jjb,jje, jmp1, llm, 2, 1, .TRUE., 1 )
       end if

       ! our work is done, exit routine
       return
    endif ! of if (llm.eq.1)

    ! General case:

    unpl2k    = 1.+ 2.* kappa

    !     -------------
    !     Calcul de pks
    !     -------------

    ijb=ij_begin
    ije=ij_end

    !$OMP DO SCHEDULE(STATIC)
    DO   ij  = ijb, ije
       pks(ij) = cpp * ( ps(ij)/preff ) ** kappa
    ENDDO
    !$OMP ENDDO
    ! Synchro OPENMP ici

    !$OMP BARRIER
    !
    !
    !    .... Calcul des coeff. alpha et beta  pour la couche l = llm ..
    !
    !$OMP DO SCHEDULE(STATIC)
    DO     ij      = ijb,ije
       alpha(ij,llm) = 0.
       beta (ij,llm) = 1./ unpl2k
    ENDDO
    !$OMP ENDDO NOWAIT
    !
    !     ... Calcul des coeff. alpha et beta  pour l = llm-1  a l = 2 ...
    !
    DO l = llm -1 , 2 , -1
       !
       !$OMP DO SCHEDULE(STATIC)
       DO ij = ijb, ije
          dellta = p(ij,l)* unpl2k + p(ij,l+1)* ( beta(ij,l+1)-unpl2k )
          alpha(ij,l)  = - p(ij,l+1) / dellta * alpha(ij,l+1)
          beta (ij,l)  =   p(ij,l  ) / dellta   
       ENDDO
       !$OMP ENDDO NOWAIT
    ENDDO

    !  ***********************************************************************
    !     .....  Calcul de pk pour la couche 1 , pres du sol  ....
    !
    !$OMP DO SCHEDULE(STATIC)
    DO   ij   = ijb, ije
       pk(ij,1) = ( p(ij,1)*pks(ij) - 0.5*alpha(ij,2)*p(ij,2) )  / &
            (  p(ij,1)* (1.+kappa) + 0.5*( beta(ij,2)-unpl2k )* p(ij,2) )
    ENDDO
    !$OMP ENDDO NOWAIT
    !
    !    ..... Calcul de pk(ij,l) , pour l = 2 a l = llm  ........
    !
    DO l = 2, llm
       !$OMP DO SCHEDULE(STATIC)
       DO   ij   = ijb, ije
          pk(ij,l) = alpha(ij,l) + beta(ij,l) * pk(ij,l-1)
       ENDDO
       !$OMP ENDDO NOWAIT        
    ENDDO

    if (present(pkf)) then
       !    calcul de pkf

       DO l = 1, llm
          !$OMP DO SCHEDULE(STATIC)
          DO   ij   = ijb, ije
             pkf(ij,l)=pk(ij,l)
          ENDDO
          !$OMP ENDDO NOWAIT             
       ENDDO

       !$OMP BARRIER

       jjb=jj_begin
       jje=jj_end
       CALL filtreg_p ( pkf,jjb,jje, jmp1, llm, 2, 1, .TRUE., 1 )
    end if

  END SUBROUTINE exner_hyb_p

end module exner_hyb_p_m

