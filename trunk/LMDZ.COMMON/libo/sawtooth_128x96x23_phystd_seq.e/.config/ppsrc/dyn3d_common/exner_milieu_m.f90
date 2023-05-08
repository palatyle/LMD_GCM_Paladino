










module exner_milieu_m

  USE comvert_mod, ONLY: preff
  USE comconst_mod, ONLY: jmp1,kappa,cpp,r

  IMPLICIT NONE

contains

  SUBROUTINE  exner_milieu ( ngrid, ps, p, pks, pk, pkf )
    !
    !     Auteurs :  F. Forget , Y. Wanherdrick
    ! P.Le Van  , Fr. Hourdin  .
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
    !     WARNING : CECI est une version speciale de exner_hyb originale
    !               Utilise dans la version martienne pour pouvoir 
    !               tourner avec des coordonnees verticales complexe
    !              => Il ne verifie PAS la condition la proportionalite en 
    !              energie totale/ interne / potentielle (F.Forget 2001)
    !    ( voir note de Fr.Hourdin )  ,
    !
    !
    include "dimensions.h"
    include "paramet.h"
    include "comgeom.h"

    INTEGER  ngrid
    REAL p(ngrid,llmp1),pk(ngrid,llm)
    real, optional:: pkf(ngrid,llm)
    REAL ps(ngrid),pks(ngrid)

    !    .... variables locales   ...

    INTEGER l, ij
    REAL dum1

    logical,save :: firstcall=.true.
    character(len=*),parameter :: modname="exner_milieu"

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

    !     -------------
    !     Calcul de pks
    !     -------------

    DO   ij  = 1, ngrid
       pks(ij) = cpp * ( ps(ij)/preff ) ** kappa
    ENDDO

    !    .... Calcul de pk  pour la couche l 
    !    --------------------------------------------
    !
    dum1 = cpp * (2*preff)**(-kappa) 
    DO l = 1, llm-1
       DO   ij   = 1, ngrid
          pk(ij,l) = dum1 * (p(ij,l) + p(ij,l+1))**kappa
       ENDDO
    ENDDO

    !    .... Calcul de pk  pour la couche l = llm ..
    !    (on met la meme distance (en log pression)  entre Pk(llm)
    !    et Pk(llm -1) qu'entre Pk(llm-1) et Pk(llm-2)

    DO   ij   = 1, ngrid
       pk(ij,llm) = pk(ij,llm-1)**2 / pk(ij,llm-2)
    ENDDO

    if (present(pkf)) then
       !    calcul de pkf
       pkf = pk
       CALL filtreg ( pkf, jmp1, llm, 2, 1, .TRUE., 1 )
    end if

  END SUBROUTINE exner_milieu

end module exner_milieu_m

