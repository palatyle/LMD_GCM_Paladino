!
! $Id: caladvtrac.F 1446 2010-10-22 09:27:25Z emillour $
!
c
c
            SUBROUTINE caladvtrac(q,pbaru,pbarv ,
     *                   p ,masse, dq ,  teta,
     *                   flxw, pk)
c
      USE infotrac, ONLY : nqtot
      USE control_mod, ONLY : iapp_tracvl,planet_type,
     &                        force_conserv_tracer
      USE comconst_mod, ONLY: dtvr
      USE planetary_operations, ONLY: planetary_tracer_amount_from_mass
      IMPLICIT NONE
c
c     Auteurs:   F.Hourdin , P.Le Van, F.Forget, F.Codron  
c
c     F.Codron (10/99) : ajout humidite specifique pour eau vapeur
c=======================================================================
c
c       Shema de  Van Leer
c
c=======================================================================


#include "dimensions.h"
#include "paramet.h"

c   Arguments:
c   ----------
      REAL pbaru( ip1jmp1,llm ),pbarv( ip1jm,llm),masse(ip1jmp1,llm)
      REAL p( ip1jmp1,llmp1),q( ip1jmp1,llm,nqtot)
      real :: dq(ip1jmp1,llm,nqtot)
      REAL teta( ip1jmp1,llm),pk( ip1jmp1,llm)
      REAL               :: flxw(ip1jmp1,llm)

c  ..................................................................
c
c  .. dq n'est utilise et dimensionne que pour l'eau  vapeur et liqu.
c
c  ..................................................................
c
c   Local:
c   ------

      EXTERNAL  advtrac,minmaxq, qminimum
      INTEGER ij,l, iq, iapptrac
      REAL finmasse(ip1jmp1,llm), dtvrtrac
      REAL :: totaltracer_old(nqtot),totaltracer_new(nqtot)
      REAL :: ratio

! Ehouarn : try to fix tracer conservation issues:
      if (force_conserv_tracer) then
        do iq=1,nqtot
          call planetary_tracer_amount_from_mass(masse,q(:,:,iq),
     &                                totaltracer_old(iq))
        enddo
      endif
cc
c
! Earth-specific stuff for the first 2 tracers (water)
      if (planet_type.eq."earth") then
C initialisation
! CRisi: il faut gérer tous les traceurs si on veut pouvoir faire des
! isotopes
!        dq(:,:,1:2)=q(:,:,1:2)
        dq(:,:,1:nqtot)=q(:,:,1:nqtot)
       
c  test des valeurs minmax
cc        CALL minmaxq(q(1,1,1),1.e33,-1.e33,'Eau vapeur (a) ')
cc        CALL minmaxq(q(1,1,2),1.e33,-1.e33,'Eau liquide(a) ')
      endif ! of if (planet_type.eq."earth")
c   advection

        CALL advtrac( pbaru,pbarv, 
     *       p,  masse,q,iapptrac, teta,
     .       flxw, pk)

c

      IF( iapptrac.EQ.iapp_tracvl ) THEN
        if (planet_type.eq."earth") then
! Earth-specific treatment for the first 2 tracers (water)
c
cc          CALL minmaxq(q(1,1,1),1.e33,-1.e33,'Eau vapeur     ')
cc          CALL minmaxq(q(1,1,2),1.e33,-1.e33,'Eau liquide    ')

cc     ....  Calcul  de deltap  qu'on stocke dans finmasse   ...
c
          DO l = 1, llm
           DO ij = 1, ip1jmp1
             finmasse(ij,l) =  p(ij,l) - p(ij,l+1) 
           ENDDO
          ENDDO
          
          !write(*,*) 'caladvtrac 87'
          CALL qminimum( q, nqtot, finmasse )
          !write(*,*) 'caladvtrac 89'

          CALL SCOPY   ( ip1jmp1*llm, masse, 1, finmasse,       1 )
          CALL filtreg ( finmasse ,  jjp1,  llm, -2, 2, .TRUE., 1 )
c
c   *****  Calcul de dq pour l'eau , pour le passer a la physique ******
c   ********************************************************************
c
          dtvrtrac = iapp_tracvl * dtvr
c
           DO iq = 1 , nqtot
            DO l = 1 , llm
             DO ij = 1,ip1jmp1
             dq(ij,l,iq) = ( q(ij,l,iq) - dq(ij,l,iq) ) * finmasse(ij,l)
     *                               /  dtvrtrac
             ENDDO
            ENDDO
           ENDDO
c
        endif ! of if (planet_type.eq."earth")
        
        ! Ehouarn : try to fix tracer conservation after tracer advection
        if (force_conserv_tracer) then
          do iq=1,nqtot
            call planetary_tracer_amount_from_mass(masse,q(:,:,iq),
     &                                  totaltracer_new(iq))
            ratio=totaltracer_old(iq)/totaltracer_new(iq)
            q(:,:,iq)=q(:,:,iq)*ratio
          enddo
        endif !of if (force_conserv_tracer)
        
      ELSE ! i.e. iapptrac.NE.iapp_tracvl
        if (planet_type.eq."earth") then
! Earth-specific treatment for the first 2 tracers (water)
          dq(:,:,1:nqtot)=0.
        endif ! of if (planet_type.eq."earth")
      ENDIF ! of IF( iapptrac.EQ.iapp_tracvl )

      END


