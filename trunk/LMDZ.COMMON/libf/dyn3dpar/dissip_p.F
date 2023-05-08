!
! $Id: $
!
      SUBROUTINE dissip_p( vcov,ucov,teta,p, dv,du,dh )
c
      USE parallel_lmdz
      USE write_field_p
      USE comconst_mod, ONLY: dtdiss
      IMPLICIT NONE


c ..  Avec nouveaux operateurs star :  gradiv2 , divgrad2, nxgraro2  ...
c                                 (  10/01/98  )

c=======================================================================
c
c   Auteur:  P. Le Van
c   -------
c
c   Objet:
c   ------
c
c   Dissipation horizontale
c
c=======================================================================
c-----------------------------------------------------------------------
c   Declarations:
c   -------------

#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"
#include "comdissnew.h"
#include "comdissipn.h"

c   Arguments:
c   ----------

      REAL,INTENT(IN) :: vcov(ip1jm,llm) ! covariant meridional wind
      REAL,INTENT(IN) :: ucov(ip1jmp1,llm) ! covariant zonal wind
      REAL,INTENT(IN) :: teta(ip1jmp1,llm) ! potentail temperature
      REAL,INTENT(IN) :: p(ip1jmp1,llmp1) ! pressure
      ! tendencies (.../s) on covariant winds and potential temperature
      REAL,INTENT(OUT) :: dv(ip1jm,llm)
      REAL,INTENT(OUT) :: du(ip1jmp1,llm)
      REAL,INTENT(OUT) :: dh(ip1jmp1,llm)

c   Local:
c   ------

      REAL gdx(ip1jmp1,llm),gdy(ip1jm,llm)
      REAL grx(ip1jmp1,llm),gry(ip1jm,llm)
      REAL te1dt(llm),te2dt(llm),te3dt(llm)
      REAL deltapres(ip1jmp1,llm)

      INTEGER l,ij

      REAL  SSUM
      integer :: ijb,ije
c-----------------------------------------------------------------------
c   initialisations:
c   ----------------

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 
      DO l=1,llm
         te1dt(l) = tetaudiv(l) * dtdiss
         te2dt(l) = tetaurot(l) * dtdiss
         te3dt(l) = tetah(l)    * dtdiss
      ENDDO
c$OMP END DO NOWAIT
c      CALL initial0( ijp1llm, du )
c      CALL initial0( ijmllm , dv )
c      CALL initial0( ijp1llm, dh )
     
      ijb=ij_begin
      ije=ij_end

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 
      DO l=1,llm
        du(ijb:ije,l)=0
        dh(ijb:ije,l)=0
      ENDDO
c$OMP END DO NOWAIT
      
      if (pole_sud) ije=ij_end-iip1

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 
      DO l=1,llm
        dv(ijb:ije,l)=0
      ENDDO
c$OMP END DO NOWAIT
     
c-----------------------------------------------------------------------
c   Calcul de la dissipation:
c   -------------------------

c   Calcul de la partie   grad  ( div ) :
c   -------------------------------------
      
     
      
      IF(lstardis) THEN
c      IF (.FALSE.) THEN
         CALL gradiv2_p( llm,ucov,vcov,nitergdiv,gdx,gdy )
      ELSE
         CALL gradiv_p ( llm,ucov,vcov,nitergdiv,gdx,gdy )
      ENDIF


      ijb=ij_begin
      ije=ij_end
      if (pole_sud) ije=ij_end-iip1

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 
      DO l=1,llm
         if (pole_nord) then
           DO ij = 1, iip1
              gdx(     ij ,l) = 0.
           ENDDO
         endif
         
         if (pole_sud) then
           DO ij = 1, iip1
              gdx(ij+ip1jm,l) = 0.
           ENDDO
         endif
         
         if (pole_nord) ijb=ij_begin+iip1
         DO ij = ijb,ije
            du(ij,l) = du(ij,l) - te1dt(l) *gdx(ij,l)
         ENDDO

         if (pole_nord) ijb=ij_begin
         DO ij = ijb,ije
            dv(ij,l) = dv(ij,l) - te1dt(l) *gdy(ij,l)
         ENDDO

       ENDDO
c$OMP END DO NOWAIT
c   calcul de la partie   n X grad ( rot ):
c   ---------------------------------------

      IF(lstardis) THEN
c      IF (.FALSE.) THEN
         CALL nxgraro2_p( llm,ucov, vcov, nitergrot,grx,gry )
      ELSE
         CALL nxgrarot_p( llm,ucov, vcov, nitergrot,grx,gry )
      ENDIF



      ijb=ij_begin
      ije=ij_end
      if (pole_sud) ije=ij_end-iip1

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)       
      DO l=1,llm
         
         if (pole_nord) then
           DO ij = 1, iip1
              grx(ij,l) = 0.
           ENDDO
         endif
         
         if (pole_nord) ijb=ij_begin+iip1
         DO ij = ijb,ije
            du(ij,l) = du(ij,l) - te2dt(l) * grx(ij,l)
         ENDDO
         
         if (pole_nord) ijb=ij_begin
         DO ij =  ijb, ije
            dv(ij,l) = dv(ij,l) - te2dt(l) * gry(ij,l)
         ENDDO
      
      ENDDO
c$OMP END DO NOWAIT

c   calcul de la partie   div ( grad ):
c   -----------------------------------

        
      IF(lstardis) THEN
c      IF (.FALSE.) THEN
    
      ijb=ij_begin
      ije=ij_end

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
       DO l = 1, llm
          DO ij = ijb, ije
            deltapres(ij,l) = AMAX1( 0.,  p(ij,l) - p(ij,l+1) )
          ENDDO
       ENDDO
c$OMP END DO NOWAIT
         CALL divgrad2_p( llm,teta, deltapres  ,niterh, gdx )
      ELSE
         CALL divgrad_p ( llm,teta, niterh, gdx        )
      ENDIF

c      call write_field3d_p('gdx2',reshape(gdx,(/iip1,jmp1,llm/)))
c      stop

      ijb=ij_begin
      ije=ij_end
      
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
      DO l = 1,llm
         DO ij = ijb,ije
            dh( ij,l ) = dh( ij,l ) - te3dt(l) * gdx( ij,l )
         ENDDO
      ENDDO
c$OMP END DO NOWAIT

      RETURN
      END
