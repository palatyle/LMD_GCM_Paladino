

      SUBROUTINE filtreg_p ( champ, ibeg, iend, nlat, nbniv, 
     &     ifiltre, iaire, griscal ,iter)
      USE parallel_lmdz, only : OMP_CHUNK
      USE mod_filtre_fft
      USE timer_filtre
      
      USE filtreg_mod
      
      IMPLICIT NONE
      
c=======================================================================
c
c   Auteur: P. Le Van        07/10/97
c   ------
c
c   Objet: filtre matriciel longitudinal ,avec les matrices precalculees
c                     pour l'operateur  Filtre    .
c   ------
c
c   Arguments:
c   ----------
c
c      
c      ibeg..iend            lattitude a filtrer
c      nlat                  nombre de latitudes du champ
c      nbniv                 nombre de niveaux verticaux a filtrer
c      champ(iip1,nblat,nbniv)  en entree : champ a filtrer
c                            en sortie : champ filtre
c      ifiltre               +1  Transformee directe
c                            -1  Transformee inverse
c                            +2  Filtre directe
c                            -2  Filtre inverse
c
c      iaire                 1   si champ intensif
c                            2   si champ extensif (pondere par les aires)
c
c      iter                  1   filtre simple
c
c=======================================================================
c
c
c                      Variable Intensive
c                ifiltre = 1     filtre directe
c                ifiltre =-1     filtre inverse
c
c                      Variable Extensive
c                ifiltre = 2     filtre directe
c                ifiltre =-2     filtre inverse
c
c
#include "dimensions.h"
#include "paramet.h"
#include "coefils.h"
c
      INTEGER ibeg,iend,nlat,nbniv,ifiltre,iter
      INTEGER i,j,l,k
      INTEGER iim2,immjm
      INTEGER jdfil1,jdfil2,jffil1,jffil2,jdfil,jffil
      
      REAL  champ( iip1,nlat,nbniv)
      
      LOGICAL    griscal
      INTEGER    hemisph, iaire
      
      REAL :: champ_fft(iip1,nlat,nbniv)
      REAL :: champ_in(iip1,nlat,nbniv)
      
      LOGICAL,SAVE     :: first=.TRUE.
c$OMP THREADPRIVATE(first) 

      REAL, DIMENSION(iip1,nlat,nbniv) :: champ_loc
      INTEGER :: ll_nb, nbniv_loc
      REAL, SAVE :: sdd12(iim,4)
c$OMP THREADPRIVATE(sdd12) 

      INTEGER, PARAMETER :: type_sddu=1
      INTEGER, PARAMETER :: type_sddv=2
      INTEGER, PARAMETER :: type_unsddu=3
      INTEGER, PARAMETER :: type_unsddv=4

      INTEGER :: sdd1_type, sdd2_type

      IF (first) THEN
         sdd12(1:iim,type_sddu) = sddu(1:iim)
         sdd12(1:iim,type_sddv) = sddv(1:iim)
         sdd12(1:iim,type_unsddu) = unsddu(1:iim)
         sdd12(1:iim,type_unsddv) = unsddv(1:iim)

         CALL Init_timer
         first=.FALSE.
      ENDIF

c$OMP MASTER      
      CALL start_timer
c$OMP END MASTER

c-------------------------------------------------------c

      IF(ifiltre.EQ.1.or.ifiltre.EQ.-1) 
     &     STOP'Pas de transformee simple dans cette version'
      
      IF( iter.EQ. 2 )  THEN
         PRINT *,' Pas d iteration du filtre dans cette version !'
     &        , ' Utiliser old_filtreg et repasser !'
         STOP
      ENDIF

      IF( ifiltre.EQ. -2 .AND..NOT.griscal )     THEN
         PRINT *,' Cette routine ne calcule le filtre inverse que '
     &        , ' sur la grille des scalaires !'
         STOP
      ENDIF

      IF( ifiltre.NE.2 .AND.ifiltre.NE. - 2 )  THEN
         PRINT *,' Probleme dans filtreg car ifiltre NE 2 et NE -2'
     &        , ' corriger et repasser !'
         STOP
      ENDIF
c

      iim2   = iim * iim
      immjm  = iim * jjm
c
c
      IF( griscal )   THEN
         IF( nlat. NE. jjp1 )  THEN
            PRINT  1111
            STOP
         ELSE
c     
            IF( iaire.EQ.1 )  THEN
               sdd1_type = type_sddv
               sdd2_type = type_unsddv
            ELSE
               sdd1_type = type_unsddv
               sdd2_type = type_sddv
            ENDIF
c
            jdfil1 = 2
            jffil1 = jfiltnu
            jdfil2 = jfiltsu
            jffil2 = jjm
         ENDIF
      ELSE
         IF( nlat.NE.jjm )  THEN
            PRINT  2222
            STOP
         ELSE
c
            IF( iaire.EQ.1 )  THEN
               sdd1_type = type_sddu
               sdd2_type = type_unsddu
            ELSE
               sdd1_type = type_unsddu
               sdd2_type = type_sddu
            ENDIF
c     
            jdfil1 = 1
            jffil1 = jfiltnv
            jdfil2 = jfiltsv
            jffil2 = jjm
         ENDIF
      ENDIF
c      
      DO hemisph = 1, 2
c     
         IF ( hemisph.EQ.1 )  THEN
cym
            jdfil = max(jdfil1,ibeg)
            jffil = min(jffil1,iend)
         ELSE
cym
            jdfil = max(jdfil2,ibeg)
            jffil = min(jffil2,iend)
         ENDIF


cccccccccccccccccccccccccccccccccccccccccccc
c Utilisation du filtre classique
cccccccccccccccccccccccccccccccccccccccccccc

         IF (.NOT. use_filtre_fft) THEN
      
c     !---------------------------------!
c     ! Agregation des niveau verticaux !
c     ! uniquement necessaire pour une  !
c     ! execution OpenMP                !
c     !---------------------------------!
            ll_nb = 0
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
            DO l = 1, nbniv
               ll_nb = ll_nb+1
               DO j = jdfil,jffil
                  DO i = 1, iim
                     champ_loc(i,j,ll_nb) = 
     &                    champ(i,j,l) * sdd12(i,sdd1_type)
                  ENDDO
               ENDDO
            ENDDO
c$OMP END DO NOWAIT

            nbniv_loc = ll_nb

            IF( hemisph.EQ.1 )      THEN
               
               IF( ifiltre.EQ.-2 )   THEN
                  DO j = jdfil,jffil
#ifdef BLAS
                     CALL DGEMM("N", "N", iim, nbniv_loc, iim, 1.0, 
     &                    matrinvn(1,1,j), iim, 
     &                    champ_loc(1,j,1), iip1*nlat, 0.0,
     &                    champ_fft(1,j-jdfil+1,1), iip1*nlat)
#else
                     champ_fft(:iim,j-jdfil+1,:)
     &                    =matmul(matrinvn(:,:,j),champ_loc(:iim,j,:))
#endif
                  ENDDO
                  
               ELSE IF ( griscal )     THEN
                  DO j = jdfil,jffil
#ifdef BLAS
                     CALL DGEMM("N", "N", iim, nbniv_loc, iim, 1.0, 
     &                    matriceun(1,1,j), iim, 
     &                    champ_loc(1,j,1), iip1*nlat, 0.0,
     &                    champ_fft(1,j-jdfil+1,1), iip1*nlat)
#else
                     champ_fft(:iim,j-jdfil+1,:)
     &                    =matmul(matriceun(:,:,j),champ_loc(:iim,j,:))
#endif
                  ENDDO
                  
               ELSE 
                  DO j = jdfil,jffil
#ifdef BLAS
                     CALL DGEMM("N", "N", iim, nbniv_loc, iim, 1.0, 
     &                    matricevn(1,1,j), iim, 
     &                    champ_loc(1,j,1), iip1*nlat, 0.0, 
     &                    champ_fft(1,j-jdfil+1,1), iip1*nlat)
#else
                     champ_fft(:iim,j-jdfil+1,:)
     &                    =matmul(matricevn(:,:,j),champ_loc(:iim,j,:))
#endif
                  ENDDO
                  
               ENDIF
               
            ELSE
               
               IF( ifiltre.EQ.-2 )   THEN
                  DO j = jdfil,jffil
#ifdef BLAS
                     CALL DGEMM("N", "N", iim, nbniv_loc, iim, 1.0, 
     &                    matrinvs(1,1,j-jfiltsu+1), iim, 
     &                    champ_loc(1,j,1), iip1*nlat, 0.0, 
     &                    champ_fft(1,j-jdfil+1,1), iip1*nlat)
#else
                     champ_fft(:iim,j-jdfil+1,:)
     &                    =matmul(matrinvs(:,:,j-jfiltsu+1),
     &                            champ_loc(:iim,j,:))
#endif
                  ENDDO
                  
               ELSE IF ( griscal )     THEN
                  
                  DO j = jdfil,jffil
#ifdef BLAS
                     CALL DGEMM("N", "N", iim, nbniv_loc, iim, 1.0, 
     &                    matriceus(1,1,j-jfiltsu+1), iim, 
     &                    champ_loc(1,j,1), iip1*nlat, 0.0, 
     &                    champ_fft(1,j-jdfil+1,1), iip1*nlat)
#else
                     champ_fft(:iim,j-jdfil+1,:)
     &                    =matmul(matriceus(:,:,j-jfiltsu+1),
     &                            champ_loc(:iim,j,:))
#endif
                  ENDDO
                  
               ELSE 
                  
                  DO j = jdfil,jffil
#ifdef BLAS
                     CALL DGEMM("N", "N", iim, nbniv_loc, iim, 1.0, 
     &                    matricevs(1,1,j-jfiltsv+1), iim, 
     &                    champ_loc(1,j,1), iip1*nlat, 0.0, 
     &                    champ_fft(1,j-jdfil+1,1), iip1*nlat)
#else
                     champ_fft(:iim,j-jdfil+1,:)
     &                    =matmul(matricevs(:,:,j-jfiltsv+1),
     &                            champ_loc(:iim,j,:))
#endif
                  ENDDO
                  
               ENDIF
               
            ENDIF
!     c     
            IF( ifiltre.EQ.2 )  THEN
               
c     !-------------------------------------!
c     ! Dés-agregation des niveau verticaux !
c     ! uniquement necessaire pour une      !
c     ! execution OpenMP                    !
c     !-------------------------------------!
               ll_nb = 0
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
               DO l = 1, nbniv
                  ll_nb = ll_nb + 1
                  DO j = jdfil,jffil
                     DO i = 1, iim
                        champ( i,j,l ) = (champ_loc(i,j,ll_nb) 
     &                       + champ_fft(i,j-jdfil+1,ll_nb))
     &                       * sdd12(i,sdd2_type)
                     ENDDO
                  ENDDO
               ENDDO
c$OMP END DO NOWAIT
               
            ELSE
               
               ll_nb = 0
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
               DO l = 1, nbniv_loc
                  ll_nb = ll_nb + 1
                  DO j = jdfil,jffil
                     DO i = 1, iim
                        champ( i,j,l ) = (champ_loc(i,j,ll_nb) 
     &                       - champ_fft(i,j-jdfil+1,ll_nb))
     &                       * sdd12(i,sdd2_type)
                     ENDDO
                  ENDDO
               ENDDO
c$OMP END DO NOWAIT
               
            ENDIF
            
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
            DO l = 1, nbniv
               DO j = jdfil,jffil
                  champ( iip1,j,l ) = champ( 1,j,l )
               ENDDO
            ENDDO
c$OMP END DO NOWAIT
            
ccccccccccccccccccccccccccccccccccccccccccccc
c Utilisation du filtre FFT
ccccccccccccccccccccccccccccccccccccccccccccc
        
         ELSE
       
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
            DO l=1,nbniv
               DO j=jdfil,jffil
                  DO  i = 1, iim
                     champ( i,j,l)= champ(i,j,l)*sdd12(i,sdd1_type)
                     champ_fft( i,j,l) = champ(i,j,l)
                  ENDDO
               ENDDO
            ENDDO
c$OMP END DO NOWAIT

            IF (jdfil<=jffil) THEN
               IF( ifiltre. EQ. -2 )   THEN
                  CALL Filtre_inv_fft(champ_fft,nlat,jdfil,jffil,nbniv) 
               ELSE IF ( griscal )     THEN
                  CALL Filtre_u_fft(champ_fft,nlat,jdfil,jffil,nbniv)
               ELSE
                  CALL Filtre_v_fft(champ_fft,nlat,jdfil,jffil,nbniv)
               ENDIF
            ENDIF


            IF( ifiltre.EQ. 2 )  THEN
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)          
               DO l=1,nbniv
                  DO j=jdfil,jffil
                     DO  i = 1, iim
                        champ( i,j,l)=(champ(i,j,l)+champ_fft(i,j,l))
     &                       *sdd12(i,sdd2_type)
                     ENDDO
                  ENDDO
               ENDDO
c$OMP END DO NOWAIT	  
            ELSE
        
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)  
               DO l=1,nbniv
                  DO j=jdfil,jffil
                     DO  i = 1, iim
                        champ(i,j,l)=(champ(i,j,l)-champ_fft(i,j,l))
     &                       *sdd12(i,sdd2_type)
                     ENDDO
                  ENDDO
               ENDDO
c$OMP END DO NOWAIT          
            ENDIF
c
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 
            DO l=1,nbniv
               DO j=jdfil,jffil
!            champ_FFT( iip1,j,l ) = champ_FFT( 1,j,l )
                  champ( iip1,j,l ) = champ( 1,j,l )
               ENDDO
            ENDDO
c$OMP END DO NOWAIT          	
         ENDIF 
c Fin de la zone de filtrage

	
      ENDDO

!      DO j=1,nlat
!     
!          PRINT *,"check FFT ----> Delta(",j,")=",
!     &            sum(champ(:,j,:)-champ_fft(:,j,:))/sum(champ(:,j,:)),
!     &            sum(champ(:,j,:)-champ_in(:,j,:))/sum(champ(:,j,:)) 
!      ENDDO
      
!          PRINT *,"check FFT ----> Delta(",j,")=",
!     &            sum(champ-champ_fft)/sum(champ)
!      
      
c
 1111 FORMAT(//20x,'ERREUR dans le dimensionnement du tableau  CHAMP a 
     &     filtrer, sur la grille des scalaires'/)
 2222 FORMAT(//20x,'ERREUR dans le dimensionnement du tableau CHAMP a fi
     &     ltrer, sur la grille de V ou de Z'/)
c$OMP MASTER      
      CALL stop_timer
c$OMP END MASTER
      RETURN
      END
