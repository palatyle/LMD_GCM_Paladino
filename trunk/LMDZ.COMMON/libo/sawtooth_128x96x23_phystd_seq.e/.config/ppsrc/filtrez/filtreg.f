










!
! $Header$
!
      SUBROUTINE filtreg ( champ, nlat, nbniv, ifiltre,iaire,
     &     griscal ,iter)
      
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
c      nblat                 nombre de latitudes a filtrer
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
!-----------------------------------------------------------------------
!   INCLUDE 'dimensions.h'
!
!   dimensions.h contient les dimensions du modele
!   ndm est tel que iim=2**ndm
!-----------------------------------------------------------------------

      INTEGER iim,jjm,llm,ndm

      PARAMETER (iim= 128,jjm=96,llm=23,ndm=1)

!-----------------------------------------------------------------------
!
! $Header$
!
!
!  ATTENTION!!!!: ce fichier include est compatible format fixe/format libre
!                 veillez  n'utiliser que des ! pour les commentaires
!                 et  bien positionner les & des lignes de continuation
!                 (les placer en colonne 6 et en colonne 73)
!
!
!-----------------------------------------------------------------------
!   INCLUDE 'paramet.h'

      INTEGER  iip1,iip2,iip3,jjp1,llmp1,llmp2,llmm1
      INTEGER  kftd,ip1jm,ip1jmp1,ip1jmi1,ijp1llm
      INTEGER  ijmllm,mvar
      INTEGER jcfil,jcfllm

      PARAMETER( iip1= iim+1,iip2=iim+2,iip3=iim+3                       &
     &    ,jjp1=jjm+1-1/jjm)
      PARAMETER( llmp1 = llm+1,  llmp2 = llm+2, llmm1 = llm-1 )
      PARAMETER( kftd  = iim/2 -ndm )
      PARAMETER( ip1jm  = iip1*jjm,  ip1jmp1= iip1*jjp1 )
      PARAMETER( ip1jmi1= ip1jm - iip1 )
      PARAMETER( ijp1llm= ip1jmp1 * llm, ijmllm= ip1jm * llm )
      PARAMETER( mvar= ip1jmp1*( 2*llm+1) + ijmllm )
      PARAMETER( jcfil=jjm/2+5, jcfllm=jcfil*llm )

!-----------------------------------------------------------------------
!
! $Header$
!
      COMMON/coefils/jfiltnu,jfiltsu,jfiltnv,jfiltsv,sddu(iim),sddv(iim)&
     & ,unsddu(iim),unsddv(iim),coefilu(iim,jjm),coefilv(iim,jjm),      &
     & modfrstu(jjm),modfrstv(jjm),eignfnu(iim,iim),eignfnv(iim,iim)    &
     & ,coefilu2(iim,jjm),coefilv2(iim,jjm)
!c
      INTEGER jfiltnu,jfiltsu,jfiltnv,jfiltsv,modfrstu,modfrstv
      REAL    sddu,sddv,unsddu,unsddv,coefilu,coefilv,eignfnu,eignfnv
      REAL    coefilu2,coefilv2

      INTEGER    nlat,nbniv,ifiltre,iter
      INTEGER    i,j,l,k
      INTEGER    iim2,immjm
      INTEGER    jdfil1,jdfil2,jffil1,jffil2,jdfil,jffil

      REAL       champ( iip1,nlat,nbniv)

      REAL       eignq(iim,nlat,nbniv), sdd1(iim),sdd2(iim)
      LOGICAL    griscal
      INTEGER    hemisph, iaire

      LOGICAL,SAVE     :: first=.TRUE.

      REAL, SAVE :: sdd12(iim,4)

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

         first=.FALSE.
      ENDIF

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
      
      iim2   = iim * iim
      immjm  = iim * jjm

      IF( griscal )   THEN
         IF( nlat. NE. jjp1 )  THEN
            PRINT  1111
            STOP
         ELSE
            
            IF( iaire.EQ.1 )  THEN
               sdd1_type = type_sddv
               sdd2_type = type_unsddv
            ELSE
               sdd1_type = type_unsddv
               sdd2_type = type_sddv
            ENDIF

c            IF( iaire.EQ.1 )  THEN
c               CALL SCOPY(  iim,    sddv, 1,  sdd1, 1 ) 
c               CALL SCOPY(  iim,  unsddv, 1,  sdd2, 1 )
c            ELSE
c               CALL SCOPY(  iim,  unsddv, 1,  sdd1, 1 )
c               CALL SCOPY(  iim,    sddv, 1,  sdd2, 1 )
c            END IF
            
            jdfil1 = 2
            jffil1 = jfiltnu
            jdfil2 = jfiltsu
            jffil2 = jjm
         END IF
      ELSE
         IF( nlat.NE.jjm )  THEN
            PRINT  2222
            STOP
         ELSE
            
            IF( iaire.EQ.1 )  THEN
               sdd1_type = type_sddu
               sdd2_type = type_unsddu
            ELSE
               sdd1_type = type_unsddu
               sdd2_type = type_sddu
            ENDIF

c            IF( iaire.EQ.1 )  THEN
c               CALL SCOPY(  iim,    sddu, 1,  sdd1, 1 ) 
c               CALL SCOPY(  iim,  unsddu, 1,  sdd2, 1 )
c            ELSE
c               CALL SCOPY(  iim,  unsddu, 1,  sdd1, 1 )
c               CALL SCOPY(  iim,    sddu, 1,  sdd2, 1 )
c            END IF
            
            jdfil1 = 1
            jffil1 = jfiltnv
            jdfil2 = jfiltsv
            jffil2 = jjm
         END IF
      END IF
      
      DO hemisph = 1, 2
         
         IF ( hemisph.EQ.1 )  THEN
            jdfil = jdfil1
            jffil = jffil1
         ELSE
            jdfil = jdfil2
            jffil = jffil2
         END IF
         
         DO l = 1, nbniv
            DO j = jdfil,jffil
               DO i = 1, iim
                  champ(i,j,l) = champ(i,j,l) * sdd12(i,sdd1_type) ! sdd1(i)
               END DO
            END DO
         END DO
         
         IF( hemisph. EQ. 1 )      THEN
            
            IF( ifiltre. EQ. -2 )   THEN
               
               DO j = jdfil,jffil
                  eignq(:,j-jdfil+1,:)
     $                 = matmul(matrinvn(:,:,j), champ(:iim,j,:))
               END DO
               
            ELSE IF ( griscal )     THEN
               
               DO j = jdfil,jffil
                  eignq(:,j-jdfil+1,:)
     $                 = matmul(matriceun(:,:,j), champ(:iim,j,:))
               END DO
               
            ELSE 
               
               DO j = jdfil,jffil
                  eignq(:,j-jdfil+1,:)
     $                 = matmul(matricevn(:,:,j), champ(:iim,j,:))
               END DO
               
            ENDIF
            
         ELSE
            
            IF( ifiltre. EQ. -2 )   THEN
               
               DO j = jdfil,jffil
                  eignq(:,j-jdfil+1,:)
     $                 = matmul(matrinvs(:,:,j-jfiltsu+1),
     $                 champ(:iim,j,:))
               END DO
               
               
            ELSE IF ( griscal )     THEN
               
               DO j = jdfil,jffil
                  eignq(:,j-jdfil+1,:)
     $                 = matmul(matriceus(:,:,j-jfiltsu+1),
     $                 champ(:iim,j,:))
               END DO
                              
            ELSE 
               
               DO j = jdfil,jffil
                  eignq(:,j-jdfil+1,:)
     $                 = matmul(matricevs(:,:,j-jfiltsv+1),
     $                 champ(:iim,j,:))
               END DO
                              
            ENDIF
            
         ENDIF
         
         IF( ifiltre.EQ. 2 )  THEN
            
            DO l = 1, nbniv
               DO j = jdfil,jffil
                  DO i = 1, iim
                     champ( i,j,l ) = 
     &                    (champ(i,j,l) + eignq(i,j-jdfil+1,l))
     &                    * sdd12(i,sdd2_type) ! sdd2(i)
                  END DO
               END DO
            END DO

         ELSE

            DO l = 1, nbniv
               DO j = jdfil,jffil
                  DO i = 1, iim
                     champ( i,j,l ) = 
     &                    (champ(i,j,l) - eignq(i,j-jdfil+1,l))
     &                    * sdd12(i,sdd2_type) ! sdd2(i)
                  END DO
               END DO
            END DO

         ENDIF

         DO l = 1, nbniv
            DO j = jdfil,jffil
               champ( iip1,j,l ) = champ( 1,j,l )
            END DO
         END DO

     
      ENDDO

1111  FORMAT(//20x,'ERREUR dans le dimensionnement du tableau  CHAMP a 
     &     filtrer, sur la grille des scalaires'/)
2222  FORMAT(//20x,'ERREUR dans le dimensionnement du tableau CHAMP a fi
     &     ltrer, sur la grille de V ou de Z'/)
      RETURN
      END
