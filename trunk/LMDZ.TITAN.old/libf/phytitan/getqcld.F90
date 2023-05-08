!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!        SUBROUTINE getoptcld(WLN,RADIUS,Q_EXT,Q_SCT,Q_ABS,Q_BAR) 
!
!        Obtention des QEXT,Q_SCT,Q_ABS,Q_BAR pour une particule de rayon RADIUS a la longueur 
!        d'onde WLN .
!
!        ARGUMENTS D'ENTREE :
!              WLN : Longueur d'onde traitee (en metres !)
!           RADIUS : Rayon de la particule (en metres !)
!
!        ARGUMENT DE SORTIE :
!           Q_EXT : section efficace d'extinction
!           Q_SCT : section efficace de diffusion
!           Q_ABS : section efficace d'absorption
!           Q_BAR : Parametre d'asymetrie
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       SUBROUTINE getoptcld(WLN,RADIUS,Q_EXT,Q_SCT,Q_ABS,Q_BAR)
         USE optcld
         IMPLICIT NONE
! ------ INPUT
         REAL   ,INTENT(in)    ::  WLN,RADIUS
! ------ OUTPUT
         REAL   ,INTENT(out)   ::  Q_EXT,Q_SCT,Q_ABS,Q_BAR
! ------ LOCAL/COMMON 
         REAL                  :: tmp
         REAL,EXTERNAL         :: get_qopt
        
!        INTERPOLATION/EXTRAPOLATION de QEXT,QABS,QBAR et CALCUL de QSCT
!        notes : 
!        Comme les indices optiques des gaz sont peu variables en ldo et qu'on
!        approxime la goutte comme etant composee d'hydrocarbone seulement, on
!        a les relations suivantes :
!           sigma(r,ldo) = sigma(r0,ldo*r0/r) * (r/r0)**2.
!           gg(r,ldo)    = gg(r0,ldo*r0/r)
!
!        La routine get_qopt calcule sigma(r0,ldo*r0/r) ou gg(r0,ldo*r0/r) (selon les inputs)      
!           ====> il ne reste plus qu'a multiplier par (r/r0)**2. les sections efficaces :)
!    
!        ------------
!        QEXT   (attention : ltq_ex car on travaille en log dans get_qopt)
!        ------------
          tmp=get_qopt(radius,wln,A_ex,B_ex,fmin_ex,ltq_ex)
          Q_EXT=tmp*(radius/r0cld)**2. 
!        ------------
!        QABS   (attention : ltq_ab car on travaille en log dans get_qopt)
!        ------------
          tmp=get_qopt(radius,wln,A_ab,B_ab,fmin_ab,ltq_ab)
          Q_ABS=tmp*(radius/r0cld)**2.
!        ------------
!        QSCT
!        ------------
          Q_SCT=Q_EXT-Q_ABS
!        ------------
!        QBAR   (attention : ltq_gg car on travaille en log dans get_qopt)
!        ------------
          tmp=get_qopt(radius,wln,A_gg,B_gg,fmin_gg,ltq_gg)
          Q_BAR=tmp

       END SUBROUTINE getoptcld 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!        REAL FUNCTION get_qopt(r,wln,A,B,fmin,)
!
!        obtention d'une propriete optique pour une particule de taille r a la longueur d'onde wln
!        a partir de la table tq. 
!        Les parametres tq,fmin,A et B definissent la propriete calculee (qext,qabs ou gg)
!
!        ARGUMENTS D'ENTREE :
!              r : Rayon de la particule (en metres !)
!            wln : Longueur d'onde traitee (en metres !)
!             tq : table de la propriete.
!           fmin : parametre pour extrapolation (debut de table)
!              A : parametre (coefficient directeur) pour extrapolation (fin de table)
!              B : parametre (ordonnee a l'origine)  pour extrapolation (fin de table)
!
!        VALEUR DE RETOUR :
!           Propriete optique recherchee a wln pour une taille r. 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       REAL FUNCTION get_qopt(r,wln,A,B,fmin,tq)
         USE optcld
         IMPLICIT NONE
! ------ INPUT
         REAL   ,INTENT(in) :: r,wln,A,B,fmin,tq(npts)
! ------ LOCAL
         REAL               :: wln_t,val
         INTEGER            :: ind,iver

!        initialisation generale 
         iver = 0
         val = 0.
         wln_t = wln * (r0cld/r)

!        Recherche du point le plus proche dans la table
         CALL locate(tq_wln,npts,wln_t,ind)  

!        Interpolation/extrapolation selon l'indice.
         IF (ind.le.0) THEN
           val = fmin
         ELSEIF(ind.ge.npts) THEN
           CALL extrapolemoi(wln_t,A,B,val,.true.)
         ELSE
           CALL interpolemoi(ind,wln_t,ltq_wln,tq,npts,val,iver,.true.)
         ENDIF
         get_qopt=val

       END FUNCTION get_qopt


