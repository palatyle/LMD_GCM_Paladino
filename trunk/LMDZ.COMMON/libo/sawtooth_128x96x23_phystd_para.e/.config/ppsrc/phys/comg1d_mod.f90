










MODULE comg1d_mod
IMPLICIT NONE  

!.......................................................................
!  le COMMON pour GRADS-1D
!  (Utilise pour les sorties format Grads dans la version 1D du modele)
!
!  on peut se dire : "on ne sauvera pas plus de 1000 variables ... hein ?"
!
      INTEGER g1d_nvarmx
      PARAMETER(g1d_nvarmx=1000)

!         * g1d_nlayer     ---> nombre de couches verticales
!         * g1d_nomfich    ---> nom du fichier grads
!         * g1d_unitfich   ---> code du fichier grads
!         * g1d_nomctl     ---> nom du fichier ctl
!         * g1d_unitctl    ---> code du fichier ctl
!         * g1d_premier    ---> variable logique pour dire si le fichier
!                               est deja ouvert
!         * g1d_irec       ---> indice de derniere ecriture
!         * g1d_nvar       ---> nombre de variables deja definies a la
!                               derniere ecriture
!         * g1d_nomvar     ---> noms des vecteurs existants
!         * g1d_dimvar     ---> taille des vecteurs
!         * g1d_titrevar   ---> titres des vecteurs
!         * g1d_tmp1       ---> caractere 
!         * g1d_tmp2       ---> caractere 

      INTEGER g1d_nlayer
      CHARACTER*100 g1d_nomfich
      INTEGER g1d_unitfich
      CHARACTER*100 g1d_nomctl
      INTEGER g1d_unitctl
      LOGICAL g1d_premier
      LOGICAL g2d_premier
      INTEGER g1d_irec
      INTEGER g2d_irec
      INTEGER g2d_appel
      INTEGER g1d_nvar
      CHARACTER*100 g1d_nomvar(0:g1d_nvarmx)
      INTEGER g1d_dimvar(0:g1d_nvarmx)
      CHARACTER*100 g1d_titrevar(0:g1d_nvarmx)
      CHARACTER*100 g1d_tmp1,g1d_tmp2

!     added by RW
      integer saveG1D
!
!.......................................................................

END MODULE comg1d_mod
