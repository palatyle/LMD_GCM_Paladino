c.......................................................................
c  le COMMON pour GRADS-1D
c  (Utilise pour les sorties format Grads dans la version 1D du modele)
c
c  on peut se dire : "on ne sauvera pas plus de 1000 variables ... hein ?"
c
      INTEGER g1d_nvarmx
      PARAMETER(g1d_nvarmx=1000)
c
c         * g1d_nlayer     ---> nombre de couches verticales
c         * g1d_nomfich    ---> nom du fichier grads
c         * g1d_unitfich   ---> code du fichier grads
c         * g1d_nomctl     ---> nom du fichier ctl
c         * g1d_unitctl    ---> code du fichier ctl
c         * g1d_premier    ---> variable logique pour dire si le fichier
c                               est deja ouvert
c         * g1d_irec       ---> indice de derniere ecriture
c         * g1d_nvar       ---> nombre de variables deja definies a la
c                               derniere ecriture
c         * g1d_nomvar     ---> noms des vecteurs existants
c         * g1d_dimvar     ---> taille des vecteurs
c         * g1d_titrevar   ---> titres des vecteurs
c         * g1d_tmp1       ---> caractere 
c         * g1d_tmp2       ---> caractere 
c
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
      CHARACTER*100 g1d_nomvar
      INTEGER g1d_dimvar
      CHARACTER*100 g1d_titrevar
      CHARACTER*100 g1d_tmp1,g1d_tmp2
c
      COMMON/COMG1DI/g1d_nlayer
     &             ,g1d_unitfich
     &             ,g1d_unitctl
     &             ,g1d_irec
     &             ,g2d_irec
     &             ,g2d_appel
     &             ,g1d_nvar
      COMMON/COMG1DC/g1d_dimvar(0:g1d_nvarmx)
     &             ,g1d_nomfich
     &             ,g1d_nomctl
     &             ,g1d_nomvar(0:g1d_nvarmx)
     &             ,g1d_titrevar(0:g1d_nvarmx)
     &             ,g1d_tmp1
     &             ,g1d_tmp2
      COMMON/COMG1DL/g1d_premier
     &             ,g2d_premier
c
c.......................................................................
