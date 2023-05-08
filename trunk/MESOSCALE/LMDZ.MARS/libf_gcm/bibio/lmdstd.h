c=======================================================================
c   INCLUDE '/usr/local/lmdgraph/libio/lmdstd.h'

      INTEGER bandmax
      PARAMETER(bandmax=24)

      COMMON/lmdstdch/gtitre,gunites,gformat,gfichier,glistfich(100),
     $                gentete,gdatedep,gdatefin
        
      CHARACTER*100 gtitre    !du champ
      CHARACTER*20  gunites   !unites du champ
      CHARACTER*20  gformat   !formats FORTRAN ou zxliN (N=1,2,3) ou ''
      CHARACTER*100 gfichier  !nom du fichier contenant le champ
      CHARACTER*100 gentete   !nom generique (exemple desript. du RUN)
      CHARACTER*8   gdatedep  !date de depart de la moy. ou evol. temp.
      CHARACTER*8   gdatefin  !date de fin '' si une seule datefin
      CHARACTER*20  glistfich !liste de champs a sortir

      COMMON/lmdstdnu/ gminimum,gmaximum,
     $                 gdeltajour(bandmax),gdeltapas(bandmax),
     $                 gnbetats(bandmax),gnbfich

      REAL    gminimum        ! \    min et max
      REAL    gmaximum        ! /    du champ
      INTEGER gnbfich         ! nombre de fichier a sortir
c  pour les evolutions temporelles:
      INTEGER gnbetats        ! nombre d'etats par bande
      INTEGER gdeltapas       ! nombre de pas de temps entre 2 sorties
      REAL    gdeltajour      ! ecart en jour des sorties sur bandes

c=======================================================================
