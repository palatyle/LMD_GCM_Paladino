      SUBROUTINE readstd(imx,champ,imax)
      IMPLICIT NONE

c=======================================================================
c
c   Auteurs:   Jan Polcher, Frederic Hourdin  20/04/91
c   --------
c
c   Objet:
c   ------
c
c Subroutine qui lit les fichiers dans le format lmdstd
c format etabli par :
c  Jan Polcher,Zhao-Xin Li, Frederic Hourdin, Emmanuelle Cohen-Solal,
c Il peut etre utilise pour sortire  des evolutions
c temporelles, des champs 2D (iim,jjm ou jjm,llm) ou des champs 3D.
c le character # place devant le format determine la fin de la 
c description du champ et permet eventulellement de rajouter
c des informations (avant le #).
c le programme de lecture lit , dans la description, les variables
c qu'il connait et saute apres jusquau #. A partir de la
c il commence a lire le champ dans le format specifie.
c
c-----------------------------------------------------------------------
c
c   Interface:
c   ----------
c
c   ARGUMENTS
c   ---------
c
c imx  : INTEGER
c        taille du tableau contenant le champ
c champ: REAL (imx)
c        champ contenant les valeurs et de longeur icount
c imax : INTEGER
c        nombre d'element dans le champ
c
c   COMMON:
c   -------
c
c=======================================================================
c  INCLUDE '/usr/local/lmdgraph/libio/lmdstd.h'
c
c     INTEGER bandmax
c     PARAMETER(bandmax=24)
c
c     COMMON/lmdstdch/gtitre,gunites,gformat,gfichier,glistfich(100),
c    $                gentete,gdatedep,gdatefin
c       
c     CHARACTER*100 gtitre    !du champ
c     CHARACTER*20  gunites   !unites du champ
c     CHARACTER*20  gformat   !formats FORTRAN ou zxliN (N=1,2,3) ou ''
c     CHARACTER*100 gfichier  !nom du fichier contenant le champ
c     CHARACTER*100 gentete   !nom generique (exemple desript. du RUN)
c     CHARACTER*8   gdatedep  !date de depart de la moy. ou evol. temp.
c     CHARACTER*8   gdatefin  !date de fin '' si une seule datefin
c     CHARACTER*20  glistfich !liste de champs a sortir
c
c     COMMON/lmdstdnu/ gminimum,gmaximum,
c    $                 gdeltajour(bandmax),gdeltapas(bandmax),
c    $                 gnbetats(bandmax),gnbfich
c
c     REAL    gminimum        ! \    min et max
c     REAL    gmaximum        ! /    du champ
c     INTEGER gnbfich         ! nombre de fichier a sortir
c pour les evolutions temporelles:
c     INTEGER gnbetats        ! nombre d'etats par bande
c     INTEGER gdeltapas       ! nombre de pas de temps entre 2 sorties
c     REAL    gdeltajour      ! ecart en jour des sorties sur bandes
c
c=======================================================================
c-----------------------------------------------------------------------
c   Declarations:
c   -------------

#include "lmdstd.h"

c   Arguments:
c   ----------

      INTEGER imx,imax
      REAL champ(imx)

c   Local:
c   ------

      INTEGER i

      CHARACTER*200 test(20)
      CHARACTER*200 helptest(20)
      CHARACTER*200 help
      INTEGER reccount

      LOGICAL check,continue

      EXTERNAL lnblnk
      INTEGER lnblnk

c-----------------------------------------------------------------------
c   Initialisations:
c   ----------------

      check = .FALSE.
      DO 9876 i=1,20
        helptest(i) = '        '
 9876 CONTINUE


c-----------------------------------------------------------------------
c   lecture le l'entete:
c   --------------------

 2222 OPEN (unit=10,file=gfichier(1:lnblnk(gfichier)), form='formatted',
     &      ERR = 9999, status='old')
      GOTO 1111

 9999 WRITE(*,*) 'Je ne peut pas lire le fichier: ',
     & gfichier(1:lnblnk(gfichier))
      CLOSE(10)
      WRITE(*,*) 'Entrez S.V.P. un autre nom de fichier (path absolu) '
      WRITE(*,*) 'ou RETURN pour stopper lexecution'
      READ (*,'(a)') gfichier
      IF(lnblnk(gfichier).EQ.0) STOP
      GOTO 2222

 1111 continue = .TRUE.
      DO 2345 i=1,20
        IF (continue) THEN
          READ (10,'(a200)') test(i)
          IF ( test(i)(1:1) .EQ. '#' .OR. test(i)(1:1) .EQ. 'f10' ) THEN
            continue = .FALSE.
            reccount = i
          ENDIF
        ENDIF
 2345 CONTINUE

      DO 1223 i=1,reccount-1
        helptest(i) = test(i)
 1223 CONTINUE
      helptest(20) = test(reccount)


      IF (check) THEN
        WRITE(*,*) 'gentete'
      ENDIF
      IF ( helptest(1) .NE. ' ' ) THEN
        read (helptest(1),'(a100)') gentete
      ELSE
c       WRITE(*,*) 'la description du run'
c       READ (*,'(a100)') gentete
        gentete = '    '
      ENDIF

      IF (check) THEN
        WRITE(*,*) 'titre'
      ENDIF
      IF ( helptest(2) .NE. ' ' ) THEN
        read (helptest(2),'(a100)') gtitre
      ELSE
c       WRITE(*,*) 'le titre du dessin'
c       READ (*,'(a100)') gtitre 
        gtitre = '    '
      ENDIF

      IF (check) THEN
        WRITE(*,*) 'datefins'
      ENDIF
      IF ( helptest(3) .NE. ' ' ) THEN
        read (helptest(3),'(a8)') gdatedep 
        read (helptest(4),'(a8)') gdatefin
      ELSE
c       WRITE(*,*)
c    &  'datefin du debut de la moyenne (ou le mois tout simplement)'
c       READ (*,'(a8)') gdatedep 
c       WRITE(*,*) 'et celle de la fin (ou rien du tout)'
c       READ (*,'(a8)') gdatefin 
        gdatedep = '       '
        gdatefin = '       '
      ENDIF

      IF (check) THEN
        WRITE(*,*) 'unitees'
      ENDIF
      IF ( helptest(5) .NE. ' ' ) THEN
        read (helptest(5),'(a20)') gunites
      ELSE
c       WRITE(*,*) 'les unitees du champ'
c       READ (*,'(a20)') gunites
        gunites = '                    '
      ENDIF

      IF (check) THEN
        WRITE(*,*) 'deltapas et nbetats'
      ENDIF
      IF ( helptest(6) .NE. ' ' ) THEN
c       READ (helptest(6),*,ERR=1212) 
c    .       (gdeltapas(i),i=1,bandmax),(gnbetats(i),i=1,bandmax)
        GOTO 1213

 1212   help =helptest(6)(1:lnblnk(helptest(6)))//
     .        helptest(7)(1:lnblnk(helptest(7)))
        helptest(6) = help
        IF (check) THEN
          WRITE(*,*) 'Problems while reading', help
        ENDIF
        DO  1214 i=7,reccount
          helptest(i) = helptest(i+1)
 1214   CONTINUE
c       READ (helptest(6),*,ERR=1212) 
c    .       (gdeltapas(i),i=1,bandmax),(gnbetats(i),i=1,bandmax)

 1213   CONTINUE
      ENDIF

      IF (check) THEN
        WRITE(*,*) 'deltajour'
      ENDIF
      IF ( helptest(7) .NE. ' ' ) THEN
c       READ (helptest(7),*,ERR=1312) (gdeltajour(i),i=1,bandmax)
        GOTO 1313

 1312   help =helptest(7)(1:lnblnk(helptest(7)))//
     .        helptest(8)(1:lnblnk(helptest(8)))
        helptest(7) = help
        IF (check) THEN
          WRITE(*,*) 'Problems while reading', help
        ENDIF
        DO  1314 i=8,reccount
          helptest(i) = helptest(i+1)
 1314   CONTINUE
c       READ (helptest(7),*,ERR=1312) (gdeltajour(i),i=1,bandmax)

 1313   CONTINUE
      ENDIF

      IF (check) THEN
        WRITE(*,*) 'min et max'
      ENDIF
      IF ( helptest(8) .NE. ' ' ) THEN
c       READ (helptest(8),*)     gminimum,gmaximum
      ENDIF

      gformat = helptest(20)(2:20)

      IF (check) THEN
        WRITE(*,*) gentete
        WRITE(*,*) gtitre
        WRITE(*,*) gdatedep
        WRITE(*,*) gdatefin 
        WRITE(*,*) gdeltapas,gnbetats
        WRITE(*,*) gdeltajour
        WRITE(*,*) gminimum,gmaximum
        WRITE(*,*) gformat
      ENDIF

c-----------------------------------------------------------------------
c   lecture du champ proprement dit:
c   --------------------------------

      READ (10,*)     imax
   
      IF(imax.GT.imx) THEN
         WRITE(*,*) 'Vous essayez de lire un champ dont la taille est'
         WRITE(*,*) 'superieure a la taille maximum autorisee'
         WRITE(*,*) 'c est a dire a la dimension du champ utilise'
         WRITE(*,*) 'pour lappel a la procedure writestd.'
         STOP
      ENDIF

      IF (gformat(1:4) .EQ. 'zxli') THEN
        CALL from36 (champ,imax,10)
      ELSE 
         IF (lnblnk(gformat).LE.1.OR.lnblnk(gformat).GE.20) THEN
            READ(10,*) (champ (i),i=1,imax)
         ELSE
            READ(10,'('//gformat(1:lnblnk(gformat))//')')
     $      (champ (i),i=1,imax)
         ENDIF
      ENDIF

      CLOSE (unit=10)


      RETURN
      END
