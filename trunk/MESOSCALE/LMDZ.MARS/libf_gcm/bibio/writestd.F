      SUBROUTINE writestd(champ,icount)

c=======================================================================
c
c   Auteurs:   Jan Polcher, Frederic Hourdin
c   --------
c
c   Objet:
c   ------
c
c Subroutine qui ecrit les fichiers dans le format de transfert
c vers le LMD.
c format etabli par :
c   Zhao-Xin Li, Frederic Hourdin, Emmanuelle Cohen-Solal,
c   Jan Polcher .
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
c champ: REAL (icount)
c        champ contenant les valeurs et de longeur icount
c icount: INTEGER
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

      IMPLICIT NONE
#include "lmdstd.h"

c   Arguments:
c   ----------

      INTEGER icount
      REAL champ(icount)

c   Local:
c   ------

      INTEGER i,imin,imax
      LOGICAL check, writeornot
      INTEGER lnblnk,ismin,ismax
      EXTERNAL lnblnk,ismin,ismax, isitor, to36

      INTEGER indexfil
      COMMON/indexfil/indexfil

c-----------------------------------------------------------------------
c   verification que le champ est dans la liste predefinie:
c   -------------------------------------------------------

      check = .FALSE.

      IF (gnbfich .EQ. 0) THEN
        writeornot = .TRUE.
      ELSE
        CALL isitor(gfichier,writeornot,gnbfich, glistfich)
      ENDIF

      IF (check) THEN
        WRITE(*,*) 'MYWRITE:: writeornot=',writeornot
      ENDIF

c-----------------------------------------------------------------------
      IF (writeornot) THEN
c-----------------------------------------------------------------------

      imin = ismin(icount,champ,1)
      imax = ismax(icount,champ,1)
      IF (check) THEN
        WRITE(*,*) gentete(1:lnblnk(gentete))
        WRITE(*,*) gtitre(1:lnblnk(gtitre))
        WRITE(*,*) gdatedep
        WRITE(*,*) gdatefin
        WRITE(*,*) gunites(1:lnblnk(gunites))
        WRITE(*,*) gdeltapas, gnbetats
        WRITE(*,*) gdeltajour
        WRITE(*,*) champ(imin),champ(imax)
        WRITE(*,*) gformat(1:lnblnk(gformat))
        WRITE(*,*) icount
      ENDIF
      WRITE(*,*) gtitre(1:lnblnk(gtitre)),'min,max:',
     &  champ(imin),champ(imax),'  ',gunites(1:lnblnk(gunites))
      
c-----------------------------------------------------------------------
c   creation eventuelle d'un nom pour le fichier d'ecriture du champ:
c   -----------------------------------------------------------------

      IF(lnblnk(gfichier).LT.1.OR.lnblnk(gfichier).GE.100) THEN
         gfichier(1:5)='champ'
         IF(indexfil.LT.1) indexfil=49
         gfichier(6:6)=char(indexfil)
	 indexfil=indexfil+1
	 PRINT*,'nom du fichier  ',gfichier
      ENDIF

c-----------------------------------------------------------------------
c   ecriture sur le fichier de sorties:
c   -----------------------------------

      OPEN (unit=55, file=gfichier(1:lnblnk(gfichier)))
      WRITE(55,'(a)') gentete(1:lnblnk(gentete))
      WRITE(55,'(a)') gtitre(1:lnblnk(gtitre))
      WRITE(55,'(a)') gdatedep
      WRITE(55,'(a)') gdatefin
      WRITE(55,'(a)') gunites(1:lnblnk(gunites))
      WRITE(55,*) gdeltapas, gnbetats
      WRITE(55,*) gdeltajour
      WRITE(55,*) champ(imin),champ(imax)
      WRITE(55,'(a)') '#'//gformat(1:lnblnk(gformat))
      WRITE(55,*) icount

      IF (check) THEN
        WRITE(*,'('//gformat(1:lnblnk(gformat))//')') champ(1)
      ENDIF

      IF (lnblnk(gformat).LE.1.OR.lnblnk(gformat).GE.20) THEN
         DO 2230 i=1,icount
            WRITE(55,*) champ(i)
2230     CONTINUE
      ELSE
         IF (gformat(1:lnblnk(gformat)) .EQ. 'zxli1') THEN
            CALL to36 (champ,icount,champ(imin),champ(imax),55, 1)
         ELSE 
            IF (gformat(1:lnblnk(gformat)) .EQ. 'zxli2') THEN
               CALL to36(champ,icount,champ(imin),champ(imax),55,2)
            ELSE
               IF (gformat(1:lnblnk(gformat)) .EQ. 'zxli3') THEN
                  CALL to36(champ,icount,champ(imin),champ(imax),55,3)
               ELSE
                  WRITE(55,'('//gformat(1:lnblnk(gformat))//')') 
     &            (champ(i),i=1,icount)
               ENDIF
            ENDIF
         ENDIF
      ENDIF

      CLOSE (unit=55)

      gfichier=' '

c-----------------------------------------------------------------------
c   si on ecrit pas le champ (writeornot.EQ..f.):
c   ---------------------------------------------

      ELSE
        WRITE(*,*) 'FILE    ',gfichier,'WAS NOT WRITEN'
      ENDIF
      
      RETURN
      END

