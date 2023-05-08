      SUBROUTINE mywrite(champ,icount,titre,unit,format,nom
     &                  , entete, date0, date, deltpas, nbeta
     &                  , deltajour, filenum, files)
C----------------------------------------------------------
C Subroutine qui ecrit les fichiers dans le format de transfert
C vers le LMD.
C format etabli par Zhao-Xin Li, Frederique Hourdin, Emmanuelle Cohen-Solal,
C Jan Polcher . Il peut etre utilise pour sortire  des evolutions
C temporelles, des champs 2D (iim,jjm ou jjm,llm) ou des champs 3D.
C le character # place devant le format determine la fin de la 
C description du champ et permet eventulellement de rajouter
C des informations (avant le #).
C le programme de lecture lit , dans la description, les variables
C qu'il connait et saute apres jusquau #. A partir de la
C il commence a lire le champ dans le format specifie.
c------------------------------------------------------------
C       ARGUMENTS
C       +++++++++
C champ: REAL (icount)
C        champ contenant les valeurs et de longeur icount
C icount: INTEGER
C         nombre d'element dans le champ
C titre: CHARACTER*100
C        nom d'ecrivant le champ
C unit: CHARACTER*20
C       unites du champ
C format: CHARACTER*20
C         format dans le quel le champ doit etre ecrit.
C         au choix:
C                  - tous les format FORTRAN
C                  - compactage de L.Li (zxliN) ou N peut valoir 1, 2 ou 3
C nom: CHARACTER*20
C      nom du fichier dans le quel le champ sera ecris, eventuellement
C      avec directory (ex: evol/tsol)
C entete: CHARACTER*100
C         une entete du fichier pouvant par exemple decrire le run
C date0: CHARACTER*8
C        date de depart de la moyenne ou de l'evolution temporelle que contient
C        le fichier
C date : CHARACTER*8
C        date de fin de la moyenne ou de l'evolution temporelle que contient
C        le fichier
C deltpas: INTEGER (bandmax)
C          decript le nombre de pas de temps entre deux sorties pour chamque
C          bande (necesaire pour les evolutions temporelles)
C nbeta: INTEGER (bandmax)
C        donne le nombre d'etas par bande (necesaire pour les evolutions 
C        temporelles)
C deltajour: REAL (bandmax)
C             l'ecart en jours des iscriptions sur chaque bande (necesaire 
C             pour les evolutions temporelles)
C filenum: INTEGER
C          nombre d'elements dans la liste files des champs desires
C          si filenum = 0 alors tous les champs sont sortis  
C files: CHARACTER*20 (100)
C        list des champ que l'on desire sortire
C
C        PARAMETER
C        +++++++++
C bandmax: INTEGER
C          le nombre maximale que peut traiter le programme
C


      IMPLICIT NONE

      INTEGER bandmax
      PARAMETER (bandmax = 24)
      INTEGER icount,i,imin,imax
      REAL champ(icount)
      CHARACTER*100 titre, entete
      CHARACTER unit*20,nom*80
      CHARACTER*20 format
      CHARACTER*8 date0, date
      INTEGER deltpas(bandmax), nbeta(bandmax)
      REAL deltajour(bandmax)
      LOGICAL check, writeornot
      INTEGER lnblnk,ismin,ismax
      INTEGER filenum
      CHARACTER*20 files(100)
      EXTERNAL lnblnk,ismin,ismax, isitor, to36

      check = .FALSE.

      IF (filenum .EQ. 0) THEN
        writeornot = .TRUE.
      ELSE
        CALL isitor(nom,writeornot,filenum, files)
      ENDIF

      IF (check) THEN
        WRITE(*,*) 'MYWRITE:: writeornot=',writeornot
      ENDIF

      IF (writeornot) THEN
      imin = ismin(icount,champ,1)
      imax = ismax(icount,champ,1)
      IF (check) THEN
        WRITE(*,*) entete(1:lnblnk(entete))
        WRITE(*,*) titre(1:lnblnk(titre))
        WRITE(*,*) date0
        WRITE(*,*) date
        WRITE(*,*) unit(1:lnblnk(unit))
        WRITE(*,*) deltpas, nbeta
        WRITE(*,*) deltajour
        WRITE(*,*) champ(imin),champ(imax)
        WRITE(*,*) format(1:lnblnk(format))
        WRITE(*,*) icount
      ENDIF
      WRITE(*,*) titre(1:lnblnk(titre)),'min,max:',
     &           champ(imin),champ(imax),'  ',unit(1:lnblnk(unit))
      
      OPEN (unit=55, file=nom(1:lnblnk(nom)))
      WRITE(55,'(a)') entete(1:lnblnk(entete))
      WRITE(55,'(a)') titre(1:lnblnk(titre))
      WRITE(55,'(a)') date0
      WRITE(55,'(a)') date
      WRITE(55,'(a)') unit(1:lnblnk(unit))
      WRITE(55,*) deltpas, nbeta
      WRITE(55,*) deltajour
      WRITE(55,*) champ(imin),champ(imax)
      WRITE(55,'(a)') '#'//format(1:lnblnk(format))
      WRITE(55,*) icount

      IF (check) THEN
        WRITE(*,'('//format(1:lnblnk(format))//')') champ(1)
      ENDIF

      IF (format(1:lnblnk(format)) .EQ. 'zxli1') THEN
        CALL to36 (champ,icount,champ(imin),champ(imax),55, 1)
      ELSE 
        IF (format(1:lnblnk(format)) .EQ. 'zxli2') THEN
          CALL to36 (champ, icount,champ(imin), champ(imax),55, 2)
        ELSE
          IF (format(1:lnblnk(format)) .EQ. 'zxli3') THEN
            CALL to36 (champ,icount,champ(imin),champ(imax),55, 3)
          ELSE
            DO 2233 i=1,icount
              WRITE(55,'('//format(1:lnblnk(format))//')') 
     &             champ(i)
 2233       CONTINUE
          ENDIF
        ENDIF
      ENDIF

      CLOSE (unit=55)

      ELSE
        WRITE(*,*) 'FILE    ',nom,'WAS NOT WRITEN'
      ENDIF
      
      RETURN
      END

