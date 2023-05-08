      SUBROUTINE readcoord(imx,jmx,lmx,imax,jmax,lmax,x,y,z,pathcoor)
      IMPLICIT NONE

c=======================================================================
c
c   Auteurs Jan Polcher, Frederic Hourdin
c   -------
c
c   Objet   Lecture du fichier contenant des coordonnees x,y et z.
c   -----
c
c   Inteface
c   --------
c
c   Arguments:
c   ----------
c
c   input
c   -----
c      imx   \
c      jmx          dimensions des coordonnees
c      lmx   /
c      pathcoor     path absolu pour le fichier des coordonnees
c                   CHARACTER*80
c
c   ouput
c   -----
c      imax         nombre de coordonnees x effectives
c      jmax         nombre de coordonnees y effectives
c      lmax         nombre de coordonnees z effectives
c      x(imx)
c      y(jmx)
c      z(lmx)
c
c=======================================================================
c-----------------------------------------------------------------------
c   declarations:
c   -------------

c   arguments:
c   ----------

      INTEGER imx,jmx,lmx
      INTEGER imax,jmax,lmax

      REAL x(imx),y(jmx),z(lmx)

      CHARACTER*80 pathcoor

c   local
c   -----

      INTEGER i,j,l

      CHARACTER*100 filecoor,messerr

      EXTERNAL lnblnk
      INTEGER lnblnk

c-----------------------------------------------------------------------
c   Initialisations:
c   ----------------

      messerr=' '
      filecoor='coord.def'

c-----------------------------------------------------------------------

100   IF(pathcoor(1:1).EQ.'/') THEN
         filecoor= pathcoor(:lnblnk(pathcoor))//filecoor
      ENDIF
      OPEN (99,file=filecoor,err=200,status='old',form='formatted')
      GO TO 300

200   WRITE(*,*) messerr
      WRITE(*,*) 'entrer le nom du fichier des coordonnees puis RETURN'
      WRITE(*,*) 'ou seulement RETURN pour arreter le programme'
      READ(*,'(a)') filecoor
      IF(filecoor(1:1).EQ.' ') STOP
      CLOSE (99)
      GOTO 100

300   CONTINUE

      READ(99,*,err=9999) imax
      IF(imax.GT.imx) THEN
        WRITE(6,1000) 'x',imx,imax
        STOP
      ENDIF
      READ(99,*) (x(i),i=1,imax)

      READ(99,*,err=9999) jmax
      IF(jmax.GT.jmx) THEN
        WRITE(6,1000) 'y',jmx,jmax
        STOP
      ENDIF
      READ(99,*,err=9999) (y(j),j=1,jmax)

      READ(99,*,err=9999) lmax
      IF(lmax.GT.lmx) THEN
        WRITE(6,1000) 'z',lmx,lmax
        STOP
      ENDIF
      READ(99,*) (z(l),l=1,lmax)

      CLOSE(99)

c-----------------------------------------------------------------------

      RETURN

9999  PRINT*,'Le fichier ',filecoor,' que vous lisez est mal structure'
      PRINT*,'Il doit contenir :'
      PRINT*,'  le nombre de points en x suivi des valeurs de x'
      PRINT*,'  puis la meme chose pour y et z'
      STOP

1000  FORMAT(/,
     $ 'Probleme avec le fichier de coordonnees',/,
     $ 'Le tableau contenant ',a3,' est dimenssionne a ',i4,/,
     $ 'et le fichier de coordonnees contient ',i4,' valeurs.')

      END
