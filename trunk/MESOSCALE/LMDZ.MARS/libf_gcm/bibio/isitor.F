
      SUBROUTINE isitor(filename,writeornot,filenum, files)

C--------------------------------------------------------------------
C Subroutine qui determine si un fichier doit etre ecrit ou non
C--------------------------------------------------------------------
C     ARGUMENTS
C     +++++++++
C     filename: CHARACTER*20 
C               nom du fichier dont on veut savoir si il doit etre ecrit
C     writeornot: LOGICAL
C                 Variable de sortie qui est vraie si le champ doit etre ecrit
C     filenum: INTEGER
C              nombre de fichier dans la list files
C     files: CHARACTER*20 files(100)
C            c'est la liste des fichier qui doivent etre ecrit et dans
C            lequel on vas voir si filename y est.
C


      IMPLICIT NONE
      CHARACTER*20 filename, checkname
      CHARACTER*1 separator
      LOGICAL writeornot, check
      INTEGER iletter, ifile
      INTEGER filenum
      CHARACTER*20 files(100)

      writeornot = .FALSE.
      check = .FALSE.
      separator = '/'

      DO 6565 iletter = 20,1,-1

        IF (filename(iletter:iletter) .EQ. separator) THEN

          checkname = filename(iletter+1:20)
    
        ELSE
        ENDIF
 6565 CONTINUE

      IF (check) THEN
        WRITE(*,*) 'ISITORNOT:: filename=',filename,
     &             'checkname=',checkname
      ENDIF

      DO 6767 ifile = 1,filenum

      IF (files(ifile) .EQ. checkname ) THEN
       
      writeornot = .TRUE. 

      ENDIF

      IF (check) THEN
        WRITE(*,*) 'ISITORNOT:: writeornot=',writeornot
      ENDIF

 6767 CONTINUE

      RETURN
      END

