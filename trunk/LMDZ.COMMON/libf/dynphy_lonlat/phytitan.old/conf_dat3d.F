!
! $Header$
!
C
C
      SUBROUTINE conf_dat3d( title, lons,lats,levs,xd,yd,zd,xf,yf,zf,
     ,                                 champd , interbar             )
c
c     Auteur : P. Le Van
c
c    Ce s-pr. configure le champ de donnees 3D 'champd' de telle facon 
c       qu'on ait     - pi    a    pi    en longitude
c       qu'on ait      pi/2.  a - pi/2.  en latitude
c      et qu'on ait les niveaux verticaux variant du sol vers le ht de l'atmos.
c           (     en Pascals   ) .
c
c      xd et yd  sont les longitudes et latitudes initiales
c      zd  les pressions initiales
c
c      xf et yf  sont les longitudes et latitudes en sortie , eventuellement
c       modifiees pour etre configurees comme ci-dessus .
c      zf  les pressions en sortie
c
c      champd   en meme temps le champ initial et  final
c
c      interbar = .TRUE.  si on appelle l'interpo. barycentrique inter_barxy
c          sinon , l'interpolation   grille_m  ( grid_atob ) .
c

      IMPLICIT NONE
 
c    ***       Arguments en  entree      ***
      CHARACTER*(*) :: title
      INTEGER lons, lats, levs
      REAL xd(lons), yd(lats), zd(levs)
      LOGICAL interbar
c
c    ***       Arguments en  sortie      ***
      REAL xf(lons), yf(lats), zf(levs)

c    ***  Arguments en entree et  sortie ***
      REAL  champd(lons,lats,levs)

c    ***  Variables locales  ***
c
      REAL pi,pis2,depi,presmax
      LOGICAL radianlon, invlon ,radianlat, invlat, invlev, alloc
      REAL rlatmin,rlatmax,oldxd1
      INTEGER i,j,ip180,ind,l

      REAL, ALLOCATABLE :: xtemp(:)
      REAL, ALLOCATABLE :: ytemp(:)
      REAL, ALLOCATABLE :: ztemp(:)
      REAL, ALLOCATABLE :: champf(:,:,:)
     

c      WRITE(6,*) '  Conf_dat3d  pour  ',title

      ALLOCATE(xtemp(lons))
      ALLOCATE(ytemp(lats))
      ALLOCATE(ztemp(levs))

      DO i = 1, lons
       xtemp(i) = xd(i)
      ENDDO
      DO j = 1, lats
       ytemp(j) = yd(j)
      ENDDO
      DO l = 1, levs
       ztemp(l) = zd(l)
      ENDDO

      pi   = 2. * ASIN(1.) 
      pis2 = pi/2.
      depi = 2. * pi

      IF( xtemp(1).GE.-pi-0.5.AND. xtemp(lons).LE.pi+0.5 )  THEN
            radianlon = .TRUE.
            invlon    = .FALSE.
      ELSE IF (xtemp(1).GE.-0.5.AND.xtemp(lons).LE.depi+0.5 ) THEN
            radianlon = .TRUE.
            invlon    = .TRUE.
      ELSE IF ( xtemp(1).GE.-180.5.AND. xtemp(lons).LE.180.5 )   THEN
            radianlon = .FALSE.
            invlon    = .FALSE.
      ELSE IF ( xtemp(1).GE.-0.5.AND.xtemp(lons).LE.360.5 )   THEN
            radianlon = .FALSE.
            invlon    = .TRUE.
      ELSE
        WRITE(6,*) 'Pbs. sur les longitudes des donnees pour le fichier'
     ,  , title
      ENDIF

      invlat = .FALSE.
      
      IF( ytemp(1).LT.ytemp(lats) ) THEN
        invlat = .TRUE.
      ENDIF

      rlatmin = MIN( ytemp(1), ytemp(lats) )
      rlatmax = MAX( ytemp(1), ytemp(lats) )
      
      IF( rlatmin.GE.-pis2-0.5.AND.rlatmax.LE.pis2+0.5)THEN
             radianlat = .TRUE.
      ELSE IF ( rlatmin.GE.-90.-0.5.AND.rlatmax.LE.90.+0.5 ) THEN
             radianlat = .FALSE.
      ELSE
        WRITE(6,*) ' Pbs. sur les latitudes des donnees pour le fichier'
     ,  , title
      ENDIF

       IF( .NOT. radianlon )  THEN
         DO i = 1, lons
          xtemp(i) = xtemp(i) * pi/180.
         ENDDO
       ENDIF

       IF( .NOT. radianlat )  THEN
         DO j = 1, lats
          ytemp(j) = ytemp(j) * pi/180.
         ENDDO   
       ENDIF


        alloc =.FALSE.

        IF ( invlon )   THEN

            ALLOCATE(champf(lons,lats,levs))
            alloc = .TRUE.

            DO i = 1 ,lons
             xf(i) = xtemp(i)
            ENDDO

            DO l = 1, levs
             DO j = 1, lats
              DO i= 1, lons
               champf (i,j,l)  = champd (i,j,l)
              ENDDO
             ENDDO
            ENDDO
c
c    ***  On tourne les longit.  pour avoir  - pi  a  +  pi  ****
c
            DO i=1,lons
             IF( xf(i).GT. pi )  THEN
              GO TO 88
             ENDIF
            ENDDO

88          CONTINUE
c
            ip180 = i

            DO i = 1,lons
             IF (xf(i).GT. pi)  THEN
              xf(i) = xf(i) - depi
             ENDIF
            ENDDO

            DO i= ip180,lons
             ind = i-ip180 +1
             xtemp(ind) = xf(i)
            ENDDO

            DO i= ind +1,lons
             xtemp(i) = xf(i-ind)
            ENDDO

c   .....    on tourne les longitudes  pour champf  ....
c
            DO l = 1,levs
              DO j = 1,lats
               DO i = ip180,lons
                ind  = i-ip180 +1
                champd (ind,j,l) = champf (i,j,l)
               ENDDO
   
               DO i= ind +1,lons
                champd (i,j,l)  = champf (i-ind,j,l)
               ENDDO
              ENDDO
            ENDDO

        ENDIF
c
c    *****   fin  de   IF(invlon)   ****
         
         IF ( invlat )    THEN

           IF(.NOT.alloc)  THEN 
            ALLOCATE(champf(lons,lats,levs))
            alloc = .TRUE.
           ENDIF

           DO j = 1, lats
            yf(j) = ytemp(j)
           ENDDO
         
           DO l = 1,levs
            DO j = 1, lats
             DO i = 1,lons
              champf(i,j,l) = champd(i,j,l)
             ENDDO
            ENDDO

            DO j = 1, lats
              ytemp( lats-j+1 ) = yf(j)
              DO i = 1, lons
               champd (i,lats-j+1,l) = champf (i,j,l)
              ENDDO
            ENDDO
          ENDDO


         ENDIF

c    *****  fin  de  IF(invlat)   ****
c
c
      IF( interbar )  THEN
        oldxd1 = xtemp(1)
        DO i = 1, lons -1
          xtemp(i) = 0.5 * ( xtemp(i) + xtemp(i+1) )
        ENDDO
          xtemp(lons) = 0.5 * ( xtemp(lons) + oldxd1 + depi )

        DO j = 1, lats -1
          ytemp(j) = 0.5 * ( ytemp(j) + ytemp(j+1) )
        ENDDO
      ENDIF
c

      invlev = .FALSE.
      IF( ztemp(1).LT.ztemp(levs) )  invlev = .TRUE.

      presmax = MAX( ztemp(1), ztemp(levs) )
      IF( presmax.LT.1200. ) THEN
         DO l = 1,levs
           ztemp(l) = ztemp(l) * 100.
         ENDDO
      ENDIF

      IF( invlev )  THEN

          IF(.NOT.alloc)  THEN
            ALLOCATE(champf(lons,lats,levs))
            alloc = .TRUE.
          ENDIF

          DO l = 1,levs
            zf(l) = ztemp(l)
          ENDDO

          DO l = 1,levs
            DO j = 1, lats
             DO i = 1,lons
              champf(i,j,l) = champd(i,j,l)
             ENDDO
            ENDDO
          ENDDO

          DO l = 1,levs
            ztemp(levs+1-l) = zf(l)
          ENDDO

          DO l = 1,levs
            DO j = 1, lats
             DO i = 1,lons
              champd(i,j,levs+1-l) = champf(i,j,l)
             ENDDO
            ENDDO
          ENDDO


      ENDIF

         IF(alloc)  DEALLOCATE(champf)

         DO i = 1, lons
           xf(i) = xtemp(i)
         ENDDO
         DO j = 1, lats
           yf(j) = ytemp(j)
         ENDDO
         DO l = 1, levs
           zf(l) = ztemp(l)
         ENDDO

      DEALLOCATE(xtemp)
      DEALLOCATE(ytemp)
      DEALLOCATE(ztemp)

      RETURN
      END
