!
! $Header$
!
C
C
      SUBROUTINE conf_dat2d( title,lons,lats,xd,yd,xf,yf,champd ,
     ,                           interbar                        )
c
c     Auteur :  P. Le Van

c    Ce s-pr. configure le champ de donnees 2D 'champd' de telle facon que
c       qu'on ait     - pi    a    pi    en longitude
c       et qu'on ait   pi/2.  a - pi/2.  en latitude
c
c      xd et yd  sont les longitudes et latitudes initiales
c      xf et yf  sont les longitudes et latitudes en sortie , eventuellement
c      modifiees pour etre configurees comme ci-dessus .

      IMPLICIT NONE
 
c    ***       Arguments en  entree      ***
      INTEGER lons,lats
      CHARACTER*25 title
      REAL xd(lons),yd(lats)
      LOGICAL interbar
c
c    ***       Arguments en  sortie      ***
      REAL xf(lons),yf(lats)
c
c    ***  Arguments en entree et  sortie ***
      REAL champd(lons,lats)

c   ***     Variables  locales  ***
c
      REAL pi,pis2,depi
      LOGICAL radianlon, invlon ,radianlat, invlat, alloc
      REAL rlatmin,rlatmax,oldxd1
      INTEGER i,j,ip180,ind

      REAL, ALLOCATABLE :: xtemp(:) 
      REAL, ALLOCATABLE :: ytemp(:) 
      REAL, ALLOCATABLE :: champf(:,:)
     
c
c      WRITE(6,*) ' conf_dat2d  pour la variable ', title

      ALLOCATE( xtemp(lons) )
      ALLOCATE( ytemp(lats) )
      ALLOCATE( champf(lons,lats) )

      DO i = 1, lons
       xtemp(i) = xd(i)
      ENDDO
      DO j = 1, lats
       ytemp(j) = yd(j)
      ENDDO

      pi   = 2. * ASIN(1.) 
      pis2 = pi/2.
      depi = 2. * pi

            radianlon = .FALSE.
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


        IF ( invlon )   THEN

           DO j = 1, lats
            DO i = 1,lons
             champf(i,j) = champd(i,j)
            ENDDO
           ENDDO

           DO i = 1 ,lons
            xf(i) = xtemp(i)
           ENDDO
c
c    ***  On tourne les longit.  pour avoir  - pi  a  +  pi  ****
c
           DO i=1,lons
            IF( xf(i).GT. pi )  THEN
            GO TO 88
            ENDIF
           ENDDO

88         CONTINUE
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

c   .....    on tourne les longitudes  pour  champf ....
c
           DO j = 1,lats

             DO i = ip180,lons
              ind  = i-ip180 +1
              champd (ind,j) = champf (i,j)
             ENDDO
   
             DO i= ind +1,lons
              champd (i,j)  = champf (i-ind,j)
             ENDDO

           ENDDO


        ENDIF
c
c    *****   fin  de   IF(invlon)   ****

         IF ( invlat )    THEN

           DO j = 1,lats
            yf(j) = ytemp(j)
           ENDDO

           DO j = 1, lats
             DO i = 1,lons
              champf(i,j) = champd(i,j)
             ENDDO
           ENDDO

           DO j = 1, lats
              ytemp( lats-j+1 ) = yf(j)
              DO i = 1, lons
               champd (i,lats-j+1) = champf (i,j)
              ENDDO
           ENDDO


         ENDIF

c    *****  fin  de  IF(invlat)   ****

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
        DEALLOCATE(champf)

       DO i = 1, lons
        xf(i) = xtemp(i)
       ENDDO
       DO j = 1, lats
        yf(j) = ytemp(j)
       ENDDO

      deallocate(xtemp)
      deallocate(ytemp)

      RETURN
      END
