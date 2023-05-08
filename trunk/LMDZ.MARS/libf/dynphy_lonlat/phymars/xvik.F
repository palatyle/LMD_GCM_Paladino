      PROGRAM xvik

      USE filtreg_mod, ONLY: inifilr
      USE comconst_mod, ONLY: dtvr,g,r,pi

      IMPLICIT NONE
c=======================================================================
c
c  Pression au site Viking
c
c=======================================================================
c-----------------------------------------------------------------------
c   declarations:
c   -------------


      include "dimensions.h"
      include "paramet.h"
      include "comdissip.h"
      include "comgeom2.h"
!#include "control.h"
      include "netcdf.inc"      


      INTEGER itau,nbpas,nbpasmx
      PARAMETER(nbpasmx=1000000)
      REAL temps(nbpasmx)
      INTEGER unitlec
      INTEGER i,j,l,jj
      REAL constR

c   Declarations NCDF:
c   -----------------
      CHARACTER*100  varname
      INTEGER ierr,nid,nvarid,dimid
      LOGICAL nc
      INTEGER start_ps(3),start_temp(4),start_co2ice(3)
      INTEGER count_ps(3),count_temp(4),count_co2ice(3)

c   declarations pour les points viking:
c   ------------------------------------
      INTEGER ivik(2),jvik(2),ifile(2),iv
      REAL lonvik(2),latvik(2),phivik(2),phisim(2)
      REAL unanj

c   variables meteo:
c   ----------------
      REAL vnat(iip1,jjm,llm),unat(iip1,jjp1,llm)
      REAL t(iip1,jjp1,llm),ps(iip1,jjp1),pstot, phis(iip1,jjp1)
      REAL co2ice(iip1,jjp1), captotN,captotS
      real t7(iip1,jjp1) ! temperature in 7th atmospheric layer

      REAL zp1,zp2,zp2_sm,zu,zv,zw(0:1,0:1,2),zalpha,zbeta

      LOGICAL firstcal,lcal,latcal,lvent,day_ls
      INTEGER*4 day0

      REAL ziceco2(iip1,jjp1)
      REAL day,zt,sollong,sol,dayw
      REAL airtot1,gh

      INTEGER ii,iyear,kyear

      CHARACTER*2 chr2

       
c   declarations de l'interface avec mywrite:
c   -----------------------------------------

      CHARACTER file*80
      CHARACTER pathchmp*80,pathsor*80,nomfich*80

c   externe:
c   --------

      EXTERNAL iniconst,inigeom,covcont,mywrite
      EXTERNAL exner,pbar
      EXTERNAL solarlong,coordij,moy2
      EXTERNAL SSUM
      REAL SSUM

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c-----------------------------------------------------------------------
c   initialisations:
c   ----------------

      chr2="0"
      unanj=669.
      print*,'WARNING!!!',unanj,'Jours/an'
      nc=.true.
      lcal=.true.
      latcal=.true.
      lvent=.false.
      day_ls=.true.

c lecture du fichier xvik.def

      phivik(1)=-3637
      phivik(2)=-4505



      OPEN(99,file='xvik.def',form='formatted')

      READ(99,*) 
      READ(99,*,iostat=ierr) phivik
      IF(ierr.NE.0) GOTO 105

      READ(99,*,END=105)
      READ(99,'(a)',END=105) pathchmp
      READ(99,*,END=105)
      READ(99,'(a)',END=105) pathsor
      READ(99,*,END=105)
c     READ(99,'(l1)',END=105) day_ls
      READ(99,'(l1)',END=105)
      READ(99,'(l1)',END=105) lcal
      READ(99,'(l1)',END=105)
      READ(99,'(l1)',END=105) lvent
      READ(99,'(l1)',END=105)
      READ(99,'(l1)',END=105) latcal
 
 105  CONTINUE
      CLOSE(99)
      write (*,*)'>>>>>>>>>>>>>>>>', phivik,g
      DO iv=1,2
         phivik(iv)=phivik(iv)*3.73
      END DO

      write(*,*) ' pathchmp:',trim(pathchmp)
      write(*,*) ' pathsor:',trim(pathsor)
      
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c   ouverture des fichiers xgraph:
c   ------------------------------

      ifile(1)=12
      ifile(2)=13
      kyear=-1
c      OPEN(77,file='xlongday',form='formatted')

      unitlec=11
c
      PRINT*,'entrer le nom du fichier NC'
      READ(5,'(a)') nomfich

      PRINT *,'nomfich : ',nomfich


c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c   grande boucle sur les fichiers histoire:
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      firstcal=.true.
      DO WHILE(len_trim(nomfich).GT.0.AND.len_trim(nomfich).LT.50)
      PRINT *,'>>>  nomfich : ',trim(nomfich)

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      file=pathchmp(1:len_trim(pathchmp))//'/'//
     s     nomfich(1:len_trim(nomfich))
      PRINT*,'file.nc: ', file(1:len_trim(file))//'.nc'
      PRINT*,'timestep ',dtvr

      IF(nc) THEN
      ierr= NF_OPEN(file(1:len_trim(file))//'.nc',NF_NOWRITE,nid)        
      ELSE
         PRINT*,'Ouverture binaire ',file
         OPEN(unitlec,file=file,status='old',form='unformatted',
     .   iostat=ierr)
      ENDIF

c----------------------------------------------------------------------
c   initialisation de la physique:
c   ------------------------------

      CALL readhead_NC(file(1:len_trim(file))//'.nc',day0,phis,constR)

      WRITE (*,*) 'day0 = ' , day0

      CALL conf_gcm( 99, .TRUE. )
      CALL iniconst
      CALL inigeom
!      CALL inifilr


c   Lecture temps :

      ierr= NF_INQ_DIMID (nid,"Time",dimid)
        IF (ierr.NE.NF_NOERR) THEN
          PRINT*, 'xvik: Le champ <Time> est absent'
          CALL abort
        ENDIF

      ierr= NF_INQ_DIMLEN (nid,dimid,nbpas)

      ierr = NF_INQ_VARID (nid, "Time", nvarid)
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, temps)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, temps)
#endif
        IF (ierr.NE.NF_NOERR) THEN
          PRINT*, 'xvik: Lecture echouee pour <Time>'
          CALL abort
        ENDIF

        PRINT*,'temps(1:10)',(temps(itau),itau=1,10)
              
c-----------------------------------------------------------------------
c   coordonnees des point Viking:
c   -----------------------------

      latvik(1)=22.27*pi/180.
      lonvik(1)=-47.95*pi/180.
      latvik(2)=47.67*pi/180.
      lonvik(2)=(360.-225.71)*pi/180.

c   ponderations pour les 4 points autour de Viking
      DO iv=1,2
        ! locate index of GCM grid points near VL
         do i=1,iim
           ! we know longitudes are ordered -180...180
           if ((lonvik(iv).ge.rlonu(i)).and.
     &         (lonvik(iv).le.rlonu(i+1))) then
             ivik(iv)=i
             exit
           endif
         enddo
         do j=1,jjm-1
           !we know tha latitudes are ordered 90...-90
           if ((latvik(iv).le.rlatv(j)).and.
     &         (latvik(iv).ge.rlatv(j+1))) then
             jvik(iv)=j
             exit
           endif
         enddo
         zalpha=(lonvik(iv)-rlonu(ivik(iv)))/
     s          (rlonu(ivik(iv)+1)-rlonu(ivik(iv)))
         zbeta=(latvik(iv)-rlatv(jvik(iv)))/
     s          (rlatv(jvik(iv)+1)-rlatv(jvik(iv)))
         zw(0,0,iv)=(1.-zalpha)*(1.-zbeta)
         zw(1,0,iv)=zalpha*(1.-zbeta)
         zw(0,1,iv)=(1.-zalpha)*zbeta
         zw(1,1,iv)=zalpha*zbeta
      ENDDO

c   altitude reelle et modele aux points Viking
      DO iv=1,2
         phisim(iv)=0.
         DO jj=0,1
            j=jvik(iv)+jj
            DO ii=0,1
               i=ivik(iv)+ii
               phisim(iv)=phisim(iv)+zw(ii,jj,iv)*phis(i,j)
            ENDDO
         ENDDO
      ENDDO
      PRINT*,'relief aux points Viking pour les sorties:',phivik

c----------------------------------------------------------------------
c   lectures des etats:
c   -------------------

       airtot1=1./(SSUM(ip1jmp1,aire,1)-SSUM(jjp1,aire,iip1))

c======================================================================
c   debut de la boucle sur les etats dans un fichier histoire:
c======================================================================
       count_ps=(/iip1,jjp1,1/)
       count_co2ice=(/iip1,jjp1,1/)
       count_temp=(/iip1,jjp1,llm,1/)
       
       DO itau=1,nbpas

       start_ps=(/1,1,itau/)
       start_co2ice=(/1,1,itau/)
       start_temp=(/1,1,1,itau/)
c   lecture drs des champs:
c   -----------------------
c         varname='u'
c         ierr=drsread (unitlec,varname,unat,itau)
c         PRINT*,'unat',unat(iip1/2,jjp1/2,llm/2)
c         varname='v'
c         ierr=drsread (unitlec,varname,vnat,itau)
c         PRINT*,'vnat',vnat(iip1/2,jjp1/2,llm/2)

ccccccccc  LECTURE Ps ccccccccccccccccccccccccccc
          ierr = NF_INQ_VARID (nid, "ps", nvarid)
#ifdef NC_DOUBLE
          ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start_ps,count_ps, ps)
#else
          ierr = NF_GET_VARA_REAL(nid, nvarid,start_ps,count_ps, ps)
#endif
          IF (ierr.NE.NF_NOERR) THEN
            PRINT*, 'xvik: Lecture echouee pour <ps>'
            CALL abort
          ENDIF
          
          PRINT*,'ps',ps(iip1/2,jjp1/2)

ccccccccc  LECTURE Temperature ccccccccccccccccccccccccccc
          ierr = NF_INQ_VARID (nid, "temp", nvarid)
#ifdef NC_DOUBLE
          ierr = NF_GET_VARA_DOUBLE(nid,nvarid,start_temp,count_temp, t)
#else
          ierr = NF_GET_VARA_REAL(nid,nvarid,start_temp,count_temp, t)
#endif
          IF (ierr.NE.NF_NOERR) THEN
            PRINT*, 'xvik: Lecture echouee pour <temp>'
            ! Ehouarn: proceed anyways
            ! CALL abort
            write(*,*)'--> Setting temperature to zero !!!'
            t(1:iip1,1:jjp1,1:llm)=0.0
            write(*,*)'--> looking for temp7 (temp in 7th layer)'
            ierr=NF_INQ_VARID(nid,"temp7", nvarid)
            if (ierr.eq.NF_NOERR) then
            write(*,*) "    OK, found temp7 variable"
#ifdef NC_DOUBLE
            ierr=NF_GET_VARA_DOUBLE(nid,nvarid,start_ps,count_ps,t7)
#else
            ierr=NF_GET_VARA_REAL(nid,nvarid,start_ps,count_ps,t7)
#endif
              if (ierr.ne.NF_NOERR) then
                write(*,*)'xvik: failed loading temp7 !'
                stop
              endif
            else ! no 'temp7' variable
              write(*,*)'  No temp7 variable either !'
              write(*,*)'  Will have to to without ...'
              t7(1:iip1,1:jjp1)=0.0
            endif
          ELSE ! t() was successfully loaded, copy 7th layer to t7()
            t7(1:iip1,1:jjp1)=t(1:iip1,1:jjp1,7)
          ENDIF

c          PRINT*,'t',t(iip1/2,jjp1/2,llm/2)

ccccccccc  LECTURE co2ice ccccccccccccccccccccccccccc
          ierr = NF_INQ_VARID (nid, "co2ice", nvarid)
#ifdef NC_DOUBLE
          ierr = NF_GET_VARA_DOUBLE(nid,nvarid,start_co2ice,
     &    count_co2ice,  co2ice)
#else
          ierr = NF_GET_VARA_REAL(nid, nvarid,start_co2ice,
     &    count_co2ice, co2ice)
#endif
          IF (ierr.NE.NF_NOERR) THEN
            PRINT*, 'xvik: Lecture echouee pour <co2ice>'
            CALL abort
          ENDIF


c Gestion du temps
c ----------------
          day=temps(itau)
          PRINT*,'day ',day
          CALL solarlong(day+day0,sollong)
          sol=day+day0+461. ! aded offset to sync with VL mission "sol 1"
          iyear=sol/unanj
          WRITE (*,*) 'iyear',iyear
          sol=sol-iyear*unanj
c
c Ouverture / fermeture des fichiers
c ----------------------------------
          IF (iyear.NE.kyear) THEN
             WRITE(chr2(1:1),'(i1)') iyear+1
             WRITE (*,*) 'iyear bis',iyear
             WRITE (*,*) 'chr2'
             WRITE (*,*)  chr2
             IF(iyear.GE.9) WRITE(chr2,'(i2)') iyear+1
             kyear=iyear
             DO ii=1,2
                CLOSE(10+ifile(ii))
                CLOSE(2+ifile(ii))
                CLOSE(4+ifile(ii))
                CLOSE(6+ifile(ii))
                CLOSE(8+ifile(ii))
                CLOSE(16+ifile(ii))
                CLOSE(12+ifile(ii))
                CLOSE(14+ifile(ii))
                CLOSE(97)
                CLOSE(98)
             ENDDO
             CLOSE(5+ifile(1))
             OPEN(ifile(1)+10,file='xpsol1'//chr2,form='formatted')
             OPEN(ifile(2)+10,file='xpsol2'//chr2,form='formatted')
c            OPEN(ifile(1)+8,file='xbpsol1'//chr2,form='formatted')
c            OPEN(ifile(2)+8,file='xbpsol2'//chr2,form='formatted')
c            OPEN(ifile(1)+2,file='xlps1'//chr2,form='formatted')
c            OPEN(ifile(2)+2,file='xlps2'//chr2,form='formatted')
             IF(lcal) THEN
c               OPEN(ifile(2)+4,file='xpressud'//chr2,form='formatted')
c               OPEN(ifile(1)+4,file='xpresnord'//chr2,form='formatted')
c               OPEN(ifile(1)+6,file='xpm2'//chr2,form='formatted')
             ENDIF
                         IF(latcal) THEN
c               OPEN(ifile(2)+14,file='xlats'//chr2,form='formatted')
c               OPEN(ifile(1)+14,file='xlatn'//chr2,form='formatted')
                         ENDIF
             IF(lvent) THEN
c               OPEN(ifile(1)+16,file='xu1'//chr2,form='formatted')
c               OPEN(ifile(2)+16,file='xu2'//chr2,form='formatted')
c               OPEN(ifile(1)+12,file='xv1'//chr2,form='formatted')
c               OPEN(ifile(2)+12,file='xv2'//chr2,form='formatted')
             ENDIF
             OPEN(97,file='xprestot'//chr2,form='formatted')
c            OPEN(98,file='xlat37_'//chr2,form='formatted')
c           WRITE(98,'(f5.1,16f7.1)') 0.,(rlonv(i)*180./pi,i=1,iim,4)
          ENDIF
 

          sollong=sollong*180./pi
          IF(day_ls) THEN
             dayw=sol
             write(*,*) 'dayw', dayw
          ELSE
             dayw=sollong
          ENDIF

c Calcul de la moyenne planetaire
c -------------------------------
          pstot=0.
          captotS=0.
          captotN=0.
          DO j=1,jjp1
             DO i=1,iim
                pstot=pstot+aire(i,j)*ps(i,j)
             ENDDO
          ENDDO
 
              DO j=1,jjp1/2
                 DO i=1,iim
                    captotN = captotN  +aire(i,j)*co2ice(i,j)
                 ENDDO
              ENDDO
              DO j=jjp1/2+1, jjp1
                 DO i=1,iim
                    captotS = captotS  +aire(i,j)*co2ice(i,j)
                 ENDDO
              ENDDO
              WRITE(97,'(4e16.6)') dayw,pstot*airtot1
     &       , captotN*g*airtot1, captotS*g*airtot1          

          IF(.NOT.firstcal) THEN
c         WRITE(98,'(f5.1,16f7.3)')
c     s    dayw,(ps(i,37),i=1,iim,4)

c boucle sur les sites vikings:
c ----------------------------

          DO iv=1,2

c interpolation de la temperature dans la 7eme couche, de la pression
c de surface et des vents aux points viking.

             zp1=0.
             zp2=0.
             zp2_sm=0.
             zt=0.
!             zu=0.
!             zv=0.
             DO jj=0,1
                j=jvik(iv)+jj
                DO ii=0,1
                   i=ivik(iv)+ii
!                   zt=zt+zw(ii,jj,iv)*t(i,j,7)
                   zt=zt+zw(ii,jj,iv)*t7(i,j)
!                   zp1=zp1+zw(ii,jj,iv)*ps(i,j)
                   zp1=zp1+zw(ii,jj,iv)*log(ps(i,j)) ! interpolate in log(P)
                    WRITE (*,*) 'ps autour iv',ps(i,j),iv
!                   zu=zu+zw(ii,jj,iv)*unat(i,j,1)/cu(i,j)
!                   zv=zv+zw(ii,jj,iv)*vnat(i,j,1)/cv(i,j)
                ENDDO
             ENDDO
             zp1=exp(zp1) ! because of the bilinear interpolation in log(P)
 
c               pression au sol extrapolee a partir de la temp. 7eme couche
           WRITE (*,*) 'constR ',constR 
           WRITE (*,*) 'zt ',zt
             gh=constR*zt
             if (gh.eq.0) then ! if we don't have temperature values
               ! assume a scale height of 10km
               zp2=zp1*exp(-(phivik(iv)-phisim(iv))/(3.73*1.e4))
             else
               zp2=zp1*exp(-(phivik(iv)-phisim(iv))/gh)
             endif
           WRITE (*,*) 'iv,pstot,zp2, zp1, phivik(iv),phisim(iv),gh'
           WRITE (*,*) iv,pstot*airtot1,zp2,zp1,phivik(iv),phisim(iv),gh
!           WRITE(ifile(iv)+10,'(2e15.5)') dayw,zp1
           WRITE(ifile(iv)+10,'(3e15.5)') dayw,zp2,zp1
           
c   sorties eventuelles de vent
!             IF(lvent) THEN
!                WRITE(ifile(iv)+16,'(2e15.5)')
!     s          dayw,zu
!                WRITE(ifile(iv)+12,'(2e15.5)')
!     s          dayw,zv
!             ENDIF
          ENDDO
c         IF (lcal) THEN
c            WRITE(ifile(1)+4,'(2e15.6)') dayw,airtot1*g*.01*
c    s       (SSUM(ip1jmp1/2,ziceco2,1)-SSUM(jjp1/2,ziceco2,iip1))
c            WRITE(ifile(2)+4,'(2e15.6)') dayw,airtot1*g*.01*
c    s       (SSUM(iip1*jjm/2,ziceco2(1,jjm/2+2),1)-
c    s       SSUM(jjm/2,ziceco2(1,jjm/2+2),iip1))
c         ENDIF
c            IF(latcal) THEN
c               CALL icelat(iim,jjm,ziceco2,rlatv,zicelat)
c               WRITE(ifile(1)+14,'(2e15.6)') dayw,zicelat(1)*180./pi
c               WRITE(ifile(2)+14,'(2e15.6)') dayw,zicelat(2)*180./pi
c            ENDIF
         ENDIF
         firstcal=.false.

c======================================================================
c   Fin de la boucle sur les etats du fichier histoire:
c======================================================================
      ENDDO

      ierr= NF_CLOSE(nid)

      PRINT*,'Fin du fichier',nomfich
      print*,'Entrer un nouveau fichier ou return pour finir'
      READ(5,'(a)',err=9999) nomfich

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c   Fin de la boucle sur les fichiers histoire:
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      ENDDO

      PRINT*,'relief du point V1',.001*phis(ivik(1),jvik(1))/g
      PRINT*,'relief du point V2',.001*phis(ivik(2),jvik(2))/g
      DO iv=1,2
         PRINT*,'Viking',iv,'   i=',ivik(iv),'j  =',jvik(iv)
         WRITE(6,7777)
     s   (rlonv(i)*180./pi,i=ivik(iv)-1,ivik(iv)+2)
         print*
         DO j=jvik(iv)-1,jvik(iv)+2
            WRITE(6,'(f8.1,10x,5f7.1)')
     s   rlatu(j)*180./pi,(phis(i,j)/(g*1000.),i=ivik(iv)-1,ivik(iv)+2)
         ENDDO
         print*
         print*,'zw'
         write(6,'(2(2f10.4/))') ((zw(ii,jj,iv),ii=0,1),jj=0,1)
         print*,'altitude interpolee (km) ',phisim(iv)/1000./g
      ENDDO
      PRINT*,'R=',r
 9999  PRINT*,'Fin '

7777  FORMAT ('latitude/longitude',4f7.1)
      END
