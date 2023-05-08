      SUBROUTINE calchim(nlon,ny,qy_c,nomqy_c,declin_rad,ls_rad,dtchim,
     .                   ctemp,cplay,cplev,czlay,czlev,
     .                   dqyc)
      
c-------------------------------------------------
c  
c     Introduction d une routine chimique
c
c     Auteur: S. Lebonnois,  01/2000 | 09/2003 
c            adaptation pour Titan 3D: 02/2009
c            adaptation pour // : 04/2013
c            extension chimie jusqu a 1300km : 10/2013
c 
c-------------------------------------------------
c
      use dimphy
      use common_mod, only:utilaer,maer,prodaer,csn,csh,psurfhaze,
     .                     NLEV,NLRT,NC,ND,NR
      USE moyzon_mod, only: tmoy,playmoy,zlaymoy,zlevmoy,klat
      use mod_grid_phy_lmdz, only: nbp_lat
      implicit none
#include "clesphys.h"
#include "YOMCST.h"

c    Arguments
c    ---------

      INTEGER      nlon                   ! nb of horiz points
      INTEGER      ny                     ! nb de composes (nqmax-nmicro) 
      REAL         qy_c(nlon,klev,NC)     ! Especes chimiques apres adv.+diss.
      character*10 nomqy_c(NC+1)          ! Noms des especes chimiques
      REAL         declin_rad,ls_rad      ! declinaison et long solaire en radians
      REAL         dtchim                 ! pas de temps chimie
      REAL         ctemp(nlon,klev)      ! Temperature
      REAL         cplay(nlon,klev)      ! pression (Pa)
      REAL         cplev(nlon,klev+1)    ! pression intercouches (Pa)
      REAL         czlay(nlon,klev)      ! altitude (m)
      REAL         czlev(nlon,klev+1)    ! altitude intercouches (m)
      
      REAL         dqyc(nlon,klev,NC)    ! Tendances especes chimiques 
      
c    Local variables :
c    -----------------

      integer i,j,l,ic,jm1

c variables envoyees dans la chimie: double precision

      REAL  temp_c(NLEV)
      REAL  press_c(NLEV)   ! T,p(mbar) a 1 lat donnee
      REAL  cqy(NLEV,NC)    ! frac mol qui seront modifiees
      REAL  cqy0(NLEV,NC)    ! frac mol avant modif
      REAL  surfhaze(NLEV)
      REAL  cprodaer(NLEV,4),cmaer(NLEV,4)
      REAL  ccsn(NLEV,4),ccsh(NLEV,4)
! rmil: milieu de couche, grille pour K,D,p,ct,T,y
! rinter: intercouche (grille RA dans la chimie)
      REAL  rmil(NLEV),rinter(NLEV),nb(NLEV)
      REAL,save :: kedd(NLEV)

      REAL  a,b
      character str1*1,str2*2
      REAL  latit
      character*20 formt1,formt2,formt3
      
c     variables locales initialisees au premier appel

      LOGICAL firstcal
      DATA    firstcal/.true./
      SAVE    firstcal
      
      integer dec,declinint,ialt
      REAL  declin_c                       ! declinaison en degres
      real  factalt,factdec,krpddec,krpddecp1,krpddecm1
      real  duree
      REAL,save :: mass(NC)
      REAL,save,allocatable :: md(:,:)
      REAL,save :: botCH4
      DATA  botCH4/0.05/
      real,save ::  r1d(131),ct1d(131),p1d(131),t1d(131) ! lecture tcp.ver
      REAL,save,allocatable :: krpd(:,:,:,:),krate(:,:)
      integer,save :: reactif(5,NR),nom_prod(NC),nom_perte(NC)
      integer,save :: prod(200,NC),perte(2,200,NC)
      character  dummy*30,fich*7,ficdec*15,curdec*15,name*10
      real  ficalt,oldq,newq,zalt
      logical okfic

c-----------------------------------------------------------------------
c***********************************************************************
c
c    Initialisations :
c    ----------------

      duree = dtchim  ! passage en real*4 pour appel a routines C

      IF (firstcal) THEN
            
        print*,'CHIMIE, premier appel'
        
c ************************************
c Au premier appel, initialisation de toutes les variables 
c pour les routines de la chimie.
c ************************************

        allocate(krpd(15,ND+1,NLRT,nbp_lat),krate(NLEV,NR),md(NLEV,NC))

c Verification du nombre de composes: coherence common_mod et nqmax-nmicro
c ----------------------------------

        if (ny.ne.NC) then
	   print*,'PROBLEME de coherence dans le nombre de composes:'
     .           ,ny,NC
           STOP
	endif

c calcul de temp_c, densites et press_c en moyenne planetaire :
c --------------------------------------------------------------

        print*,'pression, densites et temp (init chimie):'
        print*,'level, press_c, nb, temp_c'
        DO l=1,klev
         rinter(l)  = (zlevmoy(l)+RA)/1000.
         rmil(l)    = (zlaymoy(l)+RA)/1000.
c     temp_c (K):
         temp_c(l)  = tmoy(l)
c     press_c (mbar):
         press_c(l) = playmoy(l)/100.
c     nb (cm-3):
         nb(l) = 1.e-4*press_c(l) / (RKBOL*temp_c(l))
         print*,l,rmil(l)-RA/1000.,press_c(l),nb(l),temp_c(l)
        ENDDO
        rinter(klev+1)=(zlevmoy(klev+1)+RA)/1000.

c au-dessus du GCM, dr regulier et rinter(NLEV)=1290+2575 km.
       do l=klev+2,NLEV
         rinter(l) = rinter(klev+1)
     &          + (l-klev-1)*(3865.-rinter(klev+1))/(NLEV-klev-1)
         rmil(l-1) = (rinter(l-1)+rinter(l))/2.
       enddo
       rmil(NLEV) = rinter(NLEV)+(rinter(NLEV)-rinter(NLEV-1))/2.

c lecture de tcp.ver, une seule fois
c remplissage r1d,t1d,ct1d,p1d
        open(11,file='../../INPUT/tcp.ver',status='old')
        read(11,*) dummy
        do i=1,131
          read(11,*) r1d(i),t1d(i),ct1d(i),p1d(i)
c         print*,"TCP.VER ",r1d(i),t1d(i),ct1d(i),p1d(i)
        enddo 
        close(11)

c extension pour klev+1 a NLEV avec tcp.ver
c la jonction klev/klev+1 est brutale... Tant pis ?
        ialt=1
        do l=klev+1,NLEV
           zalt=rmil(l)-RA/1000.
           do i=ialt,130
            if ((zalt.ge.r1d(i)).and.(zalt.lt.r1d(i+1))) then
              ialt=i
            endif
           enddo
           factalt = (zalt-r1d(ialt))/(r1d(ialt+1)-r1d(ialt))
           press_c(l) = exp(  log(p1d(ialt))   * (1-factalt) 
     &                      + log(p1d(ialt+1)) * factalt    )
           nb(l)      = exp(  log(ct1d(ialt))   * (1-factalt) 
     &                      + log(ct1d(ialt+1)) * factalt    )
           temp_c(l)  = t1d(ialt) * (1-factalt) + t1d(ialt+1) * factalt
           print*,l,zalt,press_c(l),nb(l),temp_c(l)
        enddo
        
c caracteristiques des composes:        
c -----------------------------
        mass(:) = 0.0
        call comp(nomqy_c,nb,temp_c,mass,md)
        print*,'           Mass'
        do ic=1,NC
          print*,nomqy_c(ic),mass(ic)
c         if(nomqy_c(ic).eq.'CH4') print*,"MD=",md(:,ic)
        enddo
	

c  Stockage des composes utilises dans la prod d aerosols
c     (aerprod=1) et dans H-> H2 (htoh2=1): utilaer
c     ! decalage de 1 car utilise dans le c !

        do ic=1,NC

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!!remise de CH4 a 1.5%!!!!!!!!!!!!!!!!!!!!!!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c         if (nomqy_c(ic).eq."CH4") then
c           do l=1,llm
c            do j=1,ip1jmp1
c              if (qy_c(j,l,ic).le.0.015) qy_c(j,l,ic) = 0.015 
c            enddo
c           enddo
c         endif
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
	  if (nomqy_c(ic).eq."C4H2") then
            utilaer(10) = ic-1
	  endif
	  if (nomqy_c(ic).eq."HCN") then
            utilaer(6)  = ic-1
	  endif
	  if (nomqy_c(ic).eq."HC3N") then
            utilaer(7)  = ic-1
	  endif
	  if (nomqy_c(ic).eq."NCCN") then
            utilaer(14) = ic-1
	  endif
	  if (nomqy_c(ic).eq."CH3CN") then
            utilaer(15) = ic-1
            utilaer(16) = ic-1 ! si pas C2H3CN, CH3CN utilise, mais reac. nulle
	  endif
	  if (nomqy_c(ic).eq."H") then
            utilaer(1)  = ic-1
	  endif
	  if (nomqy_c(ic).eq."H2") then
            utilaer(2)  = ic-1
	  endif
	  if (nomqy_c(ic).eq."C2H2") then
            utilaer(3)  = ic-1
	  endif
	  if (nomqy_c(ic).eq."AC6H6") then
            utilaer(13) = ic-1
	  endif
	  if (nomqy_c(ic).eq."C2H3CN") then
            utilaer(16) = ic-1
	  endif
	  if (nomqy_c(ic).eq."C2") then
            utilaer(4)  = ic-1
	  endif
	  if (nomqy_c(ic).eq."C2H") then
            utilaer(5)  = ic-1
	  endif
	  if (nomqy_c(ic).eq."C3N") then
            utilaer(8)  = ic-1
	  endif
	  if (nomqy_c(ic).eq."H2CN") then
            utilaer(9)  = ic-1
	  endif
	  if (nomqy_c(ic).eq."C4H3") then
            utilaer(11) = ic-1
	  endif
	  if (nomqy_c(ic).eq."AC6H5") then
            utilaer(12) = ic-1
	  endif

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c         if ((nomqy_c(ic).eq."HC3N").or.
c    $        (nomqy_c(ic).eq."C3N")) then
c     DO j=1,ip1jmp1 
c       do l=1,34  ! p>~ 1 mbar
c         qy_c(j,l,ic) = 1.e-30
c       enddo
c     ENDDO
c         endif
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        enddo
	        
c taux de photodissociations:
c --------------------------
        call disso(krpd,nbp_lat) 

c reactions chimiques:
c -------------------
        call chimie(nomqy_c,nb,temp_c,krate,reactif,
     .               nom_perte,nom_prod,perte,prod)
c        print*,'nom_prod, nom_perte:'
c        do ic=1,NC
c          print*,nom_prod(ic),nom_perte(ic)
c        enddo
c        print*,'premieres prod, perte(1:reaction,2:compagnon):'
c        do ic=1,NC
c          print*,prod(1,ic),perte(1,1,ic),perte(2,1,ic)
c        enddo

c       l = klev-3
c       print*,'krate a p=',press_c(l),' reactifs et produits:'
c       do ic=1,ND+1
c         print*,ic,krpd(7,ic,l,4)*nb(l),"  ",
c    .     nomqy_c(reactif(1,ic)+1),
c    .     nomqy_c(reactif(2,ic)+1),nomqy_c(reactif(3,ic)+1),
c    .     nomqy_c(reactif(4,ic)+1),nomqy_c(reactif(5,ic)+1)
c       enddo
c       do ic=ND+2,NR
c         print*,ic,krate(l,ic),"  ",
c    .     nomqy_c(reactif(1,ic)+1),
c    .     nomqy_c(reactif(2,ic)+1),nomqy_c(reactif(3,ic)+1),
c    .     nomqy_c(reactif(4,ic)+1),nomqy_c(reactif(5,ic)+1)
c       enddo


c   init kedd
c   ---------
c   kedd stays constant with time and space 
c   => linked to pressure rather than altitude...
 
      kedd(:) = 5e2 ! valeur mise par defaut  
               ! UNITE ? doit etre ok pour gptitan
      do l=1,NLEV
       zalt=rmil(l)-RA/1000.  ! z en km
       if     ((zalt.ge.250.).and.(zalt.le.400.)) then 
         kedd(l) = 10.**(3.+(zalt-250.)/50.)
! 1E3 at 250 km
! 1E6 at 400 km
       elseif ((zalt.gt.400.).and.(zalt.le.600.)) then 
         kedd(l) = 10.**(6.+1.3*(zalt-400.)/200.)
! 2E7 at 600 km
       elseif ((zalt.gt.600.).and.(zalt.le.900.)) then 
         kedd(l) = 10.**(7.3+0.7*(zalt-600.)/300.)
! 1E8 above 900 km
       elseif ( zalt.gt.900.                    ) then 
        kedd(l) = 1.e8
       endif
      enddo

      ENDIF  ! premier appel

c***********************************************************************
c-----------------------------------------------------------------------

c   calcul declin_c (en degres)
c   ---------------------------

      declin_c = declin_rad*180./RPI
c      print*,'declinaison en degre=',declin_c
       
c***********************************************************************
c***********************************************************************
c
c                BOUCLE SUR LES LATITUDES
c
      DO j=1,nlon 
      
      if (j.eq.1) then 
         jm1=1
      else
         jm1=j-1
      endif

      if((j.eq.1).or.(klat(j).ne.klat(jm1))) then

c***********************************************************************
c***********************************************************************

c-----------------------------------------------------------------------
c
c   Distance radiale (intercouches, en km)
c   ----------------------------------------

       do l=1,klev
         rinter(l) = (RA+czlev(j,l))/1000.
         rmil(l)   = (RA+czlay(j,l))/1000.
c        print*,rinter(l)
       enddo
       rinter(klev+1)=(RA+czlev(j,klev+1))/1000.

c au-dessus du GCM, dr regulier et rinter(NLEV)=1290+2575 km.
       do l=klev+2,NLEV
         rinter(l) = rinter(klev+1)
     &          + (l-klev-1)*(3865.-rinter(klev+1))/(NLEV-klev-1)
         rmil(l-1) = (rinter(l-1)+rinter(l))/2.
       enddo
       rmil(NLEV) = rinter(NLEV)+(rinter(NLEV)-rinter(NLEV-1))/2.

c-----------------------------------------------------------------------
c
c   Temperature, pression (mbar), densite (cm-3)
c   -------------------------------------------

       DO l=1,klev
c     temp_c (K):
               temp_c(l)  = ctemp(j,l)
c     press_c (mbar):
               press_c(l) = cplay(j,l)/100.
c     nb (cm-3):
               nb(l) = 1.e-4*press_c(l) / (RKBOL*temp_c(l))
       ENDDO
c extension pour klev+1 a NLEV avec tcp.ver
       ialt=1
       do l=klev+1,NLEV
           zalt=rmil(l)-RA/1000.
           do i=ialt,130
            if ((zalt.ge.r1d(i)).and.(zalt.lt.r1d(i+1))) then
              ialt=i
            endif
           enddo
           factalt = (zalt-r1d(ialt))/(r1d(ialt+1)-r1d(ialt))
           press_c(l) = exp(  log(p1d(ialt))   * (1-factalt) 
     &                      + log(p1d(ialt+1)) * factalt    )
           nb(l)      = exp(  log(ct1d(ialt))   * (1-factalt) 
     &                      + log(ct1d(ialt+1)) * factalt    )
           temp_c(l)  = t1d(ialt) * (1-factalt) + t1d(ialt+1) * factalt
       enddo
               
c-----------------------------------------------------------------------
c
c   passage krpd => krate 
c   ---------------------
c initialisation krate pour dissociations 

      if ((declin_c*10+267).lt.14.) then
          declinint = 0
          dec       = 0
      else 
       if ((declin_c*10+267).gt.520.) then
          declinint = 14
          dec       = 534
       else 
          declinint = 1
          dec       = 27
          do while( (declin_c*10+267).ge.real(dec+20) )
            dec       = dec+40
            declinint = declinint+1
          enddo
       endif
      endif
      if ((declin_c.ge.-24.).and.(declin_c.le.24.)) then
          factdec = ( declin_c - (dec-267)/10. ) / 4.
      else
          factdec = ( declin_c - (dec-267)/10. ) / 2.7
      endif

      do l=1,NLEV

c INTERPOL EN ALT POUR k (krpd tous les deux km)
        ialt    = int((rmil(l)-RA/1000.)/2.)+1
        factalt = (rmil(l)-RA/1000.)/2.-(ialt-1)

        do i=1,ND+1
          krpddecm1 = krpd(declinint  ,i,ialt  ,klat(j)) * (1-factalt)
     &              + krpd(declinint  ,i,ialt+1,klat(j)) * factalt
          krpddec   = krpd(declinint+1,i,ialt  ,klat(j)) * (1-factalt)
     &              + krpd(declinint+1,i,ialt+1,klat(j)) * factalt
          krpddecp1 = krpd(declinint+2,i,ialt  ,klat(j)) * (1-factalt)
     &              + krpd(declinint+2,i,ialt+1,klat(j)) * factalt

                    ! ND+1 correspond a la dissociation de N2 par les GCR
          if ( factdec.lt.0. ) then 
        krate(l,i) = krpddecm1 * abs(factdec)
     &             + krpddec   * ( 1 + factdec)
          endif
          if ( factdec.gt.0. ) then
        krate(l,i) = krpddecp1 * factdec
     &             + krpddec   * ( 1 - factdec)
          endif
          if ( factdec.eq.0. ) krate(l,i) = krpddec
        enddo        
      enddo        

c-----------------------------------------------------------------------
c
c   composition 
c   ------------

       do ic=1,NC
        do l=1,klev
          cqy(l,ic)  = qy_c(j,l,ic) 
          cqy0(l,ic) = cqy(l,ic)
        enddo
       enddo

c lecture du fichier compo_klat(j) (01 à 49) pour klev+1 a NLEV

      WRITE(str2,'(i2.2)') klat(j)
      fich = "comp_"//str2//".dat"
      inquire (file=fich,exist=okfic)
      if (okfic) then
       open(11,file=fich,status='old')
c premiere ligne=declin
       read(11,'(A15)') ficdec
       write(curdec,'(E15.5)') declin_c
c si la declin est la meme: ce fichier a deja ete reecrit 
c on lit la colonne t-1 au lieu de la colonne t
c (cas d une bande de latitude partagee par 2 procs)
       do ic=1,NC
        read(11,'(A10)') name
        if (name.ne.nomqy_c(ic)) then
          print*,"probleme lecture ",fich
          print*,name,nomqy_c(ic)
          stop
        endif
        if (ficdec.eq.curdec) then
          do l=klev+1,NLEV
            read(11,*) ficalt,cqy(l,ic),newq
          enddo
        else
          do l=klev+1,NLEV
            read(11,*) ficalt,oldq,cqy(l,ic)
          enddo
        endif
       enddo
       close(11)
      else       ! le fichier n'est pas la
       do ic=1,NC
        do l=klev+1,NLEV
          cqy(l,ic)=cqy(klev,ic)
        enddo
       enddo
      endif
      cqy0 = cqy

c-----------------------------------------------------------------------
c
c   total haze area (um2/cm3)
c   -------------------------

       if (htoh2.eq.1) then
        do l=1,klev
! ATTENTION, INVERSION PAR RAPPORT A pg2.F !!!
         surfhaze(l) = psurfhaze(j,klev+1-l) 
c        if (j.eq.25)
c    .    print*,'psurfhaze en um2/cm3:',surfhaze(l)
        enddo
       endif

c-----------------------------------------------------------------------
c
c   Appel de chimietitan
c   --------------------
       
       call gptitan(rinter,temp_c,nb,
     $              nomqy_c,cqy,
     $              duree,(klat(j)-1),mass,md,
     $              kedd,botCH4,krate,reactif,
     $              nom_prod,nom_perte,prod,perte,
     $              aerprod,utilaer,cmaer,cprodaer,ccsn,ccsh,
     $              htoh2,surfhaze)
       
c   Tendances composition 
c   ---------------------
     
       do ic=1,NC
        do l=1,klev
          dqyc(j,l,ic) = (cqy(l,ic) - cqy0(l,ic))/dtchim  ! en /s
        enddo
       enddo

c-----------------------------------------------------------------------
c
c   production aer 
c   --------------
       
       if (aerprod.eq.1) then

       do ic=1,4
        do l=1,klev
          prodaer(j,l,ic) = cprodaer(l,ic) 
          maer(j,l,ic)    = cmaer(l,ic) 
          csn(j,l,ic)     = ccsn(l,ic) 
          csh(j,l,ic)     = ccsh(l,ic) 
        enddo
       enddo

       endif

c-----------------------------------------------------------------------
c
c   sauvegarde compo sur NLEV 
c   -------------------------

c dans fichier compo_klat(j) (01 à 49)
       
      open(11,file=fich,status='replace')
c premiere ligne=declin
      write(11,'(E15.5)') declin_c
      do ic=1,NC
       write(11,'(A10)') nomqy_c(ic)
       do l=klev+1,NLEV
        write(11,'(E15.5,1X,E15.5,1X,E15.5)') rmil(l)-RA/1000.,
     .                                cqy0(l,ic),cqy(l,ic)
       enddo
      enddo
      close(11)

c***********************************************************************
c***********************************************************************

c              FIN: BOUCLE SUR LES LATITUDES 

      else      ! same latitude, we don't do calculations again
        dqyc(j,:,:) = dqyc(jm1,:,:)
        if (aerprod.eq.1) then
          prodaer(j,:,:) = prodaer(jm1,:,:) 
          maer(j,:,:)    = maer(jm1,:,:)
          csn(j,:,:)     = csn(jm1,:,:)
          csh(j,:,:)     = csh(jm1,:,:)
        endif
      endif

      ENDDO
     
c***********************************************************************
c***********************************************************************


      firstcal = .false.
      RETURN
      END
