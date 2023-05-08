!
! $Id: read_pstoke0.F 1403 2010-07-01 09:02:53Z fairhead $
!
c
c
	subroutine read_pstoke0(irec,
     .   zrec,zkon,zkev,airefi,phisfi,
     .   t,mfu,mfd,en_u,de_u,en_d,de_d,coefh,
     .   fm_therm,en_therm,
     .   frac_impa,frac_nucl,pyu1,pyv1,ftsol,psrf)

C******************************************************************************
C  Frederic HOURDIN, Abderrahmane IDELKADI
C Lecture des parametres physique stockes online necessaires pour
C recalculer offline le transport des traceurs sur la meme grille que online
C A FAIRE : une seule routine au lieu de 2 (lectflux, redecoupe)!
C******************************************************************************

	use netcdf
       USE dimphy
       USE control_mod

       IMPLICIT NONE

#include "netcdf.inc"
#include "dimensions.h"
#include "paramet.h"
#include "comconst.h"
#include "comgeom.h"
#include "temps.h"
#include "ener.h"
#include "logic.h"
#include "description.h"
#include "serre.h"
#include "indicesol.h"
cccc#include "dimphy.h"
	  
	  integer kon,kev,zkon,zkev
	  parameter(kon=iim*(jjm-1)+2,kev=llm)
	  REAL phisfi(kon)
          REAL phisfi2(iim,jjm+1),airefi2(iim,jjm+1)

          REAL mfu(kon,kev), mfd(kon,kev)
          REAL en_u(kon,kev), de_u(kon,kev)
          REAL en_d(kon,kev), de_d(kon,kev)
          REAL coefh(kon,kev)

c abd 25 11 02
c Thermiques
	 REAL fm_therm(kon,kev),en_therm(kon,kev)
		REAL t(kon,kev)

          REAL mfu2(iim,jjm+1,kev), mfd2(iim,jjm+1,kev)
          REAL en_u2(iim,jjm+1,kev), de_u2(iim,jjm+1,kev)
          REAL en_d2(iim,jjm+1,kev), de_d2(iim,jjm+1,kev)
          REAL coefh2(iim,jjm+1,kev)
		REAL t2(iim,jjm+1,kev)
c Thermiques
	 REAL fm_therm2(iim,jjm+1,kev)
         REAL en_therm2(iim,jjm+1,kev)       

          REAL pl(kev)
          integer irec
          integer xid,yid,zid,tid
          integer zrec,zim,zjm
          integer ncrec,nckon,nckev,ncim,ncjm

          real airefi(kon)
          character*20 namedim

c  !! attention !!
c attention il y a aussi le pb de def kon
c dim de phis??

          REAL frac_impa(kon,kev), frac_nucl(kon,kev)
          REAL frac_impa2(iim,jjm+1,kev), 
     .     frac_nucl2(iim,jjm+1,kev)
          REAL pyu1(kon), pyv1(kon)
          REAL pyu12(iim,jjm+1), pyv12(iim,jjm+1)
          REAL ftsol(kon,nbsrf)
          REAL psrf(kon,nbsrf) 
	  REAL ftsol1(kon),ftsol2(kon),ftsol3(kon),ftsol4(kon)
          REAL psrf1(kon),psrf2(kon),psrf3(kon),psrf4(kon)
          REAL ftsol12(iim,jjm+1),ftsol22(iim,jjm+1),
     .     ftsol32(iim,jjm+1),
     .     ftsol42(iim,jjm+1)
          REAL psrf12(iim,jjm+1),psrf22(iim,jjm+1),psrf32(iim,jjm+1),
     .     psrf42(iim,jjm+1)
	
	  integer ncidp
          save ncidp
          integer varidmfu, varidmfd, varidps, varidenu, variddeu	
		integer varidt
          integer varidend,varidded,varidch,varidfi,varidfn
c therm
	  integer varidfmth,varidenth
          integer varidyu1,varidyv1,varidpl,varidai,varididvt
          integer varidfts1,varidfts2,varidfts3,varidfts4
          integer varidpsr1,varidpsr2,varidpsr3,varidpsr4
          save varidmfu, varidmfd, varidps, varidenu, variddeu
		save varidt
          save varidend,varidded,varidch,varidfi,varidfn
c therm
	   save varidfmth,varidenth
          save varidyu1,varidyv1,varidpl,varidai,varididvt
          save varidfts1,varidfts2,varidfts3,varidfts4
          save varidpsr1,varidpsr2,varidpsr3,varidpsr4

          integer l, i
          integer start(4),count(4),status
          real rcode
          logical first
          save first
          data first/.true./



c ---------------------------------------------
c   Initialisation de la lecture des fichiers
c ---------------------------------------------

      if (irec .eq. 0) then

            rcode=nf90_open('phystoke.nc',nf90_nowrite,ncidp)

            rcode = nf90_inq_varid(ncidp, 'phis', varidps)
            print*,'ncidp,varidps',ncidp,varidps

            rcode = nf90_inq_varid(ncidp, 'sig_s', varidpl)
            print*,'ncidp,varidpl',ncidp,varidpl

            rcode = nf90_inq_varid(ncidp, 'aire', varidai)
            print*,'ncidp,varidai',ncidp,varidai

                rcode = nf90_inq_varid(ncidp, 't', varidt)
                print*,'ncidp,varidt',ncidp,varidt

            rcode = nf90_inq_varid(ncidp, 'mfu', varidmfu)
            print*,'ncidp,varidmfu',ncidp,varidmfu

            rcode = nf90_inq_varid(ncidp, 'mfd', varidmfd)
            print*,'ncidp,varidmfd',ncidp,varidmfd

            rcode = nf90_inq_varid(ncidp, 'en_u', varidenu)
            print*,'ncidp,varidenu',ncidp,varidenu

            rcode = nf90_inq_varid(ncidp, 'de_u', variddeu)
            print*,'ncidp,variddeu',ncidp,variddeu

            rcode = nf90_inq_varid(ncidp, 'en_d', varidend)
            print*,'ncidp,varidend',ncidp,varidend
	
            rcode = nf90_inq_varid(ncidp, 'de_d', varidded)
            print*,'ncidp,varidded',ncidp,varidded
	
            rcode = nf90_inq_varid(ncidp, 'coefh', varidch)
            print*,'ncidp,varidch',ncidp,varidch

c Thermiques
            rcode = nf90_inq_varid(ncidp, 'fm_th', varidfmth)
            print*,'ncidp,varidfmth',ncidp,varidfmth

            rcode = nf90_inq_varid(ncidp, 'en_th', varidenth)
            print*,'ncidp,varidenth',ncidp,varidenth
	
	    rcode = nf90_inq_varid(ncidp, 'frac_impa', varidfi)
            print*,'ncidp,varidfi',ncidp,varidfi
	
	    rcode = nf90_inq_varid(ncidp, 'frac_nucl', varidfn)
            print*,'ncidp,varidfn',ncidp,varidfn
	
            rcode = nf90_inq_varid(ncidp, 'pyu1', varidyu1)
            print*,'ncidp,varidyu1',ncidp,varidyu1
	
            rcode = nf90_inq_varid(ncidp, 'pyv1', varidyv1)
            print*,'ncidp,varidyv1',ncidp,varidyv1
	
            rcode = nf90_inq_varid(ncidp, 'ftsol1', varidfts1)
            print*,'ncidp,varidfts1',ncidp,varidfts1
	
            rcode = nf90_inq_varid(ncidp, 'ftsol2', varidfts2)
            print*,'ncidp,varidfts2',ncidp,varidfts2
         
            rcode = nf90_inq_varid(ncidp, 'ftsol3', varidfts3)
            print*,'ncidp,varidfts3',ncidp,varidfts3
  
            rcode = nf90_inq_varid(ncidp, 'ftsol4', varidfts4)
            print*,'ncidp,varidfts4',ncidp,varidfts4
	
            rcode = nf90_inq_varid(ncidp, 'psrf1', varidpsr1)
            print*,'ncidp,varidpsr1',ncidp,varidpsr1
	
            rcode = nf90_inq_varid(ncidp, 'psrf2', varidpsr2)
            print*,'ncidp,varidpsr2',ncidp,varidpsr2
	
	    rcode = nf90_inq_varid(ncidp, 'psrf3', varidpsr3)
            print*,'ncidp,varidpsr3',ncidp,varidpsr3

            rcode = nf90_inq_varid(ncidp, 'psrf4', varidpsr4)
            print*,'ncidp,varidpsr4',ncidp,varidpsr4
	
c ID pour les dimensions

            status = nf_inq_dimid(ncidp,'y',yid)
            status = nf_inq_dimid(ncidp,'x',xid)
            status = nf_inq_dimid(ncidp,'sig_s',zid)
            status = nf_inq_dimid(ncidp,'time_counter',tid)

c lecture des dimensions

            status = nf_inq_dim(ncidp,yid,namedim,ncjm)
            status = nf_inq_dim(ncidp,xid,namedim,ncim)
            status = nf_inq_dim(ncidp,zid,namedim,nckev)
            status = nf_inq_dim(ncidp,tid,namedim,ncrec)
	
            zrec=ncrec
            zkev=nckev
            zim=ncim
            zjm=ncjm
	
	    zkon=zim*(zjm-2)+2
	
	    write(*,*) 'read_pstoke : zrec = ', zrec
            write(*,*) 'read_pstoke : kev = ', zkev
            write(*,*) 'read_pstoke : zim = ', zim 
            write(*,*) 'read_pstoke : zjm = ', zjm
            write(*,*) 'read_pstoke : kon = ', zkon

c niveaux de pression

            status=NF_GET_VARA_REAL(ncidp,varidpl,1,kev,pl)

c lecture de aire et phis
	
      start(1)=1
      start(2)=1
      start(3)=1
      start(4)=0

      count(1)=zim
      count(2)=zjm
      count(3)=1
      count(4)=0

c 
C**** Geopotentiel au sol ***************************************
c phis
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidps,start,count,phisfi2)
#else
        status=NF_GET_VARA_REAL(ncidp,varidps,start,count,phisfi2)
#endif
      call gr_ecrit_fi(1,kon,iim,jjm+1,phisfi2,phisfi)

C**** Aires des mails aux sol ************************************
c aire
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidai,start,count,airefi2)
#else
        status=NF_GET_VARA_REAL(ncidp,varidai,start,count,airefi2)
#endif
      call gr_ecrit_fi(1,kon,iim,jjm+1,airefi2,airefi)
      else

      print*,'ok1'

c ---------------------
c   lecture des champs
c ---------------------
	
	print*,'WARNING!!! Il n y a pas de test de coherence'
        print*,'sur le nombre de niveaux verticaux dans le fichier nc'

      start(1)=1
      start(2)=1
      start(3)=1
      start(4)=irec

      count(1)=zim
      count(2)=zjm
      count(3)=kev
      count(4)=1

C**** Temperature ********************************************
cA FAIRE : Es-ce necessaire ?

c abder t
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidt,start,count,t2)
#else
        status=NF_GET_VARA_REAL(ncidp,varidt,start,count,t2)
#endif
      call gr_ecrit_fi(kev,kon,iim,jjm+1,t2,t)

C**** Flux pour la convection (Tiedtk) ********************************************
c mfu
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidmfu,start,count,mfu2)
#else
        status=NF_GET_VARA_REAL(ncidp,varidmfu,start,count,mfu2)
#endif
      call gr_ecrit_fi(kev,kon,iim,jjm+1,mfu2,mfu)

c mfd
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidmfd,start,count,mfd2)
#else
        status=NF_GET_VARA_REAL(ncidp,varidmfd,start,count,mfd2)
#endif
      call gr_ecrit_fi(kev,kon,iim,jjm+1,mfd2,mfd)

c en_u 
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidenu,start,count,en_u2)
#else
        status=NF_GET_VARA_REAL(ncidp,varidenu,start,count,en_u2)
#endif
      call gr_ecrit_fi(kev,kon,iim,jjm+1,en_u2,en_u)

c de_u 
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,variddeu,start,count,de_u2)
#else
        status=NF_GET_VARA_REAL(ncidp,variddeu,start,count,de_u2)
#endif
      call gr_ecrit_fi(kev,kon,iim,jjm+1,de_u2,de_u)

c en_d 
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidend,start,count,en_d2)
#else
        status=NF_GET_VARA_REAL(ncidp,varidend,start,count,en_d2)
#endif
      call gr_ecrit_fi(kev,kon,iim,jjm+1,en_d2,en_d)

c de_d 
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidded,start,count,de_d2)
#else
        status=NF_GET_VARA_REAL(ncidp,varidded,start,count,de_d2)
#endif
      call gr_ecrit_fi(kev,kon,iim,jjm+1,de_d2,de_d)

C**** Coefficient de mellange turbulent *******************************************
c coefh 
	print*,'LECTURE de coefh a irec =',irec
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidch,start,count,coefh2)
#else
        status=NF_GET_VARA_REAL(ncidp,varidch,start,count,coefh2)
#endif
       call gr_ecrit_fi(kev,kon,iim,jjm+1,coefh2,coefh)
c      call dump2d(iip1,jjp1,coefh2(1,2),'COEFH2READ   ')
c      call dump2d(iim ,jjm ,coefh (2,2),'COEFH2READ   ')

C**** Flux ascendants et entrant dans le thermique **********************************
cThermiques
       print*,'LECTURE de fm_therm a irec =',irec
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidfmth,start,
     .                         count,fm_therm2)
#else
       status=NF_GET_VARA_REAL(ncidp,varidfmth,start,
     .                         count,fm_therm2)
#endif
       call gr_ecrit_fi(kev,kon,iim,jjm+1,fm_therm2,fm_therm)
       print*,'LECTURE de en_therm a irec =',irec
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidenth,start,
     .                          count,en_therm2)
#else
       status=NF_GET_VARA_REAL(ncidp,varidenth,start,
     .                          count,en_therm2)
#endif
       call gr_ecrit_fi(kev,kon,iim,jjm+1,en_therm2,en_therm)

C**** Coefficients de lessivage *******************************************
c frac_impa
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidfi,start,count,frac_impa2)
#else
        status=NF_GET_VARA_REAL(ncidp,varidfi,start,count,frac_impa2)
#endif
      call gr_ecrit_fi(kev,kon,iim,jjm+1,frac_impa2,frac_impa)

c frac_nucl

#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidfn,start,count,frac_nucl2)
#else
        status=NF_GET_VARA_REAL(ncidp,varidfn,start,count,frac_nucl2)
#endif
      call gr_ecrit_fi(kev,kon,iim,jjm+1,frac_nucl2,frac_nucl)

C**** Vents aux sol ********************************************

      start(3)=irec
      start(4)=0
      count(3)=1
      count(4)=0

c pyu1
	print*,'LECTURE de yu1 a irec =',irec
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidyu1,start,count,pyu12)
#else
        status=NF_GET_VARA_REAL(ncidp,varidyu1,start,count,pyu12)
#endif
      call gr_ecrit_fi(1,kon,iim,jjm+1,pyu12,pyu1)

c pyv1
        print*,'LECTURE de yv1 a irec =',irec
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidyv1,start,count,pyv12)
#else
        status=NF_GET_VARA_REAL(ncidp,varidyv1,start,count,pyv12)
#endif
      call gr_ecrit_fi(1,kon,iim,jjm+1,pyv12,pyv1)

C**** Temerature au sol ********************************************
c ftsol1
        print*,'LECTURE de ftsol1 a irec =',irec
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidfts1,start,count,ftsol12)
#else
      status=NF_GET_VARA_REAL(ncidp,varidfts1,start,count,ftsol12)
#endif
       call gr_ecrit_fi(1,kon,iim,jjm+1,ftsol12,ftsol1)

c ftsol2
        print*,'LECTURE de ftsol2 a irec =',irec
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidfts2,start,count,ftsol22)
#else
      status=NF_GET_VARA_REAL(ncidp,varidfts2,start,count,ftsol22)
#endif
      call gr_ecrit_fi(1,kon,iim,jjm+1,ftsol22,ftsol2)

c ftsol3
	 print*,'LECTURE de ftsol3 a irec =',irec
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidfts3,start,count,ftsol32)
#else
      status=NF_GET_VARA_REAL(ncidp,varidfts3,start,count,ftsol32)
#endif
      call gr_ecrit_fi(1,kon,iim,jjm+1,ftsol32,ftsol3)

c ftsol4
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidfts4,start,count,ftsol42)
#else
      status=NF_GET_VARA_REAL(ncidp,varidfts4,start,count,ftsol42)
#endif
      call gr_ecrit_fi(1,kon,iim,jjm+1,ftsol42,ftsol4)

C**** Nature sol ********************************************
c psrf1 
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidpsr1,start,count,psrf12)
#else
      status=NF_GET_VARA_REAL(ncidp,varidpsr1,start,count,psrf12)
#endif
c      call dump2d(iip1-1,jjm+1,psrf12,'PSRF1NC')
      call gr_ecrit_fi(1,kon,iim,jjm+1,psrf12,psrf1)

c psrf2 
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidpsr2,start,count,psrf22)
#else
      status=NF_GET_VARA_REAL(ncidp,varidpsr2,start,count,psrf22)
#endif
c      call dump2d(iip1-1,jjm+1,psrf22,'PSRF2NC')
      call gr_ecrit_fi(1,kon,iim,jjm+1,psrf22,psrf2)

c psrf3 
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidpsr3,start,count,psrf32)
#else
      status=NF_GET_VARA_REAL(ncidp,varidpsr3,start,count,psrf32)
#endif
      call gr_ecrit_fi(1,kon,iim,jjm+1,psrf32,psrf3)

c psrf4 
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidpsr4,start,count,psrf42)
#else
      status=NF_GET_VARA_REAL(ncidp,varidpsr4,start,count,psrf42)
#endif
      call gr_ecrit_fi(1,kon,iim,jjm+1,psrf42,psrf4)
	
          do i = 1,kon
	
	psrf(i,1) = psrf1(i)
        psrf(i,2) = psrf2(i)
        psrf(i,3) = psrf3(i)
c test abderr
c	print*,'Dans read_pstoke psrf3 =',psrf3(i),i
        psrf(i,4) = psrf4(i)
  
        ftsol(i,1) = ftsol1(i)
        ftsol(i,2) = ftsol2(i)
        ftsol(i,3) = ftsol3(i)
        ftsol(i,4) = ftsol4(i)
	
          enddo
	
	endif
	
	return
	
	end

