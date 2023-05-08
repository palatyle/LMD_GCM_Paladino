!
! $Id: read_pstoke.F 1403 2010-07-01 09:02:53Z fairhead $
!
c
c
	subroutine read_pstoke(irec,
     .   zrec,zklono,zklevo,airefi,phisfi,
     .   t,mfu,mfd,en_u,de_u,en_d,de_d,coefh,
     .   fm_therm,en_therm,
     .   frac_impa,frac_nucl,pyu1,pyv1,ftsol,psrf)

C******************************************************************************
C  Frederic HOURDIN, Abderrahmane IDELKADI
C Lecture des parametres physique stockes online necessaires pour
C recalculer offline le transport de traceurs sur une grille 2x plus fine que 
C celle online 
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
	
	  integer klono,klevo,imo,jmo
	  parameter (imo=iim/2,jmo=(jjm+1)/2)
	  parameter(klono=(jmo-1)*imo+2,klevo=llm)
	  REAL phisfi(klono)
          REAL phisfi2(imo,jmo+1),airefi2(imo,jmo+1)

          REAL mfu(klono,klevo), mfd(klono,klevo)
          REAL en_u(klono,klevo), de_u(klono,klevo)
          REAL en_d(klono,klevo), de_d(klono,klevo)
          REAL coefh(klono,klevo)
           REAL fm_therm(klono,klevo),en_therm(klono,klevo)

          REAL mfu2(imo,jmo+1,klevo), mfd2(imo,jmo+1,klevo)
          REAL en_u2(imo,jmo+1,klevo), de_u2(imo,jmo+1,klevo)
          REAL en_d2(imo,jmo+1,klevo), de_d2(imo,jmo+1,klevo)
          REAL coefh2(imo,jmo+1,klevo)
           REAL fm_therm2(imo,jmo+1,klevo)
           REAL en_therm2(imo,jmo+1,klevo)

          REAL pl(klevo)
          integer irec
          integer xid,yid,zid,tid
          real zrec,zklono,zklevo,zim,zjm
          integer ncrec,ncklono,ncklevo,ncim,ncjm

          real airefi(klono)
          character*20 namedim

c  !! attention !!
c attention il y a aussi le pb de def klono
c dim de phis??
	  
	 
          REAL frac_impa(klono,klevo), frac_nucl(klono,klevo)
          REAL frac_impa2(imo,jmo+1,klevo), 
     .     frac_nucl2(imo,jmo+1,klevo)
          REAL pyu1(klono), pyv1(klono)
          REAL pyu12(imo,jmo+1), pyv12(imo,jmo+1)
          REAL ftsol(klono,nbsrf)
          REAL psrf(klono,nbsrf) 
	  REAL ftsol1(klono),ftsol2(klono),ftsol3(klono),ftsol4(klono)
          REAL psrf1(klono),psrf2(klono),psrf3(klono),psrf4(klono)
          REAL ftsol12(imo,jmo+1),ftsol22(imo,jmo+1),
     .     ftsol32(imo,jmo+1),
     .     ftsol42(imo,jmo+1)
          REAL psrf12(imo,jmo+1),psrf22(imo,jmo+1),psrf32(imo,jmo+1),
     .     psrf42(imo,jmo+1)
		REAL t(klono,klevo)
		REAL t2(imo,jmo+1,klevo)	
	  integer ncidp
          save ncidp
		integer varidt
          integer varidmfu, varidmfd, varidps, varidenu, variddeu	
          integer varidend,varidded,varidch,varidfi,varidfn
           integer varidfmth,varidenth
          integer varidyu1,varidyv1,varidpl,varidai,varididvt
          integer varidfts1,varidfts2,varidfts3,varidfts4
          integer varidpsr1,varidpsr2,varidpsr3,varidpsr4
          save varidmfu, varidmfd, varidps, varidenu, variddeu
          save varidend,varidded,varidch,varidfi,varidfn
           save varidfmth,varidenth
          save varidyu1,varidyv1,varidpl,varidai,varididvt
          save varidfts1,varidfts2,varidfts3,varidfts4
          save varidpsr1,varidpsr2,varidpsr3,varidpsr4
		save varidt

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

c A FAIRE: Es-il necessaire de stocke t?
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
	
c abder (pour thermiques)
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
            status = nf_inq_dim(ncidp,zid,namedim,ncklevo)
            status = nf_inq_dim(ncidp,tid,namedim,ncrec)
	
            zrec=ncrec
            zklevo=ncklevo
            zim=ncim
            zjm=ncjm
	
	    zklono=zim*(zjm-2)+2
	
	    write(*,*) 'read_pstoke : zrec = ', zrec
            write(*,*) 'read_pstoke : zklevo = ', zklevo
            write(*,*) 'read_pstoke : zim = ', zim 
            write(*,*) 'read_pstoke : zjm = ', zjm
            write(*,*) 'read_pstoke : zklono = ', zklono

c niveaux de pression
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidpl,1,zklevo,pl)
#else
      status=NF_GET_VARA_REAL(ncidp,varidpl,1,zklevo,pl)
#endif

c lecture de aire et phis
	
      start(1)=1
      start(2)=1
      start(3)=1
      start(4)=0

      count(1)=zim
      count(2)=zjm
      count(3)=1
      count(4)=0

c phis
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidps,start,count,phisfi2)
#else
      status=NF_GET_VARA_REAL(ncidp,varidps,start,count,phisfi2)
#endif
      call gr_ecrit_fi(1,klono,imo,jmo+1,phisfi2,phisfi)

c aire
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidai,start,count,airefi2)
#else
      status=NF_GET_VARA_REAL(ncidp,varidai,start,count,airefi2)
#endif
       call gr_ecrit_fi(1,klono,imo,jmo+1,airefi2,airefi)
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
      count(3)=zklevo
      count(4)=1


C *** Lessivage******************************************************
c frac_impa 
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidfi,start,count,frac_impa2)
#else
      status=NF_GET_VARA_REAL(ncidp,varidfi,start,count,frac_impa2)
#endif
      call gr_ecrit_fi(klevo,klono,imo,jmo+1,frac_impa2,frac_impa)

c frac_nucl 
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidfn,start,count,frac_nucl2)
#else
      status=NF_GET_VARA_REAL(ncidp,varidfn,start,count,frac_nucl2)
#endif
      call gr_ecrit_fi(klevo,klono,imo,jmo+1,frac_nucl2,frac_nucl)

C*** Temperature ******************************************************
c abder t
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidt,start,count,t2)
#else
      status=NF_GET_VARA_REAL(ncidp,varidt,start,count,t2)
#endif
      call gr_ecrit_fi(klevo,klono,imo,jmo+1,t2,t)

C*** Flux pour le calcul de la convection TIEDTK ***********************
c mfu
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidmfu,start,count,mfu2)
#else
      status=NF_GET_VARA_REAL(ncidp,varidmfu,start,count,mfu2)
#endif
      call gr_ecrit_fi(klevo,klono,imo,jmo+1,mfu2,mfu)

c mfd
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidmfd,start,count,mfd2)
#else
      status=NF_GET_VARA_REAL(ncidp,varidmfd,start,count,mfd2)
#endif
      call gr_ecrit_fi(klevo,klono,imo,jmo+1,mfd2,mfd)

c en_u 
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidenu,start,count,en_u2)
#else
      status=NF_GET_VARA_REAL(ncidp,varidenu,start,count,en_u2)
#endif
      call gr_ecrit_fi(klevo,klono,imo,jmo+1,en_u2,en_u)

c de_u 
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,variddeu,start,count,de_u2)
#else
      status=NF_GET_VARA_REAL(ncidp,variddeu,start,count,de_u2)
#endif
      call gr_ecrit_fi(klevo,klono,imo,jmo+1,de_u2,de_u)

c en_d 
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidend,start,count,en_d2)
#else
      status=NF_GET_VARA_REAL(ncidp,varidend,start,count,en_d2)
#endif
      call gr_ecrit_fi(klevo,klono,imo,jmo+1,en_d2,en_d)

c de_d 
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidded,start,count,de_d2)
#else
      status=NF_GET_VARA_REAL(ncidp,varidded,start,count,de_d2)
#endif
      call gr_ecrit_fi(klevo,klono,imo,jmo+1,de_d2,de_d)

C **** Coeffecient du mellange turbulent**********************************
c coefh 
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidch,start,count,coefh2)
#else
      status=NF_GET_VARA_REAL(ncidp,varidch,start,count,coefh2)
#endif
       call gr_ecrit_fi(klevo,klono,imo,jmo+1,coefh2,coefh)

C*** Flux ascendant et entrant pour les Thermiques************************
cabder thermiques
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidfmth,start,count,fm_therm2)
#else
      status=NF_GET_VARA_REAL(ncidp,varidfmth,start,count,fm_therm2)
#endif
      call gr_ecrit_fi(klevo,klono,imo,jmo+1,fm_therm2,fm_therm)

#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidenth,start,count,en_therm2)
#else
      status=NF_GET_VARA_REAL(ncidp,varidenth,start,count,en_therm2)
#endif
      call gr_ecrit_fi(klevo,klono,imo,jmo+1,en_therm2,en_therm)

C*** Vitesses aux sol ******************************************************
      start(3)=irec
      start(4)=0
      count(3)=1
      count(4)=0
c pyu1
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidyu1,start,count,pyu12)
#else
      status=NF_GET_VARA_REAL(ncidp,varidyu1,start,count,pyu12)
#endif
      call gr_ecrit_fi(1,klono,imo,jmo+1,pyu12,pyu1)

c pyv1
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidyv1,start,count,pyv12)
#else
      status=NF_GET_VARA_REAL(ncidp,varidyv1,start,count,pyv12)
#endif
      call gr_ecrit_fi(1,klono,imo,jmo+1,pyv12,pyv1)

C*** Temperature au sol ********************************************
c ftsol1
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidfts1,start,count,ftsol12)
#else
      status=NF_GET_VARA_REAL(ncidp,varidfts1,start,count,ftsol12)
#endif
       call gr_ecrit_fi(1,klono,imo,jmo+1,ftsol12,ftsol1)

c ftsol2
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidfts2,start,count,ftsol22)
#else
      status=NF_GET_VARA_REAL(ncidp,varidfts2,start,count,ftsol22)
#endif
      call gr_ecrit_fi(1,klono,imo,jmo+1,ftsol22,ftsol2)

c ftsol3
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidfts3,start,count,ftsol32)
#else
      status=NF_GET_VARA_REAL(ncidp,varidfts3,start,count,ftsol32)
#endif
      call gr_ecrit_fi(1,klono,imo,jmo+1,ftsol32,ftsol3)

c ftsol4
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidfts4,start,count,ftsol42)
#else
      status=NF_GET_VARA_REAL(ncidp,varidfts4,start,count,ftsol42)
#endif
      call gr_ecrit_fi(1,klono,imo,jmo+1,ftsol42,ftsol4)

C*** Nature du sol **************************************************
c psrf1 
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidpsr1,start,count,psrf12)
#else
      status=NF_GET_VARA_REAL(ncidp,varidpsr1,start,count,psrf12)
#endif
      call gr_ecrit_fi(1,klono,imo,jmo+1,psrf12,psrf1)

c psrf2 
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidpsr2,start,count,psrf22)
#else
      status=NF_GET_VARA_REAL(ncidp,varidpsr2,start,count,psrf22)
#endif
      call gr_ecrit_fi(1,klono,imo,jmo+1,psrf22,psrf2)

c psrf3 
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidpsr3,start,count,psrf32)
#else
      status=NF_GET_VARA_REAL(ncidp,varidpsr3,start,count,psrf32)
#endif
      call gr_ecrit_fi(1,klono,imo,jmo+1,psrf32,psrf3)

c psrf4 
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidp,varidpsr4,start,count,psrf42)
#else
      status=NF_GET_VARA_REAL(ncidp,varidpsr4,start,count,psrf42)
#endif
      call gr_ecrit_fi(1,klono,imo,jmo+1,psrf42,psrf4)
	
          do i = 1,klono
	
	psrf(i,1) = psrf1(i)
        psrf(i,2) = psrf2(i)
        psrf(i,3) = psrf3(i)
        psrf(i,4) = psrf4(i)
  
        ftsol(i,1) = ftsol1(i)
        ftsol(i,2) = ftsol2(i)
        ftsol(i,3) = ftsol3(i)
        ftsol(i,4) = ftsol4(i)
	
          enddo
	
	endif
	
	return
	
	end

