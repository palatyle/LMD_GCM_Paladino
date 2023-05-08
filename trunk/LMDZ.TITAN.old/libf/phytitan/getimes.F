c==========================================================
c                       getimes.F
c==========================================================
c    Pack de routines pour calculer des temps d'execution.
c
c    UTILISATION :
c    L'appel a la routine initimes() doit se faire en 
c    début de programme.
c    Pour calculer le temps d'execution d'une portion de
c    code il suffit de placer la routine begintime devant
c    la portion et la routine endtime après la portion.
c
c    L'argument tps de la routine endtime retourne le 
c    temps d'execution en seconde de la portion de code.
c==========================================================

**************************************************
*  routine begintime                             *
*  Retourne dans la variable nb_ini, le nombre de*
*  cycle au moment de l'appel a begintime.       *
*                                                *
*  ARGUMENTS :                                   *
*  nb_ini (output) : nombre de cycle au lancement*
*                    de la routine.              *
**************************************************
      subroutine begintime(nb_ini)
        integer nb_ini
        CALL SYSTEM_CLOCK(COUNT=nb_ini)
        return
      end


**************************************************
*  routine endtime                               *
*  Retourne dans la variable tps, le temps entre *
*  l'appel a begintime et endtime.               *
*                                                *
*  ARGUMENTS :                                   *
*  nb_ini (input) : nombre de cycle au lancement *
*                   de la routine begintime      *
*  tps (output)   : temps en secondes ecoules    *
*                   entre begintime et endtime   *
**************************************************
       subroutine endtime(nb_ini,tps)
       implicit none
         integer nb_max,nb_sec
         common/horloge/nb_max,nb_sec
         integer nb_ini
         integer nb_fin,tmp
         real tps
         CALL SYSTEM_CLOCK(COUNT=nb_fin)
         tmp =nb_fin-nb_ini
         if (nb_fin.lt.nb_ini) tmp=tmp+nb_max
         tps = FLOAT(tmp)/ nb_sec
         return
       end

**************************************************
*  routine initimes                              *
*  Initialise le compteur                        *
*  Cette routine doit se positionnee en debut de *
*  programme                                     *
**************************************************
       subroutine initimes()
#include "itemps.h"
         integer nb_max,nb_sec
         common/horloge/nb_max,nb_sec
         CALL SYSTEM_CLOCK(COUNT_RATE=nb_sec,
     &                      COUNT_MAX=nb_max)
*        initialisation variables de stockage
         ttdynt   = 0.
         ttadvtr  = 0.
         ttphys   = 0.
         ttmuphys = 0.
         tthaze   = 0.         
         ttcclds  = 0.
         ttsclds  = 0.
         ttrad    = 0.
         ttphytra = 0.
         
         return
       end


**************************************************
*  routine printimes                             *
*  affiche des temps d'execution.                *
*  iout : sortie dans un fichier.                *
*   iout = 0 ---> sortie ecran                   *
*   iout = 1 ---> sortie dans temps.out          *
**************************************************
       subroutine printimes(iout)
#include "itemps.h"
         integer iunit
         logical ok
         ok = .true.
         iunit = 9
         if (iout.eq.0) then
           iunit = 6
         else
           do while (ok)
             iunit = iunit+1
             inquire(unit=iunit,OPENED=ok)
             if (iunit.eq.100) exit
           enddo
             if (iunit.eq.100) then
               print*,"Je n'ai pas trouve d'unite logique libre."
               print*,"J'affiche les temps a l'ecran."
               iunit = 6
             endif
         endif 
         if (iunit.ne.6) open(iunit,file="tpsprog")
         write(iunit,*) "#############################################"
         write(iunit,*) "ttdynt        :",ttdynt
         write(iunit,*) "  ttdyn       :",ttdynt-ttphys
         write(iunit,*) "  ttadvtr     :",ttdyntr
         write(iunit,*) "  ttphys      :",ttphys
         write(iunit,*) "    ttrad     :",tthaze
         write(iunit,*) "    ttphytra  :",ttphytra
         write(iunit,*) "    ttmuphys  :",ttmuphys
         write(iunit,*) "      tthaze  :",tthaze
         write(iunit,*) "      ttcclds :",ttcclds
         write(iunit,*) "      ttsclds :",ttsclds
         write(iunit,*) "#############################################"
         if (iunit.ne.6) close(iunit)
         return
       end

