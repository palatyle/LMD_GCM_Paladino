!
! $Id: ini_wake.F 1403 2010-07-01 09:02:53Z fairhead $
!
      SUBROUTINE INI_WAKE(wape,fip,it_wape_prescr,
     :     wape_prescr, fip_prescr, alp_bl_prescr, ale_bl_prescr)
***************************************************************
*                                                             *
*        INI_WAKE : variables d'initialisation de la poche    *
*                   froide, necessaires au declenchement      *
*                   de la convection.                         *
*                                                             *
*                                                             *
***************************************************************
c Arguments
c =========
c Input
c -----
c   wape           : valeur de l'energie potentielle de la poche (WAPE)
c                    dans l'etat initial
c   fip            : valeur de la puissance incidente sur le front (FIP)
c                    dans l'etat initial
c Output
c ------
c   it_wape_prescr : nombre de pas de temps pendant lesquels la WAPE
c            doit etre imposee.
c   wape_prescr    : valeur prescrite de la WAPE.
c   fip_prescr     : valeur prescrite de la FIP.
c   ale_bl_prescr  : valeur prescrite de la Ale de PBL.
c   alp_bl_prescr  : valeur prescrite de la Alp de PBL.
c
c Variables internes
c ==================
c   it = nbre de pas de temps lu
c   w  = WAPE lue
c   f  = FIP lue
c   alebl  = Ale de PBL lue
c   alpbl  = Alp de PBL lue
c
      include 'iniprint.h'
cdeclarations
      real ale_bl_prescr
      real alp_bl_prescr
      real it

! FH A mettre si besoin dans physiq.def
! FH : voir avec JYG
      it=0.
      w=4.
      f=0.1
      alebl=4.
      alpbl=0.1
c
cCR: on rajoute ale et alp de la PBL precrits
      open (99,file='ini_wake_param.data',form='formatted',
     s      status='old',err=902)
      read (99,*) it
      read (99,*) w
      read (99,*) f
      read (99,*,end=901) alebl
      read (99,*,end=901) alpbl
901   close (99)
902   continue
c
      write(lunout,*)' it,wape ',it,wape
      it_wape_prescr = it
      if (w .lt. 0) then
         wape_prescr = wape
         fip_prescr = fip
      else
         wape_prescr = w
         fip_prescr = f
      endif
c
      write(lunout,*)' alebl, alpbl ',alebl,alpbl
      ale_bl_prescr=alebl
      alp_bl_prescr=alpbl
      print *,'Initialisation de la poche : WAPE, FIP imposees ='
     $               ,wape_prescr, fip_prescr
      print *, '                   pendant ',it_wape_prescr,' steps'
c
      print *,'Initialisation de la BL: ALP, ALE imposees ='
     $               ,alp_bl_prescr, ale_bl_prescr
      return
      end
