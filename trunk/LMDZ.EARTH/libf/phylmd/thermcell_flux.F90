!
! $Id: thermcell_flux.F90 1403 2010-07-01 09:02:53Z fairhead $
!


      SUBROUTINE thermcell_flux(ngrid,klev,ptimestep,masse, &
     &       lalim,lmax,alim_star,  &
     &       entr_star,detr_star,f,rhobarz,zlev,zw2,fm,entr,  &
     &       detr,zqla,zmax,lev_out,lunout1,igout)


!---------------------------------------------------------------------------
!thermcell_flux: deduction des flux
!---------------------------------------------------------------------------

      IMPLICIT NONE
#include "iniprint.h"
      
      INTEGER ig,l
      INTEGER ngrid,klev
      
      REAL alim_star(ngrid,klev),entr_star(ngrid,klev)
      REAL detr_star(ngrid,klev)
      REAL zw2(ngrid,klev+1)
      REAL zlev(ngrid,klev+1)
      REAL masse(ngrid,klev)
      REAL ptimestep
      REAL rhobarz(ngrid,klev)
      REAL f(ngrid)
      INTEGER lmax(ngrid)
      INTEGER lalim(ngrid)
      REAL zqla(ngrid,klev)
      REAL zmax(ngrid)

      integer ncorecfm1,ncorecfm2,ncorecfm3,ncorecalpha
      integer ncorecfm4,ncorecfm5,ncorecfm6,ncorecfm7,ncorecfm8
      

      REAL entr(ngrid,klev),detr(ngrid,klev)
      REAL fm(ngrid,klev+1)
      REAL zfm

      integer igout
      integer lev_out
      integer lunout1

      REAL f_old,ddd0,eee0,ddd,eee,zzz

      REAL fomass_max,alphamax
      save fomass_max,alphamax
!$OMP THREADPRIVATE(fomass_max,alphamax)

      character (len=20) :: modname='thermcell_flux'
      character (len=80) :: abort_message

      fomass_max=0.5
      alphamax=0.7

      ncorecfm1=0
      ncorecfm2=0
      ncorecfm3=0
      ncorecfm4=0
      ncorecfm5=0
      ncorecfm6=0
      ncorecfm7=0
      ncorecfm8=0
      ncorecalpha=0

!initialisation
      fm(:,:)=0.
      
      if (prt_level.ge.10) then
         write(lunout,*) 'Dans thermcell_flux 0'
         write(lunout,*) 'flux base ',f(igout)
         write(lunout,*) 'lmax ',lmax(igout)
         write(lunout,*) 'lalim ',lalim(igout)
         write(lunout,*) 'ig= ',igout
         write(lunout,*) ' l E*    A*     D*  '
         write(lunout,'(i4,3e15.5)') (l,entr_star(igout,l),alim_star(igout,l),detr_star(igout,l) &
     &    ,l=1,lmax(igout))
      endif


!-------------------------------------------------------------------------
! Verification de la nullite des entrainement et detrainement au dessus
! de lmax(ig)
!-------------------------------------------------------------------------
     if ( prt_level > 1 ) THEN
      do l=1,klev
         do ig=1,ngrid
            if (l.le.lmax(ig)) then
               if (entr_star(ig,l).gt.1.) then
                    print*,'WARNING thermcell_flux 1 ig,l,lmax(ig)',ig,l,lmax(ig)
                    print*,'entr_star(ig,l)',entr_star(ig,l)
                    print*,'alim_star(ig,l)',alim_star(ig,l)
                    print*,'detr_star(ig,l)',detr_star(ig,l)
               endif
            else
               if (abs(entr_star(ig,l))+abs(alim_star(ig,l))+abs(detr_star(ig,l)).gt.0.) then
                    print*,'cas 1 : ig,l,lmax(ig)',ig,l,lmax(ig)
                    print*,'entr_star(ig,l)',entr_star(ig,l)
                    print*,'alim_star(ig,l)',alim_star(ig,l)
                    print*,'detr_star(ig,l)',detr_star(ig,l)
                    abort_message = ''
                    CALL abort_gcm (modname,abort_message,1)
               endif
            endif
         enddo
      enddo
     endif  !( prt_level > 1 ) THEN
!-------------------------------------------------------------------------
! Multiplication par le flux de masse issu de la femreture
!-------------------------------------------------------------------------

      do l=1,klev
         entr(:,l)=f(:)*(entr_star(:,l)+alim_star(:,l))
         detr(:,l)=f(:)*detr_star(:,l)
      enddo

      if (prt_level.ge.10) then
         write(lunout,*) 'Dans thermcell_flux 1'
         write(lunout,*) 'flux base ',f(igout)
         write(lunout,*) 'lmax ',lmax(igout)
         write(lunout,*) 'lalim ',lalim(igout)
         write(lunout,*) 'ig= ',igout
         write(lunout,*) ' l   E    D     W2'
         write(lunout,'(i4,3e15.5)') (l,entr(igout,l),detr(igout,l) &
     &    ,zw2(igout,l+1),l=1,lmax(igout))
      endif

      fm(:,1)=0.
      do l=1,klev
         do ig=1,ngrid
            if (l.lt.lmax(ig)) then
               fm(ig,l+1)=fm(ig,l)+entr(ig,l)-detr(ig,l)
            elseif(l.eq.lmax(ig)) then
               fm(ig,l+1)=0.
               detr(ig,l)=fm(ig,l)+entr(ig,l)
            else
               fm(ig,l+1)=0.
            endif
         enddo
      enddo



! Test provisoire : pour comprendre pourquoi on corrige plein de fois 
! le cas fm6, on commence par regarder une premiere fois avant les
! autres corrections.

      do l=1,klev
         do ig=1,ngrid
            if (detr(ig,l).gt.fm(ig,l)) then
               ncorecfm8=ncorecfm8+1
!              igout=ig
            endif
         enddo
      enddo

      if (prt_level.ge.10) &
    &    call printflux(ngrid,klev,lunout,igout,f,lmax,lalim, &
    &    ptimestep,masse,entr,detr,fm,'2  ')



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FH Version en cours de test;
! par rapport a thermcell_flux, on fait une grande boucle sur "l"
! et on modifie le flux avec tous les contr�les appliques d'affilee
! pour la meme couche
! Momentanement, on duplique le calcule du flux pour pouvoir comparer
! les flux avant et apres modif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do l=1,klev

         do ig=1,ngrid
            if (l.lt.lmax(ig)) then
               fm(ig,l+1)=fm(ig,l)+entr(ig,l)-detr(ig,l)
            elseif(l.eq.lmax(ig)) then
               fm(ig,l+1)=0.
               detr(ig,l)=fm(ig,l)+entr(ig,l)
            else
               fm(ig,l+1)=0.
            endif
         enddo


!-------------------------------------------------------------------------
! Verification de la positivite des flux de masse
!-------------------------------------------------------------------------

!     do l=1,klev
         do ig=1,ngrid
            if (fm(ig,l+1).lt.0.) then
!              print*,'fm1<0',l+1,lmax(ig),fm(ig,l+1)
                ncorecfm1=ncorecfm1+1
               fm(ig,l+1)=fm(ig,l)
               detr(ig,l)=entr(ig,l)
            endif
         enddo
!     enddo

      if (prt_level.ge.10) &
     &   write(lunout,'(i4,4e14.4)') l,masse(igout,l)/ptimestep, &
     &     entr(igout,l),detr(igout,l),fm(igout,l+1)

!-------------------------------------------------------------------------
!Test sur fraca croissant
!-------------------------------------------------------------------------


      if (1.eq.1) then
!     do l=1,klev
         do ig=1,ngrid
          if (l.ge.lalim(ig).and.l.le.lmax(ig) &
     &    .and.(zw2(ig,l+1).gt.1.e-10).and.(zw2(ig,l).gt.1.e-10) ) then
!  zzz est le flux en l+1 a frac constant
             zzz=fm(ig,l)*rhobarz(ig,l+1)*zw2(ig,l+1)  &
     &                          /(rhobarz(ig,l)*zw2(ig,l))
             if (fm(ig,l+1).gt.zzz) then
                detr(ig,l)=detr(ig,l)+fm(ig,l+1)-zzz
                fm(ig,l+1)=zzz
                ncorecfm4=ncorecfm4+1
             endif
          endif
        enddo
!     enddo

      if (prt_level.ge.10) &
     &   write(lunout,'(i4,4e14.4)') l,masse(igout,l)/ptimestep, &
     &     entr(igout,l),detr(igout,l),fm(igout,l+1)
      else
       if (l.eq.1) then
         print*,'Test sur les fractions croissantes inhibe dans thermcell_flux2'
       endif
      endif


!-------------------------------------------------------------------------
!test sur flux de masse croissant
!-------------------------------------------------------------------------

!     do l=1,klev
         do ig=1,ngrid
            if ((fm(ig,l+1).gt.fm(ig,l)).and.(l.gt.lalim(ig))) then
              f_old=fm(ig,l+1)
              fm(ig,l+1)=fm(ig,l)
              detr(ig,l)=detr(ig,l)+f_old-fm(ig,l+1)
               ncorecfm5=ncorecfm5+1
            endif
         enddo
!     enddo

      if (prt_level.ge.10) &
     &   write(lunout,'(i4,4e14.4)') l,masse(igout,l)/ptimestep, &
     &     entr(igout,l),detr(igout,l),fm(igout,l+1)

!-------------------------------------------------------------------------
!detr ne peut pas etre superieur a fm
!-------------------------------------------------------------------------

      if(1.eq.1) then

!     do l=1,klev
         do ig=1,ngrid
            if (entr(ig,l)<0.) then
               print*,'N1 ig,l,entr',ig,l,entr(ig,l)
               abort_message = 'entr negatif'
               CALL abort_gcm (modname,abort_message,1)
            endif
            if (detr(ig,l).gt.fm(ig,l)) then
               ncorecfm6=ncorecfm6+1
               detr(ig,l)=fm(ig,l)
!              entr(ig,l)=fm(ig,l+1)

! Dans le cas ou on est au dessus de la couche d'alimentation et que le
! detrainement est plus fort que le flux de masse, on stope le thermique.
               if (l.gt.lalim(ig)) then
                  lmax(ig)=l
                  fm(ig,l+1)=0.
                  entr(ig,l)=0.
               else
                  ncorecfm7=ncorecfm7+1
               endif
            endif

            if(l.gt.lmax(ig)) then
               detr(ig,l)=0.
               fm(ig,l+1)=0.
               entr(ig,l)=0.
            endif

            if (entr(ig,l).lt.0.) then
               print*,'ig,l,lmax(ig)',ig,l,lmax(ig)
               print*,'entr(ig,l)',entr(ig,l)
               print*,'fm(ig,l)',fm(ig,l)
               abort_message = 'probleme dans thermcell flux'
               CALL abort_gcm (modname,abort_message,1)
            endif
         enddo
!     enddo
      endif


      if (prt_level.ge.10) &
     &   write(lunout,'(i4,4e14.4)') l,masse(igout,l)/ptimestep, &
     &     entr(igout,l),detr(igout,l),fm(igout,l+1)

!-------------------------------------------------------------------------
!fm ne peut pas etre negatif
!-------------------------------------------------------------------------

!     do l=1,klev
         do ig=1,ngrid
            if (fm(ig,l+1).lt.0.) then
               detr(ig,l)=detr(ig,l)+fm(ig,l+1)
               fm(ig,l+1)=0.
!              print*,'fm2<0',l+1,lmax(ig)
               ncorecfm2=ncorecfm2+1
            endif
            if (detr(ig,l).lt.0.) then
               print*,'cas 2 : ig,l,lmax(ig)',ig,l,lmax(ig)
               print*,'detr(ig,l)',detr(ig,l)
               print*,'fm(ig,l)',fm(ig,l)
               abort_message = 'probleme dans thermcell flux'
               CALL abort_gcm (modname,abort_message,1)
            endif
        enddo
!    enddo

      if (prt_level.ge.10) &
     &   write(lunout,'(i4,4e14.4)') l,masse(igout,l)/ptimestep, &
     &     entr(igout,l),detr(igout,l),fm(igout,l+1)

!-----------------------------------------------------------------------
!la fraction couverte ne peut pas etre superieure a 1            
!-----------------------------------------------------------------------


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FH Partie a revisiter.
! Il semble qu'etaient codees ici deux optiques dans le cas
! F/ (rho *w) > 1
! soit limiter la hauteur du thermique en considerant que c'est 
! la derniere chouche, soit limiter F a rho w.
! Dans le second cas, il faut en fait limiter a un peu moins
! que ca parce qu'on a des 1 / ( 1 -alpha) un peu plus loin
! dans thermcell_main et qu'il semble de toutes facons deraisonable
! d'avoir des fractions de 1..
! Ci dessous, et dans l'etat actuel, le premier des  deux if est
! sans doute inutile.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!    do l=1,klev
        do ig=1,ngrid
           if (zw2(ig,l+1).gt.1.e-10) then
           zfm=rhobarz(ig,l+1)*zw2(ig,l+1)*alphamax
           if ( fm(ig,l+1) .gt. zfm) then
              f_old=fm(ig,l+1)
              fm(ig,l+1)=zfm
!             zw2(ig,l+1)=0.
!             zqla(ig,l+1)=0.
              detr(ig,l)=detr(ig,l)+f_old-fm(ig,l+1)
!             lmax(ig)=l+1
!             zmax(ig)=zlev(ig,lmax(ig))
!             print*,'alpha>1',l+1,lmax(ig)
              ncorecalpha=ncorecalpha+1
           endif
           endif
        enddo
!    enddo
!


      if (prt_level.ge.10) &
     &   write(lunout,'(i4,4e14.4)') l,masse(igout,l)/ptimestep, &
     &     entr(igout,l),detr(igout,l),fm(igout,l+1)

! Fin de la grande boucle sur les niveaux verticaux
      enddo

      if (prt_level.ge.10) &
    &    call printflux(ngrid,klev,lunout,igout,f,lmax,lalim, &
    &    ptimestep,masse,entr,detr,fm,'8  ')


!-----------------------------------------------------------------------
! On fait en sorte que la quantite totale d'air entraine dans le 
! panache ne soit pas trop grande comparee a la masse de la maille
!-----------------------------------------------------------------------

      if (1.eq.1) then
      do l=1,klev-1
         do ig=1,ngrid
            eee0=entr(ig,l)
            ddd0=detr(ig,l)
            eee=entr(ig,l)-masse(ig,l)*fomass_max/ptimestep
            ddd=detr(ig,l)-eee
            if (eee.gt.0.) then
                ncorecfm3=ncorecfm3+1
                entr(ig,l)=entr(ig,l)-eee
                if ( ddd.gt.0.) then
!   l'entrainement est trop fort mais l'exces peut etre compense par une
!   diminution du detrainement)
                   detr(ig,l)=ddd
                else
!   l'entrainement est trop fort mais l'exces doit etre compense en partie
!   par un entrainement plus fort dans la couche superieure
                   if(l.eq.lmax(ig)) then
                      detr(ig,l)=fm(ig,l)+entr(ig,l)
                   else
                      if(l.ge.lmax(ig).and.0.eq.1) then
                         print*,'ig,l',ig,l
                         print*,'eee0',eee0
                         print*,'ddd0',ddd0
                         print*,'eee',eee
                         print*,'ddd',ddd
                         print*,'entr',entr(ig,l)
                         print*,'detr',detr(ig,l)
                         print*,'masse',masse(ig,l)
                         print*,'fomass_max',fomass_max
                         print*,'masse(ig,l)*fomass_max/ptimestep',masse(ig,l)*fomass_max/ptimestep
                         print*,'ptimestep',ptimestep
                         print*,'lmax(ig)',lmax(ig)
                         print*,'fm(ig,l+1)',fm(ig,l+1)
                         print*,'fm(ig,l)',fm(ig,l)
                         abort_message = 'probleme dans thermcell_flux'
                         CALL abort_gcm (modname,abort_message,1)
                      endif
                      entr(ig,l+1)=entr(ig,l+1)-ddd
                      detr(ig,l)=0.
                      fm(ig,l+1)=fm(ig,l)+entr(ig,l)
                      detr(ig,l)=0.
                   endif
                endif
            endif
         enddo
      enddo
      endif
!                  
!              ddd=detr(ig)-entre
!on s assure que tout s annule bien en zmax
      do ig=1,ngrid
         fm(ig,lmax(ig)+1)=0.
         entr(ig,lmax(ig))=0.
         detr(ig,lmax(ig))=fm(ig,lmax(ig))+entr(ig,lmax(ig))
      enddo

!-----------------------------------------------------------------------
! Impression du nombre de bidouilles qui ont ete necessaires
!-----------------------------------------------------------------------

      if (ncorecfm1+ncorecfm2+ncorecfm3+ncorecfm4+ncorecfm5+ncorecalpha > 0 ) then
       if (prt_level.ge.10) then
          print*,'PB thermcell : on a du coriger ',ncorecfm1,'x fm1',&
    &     ncorecfm2,'x fm2',ncorecfm3,'x fm3 et', &
    &     ncorecfm4,'x fm4',ncorecfm5,'x fm5 et', &
    &     ncorecfm6,'x fm6', &
    &     ncorecfm7,'x fm7', &
    &     ncorecfm8,'x fm8', &
    &     ncorecalpha,'x alpha'
       endif
      endif

      if (prt_level.ge.10) &
    &    call printflux(ngrid,klev,lunout,igout,f,lmax,lalim, &
    &    ptimestep,masse,entr,detr,fm,'fin')


      return
      end

      subroutine printflux(ngrid,klev,lunout,igout,f,lmax,lalim, &
    &    ptimestep,masse,entr,detr,fm,descr)

     implicit none

      integer ngrid,klev,lunout,igout,l,lm

      integer lmax(klev),lalim(klev)
      real ptimestep,masse(ngrid,klev),entr(ngrid,klev),detr(ngrid,klev)
      real fm(ngrid,klev+1),f(ngrid)

      character*3 descr

      character (len=20) :: modname='thermcell_flux'
      character (len=80) :: abort_message

      lm=lmax(igout)+5
      if(lm.gt.klev) lm=klev

      print*,'Impression jusque lm=',lm

         write(lunout,*) 'Dans thermcell_flux '//descr
         write(lunout,*) 'flux base ',f(igout)
         write(lunout,*) 'lmax ',lmax(igout)
         write(lunout,*) 'lalim ',lalim(igout)
         write(lunout,*) 'ig= ',igout
         write(lunout,'(a3,4a14)') 'l','M','E','D','F'
         write(lunout,'(i4,4e14.4)') (l,masse(igout,l)/ptimestep, &
     &     entr(igout,l),detr(igout,l) &
     &    ,fm(igout,l+1),l=1,lm)


      do l=lmax(igout)+1,klev
          if (abs(entr(igout,l))+abs(detr(igout,l))+abs(fm(igout,l)).gt.0.) then
          print*,'cas 1 : igout,l,lmax(igout)',igout,l,lmax(igout)
          print*,'entr(igout,l)',entr(igout,l)
          print*,'detr(igout,l)',detr(igout,l)
          print*,'fm(igout,l)',fm(igout,l)
          abort_message = ''
          CALL abort_gcm (modname,abort_message,1)
          endif
      enddo

      return
      end

