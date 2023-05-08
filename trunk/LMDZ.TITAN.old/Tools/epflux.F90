subroutine  epflux(iip1,jjp1,llm,indefini,latdeg,rbar &
                  ,teta,u3d,v3d,w3d,press &
                  ,epy2d,epz2d,divep2d,vtem2d,wtem2d,acc_meridien2d &
!                 ,vpupbar,wpupbar,vpvpbar,wpvpbar,vptetapbar,wptetapbar &
                  )

      IMPLICIT NONE
!=======================================================================
!
!   Audrey Crespin  Sept 2007
!
!   source: Francois Forget    Avril 1996
!
!
! MODIF SLebonnois nov 2009
!
! Cette subroutine Calcule Les flux d'Eliassen Palm a partir
! des composantes INTERPOLE EN NIVEAU DE PRESSION
! Toutes les variables sont sur la grille de pression press(1:llm)
! On suit la methode de Peixoto et Oort
!
! On est dans le plan meridien
!
! Calcule la divergence & composantes du flux d'EP
! Calcule la deviation a la circulation meridienne moyenne
! Calcule les TEM
! Calcule les termes de Reynolds
!
!=======================================================================
!-----------------------------------------------------------------------
!   declarations:
!   -------------

include "planet.h"

!  --------
!  ARGUMENTS
!  ---------

! Inputs: 
      
      integer :: iip1,jjp1,llm
      real :: indefini,latdeg(jjp1),rbar(jjp1,llm)
      REAL :: teta(iip1,jjp1,llm)
      REAL :: u3d(iip1,jjp1,llm),v3d(iip1,jjp1,llm)
      REAL :: w3d(iip1,jjp1,llm)
      REAL :: press(llm)

! Outputs: 

      REAL :: epy2d(jjp1,llm),epz2d(jjp1,llm),divep2d(jjp1,llm)
      REAL :: vtem2d(jjp1,llm),wtem2d(jjp1,llm),acc_meridien2d(jjp1,llm)
      REAL :: vpupbar(jjp1,llm),wpupbar(jjp1,llm)
      REAL :: wpvpbar(jjp1,llm),vpvpbar(jjp1,llm)
      REAL :: vptetapbar(jjp1,llm),wptetapbar(jjp1,llm)

! Variables non zonales (vpup signifie : v^prim * u^prim)
! -------------------
      REAL :: vpup(iip1,jjp1,llm), wpup(iip1,jjp1,llm) 
      REAL :: vpvp(iip1,jjp1,llm), wpvp(iip1,jjp1,llm)
      REAL :: vptetap(iip1,jjp1,llm), wptetap(iip1,jjp1,llm)

!  Moyennes zonales
!  -------------
      REAL :: ubar(jjp1,llm),vbar(jjp1,llm),tetabar(jjp1,llm)
      REAL :: wbar(jjp1,llm)


      REAL :: dtetabardp(jjp1,llm),dubardp(jjp1,llm)
      REAL :: vz(jjp1,llm),wlat(jjp1,llm),dvzdp(jjp1,llm)
      REAL :: dwlatdlat(jjp1,llm)

      REAL :: depzdp(jjp1,llm)

!  Autres
!  ------
      real :: rlatu(jjp1) ! lat in radians (beware of rounding effects at poles)
      real :: pi
      REAL :: tetas(llm)
      logical :: peixoto

      REAL :: f(jjp1), ducosdlat(jjp1,llm)
      REAL :: depycosdlat
      INTEGER :: i,j,l, n, iim, jjm

!-----------------------------------------------------------------------

      iim = iip1-1
      jjm = jjp1-1
      pi = 2.*asin(1.)

! To avoid rounding effects at poles
      rlatu(:)=latdeg(:)*pi/180.00001

!     Calcul de moyennes zonales
!     --------------------------- 
      call moyzon(iim,jjp1,llm,indefini,u3d,ubar)
      call moyzon(iim,jjp1,llm,indefini,v3d,vbar)
      call moyzon(iim,jjp1,llm,indefini,teta,tetabar)
      call moyzon(iim,jjp1,llm,indefini,w3d,wbar)


! --*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--
      peixoto = .false.

      if (peixoto) then

!     MODIF speciale Peixoto: on utilise la moyenne 
!     globale de teta: tetas(l) a la place de tetabar
      do l=1,llm
         tetas(l) = 0
         n= 0 
         do j=2,jjm
            if (tetabar(j,l).ne.indefini) then 
               tetas(l) = tetas(l) + tetabar(j,l)
               n = n+1 
            end if
         end do
         if (n.eq.0) stop 'bug dans elias '
        tetas(l) = tetas(l) / float(n)

        do j=1,jjp1
             tetabar(j,l) = tetas(l)
        end do
      end do
!       write (*,*) 'tetas(l) ' , tetas
!       write(*,*)

      endif
! --*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--

! coriolis
! --------
      do j=1,jjp1
            f(j) = 2*omega*sin(rlatu(j))
      enddo


!     Calcul des termes non zonaux
!     ---------------------------- 
      do l=1,llm
         do j=1,jjp1
           do i=1,iip1
             if ((  v3d(i,j,l).eq.indefini).or. &
                 (   vbar(j,l).eq.indefini).or. & 
                 ( teta(i,j,l).eq.indefini).or. &
                 (tetabar(j,l).eq.indefini)) then 
               vptetap(i,j,l)= indefini
             else 
               vptetap(i,j,l)=(v3d(i,j,l)-vbar(j,l)) &
                  *(teta(i,j,l)-tetabar(j,l))
             end if
             if ((  w3d(i,j,l).eq.indefini).or. &
                 (   wbar(j,l).eq.indefini).or. &
                 ( teta(i,j,l).eq.indefini).or. &
                 (tetabar(j,l).eq.indefini)) then 
               wptetap(i,j,l)= indefini
             else 
               wptetap(i,j,l)=(w3d(i,j,l)-wbar(j,l)) &
                  *(teta(i,j,l)-tetabar(j,l))
             end if
             if ((v3d(i,j,l).eq.indefini).or. &
                 ( vbar(j,l).eq.indefini).or. &
                 (u3d(i,j,l).eq.indefini).or. &
                 ( ubar(j,l).eq.indefini)) then 
               vpup(i,j,l)= indefini
             else 
               vpup(i,j,l)=(v3d(i,j,l)-vbar(j,l))*(u3d(i,j,l)-ubar(j,l))
             end if
             if ((w3d(i,j,l).eq.indefini).or. &
                 ( wbar(j,l).eq.indefini).or. &
                 (u3d(i,j,l).eq.indefini).or. &
                 ( ubar(j,l).eq.indefini)) then 
               wpup(i,j,l)= indefini
             else 
               wpup(i,j,l)=(w3d(i,j,l)-wbar(j,l))*(u3d(i,j,l)-ubar(j,l))
             end if
             if ((v3d(i,j,l).eq.indefini).or. &
                 ( vbar(j,l).eq.indefini)) then 
               vpvp(i,j,l)= indefini
             else 
               vpvp(i,j,l)=(v3d(i,j,l)-vbar(j,l))*(v3d(i,j,l)-vbar(j,l))
             end if
             if ((w3d(i,j,l).eq.indefini).or. &
                 ( wbar(j,l).eq.indefini).or. &
                 (v3d(i,j,l).eq.indefini).or. &
                 ( vbar(j,l).eq.indefini)) then 
               wpvp(i,j,l)= indefini
             else 
               wpvp(i,j,l)=(w3d(i,j,l)-wbar(j,l))*(v3d(i,j,l)-vbar(j,l))
             end if
           end do
         end do
      end do

!     Moyennes zonales des termes non zonaux
!     -------------------------------------- 

!     Termes de Reynolds

!     flux zonaux de quantite de mouvement (vpupbar et wpupbar)
!     flux meridiens de quantite de mouvement (vpvpbar et wpvpbar)
!     flux meridien et vertical de chaleur (vptetapbar et wptetapbar)

      call moyzon(iim,jjp1,llm,indefini,vptetap,vptetapbar)
      call moyzon(iim,jjp1,llm,indefini,wptetap,wptetapbar)
      call moyzon(iim,jjp1,llm,indefini,vpup,vpupbar)
      call moyzon(iim,jjp1,llm,indefini,wpup,wpupbar)
      call moyzon(iim,jjp1,llm,indefini,vpvp,vpvpbar)
      call moyzon(iim,jjp1,llm,indefini,wpvp,wpvpbar)

!     write(*,*) 'vptetapbar(2,23) =' , vptetapbar(2,23)
!     write(*,*) 'wptetapbar(2,23) =' , wptetapbar(2,23)
!     write(*,*) 'vpupbar(2,23) =' , vpupbar(2,23)

 
!     Derivees d/dp des moyennes zonales
!     --------------------------------------

      call dx_dp(jjp1,llm,indefini,press,ubar,dubardp)
      call dx_dp(jjp1,llm,indefini,press,tetabar,dtetabardp)
!     write (*,*) 'dtetabardp(l) (K/Pa)',(dtetabardp(6,l),l=1,llm)
!       write(*,*)

!     Controle (on triche !) de dtetabardp
      do l=1,llm
         do j=1,jjp1
            if ((dtetabardp(j,l).gt.-1.e-3).and. &
               (dtetabardp(j,l).ne.indefini)) then
               dtetabardp(j,l) = -1.e-3
!              write(*,*) 'profil presque instable en j = ',j,' l= ',l
            end if
         end do
      end do
!     write(*,*)'dtetabardp(l) (K/Pa) corr ',(dtetabardp(6,l),l=1,llm)
!       write(*,*)


!     calculs intermediaires
!     ----------------------
      do l=1,llm

       if ((ubar(2,l).ne.indefini).and.(ubar(1,l).ne.indefini))then
      ducosdlat(1,l) = ( ubar(2,l)*cos(rlatu(2))-ubar(1,l)*cos(rlatu(1)) ) & 
                      /( rlatu(2) - rlatu(1)) 
       else
      ducosdlat(1,l) = indefini
       end if
        do j=2,jjm
       if ((ubar(j+1,l).ne.indefini).and.(ubar(j-1,l).ne.indefini))then
      ducosdlat(j,l) = ( ubar(j+1,l)*cos(rlatu(j+1))-ubar(j-1,l)*cos(rlatu(j-1))) &
                      /( rlatu(j+1) - rlatu(j-1))
       else
      ducosdlat(j,l) = indefini
       end if
        enddo
       if ((ubar(jjp1,l).ne.indefini).and.(ubar(jjm,l).ne.indefini))then
      ducosdlat(jjp1,l) = ( ubar(jjp1,l)*cos(rlatu(jjp1))  &
                               -ubar(jjm,l)*cos(rlatu(jjm))  ) &
                         /( rlatu(jjp1) - rlatu(jjm))
       else
      ducosdlat(jjp1,l) = indefini
       end if
           
        do j=1,jjp1
       if ((vptetapbar(j,l).ne.indefini).and. &
           (dtetabardp(j,l).ne.indefini)) then
           vz(j,l) = vptetapbar(j,l)/dtetabardp(j,l)
         wlat(j,l) = cos(rlatu(j))*vz(j,l)
       else
           vz(j,l) = indefini
         wlat(j,l) = indefini
       end if
        enddo

      enddo
             
!      Calcul des vecteurs flux Eliassen Palm et divergence
!      -----------------------------------------------------
! ref: Read 1986

      do l=1,llm
        epy2d(1,l) = indefini
        epz2d(1,l) = indefini
        do j=2,jjm
 
!     Composante y
!     ------------

            if (   (rbar(j,l).ne.indefini).and. &
                (vpupbar(j,l).ne.indefini).and. &
                (dubardp(j,l).ne.indefini).and. &
                     (vz(j,l).ne.indefini)) then
           
     epy2d(j,l) = rbar(j,l) * cos(rlatu(j))* ( vpupbar(j,l) &
                                   - dubardp(j,l) * vz(j,l) & ! terme non geo
                                             )       

!     if ((j.eq.1).and.(l.eq.23))write(*,*) 'epy2d(1,23) =' , epy2d(1,23)
!     if ((j.eq.2).and.(l.eq.23))write(*,*) 'epy2d(2,23) =' , epy2d(2,23)
!     if ((j.eq.3).and.(l.eq.23))write(*,*) 'epy2d(3,23) =' , epy2d(3,23)
           else 
     epy2d(j,l) = indefini
            end if

!     Composante z
!     ------------

            if (   (rbar(j,l).ne.indefini).and. &
                (wpupbar(j,l).ne.indefini).and. &
              (ducosdlat(j,l).ne.indefini).and. &
                     (vz(j,l).ne.indefini)) then

     epz2d(j,l) = rbar(j,l) * cos(rlatu(j))*                    &
          (  vz(j,l) * ducosdlat(j,l)/(rbar(j,l)*cos(rlatu(j))) & ! terme non geo
           - vz(j,l) * f(j)                                     &
           + wpupbar(j,l)                                       & 
          )
           else
     epz2d(j,l) = indefini
           end if

!       if ((j.eq.2).and.(l.eq.23))write(*,*) 'epz2d(2,23) =' , epz2d(2,23)
        end do
        epy2d(jjp1,l) = indefini
        epz2d(jjp1,l) = indefini
      end do

!    Moyennes euleriennes transformees (vtem2d et wtem2d)
!    -----------------------------------------------
     
      call dx_dp(jjp1,llm,indefini,press,vz,dvzdp)

      do l=1,llm

          if ( (wlat(2,l).ne.indefini) .and. &
               (wlat(1,l).ne.indefini) ) then
               dwlatdlat(1,l)=(wlat(2,l)-wlat(1,l)) &
                     /(rlatu(2)-rlatu(1))
          else
               dwlatdlat(1,l)=indefini
          end if
        do j=2,jjm
          if ( (wlat(j+1,l).ne.indefini) .and. &
               (wlat(j-1,l).ne.indefini) ) then
               dwlatdlat(j,l)=(wlat(j+1,l)-wlat(j-1,l)) &
                     /(rlatu(j+1)-rlatu(j-1))
          else
               dwlatdlat(j,l)=indefini
          end if
        enddo 
          if ( (wlat(jjp1,l).ne.indefini) .and. & 
               (wlat(jjm, l).ne.indefini) ) then
               dwlatdlat(jjp1,l)=(wlat(jjp1,l)-wlat(jjm,l)) &
                     /(rlatu(jjp1)-rlatu(jjm))
          else
               dwlatdlat(jjp1,l)=indefini
          end if

        do j=1,jjp1
          if ( (dvzdp(j,l).ne.indefini) &
           .and.(vbar(j,l).ne.indefini)) then
            vtem2d(j,l) = vbar(j,l)-dvzdp(j,l)
          else
            vtem2d(j,l) = indefini
          endif
        enddo

        do j=1,jjp1
          if (   (rbar(j,l).ne.indefini).and. &
            (dwlatdlat(j,l).ne.indefini).and. &
                 (wbar(j,l).ne.indefini)) then
           wtem2d(j,l) = wbar(j,l)+dwlatdlat(j,l)/(rbar(j,l)*cos(rlatu(j)))
          else
           wtem2d(j,l) = indefini
          endif
        enddo

      enddo

!      write(*,*) 'vtem2d(2,23) =' , vtem2d(2,23)
!      write(*,*) 'wtem2d(2,23) =' , wtem2d(2,23)

!      print*,"OK" 

!    Deviation par rapport a la circulation meridienne moyenne 
!    --------------------------------------------------------
!                    U*.Nabla ubar
      do l=1,llm
        do j=1,jjp1
          if (   (rbar(j,l).ne.indefini) &
          .and.(vtem2d(j,l).ne.indefini) &
          .and.(wtem2d(j,l).ne.indefini) &
       .and.(ducosdlat(j,l).ne.indefini) &
         .and.(dubardp(j,l).ne.indefini)) then
       acc_meridien2d(j,l) =  &
           vtem2d(j,l)*(ducosdlat(j,l)/(rbar(j,l)*cos(rlatu(j))) &
                       -f(j)) &
         + wtem2d(j,l)*dubardp(j,l)     
          else
       acc_meridien2d(j,l) = indefini
          endif
        enddo
      enddo       

!      write(*,*) 'acc_meridien2d(2,23) =' , acc_meridien2d(2,23)

!     Divergence du flux 
!     ------------------
      call dx_dp(jjp1,llm,indefini,press,epz2d,depzdp)
!     write(*,*) 'depzdp(2,23)',depzdp(j,l)

      do l=1,llm
        divep2d(1,l) =0
        do j=2,jjm
          if (    (rbar(j,l).ne.indefini) .and. &
               (epy2d(j+1,l).ne.indefini) .and. &
               (epy2d(j-1,l).ne.indefini) .and. &
                (depzdp(j,l).ne.indefini)    ) then

           depycosdlat=(epy2d(j+1,l)*cos(rlatu(j+1))   &
                      - epy2d(j-1,l)*cos(rlatu(j-1)) ) &
                       /(rlatu(j+1) - rlatu(j-1))

!      DIVERGENCE DU FLUX :

           divep2d(j,l) = depycosdlat/(rbar(j,l)*cos(rlatu(j)))+depzdp(j,l) 

!          if ((j.eq.2).and.(l.eq.23)) then
!              write(*,*) 'depycosdlat(2,23)', depycosdlat
!              write(*,*) 'depzdp(2,23)',depzdp(j,l)
!              write(*,*) 'divergence(2,23)',divep2d(j,l)
!          end if

! Wave driving (Du/Dt et non Dm/dt) :
           divep2d(j,l) = divep2d(j,l)/(rbar(j,l)*cos(rlatu(j))) 
          else
            divep2d(j,l) = indefini
          end if
        end do
        divep2d(jjp1,l) =0
      end do
 
!     write(*,*) ' fin :divep2d(2,23) =' , divep2d(2,23)

! Preparation pour sortie graphique (flux ``chapeau'' Peixoto et Oort p 391)
! A VOIR...
!      do l=1,llm
!        do j=1,jjp1
!           if (epy2d(j,l).ne.indefini) &
!           epy2d(j,l) = 2*pi*a0    *(1./g0)*cos(rlatu(j))*epy2d(j,l)
!           if (epz2d(j,l).ne.indefini) &
!           epz2d(j,l) = 2*pi*a0**2 *(1./g0)*cos(rlatu(j))*epz2d(j,l)
!        end do
!      end do

      return
      end 
