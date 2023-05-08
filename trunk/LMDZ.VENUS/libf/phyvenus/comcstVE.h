!-----------------------------------------------------------------------
! INCLUDE comcstVE.h

        integer nnuve,nbmat
        parameter (nnuve=68)  ! fichiers Vincent et Bullock
!       parameter (nnuve=598)  ! fichiers Vincent et Bullock
        parameter (nbmat=220) ! Max number of matrixes in Vincent's file

      common/comcstVE/al,bl,nlatve,indexve,nbpsve,nbszave,               &
     & psurfve,szave

      real   al(nnuve),bl(nnuve)     ! for Planck luminances calculations
! Structure of the ksi matrixes
      integer nlatve,indexve(5),nbpsve(5),nbszave(5)
      real   psurfve(16,5)           ! surface pressure in matrixes (Pa)
      real   szave(8,5)              ! converted in mu0


!-----------------------------------------------------------------------
