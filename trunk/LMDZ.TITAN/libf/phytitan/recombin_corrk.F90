SUBROUTINE recombin_corrk(q,ip,it)

  !     ===================================================================
  !     Purpose
  !     -------
  !     Recombine correlated-k of individual species in a unique corr-k.
  !
  !       -> See Lacis and Oinas (1991), Amundsen et al (2016).
  !
  !     One important hypothesis is that we assume random overlap
  !     Steps are - 1. Recombining ( Convolution spectra of two species )
  !                 2. Resorting the ngauss*ngauss obtained data
  !                 3. Rebin down to a smaller amount of bin ( here ngauss )
  !                 4. Start again if another specie have to be included
  !
  !     In each band computation is done only for species with enough optical depth.
  !     
  !     Authors
  !     -------
  !     J. Vatant d'Ollone (2018)
  !     ====================================================================

  USE radinc_h
  USE radcommon_h, only: gasi, gasv, gasi_recomb, gasv_recomb, &
       gweight, pqrold, permut_idx, tlimit, w_cum
  USE sort_mod, only: qsort, isort

  !-----------------------------------------------------------------------
  !     Declarations:
  !     -------------

  IMPLICIT NONE

  !  Arguments :
  !  -----------
  REAL, DIMENSION(L_REFVAR), INTENT(IN)  :: q

  INTEGER,                   INTENT(IN)  :: ip
  INTEGER,                   INTENT(IN)  :: it

  !  Local variables :
  !  -----------------
  INTEGER :: iw, ispec, ig, jg, ind, ibin

  REAL, DIMENSION(L_NGAUSS)             :: krecomb

  REAL, DIMENSION(L_NGAUSS*L_NGAUSS)    :: ktwospec
  REAL, DIMENSION(L_NGAUSS*L_NGAUSS)    :: ktwospec_s
  REAL, DIMENSION(L_NGAUSS*L_NGAUSS)    :: wtwospec
  REAL, DIMENSION(L_NGAUSS*L_NGAUSS)    :: wtwospec_s
  REAL, DIMENSION(L_NGAUSS*L_NGAUSS)    :: wtwospec_cum  

  REAL, DIMENSION(L_REFVAR,L_NGAUSS)    :: tmpk

  REAL :: wsplit


  ! -------------------
  ! I. INFRARED
  ! -------------------

  DO iw=1,L_NSPECTI

     krecomb(:) = q(1)*gasi(it,ip,1,iw,:) ! init for > 1 loop and works also if only one active specie

     tmpk(:,:) = gasi(it,ip,:,iw,:)

     IF ( L_REFVAR .GT. 1 ) THEN

        DO ispec=2,L_REFVAR

           ! Save ( a lot of ) CPU time, we don't add the specie if negligible absorption in the band
           IF ( q(ispec)*tmpk(ispec,L_NGAUSS-1) .LE. tlimit ) CYCLE
           !              IF ( tmpk(ispec,L_NGAUSS-1) .LE. tlimit ) CYCLE

           ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           ! 1. Recombining ~~~~~~~~~~~~~~~~~~~~~~~~~
           ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           DO ig=1,L_NGAUSS
              DO jg=1, L_NGAUSS
                 ind = jg+(ig-1)*L_NGAUSS
                 ktwospec(ind) = krecomb(ig)+q(ispec)*tmpk(ispec,jg)
                 wtwospec(ind) = gweight(ig)*gweight(jg)
              ENDDO
           ENDDO

           ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           ! 2. Resorting ~~~~~~~~~~~~~~~~~~~~~~~
           ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           
           ! Pre-sort from last step ( we have always a similar regular pattern ) to gain time for sorting
           ! NB : quite small array, quicker to permut with 2nd array than in place
           DO ind=1,L_NGAUSS*L_NGAUSS
              ktwospec_s(ind) = ktwospec(permut_idx(ind)) ! NB : won't do anything at firstcall
           ENDDO

           CALL isort(ktwospec_s,L_NGAUSS*L_NGAUSS,permut_idx)  ! Insertion sort quicker because pre-sorted
           !CALL qsort(ktwospec_s,L_NGAUSS*L_NGAUSS,permut_idx) ! Quicksort slower for pre-sorted

           ! Sort w according to permutations of k.
           ! NB : quite small array, quicker to permut with 2nd array than in place
           DO ind=1,L_NGAUSS*L_NGAUSS
              wtwospec_s(ind) = wtwospec(permut_idx(ind))
           ENDDO

           wtwospec_cum(1) = wtwospec_s(1)
           DO ind=2,L_NGAUSS*L_NGAUSS
              wtwospec_cum(ind)= wtwospec_cum(ind-1)+wtwospec_s(ind)
           ENDDO

           ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           ! 3. Rebinning on smaller amount of Gauss points ~~~~~~~~~~~
           ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           ibin=1

           krecomb(:)=0.0

           DO ig=1, L_NGAUSS-1

              DO ind=ibin,L_NGAUSS*L_NGAUSS ! rather than a while   
                 IF ( wtwospec_cum(ind) .GT. w_cum(ig) ) THEN
                    wsplit =  w_cum(ig) - wtwospec_cum(ind-1)
                    krecomb(ig)   = krecomb(ig)                                            &
                         + sum ( wtwospec_s(ibin:ind-1)*ktwospec_s(ibin:ind-1) ) &
                         + wsplit*ktwospec_s(ind)
                    krecomb(ig+1) = (wtwospec_s(ind)-wsplit)*ktwospec_s(ind)
                    ibin=ind+1
                    EXIT
                 ENDIF
              ENDDO

              krecomb(L_NGAUSS) = krecomb(L_NGAUSS) + sum ( wtwospec_s(ibin:)*ktwospec_s(ibin:) )

           ENDDO

           krecomb(1:L_NGAUSS-1) =  krecomb(1:L_NGAUSS-1) / gweight(1:L_NGAUSS-1) ! gw(L_NGAUSS)=0

        ENDDO ! ispec=2,L_REFVAR

     ENDIF ! if L_REFVAR .GT. 1

     gasi_recomb(it,ip,iw,:) = krecomb(:)

  ENDDO ! iw=1,L_NSPECTI


  ! ----------------
  ! II. VISIBLE
  ! ----------------

  DO iw=1,L_NSPECTV

     krecomb(:) = q(1)*gasv(it,ip,1,iw,:) ! init for > 1 loop and works also if only one active specie
     
     tmpk(:,:) = gasv(it,ip,:,iw,:)

     IF ( L_REFVAR .GT. 1 ) THEN

        DO ispec=2,L_REFVAR

           ! Save ( a lot of ) CPU time, we don't add the specie if negligible absorption in the band
           IF ( q(ispec)*tmpk(ispec,L_NGAUSS-1) .LE. tlimit ) CYCLE
           !              IF ( tmpk(ispec,L_NGAUSS-1) .LE. tlimit ) CYCLE

           ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           ! 1. Recombining ~~~~~~~~~~~~~~~~~~~~~~~~~
           ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           DO ig=1,L_NGAUSS
              DO jg=1, L_NGAUSS
                 ind = jg+(ig-1)*L_NGAUSS
                 ktwospec(ind) = krecomb(ig)+q(ispec)*tmpk(ispec,jg)
                 wtwospec(ind) = gweight(ig)*gweight(jg)
              ENDDO
           ENDDO

           ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           ! 2. Resorting ~~~~~~~~~~~~~~~~~~~~~~~
           ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

           ! Pre-sort from last call ( we have always a similar regular pattern ) to gain time in sorting
           ! NB : quite small array, quicker to permut with 2nd array than in place
           DO ind=1,L_NGAUSS*L_NGAUSS
              ktwospec_s(ind) = ktwospec(permut_idx(ind))
           ENDDO

           CALL isort(ktwospec_s,L_NGAUSS*L_NGAUSS,permut_idx)  ! Insertion sort quicker because pre-sorted
           !CALL qsort(ktwospec_s,L_NGAUSS*L_NGAUSS,permut_idx) ! Quicksort slower for pre-sorted

           ! Sort w according to permutations of k.
           ! NB : quite small array, quicker to permut with copy than in place
           DO ind=1,L_NGAUSS*L_NGAUSS
              wtwospec_s(ind) = wtwospec(permut_idx(ind))
           ENDDO

           wtwospec_cum(1) = wtwospec_s(1)
           DO ind=2,L_NGAUSS*L_NGAUSS
              wtwospec_cum(ind)= wtwospec_cum(ind-1)+wtwospec_s(ind)
           ENDDO

           ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           ! 3. Rebinning on smaller amount of Gauss points ~~~~~~~~~~~
           ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           ibin=1

           krecomb(:)=0.0

           DO ig=1, L_NGAUSS-1

              DO ind=ibin,L_NGAUSS*L_NGAUSS ! rather than a while   
                 IF ( wtwospec_cum(ind) .GT. w_cum(ig) ) THEN
                    wsplit =  w_cum(ig) - wtwospec_cum(ind-1)
                    krecomb(ig)   = krecomb(ig)                                             &
                         + sum ( wtwospec_s(ibin:ind-1)*ktwospec_s(ibin:ind-1) ) &
                         + wsplit*ktwospec_s(ind)
                    krecomb(ig+1) = (wtwospec_s(ind)-wsplit)*ktwospec_s(ind)
                    ibin=ind+1
                    EXIT
                 ENDIF
              ENDDO

              krecomb(L_NGAUSS) = krecomb(L_NGAUSS) + sum ( wtwospec_s(ibin:)*ktwospec_s(ibin:) )

           ENDDO

           krecomb(1:L_NGAUSS-1) =  krecomb(1:L_NGAUSS-1) / gweight(1:L_NGAUSS-1) ! gw(L_NGAUSS)=0

        ENDDO ! ispec=2,L_REFVAR

     ENDIF ! if L_REFVAR .GT. 1

     gasv_recomb(it,ip,iw,:) = krecomb(:)

  ENDDO

  ! Update saved mixing ratios
  pqrold(ip,:) = q(:)

END SUBROUTINE recombin_corrk

