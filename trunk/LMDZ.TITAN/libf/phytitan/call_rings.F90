subroutine call_rings(ngrid, ptime, pday, diurnal)
    ! A subroutine to compute the day fraction in case of rings shadowing.

    use radcommon_h, only: eclipse
    use comsaison_h, only: fract, declin
    use comcstfi_mod, only: rad, pi
    use comdiurn_h, only: coslat, sinlat, coslon, sinlon
    use callkeys_mod, only: flatten

    INTEGER, INTENT(IN) :: ngrid
    REAL, INTENT(IN) :: ptime ! "universal time", given as fraction of sol (e.g.: 0.5 for noon)
    REAL, INTENT(IN) :: pday  ! Number of days counted from the North. Spring equinoxe.	
    LOGICAL, INTENT(IN) :: diurnal
!    REAL, DIMENSION(:), INTENT(INOUT) :: fract ! day fraction for each point of the planet

!   to compute the daily average of rings shadowing
    INTEGER, PARAMETER :: nb_hours = 1536 ! set how many times per day are used
    REAL :: pas
    INTEGER :: m
    REAL :: ptime_day ! Universal time in sol fraction 
    REAL:: tmp_zls,tmp_dist_star, tmp_declin, tmp_right_ascen   ! tmp solar longitude, stellar dist, declin and RA
    REAL :: ztim1, ztim2, ztim3
    REAL, DIMENSION(:), ALLOCATABLE :: tmp_fract ! day fraction of the time interval 
    REAL, DIMENSION(:), ALLOCATABLE :: tmp_mu0 ! equivalent solar angle

!! Eclipse incoming sunlight (e.g. Saturn ring shadowing)
    ALLOCATE(eclipse(ngrid))

    write(*,*) 'Rings shadow activated'
        
    if(diurnal .eqv. .false.) then ! we need to compute the daily average insolation (day fraction) 
        pas = 1./nb_hours
        ptime_day = 0.
        fract(:) = 0.
        ALLOCATE(tmp_fract(ngrid))
        ALLOCATE(tmp_mu0(ngrid))
        tmp_fract(:) = 0.
        eclipse(:) = 0.
        tmp_mu0(:) = 0.
                    
        do m=1, nb_hours
            ptime_day = m*pas
            call stellarlong(pday+ptime_day,tmp_zls)
            call orbite(tmp_zls,tmp_dist_star,tmp_declin,tmp_right_ascen)
            
            ztim1=SIN(tmp_declin)
            ztim2=COS(tmp_declin)*COS(2.*pi*(pday+ptime_day-.5))
            ztim3=-COS(tmp_declin)*SIN(2.*pi*(pday+ptime_day-.5))

            call stelang(ngrid,sinlon,coslon,sinlat,coslat,    &
                        ztim1,ztim2,ztim3,tmp_mu0,tmp_fract, flatten)       
            call rings(ngrid, tmp_declin, ptime_day, rad, flatten, eclipse)
            fract(:) = fract(:) + (1.-eclipse(:))*tmp_fract(:) !! fract takes into account the rings shadow and the day/night alternation

        enddo        
     
        fract(:) = fract(:)/nb_hours

        DEALLOCATE(tmp_fract)
        DEALLOCATE(tmp_mu0)
                 
     else   ! instant insolation is weighted by the rings shadow 
            call rings(ngrid, declin, ptime, rad, 0., eclipse)
            fract(:) = fract(:) * (1.-eclipse)
    endif

    IF (ALLOCATED(eclipse)) DEALLOCATE(eclipse)

end subroutine call_rings
