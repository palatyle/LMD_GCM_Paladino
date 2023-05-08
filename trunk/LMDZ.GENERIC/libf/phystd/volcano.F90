    SUBROUTINE volcano(nq,ngrid,nlay,volc_loc,nline,zzlev,data_readt,ssource)

    use tracer_h
    use comgeomfi_h
    use datafile_mod, only: datadir
    use mod_phys_LMDZ_mpi_data, only: mpi_rank
    use print_control_mod, only: lunout
    IMPLICIT NONE

    ! Variable Statement

    include "dimensions.h"

       ! Local Variables
       INTEGER :: i, j, l, ii
       INTEGER :: iq        !Tracer ID
       REAL :: mmsource(nq)
       REAL :: asource(nlay,nq) !(M)ATHAM source
       ! REAL, dimension(4,173) :: data_read !173 mparts in matham, 3 tracers. 1st column is height
       ! REAL, dimension(173,4) :: data_readt !transposed (something something fortran columnwise)
       CHARACTER(LEN=20) :: tracername

       ! Inputs:
       INTEGER, intent(in) :: nq               ! Number of tracers
       INTEGER, intent(in) :: ngrid            ! Number of grid points
       INTEGER, intent(in) ::  nlay            ! Number of vertical levels
       INTEGER, intent(in) :: volc_loc
       INTEGER, intent(in) :: nline
       REAL, intent(in) :: zzlev(ngrid,nlay+1) ! height between the layers (m)
       REAL, DIMENSION(nlay+1) :: dlev
       REAL, intent(in) :: data_readt(nline,5)

       !Outputs
       REAL, intent(out) :: ssource(ngrid, nlay, nq) ! (kg.kg-1.s-1)
       ! ----------------------------------------------------------------
       write(lunout,*) "ivolc in volcano", volc_loc
          ! Loop through MATHAM levels and find closest GCM level
          ! If multiple MATHAM entries exist for one GCM level, sum together.
          ! WRITE(lunout,*) 'nlay ', nlay, 'zzlev: ', zzlev(volc_loc,:)
          asource = 0.0 ! Initilize array with zeroes
          dlev = 0.0 ! same for diff array
          ! Loop through MATHAM levels and GCM levels, find index of smallest
          ! value and use that to insert into correct loc of asource.
          DO i = 2,nline !Don't include 0 elevation
            DO j = 2,nlay+1
              dlev(j) = abs((data_readt(i,1)*1000)-zzlev(volc_loc,j)) ! *1000 for km -> m
            END DO
            l = 1
            DO iq = 1, nq
              tracername=noms(iq)
              if (tracername(1:4).eq."volc") then
                l = l + 1 ! this index coresponds to only the volcano tracers.
                          ! An assumption being made here is that the order of tracers
                          ! in the GCM = order of tracers in MATHAM
                asource(minloc(dlev(2:)),iq) = asource(minloc(dlev(2:)),iq) + data_readt(i,l)/1000 !data_readt SHOULD BE DIVIDED BY 1000 to make g/kg -> kg/kg. Keeping like for now though.
              END IF
            END DO
          END DO
          ! do ii = 1,nlay
          !   WRITE(lunout,*) "asource", mpi_rank, asource(ii,:)
          ! end do 
          ! Loop through volc tracers and levels and push MATHAM data into ssource(zdqvolc)
          DO iq=1, nq
            tracername=noms(iq)
            if (tracername(1:4).eq."volc") then ! if first 4 letters of tracer has "volc" in name
              DO l=1, nlay
                  ! write(lunout,*) volc_loc,l,iq
                  ssource(volc_loc,l,iq) = asource(l,iq)
              END DO
            END IF
          END DO
          ! do ii = 1,nlay
          !   WRITE(lunout,*) "ssource", mpi_rank, ssource(volc_loc,ii,:)
          ! end do

END SUBROUTINE volcano
