! $Id$
#ifdef CPP_XIOS
MODULE wxios
    USE xios
    USE iaxis
    USE iaxis_attr
    USE icontext_attr
    USE idate
    USE idomain_attr
    USE ifield_attr
    USE ifile_attr
    USE ixml_tree

    !Variables disponibles pendant toute l'execution du programme:
    
    INTEGER, SAVE :: g_comm
    CHARACTER(len=100), SAVE :: g_ctx_name ="LMDZ"
    TYPE(xios_context), SAVE :: g_ctx
!$OMP THREADPRIVATE(g_comm,g_cts_name,g_ctx)
    LOGICAL, SAVE :: g_flag_xml = .FALSE.
    CHARACTER(len=100) :: g_field_name = "nofield"
!$OMP THREADPRIVATE(g_flag_xml,g_field_name)
    REAL :: missing_val_omp
    REAL :: missing_val
!$OMP THREADPRIVATE(missing_val)

#ifdef XIOS1
#error "XIOS v1 no longer supported, use XIOS v2."
#endif

    CONTAINS
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !   36day => 36d etc     !!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    SUBROUTINE reformadate(odate, ndate)
        CHARACTER(len=*), INTENT(IN) :: odate
        TYPE(xios_duration) :: ndate
        
        INTEGER :: i = 0
         !!!!!!!!!!!!!!!!!!
         ! Pour XIOS:
         !  year : y
         !  month : mo
         !  day : d
         !  hour : h
         !  minute : mi
         !  second : s
         !!!!!!!!!!!!!!!!!!

        i = INDEX(odate, "day")
        IF (i > 0) THEN
            read(odate(1:i-1),*) ndate%day
        END IF

        i = INDEX(odate, "hr")
        IF (i > 0) THEN
            read(odate(1:i-1),*) ndate%hour
        END IF

        i = INDEX(odate, "mth")
        IF (i > 0) THEN
            read(odate(1:i-1),*) ndate%month
        END IF
        
        !IF (prt_level >= 10) WRITE(lunout,*) "Xios. ", odate, " => ", ndate
    END SUBROUTINE reformadate
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !   ave(X) => average etc     !!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    CHARACTER(len=7) FUNCTION reformaop(op)
        CHARACTER(len=*), INTENT(IN) :: op
        
        INTEGER :: i = 0
        reformaop = "average"
        
        IF (op.EQ."inst(X)") THEN
            reformaop = "instant"
        END IF
        
        IF (op.EQ."once") THEN
            reformaop = "once"
        END IF
        
        IF (op.EQ."t_max(X)") THEN
            reformaop = "maximum"
        END IF
        
        IF (op.EQ."t_min(X)") THEN
            reformaop = "minimum"
        END IF
        
        !IF (prt_level >= 10) WRITE(lunout,*) "Xios. ", op, " => ", reformaop
    END FUNCTION reformaop

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Routine d'initialisation      !!!!!!!!!!!!!
    !     A lancer juste après mpi_init !!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE wxios_init(xios_ctx_name, locom, outcom, type_ocean)
        USE print_control_mod, ONLY : prt_level, lunout
        IMPLICIT NONE

      CHARACTER(len=*), INTENT(IN) :: xios_ctx_name
      INTEGER, INTENT(IN), OPTIONAL :: locom
      INTEGER, INTENT(OUT), OPTIONAL :: outcom
      CHARACTER(len=6), INTENT(IN), OPTIONAL :: type_ocean

    
        TYPE(xios_context) :: xios_ctx
        INTEGER :: xios_comm

        IF (prt_level >= 10) WRITE(lunout,*) "wxios_init: Initialization"



        IF (PRESENT(locom)) THEN
          CALL xios_initialize(xios_ctx_name, local_comm = locom, return_comm = xios_comm )
          IF (prt_level >= 10) WRITE(lunout,*) "wxios_init: ctx=",xios_ctx_name," local_comm=",locom,", return_comm=",xios_comm
        ELSE
          CALL xios_initialize(xios_ctx_name, return_comm = xios_comm )
          IF (prt_level >= 10) WRITE(lunout,*) "wxios_init: ctx=",xios_ctx_name," return_comm=",xios_comm
        END IF
        
        IF (PRESENT(outcom)) THEN
          outcom = xios_comm
          IF (prt_level >= 10) WRITE(lunout,*) "wxios_init: ctx=",xios_ctx_name," outcom=",outcom
        END IF
        
        !Enregistrement des variables globales:
        g_comm = xios_comm
        g_ctx_name = xios_ctx_name
        
!        ! Si couple alors init fait dans cpl_init
!        IF (.not. PRESENT(type_ocean)) THEN
!            CALL wxios_context_init()
!        ENDIF

    END SUBROUTINE wxios_init

    SUBROUTINE wxios_context_init()
        USE print_control_mod, ONLY : prt_level, lunout
        USE mod_phys_lmdz_mpi_data, ONLY : COMM_LMDZ_PHY
        IMPLICIT NONE

        TYPE(xios_context) :: xios_ctx

!$OMP MASTER
        !Initialisation du contexte:
        CALL xios_context_initialize(g_ctx_name, COMM_LMDZ_PHY)
        CALL xios_get_handle(g_ctx_name, xios_ctx)    !Récupération
        CALL xios_set_current_context(xios_ctx)            !Activation
        g_ctx = xios_ctx

        IF (prt_level >= 10) THEN
          WRITE(lunout,*) "wxios_context_init: Current context is ",trim(g_ctx_name)
          WRITE(lunout,*) "     now call xios_solve_inheritance()"
        ENDIF
        !Une première analyse des héritages:
        CALL xios_solve_inheritance()
!$OMP END MASTER
    END SUBROUTINE wxios_context_init


    SUBROUTINE wxios_set_context()
        IMPLICIT NONE
        TYPE(xios_context) :: xios_ctx

       !$OMP MASTER
        CALL xios_get_handle(g_ctx_name, xios_ctx)    !Récupération
        CALL xios_set_current_context(xios_ctx)            !Activation
       !$OMP END MASTER

    END SUBROUTINE wxios_set_context

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Routine de paramétrisation !!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE wxios_set_cal(pasdetemps, calendrier, annee, mois, jour, heure, ini_an, ini_mois, ini_jour, ini_heure)
        USE print_control_mod, ONLY : prt_level, lunout
        IMPLICIT NONE

     !Paramètres:
     CHARACTER(len=*), INTENT(IN) :: calendrier
     INTEGER, INTENT(IN) :: annee, mois, jour, ini_an, ini_mois, ini_jour
     REAL, INTENT(IN) :: pasdetemps, heure, ini_heure
     
     !Variables:
     CHARACTER(len=80) :: abort_message
     CHARACTER(len=19) :: date
     INTEGER :: njour = 1
     
     !Variables pour xios:
     TYPE(xios_duration) :: mdtime
     !REAL(kind = 8) :: year=0, month=0, day=0, hour=0, minute=0, second=0
     
        mdtime%second=pasdetemps

        !Réglage du calendrier:
        SELECT CASE (calendrier)
            CASE('earth_360d')
                CALL xios_define_calendar("D360")
                IF (prt_level >= 10) WRITE(lunout,*) 'wxios_set_cal: Calendrier terrestre a 360 jours/an'
            CASE('earth_365d')
                CALL xios_define_calendar("NoLeap")
                IF (prt_level >= 10) WRITE(lunout,*) 'wxios_set_cal: Calendrier terrestre a 365 jours/an'
            CASE('gregorian')
                CALL xios_define_calendar("Gregorian")
                IF (prt_level >= 10) WRITE(lunout,*) 'wxios_set_cal: Calendrier gregorien'
            CASE DEFAULT
                abort_message = 'wxios_set_cal: Mauvais choix de calendrier'
                CALL abort_physic('Gcm:Xios',abort_message,1)
        END SELECT
        
        !Formatage de la date d'origine:
        WRITE(date, "(i4.4,'-',i2.2,'-',i2.2,' ',i2.2,':00:00')") annee, mois, jour, int(heure) 
        
        IF (prt_level >= 10) WRITE(lunout,*) "wxios_set_cal: Time origin: ", date

        CALL xios_set_time_origin(xios_date(annee,mois,jour,int(heure),0,0))

        !Formatage de la date de debut:

        WRITE(date, "(i4.4,'-',i2.2,'-',i2.2,' ',i2.2,':00:00')") ini_an, ini_mois, ini_jour, int(ini_heure)
        
        IF (prt_level >= 10) WRITE(lunout,*) "wxios_set_cal: Start date: ", date
        
        CALL xios_set_start_date(xios_date(ini_an,ini_mois,ini_jour,int(ini_heure),0,0))
        
        !Et enfin,le pas de temps:
        CALL xios_set_timestep(mdtime)
        IF (prt_level >= 10) WRITE(lunout,*) "wxios_set_cal: ts=",mdtime
    END SUBROUTINE wxios_set_cal

    SUBROUTINE wxios_set_timestep(ts)
        REAL, INTENT(IN) :: ts
        TYPE(xios_duration) :: mdtime     

        mdtime%timestep = ts

        CALL xios_set_timestep(mdtime)
    END SUBROUTINE wxios_set_timestep

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Pour initialiser un domaine !!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE wxios_domain_param(dom_id,flip_coordinates)
       USE dimphy, only: klon
       USE mod_phys_lmdz_transfert_para, ONLY: gather, bcast
       USE mod_phys_lmdz_para, only: jj_nb, jj_begin, jj_end, ii_begin, ii_end, &
                                     mpi_size, mpi_rank, klon_mpi, &
                                     is_sequential, is_south_pole_dyn
       USE mod_grid_phy_lmdz, only: nbp_lon, nbp_lat, klon_glo         
       USE print_control_mod, ONLY : prt_level, lunout
       USE geometry_mod

       IMPLICIT NONE
        CHARACTER(len=*),INTENT(IN) :: dom_id ! domain identifier
        LOGICAL,OPTIONAL,INTENT(IN) :: flip_coordinates ! .true. to change signs
                                                        ! of coordinates 
        LOGICAL :: flip_axes=.false. ! default; do not swap axes
        REAL   :: rlat_glo(klon_glo)
        REAL   :: rlon_glo(klon_glo)
        REAL   :: io_lat(nbp_lat)
        REAL   :: io_lon(nbp_lon)
        LOGICAL :: mask(nbp_lon,jj_nb) !Masque pour les problèmes de recouvrement MPI
        TYPE(xios_domain) :: dom
        INTEGER :: i
        LOGICAL :: boool
        

        IF (PRESENT(flip_coordinates)) flip_axes=flip_coordinates
        
        IF (flip_axes) THEN
          ! change signs of axes
          CALL gather(-latitude_deg,rlat_glo) 
          CALL bcast(rlat_glo)
          CALL gather(-longitude_deg,rlon_glo)
          CALL bcast(rlon_glo)
        ELSE
          CALL gather(latitude_deg,rlat_glo) 
          CALL bcast(rlat_glo)
          CALL gather(longitude_deg,rlon_glo)
          CALL bcast(rlon_glo)
        ENDIF
    
  !$OMP MASTER  
        io_lat(1)=rlat_glo(1)
        io_lat(nbp_lat)=rlat_glo(klon_glo)
        IF ((nbp_lon*nbp_lat) > 1) then
          DO i=2,nbp_lat-1
            io_lat(i)=rlat_glo(2+(i-2)*nbp_lon)
          ENDDO
        ENDIF

        IF (klon_glo == 1) THEN
          io_lon(1)=rlon_glo(1)
        ELSE 
          io_lon(1:nbp_lon)=rlon_glo(2:nbp_lon+1)
        ENDIF

        
        !On récupère le handle:
        CALL xios_get_domain_handle(dom_id, dom)
        
        !On parametrise le domaine:
        CALL xios_set_domain_attr_hdl(dom, ni_glo=nbp_lon, ibegin=0, ni=nbp_lon, type="rectilinear")
        CALL xios_set_domain_attr_hdl(dom, nj_glo=nbp_lat, jbegin=jj_begin-1, nj=jj_nb, data_dim=2)
        CALL xios_set_domain_attr_hdl(dom, lonvalue_1d=io_lon(1:nbp_lon), latvalue_1d=io_lat(jj_begin:jj_end))

        IF (.NOT.is_sequential) THEN
            mask(:,:)=.TRUE.
            if (ii_begin>1) mask(1:ii_begin-1,1) = .FALSE.
            if (ii_end<nbp_lon) mask(ii_end+1:nbp_lon,jj_nb) = .FALSE.
            ! special case for south pole
            if ((ii_end==1).and.(is_south_pole_dyn)) mask(1:nbp_lon,jj_nb)=.true.
            IF (prt_level >= 10) THEN
              WRITE(lunout,*) "wxios_domain_param: mpirank=",mpi_rank," mask(:,1)=",mask(:,1)
              WRITE(lunout,*) "wxios_domain_param: mpirank=",mpi_rank," mask(:,jj_nb)=",mask(:,jj_nb)
            ENDIF
            CALL xios_set_domain_attr_hdl(dom, mask_2d=mask)
        END IF

         CALL xios_is_defined_domain_attr_hdl(dom,ni_glo=boool)
        !Vérification:
        IF (xios_is_valid_domain(dom_id)) THEN
            IF (prt_level >= 10) WRITE(lunout,*) "wxios_domain_param: Domain initialized: ", trim(dom_id), boool
        ELSE
            IF (prt_level >= 10) WRITE(lunout,*) "wxios_domain_param: Invalid domain: ", trim(dom_id)
        END IF
!$OMP END MASTER
        
    END SUBROUTINE wxios_domain_param
    

    SUBROUTINE wxios_domain_param_unstructured(dom_id,flip_coordinates)
        USE geometry_mod, ONLY : longitude, latitude, boundslon, boundslat,ind_cell_glo
        USE mod_grid_phy_lmdz, ONLY : nvertex, klon_glo
        USE mod_phys_lmdz_para 
        USE nrtype, ONLY : PI
        IMPLICIT NONE
        CHARACTER(len=*),INTENT(IN) :: dom_id ! domain identifier
        LOGICAL,OPTIONAL,INTENT(IN) :: flip_coordinates ! .true. to change signs
                                                        ! of coordinates 
        LOGICAL :: flip_axes=.false. ! default; do not swap axes
        REAL :: lon_mpi(klon_mpi)
        REAL :: lat_mpi(klon_mpi)
        REAL :: boundslon_mpi(klon_mpi,nvertex)
        REAL :: boundslat_mpi(klon_mpi,nvertex)
        INTEGER :: ind_cell_glo_mpi(klon_mpi)
        TYPE(xios_domaingroup) :: dom

        IF (PRESENT(flip_coordinates)) flip_axes=flip_coordinates
        
        IF (flip_axes) THEN
          ! change signs of axes and boundaries
          CALL gather_omp(-longitude*180/PI,lon_mpi)
          CALL gather_omp(-latitude*180/PI,lat_mpi)
          CALL gather_omp(-boundslon*180/PI,boundslon_mpi)
          CALL gather_omp(-boundslat*180/PI,boundslat_mpi)
        ELSE
          CALL gather_omp(longitude*180/PI,lon_mpi)
          CALL gather_omp(latitude*180/PI,lat_mpi)
          CALL gather_omp(boundslon*180/PI,boundslon_mpi)
          CALL gather_omp(boundslat*180/PI,boundslat_mpi)
        ENDIF
        CALL gather_omp(ind_cell_glo,ind_cell_glo_mpi)
        

!$OMP MASTER
        CALL xios_get_domaingroup_handle(dom_id, dom)
        
        !On parametrise le domaine:
        CALL xios_set_attr(dom, ni_glo=klon_glo, ibegin=ij_begin-1, ni=ij_nb, type="unstructured")
        CALL xios_set_attr(dom, nvertex=nvertex, lonvalue_1d=lon_mpi, latvalue_1d=lat_mpi, &
                           bounds_lon_1d=TRANSPOSE(boundslon_mpi), bounds_lat_1d=TRANSPOSE(boundslat_mpi) )
        CALL xios_set_attr(dom, i_index=ind_cell_glo_mpi(:)-1)
!$OMP END MASTER

    END SUBROUTINE wxios_domain_param_unstructured




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Pour déclarer un axe vertical !!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE wxios_add_vaxis(axis_id, axis_size, axis_value)
        USE print_control_mod, ONLY : prt_level, lunout
        IMPLICIT NONE

        CHARACTER (len=*), INTENT(IN) :: axis_id
        INTEGER, INTENT(IN) :: axis_size
        REAL, DIMENSION(axis_size), INTENT(IN) :: axis_value
        
!        TYPE(xios_axisgroup) :: axgroup
!        TYPE(xios_axis) :: ax
!        CHARACTER(len=50) :: axis_id 
        
!        IF (len_trim(axisgroup_id).gt.len(axis_id)) THEN
!          WRITE(lunout,*) "wxios_add_vaxis: error, size of axis_id too small!!"
!          WRITE(lunout,*) "     increase it to at least ",len_trim(axisgroup_id)
!          CALL abort_gcm("wxios_add_vaxis","len(axis_id) too small",1)
!        ENDIF
!        axis_id=trim(axisgroup_id)
        
        !On récupère le groupe d'axes qui va bien:
        !CALL xios_get_axisgroup_handle(axisgroup_id, axgroup)
        
        !On ajoute l'axe correspondant à ce fichier:
        !CALL xios_add_axis(axgroup, ax, TRIM(ADJUSTL(axis_id)))
        
        !Et on le parametrise:
        !CALL xios_set_axis_attr_hdl(ax, size=axis_size, value=axis_value)
        
        ! Ehouarn: New way to declare axis, without axis_group:
        CALL xios_set_axis_attr(trim(axis_id),n_glo=axis_size,value=axis_value)

        !Vérification:
        IF (xios_is_valid_axis(TRIM(ADJUSTL(axis_id)))) THEN
            IF (prt_level >= 10) WRITE(lunout,*) "wxios_add_vaxis: Axis created: ", TRIM(ADJUSTL(axis_id))
        ELSE
            WRITE(lunout,*) "wxios_add_vaxis: Invalid axis: ", TRIM(ADJUSTL(axis_id))
        END IF

    END SUBROUTINE wxios_add_vaxis
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Pour déclarer un fichier  !!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE wxios_add_file(fname, ffreq, flvl)
        USE print_control_mod, ONLY : prt_level, lunout
        IMPLICIT NONE

        CHARACTER(len=*), INTENT(IN) :: fname
        CHARACTER(len=*), INTENT(IN) :: ffreq
        INTEGER, INTENT(IN) :: flvl
        
        TYPE(xios_file) :: x_file
        TYPE(xios_filegroup) :: x_fg
        TYPE(xios_duration) :: nffreq
        
        !On regarde si le fichier n'est pas défini par XML:
        IF (.NOT.xios_is_valid_file(fname)) THEN
            !On créé le noeud:
            CALL xios_get_filegroup_handle("defile", x_fg)
            CALL xios_add_file(x_fg, x_file, fname)
        
            !On reformate la fréquence:
            CALL reformadate(ffreq, nffreq)
        
            !On configure:
            CALL xios_set_file_attr_hdl(x_file, name="X"//fname,&
                output_freq=nffreq, output_level=flvl, enabled=.TRUE.)
       
            IF (xios_is_valid_file("X"//fname)) THEN
                IF (prt_level >= 10) THEN
                  WRITE(lunout,*) "wxios_add_file: New file: ", "X"//fname
                  WRITE(lunout,*) "wxios_add_file: output_freq=",nffreq,"; output_lvl=",flvl
                ENDIF
            ELSE
                WRITE(lunout,*) "wxios_add_file: Error, invalid file: ", "X"//trim(fname)
                WRITE(lunout,*) "wxios_add_file: output_freq=",nffreq,"; output_lvl=",flvl
            END IF
        ELSE
            IF (prt_level >= 10) THEN
              WRITE(lunout,*) "wxios_add_file: File ",trim(fname), " défined using XML."
            ENDIF
            ! Ehouarn: add an enable=.true. on top of xml definitions... why???
            CALL xios_set_file_attr(fname, enabled=.TRUE.)
        END IF
    END SUBROUTINE wxios_add_file
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Pour créer un champ      !!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE wxios_add_field(fieldname, fieldgroup, fieldlongname, fieldunit)
        USE netcdf, only: nf90_fill_real

        IMPLICIT NONE
        INCLUDE 'iniprint.h'
        
        CHARACTER(len=*), INTENT(IN) :: fieldname
        TYPE(xios_fieldgroup), INTENT(IN) :: fieldgroup
        CHARACTER(len=*), INTENT(IN) :: fieldlongname
        CHARACTER(len=*), INTENT(IN) :: fieldunit
        
        TYPE(xios_field) :: field
        CHARACTER(len=10) :: newunit
        REAL(KIND=8) :: def
        
        !La valeur par défaut des champs non définis:
        def = nf90_fill_real
        
        IF (fieldunit .EQ. " ") THEN
            newunit = "-"
        ELSE
            newunit = fieldunit
        ENDIF
        
        !On ajoute le champ:
        CALL xios_add_field(fieldgroup, field, fieldname)
        !IF (prt_level >= 10) WRITE(lunout,*) "wxios_add_field: ",fieldname,fieldgroup, fieldlongname, fieldunit
        
        !On rentre ses paramètres:
        CALL xios_set_field_attr_hdl(field, standard_name=fieldlongname, unit=newunit, default_value=def)
        IF (prt_level >= 10) WRITE(lunout,*) "wxios_add_field: Field ",trim(fieldname), "cree:"
        IF (prt_level >= 10) WRITE(lunout,*) "wxios_add_field: long_name=",trim(fieldlongname),"; unit=",trim(newunit),";  default_value=",nf90_fill_real

    END SUBROUTINE wxios_add_field
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Pour déclarer un champ      !!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE wxios_add_field_to_file(fieldname, fdim, fid, fname, fieldlongname, fieldunit, field_level, op, nam_axvert)
        USE print_control_mod, ONLY : prt_level, lunout
        IMPLICIT NONE

        CHARACTER(len=*), INTENT(IN) :: fieldname
        INTEGER, INTENT(IN)          :: fdim, fid
        CHARACTER(len=*), INTENT(IN) :: fname
        CHARACTER(len=*), INTENT(IN) :: fieldlongname
        CHARACTER(len=*), INTENT(IN) :: fieldunit
        INTEGER, INTENT(IN)          :: field_level
        CHARACTER(len=*), INTENT(IN) :: op
        
        CHARACTER(len=20) :: axis_id ! Ehouarn: dangerous...
        CHARACTER(len=20), INTENT(IN), OPTIONAL :: nam_axvert
        CHARACTER(len=100) :: operation
        TYPE(xios_file) :: f
        TYPE(xios_field) :: field
        TYPE(xios_fieldgroup) :: fieldgroup
        TYPE(xios_duration) :: freq_op

        LOGICAL :: bool=.FALSE.
        INTEGER :: lvl =0
        
        
        ! Ajout Abd pour NMC:
        IF (fid.LE.6) THEN
          axis_id="presnivs"
        ELSE
          axis_id="plev"
        ENDIF
 
        IF (PRESENT(nam_axvert)) THEN
           axis_id=nam_axvert
           print*,'nam_axvert=',axis_id
        ENDIF
        
        !on prépare le nom de l'opération:
        operation = reformaop(op)
        
        
        !On selectionne le bon groupe de champs:
        IF (fdim.EQ.2) THEN
          CALL xios_get_fieldgroup_handle("fields_2D", fieldgroup)
        ELSE
          CALL xios_get_fieldgroup_handle("fields_3D", fieldgroup)
        ENDIF
        
        !On regarde si le champ à déjà été créé ou non:
        IF (xios_is_valid_field(fieldname) .AND. .NOT. g_field_name == fieldname) THEN
            !Si ce champ existe via XML (ie, dès le premier passage, ie g_field_name != fieldname) alors rien d'autre à faire
            IF (prt_level >= 10) WRITE(lunout,*) "wxios_add_field_to_file: Field ", trim(fieldname), "exists via XML"
            g_flag_xml = .TRUE.
            g_field_name = fieldname

        ELSE IF (.NOT. g_field_name == fieldname) THEN
            !Si premier pssage et champ indéfini, alors on le créé

            IF (prt_level >= 10) WRITE(lunout,*) "wxios_add_field_to_file: Field ", trim(fieldname), "does not exist"
            
            !On le créé:
            CALL wxios_add_field(fieldname,  fieldgroup, fieldlongname, fieldunit)
            IF (xios_is_valid_field(fieldname)) THEN
                IF (prt_level >= 10) WRITE(lunout,*) "wxios_add_field_to_file: Field ", trim(fieldname), "created"
            ENDIF

            g_flag_xml = .FALSE.
            g_field_name = fieldname

        END IF

        IF (.NOT. g_flag_xml) THEN
            !Champ existe déjà, mais pas XML, alors on l'ajoute
            !On ajoute le champ:
            CALL xios_get_file_handle(fname, f)
            CALL xios_add_fieldtofile(f, field)
            
            
            !L'operation, sa frequence:
            freq_op%timestep=1
            CALL xios_set_field_attr_hdl(field, field_ref=fieldname, operation=TRIM(ADJUSTL(operation)), freq_op=freq_op, prec=4)

            
            !On rentre ses paramètres:
            CALL xios_set_field_attr_hdl(field, level=field_level, enabled=.TRUE.)
            
            IF (fdim.EQ.2) THEN
                !Si c'est un champ 2D:
                IF (prt_level >= 10) THEN
                  WRITE(lunout,*) "wxios_add_field_to_file: 2D Field ", trim(fieldname), " in ", "X"//trim(fname) ," configured with:"
                  WRITE(lunout,*) "wxios_add_field_to_file: op=", TRIM(ADJUSTL(operation))
                  WRITE(lunout,*) "wxios_add_field_to_file: freq_op=1ts","; lvl=",field_level
                ENDIF
            ELSE
                !Si 3D :
                !On ajoute l'axe vertical qui va bien:
                CALL xios_set_field_attr_hdl(field, axis_ref=TRIM(ADJUSTL(axis_id)))
                
                IF (prt_level >= 10) THEN
                  WRITE(lunout,*) "wxios_add_field_to_file: 3D Field",trim(fieldname), " in ", "X"//trim(fname), "configured with:"
                  WRITE(lunout,*) "wxios_add_field_to_file: freq_op=1ts","; lvl=",field_level
                  WRITE(lunout,*) "wxios_add_field_to_file: axis=",TRIM(ADJUSTL(axis_id))
                ENDIF
            END IF
        
        ELSE
            !Sinon on se contente de l'activer:
            CALL xios_set_field_attr(fieldname, enabled=.TRUE.)
            !NB: This will override an enable=.false. set by a user in the xml file;
            !   then the only way to not output the field is by changing its
            !   output level
        ENDIF       
        
    END SUBROUTINE wxios_add_field_to_file
    
!    SUBROUTINE wxios_update_calendar(ito)
!        INTEGER, INTENT(IN) :: ito
!        CALL xios_update_calendar(ito)
!    END SUBROUTINE wxios_update_calendar
!    
!    SUBROUTINE wxios_write_2D(fieldname, fdata)
!        CHARACTER(len=*), INTENT(IN) :: fieldname
!        REAL, DIMENSION(:,:), INTENT(IN) :: fdata
!
!        CALL xios_send_field(fieldname, fdata)
!    END SUBROUTINE wxios_write_2D
    
!    SUBROUTINE wxios_write_3D(fieldname, fdata)
!        CHARACTER(len=*), INTENT(IN) :: fieldname
!        REAL, DIMENSION(:,:,:), INTENT(IN) :: fdata
!        
!        CALL xios_send_field(fieldname, fdata)
!    END SUBROUTINE wxios_write_3D
    
    SUBROUTINE wxios_closedef()
        CALL xios_close_context_definition()
!        CALL xios_update_calendar(0)
    END SUBROUTINE wxios_closedef
    
    SUBROUTINE wxios_close()
        CALL xios_context_finalize()
         CALL xios_finalize()
     END SUBROUTINE wxios_close
END MODULE wxios
#endif
