MODULE iostart

PRIVATE
    INTEGER,SAVE :: nid_start 
    INTEGER,SAVE :: nid_restart
    
    INTEGER,SAVE :: idim1,idim2,idim3,idim4
    INTEGER,PARAMETER :: length=100
    
    INTERFACE get_field
      MODULE PROCEDURE Get_field_r1,Get_field_r2,Get_field_r3
    END INTERFACE get_field
    
    INTERFACE get_var
      MODULE PROCEDURE get_var_r0,Get_var_r1,Get_var_r2,Get_var_r3
    END INTERFACE get_var

    INTERFACE put_field
      MODULE PROCEDURE put_field_r1,put_field_r2,put_field_r3
    END INTERFACE put_field

    INTERFACE put_var
      MODULE PROCEDURE put_var_r0,put_var_r1,put_var_r2,put_var_r3
    END INTERFACE put_var

    PUBLIC get_field,get_var,put_field,put_var
    PUBLIC open_startphy,close_startphy,open_restartphy,close_restartphy
    
CONTAINS

  SUBROUTINE open_startphy(filename)
  USE netcdf
  USE mod_phys_lmdz_para
  IMPLICIT NONE
    CHARACTER(LEN=*) :: filename
    INTEGER          :: ierr

    IF (is_mpi_root .AND. is_omp_root) THEN
      ierr = NF90_OPEN (filename, NF90_NOWRITE,nid_start)
      IF (ierr.NE.NF90_NOERR) THEN
        write(6,*)' Pb d''ouverture du fichier '//filename
        write(6,*)' ierr = ', ierr
        CALL ABORT
      ENDIF
    ENDIF
   
  END SUBROUTINE open_startphy

  SUBROUTINE Close_startphy
  USE netcdf
  USE mod_phys_lmdz_para
  IMPLICIT NONE
    INTEGER          :: ierr

    IF (is_mpi_root .AND. is_omp_root) THEN
        ierr = NF90_CLOSE (nid_start)
    ENDIF

  END SUBROUTINE close_startphy


  FUNCTION Inquire_Field(Field_name)
  USE netcdf
  USE mod_phys_lmdz_para
  IMPLICIT NONE
    CHARACTER(LEN=*) :: Field_name
    LOGICAL :: inquire_field
    INTEGER :: varid
    INTEGER :: ierr
    
    IF (is_mpi_root .AND. is_omp_root) THEN
      ierr=NF90_INQ_VARID(nid_start,Field_name,varid)
      IF (ierr==NF90_NOERR) THEN
        Inquire_field=.TRUE.
      ELSE
        Inquire_field=.FALSE.
      ENDIF
    ENDIF

    CALL bcast(Inquire_field)

  END FUNCTION Inquire_Field
  
 
  SUBROUTINE Get_Field_r1(field_name,field,found)
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)    :: Field_name
    REAL,INTENT(INOUT)               :: Field(:)
    LOGICAL,INTENT(OUT),OPTIONAL   :: found 

    IF (PRESENT(found)) THEN
      CALL Get_field_rgen(field_name,field,1,found)
    ELSE
      CALL Get_field_rgen(field_name,field,1)
    ENDIF
      
  END SUBROUTINE Get_Field_r1
  
  SUBROUTINE Get_Field_r2(field_name,field,found)
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)    :: Field_name
    REAL,INTENT(INOUT)               :: Field(:,:)
    LOGICAL,INTENT(OUT),OPTIONAL   :: found 

    IF (PRESENT(found)) THEN
      CALL Get_field_rgen(field_name,field,size(field,2),found)
    ELSE
      CALL Get_field_rgen(field_name,field,size(field,2))
    ENDIF

      
  END SUBROUTINE Get_Field_r2
  
  SUBROUTINE Get_Field_r3(field_name,field,found)
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)    :: Field_name
    REAL,INTENT(INOUT)               :: Field(:,:,:)
    LOGICAL,INTENT(OUT),OPTIONAL   :: found 

    IF (PRESENT(found)) THEN
      CALL Get_field_rgen(field_name,field,size(field,2)*size(field,3),found)
    ELSE
      CALL Get_field_rgen(field_name,field,size(field,2)*size(field,3))
    ENDIF
      
  END SUBROUTINE Get_Field_r3
  
  SUBROUTINE Get_field_rgen(field_name,field,field_size,found)
  USE netcdf
  USE dimphy
  USE mod_grid_phy_lmdz
  USE mod_phys_lmdz_para
  IMPLICIT NONE
    CHARACTER(LEN=*) :: Field_name
    INTEGER          :: field_size
    REAL             :: field(klon,field_size)
    LOGICAL,OPTIONAL :: found
    
    REAL    :: field_glo(klon_glo,field_size)
    LOGICAL :: tmp_found
    INTEGER :: varid
    INTEGER :: ierr
    
    IF (is_mpi_root .AND. is_omp_root) THEN
  
      ierr=NF90_INQ_VARID(nid_start,Field_name,varid)
      
      IF (ierr==NF90_NOERR) THEN
        CALL body(field_glo)
        tmp_found=.TRUE.
      ELSE
        tmp_found=.FALSE.
      ENDIF
    
    ENDIF
    
    CALL bcast(tmp_found)

    IF (tmp_found) THEN
      CALL scatter(field_glo,field)
    ENDIF
    
    IF (PRESENT(found)) THEN
      found=tmp_found
    ELSE
      IF (.NOT. tmp_found) THEN
        PRINT*, 'phyetat0: Le champ <'//field_name//'> est absent'
        CALL abort
      ENDIF
    ENDIF
 
    
    CONTAINS
     
     SUBROUTINE body(field_glo)
       REAL :: field_glo(klon_glo*field_size)
         ierr=NF90_GET_VAR(nid_start,varid,field_glo)
         IF (ierr/=NF90_NOERR) THEN
           ! La variable exist dans le fichier mais la lecture a echouee. 
           PRINT*, 'phyetat0: Lecture echouee pour <'//field_name//'>'

           IF (field_name=='CLWCON' .OR. field_name=='RNEBCON' .OR. field_name=='RATQS') THEN
              ! Essaye de lire le variable sur surface uniqument, comme fait avant
              field_glo(:)=0.
              ierr=NF90_GET_VAR(nid_start,varid,field_glo(1:klon_glo))
              IF (ierr/=NF90_NOERR) THEN
                 PRINT*, 'phyetat0: Lecture echouee aussi en 2D pour <'//field_name//'>'
                 CALL abort
              ELSE
                 PRINT*, 'phyetat0: La variable <'//field_name//'> lu sur surface seulement'!, selon ancien format, le reste mis a zero'
              END IF
           ELSE
              CALL abort
           ENDIF
         ENDIF

     END SUBROUTINE body

  END SUBROUTINE Get_field_rgen
  

  SUBROUTINE get_var_r0(var_name,var,found)
  IMPLICIT NONE  
    CHARACTER(LEN=*),INTENT(IN)  :: var_name
    REAL,INTENT(INOUT)             :: var
    LOGICAL,OPTIONAL,INTENT(OUT) :: found

    REAL                         :: varout(1)
    
    IF (PRESENT(found)) THEN
      CALL Get_var_rgen(var_name,varout,size(varout),found)
    ELSE
      CALL Get_var_rgen(var_name,varout,size(varout))
    ENDIF
    var=varout(1)
 
  END SUBROUTINE get_var_r0

  SUBROUTINE get_var_r1(var_name,var,found)
  IMPLICIT NONE  
    CHARACTER(LEN=*),INTENT(IN)  :: var_name
    REAL,INTENT(INOUT)             :: var(:)
    LOGICAL,OPTIONAL,INTENT(OUT) :: found
    
    IF (PRESENT(found)) THEN
      CALL Get_var_rgen(var_name,var,size(var),found)
    ELSE
      CALL Get_var_rgen(var_name,var,size(var))
    ENDIF
  
  END SUBROUTINE get_var_r1

  SUBROUTINE get_var_r2(var_name,var,found)
  IMPLICIT NONE  
    CHARACTER(LEN=*),INTENT(IN)  :: var_name
    REAL,INTENT(OUT)             :: var(:,:)
    LOGICAL,OPTIONAL,INTENT(OUT) :: found
    
    IF (PRESENT(found)) THEN
      CALL Get_var_rgen(var_name,var,size(var),found)
    ELSE
      CALL Get_var_rgen(var_name,var,size(var))
    ENDIF
  
  END SUBROUTINE get_var_r2

  SUBROUTINE get_var_r3(var_name,var,found)
  IMPLICIT NONE  
    CHARACTER(LEN=*),INTENT(IN)  :: var_name
    REAL,INTENT(INOUT)             :: var(:,:,:)
    LOGICAL,OPTIONAL,INTENT(OUT) :: found
    
    IF (PRESENT(found)) THEN
      CALL Get_var_rgen(var_name,var,size(var),found)
    ELSE
      CALL Get_var_rgen(var_name,var,size(var))
    ENDIF
  
  END SUBROUTINE get_var_r3

  SUBROUTINE Get_var_rgen(var_name,var,var_size,found)
  USE netcdf
  USE dimphy
  USE mod_grid_phy_lmdz
  USE mod_phys_lmdz_para
  IMPLICIT NONE
    CHARACTER(LEN=*) :: var_name
    INTEGER          :: var_size
    REAL             :: var(var_size)
    LOGICAL,OPTIONAL :: found
    
    LOGICAL :: tmp_found
    INTEGER :: varid
    INTEGER :: ierr
    
    IF (is_mpi_root .AND. is_omp_root) THEN
  
      ierr=NF90_INQ_VARID(nid_start,var_name,varid)
      
      IF (ierr==NF90_NOERR) THEN
        ierr=NF90_GET_VAR(nid_start,varid,var)
        IF (ierr/=NF90_NOERR) THEN
          PRINT*, 'phyetat0: Lecture echouee pour <'//var_name//'>'
          CALL abort
        ENDIF
        tmp_found=.TRUE.
      ELSE
        tmp_found=.FALSE.
      ENDIF
    
    ENDIF
    
    CALL bcast(tmp_found)

    IF (tmp_found) THEN
      CALL bcast(var)
    ENDIF
    
    IF (PRESENT(found)) THEN
      found=tmp_found
    ELSE
      IF (.NOT. tmp_found) THEN
        PRINT*, 'phyetat0: La variable champ <'//var_name//'> est absente'
        CALL abort
      ENDIF
    ENDIF

  END SUBROUTINE Get_var_rgen


  SUBROUTINE open_restartphy(filename)
  USE netcdf
  USE mod_phys_lmdz_para
  USE mod_grid_phy_lmdz
  USE dimphy
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: filename
    INTEGER                     :: ierr
    
    IF (is_mpi_root .AND. is_omp_root) THEN
      ierr = NF90_CREATE(filename, NF90_CLOBBER, nid_restart)
      IF (ierr/=NF90_NOERR) THEN
        write(6,*)' Pb d''ouverture du fichier '//filename
        write(6,*)' ierr = ', ierr
        CALL ABORT
      ENDIF

      ierr = NF90_PUT_ATT (nid_restart, NF90_GLOBAL, "title","Fichier redemarrage physique")

      ierr = NF90_DEF_DIM (nid_restart, "index", length, idim1)
      ierr = NF90_DEF_DIM (nid_restart, "points_physiques", klon_glo, idim2)
      ierr = NF90_DEF_DIM (nid_restart, "horizon_vertical", klon_glo*klev, idim3)
      ierr = NF90_DEF_DIM (nid_restart, "horizon_klevp1", klon_glo*klevp1, idim4)

      ierr = NF90_ENDDEF(nid_restart)
    ENDIF

  END SUBROUTINE open_restartphy
  
  SUBROUTINE close_restartphy
  USE netcdf
  USE mod_phys_lmdz_para
  IMPLICIT NONE
    INTEGER          :: ierr

    IF (is_mpi_root .AND. is_omp_root) THEN
      ierr = NF90_CLOSE (nid_restart)
    ENDIF
 
  END SUBROUTINE close_restartphy

  
  SUBROUTINE put_field_r1(field_name,title,field)
  IMPLICIT NONE
  CHARACTER(LEN=*),INTENT(IN)    :: field_name
  CHARACTER(LEN=*),INTENT(IN)    :: title
  REAL,INTENT(IN)                :: field(:)
  
    CALL put_field_rgen(field_name,title,field,1)
  
  END SUBROUTINE put_field_r1

  SUBROUTINE put_field_r2(field_name,title,field)
  IMPLICIT NONE
  CHARACTER(LEN=*),INTENT(IN)    :: field_name
  CHARACTER(LEN=*),INTENT(IN)    :: title
  REAL,INTENT(IN)                :: field(:,:)
  
    CALL put_field_rgen(field_name,title,field,size(field,2))
  
  END SUBROUTINE put_field_r2

  SUBROUTINE put_field_r3(field_name,title,field)
  IMPLICIT NONE
  CHARACTER(LEN=*),INTENT(IN)    :: field_name
  CHARACTER(LEN=*),INTENT(IN)    :: title
  REAL,INTENT(IN)                :: field(:,:,:)
  
    CALL put_field_rgen(field_name,title,field,size(field,2)*size(field,3))
  
  END SUBROUTINE put_field_r3
  
  SUBROUTINE put_field_rgen(field_name,title,field,field_size)
  USE netcdf
  USE dimphy
  USE mod_grid_phy_lmdz
  USE mod_phys_lmdz_para
  IMPLICIT NONE
  CHARACTER(LEN=*),INTENT(IN)    :: field_name
  CHARACTER(LEN=*),INTENT(IN)    :: title
  INTEGER,INTENT(IN)             :: field_size
  REAL,INTENT(IN)                :: field(klon,field_size)
  
  REAL                           :: field_glo(klon_glo,field_size)
  INTEGER                        :: ierr
  INTEGER                        :: nvarid
  INTEGER                        :: idim
   
   
    CALL gather(field,field_glo)
    
    IF (is_mpi_root .AND. is_omp_root) THEN

      IF (field_size==1) THEN
        idim=idim2
      ELSE IF (field_size==klev) THEN
        idim=idim3
      ELSE IF (field_size==klevp1) THEN
        idim=idim4
      ELSE
        PRINT *, "erreur phyredem : probleme de dimension"
        CALL ABORT
      ENDIF
         
      ierr = NF90_REDEF (nid_restart)
#ifdef NC_DOUBLE
      ierr = NF90_DEF_VAR (nid_restart, field_name, NF90_DOUBLE,(/ idim /),nvarid)
#else
      ierr = NF90_DEF_VAR (nid_restart, field_name, NF90_FLOAT,(/ idim /),nvarid)
#endif
      IF (LEN_TRIM(title) > 0) ierr = NF90_PUT_ATT (nid_restart,nvarid,"title", title)
      ierr = NF90_ENDDEF(nid_restart)
      ierr = NF90_PUT_VAR(nid_restart,nvarid,RESHAPE(field_glo,(/klon_glo*field_size/)))
    ENDIF
    
   END SUBROUTINE put_field_rgen  
  
   SUBROUTINE put_var_r0(var_name,title,var)
   IMPLICIT NONE
     CHARACTER(LEN=*),INTENT(IN) :: var_name
     CHARACTER(LEN=*),INTENT(IN) :: title
     REAL,INTENT(IN)             :: var
     REAL                        :: varin(1)
     
     varin(1)=var
     
     CALL put_var_rgen(var_name,title,varin,size(varin))

  END SUBROUTINE put_var_r0


   SUBROUTINE put_var_r1(var_name,title,var)
   IMPLICIT NONE
     CHARACTER(LEN=*),INTENT(IN) :: var_name
     CHARACTER(LEN=*),INTENT(IN) :: title
     REAL,INTENT(IN)             :: var(:)
     
     CALL put_var_rgen(var_name,title,var,size(var))

  END SUBROUTINE put_var_r1
 
  SUBROUTINE put_var_r2(var_name,title,var)
   IMPLICIT NONE
     CHARACTER(LEN=*),INTENT(IN) :: var_name
     CHARACTER(LEN=*),INTENT(IN) :: title
     REAL,INTENT(IN)             :: var(:,:)
     
     CALL put_var_rgen(var_name,title,var,size(var))

  END SUBROUTINE put_var_r2     
  
  SUBROUTINE put_var_r3(var_name,title,var)
   IMPLICIT NONE
     CHARACTER(LEN=*),INTENT(IN) :: var_name
     CHARACTER(LEN=*),INTENT(IN) :: title
     REAL,INTENT(IN)             :: var(:,:,:)
     
     CALL put_var_rgen(var_name,title,var,size(var))

  END SUBROUTINE put_var_r3

  SUBROUTINE put_var_rgen(var_name,title,var,var_size)
  USE netcdf
  USE dimphy
  USE mod_phys_lmdz_para
  IMPLICIT NONE
     CHARACTER(LEN=*),INTENT(IN) :: var_name
     CHARACTER(LEN=*),INTENT(IN) :: title
     INTEGER,INTENT(IN)          :: var_size
     REAL,INTENT(IN)             :: var(var_size)
     
     INTEGER :: ierr
     INTEGER :: nvarid
         
    IF (is_mpi_root .AND. is_omp_root) THEN

      IF (var_size/=length) THEN
        PRINT *, "erreur phyredem : probleme de dimension"
        CALL abort
      ENDIF
      
      ierr = NF90_REDEF (nid_restart)

#ifdef NC_DOUBLE
      ierr = NF90_DEF_VAR (nid_restart, var_name, NF90_DOUBLE,(/ idim1 /),nvarid)
#else
      ierr = NF90_DEF_VAR (nid_restart, var_name, NF90_FLOAT,(/ idim1 /),nvarid)
#endif
      IF (LEN_TRIM(title)>0) ierr = NF90_PUT_ATT (nid_restart,nvarid,"title", title)
      ierr = NF90_ENDDEF(nid_restart)
     
      ierr = NF90_PUT_VAR(nid_restart,nvarid,var)

    ENDIF
    
  END SUBROUTINE put_var_rgen     
    
END MODULE iostart
