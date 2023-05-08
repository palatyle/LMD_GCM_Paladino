MODULE distrib_icosa_lmdz_mod
  
  TYPE t_distrib_physic
    INTEGER,POINTER :: index(:)    ! list of index used for thread entering in lmdz physic
    INTEGER         :: nindex      ! number of index used
    INTEGER         :: domain_ind  ! index of the related domain
  END TYPE t_distrib_physic
  
  INTEGER, SAVE :: ndomain_distrib    ! number of domain needed for thread data
!$OMP THREADPRIVATE(ndomain_distrib)

  TYPE(t_distrib_physic),ALLOCATABLE, SAVE :: distrib_physic(:) 
!$OMP THREADPRIVATE(distrib_physic)  


  INTERFACE transfer_icosa_to_lmdz
    MODULE PROCEDURE transfer_icosa_to_lmdz1d,transfer_icosa_to_lmdz2d,transfer_icosa_to_lmdz3d
  END INTERFACE transfer_icosa_to_lmdz

  INTERFACE transfer_lmdz_to_icosa
    MODULE PROCEDURE transfer_lmdz1d_to_icosa,transfer_lmdz2d_to_icosa,transfer_lmdz3d_to_icosa
  END INTERFACE transfer_lmdz_to_icosa
      
CONTAINS


  SUBROUTINE init_distrib_icosa_lmdz
  USE mod_phys_lmdz_omp_data, ONLY: klon_omp_begin, klon_omp_end
  USE domain_mod
  USE dimensions
  IMPLICIT NONE
    INTEGER :: pos,pos_tmp,nindex
    INTEGER :: ind, i,j,ij
    
    ALLOCATE(distrib_physic(ndomain))
   
    ndomain_distrib=0 
    pos=0
    DO ind=1,ndomain
      CALL swap_dimensions(ind)

! first guess to determine number of indices for this domain
      pos_tmp=pos
      nindex=0
      DO j=jj_begin,jj_end
        DO i=ii_begin,ii_end
          IF (domain(ind)%own(i,j)) THEN 
            pos_tmp=pos_tmp+1
            IF (pos_tmp >= klon_omp_begin .AND. pos_tmp <= klon_omp_end) nindex=nindex+1
          ENDIF
        ENDDO
      ENDDO

! fill the index array

      IF (nindex>0) THEN
        ndomain_distrib=ndomain_distrib+1
        ALLOCATE(distrib_physic(ndomain_distrib)%index(nindex))
        distrib_physic(ndomain_distrib)%nindex=nindex
        distrib_physic(ndomain_distrib)%domain_ind=ind
        
        nindex=0
        DO j=jj_begin,jj_end
          DO i=ii_begin,ii_end
            ij=(j-1)*iim+i
            IF (domain(ind)%own(i,j)) THEN
              pos=pos+1
              IF (pos >= klon_omp_begin .AND. pos <= klon_omp_end) THEN
                nindex=nindex+1
                distrib_physic(ndomain_distrib)%index(nindex)=ij
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ELSE
        pos=pos_tmp
      ENDIF
    
    ENDDO
       
  END SUBROUTINE init_distrib_icosa_lmdz
                  
  SUBROUTINE transfer_icosa_to_lmdz1d(f_field_icosa, field_lmdz)
  USE field_mod
  IMPLICIT NONE
    TYPE(t_field),POINTER :: f_field_icosa(:)
    REAL(rstd)         :: field_lmdz(:)
    REAL(rstd),POINTER :: field_icosa(:)
    INTEGER         :: pos, nindex,ind,i
    INTEGER,POINTER :: index(:)
    
!$OMP BARRIER
    pos=0
    DO ind=1,ndomain_distrib
      field_icosa=f_field_icosa(distrib_physic(ind)%domain_ind)
      index=>distrib_physic(ind)%index
      nindex=distrib_physic(ind)%nindex
      DO i=1,nindex
        pos=pos+1
        field_lmdz(pos)=field_icosa(index(i))
      ENDDO
    ENDDO
    
  END SUBROUTINE  transfer_icosa_to_lmdz1d
  
  SUBROUTINE transfer_icosa_to_lmdz2d(f_field_icosa, field_lmdz)
  USE field_mod
  IMPLICIT NONE
    TYPE(t_field),POINTER :: f_field_icosa(:)
    REAL(rstd)         :: field_lmdz(:,:)

    REAL(rstd),POINTER :: field_icosa(:,:)
    INTEGER         :: pos, nindex,ind,i
    INTEGER,POINTER :: index(:)
    INTEGER :: l
    
!$OMP BARRIER
    DO l=1,size(field_lmdz,2)  
      pos=0
      DO ind=1,ndomain_distrib
        field_icosa=f_field_icosa(distrib_physic(ind)%domain_ind)
        index=>distrib_physic(ind)%index
        nindex=distrib_physic(ind)%nindex
        DO i=1,nindex
          pos=pos+1
          field_lmdz(pos,l)=field_icosa(index(i),l)
        ENDDO
      ENDDO
    ENDDO
    
  END SUBROUTINE  transfer_icosa_to_lmdz2d


    
  SUBROUTINE transfer_icosa_to_lmdz3d(f_field_icosa, field_lmdz)
  USE field_mod
  IMPLICIT NONE
    TYPE(t_field),POINTER :: f_field_icosa(:)
    REAL(rstd)         :: field_lmdz(:,:,:)
    REAL(rstd),POINTER :: field_icosa(:,:,:)
    INTEGER         :: pos, nindex,ind,i
    INTEGER,POINTER :: index(:)
    INTEGER :: l,q

!$OMP BARRIER
    DO q=1,size(field_lmdz,3)  
      DO l=1,size(field_lmdz,2)  
        pos=0
        DO ind=1,ndomain_distrib
          field_icosa=f_field_icosa(distrib_physic(ind)%domain_ind)
          index=>distrib_physic(ind)%index
          nindex=distrib_physic(ind)%nindex
          DO i=1,nindex
            pos=pos+1
            field_lmdz(pos,l,q)=field_icosa(index(i),l,q)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    
  END SUBROUTINE  transfer_icosa_to_lmdz3d
      
  SUBROUTINE transfer_lmdz1d_to_icosa(field_lmdz,f_field_icosa)
  USE field_mod
  IMPLICIT NONE
    REAL(rstd)         :: field_lmdz(:)
    TYPE(t_field),POINTER :: f_field_icosa(:)
    REAL(rstd),POINTER :: field_icosa(:)
    INTEGER         :: pos, nindex,ind,i
    INTEGER,POINTER :: index(:)
    
!$OMP BARRIER
    pos=0
    DO ind=1,ndomain_distrib
      field_icosa=f_field_icosa(distrib_physic(ind)%domain_ind)
      index=>distrib_physic(ind)%index
      nindex=distrib_physic(ind)%nindex
      DO i=1,nindex
        pos=pos+1
        field_icosa(index(i))=field_lmdz(pos)
      ENDDO
    ENDDO
  END SUBROUTINE  transfer_lmdz1d_to_icosa

  SUBROUTINE transfer_lmdz2d_to_icosa(field_lmdz,f_field_icosa)
  USE field_mod
  IMPLICIT NONE
    REAL(rstd)         :: field_lmdz(:,:)
    TYPE(t_field),POINTER :: f_field_icosa(:)
    REAL(rstd),POINTER :: field_icosa(:,:)
    INTEGER         :: pos, nindex,ind,i
    INTEGER,POINTER :: index(:)
    INTEGER :: l
    
!$OMP BARRIER
    DO l=1,size(field_lmdz,2)  
      pos=0
      DO ind=1,ndomain_distrib
        field_icosa=f_field_icosa(distrib_physic(ind)%domain_ind)
        index=>distrib_physic(ind)%index
        nindex=distrib_physic(ind)%nindex
        DO i=1,nindex
          pos=pos+1
          field_icosa(index(i),l)=field_lmdz(pos,l)
        ENDDO
      ENDDO
    ENDDO
  
  END SUBROUTINE transfer_lmdz2d_to_icosa    


  SUBROUTINE transfer_lmdz3d_to_icosa(field_lmdz,f_field_icosa)
  USE field_mod
  IMPLICIT NONE
    REAL(rstd)         :: field_lmdz(:,:,:)
    TYPE(t_field),POINTER :: f_field_icosa(:)
    REAL(rstd),POINTER :: field_icosa(:,:,:)
    INTEGER         :: pos, nindex,ind,i
    INTEGER,POINTER :: index(:)
    INTEGER :: l,q
    
!$OMP BARRIER
    DO q=1,size(field_lmdz,3)  
      DO l=1,size(field_lmdz,2)  
        pos=0
        DO ind=1,ndomain_distrib
          field_icosa=f_field_icosa(distrib_physic(ind)%domain_ind)
          index=>distrib_physic(ind)%index
          nindex=distrib_physic(ind)%nindex
          DO i=1,nindex
            pos=pos+1
            field_icosa(index(i),l,q)=field_lmdz(pos,l,q)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  
  END SUBROUTINE transfer_lmdz3d_to_icosa    

END MODULE distrib_icosa_lmdz_mod
