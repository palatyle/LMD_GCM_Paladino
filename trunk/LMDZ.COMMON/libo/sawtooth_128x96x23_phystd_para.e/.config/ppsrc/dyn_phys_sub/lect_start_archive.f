










      SUBROUTINE lect_start_archive(ngrid,nlayer,
     &     date,tsurf,tsoil,emis,q2,
     &     t,ucov,vcov,ps,h,phisold_newgrid,
     &     q,qsurf,surfith,nid,
     &     rnat,pctsrf_sic,tslab,tsea_ice,sea_ice)

!      USE surfdat_h
      USE comsoil_h, ONLY: nsoilmx, layer, mlayer, volcapa, inertiedat
      USE tracer_h, ONLY: igcm_co2_ice
      USE infotrac, ONLY: tname, nqtot
      USE slab_ice_h, ONLY: noceanmx
!      USE control_mod
! to use  'getin'
      USE callkeys_mod, ONLY: ok_slab_ocean
      USE comvert_mod, ONLY: ap,bp,aps,bps,preff
      USE comconst_mod, ONLY: kappa,g,pi

c=======================================================================
c
c
c   Auteur:    05/1997 , 12/2003 : coord hybride  FF
c   ------
c
c
c   Objet:     Lecture des variables d'un fichier "start_archive"
c              Plus besoin de régler ancienne valeurs grace
c              a l'allocation dynamique de memoire (Yann Wanherdrick)
c
c
c
c=======================================================================

      implicit none

!-----------------------------------------------------------------------
!   INCLUDE 'dimensions.h'
!
!   dimensions.h contient les dimensions du modele
!   ndm est tel que iim=2**ndm
!-----------------------------------------------------------------------

      INTEGER iim,jjm,llm,ndm

      PARAMETER (iim= 128,jjm=96,llm=23,ndm=1)

!-----------------------------------------------------------------------
!#include "dimphys.h"
!#include "planete.h"
!
! $Header$
!
!
!  ATTENTION!!!!: ce fichier include est compatible format fixe/format libre
!                 veillez  n'utiliser que des ! pour les commentaires
!                 et  bien positionner les & des lignes de continuation
!                 (les placer en colonne 6 et en colonne 73)
!
!
!-----------------------------------------------------------------------
!   INCLUDE 'paramet.h'

      INTEGER  iip1,iip2,iip3,jjp1,llmp1,llmp2,llmm1
      INTEGER  kftd,ip1jm,ip1jmp1,ip1jmi1,ijp1llm
      INTEGER  ijmllm,mvar
      INTEGER jcfil,jcfllm

      PARAMETER( iip1= iim+1,iip2=iim+2,iip3=iim+3                       &
     &    ,jjp1=jjm+1-1/jjm)
      PARAMETER( llmp1 = llm+1,  llmp2 = llm+2, llmm1 = llm-1 )
      PARAMETER( kftd  = iim/2 -ndm )
      PARAMETER( ip1jm  = iip1*jjm,  ip1jmp1= iip1*jjp1 )
      PARAMETER( ip1jmi1= ip1jm - iip1 )
      PARAMETER( ijp1llm= ip1jmp1 * llm, ijmllm= ip1jm * llm )
      PARAMETER( mvar= ip1jmp1*( 2*llm+1) + ijmllm )
      PARAMETER( jcfil=jjm/2+5, jcfllm=jcfil*llm )

!-----------------------------------------------------------------------
!
! $Header$
!
!CDK comgeom2
      COMMON/comgeom/                                                   &
     & cu(iip1,jjp1),cv(iip1,jjm),unscu2(iip1,jjp1),unscv2(iip1,jjm)  , &
     & aire(iip1,jjp1),airesurg(iip1,jjp1),aireu(iip1,jjp1)           , &
     & airev(iip1,jjm),unsaire(iip1,jjp1),apoln,apols                 , &
     & unsairez(iip1,jjm),airuscv2(iip1,jjm),airvscu2(iip1,jjm)       , &
     & aireij1(iip1,jjp1),aireij2(iip1,jjp1),aireij3(iip1,jjp1)       , &
     & aireij4(iip1,jjp1),alpha1(iip1,jjp1),alpha2(iip1,jjp1)         , &
     & alpha3(iip1,jjp1),alpha4(iip1,jjp1),alpha1p2(iip1,jjp1)        , &
     & alpha1p4(iip1,jjp1),alpha2p3(iip1,jjp1),alpha3p4(iip1,jjp1)    , &
     & fext(iip1,jjm),constang(iip1,jjp1), rlatu(jjp1),rlatv(jjm),      &
     & rlonu(iip1),rlonv(iip1),cuvsurcv(iip1,jjm),cvsurcuv(iip1,jjm)  , &
     & cvusurcu(iip1,jjp1),cusurcvu(iip1,jjp1)                        , &
     & cuvscvgam1(iip1,jjm),cuvscvgam2(iip1,jjm),cvuscugam1(iip1,jjp1), &
     & cvuscugam2(iip1,jjp1),cvscuvgam(iip1,jjm),cuscvugam(iip1,jjp1) , &
     & unsapolnga1,unsapolnga2,unsapolsga1,unsapolsga2                , &
     & unsair_gam1(iip1,jjp1),unsair_gam2(iip1,jjp1)                  , &
     & unsairz_gam(iip1,jjm),aivscu2gam(iip1,jjm),aiuscv2gam(iip1,jjm)  &
     & , xprimu(iip1),xprimv(iip1)


      REAL                                                               &
     & cu,cv,unscu2,unscv2,aire,airesurg,aireu,airev,apoln,apols,unsaire &
     & ,unsairez,airuscv2,airvscu2,aireij1,aireij2,aireij3,aireij4     , &
     & alpha1,alpha2,alpha3,alpha4,alpha1p2,alpha1p4,alpha2p3,alpha3p4 , &
     & fext,constang,rlatu,rlatv,rlonu,rlonv,cuvscvgam1,cuvscvgam2     , &
     & cvuscugam1,cvuscugam2,cvscuvgam,cuscvugam,unsapolnga1           , &
     & unsapolnga2,unsapolsga1,unsapolsga2,unsair_gam1,unsair_gam2     , &
     & unsairz_gam,aivscu2gam,aiuscv2gam,cuvsurcv,cvsurcuv,cvusurcu    , &
     & cusurcvu,xprimu,xprimv
!#include "control.h"
!     NetCDF-3.
!
! netcdf version 3 fortran interface:
!

!
! external netcdf data types:
!
      integer nf_byte
      integer nf_int1
      integer nf_char
      integer nf_short
      integer nf_int2
      integer nf_int
      integer nf_float
      integer nf_real
      integer nf_double
      integer nf_ubyte
      integer nf_ushort
      integer nf_uint
      integer nf_int64
      integer nf_uint64

      parameter (nf_byte = 1)
      parameter (nf_int1 = nf_byte)
      parameter (nf_char = 2)
      parameter (nf_short = 3)
      parameter (nf_int2 = nf_short)
      parameter (nf_int = 4)
      parameter (nf_float = 5)
      parameter (nf_real = nf_float)
      parameter (nf_double = 6)
      parameter (nf_ubyte = 7)
      parameter (nf_ushort = 8)
      parameter (nf_uint = 9)
      parameter (nf_int64 = 10)
      parameter (nf_uint64 = 11)

!
! default fill values:
!
      integer           nf_fill_byte
      integer           nf_fill_int1
      integer           nf_fill_char
      integer           nf_fill_short
      integer           nf_fill_int2
      integer           nf_fill_int
      real              nf_fill_float
      real              nf_fill_real
      doubleprecision   nf_fill_double

      parameter (nf_fill_byte = -127)
      parameter (nf_fill_int1 = nf_fill_byte)
      parameter (nf_fill_char = 0)
      parameter (nf_fill_short = -32767)
      parameter (nf_fill_int2 = nf_fill_short)
      parameter (nf_fill_int = -2147483647)
      parameter (nf_fill_float = 9.9692099683868690e+36)
      parameter (nf_fill_real = nf_fill_float)
      parameter (nf_fill_double = 9.9692099683868690d+36)

!
! mode flags for opening and creating a netcdf dataset:
!
      integer nf_nowrite
      integer nf_write
      integer nf_clobber
      integer nf_noclobber
      integer nf_fill
      integer nf_nofill
      integer nf_lock
      integer nf_share
      integer nf_64bit_offset
      integer nf_64bit_data
      integer nf_cdf5
      integer nf_sizehint_default
      integer nf_align_chunk
      integer nf_format_classic
      integer nf_format_64bit
      integer nf_format_64bit_offset
      integer nf_format_64bit_data
      integer nf_format_cdf5
      integer nf_diskless
      integer nf_mmap

      parameter (nf_nowrite = 0)
      parameter (nf_write = 1)
      parameter (nf_clobber = 0)
      parameter (nf_noclobber = 4)
      parameter (nf_fill = 0)
      parameter (nf_nofill = 256)
      parameter (nf_lock = 1024)
      parameter (nf_share = 2048)
      parameter (nf_64bit_offset = 512)
      parameter (nf_64bit_data = 32)
      parameter (nf_cdf5 = nf_64bit_data)
      parameter (nf_sizehint_default = 0)
      parameter (nf_align_chunk = -1)
      parameter (nf_format_classic = 1)
      parameter (nf_format_64bit = 2)
      parameter (nf_format_64bit_offset = nf_format_64bit)
      parameter (nf_format_64bit_data = 5)
      parameter (nf_format_cdf5 = nf_format_64bit_data)
      parameter (nf_diskless = 8)
      parameter (nf_mmap = 16)

!
! size argument for defining an unlimited dimension:
!
      integer nf_unlimited
      parameter (nf_unlimited = 0)

!
! global attribute id:
!
      integer nf_global
      parameter (nf_global = 0)

!
! implementation limits:
!
      integer nf_max_dims
      integer nf_max_attrs
      integer nf_max_vars
      integer nf_max_name
      integer nf_max_var_dims

      parameter (nf_max_dims = 1024)
      parameter (nf_max_attrs = 8192)
      parameter (nf_max_vars = 8192)
      parameter (nf_max_name = 256)
      parameter (nf_max_var_dims = nf_max_dims)

!
! error codes:
!
      integer nf_noerr
      integer nf_ebadid
      integer nf_eexist
      integer nf_einval
      integer nf_eperm
      integer nf_enotindefine
      integer nf_eindefine
      integer nf_einvalcoords
      integer nf_emaxdims
      integer nf_enameinuse
      integer nf_enotatt
      integer nf_emaxatts
      integer nf_ebadtype
      integer nf_ebaddim
      integer nf_eunlimpos
      integer nf_emaxvars
      integer nf_enotvar
      integer nf_eglobal
      integer nf_enotnc
      integer nf_ests
      integer nf_emaxname
      integer nf_eunlimit
      integer nf_enorecvars
      integer nf_echar
      integer nf_eedge
      integer nf_estride
      integer nf_ebadname
      integer nf_erange
      integer nf_enomem
      integer nf_evarsize
      integer nf_edimsize
      integer nf_etrunc

      parameter (nf_noerr = 0)
      parameter (nf_ebadid = -33)
      parameter (nf_eexist = -35)
      parameter (nf_einval = -36)
      parameter (nf_eperm = -37)
      parameter (nf_enotindefine = -38)
      parameter (nf_eindefine = -39)
      parameter (nf_einvalcoords = -40)
      parameter (nf_emaxdims = -41)
      parameter (nf_enameinuse = -42)
      parameter (nf_enotatt = -43)
      parameter (nf_emaxatts = -44)
      parameter (nf_ebadtype = -45)
      parameter (nf_ebaddim = -46)
      parameter (nf_eunlimpos = -47)
      parameter (nf_emaxvars = -48)
      parameter (nf_enotvar = -49)
      parameter (nf_eglobal = -50)
      parameter (nf_enotnc = -51)
      parameter (nf_ests = -52)
      parameter (nf_emaxname = -53)
      parameter (nf_eunlimit = -54)
      parameter (nf_enorecvars = -55)
      parameter (nf_echar = -56)
      parameter (nf_eedge = -57)
      parameter (nf_estride = -58)
      parameter (nf_ebadname = -59)
      parameter (nf_erange = -60)
      parameter (nf_enomem = -61)
      parameter (nf_evarsize = -62)
      parameter (nf_edimsize = -63)
      parameter (nf_etrunc = -64)
!
! error handling modes:
!
      integer  nf_fatal
      integer nf_verbose

      parameter (nf_fatal = 1)
      parameter (nf_verbose = 2)

!
! miscellaneous routines:
!
      character*80   nf_inq_libvers
      external       nf_inq_libvers

      character*80   nf_strerror
!                         (integer             ncerr)
      external       nf_strerror

      logical        nf_issyserr
!                         (integer             ncerr)
      external       nf_issyserr

!
! control routines:
!
      integer         nf_inq_base_pe
!                         (integer             ncid,
!                          integer             pe)
      external        nf_inq_base_pe

      integer         nf_set_base_pe
!                         (integer             ncid,
!                          integer             pe)
      external        nf_set_base_pe

      integer         nf_create
!                         (character*(*)       path,
!                          integer             cmode,
!                          integer             ncid)
      external        nf_create

      integer         nf__create
!                         (character*(*)       path,
!                          integer             cmode,
!                          integer             initialsz,
!                          integer             chunksizehint,
!                          integer             ncid)
      external        nf__create

      integer         nf__create_mp
!                         (character*(*)       path,
!                          integer             cmode,
!                          integer             initialsz,
!                          integer             basepe,
!                          integer             chunksizehint,
!                          integer             ncid)
      external        nf__create_mp

      integer         nf_open
!                         (character*(*)       path,
!                          integer             mode,
!                          integer             ncid)
      external        nf_open

      integer         nf__open
!                         (character*(*)       path,
!                          integer             mode,
!                          integer             chunksizehint,
!                          integer             ncid)
      external        nf__open

      integer         nf__open_mp
!                         (character*(*)       path,
!                          integer             mode,
!                          integer             basepe,
!                          integer             chunksizehint,
!                          integer             ncid)
      external        nf__open_mp

      integer         nf_set_fill
!                         (integer             ncid,
!                          integer             fillmode,
!                          integer             old_mode)
      external        nf_set_fill

      integer         nf_set_default_format
!                          (integer             format,
!                          integer             old_format)
      external        nf_set_default_format

      integer         nf_redef
!                         (integer             ncid)
      external        nf_redef

      integer         nf_enddef
!                         (integer             ncid)
      external        nf_enddef

      integer         nf__enddef
!                         (integer             ncid,
!                          integer             h_minfree,
!                          integer             v_align,
!                          integer             v_minfree,
!                          integer             r_align)
      external        nf__enddef

      integer         nf_sync
!                         (integer             ncid)
      external        nf_sync

      integer         nf_abort
!                         (integer             ncid)
      external        nf_abort

      integer         nf_close
!                         (integer             ncid)
      external        nf_close

      integer         nf_delete
!                         (character*(*)       ncid)
      external        nf_delete

!
! general inquiry routines:
!

      integer         nf_inq
!                         (integer             ncid,
!                          integer             ndims,
!                          integer             nvars,
!                          integer             ngatts,
!                          integer             unlimdimid)
      external        nf_inq

! new inquire path

      integer nf_inq_path
      external nf_inq_path

      integer         nf_inq_ndims
!                         (integer             ncid,
!                          integer             ndims)
      external        nf_inq_ndims

      integer         nf_inq_nvars
!                         (integer             ncid,
!                          integer             nvars)
      external        nf_inq_nvars

      integer         nf_inq_natts
!                         (integer             ncid,
!                          integer             ngatts)
      external        nf_inq_natts

      integer         nf_inq_unlimdim
!                         (integer             ncid,
!                          integer             unlimdimid)
      external        nf_inq_unlimdim

      integer         nf_inq_format
!                         (integer             ncid,
!                          integer             format)
      external        nf_inq_format

!
! dimension routines:
!

      integer         nf_def_dim
!                         (integer             ncid,
!                          character(*)        name,
!                          integer             len,
!                          integer             dimid)
      external        nf_def_dim

      integer         nf_inq_dimid
!                         (integer             ncid,
!                          character(*)        name,
!                          integer             dimid)
      external        nf_inq_dimid

      integer         nf_inq_dim
!                         (integer             ncid,
!                          integer             dimid,
!                          character(*)        name,
!                          integer             len)
      external        nf_inq_dim

      integer         nf_inq_dimname
!                         (integer             ncid,
!                          integer             dimid,
!                          character(*)        name)
      external        nf_inq_dimname

      integer         nf_inq_dimlen
!                         (integer             ncid,
!                          integer             dimid,
!                          integer             len)
      external        nf_inq_dimlen

      integer         nf_rename_dim
!                         (integer             ncid,
!                          integer             dimid,
!                          character(*)        name)
      external        nf_rename_dim

!
! general attribute routines:
!

      integer         nf_inq_att
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len)
      external        nf_inq_att

      integer         nf_inq_attid
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             attnum)
      external        nf_inq_attid

      integer         nf_inq_atttype
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype)
      external        nf_inq_atttype

      integer         nf_inq_attlen
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             len)
      external        nf_inq_attlen

      integer         nf_inq_attname
!                         (integer             ncid,
!                          integer             varid,
!                          integer             attnum,
!                          character(*)        name)
      external        nf_inq_attname

      integer         nf_copy_att
!                         (integer             ncid_in,
!                          integer             varid_in,
!                          character(*)        name,
!                          integer             ncid_out,
!                          integer             varid_out)
      external        nf_copy_att

      integer         nf_rename_att
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        curname,
!                          character(*)        newname)
      external        nf_rename_att

      integer         nf_del_att
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name)
      external        nf_del_att

!
! attribute put/get routines:
!

      integer         nf_put_att_text
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             len,
!                          character(*)        text)
      external        nf_put_att_text

      integer         nf_get_att_text
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          character(*)        text)
      external        nf_get_att_text

      integer         nf_put_att_int1
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          nf_int1_t           i1vals(1))
      external        nf_put_att_int1

      integer         nf_get_att_int1
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          nf_int1_t           i1vals(1))
      external        nf_get_att_int1

      integer         nf_put_att_int2
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          nf_int2_t           i2vals(1))
      external        nf_put_att_int2

      integer         nf_get_att_int2
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          nf_int2_t           i2vals(1))
      external        nf_get_att_int2

      integer         nf_put_att_int
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          integer             ivals(1))
      external        nf_put_att_int

      integer         nf_get_att_int
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             ivals(1))
      external        nf_get_att_int

      integer         nf_put_att_int64
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          nf_int8_t           i8vals(1))
      external        nf_put_att_int64

      integer         nf_get_att_int64
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          nf_int8_t           i8vals(1))
      external        nf_get_att_int64

      integer         nf_put_att_real
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          real                rvals(1))
      external        nf_put_att_real

      integer         nf_get_att_real
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          real                rvals(1))
      external        nf_get_att_real

      integer         nf_put_att_double
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          double              dvals(1))
      external        nf_put_att_double

      integer         nf_get_att_double
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          double              dvals(1))
      external        nf_get_att_double

!
! general variable routines:
!

      integer         nf_def_var
!                         (integer             ncid,
!                          character(*)        name,
!                          integer             datatype,
!                          integer             ndims,
!                          integer             dimids(1),
!                          integer             varid)
      external        nf_def_var

      integer         nf_inq_var
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             datatype,
!                          integer             ndims,
!                          integer             dimids(1),
!                          integer             natts)
      external        nf_inq_var

      integer         nf_inq_varid
!                         (integer             ncid,
!                          character(*)        name,
!                          integer             varid)
      external        nf_inq_varid

      integer         nf_inq_varname
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name)
      external        nf_inq_varname

      integer         nf_inq_vartype
!                         (integer             ncid,
!                          integer             varid,
!                          integer             xtype)
      external        nf_inq_vartype

      integer         nf_inq_varndims
!                         (integer             ncid,
!                          integer             varid,
!                          integer             ndims)
      external        nf_inq_varndims

      integer         nf_inq_vardimid
!                         (integer             ncid,
!                          integer             varid,
!                          integer             dimids(1))
      external        nf_inq_vardimid

      integer         nf_inq_varnatts
!                         (integer             ncid,
!                          integer             varid,
!                          integer             natts)
      external        nf_inq_varnatts

      integer         nf_rename_var
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name)
      external        nf_rename_var

      integer         nf_copy_var
!                         (integer             ncid_in,
!                          integer             varid,
!                          integer             ncid_out)
      external        nf_copy_var

!
! entire variable put/get routines:
!

      integer         nf_put_var_text
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        text)
      external        nf_put_var_text

      integer         nf_get_var_text
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        text)
      external        nf_get_var_text

      integer         nf_put_var_int1
!                         (integer             ncid,
!                          integer             varid,
!                          nf_int1_t           i1vals(1))
      external        nf_put_var_int1

      integer         nf_get_var_int1
!                         (integer             ncid,
!                          integer             varid,
!                          nf_int1_t           i1vals(1))
      external        nf_get_var_int1

      integer         nf_put_var_int2
!                         (integer             ncid,
!                          integer             varid,
!                          nf_int2_t           i2vals(1))
      external        nf_put_var_int2

      integer         nf_get_var_int2
!                         (integer             ncid,
!                          integer             varid,
!                          nf_int2_t           i2vals(1))
      external        nf_get_var_int2

      integer         nf_put_var_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             ivals(1))
      external        nf_put_var_int

      integer         nf_get_var_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             ivals(1))
      external        nf_get_var_int

      integer         nf_put_var_real
!                         (integer             ncid,
!                          integer             varid,
!                          real                rvals(1))
      external        nf_put_var_real

      integer         nf_get_var_real
!                         (integer             ncid,
!                          integer             varid,
!                          real                rvals(1))
      external        nf_get_var_real

      integer         nf_put_var_double
!                         (integer             ncid,
!                          integer             varid,
!                          doubleprecision     dvals(1))
      external        nf_put_var_double

      integer         nf_get_var_double
!                         (integer             ncid,
!                          integer             varid,
!                          doubleprecision     dvals(1))
      external        nf_get_var_double

!
! single variable put/get routines:
!

      integer         nf_put_var1_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          character*1         text)
      external        nf_put_var1_text

      integer         nf_get_var1_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          character*1         text)
      external        nf_get_var1_text

      integer         nf_put_var1_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          nf_int1_t           i1val)
      external        nf_put_var1_int1

      integer         nf_get_var1_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          nf_int1_t           i1val)
      external        nf_get_var1_int1

      integer         nf_put_var1_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          nf_int2_t           i2val)
      external        nf_put_var1_int2

      integer         nf_get_var1_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          nf_int2_t           i2val)
      external        nf_get_var1_int2

      integer         nf_put_var1_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          integer             ival)
      external        nf_put_var1_int

      integer         nf_get_var1_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          integer             ival)
      external        nf_get_var1_int

      integer         nf_put_var1_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          real                rval)
      external        nf_put_var1_real

      integer         nf_get_var1_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          real                rval)
      external        nf_get_var1_real

      integer         nf_put_var1_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          doubleprecision     dval)
      external        nf_put_var1_double

      integer         nf_get_var1_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          doubleprecision     dval)
      external        nf_get_var1_double

!
! variable array put/get routines:
!

      integer         nf_put_vara_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          character(*)        text)
      external        nf_put_vara_text

      integer         nf_get_vara_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          character(*)        text)
      external        nf_get_vara_text

      integer         nf_put_vara_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          nf_int1_t           i1vals(1))
      external        nf_put_vara_int1

      integer         nf_get_vara_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          nf_int1_t           i1vals(1))
      external        nf_get_vara_int1

      integer         nf_put_vara_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          nf_int2_t           i2vals(1))
      external        nf_put_vara_int2

      integer         nf_get_vara_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          nf_int2_t           i2vals(1))
      external        nf_get_vara_int2

      integer         nf_put_vara_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             ivals(1))
      external        nf_put_vara_int

      integer         nf_get_vara_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             ivals(1))
      external        nf_get_vara_int

      integer         nf_put_vara_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          real                rvals(1))
      external        nf_put_vara_real

      integer         nf_get_vara_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          real                rvals(1))
      external        nf_get_vara_real

      integer         nf_put_vara_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          doubleprecision     dvals(1))
      external        nf_put_vara_double

      integer         nf_get_vara_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          doubleprecision     dvals(1))
      external        nf_get_vara_double

!
! strided variable put/get routines:
!

      integer         nf_put_vars_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          character(*)        text)
      external        nf_put_vars_text

      integer         nf_get_vars_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          character(*)        text)
      external        nf_get_vars_text

      integer         nf_put_vars_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          nf_int1_t           i1vals(1))
      external        nf_put_vars_int1

      integer         nf_get_vars_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          nf_int1_t           i1vals(1))
      external        nf_get_vars_int1

      integer         nf_put_vars_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          nf_int2_t           i2vals(1))
      external        nf_put_vars_int2

      integer         nf_get_vars_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          nf_int2_t           i2vals(1))
      external        nf_get_vars_int2

      integer         nf_put_vars_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             ivals(1))
      external        nf_put_vars_int

      integer         nf_get_vars_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             ivals(1))
      external        nf_get_vars_int

      integer         nf_put_vars_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          real                rvals(1))
      external        nf_put_vars_real

      integer         nf_get_vars_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          real                rvals(1))
      external        nf_get_vars_real

      integer         nf_put_vars_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          doubleprecision     dvals(1))
      external        nf_put_vars_double

      integer         nf_get_vars_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          doubleprecision     dvals(1))
      external        nf_get_vars_double

!
! mapped variable put/get routines:
!

      integer         nf_put_varm_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          character(*)        text)
      external        nf_put_varm_text

      integer         nf_get_varm_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          character(*)        text)
      external        nf_get_varm_text

      integer         nf_put_varm_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          nf_int1_t           i1vals(1))
      external        nf_put_varm_int1

      integer         nf_get_varm_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          nf_int1_t           i1vals(1))
      external        nf_get_varm_int1

      integer         nf_put_varm_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          nf_int2_t           i2vals(1))
      external        nf_put_varm_int2

      integer         nf_get_varm_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          nf_int2_t           i2vals(1))
      external        nf_get_varm_int2

      integer         nf_put_varm_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          integer             ivals(1))
      external        nf_put_varm_int

      integer         nf_get_varm_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          integer             ivals(1))
      external        nf_get_varm_int

      integer         nf_put_varm_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          real                rvals(1))
      external        nf_put_varm_real

      integer         nf_get_varm_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          real                rvals(1))
      external        nf_get_varm_real

      integer         nf_put_varm_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          doubleprecision     dvals(1))
      external        nf_put_varm_double

      integer         nf_get_varm_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          doubleprecision     dvals(1))
      external        nf_get_varm_double

!     64-bit int functions.
      integer nf_put_var1_int64
      external nf_put_var1_int64
      integer nf_put_vara_int64
      external nf_put_vara_int64
      integer nf_put_vars_int64
      external nf_put_vars_int64
      integer nf_put_varm_int64
      external nf_put_varm_int64
      integer nf_put_var_int64
      external nf_put_var_int64
      integer nf_get_var1_int64
      external nf_get_var1_int64
      integer nf_get_vara_int64
      external nf_get_vara_int64
      integer nf_get_vars_int64
      external nf_get_vars_int64
      integer nf_get_varm_int64
      external nf_get_varm_int64
      integer nf_get_var_int64
      external nf_get_var_int64


!     NetCDF-4.
!     This is part of netCDF-4. Copyright 2006, UCAR, See COPYRIGHT
!     file for distribution information.

!     Netcdf version 4 fortran interface.

!     $Id: netcdf4.inc,v 1.28 2010/05/25 13:53:02 ed Exp $

!     New netCDF-4 types.
      integer nf_string
      integer nf_vlen
      integer nf_opaque
      integer nf_enum
      integer nf_compound

      parameter (nf_string = 12)
      parameter (nf_vlen = 13)
      parameter (nf_opaque = 14)
      parameter (nf_enum = 15)
      parameter (nf_compound = 16)

!     New netCDF-4 fill values.
      integer           nf_fill_ubyte
      integer           nf_fill_ushort
!      real              nf_fill_uint
!      real              nf_fill_int64
!      real              nf_fill_uint64
      parameter (nf_fill_ubyte = 255)
      parameter (nf_fill_ushort = 65535)

!     New constants.
      integer nf_format_netcdf4
      parameter (nf_format_netcdf4 = 3)

      integer nf_format_netcdf4_classic
      parameter (nf_format_netcdf4_classic = 4)

      integer nf_netcdf4
      parameter (nf_netcdf4 = 4096)

      integer nf_classic_model
      parameter (nf_classic_model = 256)

      integer nf_chunk_seq
      parameter (nf_chunk_seq = 0)
      integer nf_chunk_sub
      parameter (nf_chunk_sub = 1)
      integer nf_chunk_sizes
      parameter (nf_chunk_sizes = 2)

      integer nf_endian_native
      parameter (nf_endian_native = 0)
      integer nf_endian_little
      parameter (nf_endian_little = 1)
      integer nf_endian_big
      parameter (nf_endian_big = 2)

!     For NF_DEF_VAR_CHUNKING
      integer nf_chunked
      parameter (nf_chunked = 0)
      integer nf_contiguous
      parameter (nf_contiguous = 1)
      integer nf_compact
      parameter (nf_compact = 2)

!     For NF_DEF_VAR_FLETCHER32
      integer nf_nochecksum
      parameter (nf_nochecksum = 0)
      integer nf_fletcher32
      parameter (nf_fletcher32 = 1)

!     For NF_DEF_VAR_DEFLATE
      integer nf_noshuffle
      parameter (nf_noshuffle = 0)
      integer nf_shuffle
      parameter (nf_shuffle = 1)

!     For NF_DEF_VAR_SZIP
      integer nf_szip_ec_option_mask
      parameter (nf_szip_ec_option_mask = 4)
      integer nf_szip_nn_option_mask
      parameter (nf_szip_nn_option_mask = 32)

!     For parallel I/O.
      integer nf_mpiio      
      parameter (nf_mpiio = 8192)
      integer nf_mpiposix
      parameter (nf_mpiposix = 16384)
      integer nf_pnetcdf
      parameter (nf_pnetcdf = 32768)

!     For NF_VAR_PAR_ACCESS.
      integer nf_independent
      parameter (nf_independent = 0)
      integer nf_collective
      parameter (nf_collective = 1)

!     New error codes.
      integer nf_ehdferr        ! Error at HDF5 layer. 
      parameter (nf_ehdferr = -101)
      integer nf_ecantread      ! Can't read. 
      parameter (nf_ecantread = -102)
      integer nf_ecantwrite     ! Can't write. 
      parameter (nf_ecantwrite = -103)
      integer nf_ecantcreate    ! Can't create. 
      parameter (nf_ecantcreate = -104)
      integer nf_efilemeta      ! Problem with file metadata. 
      parameter (nf_efilemeta = -105)
      integer nf_edimmeta       ! Problem with dimension metadata. 
      parameter (nf_edimmeta = -106)
      integer nf_eattmeta       ! Problem with attribute metadata. 
      parameter (nf_eattmeta = -107)
      integer nf_evarmeta       ! Problem with variable metadata. 
      parameter (nf_evarmeta = -108)
      integer nf_enocompound    ! Not a compound type. 
      parameter (nf_enocompound = -109)
      integer nf_eattexists     ! Attribute already exists. 
      parameter (nf_eattexists = -110)
      integer nf_enotnc4        ! Attempting netcdf-4 operation on netcdf-3 file.   
      parameter (nf_enotnc4 = -111)
      integer nf_estrictnc3     ! Attempting netcdf-4 operation on strict nc3 netcdf-4 file.   
      parameter (nf_estrictnc3 = -112)
      integer nf_enotnc3        ! Attempting netcdf-3 operation on netcdf-4 file.   
      parameter (nf_enotnc3 = -113)
      integer nf_enopar         ! Parallel operation on file opened for non-parallel access.   
      parameter (nf_enopar = -114)
      integer nf_eparinit       ! Error initializing for parallel access.   
      parameter (nf_eparinit = -115)
      integer nf_ebadgrpid      ! Bad group ID.   
      parameter (nf_ebadgrpid = -116)
      integer nf_ebadtypid      ! Bad type ID.   
      parameter (nf_ebadtypid = -117)
      integer nf_etypdefined    ! Type has already been defined and may not be edited. 
      parameter (nf_etypdefined = -118)
      integer nf_ebadfield      ! Bad field ID.   
      parameter (nf_ebadfield = -119)
      integer nf_ebadclass      ! Bad class.   
      parameter (nf_ebadclass = -120)
      integer nf_emaptype       ! Mapped access for atomic types only.   
      parameter (nf_emaptype = -121)
      integer nf_elatefill      ! Attempt to define fill value when data already exists. 
      parameter (nf_elatefill = -122)
      integer nf_elatedef       ! Attempt to define var properties, like deflate, after enddef. 
      parameter (nf_elatedef = -123)
      integer nf_edimscale      ! Probem with HDF5 dimscales. 
      parameter (nf_edimscale = -124)
      integer nf_enogrp       ! No group found.
      parameter (nf_enogrp = -125)


!     New functions.

!     Parallel I/O.
      integer nf_create_par
      external nf_create_par

      integer nf_open_par
      external nf_open_par

      integer nf_var_par_access
      external nf_var_par_access

!     Functions to handle groups.
      integer nf_inq_ncid
      external nf_inq_ncid

      integer nf_inq_grps
      external nf_inq_grps

      integer nf_inq_grpname
      external nf_inq_grpname

      integer nf_inq_grpname_full
      external nf_inq_grpname_full

      integer nf_inq_grpname_len
      external nf_inq_grpname_len

      integer nf_inq_grp_parent
      external nf_inq_grp_parent

      integer nf_inq_grp_ncid
      external nf_inq_grp_ncid

      integer nf_inq_grp_full_ncid
      external nf_inq_grp_full_ncid

      integer nf_inq_varids
      external nf_inq_varids

      integer nf_inq_dimids
      external nf_inq_dimids

      integer nf_def_grp
      external nf_def_grp

!     New rename grp function

      integer nf_rename_grp
      external nf_rename_grp

!     New options for netCDF variables.
      integer nf_def_var_deflate
      external nf_def_var_deflate

      integer nf_inq_var_deflate
      external nf_inq_var_deflate

      integer nf_def_var_szip
      external nf_def_var_szip

      integer nf_inq_var_szip
      external nf_inq_var_szip

      integer nf_def_var_fletcher32
      external nf_def_var_fletcher32

      integer nf_inq_var_fletcher32
      external nf_inq_var_fletcher32

      integer nf_def_var_chunking
      external nf_def_var_chunking

      integer nf_inq_var_chunking
      external nf_inq_var_chunking

      integer nf_def_var_fill
      external nf_def_var_fill

      integer nf_inq_var_fill
      external nf_inq_var_fill

      integer nf_def_var_endian
      external nf_def_var_endian

      integer nf_inq_var_endian
      external nf_inq_var_endian

      integer nf_def_var_filter
      external nf_def_var_filter

      integer nf_inq_var_filter
      external nf_inq_var_filter

!     User defined types.
      integer nf_inq_typeids
      external nf_inq_typeids

      integer nf_inq_typeid
      external nf_inq_typeid

      integer nf_inq_type
      external nf_inq_type

      integer nf_inq_user_type
      external nf_inq_user_type

!     User defined types - compound types.
      integer nf_def_compound
      external nf_def_compound

      integer nf_insert_compound
      external nf_insert_compound

      integer nf_insert_array_compound
      external nf_insert_array_compound

      integer nf_inq_compound
      external nf_inq_compound

      integer nf_inq_compound_name
      external nf_inq_compound_name

      integer nf_inq_compound_size
      external nf_inq_compound_size

      integer nf_inq_compound_nfields
      external nf_inq_compound_nfields

      integer nf_inq_compound_field
      external nf_inq_compound_field

      integer nf_inq_compound_fieldname
      external nf_inq_compound_fieldname

      integer nf_inq_compound_fieldindex
      external nf_inq_compound_fieldindex

      integer nf_inq_compound_fieldoffset
      external nf_inq_compound_fieldoffset

      integer nf_inq_compound_fieldtype
      external nf_inq_compound_fieldtype

      integer nf_inq_compound_fieldndims
      external nf_inq_compound_fieldndims

      integer nf_inq_compound_fielddim_sizes
      external nf_inq_compound_fielddim_sizes

!     User defined types - variable length arrays.
      integer nf_def_vlen
      external nf_def_vlen

      integer nf_inq_vlen
      external nf_inq_vlen

      integer nf_free_vlen
      external nf_free_vlen

!     User defined types - enums.
      integer nf_def_enum
      external nf_def_enum

      integer nf_insert_enum
      external nf_insert_enum

      integer nf_inq_enum
      external nf_inq_enum

      integer nf_inq_enum_member
      external nf_inq_enum_member

      integer nf_inq_enum_ident
      external nf_inq_enum_ident

!     User defined types - opaque.
      integer nf_def_opaque
      external nf_def_opaque

      integer nf_inq_opaque
      external nf_inq_opaque

!     Write and read attributes of any type, including user defined
!     types.
      integer nf_put_att
      external nf_put_att
      integer nf_get_att
      external nf_get_att

!     Write and read variables of any type, including user defined
!     types.
      integer nf_put_var
      external nf_put_var
      integer nf_put_var1
      external nf_put_var1
      integer nf_put_vara
      external nf_put_vara
      integer nf_put_vars
      external nf_put_vars
      integer nf_get_var
      external nf_get_var
      integer nf_get_var1
      external nf_get_var1
      integer nf_get_vara
      external nf_get_vara
      integer nf_get_vars
      external nf_get_vars

!     For helping F77 users with VLENs.
      integer nf_get_vlen_element
      external nf_get_vlen_element
      integer nf_put_vlen_element
      external nf_put_vlen_element

!     For dealing with file level chunk cache.
      integer nf_set_chunk_cache
      external nf_set_chunk_cache
      integer nf_get_chunk_cache
      external nf_get_chunk_cache

!     For dealing with per variable chunk cache.
      integer nf_set_var_chunk_cache
      external nf_set_var_chunk_cache
      integer nf_get_var_chunk_cache
      external nf_get_var_chunk_cache

!     NetCDF-2.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! begin netcdf 2.4 backward compatibility:
!

!      
! functions in the fortran interface
!
      integer nccre
      integer ncopn
      integer ncddef
      integer ncdid
      integer ncvdef
      integer ncvid
      integer nctlen
      integer ncsfil

      external nccre
      external ncopn
      external ncddef
      external ncdid
      external ncvdef
      external ncvid
      external nctlen
      external ncsfil


      integer ncrdwr
      integer nccreat
      integer ncexcl
      integer ncindef
      integer ncnsync
      integer nchsync
      integer ncndirty
      integer nchdirty
      integer nclink
      integer ncnowrit
      integer ncwrite
      integer ncclob
      integer ncnoclob
      integer ncglobal
      integer ncfill
      integer ncnofill
      integer maxncop
      integer maxncdim
      integer maxncatt
      integer maxncvar
      integer maxncnam
      integer maxvdims
      integer ncnoerr
      integer ncebadid
      integer ncenfile
      integer nceexist
      integer nceinval
      integer nceperm
      integer ncenotin
      integer nceindef
      integer ncecoord
      integer ncemaxds
      integer ncename
      integer ncenoatt
      integer ncemaxat
      integer ncebadty
      integer ncebadd
      integer ncests
      integer nceunlim
      integer ncemaxvs
      integer ncenotvr
      integer nceglob
      integer ncenotnc
      integer ncfoobar
      integer ncsyserr
      integer ncfatal
      integer ncverbos
      integer ncentool


!
! netcdf data types:
!
      integer ncbyte
      integer ncchar
      integer ncshort
      integer nclong
      integer ncfloat
      integer ncdouble

      parameter(ncbyte = 1)
      parameter(ncchar = 2)
      parameter(ncshort = 3)
      parameter(nclong = 4)
      parameter(ncfloat = 5)
      parameter(ncdouble = 6)

!     
!     masks for the struct nc flag field; passed in as 'mode' arg to
!     nccreate and ncopen.
!     

!     read/write, 0 => readonly 
      parameter(ncrdwr = 1)
!     in create phase, cleared by ncendef 
      parameter(nccreat = 2)
!     on create destroy existing file 
      parameter(ncexcl = 4)
!     in define mode, cleared by ncendef 
      parameter(ncindef = 8)
!     synchronise numrecs on change (x'10')
      parameter(ncnsync = 16)
!     synchronise whole header on change (x'20')
      parameter(nchsync = 32)
!     numrecs has changed (x'40')
      parameter(ncndirty = 64)  
!     header info has changed (x'80')
      parameter(nchdirty = 128)
!     prefill vars on endef and increase of record, the default behavior
      parameter(ncfill = 0)
!     do not fill vars on endef and increase of record (x'100')
      parameter(ncnofill = 256)
!     isa link (x'8000')
      parameter(nclink = 32768)

!     
!     'mode' arguments for nccreate and ncopen
!     
      parameter(ncnowrit = 0)
      parameter(ncwrite = ncrdwr)
      parameter(ncclob = nf_clobber)
      parameter(ncnoclob = nf_noclobber)

!     
!     'size' argument to ncdimdef for an unlimited dimension
!     
      integer ncunlim
      parameter(ncunlim = 0)

!     
!     attribute id to put/get a global attribute
!     
      parameter(ncglobal  = 0)

!     
!     advisory maximums:
!     
      parameter(maxncop = 64)
      parameter(maxncdim = 1024)
      parameter(maxncatt = 8192)
      parameter(maxncvar = 8192)
!     not enforced 
      parameter(maxncnam = 256)
      parameter(maxvdims = maxncdim)

!     
!     global netcdf error status variable
!     initialized in error.c
!     

!     no error 
      parameter(ncnoerr = nf_noerr)
!     not a netcdf id 
      parameter(ncebadid = nf_ebadid)
!     too many netcdfs open 
      parameter(ncenfile = -31)   ! nc_syserr
!     netcdf file exists && ncnoclob
      parameter(nceexist = nf_eexist)
!     invalid argument 
      parameter(nceinval = nf_einval)
!     write to read only 
      parameter(nceperm = nf_eperm)
!     operation not allowed in data mode 
      parameter(ncenotin = nf_enotindefine )   
!     operation not allowed in define mode 
      parameter(nceindef = nf_eindefine)   
!     coordinates out of domain 
      parameter(ncecoord = nf_einvalcoords)
!     maxncdims exceeded 
      parameter(ncemaxds = nf_emaxdims)
!     string match to name in use 
      parameter(ncename = nf_enameinuse)   
!     attribute not found 
      parameter(ncenoatt = nf_enotatt)
!     maxncattrs exceeded 
      parameter(ncemaxat = nf_emaxatts)
!     not a netcdf data type 
      parameter(ncebadty = nf_ebadtype)
!     invalid dimension id 
      parameter(ncebadd = nf_ebaddim)
!     ncunlimited in the wrong index 
      parameter(nceunlim = nf_eunlimpos)
!     maxncvars exceeded 
      parameter(ncemaxvs = nf_emaxvars)
!     variable not found 
      parameter(ncenotvr = nf_enotvar)
!     action prohibited on ncglobal varid 
      parameter(nceglob = nf_eglobal)
!     not a netcdf file 
      parameter(ncenotnc = nf_enotnc)
      parameter(ncests = nf_ests)
      parameter (ncentool = nf_emaxname) 
      parameter(ncfoobar = 32)
      parameter(ncsyserr = -31)

!     
!     global options variable. used to determine behavior of error handler.
!     initialized in lerror.c
!     
      parameter(ncfatal = 1)
      parameter(ncverbos = 2)

!
!     default fill values.  these must be the same as in the c interface.
!
      integer filbyte
      integer filchar
      integer filshort
      integer fillong
      real filfloat
      doubleprecision fildoub

      parameter (filbyte = -127)
      parameter (filchar = 0)
      parameter (filshort = -32767)
      parameter (fillong = -2147483647)
      parameter (filfloat = 9.9692099683868690e+36)
      parameter (fildoub = 9.9692099683868690e+36)

!     This is to turn on netCDF internal logging.
      integer nf_set_log_level
      external nf_set_log_level
!#include"advtrac.h"
c=======================================================================
c   Declarations
c=======================================================================

      INTEGER,INTENT(IN) :: ngrid, nlayer

c Old variables dimensions (from file)
c------------------------------------
      INTEGER   imold,jmold,lmold,nsoilold,nqold


c Variables pour les lectures des fichiers "ini" 
c--------------------------------------------------
!      INTEGER sizei,
      integer timelen,dimid
!      INTEGER length
!      parameter (length = 100)
      INTEGER tab0
      INTEGER isoil,iq,iqmax
      CHARACTER*2   str2

!      REAL dimfirst(4) ! tableau contenant les 1ers elements des dimensions

!      REAL dimlast(4) ! tableau contenant les derniers elements des dimensions

!      REAL dimcycl(4) ! tableau contenant les periodes des dimensions
!      CHARACTER*120 dimsource
!      CHARACTER*16 dimname
!      CHARACTER*80 dimtitle
!      CHARACTER*40 dimunits
!      INTEGER   dimtype

!      INTEGER dimord(4)  ! tableau contenant l''ordre
!      data dimord /1,2,3,4/ ! de sortie des dimensions

!      INTEGER vardim(4)
      REAL date
      INTEGER   memo
!      character (len=50) :: tmpname

c Variable histoire 
c------------------
      REAL vcov(iip1,jjm,llm),ucov(iip1,jjp1,llm) ! vents covariants
      REAL h(iip1,jjp1,llm),ps(iip1,jjp1)
      REAL q(iip1,jjp1,llm,nqtot),qtot(iip1,jjp1,llm)

c autre variables dynamique nouvelle grille
c------------------------------------------

c!-*-
!      integer klatdat,klongdat
!      PARAMETER (klatdat=180,klongdat=360)

c Physique sur grille scalaire 
c----------------------------

c variable physique
c------------------
      REAL tsurf(ngrid) ! surface temperature
      REAL tsoil(ngrid,nsoilmx) ! soil temperature
      REAL co2ice(ngrid) ! CO2 ice layer
      REAL emis(ngrid)
      REAL q2(ngrid,llm+1),qsurf(ngrid,nqtot)
      REAL tslab(ngrid,noceanmx)
      REAL rnat(ngrid),pctsrf_sic(ngrid)
      REAL tsea_ice(ngrid),sea_ice(ngrid)
c     REAL phisfi(ngrid)

      INTEGER i,j,l
      INTEGER nid,nvarid
c     REAL year_day,periheli,aphelie,peri_day
c     REAL obliquit,z0,emin_turb,lmixmin
c     REAL emissiv,emisice(2),albedice(2)
c     REAL iceradius(2) , dtemisice(2)

      integer ierr
      integer, dimension(4) :: start,count

c Variable nouvelle grille naturelle au point scalaire
c------------------------------------------------------
      real us(iip1,jjp1,llm),vs(iip1,jjp1,llm)
      REAL phisold_newgrid(iip1,jjp1)
      REAL t(iip1,jjp1,llm)
      real tsurfS(iip1,jjp1),tsoilS(iip1,jjp1,nsoilmx)
      real inertiedatS(iip1,jjp1,nsoilmx)
      real co2iceS(iip1,jjp1)
      real emisS(iip1,jjp1)
      REAL q2S(iip1,jjp1,llm+1),qsurfS(iip1,jjp1,nqtot)
      real tslabS(iip1,jjp1,noceanmx)
      real pctsrf_sicS(iip1,jjp1),tsea_iceS(iip1,jjp1)
      real rnatS(iip1,jjp1), sea_iceS(iip1,jjp1)

      real ptotal, co2icetotal

c Var intermediaires : vent naturel, mais pas coord scalaire
c-----------------------------------------------------------
      real vnat(iip1,jjm,llm),unat(iip1,jjp1,llm)


c Variable de l'ancienne grille 
c---------------------------------------------------------

      real, dimension(:), allocatable :: timelist
      real, dimension(:), allocatable :: rlonuold, rlatvold
      real, dimension(:), allocatable :: rlonvold, rlatuold
      real, dimension(:), allocatable :: apsold,bpsold
      real, dimension(:), allocatable :: mlayerold
      real, dimension(:,:,:), allocatable :: uold,vold,told,q2old
      real, dimension(:,:,:), allocatable :: tsoilold,qsurfold
      real, dimension(:,:,:),allocatable :: tsoiloldnew
! tsoiloldnew: old soil values, but along new subterranean grid
      real, dimension(:,:,:), allocatable :: inertiedatold
! inertiedatoldnew: old inertia values, but along new subterranean grid
      real, dimension(:,:,:), allocatable :: inertiedatoldnew
      real, dimension(:,:), allocatable :: psold,phisold
      real, dimension(:,:), allocatable :: co2iceold
      real, dimension(:,:), allocatable :: tsurfold
      real, dimension(:,:), allocatable :: emisold
      real, dimension(:,:,:,:), allocatable :: qold
      real, dimension(:,:,:), allocatable :: tslabold
      real, dimension(:,:), allocatable :: rnatold,pctsrf_sicold
      real, dimension(:,:), allocatable :: tsea_iceold,sea_iceold


      real tab_cntrl(100)

      real ptotalold, co2icetotalold

      logical :: olddepthdef=.false. ! flag
! olddepthdef=.true. if soil depths are in 'old' (unspecified) format
      logical :: depthinterpol=.false. ! flag
! depthinterpol=.true. if interpolation will be requiered
      logical :: therminertia_3D=.true. ! flag
! therminertia_3D=.true. if thermal inertia is 3D and read from datafile
c Variable intermediaires iutilise pour l'extrapolation verticale 
c----------------------------------------------------------------
      real, dimension(:,:,:), allocatable :: var,varp1 
      real, dimension(:), allocatable :: oldgrid, oldval
      real, dimension(:), allocatable :: newval

      real surfith(iip1,jjp1) ! surface thermal inertia
      ! surface thermal inertia at old horizontal grid resolution
      real, dimension(:,:), allocatable :: surfithold 

      character(len=30) :: txt ! to store some text

c=======================================================================

! 0. Preliminary stuff



!-----------------------------------------------------------------------
! 1. Read data dimensions (i.e. size and length)
!-----------------------------------------------------------------------

! 1.2 Read the various dimension lengths of data in file 

      ierr= NF_INQ_DIMID(nid,"Time",dimid)
      if (ierr.ne.NF_NOERR) then
         ierr= NF_INQ_DIMID(nid,"temps",dimid)
      endif
      ierr= NF_INQ_DIMLEN(nid,dimid,timelen)

      ierr= NF_INQ_DIMID(nid,"latitude",dimid)
      if (ierr.ne.NF_NOERR) then
         ierr= NF_INQ_DIMID(nid,"rlatu",dimid)
      endif
      ierr= NF_INQ_DIMLEN(nid,dimid,jmold)
      jmold=jmold-1

      ierr= NF_INQ_DIMID(nid,"longitude",dimid)
      if (ierr.ne.NF_NOERR) then
         ierr= NF_INQ_DIMID(nid,"rlonv",dimid)
      endif
      ierr= NF_INQ_DIMLEN(nid,dimid,imold)
      imold=imold-1

      ierr= NF_INQ_DIMID(nid,"altitude",dimid)
      if (ierr.ne.NF_NOERR) then
         ierr= NF_INQ_DIMID(nid,"sig_s",dimid)
      endif
      ierr= NF_INQ_DIMLEN(nid,dimid,lmold)

      nqold=0
      do
         write(str2,'(i2.2)') nqold+1
         ierr= NF_INQ_VARID(nid,'q'//str2,dimid)
!        write(*,*) 'q'//str2
         if (ierr.eq.NF_NOERR) then
            nqold=nqold+1
         else
            exit
         endif
      enddo

! 1.2.1 find out the # of subsurface_layers
      nsoilold=0 !dummy initialisation
      ierr=NF_INQ_DIMID(nid,"subsurface_layers",dimid)
      if (ierr.eq.NF_NOERR) then
        ierr=NF_INQ_DIMLEN(nid,dimid,nsoilold)
	if (ierr.ne.NF_NOERR) then
	 write(*,*)'lec_start_archive: ',
     &              'Failed reading subsurface_layers length'
	endif
      else
        write(*,*)"lec_start_archive: did not find subsurface_layers"
      endif

      if (nsoilold.eq.0) then ! 'old' archive format;
      ! must use Tg//str2 fields to compute nsoilold
      write(*,*)"lec_start_archive: building nsoilold from Tg fields"
        do
	 write(str2,'(i2.2)') nsoilold+1
	 ierr=NF_INQ_VARID(nid,'Tg'//str2,dimid)
	 if (ierr.eq.NF_NOERR) then
	  nsoilold=nsoilold+1
	 else
	  exit
	 endif
	enddo
      endif


      if (nsoilold.ne.nsoilmx) then ! interpolation will be required
        depthinterpol=.true.
      endif

! 1.3 Report dimensions
      
      write(*,*) "Start_archive dimensions:"
      write(*,*) "longitude: ",imold
      write(*,*) "latitude: ",jmold
      write(*,*) "altitude: ",lmold
      write(*,*) "tracers: ",nqold
      write(*,*) "subsurface_layers: ",nsoilold
      if (depthinterpol) then
      write(*,*) " => Warning, nsoilmx= ",nsoilmx
      write(*,*) '    which implies that you want subterranean interpola
     &tion.'
      write(*,*) '  Otherwise, set nsoilmx -in dimphys.h- to: ',nsoilold
      endif
      write(*,*) "time lenght: ",timelen
      write(*,*) 

!-----------------------------------------------------------------------
! 2. Allocate arrays to store datasets
!-----------------------------------------------------------------------

      allocate(timelist(timelen))
      allocate(rlonuold(imold+1), rlatvold(jmold))
      allocate(rlonvold(imold+1), rlatuold(jmold+1))
      allocate (apsold(lmold),bpsold(lmold))
      allocate(uold(imold+1,jmold+1,lmold))
      allocate(vold(imold+1,jmold+1,lmold))
      allocate(told(imold+1,jmold+1,lmold))
      allocate(psold(imold+1,jmold+1))
      allocate(phisold(imold+1,jmold+1))
      allocate(qold(imold+1,jmold+1,lmold,nqtot))
      allocate(co2iceold(imold+1,jmold+1))
      allocate(tsurfold(imold+1,jmold+1))
      allocate(emisold(imold+1,jmold+1))
      allocate(q2old(imold+1,jmold+1,lmold+1))
!      allocate(tsoilold(imold+1,jmold+1,nsoilmx))
      allocate(tsoilold(imold+1,jmold+1,nsoilold))
      allocate(tsoiloldnew(imold+1,jmold+1,nsoilmx))
      allocate(inertiedatold(imold+1,jmold+1,nsoilold)) ! soil thermal inertia
      allocate(inertiedatoldnew(imold+1,jmold+1,nsoilmx))
      ! surface thermal inertia at old horizontal grid resolution
      allocate(surfithold(imold+1,jmold+1))
      allocate(mlayerold(nsoilold))
      allocate(qsurfold(imold+1,jmold+1,nqtot))
      allocate(tslabold(imold+1,jmold+1,noceanmx))
      allocate(rnatold(imold+1,jmold+1))
      allocate(pctsrf_sicold(imold+1,jmold+1))
      allocate(tsea_iceold(imold+1,jmold+1))
      allocate(sea_iceold(imold+1,jmold+1))

      allocate(var (imold+1,jmold+1,llm))
      allocate(varp1 (imold+1,jmold+1,llm+1))

      write(*,*) 'q2',ngrid,llm+1
      write(*,*) 'q2S',iip1,jjp1,llm+1
      write(*,*) 'q2old',imold+1,jmold+1,lmold+1

!-----------------------------------------------------------------------
! 3. Read time-independent data
!-----------------------------------------------------------------------

C-----------------------------------------------------------------------
c 3.1. Lecture du tableau des parametres du run 
c     (pour  la lecture ulterieure de "ptotalold" et "co2icetotalold")
c-----------------------------------------------------------------------
c
      ierr = NF_INQ_VARID (nid, "controle", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "Lect_start_archive: champ <controle> est absent"
         CALL abort
      ENDIF
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, tab_cntrl)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echoue pour <controle>"
         CALL abort
      ENDIF
c
      tab0 = 50

c-----------------------------------------------------------------------
c 3.2 Lecture des longitudes et latitudes
c-----------------------------------------------------------------------
c
      ierr = NF_INQ_VARID (nid, "rlonv", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <rlonv> est absent"
         CALL abort
      ENDIF
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, rlonvold)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <rlonv>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "rlatu", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <rlatu> est absent"
         CALL abort
      ENDIF 
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, rlatuold)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <rlatu>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "rlonu", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <rlonu> est absent"
         CALL abort
      ENDIF
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, rlonuold)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <rlonu>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "rlatv", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <rlatv> est absent"
         CALL abort
      ENDIF
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, rlatvold)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <rlatv>"
         CALL abort
      ENDIF
c

c-----------------------------------------------------------------------
c 3.3. Lecture des niveaux verticaux
c-----------------------------------------------------------------------
c
      ierr = NF_INQ_VARID (nid, "aps", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <aps> est absent"
         apsold=0
         PRINT*, "<aps> set to 0"
      ELSE
         ierr = NF_GET_VAR_DOUBLE(nid, nvarid, apsold)
         IF (ierr .NE. NF_NOERR) THEN
            PRINT*, "lect_start_archive: Lecture echouee pour <aps>"
         ENDIF
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "bps", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <bps> est absent"
         PRINT*, "It must be an old start_archive, lets look for sig_s"
         ierr = NF_INQ_VARID (nid, "sig_s", nvarid)
         IF (ierr .NE. NF_NOERR) THEN
            PRINT*, "Nothing to do..."
            CALL abort
         ENDIF
      ENDIF
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, bpsold)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <bps>"
         CALL abort
      END IF

c-----------------------------------------------------------------------
c 3.4 Read Soil layers depths
c-----------------------------------------------------------------------
     
      ierr=NF_INQ_VARID(nid,"soildepth",nvarid)
      if (ierr.ne.NF_NOERR) then
       write(*,*)'lect_start_archive: Could not find <soildepth>'
       write(*,*)' => Assuming this is an archive in old format'
       olddepthdef=.true.
       depthinterpol=.true.
       ! this is how soil depth was defined in ye old days
	do isoil=1,nsoilold
	  mlayerold(isoil)=sqrt(887.75/3.14)*((2.**(isoil-0.5))-1.)
	enddo
      else
        ierr = NF_GET_VAR_DOUBLE(nid,nvarid,mlayerold)
       if (ierr .NE. NF_NOERR) then
         PRINT*, "lect_start_archive: Failed reading <soildepth>"
         CALL abort
       endif

      endif !of if(ierr.ne.NF_NOERR)

      ! Read (or build) mlayer()
      if (depthinterpol) then
       ! Build (default) new soil depths (mlayer(:) is in comsoil.h),
       ! as in soil_settings.F
       write(*,*)' => Building default soil depths'
       do isoil=0,nsoilmx-1
         mlayer(isoil)=2.e-4*(2.**(isoil-0.5))
       enddo
       write(*,*)' => mlayer: ',mlayer
       ! Also build (default) new soil interlayer depth layer(:)
       do isoil=1,nsoilmx
         layer(isoil)=sqrt(mlayer(0)*mlayer(1))*
     &                      ((mlayer(1)/mlayer(0))**(isoil-1))
       enddo
       write(*,*)' =>  layer: ',layer
      else ! read mlayer() from file
        ierr = NF_GET_VAR_DOUBLE(nid,nvarid,mlayer)
       if (ierr .NE. NF_NOERR) then
         PRINT*, "lect_start_archive: Failed reading <soildepth>"
         CALL abort
       endif
      endif ! of if (depthinterpol)

c-----------------------------------------------------------------------
c 3.5 Read Soil thermal inertia
c-----------------------------------------------------------------------

      ierr=NF_INQ_VARID(nid,"inertiedat",nvarid)
      if (ierr.ne.NF_NOERR) then
       write(*,*)'lect_start_archive: Could not find <inertiedat>'
       write(*,*)' => Assuming this is an archive in old format'
       therminertia_3D=.false.
       write(*,*)' => Thermal inertia will be read from reference file'
       volcapa=1.e6
       write(*,*)'    and soil volumetric heat capacity is set to ',
     &           volcapa
      else
        ierr = NF_GET_VAR_DOUBLE(nid,nvarid,inertiedatold)
       if (ierr .NE. NF_NOERR) then
         PRINT*, "lect_start_archive: Failed reading <inertiedat>"
         CALL abort
       endif
      endif

c-----------------------------------------------------------------------
c 3.6 Lecture geopotentiel au sol
c-----------------------------------------------------------------------
c
      ierr = NF_INQ_VARID (nid, "phisinit", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <phisinit> est absent"
         CALL abort
      ENDIF
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, phisold)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <phisinit>"
         CALL abort
      ENDIF

C-----------------------------------------------------------------------
c   lecture de "ptotalold" et "co2icetotalold"
c-----------------------------------------------------------------------
      ptotalold = tab_cntrl(tab0+49)
      co2icetotalold = tab_cntrl(tab0+50)
 
c-----------------------------------------------------------------------
c 4. Lecture du temps et choix
c-----------------------------------------------------------------------
 
c  lecture du temps
c
      ierr = NF_INQ_DIMID (nid, "Time", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         ierr = NF_INQ_DIMID (nid, "temps", nvarid)
         IF (ierr .NE. NF_NOERR) THEN
            PRINT*, "lect_start_archive: Le champ <Time> est absent"
            CALL abort
         endif
      ENDIF

      ierr = NF_INQ_VARID (nid, "Time", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         ierr = NF_INQ_VARID (nid, "temps", nvarid)
      endif 
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, timelist)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <Time>"
         CALL abort
      ENDIF
c
      write(*,*)
      write(*,*)
      write(*,*) 'Differentes dates des etats initiaux stockes:'
      write(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      pi=2.*ASIN(1.)
      do i=1,timelen
c       call solarlong(timelist(i),sollong(i))
c       sollong(i) = sollong(i)*180./pi
        write(*,*) 'etat initial au jour martien' ,int(timelist(i))
c       write(*,6) nint(timelist(i)),nint(mod(timelist(i),669)),
c    .    sollong(i)
      end do

   6  FORMAT(i7,i7,f9.3)
 
      write(*,*)
      write(*,*) 'Choix de la date'
 123  read(*,*,iostat=ierr) date
      if(ierr.ne.0) goto 123
      memo = 0
      do i=1,timelen
        if (date.eq.int(timelist(i))) then
            memo = i
        endif
      end do
 
      if (memo.eq.0) then
        write(*,*)
        write(*,*)
        write(*,*) 'He alors... Y sait pas lire !?!'
        write(*,*)
        write(*,*) 'Differentes dates des etats initiaux stockes:'
        write(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        do i=1,timelen
          write(*,*) 'etat initial au jour martien' ,nint(timelist(i))
c         write(*,6) nint(timelist(i)),nint(mod(timelist(i),669))
        end do
        goto 123
      endif

!-----------------------------------------------------------------------
! 5. Read (time-dependent) data from datafile
!-----------------------------------------------------------------------


c-----------------------------------------------------------------------
c 5.1 Lecture des champs 2D (co2ice, emis,ps,tsurf,Tg[10], qsurf)
c-----------------------------------------------------------------------
 
      start=(/1,1,memo,0/)
      count=(/imold+1,jmold+1,1,0/)
       
! CO2ice is now in qsurf(igcm_co2_ice) ...
!      ierr = NF_INQ_VARID (nid, "co2ice", nvarid)
!      IF (ierr .NE. NF_NOERR) THEN
!         PRINT*, "lect_start_archive: Le champ <co2ice> est absent"
!         CALL abort
!      ENDIF
!#ifdef 1
!      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,co2iceold)
!#else
!      ierr = NF_GET_VARA_REAL(nid, nvarid,start,count,co2iceold)
!#endif
!      IF (ierr .NE. NF_NOERR) THEN
!         PRINT*, "lect_start_archive: Lecture echouee pour <co2ice>"
!         PRINT*, NF_STRERROR(ierr)
!         CALL abort
!      ENDIF
c
      ierr = NF_INQ_VARID (nid, "emis", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <emis> est absent"
         CALL abort
      ENDIF
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,emisold)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <emis>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "ps", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <ps> est absent"
         CALL abort
      ENDIF
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,psold)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <ps>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "tsurf", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <tsurf> est absent"
         CALL abort
      ENDIF
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,tsurfold)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <tsurf>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "q2surf", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <q2surf> est absent"
         CALL abort
      ENDIF
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,q2old)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <q2surf>"
         CALL abort
      ENDIF
c
cc Slab ocean
      if(ok_slab_ocean) then
      start=(/1,1,1,memo/)
      count=(/imold+1,jmold+1,noceanmx,1/)

       ierr=NF_INQ_VARID(nid,"tslab",nvarid)
       IF (ierr.ne.NF_NOERR) then
          PRINT*,"lect_start_archive: Cannot find <tslab>"
       ENDIF
      ierr=NF_GET_VARA_DOUBLE(nid,nvarid,start,count,tslabold)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <tslab>"
      ENDIF


c
      start=(/1,1,memo,0/)
      count=(/imold+1,jmold+1,1,0/)

      ierr = NF_INQ_VARID (nid, "rnat", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <rnat> est absent"
      ENDIF
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,rnatold)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <rnat>"
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "pctsrf_sic", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <pctsrf_sic> est absent"
      ENDIF
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,pctsrf_sicold)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <pctsrf_sic>"
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "tsea_ice", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <tsea_ice> est absent"
      ENDIF
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,tsea_iceold)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <tsea_ice>"
      ENDIF
c
      ierr = NF_INQ_VARID (nid, "sea_ice", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <sea_ice> est absent"
      ENDIF
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,sea_iceold)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <sea_ice>"
      ENDIF
  
      ENDIF! ok_slab_ocean
c
      write(*,*)"lect_start_archive: rlonuold:"
     &           ,rlonuold," rlatvold:",rlatvold
      write(*,*)

! Surface tracers:      
      do iq=1,nqtot
        ! initialize all surface tracers to zero
        call initial0((jmold+1)*(imold+1), qsurfold(1,1,iq))
      enddo


!      print*,'tname=',tname
!      print*,'nid',nid
!      print*,'nvarid',nvarid
!      stop

      DO iq=1,nqtot
          txt=trim(tname(iq))//"_surf"
          if (txt.eq."h2o_vap") then
            ! There is no surface tracer for h2o_vap;
            ! "h2o_ice" should be loaded instead
            txt="h2o_ice_surf"
            write(*,*) 'lect_start_archive: loading surface tracer',
     &                     ' h2o_ice instead of h2o_vap'
          endif

        
        write(*,*) "lect_start_archive: loading tracer ",trim(txt)
        ierr = NF_INQ_VARID (nid,txt,nvarid)
        IF (ierr .NE. NF_NOERR) THEN
          PRINT*, "lect_start_archive: ",
     &              " Tracer <",trim(txt),"> not found"

!          print*,'RDW has added hack to let me continue...'
!          CALL abort
        ENDIF
        ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,
     &          qsurfold(1,1,iq))
        IF (ierr .NE. NF_NOERR) THEN
          PRINT*, "lect_start_archive: ",
     &             " Failed loading <",trim(txt),">"
          write (*,*) trim(txt),'    is set to 0'
!          call initial0((jmold+1)*(imold+1), qsurfold(1,1,iq))
        ENDIF

      ENDDO ! of DO iq=1,nqtot


!-----------------------------------------------------------------------
! 5.2 Read 3D subterranean fields
!-----------------------------------------------------------------------

      start=(/1,1,1,memo/)
      count=(/imold+1,jmold+1,nsoilold,1/)
!
! Read soil temperatures
!
      if (olddepthdef) then ! tsoil stored using the 'old format'
         start=(/1,1,memo,0/)
         count=(/imold+1,jmold+1,1,0/) ! because the "Tg" are 2D datasets
       do isoil=1,nsoilold
         write(str2,'(i2.2)') isoil
c
         ierr = NF_INQ_VARID (nid, "Tg"//str2, nvarid)
         IF (ierr .NE. NF_NOERR) THEN
            PRINT*, "lect_start_archive: ",
     &              "Field <","Tg"//str2,"> not found"
            CALL abort
         ENDIF
         ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,
     &          tsoilold(1,1,isoil))
         IF (ierr .NE. NF_NOERR) THEN
            PRINT*, "lect_start_archive: ",
     &            "Failed reading <","Tg"//str2,">"
            CALL abort
         ENDIF
c
       enddo ! of do isoil=1,nsoilold
      
      ! reset 'start' and 'count' to "3D" behaviour
      start=(/1,1,1,memo/)
      count=(/imold+1,jmold+1,nsoilold,1/)
      
      else
       write(*,*) "lect_start_archive: loading tsoil "
       ierr=NF_INQ_VARID(nid,"tsoil",nvarid)
       if (ierr.ne.NF_NOERR) then
        write(*,*)"lect_start_archive: Cannot find <tsoil>"
	call abort
       else
      ierr=NF_GET_VARA_DOUBLE(nid,nvarid,start,count,tsoilold)
       endif ! of if (ierr.ne.NF_NOERR)
       
      endif ! of if (olddepthdef)

!
! Read soil thermal inertias
!
!      if (.not.olddepthdef) then ! no thermal inertia data in "old" archives
!       ierr=NF_INQ_VARID(nid,"inertiedat",nvarid)
!       if (ierr.ne.NF_NOERR) then
!        write(*,*)"lect_start_archive: Cannot find <inertiedat>"
!	call abort
!       else
!#ifdef 1
!      ierr=NF_GET_VARA_DOUBLE(nid,nvarid,start,count,inertiedatold)
!#else
!      ierr=NF_GET_VARA_REAL(nid,nvarid,start,count,inertiedatold)
!#endif
!       endif ! of if (ierr.ne.NF_NOERR)
!      endif

c-----------------------------------------------------------------------
c 5.3	Lecture des champs 3D (t,u,v, q2atm,q)
c-----------------------------------------------------------------------

      start=(/1,1,1,memo/)
      count=(/imold+1,jmold+1,lmold,1/)

c
      ierr = NF_INQ_VARID (nid,"temp", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <temp> est absent"
         CALL abort
      ENDIF
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid, start, count, told)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <temp>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid,"u", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <u> est absent"
         CALL abort
      ENDIF
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,uold)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <u>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid,"v", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <v> est absent"
         CALL abort
      ENDIF
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,vold)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <v>"
         CALL abort
      ENDIF
c
      ierr = NF_INQ_VARID (nid,"q2atm", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Le champ <q2atm> est absent"
         CALL abort
      ENDIF
      ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start,count,q2old(1,1,2))
      IF (ierr .NE. NF_NOERR) THEN
         PRINT*, "lect_start_archive: Lecture echouee pour <q2atm>"
         CALL abort
      ENDIF
c

! Tracers:      
      do iq=1,nqtot
         call initial0((jmold+1)*(imold+1)*lmold,qold(1,1,1,iq) )
      enddo

      DO iq=1,nqtot
        txt=tname(iq)
        write(*,*)"lect_start_archive: loading tracer ",trim(txt)
        ierr = NF_INQ_VARID (nid,txt,nvarid)
        IF (ierr .NE. NF_NOERR) THEN
            PRINT*, "lect_start_archive: ",
     &              " Tracer <",trim(txt),"> not found"
!            CALL abort
        ENDIF
        ierr=NF_GET_VARA_DOUBLE(nid,nvarid,start,count,qold(1,1,1,iq))
        IF (ierr .NE. NF_NOERR) THEN
          PRINT*, "lect_start_archive: ",
     &             "  Failed loading <",trim(txt),">"
          write (*,*) trim(txt),'      set to 1.E-30'
          do l=1,lmold
            do j=1,jmold+1
              do i=1,imold+1
                 qold(i,j,l,iq)=1.e-30
              end do
            end do
          end do
        ENDIF

      ENDDO ! of DO iq=1,nqtot


!=======================================================================
! 6. Interpolation from old grid to new grid
!=======================================================================

c=======================================================================
c   INTERPOLATION DANS LA NOUVELLE GRILLE et initialisation des variables
c=======================================================================
c  Interpolation horizontale puis passage dans la grille physique pour 
c  les variables physique 
c  Interpolation verticale puis horizontale pour chaque variable 3D
c=======================================================================

c-----------------------------------------------------------------------
c 6.1	Variable 2d :
c-----------------------------------------------------------------------
c Relief 
      call interp_horiz (phisold,phisold_newgrid,imold,jmold,iim,jjm,1,
     &                   rlonuold,rlatvold,rlonu,rlatv)

! CO2 ice is now in qsurf(igcm_co2_ice)
!      call interp_horiz (co2iceold,co2ices,imold,jmold,iim,jjm,1,
!     &                   rlonuold,rlatvold,rlonu,rlatv)

c Temperature de surface
      call interp_horiz (tsurfold,tsurfs,imold,jmold,iim,jjm,1,
     &                   rlonuold,rlatvold,rlonu,rlatv)
      call gr_dyn_fi (1,iim+1,jjm+1,ngrid,tsurfs,tsurf)
c     write(44,*) 'tsurf', tsurf

c Temperature du sous-sol
!      call interp_horiz(tsoilold,tsoils,
!     &                  imold,jmold,iim,jjm,nsoilmx,
!     &                   rlonuold,rlatvold,rlonu,rlatv)
!      call gr_dyn_fi (nsoilmx,iim+1,jjm+1,ngrid,tsoils,tsoil)
c     write(45,*) 'tsoil',tsoil

c Emissivite de la surface
      call interp_horiz (emisold,emiss,imold,jmold,iim,jjm,1,
     &                   rlonuold,rlatvold,rlonu,rlatv)
      call gr_dyn_fi (1,iim+1,jjm+1,ngrid,emiss,emis)
c     write(46,*) 'emis',emis



c-----------------------------------------------------------------------
c 6.1.2	Traitement special de la pression au sol :
c-----------------------------------------------------------------------

c  Extrapolation la pression dans la nouvelle grille
      call interp_horiz(psold,ps,imold,jmold,iim,jjm,1,
     &                   rlonuold,rlatvold,rlonu,rlatv)

c-----------------------------------------------------------------------
c	On assure la conservation de la masse de l'atmosphere + calottes
c-----------------------------------------------------------------------

      ptotal =  0.
      DO j=1,jjp1
         DO i=1,iim
            ptotal=ptotal+ps(i,j)*aire(i,j)/g
         END DO
      END DO
      co2icetotal = 0.
      if (igcm_co2_ice.ne.0) then
        ! recast surface CO2 ice on new grid
        call interp_horiz(qsurfold(1,1,igcm_co2_ice),
     &                  qsurfs(1,1,igcm_co2_ice),
     &                  imold,jmold,iim,jjm,1,
     &                  rlonuold,rlatvold,rlonu,rlatv)
       DO j=1,jjp1
         DO i=1,iim
           !co2icetotal = co2icetotal + co2iceS(i,j)*aire(i,j)
           co2icetotal=co2icetotal+qsurfS(i,j,igcm_co2_ice)*aire(i,j)
         END DO
       END DO
      else
        write(*,*) "Warning: No co2_ice surface tracer"
      endif

      write(*,*)
      write(*,*)'Ancienne grille: masse de l atm :',ptotalold
      write(*,*)'Nouvelle grille: masse de l atm :',ptotal
      write (*,*) 'Ratio new atm./ old atm =', ptotal/ptotalold 
      write(*,*)
      write(*,*)'Ancienne grille: masse de la glace CO2:',co2icetotalold
      write(*,*)'Nouvelle grille: masse de la glace CO2:',co2icetotal
      if (co2icetotalold.ne.0.) then
      write(*,*)'Ratio new ice./old ice =',co2icetotal/co2icetotalold
      endif
      write(*,*)


      DO j=1,jjp1
         DO i=1,iip1
            ps(i,j)=ps(i,j) * ptotalold/ptotal
         END DO
      END DO

      if ( co2icetotalold.gt.0.) then 
!         DO j=1,jjp1
!            DO i=1,iip1
!               co2iceS(i,j)=co2iceS(i,j) * co2icetotalold/co2icetotal
!            END DO
!         END DO
        write(*,*) "Not enforcing conservation of surface CO2 ice"
        write(*,*) " (should be OK when CO2 ice is a tracer)"
      end if

c-----------------------------------------------------------------------
c 6.2 Subterranean 3d variables:
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c 6.2.1 Thermal Inertia
c       Note: recall that inertiedat is a common in "comsoil.h"
c-----------------------------------------------------------------------

      ! depth-wise interpolation, if required
      if (depthinterpol.and.(.not.olddepthdef)) then
        allocate(oldval(nsoilold))
	allocate(newval(nsoilmx))
        write(*,*)'lect_start_archive: WARNING: vertical interpolation o
     &f soil thermal inertia; might be wiser to reset it.'
        write(*,*)
       
        do i=1,imold+1
         do j=1,jmold+1
	   !copy old values
	   oldval(1:nsoilold)=inertiedatold(i,j,1:nsoilold)
	   !interpolate
	   call interp_line(mlayerold,oldval,nsoilold,
     &                     mlayer,newval,nsoilmx)
           !copy interpolated values
           inertiedatoldnew(i,j,1:nsoilmx)=newval(1:nsoilmx)
	 enddo
        enddo
        ! cleanup
	deallocate(oldval)
	deallocate(newval)
      endif !of if (depthinterpol)

      if (therminertia_3D) then
        ! We have inertiedatold
       if((imold.ne.iim).or.(jmold.ne.jjm)) then
       write(*,*)'lect_start_archive: WARNING: horizontal interpolation 
     &of thermal inertia; might be better to reset it.'
       write(*,*)
       endif
       
        ! Do horizontal interpolation
	if (depthinterpol) then
	  call interp_horiz(inertiedatoldnew,inertiedatS,
     &                  imold,jmold,iim,jjm,nsoilmx,
     &                   rlonuold,rlatvold,rlonu,rlatv)
	else
          call interp_horiz(inertiedatold,inertiedatS,
     &                  imold,jmold,iim,jjm,nsoilold,
     &                   rlonuold,rlatvold,rlonu,rlatv)
        endif ! of if (depthinterpol)

      else ! no 3D thermal inertia data
       write(*,*)'lect_start_archive: using reference surface inertia'
        ! Use surface inertia (and extend it to all depths)
        do i=1,nsoilmx
         inertiedatS(1:iip1,1:jjp1,i)=surfith(1:iip1,1:jjp1)
        enddo
	! Build an old resolution surface thermal inertia
	! (will be needed for tsoil interpolation)
	call interp_horiz(surfith,surfithold,
     &                    iim,jjm,imold,jmold,1,
     &                    rlonu,rlatv,rlonuold,rlatvold)
      endif


      ! Reshape inertiedatS to scalar grid as inertiedat
      call gr_dyn_fi (nsoilmx,iim+1,jjm+1,ngrid,
     &                  inertiedatS,inertiedat)
      
c-----------------------------------------------------------------------
c 6.2.2 Soil temperature
c-----------------------------------------------------------------------
!      write(*,*) 'Soil'

      !print*,'Problem in lect_start_archive interpolating'
      !print*,'to new resolution!!'

      ! Recast temperatures along soil depth, if necessary
      if (olddepthdef) then
        allocate(oldgrid(nsoilold+1))
        allocate(oldval(nsoilold+1))
	allocate(newval(nsoilmx))
        do i=1,imold+1
	 do j=1,jmold+1

            !if(i.gt.iip1 .or. j.gt.jjp1)then
               !print*,'Problem in lect_start_archive interpolating'
               !print*,'to new resolution!!'
               !call abort
            !endif

	   ! copy values
	   oldval(1)=tsurfold(i,j)
!	   oldval(1)=tsurfS(i,j)
	   oldval(2:nsoilold+1)=tsoilold(i,j,1:nsoilold)
	   ! build vertical coordinate
	   oldgrid(1)=0. ! ground
	   oldgrid(2:nsoilold+1)=mlayerold(1:nsoilold)*
     &                (surfithold(i,j)/1.e6)
          ! Note; at this stage, we impose volcapa=1.e6 above
	  ! since volcapa isn't set in old soil definitions

	  ! interpolate
	  call interp_line(oldgrid,oldval,nsoilold+1,
     &                     mlayer,newval,nsoilmx)
	 ! copy result in tsoilold
	 tsoiloldnew(i,j,1:nsoilmx)=newval(1:nsoilmx)
	 enddo
	enddo
        ! cleanup
	deallocate(oldgrid)
	deallocate(oldval)
	deallocate(newval)

      else
       if (depthinterpol) then ! if vertical interpolation is required
        allocate(oldgrid(nsoilold+1))
        allocate(oldval(nsoilold+1))
	allocate(newval(nsoilmx))
        ! build vertical coordinate
	oldgrid(1)=0. ! ground
	oldgrid(2:nsoilold+1)=mlayerold(1:nsoilold)
        do i=1,imold+1
	 do j=1,jmold+1
	   ! copy values
	   oldval(1)=tsurfold(i,j)
!	   oldval(1)=tsurfS(i,j)
	   oldval(2:nsoilold+1)=tsoilold(i,j,1:nsoilold)
	  ! interpolate
	  call interp_line(oldgrid,oldval,nsoilold+1,
     &                     mlayer,newval,nsoilmx)
	 ! copy result in tsoilold
	 tsoiloldnew(i,j,1:nsoilmx)=newval(1:nsoilmx)
	 enddo
	enddo
!	write(*,*)'tsoiloldnew(1,1,1):',tsoiloldnew(1,1,1)
        ! cleanup
	deallocate(oldgrid)
	deallocate(oldval)
	deallocate(newval)
       
       else
        tsoiloldnew(:,:,:)=tsoilold(:,:,:)
       endif ! of if (depthinterpol)
      endif ! of if (olddepthdef)

      ! Do the horizontal interpolation
       call interp_horiz(tsoiloldnew,tsoilS,
     &                  imold,jmold,iim,jjm,nsoilmx,
     &                   rlonuold,rlatvold,rlonu,rlatv)

      ! Reshape tsoilS to scalar grid as tsoil
       call gr_dyn_fi (nsoilmx,iim+1,jjm+1,ngrid,tsoilS,tsoil)

c-----------------------------------------------------------------------
c 6.3   Slab Ocean :
c-----------------------------------------------------------------------
      call interp_horiz (tslabold,tslabs,imold,jmold,iim,jjm,noceanmx,
     &                   rlonuold,rlatvold,rlonu,rlatv)
      call gr_dyn_fi (noceanmx,iim+1,jjm+1,ngrid,tslabs,tslab)

      call interp_horiz (rnatold,rnats,imold,jmold,iim,jjm,1,
     &                   rlonuold,rlatvold,rlonu,rlatv)
      call gr_dyn_fi (1,iim+1,jjm+1,ngrid,rnats,rnat)

      call interp_horiz (pctsrf_sicold,pctsrf_sics,imold,jmold,iim,
     &                   jjm,1,rlonuold,rlatvold,rlonu,rlatv)
      call gr_dyn_fi (1,iim+1,jjm+1,ngrid,pctsrf_sics,pctsrf_sic)

      call interp_horiz (tsea_iceold,tsea_ices,imold,jmold,iim,jjm,1,
     &                   rlonuold,rlatvold,rlonu,rlatv)
      call gr_dyn_fi (1,iim+1,jjm+1,ngrid,tsea_ices,tsea_ice)

      call interp_horiz (sea_iceold,sea_ices,imold,jmold,iim,jjm,1,
     &                   rlonuold,rlatvold,rlonu,rlatv)
      call gr_dyn_fi (1,iim+1,jjm+1,ngrid,sea_ices,sea_ice)

c-----------------------------------------------------------------------
c 6.4 Variable 3d :
c-----------------------------------------------------------------------
      
c temperatures atmospheriques
      write (*,*) 'lect_start_archive: told ', told (1,jmold+1,1)  ! INFO
      call interp_vert
     &    (told,var,lmold,llm,apsold,bpsold,aps,bps,
     &     psold,(imold+1)*(jmold+1))
      write (*,*) 'lect_start_archive: var ', var (1,jmold+1,1)  ! INFO
      call interp_horiz(var,t,imold,jmold,iim,jjm,llm,
     &                   rlonuold,rlatvold,rlonu,rlatv)
      write (*,*) 'lect_start_archive: t ', t(1,jjp1,1)  ! INFO

c q2 : pbl wind variance
      write (*,*) 'lect_start_archive: q2old ', q2old (1,2,1)  ! INFO
      call interp_vert (q2old,varp1,lmold+1,llm+1,
     &     apsold,bpsold,ap,bp,psold,(imold+1)*(jmold+1))
      write (*,*) 'lect_start_archive: varp1 ', varp1 (1,2,1)  ! INFO
      call interp_horiz(varp1,q2s,imold,jmold,iim,jjm,llm+1,
     &                   rlonuold,rlatvold,rlonu,rlatv)
      write (*,*) 'lect_start_archive: q2s ', q2s (1,2,1)  ! INFO
      call gr_dyn_fi (llm+1,iim+1,jjm+1,ngrid,q2s,q2)
      write (*,*) 'lect_start_archive: q2 ', q2 (1,2)  ! INFO
c     write(47,*) 'q2',q2

c calcul des champ de vent; passage en vent covariant
      write (*,*) 'lect_start_archive: uold ', uold (1,2,1)  ! INFO
      call interp_vert
     & (uold,var,lmold,llm,apsold,bpsold,aps,bps,
     &  psold,(imold+1)*(jmold+1))
      write (*,*) 'lect_start_archive: var ', var (1,2,1)  ! INFO
      call interp_horiz(var,us,imold,jmold,iim,jjm,llm,
     &                   rlonuold,rlatvold,rlonu,rlatv)
      write (*,*) 'lect_start_archive: us ', us (1,2,1)   ! INFO

      call interp_vert
     & (vold,var,lmold,llm,
     &  apsold,bpsold,aps,bps,psold,(imold+1)*(jmold+1))
      call interp_horiz(var,vs,imold,jmold,iim,jjm,llm,
     &                   rlonuold,rlatvold,rlonu,rlatv)
      call scal_wind(us,vs,unat,vnat)
      write (*,*) 'lect_start_archive: unat ', unat (1,2,1)    ! INFO
      do l=1,llm
        do j = 1, jjp1
          do i=1,iip1
            ucov( i,j,l ) = unat( i,j,l ) * cu(i,j)
c           ucov( i,j,l ) = 0
          end do
        end do
      end do 
      write (*,*) 'lect_start_archive: ucov ', ucov (1,2,1)  ! INFO
c     write(48,*) 'ucov',ucov
      do l=1,llm
        do j = 1, jjm
          do i=1,iim
            vcov( i,j,l ) = vnat( i,j,l ) * cv(i,j)
c           vcov( i,j,l ) = 0
          end do
          vcov( iip1,j,l ) = vcov( 1,j,l )
        end do
      end do
c     write(49,*) 'ucov',vcov

      if (nqtot .gt. 0) then
c traceurs surface
      do iq = 1, nqtot
            call interp_horiz(qsurfold(1,1,iq) ,qsurfs(1,1,iq),
     &                  imold,jmold,iim,jjm,1,
     &                  rlonuold,rlatvold,rlonu,rlatv)
      enddo

      call gr_dyn_fi (nqtot,iim+1,jjm+1,ngrid,qsurfs,qsurf)

c traceurs 3D
      do  iq = 1, nqtot
            call interp_vert(qold(1,1,1,iq),var,lmold,llm,
     &        apsold,bpsold,aps,bps,psold,(imold+1)*(jmold+1))
            call interp_horiz(var,q(1,1,1,iq),imold,jmold,iim,jjm,llm,
     &                  rlonuold,rlatvold,rlonu,rlatv)
      enddo
cccccccccccccccccccccccccccccc      
c  make sure that sum of q = 1      
c dominent species is = 1 - sum(all other species)      
cccccccccccccccccccccccccccccc      
c      iqmax=1
c      
c      if (nqold.gt.10) then
c       do l=1,llm
c        do j=1,jjp1
c          do i=1,iip1
c           do iq=1,nqold
c            if (q(i,j,l,iq).gt.q(i,j,l,iqmax)) then
c              iqmax=iq
c            endif
c           enddo
c           q(i,j,l,iqmax)=1.
c           qtot(i,j,l)=0
c           do iq=1,nqold
c            if (iq.ne.iqmax) then        
c              q(i,j,l,iqmax)=q(i,j,l,iqmax)-q(i,j,l,iq)        
c            endif
c           enddo !iq
c           do iq=1,nqold
c            qtot(i,j,l)=qtot(i,j,l)+q(i,j,l,iq)
c            if (i.eq.1.and.j.eq.1.and.l.Eq.1) write(*,*)' qtot(i,j,l)',
c     $    qtot(i,j,l)
c           enddo !iq
c          enddo !i   
c         enddo !j   
c       enddo !l  
c      endif
ccccccccccccccccccccccccccccccc

c     Periodicite :
      do  iq = 1, nqtot
         do l=1, llm
            do j = 1, jjp1
               q(iip1,j,l,iq) = q(1,j,l,iq)
            end do
         end do
      enddo
      
!      call gr_dyn_fi (1,iim+1,jjm+1,ngrid,co2ices,co2ice)
! no need to transfer "co2ice" any more; it is in qsurf(igcm_co2_ice)

      endif !! if nqtot .ne. 0

c-----------------------------------------------------------------------
c   Initialisation  h:	(passage de t -> h)
c-----------------------------------------------------------------------

      DO l=1,llm
         DO j=1,jjp1
            DO i=1,iim
               h(i,j,l) = t(i,j,l)*((ps(i,j)/preff)**kappa)
            END DO
            h(iip1,j,l) =  h(1,j,l)
         END DO
      END DO


c***********************************************************************
c***********************************************************************
c     Fin subroutine lecture ini
c***********************************************************************
c***********************************************************************

      deallocate(timelist)
      deallocate(rlonuold, rlatvold)
      deallocate(rlonvold, rlatuold)
      deallocate(apsold,bpsold)
      deallocate(uold)
      deallocate(vold)
      deallocate(told)
      deallocate(psold)
      deallocate(phisold)
      deallocate(qold)
      deallocate(co2iceold)
      deallocate(tsurfold)
      deallocate(emisold)
      deallocate(q2old)
      deallocate(tsoilold)
      deallocate(tsoiloldnew)
      deallocate(inertiedatold)
      deallocate(inertiedatoldnew)
      deallocate(surfithold)
      deallocate(mlayerold)
      deallocate(qsurfold)
      deallocate(var,varp1)
      deallocate(tslabold)
      deallocate(rnatold)
      deallocate(pctsrf_sicold)
      deallocate(tsea_iceold)
      deallocate(sea_iceold)

!      write(*,*)'lect_start_archive: END'
      return
      end
