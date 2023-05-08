# Automatic Makefile

FCM_BLD_CFG = /home/palatyle/LMD_gen/trunk/LMDZ.COMMON
FCM_BLD_TARGETS = liblmdz.a start2archive_128x96x23_phystd_seq.e

export PERL5LIB := /home/palatyle/LMD_gen/trunk/LMDZ.COMMON/ioipsl/modipsl/modeles/IOIPSL/tools/FCM_V1.2/lib

export FCM_ROOTDIR = /home/palatyle/LMD_gen/trunk/LMDZ.COMMON/libo/sawtooth_128x96x23_phystd_seq.e/.config
export FCM_BINDIR = /home/palatyle/LMD_gen/trunk/LMDZ.COMMON/bin
export FCM_BLDDIR = $(FCM_ROOTDIR)/bld
export FCM_CACHEDIR = $(FCM_ROOTDIR)/.cache
export FCM_CFGDIR = $(FCM_ROOTDIR)/cfg
export FCM_DONEDIR = $(FCM_ROOTDIR)/done
export FCM_ETCDIR = $(FCM_ROOTDIR)/etc
export FCM_FLAGSDIR = $(FCM_ROOTDIR)/flags
export FCM_INCDIR = $(FCM_ROOTDIR)/inc
export FCM_LIBDIR = $(FCM_ROOTDIR)/lib
export FCM_OBJDIR = $(FCM_ROOTDIR)/obj
export FCM_PPSRCDIR = $(FCM_ROOTDIR)/ppsrc
export FCM_SRCDIR = $(FCM_ROOTDIR)/src
export FCM_TMPDIR = $(FCM_ROOTDIR)/tmp

export FCM_ROOTPATH = $(FCM_ROOTDIR)
export FCM_BINPATH = $(FCM_BINDIR)
export FCM_BLDPATH = $(FCM_BLDDIR)
export FCM_CACHEPATH = $(FCM_CACHEDIR)
export FCM_CFGPATH = $(FCM_CFGDIR)
export FCM_DONEPATH = $(FCM_DONEDIR)
export FCM_ETCPATH = $(FCM_ETCDIR)
export FCM_FLAGSPATH = $(FCM_FLAGSDIR)
export FCM_INCPATH = $(FCM_INCDIR)
export FCM_LIBPATH = $(FCM_LIBDIR)
export FCM_OBJPATH = $(FCM_OBJDIR)
export FCM_PPSRCPATH = $(FCM_PPSRCDIR)
export FCM_SRCPATH = $(FCM_SRCDIR)
export FCM_TMPPATH = $(FCM_TMPDIR)

export AR = ar
export ARFLAGS = rs
export CC = %C_COMPILER
export CC_COMPILE = -c
export CC_DEFINE = -D
export CC_INCLUDE = -I
export CC_OUTPUT = -o
export CFLAGS = %C_OPTIM %INCDIR
export CPPKEYS = CPP_PHYS %CPP_DEF
export FC = mpif90
export FC_COMPILE = -c
export FC_DEFINE = -D
export FC_INCLUDE = -I
export FC_OUTPUT = -o
export FFLAGS = -c -fdefault-real-8 -fdefault-double-8 -fcray-pointer -frecursive -ffree-line-length-0 -fno-align-commons -g  -I/home/palatyle/LMD_gen/trunk/LMDZ.COMMON/libf/aeronostd -I/apps/local/netcdf/netcdff/include
export FFLAGS__muphy = -c -I/home/palatyle/LMD_gen/trunk/LMDZ.COMMON/libf/aeronostd -I/apps/local/netcdf/netcdff/include -g 
export FPPKEYS = CPP_PHYS NC_DOUBLE
export LD = mpif90
export LDFLAGS = -Wl,-rpath=/apps/local/netcdf/netcdff/lib  -L/apps/local/netcdf/netcdff/lib -lnetcdf -lnetcdff
export LD_LIBLINK = -l
export LD_LIBSEARCH = -L
export LD_OUTPUT = -o
export SHELL = /bin/bash

export FCM_VERBOSE ?= 1

vpath % $(FCM_BINPATH)
vpath %.bash $(FCM_BINPATH)
vpath %.csh $(FCM_BINPATH)
vpath %.ksh $(FCM_BINPATH)
vpath %.pl $(FCM_BINPATH)
vpath %.pm $(FCM_BINPATH)
vpath %.pro $(FCM_BINPATH)
vpath %.py $(FCM_BINPATH)
vpath %.sh $(FCM_BINPATH)
vpath %.tcl $(FCM_BINPATH)
vpath %_128x96x23_phystd_seq.e $(FCM_BINPATH)
vpath %.mk $(FCM_BLDPATH)
vpath %.done $(FCM_DONEPATH)
vpath %.etc $(FCM_DONEPATH)
vpath %.idone $(FCM_DONEPATH)
vpath %.flags $(FCM_FLAGSPATH)
vpath %.h $(FCM_INCPATH)
vpath %.h90 $(FCM_INCPATH)
vpath %.inc $(FCM_INCPATH)
vpath %.interface $(FCM_INCPATH)
vpath %.a $(FCM_LIBPATH)
vpath %.o $(FCM_OBJPATH)
vpath %.dummy $(FCM_DONEDIR)

.PHONY : all clean

all : $(FCM_BLD_TARGETS)

clean : 
	rm -rf $(FCM_BINDIR) $(FCM_LIBDIR) $(FCM_OBJDIR) $(FCM_DONEDIR)
	rm -f lib__*.a *.o

$(FCM_DONEDIR)/FCM_CP.dummy :
	touch $@

OBJECTS__strataer =

libstrataer.a : $(OBJECTS__strataer)
	fcm_internal archive $@ $(^F)

OBJECTS__sisvat =

libsisvat.a : $(OBJECTS__sisvat)
	fcm_internal archive $@ $(^F)

OBJECTS__rrtm =

librrtm.a : $(OBJECTS__rrtm)
	fcm_internal archive $@ $(^F)

OBJECTS__phys = aerave.o aerave_new.o aeropacity.o aeroptproperties.o aerosol_mod.o spectral_albedo_calc.o albedo_snow_calc.o ave_stelspec.o bilinear.o bilinearbig.o blackl.o calc_cpp3d.o calc_cpp_mugaz.o calc_rayleigh.o call_rings.o callcorrk_mod.o callkeys_mod.o callsedim.o comcstfi_mod.o comdiurn_h.o comg1d_mod.o comgeomfi_h.o comm_wrf.o comsaison_h.o comsoil_h.o conc_mod.o condense_co2.o convadj.o cp_neutral.o cpdet_phy_mod.o cvmgp.o cvmgt.o datafile_mod.o def_var.o dimphy.o dsolver.o dtridgl.o evap.o forcewcfn.o gases_h.o gfluxi.o gfluxv.o hydrol.o iniaerosol.o inifis_mod.o iniorbit.o inistats.o initracer.o iniwrite.o iniwrite_specir.o iniwrite_specvi.o iniwritesoil.o interp_line.o interpolateh2h2.o interpolateh2he.o interpolateh2ocont_ckd.o interpolateh2ocont_ppc.o interpolaten2h2.o interpolaten2n2.o iostart.o lagrange.o largescale.o mass_redistribution.o mkstats.o moistadj.o mucorr.o mufract.o multipl.o newsedim.o newtrelax.o ocean_slab_mod.o optci_mod.o optcv_mod.o orbite.o params_h.o phyetat0_mod.o phyredem.o phys_state_var_mod.o physiq_mod.o planete_mod.o planetwide_mod.o profile.o radcommon_h.o radii_mod.o radinc_h.o rain.o rings.o scopyi.o setspi.o setspv.o sfluxi.o sfluxv.o slab_ice_h.o soil.o soil_settings.o statto_mod.o stelang.o stellarlong.o stokes.o su_gases.o suaer_corrk.o sugas_corrk.o surf_heat_transp_mod.o surface_nature.o surfdat_h.o surfini.o tabfi_mod.o thermcell_closure.o thermcell_dq.o thermcell_dv2.o thermcell_env.o thermcell_flux.o thermcell_height.o thermcell_main.o thermcell_mod.o thermcell_plume.o time_phylmdz_mod.o totalcloudfrac.o tpindex.o tracer_h.o turb_mod.o turbdiff.o vdif_cd.o vdif_kc.o vdifc.o vlz_fi.o volcano.o watercommon_h.o writediagfi.o writediagsoil.o writediagspecir.o writediagspecvi.o writeg1d.o wstats.o xios_output_mod.o

libphys.a : $(OBJECTS__phys)
	fcm_internal archive $@ $(^F)

OBJECTS__phy_common = abort_physic.o geometry_mod.o init_print_control_mod.o ioipsl_getin_p_mod.o mod_grid_phy_lmdz.o mod_phys_lmdz_mpi_data.o mod_phys_lmdz_mpi_transfert.o mod_phys_lmdz_omp_data.o mod_phys_lmdz_omp_transfert.o mod_phys_lmdz_para.o mod_phys_lmdz_transfert_para.o physics_distribution_mod.o print_control_mod.o regular_lonlat_mod.o vertical_layers_mod.o write_field_phy.o xios_writefield_mod.o

libphy_common.a : $(OBJECTS__phy_common)
	fcm_internal archive $@ $(^F)

OBJECTS__muphy =

libmuphy.a : $(OBJECTS__muphy)
	fcm_internal archive $@ $(^F)

OBJECTS__misc = arth_m.o assert_eq_m.o assert_m.o cbrt.o chfev.o coefpoly_m.o scopy.o fdump.o formcoord.o handle_err_m.o i1mach.o interpolation.o ioipsl_errioipsl.o ioipsl_getincom.o ioipsl_stringop.o ismax.o ismin.o j4save.o juldate.o minmax.o minmax2.o netcdf95.o new_unit_m.o nf95_def_var_m.o nf95_get_att_m.o nf95_get_var_m.o nf95_gw_var_m.o nf95_put_att_m.o nf95_put_var_m.o nrtype.o pchdf.o pchfe.o pchfe_95_m.o pchsp.o pchsp_95_m.o pres2lev_mod.o q_sat.o ran1.o regr1_conserv_m.o regr1_lint_m.o regr3_lint_m.o simple.o slopes_m.o sort.o sort_mod.o vampir.o write_field.o xercnt.o xerhlt.o xermsg.o xerprn.o xersve.o xgetua.o

libmisc.a : $(OBJECTS__misc)
	fcm_internal archive $@ $(^F)

OBJECTS__grid =

libgrid.a : $(OBJECTS__grid)
	fcm_internal archive $@ $(^F)

OBJECTS__filtrez = acc.o eigen.o eigen_sort.o filtreg.o filtreg_mod.o inifgn.o jacobi.o mkl_dft_type.o mkl_dfti.o mod_fft.o mod_fft_fftw.o mod_fft_mathkeisan.o mod_fft_mkl.o mod_fft_wrapper.o mod_filtre_fft.o mod_filtre_fft_loc.o timer_filtre.o

libfiltrez.a : $(OBJECTS__filtrez)
	fcm_internal archive $@ $(^F)

OBJECTS__dyn_phys_sub = caldyn0.o callphysiq_mod.o datareadnc.o defrun_new.o exner_hyb.o grid_noro1.o ini_archive.o inichim_newstart.o iniphysiq_mod.o interp_vert.o lect_start_archive.o scal_wind.o wind_scal.o write_archive.o

libdyn_phys_sub.a : $(OBJECTS__dyn_phys_sub)
	fcm_internal archive $@ $(^F)

OBJECTS__dyn_phys = calfis.o calfis_p.o gr_dyn_fi.o gr_dyn_fi_p.o gr_fi_dyn.o gr_fi_dyn_p.o grid_atob_m.o inigeomphy_mod.o mod_interface_dyn_phys.o

libdyn_phys.a : $(OBJECTS__dyn_phys)
	fcm_internal archive $@ $(^F)

OBJECTS__dyn3d_common = adaptdt.o advn.o advx.o advxp.o advy.o advyp.o advz.o advzp.o bernoui.o com_io_dyn_mod.o comconst_mod.o comuforc_h.o comvert_mod.o conf_dat_m.o conf_planete.o control_mod.o convflu.o convmas.o coordij.o covcont.o covnat.o cpdet_mod.o defrun.o diagedyn.o disvert.o disvert_noterre.o diverg.o diverg_gam.o divergf.o divergst.o divgrad.o divgrad2.o dump2d.o dynetat0.o dynredem0.o dynredem_mod.o ener_mod.o enercin.o exner_hyb_m.o exner_milieu_m.o extrapol.o flumass.o fxhyp_m.o fxy.o fxyhyper.o fxysinus.o fyhyp_m.o geopot.o gr_ecrit_fi.o gr_int_dyn.o gr_u_scal.o gr_v_scal.o grad.o gradiv.o gradiv2.o grilles_gcm_netcdf_sub.o heavyside.o infotrac.o iniacademic.o iniconst.o inidissip.o inigeom.o inigrads.o initdynav.o initfluxsto.o inithist.o initial0.o inter_barxy_m.o interpost.o interpre.o invert_lat.o invert_zoom_x_m.o iso_verif_nonan_nostop.o laplacien.o laplacien_gam.o laplacien_rot.o laplacien_rotgam.o limx.o limy.o limz.o massbar.o massbarxy.o massdair.o misc_mod.o nxgrad.o nxgrad_gam.o nxgradst.o nxgraro2.o nxgrarot.o pbar.o pentes_ini.o ppm3d.o prather.o pression.o principal_cshift_m.o psextbar.o relaxetat0.o rotat.o rotat_nfil.o rotatf.o rotatst.o serre_mod.o sortvarc.o sortvarc0.o test_period.o tidal_forces.o tourpot.o traceurpole.o ugeostr.o vitvert.o writedynav.o writehist.o

libdyn3d_common.a : $(OBJECTS__dyn3d_common)
	fcm_internal archive $@ $(^F)

OBJECTS__dyn = abort_gcm.o addfi.o advect.o advtrac.o bilan_dyn.o caladvtrac.o caldyn.o check_isotopes_seq.o conf_gcm.o dissip.o dteta1.o dudv1.o dudv2.o fluxstokenc.o friction.o getparam.o groupe.o groupeun.o guide_mod.o iniinterp_horiz.o integrd.o interp_horiz.o leapfrog.o leapfrog_nogcm.o logic_mod.o mod_const_mpi.o moyzon_mod.o planetary_operations.o qminimum.o sponge_mod.o sw_case_williamson91_6.o temps_mod.o tetaleveli1j.o tetaleveli1j1.o top_bound.o vlsplt.o vlspltqs.o wrgrads.o

libdyn.a : $(OBJECTS__dyn)
	fcm_internal archive $@ $(^F)

OBJECTS__dust =

libdust.a : $(OBJECTS__dust)
	fcm_internal archive $@ $(^F)

OBJECTS__cosp =

libcosp.a : $(OBJECTS__cosp)
	fcm_internal archive $@ $(^F)

OBJECTS__cloud =

libcloud.a : $(OBJECTS__cloud)
	fcm_internal archive $@ $(^F)

OBJECTS__chem =

libchem.a : $(OBJECTS__chem)
	fcm_internal archive $@ $(^F)

OBJECTS__aerono = calchim_asis.o concentrations.o intrplf.o inv.o photochemistry_asis.o photolysis_asis.o read_phototable.o types_asis.o

libaerono.a : $(OBJECTS__aerono)
	fcm_internal archive $@ $(^F)

export OBJECTS = $(OBJECTS__aerono) $(OBJECTS__chem) $(OBJECTS__cloud) $(OBJECTS__cosp) $(OBJECTS__dust) $(OBJECTS__dyn) $(OBJECTS__dyn3d_common) $(OBJECTS__dyn_phys) $(OBJECTS__dyn_phys_sub) $(OBJECTS__filtrez) $(OBJECTS__grid) $(OBJECTS__misc) $(OBJECTS__muphy) $(OBJECTS__phy_common) $(OBJECTS__phys) $(OBJECTS__rrtm) $(OBJECTS__sisvat) $(OBJECTS__strataer)

liblmdz.a : $(OBJECTS)
	fcm_internal archive $@ $(^F)

CC.flags :
	touch $(FCM_FLAGSDIR)/$@

CFLAGS.flags : CC.flags
	touch $(FCM_FLAGSDIR)/$@

CFLAGS__aerono.flags CFLAGS__chem.flags CFLAGS__cloud.flags CFLAGS__cosp.flags CFLAGS__dust.flags CFLAGS__dyn.flags CFLAGS__dyn3d_common.flags CFLAGS__dyn_phys.flags CFLAGS__dyn_phys_sub.flags CFLAGS__filtrez.flags CFLAGS__grid.flags CFLAGS__misc.flags CFLAGS__muphy.flags CFLAGS__phy_common.flags CFLAGS__phys.flags CFLAGS__rrtm.flags CFLAGS__sisvat.flags CFLAGS__strataer.flags : CFLAGS.flags
	touch $(FCM_FLAGSDIR)/$@

CPPKEYS.flags :
	touch $(FCM_FLAGSDIR)/$@

CPPKEYS__aerono.flags CPPKEYS__chem.flags CPPKEYS__cloud.flags CPPKEYS__cosp.flags CPPKEYS__dust.flags CPPKEYS__dyn.flags CPPKEYS__dyn3d_common.flags CPPKEYS__dyn_phys.flags CPPKEYS__dyn_phys_sub.flags CPPKEYS__filtrez.flags CPPKEYS__grid.flags CPPKEYS__misc.flags CPPKEYS__muphy.flags CPPKEYS__phy_common.flags CPPKEYS__phys.flags CPPKEYS__rrtm.flags CPPKEYS__sisvat.flags CPPKEYS__strataer.flags : CPPKEYS.flags
	touch $(FCM_FLAGSDIR)/$@

FC.flags :
	touch $(FCM_FLAGSDIR)/$@

FFLAGS.flags : FC.flags
	touch $(FCM_FLAGSDIR)/$@

FFLAGS__aerono.flags FFLAGS__chem.flags FFLAGS__cloud.flags FFLAGS__cosp.flags FFLAGS__dust.flags FFLAGS__dyn.flags FFLAGS__dyn3d_common.flags FFLAGS__dyn_phys.flags FFLAGS__dyn_phys_sub.flags FFLAGS__filtrez.flags FFLAGS__grid.flags FFLAGS__misc.flags FFLAGS__muphy.flags FFLAGS__phy_common.flags FFLAGS__phys.flags FFLAGS__rrtm.flags FFLAGS__sisvat.flags FFLAGS__strataer.flags : FFLAGS.flags
	touch $(FCM_FLAGSDIR)/$@

FPPKEYS.flags :
	touch $(FCM_FLAGSDIR)/$@

FPPKEYS__aerono.flags FPPKEYS__chem.flags FPPKEYS__cloud.flags FPPKEYS__cosp.flags FPPKEYS__dust.flags FPPKEYS__dyn.flags FPPKEYS__dyn3d_common.flags FPPKEYS__dyn_phys.flags FPPKEYS__dyn_phys_sub.flags FPPKEYS__filtrez.flags FPPKEYS__grid.flags FPPKEYS__misc.flags FPPKEYS__muphy.flags FPPKEYS__phy_common.flags FPPKEYS__phys.flags FPPKEYS__rrtm.flags FPPKEYS__sisvat.flags FPPKEYS__strataer.flags : FPPKEYS.flags
	touch $(FCM_FLAGSDIR)/$@

LD.flags :
	touch $(FCM_FLAGSDIR)/$@

LD__aerono.flags LD__chem.flags LD__cloud.flags LD__cosp.flags LD__dust.flags LD__dyn.flags LD__dyn3d_common.flags LD__dyn_phys.flags LD__dyn_phys_sub.flags LD__filtrez.flags LD__grid.flags LD__misc.flags LD__muphy.flags LD__phy_common.flags LD__phys.flags LD__rrtm.flags LD__sisvat.flags LD__strataer.flags : LD.flags
	touch $(FCM_FLAGSDIR)/$@

LDFLAGS.flags :
	touch $(FCM_FLAGSDIR)/$@

LDFLAGS__aerono.flags LDFLAGS__chem.flags LDFLAGS__cloud.flags LDFLAGS__cosp.flags LDFLAGS__dust.flags LDFLAGS__dyn.flags LDFLAGS__dyn3d_common.flags LDFLAGS__dyn_phys.flags LDFLAGS__dyn_phys_sub.flags LDFLAGS__filtrez.flags LDFLAGS__grid.flags LDFLAGS__misc.flags LDFLAGS__muphy.flags LDFLAGS__phy_common.flags LDFLAGS__phys.flags LDFLAGS__rrtm.flags LDFLAGS__sisvat.flags LDFLAGS__strataer.flags : LDFLAGS.flags
	touch $(FCM_FLAGSDIR)/$@

include $(FCM_BLDDIR)/aerono.mk
include $(FCM_BLDDIR)/dyn.mk
include $(FCM_BLDDIR)/dyn3d_common.mk
include $(FCM_BLDDIR)/dyn_phys.mk
include $(FCM_BLDDIR)/dyn_phys_sub.mk
include $(FCM_BLDDIR)/filtrez.mk
include $(FCM_BLDDIR)/grid.mk
include $(FCM_BLDDIR)/misc.mk
include $(FCM_BLDDIR)/phy_common.mk
include $(FCM_BLDDIR)/phys.mk

# EOF
