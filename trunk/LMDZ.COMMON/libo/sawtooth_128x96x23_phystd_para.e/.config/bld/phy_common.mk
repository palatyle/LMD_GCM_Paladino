# Automatic Make rule for phy_common

SRCDIR0__phy_common = /home/palatyle/LMD_gen/trunk/LMDZ.COMMON/libf/phy_common

PPSRCDIR0__phy_common = /home/palatyle/LMD_gen/trunk/LMDZ.COMMON/libo/sawtooth_128x96x23_phystd_para.e/.config/ppsrc/phy_common

FFLAGS__phy_common__physics_distribution_mod.flags: \
          FFLAGS__phy_common.flags
	touch $(FCM_FLAGSDIR)/$@

physics_distribution_mod.done: \
          physics_distribution_mod.o \
          dimphy.done \
          mod_grid_phy_lmdz.done \
          mod_phys_lmdz_para.done
	touch $(FCM_DONEDIR)/$@

physics_distribution_mod.o: \
          $(PPSRCDIR0__phy_common)/physics_distribution_mod.f90 \
          FFLAGS__phy_common__physics_distribution_mod.flags \
          dimphy.o \
          mod_grid_phy_lmdz.o \
          mod_phys_lmdz_para.o
	fcm_internal compile:F phy_common $< $@

FFLAGS__phy_common__write_field_phy.flags: \
          FFLAGS__phy_common.flags
	touch $(FCM_FLAGSDIR)/$@

write_field_phy.done: \
          write_field_phy.o \
          write_field.done \
          mod_grid_phy_lmdz.done \
          mod_phys_lmdz_para.done
	touch $(FCM_DONEDIR)/$@

write_field_phy.o: \
          $(PPSRCDIR0__phy_common)/write_field_phy.f90 \
          FFLAGS__phy_common__write_field_phy.flags \
          write_field.o \
          mod_grid_phy_lmdz.o \
          mod_phys_lmdz_para.o
	fcm_internal compile:F phy_common $< $@

FFLAGS__phy_common__geometry_mod.flags: \
          FFLAGS__phy_common.flags
	touch $(FCM_FLAGSDIR)/$@

geometry_mod.done: \
          geometry_mod.o \
          mod_grid_phy_lmdz.done \
          mod_phys_lmdz_mpi_data.done \
          nrtype.done \
          print_control_mod.done
	touch $(FCM_DONEDIR)/$@

geometry_mod.o: \
          $(PPSRCDIR0__phy_common)/geometry_mod.f90 \
          FFLAGS__phy_common__geometry_mod.flags \
          mod_grid_phy_lmdz.o \
          mod_phys_lmdz_mpi_data.o \
          nrtype.o \
          print_control_mod.o
	fcm_internal compile:F phy_common $< $@

FFLAGS__phy_common__abort_physic.flags: \
          FFLAGS__phy_common.flags
	touch $(FCM_FLAGSDIR)/$@

abort_physic.done: \
          abort_physic.o \
          ioipsl_getincom.done \
          mod_phys_lmdz_para.done \
          print_control_mod.done
	touch $(FCM_DONEDIR)/$@

abort_physic.o: \
          $(PPSRCDIR0__phy_common)/abort_physic.f90 \
          FFLAGS__phy_common__abort_physic.flags \
          ioipsl_getincom.o \
          mod_phys_lmdz_para.o \
          print_control_mod.o
	fcm_internal compile:F phy_common $< $@

FFLAGS__phy_common__print_control_mod.flags: \
          FFLAGS__phy_common.flags
	touch $(FCM_FLAGSDIR)/$@

print_control_mod.done: \
          print_control_mod.o
	touch $(FCM_DONEDIR)/$@

print_control_mod.o: \
          $(PPSRCDIR0__phy_common)/print_control_mod.f90 \
          FFLAGS__phy_common__print_control_mod.flags
	fcm_internal compile:F phy_common $< $@

FFLAGS__phy_common__regular_lonlat_mod.flags: \
          FFLAGS__phy_common.flags
	touch $(FCM_FLAGSDIR)/$@

regular_lonlat_mod.done: \
          regular_lonlat_mod.o
	touch $(FCM_DONEDIR)/$@

regular_lonlat_mod.o: \
          $(PPSRCDIR0__phy_common)/regular_lonlat_mod.f90 \
          FFLAGS__phy_common__regular_lonlat_mod.flags
	fcm_internal compile:F phy_common $< $@

FFLAGS__phy_common__ioipsl_getin_p_mod.flags: \
          FFLAGS__phy_common.flags
	touch $(FCM_FLAGSDIR)/$@

ioipsl_getin_p_mod.done: \
          ioipsl_getin_p_mod.o \
          ioipsl_getincom.done \
          mod_phys_lmdz_mpi_data.done \
          mod_phys_lmdz_omp_data.done \
          mod_phys_lmdz_transfert_para.done
	touch $(FCM_DONEDIR)/$@

ioipsl_getin_p_mod.o: \
          $(PPSRCDIR0__phy_common)/ioipsl_getin_p_mod.f90 \
          FFLAGS__phy_common__ioipsl_getin_p_mod.flags \
          ioipsl_getincom.o \
          mod_phys_lmdz_mpi_data.o \
          mod_phys_lmdz_omp_data.o \
          mod_phys_lmdz_transfert_para.o
	fcm_internal compile:F phy_common $< $@

FFLAGS__phy_common__xios_writefield.flags: \
          FFLAGS__phy_common.flags
	touch $(FCM_FLAGSDIR)/$@

xios_writefield_mod.done: \
          xios_writefield_mod.o
	touch $(FCM_DONEDIR)/$@

xios_writefield_mod.o: \
          $(PPSRCDIR0__phy_common)/xios_writefield.f90 \
          FFLAGS__phy_common__xios_writefield.flags
	fcm_internal compile:F phy_common $< $@

FFLAGS__phy_common__mod_phys_lmdz_mpi_data.flags: \
          FFLAGS__phy_common.flags
	touch $(FCM_FLAGSDIR)/$@

mod_phys_lmdz_mpi_data.done: \
          mod_phys_lmdz_mpi_data.o \
          print_control_mod.done
	touch $(FCM_DONEDIR)/$@

mod_phys_lmdz_mpi_data.o: \
          $(PPSRCDIR0__phy_common)/mod_phys_lmdz_mpi_data.f90 \
          FFLAGS__phy_common__mod_phys_lmdz_mpi_data.flags \
          print_control_mod.o
	fcm_internal compile:F phy_common $< $@

FFLAGS__phy_common__mod_phys_lmdz_transfert_para.flags: \
          FFLAGS__phy_common.flags
	touch $(FCM_FLAGSDIR)/$@

mod_phys_lmdz_transfert_para.done: \
          mod_phys_lmdz_transfert_para.o \
          mod_phys_lmdz_mpi_data.done \
          mod_phys_lmdz_mpi_transfert.done \
          mod_phys_lmdz_omp_transfert.done
	touch $(FCM_DONEDIR)/$@

mod_phys_lmdz_transfert_para.o: \
          $(PPSRCDIR0__phy_common)/mod_phys_lmdz_transfert_para.f90 \
          FFLAGS__phy_common__mod_phys_lmdz_transfert_para.flags \
          mod_phys_lmdz_mpi_data.o \
          mod_phys_lmdz_mpi_transfert.o \
          mod_phys_lmdz_omp_transfert.o
	fcm_internal compile:F phy_common $< $@

FFLAGS__phy_common__mod_grid_phy_lmdz.flags: \
          FFLAGS__phy_common.flags
	touch $(FCM_FLAGSDIR)/$@

mod_grid_phy_lmdz.done: \
          mod_grid_phy_lmdz.o
	touch $(FCM_DONEDIR)/$@

mod_grid_phy_lmdz.o: \
          $(PPSRCDIR0__phy_common)/mod_grid_phy_lmdz.f90 \
          FFLAGS__phy_common__mod_grid_phy_lmdz.flags
	fcm_internal compile:F phy_common $< $@

FFLAGS__phy_common__init_print_control_mod.flags: \
          FFLAGS__phy_common.flags
	touch $(FCM_FLAGSDIR)/$@

init_print_control_mod.done: \
          init_print_control_mod.o \
          ioipsl_getin_p_mod.done \
          mod_phys_lmdz_para.done \
          print_control_mod.done
	touch $(FCM_DONEDIR)/$@

init_print_control_mod.o: \
          $(PPSRCDIR0__phy_common)/init_print_control_mod.f90 \
          FFLAGS__phy_common__init_print_control_mod.flags \
          ioipsl_getin_p_mod.o \
          mod_phys_lmdz_para.o \
          print_control_mod.o
	fcm_internal compile:F phy_common $< $@

FFLAGS__phy_common__mod_phys_lmdz_omp_transfert.flags: \
          FFLAGS__phy_common.flags
	touch $(FCM_FLAGSDIR)/$@

mod_phys_lmdz_omp_transfert.done: \
          mod_phys_lmdz_omp_transfert.o \
          mod_phys_lmdz_mpi_data.done \
          mod_phys_lmdz_omp_data.done
	touch $(FCM_DONEDIR)/$@

mod_phys_lmdz_omp_transfert.o: \
          $(PPSRCDIR0__phy_common)/mod_phys_lmdz_omp_transfert.f90 \
          FFLAGS__phy_common__mod_phys_lmdz_omp_transfert.flags \
          mod_phys_lmdz_mpi_data.o \
          mod_phys_lmdz_omp_data.o
	fcm_internal compile:F phy_common $< $@

FFLAGS__phy_common__vertical_layers_mod.flags: \
          FFLAGS__phy_common.flags
	touch $(FCM_FLAGSDIR)/$@

vertical_layers_mod.done: \
          vertical_layers_mod.o
	touch $(FCM_DONEDIR)/$@

vertical_layers_mod.o: \
          $(PPSRCDIR0__phy_common)/vertical_layers_mod.f90 \
          FFLAGS__phy_common__vertical_layers_mod.flags
	fcm_internal compile:F phy_common $< $@

FFLAGS__phy_common__mod_phys_lmdz_para.flags: \
          FFLAGS__phy_common.flags
	touch $(FCM_FLAGSDIR)/$@

mod_phys_lmdz_para.done: \
          mod_phys_lmdz_para.o \
          mod_grid_phy_lmdz.done \
          mod_phys_lmdz_mpi_data.done \
          mod_phys_lmdz_omp_data.done \
          mod_phys_lmdz_transfert_para.done \
          print_control_mod.done
	touch $(FCM_DONEDIR)/$@

mod_phys_lmdz_para.o: \
          $(PPSRCDIR0__phy_common)/mod_phys_lmdz_para.f90 \
          FFLAGS__phy_common__mod_phys_lmdz_para.flags \
          mod_grid_phy_lmdz.o \
          mod_phys_lmdz_mpi_data.o \
          mod_phys_lmdz_omp_data.o \
          mod_phys_lmdz_transfert_para.o \
          print_control_mod.o
	fcm_internal compile:F phy_common $< $@

FFLAGS__phy_common__mod_phys_lmdz_omp_data.flags: \
          FFLAGS__phy_common.flags
	touch $(FCM_FLAGSDIR)/$@

mod_phys_lmdz_omp_data.done: \
          mod_phys_lmdz_omp_data.o \
          dimphy.done \
          mod_phys_lmdz_mpi_data.done \
          print_control_mod.done
	touch $(FCM_DONEDIR)/$@

mod_phys_lmdz_omp_data.o: \
          $(PPSRCDIR0__phy_common)/mod_phys_lmdz_omp_data.f90 \
          FFLAGS__phy_common__mod_phys_lmdz_omp_data.flags \
          dimphy.o \
          mod_phys_lmdz_mpi_data.o \
          print_control_mod.o
	fcm_internal compile:F phy_common $< $@

FFLAGS__phy_common__mod_phys_lmdz_mpi_transfert.flags: \
          FFLAGS__phy_common.flags
	touch $(FCM_FLAGSDIR)/$@

mod_phys_lmdz_mpi_transfert.done: \
          mod_phys_lmdz_mpi_transfert.o \
          mod_grid_phy_lmdz.done \
          mod_phys_lmdz_mpi_data.done
	touch $(FCM_DONEDIR)/$@

mod_phys_lmdz_mpi_transfert.o: \
          $(PPSRCDIR0__phy_common)/mod_phys_lmdz_mpi_transfert.f90 \
          FFLAGS__phy_common__mod_phys_lmdz_mpi_transfert.flags \
          mod_grid_phy_lmdz.o \
          mod_phys_lmdz_mpi_data.o
	fcm_internal compile:F phy_common $< $@

