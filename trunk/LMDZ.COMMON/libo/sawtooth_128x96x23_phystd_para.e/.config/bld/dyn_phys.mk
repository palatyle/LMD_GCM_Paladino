# Automatic Make rule for dyn_phys

SRCDIR0__dyn_phys = /home/palatyle/LMD_gen/trunk/LMDZ.COMMON/libf/dynphy_lonlat

PPSRCDIR0__dyn_phys = /home/palatyle/LMD_gen/trunk/LMDZ.COMMON/libo/sawtooth_128x96x23_phystd_para.e/.config/ppsrc/dyn_phys

FFLAGS__dyn_phys__gr_dyn_fi.flags: \
          FFLAGS__dyn_phys.flags
	touch $(FCM_FLAGSDIR)/$@

gr_dyn_fi.done: \
          gr_dyn_fi.o
	touch $(FCM_DONEDIR)/$@

gr_dyn_fi.o: \
          $(PPSRCDIR0__dyn_phys)/gr_dyn_fi.f \
          FFLAGS__dyn_phys__gr_dyn_fi.flags
	fcm_internal compile:F dyn_phys $< $@

FFLAGS__dyn_phys__gr_fi_dyn_p.flags: \
          FFLAGS__dyn_phys.flags
	touch $(FCM_FLAGSDIR)/$@

gr_fi_dyn_p.done: \
          gr_fi_dyn_p.o \
          dimphy.done \
          mod_interface_dyn_phys.done \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

gr_fi_dyn_p.o: \
          $(PPSRCDIR0__dyn_phys)/gr_fi_dyn_p.f \
          FFLAGS__dyn_phys__gr_fi_dyn_p.flags \
          dimphy.o \
          mod_interface_dyn_phys.o \
          parallel_lmdz.o
	fcm_internal compile:F dyn_phys $< $@

FFLAGS__dyn_phys__inigeomphy_mod.flags: \
          FFLAGS__dyn_phys.flags
	touch $(FCM_FLAGSDIR)/$@

inigeomphy_mod.done: \
          inigeomphy_mod.o \
          comvert_mod.done \
          geometry_mod.done \
          iniprint.h.idone \
          mod_grid_phy_lmdz.done \
          mod_interface_dyn_phys.done \
          mod_phys_lmdz_para.done \
          nrtype.done \
          physics_distribution_mod.done \
          regular_lonlat_mod.done \
          vertical_layers_mod.done
	touch $(FCM_DONEDIR)/$@

inigeomphy_mod.o: \
          $(PPSRCDIR0__dyn_phys)/inigeomphy_mod.f90 \
          FFLAGS__dyn_phys__inigeomphy_mod.flags \
          comvert_mod.o \
          geometry_mod.o \
          iniprint.h \
          mod_grid_phy_lmdz.o \
          mod_interface_dyn_phys.o \
          mod_phys_lmdz_para.o \
          nrtype.o \
          physics_distribution_mod.o \
          regular_lonlat_mod.o \
          vertical_layers_mod.o
	fcm_internal compile:F dyn_phys $< $@

FFLAGS__dyn_phys__mod_interface_dyn_phys.flags: \
          FFLAGS__dyn_phys.flags
	touch $(FCM_FLAGSDIR)/$@

mod_interface_dyn_phys.done: \
          mod_interface_dyn_phys.o \
          dimensions.h.idone \
          mod_phys_lmdz_mpi_data.done
	touch $(FCM_DONEDIR)/$@

mod_interface_dyn_phys.o: \
          $(PPSRCDIR0__dyn_phys)/mod_interface_dyn_phys.f90 \
          FFLAGS__dyn_phys__mod_interface_dyn_phys.flags \
          dimensions.h \
          mod_phys_lmdz_mpi_data.o
	fcm_internal compile:F dyn_phys $< $@

FFLAGS__dyn_phys__gr_dyn_fi_p.flags: \
          FFLAGS__dyn_phys.flags
	touch $(FCM_FLAGSDIR)/$@

gr_dyn_fi_p.done: \
          gr_dyn_fi_p.o \
          dimphy.done \
          mod_interface_dyn_phys.done \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

gr_dyn_fi_p.o: \
          $(PPSRCDIR0__dyn_phys)/gr_dyn_fi_p.f \
          FFLAGS__dyn_phys__gr_dyn_fi_p.flags \
          dimphy.o \
          mod_interface_dyn_phys.o \
          parallel_lmdz.o
	fcm_internal compile:F dyn_phys $< $@

FFLAGS__dyn_phys__grid_atob_m.flags: \
          FFLAGS__dyn_phys.flags
	touch $(FCM_FLAGSDIR)/$@

grid_atob_m.done: \
          grid_atob_m.o \
          assert_eq_m.done
	touch $(FCM_DONEDIR)/$@

grid_atob_m.o: \
          $(PPSRCDIR0__dyn_phys)/grid_atob_m.f90 \
          FFLAGS__dyn_phys__grid_atob_m.flags \
          assert_eq_m.o
	fcm_internal compile:F dyn_phys $< $@

FFLAGS__dyn_phys__calfis.flags: \
          FFLAGS__dyn_phys.flags
	touch $(FCM_FLAGSDIR)/$@

calfis.done: \
          calfis.o \
          callphysiq_mod.done \
          comconst_mod.done \
          comgeom2.h.idone \
          comvert_mod.done \
          control_mod.done \
          cpdet_mod.done \
          dimensions.h.idone \
          infotrac.done \
          iniprint.h.idone \
          logic_mod.done \
          moyzon_mod.done \
          paramet.h.idone \
          write_field.done
	touch $(FCM_DONEDIR)/$@

calfis.o: \
          $(PPSRCDIR0__dyn_phys)/calfis.f \
          FFLAGS__dyn_phys__calfis.flags \
          callphysiq_mod.o \
          comconst_mod.o \
          comgeom2.h \
          comvert_mod.o \
          control_mod.o \
          cpdet_mod.o \
          dimensions.h \
          infotrac.o \
          iniprint.h \
          logic_mod.o \
          moyzon_mod.o \
          paramet.h \
          write_field.o
	fcm_internal compile:F dyn_phys $< $@

FFLAGS__dyn_phys__calfis_p.flags: \
          FFLAGS__dyn_phys.flags
	touch $(FCM_FLAGSDIR)/$@

calfis_p.done: \
          calfis_p.o \
          times.done \
          write_field.done \
          write_field_p.done \
          callphysiq_mod.done \
          comconst_mod.done \
          comgeom2.h.idone \
          comvert_mod.done \
          control_mod.done \
          cpdet_mod.done \
          dimensions.h.idone \
          dimphy.done \
          infotrac.done \
          iniprint.h.idone \
          logic_mod.done \
          mod_const_mpi.done \
          mod_interface_dyn_phys.done \
          mod_phys_lmdz_mpi_data.done \
          mod_phys_lmdz_omp_data.done \
          moyzon_mod.done \
          parallel_lmdz.done \
          paramet.h.idone
	touch $(FCM_DONEDIR)/$@

calfis_p.o: \
          $(PPSRCDIR0__dyn_phys)/calfis_p.f \
          FFLAGS__dyn_phys__calfis_p.flags \
          times.o \
          write_field.o \
          write_field_p.o \
          callphysiq_mod.o \
          comconst_mod.o \
          comgeom2.h \
          comvert_mod.o \
          control_mod.o \
          cpdet_mod.o \
          dimensions.h \
          dimphy.o \
          infotrac.o \
          iniprint.h \
          logic_mod.o \
          mod_const_mpi.o \
          mod_interface_dyn_phys.o \
          mod_phys_lmdz_mpi_data.o \
          mod_phys_lmdz_omp_data.o \
          moyzon_mod.o \
          parallel_lmdz.o \
          paramet.h
	fcm_internal compile:F dyn_phys $< $@

FFLAGS__dyn_phys__gr_fi_dyn.flags: \
          FFLAGS__dyn_phys.flags
	touch $(FCM_FLAGSDIR)/$@

gr_fi_dyn.done: \
          gr_fi_dyn.o
	touch $(FCM_DONEDIR)/$@

gr_fi_dyn.o: \
          $(PPSRCDIR0__dyn_phys)/gr_fi_dyn.f \
          FFLAGS__dyn_phys__gr_fi_dyn.flags
	fcm_internal compile:F dyn_phys $< $@

