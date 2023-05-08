# Automatic Make rule for dyn_phys_sub

SRCDIR0__dyn_phys_sub = /home/palatyle/LMD_gen/trunk/LMDZ.COMMON/libf/dynphy_lonlat/phystd

PPSRCDIR0__dyn_phys_sub = /home/palatyle/LMD_gen/trunk/LMDZ.COMMON/libo/sawtooth_128x96x23_phystd_seq.e/.config/ppsrc/dyn_phys_sub

FFLAGS__dyn_phys_sub__start2archive.flags: \
          FFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

LDFLAGS__dyn_phys_sub__start2archive.flags: \
          LDFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

LD__dyn_phys_sub__start2archive.flags: \
          LD__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

start2archive.o: \
          $(PPSRCDIR0__dyn_phys_sub)/start2archive.f \
          FFLAGS__dyn_phys_sub__start2archive.flags \
          callkeys_mod.o \
          comconst_mod.o \
          comdissip.h \
          comgeom.h \
          comsoil_h.o \
          comvert_mod.o \
          control_mod.o \
          dimensions.h \
          filtreg_mod.o \
          infotrac.o \
          iniphysiq_mod.o \
          ioipsl_getincom.o \
          mod_const_mpi.o \
          paramet.h \
          phyetat0_mod.o \
          planete_mod.o \
          slab_ice_h.o \
          temps_mod.o
	fcm_internal compile:F dyn_phys_sub $< $@

start2archive_128x96x23_phystd_seq.e: \
          start2archive.o \
          LD__dyn_phys_sub__start2archive.flags \
          LDFLAGS__dyn_phys_sub__start2archive.flags \
          $(OBJECTS) \
          callkeys_mod.done \
          comconst_mod.done \
          comdissip.h.idone \
          comgeom.h.idone \
          comsoil_h.done \
          comvert_mod.done \
          control_mod.done \
          dimensions.h.idone \
          filtreg_mod.done \
          infotrac.done \
          iniphysiq_mod.done \
          ioipsl_getincom.done \
          mod_const_mpi.done \
          paramet.h.idone \
          phyetat0_mod.done \
          planete_mod.done \
          slab_ice_h.done \
          temps_mod.done
	fcm_internal load dyn_phys_sub $< $@

FFLAGS__dyn_phys_sub__defrun_new.flags: \
          FFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

defrun_new.done: \
          defrun_new.o \
          control_mod.done \
          ioipsl_getincom.done \
          logic_mod.done \
          serre_mod.done \
          sponge_mod.done
	touch $(FCM_DONEDIR)/$@

defrun_new.o: \
          $(PPSRCDIR0__dyn_phys_sub)/defrun_new.f \
          FFLAGS__dyn_phys_sub__defrun_new.flags \
          control_mod.o \
          ioipsl_getincom.o \
          logic_mod.o \
          serre_mod.o \
          sponge_mod.o
	fcm_internal compile:F dyn_phys_sub $< $@

FFLAGS__dyn_phys_sub__exner_hyb.flags: \
          FFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

exner_hyb.done: \
          exner_hyb.o \
          comconst_mod.done \
          comvert_mod.done
	touch $(FCM_DONEDIR)/$@

exner_hyb.o: \
          $(PPSRCDIR0__dyn_phys_sub)/exner_hyb.f \
          FFLAGS__dyn_phys_sub__exner_hyb.flags \
          comconst_mod.o \
          comvert_mod.o
	fcm_internal compile:F dyn_phys_sub $< $@

FFLAGS__dyn_phys_sub__grid_noro1.flags: \
          FFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

grid_noro1.done: \
          grid_noro1.o \
          comconst_mod.done
	touch $(FCM_DONEDIR)/$@

grid_noro1.o: \
          $(PPSRCDIR0__dyn_phys_sub)/grid_noro1.f \
          FFLAGS__dyn_phys_sub__grid_noro1.flags \
          comconst_mod.o
	fcm_internal compile:F dyn_phys_sub $< $@

FFLAGS__dyn_phys_sub__interp_vert.flags: \
          FFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

interp_vert.done: \
          interp_vert.o
	touch $(FCM_DONEDIR)/$@

interp_vert.o: \
          $(PPSRCDIR0__dyn_phys_sub)/interp_vert.f \
          FFLAGS__dyn_phys_sub__interp_vert.flags
	fcm_internal compile:F dyn_phys_sub $< $@

FFLAGS__dyn_phys_sub__callphysiq_mod.flags: \
          FFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

callphysiq_mod.done: \
          callphysiq_mod.o \
          control_mod.done \
          mod_grid_phy_lmdz.done \
          physiq_mod.done
	touch $(FCM_DONEDIR)/$@

callphysiq_mod.o: \
          $(PPSRCDIR0__dyn_phys_sub)/callphysiq_mod.f90 \
          FFLAGS__dyn_phys_sub__callphysiq_mod.flags \
          control_mod.o \
          mod_grid_phy_lmdz.o \
          physiq_mod.o
	fcm_internal compile:F dyn_phys_sub $< $@

FFLAGS__dyn_phys_sub__newstart.flags: \
          FFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

LDFLAGS__dyn_phys_sub__newstart.flags: \
          LDFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

LD__dyn_phys_sub__newstart.flags: \
          LD__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

newstart.o: \
          $(PPSRCDIR0__dyn_phys_sub)/newstart.f \
          FFLAGS__dyn_phys_sub__newstart.flags \
          comconst_mod.o \
          comdissnew.h \
          comgeom2.h \
          comsoil_h.o \
          comvert_mod.o \
          control_mod.o \
          datafile_mod.o \
          dimensions.h \
          ener_mod.o \
          filtreg_mod.o \
          infotrac.o \
          iniphysiq_mod.o \
          ioipsl_getin_p_mod.o \
          iostart.o \
          mod_const_mpi.o \
          mod_phys_lmdz_para.o \
          paramet.h \
          phyetat0_mod.o \
          phyredem.o \
          serre_mod.o \
          slab_ice_h.o \
          surfdat_h.o \
          tabfi_mod.o \
          temps_mod.o \
          tracer_h.o
	fcm_internal compile:F dyn_phys_sub $< $@

newstart_128x96x23_phystd_seq.e: \
          newstart.o \
          LD__dyn_phys_sub__newstart.flags \
          LDFLAGS__dyn_phys_sub__newstart.flags \
          $(OBJECTS) \
          comconst_mod.done \
          comdissnew.h.idone \
          comgeom2.h.idone \
          comsoil_h.done \
          comvert_mod.done \
          control_mod.done \
          datafile_mod.done \
          dimensions.h.idone \
          ener_mod.done \
          filtreg_mod.done \
          infotrac.done \
          iniphysiq_mod.done \
          ioipsl_getin_p_mod.done \
          iostart.done \
          mod_const_mpi.done \
          mod_phys_lmdz_para.done \
          paramet.h.idone \
          phyetat0_mod.done \
          phyredem.done \
          serre_mod.done \
          slab_ice_h.done \
          surfdat_h.done \
          tabfi_mod.done \
          temps_mod.done \
          tracer_h.done
	fcm_internal load dyn_phys_sub $< $@

FFLAGS__dyn_phys_sub__ini_archive.flags: \
          FFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

ini_archive.done: \
          ini_archive.o \
          comconst_mod.done \
          comsoil_h.done \
          comvert_mod.done \
          slab_ice_h.done
	touch $(FCM_DONEDIR)/$@

ini_archive.o: \
          $(PPSRCDIR0__dyn_phys_sub)/ini_archive.f \
          FFLAGS__dyn_phys_sub__ini_archive.flags \
          comconst_mod.o \
          comsoil_h.o \
          comvert_mod.o \
          slab_ice_h.o
	fcm_internal compile:F dyn_phys_sub $< $@

FFLAGS__dyn_phys_sub__wind_scal.flags: \
          FFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

wind_scal.done: \
          wind_scal.o
	touch $(FCM_DONEDIR)/$@

wind_scal.o: \
          $(PPSRCDIR0__dyn_phys_sub)/wind_scal.f \
          FFLAGS__dyn_phys_sub__wind_scal.flags
	fcm_internal compile:F dyn_phys_sub $< $@

FFLAGS__dyn_phys_sub__write_archive.flags: \
          FFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

write_archive.done: \
          write_archive.o \
          comsoil_h.done \
          slab_ice_h.done
	touch $(FCM_DONEDIR)/$@

write_archive.o: \
          $(PPSRCDIR0__dyn_phys_sub)/write_archive.f \
          FFLAGS__dyn_phys_sub__write_archive.flags \
          comsoil_h.o \
          slab_ice_h.o
	fcm_internal compile:F dyn_phys_sub $< $@

FFLAGS__dyn_phys_sub__lect_start_archive.flags: \
          FFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

lect_start_archive.done: \
          lect_start_archive.o \
          callkeys_mod.done \
          comconst_mod.done \
          comsoil_h.done \
          comvert_mod.done \
          infotrac.done \
          slab_ice_h.done \
          tracer_h.done
	touch $(FCM_DONEDIR)/$@

lect_start_archive.o: \
          $(PPSRCDIR0__dyn_phys_sub)/lect_start_archive.f \
          FFLAGS__dyn_phys_sub__lect_start_archive.flags \
          callkeys_mod.o \
          comconst_mod.o \
          comsoil_h.o \
          comvert_mod.o \
          infotrac.o \
          slab_ice_h.o \
          tracer_h.o
	fcm_internal compile:F dyn_phys_sub $< $@

FFLAGS__dyn_phys_sub__datareadnc.flags: \
          FFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

datareadnc.done: \
          datareadnc.o \
          comconst_mod.done \
          datafile_mod.done \
          ioipsl_getincom.done
	touch $(FCM_DONEDIR)/$@

datareadnc.o: \
          $(PPSRCDIR0__dyn_phys_sub)/datareadnc.f \
          FFLAGS__dyn_phys_sub__datareadnc.flags \
          comconst_mod.o \
          datafile_mod.o \
          ioipsl_getincom.o
	fcm_internal compile:F dyn_phys_sub $< $@

FFLAGS__dyn_phys_sub__caldyn0.flags: \
          FFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

caldyn0.done: \
          caldyn0.o \
          comvert_mod.done
	touch $(FCM_DONEDIR)/$@

caldyn0.o: \
          $(PPSRCDIR0__dyn_phys_sub)/caldyn0.f \
          FFLAGS__dyn_phys_sub__caldyn0.flags \
          comvert_mod.o
	fcm_internal compile:F dyn_phys_sub $< $@

FFLAGS__dyn_phys_sub__scal_wind.flags: \
          FFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

scal_wind.done: \
          scal_wind.o
	touch $(FCM_DONEDIR)/$@

scal_wind.o: \
          $(PPSRCDIR0__dyn_phys_sub)/scal_wind.f \
          FFLAGS__dyn_phys_sub__scal_wind.flags
	fcm_internal compile:F dyn_phys_sub $< $@

FFLAGS__dyn_phys_sub__iniphysiq_mod.flags: \
          FFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

iniphysiq_mod.done: \
          iniphysiq_mod.o \
          comgeom.h.idone \
          comvert_mod.done \
          control_mod.done \
          dimensions.h.idone \
          dimphy.done \
          geometry_mod.done \
          infotrac.done \
          inifis_mod.done \
          inigeomphy_mod.done \
          iniprint.h.idone \
          ioipsl_getin_p_mod.done \
          mod_phys_lmdz_para.done \
          paramet.h.idone \
          planete_mod.done \
          surf_heat_transp_mod.done
	touch $(FCM_DONEDIR)/$@

iniphysiq_mod.o: \
          $(PPSRCDIR0__dyn_phys_sub)/iniphysiq_mod.f90 \
          FFLAGS__dyn_phys_sub__iniphysiq_mod.flags \
          comgeom.h \
          comvert_mod.o \
          control_mod.o \
          dimensions.h \
          dimphy.o \
          geometry_mod.o \
          infotrac.o \
          inifis_mod.o \
          inigeomphy_mod.o \
          iniprint.h \
          ioipsl_getin_p_mod.o \
          mod_phys_lmdz_para.o \
          paramet.h \
          planete_mod.o \
          surf_heat_transp_mod.o
	fcm_internal compile:F dyn_phys_sub $< $@

FFLAGS__dyn_phys_sub__inichim_newstart.flags: \
          FFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

inichim_newstart.done: \
          inichim_newstart.o \
          callkeys_mod.done \
          comvert_mod.done \
          datafile_mod.done \
          mod_grid_phy_lmdz.done \
          tracer_h.done
	touch $(FCM_DONEDIR)/$@

inichim_newstart.o: \
          $(PPSRCDIR0__dyn_phys_sub)/inichim_newstart.f90 \
          FFLAGS__dyn_phys_sub__inichim_newstart.flags \
          callkeys_mod.o \
          comvert_mod.o \
          datafile_mod.o \
          mod_grid_phy_lmdz.o \
          tracer_h.o
	fcm_internal compile:F dyn_phys_sub $< $@

