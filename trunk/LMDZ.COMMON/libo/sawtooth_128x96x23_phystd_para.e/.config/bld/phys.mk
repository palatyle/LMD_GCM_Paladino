# Automatic Make rule for phys

SRCDIR0__phys = /home/palatyle/LMD_gen/trunk/LMDZ.COMMON/libf/phystd

PPSRCDIR0__phys = /home/palatyle/LMD_gen/trunk/LMDZ.COMMON/libo/sawtooth_128x96x23_phystd_para.e/.config/ppsrc/phys

FFLAGS__phys__initracer.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

initracer.done: \
          initracer.o \
          callkeys_mod.done \
          ioipsl_getin_p_mod.done \
          surfdat_h.done \
          tracer_h.done
	touch $(FCM_DONEDIR)/$@

initracer.o: \
          $(PPSRCDIR0__phys)/initracer.f \
          FFLAGS__phys__initracer.flags \
          callkeys_mod.o \
          ioipsl_getin_p_mod.o \
          surfdat_h.o \
          tracer_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__mkstat.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

mkstats.done: \
          mkstats.o \
          mod_grid_phy_lmdz.done \
          mod_phys_lmdz_para.done \
          statto_mod.done
	touch $(FCM_DONEDIR)/$@

mkstats.o: \
          $(PPSRCDIR0__phys)/mkstat.f90 \
          FFLAGS__phys__mkstat.flags \
          mod_grid_phy_lmdz.o \
          mod_phys_lmdz_para.o \
          statto_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__cpdet_phy_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

cpdet_phy_mod.done: \
          cpdet_phy_mod.o
	touch $(FCM_DONEDIR)/$@

cpdet_phy_mod.o: \
          $(PPSRCDIR0__phys)/cpdet_phy_mod.f90 \
          FFLAGS__phys__cpdet_phy_mod.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__callcorrk.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

callcorrk_mod.done: \
          callcorrk_mod.o \
          aerosol_mod.done \
          callkeys_mod.done \
          comcstfi_mod.done \
          datafile_mod.done \
          gases_h.done \
          ioipsl_getin_p_mod.done \
          mod_phys_lmdz_para.done \
          optci_mod.done \
          optcv_mod.done \
          radcommon_h.done \
          radii_mod.done \
          radinc_h.done \
          tracer_h.done \
          watercommon_h.done
	touch $(FCM_DONEDIR)/$@

callcorrk_mod.o: \
          $(PPSRCDIR0__phys)/callcorrk.f90 \
          FFLAGS__phys__callcorrk.flags \
          aerosol_mod.o \
          callkeys_mod.o \
          comcstfi_mod.o \
          datafile_mod.o \
          gases_h.o \
          ioipsl_getin_p_mod.o \
          mod_phys_lmdz_para.o \
          optci_mod.o \
          optcv_mod.o \
          radcommon_h.o \
          radii_mod.o \
          radinc_h.o \
          tracer_h.o \
          watercommon_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__interpolateN2H2.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

interpolaten2h2.done: \
          interpolaten2h2.o \
          datafile_mod.done
	touch $(FCM_DONEDIR)/$@

interpolaten2h2.o: \
          $(PPSRCDIR0__phys)/interpolateN2H2.f90 \
          FFLAGS__phys__interpolateN2H2.flags \
          datafile_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__moistadj.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

moistadj.done: \
          moistadj.o \
          comcstfi_mod.done \
          tracer_h.done \
          watercommon_h.done
	touch $(FCM_DONEDIR)/$@

moistadj.o: \
          $(PPSRCDIR0__phys)/moistadj.f90 \
          FFLAGS__phys__moistadj.flags \
          comcstfi_mod.o \
          tracer_h.o \
          watercommon_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__cvmgp.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

cvmgp.done: \
          cvmgp.o
	touch $(FCM_DONEDIR)/$@

cvmgp.o: \
          $(PPSRCDIR0__phys)/cvmgp.f \
          FFLAGS__phys__cvmgp.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__bilinear.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

bilinear.done: \
          bilinear.o
	touch $(FCM_DONEDIR)/$@

bilinear.o: \
          $(PPSRCDIR0__phys)/bilinear.f90 \
          FFLAGS__phys__bilinear.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__radii_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

radii_mod.done: \
          radii_mod.o \
          aerosol_mod.done \
          callkeys_mod.done \
          comcstfi_mod.done \
          ioipsl_getin_p_mod.done \
          radinc_h.done \
          tracer_h.done \
          watercommon_h.done
	touch $(FCM_DONEDIR)/$@

radii_mod.o: \
          $(PPSRCDIR0__phys)/radii_mod.f90 \
          FFLAGS__phys__radii_mod.flags \
          aerosol_mod.o \
          callkeys_mod.o \
          comcstfi_mod.o \
          ioipsl_getin_p_mod.o \
          radinc_h.o \
          tracer_h.o \
          watercommon_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__stellarlong.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

stellarlong.done: \
          stellarlong.o \
          comcstfi_mod.done \
          planete_mod.done
	touch $(FCM_DONEDIR)/$@

stellarlong.o: \
          $(PPSRCDIR0__phys)/stellarlong.f \
          FFLAGS__phys__stellarlong.flags \
          comcstfi_mod.o \
          planete_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__radcommon_h.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

radcommon_h.done: \
          radcommon_h.o \
          radinc_h.done
	touch $(FCM_DONEDIR)/$@

radcommon_h.o: \
          $(PPSRCDIR0__phys)/radcommon_h.f90 \
          FFLAGS__phys__radcommon_h.flags \
          radinc_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__orbite.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

orbite.done: \
          orbite.o \
          comcstfi_mod.done \
          planete_mod.done
	touch $(FCM_DONEDIR)/$@

orbite.o: \
          $(PPSRCDIR0__phys)/orbite.f \
          FFLAGS__phys__orbite.flags \
          comcstfi_mod.o \
          planete_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__vdifc.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

vdifc.done: \
          vdifc.o \
          callkeys_mod.done \
          comcstfi_mod.done \
          radcommon_h.done \
          surfdat_h.done \
          tracer_h.done \
          watercommon_h.done
	touch $(FCM_DONEDIR)/$@

vdifc.o: \
          $(PPSRCDIR0__phys)/vdifc.f \
          FFLAGS__phys__vdifc.flags \
          callkeys_mod.o \
          comcstfi_mod.o \
          radcommon_h.o \
          surfdat_h.o \
          tracer_h.o \
          watercommon_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__planete_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

planete_mod.done: \
          planete_mod.o
	touch $(FCM_DONEDIR)/$@

planete_mod.o: \
          $(PPSRCDIR0__phys)/planete_mod.f90 \
          FFLAGS__phys__planete_mod.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__thermcell_dq.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

thermcell_dq.done: \
          thermcell_dq.o \
          print_control_mod.done \
          thermcell_mod.done
	touch $(FCM_DONEDIR)/$@

thermcell_dq.o: \
          $(PPSRCDIR0__phys)/thermcell_dq.f90 \
          FFLAGS__phys__thermcell_dq.flags \
          print_control_mod.o \
          thermcell_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__tpindex.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

tpindex.done: \
          tpindex.o \
          radinc_h.done
	touch $(FCM_DONEDIR)/$@

tpindex.o: \
          $(PPSRCDIR0__phys)/tpindex.f \
          FFLAGS__phys__tpindex.flags \
          radinc_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__multipl.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

multipl.done: \
          multipl.o
	touch $(FCM_DONEDIR)/$@

multipl.o: \
          $(PPSRCDIR0__phys)/multipl.f \
          FFLAGS__phys__multipl.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__phys_state_var_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

phys_state_var_mod.done: \
          phys_state_var_mod.o \
          comsaison_h.done \
          comsoil_h.done \
          dimphy.done \
          radcommon_h.done \
          radinc_h.done \
          slab_ice_h.done \
          surfdat_h.done \
          turb_mod.done
	touch $(FCM_DONEDIR)/$@

phys_state_var_mod.o: \
          $(PPSRCDIR0__phys)/phys_state_var_mod.f90 \
          FFLAGS__phys__phys_state_var_mod.flags \
          comsaison_h.o \
          comsoil_h.o \
          dimphy.o \
          radcommon_h.o \
          radinc_h.o \
          slab_ice_h.o \
          surfdat_h.o \
          turb_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__evap.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

evap.done: \
          evap.o \
          tracer_h.done \
          watercommon_h.done
	touch $(FCM_DONEDIR)/$@

evap.o: \
          $(PPSRCDIR0__phys)/evap.f \
          FFLAGS__phys__evap.flags \
          tracer_h.o \
          watercommon_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__profile.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

profile.done: \
          profile.o \
          ioipsl_getin_p_mod.done
	touch $(FCM_DONEDIR)/$@

profile.o: \
          $(PPSRCDIR0__phys)/profile.f \
          FFLAGS__phys__profile.flags \
          ioipsl_getin_p_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__writediagfi.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

writediagfi.done: \
          writediagfi.o \
          callkeys_mod.done \
          geometry_mod.done \
          mod_grid_phy_lmdz.done \
          mod_phys_lmdz_para.done \
          surfdat_h.done \
          time_phylmdz_mod.done
	touch $(FCM_DONEDIR)/$@

writediagfi.o: \
          $(PPSRCDIR0__phys)/writediagfi.f \
          FFLAGS__phys__writediagfi.flags \
          callkeys_mod.o \
          geometry_mod.o \
          mod_grid_phy_lmdz.o \
          mod_phys_lmdz_para.o \
          surfdat_h.o \
          time_phylmdz_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__conc_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

conc_mod.done: \
          conc_mod.o
	touch $(FCM_DONEDIR)/$@

conc_mod.o: \
          $(PPSRCDIR0__phys)/conc_mod.f90 \
          FFLAGS__phys__conc_mod.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__planetwide_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

planetwide_mod.done: \
          planetwide_mod.o \
          dimphy.done \
          mod_grid_phy_lmdz.done \
          mod_phys_lmdz_para.done
	touch $(FCM_DONEDIR)/$@

planetwide_mod.o: \
          $(PPSRCDIR0__phys)/planetwide_mod.f90 \
          FFLAGS__phys__planetwide_mod.flags \
          dimphy.o \
          mod_grid_phy_lmdz.o \
          mod_phys_lmdz_para.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__bilinearbig.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

bilinearbig.done: \
          bilinearbig.o
	touch $(FCM_DONEDIR)/$@

bilinearbig.o: \
          $(PPSRCDIR0__phys)/bilinearbig.f90 \
          FFLAGS__phys__bilinearbig.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__radinc_h.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

radinc_h.done: \
          radinc_h.o \
          bands.h.idone \
          scatterers.h.idone
	touch $(FCM_DONEDIR)/$@

radinc_h.o: \
          $(PPSRCDIR0__phys)/radinc_h.f90 \
          FFLAGS__phys__radinc_h.flags \
          bands.h \
          scatterers.h
	fcm_internal compile:F phys $< $@

FFLAGS__phys__aerosol_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

aerosol_mod.done: \
          aerosol_mod.o
	touch $(FCM_DONEDIR)/$@

aerosol_mod.o: \
          $(PPSRCDIR0__phys)/aerosol_mod.f90 \
          FFLAGS__phys__aerosol_mod.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__wstats.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

wstats.done: \
          wstats.o \
          mod_grid_phy_lmdz.done \
          mod_phys_lmdz_para.done \
          statto_mod.done
	touch $(FCM_DONEDIR)/$@

wstats.o: \
          $(PPSRCDIR0__phys)/wstats.f90 \
          FFLAGS__phys__wstats.flags \
          mod_grid_phy_lmdz.o \
          mod_phys_lmdz_para.o \
          statto_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__turbdiff.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

turbdiff.done: \
          turbdiff.o \
          callkeys_mod.done \
          comcstfi_mod.done \
          radcommon_h.done \
          surfdat_h.done \
          tracer_h.done \
          turb_mod.done \
          watercommon_h.done
	touch $(FCM_DONEDIR)/$@

turbdiff.o: \
          $(PPSRCDIR0__phys)/turbdiff.f90 \
          FFLAGS__phys__turbdiff.flags \
          callkeys_mod.o \
          comcstfi_mod.o \
          radcommon_h.o \
          surfdat_h.o \
          tracer_h.o \
          turb_mod.o \
          watercommon_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__surf_heat_transp_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

surf_heat_transp_mod.done: \
          surf_heat_transp_mod.o \
          mod_grid_phy_lmdz.done \
          slab_ice_h.done
	touch $(FCM_DONEDIR)/$@

surf_heat_transp_mod.o: \
          $(PPSRCDIR0__phys)/surf_heat_transp_mod.f90 \
          FFLAGS__phys__surf_heat_transp_mod.flags \
          mod_grid_phy_lmdz.o \
          slab_ice_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__calc_cpp3d.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

calc_cpp3d.done: \
          calc_cpp3d.o \
          comcstfi_mod.done
	touch $(FCM_DONEDIR)/$@

calc_cpp3d.o: \
          $(PPSRCDIR0__phys)/calc_cpp3d.f90 \
          FFLAGS__phys__calc_cpp3d.flags \
          comcstfi_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__iostart.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

iostart.done: \
          iostart.o \
          comsoil_h.done \
          dimphy.done \
          mod_grid_phy_lmdz.done \
          mod_phys_lmdz_para.done \
          slab_ice_h.done \
          tracer_h.done
	touch $(FCM_DONEDIR)/$@

iostart.o: \
          $(PPSRCDIR0__phys)/iostart.f90 \
          FFLAGS__phys__iostart.flags \
          comsoil_h.o \
          dimphy.o \
          mod_grid_phy_lmdz.o \
          mod_phys_lmdz_para.o \
          slab_ice_h.o \
          tracer_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__comsoil_h.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

comsoil_h.done: \
          comsoil_h.o
	touch $(FCM_DONEDIR)/$@

comsoil_h.o: \
          $(PPSRCDIR0__phys)/comsoil_h.f90 \
          FFLAGS__phys__comsoil_h.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__iniwrite_specVI.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

iniwrite_specvi.done: \
          iniwrite_specvi.o \
          comcstfi_mod.done \
          mod_grid_phy_lmdz.done \
          radcommon_h.done \
          radinc_h.done \
          regular_lonlat_mod.done \
          time_phylmdz_mod.done
	touch $(FCM_DONEDIR)/$@

iniwrite_specvi.o: \
          $(PPSRCDIR0__phys)/iniwrite_specVI.f \
          FFLAGS__phys__iniwrite_specVI.flags \
          comcstfi_mod.o \
          mod_grid_phy_lmdz.o \
          radcommon_h.o \
          radinc_h.o \
          regular_lonlat_mod.o \
          time_phylmdz_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__aerave.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

aerave.done: \
          aerave.o
	touch $(FCM_DONEDIR)/$@

aerave.o: \
          $(PPSRCDIR0__phys)/aerave.f \
          FFLAGS__phys__aerave.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__condense_co2.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

condense_co2.done: \
          condense_co2.o \
          aerosol_mod.done \
          callkeys_mod.done \
          comcstfi_mod.done \
          gases_h.done \
          geometry_mod.done \
          radii_mod.done \
          radinc_h.done \
          surfdat_h.done \
          tracer_h.done
	touch $(FCM_DONEDIR)/$@

condense_co2.o: \
          $(PPSRCDIR0__phys)/condense_co2.f90 \
          FFLAGS__phys__condense_co2.flags \
          aerosol_mod.o \
          callkeys_mod.o \
          comcstfi_mod.o \
          gases_h.o \
          geometry_mod.o \
          radii_mod.o \
          radinc_h.o \
          surfdat_h.o \
          tracer_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__hydrol.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

hydrol.done: \
          hydrol.o \
          callkeys_mod.done \
          comdiurn_h.done \
          geometry_mod.done \
          ioipsl_getin_p_mod.done \
          radinc_h.done \
          slab_ice_h.done \
          surfdat_h.done \
          tracer_h.done \
          watercommon_h.done
	touch $(FCM_DONEDIR)/$@

hydrol.o: \
          $(PPSRCDIR0__phys)/hydrol.f90 \
          FFLAGS__phys__hydrol.flags \
          callkeys_mod.o \
          comdiurn_h.o \
          geometry_mod.o \
          ioipsl_getin_p_mod.o \
          radinc_h.o \
          slab_ice_h.o \
          surfdat_h.o \
          tracer_h.o \
          watercommon_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__ocean_slab_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

ocean_slab_mod.done: \
          ocean_slab_mod.o \
          callkeys_mod.done \
          slab_ice_h.done \
          surf_heat_transp_mod.done \
          watercommon_h.done
	touch $(FCM_DONEDIR)/$@

ocean_slab_mod.o: \
          $(PPSRCDIR0__phys)/ocean_slab_mod.f90 \
          FFLAGS__phys__ocean_slab_mod.flags \
          callkeys_mod.o \
          slab_ice_h.o \
          surf_heat_transp_mod.o \
          watercommon_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__thermcell_dv2.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

thermcell_dv2.done: \
          thermcell_dv2.o
	touch $(FCM_DONEDIR)/$@

thermcell_dv2.o: \
          $(PPSRCDIR0__phys)/thermcell_dv2.f90 \
          FFLAGS__phys__thermcell_dv2.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__iniwrite.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

iniwrite.done: \
          iniwrite.o \
          comcstfi_mod.done \
          comsoil_h.done \
          mod_grid_phy_lmdz.done \
          regular_lonlat_mod.done \
          time_phylmdz_mod.done \
          vertical_layers_mod.done
	touch $(FCM_DONEDIR)/$@

iniwrite.o: \
          $(PPSRCDIR0__phys)/iniwrite.f \
          FFLAGS__phys__iniwrite.flags \
          comcstfi_mod.o \
          comsoil_h.o \
          mod_grid_phy_lmdz.o \
          regular_lonlat_mod.o \
          time_phylmdz_mod.o \
          vertical_layers_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__tabfi_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

tabfi_mod.done: \
          tabfi_mod.o \
          callkeys_mod.done \
          comcstfi_mod.done \
          comsoil_h.done \
          ioipsl_getin_p_mod.done \
          iostart.done \
          mod_phys_lmdz_para.done \
          planete_mod.done \
          surfdat_h.done \
          time_phylmdz_mod.done
	touch $(FCM_DONEDIR)/$@

tabfi_mod.o: \
          $(PPSRCDIR0__phys)/tabfi_mod.f90 \
          FFLAGS__phys__tabfi_mod.flags \
          callkeys_mod.o \
          comcstfi_mod.o \
          comsoil_h.o \
          ioipsl_getin_p_mod.o \
          iostart.o \
          mod_phys_lmdz_para.o \
          planete_mod.o \
          surfdat_h.o \
          time_phylmdz_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__calc_cpp_mugaz.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

calc_cpp_mugaz.done: \
          calc_cpp_mugaz.o \
          callkeys_mod.done \
          comcstfi_mod.done \
          gases_h.done
	touch $(FCM_DONEDIR)/$@

calc_cpp_mugaz.o: \
          $(PPSRCDIR0__phys)/calc_cpp_mugaz.f90 \
          FFLAGS__phys__calc_cpp_mugaz.flags \
          callkeys_mod.o \
          comcstfi_mod.o \
          gases_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__comg1d_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

comg1d_mod.done: \
          comg1d_mod.o
	touch $(FCM_DONEDIR)/$@

comg1d_mod.o: \
          $(PPSRCDIR0__phys)/comg1d_mod.f90 \
          FFLAGS__phys__comg1d_mod.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__comgeomfi_h.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

comgeomfi_h.done: \
          comgeomfi_h.o
	touch $(FCM_DONEDIR)/$@

comgeomfi_h.o: \
          $(PPSRCDIR0__phys)/comgeomfi_h.f90 \
          FFLAGS__phys__comgeomfi_h.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__sfluxv.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

sfluxv.done: \
          sfluxv.o \
          radcommon_h.done \
          radinc_h.done
	touch $(FCM_DONEDIR)/$@

sfluxv.o: \
          $(PPSRCDIR0__phys)/sfluxv.f \
          FFLAGS__phys__sfluxv.flags \
          radcommon_h.o \
          radinc_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__writediagspecVI.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

writediagspecvi.done: \
          writediagspecvi.o \
          callkeys_mod.done \
          geometry_mod.done \
          mod_grid_phy_lmdz.done \
          mod_phys_lmdz_para.done \
          radinc_h.done \
          time_phylmdz_mod.done
	touch $(FCM_DONEDIR)/$@

writediagspecvi.o: \
          $(PPSRCDIR0__phys)/writediagspecVI.f \
          FFLAGS__phys__writediagspecVI.flags \
          callkeys_mod.o \
          geometry_mod.o \
          mod_grid_phy_lmdz.o \
          mod_phys_lmdz_para.o \
          radinc_h.o \
          time_phylmdz_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__comcstfi_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

comcstfi_mod.done: \
          comcstfi_mod.o
	touch $(FCM_DONEDIR)/$@

comcstfi_mod.o: \
          $(PPSRCDIR0__phys)/comcstfi_mod.f90 \
          FFLAGS__phys__comcstfi_mod.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__thermcell_env.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

thermcell_env.done: \
          thermcell_env.o \
          callkeys_mod.done \
          print_control_mod.done \
          thermcell_mod.done \
          tracer_h.done \
          watercommon_h.done
	touch $(FCM_DONEDIR)/$@

thermcell_env.o: \
          $(PPSRCDIR0__phys)/thermcell_env.f90 \
          FFLAGS__phys__thermcell_env.flags \
          callkeys_mod.o \
          print_control_mod.o \
          thermcell_mod.o \
          tracer_h.o \
          watercommon_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__surfdat_h.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

surfdat_h.done: \
          surfdat_h.o
	touch $(FCM_DONEDIR)/$@

surfdat_h.o: \
          $(PPSRCDIR0__phys)/surfdat_h.f90 \
          FFLAGS__phys__surfdat_h.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__writeg1d.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

writeg1d.done: \
          writeg1d.o \
          comg1d_mod.done \
          time_phylmdz_mod.done
	touch $(FCM_DONEDIR)/$@

writeg1d.o: \
          $(PPSRCDIR0__phys)/writeg1d.f \
          FFLAGS__phys__writeg1d.flags \
          comg1d_mod.o \
          time_phylmdz_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__interpolateH2Ocont_CKD.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

interpolateh2ocont_ckd.done: \
          interpolateh2ocont_ckd.o \
          datafile_mod.done
	touch $(FCM_DONEDIR)/$@

interpolateh2ocont_ckd.o: \
          $(PPSRCDIR0__phys)/interpolateH2Ocont_CKD.f90 \
          FFLAGS__phys__interpolateH2Ocont_CKD.flags \
          datafile_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__newtrelax.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

newtrelax.done: \
          newtrelax.o \
          callkeys_mod.done \
          comcstfi_mod.done
	touch $(FCM_DONEDIR)/$@

newtrelax.o: \
          $(PPSRCDIR0__phys)/newtrelax.f90 \
          FFLAGS__phys__newtrelax.flags \
          callkeys_mod.o \
          comcstfi_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__turb_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

turb_mod.done: \
          turb_mod.o
	touch $(FCM_DONEDIR)/$@

turb_mod.o: \
          $(PPSRCDIR0__phys)/turb_mod.f90 \
          FFLAGS__phys__turb_mod.flags
	fcm_internal compile:F phys $< $@

bands.h: \
          $(SRCDIR0__phys)/bands.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

bands.h.idone: \
          $(SRCDIR0__phys)/bands.h
	touch $(FCM_DONEDIR)/$@

FFLAGS__phys__def_var.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

def_var.done: \
          def_var.o
	touch $(FCM_DONEDIR)/$@

def_var.o: \
          $(PPSRCDIR0__phys)/def_var.f90 \
          FFLAGS__phys__def_var.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__interpolateH2Ocont_PPC.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

interpolateh2ocont_ppc.done: \
          interpolateh2ocont_ppc.o \
          datafile_mod.done \
          watercommon_h.done
	touch $(FCM_DONEDIR)/$@

interpolateh2ocont_ppc.o: \
          $(PPSRCDIR0__phys)/interpolateH2Ocont_PPC.f90 \
          FFLAGS__phys__interpolateH2Ocont_PPC.flags \
          datafile_mod.o \
          watercommon_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__ave_stelspec.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

ave_stelspec.done: \
          ave_stelspec.o \
          callkeys_mod.done \
          datafile_mod.done \
          radcommon_h.done \
          radinc_h.done
	touch $(FCM_DONEDIR)/$@

ave_stelspec.o: \
          $(PPSRCDIR0__phys)/ave_stelspec.f90 \
          FFLAGS__phys__ave_stelspec.flags \
          callkeys_mod.o \
          datafile_mod.o \
          radcommon_h.o \
          radinc_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__setspv.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

setspv.done: \
          setspv.o \
          callkeys_mod.done \
          datafile_mod.done \
          radcommon_h.done \
          radinc_h.done
	touch $(FCM_DONEDIR)/$@

setspv.o: \
          $(PPSRCDIR0__phys)/setspv.f90 \
          FFLAGS__phys__setspv.flags \
          callkeys_mod.o \
          datafile_mod.o \
          radcommon_h.o \
          radinc_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__stelang.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

stelang.done: \
          stelang.o
	touch $(FCM_DONEDIR)/$@

stelang.o: \
          $(PPSRCDIR0__phys)/stelang.f \
          FFLAGS__phys__stelang.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__physiq_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

physiq_mod.done: \
          physiq_mod.o \
          aerosol_mod.done \
          callcorrk_mod.done \
          callkeys_mod.done \
          comcstfi_mod.done \
          comdiurn_h.done \
          comgeomfi_h.done \
          comsaison_h.done \
          comsoil_h.done \
          conc_mod.done \
          datafile_mod.done \
          dimensions.h.idone \
          gases_h.done \
          geometry_mod.done \
          mod_phys_lmdz_mpi_data.done \
          mod_phys_lmdz_omp_data.done \
          mod_phys_lmdz_para.done \
          ocean_slab_mod.done \
          phyetat0_mod.done \
          phyredem.done \
          phys_state_var_mod.done \
          planete_mod.done \
          planetwide_mod.done \
          print_control_mod.done \
          radcommon_h.done \
          radii_mod.done \
          radinc_h.done \
          slab_ice_h.done \
          surf_heat_transp_mod.done \
          surfdat_h.done \
          thermcell_mod.done \
          time_phylmdz_mod.done \
          tracer_h.done \
          turb_mod.done \
          vertical_layers_mod.done \
          watercommon_h.done
	touch $(FCM_DONEDIR)/$@

physiq_mod.o: \
          $(PPSRCDIR0__phys)/physiq_mod.f90 \
          FFLAGS__phys__physiq_mod.flags \
          aerosol_mod.o \
          callcorrk_mod.o \
          callkeys_mod.o \
          comcstfi_mod.o \
          comdiurn_h.o \
          comgeomfi_h.o \
          comsaison_h.o \
          comsoil_h.o \
          conc_mod.o \
          datafile_mod.o \
          dimensions.h \
          gases_h.o \
          geometry_mod.o \
          mod_phys_lmdz_mpi_data.o \
          mod_phys_lmdz_omp_data.o \
          mod_phys_lmdz_para.o \
          ocean_slab_mod.o \
          phyetat0_mod.o \
          phyredem.o \
          phys_state_var_mod.o \
          planete_mod.o \
          planetwide_mod.o \
          print_control_mod.o \
          radcommon_h.o \
          radii_mod.o \
          radinc_h.o \
          slab_ice_h.o \
          surf_heat_transp_mod.o \
          surfdat_h.o \
          thermcell_mod.o \
          time_phylmdz_mod.o \
          tracer_h.o \
          turb_mod.o \
          vertical_layers_mod.o \
          watercommon_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__scopyi.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

scopyi.done: \
          scopyi.o
	touch $(FCM_DONEDIR)/$@

scopyi.o: \
          $(PPSRCDIR0__phys)/scopyi.f \
          FFLAGS__phys__scopyi.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__statto_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

statto_mod.done: \
          statto_mod.o
	touch $(FCM_DONEDIR)/$@

statto_mod.o: \
          $(PPSRCDIR0__phys)/statto_mod.f90 \
          FFLAGS__phys__statto_mod.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__aerave_new.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

aerave_new.done: \
          aerave_new.o
	touch $(FCM_DONEDIR)/$@

aerave_new.o: \
          $(PPSRCDIR0__phys)/aerave_new.f \
          FFLAGS__phys__aerave_new.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__stokes.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

stokes.done: \
          stokes.o \
          comcstfi_mod.done
	touch $(FCM_DONEDIR)/$@

stokes.o: \
          $(PPSRCDIR0__phys)/stokes.f90 \
          FFLAGS__phys__stokes.flags \
          comcstfi_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__inistats.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

inistats.done: \
          inistats.o \
          mod_grid_phy_lmdz.done \
          mod_phys_lmdz_para.done \
          nrtype.done \
          regular_lonlat_mod.done \
          statto_mod.done \
          time_phylmdz_mod.done \
          vertical_layers_mod.done
	touch $(FCM_DONEDIR)/$@

inistats.o: \
          $(PPSRCDIR0__phys)/inistats.f \
          FFLAGS__phys__inistats.flags \
          mod_grid_phy_lmdz.o \
          mod_phys_lmdz_para.o \
          nrtype.o \
          regular_lonlat_mod.o \
          statto_mod.o \
          time_phylmdz_mod.o \
          vertical_layers_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__callsedim.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

callsedim.done: \
          callsedim.o \
          aerosol_mod.done \
          callkeys_mod.done \
          comcstfi_mod.done \
          radii_mod.done \
          radinc_h.done \
          tracer_h.done
	touch $(FCM_DONEDIR)/$@

callsedim.o: \
          $(PPSRCDIR0__phys)/callsedim.f \
          FFLAGS__phys__callsedim.flags \
          aerosol_mod.o \
          callkeys_mod.o \
          comcstfi_mod.o \
          radii_mod.o \
          radinc_h.o \
          tracer_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__thermcell_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

thermcell_mod.done: \
          thermcell_mod.o \
          ioipsl_getin_p_mod.done
	touch $(FCM_DONEDIR)/$@

thermcell_mod.o: \
          $(PPSRCDIR0__phys)/thermcell_mod.f90 \
          FFLAGS__phys__thermcell_mod.flags \
          ioipsl_getin_p_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__cp_neutral.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

cp_neutral.done: \
          cp_neutral.o \
          gases_h.done
	touch $(FCM_DONEDIR)/$@

cp_neutral.o: \
          $(PPSRCDIR0__phys)/cp_neutral.f90 \
          FFLAGS__phys__cp_neutral.flags \
          gases_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__newsedim.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

newsedim.done: \
          newsedim.o \
          comcstfi_mod.done
	touch $(FCM_DONEDIR)/$@

newsedim.o: \
          $(PPSRCDIR0__phys)/newsedim.f \
          FFLAGS__phys__newsedim.flags \
          comcstfi_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__phyredem.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

phyredem.done: \
          phyredem.o \
          callkeys_mod.done \
          comcstfi_mod.done \
          comsoil_h.done \
          geometry_mod.done \
          iostart.done \
          mod_grid_phy_lmdz.done \
          planete_mod.done \
          slab_ice_h.done \
          surfdat_h.done \
          time_phylmdz_mod.done \
          tracer_h.done
	touch $(FCM_DONEDIR)/$@

phyredem.o: \
          $(PPSRCDIR0__phys)/phyredem.f90 \
          FFLAGS__phys__phyredem.flags \
          callkeys_mod.o \
          comcstfi_mod.o \
          comsoil_h.o \
          geometry_mod.o \
          iostart.o \
          mod_grid_phy_lmdz.o \
          planete_mod.o \
          slab_ice_h.o \
          surfdat_h.o \
          time_phylmdz_mod.o \
          tracer_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__tracer_h.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

tracer_h.done: \
          tracer_h.o
	touch $(FCM_DONEDIR)/$@

tracer_h.o: \
          $(PPSRCDIR0__phys)/tracer_h.f90 \
          FFLAGS__phys__tracer_h.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__thermcell_height.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

thermcell_height.done: \
          thermcell_height.o
	touch $(FCM_DONEDIR)/$@

thermcell_height.o: \
          $(PPSRCDIR0__phys)/thermcell_height.f90 \
          FFLAGS__phys__thermcell_height.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__volcano.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

volcano.done: \
          volcano.o \
          comgeomfi_h.done \
          datafile_mod.done \
          dimensions.h.idone \
          mod_phys_lmdz_mpi_data.done \
          print_control_mod.done \
          tracer_h.done
	touch $(FCM_DONEDIR)/$@

volcano.o: \
          $(PPSRCDIR0__phys)/volcano.f90 \
          FFLAGS__phys__volcano.flags \
          comgeomfi_h.o \
          datafile_mod.o \
          dimensions.h \
          mod_phys_lmdz_mpi_data.o \
          print_control_mod.o \
          tracer_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__vdif_kc.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

vdif_kc.done: \
          vdif_kc.o
	touch $(FCM_DONEDIR)/$@

vdif_kc.o: \
          $(PPSRCDIR0__phys)/vdif_kc.f \
          FFLAGS__phys__vdif_kc.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__iniwritesoil.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

iniwritesoil.done: \
          iniwritesoil.o \
          comcstfi_mod.done \
          comsoil_h.done \
          mod_grid_phy_lmdz.done \
          regular_lonlat_mod.done
	touch $(FCM_DONEDIR)/$@

iniwritesoil.o: \
          $(PPSRCDIR0__phys)/iniwritesoil.f90 \
          FFLAGS__phys__iniwritesoil.flags \
          comcstfi_mod.o \
          comsoil_h.o \
          mod_grid_phy_lmdz.o \
          regular_lonlat_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__iniwrite_specIR.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

iniwrite_specir.done: \
          iniwrite_specir.o \
          comcstfi_mod.done \
          mod_grid_phy_lmdz.done \
          radcommon_h.done \
          radinc_h.done \
          regular_lonlat_mod.done \
          time_phylmdz_mod.done
	touch $(FCM_DONEDIR)/$@

iniwrite_specir.o: \
          $(PPSRCDIR0__phys)/iniwrite_specIR.f \
          FFLAGS__phys__iniwrite_specIR.flags \
          comcstfi_mod.o \
          mod_grid_phy_lmdz.o \
          radcommon_h.o \
          radinc_h.o \
          regular_lonlat_mod.o \
          time_phylmdz_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__mass_redistribution.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

mass_redistribution.done: \
          mass_redistribution.o \
          callkeys_mod.done \
          comcstfi_mod.done \
          planete_mod.done \
          radcommon_h.done \
          surfdat_h.done \
          tracer_h.done \
          watercommon_h.done
	touch $(FCM_DONEDIR)/$@

mass_redistribution.o: \
          $(PPSRCDIR0__phys)/mass_redistribution.f90 \
          FFLAGS__phys__mass_redistribution.flags \
          callkeys_mod.o \
          comcstfi_mod.o \
          planete_mod.o \
          radcommon_h.o \
          surfdat_h.o \
          tracer_h.o \
          watercommon_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__blackl.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

blackl.done: \
          blackl.o
	touch $(FCM_DONEDIR)/$@

blackl.o: \
          $(PPSRCDIR0__phys)/blackl.f \
          FFLAGS__phys__blackl.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__sugas_corrk.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

sugas_corrk.done: \
          sugas_corrk.o \
          callkeys_mod.done \
          comcstfi_mod.done \
          datafile_mod.done \
          gases_h.done \
          ioipsl_getin_p_mod.done \
          radcommon_h.done \
          radinc_h.done
	touch $(FCM_DONEDIR)/$@

sugas_corrk.o: \
          $(PPSRCDIR0__phys)/sugas_corrk.f90 \
          FFLAGS__phys__sugas_corrk.flags \
          callkeys_mod.o \
          comcstfi_mod.o \
          datafile_mod.o \
          gases_h.o \
          ioipsl_getin_p_mod.o \
          radcommon_h.o \
          radinc_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__interpolateN2N2.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

interpolaten2n2.done: \
          interpolaten2n2.o \
          datafile_mod.done
	touch $(FCM_DONEDIR)/$@

interpolaten2n2.o: \
          $(PPSRCDIR0__phys)/interpolateN2N2.f90 \
          FFLAGS__phys__interpolateN2N2.flags \
          datafile_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__vdif_cd.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

vdif_cd.done: \
          vdif_cd.o
	touch $(FCM_DONEDIR)/$@

vdif_cd.o: \
          $(PPSRCDIR0__phys)/vdif_cd.f \
          FFLAGS__phys__vdif_cd.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__lagrange.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

lagrange.done: \
          lagrange.o
	touch $(FCM_DONEDIR)/$@

lagrange.o: \
          $(PPSRCDIR0__phys)/lagrange.f \
          FFLAGS__phys__lagrange.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__interpolateH2He.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

interpolateh2he.done: \
          interpolateh2he.o \
          datafile_mod.done
	touch $(FCM_DONEDIR)/$@

interpolateh2he.o: \
          $(PPSRCDIR0__phys)/interpolateH2He.f90 \
          FFLAGS__phys__interpolateH2He.flags \
          datafile_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__datafile_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

datafile_mod.done: \
          datafile_mod.o
	touch $(FCM_DONEDIR)/$@

datafile_mod.o: \
          $(PPSRCDIR0__phys)/datafile_mod.f90 \
          FFLAGS__phys__datafile_mod.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__callkeys_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

callkeys_mod.done: \
          callkeys_mod.o
	touch $(FCM_DONEDIR)/$@

callkeys_mod.o: \
          $(PPSRCDIR0__phys)/callkeys_mod.f90 \
          FFLAGS__phys__callkeys_mod.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__setspi.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

setspi.done: \
          setspi.o \
          comcstfi_mod.done \
          datafile_mod.done \
          radcommon_h.done \
          radinc_h.done
	touch $(FCM_DONEDIR)/$@

setspi.o: \
          $(PPSRCDIR0__phys)/setspi.f90 \
          FFLAGS__phys__setspi.flags \
          comcstfi_mod.o \
          datafile_mod.o \
          radcommon_h.o \
          radinc_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__vlz_fi.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

vlz_fi.done: \
          vlz_fi.o
	touch $(FCM_DONEDIR)/$@

vlz_fi.o: \
          $(PPSRCDIR0__phys)/vlz_fi.f \
          FFLAGS__phys__vlz_fi.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__writediagspecIR.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

writediagspecir.done: \
          writediagspecir.o \
          callkeys_mod.done \
          geometry_mod.done \
          mod_grid_phy_lmdz.done \
          mod_phys_lmdz_para.done \
          radinc_h.done \
          time_phylmdz_mod.done
	touch $(FCM_DONEDIR)/$@

writediagspecir.o: \
          $(PPSRCDIR0__phys)/writediagspecIR.f \
          FFLAGS__phys__writediagspecIR.flags \
          callkeys_mod.o \
          geometry_mod.o \
          mod_grid_phy_lmdz.o \
          mod_phys_lmdz_para.o \
          radinc_h.o \
          time_phylmdz_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__time_phylmdz_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

time_phylmdz_mod.done: \
          time_phylmdz_mod.o
	touch $(FCM_DONEDIR)/$@

time_phylmdz_mod.o: \
          $(PPSRCDIR0__phys)/time_phylmdz_mod.f90 \
          FFLAGS__phys__time_phylmdz_mod.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__writediagsoil.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

writediagsoil.done: \
          writediagsoil.o \
          comsoil_h.done \
          geometry_mod.done \
          mod_grid_phy_lmdz.done \
          mod_phys_lmdz_para.done \
          time_phylmdz_mod.done
	touch $(FCM_DONEDIR)/$@

writediagsoil.o: \
          $(PPSRCDIR0__phys)/writediagsoil.f90 \
          FFLAGS__phys__writediagsoil.flags \
          comsoil_h.o \
          geometry_mod.o \
          mod_grid_phy_lmdz.o \
          mod_phys_lmdz_para.o \
          time_phylmdz_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__soil_settings.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

soil_settings.done: \
          soil_settings.o \
          comsoil_h.done \
          iostart.done
	touch $(FCM_DONEDIR)/$@

soil_settings.o: \
          $(PPSRCDIR0__phys)/soil_settings.f \
          FFLAGS__phys__soil_settings.flags \
          comsoil_h.o \
          iostart.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__comsaison_h.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

comsaison_h.done: \
          comsaison_h.o
	touch $(FCM_DONEDIR)/$@

comsaison_h.o: \
          $(PPSRCDIR0__phys)/comsaison_h.f90 \
          FFLAGS__phys__comsaison_h.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__interpolateH2H2.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

interpolateh2h2.done: \
          interpolateh2h2.o \
          callkeys_mod.done \
          datafile_mod.done
	touch $(FCM_DONEDIR)/$@

interpolateh2h2.o: \
          $(PPSRCDIR0__phys)/interpolateH2H2.f90 \
          FFLAGS__phys__interpolateH2H2.flags \
          callkeys_mod.o \
          datafile_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__cvmgt.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

cvmgt.done: \
          cvmgt.o
	touch $(FCM_DONEDIR)/$@

cvmgt.o: \
          $(PPSRCDIR0__phys)/cvmgt.f \
          FFLAGS__phys__cvmgt.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__albedo_calcv.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

spectral_albedo_calc.done: \
          spectral_albedo_calc.o \
          callkeys_mod.done \
          radcommon_h.done \
          radinc_h.done
	touch $(FCM_DONEDIR)/$@

spectral_albedo_calc.o: \
          $(PPSRCDIR0__phys)/albedo_calcv.f90 \
          FFLAGS__phys__albedo_calcv.flags \
          callkeys_mod.o \
          radcommon_h.o \
          radinc_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__forceWCfn.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

forcewcfn.done: \
          forcewcfn.o \
          comcstfi_mod.done \
          tracer_h.done
	touch $(FCM_DONEDIR)/$@

forcewcfn.o: \
          $(PPSRCDIR0__phys)/forceWCfn.f \
          FFLAGS__phys__forceWCfn.flags \
          comcstfi_mod.o \
          tracer_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__phyetat0_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

phyetat0_mod.done: \
          phyetat0_mod.o \
          callkeys_mod.done \
          ioipsl_getin_p_mod.done \
          iostart.done \
          slab_ice_h.done \
          surfdat_h.done \
          tabfi_mod.done \
          tracer_h.done
	touch $(FCM_DONEDIR)/$@

phyetat0_mod.o: \
          $(PPSRCDIR0__phys)/phyetat0_mod.f90 \
          FFLAGS__phys__phyetat0_mod.flags \
          callkeys_mod.o \
          ioipsl_getin_p_mod.o \
          iostart.o \
          slab_ice_h.o \
          surfdat_h.o \
          tabfi_mod.o \
          tracer_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__convadj.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

convadj.done: \
          convadj.o \
          callkeys_mod.done \
          comcstfi_mod.done \
          tracer_h.done
	touch $(FCM_DONEDIR)/$@

convadj.o: \
          $(PPSRCDIR0__phys)/convadj.f \
          FFLAGS__phys__convadj.flags \
          callkeys_mod.o \
          comcstfi_mod.o \
          tracer_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__call_rings.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

call_rings.done: \
          call_rings.o \
          callkeys_mod.done \
          comcstfi_mod.done \
          comdiurn_h.done \
          comsaison_h.done \
          radcommon_h.done
	touch $(FCM_DONEDIR)/$@

call_rings.o: \
          $(PPSRCDIR0__phys)/call_rings.f90 \
          FFLAGS__phys__call_rings.flags \
          callkeys_mod.o \
          comcstfi_mod.o \
          comdiurn_h.o \
          comsaison_h.o \
          radcommon_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__dtridgl.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

dtridgl.done: \
          dtridgl.o
	touch $(FCM_DONEDIR)/$@

dtridgl.o: \
          $(PPSRCDIR0__phys)/dtridgl.f \
          FFLAGS__phys__dtridgl.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__suaer_corrk.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

suaer_corrk.done: \
          suaer_corrk.o \
          aerosol_mod.done \
          callkeys_mod.done \
          datafile_mod.done \
          radcommon_h.done \
          radinc_h.done
	touch $(FCM_DONEDIR)/$@

suaer_corrk.o: \
          $(PPSRCDIR0__phys)/suaer_corrk.f90 \
          FFLAGS__phys__suaer_corrk.flags \
          aerosol_mod.o \
          callkeys_mod.o \
          datafile_mod.o \
          radcommon_h.o \
          radinc_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__sfluxi.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

sfluxi.done: \
          sfluxi.o \
          comcstfi_mod.done \
          radcommon_h.done \
          radinc_h.done
	touch $(FCM_DONEDIR)/$@

sfluxi.o: \
          $(PPSRCDIR0__phys)/sfluxi.f \
          FFLAGS__phys__sfluxi.flags \
          comcstfi_mod.o \
          radcommon_h.o \
          radinc_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__inifis_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

inifis_mod.done: \
          inifis_mod.o \
          callkeys_mod.done \
          comcstfi_mod.done \
          comdiurn_h.done \
          comgeomfi_h.done \
          comsoil_h.done \
          datafile_mod.done \
          init_print_control_mod.done \
          ioipsl_getin_p_mod.done \
          mod_phys_lmdz_para.done \
          planete_mod.done \
          planetwide_mod.done \
          radinc_h.done \
          time_phylmdz_mod.done
	touch $(FCM_DONEDIR)/$@

inifis_mod.o: \
          $(PPSRCDIR0__phys)/inifis_mod.f90 \
          FFLAGS__phys__inifis_mod.flags \
          callkeys_mod.o \
          comcstfi_mod.o \
          comdiurn_h.o \
          comgeomfi_h.o \
          comsoil_h.o \
          datafile_mod.o \
          init_print_control_mod.o \
          ioipsl_getin_p_mod.o \
          mod_phys_lmdz_para.o \
          planete_mod.o \
          planetwide_mod.o \
          radinc_h.o \
          time_phylmdz_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__calc_rayleigh.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

calc_rayleigh.done: \
          calc_rayleigh.o \
          comcstfi_mod.done \
          gases_h.done \
          radcommon_h.done \
          radinc_h.done
	touch $(FCM_DONEDIR)/$@

calc_rayleigh.o: \
          $(PPSRCDIR0__phys)/calc_rayleigh.f90 \
          FFLAGS__phys__calc_rayleigh.flags \
          comcstfi_mod.o \
          gases_h.o \
          radcommon_h.o \
          radinc_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__rings.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

rings.done: \
          rings.o \
          comdiurn_h.done \
          geometry_mod.done
	touch $(FCM_DONEDIR)/$@

rings.o: \
          $(PPSRCDIR0__phys)/rings.f90 \
          FFLAGS__phys__rings.flags \
          comdiurn_h.o \
          geometry_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__albedo_snow_calc.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

albedo_snow_calc.done: \
          albedo_snow_calc.o
	touch $(FCM_DONEDIR)/$@

albedo_snow_calc.o: \
          $(PPSRCDIR0__phys)/albedo_snow_calc.f90 \
          FFLAGS__phys__albedo_snow_calc.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__gfluxi.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

gfluxi.done: \
          gfluxi.o \
          comcstfi_mod.done \
          radcommon_h.done \
          radinc_h.done
	touch $(FCM_DONEDIR)/$@

gfluxi.o: \
          $(PPSRCDIR0__phys)/gfluxi.f \
          FFLAGS__phys__gfluxi.flags \
          comcstfi_mod.o \
          radcommon_h.o \
          radinc_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__aeropacity.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

aeropacity.done: \
          aeropacity.o \
          aerosol_mod.done \
          callkeys_mod.done \
          comcstfi_mod.done \
          geometry_mod.done \
          radinc_h.done \
          tracer_h.done
	touch $(FCM_DONEDIR)/$@

aeropacity.o: \
          $(PPSRCDIR0__phys)/aeropacity.f90 \
          FFLAGS__phys__aeropacity.flags \
          aerosol_mod.o \
          callkeys_mod.o \
          comcstfi_mod.o \
          geometry_mod.o \
          radinc_h.o \
          tracer_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__comm_wrf.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

comm_wrf.done: \
          comm_wrf.o
	touch $(FCM_DONEDIR)/$@

comm_wrf.o: \
          $(PPSRCDIR0__phys)/comm_wrf.f90 \
          FFLAGS__phys__comm_wrf.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__dsolver.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

dsolver.done: \
          dsolver.o
	touch $(FCM_DONEDIR)/$@

dsolver.o: \
          $(PPSRCDIR0__phys)/dsolver.f \
          FFLAGS__phys__dsolver.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__slab_ice_h.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

slab_ice_h.done: \
          slab_ice_h.o
	touch $(FCM_DONEDIR)/$@

slab_ice_h.o: \
          $(PPSRCDIR0__phys)/slab_ice_h.f90 \
          FFLAGS__phys__slab_ice_h.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__dimphy.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

dimphy.done: \
          dimphy.o
	touch $(FCM_DONEDIR)/$@

dimphy.o: \
          $(PPSRCDIR0__phys)/dimphy.f90 \
          FFLAGS__phys__dimphy.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__params_h.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

params_h.done: \
          params_h.o
	touch $(FCM_DONEDIR)/$@

params_h.o: \
          $(PPSRCDIR0__phys)/params_h.f90 \
          FFLAGS__phys__params_h.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__gfluxv.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

gfluxv.done: \
          gfluxv.o \
          radinc_h.done
	touch $(FCM_DONEDIR)/$@

gfluxv.o: \
          $(PPSRCDIR0__phys)/gfluxv.f \
          FFLAGS__phys__gfluxv.flags \
          radinc_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__soil.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

soil.done: \
          soil.o \
          comcstfi_mod.done \
          comsoil_h.done \
          geometry_mod.done \
          planete_mod.done \
          time_phylmdz_mod.done
	touch $(FCM_DONEDIR)/$@

soil.o: \
          $(PPSRCDIR0__phys)/soil.f \
          FFLAGS__phys__soil.flags \
          comcstfi_mod.o \
          comsoil_h.o \
          geometry_mod.o \
          planete_mod.o \
          time_phylmdz_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__thermcell_flux.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

thermcell_flux.done: \
          thermcell_flux.o \
          print_control_mod.done \
          thermcell_mod.done
	touch $(FCM_DONEDIR)/$@

thermcell_flux.o: \
          $(PPSRCDIR0__phys)/thermcell_flux.f90 \
          FFLAGS__phys__thermcell_flux.flags \
          print_control_mod.o \
          thermcell_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__surfini.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

surfini.done: \
          surfini.o \
          callkeys_mod.done \
          planetwide_mod.done \
          radinc_h.done \
          surfdat_h.done \
          tracer_h.done
	touch $(FCM_DONEDIR)/$@

surfini.o: \
          $(PPSRCDIR0__phys)/surfini.f \
          FFLAGS__phys__surfini.flags \
          callkeys_mod.o \
          planetwide_mod.o \
          radinc_h.o \
          surfdat_h.o \
          tracer_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__xios_output_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

xios_output_mod.done: \
          xios_output_mod.o
	touch $(FCM_DONEDIR)/$@

xios_output_mod.o: \
          $(PPSRCDIR0__phys)/xios_output_mod.f90 \
          FFLAGS__phys__xios_output_mod.flags
	fcm_internal compile:F phys $< $@

scatterers.h: \
          $(SRCDIR0__phys)/scatterers.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

scatterers.h.idone: \
          $(SRCDIR0__phys)/scatterers.h
	touch $(FCM_DONEDIR)/$@

FFLAGS__phys__mufract.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

mufract.done: \
          mufract.o
	touch $(FCM_DONEDIR)/$@

mufract.o: \
          $(PPSRCDIR0__phys)/mufract.f \
          FFLAGS__phys__mufract.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__gases_h.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

gases_h.done: \
          gases_h.o
	touch $(FCM_DONEDIR)/$@

gases_h.o: \
          $(PPSRCDIR0__phys)/gases_h.f90 \
          FFLAGS__phys__gases_h.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__surface_nature.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

surface_nature.done: \
          surface_nature.o \
          comsoil_h.done \
          geometry_mod.done \
          surfdat_h.done \
          tracer_h.done
	touch $(FCM_DONEDIR)/$@

surface_nature.o: \
          $(PPSRCDIR0__phys)/surface_nature.f \
          FFLAGS__phys__surface_nature.flags \
          comsoil_h.o \
          geometry_mod.o \
          surfdat_h.o \
          tracer_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__su_gases.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

su_gases.done: \
          su_gases.o \
          gases_h.done
	touch $(FCM_DONEDIR)/$@

su_gases.o: \
          $(PPSRCDIR0__phys)/su_gases.f90 \
          FFLAGS__phys__su_gases.flags \
          gases_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__thermcell_closure.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

thermcell_closure.done: \
          thermcell_closure.o \
          thermcell_mod.done
	touch $(FCM_DONEDIR)/$@

thermcell_closure.o: \
          $(PPSRCDIR0__phys)/thermcell_closure.f90 \
          FFLAGS__phys__thermcell_closure.flags \
          thermcell_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__rain.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

rain.done: \
          rain.o \
          comcstfi_mod.done \
          ioipsl_getin_p_mod.done \
          radii_mod.done \
          tracer_h.done \
          watercommon_h.done
	touch $(FCM_DONEDIR)/$@

rain.o: \
          $(PPSRCDIR0__phys)/rain.f90 \
          FFLAGS__phys__rain.flags \
          comcstfi_mod.o \
          ioipsl_getin_p_mod.o \
          radii_mod.o \
          tracer_h.o \
          watercommon_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__totalcloudfrac.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

totalcloudfrac.done: \
          totalcloudfrac.o \
          callkeys_mod.done \
          comdiurn_h.done \
          tracer_h.done \
          watercommon_h.done
	touch $(FCM_DONEDIR)/$@

totalcloudfrac.o: \
          $(PPSRCDIR0__phys)/totalcloudfrac.f90 \
          FFLAGS__phys__totalcloudfrac.flags \
          callkeys_mod.o \
          comdiurn_h.o \
          tracer_h.o \
          watercommon_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__optcv.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

optcv_mod.done: \
          optcv_mod.o \
          callkeys_mod.done \
          comcstfi_mod.done \
          gases_h.done \
          radcommon_h.done \
          radinc_h.done
	touch $(FCM_DONEDIR)/$@

optcv_mod.o: \
          $(PPSRCDIR0__phys)/optcv.f90 \
          FFLAGS__phys__optcv.flags \
          callkeys_mod.o \
          comcstfi_mod.o \
          gases_h.o \
          radcommon_h.o \
          radinc_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__aeroptproperties.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

aeroptproperties.done: \
          aeroptproperties.o \
          radcommon_h.done \
          radinc_h.done
	touch $(FCM_DONEDIR)/$@

aeroptproperties.o: \
          $(PPSRCDIR0__phys)/aeroptproperties.f90 \
          FFLAGS__phys__aeroptproperties.flags \
          radcommon_h.o \
          radinc_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__thermcell_plume.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

thermcell_plume.done: \
          thermcell_plume.o \
          print_control_mod.done \
          thermcell_mod.done \
          tracer_h.done \
          watercommon_h.done
	touch $(FCM_DONEDIR)/$@

thermcell_plume.o: \
          $(PPSRCDIR0__phys)/thermcell_plume.f90 \
          FFLAGS__phys__thermcell_plume.flags \
          print_control_mod.o \
          thermcell_mod.o \
          tracer_h.o \
          watercommon_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__thermcell_main.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

thermcell_main.done: \
          thermcell_main.o \
          print_control_mod.done \
          thermcell_mod.done
	touch $(FCM_DONEDIR)/$@

thermcell_main.o: \
          $(PPSRCDIR0__phys)/thermcell_main.f90 \
          FFLAGS__phys__thermcell_main.flags \
          print_control_mod.o \
          thermcell_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__mucorr.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

mucorr.done: \
          mucorr.o
	touch $(FCM_DONEDIR)/$@

mucorr.o: \
          $(PPSRCDIR0__phys)/mucorr.f \
          FFLAGS__phys__mucorr.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__watercommon_h.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

watercommon_h.done: \
          watercommon_h.o \
          comcstfi_mod.done
	touch $(FCM_DONEDIR)/$@

watercommon_h.o: \
          $(PPSRCDIR0__phys)/watercommon_h.f90 \
          FFLAGS__phys__watercommon_h.flags \
          comcstfi_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__iniaerosol.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

iniaerosol.done: \
          iniaerosol.o \
          aerosol_mod.done \
          callkeys_mod.done \
          radinc_h.done
	touch $(FCM_DONEDIR)/$@

iniaerosol.o: \
          $(PPSRCDIR0__phys)/iniaerosol.f \
          FFLAGS__phys__iniaerosol.flags \
          aerosol_mod.o \
          callkeys_mod.o \
          radinc_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__comdiurn_h.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

comdiurn_h.done: \
          comdiurn_h.o
	touch $(FCM_DONEDIR)/$@

comdiurn_h.o: \
          $(PPSRCDIR0__phys)/comdiurn_h.f90 \
          FFLAGS__phys__comdiurn_h.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__iniorbit.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

iniorbit.done: \
          iniorbit.o \
          comcstfi_mod.done \
          planete_mod.done
	touch $(FCM_DONEDIR)/$@

iniorbit.o: \
          $(PPSRCDIR0__phys)/iniorbit.f \
          FFLAGS__phys__iniorbit.flags \
          comcstfi_mod.o \
          planete_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__optci.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

optci_mod.done: \
          optci_mod.o \
          callkeys_mod.done \
          comcstfi_mod.done \
          gases_h.done \
          radcommon_h.done \
          radinc_h.done
	touch $(FCM_DONEDIR)/$@

optci_mod.o: \
          $(PPSRCDIR0__phys)/optci.f90 \
          FFLAGS__phys__optci.flags \
          callkeys_mod.o \
          comcstfi_mod.o \
          gases_h.o \
          radcommon_h.o \
          radinc_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__interp_line.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

interp_line.done: \
          interp_line.o
	touch $(FCM_DONEDIR)/$@

interp_line.o: \
          $(PPSRCDIR0__phys)/interp_line.f \
          FFLAGS__phys__interp_line.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__largescale.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

largescale.done: \
          largescale.o \
          ioipsl_getin_p_mod.done \
          tracer_h.done \
          watercommon_h.done
	touch $(FCM_DONEDIR)/$@

largescale.o: \
          $(PPSRCDIR0__phys)/largescale.f90 \
          FFLAGS__phys__largescale.flags \
          ioipsl_getin_p_mod.o \
          tracer_h.o \
          watercommon_h.o
	fcm_internal compile:F phys $< $@

