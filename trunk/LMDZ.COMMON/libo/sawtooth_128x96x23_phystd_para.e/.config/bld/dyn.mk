# Automatic Make rule for dyn

SRCDIR0__dyn = /home/palatyle/LMD_gen/trunk/LMDZ.COMMON/libf/dyn3dpar

PPSRCDIR0__dyn = /home/palatyle/LMD_gen/trunk/LMDZ.COMMON/libo/sawtooth_128x96x23_phystd_para.e/.config/ppsrc/dyn

FFLAGS__dyn__initfluxsto_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

initfluxsto_p.done: \
          initfluxsto_p.o \
          write_field.done \
          comconst_mod.done \
          comvert_mod.done \
          misc_mod.done \
          parallel_lmdz.done \
          temps_mod.done
	touch $(FCM_DONEDIR)/$@

initfluxsto_p.o: \
          $(PPSRCDIR0__dyn)/initfluxsto_p.f \
          FFLAGS__dyn__initfluxsto_p.flags \
          write_field.o \
          comconst_mod.o \
          comvert_mod.o \
          misc_mod.o \
          parallel_lmdz.o \
          temps_mod.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__sw_case_williamson91_6.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

sw_case_williamson91_6.done: \
          sw_case_williamson91_6.o \
          comconst_mod.done \
          comvert_mod.done
	touch $(FCM_DONEDIR)/$@

sw_case_williamson91_6.o: \
          $(PPSRCDIR0__dyn)/sw_case_williamson91_6.f \
          FFLAGS__dyn__sw_case_williamson91_6.flags \
          comconst_mod.o \
          comvert_mod.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__tourpot_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

tourpot_p.done: \
          tourpot_p.o \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

tourpot_p.o: \
          $(PPSRCDIR0__dyn)/tourpot_p.f \
          FFLAGS__dyn__tourpot_p.flags \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__vlsplt_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

vlsplt_p.done: \
          vlsplt_p.o \
          parallel_lmdz.done \
          vampir.done \
          comconst_mod.done \
          mod_hallo.done \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

vlsplt_p.o: \
          $(PPSRCDIR0__dyn)/vlsplt_p.f \
          FFLAGS__dyn__vlsplt_p.flags \
          parallel_lmdz.o \
          vampir.o \
          comconst_mod.o \
          mod_hallo.o \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__laplacien_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

laplacien_p.done: \
          laplacien_p.o \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

laplacien_p.o: \
          $(PPSRCDIR0__dyn)/laplacien_p.f \
          FFLAGS__dyn__laplacien_p.flags \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__laplacien_rotgam_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

laplacien_rotgam_p.done: \
          laplacien_rotgam_p.o \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

laplacien_rotgam_p.o: \
          $(PPSRCDIR0__dyn)/laplacien_rotgam_p.f \
          FFLAGS__dyn__laplacien_rotgam_p.flags \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__writehist_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

writehist_p.done: \
          writehist_p.o \
          infotrac.done \
          misc_mod.done \
          parallel_lmdz.done \
          temps_mod.done
	touch $(FCM_DONEDIR)/$@

writehist_p.o: \
          $(PPSRCDIR0__dyn)/writehist_p.f \
          FFLAGS__dyn__writehist_p.flags \
          infotrac.o \
          misc_mod.o \
          parallel_lmdz.o \
          temps_mod.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__groupe_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

groupe_p.done: \
          groupe_p.o \
          comconst_mod.done \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

groupe_p.o: \
          $(PPSRCDIR0__dyn)/groupe_p.f \
          FFLAGS__dyn__groupe_p.flags \
          comconst_mod.o \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__dissip_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

dissip_p.done: \
          dissip_p.o \
          comconst_mod.done \
          parallel_lmdz.done \
          write_field_p.done
	touch $(FCM_DONEDIR)/$@

dissip_p.o: \
          $(PPSRCDIR0__dyn)/dissip_p.f \
          FFLAGS__dyn__dissip_p.flags \
          comconst_mod.o \
          parallel_lmdz.o \
          write_field_p.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__bernoui_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

bernoui_p.done: \
          bernoui_p.o \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

bernoui_p.o: \
          $(PPSRCDIR0__dyn)/bernoui_p.f \
          FFLAGS__dyn__bernoui_p.flags \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__bilan_dyn_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

bilan_dyn_p.done: \
          bilan_dyn_p.o \
          comconst_mod.done \
          comvert_mod.done \
          misc_mod.done \
          mod_hallo.done \
          parallel_lmdz.done \
          temps_mod.done \
          write_field_p.done
	touch $(FCM_DONEDIR)/$@

bilan_dyn_p.o: \
          $(PPSRCDIR0__dyn)/bilan_dyn_p.f \
          FFLAGS__dyn__bilan_dyn_p.flags \
          comconst_mod.o \
          comvert_mod.o \
          misc_mod.o \
          mod_hallo.o \
          parallel_lmdz.o \
          temps_mod.o \
          write_field_p.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__dteta1_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

dteta1_p.done: \
          dteta1_p.o \
          parallel_lmdz.done \
          write_field_p.done
	touch $(FCM_DONEDIR)/$@

dteta1_p.o: \
          $(PPSRCDIR0__dyn)/dteta1_p.f \
          FFLAGS__dyn__dteta1_p.flags \
          parallel_lmdz.o \
          write_field_p.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__gr_u_scal_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

gr_u_scal_p.done: \
          gr_u_scal_p.o \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

gr_u_scal_p.o: \
          $(PPSRCDIR0__dyn)/gr_u_scal_p.f \
          FFLAGS__dyn__gr_u_scal_p.flags \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__dynredem_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

dynredem0_p.done: \
          dynredem0_p.o \
          comconst_mod.done \
          comgeom.h.idone \
          comgeom2.h.idone \
          comvert_mod.done \
          control_mod.done \
          dimensions.h.idone \
          dynredem_mod.done \
          ener_mod.done \
          infotrac.done \
          iniprint.h.idone \
          logic_mod.done \
          netcdf95.done \
          parallel_lmdz.done \
          paramet.h.idone \
          serre_mod.done \
          temps_mod.done
	touch $(FCM_DONEDIR)/$@

dynredem0_p.o: \
          $(PPSRCDIR0__dyn)/dynredem_p.f90 \
          FFLAGS__dyn__dynredem_p.flags \
          comconst_mod.o \
          comgeom.h \
          comgeom2.h \
          comvert_mod.o \
          control_mod.o \
          dimensions.h \
          dynredem_mod.o \
          ener_mod.o \
          infotrac.o \
          iniprint.h \
          logic_mod.o \
          netcdf95.o \
          parallel_lmdz.o \
          paramet.h \
          serre_mod.o \
          temps_mod.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__diverg_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

diverg_p.done: \
          diverg_p.o \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

diverg_p.o: \
          $(PPSRCDIR0__dyn)/diverg_p.f \
          FFLAGS__dyn__diverg_p.flags \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__rotat_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

rotat_p.done: \
          rotat_p.o \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

rotat_p.o: \
          $(PPSRCDIR0__dyn)/rotat_p.f \
          FFLAGS__dyn__rotat_p.flags \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__wrgrads.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

wrgrads.done: \
          wrgrads.o
	touch $(FCM_DONEDIR)/$@

wrgrads.o: \
          $(PPSRCDIR0__dyn)/wrgrads.f \
          FFLAGS__dyn__wrgrads.flags
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__gcm.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

LDFLAGS__dyn__gcm.flags: \
          LDFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

LD__dyn__gcm.flags: \
          LD__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

gcm.o: \
          $(PPSRCDIR0__dyn)/gcm.f \
          FFLAGS__dyn__gcm.flags \
          bands.o \
          comconst_mod.o \
          control_mod.o \
          cpdet_mod.o \
          dimphy.o \
          filtreg_mod.o \
          getparam.o \
          infotrac.o \
          iniphysiq_mod.o \
          logic_mod.o \
          mod_const_mpi.o \
          mod_grid_phy_lmdz.o \
          mod_hallo.o \
          mod_phys_lmdz_omp_data.o \
          parallel_lmdz.o \
          temps_mod.o
	fcm_internal compile:F dyn $< $@

gcm_128x96x23_phystd_para.e: \
          gcm.o \
          LD__dyn__gcm.flags \
          LDFLAGS__dyn__gcm.flags \
          $(OBJECTS) \
          bands.done \
          comconst_mod.done \
          control_mod.done \
          cpdet_mod.done \
          dimphy.done \
          filtreg_mod.done \
          getparam.done \
          infotrac.done \
          iniphysiq_mod.done \
          logic_mod.done \
          mod_const_mpi.done \
          mod_grid_phy_lmdz.done \
          mod_hallo.done \
          mod_phys_lmdz_omp_data.done \
          parallel_lmdz.done \
          temps_mod.done
	fcm_internal load dyn $< $@

FFLAGS__dyn__geopot_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

geopot_p.done: \
          geopot_p.o \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

geopot_p.o: \
          $(PPSRCDIR0__dyn)/geopot_p.f \
          FFLAGS__dyn__geopot_p.flags \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__massdair_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

massdair_p.done: \
          massdair_p.o \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

massdair_p.o: \
          $(PPSRCDIR0__dyn)/massdair_p.f \
          FFLAGS__dyn__massdair_p.flags \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__groupeun_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

groupeun_p.done: \
          groupeun_p.o \
          write_field_p.done \
          comconst_mod.done \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

groupeun_p.o: \
          $(PPSRCDIR0__dyn)/groupeun_p.f \
          FFLAGS__dyn__groupeun_p.flags \
          write_field_p.o \
          comconst_mod.o \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__temps_mod.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

temps_mod.done: \
          temps_mod.o
	touch $(FCM_DONEDIR)/$@

temps_mod.o: \
          $(PPSRCDIR0__dyn)/temps_mod.f90 \
          FFLAGS__dyn__temps_mod.flags
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__bands.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

bands.done: \
          bands.o \
          dimensions.h.idone \
          mod_phys_lmdz_para.done \
          parallel_lmdz.done \
          times.done
	touch $(FCM_DONEDIR)/$@

bands.o: \
          $(PPSRCDIR0__dyn)/bands.f90 \
          FFLAGS__dyn__bands.flags \
          dimensions.h \
          mod_phys_lmdz_para.o \
          parallel_lmdz.o \
          times.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__nxgrad_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

nxgrad_p.done: \
          nxgrad_p.o \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

nxgrad_p.o: \
          $(PPSRCDIR0__dyn)/nxgrad_p.f \
          FFLAGS__dyn__nxgrad_p.flags \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__vlspltqs_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

vlspltqs_p.done: \
          vlspltqs_p.o \
          vampir.done \
          comconst_mod.done \
          cpdet_mod.done \
          mod_hallo.done \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

vlspltqs_p.o: \
          $(PPSRCDIR0__dyn)/vlspltqs_p.f \
          FFLAGS__dyn__vlspltqs_p.flags \
          vampir.o \
          comconst_mod.o \
          cpdet_mod.o \
          mod_hallo.o \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__inithist_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

inithist_p.done: \
          inithist_p.o \
          write_field.done \
          comconst_mod.done \
          comvert_mod.done \
          infotrac.done \
          misc_mod.done \
          parallel_lmdz.done \
          temps_mod.done
	touch $(FCM_DONEDIR)/$@

inithist_p.o: \
          $(PPSRCDIR0__dyn)/inithist_p.f \
          FFLAGS__dyn__inithist_p.flags \
          write_field.o \
          comconst_mod.o \
          comvert_mod.o \
          infotrac.o \
          misc_mod.o \
          parallel_lmdz.o \
          temps_mod.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__dudv2_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

dudv2_p.done: \
          dudv2_p.o \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

dudv2_p.o: \
          $(PPSRCDIR0__dyn)/dudv2_p.f \
          FFLAGS__dyn__dudv2_p.flags \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

omp_chunk.h: \
          $(SRCDIR0__dyn)/omp_chunk.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

omp_chunk.h.idone: \
          $(SRCDIR0__dyn)/omp_chunk.h
	touch $(FCM_DONEDIR)/$@

FFLAGS__dyn__qminimum_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

qminimum_p.done: \
          qminimum_p.o \
          comvert_mod.done \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

qminimum_p.o: \
          $(PPSRCDIR0__dyn)/qminimum_p.f \
          FFLAGS__dyn__qminimum_p.flags \
          comvert_mod.o \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__covcont_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

covcont_p.done: \
          covcont_p.o \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

covcont_p.o: \
          $(PPSRCDIR0__dyn)/covcont_p.f \
          FFLAGS__dyn__covcont_p.flags \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__planetary_operations_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

planetary_operations_p.done: \
          planetary_operations_p.o \
          comgeom.h.idone \
          dimensions.h.idone \
          mod_const_mpi.done \
          parallel_lmdz.done \
          paramet.h.idone
	touch $(FCM_DONEDIR)/$@

planetary_operations_p.o: \
          $(PPSRCDIR0__dyn)/planetary_operations_p.f90 \
          FFLAGS__dyn__planetary_operations_p.flags \
          comgeom.h \
          dimensions.h \
          mod_const_mpi.o \
          parallel_lmdz.o \
          paramet.h
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__massbarxy_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

massbarxy_p.done: \
          massbarxy_p.o \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

massbarxy_p.o: \
          $(PPSRCDIR0__dyn)/massbarxy_p.f \
          FFLAGS__dyn__massbarxy_p.flags \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__mod_hallo.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

mod_hallo.done: \
          mod_hallo.o \
          mod_const_mpi.done \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

mod_hallo.o: \
          $(PPSRCDIR0__dyn)/mod_hallo.f90 \
          FFLAGS__dyn__mod_hallo.flags \
          mod_const_mpi.o \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__advect_new_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

advect_new_p.done: \
          advect_new_p.o \
          comconst_mod.done \
          comgeom.h.idone \
          dimensions.h.idone \
          logic_mod.done \
          parallel_lmdz.done \
          paramet.h.idone \
          write_field_p.done
	touch $(FCM_DONEDIR)/$@

advect_new_p.o: \
          $(PPSRCDIR0__dyn)/advect_new_p.f \
          FFLAGS__dyn__advect_new_p.flags \
          comconst_mod.o \
          comgeom.h \
          dimensions.h \
          logic_mod.o \
          parallel_lmdz.o \
          paramet.h \
          write_field_p.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__laplacien_rot_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

laplacien_rot_p.done: \
          laplacien_rot_p.o \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

laplacien_rot_p.o: \
          $(PPSRCDIR0__dyn)/laplacien_rot_p.f \
          FFLAGS__dyn__laplacien_rot_p.flags \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__sponge_mod_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

sponge_mod_p.done: \
          sponge_mod_p.o \
          write_field_p.done \
          comvert_mod.done \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

sponge_mod_p.o: \
          $(PPSRCDIR0__dyn)/sponge_mod_p.f90 \
          FFLAGS__dyn__sponge_mod_p.flags \
          write_field_p.o \
          comvert_mod.o \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__leapfrog_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

leapfrog_p.done: \
          leapfrog_p.o \
          bands.done \
          write_field.done \
          write_field_p.done \
          comconst_mod.done \
          comuforc_h.done \
          comvert_mod.done \
          control_mod.done \
          cpdet_mod.done \
          exner_hyb_m.done \
          exner_hyb_p_m.done \
          exner_milieu_m.done \
          exner_milieu_p_m.done \
          getparam.done \
          guide_p_mod.done \
          infotrac.done \
          logic_mod.done \
          misc_mod.done \
          mod_hallo.done \
          parallel_lmdz.done \
          sponge_mod_p.done \
          temps_mod.done \
          timer_filtre.done \
          times.done \
          vampir.done
	touch $(FCM_DONEDIR)/$@

leapfrog_p.o: \
          $(PPSRCDIR0__dyn)/leapfrog_p.f \
          FFLAGS__dyn__leapfrog_p.flags \
          bands.o \
          write_field.o \
          write_field_p.o \
          comconst_mod.o \
          comuforc_h.o \
          comvert_mod.o \
          control_mod.o \
          cpdet_mod.o \
          exner_hyb_m.o \
          exner_hyb_p_m.o \
          exner_milieu_m.o \
          exner_milieu_p_m.o \
          getparam.o \
          guide_p_mod.o \
          infotrac.o \
          logic_mod.o \
          misc_mod.o \
          mod_hallo.o \
          parallel_lmdz.o \
          sponge_mod_p.o \
          temps_mod.o \
          timer_filtre.o \
          times.o \
          vampir.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__convmas_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

convmas_p.done: \
          convmas_p.o \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

convmas_p.o: \
          $(PPSRCDIR0__dyn)/convmas_p.f \
          FFLAGS__dyn__convmas_p.flags \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__exner_hyb_p_m.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

exner_hyb_p_m.done: \
          exner_hyb_p_m.o \
          comconst_mod.done \
          comgeom.h.idone \
          comvert_mod.done \
          dimensions.h.idone \
          parallel_lmdz.done \
          paramet.h.idone
	touch $(FCM_DONEDIR)/$@

exner_hyb_p_m.o: \
          $(PPSRCDIR0__dyn)/exner_hyb_p_m.f90 \
          FFLAGS__dyn__exner_hyb_p_m.flags \
          comconst_mod.o \
          comgeom.h \
          comvert_mod.o \
          dimensions.h \
          parallel_lmdz.o \
          paramet.h
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__pression_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

pression_p.done: \
          pression_p.o \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

pression_p.o: \
          $(PPSRCDIR0__dyn)/pression_p.f \
          FFLAGS__dyn__pression_p.flags \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__grad_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

grad_p.done: \
          grad_p.o \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

grad_p.o: \
          $(PPSRCDIR0__dyn)/grad_p.f \
          FFLAGS__dyn__grad_p.flags \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__massbar_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

massbar_p.done: \
          massbar_p.o \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

massbar_p.o: \
          $(PPSRCDIR0__dyn)/massbar_p.f \
          FFLAGS__dyn__massbar_p.flags \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__conf_gcm.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

conf_gcm.done: \
          conf_gcm.o \
          assert_m.done \
          comconst_mod.done \
          comdissnew.h.idone \
          control_mod.done \
          dimensions.h.idone \
          infotrac.done \
          iniprint.h.idone \
          ioipsl_getincom.done \
          logic_mod.done \
          misc_mod.done \
          mod_filtre_fft.done \
          mod_hallo.done \
          paramet.h.idone \
          serre_mod.done \
          sponge_mod_p.done \
          temps_mod.done
	touch $(FCM_DONEDIR)/$@

conf_gcm.o: \
          $(PPSRCDIR0__dyn)/conf_gcm.f90 \
          FFLAGS__dyn__conf_gcm.flags \
          assert_m.o \
          comconst_mod.o \
          comdissnew.h \
          control_mod.o \
          dimensions.h \
          infotrac.o \
          iniprint.h \
          ioipsl_getincom.o \
          logic_mod.o \
          misc_mod.o \
          mod_filtre_fft.o \
          mod_hallo.o \
          paramet.h \
          serre_mod.o \
          sponge_mod_p.o \
          temps_mod.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__divgrad_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

divgrad_p.done: \
          divgrad_p.o \
          parallel_lmdz.done \
          times.done
	touch $(FCM_DONEDIR)/$@

divgrad_p.o: \
          $(PPSRCDIR0__dyn)/divgrad_p.f \
          FFLAGS__dyn__divgrad_p.flags \
          parallel_lmdz.o \
          times.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__addfi_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

addfi_p.done: \
          addfi_p.o \
          comconst_mod.done \
          control_mod.done \
          infotrac.done \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

addfi_p.o: \
          $(PPSRCDIR0__dyn)/addfi_p.f \
          FFLAGS__dyn__addfi_p.flags \
          comconst_mod.o \
          control_mod.o \
          infotrac.o \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__gr_v_scal_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

gr_v_scal_p.done: \
          gr_v_scal_p.o \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

gr_v_scal_p.o: \
          $(PPSRCDIR0__dyn)/gr_v_scal_p.f \
          FFLAGS__dyn__gr_v_scal_p.flags \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__top_bound_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

top_bound_p.done: \
          top_bound_p.o \
          comconst_mod.done \
          comvert_mod.done \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

top_bound_p.o: \
          $(PPSRCDIR0__dyn)/top_bound_p.f \
          FFLAGS__dyn__top_bound_p.flags \
          comconst_mod.o \
          comvert_mod.o \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__guide_p_mod.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

guide_p_mod.done: \
          guide_p_mod.o \
          bands.done \
          write_field_p.done \
          comconst_mod.done \
          comgeom.h.idone \
          comgeom2.h.idone \
          comvert_mod.done \
          control_mod.done \
          dimensions.h.idone \
          exner_hyb_p_m.done \
          exner_milieu_p_m.done \
          getparam.done \
          mod_hallo.done \
          netcdf95.done \
          parallel_lmdz.done \
          paramet.h.idone \
          pres2lev_mod.done \
          serre_mod.done
	touch $(FCM_DONEDIR)/$@

guide_p_mod.o: \
          $(PPSRCDIR0__dyn)/guide_p_mod.f90 \
          FFLAGS__dyn__guide_p_mod.flags \
          bands.o \
          write_field_p.o \
          comconst_mod.o \
          comgeom.h \
          comgeom2.h \
          comvert_mod.o \
          control_mod.o \
          dimensions.h \
          exner_hyb_p_m.o \
          exner_milieu_p_m.o \
          getparam.o \
          mod_hallo.o \
          netcdf95.o \
          parallel_lmdz.o \
          paramet.h \
          pres2lev_mod.o \
          serre_mod.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__enercin_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

enercin_p.done: \
          enercin_p.o \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

enercin_p.o: \
          $(PPSRCDIR0__dyn)/enercin_p.f \
          FFLAGS__dyn__enercin_p.flags \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__gradiv_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

gradiv_p.done: \
          gradiv_p.o \
          parallel_lmdz.done \
          times.done
	touch $(FCM_DONEDIR)/$@

gradiv_p.o: \
          $(PPSRCDIR0__dyn)/gradiv_p.f \
          FFLAGS__dyn__gradiv_p.flags \
          parallel_lmdz.o \
          times.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__exner_milieu_p_m.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

exner_milieu_p_m.done: \
          exner_milieu_p_m.o \
          comconst_mod.done \
          comgeom.h.idone \
          comvert_mod.done \
          dimensions.h.idone \
          parallel_lmdz.done \
          paramet.h.idone
	touch $(FCM_DONEDIR)/$@

exner_milieu_p_m.o: \
          $(PPSRCDIR0__dyn)/exner_milieu_p_m.f90 \
          FFLAGS__dyn__exner_milieu_p_m.flags \
          comconst_mod.o \
          comgeom.h \
          comvert_mod.o \
          dimensions.h \
          parallel_lmdz.o \
          paramet.h
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__initdynav_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

initdynav_p.done: \
          initdynav_p.o \
          write_field.done \
          comconst_mod.done \
          comvert_mod.done \
          infotrac.done \
          misc_mod.done \
          parallel_lmdz.done \
          temps_mod.done
	touch $(FCM_DONEDIR)/$@

initdynav_p.o: \
          $(PPSRCDIR0__dyn)/initdynav_p.f \
          FFLAGS__dyn__initdynav_p.flags \
          write_field.o \
          comconst_mod.o \
          comvert_mod.o \
          infotrac.o \
          misc_mod.o \
          parallel_lmdz.o \
          temps_mod.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__advect_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

advect_p.done: \
          advect_p.o \
          comconst_mod.done \
          logic_mod.done \
          parallel_lmdz.done \
          write_field_p.done
	touch $(FCM_DONEDIR)/$@

advect_p.o: \
          $(PPSRCDIR0__dyn)/advect_p.f \
          FFLAGS__dyn__advect_p.flags \
          comconst_mod.o \
          logic_mod.o \
          parallel_lmdz.o \
          write_field_p.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__gradiv2_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

gradiv2_p.done: \
          gradiv2_p.o \
          write_field_p.done \
          mod_hallo.done \
          parallel_lmdz.done \
          times.done
	touch $(FCM_DONEDIR)/$@

gradiv2_p.o: \
          $(PPSRCDIR0__dyn)/gradiv2_p.f \
          FFLAGS__dyn__gradiv2_p.flags \
          write_field_p.o \
          mod_hallo.o \
          parallel_lmdz.o \
          times.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__fluxstokenc_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

fluxstokenc_p.done: \
          fluxstokenc_p.o
	touch $(FCM_DONEDIR)/$@

fluxstokenc_p.o: \
          $(PPSRCDIR0__dyn)/fluxstokenc_p.f \
          FFLAGS__dyn__fluxstokenc_p.flags
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__mod_const_mpi.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

mod_const_mpi.done: \
          mod_const_mpi.o \
          ioipsl_getincom.done
	touch $(FCM_DONEDIR)/$@

mod_const_mpi.o: \
          $(PPSRCDIR0__dyn)/mod_const_mpi.f90 \
          FFLAGS__dyn__mod_const_mpi.flags \
          ioipsl_getincom.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__convmas1_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

convmas1_p.done: \
          convmas1_p.o \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

convmas1_p.o: \
          $(PPSRCDIR0__dyn)/convmas1_p.f \
          FFLAGS__dyn__convmas1_p.flags \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__write_field_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

write_field_p.done: \
          write_field_p.o \
          parallel_lmdz.done \
          write_field.done
	touch $(FCM_DONEDIR)/$@

write_field_p.o: \
          $(PPSRCDIR0__dyn)/write_field_p.f90 \
          FFLAGS__dyn__write_field_p.flags \
          parallel_lmdz.o \
          write_field.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__nxgrarot_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

nxgrarot_p.done: \
          nxgrarot_p.o \
          parallel_lmdz.done \
          times.done \
          write_field_p.done
	touch $(FCM_DONEDIR)/$@

nxgrarot_p.o: \
          $(PPSRCDIR0__dyn)/nxgrarot_p.f \
          FFLAGS__dyn__nxgrarot_p.flags \
          parallel_lmdz.o \
          times.o \
          write_field_p.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__advtrac_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

advtrac_p.done: \
          advtrac_p.o \
          bands.done \
          vampir.done \
          write_field_p.done \
          comconst_mod.done \
          comdissip.h.idone \
          comgeom2.h.idone \
          control_mod.done \
          dimensions.h.idone \
          infotrac.done \
          mod_hallo.done \
          parallel_lmdz.done \
          paramet.h.idone \
          planetary_operations_p.done \
          times.done
	touch $(FCM_DONEDIR)/$@

advtrac_p.o: \
          $(PPSRCDIR0__dyn)/advtrac_p.f90 \
          FFLAGS__dyn__advtrac_p.flags \
          bands.o \
          vampir.o \
          write_field_p.o \
          comconst_mod.o \
          comdissip.h \
          comgeom2.h \
          control_mod.o \
          dimensions.h \
          infotrac.o \
          mod_hallo.o \
          parallel_lmdz.o \
          paramet.h \
          planetary_operations_p.o \
          times.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__covnat_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

covnat_p.done: \
          covnat_p.o \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

covnat_p.o: \
          $(PPSRCDIR0__dyn)/covnat_p.f \
          FFLAGS__dyn__covnat_p.flags \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__friction_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

friction_p.done: \
          friction_p.o \
          comconst_mod.done \
          control_mod.done \
          ioipsl_getincom.done \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

friction_p.o: \
          $(PPSRCDIR0__dyn)/friction_p.f \
          FFLAGS__dyn__friction_p.flags \
          comconst_mod.o \
          control_mod.o \
          ioipsl_getincom.o \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__diverg_gam_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

diverg_gam_p.done: \
          diverg_gam_p.o \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

diverg_gam_p.o: \
          $(PPSRCDIR0__dyn)/diverg_gam_p.f \
          FFLAGS__dyn__diverg_gam_p.flags \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__convflu_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

convflu_p.done: \
          convflu_p.o \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

convflu_p.o: \
          $(PPSRCDIR0__dyn)/convflu_p.f \
          FFLAGS__dyn__convflu_p.flags \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__logic_mod.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

logic_mod.done: \
          logic_mod.o
	touch $(FCM_DONEDIR)/$@

logic_mod.o: \
          $(PPSRCDIR0__dyn)/logic_mod.f90 \
          FFLAGS__dyn__logic_mod.flags
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__flumass_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

flumass_p.done: \
          flumass_p.o \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

flumass_p.o: \
          $(PPSRCDIR0__dyn)/flumass_p.f \
          FFLAGS__dyn__flumass_p.flags \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__vlspltgen_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

vlspltgen_p.done: \
          vlspltgen_p.o \
          vampir.done \
          write_field_p.done \
          comconst_mod.done \
          infotrac.done \
          mod_hallo.done \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

vlspltgen_p.o: \
          $(PPSRCDIR0__dyn)/vlspltgen_p.f \
          FFLAGS__dyn__vlspltgen_p.flags \
          vampir.o \
          write_field_p.o \
          comconst_mod.o \
          infotrac.o \
          mod_hallo.o \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__dudv1_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

dudv1_p.done: \
          dudv1_p.o \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

dudv1_p.o: \
          $(PPSRCDIR0__dyn)/dudv1_p.f \
          FFLAGS__dyn__dudv1_p.flags \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__times.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

times.done: \
          times.o \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

times.o: \
          $(PPSRCDIR0__dyn)/times.f90 \
          FFLAGS__dyn__times.flags \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__vitvert_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

vitvert_p.done: \
          vitvert_p.o \
          comvert_mod.done \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

vitvert_p.o: \
          $(PPSRCDIR0__dyn)/vitvert_p.f \
          FFLAGS__dyn__vitvert_p.flags \
          comvert_mod.o \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__rotatf_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

rotatf_p.done: \
          rotatf_p.o \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

rotatf_p.o: \
          $(PPSRCDIR0__dyn)/rotatf_p.f \
          FFLAGS__dyn__rotatf_p.flags \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__divergf_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

divergf_p.done: \
          divergf_p.o \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

divergf_p.o: \
          $(PPSRCDIR0__dyn)/divergf_p.f \
          FFLAGS__dyn__divergf_p.flags \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__abort_gcm.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

abort_gcm.done: \
          abort_gcm.o \
          ioipsl_getincom.done \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

abort_gcm.o: \
          $(PPSRCDIR0__dyn)/abort_gcm.f \
          FFLAGS__dyn__abort_gcm.flags \
          ioipsl_getincom.o \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__rotat_nfil_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

rotat_nfil_p.done: \
          rotat_nfil_p.o \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

rotat_nfil_p.o: \
          $(PPSRCDIR0__dyn)/rotat_nfil_p.f \
          FFLAGS__dyn__rotat_nfil_p.flags \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__writedynav_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

writedynav_p.done: \
          writedynav_p.o \
          comconst_mod.done \
          infotrac.done \
          misc_mod.done \
          parallel_lmdz.done \
          temps_mod.done
	touch $(FCM_DONEDIR)/$@

writedynav_p.o: \
          $(PPSRCDIR0__dyn)/writedynav_p.f \
          FFLAGS__dyn__writedynav_p.flags \
          comconst_mod.o \
          infotrac.o \
          misc_mod.o \
          parallel_lmdz.o \
          temps_mod.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__getparam.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

getparam.done: \
          getparam.o \
          ioipsl_getincom.done \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

getparam.o: \
          $(PPSRCDIR0__dyn)/getparam.f90 \
          FFLAGS__dyn__getparam.flags \
          ioipsl_getincom.o \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__divgrad2_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

divgrad2_p.done: \
          divgrad2_p.o \
          mod_hallo.done \
          parallel_lmdz.done \
          times.done
	touch $(FCM_DONEDIR)/$@

divgrad2_p.o: \
          $(PPSRCDIR0__dyn)/divgrad2_p.f \
          FFLAGS__dyn__divgrad2_p.flags \
          mod_hallo.o \
          parallel_lmdz.o \
          times.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__integrd_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

integrd_p.done: \
          integrd_p.o \
          comconst_mod.done \
          comvert_mod.done \
          control_mod.done \
          logic_mod.done \
          parallel_lmdz.done \
          temps_mod.done
	touch $(FCM_DONEDIR)/$@

integrd_p.o: \
          $(PPSRCDIR0__dyn)/integrd_p.f \
          FFLAGS__dyn__integrd_p.flags \
          comconst_mod.o \
          comvert_mod.o \
          control_mod.o \
          logic_mod.o \
          parallel_lmdz.o \
          temps_mod.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__convmas2_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

convmas2_p.done: \
          convmas2_p.o \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

convmas2_p.o: \
          $(PPSRCDIR0__dyn)/convmas2_p.f \
          FFLAGS__dyn__convmas2_p.flags \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__parallel_lmdz.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

parallel_lmdz.done: \
          parallel_lmdz.o \
          vampir.done \
          dimensions.h.idone \
          ioipsl_getincom.done \
          mod_const_mpi.done \
          paramet.h.idone \
          vampir.done
	touch $(FCM_DONEDIR)/$@

parallel_lmdz.o: \
          $(PPSRCDIR0__dyn)/parallel_lmdz.f90 \
          FFLAGS__dyn__parallel_lmdz.flags \
          vampir.o \
          dimensions.h \
          ioipsl_getincom.o \
          mod_const_mpi.o \
          paramet.h \
          vampir.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__filtreg_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

filtreg_p.done: \
          filtreg_p.o \
          filtreg_mod.done \
          mod_filtre_fft.done \
          parallel_lmdz.done \
          timer_filtre.done
	touch $(FCM_DONEDIR)/$@

filtreg_p.o: \
          $(PPSRCDIR0__dyn)/filtreg_p.f \
          FFLAGS__dyn__filtreg_p.flags \
          filtreg_mod.o \
          mod_filtre_fft.o \
          parallel_lmdz.o \
          timer_filtre.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__nxgrad_gam_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

nxgrad_gam_p.done: \
          nxgrad_gam_p.o \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

nxgrad_gam_p.o: \
          $(PPSRCDIR0__dyn)/nxgrad_gam_p.f \
          FFLAGS__dyn__nxgrad_gam_p.flags \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__laplacien_gam_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

laplacien_gam_p.done: \
          laplacien_gam_p.o \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

laplacien_gam_p.o: \
          $(PPSRCDIR0__dyn)/laplacien_gam_p.f \
          FFLAGS__dyn__laplacien_gam_p.flags \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__nxgraro2_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

nxgraro2_p.done: \
          nxgraro2_p.o \
          mod_hallo.done \
          parallel_lmdz.done \
          times.done \
          write_field_p.done
	touch $(FCM_DONEDIR)/$@

nxgraro2_p.o: \
          $(PPSRCDIR0__dyn)/nxgraro2_p.f \
          FFLAGS__dyn__nxgraro2_p.flags \
          mod_hallo.o \
          parallel_lmdz.o \
          times.o \
          write_field_p.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__moyzon_mod.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

moyzon_mod.done: \
          moyzon_mod.o
	touch $(FCM_DONEDIR)/$@

moyzon_mod.o: \
          $(PPSRCDIR0__dyn)/moyzon_mod.f90 \
          FFLAGS__dyn__moyzon_mod.flags
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__caldyn_p.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

caldyn_p.done: \
          caldyn_p.o \
          write_field_p.done \
          comvert_mod.done \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

caldyn_p.o: \
          $(PPSRCDIR0__dyn)/caldyn_p.f \
          FFLAGS__dyn__caldyn_p.flags \
          write_field_p.o \
          comvert_mod.o \
          parallel_lmdz.o
	fcm_internal compile:F dyn $< $@

