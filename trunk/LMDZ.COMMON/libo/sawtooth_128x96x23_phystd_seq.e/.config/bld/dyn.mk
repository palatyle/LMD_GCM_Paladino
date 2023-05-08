# Automatic Make rule for dyn

SRCDIR0__dyn = /home/palatyle/LMD_gen/trunk/LMDZ.COMMON/libf/dyn3d

PPSRCDIR0__dyn = /home/palatyle/LMD_gen/trunk/LMDZ.COMMON/libo/sawtooth_128x96x23_phystd_seq.e/.config/ppsrc/dyn

FFLAGS__dyn__planetary_operations.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

planetary_operations.done: \
          planetary_operations.o \
          comconst_mod.done \
          comgeom.h.idone \
          dimensions.h.idone \
          paramet.h.idone
	touch $(FCM_DONEDIR)/$@

planetary_operations.o: \
          $(PPSRCDIR0__dyn)/planetary_operations.f90 \
          FFLAGS__dyn__planetary_operations.flags \
          comconst_mod.o \
          comgeom.h \
          dimensions.h \
          paramet.h
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__dudv1.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

dudv1.done: \
          dudv1.o
	touch $(FCM_DONEDIR)/$@

dudv1.o: \
          $(PPSRCDIR0__dyn)/dudv1.f \
          FFLAGS__dyn__dudv1.flags
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

FFLAGS__dyn__caldyn.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

caldyn.done: \
          caldyn.o \
          comvert_mod.done
	touch $(FCM_DONEDIR)/$@

caldyn.o: \
          $(PPSRCDIR0__dyn)/caldyn.f \
          FFLAGS__dyn__caldyn.flags \
          comvert_mod.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__vlsplt.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

vlsplt.done: \
          vlsplt.o \
          comconst_mod.done \
          dimensions.h.idone \
          infotrac.done \
          iniprint.h.idone \
          paramet.h.idone
	touch $(FCM_DONEDIR)/$@

vlsplt.o: \
          $(PPSRCDIR0__dyn)/vlsplt.f \
          FFLAGS__dyn__vlsplt.flags \
          comconst_mod.o \
          dimensions.h \
          infotrac.o \
          iniprint.h \
          paramet.h
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__nogcm.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

LDFLAGS__dyn__nogcm.flags: \
          LDFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

LD__dyn__nogcm.flags: \
          LD__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

nogcm.o: \
          $(PPSRCDIR0__dyn)/nogcm.f90 \
          FFLAGS__dyn__nogcm.flags \
          comconst_mod.o \
          comdissnew.h \
          comgeom.h \
          control_mod.o \
          cpdet_mod.o \
          dimensions.h \
          filtreg_mod.o \
          infotrac.o \
          iniphysiq_mod.o \
          iniprint.h \
          ioipsl_getincom.o \
          logic_mod.o \
          mod_const_mpi.o \
          paramet.h \
          temps_mod.o \
          tracstoke.h
	fcm_internal compile:F dyn $< $@

nogcm_128x96x23_phystd_seq.e: \
          nogcm.o \
          LD__dyn__nogcm.flags \
          LDFLAGS__dyn__nogcm.flags \
          $(OBJECTS) \
          comconst_mod.done \
          comdissnew.h.idone \
          comgeom.h.idone \
          control_mod.done \
          cpdet_mod.done \
          dimensions.h.idone \
          filtreg_mod.done \
          infotrac.done \
          iniphysiq_mod.done \
          iniprint.h.idone \
          ioipsl_getincom.done \
          logic_mod.done \
          mod_const_mpi.done \
          paramet.h.idone \
          temps_mod.done \
          tracstoke.h.idone
	fcm_internal load dyn $< $@

FFLAGS__dyn__advect.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

advect.done: \
          advect.o \
          comconst_mod.done \
          comgeom.h.idone \
          dimensions.h.idone \
          ener_mod.done \
          logic_mod.done \
          paramet.h.idone
	touch $(FCM_DONEDIR)/$@

advect.o: \
          $(PPSRCDIR0__dyn)/advect.f \
          FFLAGS__dyn__advect.flags \
          comconst_mod.o \
          comgeom.h \
          dimensions.h \
          ener_mod.o \
          logic_mod.o \
          paramet.h
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

FFLAGS__dyn__top_bound.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

top_bound.done: \
          top_bound.o \
          comconst_mod.done \
          comvert_mod.done
	touch $(FCM_DONEDIR)/$@

top_bound.o: \
          $(PPSRCDIR0__dyn)/top_bound.f \
          FFLAGS__dyn__top_bound.flags \
          comconst_mod.o \
          comvert_mod.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__tetaleveli1j.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

tetaleveli1j.done: \
          tetaleveli1j.o
	touch $(FCM_DONEDIR)/$@

tetaleveli1j.o: \
          $(PPSRCDIR0__dyn)/tetaleveli1j.f \
          FFLAGS__dyn__tetaleveli1j.flags
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__guide_mod.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

guide_mod.done: \
          guide_mod.o \
          write_field.done \
          comconst_mod.done \
          comgeom.h.idone \
          comgeom2.h.idone \
          comvert_mod.done \
          control_mod.done \
          dimensions.h.idone \
          exner_hyb_m.done \
          exner_milieu_m.done \
          getparam.done \
          netcdf95.done \
          paramet.h.idone \
          pres2lev_mod.done \
          serre_mod.done
	touch $(FCM_DONEDIR)/$@

guide_mod.o: \
          $(PPSRCDIR0__dyn)/guide_mod.f90 \
          FFLAGS__dyn__guide_mod.flags \
          write_field.o \
          comconst_mod.o \
          comgeom.h \
          comgeom2.h \
          comvert_mod.o \
          control_mod.o \
          dimensions.h \
          exner_hyb_m.o \
          exner_milieu_m.o \
          getparam.o \
          netcdf95.o \
          paramet.h \
          pres2lev_mod.o \
          serre_mod.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__dudv2.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

dudv2.done: \
          dudv2.o
	touch $(FCM_DONEDIR)/$@

dudv2.o: \
          $(PPSRCDIR0__dyn)/dudv2.f \
          FFLAGS__dyn__dudv2.flags
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__check_isotopes.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

check_isotopes_seq.done: \
          check_isotopes_seq.o \
          infotrac.done
	touch $(FCM_DONEDIR)/$@

check_isotopes_seq.o: \
          $(PPSRCDIR0__dyn)/check_isotopes.f \
          FFLAGS__dyn__check_isotopes.flags \
          infotrac.o
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

write_paramLMDZ_dyn.h: \
          $(SRCDIR0__dyn)/write_paramLMDZ_dyn.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

write_paramLMDZ_dyn.h.idone: \
          $(SRCDIR0__dyn)/write_paramLMDZ_dyn.h
	touch $(FCM_DONEDIR)/$@

FFLAGS__dyn__leapfrog_nogcm.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

leapfrog_nogcm.done: \
          leapfrog_nogcm.o \
          comconst_mod.done \
          comuforc_h.done \
          comvert_mod.done \
          control_mod.done \
          cpdet_mod.done \
          exner_hyb_m.done \
          exner_milieu_m.done \
          guide_mod.done \
          infotrac.done \
          logic_mod.done \
          sponge_mod.done \
          temps_mod.done \
          write_field.done
	touch $(FCM_DONEDIR)/$@

leapfrog_nogcm.o: \
          $(PPSRCDIR0__dyn)/leapfrog_nogcm.f \
          FFLAGS__dyn__leapfrog_nogcm.flags \
          comconst_mod.o \
          comuforc_h.o \
          comvert_mod.o \
          control_mod.o \
          cpdet_mod.o \
          exner_hyb_m.o \
          exner_milieu_m.o \
          guide_mod.o \
          infotrac.o \
          logic_mod.o \
          sponge_mod.o \
          temps_mod.o \
          write_field.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__groupe.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

groupe.done: \
          groupe.o \
          comconst_mod.done
	touch $(FCM_DONEDIR)/$@

groupe.o: \
          $(PPSRCDIR0__dyn)/groupe.f \
          FFLAGS__dyn__groupe.flags \
          comconst_mod.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__fluxstokenc.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

fluxstokenc.done: \
          fluxstokenc.o
	touch $(FCM_DONEDIR)/$@

fluxstokenc.o: \
          $(PPSRCDIR0__dyn)/fluxstokenc.f \
          FFLAGS__dyn__fluxstokenc.flags
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__dissip.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

dissip.done: \
          dissip.o \
          comconst_mod.done
	touch $(FCM_DONEDIR)/$@

dissip.o: \
          $(PPSRCDIR0__dyn)/dissip.f \
          FFLAGS__dyn__dissip.flags \
          comconst_mod.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__bilan_dyn.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

bilan_dyn.done: \
          bilan_dyn.o \
          comconst_mod.done \
          comvert_mod.done \
          control_mod.done \
          cpdet_mod.done \
          temps_mod.done
	touch $(FCM_DONEDIR)/$@

bilan_dyn.o: \
          $(PPSRCDIR0__dyn)/bilan_dyn.f \
          FFLAGS__dyn__bilan_dyn.flags \
          comconst_mod.o \
          comvert_mod.o \
          control_mod.o \
          cpdet_mod.o \
          temps_mod.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__integrd.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

integrd.done: \
          integrd.o \
          comconst_mod.done \
          comvert_mod.done \
          control_mod.done \
          logic_mod.done \
          temps_mod.done
	touch $(FCM_DONEDIR)/$@

integrd.o: \
          $(PPSRCDIR0__dyn)/integrd.f \
          FFLAGS__dyn__integrd.flags \
          comconst_mod.o \
          comvert_mod.o \
          control_mod.o \
          logic_mod.o \
          temps_mod.o
	fcm_internal compile:F dyn $< $@

ini_paramLMDZ_dyn.h: \
          $(SRCDIR0__dyn)/ini_paramLMDZ_dyn.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

ini_paramLMDZ_dyn.h.idone: \
          $(SRCDIR0__dyn)/ini_paramLMDZ_dyn.h
	touch $(FCM_DONEDIR)/$@

FFLAGS__dyn__tetaleveli1j1.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

tetaleveli1j1.done: \
          tetaleveli1j1.o
	touch $(FCM_DONEDIR)/$@

tetaleveli1j1.o: \
          $(PPSRCDIR0__dyn)/tetaleveli1j1.f \
          FFLAGS__dyn__tetaleveli1j1.flags
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
          paramet.h.idone \
          serre_mod.done \
          sponge_mod.done \
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
          paramet.h \
          serre_mod.o \
          sponge_mod.o \
          temps_mod.o
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
          $(PPSRCDIR0__dyn)/gcm.f90 \
          FFLAGS__dyn__gcm.flags \
          comconst_mod.o \
          comdissnew.h \
          comgeom.h \
          control_mod.o \
          cpdet_mod.o \
          dimensions.h \
          filtreg_mod.o \
          infotrac.o \
          iniphysiq_mod.o \
          iniprint.h \
          ioipsl_getincom.o \
          logic_mod.o \
          mod_const_mpi.o \
          paramet.h \
          temps_mod.o \
          tracstoke.h
	fcm_internal compile:F dyn $< $@

gcm_128x96x23_phystd_seq.e: \
          gcm.o \
          LD__dyn__gcm.flags \
          LDFLAGS__dyn__gcm.flags \
          $(OBJECTS) \
          comconst_mod.done \
          comdissnew.h.idone \
          comgeom.h.idone \
          control_mod.done \
          cpdet_mod.done \
          dimensions.h.idone \
          filtreg_mod.done \
          infotrac.done \
          iniphysiq_mod.done \
          iniprint.h.idone \
          ioipsl_getincom.done \
          logic_mod.done \
          mod_const_mpi.done \
          paramet.h.idone \
          temps_mod.done \
          tracstoke.h.idone
	fcm_internal load dyn $< $@

FFLAGS__dyn__mod_const_mpi.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

mod_const_mpi.done: \
          mod_const_mpi.o
	touch $(FCM_DONEDIR)/$@

mod_const_mpi.o: \
          $(PPSRCDIR0__dyn)/mod_const_mpi.f90 \
          FFLAGS__dyn__mod_const_mpi.flags
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__friction.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

friction.done: \
          friction.o \
          comconst_mod.done \
          control_mod.done \
          ioipsl_getincom.done
	touch $(FCM_DONEDIR)/$@

friction.o: \
          $(PPSRCDIR0__dyn)/friction.f \
          FFLAGS__dyn__friction.flags \
          comconst_mod.o \
          control_mod.o \
          ioipsl_getincom.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__qminimum.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

qminimum.done: \
          qminimum.o \
          infotrac.done
	touch $(FCM_DONEDIR)/$@

qminimum.o: \
          $(PPSRCDIR0__dyn)/qminimum.f \
          FFLAGS__dyn__qminimum.flags \
          infotrac.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__dteta1.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

dteta1.done: \
          dteta1.o
	touch $(FCM_DONEDIR)/$@

dteta1.o: \
          $(PPSRCDIR0__dyn)/dteta1.f \
          FFLAGS__dyn__dteta1.flags
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__advtrac.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

advtrac.done: \
          advtrac.o \
          comconst_mod.done \
          comdissip.h.idone \
          comgeom2.h.idone \
          control_mod.done \
          dimensions.h.idone \
          infotrac.done \
          iniprint.h.idone \
          paramet.h.idone
	touch $(FCM_DONEDIR)/$@

advtrac.o: \
          $(PPSRCDIR0__dyn)/advtrac.f90 \
          FFLAGS__dyn__advtrac.flags \
          comconst_mod.o \
          comdissip.h \
          comgeom2.h \
          control_mod.o \
          dimensions.h \
          infotrac.o \
          iniprint.h \
          paramet.h
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__leapfrog.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

leapfrog.done: \
          leapfrog.o \
          comconst_mod.done \
          comuforc_h.done \
          comvert_mod.done \
          control_mod.done \
          cpdet_mod.done \
          exner_hyb_m.done \
          exner_milieu_m.done \
          guide_mod.done \
          infotrac.done \
          logic_mod.done \
          sponge_mod.done \
          temps_mod.done \
          write_field.done
	touch $(FCM_DONEDIR)/$@

leapfrog.o: \
          $(PPSRCDIR0__dyn)/leapfrog.f \
          FFLAGS__dyn__leapfrog.flags \
          comconst_mod.o \
          comuforc_h.o \
          comvert_mod.o \
          control_mod.o \
          cpdet_mod.o \
          exner_hyb_m.o \
          exner_milieu_m.o \
          guide_mod.o \
          infotrac.o \
          logic_mod.o \
          sponge_mod.o \
          temps_mod.o \
          write_field.o
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

FFLAGS__dyn__vlspltqs.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

vlspltqs.done: \
          vlspltqs.o \
          comconst_mod.done \
          cpdet_mod.done \
          infotrac.done
	touch $(FCM_DONEDIR)/$@

vlspltqs.o: \
          $(PPSRCDIR0__dyn)/vlspltqs.f \
          FFLAGS__dyn__vlspltqs.flags \
          comconst_mod.o \
          cpdet_mod.o \
          infotrac.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__abort_gcm.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

abort_gcm.done: \
          abort_gcm.o \
          ioipsl_getincom.done
	touch $(FCM_DONEDIR)/$@

abort_gcm.o: \
          $(PPSRCDIR0__dyn)/abort_gcm.f \
          FFLAGS__dyn__abort_gcm.flags \
          ioipsl_getincom.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__sponge_mod.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

sponge_mod.done: \
          sponge_mod.o \
          comvert_mod.done
	touch $(FCM_DONEDIR)/$@

sponge_mod.o: \
          $(PPSRCDIR0__dyn)/sponge_mod.f90 \
          FFLAGS__dyn__sponge_mod.flags \
          comvert_mod.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__getparam.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

getparam.done: \
          getparam.o \
          ioipsl_getincom.done
	touch $(FCM_DONEDIR)/$@

getparam.o: \
          $(PPSRCDIR0__dyn)/getparam.f90 \
          FFLAGS__dyn__getparam.flags \
          ioipsl_getincom.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__caladvtrac.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

caladvtrac.done: \
          caladvtrac.o \
          comconst_mod.done \
          control_mod.done \
          infotrac.done \
          planetary_operations.done
	touch $(FCM_DONEDIR)/$@

caladvtrac.o: \
          $(PPSRCDIR0__dyn)/caladvtrac.f \
          FFLAGS__dyn__caladvtrac.flags \
          comconst_mod.o \
          control_mod.o \
          infotrac.o \
          planetary_operations.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__groupeun.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

groupeun.done: \
          groupeun.o \
          comconst_mod.done
	touch $(FCM_DONEDIR)/$@

groupeun.o: \
          $(PPSRCDIR0__dyn)/groupeun.f \
          FFLAGS__dyn__groupeun.flags \
          comconst_mod.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__iniinterp_horiz.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

iniinterp_horiz.done: \
          iniinterp_horiz.o
	touch $(FCM_DONEDIR)/$@

iniinterp_horiz.o: \
          $(PPSRCDIR0__dyn)/iniinterp_horiz.f \
          FFLAGS__dyn__iniinterp_horiz.flags
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

FFLAGS__dyn__addfi.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

addfi.done: \
          addfi.o \
          comconst_mod.done \
          control_mod.done \
          infotrac.done
	touch $(FCM_DONEDIR)/$@

addfi.o: \
          $(PPSRCDIR0__dyn)/addfi.f \
          FFLAGS__dyn__addfi.flags \
          comconst_mod.o \
          control_mod.o \
          infotrac.o
	fcm_internal compile:F dyn $< $@

FFLAGS__dyn__interp_horiz.flags: \
          FFLAGS__dyn.flags
	touch $(FCM_FLAGSDIR)/$@

interp_horiz.done: \
          interp_horiz.o
	touch $(FCM_DONEDIR)/$@

interp_horiz.o: \
          $(PPSRCDIR0__dyn)/interp_horiz.f \
          FFLAGS__dyn__interp_horiz.flags
	fcm_internal compile:F dyn $< $@

