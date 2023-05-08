# Automatic Make rule for dyn3d_common

SRCDIR0__dyn3d_common = /home/palatyle/LMD_gen/trunk/LMDZ.COMMON/libf/dyn3d_common

PPSRCDIR0__dyn3d_common = /home/palatyle/LMD_gen/trunk/LMDZ.COMMON/libo/sawtooth_128x96x23_phystd_seq.e/.config/ppsrc/dyn3d_common

FFLAGS__dyn3d_common__pentes_ini.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

pentes_ini.done: \
          pentes_ini.o \
          comconst_mod.done
	touch $(FCM_DONEDIR)/$@

pentes_ini.o: \
          $(PPSRCDIR0__dyn3d_common)/pentes_ini.f \
          FFLAGS__dyn3d_common__pentes_ini.flags \
          comconst_mod.o
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__vitvert.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

vitvert.done: \
          vitvert.o \
          comvert_mod.done \
          dimensions.h.idone \
          paramet.h.idone
	touch $(FCM_DONEDIR)/$@

vitvert.o: \
          $(PPSRCDIR0__dyn3d_common)/vitvert.f90 \
          FFLAGS__dyn3d_common__vitvert.flags \
          comvert_mod.o \
          dimensions.h \
          paramet.h
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__misc_mod.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

misc_mod.done: \
          misc_mod.o
	touch $(FCM_DONEDIR)/$@

misc_mod.o: \
          $(PPSRCDIR0__dyn3d_common)/misc_mod.f90 \
          FFLAGS__dyn3d_common__misc_mod.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__dynredem_mod.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

dynredem_mod.done: \
          dynredem_mod.o \
          dimensions.h.idone \
          paramet.h.idone
	touch $(FCM_DONEDIR)/$@

dynredem_mod.o: \
          $(PPSRCDIR0__dyn3d_common)/dynredem_mod.f90 \
          FFLAGS__dyn3d_common__dynredem_mod.flags \
          dimensions.h \
          paramet.h
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__nxgraro2.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

nxgraro2.done: \
          nxgraro2.o
	touch $(FCM_DONEDIR)/$@

nxgraro2.o: \
          $(PPSRCDIR0__dyn3d_common)/nxgraro2.f \
          FFLAGS__dyn3d_common__nxgraro2.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__advyp.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

advyp.done: \
          advyp.o
	touch $(FCM_DONEDIR)/$@

advyp.o: \
          $(PPSRCDIR0__dyn3d_common)/advyp.f \
          FFLAGS__dyn3d_common__advyp.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__interpost.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

interpost.done: \
          interpost.o
	touch $(FCM_DONEDIR)/$@

interpost.o: \
          $(PPSRCDIR0__dyn3d_common)/interpost.f \
          FFLAGS__dyn3d_common__interpost.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__conf_planete.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

conf_planete.done: \
          conf_planete.o \
          comconst_mod.done \
          comvert_mod.done \
          ioipsl_getincom.done
	touch $(FCM_DONEDIR)/$@

conf_planete.o: \
          $(PPSRCDIR0__dyn3d_common)/conf_planete.f90 \
          FFLAGS__dyn3d_common__conf_planete.flags \
          comconst_mod.o \
          comvert_mod.o \
          ioipsl_getincom.o
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__nxgradst.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

nxgradst.done: \
          nxgradst.o
	touch $(FCM_DONEDIR)/$@

nxgradst.o: \
          $(PPSRCDIR0__dyn3d_common)/nxgradst.f \
          FFLAGS__dyn3d_common__nxgradst.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__gr_v_scal.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

gr_v_scal.done: \
          gr_v_scal.o
	touch $(FCM_DONEDIR)/$@

gr_v_scal.o: \
          $(PPSRCDIR0__dyn3d_common)/gr_v_scal.f \
          FFLAGS__dyn3d_common__gr_v_scal.flags
	fcm_internal compile:F dyn3d_common $< $@

write_grads_dyn.h: \
          $(SRCDIR0__dyn3d_common)/write_grads_dyn.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

write_grads_dyn.h.idone: \
          $(SRCDIR0__dyn3d_common)/write_grads_dyn.h
	touch $(FCM_DONEDIR)/$@

FFLAGS__dyn3d_common__nxgrarot.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

nxgrarot.done: \
          nxgrarot.o
	touch $(FCM_DONEDIR)/$@

nxgrarot.o: \
          $(PPSRCDIR0__dyn3d_common)/nxgrarot.f \
          FFLAGS__dyn3d_common__nxgrarot.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__advzp.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

advzp.done: \
          advzp.o
	touch $(FCM_DONEDIR)/$@

advzp.o: \
          $(PPSRCDIR0__dyn3d_common)/advzp.f \
          FFLAGS__dyn3d_common__advzp.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__psextbar.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

psextbar.done: \
          psextbar.o
	touch $(FCM_DONEDIR)/$@

psextbar.o: \
          $(PPSRCDIR0__dyn3d_common)/psextbar.f \
          FFLAGS__dyn3d_common__psextbar.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__exner_milieu_m.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

exner_milieu_m.done: \
          exner_milieu_m.o \
          comconst_mod.done \
          comgeom.h.idone \
          comvert_mod.done \
          dimensions.h.idone \
          paramet.h.idone
	touch $(FCM_DONEDIR)/$@

exner_milieu_m.o: \
          $(PPSRCDIR0__dyn3d_common)/exner_milieu_m.f90 \
          FFLAGS__dyn3d_common__exner_milieu_m.flags \
          comconst_mod.o \
          comgeom.h \
          comvert_mod.o \
          dimensions.h \
          paramet.h
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__fxysinus.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

fxysinus.done: \
          fxysinus.o \
          comconst_mod.done
	touch $(FCM_DONEDIR)/$@

fxysinus.o: \
          $(PPSRCDIR0__dyn3d_common)/fxysinus.f \
          FFLAGS__dyn3d_common__fxysinus.flags \
          comconst_mod.o
	fcm_internal compile:F dyn3d_common $< $@

comdissipn.h: \
          $(SRCDIR0__dyn3d_common)/comdissipn.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

comdissipn.h.idone: \
          $(SRCDIR0__dyn3d_common)/comdissipn.h
	touch $(FCM_DONEDIR)/$@

FFLAGS__dyn3d_common__geopot.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

geopot.done: \
          geopot.o
	touch $(FCM_DONEDIR)/$@

geopot.o: \
          $(PPSRCDIR0__dyn3d_common)/geopot.f \
          FFLAGS__dyn3d_common__geopot.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__inigrads.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

inigrads.done: \
          inigrads.o
	touch $(FCM_DONEDIR)/$@

inigrads.o: \
          $(PPSRCDIR0__dyn3d_common)/inigrads.f \
          FFLAGS__dyn3d_common__inigrads.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__disvert_noterre.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

disvert_noterre.done: \
          disvert_noterre.o \
          comconst_mod.done \
          comvert_mod.done \
          ioipsl_getincom.done \
          logic_mod.done
	touch $(FCM_DONEDIR)/$@

disvert_noterre.o: \
          $(PPSRCDIR0__dyn3d_common)/disvert_noterre.f \
          FFLAGS__dyn3d_common__disvert_noterre.flags \
          comconst_mod.o \
          comvert_mod.o \
          ioipsl_getincom.o \
          logic_mod.o
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__laplacien_rotgam.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

laplacien_rotgam.done: \
          laplacien_rotgam.o
	touch $(FCM_DONEDIR)/$@

laplacien_rotgam.o: \
          $(PPSRCDIR0__dyn3d_common)/laplacien_rotgam.f \
          FFLAGS__dyn3d_common__laplacien_rotgam.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__divergst.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

divergst.done: \
          divergst.o
	touch $(FCM_DONEDIR)/$@

divergst.o: \
          $(PPSRCDIR0__dyn3d_common)/divergst.f \
          FFLAGS__dyn3d_common__divergst.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__inigeom.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

inigeom.done: \
          inigeom.o \
          comconst_mod.done \
          fxhyp_m.done \
          fyhyp_m.done \
          logic_mod.done \
          serre_mod.done
	touch $(FCM_DONEDIR)/$@

inigeom.o: \
          $(PPSRCDIR0__dyn3d_common)/inigeom.f \
          FFLAGS__dyn3d_common__inigeom.flags \
          comconst_mod.o \
          fxhyp_m.o \
          fyhyp_m.o \
          logic_mod.o \
          serre_mod.o
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__diverg.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

diverg.done: \
          diverg.o
	touch $(FCM_DONEDIR)/$@

diverg.o: \
          $(PPSRCDIR0__dyn3d_common)/diverg.f \
          FFLAGS__dyn3d_common__diverg.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__limy.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

limy.done: \
          limy.o \
          comconst_mod.done
	touch $(FCM_DONEDIR)/$@

limy.o: \
          $(PPSRCDIR0__dyn3d_common)/limy.f \
          FFLAGS__dyn3d_common__limy.flags \
          comconst_mod.o
	fcm_internal compile:F dyn3d_common $< $@

academic.h: \
          $(SRCDIR0__dyn3d_common)/academic.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

academic.h.idone: \
          $(SRCDIR0__dyn3d_common)/academic.h
	touch $(FCM_DONEDIR)/$@

FFLAGS__dyn3d_common__covnat.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

covnat.done: \
          covnat.o
	touch $(FCM_DONEDIR)/$@

covnat.o: \
          $(PPSRCDIR0__dyn3d_common)/covnat.f \
          FFLAGS__dyn3d_common__covnat.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__principal_cshift_m.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

principal_cshift_m.done: \
          principal_cshift_m.o \
          dimensions.h.idone \
          nrtype.done \
          serre_mod.done
	touch $(FCM_DONEDIR)/$@

principal_cshift_m.o: \
          $(PPSRCDIR0__dyn3d_common)/principal_cshift_m.f90 \
          FFLAGS__dyn3d_common__principal_cshift_m.flags \
          dimensions.h \
          nrtype.o \
          serre_mod.o
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__conf_dat_m.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

conf_dat_m.done: \
          conf_dat_m.o \
          assert_eq_m.done
	touch $(FCM_DONEDIR)/$@

conf_dat_m.o: \
          $(PPSRCDIR0__dyn3d_common)/conf_dat_m.f90 \
          FFLAGS__dyn3d_common__conf_dat_m.flags \
          assert_eq_m.o
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__relaxetat0.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

relaxetat0.done: \
          relaxetat0.o \
          comuforc_h.done \
          ioipsl_getincom.done
	touch $(FCM_DONEDIR)/$@

relaxetat0.o: \
          $(PPSRCDIR0__dyn3d_common)/relaxetat0.f \
          FFLAGS__dyn3d_common__relaxetat0.flags \
          comuforc_h.o \
          ioipsl_getincom.o
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__advn.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

advn.done: \
          advn.o
	touch $(FCM_DONEDIR)/$@

advn.o: \
          $(PPSRCDIR0__dyn3d_common)/advn.f \
          FFLAGS__dyn3d_common__advn.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__tourpot.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

tourpot.done: \
          tourpot.o \
          comgeom.h.idone \
          dimensions.h.idone \
          paramet.h.idone
	touch $(FCM_DONEDIR)/$@

tourpot.o: \
          $(PPSRCDIR0__dyn3d_common)/tourpot.f90 \
          FFLAGS__dyn3d_common__tourpot.flags \
          comgeom.h \
          dimensions.h \
          paramet.h
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__divgrad2.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

divgrad2.done: \
          divgrad2.o
	touch $(FCM_DONEDIR)/$@

divgrad2.o: \
          $(PPSRCDIR0__dyn3d_common)/divgrad2.f \
          FFLAGS__dyn3d_common__divgrad2.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__initdynav.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

initdynav.done: \
          initdynav.o \
          com_io_dyn_mod.done \
          comconst_mod.done \
          comgeom.h.idone \
          comvert_mod.done \
          dimensions.h.idone \
          infotrac.done \
          iniprint.h.idone \
          paramet.h.idone \
          temps_mod.done
	touch $(FCM_DONEDIR)/$@

initdynav.o: \
          $(PPSRCDIR0__dyn3d_common)/initdynav.f90 \
          FFLAGS__dyn3d_common__initdynav.flags \
          com_io_dyn_mod.o \
          comconst_mod.o \
          comgeom.h \
          comvert_mod.o \
          dimensions.h \
          infotrac.o \
          iniprint.h \
          paramet.h \
          temps_mod.o
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__massbarxy.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

massbarxy.done: \
          massbarxy.o \
          comgeom.h.idone \
          dimensions.h.idone \
          paramet.h.idone
	touch $(FCM_DONEDIR)/$@

massbarxy.o: \
          $(PPSRCDIR0__dyn3d_common)/massbarxy.f90 \
          FFLAGS__dyn3d_common__massbarxy.flags \
          comgeom.h \
          dimensions.h \
          paramet.h
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__diverg_gam.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

diverg_gam.done: \
          diverg_gam.o
	touch $(FCM_DONEDIR)/$@

diverg_gam.o: \
          $(PPSRCDIR0__dyn3d_common)/diverg_gam.f \
          FFLAGS__dyn3d_common__diverg_gam.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__disvert.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

disvert.done: \
          disvert.o \
          assert_m.done \
          comvert_mod.done \
          dimensions.h.idone \
          iniprint.h.idone \
          ioipsl_getincom.done \
          logic_mod.done \
          new_unit_m.done \
          paramet.h.idone
	touch $(FCM_DONEDIR)/$@

disvert.o: \
          $(PPSRCDIR0__dyn3d_common)/disvert.f90 \
          FFLAGS__dyn3d_common__disvert.flags \
          assert_m.o \
          comvert_mod.o \
          dimensions.h \
          iniprint.h \
          ioipsl_getincom.o \
          logic_mod.o \
          new_unit_m.o \
          paramet.h
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__test_period.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

test_period.done: \
          test_period.o \
          infotrac.done
	touch $(FCM_DONEDIR)/$@

test_period.o: \
          $(PPSRCDIR0__dyn3d_common)/test_period.f \
          FFLAGS__dyn3d_common__test_period.flags \
          infotrac.o
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__coordij.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

coordij.done: \
          coordij.o \
          comconst_mod.done
	touch $(FCM_DONEDIR)/$@

coordij.o: \
          $(PPSRCDIR0__dyn3d_common)/coordij.f \
          FFLAGS__dyn3d_common__coordij.flags \
          comconst_mod.o
	fcm_internal compile:F dyn3d_common $< $@

paramet.h: \
          $(SRCDIR0__dyn3d_common)/paramet.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

paramet.h.idone: \
          $(SRCDIR0__dyn3d_common)/paramet.h
	touch $(FCM_DONEDIR)/$@

comdissip.h: \
          $(SRCDIR0__dyn3d_common)/comdissip.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

comdissip.h.idone: \
          $(SRCDIR0__dyn3d_common)/comdissip.h
	touch $(FCM_DONEDIR)/$@

FFLAGS__dyn3d_common__laplacien.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

laplacien.done: \
          laplacien.o
	touch $(FCM_DONEDIR)/$@

laplacien.o: \
          $(PPSRCDIR0__dyn3d_common)/laplacien.f \
          FFLAGS__dyn3d_common__laplacien.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__adaptdt.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

adaptdt.done: \
          adaptdt.o \
          comconst_mod.done \
          control_mod.done
	touch $(FCM_DONEDIR)/$@

adaptdt.o: \
          $(PPSRCDIR0__dyn3d_common)/adaptdt.f \
          FFLAGS__dyn3d_common__adaptdt.flags \
          comconst_mod.o \
          control_mod.o
	fcm_internal compile:F dyn3d_common $< $@

gradsdef.h: \
          $(SRCDIR0__dyn3d_common)/gradsdef.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

gradsdef.h.idone: \
          $(SRCDIR0__dyn3d_common)/gradsdef.h
	touch $(FCM_DONEDIR)/$@

FFLAGS__dyn3d_common__grad.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

grad.done: \
          grad.o
	touch $(FCM_DONEDIR)/$@

grad.o: \
          $(PPSRCDIR0__dyn3d_common)/grad.f \
          FFLAGS__dyn3d_common__grad.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__infotrac.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

infotrac.done: \
          infotrac.o \
          control_mod.done \
          dimensions.h.idone \
          iniprint.h.idone \
          ioipsl_getincom.done
	touch $(FCM_DONEDIR)/$@

infotrac.o: \
          $(PPSRCDIR0__dyn3d_common)/infotrac.f90 \
          FFLAGS__dyn3d_common__infotrac.flags \
          control_mod.o \
          dimensions.h \
          iniprint.h \
          ioipsl_getincom.o
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__comuforc_h.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

comuforc_h.done: \
          comuforc_h.o
	touch $(FCM_DONEDIR)/$@

comuforc_h.o: \
          $(PPSRCDIR0__dyn3d_common)/comuforc_h.f90 \
          FFLAGS__dyn3d_common__comuforc_h.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__invert_lat.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

invert_lat.done: \
          invert_lat.o
	touch $(FCM_DONEDIR)/$@

invert_lat.o: \
          $(PPSRCDIR0__dyn3d_common)/invert_lat.f90 \
          FFLAGS__dyn3d_common__invert_lat.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__dynredem.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

dynredem0.done: \
          dynredem0.o \
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
          paramet.h.idone \
          serre_mod.done \
          temps_mod.done
	touch $(FCM_DONEDIR)/$@

dynredem0.o: \
          $(PPSRCDIR0__dyn3d_common)/dynredem.f90 \
          FFLAGS__dyn3d_common__dynredem.flags \
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
          paramet.h \
          serre_mod.o \
          temps_mod.o
	fcm_internal compile:F dyn3d_common $< $@

comgeom2.h: \
          $(SRCDIR0__dyn3d_common)/comgeom2.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

comgeom2.h.idone: \
          $(SRCDIR0__dyn3d_common)/comgeom2.h
	touch $(FCM_DONEDIR)/$@

FFLAGS__dyn3d_common__fxyhyper.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

fxyhyper.done: \
          fxyhyper.o
	touch $(FCM_DONEDIR)/$@

fxyhyper.o: \
          $(PPSRCDIR0__dyn3d_common)/fxyhyper.f \
          FFLAGS__dyn3d_common__fxyhyper.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__rotatst.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

rotatst.done: \
          rotatst.o
	touch $(FCM_DONEDIR)/$@

rotatst.o: \
          $(PPSRCDIR0__dyn3d_common)/rotatst.f \
          FFLAGS__dyn3d_common__rotatst.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__heavyside.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

heavyside.done: \
          heavyside.o
	touch $(FCM_DONEDIR)/$@

heavyside.o: \
          $(PPSRCDIR0__dyn3d_common)/heavyside.f \
          FFLAGS__dyn3d_common__heavyside.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__fyhyp_m.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

fyhyp_m.done: \
          fyhyp_m.o \
          coefpoly_m.done \
          dimensions.h.idone \
          nrtype.done \
          serre_mod.done
	touch $(FCM_DONEDIR)/$@

fyhyp_m.o: \
          $(PPSRCDIR0__dyn3d_common)/fyhyp_m.f90 \
          FFLAGS__dyn3d_common__fyhyp_m.flags \
          coefpoly_m.o \
          dimensions.h \
          nrtype.o \
          serre_mod.o
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__laplacien_rot.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

laplacien_rot.done: \
          laplacien_rot.o
	touch $(FCM_DONEDIR)/$@

laplacien_rot.o: \
          $(PPSRCDIR0__dyn3d_common)/laplacien_rot.f \
          FFLAGS__dyn3d_common__laplacien_rot.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__laplacien_gam.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

laplacien_gam.done: \
          laplacien_gam.o
	touch $(FCM_DONEDIR)/$@

laplacien_gam.o: \
          $(PPSRCDIR0__dyn3d_common)/laplacien_gam.f \
          FFLAGS__dyn3d_common__laplacien_gam.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__comvert_mod.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

comvert_mod.done: \
          comvert_mod.o \
          dimensions.h.idone
	touch $(FCM_DONEDIR)/$@

comvert_mod.o: \
          $(PPSRCDIR0__dyn3d_common)/comvert_mod.f90 \
          FFLAGS__dyn3d_common__comvert_mod.flags \
          dimensions.h
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__convflu.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

convflu.done: \
          convflu.o
	touch $(FCM_DONEDIR)/$@

convflu.o: \
          $(PPSRCDIR0__dyn3d_common)/convflu.f \
          FFLAGS__dyn3d_common__convflu.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__limz.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

limz.done: \
          limz.o
	touch $(FCM_DONEDIR)/$@

limz.o: \
          $(PPSRCDIR0__dyn3d_common)/limz.f \
          FFLAGS__dyn3d_common__limz.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__advz.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

advz.done: \
          advz.o
	touch $(FCM_DONEDIR)/$@

advz.o: \
          $(PPSRCDIR0__dyn3d_common)/advz.f \
          FFLAGS__dyn3d_common__advz.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__inter_barxy_m.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

inter_barxy_m.done: \
          inter_barxy_m.o \
          assert_eq_m.done \
          assert_m.done \
          comconst_mod.done \
          comgeom2.h.idone \
          dimensions.h.idone \
          paramet.h.idone
	touch $(FCM_DONEDIR)/$@

inter_barxy_m.o: \
          $(PPSRCDIR0__dyn3d_common)/inter_barxy_m.f90 \
          FFLAGS__dyn3d_common__inter_barxy_m.flags \
          assert_eq_m.o \
          assert_m.o \
          comconst_mod.o \
          comgeom2.h \
          dimensions.h \
          paramet.h
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__convmas.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

convmas.done: \
          convmas.o \
          comgeom.h.idone \
          dimensions.h.idone \
          paramet.h.idone
	touch $(FCM_DONEDIR)/$@

convmas.o: \
          $(PPSRCDIR0__dyn3d_common)/convmas.f90 \
          FFLAGS__dyn3d_common__convmas.flags \
          comgeom.h \
          dimensions.h \
          paramet.h
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__ppm3d.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

ppm3d.done: \
          ppm3d.o
	touch $(FCM_DONEDIR)/$@

ppm3d.o: \
          $(PPSRCDIR0__dyn3d_common)/ppm3d.f \
          FFLAGS__dyn3d_common__ppm3d.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__rotat_nfil.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

rotat_nfil.done: \
          rotat_nfil.o
	touch $(FCM_DONEDIR)/$@

rotat_nfil.o: \
          $(PPSRCDIR0__dyn3d_common)/rotat_nfil.f \
          FFLAGS__dyn3d_common__rotat_nfil.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__interpre.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

interpre.done: \
          interpre.o \
          comconst_mod.done \
          comvert_mod.done \
          control_mod.done
	touch $(FCM_DONEDIR)/$@

interpre.o: \
          $(PPSRCDIR0__dyn3d_common)/interpre.f \
          FFLAGS__dyn3d_common__interpre.flags \
          comconst_mod.o \
          comvert_mod.o \
          control_mod.o
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__fxy.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

fxy.done: \
          fxy.o \
          comconst_mod.done \
          serre_mod.done
	touch $(FCM_DONEDIR)/$@

fxy.o: \
          $(PPSRCDIR0__dyn3d_common)/fxy.f \
          FFLAGS__dyn3d_common__fxy.flags \
          comconst_mod.o \
          serre_mod.o
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__massdair.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

massdair.done: \
          massdair.o
	touch $(FCM_DONEDIR)/$@

massdair.o: \
          $(PPSRCDIR0__dyn3d_common)/massdair.f \
          FFLAGS__dyn3d_common__massdair.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__ugeostr.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

ugeostr.done: \
          ugeostr.o \
          comconst_mod.done \
          comgeom2.h.idone \
          dimensions.h.idone \
          paramet.h.idone
	touch $(FCM_DONEDIR)/$@

ugeostr.o: \
          $(PPSRCDIR0__dyn3d_common)/ugeostr.f90 \
          FFLAGS__dyn3d_common__ugeostr.flags \
          comconst_mod.o \
          comgeom2.h \
          dimensions.h \
          paramet.h
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__initfluxsto.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

initfluxsto.done: \
          initfluxsto.o \
          comconst_mod.done \
          comvert_mod.done \
          temps_mod.done
	touch $(FCM_DONEDIR)/$@

initfluxsto.o: \
          $(PPSRCDIR0__dyn3d_common)/initfluxsto.f \
          FFLAGS__dyn3d_common__initfluxsto.flags \
          comconst_mod.o \
          comvert_mod.o \
          temps_mod.o
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__extrapol.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

extrapol.done: \
          extrapol.o
	touch $(FCM_DONEDIR)/$@

extrapol.o: \
          $(PPSRCDIR0__dyn3d_common)/extrapol.f \
          FFLAGS__dyn3d_common__extrapol.flags
	fcm_internal compile:F dyn3d_common $< $@

description.h: \
          $(SRCDIR0__dyn3d_common)/description.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

description.h.idone: \
          $(SRCDIR0__dyn3d_common)/description.h
	touch $(FCM_DONEDIR)/$@

FFLAGS__dyn3d_common__ener_mod.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

ener_mod.done: \
          ener_mod.o \
          dimensions.h.idone \
          paramet.h.idone
	touch $(FCM_DONEDIR)/$@

ener_mod.o: \
          $(PPSRCDIR0__dyn3d_common)/ener_mod.f90 \
          FFLAGS__dyn3d_common__ener_mod.flags \
          dimensions.h \
          paramet.h
	fcm_internal compile:F dyn3d_common $< $@

comgeom.h: \
          $(SRCDIR0__dyn3d_common)/comgeom.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

comgeom.h.idone: \
          $(SRCDIR0__dyn3d_common)/comgeom.h
	touch $(FCM_DONEDIR)/$@

FFLAGS__dyn3d_common__enercin.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

enercin.done: \
          enercin.o \
          comgeom.h.idone \
          dimensions.h.idone \
          paramet.h.idone
	touch $(FCM_DONEDIR)/$@

enercin.o: \
          $(PPSRCDIR0__dyn3d_common)/enercin.f90 \
          FFLAGS__dyn3d_common__enercin.flags \
          comgeom.h \
          dimensions.h \
          paramet.h
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__advxp.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

advxp.done: \
          advxp.o
	touch $(FCM_DONEDIR)/$@

advxp.o: \
          $(PPSRCDIR0__dyn3d_common)/advxp.f \
          FFLAGS__dyn3d_common__advxp.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__flumass.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

flumass.done: \
          flumass.o \
          comgeom.h.idone \
          dimensions.h.idone \
          paramet.h.idone
	touch $(FCM_DONEDIR)/$@

flumass.o: \
          $(PPSRCDIR0__dyn3d_common)/flumass.f90 \
          FFLAGS__dyn3d_common__flumass.flags \
          comgeom.h \
          dimensions.h \
          paramet.h
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__gradiv.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

gradiv.done: \
          gradiv.o
	touch $(FCM_DONEDIR)/$@

gradiv.o: \
          $(PPSRCDIR0__dyn3d_common)/gradiv.f \
          FFLAGS__dyn3d_common__gradiv.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__writedynav.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

writedynav.done: \
          writedynav.o \
          com_io_dyn_mod.done \
          comconst_mod.done \
          comgeom.h.idone \
          dimensions.h.idone \
          infotrac.done \
          iniprint.h.idone \
          paramet.h.idone \
          temps_mod.done
	touch $(FCM_DONEDIR)/$@

writedynav.o: \
          $(PPSRCDIR0__dyn3d_common)/writedynav.f90 \
          FFLAGS__dyn3d_common__writedynav.flags \
          com_io_dyn_mod.o \
          comconst_mod.o \
          comgeom.h \
          dimensions.h \
          infotrac.o \
          iniprint.h \
          paramet.h \
          temps_mod.o
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__rotat.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

rotat.done: \
          rotat.o
	touch $(FCM_DONEDIR)/$@

rotat.o: \
          $(PPSRCDIR0__dyn3d_common)/rotat.f \
          FFLAGS__dyn3d_common__rotat.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__sortvarc.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

sortvarc.done: \
          sortvarc.o \
          comconst_mod.done \
          comgeom.h.idone \
          control_mod.done \
          dimensions.h.idone \
          ener_mod.done \
          iniprint.h.idone \
          logic_mod.done \
          paramet.h.idone
	touch $(FCM_DONEDIR)/$@

sortvarc.o: \
          $(PPSRCDIR0__dyn3d_common)/sortvarc.f \
          FFLAGS__dyn3d_common__sortvarc.flags \
          comconst_mod.o \
          comgeom.h \
          control_mod.o \
          dimensions.h \
          ener_mod.o \
          iniprint.h \
          logic_mod.o \
          paramet.h
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__dump2d.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

dump2d.done: \
          dump2d.o
	touch $(FCM_DONEDIR)/$@

dump2d.o: \
          $(PPSRCDIR0__dyn3d_common)/dump2d.f \
          FFLAGS__dyn3d_common__dump2d.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__divergf.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

divergf.done: \
          divergf.o
	touch $(FCM_DONEDIR)/$@

divergf.o: \
          $(PPSRCDIR0__dyn3d_common)/divergf.f \
          FFLAGS__dyn3d_common__divergf.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__dynetat0.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

dynetat0.done: \
          dynetat0.o \
          comconst_mod.done \
          comgeom2.h.idone \
          comvert_mod.done \
          control_mod.done \
          dimensions.h.idone \
          ener_mod.done \
          infotrac.done \
          iniprint.h.idone \
          logic_mod.done \
          paramet.h.idone \
          serre_mod.done \
          temps_mod.done
	touch $(FCM_DONEDIR)/$@

dynetat0.o: \
          $(PPSRCDIR0__dyn3d_common)/dynetat0.f90 \
          FFLAGS__dyn3d_common__dynetat0.flags \
          comconst_mod.o \
          comgeom2.h \
          comvert_mod.o \
          control_mod.o \
          dimensions.h \
          ener_mod.o \
          infotrac.o \
          iniprint.h \
          logic_mod.o \
          paramet.h \
          serre_mod.o \
          temps_mod.o
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__serre_mod.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

serre_mod.done: \
          serre_mod.o
	touch $(FCM_DONEDIR)/$@

serre_mod.o: \
          $(PPSRCDIR0__dyn3d_common)/serre_mod.f90 \
          FFLAGS__dyn3d_common__serre_mod.flags
	fcm_internal compile:F dyn3d_common $< $@

tracstoke.h: \
          $(SRCDIR0__dyn3d_common)/tracstoke.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

tracstoke.h.idone: \
          $(SRCDIR0__dyn3d_common)/tracstoke.h
	touch $(FCM_DONEDIR)/$@

FFLAGS__dyn3d_common__iniconst.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

iniconst.done: \
          iniconst.o \
          comconst_mod.done \
          comvert_mod.done \
          control_mod.done \
          dimensions.h.idone \
          iniprint.h.idone \
          ioipsl_getincom.done \
          paramet.h.idone
	touch $(FCM_DONEDIR)/$@

iniconst.o: \
          $(PPSRCDIR0__dyn3d_common)/iniconst.f90 \
          FFLAGS__dyn3d_common__iniconst.flags \
          comconst_mod.o \
          comvert_mod.o \
          control_mod.o \
          dimensions.h \
          iniprint.h \
          ioipsl_getincom.o \
          paramet.h
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__control_mod.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

control_mod.done: \
          control_mod.o
	touch $(FCM_DONEDIR)/$@

control_mod.o: \
          $(PPSRCDIR0__dyn3d_common)/control_mod.f90 \
          FFLAGS__dyn3d_common__control_mod.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__grilles_gcm_netcdf_sub.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

grilles_gcm_netcdf_sub.done: \
          grilles_gcm_netcdf_sub.o \
          comconst_mod.done \
          comgeom.h.idone \
          comvert_mod.done \
          dimensions.h.idone \
          paramet.h.idone
	touch $(FCM_DONEDIR)/$@

grilles_gcm_netcdf_sub.o: \
          $(PPSRCDIR0__dyn3d_common)/grilles_gcm_netcdf_sub.f90 \
          FFLAGS__dyn3d_common__grilles_gcm_netcdf_sub.flags \
          comconst_mod.o \
          comgeom.h \
          comvert_mod.o \
          dimensions.h \
          paramet.h
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__iniacademic.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

iniacademic.done: \
          iniacademic.o \
          write_field.done \
          academic.h.idone \
          comconst_mod.done \
          comgeom.h.idone \
          comvert_mod.done \
          control_mod.done \
          dimensions.h.idone \
          ener_mod.done \
          exner_hyb_m.done \
          exner_milieu_m.done \
          filtreg_mod.done \
          infotrac.done \
          iniprint.h.idone \
          ioipsl_getincom.done \
          logic_mod.done \
          paramet.h.idone \
          temps_mod.done
	touch $(FCM_DONEDIR)/$@

iniacademic.o: \
          $(PPSRCDIR0__dyn3d_common)/iniacademic.f90 \
          FFLAGS__dyn3d_common__iniacademic.flags \
          write_field.o \
          academic.h \
          comconst_mod.o \
          comgeom.h \
          comvert_mod.o \
          control_mod.o \
          dimensions.h \
          ener_mod.o \
          exner_hyb_m.o \
          exner_milieu_m.o \
          filtreg_mod.o \
          infotrac.o \
          iniprint.h \
          ioipsl_getincom.o \
          logic_mod.o \
          paramet.h \
          temps_mod.o
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__covcont.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

covcont.done: \
          covcont.o \
          comgeom.h.idone \
          dimensions.h.idone \
          paramet.h.idone
	touch $(FCM_DONEDIR)/$@

covcont.o: \
          $(PPSRCDIR0__dyn3d_common)/covcont.f90 \
          FFLAGS__dyn3d_common__covcont.flags \
          comgeom.h \
          dimensions.h \
          paramet.h
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__tidal_forces.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

tidal_forces.done: \
          tidal_forces.o
	touch $(FCM_DONEDIR)/$@

tidal_forces.o: \
          $(PPSRCDIR0__dyn3d_common)/tidal_forces.f \
          FFLAGS__dyn3d_common__tidal_forces.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__pbar.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

pbar.done: \
          pbar.o
	touch $(FCM_DONEDIR)/$@

pbar.o: \
          $(PPSRCDIR0__dyn3d_common)/pbar.f \
          FFLAGS__dyn3d_common__pbar.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__divgrad.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

divgrad.done: \
          divgrad.o
	touch $(FCM_DONEDIR)/$@

divgrad.o: \
          $(PPSRCDIR0__dyn3d_common)/divgrad.f \
          FFLAGS__dyn3d_common__divgrad.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__nxgrad.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

nxgrad.done: \
          nxgrad.o
	touch $(FCM_DONEDIR)/$@

nxgrad.o: \
          $(PPSRCDIR0__dyn3d_common)/nxgrad.f \
          FFLAGS__dyn3d_common__nxgrad.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__cpdet_mod.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

cpdet_mod.done: \
          cpdet_mod.o \
          comconst_mod.done \
          control_mod.done
	touch $(FCM_DONEDIR)/$@

cpdet_mod.o: \
          $(PPSRCDIR0__dyn3d_common)/cpdet_mod.f90 \
          FFLAGS__dyn3d_common__cpdet_mod.flags \
          comconst_mod.o \
          control_mod.o
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__exner_hyb_m.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

exner_hyb_m.done: \
          exner_hyb_m.o \
          comconst_mod.done \
          comgeom.h.idone \
          comvert_mod.done \
          dimensions.h.idone \
          paramet.h.idone
	touch $(FCM_DONEDIR)/$@

exner_hyb_m.o: \
          $(PPSRCDIR0__dyn3d_common)/exner_hyb_m.f90 \
          FFLAGS__dyn3d_common__exner_hyb_m.flags \
          comconst_mod.o \
          comgeom.h \
          comvert_mod.o \
          dimensions.h \
          paramet.h
	fcm_internal compile:F dyn3d_common $< $@

clesph0.h: \
          $(SRCDIR0__dyn3d_common)/clesph0.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

clesph0.h.idone: \
          $(SRCDIR0__dyn3d_common)/clesph0.h
	touch $(FCM_DONEDIR)/$@

FFLAGS__dyn3d_common__gr_ecrit_fi.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

gr_ecrit_fi.done: \
          gr_ecrit_fi.o
	touch $(FCM_DONEDIR)/$@

gr_ecrit_fi.o: \
          $(PPSRCDIR0__dyn3d_common)/gr_ecrit_fi.f \
          FFLAGS__dyn3d_common__gr_ecrit_fi.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__limx.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

limx.done: \
          limx.o
	touch $(FCM_DONEDIR)/$@

limx.o: \
          $(PPSRCDIR0__dyn3d_common)/limx.f \
          FFLAGS__dyn3d_common__limx.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__gr_u_scal.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

gr_u_scal.done: \
          gr_u_scal.o
	touch $(FCM_DONEDIR)/$@

gr_u_scal.o: \
          $(PPSRCDIR0__dyn3d_common)/gr_u_scal.f \
          FFLAGS__dyn3d_common__gr_u_scal.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__defrun.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

defrun.done: \
          defrun.o \
          control_mod.done \
          logic_mod.done \
          serre_mod.done
	touch $(FCM_DONEDIR)/$@

defrun.o: \
          $(PPSRCDIR0__dyn3d_common)/defrun.f \
          FFLAGS__dyn3d_common__defrun.flags \
          control_mod.o \
          logic_mod.o \
          serre_mod.o
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__inidissip.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

inidissip.done: \
          inidissip.o \
          comconst_mod.done \
          comdissipn.h.idone \
          comvert_mod.done \
          control_mod.done \
          dimensions.h.idone \
          iniprint.h.idone \
          logic_mod.done \
          paramet.h.idone
	touch $(FCM_DONEDIR)/$@

inidissip.o: \
          $(PPSRCDIR0__dyn3d_common)/inidissip.f90 \
          FFLAGS__dyn3d_common__inidissip.flags \
          comconst_mod.o \
          comdissipn.h \
          comvert_mod.o \
          control_mod.o \
          dimensions.h \
          iniprint.h \
          logic_mod.o \
          paramet.h
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__advy.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

advy.done: \
          advy.o
	touch $(FCM_DONEDIR)/$@

advy.o: \
          $(PPSRCDIR0__dyn3d_common)/advy.f \
          FFLAGS__dyn3d_common__advy.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__initial0.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

initial0.done: \
          initial0.o
	touch $(FCM_DONEDIR)/$@

initial0.o: \
          $(PPSRCDIR0__dyn3d_common)/initial0.f \
          FFLAGS__dyn3d_common__initial0.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__bernoui.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

bernoui.done: \
          bernoui.o
	touch $(FCM_DONEDIR)/$@

bernoui.o: \
          $(PPSRCDIR0__dyn3d_common)/bernoui.f \
          FFLAGS__dyn3d_common__bernoui.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__inithist.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

inithist.done: \
          inithist.o \
          com_io_dyn_mod.done \
          comconst_mod.done \
          comvert_mod.done \
          infotrac.done \
          temps_mod.done
	touch $(FCM_DONEDIR)/$@

inithist.o: \
          $(PPSRCDIR0__dyn3d_common)/inithist.f \
          FFLAGS__dyn3d_common__inithist.flags \
          com_io_dyn_mod.o \
          comconst_mod.o \
          comvert_mod.o \
          infotrac.o \
          temps_mod.o
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__gr_int_dyn.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

gr_int_dyn.done: \
          gr_int_dyn.o
	touch $(FCM_DONEDIR)/$@

gr_int_dyn.o: \
          $(PPSRCDIR0__dyn3d_common)/gr_int_dyn.f \
          FFLAGS__dyn3d_common__gr_int_dyn.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__com_io_dyn_mod.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

com_io_dyn_mod.done: \
          com_io_dyn_mod.o
	touch $(FCM_DONEDIR)/$@

com_io_dyn_mod.o: \
          $(PPSRCDIR0__dyn3d_common)/com_io_dyn_mod.f90 \
          FFLAGS__dyn3d_common__com_io_dyn_mod.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__gradiv2.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

gradiv2.done: \
          gradiv2.o
	touch $(FCM_DONEDIR)/$@

gradiv2.o: \
          $(PPSRCDIR0__dyn3d_common)/gradiv2.f \
          FFLAGS__dyn3d_common__gradiv2.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__fxhyp_m.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

fxhyp_m.done: \
          fxhyp_m.o \
          arth_m.done \
          dimensions.h.idone \
          invert_zoom_x_m.done \
          nrtype.done \
          principal_cshift_m.done \
          serre_mod.done
	touch $(FCM_DONEDIR)/$@

fxhyp_m.o: \
          $(PPSRCDIR0__dyn3d_common)/fxhyp_m.f90 \
          FFLAGS__dyn3d_common__fxhyp_m.flags \
          arth_m.o \
          dimensions.h \
          invert_zoom_x_m.o \
          nrtype.o \
          principal_cshift_m.o \
          serre_mod.o
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__invert_zoom_x_m.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

invert_zoom_x_m.done: \
          invert_zoom_x_m.o \
          coefpoly_m.done \
          dimensions.h.idone \
          nrtype.done \
          serre_mod.done
	touch $(FCM_DONEDIR)/$@

invert_zoom_x_m.o: \
          $(PPSRCDIR0__dyn3d_common)/invert_zoom_x_m.f90 \
          FFLAGS__dyn3d_common__invert_zoom_x_m.flags \
          coefpoly_m.o \
          dimensions.h \
          nrtype.o \
          serre_mod.o
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__diagedyn.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

diagedyn.done: \
          diagedyn.o \
          control_mod.done \
          cpdet_mod.done
	touch $(FCM_DONEDIR)/$@

diagedyn.o: \
          $(PPSRCDIR0__dyn3d_common)/diagedyn.f \
          FFLAGS__dyn3d_common__diagedyn.flags \
          control_mod.o \
          cpdet_mod.o
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__rotatf.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

rotatf.done: \
          rotatf.o
	touch $(FCM_DONEDIR)/$@

rotatf.o: \
          $(PPSRCDIR0__dyn3d_common)/rotatf.f \
          FFLAGS__dyn3d_common__rotatf.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__pression.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

pression.done: \
          pression.o \
          dimensions.h.idone \
          paramet.h.idone
	touch $(FCM_DONEDIR)/$@

pression.o: \
          $(PPSRCDIR0__dyn3d_common)/pression.f90 \
          FFLAGS__dyn3d_common__pression.flags \
          dimensions.h \
          paramet.h
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__iso_verif_dyn.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

iso_verif_nonan_nostop.done: \
          iso_verif_nonan_nostop.o \
          infotrac.done
	touch $(FCM_DONEDIR)/$@

iso_verif_nonan_nostop.o: \
          $(PPSRCDIR0__dyn3d_common)/iso_verif_dyn.f \
          FFLAGS__dyn3d_common__iso_verif_dyn.flags \
          infotrac.o
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__advx.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

advx.done: \
          advx.o
	touch $(FCM_DONEDIR)/$@

advx.o: \
          $(PPSRCDIR0__dyn3d_common)/advx.f \
          FFLAGS__dyn3d_common__advx.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__sortvarc0.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

sortvarc0.done: \
          sortvarc0.o \
          comconst_mod.done \
          ener_mod.done
	touch $(FCM_DONEDIR)/$@

sortvarc0.o: \
          $(PPSRCDIR0__dyn3d_common)/sortvarc0.f \
          FFLAGS__dyn3d_common__sortvarc0.flags \
          comconst_mod.o \
          ener_mod.o
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__writehist.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

writehist.done: \
          writehist.o \
          com_io_dyn_mod.done \
          infotrac.done \
          temps_mod.done
	touch $(FCM_DONEDIR)/$@

writehist.o: \
          $(PPSRCDIR0__dyn3d_common)/writehist.f \
          FFLAGS__dyn3d_common__writehist.flags \
          com_io_dyn_mod.o \
          infotrac.o \
          temps_mod.o
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__prather.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

prather.done: \
          prather.o \
          comconst_mod.done
	touch $(FCM_DONEDIR)/$@

prather.o: \
          $(PPSRCDIR0__dyn3d_common)/prather.f \
          FFLAGS__dyn3d_common__prather.flags \
          comconst_mod.o
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__nxgrad_gam.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

nxgrad_gam.done: \
          nxgrad_gam.o
	touch $(FCM_DONEDIR)/$@

nxgrad_gam.o: \
          $(PPSRCDIR0__dyn3d_common)/nxgrad_gam.f \
          FFLAGS__dyn3d_common__nxgrad_gam.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__massbar.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

massbar.done: \
          massbar.o \
          comgeom.h.idone \
          dimensions.h.idone \
          paramet.h.idone
	touch $(FCM_DONEDIR)/$@

massbar.o: \
          $(PPSRCDIR0__dyn3d_common)/massbar.f90 \
          FFLAGS__dyn3d_common__massbar.flags \
          comgeom.h \
          dimensions.h \
          paramet.h
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__comconst_mod.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

comconst_mod.done: \
          comconst_mod.o
	touch $(FCM_DONEDIR)/$@

comconst_mod.o: \
          $(PPSRCDIR0__dyn3d_common)/comconst_mod.f90 \
          FFLAGS__dyn3d_common__comconst_mod.flags
	fcm_internal compile:F dyn3d_common $< $@

FFLAGS__dyn3d_common__traceurpole.flags: \
          FFLAGS__dyn3d_common.flags
	touch $(FCM_FLAGSDIR)/$@

traceurpole.done: \
          traceurpole.o \
          control_mod.done
	touch $(FCM_DONEDIR)/$@

traceurpole.o: \
          $(PPSRCDIR0__dyn3d_common)/traceurpole.f \
          FFLAGS__dyn3d_common__traceurpole.flags \
          control_mod.o
	fcm_internal compile:F dyn3d_common $< $@

comdissnew.h: \
          $(SRCDIR0__dyn3d_common)/comdissnew.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

comdissnew.h.idone: \
          $(SRCDIR0__dyn3d_common)/comdissnew.h
	touch $(FCM_DONEDIR)/$@

