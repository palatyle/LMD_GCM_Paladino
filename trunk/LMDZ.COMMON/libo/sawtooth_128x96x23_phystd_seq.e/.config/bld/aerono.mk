# Automatic Make rule for aerono

SRCDIR0__aerono = /home/palatyle/LMD_gen/trunk/LMDZ.COMMON/libf/aeronostd

PPSRCDIR0__aerono = /home/palatyle/LMD_gen/trunk/LMDZ.COMMON/libo/sawtooth_128x96x23_phystd_seq.e/.config/ppsrc/aerono

chimiedata.h: \
          $(SRCDIR0__aerono)/chimiedata.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

chimiedata.h.idone: \
          $(SRCDIR0__aerono)/chimiedata.h
	touch $(FCM_DONEDIR)/$@

FFLAGS__aerono__calchim_asis.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

calchim_asis.done: \
          calchim_asis.o \
          callkeys_mod.done \
          comcstfi_mod.done \
          conc_mod.done \
          tracer_h.done
	touch $(FCM_DONEDIR)/$@

calchim_asis.o: \
          $(PPSRCDIR0__aerono)/calchim_asis.f90 \
          FFLAGS__aerono__calchim_asis.flags \
          callkeys_mod.o \
          comcstfi_mod.o \
          conc_mod.o \
          tracer_h.o
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__photolysis_asis.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

photolysis_asis.done: \
          photolysis_asis.o \
          callkeys_mod.done \
          comcstfi_mod.done
	touch $(FCM_DONEDIR)/$@

photolysis_asis.o: \
          $(PPSRCDIR0__aerono)/photolysis_asis.f90 \
          FFLAGS__aerono__photolysis_asis.flags \
          callkeys_mod.o \
          comcstfi_mod.o
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__photochemistry_asis.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

photochemistry_asis.done: \
          photochemistry_asis.o \
          callkeys_mod.done \
          comcstfi_mod.done \
          tracer_h.done \
          types_asis.done
	touch $(FCM_DONEDIR)/$@

photochemistry_asis.o: \
          $(PPSRCDIR0__aerono)/photochemistry_asis.f90 \
          FFLAGS__aerono__photochemistry_asis.flags \
          callkeys_mod.o \
          comcstfi_mod.o \
          tracer_h.o \
          types_asis.o
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__inv.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

inv.done: \
          inv.o
	touch $(FCM_DONEDIR)/$@

inv.o: \
          $(PPSRCDIR0__aerono)/inv.f \
          FFLAGS__aerono__inv.flags
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__intrplf.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

intrplf.done: \
          intrplf.o
	touch $(FCM_DONEDIR)/$@

intrplf.o: \
          $(PPSRCDIR0__aerono)/intrplf.f \
          FFLAGS__aerono__intrplf.flags
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__read_phototable.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

read_phototable.done: \
          read_phototable.o \
          datafile_mod.done \
          ioipsl_getin_p_mod.done \
          ioipsl_getincom.done
	touch $(FCM_DONEDIR)/$@

read_phototable.o: \
          $(PPSRCDIR0__aerono)/read_phototable.f90 \
          FFLAGS__aerono__read_phototable.flags \
          datafile_mod.o \
          ioipsl_getin_p_mod.o \
          ioipsl_getincom.o
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__types_asis.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

types_asis.done: \
          types_asis.o
	touch $(FCM_DONEDIR)/$@

types_asis.o: \
          $(PPSRCDIR0__aerono)/types_asis.f90 \
          FFLAGS__aerono__types_asis.flags
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__concentrations.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

concentrations.done: \
          concentrations.o \
          callkeys_mod.done \
          comcstfi_mod.done \
          conc_mod.done \
          tracer_h.done
	touch $(FCM_DONEDIR)/$@

concentrations.o: \
          $(PPSRCDIR0__aerono)/concentrations.f \
          FFLAGS__aerono__concentrations.flags \
          callkeys_mod.o \
          comcstfi_mod.o \
          conc_mod.o \
          tracer_h.o
	fcm_internal compile:F aerono $< $@

