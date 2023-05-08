# Automatic Make rule for filtrez

SRCDIR0__filtrez = /home/palatyle/LMD_gen/trunk/LMDZ.COMMON/libf/filtrez

PPSRCDIR0__filtrez = /home/palatyle/LMD_gen/trunk/LMDZ.COMMON/libo/sawtooth_128x96x23_phystd_para.e/.config/ppsrc/filtrez

filtrez.etc : \
          $(SRCDIR0__filtrez)/parafilt.h_96x71x19 \
          $(SRCDIR0__filtrez)/parafilt.h_192x142x29 \
          $(FCM_DONEDIR)/FCM_CP.dummy
	cp $^ $(FCM_ETCDIR)
	touch $(FCM_DONEDIR)/$@

FFLAGS__filtrez__mod_fft_mkl.flags: \
          FFLAGS__filtrez.flags
	touch $(FCM_FLAGSDIR)/$@

mod_fft_mkl.done: \
          mod_fft_mkl.o
	touch $(FCM_DONEDIR)/$@

mod_fft_mkl.o: \
          $(PPSRCDIR0__filtrez)/mod_fft_mkl.f90 \
          FFLAGS__filtrez__mod_fft_mkl.flags
	fcm_internal compile:F filtrez $< $@

FFLAGS__filtrez__mod_fft_wrapper.flags: \
          FFLAGS__filtrez.flags
	touch $(FCM_FLAGSDIR)/$@

mod_fft_wrapper.done: \
          mod_fft_wrapper.o
	touch $(FCM_DONEDIR)/$@

mod_fft_wrapper.o: \
          $(PPSRCDIR0__filtrez)/mod_fft_wrapper.f90 \
          FFLAGS__filtrez__mod_fft_wrapper.flags
	fcm_internal compile:F filtrez $< $@

FFLAGS__filtrez__mod_fft_mathkeisan.flags: \
          FFLAGS__filtrez.flags
	touch $(FCM_FLAGSDIR)/$@

mod_fft_mathkeisan.done: \
          mod_fft_mathkeisan.o
	touch $(FCM_DONEDIR)/$@

mod_fft_mathkeisan.o: \
          $(PPSRCDIR0__filtrez)/mod_fft_mathkeisan.f90 \
          FFLAGS__filtrez__mod_fft_mathkeisan.flags
	fcm_internal compile:F filtrez $< $@

FFLAGS__filtrez__filtreg.flags: \
          FFLAGS__filtrez.flags
	touch $(FCM_FLAGSDIR)/$@

filtreg.done: \
          filtreg.o \
          filtreg_mod.done
	touch $(FCM_DONEDIR)/$@

filtreg.o: \
          $(PPSRCDIR0__filtrez)/filtreg.f \
          FFLAGS__filtrez__filtreg.flags \
          filtreg_mod.o
	fcm_internal compile:F filtrez $< $@

FFLAGS__filtrez__inifgn.flags: \
          FFLAGS__filtrez.flags
	touch $(FCM_FLAGSDIR)/$@

inifgn.done: \
          inifgn.o
	touch $(FCM_DONEDIR)/$@

inifgn.o: \
          $(PPSRCDIR0__filtrez)/inifgn.f \
          FFLAGS__filtrez__inifgn.flags
	fcm_internal compile:F filtrez $< $@

FFLAGS__filtrez__mod_fft_fftw.flags: \
          FFLAGS__filtrez.flags
	touch $(FCM_FLAGSDIR)/$@

mod_fft_fftw.done: \
          mod_fft_fftw.o
	touch $(FCM_DONEDIR)/$@

mod_fft_fftw.o: \
          $(PPSRCDIR0__filtrez)/mod_fft_fftw.f90 \
          FFLAGS__filtrez__mod_fft_fftw.flags
	fcm_internal compile:F filtrez $< $@

FFLAGS__filtrez__filtreg_mod.flags: \
          FFLAGS__filtrez.flags
	touch $(FCM_FLAGSDIR)/$@

filtreg_mod.done: \
          filtreg_mod.o \
          logic_mod.done \
          mod_filtre_fft.done \
          mod_filtre_fft_loc.done \
          serre_mod.done
	touch $(FCM_DONEDIR)/$@

filtreg_mod.o: \
          $(PPSRCDIR0__filtrez)/filtreg_mod.f90 \
          FFLAGS__filtrez__filtreg_mod.flags \
          logic_mod.o \
          mod_filtre_fft.o \
          mod_filtre_fft_loc.o \
          serre_mod.o
	fcm_internal compile:F filtrez $< $@

FFLAGS__filtrez__acc.flags: \
          FFLAGS__filtrez.flags
	touch $(FCM_FLAGSDIR)/$@

acc.done: \
          acc.o
	touch $(FCM_DONEDIR)/$@

acc.o: \
          $(PPSRCDIR0__filtrez)/acc.f \
          FFLAGS__filtrez__acc.flags
	fcm_internal compile:F filtrez $< $@

FFLAGS__filtrez__eigen_sort.flags: \
          FFLAGS__filtrez.flags
	touch $(FCM_FLAGSDIR)/$@

eigen_sort.done: \
          eigen_sort.o
	touch $(FCM_DONEDIR)/$@

eigen_sort.o: \
          $(PPSRCDIR0__filtrez)/eigen_sort.f \
          FFLAGS__filtrez__eigen_sort.flags
	fcm_internal compile:F filtrez $< $@

FFLAGS__filtrez__mkl_dft_type.flags: \
          FFLAGS__filtrez.flags
	touch $(FCM_FLAGSDIR)/$@

mkl_dft_type.done: \
          mkl_dft_type.o
	touch $(FCM_DONEDIR)/$@

mkl_dft_type.o: \
          $(SRCDIR0__filtrez)/mkl_dft_type.f90 \
          FFLAGS__filtrez__mkl_dft_type.flags
	fcm_internal compile:F filtrez $< $@

FFLAGS__filtrez__jacobi.flags: \
          FFLAGS__filtrez.flags
	touch $(FCM_FLAGSDIR)/$@

jacobi.done: \
          jacobi.o
	touch $(FCM_DONEDIR)/$@

jacobi.o: \
          $(PPSRCDIR0__filtrez)/jacobi.f90 \
          FFLAGS__filtrez__jacobi.flags
	fcm_internal compile:F filtrez $< $@

FFLAGS__filtrez__mkl_dfti.flags: \
          FFLAGS__filtrez.flags
	touch $(FCM_FLAGSDIR)/$@

mkl_dfti.done: \
          mkl_dfti.o \
          mkl_dft_type.done
	touch $(FCM_DONEDIR)/$@

mkl_dfti.o: \
          $(SRCDIR0__filtrez)/mkl_dfti.f90 \
          FFLAGS__filtrez__mkl_dfti.flags \
          mkl_dft_type.o
	fcm_internal compile:F filtrez $< $@

FFLAGS__filtrez__mod_fft.flags: \
          FFLAGS__filtrez.flags
	touch $(FCM_FLAGSDIR)/$@

mod_fft.done: \
          mod_fft.o \
          mod_fft_wrapper.done
	touch $(FCM_DONEDIR)/$@

mod_fft.o: \
          $(PPSRCDIR0__filtrez)/mod_fft.f90 \
          FFLAGS__filtrez__mod_fft.flags \
          mod_fft_wrapper.o
	fcm_internal compile:F filtrez $< $@

FFLAGS__filtrez__eigen.flags: \
          FFLAGS__filtrez.flags
	touch $(FCM_FLAGSDIR)/$@

eigen.done: \
          eigen.o
	touch $(FCM_DONEDIR)/$@

eigen.o: \
          $(PPSRCDIR0__filtrez)/eigen.f \
          FFLAGS__filtrez__eigen.flags
	fcm_internal compile:F filtrez $< $@

FFLAGS__filtrez__mod_filtre_fft_loc.flags: \
          FFLAGS__filtrez.flags
	touch $(FCM_FLAGSDIR)/$@

mod_filtre_fft_loc.done: \
          mod_filtre_fft_loc.o \
          dimensions.h.idone \
          mod_fft.done \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

mod_filtre_fft_loc.o: \
          $(PPSRCDIR0__filtrez)/mod_filtre_fft_loc.f90 \
          FFLAGS__filtrez__mod_filtre_fft_loc.flags \
          dimensions.h \
          mod_fft.o \
          parallel_lmdz.o
	fcm_internal compile:F filtrez $< $@

FFLAGS__filtrez__timer_filtre.flags: \
          FFLAGS__filtrez.flags
	touch $(FCM_FLAGSDIR)/$@

timer_filtre.done: \
          timer_filtre.o
	touch $(FCM_DONEDIR)/$@

timer_filtre.o: \
          $(PPSRCDIR0__filtrez)/timer_filtre.f90 \
          FFLAGS__filtrez__timer_filtre.flags
	fcm_internal compile:F filtrez $< $@

FFLAGS__filtrez__mod_filtre_fft.flags: \
          FFLAGS__filtrez.flags
	touch $(FCM_FLAGSDIR)/$@

mod_filtre_fft.done: \
          mod_filtre_fft.o \
          dimensions.h.idone \
          mod_fft.done \
          parallel_lmdz.done
	touch $(FCM_DONEDIR)/$@

mod_filtre_fft.o: \
          $(PPSRCDIR0__filtrez)/mod_filtre_fft.f90 \
          FFLAGS__filtrez__mod_filtre_fft.flags \
          dimensions.h \
          mod_fft.o \
          parallel_lmdz.o
	fcm_internal compile:F filtrez $< $@

coefils.h: \
          $(SRCDIR0__filtrez)/coefils.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

coefils.h.idone: \
          $(SRCDIR0__filtrez)/coefils.h
	touch $(FCM_DONEDIR)/$@

parafilt.h: \
          $(SRCDIR0__filtrez)/parafilt.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

parafilt.h.idone: \
          $(SRCDIR0__filtrez)/parafilt.h
	touch $(FCM_DONEDIR)/$@

