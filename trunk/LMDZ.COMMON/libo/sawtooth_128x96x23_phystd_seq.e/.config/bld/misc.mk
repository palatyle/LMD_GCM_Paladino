# Automatic Make rule for misc

SRCDIR0__misc = /home/palatyle/LMD_gen/trunk/LMDZ.COMMON/libf/misc

PPSRCDIR0__misc = /home/palatyle/LMD_gen/trunk/LMDZ.COMMON/libo/sawtooth_128x96x23_phystd_seq.e/.config/ppsrc/misc

FFLAGS__misc__xgetua.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

xgetua.done: \
          xgetua.o
	touch $(FCM_DONEDIR)/$@

xgetua.o: \
          $(PPSRCDIR0__misc)/xgetua.f \
          FFLAGS__misc__xgetua.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__ran1.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

ran1.done: \
          ran1.o
	touch $(FCM_DONEDIR)/$@

ran1.o: \
          $(PPSRCDIR0__misc)/ran1.f \
          FFLAGS__misc__ran1.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__nf95_def_var_m.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

nf95_def_var_m.done: \
          nf95_def_var_m.o \
          handle_err_m.done
	touch $(FCM_DONEDIR)/$@

nf95_def_var_m.o: \
          $(PPSRCDIR0__misc)/nf95_def_var_m.f90 \
          FFLAGS__misc__nf95_def_var_m.flags \
          handle_err_m.o
	fcm_internal compile:F misc $< $@

FFLAGS__misc__regr3_lint_m.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

regr3_lint_m.done: \
          regr3_lint_m.o \
          assert_eq_m.done \
          interpolation.done
	touch $(FCM_DONEDIR)/$@

regr3_lint_m.o: \
          $(PPSRCDIR0__misc)/regr3_lint_m.f90 \
          FFLAGS__misc__regr3_lint_m.flags \
          assert_eq_m.o \
          interpolation.o
	fcm_internal compile:F misc $< $@

FFLAGS__misc__xerhlt.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

xerhlt.done: \
          xerhlt.o
	touch $(FCM_DONEDIR)/$@

xerhlt.o: \
          $(PPSRCDIR0__misc)/xerhlt.f \
          FFLAGS__misc__xerhlt.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__xersve.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

xersve.done: \
          xersve.o
	touch $(FCM_DONEDIR)/$@

xersve.o: \
          $(PPSRCDIR0__misc)/xersve.f \
          FFLAGS__misc__xersve.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__arth_m.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

arth_m.done: \
          arth_m.o
	touch $(FCM_DONEDIR)/$@

arth_m.o: \
          $(PPSRCDIR0__misc)/arth_m.f90 \
          FFLAGS__misc__arth_m.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__nf95_get_att_m.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

nf95_get_att_m.done: \
          nf95_get_att_m.o \
          handle_err_m.done \
          simple.done
	touch $(FCM_DONEDIR)/$@

nf95_get_att_m.o: \
          $(PPSRCDIR0__misc)/nf95_get_att_m.f90 \
          FFLAGS__misc__nf95_get_att_m.flags \
          handle_err_m.o \
          simple.o
	fcm_internal compile:F misc $< $@

FFLAGS__misc__regr1_conserv_m.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

regr1_conserv_m.done: \
          regr1_conserv_m.o \
          assert_eq_m.done \
          assert_m.done \
          interpolation.done
	touch $(FCM_DONEDIR)/$@

regr1_conserv_m.o: \
          $(PPSRCDIR0__misc)/regr1_conserv_m.f90 \
          FFLAGS__misc__regr1_conserv_m.flags \
          assert_eq_m.o \
          assert_m.o \
          interpolation.o
	fcm_internal compile:F misc $< $@

FFLAGS__misc__j4save.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

j4save.done: \
          j4save.o
	touch $(FCM_DONEDIR)/$@

j4save.o: \
          $(PPSRCDIR0__misc)/j4save.f \
          FFLAGS__misc__j4save.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__xercnt.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

xercnt.done: \
          xercnt.o
	touch $(FCM_DONEDIR)/$@

xercnt.o: \
          $(PPSRCDIR0__misc)/xercnt.f \
          FFLAGS__misc__xercnt.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__ioipsl_stringop.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

ioipsl_stringop.done: \
          ioipsl_stringop.o
	touch $(FCM_DONEDIR)/$@

ioipsl_stringop.o: \
          $(PPSRCDIR0__misc)/ioipsl_stringop.f90 \
          FFLAGS__misc__ioipsl_stringop.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__cray.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

scopy.done: \
          scopy.o
	touch $(FCM_DONEDIR)/$@

scopy.o: \
          $(PPSRCDIR0__misc)/cray.f \
          FFLAGS__misc__cray.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__new_unit_m.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

new_unit_m.done: \
          new_unit_m.o
	touch $(FCM_DONEDIR)/$@

new_unit_m.o: \
          $(PPSRCDIR0__misc)/new_unit_m.f90 \
          FFLAGS__misc__new_unit_m.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__minmax.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

minmax.done: \
          minmax.o
	touch $(FCM_DONEDIR)/$@

minmax.o: \
          $(PPSRCDIR0__misc)/minmax.f \
          FFLAGS__misc__minmax.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__pres2lev_mod.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

pres2lev_mod.done: \
          pres2lev_mod.o
	touch $(FCM_DONEDIR)/$@

pres2lev_mod.o: \
          $(PPSRCDIR0__misc)/pres2lev_mod.f90 \
          FFLAGS__misc__pres2lev_mod.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__slopes_m.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

slopes_m.done: \
          slopes_m.o
	touch $(FCM_DONEDIR)/$@

slopes_m.o: \
          $(PPSRCDIR0__misc)/slopes_m.f90 \
          FFLAGS__misc__slopes_m.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__i1mach.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

i1mach.done: \
          i1mach.o
	touch $(FCM_DONEDIR)/$@

i1mach.o: \
          $(PPSRCDIR0__misc)/i1mach.f \
          FFLAGS__misc__i1mach.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__netcdf95.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

netcdf95.done: \
          netcdf95.o \
          handle_err_m.done \
          nf95_def_var_m.done \
          nf95_get_att_m.done \
          nf95_get_var_m.done \
          nf95_gw_var_m.done \
          nf95_put_att_m.done \
          nf95_put_var_m.done \
          simple.done
	touch $(FCM_DONEDIR)/$@

netcdf95.o: \
          $(PPSRCDIR0__misc)/netcdf95.f90 \
          FFLAGS__misc__netcdf95.flags \
          handle_err_m.o \
          nf95_def_var_m.o \
          nf95_get_att_m.o \
          nf95_get_var_m.o \
          nf95_gw_var_m.o \
          nf95_put_att_m.o \
          nf95_put_var_m.o \
          simple.o
	fcm_internal compile:F misc $< $@

FFLAGS__misc__coefpoly_m.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

coefpoly_m.done: \
          coefpoly_m.o \
          nrtype.done
	touch $(FCM_DONEDIR)/$@

coefpoly_m.o: \
          $(PPSRCDIR0__misc)/coefpoly_m.f90 \
          FFLAGS__misc__coefpoly_m.flags \
          nrtype.o
	fcm_internal compile:F misc $< $@

FFLAGS__misc__vampir.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

vampir.done: \
          vampir.o
	touch $(FCM_DONEDIR)/$@

vampir.o: \
          $(PPSRCDIR0__misc)/vampir.f90 \
          FFLAGS__misc__vampir.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__nf95_put_att_m.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

nf95_put_att_m.done: \
          nf95_put_att_m.o \
          handle_err_m.done
	touch $(FCM_DONEDIR)/$@

nf95_put_att_m.o: \
          $(PPSRCDIR0__misc)/nf95_put_att_m.f90 \
          FFLAGS__misc__nf95_put_att_m.flags \
          handle_err_m.o
	fcm_internal compile:F misc $< $@

FFLAGS__misc__interpolation.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

interpolation.done: \
          interpolation.o
	touch $(FCM_DONEDIR)/$@

interpolation.o: \
          $(PPSRCDIR0__misc)/interpolation.f90 \
          FFLAGS__misc__interpolation.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__xermsg.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

xermsg.done: \
          xermsg.o
	touch $(FCM_DONEDIR)/$@

xermsg.o: \
          $(PPSRCDIR0__misc)/xermsg.f \
          FFLAGS__misc__xermsg.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__chfev.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

chfev.done: \
          chfev.o
	touch $(FCM_DONEDIR)/$@

chfev.o: \
          $(PPSRCDIR0__misc)/chfev.f \
          FFLAGS__misc__chfev.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__ioipsl_errioipsl.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

ioipsl_errioipsl.done: \
          ioipsl_errioipsl.o
	touch $(FCM_DONEDIR)/$@

ioipsl_errioipsl.o: \
          $(PPSRCDIR0__misc)/ioipsl_errioipsl.f90 \
          FFLAGS__misc__ioipsl_errioipsl.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__assert_m.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

assert_m.done: \
          assert_m.o
	touch $(FCM_DONEDIR)/$@

assert_m.o: \
          $(PPSRCDIR0__misc)/assert_m.f90 \
          FFLAGS__misc__assert_m.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__pchfe.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

pchfe.done: \
          pchfe.o
	touch $(FCM_DONEDIR)/$@

pchfe.o: \
          $(PPSRCDIR0__misc)/pchfe.f \
          FFLAGS__misc__pchfe.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__nf95_get_var_m.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

nf95_get_var_m.done: \
          nf95_get_var_m.o \
          handle_err_m.done
	touch $(FCM_DONEDIR)/$@

nf95_get_var_m.o: \
          $(PPSRCDIR0__misc)/nf95_get_var_m.f90 \
          FFLAGS__misc__nf95_get_var_m.flags \
          handle_err_m.o
	fcm_internal compile:F misc $< $@

FFLAGS__misc__regr1_lint_m.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

regr1_lint_m.done: \
          regr1_lint_m.o \
          assert_eq_m.done \
          interpolation.done
	touch $(FCM_DONEDIR)/$@

regr1_lint_m.o: \
          $(PPSRCDIR0__misc)/regr1_lint_m.f90 \
          FFLAGS__misc__regr1_lint_m.flags \
          assert_eq_m.o \
          interpolation.o
	fcm_internal compile:F misc $< $@

FFLAGS__misc__ismin.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

ismin.done: \
          ismin.o
	touch $(FCM_DONEDIR)/$@

ismin.o: \
          $(PPSRCDIR0__misc)/ismin.f \
          FFLAGS__misc__ismin.flags
	fcm_internal compile:F misc $< $@

iniprint.h: \
          $(SRCDIR0__misc)/iniprint.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

iniprint.h.idone: \
          $(SRCDIR0__misc)/iniprint.h
	touch $(FCM_DONEDIR)/$@

FFLAGS__misc__sort.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

sort.done: \
          sort.o
	touch $(FCM_DONEDIR)/$@

sort.o: \
          $(PPSRCDIR0__misc)/sort.f \
          FFLAGS__misc__sort.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__pchsp_95_m.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

pchsp_95_m.done: \
          pchsp_95_m.o \
          assert_eq_m.done
	touch $(FCM_DONEDIR)/$@

pchsp_95_m.o: \
          $(PPSRCDIR0__misc)/pchsp_95_m.f90 \
          FFLAGS__misc__pchsp_95_m.flags \
          assert_eq_m.o
	fcm_internal compile:F misc $< $@

FFLAGS__misc__nf95_gw_var_m.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

nf95_gw_var_m.done: \
          nf95_gw_var_m.o \
          nf95_get_var_m.done \
          simple.done
	touch $(FCM_DONEDIR)/$@

nf95_gw_var_m.o: \
          $(PPSRCDIR0__misc)/nf95_gw_var_m.f90 \
          FFLAGS__misc__nf95_gw_var_m.flags \
          nf95_get_var_m.o \
          simple.o
	fcm_internal compile:F misc $< $@

FFLAGS__misc__ioipsl_getincom.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

ioipsl_getincom.done: \
          ioipsl_getincom.o \
          ioipsl_errioipsl.done \
          ioipsl_stringop.done
	touch $(FCM_DONEDIR)/$@

ioipsl_getincom.o: \
          $(PPSRCDIR0__misc)/ioipsl_getincom.f90 \
          FFLAGS__misc__ioipsl_getincom.flags \
          ioipsl_errioipsl.o \
          ioipsl_stringop.o
	fcm_internal compile:F misc $< $@

FFLAGS__misc__q_sat.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

q_sat.done: \
          q_sat.o
	touch $(FCM_DONEDIR)/$@

q_sat.o: \
          $(PPSRCDIR0__misc)/q_sat.f \
          FFLAGS__misc__q_sat.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__xerprn.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

xerprn.done: \
          xerprn.o
	touch $(FCM_DONEDIR)/$@

xerprn.o: \
          $(PPSRCDIR0__misc)/xerprn.f \
          FFLAGS__misc__xerprn.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__write_field.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

write_field.done: \
          write_field.o
	touch $(FCM_DONEDIR)/$@

write_field.o: \
          $(PPSRCDIR0__misc)/write_field.f90 \
          FFLAGS__misc__write_field.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__ismax.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

ismax.done: \
          ismax.o
	touch $(FCM_DONEDIR)/$@

ismax.o: \
          $(PPSRCDIR0__misc)/ismax.f \
          FFLAGS__misc__ismax.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__pchfe_95_m.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

pchfe_95_m.done: \
          pchfe_95_m.o \
          assert_eq_m.done
	touch $(FCM_DONEDIR)/$@

pchfe_95_m.o: \
          $(PPSRCDIR0__misc)/pchfe_95_m.f90 \
          FFLAGS__misc__pchfe_95_m.flags \
          assert_eq_m.o
	fcm_internal compile:F misc $< $@

FFLAGS__misc__minmax2.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

minmax2.done: \
          minmax2.o
	touch $(FCM_DONEDIR)/$@

minmax2.o: \
          $(PPSRCDIR0__misc)/minmax2.f \
          FFLAGS__misc__minmax2.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__pchdf.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

pchdf.done: \
          pchdf.o
	touch $(FCM_DONEDIR)/$@

pchdf.o: \
          $(PPSRCDIR0__misc)/pchdf.f \
          FFLAGS__misc__pchdf.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__pchsp.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

pchsp.done: \
          pchsp.o
	touch $(FCM_DONEDIR)/$@

pchsp.o: \
          $(PPSRCDIR0__misc)/pchsp.f \
          FFLAGS__misc__pchsp.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__nrtype.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

nrtype.done: \
          nrtype.o
	touch $(FCM_DONEDIR)/$@

nrtype.o: \
          $(PPSRCDIR0__misc)/nrtype.f90 \
          FFLAGS__misc__nrtype.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__fdump.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

fdump.done: \
          fdump.o
	touch $(FCM_DONEDIR)/$@

fdump.o: \
          $(PPSRCDIR0__misc)/fdump.f \
          FFLAGS__misc__fdump.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__handle_err_m.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

handle_err_m.done: \
          handle_err_m.o
	touch $(FCM_DONEDIR)/$@

handle_err_m.o: \
          $(PPSRCDIR0__misc)/handle_err_m.f90 \
          FFLAGS__misc__handle_err_m.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__nf95_put_var_m.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

nf95_put_var_m.done: \
          nf95_put_var_m.o \
          handle_err_m.done
	touch $(FCM_DONEDIR)/$@

nf95_put_var_m.o: \
          $(PPSRCDIR0__misc)/nf95_put_var_m.f90 \
          FFLAGS__misc__nf95_put_var_m.flags \
          handle_err_m.o
	fcm_internal compile:F misc $< $@

FFLAGS__misc__cbrt.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

cbrt.done: \
          cbrt.o
	touch $(FCM_DONEDIR)/$@

cbrt.o: \
          $(PPSRCDIR0__misc)/cbrt.f \
          FFLAGS__misc__cbrt.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__assert_eq_m.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

assert_eq_m.done: \
          assert_eq_m.o
	touch $(FCM_DONEDIR)/$@

assert_eq_m.o: \
          $(PPSRCDIR0__misc)/assert_eq_m.f90 \
          FFLAGS__misc__assert_eq_m.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__formcoord.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

formcoord.done: \
          formcoord.o
	touch $(FCM_DONEDIR)/$@

formcoord.o: \
          $(PPSRCDIR0__misc)/formcoord.f \
          FFLAGS__misc__formcoord.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__juldate.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

juldate.done: \
          juldate.o
	touch $(FCM_DONEDIR)/$@

juldate.o: \
          $(PPSRCDIR0__misc)/juldate.f \
          FFLAGS__misc__juldate.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__sort_mod.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

sort_mod.done: \
          sort_mod.o
	touch $(FCM_DONEDIR)/$@

sort_mod.o: \
          $(PPSRCDIR0__misc)/sort_mod.f90 \
          FFLAGS__misc__sort_mod.flags
	fcm_internal compile:F misc $< $@

FFLAGS__misc__simple.flags: \
          FFLAGS__misc.flags
	touch $(FCM_FLAGSDIR)/$@

simple.done: \
          simple.o \
          handle_err_m.done
	touch $(FCM_DONEDIR)/$@

simple.o: \
          $(PPSRCDIR0__misc)/simple.f90 \
          FFLAGS__misc__simple.flags \
          handle_err_m.o
	fcm_internal compile:F misc $< $@

