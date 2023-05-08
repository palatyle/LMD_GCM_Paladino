#

LN      =       ln -sf
MAKE    =       make -i -r
RM      =       rm -f


MODULES = \
	module_physics_addtendc.o \
	module_physics_init.o \
	variables_mod.o \
	update_inputs_physiq_mod.o \
	update_outputs_physiq_mod.o \
	iniphysiq_mod.o \
	callphysiq_mod.o \
	module_lmd_driver.o 
 
OBJS    =

NMM_MODULES =  

LIBTARGET    =  physics
TARGETDIR    =  ./

$(LIBTARGET) :
		if [ $(WRF_NMM_CORE) -eq 1 ] ; then \
                  $(MAKE) nmm_contrib                                      ; \
                  $(AR) ../main/libwrflib.a $(MODULES) $(OBJS) $(NMM_OBJS) $(NMM_MODULES) ; \
                else                                                         \
                  $(MAKE) non_nmm                                          ; \
                  $(AR) ../main/libwrflib.a $(MODULES) $(OBJS)             ; \
                fi

include ../configure.wrf

nmm_contrib : $(NMM_OBJS) $(NMM_MODULES) $(MODULES) $(OBJS)

non_nmm : $(MODULES) $(OBJS)

clean:
	@ echo 'use the clean script'

# DEPENDENCIES : only dependencies after this line (don't remove the word DEPENDENCIES)
# si ce qui est à droite a changé, on recompile à gauche

module_physics_addtendc.o: \
		../frame/module_state_description.o \
		../frame/module_configure.o

module_physics_init.o : \
		../frame/module_state_description.o \
		../frame/module_configure.o \
		../frame/module_wrf_error.o \
		../frame/module_dm.o \
		../share/module_model_constants.o 

module_lmd_driver.o: \
                ../frame/module_state_description.o \
                ../frame/module_wrf_error.o \
                ../frame/module_configure.o \
                ../share/module_model_constants.o \

#               ../frame/module_wrf_error.o \
#               ../share/module_model_constants.o

# DO NOT DELETE
