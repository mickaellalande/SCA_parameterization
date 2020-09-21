# Automatic Make rule for sechiba

PPSRCDIR0__sechiba = /gpfswork/rech/goe/ufz23bm/SCA_parameterization/modipsl/modeles/ORCHIDEE/.config/ppsrc/sechiba

SRCDIR0__sechiba = /gpfswork/rech/goe/ufz23bm/SCA_parameterization/modipsl/modeles/ORCHIDEE/src_sechiba

sechiba.etc : \
          $(SRCDIR0__sechiba)/Makefile \
          $(SRCDIR0__sechiba)/AA_make.ldef \
          $(SRCDIR0__sechiba)/AA_make \
          $(FCM_DONEDIR)/FCM_CP.dummy
	cp $^ $(FCM_ETCDIR)
	touch $(FCM_DONEDIR)/$@

FFLAGS__sechiba__routing.flags: \
          FFLAGS__sechiba.flags
	touch $(FCM_FLAGSDIR)/$@

routing.done: \
          routing.o \
          constantes.done \
          constantes_soil.done \
          grid.done \
          interpol_help.done \
          ioipsl_para.done \
          mod_orchidee_para.done \
          pft_parameters.done \
          sechiba_io_p.done \
          time.done \
          xios_orchidee.done
	touch $(FCM_DONEDIR)/$@

routing.o: \
          $(PPSRCDIR0__sechiba)/routing.f90 \
          FFLAGS__sechiba__routing.flags \
          constantes.o \
          constantes_soil.o \
          grid.o \
          interpol_help.o \
          ioipsl_para.o \
          mod_orchidee_para.o \
          pft_parameters.o \
          sechiba_io_p.o \
          time.o \
          xios_orchidee.o
	fcm_internal compile:F sechiba $< $@

FFLAGS__sechiba__hydrolc.flags: \
          FFLAGS__sechiba.flags
	touch $(FCM_FLAGSDIR)/$@

hydrolc.done: \
          hydrolc.o \
          constantes.done \
          constantes_soil.done \
          explicitsnow.done \
          grid.done \
          ioipsl_para.done \
          mod_orchidee_para.done \
          pft_parameters.done \
          sechiba_io_p.done \
          time.done \
          xios_orchidee.done
	touch $(FCM_DONEDIR)/$@

hydrolc.o: \
          $(PPSRCDIR0__sechiba)/hydrolc.f90 \
          FFLAGS__sechiba__hydrolc.flags \
          constantes.o \
          constantes_soil.o \
          explicitsnow.o \
          grid.o \
          ioipsl_para.o \
          mod_orchidee_para.o \
          pft_parameters.o \
          sechiba_io_p.o \
          time.o \
          xios_orchidee.o
	fcm_internal compile:F sechiba $< $@

FFLAGS__sechiba__slowproc.flags: \
          FFLAGS__sechiba.flags
	touch $(FCM_FLAGSDIR)/$@

slowproc.done: \
          slowproc.o \
          constantes.done \
          constantes_soil.done \
          grid.done \
          interpol_help.done \
          interpweight.done \
          ioipsl_para.done \
          mod_orchidee_para.done \
          pft_parameters.done \
          sechiba_io_p.done \
          stomate.done \
          stomate_data.done \
          time.done \
          xios_orchidee.done
	touch $(FCM_DONEDIR)/$@

slowproc.o: \
          $(PPSRCDIR0__sechiba)/slowproc.f90 \
          FFLAGS__sechiba__slowproc.flags \
          constantes.o \
          constantes_soil.o \
          grid.o \
          interpol_help.o \
          interpweight.o \
          ioipsl_para.o \
          mod_orchidee_para.o \
          pft_parameters.o \
          sechiba_io_p.o \
          stomate.o \
          stomate_data.o \
          time.o \
          xios_orchidee.o
	fcm_internal compile:F sechiba $< $@

FFLAGS__sechiba__explicitsnow.flags: \
          FFLAGS__sechiba.flags
	touch $(FCM_FLAGSDIR)/$@

explicitsnow.done: \
          explicitsnow.o \
          constantes.done \
          constantes_soil.done \
          ioipsl_para.done \
          pft_parameters.done \
          qsat_moisture.done \
          sechiba_io_p.done \
          time.done \
          xios_orchidee.done
	touch $(FCM_DONEDIR)/$@

explicitsnow.o: \
          $(PPSRCDIR0__sechiba)/explicitsnow.f90 \
          FFLAGS__sechiba__explicitsnow.flags \
          constantes.o \
          constantes_soil.o \
          ioipsl_para.o \
          pft_parameters.o \
          qsat_moisture.o \
          sechiba_io_p.o \
          time.o \
          xios_orchidee.o
	fcm_internal compile:F sechiba $< $@

FFLAGS__sechiba__intersurf.flags: \
          FFLAGS__sechiba.flags
	touch $(FCM_FLAGSDIR)/$@

intersurf.done: \
          intersurf.o \
          constantes.done \
          constantes_soil.done \
          control.done \
          grid.done \
          ioipsl_para.done \
          ioipslctrl.done \
          mod_orchidee_para.done \
          pft_parameters.done \
          sechiba.done \
          solar.done \
          thermosoilc.done \
          time.done \
          timer.done \
          xios_orchidee.done
	touch $(FCM_DONEDIR)/$@

intersurf.o: \
          $(PPSRCDIR0__sechiba)/intersurf.f90 \
          FFLAGS__sechiba__intersurf.flags \
          constantes.o \
          constantes_soil.o \
          control.o \
          grid.o \
          ioipsl_para.o \
          ioipslctrl.o \
          mod_orchidee_para.o \
          pft_parameters.o \
          sechiba.o \
          solar.o \
          thermosoilc.o \
          time.o \
          timer.o \
          xios_orchidee.o
	fcm_internal compile:F sechiba $< $@

FFLAGS__sechiba__sechiba_io_p.flags: \
          FFLAGS__sechiba.flags
	touch $(FCM_FLAGSDIR)/$@

sechiba_io_p.done: \
          sechiba_io_p.o \
          constantes.done \
          ioipsl_para.done \
          mod_orchidee_para.done
	touch $(FCM_DONEDIR)/$@

sechiba_io_p.o: \
          $(PPSRCDIR0__sechiba)/sechiba_io_p.f90 \
          FFLAGS__sechiba__sechiba_io_p.flags \
          constantes.o \
          ioipsl_para.o \
          mod_orchidee_para.o
	fcm_internal compile:F sechiba $< $@

FFLAGS__sechiba__chemistry.flags: \
          FFLAGS__sechiba.flags
	touch $(FCM_FLAGSDIR)/$@

chemistry.done: \
          chemistry.o \
          constantes.done \
          grid.done \
          interpweight.done \
          ioipsl_para.done \
          pft_parameters.done \
          qsat_moisture.done \
          sechiba_io_p.done \
          time.done \
          xios_orchidee.done
	touch $(FCM_DONEDIR)/$@

chemistry.o: \
          $(PPSRCDIR0__sechiba)/chemistry.f90 \
          FFLAGS__sechiba__chemistry.flags \
          constantes.o \
          grid.o \
          interpweight.o \
          ioipsl_para.o \
          pft_parameters.o \
          qsat_moisture.o \
          sechiba_io_p.o \
          time.o \
          xios_orchidee.o
	fcm_internal compile:F sechiba $< $@

FFLAGS__sechiba__enerbil.flags: \
          FFLAGS__sechiba.flags
	touch $(FCM_FLAGSDIR)/$@

enerbil.done: \
          enerbil.o \
          constantes.done \
          constantes_soil.done \
          explicitsnow.done \
          ioipsl_para.done \
          pft_parameters.done \
          qsat_moisture.done \
          sechiba_io_p.done \
          time.done \
          xios_orchidee.done
	touch $(FCM_DONEDIR)/$@

enerbil.o: \
          $(PPSRCDIR0__sechiba)/enerbil.f90 \
          FFLAGS__sechiba__enerbil.flags \
          constantes.o \
          constantes_soil.o \
          explicitsnow.o \
          ioipsl_para.o \
          pft_parameters.o \
          qsat_moisture.o \
          sechiba_io_p.o \
          time.o \
          xios_orchidee.o
	fcm_internal compile:F sechiba $< $@

FFLAGS__sechiba__condveg.flags: \
          FFLAGS__sechiba.flags
	touch $(FCM_FLAGSDIR)/$@

condveg.done: \
          condveg.o \
          constantes.done \
          constantes_soil.done \
          grid.done \
          interpol_help.done \
          interpweight.done \
          ioipsl_para.done \
          mod_orchidee_para.done \
          pft_parameters.done \
          qsat_moisture.done \
          sechiba_io_p.done \
          xios_orchidee.done
	touch $(FCM_DONEDIR)/$@

condveg.o: \
          $(PPSRCDIR0__sechiba)/condveg.f90 \
          FFLAGS__sechiba__condveg.flags \
          constantes.o \
          constantes_soil.o \
          grid.o \
          interpol_help.o \
          interpweight.o \
          ioipsl_para.o \
          mod_orchidee_para.o \
          pft_parameters.o \
          qsat_moisture.o \
          sechiba_io_p.o \
          xios_orchidee.o
	fcm_internal compile:F sechiba $< $@

FFLAGS__sechiba__ioipslctrl.flags: \
          FFLAGS__sechiba.flags
	touch $(FCM_FLAGSDIR)/$@

ioipslctrl.done: \
          ioipslctrl.o \
          constantes.done \
          constantes_soil.done \
          grid.done \
          ioipsl_para.done \
          mod_orchidee_para.done \
          pft_parameters.done \
          thermosoilc.done \
          time.done
	touch $(FCM_DONEDIR)/$@

ioipslctrl.o: \
          $(PPSRCDIR0__sechiba)/ioipslctrl.f90 \
          FFLAGS__sechiba__ioipslctrl.flags \
          constantes.o \
          constantes_soil.o \
          grid.o \
          ioipsl_para.o \
          mod_orchidee_para.o \
          pft_parameters.o \
          thermosoilc.o \
          time.o
	fcm_internal compile:F sechiba $< $@

FFLAGS__sechiba__qsat_moisture.flags: \
          FFLAGS__sechiba.flags
	touch $(FCM_FLAGSDIR)/$@

qsat_moisture.done: \
          qsat_moisture.o \
          constantes.done \
          constantes_soil.done
	touch $(FCM_DONEDIR)/$@

qsat_moisture.o: \
          $(PPSRCDIR0__sechiba)/qsat_moisture.f90 \
          FFLAGS__sechiba__qsat_moisture.flags \
          constantes.o \
          constantes_soil.o
	fcm_internal compile:F sechiba $< $@

FFLAGS__sechiba__diffuco.flags: \
          FFLAGS__sechiba.flags
	touch $(FCM_FLAGSDIR)/$@

diffuco.done: \
          diffuco.o \
          chemistry.done \
          constantes.done \
          constantes_soil.done \
          grid.done \
          ioipsl_para.done \
          pft_parameters.done \
          qsat_moisture.done \
          sechiba_io_p.done \
          time.done \
          xios_orchidee.done
	touch $(FCM_DONEDIR)/$@

diffuco.o: \
          $(PPSRCDIR0__sechiba)/diffuco.f90 \
          FFLAGS__sechiba__diffuco.flags \
          chemistry.o \
          constantes.o \
          constantes_soil.o \
          grid.o \
          ioipsl_para.o \
          pft_parameters.o \
          qsat_moisture.o \
          sechiba_io_p.o \
          time.o \
          xios_orchidee.o
	fcm_internal compile:F sechiba $< $@

FFLAGS__sechiba__thermosoilc.flags: \
          FFLAGS__sechiba.flags
	touch $(FCM_FLAGSDIR)/$@

thermosoilc.done: \
          thermosoilc.o \
          constantes.done \
          constantes_soil.done \
          grid.done \
          interpweight.done \
          ioipsl_para.done \
          sechiba_io_p.done \
          time.done \
          xios_orchidee.done
	touch $(FCM_DONEDIR)/$@

thermosoilc.o: \
          $(PPSRCDIR0__sechiba)/thermosoilc.f90 \
          FFLAGS__sechiba__thermosoilc.flags \
          constantes.o \
          constantes_soil.o \
          grid.o \
          interpweight.o \
          ioipsl_para.o \
          sechiba_io_p.o \
          time.o \
          xios_orchidee.o
	fcm_internal compile:F sechiba $< $@

FFLAGS__sechiba__hydrol.flags: \
          FFLAGS__sechiba.flags
	touch $(FCM_FLAGSDIR)/$@

hydrol.done: \
          hydrol.o \
          constantes.done \
          constantes_soil.done \
          explicitsnow.done \
          grid.done \
          pft_parameters.done \
          sechiba_io_p.done \
          time.done \
          xios_orchidee.done
	touch $(FCM_DONEDIR)/$@

hydrol.o: \
          $(PPSRCDIR0__sechiba)/hydrol.f90 \
          FFLAGS__sechiba__hydrol.flags \
          constantes.o \
          constantes_soil.o \
          explicitsnow.o \
          grid.o \
          pft_parameters.o \
          sechiba_io_p.o \
          time.o \
          xios_orchidee.o
	fcm_internal compile:F sechiba $< $@

FFLAGS__sechiba__sechiba.flags: \
          FFLAGS__sechiba.flags
	touch $(FCM_FLAGSDIR)/$@

sechiba.done: \
          sechiba.o \
          chemistry.done \
          condveg.done \
          constantes.done \
          constantes_soil.done \
          diffuco.done \
          enerbil.done \
          grid.done \
          hydrol.done \
          hydrolc.done \
          ioipsl_para.done \
          pft_parameters.done \
          routing.done \
          sechiba_io_p.done \
          slowproc.done \
          thermosoil.done \
          thermosoilc.done \
          time.done \
          xios_orchidee.done
	touch $(FCM_DONEDIR)/$@

sechiba.o: \
          $(PPSRCDIR0__sechiba)/sechiba.f90 \
          FFLAGS__sechiba__sechiba.flags \
          chemistry.o \
          condveg.o \
          constantes.o \
          constantes_soil.o \
          diffuco.o \
          enerbil.o \
          grid.o \
          hydrol.o \
          hydrolc.o \
          ioipsl_para.o \
          pft_parameters.o \
          routing.o \
          sechiba_io_p.o \
          slowproc.o \
          thermosoil.o \
          thermosoilc.o \
          time.o \
          xios_orchidee.o
	fcm_internal compile:F sechiba $< $@

FFLAGS__sechiba__thermosoil.flags: \
          FFLAGS__sechiba.flags
	touch $(FCM_FLAGSDIR)/$@

thermosoil.done: \
          thermosoil.o \
          constantes.done \
          constantes_soil.done \
          grid.done \
          interpweight.done \
          ioipsl_para.done \
          sechiba_io_p.done \
          time.done \
          xios_orchidee.done
	touch $(FCM_DONEDIR)/$@

thermosoil.o: \
          $(PPSRCDIR0__sechiba)/thermosoil.f90 \
          FFLAGS__sechiba__thermosoil.flags \
          constantes.o \
          constantes_soil.o \
          grid.o \
          interpweight.o \
          ioipsl_para.o \
          sechiba_io_p.o \
          time.o \
          xios_orchidee.o
	fcm_internal compile:F sechiba $< $@

