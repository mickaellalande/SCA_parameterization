# Automatic Make rule for orchidee_ol

PPSRCDIR0__orchidee_ol = /gpfswork/rech/goe/ufz23bm/SCA_parameterization/modipsl/modeles/ORCHIDEE/.config/ppsrc/orchidee_ol

SRCDIR0__orchidee_ol = /gpfswork/rech/goe/ufz23bm/SCA_parameterization/modipsl/modeles/ORCHIDEE/src_driver

orchidee_ol.etc : \
          $(SRCDIR0__orchidee_ol)/Makefile \
          $(SRCDIR0__orchidee_ol)/AA_make.ldef \
          $(SRCDIR0__orchidee_ol)/AA_make \
          $(FCM_DONEDIR)/FCM_CP.dummy
	cp $^ $(FCM_ETCDIR)
	touch $(FCM_DONEDIR)/$@

FFLAGS__orchidee_ol__getlandseamask.flags: \
          FFLAGS__orchidee_ol.flags
	touch $(FCM_FLAGSDIR)/$@

getlandseamask.done: \
          getlandseamask.o \
          constantes.done \
          constantes_soil_var.done \
          constantes_var.done \
          control.done \
          ioipsl_para.done \
          mod_orchidee_para.done
	touch $(FCM_DONEDIR)/$@

getlandseamask.o: \
          $(PPSRCDIR0__orchidee_ol)/getlandseamask.f90 \
          FFLAGS__orchidee_ol__getlandseamask.flags \
          constantes.o \
          constantes_soil_var.o \
          constantes_var.o \
          control.o \
          ioipsl_para.o \
          mod_orchidee_para.o
	fcm_internal compile:F orchidee_ol $< $@

FFLAGS__orchidee_ol__forcesoil.flags: \
          FFLAGS__orchidee_ol.flags
	touch $(FCM_FLAGSDIR)/$@

LDFLAGS__orchidee_ol__forcesoil.flags: \
          LDFLAGS__orchidee_ol.flags
	touch $(FCM_FLAGSDIR)/$@

LD__orchidee_ol__forcesoil.flags: \
          LD__orchidee_ol.flags
	touch $(FCM_FLAGSDIR)/$@

forcesoil.exe: \
          forcesoil.o \
          LD__orchidee_ol__forcesoil.flags \
          LDFLAGS__orchidee_ol__forcesoil.flags \
          $(OBJECTS__parallel) \
          $(OBJECTS__orchidee_ol) \
          $(OBJECTS__orchoasis) \
          $(OBJECTS__global) \
          $(OBJECTS__parameters) \
          $(OBJECTS__ext_src) \
          $(OBJECTS__stomate) \
          $(OBJECTS__sechiba) \
          constantes.done \
          constantes_mtc.done \
          ioipsl_para.done \
          mod_orchidee_para.done \
          pft_parameters.done \
          stomate_data.done \
          stomate_soilcarbon.done
	fcm_internal load orchidee_ol $< $@

forcesoil.o: \
          $(PPSRCDIR0__orchidee_ol)/forcesoil.f90 \
          FFLAGS__orchidee_ol__forcesoil.flags \
          constantes.o \
          constantes_mtc.o \
          ioipsl_para.o \
          mod_orchidee_para.o \
          pft_parameters.o \
          stomate_data.o \
          stomate_soilcarbon.o
	fcm_internal compile:F orchidee_ol $< $@

FFLAGS__orchidee_ol__teststomate.flags: \
          FFLAGS__orchidee_ol.flags
	touch $(FCM_FLAGSDIR)/$@

LDFLAGS__orchidee_ol__teststomate.flags: \
          LDFLAGS__orchidee_ol.flags
	touch $(FCM_FLAGSDIR)/$@

LD__orchidee_ol__teststomate.flags: \
          LD__orchidee_ol.flags
	touch $(FCM_FLAGSDIR)/$@

teststomate.exe: \
          teststomate.o \
          LD__orchidee_ol__teststomate.flags \
          LDFLAGS__orchidee_ol__teststomate.flags \
          $(OBJECTS__parallel) \
          $(OBJECTS__orchidee_ol) \
          $(OBJECTS__orchoasis) \
          $(OBJECTS__global) \
          $(OBJECTS__parameters) \
          $(OBJECTS__ext_src) \
          $(OBJECTS__stomate) \
          $(OBJECTS__sechiba) \
          constantes.done \
          constantes_soil.done \
          grid.done \
          intersurf.done \
          ioipsl_para.done \
          ioipslctrl.done \
          mod_orchidee_para.done \
          mod_orchidee_para_var.done \
          pft_parameters.done \
          slowproc.done \
          stomate.done \
          stomate_data.done \
          time.done
	fcm_internal load orchidee_ol $< $@

teststomate.o: \
          $(PPSRCDIR0__orchidee_ol)/teststomate.f90 \
          FFLAGS__orchidee_ol__teststomate.flags \
          constantes.o \
          constantes_soil.o \
          grid.o \
          intersurf.o \
          ioipsl_para.o \
          ioipslctrl.o \
          mod_orchidee_para.o \
          mod_orchidee_para_var.o \
          pft_parameters.o \
          slowproc.o \
          stomate.o \
          stomate_data.o \
          time.o
	fcm_internal compile:F orchidee_ol $< $@

FFLAGS__orchidee_ol__orchideedriver.flags: \
          FFLAGS__orchidee_ol.flags
	touch $(FCM_FLAGSDIR)/$@

LDFLAGS__orchidee_ol__orchideedriver.flags: \
          LDFLAGS__orchidee_ol.flags
	touch $(FCM_FLAGSDIR)/$@

LD__orchidee_ol__orchideedriver.flags: \
          LD__orchidee_ol.flags
	touch $(FCM_FLAGSDIR)/$@

orchidedriver.o: \
          $(PPSRCDIR0__orchidee_ol)/orchideedriver.f90 \
          FFLAGS__orchidee_ol__orchideedriver.flags \
          constantes.o \
          control.o \
          forcing_tools.o \
          globgrd.o \
          grid.o \
          ioipsl_para.o \
          ioipslctrl.o \
          mod_orchidee_para.o \
          sechiba.o \
          thermosoilc.o \
          time.o \
          timer.o
	fcm_internal compile:F orchidee_ol $< $@

orchideedriver.exe: \
          orchidedriver.o \
          LD__orchidee_ol__orchideedriver.flags \
          LDFLAGS__orchidee_ol__orchideedriver.flags \
          $(OBJECTS__parallel) \
          $(OBJECTS__orchidee_ol) \
          $(OBJECTS__orchoasis) \
          $(OBJECTS__global) \
          $(OBJECTS__parameters) \
          $(OBJECTS__ext_src) \
          $(OBJECTS__stomate) \
          $(OBJECTS__sechiba) \
          constantes.done \
          control.done \
          forcing_tools.done \
          globgrd.done \
          grid.done \
          ioipsl_para.done \
          ioipslctrl.done \
          mod_orchidee_para.done \
          sechiba.done \
          thermosoilc.done \
          time.done \
          timer.done
	fcm_internal load orchidee_ol $< $@

FFLAGS__orchidee_ol__readdim2.flags: \
          FFLAGS__orchidee_ol.flags
	touch $(FCM_FLAGSDIR)/$@

readdim2.done: \
          readdim2.o \
          timer.done \
          constantes.done \
          grid.done \
          ioipsl_para.done \
          mod_orchidee_para.done \
          solar.done \
          time.done \
          weather.done
	touch $(FCM_DONEDIR)/$@

readdim2.o: \
          $(PPSRCDIR0__orchidee_ol)/readdim2.f90 \
          FFLAGS__orchidee_ol__readdim2.flags \
          timer.o \
          constantes.o \
          grid.o \
          ioipsl_para.o \
          mod_orchidee_para.o \
          solar.o \
          time.o \
          weather.o
	fcm_internal compile:F orchidee_ol $< $@

FFLAGS__orchidee_ol__getprec.flags: \
          FFLAGS__orchidee_ol.flags
	touch $(FCM_FLAGSDIR)/$@

LDFLAGS__orchidee_ol__getprec.flags: \
          LDFLAGS__orchidee_ol.flags
	touch $(FCM_FLAGSDIR)/$@

LD__orchidee_ol__getprec.flags: \
          LD__orchidee_ol.flags
	touch $(FCM_FLAGSDIR)/$@

getprec.exe: \
          getprec.o \
          LD__orchidee_ol__getprec.flags \
          LDFLAGS__orchidee_ol__getprec.flags \
          $(OBJECTS__parallel) \
          $(OBJECTS__orchidee_ol) \
          $(OBJECTS__orchoasis) \
          $(OBJECTS__global) \
          $(OBJECTS__parameters) \
          $(OBJECTS__ext_src) \
          $(OBJECTS__stomate) \
          $(OBJECTS__sechiba)
	fcm_internal load orchidee_ol $< $@

getprec.o: \
          $(PPSRCDIR0__orchidee_ol)/getprec.f90 \
          FFLAGS__orchidee_ol__getprec.flags
	fcm_internal compile:F orchidee_ol $< $@

FFLAGS__orchidee_ol__globgrd.flags: \
          FFLAGS__orchidee_ol.flags
	touch $(FCM_FLAGSDIR)/$@

globgrd.done: \
          globgrd.o \
          forcing_tools.done \
          grid.done
	touch $(FCM_DONEDIR)/$@

globgrd.o: \
          $(PPSRCDIR0__orchidee_ol)/globgrd.f90 \
          FFLAGS__orchidee_ol__globgrd.flags \
          forcing_tools.o \
          grid.o
	fcm_internal compile:F orchidee_ol $< $@

FFLAGS__orchidee_ol__dim2_driver.flags: \
          FFLAGS__orchidee_ol.flags
	touch $(FCM_FLAGSDIR)/$@

LDFLAGS__orchidee_ol__dim2_driver.flags: \
          LDFLAGS__orchidee_ol.flags
	touch $(FCM_FLAGSDIR)/$@

LD__orchidee_ol__dim2_driver.flags: \
          LD__orchidee_ol.flags
	touch $(FCM_FLAGSDIR)/$@

dim2_driver.exe: \
          driver.o \
          LD__orchidee_ol__dim2_driver.flags \
          LDFLAGS__orchidee_ol__dim2_driver.flags \
          $(OBJECTS__parallel) \
          $(OBJECTS__orchidee_ol) \
          $(OBJECTS__orchoasis) \
          $(OBJECTS__global) \
          $(OBJECTS__parameters) \
          $(OBJECTS__ext_src) \
          $(OBJECTS__stomate) \
          $(OBJECTS__sechiba) \
          constantes.done \
          grid.done \
          intersurf.done \
          ioipsl_para.done \
          mod_orchidee_para.done \
          readdim2.done \
          time.done \
          timer.done
	fcm_internal load orchidee_ol $< $@

driver.o: \
          $(PPSRCDIR0__orchidee_ol)/dim2_driver.f90 \
          FFLAGS__orchidee_ol__dim2_driver.flags \
          constantes.o \
          grid.o \
          intersurf.o \
          ioipsl_para.o \
          mod_orchidee_para.o \
          readdim2.o \
          time.o \
          timer.o
	fcm_internal compile:F orchidee_ol $< $@

FFLAGS__orchidee_ol__forcing_tools.flags: \
          FFLAGS__orchidee_ol.flags
	touch $(FCM_FLAGSDIR)/$@

forcing_tools.done: \
          forcing_tools.o \
          constantes.done \
          mod_orchidee_para.done \
          solar.done \
          time.done
	touch $(FCM_DONEDIR)/$@

forcing_tools.o: \
          $(PPSRCDIR0__orchidee_ol)/forcing_tools.f90 \
          FFLAGS__orchidee_ol__forcing_tools.flags \
          constantes.o \
          mod_orchidee_para.o \
          solar.o \
          time.o
	fcm_internal compile:F orchidee_ol $< $@

FFLAGS__orchidee_ol__testrouting.flags: \
          FFLAGS__orchidee_ol.flags
	touch $(FCM_FLAGSDIR)/$@

LDFLAGS__orchidee_ol__testrouting.flags: \
          LDFLAGS__orchidee_ol.flags
	touch $(FCM_FLAGSDIR)/$@

LD__orchidee_ol__testrouting.flags: \
          LD__orchidee_ol.flags
	touch $(FCM_FLAGSDIR)/$@

testrouting.exe: \
          testrouting.o \
          LD__orchidee_ol__testrouting.flags \
          LDFLAGS__orchidee_ol__testrouting.flags \
          $(OBJECTS__parallel) \
          $(OBJECTS__orchidee_ol) \
          $(OBJECTS__orchoasis) \
          $(OBJECTS__global) \
          $(OBJECTS__parameters) \
          $(OBJECTS__ext_src) \
          $(OBJECTS__stomate) \
          $(OBJECTS__sechiba) \
          constantes.done \
          constantes_soil_var.done \
          constantes_var.done \
          control.done \
          getlandseamask.done \
          grid.done \
          ioipsl_para.done \
          mod_orchidee_para.done \
          pft_parameters.done \
          routing.done \
          time.done \
          timer.done
	fcm_internal load orchidee_ol $< $@

testrouting.o: \
          $(PPSRCDIR0__orchidee_ol)/testrouting.f90 \
          FFLAGS__orchidee_ol__testrouting.flags \
          constantes.o \
          constantes_soil_var.o \
          constantes_var.o \
          control.o \
          getlandseamask.o \
          grid.o \
          ioipsl_para.o \
          mod_orchidee_para.o \
          pft_parameters.o \
          routing.o \
          time.o \
          timer.o
	fcm_internal compile:F orchidee_ol $< $@

FFLAGS__orchidee_ol__weather.flags: \
          FFLAGS__orchidee_ol.flags
	touch $(FCM_FLAGSDIR)/$@

weather.done: \
          weather.o \
          constantes.done \
          ioipsl_para.done \
          mod_orchidee_para.done \
          solar.done \
          time.done
	touch $(FCM_DONEDIR)/$@

weather.o: \
          $(PPSRCDIR0__orchidee_ol)/weather.f90 \
          FFLAGS__orchidee_ol__weather.flags \
          constantes.o \
          ioipsl_para.o \
          mod_orchidee_para.o \
          solar.o \
          time.o
	fcm_internal compile:F orchidee_ol $< $@

