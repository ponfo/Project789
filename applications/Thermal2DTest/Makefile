# Makefile

COMPILER	:=  ifort
VPATH		:=  $(PWD)/src
BINDIR		:=  $(PWD)/Bin
OBJECTDIR	:=  $(PWD)/Objects
FRAMEWORK       :=  ../lib/project789.a
FFLAGS		:=  -Ofast -free -static -check bounds -traceback -I $(FRAMEWORK) -module $(OBJECTDIR)
FFLAGSDebug 	:=  -O0 -fpp -check bounds -traceback -warn nounused -I $(FRAMEWORK) -module $(OBJECTDIR)

OBJECTS := $(BINDIR)/ThermalMaterial.o     \
	$(BINDIR)/ThermalElement.o         \
	$(BINDIR)/ConvectionOnLine.o       \
	$(BINDIR)/FluxOnLine.o             \
	$(BINDIR)/ThermalModel.o           \
	$(BINDIR)/DirectBuilderAndSolver.o \
	$(BINDIR)/DirectStrategy.o         \
	$(BINDIR)/Thermal2DApplication.o   \
	$(BINDIR)/DataInput.o              \
	$(BINDIR)/main.o

main: $(OBJECTS)
	$(COMPILER) $(FFLAGS) $(OBJECTS) -o main

debug: $(OBJECTS)
	$(COMPILER) $(FFLAGSDebug) $(OBJECTS) -o main

$(BINDIR)/main.o : main.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@ 

$(BINDIR)/ThermalMaterial.o : $(VPATH)/ThermalMaterial.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ThermalElement.o : $(VPATH)/CustomElement/ThermalElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ConvectionOnLine.o : $(VPATH)/CustomCondition/ConvectionOnLine.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/FluxOnLine.o : $(VPATH)/CustomCondition/FluxOnLine.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ThermalModel.o : $(VPATH)/Model/ThermalModel.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/DirectBuilderAndSolver.o : $(VPATH)/SolvingStrategy/DirectBuilderAndSolver.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/DirectStrategy.o : $(VPATH)/SolvingStrategy/DirectStrategy.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Thermal2DApplication.o : $(VPATH)/Thermal2DApplication.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/DataInput.o : $(VPATH)/CustomIO/DataInput.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@


clean:
	rm -f $(BINDIR)/*.o main $(OBJECTDIR)/*.mod

help:
	@echo ""
	@echo " main: compila la aplicacion"
	@echo ""
	@echo " debug: busca errores"
	@echo ""
	@echo " clean: borra los *.o y *.mod"
	@echo ""
#------------------------------------------------------------------------------

