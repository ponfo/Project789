# Makefile

VPATH		:=  src
BINDIR		:=  Bin
OBJECTDIR	:=  Objects
LIB		:=  ../../lib
LIBOBJECTS	:=  $(LIB)/Objects
COMPILER	:=  ifort
FFLAGS		:=  -Ofast -qopenmp -free -check bounds -mkl -liomp5 -lpthread -ldl -traceback -module $(OBJECTDIR)
FFLAGSDebug 	:=  -O0 -fpp -check bounds -traceback -warn nounused -module $(OBJECTDIR)

OBJECTS := $(BINDIR)/ThermalMaterial.o     \
	$(BINDIR)/ThermalElement.o         \
	$(BINDIR)/ConvectionOnLine.o       \
	$(BINDIR)/FluxOnLine.o             \
	$(BINDIR)/HeatFlux.o               \
	$(BINDIR)/ThermalModel.o           \
	$(BINDIR)/DirectBuilderAndSolver.o \
	$(BINDIR)/DirectScheme.o           \
	$(BINDIR)/DirectStrategy.o         \
	$(BINDIR)/Thermal2DApplication.o   \
	$(BINDIR)/DataInput.o              \
	$(BINDIR)/main.o

LIBRARIES := $(LIB)/project789.a

#========================================================================
main: $(OBJECTS)
	$(COMPILER) $(FFLAGS) $^ -I$(LIBOBJECTS) -L$(LIB) $(LIBRARIES) -o $@

clean:
	rm -f $(BINDIR)/*.o main $(OBJECTDIR)/*.mod
#========================================================================

$(BINDIR)/main.o : main.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@ 

$(BINDIR)/ThermalMaterial.o : $(VPATH)/ThermalMaterial.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/ThermalElement.o : $(VPATH)/CustomElements/ThermalElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/ConvectionOnLine.o : $(VPATH)/CustomConditions/ConvectionOnLine.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/FluxOnLine.o : $(VPATH)/CustomConditions/FluxOnLine.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/HeatFlux.o : $(VPATH)/CustomResult/HeatFlux.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/ThermalModel.o : $(VPATH)/Model/ThermalModel.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/DirectBuilderAndSolver.o : $(VPATH)/SolvingStrategy/DirectBuilderAndSolver.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/DirectScheme.o : $(VPATH)/SolvingStrategy/DirectScheme.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/DirectStrategy.o : $(VPATH)/SolvingStrategy/DirectStrategy.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/Thermal2DApplication.o : $(VPATH)/Thermal2DApplication.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/DataInput.o : $(VPATH)/CustomIO/DataInput.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@



help:
	@echo ""
	@echo " make      : compila la aplicacion"
	@echo ""
	@echo " make clean: borra los *.o y *.mod"
	@echo ""
#------------------------------------------------------------------------------

