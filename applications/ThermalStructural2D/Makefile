# Makefile

VPATH		:=  src
BINDIR		:=  Bin
OBJECTDIR	:=  Objects
LIB		:=  ../../lib
THERMALSRC      :=  ../Thermal2D/src
STRUCTURALSRC   :=  ../Structural2D/src
LIBOBJECTS	:=  $(LIB)/Objects
COMPILER	:=  ifort
FFLAGS		:=  -Ofast -qopenmp -free -check bounds -mkl -liomp5 -lpthread -ldl -traceback -module $(OBJECTDIR)
FFLAGSDebug 	:=  -O0 -fpp -check bounds -traceback -warn nounused -module $(OBJECTDIR)

OBJECTS := $(BINDIR)/ThermalMaterial.o          \
	$(BINDIR)/StructuralMaterial.o          \
	$(BINDIR)/ThermalElement.o              \
	$(BINDIR)/ConvectionOnLine.o            \
	$(BINDIR)/FluxOnLine.o                  \
	$(BINDIR)/HeatFlux.o                    \
	$(BINDIR)/ThermalModel.o                \
	$(BINDIR)/ThermalBuilderAndSolver.o     \
	$(BINDIR)/ThermalScheme.o               \
	$(BINDIR)/ThermalStrategy.o             \
	                                        \
	$(BINDIR)/Pressure.o                    \
	$(BINDIR)/NormalStress.o                \
	$(BINDIR)/ShearStress.o                 \
	$(BINDIR)/Strain.o                      \
	$(BINDIR)/StructuralModel.o             \
	$(BINDIR)/StructuralBuilderAndSolver.o  \
	$(BINDIR)/StructuralScheme.o            \
	$(BINDIR)/StructuralStrategy.o          \
                                                \
	$(BINDIR)/ThermalStructuralElement.o    \
	$(BINDIR)/ThermalStruct2DApplication.o  \
	$(BINDIR)/DataInput.o                   \
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

$(BINDIR)/ThermalMaterial.o : $(THERMALSRC)/ThermalMaterial.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/ThermalElement.o : $(THERMALSRC)/CustomElements/ThermalElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/ConvectionOnLine.o : $(THERMALSRC)/CustomConditions/ConvectionOnLine.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/FluxOnLine.o : $(THERMALSRC)/CustomConditions/FluxOnLine.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/HeatFlux.o : $(THERMALSRC)/CustomResult/HeatFlux.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/ThermalModel.o : $(THERMALSRC)/Model/ThermalModel.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/ThermalBuilderAndSolver.o : $(THERMALSRC)/SolvingStrategy/ThermalBuilderAndSolver.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/ThermalScheme.o : $(THERMALSRC)/SolvingStrategy/ThermalScheme.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/ThermalStrategy.o : $(THERMALSRC)/SolvingStrategy/ThermalStrategy.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@


$(BINDIR)/Pressure.o : $(STRUCTURALSRC)/CustomConditions/Pressure.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/NormalStress.o : $(STRUCTURALSRC)/CustomResult/NormalStress.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/ShearStress.o : $(STRUCTURALSRC)/CustomResult/ShearStress.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/Strain.o : $(STRUCTURALSRC)/CustomResult/Strain.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/StructuralModel.o : $(STRUCTURALSRC)/Model/StructuralModel.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/StructuralBuilderAndSolver.o : $(STRUCTURALSRC)/SolvingStrategy/StructuralBuilderAndSolver.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/StructuralScheme.o : $(STRUCTURALSRC)/SolvingStrategy/StructuralScheme.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/StructuralStrategy.o : $(STRUCTURALSRC)/SolvingStrategy/StructuralStrategy.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@


$(BINDIR)/StructuralMaterial.o : $(VPATH)/ThermalStructuralMaterial.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/ThermalStructuralElement.o : $(VPATH)/CustomElements/ThermalStructuralElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/ThermalStruct2DApplication.o : $(VPATH)/ThermalStruct2DApplication.f90
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

