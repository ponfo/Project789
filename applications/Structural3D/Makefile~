# Makefile

VPATH		:=  src
BINDIR		:=  Bin
OBJECTDIR	:=  Objects
LIB		:=  ../../lib
LIBOBJECTS	:=  $(LIB)/Objects
COMPILER	:=  ifort
FFLAGS		:=  -Ofast -qopenmp -free -check bounds -mkl -liomp5 -lpthread -ldl -traceback -module $(OBJECTDIR)
FFLAGSDebug 	:=  -O0 -fpp -check bounds -traceback -warn nounused -module $(OBJECTDIR)

OBJECTS := $(BINDIR)/StructuralMaterial.o      \
	$(BINDIR)/StructuralElement.o          \
	$(BINDIR)/Pressure.o                   \
	$(BINDIR)/NormalStress.o               \
	$(BINDIR)/ShearStress.o                \
	$(BINDIR)/Strain.o                     \
	$(BINDIR)/StructuralModel.o            \
	$(BINDIR)/StructuralBuilderAndSolver.o \
	$(BINDIR)/StructuralScheme.o           \
	$(BINDIR)/StructuralStrategy.o         \
	$(BINDIR)/Structural2DApplication.o    \
	$(BINDIR)/DataInput.o                  \
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

$(BINDIR)/StructuralMaterial.o : $(VPATH)/StructuralMaterial.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/StructuralElement.o : $(VPATH)/CustomElements/StructuralElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/Pressure.o : $(VPATH)/CustomConditions/Pressure.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/NormalStress.o : $(VPATH)/CustomResult/NormalStress.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/ShearStress.o : $(VPATH)/CustomResult/ShearStress.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/Strain.o : $(VPATH)/CustomResult/Strain.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/StructuralModel.o : $(VPATH)/Model/StructuralModel.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/StructuralBuilderAndSolver.o : $(VPATH)/SolvingStrategy/StructuralBuilderAndSolver.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/StructuralScheme.o : $(VPATH)/SolvingStrategy/StructuralScheme.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/StructuralStrategy.o : $(VPATH)/SolvingStrategy/StructuralStrategy.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/Structural2DApplication.o : $(VPATH)/Structural2DApplication.f90
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

