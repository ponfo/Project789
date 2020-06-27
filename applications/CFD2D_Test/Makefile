# Makefile

VPATH		:=  src
BINDIR		:=  Bin
OBJECTDIR	:=  Objects
LIB		:=  ../../lib
LIBOBJECTS	:=  $(LIB)/Objects
COMPILER	:=  ifort
FFLAGS		:=  -Ofast -qopenmp -free -mkl -liomp5 -lpthread -ldl -module $(OBJECTDIR)
FFLAGSDebug 	:=  -O0 -fpp -check bounds -traceback -warn nounused -module $(OBJECTDIR)

OBJECTS := $(BINDIR)/CFDMaterial.o      \
	$(BINDIR)/CFDElement.o      	\
	$(BINDIR)/Print.o               \
	$(BINDIR)/Results.o             \
	$(BINDIR)/NormalVelocity.o      \
	$(BINDIR)/CFDModel.o            \
	$(BINDIR)/CFDApplication.o      \
	$(BINDIR)/CFDBuilderAndSolver.o \
	$(BINDIR)/CFDScheme.o           \
	$(BINDIR)/NavierStokes2D.o      \
	$(BINDIR)/DataInput.o           \
	$(BINDIR)/SolvingStrategy.o     \
	$(BINDIR)/CFDStrategy.o         \
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
$(BINDIR)/CFDElement.o : $(VPATH)/CustomElements/CFDElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/CFDModel.o : $(VPATH)/Model/CFDModel.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/CFDBuilderAndSolver.o : $(VPATH)/SolvingStrategy/CFDBuilderAndSolver.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/CFDScheme.o : $(VPATH)/SolvingStrategy/CFDScheme.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/CFDStrategy.o : $(VPATH)/SolvingStrategy/CFDStrategy.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/SolvingStrategy.o : $(VPATH)/SolvingStrategy/SolvingStrategy.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/CFDApplication.o : $(VPATH)/CFDApplication.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/DataInput.o : $(VPATH)/CustomIO/DataInput.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/CFDMaterial.o : $(VPATH)/CFDMaterial.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/Print.o : $(VPATH)/Process/Print.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/Results.o : $(VPATH)/CustomResults/Results.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/NormalVelocity.o : $(VPATH)/CustomConditions/NormalVelocity.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/NavierStokes2D.o : $(VPATH)/SolvingStrategy/NavierStokes2D.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
#========================================================================
help:
	@echo ""
	@echo " make      : compila la aplicacion"
	@echo ""
	@echo " make clean: borra los *.o y *.mod"
	@echo ""
#========================================================================

