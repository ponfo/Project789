# Makefile
# Defaults
VPATH		:=  Source
BINDIR		:=  Bin
OBJECTDIR	:=  Objects
LIB		:=  ../../lib
LIBOBJECTS	:=  $(LIB)/Objects
COMPILER	:=  ifort
FFLAGS		:=  -Ofast -qopenmp -free -check bounds -mkl -liomp5 -lpthread -ldl -traceback -module $(OBJECTDIR)
FFLAGSDebug 	:=  -O0 -fpp -check bounds -traceback -warn nounused -module $(OBJECTDIR)

OBJECTS :=$(BINDIR)/				\
						\
	$(BINDIR)/ThermalModel.f90		\
						\
	$(BINDIR)/DirectBuilderAndSolver.f90	\
	$(BINDIR)/DirectScheme.f90		\
	$(BINDIR)/DirectStrategy.f90		\
						\
	$(BINDIR)/Thermal2DGid.o                \
						\
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

$(BINDIR)/Thermal2DGid.o : $(VPATH)/Thermal2DGid.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@

$(BINDIR)/DirectStrategy.o : $(VPATH)/SolvingStrategy/DirectStrategy.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/DirectScheme.o : $(VPATH)/SolvingStrategy/DirectScheme.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/DirectBuilderAndSolver.o : $(VPATH)/SolvingStrategy/DirectBuilderAndSolver.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@

$(BINDIR)/Thermalmodel.o : $(VPATH)/Model/ThermalModel.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@

#------------------------------------------------------------------------------
help:
	@echo "                                                 "
	@echo " make              : compila la aplicacion       "
	@echo "                                                 "
	@echo " make clean        : borra los *.o, *.mod y main "
	@echo "                                                 "
#------------------------------------------------------------------------------
