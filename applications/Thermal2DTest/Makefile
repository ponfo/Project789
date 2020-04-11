# Makefile

help:
	@echo ""
	@echo " main: compila la aplicacion"
	@echo ""
	@echo " debug: busca errores"
	@echo ""
	@echo " clean: borra los *.o y *.mod"
	@echo ""
#------------------------------------------------------------------------------

# Defaults
COMPILER	:=  ifort
FFLAGS		:=  -Ofast -qopenmp -free -check bounds -mkl -liomp5 -lpthread -ldl -traceback -module $(OBJECTDIR)
FFLAGSDebug 	:=  -O0 -fpp -check bounds -traceback -warn nounused -module $(OBJECTDIR)
VPATH		:=  Source
BINDIR		:=  Bin
OBJECTDIR	:=  Objects
LIB		:=  ../../lib
LIBOBJECTS	:=  $(LIB)/Objects

OBJECTS := $(BINDIR)/ThermalMaterial.o          \
                                                \
	$(BINDIR)/ThermalElement1D.o            \
                                                \
	$(BINDIR)/ThermalElement2D.o            \
	$(BINDIR)/ThermalLinTriangElement.o     \
	$(BINDIR)/ThermalQuadTriangElement.o    \
	$(BINDIR)/ThermalLinQuadElement.o       \
	$(BINDIR)/ThermalQuadQuadElement.o      \
	$(BINDIR)/ThermalElementList2D.o        \
                                                \
	$(BINDIR)/DirichletPoint.o              \
	$(BINDIR)/NormalFluxPoint.o             \
	$(BINDIR)/ConvectionPoint.o             \
	$(BINDIR)/ThermalBoundaryCondition1D.o  \
	$(BINDIR)/NormalFluxLine.o              \
	$(BINDIR)/ConvectionLine.o              \
	$(BINDIR)/ThermalBoundaryCondition2D.o  \
                                                \
	$(BINDIR)/ThermalDomain.o               \
                                                \
	$(BINDIR)/ThermalProblem.o              \
	                                        \
	$(BINDIR)/HeatFlux.o                    \
	                                        \
	$(BINDIR)/DataInput.o                   \
	$(BINDIR)/DataOutput.o                  \
                                                \
	$(BINDIR)/Solver.o                      \
	$(BINDIR)/main.o

LIBRARIES := $(LIB)/project789.a

main: $(OBJECTS)
	$(COMPILER) $(FFLAGS) $(OBJECTS) -I$(LIBOBJECTS) -L$(LIB) $(LIBRARIES) -o main

debug: $(OBJECTS)
	$(COMPILER) $(FFLAGSDebug) $(OBJECTS) -o main

$(BINDIR)/main.o : main.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@ 

$(BINDIR)/ThermalMaterial.o : $(VPATH)/Material/ThermalMaterial.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@

$(BINDIR)/ThermalElement1D.o : $(VPATH)/Element/1D/ThermalElement1D.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/ThermalLinLineElement.o : $(VPATH)/Element/1D/ThermalLinLineElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/ThermalQuadLineElement.o : $(VPATH)/Element/1D/ThermalQuadLineElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/ThermalElementList1D.o : $(VPATH)/Element/1D/ThermalElementList1D.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@

$(BINDIR)/ThermalElement2D.o : $(VPATH)/Element/2D/ThermalElement2D.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/ThermalLinTriangElement.o : $(VPATH)/Element/2D/ThermalLinTriangElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/ThermalQuadTriangElement.o : $(VPATH)/Element/2D/ThermalQuadTriangElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/ThermalLinQuadElement.o : $(VPATH)/Element/2D/ThermalLinQuadElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/ThermalQuadQuadElement.o : $(VPATH)/Element/2D/ThermalQuadQuadElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/ThermalElementList2D.o : $(VPATH)/Element/2D/ThermalElementList2D.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@

$(BINDIR)/DirichletPoint.o : $(VPATH)/BoundaryCondition/DirichletPoint.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/NormalFluxPoint.o : $(VPATH)/BoundaryCondition/NormalFluxPoint.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/ConvectionPoint.o : $(VPATH)/BoundaryCondition/ConvectionPoint.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/ThermalBoundaryCondition1D.o : $(VPATH)/BoundaryCondition/ThermalBoundaryCondition1D.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/NormalFluxLine.o : $(VPATH)/BoundaryCondition/NormalFluxLine.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/ConvectionLine.o : $(VPATH)/BoundaryCondition/ConvectionLine.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/ThermalBoundaryCondition2D.o : $(VPATH)/BoundaryCondition/ThermalBoundaryCondition2D.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@

$(BINDIR)/ThermalDomain.o : $(VPATH)/Domain/ThermalDomain.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@

$(BINDIR)/ThermalProblem.o : $(VPATH)/Problem/ThermalProblem.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@

$(BINDIR)/Solver.o : $(VPATH)/Solver/Solver.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@

$(BINDIR)/HeatFlux.o : $(VPATH)/PostProcess/HeatFlux.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@

$(BINDIR)/DataInput.o : $(VPATH)/DataIO/DataInput.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@
$(BINDIR)/DataOutput.o : $(VPATH)/DataIO/DataOutput.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@

clean:
	rm -f $(BINDIR)/*.o main $(OBJECTDIR)/*.mod
