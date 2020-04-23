# Makefile

# Defaults
VPATH		:=  src
BINDIR		:=  lib/Bin
OBJECTDIR	:=  lib/Objects
LIBDIR		:=  lib
COMPILER	:=  ifort
FFLAGS		:=  -Ofast -qopenmp -free -check bounds -mkl -liomp5 -lpthread -ldl -traceback -module $(OBJECTDIR)
FFLAGSDebug 	:=  -O0 -fpp -check bounds -traceback -warn nounused -module $(OBJECTDIR)

OBJECTS := $(BINDIR)/Debugger.o                 \
	$(BINDIR)/utilities.o                   \
	$(BINDIR)/quicksort.o                   \
	$(BINDIR)/SparseKit.o                   \
                                                \
	$(BINDIR)/Point.o                       \
	$(BINDIR)/PointPtr.o                    \
                                                \
	$(BINDIR)/Integrator.o                  \
	$(BINDIR)/IntegratorPtr.o               \
                                                \
	$(BINDIR)/Material.o                    \
	$(BINDIR)/MaterialPtr.o                 \
                                                \
	$(BINDIR)/Geometry.o                    \
	$(BINDIR)/GeometryPtr.o                 \
                                                \
	$(BINDIR)/Element.o                     \
                                                \
	$(BINDIR)/Element1D.o                   \
	$(BINDIR)/Element1DPtr.o		\
						\
	$(BINDIR)/Element2D.o                   \
	$(BINDIR)/Element2DPtr.o                \
	$(BINDIR)/TriangElement.o               \
	$(BINDIR)/QuadElement.o                 \
                                                \
	$(BINDIR)/functionOnPoints.o            \
	$(BINDIR)/functionOnLines.o             \
	$(BINDIR)/functionOnSurfaces.o          \
	$(BINDIR)/PointSource.o                 \
	$(BINDIR)/LineSource.o                  \
	$(BINDIR)/SurfaceSource.o               \
	$(BINDIR)/Source.o                      \
                                                \
	$(BINDIR)/Domain.o                      \
                                                \
	$(BINDIR)/Problem.o                     \
	                                        \
	$(BINDIR)/IOData.o                      \
                                                \
	$(BINDIR)/Solver.o  

#=============================================================
main: $(OBJECTS)
	ar rcv $(LIBDIR)/project789.a $(OBJECTS)

clean:
	rm -f $(BINDIR)/*.o $(LIBDIR)/*.a $(OBJECTDIR)/*.mod
#=============================================================

%*.a : $(BINDIR)/%.o
	ar rcv $^ $@

$(BINDIR)/main.o : main.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@ 

$(BINDIR)/Debugger.o : $(VPATH)/lib/Debugger.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/utilities.o : $(VPATH)/lib/Utilities.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/quicksort.o : $(VPATH)/lib/Quicksort.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/SparseKit.o : $(VPATH)/lib/SparseKit.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/Point.o : $(VPATH)/Point/Point.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/PointPtr.o : $(VPATH)/Point/PointPtr.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/Integrator.o : $(VPATH)/Integrator/Integrator.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/IntegratorPtr.o : $(VPATH)/Integrator/IntegratorPtr.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/Material.o : $(VPATH)/Material/Material.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/MaterialPtr.o : $(VPATH)/Material/MaterialPtr.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/Geometry.o : $(VPATH)/Geometry/Geometry.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/GeometryPtr.o : $(VPATH)/Geometry/GeometryPtr.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/Element.o : $(VPATH)/Element/Element.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/Element1D.o : $(VPATH)/Element/1D/Element1D.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Element1DPtr.o : $(VPATH)/Element/1D/Element1DPtr.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/Element2D.o : $(VPATH)/Element/2D/Element2D.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Element2DPtr.o : $(VPATH)/Element/2D/Element2DPtr.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/TriangElement.o : $(VPATH)/Element/2D/TriangElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/QuadElement.o : $(VPATH)/Element/2D/QuadElement.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/functionOnPoints.o : $(VPATH)/Source/functionOnPoints.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/functionOnLines.o : $(VPATH)/Source/functionOnLines.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/functionOnSurfaces.o : $(VPATH)/Source/functionOnSurfaces.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/PointSource.o : $(VPATH)/Source/PointSource.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/LineSource.o : $(VPATH)/Source/LineSource.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/SurfaceSource.o : $(VPATH)/Source/SurfaceSource.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Source.o : $(VPATH)/Source/Source.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/Domain.o : $(VPATH)/Domain/Domain.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/Problem.o : $(VPATH)/Problem/Problem.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/Solver.o : $(VPATH)/Solver/Solver.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/IOData.o : $(VPATH)/DataIO/IOData.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
#------------------------------------------------------------------------------
help:
	@echo ""
	@echo " make: compila la librería"
	@echo ""
	@echo ""
	@echo " clean: borra los *.o y *.mod"
	@echo ""
#------------------------------------------------------------------------------
