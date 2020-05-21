# Makefile

# Defaults
VPATH		:=  src
BINDIR		:=  lib/Bin
OBJECTDIR	:=  lib/Objects
LIBDIR		:=  lib
COMPILER	:=  ifort
FFLAGS		:=  -Ofast -qopenmp -free -check bounds -mkl -liomp5 -lpthread -ldl -traceback -module $(OBJECTDIR)
FFLAGSDebug 	:=  -O0 -fpp -check bounds -traceback -warn nounused -module $(OBJECTDIR)

OBJECTS := $(BINDIR)/Utilities.o                \
	$(BINDIR)/Debugger.o                    \
	$(BINDIR)/Quicksort.o                   \
	$(BINDIR)/SparseKit.o                   \
	                                        \
	$(BINDIR)/FortranParser.o               \
	$(BINDIR)/FortranParserPtr.o            \
	$(BINDIR)/Property.o                    \
	                                        \
	$(BINDIR)/Source.o                      \
	$(BINDIR)/SourcePtr.o                   \
                                                \
	$(BINDIR)/LeftHandSide.o                \
                                                \
	$(BINDIR)/Point.o                       \
	$(BINDIR)/Dof.o                         \
	$(BINDIR)/Node.o                        \
	$(BINDIR)/NodePtr.o                     \
                                                \
	$(BINDIR)/Integrator.o                  \
	$(BINDIR)/IntegratorPtr.o               \
	                                        \
	$(BINDIR)/Geometry.o                    \
	$(BINDIR)/Line1D2Node.o                 \
	$(BINDIR)/Line1D3Node.o                 \
	$(BINDIR)/Line2D2Node.o                 \
	$(BINDIR)/Line2D3Node.o                 \
	$(BINDIR)/Line3D2Node.o                 \
	$(BINDIR)/Line3D3Node.o                 \
	$(BINDIR)/Triangle2D3Node.o             \
	$(BINDIR)/Triangle2D6Node.o             \
	$(BINDIR)/Triangle3D3Node.o             \
	$(BINDIR)/Triangle3D6Node.o             \
	$(BINDIR)/Quadrilateral2D4Node.o        \
	$(BINDIR)/Quadrilateral2D8Node.o        \
	$(BINDIR)/Quadrilateral3D4Node.o        \
	$(BINDIR)/Quadrilateral3D8Node.o        \
	$(BINDIR)/Tetrahedron3D4Node.o          \
	$(BINDIR)/Tetrahedron3D10Node.o         \
	$(BINDIR)/Hexahedron3D8Node.o           \
	$(BINDIR)/Hexahedron3D20Node.o          \
	                                        \
	$(BINDIR)/ProcessInfo.o			\
                                                \
	$(BINDIR)/Element.o                     \
	$(BINDIR)/ElementPtr.o                  \
	                                        \
	$(BINDIR)/Condition.o                   \
	$(BINDIR)/ConditionPtr.o                \
                                                \
	$(BINDIR)/Mesh.o                        \
	$(BINDIR)/Model.o			\
						\
	$(BINDIR)/Preconditioner.o		\
	$(BINDIR)/PreconditionerMethod.o	\
	$(BINDIR)/UsePreconditioner.o		\
						\
	$(BINDIR)/ReorderSystem.o		\
	$(BINDIR)/ReorderSystemMethod.o		\
	$(BINDIR)/UseReorderSystem.o		\
                                                \
	$(BINDIR)/DirectLinearSolver.o		\
	$(BINDIR)/IterativeLinearSolver.o	\
	$(BINDIR)/LinearSolver.o		\
	$(BINDIR)/mklPardiso.o                  \
	$(BINDIR)/IterativeLinearSolverMethod.o \
						\
	$(BINDIR)/NonLinearSolvers.o		\
	$(BINDIR)/NonLinearSolver.o		\
	$(BINDIR)/Method.o			\
						\
	$(BINDIR)/Process.o			\
						\
	$(BINDIR)/BuilderAndSolver.o		\
	$(BINDIR)/Scheme.o			\
	$(BINDIR)/Integrand.o			\
	$(BINDIR)/BaseModel.o			\
	$(BINDIR)/RK2.o				\
	$(BINDIR)/RK4.o				\
	$(BINDIR)/AdamsB4.o			\
	$(BINDIR)/ExplicitEuler.o		\
	$(BINDIR)/SolvingStrategy.o             \
                                                \
	$(BINDIR)/GIDDataOutput.o               

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
$(BINDIR)/Utilities.o : $(VPATH)/lib/Utilities.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Quicksort.o : $(VPATH)/lib/Quicksort.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/SparseKit.o : $(VPATH)/lib/SparseKit.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/Point.o : $(VPATH)/node/Point.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Dof.o : $(VPATH)/node/Dof.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Node.o : $(VPATH)/node/Node.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/NodePtr.o : $(VPATH)/node/NodePtr.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/Integrator.o : $(VPATH)/integrator/Integrator.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/IntegratorPtr.o : $(VPATH)/integrator/IntegratorPtr.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/Geometry.o : $(VPATH)/geometry/Geometry.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Line1D2Node.o : $(VPATH)/geometry/Line1D2Node.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Line1D3Node.o : $(VPATH)/geometry/Line1D3Node.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Line2D2Node.o : $(VPATH)/geometry/Line2D2Node.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Line2D3Node.o : $(VPATH)/geometry/Line2D3Node.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Triangle2D3Node.o : $(VPATH)/geometry/Triangle2D3Node.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Triangle2D6Node.o : $(VPATH)/geometry/Triangle2D6Node.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Quadrilateral2D4Node.o : $(VPATH)/geometry/Quadrilateral2D4Node.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Quadrilateral2D8Node.o : $(VPATH)/geometry/Quadrilateral2D8Node.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Line3D2Node.o : $(VPATH)/geometry/Line3D2Node.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Line3D3Node.o : $(VPATH)/geometry/Line3D3Node.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Triangle3D3Node.o : $(VPATH)/geometry/Triangle3D3Node.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Triangle3D6Node.o : $(VPATH)/geometry/Triangle3D6Node.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Quadrilateral3D4Node.o : $(VPATH)/geometry/Quadrilateral3D4Node.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Quadrilateral3D8Node.o : $(VPATH)/geometry/Quadrilateral3D8Node.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Tetrahedron3D4Node.o : $(VPATH)/geometry/Tetrahedron3D4Node.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Tetrahedron3D10Node.o : $(VPATH)/geometry/Tetrahedron3D10Node.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Hexahedron3D8Node.o : $(VPATH)/geometry/Hexahedron3D8Node.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Hexahedron3D20Node.o : $(VPATH)/geometry/Hexahedron3D20Node.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/FortranParser.o : $(VPATH)/property/FortranParser.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/FortranParserPtr.o : $(VPATH)/property/FortranParserPtr.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Property.o : $(VPATH)/property/Property.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/Element.o : $(VPATH)/element/Element.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ElementPtr.o : $(VPATH)/element/ElementPtr.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/Condition.o : $(VPATH)/condition/Condition.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ConditionPtr.o : $(VPATH)/condition/ConditionPtr.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/Source.o : $(VPATH)/sources/Source.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/SourcePtr.o : $(VPATH)/sources/SourcePtr.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/LeftHandSide.o : $(VPATH)/leftHandSide/LeftHandSide.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/Mesh.o : $(VPATH)/model/Mesh.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ProcessInfo.o : $(VPATH)/model/ProcessInfo.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Model.o : $(VPATH)/model/Model.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/mklPardiso.o : $(VPATH)/solvers/Linear/Direct/Solvers/mklPardiso.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/DirectLinearSolver.o : $(VPATH)/solvers/Linear/Direct/DirectLinearSolver.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Preconditioner.o : $(VPATH)/solvers/Linear/Preconditioner/Preconditioner.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/PreconditionerMethod.o : $(VPATH)/solvers/Linear/Preconditioner/PreconditionerMethod.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/UsePreconditioner.o : $(VPATH)/solvers/Linear/Preconditioner/UsePreconditioner.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/IterativeLinearSolverMethod.o : $(VPATH)/solvers/Linear/Iterative/Solvers/IterativeLinearSolverMethod.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/IterativeLinearSolver.o : $(VPATH)/solvers/Linear/Iterative/IterativeLinearSolver.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ReorderSystem.o : $(VPATH)/solvers/Linear/Reorder/ReorderSystem.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ReorderSystemMethod.o : $(VPATH)/solvers/Linear/Reorder/ReorderSystemMethod.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/UseReorderSystem.o : $(VPATH)/solvers/Linear/Reorder/UseReorderSystem.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/LinearSolver.o : $(VPATH)/solvers/Linear/LinearSolver.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/NonLinearSolvers.o : $(VPATH)/solvers/NonLinear/NonLinearSolvers.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Method.o : $(VPATH)/solvers/NonLinear/Method.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/NonLinearSolver.o : $(VPATH)/solvers/NonLinear/NonLinearSolver.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/Process.o : $(VPATH)/Process/Process.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/BuilderAndSolver.o : $(VPATH)/SolvingStrategy/BuilderAndSolver.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Scheme.o : $(VPATH)/SolvingStrategy/Scheme.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/SolvingStrategy.o : $(VPATH)/SolvingStrategy/SolvingStrategy.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/BaseModel.o : $(VPATH)/IntegrationSchemes/BaseModel.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/Integrand.o : $(VPATH)/IntegrationSchemes/Integrand.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/ExplicitEuler.o : $(VPATH)/IntegrationSchemes/ExplicitEuler.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/RK2.o : $(VPATH)/IntegrationSchemes/RK2.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/RK4.o : $(VPATH)/IntegrationSchemes/RK4.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@
$(BINDIR)/AdamsB4.o : $(VPATH)/IntegrationSchemes/AdamsB4.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@

$(BINDIR)/GIDDataOutput.o : $(VPATH)/io/GIDDataOutput.f90
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
