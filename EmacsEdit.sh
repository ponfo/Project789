#!/bin/sh

emacs   src/lib/Utilities.f90    \
	src/lib/Debugger.f90     \
	src/lib/Quicksort.f90    \
	src/lib/SparseKit.f90    \
	\
	src/solvers/mklPardiso.f90 \
	\
	src/integrator/Integrator.f90  \
	src/integrator/IntegratorPtr.f90 \
	\
	src/property/FortranParser.f90  \
	src/property/Property.f90       \
	\
	src/node/Point.f90 \
	src/node/Dof.f90   \
	src/node/Node.f90  \
	src/node/NodePtr.f90 \
	\
	src/geometry/Geometry.f90 \
	src/geometry/Triangle2D3Node.f90 \
	src/geometry/Triangle2D6Node.f90 \
	src/geometry/Quadrilateral2D4Node.f90 \
	src/geometry/Quadrilateral2D8Node.f90 \
	\
	src/includes/GeometryObject.f90 \
	src/includes/Element.f90 \
	src/includes/Condition.f90\
	\
	src/model/Mesh.f90 \
	src/model/ProcessInfo.f90 \
	src/model/Model.f90       &
