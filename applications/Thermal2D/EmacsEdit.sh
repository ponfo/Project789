#!/bin/sh

emacs ThermalMaterial.f90                   \
      CustomConditions/FluxOnLine.f90       \
      CustomConditions/ConvectionOnLine.f90 \
      CustomElements/ThermalElement.f90     \
      CustomIO/DataInput.f90                \
      CustomIO/DataOutput.f90               \
      Thermal2DApplication.f90              \
      main.f90                              &
