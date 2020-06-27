module HeatFluxM
  
  use UtilitiesM

  private
  public :: HeatFluxDT

  type :: HeatFluxDT
     integer(ikind), dimension(:)  , allocatable :: triangElemID
     real(rkind)   , dimension(:,:), allocatable :: triangGPoint
     real(rkind)   , dimension(:,:), allocatable :: triangFlux
     integer(ikind), dimension(:)  , allocatable :: quadElemID
     real(rkind)   , dimension(:,:), allocatable :: quadGPoint
     real(rkind)   , dimension(:,:), allocatable :: quadFlux
  end type HeatFluxDT
  
end module HeatFluxM
