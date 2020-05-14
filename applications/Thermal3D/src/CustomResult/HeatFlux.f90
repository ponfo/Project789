module HeatFluxM
  use UtilitiesM

  private
  public :: HeatFluxDT

  type :: HeatFluxDT
     integer(ikind), dimension(:)  , allocatable :: tetraElemID
     real(rkind)   , dimension(:,:), allocatable :: tetraGPoint
     real(rkind)   , dimension(:,:), allocatable :: tetraFlux
     integer(ikind), dimension(:)  , allocatable :: hexaElemID
     real(rkind)   , dimension(:,:), allocatable :: hexaGPoint
     real(rkind)   , dimension(:,:), allocatable :: hexaFlux
  end type HeatFluxDT
end module HeatFluxM
