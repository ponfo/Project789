module HeatFluxM
  use UtilitiesM

  private
  public :: HeatFluxDT

  type :: HeatFluxDT
     integer(ikind), dimension(:)  , allocatable :: tetraElemID
     real(rkind)   , dimension(:,:), allocatable :: tetraGPoint
     real(rkind)   , dimension(:,:), allocatable :: tetraFlux
  end type HeatFluxDT
end module HeatFluxM
