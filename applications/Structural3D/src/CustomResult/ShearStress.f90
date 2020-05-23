module ShearStressM
  use UtilitiesM
  
  implicit none
  private
  
  public :: ShearStressDT
  type ShearStressDT
     integer(ikind) , dimension(:)  , allocatable :: tetraElemID
     real(rkind)    , dimension(:,:), allocatable :: tetraGPoint
     real(rkind)    , dimension(:,:), allocatable :: tetraShS
     integer(ikind) , dimension(:)  , allocatable :: hexaElemID
     real(rkind)    , dimension(:,:), allocatable :: hexaGPoint
     real(rkind)    , dimension(:,:), allocatable :: hexaShS
  end type ShearStressDT
  
end module ShearStressM
