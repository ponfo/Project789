module NormalStressM
  use UtilitiesM
  
  implicit none
  private
  
  public :: NormalStressDT
  type NormalStressDT
     integer(ikind) , dimension(:)  , allocatable :: tetraElemID
     real(rkind)    , dimension(:,:), allocatable :: tetraGPoint
     real(rkind)    , dimension(:,:), allocatable :: tetraNS
     integer(ikind) , dimension(:)  , allocatable :: hexaElemID
     real(rkind)    , dimension(:,:), allocatable :: hexaGPoint
     real(rkind)    , dimension(:,:), allocatable :: hexaNS
  end type NormalStressDT

end module NormalStressM
