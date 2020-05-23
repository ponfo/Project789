module StrainM
  use UtilitiesM
  
  implicit none
  private
  
  public :: StrainDT
  type StrainDT
     integer(ikind) , dimension(:)  , allocatable :: tetraElemID
     real(rkind)    , dimension(:,:), allocatable :: tetraGPoint
     real(rkind)    , dimension(:,:), allocatable :: tetraEp
     integer(ikind) , dimension(:)  , allocatable :: hexaElemID
     real(rkind)    , dimension(:,:), allocatable :: hexaGPoint
     real(rkind)    , dimension(:,:), allocatable :: hexaEp
  end type StrainDT

end module StrainM
