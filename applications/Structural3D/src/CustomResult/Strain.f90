module StrainM
  use UtilitiesM
  
  implicit none
  private
  
  public :: StrainDT
  type StrainDT
     integer(ikind) , dimension(:)  , allocatable :: triangElemID
     real(rkind)    , dimension(:,:), allocatable :: triangGPoint
     real(rkind)    , dimension(:,:), allocatable :: triangEp
     integer(ikind) , dimension(:)  , allocatable :: quadElemID
     real(rkind)    , dimension(:,:), allocatable :: quadGPoint
     real(rkind)    , dimension(:,:), allocatable :: quadEp
  end type StrainDT

end module StrainM
