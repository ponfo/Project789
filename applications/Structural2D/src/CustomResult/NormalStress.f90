module NormalStressM
  
  use UtilitiesM
  
  implicit none
  private
  
  public :: NormalStressDT
  type NormalStressDT
     integer(ikind) , dimension(:)  , allocatable :: triangElemID
     real(rkind)    , dimension(:,:), allocatable :: triangGPoint
     real(rkind)    , dimension(:,:), allocatable :: triangNS
     integer(ikind) , dimension(:)  , allocatable :: quadElemID
     real(rkind)    , dimension(:,:), allocatable :: quadGPoint
     real(rkind)    , dimension(:,:), allocatable :: quadNS
  end type NormalStressDT

end module NormalStressM
