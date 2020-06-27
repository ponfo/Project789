module ShearStressM
  
  use UtilitiesM
  
  implicit none
  private
  
  public :: ShearStressDT
  type ShearStressDT
     integer(ikind) , dimension(:)  , allocatable :: triangElemID
     real(rkind)    , dimension(:,:), allocatable :: triangGPoint
     real(rkind)    , dimension(:  ), allocatable :: triangShS
     integer(ikind) , dimension(:)  , allocatable :: quadElemID
     real(rkind)    , dimension(:,:), allocatable :: quadGPoint
     real(rkind)    , dimension(:  ), allocatable :: quadShS
  end type ShearStressDT
  
end module ShearStressM
