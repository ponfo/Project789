module ResultsM
  
  use UtilitiesM
  
  implicit none
  private
  
  public :: ResultsDT
  type ResultsDT
     real(rkind) , dimension(:)  , allocatable :: density
     real(rkind) , dimension(:)  , allocatable :: internalEnergy
     real(rkind) , dimension(:)  , allocatable :: pressure
     real(rkind) , dimension(:)  , allocatable :: mach
     real(rkind) , dimension(:)  , allocatable :: temperature
     real(rkind) , dimension(:,:), allocatable :: velocity
  end type ResultsDT

end module ResultsM
