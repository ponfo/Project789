module CFDMaterialM
  use UtilitiesM
  use PropertyM
  
  implicit none
  
  private
  public :: CFDMaterialDT, cfdMaterial
  
  type :: CFDMaterialDT
     real(rkind) :: R
     real(rkind) :: gamma
     real(rkind) :: mu
     real(rkind) :: k
   contains
     procedure :: init
  end type CFDMaterialDT
  
  interface cfdMaterial
     procedure :: constructor
  end interface cfdMaterial
  
contains
  
  type(CFDMaterialDT) function constructor(R, gamma, mu, k)
    implicit none
    real(rkind), intent(in) :: R
    real(rkind), intent(in) :: gamma
    real(rkind), intent(in) :: mu
    real(rkind), intent(in) :: k
    call constructor%init(R, gamma, mu, k)
  end function constructor

  subroutine init(this, R, gamma, mu, k)
    implicit none
    class(CFDMaterialDT), intent(inout) :: this
    real(rkind)         , intent(in)    :: R
    real(rkind)         , intent(in)    :: gamaa
    real(rkind)         , intent(in)    :: mu
    real(rkind)         , intent(in)    :: k
    this%R     = R
    this%gamma = gamma
    this%mu    = mu
    this%k     = k
  end subroutine init
  
end module CFDMaterialM
