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
     real(rkind) :: V_inf
     real(rkind) :: T_inf
     real(rkind) :: P_inf
     real(rkind) :: Vc
   contains
     procedure :: init
  end type CFDMaterialDT
  
  interface cfdMaterial
     procedure :: constructor
  end interface cfdMaterial
  
contains
  
  type(CFDMaterialDT) function constructor(R, gamma, mu, k, V, T, P)
    implicit none
    real(rkind), intent(in) :: R
    real(rkind), intent(in) :: gamma
    real(rkind), intent(in) :: mu
    real(rkind), intent(in) :: k
    real(rkind), intent(in) :: V
    real(rkind), intent(in) :: T
    real(rkind), intent(in) :: P
    call constructor%init(R, gamma, mu, k, V, T, P)
  end function constructor

  subroutine init(this, R, gamma, mu, k, V, T, P)
    implicit none
    class(CFDMaterialDT), intent(inout) :: this
    real(rkind)         , intent(in)    :: R
    real(rkind)         , intent(in)    :: gamma
    real(rkind)         , intent(in)    :: mu
    real(rkind)         , intent(in)    :: k
    real(rkind)         , intent(in)    :: V
    real(rkind)         , intent(in)    :: T
    real(rkind)         , intent(in)    :: P
    this%R     = R
    this%gamma = gamma
    this%mu    = mu
    this%k     = k
    this%V_inf = V
    this%T_inf = T
    this%P_inf = P
    this%Vc    = sqrt(gamma*R*T)
  end subroutine init
  
end module CFDMaterialM
