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
     real(rkind) :: Vx_inf
     real(rkind) :: Vy_inf
     real(rkind) :: T_inf
     real(rkind) :: rho
     real(rkind) :: P_inf
     real(rkind) :: Vc
   contains
     procedure :: init
  end type CFDMaterialDT
  
  interface cfdMaterial
     procedure :: constructor
  end interface cfdMaterial
  
contains
  
  type(CFDMaterialDT) function constructor(R, gamma, mu, k, Vx, Vy, T, rho)
    implicit none
    real(rkind), intent(in) :: R
    real(rkind), intent(in) :: gamma
    real(rkind), intent(in) :: mu
    real(rkind), intent(in) :: k
    real(rkind), intent(in) :: Vx
    real(rkind), intent(in) :: Vy
    real(rkind), intent(in) :: T
    real(rkind), intent(in) :: rho
    call constructor%init(R, gamma, mu, k, Vx, Vy, T, rho)
  end function constructor

  subroutine init(this, R, gamma, mu, k, Vx, Vy, T, rho)
    implicit none
    class(CFDMaterialDT), intent(inout) :: this
    real(rkind)         , intent(in)    :: R
    real(rkind)         , intent(in)    :: gamma
    real(rkind)         , intent(in)    :: mu
    real(rkind)         , intent(in)    :: k
    real(rkind)         , intent(in)    :: Vx
    real(rkind)         , intent(in)    :: Vy
    real(rkind)         , intent(in)    :: T
    real(rkind)         , intent(in)    :: rho
    this%R      = R
    this%gamma  = gamma
    this%mu     = mu
    this%k      = k
    this%Vx_inf = Vx
    this%Vy_inf = Vy
    this%T_inf  = T
    this%rho    = rho
    this%P_inf  = rho*R*T
    this%Vc     = sqrt(gamma*R*T)
  end subroutine init
  
end module CFDMaterialM
