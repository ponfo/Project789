module ThermalMaterialM
  
  use UtilitiesM
  use PropertyM
  
  implicit none
  
  private
  public :: ThermalMaterialDT, thermalMaterial
  
  type :: ThermalMaterialDT
     real(rkind), dimension(3) :: conductivity
   contains
     procedure :: init
  end type ThermalMaterialDT
  
  interface thermalMaterial
     procedure :: constructor
  end interface thermalMaterial
  
contains
  
  type(ThermalMaterialDT) function constructor(kx,ky,kz)
    implicit none
    real(rkind), intent(in) :: kx, ky, kz
    call constructor%init(kx, ky, kz)
  end function constructor
  
  subroutine init(this, kx, ky, kz)
    implicit none
    real(rkind), intent(in) :: kx, ky, kz
    class(ThermalMaterialDT), intent(inout) :: this
    this%conductivity = (/kx, ky, kz/)
  end subroutine init
  
end module ThermalMaterialM
