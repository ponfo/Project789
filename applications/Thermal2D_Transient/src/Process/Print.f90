module PrintM

  use UtilitiesM

  use HeatFluxM
  
  use GIDDataOutputM

  use ProcessM
  
  implicit none

  private
  public :: PrintDT

  type, extends(NewProcessDT) :: PrintDT
   contains
     procedure :: initPrint
     procedure :: useProcess => use
     procedure :: print
  end type PrintDT

contains

  subroutine initPrint(this)
    implicit none
    class(PrintDT), intent(inout) :: this
    call initDataOutput()
  end subroutine initPrint
  
  subroutine print(this, dof, heatFlux, step)
    implicit none
    class(PrintDT)           , intent(inout) :: this
    class(HeatFluxDT)        , intent(inout) :: heatFlux
    real(rkind), dimension(:), intent(inout) :: dof
    integer(ikind)           , intent(inout) :: step
    call printResults(resultName = 'Temperature'      &
         , step         = step                        &
         , graphType    = 'Scalar'                    &
         , locationName = 'onNodes'                   &
         , resultNumber = size(dof)                   &
         , component1   = dof                         )
    call printResults(resultName = 'FluxOnTriangs'    &
         , type         = 'Triangle'                  &
         , step         = step                        &
         , graphType    = 'Vector'                    &
         , locationName = 'onGaussPoints'             &
         , gaussPoints  = heatFlux%triangGPoint       &
         , resultNumber = size(heatFlux%triangElemID) &
         , elemID       = heatFlux%triangElemID       &
         , component1   = heatFlux%triangFlux(:,1)    &
         , component2   = heatFlux%triangFlux(:,2)    )
    call printResults(resultName = 'FluxOnQuads'      &
         , type         = 'Quadrilateral'             &
         , step         = step                        &
         , graphType    = 'Vector'                    &
         , locationName = 'onGaussPoints'             &
         , gaussPoints  = heatFlux%quadGPoint         &
         , resultNumber = size(heatFlux%quadElemID)   &
         , elemID       = heatFlux%quadElemID         &
         , component1   = heatFlux%quadFlux(:,1)      &
         , component2   = heatFlux%quadFlux(:,2)      )
  end subroutine print

  subroutine use(this)
    implicit none
    class(PrintDT), intent(inout) :: this
  end subroutine use
  
end module PrintM
