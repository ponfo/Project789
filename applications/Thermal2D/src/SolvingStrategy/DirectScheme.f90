module DirectSchemeM

  use UtilitiesM
  use DebuggerM

  use IntegratorPtrM

  use ElementPtrM
  
  use ThermalmodelM
  
  use SchemeM

  implicit none

  private
  public :: DirectSchemeDT

  type, extends(NewSchemeDT) :: DirectSchemeDT
   contains
     procedure, public              :: calculateFlux
  end type DirectSchemeDT
  
contains

  subroutine calculateFlux(this, model)
    implicit none
    class(DirectSchemeDT), intent(inout)               :: this
    class(ThermalModelDT), intent(inout)               :: model
    logical                                            :: firstTriang = .true.
    logical                                            :: firstQuad = .true.
    integer(ikind)                                     :: iElem, iGauss, nNode, nTriang, nQuad
    integer(ikind)                                     :: triangCounter, quadCounter
    integer(ikind)                                     :: triangPointCounter, quadPointCounter
    integer(ikind)                                     :: nElem, nGauss, nPointsTriang, nPointsQuad
    real(rkind)          , dimension(:,:), allocatable :: localResultMat
    type(ElementPtrDT)                                 :: element
    type(IntegratorPtrDT)                              :: integrator
    write(*,*) '***  Direct Scheme ***'
    write(*,*) '*** Calculate Flux ***'
    nElem = model%getnElement()
    nPointsTriang = 0
    nTriang = 0
    nPointsQuad = 0
    nQuad = 0
    do iElem = 1, nElem
       element = model%getElement(iElem)
       nNode = element%getnNode()
       integrator = element%getIntegrator()
       nGauss = integrator%ptr%integTerms
       if(nNode == 3 .or. nNode == 6) then
          nPointsTriang = nPointsTriang + nGauss
          nTriang = nTriang + 1
          if(firstTriang) then
             allocate(model%heatFlux%triangGPoint(nGauss,2))
             model%heatFlux%triangGPoint = integrator%ptr%gPoint
             firstTriang = .false.
          end if
       else if(nNode == 4 .or. nNode == 8) then
          nPointsQuad = nPointsQuad + nGauss
          nQuad = nQuad + 1
          if(firstQuad) then
             allocate(model%heatFlux%quadGPoint(nGauss,2))
             model%heatFlux%quadGPoint = integrator%ptr%gPoint
             firstQuad = .false.
          end if
       end if
    end do
    allocate(model%heatFlux%triangElemID(nTriang))
    allocate(model%heatFlux%triangFlux(nPointsTriang,2))
    allocate(model%heatFlux%quadElemID(nQuad))
    allocate(model%heatFlux%quadFlux(nPointsQuad,2))
    triangCounter = 0
    triangPointCounter = 0
    quadCounter = 0
    quadPointCounter = 0
    do iElem = 1, nElem
       element = model%getElement(iElem)
       nNode = element%getnNode()
       call element%calculateResults(localResultMat)
       nGauss = size(localResultMat,1)
       if(nNode == 3 .or. nNode == 6) then
          triangCounter = triangCounter + 1
          model%heatFlux%triangElemID(triangCounter) = iElem
          do iGauss = 1, nGauss
             triangPointCounter = triangPointCounter + 1
             model%heatFlux%triangFlux(triangPointCounter,1) = localResultMat(iGauss,1)
             model%heatFlux%triangFlux(triangPointCounter,2) = localResultMat(iGauss,2)
          end do
       else if(nNode == 4 .or. nNode == 8) then
          quadCounter = quadCounter + 1
          model%heatFlux%quadElemID(quadCounter) = iElem
          do iGauss = 1, nGauss
             quadPointCounter = quadPointCounter + 1
             model%heatFlux%quadFlux(quadPointCounter,1) = localResultMat(iGauss,1)
             model%heatFlux%quadFlux(quadPointCounter,2) = localResultMat(iGauss,2)
          end do
       end if
       deallocate(localResultMat)
    end do
  end subroutine calculateFlux

end module DirectSchemeM