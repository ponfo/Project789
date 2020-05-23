module StructuralSchemeM

  use UtilitiesM
  use DebuggerM

  use IntegratorPtrM

  use ElementPtrM
  
  use StructuralmodelM
  
  use SchemeM

  use ProcessM

  implicit none

  private
  public :: StructuralSchemeDT

  type, extends(NewSchemeDT) :: StructuralSchemeDT
   contains
     procedure, public              :: calculatePost
     procedure, nopass              :: integrate => integrator
  end type StructuralSchemeDT
  
contains

  subroutine calculatePost(this, model)
    implicit none
    class(StructuralSchemeDT), intent(inout)             :: this
    class(StructuralModelDT) , intent(inout)             :: model
    logical                                              :: firstTriang = .true.
    logical                                              :: firstQuad = .true.
    integer(ikind)                                       :: iElem, iGauss, nNode, nTriang, nQuad
    integer(ikind)                                       :: triangCounter, quadCounter
    integer(ikind)                                       :: triangPointCounter, quadPointCounter
    integer(ikind)                                       :: nElem, nGauss, nPointsTriang, nPointsQuad
    real(rkind)          , dimension(:,:,:), allocatable :: localResultMat
    type(ElementPtrDT)                                   :: element
    type(IntegratorPtrDT)                                :: integrator
    write(*,*) '***  Structural Scheme ***'
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
       nGauss = integrator%getIntegTerms()
       if(nNode == 3 .or. nNode == 6) then
          nPointsTriang = nPointsTriang + nGauss
          nTriang = nTriang + 1
          if(firstTriang) then
             model%normalStress%triangGPoint = integrator%getGPointFull()
             model%shearStress%triangGPoint = integrator%getGPointFull()
             model%strain%triangGPoint = integrator%getGPointFull()
             firstTriang = .false.
          end if
       else if(nNode == 4 .or. nNode == 8) then
          nPointsQuad = nPointsQuad + nGauss
          nQuad = nQuad + 1
          if(firstQuad) then
             model%normalStress%quadGPoint = integrator%getGPointFull()
             model%shearStress%quadGPoint = integrator%getGPointFull()
             model%strain%quadGPoint = integrator%getGPointFull()
             firstQuad = .false.
          end if
       end if
    end do
    allocate(model%normalStress%triangElemID(nTriang))
    allocate(model%normalStress%triangNS(nPointsTriang,2))
    allocate(model%normalStress%quadElemID(nQuad))
    allocate(model%normalStress%quadNS(nPointsQuad,2))
    allocate(model%shearStress%triangElemID(nTriang))
    allocate(model%shearStress%triangShS(nPointsTriang))
    allocate(model%shearStress%quadElemID(nQuad))
    allocate(model%shearStress%quadShS(nPointsQuad))
    allocate(model%strain%triangElemID(nTriang))
    allocate(model%strain%triangEp(nPointsTriang,2))
    allocate(model%strain%quadElemID(nQuad))
    allocate(model%strain%quadEp(nPointsQuad,2))
    triangCounter = 0
    triangPointCounter = 0
    quadCounter = 0
    quadPointCounter = 0
    do iElem = 1, nElem
       element = model%getElement(iElem)
       nNode = element%getnNode()
       call element%calculateResults(model%processInfo, localResultMat)
       nGauss = size(localResultMat,2)
       if(nNode == 3 .or. nNode == 6) then
          triangCounter = triangCounter + 1
          model%normalStress%triangElemID(triangCounter) = iElem
          model%shearStress%triangElemID(triangCounter) = iElem
          model%strain%triangElemID(triangCounter) = iElem
          do iGauss = 1, nGauss
             triangPointCounter = triangPointCounter + 1
             model%normalStress%triangNS(triangPointCounter,1) = localResultMat(1,iGauss,1)
             model%normalStress%triangNS(triangPointCounter,2) = localResultMat(1,iGauss,2)
             model%shearStress%triangShS(triangPointCounter)   = localResultMat(2,iGauss,1)
             model%strain%triangEp(triangPointCounter,1)       = localResultMat(3,iGauss,1)
             model%strain%triangEp(triangPointCounter,2)       = localResultMat(3,iGauss,2)
          end do
       else if(nNode == 4 .or. nNode == 8) then
          quadCounter = quadCounter + 1
          model%normalStress%quadElemID(quadCounter) = iElem
          model%shearStress%quadElemID(quadCounter) = iElem
          model%strain%quadElemID(quadCounter) = iElem
          do iGauss = 1, nGauss
             quadPointCounter = quadPointCounter + 1
             model%normalStress%quadNS(quadPointCounter,1) = localResultMat(1,iGauss,1)
             model%normalStress%quadNS(quadPointCounter,2) = localResultMat(1,iGauss,2)
             model%shearStress%quadShS(quadPointCounter)   = localResultMat(2,iGauss,1)
             model%strain%quadEp(quadPointCounter,1)       = localResultMat(3,iGauss,1)
             model%strain%quadEp(quadPointCounter,2)       = localResultMat(3,iGauss,2)
          end do
       end if
       deallocate(localResultMat)
    end do
  end subroutine calculatePost
  
  subroutine integrator(this, dt)
    implicit none
    class(NewProcessDT), intent(inout) :: this
    real(rkind)        , intent(in)    :: dt
  end subroutine integrator

end module StructuralSchemeM
