module StructuralSchemeM

  use UtilitiesM

  use IntegratorPtrM

  use ElementPtrM
  
  use StructuralmodelM
  
  use SchemeM

  implicit none

  private
  public :: StructuralSchemeDT

  type, extends(NewSchemeDT) :: StructuralSchemeDT
   contains
     procedure, public              :: calculatePost
  end type StructuralSchemeDT
  
contains

  subroutine calculatePost(this, model)
    implicit none
    class(StructuralSchemeDT), intent(inout)             :: this
    class(StructuralModelDT) , intent(inout)             :: model
    logical                                              :: firstTetra = .true.
    logical                                              :: firstHexa = .true.
    integer(ikind)                                       :: iElem, iGauss, nNode, nTetra, nHexa
    integer(ikind)                                       :: tetraCounter, hexaCounter
    integer(ikind)                                       :: tetraPointCounter, hexaPointCounter
    integer(ikind)                                       :: nElem, nGauss, nPointsTetra, nPointsHexa
    real(rkind)          , dimension(:,:,:), allocatable :: localResultMat
    type(ElementPtrDT)                                   :: element
    type(IntegratorPtrDT)                                :: integrator
    write(*,*) '***  Structural Scheme ***'
    write(*,*) '*** Calculate Flux ***'
    nElem = model%getnElement()
    nPointsTetra = 0
    nTetra = 0
    nPointsHexa = 0
    nHexa = 0
    do iElem = 1, nElem
       element = model%getElement(iElem)
       nNode = element%getnNode()
       integrator = element%getIntegrator()
       nGauss = integrator%getIntegTerms()
       if(nNode == 4 .or. nNode == 10) then
          nPointsTetra = nPointsTetra + nGauss
          nTetra = nTetra + 1
          if(firstTetra) then
             model%normalStress%tetraGPoint = integrator%getGPointFull()
             model%shearStress%tetraGPoint = integrator%getGPointFull()
             model%strain%tetraGPoint = integrator%getGPointFull()
             firstTetra = .false.
          end if
       else if(nNode == 8 .or. nNode == 20) then
          nPointsHexa = nPointsHexa + nGauss
          nHexa = nHexa + 1
          if(firstHexa) then
             model%normalStress%hexaGPoint = integrator%getGPointFull()
             model%shearStress%hexaGPoint = integrator%getGPointFull()
             model%strain%hexaGPoint = integrator%getGPointFull()
             firstHexa = .false.
          end if
       end if
    end do
    allocate(model%normalStress%tetraElemID(nTetra))
    allocate(model%normalStress%tetraNS(nPointsTetra,3))
    allocate(model%normalStress%hexaElemID(nHexa))
    allocate(model%normalStress%hexaNS(nPointsHexa,3))
    allocate(model%shearStress%tetraElemID(nTetra))
    allocate(model%shearStress%tetraShS(nPointsTetra,3))
    allocate(model%shearStress%hexaElemID(nHexa))
    allocate(model%shearStress%hexaShS(nPointsHexa,3))
    allocate(model%strain%tetraElemID(nTetra))
    allocate(model%strain%tetraEp(nPointsTetra,3))
    allocate(model%strain%hexaElemID(nHexa))
    allocate(model%strain%hexaEp(nPointsHexa,3))
    tetraCounter = 0
    tetraPointCounter = 0
    hexaCounter = 0
    hexaPointCounter = 0
    do iElem = 1, nElem
       element = model%getElement(iElem)
       nNode = element%getnNode()
       call element%calculateResults(model%processInfo, localResultMat)
       nGauss = size(localResultMat,2)
       if(nNode == 4 .or. nNode == 10) then
          tetraCounter = tetraCounter + 1
          model%normalStress%tetraElemID(tetraCounter) = iElem
          model%shearStress%tetraElemID(tetraCounter) = iElem
          model%strain%tetraElemID(tetraCounter) = iElem
          do iGauss = 1, nGauss
             tetraPointCounter = tetraPointCounter + 1
             model%normalStress%tetraNS(tetraPointCounter,1) = localResultMat(1,iGauss,1)
             model%normalStress%tetraNS(tetraPointCounter,2) = localResultMat(1,iGauss,2)
             model%normalStress%tetraNS(tetraPointCounter,3) = localResultMat(1,iGauss,3)
             model%shearStress%tetraShS(tetraPointCounter,1) = localResultMat(2,iGauss,1)
             model%shearStress%tetraShS(tetraPointCounter,2) = localResultMat(2,iGauss,2)
             model%shearStress%tetraShS(tetraPointCounter,3) = localResultMat(2,iGauss,3)
             model%strain%tetraEp(tetraPointCounter,1)       = localResultMat(3,iGauss,1)
             model%strain%tetraEp(tetraPointCounter,2)       = localResultMat(3,iGauss,2)
             model%strain%tetraEp(tetraPointCounter,3)       = localResultMat(3,iGauss,3)
          end do
       else if(nNode == 8 .or. nNode == 20) then
          hexaCounter = hexaCounter + 1
          model%normalStress%hexaElemID(hexaCounter) = iElem
          model%shearStress%hexaElemID(hexaCounter) = iElem
          model%strain%hexaElemID(hexaCounter) = iElem
          do iGauss = 1, nGauss
             hexaPointCounter = hexaPointCounter + 1
             model%normalStress%hexaNS(hexaPointCounter,1) = localResultMat(1,iGauss,1)
             model%normalStress%hexaNS(hexaPointCounter,2) = localResultMat(1,iGauss,2)
             model%normalStress%hexaNS(hexaPointCounter,3) = localResultMat(1,iGauss,3)
             model%shearStress%hexaShS(hexaPointCounter,1) = localResultMat(2,iGauss,1)
             model%shearStress%hexaShS(hexaPointCounter,2) = localResultMat(2,iGauss,2)
             model%shearStress%hexaShS(hexaPointCounter,3) = localResultMat(2,iGauss,3)
             model%strain%hexaEp(hexaPointCounter,1)       = localResultMat(3,iGauss,1)
             model%strain%hexaEp(hexaPointCounter,2)       = localResultMat(3,iGauss,2)
             model%strain%hexaEp(hexaPointCounter,3)       = localResultMat(3,iGauss,3)
          end do
       end if
       deallocate(localResultMat)
    end do
  end subroutine calculatePost

end module StructuralSchemeM
