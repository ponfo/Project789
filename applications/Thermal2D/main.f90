program main
  
  use DataInputM
  use Thermal2DApplicationM
  use ThermalStrategyM
  use SolvingStrategyM
  use GIDDataOutputM

  implicit none

  type(Thermal2DApplicationDT)   :: application
  type(ThermalStrategyDT)        :: strategy
  type(SolvingStrategyDT)        :: solvingStrategy
 
  call initFEM2D(application)
  solvingStrategy = InitSolvingStrategy(strategy, application%model)
  call solvingStrategy%useStrategy()
  call initDataOutput()
  call printResults(resultName = 'Temperature'       &
       , step         = 1                            &
       , graphType    = 'Scalar'                     &
       , locationName = 'onNodes'                    &
       , resultNumber = application%model%getnNode()  &
       , component1   = application%model%dof         )
  call printResults(resultName = 'FluxOnTriangs'                     &
       , type         = 'Triangle'                                   &
       , step         = 1                                            &
       , graphType    = 'Vector'                                     &
       , locationName = 'onGaussPoints'                              &
       , gaussPoints  = application%model%heatFlux%triangGPoint       &
       , resultNumber = size(application%model%heatFlux%triangElemID) &
       , elemID       = application%model%heatFlux%triangElemID       &
       , component1   = application%model%heatFlux%triangFlux(:,1)    &
       , component2   = application%model%heatFlux%triangFlux(:,2)    )
  call printResults(resultName = 'FluxOnQuads'                      &
       , type         = 'Quadrilateral'                             &
       , step         = 1                                           &
       , graphType    = 'Vector'                                    &
       , locationName = 'onGaussPoints'                             &
       , gaussPoints  = application%model%heatFlux%quadGPoint        &
       , resultNumber = size(application%model%heatFlux%quadElemID)  &
       , elemID       = application%model%heatFlux%quadElemID        &
       , component1   = application%model%heatFlux%quadFlux(:,1)     &
       , component2   = application%model%heatFlux%quadFlux(:,2)    )
 
end program main
