program main
  
  use DataInputM
  use Thermal3DApplicationM
  use ThermalStrategyM
  use GIDDataOutputM

  implicit none

  type(Thermal3DApplicationDT)   :: application
  type(ThermalStrategyDT)        :: thermalStrategy

  call initFEM3D(application)
  call thermalStrategy%buildStrategyAndSolve(application%model)
  call initDataOutput()
  call printResults(resultName = 'Temperature'       &
       , step         = 1                            &
       , graphType    = 'Scalar'                     &
       , locationName = 'onNodes'                    &
       , resultNumber = application%model%getnNode()  &
       , component1   = application%model%dof         )
  call printResults(resultName = 'FluxOnTetras'                     &
       , type         = 'Tetrahedron'                               &
       , step         = 1                                           &
       , graphType    = 'Vector'                                    &
       , locationName = 'onGaussPoints'                             &
       , gaussPoints  = application%model%heatFlux%tetraGPoint       &
       , resultNumber = size(application%model%heatFlux%tetraElemID) &
       , elemID       = application%model%heatFlux%tetraElemID       &
       , component1   = application%model%heatFlux%tetraFlux(:,1)    &
       , component2   = application%model%heatFlux%tetraFlux(:,2)    &
       , component3   = application%model%heatFlux%tetraFlux(:,3)    )
 
end program main