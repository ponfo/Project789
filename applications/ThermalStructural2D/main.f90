program main
  
  use DataInputM
  use ThermalStruct2DApplicationM
  use SolvingStrategyM
  use ThermalStrategyM
  use StructuralStrategyM
  use GIDDataOutputM

  implicit none

  type(ThermalStruct2DApplicationDT) :: application
  type(SolvingStrategyDT)            :: solvingStrategy
  type(ThermalStrategyDT)            :: thermalStrategy
  type(StructuralStrategyDT)         :: structuralStrategy

  call initFEM2D(application)
  solvingStrategy = InitSolvingStrategy(thermalStrategy, application%thermalModel)
  call solvingStrategy%useStrategy()
  call initDataOutput()
  call printResults(resultName = 'Temperature'              &
       , step         = 1                                   &
       , graphType    = 'Scalar'                            &
       , locationName = 'onNodes'                           &
       , resultNumber = application%thermalModel%getnNode()  &
       , component1   = application%thermalModel%dof         )
  call printResults(resultName = 'FluxOnTriangs'                            &
       , type         = 'Triangle'                                          &
       , step         = 1                                                   &
       , graphType    = 'Vector'                                            &
       , locationName = 'onGaussPoints'                                     &
       , gaussPoints  = application%thermalModel%heatFlux%triangGPoint       &
       , resultNumber = size(application%thermalModel%heatFlux%triangElemID) &
       , elemID       = application%thermalModel%heatFlux%triangElemID       &
       , component1   = application%thermalModel%heatFlux%triangFlux(:,1)    &
       , component2   = application%thermalModel%heatFlux%triangFlux(:,2)    )
  call printResults(resultName = 'FluxOnQuads'                             &
       , type         = 'Quadrilateral'                                    &
       , step         = 1                                                  &
       , graphType    = 'Vector'                                           &
       , locationName = 'onGaussPoints'                                    &
       , gaussPoints  = application%thermalModel%heatFlux%quadGPoint        &
       , resultNumber = size(application%thermalModel%heatFlux%quadElemID)  &
       , elemID       = application%thermalModel%heatFlux%quadElemID        &
       , component1   = application%thermalModel%heatFlux%quadFlux(:,1)     &
       , component2   = application%thermalModel%heatFlux%quadFlux(:,2)    )

  call application%transitionToStructural()
  
  solvingStrategy = InitSolvingStrategy(structuralStrategy, application%structuralModel)
  call solvingStrategy%useStrategy()
  
  call printResults(resultName = 'Displacement'                &
       , step         = 1                                      &
       , graphType    = 'Vector'                               &
       , locationName = 'onNodes'                              &
       , resultNumber = application%structuralModel%getnNode()  &
       , component1   = application%structuralModel%dof         )
  
  call printResults(resultName = 'NormalStressOnTriangs'                            &
       , type         = 'Triangle'                                                  &
       , step         = 1                                                           &
       , graphType    = 'Vector'                                                    &
       , locationName = 'onGaussPoints'                                             &
       , gaussPoints  = application%structuralModel%normalStress%triangGPoint        &
       , resultNumber = size(application%structuralModel%normalStress%triangElemID)  &
       , elemID       = application%structuralModel%normalStress%triangElemID        &
       , component1   = application%structuralModel%normalStress%triangNS(:,1)       &
       , component2   = application%structuralModel%normalStress%triangNS(:,2)       )
  call printResults(resultName = 'NormalStressOnQuads'                            &
       , type         = 'Quadrilateral'                                           &
       , step         = 1                                                         &
       , graphType    = 'Vector'                                                  &
       , locationName = 'onGaussPoints'                                           &
       , gaussPoints  = application%structuralModel%normalStress%quadGPoint        &
       , resultNumber = size(application%structuralModel%normalStress%quadElemID)  &
       , elemID       = application%structuralModel%normalStress%quadElemID        &
       , component1   = application%structuralModel%normalStress%quadNS(:,1)       &
       , component2   = application%structuralModel%normalStress%quadNS(:,2)       )

  call printResults(resultName = 'ShearStressOnTriangs'                           &
       , type         = 'Triangle'                                                &
       , step         = 1                                                         &
       , graphType    = 'Scalar'                                                  &
       , locationName = 'onGaussPoints'                                           &
       , gaussPoints  = application%structuralModel%shearStress%triangGPoint       &
       , resultNumber = size(application%structuralModel%shearStress%triangElemID) &
       , elemID       = application%structuralModel%shearStress%triangElemID       &
       , component1   = application%structuralModel%shearStress%triangShS          )
  call printResults(resultName = 'ShearStressOnQuads'                            &
       , type         = 'Quadrilateral'                                          &
       , step         = 1                                                        &
       , graphType    = 'Scalar'                                                 &
       , locationName = 'onGaussPoints'                                          &
       , gaussPoints  = application%structuralModel%shearStress%quadGPoint        &
       , resultNumber = size(application%structuralModel%shearStress%quadElemID)  &
       , elemID       = application%structuralModel%shearStress%quadElemID        &
       , component1   = application%structuralModel%shearStress%quadShS           )

  call printResults(resultName = 'StrainStressOnTriangs'                     &
       , type         = 'Triangle'                                           &
       , step         = 1                                                    &
       , graphType    = 'Vector'                                             &
       , locationName = 'onGaussPoints'                                      &
       , gaussPoints  = application%structuralModel%strain%triangGPoint       &
       , resultNumber = size(application%structuralModel%strain%triangElemID) &
       , elemID       = application%structuralModel%strain%triangElemID       &
       , component1   = application%structuralModel%strain%triangEp(:,1)      &
       , component2   = application%structuralModel%strain%triangEp(:,2)      )
  call printResults(resultName = 'StrainStressOnQuads'                      &
       , type         = 'Quadrilateral'                                     &
       , step         = 1                                                   &
       , graphType    = 'Vector'                                            &
       , locationName = 'onGaussPoints'                                     &
       , gaussPoints  = application%structuralModel%strain%quadGPoint        &
       , resultNumber = size(application%structuralModel%strain%quadElemID)  &
       , elemID       = application%structuralModel%strain%quadElemID        &
       , component1   = application%structuralModel%strain%quadEp(:,1)       &
       , component2   = application%structuralModel%strain%quadEp(:,2)       )

  call finishProgram()
 
end program main
