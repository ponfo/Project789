program main
  
  use DataInputM
  use Structural3DApplicationM
  use StructuralStrategyM
  use GIDDataOutputM

  implicit none

  type(Structural3DApplicationDT)   :: application
  type(StructuralStrategyDT)        :: structuralStrategy

  call initFEM3D(application)
  call structuralStrategy%buildStrategyAndSolve(application%model)
  call initDataOutput()
  call printResults(resultName = 'Displacement'      &
       , step         = 1                            &
       , graphType    = 'Vector'                     &
       , locationName = 'onNodes'                    &
       , resultNumber = application%model%getnNode()  &
       , nDof         = 3                             &
       , component1   = application%model%dof         )
  
  call printResults(resultName = 'NormalStressOnTetras'                  &
       , type         = 'Tetrahedra'                                      &
       , step         = 1                                                 &
       , graphType    = 'Vector'                                          &
       , locationName = 'onGaussPoints'                                   &
       , gaussPoints  = application%model%normalStress%tetraGPoint       &
       , resultNumber = size(application%model%normalStress%tetraElemID) &
       , elemID       = application%model%normalStress%tetraElemID       &
       , component1   = application%model%normalStress%tetraNS(:,1)      &
       , component2   = application%model%normalStress%tetraNS(:,2)      &
       , component3   = application%model%normalStress%tetraNS(:,3)      )
  call printResults(resultName = 'NormalStressOnHexas'                   &
       , type         = 'Hexahedra'                                      &
       , step         = 1                                                &
       , graphType    = 'Vector'                                         &
       , locationName = 'onGaussPoints'                                  &
       , gaussPoints  = application%model%normalStress%hexaGPoint        &
       , resultNumber = size(application%model%normalStress%hexaElemID)  &
       , elemID       = application%model%normalStress%hexaElemID        &
       , component1   = application%model%normalStress%hexaNS(:,1)       &
       , component2   = application%model%normalStress%hexaNS(:,2)       &
       , component3   = application%model%normalStress%hexaNS(:,3)       )

  call printResults(resultName = 'ShearStressOnTetras'                  &
       , type         = 'Tetrahedra'                                    &
       , step         = 1                                               &
       , graphType    = 'Vector'                                        &
       , locationName = 'onGaussPoints'                                 &
       , gaussPoints  = application%model%shearStress%tetraGPoint       &
       , resultNumber = size(application%model%shearStress%tetraElemID) &
       , elemID       = application%model%shearStress%tetraElemID       &
       , component1   = application%model%shearStress%tetraShS(:,1)     &
       , component2   = application%model%shearStress%tetraShS(:,2)     &
       , component3   = application%model%shearStress%tetraShS(:,3)     )
  call printResults(resultName = 'ShearStressOnHexas'                  &
       , type         = 'Hexahedra'                                    &
       , step         = 1                                              &
       , graphType    = 'Vector'                                       &
       , locationName = 'onGaussPoints'                                &
       , gaussPoints  = application%model%shearStress%hexaGPoint        &
       , resultNumber = size(application%model%shearStress%hexaElemID)  &
       , elemID       = application%model%shearStress%hexaElemID        &
       , component1   = application%model%shearStress%hexaShS(:,1)      &
       , component2   = application%model%shearStress%hexaShS(:,2)      &
       , component3   = application%model%shearStress%hexaShS(:,3)      )

  call printResults(resultName = 'StrainStressOnTetras'           &
       , type         = 'Tetrahedra'                               &
       , step         = 1                                          &
       , graphType    = 'Vector'                                   &
       , locationName = 'onGaussPoints'                            &
       , gaussPoints  = application%model%strain%tetraGPoint       &
       , resultNumber = size(application%model%strain%tetraElemID) &
       , elemID       = application%model%strain%tetraElemID       &
       , component1   = application%model%strain%tetraEp(:,1)      &
       , component2   = application%model%strain%tetraEp(:,2)      &
       , component3   = application%model%strain%tetraEp(:,3)      )
  call printResults(resultName = 'StrainStressOnHexas'            &
       , type         = 'Hexahedra'                               &
       , step         = 1                                         &
       , graphType    = 'Vector'                                  &
       , locationName = 'onGaussPoints'                           &
       , gaussPoints  = application%model%strain%hexaGPoint        &
       , resultNumber = size(application%model%strain%hexaElemID)  &
       , elemID       = application%model%strain%hexaElemID        &
       , component1   = application%model%strain%hexaEp(:,1)       &
       , component2   = application%model%strain%hexaEp(:,2)       &
       , component3   = application%model%strain%hexaEp(:,3)       )
 
end program main
