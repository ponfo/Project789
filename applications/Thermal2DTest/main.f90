program main
  
  use DataInputM
  use Thermal2DApplicationM
  use DirectBuilderAndSolverM
  use GIDDataOutputM

  implicit none

  type(Thermal2DApplicationDT)   :: application
  type(DirectBuilderAndSolverDT) :: builderAndSolver

  call initFEM2D(application)
  call builderAndSolver%buildAndSolve(application%model)
  call printResults(resultName = 'Temperature'  &
       , step         = 1                            &
       , graphType    = 'Scalar'                     &
       , locationName = 'onNodes'                    &
       , resultNumber = application%model%getnNode()  &
       , component1   = application%model%dof         )

end program main
