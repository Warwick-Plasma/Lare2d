PROGRAM lare2d

  USE shared_data
  USE deck
  USE initial_conditions
  USE setup
  USE boundary
  USE openboundary
  USE diagnostics
  USE lagran
  USE remap
  USE mpi_routines
  USE welcome

  IMPLICIT NONE

  INTEGER :: i = 0

  CALL minimal_init
  CALL mpi_minimal_init
  CALL welcome_message
 
  CALL read_deck 
  CALL set_boundary_conditions
  CALL mpi_initialise
  CALL after_control
  CALL open_files
  CALL grid
  time = 0.0_num
  CALL equilibrium    ! define initial profiles
  IF (restart) CALL restart_data
  CALL set_boundary_conditions
  CALL boundary_conditions

  CALL output_routines(i)
  CALL eta_calc


  DO
     IF ((i >= nsteps .AND. nsteps >=0) .OR. (time >= t_end)) EXIT
     i = i + 1
     CALL eta_calc
     CALL set_dt
     CALL lagrangian_step
     CALL eulerian_remap(i)
     IF (rke) CALL energy_correction
     IF (any_open) THEN
        CALL open_bcs
     ENDIF
     CALL boundary_conditions  
     CALL output_routines(i)
  END DO

  CALL mpi_close
  CALL close_files
  CALL MPI_FINALIZE(errcode)

END PROGRAM lare2d
