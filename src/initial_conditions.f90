MODULE initial_conditions

  USE shared_data
  USE neutral
  USE diagnostics

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: set_initial_conditions

CONTAINS

  !****************************************************************************
  ! This function sets up the initial condition for the code
  ! The variables which must be set are:
  !   rho - density
  !   v{x,y,z} - Velocities in x, y, z
  !   b{x,y,z} - Magnetic fields in x, y, z
  !   energy - Specific internal energy
  !   grav - Gravity
  !
  ! You may also need the neutral fraction. This can be calculated by a
  ! function call to get_neutral(temperature, rho). This routine is in
  ! core/neutral.f90 and requires the local temperature and mass density.
  ! For example to set xi_n to the neutral fraction use:
  !   xi_n = get_neutral(temperature, rho)
  ! 
  ! If using Hall_MHD then you must specific lambda_i in this routine
  !****************************************************************************


  SUBROUTINE set_initial_conditions

    REAL(num), DIMENSION(:,:), ALLOCATABLE :: temperature

    ! Below are all the variables which must be defined and their sizes

    vx(-2:nx+2, -2:ny+2) = 0.0_num
    vy(-2:nx+2, -2:ny+2) = 0.0_num
    vz(-2:nx+2, -2:ny+2) = 0.0_num

    bx(-2:nx+2, -1:ny+2) = 0.0_num
    by(-1:nx+2, -2:ny+2) = 0.0_num
    bz(-1:nx+2, -1:ny+2) = 0.0_num

    rho(-1:nx+2, -1:ny+2) = 1.0_num
    energy(-1:nx+2, -1:ny+2) = 0.1_num

    grav(-1:ny+2) = 0.0_num

    ! If defining the initial conditions using temperature then use
    ALLOCATE(temperature(-1:nx+2, -1:ny+2))
    temperature(-1:nx+2, -1:ny+2) = 0.5_num
    ! Then fix the energy, for a fully ionised plasma, from
    energy(:,:) = 2.0_num * temperature(:,:) / (gamma - 1.0_num) 

    ! If neutrals included xi_n is a function of temperature so iteration required
    ! Iteration not shown in this example - see examples in Old directory
    ! Set the neutral fraction if needed
    IF (eos_number /= EOS_IDEAL) THEN
      DO iy = -1, ny+2
        DO ix = -1, nx+2
          xi_n(ix,iy) = get_neutral(temperature(ix,iy), rho(ix,iy))
        END DO
      END DO
    END IF

    IF (hall_mhd) THEN
      lambda_i(0:nx, 0:ny) = 1.0_num
    END IF

    ! If probe points needed add them here
    CALL add_probe(0.0_num, 0.0_num)

    DEALLOCATE(temperature)

  END SUBROUTINE set_initial_conditions

END MODULE initial_conditions
