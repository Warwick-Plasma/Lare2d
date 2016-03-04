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
  !
  ! You may also need the neutral fraction. This can be calculated by a
  ! function call to get_neutral(temperature, rho, z). This routine is in
  ! core/neutral.f90 and requires the local temperature and mass density.
  ! For example to set xi_n to the neutral fraction use:
  !   xi_n = get_neutral(temperature, rho, z)
  !****************************************************************************

  SUBROUTINE set_initial_conditions

    INTEGER:: ix

    ! 1d MHD Brio-Wu
    ! typically run with 800 points until t=0.1
    !normalised equations

    gamma = 5.0_num / 3.0_num

    vx = 0.0_num
    vy = 0.0_num
    vz = 0.0_num
    bx = 0.75_num
    by = 0.0_num
    bz = 0.0_num

    DO ix = -1, nx+2
       IF (xc(ix) >= 0.0_num) THEN
          rho(ix,:) = 0.125_num
          energy(ix,:) = 0.1_num / (gamma - 1.0_num)
          by(ix,:) = -1.0_num
       ELSE
          rho(ix,:) = 1.0_num
          energy(ix,:) = 1.0_num / (gamma - 1.0_num)
          by(ix,:) = 1.0_num
       END IF
    END DO

    energy = energy / rho

  END SUBROUTINE set_initial_conditions



END MODULE initial_conditions
