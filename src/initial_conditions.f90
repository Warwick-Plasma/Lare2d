MODULE initial_conditions

  USE shared_data
  USE eos
  USE neutral

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: set_initial_conditions

CONTAINS

  !---------------------------------------------------------------------------
  ! This function sets up the initial condition for the code
  ! The variables which must be set are
  ! Rho - density
  ! V{x, y, z} - Velocities in x, y, z
  ! B{x, y, z} - Magnetic fields in x, y, z
  ! Energy - Specific internal energy
  ! Since temperature is a more intuitive quantity than specific internal energy
  ! There is a helper function Get_Energy which converts temperature to energy
  ! The syntax for this function is
  !
  ! CALL Get_Energy(density, temperature, equation_of_state,  &
  !     output_energy)
  !
  ! REAL(num) :: density - The density at point (ix, iy) on the grid
  ! REAL(num) :: temperature - The temperature at point (ix, iy) on the grid
  ! INTEGER :: equation_of_state - The code for the equation of state to use.
  !            The global equation of state for the code is eos_number
  ! REAL(num) :: output_energy - The specific internal energy returned by
  !              the routine
  !
  ! You may also need the neutral fraction. This can be calculated by a function
  ! call to  get_neutral(temperature, rho). This routine is in core/neutral.f90
  ! and requires the local temperature and mass density. For example to set
  ! xi_n to the neutral fraction use
  ! xi_n = get_neutral(temperature, rho)
  !---------------------------------------------------------------------------
  SUBROUTINE set_initial_conditions

    vx = 0.0_num
    vy = 0.0_num
    vz = 0.0_num
    bx = 0.0_num  
    by = 0.0_num 
    bz = 0.0_num
    energy = 1.e-5_num
    rho = 1.0_num

    grav = 0.0_num  

    DO iy = -1, ny+2
      DO ix = -1, nx+2 
        energy(ix,iy) = 1.e-4_num + 1.0_num * EXP(-((xc(ix)+1.0_num)**2+(yc(iy)-1.0_num)**2)/0.001_num)   
        bx(ix,iy) = -xc(ix) 
        by(ix,iy) = yc(iy) 
      END DO
    END DO

  END SUBROUTINE set_initial_conditions

END MODULE initial_conditions
