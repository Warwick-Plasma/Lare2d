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

    INTEGER :: test
    LOGICAL :: second_call = .FALSE.

    IF (second_call) THEN
      vx = 0.0_num
      vy = 0.0_num
      vz = 0.0_num
      bx = 0.0_num
      by = 0.0_num
      bz = 0.0_num
      grav = 0.0_num
      rho = 1.0_num
      energy = 0.1_num
    END IF

    test = 1

    SELECT CASE(test)

!--------------------------------------------------------
      CASE(20) !2D X-point, beta_max ~ 1

        bc1 = 2
        bc2 = 2
        bc3 = 2
        visc1 = 0.0_num
        visc2 = 0.0_num    
        xbc_min = BC_OPEN
        xbc_max = BC_OPEN
        ybc_min = BC_OPEN
        ybc_max = BC_OPEN

        IF (second_call) THEN
          energy = 1.0_num
          DO iy = 0, ny
            DO ix = 0, nx
              energy(ix,iy) = energy(ix,iy)  &
                 * (1.0_num + 0.01_num * EXP(-(xc(ix)**2+yc(iy)**2)/0.01_num))
            END DO
          END DO
          DO iy = -1, ny+2
            DO ix = -2, nx+2
              bx(ix,iy) = xb(ix)
            END DO
          END DO

          DO iy = -2, ny+2
            DO ix = -1, nx+2
              by(ix,iy) = -yb(iy)
            END DO
          END DO
        END IF

!--------------------------------------------------------
      CASE (2) !1D Alfven wave test, beta ~1

        bc1 = 2
        bc2 = 2
        bc3 = 2
        visc1 = 0.0_num
        visc2 = 0.0_num    
        xbc_min = BC_OPEN
        xbc_max = BC_OPEN
        ybc_min = BC_PERIODIC
        ybc_max = BC_PERIODIC
        IF (second_call) THEN
         bx = 1.0_num
         energy = 1.0_num      
         DO iy = 0, ny
           DO ix = 0, nx
             vy(ix,iy) = 0.01_num * EXP(-(xb(ix)**2)/0.01_num)
           END DO
         END DO
        END IF

!--------------------------------------------------------
      CASE (1) !1D acoustic test

        bc1 = 1
        bc2 = 1
        bc3 = 1
        visc1 = 0.0_num
        visc2 = 0.0_num    
        xbc_min = BC_OPEN
        xbc_max = BC_OPEN
        ybc_min = BC_PERIODIC
        ybc_max = BC_PERIODIC
        IF (second_call) THEN
         energy = 1.0_num      
         DO iy = 0, ny
           DO ix = 0, nx
             vx(ix,iy) = 0.01_num * EXP(-(xb(ix)**2)/0.01_num)
           END DO
         END DO
        END IF

    END SELECT

    second_call = .TRUE.

  END SUBROUTINE set_initial_conditions

END MODULE initial_conditions
