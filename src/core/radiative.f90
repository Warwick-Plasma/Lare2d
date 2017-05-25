MODULE radiative

  USE shared_data
  USE boundary
  USE neutral

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rad_losses

  INTEGER, PARAMETER :: n = 7   !set the number of temperature boundaries used in Q(T)
  INTEGER, PARAMETER  :: kmax = n - 7
  REAL(num), DIMENSION(:), ALLOCATABLE :: t_boundary, alpha, psi

CONTAINS


  SUBROUTINE rad_losses

    LOGICAL :: first_call = .TRUE.

    IF (first_call) THEN
      first_call = .FALSE.
      ALLOCATE (t_boundary(1:n), alpha(1:kmax), psi(1:kmax))
      CALL setup_loss_function
    END IF

  END SUBROUTINE rad_losses



  SUBROUTINE setup_loss_function

    REAL(num) :: frac

    t_boundary = (/0.02_num, 0.0398_num, 0.0794_num, 0.251_num, 0.562_num, 1.995_num, 10.0_num/) * 1e6_num

    frac = 2.0_num / 3.0_num
    alpha = (/0.0_num, 2.0_num, 0.0_num, -2.0_num, 0.0_num, frac/)

    psi = (/-21.85_num, -31.0_num, -21.2_num, -10.4_num, -21.94_num, -17.73_num/)
    psi = 10**(psi)

  END SUBROUTINE setup_loss_function



  FUNCTION heating(density, temperature)

    ! For a given density and temperature returns a user specific
    ! heating function.
    ! Input density is the mass density. 
    ! Both density and temperature are in Lare2d normalised units

    REAL(num), INTENT(IN) :: density, temperature
    REAL(num) :: heating

    ! First specify the heating in S.I. units
    heating = 0.0_num

    ! Convert to LareXd normalised units
    heating = heating * h_star

  END FUNCTION heating


END MODULE radiative
