MODULE radiative

  USE shared_data
  USE boundary
  USE neutral

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rad_losses

CONTAINS

  SUBROUTINE rad_losses

  END SUBROUTINE


  SUBROUTINE loss_function(density, temperature, xi, height, rad)

    ! Returns the normalised RTV losses
    ! Input density is the mass density 
    ! Both density and temperature are in Lare2d normalised units
    ! To replace with a different loss function note that this
    ! routine first calculates the S.I. loss rate as defined in the manual.
    ! This is then converted to LareXd units by multiplying by h_star

    REAL(num), INTENT(IN) :: density, temperature, xi, height
    REAL(num), INTENT(OUT) :: rad
    REAL(num), DIMENSION(7) :: trange = (/0.02_num, 0.0398_num, 0.0794_num, &
        0.251_num, 0.562_num, 1.995_num, 10.0_num/)
    REAL(num), DIMENSION(6) :: psi = (/1.2303_num, 870.96_num, 5.496_num, &
        0.3467_num, 1.0_num, 1.6218_num/)
    REAL(num), DIMENSION(6) :: alpha = (/0.0_num, 2.0_num, 0.0_num, &
        -2.0_num, 0.0_num, -2.0_num / 3.0_num/)
    REAL(num) :: tmk, factor
    INTEGER :: i

    rad = 0.0_num

    IF (.NOT. radiation) RETURN
    IF (height < 1.0_num) RETURN

    tmk = temperature * t2tmk
    IF (tmk < trange(1) .OR. tmk > trange(7)) RETURN

    DO i = 1, 6
      IF (tmk >= trange(i) .AND. tmk <= trange(i+1)) EXIT
    END DO

    ! Account for reduced electron number density due to neutrals
    factor = (1.0_num - xi)**2
    IF (eos_number == EOS_IDEAL) factor = 1.0_num

    rad = factor * density**2 * psi(i) * tmk**alpha(i)
    ! Determine the radiative losses in S.I. units
    rad = rad * lr_star

    ! Convert S.I. unit radiative losses to LareXd normalised value
    rad = rad * h_star 

  END SUBROUTINE loss_function



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
