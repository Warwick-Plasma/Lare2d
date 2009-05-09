! All the subroutines in this module are for the partially ionised flux
! emergence simulations; see Leake & Arber, 2006

MODULE neutral

  USE shared_data
  USE boundary
  USE eos

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: perpendicular_resistivity, newton_relax, &
      neutral_fraction, setup_neutral, get_neutral

CONTAINS

  SUBROUTINE setup_neutral
    ! Ion neutral collision cross section(m^2)
    REAL(num) :: Sigma_in = 5.0e-19_num
    REAL(num) :: Tr

    IF (include_neutrals) THEN
      ALLOCATE(xi_n(-1:nx+2, -1:ny+2))
      xi_n = 0.0_num
    END IF

    IF (cowling_resistivity) THEN
      ALLOCATE(eta_perp(-1:nx+2, -1:ny+2))
      ALLOCATE(parallel_current(0:nx, 0:ny))
      ALLOCATE(perp_current(0:nx, 0:ny))
    END IF

    ! Ionisation potential of hydrogen(J)
    ionise_pot = ionise_pot_0

    ! Temperature of the photospheric radiation field
    Tr = 7230.85_num
    Tr_bar = 1.0_num / Tr

    ! Calculate fbar^(2 / 3) in (k^-1 m^-2)
    f_bar = (pi * me_0 * kb_0) / h_0**2
    f_bar = f_bar**(3.0_num / 2.0_num)

    ! Calculate tbar in (K)
    t_bar = ionise_pot / kb

    ! Calculate rbar in (kg^-1)
    r_bar = 4.0_num / MBAR

    ! Calculate eta_bar in (m^4 / (k s kg^2))
    eta_bar = (0.5_num / MBAR &
        * SQRT(16.0_num * kb_0 / (pi * MBAR)) * sigma_in)**(-1)

  END SUBROUTINE setup_neutral



  SUBROUTINE perpendicular_resistivity

    ! This subroutine calculates the cross field resistivity at the current
    ! temperature. If you're not using the "neutral_fraction" subroutine and
    ! the "Saha" equation of state, then this routine will give VERY strange
    ! results which are probably meaningless.

    ! normalising values are L0 = 150km, v0 = 6.5km / s, rho0 = 2.7e-4 kg / m3
    ! t0 = 23s, T0 = 6420K, P0 = 1.2e4 Pa, B0 = 1200G (0.12T)

    REAL(num) :: f, xi_v, bxv, byv, bzv, bfieldsq, rho_v, T_v, T
    INTEGER :: ixp, iyp

    DO iy = 0, ny
      DO ix = 0, nx
        ixp = ix + 1
        iyp = iy + 1

        ! Get the vertex density
        rho_v = rho(ix, iy) * cv(ix, iy) + rho(ixp, iy) * cv(ixp, iy) &
            + rho(ix, iyp) * cv(ix, iyp) + rho(ixp, iyp) * cv(ixp, iyp)
        rho_v = rho_v / (cv(ix, iy) + cv(ixp, iy) + cv(ix, iyp) + cv(ixp, iyp))

        ! Get the vertex magnetic field
        bxv = (bx(ix, iy) + bx(ix, iyp)) / 2.0_num
        byv = (by(ix, iy) + by(ixp, iy)) / 2.0_num
        bzv = (bz(ix, iy) + bz(ixp, iy) + bz(ix, iyp) &
            + bz(ixp, iyp)) / 4.0_num
        bfieldsq = bxv**2 + byv**2 + bzv**2

        ! Get the vertex temperature
        CALL Get_Temp(rho(ix, iy), energy(ix, iy), eos_number, ix, iy, T_v)
        CALL Get_Temp(rho(ixp, iy), energy(ixp, iy), eos_number, ix, iy, T)

        T_v = T_v + T
        CALL Get_Temp(rho(ix, iyp), energy(ix, iyp), eos_number, ix, iy, T)

        T_v = T_v + T
        CALL Get_Temp(rho(ixp, iyp), energy(ixp, iyp), eos_number, ix, iy, T)

        T_v = T_v + T
        T_v = T_v / 4.0_num

        xi_v = Get_Neutral(T_v, rho_v)

        f = MAX(1.0_num - xi_v, none_zero)
        IF (f .GT. 0) THEN
          eta_perp(ix, iy) = eta_bar * xi_v / f * bfieldsq &
              / rho_v**2 / SQRT(T_v)
        ELSE
          eta_perp(ix, iy) = 0.0_num
        END IF

        eta_perp(ix, iy) = MIN(eta_perp(ix, iy), 4.0_num)
      END DO
    END DO

  END SUBROUTINE perpendicular_resistivity



  FUNCTION Get_Neutral(T_v, rho_v)

    REAL(num), INTENT(IN) :: T_V, rho_v
    REAL(num) :: Get_Neutral
    REAL(num) :: f, b, r

    f = f_bar * SQRT(T_v) * EXP(-T_bar / T_v) ! T from b has been cancelled
    b = Tr_bar * EXP(0.25_num * T_bar / T_v * (Tr_bar * T_v - 1.0_num))
    r = 0.5_num * (-1.0_num + SQRT(1.0_num + r_bar * rho_v * b / f))
    Get_Neutral = r / (1.0_num + r)

  END FUNCTION  Get_Neutral



  SUBROUTINE neutral_fraction(material)

    INTEGER, INTENT(IN) :: material
    REAL(num) :: bof, r, T, rho0, e0, dx, x
    REAL(num), DIMENSION(2) :: Ta, fa, xi_a
    REAL(num) :: ionise_pot_local
    INTEGER :: loop

    IF (material .EQ. EOS_ION) THEN
      ionise_pot_local = ionise_pot
    ELSE
      ionise_pot_local = 0.0_num
    END IF

    ! Variable bof is b / f in the original version
    DO iy = -1, ny + 2
      DO ix = -1, nx + 2
        rho0 = rho(ix, iy)
        e0 = energy(ix, iy)
        Ta = (gamma - 1.0_num) &
            * (/ MAX((e0 - ionise_pot_local) / 2.0_num, none_zero), e0 /)

        IF (Ta(1) > Ta(2)) THEN
          PRINT * , "Temperature bounds problem", Ta
          STOP
        END IF

        dx = Ta(2) - Ta(1)
        T = Ta(1)

        DO loop = 1, 100
          dx = dx / 2.0_num
          x = T  + dx
          bof = Tr_bar / (f_bar * SQRT(x)) &
              * EXP((0.25_num * (Tr_bar * x - 1.0_num) + 1.0_num) * T_bar / x)
          r = 0.5_num * (-1.0_num + SQRT(1.0_num + r_bar * rho0 * bof))
          xi_a(1) = r / (1.0_num + r)
          fa(1) = x - (gamma - 1.0_num) * (e0 &
              - (1.0_num - xi_a(1)) * ionise_pot_local) / (2.0_num - xi_a(1))
          IF (fa(1) <= 0.0_num) T = x
          IF (ABS(dx) < 1.e-8_num .OR. fa(1) == 0.0_num) EXIT
        END DO

        bof = Tr_bar / (f_bar * SQRT(T)) &
            * EXP((0.25_num * (Tr_bar * T - 1.0_num) + 1.0_num) * T_bar / T)
        r = 0.5_num * (-1.0_num + SQRT(1.0_num + r_bar * rho0 * bof))
        xi_n(ix, iy) = r / (1.0_num + r)
      END DO
    END DO

  END SUBROUTINE neutral_fraction



  SUBROUTINE newton_relax

!!$    INTEGER, DIMENSION(1) :: ref_index, z0(1) = 1
!!$    LOGICAL :: first_call = .TRUE., run_loop = .TRUE.
!!$
!!$    ! This should only be run above the photosphere so first call sets up
!!$    ! the lowest value of iz to use if at all, the -2 is due to zc starting
!!$    ! at -1
!!$    IF (first_call) THEN
!!$      z0 = MINLOC(ABS(zc - 0.0_num)) - 2
!!$      ! This process doesn't have any cells in the corona
!!$      IF (z0(1) > nz) run_loop = .FALSE.
!!$      IF (z0(1) < 1) z0(1) = 1 ! only need to run over the internal domain
!!$      first_call = .FALSE.
!!$    END IF
!!$
!!$    ! For every point need to find the reference density value and hence
!!$    ! the tau and temperature
!!$    IF (run_loop) THEN
!!$      DO iz = z0(1), nz
!!$        DO iy = 1, ny
!!$          DO ix = 1, nx
!!$            ! the 2 is subtracted due to rho_ref starting at -1
!!$            ref_index = MINLOC(ABS(rho(ix, iy) - rho_ref)) - 2
!!$            energy(ix, iy) = (energy(ix, iy) + dt / tau_ref(ref_index(1)) * &
!!$                T_ref(ref_index(1)) / (gamma - 1.0_num)) &
!!$                / (1.0_num + dt / tau_ref(ref_index(1)))
!!$          END DO
!!$        END DO
!!$      END DO
!!$    END IF

    CALL energy_bcs

  END SUBROUTINE newton_relax

END MODULE neutral
