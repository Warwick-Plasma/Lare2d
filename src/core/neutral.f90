! All the subroutines in this module are for the partially ionised flux
! emergence simulations; see Leake & Arber, 2006

MODULE neutral

  USE shared_data
  USE boundary

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: perpendicular_resistivity, neutral_fraction, setup_neutral, get_neutral

CONTAINS

  SUBROUTINE setup_neutral
    ! Ion neutral collision cross section(m^2)
    REAL(num) :: sigma_in = 5.0e-19_num
    REAL(num) :: mbar

    ALLOCATE(xi_n(-1:nx+2, -1:ny+2))
    xi_n = 0.0_num       

    IF (cowling_resistivity) THEN
      ALLOCATE(eta_perp(-1:nx+2, -1:ny+2))
      ALLOCATE(parallel_current(0:nx, 0:ny))
      ALLOCATE(perp_current(0:nx, 0:ny))
    END IF
    
    ! Temperature of the photospheric radiation field
    tr = 6420.0_num
    
    ! Calculate fbar^(2 / 3) in (k^-1 m^-2)
<<<<<<< HEAD:src/core/neutral.f90
    f_bar = 2.0_num * pi * (me_si / hp_si) * (kb_si / hp_si)
    f_bar = f_bar**(3.0_num / 2.0_num)
=======
    f_bar = pi * (me_si / hp_si) * (kb_si / hp_si)
    f_bar = SQRT(2.0_num) * f_bar**(3.0_num / 2.0_num)
>>>>>>> Updated to match local version:src/core/neutral.f90
    
    ! Calculate tbar in (K)
    t_bar = ionise_pot_si / kb_si
    
    ! Calculate rbar in (kg^-1)
    mbar = mh_si * mf                     
    r_bar = 4.0_num / mbar
    
    ! Calculate eta_bar in (m^4 / (k s kg^2))
    eta_bar = 2.0_num * mbar &
        / (SQRT(16.0_num * kb_si / (pi * mbar)) * sigma_in)
    
  END SUBROUTINE setup_neutral


                                        
  SUBROUTINE perpendicular_resistivity

    ! This subroutine calculates the cross field resistivity at the current
    ! temperature. 

    REAL(num) :: f, xi_v, bxv, byv, bzv, bfieldsq, rho_v, t_v, T
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
        t_v = (gamma - 1.0_num) &
            * (energy(ix,iy) - (1.0_num - xi_n(ix, iy)) * ionise_pot) &
            / ((2.0_num - xi_n(ix, iy)))
        T = (gamma - 1.0_num) &
            * (energy(ixp,iy) - (1.0_num - xi_n(ixp, iy)) * ionise_pot) &
            / ((2.0_num - xi_n(ixp, iy)))
        t_v = t_v + T

        T = (gamma - 1.0_num) &
            * (energy(ix,iyp) - (1.0_num - xi_n(ix, iyp)) * ionise_pot) &
            / ((2.0_num - xi_n(ix, iyp)))
        t_v = t_v + T

        T = (gamma - 1.0_num) &
            * (energy(ixp,iyp) - (1.0_num - xi_n(ixp, iyp)) * ionise_pot) &
            / ((2.0_num - xi_n(ixp, iyp)))
        t_v = t_v + T

        t_v = t_v / 4.0_num

        xi_v = get_neutral(t_v, rho_v, yb(iy))

        f = MAX(1.0_num - xi_v, none_zero)
        IF (f .GT. 0) THEN
          eta_perp(ix, iy) = eta_bar * xi_v / f * bfieldsq &
              / rho_v**2 / SQRT(t_v)
        ELSE
          eta_perp(ix, iy) = 0.0_num
        END IF

!        eta_perp(ix, iy) = MIN(eta_perp(ix, iy), 4.0_num)
      END DO
    END DO

  END SUBROUTINE perpendicular_resistivity



  FUNCTION get_neutral(t_v, rho_v, height)
    
    REAL(num), INTENT(IN) :: t_v, rho_v, height
    REAL(num) :: get_neutral
    REAL(num) :: bof, r, t_rad, dilution
    
    t_rad = tr
    dilution = 0.5_num
    IF (height <= 0.0_num) THEN
      t_rad = t_v
      dilution = 1.0_num  
    END IF
      
    bof = 1.0_num / (dilution * f_bar * t_rad * SQRT(t_v)) &
        * EXP((0.25_num * (t_v / t_rad - 1.0_num) + 1.0_num) &
        * T_bar / t_v)
    r = 0.5_num * (-1.0_num + SQRT(1.0_num + r_bar * rho_v * bof))
    get_neutral = r / (1.0_num + r)
    
  END FUNCTION  get_neutral



  SUBROUTINE neutral_fraction

    REAL(num) :: bof, r, T, rho0, e0, dx, x
    REAL(num), DIMENSION(2) :: ta, fa, xi_a
    REAL(num) :: ionise_pot_local
    INTEGER :: loop

    ! Variable bof is b / f in the original version
    DO iy = -1, ny + 2
      DO ix = -1, nx + 2
        rho0 = rho(ix, iy)
        e0 = energy(ix, iy)
        ta = (gamma - 1.0_num) &
            * (/ MAX((e0 - ionise_pot_local) / 2.0_num, none_zero), e0 /)

        IF (ta(1) > ta(2)) THEN
          PRINT * , "Temperature bounds problem", ta
          STOP
        END IF

        dx = ta(2) - ta(1)
        t = ta(1)

        DO loop = 1, 100
          dx = dx / 2.0_num
          x = t  + dx
          xi_a(1) = get_neutral(x, rho0, yb(iy))   
          fa(1) = x - (gamma - 1.0_num) * (e0 &
              - (1.0_num - xi_a(1)) * ionise_pot_local) / (2.0_num - xi_a(1))
          IF (fa(1) <= 0.0_num) t = x
          IF (ABS(dx) < 1.e-8_num .OR. fa(1) == 0.0_num) EXIT
        END DO

        xi_n(ix, iy) = get_neutral(x, rho0, yb(iy))   
      END DO
    END DO

  END SUBROUTINE neutral_fraction



!   SUBROUTINE newton_relax
! 
! !this routine currently never used
! 
!     INTEGER, DIMENSION(1) :: ref_index, z0(1) = 1
!     LOGICAL :: first_call = .TRUE., run_loop = .TRUE.
! 
!     ! This should only be run above the photosphere so first call sets up
!     ! the lowest value of iz to use if at all, the -2 is due to zc starting
!     ! at -1
!     IF (first_call) THEN
!       z0 = MINLOC(ABS(zc - 0.0_num)) - 2
!       ! This process doesn't have any cells in the corona
!       IF (z0(1) > nz) run_loop = .FALSE.
!       IF (z0(1) < 1) z0(1) = 1 ! only need to run over the internal domain
!       first_call = .FALSE.
!     END IF
! 
!     ! For every point need to find the reference density value and hence
!     ! the tau and temperature
!     IF (run_loop) THEN
!       DO iz = z0(1), nz
!         DO iy = 1, ny
!           DO ix = 1, nx
!             ! the 2 is subtracted due to rho_ref starting at -1
!             ref_index = MINLOC(ABS(rho(ix, iy) - rho_ref)) - 2
!             energy(ix, iy) = (energy(ix, iy) + dt / tau_ref(ref_index(1)) * &
!                 T_ref(ref_index(1)) / (gamma - 1.0_num)) &
!                 / (1.0_num + dt / tau_ref(ref_index(1)))
!           END DO
!         END DO
!       END DO
!     END IF
! 
!     CALL energy_bcs   
! 
!   END SUBROUTINE newton_relax
!                

END MODULE neutral
