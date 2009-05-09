!-------------------------------------------------------------------------
! mass coordinate based Van Leer limited remap.
! See Bram van Leer, JCP, vol 135, p229, (1997)
!-------------------------------------------------------------------------
MODULE yremap

  USE shared_data; USE boundary

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: remap_y

  REAL(num), DIMENSION(:, :), ALLOCATABLE :: rho1, dm, cv2, flux
  REAL(num), DIMENSION(:, :), ALLOCATABLE :: dyb1, rho_v, rho_v1, dyc1

CONTAINS

  SUBROUTINE remap_y ! remap onto original Eulerian grid

    REAL(num) :: vxb, vxbm, vyb, vybm, dv

    ALLOCATE (rho1(-1:nx+2, -1:ny+2), dm(-1:nx+2, -1:ny+2), &
        cv2(-1:nx+2, -1:ny+2), flux(-1:nx+2, -2:ny+2), &
        dyb1(-1:nx+2, -1:ny+2), rho_v(-1:nx+2, -1:ny+2), &
        rho_v1(-1:nx+2, -1:ny+2), dyc1(-1:nx+1, -1:ny+1))

    dm = 0.0_num
    rho1 = rho ! store initial density in rho1

    DO iy = -1, ny+2
      DO ix = -1, nx+2
        ixm = ix - 1
        iym = iy - 1

        vxb = (vx1(ix, iy) + vx1(ix, iym)) / 2.0_num     ! vx at Ex(i, j)
        vxbm = (vx1(ixm, iy) + vx1(ixm, iym)) / 2.0_num  ! vx at Ex(i-1, j)
        vyb = (vy1(ix, iy) + vy1(ixm, iy)) / 2.0_num     ! vy at Ey(i, j)
        vybm = (vy1(ix, iym) + vy1(ixm, iym)) / 2.0_num  ! vy at Ey(i, j-1)

        dv = (REAL(xpass, num) * (vxb - vxbm) / dxb(ix) &
            + (vyb - vybm) / dyb(iy)) * dt

        ! control volume before remap
        cv1(ix, iy) = cv(ix, iy) * (1.0_num + dv)

        dv = REAL(xpass, num) * (vxb - vxbm) / dxb(ix) * dt

        ! control volume after remap
        cv2(ix, iy) = cv(ix, iy) * (1.0_num + dv)

        ! dyb before remap
        dyb1(ix, iy) = dyb(iy) + (vyb - vybm) * dt
      END DO
    END DO

    DO iy = -1, ny+1
      DO ix = -1, nx+1
        dyc1(ix, iy) = 0.5_num * (dyb1(ix, iy) + dyb1(ix, iy+1))
      END DO
    END DO

    ! Evans and Hawley (ApJ, vol 332, p650, (1988))
    ! constrained transport remap of magnetic fluxes
    CALL vy_bx_flux
    DO iy = 1, ny
      DO ix = 0, nx
        iym = iy - 1
        bx(ix, iy) = bx(ix, iy) - flux(ix, iy) + flux(ix, iym)
      END DO
    END DO

    DO iy = 0, ny
      DO ix = 1, nx
        ixm = ix - 1
        by(ix, iy) = by(ix, iy) + flux(ix, iy) - flux(ixm, iy)
      END DO
    END DO

    CALL vy_bz_flux
    DO iy = 1, ny
      DO ix = 1, nx
        iym = iy - 1
        bz(ix, iy) = bz(ix, iy) - flux(ix, iy) + flux(ix, iym)
      END DO
    END DO

    ! remap of mass + calculation of mass fluxes (dm) needed for later remaps
    CALL y_mass_flux ! calculates dm(0:nx, 0:ny+1)
    CALL dm_y_bcs    ! need dm(-1:nx+1, 0:ny+1) for velocity remap
    DO iy = 1, ny
      DO ix = 1, nx
        iym = iy - 1
        rho(ix, iy) = (rho1(ix, iy) * cv1(ix, iy) &
            + dm(ix, iym) - dm(ix, iy)) / cv2(ix, iy)
      END DO
    END DO

    ! remap specific energy density using mass coordinates
    CALL y_energy_flux
    DO iy = 1, ny
      DO ix = 1, nx
        iym = iy - 1
        energy(ix, iy) = (energy(ix, iy) * cv1(ix, iy) * rho1(ix, iy) &
            + flux(ix, iym) - flux(ix, iy)) / (cv2(ix, iy) * rho(ix, iy))
      END DO
    END DO

    ! redefine dyb1, cv1, cv2, dm and vy1 for velocity (vertex) cells
    ! in some of these calculations the flux variable is used as a
    ! temporary array
    DO iy = -1, ny+1
      DO ix = 0, nx
        ixp = ix + 1
        iyp = iy + 1

        ! vertex density before remap
        rho_v(ix, iy) = rho1(ix, iy) * cv1(ix, iy) &
            + rho1(ixp, iy ) * cv1(ixp, iy ) &
            + rho1(ix , iyp) * cv1(ix , iyp) &
            + rho1(ixp, iyp) * cv1(ixp, iyp)

        rho_v(ix, iy) = rho_v(ix, iy) / (cv1(ix, iy) + cv1(ixp, iy) &
            + cv1(ix, iyp) + cv1(ixp, iyp))
      END DO
    END DO

    DO iy = 0, ny
      DO ix = 0, nx
        ixp = ix + 1
        iyp = iy + 1
        flux(ix, iy) = cv1(ix, iy) + cv1(ixp, iy) + cv1(ix, iyp) + cv1(ixp, iyp)
      END DO
    END DO
    ! cv1 = vertex CV before remap
    cv1(0:nx, 0:ny) = flux(0:nx, 0:ny) / 4.0_num

    DO iy = 0, ny
      DO ix = 0, nx
        ixp = ix + 1
        iyp = iy + 1
        flux(ix, iy) = cv2(ix, iy) + cv2(ixp, iy) + cv2(ix, iyp) + cv2(ixp, iyp)
      END DO
    END DO
    ! cv2 = vertex CV after remap
    cv2(0:nx, 0:ny) = flux(0:nx, 0:ny) / 4.0_num

    DO iy = -2, ny+1
      DO ix = 0, nx
        iyp = iy + 1
        flux(ix, iy) = (vy1(ix, iy) + vy1(ix, iyp)) / 2.0_num
      END DO
    END DO
    ! vertex velocity used in remap
    vy1(0:nx, -2:ny+1) = flux(0:nx, -2:ny+1)

    DO iy = -1, ny+1
      DO ix = 0, nx
        iym = iy - 1
        dyb1(ix, iy) = dyc(iy) + (vy1(ix, iy) - vy1(ix, iym)) * dt
      END DO
    END DO

    DO iy = -1, ny
      DO ix = 0, nx
        ixp = ix + 1
        iyp = iy + 1
        flux(ix, iy) = dm(ix, iy) + dm(ixp, iy) + dm(ix, iyp) + dm(ixp, iyp)
      END DO
    END DO
    ! mass flux out of vertex CV
    dm(0:nx, -1:ny) = flux(0:nx, -1:ny) / 4.0_num

    DO iy = 0, ny
      DO ix = 0, nx
        iym = iy - 1
        ! vertex density after remap
        rho_v1(ix, iy) = (rho_v(ix, iy) * cv1(ix, iy) &
            + dm(ix, iym) - dm(ix, iy)) / cv2(ix, iy)
      END DO
    END DO

    CALL y_momz_flux
    DO iy = 0, ny
      DO ix = 0, nx
        iym = iy - 1
        vz(ix, iy) = (rho_v(ix, iy) * vz(ix, iy) * cv1(ix, iy) &
            + flux(ix, iym) - flux(ix, iy)) / (cv2(ix, iy) * rho_v1(ix, iy))
      END DO
    END DO

    CALL y_momx_flux
    DO iy = 0, ny
      DO ix = 0, nx
        iym = iy - 1
        vx(ix, iy) = (rho_v(ix, iy) * vx(ix, iy) * cv1(ix, iy) &
            + flux(ix, iym) - flux(ix, iy)) / (cv2(ix, iy) * rho_v1(ix, iy))
      END DO
    END DO

    CALL y_momy_flux
    DO iy = 0, ny
      DO ix = 0, nx
        iym = iy - 1
        vy(ix, iy) = (rho_v(ix, iy) * vy(ix, iy) * cv1(ix, iy) &
            + flux(ix, iym) - flux(ix, iy)) / (cv2(ix, iy) * rho_v1(ix, iy))
      END DO
    END DO

    CALL boundary_conditions

    DEALLOCATE (rho1, dm, cv2, flux, dyb1, rho_v, rho_v1, dyc1)
    ypass = 0

  END SUBROUTINE remap_y



  SUBROUTINE vy_bx_flux

    REAL(num) :: v_advect
    REAL(num) :: db, dbyp, dbyp2, dbym
    INTEGER :: iyp2

    DO iy = 0, ny
      DO ix = 0, nx
        iym  = iy - 1
        iyp  = iy + 1
        iyp2 = iy + 2
        ixp = ix + 1

        v_advect = vy1(ix, iy)

        db    = (dyb1(ix, iy  ) + dyb1(ixp, iy  )) / 2.0_num
        dbyp  = (dyb1(ix, iyp ) + dyb1(ixp, iyp )) / 2.0_num
        dbyp2 = (dyb1(ix, iyp2) + dyb1(ixp, iyp2)) / 2.0_num
        dbym  = (dyb1(ix, iym ) + dyb1(ixp, iym )) / 2.0_num

        w4 = bx(ix, iy ) / db
        w5 = bx(ix, iyp) / dbyp

        flux(ix, iy) = (MAX(0.0_num, v_advect) * w4 &
            + MIN(0.0_num, v_advect) * w5) * dt

        w1 = bx(ix, iyp ) / dbyp  - bx(ix, iy ) / db
        w2 = bx(ix, iy  ) / db    - bx(ix, iym) / dbym
        w3 = bx(ix, iyp2) / dbyp2 - bx(ix, iyp) / dbyp

        IF (v_advect > 0.0) THEN
          w5 = ABS(v_advect) * dt / db
          w4 = (2.0_num - w5) * ABS(w1) / dyc1(ix, iy) &
              + (1.0_num + w5) * ABS(w2) / dyc1(ix, iym)
          w4 = w4 / 6.0_num
          w8 = 0.5_num * (SIGN(1.0_num, w1) + SIGN(1.0_num, w2))
          w6 = w8 * MIN(ABS(w4) * dyb1(ix, iy), ABS(w1), ABS(w2))
          flux(ix, iy) = flux(ix, iy) &
              + bzone(ix, iy) * v_advect * dt * w6 * (1.0_num - w5)
        ELSE
          w5 = ABS(v_advect) * dt / dbyp
          w4 = (2.0_num - w5) * ABS(w1) / dyc1(ix, iy) &
              + (1.0_num + w5) * ABS(w3) / dyc1(ix, iyp)
          w4 = w4 / 6.0_num
          w8 = 0.5_num * (SIGN(1.0_num, w1) + SIGN(1.0_num, w3))
          w6 = -w8 * MIN(ABS(w4) * dyb1(ix, iyp), ABS(w1), ABS(w3))
          flux(ix, iy) = flux(ix, iy) &
              + bzone(ix, iy) * v_advect * dt * w6 * (1.0_num - w5)
        END IF

      END DO
    END DO

  END SUBROUTINE vy_bx_flux



  SUBROUTINE vy_bz_flux

    REAL(num) :: v_advect
    INTEGER :: iyp2

    DO iy = 0, ny
      DO ix = 1, nx
        iym  = iy - 1
        iyp  = iy + 1
        iyp2 = iy + 2
        ixm = ix - 1

        v_advect = (vy1(ix, iy) + vy1(ixm, iy)) / 2.0_num

        w4 = bz(ix, iy ) / dyb1(ix, iy )
        w5 = bz(ix, iyp) / dyb1(ix, iyp)

        flux(ix, iy) = (MAX(0.0_num, v_advect) * w4 &
            + MIN(0.0_num, v_advect) * w5) * dt

        w1 = bz(ix, iyp ) / dyb1(ix, iyp ) - bz(ix, iy ) / dyb1(ix, iy )
        w2 = bz(ix, iy  ) / dyb1(ix, iy  ) - bz(ix, iym) / dyb1(ix, iym)
        w3 = bz(ix, iyp2) / dyb1(ix, iyp2) - bz(ix, iyp) / dyb1(ix, iyp)

        IF (v_advect > 0.0) THEN
          w5 = ABS(v_advect) * dt / dyb1(ix, iy)
          w4 = (2.0_num - w5) * ABS(w1) / dyc1(ix, iy) &
              + (1.0_num + w5) * ABS(w2) / dyc1(ix, iym)
          w4 = w4 / 6.0_num
          w8 = 0.5_num * (SIGN(1.0_num, w1) + SIGN(1.0_num, w2))
          w6 = w8 * MIN(ABS(w4) * dyb1(ix, iy), ABS(w1), ABS(w2))
          flux(ix, iy) = flux(ix, iy) &
              + bzone(ix, iy) * v_advect * dt * w6 * (1.0_num - w5)
        ELSE
          w5 = ABS(v_advect) * dt / dyb1(ix, iyp)
          w4 = (2.0_num - w5) * ABS(w1) / dyc1(ix, iy) &
              + (1.0_num + w5) * ABS(w3) / dyc1(ix, iyp)
          w4 = w4 / 6.0_num
          w8 = 0.5_num * (SIGN(1.0_num, w1) + SIGN(1.0_num, w3))
          w6 = -w8 * MIN(ABS(w4) * dyb1(ix, iyp), ABS(w1), ABS(w3))
          flux(ix, iy) = flux(ix, iy) &
              + bzone(ix, iy) * v_advect * dt * w6 * (1.0_num - w5)
        END IF

      END DO
    END DO

  END SUBROUTINE vy_bz_flux



  SUBROUTINE y_mass_flux

    REAL(num) :: v_advect, flux_rho
    INTEGER :: iyp2

    DO iy = 0, ny
      DO ix = 0, nx+1
        iym  = iy - 1
        iyp  = iy + 1
        iyp2 = iy + 2
        ixm = ix - 1

        v_advect = (vy1(ix, iy) + vy1(ixm, iy)) / 2.0_num

        dm(ix, iy) = (MAX(0.0_num, v_advect) * rho(ix, iy) &
            + MIN(0.0_num, v_advect) * rho(ix, iyp)) * dt

        w1 = rho(ix, iyp ) - rho(ix, iy )
        w2 = rho(ix, iy  ) - rho(ix, iym)
        w3 = rho(ix, iyp2) - rho(ix, iyp)

        IF (v_advect > 0.0) THEN
          w5 = ABS(v_advect) * dt / dyb1(ix, iy)
          w4 = (2.0_num - w5) * ABS(w1) / dyc1(ix, iy) &
              + (1.0_num + w5) * ABS(w2) / dyc1(ix, iym)
          w4 = w4 / 6.0_num
          w8 = 0.5_num * (SIGN(1.0_num, w1) + SIGN(1.0_num, w2))
          w6 = w8 * MIN(ABS(w4) * dyb1(ix, iy), ABS(w1), ABS(w2))
          flux_rho = v_advect * dt * w6 * (1.0_num - w5)
        ELSE
          w5 = ABS(v_advect) * dt / dyb1(ix, iyp)
          w4 = (2.0_num - w5) * ABS(w1) / dyc1(ix, iy) &
              + (1.0_num + w5) * ABS(w3) / dyc1(ix, iyp)
          w4 = w4 / 6.0_num
          w8 = 0.5_num * (SIGN(1.0_num, w1) + SIGN(1.0_num, w3))
          w6 = -w8 * MIN(ABS(w4) * dyb1(ix, iyp), ABS(w1), ABS(w3))
          flux_rho = v_advect * dt * w6 * (1.0_num - w5)
        END IF
        dm(ix, iy) = (bzone(ix, iy) * flux_rho + dm(ix, iy)) * dxb(ix)
      END DO
    END DO

  END SUBROUTINE y_mass_flux



  SUBROUTINE y_energy_flux ! energy remap in mass coordinates

    REAL(num) :: v_advect
    INTEGER :: iyp2

    DO iy = 0, ny
      DO ix = 0, nx
        iym  = iy - 1
        iyp  = iy + 1
        iyp2 = iy + 2
        ixm = ix - 1

        v_advect = (vy1(ix, iy) + vy1(ixm, iy)) / 2.0_num

        w1 = energy(ix, iyp ) - energy(ix, iy )
        w2 = energy(ix, iy  ) - energy(ix, iym)
        w3 = energy(ix, iyp2) - energy(ix, iyp)

        IF (v_advect > 0.0) THEN
          w5 = ABS(v_advect) * dt / dyb1(ix, iy)
          w7 = energy(ix, iy)
          w6 = ABS(dm(ix, iy)) / (rho1(ix, iy) * dyb1(ix, iy) * dxb(ix))
          w4 = (2.0_num - w5) * ABS(w1) / dyc(iy) &
              + (1.0_num + w5) * ABS(w2) / dyc(iym)
          w4 = w4 / 6.0_num
          w8 = 0.5_num * (SIGN(1.0_num, w1) + SIGN(1.0_num, w2))
          w5 = w8 * MIN(ABS(w4) * dyb(iy), ABS(w1), ABS(w2))
        ELSE
          w5 = ABS(v_advect) * dt / dyb1(ix, iyp)
          w7 = energy(ix, iyp)
          w6 = ABS(dm(ix, iy)) / (rho1(ix, iyp) * dyb1(ix, iyp) * dxb(ix))
          w4 = (2.0_num - w5) * ABS(w1) / dyc(iy) &
              + (1.0_num + w5) * ABS(w3) / dyc(iyp)
          w4 = w4 / 6.0_num
          w8 = 0.5_num * (SIGN(1.0_num, w1) + SIGN(1.0_num, w3))
          w5 = -w8 * MIN(ABS(w4) * dyb(iyp), ABS(w1), ABS(w3))
        END IF

        flux(ix, iy) = dm(ix, iy) * (w7 + bzone(ix, iy) * w5 * (1.0_num - w6))
      END DO
    END DO

  END SUBROUTINE y_energy_flux



  SUBROUTINE y_momx_flux ! energy remap in mass coordinates

    REAL(num) :: v_advect, m, mp, ai, aip, dk
    INTEGER :: iyp2

    DO iy = -1, ny
      DO ix = 0, nx
        iym  = iy - 1
        iyp  = iy + 1
        iyp2 = iy + 2

        v_advect = vy1(ix, iy)

        w1 = vx(ix, iyp ) - vx(ix, iy )
        w2 = vx(ix, iy  ) - vx(ix, iym)
        w3 = vx(ix, iyp2) - vx(ix, iyp)

        IF (v_advect > 0.0) THEN
          w5 = ABS(v_advect) * dt / dyb1(ix, iy)
          w7 = vx(ix, iy)
          w6 = ABS(dm(ix, iy)) / (rho_v(ix, iy) * dyb1(ix, iy) * dxc(ix))
          w4 = (2.0_num - w5) * ABS(w1) / dyb(iyp) &
              + (1.0_num + w5) * ABS(w2) / dyb(iy)
          w4 = w4 / 6.0_num
          w8 = 0.5_num * (SIGN(1.0_num, w1) + SIGN(1.0_num, w2))
          w5 = w8 * MIN(ABS(w4) * dyc(iy), ABS(w1), ABS(w2))
        ELSE
          w5 = ABS(v_advect) * dt / dyb1(ix, iyp)
          w7 = vx(ix, iyp)
          w6 = ABS(dm(ix, iy)) / (rho_v(ix, iyp) * dyb1(ix, iyp) * dxc(ix))
          w4 = (2.0_num - w5) * ABS(w1) / dyb(iyp) &
              + (1.0_num + w5) * ABS(w3) / dyb(iyp2)
          w4 = w4 / 6.0_num
          w8 = 0.5_num * (SIGN(1.0_num, w1) + SIGN(1.0_num, w3))
          w5 = -w8 * MIN(ABS(w4) * dyc(iyp), ABS(w1), ABS(w3))
        END IF

        flux(ix, iy) = w7 + bzone(ix, iy) * w5 * (1.0_num - w6)
      END DO
    END DO

    IF (rke) THEN
      DO iy = 0, ny-1
        DO ix = 0, nx
          iym = iy - 1
          iyp = iy + 1
          ixp = ix + 1

          m = rho_v1(ix, iy) * cv2(ix, iy)
          mp = rho_v1(ix, iyp) * cv2(ix, iyp)

          ai = (vx(ix, iy) - flux(ix, iym)) * dm(ix, iym) / m &
              + (flux(ix, iy) - vx(ix, iy)) * dm(ix, iy) / m

          aip = (vx(ix, iyp) - flux(ix, iy)) * dm(ix, iy) / mp &
              + (flux(ix, iyp) - vx(ix, iyp)) * dm(ix, iyp) / mp

          dk = (vx(ix, iyp) - vx(ix, iy)) * (flux(ix, iy) &
              - 0.5_num * (vx(ix, iyp) + vx(ix, iy))) &
              + 0.5_num * ai * (flux(ix, iy) - vx(ix, iy)) &
              + 0.5_num * aip * (vx(ix, iyp) - flux(ix, iy))

          dk = dk * dm(ix, iy) / 2.0_num
          delta_ke(ix , iyp) = delta_ke(ix , iyp) + dk
          delta_ke(ixp, iyp) = delta_ke(ixp, iyp) + dk
        END DO
      END DO
    END IF

    flux(0:nx, -1:ny) = flux(0:nx, -1:ny) * dm(0:nx, -1:ny)

  END SUBROUTINE y_momx_flux



  SUBROUTINE y_momz_flux ! energy remap in mass coordinates

    REAL(num) :: v_advect, m, mp, ai, aip, dk
    INTEGER :: iyp2

    DO iy = -1, ny
      DO ix = 0, nx
        iym  = iy - 1
        iyp  = iy + 1
        iyp2 = iy + 2

        v_advect = vy1(ix, iy)

        w1 = vz(ix, iyp ) - vz(ix, iy )
        w2 = vz(ix, iy  ) - vz(ix, iym)
        w3 = vz(ix, iyp2) - vz(ix, iyp)

        IF (v_advect > 0.0) THEN
          w5 = ABS(v_advect) * dt / dyb1(ix, iy)
          w7 = vz(ix, iy)
          w6 = ABS(dm(ix, iy)) / (rho_v(ix, iy) * dyb1(ix, iy) * dxc(ix))
          w4 = (2.0_num - w5) * ABS(w1) / dyb(iyp) &
              + (1.0_num + w5) * ABS(w2) / dyb(iy)
          w4 = w4 / 6.0_num
          w8 = 0.5_num * (SIGN(1.0_num, w1) + SIGN(1.0_num, w2))
          w5 = w8 * MIN(ABS(w4) * dyc(iy), ABS(w1), ABS(w2))
        ELSE
          w5 = ABS(v_advect) * dt / dyb1(ix, iyp)
          w7 = vz(ix, iyp)
          w6 = ABS(dm(ix, iy)) / (rho_v(ix, iyp) * dyb1(ix, iyp) * dxc(ix))
          w4 = (2.0_num - w5) * ABS(w1) / dyb(iyp) &
              + (1.0_num + w5) * ABS(w3) / dyb(iyp2)
          w4 = w4 / 6.0_num
          w8 = 0.5_num * (SIGN(1.0_num, w1) + SIGN(1.0_num, w3))
          w5 = -w8 * MIN(ABS(w4) * dyc(iyp), ABS(w1), ABS(w3))
        END IF

        flux(ix, iy) = w7 + bzone(ix, iy) * w5 * (1.0_num - w6)
      END DO
    END DO

    IF (rke) THEN
      DO iy = 0, ny-1
        DO ix = 0, nx
          iym = iy - 1
          iyp = iy + 1
          ixp = ix + 1

          m = rho_v1(ix, iy) * cv2(ix, iy)
          mp = rho_v1(ix, iyp) * cv2(ix, iyp)

          ai = (vz(ix, iy) - flux(ix, iym)) * dm(ix, iym) / m &
              + (flux(ix, iy) - vz(ix, iy)) * dm(ix, iy) / m

          aip = (vz(ix, iyp) - flux(ix, iy)) * dm(ix, iy) / mp &
              + (flux(ix, iyp) - vz(ix, iyp)) * dm(ix, iyp) / mp

          dk = (vz(ix, iyp) - vz(ix, iy)) * (flux(ix, iy) &
              - 0.5_num * (vz(ix, iyp) + vz(ix, iy))) &
              + 0.5_num * ai * (flux(ix, iy) - vz(ix, iy)) &
              + 0.5_num * aip * (vz(ix, iyp) - flux(ix, iy))

          dk = dk * dm(ix, iy) / 2.0_num
          delta_ke(ix , iyp) = delta_ke(ix , iyp) + dk
          delta_ke(ixp, iyp) = delta_ke(ixp, iyp) + dk
        END DO
      END DO
    END IF

    flux(0:nx, -1:ny) = flux(0:nx, -1:ny) * dm(0:nx, -1:ny)

  END SUBROUTINE y_momz_flux



  SUBROUTINE y_momy_flux

    REAL(num) :: v_advect, m, mp, ai, aip, dk
    INTEGER :: iyp2

    DO iy = -1, ny
      DO ix = 0, nx
        iym  = iy - 1
        iyp  = iy + 1
        iyp2 = iy + 2

        v_advect = vy1(ix, iy)

        w1 = vy(ix, iyp ) - vy(ix, iy )
        w2 = vy(ix, iy  ) - vy(ix, iym)
        w3 = vy(ix, iyp2) - vy(ix, iyp)

        IF (v_advect > 0.0) THEN
          w5 = ABS(v_advect) * dt / dyb1(ix, iy)
          w7 = vy(ix, iy)
          w6 = ABS(dm(ix, iy)) / (rho_v(ix, iy) * dyb1(ix, iy) * dxc(ix))
          w4 = (2.0_num - w5) * ABS(w1) / dyb(iyp)&
              + (1.0_num + w5) * ABS(w2) / dyb(iy)
          w4 = w4 / 6.0_num
          w8 = 0.5_num * (SIGN(1.0_num, w1) + SIGN(1.0_num, w2))
          w5 = w8 * MIN(ABS(w4) * dyc(iy), ABS(w1), ABS(w2))
        ELSE
          w5 = ABS(v_advect) * dt / dyb1(ix, iyp)
          w7 = vy(ix, iyp)
          w6 = ABS(dm(ix, iy)) / (rho_v(ix, iyp) * dyb1(ix, iyp) * dxc(ix))
          w4 = (2.0_num - w5) * ABS(w1) / dyb(iyp) &
              + (1.0_num + w5) * ABS(w3) / dyb(iyp2)
          w4 = w4 / 6.0_num
          w8 = 0.5_num * (SIGN(1.0_num, w1) + SIGN(1.0_num, w3))
          w5 = -w8 * MIN(ABS(w4) * dyc(iyp), ABS(w1), ABS(w3))
        END IF

        flux(ix, iy) = w7 + bzone(ix, iy) * w5 * (1.0_num - w6)
      END DO
    END DO

    IF (rke) THEN
      DO iy = 0, ny-1
        DO ix = 0, nx
          iym = iy - 1
          iyp = iy + 1
          ixp = ix + 1

          m = rho_v1(ix, iy) * cv2(ix, iy)
          mp = rho_v1(ix, iyp) * cv2(ix, iyp)

          ai = (vy(ix, iy) - flux(ix, iym)) * dm(ix, iym) / m &
              + (flux(ix, iy) - vy(ix, iy)) * dm(ix, iy) / m

          aip = (vy(ix, iyp) - flux(ix, iy)) * dm(ix, iy) / mp &
              + (flux(ix, iyp) - vy(ix, iyp)) * dm(ix, iyp) / mp

          dk = (vy(ix, iyp) - vy(ix, iy)) * (flux(ix, iy) &
              - 0.5_num * (vy(ix, iyp) + vy(ix, iy))) &
              + 0.5_num * ai * (flux(ix, iy) - vy(ix, iy)) &
              + 0.5_num * aip * (vy(ix, iyp) - flux(ix, iy))

          dk = dk * dm(ix, iy) / 2.0_num
          delta_ke(ix , iyp) = delta_ke(ix , iyp) + dk
          delta_ke(ixp, iyp) = delta_ke(ixp, iyp) + dk
        END DO
      END DO
    END IF

    flux(0:nx, -1:ny) = flux(0:nx, -1:ny) * dm(0:nx, -1:ny)

  END SUBROUTINE y_momy_flux



  SUBROUTINE dm_y_bcs

    CALL MPI_SENDRECV(dm(0:nx+1, 1), nx+2, mpireal, &
        down, tag, dm(0:nx+1, ny+1), nx+2, mpireal, &
        up, tag, comm, status, errcode)

    IF (up == MPI_PROC_NULL) &
        dm(0:nx+1, ny+1) = dm(0:nx+1, ny)

    CALL MPI_SENDRECV(dm(0:nx+1, ny-1), nx+2, mpireal, &
        up, tag, dm(0:nx+1, -1), nx+2, mpireal, &
        down, tag, comm, status, errcode)

    IF (down == MPI_PROC_NULL) &
        dm(0:nx+1, -1) = dm(0:nx+1, 0)

  END SUBROUTINE dm_y_bcs

END MODULE yremap
