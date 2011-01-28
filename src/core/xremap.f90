!-------------------------------------------------------------------------
! mass coordinate based Van Leer limited remap.
! See Bram van Leer, JCP, vol 135, p229, (1997)
!-------------------------------------------------------------------------
MODULE xremap

  USE shared_data; USE boundary

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: remap_x

  REAL(num), DIMENSION(:, :), ALLOCATABLE :: dxb1, rho_v, rho_v1
  REAL(num), DIMENSION(:, :), ALLOCATABLE :: rho1, dm, cv2, flux, dxc1

CONTAINS

  SUBROUTINE remap_x ! remap onto original Eulerian grid

    REAL(num) :: vxb, vxbm, vyb, vybm, dv

    ALLOCATE (rho1(-1:nx+2, -1:ny+2), dm(-1:nx+2, -1:ny+2), &
        cv2(-1:nx+2, -1:ny+2), flux(-2:nx+2, -1:ny+2), &
        dxb1(-1:nx+2, -1:ny+2), rho_v(-1:nx+2, -1:ny+2), &
        rho_v1(-1:nx+2, -1:ny+2), dxc1(-1:nx+1, -1:ny+1))

    dm = 0.0_num
    rho1 = rho ! store initial density in rho1

    DO iy = -1, ny+2
      iym = iy - 1
      DO ix = -1, nx+2
        ixm = ix - 1

        vxb = (vx1(ix, iy) + vx1(ix, iym)) * 0.5_num     ! vx at Ex(i, j)
        vxbm = (vx1(ixm, iy) + vx1(ixm, iym)) * 0.5_num  ! vx at Ex(i-1, j)
        vyb = (vy1(ix, iy) + vy1(ixm, iy)) * 0.5_num     ! vy at Ey(i, j)
        vybm = (vy1(ix, iym) + vy1(ixm, iym)) * 0.5_num  ! vy at Ey(i, j-1)

        dv = (REAL(ypass, num) * (vyb - vybm) / dyb(iy) &
            + (vxb - vxbm) / dxb(ix)) * dt

        ! control volume before remap
        cv1(ix, iy) = cv(ix, iy) * (1.0_num + dv)

        dv = REAL(ypass, num) * (vyb - vybm) / dyb(iy) * dt

        ! control volume after remap
        cv2(ix, iy) = cv(ix, iy) * (1.0_num + dv)

        ! dxb before remap
        dxb1(ix, iy) = dxb(ix) + (vxb - vxbm) * dt
      END DO
    END DO

    DO iy = -1, ny+1
      DO ix = -1, nx+1
        dxc1(ix, iy) = 0.5_num * (dxb1(ix, iy) + dxb1(ix+1, iy))
      END DO
    END DO

    ! Evans and Hawley (ApJ, vol 332, p650, (1988))
    ! constrained transport remap of magnetic fluxes
    CALL vx_by_flux
    DO iy = 0, ny
      DO ix = 1, nx
        ixm = ix - 1
        by(ix, iy) = by(ix, iy) - flux(ix, iy) + flux(ixm, iy)
      END DO
    END DO

    DO iy = 1, ny
      iym = iy - 1
      DO ix = 0, nx
        bx(ix, iy) = bx(ix, iy) + flux(ix, iy) - flux(ix, iym)
      END DO
    END DO

    CALL vx_bz_flux
    DO iy = 1, ny
      DO ix = 1, nx
        ixm = ix - 1
        bz(ix, iy) = bz(ix, iy) - flux(ix, iy) + flux(ixm, iy)
      END DO
    END DO

    ! remap of mass + calculation of mass fluxes (dm) needed for later remaps
    CALL x_mass_flux ! calculates dm(0:nx, 0:ny+1)
    CALL dm_x_bcs    ! need dm(-1:nx+1, 0:ny+1) for velocity remap
    DO iy = 1, ny
      DO ix = 1, nx
        ixm = ix - 1
        rho(ix, iy) = (rho1(ix, iy) * cv1(ix, iy) &
            + dm(ixm, iy) - dm(ix, iy)) / cv2(ix, iy)
      END DO
    END DO

    ! remap specific energy density using mass coordinates
    CALL x_energy_flux
    DO iy = 1, ny
      DO ix = 1, nx
        ixm = ix - 1
        energy(ix, iy) = (energy(ix, iy) * cv1(ix, iy) * rho1(ix, iy) &
            + flux(ixm, iy) - flux(ix, iy)) / (cv2(ix, iy) * rho(ix, iy))
      END DO
    END DO

    ! redefine dxb1, cv1, cv2, dm and vx1 for velocity (vertex) cells
    ! in some of these calculations the flux variable is used as a
    ! temporary array
    DO iy = 0, ny
      iyp = iy + 1
      DO ix = -1, nx+1
        ixp = ix + 1

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
      iyp = iy + 1
      DO ix = 0, nx
        ixp = ix + 1
        flux(ix, iy) = cv1(ix, iy) + cv1(ixp, iy) + cv1(ix, iyp) + cv1(ixp, iyp)
      END DO
    END DO
    ! cv1 = vertex CV before remap
    cv1(0:nx, 0:ny) = flux(0:nx, 0:ny) * 0.25_num

    DO iy = 0, ny
      iyp = iy + 1
      DO ix = 0, nx
        ixp = ix + 1
        flux(ix, iy) = cv2(ix, iy) + cv2(ixp, iy) + cv2(ix, iyp) + cv2(ixp, iyp)
      END DO
    END DO
    ! cv2 = vertex CV after remap
    cv2(0:nx, 0:ny) = flux(0:nx, 0:ny) * 0.25_num

    DO iy = 0, ny
      DO ix = -2, nx+1
        ixp = ix + 1
        flux(ix, iy) = (vx1(ix, iy) + vx1(ixp, iy)) * 0.5_num
      END DO
    END DO
    ! vertex velocity used in remap
    vx1(-2:nx+1, 0:ny) = flux(-2:nx+1, 0:ny)

    DO iy = 0, ny
      DO ix = -1, nx+1
        ixm = ix - 1
        ! dxb1 = width of vertex CV before remap
        dxb1(ix, iy) = dxc(ix) + (vx1(ix, iy) - vx1(ixm, iy)) * dt
      END DO
    END DO

    DO iy = 0, ny
      iyp = iy + 1
      DO ix = -1, nx
        ixp = ix + 1
        flux(ix, iy) = dm(ix, iy) + dm(ixp, iy) + dm(ix, iyp) + dm(ixp, iyp)
      END DO
    END DO
    ! mass flux out of vertex CV
    dm(-1:nx, 0:ny) = flux(-1:nx, 0:ny) * 0.25_num

    DO iy = 0, ny
      DO ix = 0, nx
        ixm = ix - 1
        ! vertex density after remap
        rho_v1(ix, iy) = (rho_v(ix, iy) * cv1(ix, iy) &
            + dm(ixm, iy) - dm(ix, iy)) / cv2(ix, iy)
      END DO
    END DO

    CALL x_momy_flux
    DO iy = 0, ny
      DO ix = 0, nx
        ixm = ix - 1
        vy(ix, iy) = (rho_v(ix, iy) * vy(ix, iy) * cv1(ix, iy) &
            + flux(ixm, iy) - flux(ix, iy)) / (cv2(ix, iy) * rho_v1(ix, iy))
      END DO
    END DO

    CALL x_momz_flux
    DO iy = 0, ny
      DO ix = 0, nx
        ixm = ix - 1
        vz(ix, iy) = (rho_v(ix, iy) * vz(ix, iy) * cv1(ix, iy) &
            + flux(ixm, iy) - flux(ix, iy)) / (cv2(ix, iy) * rho_v1(ix, iy))
      END DO
    END DO

    CALL x_momx_flux
    DO iy = 0, ny
      DO ix = 0, nx
        ixm = ix - 1
        vx(ix, iy) = (rho_v(ix, iy) * vx(ix, iy) * cv1(ix, iy) &
            + flux(ixm, iy) - flux(ix, iy)) / (cv2(ix, iy) * rho_v1(ix, iy))
      END DO
    END DO

    CALL boundary_conditions

    DEALLOCATE (rho1, dm, cv2, flux, dxb1, rho_v, rho_v1, dxc1)
    xpass = 0

  END SUBROUTINE remap_x



  SUBROUTINE vx_by_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: db, dbxp, dbxp2, dbxm
    INTEGER :: ixp2

    DO iy = 0, ny
      iyp = iy + 1
      DO ix = 0, nx
        ixm  = ix - 1
        ixp  = ix + 1
        ixp2 = ix + 2

        v_advect = vx1(ix, iy)

        db    = (dxb1(ix  , iy) + dxb1(ix  , iyp)) * 0.5_num
        dbxp  = (dxb1(ixp , iy) + dxb1(ixp , iyp)) * 0.5_num
        dbxp2 = (dxb1(ixp2, iy) + dxb1(ixp2, iyp)) * 0.5_num
        dbxm  = (dxb1(ixm , iy) + dxb1(ixm , iyp)) * 0.5_num

        w4 = by(ix , iy) / db
        w5 = by(ixp, iy) / dbxp

        flux(ix, iy) = (MAX(0.0_num, v_advect) * w4 &
            + MIN(0.0_num, v_advect) * w5) * dt

        w1 = by(ixp , iy) / dbxp  - by(ix , iy) / db
        w2 = by(ix  , iy) / db    - by(ixm, iy) / dbxm
        w3 = by(ixp2, iy) / dbxp2 - by(ixp, iy) / dbxp

        ! vad_p and vad_m are logical switches which determine v_advect>=0
        ! and v_advect<0 respectively. It's written this way to allow vector
        ! optimization

        vad_p = MAX(SIGN(1.0_num, v_advect  ), 0.0_num)
        vad_m = -MIN(SIGN(1.0_num, v_advect  ), 0.0_num)

        w7 = ABS(v_advect) * dt / (db * vad_p + dbxp * vad_m)

        w9 = (2.0_num - w7) * ABS(w1) / dxc1(ix, iy) &
            + (1.0_num + w7) * (ABS(w2) / dxc1(ixm, iy) * vad_p &
            + ABS(w3) / dxc1(ixp, iy) * vad_m)

        w9 = w9 * sixth

        w8 = 0.5_num * (SIGN(1.0_num, w1) &
            + SIGN(1.0_num, w2 * vad_p + w3 * vad_m))

        w6 = SIGN(1.0_num, v_advect) * w8 &
            * MIN(ABS(w9) * (dxb1(ix, iy) * vad_p + dxb1(ixp, iy) * vad_m), &
            ABS(w1), ABS(w2 * vad_p + w3 * vad_m))

        flux(ix, iy) = flux(ix, iy) &
            + v_advect * dt * w6 * (1.0_num - w7)
      END DO
    END DO

  END SUBROUTINE vx_by_flux



  SUBROUTINE vx_bz_flux

    REAL(num) :: v_advect, vad_p, vad_m
    INTEGER :: ixp2

    DO iy = 1, ny
      iym = iy - 1
      DO ix = 0, nx
        ixm  = ix - 1
        ixp  = ix + 1
        ixp2 = ix + 2

        v_advect = (vx1(ix, iy) + vx1(ix, iym)) * 0.5_num

        w4 = bz(ix , iy) / dxb1(ix , iy)
        w5 = bz(ixp, iy) / dxb1(ixp, iy)

        flux(ix, iy) = (MAX(0.0_num, v_advect) * w4 &
            + MIN(0.0_num, v_advect) * w5) * dt

        w1 = bz(ixp , iy) / dxb1(ixp , iy) - bz(ix , iy) / dxb1(ix , iy)
        w2 = bz(ix  , iy) / dxb1(ix  , iy) - bz(ixm, iy) / dxb1(ixm, iy)
        w3 = bz(ixp2, iy) / dxb1(ixp2, iy) - bz(ixp, iy) / dxb1(ixp, iy)

        ! vad_p and vad_m are logical switches which determine v_advect>=0
        ! and v_advect<0 respectively. It's written this way to allow vector
        ! optimization

        vad_p = MAX(SIGN(1.0_num, v_advect  ), 0.0_num)
        vad_m = - MIN(SIGN(1.0_num, v_advect  ), 0.0_num)

        w7 = ABS(v_advect) * dt / (dxb1(ix, iy) * vad_p + dxb1(ixp, iy) * vad_m)

        w9 = (2.0_num - w7) * ABS(w1) / dxc1(ix, iy) &
            + (1.0_num + w7) * (ABS(w2) / dxc1(ixm, iy) * vad_p &
            + ABS(w3) / dxc1(ixp, iy) * vad_m)

        w9 = w9 * sixth
        w8 = 0.5_num * (SIGN(1.0_num, w1) &
            + SIGN(1.0_num, w2 * vad_p + w3 * vad_m))

        w6 = SIGN(1.0_num, v_advect) * w8 &
            * MIN(ABS(w9) * (dxb1(ix, iy) * vad_p + dxb1(ixp, iy) * vad_m), &
            ABS(w1), ABS(w2 * vad_p + w3 * vad_m))

        flux(ix, iy) = flux(ix, iy) &
            +  v_advect * dt * w6 * (1.0_num - w7)
      END DO
    END DO

  END SUBROUTINE vx_bz_flux



  SUBROUTINE x_mass_flux

    REAL(num) :: v_advect, flux_rho, vad_p, vad_m
    INTEGER :: ixp2

    DO iy = 0, ny+1
      iym = iy - 1
      DO ix = 0, nx
        ixm  = ix - 1
        ixp  = ix + 1
        ixp2 = ix + 2

        v_advect = (vx1(ix, iy) + vx1(ix, iym)) * 0.5_num

        dm(ix, iy) = (MAX(0.0_num, v_advect) * rho(ix, iy) &
            + MIN(0.0_num, v_advect) * rho(ixp, iy)) * dt

        w1 = rho(ixp , iy) - rho(ix , iy)
        w2 = rho(ix  , iy) - rho(ixm, iy)
        w3 = rho(ixp2, iy) - rho(ixp, iy)

        ! vad_p and vad_m are logical switches which determine v_advect>=0
        ! and v_advect<0 respectively. It's written this way to allow vector
        ! optimization

        vad_p = MAX(SIGN(1.0_num, v_advect  ), 0.0_num)
        vad_m = - MIN(SIGN(1.0_num, v_advect  ), 0.0_num)

        w5 = ABS(v_advect) * dt / (dxb1(ix, iy) * vad_p + dxb1(ixp, iy) * vad_m)

        w4 = (2.0_num - w5) * ABS(w1) / dxc1(ix, iy) &
            + (1.0_num + w5) * (ABS(w2) / dxc1(ixm, iy) * vad_p &
            + ABS(w3) / dxc1(ixp, iy) * vad_m)

        w4 = w4 * sixth
        w8 = 0.5_num * (SIGN(1.0_num, w1) &
            + SIGN(1.0_num, w2 * vad_p + w3 * vad_m))

        w6 = SIGN(1.0_num, v_advect) * w8 &
            * MIN(ABS(w4) * (dxb1(ix, iy) * vad_p + dxb1(ixp, iy) * vad_m), &
            ABS(w1), ABS(w2 * vad_p + w3 * vad_m))

        flux_rho = v_advect * dt * w6 * (1.0_num - w5)

        dm(ix, iy) = (flux_rho + dm(ix, iy)) * dyb(iy)
      END DO
    END DO

  END SUBROUTINE x_mass_flux



  SUBROUTINE x_energy_flux ! energy remap in mass coordinates

    REAL(num) :: v_advect, vad_p, vad_m
    INTEGER :: ixp2

    DO iy = 0, ny
      iym = iy - 1
      DO ix = 0, nx
        ixm  = ix - 1
        ixp  = ix + 1
        ixp2 = ix + 2

        v_advect = (vx1(ix, iy) + vx1(ix, iym)) * 0.5_num

        w1 = energy(ixp , iy) - energy(ix , iy)
        w2 = energy(ix  , iy) - energy(ixm, iy)
        w3 = energy(ixp2, iy) - energy(ixp, iy)

        ! vad_p and vad_m are logical switches which determine v_advect>=0
        ! and v_advect<0 respectively. It's written this way to allow vector
        ! optimization

        vad_p = MAX(SIGN(1.0_num, v_advect  ), 0.0_num)
        vad_m = - MIN(SIGN(1.0_num, v_advect  ), 0.0_num)

        w5 = ABS(v_advect) * dt / (dxb1(ix, iy) * vad_p + dxb1(ixp, iy) * vad_m)

        w7 = energy(ix, iy) * vad_p + energy(ixp, iy) * vad_m

        w6 = ABS(dm(ix, iy)) / ((rho1(ix, iy) * dxb1(ix, iy) * vad_p + &
            rho1(ixp, iy) * dxb1(ixp, iy) * vad_m) * dyb(iy))

        w4 = (2.0_num - w5) * ABS(w1) / dxc(ix) &
            + (1.0_num + w5) * (ABS(w2) / dxc(ixm) * vad_p &
            + ABS(w3) / dxc(ixp) * vad_m)

        w4 = w4 * sixth
        w8 = 0.5_num * (SIGN(1.0_num, w1) &
            + SIGN(1.0_num, w2 * vad_p + w3 * vad_m))

        w9 = SIGN(1.0_num, v_advect) * w8 &
            * MIN(ABS(w4) * (dxb(ix) * vad_p + dxb(ixp) * vad_m), &
            ABS(w1), ABS(w2 * vad_p + w3 * vad_m))

        flux(ix, iy) = dm(ix, iy) * (w7 + w9 * (1.0_num - w6))
      END DO
    END DO

  END SUBROUTINE x_energy_flux



  SUBROUTINE x_momy_flux ! energy remap in mass coordinates

    REAL(num) :: v_advect, m, mp, ai, aip, dk, vad_p, vad_m
    INTEGER :: ixp2

    DO iy = 0, ny
      DO ix = -1, nx
        ixm  = ix - 1
        ixp  = ix + 1
        ixp2 = ix + 2

        v_advect = vx1(ix, iy)

        w1 = vy(ixp , iy) - vy(ix , iy)
        w2 = vy(ix  , iy) - vy(ixm, iy)
        w3 = vy(ixp2, iy) - vy(ixp, iy)

        ! vad_p and vad_m are logical switches which determine v_advect>=0
        ! and v_advect<0 respectively. It's written this way to allow vector
        ! optimization

        vad_p = MAX(SIGN(1.0_num, v_advect  ), 0.0_num)
        vad_m = - MIN(SIGN(1.0_num, v_advect  ), 0.0_num)

        w5 = ABS(v_advect) * dt / (dxb1(ix, iy) * vad_p + dxb1(ixp, iy) * vad_m)

        w7 = vy(ix, iy) * vad_p + vy(ixp, iy) * vad_m

        w6 = ABS(dm(ix, iy)) / ((rho_v(ix, iy) * dxb1(ix, iy) * vad_p + &
            rho_v(ixp, iy) * dxb1(ixp, iy) * vad_m) * dyc(iy))

        w4 = (2.0_num - w5) * ABS(w1) / dxb(ixp) &
            + (1.0_num + w5) * (ABS(w2) / dxb(ix) * vad_p &
            + ABS(w3) / dxb(ixp2) * vad_m)

        w4 = w4 * sixth
        w8 = 0.5_num * (SIGN(1.0_num, w1) &
            + SIGN(1.0_num, w2 * vad_p + w3 * vad_m))

        w9 = SIGN(1.0_num, v_advect) * w8 &
            * MIN(ABS(w4) * (dxc(ix) * vad_p + dxc(ixp) * vad_m), &
            ABS(w1), ABS(w2 * vad_p + w3 * vad_m))

        flux(ix, iy) = w7 + w9 * (1.0_num - w6)
      END DO
    END DO

    IF (rke) THEN
      DO iy = 0, ny
        iyp = iy + 1
        DO ix = 0, nx-1
          ixm = ix - 1
          ixp = ix + 1

          m = rho_v1(ix, iy) * cv2(ix, iy)
          mp = rho_v1(ixp, iy) * cv2(ixp, iy)

          ai = (vy(ix, iy) - flux(ixm, iy)) * dm(ixm, iy) / m &
              + (flux(ix, iy) - vy(ix, iy)) * dm(ix, iy) / m

          aip = (vy(ixp, iy) - flux(ix, iy)) * dm(ix, iy) / mp &
              + (flux(ixp, iy) - vy(ixp, iy)) * dm(ixp, iy) / mp

          dk = (vy(ixp, iy) - vy(ix, iy)) * (flux(ix, iy) &
              - 0.5_num * (vy(ixp, iy) + vy(ix, iy))) &
              + 0.5_num * ai * (flux(ix, iy) - vy(ix, iy)) &
              + 0.5_num * aip * (vy(ixp, iy) - flux(ix, iy))

          dk = dk * dm(ix, iy) * 0.5_num
          delta_ke(ixp, iy ) = delta_ke(ixp, iy ) + dk
          delta_ke(ixp, iyp) = delta_ke(ixp, iyp) + dk
        END DO
      END DO
    END IF

    flux(-1:nx, 0:ny) = flux(-1:nx, 0:ny) * dm(-1:nx, 0:ny)

  END SUBROUTINE x_momy_flux



  SUBROUTINE x_momz_flux ! energy remap in mass coordinates

    REAL(num) :: v_advect, m, mp, ai, aip, dk, vad_p, vad_m
    INTEGER :: ixp2

    DO iy = 0, ny
      DO ix = -1, nx
        ixm  = ix - 1
        ixp  = ix + 1
        ixp2 = ix + 2

        v_advect = vx1(ix, iy)

        w1 = vz(ixp , iy) - vz(ix , iy)
        w2 = vz(ix  , iy) - vz(ixm, iy)
        w3 = vz(ixp2, iy) - vz(ixp, iy)

        ! vad_p and vad_m are logical switches which determine v_advect>=0
        ! and v_advect<0 respectively. It's written this way to allow vector
        ! optimization

        vad_p = MAX(SIGN(1.0_num, v_advect  ), 0.0_num)
        vad_m = - MIN(SIGN(1.0_num, v_advect  ), 0.0_num)

        w5 = ABS(v_advect) * dt / (dxb1(ix, iy) * vad_p + dxb1(ixp, iy) * vad_m)

        w7 = vz(ix, iy) * vad_p + vz(ixp, iy) * vad_m

        w6 = ABS(dm(ix, iy)) / ((rho_v(ix, iy) * dxb1(ix, iy) * vad_p + &
            rho_v(ixp, iy) * dxb1(ixp, iy) * vad_m) * dyc(iy))

        w4 = (2.0_num - w5) * ABS(w1) / dxb(ixp) &
            + (1.0_num + w5) * (ABS(w2) / dxb(ix) * vad_p &
            + ABS(w3) / dxb(ixp2) * vad_m)

        w4 = w4 * sixth
        w8 = 0.5_num * (SIGN(1.0_num, w1) &
            + SIGN(1.0_num, w2 * vad_p + w3 * vad_m))

        w9 = SIGN(1.0_num, v_advect) * w8 &
            * MIN(ABS(w4) * (dxc(ix) * vad_p + dxc(ixp) * vad_m), &
            ABS(w1), ABS(w2 * vad_p + w3 * vad_m))

        flux(ix, iy) = w7 + w9 * (1.0_num - w6)
      END DO
    END DO

    IF (rke) THEN
      DO iy = 0, ny
        iyp = iy + 1
        DO ix = 0, nx-1
          ixm = ix - 1
          ixp = ix + 1

          m = rho_v1(ix, iy) * cv2(ix, iy)
          mp = rho_v1(ixp, iy) * cv2(ixp, iy)

          ai = (vz(ix, iy) - flux(ixm, iy)) * dm(ixm, iy) / m &
              + (flux(ix, iy) - vz(ix, iy)) * dm(ix, iy) / m

          aip = (vz(ixp, iy) - flux(ix, iy)) * dm(ix, iy) / mp &
              + (flux(ixp, iy) - vz(ixp, iy)) * dm(ixp, iy) / mp

          dk = (vz(ixp, iy) - vz(ix, iy)) * (flux(ix, iy) &
              - 0.5_num * (vz(ixp, iy) + vz(ix, iy))) &
              + 0.5_num * ai * (flux(ix, iy) - vz(ix, iy)) &
              + 0.5_num * aip * (vz(ixp, iy) - flux(ix, iy))
          dk = dk * dm(ix, iy) * 0.5_num
          delta_ke(ixp, iy ) = delta_ke(ixp, iy ) + dk
          delta_ke(ixp, iyp) = delta_ke(ixp, iyp) + dk
        END DO
      END DO
    END IF

    flux(-1:nx, 0:ny) = flux(-1:nx, 0:ny) * dm(-1:nx, 0:ny)

  END SUBROUTINE x_momz_flux



  SUBROUTINE x_momx_flux

    REAL(num) :: v_advect, m, mp, ai, aip, dk, vad_p, vad_m
    INTEGER :: ixp2

    DO iy = 0, ny
      DO ix = -1, nx
        ixm  = ix - 1
        ixp  = ix + 1
        ixp2 = ix + 2

        v_advect = vx1(ix, iy)

        w1 = vx(ixp , iy) - vx(ix , iy)
        w2 = vx(ix  , iy) - vx(ixm, iy)
        w3 = vx(ixp2, iy) - vx(ixp, iy)

        ! vad_p and vad_m are logical switches which determine v_advect>=0
        ! and v_advect<0 respectively. It's written this way to allow vector
        ! optimization

        vad_p = MAX(SIGN(1.0_num, v_advect  ), 0.0_num)
        vad_m = - MIN(SIGN(1.0_num, v_advect  ), 0.0_num)

        w5 = ABS(v_advect) * dt / (dxb1(ix, iy) * vad_p + dxb1(ixp, iy) * vad_m)

        w7 = vx(ix, iy) * vad_p + vx(ixp, iy) * vad_m

        w6 = ABS(dm(ix, iy)) / ((rho_v(ix, iy) * dxb1(ix, iy) * vad_p &
            + rho_v(ixp, iy) * dxb1(ixp, iy) * vad_m) * dyc(iy))

        w4 = (2.0_num - w5) * ABS(w1) / dxb(ixp) &
            + (1.0_num + w5) * (ABS(w2) / dxb(ix) * vad_p &
            + ABS(w3) / dxb(ixp2) * vad_m)

        w4 = w4 * sixth
        w8 = 0.5_num * (SIGN(1.0_num, w1) &
            + SIGN(1.0_num, w2 * vad_p + w3 * vad_m))

        w9 = SIGN(1.0_num, v_advect) * w8 &
            * MIN(ABS(w4) * (dxc(ix) * vad_p + dxc(ixp) * vad_m), &
            ABS(w1), ABS(w2 * vad_p + w3 * vad_m))

        flux(ix, iy) = w7 + w9 * (1.0_num - w6)
      END DO
    END DO

    IF (rke) THEN
      DO iy = 0, ny
        iyp = iy + 1
        DO ix = 0, nx-1
          ixm = ix - 1
          ixp = ix + 1

          m = rho_v1(ix, iy) * cv2(ix, iy)
          mp = rho_v1(ixp, iy) * cv2(ixp, iy)

          ai = (vx(ix, iy) - flux(ixm, iy)) * dm(ixm, iy) / m &
              + (flux(ix, iy) - vx(ix, iy)) * dm(ix, iy) / m

          aip = (vx(ixp, iy) - flux(ix, iy)) * dm(ix, iy) / mp &
              + (flux(ixp, iy) - vx(ixp, iy)) * dm(ixp, iy) / mp

          dk = (vx(ixp, iy) - vx(ix, iy)) * (flux(ix, iy) &
              - 0.5_num * (vx(ixp, iy) + vx(ix, iy))) &
              + 0.5_num * ai * (flux(ix, iy) - vx(ix, iy)) &
              + 0.5_num * aip * (vx(ixp, iy) - flux(ix, iy))

          dk = dk * dm(ix, iy) * 0.5_num
          delta_ke(ixp, iy ) = delta_ke(ixp, iy ) + dk
          delta_ke(ixp, iyp) = delta_ke(ixp, iyp) + dk
        END DO
      END DO
    END IF

    flux(-1:nx, 0:ny) = flux(-1:nx, 0:ny) * dm(-1:nx, 0:ny)

  END SUBROUTINE x_momx_flux



  SUBROUTINE dm_x_bcs

    CALL MPI_SENDRECV(dm(1, 0:ny+1), ny+2, mpireal, &
        left, tag, dm(nx+1, 0:ny+1), ny+2, mpireal, &
        right, tag, comm, status, errcode)

    IF (right == MPI_PROC_NULL) &
        dm(nx+1, 0:ny+1) = dm(nx, 0:ny+1)

    CALL MPI_SENDRECV(dm(nx-1, 0:ny+1), ny+2, mpireal, &
        right, tag, dm(-1, 0:ny+1), ny+2, mpireal, &
        left, tag, comm, status, errcode)

    IF (left == MPI_PROC_NULL) &
        dm(-1, 0:ny+1) = dm(0, 0:ny+1)

  END SUBROUTINE dm_x_bcs

END MODULE xremap
