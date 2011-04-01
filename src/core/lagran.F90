MODULE lagran

  USE shared_data
  USE boundary
  USE neutral
  USE eos 
  USE conduct

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: lagrangian_step, eta_calc

  ! only used inside lagran.f90
  REAL(num), DIMENSION(:, :), ALLOCATABLE :: qxy, qxz, qyz
  REAL(num), DIMENSION(:, :), ALLOCATABLE :: qxx, qyy, visc_heat, pressure
  REAL(num), DIMENSION(:, :), ALLOCATABLE :: flux_x, flux_y, flux_z, curlb

CONTAINS

  ! This subroutine manages the progress of the lagrangian step
  SUBROUTINE lagrangian_step

    INTEGER :: substeps, subcycle
    REAL(num) :: actual_dt, dt_sub

    ALLOCATE (bx1(0:nx+1, 0:ny+1), by1(0:nx+1, 0:ny+1), bz1(0:nx+1, 0:ny+1), &
        qxy(0:nx+1, 0:ny+1), qxz(0:nx+1, 0:ny+1), qyz(0:nx+1, 0:ny+1), &
        qxx(0:nx+1, 0:ny+1), qyy(0:nx+1, 0:ny+1), &
        visc_heat(-1:nx+2, -1:ny+2), pressure(-1:nx+2, -1:ny+2), &
        flux_x(0:nx, 0:ny), flux_y(0:nx, 0:ny), flux_z(0:nx, 0:ny), &
        curlb(0:nx, 0:ny))
      
    IF (resistive_mhd .OR. hall_mhd) THEN
      ! if subcycling isn't wanted set dt = dtr in set_dt, don't just
      ! set substeps to 1.
      IF (resistive_mhd) THEN
        dt_sub = dtr
      ELSE
        dt_sub = dth
      END IF

      IF (resistive_mhd .AND. hall_mhd) dt_sub = MIN(dtr, dth)
      substeps = INT(dt / dt_sub) + 1

      IF (substeps > peak_substeps) peak_substeps = substeps
      actual_dt = dt
      dt = dt / REAL(substeps, num)

      DO subcycle = 1, substeps
        CALL eta_calc
        IF (include_neutrals) CALL neutral_fraction(eos_number)
        IF (cowling_resistivity) CALL perpendicular_resistivity
        IF (hall_mhd) CALL hall_effects
        IF (resistive_mhd) CALL resistive_effects
      END DO

      dt = actual_dt
    END IF               
    
    IF (conduction) CALL conduct_heat  
     
    DO iy = 0, ny + 1
       DO ix = 0, nx + 1
         ixm = ix - 1
         iym = iy - 1
         bx1(ix, iy) = (bx(ix, iy) + bx(ixm, iy)) * 0.5_num
         by1(ix, iy) = (by(ix, iy) + by(ix, iym)) * 0.5_num
       END DO
     END DO
     bz1 = bz(0:nx+1, 0:ny+1) ! bz1 = bz at C
 
     CALL predictor_corrector_step
 
     DEALLOCATE (bx1, by1, bz1, qxy, qxz, qyz, qxx, qyy, visc_heat, pressure, &
         flux_x, flux_y, flux_z, curlb)
 
      CALL energy_bcs
      CALL density_bcs
      CALL velocity_bcs  

  END SUBROUTINE lagrangian_step



  ! The main predictor / corrector step which advances the momentum equation
  SUBROUTINE predictor_corrector_step

    REAL(num) :: e1, rho_v
    REAL(num) :: fx, fy, fz
    REAL(num) :: vxb, vxbm, vyb, vybm
    REAL(num) :: dv
    REAL(num) :: bxv, byv, bzv, jx, jy, jz
    REAL(num) :: cvx, cvxp, cvy, cvyp

    dt2 = dt * 0.5_num
    CALL viscosity_and_b_update

    bx1 = bx1 * cv1(0:nx+1, 0:ny+1)
    by1 = by1 * cv1(0:nx+1, 0:ny+1)
    bz1 = bz1 * cv1(0:nx+1, 0:ny+1)

    DO iy = 0, ny + 1
      DO ix = 0, nx + 1
        dv = cv1(ix, iy) / cv(ix, iy) - 1.0_num
        ! predictor energy
#ifdef Q_MONO
        ! add shock viscosity
        pressure(ix, iy) = pressure(ix, iy) + p_visc(ix, iy)
#endif
        e1 = energy(ix, iy) - pressure(ix, iy) * dv / rho(ix, iy)
        e1 = e1 + visc_heat(ix, iy) * dt2 / rho(ix, iy)

        ! now define the predictor step pressures
        CALL get_pressure(rho(ix, iy) * cv(ix, iy) / cv1(ix, iy), &
            e1, eos_number, ix, iy, pressure(ix, iy))
#ifdef Q_MONO
        ! add shock viscosity
        pressure(ix, iy) = pressure(ix, iy) + p_visc(ix, iy)
#endif
      END DO
    END DO

    DO iy = 0, ny
      DO ix = 0, nx
        ixp = ix + 1
        iyp = iy + 1

        ! P total at Ey(i, j)
        w1 = (pressure(ix , iy) + pressure(ix , iyp)) * 0.5_num
        ! P total at Ey(i+1, j)
        w2 = (pressure(ixp, iy) + pressure(ixp, iyp)) * 0.5_num
        fx = -(w2 - w1) / dxc(ix)
        ! P total at Ex(i, j)
        w1 = (pressure(ix, iy ) + pressure(ixp, iy )) * 0.5_num
        ! P total at Ex(i, j+1)
        w2 = (pressure(ix, iyp) + pressure(ixp, iyp)) * 0.5_num
        fy = -(w2 - w1) / dyc(iy)
        fz = 0.0_num

        ! add parallel component of viscosity
        w1 = (qxx(ix, iy) + qxx(ix, iyp)) * 0.5_num
        w2 = (qxx(ixp, iy) + qxx(ixp, iyp)) * 0.5_num
        fx = fx + (w2 - w1) / dxc(ix)
        w1 = (qyy(ix, iy) + qyy(ixp, iy)) * 0.5_num
        w2 = (qyy(ix, iyp) + qyy(ixp, iyp)) * 0.5_num
        fy = fy + (w2 - w1) / dyc(iy)

        ! add shear forces
        w1 = (qxy(ix, iy) + qxy(ixp, iy)) * 0.5_num
        w2 = (qxy(ix, iyp) + qxy(ixp, iyp)) * 0.5_num
        fx = fx + (w2 - w1) / dyc(iy)
        w1 = (qxy(ix, iy) + qxy(ix, iyp)) * 0.5_num
        w2 = (qxy(ixp, iy) + qxy(ixp, iyp)) * 0.5_num
        fy = fy + (w2 - w1) / dxc(ix)
        w1 = (qxz(ix, iy) + qxz(ix, iyp)) * 0.5_num
        w2 = (qxz(ixp, iy) + qxz(ixp, iyp)) * 0.5_num
        fz = (w2 - w1) / dxc(ix)
        w1 = (qyz(ix, iy) + qyz(ixp, iy)) * 0.5_num
        w2 = (qyz(ix, iyp) + qyz(ixp, iyp)) * 0.5_num
        fz = fz + (w2 - w1) / dyc(iy)

        cvx = cv1(ix, iy) + cv1(ix, iyp)
        cvxp = cv1(ixp, iy) + cv1(ixp, iyp)
        cvy = cv1(ix, iy) + cv1(ixp, iy)
        cvyp = cv1(ix, iyp) + cv1(ixp, iyp)

        w1 = (bz1(ix, iy) + bz1(ixp, iy)) / cvy
        w2 = (bz1(ix, iyp) + bz1(ixp, iyp)) / cvyp
        jx = (w2 - w1) / dyc(iy)

        w1 = (bz1(ix, iy) + bz1(ix, iyp)) / cvx
        w2 = (bz1(ixp, iy) + bz1(ixp, iyp)) / cvxp
        jy = -(w2 - w1) / dxc(ix)

        w1 = (by1(ix, iy) + by1(ix, iyp)) / cvx
        w2 = (by1(ixp, iy) + by1(ixp, iyp)) / cvxp
        jz = (w2 - w1) / dxc(ix)
        w1 = (bx1(ix, iy) + bx1(ixp, iy)) / cvy
        w2 = (bx1(ix, iyp) + bx1(ixp, iyp)) / cvyp
        jz = jz - (w2 - w1) / dyc(iy)

        bxv = (bx1(ix, iy) + bx1(ixp, iy) + bx1(ix, iyp) + bx1(ixp, iyp)) &
            / (cvx + cvxp)
        byv = (by1(ix, iy) + by1(ixp, iy) + by1(ix, iyp) + by1(ixp, iyp)) &
            / (cvx + cvxp)
        bzv = (bz1(ix, iy) + bz1(ixp, iy) + bz1(ix, iyp) + bz1(ixp, iyp)) &
            / (cvx + cvxp)

        fx = fx + (jy * bzv - jz * byv)
        fy = fy + (jz * bxv - jx * bzv)
        fz = fz + (jx * byv - jy * bxv)

        rho_v = (rho(ix, iy) * cv(ix, iy) + rho(ixp, iy) * cv(ixp, iy) &
            + rho(ix, iyp) * cv(ix, iyp) + rho(ixp, iyp) * cv(ixp, iyp))
        rho_v = rho_v / (cv(ix, iy) + cv(ixp, iy) + cv(ix, iyp) + cv(ixp, iyp))

        fy = fy - rho_v * grav(iy)

        ! find half step velocity needed for remap
        vx1(ix, iy) = vx(ix, iy) + dt2 * fx / rho_v
        vy1(ix, iy) = vy(ix, iy) + dt2 * fy / rho_v
        vz1(ix, iy) = vz(ix, iy) + dt2 * fz / rho_v

        ! velocity at the end of the Lagrangian step
        vx(ix, iy) = vx(ix, iy) + dt * fx / rho_v
        vy(ix, iy) = vy(ix, iy) + dt * fy / rho_v
        vz(ix, iy) = vz(ix, iy) + dt * fz / rho_v

      END DO
    END DO

    CALL remap_v_bcs

    CALL visc_heating

    ! finally correct density and energy to final values
    DO iy = 1, ny
      DO ix = 1, nx
        ixm = ix - 1
        iym = iy - 1

        vxb = (vx1(ix, iy) + vx1(ix, iym)) * 0.5_num     ! vx1 at Ex(i, j)
        vxbm = (vx1(ixm, iy) + vx1(ixm, iym)) * 0.5_num  ! vx1 at Ex(i-1, j)
        vyb = (vy1(ix, iy) + vy1(ixm, iy)) * 0.5_num     ! vy1 at Ey(i, j)
        vybm = (vy1(ix, iym) + vy1(ixm, iym)) * 0.5_num  ! vy1 at Ey(i-1, j)
        dv = ((vxb - vxbm) / dxb(ix) + (vyb - vybm) / dyb(iy)) * dt
        cv1(ix, iy) = cv(ix, iy) * (1.0_num + dv)

        ! energy at end of Lagrangian step
        energy(ix, iy) = energy(ix, iy) &
            - pressure(ix, iy) * dv / rho(ix, iy)

        energy(ix, iy) = energy(ix, iy) &
              + dt * visc_heat(ix, iy) / rho(ix, iy)

        rho(ix, iy) = rho(ix, iy) / (1.0_num + dv)

        total_visc_heating = total_visc_heating &
              + dt * visc_heat(ix, iy) * cv(ix, iy)

#ifdef Q_MONO
        total_visc_heating = total_visc_heating &
            - p_visc(ix, iy) * dv * cv(ix, iy)
#endif

      END DO
    END DO

  END SUBROUTINE predictor_corrector_step



  ! This subroutine calculates the viscous effects and updates the
  ! magnetic field
  SUBROUTINE viscosity_and_b_update

    REAL(num) :: vxb, vxbm, vyb, vybm, vzb, vzbm, dvxdy, dvydx
    REAL(num) :: p, pxp, pxm, pyp, pym, fx, fy, dv
    REAL(num) :: dvxdx, dvydy, dvxy, dvzdx, dvzdy, s, L, cf, L2
    REAL(num) :: sxx, syy, sxy, sxz, syz
    REAL(num) :: case1, case2, flip, cs

    p_visc = 0.0_num

    DO iy = -1, ny + 2
      DO ix = -1, nx + 2
        CALL get_pressure(rho(ix, iy), energy(ix, iy), eos_number, &
            ix, iy, pressure(ix, iy))
      END DO
    END DO

    DO iy = 0, ny + 1
      DO ix = 0, nx + 1
        ixp = ix + 1
        iyp = iy + 1
        ixm = ix - 1
        iym = iy - 1

        vxb = (vx(ix, iy) + vx(ix, iym)) * 0.5_num     ! vx at Ex(i, j)
        vxbm = (vx(ixm, iy) + vx(ixm, iym)) * 0.5_num  ! vx at Ex(i-1, j)
        vyb = (vy(ix, iy) + vy(ixm, iy)) * 0.5_num     ! vy at Ey(i, j)
        vybm = (vy(ix, iym) + vy(ixm, iym)) * 0.5_num  ! vy at Ey(i-1, j)
        dv = ((vxb - vxbm) / dxb(ix) + (vyb - vybm) / dyb(iy)) * dt2
        cv1(ix, iy) = cv(ix, iy) * (1.0_num + dv)

        dvxdx = (vxb - vxbm) / dxb(ix)
        dvydy = (vyb - vybm) / dyb(iy)

        vxb = (vx(ix, iy) + vx(ixm, iy)) * 0.5_num     ! vx at Ey(i, j)
        vxbm = (vx(ix, iym) + vx(ixm, iym)) * 0.5_num  ! vx at Ey(i, j-1)
        vyb = (vy(ix, iy) + vy(ix, iym)) * 0.5_num     ! vy at Ex(i, j)
        vybm = (vy(ixm, iy) + vy(ixm, iym)) * 0.5_num  ! vy at Ex(i-1, j)

        dvxdy = (vxb - vxbm) / dyb(iy)
        dvydx = (vyb - vybm) / dxb(ix)
        dvxy = dvxdy + dvydx

        sxy = dvxy * 0.5_num
        sxx = 2.0_num * dvxdx * third - dvydy * third
        syy = 2.0_num * dvydy * third - dvxdx * third

        vzb = (vz(ix, iy) + vz(ix, iym)) * 0.5_num
        vzbm = (vz(ixm, iy) + vz(ixm, iym)) * 0.5_num
        dvzdx = (vzb - vzbm) / dxb(ix)
        sxz = dvzdx * 0.5_num

        vzb = (vz(ix, iy) + vz(ixm, iy)) * 0.5_num
        vzbm = (vz(ix, iym) + vz(ixm, iym)) * 0.5_num
        dvzdy = (vzb - vzbm) / dyb(iy)
        syz = dvzdy * 0.5_num

        p = pressure(ix, iy)
        pxp = pressure(ixp, iy)
        pxm = pressure(ixm, iy)
        pyp = pressure(ix, iyp)
        pym = pressure(ix, iym)

        fx = -(pxp - pxm) / dxb(ix)
        fy = -(pyp - pym) / dyb(iy)

        w1 = fx**2 + fy**2
        s = (dvxdx * fx**2 + dvydy * fy**2 + dvxy * fx * fy) &
            / MAX(w1, none_zero)

        ! (dxb(ix) * ABS(fy) > dyb(iy) * ABS(fx)) ? 1: - 1
        flip = SIGN(1.0_num, dxb(ix) * ABS(fy) - dyb(iy) * ABS(fx))
        case1 = -MIN(flip, 0.0_num) ! (flip < 0) ? 1:0
        case2 =  MAX(flip, 0.0_num) ! (flip > 0) ? 1:0

        w2 =  dxb(ix)**2 * w1 * case1 / MAX(fx**2, none_zero) &
            + dyb(iy)**2 * w1 * case2 / MAX(fy**2, none_zero) &
            + (dxb(ix)**2 + dyb(iy)**2) * (1.0_num - case1) * (1.0_num - case2)

        flip = -MIN(SIGN(1.0_num, w1 - 1.e-6_num), 0.0_num) ! (w1 < 1.e-6) ? 1:0

        w2 = w2 * (1.0_num - flip) + flip * MIN(dxb(ix), dyb(iy))**2
        L = SQRT(w2)

        L2 = L
        IF (s > 0.0_num .OR. dv > 0.0_num) L = 0.0_num

        w1 = (bx1(ix, iy)**2 + by1(ix, iy)**2 + bz1(ix, iy)**2) / rho(ix, iy)
        cs = SQRT(gamma*(gamma-1.0_num)*energy(ix,iy))
        cf = SQRT(cs**2 + w1)

        p_visc(ix, iy) = visc1 * ABS(s) * L * cf * rho(ix, iy) &
            + visc2 * (s * L)**2 * rho(ix, iy)

        qxy(ix, iy) = 0.0_num
        qxz(ix, iy) = 0.0_num
        qyz(ix, iy) = 0.0_num
        qxx(ix, iy) = 0.0_num
        qyy(ix, iy) = 0.0_num
#ifndef Q_MONO
        qxy(ix, iy) = sxy * (L2 * rho(ix, iy) &
            * (visc1 * cf + L2 * visc2 * ABS(sxy))) 
        qxz(ix, iy) = sxz * (L2 * rho(ix, iy) &   
            * (visc1 * cf + L2 * visc2 * ABS(sxz))) 
        qyz(ix, iy) = syz * (L2 * rho(ix, iy) &
            * (visc1 * cf + L2 * visc2 * ABS(syz))) 
        qxx(ix, iy) = sxx * (L2 * rho(ix, iy) &
            * (visc1 * cf + L2 * visc2 * ABS(sxx))) 
        qyy(ix, iy) = syy * (L2 * rho(ix, iy) &
            * (visc1 * cf + L2 * visc2 * ABS(syy))) 
#endif
        qxy(ix, iy) = qxy(ix, iy) + 2.0_num * sxy * rho(ix, iy) * visc3 
        qxz(ix, iy) = qxz(ix, iy) + 2.0_num * sxz * rho(ix, iy) * visc3 
        qyz(ix, iy) = qyz(ix, iy) + 2.0_num * syz * rho(ix, iy) * visc3 
        qxx(ix, iy) = qxx(ix, iy) + 2.0_num * sxx * rho(ix, iy) * visc3 
        qyy(ix, iy) = qyy(ix, iy) + 2.0_num * syy * rho(ix, iy) * visc3 

        visc_heat(ix, iy) = qxy(ix, iy) * dvxy + qxz(ix, iy) * dvzdx &
              + qyz(ix, iy) * dvzdy + qxx(ix, iy) * dvxdx + qyy(ix, iy) * dvydy

        w3 = bx1(ix, iy) * dvxdx + by1(ix, iy) * dvxdy
        w4 = bx1(ix, iy) * dvydx + by1(ix, iy) * dvydy
        w5 = bx1(ix, iy) * dvzdx + by1(ix, iy) * dvzdy

        bx1(ix, iy) = (bx1(ix, iy) + w3 * dt2) / (1.0_num + dv)
        by1(ix, iy) = (by1(ix, iy) + w4 * dt2) / (1.0_num + dv)
        bz1(ix, iy) = (bz1(ix, iy) + w5 * dt2) / (1.0_num + dv)

      END DO
    END DO

  END SUBROUTINE viscosity_and_b_update



  ! Calculate the viscous heating
  SUBROUTINE visc_heating

    REAL(num) :: vxb, vxbm, vyb, vybm, vzb, vzbm
    REAL(num) :: dvxdx, dvydy, dvxy, dvxz, dvyz

    DO iy = 0, ny + 1
      DO ix = 0, nx + 1
        ixp = ix + 1
        iyp = iy + 1
        ixm = ix - 1
        iym = iy - 1
        vxb = (vx1(ix, iy) + vx1(ix, iym)) * 0.5_num     ! vx at Ex(i, j)
        vxbm = (vx1(ixm, iy) + vx1(ixm, iym)) * 0.5_num  ! vx at Ex(i-1, j)
        vyb = (vy1(ix, iy) + vy1(ixm, iy)) * 0.5_num     ! vy at Ey(i, j)
        vybm = (vy1(ix, iym) + vy1(ixm, iym)) * 0.5_num  ! vy at Ey(i-1, j)

        dvxdx = (vxb - vxbm) / dxb(ix)
        dvydy = (vyb - vybm) / dyb(iy)

        vxb = (vx1(ix, iy) + vx1(ixm, iy)) * 0.5_num     ! vx at Ey(i, j)
        vxbm = (vx1(ix, iym) + vx1(ixm, iym)) * 0.5_num  ! vx at Ey(i, j-1)
        vyb = (vy1(ix, iy) + vy1(ix, iym)) * 0.5_num     ! vy at Ex(i, j)
        vybm = (vy1(ixm, iy) + vy1(ixm, iym)) * 0.5_num  ! vy at Ex(i-1, j)

        dvxy = (vxb - vxbm) / dyb(iy) + (vyb - vybm) / dxb(ix)

        vzb = (vz1(ix, iy) + vz1(ix, iym)) * 0.5_num
        vzbm = (vz1(ixm, iy) + vz1(ixm, iym)) * 0.5_num
        dvxz = (vzb - vzbm) / dxb(ix)

        vzb = (vz1(ix, iy) + vz1(ixm, iy)) * 0.5_num
        vzbm = (vz1(ix, iym) + vz1(ixm, iym)) * 0.5_num
        dvyz = (vzb - vzbm) / dyb(iy)

        visc_heat(ix, iy) = qxy(ix, iy) * dvxy + qxz(ix, iy) * dvxz &
            + qyz(ix, iy) * dvyz + qxx(ix, iy) * dvxdx + qyy(ix, iy) * dvydy
      END DO
    END DO

    visc_heat = MAX(visc_heat, 0.0_num)

  END SUBROUTINE visc_heating



  ! Calculate the spatial profile of the resistivity at the current timestep
  ! Note that this is a core routine so it works in normalised units
  ! This includes lengths etc.
  SUBROUTINE eta_calc

    REAL(num) :: jx, jy, jz, jx1, jx2, jy1, jy2
    INTEGER :: ixp, iyp

    IF (resistive_mhd) THEN
      DO iy = -1, ny + 1
        DO ix = -1, nx + 1
          ixp = ix + 1
          iyp = iy + 1
  
          jx1 = (bz(ix, iyp) - bz(ix, iy)) / dyc(iy)
          jx2 = (bz(ixp, iyp) - bz(ixp, iy)) / dyc(iy)
          jy1 = -(bz(ixp, iy) - bz(ix, iy)) / dxc(ix)
          jy2 = -(bz(ixp, iyp) - bz(ix, iyp)) / dxc(ix)
          jx = (jx1 + jx2) * 0.5_num
          jy = (jy1 + jy2) * 0.5_num
          jz = (by(ixp, iy) - by(ix, iy)) / dxc(ix) &
              - (bx(ix, iyp) - bx(ix, iy)) / dyc(iy)
  
          IF (SQRT(jx**2 + jy**2 + jz**2) .GT. j_max) THEN
            eta(ix, iy) = eta_background + eta0
          ELSE
            eta(ix, iy) = eta_background
          END IF
  
        END DO
      END DO
    ELSE
        eta = 0.0_num        
    END IF

  END SUBROUTINE eta_calc



  ! Calculate the effect of resistivity on the magnetic field and Ohmic heating
  ! Use the subroutine rkstep
  SUBROUTINE resistive_effects

    REAL(num) :: dt6, half_dt
    REAL(num) :: jx1, jx2, jy1, jy2
#ifdef Q_FOURTHORDER    
    REAL(num), DIMENSION(:, :), ALLOCATABLE :: k1x, k2x, k3x, k4x
    REAL(num), DIMENSION(:, :), ALLOCATABLE :: k1y, k2y, k3y, k4y
    REAL(num), DIMENSION(:, :), ALLOCATABLE :: k1z, k2z, k3z, k4z
    REAL(num), DIMENSION(:, :), ALLOCATABLE :: c1, c2, c3, c4

    ALLOCATE(k1x(0:nx, 0:ny), k2x(0:nx, 0:ny), k3x(0:nx, 0:ny), k4x(0:nx, 0:ny))
    ALLOCATE(k1y(0:nx, 0:ny), k2y(0:nx, 0:ny), k3y(0:nx, 0:ny), k4y(0:nx, 0:ny))
    ALLOCATE(k1z(0:nx, 0:ny), k2z(0:nx, 0:ny), k3z(0:nx, 0:ny), k4z(0:nx, 0:ny))
    ALLOCATE(c1(0:nx, 0:ny), c2(0:nx, 0:ny), c3(0:nx, 0:ny), c4(0:nx, 0:ny))
#endif

    bx1 = bx(0:nx+1, 0:ny+1)
    by1 = by(0:nx+1, 0:ny+1)
    bz1 = bz(0:nx+1, 0:ny+1)

    ! step 1
    CALL rkstep

#ifndef Q_FOURTHORDER
    DO iy = 0, ny
      DO ix = 1, nx
        by(ix, iy) = by1(ix, iy) &
            + (-flux_z(ix, iy) + flux_z(ix-1, iy)) * dt / dxb(ix)
      END DO
    END DO
    
    DO iy = 1, ny
      DO ix = 0, nx
        bx(ix, iy) = bx1(ix, iy) &
            + (flux_z(ix, iy) - flux_z(ix, iy-1)) * dt / dyb(iy)
      END DO
    END DO
    
    DO iy = 1, ny
      DO ix = 1, nx
        ixm = ix - 1
        iym = iy - 1
        w1 = -flux_x(ix, iy) - flux_x(ixm, iy) + flux_x(ixm, iym) + flux_x(ix, iym)
        w2 = w1 + flux_y(ix, iy) + flux_y(ix, iym) - flux_y(ixm, iy) - flux_y(ixm, iym)
        bz(ix, iy) = bz1(ix, iy) + w2 * dt / cv(ix, iy)
      END DO
    END DO
    
    CALL bfield_bcs
    
    DO iy = 1, ny
      DO ix = 1, nx
        ixm = ix - 1
        iym = iy - 1
        energy(ix, iy) = energy(ix, iy) &
            + (curlb(ix, iy) + curlb(ixm, iy) + curlb(ix, iym) + curlb(ixm, iym)) &
            * dt / (4.0_num * rho(ix, iy))
      END DO
    END DO  
    CALL energy_bcs

    DO iy = 0, ny
      DO ix = 0, nx
        w1 = dt * dxc(ix) * dyc(iy) * curlb(ix, iy)
        IF ((ix == 0) .OR. (ix == nx)) THEN
          w1 = w1 * 0.5_num
        END IF
        IF ((iy == 0) .OR. (iy == ny)) THEN
          w1 = w1 * 0.5_num
        END IF
        total_ohmic_heating = total_ohmic_heating + w1
      END DO
    END DO    
#else
    half_dt = dt * 0.5_num 
    dt6 = dt * sixth
          
    k1x = flux_x
    k1y = flux_y
    k1z = flux_z
    c1 = curlb
    ! step 2
    DO iy = 1, ny
      DO ix = 0, nx
        bx(ix, iy) = bx1(ix, iy) + (k1z(ix, iy) - k1z(ix, iy-1)) * half_dt / dyb(iy)
      END DO
    END DO

    DO iy = 0, ny
      DO ix = 1, nx
        by(ix, iy) = by1(ix, iy) + (-k1z(ix, iy) + k1z(ix-1, iy)) * half_dt / dxb(ix)
      END DO
    END DO

    DO iy = 1, ny
      DO ix = 1, nx
        ixm = ix - 1
        iym = iy - 1
        w1 = -k1x(ix, iy) - k1x(ixm, iy) + k1x(ixm, iym) + k1x(ix, iym)
        w2 = w1 + k1y(ix, iy) + k1y(ix, iym) - k1y(ixm, iy) - k1y(ixm, iym)
        bz(ix, iy) = bz1(ix, iy) + w2 * half_dt / cv(ix, iy)
      END DO
    END DO

    CALL bfield_bcs
    CALL rkstep
    k2x = flux_x
    k2y = flux_y
    k2z = flux_z
    c2 = curlb

    ! step 3
    DO iy = 1, ny
      DO ix = 0, nx
        bx(ix, iy) = bx1(ix, iy) + (k2z(ix, iy) - k2z(ix, iy-1)) * half_dt / dyb(iy)
      END DO
    END DO

    DO iy = 0, ny
      DO ix = 1, nx
        by(ix, iy) = by1(ix, iy) + (-k2z(ix, iy) + k2z(ix-1, iy)) * half_dt / dxb(ix)
      END DO
    END DO

    DO iy = 1, ny
      DO ix = 1, nx
        ixm = ix - 1
        iym = iy - 1
        w1 = -k2x(ix, iy) - k2x(ixm, iy) + k2x(ixm, iym) + k2x(ix, iym)
        w2 = w1 + k2y(ix, iy) + k2y(ix, iym) - k2y(ixm, iy) - k2y(ixm, iym)
        bz(ix, iy) = bz1(ix, iy) + w2 * half_dt / cv(ix, iy)
      END DO
    END DO

    CALL bfield_bcs
    CALL rkstep
    k3x = flux_x
    k3y = flux_y
    k3z = flux_z
    c3 = curlb

    ! step 4
    DO iy = 1, ny
      DO ix = 0, nx
        bx(ix, iy) = bx1(ix, iy) + (k3z(ix, iy) - k3z(ix, iy-1)) * dt / dyb(iy)
      END DO
    END DO

    DO iy = 0, ny
      DO ix = 1, nx
        by(ix, iy) = by1(ix, iy) + (-k3z(ix, iy) + k3z(ix-1, iy)) * dt / dxb(ix)
      END DO
    END DO

    DO iy = 1, ny
      DO ix = 1, nx
        ixm = ix - 1
        iym = iy - 1
        w1 = -k3x(ix, iy) - k3x(ixm, iy) + k3x(ixm, iym) + k3x(ix, iym)
        w2 = w1 + k3y(ix, iy) + k3y(ix, iym) - k3y(ixm, iy) - k3y(ixm, iym)
        bz(ix, iy) = bz1(ix, iy) + w2 * dt / cv(ix, iy)
      END DO
    END DO

    CALL bfield_bcs
    CALL rkstep
    k4x = flux_x
    k4y = flux_y
    k4z = flux_z
    c4 = curlb

    ! full update
    k3x = k1x + 2.0_num * k2x + 2.0_num * k3x + k4x
    k3y = k1y + 2.0_num * k2y + 2.0_num * k3y + k4y
    k3z = k1z + 2.0_num * k2z + 2.0_num * k3z + k4z
    c1 = c1 + 2.0_num * c2 + 2.0_num * c3 + c4

    DO iy = 0, ny
      DO ix = 1, nx
        by(ix, iy) = by1(ix, iy) &
            + (-k3z(ix, iy) + k3z(ix-1, iy)) * dt6 / dxb(ix)
      END DO
    END DO

    DO iy = 1, ny
      DO ix = 0, nx
        bx(ix, iy) = bx1(ix, iy) &
            + (k3z(ix, iy) - k3z(ix, iy-1)) * dt6 / dyb(iy)
      END DO
    END DO

    DO iy = 1, ny
      DO ix = 1, nx
        ixm = ix - 1
        iym = iy - 1
        w1 = -k3x(ix, iy) - k3x(ixm, iy) + k3x(ixm, iym) + k3x(ix, iym)
        w2 = w1 + k3y(ix, iy) + k3y(ix, iym) - k3y(ixm, iy) - k3y(ixm, iym)
        bz(ix, iy) = bz1(ix, iy) + w2 * dt6 / cv(ix, iy)
      END DO
    END DO

    CALL bfield_bcs

    DO iy = 1, ny
      DO ix = 1, nx
        ixm = ix - 1
        iym = iy - 1
        energy(ix, iy) = energy(ix, iy) &
            + (c1(ix, iy) + c1(ixm, iy) + c1(ix, iym) + c1(ixm, iym)) &
            * dt6 / (4.0_num * rho(ix, iy))
      END DO
    END DO

    CALL energy_bcs

    DO iy = 0, ny
      DO ix = 0, nx
        w1 = dt6 * dxc(ix) * dyc(iy) * c1(ix, iy)
        IF ((ix == 0) .OR. (ix == nx)) THEN
          w1 = w1 * 0.5_num
        END IF
        IF ((iy == 0) .OR. (iy == ny)) THEN
          w1 = w1 * 0.5_num
        END IF
        total_ohmic_heating = total_ohmic_heating + w1
      END DO
    END DO
#endif

    DO iy = 0, ny + 1
      DO ix = 0, nx + 1
        ixp = ix + 1
        iyp = iy + 1

        jx1 = (bz(ix, iyp) - bz(ix, iy)) / dyc(iy)
        jx2 = (bz(ixp, iyp) - bz(ixp, iy)) / dyc(iy)
        jy1 = -(bz(ixp, iy) - bz(ix, iy)) / dxc(ix)
        jy2 = -(bz(ixp, iyp) - bz(ix, iyp)) / dxc(ix)
        jx_r(ix, iy) = (jx1 + jx2) * 0.5_num
        jy_r(ix, iy) = (jy1 + jy2) * 0.5_num
        jz_r(ix, iy) = (by(ixp, iy) - by(ix, iy)) / dxc(ix) &
            - (bx(ix, iyp) - bx(ix, iy)) / dyc(iy)
      END DO
    END DO

    ! Once more to get j_perp and j_par correct
    CALL rkstep

#ifdef Q_FOURTHORDER
    DEALLOCATE(k1x, k2x, k3x, k4x, k1y, k2y, k3y, k4y, k1z, k2z, k3z, k4z)
    DEALLOCATE(c1, c2, c3, c4)
#endif

  END SUBROUTINE resistive_effects



  ! calculates 'k' values from b[xyz]1 values
  SUBROUTINE rkstep

    REAL(num) :: jx, jy, jz
    REAL(num) :: jx1, jy1, jx2, jy2
    REAL(num) :: bxv, byv, bzv
    REAL(num) :: magn_b
    REAL(num) :: j_par_x, j_par_y, j_par_z
    REAL(num) :: j_perp_x, j_perp_y, j_perp_z
    REAL(num) :: magn_j_perp, magn_j_par

    IF (.NOT. cowling_resistivity) THEN
      DO iy = 0, ny 
        iyp = iy + 1
        DO ix = 0, nx 
          ixp = ix + 1

          jx1 = (bz(ix, iyp) - bz(ix, iy)) / dyc(iy)
          jx2 = (bz(ixp, iyp) - bz(ixp, iy)) / dyc(iy)
          jy1 = -(bz(ixp, iy) - bz(ix, iy)) / dxc(ix)
          jy2 = -(bz(ixp, iyp) - bz(ix, iyp)) / dxc(ix)
          jx = (jx1 + jx2) * 0.5_num
          jy = (jy1 + jy2) * 0.5_num
          jz = (by(ixp, iy) - by(ix, iy)) / dxc(ix) &
              - (bx(ix, iyp) - bx(ix, iy)) / dyc(iy) 
           
          flux_x(ix, iy) = -jx * eta(ix, iy) * dxc(ix) * 0.5_num
          flux_y(ix, iy) = -jy * eta(ix, iy) * dyc(iy) * 0.5_num
          flux_z(ix, iy) = -jz * eta(ix, iy)
          curlb(ix, iy) = eta(ix, iy) * (jx**2 + jy**2 + jz**2)                        
        END DO
      END DO
    ELSE
      ! Use partially ionised flux calculation
      DO iy = 0, ny
        iyp = iy + 1
        DO ix = 0, nx
          ixp = ix + 1

          jx1 = (bz(ix, iyp) - bz(ix, iy)) / dyc(iy)
          jx2 = (bz(ixp, iyp) - bz(ixp, iy)) / dyc(iy)
          jy1 = -(bz(ixp, iy) - bz(ix, iy)) / dxc(ix)
          jy2 = -(bz(ixp, iyp) - bz(ix, iyp)) / dxc(ix)
          jx = (jx1 + jx2) * 0.5_num
          jy = (jy1 + jy2) * 0.5_num
          jz = (by(ixp, iy) - by(ix, iy)) / dxc(ix) &
              - (bx(ix, iyp) - bx(ix, iy)) / dyc(iy)
          
          ! B at vertices
          bxv = (bx(ix, iy) + bx(ix, iyp)) * 0.5_num
          byv = (by(ix, iy) + by(ixp, iy)) * 0.5_num
          bzv = (bz(ix, iy) + bz(ixp, iy) &
              + bz(ix, iyp) + bz(ixp, iyp)) * 0.25_num

          magn_b = bxv**2 + byv**2 + bzv**2

          ! Calculate parallel and perpendicular currents
          j_par_x = (jx * bxv + jy * byv &
              + jz * bzv) * bxv / MAX(magn_b, none_zero)
          j_par_y = (jx * bxv + jy * byv &
              + jz * bzv) * byv / MAX(magn_b, none_zero)
          j_par_z = (jx * bxv + jy * byv &
              + jz * bzv) * bzv / MAX(magn_b, none_zero)

          ! If b = 0 then there is no parallel current
          IF (magn_b .LT. none_zero) THEN
            j_par_x = 0.0_num
            j_par_y = 0.0_num
            j_par_z = 0.0_num
          END IF

          ! Calculate perpendicular current
          j_perp_x = jx - j_par_x
          j_perp_y = jy - j_par_y
          j_perp_z = jz - j_par_z

          magn_j_par = SQRT(j_par_x**2 + j_par_y**2 + j_par_z**2)
          magn_j_perp = SQRT(j_perp_x**2 + j_perp_y**2 + j_perp_z**2)

          parallel_current(ix, iy) = magn_j_par
          perp_current(ix, iy) = magn_j_perp

          ! This isn't really curlb. It's actually heat flux
          curlb(ix, iy) = eta(ix, iy) * magn_j_par**2 &
              + (eta_perp(ix, iy) + eta(ix, iy)) * magn_j_perp**2

          flux_x(ix, iy) = -((j_par_x * eta(ix, iy) &
              + j_perp_x * (eta_perp(ix, iy) + eta(ix, iy)))&
              * dxc(ix) * 0.5_num)
          flux_y(ix, iy) = -((j_par_y * eta(ix, iy) &
              + j_perp_y * (eta_perp(ix, iy) + eta(ix, iy)))&
              * dyc(iy) * 0.5_num)
          flux_z(ix, iy) = -((j_par_z * eta(ix, iy) &
              + j_perp_z * (eta_perp(ix, iy) + eta(ix, iy))))
        END DO
      END DO
    END IF

  END SUBROUTINE rkstep



  ! Calculate the effect of the Hall term on the magnetic field
  ! Uses the subroutine rkstep1
  SUBROUTINE hall_effects

    REAL(num) :: dt6
    REAL(num), DIMENSION(:, :), ALLOCATABLE :: k1x, k2x, k3x, k4x
    REAL(num), DIMENSION(:, :), ALLOCATABLE :: k1y, k2y, k3y, k4y
    REAL(num), DIMENSION(:, :), ALLOCATABLE :: k1z, k2z, k3z, k4z

    ALLOCATE(k1x(0:nx, 0:ny), k2x(0:nx, 0:ny), k3x(0:nx, 0:ny), k4x(0:nx, 0:ny))
    ALLOCATE(k1y(0:nx, 0:ny), k2y(0:nx, 0:ny), k3y(0:nx, 0:ny), k4y(0:nx, 0:ny))
    ALLOCATE(k1z(0:nx, 0:ny), k2z(0:nx, 0:ny), k3z(0:nx, 0:ny), k4z(0:nx, 0:ny))

    dt = dt * 0.5_num

    bx1 = bx(0:nx+1, 0:ny+1)
    by1 = by(0:nx+1, 0:ny+1)
    bz1 = bz(0:nx+1, 0:ny+1)

    ! step 1
    CALL rkstep1
    k1x = flux_x
    k1y = flux_y
    k1z = flux_z

    ! step 2
    DO iy = 1, ny
      DO ix = 0, nx
        bx(ix, iy) = bx1(ix, iy) + (k1z(ix, iy) - k1z(ix, iy-1)) * dt / dyb(iy)
      END DO
    END DO

    DO iy = 0, ny
      DO ix = 1, nx
        by(ix, iy) = by1(ix, iy) + (-k1z(ix, iy) + k1z(ix-1, iy)) * dt / dxb(ix)
      END DO
    END DO

    DO iy = 1, ny
      DO ix = 1, nx
        ixm = ix - 1
        iym = iy - 1
        w1 = -k1x(ix, iy) - k1x(ixm, iy) + k1x(ixm, iym) + k1x(ix, iym)
        w2 = w1 + k1y(ix, iy) + k1y(ix, iym) - k1y(ixm, iy) - k1y(ixm, iym)
        bz(ix, iy) = bz1(ix, iy) + w2 * dt / cv(ix, iy)
      END DO
    END DO

    CALL bfield_bcs
    CALL rkstep1
    k2x = flux_x
    k2y = flux_y
    k2z = flux_z

    ! step 3
    DO iy = 1, ny
      DO ix = 0, nx
        bx(ix, iy) = bx1(ix, iy) + (k2z(ix, iy) - k2z(ix, iy-1)) * dt / dyb(iy)
      END DO
    END DO

    DO iy = 0, ny
      DO ix = 1, nx
        by(ix, iy) = by1(ix, iy) + (-k2z(ix, iy) + k2z(ix-1, iy)) * dt / dxb(ix)
      END DO
    END DO

    DO iy = 1, ny
      DO ix = 1, nx
        ixm = ix - 1
        iym = iy - 1
        w1 = -k2x(ix, iy) - k2x(ixm, iy) + k2x(ixm, iym) + k2x(ix, iym)
        w2 = w1 + k2y(ix, iy) + k2y(ix, iym) - k2y(ixm, iy) - k2y(ixm, iym)
        bz(ix, iy) = bz1(ix, iy) + w2 * dt / cv(ix, iy)
      END DO
    END DO

    CALL bfield_bcs
    CALL rkstep1
    k3x = flux_x
    k3y = flux_y
    k3z = flux_z

    dt = 2.0_num * dt

    ! step 4
    DO iy = 1, ny
      DO ix = 0, nx
        bx(ix, iy) = bx1(ix, iy) + (k3z(ix, iy) - k3z(ix, iy-1)) * dt / dyb(iy)
      END DO
    END DO

    DO iy = 0, ny
      DO ix = 1, nx
        by(ix, iy) = by1(ix, iy) + (-k3z(ix, iy) + k3z(ix-1, iy)) * dt / dxb(ix)
      END DO
    END DO

    DO iy = 1, ny
      DO ix = 1, nx
        ixm = ix - 1
        iym = iy - 1
        w1 = -k3x(ix, iy) - k3x(ixm, iy) + k3x(ixm, iym) + k3x(ix, iym)
        w2 = w1 + k3y(ix, iy) + k3y(ix, iym) - k3y(ixm, iy) - k3y(ixm, iym)
        bz(ix, iy) = bz1(ix, iy) + w2 * dt / cv(ix, iy)
      END DO
    END DO

    CALL bfield_bcs
    CALL rkstep1
    k4x = flux_x
    k4y = flux_y
    k4z = flux_z

    ! full update
    dt6 = dt * sixth
    k3x = k1x + 2.0_num * k2x + 2.0_num * k3x + k4x
    k3y = k1y + 2.0_num * k2y + 2.0_num * k3y + k4y
    k3z = k1z + 2.0_num * k2z + 2.0_num * k3z + k4z

    DO iy = 0, ny
      DO ix = 1, nx
        by(ix, iy) = by1(ix, iy) &
            + (-k3z(ix, iy) + k3z(ix-1, iy)) * dt6 / dxb(ix)
      END DO
    END DO

    DO iy = 1, ny
      DO ix = 0, nx
        bx(ix, iy) = bx1(ix, iy) &
            + (k3z(ix, iy) - k3z(ix, iy-1)) * dt6 / dyb(iy)
      END DO
    END DO

    DO iy = 1, ny
      DO ix = 1, nx
        ixm = ix - 1
        iym = iy - 1
        w1 = -k3x(ix, iy) - k3x(ixm, iy) + k3x(ixm, iym) + k3x(ix, iym)
        w2 = w1 + k3y(ix, iy) + k3y(ix, iym) - k3y(ixm, iy) - k3y(ixm, iym)
        bz(ix, iy) = bz1(ix, iy) + w2 * dt6 / cv(ix, iy)
      END DO
    END DO

    CALL bfield_bcs

    DEALLOCATE(k1x, k2x, k3x, k4x, k1y, k2y, k3y, k4y, k1z, k2z, k3z, k4z)

  END SUBROUTINE hall_effects



  SUBROUTINE rkstep1

    REAL(num), DIMENSION(:, :), ALLOCATABLE :: jx, jy, jz
    REAL(num), DIMENSION(:, :), ALLOCATABLE :: hxflux, hyflux, hzflux
    REAL(num) :: jx1, jy1, jx2, jy2
    REAL(num) :: bxv, byv, bzv, rho_v
    REAL(num) :: f1, f2, area
    REAL(num) :: jadp, jadm
    INTEGER :: ixp2, iyp2

    ALLOCATE(jx(-1:nx+1, -1:ny+1), jy(-1:nx+1, -1:ny+1), jz(-1:nx+1, -1:ny+1))
    ALLOCATE(hxflux(-1:nx+1, -1:ny+1), hyflux(-1:nx+1, -1:ny+1), &
        hzflux(-1:nx+1, -1:ny+1))

    DO iy = -1, ny + 1
      DO ix = -1, nx + 1
        ixp = ix + 1
        iyp = iy + 1

        jx1 = (bz(ix, iyp) - bz(ix, iy)) / dyc(iy)
        jx2 = (bz(ixp, iyp) - bz(ixp, iy)) / dyc(iy)
        jy1 = -(bz(ixp, iy) - bz(ix, iy)) / dxc(ix)
        jy2 = -(bz(ixp, iyp) - bz(ix, iyp)) / dxc(ix)
        jx(ix, iy) = (jx1 + jx2) * 0.5_num
        jy(ix, iy) = (jy1 + jy2) * 0.5_num
        jz(ix, iy) = (by(ixp, iy) - by(ix, iy)) / dxc(ix) &
            - (bx(ix, iyp) - bx(ix, iy)) / dyc(iy)
      END DO
    END DO

    hxflux = 0.0_num
    hyflux = 0.0_num
    hzflux = 0.0_num

    DO iy = 0, ny
      DO ix = 0, nx
        ixp = ix + 1
        ixp2 = ix + 2
        ixm = ix - 1
        iyp = iy + 1
        iyp2 = iy + 2
        iym = iy - 1

        bxv = (bx(ix, iy) * dyb(iy) + bx(ix, iyp) * dyb(iyp)) &
            / (dyb(iy) + dyb(iyp))
        byv = (by(ix, iy) * dxb(ix) + by(ixp, iy) * dxb(ixp)) &
            / (dxb(ix) + dxb(ixp))
        bzv = (bz(ix, iy) * cv(ix, iy) + bz(ixp, iy) * cv(ixp, iy) &
            + bz(ix, iyp) * cv(ix, iyp) + bz(ixp, iyp) * cv(ixp, iyp))

        rho_v = (rho(ix, iy) * cv(ix, iy) + rho(ixp, iy) * cv(ixp, iy) &
            + rho(ix, iyp) * cv(ix, iyp) + rho(ixp, iyp) * cv(ixp, iyp))

        area = (cv(ix, iy) + cv(ixp, iy) + cv(ix, iyp) + cv(ixp, iyp))
        bzv = bzv / area
        rho_v = rho_v / area

        ! Evaluate the flux due to the Hall term but upwind in direction of
        ! negative current density i.e. upwind in the direction of the
        ! effective advection velocity
        f1 = MAX(0.0_num, jy(ix, iy)) * (bz(ix, iyp) + bz(ixp, iyp)) * 0.5_num &
            + MIN(0.0_num, jy(ix, iy)) * (bz(ix, iy) + bz(ixp, iy)) * 0.5_num

        w1 = bz(ix, iyp) - bz(ix, iy) + bz(ixp, iyp) - bz(ixp, iy)
        w2 = bz(ix, iy) - bz(ix, iym) + bz(ixp, iy) - bz(ixp, iym)
        w3 = bz(ix, iyp2) - bz(ix, iyp) + bz(ixp, iyp2) - bz(ixp, iyp)
        w1 = 0.5_num * w1
        w2 = 0.5_num * w2
        w3 = 0.5_num * w3

        jadm = -MIN(SIGN(1.0_num, jy(ix, iy)), 0.0_num)
        jadp = 1.0_num - jadm

        w5 = lambda_i(ix, iy) * ABS(jy(ix, iy)) * dt / dyb(iy) / rho_v * jadm &
            + lambda_i(ix, iy) * ABS(jy(ix, iy)) * dt / dyb(iyp) / rho_v * jadp
        w4 = ((2.0_num - w5) * ABS(w1) / dyc(iy) &
            + (1.0_num + w5) * ABS(w2) / dyc(iym)) * jadm + &
            ((2.0_num - w5) * ABS(w1) / dyc(iy) &
            + (1.0_num + w5) * ABS(w3) / dyc(iyp)) * jadp

        w4 = w4 * sixth
        w8 = 0.5_num * (SIGN(1.0_num, w1) &
            + SIGN(1.0_num, w2 * jadm + w3 * jadp))

        w6 = -SIGN(1.0_num, jy(ix, iy)) * w8 &
            * MIN(ABS(w4) * (dyb(iy) * jadm + dyb(iyp) * jadp), &
            ABS(w1), ABS(w2 * jadm + w3 * jadp))

        f1 = f1 + jy(ix, iy) * w6 * (1.0_num - w5)
        f2 = jz(ix, iy) * byv
        hxflux(ix, iy) = lambda_i(ix, iy) * (f1 - f2) / rho_v

        f1 = jz(ix, iy) * bxv
        f2 = MAX(0.0_num, jx(ix, iy)) * (bz(ixp, iyp) + bz(ixp, iy)) * 0.5_num &
            + MIN(0.0_num, jx(ix, iy)) * (bz(ix, iy) + bz(ix, iyp)) * 0.5_num

        w1 = bz(ixp, iy) - bz(ix, iy) + bz(ixp, iyp) - bz(ix, iyp)
        w2 = bz(ix, iy) - bz(ixm, iy) + bz(ix, iyp) - bz(ixm, iyp)
        w3 = bz(ixp2, iy) - bz(ixp, iy) + bz(ixp2, iyp) - bz(ixp, iyp)
        w1 = 0.5_num * w1
        w2 = 0.5_num * w2
        w3 = 0.5_num * w3

        jadm = -MIN(SIGN(1.0_num, jy(ix, iy)), 0.0_num)
        jadp = 1.0_num - jadm

        w5 = lambda_i(ix, iy) * ABS(jx(ix, iy)) * dt / &
            (dxb(ix) * jadm + dxb(ixp) * jadp) / rho_v
        w4 = (2.0_num - w5) * ABS(w1) / dxc(ix) &
            + (1.0_num + w5) * ABS(w2 * jadm + w3 * jadp) &
            / (dxc(ixm) * jadm + dxc(ixp) * jadp)

        w4 = w4 * sixth
        w8 = 0.5_num * (SIGN(1.0_num, w1) &
            + SIGN(1.0_num, w2 * jadm + w3 * jadm))
        w6 = -SIGN(1.0_num, jy(ix, iy)) * w8 &
            * MIN(ABS(w4) * (dxb(ix) * jadm + dxb(ixp) * jadp), &
            ABS(w1), ABS(w2 * jadm + w3 * jadp))

        f2 = f2 + jx(ix, iy) * w6 * (1.0_num - w5)

        hyflux(ix, iy) = lambda_i(ix, iy) * (f1 - f2) / rho_v

        f1 = MAX(0.0_num, jx(ix, iy)) * by(ixp, iy) &
            + MIN(0.0_num, jx(ix, iy)) * by(ix, iy)

        w1 = by(ixp, iy) - by(ix, iy)
        w2 = by(ix, iy) - by(ixm, iy)
        w3 = by(ixp2, iy) - by(ixp, iy)

        jadm = -MIN(SIGN(1.0_num, jy(ix, iy)), 0.0_num)
        jadp = 1.0_num - jadm

        w5 = lambda_i(ix, iy) * ABS(jx(ix, iy)) * dt &
            / (dxb(ix) * jadm + dxb(ixp) * jadp) / rho_v

        w4 = (2.0_num - w5) * ABS(w1) / dxc(ix) &
            + (1.0_num + w5) * (ABS(w2) / dxc(ixm) * jadm + &
            ABS(w3) / dxc(ixp) * jadp)

        w4 = w4 * sixth
        w8 = 0.5_num * (SIGN(1.0_num, w1) &
            + SIGN(1.0_num, w2 * jadm + w3 * jadp))

        w6 = -SIGN(1.0_num, jy(ix, iy)) * w8 &
            * MIN(ABS(w4) * (dxb(ix) * jadm + dxb(ixp) * jadp), &
            ABS(w1), ABS(w2 * jadm + w3 * jadp))

        f1 = f1 + jx(ix, iy) * w6 * (1.0_num - w5)

        f2 = MAX(0.0_num, jy(ix, iy)) * bx(ix, iyp) &
            + MIN(0.0_num, jy(ix, iy)) * bx(ix, iy)

        w1 = bx(ix, iyp) - bx(ix, iy)
        w2 = bx(ix, iy) - bx(ix, iym)
        w3 = bx(ix, iyp2) - bx(ix, iyp)

        jadm = -MIN(SIGN(1.0_num, jy(ix, iy)), 0.0_num)
        jadp = 1.0_num - jadm

        w5 = lambda_i(ix, iy) * ABS(jy(ix, iy)) * dt &
            / (dyb(iy) * jadm + dyb(iyp) * jadp) / rho_v

        w4 = (2.0_num - w5) * ABS(w1) / dyc(iy) &
            + (1.0_num + w5) * (ABS(w2) / dyc(iym) * jadm + &
            ABS(w3) / dyc(iy) * jadp)

        w4 = w4 * sixth
        w8 = 0.5_num * (SIGN(1.0_num, w1) &
            + SIGN(1.0_num, w2 * jadm + w3 * jadp))

        w6 = -SIGN(1.0_num, jy(ix, iy)) * w8 * MIN(ABS(w4) * dyb(iy), &
            ABS(w1), ABS(w2 * jadm + w3 * jadp))
        f2 = f2 + jy(ix, iy) * w6 * (1.0_num - w5)

        hzflux(ix, iy) = lambda_i(ix, iy) * (f1 - f2) / rho_v

      END DO
    END DO

    DO iy = 0, ny
      DO ix = 0, nx
        flux_x(ix, iy) = -hxflux(ix, iy) * dxc(ix) * 0.5_num
        flux_y(ix, iy) = -hyflux(ix, iy) * dyc(iy) * 0.5_num
        flux_z(ix, iy) = -hzflux(ix, iy)
      END DO
    END DO

    DEALLOCATE (jx, jy, jz, hxflux, hyflux, hzflux)

  END SUBROUTINE rkstep1

END MODULE lagran
