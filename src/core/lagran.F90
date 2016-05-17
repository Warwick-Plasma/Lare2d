!******************************************************************************
! Lagrangian step routines
!******************************************************************************

MODULE lagran

  USE shared_data
  USE boundary
  USE neutral
  USE conduct

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: lagrangian_step, eta_calc

  ! Only used inside lagran.f90
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: alpha1, alpha2
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: visc_heat, pressure, rho_v, cv_v
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: flux_x, flux_y, flux_z, curlb
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: fx_visc, fy_visc, fz_visc
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: fx, fy, fz

CONTAINS

  !****************************************************************************
  ! This subroutine manages the progress of the lagrangian step
  !****************************************************************************

  SUBROUTINE lagrangian_step

    INTEGER :: substeps, subcycle
    REAL(num) :: actual_dt, dt_sub

    ALLOCATE(bx1(-1:nx+2,-1:ny+2))
    ALLOCATE(by1(-1:nx+2,-1:ny+2))
    ALLOCATE(bz1(-1:nx+2,-1:ny+2))
    ALLOCATE(alpha1(0:nx+1,0:ny+2))
    ALLOCATE(alpha2(-1:nx+1,0:ny+1))
    ALLOCATE(visc_heat(0:nx+1,0:ny+1))
    ALLOCATE(pressure(-1:nx+2,-1:ny+2))
    ALLOCATE(rho_v(-1:nx+1,-1:ny+1))
    ALLOCATE(cv_v(-1:nx+1,-1:ny+1))
    ALLOCATE(fx(0:nx,0:ny))
    ALLOCATE(fy(0:nx,0:ny))
    ALLOCATE(fz(0:nx,0:ny))
    ALLOCATE(fx_visc(0:nx,0:ny))
    ALLOCATE(fy_visc(0:nx,0:ny))
    ALLOCATE(fz_visc(0:nx,0:ny))
    ALLOCATE(flux_x(0:nx,0:ny))
    ALLOCATE(flux_y(0:nx,0:ny))
    ALLOCATE(flux_z(0:nx,0:ny))
    ALLOCATE(curlb (0:nx,0:ny))

    DO iy = -1, ny + 2
      iym = iy - 1
      DO ix = -1, nx + 2
        ixm = ix - 1
        bx1(ix,iy) = (bx(ix,iy) + bx(ixm,iy )) * 0.5_num
        by1(ix,iy) = (by(ix,iy) + by(ix ,iym)) * 0.5_num
        bz1(ix,iy) = bz(ix,iy)

        pressure(ix,iy) = (gamma - 1.0_num) * rho(ix,iy) &
            * (energy(ix,iy) - (1.0_num - xi_n(ix,iy)) * ionise_pot)
      END DO
    END DO

    DO iy = -1, ny + 1
      iyp = iy + 1
      DO ix = -1, nx + 1
        ixp = ix + 1
        rho_v(ix,iy) = rho(ix,iy) * cv(ix,iy) + rho(ixp,iy) * cv(ixp,iy) &
            +   rho(ix,iyp) * cv(ix,iyp) + rho(ixp,iyp) * cv(ixp,iyp)
        cv_v(ix,iy) = cv(ix,iy) + cv(ixp,iy) + cv(ix,iyp) + cv(ixp,iyp)
        rho_v(ix,iy) = rho_v(ix,iy) / cv_v(ix,iy)

        cv_v(ix,iy) = 0.25_num * cv_v(ix,iy)
      END DO
    END DO

    CALL shock_viscosity
    CALL set_dt
    dt2 = dt * 0.5_num

    IF (resistive_mhd .OR. hall_mhd) THEN
      ! If subcycling isn't wanted set dt = dtr in set_dt, don't just
      ! set substeps to 1.
      IF (resistive_mhd) THEN
        dt_sub = dtr
        IF (hall_mhd) dt_sub = MIN(dtr, dth)
      ELSE
        dt_sub = dth
      END IF

      substeps = INT(dt / dt_sub) + 1

      IF (substeps > peak_substeps) peak_substeps = substeps
      actual_dt = dt
      dt = dt / REAL(substeps, num)

      DO subcycle = 1, substeps
        CALL eta_calc
        IF (eos_number /= EOS_IDEAL) CALL neutral_fraction
        IF (cowling_resistivity) CALL perpendicular_resistivity
        IF (hall_mhd) CALL hall_effects
        IF (resistive_mhd) CALL resistive_effects
      END DO

      dt = actual_dt
    END IF

    IF (conduction .OR. coronal_heating .OR. radiation) CALL conduct_heat

    CALL predictor_corrector_step

    DEALLOCATE(bx1, by1, bz1, alpha1, alpha2)
    DEALLOCATE(visc_heat, pressure, rho_v, cv_v, flux_x, flux_y, flux_z, curlb)
    DEALLOCATE(fx_visc, fy_visc, fz_visc)
    DEALLOCATE(fx, fy, fz)
    DEALLOCATE(flux_x, flux_y, flux_z)
    DEALLOCATE(curlb)

    CALL energy_bcs
    CALL density_bcs
    CALL velocity_bcs

  END SUBROUTINE lagrangian_step



  !****************************************************************************
  ! The main predictor / corrector step which advances the momentum equation
  !****************************************************************************

  SUBROUTINE predictor_corrector_step

    REAL(num) :: pp, ppx, ppy, ppxy
    REAL(num) :: e1
    REAL(num) :: vxb, vxbm, vyb, vybm
    REAL(num) :: bxv, byv, bzv, jx, jy, jz
    REAL(num) :: cvx, cvxp, cvy, cvyp
    REAL(num) :: dv

    CALL b_field_and_cv1_update

    bx1 = bx1 * cv1
    by1 = by1 * cv1
    bz1 = bz1 * cv1

    DO iy = 0, ny + 1
      DO ix = 0, nx + 1
        dv = cv1(ix,iy) / cv(ix,iy) - 1.0_num
        ! Predictor energy
        e1 = energy(ix,iy) - pressure(ix,iy) * dv / rho(ix,iy)
        e1 = e1 + visc_heat(ix,iy) * dt2 / rho(ix,iy)

        ! Now define the predictor step pressures
        pressure(ix,iy) = (e1 - (1.0_num - xi_n(ix,iy)) * ionise_pot) &
            * (gamma - 1.0_num) * rho(ix,iy) * cv(ix,iy) / cv1(ix,iy)
      END DO
    END DO

    DO iy = 0, ny
      iyp = iy + 1
      iym = iy - 1
      DO ix = 0, nx
        ixp = ix + 1
        ixm = ix - 1

        pp    = pressure(ix ,iy )
        ppx   = pressure(ixp,iy )
        ppy   = pressure(ix ,iyp)
        ppxy  = pressure(ixp,iyp)

        ! P total at Ex(i,j)
        w1 = (pp + ppy) * 0.5_num
        ! P total at Ex(i+1,j)
        w2 = (ppx + ppxy) * 0.5_num
        fx(ix,iy) = -(w2 - w1) / dxc(ix)

        ! P total at Ey(i,j)
        w1 = (pp + ppx) * 0.5_num
        ! P total at Ey(i,j+1)
        w2 = (ppy + ppxy) * 0.5_num
        fy(ix,iy) = -(w2 - w1) / dyc(iy)

        fz(ix,iy) = 0.0_num

        cvx  = cv1(ix ,iy ) + cv1(ix ,iyp)
        cvxp = cv1(ixp,iy ) + cv1(ixp,iyp)
        cvy  = cv1(ix ,iy ) + cv1(ixp,iy )
        cvyp = cv1(ix ,iyp) + cv1(ixp,iyp)

        w1 = (bz1(ix ,iy ) + bz1(ixp,iy )) / cvy
        w2 = (bz1(ix ,iyp) + bz1(ixp,iyp)) / cvyp
        jx = (w2 - w1) / dyc(iy)

        w1 = (bz1(ix ,iy ) + bz1(ix ,iyp)) / cvx
        w2 = (bz1(ixp,iy ) + bz1(ixp,iyp)) / cvxp
        jy = -(w2 - w1) / dxc(ix)

        w1 = (by1(ix ,iy ) + by1(ix ,iyp)) / cvx
        w2 = (by1(ixp,iy ) + by1(ixp,iyp)) / cvxp
        jz = (w2 - w1) / dxc(ix)

        w1 = (bx1(ix ,iy ) + bx1(ixp,iy )) / cvy
        w2 = (bx1(ix ,iyp) + bx1(ixp,iyp)) / cvyp
        jz = jz - (w2 - w1) / dyc(iy)

        bxv = (bx1(ix,iy ) + bx1(ixp,iy ) + bx1(ix,iyp) + bx1(ixp,iyp)) &
            / (cvx + cvxp)

        byv = (by1(ix,iy ) + by1(ixp,iy ) + by1(ix,iyp) + by1(ixp,iyp)) &
            / (cvx + cvxp)

        bzv = (bz1(ix,iy ) + bz1(ixp,iy ) + bz1(ix,iyp) + bz1(ixp,iyp)) &
            / (cvx + cvxp)

        fx(ix,iy) = fx(ix,iy) + (jy * bzv - jz * byv)
        fy(ix,iy) = fy(ix,iy) + (jz * bxv - jx * bzv)
        fz(ix,iy) = fz(ix,iy) + (jx * byv - jy * bxv)

        fy(ix,iy) = fy(ix,iy) - rho_v(ix,iy) * grav(iy)

        ! Find half step velocity needed for remap
        vx1(ix,iy) = vx(ix,iy) + dt2 * (fx_visc(ix,iy) + fx(ix,iy)) / rho_v(ix,iy)
        vy1(ix,iy) = vy(ix,iy) + dt2 * (fy_visc(ix,iy) + fy(ix,iy)) / rho_v(ix,iy)
        vz1(ix,iy) = vz(ix,iy) + dt2 * (fz_visc(ix,iy) + fz(ix,iy)) / rho_v(ix,iy)
      END DO
    END DO

    CALL remap_v_bcs

    bx1 = bx1 / cv1
    by1 = by1 / cv1
    bz1 = bz1 / cv1
    CALL shock_heating

    DO iy = 0, ny
      DO ix = 0, nx
        ! Velocity at the end of the Lagrangian step
        vx(ix,iy) = vx(ix,iy) + dt * (fx_visc(ix,iy) + fx(ix,iy)) / rho_v(ix,iy)
        vy(ix,iy) = vy(ix,iy) + dt * (fy_visc(ix,iy) + fy(ix,iy)) / rho_v(ix,iy)
        vz(ix,iy) = vz(ix,iy) + dt * (fz_visc(ix,iy) + fz(ix,iy)) / rho_v(ix,iy)
      END DO
    END DO
    
    CALL velocity_bcs

    ! Finally correct density and energy to final values
    DO iy = 1, ny
      iym = iy - 1
      DO ix = 1, nx
        ixm = ix - 1

        ! vx1 at Bx(i,j)
        vxb  = (vx1(ix ,iy ) + vx1(ix ,iym)) * 0.5_num
        ! vx1 at Bx(i-1,j)
        vxbm = (vx1(ixm,iy ) + vx1(ixm,iym)) * 0.5_num
        ! vy1 at By(i,j)
        vyb  = (vy1(ix ,iy ) + vy1(ixm,iy )) * 0.5_num
        ! vy1 at By(i,j-1)
        vybm = (vy1(ix ,iym) + vy1(ixm,iym)) * 0.5_num

        dv = ((vxb - vxbm) / dxb(ix) + (vyb - vybm) / dyb(iy)) * dt

        cv1(ix,iy) = cv(ix,iy) * (1.0_num + dv)

        ! Energy at end of Lagrangian step
        energy(ix,iy) = energy(ix,iy) &
            + (dt * visc_heat(ix,iy) - dv * pressure(ix,iy)) &
            / rho(ix,iy)

        rho(ix,iy) = rho(ix,iy) / (1.0_num + dv)

        total_visc_heating = total_visc_heating &
            + dt * visc_heat(ix,iy) * cv(ix,iy)

      END DO
    END DO

  END SUBROUTINE predictor_corrector_step



  !****************************************************************************
  ! This subroutine calculates the viscous effects and updates the
  ! magnetic field
  !****************************************************************************

  SUBROUTINE shock_viscosity

    REAL(num) :: dvdots, dx, dxm, dxp
    REAL(num) :: b2, rmin
    REAL(num) :: a1, a2, a3, a4
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: cs, cs_v
    INTEGER :: i0, i1, i2, i3, j0, j1, j2, j3
    LOGICAL, SAVE :: first_call = .TRUE.

    ALLOCATE(cs(-1:nx+2,-1:ny+2), cs_v(-1:nx+1,-1:ny+1))

    IF (first_call) THEN
      first_call = .FALSE.
      visc2_norm = 0.25_num * (gamma + 1.0_num) * visc2
    END IF

    p_visc = 0.0_num
    visc_heat = 0.0_num

    DO ix = -1, nx + 2
      DO iy = -1, ny + 2
        rmin = MAX(rho(ix,iy), none_zero)
        b2 = bx1(ix,iy)**2 + by1(ix,iy)**2 + bz1(ix,iy)**2
        cs(ix,iy) = SQRT((gamma * pressure(ix,iy) + b2) / rmin)
      END DO
    END DO

    DO iy = -1, ny + 1
      iyp = iy + 1
      DO ix = -1, nx + 1
        ixp = ix + 1

        cs_v(ix,iy) = cs(ix,iy) * cv(ix,iy) + cs(ixp,iy) * cv(ixp,iy) &
            +   cs(ix,iyp) * cv(ix,iyp) + cs(ixp,iyp) * cv(ixp,iyp)
        cs_v(ix,iy) = 0.25_num * cs_v(ix,iy) / cv_v(ix,iy)

      END DO
    END DO

    DO iy = 0, ny + 2
      iym = iy - 1
      iyp = iy + 1
      DO ix = 0, nx + 1
        ixm = ix - 1
        ixp = ix + 1

        ! Edge viscosities from Caramana
        ! Triangles numbered as in Goffrey thesis

        ! Edge viscosity for triangle 1
        i1 = ixm
        j1 = iym
        i2 = ix
        j2 = iym
        i0 = i1 - 1
        j0 = j1
        i3 = i2 + 1
        j3 = j1
        dx = dxb(ix)
        dxp = dxb(ixp)
        dxm = dxb(ixm)
        ! dv in direction of dS, i.e. dv.dS / abs(dS)
        dvdots = - (vx(i1,j1) - vx(i2,j2))
        ! Force on node is alpha*dv*ds but store only alpha and convert to force
        ! when needed.  
        alpha1(ix,iy) = edge_viscosity()
      END DO
    END DO

      ! Edge viscosity for triangle 2
    DO iy = 0, ny + 1
      iym = iy - 1
      iyp = iy + 1
      DO ix = -1, nx + 1
        ixm = ix - 1
        ixp = ix + 1

        i1 = ix
        j1 = iym
        i2 = ix
        j2 = iy
        i0 = i1
        j0 = j1 - 1
        i3 = i2
        j3 = j1 + 1
        dx = dyb(iy)
        dxp = dyb(iyp)
        dxm = dyb(iym)
        dvdots = - (vy(i1,j1) - vy(i2,j2))
        alpha2(ix,iy) = edge_viscosity()
      END DO
    END DO

    DO iy = 0, ny + 1 
      iym = iy - 1
      iyp = iy + 1
      DO ix = 0, nx + 1 
        ixm = ix - 1
        ixp = ix + 1
        ! Estimate p_visc based on alpha * dv, for timestep control
        a1 = ((vx(ixm,iym) - vx(ix ,iym))**2  &
              + (vy(ixm,iym) - vy(ix ,iym))**2 + (vz(ixm,iym) - vz(ix ,iym))**2) 
        a2 = ((vx(ix ,iym) - vx(ix ,iy ))**2  &
              + (vy(ix ,iym) - vy(ix ,iy ))**2 + (vz(ix ,iym) - vz(ix ,iy ))**2)
        a3 = ((vx(ix ,iy ) - vx(ixm,iy ))**2  &
              + (vy(ix ,iy ) - vy(ixm,iy ))**2 + (vz(ix ,iy ) - vz(ixm,iy ))**2) 
        a4 = ((vx(ixm,iy ) - vx(ixm,iym))**2  &
              + (vy(ixm,iy ) - vy(ixm,iym))**2 + (vz(ixm,iy ) - vz(ixm,iym))**2)

        p_visc(ix,iy) = MAX(p_visc(ix,iy), - alpha1(ix,iy)*SQRT(a1)) 
        p_visc(ix,iy) = MAX(p_visc(ix,iy), - alpha2(ix,iy)*SQRT(a2)) 

        visc_heat(ix,iy) = &
            - 0.5_num * dyb(iy) * alpha1(ix ,iy ) * a1 &
            - 0.5_num * dxb(ix) * alpha2(ix ,iy ) * a2 &
            - 0.5_num * dyb(iy) * alpha1(ix ,iyp) * a3 &
            - 0.5_num * dxb(ix) * alpha2(ixm,iy ) * a4

        visc_heat(ix,iy) = visc_heat(ix,iy) / cv(ix,iy)
      END DO
    END DO

    fx_visc = 0.0_num
    fy_visc = 0.0_num
    fz_visc = 0.0_num
    DO iy = 0, ny
      iym = iy - 1
      iyp = iy + 1
      DO ix = 0, nx
        ixm = ix - 1
        ixp = ix + 1

        a1 = alpha1(ix ,iyp) * dyc(iy)
        a2 = alpha1(ixp,iyp) * dyc(iy)
        a3 = alpha2(ix ,iy ) * dxc(ix)
        a4 = alpha2(ix ,iyp) * dxc(ix)

        fx_visc(ix,iy) = (a1 * (vx(ix,iy) - vx(ixm,iy )) &
                        + a2 * (vx(ix,iy) - vx(ixp,iy )) &
                        + a3 * (vx(ix,iy) - vx(ix ,iym)) &
                        + a4 * (vx(ix,iy) - vx(ix ,iyp)) ) / cv_v(ix,iy)

        fy_visc(ix,iy) = (a1 * (vy(ix,iy) - vy(ixm,iy )) &
                        + a2 * (vy(ix,iy) - vy(ixp,iy )) &
                        + a3 * (vy(ix,iy) - vy(ix ,iym)) &
                        + a4 * (vy(ix,iy) - vy(ix ,iyp)) ) / cv_v(ix,iy)

        fz_visc(ix,iy) = (a1 * (vz(ix,iy) - vz(ixm,iy )) &
                        + a2 * (vz(ix,iy) - vz(ixp,iy )) &
                        + a3 * (vz(ix,iy) - vz(ix ,iym)) &
                        + a4 * (vz(ix,iy) - vz(ix ,iyp)) ) / cv_v(ix,iy)

      END DO
    END DO

    DEALLOCATE(cs, cs_v)

  CONTAINS

    DOUBLE PRECISION FUNCTION edge_viscosity()

      ! Actually returns q_k_bar = q_kur*(1-psi) / abs(dv)
      ! Other symbols follow notation in Caramana

      REAL(num) :: dvx, dvy, dvz, dv, dv2
#ifdef SHOCKLIMITER      
      REAL(num) :: dvxm, dvxp, dvym, dvyp, dvzm, dvzp
      REAL(num) :: rl, rr
#endif
      REAL(num) :: psi, rho_edge, cs_edge, q_k_bar

#ifdef SHOCKLIMITER 
      ! Turn off shock viscosity if cell edge expanding
      dvdots = MIN(0.0_num, dvdots)
#else
      ! Allow shock viscoity on expanding edge
      dvdots = -ABS(dvdots)
#endif

      rho_edge = 2.0_num * rho_v(i1,j1) * rho_v(i2,j2) &
          / (rho_v(i1,j1) + rho_v(i2,j2))
      cs_edge = MIN(cs_v(i1,j1), cs_v(i2,j2))

      dvx = vx(i1,j1) - vx(i2,j2)
      dvy = vy(i1,j1) - vy(i2,j2)
      dvz = vz(i1,j1) - vz(i2,j2)
      dv2 = dvx**2 + dvy**2 + dvz**2
      dv = SQRT(dv2)
      psi = 0.0_num
      IF (dv * dt / dx < 1.e-14_num) THEN
        dvdots = 0.0_num
      ELSE
        dvdots = dvdots / dv
      END IF

#ifdef SHOCKLIMITER
      dvxm = vx(i0,j0) - vx(i1,j1)
      dvxp = vx(i2,j2) - vx(i3,j3)
      dvym = vy(i0,j0) - vy(i1,j1)
      dvyp = vy(i2,j2) - vy(i3,j3)
      dvzm = vz(i0,j0) - vz(i1,j1)
      dvzp = vz(i2,j2) - vz(i3,j3)
      IF (dv * dt / dx < 1.e-14_num) THEN
        rl = 1.0_num
        rr = 1.0_num
      ELSE
        rl = (dvxp * dvx + dvyp * dvy + dvzp * dvz) * dx / (dxp * dv2)
        rr = (dvxm * dvx + dvym * dvy + dvzm * dvz) * dx / (dxm * dv2)
      END IF
      psi = MIN(0.5_num * (rr + rl), 2.0_num * rl, 2.0_num * rr, 1.0_num)
      psi = MAX(0.0_num, psi)
#endif

      ! Find q_kur / abs(dv)
      q_k_bar = rho_edge &
          * (visc2_norm * dv + SQRT(visc2_norm**2 * dv2 + (visc1 * cs_edge)**2))

      edge_viscosity = q_k_bar * (1.0_num - psi) * dvdots

    END FUNCTION edge_viscosity

  END SUBROUTINE shock_viscosity



  SUBROUTINE shock_heating

    REAL(num) :: a1, a2, a3, a4

    visc_heat = 0.0_num

    DO iy = 0, ny + 1 
      iym = iy - 1
      iyp = iy + 1
      DO ix = 0, nx + 1 
        ixm = ix - 1
        ixp = ix + 1

        a1 =  (vx(ixm,iym) - vx(ix ,iym))*(vx1(ixm,iym) - vx1(ix ,iym)) &
            + (vy(ixm,iym) - vy(ix ,iym))*(vy1(ixm,iym) - vy1(ix ,iym)) &
            + (vz(ixm,iym) - vz(ix ,iym))*(vz1(ixm,iym) - vz1(ix ,iym)) 
        a2 =  (vx(ix ,iym) - vx(ix ,iy ))*(vx1(ix ,iym) - vx1(ix ,iy )) &
            + (vy(ix ,iym) - vy(ix ,iy ))*(vy1(ix ,iym) - vy1(ix ,iy )) &
            + (vz(ix ,iym) - vz(ix ,iy ))*(vz1(ix ,iym) - vz1(ix ,iy ))
        a3 =  (vx(ix ,iy ) - vx(ixm,iy ))*(vx1(ix ,iy ) - vx1(ixm,iy )) &
            + (vy(ix ,iy ) - vy(ixm,iy ))*(vy1(ix ,iy ) - vy1(ixm,iy )) &
            + (vz(ix ,iy ) - vz(ixm,iy ))*(vz1(ix ,iy ) - vz1(ixm,iy ))
        a4 =  (vx(ixm,iy ) - vx(ixm,iym))*(vx1(ixm,iy ) - vx1(ixm,iym)) &
            + (vy(ixm,iy ) - vy(ixm,iym))*(vy1(ixm,iy ) - vy1(ixm,iym)) &
            + (vz(ixm,iy ) - vz(ixm,iym))*(vz1(ixm,iy ) - vz1(ixm,iym))

        visc_heat(ix,iy) = &
            - 0.5_num * dyb(iy) * alpha1(ix,iy) * a1 &
            - 0.5_num * dxb(ix) * alpha2(ix,iy) * a2 &
            - 0.5_num * dyb(iy) * alpha1(ix,iyp) * a3 &
            - 0.5_num * dxb(ix) * alpha2(ixm,iy) * a4

        visc_heat(ix,iy) = visc_heat(ix,iy) / cv(ix,iy)
      END DO
    END DO

    visc_heat = MAX(visc_heat, 0.0_num)

  END SUBROUTINE shock_heating



  SUBROUTINE b_field_and_cv1_update

    REAL(num) :: vxb, vxbm, vyb, vybm, vzb, vzbm
    REAL(num) :: dvxdx, dvydx, dvzdx
    REAL(num) :: dvxdy, dvydy, dvzdy
    REAL(num) :: dv

    DO iy = -1, ny + 2
      iym = iy - 1
      DO ix = -1, nx + 2
        ixm = ix - 1

        ! vx at Bx(i,j)
        vxb  = (vx(ix ,iy ) + vx(ix ,iym)) * 0.5_num
        ! vx at Bx(i-1,j)
        vxbm = (vx(ixm,iy ) + vx(ixm,iym)) * 0.5_num
        ! vy at By(i,j)
        vyb  = (vy(ix ,iy ) + vy(ixm,iy )) * 0.5_num
        ! vy at By(i,j-1)
        vybm = (vy(ix ,iym) + vy(ixm,iym)) * 0.5_num

        dvxdx = (vxb - vxbm) / dxb(ix)
        dvydy = (vyb - vybm) / dyb(iy)

        dv = (dvxdx + dvydy) * dt2
        cv1(ix,iy) = cv(ix,iy) * (1.0_num + dv)

        ! vx at By(i,j)
        vxb  = (vx(ix ,iy ) + vx(ixm,iy )) * 0.5_num
        ! vx at By(i,j-1)
        vxbm = (vx(ix ,iym) + vx(ixm,iym)) * 0.5_num
        ! vy at Bx(i,j)
        vyb  = (vy(ix ,iy ) + vy(ix ,iym)) * 0.5_num
        ! vy at Bx(i-1,j)
        vybm = (vy(ixm,iy ) + vy(ixm,iym)) * 0.5_num

        dvxdy = (vxb - vxbm) / dyb(iy)
        dvydx = (vyb - vybm) / dxb(ix)

        ! vz at Bx(i,j)
        vzb  = (vz(ix ,iy ) + vz(ix ,iym)) * 0.5_num
        ! vz at Bx(i-1,j)
        vzbm = (vz(ixm,iy ) + vz(ixm,iym)) * 0.5_num
        ! vz at By(i,j)
        dvzdx = (vzb - vzbm) / dxb(ix)

        vzb  = (vz(ix ,iy ) + vz(ixm,iy )) * 0.5_num
        ! vz at By(i,j-1)
        vzbm = (vz(ix ,iym) + vz(ixm,iym)) * 0.5_num
        dvzdy = (vzb - vzbm) / dyb(iy)

        w3 =  bx1(ix,iy) * dvxdx + by1(ix,iy) * dvxdy
        w4 =  bx1(ix,iy) * dvydx + by1(ix,iy) * dvydy
        w5 =  bx1(ix,iy) * dvzdx + by1(ix,iy) * dvzdy

        bx1(ix,iy) = (bx1(ix,iy) + w3 * dt2) / (1.0_num + dv)
        by1(ix,iy) = (by1(ix,iy) + w4 * dt2) / (1.0_num + dv)
        bz1(ix,iy) = (bz1(ix,iy) + w5 * dt2) / (1.0_num + dv)
      END DO
    END DO

  END SUBROUTINE b_field_and_cv1_update



  !****************************************************************************
  ! Sets CFL limited step
  !****************************************************************************

  SUBROUTINE set_dt

    ! Assumes all variables are defined at the same point. Be careful with
    ! setting 'dt_multiplier' if you expect massive changes across cells.

    REAL(num) :: vxbm, vxbp, avxm, avxp, dvx, ax
    REAL(num) :: vybm, vybp, avym, avyp, dvy, ay
    REAL(num) :: cs2, c_visc2, area, rho0, length
    REAL(num) :: dxlocal, dt_local, dtr_local, dt1, dt2, dth_local
    REAL(num) :: dt_locals(3), dt_min(3)
    REAL(num) :: dt0, time_dump, time_rem
    REAL(num) :: dt_fudge = 1e-4_num
    CHARACTER(LEN=1) :: dt_reason
    LOGICAL :: is_restart = .FALSE.
    LOGICAL, SAVE :: first = .TRUE.

    IF (first) THEN
      first = .FALSE.
      IF (restart) THEN
        dt = dt_from_restart
        RETURN
      END IF
    END IF

    dt_local = largest_number
    dtr_local = largest_number
    dth_local = largest_number

    DO iy = 0, ny
      iym = iy - 1
      DO ix = 0, nx
        ixm = ix - 1

        ! Fix dt for Lagrangian step
        w1 = bx1(ix,iy)**2 + by1(ix,iy)**2 + bz1(ix,iy)**2
        ! Sound speed squared
        rho0 = MAX(rho(ix,iy), none_zero)
        cs2 = (gamma * pressure(ix,iy) + w1)/ rho0

        !effective speed from viscous pressure
        c_visc2 = p_visc(ix,iy) / rho0

        ! length based on simple DYNA2D estimates
        length = dxb(ix) * dyb(iy) / SQRT(dxb(ix)**2 + dyb(iy)**2)

        ! Find ideal MHD CFL limit for Lagrangian step
        dt1 = length / (SQRT(c_visc2) + SQRT(cs2 + c_visc2))
        dt_local = MIN(dt_local, dt1)

        ! Now find dt for remap step
        ax = 0.5_num * dyb(iy)
        ay = 0.5_num * dxb(ix)
        vxbm = (vx(ixm,iy ) + vx(ixm,iym)) * ax
        vxbp = (vx(ix ,iy ) + vx(ix ,iym)) * ax
        vybm = (vy(ix ,iym) + vy(ixm,iym)) * ay
        vybp = (vy(ix ,iy ) + vy(ixm,iy )) * ay

        dvx = ABS(vxbp - vxbm)
        dvy = ABS(vybp - vybm)
        avxm = ABS(vxbm)
        avxp = ABS(vxbp)
        avym = ABS(vybm)
        avyp = ABS(vybp)

        area = dxb(ix) * dyb(iy)
        dt1 = area / MAX(avxm, avxp, dvx, 1.e-10_num * area)
        dt2 = area / MAX(avym, avyp, dvy, 1.e-10_num * area)

        ! Fix dt for remap step
        dt_local = MIN(dt_local, dt1, dt2)

        ! Note resistive limits assumes uniform resistivity hence cautious
        ! factor 0.2
        dxlocal = 1.0_num / (1.0_num / dxb(ix)**2 + 1.0_num / dyb(iy)**2)

        IF (cowling_resistivity) THEN
          dt1 = 0.2_num * dxlocal &
              / MAX(MAX(eta(ix,iy), eta_perp(ix,iy)), none_zero)
        ELSE
          dt1 = 0.2_num * dxlocal / MAX(eta(ix,iy), none_zero)
        END IF

        ! Adjust to accomodate resistive effects
        dtr_local = MIN(dtr_local, dt1)

        ! Hall MHD CFL limit
        dt1 = 0.75_num * rho(ix,iy) * MIN(dxb(ix), dyb(iy))**2 &
            / MAX(lambda_i(ix,iy) * SQRT(w1), none_zero)

        dth_local = MIN(dth_local, dt1)
      END DO
    END DO

    dt_locals(1) = dt_local
    dt_locals(2) = dtr_local
    dt_locals(3) = dth_local

    CALL MPI_ALLREDUCE(dt_locals, dt_min, 3, mpireal, MPI_MIN, comm, errcode)

    dt  = dt_multiplier * dt_min(1)
    dtr = dt_multiplier * dt_min(2)
    dth = dt_multiplier * dt_min(3)

    time_dump = time_prev + dt_snapshots
    IF (t_end < time_dump) time_dump = t_end

    IF (step < nramp_start) THEN
      ! At start of simulation, slowly ramp up from zero dt
      dt_reason = 's'
      dt = (dt_fudge + (step / REAL(nramp_start,num))**2) * dt
    ELSE IF (nramp > 0) THEN
      ! After output dump, slowly ramp up from old dt
      dt_reason = 'u'
      dt0 = dt_previous + nramp * dt_factor
      nramp = nramp - 1
      dt  = MIN(dt0, dt)
    ELSE IF (nramp < 0) THEN
      ! Slowly ramp down dt until next output dump
      dt_reason = 'd'
      nramp = nramp - 1
      dt0 = dt
      dt  = dt_previous - nramp * dt_factor
      dt  = MIN(dt0, dt)
      IF (nramp == -nrsteps .OR. dt > (time_dump - time)) THEN
        dt = time_dump - time
        nramp = -nramp
      END IF
      IF (dt < t_end * 1.e-10_num) dt = dt0
    ELSE IF (nramp_steps > 0) THEN
      ! Less than nramp steps until next output or end of simulation
      time_rem = time_dump - time
      IF (time_rem > 0.1_num * dt .AND. nramp_steps * dt > time_rem) THEN
        nrsteps = FLOOR(time_rem / dt) + 1
        IF (nramp_steps < nrsteps) nrsteps = nramp_steps
        dt_factor = 2.0_num * (time_rem - nrsteps * dt) &
            / (nrsteps**2 + nrsteps)
        dt_reason = 'd'
        nramp = -1
        IF (nrsteps == 1) nramp = 0
        dt0 = dt
        dt_previous = dt0
        dt  = dt_previous + dt_factor
        dt  = MIN(dt0, dt)
        IF (dt < t_end * 1.e-10_num) dt = dt0
      END IF
    END IF

    IF (is_restart) THEN
      dt_reason = 'r'
      dt = dt_from_restart
      nramp = nramp + 1
    END IF

    time = time + dt

  END SUBROUTINE set_dt




  !****************************************************************************
  ! Calculate the spatial profile of the resistivity at the current timestep
  ! Note that this is a core routine so it works in normalised units
  ! This includes lengths etc.
  !****************************************************************************

  SUBROUTINE eta_calc

    REAL(num) :: jx, jy, jz, jxp, jyp
    INTEGER :: ixp, iyp

    IF (resistive_mhd) THEN
      DO iy = -1, ny + 1
        iyp = iy + 1
        DO ix = -1, nx + 1
          ixp = ix + 1

          ! jx at Ex(i,j)
          jx  = (bz(ix ,iyp) - bz(ix ,iy )) / dyc(iy)

          ! jx at Ex(i+1,j)
          jxp = (bz(ixp,iyp) - bz(ixp,iy )) / dyc(iy)

          ! jy at Ey(i,j)
          jy  = (bz(ix ,iy ) - bz(ixp,iy )) / dxc(ix)

          ! jy at Ey(i,j+1)
          jyp = (bz(ix ,iyp) - bz(ixp,iyp)) / dxc(ix)

          ! jz at Ez(i,j)
          jz  = (by(ixp,iy ) - by(ix ,iy )) / dxc(ix) &
              - (bx(ix ,iyp) - bx(ix ,iy )) / dyc(iy)

          ! Current at V
          jx = (jx + jxp) * 0.5_num
          jy = (jy + jyp) * 0.5_num

          IF (SQRT(jx**2 + jy**2 + jz**2) > j_max) THEN
            eta(ix,iy) = eta_background + eta0
          ELSE
            eta(ix,iy) = eta_background
          END IF
        END DO
      END DO
    ELSE
      eta = 0.0_num
    END IF

  END SUBROUTINE eta_calc



  !****************************************************************************
  ! Calculate the effect of resistivity on the magnetic field and Ohmic heating
  ! Use the subroutine rkstep
  !****************************************************************************

  SUBROUTINE resistive_effects

    REAL(num) :: jx1, jx2, jy1, jy2, jz1
#ifdef FOURTHORDER
    REAL(num) :: dt6, half_dt
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: k1x, k2x, k3x, k4x
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: k1y, k2y, k3y, k4y
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: k1z, k2z, k3z, k4z
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: c1, c2, c3, c4

    ALLOCATE(k1x(0:nx,0:ny), k2x(0:nx,0:ny))
    ALLOCATE(k3x(0:nx,0:ny), k4x(0:nx,0:ny))
    ALLOCATE(k1y(0:nx,0:ny), k2y(0:nx,0:ny))
    ALLOCATE(k3y(0:nx,0:ny), k4y(0:nx,0:ny))
    ALLOCATE(k1z(0:nx,0:ny), k2z(0:nx,0:ny))
    ALLOCATE(k3z(0:nx,0:ny), k4z(0:nx,0:ny))
    ALLOCATE( c1(0:nx,0:ny),  c2(0:nx,0:ny))
    ALLOCATE( c3(0:nx,0:ny),  c4(0:nx,0:ny))
#endif

    bx1 = bx(-1:nx+2,-1:ny+2)
    by1 = by(-1:nx+2,-1:ny+2)
    bz1 = bz(-1:nx+2,-1:ny+2)

    ! Step 1
    CALL rkstep

    ! Default is first order in time
#ifndef FOURTHORDER
    CALL bstep(flux_x, flux_y, flux_z, dt)

    DO iy = 1, ny
      iym = iy - 1
      DO ix = 1, nx
        ixm = ix - 1
        energy(ix,iy) = energy(ix,iy) &
            + (curlb(ix ,iy ) + curlb(ixm,iy )  &
            +  curlb(ix ,iym) + curlb(ixm,iym)) &
            * dt / (4.0_num * rho(ix,iy))
      END DO
    END DO

    CALL energy_bcs

    DO iy = 0, ny
      DO ix = 0, nx
        w1 = dt * dxc(ix) * dyc(iy) * curlb(ix,iy)
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
    ! If complier flag set then use 4th order Runge-Kutta
    half_dt = dt * 0.5_num
    dt6 = dt * sixth

    k1x = flux_x
    k1y = flux_y
    k1z = flux_z
    c1 = curlb

    ! Step 2
    CALL bstep(k1x, k1y, k1z, half_dt)

    CALL rkstep

    k2x = flux_x
    k2y = flux_y
    k2z = flux_z
    c2 = curlb

    ! Step 3
    CALL bstep(k2x, k2y, k2z, half_dt)

    CALL rkstep

    k3x = flux_x
    k3y = flux_y
    k3z = flux_z
    c3 = curlb

    ! Step 4
    CALL bstep(k3x, k3y, k3z, dt)

    CALL rkstep

    k4x = flux_x
    k4y = flux_y
    k4z = flux_z
    c4 = curlb

    ! Full update
    k1x = k1x + 2.0_num * k2x + 2.0_num * k3x + k4x
    k1y = k1y + 2.0_num * k2y + 2.0_num * k3y + k4y
    k1z = k1z + 2.0_num * k2z + 2.0_num * k3z + k4z
    c1 = c1 + 2.0_num * c2 + 2.0_num * c3 + c4

    CALL bstep(k1x, k1y, k1z, dt6)

    DO iy = 1, ny
      iym = iy - 1
      DO ix = 1, nx
        ixm = ix - 1
        energy(ix,iy) = energy(ix,iy) &
            + (c1(ix ,iy ) + c1(ixm,iy ) + c1(ix ,iym) + c1(ixm,iym)) &
            * dt6 / (4.0_num * rho(ix,iy))
      END DO
    END DO

    CALL energy_bcs

    DO iy = 0, ny
      DO ix = 0, nx
        w1 = dt6 * dxc(ix) * dyc(iy) * c1(ix,iy)
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

    DO iy = 0, ny
      iyp = iy + 1
      DO ix = 0, nx
        ixp = ix + 1

        jx1 = (bz(ix ,iyp) - bz(ix ,iy )) / dyc(iy)
        jx2 = (bz(ixp,iyp) - bz(ixp,iy )) / dyc(iy)
        jy1 = (bz(ix ,iy ) - bz(ixp,iy )) / dxc(ix)
        jy2 = (bz(ix ,iyp) - bz(ixp,iyp)) / dxc(ix)
        jz1 = (by(ixp,iy ) - by(ix ,iy )) / dxc(ix) &
            - (bx(ix ,iyp) - bx(ix ,iy )) / dyc(iy)

        jx_r(ix,iy) = (jx1 + jx2) * 0.5_num
        jy_r(ix,iy) = (jy1 + jy2) * 0.5_num
        jz_r(ix,iy) = jz1
      END DO
    END DO

    ! Once more to get j_perp and j_par correct
    CALL rkstep

#ifdef FOURTHORDER
    DEALLOCATE(k1x, k2x, k3x, k4x, k1y, k2y, k3y, k4y, k1z, k2z, k3z, k4z)
    DEALLOCATE(c1, c2, c3, c4)
#endif

    DO iy = -1, ny + 2
      iym = iy - 1
      DO ix = -1, nx + 2
        ixm = ix - 1
        bx1(ix,iy) = (bx(ix,iy) + bx(ixm,iy )) * 0.5_num
        by1(ix,iy) = (by(ix,iy) + by(ix ,iym)) * 0.5_num
        bz1(ix,iy) = bz(ix,iy)

        pressure(ix,iy) = (gamma - 1.0_num) * rho(ix,iy) &
            * (energy(ix,iy) - (1.0_num - xi_n(ix,iy)) * ionise_pot)
      END DO
    END DO

  END SUBROUTINE resistive_effects



  !****************************************************************************
  ! Calculates 'k' values from b[xyz]1 values
  !****************************************************************************

  SUBROUTINE rkstep

    REAL(num) :: jx, jy, jz
    REAL(num) :: jx1, jy1, jx2, jy2
    REAL(num) :: bxv, byv, bzv
    REAL(num) :: magn_b
    REAL(num) :: j_par_x, j_par_y, j_par_z
    REAL(num) :: j_perp_x, j_perp_y, j_perp_z
    REAL(num) :: magn_j_perp, magn_j_par

    IF (.NOT. cowling_resistivity) THEN
      ! Use simple flux calculation
      DO iy = 0, ny
        iyp = iy + 1
        DO ix = 0, nx
          ixp = ix + 1

          jx1 = (bz(ix ,iyp) - bz(ix ,iy )) / dyc(iy)
          jx2 = (bz(ixp,iyp) - bz(ixp,iy )) / dyc(iy)
          jy1 = (bz(ix ,iy ) - bz(ixp,iy )) / dxc(ix)
          jy2 = (bz(ix ,iyp) - bz(ixp,iyp)) / dxc(ix)
          jz  = (by(ixp,iy ) - by(ix ,iy )) / dxc(ix) &
              - (bx(ix ,iyp) - bx(ix ,iy )) / dyc(iy)

          jx = (jx1 + jx2) * 0.5_num
          jy = (jy1 + jy2) * 0.5_num

          flux_x(ix,iy) = -jx * eta(ix,iy) * dxc(ix) * 0.5_num
          flux_y(ix,iy) = -jy * eta(ix,iy) * dyc(iy) * 0.5_num
          flux_z(ix,iy) = -jz * eta(ix,iy)
          ! This isn't really curlb. It's actually heat flux
          curlb(ix,iy) = eta(ix,iy) * (jx**2 + jy**2 + jz**2)
        END DO
      END DO
    ELSE
      ! Use partially ionised flux calculation
      DO iy = 0, ny
        iyp = iy + 1
        DO ix = 0, nx
          ixp = ix + 1

          jx1 = (bz(ix ,iyp) - bz(ix ,iy )) / dyc(iy)
          jx2 = (bz(ixp,iyp) - bz(ixp,iy )) / dyc(iy)
          jy1 = (bz(ix ,iy ) - bz(ixp,iy )) / dxc(ix)
          jy2 = (bz(ix ,iyp) - bz(ixp,iyp)) / dxc(ix)
          jz  = (by(ixp,iy ) - by(ix ,iy )) / dxc(ix) &
              - (bx(ix ,iyp) - bx(ix ,iy )) / dyc(iy)

          jx = (jx1 + jx2) * 0.5_num
          jy = (jy1 + jy2) * 0.5_num

          ! B at vertices
          bxv = (bx(ix,iy ) + bx(ix ,iyp)) * 0.5_num
          byv = (by(ix,iy ) + by(ixp,iy )) * 0.5_num
          bzv = (bz(ix,iy ) + bz(ixp,iy ) &
              +  bz(ix,iyp) + bz(ixp,iyp)) * 0.25_num

          magn_b = bxv**2 + byv**2 + bzv**2

          ! Calculate parallel and perpendicular currents
          IF (magn_b > none_zero) THEN
            j_par_x = (jx * bxv + jy * byv + jz * bzv) * bxv / magn_b
            j_par_y = (jx * bxv + jy * byv + jz * bzv) * byv / magn_b
            j_par_z = (jx * bxv + jy * byv + jz * bzv) * bzv / magn_b
          ELSE
            ! If b = 0 then there is no parallel current
            j_par_x = 0.0_num
            j_par_y = 0.0_num
            j_par_z = 0.0_num
          END IF

          ! Calculate perpendicular current
          j_perp_x = jx - j_par_x
          j_perp_y = jy - j_par_y
          j_perp_z = jz - j_par_z

          magn_j_par  = SQRT(j_par_x**2 + j_par_y**2 + j_par_z**2)
          magn_j_perp = SQRT(j_perp_x**2 + j_perp_y**2 + j_perp_z**2)

          parallel_current(ix,iy) = magn_j_par
          perp_current(ix,iy) = magn_j_perp

          ! This isn't really curlb. It's actually heat flux
          curlb(ix,iy) = eta(ix,iy) * magn_j_par**2 &
              + (eta_perp(ix,iy) + eta(ix,iy)) * magn_j_perp**2

          flux_x(ix,iy) = -((j_par_x + j_perp_x) * eta(ix,iy) &
              + j_perp_x * eta_perp(ix,iy)) * dxc(ix) * 0.5_num
          flux_y(ix,iy) = -((j_par_y + j_perp_y) * eta(ix,iy) &
              + j_perp_y * eta_perp(ix,iy)) * dyc(iy) * 0.5_num
          flux_z(ix,iy) = -((j_par_z + j_perp_z) * eta(ix,iy) &
              + j_perp_z * eta_perp(ix,iy))
        END DO
      END DO
    END IF

  END SUBROUTINE rkstep



  !****************************************************************************
  ! Calculate the effect of the Hall term on the magnetic field
  ! Uses the subroutine rkstep1
  !****************************************************************************

  SUBROUTINE hall_effects

    REAL(num) :: dt6, half_dt
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: k1x, k2x, k3x, k4x
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: k1y, k2y, k3y, k4y
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: k1z, k2z, k3z, k4z

    ALLOCATE(k1x(0:nx,0:ny), k2x(0:nx,0:ny), k3x(0:nx,0:ny), k4x(0:nx,0:ny))
    ALLOCATE(k1y(0:nx,0:ny), k2y(0:nx,0:ny), k3y(0:nx,0:ny), k4y(0:nx,0:ny))
    ALLOCATE(k1z(0:nx,0:ny), k2z(0:nx,0:ny), k3z(0:nx,0:ny), k4z(0:nx,0:ny))

    half_dt = dt * 0.5_num
    dt6 = dt * sixth

    bx1 = bx(0:nx+1,0:ny+1)
    by1 = by(0:nx+1,0:ny+1)
    bz1 = bz(0:nx+1,0:ny+1)

    ! Step 1
    CALL rkstep1(half_dt)

    k1x = flux_x
    k1y = flux_y
    k1z = flux_z

    ! Step 2
    CALL bstep(k1x, k1y, k1z, half_dt)

    CALL rkstep1(half_dt)

    k2x = flux_x
    k2y = flux_y
    k2z = flux_z

    ! Step 3
    CALL bstep(k2x, k2y, k2z, half_dt)

    CALL rkstep1(half_dt)

    k3x = flux_x
    k3y = flux_y
    k3z = flux_z

    ! Step 4
    CALL bstep(k3x, k3y, k3z, dt)

    CALL rkstep1(dt)

    k4x = flux_x
    k4y = flux_y
    k4z = flux_z

    ! Full update
    k1x = k1x + 2.0_num * k2x + 2.0_num * k3x + k4x
    k1y = k1y + 2.0_num * k2y + 2.0_num * k3y + k4y
    k1z = k1z + 2.0_num * k2z + 2.0_num * k3z + k4z

    CALL bstep(k1x, k1y, k1z, dt6)

    DEALLOCATE(k1x, k2x, k3x, k4x, k1y, k2y, k3y, k4y, k1z, k2z, k3z, k4z)

  END SUBROUTINE hall_effects



  SUBROUTINE rkstep1(dt)

    REAL(num), INTENT(IN) :: dt
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: jx, jy, jz
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: hxflux, hyflux, hzflux
    REAL(num) :: jx1, jy1, jz1, jx2, jy2
    REAL(num) :: bxv, byv, bzv, rho_v
    REAL(num) :: f1, f2, area
    REAL(num) :: j_advect, jad_p, jad_m
    REAL(num) :: fm, fi, fp, fp2, fu, dfm, dfi, dfp, dfu, sign_j
    REAL(num) :: dxci, dxcu, dxbu, phi, ss, Da, Di
    INTEGER :: ixp2, iyp2

    ALLOCATE(jx(-1:nx+1,-1:ny+1))
    ALLOCATE(jy(-1:nx+1,-1:ny+1))
    ALLOCATE(jz(-1:nx+1,-1:ny+1))
    ALLOCATE(hxflux(-1:nx+1,-1:ny+1))
    ALLOCATE(hyflux(-1:nx+1,-1:ny+1))
    ALLOCATE(hzflux(-1:nx+1,-1:ny+1))

    DO iy = -1, ny + 1
      iyp = iy + 1
      DO ix = -1, nx + 1
        ixp = ix + 1

        jx1 = (bz(ix ,iyp) - bz(ix ,iy )) / dyc(iy)
        jx2 = (bz(ixp,iyp) - bz(ixp,iy )) / dyc(iy)
        jy1 = (bz(ix ,iy ) - bz(ixp,iy )) / dxc(ix)
        jy2 = (bz(ix ,iyp) - bz(ixp,iyp)) / dxc(ix)
        jz1 = (by(ixp,iy ) - by(ix ,iy )) / dxc(ix) &
            - (bx(ix ,iyp) - bx(ix ,iy )) / dyc(iy)

        jx(ix,iy) = (jx1 + jx2) * 0.5_num
        jy(ix,iy) = (jy1 + jy2) * 0.5_num
        jz(ix,iy) = jz1
      END DO
    END DO

    hxflux = 0.0_num
    hyflux = 0.0_num
    hzflux = 0.0_num

    DO iy = 0, ny
      iym  = iy - 1
      iyp  = iy + 1
      iyp2 = iy + 2
      DO ix = 0, nx
        ixm  = ix - 1
        ixp  = ix + 1
        ixp2 = ix + 2

        bxv = (bx(ix,iy) * dyb(iy) + bx(ix,iyp) * dyb(iyp)) &
            / (dyb(iy) + dyb(iyp))
        byv = (by(ix,iy) * dxb(ix) + by(ixp,iy) * dxb(ixp)) &
            / (dxb(ix) + dxb(ixp))
        bzv = (bz(ix,iy ) * cv(ix,iy ) + bz(ixp,iy ) * cv(ixp,iy ) &
            +  bz(ix,iyp) * cv(ix,iyp) + bz(ixp,iyp) * cv(ixp,iyp))

        rho_v = (rho(ix,iy ) * cv(ix,iy ) + rho(ixp,iy ) * cv(ixp,iy ) &
            +    rho(ix,iyp) * cv(ix,iyp) + rho(ixp,iyp) * cv(ixp,iyp))

        area = (cv(ix,iy) + cv(ixp,iy) + cv(ix,iyp) + cv(ixp,iyp))
        bzv = bzv / area
        rho_v = rho_v / area

        ! Evaluate the flux due to the Hall term but upwind in direction of
        ! negative current density i.e. upwind in the direction of the
        ! effective advection velocity
        fm  = (bz(ix,iym ) + bz(ixp,iym )) * 0.5_num
        fi  = (bz(ix,iy  ) + bz(ixp,iy  )) * 0.5_num
        fp  = (bz(ix,iyp ) + bz(ixp,iyp )) * 0.5_num
        fp2 = (bz(ix,iyp2) + bz(ixp,iyp2)) * 0.5_num

        dfm = fi - fm
        dfi = fp - fi
        dfp = fp2 - fp

        j_advect = jy(ix,iy)
        sign_j = -SIGN(1.0_num, j_advect)
        jad_p = (sign_j + 1.0_num) * 0.5_num
        jad_m = 1.0_num - jad_p

        fu = fi * jad_p + fp * jad_m
        dfu = dfm * jad_p + dfp * jad_m
        dxci = dyc(iy)
        dxcu = dyc(iym) * jad_p + dyc(iyp) * jad_m
        dxbu = dyb(iy ) * jad_p + dyb(iyp) * jad_m

        phi = lambda_i(ix,iy) * ABS(j_advect) * dt / rho_v / dxbu

        Da =  (2.0_num - phi) * ABS(dfi) / dxci &
            + (1.0_num + phi) * ABS(dfu) / dxcu
        Da = Da * sixth

        ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

        Di = sign_j * ss * MIN(ABS(Da) * dxbu, ABS(dfi), ABS(dfu))

        f1 = (fu + Di * (1.0_num - phi)) * j_advect
        f2 = jz(ix,iy) * byv
        hxflux(ix,iy) = lambda_i(ix,iy) * (f1 - f2) / rho_v

        fm  = (bz(ixm ,iy) + bz(ixm ,iyp)) * 0.5_num
        fi  = (bz(ix  ,iy) + bz(ix  ,iyp)) * 0.5_num
        fp  = (bz(ixp ,iy) + bz(ixp ,iyp)) * 0.5_num
        fp2 = (bz(ixp2,iy) + bz(ixp2,iyp)) * 0.5_num

        dfm = fi - fm
        dfi = fp - fi
        dfp = fp2 - fp

        j_advect = jx(ix,iy)
        sign_j = -SIGN(1.0_num, j_advect)
        jad_p = (sign_j + 1.0_num) * 0.5_num
        jad_m = 1.0_num - jad_p

        fu = fi * jad_p + fp * jad_m
        dfu = dfm * jad_p + dfp * jad_m
        dxci = dxc(ix)
        dxcu = dxc(ixm) * jad_p + dxc(ixp) * jad_m
        dxbu = dxb(ix ) * jad_p + dxb(ixp) * jad_m

        phi = lambda_i(ix,iy) * ABS(j_advect) * dt / rho_v / dxbu

        Da =  (2.0_num - phi) * ABS(dfi) / dxci &
            + (1.0_num + phi) * ABS(dfu) / dxcu
        Da = Da * sixth

        ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

        Di = sign_j * ss * MIN(ABS(Da) * dxbu, ABS(dfi), ABS(dfu))

        f1 = jz(ix,iy) * bxv
        f2 = (fu + Di * (1.0_num - phi)) * j_advect

        hyflux(ix,iy) = lambda_i(ix,iy) * (f1 - f2) / rho_v

        fm  = by(ixm ,iy)
        fi  = by(ix  ,iy)
        fp  = by(ixp ,iy)
        fp2 = by(ixp2,iy)

        dfm = fi - fm
        dfi = fp - fi
        dfp = fp2 - fp

        j_advect = jx(ix,iy)
        sign_j = -SIGN(1.0_num, j_advect)
        jad_p = (sign_j + 1.0_num) * 0.5_num
        jad_m = 1.0_num - jad_p

        fu = fi * jad_p + fp * jad_m
        dfu = dfm * jad_p + dfp * jad_m
        dxci = dxc(ix)
        dxcu = dxc(ixm) * jad_p + dxc(ixp) * jad_m
        dxbu = dxb(ix ) * jad_p + dxb(ixp) * jad_m

        phi = lambda_i(ix,iy) * ABS(j_advect) * dt / rho_v / dxbu

        Da =  (2.0_num - phi) * ABS(dfi) / dxci &
            + (1.0_num + phi) * ABS(dfu) / dxcu
        Da = Da * sixth

        ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

        Di = sign_j * ss * MIN(ABS(Da) * dxbu, ABS(dfi), ABS(dfu))

        f1 = (fu + Di * (1.0_num - phi)) * j_advect

        fm  = bx(ix,iym )
        fi  = bx(ix,iy  )
        fp  = bx(ix,iyp )
        fp2 = bx(ix,iyp2)

        dfm = fi - fm
        dfi = fp - fi
        dfp = fp2 - fp

        j_advect = jy(ix,iy)
        sign_j = -SIGN(1.0_num, j_advect)
        jad_p = (sign_j + 1.0_num) * 0.5_num
        jad_m = 1.0_num - jad_p

        fu = fi * jad_p + fp * jad_m
        dfu = dfm * jad_p + dfp * jad_m
        dxci = dyc(iy)
        dxcu = dyc(iym) * jad_p + dyc(iy) * jad_m
        dxbu = dyb(iy ) * jad_p + dyb(iyp) * jad_m

        phi = lambda_i(ix,iy) * ABS(j_advect) * dt / rho_v / dxbu

        Da =  (2.0_num - phi) * ABS(dfi) / dxci &
            + (1.0_num + phi) * ABS(dfu) / dxcu
        Da = Da * sixth

        ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

        Di = sign_j * ss * MIN(ABS(Da) * dxbu, ABS(dfi), ABS(dfu))

        f2 = (fu + Di * (1.0_num - phi)) * j_advect

        hzflux(ix,iy) = lambda_i(ix,iy) * (f1 - f2) / rho_v
      END DO
    END DO

    DO iy = 0, ny
      DO ix = 0, nx
        flux_x(ix,iy) = -hxflux(ix,iy) * dxc(ix) * 0.5_num
        flux_y(ix,iy) = -hyflux(ix,iy) * dyc(iy) * 0.5_num
        flux_z(ix,iy) = -hzflux(ix,iy)
      END DO
    END DO

    DEALLOCATE (jx, jy, jz, hxflux, hyflux, hzflux)

  END SUBROUTINE rkstep1



  SUBROUTINE bstep(kx, ky, kz, dt)

    REAL(num), DIMENSION(0:,0:), INTENT(IN) :: kx, ky, kz
    REAL(num), INTENT(IN) :: dt

    DO iy = 1, ny
      iym = iy - 1
      DO ix = 0, nx
        bx(ix,iy) = bx1(ix,iy) &
            + (kz(ix,iy ) - kz(ix,iym)) * dt / dyb(iy)
      END DO
    END DO

    DO iy = 0, ny
      DO ix = 1, nx
        ixm = ix - 1
        by(ix,iy) = by1(ix,iy) &
            - (kz(ix ,iy) - kz(ixm,iy)) * dt / dxb(ix)
      END DO
    END DO

    DO iy = 1, ny
      iym = iy - 1
      DO ix = 1, nx
        ixm = ix - 1
        bz(ix,iy) = bz1(ix,iy) &
            + (ky(ix ,iy ) - ky(ixm,iy ) &
            +  ky(ix ,iym) - ky(ixm,iym) &
            -  kx(ix ,iy ) + kx(ix ,iym) &
            -  kx(ixm,iy ) + kx(ixm,iym)) * dt / cv(ix,iy)
      END DO
    END DO

    CALL bfield_bcs

  END SUBROUTINE bstep

END MODULE lagran
