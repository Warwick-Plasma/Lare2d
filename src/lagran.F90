MODULE lagran

  USE shared_data; USE boundary
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: lagrangian_step, eta_calc

  ! only used inside lagran.f90
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: qxy, qxz, qyz
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: qxx, qyy, visc_heat, pressure
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: flux_x, flux_y, flux_z, curlb

CONTAINS


  SUBROUTINE lagrangian_step

    INTEGER :: substeps, subcycle
    REAL(num) :: actual_dt, dt_sub

    ALLOCATE (bx1(0:nx+1,0:ny+1),by1(0:nx+1,0:ny+1),bz1(0:nx+1,0:ny+1), &
         qxy(0:nx+1,0:ny+1),qxz(0:nx+1,0:ny+1),qyz(0:nx+1,0:ny+1), &
         qxx(0:nx+1,0:ny+1),qyy(0:nx+1,0:ny+1), &
         visc_heat(0:nx+1,0:ny+1),pressure(0:nx+1,0:ny+1), &
         flux_x(0:nx,0:ny),flux_y(0:nx,0:ny),flux_z(0:nx,0:ny),curlb(0:nx,0:ny))

    IF (resistiveMHD) THEN
       ! if subcycling isn't wanted set dt = dtr in set_dt, don't just
       ! set substeps to 1.
       dt_sub = dtr
       substeps = INT(dt/dt_sub) + 1
       IF (substeps > peak_substeps) peak_substeps = substeps
       actual_dt = dt
       dt = dt / REAL(substeps,num)
       DO subcycle = 1, substeps
          CALL eta_calc
          IF (resistiveMHD) CALL resistive_effects  
       ENDDO
       dt = actual_dt
    ENDIF

    DO iy = 0, ny+1
       DO ix = 0, nx+1
          ixm = ix - 1
          iym = iy - 1
          bx1(ix,iy) = (bx(ix,iy) + bx(ixm,iy)) / 2.0_num
          by1(ix,iy) = (by(ix,iy) + by(ix,iym)) / 2.0_num
       END DO
    END DO
    bz1 = bz(0:nx+1,0:ny+1)      ! bz1 = bz at C

    CALL predictor_corrector_step

    DEALLOCATE (bx1, by1, bz1, qxy, qxz, qyz, qxx, qyy, visc_heat, pressure, &
         flux_x, flux_y, flux_z, curlb)

    CALL energy_bcs
    CALL density_bcs
    CALL velocity_bcs

  END SUBROUTINE lagrangian_step



  SUBROUTINE predictor_corrector_step

    REAL(num) :: p, pxp, pyp, pxpyp, p1, p2, p3, p4
    REAL(num) :: e1xp, e1yp, e1xpyp

    REAL(num) :: e1, rho_v
    REAL(num) :: fx, fy, fz
    REAL(num) :: vxb, vxbm, vyb, vybm
    REAL(num) :: dv, dvxp, dvyp, dvxpyp
    REAL(num) :: bxv, byv, bzv, jx, jy, jz
    REAL(num) :: cvx, cvxp, cvy, cvyp

    REAL(num) :: mx,mn,momentum

    dt2 = dt / 2.0_num
    CALL viscosity_and_b_update

    bx1 = bx1*cv1(0:nx+1,0:ny+1)
    by1 = by1*cv1(0:nx+1,0:ny+1)
    bz1 = bz1*cv1(0:nx+1,0:ny+1)
    momentum=0.0_num

    DO iy = 0, ny+1
          !DEC$ IVDEP
          !DEC$ VECTOR ALWAYS
       DO ix = 0, nx+1
          dv = cv1(ix,iy) / cv(ix,iy) - 1.0_num
          e1 = energy(ix,iy) - pressure(ix,iy) * dv/rho(ix,iy)   !predictor energy
#ifdef Q_MONO
          IF (visc3 > 1.e-6_num) THEN
#endif
             e1 = e1  + visc_heat(ix,iy)*dt2/rho(ix,iy)   
#ifdef Q_MONO
          END IF
#endif
#ifdef Q_MONO
          e1 = e1 - p_visc(ix,iy) * dv/rho(ix,iy) 
#endif

          ! now define the predictor step pressures
          pressure(ix,iy) = e1 * (gamma - 1.0_num) * rho(ix,iy) * cv(ix,iy) / cv1(ix,iy)
          ! add shock viscosity
#ifdef Q_MONO
          pressure(ix,iy) = pressure(ix,iy) + p_visc(ix,iy) 
#endif
       END DO
    END DO

    DO iy = 0, ny
          !DEC$ IVDEP
          !DEC$ VECTOR ALWAYS
       DO ix = 0, nx
          ixp = ix + 1
          iyp = iy + 1
          w1 = (pressure(ix,iy) + pressure(ix,iyp)) / 2.0_num     ! P total at Ey(i,j)
          w2 = (pressure(ixp,iy) + pressure(ixp,iyp)) / 2.0_num        ! P total at Ey(i+1,j)
          fx = - (w2 - w1) / dxc(ix)
          w1 = (pressure(ix,iy) + pressure(ixp,iy)) / 2.0_num     ! P total at Ex(i,j)
          w2 = (pressure(ix,iyp) + pressure(ixp,iyp)) / 2.0_num        ! P total at Ex(i,j+1)
          fy = - (w2 - w1) / dyc(iy)
          fz = 0.0_num
#ifdef Q_MONO
          IF (visc3 > 1.e-6_num) THEN
#endif
             ! add parallel component of viscosity
             w1 = (qxx(ix,iy) + qxx(ix,iyp)) / 2.0_num     
             w2 = (qxx(ixp,iy) + qxx(ixp,iyp)) / 2.0_num        
             fx = fx + (w2 - w1) / dxc(ix)
             w1 = (qyy(ix,iy) + qyy(ixp,iy)) / 2.0_num     
             w2 = (qyy(ix,iyp) + qyy(ixp,iyp)) / 2.0_num        
             fy = fy + (w2 - w1) / dyc(iy)

             ! add shear forces
             w1 = (qxy(ix,iy) + qxy(ixp,iy)) / 2.0_num
             w2 = (qxy(ix,iyp) + qxy(ixp,iyp)) / 2.0_num
             fx = fx + (w2 - w1) / dyc(iy)
             w1 = (qxy(ix,iy) + qxy(ix,iyp)) / 2.0_num
             w2 = (qxy(ixp,iy) + qxy(ixp,iyp)) / 2.0_num
             fy = fy + (w2 - w1) / dxc(ix)
             w1 = (qxz(ix,iy)+qxz(ix,iyp)) / 2.0_num
             w2 = (qxz(ixp,iy)+qxz(ixp,iyp)) / 2.0_num
             fz = (w2 - w1) / dxc(ix)
             w1 = (qyz(ix,iy)+qyz(ixp,iy)) / 2.0_num
             w2 = (qyz(ix,iyp)+qyz(ixp,iyp)) / 2.0_num
             fz = fz + (w2 - w1) / dyc(iy)
#ifdef Q_MONO
          END IF
#endif

          cvx = cv1(ix,iy)+cv1(ix,iyp)
          cvxp = cv1(ixp,iy)+cv1(ixp,iyp)
          cvy = cv1(ix,iy)+cv1(ixp,iy)
          cvyp = cv1(ix,iyp)+cv1(ixp,iyp)

          w1 = (bz1(ix,iy)+bz1(ixp,iy)) / cvy
          w2 = (bz1(ix,iyp)+bz1(ixp,iyp)) / cvyp
          jx = (w2 - w1) / dyc(iy)

          w1 = (bz1(ix,iy)+bz1(ix,iyp)) / cvx
          w2 = (bz1(ixp,iy)+bz1(ixp,iyp)) / cvxp
          jy = - (w2 - w1) / dxc(ix)

          w1 = (by1(ix,iy)+by1(ix,iyp)) / cvx
          w2 = (by1(ixp,iy)+by1(ixp,iyp)) / cvxp
          jz = (w2 - w1) / dxc(ix)
          w1 = (bx1(ix,iy)+bx1(ixp,iy)) / cvy
          w2 = (bx1(ix,iyp)+bx1(ixp,iyp)) / cvyp
          jz = jz - (w2 - w1) / dyc(iy)

          bxv = (bx1(ix,iy)+bx1(ixp,iy)+bx1(ix,iyp)+bx1(ixp,iyp)) &
               / (cvx + cvxp)
          byv = (by1(ix,iy)+by1(ixp,iy)+by1(ix,iyp)+by1(ixp,iyp)) &
               / (cvx + cvxp)
          bzv = (bz1(ix,iy)+bz1(ixp,iy)+bz1(ix,iyp)+bz1(ixp,iyp)) &
               / (cvx + cvxp)

          fx = fx + (jy*bzv - jz*byv) 
          fy = fy + (jz*bxv - jx*bzv)
          fz = fz + (jx*byv - jy*bxv) 

          rho_v = (rho(ix,iy)*cv(ix,iy) + rho(ixp,iy)*cv(ixp,iy) &
               + rho(ix,iyp)*cv(ix,iyp) +rho(ixp,iyp)*cv(ixp,iyp))
          rho_v = rho_v / (cv(ix,iy) + cv(ixp,iy) + cv(ix,iyp) + cv(ixp,iyp))

          fy = fy - rho_v * grav(iy)

          ! find half step velocity needed for remap
          vx1(ix,iy) = vx(ix,iy) + dt2 * fx / rho_v
          vy1(ix,iy) = vy(ix,iy) + dt2 * fy / rho_v
          vz1(ix,iy) = vz(ix,iy) + dt2 * fz / rho_v

          ! velocity at the end of the Lagrangian step
          vx(ix,iy) = vx(ix,iy) + dt * fx / rho_v
          vy(ix,iy) = vy(ix,iy) + dt * fy / rho_v
          vz(ix,iy) = vz(ix,iy) + dt * fz / rho_v

       END DO
    END DO

    CALL store_boundary_dv
    CALL remap_v_bcs
#ifdef Q_MONO
    IF (visc3 > 1.e-6_num) THEN
#endif
       CALL visc_heating
#ifdef Q_MONO
    ENDIF
#endif

    ! finally correct density and energy to final values
    DO iy = 1, ny
          !DEC$ IVDEP
          !DEC$ VECTOR ALWAYS
       DO ix = 1, nx
          ixm = ix - 1
          iym = iy - 1

          vxb = (vx1(ix,iy) + vx1(ix,iym)) / 2.0_num     !vx1 at Ex(i,j)
          vxbm = (vx1(ixm,iy) + vx1(ixm,iym)) / 2.0_num  !vx1 at Ex(i-1,j)
          vyb = (vy1(ix,iy) + vy1(ixm,iy)) / 2.0_num     !vy1 at Ey(i,j)
          vybm = (vy1(ix,iym) + vy1(ixm,iym)) / 2.0_num  !vy1 at Ey(i-1,j)
          dv = ((vxb - vxbm)/dxb(ix) + (vyb - vybm)/dyb(iy)) * dt
          cv1(ix,iy) = cv(ix,iy) * (1.0_num + dv)

          !(dv > 0) ? pvisc:0
          p_visc(ix,iy) = p_visc(ix,iy) * SIGN(1.0_num,MAX(dv,0.0_num)) 
          !energy at end of Lagrangian step
#ifdef Q_MONO
          pressure(ix,iy) = pressure(ix,iy) - p_visc(ix,iy) 
#endif
          energy(ix,iy) = energy(ix,iy)  &
               - pressure(ix,iy) * dv / rho(ix,iy) 
#ifdef Q_MONO
          IF (visc3 > 1.e-6_num) THEN
#endif
             energy(ix,iy) = energy(ix,iy)  &
                  + dt * visc_heat(ix,iy) / rho(ix,iy)
#ifdef Q_MONO
          END IF
          energy(ix,iy) = energy(ix,iy)  &
               - p_visc(ix,iy) * dv / rho(ix,iy) 
#endif

          rho(ix,iy) = rho(ix,iy) / (1.0_num + dv) 
#ifdef Q_MONO
          IF (visc3 > 1.e-6_num) THEN
#endif
             total_visc_heating = total_visc_heating  &
                  + dt * visc_heat(ix,iy) * cv(ix,iy)
#ifdef Q_MONO
          END IF
          total_visc_heating = total_visc_heating  &
               - p_visc(ix,iy) * dv * cv(ix,iy) 
#endif

       END DO
    END DO

  END SUBROUTINE predictor_corrector_step



  SUBROUTINE viscosity_and_b_update

    REAL(num) :: vxb, vxbm, vyb, vybm, vzb, vzbm, dvxdy, dvydx
    REAL(num) :: p, pxp, pxm, pyp, pym, fx, fy, dv
    REAL(num) :: dvxdx, dvydy, dvxy, dvzdx, dvzdy, s, L, cf, L2
    REAL(num) :: sxx, syy, sxy, sxz, syz
    REAL(num) :: case1,case2,flip

    p_visc = 0.0_num
    pressure = energy(0:nx+1,0:ny+1)*(gamma-1.0_num)*rho(0:nx+1,0:ny+1)

    DO iy = 0, ny+1
          !DEC$ IVDEP
          !DEC$ VECTOR ALWAYS
       DO ix = 0, nx+1
          ixp = ix + 1
          iyp = iy + 1
          ixm = ix - 1
          iym = iy - 1
          vxb = (vx(ix,iy) + vx(ix,iym)) / 2.0_num     !vx at Ex(i,j)
          vxbm = (vx(ixm,iy) + vx(ixm,iym)) / 2.0_num  !vx at Ex(i-1,j)
          vyb = (vy(ix,iy) + vy(ixm,iy)) / 2.0_num     !vy at Ey(i,j)
          vybm = (vy(ix,iym) + vy(ixm,iym)) / 2.0_num  !vy at Ey(i-1,j)
          dv = ((vxb - vxbm)/dxb(ix) + (vyb - vybm)/dyb(iy)) * dt2
          cv1(ix,iy) = cv(ix,iy) * (1.0_num + dv)

          dvxdx = (vxb - vxbm) / dxb(ix)
          dvydy = (vyb - vybm) / dyb(iy)

          vxb = (vx(ix,iy) + vx(ixm,iy)) / 2.0_num     !vx at Ey(i,j)
          vxbm = (vx(ix,iym) + vx(ixm,iym)) / 2.0_num  !vx at Ey(i,j-1)
          vyb = (vy(ix,iy) + vy(ix,iym)) / 2.0_num     !vy at Ex(i,j)
          vybm = (vy(ixm,iy) + vy(ixm,iym)) / 2.0_num  !vy at Ex(i-1,j)

          dvxdy = (vxb - vxbm)/dyb(iy)
          dvydx = (vyb - vybm)/dxb(ix)
          dvxy = dvxdy + dvydx

          sxy = dvxy / 2.0_num
          sxx = 2.0_num * dvxdx / 3.0_num - dvydy / 3.0_num
          syy = 2.0_num * dvydy / 3.0_num - dvxdx / 3.0_num

          vzb = (vz(ix,iy) + vz(ix,iym)) / 2.0_num
          vzbm = (vz(ixm,iy) + vz(ixm,iym)) / 2.0_num
          dvzdx = (vzb - vzbm) / dxb(ix) 
          sxz = dvzdx / 2.0_num

          vzb = (vz(ix,iy) + vz(ixm,iy)) / 2.0_num
          vzbm = (vz(ix,iym) + vz(ixm,iym)) / 2.0_num
          dvzdy = (vzb - vzbm) / dyb(iy) 
          syz = dvzdy / 2.0_num

          p = energy(ix,iy) * (gamma - 1.0_num) * rho(ix,iy)
          pxp = energy(ixp,iy) * (gamma - 1.0_num) * rho(ixp,iy)
          pxm = energy(ixm,iy) * (gamma - 1.0_num) * rho(ixm,iy)
          pyp = energy(ix,iyp) * (gamma - 1.0_num) * rho(ix,iyp)
          pym = energy(ix,iym) * (gamma - 1.0_num) * rho(ix,iym)

          fx = - (pxp - pxm) / dxb(ix) 
          fy = - (pyp - pym) / dyb(iy)  

          w1 = fx**2 + fy**2
          s = (dvxdx*fx**2 + dvydy*fy**2 + dvxy*fx*fy)   &
               / MAX(w1, none_zero)

          !(dxb(ix)*ABS(fy) > dyb(iy)*ABS(fx)) ? 1:-1
          flip=SIGN(1.0_num,dxb(ix)*ABS(fy)-dyb(iy)*ABS(fx))  
          case1 = -MIN(flip,0.0_num) !(flip < 0) ? 1:0
          case2 = MAX(flip,0.0_num)  ! (flip > 0) ? 1:0

          w2 = dxb(ix)**2*w1 / MAX(fx**2, none_zero) * case1&
               + dyb(iy)**2*w1 / MAX(fy**2, none_zero) *case2&
               + dxb(ix)**2+dyb(iy)**2 * (1.0_num-case1) * (1.0_num-case2)

          flip = -MIN(SIGN(1.0_num,w1-1.e-6_num),0.0_num) !(w1 < 1.e-6) ? 1:0

          w2 = w2 * (1.0_num-flip) + flip*MIN(dxb(ix),dyb(iy))**2
          L = SQRT(w2)

          L2 = L
          IF (s > 0.0_num .OR. dv > 0.0_num) L = 0.0_num

          w1 = bx1(ix,iy)**2 + by1(ix,iy)**2 + bz1(ix,iy)**2
          cf = SQRT((w1+gamma*p)/ rho(ix,iy))

          p_visc(ix,iy) = visc1*ABS(s)*L*cf*rho(ix,iy) &
               + visc2*(s*L)**2*rho(ix,iy)

          qxy(ix,iy) = 0.0_num
          qxz(ix,iy) = 0.0_num
          qyz(ix,iy) = 0.0_num
          qxx(ix,iy) = 0.0_num
          qyy(ix,iy) = 0.0_num
#ifndef Q_MONO      
          qxy(ix,iy) = sxy * (L2 * rho(ix,iy)  &
               * (visc1 * cf + L2 * visc2 * ABS(s)))
          qxz(ix,iy) = sxz * (L2 * rho(ix,iy)  &
               * (visc1 * cf + L2 * visc2 * ABS(s)))
          qyz(ix,iy) = syz * (L2 * rho(ix,iy)  &
               * (visc1 * cf + L2 * visc2 * ABS(s)))
          qxx(ix,iy) = sxx * (L2 * rho(ix,iy)  &
               * (visc1 * cf + L2 * visc2 * ABS(s)))
          qyy(ix,iy) = syy * (L2 * rho(ix,iy)  &
               * (visc1 * cf + L2 * visc2 * ABS(s)))
#endif

          flip= - MIN(SIGN(1.0_num,visc3-1.e-6_num),0.0_num)!(visc3 < 1.e-6) ? 1:0
          qxy(ix,iy) = qxy(ix,iy) + sxy * rho(ix,iy) * visc3 * flip 
          qxz(ix,iy) = qxz(ix,iy) + sxz * rho(ix,iy) * visc3 * flip
          qyz(ix,iy) = qyz(ix,iy) + syz * rho(ix,iy) * visc3 * flip
          qxx(ix,iy) = qxx(ix,iy) + sxx * rho(ix,iy) * visc3 * flip  
          qyy(ix,iy) = qyy(ix,iy) + syy * rho(ix,iy) * visc3 * flip  

#ifdef Q_MONO
          IF (visc3 > 1.e-6_num) THEN
#endif
             visc_heat(ix,iy) = qxy(ix,iy)*dvxy + qxz(ix,iy)*dvzdx &
                  + qyz(ix,iy)*dvzdy + qxx(ix,iy)*dvxdx + qyy(ix,iy)*dvydy 
#ifdef Q_MONO
          END IF
#endif

          w3 = bx1(ix,iy)*dvxdx + by1(ix,iy)*dvxdy
          w4 = bx1(ix,iy)*dvydx + by1(ix,iy)*dvydy
          w5 = bx1(ix,iy)*dvzdx + by1(ix,iy)*dvzdy

          bx1(ix,iy) = (bx1(ix,iy) + w3 * dt2) / (1.0_num + dv)
          by1(ix,iy) = (by1(ix,iy) + w4 * dt2) / (1.0_num + dv)
          bz1(ix,iy) = (bz1(ix,iy) + w5 * dt2) / (1.0_num + dv)

       END DO
    END DO

  END SUBROUTINE viscosity_and_b_update




  SUBROUTINE visc_heating

    REAL(num) :: vxb, vxbm, vyb, vybm, vzb, vzbm
    REAL(num) :: dvxdx, dvydy, dvxy, dvxz, dvyz

    DO iy = 0, ny+1
          !DEC$ IVDEP
          !DEC$ VECTOR ALWAYS
       DO ix = 0, nx+1
          ixp = ix + 1
          iyp = iy + 1
          ixm = ix - 1
          iym = iy - 1
          vxb = (vx1(ix,iy) + vx1(ix,iym)) / 2.0_num     !vx at Ex(i,j)
          vxbm = (vx1(ixm,iy) + vx1(ixm,iym)) / 2.0_num  !vx at Ex(i-1,j)
          vyb = (vy1(ix,iy) + vy1(ixm,iy)) / 2.0_num     !vy at Ey(i,j)
          vybm = (vy1(ix,iym) + vy1(ixm,iym)) / 2.0_num  !vy at Ey(i-1,j)

          dvxdx = (vxb - vxbm) / dxb(ix)
          dvydy = (vyb - vybm) / dyb(iy)

          vxb = (vx1(ix,iy) + vx1(ixm,iy)) / 2.0_num     !vx at Ey(i,j)
          vxbm = (vx1(ix,iym) + vx1(ixm,iym)) / 2.0_num  !vx at Ey(i,j-1)
          vyb = (vy1(ix,iy) + vy1(ix,iym)) / 2.0_num     !vy at Ex(i,j)
          vybm = (vy1(ixm,iy) + vy1(ixm,iym)) / 2.0_num  !vy at Ex(i-1,j)

          dvxy = (vxb - vxbm)/dyb(iy) + (vyb - vybm)/dxb(ix)

          vzb = (vz1(ix,iy) + vz1(ix,iym)) / 2.0_num
          vzbm = (vz1(ixm,iy) + vz1(ixm,iym)) / 2.0_num
          dvxz = (vzb - vzbm) / dxb(ix) 

          vzb = (vz1(ix,iy) + vz1(ixm,iy)) / 2.0_num
          vzbm = (vz1(ix,iym) + vz1(ixm,iym)) / 2.0_num
          dvyz = (vzb - vzbm) / dyb(iy) 

          visc_heat(ix,iy) = qxy(ix,iy)*dvxy + qxz(ix,iy)*dvxz &
               + qyz(ix,iy)*dvyz + qxx(ix,iy) * dvxdx + qyy(ix,iy) * dvydy 
       END DO
    END DO

    visc_heat = MAX(visc_heat, 0.0_num)

  END SUBROUTINE visc_heating



  SUBROUTINE eta_calc

    REAL(num) :: jx, jy, jz, jxxp, jyyp, flux
    REAL(num) :: current, cs, d, r2

    eta = 0.0_num
    lambda_i = lambda0

    DO iy = -1, ny+1
       DO ix = -1, nx+1
          ixp = ix + 1
          iyp = iy + 1
          jx = (bz(ix,iyp) - bz(ix,iy)) / dyc(iy)         ! jx at Ey
          jxxp = (bz(ixp,iyp) - bz(ixp,iy)) / dyc(iy)    ! jx at Ey(i+1,j)
          jy = - (bz(ixp,iy) - bz(ix,iy)) / dxc(ix)       ! jy at Ex
          jyyp = - (bz(ixp,iyp) - bz(ix,iyp)) / dxc(ix)  ! jy at Ex(i,j+1)
          jz = (by(ixp,iy) - by(ix,iy)) / dxc(ix)  &      ! jz at V
               - (bx(ix,iyp) - bx(ix,iy)) / dyc(iy)

          !current at V
          w4 = (jx+jxxp) / 2.0_num
          w5 = (jy+jyyp) / 2.0_num
          w6 = jz
          current = SQRT(w4**2 + w5**2 + w6**2) / rho(ix,iy)
          cs = SQRT(gamma * (gamma - 1.0_num) * energy(ix,iy)) * j_max * 4.0_num
          IF (current > j_max) THEN
             eta(ix,iy) = eta0 
             lambda_i(ix,iy) = lambda0
          ELSE
             eta(ix,iy) = eta_background
             lambda_i(ix,iy) = 0.0_num
          ENDIF
       END DO
    END DO


    IF (.NOT. resistiveMHD) eta = 0.0_num

  END SUBROUTINE eta_calc



  SUBROUTINE resistive_effects

    REAL(num) :: jx, jy, jz, jxxp, jyyp, flux
    REAL(num) :: rho_v, half_dt, dt6
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: k1x,k2x,k3x,k4x
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: k1y,k2y,k3y,k4y
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: k1z,k2z,k3z,k4z
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: c1,c2,c3,c4

    ALLOCATE(k1x(0:nx,0:ny),k2x(0:nx,0:ny),k3x(0:nx,0:ny),k4x(0:nx,0:ny))
    ALLOCATE(k1y(0:nx,0:ny),k2y(0:nx,0:ny),k3y(0:nx,0:ny),k4y(0:nx,0:ny))
    ALLOCATE(k1z(0:nx,0:ny),k2z(0:nx,0:ny),k3z(0:nx,0:ny),k4z(0:nx,0:ny))
    ALLOCATE(c1(0:nx,0:ny),c2(0:nx,0:ny),c3(0:nx,0:ny),c4(0:nx,0:ny))

    dt = dt / 2.0_num

    bx1 = bx(0:nx+1,0:ny+1)
    by1 = by(0:nx+1,0:ny+1)
    bz1 = bz(0:nx+1,0:ny+1)

    ! step 1
    CALL rkstep
    k1x = flux_x
    k1y = flux_y
    k1z = flux_z
    c1 = curlb

    ! step 2
    DO iy = 1, ny 
          !DEC$ IVDEP
          !DEC$ VECTOR ALWAYS
       DO ix = 0, nx
          bx(ix,iy) = bx1(ix,iy) + (k1z(ix,iy) - k1z(ix,iy-1)) * dt / dyb(iy)
       END DO
    END DO
    DO iy = 0, ny
          !DEC$ IVDEP
          !DEC$ VECTOR ALWAYS
       DO ix = 1, nx
          by(ix,iy) = by1(ix,iy) + (- k1z(ix,iy) + k1z(ix-1,iy)) * dt / dxb(ix)
       END DO
    END DO
    DO iy = 1, ny
          !DEC$ IVDEP
          !DEC$ VECTOR ALWAYS
       DO ix = 1, nx
          ixm = ix - 1
          iym = iy - 1 
          w1 = - k1x(ix,iy) - k1x(ixm,iy) + k1x(ixm,iym) + k1x(ix,iym)
          w2 = w1 + k1y(ix,iy) + k1y(ix,iym) - k1y(ixm,iy) - k1y(ixm,iym)
          bz(ix,iy) = bz1(ix,iy) + w2 * dt / cv(ix,iy)
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
          !DEC$ IVDEP
          !DEC$ VECTOR ALWAYS
       DO ix = 0, nx
          bx(ix,iy) = bx1(ix,iy) + (k2z(ix,iy) - k2z(ix,iy-1)) * dt / dyb(iy)
       END DO
    END DO
    DO iy = 0, ny
          !DEC$ IVDEP
          !DEC$ VECTOR ALWAYS
       DO ix = 1, nx
          by(ix,iy) = by1(ix,iy) + (- k2z(ix,iy) + k2z(ix-1,iy)) * dt / dxb(ix)
       END DO
    END DO
    DO iy = 1, ny
          !DEC$ IVDEP
          !DEC$ VECTOR ALWAYS
       DO ix = 1, nx
          ixm = ix - 1
          iym = iy - 1 
          w1 = - k2x(ix,iy) - k2x(ixm,iy) + k2x(ixm,iym) + k2x(ix,iym)
          w2 = w1 + k2y(ix,iy) + k2y(ix,iym) - k2y(ixm,iy) - k2y(ixm,iym)
          bz(ix,iy) = bz1(ix,iy) + w2 * dt / cv(ix,iy)
       END DO
    END DO
    CALL bfield_bcs
    CALL rkstep
    k3x = flux_x
    k3y = flux_y
    k3z = flux_z
    c3 = curlb

    dt = 2.0_num * dt

    !step 4
    DO iy = 1, ny 
          !DEC$ IVDEP
          !DEC$ VECTOR ALWAYS
       DO ix = 0, nx
          bx(ix,iy) = bx1(ix,iy) + (k3z(ix,iy) - k3z(ix,iy-1)) * dt / dyb(iy)
       END DO
    END DO
    DO iy = 0, ny
          !DEC$ IVDEP
          !DEC$ VECTOR ALWAYS
       DO ix = 1, nx
          by(ix,iy) = by1(ix,iy) + (- k3z(ix,iy) + k3z(ix-1,iy)) * dt / dxb(ix)
       END DO
    END DO
    DO iy = 1, ny
          !DEC$ IVDEP
          !DEC$ VECTOR ALWAYS
       DO ix = 1, nx
          ixm = ix - 1
          iym = iy - 1 
          w1 = - k3x(ix,iy) - k3x(ixm,iy) + k3x(ixm,iym) + k3x(ix,iym)
          w2 = w1 + k3y(ix,iy) + k3y(ix,iym) - k3y(ixm,iy) - k3y(ixm,iym)
          bz(ix,iy) = bz1(ix,iy) + w2 * dt / cv(ix,iy)
       END DO
    END DO
    CALL bfield_bcs
    CALL rkstep
    k4x = flux_x
    k4y = flux_y
    k4z = flux_z
    c4 = curlb

    !full update
    dt6 = dt / 6.0_num
    k3x = k1x + 2.0_num * k2x + 2.0_num * k3x + k4x
    k3y = k1y + 2.0_num * k2y + 2.0_num * k3y + k4y
    k3z = k1z + 2.0_num * k2z + 2.0_num * k3z + k4z
    c1 = c1 + 2.0_num * c2 + 2.0_num * c3 + c4
    DO iy = 0, ny
          !DEC$ IVDEP
          !DEC$ VECTOR ALWAYS
       DO ix = 1, nx
          by(ix,iy) = by1(ix,iy) + (- k3z(ix,iy) + k3z(ix-1,iy)) * dt6 / dxb(ix)
       END DO
    END DO
    DO iy = 1, ny 
          !DEC$ IVDEP
          !DEC$ VECTOR ALWAYS
       DO ix = 0, nx
          bx(ix,iy) = bx1(ix,iy) + (k3z(ix,iy) - k3z(ix,iy-1)) * dt6 / dyb(iy)
       END DO
    END DO
    DO iy = 1, ny
          !DEC$ IVDEP
          !DEC$ VECTOR ALWAYS
       DO ix = 1, nx
          ixm = ix - 1
          iym = iy - 1 
          w1 = - k3x(ix,iy) - k3x(ixm,iy) + k3x(ixm,iym) + k3x(ix,iym)
          w2 = w1 + k3y(ix,iy) + k3y(ix,iym) - k3y(ixm,iy) - k3y(ixm,iym)
          bz(ix,iy) = bz1(ix,iy) + w2 * dt6 / cv(ix,iy)
       END DO
    END DO
    CALL bfield_bcs

    DO iy = 1, ny
       DO ix = 1, nx
          ixm = ix - 1
          iym = iy - 1
          energy(ix,iy) = energy(ix,iy) + (c1(ix,iy) + c1(ixm,iy) + c1(ix,iym) + c1(ixm,iym)) &
               * dt6 / (4.0_num * rho(ix,iy))
       END DO
    END DO

    CALL energy_bcs

    DO iy = 0, ny
       DO ix = 0, nx
          w1 = dt6 * dxc(ix) * dyc(iy) * c1(ix,iy) 
          IF ((ix == 0) .OR. (ix == nx)) THEN
             w1 = w1 * 0.5_num
          ENDIF
          IF ((iy == 0) .OR. (iy == ny)) THEN
             w1 = w1 * 0.5_num
          ENDIF
          total_ohmic_heating = total_ohmic_heating + w1
       END DO
    END DO

    DEALLOCATE(k1x,k2x,k3x,k4x,k1y,k2y,k3y,k4y,k1z,k2z,k3z,k4z)
    DEALLOCATE(c1,c2,c3,c4)

  END SUBROUTINE resistive_effects



  ! calculates 'k' values from b[xyz]1 values
  SUBROUTINE rkstep

    REAL(num), DIMENSION(:,:), ALLOCATABLE :: jx, jy, jz
    REAL(num), DIMENSION(:,:), ALLOCATABLE  :: hxflux, hyflux, hzflux
    REAL(num) :: jx1, jy1, jz1, jx2, jy2
    REAL(num) :: bxv, byv, bzv, rho_v
    REAL(num) :: f1, f2, area
    INTEGER :: ixp2, iyp2

    ALLOCATE(jx(-1:nx+1,-1:ny+1),jy(-1:nx+1,-1:ny+1),jz(-1:nx+1,-1:ny+1))

    DO iy = -1, ny+1
          !DEC$ IVDEP
          !DEC$ VECTOR ALWAYS
       DO ix = -1, nx+1
          ixp = ix + 1
          iyp = iy + 1
          jx1 = (bz(ix,iyp) - bz(ix,iy)) / dyc(iy)      
          jx2 = (bz(ixp,iyp) - bz(ixp,iy)) / dyc(iy) 
          jy1 = - (bz(ixp,iy) - bz(ix,iy)) / dxc(ix)   
          jy2 = - (bz(ixp,iyp) - bz(ix,iyp)) / dxc(ix)
          jx(ix,iy) = (jx1 + jx2) / 2.0_num    
          jy(ix,iy) = (jy1 + jy2) / 2.0_num    
          jz(ix,iy) = (by(ixp,iy) - by(ix,iy)) / dxc(ix)  &  
               - (bx(ix,iyp) - bx(ix,iy)) / dyc(iy)
       END DO
    END DO

    DO iy = 0, ny
          !DEC$ IVDEP
          !DEC$ VECTOR ALWAYS
       DO ix = 0, nx
          flux_x(ix,iy) = - jx(ix,iy) * eta(ix,iy) * dxc(ix) / 2.0_num
          flux_y(ix,iy) = - jy(ix,iy) * eta(ix,iy) * dyc(iy) / 2.0_num
          flux_z(ix,iy) = - jz(ix,iy) * eta(ix,iy)
          curlb(ix,iy) = eta(ix,iy) * (jx(ix,iy)**2 + jy(ix,iy)**2 + jz(ix,iy)**2)
       END DO
    END DO

    DEALLOCATE (jx, jy, jz)

  END SUBROUTINE rkstep



  SUBROUTINE store_boundary_dv

    REAL(num) :: dvx, dvy, dvz

    IF (xbc_right == open .AND. right == MPI_PROC_NULL) THEN
       DO iy = -2,ny+2
          dvx=2.0_num*(vx(nx,iy)-vx1(nx,iy))
          dvy=2.0_num*(vy(nx,iy)-vy1(nx,iy))
          dvz=2.0_num*(vz(nx,iy)-vz1(nx,iy))
          dv_right(iy)=SQRT(dvx**2+dvy**2+dvz**2)
       END DO
    END IF
    IF (xbc_left == open .AND. left == MPI_PROC_NULL) THEN
       DO iy = -2,ny+2
          dvx=2.0_num*(vx(0,iy)-vx1(0,iy))
          dvy=2.0_num*(vy(0,iy)-vy1(0,iy))
          dvz=2.0_num*(vz(0,iy)-vz1(0,iy))
          dv_left(iy)=SQRT(dvx**2+dvy**2+dvz**2)
       END DO
    END IF
    IF (ybc_up == open .AND. up == MPI_PROC_NULL) THEN
       DO ix = -2,nx+2
          dvx=2.0_num*(vx(ix,ny)-vx1(ix,ny))
          dvy=2.0_num*(vy(ix,ny)-vy1(ix,ny))
          dvz=2.0_num*(vz(ix,ny)-vz1(ix,ny))
          dv_up(ix)=SQRT(dvx**2+dvy**2+dvz**2)
       END DO
    END IF
    IF (ybc_down == open .AND. down == MPI_PROC_NULL) THEN
       DO ix = -2,nx+2
          dvx=2.0_num*(vx(ix,0)-vx1(ix,0))
          dvy=2.0_num*(vy(ix,0)-vy1(ix,0))
          dvz=2.0_num*(vz(ix,0)-vz1(ix,0))
          dv_down(ix)=SQRT(dvx**2+dvy**2+dvz**2)
       END DO
    END IF

  END SUBROUTINE store_boundary_dv


END MODULE lagran
  
