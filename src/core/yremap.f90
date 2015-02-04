!******************************************************************************
! Mass coordinate based Van Leer limited remap.
! See Bram van Leer, JCP, vol 135, p229, (1997)
!******************************************************************************

MODULE yremap

  USE shared_data; USE boundary

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: remap_y

  REAL(num), DIMENSION(:,:), ALLOCATABLE :: dyb1, dyc1, rho_v, rho_v1
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: rho1, dm, cv2, flux

CONTAINS

  ! Remap onto original Eulerian grid

  SUBROUTINE remap_y

    REAL(num) :: vxb, vxbm, vyb, vybm, dv

    ALLOCATE(rho1  (-1:nx+2, -1:ny+2))
    ALLOCATE(dm    (-1:nx+2, -1:ny+2))
    ALLOCATE(cv2   (-1:nx+2, -1:ny+2))
    ALLOCATE(flux  (-1:nx+2, -2:ny+2))
    ALLOCATE(dyb1  (-1:nx+2, -1:ny+2))
    ALLOCATE(dyc1  (-1:nx+1, -1:ny+1))
    ALLOCATE(rho_v (-1:nx+2, -1:ny+2))
    ALLOCATE(rho_v1(-1:nx+2, -1:ny+2))

    dm = 0.0_num
    ! Store initial density in rho1
    rho1 = rho

    DO iy = -1, ny + 2
      iym = iy - 1
      DO ix = -1, nx + 2
        ixm = ix - 1

        ! vx at Bx(i,j)
        vxb  = (vx1(ix ,iy ) + vx1(ix ,iym)) * 0.5_num

        ! vx at Bx(i-1,j)
        vxbm = (vx1(ixm,iy ) + vx1(ixm,iym)) * 0.5_num

        ! vy at By(i,j)
        vyb  = (vy1(ix ,iy ) + vy1(ixm,iy )) * 0.5_num

        ! vy at By(i,j-1)
        vybm = (vy1(ix ,iym) + vy1(ixm,iym)) * 0.5_num

        dv = (REAL(xpass, num) * (vxb - vxbm) / dxb(ix) &
            + (vyb - vybm) / dyb(iy)) * dt

        ! Control volume before remap
        cv1(ix,iy) = cv(ix,iy) * (1.0_num + dv)

        dv = REAL(xpass, num) * (vxb - vxbm) / dxb(ix) * dt

        ! Control volume after remap
        cv2(ix,iy) = cv(ix,iy) * (1.0_num + dv)

        ! dyb before remap
        dyb1(ix,iy) = dyb(iy) + (vyb - vybm) * dt
      END DO
    END DO

    DO iy = -1, ny + 1
      iyp = iy + 1
      DO ix = -1, nx + 1
        ! dyc before remap
        dyc1(ix,iy) = 0.5_num * (dyb1(ix,iy) + dyb1(ix,iyp))
      END DO
    END DO

    ! Evans and Hawley (ApJ, vol 332, p650, (1988))
    ! constrained transport remap of magnetic fluxes
    CALL vy_bx_flux

    DO iy = 1, ny
      iym = iy - 1
      DO ix = 0, nx
        bx(ix,iy) = bx(ix,iy) - flux(ix,iy) + flux(ix,iym)
      END DO
    END DO

    DO iy = 0, ny
      DO ix = 1, nx
        ixm = ix - 1
        by(ix,iy) = by(ix,iy) + flux(ix,iy) - flux(ixm,iy)
      END DO
    END DO

    CALL vy_bz_flux

    DO iy = 1, ny
      iym = iy - 1
      DO ix = 1, nx
        bz(ix,iy) = bz(ix,iy) - flux(ix,iy) + flux(ix,iym)
      END DO
    END DO

    ! Remap of mass + calculation of mass fluxes (dm) needed for later remaps
    ! Calculates dm(0:nx+1,0:ny)
    CALL y_mass_flux
    ! Need dm(0:nx+1,-1:ny+1) for velocity remap
    CALL dm_y_bcs

    DO iy = 1, ny
      iym = iy - 1
      DO ix = 1, nx
        rho(ix,iy) = (rho1(ix,iy) * cv1(ix,iy) &
            + dm(ix,iym) - dm(ix,iy)) / cv2(ix,iy)
      END DO
    END DO

    ! Remap specific energy density using mass coordinates
    CALL y_energy_flux

    DO iy = 1, ny
      iym = iy - 1
      DO ix = 1, nx
        energy(ix,iy) = (energy(ix,iy) * cv1(ix,iy) &
            * rho1(ix,iy) + flux(ix,iym) - flux(ix,iy)) &
            / (cv2(ix,iy) * rho(ix,iy))
      END DO
    END DO

    ! Redefine dyb1, cv1, cv2, dm and vy1 for velocity (vertex) cells.
    ! In some of these calculations the flux variable is used as a
    ! temporary array
    DO iy = -1, ny + 1
      iyp = iy + 1
      DO ix = 0, nx
        ixp = ix + 1

        ! Vertex density before remap
        rho_v(ix,iy) = rho1(ix,iy) * cv1(ix,iy) &
            + rho1(ixp,iy ) * cv1(ixp,iy ) &
            + rho1(ix ,iyp) * cv1(ix ,iyp) &
            + rho1(ixp,iyp) * cv1(ixp,iyp)

        rho_v(ix,iy) = rho_v(ix,iy) &
            / (cv1(ix,iy ) + cv1(ixp,iy ) + cv1(ix,iyp) + cv1(ixp,iyp))
      END DO
    END DO

    DO iy = 0, ny
      iyp = iy + 1
      DO ix = 0, nx
        ixp = ix + 1
        flux(ix,iy) = cv1(ix,iy) + cv1(ixp,iy) + cv1(ix,iyp) + cv1(ixp,iyp)
      END DO
    END DO

    ! cv1 = vertex CV before remap
    cv1(0:nx,0:ny) = flux(0:nx,0:ny) * 0.25_num

    DO iy = 0, ny
      iyp = iy + 1
      DO ix = 0, nx
        ixp = ix + 1
        flux(ix,iy) = cv2(ix,iy) + cv2(ixp,iy) + cv2(ix,iyp) + cv2(ixp,iyp)
      END DO
    END DO

    ! cv2 = vertex CV after remap
    cv2(0:nx,0:ny) = flux(0:nx,0:ny) * 0.25_num

    DO iy = -2, ny + 1
      iyp = iy + 1
      DO ix = 0, nx
        flux(ix,iy) = (vy1(ix,iy) + vy1(ix,iyp)) * 0.5_num
      END DO
    END DO

    ! Vertex boundary velocity used in remap
    vy1(0:nx,-2:ny+1) = flux(0:nx,-2:ny+1)

    ! Calculate vertex-centred lengths

    DO iy = -1, ny + 2
      iym = iy - 1
      DO ix = -1, nx + 2
        ! dyb before remap
        dyb1(ix,iy) = dyb(iy) + (vy1(ix,iy) - vy1(ix,iym)) * dt
      END DO
    END DO

    DO iy = -1, ny + 1
      iyp = iy + 1
      DO ix = -1, nx + 1
        ! dyc before remap
        dyc1(ix,iy) = 0.5_num * (dyb1(ix,iy) + dyb1(ix,iyp))
      END DO
    END DO

    DO iy = -1, ny
      iyp = iy + 1
      DO ix = 0, nx
        ixp = ix + 1
        flux(ix,iy) = dm(ix,iy) + dm(ixp,iy) + dm(ix,iyp) + dm(ixp,iyp)
      END DO
    END DO

    ! Mass flux out of vertex CV
    dm(0:nx,-1:ny) = flux(0:nx,-1:ny) * 0.25_num

    DO iy = 0, ny
      iym = iy - 1
      DO ix = 0, nx
        ! Vertex density after remap
        rho_v1(ix,iy) = (rho_v(ix,iy) * cv1(ix,iy) &
            + dm(ix,iym) - dm(ix,iy)) / cv2(ix,iy)
      END DO
    END DO

    CALL y_momx_flux

    DO iy = 0, ny
      iym = iy - 1
      DO ix = 0, nx
        vx(ix,iy) = (rho_v(ix,iy) * vx(ix,iy) * cv1(ix,iy) &
            + flux(ix,iym) - flux(ix,iy)) / (cv2(ix,iy) * rho_v1(ix,iy))
      END DO
    END DO

    CALL y_momy_flux

    DO iy = 0, ny
      iym = iy - 1
      DO ix = 0, nx
        vy(ix,iy) = (rho_v(ix,iy) * vy(ix,iy) * cv1(ix,iy) &
            + flux(ix,iym) - flux(ix,iy)) / (cv2(ix,iy) * rho_v1(ix,iy))
      END DO
    END DO

    CALL y_momz_flux

    DO iy = 0, ny
      iym = iy - 1
      DO ix = 0, nx
        vz(ix,iy) = (rho_v(ix,iy) * vz(ix,iy) * cv1(ix,iy) &
            + flux(ix,iym) - flux(ix,iy)) / (cv2(ix,iy) * rho_v1(ix,iy))
      END DO
    END DO

    CALL boundary_conditions

    DEALLOCATE (rho1, dm, cv2, flux, dyb1, dyc1, rho_v, rho_v1)
    ypass = 0

  END SUBROUTINE remap_y



  SUBROUTINE vy_bx_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: fm, fi, fp, fp2, fu, dfm, dfi, dfp, dfu, sign_v
    REAL(num) :: dyci, dycu, dybu, phi, ss, Da, Di
    REAL(num) :: dyu, dby, dbyp, dbyp2, dbym
    INTEGER :: iyp2

    DO iy = 0, ny
      iym  = iy - 1
      iyp  = iy + 1
      iyp2 = iy + 2
      DO ix = 0, nx
        ixp = ix + 1

        v_advect = vy1(ix,iy)

        dby   = (dyb1(ix,iy  ) + dyb1(ixp,iy  )) * 0.5_num
        dbyp  = (dyb1(ix,iyp ) + dyb1(ixp,iyp )) * 0.5_num
        dbyp2 = (dyb1(ix,iyp2) + dyb1(ixp,iyp2)) * 0.5_num
        dbym  = (dyb1(ix,iym ) + dyb1(ixp,iym )) * 0.5_num

        fm  = bx(ix,iym ) / dbym
        fi  = bx(ix,iy  ) / dby
        fp  = bx(ix,iyp ) / dbyp
        fp2 = bx(ix,iyp2) / dbyp2

        dfm = fi - fm
        dfi = fp - fi
        dfp = fp2 - fp

        ! vad_p and vad_m are logical switches which determine v_advect>=0
        ! and v_advect<0 respectively. It's written this way to allow vector
        ! optimization

        sign_v = SIGN(1.0_num, v_advect)
        vad_p = (sign_v + 1.0_num) * 0.5_num
        vad_m = 1.0_num - vad_p

        fu = fi * vad_p + fp * vad_m
        dfu = dfm * vad_p + dfp * vad_m
        dyci = dyc1(ix,iy )
        dycu = dyc1(ix,iym) * vad_p + dyc1(ix,iyp) * vad_m
        dybu = dyb1(ix,iy ) * vad_p + dyb1(ix,iyp) * vad_m

        dyu = dby * vad_p + dbyp * vad_m
        phi = ABS(v_advect) * dt / dyu

        Da =  (2.0_num - phi) * ABS(dfi) / dyci &
            + (1.0_num + phi) * ABS(dfu) / dycu
        Da = Da * sixth

        ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

        Di = sign_v * ss * MIN(ABS(Da) * dybu, ABS(dfi), ABS(dfu))

        flux(ix,iy) = (fu + Di * (1.0_num - phi)) * v_advect * dt
      END DO
    END DO

  END SUBROUTINE vy_bx_flux



  SUBROUTINE vy_bz_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: fm, fi, fp, fp2, fu, dfm, dfi, dfp, dfu, sign_v
    REAL(num) :: dyci, dycu, dybu, phi, ss, Da, Di
    REAL(num) :: dyu, dby, dbyp, dbyp2, dbym
    INTEGER :: iyp2

    DO iy = 0, ny
      iym  = iy - 1
      iyp  = iy + 1
      iyp2 = iy + 2
      DO ix = 0, nx
        ixm = ix - 1

        v_advect = (vy1(ix,iy) + vy1(ixm,iy)) * 0.5_num

        dby   = dyb1(ix,iy  )
        dbyp  = dyb1(ix,iyp )
        dbyp2 = dyb1(ix,iyp2)
        dbym  = dyb1(ix,iym )

        fm  = bz(ix,iym ) / dbym
        fi  = bz(ix,iy  ) / dby
        fp  = bz(ix,iyp ) / dbyp
        fp2 = bz(ix,iyp2) / dbyp2

        dfm = fi - fm
        dfi = fp - fi
        dfp = fp2 - fp

        ! vad_p and vad_m are logical switches which determine v_advect>=0
        ! and v_advect<0 respectively. It's written this way to allow vector
        ! optimization

        sign_v = SIGN(1.0_num, v_advect)
        vad_p = (sign_v + 1.0_num) * 0.5_num
        vad_m = 1.0_num - vad_p

        fu = fi * vad_p + fp * vad_m
        dfu = dfm * vad_p + dfp * vad_m
        dyci = dyc1(ix,iy )
        dycu = dyc1(ix,iym) * vad_p + dyc1(ix,iyp) * vad_m
        dybu = dyb1(ix,iy ) * vad_p + dyb1(ix,iyp) * vad_m

        dyu = dby * vad_p + dbyp * vad_m
        phi = ABS(v_advect) * dt / dyu

        Da =  (2.0_num - phi) * ABS(dfi) / dyci &
            + (1.0_num + phi) * ABS(dfu) / dycu
        Da = Da * sixth

        ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

        Di = sign_v * ss * MIN(ABS(Da) * dybu, ABS(dfi), ABS(dfu))

        flux(ix,iy) = (fu + Di * (1.0_num - phi)) * v_advect * dt
      END DO
    END DO

  END SUBROUTINE vy_bz_flux



  SUBROUTINE y_mass_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: fm, fi, fp, fp2, fu, dfm, dfi, dfp, dfu, sign_v
    REAL(num) :: dyci, dycu, dybu, phi, ss, Da, Di
    REAL(num) :: area
    INTEGER :: iyp2

    DO iy = 0, ny
      iym  = iy - 1
      iyp  = iy + 1
      iyp2 = iy + 2
      DO ix = 0, nx + 1
        ixm = ix - 1
        area = dxb(ix)

        v_advect = (vy1(ix,iy) + vy1(ixm,iy)) * 0.5_num

        fm  = rho(ix,iym )
        fi  = rho(ix,iy  )
        fp  = rho(ix,iyp )
        fp2 = rho(ix,iyp2)

        dfm = fi - fm
        dfi = fp - fi
        dfp = fp2 - fp

        ! vad_p and vad_m are logical switches which determine v_advect>=0
        ! and v_advect<0 respectively. It's written this way to allow vector
        ! optimization

        sign_v = SIGN(1.0_num, v_advect)
        vad_p = (sign_v + 1.0_num) * 0.5_num
        vad_m = 1.0_num - vad_p

        fu = fi * vad_p + fp * vad_m
        dfu = dfm * vad_p + dfp * vad_m
        dyci = dyc1(ix,iy )
        dycu = dyc1(ix,iym) * vad_p + dyc1(ix,iyp) * vad_m
        dybu = dyb1(ix,iy ) * vad_p + dyb1(ix,iyp) * vad_m

        phi = ABS(v_advect) * dt / dybu

        Da =  (2.0_num - phi) * ABS(dfi) / dyci &
            + (1.0_num + phi) * ABS(dfu) / dycu
        Da = Da * sixth

        ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

        Di = sign_v * ss * MIN(ABS(Da) * dybu, ABS(dfi), ABS(dfu))

        dm(ix,iy) = (fu + Di * (1.0_num - phi)) * v_advect * dt * area
      END DO
    END DO

  END SUBROUTINE y_mass_flux



  ! Energy remap in mass coordinates

  SUBROUTINE y_energy_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: fm, fi, fp, fp2, fu, dfm, dfi, dfp, dfu, sign_v
    REAL(num) :: dyci, dycu, dybu, phi, ss, Da, Di
    REAL(num) :: area, rhou, dmu
    INTEGER :: iyp2

    DO iy = 0, ny
      iym  = iy - 1
      iyp  = iy + 1
      iyp2 = iy + 2
      DO ix = 0, nx
        ixm = ix - 1
        area = dxb(ix)

        v_advect = (vy1(ix,iy) + vy1(ixm,iy)) * 0.5_num

        fm  = energy(ix,iym )
        fi  = energy(ix,iy  )
        fp  = energy(ix,iyp )
        fp2 = energy(ix,iyp2)

        dfm = fi - fm
        dfi = fp - fi
        dfp = fp2 - fp

        ! vad_p and vad_m are logical switches which determine v_advect>=0
        ! and v_advect<0 respectively. It's written this way to allow vector
        ! optimization

        sign_v = SIGN(1.0_num, v_advect)
        vad_p = (sign_v + 1.0_num) * 0.5_num
        vad_m = 1.0_num - vad_p

        fu = fi * vad_p + fp * vad_m
        dfu = dfm * vad_p + dfp * vad_m
        dyci = dyc1(ix,iy )
        dycu = dyc1(ix,iym) * vad_p + dyc1(ix,iyp) * vad_m
        dybu = dyb1(ix,iy ) * vad_p + dyb1(ix,iyp) * vad_m

        phi = ABS(v_advect) * dt / dybu

        Da =  (2.0_num - phi) * ABS(dfi) / dyci &
            + (1.0_num + phi) * ABS(dfu) / dycu
        Da = Da * sixth

        ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

        Di = sign_v * ss * MIN(ABS(Da) * dybu, ABS(dfi), ABS(dfu))

        rhou = rho1(ix,iy) * vad_p + rho1(ix,iyp) * vad_m
        dmu = ABS(dm(ix,iy)) / area / dybu / rhou

        flux(ix,iy) = (fu + Di * (1.0_num - dmu)) * dm(ix,iy)
      END DO
    END DO

  END SUBROUTINE y_energy_flux



  SUBROUTINE y_momx_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: fm, fi, fp, fp2, fu, dfm, dfi, dfp, dfu, sign_v
    REAL(num) :: dyci, dycu, dybu, phi, ss, Da, Di
    REAL(num) :: area, rhou, dmu, m, mp, ai, aip, dk
    INTEGER :: iyp2

    DO iy = -1, ny
      iym  = iy - 1
      iyp  = iy + 1
      iyp2 = iy + 2
      DO ix = 0, nx
        area = dxc(ix)

        v_advect = vy1(ix,iy)

        fm  = vx(ix,iym )
        fi  = vx(ix,iy  )
        fp  = vx(ix,iyp )
        fp2 = vx(ix,iyp2)

        dfm = fi - fm
        dfi = fp - fi
        dfp = fp2 - fp

        ! vad_p and vad_m are logical switches which determine v_advect>=0
        ! and v_advect<0 respectively. It's written this way to allow vector
        ! optimization

        sign_v = SIGN(1.0_num, v_advect)
        vad_p = (sign_v + 1.0_num) * 0.5_num
        vad_m = 1.0_num - vad_p

        fu = fi * vad_p + fp * vad_m
        dfu = dfm * vad_p + dfp * vad_m
        dyci = dyb1(ix,iyp)
        dycu = dyb1(ix,iy ) * vad_p + dyb1(ix,iyp2) * vad_m
        dybu = dyc1(ix,iy ) * vad_p + dyc1(ix,iyp ) * vad_m

        phi = ABS(v_advect) * dt / dybu

        Da =  (2.0_num - phi) * ABS(dfi) / dyci &
            + (1.0_num + phi) * ABS(dfu) / dycu
        Da = Da * sixth

        ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

        Di = sign_v * ss * MIN(ABS(Da) * dybu, ABS(dfi), ABS(dfu))

        rhou = rho_v(ix,iy) * vad_p + rho_v(ix,iyp) * vad_m
        dmu = ABS(dm(ix,iy)) / area / dybu / rhou

        flux(ix,iy) = fu + Di * (1.0_num - dmu)
      END DO
    END DO

    IF (rke) THEN
      DO iy = 0, ny - 1
        iym = iy - 1
        iyp = iy + 1
        DO ix = 0, nx
          ixp = ix + 1

          m  = rho_v1(ix,iy ) * cv2(ix,iy )
          mp = rho_v1(ix,iyp) * cv2(ix,iyp)

          ai =  (vx(ix,iy ) - flux(ix,iym)) * dm(ix,iym) / m &
              - (vx(ix,iy ) - flux(ix,iy )) * dm(ix,iy ) / m

          aip = (vx(ix,iyp) - flux(ix,iy )) * dm(ix,iy ) / mp &
              - (vx(ix,iyp) - flux(ix,iyp)) * dm(ix,iyp) / mp

          dk = (vx(ix,iyp) - vx(ix,iy)) * (flux(ix,iy) &
              - 0.5_num * (vx(ix,iyp) + vx(ix,iy))) &
              - 0.5_num * ai  * (vx(ix,iy ) - flux(ix,iy)) &
              + 0.5_num * aip * (vx(ix,iyp) - flux(ix,iy))

          dk = dk * dm(ix,iy) * 0.5_num
          delta_ke(ix ,iyp) = delta_ke(ix ,iyp) + dk
          delta_ke(ixp,iyp) = delta_ke(ixp,iyp) + dk
        END DO
      END DO
    END IF

    flux(0:nx,-1:ny) = flux(0:nx,-1:ny) * dm(0:nx,-1:ny)

  END SUBROUTINE y_momx_flux



  SUBROUTINE y_momy_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: fm, fi, fp, fp2, fu, dfm, dfi, dfp, dfu, sign_v
    REAL(num) :: dyci, dycu, dybu, phi, ss, Da, Di
    REAL(num) :: area, rhou, dmu, m, mp, ai, aip, dk
    INTEGER :: iyp2

    DO iy = -1, ny
      iym  = iy - 1
      iyp  = iy + 1
      iyp2 = iy + 2
      DO ix = 0, nx
        area = dxc(ix)

        v_advect = vy1(ix,iy)

        fm  = vy(ix,iym )
        fi  = vy(ix,iy  )
        fp  = vy(ix,iyp )
        fp2 = vy(ix,iyp2)

        dfm = fi - fm
        dfi = fp - fi
        dfp = fp2 - fp

        ! vad_p and vad_m are logical switches which determine v_advect>=0
        ! and v_advect<0 respectively. It's written this way to allow vector
        ! optimization

        sign_v = SIGN(1.0_num, v_advect)
        vad_p = (sign_v + 1.0_num) * 0.5_num
        vad_m = 1.0_num - vad_p

        fu = fi * vad_p + fp * vad_m
        dfu = dfm * vad_p + dfp * vad_m
        dyci = dyb1(ix,iyp)
        dycu = dyb1(ix,iy ) * vad_p + dyb1(ix,iyp2) * vad_m
        dybu = dyc1(ix,iy ) * vad_p + dyc1(ix,iyp ) * vad_m

        phi = ABS(v_advect) * dt / dybu

        Da =  (2.0_num - phi) * ABS(dfi) / dyci &
            + (1.0_num + phi) * ABS(dfu) / dycu
        Da = Da * sixth

        ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

        Di = sign_v * ss * MIN(ABS(Da) * dybu, ABS(dfi), ABS(dfu))

        rhou = rho_v(ix,iy) * vad_p + rho_v(ix,iyp) * vad_m
        dmu = ABS(dm(ix,iy)) / area / dybu / rhou

        flux(ix,iy) = fu + Di * (1.0_num - dmu)
      END DO
    END DO

    IF (rke) THEN
      DO iy = 0, ny - 1
        iym = iy - 1
        iyp = iy + 1
        DO ix = 0, nx
          ixp = ix + 1

          m  = rho_v1(ix,iy ) * cv2(ix,iy )
          mp = rho_v1(ix,iyp) * cv2(ix,iyp)

          ai =  (vy(ix,iy ) - flux(ix,iym)) * dm(ix,iym) / m &
              - (vy(ix,iy ) - flux(ix,iy )) * dm(ix,iy ) / m

          aip = (vy(ix,iyp) - flux(ix,iy )) * dm(ix,iy ) / mp &
              - (vy(ix,iyp) - flux(ix,iyp)) * dm(ix,iyp) / mp

          dk = (vy(ix,iyp) - vy(ix,iy)) * (flux(ix,iy) &
              - 0.5_num * (vy(ix,iyp) + vy(ix,iy))) &
              - 0.5_num * ai  * (vy(ix,iy ) - flux(ix,iy)) &
              + 0.5_num * aip * (vy(ix,iyp) - flux(ix,iy))

          dk = dk * dm(ix,iy) * 0.5_num
          delta_ke(ix ,iyp) = delta_ke(ix ,iyp) + dk
          delta_ke(ixp,iyp) = delta_ke(ixp,iyp) + dk
        END DO
      END DO
    END IF

    flux(0:nx,-1:ny) = flux(0:nx,-1:ny) * dm(0:nx,-1:ny)

  END SUBROUTINE y_momy_flux



  SUBROUTINE y_momz_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: fm, fi, fp, fp2, fu, dfm, dfi, dfp, dfu, sign_v
    REAL(num) :: dyci, dycu, dybu, phi, ss, Da, Di
    REAL(num) :: area, rhou, dmu, m, mp, ai, aip, dk
    INTEGER :: iyp2

    DO iy = -1, ny
      iym  = iy - 1
      iyp  = iy + 1
      iyp2 = iy + 2
      DO ix = 0, nx
        area = dxc(ix)

        v_advect = vy1(ix,iy)

        fm  = vz(ix,iym )
        fi  = vz(ix,iy  )
        fp  = vz(ix,iyp )
        fp2 = vz(ix,iyp2)

        dfm = fi - fm
        dfi = fp - fi
        dfp = fp2 - fp

        ! vad_p and vad_m are logical switches which determine v_advect>=0
        ! and v_advect<0 respectively. It's written this way to allow vector
        ! optimization

        sign_v = SIGN(1.0_num, v_advect)
        vad_p = (sign_v + 1.0_num) * 0.5_num
        vad_m = 1.0_num - vad_p

        fu = fi * vad_p + fp * vad_m
        dfu = dfm * vad_p + dfp * vad_m
        dyci = dyb1(ix,iyp)
        dycu = dyb1(ix,iy ) * vad_p + dyb1(ix,iyp2) * vad_m
        dybu = dyc1(ix,iy ) * vad_p + dyc1(ix,iyp ) * vad_m

        phi = ABS(v_advect) * dt / dybu

        Da =  (2.0_num - phi) * ABS(dfi) / dyci &
            + (1.0_num + phi) * ABS(dfu) / dycu
        Da = Da * sixth

        ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

        Di = sign_v * ss * MIN(ABS(Da) * dybu, ABS(dfi), ABS(dfu))

        rhou = rho_v(ix,iy) * vad_p + rho_v(ix,iyp) * vad_m
        dmu = ABS(dm(ix,iy)) / area / dybu / rhou

        flux(ix,iy) = fu + Di * (1.0_num - dmu)
      END DO
    END DO

    IF (rke) THEN
      DO iy = 0, ny - 1
        iym = iy - 1
        iyp = iy + 1
        DO ix = 0, nx
          ixp = ix + 1

          m  = rho_v1(ix,iy ) * cv2(ix,iy )
          mp = rho_v1(ix,iyp) * cv2(ix,iyp)

          ai =  (vz(ix,iy ) - flux(ix,iym)) * dm(ix,iym) / m &
              - (vz(ix,iy ) - flux(ix,iy )) * dm(ix,iy ) / m

          aip = (vz(ix,iyp) - flux(ix,iy )) * dm(ix,iy ) / mp &
              - (vz(ix,iyp) - flux(ix,iyp)) * dm(ix,iyp) / mp

          dk = (vz(ix,iyp) - vz(ix,iy)) * (flux(ix,iy) &
              - 0.5_num * (vz(ix,iyp) + vz(ix,iy))) &
              - 0.5_num * ai  * (vz(ix,iy ) - flux(ix,iy)) &
              + 0.5_num * aip * (vz(ix,iyp) - flux(ix,iy))

          dk = dk * dm(ix,iy) * 0.5_num
          delta_ke(ix ,iyp) = delta_ke(ix ,iyp) + dk
          delta_ke(ixp,iyp) = delta_ke(ixp,iyp) + dk
        END DO
      END DO
    END IF

    flux(0:nx,-1:ny) = flux(0:nx,-1:ny) * dm(0:nx,-1:ny)

  END SUBROUTINE y_momz_flux



  SUBROUTINE dm_y_bcs

    CALL MPI_SENDRECV(&
        dm(-1,1   ), 1, cell_yface, proc_y_min, tag, &
        dm(-1,ny+1), 1, cell_yface, proc_y_max, tag, &
        comm, status, errcode)

    IF (proc_y_max == MPI_PROC_NULL) &
        dm(0:nx+1,ny+1) = dm(0:nx+1,ny)

    CALL MPI_SENDRECV(&
        dm(-1,ny-1), 1, cell_yface, proc_y_max, tag, &
        dm(-1,-1  ), 1, cell_yface, proc_y_min, tag, &
        comm, status, errcode)

    IF (proc_y_min == MPI_PROC_NULL) &
        dm(0:nx+1,-1) = dm(0:nx+1,0)

  END SUBROUTINE dm_y_bcs

END MODULE yremap
