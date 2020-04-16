  ! Copyright 2020 University of Warwick

  ! Licensed under the Apache License, Version 2.0 (the "License");
  ! you may not use this file except in compliance with the License.
  ! You may obtain a copy of the License at

  !    http://www.apache.org/licenses/LICENSE-2.0

  ! Unless required by applicable law or agreed to in writing, software
  ! distributed under the License is distributed on an "AS IS" BASIS,
  ! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  ! See the License for the specific language governing permissions and
  ! limitations under the License.
  
!******************************************************************************
! Mass coordinate based Van Leer limited remap.
! See Bram van Leer, JCP, vol 135, p229, (1997)
!******************************************************************************

MODULE xremap

  USE shared_data; USE boundary

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: remap_x

  REAL(num), DIMENSION(:,:), ALLOCATABLE :: dxb1, dxc1, rho_v, rho_v1
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: rho1, dm, cv2, flux

CONTAINS

  ! Remap onto original Eulerian grid

  SUBROUTINE remap_x

    REAL(num) :: vxb, vxbm, vyb, vybm, dv

    IF (predictor_step) THEN
      ALLOCATE(cv2   (-1:nx+2, -1:ny+2))
      ALLOCATE(flux  (-2:nx+2, -1:ny+2))
      ALLOCATE(dxb1  (-1:nx+2, -1:ny+2))
      ALLOCATE(dxc1  (-1:nx+2, -1:ny+2))
    ELSE
      ALLOCATE(rho1  (-1:nx+2, -1:ny+2))
      ALLOCATE(dm    (-1:nx+2, -1:ny+2))
      ALLOCATE(cv2   (-1:nx+2, -1:ny+2))
      ALLOCATE(flux  (-2:nx+2, -1:ny+2))
      ALLOCATE(dxb1  (-1:nx+2, -1:ny+2))
      ALLOCATE(dxc1  (-1:nx+2, -1:ny+2))
      ALLOCATE(rho_v (-1:nx+2, -1:ny+2))
      ALLOCATE(rho_v1(-1:nx+2, -1:ny+2))
    END IF

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

        dv = (REAL(ypass, num) * (vyb - vybm) / dyb(iy) &
            + (vxb - vxbm) / dxb(ix)) * dt

        ! Control volume before remap
        cv1(ix,iy) = cv(ix,iy) * (1.0_num + dv)

        dv = REAL(ypass, num) * (vyb - vybm) / dyb(iy) * dt

        ! Control volume after remap
        cv2(ix,iy) = cv(ix,iy) * (1.0_num + dv)

        ! dxb before remap
        dxb1(ix,iy) = dxb(ix) + (vxb - vxbm) * dt
      END DO
    END DO

    DO iy = -1, ny + 2
      DO ix = -1, nx + 2
        ixp = MIN(ix + 1, nx + 2)
        ! dxc before remap
        dxc1(ix,iy) = 0.5_num * (dxb1(ix,iy) + dxb1(ixp,iy))
      END DO
    END DO

    ! Evans and Hawley (ApJ, vol 332, p650, (1988))
    ! constrained transport remap of magnetic fluxes
    CALL vx_by_flux

    DO iy = 0, ny + 1
      DO ix = 0, nx + 1
        ixm = ix - 1
        by(ix,iy) = by(ix,iy) - flux(ix,iy) + flux(ixm,iy)
      END DO
    END DO

    DO iy = 0, ny + 1
      iym = iy - 1
      DO ix = 0, nx + 1
        bx(ix,iy) = bx(ix,iy) + flux(ix,iy) - flux(ix,iym)
      END DO
    END DO

    CALL vx_bz_flux

    DO iy = 0, ny + 1
      DO ix = 0, nx +1
        ixm = ix - 1
        bz(ix,iy) = bz(ix,iy) - flux(ix,iy) + flux(ixm,iy)
      END DO
    END DO

    IF (predictor_step) THEN
      DEALLOCATE(cv2, flux, dxb1, dxc1)
      RETURN
    END IF

    dm = 0.0_num
    ! Store initial density in rho1
    rho1(:,:) = rho(:,:)

    ! Remap of mass + calculation of mass fluxes (dm) needed for later remaps
    ! Calculates dm(0:nx,0:ny+1)
    CALL x_mass_flux
    ! Need dm(-1:nx+1,0:ny+1) for velocity remap
    CALL dm_x_bcs

    DO iy = 1, ny
      DO ix = 1, nx
        ixm = ix - 1
        rho(ix,iy) = (rho1(ix,iy) * cv1(ix,iy) &
            + dm(ixm,iy) - dm(ix,iy)) / cv2(ix,iy)
      END DO
    END DO

    ! Remap specific energy density using mass coordinates
    CALL x_energy_flux

    DO iy = 1, ny
      DO ix = 1, nx
        ixm = ix - 1
        energy(ix,iy) = (energy(ix,iy) * cv1(ix,iy) &
            * rho1(ix,iy) + flux(ixm,iy) - flux(ix,iy)) &
            / (cv2(ix,iy) * rho(ix,iy))
      END DO
    END DO

    ! Redefine dxb1, cv1, cv2, dm and vx1 for velocity (vertex) cells.
    ! In some of these calculations the flux variable is used as a
    ! temporary array
    DO iy = 0, ny
      iyp = iy + 1
      DO ix = -1, nx + 1
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

    ! cv1 = vertex CV before remap
    DO iy = 0, ny
      iyp = iy + 1
      DO ix = 0, nx
        ixp = ix + 1
        flux(ix,iy) = cv1(ix,iy) + cv1(ixp,iy) + cv1(ix,iyp) + cv1(ixp,iyp)
      END DO
    END DO
    cv1(0:nx,0:ny) = flux(0:nx,0:ny) * 0.25_num

    ! cv2 = vertex CV after remap
    DO iy = 0, ny
      iyp = iy + 1
      DO ix = 0, nx
        ixp = ix + 1
        flux(ix,iy) = cv2(ix,iy) + cv2(ixp,iy) + cv2(ix,iyp) + cv2(ixp,iyp)
      END DO
    END DO
    cv2(0:nx,0:ny) = flux(0:nx,0:ny) * 0.25_num

    ! Vertex boundary velocity used in remap
    DO iy = 0, ny
      DO ix = -2, nx + 1
        ixp = ix + 1
        flux(ix,iy) = (vx1(ix,iy) + vx1(ixp,iy)) * 0.5_num
      END DO
    END DO
    vx1(-2:nx+1,0:ny) = flux(-2:nx+1,0:ny)

    ! Mass flux out of vertex CV
    DO iy = 0, ny
      iyp = iy + 1
      DO ix = -1, nx
        ixp = ix + 1
        flux(ix,iy) = dm(ix,iy) + dm(ixp,iy) + dm(ix,iyp) + dm(ixp,iyp)
      END DO
    END DO
    dm(-1:nx,0:ny) = flux(-1:nx,0:ny) * 0.25_num

    ! Vertex density after remap
    DO iy = 0, ny
      DO ix = 0, nx
        ixm = ix - 1
        rho_v1(ix,iy) = (rho_v(ix,iy) * cv1(ix,iy) &
            + dm(ixm,iy) - dm(ix,iy)) / cv2(ix,iy)
      END DO
    END DO

    CALL x_momx_flux

    DO iy = 0, ny
      DO ix = 0, nx
        ixm = ix - 1
        vx(ix,iy) = (rho_v(ix,iy) * vx(ix,iy) * cv1(ix,iy) &
            + flux(ixm,iy) - flux(ix,iy)) / (cv2(ix,iy) * rho_v1(ix,iy))
      END DO
    END DO

    CALL x_momy_flux

    DO iy = 0, ny
      DO ix = 0, nx
        ixm = ix - 1
        vy(ix,iy) = (rho_v(ix,iy) * vy(ix,iy) * cv1(ix,iy) &
            + flux(ixm,iy) - flux(ix,iy)) / (cv2(ix,iy) * rho_v1(ix,iy))
      END DO
    END DO

    CALL x_momz_flux

    DO iy = 0, ny
      DO ix = 0, nx
        ixm = ix - 1
        vz(ix,iy) = (rho_v(ix,iy) * vz(ix,iy) * cv1(ix,iy) &
            + flux(ixm,iy) - flux(ix,iy)) / (cv2(ix,iy) * rho_v1(ix,iy))
      END DO
    END DO

    CALL boundary_conditions
      
    DEALLOCATE (rho1, dm, cv2, flux, dxb1, dxc1, rho_v, rho_v1)
    xpass = 0

  END SUBROUTINE remap_x



  SUBROUTINE vx_by_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: fm, fi, fp, fp2, fu, dfm, dfi, dfp, dfu, sign_v
    REAL(num) :: dxci, dxcu, dxbu, phi, ss, Da, Di
    REAL(num) :: dxu, dbx, dbxp, dbxp2, dbxm
    INTEGER :: ixp2

    DO iy = -1, ny+1
      iyp = iy + 1
      DO ix = -1, nx+1
        ixm  = MAX(ix - 1, -1)
        ixp  = ix + 1
        ixp2 = MIN(ix + 2, nx + 2)

        v_advect = vx1(ix,iy)

        dbx   = (dxb1(ix  ,iy) + dxb1(ix  ,iyp)) * 0.5_num
        dbxp  = (dxb1(ixp ,iy) + dxb1(ixp ,iyp)) * 0.5_num
        dbxp2 = (dxb1(ixp2,iy) + dxb1(ixp2,iyp)) * 0.5_num
        dbxm  = (dxb1(ixm ,iy) + dxb1(ixm ,iyp)) * 0.5_num

        fm  = by(ixm,iy) / dbxm
        fi  = by(ix ,iy) / dbx
        fp  = by(ixp,iy) / dbxp
        fp2 = by(ixp2,iy) / dbxp2

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
        dxci = dxc1(ix ,iy)
        dxcu = dxc1(ixm,iy) * vad_p + dxc1(ixp,iy) * vad_m
        dxbu = dxb1(ix ,iy) * vad_p + dxb1(ixp,iy) * vad_m

        dxu = dbx * vad_p + dbxp * vad_m
        phi = ABS(v_advect) * dt / dxu

        Da =  (2.0_num - phi) * ABS(dfi) / dxci &
            + (1.0_num + phi) * ABS(dfu) / dxcu
        Da = Da * sixth

        ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

        Di = sign_v * ss * MIN(ABS(Da) * dxbu, ABS(dfi), ABS(dfu))

        flux(ix,iy) = (fu + Di * (1.0_num - phi)) * v_advect * dt
      END DO
    END DO

  END SUBROUTINE vx_by_flux



  SUBROUTINE vx_bz_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: fm, fi, fp, fp2, fu, dfm, dfi, dfp, dfu, sign_v
    REAL(num) :: dxci, dxcu, dxbu, phi, ss, Da, Di
    REAL(num) :: dxu, dbx, dbxp, dbxp2, dbxm
    INTEGER :: ixp2

    DO iy = -1, ny + 1
      iym = iy - 1
      DO ix = -1, nx + 1
        ixm  = MAX(ix - 1, -1)
        ixp  = ix + 1
        ixp2 = MIN(ix + 2, nx+2)

        v_advect = (vx1(ix,iy) + vx1(ix,iym)) * 0.5_num

        dbx   = dxb1(ix  ,iy)
        dbxp  = dxb1(ixp ,iy)
        dbxp2 = dxb1(ixp2,iy)
        dbxm  = dxb1(ixm ,iy)

        fm  = bz(ixm,iy) / dbxm
        fi  = bz(ix ,iy) / dbx
        fp  = bz(ixp,iy) / dbxp
        fp2 = bz(ixp2,iy) / dbxp2

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
        dxci = dxc1(ix ,iy)
        dxcu = dxc1(ixm,iy) * vad_p + dxc1(ixp,iy) * vad_m
        dxbu = dxb1(ix ,iy) * vad_p + dxb1(ixp,iy) * vad_m

        dxu = dbx * vad_p + dbxp * vad_m
        phi = ABS(v_advect) * dt / dxu

        Da =  (2.0_num - phi) * ABS(dfi) / dxci &
            + (1.0_num + phi) * ABS(dfu) / dxcu
        Da = Da * sixth

        ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

        Di = sign_v * ss * MIN(ABS(Da) * dxbu, ABS(dfi), ABS(dfu))

        flux(ix,iy) = (fu + Di * (1.0_num - phi)) * v_advect * dt
      END DO
    END DO

  END SUBROUTINE vx_bz_flux



  SUBROUTINE x_mass_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: fm, fi, fp, fp2, fu, dfm, dfi, dfp, dfu, sign_v
    REAL(num) :: dxci, dxcu, dxbu, phi, ss, Da, Di
    REAL(num) :: area
    INTEGER :: ixp2

    DO iy = 0, ny + 1
      iym = iy - 1
      area = dyb(iy)
      DO ix = 0, nx
        ixm  = ix - 1
        ixp  = ix + 1
        ixp2 = ix + 2

        v_advect = (vx1(ix,iy) + vx1(ix,iym)) * 0.5_num

        fm  = rho(ixm ,iy)
        fi  = rho(ix  ,iy)
        fp  = rho(ixp ,iy)
        fp2 = rho(ixp2,iy)

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
        dxci = dxc1(ix ,iy)
        dxcu = dxc1(ixm,iy) * vad_p + dxc1(ixp,iy) * vad_m
        dxbu = dxb1(ix ,iy) * vad_p + dxb1(ixp,iy) * vad_m

        phi = ABS(v_advect) * dt / dxbu

        Da =  (2.0_num - phi) * ABS(dfi) / dxci &
            + (1.0_num + phi) * ABS(dfu) / dxcu
        Da = Da * sixth

        ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

        Di = sign_v * ss * MIN(ABS(Da) * dxbu, ABS(dfi), ABS(dfu))

        dm(ix,iy) = (fu + Di * (1.0_num - phi)) * v_advect * dt * area
      END DO
    END DO

  END SUBROUTINE x_mass_flux



  ! Energy remap in mass coordinates

  SUBROUTINE x_energy_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: fm, fi, fp, fp2, fu, dfm, dfi, dfp, dfu, sign_v
    REAL(num) :: dxci, dxcu, dxbu, phi, ss, Da, Di
    REAL(num) :: area, rhou, dmu
    INTEGER :: ixp2

    DO iy = 0, ny
      iym = iy - 1
      area = dyb(iy)
      DO ix = 0, nx
        ixm  = ix - 1
        ixp  = ix + 1
        ixp2 = ix + 2

        v_advect = (vx1(ix,iy) + vx1(ix,iym)) * 0.5_num

        fm  = energy(ixm ,iy)
        fi  = energy(ix  ,iy)
        fp  = energy(ixp ,iy)
        fp2 = energy(ixp2,iy)

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
        dxci = dxc1(ix ,iy)
        dxcu = dxc1(ixm,iy) * vad_p + dxc1(ixp,iy) * vad_m
        dxbu = dxb1(ix ,iy) * vad_p + dxb1(ixp,iy) * vad_m

        phi = ABS(v_advect) * dt / dxbu

        Da =  (2.0_num - phi) * ABS(dfi) / dxci &
            + (1.0_num + phi) * ABS(dfu) / dxcu
        Da = Da * sixth

        ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

        Di = sign_v * ss * MIN(ABS(Da) * dxbu, ABS(dfi), ABS(dfu))

        rhou = rho1(ix,iy) * vad_p + rho1(ixp,iy) * vad_m
        dmu = ABS(dm(ix,iy)) / area / dxbu / rhou

        flux(ix,iy) = (fu + Di * (1.0_num - dmu)) * dm(ix,iy)
      END DO
    END DO

  END SUBROUTINE x_energy_flux



  SUBROUTINE x_momx_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: fm, fi, fp, fp2, fu, dfm, dfi, dfp, dfu, sign_v
    REAL(num) :: dxci, dxcu, dxbu, phi, ss, Da, Di
    REAL(num) :: area, rhou, dmu, m, mp, ai, aip, dk
    INTEGER :: ixp2

    DO iy = 0, ny
      area = dyc(iy)
      DO ix = -1, nx
        ixm  = ix - 1
        ixp  = ix + 1
        ixp2 = ix + 2

        v_advect = vx1(ix,iy)

        fm  = vx(ixm ,iy)
        fi  = vx(ix  ,iy)
        fp  = vx(ixp ,iy)
        fp2 = vx(ixp2,iy)

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
        dxci = dxb1(ixp,iy)
        dxcu = dxb1(ix ,iy) * vad_p + dxb1(ixp2,iy) * vad_m
        dxbu = dxc1(ix ,iy) * vad_p + dxc1(ixp ,iy) * vad_m

        phi = ABS(v_advect) * dt / dxbu

        Da =  (2.0_num - phi) * ABS(dfi) / dxci &
            + (1.0_num + phi) * ABS(dfu) / dxcu
        Da = Da * sixth

        ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

        Di = sign_v * ss * MIN(ABS(Da) * dxbu, ABS(dfi), ABS(dfu))

        rhou = rho_v(ix,iy) * vad_p + rho_v(ixp,iy) * vad_m
        dmu = ABS(dm(ix,iy)) / area / dxbu / rhou

        flux(ix,iy) = fu + Di * (1.0_num - dmu)
      END DO
    END DO

    IF (rke) THEN
      DO iy = 0, ny
        iyp = iy + 1
        DO ix = 0, nx - 1
          ixm = ix - 1
          ixp = ix + 1

          m  = rho_v1(ix ,iy) * cv2(ix ,iy)
          mp = rho_v1(ixp,iy) * cv2(ixp,iy)

          ai =  (vx(ix ,iy) - flux(ixm,iy)) * dm(ixm,iy) / m &
              - (vx(ix ,iy) - flux(ix ,iy)) * dm(ix ,iy) / m

          aip = (vx(ixp,iy) - flux(ix ,iy)) * dm(ix ,iy) / mp &
              - (vx(ixp,iy) - flux(ixp,iy)) * dm(ixp,iy) / mp

          dk = (vx(ixp,iy) - vx(ix,iy)) * (flux(ix,iy) &
              - 0.5_num * (vx(ixp,iy) + vx(ix,iy))) &
              - 0.5_num * ai  * (vx(ix ,iy) - flux(ix,iy)) &
              + 0.5_num * aip * (vx(ixp,iy) - flux(ix,iy))

          dk = dk * dm(ix,iy) * 0.5_num
          delta_ke(ixp,iy ) = delta_ke(ixp,iy ) + dk
          delta_ke(ixp,iyp) = delta_ke(ixp,iyp) + dk
        END DO
      END DO
    END IF

    flux(-1:nx,0:ny) = flux(-1:nx,0:ny) * dm(-1:nx,0:ny)

  END SUBROUTINE x_momx_flux



  SUBROUTINE x_momy_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: fm, fi, fp, fp2, fu, dfm, dfi, dfp, dfu, sign_v
    REAL(num) :: dxci, dxcu, dxbu, phi, ss, Da, Di
    REAL(num) :: area, rhou, dmu, m, mp, ai, aip, dk
    INTEGER :: ixp2

    DO iy = 0, ny
      area = dyc(iy)
      DO ix = -1, nx
        ixm  = ix - 1
        ixp  = ix + 1
        ixp2 = ix + 2

        v_advect = vx1(ix,iy)

        fm  = vy(ixm ,iy)
        fi  = vy(ix  ,iy)
        fp  = vy(ixp ,iy)
        fp2 = vy(ixp2,iy)

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
        dxci = dxb1(ixp,iy)
        dxcu = dxb1(ix ,iy) * vad_p + dxb1(ixp2,iy) * vad_m
        dxbu = dxc1(ix ,iy) * vad_p + dxc1(ixp ,iy) * vad_m

        phi = ABS(v_advect) * dt / dxbu

        Da =  (2.0_num - phi) * ABS(dfi) / dxci &
            + (1.0_num + phi) * ABS(dfu) / dxcu
        Da = Da * sixth

        ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

        Di = sign_v * ss * MIN(ABS(Da) * dxbu, ABS(dfi), ABS(dfu))

        rhou = rho_v(ix,iy) * vad_p + rho_v(ixp,iy) * vad_m
        dmu = ABS(dm(ix,iy)) / area / dxbu / rhou

        flux(ix,iy) = fu + Di * (1.0_num - dmu)
      END DO
    END DO

    IF (rke) THEN
      DO iy = 0, ny
        iyp = iy + 1
        DO ix = 0, nx - 1
          ixm = ix - 1
          ixp = ix + 1

          m  = rho_v1(ix ,iy) * cv2(ix ,iy)
          mp = rho_v1(ixp,iy) * cv2(ixp,iy)

          ai =  (vy(ix ,iy) - flux(ixm,iy)) * dm(ixm,iy) / m &
              - (vy(ix ,iy) - flux(ix ,iy)) * dm(ix ,iy) / m

          aip = (vy(ixp,iy) - flux(ix ,iy)) * dm(ix ,iy) / mp &
              - (vy(ixp,iy) - flux(ixp,iy)) * dm(ixp,iy) / mp

          dk = (vy(ixp,iy) - vy(ix,iy)) * (flux(ix,iy) &
              - 0.5_num * (vy(ixp,iy) + vy(ix,iy))) &
              - 0.5_num * ai  * (vy(ix ,iy) - flux(ix,iy)) &
              + 0.5_num * aip * (vy(ixp,iy) - flux(ix,iy))

          dk = dk * dm(ix,iy) * 0.5_num
          delta_ke(ixp,iy ) = delta_ke(ixp,iy ) + dk
          delta_ke(ixp,iyp) = delta_ke(ixp,iyp) + dk
        END DO
      END DO
    END IF

    flux(-1:nx,0:ny) = flux(-1:nx,0:ny) * dm(-1:nx,0:ny)

  END SUBROUTINE x_momy_flux



  SUBROUTINE x_momz_flux

    REAL(num) :: v_advect, vad_p, vad_m
    REAL(num) :: fm, fi, fp, fp2, fu, dfm, dfi, dfp, dfu, sign_v
    REAL(num) :: dxci, dxcu, dxbu, phi, ss, Da, Di
    REAL(num) :: area, rhou, dmu, m, mp, ai, aip, dk
    INTEGER :: ixp2

    DO iy = 0, ny
      area = dyc(iy)
      DO ix = -1, nx
        ixm  = ix - 1
        ixp  = ix + 1
        ixp2 = ix + 2

        v_advect = vx1(ix,iy)

        fm  = vz(ixm ,iy)
        fi  = vz(ix  ,iy)
        fp  = vz(ixp ,iy)
        fp2 = vz(ixp2,iy)

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
        dxci = dxb1(ixp,iy)
        dxcu = dxb1(ix ,iy) * vad_p + dxb1(ixp2,iy) * vad_m
        dxbu = dxc1(ix ,iy) * vad_p + dxc1(ixp ,iy) * vad_m

        phi = ABS(v_advect) * dt / dxbu

        Da =  (2.0_num - phi) * ABS(dfi) / dxci &
            + (1.0_num + phi) * ABS(dfu) / dxcu
        Da = Da * sixth

        ss = 0.5_num * (SIGN(1.0_num, dfi) + SIGN(1.0_num, dfu))

        Di = sign_v * ss * MIN(ABS(Da) * dxbu, ABS(dfi), ABS(dfu))

        rhou = rho_v(ix,iy) * vad_p + rho_v(ixp,iy) * vad_m
        dmu = ABS(dm(ix,iy)) / area / dxbu / rhou

        flux(ix,iy) = fu + Di * (1.0_num - dmu)
      END DO
    END DO

    IF (rke) THEN
      DO iy = 0, ny
        iyp = iy + 1
        DO ix = 0, nx - 1
          ixm = ix - 1
          ixp = ix + 1

          m  = rho_v1(ix ,iy) * cv2(ix ,iy)
          mp = rho_v1(ixp,iy) * cv2(ixp,iy)

          ai =  (vz(ix ,iy) - flux(ixm,iy)) * dm(ixm,iy) / m &
              - (vz(ix ,iy) - flux(ix ,iy)) * dm(ix ,iy) / m

          aip = (vz(ixp,iy) - flux(ix ,iy)) * dm(ix ,iy) / mp &
              - (vz(ixp,iy) - flux(ixp,iy)) * dm(ixp,iy) / mp

          dk = (vz(ixp,iy) - vz(ix,iy)) * (flux(ix,iy) &
              - 0.5_num * (vz(ixp,iy) + vz(ix,iy))) &
              - 0.5_num * ai  * (vz(ix ,iy) - flux(ix,iy)) &
              + 0.5_num * aip * (vz(ixp,iy) - flux(ix,iy))

          dk = dk * dm(ix,iy) * 0.5_num
          delta_ke(ixp,iy ) = delta_ke(ixp,iy ) + dk
          delta_ke(ixp,iyp) = delta_ke(ixp,iyp) + dk
        END DO
      END DO
    END IF

    flux(-1:nx,0:ny) = flux(-1:nx,0:ny) * dm(-1:nx,0:ny)

  END SUBROUTINE x_momz_flux



  SUBROUTINE dm_x_bcs

    CALL MPI_SENDRECV(&
        dm(1   ,-1), 1, cell_xface, proc_x_min, tag, &
        dm(nx+1,-1), 1, cell_xface, proc_x_max, tag, &
        comm, status, errcode)

    IF (proc_x_max == MPI_PROC_NULL) &
        dm(nx+1,0:ny+1) = dm(nx,0:ny+1)

    CALL MPI_SENDRECV(&
        dm(nx-1,-1), 1, cell_xface, proc_x_max, tag, &
        dm(-1  ,-1), 1, cell_xface, proc_x_min, tag, &
        comm, status, errcode)

    IF (proc_x_min == MPI_PROC_NULL) &
        dm(-1,0:ny+1) = dm(0,0:ny+1)

  END SUBROUTINE dm_x_bcs

END MODULE xremap
