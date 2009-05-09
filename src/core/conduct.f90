MODULE conduct

  USE shared_data
  USE boundary
  USE normalise
  USE eos

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Conduct_Heat

CONTAINS

  !  Subroutine implementing Braginskii parallel thermal conduction
  ! Note that this subroutine assumes that it is possible to
  ! Convert from specific internal energy to temperature by a constant
  ! factor.

  ! For this as well as other reasons, this routine doesn't work properly
  ! for partially ionised plasmas.
  SUBROUTINE Conduct_Heat

    REAL(num), DIMENSION(:, :), ALLOCATABLE :: kx, ky, ux, uy
    REAL(num), DIMENSION(:, :), ALLOCATABLE :: energy0, EnergyToT
    REAL(num) :: e, T
    REAL(num) :: B, bxc, byc, bzc
    REAL(num) :: pow = 5.0_num / 2.0_num
    REAL(num) :: QpX, QmX, Q0X
    REAL(num) :: QpY, QmY, Q0Y

    REAL(num) :: MpX, MmX, M0X
    REAL(num) :: MpY, MmY, M0Y

    REAL(num) :: Rx, Ry
    REAL(num) :: Rxx, Ryy
    REAL(num) :: Rxy

    REAL(num) :: KxX, KyX
    REAL(num) :: KxY, KyY

    REAL(num) :: uxX, uyY

    REAL(num) :: A1, A2, mx, Q, errtot, mx1, errtot_prev
    REAL(num), SAVE :: w = 1.9_num

    REAL(dbl) :: runtime = 0.0_dbl, starttime, endtime

    INTEGER :: CYCLE, sweep, mx_x, mx_y, mx_z, start_index

    LOGICAL :: converged

    REAL(num), PARAMETER :: b_min = 1.0e-4_num

    w = (w + 2.0_num) / 2.0_num
    mx = 0.0_num

    ALLOCATE(kx(0:nx+1, 0:ny+1), ky(0:nx+1, 0:ny+1))
    ALLOCATE(ux(0:nx+1, 0:ny+1), uy(0:nx+1, 0:ny+1))
    ALLOCATE(energy0(-1:nx+1, -1:ny+1), EnergyToT(-1:nx+1, -1:ny+1))

    energy0 = energy(-1:nx+1, -1:ny+1)
    converged = .FALSE.
    DO iy = -1, ny + 1
      DO ix = -1, nx + 1
        e = energy(ix, iy)
        CALL Get_Temp(rho(ix, iy), e, eos_number, ix, iy, T)
        EnergyToT(ix, iy) = T / e
      END DO
    END DO

    DO iy = 0, ny + 1
      DO ix = 0, nx + 1
        bxc = (bx(ix, iy) * dxb(ix) + bx(ix-1, iy) * dxb(ix-1)) &
            / (dxb(ix) + dxb(ix-1))
        byc = (by(ix, iy) * dyb(iy) + by(ix, iy-1) * dyb(iy-1)) &
            / (dyb(iy) + dyb(iy-1))
        B = SQRT(bxc**2 + byc**2 + bz(ix, iy)**2)

        ! Direction of magnetic field
        ux(ix, iy) = bxc**2 / (B**2 + b_min**2)
        uy(ix, iy) = byc**2 / (B**2 + b_min**2)

        ! Kappa along magnetic field
        kx(ix, iy) = kappa * (EnergyToT(ix, iy) * energy0(ix, iy))**pow &
            * (ABS(ux(ix, iy)) + b_min**2 / (B**2 + b_min**2))
        ky(ix, iy) = kappa * (EnergyToT(ix, iy) * energy0(ix, iy))**pow &
            * (ABS(uy(ix, iy)) + b_min**2 / (B**2 + b_min**2))
      END DO
    END DO

    DO CYCLE = 0, 500
      errtot = 0.0_num
      mx = 0.0_num
      mx1 = 0.0_num
      DO iy = 1, ny
        start_index = 1
        DO ix = start_index, nx! , 2
          QpX = dxc(ix-1) / (dxc(ix) * (dxc(ix) + dxc(ix-1)))
          QmX = dxc(ix) / (dxc(ix-1) * (dxc(ix) + dxc(ix-1)))
          Q0X = (dxc(ix)**2 - dxc(ix-1)**2) &
              / (dxc(ix) * dxc(ix-1) * (dxc(ix) + dxc(ix-1)))

          QpY = dyc(iy-1) / (dyc(iy) * (dyc(iy) + dyc(iy-1)))
          QmY = dyc(iy) / (dyc(iy-1) * (dyc(iy) + dyc(iy-1)))
          Q0Y = (dyc(iy)**2 - dyc(iy-1)**2) &
              / (dyc(iy) * dyc(iy-1) * (dyc(iy) + dyc(iy-1)))

          MpX = 1.0_num / (dxc(ix) * dxb(ix))
          MmX = 1.0_num / (dxc(ix-1) * dxb(ix))
          M0X = (dxc(ix) + dxc(ix-1)) / (dxc(ix) * dxc(ix-1) * dxb(ix))

          MpY = 1.0_num / (dyc(iy) * dyb(iy))
          MmY = 1.0_num / (dyc(iy-1) * dyb(iy))
          M0Y = (dyc(iy) + dyc(iy-1)) / (dyc(iy) * dyc(iy-1) * dyb(iy))

          Rx = QpX * EnergyToT(ix+1, iy) * energy(ix+1, iy) &
              - QmX * EnergyToT(ix-1, iy) * energy(ix-1, iy)
          Ry = QpY * EnergyToT(ix, iy+1) * energy(ix, iy+1) &
              - QmY * EnergyToT(ix, iy-1) * energy(ix, iy-1)

          Rxx = MpX * EnergyToT(ix+1, iy) * energy(ix+1, iy) &
              + MmX * EnergyToT(ix-1, iy) * energy(ix-1, iy)
          Ryy = MpY * EnergyToT(ix, iy+1) * energy(ix, iy+1) &
              + MmY * EnergyToT(ix, iy-1) * energy(ix, iy-1)

          Rxy = QpY / 16.0_num &
              * (QpX * EnergyToT(ix+1, iy+1) * energy(ix+1, iy+1) &
              - QmX * EnergyToT(ix-1, iy+1) * energy(ix-1, iy+1) &
              + Q0X * EnergyToT(ix, iy+1) * energy(ix, iy+1)) &
              - QmY / 16.0_num &
              * (QpX * EnergyToT(ix+1, iy-1) * energy(ix+1, iy-1) &
              - QmX * EnergyToT(ix-1, iy-1) * energy(ix-1, iy-1) &
              + Q0X * EnergyToT(ix , iy-1) * energy(ix, iy-1)) &
              + Q0Y / 16.0_num &
              * (QpX * EnergyToT(ix+1, iy) * energy(ix+1, iy) &
              - QmX * EnergyToT(ix-1, iy) * energy(ix-1, iy))

          KxX = QpX * kx(ix+1, iy) - QmX * kx(ix-1, iy) + Q0X * kx(ix, iy)
          KyX = QpX * ky(ix+1, iy) - QmX * ky(ix-1, iy) + Q0X * ky(ix, iy)

          KxY = QpY * kx(ix, iy+1) - QmY * kx(ix, iy-1) + Q0Y * kx(ix, iy)
          KyY = QpY * ky(ix, iy+1) - QmY * ky(ix, iy-1) + Q0Y * ky(ix, iy)

          uxX = QpX * ux(ix+1, iy) - QmX * ux(ix-1, iy) + Q0X * ux(ix, iy)
          uyY = QpY * uy(ix, iy+1) - QmY * uy(ix, iy-1) + Q0Y * uy(ix, iy)

          bxc = (bx(ix, iy) + bx(ix-1, iy)) / 2.0_num
          byc = (by(ix, iy) + by(ix, iy-1)) / 2.0_num
          bzc = bz(ix, iy)

          B = SQRT(bxc**2 + byc**2 + bzc**2)
          ! Second differentials in T
          A1 = M0X * ux(ix, iy) * kx(ix, iy) &
              + M0Y * uy(ix, iy) * ky(ix, iy) &
              + Q0X * Q0Y / 16.0_num &
              * (ux(ix, iy) * ky(ix, iy) + uy(ix, iy) * kx(ix, iy))
          ! Differentials in kx, ky
          A1 = A1 + ux(ix, iy) * (KxX * Q0X + KyX * Q0Y) &
              + uy(ix, iy) * (KxY * Q0X + KyY * Q0Y)
          ! Differentials in ux, uy
          A1 = A1 + Q0X * kx(ix, iy) * (uxX + uyY) &
              + Q0Y * ky(ix, iy) * (uxX + uyY)

          ! Second differentials in T
          A2 = Rxx * ux(ix, iy) * kx(ix, iy) + Ryy * uy(ix, iy) * ky(ix, iy) &
              + Rxy * (ux(ix, iy) * ky(ix, iy) + uy(ix, iy) * kx(ix, iy))
          ! Differentials in kx, ky
          A2 = A2 + ux(ix, iy) * (KxX * Rx + KyX * Ry) &
              + uy(ix, iy) * (KxY * Rx + KyY * Ry)
          ! Differentials in ux, uy
          A2 = A2 + uxX * (kx(ix, iy) * Rx + ky(ix, iy) * Ry) &
              + uyY * (kx(ix, iy) * Rx + ky(ix, iy) * Ry)

          A1 = A1 + (kx(ix, iy) * M0X + ky(ix, iy) * M0Y &
              + KxX * Q0X + KyY * Q0Y) * b_min**2 / (b**2 + b_min**2)
          A2 = A2 + (kx(ix, iy) * Rxx + ky(ix, iy) * Ryy &
              + KxX * Rx + KyY * Ry) * b_min**2 / (b**2 + b_min**2)
          A1 = A1 * dt / rho(ix, iy)
          A2 = A2 * dt / rho(ix, iy)

          Q = energy(ix, iy)
          energy(ix, iy) = (1.0_num-w) * energy(ix, iy) &
              + w / (1.0_num / EnergyToT(ix, iy) + A1) &
              * (energy0(ix, iy) / EnergyToT(ix, iy) + A2 / EnergyToT(ix, iy))
          Q = (Q - energy(ix, iy)) / Q

          errtot = errtot + ABS(Q) * dxc(ix) * dyc(iy)
        END DO
      END DO
      CALL energy_bcs

      CALL MPI_ALLREDUCE(errtot, mx, 1, mpireal, MPI_SUM, comm, errcode)
      errtot = mx

      IF (errtot .GT. errtot_prev) w = (1.0_num+w) / 2.0_num
      errtot_prev = errtot

      IF (errtot .LT. 1e-4_num) THEN
        converged = .TRUE.
        EXIT
      END IF
    END DO

    IF (.NOT. converged) &
        PRINT * , "Solution failed to converge during heat conduction", errtot

    ! PRINT * , "Cycle Over", w
    DEALLOCATE(kx, ky)
    DEALLOCATE(ux, uy)
    DEALLOCATE(energy0, EnergyToT)

  END SUBROUTINE Conduct_Heat

END MODULE conduct
