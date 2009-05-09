MODULE conduct

  USE shared_data
  USE boundary
  USE normalise
  USE eos

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: conduct_heat

CONTAINS

  !  Subroutine implementing Braginskii parallel thermal conduction
  ! Note that this subroutine assumes that it is possible to
  ! Convert from specific internal energy to temperature by a constant
  ! factor.

  ! For this as well as other reasons, this routine doesn't work properly
  ! for partially ionised plasmas.
  SUBROUTINE conduct_heat

    REAL(num), DIMENSION(:, :), ALLOCATABLE :: kx, ky, ux, uy
    REAL(num), DIMENSION(:, :), ALLOCATABLE :: energy0, energytot
    REAL(num) :: e, T
    REAL(num) :: B, bxc, byc, bzc
    REAL(num) :: pow = 5.0_num / 2.0_num
    REAL(num) :: qpx, qmx, q0x
    REAL(num) :: qpy, qmy, q0y

    REAL(num) :: mpx, mmx, m0x
    REAL(num) :: mpy, mmy, m0y

    REAL(num) :: rx, ry
    REAL(num) :: rxx, ryy
    REAL(num) :: rxy

    REAL(num) :: kxx, kyx
    REAL(num) :: kxy, kyy

    REAL(num) :: uxx, uyy

    REAL(num) :: a1, a2, mx, Q, errtot, mx1, errtot_prev = 0.0_num
    REAL(num), SAVE :: w = 1.9_num

    INTEGER :: CYCLE, start_index

    LOGICAL :: converged

    REAL(num), PARAMETER :: b_min = 1.0e-4_num

    w = (w + 2.0_num) / 2.0_num
    mx = 0.0_num

    ALLOCATE(kx(0:nx+1, 0:ny+1), ky(0:nx+1, 0:ny+1))
    ALLOCATE(ux(0:nx+1, 0:ny+1), uy(0:nx+1, 0:ny+1))
    ALLOCATE(energy0(-1:nx+1, -1:ny+1), energytot(-1:nx+1, -1:ny+1))

    energy0 = energy(-1:nx+1, -1:ny+1)
    converged = .FALSE.
    DO iy = -1, ny + 1
      DO ix = -1, nx + 1
        e = energy(ix, iy)
        CALL get_temp(rho(ix, iy), e, eos_number, ix, iy, T)
        energytot(ix, iy) = T / e
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
        kx(ix, iy) = kappa * (energytot(ix, iy) * energy0(ix, iy))**pow &
            * (ABS(ux(ix, iy)) + b_min**2 / (B**2 + b_min**2))
        ky(ix, iy) = kappa * (energytot(ix, iy) * energy0(ix, iy))**pow &
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
          qpx = dxc(ix-1) / (dxc(ix) * (dxc(ix) + dxc(ix-1)))
          qmx = dxc(ix) / (dxc(ix-1) * (dxc(ix) + dxc(ix-1)))
          q0x = (dxc(ix)**2 - dxc(ix-1)**2) &
              / (dxc(ix) * dxc(ix-1) * (dxc(ix) + dxc(ix-1)))

          qpy = dyc(iy-1) / (dyc(iy) * (dyc(iy) + dyc(iy-1)))
          qmy = dyc(iy) / (dyc(iy-1) * (dyc(iy) + dyc(iy-1)))
          q0y = (dyc(iy)**2 - dyc(iy-1)**2) &
              / (dyc(iy) * dyc(iy-1) * (dyc(iy) + dyc(iy-1)))

          mpx = 1.0_num / (dxc(ix) * dxb(ix))
          mmx = 1.0_num / (dxc(ix-1) * dxb(ix))
          m0x = (dxc(ix) + dxc(ix-1)) / (dxc(ix) * dxc(ix-1) * dxb(ix))

          mpy = 1.0_num / (dyc(iy) * dyb(iy))
          mmy = 1.0_num / (dyc(iy-1) * dyb(iy))
          m0y = (dyc(iy) + dyc(iy-1)) / (dyc(iy) * dyc(iy-1) * dyb(iy))

          rx = qpx * energytot(ix+1, iy) * energy(ix+1, iy) &
              - qmx * energytot(ix-1, iy) * energy(ix-1, iy)
          ry = qpy * energytot(ix, iy+1) * energy(ix, iy+1) &
              - qmy * energytot(ix, iy-1) * energy(ix, iy-1)

          rxx = mpx * energytot(ix+1, iy) * energy(ix+1, iy) &
              + mmx * energytot(ix-1, iy) * energy(ix-1, iy)
          ryy = mpy * energytot(ix, iy+1) * energy(ix, iy+1) &
              + mmy * energytot(ix, iy-1) * energy(ix, iy-1)

          rxy = qpy / 16.0_num &
              * (qpx * energytot(ix+1, iy+1) * energy(ix+1, iy+1) &
              - qmx * energytot(ix-1, iy+1) * energy(ix-1, iy+1) &
              + q0x * energytot(ix, iy+1) * energy(ix, iy+1)) &
              - qmy / 16.0_num &
              * (qpx * energytot(ix+1, iy-1) * energy(ix+1, iy-1) &
              - qmx * energytot(ix-1, iy-1) * energy(ix-1, iy-1) &
              + q0x * energytot(ix , iy-1) * energy(ix, iy-1)) &
              + q0y / 16.0_num &
              * (qpx * energytot(ix+1, iy) * energy(ix+1, iy) &
              - qmx * energytot(ix-1, iy) * energy(ix-1, iy))

          kxx = qpx * kx(ix+1, iy) - qmx * kx(ix-1, iy) + q0x * kx(ix, iy)
          kyx = qpx * ky(ix+1, iy) - qmx * ky(ix-1, iy) + q0x * ky(ix, iy)

          kxy = qpy * kx(ix, iy+1) - qmy * kx(ix, iy-1) + q0y * kx(ix, iy)
          kyy = qpy * ky(ix, iy+1) - qmy * ky(ix, iy-1) + q0y * ky(ix, iy)

          uxx = qpx * ux(ix+1, iy) - qmx * ux(ix-1, iy) + q0x * ux(ix, iy)
          uyy = qpy * uy(ix, iy+1) - qmy * uy(ix, iy-1) + q0y * uy(ix, iy)

          bxc = (bx(ix, iy) + bx(ix-1, iy)) / 2.0_num
          byc = (by(ix, iy) + by(ix, iy-1)) / 2.0_num
          bzc = bz(ix, iy)

          B = SQRT(bxc**2 + byc**2 + bzc**2)
          ! Second differentials in T
          a1 = m0x * ux(ix, iy) * kx(ix, iy) &
              + m0y * uy(ix, iy) * ky(ix, iy) &
              + q0x * q0y / 16.0_num &
              * (ux(ix, iy) * ky(ix, iy) + uy(ix, iy) * kx(ix, iy))
          ! Differentials in kx, ky
          a1 = a1 + ux(ix, iy) * (kxx * q0x + kyx * q0y) &
              + uy(ix, iy) * (kxy * q0x + kyy * q0y)
          ! Differentials in ux, uy
          a1 = a1 + q0x * kx(ix, iy) * (uxx + uyy) &
              + q0y * ky(ix, iy) * (uxx + uyy)

          ! Second differentials in T
          a2 = rxx * ux(ix, iy) * kx(ix, iy) + ryy * uy(ix, iy) * ky(ix, iy) &
              + rxy * (ux(ix, iy) * ky(ix, iy) + uy(ix, iy) * kx(ix, iy))
          ! Differentials in kx, ky
          a2 = a2 + ux(ix, iy) * (kxx * rx + kyx * ry) &
              + uy(ix, iy) * (kxy * rx + kyy * ry)
          ! Differentials in ux, uy
          a2 = a2 + uxx * (kx(ix, iy) * rx + ky(ix, iy) * ry) &
              + uyy * (kx(ix, iy) * rx + ky(ix, iy) * ry)

          a1 = a1 + (kx(ix, iy) * m0x + ky(ix, iy) * m0y &
              + kxx * q0x + kyy * q0y) * b_min**2 / (b**2 + b_min**2)
          a2 = a2 + (kx(ix, iy) * rxx + ky(ix, iy) * ryy &
              + kxx * rx + kyy * ry) * b_min**2 / (b**2 + b_min**2)
          a1 = a1 * dt / rho(ix, iy)
          a2 = a2 * dt / rho(ix, iy)

          Q = energy(ix, iy)
          energy(ix, iy) = (1.0_num-w) * energy(ix, iy) &
              + w / (1.0_num / energytot(ix, iy) + a1) &
              * (energy0(ix, iy) / energytot(ix, iy) + a2 / energytot(ix, iy))
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
    DEALLOCATE(energy0, energytot)

  END SUBROUTINE conduct_heat

END MODULE conduct
