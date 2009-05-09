MODULE conduct

  USE shared_data
  USE boundary
  USE eos

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: conduct_heat

CONTAINS

  ! Subroutine implementing Braginskii parallel thermal conduction
  ! Note that this subroutine assumes that it is possible to
  ! Convert from specific internal energy to temperature by a constant
  ! factor.
  ! For this as well as other reasons, this routine doesn't work properly
  ! for partially ionised plasmas. 
	! Notation and algorithm in Appendix of Manual
  SUBROUTINE conduct_heat

    REAL(num), DIMENSION(:, :), ALLOCATABLE :: kx, ky, ux, uy
    REAL(num), DIMENSION(:, :), ALLOCATABLE :: e2temp, kp, energy0
    REAL(num) :: e, T
    REAL(num) :: b, bxc, byc, bzc
    REAL(num) :: pow = 5.0_num / 2.0_num
    REAL(num) :: qpx, qmx, q0x
    REAL(num) :: qpy, qmy, q0y
    REAL(num) :: mpx, mmx, m0x
    REAL(num) :: mpy, mmy, m0y
    REAL(num) :: rx, ry
    REAL(num) :: rxx, ryy
    REAL(num) :: rxy
    REAL(num) :: kxx, kyx, kpx
    REAL(num) :: kxy, kyy, kpy
    REAL(num) :: uxx, uyy
    REAL(num) :: a1, a2, error, q, errmax, errmax_prev = 0.0_num
    REAL(num) :: w 

    INTEGER :: loop

    LOGICAL :: converged
    REAL(num), PARAMETER :: b_min = 1.0e-4_num

    ALLOCATE(kx(0:nx+1, 0:ny+1), ky(0:nx+1, 0:ny+1))
    ALLOCATE(ux(0:nx+1, 0:ny+1), uy(0:nx+1, 0:ny+1))
    ALLOCATE(e2temp(-1:nx+1, -1:ny+1), kp(-1:nx+1, -1:ny+1), energy0(-1:nx+2, -1:ny+2))
            
kappa_0 = 1.0_num
energy = energy + 4.0_num * kappa_0 * dt / 7.0_num
CALL energy_bcs

		! find factor reuired to convert between energy and temperature
    DO iy = -1, ny + 1
      DO ix = -1, nx + 1
        e = energy(ix, iy)
        CALL get_temp(rho(ix, iy), e, eos_number, ix, iy, T)
        e2temp(ix, iy) = T / e
      END DO
    END DO

    DO iy = 0, ny + 1
      DO ix = 0, nx + 1
        bxc = (bx(ix, iy) + bx(ix-1, iy)) / 2.0_num
        byc = (by(ix, iy) + by(ix, iy-1)) / 2.0_num
        b = SQRT(bxc**2 + byc**2 + bz(ix, iy)**2)

        ! Direction of magnetic field
        ux(ix, iy) = bxc / SQRT(B**2 + b_min**2)
        uy(ix, iy) = byc / SQRT(B**2 + b_min**2)

        ! Kappa along magnetic field, now a vector 
				T = (e2temp(ix, iy) * energy(ix, iy))**pow
        kx(ix, iy) = kappa_0 * T * ux(ix, iy) 
        ky(ix, iy) = kappa_0 * T * uy(ix, iy) 

				! Kappa isotropic
				kp(ix,iy) = kappa_0 * T * b_min**2 / (B**2 + b_min**2)
      END DO
    END DO
     
    converged = .FALSE. 
    w = 1.9_num       ! initial over-relaxation parameter
		! store energy^{n} 
		energy0 = energy   
		! interate to get energy^{n+1} by SOR Guass-Seidel
    DO loop = 0, 100
      errmax = 0.0_num
      error = 0.0_num
      DO iy = 1, ny
        DO ix = 1, nx
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

          rx = qpx * e2temp(ix+1, iy) * energy(ix+1, iy) &
              - qmx * e2temp(ix-1, iy) * energy(ix-1, iy)
          ry = qpy * e2temp(ix, iy+1) * energy(ix, iy+1) &
              - qmy * e2temp(ix, iy-1) * energy(ix, iy-1)

          rxx = mpx * e2temp(ix+1, iy) * energy(ix+1, iy) &
              + mmx * e2temp(ix-1, iy) * energy(ix-1, iy)
          ryy = mpy * e2temp(ix, iy+1) * energy(ix, iy+1) &
              + mmy * e2temp(ix, iy-1) * energy(ix, iy-1)

          rxy = qpy * &
									(qpx * e2temp(ix+1, iy+1) * energy(ix+1, iy+1) &
		              - qmx * e2temp(ix-1, iy+1) * energy(ix-1, iy+1) &
		              + q0x * e2temp(ix, iy+1) * energy(ix, iy+1)) &
              - qmy * &
									(qpx * e2temp(ix+1, iy-1) * energy(ix+1, iy-1) &
		              - qmx * e2temp(ix-1, iy-1) * energy(ix-1, iy-1) &
		              + q0x * e2temp(ix , iy-1) * energy(ix, iy-1)) &
              + q0y * &
									(qpx * e2temp(ix+1, iy) * energy(ix+1, iy) &
		              - qmx * e2temp(ix-1, iy) * energy(ix-1, iy))

          kxx = qpx * kx(ix+1, iy) - qmx * kx(ix-1, iy) + q0x * kx(ix, iy)
          kyx = qpx * ky(ix+1, iy) - qmx * ky(ix-1, iy) + q0x * ky(ix, iy)
          kpx = qpx * kp(ix+1, iy) - qmx * kp(ix-1, iy) + q0x * kx(ix, iy)

          kxy = qpy * kx(ix, iy+1) - qmy * kx(ix, iy-1) + q0y * kx(ix, iy)
          kyy = qpy * ky(ix, iy+1) - qmy * ky(ix, iy-1) + q0y * ky(ix, iy)
          kpy = qpy * kp(ix, iy+1) - qmy * kp(ix, iy-1) + q0y * kp(ix, iy)

          uxx = qpx * ux(ix+1, iy) - qmx * ux(ix-1, iy) + q0x * ux(ix, iy)
          uyy = qpy * uy(ix, iy+1) - qmy * uy(ix, iy-1) + q0y * uy(ix, iy)

          bxc = (bx(ix, iy) + bx(ix-1, iy)) / 2.0_num
          byc = (by(ix, iy) + by(ix, iy-1)) / 2.0_num
          bzc = bz(ix, iy)

          B = SQRT(bxc**2 + byc**2 + bzc**2)
          ! Second differentials in T
          a1 = m0x * ux(ix, iy) * kx(ix, iy) &
              + m0y * uy(ix, iy) * ky(ix, iy) &
              + q0x * q0y * (ux(ix, iy) * ky(ix, iy) + uy(ix, iy) * kx(ix, iy))
          ! Differentials in kx, ky
          a1 = a1 - ux(ix, iy) * (kxx * q0x + kyx * q0y) &
              - uy(ix, iy) * (kxy * q0x + kyy * q0y)
          ! Differentials in ux, uy
          a1 = a1 - q0x * kx(ix, iy) * (uxx + uyy) &
              - q0y * ky(ix, iy) * (uxx + uyy)

          ! Second differentials in T
          a2 = rxx * ux(ix, iy) * kx(ix, iy) + ryy * uy(ix, iy) * ky(ix, iy) &
              + rxy * (ux(ix, iy) * ky(ix, iy) + uy(ix, iy) * kx(ix, iy))
          ! Differentials in kx, ky
          a2 = a2 + ux(ix, iy) * (kxx * rx + kyx * ry) &
              + uy(ix, iy) * (kxy * rx + kyy * ry)
          ! Differentials in ux, uy
          a2 = a2 + uxx * (kx(ix, iy) * rx + ky(ix, iy) * ry) &
              + uyy * (kx(ix, iy) * rx + ky(ix, iy) * ry)
          ! add isotropic elements
          a1 = a1 + kp(ix, iy) * (m0x + m0y) &
                    - kpx * q0x - kpy * q0y 
          a2 = a2 + kp(ix, iy) * (rxx + ryy) &
                    + kpx * rx + kpy * ry  

          a1 = a1 * dt * e2temp(ix, iy) / rho(ix, iy) 
          a2 = a2 * dt / rho(ix, iy) 

          Q = energy(ix, iy)
          energy(ix, iy) = (1.0_num-w) * energy(ix, iy) &
              + w * (energy0(ix, iy)  + a2) / (1.0_num + a1)
          Q = (Q - energy(ix, iy)) / Q

          errmax = MAX(errmax, ABS(Q))
        END DO
      END DO 

      CALL energy_bcs
  
      CALL MPI_ALLREDUCE(errmax, error, 1, mpireal, MPI_MAX, comm, errcode)
      errmax = error

      IF (errmax .GT. errmax_prev) w = (1.0_num + w) / 2.0_num
      errmax_prev = errmax

      IF (errmax .LT. 1.e-5_num) THEN
        converged = .TRUE.  
        EXIT
      END IF
    END DO

    IF (.NOT. converged) &
        PRINT * , "Solution failed to converge during heat conduction", errmax

    DEALLOCATE(kx, ky)
    DEALLOCATE(ux, uy)
    DEALLOCATE(kp, e2temp, energy0)

  END SUBROUTINE conduct_heat

END MODULE conduct
