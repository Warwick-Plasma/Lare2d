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

<<<<<<< HEAD:src/core/conduct.f90
    REAL(num), DIMENSION(:, :), ALLOCATABLE :: kx, ky, ux, uy
    REAL(num), DIMENSION(:, :), ALLOCATABLE :: e2temp, kp, energy0, energy1
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
    REAL(num) :: a1, a2, error, errmax, errmax_prev = 0.0_num
    REAL(num) :: w, residual
=======
    REAL(num), DIMENSION(:, :), ALLOCATABLE :: uxkx, uxky, uykx, uyky 
    REAL(num), DIMENSION(:, :), ALLOCATABLE :: energy0, limiter
    REAL(num) :: e2t, exb, eyb 
    REAL(num) :: b, bxc, byc, bzc, bpx, bpy 
    REAL(num) :: ux, uy
    REAL(num) :: pow = 5.0_num / 2.0_num  
    REAL(num) :: a1, a2, a3, error, errmax, errmax_prev = 0.0_num
    REAL(num) :: w, residual, q_shx, q_shy, q_sh, q_f, q_nl
    REAL(num) :: initial_energy, final_energy
>>>>>>> conduction:src/core/conduct.f90

    INTEGER :: loop, redblack, x1, y1

    LOGICAL :: converged
<<<<<<< HEAD:src/core/conduct.f90
    REAL(num), PARAMETER :: fractional_error = 1.0e-3_num
    REAL(num), PARAMETER :: b_min = 1.0e-4_num

    ALLOCATE(kx(-1:nx+2, -1:ny+2), ky(-1:nx+2, -1:ny+2))
    ALLOCATE(ux(-1:nx+2, -1:ny+2), uy(-1:nx+2, -1:ny+2))
    ALLOCATE(e2temp(-1:nx+2, -1:ny+2), kp(-1:nx+2, -1:ny+2), energy0(-1:nx+2, -1:ny+2))
    ALLOCATE(energy1(-1:nx+2, -1:ny+2))
            
		! find factor reuired to convert between energy and temperature
    DO iy = -1, ny + 2
      DO ix = -1, nx + 2
        e = energy(ix, iy)
        CALL get_temp(rho(ix, iy), e, eos_number, ix, iy, T)
        e2temp(ix, iy) = T / MAX(e, none_zero)
=======
    REAL(num), PARAMETER :: fractional_error = 1.0e-2_num
    REAL(num), PARAMETER :: b_min = 1.0e-3_num

    ALLOCATE(uxkx(-1:nx+1, 0:ny+1), uxky(-1:nx+1, 0:ny+1))
    ALLOCATE(uykx(-1:nx+1, 0:ny+1), uyky(-1:nx+1, 0:ny+1))
    ALLOCATE(energy0(-1:nx+2, -1:ny+2))
    ALLOCATE(limiter(-1:nx+2, -1:ny+2))
            
    a1 = 0.0_num
    DO iy = 1, ny 
      DO ix = 1, nx
        a1 = a1 + energy(ix,iy) * dxb(ix) * dyb(iy)
>>>>>>> conduction:src/core/conduct.f90
      END DO
    END DO
    CALL MPI_ALLREDUCE(a1, initial_energy, 1, mpireal, MPI_SUM, comm, errcode) 

<<<<<<< HEAD:src/core/conduct.f90
    DO iy = -1, ny + 2
      DO ix = -1, nx + 2
        bxc = (bx(ix, iy) + bx(ix-1, iy)) / 2.0_num
        byc = (by(ix, iy) + by(ix, iy-1)) / 2.0_num
        b = SQRT(bxc**2 + byc**2 + bz(ix, iy)**2 + b_min**2)

        ! Direction of magnetic field
        ux(ix, iy) = bxc / b
        uy(ix, iy) = byc / b
=======
		! find factor required to convert between energy and temperature
		! N.B. only works for simple EOS
    e2t = (gamma - 1.0_num) / 2.0_num  
>>>>>>> conduction:src/core/conduct.f90

    DO iy = -1, ny + 1
      DO ix = -1, nx + 1 

        ! x face centred B field
        bxc = bx(ix, iy) 
        byc = (by(ix, iy) + by(ix+1, iy) + by(ix, iy-1) + by(ix+1, iy-1)) / 4.0_num
        bzc = (bz(ix, iy) + bz(ix+1,iy)) / 2.0_num
        bpx = SQRT(bxc**2 + byc**2 + bzc**2 + b_min**2) 

        exb = (energy(ix, iy) + energy(ix+1,iy)) / 2.0_num

        ! Direction of magnetic field on x face 
        ux = bxc / bpx       
        uy = byc / bpx       
        ! Kappa along magnetic field, now a vector  
        uxkx(ix, iy) = ux * ux * kappa_0 * (e2t * exb)**pow 
        uxky(ix, iy) = ux * uy * kappa_0 * (e2t * exb)**pow 
        ! add symmetic conduction near b=0 points 
        uxkx(ix,iy) = uxkx(ix,iy) + b_min**2 / (bpx**2 + b_min**2) * kappa_0 * (e2t * exb)**pow

        ! y face centred B field
        bxc = (bx(ix, iy) + bx(ix, iy+1) + bx(ix-1, iy) + bx(ix-1, iy+1) ) / 4.0_num
        byc = by(ix, iy)      
        bzc = (bz(ix,iy) + bz(ix,iy+1)) / 2.0_num
        bpy = SQRT(bxc**2 + byc**2 + bzc**2 + b_min**2)

        eyb = (energy(ix, iy) + energy(ix,iy+1)) / 2.0_num

        ! Direction of magnetic field on y face 
        uy = byc / bpy       
        ux = bxc / bpy      
        ! Kappa along magnetic field, now a vector  
        uykx(ix, iy) = uy * ux * kappa_0 * (e2t * eyb)**pow 
        uyky(ix, iy) = uy * uy * kappa_0 * (e2t * eyb)**pow    
        ! add symmetic conduction near b=0 points 
        uyky(ix,iy) = uyky(ix,iy) + b_min**2 / (bpy**2 + b_min**2) * kappa_0 * (e2t * eyb)**pow   

<<<<<<< HEAD:src/core/conduct.f90
				! Kappa isotropic
				kp(ix,iy) = kappa_0 * T * b_min**2 / b**2 
      END DO
    END DO 
    
    converged = .FALSE. 
    w = 1.6_num       ! initial over-relaxation parameter
		! store energy^{n} 
		energy0 = energy  
		! interate to get energy^{n+1} by SOR Guass-Seidel
    DO loop = 0, 100 
      energy1 = energy
=======
      END DO
    END DO  

    IF (heat_flux_limiter) THEN
       DO iy = 0, ny + 1
         DO ix = 0, nx + 1  
           ! estimate the parallel heat flux and the centre of a cell
           q_shx = (uxkx(ix,iy) + uxkx(ix-1,iy)) * e2t &
                * (energy(ix+1,iy) - energy(ix-1,iy)) / dxc(ix) 
           q_shy = (uyky(ix,iy) + uyky(ix,iy-1)) * e2t &
                * (energy(ix,iy+1) - energy(ix,iy-1)) / dyc(iy) 
           q_sh = SQRT(q_shx**2 + q_shy**2) / 16.0_num  
           ! estimate the free streaming limit
           ! 42.85 = SRQT(m_p/m_e)    
           q_f = 42.85_num * flux_limiter * rho(ix,iy) * (e2t * MIN(energy(ix,iy), temperature_100mk))**1.5_num            
           q_nl = 1.0_num / (1.0_num / MAX(q_sh, none_zero) + 1.0_num / MAX(q_f, none_zero)) 
           limiter(ix,iy) = q_nl / MAX(q_sh, none_zero) 
         END DO
      END DO  
      DO iy = 0, ny+1
        DO ix = 0, nx+1  
          uxkx(ix,iy) = uxkx(ix,iy) * (limiter(ix,iy) + limiter(ix+1,iy)) 
          uxky(ix,iy) = uxky(ix,iy) * (limiter(ix,iy) + limiter(ix+1,iy))
          uyky(ix,iy) = uyky(ix,iy) * (limiter(ix,iy) + limiter(ix,iy+1)) 
          uykx(ix,iy) = uykx(ix,iy) * (limiter(ix,iy) + limiter(ix,iy+1))
        END DO
      END DO
    END IF  
    
    converged = .FALSE. 
    w = 1.5_num       ! initial over-relaxation parameter  
		! store energy^{n} 
		energy0 = energy  
		! interate to get energy^{n+1} by SOR Guass-Seidel
    iterate: DO loop = 1, 100
>>>>>>> conduction:src/core/conduct.f90
      errmax = 0.0_num
      error = 0.0_num
      y1 = 1 
      DO redblack = 1, 2
        x1 = y1 
<<<<<<< HEAD:src/core/conduct.f90
        DO iy = 1, ny
          qpy = dyc(iy-1) / (dyc(iy) * (dyc(iy) + dyc(iy-1)))
          qmy = dyc(iy) / (dyc(iy-1) * (dyc(iy) + dyc(iy-1)))
          q0y = (dyc(iy)**2 - dyc(iy-1)**2) &
              / (dyc(iy) * dyc(iy-1) * (dyc(iy) + dyc(iy-1)))
          mpy = 1.0_num / (dyc(iy) * dyb(iy))
          mmy = 1.0_num / (dyc(iy-1) * dyb(iy))
          m0y = (dyc(iy) + dyc(iy-1)) / (dyc(iy) * dyc(iy-1) * dyb(iy))          
          DO ix = x1, nx, 2
            qpx = dxc(ix-1) / (dxc(ix) * (dxc(ix) + dxc(ix-1)))
            qmx = dxc(ix) / (dxc(ix-1) * (dxc(ix) + dxc(ix-1)))
            q0x = (dxc(ix)**2 - dxc(ix-1)**2) &
                / (dxc(ix) * dxc(ix-1) * (dxc(ix) + dxc(ix-1)))
            mpx = 1.0_num / (dxc(ix) * dxb(ix))
            mmx = 1.0_num / (dxc(ix-1) * dxb(ix))
            m0x = (dxc(ix) + dxc(ix-1)) / (dxc(ix) * dxc(ix-1) * dxb(ix))

            rx = qpx * e2temp(ix+1, iy) * energy(ix+1, iy) &
                - qmx * e2temp(ix-1, iy) * energy(ix-1, iy)
            ry = qpy * e2temp(ix, iy+1) * energy(ix, iy+1) &
                - qmy * e2temp(ix, iy-1) * energy(ix, iy-1)

            rxx = mpx * e2temp(ix+1, iy) * energy(ix+1, iy) &
                + mmx * e2temp(ix-1, iy) * energy(ix-1, iy)
            ryy = mpy * e2temp(ix, iy+1) * energy(ix, iy+1) &
                + mmy * e2temp(ix, iy-1) * energy(ix, iy-1)

            rxy = qpy * &
  									(qpx * e2temp(ix+1, iy+1) * energy1(ix+1, iy+1) &
  		              - qmx * e2temp(ix-1, iy+1) * energy1(ix-1, iy+1) &
  		              + q0x * e2temp(ix, iy+1) * energy(ix, iy+1)) &
                - qmy * &
  									(qpx * e2temp(ix+1, iy-1) * energy1(ix+1, iy-1) &
  		              - qmx * e2temp(ix-1, iy-1) * energy1(ix-1, iy-1) &
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

            residual = energy(ix, iy) &
                  - (energy0(ix, iy)  + a2) / (1.0_num + a1) 
            energy(ix, iy) = MAX(energy(ix, iy) - w * residual, 0.0_num) 
            error = ABS(residual / MAX(energy(ix, iy), energy0(ix, iy), none_zero))     
            errmax = MAX(errmax, error)
                
=======
        DO iy = 1, ny               
          DO ix = x1, nx, 2
            
            ! terms containing energy(ix,iy) resulting from 
            ! d^2/dx^2 and d^2/dy^2 derivatives
            a1 = uxkx(ix,iy)/(dxc(ix)*dxb(ix)) + uxkx(ix-1,iy)/(dxc(ix-1)*dxb(ix)) &
               + uyky(ix,iy)/(dyc(iy)*dyb(iy)) + uyky(ix,iy-1)/(dyc(iy-1)*dyb(iy))  
            
            ! terms not containing energy(ix,iy) resulting from 
            ! d^2/dx^2 and d^2/dy^2 derivatives
            a2 = uxkx(ix,iy)*e2t*energy(ix+1,iy)/(dxc(ix)*dxb(ix)) &
              + uxkx(ix-1,iy)*e2t*energy(ix-1,iy)/(dxc(ix-1)*dxb(ix))  
            a2 = a2 + uyky(ix,iy)*e2t*energy(ix,iy+1)/(dyc(iy)*dyb(iy)) &
              + uyky(ix,iy-1)*e2t*energy(ix,iy-1)/(dyc(iy-1)*dyb(iy))              
              
            ! terms not containing energy(ix,iy) resulting from 
            ! d^2/dxdy cross derivatives                  
            a2 = a2 + uxky(ix,iy) * e2t &
              * (energy(ix+1,iy+1) + energy(ix,iy+1) - energy(ix+1,iy-1) - energy(ix,iy-1)) &
              / (2.0_num * dxb(ix) * (dyc(iy) + dyc(iy-1)))  
            a2 = a2 - uxky(ix-1,iy) * e2t &
              * (energy(ix,iy+1) + energy(ix-1,iy+1) - energy(ix,iy-1) - energy(ix-1,iy-1)) &
              / (2.0_num * dxb(ix) * (dyc(iy) + dyc(iy-1)))  
  
            ! terms not containing energy(ix,iy) resulting from 
            ! d^2/dydx cross derivatives
            a2 = a2 + uykx(ix,iy) * e2t &
              * (energy(ix+1,iy+1) + energy(ix+1,iy) - energy(ix-1,iy+1) - energy(ix-1,iy)) &
              / (2.0_num * dyb(iy) * (dxc(ix) + dxc(ix-1)))  
            a2 = a2 - uykx(ix,iy-1) * e2t &
              * (energy(ix+1,iy) + energy(ix+1,iy-1) - energy(ix-1,iy) - energy(ix-1,iy-1)) &
              / (2.0_num * dyb(iy) * (dxc(ix) + dxc(ix-1)))  
  
            a1 = a1 * dt * e2t / rho(ix, iy)     
            a2 = a2 * dt / rho(ix, iy)         

            residual = energy(ix, iy) &
                  - (energy0(ix, iy)  + a2) / (1.0_num + a1) 
            energy(ix, iy) = MAX(energy(ix, iy) - w * residual, 0.0_num)
            error = ABS(residual / MAX(energy(ix, iy), energy0(ix,iy), none_zero))     
            errmax = MAX(errmax, error)
            
>>>>>>> conduction:src/core/conduct.f90
          END DO 
          x1 = 3 - x1
        END DO 
        y1 = 3 - y1
<<<<<<< HEAD:src/core/conduct.f90

        CALL energy_bcs

      END DO

      CALL MPI_ALLREDUCE(errmax, error, 1, mpireal, MPI_MAX, comm, errcode)
      errmax = error

=======
        
        CALL energy_bcs

      END DO
       
      CALL MPI_ALLREDUCE(errmax, error, 1, mpireal, MPI_MAX, comm, errcode)
      errmax = error      

>>>>>>> conduction:src/core/conduct.f90
      IF (errmax .LT. fractional_error) THEN
        converged = .TRUE.  
        EXIT iterate
      END IF
    END DO iterate
   
    IF (rank == 0 .AND. .NOT. converged) PRINT * , "Conduction failed at t = ", time

    a1 = 0.0_num
    DO iy = 1, ny 
      DO ix = 1, nx 
        a1 = a1 + energy(ix,iy) * dxb(ix) * dyb(iy)
      END DO
    END DO 
    CALL MPI_ALLREDUCE(a1, final_energy, 1, mpireal, MPI_SUM, comm, errcode)
    IF (initial_energy < final_energy) energy = energy * initial_energy / final_energy

<<<<<<< HEAD:src/core/conduct.f90
    DEALLOCATE(kx, ky)
    DEALLOCATE(ux, uy)
    DEALLOCATE(kp, e2temp, energy0, energy1)

  END SUBROUTINE conduct_heat
          


=======
    DEALLOCATE(uxkx, uxky)
    DEALLOCATE(uykx, uyky)
    DEALLOCATE(energy0)
    DEALLOCATE(limiter)

  END SUBROUTINE conduct_heat
          
>>>>>>> conduction:src/core/conduct.f90

END MODULE conduct
