MODULE conduct

  USE shared_data
  USE boundary    

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: conduct_heat    
  
  REAL(num) :: heat0 = 0.0_num
  REAL(num) :: rho_cor = 0.0_num 
  REAL(num), DIMENSION(:, :), ALLOCATABLE  :: temperature     

CONTAINS

  ! Subroutine implementing Braginskii parallel thermal conduction
  ! Note that this subroutine assumes that it is possible to
  ! Convert from specific internal energy to temperature by a constant
  ! factor.
  ! For this as well as other reasons, this routine doesn't work properly
  ! for partially ionised plasmas. 
	! Notation and algorithm in Appendix of Manual
  SUBROUTINE conduct_heat

    REAL(num), DIMENSION(:, :), ALLOCATABLE :: uxkx, uxky, uykx, uyky 
    REAL(num), DIMENSION(:, :), ALLOCATABLE :: energy0, limiter, temperature0
    REAL(num), DIMENSION(:, :), ALLOCATABLE :: radiation
    REAL(num), DIMENSION(:, :), ALLOCATABLE  :: heat_in 
    REAL(num), DIMENSION(:, :), ALLOCATABLE  :: alpha
    REAL(num) :: txb, tyb 
    REAL(num) :: b, bxc, byc, bzc, bpx, bpy 
    REAL(num) :: ux, uy
    REAL(num) :: pow = 5.0_num / 2.0_num  
    REAL(num) :: a1, a2, a3, a4, a5, error, errmax
    REAL(num) :: w, residual, q_shx, q_shy, q_sh, q_f, q_nl 
    REAL(num) :: rad_max, rad, alf

    INTEGER :: loop, redblack, x1, y1

    LOGICAL :: converged
    REAL(num), PARAMETER :: fractional_error = 1.0e-4_num
    REAL(num), PARAMETER :: b_min = 1.0e-3_num

    ALLOCATE(uxkx(-1:nx+1, -1:ny+1), uxky(-1:nx+1, -1:ny+1))
    ALLOCATE(uykx(-1:nx+1, -1:ny+1), uyky(-1:nx+1, -1:ny+1))
    ALLOCATE(energy0(-1:nx+2, -1:ny+2))
    ALLOCATE(limiter(-1:nx+2, -1:ny+2))
    ALLOCATE(radiation(1:nx, 1:ny))
    ALLOCATE(heat_in(1:nx, 1:ny))
    ALLOCATE(alpha(1:nx, 1:ny))
    ALLOCATE(temperature(-1:nx+2, -1:ny+2))
    ALLOCATE(temperature0(-1:nx+2, -1:ny+2))
            
    temperature = (gamma - 1.0_num) * (energy - (1.0_num - xi_n) * ionise_pot) &
                 / (2.0_num - xi_n)            

    DO iy = -1, ny + 1
      DO ix = -1, nx + 1 

        ! x face centred B field
        bxc = bx(ix, iy) 
        byc = (by(ix, iy) + by(ix+1, iy) + by(ix, iy-1) + by(ix+1, iy-1)) / 4.0_num
        bzc = (bz(ix, iy) + bz(ix+1,iy)) / 2.0_num
        bpx = SQRT(bxc**2 + byc**2 + bzc**2 + b_min**2)
        bpx = MAX(bpx, none_zero) 

        txb = (temperature(ix, iy) + temperature(ix+1,iy)) / 2.0_num

        ! Direction of magnetic field on x face 
        ux = bxc / bpx       
        uy = byc / bpx       
        ! Kappa along magnetic field, now a vector  
        uxkx(ix, iy) = ux * ux * kappa_0 * txb**pow 
        uxky(ix, iy) = ux * uy * kappa_0 * txb**pow 
        ! add symmetic conduction near b=0 points 
        uxkx(ix,iy) = uxkx(ix,iy) + b_min**2 * kappa_0 * txb**pow &
              / (bpx**2 + b_min**2)

        ! y face centred B field
        bxc = (bx(ix, iy) + bx(ix, iy+1) + bx(ix-1, iy) + bx(ix-1, iy+1) ) / 4.0_num
        byc = by(ix, iy)      
        bzc = (bz(ix,iy) + bz(ix,iy+1)) / 2.0_num
        bpy = SQRT(bxc**2 + byc**2 + bzc**2 + b_min**2)
        bpy = MAX(bpy, none_zero)

        tyb = (temperature(ix, iy) + temperature(ix,iy+1)) / 2.0_num

        ! Direction of magnetic field on y face 
        uy = byc / bpy       
        ux = bxc / bpy      
        ! Kappa along magnetic field, now a vector  
        uykx(ix, iy) = uy * ux * kappa_0 * tyb**pow 
        uyky(ix, iy) = uy * uy * kappa_0 * tyb**pow    
        ! add symmetic conduction near b=0 points 
        uyky(ix,iy) = uyky(ix,iy) + b_min**2 * kappa_0 * tyb**pow &
            / (bpy**2 + b_min**2)  

      END DO
    END DO  

    IF (heat_flux_limiter) THEN
       DO iy = 0, ny + 1
         DO ix = 0, nx + 1  
           ! estimate the parallel heat flux and the centre of a cell
           q_shx = &
                (uxkx(ix,iy) + uxkx(ix-1,iy))   &
                * (temperature(ix+1,iy) - temperature(ix-1,iy)) / dxc(ix) &
              + (uxky(ix,iy) + uxky(ix,iy-1))  &
                * (temperature(ix,iy+1) - temperature(ix,iy-1)) / dyc(iy)
           q_shy = &
                (uykx(ix,iy) + uykx(ix-1,iy))   &
                * (temperature(ix+1,iy) - temperature(ix-1,iy)) / dxc(ix) &
             +  (uyky(ix,iy) + uyky(ix,iy-1))   &
                * (temperature(ix,iy+1) - temperature(ix,iy-1)) / dyc(iy) 
           q_sh = SQRT(q_shx**2 + q_shy**2) / 16.0_num  
           ! estimate the free streaming limit
           ! 42.85 = SRQT(m_p/m_e)    
           q_f = 42.85_num * flux_limiter * rho(ix,iy) &
            * MIN(temperature(ix,iy), temperature_100mk)**1.5_num            
           q_nl = 1.0_num / (1.0_num / MAX(q_sh, none_zero) + 1.0_num / MAX(q_f, none_zero)) 
           limiter(ix,iy) = q_nl / MAX(q_sh, none_zero) / 2.0_num 
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
    w = 1.6_num       ! initial over-relaxation parameter  
		! store energy^{n} 
		energy0 = energy  
		temperature0 = temperature  

    radiation = 0.0_num
    heat_in = 0.0_num 
    alpha = 0.0_num
    DO iy = 1, ny
      DO ix = 1, nx  
        heat_in(ix,iy) = heating(rho(ix,iy), temperature0(ix,iy))  
        CALL rad_losses(rho(ix,iy), temperature0(ix,iy), xi_n(ix,iy), rad, alf)
        alpha(ix,iy) = alf  
        IF (yc(iy) > 1.0_num) THEN  
          radiation(ix,iy) =  rad 
        END IF    
      END DO
    END DO      

		! interate to get energy^{n+1} by SOR Guass-Seidel
    iterate: DO loop = 1, 100
      errmax = 0.0_num 
      error = 0.0_num
      y1 = 1 
      DO redblack = 1, 2
        x1 = y1 
        DO iy = 1, ny               
          DO ix = x1, nx, 2
            
            ! terms containing energy(ix,iy) resulting from 
            ! d^2/dx^2 and d^2/dy^2 derivatives
            a1 = uxkx(ix,iy)/(dxc(ix)*dxb(ix)) + uxkx(ix-1,iy)/(dxc(ix-1)*dxb(ix)) &
               + uyky(ix,iy)/(dyc(iy)*dyb(iy)) + uyky(ix,iy-1)/(dyc(iy-1)*dyb(iy))   
                           
            ! terms not containing temperature(ix,iy) resulting from 
            ! d^2/dx^2 and d^2/dy^2 derivatives
            a2 = uxkx(ix,iy)*temperature(ix+1,iy)/(dxc(ix)*dxb(ix)) &
              + uxkx(ix-1,iy)*temperature(ix-1,iy)/(dxc(ix-1)*dxb(ix))  
            a2 = a2 + uyky(ix,iy)*temperature(ix,iy+1)/(dyc(iy)*dyb(iy)) &
              + uyky(ix,iy-1)*temperature(ix,iy-1)/(dyc(iy-1)*dyb(iy))              
              
            ! terms not containing temperature(ix,iy) resulting from 
            ! d^2/dxdy cross derivatives                  
            a2 = a2 + uxky(ix,iy)  &
              * (temperature(ix+1,iy+1) + temperature(ix,iy+1) - temperature(ix+1,iy-1) - temperature(ix,iy-1)) &
              / (2.0_num * dxb(ix) * (dyc(iy) + dyc(iy-1)))  
            a2 = a2 - uxky(ix-1,iy)  &
              * (temperature(ix,iy+1) + temperature(ix-1,iy+1) - temperature(ix,iy-1) - temperature(ix-1,iy-1)) &
              / (2.0_num * dxb(ix) * (dyc(iy) + dyc(iy-1)))  
  
            ! terms not containing temperature(ix,iy) resulting from 
            ! d^2/dydx cross derivatives
            a2 = a2 + uykx(ix,iy)  &
              * (temperature(ix+1,iy+1) + temperature(ix+1,iy) - temperature(ix-1,iy+1) - temperature(ix-1,iy)) &
              / (2.0_num * dyb(iy) * (dxc(ix) + dxc(ix-1)))  
            a2 = a2 - uykx(ix,iy-1)  &
              * (temperature(ix+1,iy) + temperature(ix+1,iy-1) - temperature(ix-1,iy) - temperature(ix-1,iy-1)) &
              / (2.0_num * dyb(iy) * (dxc(ix) + dxc(ix-1))) 

            a3 = (a1 + radiation(ix,iy) * alpha(ix,iy)) * temperature(ix,iy) / energy(ix,iy) 

            a4 = a2 + heat_in(ix,iy)  - (1.0_num - alpha(ix,iy)) * radiation(ix,iy) * temperature0(ix,iy)
  
            a3 = a3 * dt / rho(ix, iy)     
            a4 = a4 * dt / rho(ix, iy)         

            residual = energy(ix, iy) &
                  - (energy0(ix, iy)  + a4) / (1.0_num + a3) 
            energy(ix, iy) = MAX(energy(ix, iy) - w * residual, (1.0_num - xi_n(ix,iy)) * ionise_pot) 
            error = ABS(residual) / energy0(ix,iy)
            errmax = MAX(errmax, error)
            
          END DO 
          x1 = 3 - x1
        END DO 
        y1 = 3 - y1    
        
        CALL energy_bcs  
        
        temperature = (gamma - 1.0_num) * (energy - (1.0_num - xi_n) * ionise_pot) &
                 / (2.0_num - xi_n)            
        
      END DO       
             
      CALL MPI_ALLREDUCE(errmax, error, 1, mpireal, MPI_MAX, comm, errcode)
      errmax = error      

      IF (errmax .LT. fractional_error) THEN
        converged = .TRUE.  
        EXIT iterate
      END IF

    END DO iterate
  
    IF (rank == 0 .AND. .NOT. converged) PRINT * , "Conduction failed at t = ", time

    DEALLOCATE(uxkx, uxky)
    DEALLOCATE(uykx, uyky)
    DEALLOCATE(energy0)
    DEALLOCATE(limiter)
    DEALLOCATE(radiation)
    DEALLOCATE(heat_in)
    DEALLOCATE(alpha)
    DEALLOCATE(temperature)
    DEALLOCATE(temperature0)

  END SUBROUTINE conduct_heat
               
               

  SUBROUTINE rad_losses(density, temperature, xi, rad, alf)  
    ! returns the normalised RTV losses divided by normalised temperature
    REAL(num), INTENT(IN) :: density, temperature, xi  
    REAL(num), INTENT(OUT) :: rad, alf 

    REAL(num), DIMENSION(7) :: trange = (/0.02_num,0.0398_num,0.0794_num,0.251_num, 0.562_num,1.995_num,10.0_num /)
    REAL(num), DIMENSION(6) :: psi = (/1.2303_num, 870.96_num, 5.496_num, 0.3467_num, 1.0_num, 1.6218_num /)
    REAL(num), DIMENSION(6) :: alpha = (/0.0_num, 2.0_num, 0.0_num, -2.0_num, 0.0_num, -2.0_num/3.0_num /) 
    REAL(num) :: tmk, factor
    INTEGER :: i
                      
    rad = 0.0_num     
    alf = 0.0_num                           
    IF(.NOT. radiation) RETURN

    tmk = temperature * t2tmk   
    IF (tmk < trange(1) .OR. tmk > trange(7)) RETURN
                                                      
    DO i = 1, 6
      IF (tmk >= trange(i) .AND. tmk <= trange(i+1)) EXIT
    END DO 
                                                
    factor = (1.0_num - xi)**2
    IF (eos_number == EOS_IDEAL) factor = 1.0_num
    rad = factor * density**2 * psi(i) * tmk**(alpha(i)-1.0_num)  
    rad = rad * h_star * lr_star * t2tmk   
    alf = alpha(i)

  END SUBROUTINE rad_losses  
       
  
  
  FUNCTION heating(density, t0) 
    REAL(num), INTENT(IN) :: density, t0
    REAL(num) :: heating
    REAL(num) :: tmk, a1, a2, rad, alf, height
    LOGICAL :: first_call = .TRUE.
    REAL(num) :: heat0 = 0.0_num      
    REAL(num) :: rho_coronal = 0.0_num   
    INTEGER :: loop
    
    IF (first_call) THEN
      a1 = 0.0_num     
      IF (proc_y_max == MPI_PROC_NULL) THEN    
         CALL rad_losses(rho(1,ny), temperature(1,ny), xi_n(1,ny), rad, alf)
         a1 = rad * temperature(1,ny) / rho(1,ny)**2 
      END IF               
      CALL MPI_ALLREDUCE(a1, heat0, 1, mpireal, MPI_MAX, comm, errcode)  
       
      ! choose a reference density based on height   
      height = 15.0_num  
      a2 = 0.0_num
      DO loop = 1, ny
         IF (yb(loop) >= height .AND. yb(loop-1) < height) a2 = rho(1,loop)  
      END DO 
      CALL MPI_ALLREDUCE(a2, rho_coronal, 1, mpireal, MPI_MAX, comm, errcode)    

      first_call = .FALSE.               
    END IF     

    heating = 0.0_num
    IF (.NOT. coronal_heating) RETURN

    tmk = t0 * t2tmk                       
    IF(density < rho_coronal .AND. tmk > 0.02_num) heating = 100.0_num * heat0 * density**2

  END FUNCTION heating
  
                

!  SUBROUTINE get_temperature                
! Full iterative solution for T and xi_n below has been tested and gives same answer
! as simple approach above for flux emergence to within 0.2%
! Should routinely test on new problems to be sure.
! 
!     REAL(num) :: bof, r, T, rho0, e0, dx, x
!     REAL(num), DIMENSION(2) :: ta, fa, xi_a
!     INTEGER :: loop
!                                                                      
!     xi_n = 1.0_num                                                                 
!     IF (eos_number == EOS_IDEAL .AND. (.NOT. neutral_gas)) xi_n = 0.0_num    
!     temperature = (gamma - 1.0_num) * energy / (2.0_num - xi_n)
!     IF (eos_number == EOS_IDEAL) RETURN
! 
!     ! Variable bof is b / f in the original version
!     DO iy = -1, ny + 2
!       DO ix = -1, nx + 2
!         rho0 = rho(ix, iy)
!         e0 = energy(ix, iy)
!         ta = (gamma - 1.0_num) &
!             * (/ MAX((e0 - ionise_pot) / 2.0_num, none_zero), e0 /)
! 
!         IF (ta(1) > ta(2)) THEN
!           PRINT * , "Temperature bounds problem", ta
!           STOP
!         END IF
! 
!         dx = ta(2) - ta(1)
!         t = ta(1)
! 
!         DO loop = 1, 100
!           dx = dx / 2.0_num
!           x = t  + dx
!           xi_a(1) = get_neutral(x, rho0, yb(iy))   
!           fa(1) = x - (gamma - 1.0_num) * (e0 &
!               - (1.0_num - xi_a(1)) * ionise_pot) / (2.0_num - xi_a(1))
!           IF (fa(1) <= 0.0_num) t = x
!           IF (ABS(dx) < 1.e-8_num .OR. fa(1) == 0.0_num) EXIT
!         END DO
! 
!         xi_n(ix, iy) = get_neutral(x, rho0, yb(iy))  
!         temperature(ix,iy) = x 
!       END DO
!      END DO                   
!
!  END SUBROUTINE get_temperature
                                     


END MODULE conduct

