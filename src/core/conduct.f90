MODULE conduct

  USE shared_data
  USE boundary

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: conduct_heat
  
  INTEGER :: n_s_stages
  REAL(num),PARAMETER :: pow = 5.0_num/2.0_num
  REAL(num),PARAMETER :: min_b = 1.0e-6_num

CONTAINS


  !****************************************************************************
  ! Subroutine calculating the number of stages needed for the RKL2
  !****************************************************************************

  SUBROUTINE s_stages

    REAL (num)  ::  stages, dt_parab, dt1, dt2
    REAL(num) :: gm1
    INTEGER :: n_s_stages_local

    dt_parab = 1.e10_num
    gm1 = 0.5_num * (gamma - 1.0_num)

    DO iy = 1, ny
      DO ix = 1, nx
        ! estimate explicit thermal conduction time-step
        dt1 = gm1 * rho(ix,iy)*dxb(ix)**2 / (kappa_0 * (gm1 * energy(ix,iy))**pow)
        dt2 = gm1 * rho(ix,iy)*dyb(iy)**2 / (kappa_0 * (gm1 * energy(ix,iy))**pow)
        dt_parab = MIN(dt_parab, dt1, dt2)
      END DO
    END DO

    dt_parab = dt_multiplier * dt_parab

    stages = 0.5_num*(SQRT(9.0_num + 16.0_num * (dt/dt_parab))-1.0_num)
    n_s_stages = CEILING(stages)
    IF (MODULO(n_s_stages,2) .EQ. 0) THEN 
      n_s_stages = n_s_stages + 1
    ENDIF
    n_s_stages_local = n_s_stages
    CALL MPI_ALLREDUCE(n_s_stages_local, n_s_stages, 1, &
        MPI_INTEGER, MPI_MAX, comm, errcode)

  END SUBROUTINE s_stages



    !****************************************************************************
    ! Subroutine to calculate the heat flux
    !****************************************************************************

    SUBROUTINE heat_flux(temperature, flux)
      REAL(num),INTENT(IN),DIMENSION(-1:,-1:) :: temperature
      REAL(num),INTENT(OUT),DIMENSION(-1:,-1:) :: flux
      INTEGER :: ix, ixp, ixm
      REAL(num) :: tb1, tb2, tg1, tg2, fc_sp1, fc_sp2, rho_b1, rho_b2
      REAL(num) :: tg_a1, tg_a2, tb_p2, tb_p1, tb_m1, tb_m2
      REAL(num) :: fc_sa1, fc_sa2, fc1, fc2, modb1, modb2
      REAL(num) :: byf1, byf2, bxf1, bxf2, bzf1, bzf2

      DO iy = 1,ny
        iyp = iy + 1
        iym = iy - 1
        DO ix = 1,nx
          ixp = ix + 1
          ixm = ix - 1
  
          !X flux
          byf1 = 0.25_num*(by(ixm,iy)+by(ix,iy)+by(ixm,iym)+by(ix,iym))
          byf2 = 0.25_num*(by(ix,iy)+by(ixp,iy)+by(ix,iym)+by(ixp,iym))
          bzf1 = 0.5_num*(bz(ix,iy)+bz(ixm,iy))
          bzf2 = 0.5_num*(bz(ix,iy)+bz(ixp,iy))

          modb2 = SQRT(bx(ix,iy)**2+byf2**2+bzf2**2)
          modb1 = SQRT(bx(ixm,iy)**2+byf1**2+bzf1**2)
          
          !Braginskii Conductive Flux
          !Temperature at the x boundaries in the current cell
          tb2 = (temperature(ixp,iy) + temperature(ix,iy))/2.0_num
          tb1 = (temperature(ix,iy) + temperature(ixm,iy))/2.0_num
  
          !Temperature at the x boundaries in the cell above		  
          tb_p2 = (temperature(ixp,iyp)+temperature(ix,iyp))/2.0_num
          tb_p1 = (temperature(ix,iyp)+temperature(ixm,iyp))/2.0_num
          !Temperature at the x boundaries in the cell below		  
          tb_m2 = (temperature(ixp,iym)+temperature(ix,iym))/2.0_num
          tb_m1 = (temperature(ix,iym)+temperature(ixm,iym))/2.0_num
          !X temperature gradient at the x boundaries of the current cell		  
          tg2 = (temperature(ixp,iy) - temperature(ix,iy))/dxc(ix)
          tg1 = (temperature(ix,iy) - temperature(ixm,iy))/dxc(ixm)
          !Y temperature gradient at the x boundaries of the current cell
          !Uses centred difference on averaged values, so likely very smoothed
          tg_a2 = (tb_p2 - tb2) * dyc(iym) / dyc(iy) + (tb2 - tb_m2) * dyc(iy) / dyc(iym) &
                / (dyc(iy) + dyc(iym))
          tg_a1 = (tb_p1 - tb1) * dyc(iym) / dyc(iy) + (tb1 - tb_m1) * dyc(iy) / dyc(iym) &
                / (dyc(iy) + dyc(iym))

          fc_sp2 = kappa_0 * tb2**pow * (bx(ix,iy) * (tg2 * bx(ix,iy) + &
              tg_a2 * byf2)+tg2*min_b)/(modb2**2+min_b) 
          fc_sp1 = kappa_0 * tb1**pow * (bx(ixm,iy) * (tg1 * bx(ixm,iy) + &
              tg_a1 * byf1)+tg1*min_b)/(modb1**2+min_b)

          ! Saturated Conductive Flux
          rho_b2 = (rho(ixp,iy)+rho(ix,iy))/2.0_num
          rho_b1 = (rho(ix,iy)+rho(ixm,iy))/2.0_num
          fc_sa2 = 64.27_num * rho_b2 * tb2**(3.0_num/2.0_num)  !64.27 = (3/2)*SRQT(m_p/m_e)
          fc_sa1 = 64.27_num * rho_b1 * tb1**(3.0_num/2.0_num)

          ! Conductive Flux Limiter. Note flux_limiter is inverse of usual definition here
          fc2 = 1.0_num / (1.0_num/fc_sp2 + flux_limiter/fc_sa2)
          fc1 = 1.0_num / (1.0_num/fc_sp1 + flux_limiter/fc_sa1)

          flux(ix,iy) = (fc2-fc1)/dxb(ix)

          !Y flux
          bxf1 = 0.25_num*(bx(ix,iy)+bx(ix,iyp)+bx(ixm,iyp)+bx(ixm,iy))
          bxf2 = 0.25_num*(bx(ix,iym)+bx(ix,iy)+bx(ixm,iy)+bx(ixm,iym))
          bzf1 = 0.5_num*(bz(ix,iy)+bz(ix,iym))
          bzf2 = 0.5_num*(bz(ix,iy)+bz(ix,iyp))
          modb1 = SQRT(by(ix,iy)**2+bxf1**2+bzf1**2)
          modb2 = SQRT(by(ix,iym)**2+bxf2**2+bzf2**2)

          ! Braginskii Conductive Flux
          tb2 = (temperature(ix,iyp) + temperature(ix,iy))/2.0_num
          tb1 = (temperature(ix,iy) + temperature(ix,iym))/2.0_num

          !Temperature at the y boundaries in the cell right		  
          tb_p2 = (temperature(ixp,iyp)+temperature(ixp,iy))/2.0_num
          tb_p1 = (temperature(ixp,iy)+temperature(ixp,iym))/2.0_num
          !Temperature at the y boundaries in the cell left		  
          tb_m2 = (temperature(ixm,iyp)+temperature(ixm,iy))/2.0_num
          tb_m1 = (temperature(ixm,iy)+temperature(ixm,iym))/2.0_num
          !Y temperature gradient at the y boundaries of the current cell		  
          tg2 = (temperature(ix,iyp) - temperature(ix,iy))/dyc(iy)
          tg1 = (temperature(ix,iy) - temperature(ix,iym))/dyc(iym)
          !X temperature gradient at the y boundaries of the current cell
          !Uses centred difference on averaged values, so likely very smoothed
          tg_a2 = (tb_p2 - tb2) * dxc(ixm) / dxc(ix) + (tb2 - tb_m2) * dxc(ix) / dxc(ixm) &
                / (dxc(ix) + dyc(ixm))
          tg_a1 = (tb_p1 - tb1) * dxc(ixm) / dxc(ix) + (tb1 - tb_m1) * dxc(ix) / dxc(ixm) &
                / (dxc(ix) + dxc(ixm))              
 
          fc_sp2 = kappa_0 * tb2**pow * (by(ix,iy) * (tg2 * by(ix,iy)&
              + tg_a2 * bxf2)+min_b*tg2)/(modb2**2+min_b)
          fc_sp1 = kappa_0 * tb1**pow * (by(ix,iym) * (tg1 * by(ix,iym)&
              + tg_a1 * bxf1)+min_b*tg1)/(modb1**2+min_b)

          ! Saturated Conductive Flux     
          rho_b2 = (rho(ix,iyp)+rho(ix,iy))/2.0_num
          rho_b1 = (rho(ix,iy)+rho(ix,iym))/2.0_num
          fc_sa2 = 64.27_num * rho_b2 * tb2**(3.0_num/2.0_num)  !64.27 = (3/2)*SRQT(m_p/m_e)
          fc_sa1 = 64.27_num * rho_b1 * tb1**(3.0_num/2.0_num)

          ! Conductive Flux Limiter
          fc2 = 1.0_num / (1.0_num/fc_sp2 + flux_limiter/fc_sa2)
          fc1 = 1.0_num / (1.0_num/fc_sp1 + flux_limiter/fc_sa1)

          flux(ix,iy) = flux(ix,iy) + (fc2-fc1)/dyb(iy)
        END DO
      END DO

  END SUBROUTINE heat_flux



  !****************************************************************************
  ! Subroutine implementing Braginskii parallel thermal conduction.
  ! Notation and algorithm in Appendix of Manual
  !****************************************************************************

  SUBROUTINE conduct_heat
  
    CALL s_stages
    CALL heat_conduct_sts2

  END SUBROUTINE conduct_heat


  !****************************************************************************
  ! Implementation of the RKL2 scheme
  !****************************************************************************

  SUBROUTINE heat_conduct_sts2

    !Superstepping based conduction code
    !2nd order Runge-Kutta-Lagrange (RKL2) scheme
    !Based on Meyer at al. 2012 variables named as in that paper
    REAL(num) :: pow=5.0_num/2.0_num
    REAL(num), DIMENSION(:,:), ALLOCATABLE  :: flux
    REAL(num)  ::  omega_1
    REAL(num), DIMENSION(0:n_s_stages)  :: a, b
    REAL(num), DIMENSION(1:n_s_stages)  :: mu_tilde
    REAL(num), DIMENSION(2:n_s_stages)  :: mu, nu, gamma_tilde
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE  :: Y     ! intermediate solutions Y=[Y_0, Y_j-2, Y_j-1 , Y_j]
    REAL(num)  ::  tb1, tb2, tg1, tg2                  ! temperature gradient and temperature on the boundary
    REAL(num)  ::  rho_b1, rho_b2                      ! density on the boundary
    REAL(num)  ::  c0, c1
    REAL(num)  ::  Lc_Yj_1                  ! L^c(Y_j-1)
    REAL(num), DIMENSION(1:nx,1:ny)  :: Lc_Y0    ! L^c(Y_0)
    REAL(num)  ::  Fc_sp1, Fc_sp2, Fc_sa1, Fc_sa2, Fc1, Fc2 ! spitzer and saturated conductive flux terms
    INTEGER :: j

    ALLOCATE(flux(-1:nx+2,-1:ny+2))
    ALLOCATE(Y(0:3,-1:nx+2,-1:ny+2))
    flux = 0.0_num
    Y = 0.0_num
    
    omega_1 = 4.0_num/(DBLE(n_s_stages)**2.0_num + DBLE(n_s_stages) - 2.0_num)
    b(0:2) = 1.0_num/3.0_num
    DO j = 3, n_s_stages
       b(j) = (DBLE(j)**2.0_num + DBLE(j)- 2.0_num) / &
           (2.0_num * DBLE(j) *(DBLE(j) + 1.0_num))
    ENDDO 
    a = 1.0_num - b
    mu_tilde(1) = omega_1/3.0_num

    DO j = 2,n_s_stages
      mu_tilde(j) = ((2.0_num *DBLE(j) - 1.0_num) /DBLE(j))* &
          omega_1 *(b(j)/b(j-1))
      mu(j) = ((2.0_num *DBLE(j) - 1.0_num) /DBLE(j)) *(b(j)/b(j-1))
      nu(j) = -1.0_num * ((DBLE(j) - 1.0_num) /DBLE(j))*(b(j)/b(j-2))
      gamma_tilde(j) = -1.0_num*a(j-1)*mu_tilde(j)
    ENDDO

    !!!!! RKL2 s stage scheme !!!!!
    ! initial condition
    Y(0,:,:) = energy
    !! boundary conditions 
    DO j = 1,3
      Y(j,:,:) = energy
    END DO


    !! First STS stage
    CALL heat_flux(Y(0,:,:) * ((gamma - 1.0_num)/2.0_num),flux)
    DO iy = 1,ny
      DO ix = 1,nx 
        ! Store L^c(Y_0) to use in stages 2-s
        Lc_Y0(ix,iy) = (1.0_num / rho(ix,iy)) * flux(ix,iy)
        c0 =  mu_tilde(1) * dt * Lc_Y0(ix,iy)
        Y(1,ix,iy) = Y(0,ix,iy) + c0
      ENDDO
    ENDDO

    Y(2,1:nx,1:ny) = Y(1,1:nx,1:ny)
    Y(1,1:nx,1:ny) = Y(0,1:nx,1:ny)

    DO j = 2,n_s_stages
      CALL heat_flux(Y(2,:,:) * ((gamma - 1.0_num)/2.0_num),flux)
      DO iy = 1,ny
        DO ix = 1,nx
          c0 = gamma_tilde(j) * dt * Lc_Y0(ix,iy)
          Lc_Yj_1 = (1.0_num / rho(ix,iy)) * flux(ix,iy)
          c1 = mu_tilde(j) * dt * Lc_Yj_1
        !Y_j
          Y(3,ix,iy) = mu(j)*Y(2,ix,iy) + nu(j)*Y(1,ix,iy) + (1.0_num -mu(j)-nu(j))* Y(0,ix,iy)
          Y(3,ix,iy) = Y(3,ix,iy) + c1 + c0
        END DO
      END DO

      IF ( j .LT. n_s_stages ) THEN
        !for jth stage  Y=[Y_0, Y_j-2, Y_j-1 , Y_j]
        Y(1,:,:) = Y(2,:,:)
        !This is not ideal, but it allows you to not have special boundary conditions
        energy = Y(3,:,:)
        CALL energy_bcs
        Y(2,:,:) = energy
        Y(3,1:nx,1:ny) = 0.0_num
      END IF
    ENDDO

    energy(1:nx,1:ny) = Y(3,1:nx,1:ny)
    CALL energy_bcs

  END SUBROUTINE heat_conduct_sts2




  SUBROUTINE rad_losses(density, temperature, xi, height, rad, alf)

    ! Returns the normalised RTV losses

    REAL(num), INTENT(IN) :: density, temperature, xi, height
    REAL(num), INTENT(OUT) :: rad, alf

    REAL(num), DIMENSION(7) :: trange = (/0.02_num, 0.0398_num, 0.0794_num, &
        0.251_num, 0.562_num, 1.995_num, 10.0_num/)
    REAL(num), DIMENSION(6) :: psi = (/1.2303_num, 870.96_num, 5.496_num, &
        0.3467_num, 1.0_num, 1.6218_num/)
    REAL(num), DIMENSION(6) :: alpha = (/0.0_num, 2.0_num, 0.0_num, &
        -2.0_num, 0.0_num, -2.0_num / 3.0_num/)
    REAL(num) :: tmk, factor
    INTEGER :: i

    rad = 0.0_num
    alf = 0.0_num

    IF (.NOT. radiation) RETURN
    IF (height < 1.0_num) RETURN

    tmk = temperature * t2tmk
    IF (tmk < trange(1) .OR. tmk > trange(7)) RETURN

    DO i = 1, 6
      IF (tmk >= trange(i) .AND. tmk <= trange(i+1)) EXIT
    END DO

    ! Account for reduced electron number density due to neutrals
    factor = (1.0_num - xi)**2
    IF (eos_number == EOS_IDEAL) factor = 1.0_num

    rad = factor * density**2 * psi(i) * tmk**alpha(i)
    rad = rad * h_star * lr_star
    alf = alpha(i)

  END SUBROUTINE rad_losses



  FUNCTION heating(density, t0)

    ! For a given density and temperature returns a user specific
    ! heating function

    REAL(num), INTENT(IN) :: density, t0
    REAL(num) :: heating
    REAL(num) :: tmk, a1, a2, rad, alf, height
    LOGICAL, SAVE :: first_call = .TRUE.
    REAL(num) :: heat0 = 0.0_num
    REAL(num) :: rho_coronal = 0.0_num
    LOGICAL :: store_state
    INTEGER :: loop

    IF (first_call) THEN
      a1 = 0.0_num
      IF (proc_y_max == MPI_PROC_NULL) THEN
        store_state = radiation
        radiation = .TRUE.
!        CALL rad_losses(rho(1,ny), temperature(1,ny), &
!                        xi_n(1,ny), yc(ny), rad, alf)
        radiation = store_state
        a1 = rad / rho(1,ny)**2
      END IF

      CALL MPI_ALLREDUCE(a1, heat0, 1, mpireal, MPI_MAX, comm, errcode)

      ! Choose a reference density based on height
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
    ! For low density and high temeprature define a heating course term
    IF (density < rho_coronal .AND. tmk > 0.02_num) &
        heating = 100.0_num * heat0 * density**2

  END FUNCTION heating

END MODULE conduct
