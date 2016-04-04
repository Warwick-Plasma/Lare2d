MODULE conduct

  USE shared_data
  USE boundary

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: conduct_heat

  INTEGER :: n_s_stages
  REAL(num), PARAMETER :: pow = 2.5_num
  REAL(num), PARAMETER :: min_b = 1.0e-6_num

CONTAINS

  !****************************************************************************
  ! Subroutine calculating the number of stages needed for the RKL2
  !****************************************************************************

  SUBROUTINE s_stages

    REAL(num) :: stages, dt_parab, dt1, dt2
    REAL(num) :: gm1, temp
    INTEGER :: n_s_stages_local

    dt_parab = 1.e10_num
    gm1 = 0.5_num * (gamma - 1.0_num)

    DO iy = 1, ny
      DO ix = 1, nx
        ! Estimate explicit thermal conduction time-step
        temp = gm1 * energy(ix,iy)
        temp = gm1 * rho(ix,iy) / (kappa_0 * temp**pow)
        dt1 = temp * dxb(ix)**2
        dt2 = temp * dyb(iy)**2
        dt_parab = MIN(dt_parab, dt1, dt2)
      END DO
    END DO

    dt_parab = dt_multiplier * dt_parab

    stages = 0.5_num * (SQRT(9.0_num + 16.0_num * (dt / dt_parab)) - 1.0_num)
    n_s_stages_local = CEILING(stages)

    IF (MODULO(n_s_stages_local,2) == 0) THEN
      n_s_stages_local = n_s_stages_local + 1
    END IF

    CALL MPI_ALLREDUCE(n_s_stages_local, n_s_stages, 1, MPI_INTEGER, MPI_MAX, &
        comm, errcode)

  END SUBROUTINE s_stages



  !****************************************************************************
  ! Subroutine to calculate the heat flux
  !****************************************************************************

  SUBROUTINE heat_flux(temperature, flux)

    REAL(num), INTENT(IN), DIMENSION(-1:,-1:) :: temperature
    REAL(num), INTENT(OUT), DIMENSION(-1:,-1:) :: flux
    INTEGER :: ix, ixp, ixm
    REAL(num) :: tb1, tb2, tg1, tg2, fc_sp1, fc_sp2, rho_b1, rho_b2
    REAL(num) :: tg_a1, tg_a2, tb_p2, tb_p1, tb_m1, tb_m2
    REAL(num) :: fc_sa1, fc_sa2, fc1, fc2, modb1, modb2
    REAL(num) :: byf1, byf2, bxf1, bxf2, bzf1, bzf2

    DO iy = 1, ny
      iyp = iy + 1
      iym = iy - 1
      DO ix = 1, nx
        ixp = ix + 1
        ixm = ix - 1

        ! X flux
        byf1 = 0.25_num * (by(ix,iy) + by(ix,iym) + by(ixm,iy) + by(ixm,iym))
        byf2 = 0.25_num * (by(ix,iy) + by(ix,iym) + by(ixp,iy) + by(ixp,iym))
        bzf1 = 0.5_num * (bz(ix,iy) + bz(ixm,iy))
        bzf2 = 0.5_num * (bz(ix,iy) + bz(ixp,iy))

        modb1 = SQRT(bx(ixm,iy)**2 + byf1**2 + bzf1**2)
        modb2 = SQRT(bx(ix ,iy)**2 + byf2**2 + bzf2**2)

        ! Braginskii Conductive Flux
        ! Temperature at the x boundaries in the current cell
        tb1 = 0.5_num * (temperature(ix,iy) + temperature(ixm,iy))
        tb2 = 0.5_num * (temperature(ix,iy) + temperature(ixp,iy))

        ! Temperature at the x boundaries in the cell above		
        tb_p1 = 0.5_num * (temperature(ix,iyp) + temperature(ixm,iyp))
        tb_p2 = 0.5_num * (temperature(ix,iyp) + temperature(ixp,iyp))
        ! Temperature at the x boundaries in the cell below		
        tb_m1 = 0.5_num * (temperature(ix,iym) + temperature(ixm,iym))
        tb_m2 = 0.5_num * (temperature(ix,iym) + temperature(ixp,iym))
        ! X temperature gradient at the x boundaries of the current cell		
        tg1 = (temperature(ix,iy) - temperature(ixm,iy)) / dxc(ixm)
        tg2 = (temperature(ix,iy) - temperature(ixp,iy)) / dxc(ix)
        ! Y temperature gradient at the x boundaries of the current cell
        ! Uses centred difference on averaged values, so likely very smoothed
        tg_a1 = (tb_p1 - tb_m1) / (dyc(iy) + dyc(iym))
        tg_a2 = (tb_p2 - tb_m2) / (dyc(iy) + dyc(iym))

        fc_sp1 = kappa_0 * tb1**pow / (modb1**2 + min_b) &
            * (bx(ixm,iy) * (tg1 * bx(ixm,iy) + tg_a1 * byf1) + tg1 * min_b)
        fc_sp2 = kappa_0 * tb2**pow / (modb2**2 + min_b) &
            * (bx(ix ,iy) * (tg2 * bx(ix ,iy) + tg_a2 * byf2) + tg2 * min_b)

        ! Saturated Conductive Flux
        rho_b1 = 0.5_num * (rho(ix,iy) + rho(ixm,iy))
        rho_b2 = 0.5_num * (rho(ix,iy) + rho(ixp,iy))
        fc_sa1 = 42.85_num * rho_b1 * tb1**1.5_num
        fc_sa2 = 42.85_num * rho_b2 * tb2**1.5_num  ! 42.85 = SRQT(m_p/m_e)

        ! Conductive Flux Limiter. Note flux_limiter is inverse of usual
        ! definition here
        fc1 = 1.0_num / (1.0_num / fc_sp1 + flux_limiter / fc_sa1)
        fc2 = 1.0_num / (1.0_num / fc_sp2 + flux_limiter / fc_sa2)

        flux(ix,iy) = (fc2 - fc1) / dxb(ix)

        ! Y flux
        bxf1 = 0.25_num * (bx(ix,iy) + bx(ixm,iy) + bx(ix,iyp) + bx(ixm,iyp))
        bxf2 = 0.25_num * (bx(ix,iy) + bx(ixm,iy) + bx(ix,iym) + bx(ixm,iym))
        bzf1 = 0.5_num * (bz(ix,iy) + bz(ix,iym))
        bzf2 = 0.5_num * (bz(ix,iy) + bz(ix,iyp))
        modb1 = SQRT(by(ix,iy )**2 + bxf1**2 + bzf1**2)
        modb2 = SQRT(by(ix,iym)**2 + bxf2**2 + bzf2**2)

        ! Braginskii Conductive Flux
        tb1 = 0.5_num * (temperature(ix,iy) + temperature(ix,iym))
        tb2 = 0.5_num * (temperature(ix,iy) + temperature(ix,iyp))

        ! Temperature at the y boundaries in the cell right		
        tb_p1 = 0.5_num * (temperature(ixp,iy) + temperature(ixp,iym))
        tb_p2 = 0.5_num * (temperature(ixp,iy) + temperature(ixp,iyp))
        ! Temperature at the y boundaries in the cell left		
        tb_m1 = 0.5_num * (temperature(ixm,iy) + temperature(ixm,iym))
        tb_m2 = 0.5_num * (temperature(ixm,iy) + temperature(ixm,iyp))
        ! Y temperature gradient at the y boundaries of the current cell		
        tg1 = (temperature(ix,iy) - temperature(ix,iym)) / dyc(iym)
        tg2 = (temperature(ix,iy) - temperature(ix,iyp)) / dyc(iy )
        ! X temperature gradient at the y boundaries of the current cell
        ! Uses centred difference on averaged values, so likely very smoothed
        tg_a1 = (tb_p1 - tb_m1) / (dxc(ix) + dxc(ixm))
        tg_a2 = (tb_p2 - tb_m2) / (dxc(ix) + dxc(ixm))

        fc_sp1 = kappa_0 * tb1**pow / (modb1**2 + min_b) &
            * (by(ix,iym) * (tg1 * by(ix,iym) + tg_a1 * bxf1) + min_b * tg1)
        fc_sp2 = kappa_0 * tb2**pow / (modb2**2 + min_b) &
            * (by(ix,iy ) * (tg2 * by(ix,iy ) + tg_a2 * bxf2) + min_b * tg2)

        ! Saturated Conductive Flux
        rho_b1 = 0.5_num * (rho(ix,iy) + rho(ix,iym))
        rho_b2 = 0.5_num * (rho(ix,iy) + rho(ix,iyp))
        fc_sa1 = 42.85_num * rho_b1 * tb1**1.5_num
        fc_sa2 = 42.85_num * rho_b2 * tb2**1.5_num  ! 42.85 = SRQT(m_p/m_e)

        ! Conductive Flux Limiter
        fc1 = 1.0_num / (1.0_num / fc_sp1 + flux_limiter / fc_sa1)
        fc2 = 1.0_num / (1.0_num / fc_sp2 + flux_limiter / fc_sa2)

        flux(ix,iy) = flux(ix,iy) + (fc2 - fc1) / dyb(iy)
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

    ! Superstepping based conduction code
    ! 2nd order Runge-Kutta-Lagrange (RKL2) scheme
    ! Based on Meyer at al. 2012 variables named as in that paper

    REAL(num), DIMENSION(:,:), ALLOCATABLE :: flux
    REAL(num) :: omega_1
    REAL(num), DIMENSION(0:n_s_stages) :: a, b
    REAL(num), DIMENSION(1:n_s_stages) :: mu_tilde
    REAL(num), DIMENSION(2:n_s_stages) :: mu, nu, gamma_tilde
    ! intermediate solutions Y=[Y_0, Y_j-2, Y_j-1, Y_j]
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: Y
    REAL(num) :: c0, c1, dj, fac
    REAL(num) :: Lc_Yj_1                  ! L^c(Y_j-1)
    REAL(num), DIMENSION(1:nx,1:ny) :: Lc_Y0    ! L^c(Y_0)
    INTEGER :: j

    ALLOCATE(flux(-1:nx+2,-1:ny+2))
    ALLOCATE(Y(0:3,-1:nx+2,-1:ny+2))
    flux = 0.0_num
    Y = 0.0_num

    dj = REAL(n_s_stages, num)
    omega_1 = 4.0_num / (dj**2 + dj - 2.0_num)
    b(0:2) = 1.0_num / 3.0_num
    DO j = 3, n_s_stages
      dj = REAL(j, num)
      b(j) = (dj**2 + dj - 2.0_num) / (2.0_num * dj * (dj + 1.0_num))
    END DO
    a = 1.0_num - b
    mu_tilde(1) = omega_1 / 3.0_num

    DO j = 2, n_s_stages
      dj = REAL(j, num)
      fac = (2.0_num * dj - 1.0_num) / dj
      mu_tilde(j) = fac * omega_1 * (b(j) / b(j-1))
      mu(j) = fac * (b(j) / b(j-1))
      nu(j) = (1.0_num - dj) / dj * (b(j) / b(j-2))
      gamma_tilde(j) = -1.0_num * a(j-1) * mu_tilde(j)
    END DO

    !!!!! RKL2 s stage scheme !!!!!
    ! initial condition
    Y(0,:,:) = energy
    !! boundary conditions
    DO j = 1, 3
      Y(j,:,:) = energy
    END DO


    !! First STS stage
    CALL heat_flux(Y(0,:,:) * 0.5_num * (gamma - 1.0_num), flux)

    DO iy = 1, ny
      DO ix = 1, nx
        ! Store L^c(Y_0) to use in stages 2-s
        Lc_Y0(ix,iy) = flux(ix,iy) / rho(ix,iy)
        c0 =  mu_tilde(1) * dt * Lc_Y0(ix,iy)
        Y(1,ix,iy) = Y(0,ix,iy) + c0
      END DO
    END DO

    Y(2,1:nx,1:ny) = Y(1,1:nx,1:ny)
    Y(1,1:nx,1:ny) = Y(0,1:nx,1:ny)

    DO j = 2, n_s_stages
      CALL heat_flux(Y(2,:,:) * 0.5_num * (gamma - 1.0_num), flux)

      DO iy = 1, ny
        DO ix = 1, nx
          c0 = gamma_tilde(j) * dt * Lc_Y0(ix,iy)
          Lc_Yj_1 = flux(ix,iy) / rho(ix,iy)
          c1 = mu_tilde(j) * dt * Lc_Yj_1
          ! Y_j
          Y(3,ix,iy) = mu(j) * Y(2,ix,iy) + nu(j) * Y(1,ix,iy) &
              + (1.0_num - mu(j) - nu(j)) * Y(0,ix,iy)
          Y(3,ix,iy) = Y(3,ix,iy) + c1 + c0
        END DO
      END DO

      IF (j < n_s_stages) THEN
        ! For jth stage  Y=[Y_0, Y_j-2, Y_j-1, Y_j]
        Y(1,:,:) = Y(2,:,:)
        ! This is not ideal, but it allows you to not have special boundary
        ! conditions
        energy = Y(3,:,:)
        CALL energy_bcs
        Y(2,:,:) = energy
        Y(3,1:nx,1:ny) = 0.0_num
      END IF
    END DO

    energy(1:nx,1:ny) = Y(3,1:nx,1:ny)

    CALL energy_bcs

  END SUBROUTINE heat_conduct_sts2



  SUBROUTINE rad_losses(density, temperature, xi, height, rad)

    ! Returns the normalised RTV losses

    REAL(num), INTENT(IN) :: density, temperature, xi, height
    REAL(num), INTENT(OUT) :: rad
    REAL(num), DIMENSION(7) :: trange = (/0.02_num, 0.0398_num, 0.0794_num, &
        0.251_num, 0.562_num, 1.995_num, 10.0_num/)
    REAL(num), DIMENSION(6) :: psi = (/1.2303_num, 870.96_num, 5.496_num, &
        0.3467_num, 1.0_num, 1.6218_num/)
    REAL(num), DIMENSION(6) :: alpha = (/0.0_num, 2.0_num, 0.0_num, &
        -2.0_num, 0.0_num, -2.0_num / 3.0_num/)
    REAL(num) :: tmk, factor
    INTEGER :: i

    rad = 0.0_num

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

  END SUBROUTINE rad_losses



  FUNCTION heating(density, t0)

    ! For a given density and temperature returns a user specific
    ! heating function

    REAL(num), INTENT(IN) :: density, t0
    REAL(num) :: heating
    REAL(num) :: tmk
    REAL(num) :: heat0 = 0.0_num
    REAL(num) :: rho_coronal = 0.0_num

    heating = 0.0_num
    IF (.NOT. coronal_heating) RETURN

    tmk = t0 * t2tmk
    ! For low density and high temeprature define a heating course term
    IF (density < rho_coronal .AND. tmk > 0.02_num) &
        heating = 100.0_num * heat0 * density**2

  END FUNCTION heating

END MODULE conduct
