MODULE conduct

  USE shared_data
  USE boundary
  USE neutral

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
        temp = gm1 * (2.0_num - xi_n(ix,iy)) &
            *(energy(ix,iy)-(1.0_num - xi_n(ix,iy))&
            *ionise_pot)
        temp = gm1 * rho(ix,iy) / (kappa_0 * temp**pow)
        dt1 = temp * dxb(ix)**2
        dt2 = temp * dyb(iy)**2
        dt_parab = MIN(dt_parab, dt1, dt2)
      END DO
    END DO

    dt_parab = dt_multiplier * dt_parab/SQRT(2.0_num)

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
    REAL(num) :: tb, tg, fc_sp, rho_b
    REAL(num) :: tg_a, tb_p, tb_m
    REAL(num) :: fc_sa, fc, modb, hfl
    REAL(num) :: byf, bxf, bzf

    hfl = 0.0_num
    IF (heat_flux_limiter) hfl = 1.0_num

    flux=0.0_num
    DO iy = 0, ny
      iyp = iy + 1
      iym = iy - 1
      DO ix = 0, nx
        ixp = ix + 1
        ixm = ix - 1

        ! X flux
        byf = 0.25_num * (by(ix,iy) + by(ix,iym) + by(ixp,iy) + by(ixp,iym))
        bzf = 0.5_num * (bz(ix,iy) + bz(ixp,iy))

        modb = SQRT(bx(ix ,iy)**2 + byf**2 + bzf**2)

        ! Braginskii Conductive Flux
        ! Temperature at the x boundaries in the current cell
        tb = 0.5_num * (temperature(ix,iy) + temperature(ixp,iy))

        ! Temperature at the x boundaries in the cell above
        tb_p = 0.5_num * (temperature(ix,iyp) + temperature(ixp,iyp))
        ! Temperature at the x boundaries in the cell below
        tb_m = 0.5_num * (temperature(ix,iym) + temperature(ixp,iym))
        ! X temperature gradient at the x boundaries of the current cell
        tg = (temperature(ixp,iy) - temperature(ix,iy)) / dxc(ix)
        ! Y temperature gradient at the x boundaries of the current cell
        ! Uses centred difference on averaged values, so likely very smoothed
        tg_a = (tb_p - tb_m) / (dyc(iy) + dyc(iym))

        fc_sp = - kappa_0 * tb**pow / (modb**2 + min_b) &
            * (bx(ix ,iy) * (tg * bx(ix ,iy) + tg_a * byf) + tg * min_b)

        ! Saturated Conductive Flux
        rho_b = 0.5_num * (rho(ix,iy) + rho(ixp,iy))
        fc_sa = 42.85_num * flux_limiter * rho_b * tb**1.5_num  ! 42.85 = SRQT(m_p/m_e)

        ! Conductive Flux Limiter. 
        fc = (1.0_num - hfl) * fc_sp + hfl * fc_sp * fc_sa / MAX(ABS(fc_sp) + fc_sa, none_zero)

        flux(ix,iy) = flux(ix,iy) - fc / dxb(ix)
        flux(ixp,iy) = flux(ixp,iy) + fc / dxb(ix)

        ! Y flux
        bxf = 0.25_num * (bx(ix,iy) + bx(ixm,iy) + bx(ix,iym) + bx(ixm,iym))
        bzf = 0.5_num * (bz(ix,iy) + bz(ix,iyp))
        modb = SQRT(by(ix,iym)**2 + bxf**2 + bzf**2)

        ! Braginskii Conductive Flux
        tb = 0.5_num * (temperature(ix,iy) + temperature(ix,iyp))

        ! Temperature at the y boundaries in the cell right
        tb_p = 0.5_num * (temperature(ixp,iy) + temperature(ixp,iyp))
        ! Temperature at the y boundaries in the cell left
        tb_m = 0.5_num * (temperature(ixm,iy) + temperature(ixm,iyp))
        ! Y temperature gradient at the y boundaries of the current cell
        tg = (temperature(ix,iyp) - temperature(ix,iy)) / dyc(iy)
        ! X temperature gradient at the y boundaries of the current cell
        ! Uses centred difference on averaged values, so likely very smoothed
        tg_a = (tb_p - tb_m) / (dxc(ix) + dxc(ixm))

        fc_sp = - kappa_0 * tb**pow / (modb**2 + min_b) &
            * (by(ix,iy ) * (tg * by(ix,iy ) + tg_a * bxf) + min_b * tg)

        ! Saturated Conductive Flux
        rho_b = 0.5_num * (rho(ix,iy) + rho(ix,iyp))
        fc_sa = 42.85_num * flux_limiter * rho_b * tb**1.5_num  ! 42.85 = SRQT(m_p/m_e)

        ! Conductive Flux Limiter
        fc = (1.0_num - hfl) * fc_sp + hfl * fc_sp * fc_sa / MAX(ABS(fc_sp) + fc_sa, none_zero)

        flux(ix,iy) = flux(ix,iy) - fc / dyb(iy)
        flux(ix,iyp) = flux(ix,iyp) + fc / dyb(iy)
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

    REAL(num), DIMENSION(:,:), ALLOCATABLE :: flux, Lc_Y0
    REAL(num) :: omega_1
    REAL(num), DIMENSION(0:n_s_stages) :: a, b
    REAL(num), DIMENSION(1:n_s_stages) :: mu_tilde
    REAL(num), DIMENSION(2:n_s_stages) :: mu, nu, gamma_tilde
    ! intermediate solutions Y=[Y_0, Y_j-2, Y_j-1, Y_j]
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: Y
    REAL(num) :: c0, c1, dj, fac
    REAL(num) :: Lc_Yj_1                  ! L^c(Y_j-1)
    INTEGER :: j

    ALLOCATE(flux(-1:nx+2,-1:ny+2))
    ALLOCATE(Lc_Y0(1:nx,1:ny))
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
    temperature(:,:) = (gamma-1.0_num) / (2.0_num - xi_n(:,:)) &
        *(Y(0,:,:)-(1.0_num - xi_n(:,:)) * ionise_pot)

    CALL temperature_bcs
    CALL heat_flux(temperature, flux)

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
      temperature(:,:) = (gamma-1.0_num) / (2.0_num - xi_n(:,:)) &
         * (Y(2,:,:)-(1.0_num - xi_n(:,:)) * ionise_pot)
		  
      CALL temperature_bcs
      CALL heat_flux(temperature, flux)

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
        energy(:,:) = Y(3,:,:)
        CALL energy_bcs
        IF (eos_number /= EOS_IDEAL) CALL neutral_fraction
        Y(2,:,:) = energy
        Y(3,1:nx,1:ny) = 0.0_num
      END IF
    END DO

    energy(1:nx,1:ny) = Y(3,1:nx,1:ny)
    CALL energy_bcs

    DEALLOCATE(flux)
    DEALLOCATE(Y)
    DEALLOCATE(Lc_Y0)

  END SUBROUTINE heat_conduct_sts2



END MODULE conduct
