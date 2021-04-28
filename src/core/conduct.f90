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
  
MODULE conduct

  USE shared_data
  USE boundary
  USE neutral

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: conduct_heat, calc_s_stages

  REAL(num), PARAMETER :: pow = 2.5_num
  REAL(num), PARAMETER :: min_b = 1.0e-6_num

CONTAINS

  !****************************************************************************
  ! Subroutine implementing Braginskii parallel thermal conduction.
  ! Notation and algorithm in Appendix of Manual
  !****************************************************************************

  SUBROUTINE conduct_heat

    ALLOCATE(larsen_factor(0:nx,0:ny))
    CALL calc_s_stages(.FALSE.)
    CALL heat_conduct_sts2
    DEALLOCATE(larsen_factor)

  END SUBROUTINE conduct_heat


  !****************************************************************************
  ! Subroutine calculating the number of stages needed for the RKL2
  !****************************************************************************

  SUBROUTINE calc_s_stages(lagrangian_call)

    LOGICAL, INTENT(IN) :: lagrangian_call
    REAL(num) :: stages, dt_parab, dt1, dt2
    REAL(num) :: q_fs, q_fs2, q_spx, q_spy, q_sp2
    REAL(num) :: temp, kappa1
    INTEGER :: n_s_stages_local

    ! Make sure arrays are allocated if calling this routine just to determine
    ! dt from the lagrangian setp
    IF (lagrangian_call) ALLOCATE(larsen_factor(0:nx,0:ny))

    DO iy=-1,ny+2
      DO ix=-1,nx+2
        temperature(ix,iy) = (gamma - 1.0_num) / (2.0_num - xi_n(ix,iy)) &
            * (energy(ix,iy) - (1.0_num - xi_n(ix,iy)) * ionise_pot)
      ENDDO
    ENDDO

    ! Include flux limiting through a larson factor correction to the
    ! conductivity
      DO iy = 0, ny
        iym = iy - 1
        iyp = iy + 1
        DO ix = 0, nx
          ixm = ix - 1
          ixp = ix + 1
          temp = temperature(ix,iy)**pow
          q_fs = flux_limiter * 42.85_num * rho(ix,iy) &  ! 42.85 = SQRT(m_i/m_e)
              * temperature(ix,iy)**1.5_num
          q_fs2 = q_fs**2
          q_spx = - kappa_0 * temp &
              * (temperature(ixp,iy) - temperature(ixm,iy)) &
              * 0.5_num / dxb(ix)
          q_spy = - kappa_0 * temp &
              * (temperature(ix,iyp) - temperature(ix,iym)) &
              * 0.5_num / dyb(iy)
          q_sp2 = q_spx**2 + q_spy**2
          larsen_factor(ix,iy) = q_fs / SQRT(q_fs2 + q_sp2)
        END DO
      END DO

    IF (.NOT. heat_flux_limiter) larsen_factor = 1.0_num

    dt_parab = 1.e10_num

    DO iy = 1, ny
      DO ix = 1, nx
        ! Estimate explicit thermal conduction time-step
        temp = (gamma - 1.0_num) / (2.0_num - xi_n(ix,iy)) &
             * (energy(ix,iy) - (1.0_num - xi_n(ix,iy)) * ionise_pot)
        kappa1 = kappa_0 * larsen_factor(ix,iy)
        temp = gm1 * rho(ix,iy) / (kappa1 * temp**pow)
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

    IF (lagrangian_call) DEALLOCATE(larsen_factor)

  END SUBROUTINE calc_s_stages



  !****************************************************************************
  ! Subroutine to calculate the heat flux
  !****************************************************************************

  SUBROUTINE heat_flux(temperature, flux)

    REAL(num), INTENT(IN), DIMENSION(-1:,-1:) :: temperature
    REAL(num), INTENT(OUT), DIMENSION(-1:,-1:) :: flux
    INTEGER :: ix, ixp, ixm
    REAL(num) :: tb, tg, fc_sp
    REAL(num) :: tg_a, tb_p, tb_m
    REAL(num) :: modb
    REAL(num) :: byf, bxf, bzf

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

        modb = bx(ix ,iy)**2 + byf**2 + bzf**2 + min_b

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

        fc_sp = - larsen_factor(ix,iy) * kappa_0 * tb**pow / modb &
            * (bx(ix ,iy) * (tg * bx(ix ,iy) + tg_a * byf) + tg * min_b)

        flux(ix,iy) = flux(ix,iy) - fc_sp / dxb(ix)
        flux(ixp,iy) = flux(ixp,iy) + fc_sp / dxb(ix)

        ! Y flux
        bxf = 0.25_num * (bx(ix,iy) + bx(ixm,iy) + bx(ix,iyp) + bx(ixm,iyp))
        bzf = 0.5_num * (bz(ix,iy) + bz(ix,iyp))
        modb = by(ix,iy)**2 + bxf**2 + bzf**2 + min_b

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

        fc_sp = - larsen_factor(ix,iy) * kappa_0 * tb**pow / modb &
            * (by(ix,iy ) * (tg * by(ix,iy ) + tg_a * bxf) + min_b * tg)

        flux(ix,iy) = flux(ix,iy) - fc_sp / dyb(iy)
        flux(ix,iyp) = flux(ix,iyp) + fc_sp / dyb(iy)
      END DO
    END DO

  END SUBROUTINE heat_flux




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
