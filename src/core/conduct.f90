MODULE conduct

  USE shared_data
  USE boundary

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: conduct_heat

  REAL(num), DIMENSION(:,:), ALLOCATABLE :: temperature

CONTAINS

  !****************************************************************************
  ! Subroutine implementing Braginskii parallel thermal conduction.
  ! Notation and algorithm in Appendix of Manual
  !****************************************************************************

  SUBROUTINE conduct_heat

    REAL(num), DIMENSION(:,:), ALLOCATABLE :: uxkx, uxky
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: uykx, uyky
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: energy0, limiter, temperature0
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: radiation, heat_in, alpha
    REAL(num) :: txb, tyb
    REAL(num) :: bxc, byc, bzc, bpx, bpy
    REAL(num) :: ux, uy
    REAL(num) :: pow = 5.0_num / 2.0_num
    REAL(num) :: a1, a2, a3, a4, error, errmax, rad, alf, limit
    REAL(num) :: w, residual, q_f, q_nl, q_sh, q_shx, q_shy
    REAL(num) :: duxdxi2, duxdxm2, duydyi2, duydym2
    INTEGER :: loop, redblack, x1, y1
    LOGICAL :: converged
    REAL(num), PARAMETER :: fractional_error = 1.0e-4_num
    REAL(num), PARAMETER :: b_min = 1.0e-3_num

    ALLOCATE(uxkx(-1:nx+1,-1:ny+1))
    ALLOCATE(uxky(-1:nx+1,-1:ny+1))
    ALLOCATE(uykx(-1:nx+1,-1:ny+1))
    ALLOCATE(uyky(-1:nx+1,-1:ny+1))
    ALLOCATE(energy0     (-1:nx+2,-1:ny+2))
    ALLOCATE(limiter     (-1:nx+2,-1:ny+2))
    ALLOCATE(temperature0(-1:nx+2,-1:ny+2))
    ALLOCATE(temperature (-1:nx+2,-1:ny+2))
    ALLOCATE(radiation   (-1:nx+2,-1:ny+2))
    ALLOCATE(heat_in     (-1:nx+2,-1:ny+2))
    ALLOCATE(alpha       (-1:nx+2,-1:ny+2))

    temperature = (gamma - 1.0_num) &
        * (energy - (1.0_num - xi_n) * ionise_pot) / (2.0_num - xi_n)

    ! A wasteful but simple way to turn off conduction.
    ! Mostly not needed so only needs to be improved if heating/radiation
    ! needed without conduction
    IF (.NOT. conduction) kappa_0 = 0.0_num

    DO iy = -1, ny + 1
      iym = iy - 1
      iyp = iy + 1
      DO ix = -1, nx + 1
        ixm = ix - 1
        ixp = ix + 1

        ! x face centred B field
        bxc = bx(ix,iy)
        byc = by(ix,iy) + by(ixp,iy) + by(ix,iym) + by(ixp,iym)
        bzc = bz(ix,iy) + bz(ixp,iy)
        byc = 0.25_num * byc
        bzc = 0.5_num  * bzc
        bpx = SQRT(bxc**2 + byc**2 + bzc**2 + b_min**2)
        bpx = MAX(bpx, none_zero)

        txb = 0.5_num * (temperature(ix,iy) + temperature(ixp,iy))

        ! Direction of magnetic field on x face
        ux = bxc / bpx
        uy = byc / bpx

        ! Kappa along magnetic field, now a vector
        uxkx(ix,iy) = ux * ux * kappa_0 * txb**pow
        uxky(ix,iy) = ux * uy * kappa_0 * txb**pow

        ! Add symmetic conduction near b=0 points
        uxkx(ix,iy) = uxkx(ix,iy) &
            + b_min**2 * kappa_0 * txb**pow / (bpx**2 + b_min**2)

        ! y face centred B field
        bxc = bx(ix,iy) + bx(ix,iyp) + bx(ixm,iy) + bx(ixm,iyp)
        byc = by(ix,iy)
        bzc = bz(ix,iy) + bz(ix,iyp)
        bxc = 0.25_num * bxc
        bzc = 0.5_num  * bzc
        bpy = SQRT(bxc**2 + byc**2 + bzc**2 + b_min**2)
        bpy = MAX(bpy, none_zero)

        tyb = 0.5_num * (temperature(ix,iy) + temperature(ix,iyp))

        ! Direction of magnetic field on y face
        ux = bxc / bpy
        uy = byc / bpy

        ! Kappa along magnetic field, now a vector
        uykx(ix,iy) = uy * ux * kappa_0 * tyb**pow
        uyky(ix,iy) = uy * uy * kappa_0 * tyb**pow

        ! Add symmetic conduction near b=0 points
        uyky(ix,iy) = uyky(ix,iy) &
            + b_min**2 * kappa_0 * tyb**pow / (bpy**2 + b_min**2)
      END DO
    END DO

    IF (heat_flux_limiter) THEN
      DO iy = 0, ny + 1
        iym = iy - 1
        iyp = iy + 1
        DO ix = 0, nx + 1
          ixm = ix - 1
          ixp = ix + 1

          ! Estimate the parallel heat flux and the centre of a cell
          q_shx = (uxkx(ix,iy) + uxkx(ixm,iy)) / dxc(ix) &
              *   (temperature(ixp,iy) - temperature(ixm,iy)) &
              +   (uxky(ix,iy) + uxky(ix,iym)) / dyc(iy) &
              *   (temperature(ix,iyp) - temperature(ix,iym))

          q_shy = (uykx(ix,iy) + uykx(ixm,iy)) / dxc(ix) &
              *   (temperature(ixp,iy) - temperature(ixm,iy)) &
              +   (uyky(ix,iy) + uyky(ix,iym)) / dyc(iy) &
              *   (temperature(ix,iyp) - temperature(ix,iym))

          q_sh = SQRT(q_shx**2 + q_shy**2) / 16.0_num

          ! Estimate the free streaming limit
          ! 42.85 = SRQT(m_p/m_e)
          q_f = 42.85_num * flux_limiter * rho(ix,iy) &
              * MIN(temperature(ix,iy), temperature_100mk)**1.5_num
          q_nl = 1.0_num / (1.0_num / MAX(q_sh, none_zero) &
              +  1.0_num / MAX(q_f, none_zero))
          limiter(ix,iy) = q_nl / MAX(q_sh, none_zero) / 2.0_num
        END DO
      END DO

      DO iy = 0, ny + 1
        DO ix = 0, nx + 1
          limit = limiter(ix,iy) + limiter(ix+1,iy)
          uxkx(ix,iy) = uxkx(ix,iy) * limit
          uxky(ix,iy) = uxky(ix,iy) * limit

          limit = limiter(ix,iy) + limiter(ix,iy+1)
          uykx(ix,iy) = uykx(ix,iy) * limit
          uyky(ix,iy) = uyky(ix,iy) * limit
        END DO
      END DO
    END IF

    converged = .FALSE.
    w = 1.6_num       ! Initial over-relaxation parameter
    ! Store energy^{n}
    energy0 = energy
    temperature0 = temperature

    radiation = 0.0_num
    heat_in = 0.0_num
    alpha = 0.0_num
    DO iy = 1, ny
      DO ix = 1, nx
        heat_in(ix,iy) = heating(rho(ix,iy), temperature0(ix,iy))
        CALL rad_losses(rho(ix,iy), temperature0(ix,iy), &
                        xi_n(ix,iy), yc(iy), rad, alf)
        alpha(ix,iy) = alf
        radiation(ix,iy) = rad
      END DO
    END DO

    ! Iterate to get energy^{n+1} by SOR Guass-Seidel
    iterate: DO loop = 1, 100
      errmax = 0.0_num
      error = 0.0_num
      y1 = 1
      DO redblack = 1, 2
        x1 = y1
        DO iy = 1, ny
          iym = iy - 1
          iyp = iy + 1
          DO ix = x1, nx, 2
            ixm = ix - 1
            ixp = ix + 1

            duxdxi2 = uxkx(ix ,iy ) / dxc(ix ) / dxb(ix)
            duxdxm2 = uxkx(ixm,iy ) / dxc(ixm) / dxb(ix)
            duydyi2 = uyky(ix ,iy ) / dyc(iy ) / dyb(iy)
            duydym2 = uyky(ix ,iym) / dyc(iym) / dyb(iy)

            ! Terms containing energy(ix,iy) resulting from
            ! d^2/dx^2 and d^2/dy^2 derivatives
            a1 = duxdxi2 + duxdxm2 + duydyi2 + duydym2

            ! Terms not containing temperature(ix,iy) resulting from
            ! d^2/dx^2 and d^2/dy^2 derivatives
            a2 =  duxdxi2 * temperature(ixp,iy ) &
                + duxdxm2 * temperature(ixm,iy ) &
                + duydyi2 * temperature(ix ,iyp) &
                + duydym2 * temperature(ix ,iym)

            ! Terms not containing temperature(ix,iy) resulting from
            ! d^2/dxdy cross derivatives
            a2 = a2 + uxky(ix ,iy ) &
                * (temperature(ixp,iyp) + temperature(ix ,iyp) &
                -  temperature(ixp,iym) - temperature(ix ,iym)) &
                / (2.0_num * dxb(ix) * (dyc(iy) + dyc(iym)))
            a2 = a2 - uxky(ixm,iy ) &
                * (temperature(ix ,iyp) + temperature(ixm,iyp) &
                -  temperature(ix ,iym) - temperature(ixm,iym)) &
                / (2.0_num * dxb(ix) * (dyc(iy) + dyc(iym)))

            ! Terms not containing temperature(ix,iy) resulting from
            ! d^2/dydx cross derivatives
            a2 = a2 + uykx(ix ,iy ) &
                * (temperature(ixp,iyp) + temperature(ixp,iy ) &
                -  temperature(ixm,iyp) - temperature(ixm,iy )) &
                / (2.0_num * dyb(iy) * (dxc(ix) + dxc(ixm)))
            a2 = a2 - uykx(ix ,iym) &
                * (temperature(ixp,iy ) + temperature(ixp,iym) &
                -  temperature(ixm,iy ) - temperature(ixm,iym)) &
                / (2.0_num * dyb(iy) * (dxc(ix) + dxc(ixm)))

            a3 = (a1 + radiation(ix,iy) * alpha(ix,iy) &
                / temperature0(ix,iy)) * temperature(ix,iy) &
                / energy(ix,iy)

            a4 = a2 + heat_in(ix,iy) &
                - (1.0_num - alpha(ix,iy)) * radiation(ix,iy)

            a3 = a3 * dt / rho(ix,iy)
            a4 = a4 * dt / rho(ix,iy)

            residual = energy(ix,iy) &
                - (energy0(ix,iy) + a4) / (1.0_num + a3)
            energy(ix,iy) = MAX(energy(ix,iy) - w * residual, &
                (1.0_num - xi_n(ix,iy)) * ionise_pot)
            error = ABS(residual) / energy0(ix,iy)
            errmax = MAX(errmax, error)
          END DO
          x1 = 3 - x1
        END DO
        y1 = 3 - y1

        CALL energy_bcs

        temperature = (gamma - 1.0_num) &
            * (energy - (1.0_num - xi_n) * ionise_pot) / (2.0_num - xi_n)
      END DO

      CALL MPI_ALLREDUCE(errmax, error, 1, mpireal, MPI_MAX, comm, errcode)
      errmax = error

      IF (errmax .LT. fractional_error) THEN
        converged = .TRUE.
        EXIT iterate
      END IF
    END DO iterate

    IF (rank == 0 .AND. .NOT. converged) &
        PRINT*, 'Conduction failed at t = ', time

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
        CALL rad_losses(rho(1,ny), temperature(1,ny), &
                        xi_n(1,ny), yc(ny), rad, alf)
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
