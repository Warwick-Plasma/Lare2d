!******************************************************************************
! Controls all I/O and diagnostics. Output files are 'lare2d.dat',
! 'control.dat', 'en.dat' and a series of snapshots in 'nnnn.sdf'
!******************************************************************************

MODULE diagnostics

  USE shared_data
  USE boundary
  USE conduct
  USE output_cartesian
  USE output
  USE iocontrol
  USE version_data

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: set_dt, output_routines, energy_correction

CONTAINS

  !****************************************************************************
  ! Call the output routines
  !****************************************************************************

  SUBROUTINE output_routines(i) ! i = step index

    INTEGER, INTENT(IN) :: i

    INTEGER, PARAMETER :: out = 1000
    INTEGER, SAVE :: index = 1, step = 1
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: array
    LOGICAL :: print_arrays, last_call
    REAL(num), DIMENSION(c_ndims) :: stagger = 0.0_num
    INTEGER, DIMENSION(c_ndims) :: dims

    ! This output routine uses the same structure as needed for MPI output.
    ! This is more complicated than need for the serial code
    ! rank always equals zero in this serial code
    CHARACTER(LEN=9+data_dir_max_length+n_zeros) :: filename
    CHARACTER(LEN=35) :: filename_desc

    REAL(num) :: t_out = 0.0_num
    REAL(num) :: en_ke = 0.0_num, en_int = 0.0_num
    REAL(num) :: en_b = 0.0_num, heating_visc = 0.0_num
    REAL(num) :: heating_ohmic = 0.0_num
    REAL(num) :: total

    dims = (/nx_global+1, ny_global+1/)

    ! Make sure output fits arrays
    IF (nsteps >= out) step = nsteps / out + 1

    ! Done just once at the start
    IF (i == 0 .AND. rank == 0) THEN
      CALL output_log
      IF (.NOT. restart) WRITE(en_unit) num, 6
    END IF

    ! Do every (step) steps
    IF (MOD(i, step) == 0 .OR. last_call) THEN
      t_out = time
      CALL energy_account(en_b, en_ke, en_int)

      CALL MPI_ALLREDUCE(total_visc_heating, total, 1, mpireal, MPI_SUM, &
          comm, errcode)

      heating_visc = total

      CALL MPI_ALLREDUCE(total_ohmic_heating, total, 1, mpireal, MPI_SUM, &
          comm, errcode)

      heating_ohmic = total

      IF (rank == 0) THEN
        WRITE(en_unit) t_out, en_b, en_ke, en_int
        WRITE(en_unit) heating_visc, heating_ohmic
      END IF

      index = index + 1
    END IF

    ! Check if snapshot is needed
    CALL io_test(i, print_arrays, last_call)

    ! Output a snapshot of arrays
    IF (print_arrays) THEN
      IF (rank == 0) THEN
        WRITE(stat_unit,*) 'Dumping ', file_number, ' at time', time
        CALL FLUSH(stat_unit)
      END IF

      ! Set the filename
      WRITE(filename_desc, '("(''nfs:'', a, ''/'', i", i3.3, ".", i3.3, &
          & ", ''.cfd'')")') n_zeros, n_zeros
      WRITE(filename, filename_desc) TRIM(data_dir), file_number

      CALL cfd_open(filename, rank, comm, MPI_MODE_CREATE + MPI_MODE_WRONLY)
      CALL cfd_write_snapshot_data(REAL(time, dbl), i, 0)

      ALLOCATE(array(0:nx,0:ny))

      CALL cfd_write_2d_cartesian_grid('Grid', 'Grid', &
          xb_global(0:nx_global), yb_global(0:ny_global), 0)

      IF (dump_mask(1)) THEN
        array = rho(0:nx,0:ny)
        CALL cfd_write_2d_cartesian_variable_parallel('Rho', 'Fluid', &
            dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      IF (dump_mask(2)) THEN
        array = energy(0:nx,0:ny)
        CALL cfd_write_2d_cartesian_variable_parallel('Energy', 'Fluid', &
            dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      IF (dump_mask(3)) THEN
        array = vx(0:nx,0:ny)
        CALL cfd_write_2d_cartesian_variable_parallel('Vx', 'Velocity', &
            dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      IF (dump_mask(4)) THEN
        array = vy(0:nx,0:ny)
        CALL cfd_write_2d_cartesian_variable_parallel('Vy', 'Velocity', &
            dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      IF (dump_mask(5)) THEN
        array = vz(0:nx,0:ny)
        CALL cfd_write_2d_cartesian_variable_parallel('Vz', 'Velocity', &
            dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      IF (dump_mask(6)) THEN
        array = bx(0:nx,0:ny)
        CALL cfd_write_2d_cartesian_variable_parallel('Bx', 'Magnetic_Field', &
            dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      IF (dump_mask(7)) THEN
        array = by(0:nx,0:ny)
        CALL cfd_write_2d_cartesian_variable_parallel('By', 'Magnetic_Field', &
            dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      IF (dump_mask(8)) THEN
        array = bz(0:nx,0:ny)
        CALL cfd_write_2d_cartesian_variable_parallel('Bz', 'Magnetic_Field', &
            dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      IF (dump_mask(9)) THEN
        DO iy = 0, ny
          DO ix = 0, nx
            array(ix,iy) = (gamma - 1.0_num) / (2.0_num - xi_n(ix,iy)) &
                * (energy(ix,iy) - (1.0_num - xi_n(ix,iy)) * ionise_pot)
          END DO
        END DO
        CALL cfd_write_2d_cartesian_variable_parallel('Temperature', 'Fluid', &
            dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      IF (dump_mask(10)) THEN
        DO iy = 0, ny
          DO ix = 0, nx
            array(ix,iy) = (gamma - 1.0_num) * rho(ix,iy) &
                * (energy(ix,iy) - (1.0_num - xi_n(ix,iy)) * ionise_pot)
          END DO
        END DO
        CALL cfd_write_2d_cartesian_variable_parallel('Pressure', 'Fluid', &
            dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      IF (dump_mask(11)) THEN
        array = SQRT(gamma * (gamma - 1.0_num) * energy(1:nx,1:ny))
        CALL cfd_write_2d_cartesian_variable_parallel('cs', 'Fluid', &
            dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      IF (dump_mask(12)) THEN
        array = parallel_current(0:nx,0:ny)
        CALL cfd_write_2d_cartesian_variable_parallel('j_par', 'PIP', &
            dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      IF (dump_mask(13)) THEN
        array = perp_current(0:nx,0:ny)
        CALL cfd_write_2d_cartesian_variable_parallel('j_perp', 'PIP', &
            dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      IF (dump_mask(14)) THEN
        array = xi_n(0:nx,0:ny)
        CALL cfd_write_2d_cartesian_variable_parallel('neutral_fraction', &
            'PIP', dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      IF (dump_mask(15)) THEN
        array = eta_perp(0:nx,0:ny)
        CALL cfd_write_2d_cartesian_variable_parallel('eta_perp', 'PIP', &
            dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      IF (dump_mask(16)) THEN
        array = eta(0:nx,0:ny)
        CALL cfd_write_2d_cartesian_variable_parallel('eta', 'PIP', &
            dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      IF (dump_mask(17)) THEN
        array = jx_r(0:nx,0:ny)
        CALL cfd_write_2d_cartesian_variable_parallel('jx', 'current', &
            dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      IF (dump_mask(18)) THEN
        array = jy_r(0:nx,0:ny)
        CALL cfd_write_2d_cartesian_variable_parallel('jy', 'current', &
            dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      IF (dump_mask(19)) THEN
        array = jz_r(0:nx,0:ny)
        CALL cfd_write_2d_cartesian_variable_parallel('jz', 'current', &
            dims, stagger, 'Grid', 'Grid', array, subtype)
      END IF

      DEALLOCATE(array)

      ! Close the file
      CALL cfd_close

      file_number = file_number + 1
    END IF

    ! Output energy diagnostics etc
    IF (last_call .AND. rank == 0) THEN
      WRITE(stat_unit,*) 'final nsteps / time = ', i, time
    END IF

  END SUBROUTINE output_routines



  !****************************************************************************
  ! Test whether any of the conditions for doing output on the current
  ! iteration are met
  !****************************************************************************

  SUBROUTINE io_test(i, print_arrays, last_call)

    INTEGER, INTENT(IN) :: i
    LOGICAL, INTENT(OUT) :: print_arrays, last_call

    REAL(num), SAVE :: t1 = 0.0_num

    IF (restart) THEN
      t1 = time
      restart = .FALSE.
    END IF

    print_arrays = .FALSE.
    last_call = .FALSE.

    IF (time >= t1) THEN
      print_arrays = .TRUE.
      t1 = t1 + dt_snapshots
    END IF

    IF (time >= t_end .OR. i == nsteps) THEN
      last_call = .TRUE.
      print_arrays = .TRUE.
    END IF

  END SUBROUTINE io_test



  SUBROUTINE set_dt ! sets CFL limited step

    ! Assumes all variables are defined at the same point. Be careful with
    ! setting 'dt_multiplier' if you expect massive changes across cells.

    REAL(num) :: vxbm, vxbp, avxm, avxp, dvx, ax
    REAL(num) :: vybm, vybp, avym, avyp, dvy, ay
    REAL(num) :: cons, cs, area
    REAL(num) :: dxlocal, dt_local, dtr_local, dt1, dt2, dth_local, dt_rad
    REAL(num) :: dt_locals(3), dt_min(3)

    dt_local = largest_number
    dtr_local = largest_number
    dth_local = largest_number
    dt_rad = largest_number
    cons = gamma * (gamma - 1.0_num)

    DO iy = 0, ny
      iym = iy - 1
      DO ix = 0, nx
        ixm = ix - 1

        ! Fix dt for Lagrangian step
        w1 = bx(ix,iy)**2 + by(ix,iy)**2 + bz(ix,iy)**2
        ! Sound speed squared
        cs = cons * energy(ix,iy)

        w2 = SQRT(cs + w1 / MAX(rho(ix,iy), none_zero) &
            + 2.0_num * p_visc(ix,iy) / MAX(rho(ix,iy), none_zero))

        ! Find ideal MHD CFL limit for Lagrangian step
        dt1 = MIN(dxb(ix), dyb(iy)) / w2
        dt_local = MIN(dt_local, dt1)

        ! Now find dt for remap step
        ax = 0.5_num * dyb(iy)
        ay = 0.5_num * dxb(ix)
        vxbm = (vx(ixm,iy ) + vx(ixm,iym)) * ax
        vxbp = (vx(ix ,iy ) + vx(ix ,iym)) * ax
        vybm = (vy(ix ,iym) + vy(ixm,iym)) * ay
        vybp = (vy(ix ,iy ) + vy(ixm,iy )) * ay

        dvx = ABS(vxbp - vxbm)
        dvy = ABS(vybp - vybm)
        avxm = ABS(vxbm)
        avxp = ABS(vxbp)
        avym = ABS(vybm)
        avyp = ABS(vybp)

        area = dxb(ix) * dyb(iy)
        dt1 = area / MAX(avxm, avxp, dvx, 1.e-10_num * area)
        dt2 = area / MAX(avym, avyp, dvy, 1.e-10_num * area)

        ! Fix dt for remap step
        dt_local = MIN(dt_local, dt1, dt2, dt_rad)

        ! Note resistive limits assumes uniform resistivity hence cautious
        ! factor 0.2
        dxlocal = 1.0_num / (1.0_num / dxb(ix)**2 + 1.0_num / dyb(iy)**2)

        IF (cowling_resistivity) THEN
          dt1 = 0.2_num * dxlocal &
              / MAX(MAX(eta(ix,iy), eta_perp(ix,iy)), none_zero)
        ELSE
          dt1 = 0.2_num * dxlocal / MAX(eta(ix,iy), none_zero)
        END IF

        ! Adjust to accomodate resistive effects
        dtr_local = MIN(dtr_local, dt1)

        ! Hall MHD CFL limit
        dt1 = 0.75_num * rho(ix,iy) * MIN(dxb(ix), dyb(iy))**2 &
            / MAX(lambda_i(ix,iy) * SQRT(w1), none_zero)

        dth_local = MIN(dth_local, dt1)
      END DO
    END DO

    dt_locals(1) = dt_local
    dt_locals(2) = dtr_local
    dt_locals(3) = dth_local

    CALL MPI_ALLREDUCE(dt_locals, dt_min, 3, mpireal, MPI_MIN, comm, errcode)

    dt  = dt_multiplier * dt_min(1)
    dtr = dt_multiplier * dt_min(2)
    dth = dt_multiplier * dt_min(3)

    time = time + dt

  END SUBROUTINE set_dt



  SUBROUTINE energy_account(energy_b, energy_ke, energy_int)

    REAL(num), INTENT(OUT) :: energy_b, energy_ke, energy_int
    REAL(dbl) :: energy_b_local, energy_ke_local, energy_int_local
    REAL(dbl) :: energy_local(3), energy_sum(3)
    REAL(dbl) :: cv_v, rho_v, w1, w2, w3

    energy_b_local   = 0.0_dbl
    energy_ke_local  = 0.0_dbl
    energy_int_local = 0.0_dbl

    DO iy = 1, ny
      iym = iy - 1
      DO ix = 1, nx
        ixm = ix - 1

        w1 = (bx(ix,iy)**2 + bx(ixm,iy )**2) * 0.5_num
        w2 = (by(ix,iy)**2 + by(ix ,iym)**2) * 0.5_num
        w3 = bz(ix,iy)**2
        w1 = (w1 + w2 + w3) * 0.5_dbl
        energy_b_local = energy_b_local + w1 * cv(ix,iy)

        energy_int_local = energy_int_local &
            + energy(ix,iy) * rho(ix,iy) * cv(ix,iy)
      END DO
    END DO

    DO iy = 0, ny
      iyp = iy + 1
      DO ix = 0, nx
        ixp = ix + 1

        ! WARNING the KE is summed on the vertices
        rho_v = rho(ix,iy ) * cv(ix,iy ) + rho(ixp,iy ) * cv(ixp,iy ) &
              + rho(ix,iyp) * cv(ix,iyp) + rho(ixp,iyp) * cv(ixp,iyp)

        cv_v = cv(ix,iy) + cv(ixp,iy) + cv(ix,iyp) + cv(ixp,iyp)

        rho_v = rho_v / cv_v
        cv_v = cv_v * 0.25_dbl
        w1 = rho_v * cv_v * (vx(ix,iy)**2 + vy(ix,iy)**2 + vz(ix,iy)**2)

        IF (ix == 0 .OR. ix == nx) THEN
          w1 = w1 * 0.5_dbl
        END IF

        IF (iy == 0 .OR. iy == ny) THEN
          w1 = w1 * 0.5_dbl
        END IF

        energy_ke_local = energy_ke_local + w1 * 0.5_dbl
      END DO
    END DO

    energy_local(1) = energy_b_local
    energy_local(2) = energy_ke_local
    energy_local(3) = energy_int_local

    CALL MPI_ALLREDUCE(energy_local, energy_sum, 3, MPI_DOUBLE_PRECISION, &
        MPI_SUM, comm, errcode)

    energy_b   = REAL(energy_sum(1), num)
    energy_ke  = REAL(energy_sum(2), num)
    energy_int = REAL(energy_sum(3), num)

  END SUBROUTINE energy_account



  SUBROUTINE energy_correction

    delta_ke = -delta_ke
    WHERE (delta_ke < 0.0_num) delta_ke = 0.0_num
    delta_ke = delta_ke / (rho * cv)

    DO iy = 1, ny
      DO ix = 1, nx
        energy(ix,iy) = energy(ix,iy) + delta_ke(ix,iy)
      END DO
    END DO

    CALL energy_bcs

  END SUBROUTINE energy_correction



  SUBROUTINE output_log

    ! Writes basic data to 'lare2d.dat'

    WRITE(stat_unit,*) ascii_header
    WRITE(stat_unit,*)
    WRITE(stat_unit,*) 'nprocx, nprocy = ', nprocx, nprocy
    WRITE(stat_unit,*) 'nx, ny = ', nx, ny
    WRITE(stat_unit,*)
    WRITE(stat_unit,*) 'length_x = ', length_x
    WRITE(stat_unit,*) 'length_y = ', length_y
    WRITE(stat_unit,*)
#ifdef QMONO
    WRITE(stat_unit,*) 'q_mono viscosity'
#else
    WRITE(stat_unit,*) 'tensor shock viscosity'
#endif
    WRITE(stat_unit,*) 'linear viscosity coeff = ', visc1
    WRITE(stat_unit,*) 'quadratic viscosity coeff = ', visc2
    WRITE(stat_unit,*) 'uniform tensor viscosity coeff = ', visc3
    WRITE(stat_unit,*) 'j_max = ', j_max
    WRITE(stat_unit,*) 'eta0 = ', eta0
    WRITE(stat_unit,*) 'eta_background = ', eta_background
    WRITE(stat_unit,*) 'kappa = ', kappa_0
    WRITE(stat_unit,*)
    WRITE(stat_unit,*) 't_start, t_end = ', time, t_end
    WRITE(stat_unit,*) 'nsteps =', nsteps
    WRITE(stat_unit,*)

  END SUBROUTINE output_log

END MODULE diagnostics
