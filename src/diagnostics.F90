!******************************************************************************
! Controls all I/O and diagnostics. Output files are 'lare2d.dat',
! 'control.dat', 'en.dat' and a series of snapshots in 'nnnn.sdf'
!******************************************************************************

MODULE diagnostics

  USE shared_data
  USE boundary
  USE conduct
  USE sdf
  USE version_data

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: set_dt, output_routines, energy_correction, write_file, output_log

  REAL(dbl) :: visc_heating
  LOGICAL, SAVE :: visc_heating_updated = .FALSE.

CONTAINS

  !****************************************************************************
  ! Call the output routines
  !****************************************************************************

  SUBROUTINE output_routines(i) ! i = step index

    INTEGER, INTENT(IN) :: i
    INTEGER, PARAMETER :: outstep = 1
    LOGICAL :: print_arrays, last_call
    REAL(num) :: t_out = 0.0_num
    REAL(num) :: en_ke = 0.0_num, en_int = 0.0_num, en_b = 0.0_num
    REAL(dbl) :: var_local(en_nvars-1), var_sum(en_nvars-1)

#ifdef NO_IO
    RETURN
#endif

    visc_heating_updated = .FALSE.

    ! Do every (outstep) steps
    IF (MOD(i, outstep) == 0 .OR. last_call) THEN
      t_out = time
      CALL energy_account(en_b, en_ke, en_int, .FALSE.)

      var_local(1) = en_b
      var_local(2) = en_ke
      var_local(3) = en_int
      var_local(4) = total_visc_heating
      var_local(5) = total_ohmic_heating

      CALL MPI_ALLREDUCE(var_local, var_sum, en_nvars-1, MPI_DOUBLE_PRECISION, &
          MPI_SUM, comm, errcode)

      visc_heating = var_sum(4)
      visc_heating_updated = .TRUE.

      IF (rank == 0) THEN
        WRITE(en_unit) t_out, REAL(var_sum, num)
      END IF
    END IF

    ! Check if snapshot is needed
    CALL io_test(i, print_arrays, last_call)

    IF (print_arrays) CALL write_file(i)

    ! Output energy diagnostics etc
    IF (last_call .AND. rank == 0) THEN
      WRITE(stat_unit,*) 'final nsteps / time = ', i, time
    END IF

  END SUBROUTINE output_routines



  !****************************************************************************
  ! Write SDF file
  !****************************************************************************

  SUBROUTINE write_file(i) ! i = step index

    INTEGER, INTENT(IN) :: i
    REAL(num), DIMENSION(:), ALLOCATABLE :: work
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: array
    LOGICAL :: restart_flag, convert
    INTEGER, DIMENSION(c_ndims) :: global_dims, dims
    CHARACTER(LEN=22) :: filename_fmt
    CHARACTER(LEN=5+n_zeros+c_id_length) :: filename
    CHARACTER(LEN=6+data_dir_max_length+n_zeros+c_id_length) :: full_filename
    CHARACTER(LEN=c_id_length) :: varname, units
    TYPE(sdf_file_handle) :: sdf_handle
    LOGICAL, SAVE :: first = .TRUE.

    global_dims = (/ nx_global, ny_global /)

    IF (first) THEN
      ! Resize the {x,y}b_global to be the correct size for output
      ALLOCATE(work(-2:MAX(nx_global,ny_global)+2))

      work(-2:nx_global+2) = xb_global
      DEALLOCATE(xb_global)
      ALLOCATE(xb_global(0:nx_global))
      xb_global = work(0:nx_global)

      work(-2:ny_global+2) = yb_global
      DEALLOCATE(yb_global)
      ALLOCATE(yb_global(0:ny_global))
      yb_global = work(0:ny_global)

      DEALLOCATE(work)
      first = .FALSE.
    END IF

    ! Output a snapshot of arrays
    IF (rank == 0) THEN
      WRITE(stat_unit,*) 'Dumping ', file_number, ' at time', time
      CALL FLUSH(stat_unit)
    END IF

    ! Set the filename. Allows a maximum of 10^999 output dumps.
    WRITE(filename_fmt, '(''(a, i'', i3.3, ''.'', i3.3, '', ".sdf")'')') &
        n_zeros, n_zeros
    WRITE(filename, filename_fmt) TRIM(file_prefix), file_number
    full_filename = TRIM(filesystem) // TRIM(data_dir) // '/' // TRIM(filename)

    ! If dump_mask(1:8) are true then this file can be used for restarting
    restart_flag = ALL(dump_mask(1:8))

    convert = .FALSE.

    IF (.NOT.visc_heating_updated) THEN
      CALL MPI_ALLREDUCE(total_visc_heating, visc_heating, 1, &
          MPI_DOUBLE_PRECISION, MPI_SUM, comm, errcode)
      visc_heating_updated = .TRUE.
    END IF

    CALL sdf_open(sdf_handle, full_filename, comm, c_sdf_write)
    CALL sdf_write_header(sdf_handle, TRIM(c_code_name), 1, i, time, &
        restart_flag, jobid)
    CALL sdf_write_run_info(sdf_handle, c_version, c_revision, c_minor_rev, &
        c_commit_id, '', c_compile_machine, c_compile_flags, 0_8, &
        c_compile_date, run_date)
    CALL sdf_write_cpu_split(sdf_handle, 'cpu_rank', 'CPUs/Original rank', &
        cell_nx_maxs, cell_ny_maxs)
    CALL sdf_write_srl(sdf_handle, 'dt', 'Time increment', dt)
    CALL sdf_write_srl(sdf_handle, 'time_prev', 'Last dump time requested', &
        time_prev)
    CALL sdf_write_srl(sdf_handle, 'visc_heating', 'Viscous heating total', &
        visc_heating)

    CALL sdf_write_srl_plain_mesh(sdf_handle, 'grid', 'Grid/Grid', &
        xb_global, yb_global, convert)

    IF (dump_mask(1)) THEN
      varname = 'Rho'
      units = 'kg/m^3'
      dims = global_dims

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'Fluid/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_cell_centre, 'grid', rho, &
          cell_distribution, cell_subarray, convert)
    END IF

    IF (dump_mask(2)) THEN
      varname = 'Energy'
      units = 'J'
      dims = global_dims

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'Fluid/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_cell_centre, 'grid', energy, &
          cell_distribution, cell_subarray, convert)
    END IF

    IF (dump_mask(3)) THEN
      varname = 'Vx'
      units = 'm/s'
      dims = global_dims + 1

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'Velocity/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_vertex, 'grid', vx, &
          node_distribution, node_subarray, convert)
    END IF

    IF (dump_mask(4)) THEN
      varname = 'Vy'
      units = 'm/s'
      dims = global_dims + 1

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'Velocity/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_vertex, 'grid', vy, &
          node_distribution, node_subarray, convert)
    END IF

    IF (dump_mask(5)) THEN
      varname = 'Vz'
      units = 'm/s'
      dims = global_dims + 1

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'Velocity/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_vertex, 'grid', vz, &
          node_distribution, node_subarray, convert)
    END IF

    IF (dump_mask(6)) THEN
      varname = 'Bx'
      units = 'T'
      dims = global_dims
      dims(1) = dims(1) + 1

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'Magnetic_Field/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_bx, 'grid', bx, &
          bx_distribution, bx_subarray, convert)
    END IF

    IF (dump_mask(7)) THEN
      varname = 'By'
      units = 'T'
      dims = global_dims
      dims(2) = dims(2) + 1

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'Magnetic_Field/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_by, 'grid', by, &
          by_distribution, by_subarray, convert)
    END IF

    IF (dump_mask(8)) THEN
      varname = 'Bz'
      units = 'T'
      dims = global_dims

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'Magnetic_Field/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_bz, 'grid', bz, &
          bz_distribution, bz_subarray, convert)
    END IF

    IF (dump_mask(9)) THEN
      varname = 'Temperature'
      units = 'K'
      dims = global_dims

      IF (.NOT.ALLOCATED(array)) ALLOCATE(array(nx,ny))

      DO iy = 1, ny
        DO ix = 1, nx
          array(ix,iy) = (gamma - 1.0_num) / (2.0_num - xi_n(ix,iy)) &
              * (energy(ix,iy) - (1.0_num - xi_n(ix,iy)) * ionise_pot)
        END DO
      END DO

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'Fluid/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_cell_centre, 'grid', array, &
          cell_distribution, cellng_subarray, convert)
    END IF

    IF (dump_mask(10)) THEN
      varname = 'Pressure'
      units = 'Pa'
      dims = global_dims

      IF (.NOT.ALLOCATED(array)) ALLOCATE(array(nx,ny))

      DO iy = 1, ny
        DO ix = 1, nx
          array(ix,iy) = (gamma - 1.0_num) * rho(ix,iy) &
              * (energy(ix,iy) - (1.0_num - xi_n(ix,iy)) * ionise_pot)
        END DO
      END DO

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'Fluid/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_cell_centre, 'grid', array, &
          cell_distribution, cellng_subarray, convert)
    END IF

    IF (dump_mask(11)) THEN
      varname = 'Cs'
      units = 'm/s'
      dims = global_dims

      IF (.NOT.ALLOCATED(array)) ALLOCATE(array(nx,ny))

      array = SQRT(gamma * (gamma - 1.0_num) * energy(1:nx,1:ny))

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'Fluid/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_cell_centre, 'grid', array, &
          cell_distribution, cellng_subarray, convert)
    END IF

    IF (dump_mask(12)) THEN
      varname = 'j_par'
      units = 'A/m^2'
      dims = global_dims + 1

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'PIP/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_vertex, 'grid', parallel_current, &
          node_distribution, nodeng_subarray, convert)
    END IF

    IF (dump_mask(13)) THEN
      varname = 'j_perp'
      units = 'A/m^2'
      dims = global_dims + 1

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'PIP/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_vertex, 'grid', perp_current, &
          node_distribution, nodeng_subarray, convert)
    END IF

    IF (dump_mask(14)) THEN
      varname = 'neutral_fraction'
      units = '%'
      dims = global_dims + 1

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'PIP/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_vertex, 'grid', xi_n, &
          node_distribution, nodeng_subarray, convert)
    END IF

    IF (dump_mask(15)) THEN
      varname = 'eta_perp'
      units = ''
      dims = global_dims + 1

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'PIP/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_vertex, 'grid', eta_perp, &
          node_distribution, nodeng_subarray, convert)
    END IF

    IF (dump_mask(16)) THEN
      varname = 'eta'
      units = ''
      dims = global_dims + 1

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'PIP/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_vertex, 'grid', eta, &
          node_distribution, nodeng_subarray, convert)
    END IF

    IF (dump_mask(17)) THEN
      varname = 'Jx'
      units = 'A/m^2'
      dims = global_dims + 1

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'Current/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_vertex, 'grid', jx_r, &
          node_distribution, nodeng_subarray, convert)
    END IF

    IF (dump_mask(18)) THEN
      varname = 'Jy'
      units = 'A/m^2'
      dims = global_dims + 1

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'Current/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_vertex, 'grid', jy_r, &
          node_distribution, nodeng_subarray, convert)
    END IF

    IF (dump_mask(19)) THEN
      varname = 'Jz'
      units = 'A/m^2'
      dims = global_dims + 1

      CALL sdf_write_plain_variable(sdf_handle, TRIM(varname), &
          'Current/' // TRIM(varname), TRIM(units), dims, &
          c_stagger_vertex, 'grid', jz_r, &
          node_distribution, nodeng_subarray, convert)
    END IF

    IF (ALLOCATED(array)) DEALLOCATE(array)

    ! Close the file
    CALL sdf_close(sdf_handle)

    file_number = file_number + 1

  END SUBROUTINE write_file



  !****************************************************************************
  ! Test whether any of the conditions for doing output on the current
  ! iteration are met
  !****************************************************************************

  SUBROUTINE io_test(i, print_arrays, last_call)

    INTEGER, INTENT(IN) :: i
    LOGICAL, INTENT(OUT) :: print_arrays, last_call

    REAL(num), SAVE :: t1 = 0.0_num

    IF (restart) THEN
      t1 = time_prev + dt_snapshots
    END IF

    print_arrays = .FALSE.
    last_call = .FALSE.

    IF (time >= t1) THEN
      print_arrays = .TRUE.
      time_prev = t1
      t1 = t1 + dt_snapshots
    END IF

    IF (time >= t_end .OR. i == nsteps) THEN
      last_call = .TRUE.
      print_arrays = .TRUE.
    END IF

    IF (restart) THEN
      print_arrays = .FALSE.
      file_number = file_number + 1
      restart = .FALSE.
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
    LOGICAL, SAVE :: first = .TRUE.

    IF (first) THEN
      first = .FALSE.
      IF (restart) THEN
        dt = dt_from_restart
        RETURN
      END IF
    END IF

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



  SUBROUTINE energy_account(energy_b, energy_ke, energy_int, do_sum)

    REAL(dbl), INTENT(OUT) :: energy_b, energy_ke, energy_int
    LOGICAL, INTENT(IN) :: do_sum
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

    IF (do_sum) THEN
      energy_local(1) = energy_b_local
      energy_local(2) = energy_ke_local
      energy_local(3) = energy_int_local

      CALL MPI_ALLREDUCE(energy_local, energy_sum, 3, MPI_DOUBLE_PRECISION, &
          MPI_SUM, comm, errcode)

      energy_b   = energy_sum(1)
      energy_ke  = energy_sum(2)
      energy_int = energy_sum(3)
    ELSE
      energy_b   = energy_b_local
      energy_ke  = energy_ke_local
      energy_int = energy_int_local
    END IF

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

    IF (restart) THEN
      WRITE(stat_unit,*)
      WRITE(stat_unit,*) '#####################################################'
      WRITE(stat_unit,*)
      WRITE(stat_unit,*) 'Restarting from step ', step, ' and time ', time
      WRITE(stat_unit,*)
    END IF

    WRITE(stat_unit,*) ascii_header
    WRITE(stat_unit,*)
    WRITE(stat_unit,*) 'nprocx, nprocy = ', nprocx, nprocy
    WRITE(stat_unit,*) 'nx, ny = ', nx_global, ny_global
    WRITE(stat_unit,*)
    WRITE(stat_unit,*) 'length_x = ', length_x
    WRITE(stat_unit,*) 'length_y = ', length_y
    WRITE(stat_unit,*)
#ifdef QMONO
    WRITE(stat_unit,*) 'q_mono viscosity (-DQMONO)'
#else
    WRITE(stat_unit,*) 'tensor shock viscosity'
#endif
#ifdef FOURTHORDER
    WRITE(stat_unit,*) '4th-order resistive update (-DFOURTHORDER)'
#endif
#ifdef SINGLE
    WRITE(stat_unit,*) 'single precision (-DSINGLE)'
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
