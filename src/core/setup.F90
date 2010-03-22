MODULE setup

  USE shared_data
  USE normalise
  USE iocommon
  USE iocontrol
  USE input
  USE input_cartesian

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: before_control, after_control
  PUBLIC :: grid
  PUBLIC :: open_files, close_files, restart_data
  PUBLIC :: normalise_code, normalise_neutral

  REAL(num), DIMENSION(:), ALLOCATABLE :: dxnew, dynew

CONTAINS

  SUBROUTINE before_control
    ! Setup basic variables which have to have default values

    nprocx = 0
    nprocy = 0

    time = 0.0_num
    gamma = 5.0_num / 3.0_num

    IF (num .EQ. 4) mpireal = MPI_REAL

  END SUBROUTINE before_control



  SUBROUTINE after_control
    ! Setup arrays and other variables which can only be set after
    ! user input

    IF (IAND(initial, IC_RESTART) .EQ. 0) restart_snapshot = 0

    p_visc = 0.0_num
    eta = 0.0_num
    grav = 0.0_num

    rho = 0.0_num
    energy = 0.0_num
    bx = 0.0_num
    by = 0.0_num
    bz = 0.0_num
    vx = 0.0_num
    vy = 0.0_num
    vz = 0.0_num

  END SUBROUTINE after_control



  SUBROUTINE normalise_code

    ! Normalise physical parameters (eta, grav, visc3 etc.)
    CALL normalise_constants    ! setup.f90
    ! Normalise the grid
    CALL normalise_grid         ! setup.f90
    ! Normalise the actual initial conditions
    CALL normalise_eqm          ! setup.f90
    ! Normalise the constants needed for the
    ! partially ionised plasma routines
    IF (include_neutrals) CALL normalise_neutral ! setup.f90

  END SUBROUTINE normalise_code



  SUBROUTINE normalise_grid

    ! Normalise the grid, including control volumes

    xb = xb / L0
    xc = xc / L0
    xb_global = xb_global / L0
    dxc = dxc / L0
    dxb = dxb / L0

    yb = yb / L0
    yc = yc / L0
    yb_global = yb_global / L0
    dyc = dyc / L0
    dyb = dyb / L0

    cv = cv / L0**2

  END SUBROUTINE normalise_grid



  SUBROUTINE normalise_eqm

    ! Normalise the initial conditions

    rho = rho / rho0
    energy = energy / ENERGY0
    vx = vx / vel0
    vy = vy / vel0
    vz = vz / vel0
    bx = bx / B0
    by = by / B0
    bz = bz / B0

  END SUBROUTINE normalise_eqm



  SUBROUTINE normalise_constants
    
    ! Normalise gravity etc.
    
    grav = grav / grav0
    visc3 = visc3 / visc0
    eta0 = eta0 / res0
    eta_background = eta_background / res0
    
    ! the SI value for the constant in the conductivity
    ! assuming ln(Lambda) = 18.4
    kappa_0 = 1.0e-11_num
    kappa_0 = kappa_0 / kappa0
    
    time = time / T0
    t_end = t_end / T0
    dt_snapshots = dt_snapshots / T0
    lambda_i = lambda_i / L0
    
  END SUBROUTINE normalise_constants



  SUBROUTINE normalise_neutral 
  
    ! Normalise constants used in the calculation of properties
    ! Of partially ionised plasmas

    ! Normalised mass

    REAL(num) :: eta_bar_0

    pressure0 = B0**2 / mu0 ! Pressure 
    energy0 = B0**2 / (mu0 * rho0)   
    res0 = 1.0_num
    temp0 = mbar * pressure0 / (kb * rho0) ! Temperature in K
    
    ! Normalise tbar
    t_bar = t_bar / temp0

    ! Redefine rbar to incluse normalisation from rho and b/f
    r_bar = r_bar * rho0 / temp0**(3.0_num / 2.0_num)

    ! Normalise eta_bar
    eta_bar_0 = rho0**2 * SQRT(temp0) * res0 / B0**2
    eta_bar = eta_bar / eta_bar_0

    ! Finally normalise ion_mass and ionise_pot which are needed in the code
    ionise_pot = ionise_pot / (energy0 * mbar)

    tr = tr / temp0

  END SUBROUTINE normalise_neutral



  SUBROUTINE grid ! stretched and staggered grid

    REAL(num) :: dx, dy, xcstar, ycstar
    INTEGER :: ix, iy

    ALLOCATE (xb_global(-2:nx_global+2))
    ALLOCATE (yb_global(-2:ny_global+2))
    ALLOCATE (dxnew(-2:nx_global+2))
    ALLOCATE (dynew(-2:ny_global+2))

    ! initially assume uniform grid
    dx = 1.0_num / REAL(nx_global, num)
    dy = 1.0_num / REAL(ny_global, num)

    length_x = x_end - x_start
    length_y = y_end - y_start

    ! grid cell boundary for x coordinates
    ! set to - 0.5_num to have x = 0 in centre of domain
    xb_global(0) = 0.0_num
    DO ix = -2, nx_global+2
      xb_global(ix) = xb_global(0) + REAL(ix, num) * dx
    END DO
    xb_global = xb_global * length_x + x_start

    IF (x_stretch) CALL stretch_x ! stretch grid ?

    ! define position of ghost cells using sizes of adjacent cells
    IF (xbc_right == BC_PERIODIC) THEN
      xb_global(nx_global+1) = xb_global(nx_global) &
          + (xb_global(1) - xb_global(0))
      xb_global(nx_global+2) = xb_global(nx_global) &
          + (xb_global(2) - xb_global(0))
      xb_global(-1) = xb_global(0) &
          - (xb_global(nx) - xb_global(nx-1))
      xb_global(-2) = xb_global(0) &
          - (xb_global(nx) - xb_global(nx-2))
    ELSE
      xb_global(nx_global+1) = 2.0_num * xb_global(nx_global) &
          - xb_global(nx_global-1)
      xb_global(nx_global+2) = 2.0_num * xb_global(nx_global) &
          - xb_global(nx_global-2)
      xb_global(-1) = 2.0_num * xb_global(0) - xb_global(1)
      xb_global(-2) = 2.0_num * xb_global(0) - xb_global(2)
    END IF

    xb = xb_global(coordinates(2)*nx-2:coordinates(2)*nx+nx+2)

    DO ix = -1, nx+2
      ixm = ix - 1
      xc(ix) = 0.5_num * (xb(ixm) + xb(ix)) ! cell centre
    END DO

    DO ix = -1, nx+1
      ixp = ix + 1
      dxc(ix) = xc(ixp) - xc(ix) ! distance between centres
    END DO

    IF (coordinates(2) == nprocx - 1) THEN
      dxc(nx+2) = dxc(nx+1)
    ELSE
      xcstar = 0.5_num * (xb(nx+2) + xb_global(coordinates(2)*nx+nx+3))
      dxc(nx+2) = xcstar - xc(nx+2)
    END IF

    DO ix = -1, nx+2
      ixm = ix - 1
      dxb(ix) = xb(ix) - xb(ixm) ! cell width
    END DO

    yb_global(0) = 0.0_num ! repeat for y
    DO iy = -2, ny_global+2
      yb_global(iy) = yb_global(0) + REAL(iy, num) * dy
    END DO
    yb_global = yb_global * length_y + y_start

    IF (y_stretch) CALL stretch_y

    IF (ybc_up == BC_PERIODIC) THEN
      yb_global(ny_global+1) = yb_global(ny_global) &
          + (yb_global(1) - yb_global(0))
      yb_global(ny_global+2) = yb_global(ny_global) &
          + (yb_global(2) - yb_global(0))
      yb_global(-1) = yb_global(0) &
          - (yb_global(ny) - yb_global(ny-1))
      yb_global(-2) = yb_global(0) &
          - (yb_global(ny) - yb_global(ny-2))
    ELSE
      yb_global(ny_global+1) = 2.0_num * yb_global(ny_global) &
          - yb_global(ny_global-1)
      yb_global(ny_global+2) = 2.0_num * yb_global(ny_global) &
          - yb_global(ny_global-2)
      yb_global(-1) = 2.0_num * yb_global(0) - yb_global(1)
      yb_global(-2) = 2.0_num * yb_global(0) - yb_global(2)
    END IF

    yb = yb_global(coordinates(1)*ny-2:coordinates(1)*ny+ny+2)

    DO iy = -1, ny+2
      iym = iy - 1
      yc(iy) = 0.5_num * (yb(iym) + yb(iy))
    END DO

    DO iy = -1, ny+1
      iyp = iy + 1
      dyc(iy) = yc(iyp) - yc(iy)
    END DO

    IF (coordinates(1) == nprocy - 1) THEN
      dyc(ny+2) = dyc(ny+1)
    ELSE
      ycstar = 0.5_num * (yb(ny+2) + yb_global(coordinates(1)*ny+ny+3))
      dyc(ny+2) = ycstar - yc(ny+2)
    END IF

    DO iy = -1, ny+2
      iym = iy - 1
      dyb(iy) = yb(iy) - yb(iym)
    END DO

    DO ix = -1, nx+2
      DO iy = -1, ny+2
        cv(ix, iy) = dxb(ix) * dyb(iy) ! define the cell area
      END DO
    END DO

    DEALLOCATE(dxnew, dynew)

  END SUBROUTINE grid



  ! Subroutine stretches the grid in the x direction
  SUBROUTINE stretch_x ! replace with any stretching algorithm as needed

    REAL(num) :: width, dx, L, f, lx_new

    ! new total length
    lx_new = 200.0_num

    ! centre of tanh stretching in unstretched coordinates
    L = length_x / 1.5_num

    ! width of tanh stretching in unstretched coordinates
    width = length_x / 10.0_num

    f = (lx_new - length_x) / (length_x - L) / 2.0_num

    dx = length_x / REAL(nx_global, num)
    dxnew = dx + f * (1.0_num + TANH((ABS(xb_global) - L) / width)) * dx

!!$    DO ix = nx_global/2+1, nx_global+2
!!$      xb_global(ix) = xb_global(ix-1) + dxnew(ix)
!!$    END DO
!!$    DO ix = nx_global/2-1, -2, -1
!!$      xb_global(ix) = xb_global(ix+1) - dxnew(ix)
!!$    END DO

    DO ix = 1, nx_global+2
      xb_global(ix) = xb_global(ix-1) + dxnew(ix)
    END DO

    length_x = lx_new

  END SUBROUTINE stretch_x



  ! Subroutine stretches the domain in the y direction
  SUBROUTINE stretch_y ! stretch domain upwards only

    REAL(num) :: width, dy, L, f, ly_new

    ! new tolal length
    ly_new = 100.0_num

    ! centre of tanh stretching in unstretched coordinates
    L = length_y / 1.5_num

    ! width of tanh stretching in unstretched coordinates
    width = length_y / 10.0_num

    f = (ly_new - length_y) / (length_y - L) / 2.0_num

    dy = length_y / REAL(ny_global, num)
    dynew = dy + f * (1.0_num + TANH((ABS(yb_global) - L) / width)) * dy

!!$    DO iy = ny_global/2+1, ny_global+2
!!$      yb_global(iy) = yb_global(iy-1) + dynew(iy)
!!$    END DO
!!$    DO iy = ny_global/2-1, -2, -1
!!$      yb_global(iy) = yb_global(iy+1) - dynew(iy)
!!$    END DO

    DO iy = 1, ny_global+2
      yb_global(iy) = yb_global(iy-1) + dynew(iy)
    END DO

    length_y = ly_new

  END SUBROUTINE stretch_y



  ! Open the output diagnostic files
  SUBROUTINE open_files

    CHARACTER(LEN = 11+data_dir_max_length) :: file2
    CHARACTER(LEN = 7+data_dir_max_length) :: file3
    INTEGER :: ios

    IF (rank == 0) THEN
      WRITE(file2, '(a, "/lare2d.dat")') TRIM(data_dir)
      OPEN(unit = 20, STATUS = 'REPLACE', FILE = file2, iostat = ios)

      IF (ios .NE. 0) THEN
        PRINT *, "Unable to open file lare2d.dat for writing. This is ", &
                 "most commonly caused by the output directory not existing"
        PRINT *, " "
        PRINT *, " "
        CALL MPI_ABORT(comm, errcode)
      END IF

      WRITE(file3, '(a, "/en.dat")') TRIM(data_dir)
      OPEN(unit = 30, STATUS = 'REPLACE', FILE = file3, &
          FORM = "binary", iostat = ios)

      IF (ios .NE. 0) THEN
        PRINT *, "Unable to open file en.dat for writing. This is ", &
                 "most commonly caused by the output directory not existing"
        PRINT *, " "
        PRINT *, " "
        CALL MPI_ABORT(comm, errcode)
      END IF
    END IF

  END SUBROUTINE open_files



  ! Close the output diagnostic files
  SUBROUTINE close_files

    IF (rank == 0) THEN
      CLOSE(unit = 20)
      CLOSE(unit = 30)
    END IF

  END SUBROUTINE close_files



  ! Subroutine to perform string comparisons
  FUNCTION str_cmp(str_in, str_test)

    CHARACTER(*), INTENT(IN) :: str_in, str_test
    CHARACTER(30) :: str_trim
    LOGICAL :: str_cmp

    str_trim = TRIM(ADJUSTL(str_in))

    IF (LEN(str_test) .GT. LEN(str_in)) THEN
      str_cmp = .FALSE.
      RETURN
    END IF

    IF (str_trim(LEN(str_test)+1:LEN(str_test)+1) .NE. " ") THEN
      str_cmp = .FALSE.
      RETURN
    END IF

    str_cmp = str_trim(1:LEN(str_test)) == str_test

  END FUNCTION str_cmp



  ! Restart from previous output dumps
  SUBROUTINE restart_data

    CHARACTER(LEN = 20+data_dir_max_length) :: filename
    CHARACTER(LEN = 20) :: name, class, mesh_name, mesh_class
    INTEGER :: nblocks, type, nd, sof, snap
    INTEGER, DIMENSION(2) :: dims 
    REAL(dbl) :: time_d
    REAL(num), DIMENSION(2) :: extent
    REAL(num), DIMENSION(2) :: stagger
    REAL(num), DIMENSION(:, :), ALLOCATABLE :: data

    ! Create the filename for the last snapshot
#ifdef MHDCLUSTER
    WRITE(filename, '("nfs:", a, "/", i4.4, ".cfd")') &
        TRIM(data_dir), restart_snapshot
#else
    WRITE(filename, '(a, "/", i4.4, ".cfd")') &
        TRIM(data_dir), restart_snapshot
#endif

    output_file = restart_snapshot

    ALLOCATE(data(0:nx, 0:ny))
    CALL cfd_open(filename, rank, comm, MPI_MODE_RDONLY)
    ! Open the file
    nblocks = cfd_get_nblocks()

    DO ix = 1, nblocks
      CALL cfd_get_next_block_info_all(name, class, type)
      IF (rank == 0) PRINT *, ix, name, class, type

      IF (type == TYPE_SNAPSHOT) THEN
        CALL cfd_get_snapshot(time_d, snap)
        time = time_d
      END IF

      IF (type == TYPE_MESH) THEN
        ! Strangely, LARE doesn't actually read in the grid from a file
        ! This can be fixed, but for the moment, just go with the flow and
        ! Replicate the old behaviour

        CALL cfd_skip_block()
      ELSE IF (type == TYPE_MESH_VARIABLE) THEN
        CALL cfd_get_common_meshtype_metadata_all(type, nd, sof)

        IF (nd /= DIMENSION_2D) THEN
          IF (rank == 0) PRINT *, "Non 2D Dataset found in input file, ", &
              "ignoring and continuting."
          CALL cfd_skip_block()
          CYCLE
        END IF

        IF (type /= VAR_CARTESIAN) THEN
          IF (rank == 0) PRINT *, "Non - Cartesian variable block found ", &
              "in file, ignoring and continuing"
          CALL cfd_skip_block()
          CYCLE
        END IF

        ! We now have a valid variable, let's load it up
        ! First error trapping
        CALL cfd_get_nd_cartesian_variable_metadata_all(nd, dims, extent, &
            stagger, mesh_name, mesh_class)

        IF (dims(1) /= nx_global+1 .OR. dims(2) /= ny_global+1) THEN
          IF (rank == 0) PRINT *, "Size of grid represented by one more ", &
              "variables invalid. Continuing"
          CALL cfd_skip_block
          CYCLE
        END IF

        IF (sof /= num) THEN
          IF (rank == 0) PRINT *, "Precision of data does not match ", &
              "precision of code. Continuing."
          CALL cfd_skip_block
        END IF

        ! We're not interested in the other parameters, so if we're here,
        ! load up the data

        CALL cfd_get_2d_cartesian_variable_parallel(data, subtype)

        ! Now have the data, just copy it to correct place

        IF (str_cmp(name(1:3), "Rho")) THEN
          rho(0:nx, 0:ny) = data
        END IF

        IF (str_cmp(name(1:6), "Energy")) THEN
          energy(0:nx, 0:ny) = data
        END IF

        IF (str_cmp(name(1:2), "Vx")) THEN
          vx(0:nx, 0:ny) = data
        END IF

        IF (str_cmp(name(1:2), "Vy")) THEN
          vy(0:nx, 0:ny) = data
        END IF

        IF (str_cmp(name(1:2), "Vz")) THEN
          vz(0:nx, 0:ny) = data
        END IF

        IF (str_cmp(name(1:2), "Bx")) THEN
          bx(0:nx, 0:ny) = data
        END IF

        IF (str_cmp(name(1:2), "By")) THEN
          by(0:nx, 0:ny) = data
        END IF

        IF (str_cmp(name(1:2), "Bz")) THEN
          bz(0:nx, 0:ny) = data
        END IF

        ! Should be at end of block, but force the point anyway
        CALL cfd_skip_block()
      ELSE
        ! Unknown block, just skip it
        CALL cfd_skip_block()
      END IF
    END DO

    DEALLOCATE(data)

    CALL cfd_close()

    CALL MPI_BARRIER(comm, errcode)

  END SUBROUTINE restart_data

END MODULE setup
