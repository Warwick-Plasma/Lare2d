MODULE mpiboundary

  USE shared_data
  IMPLICIT NONE

CONTAINS

  SUBROUTINE bfield_MPI

    CALL MPI_SENDRECV(bx(1:2, :), 2*(ny+4), mpireal, proc_x_min, tag, &
        bx(nx+1:nx+2, :), 2*(ny+4), mpireal, proc_x_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(bx(nx-2:nx, :), 3*(ny+4), mpireal, proc_x_max, tag, &
        bx(-2:0, :), 3*(ny+4), mpireal, proc_x_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(by(1:2, :), 2*(ny+5), mpireal, proc_x_min, tag, &
        by(nx+1:nx+2, :), 2*(ny+5), mpireal, proc_x_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(by(nx-1:nx, :), 2*(ny+5), mpireal, proc_x_max, tag, &
        by(-1:0, :), 2*(ny+5), mpireal, proc_x_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(bz(1:2, :), 2*(ny+4), mpireal, proc_x_min, tag, &
        bz(nx+1:nx+2, :), 2*(ny+4), mpireal, proc_x_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(bz(nx-1:nx, :), 2*(ny+4), mpireal, proc_x_max, tag, &
        bz(-1:0, :), 2*(ny+4), mpireal, proc_x_min, tag, comm, &
        status, errcode)

    CALL MPI_SENDRECV(bx(:, ny-1:ny), 2*(nx+5), mpireal, proc_y_max, tag, &
        bx(:, -1:0), 2*(nx+5), mpireal, proc_y_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(bx(:, 1:2), 2*(nx+5), mpireal, proc_y_min, tag, &
        bx(:, ny+1:ny+2), 2*(nx+5), mpireal, proc_y_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(by(:, ny-2:ny), 3*(nx+4), mpireal, proc_y_max, tag, &
        by(:, -2:0), 3*(nx+4), mpireal, proc_y_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(by(:, 1:2), 2*(nx+4), mpireal, proc_y_min, tag, &
        by(:, ny+1:ny+2), 2*(nx+4), mpireal, proc_y_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(bz(:, ny-1:ny), 2*(nx+4), mpireal, proc_y_max, tag, &
        bz(:, -1:0), 2*(nx+4), mpireal, proc_y_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(bz(:, 1:2), 2*(nx+4), mpireal, proc_y_min, tag, &
        bz(:, ny+1:ny+2), 2*(nx+4), mpireal, proc_y_max, tag, comm, &
        status, errcode)

  END SUBROUTINE bfield_MPI



  SUBROUTINE bz_MPI

    CALL MPI_SENDRECV(bz(1:2, :), 2*(ny+4), mpireal, proc_x_min, tag, &
        bz(nx+1:nx+2, :), 2*(ny+4), mpireal, proc_x_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(bz(nx-1:nx, :), 2*(ny+4), mpireal, proc_x_max, tag, &
        bz(-1:0, :), 2*(ny+4), mpireal, proc_x_min, tag, comm, &
        status, errcode)

    CALL MPI_SENDRECV(bz(:, ny-1:ny), 2*(nx+4), mpireal, proc_y_max, tag, &
        bz(:, -1:0), 2*(nx+4), mpireal, proc_y_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(bz(:, 1:2), 2*(nx+4), mpireal, proc_y_min, tag, &
        bz(:, ny+1:ny+2), 2*(nx+4), mpireal, proc_y_max, tag, comm, &
        status, errcode)

  END SUBROUTINE bz_MPI



  SUBROUTINE energy_MPI

    CALL MPI_SENDRECV(energy(1:2, :), 2*(ny+4), mpireal, proc_x_min, tag, &
        energy(nx+1:nx+2, :), 2*(ny+4), mpireal, proc_x_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(energy(nx-1:nx, :), 2*(ny+4), mpireal, proc_x_max, tag, &
        energy(-1:0, :), 2*(ny+4), mpireal, proc_x_min, tag, comm, &
        status, errcode)

    CALL MPI_SENDRECV(energy(:, ny-1:ny), 2*(nx+4), mpireal, proc_y_max, tag, &
        energy(:, -1:0), 2*(nx+4), mpireal, proc_y_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(energy(:, 1:2), 2*(nx+4), mpireal, proc_y_min, tag, &
        energy(:, ny+1:ny+2), 2*(nx+4), mpireal, proc_y_max, tag, comm, &
        status, errcode)

  END SUBROUTINE energy_MPI



  SUBROUTINE velocity_MPI

    CALL MPI_SENDRECV(vx(1:2, :), 2*(ny+5), mpireal, proc_x_min, tag, &
        vx(nx+1:nx+2, :), 2*(ny+5), mpireal, proc_x_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vx(nx-2:nx, :), 3*(ny+5), mpireal, proc_x_max, tag, &
        vx(-2:0, :), 3*(ny+5), mpireal, proc_x_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vy(1:2, :), 2*(ny+5), mpireal, proc_x_min, tag, &
        vy(nx+1:nx+2, :), 2*(ny+5), mpireal, proc_x_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vy(nx-2:nx, :), 3*(ny+5), mpireal, proc_x_max, tag, &
        vy(-2:0, :), 3*(ny+5), mpireal, proc_x_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vz(1:2, :), 2*(ny+5), mpireal, proc_x_min, tag, &
        vz(nx+1:nx+2, :), 2*(ny+5), mpireal, proc_x_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vz(nx-2:nx, :), 3*(ny+5), mpireal, proc_x_max, tag, &
        vz(-2:0, :), 3*(ny+5), mpireal, proc_x_min, tag, comm, &
        status, errcode)

    CALL MPI_SENDRECV(vx(:, ny-2:ny), 3*(nx+5), mpireal, proc_y_max, tag, &
        vx(:, -2:0), 3*(nx+5), mpireal, proc_y_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vx(:, 1:2), 2*(nx+5), mpireal, proc_y_min, tag, &
        vx(:, ny+1:ny+2), 2*(nx+5), mpireal, proc_y_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vy(:, ny-2:ny), 3*(nx+5), mpireal, proc_y_max, tag, &
        vy(:, -2:0), 3*(nx+5), mpireal, proc_y_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vy(:, 1:2), 2*(nx+5), mpireal, proc_y_min, tag, &
        vy(:, ny+1:ny+2), 2*(nx+5), mpireal, proc_y_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vz(:, ny-2:ny), 3*(nx+5), mpireal, proc_y_max, tag, &
        vz(:, -2:0), 3*(nx+5), mpireal, proc_y_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vz(:, 1:2), 2*(nx+5), mpireal, proc_y_min, tag, &
        vz(:, ny+1:ny+2), 2*(nx+5), mpireal, proc_y_max, tag, comm, &
        status, errcode)

  END SUBROUTINE velocity_MPI



  SUBROUTINE remap_v_MPI

    CALL MPI_SENDRECV(vx1(1:2, :), 2*(ny+5), mpireal, proc_x_min, tag, &
        vx1(nx+1:nx+2, :), 2*(ny+5), mpireal, proc_x_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vx1(nx-2:nx, :), 3*(ny+5), mpireal, proc_x_max, tag, &
        vx1(-2:0, :), 3*(ny+5), mpireal, proc_x_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vy1(1:2, :), 2*(ny+5), mpireal, proc_x_min, tag, &
        vy1(nx+1:nx+2, :), 2*(ny+5), mpireal, proc_x_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vy1(nx-2:nx, :), 3*(ny+5), mpireal, proc_x_max, tag, &
        vy1(-2:0, :), 3*(ny+5), mpireal, proc_x_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vz1(1:2, :), 2*(ny+5), mpireal, proc_x_min, tag, &
        vz1(nx+1:nx+2, :), 2*(ny+5), mpireal, proc_x_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vz1(nx-2:nx, :), 3*(ny+5), mpireal, proc_x_max, tag, &
        vz1(-2:0, :), 3*(ny+5), mpireal, proc_x_min, tag, comm, &
        status, errcode)

    CALL MPI_SENDRECV(vx1(:, ny-2:ny), 3*(nx+5), mpireal, proc_y_max, tag, &
        vx1(:, -2:0), 3*(nx+5), mpireal, proc_y_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vx1(:, 1:2), 2*(nx+5), mpireal, proc_y_min, tag, &
        vx1(:, ny+1:ny+2), 2*(nx+5), mpireal, proc_y_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vy1(:, ny-2:ny), 3*(nx+5), mpireal, proc_y_max, tag, &
        vy1(:, -2:0), 3*(nx+5), mpireal, proc_y_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vy1(:, 1:2), 2*(nx+5), mpireal, proc_y_min, tag, &
        vy1(:, ny+1:ny+2), 2*(nx+5), mpireal, proc_y_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vz1(:, ny-2:ny), 3*(nx+5), mpireal, proc_y_max, tag, &
        vz1(:, -2:0), 3*(nx+5), mpireal, proc_y_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(vz1(:, 1:2), 2*(nx+5), mpireal, proc_y_min, tag, &
        vz1(:, ny+1:ny+2), 2*(nx+5), mpireal, proc_y_max, tag, comm, &
        status, errcode)

  END SUBROUTINE remap_v_MPI



  SUBROUTINE density_MPI

    CALL MPI_SENDRECV(rho(1:2, :), 2*(ny+4), mpireal, proc_x_min, tag, &
        rho(nx+1:nx+2, :), 2*(ny+4), mpireal, proc_x_max, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(rho(nx-1:nx, :), 2*(ny+4), mpireal, proc_x_max, tag, &
        rho(-1:0, :), 2*(ny+4), mpireal, proc_x_min, tag, comm, &
        status, errcode)

    CALL MPI_SENDRECV(rho(:, ny-1:ny), 2*(nx+4), mpireal, proc_y_max, tag, &
        rho(:, -1:0), 2*(nx+4), mpireal, proc_y_min, tag, comm, &
        status, errcode)
    CALL MPI_SENDRECV(rho(:, 1:2), 2*(nx+4), mpireal, proc_y_min, tag, &
        rho(:, ny+1:ny+2), 2*(nx+4), mpireal, proc_y_max, tag, comm, &
        status, errcode)

  END SUBROUTINE density_MPI

END MODULE mpiboundary
