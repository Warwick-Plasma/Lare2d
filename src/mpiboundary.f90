MODULE mpiboundary

  USE shared_data
  IMPLICIT NONE

CONTAINS

  SUBROUTINE bfield_MPI

    CALL MPI_SENDRECV(bx(1:2,:), 2*(ny+4), mpireal, left, tag,  &
         bx(nx+1:nx+2,:), 2*(ny+4), mpireal, right, tag,comm,status,&
         errcode)
    CALL MPI_SENDRECV(bx(nx-2:nx,:), 3*(ny+4), mpireal, right, tag,  &
         bx(-2:0,:), 3*(ny+4), mpireal, left, tag, comm, status, errcode)
    CALL MPI_SENDRECV(by(1:2,:), 2*(ny+5), mpireal, left, tag,  &
         by(nx+1:nx+2,:), 2*(ny+5), mpireal, right, tag,comm,status,&
         errcode)
    CALL MPI_SENDRECV(by(nx-1:nx,:), 2*(ny+5), mpireal, right, tag,  &
         by(-1:0,:), 2*(ny+5), mpireal, left, tag, comm, status, errcode)
    CALL MPI_SENDRECV(bz(1:2,:), 2*(ny+4), mpireal, left, tag,  &
         bz(nx+1:nx+2,:), 2*(ny+4), mpireal, right, tag,comm,status,&
         errcode)
    CALL MPI_SENDRECV(bz(nx-1:nx,:), 2*(ny+4), mpireal, right, tag,  &
         bz(-1:0,:), 2*(ny+4), mpireal, left, tag, comm, status, errcode)

    CALL MPI_SENDRECV(bx(:,ny-1:ny), 2*(nx+5), mpireal, up, tag,  &
         bx(:,-1:0), 2*(nx+5), mpireal, down, tag, comm, status, errcode)
    CALL MPI_SENDRECV(bx(:,1:2), 2*(nx+5), mpireal, down, tag,  &
         bx(:,ny+1:ny+2), 2*(nx+5), mpireal, up, tag, comm, status, errcode)
    CALL MPI_SENDRECV(by(:,ny-2:ny), 3*(nx+4), mpireal, up, tag,  &
         by(:,-2:0), 3*(nx+4), mpireal, down, tag, comm, status, errcode)
    CALL MPI_SENDRECV(by(:,1:2), 2*(nx+4), mpireal, down, tag,  &
         by(:,ny+1:ny+2), 2*(nx+4), mpireal, up, tag, comm, status, errcode)
    CALL MPI_SENDRECV(bz(:,ny-1:ny), 2*(nx+4), mpireal, up, tag,  &
         bz(:,-1:0), 2*(nx+4), mpireal, down, tag, comm, status, errcode)
    CALL MPI_SENDRECV(bz(:,1:2), 2*(nx+4), mpireal, down, tag,  &
         bz(:,ny+1:ny+2), 2*(nx+4), mpireal, up, tag, comm, status, errcode)

  END SUBROUTINE bfield_MPI



  SUBROUTINE bz_MPI

    CALL MPI_SENDRECV(bz(1:2,:), 2*(ny+4), mpireal, left, tag,  &
         bz(nx+1:nx+2,:), 2*(ny+4), mpireal, right,tag,comm,status, &
         errcode)
    CALL MPI_SENDRECV(bz(nx-1:nx,:), 2*(ny+4), mpireal, right, tag,  &
         bz(-1:0,:), 2*(ny+4), mpireal, left, tag, comm, status, errcode)

    CALL MPI_SENDRECV(bz(:,ny-1:ny), 2*(nx+4), mpireal, up, tag,  &
         bz(:,-1:0), 2*(nx+4), mpireal, down, tag, comm, status, errcode)
    CALL MPI_SENDRECV(bz(:,1:2), 2*(nx+4), mpireal, down, tag,  &
         bz(:,ny+1:ny+2), 2*(nx+4), mpireal, up, tag, comm, status, errcode)

  END SUBROUTINE bz_MPI



  SUBROUTINE energy_MPI

    CALL MPI_SENDRECV(energy(1:2,:), 2*(ny+4), mpireal, left, tag,  &
         energy(nx+1:nx+2,:), 2*(ny+4), mpireal, right, tag, comm, &
         status, &
         errcode)
    CALL MPI_SENDRECV(energy(nx-1:nx,:), 2*(ny+4), mpireal, right, tag, &
         energy(-1:0,:), 2*(ny+4), mpireal, left, tag, comm, status, &
         errcode)

    CALL MPI_SENDRECV(energy(:,ny-1:ny), 2*(nx+4), mpireal, up, tag,  &
         energy(:,-1:0), 2*(nx+4), mpireal, down, tag, comm, status, errcode)
    CALL MPI_SENDRECV(energy(:,1:2), 2*(nx+4), mpireal, down, tag,  &
         energy(:,ny+1:ny+2), 2*(nx+4), mpireal, up, tag, comm, status, &
         errcode)

  END SUBROUTINE energy_MPI



  SUBROUTINE velocity_MPI

    CALL MPI_SENDRECV(vx(1:2,:), 2*(ny+5), mpireal, left, tag,  &
         vx(nx+1:nx+2,:), 2*(ny+5), mpireal, right, tag,comm,status,&
         errcode)
    CALL MPI_SENDRECV(vx(nx-2:nx,:), 3*(ny+5), mpireal, right, tag,  &
         vx(-2:0,:), 3*(ny+5), mpireal, left, tag, comm, status, errcode)
    CALL MPI_SENDRECV(vy(1:2,:), 2*(ny+5), mpireal, left, tag,  &
         vy(nx+1:nx+2,:), 2*(ny+5), mpireal, right, tag,comm,status,&
         errcode)
    CALL MPI_SENDRECV(vy(nx-2:nx,:), 3*(ny+5), mpireal, right, tag,  &
         vy(-2:0,:), 3*(ny+5), mpireal, left, tag, comm, status, errcode)
    CALL MPI_SENDRECV(vz(1:2,:), 2*(ny+5), mpireal, left, tag,  &
         vz(nx+1:nx+2,:), 2*(ny+5),mpireal, right, tag,comm,status, &
         errcode)
    CALL MPI_SENDRECV(vz(nx-2:nx,:), 3*(ny+5), mpireal, right, tag,  &
         vz(-2:0,:), 3*(ny+5), mpireal, left, tag, comm, status, errcode)

    CALL MPI_SENDRECV(vx(:,ny-2:ny), 3*(nx+5), mpireal, up, tag,  &
         vx(:,-2:0), 3*(nx+5), mpireal, down, tag, comm, status, errcode)
    CALL MPI_SENDRECV(vx(:,1:2), 2*(nx+5), mpireal, down, tag,  &
         vx(:,ny+1:ny+2), 2*(nx+5), mpireal, up, tag, comm, status, errcode)
    CALL MPI_SENDRECV(vy(:,ny-2:ny), 3*(nx+5), mpireal, up, tag,  &
         vy(:,-2:0), 3*(nx+5), mpireal, down, tag, comm, status, errcode)
    CALL MPI_SENDRECV(vy(:,1:2), 2*(nx+5), mpireal, down, tag,  &
         vy(:,ny+1:ny+2), 2*(nx+5), mpireal, up, tag, comm, status, errcode)
    CALL MPI_SENDRECV(vz(:,ny-2:ny), 3*(nx+5), mpireal, up, tag,  &
         vz(:,-2:0), 3*(nx+5), mpireal, down, tag, comm, status, errcode)
    CALL MPI_SENDRECV(vz(:,1:2), 2*(nx+5), mpireal, down, tag,  &
         vz(:,ny+1:ny+2), 2*(nx+5), mpireal, up, tag, comm, status, errcode)

  END SUBROUTINE velocity_MPI



  SUBROUTINE remap_v_MPI
    CALL MPI_SENDRECV(vx1(1:2,:), 2*(ny+5), mpireal, left, tag,  &
         vx1(nx+1:nx+2,:),2*(ny+5),mpireal,right,tag,comm,status,errcode)
    CALL MPI_SENDRECV(vx1(nx-2:nx,:), 3*(ny+5), mpireal, right, tag,  &
         vx1(-2:0,:), 3*(ny+5), mpireal, left, tag, comm, status,errcode)
    CALL MPI_SENDRECV(vy1(1:2,:), 2*(ny+5), mpireal, left, tag,  &
         vy1(nx+1:nx+2,:),2*(ny+5),mpireal,right,tag,comm,status,errcode)
    CALL MPI_SENDRECV(vy1(nx-2:nx,:), 3*(ny+5), mpireal, right, tag,  &
         vy1(-2:0,:), 3*(ny+5), mpireal, left, tag, comm,status,errcode)
    CALL MPI_SENDRECV(vz1(1:2,:), 2*(ny+5), mpireal, left, tag,  &
         vz1(nx+1:nx+2,:),2*(ny+5),mpireal,right,tag,comm,status,errcode)
    CALL MPI_SENDRECV(vz1(nx-2:nx,:), 3*(ny+5), mpireal, right, tag,  &
         vz1(-2:0,:), 3*(ny+5), mpireal, left, tag, comm,status,errcode)

    CALL MPI_SENDRECV(vx1(:,ny-2:ny), 3*(nx+5), mpireal, up, tag,  &
         vx1(:,-2:0), 3*(nx+5), mpireal, down, tag, comm, status, errcode)
    CALL MPI_SENDRECV(vx1(:,1:2), 2*(nx+5), mpireal, down, tag,  &
         vx1(:,ny+1:ny+2), 2*(nx+5), mpireal, up, tag, comm, status, errcode)
    CALL MPI_SENDRECV(vy1(:,ny-2:ny), 3*(nx+5), mpireal, up, tag,  &
         vy1(:,-2:0), 3*(nx+5), mpireal, down, tag, comm, status, errcode)
    CALL MPI_SENDRECV(vy1(:,1:2), 2*(nx+5), mpireal, down, tag,  &
         vy1(:,ny+1:ny+2), 2*(nx+5), mpireal, up, tag, comm, status, errcode)
    CALL MPI_SENDRECV(vz1(:,ny-2:ny), 3*(nx+5), mpireal, up, tag,  &
         vz1(:,-2:0), 3*(nx+5), mpireal, down, tag, comm, status, errcode)
    CALL MPI_SENDRECV(vz1(:,1:2), 2*(nx+5), mpireal, down, tag,  &
         vz1(:,ny+1:ny+2), 2*(nx+5), mpireal, up, tag, comm, status, errcode)
  END SUBROUTINE remap_v_MPI



  SUBROUTINE density_MPI

    CALL MPI_SENDRECV(rho(1:2,:), 2*(ny+4), mpireal, left, tag,  &
         rho(nx+1:nx+2,:), 2*(ny+4), mpireal, right, tag, comm, status, &
         errcode)
    CALL MPI_SENDRECV(rho(nx-1:nx,:), 2*(ny+4), mpireal, right, tag,  &
         rho(-1:0,:), 2*(ny+4), mpireal, left, tag, comm,status,errcode)

    CALL MPI_SENDRECV(rho(:,ny-1:ny), 2*(nx+4), mpireal, up, tag,  &
         rho(:,-1:0), 2*(nx+4), mpireal, down, tag, comm, status, errcode)
    CALL MPI_SENDRECV(rho(:,1:2), 2*(nx+4), mpireal, down, tag,  &
         rho(:,ny+1:ny+2), 2*(nx+4), mpireal, up, tag, comm, status, errcode)

  END SUBROUTINE density_MPI

END MODULE
