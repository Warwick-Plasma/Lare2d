!*******************************************************************
! All the real ghost cell values are controlled by these routines.
! To speed things up it may be worth having this routine hard coded
! for each particular run, i.e. remove all if statements.
!*******************************************************************
MODULE boundary

  USE shared_data
  USE mpiboundary
  USE openboundary
  IMPLICIT NONE

CONTAINS


  SUBROUTINE set_boundary_conditions

    REAL(num) :: a
    LOGICAL :: first_call =.TRUE.

    IF (first_call) THEN
       any_open = .FALSE.
       IF ((xbc_right == open) .OR. (xbc_left == open)   &
            .OR. (ybc_up == open) .OR. (ybc_down == open)) any_open = .TRUE.
       first_call = .FALSE.
    ELSE
       ! when bzone=0 uses first order remap scheme so that farfield not used in remap
       bzone = 1.0_num

       IF (xbc_right == other .AND. right == MPI_PROC_NULL) bzone(nx-4:nx+2,:) = 0.0_num 
       IF (ybc_up == other .AND. up == MPI_PROC_NULL) bzone(:,ny-4:ny+2) = 0.0_num

       IF (xbc_right == open .AND. right == MPI_PROC_NULL) bzone(nx-4:nx+2,:) = 0.0_num 
       IF (xbc_left == open .AND. left == MPI_PROC_NULL) bzone(-1:4,:) = 0.0_num
       IF (ybc_up == open .AND. up == MPI_PROC_NULL) bzone(:,ny-4:ny+2) = 0.0_num
       IF (ybc_down == open .AND. down == MPI_PROC_NULL) bzone(:,-1:4) = 0.0_num

    END IF

  END SUBROUTINE set_boundary_conditions



  SUBROUTINE boundary_conditions

    REAL(num) :: a, d

    CALL bfield_bcs
    CALL energy_bcs
    CALL density_bcs
    CALL velocity_bcs

    CALL damp_boundaries

  END SUBROUTINE boundary_conditions




  SUBROUTINE bfield_bcs

    CALL bfield_MPI

    IF (right == MPI_PROC_NULL .AND. xbc_right == other) THEN
       bx(nx+1,:) = bx(nx-1,:)
       bx(nx+2,:) = bx(nx-2,:)
       by(nx+1,:) = by(nx,:)
       by(nx+2,:) = by(nx-1,:)
       bz(nx+1,:) = bz(nx,:)
       bz(nx+2,:) = bz(nx-1,:)
    END IF
    IF (left == MPI_PROC_NULL .AND. xbc_left == other) THEN
       bx(-1,:) = bx(1,:)
       bx(-2,:) = bx(2,:)
       by(0,:)  = by(1,:)
       by(-1,:) = by(2,:)
       bz(0,:)  = bz(1,:)
       bz(-1,:) = bz(2,:)
    END IF

    IF (up == MPI_PROC_NULL .AND. ybc_up == other) THEN
       bx(:,ny+1) = bx(:,ny)
       bx(:,ny+2) = bx(:,ny-1)
       by(:,ny+1) = by(:,ny-1)
       by(:,ny+2) = by(:,ny-2)
       bz(:,ny+1) = bz(:,ny)
       bz(:,ny+2) = bz(:,ny-1)
    END IF
    IF (down == MPI_PROC_NULL .AND. ybc_down == other) THEN
       bx(:,0) = bx(:,1)
       bx(:,-1) = bx(:,2)
       by(:,-1) = by(:,1)
       by(:,-2) = by(:,2)
       bz(:,0) = bz(:,1)
       bz(:,-1) = bz(:,2)
    END IF

  END SUBROUTINE bfield_bcs



  SUBROUTINE bz_bcs

    CALL bz_MPI

    IF (right == MPI_PROC_NULL .AND. xbc_right == other) THEN
       bz(nx+1,:) = bz(nx,:)
       bz(nx+2,:) = bz(nx-1,:)
    END IF
    IF (left == MPI_PROC_NULL .AND. xbc_left == other) THEN
       bz(0,:) = bz(1,:)
       bz(-1,:) = bz(2,:)
    END IF

    IF (up == MPI_PROC_NULL .AND. ybc_up == other) THEN
       bz(:,ny+1) = bz(:,ny)
       bz(:,ny+2) = bz(:,ny-1)
    END IF
    IF (down == MPI_PROC_NULL .AND. ybc_down == other) THEN
       bz(:,0) = bz(:,1)
       bz(:,-1) = bz(:,2)
    END IF

  END SUBROUTINE bz_bcs



  SUBROUTINE energy_bcs

    CALL energy_MPI

    IF (right == MPI_PROC_NULL .AND. xbc_right == other) THEN
       energy(nx+1,:) = energy(nx,:)
       energy(nx+2,:) = energy(nx-1,:)
    END IF
    IF (left == MPI_PROC_NULL .AND. xbc_left == other) THEN
       energy(0,:) = energy(1,:)
       energy(-1,:) = energy(2,:)
    END IF

    IF (up == MPI_PROC_NULL .AND. ybc_up == other) THEN
       energy(:,ny+1) = energy(:,ny)
       energy(:,ny+2) = energy(:,ny-1)
    END IF
    IF (down == MPI_PROC_NULL .AND. ybc_down == other) THEN
       energy(:,0) = energy(:,1)
       energy(:,-1) = energy(:,2)
    END IF

  END SUBROUTINE energy_bcs



  SUBROUTINE velocity_bcs

    CALL velocity_MPI

    IF (right == MPI_PROC_NULL .AND. xbc_right == other) THEN
       vx(nx,:)   = 0.0_num 
       vx(nx+1,:) = 0.0_num
       vx(nx+2,:) = 0.0_num
       vy(nx+1,:) = 0.0_num
       vy(nx+2,:) = 0.0_num
       vz(nx+1,:) = 0.0_num
       vz(nx+2,:) = 0.0_num
    END IF
    IF (left == MPI_PROC_NULL .AND. xbc_left == other) THEN
       vx(-2,:) = 0.0_num
       vx(-1,:) = 0.0_num
       vy(-1,:) = 0.0_num
       vy(-2,:) = 0.0_num
       vz(-1,:) = 0.0_num
       vz(-2,:) = 0.0_num
    END IF

    IF (up == MPI_PROC_NULL .AND. ybc_up == other) THEN
       vx(:,ny+1) = 0.0_num
       vx(:,ny+2) = 0.0_num
       vy(:,ny+1) = 0.0_num
       vy(:,ny+2) = 0.0_num
       vz(:,ny+1) = 0.0_num
       vz(:,ny+2) = 0.0_num
    END IF
    IF (down == MPI_PROC_NULL .AND. ybc_down == other) THEN
       vx(:,-1) = 0.0_num
       vx(:,-2) = 0.0_num
       vy(:,-1) = 0.0_num
       vy(:,-2) = 0.0_num
       vz(:,-1) = 0.0_num
       vz(:,-2) = 0.0_num
    END IF

  END SUBROUTINE velocity_bcs



  SUBROUTINE remap_v_bcs

    REAL(num) :: theta

    CALL remap_v_MPI

    IF (right == MPI_PROC_NULL) THEN
       IF (xbc_right == other) THEN
          vx1(nx,:)   = 0.0_num 
          vx1(nx+1,:) = 0.0_num
          vx1(nx+2,:) = 0.0_num
          vy1(nx+1,:) = 0.0_num
          vy1(nx+2,:) = 0.0_num
          vz1(nx+1,:) = 0.0_num
          vz1(nx+2,:) = 0.0_num
       ELSE IF (xbc_right == open) THEN
          vx1(nx+1, 1:ny) = 0.0_num
          vx1(nx+2, 1:ny) = 0.0_num
          vy1(nx+1, 1:ny) = 0.0_num
          vy1(nx+2, 1:ny) = 0.0_num
          vz1(nx+1, 1:ny) = 0.0_num
          vz1(nx+2, 1:ny) = 0.0_num
       END IF
    END IF
    IF (left == MPI_PROC_NULL) THEN
       IF (xbc_left == other) THEN
          vx1(-2,:) = 0.0_num
          vx1(-1,:) = 0.0_num
          vy1(-1,:) = 0.0_num
          vy1(-2,:) = 0.0_num
          vz1(-1,:) = 0.0_num
          vz1(-2,:) = 0.0_num
       ELSE IF (xbc_left == open) THEN
          vx1(-2,:) = 0.0_num
          vy1(-2,:) = 0.0_num
          vz1(-2,:) = 0.0_num
          vx1(-1,:) = 0.0_num
          vy1(-1,:) = 0.0_num
          vz1(-1,:) = 0.0_num
       END IF
    END IF

    IF (up == MPI_PROC_NULL) THEN
       IF (ybc_up == other) THEN
          vx1(:,ny+1) = 0.0_num
          vx1(:,ny+2) = 0.0_num
          vy1(:,ny+1) = 0.0_num
          vy1(:,ny+2) = 0.0_num
          vz1(:,ny+1) = 0.0_num
          vz1(:,ny+2) = 0.0_num
       ELSE IF (ybc_up == open) THEN
          vx1(:,ny+1) = 0.0_num
          vx1(:,ny+2) = 0.0_num
          vy1(:,ny+1) = 0.0_num
          vy1(:,ny+2) = 0.0_num
          vz1(:,ny+1) = 0.0_num
          vz1(:,ny+2) = 0.0_num
       END IF
    END IF
    IF (down == MPI_PROC_NULL) THEN
       IF (ybc_down == other) THEN
          vx1(:,-1) = 0.0_num
          vx1(:,-2) = 0.0_num
          vy1(:,-1) = 0.0_num
          vy1(:,-2) = 0.0_num
          vz1(:,-1) = 0.0_num
          vz1(:,-2) = 0.0_num
       ELSE IF (ybc_down == open) THEN
          vx1(:,-2) = 0.0_num
          vy1(:,-2) = 0.0_num
          vz1(:,-2) = 0.0_num
          vx1(:,-1) = 0.0_num
          vy1(:,-1) = 0.0_num
          vz1(:,-1) = 0.0_num
       END IF
    END IF

  END SUBROUTINE remap_v_bcs



  SUBROUTINE density_bcs

    CALL density_MPI

    IF (right == MPI_PROC_NULL .AND. xbc_right == other) THEN
       rho(nx+1,:) = rho(nx,:)
       rho(nx+2,:) = rho(nx-1,:)
    END IF
    IF (left == MPI_PROC_NULL .AND. xbc_left == other) THEN
       rho(0,:) = rho(1,:)
       rho(-1,:) = rho(2,:)
    END IF

    IF (up == MPI_PROC_NULL .AND. ybc_up == other) THEN
       rho(:,ny+1) = rho(:,ny)
       rho(:,ny+2) = rho(:,ny-1)
    END IF
    IF (down == MPI_PROC_NULL .AND. ybc_down == other) THEN
       rho(:,0) = rho(:,1)
       rho(:,-1) = rho(:,2)
    END IF

  END SUBROUTINE density_bcs




END MODULE boundary
