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
       IF ((xbc_right == BC_OPEN) .OR. (xbc_left == BC_OPEN)   &
            .OR. (ybc_up == BC_OPEN) .OR. (ybc_down == BC_OPEN)) any_open = .TRUE.
       first_call = .FALSE.
    ELSE
       ! when bzone=0 uses first order remap scheme so that farfield not used in remap
       bzone = 1.0_num

       IF (xbc_right == BC_OPEN .AND. right == MPI_PROC_NULL) bzone(nx-4:nx+2,:) = 0.0_num 
       IF (xbc_left == BC_OPEN .AND. left == MPI_PROC_NULL) bzone(-1:4,:) = 0.0_num
       IF (ybc_up == BC_OPEN .AND. up == MPI_PROC_NULL) bzone(:,ny-4:ny+2) = 0.0_num
       IF (ybc_down == BC_OPEN .AND. down == MPI_PROC_NULL) bzone(:,-1:4) = 0.0_num

    END IF

  END SUBROUTINE set_boundary_conditions



  SUBROUTINE boundary_conditions

    REAL(num) :: a, d

    CALL bfield_bcs
    CALL energy_bcs
    CALL density_bcs
    CALL velocity_bcs

    CALL damp_boundaries    ! in openboundary.f90

  END SUBROUTINE boundary_conditions



  SUBROUTINE bfield_bcs

    CALL bfield_MPI

    IF (right == MPI_PROC_NULL .AND. xbc_right == BC_OTHER) THEN
       bx(nx+1,:) = bx(nx-1,:)
       bx(nx+2,:) = bx(nx-2,:)
       by(nx+1,:) = by(nx,:)
       by(nx+2,:) = by(nx-1,:)
       bz(nx+1,:) = bz(nx,:)
       bz(nx+2,:) = bz(nx-1,:)
    END IF
    IF (left == MPI_PROC_NULL .AND. xbc_left == BC_OTHER) THEN
       bx(-1,:) = bx(1,:)
       bx(-2,:) = bx(2,:)
       by(0,:)  = by(1,:)
       by(-1,:) = by(2,:)
       bz(0,:)  = bz(1,:)
       bz(-1,:) = bz(2,:)
    END IF

    IF (up == MPI_PROC_NULL .AND. ybc_up == BC_OTHER) THEN
       bx(:,ny+1) = bx(:,ny)
       bx(:,ny+2) = bx(:,ny-1)
       by(:,ny+1) = by(:,ny-1)
       by(:,ny+2) = by(:,ny-2)
       bz(:,ny+1) = bz(:,ny)
       bz(:,ny+2) = bz(:,ny-1)
    END IF
    IF (down == MPI_PROC_NULL .AND. ybc_down == BC_OTHER) THEN
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

    IF (right == MPI_PROC_NULL .AND. xbc_right == BC_OTHER) THEN
       bz(nx+1,:) = bz(nx,:)
       bz(nx+2,:) = bz(nx-1,:)
    END IF
    IF (left == MPI_PROC_NULL .AND. xbc_left == BC_OTHER) THEN
       bz(0,:) = bz(1,:)
       bz(-1,:) = bz(2,:)
    END IF


    IF (up == MPI_PROC_NULL .AND. ybc_up == BC_OTHER) THEN
       bz(:,ny+1) = bz(:,ny)
       bz(:,ny+2) = bz(:,ny-1)
    END IF
    IF (down == MPI_PROC_NULL .AND. ybc_down == BC_OTHER) THEN
       bz(:,0) = bz(:,1)
       bz(:,-1) = bz(:,2)
    END IF

  END SUBROUTINE bz_bcs



  SUBROUTINE energy_bcs

    CALL energy_MPI

    IF (right == MPI_PROC_NULL .AND. xbc_right == BC_OTHER) THEN
       energy(nx+1,:) = energy(nx,:)
       energy(nx+2,:) = energy(nx-1,:)
    END IF
    IF (left == MPI_PROC_NULL .AND. xbc_left == BC_OTHER) THEN
       energy(0,:) = energy(1,:)
       energy(-1,:) = energy(2,:)
    END IF

    IF (up == MPI_PROC_NULL .AND. ybc_up == BC_OTHER) THEN
       energy(:,ny+1) = energy(:,ny)
       energy(:,ny+2) = energy(:,ny-1)
    END IF
    IF (down == MPI_PROC_NULL .AND. ybc_down == BC_OTHER) THEN
       energy(:,0) = energy(:,1)
       energy(:,-1) = energy(:,2)
    END IF

  END SUBROUTINE energy_bcs



  SUBROUTINE velocity_bcs

    CALL velocity_MPI

    IF (right == MPI_PROC_NULL .AND. xbc_right == BC_OTHER) THEN
       vx(nx+1,:) = vx(nx-1,:)
       vx(nx+2,:) = vx(nx-2,:)
       vy(nx+1,:) = vy(nx-1,:)
       vy(nx+2,:) = vy(nx-2,:)
       vz(nx+1,:) = vz(nx-1,:)
       vz(nx+2,:) = vz(nx-2,:)
    END IF
    IF (left == MPI_PROC_NULL .AND. xbc_left == BC_OTHER) THEN
       vx(-2,:) = vx(2,:)
       vx(-1,:) = vx(1,:)
       vy(-1,:) = vy(1,:)
       vy(-2,:) = vy(2,:)
       vz(-1,:) = vz(1,:)
       vz(-2,:) = vz(2,:)
    END IF

    IF (up == MPI_PROC_NULL .AND. ybc_up == BC_OTHER) THEN
       vx(:,ny+1) = vx(:,ny-1)
       vx(:,ny+2) = vx(:,ny-2)
       vy(:,ny+1) = vy(:,ny-1)
       vy(:,ny+2) = vy(:,ny-2)
       vz(:,ny+1) = vz(:,ny-1)
       vz(:,ny+2) = vz(:,ny-2)
    END IF
    IF (down == MPI_PROC_NULL .AND. ybc_down == BC_OTHER) THEN
       vx(:,-1) = vx(:,1)
       vx(:,-2) = vx(:,2)
       vy(:,-1) = vy(:,1)
       vy(:,-2) = vy(:,2)
       vz(:,-1) = vz(:,1)
       vz(:,-2) = vz(:,2)
    END IF

  END SUBROUTINE velocity_bcs



  SUBROUTINE remap_v_bcs

    REAL(num) :: theta

    CALL remap_v_MPI

    IF (right == MPI_PROC_NULL) THEN
       IF (xbc_right == BC_OTHER) THEN
          vx1(nx+1,:) = vx1(nx-1,:)
          vx1(nx+2,:) = vx1(nx-2,:)
          vy1(nx+1,:) = vy1(nx-1,:)
          vy1(nx+2,:) = vy1(nx-2,:)
          vz1(nx+1,:) = vz1(nx-1,:)
          vz1(nx+2,:) = vz1(nx-2,:)
       ELSE IF (xbc_right == BC_OPEN) THEN
          vx1(nx+1, 1:ny) = vx1(nx, 1:ny)
          vx1(nx+2, 1:ny) = vx1(nx, 1:ny)
          vy1(nx+1, 1:ny) = vy1(nx, 1:ny)
          vy1(nx+2, 1:ny) = vy1(nx, 1:ny)
          vz1(nx+1, 1:ny) = vz1(nx, 1:ny)
          vz1(nx+2, 1:ny) = vz1(nx, 1:ny)
       END IF
    END IF
    IF (left == MPI_PROC_NULL) THEN
       IF (xbc_left == BC_OTHER) THEN
          vx1(-2,:) = vx1(2,:)
          vx1(-1,:) = vx1(1,:)
          vy1(-1,:) = vy1(1,:)
          vy1(-2,:) = vy1(2,:)
          vz1(-1,:) = vz1(1,:)
          vz1(-2,:) = vz1(2,:)
       ELSE IF (xbc_left == BC_OPEN) THEN
          vx1(-2,:) = vx1(0,:) 
          vy1(-2,:) = vy1(0,:) 
          vz1(-2,:) = vz1(0,:) 
          vx1(-1,:) = vx1(0,:) 
          vy1(-1,:) = vy1(0,:) 
          vz1(-1,:) = vz1(0,:) 
       END IF
    END IF

    IF (up == MPI_PROC_NULL) THEN
       IF (ybc_up == BC_OTHER) THEN
          vx1(:,ny+1) = vx1(:,ny-1)
          vx1(:,ny+2) = vx1(:,ny-2)
          vy1(:,ny+1) = vy1(:,ny-1)
          vy1(:,ny+2) = vy1(:,ny-2)
          vz1(:,ny+1) = vz1(:,ny-1)
          vz1(:,ny+2) = vz1(:,ny-2)
       ELSE IF (ybc_up == BC_OPEN) THEN
          vx1(:,ny+1) = vx1(:,ny)
          vx1(:,ny+2) = vx1(:,ny)
          vy1(:,ny+1) = vy1(:,ny)
          vy1(:,ny+2) = vy1(:,ny)
          vz1(:,ny+1) = vz1(:,ny)
          vz1(:,ny+2) = vz1(:,ny)
       END IF
    END IF
    IF (down == MPI_PROC_NULL) THEN
       IF (ybc_down == BC_OTHER) THEN
          vx1(:,-1) = vx1(:,1)
          vx1(:,-2) = vx1(:,2)
          vy1(:,-1) = vy1(:,1)
          vy1(:,-2) = vy1(:,2)
          vz1(:,-1) = vz1(:,1)
          vz1(:,-2) = vz1(:,2)
       ELSE IF (ybc_down == BC_OPEN) THEN
          vx1(:,-2) = vx1(:,0)
          vy1(:,-2) = vy1(:,0)
          vz1(:,-2) = vz1(:,0)
          vx1(:,-1) = vx1(:,0)
          vy1(:,-1) = vy1(:,0)
          vz1(:,-1) = vz1(:,0)
       END IF
    END IF

  END SUBROUTINE remap_v_bcs



  SUBROUTINE density_bcs

    REAL(num) :: a,b

    CALL density_MPI

    IF (right == MPI_PROC_NULL .AND. xbc_right == BC_OTHER) THEN
       rho(nx+1,:) = rho(nx,:)
       rho(nx+2,:) = rho(nx-1,:)
    END IF
    IF (left == MPI_PROC_NULL .AND. xbc_left == BC_OTHER) THEN
       rho(0,:) = rho(1,:)
       rho(-1,:) = rho(2,:)
    END IF

    IF (up == MPI_PROC_NULL .AND. ybc_up == BC_OTHER) THEN
       rho(:,ny+1) = rho(:,ny)
       rho(:,ny+2) = rho(:,ny-1)
    END IF
    a=1.0_num/dyc(0)
    b=grav(0)/(gamma-1.0_num)
    IF (down == MPI_PROC_NULL .AND. ybc_down == BC_OTHER) THEN
       rho(:,0) = rho(:,1) * (a+2.0_num * b/(energy(:,0)+energy(:,1)))/&
            (a-2.0_num * b/(energy(:,0)+energy(:,1)))
       rho(:,-1) = rho(:,0)
    END IF

  END SUBROUTINE density_bcs




END MODULE boundary