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
  
!******************************************************************************
! This module contains the boundary conditions for the entire code
! Any new boundary conditions should be added here
!******************************************************************************

MODULE boundary

  USE shared_data
  USE mpiboundary
  USE random_generator

  IMPLICIT NONE

  REAL(num), DIMENSION(:), ALLOCATABLE :: drive_axis
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: drive_phase, drive_amp
  INTEGER :: drive_nel

  SAVE

CONTAINS

  !****************************************************************************
  ! Set up any necessary variables for the chosen boundary conditions
  !****************************************************************************

  SUBROUTINE set_boundary_conditions

    any_open = .FALSE.
    IF (xbc_min == BC_OPEN .OR. xbc_max == BC_OPEN &
        .OR. ybc_min == BC_OPEN .OR. ybc_max == BC_OPEN) any_open = .TRUE.

    IF (driven_boundary) CALL setup_driver_spectrum

  END SUBROUTINE set_boundary_conditions


  !****************************************************************************
  ! Routines to produce a spectral driver from a set of sine waves
  ! These are example driver routines setup for ybc_min only
  !****************************************************************************

  SUBROUTINE setup_driver_spectrum

    REAL(num) :: min_omega, max_omega
    INTEGER :: iel

    ! Set up a driver with 1000 elements
    drive_nel = 1000
    ALLOCATE(drive_axis(drive_nel))
    ALLOCATE(drive_amp(-2:nx+2,drive_nel))
    ALLOCATE(drive_phase(-2:nx+2,drive_nel))

    min_omega = 0.01_num
    max_omega = 10.0_num

    ! Initialize the random number generator. Change the seed to get
    ! different results
    CALL random_init_local(76783467)

    DO iel = 1, drive_nel
      ! Uniformly spaced frequency bins
      drive_axis(iel) = REAL(iel - 1, num) / REAL(drive_nel - 1, num) &
          * (max_omega - min_omega) + min_omega
      ! Random phase
      drive_phase(:,iel) = random() * 2.0_num * pi
      ! Kolmogorov amplitude
      drive_amp(:,iel) = 1.0e-4_num * drive_axis(iel)**(-2.5_num / 3.0_num)
    END DO

  END SUBROUTINE setup_driver_spectrum



  SUBROUTINE produce_spectrum(dat, time, rise_time)

    REAL(num), DIMENSION(-2:,:), INTENT(INOUT) :: dat
    REAL(num), INTENT(IN) :: time, rise_time
    REAL(num) :: val
    INTEGER :: iel, ix

    DO ix = -2, nx + 2
      val = 0.0_num
      DO iel = 1, drive_nel
        val = val + drive_amp(ix,iel) &
            * SIN(drive_axis(iel) * time + drive_phase(ix,iel))
      END DO
      dat(ix,:) = val
    END DO

    IF (time < rise_time) THEN
      dat = dat * 0.5_num * (1.0_num - COS(time * pi / rise_time))
    END IF

  END SUBROUTINE produce_spectrum


  !****************************************************************************
  ! Call all of the boundaries needed by the core Lagrangian solver
  !****************************************************************************

  SUBROUTINE boundary_conditions

    CALL bfield_bcs
    CALL energy_bcs
    CALL density_bcs
    CALL velocity_bcs

  END SUBROUTINE boundary_conditions



  !****************************************************************************
  ! Boundary conditions for magnetic field through plane
  !****************************************************************************

  SUBROUTINE bfield_bcs

    CALL bfield_mpi

    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_USER) THEN
      bx(-1,:) = bx(1,:)
      bx(-2,:) = bx(2,:)
      by( 0,:) = by(1,:)
      by(-1,:) = by(2,:)
      bz( 0,:) = bz(1,:)
      bz(-1,:) = bz(2,:)
    END IF

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_USER) THEN
      bx(nx+1,:) = bx(nx-1,:)
      bx(nx+2,:) = bx(nx-2,:)
      by(nx+1,:) = by(nx  ,:)
      by(nx+2,:) = by(nx-1,:)
      bz(nx+1,:) = bz(nx  ,:)
      bz(nx+2,:) = bz(nx-1,:)
    END IF

    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_USER) THEN
      bx(:, 0) = bx(:,1)
      bx(:,-1) = bx(:,2)
      by(:,-1) = by(:,1)
      by(:,-2) = by(:,2)
      bz(:, 0) = bz(:,1)
      bz(:,-1) = bz(:,2)
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_USER) THEN
      bx(:,ny+1) = bx(:,ny  )
      bx(:,ny+2) = bx(:,ny-1)
      by(:,ny+1) = by(:,ny-1)
      by(:,ny+2) = by(:,ny-2)
      bz(:,ny+1) = bz(:,ny  )
      bz(:,ny+2) = bz(:,ny-1)
    END IF

  END SUBROUTINE bfield_bcs



  SUBROUTINE bz_bcs

    CALL bz_mpi

    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_USER) THEN
      bz( 0,:) = bz(1,:)
      bz(-1,:) = bz(2,:)
    END IF

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_USER) THEN
      bz(nx+1,:) = bz(nx  ,:)
      bz(nx+2,:) = bz(nx-1,:)
    END IF

    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_USER) THEN
      bz(:, 0) = bz(:,1)
      bz(:,-1) = bz(:,2)
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_USER) THEN
      bz(:,ny+1) = bz(:,ny  )
      bz(:,ny+2) = bz(:,ny-1)
    END IF

  END SUBROUTINE bz_bcs



  !****************************************************************************
  ! Boundary conditions for specific internal energy
  !****************************************************************************

  SUBROUTINE energy_bcs

    CALL energy_mpi

    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_USER) THEN
      energy( 0,:) = energy(1,:)
      energy(-1,:) = energy(2,:)
    END IF

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_USER) THEN
      energy(nx+1,:) = energy(nx  ,:)
      energy(nx+2,:) = energy(nx-1,:)
    END IF

    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_USER) THEN
      energy(:, 0) = energy(:,1)
      energy(:,-1) = energy(:,2)
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_USER) THEN
      energy(:,ny+1) = energy(:,ny  )
      energy(:,ny+2) = energy(:,ny-1)
    END IF

  END SUBROUTINE energy_bcs



  SUBROUTINE temperature_bcs

    CALL temperature_mpi

    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_USER) THEN
      temperature( 0,:) = temperature(1,:)
      temperature(-1,:) = temperature(2,:)
    END IF

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_USER) THEN
      temperature(nx+1,:) = temperature(nx  ,:)
      temperature(nx+2,:) = temperature(nx-1,:)
    END IF

    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_USER) THEN
      temperature(:, 0) = temperature(:,1)
      temperature(:,-1) = temperature(:,2)
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_USER) THEN
      temperature(:,ny+1) = temperature(:,ny  )
      temperature(:,ny+2) = temperature(:,ny-1)
    END IF

  END SUBROUTINE temperature_bcs



  !****************************************************************************
  ! Boundary conditions for density
  !****************************************************************************

  SUBROUTINE density_bcs

    CALL density_mpi

    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_USER) THEN
      rho( 0,:) = rho(1,:)
      rho(-1,:) = rho(2,:)
    END IF

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_USER) THEN
      rho(nx+1,:) = rho(nx  ,:)
      rho(nx+2,:) = rho(nx-1,:)
    END IF

    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_USER) THEN
      rho(:, 0) = rho(:,1)
      rho(:,-1) = rho(:,2)
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_USER) THEN
      rho(:,ny+1) = rho(:,ny  )
      rho(:,ny+2) = rho(:,ny-1)
    END IF

  END SUBROUTINE density_bcs



  !****************************************************************************
  ! Full timestep velocity boundary conditions
  !****************************************************************************

  SUBROUTINE velocity_bcs

    CALL velocity_mpi

    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_USER) THEN
      vx(-2:0,:) = 0.0_num
      vy(-2:0,:) = 0.0_num
      vz(-2:0,:) = 0.0_num
    END IF

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_USER) THEN
      vx(nx:nx+2,:) = 0.0_num
      vy(nx:nx+2,:) = 0.0_num
      vz(nx:nx+2,:) = 0.0_num
    END IF

    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_USER) THEN
      vx(:,-2:0) = 0.0_num
      vy(:,-2:0) = 0.0_num
      CALL produce_spectrum(vz(:,-2:0), time, 1.0_num)
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_USER) THEN
      vx(:,ny:ny+2) = 0.0_num
      vy(:,ny:ny+2) = 0.0_num
      vz(:,ny:ny+2) = 0.0_num
    END IF

  END SUBROUTINE velocity_bcs



  !****************************************************************************
  ! Half timestep velocity boundary conditions
  !****************************************************************************

  SUBROUTINE remap_v_bcs

    CALL remap_v_mpi

    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_USER) THEN
      vx1(-2:0,:) = 0.0_num
      vy1(-2:0,:) = 0.0_num
      vz1(-2:0,:) = 0.0_num
    END IF

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_USER) THEN
      vx1(nx:nx+2,:) = 0.0_num
      vy1(nx:nx+2,:) = 0.0_num
      vz1(nx:nx+2,:) = 0.0_num
    END IF

    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_USER) THEN
      vx1(:,-2:0) = 0.0_num
      vy1(:,-2:0) = 0.0_num
      CALL produce_spectrum(vz1(:,-2:0), time - 0.5_num * dt, 1.0_num)
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_USER) THEN
      vx1(:,ny:ny+2) = 0.0_num
      vy1(:,ny:ny+2) = 0.0_num
      vz1(:,ny:ny+2) = 0.0_num
    END IF

  END SUBROUTINE remap_v_bcs



END MODULE boundary
