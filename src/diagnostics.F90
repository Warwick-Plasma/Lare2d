!*************************************************************************
! Controls al I/O and diagnostics. Output files are 'lare2d.dat',
! 'control.dat', 'en.dat' and a series of snapshots in 'fort.5x'
! The idl package in 'plot.pro' gives simple loading and surface
! plotting based on these files. This isn't documented but is very simple!
!*************************************************************************
MODULE diagnostics

  USE shared_data ; USE boundary
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: set_dt, output_routines, energy_correction

CONTAINS


  SUBROUTINE output_routines(i)   ! i=step index
    ! if halt set to false then code stops
    INTEGER, INTENT(IN) :: i

    INTEGER, SAVE :: output_file = 0
    INTEGER, PARAMETER :: out = 1000
    INTEGER, SAVE :: index = 1, step = 1
    INTEGER :: filehandle, localcellcount
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: data
    LOGICAL :: print_arrays, last_call

    !this output routine uses the same sturcture as needed for mpi output
    !this is more complicated than need for the serial code
    !rank always equals zero in this serial code
    CHARACTER(LEN=12+Data_Dir_Max_Length) :: filename

    REAL(num) :: t_out = 0.0_num
    REAL(num) :: en_ke = 0.0_num, en_int = 0.0_num
    REAL(num) :: en_b = 0.0_num, heating_visc = 0.0_num
    REAL(num) :: heating_ohmic = 0.0_num
    REAL(num) :: j_max_local, total

    IF (nsteps >= out) step = nsteps/out + 1 ! make sure output fits arrays
    IF (i == 0 .AND. rank == 0) THEN         ! done just once at the start
       CALL output_log
       IF (.NOT. restart) WRITE(30) num, 6
    END IF
    IF (i == 0) output_file = output_file + restart_snapshot

    IF (MOD(i,step) .EQ. 0 .OR. last_call) THEN       ! do every (step) steps
       t_out = time

       CALL energy_account(en_b, en_ke, en_int)

       CALL MPI_ALLREDUCE(total_visc_heating, total, 1, mpireal, MPI_SUM, &
            comm, errcode)
       heating_visc = total

       CALL MPI_ALLREDUCE(total_ohmic_heating, total, 1, mpireal, MPI_SUM, &
            comm, errcode)
       heating_ohmic = total
       IF (rank .EQ. 0) THEN
          WRITE(30) t_out, en_b, en_ke, en_int
          WRITE(30) heating_visc, heating_ohmic
       ENDIF

       index = index + 1
    END IF

    CALL io_test(i, print_arrays, last_call) !check if snapshot is needed

    IF (print_arrays) THEN   ! output a snapshot of arrays
 
#ifdef NONMPIIO
       ! The old output routines
       WRITE(filename, '(a,"/",i3.3,i4.4,".dat")') TRIM(data_dir), rank, output_file
       !after all of this charcter manipulation the filename for output
       !is of the form rrrnnnn.dat where rrr is the rank of the process
       !and nnnn is the snapshot number

       OPEN(unit=50, FORM = 'UNFORMATTED',STATUS = 'REPLACE', &
            FILE = TRIM(filename))
       WRITE(50) REAL(nx_global,num)-0.1_num,REAL(ny_global,num)-0.1_num
       WRITE(50)  REAL(nx,num)-0.1_num,REAL(ny,num)-0.1_num
       WRITE(50) time
       WRITE(50) rho(1:nx,1:ny),(gamma - 1.0_num) * energy(1:nx,1:ny)
       WRITE(50) vx(1:nx,1:ny), vy(1:nx,1:ny), &
            vz(1:nx,1:ny)
       WRITE(50) bx(1:nx,1:ny), by(1:nx,1:ny), &
            bz(1:nx,1:ny)
       WRITE(50) xb(1:nx), yb(1:ny)
       WRITE(50) xb(1:nx+1), yb(1:ny+1)
       CLOSE(50)
#else       
       ! MPI file output routines  --------------------------------------

       ! Create the filename
#ifdef MHDCLUSTER
       WRITE(filename, '("nfs:",a,"/",i4.4,".lld")') TRIM(data_dir), output_file
#else
       WRITE(filename, '(a,"/",i4.4,".lld")') TRIM(data_dir)	, output_file
#endif
       ! Rank 0 deletes the file as there is no replace option on the file open,
       ! not deleting causes literal overwriting and potential confusion.
       IF (rank == 0) CALL MPI_FILE_DELETE(TRIM(filename), MPI_INFO_NULL, errcode)

       ! Wait for the deletetion before we open the file...
       CALL MPI_BARRIER(comm,errcode)

       ! Open the file, filehandle is the MPI file unit number.
       CALL MPI_FILE_OPEN(comm, TRIM(filename), MPI_MODE_CREATE + MPI_MODE_WRONLY, &
            MPI_INFO_NULL, filehandle, errcode)

       ! If rank 0 then write the file header.
       IF (rank == 0) THEN
          CALL MPI_FILE_WRITE(filehandle,(/nx_global, ny_global/), 2, MPI_INTEGER, &
               status, errcode)
          CALL MPI_FILE_WRITE(filehandle,num,1,MPI_INTEGER,status,errcode)
          CALL MPI_FILE_WRITE(filehandle,time,1,mpireal,status,errcode)
          CALL MPI_FILE_WRITE(filehandle,xb_global(0:nx_global),nx_global+1,mpireal,status,errcode)
          CALL MPI_FILE_WRITE(filehandle,yb_global(0:ny_global),ny_global+1,mpireal,status,errcode)
       ENDIF

       ! Set my view of the file
       CALL MPI_FILE_SET_VIEW(filehandle, initialdisp, mpireal, subtype,&
            "native", MPI_INFO_NULL, errcode)

       localcellcount = (nx+1) * (ny+1)
       ALLOCATE(data(0:nx,0:ny))

       data = rho(0:nx,0:ny) 
       CALL MPI_FILE_WRITE_ALL(filehandle,data, localcellcount, mpireal, status, errcode)
       data = (gamma-1.0_num) * energy(0:nx,0:ny) 
       CALL MPI_FILE_WRITE_ALL(filehandle,data, localcellcount, mpireal, status, errcode)
       data = vx(0:nx,0:ny) 
       CALL MPI_FILE_WRITE_ALL(filehandle,data, localcellcount, mpireal, status, errcode)
       data = vy(0:nx,0:ny)
       CALL MPI_FILE_WRITE_ALL(filehandle,data, localcellcount, mpireal, status, errcode)
       data = vz(0:nx,0:ny) 
       CALL MPI_FILE_WRITE_ALL(filehandle,data, localcellcount, mpireal, status, errcode)
       data = bx(0:nx,0:ny) 
       CALL MPI_FILE_WRITE_ALL(filehandle,data, localcellcount, mpireal, status, errcode)
       data = by(0:nx,0:ny)
       CALL MPI_FILE_WRITE_ALL(filehandle,data, localcellcount, mpireal, status, errcode)
       data = bz(0:nx,0:ny)
       CALL MPI_FILE_WRITE_ALL(filehandle,data, localcellcount, mpireal, status, errcode)

       CALL MPI_FILE_CLOSE(filehandle, errcode)
       DEALLOCATE(data)
       IF (rank == 0) THEN
          WRITE(20,*) 'Succesfully dumped data at time ',time,' seconds after ',i,' iterations'
          !FLUSH(20)
       ENDIF
#endif

       output_file = output_file + 1

    END IF

    IF (last_call .AND. rank == 0) THEN   ! output energy diagnostics etc
       WRITE(20,*) 'final nsteps / time =', i, time
    END IF

  END SUBROUTINE output_routines



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



  SUBROUTINE set_dt        ! sets CFL limited step
! Assumes all variables are defined at the same point. Be careful
! with setting 'dt_multiplier' if you expect massive changes across
! cells.

    REAL(num) :: cons, dt1, dt2, dt3, dt4, dt5, dt_local, dxlocal
    REAL(num) :: dtr_local, dth_local

    dt_local = largest_number
    dtr_local = largest_number
    dth_local = largest_number
    cons = gamma * (gamma - 1.0_num)

    DO iy = 0, ny
       DO ix = 0, nx
          ixm = ix - 1
          iym = iy - 1

          w1 = (bx(ix,iy)**2 + by(ix,iy)**2 + bz(ix,iy)**2) 
          w2 = SQRT((w1 + energy(ix,iy)*rho(ix,iy)*cons) / MAX(rho(ix,iy), none_zero))  &
               + 2.0_num * SQRT(p_visc(ix,iy) / MAX(rho(ix,iy), none_zero))
          w2 = w2 * (1.0_num + visc1)
          dt1 = dxb(ix) / (w2 + ABS(vx(ix,iy)))
          dt2 = dyb(iy) / (w2 + ABS(vy(ix,iy)))
! find ideal MHD CFL limit
          dt_local = MIN(dt_local, dt1, dt2)
! note resistive limits assumes uniform resistivity hence cautious factor 0.2
          dxlocal = 1.0_num / (1.0_num/dxb(ix)**2 + 1.0_num/dyb(iy)**2)
          dt3 = 0.2_num * dxlocal / MAX(eta(ix,iy), none_zero)
! Hall MHD CFL limit
          dt4 = 0.75_num * rho(ix,iy) * MIN(dxb(ix),dyb(iy))**2 &
               / MAX(lambda_i(ix,iy), none_zero) /  MAX(SQRT(w1), none_zero)
! adjust to accomodate resistive effects
          dtr_local = MIN(dtr_local, dt3)
          dth_local = MIN(dth_local, dt4)
! correct to stop overlapping of Lagrangian cells	
          w1 = ABS(vx(ix,iy) - vx(ixm,iy)) / dxb(ix)   &
               + ABS(vy(ix,iy) - vy(ix,iym)) / dyb(iy)
          dt5 = 1.0_num / MAX(w1, none_zero)
          dt_local = MIN(dt_local, dt5)
       END DO
    END DO

    CALL MPI_ALLREDUCE(dt_local, dt, 1, mpireal, MPI_MIN, comm, errcode)
    CALL MPI_ALLREDUCE(dtr_local, dtr, 1, mpireal, MPI_MIN, comm, errcode)
    CALL MPI_ALLREDUCE(dth_local, dth, 1, mpireal, MPI_MIN, comm, errcode)

    dtr = dt_multiplier * dtr
    dth = dt_multiplier * dth
    dt = dt_multiplier * dt

    time = time + dt

  END SUBROUTINE set_dt



  SUBROUTINE energy_account(energy_b, energy_ke, energy_int)

    REAL(num), INTENT(OUT) :: energy_b, energy_ke, energy_int
    REAL(num) :: energy_b_local, energy_ke_local, energy_int_local
    REAL(num) :: cv_v, rho_v

    energy_b_local = 0.0_num
    energy_ke_local= 0.0_num
    energy_int_local = 0.0_num

    DO iy = 1, ny
       DO ix = 1, nx
          ixm = ix - 1
          iym = iy - 1
          ixp = ix + 1
          iyp = iy + 1
          w2 = (bx(ix,iy)**2 + bx(ixm,iy)**2)/2.0_num
          w3 = (by(ix,iy)**2 + by(ix,iym)**2)/2.0_num
          w4 = bz(ix,iy)**2
          w1 = (w2 + w3 + w4) / 2.0_num
          energy_b_local = energy_b_local + w1*cv(ix,iy)

          energy_int_local = energy_int_local &
               + energy(ix,iy)*rho(ix,iy)*cv(ix,iy)

!WARNING the KE is summed on the vertices
          rho_v = (rho(ix,iy)*cv(ix,iy) + rho(ixp,iy)*cv(ixp,iy) &
               + rho(ix,iyp)*cv(ix,iyp) +rho(ixp,iyp)*cv(ixp,iyp))
          rho_v = rho_v / (cv(ix,iy) + cv(ixp,iy) + cv(ix,iyp) + cv(ixp,iyp))
          cv_v = (cv(ix,iy) + cv(ixp,iy) + cv(ix,iyp)+ cv(ixp,iyp))/4.0_num
          w1 = rho_v * (vx(ix,iy)**2 + vy(ix,iy)**2 + vz(ix,iy)**2) * cv_v
          energy_ke_local = energy_ke_local + w1 / 2.0_num
       END DO
    END DO

    energy_b_local = energy_b_local

    CALL MPI_ALLREDUCE(energy_ke_local, energy_ke, 1, mpireal, MPI_SUM, &
         comm, errcode)
    CALL MPI_ALLREDUCE(energy_b_local, energy_b, 1, mpireal, MPI_SUM, &
         comm, errcode)
    CALL MPI_ALLREDUCE(energy_int_local, energy_int, 1, mpireal, MPI_SUM, &
         comm, errcode)

  END SUBROUTINE energy_account



  SUBROUTINE energy_correction

    delta_ke = - delta_ke
    WHERE (delta_ke < 0.0_num) delta_ke = 0.0_num
    delta_ke = delta_ke / (rho * cv)

    DO iy = 1, ny
       DO ix = 1, nx
          energy(ix,iy) = energy(ix,iy) + delta_ke(ix,iy) 
       END DO
    END DO
    CALL energy_bcs

  END SUBROUTINE energy_correction



  SUBROUTINE output_log  !writes basic data to 'lare2d.dat'

    REAL(num) :: temp
    
    WRITE(20,*) 'nprocx, nprocy = ', nprocx, nprocy
    WRITE(20,*) 'nx, ny = ', nx, ny
    WRITE(20,*)
    WRITE(20,*) 'length_x = ', length_x 
    WRITE(20,*) 'length_y = ', length_y 
    WRITE(20,*)
#ifndef Q_MONO
    WRITE(20,*) 'tensor shock viscosity'
#else
    WRITE(20,*) 'q_mono viscosity'
#endif
    WRITE(20,*) 'linear viscosity coeff = ', visc1
    WRITE(20,*) 'quadratic viscosity coeff = ', visc2
    WRITE(20,*) 'uniform tensor viscosity coeff = ', visc3 
    WRITE(20,*) 'j_max = ', j_max 
    WRITE(20,*) 'vc = ', vc 
    WRITE(20,*) 'eta0 = ', eta0 
    WRITE(20,*) 'eta_background = ', eta_background 
    WRITE(20,*)
    WRITE(20,*) 't_start, t_end = ', time, t_end
    WRITE(20,*) 'nsteps =',nsteps
    WRITE(20,*)

  END SUBROUTINE output_log

END MODULE diagnostics
