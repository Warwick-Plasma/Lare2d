MODULE setup

  USE shared_data
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: minimal_init,after_control
  PUBLIC ::  grid
  PUBLIC :: open_files, close_files, restart_data

  REAL(num), DIMENSION(:), ALLOCATABLE:: dxnew, dynew

CONTAINS

  SUBROUTINE minimal_init

    nprocx=0
    nprocy=0

    IF (num .EQ. 4) mpireal=MPI_REAL

  END SUBROUTINE



  SUBROUTINE after_control

    IF (.NOT. restart)  restart_snapshot = 0

    p_visc = 0.0_num
    eta = 0.0_num

    grav = 0.0_num
  END SUBROUTINE after_control


  
  SUBROUTINE grid                 ! stretched and staggered grid

    REAL(num) :: dx, dy, dxmin, dymin,xcstar,ycstar
    INTEGER :: ix, iy

    ALLOCATE (xb_global(-2:nx_global+2))
    ALLOCATE (yb_global(-2:ny_global+2))
    ALLOCATE (dxnew(-2:nx_global+2))
    ALLOCATE (dynew(-2:ny_global+2))

    dx = 1.0_num / REAL(nx_global, num)       ! initially assume uniform grid
    dy = 1.0_num / REAL(ny_global, num)

    length_x=x_end-x_start
    length_y=y_end-y_start

    ! grid cell boundary for x coordinates
    ! set to -0.5_num to have x=0 in centre of domain
    xb_global(0) = 0.0_num
    DO ix = -2, nx_global+2
       xb_global(ix) =  xb_global(0) + REAL(ix, num)*dx
    END DO
    xb_global = xb_global * length_x + x_start
    IF (x_stretch) CALL stretch_x     ! stretch grid ?

    !define position of ghost cells using sizes of adjacent cells
    IF (xbc_right == periodic) THEN
       xb_global(nx_global+1) = xb_global(nx_global) &
            + (xb_global(1) - xb_global(0))
       xb_global(nx_global+2) = xb_global(nx_global) &
            + (xb_global(2) - xb_global(0))
       xb_global(-1) = xb_global(0) &
            - (xb_global(nx) - xb_global(nx-1))
       xb_global(-2) = xb_global(0) &
            - (xb_global(nx) - xb_global(nx-2))
    ELSE
       xb_global(nx_global+1) = 2.0_num * xb_global(nx_global) - xb_global(nx_global-1)
       xb_global(nx_global+2) = 2.0_num * xb_global(nx_global) - xb_global(nx_global-2)
       xb_global(-1) = 2.0_num * xb_global(0) - xb_global(1)
       xb_global(-2) = 2.0_num * xb_global(0) - xb_global(2)
    END IF

    xb = xb_global(coordinates(2)*nx-2:coordinates(2)*nx+nx+2)

    DO ix = -1, nx+2
       ixm = ix - 1
       xc(ix) = 0.5_num*(xb(ixm) + xb(ix))     ! cell centre
    END DO
    DO ix = -1, nx+1
       ixp = ix + 1
       dxc(ix) = xc(ixp) - xc(ix)    ! distance between centres
    END DO
    IF (coordinates(2)==nprocx-1) THEN
       dxc(nx+2) = dxc(nx+1)
    ELSE
       xcstar = 0.5_num*(xb(nx+2) + xb_global(coordinates(2)*nx+nx+3))
       dxc(nx+2) = xcstar - xc(nx+2)
    END IF
    DO ix = -1, nx+2
       ixm = ix - 1
       dxb(ix) = xb(ix) - xb(ixm)    ! cell width
    END DO

    yb_global(0) = 0.0_num                   ! repeat for y
    DO iy = -2, ny_global+2
       yb_global(iy) =  yb_global(0) + REAL(iy, num)*dy
    END DO
    yb_global = yb_global * length_y + y_start
    IF (y_stretch) CALL stretch_y
    IF (ybc_up == periodic) THEN
       yb_global(ny_global+1) = yb_global(ny_global) &
            + (yb_global(1) - yb_global(0))
       yb_global(ny_global+2) = yb_global(ny_global) &
            + (yb_global(2) - yb_global(0))
       yb_global(-1) = yb_global(0) &
            - (yb_global(ny) - yb_global(ny-1))
       yb_global(-2) = yb_global(0) &
            - (yb_global(ny) - yb_global(ny-2))
    ELSE
       yb_global(ny_global+1) = 2.0_num * yb_global(ny_global) - yb_global(ny_global-1)
       yb_global(ny_global+2) = 2.0_num * yb_global(ny_global) - yb_global(ny_global-2)
       yb_global(-1) = 2.0_num * yb_global(0) - yb_global(1)
       yb_global(-2) = 2.0_num * yb_global(0) - yb_global(2)
    END IF

    yb = yb_global(coordinates(1)*ny-2:coordinates(1)*ny+ny+2)

    DO iy = -1, ny+2
       iym = iy - 1
       yc(iy) = 0.5_num*(yb(iym) + yb(iy))
    END DO
    DO iy = -1, ny+1
       iyp = iy + 1
       dyc(iy) = yc(iyp) - yc(iy)
    END DO
    IF (coordinates(1)==nprocy-1) THEN
       dyc(ny+2) = dyc(ny+1)
    ELSE
       ycstar = 0.5_num*(yb(ny+2) + yb_global(coordinates(1)*ny+ny+3))
       dyc(ny+2) = ycstar - yc(ny+2)
    END IF
    DO iy = -1, ny+2
       iym = iy - 1
       dyb(iy) = yb(iy) - yb(iym)
    END DO
    DO ix = -1, nx+2
       DO iy = -1, ny+2
          cv(ix,iy) = dxb(ix) * dyb(iy)     ! define the cell area
       END DO
    END DO

    CALL MPI_ALLREDUCE(MINVAL(dxb(1:nx)), dxmin, 1, mpireal, MPI_MIN, &
         comm, errcode)
    CALL MPI_ALLREDUCE(MINVAL(dyb(1:ny)), dymin, 1, mpireal, MPI_MIN, &
         comm, errcode)
    min_grid_spacing = MIN(dxmin, dymin)

    DEALLOCATE( dxnew, dynew)

  END SUBROUTINE grid



  SUBROUTINE stretch_x   ! replace with any stretching algorithm as needed

    REAL(num) :: width, dx, L, f, lx_new

    lx_new = 200.0_num                ! new tolal length
    L = length_x / 1.5_num       ! centre of tanh stretching in unstretched coordinates
    width = length_x / 10.0_num       ! width of tanh stretching in unstretched coordinates

    f = (lx_new - length_x)/(length_x - L)/2.0_num

    dx = length_x / REAL(nx_global,num)  
    dxnew = dx + f*(1.0_num+TANH((ABS(xb_global)-L)/width))*dx
    
!!$    DO ix = nx_global/2+1, nx_global+2
!!$       xb_global(ix) = xb_global(ix-1) + dxnew(ix)
!!$    ENDDO
!!$    DO ix = nx_global/2-1, -2, -1
!!$       xb_global(ix) = xb_global(ix+1) - dxnew(ix)
!!$    ENDDO

    DO ix = 1, nx_global+2
       xb_global(ix) = xb_global(ix-1) + dxnew(ix)
    ENDDO

    length_x = lx_new
    
  END SUBROUTINE stretch_x



  SUBROUTINE stretch_y !stretch domain upwards only

    REAL(num) :: width, dy, L, f, ly_new

    ly_new = 100.0_num                ! new tolal length
    L = length_y / 1.5_num       ! centre of tanh stretching in unstretched coordinates
    width = length_y / 10.0_num       ! width of tanh stretching in unstretched coordinates

    f = (ly_new - length_y)/(length_y - L)/2.0_num

    dy = length_y / REAL(ny_global,num)  
    dynew = dy + f*(1.0_num+TANH((ABS(yb_global)-L)/width))*dy

!!$    DO iy = ny_global/2+1, ny_global+2
!!$       yb_global(iy) = yb_global(iy-1) + dynew(iy)
!!$    ENDDO
!!$    DO iy = ny_global/2-1, -2, -1
!!$       yb_global(iy) = yb_global(iy+1) - dynew(iy)
!!$    ENDDO

    DO iy = 1, ny_global+2
       yb_global(iy) = yb_global(iy-1) + dynew(iy)
    ENDDO
    
    length_y = ly_new

  END SUBROUTINE stretch_y



  SUBROUTINE open_files

    CHARACTER(LEN=11+Data_Dir_Max_Length) :: file2
    CHARACTER(LEN=7+Data_Dir_Max_Length) :: file3

    IF (rank == 0) THEN
       WRITE(file2, '(a,"/lare2d.dat")') TRIM(data_dir)
       OPEN(unit=20, STATUS = 'REPLACE',FILE = file2)
       WRITE(file3, '(a,"/en.dat")') TRIM(data_dir)
       OPEN(unit=30, STATUS = 'REPLACE',FILE = file3,FORM="binary")
    END IF

  END SUBROUTINE open_files



  SUBROUTINE close_files

    IF (rank == 0) THEN
       CLOSE(unit=20)
       CLOSE(unit=30)
    END IF

  END SUBROUTINE close_files



  SUBROUTINE restart_data
    CHARACTER(LEN=20+Data_Dir_Max_Length) :: filename
    INTEGER :: filehandle, localcellcount, tn(2),sof
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: data

#ifdef NONMPIIO
    print*, "Cannot restart from non-MPI io yet"
    STOP
#else

    ! Create the filename for the last snapshot
#ifdef MHDCLUSTER
    WRITE(filename, '("nfs:",a,"/",i4.4,".lld")') TRIM(data_dir), restart_snapshot
#else
    WRITE(filename, '(a,"/",i4.4,".lld")') TRIM(data_dir), restart_snapshot
#endif
    ! Open the file
    CALL MPI_FILE_OPEN(comm, filename, MPI_MODE_RDONLY, &
         MPI_INFO_NULL, filehandle, errcode)

    CALL MPI_FILE_READ_ALL(filehandle,tn,2,MPI_INTEGER,status,errcode)
    CALL MPI_FILE_READ_ALL(filehandle,sof,1,MPI_INTEGER,status,errcode)
    CALL MPI_FILE_READ_ALL(filehandle,time,1,mpireal,status,errcode)


    IF ((tn(1) /= nx_global).OR.(tn(2)/=ny_global) .OR. (sof/=num)) THEN
       IF (rank ==0) THEN
          WRITE(20,*) 'Error: global dimensions do not match restart file.'
       ENDIF
       CALL MPI_FILE_CLOSE(filehandle, errcode)
       CALL close_files
       CALL MPI_FINALIZE(errcode)
       STOP
    ENDIF

    ! Set my view of the file
    CALL MPI_FILE_SET_VIEW(filehandle, initialdisp, mpireal, subtype,&
         "native", MPI_INFO_NULL, errcode)

    localcellcount = (nx+1) * (ny+1)
    ALLOCATE(data(0:nx,0:ny))
    CALL MPI_FILE_READ_ALL(filehandle,data, localcellcount, mpireal, status, errcode)
    rho(0:nx,0:ny) = data
    CALL MPI_FILE_READ_ALL(filehandle,data, localcellcount, mpireal, status, errcode)
    energy(0:nx,0:ny) = data / (gamma - 1.0_num)
    CALL MPI_FILE_READ_ALL(filehandle,data, localcellcount, mpireal, status, errcode)
    vx(0:nx,0:ny) = data
    CALL MPI_FILE_READ_ALL(filehandle,data, localcellcount, mpireal, status, errcode)
    vy(0:nx,0:ny) = data
    CALL MPI_FILE_READ_ALL(filehandle,data, localcellcount, mpireal, status, errcode)
    vz(0:nx,0:ny) = data
    CALL MPI_FILE_READ_ALL(filehandle,data, localcellcount, mpireal, status, errcode)
    bx(0:nx,0:ny) = data
    CALL MPI_FILE_READ_ALL(filehandle,data, localcellcount, mpireal, status, errcode)
    by(0:nx,0:ny) = data
    CALL MPI_FILE_READ_ALL(filehandle,data, localcellcount, mpireal, status, errcode)
    bz(0:nx,0:ny) = data

    DEALLOCATE(data)

#endif

  END SUBROUTINE restart_data



END MODULE setup
