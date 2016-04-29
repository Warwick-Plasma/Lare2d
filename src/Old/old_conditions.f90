! This file contains example initial conditions used in previous simulations



    SUBROUTINE set_initial_conditions
    ! Simple Harris current sheet
    ! Normalise equations
    INTEGER:: ix, iy
    REAL(num) :: x0,lamda,rhoinf,rho0,b0

    bx = 0.0_num
    by = 0.0_num
    bz = 0.0_num
    vx = 0.0_num
    vy = 0.0_num
    vz = 0.0_num

    lamda = 1.0_num
    rho0 = 1.0_num
    rhoinf = rho0 * 0.2_num
    b0 = 1.0_num

    DO iy = -1, ny+2
       bx(:,iy) = b0* tanh(yc(iy) / lamda)
       rho(:,iy) = rho0 * (1.0_num / COSH(yc(iy) / lamda))**2 + rhoinf
    END DO

    energy = 0.75_num

  END SUBROUTINE set_initial_conditions



SUBROUTINE set_initial_conditions
    INTEGER:: ix, iy

    ! strong shock test
    ! standard test with nx=1200 and Lx=1
    ! run to t=0.038
    !Normalised equations

    xc = xc + length_x / 2.0_num
    xb = xb + length_x / 2.0_num
    yc = yc + length_y / 2.0_num
    yb = yb + length_y / 2.0_num

    vx = 0.0_num
    vy = 0.0_num
    vz = 0.0_num
    bx = 0.0_num
    by = 0.0_num
    bz = 0.0_num
    rho = 1.0_num

    gamma = 1.4_num

    DO ix = -1, nx+2
       IF (xc(ix) < 0.1_num) THEN
          energy(ix,:) = 1000.0_num / (gamma - 1.0_num)
       ELSE IF (xc(ix) < 0.9_num) THEN
          energy(ix,:) = 0.01_num / (gamma - 1.0_num)
       ELSE
          energy(ix,:) = 100.0_num / (gamma - 1.0_num)
       END IF
    END DO

  END SUBROUTINE set_initial_conditions



  SUBROUTINE set_initial_conditions

    INTEGER:: ix, iy

    ! 1d MHD Brio-Wu
    ! typically run with 800 points until t=0.1
    !normalised equations

    gamma = 1.4_num

    vx = 0.0_num
    vy = 0.0_num
    vz = 0.0_num
    bx = 0.75_num
    by = 0.0_num
    bz = 0.0_num

    DO ix = -1, nx+2
       IF (xc(ix) >= 0.0_num) THEN
          rho(ix,:) = 0.125_num
          energy(ix,:) = 0.1_num / (gamma - 1.0_num)
          by(ix,:) = -1.0_num
       ELSE
          rho(ix,:) = 1.0_num
          energy(ix,:) = 1.0_num / (gamma - 1.0_num)
          by(ix,:) = 1.0_num
       END IF
    END DO

    energy = energy / rho

  END SUBROUTINE set_initial_conditions



  SUBROUTINE set_initial_conditions

    INTEGER:: ix, iy
    REAL(num) :: b1, b2, x, y, s

    ! 2d MHD Brio-Wu
    ! needs open bcs. OK along diagonal with periodic bcs
    ! with Lx=Ly=1/SQRT(2) and t_end=0.1
    ! normalised equations
    vx = 0.0_num
    vy = 0.0_num
    vz = 0.0_num
    bz = 0.0_num
    b1 = 0.75_num
    b2 = 1.0_num

    DO ix = -1, nx+2
       DO iy = -1, ny+2
          x = xc(ix) + 0.5_num / SQRT(2.0_num)
          y = yc(iy) + 0.5_num / SQRT(2.0_num)
          s = (x + y) / SQRT(2.0_num)
          IF (s > 0.5_num + 0.1_num*dxb(1)) THEN
             rho(ix,iy) = 0.125_num
             energy(ix,iy) = 0.1_num / (gamma - 1.0_num)
             bx(ix,iy) = (b1 + b2) / SQRT(2.0_num)
             by(ix,iy) = (b1 - b2) / SQRT(2.0_num)
          ELSE
             rho(ix,iy) = 1.0_num
             energy(ix,iy) = 1.0_num / (gamma - 1.0_num)
             bx(ix,iy) = (b1 - b2) / SQRT(2.0_num)
             by(ix,iy) = (b1 + b2) / SQRT(2.0_num)
          END IF
       END DO
    END DO
    bx(-2,:) = bx(-1,:)
    by(:,-2) = by(:,-1)
    energy = energy / rho

  END SUBROUTINE set_initial_conditions



  SUBROUTINE set_initial_conditions

    INTEGER:: ix, iy
    ! Ortzang-Tang
    ! needs Lx=Ly=2*pi periodic bcs
    ! run to t=pi
    !normalise equations

    xc = xc + length_x / 2.0_num
    xb = xb + length_x / 2.0_num
    yc = yc + length_y / 2.0_num
    yb = yb + length_y / 2.0_num

    gamma = 1.667_num

    pressure = gamma
    rho = 2.778_num
    energy = pressure / ((gamma - 1.0_num)*rho)
    vz = 0.0_num
    bz = 0.0_num

    DO ix = -1, nx+2
       DO iy = -1, ny+2
          bx(ix,iy) = - SIN(yc(iy))
          by(ix,iy) = SIN(2.0_num*xc(ix))
          vx(ix,iy) = - SIN(yb(iy))
          vy(ix,iy) = SIN(xb(ix))
       END DO
    END DO
    bx(-2,:) = bx(nx-2,:)
    by(:,-2) = by(:,ny-2)
    vx(-2,:) = vx(nx-2,:)
    vy(:,-2) = vy(:,ny-2)

  END SUBROUTINE set_initial_conditions



  SUBROUTINE set_initial_conditions

    INTEGER:: ix, iy
    REAL(num) :: angle
    REAL(num) :: k, a, cs, r, vpar, vperp, bpar, bperp

    ! 2d spherical Alfven wave
    ! set angle=pi/6 to get Alfven wave test in JCP
    ! assumes periodic bcs. Lx=1/cos(angle). Ly=1/sin(angle)
    !normalised equations

    angle = pi / 6.0_num

    rho = 1.0_num
    energy = 1.0_num
    vx = 0.0_num
    vy = 0.0_num
    vz = 0.0_num
    bx = 0.0_num
    by = 0.0_num
    bz = 0.0_num
    a = 0.1_num
    cs = SQRT(gamma*(gamma-1.0_num))
    k = 2.0_num * pi
    xb = xb + 0.5_num * length_x
    xc = xc + 0.5_num * length_x
    yb = yb + 0.5_num * length_y
    yc = yc + 0.5_num * length_y
    DO ix = -1, nx+2
       DO iy = -1, ny+2
          r = xb(ix)*COS(angle) + yb(iy)*SIN(angle)
          vpar = 0.0_num
          vperp = a * SIN(k*r)
          vx(ix,iy) = vpar * COS(angle) - vperp * SIN(angle)
          vy(ix,iy) = vpar * SIN(angle) + vperp * COS(angle)
          vz(ix,iy) = a * COS(k*r)

          r = xb(ix)*COS(angle) + yc(iy)*SIN(angle)
          bpar = 1.0_num
          bperp = a * SIN(k*r)
          bx(ix,iy) = bpar * COS(angle) - bperp * SIN(angle)

          r = xc(ix)*COS(angle) + yb(iy)*SIN(angle)
          bpar = 1.0_num
          bperp = a * SIN(k*r)
          by(ix,iy) = bpar * SIN(angle) + bperp * COS(angle)

          r = xc(ix)*COS(angle) + yc(iy)*SIN(angle)
          bz(ix,iy) = a * COS(k*r)
       END DO
    END DO
    bx(-2,:) = bx(nx-2,:)
    by(:,-2) = by(:,ny-2)
    vx(-2,:) = vx(nx-2,:)
    vy(:,-2) = vy(:,ny-2)

    energy = 1.0_num / (gamma - 1.0_num)

  END SUBROUTINE set_initial_conditions



  SUBROUTINE set_initial_conditions
! Startified model atmosphere that uses the neutral fraction to correctly
! include the effects of pratial pressures
! S.I. units
    INTEGER:: ix, iy, icycle,flipy

    REAL(num)::dg,m=1.5_num,a=1.0_num,y_cor=3.75e6,t_ph=6420.0_num
    REAL(num) :: t_cor=963000.0_num,wtr=7.5e5_num,T
    REAL(num) :: yt,q,r,bphi,b1,pb,b,b_0
    INTEGER :: eos_this_cycle,max_cycles

    vx = 0.0_num
    vy = 0.0_num
    vz = 0.0_num
    bx = 0.0_num
    by = 0.0_num
    bz = 0.0_num
    energy = 1.0_num

    gamma = 5.0_num / 3.0_num

    grav = 274.0_num
    IF (coordinates(1) == 0) grav(0) = 0.0_num
    IF (coordinates(1) == nprocy-1) grav(ny) = 0.0_num

    max_cycles=100
    IF (eos_number == EOS_IDEAL) max_cycles=1

    !First cycle always use ideal gas equation of state for first guess
    eos_this_cycle=EOS_IDEAL
    DO icycle=1,max_cycles
     !define energy profile
      DO ix = -1,nx+2
          DO iy = -1,ny+2
            IF (yc(iy) < 0.0_num) THEN
                T = t_ph - (a * grav(iy) * yc(iy) * mbar/kb/ (m+1.0_num) )
            ELSE
              T = t_ph + &
                   ((t_cor-t_ph)*0.5_num * (TANH((yc(iy)-y_cor)/wtr)+1.0_num))
            END IF
            CALL Get_Energy(rho(ix,iy),T,eos_this_cycle,ix,iy,energy(ix,iy))
          END DO
        END DO

        CALL neutral_fraction(eos_number)
        !solve HS eqn to get rho profile
        rho = 2.7e-4_num

     !wait for lower process
     IF (coordinates(1)/=0) THEN
        CALL MPI_RECV(rho(:,-1),nx+4,mpireal,proc_y_min,tag,comm,status,errcode)
     END IF
     !calculate the density profile in y for the specified energy profile
     DO ix = -1,nx+2
        DO iy = 1,ny+2
           dg = 1.0_num/(dyb(iy)+dyb(iy-1))
           rho(ix,iy) = rho(ix,iy-1)*(energy(ix,iy-1)/dyc(iy-1)*(gamma-1.0_num)&
                -grav(iy-1)*dyb(iy-1)*dg)
           rho(ix,iy) = rho(ix,iy)/(energy(ix,iy)/dyc(iy-1)*(gamma-1.0_num)&
                +grav(iy-1)*dyb(iy)*dg)
        END DO
     END DO
     !now send to higher process
     IF (coordinates(1)/=nprocy-1) THEN
        CALL MPI_SEND(rho(:,ny-1),nx+4,mpireal,proc_y_max,tag,comm,errcode)
     END IF
     !If there are any other cycles the use the real equation of state
     eos_this_cycle=eos_number
  END DO

  by = 0.001_num

END SUBROUTINE set_initial_conditions


SUBROUTINE set_initial_conditions

  ! This is about the most complicated example for initial conditions
  ! used here as it covers including gravity and neutrals.
  ! The normalisation assumed is that from the defauls control.f90

  INTEGER :: loop
  INTEGER :: ix, iy, iy1
  REAL(num) :: a1, a2, dg
  REAL(num) :: a = 2.0_num, Tph = 9.8_num, Tcor = 980.0_num, ycor = 11.78_num
  REAL(num) :: betafs = 0.25_num, yfsl = -5.0_num, yfsu = 0.0_num
  REAL(num) :: wtr = 0.4_num, wfsl = 0.5_num, wfsu = 0.5_num
  REAL(num) :: r1, maxerr, xi_v
  REAL(num), DIMENSION(:), ALLOCATABLE :: yc_global, dyb_global, dyc_global
  REAL(num), DIMENSION(:), ALLOCATABLE :: grav_ref, temp_ref, rho_ref
  REAL(num), DIMENSION(:), ALLOCATABLE :: beta_ref, mag_ref, mu_m

  ALLOCATE( yc_global(-1:ny_global+1))
  ALLOCATE(dyb_global(-1:ny_global+1))
  ALLOCATE(dyc_global(-1:ny_global+1))
  ALLOCATE(  grav_ref(-1:ny_global+2))
  ALLOCATE(  temp_ref(-1:ny_global+2))
  ALLOCATE(   rho_ref(-1:ny_global+2))
  ALLOCATE(   mag_ref(-1:ny_global+2))
  ALLOCATE(  beta_ref(-1:ny_global+2))
  ALLOCATE(      mu_m(-1:ny_global+2))

  vx = 0.0_num
  vy = 0.0_num
  vz = 0.0_num
  bx = 0.0_num
  by = 0.0_num
  bz = 0.0_num

  ! Fill in yc_global with the positions central to the yb_global points
  DO iy = -1, ny_global + 1
    yc_global(iy) = 0.5_num * (yb_global(iy-1) + yb_global(iy))
  END DO

  ! Fill in dyb_global and dyc_global
  DO iy = -1, ny_global
    dyb_global(iy) = yb_global(iy) - yb_global(iy-1)
    dyc_global(iy) = yc_global(iy+1) - yc_global(iy)
  END DO

  ! Fill in the reference gravity array - lowering grav to zero at the top
  ! of the corona smoothly from a1 to grav = 0 at a2 and above
  grav_ref = 11.78_num
  a1 = yb_global(ny_global) - 20.0_num
  a2 = yb_global(ny_global) - 5.0_num
  DO iy = 0, ny_global + 2
    IF (yb_global(iy) > a1) THEN
      grav_ref(iy) = 11.78_num &
          * (1.0_num + COS(pi * (yb_global(iy) - a1) / (a2 - a1))) / 2.0_num
    END IF
    IF (yb_global(iy) > a2) THEN
      grav_ref(iy) = 0.0_num
    END IF
  END DO
  grav_ref(-1) = grav_ref(0)
  grav_ref(ny_global+1:ny_global+2) = grav_ref(ny_global)

  ! Beta profile from Archontis 2009 but in 2D
  ! Similar to that of Nozawa 1991
  ! NB : The variable beta used here is actually 1/beta
  beta_ref = 0.0_num
  DO iy = -1, ny_global + 1
    IF ((yc_global(iy) > yfsl) .AND. (yc_global(iy) < yfsu)) THEN
      beta_ref(iy) = betafs &
          * (0.5_num * (TANH((yc_global(iy) - yfsl) / wfsl) + 1.0_num)) &
          * (0.5_num * (1.0_num - TANH((yc_global(iy) - yfsu) / wfsu)))
    END IF
  END DO

  ! Calculate the density profile, starting from the refence density at the
  ! photosphere and calculating up and down from there including beta
  rho_ref = 1.0_num
  mu_m = 1.0_num
  IF (eos_number == EOS_IDEAL .AND. (.NOT. neutral_gas)) mu_m = 0.5_num

  DO loop = 1, 1000
    maxerr = 0.0_num
    ! Go from photosphere down
    DO iy = -1, ny_global + 1
      IF (yc_global(iy) < 0.0_num) THEN
        temp_ref(iy) = Tph - a * (gamma - 1.0_num) &
            * yc_global(iy) * grav_ref(iy) * mu_m(iy) / gamma
      END IF
      IF (yc_global(iy) >= 0.0_num) THEN
        a1 = 0.5_num * (TANH((yc_global(iy) - ycor) / wtr) + 1.0_num)
        temp_ref(iy) = Tph - 1.0 + (Tcor - Tph)**a1
      END IF
    END DO
    temp_ref(ny_global+1:ny_global+2) = temp_ref(ny_global)

    DO iy = ny_global, 0, -1
      IF (yc_global(iy) < 0.0_num) THEN
        iym = iy - 1
        dg = 1.0_num / (dyb_global(iy) + dyb_global(iym))

        rho_ref(iym) = rho_ref(iy ) * (temp_ref(iy ) &
            * (1.0_num + beta_ref(iy )) / dyc_global(iym) / mu_m(iy ) &
            + grav_ref(iym) * dyb_global(iy ) * dg)

        rho_ref(iym) = rho_ref(iym) / (temp_ref(iym) &
            * (1.0_num + beta_ref(iym)) / dyc_global(iym) / mu_m(iym) &
            - grav_ref(iym) * dyb_global(iym) * dg)
      END IF
    END DO

    ! Now move from the photosphere up to the corona
    DO iy = 0, ny_global
      IF (yc_global(iy) >= 0.0_num) THEN
        iym = iy - 1
        dg = 1.0_num / (dyb_global(iy) + dyb_global(iym))

        rho_ref(iy)  = rho_ref(iym) * (temp_ref(iym) &
            * (1.0_num + beta_ref(iym)) / dyc_global(iym) / mu_m(iym) &
            - grav_ref(iym) * dyb_global(iym) * dg)

        rho_ref(iy)  = rho_ref(iy ) / (temp_ref(iy ) &
            * (1.0_num + beta_ref(iy )) / dyc_global(iym) / mu_m(iy ) &
            + grav_ref(iym) * dyb_global(iy ) * dg)
      END IF
    END DO

    IF (eos_number /= EOS_IDEAL) THEN
      DO iy = 0, ny_global
        xi_v = get_neutral(temp_ref(iy), rho_ref(iy))
        r1 = mu_m(iy)
        mu_m(iy) = 1.0_num / (2.0_num - xi_v)
        maxerr = MAX(maxerr, ABS(mu_m(iy) - r1))
      END DO
    END IF

    IF (maxerr < 1.e-16_num) EXIT
  END DO

  rho_ref(ny_global+1:ny_global+2) = rho_ref(ny_global)

  ! Magnetic flux sheet profile from Archontis2009
  ! Similar structure to the 2D version used in Nozawa1991 and Isobe2006
  DO iy= -1, ny_global + 2
    mag_ref(iy) = SQRT(2.0_num * beta_ref(iy) * temp_ref(iy) * rho_ref(iy) &
        / mu_m(iy))
  END DO

  ! Fill in all the final arrays from the ref arrays
  iy1 = n_global_min(2) - 1

  DO iy = -1, ny + 2
    grav(iy) = grav_ref(iy1)
    DO ix = -1, nx + 2
      rho(ix,iy) = rho_ref(iy1)
      energy(ix,iy) = temp_ref(iy1)
      bx(ix,iy) = mag_ref(iy1)

      IF (eos_number /= EOS_IDEAL) THEN
        xi_v = get_neutral(energy(ix,iy), rho(ix,iy))
      ELSE
        IF (neutral_gas) THEN
          xi_v = 1.0_num
        ELSE
          xi_v = 0.0_num
        END IF
      END IF

      energy(ix,iy) = (energy(ix,iy) * (2.0_num - xi_v) &
          + (1.0_num - xi_v) * ionise_pot * (gamma - 1.0_num)) &
          / (gamma - 1.0_num)
    END DO
    iy1 = iy1 + 1
  END DO

  DO ix = -1, nx + 2
    energy(ix,ny+2) = energy(ix,ny+1)
  END DO

  DEALLOCATE(yc_global, dyb_global, dyc_global, mu_m)
  DEALLOCATE(grav_ref, temp_ref, rho_ref, beta_ref, mag_ref)

  CALL add_probe(0.0_num, 0.0_num)

END SUBROUTINE set_initial_conditions
