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
    ENDDO

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

    gamma = 5.0_num / 3.0_num

    vx = 0.0_num
    vy = 0.0_num
    vz = 0.0_num
    bx =  0.75_num
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
    IF (coordinates(1) .EQ. 0) grav(0) = 0.0_num
    IF (coordinates(1) .EQ. nprocy-1) grav(ny) = 0.0_num

    max_cycles=100
    IF (eos_number .EQ. EOS_IDEAL) max_cycles=1

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
        CALL MPI_RECV(rho(:,-1),nx+4,mpireal,down,tag,comm,status,errcode)
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
        CALL MPI_SEND(rho(:,ny-1),nx+4,mpireal,up,tag,comm,errcode)
     END IF
     !If there are any other cycles the use the real equation of state
     eos_this_cycle=eos_number
  ENDDO

  by = 0.001_num

END SUBROUTINE set_initial_conditions



  SUBROUTINE set_initial_conditions
  ! 2D flux emergence based on Fan atmosphere with flux tube
  ! as used (or at least close to) in Leake's 2D flux emergence paper
  !S.I. units
    INTEGER:: ix, iy, icycle,flipy

    REAL(num)::dg,m=1.5_num,a=1.0_num,y_cor=3.75e6,t_ph=6420.0_num
    REAL(num) :: t_cor=963000.0_num,wtr=7.5e5_num,T
    REAL(num) :: yt,q,r,bphi,b1,pb,b,b_0

    vx = 0.0_num
    vy = 0.0_num
    vz = 0.0_num
    bx = 0.0_num
    by = 0.0_num
    bz = 0.0_num
    energy = 1.0_num

    gamma = 5.0_num / 3.0_num

    grav = 274.0_num
    IF (coordinates(1) .EQ. 0) grav(0) = 0.0_num
    IF (coordinates(1) .EQ. nprocy-1) grav(ny) = 0.0_num

    !define energy profile
    DO ix = -1,nx+2
       DO iy = -1,ny+2
          IF (yc(iy) < 0.0_num) THEN
             T = t_ph - (a * grav(iy) * yc(iy) * mbar/kb/ (m+1.0_num) )
          ELSE 
             T = t_ph + &
                  ((t_cor-t_ph)*0.5_num * (TANH((yc(iy)-y_cor)/wtr)+1.0_num))
          END IF
          CALL Get_Energy(rho(ix,iy),T,EOS_IDEAL,ix,iy,energy(ix,iy))
       END DO
    END DO

    !solve HS eqn to get rho profile
    rho = 0.0135_num

    !wait for lower process
    IF (coordinates(1)/=0) THEN
       CALL MPI_RECV(rho(:,-1),nx+4,mpireal,down,tag,comm,status,errcode)
    END IF
    !calculate the density profile in y for the specified energy profile
    DO ix = -1,nx+2
       DO iy = 0,ny+2
          dg = 1.0_num/(dyb(iy)+dyb(iy-1))
          rho(ix,iy) = rho(ix,iy-1)*(energy(ix,iy-1)/dyc(iy-1)*(gamma-1.0_num)&
               -grav(iy-1)*dyb(iy-1)*dg)
          rho(ix,iy) = rho(ix,iy)/(energy(ix,iy)/dyc(iy-1)*(gamma-1.0_num)&
               +grav(iy-1)*dyb(iy)*dg)
       END DO
    END DO
    !now send to higher process
    IF (coordinates(1)/=nprocy-1) THEN
       CALL MPI_SEND(rho(:,ny-1),nx+4,mpireal,up,tag,comm,errcode)
    END IF

    yt = -12.0_num * 150.0e3_num
    b_0 =  5.0_num * 0.12_num
    a  =  2.0_num * 150.0e3_num
    q  = -1.0_num/a
    DO ix = -1,nx+2
       DO iy = -1,ny+2
          r=SQRT(xc(ix)**2+(yc(iy)-yt)**2)
          bz(ix,iy) =  b_0 * EXP(-(r/a)**2)
          bphi = bz(ix,iy) * q * r

          r=SQRT(xb(ix)**2+(yc(iy)-yt)**2)
          b1 =  b_0 * EXP(-(r/a)**2)
          bx(ix,iy) = -b1 * q * (yc(iy)-yt) 

          r=SQRT(xc(ix)**2+(yb(iy)-yt)**2)
          b1 =  b_0 * EXP(-(r/a)**2)
          by(ix,iy) =  b1 * q *  xc(ix)

          pb = (-0.25_num * (bz(ix,iy)**2))/mu0 - (0.5_num * (bphi**2))/mu0
          rho(ix,iy) = rho(ix,iy) + (pb/(energy(ix,iy)*(gamma-1.0_num)))
       END DO
    END DO

  END SUBROUTINE set_initial_conditions



  SUBROUTINE set_initial_conditions
  ! 2D flux emergence based on Fan atmosphere with flux tube
  ! as used (or at least close to) in Leake's 2D flux emergence paper
  ! but now with iteration to include effects of neutrals in equilibrium
  !S.I.units
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
    IF (coordinates(1) .EQ. 0) grav(0) = 0.0_num
    IF (coordinates(1) .EQ. nprocy-1) grav(ny) = 0.0_num

    max_cycles=100
    IF (eos_number .EQ. EOS_IDEAL) max_cycles=1

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

       !solve HS eqn to get rho profile
       rho = 0.0135_num

       !wait for lower process
       IF (coordinates(1)/=0) THEN
          CALL MPI_RECV(rho(:,-1),nx+4,mpireal,down,tag,comm,status,errcode)
       END IF
       !calculate the density profile in y for the specified energy profile
       DO ix = -1,nx+2
          DO iy = 0,ny+2
             dg = 1.0_num/(dyb(iy)+dyb(iy-1))
             rho(ix,iy) = rho(ix,iy-1)*(energy(ix,iy-1)/dyc(iy-1)*(gamma-1.0_num)&
                  -grav(iy-1)*dyb(iy-1)*dg)
             rho(ix,iy) = rho(ix,iy)/(energy(ix,iy)/dyc(iy-1)*(gamma-1.0_num)&
                  +grav(iy-1)*dyb(iy)*dg)
          END DO
       END DO
       !now send to higher process
       IF (coordinates(1)/=nprocy-1) THEN
          CALL MPI_SEND(rho(:,ny-1),nx+4,mpireal,up,tag,comm,errcode)
       END IF
       !If there are any other cycles the use the real equation of state
       eos_this_cycle=eos_number
    ENDDO

    yt = -12.0_num * 150.0e3_num
    b_0 =  5.0_num * 0.12_num
    a  =  2.0_num * 150.0e3_num
    q  = -1.0_num/a
    DO ix = -1,nx+2
       DO iy = -1,ny+2
          r=SQRT(xc(ix)**2+(yc(iy)-yt)**2)
          bz(ix,iy) =  b_0 * EXP(-(r/a)**2)
          bphi = bz(ix,iy) * q * r

          r=SQRT(xb(ix)**2+(yc(iy)-yt)**2)
          b1 =  b_0 * EXP(-(r/a)**2)
          bx(ix,iy) = -b1 * q * (yc(iy)-yt) 

          r=SQRT(xc(ix)**2+(yb(iy)-yt)**2)
          b1 =  b_0 * EXP(-(r/a)**2)
          by(ix,iy) =  b1 * q *  xc(ix)

          pb = (-0.25_num * (bz(ix,iy)**2))/mu0 - (0.5_num * (bphi**2))/mu0
          rho(ix,iy) = rho(ix,iy) + (pb/(energy(ix,iy)*(gamma-1.0_num)))
       END DO
    END DO

  END SUBROUTINE set_initial_conditions

