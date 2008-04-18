

    SUBROUTINE equilibrium

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

  END SUBROUTINE equilibrium






SUBROUTINE equilibrium

    INTEGER:: ix, iy

    ! strong shock test
    ! standard test with nx=1200 and Lx=1
    ! run to t=0.038

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

  END SUBROUTINE equilibrium



  SUBROUTINE equilibrium

    INTEGER:: ix, iy

    ! 1d MHD Brio-Wu
    ! typically run with 800 points until t=0.1

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

  END SUBROUTINE equilibrium



  SUBROUTINE equilibrium

    INTEGER:: ix, iy
    REAL(num) :: b1, b2, x, y, s

    ! 2d MHD Brio-Wu
    ! needs open bcs. OK along diagonal with periodic bcs
    ! with Lx=Ly=1/SQRT(2) and t_end=0.1
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

  END SUBROUTINE equilibrium



  SUBROUTINE equilibrium

    INTEGER:: ix, iy
    ! Ortzang-Tang
    ! needs Lx=Ly=2*pi periodic bcs
    ! run to t=pi

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

  END SUBROUTINE equilibrium



  SUBROUTINE equilibrium

    INTEGER:: ix, iy
    REAL(num) :: angle
    REAL(num) :: k, a, cs, r, vpar, vperp, bpar, bperp

    ! 2d spherical Alfven wave
    ! set angle=pi/6 to get Alfven wave test in JCP
    ! assumes periodic bcs. Lx=1/cos(angle). Ly=1/sin(angle)

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

  END SUBROUTINE equilibrium



  SUBROUTINE equilibrium

    INTEGER:: ix, iy
    REAL(num) :: x0
    REAL(num) :: a, cs, r

    !  MHD pulse
    rho = 1.0_num
    energy = 1.0_num
    vx = 0.0_num
    vy = 0.0_num
    vz = 0.0_num
    bx = 1.0_num
    by = 0.0_num
    bz = 0.0_num
    a = 0.1_num
    x0 = 0.04_num

    cs = SQRT(gamma*(gamma-1.0_num))
    DO ix = 0, nx
       DO iy = 0, ny
          r = xb(ix)**2 !+ yb(iy)**2 
          vy(ix,iy) = a * EXP(-r/x0)/cs
       END DO
    END DO

  END SUBROUTINE equilibrium



  SUBROUTINE equilibrium

    INTEGER:: ix, iy, loop
    REAL(num):: b0, pressure
    REAL(num) :: b1, b2, x, y, s, s1, x0
    REAL(num) :: k, a, cs, r, vpar, vperp, bpar, bperp
    REAL(num) :: yt, q, bphi, pb
    REAL(num) :: ys, v0, lambda, a1, a2

    !  Fan stratified atmosphere
    IF (y_stretch) THEN
       yb = yb -20.0_num
       yc = yc -20.0_num
    ELSE
       yb = yb + (length_y/2.0_num) -20.0_num
       yc = yc + (length_y/2.0_num) -20.0_num
    END IF

    a1 = 40.0_num
    a2 = 60.0_num
    WHERE (yb > a1) grav = grav(0)*(1.0_num+COS(pi*(yb-a1)/(a2-a1)))/2.0_num
    WHERE (yb > a2) grav = 0.0_num

    CALL fan

    yt = -12.0_num
    b0 =  5.0_num
    a  =  4.0_num
    q  = -1.0_num/a
    DO ix = -1,nx+2
       DO iy = -1,ny+2
          r=SQRT(xc(ix)**2+(yc(iy)-yt)**2)
          bz(ix,iy) =  b0 * EXP(-(r/a)**2)
          bphi = bz(ix,iy) * q * r

          r=SQRT(xb(ix)**2+(yc(iy)-yt)**2)
          b1 =  b0 * EXP(-(r/a)**2)
          bx(ix,iy) = -b1 * q * (yc(iy)-yt) 

          r=SQRT(xc(ix)**2+(yb(iy)-yt)**2)
          b1 =  b0 * EXP(-(r/a)**2)
          by(ix,iy) =  b1 * q *  xc(ix)

          pb = (-0.25_num * (bz(ix,iy)**2)) - (0.5_num * (bphi**2))
          rho(ix,iy) = rho(ix,iy) + (pb/(energy(ix,iy)*(gamma-1.0_num)))
       END DO
    END DO

    bx(-2,:) = bx(-1,:)
    by(:,-2) = by(:,-1)
    vx(-2,:) = vx(-1,:)
    vy(:,-2) = vy(:,-1)

  END SUBROUTINE equilibrium



  SUBROUTINE fan !includes convection zone photosphere and corona

    REAL(num)::dg,m=1.5_num,a=1.1_num,y_cor=25.0_num,t_ph=1.0_num
    REAL(num) :: t_cor=150.0_num,wtr=5.0_num
    
    !define energy profile
    DO ix=-1,nx+2
       DO iy=-1,ny+2
          IF (yc(iy)<0.0_num) THEN
             energy(ix,iy) = t_ph - (a * grav(iy) * yc(iy) / (m+1.0_num) )
          ELSE 
             energy(ix,iy) = t_ph + &
                  ((t_cor-t_ph)*0.5_num * (TANH((yc(iy)-y_cor)/wtr)+1.0_num))
          END IF
       END DO
    END DO
    energy = energy / (gamma-1.0_num)

    !solve HS eqn to get rho profile
    rho = 50.0_num

    !wait for lower process
    IF (coordinates(1)/=0) THEN
       CALL MPI_RECV(rho(:,-1),nx+4,mpireal,down,tag,comm,status,errcode)
    END IF
    !calculate the density profile in y for the specified energy profile
    DO ix= -1,nx+2
       DO iy=0,ny+2
          dg=1.0_num/(dyb(iy)+dyb(iy-1))
          rho(ix,iy)=rho(ix,iy-1)*(energy(ix,iy-1)/dyc(iy-1)*(gamma-1.0_num)&
               -grav(iy-1)*dyb(iy-1)*dg)
          rho(ix,iy)=rho(ix,iy)/(energy(ix,iy)/dyc(iy-1)*(gamma-1.0_num)&
               +grav(iy-1)*dyb(iy)*dg)
       END DO
    END DO
    !now send to higher process
    IF (coordinates(1)/=nprocy-1) THEN
       CALL MPI_SEND(rho(:,ny-1),nx+4,mpireal,up,tag,comm,errcode)
    END IF

  END SUBROUTINE fan



  SUBROUTINE equilibrium

    INTEGER:: ix, iy, ivp, loop
    REAL(num) :: angle
    REAL(num):: b0, pressure
    REAL(num) :: b1, b2, x, y, s, s1, x0
    REAL(num) :: k, a, cs, r, vpar, vperp, bpar, bperp
    REAL(num) :: yt, q, bphi, pb
    REAL(num) :: ys, v0, lambda, a1, a2

    ! Shibata stratified atmosphere

    IF (y_stretch) THEN
       yb = yb - 20.0_num 
       yc = yc - 20.0_num
    ELSE
       yb = yb + (length_y/2.0_num) - 20.0_num
       yc = yc + (length_y/2.0_num) - 20.0_num
    END IF

    CALL fan

    !define flux sheet and hold in MEQ
    ys = -10.0_num
    a = 4.0_num
    b0 = 0.4_num
    DO ix=-1,nx+2
       DO iy=-1,ny+2
          r = yb(iy)-ys
          bx(ix,iy) = b0 * EXP(-(r/a)**2)
          pb = 0.5_num * (bx(ix,iy)**2)
          energy(ix,iy) = energy(ix,iy) - (pb/(rho(ix,iy)*(gamma-1.0_num)))
       END DO
    END DO

    !add non linear velocity pertubation
    v0=0.01_num
    lambda = 25.0_num
    DO ix=-1,nx+2
       DO iy=-1,ny+2
          r = yb(iy)-ys
          IF (ABS(xc(ix))<lambda/4.0_num) THEN
             vy(ix,iy) = v0 * cos(2.0_num*pi*xb(ix)/lambda) * EXP(-(r/a)**2)
          END IF
       END DO
    END DO

  END SUBROUTINE equilibrium




    
  SUBROUTINE equilibrium  !2D model flare loop interaction
    
    INTEGER :: ix, iy
    REAL(num) :: plasma_beta, total_pressure
    REAL(num) :: rho_e
    REAL(num) :: radius, rho_0
    REAL(num) :: loop_y, skin, gradrho, drho, signr
    REAL(num) :: flare_x, flare_y, a, r, r0

    !external values
    plasma_beta = 0.001_num
    rho_e = 1.0_num
    rho = rho_e
    vx = 0.0_num
    vy = 0.0_num
    vz = 0.0_num
    bx = 1.0_num
    by = 0.0_num
    bz = 0.0_num

    !calculate energy from plasma_beta
    energy = plasma_beta * (bx**2 + by**2 + bz**2) &
             / (2.0_num * rho * (gamma - 1.0_num))
    total_pressure = energy(0,0) * rho_e * (gamma - 1.0_num) + &
                     0.5_num * (bx(0,0)**2 + by(0,0)**2 + bz(0,0)**2)

    !coronal loop values
    loop_y = 0.0_num
    radius = 20.0_num
    skin = 0.0_num
    rho_0 = 100.0_num

    !create loop
    gradrho = 6.5_num / skin ! > atanh(0.99999) / skin
    drho = 0.5_num * (rho_0 - rho_e)
    DO iy = -1, ny+2
        r = yc(iy) - loop_y
        IF (r.LT.0.0_num) THEN
            signr = -1.0_num
        ELSE
            signr = 1.0_num
        END IF
        IF (ABS(r).GT.radius) THEN
            r0 = loop_y + signr * (radius + skin)
            rho(:,iy) = drho * tanh(gradrho*(-signr)*(yc(iy)-r0)) + drho + rho_e
        ELSE
            rho(:,iy) = rho_0
        END IF
    END DO

    !total pressure balance
    energy = (total_pressure - 0.5_num * (bx**2 + by**2 + bz**2)) &
             / (rho * (gamma - 1.0_num))

    !apply flare
    flare_x = 0.0_num
    flare_y = 20.0_num
    a = 10.0_num
    r0 = 10.0_num
    DO ix = -1, nx+2
        DO iy = -1, ny+2
            r = (xc(ix)-flare_x)**2 + (yc(iy)-flare_y)**2
            energy(ix,iy) = energy(ix,iy) + a * EXP(-r/r0)
        END DO
    END DO

  END SUBROUTINE equilibrium



  SUBROUTINE equilibrium

    ! NB: in this routine the variable energy is set to the temperature it
    ! is converted to specific internal energy density at the end

    INTEGER:: ix, iy
    REAL(num) :: x0, pg0
    REAL(num) :: a, r, b0

    rho = 1.e15_num * m_proton   ! n = 1.e15 m^(-3)
    energy = 2.e6_num            ! T = 2 MK
    vx = 0.0_num
    vy = 0.0_num
    vz = 0.0_num
    bx = 0.0_num                ! B = 100 G
    by = 0.0_num
    bz = 0.0_num
    gamma = 5.0_num / 3.0_num
    b0 = 1.737e-3_num
    pg0 =  1.e15_num * m_proton * 1.5e5_num**2 / gamma

    DO iy = -1, ny + 2
       bx(:,iy) = b0 * tanh(yc(iy) / 3.0e6_num)
       rho(:,iy) = 1.e15_num * m_proton + (mean_mass / 2.e6_num / k_b) * (b0**2 / 2.0_num / mu_0) * (1.0_num - tanh(yc(iy) / 3.e6_num)**2)
    END DO
    
    energy = energy * k_B / mean_mass / (gamma - 1.0_num)

  END SUBROUTINE equilibrium
