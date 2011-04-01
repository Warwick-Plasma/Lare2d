MODULE initial_conditions

  USE shared_data
  USE neutral

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: set_initial_conditions

CONTAINS


  SUBROUTINE set_initial_conditions
    ! This is about the most complicated example for initial conditions
    ! Used here as it covers including gravity and neutrals
    ! The normalisation assumed is that from the defauls control.f90    

    INTEGER :: loop, i, j, k
    INTEGER :: ix, iy, iz
    REAL(num) :: a1, a2, blah, dg
    REAL(num) :: a=1.0_num, Tph=9.8_num, Tcor=980.0_num, ycor=11.0_num, wtr=0.6_num
    REAL(num) :: betafs=0.25_num, yfsl=-10.0_num, yfsu=-1.0_num, wfsl=0.5_num, wfsu=0.5_num
    REAL(num) :: r1, maxerr, xi_v, ran1
    REAL(num) :: amp, wptb, yptb1, yptb2, yptb, np
    REAL(num), DIMENSION(:), ALLOCATABLE :: yc_global, dyb_global, dyc_global
    REAL(num), DIMENSION(:), ALLOCATABLE :: grav_ref, temp_ref, rho_ref
    REAL(num), DIMENSION(:), ALLOCATABLE :: beta_ref, mag_ref, mu_m
    
    ALLOCATE(yc_global(-1:ny_global+1))
    ALLOCATE(dyb_global(-1:ny_global+1), dyc_global(-1:ny_global))
    ALLOCATE(grav_ref(-1:ny_global+2), temp_ref(-1:ny_global+2))
    ALLOCATE(rho_ref(-1:ny_global+2), mag_ref(-1:ny_global+2))
    ALLOCATE(beta_ref(-1:ny_global+2), mu_m(-1:ny_global+2))

    vx = 0.0_num
    vy = 0.0_num
    vz = 0.0_num
    bx = 0.0_num
    by = 0.0_num
    bz = 0.0_num

    !fill in yc_global with the positions central to the yb_global points
    DO iy = -1,ny_global+1
       yc_global(iy) = 0.5_num * (yb_global(iy-1) + yb_global(iy))
    END DO

    !fill in dyb_global and dyc_global
    DO iy = -1,ny_global
       dyb_global(iy) = yb_global(iy) - yb_global(iy-1)
       dyc_global(iy) = yc_global(iy+1) - yc_global(iy)
    END DO

    !fill in the reference gravity array - lowering grav to zero at the top 
    !of the corona smoothly from a1 to grav=0 at a2 and above
    grav_ref = 11.78_num
    a1 = yb_global(ny_global) - 20.0_num
    a2 = yb_global(ny_global) - 5.0_num
    DO iy = 0,ny_global+2
       IF (yb_global(iy) > a1) THEN
          grav_ref(iy) = 11.78_num * (1.0_num + COS(pi * (yb_global(iy) - a1) &
               / (a2-a1))) / 2.0_num
       END IF
       IF (yb_global(iy) > a2) THEN
          grav_ref(iy) = 0.0_num
       END IF
    END DO    
    grav_ref(-1) = grav_ref(0)
    grav_ref(ny_global+1:ny_global+2) = grav_ref(ny_global)

    !beta profile from Archontis 2009 but in 2D
    !similar to that of Nozawa 1991
    !NB : The variable beta used here is actually 1/beta
    beta_ref = 0.0_num
    DO iy = -1,ny_global+1
      IF ((yc_global(iy) .GT. yfsl) .AND. (yc_global(iy) .LT. yfsu)) THEN
        beta_ref(iy) = betafs * &
              (0.5_num * (TANH((yc_global(iy) - yfsl) / wfsl) + 1.0_num)) * &
              (0.5_num * (1.0_num - TANH((yc_global(iy) - yfsu) / wfsu)))
      END IF
    END DO

    !calculate the density profile, starting from the refence density at the
    !photosphere and calculating up and down from there including beta
    rho_ref = 1.0_num
    mu_m = 1.0_num
    IF (eos_number == EOS_IDEAL .AND. (.NOT. neutral_gas)) mu_m = 0.5_num
    DO loop = 1,1000
       maxerr = 0.0_num
       !Go from photosphere down
       DO iy = -1,ny_global+1
          IF (yc_global(iy) < 0.0_num) THEN
             temp_ref(iy) = Tph - a * (gamma - 1.0_num) &
                  * yc_global(iy) * grav_ref(iy) * mu_m(iy) / gamma 
          END IF
          IF (yc_global(iy) >= 0.0_num) THEN
             temp_ref(iy) = Tph + (Tcor - Tph) * 0.5_num &
                  * (TANH((yc_global(iy) - ycor) / wtr) + 1.0_num)
          END IF
       END DO
       temp_ref(ny_global+1:ny_global+2) = temp_ref(ny_global)
       
       DO iy = ny_global,0,-1
          IF (yc_global(iy) < 0.0_num) THEN  
             dg = 1.0_num / (dyb_global(iy) + dyb_global(iy-1))
             rho_ref(iy-1) = rho_ref(iy) * (temp_ref(iy)*(1.0_num+beta_ref(iy)) &
                  /dyc_global(iy-1)/mu_m(iy)+grav_ref(iy-1)*dyb_global(iy)*dg)
             rho_ref(iy-1) = rho_ref(iy-1) / (temp_ref(iy-1)*(1.0_num+beta_ref(iy-1)) &
                  /dyc_global(iy-1)/mu_m(iy-1)-grav_ref(iy-1)*dyb_global(iy-1)*dg)
          END IF
       END DO
       !Now move from the photosphere up to the corona
       DO iy = 0,ny_global
          IF (yc_global(iy) >= 0.0_num) THEN
             dg = 1.0_num / (dyb_global(iy)+dyb_global(iy-1))
             rho_ref(iy) = rho_ref(iy-1) * (temp_ref(iy-1)*(1.0_num+beta_ref(iy-1)) &
                  /dyc_global(iy-1)/mu_m(iy-1)-grav_ref(iy-1)*dyb_global(iy-1)*dg)
             rho_ref(iy) = rho_ref(iy) / (temp_ref(iy)*(1.0_num+beta_ref(iy)) &
                  /dyc_global(iy-1)/mu_m(iy)+grav_ref(iy-1)*dyb_global(iy)*dg)
          END IF
       END DO
       IF (eos_number /= EOS_IDEAL) THEN
          DO iy=0,ny_global,1
             xi_v = get_neutral(temp_ref(iy),rho_ref(iy),yb(iy))
             r1 = mu_m(iy)
             mu_m(iy) = 1.0_num / (2.0_num-xi_v)
             maxerr = MAX(maxerr, ABS(mu_m(iy) - r1))
          END DO
       END IF
       IF (maxerr < 1.e-16_num) EXIT
    END DO 

    rho_ref(ny_global+1:ny_global+2) = rho_ref(ny_global)
                                  
    !magnetic flux sheet profile from Archontis2009
    !similar structure to the 2D version used in Nozawa1991 and Isobe2006
    DO iy= -1,ny_global+2,1
       mag_ref(iy) = SQRT(2.0_num * beta_ref(iy) * temp_ref(iy) * rho_ref(iy) / mu_m(iy))
    END DO

    !fill in all the final arrays from the ref arrays
    grav(:) = grav_ref(coordinates(1)*ny-1:coordinates(1)*ny+ny+2)
    DO ix = -1,nx+2,1
       rho(ix,:) = rho_ref(coordinates(1)*ny-1:coordinates(1)*ny+ny+2)
       energy(ix,:) = temp_ref(coordinates(1)*ny-1:coordinates(1)*ny+ny+2)
       bx(ix,:) = mag_ref(coordinates(1)*ny-1:coordinates(1)*ny+ny+2)
    END DO
    DO ix = -1,nx+2,1
       DO iy = -1,ny+2,1 
         IF (eos_number /= EOS_IDEAL) THEN         
           xi_v = get_neutral(energy(ix,iy), rho(ix,iy), yb(iy))
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
    END DO
    DO ix=-1,nx+2,1
        energy(ix,ny+2) = energy(ix,ny+1)
    END DO

    amp = 0.01_num
    wptb = 20.0_num

    DO iy=1,ny,1
      IF ((yc_global(iy) .GT. -10.0_num) .AND. (yc_global(iy) .LT. -1.0_num)) THEN
        DO ix=1,nx,1
          vy(ix,iy) = (amp / 4.0_num) * COS(2.0_num*pi*xb(ix)/wptb) &
              * (TANH((yb(iy)-yfsl)/0.5_num)-TANH((yb(iy)-yfsu)/0.5_num))
        END DO
      END IF
    END DO  


    DEALLOCATE(yc_global, dyb_global, dyc_global, mu_m)
    DEALLOCATE(grav_ref, temp_ref, rho_ref, beta_ref, mag_ref)


  END SUBROUTINE set_initial_conditions

END MODULE initial_conditions
