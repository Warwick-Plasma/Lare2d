MODULE openboundary

  USE shared_data

  IMPLICIT NONE

  REAL(num) :: pfar, rhofar, uxfar, uyfar, uzfar
  REAL(num) :: bxfar, byfar, bzfar, pbc, vnorm, fraction

  REAL(num), DIMENSION(0:1) :: vxbc, vybc, vzbc
  REAL(num), DIMENSION(0:1) :: bxbc, bybc, bzbc
  REAL(num), DIMENSION(0:1) :: rbc, ebc

CONTAINS



  SUBROUTINE open_bcs

    REAL(num) :: bperp

    ! update ghost cells based on Riemann problem with farfield
    ! only expected to work perfectly for prblems with straight B field
    ! through boundaries which do not drastically change shape during the
    ! simulation.    
    
    fraction = 0.9_num

    ! right boundary
    IF (xbc_right == BC_OPEN .AND. right == MPI_PROC_NULL) THEN
      DO iy = -1, ny + 1
        ! variables carried out of domain by Riemann invariants
        vxbc(1) = vx(nx, iy)
        vybc(1) = vy(nx, iy)
        vzbc(1) = vz(nx, iy)
        bxbc(1) = bx(nx, iy)
        bybc(1) = by(nx, iy)
        bzbc(1) = bz(nx, iy)
        rbc(1) = rho(nx, iy)
        ebc(1) = energy(nx, iy)

        pbc = (gamma - 1.0_num) * energy(nx, iy) * rho(nx, iy)

        ! farfield values carried into domain
        pfar = (gamma - 1.0_num) * energy(nx+2, iy) * rho(nx+2, iy)

        rhofar = rho(nx+2, iy)
        uxfar = vx(nx+2, iy)
        uyfar = vy(nx+2, iy)
        uzfar = vz(nx+2, iy)
        bxfar = bx(nx+2, iy)
        byfar = by(nx+2, iy)
        bzfar = bz(nx+2, iy)

        vnorm = vx(nx, iy)
        
        ! select correct open bc solver
        bperp = SQRT(byfar**2 + bzfar**2)
        IF (ABS(bxfar) <= fraction * bperp) THEN
          CALL open_bcs_1
        ELSE 
          CALL open_bcs_2
        END IF

        rho(nx+1, iy) = rbc(0)
        energy(nx+1, iy) = ebc(0)
        bz(nx+1, iy) = bzbc(0)
        by(nx+1, iy) = bybc(0)
        bx(nx+1, iy) = bxbc(0)
        vx(nx+1, iy) = vxbc(0)
        vy(nx+1, iy) = vybc(0)
        vz(nx+1, iy) = vzbc(0)
      END DO
    END IF

    ! left bounday
    IF (xbc_left == BC_OPEN .AND. left == MPI_PROC_NULL) THEN
      DO iy = -1, ny + 1
        vxbc(1) = - vx(0, iy)
        vybc(1) = vy(0, iy)
        vzbc(1) = vz(0, iy)
        bxbc(1) = - bx(0, iy)
        bybc(1) = by(1, iy)
        bzbc(1) = bz(1, iy)
        rbc(1) = rho(1, iy)
        ebc(1) = energy(1, iy)

        pbc = (gamma - 1.0_num) * energy(1, iy) * rho(1, iy)

        pfar = (gamma - 1.0_num) * energy(-1, iy) * rho(-1, iy)

        rhofar = rho(-1, iy)
        uxfar = - vx(-2, iy)
        uyfar = vy(-2, iy)
        uzfar = vz(-2, iy)
        bxfar = - bx(-2, iy)
        byfar = by(-1, iy)
        bzfar = bz(-1, iy)

        vnorm = - vx(0, iy)
        bperp = SQRT(byfar**2 + bzfar**2)

        IF (ABS(bxfar) <= fraction * bperp) THEN
          CALL open_bcs_1
        ELSE 
          CALL open_bcs_2
        END IF            

        rho(0, iy) = rbc(0)
        energy(0, iy) = ebc(0)
        bz(0, iy) = bzbc(0)
        bx(-1, iy) = - bxbc(0)
        vx(-1, iy) = - vxbc(0)
        vy(-1, iy) = vybc(0)
        vz(-1, iy) = vzbc(0)
        by(0, iy) = bybc(0)
      END DO
    END IF

    ! top boundary
    IF (ybc_up == BC_OPEN .AND. up == MPI_PROC_NULL) THEN
      DO ix = -1, nx + 1
        vxbc(1) = vy(ix, ny)
        vybc(1) = vx(ix, ny)
        vzbc(1) = vz(ix, ny)
        bxbc(1) = by(ix, ny)
        bybc(1) = bx(ix, ny)
        bzbc(1) = bz(ix, ny)
        rbc(1) = rho(ix, ny)
        ebc(1) = energy(ix, ny)

        pbc = (gamma - 1.0_num) * energy(ix, ny) * rho(ix, ny)

        pfar = (gamma - 1.0_num) * energy(ix, ny+2) * rho(ix, ny+2)

        rhofar = rho(ix, ny+2)
        uxfar = vy(ix, ny+2)
        uyfar = vx(ix, ny+2)
        uzfar = vz(ix, ny+2)
        bxfar = by(ix, ny+2)
        byfar = bx(ix, ny+2)
        bzfar = bz(ix, ny+2)

        vnorm = vy(ix, ny)
        bperp = SQRT(byfar**2 + bzfar**2)

        IF (ABS(bxfar) <= fraction * bperp) THEN
          CALL open_bcs_1
        ELSE 
          CALL open_bcs_2
        END IF  
        
        rho(ix, ny+1) = rbc(0)
        energy(ix, ny+1) = ebc(0)
        bz(ix, ny+1) = bzbc(0)
        by(ix, ny+1) = bxbc(0)
        vx(ix, ny+1) = vybc(0)
        vy(ix, ny+1) = vxbc(0)
        vz(ix, ny+1) = vzbc(0)
        bx(ix, ny+1) = bybc(0)
      END DO
    END IF

    ! bottom boundary
    IF (ybc_down == BC_OPEN .AND. down == MPI_PROC_NULL) THEN
      DO ix = -1, nx + 1
        vxbc(1) = - vy(ix, 0)
        vybc(1) = vx(ix, 0)
        vzbc(1) = vz(ix, 0)
        bxbc(1) = - by(ix, 0)
        bybc(1) = bx(ix, 1)
        bzbc(1) = bz(ix, 1)
        rbc(1) = rho(ix, 1)
        ebc(1) = energy(ix, 1)

        pbc = (gamma - 1.0_num) * energy(ix, 1) * rho(ix, 1)

        pfar = (gamma - 1.0_num) * energy(ix, -1) * rho(ix, -1)

        rhofar = rho(ix, -1)
        uxfar = - vy(ix, -2)
        uyfar = vx(ix, -2)
        uzfar = vz(ix, -2)
        bxfar = - by(ix, -2)
        byfar = bx(ix, -1)
        bzfar = bz(ix, -1)

        vnorm = -vy(ix, 0)
        bperp = SQRT(byfar**2 + bzfar**2)

        IF (ABS(bxfar) <= fraction * bperp) THEN
          CALL open_bcs_1
        ELSE 
          CALL open_bcs_2
        END IF

        rho(ix, 0) = rbc(0)
        energy(ix, 0) = ebc(0)
        bz(ix, 0) = bzbc(0)
        by(ix, -1) = - bxbc(0)
        vx(ix, -1) = vybc(0)
        vy(ix, -1) = - vxbc(0)
        vz(ix, -1) = vzbc(0)
        bx(ix, 0) = bybc(0)
      END DO
    END IF

  END SUBROUTINE open_bcs



  SUBROUTINE open_bcs_1

    ! open bc when bx = 0
    REAL(num) :: c0, ct, cf
    REAL(num) :: pg, rhog, cffar, c0far, ctfar
    REAL(num) :: pmagg
    REAL(num), DIMENSION(3) :: vtest, pstar, vstar, rhostar, pmagstar, bystar, bzstar
    INTEGER :: i

    c0far = SQRT(gamma * pfar / rhofar)
    ctfar = SQRT((byfar**2 + bzfar**2) / rhofar)
    cffar = SQRT(c0far**2 + ctfar**2)   

    c0 = SQRT(gamma * (gamma - 1.0_num) * ebc(1))
    ct = SQRT((bybc(1)**2 + bzbc(1)**2) / rbc(1))
    cf = SQRT(c0**2 + ct**2)

    ! Define the speeds of the characteristics to be checked along
    vtest(1) = vnorm + cf
    vtest(2) = vnorm - cf
    vtest(3) = vnorm

    DO i = 1, 3
      IF (vtest(i) >= 0.0_num) THEN
        pstar(i) =  pbc + 0.5_num * (bybc(1)**2 + bzbc(1)**2)
        pmagstar(i) =  0.5_num * (bybc(1)**2 + bzbc(1)**2)
        vstar(i) = vxbc(1)
        rhostar(i) = rbc(1) 
        bystar(i) = bybc(1)
        bzstar(i) = bzbc(1)
      ELSE
        pstar(i) = pfar + 0.5_num * (byfar**2 + bzfar**2)
        pmagstar(i) = 0.5_num * (byfar**2 + bzfar**2)
        vstar(i) = uxfar
        rhostar(i) = rhofar
        bystar(i) = byfar
        bzstar(i) = bzfar
      END IF
    END DO

    bxbc(0) = bxbc(1)
    bybc(0) = 0.5_num * (bystar(1) + bystar(2))    
    bzbc(0) = 0.5_num * (bzstar(1) + bzstar(2)) 

    pg = 0.5_num &
        * (pstar(1) + pstar(2) + rhofar * cffar * (vstar(1) - vstar(2)))
    pmagg = 0.5_num * (bybc(0)**2 + bzbc(0)**2)  
    
    rhog = (pg - pmagg - pstar(3) + pmagstar(3)) / cffar**2 + rhostar(3)
    rbc(0) = MAX(rhog, none_zero)    

    ebc(0) = MAX(pg - pmagg, none_zero) / (gamma - 1.0_num) / rbc(0)  
    
    vxbc(0) = 0.5_num &
        * (vstar(1) + vstar(2) + (pstar(1) - pstar(2)) / (rhofar * cffar))
    vybc(0) = vybc(1)
    vzbc(0) = vzbc(1)  

  END SUBROUTINE open_bcs_1



  SUBROUTINE open_bcs_2

    ! open bc when bperp = 0
    REAL(num) :: lambdayfar, lambdazfar
    REAL(num) :: c0, cx
    REAL(num) :: pg, rhog, c0far, cxfar
    REAL(num) :: pmagg, lambdag, beta
    REAL(num), DIMENSION(5) :: vtest, pstar, uxstar, rhostar, pmagstar
    REAL(num), DIMENSION(5) :: uystar, lambdaystar, lambdazstar, uzstar
    INTEGER :: i

    lambdayfar = -bxfar * byfar
    lambdazfar = -bxfar * bzfar
    c0far = SQRT(gamma * pfar / rhofar)
    cxfar = SQRT(bxfar**2 / rhofar)  
    beta = (c0far / cxfar)**2

    c0 = SQRT(gamma * (gamma - 1.0_num) * ebc(1))
    cx = SQRT(bxbc(1)**2 / rbc(1))

    ! Define the speeds of the characteristics to be checked along
    vtest(1) = vnorm + c0
    vtest(2) = vnorm - c0
    vtest(3) = vnorm + cx
    vtest(4) = vnorm - cx
    vtest(5) = vnorm

    DO i = 1, 5
      IF (vtest(i) >= 0.0_num) THEN
        pstar(i) = pbc + 0.5_num * (bybc(1)**2 + bzbc(1)**2 - bxbc(1)**2)
        pmagstar(i) = 0.5_num * (bybc(1)**2 + bzbc(1)**2 - bxbc(1)**2)
        lambdaystar(i) = -bybc(1) * bxbc(1)
        lambdazstar(i) = -bzbc(1) * bxbc(1)
        uxstar(i) = vxbc(1)
        uystar(i) = vybc(1)
        uzstar(i) = vzbc(1)
        rhostar(i) = rbc(1)
      ELSE
        pstar(i) = pfar + 0.5_num * (byfar**2 + bzfar**2 - bxfar**2)
        pmagstar(i) = 0.5_num * (byfar**2 + bzfar**2 - bxfar**2)
        lambdaystar(i) = lambdayfar
        lambdazstar(i) = lambdazfar
        uxstar(i) = uxfar
        uystar(i) = uyfar
        uzstar(i) = uzfar
        rhostar(i) = rhofar
      END IF
    END DO

    bxbc(0) = bxbc(1)  
    
    lambdag = 0.5_num * (lambdaystar(3) + lambdaystar(4) &
        + rhofar * cxfar * (uystar(3) - uystar(4)))
    bybc(0) = -lambdag / bxbc(0)                   
    
    lambdag = 0.5_num * (lambdazstar(3) + lambdazstar(4) &
        + rhofar * cxfar * (uzstar(3) - uzstar(4)))
    bzbc(0) = -lambdag  / bxbc(0)                  
    
 
    IF (beta > 0.01_num .AND. beta < 10.0_num) THEN
       pmagg = 0.5_num * (pmagstar(1) + pmagstar(2))
       pg = 0.5_num * (pstar(1) + pstar(2) + rhofar * c0far * (uxstar(1) - uxstar(2)))
  
       rhog = (ABS(pg - pmagg) - ABS(pstar(5) - pmagstar(5))) / c0far**2 + rhostar(5) 
  
       rbc(0) = MAX(rhog, none_zero)
       ebc(0) = MAX(pg - pmagg, none_zero) / (gamma - 1.0_num) / rbc(0)
  
       vxbc(0) = 0.5_num &
              * (uxstar(1) + uxstar(2) + (pstar(1) - pstar(2)) / (rhofar * c0far))
       vybc(0) = 0.5_num * (uystar(3) + uystar(4) &
              + (lambdaystar(3) - lambdaystar(4)) / (rhofar * cxfar))
       vzbc(0) = 0.5_num * (uzstar(3) + uzstar(4) &
              + (lambdazstar(3) - lambdazstar(4)) / (rhofar * cxfar))  
    ELSE IF (beta > 10.0_num) THEN
       lambdag = 0.5_num * (lambdaystar(3) + lambdaystar(4))
       bybc(0) = -lambdag / bxbc(0)                   
    
       lambdag = 0.5_num * (lambdazstar(3) + lambdazstar(4))
       bzbc(0) = -lambdag  / bxbc(0)
       pmagg = 0.5_num * (pmagstar(1) + pmagstar(2))
       pg = 0.5_num * (pstar(1) + pstar(2) + rhofar * c0far * (uxstar(1) - uxstar(2)))
  
       rhog = (ABS(pg - pmagg) - ABS(pstar(5) - pmagstar(5))) / c0far**2 + rhostar(5) 
  
       rbc(0) = MAX(rhog, none_zero)
       ebc(0) = MAX(pg - pmagg, none_zero) / (gamma - 1.0_num) / rbc(0)
  
       vxbc(0) = 0.5_num &
              * (uxstar(1) + uxstar(2) + (pstar(1) - pstar(2)) / (rhofar * c0far))
       vybc(0) = 0.5_num * (uystar(3) + uystar(4)) 
       vzbc(0) = 0.5_num * (uzstar(3) + uzstar(4)) 
    ELSE 
       pmagg = 0.5_num * (pmagstar(1) + pmagstar(2))
       pg = 0.5_num * (pstar(1) + pstar(2)) 
  
       rhog = rhostar(5) 
  
       rbc(0) = MAX(rhog, none_zero)
       ebc(0) = MAX(pg - pmagg, none_zero) / (gamma - 1.0_num) / rbc(0)
  
       vxbc(0) = 0.5_num * (uxstar(1) + uxstar(2))
       vybc(0) = 0.5_num * (uystar(3) + uystar(4) &
              + (lambdaystar(3) - lambdaystar(4)) / (rhofar * cxfar))
       vzbc(0) = 0.5_num * (uzstar(3) + uzstar(4) &
              + (lambdazstar(3) - lambdazstar(4)) / (rhofar * cxfar))   
    ENDIF   


  END SUBROUTINE open_bcs_2



! This routine is currently not used and is here for later extensions and
! improvements to the open boundary conditions if they are needed.
  SUBROUTINE open_bcs_3
    ! Solve for when bx and bperp are non zero. Solves in the coordinate system
    ! such that y-axis points along by_farfield  
    REAL(num), DIMENSION(7) :: vtest
    INTEGER :: i
    REAL(num) :: a, b, c, d, e, f, g
    REAL(num) :: pmagg, pmagfar, beta
    REAL(num) :: c0, cx, ct, cf, cs
    REAL(num) :: lambdafar, byfar2
    REAL(num) :: c0far, cxfar, ctfar, cffar, csfar
    REAL(num) :: pg, rhog, uxg, uyg, uzg, lambdag, byg, bxg, bzg
    REAL(num), DIMENSION(7) :: pstar, uxstar, uystar, uzstar, rhostar
    REAL(num), DIMENSION(7) :: lambdastar, pmagstar, bzstar

    ! Setup the far field variables
    c0far = SQRT(gamma * pfar / rhofar)
    pmagfar = 0.5_num * (byfar**2 + bzfar**2 - bxfar**2)
    pfar = pfar + pmagfar
    lambdafar = -bxfar * byfar
    cxfar = ABS(bxfar) / SQRT(rhofar)
    ctfar = SQRT((byfar**2 + bzfar**2) / rhofar)
    a =  c0far**2 + cxfar**2 + ctfar**2
    b =  4.0_num * c0far**2 * cxfar**2
    cffar = SQRT(0.5_num * (a + SQRT(a**2 - b)))
    csfar = SQRT(0.5_num * (a - SQRT(a**2 - b))) 
    beta = c0far**2 / (cxfar**2 + ctfar**2)

    ! Setup the speeds
    c0 = SQRT(gamma * (gamma - 1.0_num) * ebc(1))
    cx = ABS(bxbc(1)) / SQRT(rbc(1))
    ct = SQRT((bybc(1)**2 + bzbc(1)**2) / rbc(1))
    a =  c0**2 + cx**2 + ct**2
    b =  4.0_num * c0**2 * cx**2
    cf = SQRT(0.5_num * (a + SQRT(a**2 - b)))
    cs = SQRT(0.5_num * (a - SQRT(a**2 - b)))

    ! Define the speeds of the characteristics to be checked along
    vtest(1) = vnorm + cf
    vtest(2) = vnorm - cf
    vtest(3) = vnorm - cs
    vtest(4) = vnorm + cs
    vtest(5) = vnorm
    vtest(6) = vnorm + cx
    vtest(7) = vnorm - cx

    ! Now check which characteristics are inflowing, outflowing, non-moving
    DO i = 1, 7
      IF (vtest(i)  >= 0.0_num) THEN
        pstar(i) = pbc + 0.5_num * (bybc(1)**2 + bzbc(1)**2 - bxbc(1)**2)
        pmagstar(i) = 0.5_num * (bybc(1)**2 + bzbc(1)**2 - bxbc(1)**2)
        rhostar(i) = rbc(1)
        lambdastar(i) = -bxbc(1) * bybc(1) 
        bzstar(i) = bzbc(1)
        uystar(i) = vybc(1) 
        uzstar(i) = vzbc(1)
        uxstar(i) = vxbc(1)
      ELSE
        pstar(i) = pfar
        pmagstar(i) = pmagfar
        rhostar(i) = rhofar
        lambdastar(i) = lambdafar
        bzstar(i) = bzfar
        uystar(i) = uyfar
        uzstar(i) = uzfar
        uxstar(i) = uxfar
      END IF

    END DO

    ! Now setup the constants that are defined in the solution
    a = cffar**2 - cxfar**2
    b = (lambdafar / rhofar)
    c = csfar**2 - cxfar**2

    d = pstar(1) + pstar(2) + rhofar * cffar * (uxstar(1) - uxstar(2))
    e = lambdastar(1) + lambdastar(2) + rhofar * cffar * (uystar(1) - uystar(2))
    f = pstar(3) + pstar(4) - rhofar * csfar * (uxstar(3) - uxstar(4))
    g = lambdastar(3) + lambdastar(4) - rhofar * csfar * (uystar(3) - uystar(4))
    pg = 0.5_num * (a*d + b*e - c*f - b*g) / (a - c)
    lambdag = 0.5_num * (c * (a*d + b*e) - a * (c*f + b*g)) / (b * (c - a))

    d = (pstar(1) - pstar(2)) / (rhofar * cffar) + (uxstar(1) + uxstar(2))
    e = (lambdastar(1) - lambdastar(2)) &
        / (rhofar * cffar) + (uystar(1) + uystar(2))
    f = (pstar(4) - pstar(3)) / (rhofar * csfar) + (uxstar(3) + uxstar(4))
    g = (lambdastar(4) - lambdastar(3)) &
        / (rhofar * csfar) + (uystar(3) + uystar(4))
    uxg = 0.5_num * (a*d + b*e - c*f - b*g) / (a - c)
    uyg = 0.5_num * (c * (a*d + b*e) - a * (c*f + b*g)) / (b * (c - a))

    a = cxfar * rhofar / bxfar
    bzg = 0.0_num
    uzg = 0.0_num

    bxg = bxbc(1)
    byg = -lambdag / bxg

    IF (beta > 0.01_num) THEN
      pmagg = 0.5_num * (byg**2 + bzg**2 - bxg**2)

      rhog = (ABS(pg - pmagg) - ABS(pstar(5) - pmagstar(5))) / c0**2 + rhostar(5)
      rhog = MAX(rhog, none_zero)
      rbc(0) = rhog
      ebc(0) = MAX(pg - pmagg, none_zero) / ((gamma - 1.0_num) * rhog)
      vxbc(0) = uxg
      vybc(0) = uyg 
    ELSE
      rbc(0) = rhostar(5)
      ebc(0) = ebc(1)
      vxbc(0) = 0.5_num * (uxstar(1) + uxstar(2))
      vybc(0) = 0.5_num * (uystar(1) + uystar(2))
    ENDIF

    ! rotate back to grid coordinate system
    bxbc(0) = bxg
    bybc(0) = byg 
    bzbc(0) = bzg

  END SUBROUTINE open_bcs_3

END MODULE openboundary
