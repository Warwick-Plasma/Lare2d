MODULE remap

  USE shared_data; USE xremap ;  USE yremap ; USE zremap
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: eulerian_remap

CONTAINS

  SUBROUTINE eulerian_remap(i)  ! Strang splitting

    INTEGER, INTENT(IN) :: i
    INTEGER :: case_test

    delta_ke = 0.0_num
    xpass = 1
    ypass = 1

    DO ix = -1, nx+2
       DO iy = -1, ny+2
          bx(ix,iy) = bx(ix,iy) * dyb(iy)  
          by(ix,iy) = by(ix,iy) * dxb(ix)
          bz(ix,iy) = bz(ix,iy) * cv(ix,iy)
       END DO
    END DO
    bx(-2,:) = bx(-2,:) * dyb
    by(:,-2) = by(:,-2) * dxb

    case_test = MODULO(i,6)

    SELECT CASE(case_test)  ! Strang ordering
    CASE (0)
       CALL remap_x
       CALL remap_y
       CALL remap_z
    CASE (1)
       CALL remap_y
       CALL remap_z
       CALL remap_x
    CASE (2)
       CALL remap_z
       CALL remap_x
       CALL remap_y
    CASE (3)
       CALL remap_x
       CALL remap_z
       CALL remap_y
    CASE (4)
       CALL remap_z
       CALL remap_y
       CALL remap_x
    CASE (5)
       CALL remap_y
       CALL remap_x
       CALL remap_z
    END SELECT

    DO ix = -1, nx+2
       DO iy = -1, ny+2
          bx(ix,iy) = bx(ix,iy) / dyb(iy)  
          by(ix,iy) = by(ix,iy) / dxb(ix)
          bz(ix,iy) = bz(ix,iy) / cv(ix,iy)
       END DO
    END DO
    bx(-2,:) = bx(-2,:) / dyb
    by(:,-2) = by(:,-2) / dxb

  END SUBROUTINE eulerian_remap

END MODULE remap
