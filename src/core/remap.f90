MODULE remap

  USE shared_data
  USE xremap
  USE yremap
  USE zremap

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: eulerian_remap

CONTAINS

  ! Strang splitting

  SUBROUTINE eulerian_remap(i)

    INTEGER, INTENT(IN) :: i
    INTEGER :: case_test

    IF (rke) delta_ke = 0.0_num
    xpass = 1
    ypass = 1

    DO ix = -1, nx + 2
      DO iy = -1, ny + 2
        bx(ix,iy) = bx(ix,iy) * dyb(iy)
        by(ix,iy) = by(ix,iy) * dxb(ix)
        bz(ix,iy) = bz(ix,iy) * cv(ix,iy)
      END DO
    END DO

    DO iy = -1, ny + 2
      bx(-2,iy) = bx(-2,iy) * dyb(iy)
    END DO

    DO ix = -1, nx + 2
      by(ix,-2) = by(ix,-2) * dxb(ix)
    END DO

    case_test = MODULO(i, 6)

    ! Strang ordering
    SELECT CASE(case_test)
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

    DO iy = -1, ny + 2
      DO ix = -1, nx + 2
        bx(ix,iy) = bx(ix,iy) / dyb(iy)
        by(ix,iy) = by(ix,iy) / dxb(ix)
        bz(ix,iy) = bz(ix,iy) / cv(ix,iy)
      END DO
    END DO

    DO iy = -1, ny + 2
      bx(-2,iy) = bx(-2,iy) / dyb(iy)
    END DO

    DO ix = -1, nx + 2
      by(ix,-2) = by(ix,-2) / dxb(ix)
    END DO

  END SUBROUTINE eulerian_remap

END MODULE remap
