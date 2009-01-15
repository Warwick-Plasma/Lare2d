MODULE zremap

  USE shared_data; USE boundary
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: remap_z

CONTAINS

  SUBROUTINE remap_z

    REAL(num) :: v_advect, flux1, flux2

    DO iy = 1, ny
       DO ix = 1, nx
          ixm = ix - 1
          iym = iy - 1
          v_advect = (vz1(ix,iy) + vz1(ix,iym)) / 2.0_num
          flux1 = v_advect * bx(ix,iy) * dt
          v_advect = (vz1(ixm,iy) + vz1(ixm,iym)) / 2.0_num
          flux2 = v_advect * bx(ixm,iy) * dt
          bz(ix,iy) = bz(ix,iy) + flux1 - flux2
       END DO
    END DO
    DO iy = 1, ny
       DO ix = 1, nx
          ixm = ix - 1
          iym = iy - 1
          v_advect = (vz1(ix,iy) + vz1(ixm,iy)) / 2.0_num
          flux1 = v_advect * by(ix,iy) * dt
          v_advect = (vz1(ix,iym) + vz1(ixm,iym)) / 2.0_num
          flux2 = v_advect * by(ix,iym) * dt
          bz(ix,iy) = bz(ix,iy) + flux1 - flux2
       END DO
    END DO

    CALL bz_bcs

  END SUBROUTINE remap_z

END MODULE zremap
