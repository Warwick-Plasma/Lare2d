PRO addaz,data

on_error, 2

close,1

IF N_ELEMENTS(data) NE 0 THEN BEGIN

    nx=data.nx
    ny=data.ny

    dyb = shift(data.y,-1) - data.y
    dxb = shift(data.x,-1) - data.x

    unit = (data.prec EQ 4) ? 0.0 : 0.0D

    Az=data.bx

    Az(0,0) = unit

    FOR iy = 1, ny DO BEGIN
        Az(0,iy) = data.bx(0,iy) * dyb(iy) + Az(0,iy-1)
    END

    FOR ix = 1, nx DO BEGIN
        Az(ix,0) = Az(ix-1,0) - data.by(ix,0) * dxb(ix) 
    END

    FOR iy = 1, ny DO BEGIN
        FOR ix = 1, nx DO BEGIN
            Az(ix,iy) = data.bx(ix,iy) * dyb(iy) + Az(ix,iy-1)
        END
    END

data=CREATE_STRUCT(data,{az: az})

ENDIF ELSE BEGIN
    print, "Invalid variable passed"
    print, "Use: addcurrents, <data sturcture>"
ENDELSE

END


