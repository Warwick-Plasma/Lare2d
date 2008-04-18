PRO addcurrents, data

on_error, 2

close,1

IF N_ELEMENTS(data) NE 0 THEN BEGIN
    nx = data.nx+1
    ny = data.ny+1
    prec = data.prec

    dy = shift(data.y,-1) - data.y
    dx = shift(data.x,-1) - data.x
    
    unit = (prec EQ 4) ? 1.0 : 1.0D
    
    adx = dx # replicate(unit,ny)
    ady = replicate(unit,nx) # dy
    
    jx = (shift(data.bz,0,-1) - data.bz) / ady
    jy = - (shift(data.bz,-1,0) - data.bz) / adx
    jz = (shift(data.by,-1,0) - data.by) / adx - (shift(data.bx,0,-1) - data.bx) / ady

; Sadly this gets the n boundary wrong so set to zero
    zero = (prec EQ 4) ? 0.0 : 0.0D
    jx(nx-1,*) = zero
    jx(*,ny-1) = zero
    jy(nx-1,*) = zero
    jy(*,ny-1) = zero
    jz(nx-1,*) = zero
    jz(*,ny-1) = zero
    
    curlb = SQRT(jx^2 + jy^2 + jz^2)
    data = CREATE_STRUCT(data,{jx: jx,jy: jy,jz: jz, curlb: curlb})
    
ENDIF ELSE BEGIN
    print, "Invalid variable passed"
    print, "Use: addcurrents, <data sturcture>"
ENDELSE


END
