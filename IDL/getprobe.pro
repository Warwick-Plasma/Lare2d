FUNCTION getprobe,probe,wkdir=wkdir,single=single

  COMMON background, wkdir_global, retro_global
  IF NOT KEYWORD_SET(wkdir) THEN wkdir=wkdir_global
  file = wkdir + string(probe,format='("/probe",I03,".dat")')
  n=1ll
  close,10
  openr,10,file
  readu,10,n
  nvars=3
  IF NOT KEYWORD_SET(oldstyle) THEN nvars=7
  IF NOT KEYWORD_SET(single) THEN BEGIN
    data=dblarr(nvars,n)
    posx=1.0d
    posy=1.0d
  ENDIF ELSE BEGIN
    data=fltarr(nvars,n)
    posx=1.0
    posy=1.0
  ENDELSE
  readu,10,posx,posy
  readu,10,data
  close,10
  IF (KEYWORD_SET(oldstyle)) THEN BEGIN
    probe = {filename:file,time:reform(data(0,*)),x_pos:posx, y_pos:posy, vz:reform(data(1,*)), bz:reform(data(2,*))}
  ENDIF ELSE BEGIN
    probe = {filename:file,time:reform(data(0,*)),x_pos:posx, y_pos:posy, vx:reform(data(1,*)), vy:reform(data(2,*)), vz:reform(data(3,*)),bx:reform(data(4,*)), by:reform(data(5,*)), bz:reform(data(6,*))}
  ENDELSE
  return,probe
END
