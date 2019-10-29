pro getspectrum, wkdir_in

 
  rho0 = 1.67e-4
  l0 = 180.e3
  b0 = 0.03

  Q0 = 1.60217646d-19 ; proton charge [C]
  M0 = 9.10938188d-31 ; electron mass [kg]
  kb = 1.3806503d-23  ; Boltzmann's constant [J/K]
  epsilon0=8.8541878176203e-12
  mu0=4.0e-7*!dpi

  IF (N_ELEMENTS(mbar) EQ 0) THEN mbar=1.2*m0*1836.2

  v0=b0/SQRT(mu0*rho0)    ;m/s
  t0=l0/v0                ;s

  wkdir = ''

  IF (N_ELEMENTS(wkdir_in) NE 0) THEN BEGIN
    wkdir = wkdir_in
  ENDIF

  IF (wkdir EQ '') THEN wkdir = 'Data'

  file = wkdir + '/spectrum.dat'

  OPENR, lun, file, /GET_LUN ;, /F77_UNFORMATTED

  npts = 1000
  drive_axis = dblarr(npts)
  drive_amp = dblarr(npts)

  POINT_LUN, lun, 4
  READU, lun, drive_axis
  READU, lun, drive_amp

  pi = 3.1415926535
  plot, drive_axis[0:npts-1] / (2 * pi * t0), drive_amp[0:npts-1]  * v0, /xlog,/ylog $
    ,ytitle='Driver amplitude [m/s]' , xtitle='Frequency [Hz]'
 
  FREE_LUN, lun
  CLOSE, lun

  
END
