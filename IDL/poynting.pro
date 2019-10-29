 pro poynting_flux, i_start, i_final

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
 	p0=b0^2/mu0             ;Pa
 	t0=l0/v0                ;s
 	j0=b0/(mu0*l0)          ;A/m
 	e0=v0*b0                ;V/m
 	energy0=v0^2            ;J/kg
 	temp0=energy0*mbar/kb   ;K
 	eta0=mu0*l0*v0          ;???

   wkdir = 'Data'
   a=getdata(wkdir,0)
   it = i_final-i_start+1
   nx = a.grid.npts[0]
   flux = dblarr(it, nx)

   for i = 0, i_final - i_start do begin
      a = getdata(wkdir,i + i_start)
      for ix = 0, nx - 2 do begin 
         flux(i,ix) = - a.by[ix,0] * (a.vz[ix,0] * a.bz[ix,0] + a.vx[ix,0] * a.bx[ix,0])
      end
   end
   flux = flux * v0 * b0^2 / mu0 * 1.e3 / 1.e7
   shade_surf,flux, ztitle='Poynting flux in 10^7 ergs/cm^2',charsize=3

   print, 'mean poynting flux in 10^7 erg/cm^2 = ', mean(flux) 
   print, 'peak poynting flux in 10^7 erg/cm^2 = ', max(flux) 

 end