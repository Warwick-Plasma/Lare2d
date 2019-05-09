 pro figure4_minimal, i_start, i_final

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
   final = getdata(wkdir,i_final, /viscous_heat,/grid)
   start = getdata(wkdir,i_start, /viscous_heat,/grid)

   sz=size(final.viscous_heat)

   ; * 10 to get from J/m^3 to ergs/cm^3
   plot,final.grid.y*l0/1.0e3,total(final.viscous_heat - start.viscous_heat,1)/sz(1) * rho0 * energy0 / (final.time - start.time)/t0 * 10.0,xrange=[0,2500],yrange=[1e-4,10],/xsty,/ysty,/ylog
 

   wkdir = '../Run2/Data'
   final = getdata(wkdir,i_final, /viscous_heat,/grid)
   start = getdata(wkdir,i_start, /viscous_heat,/grid)

   sz=size(final.viscous_heat)

   ; * 10 to get from J/m^3 to ergs/cm^3
   oplot,final.grid.y*l0/1.0e3,total(final.viscous_heat - start.viscous_heat,1)/sz(1) * rho0 * energy0 / (final.time - start.time)/t0 * 10.0,l=2
 

   wkdir = '../Run3/Data'
   final = getdata(wkdir,i_final, /viscous_heat,/grid)
   start = getdata(wkdir,i_start, /viscous_heat,/grid)

   sz=size(final.viscous_heat)

   ; * 10 to get from J/m^3 to ergs/cm^3
   oplot,final.grid.y*l0/1.0e3,total(final.viscous_heat - start.viscous_heat,1)/sz(1) * rho0 * energy0 / (final.time - start.time)/t0 * 10.0,l=4


 end