!**************************************************************
! This module contains the conversion factors between normalised
! and internal code units
!**************************************************************

MODULE normalise
  USE constants
  USE shared_data
  IMPLICIT NONE
  
CONTAINS

  
  SUBROUTINE normalise_transport
  
    REAL(num) :: eta_bar_0, pressure0, energy0, temp0, mbar, kappa0
  
    pressure0 = B0**2 / mu0_si ! Pressure 
    energy0 = B0**2 / (mu0_si * rho0)   
    mbar = mf * mh_si
    temp0 = mbar * pressure0 / (kb_si * rho0) ! Temperature in K
    
    ! Normalise tbar, r_bar and eta_bar for including Cowling resistivity and neutrals
    t_bar = t_bar / temp0
    r_bar = r_bar * rho0 / temp0**(3.0_num / 2.0_num) 
    eta_bar_0 = rho0**2 * SQRT(temp0) / B0**2
    eta_bar = eta_bar / eta_bar_0
  
    ! Normalise ionise_pot 
    ionise_pot = ionise_pot_si / (energy0 * mbar)
     
    ! normalise tr required for get_neutral etc.
    tr = tr / temp0                             
    
    ! normalise paralle thermal conductivity
    kappa0 = energy0**(3.0_num / 2.0_num) * rho0 * L0 &
        / (mbar / kb_si * energy0)**(7.0_num / 2.0_num)	
  	kappa_0 = 1.e-11_num / kappa0   
  	
  	!find the normalised temperature corresponding to 100MK
    temp0 = (mbar / kb_si) * B0**2 / (mu0_si * rho0) ! Temperature in K
  	temperature_100mk = 1.e8_num / temp0  
  
  END SUBROUTINE normalise_transport
  


END MODULE normalise
