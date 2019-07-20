module phys_constants

 real    :: rho_wat
       !|density of water [kg/m**3]
 real    :: rho_air
       !|density of air [kg/m**3]
 real    :: Cp_W
       !|specific heat of water [J/(kg*K)]  
 real    :: Cp_A
       !|specific heat of air at constant pressure [J/(kg*K)]
 real    :: Lwv 
       !|latent heat of vaporization [J/kg]
 real    :: sigma
       !|Stefan-Boltzmann constant [W/(m**2*K**4)] 
 real    :: eps
       !|water surface emissivity [n/d]
 real    :: eps_a
       !|atmoshere emissivity [n/d]
 real    :: x1
       !|Brent's formula first parameter [n/d]
 real    :: x2
       !|Brent's formula second parameter [n/d]
 real    :: Rd
       !|gas constant for dry air [J/(kg*K)]
 real    :: Rwv
       !|gas constant for water vapor [J/(kg*K)]
 real    :: R_univ
       !|universal gas constant [J/(mol*K)]
 real    :: Atm_P
       !|reference atmospheric pressure
 real    :: albedo_water
       !|albedo of water surface [n/d]
 real    :: S_const
       !|The solar constant [W/m2]
 real    :: pi
       !|The Pi number
 
 contains

 subroutine define_phys_constants

  implicit none

  rho_wat = 1000
  rho_air = 1.273
  Cp_W = 3990
  Cp_A = 1005
  Lwv = 2501000
  sigma = 5.670367E-08
  eps = 0.98 
  eps_a = 0.68
  x1 = 0.526
  x2 = 0.065
  Rd = 287.05
  Rwv = 461.91
  R_univ = 8.31441 
  Atm_P = 1.E+05 
  albedo_water = 0.07
  S_const = 1367 
  pi = 3.14

 end subroutine define_phys_constants

end module phys_constants
