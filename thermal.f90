module thermal

 use phys_constants, only: define_phys_constants 
 use surf_scheme1_mod

 !|contains only "river surface -- atmosphere" interface processes! 

!*                                                                   =
!*   THERMAL:                                                        =
!*   Purpose: calculating of shortwave and longwave radiation fluxes =
!*   Author: Ivan Korpushenkov                                       =
!*   Data last revised: 18.05.2018                                   =
!*                                                                   =
!*====================================================================
 
 contains

 real function TetaT(Ta, pressure_a) 
    !|for calculating of potential temperature                  =
    !*simple formula                                            =
    !* Input: Ta - air temperature [cel.]                       =
    !*                                                          =
    !*===========================================================
  
  implicit none
   
  real, intent(in) :: Ta, pressure_a

  TetaT = Ta * ((1000 / pressure_a)**0.286)

 end function TetaT
 
 real function U_2(U_10)
    !|for calculating of wind speed at z = 2 m                  =
    !* Input: U_10 -- wind speed at z=10 m                      =
    !*                                                          =
    !*===========================================================

  implicit none

  real, intent(in) :: U_10
  real :: z_0 = 1.E-04 

  U_2 = U_10 * (log(real(2) / z_0) / log(real(10) / z_0))

 end function U_2

 
 real function sin_h0(pi, lat, delta, HRA)
   !|for calculating of sun raditaion incoming                  = 
   !*
   !* Input: S_const - the solar constant                       =
   !*        lat - latitude                                     =
   !*        del - angle of declination of the Sun              =
   !*        alpha - hour angle of the sun                      =
   !*============================================================

  implicit none

  real, intent(in) :: pi, lat, delta, HRA

  sin_h0 = sin(lat*(pi/180))*sin(delta*(pi/180)       ) + &
  & cos(delta*(pi/180))*cos(lat*(pi/180))*cos(HRA*(pi/180))
 
 end function sin_h0


 real function S_wave(pi, S_constant, sin_h0, nCloud)

  implicit none

  real, intent(in) :: pi, S_constant, sin_h0, nCloud

  real, parameter :: C_sh = 0.5607
  real, parameter :: tet = 0.105

  real :: nClouds, h0, ep 
 
  h0 = asin(sin_h0)*(180/pi)
  ep = -0.0022*h0 + 0.27 
  nClouds = nCloud/10

  if ((asin(sin_h0)*(180/pi)) < 0) then
     S_wave = 0
  else
     S_wave = ((S_constant*sin_h0) / &
     & (1 + ((ep*tet) / sin_h0)))  * &
     & (1 -              C_sh*nClouds)
  endif

 end function S_wave

 real function S_h(water_wap, pressure)
    !|for calculating of a specific humidity                    =
    !* Input: water_wap -- water wapor pressure(hPa)            =
    !*                                                          =
    !*                                                          =
    !*===========================================================
  implicit none
  
  real, intent(in) :: water_wap, pressure

  S_h = (0.623*water_wap) / (pressure - 0.377*water_wap)

 end function S_h


 real function S_0(Tau, S, pressure)
    !|for calculating of a specific humidity at water surf.
    !* Input: Tau = T * S (m**2 K)
    !*        pressure - air pressure (hPa)
    !*        S - cross-section area, m**2
    !*===========================================================
  implicit none

  real, intent(in) :: Tau, S, pressure
  real :: T_s, water_wap, E
  real, parameter :: f_s = 0.9, E_0 = 6.107
  real, parameter :: a = 7.6326, b = 265.5

  T_s = Tau / S

  E = E_0*10**((a*T_s)/(b+T_s))
  water_wap = E*f_s

  S_0 = (0.623*water_wap) / (pressure - 0.377*water_wap)
 
 end function S_0

! real function Esum(Tau, S, Ta, water_wap, x1, x2, eps, eps_a,sigma)
 real function Esum(Ta, water_wap, x1, x2, eps, eps_a, sigma)
    !* for calculating longwave emissivity                      =
    !* Input: T - water temperature (in i-1 step of time)[K.]   =  
    !*        Ta - air temperature [K.]                         =
    !*        e - water vapor pressure [hPa]                    =
    !*===========================================================   

  implicit none
  
  real, intent(in) :: Ta, water_wap !, Tau, S
  real, intent(in) :: x1, x2, eps, eps_a, sigma

  real :: T_s

!  T_s = Tau / S

  !1. Brent's(1932) formula (for clear sky )
  Esum = sigma * eps_a * (Ta**4) * (x1 + x2*sqrt(water_wap)     )

  !2. Prata's(1996) formula (for clear sky)
  !Esum = (1 - (1 + a_5 * (water_wap / Ta))                   * &
  !& exp(-(b_5 + a_5 * c_5*(water_wap / Ta))**0.5))*sigma*(Ta**4)

  !3. Duerte's et al.(2006) formula (cloudly sky)
  !Esum = Esum * (1 - n**(0.671)) + 0.990 * (n**0.671)*sigma*(Ta**4)

 end function Esum

!  E_up = sigma*eps*(T_s**4)

 real function H_D(Cp_A, rho_air, C_h, U_2, pressure, Tau, S, T_2)
    !|aerodynamic formula for calculating of sensible heat      =      
    !*Input: Cp_A - specific heat of air                        =
    !*       rhoA - density of dry air                          =
    !*       bix2(9) - transfer coefficient (Monin-Obukhov)     =
    !*       U_2 - wind speed at z = 2 m.                       =
    !*       T - water surface temperature                      =
    !*       TetaT - potentinal temperature                     =
    !*                                                          =
    !*===========================================================

  implicit none

  real, intent(in) :: Cp_A, rho_air, U_2, & 
  & Tau, S,  T_2, pressure
  real :: Teta_0, Teta_2, T_s
  real(kind=8) :: C_h

  T_s = Tau / S
 
  Teta_0 = T_s*((1000/pressure)**0.286)
  Teta_2 = T_2*((1000/pressure)**0.286)
 
  H_D = Cp_A*rho_air*C_h*abs(U_2)*(Teta_2 - Teta_0)

 end function H_D


 real function L_E(rho_air, Lwv, C_e, U_2, Sp_humidity0, &
 &   Sp_humidity2)
    !|aerodynamic formula for calculating of latent heat        =      
    !* Input: rho_air -- density of air                         =
    !*        Lwv -- latent heat of vaporization                =                             
    !*        bix2(9) -- exchange coefficient (Monin-Obukhov)   = 
    !*        U_2 -- wind speed at z= 2 m.                      =
    !*        Sp_humidity0 -- specific humidity at z = 0m       =
    !*        Sp_humidity2 -- specific humidity at z = 2m       =
    !*===========================================================

  implicit none

  real, intent(in) :: rho_air, Lwv,  U_2, &
  & Sp_humidity0, Sp_humidity2
  
  real(kind=8) :: C_e

  L_E = rho_air*Lwv*C_e*abs(U_2)*(Sp_humidity2 - Sp_humidity0)


 end function L_E


end module
