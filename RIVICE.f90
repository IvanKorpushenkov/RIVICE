subroutine rivice(ntime, temp_air, temp_water, solar_wave, thickness_ice)
        !|ntime - number of time steps
        !|temp_air - air temperature [K]
        !|temp_water - water temperature [K]
        !|solar_wave - solar radiation [W/m**2]
        !|thickness_ice0 - thickness of river ice before [m]
        !|thickness_ice - actually thickness of river ice [m]

 implicit none

 !#in
 integer, intent(in) :: ntime
 real :: temp_air(1:ntime)
 real :: temp_water(1:ntime)
 real :: solar_wave(1:ntime)

 !#local
 integer, parameter :: time_step = 1
          !|time step [s]
 integer :: i
          !loop index
 real :: temp_ice(1:ntime)
       !|temp_ice - river ice body temperature [K]
 real :: temp_botice
       !|temp_botice - river ice bottom temperature [K]

 !#out
 real, intent(out) :: thickness_ice(1:ntime)

 !#phys_constants
 real, parameter :: CONDUCTIVITY_ICE  = 2.23
                  !|Thermal conductivity of river ice [W/m*K]
 real, parameter :: DENSITY_RIVER_ICE = 916.8
                  !|Density of river ice [kg/m**3]
 real, parameter :: HEAT_FUSION_ICE = 330*(10**3)
                  !|[joule/kg]

 !#simlification
 temp_botice = 0

 !#ordinary_differential_equation_of_river_ice_thickness_dynamics
 !#Donchenko R. V. Ice regime of the USSR rivers, 1987
 !#eq.4.3
 !#Euler_method

 !#first_step
 temp_ice(1) = min(temp_air(1), 0.0)
 thickness_ice(1) = 1.E-03 + (time_step * HEAT_FUSION_ICE**(-1) * DENSITY_RIVER_ICE**(-1)) * &
         (temp_ice(1) * CONDUCTIVITY_ICE * 1.E-03**(-1) - solar_wave(1))

 do i = 2, ntime
 temp_ice(i-1) = min(temp_air(i-1), temp_water(i-1))
 thickness_ice(i) = thickness_ice(i-1) + (time_step * HEAT_FUSION_ICE**(-1) * DENSITY_RIVER_ICE**(-1)) * &
         (temp_ice(i-1) * CONDUCTIVITY_ICE * thickness_ice(i-1)**(-1) - solar_wave(i-1))
 if (thickness_ice(i).lt.0) then
         thickness_ice(i) = 0
 endif
 enddo

 return

end subroutine

