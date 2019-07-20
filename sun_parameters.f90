module sun_parameters

 
 contains

 subroutine astronomical_algorithm(pi, LC, day_i, month, lon, astro)
  
  implicit none

  !in
  real, intent(in) :: pi
  real, intent(in) :: LC,day_i, month
  real, intent(in) :: lon

  real :: N_zone, B_d, EoT, LSTM, TC    
  real :: tau_o, tau_m, HRA, delta  
  real :: day

  !out
  real :: ast(2)
  real, intent(out) :: astro(2)  
  
  
  if     (month == 1)  then
      day = day_i
  elseif (month == 2)  then
      day = day_i + 31
  elseif (month == 3)  then
      day = day_i + 59
  elseif (month == 4)  then
      day = day_i + 90
  elseif (month == 5)  then
      day = day_i + 120
  elseif (month == 6)  then
      day = day_i + 151
  elseif (month == 7)  then
      day = day_i + 181
  elseif (month == 8)  then
      day = day_i + 212
  elseif (month == 9)  then
      day = day_i + 243
  elseif (month == 10) then
      day = day_i + 273
  elseif (month == 11) then
      day = day_i + 304
  else
      day = day_i + 334
  endif

  !correction
!  B_d = (real(360) / real(365))*(day - 81)  
!  B_d = (360*(day - 81))/ 365
   B_d = 2*pi*(day-81)/365
 
  !Equation of time(minutes!)
!  EoT = 7.53*cos(B_d) + 1.5*sin(B_d) - 9.87*sin(2*B_d)
  EoT = 9.87*sin(2*B_d) - 7.67*sin(B_d + 78.7)

  !define number of zone
  N_zone = int(abs(lon) / 15)

  !define tau_m (average time of zone)
  if (lon < 0) then
      tau_m = 0 + N_zone
  else
      tau_m = 12 + N_zone
  end if

  LSTM = 15*N_zone 
  TC = (4*(LSTM - lon))/real(60) + EoT/real(60)   
 
  !local solar time(LST)
!  tau_o = LC + (4*(lon - 15*N_zone)-EoT) / 60
!  tau_o = LC + lon + EoT
  tau_o =  LC + TC/real(60)

  !Hour angle(HRA) 
  HRA = 15*(tau_o - 12)

  !angle of declination of the sun
  delta = 23.45*sin(B_d)   

  ast(1) = delta !degree
  ast(2) = HRA !degree

  astro = ast

  return
 end subroutine astronomical_algorithm


end module sun_parameters
