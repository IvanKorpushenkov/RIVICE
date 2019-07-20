module functions

contains

  real function TKel(tCel)
  !=       !|from Cel. degree to Kelvin degree          =
  !=       Input: tCel -- temperature at Celsius degrees=
  !=                                                    =
  !======================================================

  real, intent(in ) :: tCel

  TKel = 273.15 + tCel

 end function TKel

 real function Sp_humidity(water_wap, pressure_a) 
  != calculating of specific humidity                 =
  != Input: water_wap -- water wapor pressure         = 
  !=        pressure_a -- current atmospheric pressure=
  !=                                                  =
  !====================================================

  real, intent(in) :: water_wap
  real, intent(in) :: pressure_a

  Sp_humidity = (0.623*water_wap) / &
  & (pressure_a - 0.377*water_wap)


 end function Sp_humidity



end module functions
