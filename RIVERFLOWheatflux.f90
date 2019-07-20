PROGRAM RIVERFLOW1D

!*    RIVERFLOW1D: river hydrothermodynamic model for INM-MSU                       =
!*    Filename: RIVERFLOW1Dheatflux.f90                                             =
!*    Purpose: modeling of heatflux on ideal river channel                          =
!*    Author: Ivan Korpushenkov, Moscow, Lomonosov MSU, department of land hydrology=
!*    Created:                                                                      =
!*    Copyright: (c) Ivan Korpushenkov 2018                                         =
!*    License:                                                                      =
!*    Date last revised: 23.05.2019                                                 =
!*===================================================================================


!
!*    MAIN PROGRAM                                                                  =
!*    FOR TRAPEZOID CROSS-SECTION AREA OF RIVER                                     =
!*    Trapezoid -- slope angle = 30 degrees!                                        =
!*===================================================================================
use datalist
use data_assimilation
use phys_constants
use hydraulic
use sun_parameters
use surf_scheme1_mod
use thermal
use functions
    !|datalist -
    !|phys_constants -
    !|sun_parameters -
    !|thermal -
    !|functions -
    !|data_assimilation - i/o procedures, boundary and initial conditions
    !|hydraulic - hydraulic radius calculate, manning's roughness coef. calculate
    !|SURF_SCHEME1_MOD - procedure for calculates exchange coefficients in aerodynamic
    !                     formulas (Monin-Obukhov equations)

implicit none


call PBLDAT()
    !|Initializing data for Monin-Obukhov procedure

call FILEREADING(lmax, qmax, tmax, hmax, amax, tribmax)
    !|initializing data reading procedure

ntime = int(qmax - 1)

allocate(sun_h(1:ntime, 2:N-2), astr1(1:ntime, 2:N-2   ),&
&         astr2(1:ntime, 2:N-2), H_sun(1:ntime, 2:N-2),  &
& H_down(1:ntime, 1:N-2), H_hd(1:ntime, 2:N-2),         &
& H_LE(1:ntime, 2:N-2))

allocate(  Ta(1:ntime), water_wap(1:ntime)                   ,&
& nCloud(  1:ntime), T_air(1:ntime, 1:N-2)                   ,&
& wap_water  (1:ntime, 1:N-2), clouds_n(1:ntime, 1:N-2)      ,&
& pressure(1:ntime), w_speed(1:ntime), speed_w(1:ntime,1:N-2),&
& Press(1:ntime, 1:N-2))

allocate(date(1:ntime, 1:3))
allocate(T_int(1:ntime))
allocate(length(2:lmax),                dtime(2:lmax), &
&                              dxx(2:lmax), xii(2:lmax))
allocate(discharge0(1:ntime))
allocate(hriver0(2:tmax), hriver00(1:tmax-1))
allocate(Hriver(2:hmax), Scross(2:hmax))
allocate(No1(2:tribmax),               No2(2:tribmax), &
&                        No3(2:tribmax), No4(2:tribmax))
allocate(trib(1:N, 2:tribmax))

call MODELPARAMETERS(lmax, length, dtime, dxx, xii)

!>River length
Le = length(2)

!>Time step
dt = dtime(2)

!>Lentgh step(meters)
dx = dxx(2)

!>River slope
xi(1:N) = xii(2)

call ACOEF(amax, a1)

!Initial conditions
! Set h, S
call INITIALCONDITIONS(tmax, hriver00)

!>Cross-section initial====================================================

!*Initialize S                                                            =

do j = 1, N
   S(j) = a1 * hriver00(j)
enddo
!>=========================================================================

call BOUNDARYCONDITIONS(qmax, discharge0, T_int)
T_init = Tkel(T_int(1))

!>Temperature initial conditions===========================================

!*Define initial conditions of Tau (delta Tau/delta t)
!*where Tau is T*S

do j = 1, N
   Tau(j) = T_init * S(j)
enddo
!>=========================================================================

!write(0,*) S(1:N), Tau(1:N)

!>Export of results for a .txt file

!do j = 1, N
!   T_s(j) = Tau(j)/S(j)
!enddo

!write(0,*) T_s(1:N)

do i = 1, ntime

  call BOUNDARYCONDITIONS(qmax,               &
  &                          discharge0, T_int)
  call SCROSSAREA(hmax, a1,        medhriver, &
  & medarea,                    Hriver, Scross)
  call HYDRAULICRADIUStime(hmax, Hriver,      &
  &                          Scross, medradius)
  call MANNCOEF(lmax, qmax, xii, discharge0,  &
  & medarea, medradius, xnn)

  xn(1:N) = xnn

  uS(1) = discharge0(i)
  T(1) = Tkel(T_int(i))

  !>thermal conductivity boundary condition
!  T(1) = 11
!  Tau(1) = T(1) * S(1)
!  uTau(1) = Tau(1) * u(1)


  call meteo_data(qmax, Ta, water_wap, nCloud,&
  & pressure, w_speed                         )
  call data_time(                 ntime, date )

  call astronomical_algorithm(pi, date(i,1),  &
  &          date(i,2), date(i,3), lon, astro )
  call define_phys_constants()

  Ta(i) = Tkel(Ta(i))

  astr1(i,2:N-2)  = astro(1)
  astr2(i,2:N-2)  = astro(2)
  sun_h(i, 2:N-2) = sin_h0(pi, lat,   &
  &              astro(1), astro(2)   )

  T_air(i,    2:N-2) =            Ta(i)
  wap_water(i,2:N-2) =     water_wap(i)
  clouds_n(i,2:N-2 ) =        nCloud(i)
  Press(i, 2:N-2   ) =      pressure(i)
  speed_w(i, 2:N-2 ) =  U_2(w_speed(i))

  !>Boundary conditions
  call HCALC(N, a1, S, h)
  call RADIUS(N, S, h, RADIUS_)
!  RADIUS_(:) = sum(RADIUS_(1:N-1))/real(N-1) !do not use for MacCormack scheme

  !>Manning_equation------------------------------------------------------------
  do j = 2, N-1
    u(j) = (1./xn(j)) * (0.5 * (RADIUS_(j-1) + RADIUS_(j))) ** x23 * sqrt(xi(j))
  enddo
  !>----------------------------------------------------------------------------

  u(1) = (1./xn(1))*(RADIUS_(1)) ** x23*sqrt(xi(1))
  u(N) = (1./xn(N))*(RADIUS_(N-1))**x23*sqrt(xi(N))

  !>FOR_DISCHARGE---------------------------------------------------------------
  !>leftmost_point
  uS(2) = u(2) * 0.5 * (S(1) +      S(2))
  xx = -(dt/dx) * (uS(2) - uS(1)) + S(1 )  ! S(1)

  !rightmost point
  uS(N-1)   = u(N-1)  *0.5*(S(N-2) + S (N-1 ))

  !extrapolate to S(N) from S(N-1) and  S(N-2)
  S(N) = S(N-2) + ((S(N-1) - S(N-2))/((N-1) - (N-2)))*(N - (N-2))
  uS(N) = u(N)*0.5*(S(N-1)+S(N))
  xxx = -(dt/dx)*(uS(N) - uS(N-1)) + S(N-1) !S(N-1)
  !>----------------------------------------------------------------------------

  T(1) = 11
  Tau(1) = T(1) * S(1)
  uTau(1) = Tau(1) * u(1)

  !-----------for heat flux-----------------------------------------------------
  !>rightmost point

  uTau(N-2) = u(N-2) * Tau(N-2)
  uTau(N-1) = u(N-1) * Tau(N-1)
  Tau(N-1) = - (dt/dx)*(uTau(N-1) - uTau(N-2)) + Tau(N-1)

  !>extrapolate_to_Tau(N)_from_Tau(N-1)_and_Tau(N-2)
  Tau(N) = Tau(N-2) + ((Tau(N-1) - Tau(N-2))/((N-1) - (N-2)))*(N - (N-2))

  !-----------------------------------------------------------------------------

!  do j = 2, N-2
!
!   uS(j)   = u(j)  *0.5*(S(j-1) +S (j)  )
!   uS(j+1) = u(j+1)*0.5*(S(j)   +S (j+1))
!
!   uTau(j) = u(j)*Tau(j)
!
!  !>Continuity equation(forward scheme)#################################################
!    S(j) = -(dt/dx)*(uS(j+1) - uS(j)) + S(j) ! + ((trib(j, i+1) * dt)/dx)
!  !>####################################################################################
!

  !MacCormack_method====================================================================

  !1. Predictor step
  do j = 2, N-2

      uS(j   ) = u(j) *    S(j)
      uS(j+1 ) = u(j+1)* S(j+1)

      S_pre(j) = S(j) - (dt/dx) * (uS(j+1) - uS(j))

  enddo

  call HCALC(N, a1, S_pre, h)
  call RADIUS(N, S_pre, h, RADIUS_)

  do j = 2, N-1
    u(j) = (1./xn(j)) * (0.5 * (RADIUS_(j-1) + RADIUS_(j))) ** x23 * sqrt(xi(j))
  enddo

  !2. Corrector step
  do j = 2, N-2

      uS(j)   = u(j)*  S_pre(j  )
      uS(j-1) = u(j-1)*S_pre(j-1)

      S(j) = 0.5*((S(j) + S_pre(j)) - 0.5*(dt/dx)*(uS(j) - uS(j-1)))

  enddo
  !=====================================================================================

  !input variables for Monin-Obukhov procedure
  bx2(1) = speed_w(i,j)
  bx2(2) = TetaT(T_air(i,j),    Press(i,j))
  bx2(3) = TetaT((Tau(j)/S(j)), Press(i,j))
  bx2(4) = S_0(Tau(j), S(j), Press(i,j))
  bx2(5) = S_h(wap_water(i,j),  Press(i,j))
  bx2(6) = 2.
  bx2(7) = 1.E-04

  call DRAGVL(bx2, bix2, itdrag)
    !|Initializing Monin-Obukhov procedure for aerodynamic formulas


  uTau(j-1) = u(j-1)*Tau(j-1)
  uTau(j  ) = u(j)  * Tau(j )

  C_1 = dt/(Cp_W*rho_wat)

!  H_sun  = S_wave(pi, S_const, sun_h(i, j), clouds_n(i,j))*(1-albedo_water)  !*C_1
!  H_down = Esum(T_air(i,j), wap_water(i,j),x1,x2,eps,eps_a,sigma)            !*C_1
!  H_hd   = H_D(Cp_A, rho_air, bix2(9), speed_w(i,j), Press(i,j), Tau(j), S(j), T_air(i,j))!*C_1
!  H_LE   = L_E(rho_air, Lwv, bix2(9), speed_w(i,j), S_0(Tau(j), S(j), Press(i,j)),           &
!  & S_h(wap_water(i,j), Press(i,j)))*C_1

  !>Heat_flux_scheme####################################################################
  Tau(j) = - (dt/dx)*(uTau(j) - uTau(j-1)) + Tau(j) ! + H_sun(i,j) + H_long(i,j) !+ H_LE(i,j)   ! + &
  & S_wave(pi, S_const, sun_h(i, j), clouds_n(i,j))*(1-albedo_water)*C_1                    + &
  & Esum(Tau(j), S(j), T_air(i,j), wap_water(i,j),x1,x2,eps,eps_a,sigma)*C_1                - &
  & H_D(Cp_A, rho_air, bix2(9), speed_w(i,j), Press(i,j), Tau(j), S(j), T_air(i,j))*C_1     - &
  & L_E(rho_air, Lwv, bix2(9), speed_w(i,j), S_0(Tau(j), S(j), Press(i,j)),                   &
  & S_h(wap_water(i,j), Press(i,j)))*C_1
  !>#############################################################################################

  enddo
  S(1  ) = xx
  S(N-1) = xxx



enddo

END PROGRAM RIVERFLOW1D

