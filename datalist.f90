module datalist

 implicit none

 !>Number of mesh cells
 integer, parameter :: N = 100 
 integer :: ntime 

 real, parameter :: x23 = 2./3.
 real, parameter :: dampcoef = 5.E-02 !FIltering coef. var1=5.E-02 var2=44.E-03 

 !>Loop indices
 integer :: i, j

 !>The numerical scheme dimensions
 real :: S(1:N), h(1:N-1), u(1:N), uS(1:N),       xi(1:N),      xn(1:N),             & 
 &       Q(1:N-1), RADIUS_(1:N), work(1:N), T(1:N), Tau(1:N),    uTau(1:N), E_e(1:N),&
 & T_s(1:N), S_pre(1:N)

 !atmospheris parameters
 real :: pressure_a
 real, parameter :: Sf = 0
! real :: sin_h0
 
 !sun parameters
 real :: HRA
       !|hour angle of the sun
 real :: delta
       !|declination of the sun
 real, parameter :: lat = 44
       !|latitude
 real, parameter :: lon = -94
       !|longitude (if minus - its west hemisphere) 
 real :: S_constant
       !|insolation
 real :: astro(2)
       !|incoming radiation fluxes

 !>PROGONKA coefs
 real :: a(1:N), b(1:N), c(1:N), f(1:N)

 !>River length
 real :: Le
 
 !>Steps
 integer :: dt
 real ::  dx
 real :: xx, xxx
 real :: xnn

 real :: medarea, medhriver, medradius
 real :: a1, a2

 integer :: lmax
 integer :: qmax
 integer :: tmax
 integer :: hmax
 integer :: amax
 integer :: tribmax
 integer :: hmaxxx

 !----COS(X) parameters-------------------------------------------
! integer, parameter :: pi = 3.14                                 
 real, parameter :: amp = 30                                     
 real, parameter :: Am = 11
 real :: omega                                                  
 real :: teta
 !----------------------------------------------------------------

 !SUBROUTINE MODELPARAMETERS dimensions
 integer, allocatable :: dtime(:)
 real, allocatable :: length(:),  dxx(:), xii(:)

 !SUBROUTINE BOUNDARYCONDITIONS dimensions
 real, allocatable :: discharge0(:)

 !SUBROUTINE TRIBUTARIES dimensions
 real, allocatable :: trib(:,:)
 real, allocatable :: No1(:)
 real, allocatable :: No2(:)
 real, allocatable :: No3(:)
 real, allocatable :: No4(:)

 !SUBROUTINE INITIALCONDITIONS dimensions
 real, allocatable :: hriver0(:), hriver00(:)

 !SUBROUTINE SCROSSAREA dimensions
 real, allocatable :: Hriver(:), Scross(:)

 !subroutine meteo_data
 real, allocatable :: Ta(:), water_wap(:), nCloud(:), &
 & T_air(:,:), wap_water(:,:), clouds_n(:,:),         &
 & pressure(:), w_speed(:), speed_w(:,:),    Press(:,:)

 !astronomical and time data
 real, allocatable :: date(:,:)
 real, allocatable :: sun_h(:,:)
 real, allocatable :: astr1(:,:),astr2(:,:)
 real, allocatable :: H_sun(:,:), H_down(:,:), H_hd(:,:), &
 & H_LE(:,:)
 
 !input temp (temp measurement by USGS)
 real, allocatable :: T_int(:)
 real :: T_init

 !dt/Cp_W*1000 constant value
 real :: C_1 

 !Some parameters of Monin-Obukhov procedure
 real(kind=8) :: bx2(7)
 real(kind=8) :: bix2(11)
 integer(kind=4), parameter :: itdrag = 10
    !|number of iterations in Monin-Obukhov procedure
      
end module datalist
