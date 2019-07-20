module data_assimilation

 contains


 subroutine FILEREADING(lmax, qmax, tmax, hmax, amax, tribmax)

  implicit none
  integer :: IOS

  integer :: Lmaxx
  integer :: lmax

  integer :: Qmaxx
  integer, intent(out) :: qmax

  integer :: Tmaxx
  integer :: tmax

  integer :: Hmaxx
  integer :: hmax

  integer :: Amaxx
  integer :: amax

  integer :: Tribmaxx
  integer :: tribmax

  integer :: metmaxx
  integer :: metmax

  !Modelparameters.txt   reading---------------------------------------------------------------------------------
  open(1, file = 'D:\Магистерская\Модель\RIVERFLOW1D\Modelinfo\Modelparameters.txt', status = 'old')
  read(1, *)
  Lmaxx = 0
  do
     Lmaxx = Lmaxx + 1
     read(1,*, iostat = IOS)
  if (IOS.lt.0) exit
  enddo
  close(1)
  lmax = Lmaxx

  !------------------------------------------------------------------------------------------------------------

  !Boundary.txt reading----------------------------------------------------------------------------------------
  open(2, file = 'D:\Магистерская\Модель\RIVERFLOW1D\Modelinfo\boundary2018.txt', status = 'old')
  read(2, *)
  Qmaxx = 0
  do
     Qmaxx = Qmaxx + 1
     read(2,*,iostat = IOS)
  if (IOS.lt.0) exit
  enddo
  close(2)
  qmax = Qmaxx
  !------------------------------------------------------------------------------------------------------------

  !Initial.txt reading-----------------------------------------------------------------------------------------
  open(3, file = 'D:\Магистерская\Модель\RIVERFLOW1D\Modelinfo\initial2018.txt', status = 'old')
  read(3, *)
  Tmaxx = 0
  do
     Tmaxx = Tmaxx + 1
     read(3,*,iostat = IOS)
  if (IOS.lt.0) exit
  enddo
  close(3)
  tmax = Tmaxx
  !------------------------------------------------------------------------------------------------------------

  !hriver.txt reading------------------------------------------------------------------------------------------
  open(4, file = 'D:\Магистерская\Модель\RIVERFLOW1D\Modelinfo\hriver2018.txt', status = 'old')
  read(4, *)
  Hmaxx = 0
  do
     Hmaxx = Hmaxx + 1
     read(4,*, iostat = IOS)
  if (IOS.lt.0) exit
  enddo
  close(4)
  hmax = Hmaxx
  !------------------------------------------------------------------------------------------------------------

  !BandH.txt reading-------------------------------------------------------------------------------------------
  open(5, file = 'D:\Магистерская\Модель\RIVERFLOW1D\Modelinfo\BandHMinnesota.txt', status = 'old')
  read(5, *)
  Amaxx = 0
  do
     Amaxx = Amaxx + 1
     read(5, *, iostat = IOS)
  if (IOS.lt.0) exit
  enddo
  close(5)
  amax = Amaxx
  !------------------------------------------------------------------------------------------------------------

  !tributaries.txt reading-------------------------------------------------------------------------------------
  !open(6, file = '../Modelinfo/tributaries1982.txt', status = 'old')
  !read(6, *)
  !Tribmaxx = 0
  !do
  !   Tribmaxx = Tribmaxx + 1
  !   read(6, *, iostat = IOS)
  !if (IOS.lt.0) exit
  !enddo
  !close(6)
  !tribmax = Tribmaxx


  !------------------------------------------------------------------------------------------------------------


  !>air temperature reading------------------------------------------------------------------------------------
  !open(7, file = '../Modelinfo/', status = 'old')
  !read(7,*)
  !metmax = 0
  !do
  !  metmax = metmax + 1
  !   read(7, *, iostat = IOS)
  !if (IOS.lt.0) exit
  !enddo
  !close(7)
  !metmax = metmax

 end subroutine


subroutine MODELPARAMETERS(lmax, length, dtime, dxx, xii)

implicit none
integer :: IOS


integer :: lmax
real :: length(2:lmax)
integer :: dtime(2:lmax)
real :: dxx(2:lmax)
real :: xii(2:lmax)
real :: Length1(2:lmax)
integer :: Dtime1(2:lmax)
real :: Dxx1(2:lmax)
real :: Xii1(2:lmax)

!Loop indices
integer :: i


open(1, file = 'D:\Магистерская\Модель\RIVERFLOW1D\Modelinfo\Modelparameters.txt', status = 'old')
read(1, *)
do i = 2, lmax
   read(1,*) Length1(i), Dtime1(i), Dxx1(i), Xii1(i)
if (IOS.lt.0) exit
enddo
close(1)

length(2:lmax) = Length1(2:lmax)
dtime(2:lmax) = Dtime1(2:lmax)
dxx(2:lmax) = Dxx1(2:lmax)
xii(2:lmax) = Xii1(2:lmax)




end subroutine MODELPARAMETERS


subroutine BOUNDARYCONDITIONS(qmax, discharge0, T_int)

implicit none
integer :: IOS

!Input variables
integer, intent(in) :: qmax

!Local variables
real :: time(2:qmax)
real, intent(out) :: discharge0(2:qmax)

!real :: Time11(2:qmax)
real :: Discharge11(2:qmax)
character :: workchar*8

real :: T_xxx(1:qmax)

real, intent(out) :: T_int(1:qmax)

!Loop indices
integer :: i

open(2, file = 'D:\Магистерская\Модель\RIVERFLOW1D\Modelinfo\boundary2018.txt', status = 'old')
read(2, *, iostat = IOS)
do i = 2, qmax
    read(2,*, iostat = IOS) Discharge11(i)
if(IOS.lt.0) exit
enddo
close(2)
discharge0(2:qmax) = Discharge11(2:qmax)


 open(11, file = 'D:\Магистерская\Модель\RIVERFLOW1D\Modelinfo\tint.txt', status='old')
 read(11,*)
 do i = 1, qmax
     read(11,*, iostat=IOS) T_xxx(i)
 if (IOS.lt.0) exit
 enddo
 close(11)
 T_int(1:qmax) = T_xxx(1:qmax)


end subroutine BOUNDARYCONDITIONS


!subroutine TRIBUTARIES(tribmax, No1, No2, No3, No4)
!
!implicit none
!integer :: IOS
!
!Input variables
!integer :: tribmax
!
!Local variables
!real :: No1(2:tribmax)
!real :: No2(2:tribmax)
!real :: No3(2:tribmax)
!real :: No4(2:tribmax)
!
!Loop indices
!integer :: i
!
!open(6, file = '../Modelinfo/tributaries1982.txt', status = 'old')
!read(6, *, iostat = IOS)
!do i = 6, tribmax
!    read(6,*, iostat = IOS) No1(i), No2(i), No3(i), No4(i)
!if(IOS.lt.0) exit
!enddo
!close(6)
!
!No1(2:tribmax) = No1(2:tribmax)
!No2(2:tribmax) = No2(2:tribmax)
!No3(2:tribmax) = No3(2:tribmax)
!No4(2:tribmax) = No4(2:tribmax)
!
!>open(17, file = '../Output/test.txt' ,status = 'unknown')
!>close(17)
!
!
!
!end subroutine TRIBUTARIES

subroutine INITIALCONDITIONS(tmax, hriver00)


implicit none
integer ::  IOS

integer :: tmax
real :: lx(2:tmax)
real :: hriver0(2:tmax)
real :: hriver00(1:tmax-1)
real :: Lxx(2:tmax)
real :: Hriver011(2:tmax)

integer :: i

open(3, file = 'D:\Магистерская\Модель\RIVERFLOW1D\Modelinfo\initial2018.txt', status = 'old')
read(3,*, iostat = IOS)
do i = 2, tmax
   read(3,*,iostat = IOS) Lxx(i), Hriver011(i)
if(IOS.lt.0) exit
enddo
close(3)

lx(2:tmax) = Lxx(2:tmax)
hriver00(1:tmax-1) = Hriver011(2:tmax)



end subroutine INITIALCONDITIONS

subroutine ACOEF(amax, a1)

!Isosceles trapezoid river channel cross-section area

implicit none
integer :: IOS
integer :: amax

!Input variables
real :: B(2:     amax)
real :: Briv(2:  amax) !river water width
real :: crossar(1:  2) !cross-section area
real :: H0(2:    amax)
real :: hriver(2:amax) !river depth

real :: dS1, dS2, dh1, dh2

!linear function coefficient's
real :: a11
real :: a1

!Loop incidens
integer ::  i


open(5, file = 'D:\Магистерская\Модель\RIVERFLOW1D\Modelinfo\BandHMinnesota.txt', status = 'old')
read(5, *, iostat = IOS)
do i = 2, amax
   read(5, *, iostat = IOS) H0(i), B(i)
if(IOS.lt.0) exit
enddo
close(5)

hriver(2:amax) = H0(2:amax)
Briv(2:amax) = B(2:amax)

!a1 coef calculation
a1 = Briv(2) -(hriver(2) * sqrt(3.))


end subroutine ACOEF

subroutine SCROSSAREA(hmax, a1, medhriver, medarea, Hriver, Scross)

implicit none


integer :: IOS
integer :: hmax



!Local variables
real :: medhriver, medarea
real :: medhriver11
real :: Hriver(2:hmax)
real :: Scross(2:hmax)
real :: hriver011(2:hmax)
real :: scrossxx(2:hmax)
real :: a1, a2


!Loop indices
integer :: i

open(4, file = 'D:\Магистерская\Модель\RIVERFLOW1D\Modelinfo\hriver2018.txt', status = 'old')
read(4,*)
do i = 2, hmax
   read(4,*,iostat = IOS) hriver011(i)
if(IOS.lt.0) exit
enddo
close(4)

!Calculate of hriver average value
medhriver = 0
do i = 2, hmax
   medhriver11 = medhriver11 + hriver011(i)
enddo
medhriver = medhriver11/hmax
medarea = a1 * medhriver


do i = 2, hmax
   scrossxx(i) = a1 * hriver011(i)
enddo

Hriver(2:hmax) = hriver011(2:hmax)
Scross(2:hmax) = scrossxx(2:hmax)



end subroutine SCROSSAREA


subroutine meteo_data(qmax, air_T, w_W, n_C, pressure, w_speed)

 implicit none


 integer,intent(in) :: qmax

 integer :: k, IOS
 real :: air_xxxx(1:qmax), wW_xxxx(1:qmax), nC_xxxx(1:qmax),   &
 &                                p_xxxx(1:qmax), w_xxxx(1:qmax)

 real, intent(out) :: air_T(1:qmax), w_W(1:qmax), n_C(1:qmax), &
 &                             pressure(1:qmax), w_speed(1:qmax)


 open(7, file = 'D:\Магистерская\Модель\RIVERFLOW1D\Modelinfo\air.txt', status='old')
 read(7,*)
 do k = 1, qmax
     read(7,*, iostat=IOS) air_xxxx(k)
 if (IOS.lt.0) exit
 enddo
 close(7)
 air_T(1:qmax) = air_xxxx(1:qmax)

 open(8, file = 'D:\Магистерская\Модель\RIVERFLOW1D\Modelinfo\waterw.txt', status='old')
 read(8,*)
 do k = 1, qmax
     read(8,*, iostat=IOS) wW_xxxx(k)
 if (IOS.lt.0) exit
 enddo
 close(8)
 w_W(1:qmax) = wW_xxxx(1:qmax)

 open(9, file = 'D:\Магистерская\Модель\RIVERFLOW1D\Modelinfo\clouds.txt', status='old')
 read(9,*)
 do k = 1, qmax
     read(9,*, iostat=IOS) nC_xxxx(k)
 if (IOS.lt.0) exit
 enddo
 close(9)
 n_C(1:qmax) = nC_xxxx(1:qmax)

 open(12, file = 'D:\Магистерская\Модель\RIVERFLOW1D\Modelinfo\pressurem.txt', status='old')
 read(12,*)
 do k = 1, qmax
     read(12,*, iostat=IOS) p_xxxx(k)
 if (IOS.lt.0) exit
 enddo
 close(12)
 pressure(1:qmax) = p_xxxx(1:qmax)

 open(13, file = 'D:\Магистерская\Модель\RIVERFLOW1D\Modelinfo\Wspeed.txt', status='old')
 read(13,*)
 do k = 1, qmax
     read(13,*, iostat=IOS) w_xxxx(k)
 if (IOS.lt.0) exit
 enddo
 close(13)
 w_speed(1:qmax) = w_xxxx(1:qmax)

end subroutine meteo_data

subroutine data_time(ntime, date)

 implicit none

 integer, intent(in) :: ntime
 integer :: n = 3
 integer :: k_1, k_2, IOS
 real :: dat(1:ntime, 1:3)
 real, intent(out) :: date(1:ntime, 1:3)

 open(10, file = 'D:\Магистерская\Модель\RIVERFLOW1D\Modelinfo\dates.txt', status='old')
 do k_1 = 1, ntime
     read(10,*, iostat=IOS) (dat(k_1, k_2), k_2 = 1,3)
 enddo
 close(10)
 date(1:ntime,1:3) = dat(1:ntime,1:3)

end subroutine data_time

!SUBROUTINE OUTPUTFILE()
!
!character(len=10) :: date_char
!character(len=2) :: dd, mm, yy
!integer :: day, month, year
!
!
!write(dd,'(a2)') day
!write(mm,'(a2)') month
!write(yy,'(a2)') year
!date_char = dd//'/'//mm//'/'//yy
!
!write(1,*) date_char, disccharge(1)
!
!END SUBROUTINE OUTPUTFILE

subroutine PROGONKA(a, b, c, f, y, K, N)
implicit none
!FACTORIZATION METHOD FOR THE FOLLOWING SYSTEM OF LINEAR EQUATIONS:
!-a(i)*y(i-1)+c(i)*y(i)-b(i)*y(i+1)=f(i) i=K+1,N-1
! c(K)*y(K)-b(K)*y(K+1)=f(K)
!-a(N)*y(N-1)+c(N)*y(N)=f(N)
!
!
integer, intent(in) :: K, N
real, intent(in) :: a(K:N), b(K:N), c(K:N), f(K:N)
real, intent(out) :: y(K:N)

real :: alpha(K+1:N+1), beta(K+1:N+1)
integer :: i

SAVE

alpha(K+1) = b(K)/c(K)
beta(K+1) = f(K)/c(K)
do i = K+2, N+1
   alpha(i) = b(i-1)/(c(i-1)-a(i-1)*alpha(i-1))
   beta(i) = (f(i-1)+a(i-1)*beta(i-1))/ &
   & (c(i-1)-a(i-1)*alpha(i-1))
end do
y(N) = beta(N+1)
do i = N-1, K, -1
   y(i) = alpha(i+1)*y(i+1)+beta(i+1)
end do

end subroutine PROGONKA





end module data_assimilation
