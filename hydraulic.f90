module hydraulic

 contains

 subroutine HYDRAULICRADIUStime(hmax, Hriver, Scross, medradius)

 implicit none
 !integer :: IOS

 integer :: hmax 

 !Input variables
 real :: Hriver(2:hmax)
 real :: d(2:hmax)
 real :: A(2:hmax)
 real :: Scross(2:hmax)
 real :: S(2:hmax)
 real :: B(2:hmax)
 real :: m(2:hmax)
 real :: Pwater(2:hmax)
 real :: radius_(2:hmax)
 real :: medradius

 !Loop indices
 integer :: i, k

 do k = 2, hmax 
    d(k) = Hriver( k) / real ( sqrt(3.  ))
    A(k) = (Hriver(k) ** 2)/real(sqrt(3.))
    S(k) = Scross(k) - A(k)
    B(k) = S(k)/Hriver(k)
    Pwater(k) = B(k) + 2 * sqrt(Hriver(k) ** 2  + d(k) ** 2)
    radius_(k) = Scross(k)/Pwater(k)
 enddo

 !Calculate of hriver average value
 medradius = 0
 do i = 2, hmax
    medradius = medradius + radius_(i)
 enddo
 medradius = medradius/hmax

 end subroutine HYDRAULICRADIUStime

 subroutine MANNCOEF(lmax, qmax, xii, discharge0, medarea, medradius, xnn)

 implicit none
 
 !Input variables
 integer :: lmax
 integer :: qmax

 !real :: disch(1:tmax)
 real :: discharge0(1:qmax)
 real, parameter :: x23 = 2./3 
 real :: medarea
 real :: xii(2:lmax)
 real :: medradius

 !Local variables
 real :: meddisch
 real :: xnn

 integer :: i

 !Calculate of average value
 meddisch = 0
 do i = 2, qmax
    meddisch = meddisch + discharge0(i)
 enddo
 meddisch = meddisch/qmax


 xnn = (meddisch/medarea) * sqrt(xii(2)) * (medradius ** x23) 

 end subroutine MANNCOEF

SUBROUTINE HCALC(N, a1, S, h)

 implicit none

 integer ::  N 

 real :: a1
 real :: S(1:N-1)
 real :: h(1:N-1)
 integer :: j

 do j = 1, N-1
    h(j) = (S(j)) / (a1)
 enddo

END SUBROUTINE HCALC

subroutine RADIUS(N, S, h, RADIUS_)

!--------------------------------------------------------------------------------------
!
!SUBROUTINE for a calculate of hydraulic radius
!
!Input data: N(number of mesh cells), S(cross-section area of river), h(depth of river) 
!
!Output data: RADIUS_(hydraulic radius of river)
!
!--------------------------------------------------------------------------------------

 implicit none

 integer N 

 !Input variables
 real :: S(1:N-1), h(1:N-1), d(1:N-1), A(1:N-1), Sc(1:N-1), B(1:N-1), m(1:N-1),      &
 & Pw(1:N-1), RAD(1:N-1), RADIUS_(1:N-1)

! real :: h(1:N-1)
! real :: d(1:N-1)
! real :: A(1:N-1)
! real :: Sc(1:N-1)
! real :: B(1:N-1)
! real :: m(1:N-1)
! real :: Pw(1:N-1)
! real :: RAD(1:N-1)
! real :: RADIUS_(1:N-1)

 integer :: j

 do j = 1, N-1 
  
    d(j  )  = h(j)        /            sqrt(3.)
    A(j  )  = (h(j) ** 2) /            sqrt(3.)
    Sc(j )  = S(j)        -                A(j)
    B(j  )  = Sc(j)       /                h(j)
    Pw(j )  = B(j) + 2.*sqrt(h(j)**2 + d(j)**2)
    RAD(j)  = S(j)        /               Pw(j)

 enddo 

RADIUS_(1:N-1) = RAD(1:N-1)

!implicit none

!integer :: N
!  integer,intent(in) :: N 
!  real   ,intent(in) :: S(1:N-1), h(1:N-1)

!  real :: d(1:N-1), A(1:N-1), Sc(1:N-1), B(1:N-1), &
!  & m(1:N-1), Pw(1:N-1),                  RAD(1:N-1)
!  integer :: j

!  real, intent(out) :: RADIUS_(1:N-1)

!real :: S(1:N-1)
!real :: h(1:N-1)
!real :: d(1:N-1)
!real :: A(1:N-1)
!real :: Sc(1:N-1)
!real :: B(1:N-1)
!real :: m(1:N-1)
!real :: Pw(1:N-1)
!real :: RAD(1:N-1)
!real :: RADIUS_(1:N-1)

!  integer :: j
!
!  do j = 1, N-1 
!     d(j) = h(j) / sqrt(3.)
!     A(j) = (h(j) ** 2) / sqrt(3.)
!     Sc(j) = S(j) - A(j)
!     B(j) = Sc(j) / h(j)
!     Pw(j) = B(j) + 2.*sqrt(h(j)**2 + d(j)**2)
!     RAD(j) = S(j) / Pw(j)
!  enddo 

!open(12,file = '../Output/tempTest.txt' , status = 'unknown')
!do i = 1, ntime
!  RADIUS_(1:N-1) = RAD(1:N-1)
!  write(12,*) S(1)
!enddo
!close(12)

end subroutine RADIUS

end module
