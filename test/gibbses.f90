module temperature
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    private
    public :: T
    integer :: i
    real(kind=real64), parameter :: T(198) = [298.15d0, (300.0d0 + (i - 1) * 100.0d0, i=1, 197)]
end module temperature

module gfe
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    private
    public :: Species, gibbs, polynomy, enthalpy, entropy!, state

    TYPE polynomy
        real(real64) :: tmin, tmax         !minimum & maximum temperature in the segement.
        real(real64), dimension(7) :: a
        real(real64), dimension(2) :: b
    end type polynomy

    TYPE Species
        character(LEN=:), allocatable :: name
        integer :: num          !number of segments
        type(polynomy), allocatable :: c(:)
        real(real64), allocatable :: G(:)
    end type Species

contains
    real(real64) elemental function enthalpy(sp, T) !H/(RT)
        type(Species), INTENT(in) :: sp
        real(real64), INTENT(IN) :: T
        integer :: i
        real(real64) :: a(7), b(2)
        do i = 1, sp%num
            if (T >= sp%c(i)%tmin .and. T <= sp%c(i)%tmax) then
                a = sp%c(i)%a
                b = sp%c(i)%b
                enthalpy = (-a(1) / T + a(2) * log(T) + b(1)) / T + a(3) + &
                     &T * (a(4) / 2 + T * (a(5) / 3 + T * (a(6) / 4 + T * (a(7) / 5))))
            else
                CYCLE
            end if
        end do
    end function enthalpy

    real(real64) elemental function entropy(sp, T) !S/R
        type(Species), intent(in) :: sp
        real(real64), INTENT(in) :: T
        integer :: i
        real(real64) :: a(7), b(2)
        do i = 1, sp%num
            if (T >= sp%c(i)%tmin .and. T <= sp%c(i)%tmax) then
                a = sp%c(i)%a
                b = sp%c(i)%b
                entropy = ((-a(1) / 2) / T - a(2)) / T + a(3) * log(T) + b(2) + &
                     &T * (a(4) + T * (a(5) / 2 + T * (a(6) / 3 + T * (a(7) / 4))))
            else
                CYCLE
            end if
        end do
    end function entropy

    real(real64) elemental function gibbs(sp, T)
        type(Species), intent(in) :: sp
        real(real64), intent(in) :: T
        gibbs = T * (enthalpy(sp, T) - entropy(sp, T))
    end function gibbs

end module gfe

program main
    use, intrinsic :: iso_fortran_env, only: real64
    use temperature, only: T
    use gfe, only: Species, gibbs, polynomy, enthalpy, entropy!, state

    implicit none

    integer :: i
    real(real64), parameter :: R = 8.314462618
    real(real64), allocatable :: H(:)

    type(species) :: n
    type(Species) :: n2

    n%name = 'N'
    n%num = 3
    allocate (n%c(n%num))
    n%c(1) = polynomy(200.0_real64, 1000.0_real64,&
         & [real(real64) :: 0, 0, 2.5, 0, 0, 0, 0], [5.610463780D+04, 4.193905026D+00])
    n%c(2) = polynomy(1000.0_real64, 6000.0_real64,&
         & [8.876501380D+04, -1.071231500D+02, 2.362188287D+00, 2.916720081D-04,&
         &-1.729515100D-07, 4.012657880D-11, -2.677227571D-15], [5.697351330D+04, 4.865231506D+00])
    n%c(3) = polynomy(6000.0_real64, 20000.0_real64,&
         & [5.475181050D+08, -3.107574980D+05, 6.916782740D+01, -6.847988130D-03,&
         &3.827572400D-07, -1.098367709D-11, 1.277986024D-16], [2.550585618D+06, -5.848769753D+02])

    n2%name = 'N2'
    n2%num = 3
    allocate (n2%c(n2%num))
    n2%c(1) = polynomy(200.0_real64, 1000.0_real64,&
         & [2.210371497D+04, -3.818461820D+02, 6.082738360D+00, -8.530914410D-03,&
         &1.384646189D-05, -9.625793620D-09, 2.519705809D-12], [7.108460860D+02, -1.076003744D+01])
    n2%c(2) = polynomy(1000.0_real64, 6000.0_real64,&
         & [5.877124060D+05, -2.239249073D+03, 6.066949220D+00, -6.139685500D-04,&
         &1.491806679D-07, -1.923105485D-11, 1.061954386D-15], [1.283210415D+04, -1.586640027D+01])
    n2%c(3) = polynomy(6000.0_real64, 20000.0_real64,&
         & [8.310139160D+08, -6.420733540D+05, 2.020264635D+02, -3.065092046D-02,&
         &2.486903333D-06, -9.705954110D-11, 1.437538881D-15], [4.938707040D+06, -1.672099740D+03])

    allocate (n%G(size(T)))
    allocate (n2%G(size(T)))
    allocate (H(size(T)))

    H = R * T * (enthalpy(n, T) - enthalpy(n2, T) / 2) ! 窒素原子の生成エンタルピー
    !n%G = H - R * T * (entropy(n, T) - entropy(n2, T) / 2) ! 窒素原子の標準生成ギブスエネルギー
    n%G = 2 * R * gibbs(n, T) - R * gibbs(n2, T)
    do i = 1, size(T)
        print *, T(i), H(i) / 1000_real64, n%G(i) / 1000_real64
    end do
end program main

