program main
    use, intrinsic :: iso_fortran_env, only: dp => real128
    implicit none (type, external)

    integer, parameter :: n = 11
    integer :: i
    integer, parameter :: max_iter = 100000
    real(dp), parameter :: err = 1.0e-33_dp

    real(dp), parameter :: PP = 1.0_dp
    real(dp), parameter :: nn = 78.0_dp, no = 21.0_dp

    real(dp) :: x(n)
    real(dp) :: dx(2)
    real(dp) :: f(2)
    real(dp) :: k(n-3) !これは別の場所で代入されている前提とする。

    real(dp) :: max_step, scale_factor
    real(dp) :: max_dx

    !k(:) = [0.11006137d+03, 0.22127333d+01, 0.31448583d-05, 0.21392372d-14, 0.71661612d-15,&
    !      & 0.17097567d-12, 0.11950020d-10, 0.29310428d-07]
    k(:) = [0.11500594q+106, 0.63556694q-81, 0.24410642q-159, 0.39124095q-230, 0.32510580q-245, 0.10257573q-122, 0.37281879q-103, 0.81561706q-52]
    !x(:) = [0.001_dp, 0.001_dp, 0.001_dp, 1.0_dp, 1.0_dp, 0.00001_dp, 0.00001_dp, 0.001_dp, 0.001_dp, 0.0001_dp, 0.0001_dp]

    x(:) = 0.01_dp
    !f = 0.0_dp
    max_step = 2.0_dp

    newrap_loop: do i=1, max_iter

	call two_variables(n, x(1:2), dx, k, f)

	!print *, l
	!print *, f

        if (maxval(abs(f)) < err) then
	    !print *, l
	    call substitution(n, x, k)
	    print *, "--- loop = ", i, " ---"
	    print *, 'Maybe converged?'
            print*, x
            stop
        end if

	!max_dx = maxval(abs(f))

        ! if (max_dx > max_step) then
   !     scale_factor = max_step / max_dx
   ! else
       scale_factor = 1.0_dp
	! end if

        x(1:2) = x(1:2) + scale_factor * dx(1:2)

	!print*, x
    end do newrap_loop

  ! 1.161108219318809396051326940527315E-0041
  ! 1.386817472860182497108565377155488E-0080
  ! 1.851877589960718556744889064738836E-0016
  ! 0.212121212121212028618241714085284
  ! 0.787878787878787873312553524859743
  ! 1.253510140071968045419621872097713E-0185
  ! 1.244098795583205891871379531010384E-0239
  ! 3.815934069178347102227887702156583E-0119
  ! 1.978551000919715250594289050684085E-0177
  ! 3.624008041554558652764309209140806E-0087
  ! 3.624008041554558652764309209140844E-0087
    print *, "Loop reached maximam iteration without convergence."
    
    call substitution(n, x, k)
    print*,x
    !ifx -Warn all -O3 -xHOST -qmkl newrap_two_vars.f90 && ./a.out
    !ifx -O0 -g -traceback -fpe0 -qmkl newrap_two_vars.f90 && ./a.out
contains
    subroutine two_variables(num, p, ddp, KK, ff)
	integer, intent(in) :: num
	real(dp), intent(in) :: p(2)
	real(dp), intent(out) :: ddp(2)
	real(dp), intent(in) :: KK(3:num-1)
	real(dp) :: f11, A, B, f21, f22, C, D
	real(dp), intent(out) :: ff(2)
	f11 = P(1)**2 / KK(4) + P(1) + P(2)**2 / KK(5) + P(2) + P(1) * P(2) * KK(3) - PP
        ff(1) = 4 * (P(1)**2 * KK(8) + P(1) * KK(6) + P(2)**2 * KK(9) + P(2) * KK(7) + P(1) * P(2) * KK(10)) - f11**2
        A = 2 * (2 * (2 * P(1) * KK(8) + KK(6) + P(2) * KK(10)) - (2 * P(1) / KK(4) + 1 + P(2) * KK(3)) * f11)
        B = 2 * (2 * (2 * P(2) * KK(9) + KK(7) + P(1) * KK(10)) - (2 * P(2) / KK(5) + 1 + P(1) * KK(3)) * f11)
        f21 = PP - (P(1)**2 / KK(4) + P(1) + P(2)**2 / KK(5) + P(2) + P(1) * P(2) * KK(3))
        f22 = 2 * nn * P(1)**2 / KK(4) + nn * P(1) - 2 * no * P(2)**2 / KK(5) - no * P(2) + (nn - no) * P(1) * P(2) * KK(3)
        ff(2)=f21*f22+4*nn*P(1)**2*KK(8)+2*nn*P(1)*KK(6)-4*no*P(2)**2*KK(9)-2*no*P(2)*KK(7)+2*(nn-no)*P(1)*P(2)*KK(10)
        C=-(2*P(1)/KK(4)+1+P(2)*KK(3))*f22+f21*(4*nn*P(1)/KK(4)+nn+(nn-no)*P(2)*KK(3))+8*nn*P(1)*KK(8)+2*nn*KK(6)+&
        &2 * (nn - no) * P(2) * KK(10)
        D=-(2*P(2)/KK(5)+1+P(1)*KK(3))*f22-f21*(4*no*P(2)/KK(5)+no-(nn-no)*P(1)*KK(3))-8*no*P(2)*KK(9)-2*no*KK(7)+&
        &2 * (nn - no) * P(1) * KK(10)
        ddp(1) = (-D * ff(1) + B * ff(2)) / (A * D - B * C)
        ddp(2) = (C * ff(1) - A * ff(2)) / (A * D - B * C)
    end subroutine two_variables

    subroutine substitution(num, p, KK)
	integer, intent(in) :: num
	real(dp), intent(inout) :: p(num)
	real(dp), intent(in) :: KK(3:num-1)
        p(3) = p(1) * p(2) * KK(3)
        p(4) = P(1)**2 / KK(4)
        p(5) = p(2)**2 / KK(5)
        p(11)=sqrt(p(1)*KK(6) + p(2)*2*KK(7) + p(1)**2*KK(8) + p(2)**2*KK(9) + p(1)*p(2)*KK(10))
        p(6) = p(1) / p(11) * KK(6)
        p(7) = p(2) / p(11) * KK(7)
        p(8) = p(1)**2 / p(11) * KK(8)
        p(9) = p(2)**2 / p(11) * KK(9)
        p(10) = p(1) * p(2) / p(11) * KK(10)
    end subroutine substitution
end program main
