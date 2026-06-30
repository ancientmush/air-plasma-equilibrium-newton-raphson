module newton_raphson
    use, intrinsic :: iso_fortran_env, only: dp => real64, stderr => error_unit
    use lapack95, only: gesv
    implicit none(type, external)
    private
    public :: main_loop
    external dgesv

    !Backtracking constants and a variable.
    real(dp), parameter :: rate = 0.5_dp
    real(dp) :: alpha
    real(dp), parameter :: c1 = 1.0e-4_dp, c2 = 0.9_dp

contains

    function jacob(num, nn, no, p, KK) result(j)
        integer, intent(in)  :: num
        real(dp), intent(in) :: nn, no
        real(dp), intent(in) :: p(num)
        real(dp), intent(in) :: KK(3:num - 1)
        real(dp) :: j(num, num)
        real(dp) :: ep(num)

        ep = exp(p)

        j(:, 1) = [-2.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, 1.0_dp, -ep(1), 0.0_dp, nn * ep(1)]
        j(:, 2) = [0.0_dp, -2.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 2.0_dp, 1.0_dp, -ep(2), 0.0_dp, -no * ep(2)]
        j(:, 3) = [0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -ep(3), 0.0_dp, (nn - no) * ep(3)]
        j(:, 4) = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -ep(4), 0.0_dp, nn * 2 * ep(4)]
        j(:, 5) = [0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -ep(5), 0.0_dp, -no * 2 * ep(5)]
        j(:, 6) = [0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -ep(6), -ep(6), nn * ep(6)]
        j(:, 7) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -ep(7), -ep(7), -no * ep(7)]
        j(:, 8) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, -ep(8), -ep(8), nn * 2 * ep(8)]
        j(:, 9) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, -ep(9), -ep(9), -no * 2 * ep(9)]
        j(:, 10) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, -ep(10), -ep(10), (nn - no) * ep(10)]
        j(:, 11) = [0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -ep(11), ep(11), 0.0_dp]
    end function jacob

    function func(num, PP, nn, no, p, KK) result(res)
        integer, intent(in) :: num
        real(dp), intent(in) :: PP
        real(dp), intent(in) :: nn, no
        real(dp), intent(in) :: p(num)
        real(dp), intent(in) :: KK(3:num - 1)
        real(dp) :: res(num)
        real(dp) :: ep(num)

        ep = exp(p)

        res(1) = log(KK(4)) - 2 * p(1) + p(4)
        res(2) = log(KK(5)) - 2 * p(2) + p(5)
        res(3) = log(KK(3)) - p(3) + p(2) + p(1)
        res(4) = log(KK(6)) - p(6) - p(11) + p(1)
        res(5) = log(KK(7)) - p(7) - p(11) + p(2)
        res(6) = log(KK(8)) - p(8) - p(11) + 2 * p(1)
        res(7) = log(KK(9)) - p(9) - p(11) + 2 * p(2)
        res(8) = log(KK(10)) - p(10) - p(11) + p(1) + p(2)
        res(9) = PP - sum(ep)
        res(10) = ep(11) - sum(ep(6:10))
        res(11) = nn * (ep(1) + ep(3) + 2 * ep(4) + ep(6) + 2 * ep(8) + ep(10)) &
             & - no * (ep(2) + ep(3) + 2 * ep(5) + ep(7) + 2 * ep(9) + ep(10))
    end function func

    subroutine backtracking(num, PP, nn, no, p_old, KK, delta_p, fu, beta, c)
        integer, intent(in) :: num
        real(dp), intent(in) :: PP
        real(dp), intent(in) :: nn, no
        real(dp), intent(in) :: p_old(num)
        real(dp), intent(in) :: KK(3:num - 1)
        real(dp), intent(in) :: delta_p(num)
        real(dp), intent(in) :: fu(num)
        real(dp), intent(in) :: beta
        real(dp), intent(out) :: c

        real(dp) :: p_new(num), f_new(num)
        real(dp) :: rss_old, rss_new, dir_deriv

        rss_old = 0.5_dp * norm2(fu)**2
        dir_deriv = -2.0_dp * rss_old

        c = 1.0_dp
        do
            p_new = p_old + c * delta_p
            f_new = func(num, PP, nn, no, p_new, KK)
            rss_new = 0.5_dp * norm2(f_new)**2

            if (rss_new <= rss_old + c1 * c * dir_deriv) exit

            c = c * beta
            if (c < 1.0e-10_dp) exit
        end do
    end subroutine backtracking

    subroutine main_loop(n, max_iter, err, P_atm, Nn, No, x, K)
        integer, intent(in) :: n
        integer, intent(in) :: max_iter
        real(dp), intent(in) :: err

        real(dp), intent(in) :: P_atm
        real(dp), intent(in) :: Nn, No

        real(dp), intent(inout) :: x(n)
        real(dp) :: dx(n)
        real(dp) :: jacb(n, n)
        real(dp) :: f(n)
        real(dp), intent(in) :: K(n - 3)

        ! LAPACK用変数
        integer :: ipiv(n)     ! ピボット情報
        integer :: info        ! ステータス

        integer :: i

        newrap_loop: do i = 1, max_iter
            f = func(n, P_atm, Nn, No, x, K)

            if (maxval(abs(f)) < err) then
                ! print *, "Converged at loop = ", i - 1
                return
            end if

            jacb(:, :) = jacob(n, Nn, No, x, K)
            dx = -f

            call dgesv(n, 1, jacb, n, ipiv, dx, n, info)

            if (info /= 0) then
                write (stderr, *) "Error: LAPACK DGESV failed with info = ", info
                return
            end if

            if (maxval(abs(dx)) > 10.0_dp) then
                dx = dx * (10.0_dp / maxval(abs(dx)))
            end if

            call backtracking(n, P_atm, Nn, No, x, K, dx, f, rate, alpha)

            x = x + alpha * dx
        end do newrap_loop

        print *, "Loop reached maximam iteration without convergence."
    end subroutine main_loop
end module newton_raphson
