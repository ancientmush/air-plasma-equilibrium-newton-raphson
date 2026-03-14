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
    subroutine calc_f(num, PP, nn, no, p, ff, KK)
        integer, intent(in) :: num
        real(dp), intent(in) :: PP
        real(dp), intent(in) :: nn, no
        real(dp), intent(in) :: p(num)
        real(dp), intent(out) :: ff(num)
        real(dp), intent(in) :: KK(3:num - 1)

        ff(1) = KK(4) * p(4) - p(1)**2
        ff(2) = KK(5) * p(5) - p(2)**2
        ff(3) = KK(3) * p(1) * p(2) - p(3)
        ff(4) = KK(6) * p(1) - p(6) * p(11)
        ff(5) = KK(7) * p(2) - p(7) * p(11)
        ff(6) = KK(8) * p(1)**2 - p(8) * p(11)
        ff(7) = KK(9) * p(2)**2 - p(9) * p(11)
        ff(8) = KK(10) * p(1) * p(2) - p(10) * p(11)
        ff(9) = PP - sum(p)
        ff(10) = p(11) - sum(p(6:10))
        ff(11) = nn * (p(1) + p(3) + 2 * p(4) + p(6) + 2 * p(8) + p(10)) &
             &- no * (p(2) + p(3) + 2 * p(5) + p(7) + 2 * p(9) + p(10))
    end subroutine calc_f

    subroutine jacub(num, nn, no, p, j, KK)
        integer, intent(in)  :: num
        real(dp), intent(in) :: nn, no
        real(dp), intent(in) ::p(num)
        real(dp), intent(out) :: j(num, num)
        real(dp), intent(in) :: KK(3:num - 1)

        !j(関数i, 変数j)
        j(:, 1) = [-2 * p(1), 0.0_dp, KK(3) * p(2), KK(6), 0.0_dp, 2 * KK(8) * p(1), 0.0_dp, KK(10) * p(2), -1.0_dp, 0.0_dp, nn]
        j(:, 2) = [0.0_dp, -2 * p(2), KK(3) * p(1), 0.0_dp, KK(7), 0.0_dp, 2 * KK(9) * p(2), KK(10) * p(1), -1.0_dp, 0.0_dp, -nn]
        j(:, 3) = [0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, nn - no]
        j(:, 4) = [KK(4), 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 2 * nn]
        j(:, 5) = [0.0_dp, KK(5), 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, -2 * no]
        j(:, 6) = [0.0_dp, 0.0_dp, 0.0_dp, -p(11), 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, -1.0_dp, nn]
        j(:, 7) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -p(11), 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, -1.0_dp, -no]
        j(:, 8) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -p(11), 0.0_dp, 0.0_dp, -1.0_dp, -1.0_dp, 2 * nn]
        j(:, 9) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -p(11), 0.0_dp, -1.0_dp, -1.0_dp, -2 * no]
        j(:, 10) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -p(11), -1.0_dp, -1.0_dp, nn - no]
        j(:, 11) = [0.0_dp, 0.0_dp, 0.0_dp, -p(6), -p(7), -p(8), -p(9), -p(10), -1.0_dp, 1.0_dp, 0.0_dp]
    end subroutine jacub

    function jacob(num, nn, no, p, KK) result(j)
        integer, intent(in)  :: num
        real(dp), intent(in) :: nn, no
        real(dp), intent(in) :: p(num)
        real(dp), intent(in) :: KK(3:num - 1)
        real(dp) :: j(num, num)

        j(:, 1) = [-2.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, 1.0_dp, -exp(p(1)), 0.0_dp, nn * exp(p(1))]
        j(:, 2) = [0.0_dp, -2.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 2.0_dp, 1.0_dp, -exp(p(2)), 0.0_dp, -no * exp(p(2))]
        j(:, 3) = [0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -exp(p(3)), 0.0_dp, (nn - no) * exp(p(3))]
        j(:, 4) = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -exp(p(4)), 0.0_dp, nn * 2 * exp(p(4))]
        j(:, 5) = [0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -exp(p(5)), 0.0_dp, -no * 2 * exp(p(5))]
        j(:, 6) = [0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -exp(p(6)), -exp(p(6)), nn * exp(p(6))]
        j(:, 7) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -exp(p(7)), -exp(p(7)), -no * exp(p(7))]
        j(:, 8) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, -exp(p(8)), -exp(p(8)), nn * 2 * exp(p(8))]
        j(:, 9) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, -exp(p(9)), -exp(p(9)), -no * 2 * exp(p(9))]
      j(:, 10) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, -exp(p(10)), -exp(p(10)), (nn - no) * exp(p(10))]
        j(:, 11) = [0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -exp(p(11)), exp(p(11)), 0.0_dp]
    end function jacob

    function func(num, PP, nn, no, p, KK) result(res)
        integer, intent(in) :: num
        real(dp), intent(in) :: PP
        real(dp), intent(in) :: nn, no
        real(dp), intent(in) :: p(num)
        real(dp), intent(in) :: KK(3:num - 1)
        real(dp) :: res(num)

        res(1) = log(KK(4)) - 2 * p(1) + p(4)
        res(2) = log(KK(5)) - 2 * p(2) + p(5)
        res(3) = log(KK(3)) - p(3) + p(2) + p(1)
        res(4) = log(KK(6)) - p(6) - p(11) + p(1)
        res(5) = log(KK(7)) - p(7) - p(11) + p(2)
        res(6) = log(KK(8)) - p(8) - p(11) + 2 * p(1)
        res(7) = log(KK(5)) - p(9) - p(11) + 2 * p(2)
        res(8) = log(KK(10)) - p(10) - p(11) + p(1) + p(2)
        res(9) = PP - sum(exp(p))
        res(10) = exp(p(11)) - sum(exp(p(6:10)))
        res(11) = nn * (exp(p(1)) + exp(p(3)) + 2 * exp(p(4)) + exp(p(6)) + 2 * exp(p(8)) + exp(p(10))) &
             & - no * (exp(p(2)) + exp(p(3)) + 2 * exp(p(5)) + exp(p(7)) + 2 * exp(p(9)) + exp(p(10)))
    end function func

    function func_rss(num, PP, nn, no, p, KK, fun) result(r)
        integer, intent(in) :: num
        real(dp), intent(in) :: PP
        real(dp), intent(in) :: nn, no
        real(dp), intent(in) :: p(num)
        real(dp), intent(in) :: KK(3:num - 1)
        interface
            function fun(number, d, e1, e2, g, h) result(res)
                import :: dp
                implicit none(type, external)
                integer, intent(in) :: number
                real(dp), intent(in) :: d
                real(dp), intent(in) :: e1, e2
                real(dp), intent(in) :: g(number)
                real(dp), intent(in) :: h(number - 3)
                real(dp) ::res(number)
            end function fun
        end interface
        real(dp) :: r
        real(dp) :: fu(num)

        !fu(:) = func(num, p, KK)
        r = 1 / 2 * (norm2(fun(num, PP, nn, no, p, KK)))**2
    end function func_rss

    subroutine backtracking(num, PP, nn, no, p_old, KK, delta_p, fu, beta, c, fun, fun_r, jac) !rateを返す
        integer, intent(in) :: num
        real(dp), intent(in) :: PP
        real(dp), intent(in) :: nn, no
        real(dp), intent(in) :: p_old(num)
        real(dp), intent(in) :: KK(3:num - 1)
        real(dp), intent(in) :: delta_p(num)
        real(dp), intent(in) :: fu(num)
        real(dp), intent(in) :: beta
        real(dp), intent(inout) :: c
        interface
            function fun(number, d, e1, e2, g, h) result(res)
                import :: dp
                implicit none(type, external)
                integer, intent(in) :: number
                real(dp), intent(in) :: d
                real(dp), intent(in) :: e1, e2
                real(dp), intent(in) :: g(number)
                real(dp), intent(in) :: h(number - 3)
                real(dp) :: res(number)
            end function fun
            function fun_r(numberr, dd, ee1, ee2, gg, hh, fun) result(res1)
                import :: dp
                implicit none(type, external)
                integer, intent(in) :: numberr
                real(dp), intent(in) :: dd
                real(dp), intent(in) :: ee1, ee2
                real(dp), intent(in) :: gg(numberr)
                real(dp), intent(in) :: hh(numberr - 3)
                interface
                    function fun(nummer, u, v1, v2, w, z) result(resl)
                        import :: dp
                        implicit none(type, external)
                        integer, intent(in) :: nummer
                        real(dp), intent(in) :: u
                        real(dp), intent(in) :: v1, v2
                        real(dp), intent(in) :: w(nummer)
                        real(dp), intent(in) :: z(nummer - 3)
                        real(dp) :: resl(nummer)
                    end function fun
                end interface
                real(dp) :: res1
            end function fun_r
            function jac(numel, eee1, eee2, ggg, hhh) result(res2)
                import :: dp
                implicit none(type, external)
                integer, intent(in) :: numel
                real(dp), intent(in) :: eee1, eee2
                real(dp), intent(in) :: ggg(numel)
                real(dp), intent(in) :: hhh(numel - 3)
                real(dp) :: res2(numel, numel)
            end function jac
        end interface
        integer :: l
        real(dp) :: p_new(num), f_new(num), j_new(num, num)
        real(dp) :: rss_old, rss_new

        rss_old = fun_r(num, PP, nn, no, p_old, KK, fun)
        l = 0
        c = beta**0 !Start with 1.
        do while (.not. rss_new <= rss_old - c1 * c * (norm2(fu))**2)
            !print *, c
            p_new = p_old + c * delta_p
            rss_new = fun_r(num, PP, nn, no, p_new, KK, fun)
            f_new = fun(num, PP, nn, no, p_new, KK)
            j_new = jac(num, nn, no, p_new, KK)
            if (dot_product(matmul(f_new, j_new), delta_p) <= c2 * 2.0_dp * rss_old) exit
            l = l + 1
            c = beta**l
        end do
    end subroutine backtracking

    subroutine main_loop(n, max_iter, err, P_atm, Nn, No, x, K)
        integer, intent(in) :: n
        integer, intent(in) :: max_iter
        real(dp), intent(in) :: err

        real(dp), intent(in) :: P_atm
        real(dp), intent(in) :: Nn, No

        real(dp), intent(out) :: x(n)
        real(dp) :: dx(n)
        real(dp) :: jacb(n, n)
        real(dp) :: f(n)
        real(dp), intent(in) :: K(n - 3)

        ! LAPACK用変数
        integer :: ipiv(n)     ! ピボット情報
        integer :: info        ! ステータス

        integer :: i

        newrap_loop: do i = 1, max_iter

            !call calc_f(n, P_atm, Nn, No, x, f, K)
            f = func(n, P_atm, Nn, No, x, K)
            !print *, i
            !print *, f

            if (maxval(abs(f)) < err) then
                !print *, i
                !print *, x
                print *, "Converged at loop = ", i - 1
                print *, exp(x)
                return
            end if

            !call jacub(n, Nn, No, x, jacb, K)
            jacb(:, :) = jacob(n, Nn, No, x, K)

            dx = -f
            !print *, "--- loop = ", i, " ---"
            !print *, "Max F   = ", maxval(abs(f))
            !print *, "Max Jac = ", maxval(abs(jacb))
            call dgesv(n, 1, jacb, n, ipiv, dx, n, info)

            if (info /= 0) then
                write (stderr, *) "Error: LAPACK DGESV failed with info = ", info
                !exit newrap_loop
                return
            end if

            call backtracking(n, P_atm, Nn, No, x, K, dx, f, rate, alpha, func, func_rss, jacob)

            x = x + alpha * dx

            f(:) = 0.0_dp
            jacb(:, :) = 0.0_dp
        end do newrap_loop
        !print *, "loop = ", i - 1
        !print *, sum(exp(x))

        print *, "Loop reached maximam iteration without convergence."
    end subroutine main_loop
end module newton_raphson
