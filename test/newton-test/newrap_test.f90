program main
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use lapack95, only: gesv
    implicit none(type, external)

    external :: dgesv

    integer, parameter :: n = 11
    integer :: i
    integer, parameter :: max_iter = 5000
    real(dp), parameter :: err = 1.0e-15_dp

    real(dp), parameter :: PP = 1.0_dp
    real(dp), parameter :: nn = 78.0_dp, no = 21.0_dp

    real(dp) :: x(n)
    real(dp) :: dx(n)
    real(dp) :: jacb(n, n)
    real(dp) :: f(n)
    real(dp) :: k(n - 3) !これは別の場所で代入されている前提とする。

    ! LAPACK用変数
    integer :: ipiv(n)     ! ピボット情報
    integer :: info        ! ステータス

    real(dp), parameter :: rate = 0.5_dp
    real(dp) :: alpha
    real(dp), parameter :: c1 = 1.0e-4_dp, c2 = 0.9_dp

    !k(:) = [0.11006137d+03, 0.22127333d+01, 0.31448583d-05, 0.21392372d-14, 0.71661612d-15, 0.17097567d-12, 0.11950020d-10, 0.29310428d-07]
    k(:) = [0.11500594d+106, 0.63556694d-81, 0.24410642d-159, 0.39124095d-230, 0.32510580d-245, 0.10257573d-122, 0.37281879d-103, 0.81561706d-52]
    !x(:) = [0.001_dp, 0.001_dp, 0.001_dp, 1.0_dp, 1.0_dp, 0.00001_dp, 0.00001_dp, 0.001_dp, 0.001_dp, 0.0001_dp, 0.0001_dp]

    x(:) = 0.1_dp

    newrap_loop: do i = 1, max_iter
        !alpha=0.0_dp

        !call calc_f(n, x, f, k)
        f = func(n, x, k)

        !print *, f

        if (maxval(abs(f)) < err) then
            print *, "loop=", i
            print *, "Maybe converged?"
            print *, exp(x)
            !print*, x
            stop
        end if

        !call jacub(n, x, jacb, k)
        jacb(:, :) = jacob(n, x, k)

        !print *, jacb

        dx = -f

        !print *, "--- loop = ", i, " ---"
        !print *, "Max F   = ", maxval(abs(f))
        !print *, "Max Jac = ", maxval(abs(jacb))
        !print *, jacb
        call dgesv(n, 1, jacb, n, ipiv, dx, n, info)

        !print *, maxval(abs(dx))

        if (info /= 0) then
            print *, "loop =", i
            print *, exp(x)
            print *, "Error: LAPACK DGESV failed with info = ", info
            stop
        end if

        !print *, func_rss(n, dx, K, func)
        call backtracking(n, x, k, dx, f, rate, alpha, func_rss) !returns lambda

        x = x + alpha * dx

        !print*, x
    end do newrap_loop

    !!$ 1.161108219318815E-041
    !!$ 1.386817472860219E-080
    !!$ 1.851877589960801E-016
    !!$ 0.212121212121212
    !!$ 0.787878787878788
    !!$ 1.253510140071907E-185
    !!$ 1.244098795583246E-239
    !!$ 3.815934069178506E-119
    !!$ 1.295473872499711E-233
    !!$ 3.624008041554619E-087
    !!$ 3.624008041554619E-087

    print *, "Loop reached maximam iteration without convergence."
    print *, exp(x)
    !print*,x
    !ifx -Warn all -O3 -xHOST -qmkl newrap_test.f90 -o newrap && ./newrap
    !ifx -O0 -g -traceback -fpe0 -qmkl newrap_test.f90 -o newrap && ./newrap
contains
    function func(num, p, KK) result(res)
        integer, intent(in) :: num
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

    function func_rss(num, p, KK) result(r)
        integer, intent(in) :: num
        real(dp), intent(in) :: p(num)
        real(dp), intent(in) :: KK(3:num - 1)
        real(dp) :: r
        real(dp) :: fu(num)

        !fu(:) = func(num, p, KK)
        r = 1 / 2 * (norm2(func(num, p, KK)))**2
    end function func_rss

    subroutine backtracking(num, p_old, KK, delta_p, fu, beta, c, fun_r) !rateを返す
        integer, intent(in) :: num
        real(dp), intent(in) :: p_old(num)
        real(dp), intent(in) :: KK(3:num - 1)
        real(dp), intent(in) :: delta_p(num)
        real(dp), intent(in) :: fu(num)
        real(dp), intent(in) :: beta
        real(dp), intent(inout) :: c
        interface
            function fun_r(number, g, h) result(res1)
                import :: dp
                implicit none(type, external)
                integer, intent(in) :: number
                real(dp), intent(in) :: g(number)
                real(dp), intent(in) :: h(number - 3)
                real(dp) :: res1
            end function fun_r
        end interface
        integer :: l
        real(dp) :: p_new(num), f_new(num), j_new(num, num)
        real(dp) :: rss_old, rss_new

        rss_old = fun_r(num, p_old, KK)
        l = 0
        c = beta**0 !Start with 1.
        do while (.not. rss_new <= rss_old - c1 * c * (norm2(fu))**2)
            !print *, c
            p_new = p_old + c * delta_p
            rss_new = fun_r(num, p_new, KK)
            f_new = func(num, p_new, KK)
            j_new = jacob(num, p_new, KK)
            if (dot_product(matmul(f_new, j_new), delta_p) <= c2 * 2.0_dp * rss_old) exit
	    l = l + 1
            c = beta**l
        end do
    end subroutine backtracking

    subroutine calc_f(num, p, ff, kk)
        integer, intent(in) :: num
        real(dp), intent(in) :: p(num)
        real(dp), intent(out) :: ff(num)
        real(dp), intent(in) :: KK(3:num - 1)

        !!$ ff(1) = KK(4)*p(4) - p(1)**2
        !!$ ff(2) = KK(5)*p(5) - p(2)**2
        !!$ ff(3) = KK(3)*p(1)*p(2) - p(3)
        !!$ ff(4) = KK(6)*p(1) - p(6)*p(11)
        !!$ ff(5) = KK(7)*p(2) - p(7)*p(11)
        !!$ ff(6) = KK(8)*p(1)**2 - p(8)*p(11)
        !!$ ff(7) = KK(9)*p(2)**2 - p(9)*p(11)
        !!$ ff(8) = KK(10)*p(1)*p(2) - p(10)*p(11)
        !!$ ff(9) = PP - sum(p)
        !!$ ff(10) = p(11) - sum(p(6:10))
        !!$ ff(11) = nn*(p(1) + p(3) + 2*p(4) + p(6) + 2*p(8) + p(10)) &
        !!$      &- no*(p(2) + p(3) + 2*p(5) + p(7) + 2*p(9) + p(10))
        ff(1) = log(KK(4)) - 2 * p(1) + p(4)
        ff(2) = log(KK(5)) - 2 * p(2) + p(5)
        ff(3) = log(KK(3)) - p(3) + p(2) + p(1)
        ff(4) = log(KK(6)) - p(6) - p(11) + p(1)
        ff(5) = log(KK(7)) - p(7) - p(11) + p(2)
        ff(6) = log(KK(8)) - p(8) - p(11) + 2 * p(1)
        ff(7) = log(KK(5)) - p(9) - p(11) + 2 * p(2)
        ff(8) = log(KK(10)) - p(10) - p(11) + p(1) + p(2)
        ff(9) = PP - sum(exp(p))
        ff(10) = exp(p(11)) - sum(exp(p(6:10)))
        ff(11) = nn * (exp(p(1)) + exp(p(3)) + 2 * exp(p(4)) + exp(p(6)) + 2 * exp(p(8)) + exp(p(10))) &
             & - no * (exp(p(2)) + exp(p(3)) + 2 * exp(p(5)) + exp(p(7)) + 2 * exp(p(9)) + exp(p(10)))
    end subroutine calc_f

    subroutine jacub(num, p, j, KK)
        integer, intent(in)  :: num
        real(dp), intent(in) ::p(num)
        real(dp), intent(out) :: j(num, num)
        real(dp), intent(in) :: KK(3:num - 1)

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
        !j(変数i, 関数j)
        ! j(1,:)  = [  -2*p(1)         ,    0.0_dp   , 0.0_dp ,  KK(4) ,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp]
        ! j(2,:)  = [   0.0_dp         ,   -2*p(2)   , 0.0_dp ,  0.0_dp,  KK(5) ,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp]
        ! j(3,:)  = [ KK(3)*p(2) ,  KK(3)*p(1) , -1.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp]
        ! j(4,:)  = [    KK(6)         ,    0.0_dp   , 0.0_dp ,  0.0_dp,  0.0_dp,  -p(11),  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  -p(6) ]
        ! j(5,:)  = [   0.0_dp         ,     KK(7)   , 0.0_dp ,  0.0_dp,  0.0_dp,  0.0_dp,  -p(11),  0.0_dp,  0.0_dp,  0.0_dp,  -p(7) ]
        ! j(6,:)  = [2*KK(8)*p(1),    0.0_dp   , 0.0_dp ,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  -p(11),  0.0_dp,  0.0_dp,  -p(8) ]
        ! j(7,:)  = [   0.0_dp         , 2*KK(9)*p(2), 0.0_dp ,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  -p(11),  0.0_dp,  -p(9) ]
        ! j(8,:)  = [KK(10)*p(2) , KK(10)*p(1) , 0.0_dp ,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  -p(11),  -p(10)]
        ! j(9,:)  = [  -1.0_dp         ,   -1.0_dp   , -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp]
        ! j(10,:) = [   0.0_dp         ,    0.0_dp   , 0.0_dp ,  0.0_dp,  0.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp,  1.0_dp]
        ! j(11,:) = [     nn         ,      -no    ,  nn-no ,   2*nn ,   -2*no,    nn  ,   -no  ,   2*nn ,  -2*no ,   nn-no,  0.0_dp]

        !j(関数i, 変数j)
        !!$ j(:,1) = [-2*p(1), 0.0_dp , KK(3)*p(2), KK(6), 0.0_dp, 2*KK(8)*p(1), 0.0_dp, KK(10)*p(2), -1.0_dp, 0.0_dp, nn]
        !!$ j(:,2) = [0.0_dp , -2*p(2), KK(3)*p(1), 0.0_dp, KK(7), 0.0_dp, 2*KK(9)*p(2), KK(10)*p(1), -1.0_dp, 0.0_dp, -nn]
        !!$ j(:,3) = [0.0_dp , 0.0_dp ,  -1.0_dp  , 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, nn-no]
        !!$ j(:,4) = [ KK(4) , 0.0_dp ,   0.0_dp  , 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 2*nn]
        !!$ j(:,5) = [0.0_dp ,  KK(5) ,   0.0_dp  , 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, -2*no ]
        !!$ j(:,6) = [0.0_dp , 0.0_dp ,   0.0_dp  , -p(11), 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, -1.0_dp, nn]
        !!$ j(:,7) = [0.0_dp , 0.0_dp ,   0.0_dp  , 0.0_dp, -p(11), 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, -1.0_dp, -no]
        !!$ j(:,8) = [0.0_dp , 0.0_dp ,   0.0_dp  , 0.0_dp, 0.0_dp, -p(11), 0.0_dp, 0.0_dp, -1.0_dp, -1.0_dp, 2*nn]
        !!$ j(:,9) = [0.0_dp , 0.0_dp ,   0.0_dp  , 0.0_dp, 0.0_dp, 0.0_dp, -p(11), 0.0_dp, -1.0_dp, -1.0_dp, -2*no]
        !!$ j(:,10) =[0.0_dp , 0.0_dp ,   0.0_dp  , 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -p(11), -1.0_dp, -1.0_dp, nn-no]
        !!$ j(:,11) =[0.0_dp , 0.0_dp ,   0.0_dp  , -p(6) , -p(7) , -p(8) , -p(9) , -p(10), -1.0_dp, 1.0_dp, 0.0_dp]
    end subroutine jacub

    function jacob(num, p, KK) result(j)
        integer, intent(in)  :: num
        real(dp), intent(in) ::p(num)
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
        !j(変数i, 関数j)
        ! j(1,:)  = [  -2*p(1)         ,    0.0_dp   , 0.0_dp ,  KK(4) ,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp]
        ! j(2,:)  = [   0.0_dp         ,   -2*p(2)   , 0.0_dp ,  0.0_dp,  KK(5) ,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp]
        ! j(3,:)  = [ KK(3)*p(2) ,  KK(3)*p(1) , -1.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp]
        ! j(4,:)  = [    KK(6)         ,    0.0_dp   , 0.0_dp ,  0.0_dp,  0.0_dp,  -p(11),  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  -p(6) ]
        ! j(5,:)  = [   0.0_dp         ,     KK(7)   , 0.0_dp ,  0.0_dp,  0.0_dp,  0.0_dp,  -p(11),  0.0_dp,  0.0_dp,  0.0_dp,  -p(7) ]
        ! j(6,:)  = [2*KK(8)*p(1),    0.0_dp   , 0.0_dp ,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  -p(11),  0.0_dp,  0.0_dp,  -p(8) ]
        ! j(7,:)  = [   0.0_dp         , 2*KK(9)*p(2), 0.0_dp ,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  -p(11),  0.0_dp,  -p(9) ]
        ! j(8,:)  = [KK(10)*p(2) , KK(10)*p(1) , 0.0_dp ,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  -p(11),  -p(10)]
        ! j(9,:)  = [  -1.0_dp         ,   -1.0_dp   , -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp]
        ! j(10,:) = [   0.0_dp         ,    0.0_dp   , 0.0_dp ,  0.0_dp,  0.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp,  1.0_dp]
        ! j(11,:) = [     nn         ,      -no    ,  nn-no ,   2*nn ,   -2*no,    nn  ,   -no  ,   2*nn ,  -2*no ,   nn-no,  0.0_dp]

        !j(関数i, 変数j)
        !!$ j(:,1) = [-2*p(1), 0.0_dp , KK(3)*p(2), KK(6), 0.0_dp, 2*KK(8)*p(1), 0.0_dp, KK(10)*p(2), -1.0_dp, 0.0_dp, nn]
        !!$ j(:,2) = [0.0_dp , -2*p(2), KK(3)*p(1), 0.0_dp, KK(7), 0.0_dp, 2*KK(9)*p(2), KK(10)*p(1), -1.0_dp, 0.0_dp, -nn]
        !!$ j(:,3) = [0.0_dp , 0.0_dp ,  -1.0_dp  , 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, nn-no]
        !!$ j(:,4) = [ KK(4) , 0.0_dp ,   0.0_dp  , 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 2*nn]
        !!$ j(:,5) = [0.0_dp ,  KK(5) ,   0.0_dp  , 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, -2*no ]
        !!$ j(:,6) = [0.0_dp , 0.0_dp ,   0.0_dp  , -p(11), 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, -1.0_dp, nn]
        !!$ j(:,7) = [0.0_dp , 0.0_dp ,   0.0_dp  , 0.0_dp, -p(11), 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, -1.0_dp, -no]
        !!$ j(:,8) = [0.0_dp , 0.0_dp ,   0.0_dp  , 0.0_dp, 0.0_dp, -p(11), 0.0_dp, 0.0_dp, -1.0_dp, -1.0_dp, 2*nn]
        !!$ j(:,9) = [0.0_dp , 0.0_dp ,   0.0_dp  , 0.0_dp, 0.0_dp, 0.0_dp, -p(11), 0.0_dp, -1.0_dp, -1.0_dp, -2*no]
        !!$ j(:,10) =[0.0_dp , 0.0_dp ,   0.0_dp  , 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -p(11), -1.0_dp, -1.0_dp, nn-no]
        !!$ j(:,11) =[0.0_dp , 0.0_dp ,   0.0_dp  , -p(6) , -p(7) , -p(8) , -p(9) , -p(10), -1.0_dp, 1.0_dp, 0.0_dp]
    end function jacob

end program main
