program plot
    use, intrinsic :: iso_fortran_env
    implicit none(type, external)

    character(len=64) :: pressure(3)
    real(real64), parameter :: P_atm(3) = [real(real64) :: 1.0, 0.1, 0.01]
    integer, parameter :: data = 60
    integer :: k, l
#ifdef USE_CUSTOM_VALUE
    real(real64), parameter :: no = real(VAL_NO, kind=real64), nn = real(VAL_NN, kind=real64)
#else
    real(real64), parameter :: no = 78.0_real64, nn = 21.0_real64
#endif
    character(len=64) :: ratio

    write (ratio, '(I0, A, I0)') nint(no), "v", nint(nn)

    pressure = period_to_p(P_atm)

    do l = 1, 3
        pressure(l) = trim(pressure(l))//"atm"
        open (unit=10, file='../output/'//trim(pressure(l))//trim(ratio)//'.dat', status='replace')
        write (10, '(I2)') data
        close (10)
    end do

contains
    function remove_char(input_char, char_to_remove) result(output_char)
        character(len=*), intent(in) :: input_char
        character(len=1), intent(in) :: char_to_remove
        character(len=:), allocatable :: output_char
        integer :: len_input, len_output
        integer :: i, j

        len_input = LEN_TRIM(input_char)
        len_output = 0

        do i = 1, len_input
            if (input_char(i:i) /= char_to_remove) then
                len_output = len_output + 1
            end if
        end do

        allocate (character(len=len_output) :: output_char)

        j = 1
        do i = 1, len_input
            if (input_char(i:i) /= char_to_remove) then
                output_char(j:j) = input_char(i:i)
                j = j + 1
            end if
        end do
    end function remove_char

    elemental function period_to_p(val) result(res)
        real(real64), intent(in) :: val
        character(len=10) :: res
        integer :: pos

        write (res, '(F10.2)') val
        res = adjustl(res)

        pos = index(res, ".")
        if (pos > 0) then
            res(pos:pos) = "p"
        end if

    end function period_to_p
end program plot

