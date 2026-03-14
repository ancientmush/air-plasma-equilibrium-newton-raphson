program plot
    use, intrinsic :: iso_fortran_env
    implicit none(type, external)

    character(len=64) :: pressure(3)
    real(real64), parameter :: P_atm(3) = [real(real64) :: 1.0, 0.1, 0.01]
    integer, parameter :: data = 60
    integer :: k, l
#ifdef USE_CUSTOM_VALUE
    real(real64), parameter :: no = VAL_NO, nn = VAL_NN !どうやって明示的に倍精度として定義する?
#else
    real(real64), parameter :: no = 78.0_real64, nn = 21.0_real64
#endif
    character(len=64) :: ratio

    write (pressure, '(F4.2, "atm")') P_atm
    write (ratio, '(I0, A, I0)') nint(no), "v", nint(nn)
    print *, ratio

    !print *, pressure
    do k = 1, 3
        call replace_deci_point(pressure(k))
    end do
    !print *, pressure

    ! do k = 1, size(pressure)
    !     pressure(k) = remove_char(pressure(k), '.')
    ! end do

    do l = 1, 3
        open (unit=10, file=trim(pressure(l))//trim(ratio)//'.dat', status='replace')
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

    subroutine replace_deci_point(input_char)
        character(len=*), intent(inout) :: input_char
        integer :: pos

        pos = index(trim(input_char), ".")

        if (pos > 0) then
            input_char(pos:pos) = "p"
        end if

        input_char = trim(input_char)

    end subroutine replace_deci_point
end program plot

