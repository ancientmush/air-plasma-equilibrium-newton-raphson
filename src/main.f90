program main
    use mod_constants, only: real64, R, T, init
    use mod_types, only: atom, molecule, ion, compound, imolecule
    use newton_raphson, only: main_loop

    implicit none(type, external)
    integer :: i, j, l
    integer, parameter :: k = 11
    integer, parameter :: max_iter = 10000
    real(real64), parameter :: erro = 1.0d-10
    real(real64), parameter :: Patm(3) = [real(real64) :: 1.0, 0.1, 0.01]
    real(real64), parameter :: ratio_N = 78.0_real64, ratio_O = 21.0_real64
    real(real64), allocatable :: kp(:, :)
    real(real64), allocatable :: p(:, :)

    type(atom) :: n, o, e
    type(molecule) :: n2, o2
    type(ion) :: np, op, nop
    type(imolecule) :: op2, np2
    type(compound) :: no

    ! Storing data from mod_constants to types defined in mod_types.
    call init()
    call data_storing()

    ! Calculating equilibrium constants and write into a file "k_test.dat".
    allocate (kp(size(T), k - 3))
    kp(:, 1) = no%kps(T)
    kp(:, 2) = o2%kps(T)
    kp(:, 3) = n2%kps(T)
    kp(:, 4) = op%kps(T)
    kp(:, 5) = np%kps(T)
    kp(:, 6) = op2%kps(T)
    kp(:, 7) = np2%kps(T)
    kp(:, 8) = nop%kps(T)

    open (10, file='output/k_test.dat', status='replace')
    do i = 1, size(T)
        write (10, '(f8.2,11E16.8)') T(i), (kp(i, j), j=1, 8)
    end do
    close (10)

    do l = 1, size(Patm)
        allocate (p(size(T), k))
        open (30, file='output/pressure_test.dat', status='replace')
        do i = 1, size(T)
            p(i, :) = Patm(l)
            print *, 'T=', T(i)
            call main_loop(k, max_iter, erro, Patm(l), ratio_N, ratio_O, p(i, :), Kp(i, :))
            write (30, '(f8.2, 1x, 11E16.8E3)') T(i), (exp(p(i, j)) / Patm(l), j=1, 11)
        end do
        close (30)
        deallocate (p)
    end do

    deallocate (kp)

contains
    subroutine data_storing()
        n = atom(name='N')
        o = atom(name='O')
        n2 = molecule(name='N2', atendee=1)
        o2 = molecule(name='O2', atendee=1)
        e = atom(name='e-')
        np = ion(name='N+', atendee=2)
        op = ion(name='O+', atendee=2)
        op2 = imolecule(name='O2+', atendee=2)
        np2 = imolecule(name='N2+', atendee=2)
        no = compound(name='NO', atendee=2)
        nop = ion(name='NO+', atendee=3)
        call n%data()
        call o%data()
        call n2%data()
        call o2%data()
        call e%data()
        call np%data()
        call op%data()
        call np2%data()
        call op2%data()
        call no%data()
        call nop%data()
        n2%other_used_species(1)%g_used = n%g
        o2%other_used_species(1)%g_used = o%g
        no%other_used_species(1)%g_used = n%g
        no%other_used_species(2)%g_used = o%g
        !(1) must be of e- for type(ion)
        op%other_used_species(1)%g_used = e%g
        op%other_used_species(2)%g_used = o%g
        np%other_used_species(1)%g_used = e%g
        np%other_used_species(2)%g_used = n%g
        op2%other_used_species(1)%g_used = e%g
        op2%other_used_species(2)%g_used = o%g
        np2%other_used_species(1)%g_used = e%g
        np2%other_used_species(2)%g_used = n%g
        nop%other_used_species(1)%g_used = e%g
        nop%other_used_species(2)%g_used = o%g
        nop%other_used_species(3)%g_used = n%g
    end subroutine data_storing
end program main
