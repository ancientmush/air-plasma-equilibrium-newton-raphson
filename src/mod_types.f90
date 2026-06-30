module mod_types
    use mod_constants, only: real64, R
    implicit none(type, external)
    private
    public :: atom, molecule, compound, ion, imolecule

    TYPE polynomy
        real(real64) :: tmin, tmax         !minimum & maximum temperature in the segement.
        real(real64), dimension(7) :: a
        real(real64), dimension(2) :: b
    END TYPE polynomy

    type reactants
        real(real64), allocatable :: g_used(:)
    end type reactants

    type atom ! for electron as well
        character(LEN=:), allocatable :: name
        integer :: num_seg = 3 !Default value as it's general.
        type(polynomy), allocatable :: c(:)
        real(real64), allocatable :: g(:)
    contains
        procedure, pass :: abs_gibbs => gibbs
        procedure, pass :: data => storing
    end type atom

    type, abstract, extends(atom) :: intermediate
        !intermediate class to implement compound species
        integer :: atendee
        type(reactants), allocatable :: other_used_species(:)
        real(real64), allocatable :: dg(:)
        real(real64), allocatable :: kp(:)
    contains
        procedure, pass :: kps => equilibrium_constant
        procedure(gibbs_change), deferred, pass :: delta_g
    end type intermediate

    abstract interface
        elemental subroutine gibbs_change(this)
            import :: intermediate
            implicit none(type, external)
            class(intermediate), intent(inout) :: this
        end subroutine gibbs_change
    end interface

    interface
        module elemental real(real64) function enthalpy(this, T) ! Calculating Gibbs
            implicit none(type, external)
            class(atom), INTENT(IN) :: this
            real(real64), INTENT(IN) :: T
        end function enthalpy
        module elemental real(real64) function entropy(this, T) ! Calculating Gibbs
            implicit none(type, external)
            class(atom), INTENT(IN) :: this
            real(real64), INTENT(IN) :: T
        end function entropy
        module elemental real(real64) function gibbs(this, T)
            implicit none(type, external)
            class(atom), intent(in) :: this
            real(real64), intent(in) :: T
        end function gibbs
        module function equilibrium_constant(this, T)
            implicit none(type, external)
            class(intermediate), intent(inout) :: this
            real(real64), intent(in) :: T(:)
            real(real64) :: equilibrium_constant(size(T))
        end function equilibrium_constant
        module subroutine storing(this) ! store each species' values
            use mod_constants, only: a, b, T
            implicit none(type, external)
            class(atom), intent(inout) :: this
        end subroutine storing
    end interface

    type, extends(intermediate) :: molecule !for the likes of N2
    contains
        procedure :: delta_g => clevage
    end type molecule

    type, extends(intermediate) :: compound !NO
    contains
        procedure :: delta_g => dissociation
    end type compound

    type, extends(intermediate) :: ion ! for the likes of O+
    contains
        procedure :: delta_g => ionization
    end type ion

    type, extends(molecule) :: imolecule ! for the likes of O2+
    contains
        procedure :: delta_g => ion_clevage
    end type imolecule

contains
    elemental subroutine clevage(this) !N2, O2
        class(molecule), intent(inout) :: this
        this%dg = 2 * this%other_used_species(1)%g_used - this%g
    end subroutine clevage
    elemental subroutine dissociation(this) !NO
        class(compound), intent(inout) :: this
        integer :: i
        this%dg = this%g
        do i = 1, this%atendee
            this%dg = this%dg - this%other_used_species(i)%g_used
        end do
    end subroutine dissociation
    elemental subroutine ionization(this) !O+, N+, NO+
        class(ion), intent(inout) :: this
        integer :: i
        this%dg = this%g + this%other_used_species(1)%g_used
        do i = 2, this%atendee
            this%dg = this%dg - this%other_used_species(i)%g_used
        end do
    end subroutine ionization
    elemental subroutine ion_clevage(this) !N2+, O2+
        class(imolecule), intent(inout) :: this
        this%dg = this%g + this%other_used_species(1)%g_used - 2 * this%other_used_species(2)%g_used
    end subroutine ion_clevage
end module mod_types

