submodule(mod_types) subroutines
    use mod_constants, only: a,b,tsec,T
    implicit none (type, external)
    INTEGER :: i
contains
    module procedure storing
        select type(this)
        type is (atom)
            select case (trim(this%name))
            case("N")
                ALLOCATE (this%c(this%num_seg))
                this%c(1) = polynomy(tsec(1,0), tsec(1,2), a(1,1,:), b(1,1,:))
                this%c(2:3) =[(polynomy(tsec(i,1), tsec(i,2), a(1,i,:), b(1,i,:)), i=2,3)]
            case("O")
                ALLOCATE (this%c(this%num_seg))
                this%c(1) = polynomy(tsec(1,0), tsec(1,2), a(3,1,:), b(3,1,:))
                this%c(2:3) =[(polynomy(tsec(i,1), tsec(i,2), a(3,i,:), b(3,i,:)), i=2,3)]
            case("e-")
                allocate(this%c(this%num_seg))
                this%c = [(polynomy(tsec(i,1), tsec(i,2), a(11,i,:), b(11,i,:)), i=1,3)]
                case default
                print*,"there is no such physical species."
            end select
        class is (intermediate)
            select case (trim(this%name))
            case("N2")
                ALLOCATE (this%c(this%num_seg))
                this%c(1) = polynomy(tsec(1,0), tsec(1,2), a(2,1,:), b(2,1,:))
                this%c(2:3) =[(polynomy(tsec(i,1), tsec(i,2), a(2,i,:), b(2,i,:)), i=2,3)]
                allocate(this%other_used_species(this%atendee))
            case("O2")
                ALLOCATE (this%c(this%num_seg))
                this%c(1) = polynomy(tsec(1,0), tsec(1,2), a(4,1,:), b(4,1,:))
                this%c(2:3) =[(polynomy(tsec(i,1), tsec(i,2), a(4,i,:), b(4,i,:)), i=2,3)]
                allocate(this%other_used_species(this%atendee))
            case("NO")
                ALLOCATE (this%c(this%num_seg))
                this%c(1) = polynomy(tsec(1,0), tsec(1,2), a(5,1,:), b(5,1,:))
                this%c(2:3) =[(polynomy(tsec(i,1), tsec(i,2), a(5,i,:), b(5,i,:)), i=2,3)]
                allocate(this%other_used_species(this%atendee))
            case("NO+")
                ALLOCATE (this%c(this%num_seg))
                this%c = [(polynomy(tsec(i,1), tsec(i,2), a(6,i,:), b(6,i,:)), i=1,3)]
                allocate(this%other_used_species(this%atendee))
            case("N+")
                ALLOCATE (this%c(this%num_seg))
                this%c = [(polynomy(tsec(i,1), tsec(i,2), a(7,i,:), b(7,i,:)), i=1,3)]
                allocate(this%other_used_species(this%atendee))
            case("N2+")
                ALLOCATE (this%c(this%num_seg))
                this%c = [(polynomy(tsec(i,1), tsec(i,2), a(8,i,:), b(8,i,:)), i=1,3)]
                allocate(this%other_used_species(this%atendee))
            case("O+")
                ALLOCATE (this%c(this%num_seg))
                this%c = [(polynomy(tsec(i,1), tsec(i,2), a(9,i,:), b(9,i,:)), i=1,3)]
                allocate(this%other_used_species(this%atendee))
            case("O2+")
                ALLOCATE (this%c(this%num_seg))
                this%c = [(polynomy(tsec(i,1), tsec(i,2), a(10,i,:), b(10,i,:)), i=1,3)]
                allocate(this%other_used_species(this%atendee))
            case default
                print*, "there is no such phisical species."
            end select
        class default
            print*, "there is no such class."
        end select
        this%g = this%abs_gibbs(T)
        end procedure storing
end submodule subroutines
