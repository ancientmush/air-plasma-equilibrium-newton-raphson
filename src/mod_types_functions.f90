submodule(mod_types) functions
    implicit none (type, external)
contains
    module procedure enthalpy !H/(RT)
        integer :: i
        do i = 1, this%num_seg
            if (T >= this%c(i)%tmin .and. T <= this%c(i)%tmax) then
                associate(a=>this%c(i)%a, b=>this%c(i)%b)
                    enthalpy = (-a(1) / T + a(2) * log(T) + b(1)) / T + a(3) + &
                         &T * (a(4) / 2 + T * (a(5) / 3 + T * (a(6) / 4 + T * (a(7) / 5))))
                end associate
            else
                CYCLE
            end if
        end do
        end procedure enthalpy

    module procedure entropy !S/R
        integer :: i
        do i = 1, this%num_seg
            if (T >= this%c(i)%tmin .and. T <= this%c(i)%tmax) then
                associate(a=>this%c(i)%a, b=>this%c(i)%b)
                    entropy = ((-a(1) / 2) / T - a(2)) / T + a(3) * log(T) + b(2) + &
                         &T * (a(4) + T * (a(5) / 2 + T * (a(6) / 3 + T * (a(7) / 4))))
                end associate
            else
                CYCLE
            end if
        end do
    end procedure entropy

    module procedure gibbs
        gibbs = R * T * (enthalpy(this, T) - entropy(this, T))
    end procedure gibbs

    module procedure equilibrium_constant
        call this%delta_g()
        equilibrium_constant = EXP(-this%dg/(R*T))
    end procedure equilibrium_constant
end submodule functions
