
subroutine bulid_kit
    use def
    implicit real*8(a-h,o-z)

    T=0
    do istate=0,Nstate-1
        do i=1,Ngrid
            T(i+istate,i+istate)=-2
            if(i>1)then
                T(i+istate,i-1+istate)=1
                T(i-1+istate,i+istate)=1 
            end if 
        end do
    end do

    T=-T/(2*mass*dx**2)


   


end subroutine


subroutine build_pot
    use def
    implicit real*8(a-h,o-z)

    select case(trim(adjustl(potential)))
    case("HO")
        call build_HO
    end select

end subroutine




subroutine eigensolver
    use def
    use math
    implicit real*8(a-h,o-z)
    real*8,allocatable :: c(:,:)

    ! allocate(c(Ngrid*Nstate,Ngrid*Nstate))

    call bulid_kit
    call build_pot

    call dia_symmat(Ngrid*Nstate,T+V,E,psi_real)

    


end subroutine