
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
    case("QP")
        call build_QP
    case("dualHO")
        call build_dualHO
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

    call dia_symmat(Ngrid*Nstate,T+V,E,eigenwf)


end subroutine



subroutine initial_wf
    use def
    use constant
    implicit real*8(a-h,o-z)

    select case(trim(adjustl(potential)))
    case("dualHO")
        psi=0
        do i=1,Ngrid
            psi(i)=(omega1_dualHO/pi)**0.25*exp(-omega1_dualHO*(R(i)+q_dualHO)**2/2)*exp(im*P0*(R(i)+q_dualHO))
        end do
    end select


end subroutine



subroutine DVRpropagator
    use def
    use math
    use constant
    implicit real*8(a-h,o-z)
    real*8,allocatable :: LDA(:,:)



    call eigensolver
    allocate(LDA(Ngrid*Nstate,Ngrid*Nstate))
    LDA=0
    do i=1,Ngrid*Nstate 
        LDA(i,i)=E(i)
    end do

    propagator=matmul(eigenwf,matmul(exp(-im*LDA*dt),transpose(eigenwf)))

    deallocate(LDA)

    call initial_wf

    nstep=ttot/dt 
    allocate(rho(Nstate,nstep))
    allocate(time(nstep))

    do i=1,Nstate
        rho(i,1)=sum(real(psi(1+(i-1)*Ngrid:Ngrid+(i-1)*Ngrid))**2+imag(psi(1+(i-1)*Ngrid:Ngrid+(i-1)*Ngrid))**2)*dx
    end do
    time(1)=0

    do istep=2,Nstep
        time(istep)=(istep-1)*dt
        psi=matmul(propagator,psi)
        do i=1,Nstate
            rho(i,istep)=sum(real(psi(1+(i-1)*Ngrid:Ngrid+(i-1)*Ngrid))**2+imag(psi(1+(i-1)*Ngrid:Ngrid+(i-1)*Ngrid))**2)*dx
        end do
    end do
        




    

end subroutine