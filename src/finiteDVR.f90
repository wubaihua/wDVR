
subroutine bulid_kit
    use def
    implicit real*8(a-h,o-z)

    T=0
    do istate=0,Nstate-1
        do i=1,Ngrid
            T(i+istate*Ngrid,i+istate*Ngrid)=-2
            if(i>1)then
                T(i+istate*Ngrid,i-1+istate*Ngrid)=1
                T(i-1+istate*Ngrid,i+istate*Ngrid)=1 
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
    case("AHO")
        call build_AHO
    case("dualHO")
        call build_dualHO
    case("3morse")
        call build_3morse
    case("tully")
        call build_tully
    case("ivp")
        call build_ivp
    case("15s")
        call build_15s
    end select

end subroutine




subroutine eigensolver
    use def
    use math
    implicit real*8(a-h,o-z)
    real*8,allocatable :: c(:,:)
    REAL*8 V_inte(Ngrid*Nstate)
    ! real*8 parti_fun

    ! allocate(c(Ngrid*Nstate,Ngrid*Nstate))

    call bulid_kit
    call build_pot

    call dia_symmat(Ngrid*Nstate,T+V,E,eigenwf)

    ! allocate(c(Ngrid*Nstate,Ngrid*Nstate))
    ! c=matmul(eigenwf,transpose(eigenwf))
    ! do i=1,Ngrid*Nstate
    !     print*, c(i,i)
    ! end do
    ! print*, sum(c)
    
    do i=1,Ngrid*Nstate
        V_inte(i)=V(i,i)
    end do

    Vpot=0
    kinetic=0

    parti_fun=0

    ! write(*,*) sum(eigenwf(:,1)**2)*dx
    ! print*, beta

    do i=1,Ngrid*Nstate
        parti_fun=parti_fun+exp(-beta*E(i))
        ! eigenwf(:,i)=eigenwf(:,i)/(sqrt(sum(eigenwf(:,i)**2)*dx))
        ! write(*,*) sum(eigenwf(:,i)**2)*dx
        Vpot=Vpot+exp(-beta*E(i))*sum(eigenwf(:,i)**2*V_inte(:))*dx/(sum(eigenwf(:,i)**2)*dx)
        kinetic=kinetic+exp(-beta*E(i))*sum(eigenwf(:,i)*matmul(T,eigenwf(:,i)))*dx/(sum(eigenwf(:,i)**2)*dx)
    end do

    Vpot=Vpot/parti_fun
    kinetic=kinetic/parti_fun
    


end subroutine



subroutine initial_wf
    use def
    use constant
    implicit real*8(a-h,o-z)
    real*8 re_3morse

    select case(trim(adjustl(potential)))
    case("dualHO")
        psi=0
        do i=1,Ngrid
            psi(i)=(omega1_dualHO/pi)**0.25*exp(-omega1_dualHO*(R(i)+q_dualHO)**2/2)*exp(im*P0*(R(i)+q_dualHO))
        end do
    case("3morse")
        if(type_morse==1)then
            re_3morse=2.9
        else if(type_morse==2)then
            re_3morse=3.3
        else if(type_morse==3)then
            re_3morse=2.1
        end if  
        psi=0
        do i=1,Ngrid
            psi(i)=(mass*5E-3/pi)**0.25*exp(-mass*5E-3*(R(i)-re_3morse)**2/2)*exp(im*P0*R(i))
        end do
    case("tully")
        gamma_tully=1
        psi=0
        do i=1,Ngrid
            psi(i)=(gamma_tully/pi)**0.25*exp(-gamma_tully/2*(R(i)-R0_tully)**2)*exp(im*P0*(R(i)-R0_tully))
        end do
    case("ivp")
        x0_ivp=0
        omega_ivp=0.02
        psi=0
        do i=1,Ngrid
            psi(i)=(mass*omega_ivp/pi)**0.25*exp(-mass*omega_ivp/2*(R(i)-x0_ivp)**2)*exp(im*P0*(R(i)-x0_ivp))
        end do
    case("15s")
        

        sigma_15s=sqrt(mass*2*Dg_15s*alpha_15s**2)
        do i=1,Ngrid
            psi(i)=(sigma_15s/pi)**0.25*exp(-sigma_15s/2*(R(i)-13)**2)*exp(im*P0*(R(i)-13))
        end do
    end select


end subroutine



subroutine DVRpropagator
    use def
    use math
    use constant
    implicit real*8(a-h,o-z)
    real*8,allocatable :: LDA(:,:)
    complex*16,allocatable :: elda(:,:)



    call eigensolver
    ! allocate(LDA(Ngrid*Nstate,Ngrid*Nstate))
    allocate(eLDA(Ngrid*Nstate,Ngrid*Nstate))
    ! LDA=0
    elda=(0.0_8,0.0_8)
    do i=1,Ngrid*Nstate 
        ! LDA(i,i)=E(i)
        elda(i,i)=exp(-im*E(i)*dt)
    end do

    

    propagator=matmul(eigenwf,matmul(eLDA,transpose(eigenwf)))

    deallocate(eLDA)

    call initial_wf

    nstep=ttot/dt 
    allocate(rho(Nstate,nstep))
    allocate(time(nstep))
    ! if(trim(adjustl(potential))=="tully")then
    !     allocate(pop_f(Nstate,nstep))
    !     allocate(pop_b(Nstate,nstep))
    !     ! do i=1,Ngrid
    !     !     if(R(i)<0 .and. R(i+1)>0)then
    !     !         i_zero=i 
    !     !         exit
    !     !     end if
    !     ! end do
    ! end if
    do i=1,Nstate
        rho(i,1)=sum(real(psi(1+(i-1)*Ngrid:Ngrid+(i-1)*Ngrid))**2+imag(psi(1+(i-1)*Ngrid:Ngrid+(i-1)*Ngrid))**2)*dx
    end do
    ! if(trim(adjustl(potential))=="tully")then
    !     ! pop_b(1,1)=sum(real(psi(1:i_zero))**2+imag(psi(1:i_zero))**2)*dx
    !     ! pop_b(2,1)=sum(real(psi(1+Ngrid:i_zero+Ngrid))**2+imag(psi(1:i_zero))**2)*dx
    !     ! pop_f(1,1)=sum(real(psi(i_zero+1:Ngrid))**2+imag(psi(i_zero+1:Ngrid))**2)*dx
    !     ! pop_f(2,1)=sum(real(psi(i_zero+1+Ngrid:2*Ngrid))**2+imag(psi(i_zero+1+Ngrid:2*Ngrid))**2)*dx
    !     call cal_fb_pop(dx,Ngrid,R,psi(1:Ngrid),pop_f(1,1),pop_b(1,1))
    !     call cal_fb_pop(dx,Ngrid,R,psi(1+Ngrid:2*Ngrid),pop_f(2,1),pop_b(2,1))
    ! end if
    time(1)=0

    do istep=2,Nstep
        time(istep)=(istep-1)*dt
        psi=matmul(propagator,psi)
        do i=1,Nstate
            rho(i,istep)=sum(real(psi((1+(i-1)*Ngrid):(Ngrid+(i-1)*Ngrid)))**2+imag(psi(1+(i-1)*Ngrid:Ngrid+(i-1)*Ngrid))**2)*dx
        end do
        ! if(trim(adjustl(potential))=="tully")then
        !     ! pop_b(1,istep)=sum(real(psi(1:i_zero))**2+imag(psi(1:i_zero))**2)*dx
        !     ! pop_b(2,istep)=sum(real(psi(1+Ngrid:i_zero+Ngrid))**2+imag(psi(1:i_zero))**2)*dx
        !     ! pop_f(1,istep)=sum(real(psi(i_zero+1:Ngrid))**2+imag(psi(i_zero+1:Ngrid))**2)*dx
        !     ! pop_f(2,istep)=sum(real(psi(i_zero+1+Ngrid:2*Ngrid))**2+imag(psi(i_zero+1+Ngrid:2*Ngrid))**2)*dx
        !     call cal_fb_pop(dx,Ngrid,R,psi(1:Ngrid),pop_f(1,istep),pop_b(1,istep))
        !     call cal_fb_pop(dx,Ngrid,R,psi(1+Ngrid:2*Ngrid),pop_f(2,istep),pop_b(2,istep))
        ! end if
    end do
        
    



    

end subroutine




subroutine cal_fb_pop(dx,Ngrid,R,psi,pop_f,pop_b)
    use constant
    implicit real*8(a-h,o-z)
    integer Ngrid
    real*8 R(Ngrid),pop_f,pop_b,dx,dp,P_start,P_end
    complex*16 psi(Ngrid)
    real*8,allocatable :: P(:)
    complex*16,allocatable :: psi_p(:)

    dp=dx
    P_start=R(1)
    P_end=R(Ngrid)
    ngrid_p=(P_end-P_start)/dp
    allocate(P(ngrid_p))
    allocate(psi_p(ngrid_p))

    do i=1,ngrid_p
        P(i)=P_start+i*dp
        psi_p(i)=sum(psi(:)*exp(-im*P(i)*R(:)))*dx/sqrt(2*pi)
    end do

    ! open(34,file="psi_p.dat")
    ! do i=1,ngrid_p
    !     write(34,*) p(i),real(psi_p(i))**2+imag(psi_p(i))**2
    ! end do
    ! stop

    do i=1,Ngrid_p
        if(P(i)<0 .and. P(i+1)>0)then
            i_zero=i 
            exit
        end if
    end do

    pop_b=sum(real(psi_p(1:i_zero))**2+imag(psi_p(1:i_zero))**2)*dp
    pop_f=sum(real(psi_p(1+i_zero:Ngrid_p))**2+imag(psi_p(1+i_zero:Ngrid_p))**2)*dp

    deallocate(P)
    deallocate(psi_p)


end subroutine





subroutine DVRcorrefun
    use def
    use math
    use constant
    use omp_lib
    implicit real*8(a-h,o-z)
    real*8,allocatable :: x2_grid(:)
    complex*16,allocatable :: p_grid(:,:)

    allocate(x2_grid(Ngrid))
    allocate(p_grid(Ngrid,Ngrid))

    x2_grid=R**2
    do i=1,Ngrid 
        p_grid(i,i)=-1.0
        if(i<Ngrid)p_grid(i,i+1)=1.0
    end do
    p_grid=p_grid*(-im)/dx 

    call eigensolver
    nstep=ttot/dt 
    allocate(corre_fun(nstep,20))
    allocate(time(nstep))
    corre_fun=(0.0_8,0.0_8)
    time(1)=0

    !$OMP PARALLEL DO PRIVATE(n,m) SHARED(corre_fun,E,Ngrid,beta,eigenwf,p_grid,dx,x2_grid)
    do n=1,Ngrid
        do m=1,Ngrid
            if(n==m)then
                corre_fun(1,1)=corre_fun(1,1)+exp(-beta*E(n))*(sum(eigenwf(:,n)*matmul(p_grid,eigenwf(:,m)))*dx)*(sum(eigenwf(:,m)*matmul(p_grid,eigenwf(:,n)))*dx)/(sum(eigenwf(:,n)**2)*dx*sum(eigenwf(:,m)**2)*dx)
            else
                corre_fun(1,1)=corre_fun(1,1)+(exp(-beta*E(m))-exp(-beta*E(n)))/(beta*(E(n)-E(m)))*(sum(eigenwf(:,n)*matmul(p_grid,eigenwf(:,m)))*dx)*(sum(eigenwf(:,m)*matmul(p_grid,eigenwf(:,n)))*dx)/(sum(eigenwf(:,n)**2)*dx*sum(eigenwf(:,m)**2)*dx)
            end if
            corre_fun(1,2)=corre_fun(1,2)+exp(-beta*E(n))*(sum(eigenwf(:,n)*eigenwf(:,m)*x2_grid)*dx)**2/(sum(eigenwf(:,n)**2)*dx*sum(eigenwf(:,m)**2)*dx)
            ! write(*,*) 'test1'
        end do
    end do
    !$OMP END PARALLEL DO
    
    ! corre_fun=corre_fun/parti_fun 
    ! write(*,*) corre_fun(1,2)

    !$OMP PARALLEL DO PRIVATE(n,m,istep) SHARED(corre_fun,E,Ngrid,beta,eigenwf,p_grid,dx,x2_grid,time,dt)
    do istep=2,Nstep
        time(istep)=(istep-1)*dt
        ! psi=matmul(propagator,psi)
        do n=1,Ngrid
            do m=1,Ngrid
                if(n==m)then
                    corre_fun(istep,1)=corre_fun(istep,1)+exp(-beta*E(n))*(sum(eigenwf(:,n)*matmul(p_grid,eigenwf(:,m)))*dx)*(sum(eigenwf(:,m)*matmul(p_grid,eigenwf(:,n)))*dx)/(sum(eigenwf(:,n)**2)*dx*sum(eigenwf(:,m)**2)*dx)*exp(im*(E(m)-E(n))*time(istep))
                else
                    corre_fun(istep,1)=corre_fun(istep,1)+(exp(-beta*E(m))-exp(-beta*E(n)))/(beta*(E(n)-E(m)))*(sum(eigenwf(:,n)*matmul(p_grid,eigenwf(:,m)))*dx)*(sum(eigenwf(:,m)*matmul(p_grid,eigenwf(:,n)))*dx)/(sum(eigenwf(:,n)**2)*dx*sum(eigenwf(:,m)**2)*dx)*exp(im*(E(m)-E(n))*time(istep))
                end if
                corre_fun(istep,2)=corre_fun(istep,2)+exp(-beta*E(n)+im*(E(m)-E(n))*time(istep))*(sum(eigenwf(:,n)*eigenwf(:,m)*x2_grid)*dx)**2/(sum(eigenwf(:,n)**2)*dx*sum(eigenwf(:,m)**2)*dx)
            end do
        end do 
        
    end do
    !$OMP END PARALLEL DO
        
    corre_fun=corre_fun/parti_fun
   
end subroutine