! potential of 1-dim 1-state system



subroutine build_HO
    use def
    implicit real*8(a-h,o-z)


    V=0
    do i=1,Ngrid
        V(i,i)=0.5*mass*omega_HO**2*(R(i)-R0_HO)**2
    end do


end subroutine



subroutine build_QP
    use def
    implicit real*8(a-h,o-z)


    V=0
    do i=1,Ngrid
        V(i,i)=0.25*R(i)**4
    end do


end subroutine



subroutine build_AHO
    use def
    implicit real*8(a-h,o-z)


    V=0
    do i=1,Ngrid
        V(i,i)=R(i)**2-0.1*R(i)**3+0.1*R(i)**4
    end do


end subroutine




subroutine build_dualHO
    use def
    implicit real*8(a-h,o-z)


    V=0
    do i=1,Ngrid
        V(i,i)=0.5*mass*omega1_dualHO**2*(R(i)+q_dualHO)**2+eps_dualHO
        V(i+Ngrid,i+Ngrid)=0.5*mass*omega2_dualHO**2*R(i)**2-eps_dualHO
        V(i,i+Ngrid)=delta_dualHO
        V(i+Ngrid,i)=delta_dualHO
    end do


end subroutine





subroutine build_3morse
    use def
    implicit real*8(a-h,o-z)
    real*8 A_morse(3,3),beta_morse(3),D_morse(3),C_morse(3),R0_morse(3),R00_morse(3,3),alpha_morse(3,3)



    select case(type_morse)
    case(1)
        D_morse=(/0.003_8, 0.004_8, 0.003_8 /)
        beta_morse=(/0.65_8, 0.6_8, 0.65_8 /)
        R0_morse=(/5.0_8, 4.0_8, 6.0_8 /)
        C_morse=(/0.0_8, 0.01_8, 0.006_8 /)
        A_morse(1,2)=0.002
        A_morse(2,3)=0.002
        A_morse(1,3)=0
        R00_morse(1,2)=3.4
        R00_morse(2,3)=4.8
        R00_morse(1,3)=0
        alpha_morse(1,2)=16
        alpha_morse(2,3)=16
        alpha_morse(1,3)=0
    case(2)
        D_morse=(/0.02_8, 0.01_8, 0.003_8 /)
        beta_morse=(/0.65_8, 0.4_8, 0.65_8 /)
        R0_morse=(/4.5_8, 4.0_8, 4.4_8 /)
        C_morse=(/0.0_8, 0.01_8, 0.02_8 /)
        A_morse(1,2)=0.005
        A_morse(2,3)=0
        A_morse(1,3)=0.005
        R00_morse(1,2)=3.66
        R00_morse(2,3)=0
        R00_morse(1,3)=3.34
        alpha_morse(1,2)=32
        alpha_morse(2,3)=0
        alpha_morse(1,3)=32
    case(3)
        D_morse=(/0.02_8, 0.02_8, 0.003_8 /)
        beta_morse=(/0.4_8, 0.65_8, 0.65_8 /)
        R0_morse=(/4.0_8, 4.5_8, 6.0_8 /)
        C_morse=(/0.02_8, 0.0_8, 0.02_8 /)
        A_morse(1,2)=0.005
        A_morse(2,3)=0
        A_morse(1,3)=0.005
        R00_morse(1,2)=3.4
        R00_morse(2,3)=0
        R00_morse(1,3)=4.97
        alpha_morse(1,2)=32
        alpha_morse(2,3)=0
        alpha_morse(1,3)=32
    end select

    do i=1,3
        do j=1,i-1
            A_morse(i,j)=A_morse(j,i)
            R00_morse(i,j)=R00_morse(j,i)
            alpha_morse(i,j)=alpha_morse(j,i)
        end do
    end do


    V=0

    do m=1,3
        do n=1,3
            do i=1,Ngrid
                if(m==n)then
                    V(i+(m-1)*Ngrid,i+(m-1)*Ngrid)=D_morse(m)*(1-exp(-beta_morse(m)*(R(i)-R0_morse(m))))**2+C_morse(m)
                else
                    V(i+(m-1)*Ngrid,i+(n-1)*Ngrid)=A_morse(m,n)*exp(-alpha_morse(m,n)*(R(i)-R00_morse(m,n))**2)
                end if      
            end do
        end do
    end do
    


end subroutine




subroutine build_tully
    use def
    implicit real*8(a-h,o-z)
    real*8 A_tully,B_tully,C_tully,D_tully,E_tully



    select case(type_tully)
    case(1)
        A_tully=0.01
        B_tully=1.6
        C_tully=0.005
        D_tully=1
        mass=2000
        R0_tully=-3.8
    case(2)
        A_tully=0.1
        B_tully=0.28
        C_tully=0.015
        D_tully=0.06
        E_tully=0.05
        mass=2000
        R0_tully=-10
        
    case(3)
        A_tully=6E-4
        B_tully=0.1
        C_tully=0.9
        mass=2000
        R0_tully=-13
       
    end select



    V=0


    select case(type_tully)
    case(1)  
        do i=1,Ngrid
            if(R(i)>0)then
                V(i,i)=A_tully*(1-exp(-B_tully*R(i)))
            else
                V(i,i)=-A_tully*(1-exp(B_tully*R(i)))
            end if  
            V(i+Ngrid,i+Ngrid)=-V(i,i)
            V(i,i+Ngrid)=C_tully*exp(-D_tully*R(i)**2)
            V(i+Ngrid,i)=V(i,i+Ngrid)
        end do            
    case(2)
        do i=1,Ngrid
            V(i,i)=0
            V(i+Ngrid,i+Ngrid)=-A_tully*exp(-B_tully*R(i)**2)+E_tully
            V(i,i+Ngrid)=C_tully*exp(-D_tully*R(i)**2)
            V(i+Ngrid,i)=V(i,i+Ngrid)
        end do   
    case(3)
        do i=1,Ngrid 
            V(i,i)=-A_tully
            V(i+Ngrid,i+Ngrid)=A_tully
            if(R(i)>0)then
                V(i,i+Ngrid)=B_tully*(1-(exp(-C_tully*R(i))-1))
            else
                V(i,i+Ngrid)=B_tully*(1+(exp(C_tully*R(i))-1))
            end if  
            V(i+Ngrid,i)=V(i,i+Ngrid)
        end do  

    end select
    


end subroutine


subroutine build_ivp
    use def
    implicit real*8(a-h,o-z)
    
    

    omega_ivp=0.02
    ! omega_ivp=0.004472
    x0_ivp=0
    v01_ivp=0.4
    de_ivp=0.1
    a_ivp=0.6
    re_ivp=2.5
    d0_ivp=0

    delta_ivp=0.2

    k_ivp=0.4
    R0_ivp=0.5
    sg_ivp=0.05
    
    

    do i=1,Ngrid
        
        V(i,i)=0.5*mass*omega_ivp**2*(R(i)-x0_ivp)**2+v01_ivp
          
        V(i+Ngrid,i+Ngrid)=de_ivp*(1-exp(-a_ivp*(R(i)-re_ivp)))**2+d0_ivp
        V(i,i+Ngrid)=k_ivp*exp(-((R(i)-R0_ivp)**2)/sg_ivp)
        V(i+Ngrid,i)=V(i,i+Ngrid)
    end do    
    


end subroutine