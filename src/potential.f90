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
        D_morse=(/0.003_8, 0.004_8, 0.005_8 /)
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