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