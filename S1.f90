! potential of 1-dim 1-state system



subroutine build_HO
    use def
    implicit real*8(a-h,o-z)


    V=0
    do i=1,Ngrid
        V(i,i)=0.5*mass*omega_HO**2*R(i)**2
    end do


end subroutine