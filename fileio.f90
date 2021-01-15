subroutine read_input(idinp)
    use def
    implicit real*8(a-h,o-z)
    integer idinp

    read(idinp,*) 
    read(idinp,*) potential
    select case(trim(adjustl(potential)))
    case("HO")
        Nstate=1
        read(idinp,*) 
        read(idinp,*) omega_HO
    end select
    read(idinp,*) 
    read(idinp,*) mass
    read(idinp,*) 
    read(idinp,*) iop
    read(idinp,*) 
    read(idinp,*) jop
    read(idinp,*) 
    read(idinp,*) grid_start
    read(idinp,*) 
    read(idinp,*) grid_end
    read(idinp,*) 
    read(idinp,*) dx



end subroutine



subroutine output
    use def
    implicit real*8(a-h,o-z)


    100 format(<Ngrid+1>E18.8E3)

    select case(jop)
    case(1)

        open(15,file=filepath(1:len(trim(filepath))-4)//"_energy.out")
        open(16,file=filepath(1:len(trim(filepath))-4)//"_wavefun.out")

        do i=1,Nstate*Ngrid
            write(15,*) E(i)
        end do
        do i=1,Ngrid
            write(16,100) R(i),psi_real(i,:)
        end do



    end select




end subroutine