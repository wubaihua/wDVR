subroutine read_input(idinp)
    use def
    implicit real*8(a-h,o-z)
    integer idinp

    read(idinp,*) 
    read(idinp,*) potential
    read(idinp,*) 
    select case(trim(adjustl(potential)))
    case("HO")
        Nstate=1
        read(idinp,*) 
        read(idinp,*) omega_HO
        read(idinp,*) 
        read(idinp,*) R0_HO
    case("QP")
        Nstate=1
    case("dualHO")
        Nstate=2
        read(idinp,*) 
        read(idinp,*) omega1_dualHO,omega2_dualHO
        read(idinp,*)
        read(idinp,*) eps_dualHO,delta_dualHO
        read(idinp,*) 
        read(idinp,*) q_dualHO
    case("3morse")
        Nstate=3
        read(idinp,*) 
        read(idinp,*) type_morse
    end select
    read(idinp,*) 
    read(idinp,*) mass
    read(idinp,*) 
    read(idinp,*) iop
    read(idinp,*) 
    read(idinp,*) jop
    select case(jop)
    case(2)
        read(idinp,*) 
        read(idinp,*) dt
        read(idinp,*) 
        read(idinp,*) ttot
        read(idinp,*) 
        read(idinp,*) P0
    end select
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
    101 format(<Nstate+1>E18.8E3)

    select case(jop)
    case(1)

        open(15,file=filepath(1:len(trim(filepath))-4)//"_energy.out")
        open(16,file=filepath(1:len(trim(filepath))-4)//"_wavefun.out")

        do i=1,Nstate*Ngrid
            write(15,*) E(i)
        end do
        do i=1,Ngrid
            write(16,100) R(i),eigenwf(i,:)
        end do

    case(2)
        open(20,file=filepath(1:len(trim(filepath))-4)//"_population.out")

        do i=1,nstep
            write(20,101) time(i),rho(:,i)
        end do
    end select




end subroutine