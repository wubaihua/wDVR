module def
    implicit real*8(a-h,o-z)

    character*200 filepath
    character*20 potential

    integer iop
    integer jop ! index of job type

    real*8 grid_start,grid_end,dx
    integer Ngrid
    real*8,allocatable :: R(:)

    integer Nstate

    real*8 beta 
    real*8 vpot,kinetic

    real*8,allocatable :: T(:,:),V(:,:) ! kinectic & potential matrix for DVR
    real*8,allocatable :: eigenwf(:,:) ! DVR eigen function
    complex*16,allocatable :: psi(:) ! DVR time-dependent wave function
    real*8,allocatable :: E(:)

    real*8,allocatable :: pot(:,:) ! multi-state potential matrix


    real*8 parti_fun
   

    complex*16,allocatable :: propagator(:,:) ! DVR wave function
    real*8 ttot,dt
    real*8 P0

    real*8,allocatable :: rho(:,:)

    real*8,allocatable :: time(:)

    complex*16,allocatable :: corre_fun(:,:)

    integer nstep

    real*8 mass

    real*8 unittrans

    ! parameters for harmonic oscillator
    real*8 omega_HO,R0_HO




    ! parameters for 2-state dual harmonic oscillators
    real*8 omega1_dualHO,omega2_dualHO,eps_dualHO,delta_dualHO, q_dualHO

    ! parameters for 3-state morse potential
    integer type_morse

    ! parameters for 2-state tully model
    integer type_tully
    real*8 R0_tully,gamma_tully
    real*8,allocatable :: pop_f(:,:),pop_b(:,:)

    ! parameters for 15-state model
    real*8 Re_15s,Dg_15s,alpha_15s,eta_15s,De_15s,Dc_15s,lambda_15s,R_cross_15s(15),sigma_15s


contains

    subroutine initialization
        implicit real*8(a-h,o-z)

        Ngrid=(grid_end-grid_start)/dx
        allocate(R(Ngrid))
        do i=1,Ngrid
            R(i)=grid_start+i*dx
        end do

        allocate(T(Nstate*Ngrid,Nstate*Ngrid))
        allocate(V(Nstate*Ngrid,Nstate*Ngrid))

        allocate(psi(Nstate*Ngrid))
        allocate(eigenwf(Nstate*Ngrid,Nstate*Ngrid))
        allocate(propagator(Nstate*Ngrid,Nstate*Ngrid))
        allocate(E(Nstate*Ngrid))

    end subroutine













end module