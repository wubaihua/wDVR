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

    real*8,allocatable :: T(:,:),V(:,:) ! kinectic & potential matrix for DVR
    real*8,allocatable :: psi_real(:,:) ! DVR real wave function
    complex*16,allocatable :: psi(:,:) ! DVR wave function
    real*8,allocatable :: E(:)

    real*8,allocatable :: pot(:,:) ! multi-state potential matrix
   

    complex*16,allocatable :: propagator(:,:) ! DVR wave function


    real*8 mass



    real*8 omega_HO



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

        allocate(psi(Nstate*Ngrid,Nstate*Ngrid))
        allocate(psi_real(Nstate*Ngrid,Nstate*Ngrid))
        allocate(propagator(Nstate*Ngrid,Nstate*Ngrid))
        allocate(E(Nstate*Ngrid))

    end subroutine













end module