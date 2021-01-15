module constant
    
    real*8,parameter :: pi=3.1415926535897932

    complex*16,parameter :: im=(0.0_8 , 1.0_8)
    
    
    
    real*8,parameter :: hbar=1.0
    
    real*8,parameter :: kb=1.3806E-5/4.35974
    
    
    
    real*8,parameter :: au2ev=27.211386245988
    
    
     !-- unit conversion constant
    !-- a_2_b : 1(a) == a_2_b (b)
    !-- energy convertors
    !> energy convertor: au to wavenumber
    real*8, parameter :: au_2_wn = 219474.6313702
    !> energy convertor: au to electron Volt
    real*8, parameter :: au_2_eV = 27.21138602
    !> energy convertor: au to kcal per mole
    real*8, parameter :: au_2_kcalpmol = 627.509474
    !> energy convertor: au to kJ per mole
    real*8, parameter :: au_2_kJpmol = 2625.499638
    
    !-- time convertors
    !> time convertor: au to femtosecond
    real*8, parameter :: au_2_fs = 2.418884326505E-02      ! (fs)
    real*8, parameter :: au_2_ps = 2.418884326505E-05      ! (fs)
    
    
    
    !> mass convertor: amu to au 
    real*8, parameter :: amu_2_au = 1.82289E3
    
    !> length convertor: au to Angstrom 
    real*8, parameter :: au_2_angstrom = 0.529177
    
    
    
end module
