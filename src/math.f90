module math
    use constant
    use LAPACK95
contains

!!Calculate the distance between two points r1 and r2.
subroutine math_distance(r1,r2,d)
    implicit none
    real*8,intent(in) :: r1(:),r2(:)
    real*8 d
    integer nsize,i
    
    nsize=size(r1)
    d=0
    do i=1,nsize
        d=d+(r1(i)-r2(i))**2
    end do
    
    d=sqrt(d)
end subroutine
!calcualte the norm or \vec x, i.e. |\vec x|.
subroutine math_vecnorm(x,norm)
    implicit none
    real*8,intent(in) :: x(:)
    real*8 norm
    integer nsize,i
    norm=0
    do i=1,nsize
        norm=norm+x(i)**2
    end do
    norm =sqrt(norm)


end subroutine

!calculate \vec x \cdot \vec y 
subroutine math_dotprod(x,y,p)
    implicit none
    real*8,intent(in):: x(:),y(:)
    real*8 p
    integer nsize,i
    
    p=0
    nsize=size(x)
    do i=1,nsize
        p=p+x(i)*y(i)
    end do
    
end subroutine

!calculate dot product between two matrices. 
function math_matdotprod(x,y)
    implicit none
    real*8,intent(in):: x(:,:),y(:,:)
    real*8 p,math_matdotprod
    integer nsize1,nsize2,i,j
    
    p=0
    nsize1=size(x(:,1))
    nsize2=size(x(1,:))
    do i=1,nsize1
        do j=1,nsize2
            p=p+x(i,j)*y(i,j)
        end do
    end do

    math_matdotprod=p
    
end function

!calculate \vec x \times \vec y,only in dim=3
subroutine math_timeprod(x,y,p)
    implicit none
    real*8 x(3),y(3)
    real*8 p(3)
    integer nsize,i
    
    p(1)=x(2)*y(3)-x(3)*y(2)
    p(2)=-x(1)*y(3)+x(3)*y(1)
    p(3)=x(1)*y(2)-x(2)*y(1)
    
end subroutine



!Calculate the angle of x-y-z.
subroutine math_angle(x,y,z,theta)
    implicit none
    real*8,intent(in) :: x(:),y(:),z(:)
    real*8 theta,L1,L2,prod
    integer nsize,i
    
    nsize=size(x)
    prod=0
    
    do i=1,nsize
        prod=prod+(x(i)-y(i))*(z(i)-y(i))
    end do
    call math_distance(x,y,L1)
    call math_distance(z,y,L2)
    
    theta=acos(prod/(L1*L2))*180/pi   
end subroutine

!Calculate the dihedral angle between a-b-c-d,only in dim=3.
subroutine math_diangle(a,b,c,d,dia)    
    implicit none
    real*8 a(3),b(3),c(3),d(3),dia
    real*8 ab(3),bc(3),cd(3),g(3),h(3)
    real*8 p,norm1,norm2
    ab=a-b
    bc=b-c
    cd=c-d
    
    call math_timeprod(bc,ab,g)
    call math_timeprod(bc,cd,h)
    call math_dotprod(g,h,p)
    call math_vecnorm(g,norm1)
    call math_vecnorm(h,norm2)
    
    dia=acos(p/(norm1*norm2))*180/pi
end subroutine

!Generate two random number x1,x2 satisfy Gaussian distribution N(miu,sigma).
subroutine box_muller(x1,x2,sigma,miu)
    implicit none
    !real*8, parameter :: pi = 3.14159265358979323846
    !real*8, parameter :: hbar    = 1.0
    real*8 x1,x2,u1,u2,sigma,miu
    call RANDOM_NUMBER(u1)
    call RANDOM_NUMBER(u2)
    x1=sqrt(-2*log(u1))*cos(2*pi*u2)
    x2=sqrt(-2*log(u1))*sin(2*pi*u2)
    !x1=(x1-miu)/sigma
    !x2=(x2-miu)/sigma
    x1=x1*sigma+miu
    x2=x2*sigma+miu
end subroutine
   
function Kronecker_delta(i,j)
    integer Kronecker_delta,i,j
    if(i==j)then
        Kronecker_delta=1
    else
        Kronecker_delta=0
    end if
    
end function



function heaviside(x)
    implicit none
    real*8 x
    integer heaviside
    if(x<0)then
        heaviside=0
    else
        heaviside=1
    end if
    

end function


function squa_window(n,N0,gamma)
    implicit none
    real*8 n,gamma,squa_window
    integer N0
    !integer,external :: heaviside

    squa_window=heaviside(gamma-abs(n-N0))/(2*gamma)

end function

function tria_window(n1,n2,N,gamma)
    implicit none
    real*8 n1,n2,gamma,tria_window
    integer N

    if(N==1)then
        tria_window=2*heaviside(n1+gamma-1)*heaviside(n2+gamma)*heaviside(2-2*gamma-n1-n2)
    elseif(N==2)then
        tria_window=2*heaviside(n1+gamma)*heaviside(n2+gamma-1)*heaviside(2-2*gamma-n1-n2)
    end if

end function

function sqc_w1(n,gamma,F)
    real*8 sqc_w1,n,gamma
    integer F

    if(n>1-gamma .and. n<2-gamma)then
        sqc_w1=(2-gamma-n)**(2-F)
    else
        sqc_w1=0
    end if

end function


function sqc_w0(n,np,gamma)
    real*8 np,n,gamma
    integer F,sqc_w0

    if(np<2-2*gamma-n)then
        sqc_w0=1
    else
        sqc_w0=0
    end if

end function

function sqc_w(n,k,F,gamma)
    integer k,F,i
    real*8 sqc_w,n(F),gamma

    sqc_w=1

    do i=1,F
        if(i==k)then
            sqc_w=sqc_w*sqc_w1(n(k),gamma,F)
        else
            sqc_w=sqc_w*sqc_w0(n(k),n(i),gamma)
        end if
    end do

end function













subroutine math_sort(n,a)
    implicit none
    integer n,i,j(1)
    real*8 a(n),b(n)
    


    do i=1,n
        j=maxloc(a)
        b(n-i+1)=a(j(1))
        a(j(1))=minval(a)
    end do

    a=b

end subroutine


subroutine random_prob(n,p)
    implicit none
    integer n,j
    integer*2 compar
    real*8 p(n),p0(n-1)

    
    call random_number(p0)
    !write(*,*) p0
    call math_sort(n-1,p0)
    !write(*,*) p0
    do j=1,n
        if(j==1)then
            p(j)=p0(1)
        else if(j==n)then
            p(j)=1-p0(n-1)
        else
            p(j)=p0(j)-p0(j-1)
        end if
    end do
    
    ! call random_number(p0)
    ! P=p0/sum(p0)



end subroutine

subroutine dia_symmat(n,A,E,C)
    implicit none
    integer n,info
    real*8 A(n,n),E(n),C(n,n),work(3*n)
    
    C=A
    call dsyev('V','L',n,C,n,E,work,3*n,info)

end subroutine


subroutine dia_hermitemat(n,A,E,C)
    implicit none
    integer n,info
    complex*16 A(n,n),C(n,n),work(3*n)
    real*8 rwork(3*n),E(n)
    
    C=A
    call zheev('V','L',n,C,n,E,work,3*n,rwork,info)

end subroutine



function trace(n,A)
    implicit none
    integer n,i
    real*8 A(n,n),tA,trace

    tA=0
    do i=1,n
        tA=tA+A(i,i)
    end do

    trace=tA

end function



    
end module



! program test
!     use math
!     real*8 a(7),p2(2),p7(7),x
!     integer i

!     ! a=(/1.0_8, 4.7_8, 2.7_8, 8.3_8, 92.0_8, 0.3_8, 7.4_8/)

!     ! call math_sort(7,a) 

!     ! write(*,*) a

!     open(11,file="rp.dat")

!     do i=1,10000

!     call random_prob(2,p2)
!     call random_number(x)
!     write(11,"(4F14.8)") p2,x,1-x

!     end do




! end program