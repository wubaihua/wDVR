program wDVR
    use def
    implicit real*8(a-h,o-z)
    



    call getarg(1,filepath)
    
    open(10, file=trim(filepath), status='old')  

    call read_input(10)

    call initialization


    select case(jop)
    case(1)
        call eigensolver
    case(2)
        if(Nstate>1)then
            call DVRpropagator
        else
            call DVRcorrefun
        end if
    end select

    call output




end program