program main
    use globalFileName
    use ga_module
    implicit none
    real(8) :: alpha, beta, gamma, eta
    real(8), allocatable :: array(:)
    integer :: arraySize, i

    allocate( array(indv_size) )
    call mutation(array)

    do i = 1, indv_size
        print *, array(i)
    end do

    !call heat_conduction(1.0d-5, 1.0d-5, 1.0d-4, 5.0d2, highT, lowT)
    !arraySize = size( highT )


contains


end program main
