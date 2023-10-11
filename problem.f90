module globalFileName
  implicit none
  character(len=30) :: csvName, binName
end module globalFileName

module myProblem
  implicit none
  interface is_inf
      module procedure is_inf4
      module procedure is_inf8
  end interface

contains
  subroutine heat_conduction(alpha, beta, gamma, eta, highT, lowT)
    implicit none
    real(8), intent(in)  :: alpha, beta, gamma, eta ! fitting parameters
    real(8), allocatable, intent(out) :: highT(:), lowT(:) ! temperature array for fitting
    integer, parameter :: atLowT = 123 ! measurement point of lowT
    integer, parameter :: n = 200 ! array size
    integer, parameter :: max_steps = 800000
    real(8), parameter :: L = 2.45d-2, dx = L/n, dt = 0.0001
    real(8), parameter :: T0 = 300       ! ambient temperature [K]
    real(8), parameter :: sigma  = 5.670374419d-8 ! Stefan-Boltzmann constant [W/m2/K4]
    real(8), parameter :: Wmax  = 100 ! Maximum heater power [W]
    real(8), parameter :: t_Vincr  = 20 ! duration of heater power increase [s]
    real(8), parameter :: t_heat  = 60 ! Heating duration [s]
    integer, parameter :: nOut = 10000 ! Output the results every nOut time step
    real(8) :: ht
    real(8) :: Th_coef ! gradient of Th [K/s]
    real(8) :: Th ! The highest temperature [K]
    real(8) :: T04, u4 ! Temperature squared
    real(8) :: conductivity(n), HeatTrans(n), radiation(n), HeatGain(n)
    real(8) :: u(n), new_u(n)
    integer :: i, t

    allocate( highT(max_steps/nOut) )
    allocate( lowT(max_steps/nOut) )
    do i = 1, max_steps/nOut
        highT(i) = T0
        lowT(i)  = T0
    end do

    ! Setting initial conditions
    Th_coef = Wmax / t_Vincr / 2.0
    Th = Wmax / 2.0
    do i = 1, n
        u(i) = T0
        conductivity(i) = 0
        HeatTrans(i)    = 0
        HeatGain(i)     = 0
        radiation(i)    = 0
    end do

    ! Time stepping loop
    do t = 1, max_steps
        if ( Th < Wmax ) then
            Th = Th_coef * dt * t + ( Wmax / 2.0 )
        else
            Th = Wmax
        end if
        ! Solving the heat conduction differential equation
        do i = 2, n-1
            ! conductivity term
            conductivity(i) = alpha / dx**2 * ( u(i+1) - 2*u(i) + u(i-1) )
            ! heat transfer term
            ht = heat_trans_coef(u(i)-T0, beta)
            HeatTrans(i) = ht * ( u(i) - T0 )
            ! radiatoin term
            u4  = u(i) * u(i) * u(i) * u(i)
            T04 = T0 * T0 * T0 * T0
            radiation(i) = sigma * gamma * ( u4-T04 )

            new_u(i) = u(i) + dt * ( conductivity(i) - radiation(i) - HeatTrans(i) )
        end do
        ! adiabatic boundary
        new_u(1) = new_u(2)
        if ( dt*t <= t_heat ) then
            ! Heat flux (heat source) boundary
            HeatGain(n) = eta * Th
            new_u(n) = new_u(n-1) + eta * Th * dx / alpha
        else
            ! adiabatic boundary
            HeatGain(n) = 0
            new_u(n) = new_u(n-1)
        end if

        ! Updating the temperature array
        u(:) = new_u(:)

        ! Output temperature distribution every nOut steps to a CSV file
        if ( mod(t, nOut) == 0 ) then
            !call write_csv(u, conductivity, radiation, HeatTrans, HeatGain, n, t, dx)
            highT(t/nOut) = u(n)
            lowT(t/nOut) = u(atLowT)
        end if
    end do
    ! save fitting parameters
    call save_data(alpha, beta, gamma, eta, highT, lowT, max_steps/nOut, dt)

  end subroutine heat_conduction

  real(8) function heat_trans_coef(dT, C)
    implicit none
    real(8), parameter :: L = 4.5d-3 ! heater diameter [m]
    real(8) :: dT
    real(8) :: C

    heat_trans_coef = 2.51 * C * (dT/L)**(0.25)
  end function heat_trans_coef

  subroutine save_data(alpha, beta, gamma, eta, highT, lowT, len, dt)
    use globalFileName
    implicit none
    real(8) :: alpha, beta, gamma, eta
    real(8) :: dt
    real(8), dimension(len) :: highT, lowT
    integer :: i, len
    integer, parameter :: iounit = 11
    integer :: ierr
    real(8), dimension(len) :: t

    t(1) = 0
    do i = 1, len
        highT(i) = highT(i) - 273.15
        lowT(i)  = lowT(i) - 273.15
        if (i /= 1) then
            t(i) = t(i-1) + dt
        end if
    end do

    ! Creating a filename using the timestep
    write(binName, '(A, A,I8.8,A)') 'BIN/output_.bin'
    binName = trim( adjustl(binName) )

    open(unit=iounit, file=binName, form='unformatted')
    if (ierr /= 0) then
        print *, 'Error opening file:', binName
        return
    end if
    write(iounit) alpha
    write(iounit) beta
    write(iounit) gamma
    write(iounit) eta
    write(iounit) t(1:len)
    write(iounit) lowT(1:len)
    write(iounit) highT(1:len)

    close(iounit)
  end subroutine save_data
  subroutine read_data(alpha, beta, gamma, eta, len, t, lowT, highT)
    use globalFileName
    implicit none
    real(8) :: alpha, beta, gamma, eta
    integer :: len
    real(8), dimension(len) :: t, lowT, highT
    integer, parameter :: iounit = 11
    integer :: ierr

    open(unit=iounit, file=binName, form='unformatted')
    if (ierr /= 0) then
        print *, 'Error opening file:', binName
        return
    end if

    read(iounit) alpha
    read(iounit) beta
    read(iounit) gamma
    read(iounit) eta
    read(iounit) t(1:len)
    read(iounit) lowT(1:len)
    read(iounit) highT(1:len)
  end subroutine read_data

  subroutine write_csv(temp, cond, rad, HT, HG, len, timestep, dx)
    use globalFileName
    implicit none
    real(8), dimension(len) :: temp, cond, rad, HT, HG
    integer :: i, len, timestep
    real(8) :: dx
    character(len=30) :: tmp_timestep

    ! Creating a filename using the timestep
    write(csvName, '(A,I8.8,A)') 'DATA/output_', timestep,  '.csv'
    csvName = trim( adjustl(csvName) )

    ! Opening the file and writing the data
    open(unit=10, file=csvName, status='unknown')
    write(10, *) 'x[m],T[degC],conductivity[-],radiation[-],HeatTransfer[-],HeatGain[-]'
    do i = 1, len
        write(10, '(E16.8, A, E16.8, A, E16.8, A, E16.8, A, E16.8, A, E16.8)')&
        &   dx * (i-1), ',', temp(i)-273.15, ',', cond(i), ',', rad(i), ',',  HT(i), ',', HG(i)
    end do
    close(10)
  end subroutine write_csv

  subroutine correlationFunction(x, y, n, r)
    implicit none
    integer, intent(in)  :: n ! array size
    real(8), intent(in)  :: x(n), y(n)
    real(8), intent(out) :: r ! correlation factor
    real(8) :: sum_x, sum_y, sum_xy, sum_x2, sum_y2
    integer :: i

    sum_x = sum(x)
    sum_y = sum(y)
    sum_xy = sum(x * y)
    sum_x2 = sum(x * x)
    sum_y2 = sum(y * y)

    r = ( n*sum_xy - sum_x*sum_y ) / &
        sqrt( (n*sum_x2 - sum_x*sum_x) * (n*sum_y2 - sum_y*sum_y) )

    if ( isnan(r) ) then
        r = 0.0
    end if
    if (is_inf( abs(r) )) then
        r = 0.0
    end if
  end subroutine correlationFunction

  real(8) function gauss(x, x0, sigma)
    implicit none
    real(8) :: x     ! position
    real(8) :: x0    ! center of the gauss distribution
    real(8) :: sigma ! width of gauss distribution
    gauss = exp( -(x-x0)*(x-x0) * 0.5 / sigma / sigma )
  end function gauss
  logical function is_inf4(val) result(is_inf)
    real(4), intent(in) :: val
    integer(4) :: n
    integer :: i

    n = transfer(val, 0_4)

    ! Exponent part
    is_inf = .true.
    do i = 23, 30
        if( .not. btest(n,i) ) then
            is_inf = .false.
            return
        end if
    end do

    ! Fraction part
    is_inf = .true.
    do i = 0, 22
        if( btest(n,i) ) then
            is_inf = .false.
            return
        end if
    end do
  end function is_inf4
  logical function is_inf8(val) result(is_inf)
    real(8), intent(in) :: val
    integer(8) :: n
    integer :: i

    n = transfer(val, 0_8)

    ! Exponent part
    is_inf = .true.
    do i = 52, 62
        if( .not. btest(n,i) ) then
            is_inf = .false.
            return
        end if
    end do

    ! Fraction part
    is_inf = .true.
    do i = 0, 51
        if( btest(n,i) ) then
            is_inf = .false.
            return
        end if
    end do
  end function is_inf8
end module myProblem
