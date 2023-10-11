program heat_conduction
    implicit none
    integer, parameter :: n = 200
    real(8), parameter :: L = 2.45d-2, dx = L / n, dt = 0.0001
    real(8), parameter :: max_time = 80 ! second
    integer, parameter :: max_steps = max_time / dt
    real(8), parameter :: alpha = 8.5194d-6 ! mimicing thermal diffusivity [m2/s]
    real(8), parameter :: beta  = 6.129d-4 ! mimicing heat transfer coefficient [W/m2/K]
    real(8), parameter :: gamma = 1.916d-6 ! mimicing emissivity [-]
    real(8), parameter :: eta   = 643.06  ! mimicing heat generation coefficient [W/m2/K]
    real(8), parameter :: T0 = 300       ! ambient temperature [K]
    real(8), parameter :: sigma  = 5.670374419d-8 ! Stefan-Boltzmann constant [W/m2/K4]
    real(8), parameter :: Wmax  = 100 ! Maximum heater power [W]
    real(8), parameter :: t_Vincr  = 20 ! duration of heater power increase [s]
    real(8), parameter :: t_heat  = 60 ! Heating duration [s]
    integer, parameter :: nOut = 5000 ! Output the results every nOut time step
    real(8) :: ht
    real(8) :: Th_coef ! gradient of Th [K/s]
    real(8) :: Th ! The highest temperature [K]
    real(8) :: T04, u4 ! Temperature squared
    real(8) :: conductivity(n), HeatTrans(n), radiation(n), HeatGain(n)
    real(8) :: u(n), new_u(n)
    integer :: i, t


    ! Setting initial conditions
    Th_coef = Wmax / t_Vincr
    Th = 0.0
    do i = 1, n
        u(i) = T0
        conductivity(i) = 0.0
        HeatTrans(i)    = 0.0
        HeatGain(i)     = 0.0
        radiation(i)    = 0.0
    end do
    print *, alpha * dt / dx /dx

    ! Time stepping loop
    do t = 1, max_steps
        if ( Th < Wmax ) then
            Th = Th_coef * dt * t
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
        new_u(1) = new_u(2)
        if ( dt*t <= t_heat ) then
            HeatGain(n) = eta * Th
            new_u(n) = new_u(n-1) + eta * Th * dx
        else
            HeatGain(n) = 0
            new_u(n) = new_u(n-1)
        end if

        ! Updating the temperature array
        u(:) = new_u(:)

        ! Output temperature distribution every nOut steps to a CSV file
        if ( mod(t, nOut) == 0 ) then
            call write_csv(u, conductivity, radiation, HeatTrans, HeatGain, n, t*dt, dx)
        end if
    end do
contains
    ! Subroutine to write temperature data to a CSV file
    subroutine write_csv(temp, cond, rad, HT, HG, len, timestep, dx)
        real(8), dimension(len) :: temp, cond, rad, HT, HG
        integer :: i, len
        real(8) :: timestep
        real(8) :: dx
        character(len=30) :: tmp_timestep
        character(len=30) :: filename

        ! Creating a filename using the timestep
        !write(tmp_timestep, '(I8)') timestep
        write(filename, '(A,e13.6,A)') 'DATA/output_', timestep, '.csv'
        filename = trim( adjustl(filename) )

        ! Opening the file and writing the data
        open(unit=10, file=filename, status='unknown')
        write(10, *) 'x[m],T[degC],conductivity[-],radiation[-],HeatTransfer[-],HeatGain[-]'
        do i = 1, len
            write(10, '(E16.8, A, E16.8, A, E16.8, A, E16.8, A, E16.8, A, E16.8)')&
            &   dx * (i-1), ',', temp(i)-273.15, ',', cond(i), ',', rad(i), ',',  HT(i), ',', HG(i)
        end do
        close(10)
    end subroutine write_csv
    real(8) function heat_trans_coef(dT, C)
        real(8), parameter :: L = 4.5d-3 ! heater diameter [m]
        real(8) :: dT
        real(8) :: C

        heat_trans_coef = 2.51 * C * (dT/L)**(0.25)
    end function heat_trans_coef
    real(8) function gauss(x, x0, sigma)
        real(8) :: x     ! position
        real(8) :: x0    ! center of the gauss distribution
        real(8) :: sigma ! width of gauss distribution
        gauss = exp( -(x-x0)*(x-x0) * 0.5 / sigma / sigma )
    end function gauss
end program heat_conduction
