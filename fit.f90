program heat_conduction
    implicit none
    integer, parameter :: n = 200, max_steps = 600000
    real(8), parameter :: L = 2.85d-2, dx = L / n, dt = 0.0001
    real(8), parameter :: x_heatMax = 2.45d-2 ! Maximum heating position [m]
    real(8), parameter :: alpha = 2.5d-5 ! mimicing thermal diffusivity [m2/s]
    real(8), parameter :: beta  = 0.52 ! mimicing heat transfer coefficient [W/m2/K]
    real(8), parameter :: gamma = 1.0 ! mimicing emissivity [-]
    real(8), parameter :: eta   = 1.0 ! mimicing heat generation coefficient [W/m2/K]
    real(8), parameter :: T0 = 300       ! ambient temperature [K]
    real(8), parameter :: sigma  = 5.670374419d-8 ! Stefan-Boltzmann constant [W/m2/K4]
    real(8), parameter :: Wmax  = 100 ! Maximum heater power [W]
    real(8), parameter :: t_Vincr  = 10 ! duration of heater power increase [s]
    real(8), parameter :: t_heat  = 30 ! Heating duration [s]
    real(8) :: ht
    real(8) :: Th_coef ! gradient of Th [K/s]
    real(8) :: Th ! The highest temperature [K]
    real(8) :: T04, u4 ! Temperature squared
    real(8) :: conductivity(n), HeatTrans(n), HeatGain(n), radiation(n)
    real(8) :: u(n), new_u(n), x
    integer :: i, t

    ! Setting initial conditions
    Th_coef = Wmax / t_Vincr
    Th = 0
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
            Th = Th_coef * dt * t
        else
            Th = Wmax
        end if
        ! Solving the heat conduction differential equation
        do i = 2, n-1
            ht = heat_trans_coef(u(i)-T0, beta)
            conductivity(i) = alpha / dx**2 * ( u(i+1) - 2*u(i) + u(i-1) )
            HeatTrans(i) = ht * ( u(i) - T0 )
            u4  = u(i) * u(i) * u(i) * u(i)
            T04 = T0 * T0 * T0 * T0
            radiation(i) = sigma * gamma * ( u4-T04 )
            if ( dt*t <= t_heat ) then
                HeatGain(i) = eta * Th * gauss( (i-1)*dx, x_heatMax, L/20.0 )
            else
                HeatGain(i) = 0
            end if
            new_u(i) = u(i) + dt * ( conductivity(i) - radiation(i) - HeatTrans(i) + HeatGain(i) )
        end do
        new_u(1) = T0
        new_u(n) = ( ht * T0 * dx + alpha * u(n-1) ) / ( ht*dx + alpha )

        ! Updating the temperature array
        u(:) = new_u(:)

        ! Output temperature distribution every 10 steps to a CSV file
        if ( mod(t, 10000) == 0 ) then
            call write_csv(u, conductivity, radiation, HeatTrans, HeatGain, n, t, dx)
        end if
    end do
contains
    ! Subroutine to write temperature data to a CSV file
    subroutine write_csv(temp, cond, rad, HT, HG, len, timestep, dx)
        real(8), dimension(len) :: temp, cond, rad, HT, HG
        integer :: i, len, timestep
        real(8) :: dx
        character(len=30) :: tmp_timestep
        character(len=30) :: filename

        ! Creating a filename using the timestep
        !write(tmp_timestep, '(I8)') timestep
        write(filename, '(A,I8.8,A)') 'DATA/output_', timestep,  '.csv'
        filename = trim( adjustl(filename) )

        ! Opening the file and writing the data
        open(unit=10, file=filename, status='unknown')
        write(10, *) 'x[m],T[K],conductivity[-],radiation[-],HeatTransfer[-],HeatGain[-]'
        do i = 1, len
            write(10, '(E16.8, A, E16.8, A, E16.8, A, E16.8, A, E16.8, A, E16.8)')&
            &   dx * (i-1), ',', temp(i), ',', cond(i), ',', rad(i), ',',  HT(i), ',', HG(i)
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
