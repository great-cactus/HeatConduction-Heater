program heat_conduction
    implicit none
    integer, parameter :: n = 200, max_steps = 10000
    real, parameter :: L = 0.1, dx = L / n, dt = 0.0001
    real, parameter :: alpha = 2.5e-5 ! thermal diffusivity [m2/s]
    real, parameter :: h = 1 ! heat transfer coefficient [W/m2/K]
    real, parameter :: T0 = 300 ! ambient temperature [K]
    real, parameter :: t_heat  = 0.05 ! Heating duration [s]
    real, parameter :: Temp_heat  = 1000 ! The highest heating temperature [K]
    real :: Th_coef ! gradient of Th [K/s]
    real :: Th ! The highest temperature [K]
    real :: gauss
    real :: conductivity(n), HeatTrans(n), HeatGain(n)
    real :: u(n), new_u(n), x
    integer :: i, t

    ! Setting initial conditions
    Th_coef = Temp_heat / t_heat
    Th = T0
    do i = 1, n
        x = (i - 1) * dx
        gauss = exp( -(x - L/2.) / (L/10.)*(x - L/2.) / (L/10.) )
        u(i) = T0 + ( Th-T0  )* gauss
        conductivity(i) = 0
        HeatTrans(i) = 0
        HeatGain(i) = 0
    end do

    ! Time stepping loop
    do t = 1, max_steps
        Th = Th_coef * dt * i
        ! Solving the heat conduction differential equation
        do i = 2, n-1
            conductivity(i) = alpha / dx**2 * ( u(i+1) - 2*u(i) + u(i-1) )
            HeatTrans(i) = h * ( u(i) - T0 )
            gauss = exp( -((i-1)*dx - L/2.) / (L/10.) * ((i-1)*dx - L/2.)/ (L/10.) )
            if ( dt * t <= t_heat ) then
                write(*,*)i, dt*i, t_heat
                HeatGain(i) = 5e3 * Th * gauss
            else
                HeatGain(i) = 0
            end if
            new_u(i) = u(i) + dt * ( conductivity(i) - HeatTrans(i) + HeatGain(i) )
        end do
        new_u(1) = T0
        new_u(n) = ( h * T0 * dx + alpha * u(n-1) ) / ( h*dx + alpha )

        ! Updating the temperature array
        u(:) = new_u(:)

        ! Output temperature distribution every 10 steps to a CSV file
        if ( mod(t, 100) == 0 ) then
            call write_csv(u, conductivity, HeatTrans, HeatGain, n, t, dx)
        end if
    end do
contains
    ! Subroutine to write temperature data to a CSV file
    subroutine write_csv(data, cond, HT, HG, len, timestep, dx)
        real, dimension(len) :: data, cond, HT, HG
        integer :: i, len, timestep
        real :: dx
        character(len=20) :: filename

        ! Creating a filename using the timestep
        write(filename, '(A,I4.4,A)') 'DATA/output_', timestep, '.csv'

        ! Opening the file and writing the data
        open(unit=10, file=filename, status='unknown')
        write(10, *) 'x[cm],T[K],conductivity[-],HeatTransfer[-],HeatGain[-]'
        do i = 1, len
            write(10, '(E16.8, A, E16.8, A, E16.8, A, E16.8, A, E16.8 )') dx * (i-1), ',', data(i), ',', cond(i), ',', HT(i), ',', HG(i)
        end do
        close(10)
    end subroutine write_csv
end program heat_conduction
