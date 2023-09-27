program heat_conduction
    implicit none
    integer, parameter :: n = 200, max_steps = 10000
    real(8), parameter :: L = 0.1, dx = L / n, dt = 0.0001
    real(8), parameter :: alpha = 2.5d-5 ! thermal diffusivity [m2/s]
    real(8), parameter :: ht = 40 ! heat transfer coefficient [W/m2/K]
    real(8), parameter :: hg = 60 ! heat gain coefficient [W/m2]
    real(8), parameter :: T0 = 300 ! ambient temperature [K]
    real(8), parameter :: t_heat  = 0.05 ! Heating duration [s]
    real(8), parameter :: Temp_heat  = 1000 ! The highest heating temperature [K]
    real(8), parameter :: sigma  = 5.670374419d-8 ! Stefan-Boltzmann constant [W/m2/K4]
    real(8), parameter :: epsilon  = 1 ! emissivity
    real(8) :: Th_coef ! gradient of Th [K/s]
    real(8) :: Th ! The highest temperature [K]
    real(8) :: T04, u4 ! Temperature squared
    real(8) :: gauss
    real(8) :: conductivity(n), HeatTrans(n), HeatGain(n), radiation(n)
    real(8) :: u(n), new_u(n), x
    integer :: i, t

    ! Setting initial conditions
    Th_coef = Temp_heat / t_heat
    Th = T0
    do i = 1, n
        x = (i - 1) * dx
        gauss = exp( -(x - L/2.) / (L/10.)*(x - L/2.) / (L/10.) )
        u(i) = T0 + ( Th-T0  )* gauss
        conductivity(i) = 0
        HeatTrans(i)    = 0
        HeatGain(i)     = 0
        radiation(i)    = 0
    end do

    ! Time stepping loop
    do t = 1, max_steps
        Th = Th_coef * dt * i
        ! Solving the heat conduction differential equation
        do i = 2, n-1
            conductivity(i) = alpha / dx**2 * ( u(i+1) - 2*u(i) + u(i-1) )
            HeatTrans(i) = ht * ( u(i) - T0 )
            u4  = u(i) * u(i) * u(i) * u(i)
            T04 = T0 * T0 * T0 * T0
            radiation(i) = sigma * epsilon * ( u4-T04 )
            gauss = exp( -((i-1)*dx - L/2.) / (L/10.) * ((i-1)*dx - L/2.)/ (L/10.) )
            if ( dt*t <= t_heat ) then
                HeatGain(i) = hg * Th * gauss
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
        if ( mod(t, 10) == 0 ) then
            call write_csv(u, conductivity, radiation, HeatTrans, HeatGain, n, t, dx)
        end if
    end do
contains
    ! Subroutine to write temperature data to a CSV file
    subroutine write_csv(temp, cond, rad, HT, HG, len, timestep, dx)
        real(8), dimension(len) :: temp, cond, rad, HT, HG
        integer :: i, len, timestep
        real(8) :: dx
        character(len=20) :: filename

        ! Creating a filename using the timestep
        write(filename, '(A,I4.4,A)') 'DATA/output_', timestep, '.csv'

        ! Opening the file and writing the data
        open(unit=10, file=filename, status='unknown')
        write(10, *) 'x[cm],T[K],conductivity[-],radiation[-],HeatTransfer[-],HeatGain[-]'
        do i = 1, len
            write(10, '(E16.8, A, E16.8, A, E16.8, A, E16.8, A, E16.8, A, E16.8)')&
            &   dx * (i-1), ',', temp(i), ',', cond(i), ',', rad(i), ',',  HT(i), ',', HG(i)
        end do
        close(10)
    end subroutine write_csv
end program heat_conduction
