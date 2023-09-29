program main
    implicit none

contains
    subroutine correlationFunction(x, y, n, r)
        integer, intent(in)  :: n ! array size
        real(8), intent(in)  :: x(n), y(n)
        real(8), intent(out) :: r ! correlation factor
        real(8) :: sum_x, sum_y, sum_xy, sum_x2, sum_y2
        real(8) :: numerator, denominator
        integer :: i

        sum_x = 0.0
        sum_y = 0.0
        sum_xy = 0.0
        sum_x2 = 0.0
        sum_y2 = 0.0

        do i = 1, n
            sum_x = sum_x + x(i)
            sum_y = sum_y + y(i)
            sum_xy = sum_xy + x(i) * y(i)
            sum_x2 = sum_x2 + x(i) * x(i)
            sum_y2 = sum_y2 + y(i) * y(i)
        end do

        numerator = n*sum_xy - sum_x*sum_y
        denominator = sqrt((n*sum_x2 - sum_x*sum_x)*(n*sum_y2 - sum_y*sum_y))
        r = numerator / denominator
    end subroutine correlationFunction

    subroutine heat_conduction(alpha, beta, gamma, eta, n, highT, lowT)
        real(8), intent(in)  :: alpha, beta, gamma, eta ! fitting parameters
        integer, intelt(in)  :: n ! array size
        real(8), intent(out) :: highT(n), lowT(n) ! temperature array for fitting
        integer, parameter :: max_steps = 800000
        real(8), parameter :: L = 2.45d-2, dx = L / n, dt = 0.0001
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
            new_u(1) = T0
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
                call write_csv(u, conductivity, radiation, HeatTrans, HeatGain, n, t, dx, dir)
            end if
        end do
    end subroutine heat_conduction

    subroutine save_data(temp, cond, rad, HT, HG, len, timestep, dx, dir)
        real(8), dimension(len) :: temp, cond, rad, HT, HG
        character(*) :: dir
        integer :: i, len, timestep
        real(8) :: dx
        character(len=30) :: tmp_timestep
        character(len=30) :: filename
        integer :: iounit
        integer :: ierr

        ! Creating a filename using the timestep
        write(filename, '(A, A,I8.8,A)') dir, '/output_', timestep,  '.bin'
        filename = trim( adjustl(filename) )

        open(newunit=iounit, file=filename, form='unformatted')
        if (ierr /= 0) then
            print *, 'Error opening file:', filename
            return
        end if
        write(iounit) cond, rad, HT, HG
    end subroutine save_data

    subroutine write_csv(temp, cond, rad, HT, HG, len, timestep, dx, dir)
        real(8), dimension(len) :: temp, cond, rad, HT, HG
        character(*) :: dir
        integer :: i, len, timestep
        real(8) :: dx
        character(len=30) :: tmp_timestep
        character(len=30) :: filename

        ! Creating a filename using the timestep
        write(filename, '(A, A,I8.8,A)') dir, '/output_', timestep,  '.csv'
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

    real(8) function gauss(x, x0, sigma)
        real(8) :: x     ! position
        real(8) :: x0    ! center of the gauss distribution
        real(8) :: sigma ! width of gauss distribution
        gauss = exp( -(x-x0)*(x-x0) * 0.5 / sigma / sigma )
    end function gauss
end program main
