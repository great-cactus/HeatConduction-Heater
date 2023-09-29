subroutine correlationFunction(x, y, n, r)
    implicit none
    integer, intent(in) :: n ! array size
    real(8), intent(in), dimension :: x(n), y(n)
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
