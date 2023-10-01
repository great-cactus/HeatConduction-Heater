program main
  use globalFileName
  use ga_module
  implicit none
  integer :: i

  call runGA()
  call write_csv(evals, population, eval_func_len, "evalL,evalH,evalTot", "eval.csv")
contains
  subroutine write_csv(array, numRows, numCols, colName, fileName)
    integer     , intent(in) :: numRows, numCols
    real(8)     , intent(in) :: array(numRows, numCols)
    character(*), intent(in) :: colName, fileName
    integer, parameter :: ioUnit = 15
    integer :: i, j
    character(20) :: fmt

    write(fmt, '("(F0.",I0,"A)")') 9, ','

    open(unit=ioUnit, file=fileName, action='write', status='replace')
    write(ioUnit, *) colName

    do i = 1, numRows
        do j = 1, numCols
            if ( j > numCols ) then
                write(ioUnit, fmt, advance='no') array(i, j)
            else
                write(ioUnit, '(F0.9)') array(i, j)
            end if
        end do
    end do

    close(ioUnit)
  end subroutine write_csv
end program main
