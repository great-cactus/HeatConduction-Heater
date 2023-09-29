module ga_module
  implicit none
  private
  public :: GA_init, GA_step, arithmetic_sequence_sum, arithmetic_sequence_sum_inverse
  real(8), parameter :: mutation_rate = 0.1
  real(8), parameter :: elite_fraction = 0.2
  integer, parameter :: population = 200
  integer, parameter :: indv_size = 5

contains

  ! Arithmetic sequence sum
  real function arithmetic_sequence_sum(size, start, diff)
    integer, intent(in) :: size
    real, intent(in), optional :: start, diff
    real :: start_, diff_
    start_ = 1
    diff_ = 1
    if (present(start)) start_ = start
    if (present(diff)) diff_ = diff
    arithmetic_sequence_sum = size * (2 * start_ + (size - 1) * diff_) / 2
  end function arithmetic_sequence_sum

  ! Arithmetic sequence sum inverse
  real function arithmetic_sequence_sum_inverse(val, start, diff)
    real, intent(in) :: val
    real, intent(in), optional :: start, diff
    real :: start_, diff_, t
    start_ = 1
    diff_ = 1
    if (present(start)) start_ = start
    if (present(diff)) diff_ = diff
    if (diff_ == 0) then
      arithmetic_sequence_sum_inverse = val
      return
    end if
    t = diff_ - 2 * start_ + sqrt((2 * start_ - diff_) ** 2 + 8 * diff_ * val)
    arithmetic_sequence_sum_inverse = t / (2 * diff_)
  end function arithmetic_sequence_sum_inverse

  subroutine selectRoulette(individuals)
      integer :: i
  end subroutine selectRoulette

  ! Initialization of GA
  subroutine GA_init(individual_max, save_elite, select_method, mutation, individuals, best_individual)
    integer, intent(in) :: individual_max
    logical, intent(in), optional :: save_elite
    character(len=*), intent(in), optional :: select_method
    real, intent(in), optional :: mutation
    ! Rest of the code for initialization
  end subroutine GA_init

  subroutine cross(parent1, parent2, child)
    real(8), intent(in)  :: parent1(indv_size), parent2(indv_size)
    real(8), intent(out) :: child(indv_size)
    integer :: i

    do i = 1, indv_size
        if rnd() < 0.5 then
            child(i) = parent1(i)
        else
            child(i) = parent2(i)
        end if
    end do
  end subroutine cross
  subroutine mutation(mutated_child)
      real(8) :: mutated_child(indv_size)

    do i = 1, indv_size
        mutated_child(i) = rnd()
    end do

  end subroutine mutation

  ! One step of GA
  subroutine GA_step(individuals, best_individual, new_individuals)
    real(8), intent(in)  :: individuals(population)
    real(8), intent(out) :: best_individual
    real(8), intent(out) :: new_individuals(population)
    real(8) :: parent1, parent2, child, mutated_child
    integer :: i

    parent1 = selectRoulette(individuals)
    parent2 = selectRoulette(individuals)
    child = cross(individuals)

    do i = 1, population
        if i <= indv_size * mutation_rate then
            new_individuals(i) = 0
        else if ( i <= indv_size * (elite_fraction + mutation_rate)&
            & .and. i > indv_size * mutation_rate ) then
            new_individuals(i) = best_individual
        else
            new_individuals(i) = child
        end if
    end do
  end subroutine GA_step

end module ga_module

