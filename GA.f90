module ga_module
  implicit none
  public
  !public :: GA_init, GA_step, arithmetic_sequence_sum, arithmetic_sequence_sum_inverse
  real(8), parameter :: mutation_rate = 0.1
  real(8), parameter :: elite_fraction = 0.2
  integer, parameter :: population = 200
  integer, parameter :: indv_size = 4

contains
  subroutine selectRoulette(individuals, eval, selected_invd)
    real(8), dimension(population, indv_size), intent(in) :: individuals
    real(8), dimension(population)           , intent(in) :: eval
    real(8), dimension(indv_size)            , intent(out):: selected_invd
    real(8), dimension(population) :: weight
    real(8) :: w_min
    real(8) :: w_tot
    real(8) :: rnd_num, cum_weight
    integer :: i, j

    ! Create weight array
    w_min = minval( eval )
    if ( w_min < 0 ) then
        do i = 1, population
            weight(i) = eval(i) + w_min
        end do
    else
        do i = 1, population
            weight(i) = eval(i)
        end do
    end if

    ! Generate a random number between 0 and total_weight
    w_tot = sum( weight )
    rnd_num = rnd(1) * w_tot

    ! Select the individual based on the random number and weights
    cum_weight = 0.0
    do i = 1, population
        cum_weight = cum_weight + eval(i)
        if (rnd_num <= cum_weight) then
            do j = 1, indv_size
                selected_invd(j) = individuals(i, j)
            end do
            exit
        end if
    end do
   end subroutine selectRoulette

  ! Initialization of GA
!  subroutine GA_init(individual_max, save_elite, select_method, mutation, individuals, best_individual)
!    integer, intent(in) :: individual_max
!    logical, intent(in), optional :: save_elite
!    character(len=*), intent(in), optional :: select_method
!    real, intent(in), optional :: mutation
!    ! Rest of the code for initialization
!  end subroutine GA_init

  subroutine cross(parent1, parent2, child)
    real(8), intent(in)  :: parent1(indv_size), parent2(indv_size)
    real(8), intent(out) :: child(indv_size)
    integer :: i

    do i = 1, indv_size
        if ( rnd(indv_size) < 0.5 ) then
            child(i) = parent1(i)
        else
            child(i) = parent2(i)
        end if
    end do
  end subroutine cross
  subroutine mutation(mutated_child)
    real(8), dimension(indv_size) :: mutated_child
    integer :: i, rnd_int

    ! alpha
    rnd_int = rnd(1) * 10
    mutated_child(1) = rnd(2) * 0.1 ** rnd_int
    ! beta
    rnd_int = rnd(3) * 10
    mutated_child(2) = rnd(4) * 0.1 ** rnd_int
    ! gamma
    rnd_int = rnd(5) * 10
    mutated_child(3) = rnd(6) * 0.1 ** rnd_int
    ! eta
    rnd_int = rnd(7) * 10
    mutated_child(4) = rnd(8) * 10 ** rnd_int

  end subroutine mutation

  subroutine get_score(cur_indv, eval)
    use problem
    real(8), dimension(indv_size), intent(in) :: cur_indv
    real(8), intent(out) :: eval
    real(8), allocatable :: highT(:), lowT(:), highTRef(:), lowTRef(:)
    integer :: arraySize
    real(8) :: evalL, evalH

    call heat_conduction(cur_indv(1), cur_indv(2), cur_indv(3), cur_indv(4), highT, lowT)
    arraySize = size( highT )
    allocate( highTRef(arraySize) )
    allocate( lowTRef(arraySize) )
    highTRef(:) = 0.0
    lowTRef(:) = 0.0
    call correlationFunction(highTRef, highT, arraySize, evalH)
    call correlationFunction(lowTRef , lowT , arraySize, evalL)

    eval = ( evalH + evalL ) /2.0
  end subroutine get_score

  ! One step of GA
  subroutine GA_step(individuals, eval, best_individual, new_individuals)
    real(8), dimension(population, indv_size), intent(in)  :: individuals, eval
    real(8), dimension(population, indv_size), intent(out) :: new_individuals
    real(8), dimension(indv_size)            , intent(out) :: best_individual
    real(8), dimension(indv_size) :: parent1, parent2, child, mutated_child
    integer :: i, j

    call selectRoulette(individuals, eval, parent1)
    call selectRoulette(individuals, eval, parent2)
    call cross(parent1, parent2, child)

    do i = 1, population
        if ( i <= indv_size * mutation_rate ) then
            do j = 1, indv_size
                new_individuals(i, j) = 0.0
            end do
        else if ( i <= indv_size * (elite_fraction + mutation_rate)&
            & .and. i > indv_size * mutation_rate ) then
            do j = 1, indv_size
                new_individuals(i, j) = best_individual(j)
            end do
        else
            do j = 1, indv_size
                new_individuals(i, j) = child(j)
            end do
        end if
    end do
  end subroutine GA_step

  real(8) function rnd(mixer)
    integer :: mixer, seedsize, i
    integer, allocatable :: seed(:)

    call random_seed(size=seedsize)
    allocate( seed(seedsize) )

    do i = 1, seedsize
        call system_clock( count=seed(i) )
        seed(i) = seed(i) * mixer
    end do
    call random_seed( put=seed(:) )
    call random_number(rnd)
  end function rnd

end module ga_module

