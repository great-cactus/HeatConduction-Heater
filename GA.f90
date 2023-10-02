module ga_module
  implicit none
  public
  !....
  !    Here is the user definition parameters
  !....
  real(8), parameter :: mutation_rate = 0.5  ! mutation rate
  real(8), parameter :: elite_fraction = 0.1 ! elite fraction
  integer, parameter :: indv_size = 4        ! size of each gene
  integer, parameter :: population = 10      ! population
  integer, parameter :: max_gen = 5          ! maximum generation
  integer, parameter :: eval_func_len = 3    ! length of evaluation function array
  !....
  !....
  real(8), dimension(population, eval_func_len) :: evals

contains
  subroutine runGA()
    integer :: i
    real(8), dimension(max_gen, indv_size)    :: saved_elites
    real(8), dimension(population, indv_size) :: individuals, new_individuals
    real(8), dimension(indv_size) :: best_individual, king_of_kings
    real(8), dimension(eval_func_len) :: cur_best_eval
    real(8) :: eval_best

    eval_best = 0.0
    call GA_init(individuals, saved_elites)
    do i = 1, max_gen
        print *, i, 'generation'
        call GA_step(individuals, best_individual, new_individuals, cur_best_eval)
        saved_elites(i, :) = best_individual(:)
        evals(i, :) = cur_best_eval(:)
        individuals(:, :) = new_individuals(:, :)

        if ( eval_best < cur_best_eval(eval_func_len) ) then
            eval_best = cur_best_eval(eval_func_len)
            king_of_kings(:) = best_individual(:)
        end if
    end do

    print *, 'Finished running'
    print *, 'The king of kings is'
    print *, king_of_kings(:)
    do i = 1, max_gen
        print *, evals(i, :)
    end do

  end subroutine runGA

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
    rnd_num = rnd() * w_tot

    ! Select the individual based on the random number and weights
    cum_weight = 0.0
    do i = 1, population
        cum_weight = cum_weight + eval(i)
        if (rnd_num <= cum_weight) then
            selected_invd(:) = individuals(i, :)
            exit
        end if
    end do
   end subroutine selectRoulette

  ! Initialization of GA
  subroutine GA_init(individuals, saved_elites)
    real(8), dimension(max_gen, indv_size)   , intent(out) :: saved_elites
    real(8), dimension(population, indv_size), intent(out) :: individuals
    real(8) :: rnd_num, diff_ratio
    integer :: i, j

    individuals(1, 1) = 1.0d-5
    individuals(1, 2) = 1.0d-5
    individuals(1, 3) = 1.0d-4
    individuals(1, 4) = 5.0d2
    do i = 1, population - 10
        rnd_num = rnd()
        diff_ratio = rnd_num / ( 1.0 - rnd_num )
        individuals(i, 1) = individuals(1, 1) * diff_ratio
        rnd_num = rnd()
        diff_ratio = rnd_num / ( 1.0 - rnd_num )
        individuals(i, 2) = individuals(1, 2) * diff_ratio
        rnd_num = rnd()
        diff_ratio = rnd_num / ( 1.0 - rnd_num )
        individuals(i, 3) = individuals(1, 3) * diff_ratio
        rnd_num = rnd()
        diff_ratio = rnd_num / ( 1.0 - rnd_num )
        individuals(i, 4) = individuals(1, 4) * diff_ratio
    end do
    do i = ( population - 10), population
        call mutation(individuals(i, :))
    end do

    do i = 1, max_gen
        do j = 1, indv_size
            saved_elites(i, j) = 0.0
        end do
    end do
  end subroutine GA_init

  subroutine cross(parent1, parent2, child)
    real(8), intent(in)  :: parent1(indv_size), parent2(indv_size)
    real(8), intent(out) :: child(indv_size)
    integer :: i

    do i = 1, indv_size
        if ( rnd() < 0.5 ) then
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
    rnd_int = rnd() * 10
    mutated_child(1) = rnd() * 0.1 ** rnd_int
    ! beta
    rnd_int = rnd() * 10
    mutated_child(2) = rnd() * 0.1 ** rnd_int
    ! gamma
    rnd_int = rnd() * 10
    mutated_child(3) = rnd() * 0.1 ** rnd_int
    ! eta
    rnd_int = rnd() * 10
    mutated_child(4) = rnd() * 10 ** rnd_int

  end subroutine mutation

  subroutine get_score(cur_indv, eval)
    use problem
    real(8), dimension(indv_size)    , intent(in)  :: cur_indv
    real(8), dimension(eval_func_len), intent(out) :: eval
    real(8), allocatable :: highT(:), lowT(:)
    real(8), dimension(80) :: highTRef = (/39.9,106.6,174.0,236.2 ,293.8 ,348.2 ,399.4 ,&
        &447.7 ,494.7 ,540.4 ,585.5 ,629.8 ,674.1 ,714.3 ,751.4 ,786.7 ,821.9 ,852.5 ,878.3 ,897.6,&
        &916.7 ,936.5 ,955.9 ,974.8 ,989.9 ,1004.3,1014.0,1023.8,1032.4,1038.5,1042.5,1047.6,1051.4,&
        &1056.0,1055.3,1056.1,1057.3,1058.6,1061.0,1060.5,1060.7,1061.1,1061.9,1062.4,1063.6,&
        &1063.3,1062.7,1063.4,1063.4,1062.8,1062.1,1061.4,1056.2,1062.7,1066.7,1065.8,1063.7,&
        &1063.7,1063.5,1065.4,1064.4,1024.1,971.8 ,922.7 ,877.5 ,837.5 ,799.7 ,765.8 ,737.5 ,710.4,&
        &684.1 ,658.6 ,634.0 ,611.6 ,590.2 ,570.4 ,551.3 ,532.6 ,516.1 ,499.8/) !,484.3/)
    real(8), dimension(80) :: lowTRef = (/27.5,35.8,51.3,71.0,93.1,116.8,141.3,166.0,190.3,&
        &214.5,238.6,262.3,285.1,307.9,329.9,351.5,372.4,392.6,412.2,430.6,448.3,465.3,&
        &481.5,497.0,511.6,524.6,536.7,547.9,558.0,567.0,575.0,582.2,588.6,594.1,599.0,603.9,608.0,&
        &611.9,615.3,618.6,621.7,623.9,624.6,626.6,628.7,630.9,633.0,635.2,637.3,639.0,640.6,&
        &642.1,643.1,644.3,644.6,645.7,647.7,649.2,650.4,651.4,652.5,647.3,634.4,617.6,&
        &599.0,579.0,559.2,539.8,521.3,503.6,486.9,471.1,456.3,442.2,428.9,416.3,404.2,392.7,&
        &381.7,371.3/)!,361.6/)
    integer :: arraySize, max_idx
    real(8) :: evalL, evalH
    real(8) :: highTMaxRef, lowTMaxRef
    real(8) :: highTMax, lowTMax

    call heat_conduction(cur_indv(1), cur_indv(2), cur_indv(3), cur_indv(4), highT, lowT)
    arraySize = size( highT )
    call correlationFunction(highTRef, highT, arraySize, evalH)
    call correlationFunction(lowTRef , lowT , arraySize, evalL)
    deallocate( highT )
    deallocate( lowT )

    !call find_maxIdx(highTRef, arraySize, max_idx)
    !highTMaxRef = highTRef(max_idx)
    !call find_maxIdx(lowTRef, arraySize, max_idx)
    !lowTMaxRef = lowTRef(max_idx)

    eval(1) = evalL
    eval(2) = evalH
    eval(3) = ( evalH + evalL ) /2.0
  end subroutine get_score

  ! One step of GA
  subroutine GA_step(individuals, best_individual, new_individuals, eval_best)
    real(8), dimension(population, indv_size), intent(in)  :: individuals
    real(8), dimension(population, indv_size), intent(out) :: new_individuals
    real(8), dimension(indv_size)            , intent(out) :: best_individual
    real(8), dimension(eval_func_len)        , intent(out) :: eval_best
    real(8), dimension(population, eval_func_len):: eval
    real(8), dimension(indv_size) :: parent1, parent2, child, mutated_child, tmp_indv
    integer :: i, j, best_idx

    ! Get score
    do i = 1, population
        print *, 'population', i
        call get_score(individuals(i,:), eval(i, :))
    end do
    ! Find best individual
    call find_maxIdx(eval(:, eval_func_len), population, best_idx)
    best_individual(:) = individuals(best_idx, :)
    eval_best(:) = eval(best_idx, :)

    ! Create crossed individual
    call selectRoulette(individuals, eval(:,eval_func_len), parent1)
    call selectRoulette(individuals, eval(:,eval_func_len), parent2)
    call cross(parent1, parent2, child)

    ! Next generation
    do i = 1, population
        if ( i <= population * mutation_rate ) then
            call mutation( new_individuals(i, :))
        else if ( i <= population * (elite_fraction + mutation_rate)&
            & .and. i > population * mutation_rate ) then
            new_individuals(i, :) = best_individual(:)
        else
            new_individuals(i, :) = child(:)
        end if
    end do
  end subroutine GA_step

  real(8) function rnd()
    integer :: seedsize
    integer, allocatable :: seed(:)

    ! Setup the seed
    call random_seed( size=seedsize )
    allocate( seed(seedsize) )
    call random_seed( get=seed(:) )
    call random_seed( put=seed(:) )
    deallocate( seed )

    call random_number(rnd)
  end function rnd

  subroutine find_maxIdx(array, len, max_idx)
    real(8), dimension(len), intent(in) :: array
    integer, intent(in) :: len
    integer, intent(out) :: max_idx
    real(8) :: max_val
    integer :: i

    max_val = array(1)
    max_idx = 1

    do i = 2, len
        if ( array(i) > max_idx ) then
            max_val = array(i)
            max_idx = i
        end if
    end do
  end subroutine find_maxIdx

end module ga_module
