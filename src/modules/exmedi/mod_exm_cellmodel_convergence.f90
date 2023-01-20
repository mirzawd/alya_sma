!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_exm_cellmodel_convergence
    use def_master, only : rp, ip, lg
    implicit none

    private

    integer(ip), parameter, public :: &
        METRIC_RMSE = 1_ip


    type ConvergenceChecker
        private 
        real(rp), dimension(:), pointer, public :: array             ! fill this array to check for convergence
        real(rp), dimension(:), pointer         :: array_prev        ! array at  previous iteration to compare
        real(rp)                                :: toler             ! tolerance
        integer(ip)                             :: steps_max         ! number of consecutive steps the difference has to be < toler
        integer(ip)                             :: steps             ! counter from 1 to steps_max
        integer(4)                              :: metric            ! comparison metric 
        real(rp)                                :: metric_value      ! value

        contains 
        private
        procedure, public :: Init             ! Initialize all variables
        procedure, public :: SetMaxSteps      ! set steps_max
        procedure, public :: GetMaxSteps      ! Get steps_max
        procedure, public :: SetTolerance     ! change tolerance
        procedure, public :: GetTolerance     ! get tolerance
        procedure, public :: EndIteration     ! save the current array as previous
        procedure, public :: Converged        ! return if the difference is < tolerance and steps==steps_max
        procedure, public :: GetMetricValue   ! get latest result of metric evaluation

        final :: destructor
    end type ConvergenceChecker

    public :: ConvergenceChecker

contains

real(rp) function GetMetricValue(this)
    implicit none
    class(ConvergenceChecker) :: this

    GetMetricValue = this % metric_value
end function GetMetricValue



subroutine EndIteration(this)
    use mod_memory, only : memory_size
    implicit none
    class(ConvergenceChecker) :: this
    integer(ip) :: i

    do i = 1,memory_size(this % array_prev, 1_ip)
        this % array_prev(i) = this % array(i)
    end do

end subroutine EndIteration



logical(lg) function Converged(this)
    use mod_memory, only : memory_size
    implicit none
    class(ConvergenceChecker) :: this
    integer(ip) :: n

    select case ( this % metric )
    case ( METRIC_RMSE )
        n = memory_size(this % array, 1_ip)

        this % metric_value = sqrt(sum( (this % array - this % array_prev)**2/n ))
        if ( this % metric_value < this % toler ) then
            this % steps = this % steps + 1
        else
            this % steps = 0_ip
        end if
    case default
        call runend("ConvergenceChecker: Unknown metric.")
    end select

    Converged = (this % steps) .GE. (this % steps_max)
end function Converged


real(rp) function GetTolerance(this)
    implicit none
    class(ConvergenceChecker) :: this

    GetTolerance = this % toler
end function GetTolerance


subroutine SetTolerance(this, toler)
    implicit none
    class(ConvergenceChecker) :: this
    real(rp), intent(in) :: toler

    this % toler = toler
end subroutine SetTolerance




integer(ip) function GetMaxSteps(this)
    implicit none
    class(ConvergenceChecker) :: this

    GetMaxSteps = this % steps_max
end function GetMaxSteps


subroutine SetMaxSteps(this, max_steps)
    implicit none
    class(ConvergenceChecker) :: this
    integer(ip), intent(in) :: max_steps

    this % steps_max = max_steps
end subroutine SetMaxSteps


subroutine destructor(this)
    use mod_memory, only : memory_deallo
    use def_master, only : mem_modul, modul
    implicit none
    type(ConvergenceChecker) :: this

    call memory_deallo(mem_modul(1:2, modul), 'ConvergenceChecker % array', 'mod_exm_cellmodel_convergence', this % array)
    call memory_deallo(mem_modul(1:2, modul), 'ConvergenceChecker % array_prev', 'mod_exm_cellmodel_convergence', this % array_prev)
end subroutine destructor


subroutine Init(this, toler, array_length)
    use mod_memory, only : memory_alloca
    use def_master, only : mem_modul, modul

    implicit none
    class(ConvergenceChecker) this
    real(rp)    :: toler
    integer(ip) :: array_length
    integer(ip) :: i

    nullify(this % array)
    nullify(this % array_prev)

    this % toler        = toler    
    this % steps        = 0_ip
    this % steps_max    = 3_ip         
    this % metric       = METRIC_RMSE
    this % metric_value = 0.0_rp
    
    call memory_alloca(mem_modul(1:2, modul), 'ConvergenceChecker % array', 'mod_exm_cellmodel_convergence', this % array, array_length)
    call memory_alloca(mem_modul(1:2, modul), 'ConvergenceChecker % array_prev', 'mod_exm_cellmodel_convergence', this % array_prev, array_length)

    do i = 1,array_length
        this % array_prev(i) = 1.0e+30_rp ! some large numbers
        this % array(i) = 0.0_rp
    end do


end subroutine Init



end module mod_exm_cellmodel_convergence
