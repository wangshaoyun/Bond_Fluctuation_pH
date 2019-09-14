program main
use global_variables
use initialize_update
use compute_energy
use input_output
implicit none

  !##########Data Dictionary############!
    integer :: i, j, k
    real*8  :: EE=0, EE1=0, st, fn
    real*8  :: DeltaE, time(3)
  !#####################################!

  !################begin#################!
  call cpu_time(started)
  call random_seed()
  !
  !input and initialize system, timing and histogram parameters.
  call initialize_parameters
!   !
!   !initialize energy and parameters of potential
!   call initialize_energy_parameters
  !
  !
!   if (restart_or_continue == 0) then
    i=1
    !
    !initialize position  
    call initialize_position
    call write_pos
    call write_pos1(i)
    !
    !error analysis
    call error_analysis
    call write_hist    
  else if (restart_or_continue /= 0) then
    !
    !read position and histogram data
    call continue_read_data(i)
    !
    !error analysis
    call error_analysis
  end if
  !##############preheating##############!
  if ( i <= StepNum0 ) then
    do step=i, StepNum0
      call monte_carlo_move(EE, DeltaE)
      if ( mod(step,DeltaStep1) == 0 ) then
        call compute_physical_quantities
        call write_physical_quantities( step )
      end if
      if ( mod(step,DeltaStep3) == 0 ) then
        call error_analysis
        call write_pos1(step)
      end if
    end do
    i=step
  end if

  !################running###############!
  do step=i,StepNum0+StepNum            
    call monte_carlo_move(EE, DeltaE)
    if ( mod(step,DeltaStep1) == 0 ) then
      call compute_physical_quantities
      call write_physical_quantities( step )
    end if
    if ( mod(step,DeltaStep2) == 0 ) then     
      call histogram
    end if
    if ( mod(step,DeltaStep3) == 0 ) then
      call write_pos1(step)               
      call write_hist
    end if
  end do

  !##################end#################!
  call cpu_time(finished)
  total_time=finished-started+total_time
  call write_time(total_time)
  write(*,*) 'Finished!'

end program main



