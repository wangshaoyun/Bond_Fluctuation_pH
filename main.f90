program main
use global_variables
use initialize_update
use compute_energy
use input_output
implicit none

  !#################data#################!
  integer :: i
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
    call write_pos1
!     !
!     !Construct the real verlet list and real_point vector
!     call construct_cell_list
!     !
!     !compute energy
!     call compute_energy
!     call write_hist    
!   else if (restart_or_continue /= 0) then
!     !
!     !read position and histogram data
!     call continue_read_data(i)
!     !
!     !Construct the real verlet list and real_point vector
!     call real_verlet_list
!     !
!     !compute energy
!     call compute_energy
!   end if
  !##############preheating##############!
  if ( i <= StepNum0 ) then
    do step=i, StepNum0
      call monte_carlo_move
!       call update_cell_list
!       if ( mod(step,DeltaStep1) == 0 ) then
!         call height
!         call write_height(step)
!       end if
!       if ( mod(step,DeltaStep3) == 0 ) then
!         call write_pos1
!       end if
    end do
    i=step
  end if

!   !################running###############!
!   do step=i,StepNum0+StepNum            
!     call new_position
!     call update_verlet_list
!     if ( mod(step,DeltaStep1) == 0 ) then
!       call height
!       call write_height(step)
!     end if
!     if ( mod(step,DeltaStep2) == 0 ) then     
!       call histogram
!     end if
!     if ( mod(step,DeltaStep3) == 0 ) then
!       call write_pos1                   
!       call write_hist
!     end if
!   end do

!   !##################end#################!
!   call cpu_time(finished)
!   total_time=finished-started+total_time
!   call write_time(total_time)
!   write(*,*) 'Finished!'

end program main

