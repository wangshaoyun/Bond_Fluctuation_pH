module global_variables
implicit none
save
!########################constants#########################!
  real*8, parameter:: pi=3.141592653589793D0      
  real*8, parameter:: gamma=.5772156649015329D0   !Euler Gamma
!########################constants#########################!

!#####################system parameters####################!
  integer :: Npe      !Total monomers in Polyelectrolytes(PE)
  integer :: arm      !Arms of star brushes
                      !including the chains anchored to the plate
  integer :: Nma      !Monomers of each arm
  integer :: Ns       !Monomers of each star
  integer :: Nga      !chains grafted on plate, must square of an integer
  integer :: Nq       !Total charge in the system, ions + aions
  integer :: Nq_PE    !Charged monomers of PE
  integer :: Nq_net   !Net charged monomers of PE, or protonated monomers
  integer :: NN       !Total particles in the system
  integer :: NN_net   !Net particles in system
  integer :: Nq_salt_ions !Charged salt ions, which not include anions.
  integer :: man_s    !Manning effect: star chains
  real*8  :: Lx       !length in x direction
  real*8  :: Ly       !length number in y direction
  real*8  :: Lz       !length number between two plate
  integer :: Lx2      !Lattice number in x direction
  integer :: Ly2      !Lattice number in y direction
  integer :: Lz2      !Lattice number between two plate
  integer :: N_bond   !Number of all bonds in system
  integer :: qq       !Charge of charged monomers
  integer :: qqi      !Charge of salt ions
  real*8  :: ion_ratio!Ratio of salt ions to the charge quantites of PE
  real*8  :: Z_empty  !Empty space ratio of height and length in slab geometry
  real*8  :: sigmag   !Grafting density of brushes on the plate (sigma unit)
  real*8  :: sigmag1  !True sigma of the lattice system
  real*8  :: Beta     !Beta=1/(kB*T), T is temperature, 
                      !kB is Boltzmann constant
  real*8  :: pH_pKa   !pH - pKa
  integer :: accept   !accepted move
  real*8  :: accept_ratio 
!##########################################################!

!##################running and Histogram###################!
  integer :: restart_or_continue  !restart or continue after breaking off
  integer :: StepNum0             !Steps of preheating
  integer :: StepNum              !Steps of running
  integer :: Deltastep            !Each deltastep of pH move
  integer :: DeltaStep1           !step inteval, height
  integer :: DeltaStep2           !step inteval, histogram
  integer :: DeltaStep3           !step inteval, write data
  integer :: step                 !ordering number
  !
  !timing
  real*8  :: started              !time at starting
  real*8  :: finished             !time at finishing
  real*8  :: total_time=0         !total time of the simulation
  !
  !histogram
  integer :: SizeHist=500         !number of histogram which is equally divided
!################end running and Histogram#################!

!##########################arrays##########################!
  integer, allocatable, dimension(:,:) :: pos    !array of position
  integer, dimension(4) :: pos_ip0                !old position of ip
  integer, dimension(4) :: pos_ip1                !new position of ip
  integer, dimension(4) :: pos_ip0i               !old position of ions
  integer, dimension(4) :: pos_ip1i               !new position of ions
  integer :: ip                                   !The particle that is choosed
  integer :: ip1                                  !The ions choosed
  integer, dimension(6,3) :: new_direction
  !
  ! bonds connected to each monomer, the direction of first one is from before 
  ! to itself, ohters if from itself to the next ones 
  integer, allocatable, dimension(:,:) :: monbd  
  !       
  ! bonds number of each bonds
  integer, allocatable, dimension(:)   :: bond_numb 
  !  
  ! status of occupation of lattice vertexes  
  integer(kind=1), allocatable, dimension(:,:,:):: latt 
  !                                                      
  ! i is lattice number, or the 
  ! left lattice point, i_plus is the right lattice point
  integer, allocatable, dimension(:) :: ipx   
  integer, allocatable, dimension(:) :: ipy
  integer, allocatable, dimension(:) :: ipz
  !
  ! ip2x equals current lattice vertex plus 2, which means the right lattice 
  ! vertex after forward move  
  integer, allocatable, dimension(:) :: ip2x  
  integer, allocatable, dimension(:) :: ip2y
  integer, allocatable, dimension(:) :: ip2z
  !
  ! imx equals current lattice vertex minus 1, which means the left lattice 
  ! vertex after backward move
  integer, allocatable, dimension(:) :: imx  
  integer, allocatable, dimension(:) :: imy
  integer, allocatable, dimension(:) :: imz
  !
  ! record 108 probable bond vector
  integer, allocatable, dimension(:,:) :: bonds  
  !
  ! bond number after 6 kinds of move
  integer, allocatable, dimension(:,:) :: move    
  ! 
  ! bonds length array
  real*8, allocatable, dimension(:) :: b1   
  !
  ! square of bonds length array
  real*8, allocatable, dimension(:) :: b12  
!########################end arrays########################!


contains

subroutine periodic_condition(rr)
  !--------------------------------------!
  !Peridodic condition of position vector
  !rr(2) in slab geometry.
  !   
  !Input
  !   rr
  !Output
  !   rr
  !External Variables
  !   Lx, Ly
  !Routine Referenced:
  !1.
  !--------------------------------------!
  implicit none
  integer, intent(inout) :: rr(2)

  if ( rr(1) > Lx2 ) then
    rr(1) = rr(1) - Lx2
  elseif( rr(1) <= 0 ) then
    rr(1) = rr(1) + Lx2
  end if
  if ( rr(2) > Ly2 ) then
    rr(2) = rr(2) - Ly2
  elseif( rr(2) <= 0 ) then
    rr(2) = rr(2) + Ly2
  end if

end subroutine periodic_condition

end module global_variables 


