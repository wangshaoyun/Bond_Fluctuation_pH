module compute_energy
implicit none

save
  !
  !coulomb potential
  real*8           :: lb          !Bjerrum length
  real*8,  private :: EF          !electric field
  real*8,  private :: tol         !tolerance
  real*8,  private :: tau_rf      !ratio of time between real and fourier space
  real*8,  private :: alpha       !Ewald screening parameter alpha
  real*8,  private :: alpha2      !alpha2=alpha*alpha
  !
  !real space
  real*8,  private :: rcc         !Cut off radius of real space
  real*8,  private :: rcc2        !rcc2=rcc*rcc  !
  !reciprocal space
  integer, private :: Kmax1       !max wave number of x direction
  integer, private :: Kmax2       !max wave number of y direction 
  integer, private :: Kmax3       !max wave number of z direction
  integer, private :: K_total     !Total wave number in reciprocal space
  !
  !arrays
  real*8,  allocatable, dimension(:,:), private :: posq        
                        !array of position of charged particle
  integer, allocatable, dimension(:,:), private :: real_pair_list
                        !real space verlet list
  integer, allocatable, dimension( : )          :: charge
                        !charge number to monomer number
  real*8,  allocatable, dimension( : ), private :: exp_ksqr 
                        !coefficients in Fourier space
  real*8,  allocatable, dimension( : ), private :: real_fun 
                        !Function list of erfc and coefficients in real space

contains 

subroutine initialize_energy_parameter
  use global_variables
  implicit none

  !
  !read force parameters from file
  call read_force_parameters

  !
  !
  if ( qq /= 0 ) then
    !
    !Initialize ewald parameters and array allocate.
    call Initialize_ewald_parameters
  end if


end subroutine initialize_energy_parameter


subroutine read_force_parameters
  use global_variables
  implicit none

  open(10,file='./force_data.txt')
    read(10,*) epsilon
    read(10,*) sigma  
    read(10,*) rcl    
    read(10,*) rvl        
    read(10,*) rsk
    read(10,*) R0_2
    read(10,*) kFENE
    read(10,*) ordr(1)
    read(10,*) ordr(2)
    read(10,*) ordr(3)
    read(10,*) lb
    read(10,*) xi
    read(10,*) EF           
    read(10,*) tol
    read(10,*) tau_rf
  close(10)

end subroutine read_force_parameters




end module compute_energy


