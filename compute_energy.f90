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
  integer, allocatable, dimension( : )          :: charge
                        !charge number to monomer number
  integer, allocatable, dimension( : )          :: cell_list_q
                        !cell list of charge
  integer, allocatable, dimension( : )          :: inv_cell_list_q
                        !inverse cell list of charge
  integer, allocatable, dimension( : )          :: cell_list_r
                        !cell list in real space
  integer, allocatable, dimension( : )          :: inv_cell_list_r
                        !inverse cell list in real space
  integer, allocatable, dimension( : )          :: hoc_r     ! head of chains
  integer, allocatable, dimension( : )          :: inv_hoc_r ! head of chains
  real*8,  allocatable, dimension( : ), private :: exp_ksqr 
                        !coefficients in Fourier space
  real*8,  allocatable, dimension( : ), private :: real_fun 
                        !Function list of erfc and coefficients in real space
  real,    allocatable, dimension( : ), private :: fourier_ij 
                        !Coulomb energy of i,j in fourier space
  real,    allocatable, dimension( : ), private :: real_ij 
                        !Coulomb energy of i,j in real space

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
    !Initialize charge with lind list
    call Initialize_Charge
    !
    !Initialize cell list of charge
    call Initialize_cell_list_q
    !
    !Initialize real cell list
    call Initialize_real_cell_list
  end if

end subroutine initialize_energy_parameter


subroutine read_force_parameters
  use global_variables
  implicit none
  integer :: i, j, k, l, m, n
  real :: rl

  open(10,file='./force_data.txt')
    read(10,*) rcc 
    read(10,*) lb
    read(10,*) EF           
    read(10,*) tol
    read(10,*) tau_rf
  close(10)

  alpha = tol / (1.*rcc/2)                    !alpha = 2.5/3
  alpha2 = alpha * alpha 

!   allocate(fourier_ij(-260:260,-260:260,1:Lz))
!   m = 521*521*Lz
!   open(100, file='./data/fourier.txt') 
!     do i = 1, m
!       read(100,*) j, k, l, rl
!       fourier_ij(j,k,l) = rl
!     end do
!   close(100)

!   n = 2*nint(rcc)
!   allocate(real_ij(-n:n,-n:n,-n:n))
!   n = (n+1)*(n+1)*(n+1)
!   open(100, file ='./data/real.txt')
!     do i = 1, n
!       read(100,*) j,k,l,rl
!       real_ij(j,k,l) = rl
!     end do
!   close(100)

end subroutine read_force_parameters


subroutine Initialize_Charge
  use global_variables
  implicit none
  integer :: i, j

  allocate(charge(Nq))

  j = 0
  do i = 1, NN
    if ( pos(i,4) /= 0 ) then
      j = j + 1
      charge(j) = i
    end if
  end do

end subroutine Initialize_Charge


subroutine Initialize_cell_list_q
  use global_variables
  implicit none
  integer :: i, j, k

  allocate(cell_list_q(Nq+1)) ! the last one is head of the list
  allocate(inv_cell_list_q(Nq+1)) ! the last one is the head of the list

  cell_list_q(Nq+1) = 0

  do i = 1, Nq
    cell_list_q(i) = cell_list_q(Nq+1)
    cell_list_q(Nq+1) = i
  end do

  j = cell_list_q(Nq+1)
  k = 0
  do while( j /= 0 )
    k = k+1
    inv_cell_list_q(k) = j
    inv_cell_list_q(Nq+1) = k
    j = cell_list_q(j)
  end do

  open(100,file='cell_list_q.txt')
    do i = 1, Nq + 1
      write(100,*) cell_list_q(i), inv_cell_list_q(i)
    end do
  close(100)

end subroutine Initialize_cell_list_q


subroutine Initialize_real_cell_list
  use global_variables
  implicit none
  integer :: i, j, k


end subroutine Initialize_real_cell_list


subroutine Delta_Energy(DeltaE)
  use global_variables
  implicit none
  real*8, intent(out) :: DeltaE

  DeltaE = -1

end subroutine Delta_Energy


end module compute_energy


