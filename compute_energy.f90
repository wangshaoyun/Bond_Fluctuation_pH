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
  integer, private :: K1_half
  integer, private :: K2_half  
  integer, private :: K3_half   
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
  real,  allocatable, dimension(:,:,:), private :: fourier_ij 
                        !Coulomb energy of i,j in fourier space
  real,  allocatable, dimension(:,:,:), private :: real_ij 
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
    !pre-calculate arrays in fourier space
    call pre_calculte_fourier_space
    !
    !pre-calculate arrays in real space
    call pre_calculate_real_space
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
    read(10,*) lb
    read(10,*) EF           
    read(10,*) tol
    read(10,*) tau_rf
  close(10)



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


subroutine pre_calculate_fourier_space
  use global_variables
  implicit none
  include "fftw3.f90"
  integer :: i, j, k, k2, p, q, r
  complex(kind=8) :: exp_x,exp_y,exp_z
  real(kind=8),allocatable,dimension(:,:,:)::exp_k2
  complex(kind=8),allocatable,dimension(:,:,:)::FT_exp_k2
  integer ( kind = 8 ) plan_backward
  integer ( kind = 8 ) plan_forward
  integer ( kind = 8 ) plan

  Kmax1 = Lx
  Kmax2 = Ly
  Kmax3 = nint(Lz*Z_empty)
  if ( mod(Kmax1,2) == 0 ) then
    K1_half = Kmax1 / 2
    exp_x = cmplx(cos(pi),-sin(pi))
  else
    K1_half = (Kmax1-1) / 2
    exp_x = cmplx(cos(pi*(Kmax1-1)/Kmax1),sin(pi*(Kmax1-1)/Kmax1))
  end if
  if ( mod(Kmax2,2) == 0 ) then
    K2_half = Kmax2 / 2
    exp_y = cmplx(cos(pi),-sin(pi))
  else
    K2_half = (Kmax2-1) / 2
    exp_y = cmplx(cos(pi*(Kmax2-1)/Kmax2),sin(pi*(Kmax2-1)/Kmax2))
  end if
  if ( mod(Kmax3,2) == 0 ) then
    K3_half = Kmax3 / 2
    exp_z = cmplx(cos(pi),-sin(pi))
  else
    K3_half = (Kmax3-1) / 2
    exp_z = cmplx(cos(pi*(Kmax3-1)/Kmax3),sin(pi*(Kmax3-1)/Kmax3))
  end if

  rcc = tol*tol/pi*2
  alpha = tol / (rcc/2)                    !alpha = 2.5/3
  alpha2 = alpha * alpha 

  if (allocated(exp_k2)) deallocate( exp_k2 )
  if (allocated(FT_exp_k2)) deallocate( FT_exp_k2 )
  allocate( exp_k2(Kmax1,Kmax2,Kmax3) )
  allocate( FT_exp_k2(Kmax1,Kmax2,Kmax3) )

  do i = 1, Kmax1
    do j = 1, Kmax2
      do k = 1, Kmax3
        k2=(i-K1_half-1)**2+(i-K1_half-1)**2+(i-K1_half-1)**2
        if (k2 == 0) then
          exp_k2(i,j,k) == 0
        else
          exp_k2(i,j,k) = cmplx(exp(-1.D0*k2/4/alpha)/k2,0)
        end if
      end do
    end do
  end do

  call dfftw_plan_dft_3d_ ( plan_forward, Kmax1, Kmax2, Kmax3, exp_k2, FT_exp_k2, FFTW_FORWARD, FFTW_Estimate ) 

  call dfftw_execute_ ( plan_forward )
  
  call dfftw_destroy_plan_ ( plan_forward )

  deallocate(exp_k2)
  if (allocated(fourier_ij)) deallocate( fourier_ij )
  allocate( fourier_ij( -K1_half:(Kmax1-K1_half-1), -K2_half:(Kmax2-K2_half-1),-K3_half:(Kmax3-K3_half-1) )
  do i = 1, Kmax1
    do j = 1, Kmax2 
      do k = 1, Kmax3
        if ( i > Kmax1-K1_half-1 ) then
          p = i - Kmax1
        end if
        if ( j > Kmax2-K2_half-1 ) then
          q = j - Kmax2
        end if
        if ( k > Kmax3-K3_half-1 ) then
          r = k - Kmax3
        end if
        fourier_ij(p,q,r) = real( (exp_x**p) * (exp_y**q) * (exp_z**r) * &
          &                      FT_exp_k2(i,j,k) )
      end do
    end do
  end do
  fourier_ij = 2*pi*8/(Lx*Ly*Lz*Z_empty)*fourier_ij
  deallocate(FT_exp_k2)

end subroutine pre_calculate_fourier_space


subroutine pre_calculate_real_space
  use global_variables
  implicit none
  integer :: i, j, k
  integer :: nn, nn_half
  real*8 :: rr

  nn = floor(rcc)
  if ( mod(nn,2) == 0 ) then
    nn_half = nn / 2
  else
    nn_half = (nn-1) / 2
  end if 
  allocate( real_ij( -nn_half:nn-nn_half-1, -nn_half:nn-nn_half-1, -nn_half:nn-nn_half-1 ) )

  do i = 1, nn
    do j = 1, nn
      do k = 1, nn
        if ( i > nn-nn_half-1 ) then
          p = i - nn
        end if
        if ( j > nn-nn_half-1 ) then
          q = j - nn
        end if
        if ( k > nn-nn_half-1 ) then
          r = k - nn
        end if
        rr = sqrt(1.*p**2 + 1.*q**2 + 1.*r**2)
        real_ij(p,q,r) = erfc(alpha*rr)/rr/2
      end do
    end do
  end do

end subroutine pre_calculate_real_space



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

!   open(100,file='cell_list_q.txt')
!     do i = 1, Nq + 1
!       write(100,*) cell_list_q(i), inv_cell_list_q(i)
!     end do
!   close(100)

end subroutine Initialize_cell_list_q


subroutine Initialize_real_cell_list
  use global_variables
  implicit none
  integer :: i, j, k, l, m
  integer :: nx, ny, nz
  integer :: icelx, icely, icelz
  real*8 :: rnx,rny,rnz

  nx = int(1.*Lx/rcc)
  ny = int(1.*Ly/rcc)
  nz = int(1.*Lz/rcc)
  rnx = 1.*Lx/nx
  rny = 1.*Ly/ny
  rnz = 1.*Lz/nz
  ncel = nx*ny*nz

  allocate(hoc_r(nx,ny,nz))
  allocate(inv_hoc_r(nx,ny,nz))
  hoc_r = 0
  inv_hoc_r = 0

  allocate(cell_list_r(NN))
  allocate(inv_cell_list_r(NN))

  do i = 1, NN
    icelx = int(pos(i,1)/rnx)
    icely = int(pos(i,2)/rny)
    icelz = int(pos(i,3)/rnz)
    cell_list_r(i) = hoc(icelx,icely,icelz)
    hoc(icelx,icely,icelz) = i
  end do

  do i = 1, nx
    do j = 1, ny
      do k = 1, nz
        l = hoc_r(i,j,k)
        m = 0
        do while( l /= 0 )
          m = m + 1
          inv_cell_list_r(k) = l
          inv_hoc_r(icelx,icely,icelz) = k
          l = cell_list_r(l)
        end do
      end do
    end do
  end do

!   open(100,file='cell_list_q.txt')
!     do i = 1, NN
!       write(100,*) cell_list_q(i), inv_cell_list_q(i)
!     end do
!   close(100)

!   open(100,file='hoc_r.txt')
!     do i = 1, nx
!       do j = 1, ny
!         do k = 1, nz
!          write(100,*) hoc_r(icelx,icely,icelz),inv_hoc_r(icelx,icely,icelz)
!         end do
!       end do
!     end do
!   close(100)

end subroutine Initialize_real_cell_list


subroutine Delta_Energy(DeltaE)
  use global_variables
  implicit none
  real*8, intent(out) :: DeltaE

  DeltaE = -1

end subroutine Delta_Energy


end module compute_energy


