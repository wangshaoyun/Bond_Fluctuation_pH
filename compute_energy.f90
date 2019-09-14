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
  real*8,  private :: clx         !length of cell in x direction
  real*8,  private :: cly         !length of cell in y direction
  real*8,  private :: clz         !length of cell in z direction
  integer, private :: nclx        !number of cell in x direction
  integer, private :: ncly        !number of cell in y direction
  integer, private :: nclz        !number of cell in z direction 
  !reciprocal space
  integer, private :: Kmax1       !max wave number of x direction
  integer, private :: Kmax2       !max wave number of y direction 
  integer, private :: Kmax3       !max wave number of z direction
  integer, private :: K1_half
  integer, private :: K2_half  
  integer, private :: K3_half   
  integer, private :: K_total     !Total wave number in reciprocal space
  !
  !array of position of charged particle
  integer, allocatable, dimension(:,:,:), private :: cell_near_list 
  !
  !charge number to monomer number        
  integer, allocatable, dimension( : )          :: charge
  !  
  !charge number to monomer number                    
  integer, allocatable, dimension( : ), private :: inv_charge
  !
  !cell list of charge
  integer, allocatable, dimension( : ), private :: cell_list_q
  !
  !inverse cell list of charge
  integer, allocatable, dimension( : ), private :: inv_cell_list_q
  !
  !cell list in real space
  integer, allocatable, dimension( : ), private :: cell_list_r
  !
  !inverse cell list in real space
  integer, allocatable, dimension( : ), private :: inv_cell_list_r
  !
  !inverse cell list in real space                      
  integer, allocatable, dimension( : ), private :: inv_cell_list_r
  !
  ! head of chains
  integer, allocatable, dimension( : ), private :: hoc_r     
  !
  ! head of chains
  integer, allocatable, dimension( : ), private :: inv_hoc_r
  !
  ! Periodic condition
  integer, allocatable, dimension( : ), private :: periodic_x
  !
  !Periodic condition  
  integer, allocatable, dimension( : ), private :: periodic_y
  !
  !Coulomb energy of i,j in fourier space
  real,  allocatable, dimension(:,:,:), private :: fourier_ij 
  !
  !Coulomb energy of i,j in real space
  real,  allocatable, dimension(:,:,:), private :: real_ij 


contains 

subroutine initialize_energy_parameter
  use global_variables
  implicit none

  !
  !read force parameters from file
  call read_force_parameters
  !
  ! periodic condition array
  call Periodic_array
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


subroutine error_analysis
  use global_variables
  implicit none
  real*8 :: EE1, EE2, absolute_error, relative_error
  real*8 :: real_time, fourier_time, time_ewald
  real*8 :: st, fn

  call cpu_time(st)
  call compute_energy_Ewald(EE1)
  call cpu_time(fn)
  time_ewald = fn - st

  call compute_energy_lookup_table(EE2, real_time, fourier_time)

  absolute_error = abs(EE2-EE1)

  relative_error = absolute_error / EE1

  write(*,*) 
  write(*,*) '******************error_analysis********************'
  write(*,*) 'absolute error         absolute_error:', absolute_error
  write(*,*) 'relative error         relative_error:', relative_error
  write(*,*) 'real time                   real_time:', real_time
  write(*,*) 'fourier time             fourier_time:', fourier_time
  write(*,*) 'time ewald                 time_ewald:', time_ewald
  write(*,*) '****************************************************'
  write(*,*) 
  write(*,*) 

end subroutine error_analysis


subroutine compute_energy_Ewald(EE)
  use global_variables
  implicit none
  real*8, intent(out) :: EE

  EE = 0

end subroutine compute_energy_Ewald



subroutine Periodic_array
  use global_variables
  implicit none
  integer :: i

  allocate(periodic_x(-Lx:Lx))
  allocate(periodic_y(-Ly:Ly))
  Periodic_x = 0
  Periodic_y = 0

  if (mod(Lx,2) == 0) then
    do i = -Lx, Lx
      if (i<-Lx/2) then
        periodic_x(i) = i + Lx
      elseif (i>=Lx/2) then
        periodic_x(i) = i - Lx
      end if
    end do
  else
    do i = -Lx, Lx
      if (i<-Lx/2) then
        periodic_x(i) = i + Lx
      elseif (i>Lx/2) then
        periodic_x(i) = i - Lx
      end if
    end do
  end if 

  if (mod(Ly,2) == 0) then
    do i = -Ly, Ly
      if (i<-Ly/2) then
        periodic_y(i) = i + Ly
      elseif (i>=Ly/2) then
        periodic_y(i) = i - Ly
      end if
    end do
  else
    do i = -Ly, Ly
      if (i<-Ly/2) then
        periodic_y(i) = i + Ly
      elseif (i>Ly/2) then
        periodic_y(i) = i - Ly
      end if
    end do
  end if 

end subroutine Periodic_array


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
  allocate(inv_charge(NN))

  j = 0
  do i = 1, NN
    if ( pos(i,4) /= 0 ) then
      j = j + 1
      charge(j) = i
    end if
  end do

  j = 0
  do i = 1, NN
    if ( pos(i,4) /=0 ) then
      j = j + 1
      inv_charge(i) = j
    else
      inv_charge(i) = 0
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

  inv_cell_list_q(Nq+1) = 0
  do i = Nq, 1, -1
    inv_cell_list_q(i) = inv_cell_list_q(N+1)
    inv_cell_list_q(Nq+1) = i
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
  integer :: i, j, k, l, m, n, p, q, r, x, y, z
  integer :: icelx, icely, icelz

  nclx = int(1.*Lx/rcc)
  ncly = int(1.*Ly/rcc)
  nclz = int(1.*Lz/rcc)
  clx = 1.*Lx/nclx
  cly = 1.*Ly/ncly
  clz = 1.*Lz/nclz
  ncel = nclx*ncly*nclz

  allocate(hoc_r(nclx,ncly,nclz))
  allocate(inv_hoc_r(nclx,ncly,nclz))
  hoc_r = 0
  inv_hoc_r = 0

  allocate(cell_list_r(Nq))
  allocate(inv_cell_list_r(Nq))

  do i = 1, Nq
    j = charge(i)
    icelx = int(pos(j,1)/clx)
    icely = int(pos(j,2)/cly)
    icelz = int(pos(j,3)/clz)
    cell_list_r(i) = hoc(icelx,icely,icelz)
    hoc(icelx,icely,icelz) = i
  end do

  do i = Nq, 1, -1
    j = charge(i)
    icelx = int(pos(j,1)/clx)
    icely = int(pos(j,2)/cly)
    icelz = int(pos(j,3)/clz)
    inv_cell_list_r(i) = inv_hoc_r(icelx,icely,icelz)
    inv_hoc_r(icelx,icely,icelz) = i
  end do

  allocate(cell_near_list(nclx*ncly*nclz,28,3)))
  cell_near_list = 0
  m = 0
  do i = 1, nclx
    do j = 1, ncly
      do k = 1, nclz
        m = m + 1
        n = 0
        do p = -1, 1
          do q = -1, 1
            do r = -1, 1
              x = i + p
              y = j + q
              z = k + r
              if (z>0 .and. z<=nclz) then
                n = n + 1
                if (x>nclx) then
                  x = x - nclx
                elseif (x<=0) then
                  x = x + nclx
                end if
                if (y>ncly) then
                  y = y - ncly
                elseif (y<=0) then
                  y = y + ncly
                end if
                cell_near_list(m,n,1) = x
                cell_near_list(m,n,2) = y
                cell_near_list(m,n,3) = z
              end if
            end do
          end do
        end do
        cell_near_list(m,28,1) = n
      end do
    end do
  end do

!   open(100,file='cell_list_q.txt')
!     do i = 1, NN
!       write(100,*) cell_list_q(i), inv_cell_list_q(i)
!     end do
!   close(100)

!   open(100,file='hoc_r.txt')
!     do i = 1, nclx
!       do j = 1, ncly
!         do k = 1, nclz
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
  integer :: i,j,k,x,y,z,icelx,icey,icelz,ncel
  real*8 :: EE1,EE2,rr
  real*8, dimension(3) :: rij

  DeltaE = 0

  !
  ! Fourier Space
  EE1 = 0
  EE2 = 0
  j = cell_list_q(Nq+1)
  do while (j/=0)
    k = charge(j)
    if (k/=ip) then
      (/x,y,z/) = pos_ip0(1:3) - pos(k,1:3)
      x = Periodic_x(x)
      y = Periodic_y(y)
      EE1 = EE1 + pos(k,4)*fourier_ij(x,y,z)

      (/x,y,z/) = pos_ip1(1:3) - pos(k,1:3)
      x = Periodic_x(x)
      y = Periodic_y(y)
      EE2 = EE2 + pos(k,4)*fourier_ij(x,y,z)
    end if
  end do

  DeltaE = pos_ip1(4) * (EE2-EE1)

  !
  ! Real Space
  EE1 = 0
  icelx = int(pos_ip0(1)/clx)
  icely = int(pos_ip0(2)/cly)
  icelz = int(pos_ip0(3)/clz) 
  ncel = (icelx-1)*ncly*nclz+(icely-1)*nclz+icelz
  do i = 1, cell_near_list(ncel,28,1)
    icelx = cell_near_list(ncel,i,1)
    icely = cell_near_list(ncel,i,2)
    icelz = cell_near_list(ncel,i,3)   
    j = hoc_r(icelx,icely,icelz) 
    do while (j/=0) 
      k = charge(j)
      if (k/=ip) then
        (/x,y,z/) = pos_ip0(1:3) - pos(k,1:3)
        x = Periodic_x(x)
        y = Periodic_y(y)
        if ((x*x+y*y+z*z)<rcc) then
          EE1 = EE1 + pos(k,4)*real_ij(x,y,z)
        end if
      end if
    end do
  end do
  EE1 = EE1 * pos_ip1(4)

  EE2 = 0
  icelx = int(pos_ip1(1)/clx)
  icely = int(pos_ip1(2)/cly)
  icelz = int(pos_ip1(3)/clz) 
  ncel = (icelx-1)*ncly*nclz+(icely-1)*nclz+icelz
  do i = 1, cell_near_list(ncel,28,1)
    icelx = cell_near_list(ncel,i,1)
    icely = cell_near_list(ncel,i,2)
    icelz = cell_near_list(ncel,i,3)   
    j = hoc_r(icelx,icely,icelz) 
    do while (j/=0) 
      k = charge(j)
      if (k/=ip) then
        (/x,y,z/) = pos_ip1(1:3) - pos(k,1:3)
        x = Periodic_x(x)
        y = Periodic_y(y)
        if ((x*x+y*y+z*z)<rcc) then
          EE2 = EE2 + pos(k,4)*real_ij(x,y,z)
        end if
      end if
    end do
  end do
  EE2 = EE2 * pos_ip1(4)
  DeltaE = DeltaE + EE2 - EE1

end subroutine Delta_Energy


subroutine Delta_Energy_add(DeltaE)
  use global_variables
  implicit none
  real*8, intent(out) :: DeltaE
  integer :: i,j,k,x,y,z,qq1,qq2,icelx,icey,icelz,ncel
  real*8 :: EE1,EE2,rr
  real*8, dimension(3) :: rij

  DeltaE = 0

  !
  ! Fourier Space
  EE1 = 0
  EE2 = 0
  qq1 = pos_ip1(4)
  qq2 = pos_ip1i(4)
  j = cell_list_q(Nq+1)
  do while (j/=0)
    k = charge(j)
    (/x,y,z/) = pos_ip1(1:3) - pos(k,1:3)
    x = Periodic_x(x)
    y = Periodic_y(y)
    EE1 = EE1 + pos(k,4)*fourier_ij(x,y,z)

    (/x,y,z/) = pos_ip1i(1:3) - pos(k,1:3)
    x = Periodic_x(x)
    y = Periodic_y(y)
    EE2 = EE2 + pos(k,4)*fourier_ij(x,y,z)
  end do

  DeltaE = EE1*qq1 + EE2*qq2

  !
  ! Real Space
  EE1 = 0
  icelx = int(pos_ip1(1)/clx)
  icely = int(pos_ip1(2)/cly)
  icelz = int(pos_ip1(3)/clz) 
  ncel = (icelx-1)*ncly*nclz+(icely-1)*nclz+icelz
  do i = 1, cell_near_list(ncel,28,1)
    icelx = cell_near_list(ncel,i,1)
    icely = cell_near_list(ncel,i,2)
    icelz = cell_near_list(ncel,i,3)   
    j = hoc_r(icelx,icely,icelz) 
    do while (j/=0) 
      k = charge(j)
      (/x,y,z/) = pos_ip1(1:3) - pos(k,1:3)
      x = Periodic_x(x)
      y = Periodic_y(y)
      if ((x*x+y*y+z*z)<rcc) then
        EE1 = EE1 + pos(k,4)*real_ij(x,y,z)
      end if
    end do
  end do
  EE1 = EE1 * qq1

  EE2 = 0
  icelx = int(pos_ip1i(1)/clx)
  icely = int(pos_ip1i(2)/cly)
  icelz = int(pos_ip1i(3)/clz) 
  ncel = (icelx-1)*ncly*nclz+(icely-1)*nclz+icelz
  do i = 1, cell_near_list(ncel,28,1)
    icelx = cell_near_list(ncel,i,1)
    icely = cell_near_list(ncel,i,2)
    icelz = cell_near_list(ncel,i,3)   
    j = hoc_r(icelx,icely,icelz) 
    do while (j/=0) 
      k = charge(j)
      (/x,y,z/) = pos_ip1i(1:3) - pos(k,1:3)
      x = Periodic_x(x)
      y = Periodic_y(y)
      if ((x*x+y*y+z*z)<rcc) then
        EE2 = EE2 + pos(k,4)*real_ij(x,y,z)
      end if
    end do
  end do
  EE2 = EE2 * qq1
  DeltaE = DeltaE + EE1 + EE2

  !
  !interaction of the added two particles
  (/x,y,z/) = pos_ip1i(1:3) - pos_ip1(1:3)
  x = Periodic_x(x)
  y = Periodic_y(y)
  DeltaE = DeltaE + qq1*qq2*(fourier_ij(x,y,z)+real_ij(x,y,z))

end subroutine Delta_Energy_add


subroutine Delta_Energy_delete(DeltaE)
  use global_variables
  implicit none
  real*8, intent(out) :: DeltaE
  integer :: i,j,k,x,y,z,qq1,qq2,icelx,icey,icelz,ncel
  real*8 :: EE1,EE2,rr
  real*8, dimension(3) :: rij

  DeltaE = 0

  !
  ! Fourier Space
  EE1 = 0
  EE2 = 0
  qq1 = pos_ip1(4)
  qq2 = pos_ip1i(4)
  j = cell_list_q(Nq+1)
  do while (j/=0)
    k = charge(j)
    if (k/=ip) then
      (/x,y,z/) = pos_ip1(1:3) - pos(k,1:3)
      x = Periodic_x(x)
      y = Periodic_y(y)
      EE1 = EE1 + pos(k,4)*fourier_ij(x,y,z)
    end if
    if (k/=ip1) then
      (/x,y,z/) = pos_ip1i(1:3) - pos(k,1:3)
      x = Periodic_x(x)
      y = Periodic_y(y)
      EE2 = EE2 + pos(k,4)*fourier_ij(x,y,z)
    end if
  end do

  DeltaE = EE1*qq1 + EE2*qq2

  !
  ! Real Space
  EE1 = 0
  icelx = int(pos_ip1(1)/clx)
  icely = int(pos_ip1(2)/cly)
  icelz = int(pos_ip1(3)/clz) 
  ncel = (icelx-1)*ncly*nclz+(icely-1)*nclz+icelz
  do i = 1, cell_near_list(ncel,28,1)
    icelx = cell_near_list(ncel,i,1)
    icely = cell_near_list(ncel,i,2)
    icelz = cell_near_list(ncel,i,3)   
    j = hoc_r(icelx,icely,icelz) 
    do while (j/=0) 
      k = charge(j)
      if (k/=ip) then
        (/x,y,z/) = pos_ip1(1:3) - pos(k,1:3)
        x = Periodic_x(x)
        y = Periodic_y(y)
        if ((x*x+y*y+z*z)<rcc) then
          EE1 = EE1 + pos(k,4)*real_ij(x,y,z)
        end if
      end if
    end do
  end do
  EE1 = EE1 * qq1

  EE2 = 0
  icelx = int(pos_ip1i(1)/clx)
  icely = int(pos_ip1i(2)/cly)
  icelz = int(pos_ip1i(3)/clz) 
  ncel = (icelx-1)*ncly*nclz+(icely-1)*nclz+icelz
  do i = 1, cell_near_list(ncel,28,1)
    icelx = cell_near_list(ncel,i,1)
    icely = cell_near_list(ncel,i,2)
    icelz = cell_near_list(ncel,i,3)   
    j = hoc_r(icelx,icely,icelz) 
    do while (j/=0) 
      k = charge(j)
      if (k/=ip1) then
        (/x,y,z/) = pos_ip1i(1:3) - pos(k,1:3)
        x = Periodic_x(x)
        y = Periodic_y(y)
        if ((x*x+y*y+z*z)<rcc) then
          EE2 = EE2 + pos(k,4)*real_ij(x,y,z)
        end if
      end if
    end do
  end do
  EE2 = EE2 * qq1
  DeltaE = DeltaE + EE1 + EE2
  !
  !interaction of the added two particles
  (/x,y,z/) = pos_ip1i(1:3) - pos_ip1(1:3)
  x = Periodic_x(x)
  y = Periodic_y(y)
  DeltaE = DeltaE - qq1*qq2*(fourier_ij(x,y,z)+real_ij(x,y,z))

  DeltaE = -DeltaE

end subroutine Delta_Energy_delete


subroutine update_real_cell_list
  use global_variables
  implicit none

  call update_real_cell_list_delete
  call update_real_cell_list_add

end subroutine update_real_cell_list


subroutine update_real_cell_list_add
  use global_variables
  implicit none
  integer :: icelx, icely, icelz
  integer :: nti,bfi,ii       ! next particle of ii, before particle of ii
  integer :: ed, st 

  icelx = int(pos_ip1(1)/clx)
  icely = int(pos_ip1(2)/cly)
  icelz = int(pos_ip1(3)/clz)  

  ii = inv_charge(ip)   !ii belongs to [1,Nq]

  cell_list_r(ii) = hoc_r(icelx,icely,icelz)
  hoc_r(icelx,icely,icelz) = ii

  inv_cell_list_r(ii) = 0
  if ( inv_hoc_r(icelx,icely,icelz) /=0 ) then
    inv_cell_list_r( hoc_r(icelx,icely,icelz) ) = ii
  else
    inv_hoc_r(icelx,icely,icelz) = ii
  end if

end subroutine update_real_cell_list_add


subroutine update_real_cell_list_delete
  use global_variables
  implicit none
  integer :: icelx, icely, icelz
  integer :: nti,bfi,ii       ! next particle of ii, before particle of ii
  integer :: ed, st 

  icelx = int(pos_ip0(1)/clx)
  icely = int(pos_ip0(2)/cly)
  icelz = int(pos_ip0(3)/clz) 

  ii = inv_charge(ip)   !ii belongs to [1,Nq]

  nti = cell_list_r(ii)
  bfi = inv_cell_list_r(ii)

  if ( bfi/=0 .and. nti/=0 ) then        !middle
    cell_list_r(nti) = bfi
    inv_cell_list_r(bfi) = nti
  elseif ( bfi==0 .and. nti/=0 ) then    !the first one
    cell_list_r(nti) = bfi
    inv_hoc_r(icelx,icely,icelz) = nti
  elseif ( bfi/=0 .and. nti==0 ) then
    hoc_r(icelx,icely,icelz) = bfi
    inv_cell_list_q(bfi) = nti
  else
    hoc_r(icelx,icely,icelz) = bfi
    inv_hoc_r(icelx,icely,icelz) = nti
  end if

end subroutine update_real_cell_list_delete


subroutine update_charge_cell_list_add
  use global_variables
  implicit none
  integer :: ii       

  ii = inv_charge(ip)         ! ii belongs to [1,Nq]

  cell_list_q(ii) = cell_list_q(Nq+1)
  cell_list_q(Nq+1) = ii

  inv_cell_list_q(ii) = 0
  if ( cell_list_q(Nq+1)/=0 ) then
    inv_cell_list_q(cell_list_q(Nq+1)) = ii
  else
    inv_cell_list_q(Nq+1) = ii
  end if

end subroutine update_charge_cell_list_add


subroutine update_charge_cell_list_delete
  use global_variables
  implicit none
  integer :: nti,bfi,ii       ! next particle of ii, before particle of ii

  ii = inv_charge(ip)         !ii belongs to [1,Nq]
  
  bfi = cell_list_q(ii)
  nti = inv_cell_list_q(ii)

  if ( bfi/=0 .and. nti/=0 ) then        !middle
    cell_list_q(nti) = bfi
    inv_cell_list_q(bfi) = nti
  elseif ( bfi==0 .and. nti/=0 ) then    !the first one
    cell_list_q(nti) = bfi
    inv_cell_list_q(Nq+1) = nti
  elseif ( bfi/=0 .and. nti==0 ) then
    cell_list_q(Nq+1) = bfi
    inv_cell_list_q(bfi) = nti
  else
    cell_list_q(Nq+1) = bfi
    inv_cell_list_q(Nq+1) = nti
  end if

end subroutine update_charge_cell_list_delete



end module compute_energy


