module compute_energy
implicit none

save
  !
  !coulomb potential
  real*8,  private :: lb          !Bjerrum length
  real*8,  private :: EF          !electric field
  real*8,  private :: tol         !tolerance
  real*8, private  :: tau_rf      !time ratio of real space and fourier space
  real*8,  private :: alpha       !Ewald screening parameter alpha
  real*8,  private :: alpha2      !alpha2=alpha*alpha
  real*8,  private :: Mz          !total dipole moment
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
  !neighbor cells of the center cell
  integer, allocatable, dimension(:,:,:), private :: cell_near_list 
  !
  !charge number to monomer number        
  integer, allocatable, dimension( : )          :: charge
  !  
  !monomer number to charge number
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
  ! head of chains, cell list
  integer, allocatable, dimension(:,:,:), private :: hoc_r     
  !
  ! head of chains, inverse cell list
  integer, allocatable, dimension(:,:,:), private :: inv_hoc_r
  !
  ! Periodic condition
  integer, allocatable, dimension( : ), private :: periodic_x
  !
  !Periodic condition  
  integer, allocatable, dimension( : ), private :: periodic_y
  !
  !Periodic condition  
  integer, allocatable, dimension( : ), private :: periodic_z
  !
  !Periodic condition  
  integer, allocatable, dimension( : ), private :: fourier_x
  !
  !Periodic condition  
  integer, allocatable, dimension( : ), private :: fourier_y
  !
  !Periodic condition  
  integer, allocatable, dimension( : ), private :: fourier_z
  !
  !Coulomb energy of i,j in fourier space
  real,  allocatable, dimension(:), private :: fourier_ij 
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
  if ( qq /= 0 ) then
    !
    !pre-calculate arrays in fourier space
    call pre_calculate_fourier_space
    !
    ! periodic condition array
    call Periodic_array
    !
    !pre-calculate arrays in real space
    call pre_calculate_real_space
  end if

end subroutine initialize_energy_parameter


subroutine initialize_energy_arrays
  use global_variables
  implicit none
  
  if ( qq /= 0 ) then
    !
    !Initialize charge with lind list. From this subroutine, pos array is needed.
    call Initialize_Charge
    !
    !Initialize cell list of charge
    call Initialize_cell_list_q
    !
    !Initialize real cell list
    call Initialize_real_cell_list
  end if

  call write_energy_parameters

end subroutine initialize_energy_arrays


subroutine read_force_parameters
  use global_variables
  implicit none

  open(10,file='./energy_data.txt')
    read(10,*) lb
    read(10,*) EF           
    read(10,*) tol
    read(10,*) tau_rf
  close(10)

end subroutine read_force_parameters


subroutine write_energy_parameters

  write(*,*) '****************************************************************'
  write(*,*) '*************************system_data****************************'
  write(*,*) 'rcc           :', rcc
  write(*,*) 'nclx          :', nclx
  write(*,*) 'ncly          :', ncly
  write(*,*) 'nclz          :', nclz
  write(*,*) 'clx           :', clx
  write(*,*) 'cly           :', cly
  write(*,*) 'clz           :', clz
  write(*,*) '****************************************************************'
  write(*,*)
  write(*,*)
end subroutine write_energy_parameters


subroutine Periodic_array
  use global_variables
  implicit none
  integer :: i

  allocate(periodic_x(-Lx2:Lx2))
  allocate(periodic_y(-Ly2:Ly2))
  allocate(periodic_z(-Lz2:Lz2))
  allocate(fourier_x(-Lx2:Lx2))
  allocate(fourier_y(-Ly2:Ly2))
  allocate(fourier_z(-Lz2:Lz2))
  Periodic_x = 0
  Periodic_y = 0
  fourier_x = 0
  fourier_y = 0
  fourier_z = 0

  if (mod(Lx2,2) == 0) then
    do i = -Lx2, Lx2
      if (i<-Lx2/2) then
        periodic_x(i) = i + Lx2
      elseif (i>=Lx2/2) then
        periodic_x(i) = i - Lx2
      else
        periodic_x(i) = i
      end if
    end do
  else
    do i = -Lx2, Lx2
      if (i<-Lx2/2) then
        periodic_x(i) = i + Lx2
      elseif (i>Lx2/2) then
        periodic_x(i) = i - Lx2
      else
        periodic_x(i) = i
      end if
    end do
  end if
  do i = -Lx2, Lx2
    if (Periodic_x(i)<0) then
      periodic_x(i) = -Periodic_x(i)
    end if
  end do

  if (mod(Ly2,2) == 0) then
    do i = -Ly2, Ly2
      if (i<-Ly2/2) then
        periodic_y(i) = i + Ly2
      elseif (i>=Ly2/2) then
        periodic_y(i) = i - Ly2
      else
        Periodic_y(i) = i
      end if
    end do
  else
    do i = -Ly2, Ly2
      if (i<-Ly2/2) then
        periodic_y(i) = i + Ly2
      elseif (i>Ly2/2) then
        periodic_y(i) = i - Ly2
      else
        Periodic_y(i) = i
      end if
    end do
  end if 
  do i = -Ly2, Ly2
    if (Periodic_y(i)<0) then
      periodic_y(i) = -Periodic_y(i)
    end if
  end do

  do i = -Lz2, Lz2
    if (i>=0) then
      Periodic_z(i) = i;
    else
      Periodic_z(i) = -i;
    end if 
  end do

  do i = -Lx2, Lx2
    fourier_x(i) = Periodic_x(i) + 1
    fourier_x(i) = (fourier_x(i)-1)*(K2_half+1)*(Lz2+1)
  end do

  do i = -Ly2, Ly2
    fourier_y(i) = Periodic_y(i) + 1
    fourier_y(i) = (fourier_y(i)-1)*(Lz2+1)
  end do

  do i = -Lz2, Lz2
    fourier_z(i) = Periodic_z(i) + 1
  end do


!   open(100,file='./data/periodic_xy.txt')
!     do i = -Ly2, Ly2
!       write(100,*) i,periodic_x(i),Periodic_y(i)
!     end do
!   close(100)

end subroutine Periodic_array


subroutine energy_lookup_table(EE, rt, ft)
  use global_variables
  implicit none
  real*8, intent(out) :: EE
  real*8, intent(out) :: rt
  real*8, intent(out) :: ft
  integer :: i, j, k, l, m, n, x1, y1, z1, x, y, z, t
  integer :: icelx, icely, icelz, ncel
  real*8 :: st, fn, EE1, EE2,rr(4)

  EE = 0
  !
  !real space
  call cpu_time(st)
  do i = 1, Nq
    m = charge(i)
    icelx = int(pos(m,1)/clx)+1
    icely = int(pos(m,2)/cly)+1
    icelz = int(pos(m,3)/clz)+1
    ncel = (icelx-1)*ncly*nclz+(icely-1)*nclz+icelz
    do j = 1, cell_near_list(ncel,28,1)
      icelx = cell_near_list(ncel,j,1)
      icely = cell_near_list(ncel,j,2)
      icelz = cell_near_list(ncel,j,3)
      k = hoc_r(icelx,icely,icelz)
      do while(k/=0)
        l = charge(k)
        if (l/=m) then
          x = pos(m,1) - pos(l,1)
          y = pos(m,2) - pos(l,2)
          z = pos(m,3) - pos(l,3)
          if ((x*x+y*y+z*z)<rcc2) then
            EE = EE + pos(m,4)*pos(l,4)*real_ij(periodic_x(x),periodic_y(y),Periodic_z(z))
          end if
        end if
        k = cell_list_r(k)
      end do
    end do
  end do
  call cpu_time(fn)
  rt = fn - st 

  EE1 = EE
  write(*,*) 'real energy', EE, rt
  !
  !fourier space
  call cpu_time(st)
  write(*,*) size(fourier_ij)
  do i = 1, Nq
    m = charge(i)
    EE2=0
    x = pos(m,1)
    y = pos(m,2)
    z = pos(m,3)
    do j = 1, Nq
      n = charge(j)
      if (m/=n) then
        x1 = x - pos(n,1)
        y1 = y - pos(n,2)
        z1 = z - pos(n,3)  
        t = fourier_x(x1)+fourier_y(y1)+fourier_z(z1)
        if (t<1 .or. t>(K1_half+1)*(K2_half+1)*(Lz2+1)) then
          write(*,*) t,x1,y1,z1,fourier_x(x1),fourier_y(y1),fourier_z(z1)
        end if
        EE2=EE2+pos(n,4)*fourier_ij(fourier_x(x1)+fourier_y(y1)+fourier_z(z1))
      end if
    end do
    EE2= pos(m,4)*EE2
    EE = EE + EE2
  end do
  call cpu_time(fn)
  ft = fn - st
  write(*,*) 'fourier energy', EE-EE1,ft
  write(*,*) 'coulomb energy', EE
  !
  !modified term of slab geometry
  EE = EE/2 + 2*pi / (Lx*Ly*Lz*Z_empty) * lb/Beta * Mz**2
  write(*,*) 'coulomb energy with modification', EE
  stop
end subroutine energy_lookup_table


subroutine pre_calculate_fourier_space
  use global_variables
  implicit none
  include "fftw3.f90"
  integer :: i, j, k, p, q, r, t
  real*8 :: kx,ky,kz,krx,kry,krz,k2
  complex(kind=8) :: exp_x,exp_y,exp_z,test,a1,a2,a3
  complex(kind=8), allocatable, dimension(:,:,:) :: exp_k2
  complex(kind=8), allocatable, dimension(:,:,:) :: FT_exp_k2
  real*8, allocatable,dimension(:,:,:) :: fourier_ij1
  integer ( kind = 8 ) plan_backward
  integer ( kind = 8 ) plan_forward
  integer ( kind = 8 ) plan
!   complex(kind=8) :: fij
!   integer :: rr(3)

  alpha = pi/2/tol
  alpha2 = alpha*alpha
  rcc = 1.D0*floor(2*tol*tol/pi*2)       ! 2*2.5*2.5/pi*2=7.96, lattice unit
  rcc2 = rcc*rcc                         ! rcc is dependent on tol only.
  write(*,*) 'rcc=',rcc

  Kmax1 = Lx2
  Kmax2 = Ly2
  Kmax3 = nint(Lz2*Z_empty)
  if ( mod(Kmax1,2) == 0 ) then
    K1_half = Kmax1 / 2
    exp_x = cmplx(-1,0)
  else
    K1_half = (Kmax1-1) / 2
    exp_x = cmplx(cos(pi*(Kmax1-1)/Kmax1),sin(pi*(Kmax1-1)/Kmax1))
  end if
  if ( mod(Kmax2,2) == 0 ) then
    K2_half = Kmax2 / 2
    exp_y = cmplx(-1,0)
  else
    K2_half = (Kmax2-1) / 2
    exp_y = cmplx(cos(pi*(Kmax2-1)/Kmax2),sin(pi*(Kmax2-1)/Kmax2))
  end if
  if ( mod(Kmax3,2) == 0 ) then
    K3_half = Kmax3 / 2
    exp_z = cmplx(-1,0)
  else
    K3_half = (Kmax3-1) / 2
    exp_z = cmplx(cos(pi*(Kmax3-1)/Kmax3),sin(pi*(Kmax3-1)/Kmax3))
  end if

  !
  !for the maxium situation sigmag=1e-3, (500,500,1200), it will need 9.6G RAM
  if (allocated(exp_k2)) deallocate( exp_k2 )
  if (allocated(FT_exp_k2)) deallocate( FT_exp_k2 )
  allocate( exp_k2(Kmax1,Kmax2,Kmax3) )
  allocate( FT_exp_k2(Kmax1,Kmax2,Kmax3) )

!   fij=0
!   rr=(/-60,-60,-200/)
  do i = 1, Kmax1
    do j = 1, Kmax2
      do k = 1, Kmax3
        kx = 2*pi*(i-K1_half-1)/Kmax1
        ky = 2*pi*(j-K2_half-1)/Kmax2
        kz = 2*pi*(k-K3_half-1)/Kmax3
        k2 = kx*kx + ky*ky + kz*kz
        if (k2 == 0) then
          exp_k2(i,j,k) = 0
        else
          exp_k2(i,j,k) = cmplx(exp(-k2/4/alpha2)/k2,0)
        end if
!         krx = kx*rr(1)
!         kry = ky*rr(2)
!         krz = kz*rr(3)
!         fij = fij + exp_k2(i,j,k)*cmplx(cos(krx+kry+krz),sin(krx+kry+krz))
      end do
    end do
  end do
!   fij = (lb/beta*4*pi/(Lx*Ly*Lz*Z_empty))*fij

  call dfftw_plan_dft_3d_ ( plan_forward, Kmax1, Kmax2, Kmax3, exp_k2, FT_exp_k2, FFTW_FORWARD, FFTW_Estimate ) 

  call dfftw_execute_ ( plan_forward )
  
  call dfftw_destroy_plan_ ( plan_forward )

  deallocate(exp_k2)

  if (allocated(fourier_ij1)) deallocate( fourier_ij1 )
  allocate( fourier_ij1( -K1_half:(Kmax1-K1_half-1), -K2_half:(Kmax2-K2_half-1),-K3_half:(Kmax3-K3_half-1) ) )

  do i = 1, Kmax1
    do j = 1, Kmax2 
      do k = 1, Kmax3
        p=i-1
        q=j-1
        r=k-1
        if ( p > Kmax1-K1_half-1 ) then
          p = p - Kmax1
        end if
        if ( q > Kmax2-K2_half-1 ) then
          q = q - Kmax2
        end if
        if ( r > Kmax3-K3_half-1 ) then
          r = r - Kmax3
        end if
        fourier_ij1(p,q,r)=real(exp_x**(-p)*exp_y**(-q)*exp_z**(-r)* &
          FT_exp_k2(i,j,k))
      end do
    end do
  end do

  deallocate(FT_exp_k2)

  fourier_ij1 = (lb/beta*4*pi/(Lx*Ly*Lz*Z_empty))*fourier_ij1

  if (allocated(fourier_ij)) deallocate( fourier_ij )
  allocate( fourier_ij( (K1_half+1)*(K2_half+1)*(Lz2+1) ) )  

  do i = 0, K1_half
    do j = 0, K2_half
      do k = 0, Lz2
        p = i + 1
        q = j + 1
        r = k + 1
        t = (p-1)*(K2_half+1)*(Lz2+1) + (q-1)*(Lz2+1) + r
        fourier_ij(t) = fourier_ij1(-i,-j,-k)
      end do
    end do
  end do
  deallocate(fourier_ij1)

end subroutine pre_calculate_fourier_space


subroutine pre_calculate_real_space
  use global_variables
  implicit none
  integer :: i, j, k 
  integer :: nnl, nn_half
  real*8 :: rr, x, y, z

  nnl = nint(rcc)   ! cut off length, lattice unit
  if (allocated(real_ij)) deallocate(real_ij)
  allocate( real_ij(0:nnl, 0:nnl, 0:nnl) )

  do i = 0, nnl
    do j = 0, nnl
      do k = 0, nnl
        x = i/2.D0                      !sigma unit
        y = j/2.D0
        z = k/2.D0
        rr = sqrt(x*x+y*y+z*z)
        real_ij(i,j,k) = erfc(alpha*rr)/rr
      end do
    end do
  end do

  real_ij = real_ij * lb / beta

end subroutine pre_calculate_real_space



subroutine Initialize_Charge
  use global_variables
  implicit none
  integer :: i, j

  Mz = 0
  allocate(charge(Nq))
  allocate(inv_charge(NN))

  j = 0
  do i = 1, NN
    if ( pos(i,4) /= 0 ) then
      j = j + 1
      charge(j) = i
      Mz = Mz + pos(i,4)*pos(i,3)
    end if
  end do

  j = 0
  do i = 1, NN
    if ( pos(i,4) /= 0 ) then
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

  allocate(cell_list_q(Nq+1))     ! the last one is head of the list
  allocate(inv_cell_list_q(Nq+1)) ! the last one is the head of the list

  !assume initial state, all particles are charged.
  cell_list_q(Nq+1) = 0
  do i = 1, Nq
    cell_list_q(i) = cell_list_q(Nq+1)
    cell_list_q(Nq+1) = i
  end do

  inv_cell_list_q(Nq+1) = 0
  do i = Nq, 1, -1
    inv_cell_list_q(i) = inv_cell_list_q(Nq+1)
    inv_cell_list_q(Nq+1) = i
  end do

!   open(100,file='./data/cell_list_q.txt')
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

  nclx = int(Lx2/(rcc+1))     !cell numbers in x direction
  ncly = int(Ly2/(rcc+1))
  nclz = int(Lz2/(rcc+1))
  clx = 1.D0*Lx2/nclx         !cell length    
  cly = 1.D0*Ly2/ncly
  clz = 1.D0*Lz2/nclz

  !
  ! maxium situation, (125,125,100), 6.2Mb RAM is needed.
  allocate(hoc_r(nclx,ncly,nclz))
  allocate(inv_hoc_r(nclx,ncly,nclz))
  hoc_r = 0
  inv_hoc_r = 0

  allocate(cell_list_r(Nq))
  allocate(inv_cell_list_r(Nq))

  do i = 1, Nq
    j = charge(i)
    icelx = int(pos(j,1)/clx) + 1
    icely = int(pos(j,2)/cly) + 1
    icelz = int(pos(j,3)/clz) + 1
    cell_list_r(i) = hoc_r(icelx,icely,icelz)
    hoc_r(icelx,icely,icelz) = i
  end do

  do i = Nq, 1, -1
    j = charge(i)
    icelx = int(pos(j,1)/clx) + 1
    icely = int(pos(j,2)/cly) + 1
    icelz = int(pos(j,3)/clz) + 1
    inv_cell_list_r(i) = inv_hoc_r(icelx,icely,icelz)
    inv_hoc_r(icelx,icely,icelz) = i
  end do

  !
  ! maxium situation, (125*125*100,28,3), 500Mb RAM is needed.
  allocate(cell_near_list(nclx*ncly*nclz,28,3))
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

!   open(100,file='./data/cell_list_r.txt')
!     do i = 1, NN
!       write(100,*) i,cell_list_r(i), inv_cell_list_r(i)
!     end do
!   close(100)

!   open(100,file='./data/hoc_r.txt')
!     do i = 1, nclx
!       do j = 1, ncly
!         do k = 1, nclz
!          write(100,*) i,j,k,hoc_r(i,j,k),inv_hoc_r(i,j,k)
!         end do
!       end do
!     end do
!   close(100)

end subroutine Initialize_real_cell_list


subroutine Delta_Energy(DeltaE)
  use global_variables
  implicit none
  real*8, intent(out) :: DeltaE
  integer :: i,j,k,x,y,z,x1,y1,z1,t,icelx,icely,icelz,ncel
  real*8 :: EE1,EE2,rr,dMz
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
      x1 = pos_ip0(1) - pos(k,1)
      y1 = pos_ip0(2) - pos(k,2)
      z1 = pos_ip0(3) - pos(k,3)     
      x1 = fourier_x(x1)
      y1 = fourier_y(y1)
      z1 = fourier_z(z1)
      t = x1 + y1 + z1
      EE1 = EE1 + pos(k,4)*fourier_ij(t)

      x1 = pos_ip1(1) - pos(k,1)
      y1 = pos_ip1(2) - pos(k,2)
      z1 = pos_ip1(3) - pos(k,3)
      x1 = fourier_x(x1)
      y1 = fourier_y(y1)
      z1 = fourier_z(z1)
      t = x1 + y1 + z1
      EE2 = EE2 + pos(k,4)*fourier_ij(t)
    end if
  end do

  DeltaE = pos_ip1(4) * (EE2-EE1)

  !
  ! Real Space
  EE1 = 0
  icelx = int(pos_ip0(1)/clx)+1
  icely = int(pos_ip0(2)/cly)+1
  icelz = int(pos_ip0(3)/clz)+1
  ncel = (icelx-1)*ncly*nclz+(icely-1)*nclz+icelz
  do i = 1, cell_near_list(ncel,28,1)
    icelx = cell_near_list(ncel,i,1)
    icely = cell_near_list(ncel,i,2)
    icelz = cell_near_list(ncel,i,3)   
    j = hoc_r(icelx,icely,icelz) 
    do while (j/=0) 
      k = charge(j)
      if (k/=ip) then
        x = pos_ip0(1) - pos(k,1)
        y = pos_ip0(2) - pos(k,2)
        z = pos_ip0(3) - pos(k,3)
        if ((x*x+y*y+z*z)<rcc2) then
          EE1 = EE1 + pos_ip0(4)*pos(k,4)*real_ij(periodic_x(x),periodic_y(y),Periodic_z(z))
        end if
      end if
    end do
  end do
  EE1 = EE1 * pos_ip1(4)

  EE2 = 0
  icelx = int(pos_ip1(1)/clx)+1
  icely = int(pos_ip1(2)/cly)+1
  icelz = int(pos_ip1(3)/clz)+1
  ncel = (icelx-1)*ncly*nclz+(icely-1)*nclz+icelz
  do i = 1, cell_near_list(ncel,28,1)
    icelx = cell_near_list(ncel,i,1)
    icely = cell_near_list(ncel,i,2)
    icelz = cell_near_list(ncel,i,3)   
    j = hoc_r(icelx,icely,icelz) 
    do while (j/=0) 
      k = charge(j)
      if (k/=ip) then
        x = pos_ip1(1) - pos(k,1)
        y = pos_ip1(2) - pos(k,2)
        z = pos_ip1(3) - pos(k,3)
        if ((x*x+y*y+z*z)<rcc2) then
          EE2 = EE2 + pos_ip1(4)*pos(k,4)*real_ij(periodic_x(x),periodic_y(y),Periodic_z(z))
        end if
      end if
    end do
  end do
  EE2 = EE2 * pos_ip1(4)
  DeltaE = DeltaE + EE2 - EE1

  !
  ! Modified term of slab geometry
  dMz = pos_ip0(4) * (pos_ip1(3) - pos_ip0(3))/2.D0         !sigma unit

  DeltaE = DeltaE + 2*pi / (Lx*Ly*Lz*Z_empty) * lb/Beta * (2*Mz*dMz + dMz*dMz)

  Mz = Mz + dMz
  !
  ! External field energy
  DeltaE=DeltaE-EF*pos_ip0(4)*(pos_ip1(3)-pos_ip0(3))/2.D0   !simga unit

end subroutine Delta_Energy


subroutine Delta_Energy_add(DeltaE)
  use global_variables
  implicit none
  real*8, intent(out) :: DeltaE
  integer :: i,j,k,x,y,z,x1,y1,z1,t,qq1,qq2,icelx,icely,icelz,ncel
  real*8 :: EE1,EE2,rr,dMz
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
    x1 = pos_ip0(1) - pos(k,1)
    y1 = pos_ip0(2) - pos(k,2)
    z1 = pos_ip0(3) - pos(k,3)     
    x1 = fourier_x(x1)
    y1 = fourier_y(y1)
    z1 = fourier_z(z1)
    t = x1 + y1 + z1
    EE1 = EE1 + pos(k,4)*fourier_ij(t)

    x1 = pos_ip1(1) - pos(k,1)
    y1 = pos_ip1(2) - pos(k,2)
    z1 = pos_ip1(3) - pos(k,3)
    x1 = fourier_x(x1)
    y1 = fourier_y(y1)
    z1 = fourier_z(z1)
    t = x1 + y1 + z1
    EE2 = EE2 + pos(k,4)*fourier_ij(t)
  end do

  DeltaE = EE1*qq1 + EE2*qq2

  !
  ! Real Space
  EE1 = 0
  icelx = int(pos_ip1(1)/clx)+1
  icely = int(pos_ip1(2)/cly)+1
  icelz = int(pos_ip1(3)/clz)+1
  ncel = (icelx-1)*ncly*nclz+(icely-1)*nclz+icelz
  do i = 1, cell_near_list(ncel,28,1)
    icelx = cell_near_list(ncel,i,1)
    icely = cell_near_list(ncel,i,2)
    icelz = cell_near_list(ncel,i,3)   
    j = hoc_r(icelx,icely,icelz) 
    do while (j/=0) 
      k = charge(j)
      x = pos_ip1(1) - pos(k,1)
      y = pos_ip1(2) - pos(k,2)
      z = pos_ip1(3) - pos(k,3)
      if ((x*x+y*y+z*z)<rcc2) then
        EE1 = EE1 + pos_ip1(4)*pos(k,4)*real_ij(periodic_x(x),periodic_y(y),Periodic_z(z))
      end if
    end do
  end do
  EE1 = EE1 * qq1

  EE2 = 0
  icelx = int(pos_ip1i(1)/clx)+1
  icely = int(pos_ip1i(2)/cly)+1
  icelz = int(pos_ip1i(3)/clz)+1
  ncel = (icelx-1)*ncly*nclz+(icely-1)*nclz+icelz
  do i = 1, cell_near_list(ncel,28,1)
    icelx = cell_near_list(ncel,i,1)
    icely = cell_near_list(ncel,i,2)
    icelz = cell_near_list(ncel,i,3)   
    j = hoc_r(icelx,icely,icelz) 
    do while (j/=0) 
      k = charge(j)
      x = pos_ip1i(1) - pos(k,1)
      y = pos_ip1i(2) - pos(k,2)
      z = pos_ip1i(3) - pos(k,3)
      if ((x*x+y*y+z*z)<rcc2) then
        EE2=EE2+pos_ip1i(4)*pos(k,4)*real_ij(periodic_x(x),periodic_y(y),Periodic_z(z))
      end if
    end do
  end do
  EE2 = EE2 * qq1
  DeltaE = DeltaE + EE1 + EE2

  !
  !interaction of the added two particles
  x = pos_ip1i(1) - pos_ip1(1)
  y = pos_ip1i(2) - pos_ip1(2)
  z = pos_ip1i(3) - pos_ip1(3)
  x1 = fourier_x(x)
  y1 = fourier_y(y)
  z1 = fourier_z(z)
  x = Periodic_x(x)
  y = Periodic_y(y)
  z = Periodic_z(z)
  t = x1 + y1 + z1
  DeltaE = DeltaE + qq1*qq2*(fourier_ij(t)+real_ij(x,y,z))

  !
  !modified term of slab geometry
  dMz = pos_ip1i(4)*pos_ip1i(3) + pos_ip1(4)*pos_ip1(3)
  DeltaE = DeltaE + 2*pi / (Lx*Ly*Lz*Z_empty) * lb/Beta * (2*Mz*dMz + dMz*dMz)

  !
  !External field energy
  DeltaE = DeltaE - EF * pos_ip1i(4)*pos_ip1i(3) - EF * pos_ip1(4)*pos_ip1(3)

end subroutine Delta_Energy_add


subroutine Delta_Energy_delete(DeltaE)
  use global_variables
  implicit none
  real*8, intent(out) :: DeltaE
  integer :: i,j,k,x,y,z,x1,y1,z1,t,qq1,qq2,icelx,icely,icelz,ncel
  real*8 :: EE1,EE2,rr,dMz
  real*8, dimension(3) :: rij

  DeltaE = 0

  !
  ! Fourier Space
  EE1 = 0
  EE2 = 0
  qq1 = pos_ip0(4)
  qq2 = pos_ip0i(4)
  j = cell_list_q(Nq+1)
  do while (j/=0)
    k = charge(j)
    if (k/=ip) then
      x1 = pos_ip0(1) - pos(k,1)
      y1 = pos_ip0(2) - pos(k,2)
      z1 = pos_ip0(3) - pos(k,3)     
      x1 = fourier_x(x1)
      y1 = fourier_y(y1)
      z1 = fourier_z(z1)
      t = x1 + y1 + z1
      EE1 = EE1 + pos(k,4)*fourier_ij(t)
    end if
    if (k/=ip1) then
      x1 = pos_ip1(1) - pos(k,1)
      y1 = pos_ip1(2) - pos(k,2)
      z1 = pos_ip1(3) - pos(k,3)
      x1 = fourier_x(x1)
      y1 = fourier_y(y1)
      z1 = fourier_z(z1)
      t = x1 + y1 + z1
      EE2 = EE2 + pos(k,4)*fourier_ij(t)
    end if
  end do

  DeltaE = EE1*qq1 + EE2*qq2

  !
  ! Real Space
  EE1 = 0
  icelx = int(pos_ip1(1)/clx)+1
  icely = int(pos_ip1(2)/cly)+1
  icelz = int(pos_ip1(3)/clz)+1
  ncel = (icelx-1)*ncly*nclz+(icely-1)*nclz+icelz
  do i = 1, cell_near_list(ncel,28,1)
    icelx = cell_near_list(ncel,i,1)
    icely = cell_near_list(ncel,i,2)
    icelz = cell_near_list(ncel,i,3)   
    j = hoc_r(icelx,icely,icelz) 
    do while (j/=0) 
      k = charge(j)
      if (k/=ip) then
        x = pos_ip1(1) - pos(k,1)
        y = pos_ip1(2) - pos(k,2)
        z = pos_ip1(3) - pos(k,3)
        if ((x*x+y*y+z*z)<rcc2) then
          EE1 = EE1 + pos_ip1(4)*pos(k,4)*real_ij(periodic_x(x),periodic_y(y),Periodic_z(z))
        end if
      end if
    end do
  end do
  EE1 = EE1 * qq1

  EE2 = 0
  icelx = int(pos_ip1i(1)/clx)+1
  icely = int(pos_ip1i(2)/cly)+1
  icelz = int(pos_ip1i(3)/clz)+1 
  ncel = (icelx-1)*ncly*nclz+(icely-1)*nclz+icelz
  do i = 1, cell_near_list(ncel,28,1)
    icelx = cell_near_list(ncel,i,1)
    icely = cell_near_list(ncel,i,2)
    icelz = cell_near_list(ncel,i,3)   
    j = hoc_r(icelx,icely,icelz) 
    do while (j/=0) 
      k = charge(j)
      if (k/=ip1) then
        x = pos_ip1i(1) - pos(k,1)
        y = pos_ip1i(2) - pos(k,2)
        z = pos_ip1i(3) - pos(k,3)
        if ((x*x+y*y+z*z)<rcc2) then
          EE2 = EE2 + pos_ip1i(3)*pos(k,4)*real_ij(periodic_x(x),periodic_y(y),Periodic_z(z))
        end if
      end if
    end do
  end do
  EE2 = EE2 * qq1
  DeltaE = DeltaE + EE1 + EE2
  !
  !interaction of the added two particles
  x = pos_ip1i(1) - pos_ip1(1)
  y = pos_ip1i(2) - pos_ip1(2)
  z = pos_ip1i(3) - pos_ip1(3)
  x1 = fourier_x(x)
  y1 = fourier_y(y)
  z1 = fourier_z(z)
  x = Periodic_x(x)
  y = Periodic_y(y)
  z = Periodic_z(z)
  t = x1 + y1 + z1
  DeltaE = DeltaE - qq1*qq2*(fourier_ij(t)+real_ij(x,y,z))
  !
  !modified term of slab geometry
  dMz = pos_ip0i(4)*pos_ip0i(3) + pos_ip0(4)*pos_ip0(3)
  DeltaE = DeltaE + 2*pi / (Lx*Ly*Lz*Z_empty) * lb/Beta * (2*Mz*dMz + dMz*dMz)
  !
  !External field energy
  DeltaE = DeltaE - EF * pos_ip1i(4)*pos_ip1i(3) - EF * pos_ip1(4)*pos_ip1(3)

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

  icelx = int(pos_ip1(1)/clx)+1
  icely = int(pos_ip1(2)/cly)+1
  icelz = int(pos_ip1(3)/clz)+1

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

  icelx = int(pos_ip0(1)/clx)+1
  icely = int(pos_ip0(2)/cly)+1
  icelz = int(pos_ip0(3)/clz)+1

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


