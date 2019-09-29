module compute_energy_ewald
  !--------------------------------------!
  ! 
  !--------------------------------------!
  implicit none

  save

!############coefficient in potential function#############!
  !
  !coulomb
  !
  !coulomb potential
  real*8,  private :: lb          !Bjerrum length
  real*8,  private :: EF          !electric field
  real*8,  private :: tol         !tolerance
  real*8,  private :: tau_rf      !time ratio of real space and fourier space
  real*8,  private :: alpha       !Ewald screening parameter alpha
  real*8,  private :: alpha2      !alpha2=alpha*alpha
  real*8,  private :: Mz          !total dipole moment
  real*8,  private :: dMz         !delta Mz
  real*8,  private :: Mz_coef     !Coefficients in correction coulomb energy
  real*8,  private :: err
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
  integer, private :: K_total     !Total wave number in reciprocal space
!##########end coefficient in potential function###########!


!##########################arrays##########################!
  !
  !charge number to monomer number        
  integer, allocatable, dimension( : )          :: charge
  !  
  !monomer number to charge number
  integer, allocatable, dimension( : ), private :: inv_charge
  !
  !neighbor cells of the center cell
  integer, allocatable, dimension(:,:,:), private :: cell_near_list 
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
  !Coulomb energy of i,j in real space
  real,  allocatable, dimension(:,:,:), private :: real_ij 
  !
  !coefficients in Fourier space
  real*8,  allocatable, dimension( : ), private :: exp_ksqr
  !
  !structure factor
  complex(kind=8), allocatable, dimension( : ), private :: rho_k
  !
  !difference of structure factor
  complex(kind=8), allocatable, dimension( : ), private :: delta_rhok
  !
  !difference of structure factor
  complex(kind=8), allocatable, dimension( : ), private :: delta_cosk
  !
  !wave vector ordinal number
  integer, allocatable, dimension(:,:), private :: totk_vectk
  !
  !same as lj_point and real_point
  real*8,  allocatable, dimension(:,:), private :: posq
!########################end arrays########################!


contains


subroutine initialize_energy_parameters_Ewald
  !--------------------------------------!
  !Initial parameters are not inputted from file and compute
  !the total energy of the system.
  !Input
  !   
  !Output
  !   
  !External Variables
  !   
  !Routine Referenced:
  !1.
  !Reference:
  !The computation of alpha, rc_real et al are refered to
  !Frenkel, Smit, 'Understanding molecular simulation: from
  !algorithm to applications', Elsevier, 2002, pp.304-306.
  !--------------------------------------!
  use global_variables
  implicit none
  !
  !read energy parameters from file
  call read_energy_parameters_Ewald
  !
  !
  if ( qq /= 0 ) then
    !
    !Initialize ewald parameters and array allocate.
    call Initialize_ewald_parameters
    !
    !Construct the array totk_vectk(K_total,3), and allocate
    !rho_k(K_total), delta_rhok(K_total).
    call build_totk_vectk
    !
    !Construct the coefficients vector in Fourier space
    call build_exp_ksqr
    !
    !
    call Periodic_array
  end if

end subroutine initialize_energy_parameters_Ewald


subroutine initialize_energy_arrays_ewald
  use global_variables
  implicit none
  !
  !
  if ( qq /= 0 ) then
    !
    !Initialize charge with lind list. From this subroutine, pos array is needed.
    call Build_Charge_Ewald
    !
    !Initialize cell list of charge
    call Initialize_cell_list_q_Ewald
    !
    !Initialize real cell list
    call Initialize_real_cell_list_Ewald
    !
    !Construct the structure factor rho_k
    call build_rho_k
  end if

  call write_energy_parameters_Ewald

end subroutine initialize_energy_arrays_ewald


subroutine Periodic_array
  use global_variables
  implicit none
  integer :: i

  allocate(periodic_x(-Lx2:Lx2))
  allocate(periodic_y(-Ly2:Ly2))
  allocate(periodic_z(-Lz2:Lz2))
  Periodic_x = 0
  Periodic_y = 0
  periodic_z = 0

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


  open(120,file='./data/periodic_xy.txt')
    do i = -Lx2, Lx2
      write(120,*) i,periodic_x(i),Periodic_y(i)
    end do
  close(120)

  open(121,file='./data/periodic_z.txt')
    do i = -Lz2, Lz2
      write(121,*) i,periodic_z(i)
    end do
  close(121)

end subroutine Periodic_array


subroutine total_energy_ewald(EE, rt, ft)
  !--------------------------------------!
  !
  !   
  !Input
  !   
  !Output
  !   
  !External Variables
  !   
  !Routine Referenced:
  !1. 
  !--------------------------------------!
  use global_variables
  implicit none
  real*8, intent(out) :: EE
  real*8, intent(out) :: rt
  real*8, intent(out) :: ft
  integer :: i, j, k, l, m, n, x1, y1, z1, x, y, z, t
  integer :: icelx, icely, icelz, ncel
  real*8 :: st, fn, EE1, EE2,rr(4),q_total,Ec

  EE = 0
  Ec = 0
  !
  !real space
  call cpu_time(st)
  do i = 1, Nq
    m = charge(i)
    EE1 = 0
    icelx = int((pos(m,1)-1)/clx)+1
    icely = int((pos(m,2)-1)/cly)+1
    icelz = int((pos(m,3)-1)/clz)+1
    ncel = (icelx-1)*ncly*nclz+(icely-1)*nclz+icelz
    do j = 1, cell_near_list(ncel,28,1)
      icelx = cell_near_list(ncel,j,1)
      icely = cell_near_list(ncel,j,2)
      icelz = cell_near_list(ncel,j,3)
      k = hoc_r(icelx,icely,icelz)
      do while(k/=0)
        l = charge(k)
        if (l/=m) then
          x = pos(m,1)-pos(l,1)
          y = pos(m,2)-pos(l,2)
          z = pos(m,3)-pos(l,3)
          x = periodic_x(x)
          y = periodic_y(y)
          z = Periodic_z(z)
          if ((x*x+y*y+z*z)<rcc2) then
            EE1=EE1+pos(l,4)*real_ij(x,y,z)
          end if
        end if
        k = cell_list_r(k)
      end do
    end do
    EE = EE + EE1*pos(m,4)/2.D0 
    !
    !external field
    EE = EE - EF*pos(m,4)*pos(m,3)/2.D0
    !
    !self corrected term
    Ec = Ec - sqrt(alpha2/pi) * posq(i,4) * posq(i,4)
  end do
  call cpu_time(fn)
  rt = fn - st
  !
  !fourier space
  EE = EE + Ec/Beta*lb + sum( exp_ksqr * real( conjg(rho_k) * rho_k ) )/2.D0
  !
  !modified term of slab geometry
  EE = EE + Mz_coef * Mz**2
  
end subroutine total_energy_ewald


subroutine Delta_Energy_Ewald(DeltaE)
  !--------------------------------------!
  !Compute change of energy.
  !   
  !Input
  !   
  !Output
  !   DeltaE
  !External Variables
  !   pos_ip0, pos_ip1, ip
  !   inv_charge, DeltaE, EF
  !Routine Referenced:
  !1.Delta_LJ_Energy(DeltaE)
  !2.Delta_FENE_Energy(DeltaE)
  !3.Delta_real_Energy(DeltaE)
  !4.Delta_Reciprocal_Energy(DeltaE)
  !--------------------------------------!
  use global_variables
  implicit none
	real*8,  intent(out) :: DeltaE

  DeltaE = 0
  !
  !Compute Coulomb energy
  !
  !Compute coulomb energy change in real space
  call Delta_real_Energy(DeltaE)
  !
  !Compute coulomb energy change in reciprocal space
  call Delta_Reciprocal_Energy(DeltaE)

end subroutine Delta_Energy_Ewald


subroutine Delta_Energy_Ewald_add(DeltaE)
  !--------------------------------------!
  !Compute change of energy.
  !   
  !Input
  !   
  !Output
  !   DeltaE
  !External Variables
  !   pos_ip0, pos_ip1, ip
  !   inv_charge, DeltaE, EF
  !Routine Referenced:
  !1.Delta_LJ_Energy(DeltaE)
  !2.Delta_FENE_Energy(DeltaE)
  !3.Delta_real_Energy(DeltaE)
  !4.Delta_Reciprocal_Energy(DeltaE)
  !--------------------------------------!
  use global_variables
  implicit none
  real*8,  intent(out) :: DeltaE

  DeltaE = 0
  !
  !Compute Coulomb energy
  !
  !Compute coulomb energy change in real space
  call Delta_real_Energy_add(DeltaE)
  !
  !Compute coulomb energy change in reciprocal space
  call Delta_Reciprocal_Energy_add(DeltaE)

end subroutine Delta_Energy_Ewald_add


subroutine Delta_Energy_Ewald_delete(DeltaE)
  !--------------------------------------!
  !Compute change of energy.
  !   
  !Input
  !   
  !Output
  !   DeltaE
  !External Variables
  !   pos_ip0, pos_ip1, ip
  !   inv_charge, DeltaE, EF
  !Routine Referenced:
  !1.Delta_LJ_Energy(DeltaE)
  !2.Delta_FENE_Energy(DeltaE)
  !3.Delta_real_Energy(DeltaE)
  !4.Delta_Reciprocal_Energy(DeltaE)
  !--------------------------------------!
  use global_variables
  implicit none
  real*8,  intent(out) :: DeltaE

  DeltaE = 0
  !
  !Compute coulomb energy change in real space
  call Delta_real_Energy_delete(DeltaE)
  !
  !Compute coulomb energy change in reciprocal space
  call Delta_Reciprocal_Energy_delete(DeltaE)

end subroutine Delta_Energy_Ewald_delete


subroutine Delta_real_Energy(DeltaE)
  !--------------------------------------!
  !
  !   
  !Input
  !   DeltaE
  !Output
  !   DeltaE
  !External Variables
  !   real_pair_list, real_point
  !   pos_ip0, pos_ip1, 
  !   Lx, Ly, alpha, lb, Beta
  !Reference:
  !1.The correction of Coulomb energy in slab geometry is
  !  referred to:
  !  In-Chul Yeh, Max L. Berkowitz, 'Ewald summation for 
  !  systems with slab geometry', Journal of Chemical 
  !  Physics, vol 111, number 7, (1999).
  !--------------------------------------!
  use global_variables
  implicit none
  real*8, intent(inout) :: DeltaE
  real*8  :: rij(3)
  real*8  :: EE1, EE2
  integer :: i,j,k,x,y,z,x1,y1,z1,t,icelx,icely,icelz,ncel

  EE1 = 0
  icelx = int((pos_ip0(1)-1)/clx)+1
  icely = int((pos_ip0(2)-1)/cly)+1
  icelz = int((pos_ip0(3)-1)/clz)+1
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
        x = periodic_x(x)
        y = periodic_y(y)
        z = Periodic_z(z)
        if ((x*x+y*y+z*z)<rcc2) then
          EE1=EE1+pos(k,4)*real_ij(x,y,z)
        end if
      end if
      j = cell_list_r(j)
    end do
  end do
  EE1 = EE1 * pos_ip1(4)

  EE2 = 0
  icelx = int((pos_ip1(1)-1)/clx)+1
  icely = int((pos_ip1(2)-1)/cly)+1
  icelz = int((pos_ip1(3)-1)/clz)+1
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
        x = periodic_x(x)
        y = periodic_y(y)
        z = Periodic_z(z)
        if ((x*x+y*y+z*z)<rcc2) then
          EE2=EE2+pos(k,4)*real_ij(x,y,z)
        end if
      end if
      j = cell_list_r(j)
    end do
  end do
  EE2 = EE2 * pos_ip1(4)
  DeltaE = DeltaE + EE2 - EE1

  !
  !Change of correction energy in slab geometry
  dMz = pos_ip0(4) * (pos_ip1(3) - pos_ip0(3))/2.D0
  DeltaE = DeltaE + Mz_coef * (2*Mz*dMz + dMz*dMz)
  !
  ! External field energy
  DeltaE=DeltaE-EF*pos_ip0(4)*(pos_ip1(3)-pos_ip0(3))/2.D0   !simga unit

end subroutine Delta_real_Energy


subroutine Delta_real_Energy_add(DeltaE)
  use global_variables
  implicit none
  real*8, intent(out) :: DeltaE
  integer :: i,j,k,x,y,z,x1,y1,z1,t,qq1,qq2,icelx,icely,icelz,ncel
  real*8 :: EE1,EE2
  real*8, dimension(3) :: rij

  DeltaE = 0
  EE1 = 0
  icelx = int((pos_ip1(1)-1)/clx)+1
  icely = int((pos_ip1(2)-1)/cly)+1
  icelz = int((pos_ip1(3)-1)/clz)+1
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
      x = periodic_x(x)
      y = periodic_y(y)
      z = Periodic_z(z)
      if ((x*x+y*y+z*z)<rcc2) then
        EE1=EE1+pos(k,4)*real_ij(x,y,z)
      end if
      j = cell_list_r(j)
    end do
  end do
  EE1 = EE1 * qq1

  EE2 = 0
  icelx = int((pos_ip1i(1)-1)/clx)+1
  icely = int((pos_ip1i(2)-1)/cly)+1
  icelz = int((pos_ip1i(3)-1)/clz)+1
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
      x = periodic_x(x)
      y = periodic_y(y)
      z = Periodic_z(z)
      if ((x*x+y*y+z*z)<rcc2) then
        EE2=EE2+pos(k,4)*real_ij(x,y,z)
      end if
      j = cell_list_r(j)
    end do
  end do
  EE2 = EE2 * qq2
  DeltaE = DeltaE + EE1 + EE2

  !
  !interaction of the added two particles
  x = pos_ip1i(1) - pos_ip1(1)
  y = pos_ip1i(2) - pos_ip1(2)
  z = pos_ip1i(3) - pos_ip1(3)
  x = Periodic_x(x)
  y = Periodic_y(y)
  z = Periodic_z(z)
  if (x*x+y*y+z*z<rcc2) then
    DeltaE = DeltaE + qq1*qq2*real_ij(x,y,z)
  end if
  !
  !modified term of slab geometry
  dMz = pos_ip1i(4)*pos_ip1i(3)/2.D0 + pos_ip1(4)*pos_ip1(3)/2.D0
  DeltaE = DeltaE + 2*pi / (Lx*Ly*Lz*Z_empty) * lb/Beta * (2*Mz*dMz + dMz*dMz)
  !
  !External field energy
  DeltaE=DeltaE-EF*pos_ip1i(4)*pos_ip1i(3)/2.D0-EF*pos_ip1(4)*pos_ip1(3)/2.D0
  !
  !
  DeltaE = DeltaE - sqrt(alpha2/pi)*(pos_ip1(4)**2+pos_ip1i(4)**2)

end subroutine Delta_real_Energy_add


subroutine Delta_real_Energy_delete(DeltaE)
  use global_variables
  implicit none
  real*8, intent(out) :: DeltaE
  integer :: i,j,k,x,y,z,x1,y1,z1,t,qq1,qq2,icelx,icely,icelz,ncel
  real*8 :: EE1,EE2
  real*8, dimension(3) :: rij

  DeltaE = 0
  !
  ! Real Space
  EE1 = 0
  icelx = int((pos_ip1(1)-1)/clx)+1
  icely = int((pos_ip1(2)-1)/cly)+1
  icelz = int((pos_ip1(3)-1)/clz)+1
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
        x = periodic_x(x)
        y = periodic_y(y)
        z = Periodic_z(z)
        if ((x*x+y*y+z*z)<rcc2) then
          EE1=EE1+pos(k,4)*real_ij(x,y,z)
        end if
      end if
      j = cell_list_r(j)
    end do
  end do
  EE1 = -EE1 * qq1

  EE2 = 0
  icelx = int((pos_ip1i(1)-1)/clx)+1
  icely = int((pos_ip1i(2)-1)/cly)+1
  icelz = int((pos_ip1i(3)-1)/clz)+1 
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
        x = periodic_x(x)
        y = periodic_y(y)
        z = Periodic_z(z)
        if ((x*x+y*y+z*z)<rcc2) then
          EE2=EE2+pos(k,4)*real_ij(x,y,z)
        end if
      end if
      j = cell_list_r(j)
    end do
  end do
  EE2 = -EE2 * qq2
  DeltaE = DeltaE + EE1 + EE2
  if (abs(DeltaE)>1e10) then
    write(*,*) 'real'
    stop
  end if
  !
  !interaction of the added two particles
  x = pos_ip1i(1) - pos_ip1(1)
  y = pos_ip1i(2) - pos_ip1(2)
  z = pos_ip1i(3) - pos_ip1(3)
  x = Periodic_x(x)
  y = Periodic_y(y)
  z = Periodic_z(z)
  if ((x*x+y*y+z*z)<rcc2) then
    DeltaE = DeltaE + qq1*qq2*real_ij(x,y,z)
  end if 
  !
  !modified term of slab geometry
  dMz = - pos_ip0i(4)*pos_ip0i(3)/2.D0 - pos_ip0(4)*pos_ip0(3)/2.D0
  DeltaE = DeltaE + 2*pi / (Lx*Ly*Lz*Z_empty) * lb/Beta * (2*Mz*dMz + dMz*dMz)
  !
  !External field energy
  DeltaE=DeltaE+EF * pos_ip0i(4)*pos_ip0i(3)/2.D0+EF*pos_ip0(4)*pos_ip0(3)/2.D0

  DeltaE = DeltaE + sqrt(alpha2/pi)*(qq1**2+qq2**2)

end subroutine Delta_real_Energy_delete


subroutine Delta_Reciprocal_Energy(DeltaE)  
  !------------------------------------!
  !This program is used to calculate the difference of the 
  !electrical energy in reciprocal space between new and old
  !position.
  !   
  !Input
  !   
  !Output
  !   Del_Recip_erg
  !External Variables
  !   exp_ksqr, rho_k, delta_rhok, totk_vectk
  !   pos_ip0, pos_ip1, ip
  !   Kmax1, Kmax2, Kmax3, Lx, Ly, Lz, Z_empty
  !Routine Referenced:
  !1.
  !------------------------------------!
  use global_variables
  implicit none
  real*8, intent(inout) :: DeltaE
  real*8  :: Del_Recip_erg
  complex(kind=8) :: eikx0( -Kmax1:Kmax1 ), eikx1( -Kmax1:Kmax1 )
  complex(kind=8) :: eiky0( -Kmax2:Kmax2 ), eiky1( -Kmax2:Kmax2 )
  complex(kind=8) :: eikz0( 0:Kmax3 ), eikz1( 0:Kmax3 )
  complex(kind=8) :: eikr0, eikr1
  real*8  :: c1, c2, c3
  integer :: ord(3), i, p, q, r

  c1 = 2*pi / Lx
  c2 = 2*pi / Ly
  c3 = 2*pi / (Lz*Z_empty)

  eikx0(0)  = ( 1,0 )
  eiky0(0)  = ( 1,0 )
  eikz0(0)  = ( 1,0 )
  eikx0(1)  = cmplx( cos(c1*pos_ip0(1)/2.D0), sin(-c1*pos_ip0(1)/2.D0), 8 )
  eiky0(1)  = cmplx( cos(c2*pos_ip0(2)/2.D0), sin(-c2*pos_ip0(2)/2.D0), 8 )
  eikz0(1)  = cmplx( cos(c3*pos_ip0(3)/2.D0), sin(-c3*pos_ip0(3)/2.D0), 8 )
  eikx0(-1) = conjg( eikx0(1) )
  eiky0(-1) = conjg( eiky0(1) )

  do p = 2, Kmax1
    eikx0(p)  = eikx0(p-1) * eikx0(1)
    eikx0(-p) = conjg( eikx0(p) )
  end do
  do q = 2, Kmax2
    eiky0(q)  = eiky0(q-1) * eiky0(1)
    eiky0(-q) = conjg(eiky0(q))
  end do
  do r = 2, Kmax3
    eikz0(r)  = eikz0(r-1) * eikz0(1)
  end do

  eikx1(0)  = ( 1,0 )
  eiky1(0)  = ( 1,0 )
  eikz1(0)  = ( 1,0 )
  eikx1(1)  = cmplx( cos(c1*pos_ip1(1)/2.D0), sin(-c1*pos_ip1(1)/2.D0), 8 )
  eiky1(1)  = cmplx( cos(c2*pos_ip1(2)/2.D0), sin(-c2*pos_ip1(2)/2.D0), 8 )
  eikz1(1)  = cmplx( cos(c3*pos_ip1(3)/2.D0), sin(-c3*pos_ip1(3)/2.D0), 8 )
  eikx1(-1) = conjg( eikx1(1) )
  eiky1(-1) = conjg( eiky1(1) )

  do p=2, Kmax1
    eikx1(p)  = eikx1(p-1) * eikx1(1)
    eikx1(-p) = conjg( eikx1(p) )
  end do
  do q=2, Kmax2
    eiky1(q)  = eiky1(q-1) * eiky1(1)
    eiky1(-q) = conjg(eiky1(q))
  end do
  do r=2, Kmax3
    eikz1(r)  = eikz1(r-1) * eikz1(1)
  end do

  do i=1, K_total
    ord = totk_vectk(i,1:3)
    eikr0 = eikx0(ord(1)) * eiky0(ord(2)) * eikz0(ord(3))
    eikr1 = eikx1(ord(1)) * eiky1(ord(2)) * eikz1(ord(3))
    delta_rhok(i) = eikr1 - eikr0
    delta_cosk(i) = 1 - real( conjg(eikr1) * eikr0 )
  end do

  delta_rhok = delta_rhok * pos_ip0(4)

  delta_cosk = delta_cosk * ( pos_ip0(4) * pos_ip0(4) )

  Del_Recip_erg = sum( exp_ksqr * ( Real( rho_k * delta_rhok ) + delta_cosk ) )

  DeltaE = DeltaE + Del_Recip_erg

end subroutine Delta_Reciprocal_Energy


subroutine Delta_Reciprocal_Energy_add(DeltaE)
  use global_variables
  implicit none
  real*8, intent(inout) :: DeltaE
  real*8  :: Del_Recip_erg
  complex(kind=8) :: eikx0( -Kmax1:Kmax1 ), eikx1( -Kmax1:Kmax1 )
  complex(kind=8) :: eiky0( -Kmax2:Kmax2 ), eiky1( -Kmax2:Kmax2 )
  complex(kind=8) :: eikz0( 0:Kmax3 ), eikz1( 0:Kmax3 )
  complex(kind=8) :: eikr0, eikr1
  real*8  :: c1, c2, c3
  integer :: ord(3), i, p, q, r

  c1 = 2*pi / Lx
  c2 = 2*pi / Ly
  c3 = 2*pi / (Lz*Z_empty)

  eikx0(0)  = ( 1,0 )
  eiky0(0)  = ( 1,0 )
  eikz0(0)  = ( 1,0 )
  eikx0(1)  = cmplx( cos(c1*pos_ip1(1)/2.D0), sin(-c1*pos_ip1(1)/2.D0), 8 )
  eiky0(1)  = cmplx( cos(c2*pos_ip1(2)/2.D0), sin(-c2*pos_ip1(2)/2.D0), 8 )
  eikz0(1)  = cmplx( cos(c3*pos_ip1(3)/2.D0), sin(-c3*pos_ip1(3)/2.D0), 8 )
  eikx0(-1) = conjg( eikx0(1) )
  eiky0(-1) = conjg( eiky0(1) )

  do p = 2, Kmax1
    eikx0(p)  = eikx0(p-1) * eikx0(1)
    eikx0(-p) = conjg( eikx0(p) )
  end do
  do q = 2, Kmax2
    eiky0(q)  = eiky0(q-1) * eiky0(1)
    eiky0(-q) = conjg(eiky0(q))
  end do
  do r = 2, Kmax3
    eikz0(r)  = eikz0(r-1) * eikz0(1)
  end do

  eikx1(0)  = ( 1,0 )
  eiky1(0)  = ( 1,0 )
  eikz1(0)  = ( 1,0 )
  eikx1(1)  = cmplx( cos(c1*pos_ip1i(1)/2.D0), sin(-c1*pos_ip1i(1)/2.D0), 8 )
  eiky1(1)  = cmplx( cos(c2*pos_ip1i(2)/2.D0), sin(-c2*pos_ip1i(2)/2.D0), 8 )
  eikz1(1)  = cmplx( cos(c3*pos_ip1i(3)/2.D0), sin(-c3*pos_ip1i(3)/2.D0), 8 )
  eikx1(-1) = conjg( eikx1(1) )
  eiky1(-1) = conjg( eiky1(1) )

  do p=2, Kmax1
    eikx1(p)  = eikx1(p-1) * eikx1(1)
    eikx1(-p) = conjg( eikx1(p) )
  end do
  do q=2, Kmax2
    eiky1(q)  = eiky1(q-1) * eiky1(1)
    eiky1(-q) = conjg(eiky1(q))
  end do
  do r=2, Kmax3
    eikz1(r)  = eikz1(r-1) * eikz1(1)
  end do

  do i=1, K_total
    ord = totk_vectk(i,1:3)
    eikr0 = eikx0(ord(1)) * eiky0(ord(2)) * eikz0(ord(3))
    eikr1 = eikx1(ord(1)) * eiky1(ord(2)) * eikz1(ord(3))
    delta_rhok(i) = eikr1 - eikr0
    delta_cosk(i) = 1 - real( conjg(eikr1) * eikr0 )
  end do

  Del_Recip_erg = sum( exp_ksqr * ( Real( rho_k * delta_rhok ) + delta_cosk ) )

  DeltaE = DeltaE + Del_Recip_erg

end subroutine Delta_Reciprocal_Energy_add


subroutine Delta_Reciprocal_Energy_delete(DeltaE)
  use global_variables
  real*8, intent(inout) :: DeltaE
  real*8  :: Del_Recip_erg
  complex(kind=8) :: eikx0( -Kmax1:Kmax1 ), eikx1( -Kmax1:Kmax1 )
  complex(kind=8) :: eiky0( -Kmax2:Kmax2 ), eiky1( -Kmax2:Kmax2 )
  complex(kind=8) :: eikz0( 0:Kmax3 ), eikz1( 0:Kmax3 )
  complex(kind=8) :: eikr0, eikr1
  real*8  :: c1, c2, c3
  integer :: ord(3), i, p, q, r

  c1 = 2*pi / Lx
  c2 = 2*pi / Ly
  c3 = 2*pi / (Lz*Z_empty)

  eikx0(0)  = ( 1,0 )
  eiky0(0)  = ( 1,0 )
  eikz0(0)  = ( 1,0 )
  eikx0(1)  = cmplx( cos(c1*pos_ip0(1)/2.D0), sin(-c1*pos_ip0(1)/2.D0), 8 )
  eiky0(1)  = cmplx( cos(c2*pos_ip0(2)/2.D0), sin(-c2*pos_ip0(2)/2.D0), 8 )
  eikz0(1)  = cmplx( cos(c3*pos_ip0(3)/2.D0), sin(-c3*pos_ip0(3)/2.D0), 8 )
  eikx0(-1) = conjg( eikx0(1) )
  eiky0(-1) = conjg( eiky0(1) )

  do p = 2, Kmax1
    eikx0(p)  = eikx0(p-1) * eikx0(1)
    eikx0(-p) = conjg( eikx0(p) )
  end do
  do q = 2, Kmax2
    eiky0(q)  = eiky0(q-1) * eiky0(1)
    eiky0(-q) = conjg(eiky0(q))
  end do
  do r = 2, Kmax3
    eikz0(r)  = eikz0(r-1) * eikz0(1)
  end do

  eikx1(0)  = ( 1,0 )
  eiky1(0)  = ( 1,0 )
  eikz1(0)  = ( 1,0 )
  eikx1(1)  = cmplx( cos(c1*pos_ip0i(1)/2.D0), sin(-c1*pos_ip0i(1)/2.D0), 8 )
  eiky1(1)  = cmplx( cos(c2*pos_ip0i(2)/2.D0), sin(-c2*pos_ip0i(2)/2.D0), 8 )
  eikz1(1)  = cmplx( cos(c3*pos_ip0i(3)/2.D0), sin(-c3*pos_ip0i(3)/2.D0), 8 )
  eikx1(-1) = conjg( eikx1(1) )
  eiky1(-1) = conjg( eiky1(1) )

  do p=2, Kmax1
    eikx1(p)  = eikx1(p-1) * eikx1(1)
    eikx1(-p) = conjg( eikx1(p) )
  end do
  do q=2, Kmax2
    eiky1(q)  = eiky1(q-1) * eiky1(1)
    eiky1(-q) = conjg(eiky1(q))
  end do
  do r=2, Kmax3
    eikz1(r)  = eikz1(r-1) * eikz1(1)
  end do

  do i=1, K_total
    ord = totk_vectk(i,1:3)
    eikr0 = eikx0(ord(1)) * eiky0(ord(2)) * eikz0(ord(3))
    eikr1 = eikx1(ord(1)) * eiky1(ord(2)) * eikz1(ord(3))
    delta_rhok(i) = - eikr1 + eikr0
    delta_cosk(i) = 1 - real( conjg(eikr1) * eikr0 )
  end do

  Del_Recip_erg = sum( exp_ksqr * ( Real( rho_k * delta_rhok ) + delta_cosk ) )

  DeltaE = DeltaE + Del_Recip_erg

end subroutine Delta_Reciprocal_Energy_delete


subroutine update_real_cell_list_Ewald
  use global_variables
  implicit none
  integer :: icelx,icely,icelz

  icelx = int((pos_ip0(1)-1)/clx)+1
  icely = int((pos_ip0(2)-1)/cly)+1
  icelz = int((pos_ip0(3)-1)/clz)+1
  call update_real_cell_list_delete_Ewald(ip,icelx,icely,icelz)

  icelx = int((pos_ip1(1)-1)/clx)+1
  icely = int((pos_ip1(2)-1)/cly)+1
  icelz = int((pos_ip1(3)-1)/clz)+1
  call update_real_cell_list_add_Ewald(ip,icelx,icely,icelz)

  Mz = Mz + dMz

end subroutine update_real_cell_list_Ewald


subroutine update_cell_list_pH_add_Ewald
  use global_variables
  implicit none
  integer :: icelx,icely,icelz

  icelx = int((pos_ip1(1)-1)/clx)+1
  icely = int((pos_ip1(2)-1)/cly)+1
  icelz = int((pos_ip1(3)-1)/clz)+1
  call update_real_cell_list_add_Ewald(ip,icelx,icely,icelz)
  call update_charge_cell_list_add_Ewald(ip)

  icelx = int((pos_ip1i(1)-1)/clx)+1
  icely = int((pos_ip1i(2)-1)/cly)+1
  icelz = int((pos_ip1i(3)-1)/clz)+1
  call update_real_cell_list_add_Ewald(ip1,icelx,icely,icelz)
  call update_charge_cell_list_add_Ewald(ip1)

  Mz = Mz + dMz

  Nq_net = Nq_net + 1

end subroutine update_cell_list_pH_add_Ewald


subroutine update_cell_list_pH_delete_Ewald
  use global_variables
  implicit none
  integer :: icelx,icely,icelz

  icelx = int((pos_ip0(1)-1)/clx)+1
  icely = int((pos_ip0(2)-1)/cly)+1
  icelz = int((pos_ip0(3)-1)/clz)+1
  call update_real_cell_list_delete_Ewald(ip,icelx,icely,icelz)
  call update_charge_cell_list_delete_Ewald(ip)

  icelx = int((pos_ip0i(1)-1)/clx)+1
  icely = int((pos_ip0i(2)-1)/cly)+1
  icelz = int((pos_ip0i(3)-1)/clz)+1 
  call update_real_cell_list_delete_Ewald(ip1,icelx,icely,icelz)
  call update_charge_cell_list_delete_Ewald(ip1)

  Mz = Mz + dMz

  Nq_net = Nq_net - 1

end subroutine update_cell_list_pH_delete_Ewald


subroutine update_real_cell_list_add_Ewald(iq,icelx,icely,icelz)
  use global_variables
  implicit none
  integer, intent(in) :: iq
  integer, intent(in) :: icelx, icely, icelz
  integer :: nti,bfi,ii,j       ! next particle of ii, before particle of ii
  integer :: ed, st 

  ii = inv_charge(iq)   !ii belongs to [1,Nq]

  inv_cell_list_r(ii) = 0
  if ( inv_hoc_r(icelx,icely,icelz) /=0 ) then
    inv_cell_list_r( hoc_r(icelx,icely,icelz) ) = ii
  else
    inv_hoc_r(icelx,icely,icelz) = ii
  end if

  cell_list_r(ii) = hoc_r(icelx,icely,icelz)
  hoc_r(icelx,icely,icelz) = ii

!   j = hoc_r(icelx,icely,icelz)  
!   do while (j/=0)
!     j = cell_list_r(j)
!     if ((j==cell_list_r(j) .and. j/=0).or. j>20000) then
!       write(*,*) j
!       write(*,*) 'add'
!       stop
!     end if
!   end do

end subroutine update_real_cell_list_add_Ewald


subroutine update_real_cell_list_delete_Ewald(iq,icelx,icely,icelz)
  use global_variables
  implicit none
  integer, intent(in) :: iq
  integer, intent(in) :: icelx, icely, icelz
  integer :: nti,bfi,ii,j       ! next particle of ii, before particle of ii
  integer :: ed, st 

  ii = inv_charge(iq)   !ii belongs to [1,Nq]

  bfi = cell_list_r(ii)
  nti = inv_cell_list_r(ii)

  if ( bfi/=0 .and. nti/=0 ) then        !middle
    cell_list_r(nti) = bfi
    inv_cell_list_r(bfi) = nti
  elseif ( bfi==0 .and. nti/=0 ) then    !the first one
    cell_list_r(nti) = bfi
    inv_hoc_r(icelx,icely,icelz) = nti
  elseif ( bfi/=0 .and. nti==0 ) then    !the last one
    hoc_r(icelx,icely,icelz) = bfi
    inv_cell_list_r(bfi) = nti
  else                                   !only one
    hoc_r(icelx,icely,icelz) = nti
    inv_hoc_r(icelx,icely,icelz) = bfi
  end if

!   j = hoc_r(icelx,icely,icelz) 
!   do while (j/=0)
!     j = cell_list_r(j)
!     if ((j==cell_list_r(j) .and. j/=0).or.j>20000) then
!       write(*,*) 'delete'
!       write(*,*) j,bfi,nti,ii,icelx,icely,icelz
!       stop
!     end if
!   end do

end subroutine update_real_cell_list_delete_Ewald


subroutine update_charge_cell_list_add_Ewald(iq)
  use global_variables
  implicit none
  integer, intent(in) :: iq
  integer :: ii, j     

  ii = inv_charge(iq)         ! ii belongs to [1,Nq]

  inv_cell_list_q(ii) = 0
  if ( cell_list_q(Nq+1)/=0 ) then
    inv_cell_list_q(cell_list_q(Nq+1)) = ii
  else
    inv_cell_list_q(Nq+1) = ii
  end if

  cell_list_q(ii) = cell_list_q(Nq+1)
  cell_list_q(Nq+1) = ii

!   j = cell_list_q(Nq+1) 
!   do while (j/=0)
!     j = cell_list_q(j)
!     if ((j==cell_list_q(j) .and. j/=0) .or. j>20000) then
!       write(*,*) 'add pH'
!       write(*,*) j
!       stop
!     end if
!   end do

end subroutine update_charge_cell_list_add_Ewald


subroutine update_charge_cell_list_delete_Ewald(iq)
  use global_variables
  implicit none
  integer, intent(in) :: iq
  integer :: nti,bfi,ii,j       ! next particle of ii, before particle of ii

  ii = inv_charge(iq)         !ii belongs to [1,Nq]
  
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
    cell_list_q(Nq+1) = nti
    inv_cell_list_q(Nq+1) = bfi
  end if

!   j = cell_list_q(Nq+1) 
!   do while (j/=0)
!     j = cell_list_q(j)
!     if ((j==cell_list_q(j) .and. j/=0).or.j>20000) then
!       write(*,*) 'delete pH'
!       write(*,*) j,ii
!       stop
!     end if
!   end do

end subroutine update_charge_cell_list_delete_Ewald


subroutine error_analysis_ewald(EE1)
  !--------------------------------------!
  !
  !Reference:
  !1.The error is estimated by Kolafa1992,
  !  JIRI KOLAFA, JOHN W. PERRAM, 'CUTOFF ERRORS IN THE 
  !  EWALD SUMMATION FORMULAE FOR POINT CHARGE SYSTEMS',
  !  Molecular Simulation, Vol. 9(5), pp.351-368, (1992).
  !2.The true error e.g. RMS energy is defined at:
  !  Ulrich Essmann, Lalith Perera et al, 'A smooth particle
  !  mesh Ewald method', Journal of Chemical Physics, Vol. 103(9),
  !  (1995).
  !--------------------------------------!
  use global_variables
  implicit none
  real*8 :: e_r, e_k, q_tot, rmse, tol1, sumf1, sumf2, tm1, tm2
  real*8 :: EE0, EE2
  real*8, intent(out) :: EE1
  integer i,j,m

!   tol = 5                

  call total_energy_ewald(EE1,tm1,tm2)

end subroutine error_analysis_ewald


subroutine read_energy_parameters_Ewald
  !--------------------------------------!
  !
  !--------------------------------------!
  use global_variables
  implicit none

  open(unit=100, file='energy_data.txt')
    read(100,*) lb
    read(100,*) Ef
    read(100,*) tol
    read(100,*) tau_rf
  close(100)

end subroutine read_energy_parameters_Ewald


subroutine write_energy_parameters_Ewald
  !--------------------------------------!
  !
  !--------------------------------------!
  use global_variables
  implicit none

  write(*,*) '******************  Potential  *********************'
  write(*,*) 'Kmax1      :', Kmax1
  write(*,*) 'Kmax2      :', Kmax2
  write(*,*) 'Kmax3      :', Kmax3
  write(*,*) 'K_total    :', K_total
  write(*,*) 'alpha      :', alpha
  write(*,*) 'tol        :', tol
  write(*,*) 'tau_rf     :', tau_rf
  write(*,*) 'rcc        :', rcc
  write(*,*) 'nclx       :', nclx
  write(*,*) 'ncly       :', ncly
  write(*,*) 'nclz       :', nclz
  write(*,*) 'clx        :', clx
  write(*,*) 'cly        :', cly
  write(*,*) 'clz        :', clz
  write(*,*) '****************************************************'

end subroutine write_energy_parameters_Ewald


subroutine Initialize_ewald_parameters
  !--------------------------------------!
  !
  !--------------------------------------!
  use global_variables
  implicit none
  real*8 :: rho, v_verlet

  !
  ! alpha and rcc
  alpha    = ( tau_rf * pi**3 * Nq / (Lx*Ly)**2/Lz/Z_empty ) ** (1.D0/6)
  alpha2   = alpha * alpha
  rcc  = tol / alpha * 2
  rcc2 = rcc * rcc

  !
  !use verlet list in real space
  Kmax1 = ceiling(tol*Lx*alpha/pi)
  Kmax2 = ceiling(tol*Ly*alpha/pi)
  Kmax3 = ceiling(tol*Lz*Z_empty*alpha/pi)

  !
  !Cell list parameters
  nclx = int(Lx2/(rcc+1))     !cell numbers in x direction
  ncly = int(Ly2/(rcc+1))
  nclz = int(Lz2/(rcc+1))
  clx = 1.D0*Lx2/nclx         !cell length    
  cly = 1.D0*Ly2/ncly
  clz = 1.D0*Lz2/nclz
  Mz_coef = 2*pi / (Lx*Ly*Lz*Z_empty) * lb/Beta

end subroutine Initialize_ewald_parameters


subroutine update_rhok
  !--------------------------------------!
  !
  !--------------------------------------!
  use global_variables
  implicit none

  rho_k = rho_k + Conjg( delta_rhok )

  Mz = Mz + dMz
  
end subroutine update_rhok


subroutine pre_calculate_real_space
  use global_variables
  implicit none
  integer :: i, j, k 
  integer :: nnl, nn_half
  real*8 :: rr, x, y, z

  nnl = nint(rcc+1)    ! cut off length, lattice unit
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


subroutine Build_Charge_Ewald
  !--------------------------------------!
  !Initialize and charge.
  !   
  !Input
  !   pos
  !Output
  !   charge, inv_charge, Mz
  !External Variables
  !   pos, charge, NN, Nq
  !--------------------------------------!
  use global_variables
  implicit none
  integer :: i, j

  Mz = 0
  if (allocated(charge)) deallocate(charge)
  allocate(charge(Nq))
  if (allocated(inv_charge)) deallocate(inv_charge)
  allocate(inv_charge(NN))

  j = 0
  do i = 1, NN
    if ( pos(i,4) /= 0 ) then
      j = j + 1
      charge(j) = i
      Mz = Mz + pos(i,4)*pos(i,3)/2.D0 !sigma unit
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

end subroutine Build_Charge_Ewald


subroutine Initialize_cell_list_q_Ewald
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

  open(112,file='./data/cell_list_q.txt')
    do i = 1, Nq + 1
      write(112,*) cell_list_q(i), inv_cell_list_q(i)
    end do
  close(112)

end subroutine Initialize_cell_list_q_Ewald


subroutine Initialize_real_cell_list_Ewald
  use global_variables
  implicit none
  integer :: i, j, k, l, m, n, p, q, r, x, y, z
  integer :: icelx, icely, icelz

  !
  ! maxium situation, (125,125,100), 6.2Mb RAM is needed.
  allocate(hoc_r(nclx,ncly,nclz))
  allocate(inv_hoc_r(nclx,ncly,nclz))
  hoc_r = 0
  inv_hoc_r = 0

  allocate(cell_list_r(Nq))
  allocate(inv_cell_list_r(Nq))
  cell_list_r = 0
  inv_cell_list_r = 0

  do i = 1, Nq
    j = charge(i)
    icelx = int((pos(j,1)-1)/clx) + 1
    icely = int((pos(j,2)-1)/cly) + 1
    icelz = int((pos(j,3)-1)/clz) + 1
    cell_list_r(i) = hoc_r(icelx,icely,icelz)
    hoc_r(icelx,icely,icelz) = i
  end do

  do i = Nq, 1, -1
    j = charge(i)
    icelx = int((pos(j,1)-1)/clx) + 1
    icely = int((pos(j,2)-1)/cly) + 1
    icelz = int((pos(j,3)-1)/clz) + 1
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

!   open(113,file='./data/cell_list_r.txt')
!     do i = 1, Nq
!       write(113,*) i, cell_list_r(i), inv_cell_list_r(i)
!     end do
!   close(113)


!   open(100,file='./data/hoc_r.txt')
!     do i = 1, nclx
!       do j = 1, ncly
!         do k = 1, nclz
!          write(100,*) i,j,k,hoc_r(i,j,k),inv_hoc_r(i,j,k)
!         end do
!       end do
!     end do
!   close(100)

end subroutine Initialize_real_cell_list_Ewald


subroutine build_totk_vectk
  !--------------------------------------!
  !exp_ksqr, rho_k, delta_rhok are all vectors with size of
  !K_total. For i = 1 to K_total, we often need to know 
  !corresponding wave number kx,ky,kz. This progam build a 
  !array totk_vectk(1:K_total,3) to store kx,ky,kz.
  !What's more, rho_k and delta_rhok are allocated here.
  !   
  !Input
  !   
  !Output
  !   totk_vectk
  !   K_total
  !External Variables
  !   K_total, Kmax1, Kmax2, Kmax3
  !   totk_vectk
  !Routine Referenced:
  !   
  !--------------------------------------!
  use global_variables
  implicit none
  integer :: i, j, k, l
  real*8  :: ksqr, k1, k2, k3, factor, kcut

  K_total=0
  do k = 0, Kmax3
    do i = -Kmax1, Kmax1
      do j = -Kmax2, Kmax2
        kcut = (1.D0*i/Kmax1) * (1.D0*i/Kmax1) &
             + (1.D0*j/Kmax2) * (1.D0*j/Kmax2) &
             + (1.D0*k/Kmax3) * (1.D0*k/Kmax3)
        if ( kcut>1 .or. kcut==0 ) cycle
        K_total = K_total + 1
      end do
    end do
  end do

  if ( allocated(totk_vectk) ) deallocate(totk_vectk)
  if ( allocated(rho_k)      ) deallocate(rho_k)
  if ( allocated(delta_rhok) ) deallocate(delta_rhok)
  if ( allocated(delta_cosk) ) deallocate(delta_cosk)
  allocate( totk_vectk( K_total, 3 ) )
  allocate( rho_k( K_total )         )
  allocate( delta_rhok( K_total )    )
  allocate( delta_cosk( K_total )    )
  totk_vectk = 0
  rho_k      = 0
  delta_rhok = 0
  delta_cosk = 0

  l=0
  do k = 0, Kmax3
    do i = -Kmax1, Kmax1
      do j = -Kmax2, Kmax2
        kcut = (1.D0*i/Kmax1) * (1.D0*i/Kmax1) &
             + (1.D0*j/Kmax2) * (1.D0*j/Kmax2) &
             + (1.D0*k/Kmax3) * (1.D0*k/Kmax3)
        if ( kcut>1 .or. kcut==0 ) cycle
        l = l + 1
        totk_vectk( l, 1 ) = i
        totk_vectk( l, 2 ) = j
        totk_vectk( l, 3 ) = k
      end do
    end do
  end do

end subroutine build_totk_vectk


subroutine build_exp_ksqr
  !--------------------------------------!
  !Reciprocal energy is divided to three parts: 
  !1.structrure factor is referred to rho_k.
  !2.difference of structure factor between new and old
  !position is referred to delta_rhok.
  !3.the other which includes exp(k^2/4/alpha) is referred 
  !to exp_ksqr.
  !This program is used to bulid the third part.
  !
  !Input
  !   
  !Output
  !   exp_ksqr
  !External Variables
  !   K_total
  !   Kmax1, Kmax2, Kmax3
  !   alpha2, lb
  !   Lx, Ly, Lz, Z_empty
  !   Beta
  !Reference:
  !Frenkel, Smit, 'Understanding molecular simulation: from
  !algorithm to applications', Elsevier, 2002, pp.300(12.1.25),
  !however his alpha is alpha^2 in this program.b 
  !--------------------------------------!
  use global_variables
  implicit none
  integer :: i, j, k, l, ord(3)
  real*8  :: ksqr, k1, k2, k3, factor

  if ( allocated(exp_ksqr) ) deallocate(exp_ksqr)
  allocate( exp_ksqr(K_total) )
  exp_ksqr = 0

  l = 0
  do i = 1, K_total
    ord = totk_vectk(i,:)
    if ( ord(3) == 0 ) then
      factor = 1
    else
      factor = 2
    end if
    k1   = 2*pi*ord(1) / Lx
    k2   = 2*pi*ord(2) / Ly
    k3   = 2*pi*ord(3) / Lz/Z_empty 
    ksqr = k1*k1 + k2*k2 + k3*k3 
    exp_ksqr(i) = factor * 4*pi / (Lx*Ly*Lz*Z_empty) *  &
                  exp(-ksqr/4/alpha2) / ksqr * lb / Beta     
  end do

end subroutine build_exp_ksqr


subroutine build_rho_k
  !--------------------------------------!
  !Calculate the structure factor array.
  !   
  !Input
  !   
  !Output
  !   rho_k
  !External Variables
  !   pos, charge
  !   Nq, Lx, Ly, Lz, Z_empty, K_total
  !Routine Referenced:
  !1.
  !--------------------------------------!
  use global_variables
  implicit none
  complex(kind=8) :: eikx(1:Nq, -Kmax1:Kmax1)
  complex(kind=8) :: eiky(1:Nq, -Kmax2:Kmax2)
  complex(kind=8) :: eikz(1:Nq, 0:Kmax3)
  integer i,j,l,m,n,p,q,r,ord(3)
  real*8 :: c1, c2, c3
  real*8 :: zq(Nq)
  rho_k = 0

  do m = 1, Nq
    n = charge(m)
    zq(m) = pos(n,4)
  end do

  c1 = 2*pi/Lx
  c2 = 2*pi/Ly
  c3 = 2*pi/Lz/Z_empty
  do m = 1, Nq
    i = charge(m)
    eikx(m,0)  = (1,0)
    eiky(m,0)  = (1,0)
    eikz(m,0)  = (1,0)

    eikx(m,1)  = cmplx( cos(c1*pos(i,1)/2.D0), sin(c1*pos(i,1)/2.D0), 8 )
    eiky(m,1)  = cmplx( cos(c2*pos(i,2)/2.D0), sin(c2*pos(i,2)/2.D0), 8 )
    eikz(m,1)  = cmplx( cos(c3*pos(i,3)/2.D0), sin(c3*pos(i,3)/2.D0), 8 )

    eikx(m,-1) = conjg(eikx(m,1))
    eiky(m,-1) = conjg(eiky(m,1))
  end do

  do p=2, Kmax1
    do m=1, Nq
      eikx(m,p)=eikx(m,p-1)*eikx(m,1)
      eikx(m,-p)=conjg(eikx(m,p))
    end do
  end do
  do q=2, Kmax2
    do m=1, Nq
      eiky(m,q)=eiky(m,q-1)*eiky(m,1)
      eiky(m,-q)=conjg(eiky(m,q))
    end do
  end do
  do r=2, Kmax3
    do m=1, Nq
      eikz(m,r)=eikz(m,r-1)*eikz(m,1)
    end do
  end do

  do i = 1, K_total
    ord = totk_vectk(i,:)
    do m = 1, Nq
      rho_k(i) = rho_k(i) + &
                 zq(m) * eikx(m,ord(1)) * eiky(m,ord(2)) * eikz(m,ord(3))
    end do
  end do

end subroutine build_rho_k

end module compute_energy_ewald


















