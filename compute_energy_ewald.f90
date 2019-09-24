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
!##########end coefficient in potential function###########!


!##########################arrays##########################!
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


subroutine initialize_energy_parameters
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
  call read_energy_parameters
  !
  !
  if ( qq /= 0 ) then
    !
    !Initialize ewald parameters and array allocate.
    call Initialize_ewald_parameters
    !
    !Construct the relation vector charge(Nq) of pos(NN,4)
    !and posq(Nq,4). pos(NN,4) are known.
    call build_charge
    !
    !Construct the real verlet list and real_point vector
    if ( real_verlet == 1 ) then
      call build_real_verlet_list
    end if
    !
    !Construct the array totk_vectk(K_total,3), and allocate
    !rho_k(K_total), delta_rhok(K_total).
    call build_totk_vectk
    !
    !Construct the coefficients vector in Fourier space
    call build_exp_ksqr
    !
    !Construct the structure factor rho_k
    call build_rho_k
  end if
  !
  !write energy parameters
  call write_energy_parameters

end subroutine initialize_energy_parameters


subroutine total_energy (EE)
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

  EE=0

  if ( qq /= 0 ) then
 
    call Coulomb_Energy(EE)

    call Electrical_Energy(EE)

  end if

end subroutine total_energy


subroutine Coulomb_energy ( EE )
  !--------------------------------------!
  !Compute Coulomb energy, including real space,
  !reciprocal space, self energy and slab energy.
  !   
  !Input
  !   EE
  !Output
  !   EE
  !External Variables
  !   real_point, real_pair_list, exp_ksqr, rho_k
  !   Nq, alpha, Beta, lb
  !Routine Referenced:
  !1. rij_and_rr
  !--------------------------------------!
  use global_variables
  implicit none
  real*8, intent(inout) :: EE
  integer :: i, j, k, l, m, n
  real*8  :: Ec, rij(3), rr2, rr, ord(3)

  Ec = 0
  Mz = 0
  if ( real_verlet == 1 ) then
    do n = 1, Nq
      i = charge(n)
      if ( n == 1 ) then
        k = 1
        l = real_point(1)
      else
        k = real_point(n-1) + 1
        l = real_point(n)
      end if
      do m = k, l
        j = real_pair_list(m)
        call rij_and_rr( rij, rr2, i, j )
        if ( rr2 < rc_real*rc_real ) then
          rr = sqrt( rr2 )
          !
          !Real space energy
          Ec = Ec + posq(i,4) * posq(j,4) * erfc(alpha*rr) / rr / 2 
        end if
      end do
      !z component of dipole moment
      Mz = Mz + posq(i,4)*posq(i,3)
      !
      !Self energy
      Ec = Ec - sqrt(alpha2/pi) * posq(i,4) * posq(i,4)
    end do
  else
    do m = 1, Nq-1
      i = charge(m)
      do n = m+1, Nq
        j = charge(n)
        call rij_and_rr( rij, rr2, i, j )
        if ( rr2 < rc_real*rc_real ) then
          rr = sqrt( rr2 )
          !
          !Real space energy
          Ec = Ec + posq(i,4) * posq(j,4) * erfc(alpha*rr) / rr
        end if
      end do
      !z component of dipole moment
      Mz = Mz + posq(i,4)*posq(i,3)
      !
      !Self energy
      Ec = Ec - sqrt(alpha2/pi) * posq(i,4) * posq(i,4)
    end do
    i  = charge(Nq)
    Mz = Mz + posq(i,4)*posq(i,3)
    !
    !Self energy
    Ec = Ec - sqrt(alpha2/pi) * posq(i,4) * posq(i,4)
  end if 
  !
  !Reciprocal space energy
  EE = EE + Ec / Beta * lb + sum( exp_ksqr * real( conjg(rho_k) * rho_k ) ) / 2
  !
  !Correction energy of slab energy
  Mz_coef = 2*pi / (Lx*Ly*Lz*Z_empty) * lb/Beta
  EE = EE + Mz_coef * Mz**2
  
end subroutine Coulomb_energy


subroutine Electrical_energy(EE)
  !--------------------------------------!
  !Compute electric energy of external field.
  !   
  !Input
  !   EE
  !Output
  !   EE
  !External Variables
  !   Nq
  !Routine Referenced:
  !1.
  !--------------------------------------!
  use global_variables
  implicit none
  real*8, intent(inout) :: EE
  integer :: i, j

  do j = 1, Nq

    !
    !EF is positive from bottom to top
    !The potential at the bottom plate is zero.
    i = charge(j)
    EE = EE - EF * posq(i,3) * posq(i,4) 

  end do

end subroutine Electrical_energy


subroutine Delta_Energy(DeltaE)
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
  !Compute energy of LJ potential
  call Delta_LJ_Energy(DeltaE)
  !
  !Compute Delta energy of FENE potential
  if ( ip <= Npe ) then
    call Delta_FENE_Energy(DeltaE)
  end if
  !
  !Compute Coulomb energy
  if ( pos_ip0(4) /= 0 ) then
    !
    !Compute coulomb energy change in real space
    call Delta_real_Energy(DeltaE)
    !
    !Compute coulomb energy change in reciprocal space
    call Delta_Reciprocal_Energy(DeltaE)
    !
    !Compute electric energy change of external electric field
    DeltaE = DeltaE - EF * pos_ip0(4) * ( pos_ip1(3) - pos_ip0(3) )
  end if 

end subroutine Delta_Energy


subroutine Delta_Energy_time( DeltaE, time )
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
  real*8, intent(out) :: DeltaE
  real*8, dimension(3), intent(inout) :: time
  real*8 :: st_lj, fn_lj, st_real, fn_real
  real*8 :: st_Fourier, fn_Fourier

  DeltaE = 0
  !
  !Compute energy of LJ potential
  call cpu_time(st_lj)
  call Delta_LJ_Energy(DeltaE)
  call cpu_time(fn_lj)
  time(1) = time(1) + fn_lj - st_lj
  !
  !Compute Delta energy of FENE potential
  if ( ip <= Npe ) then
    call Delta_FENE_Energy(DeltaE)
  end if
  !
  !Compute Coulomb energy
  if ( pos_ip0(4) /= 0 ) then
    !
    !Compute coulomb energy change in real space
    call cpu_time(st_real)
    call Delta_real_Energy(DeltaE)
    call cpu_time(fn_real)
    time(2) = time(2) + fn_real - st_real
    !
    !Compute coulomb energy change in reciprocal space
    call cpu_time(st_Fourier)
    call Delta_Reciprocal_Energy(DeltaE)
    call cpu_time(fn_Fourier)
    time(3) = time(3) + fn_Fourier - st_Fourier
    !
    !Compute electric energy change of external electric field
    DeltaE = DeltaE - EF * pos_ip0(4) * ( pos_ip1(3) - pos_ip0(3) )
  end if 

end subroutine Delta_Energy_time


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
  integer :: i, j, k, l, ipq
  real*8  :: rr, rij(3), h_lx, h_ly, nh_lx, nh_ly
  real*8  :: EE, dMz

  EE=0
  h_lx = Lx / 2
  h_ly = Ly / 2
  nh_lx = - h_lx
  nh_ly = - h_ly
  if ( real_verlet == 1 ) then
    ipq = inv_charge(ip)
    if ( ipq == 1 ) then
      k = 1
      l = real_point(1)
    else
      k = real_point(ipq-1) + 1
      l = real_point(ipq)
    end if
    do j = k, l
      i = real_pair_list(j)
      !
      !Energy of Coulomb potential of old configuration in real space
      !
      rij = posq(i,1:3) - pos_ip0(1:3)
      !
      !Periodic condition
      if ( rij(1) > h_lx ) then
        rij(1) = rij(1) - Lx
      elseif ( rij(1) < nh_lx ) then
        rij(1) = rij(1) + Lx
      end if
      if ( rij(2) > h_ly ) then
        rij(2) = rij(2) - Ly
      elseif ( rij(2) < nh_ly ) then
        rij(2) = rij(2) + Ly
      end if
      rr = rij(1) * rij(1) + rij(2) * rij(2) + rij(3) * rij(3)
      if ( rr < rc_real2 ) then
        rr = sqrt(rr)
        EE = EE - posq(i,4) * pos_ip0(4) * erfc(alpha * rr) / rr
      end if
      !
      !Energy of Coulomb potential of new configuration in real space
      !
      rij = posq(i,1:3) - pos_ip1(1:3)
      !
      !Periodic condition
      if ( rij(1) > h_lx ) then
        rij(1) = rij(1) - Lx
      elseif ( rij(1) < nh_lx ) then
        rij(1) = rij(1) + Lx
      end if
      if ( rij(2) > h_ly ) then
        rij(2) = rij(2) - Ly
      elseif ( rij(2) < nh_ly ) then
        rij(2) = rij(2) + Ly
      end if
      rr = rij(1) * rij(1) + rij(2) * rij(2) + rij(3) * rij(3)
      if ( rr < rc_real2 ) then
        rr = sqrt(rr)
        EE = EE + posq(i,4) * pos_ip0(4) * erfc(alpha * rr) / rr
      end if
    end do
  else 
    do j = 1, Nq
      i = charge(j)
      if ( i == ip ) cycle
      !
      !Energy of Coulomb potential of old configuration in real space
      !
      rij = posq(i,1:3) - pos_ip0(1:3)
      !
      !Periodic condition
      if ( rij(1) > h_lx ) then
        rij(1) = rij(1) - Lx
      elseif ( rij(1) < nh_lx ) then
        rij(1) = rij(1) + Lx
      end if
      if ( rij(2) > h_ly ) then
        rij(2) = rij(2) - Ly
      elseif ( rij(2) < nh_ly ) then
        rij(2) = rij(2) + Ly
      end if
      rr = sqrt( rij(1) * rij(1) + rij(2) * rij(2) + rij(3) * rij(3) )
      EE = EE - posq(i,4) * pos_ip0(4) * erfc(alpha * rr) / rr
      !
      !Energy of Coulomb potential of new configuration in real space
      !
      rij = posq(i,1:3) - pos_ip1(1:3)
      !
      !Periodic condition
      if ( rij(1) > h_lx ) then
        rij(1) = rij(1) - Lx
      elseif ( rij(1) < nh_lx ) then
        rij(1) = rij(1) + Lx
      end if
      if ( rij(2) > h_ly ) then
        rij(2) = rij(2) - Ly
      elseif ( rij(2) < nh_ly ) then
        rij(2) = rij(2) + Ly
      end if
      rr = sqrt( rij(1) * rij(1) + rij(2) * rij(2) + rij(3) * rij(3) )
      EE = EE + posq(i,4) * pos_ip0(4) * erfc(alpha * rr) / rr
    end do
  end if

  dMz = pos_ip0(4) * (pos_ip1(3) - pos_ip0(3))
  !
  !Change of correction energy in slab geometry
  DeltaE = DeltaE +  lB/Beta * EE + Mz_coef * (2*Mz*dMz + dMz*dMz)

  Mz = Mz + dMz

end subroutine Delta_real_Energy


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
  eikx0(1)  = cmplx( cos( c1 * pos_ip0(1) ), sin( -c1 * pos_ip0(1) ), 8 )
  eiky0(1)  = cmplx( cos( c2 * pos_ip0(2) ), sin( -c2 * pos_ip0(2) ), 8 )
  eikz0(1)  = cmplx( cos( c3 * pos_ip0(3) ), sin( -c3 * pos_ip0(3) ), 8 )
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
  eikx1(1)  = cmplx( cos( c1 * pos_ip1(1) ), sin( -c1 * pos_ip1(1) ), 8 )
  eiky1(1)  = cmplx( cos( c2 * pos_ip1(2) ), sin( -c2 * pos_ip1(2) ), 8 )
  eikz1(1)  = cmplx( cos( c3 * pos_ip1(3) ), sin( -c3 * pos_ip1(3) ), 8 )
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

  tol = 5                

  do i = 1, NN
    if (mod(Lx2,2)==0) then
      if (pos(i,1)>=Lx2/2) then
        posq(i,1) = ( pos(i,1) - Lx2 ) / 2.D0
      elseif(pos(i,1)<-Lx2/2) then
        posq(i,1) = ( pos(i,1) + Lx2 ) / 2.D0
      end if
    else
      if (pos(i,2)>=Ly2/2) then
        posq(i,2) = ( pos(i,2) - Ly2 ) / 2.D0
      elseif(pos(i,2)<-Ly2/2) then
        posq(i,2) = ( pos(i,2) + Ly2 ) / 2.D0
      end if
    end if
    if (mod(Ly2,2)==0) then
      if (pos(i,2)>=Ly2/2) then
        posq(i,2) = ( pos(i,2) - Ly2 ) / 2.D0
      elseif(pos(i,2)<-Ly2/2) then
        posq(i,2) = ( pos(i,2) + Ly2 ) / 2.D0
      end if
    else
      if (pos(i,2)>Ly2/2) then
        posq(i,2) = ( pos(i,2) - Ly2 ) / 2.D0
      elseif(pos(i,2)<-Ly2/2) then
        posq(i,2) = ( pos(i,2) + Ly2 ) / 2.D0
      end if
    end if
    posq(i,3) = pos(i,3) / 2.D0
    posq(i,4) = pos(i,4)
  end do

  call Coulomb_energy(EE1)

end subroutine error_analysis_ewald


subroutine read_energy_parameters
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

  epsilon = 1
  sigma = 1
  rc_lj = 1.12
  rv_lj = 2
  rsk_lj = 1
  rsk_real = 1
  R0_2 = 2.25
  kFENE = 30

  allocate(posq(NN,4))

end subroutine read_energy_parameters


subroutine write_energy_parameters
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
  write(*,*) 'real_verlet:', real_verlet
  write(*,*) 'rc_real    :', rc_real
  write(*,*) 'tol        :', tol
  write(*,*) 'tau_rf     :', tau_rf
  write(*,*) 'rc_lj      :', rc_lj
  write(*,*) 'rv_lj      :', rv_lj
  write(*,*) '****************************************************'

end subroutine write_energy_parameters


subroutine Initialize_ewald_parameters
  !--------------------------------------!
  !
  !--------------------------------------!
  use global_variables
  implicit none
  real*8 :: rho, v_verlet

  alpha    = ( tau_rf * pi**3 * Nq / (Lx*Ly)**2/Lz/Z_empty ) ** (1.D0/6)
  alpha2   = alpha * alpha
  rc_real  = tol / alpha
  rc_real2 = rc_real * rc_real
  rv_real  = rc_real + rsk_real
  !
  !use verlet list in real space
  if ( (int(Lx/rv_real) * int(Ly/rv_real) * int(Lz/rv_real)) > 27 ) then 
    Kmax1 = ceiling(tol*Lx*alpha/pi)
    Kmax2 = ceiling(tol*Ly*alpha/pi)
    Kmax3 = ceiling(tol*Lz*Z_empty*alpha/pi)
    real_verlet = 1
  !
  !don't use verlet list in real space
  else
    if ( Lx > Ly ) then
      rc_real = Ly/2
    else
      rc_real = Lx/2
    end if
    rc_real2 = rc_real * rc_real
    Kmax1    = ceiling(tol*tol/pi*Lx/rc_real)
    Kmax2    = ceiling(tol*tol/pi*Ly/rc_real)
    Kmax3    = ceiling(tol*tol/pi*Lz*Z_empty/rc_real)
    alpha    = tol / rc_real
    alpha2   = alpha * alpha
    rv_real  = rc_real + rsk_real
    real_verlet = 0
  end if
  !
  !allocate verlet list of real space
  if ( allocated(real_point) ) deallocate(real_point)
  allocate( real_point(Nq) )
  real_point = 0
  rho = Nq / (Lx * Ly * Lz)
  v_verlet = 8.D0/3 * pi * rv_real**3
  if ( allocated(real_pair_list) ) deallocate(real_pair_list)
  allocate( real_pair_list(25*Nq*ceiling(rho*v_verlet)) )
  real_pair_list = 0

end subroutine Initialize_ewald_parameters


subroutine update_rhok
  !--------------------------------------!
  !
  !--------------------------------------!
  use global_variables
  implicit none

  rho_k = rho_k + Conjg( delta_rhok )
  
end subroutine update_rhok


subroutine build_charge
  !--------------------------------------!
  !Initialize and charge.
  !   
  !Input
  !   pos
  !Output
  !   charge
  !External Variables
  !   pos, charge, NN, Nq
  !--------------------------------------!
  use global_variables
  implicit none
  integer :: i, j

  allocate( charge(Nq) )
  allocate( inv_charge(NN) )
  inv_charge = 0
  
  j=0
  do i = 1, NN
    if ( pos( i, 4 ) /= 0 ) then
      j             = j + 1
      charge(j)     = i
      inv_charge(i) = j
    end if 
  end do

end subroutine build_charge


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
    zq(m) = posq(n,4)
  end do

  c1 = 2*pi/Lx
  c2 = 2*pi/Ly
  c3 = 2*pi/Lz/Z_empty
  do m = 1, Nq
    i = charge(m)
    eikx(m,0)  = (1,0)
    eiky(m,0)  = (1,0)
    eikz(m,0)  = (1,0)

    eikx(m,1)  = cmplx( cos(c1*posq(i,1)), sin(c1*posq(i,1)), 8 )
    eiky(m,1)  = cmplx( cos(c2*posq(i,2)), sin(c2*posq(i,2)), 8 )
    eikz(m,1)  = cmplx( cos(c3*posq(i,3)), sin(c3*posq(i,3)), 8 )

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


subroutine build_fene_list
  !--------------------------------------!
  !Construct the fene_list array.
  !   
  !Input
  !   
  !Output
  !   fene_list, fene_point
  !External Variables
  !   N_bond, Npe, Nml, Ngl
  !--------------------------------------!
  use global_variables
  implicit none
  integer :: i, j, k, l

  N_bond = Ngl * (Nml-1) * 2
  allocate( fene_list(N_bond) )
  allocate( fene_point(Npe) )

  l = 0
  do i = 1, Npe
    if ( mod( i, Nml ) == 1 ) then
      l = l + 1
      fene_list(l)  = i + 1 
      fene_point(i) = l
    elseif( mod( i, Nml ) == 0 ) then
      l = l + 1
      fene_list(l)  = i - 1
      fene_point(i) = l
    else
      l = l + 1
      fene_list(l)  = i - 1
      l = l + 1
      fene_list(l)  = i + 1
      fene_point(i) = l
    end if
  end do

end subroutine build_fene_list


subroutine rij_and_rr(rij, rsqr, i, j)
  !-----------------------------------------!
  !compute displacement vector and displacement of two particles
  !input:
  !  post(pos or pos1), i, j(particle number) 
  !output:
  !  rij(displacement vecter), rr(square of displacement)
  !External Variant:
  !  Lz(used in period condition)
  !note:
  !  including period condition
  !-----------------------------------------!
  use global_variables
  implicit none
  real*8, dimension(3), intent(out) :: rij
  real*8, intent(out) :: rsqr
  integer, intent(in) :: i
  integer, intent(in) :: j

  rij = posq(i,1:3) - posq(j,1:3)

  if ( rij(1) > Lx/2 ) then
    rij(1) = rij(1) - Lx
  elseif( rij(1) <= -Lx/2 ) then
    rij(1) = rij(1) + Lx
  end if
  if ( rij(2) > Ly/2 ) then
    rij(2) = rij(2) - Ly
  elseif( rij(2) <= -Ly/2 ) then
    rij(2) = rij(2) + Ly
  end if

  rsqr = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)

end subroutine rij_and_rr


end module compute_energy_ewald


















