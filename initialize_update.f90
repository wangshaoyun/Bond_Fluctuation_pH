module initialize_update
implicit none

contains

subroutine initialize_position
  use global_variables
  implicit none

  !
  !Graft star brushes, the arms are all vertical to the plate.
  !The start point on plate can be uniform of random
  call uniform_star_brushes
  !
  !random distribution in the box
  if ( qq /= 0 ) then
    call initialize_ions
  end if

!   !
!   !initialize chargen on PE
!   do i=1, Nq_PE
!     pos(charge(i),4) = qq                 
!   end do

end subroutine initialize_position


subroutine uniform_star_brushes
  use global_variables
  implicit none
  integer :: i, j, k, l, m, n
  integer :: nx, dLx
  integer :: base, base1, base2
  integer :: xi, yi, zi, xp, yp, zp
  integer, dimension(5,4) :: bond_vector

  nx = nint( sqrt(1.*Nga) )
  dLx = Lx / nx
  !up, bond_vector is bond number in bonds array by looking up bonds.txt
  bond_vector(1,1) = 0
  bond_vector(1,2) = 0
  bond_vector(1,3) = 2
  bond_vector(1,4) = 3     !bond number in bonds.txt
  !right
  bond_vector(2,1) = 2
  bond_vector(2,2) = 0
  bond_vector(2,3) = 0 
  bond_vector(2,4) = 1     
  !forward
  bond_vector(3,1) = 0
  bond_vector(3,2) = 2
  bond_vector(3,3) = 0 
  bond_vector(3,4) = 5     
  !left
  bond_vector(4,1) = -2
  bond_vector(4,2) = 0
  bond_vector(4,3) = 0   
  bond_vector(4,4) = 2     
  !backward
  bond_vector(5,1) = 0
  bond_vector(5,2) = -2
  bond_vector(5,3) = 0
  bond_vector(5,4) = 6     

  do i = 1, nx
    do j = 1, nx
      l = (i-1)*nx + j
      base = (l-1)*(arm*Nma+1)
      pos(base+1, 1) = (i-1)*dLx + 3
      pos(base+1, 2) = (j-1)*dLx + 3
      pos(base+1, 3) = 1
      monbd(base+1,1) = base + 2 - l    
      monbd(base+1,arm+1) = 1
      do k = 2, Nma + 1
        pos(base+k, :) = pos(base+k-1,:) + bond_vector(1,1:3)
        bond_numb(base+k-l) = bond_vector(1,4)
        monbd(base+k,1) = base + k - l       
        monbd(base+k,2) = base + k + 1 - l 
        monbd(base+k,arm+1) = 2
      end do
      do m = 2, arm
        base1 = base + Nma + 1      
        base2 = base1 + (m-2) * Nma 
        pos(base2+1, :) = pos(base1, :) + bond_vector(m-1,1:3)
        if (pos(base2+1, 1)>Lx) then
          pos(base2+1,1) = pos(base2+1,1)-Lx
        end if
        if (pos(base2+1, 2)>Ly) then
          pos(base2+1,2) = pos(base2+1,2)-Ly
        end if
        bond_numb(base2+1-l) = bond_vector(m,4)
        monbd(base1,m) = base2 + 1 - l     
        monbd(base1,arm+1) = arm
        monbd(base2+1,1) = base2 + 1 - l   
        monbd(base2+1,2) = base2 + 2 - l    
        monbd(base2+1,arm+1) = 2
        do n = 2, Nma
          pos(base2+n,:) = pos(base2+n-1,:) + bond_vector(1,1:3)
          bond_numb(base2+n-l) = bond_vector(1,4)
          monbd(base2+n,1) = base2 + n - l                
          monbd(base2+n,2) = base2 + n + 1 - l            
          monbd(base2+n,arm+1) = 2
        end do
        monbd(base2+Nma,2) = 0               
        monbd(base2+Nma,arm+1) = 1
      end do   
    end do
  end do

  !the vortexes of the lattice with particle are occupied
  do i = 1, Npe
    xi = pos(i,1)
    yi = pos(i,2)
    zi = pos(i,3)
    if (zi == 0) then
      cycle
    end if
    xp = ipx(xi)
    yp = ipy(yi)
    zp = ipz(zi)
    latt(xi,yi,zi) = 1
    latt(xi,yi,zp) = 1
    latt(xi,yp,zi) = 1
    latt(xi,yp,zp) = 1
    latt(xp,yi,zi) = 1
    latt(xp,yi,zp) = 1
    latt(xp,yp,zi) = 1
    latt(xp,yp,zp) = 1
  end do

  !the above and bottom plates are occupied
  do i = 1, Lx
    do j = 1, Ly
      latt(i,j,1) = 1
      latt(i,j,Lz+1) = 1
    end do
  end do

end subroutine uniform_star_brushes


subroutine initialize_ions
  use global_variables
  implicit none
  integer :: i, j, k
  integer :: total
  integer :: xi, yi, zi, xp, yp, zp
  real*8, dimension(3) :: rnd
  logical :: test

  do i = Npe+1, NN
    test=.true.
    do while (test)
      call random_number(rnd)
      xi = floor(rnd(1)*Lx) + 1
      yi = floor(rnd(2)*Ly) + 1
      zi = floor(rnd(3)*Lz) + 1
      xp = ipx(xi)
      yp = ipy(yi)
      zp = ipz(zi)
      if (zp==0) then
        write(*,*) zp
      end if
      total = latt(xi,yi,zi)+latt(xi,yi,zp)+latt(xi,yp,zi)+latt(xi,yp,zp) +   &
              latt(xp,yi,zi)+latt(xp,yi,zp)+latt(xp,yp,zi)+latt(xp,yp,zp)
      if ( total == 0 ) then
        pos(i,1) = xi
        pos(i,2) = yi
        pos(i,3) = zi
        pos(i,4) = -qq
        latt(xi,yi,zi) = 1
        latt(xi,yi,zp) = 1
        latt(xi,yp,zi) = 1
        latt(xi,yp,zp) = 1
        latt(xp,yi,zi) = 1
        latt(xp,yi,zp) = 1
        latt(xp,yp,zi) = 1
        latt(xp,yp,zp) = 1        
        test = .false.
      end if
    end do
  end do

end subroutine initialize_ions


subroutine monte_carlo_move
  use global_variables
  implicit none
  integer :: i

  do i = 1 : NN - Nga
    call choose_particle
    call new_position
    call delta_energy
    call move_or_not
  end do

end subroutine monte_carlo_move


subroutine choose_particle
  use global_variables
  implicit none
  real*8 :: rnd
  integer :: Ns

  call random_number(rnd)
  ip = floor(rnd*NN) + 1

  Ns = arm*Nma+1
  !
  !The monomer anchored on the plate can't move, so we need to choose again.
  do while( mod(ip,Ns) == 1 .and. ip <= Npe )
    call random_number(rnd)
    ip = int(rnd*NN) + 1
  end do

end subroutine choose_particle


subroutine new_position
  use global_variables
  implicit none
  integer :: dir  !direction, with 6 choice
  integer :: bn   ! number of bonds connect to the particle
  integer, allocatable, dimension(:) :: new_bonds
  integer :: i
  real*8  :: rnd

  call random_number(rnd)
  dir = floor(6*rnd)+1

  bn = monbd(ip,arm+1)
  allocate(new_bonds(bn))

  new_bonds(1) = move( bond_number( monbd( ip, 1 ) ), dir )
  do i = 2, bn
    new_bonds(i) = move( bond_number( monbd( ip, i ) ), 7-dir )
  end do

  pos_ip0 = pos(ip,:)
  pos_ip1(1:3) = pos_ip0(1:3) + new_direction(dir,:)
  pos_ip1(4)   = pos_ip0(4)
  call periodic_condition( pos_ip1(1:2) )

end subroutine new_position


end module initialize_update

