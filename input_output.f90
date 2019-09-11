module input_output
implicit none 

save

    integer, allocatable, dimension(:,:)  :: hist1   !brushes distribution
!   real*8, allocatable, dimension(:,:)  :: hist2   !height distribution
!   real*8, allocatable, dimension(:,:)  :: hist3   !ions distribution
!   real*8, allocatable, dimension(:,:)  :: hist4   !cos distribution
!   real*8, allocatable, dimension(:,:)  :: hist5   !max height distribution
!   real*8, allocatable, dimension(:,:)  :: hist6   !polar angle distribution
!   integer, allocatable, dimension(:,:) :: hist7   !top view of monomer
!   integer, allocatable, dimension(:,:) :: hist8   !side view of monomer
!   integer, allocatable, dimension(:,:) :: hist9   !side view of counterion

contains


subroutine initialize_parameters
  !--------------------------------------!
  !Initialize system parameters
  !
  !Input
  !  none
  !Output
  !  none
  !External Variables
  !  none
  !Routine Referenced:
  !  none
  !--------------------------------------!
  implicit none

  call read_data

  call data_operation

  call write_data

  call data_allocate

  call initialize_bond_vector

  call initialize_move

end subroutine initialize_parameters


subroutine read_data
  !--------------------------------------!
  !read system parameters
  !
  !Input
  !  none
  !Output
  !  none
  !External Variables
  !  parameters in global variables
  !Routine Referenced:
  !  none
  !--------------------------------------!
  use global_variables
  implicit none
  !
  !read system data
  open(10,file='./system_data.txt')
    read(10,*) Lz
    read(10,*) Z_empty  
    read(10,*) sigmag   
    read(10,*) Beta       
    read(10,*) qq
    read(10,*) arm
    read(10,*) Nma
    read(10,*) Nga
    read(10,*) man_s
    read(10,*) StepNum0           
    read(10,*) StepNum            
    read(10,*) DeltaStep1
    read(10,*) DeltaStep2     
    read(10,*) DeltaStep3 
    read(10,*) pH_pKa        
  close(10)
end subroutine read_data


subroutine data_operation
  !--------------------------------------!
  !Initialize system parameters
  !and judge whether restarted or continue
  !
  !Input
  !  
  !Output
  !   
  !External Variables
  !  
  !Routine Referenced:
  !
  !--------------------------------------!
  use global_variables
  implicit none
  logical alive
  integer :: i, nx, dLx

  !
  !The total monomers of star brushes
  Npe = ( Nma*arm + 1 ) * Nga
  !
  !The total charged monomers on PE
  if ( abs(qq) == 0 ) then
    Nq_PE = 0
  else
    if ( man_s /= 0) then
      Nq_PE = Nma / man_s * arm * Nga
    else
      Nq_PE = 0
    end if
  end if

  !The total charged particles Nq and total particles NN in system
  Nq = Nq_PE * 2
  NN = Npe + Nq_PE
  !
  !System size, keep mod(Lx,nx)=0
  Lx = nint(sqrt( Nga / sigmag )) * 2
  nx = nint( sqrt(1.*Nga) )
  Lx = Lx - mod( Lx, nx )
  Ly = Lx
  Z_empty = ( 1.*Lz + 1.*Lz * Z_empty ) / (1.*Lz)
  !
  !number of bonds in system
  N_bond = Nma * arm * Nga
  !
  !whether continue or restart
  Inquire( file='start_time.txt', exist=alive )
  if (alive) then
    open(11,file='./start_time.txt')
      read(11,*) restart_or_continue
    close(11)
  else
    restart_or_continue = 0
  end if
end subroutine data_operation


subroutine write_data
  !--------------------------------------!
  !Write parameters to screen
  !
  !Input
  !  
  !Output
  !   
  !External Variables
  !  
  !Routine Referenced:
  !
  !--------------------------------------!
  use global_variables
  implicit none

  write(*,*)
  write(*,*)
  write(*,*) '****************************************************************'
  write(*,*) '*************************system_data****************************'
  write(*,*) 'Incert space in z direction,       Z_empty:', Z_empty-1
  write(*,*) 'Grafting density,                   sigmag:', sigmag
  write(*,*) 'Beta=1/kT,                            Beta:', Beta
  write(*,*) 'Charges of polymer,                     qq:', qq
  write(*,*) 'Arms of star brushes                   Arm:', arm
  write(*,*) 'Monomers of each arm                   Nma:', Nma
  write(*,*) 'Number of grafted star chains          Nga:', Nga
  write(*,*) 'Each man_s monomers with one charge, man_s:', man_s
  write(*,*) 'total particles,                        NN:', NN
  write(*,*) 'total charged particles,                Nq:', Nq
  write(*,*) 'total brushes particles,               Npe:', Npe
  write(*,*) 'Lattice number in x direction,          Lx:', Lx
  write(*,*) 'Lattice number in y direction,          Ly:', Ly
  write(*,*) 'Lattice number in z direction,          Lz:', Lz
  write(*,*) 'Number of bonds in system,          N_bond:', N_bond
  write(*,*) 'pH-pKa,                             pH-pKa:', pH_pKa
  write(*,*) '****************************************************************'

  write(*,*)
  write(*,*) '****************************************************************'
  write(*,*) '************************running_steps***************************'
  write(*,*) 'restart (0), continue (0), restart_continue:',restart_or_continue
  write(*,*) 'Preheating steps                   StepNum0:', StepNum0
  write(*,*) 'Running steps                       StepNum:', StepNum
  write(*,*) 'Total steps                StepNum0+StepNum:', (StepNum0+StepNum)
  write(*,*) 'Step inteval, physical quantum,  DeltaStep1:', DeltaStep1
  write(*,*) 'Step inteval, histogram,         DeltaStep2:', DeltaStep2
  write(*,*) 'Step inteval, output data,       DeltaStep3:', DeltaStep3
  write(*,*) '****************************************************************'
  write(*,*)
  write(*,*)

end subroutine write_data

subroutine data_allocate
  !--------------------------------------!
  !Allocate pos, vel, acc and histogram arrays
  !Initialize histogram arrays
  !
  !Input
  !  
  !Output
  !   
  !External Variables
  !  Arrays in Global Variables
  !Routine Referenced:
  !  David P. Landau, Kurt Binder. A Guide to Monte Carlo Simulation in
  !  Statistical Physics, 3rd (Cambridge, UK: Cambridge University Press, 
  !  2009), 453.
  !--------------------------------------!
  use global_variables
  implicit none
  integer :: i
  !
  !position, velocity, and acceleration
  allocate( pos(NN,4) )
  allocate( latt(Lx,Ly,Lz+1) )
  pos = 0
  latt = 0

  allocate( ipx(Lx) )
  allocate( ipy(Ly) )
  allocate( ipz(Lz) )

  allocate( ip2x(Lx) )
  allocate( ip2y(Ly) )
  allocate( ip2z(Lz) )

  allocate( imx(Lx) )
  allocate( imy(Ly) ) 
  allocate( imz(Lz) )

  allocate( bonds(110,3) )
  allocate( move(109,6) )
  allocate( b1(108) )
  allocate( b12(108) )
  allocate( monbd(Npe, arm+1) )
  allocate( bond_numb(N_Bond) )
  monbd = 0

  do i = 1, Lx
    ipx(i) = i+1
    ip2x(i) = i+2
    imx(i) = i-1
  end do
  ipx(Lx) = 1       !periodical boundary condition
  ip2x(Lx-1) = 1
  ip2x(Lx) = 2
  imx(1) = Lx

  do i = 1, Ly
    ipy(i) = i+1
    ip2y(i) = i+2
    imy(i) = i-1
  end do
  ipy(Ly) = 1       !periodical boundary condition
  ip2y(Ly-1) = 1
  ip2y(Ly) = 2
  imy(1) = Ly

  do i = 1, Lz      !finite boundary condition
    ipz(i) = i+1
    ip2z(i) = i+2
    imz(i) = i-1
  end do 
  ip2z(Lz) = Lz + 1
  imz(1) = 1   

  !
  !Allocate arrays and initialize them
  allocate( hist1(Lz, 2)        )
!   allocate( hist2(SizeHist, 2)        )
!   allocate( hist3(SizeHist, 2)        )
!   allocate( hist4(Nml, 2)             )
!   allocate( hist5(SizeHist, 2)        )
!   allocate( hist6(SizeHist, 2)        )
!   allocate( hist7(SizeHist, SizeHist) )
!   allocate( hist8(SizeHist, SizeHist) )
!   allocate( hist9(SizeHist, SizeHist) )
  hist1=0
!   hist2=0
!   hist3=0
!   hist4=0
!   hist5=0
!   hist6=0
!   hist7=0
!   hist8=0
!   hist9=0        

end subroutine data_allocate


subroutine continue_read_data(l)
  !------------------------------------!
  !
  !------------------------------------!
  use global_variables
  implicit none
  integer, intent(out) :: l
  integer :: i, j 

!   open(20,file='./data/pos1.txt')
!     read(20,*) ((pos(i,j),j=1,4),i=1,NN)
!   close(20)
!   open(19,file='./start_time.txt')
!     read(19,*)
!     read(19,*) l
!     read(19,*) total_time
!   close(19)
  open(22,file='./data/hist1.txt')
!   open(23,file='./data/hist2.txt')
!   open(24,file='./data/hist3.txt')
!   open(25,file='./data/hist4.txt')
!   open(26,file='./data/hist5.txt')
!   open(27,file='./data/hist6.txt')
!   open(28,file='./data/hist7.txt')
!   open(29,file='./data/hist8.txt')
!   open(30,file='./data/hist9.txt')
    read(22,*) ((hist1(i,j),j=1,2),i=1,Lz)
!     read(23,*) ((hist2(i,j),j=1,2),i=1,SizeHist)
!     read(24,*) ((hist3(i,j),j=1,2),i=1,SizeHist)
!     read(25,*) ((hist4(i,j),j=1,2),i=1,Nml)
!     read(26,*) ((hist5(i,j),j=1,2),i=1,SizeHist)
!     read(27,*) ((hist6(i,j),j=1,2),i=1,SizeHist)
!     read(28,*) ((hist7(i,j),j=1,SizeHist),i=1,SizeHist)
!     read(29,*) ((hist8(i,j),j=1,SizeHist),i=1,SizeHist)
!     read(30,*) ((hist9(i,j),j=1,SizeHist),i=1,SizeHist)
!   close(30)
!   close(29)
!   close(28)
!   close(27)
!   close(26)
!   close(25)
!   close(24)
!   close(23)
  close(22)
end subroutine continue_read_data


subroutine initialize_bond_vector
  !--------------------------------------!
  ! The bonds.txt file are generated by MATLAB according to the subroutine
  ! bdibfl in pp. 450-452 in Landau's book.
  !
  !Input
  !  
  !Output
  !   
  !External Variables
  !  Arrays in Global Variables
  !Routine Referenced:
  !  David P. Landau, Kurt Binder. A Guide to Monte Carlo Simulation in
  !  Statistical Physics, 3rd (Cambridge, UK: Cambridge University Press, 
  !  2009), see Subroutine bdibfl in pp. 450-452.
  !--------------------------------------!
  use global_variables
  implicit none
  integer :: i, j

  open(100,file='bonds.txt')
    read(100,*) ((bonds(i,j),j=1,3),i=1,108)
  close(100)
  !
  !testify bonds array
  !   do i = 1, 108
  !    write(*,*) bonds(i,:)
  !   end do

  do i = 1, 108
    b12(i) = dot_product(bonds(i,:), bonds(i,:))
    b1(i)  = sqrt(b12(i))
  end do

  !testify b1, b12 array
  !   do i = 1, 108
  !     write(*,*) b1(i), b12(i)
  !   end do

  new_direction(1,1) = 1
  new_direction(1,2) = 0
  new_direction(1,3) = 0

  new_direction(2,1) = 0
  new_direction(2,2) = 1
  new_direction(2,3) = 0

  new_direction(3,1) = 0
  new_direction(3,2) = 0
  new_direction(3,3) = 1  

  new_direction(4,1) = 0
  new_direction(4,2) = 0
  new_direction(4,3) = -1

  new_direction(5,1) = 0
  new_direction(5,2) = -1 
  new_direction(5,3) = 0   

  new_direction(6,1) = -1
  new_direction(6,2) = 0
  new_direction(6,3) = 0

end subroutine initialize_bond_vector


subroutine initialize_move
  !--------------------------------------!
  ! To generate the array move(109,6) whose element move(i,j) is a bond number 
  ! belonging to bonds array after the original bond, i-th bond in bonds array,
  ! with the j-th move that have 6 probable choice.
  !
  !Input
  !  
  !Output
  !   
  !External Variables
  !  
  !Routine Referenced:
  ! David P. Landau, Kurt Binder. A Guide to Monte Carlo Simulation in
  ! Statistical Physics, 3rd (Cambridge, UK: Cambridge University Press, 
  ! 2009), see Subroutine inimove in pp. 461.
  !Note:
  ! There are little errors in Landau's book, and are corrected here.
  !--------------------------------------!
  use global_variables
  implicit none
  integer :: i, j, k
  integer, dimension(6,3) :: new
  logical :: test

  move = 0     ! move(i,j)==0, means the new bond is broken
  do i = 1, 108
    new(1,1) = bonds(i,1) + 1
    new(1,2) = bonds(i,2)
    new(1,3) = bonds(i,3)

    new(2,1) = bonds(i,1)
    new(2,2) = bonds(i,2) + 1
    new(2,3) = bonds(i,3)

    new(3,1) = bonds(i,1)
    new(3,2) = bonds(i,2)
    new(3,3) = bonds(i,3) + 1  

    new(4,1) = bonds(i,1)
    new(4,2) = bonds(i,2)
    new(4,3) = bonds(i,3) - 1

    new(5,1) = bonds(i,1)
    new(5,2) = bonds(i,2) - 1 
    new(5,3) = bonds(i,3)    

    new(6,1) = bonds(i,1) - 1
    new(6,2) = bonds(i,2)
    new(6,3) = bonds(i,3)

    do j = 1, 6
      test = .false.
      do k = 1, 108 
        test = ( new(j,1) .eq. bonds(k,1) ) .and. &
  &            ( new(j,2) .eq. bonds(k,2) ) .and. &
  &            ( new(j,3) .eq. bonds(k,3) )
        if (test) then
          move(i,j) = k ! move(i,j)==k, means the new bond 
                        ! is k bond in bonds array
          exit
        end if
      end do
    end do
  end do

  do i = 1, 6
    move(109,i) = 109   ! move(i,j)==109, means the end bond of one chain
  end do

  !testify the move array
  !   do i = 1, 109
  !     write(*,*) move(i,:)
  !   end do

end subroutine initialize_move


subroutine histogram
  !----------------------------------------!
  !input:
  !  pos
  !output:
  !  hist1(distribution hisotgram from PE to rod)
  !External Variants:
  !  Npe, Ngl, NN, Nml, Sizehist, Lz 
  !----------------------------------------!
  use global_variables
  implicit none
  integer :: i, j, k
  real*8, dimension(3) :: rij
  real*8 rsqr, max_h, theta, rr
  
  !
  !hist1
  do i=1, Npe
    if ( mod(i,(arm*Nma+1))==1 .and. i<=Npe ) cycle
    k = pos(i,3)
    if ( k<=0 .or. k>Lz ) then
      write(*,*) 'Wrong in histogram1'
      cycle
    end if
    hist1(k,2) = hist1(k,2) + 1
  end do
!   !
!   !hist2
!   do i = 1, Nga
!     k = ceiling( pos((i-1)*(Nma*arm+1)+Nma+1,3) / (1.*Lz/SizeHist) )
!     if ( k<=0 .or. k>SizeHist ) then
!       write(*,*) 'Wrong in histogram2'
!       cycle
!     end if
!     hist2(k,2) = hist2(k,2) + 1
!   end do
!   !
!   !hist3
!   do i = 1, Ngl*nint(abs(qq))
!     k = ceiling( pos(Ngl*Nml+i,3) / (Lz/SizeHist) )
!     if ( k<=0 .or. k>SizeHist ) then
!       write(*,*) 'Wrong in histogram3'
!       cycle
!     end if
!     hist3(k,2) = hist3(k,2) + 1
!   end do
!   !
!   !hist4
!   do i=1, Ngl
!     do j=2, Nml
!       call rij_and_rr(rij, rsqr, (i-1)*Nml+j, (i-1)*Nml+j-1)
!       hist4(j,2) = hist4(j,2) + rij(3)/sqrt(rsqr)/Ngl
!     end do
!   end do
!   !
!   !hist5
!   do i = 1, Ngl
!     max_h = 0
!     do j = 1, Nml
!       if ( max_h < pos((i-1)*Nml+j,3) ) then
!         max_h = pos((i-1)*Nml+j,3)
!       end if
!     end do
!     k = ceiling( max_h / (Lz/SizeHist) )
!     if ( k<=0 .or. k>SizeHist ) then
!       write(*,*) 'Wrong in histogram5'
!       cycle
!     end if
!     hist5(k,2) = hist5(k,2) + 1
!   end do
!   !
!   !hist6
!   do i = 1, Ngl
!     call rij_and_rr(rij,rr,i*Nml,(i-1)*Nml+1)
!     theta = acos( pos(i*Nml,3) / sqrt(rr) )
!     k = ceiling( theta / (pi/2/SizeHist) )
!     if ( k<=0 .or. k>SizeHist ) then
!       write(*,*) 'Wrong in histogram6'
!       cycle
!     end if
!     hist6(k,2) = hist6(k,2) + 1
!   end do
!   !
!   !hist7
!   do k = 1, Npe
!     if ( mod(k,Nml)==1 .and. k<=Npe ) cycle
!     i = ceiling( (pos(k,1)+Lx/2) / (Lx/SizeHist) )
!     j = ceiling( (pos(k,2)+Ly/2) / (Ly/SizeHist) )
!     if (i<=0 .or. i>SizeHist .or. j<=0 .or. j>SizeHist) then
!       write(*,*) 'Wrong in histogram7'
!       cycle
!     end if
!     hist7(i,j)=hist7(i,j)+1
!   end do
!   !
!   !hist8
!   do k = 1, Npe
!     if ( mod(k,Nml)==1 .and. k<=Npe ) cycle
!     i = ceiling( (pos(k,2)+Ly/2) / (Ly/SizeHist) )
!     j = ceiling( pos(k,3) / (Lz/SizeHist) )
!     if ( i<=0 .or. i>SizeHist .or. j<=0 .or. j>SizeHist ) then
!       write(*,*) 'Wrong in histogram8'
!       write(*,*) i, j, k, pos(k,1:3)
!       cycle
!     end if
!     hist8(i,j) = hist8(i,j) + 1
!   end do
!   !
!   !hist9
!   do k = Npe+1, NN
!     i = ceiling( (pos(k,2)+Ly/2) / (Ly/SizeHist) )
!     j = ceiling( pos(k,3) / (Lz/SizeHist) )
!     if ( i<=0 .or. i>SizeHist .or. j<=0 .or. j>SizeHist ) then
!       write(*,*) 'Wrong in histogram9'
!       cycle
!     end if
!     hist9(i,j) = hist9(i,j) + 1
!   end do

end subroutine histogram


subroutine write_pos
  use global_variables
  implicit none
  integer :: i, j, k

  open(100,file='./data/pos.txt')
    do i = 1, NN
      write(100,'(4I10)') (pos(i,j),j=1,4)
    end do
  close(100)

  open(100,file='./data/bond_numb.txt')
    do i = 1, N_Bond
      write(100,*) bond_numb(i)
    end do
  close(100)

  open(100,file='./data/monbd.txt')
    do i = 1, Npe
      write(100,*) monbd(i,:)
    end do
  close(100)

  open(100,file='./data/move.txt')
    do i = 1, 109
      write(100,*) move(i,:)
    end do
  close(100)  

  open(100, file='./data/ipx.txt')
    do i = 1, Lx
      write(100,*) imx(i),ipx(i),ip2x(i)
    end do
  close(100)

  open(100, file='./data/ipy.txt')
    do i = 1, Ly
      write(100,*) imy(i),ipy(i),ip2y(i)
    end do
  close(100)

  open(100, file='./data/ipz.txt')
    do i = 1, Lz
      write(100,*) imz(i),ipz(i),ip2z(i)
    end do
  close(100)

  open(100,file='./data/latt.txt')
    do i = 1, Lx
      do j = 1, Ly
        do k = 1, Lz+1
          write(100,'(4I6)') i,j,k,latt(i,j,k)
        end do
      end do
    end do
  close(100)

end subroutine write_pos

subroutine write_pos1
  use global_variables
  implicit none
  integer :: i, j, k

  open(100,file='./data/pos1.txt')
    do i = 1, NN
      write(100,'(4I10)') (pos(i,j),j=1,4)
    end do
  close(100)

  open(100,file='./data/bond_numb1.txt')
    do i = 1, N_Bond
      write(100,*) bond_numb(i)
    end do
  close(100)

  open(100,file='./data/latt1.txt')
    do i = 1, Lx
      do j = 1, Ly
        do k = 1, Lz+1
          write(100,'(4I6)') i,j,k,latt(i,j,k)
        end do
      end do
    end do
  close(100)

end subroutine write_pos1


subroutine write_hist
  !----------------------------------------!
  !Write distribution histogram to the file hist1.txt ... hist4.txt
  !input:
  !  hist1, hist2, hist3, hist4,...hist9
  !External Variants: 
  !  SizeHist, Lz, 
  !----------------------------------------!
  use global_variables
  implicit none
  integer i,j
  
  open(34,file='./data/hist1.txt')
!   open(35,file='./data/hist2.txt')
!   open(36,file='./data/hist3.txt')
!   open(37,file='./data/hist5.txt')
!   open(38,file='./data/hist6.txt')
    do i=1,Lz
      write(34,340) i, hist1(i,2)
!       hist2(i,1)=i*Lz/SizeHist
!       write(35,340) hist2(i,1), hist2(i,2)
!       hist3(i,1)=i*Lz/SizeHist
!       write(36,340) hist3(i,1), hist3(i,2)
!       hist5(i,1)=i*Lz/SizeHist
!       write(37,340) hist5(i,1), hist5(i,2)
!       hist6(i,1)=i*pi/2/SizeHist
!       write(38,340) hist6(i,1), hist6(i,2)
    end do
    340 format(2I10)
  close(34)
!   close(35)
!   close(36)
!   close(37)
!   close(38)

!   open(39,file='./data/hist4.txt')
!     do i=1,Nml
!       hist4(i,1)=i*1.
!       write(39,340) hist4(i,1), hist4(i,2)
!     end do
!   close(39)

!   open(40,file='./data/hist7.txt')
!   open(41,file='./data/hist8.txt')
!   open(42,file='./data/hist9.txt')
!     do i=1,SizeHist
!       write(40,'(500I10)') (hist7(i,j),j=1,SizeHist) 
!       write(41,'(500I10)') (hist8(i,j),j=1,SizeHist)  
!       write(42,'(500I10)') (hist9(i,j),j=1,SizeHist)  
!     end do
!   close(40)
!   close(41)
!   close(42)

end subroutine write_hist


subroutine write_time(time)
  !------------------------------------!
  !
  !------------------------------------!
  use global_variables
  implicit none

  real*8, intent(in) :: time
  open(10,file='./data/time.txt')
    write(10,*) 'time:(seconds)', real(total_time)
    write(10,*) 'time:(hours)  ', real(total_time/3600)
    write(10,*) 'time:(days)   ', real(total_time/86400)
    write(10,*) 'Lx:           ', real(Lx)
    write(10,*) 'Ly:           ', real(Ly)
    write(10,*) 'Lz:           ', real(Lz)
    write(10,*) 'Nga:          ', Nga
    write(10,*) 'Nma:          ', Nma
    write(10,*) 'Nq:           ', Nq
    write(10,*) 'NN:           ', NN
  close(10)

end subroutine write_time


end module input_output


