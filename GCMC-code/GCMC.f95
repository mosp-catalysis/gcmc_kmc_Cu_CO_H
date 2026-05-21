! ! GCMC code, for Cu100-eCO2RR work
! ! by Shuoqi Zhang, written in 2024.9.19-2024.9.29
! ! hope it can work :> 

subroutine init_random_seed(seed_flag)   ! random seed
  implicit none
  
  integer(kind=1), intent(IN) :: seed_flag
  integer(kind=4) :: i, n, clock
  integer, dimension(:), allocatable :: seed  ! ArrayList

  ! seed_flag = 0：generate new seed and record
  if(seed_flag == 0) then
    call random_seed(size=n)
    allocate(seed(n))
    call system_clock(count=clock)
    seed = clock + 37*(/(i - 1, i = 1, n)/)
    call random_seed(put=seed)
    open (1, file='seed.txt', status='replace', action='write')
      write (1, *) n
      write (1, *) clock
    close (1)
  ! seed_flag = 1：use old seed
  elseif(seed_flag == 1) then
    open (1, file='seed.txt', status='old', action='read')
      read (1, *) n
      read (1, *) clock
    close (1)
    allocate(seed(n))
    seed = clock + 37*(/(i - 1, i=1, n)/)
    call random_seed(put=seed)
  end if
  deallocate(seed)
end subroutine

subroutine random_integer(r_int, low, high)
  implicit none

  integer(kind=4), intent(out) :: r_int
  integer(kind=4), intent(in) :: high, low
  real(kind=4) :: r_float

  call random_number(r_float)
  r_int = floor(r_float * (high - low + 1)) + low
end subroutine

  

program main
  implicit none

  !!! 结构部分的参数
  ! flag变量
  integer(kind=1) :: seed_flag
  integer(kind=1), allocatable :: cov_type(:)     ! cov_type的flag意义：0: 空位点。1: 有CO吸附。2: 有H吸附。5：体相
  ! 原子编号计数等变量
  integer(kind=4) :: i, j, temp_int4, num_atoms, num_surfatom, num_H, num_CO
  real(kind=4) :: cov_H, cov_CO
  ! 相应的数组：cn、cn_2、近邻、次近邻列表
  integer(kind=4), allocatable :: cn(:), nnsite(:, :), surf_index(:)
  ! 各种长度；坐标；gcn
  real(kind=4) :: ra, rx, ry, rz, z_surf
  real(kind=4) :: rx_slab, ry_slab, rz_slab
  real(kind=4), allocatable :: x(:), y(:), z(:), gcn(:) 
  ! ASCII型2字符长度字符串：用于存元素
  character(kind=1, len=2), allocatable :: element(:)


  !!! 能量部分的参数
  ! M-adsorbate作用
  ! real(kind=8) :: M_CO_slope, M_CO_intercept, M_H_slope, M_H_intercept
  ! AD-AD作用
  real(kind=8) :: E_COrep, E_Hrep, E_CO_Hrep
  ! 环境
  real(kind=8) :: pCO, Temp, pH, U_SHE
  real(kind=8) :: GT_COgas, GT_COads, GT_H2gas, GT_Hads
  ! 能量相关变量
  real(kind=8) :: E_now, E_test, r_accept, u_perH, u_perCO
  ! 参数
  real(kind=8), parameter :: kb=8.6173324D-05     ! 单位eV/K
  real(kind=8), parameter :: p0 = 100000          ! 单位Pa


  !!! 步数等变量
  ! 步数记录&文件名
  integer(kind=4) :: nLoop, nLoop_peroutput, nstep_tot, n_file
  character(kind=1, len=100) filename_struc, filename_proc, char_seq
  ! 确保随机事件位点的不重复 & 每位点的变化次数计数
  ! integer(kind=1), allocatable :: chosen_site(:)
  integer(kind=4), allocatable :: Nevent_persite(:)
  ! 事件接受相关变量
  integer(kind=1) :: accept_flag, last_steps(100000)  ! flag数组，保存当前步的拒绝状态
  integer(kind=4) :: total_rej
  real(kind=4) :: ratio_rej
  ! 随机数存储
  integer(kind=4) :: site_seq, type_event
  real(kind=8) :: R_random


  ! 声明一个派生类型，来记录发生的事件
  type :: event_record
    integer(kind=4) :: position
    integer(kind=1) :: cov_before
  end type event_record
  type(event_record) :: steps






  !!! initialize from input files
  open(3, file="input", status="old", action='read')
    read(3, *) seed_flag
    read(3, *) rx_slab, ry_slab, rz_slab
    read(3, *) z_surf
    ! read(3, *) M_CO_slope, M_CO_intercept, M_H_slope, M_H_intercept
    read(3, *) E_COrep, E_CO_Hrep, E_Hrep
    read(3, *) GT_COgas, GT_COads, GT_H2gas, GT_Hads
    read(3, *) pCO, Temp, pH, U_SHE
    read(3, *) nLoop, nLoop_peroutput
  close(3)



  ! 化学势值    ! 单位eV
  u_perCO =  kb * Temp * log(pCO / p0) + GT_COgas - GT_COads
  u_perH = -kb * Temp * pH * log(10.0) - U_SHE + 0.5 * GT_H2gas - GT_Hads
  ! debug段
  write(*, *) 'u_perCO:', u_perCO, 'u_perH:', u_perH
  write(*, *)


  ! 算整个slab的rx, ry, rz（用于周期性边界条件）
  ! rx_slab = rx_peratom * tot_Natom_x
  ! ry_slab = ry_peratom * tot_Natom_y
  ! rz_slab = rz_peratom * tot_Natom_z

  open(9, file='ini.xyz', status="old", action='read')
    read(9, *) num_atoms

    ! 初始化待分配数组
    allocate(element(num_atoms))
    element = ''
    allocate(x(num_atoms))
    x = 0.0
    allocate(y(num_atoms))
    y = 0.0
    allocate(z(num_atoms))
    z = 0.0
    allocate(cov_type(num_atoms))
    cov_type = 0
    allocate(cn(num_atoms))
    cn = 0
    allocate(nnsite(12, num_atoms))
    nnsite = 0
    allocate(gcn(num_atoms))
    gcn = 0.0

    read(9, *)
    do i = 1, num_atoms
      read(9, *) element(i), x(i), y(i), z(i), cov_type(i)
    end do
  close(9)

  ! z_surf = (maxval(z)-minval(z)) / 2.0 + minval(z)


  num_surfatom = 0
  write(*,*) "Get CN, NNsite, surface index"
  ! 生成CN、nnsite列表
  do i = 1, num_atoms
    do j = i + 1, num_atoms
      rx = abs(x(i)-x(j))
      ry = abs(y(i)-y(j))
      rz = abs(z(i)-z(j))
      
      ! 处理周期性边界条件
      if(rx * 2 > rx_slab) rx = (rx_slab - rx)
      if(ry * 2 > ry_slab) ry = (ry_slab - ry)
      if(rz * 2 > rz_slab) rz = (rz_slab - rz)
      
      ra = sqrt(rx ** 2 + ry ** 2 + rz ** 2)
      if(ra < 3.0) then
        cn(i) = cn(i) + 1
        cn(j) = cn(j) + 1
        nnsite(cn(i), i) = j
        nnsite(cn(j), j) = i
      end if
    end do


    ! 体相、表面分开，并统计表面原子数
    if(cn(i) == 12 .or. z(i) < z_surf) then
      cov_type(i) = 5
    else
      num_surfatom = num_surfatom + 1
      if(cov_type(i) >= 3) then     ! debug段：表面原子如果被定义为3-5，则一律设为0
        cov_type(i) = 0
        write(*, *) 'Some surface atoms were defined as bulk in ini.xyz! Already put them to empty atoms'
      end if
    end if
  end do

  ! 生成表面原子序号的索引数组
  allocate(surf_index(num_surfatom))
  surf_index = 0
  j = 0
  do i = 1, num_atoms
    if(cov_type(i) < 5) then
      j = j + 1
      surf_index(j) = i
    end if
  end do
  ! debug段：用于检查表面序号索引有无错误
  if(j /= num_surfatom) then
    write(*, *) "Number of surface atoms occurs error! Please check the code"
    stop 
  end if


  allocate(Nevent_persite(num_surfatom))
  Nevent_persite = 0


  ! 计算GCN
  do i = 1, num_atoms
    do j = 1, 12
      temp_int4 = nnsite(j, i)
      if(temp_int4 /= 0) gcn(i) = gcn(i) + cn(temp_int4)
      ! debug段
      ! write(*, *) 'cn(temp_int4):', cn(temp_int4), 'temp_int4', temp_int4
      ! write(*, *) 'gcn(i)', gcn(i)
    end do
    gcn(i) = gcn(i) / 12.0
  end do



  ratio_rej = 0.0

  ! 计算本结构的覆盖度
  call cov_H_CO(num_CO, num_H, cov_CO, cov_H)
  ! 计算初始结构的吸附能
  call calc_Eads(E_now)
  ! 输出first_one.xyz
  open(10, file = 'first_one.xyz', status='replace', action='write')     ! 10为结构轨迹编号
    write(10, *) num_atoms
    ! 步数、吸附能量、CO个数、CO覆盖度、H个数、H覆盖度、最近10w步拒绝比例、每原子变化数的最小值
    write(10, '(A, I10, F24.14, I10, F15.9, I10, F15.9, F15.9, I10, I10)') &
      'first_one', nstep_tot, E_now, num_CO, cov_CO, num_H, cov_H, ratio_rej, &
      minval(Nevent_persite), maxval(Nevent_persite) 
    do i = 1, num_atoms
      if(cov_type(i) == 0) write(10, *) 'Pt', x(i), y(i), z(i), cov_type(i)
    end do 
    do i = 1, num_atoms
      if(cov_type(i) == 1) write(10, *) 'Au', x(i), y(i), z(i), cov_type(i)
    end do
    do i = 1, num_atoms
      if(cov_type(i) == 2) write(10, *) 'Ir', x(i), y(i), z(i), cov_type(i)
    end do
    do i = 1, num_atoms
      if(cov_type(i) == 5) write(10, *) 'Co', x(i), y(i), z(i), cov_type(i)         
    end do
  close(10)


  call init_random_seed(seed_flag)
  n_file = 0
  nstep_tot = 0
  last_steps = 0

  ! 写主体部分的输出文件
  write(char_seq, '(i2.2)') n_file
  ! write(filename_traj, '(A, A, A)') 'atom_str_',trim(char_seq),'.xyz'  ! real atoms coordination
  write(filename_struc, '(A, A, A)') 'struc_rec_',trim(char_seq),'.dat'  !structure information
  write(filename_proc, '(A, A, A)') 'proc_rec_',trim(char_seq),'.dat'  ! process information

  ! open(10, file=filename_traj, status='new', action='write')
  open(11, file=filename_struc, status='new', action='write')
  open(12, file=filename_proc, status='new', action='write')
    call output_traj


    write(*, *) 'GCMC cycle start!'
    do nstep_tot = 1, nLoop

      ! 输出文件间隔控制
      if(nstep_tot > 1 .and. mod(nstep_tot, nLoop_peroutput) == 1) then
        n_file = n_file + 1
        write(char_seq, '(i2.2)') n_file
        ! write(filename_traj, '(A, A, A)') 'atom_str_',trim(char_seq),'.xyz'  !real atoms coordination
        write(filename_struc, '(A, A, A)') 'struc_rec_',trim(char_seq),'.dat'  !structure information
        write(filename_proc, '(A, A, A)') 'proc_rec_',trim(char_seq),'.dat'  ! process information

        ! close(10)
        close(11)
        close(12)

        ! open(10, file=filename_traj, status='new', action='write')
        open(11, file=filename_struc, status='new', action='write')
        open(12, file=filename_proc, status='new', action='write')
      end if

      
      ! 初始化steps记录
      steps%position = 0
      steps%cov_before = 0

      ! 随机决定事件发生位点
      call random_integer(site_seq, 1, num_surfatom)
      steps%position = site_seq
      steps%cov_before = cov_type(surf_index(site_seq))

      ! 该事件位点计数一次
      Nevent_persite(site_seq) = Nevent_persite(site_seq) + 1

      ! 随机决定事件种类（两种事件：将这个原子的cov_type随机更换为另外两种）
      call random_integer(type_event, 1, 2)

      if(cov_type(surf_index(site_seq)) == 0) then       ! 对空原子：换成1 或 换成2 
        if(type_event == 1) then
          cov_type(surf_index(site_seq)) = 1
        else
          cov_type(surf_index(site_seq)) = 2
        end if

      elseif(cov_type(surf_index(site_seq)) == 1) then       ! 对CO吸附原子：换成0 或 换成2 
        if(type_event == 1) then
          cov_type(surf_index(site_seq)) = 0
        else
          cov_type(surf_index(site_seq)) = 2
        end if

      elseif(cov_type(surf_index(site_seq)) == 2) then       ! 对H吸附原子：换成0 或 换成1
        if(type_event == 1) then
          cov_type(surf_index(site_seq)) = 0
        else
          cov_type(surf_index(site_seq)) = 1
        end if

      else     ! debug段：表面原子只能有0, 1, 2三种cov_type
        write(*, *) "Coverage type of surface atom occurs error! Please check the code"
        stop
      end if

      ! 计算本结构的覆盖度
      call cov_H_CO(num_CO, num_H, cov_CO, cov_H)
      ! 计算试探事件后的吸附能
      call calc_Eads(E_test)
      ! debug段
      ! write(*, *) 'E_now:', E_now, 'E_test:', E_test


      ! 判断事件是否接受
      if(E_now > E_test) then    ! 如果E_now > E_test，则事件接受，保持状态更新，更新能量
        accept_flag = 1
        write(12, '(I10, F24.14, F24.14, F24.14, I5, F15.9, F15.9)') nstep_tot, E_now, &
          E_test, (E_test - E_now), accept_flag, 0.0, 0.0

        E_now = E_test
        last_steps(mod(nstep_tot, 100000)) = 0    ! last_steps中，被接受的步记为0

        
      else                       ! 如果E_now <= E_test：

        ! debug段 
        ! write(*, *) 'generating r_accept! E_now - E_test = ', (E_now - E_test)
        !按照能量差的玻尔兹曼因子，生成接受事件的概率
        r_accept = exp((E_now - E_test) / (kb * Temp))
        ! debug段：输出r_accept，如果大于1或小于0则停止代码
        ! write(*, *) 'step:', nstep_tot, 'r_accept: ', r_accept
        if(r_accept > 1 .or. r_accept < 0) then
          write(*, *) 'Error: r_accept = ', r_accept, '! Please check the code'
          stop
        end if

        call random_number(R_random)

        if(R_random < r_accept) then   ! 如果随机数<概率则接受，保持cov状态更新，更新能量
          accept_flag = 1
          write(12, '(I10, F24.14, F24.14, F24.14, I5, F15.9, F15.9)') nstep_tot, E_now, &
            E_test, (E_test - E_now), accept_flag, r_accept, R_random

          E_now = E_test
          last_steps(mod(nstep_tot, 100000)) = 0

        else                           ! 如果随机数<概率，则拒绝事件，回退cov_type更新，回退事件计数，能量不变
          accept_flag = 0
          ! write(12, '(I10, F23.14, F23.14, F23.14, I5, F15.9, F15.9)') nstep_tot, E_now, &
          !   E_test, (E_test - E_now), accept_flag, r_accept, R_random

          cov_type(surf_index(steps%position)) = steps%cov_before
          Nevent_persite(steps%position) = Nevent_persite(steps%position) - 1

          call cov_H_CO(num_CO, num_H, cov_CO, cov_H)      ! 计算回退后的覆盖度
          ! call calc_Eads(E_now)    

          last_steps(mod(nstep_tot, 100000)) = 1    ! last_steps中，被拒绝的步记为1
          
        end if

      end if

      ! 拒绝率计算
      total_rej = 0
      do j = 1, 100000
        total_rej = total_rej + last_steps(j) 
      end do
      ! debug段：输出last_steps总和
      ! write(*, *) 'sum of last_steps:', total_rej
      if(nstep_tot < 100000) then
        ratio_rej = REAL(total_rej) / REAL(nstep_tot)
      else
        ratio_rej = REAL(total_rej) / 100000
      end if


      ! 假如事件被接受，则输出轨迹
      if(accept_flag == 1) then
        call output_traj
      end if


    end do


  ! close(10)
  close(11)
  close(12)

  write (*, *) "====== END ======="

  !输出last_one.xyz，加上体相原子
  open(10, file='last_one.xyz', status='replace', action='write')
    write(10, *) num_atoms
    ! 步数、吸附能量、CO个数、CO覆盖度、H个数、H覆盖度、最近10w步拒绝比例、每原子变化数的最小值
    write(10, '(A, I10, F24.14, I10, F15.9, I10, F15.9, F15.9, I10, I10)') &
      'last_one', nstep_tot, E_now, num_CO, cov_CO, num_H, cov_H, ratio_rej, &
      minval(Nevent_persite), maxval(Nevent_persite) 

    do i = 1, num_atoms
      if(cov_type(i) == 0) write(10, *) 'Pt', x(i), y(i), z(i), cov_type(i)
    end do 
    do i = 1, num_atoms
      if(cov_type(i) == 1) write(10, *) 'Au', x(i), y(i), z(i), cov_type(i)
    end do
    do i = 1, num_atoms
      if(cov_type(i) == 2) write(10, *) 'Ir', x(i), y(i), z(i), cov_type(i)
    end do
    do i = 1, num_atoms
      if(cov_type(i) == 5) write(10, *) 'Co', x(i), y(i), z(i), cov_type(i)         
    end do
  close(10)


  deallocate(element)
  deallocate(x)
  deallocate(y)
  deallocate(z)
  deallocate(cov_type)
  deallocate(cn)
  deallocate(nnsite)
  deallocate(gcn)
  deallocate(surf_index)
  deallocate(Nevent_persite)


contains
  subroutine cov_H_CO(N_CO, N_H, C_CO, C_H)   ! 计算CO、H数量，及覆盖度
    implicit none

    integer(kind=4), intent(out) :: N_CO, N_H
    real(kind=4), intent(out) :: C_CO, C_H

    N_CO = 0
    N_H = 0

    do i = 1, num_atoms
      if(cov_type(i) == 1) then
        N_CO = N_CO + 1
      elseif(cov_type(i) == 2) then
        N_H = N_H + 1
      end if
    end do 

    C_CO = REAL(N_CO) / REAL(num_surfatom)
    C_H = REAL(N_H) / REAL(num_surfatom)
    ! debug段
    ! write(*, *) 'num_CO:', N_CO, 'num_H:', N_H
    ! write(*, *) 'num_surf:', num_surfatom
    ! write(*, *) 'C_CO:', C_CO, 'C_H:', C_H
  end subroutine

  subroutine output_traj
    implicit none

    ! 输出坐标 
    ! write(10, *) num_surfatom
    ! write(10, *) nstep_tot, E_now, num_CO, cov_CO, num_H, cov_H, ratio_rej    ! 吸附能量、CO个数、CO覆盖度、H个数、H覆盖度、拒绝比例（拒绝次数/最近10w步）
    ! do i = 1, num_atoms
    !   if(cov_type(i) == 0) write(10, *) 'Pt', x(i), y(i), z(i), cov_type(i)
    ! end do 
    ! do i = 1, num_atoms
    !   if(cov_type(i) == 1) write(10, *) 'Au', x(i), y(i), z(i), cov_type(i)
    ! end do
    ! do i = 1, num_atoms
    !   if(cov_type(i) == 2) write(10, *) 'Ir', x(i), y(i), z(i), cov_type(i)
    ! end do

    ! 输出能量等
    write(11, '(I10, F24.14, I10, F15.9, I10, F15.9, F15.9, I10, I10)') &
      nstep_tot, E_now, num_CO, cov_CO, num_H, cov_H, ratio_rej, &
      minval(Nevent_persite), maxval(Nevent_persite) 

  end subroutine


  subroutine calc_Eads(E_ads)
    implicit none

    real(kind=8), intent(out) :: E_ads
    real(kind=8) :: E_tot

    E_tot = 0.0

    do i = 1, num_atoms
      if(cov_type(i) < 5) then
        ! MM能量被当作0点

        ! M-AD能量（使用SISSO拟合的吸附能公式）
        if(cov_type(i) == 1) then
          E_tot = E_tot - 0.0245*(cn(i)-gcn(i))*U_SHE*U_SHE - 0.00215*gcn(i)*gcn(i)*U_SHE &
                  + 0.000771*cn(i)*gcn(i)*gcn(i) - 1.103
          ! E_tot = E_tot + (-0.243*gcn(i) + 0.168*cn(i) + 0.168) * exp(U_SHE) &
                  ! + 0.00119*cn(i)*gcn(i)*gcn(i) - 1.17
          ! E_tot = E_tot + M_CO_slope * gcn(i) + M_CO_intercept
          ! debug段
          ! write(*, *) 'gcn(i):', gcn(i)
          ! write(*, *) 'E_tot1:', E_tot
        elseif(cov_type(i) == 2) then ! H使用固定值吸附能
          E_tot = E_tot - 0.1
          ! E_tot = E_tot + 0.892*U_SHE*U_SHE/(cn(i)*gcn(i)) + (0.820/gcn(i)-0.00179*cn(i)*cn(i))*U_SHE &
          !         - 0.140
          ! E_tot = E_tot - 0.0111*U_SHE*U_SHE/(cn(i)-gcn(i)) - 0.204*U_SHE/cn(i) &
          !         + 0.0172*(cn(i)-2*gcn(i)) - 0.0709
          ! E_tot = E_tot + (-0.0112*gcn(i) + 0.0112) * exp(U_SHE) &
                  ! - 0.198*cn(i)/gcn(i) - 0.198*gcn(i) + 0.013*cn(i)*gcn(i) + 0.729
          ! E_tot = E_tot + M_H_slope * gcn(i) + M_H_intercept
          ! debug段
          ! write(*, *) 'E_tot2:', E_tot
        end if

        ! AD-AD能量
        if(cn(i) > 6 .and. cov_type(i) == 1) then     ! CN>6且CO吸附位点：累加紧邻的CO-CO与CO-H排斥
          do j = 1, 12
            if(cov_type(nnsite(j, i)) == 1) E_tot = E_tot + E_COrep
            if(cov_type(nnsite(j, i)) == 2) E_tot = E_tot + E_CO_Hrep
          end do
        elseif(cn(i) > 6 .and. cov_type(i) == 2) then ! CN>6且H吸附位点：累加近邻的H-CO与H-H排斥
          do j = 1, 12
            if(cov_type(nnsite(j, i)) == 1) E_tot = E_tot + E_CO_Hrep
            if(cov_type(nnsite(j, i)) == 2) E_tot = E_tot + E_Hrep
          end do
        end if

      end if
    end do

    ! 计算吸附能
    E_ads = E_tot - num_CO * u_perCO - num_H * u_perH
    ! debug段
    ! write(*, *) 'E_tot:', E_tot, 'num_CO:', num_CO, 'num_H:', num_H, 'u_perCO:', u_perCO, 'u_perH', u_perH, 'E_ads:', E_ads
  end subroutine

end program
