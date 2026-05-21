! ! random cov_type after GCMC


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
        open (1, file='seed-gencov.txt', status='new', action='write')
            write (1, *) n
            write (1, *) clock
        close (1)
    ! seed_flag = 1：use old seed
    elseif(seed_flag == 1) then
        open (1, file='seed-gencov.txt', status='old', action='read')
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
    logical :: flag_is_found

    ! 原子编号计数等变量
    integer(kind=4) :: i, j, temp_int4, num_atoms, num_surfatom, num_H, num_CO
    real(kind=4) :: cov_H, cov_CO, ratio_rej
    ! 坐标
    real(kind=4), allocatable :: x(:), y(:), z(:)
    real(kind=4) :: ra, rx, ry, rz, z_surf
    real(kind=4) :: rx_slab, ry_slab, rz_slab  
    ! 相应的数组：cn、cn_2、近邻、次近邻列表
    integer(kind=4), allocatable :: cn(:), surf_index(:)
    character(len=50) :: char, elem
    ! 步数
    integer(kind=4) :: nstep_tot, min_event, max_event
    ! 能量
    real(kind=8) :: E_tot

    ! 为了随机撒点建随机结构的变量
    integer,allocatable :: indices(:)   ! 存储已选择的索引


    open(1, file='last_one.xyz', status='old', action='read')
        read(1, *) num_atoms
        read(1, *) char, nstep_tot, E_tot, num_CO, cov_CO, num_H, cov_H, ratio_rej, min_event, max_event

        ! 初始化待分配数组
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

        do i = 1, num_atoms
            read(1, *) elem, x(i), y(i), z(i)
        end do
    close(1)

    ! z_surf = (maxval(z)-minval(z)) / 2.0 + minval(z)


    !!! initialize from input files
    open(3, file="input", status="old", action='read')
        read(3, *) seed_flag
        read(3, *) rx_slab, ry_slab, rz_slab
        read(3, *) z_surf
    close(3)

    ! 算整个slab的rx, ry, rz（用于周期性边界条件）
    ! rx_slab = rx_peratom * tot_Natom_x
    ! ry_slab = ry_peratom * tot_Natom_y
    ! rz_slab = rz_peratom * tot_Natom_z

    num_surfatom = 0
    ! 生成CN
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
    write(*, *) num_surfatom  ! debug


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
  
    call init_random_seed(seed_flag)

    ! 将num_CO和num_H随机撒点到表面原子上，输出结构“last_one_random.xyz”
    allocate(indices(num_CO + num_H))
    
    ! 替换为CO
    do i = 1, num_CO
        flag_is_found = .false.
        do while(.not. flag_is_found)
            call random_integer(temp_int4, 1, num_surfatom)
            if(.not. any(indices == temp_int4)) then
                cov_type(surf_index(temp_int4)) = 1
                indices(i) = temp_int4 
                flag_is_found = .true.
            end if
        end do
    end do

    ! 替换为H
    do i = 1, num_H
        flag_is_found = .false.
        do while(.not. flag_is_found)
            call random_integer(temp_int4, 1, num_surfatom)
            if(.not. any(indices == temp_int4)) then
                cov_type(surf_index(temp_int4)) = 2
                indices(i) = temp_int4 
                flag_is_found = .true.
            end if
        end do
    end do

    open(10, file="last_one_random.xyz", status="new", action="write")
        write(10, *) num_atoms
        write(10, '(A, I10, F24.14, I10, F15.9, I10, F15.9, F15.9, I10, I10)') &
            'last_one', nstep_tot, E_tot, num_CO, cov_CO, num_H, cov_H, ratio_rej, &
            min_event, max_event
   
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

    
    deallocate(x)
    deallocate(y)
    deallocate(z)
    deallocate(cov_type)
    deallocate(cn)
    deallocate(surf_index)
    deallocate(indices)
  

end program main