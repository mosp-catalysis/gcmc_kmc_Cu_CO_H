module struct_data
    use outer_data, only: flag_restart, elem, rx_slab, ry_slab, rz_slab, specie
    implicit none
    ! This module is used to read and record structure data

    integer(kind=4), private :: i, j

    integer(kind=4), protected :: nbulk ! number of atoms in bulk
    real(kind=4), protected :: xxx(300000), yyy(300000), zzz(300000)  ! coordinate of the grid
    
    integer(kind=4), protected :: natoms, last_ADS1_cov, last_ADS2_cov ! number of atoms in NPs
    real(kind=4), protected :: x(300000), y(300000), z(300000)  ! coordinate of the NP
    integer(kind=4), protected :: nnsite(12, 300000), nnnsite(55, 300000), atom_num(300000) ! # of CN and SCN
    real(kind=4), protected :: z_surf

    real(kind=8), public :: gcn(300000) ! record generalized coordination number of each site
    integer(kind=4), public :: cn(300000), cn_2(300000), effcn(300000) ! CN(bulk site include), SCN, EFFCN(NP site only)
    integer(kind=4), public :: cov_type(300000) ! 0-none; 1-ADS1; 2-ADS2; 5-bulk
    integer(kind=4), public :: site_type(300000) ! 1-bulk; 2-surface site; 3-outer grid
    integer(kind=4), public :: R_TON(300000), ADS_ADS1_TON(300000), ADS_ADS2_TON(300000)
    
    integer(kind=4), public :: ijpick_cn(2, 12, 300000)
    real(kind=8), public :: ijpick_gcn(2, 12, 300000)


    contains

    subroutine calc_para
        if(specie == 'P') then
            continue    !!! 还没写！
            ! clx = latt_para
            ! cly = latt_para
            ! clz = latt_para
            ! nlx = int(dimx/clx)
            ! nly = int(dimy/cly)
            ! nlz = int(dimz/clz)
            ! bx = (nlx-0.5)*clx
            ! by = (nly-0.5)*cly
            ! bz = (nlz-0.5)*clz
        elseif(specie == 'S') then
            ! clz = latt_para
            ! rx_slab = rx_peratom * nlx
            ! ry_slab = ry_peratom * nly
            ! rz_slab = rz_peratom * nlz          
        end if
    end subroutine calc_para


    subroutine read_str
        last_ADS1_cov = 0
        last_ADS2_cov = 0
        if(flag_restart == 0) then     ! 产生新结构
            if(specie == 'P') call read_str_new_p
            if(specie == 'S') call read_str_new_s
        end if                    ! 读取旧结构
        if(flag_restart == 1) call read_str_last
    end subroutine read_str


    subroutine read_str_new_p
        continue       !!! 还没写！
    end subroutine read_str_new_p


    subroutine read_str_new_s   !!! 还没写！这段函数尚未验证
        real(kind=4), parameter :: Z_INI = 3.96170   ! Total_bulk.xyz里面最小的z
        real(kind=4) :: rz_max

        call calc_para
        rz_max = Z_INI + rz_slab / 2.0 - 0.1    ! 计算新结构的z最大值
        z_surf = rz_max / 2.0
        write(*, *) rz_slab, rz_max, z_surf

        open(2, file="Total_bulk.xyz", status="old", action="read")
            read(2, *) nbulk
            read(2, *)
            natoms = 0
            do i = 1, nbulk
                atom_num(i) = i
                read(2, *) elem, xxx(i), yyy(i), zzz(i)
                site_type(i) = 3
                cov_type(i) = 5
                if(zzz(i) < rz_max) then    ! 产生实际结构
                    natoms = natoms + 1
                    site_type(i) = 1
                    cov_type(i) = 0
                    x(natoms) = xxx(i)
                    y(natoms) = yyy(i)
                    z(natoms) = zzz(i)
                end if
            end do
        close(2)
        ! creat ini file
        open(3, file="ini_new.xyz", status="new", action="write")
            write(3, *) natoms
            write(3, *)
            do i = 1, natoms
                write(3, *) elem, x(i), y(i), z(i), atom_num(i)
            end do
        close(3)
    end subroutine read_str_new_s


    subroutine read_str_last
        real(kind=4) :: dx, dy, dz
        integer(kind=4) :: cov_type_ini(300000)

        call calc_para
        ! read grid
        open(2, file="Total_bulk.xyz", status="old", action="read")
            read(2, *) nbulk
            read(2, *)
            do i = 1, nbulk
                atom_num(i) = i
                read (2, *) elem, xxx(i), yyy(i), zzz(i)
                site_type(i) = 3
                cov_type(i) = 5
            end do
        close(2)
        z_surf = (maxval(zzz)-minval(zzz)) / 4.0 + minval(zzz)

        ! read structure
        open(3, file="ini.xyz", status="old", action="read")
            read(3, *) natoms
            read(3, *)
            do i = 1, natoms
                read (3, *) elem, x(i), y(i), z(i), cov_type_ini(i)
            end do
        close(3)

        ! bulk sites and covs classification
        do i = 1, natoms
            do j = 1, nbulk
                dx = abs(x(i) - xxx(j))
                dy = abs(y(i) - yyy(j))
                dz = abs(z(i) - zzz(j))
                if (dx <= 0.001 .and. dy <= 0.001 .and. dz <= 0.001) then
                    site_type(j) = 1
                    cov_type(j) = cov_type_ini(i)
                    exit
                end if
            end do
            if (j > nbulk) then
                write(*, *) "Error: this atom not match with grid! x=", x(i), "y=", y(i), "z=", z(i)
                stop
            end if
        end do
    end subroutine read_str_last
      

    subroutine count_cn
        ! record cn (nnsite) and scn (nnnsite) of each atom
        real(kind=4) :: rx, ry, rz, dr
        integer(kind=4) :: nn_j

        cn = 0
        cn_2 = 0
        effcn = 0
        nnsite = 0
        nnnsite = 0
        open(4, file="ini_update_cn_gcn.xyz")

        write(4, *) natoms
        write(4, *)
        do i = 1, nbulk
            do j = i + 1, nbulk
                rx = abs(xxx(i) - xxx(j))
                ry = abs(yyy(i) - yyy(j))
                rz = abs(zzz(i) - zzz(j))
                ! extended boundary
                if(rx*2 > rx_slab) rx = rx_slab - rx
                if(ry*2 > ry_slab) ry = ry_slab - ry
                if(rz*2 > rz_slab) rz = rz_slab - rz
        
                dr = sqrt(rx**2 + ry**2 + rz**2)
                ! if ( (i.eq.58 .or. j.eq.58) .and. (dr.gt.2.0 .and. dr.lt.3.5)) then
                !     write(*, *) dr, i, j, xxx(i) - xxx(j), yyy(i) - yyy(j), zzz(i) - zzz(j)
                ! end if

                if(dr < 2.7) then
                    cn(i) = cn(i) + 1
                    cn(j) = cn(j) + 1
                    nnsite(cn(i), i) = j
                    nnsite(cn(j), j) = i
                    ! 用于debug
                    if (dr < 2) write(*, "(F5.3, 2I6)") dr, i, j
                    if (cn(i) > 12 .or. cn(j) > 12) then
                        write(*, *) "Error: too many nn. dr=", dr, "i=", i, "j=", j, "cn(i)=", cn(i), "cn(j)=", cn(j)
                        if (cn(i) > 12) write(*, "(F5.3, 15I6)") dr, i, j, cn(i), nnsite(:, i)
                        if (cn(j) > 12) write(*, "(F5.3, 15I6)") dr, i, j, cn(j), nnsite(:, j)
                        stop
                    end if

                    if(site_type(i) < 3 .and. site_type(j) < 3) then
                        effcn(i) = effcn(i) + 1
                        effcn(j) = effcn(j) + 1
                    end if
                end if

                if(dr < 5.2) then
                    cn_2(i) = cn_2(i) + 1
                    cn_2(j) = cn_2(j) + 1
                    ! 用于debug
                    if (cn_2(i) > 54 .or. cn_2(j) > 54) then
                        write(*, *) "Error: too many s-nn. dr=", dr, "i=", i, "j=", j
                        stop
                    end if
                    ! if (cn_2(i).gt.54) write(*, "(F5.3, 15I6)") dr, i, j, cn_2(i)
                    ! if (cn_2(j).gt.54) write(*, "(F5.3, 15I6)") dr, i, j, cn_2(j)
                    nnnsite(cn_2(i), i) = j
                    nnnsite(cn_2(j), j) = i
                end if
            end do
        end do

        ! calculate and record gcn of each atom
        gcn = 0.0
        do i = 1, nbulk
            if(site_type(i) < 3) then
                do j = 1, 12
                    nn_j = nnsite(j, i)
                    gcn(i) = gcn(i) + effcn(nn_j)
                end do
                gcn(i) = gcn(i) / 12.0
                write(4, *) elem, xxx(i), yyy(i), zzz(i), atom_num(i), effcn(i), gcn(i), cov_type(i)
            end if
        end do

        close(4)
    end subroutine count_cn


    subroutine define_site   !!!还没写！NP部分未经验证
        if(specie == 'P') then   ! 对于NP：表面原子定义site为2，新生成的ini结构则再定义cov为0
            do i = 1, nbulk
                if(effcn(i) > 0 .and. effcn(i) < 12) then
                    site_type(i) = 2
                    if(flag_restart == 0) cov_type(i) = 0
                end if
            end do
        else                     ! 对于slab：所有effcn=12的设置为site1，cov5。大于z_surf的原子再区分表面出来，小于z_surf的原子则不改
            do i = 1, nbulk
                if(site_type(i) < 3) then
                    if(zzz(i) > z_surf) then
                        if(effcn(i) >= 0 .and. effcn(i) < 12) then
                            site_type(i) = 2
                            if(flag_restart == 0) cov_type(i) = 0
                        elseif(effcn(i) == 12) then
                            site_type(i) = 1
                            cov_type(i) = 5
                        end if
                    else
                        site_type(i) = 1
                        cov_type(i) = 5
                    end if
                end if
            end do
        end if
    end subroutine define_site
end module struct_data
