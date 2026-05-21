! ! KMC code, by myself!
! ! by Shuoqi Zhang, written in 2024.10.12-
! ! hope it can work :>

program main
    use outer_data
    use struct_data
    use Array_oper
    implicit none

    integer,parameter :: r400 = selected_real_kind(r=400)

    integer(kind=4), parameter :: N_SINGLE_EVENT = 4 ! number of single-site-process
    integer(kind=4), DIMENSION(N_SINGLE_EVENT) :: single_event = (/2, 3, 4, 5/) ! # of single-site-process

    integer(kind=8), parameter :: refresh_int = 1000000 ! inteval to update rtot, to avoid error accumulation

    character(len=120) :: filename1, filename2, filename3, char_seq, fmt_string1, fmt_string2
    character(len=120) :: filename4, filename5, filename6

    integer(kind=4) :: i, j, k, r_flag, flag_record   ! r_flag：记录结构的flag（用于防止每步aj都输出结构时取余全部为0）。flag_record：记录step信息的flag（以防重复）

    real(kind=r400) :: r_site, rpoint, rtot, rsite(300000), revent(N_EVENT_TYPE, 300000), rneis(12, N_EVENT_TYPE, 300000)
    real(kind=r400) :: revent_tot(N_EVENT_TYPE)
    real(kind=8) :: Ea_r(2, 12, 300000)

    integer(kind=4) :: j_cn, i_nn ! # of nnsite 
    integer(kind=4) :: ipick, jpick, kpick, nn_i, nn_j, nn_k, nf_1, neis_i, neis_j, cov_old
    integer(kind=8) :: n_step_tot, n_step_each(N_EVENT_TYPE), ii_step
    integer(kind=4) :: natoms_ij(108), ni_effcn, i_nADS1, i_nADS2, j_nADS1, j_nADS2
    ! integer(kind=8) :: r_last_step(100000) ! last step of coverage state changes
    real(kind=8) :: gcn_ni
    real(kind=8) :: randum
    real(kind=r400) :: ctime, dtime, sum_rsite, sum_revent, sum_rneis
    integer(kind=4) :: update_list_i(108)

    ! get parameter from 'input'
    call read_data 
    ! generate grid and read initial structure
    call read_str
    ! calculate cn, cnn, gcn
    call count_cn
    ! define surface atom
    call define_site

    ! calculate rtot, rsite, rneis, revent
    rtot = 0.0
    rsite = 0.0
    revent = 0.0
    rneis = 0.0
    revent_tot = 0.0
    Ea_r = 0.0
    r_flag = 1 ! flag to do structure record

    do i = 1, nbulk
        if(site_type(i) == 2) then
            do k = 1, N_EVENT_TYPE
                if (enable_evt(k)) then
                    if(any(k == single_event)) then   ! 单位点事件
                        r_site = 0
                        call rijk(r_site, i, k, 0, Ea_r)
                        revent(k, i) = r_site
                    else    ! 双位点事件，再迭代j
                        do j = 1, 12
                            r_site = 0
                            call rijk(r_site, i, k, j, Ea_r)
                            rneis(j, k, i) = r_site
                            revent(k, i) = revent(k, i) + r_site
                        end do
                    end if
                    revent_tot(k) = revent_tot(k) + revent(k, i) 
                    rsite(i) = rsite(i) + revent(k, i)
                end if
            end do
            rtot = rtot + rsite(i)
        end if
    end do
    write(*, *) 'inital rtot=', rtot

    call init_random_seed(flag_reseed)
    n_step_tot = 0
    n_step_each = 0
    ii_step = 0
    ctime = 0.0
    ! new files
    nf_1 = 0
    write(char_seq,'(i2.2)') nf_1
    write(filename1, '(3A)') 'atom_inf_',trim(char_seq),'.dat'  ! atoms coordination information (site_type cov_type) in grid space
    ! write(filename2, '(3A)') 'cn_gcn_', trim(char_seq), '.dat'  ! gcn & cn data (in real atomistic structure)
    write(filename3, '(3A)') 'Ea_react_',trim(char_seq),'.dat'   ! ctime, n_tot_step of reaction
    write(filename4, '(3A)') 'Ea_atom_',trim(char_seq),'.dat'   ! energy barrier of atom jumping, cn, n_tot_step
    write(filename5, '(3A)') 'step_rec_',trim(char_seq),'.dat'  ! steps of every event
    write(filename6, '(3A)') 'step_rate_',trim(char_seq),'.dat'   

    open(11, file=filename1, status='replace', action='write')
    ! open(12, file=filename2, status='replace', action='write')
    open(13, file=filename3, status='replace', action='write')
    open(14, file=filename4, status='replace', action='write')
    open(15, file=filename5, status='replace', action='write')
    open(16, file=filename6, status='replace', action='write')

    ! record first structural coordination
    k = 0
    write(11, '(I10, I10, e22.12)') ii_step, n_step_tot, ctime
    do i = 1, nbulk
        if(site_type(i) < 3) then
            write(11, '(I10, I5, I5, F10.4)') atom_num(i), cov_type(i), effcn(i), gcn(i)
            if(cov_type(i) < 5) k = k + 1
        end if
    end do
    write(fmt_string1, '(A, I0, A)') "(e22.12, 2I13, ", N_EVENT_TYPE, "I13)"       ! 用于生成输出n_step_each格式的字符串
    write(fmt_string2, '(A, I0, A)') "(e22.12, I13, e22.12, ", N_EVENT_TYPE, "e22.12)"
    write(15, fmt_string1) ctime, k, n_step_tot, n_step_each


    ! KMC cycle
    write(*, *) 'KMC cycle starts!'
    do while (n_step_tot < nLoop)
    ! do while (ctime <= tend)
        ! Circular recording
        if(n_step_tot > 0 .and. mod(n_step_tot, 20000000) == 0) then
            nf_1 = nf_1 + 1
            write(char_seq,'(i2.2)') nf_1
            write(filename1, '(3A)') 'atom_inf_',trim(char_seq),'.dat'  ! atoms coordination information (site_type cov_type) in grid space
            ! write(filename2, '(3A)') 'cn_gcn_', trim(char_seq), '.dat'  ! gcn & cn data (in real atomistic structure)
            write(filename3, '(3A)') 'Ea_react_',trim(char_seq),'.dat'   ! ctime, n_tot_step of reaction
            write(filename4, '(3A)') 'Ea_atom_',trim(char_seq),'.dat'   ! energy barrier of atom jumping, cn, n_tot_step
            write(filename5, '(3A)') 'step_rec_',trim(char_seq),'.dat'  ! steps of every event
            write(filename6, '(3A)') 'step_rate_',trim(char_seq),'.dat'   
        
            open(11, file=filename1, status='replace', action='write')
            ! open(12, file=filename2, status='replace', action='write')
            open(13, file=filename3, status='replace', action='write')
            open(14, file=filename4, status='replace', action='write')
            open(15, file=filename5, status='replace', action='write')
            open(16, file=filename6, status='replace', action='write')
        end if

        ! refresh rtot to avoid error accumulate
        if(n_step_tot > 0 .and. mod(n_step_tot, refresh_int) == 0) then
            sum_rsite = 0.0
            do i = 1, nbulk
                if(site_type(i) == 2) then
                    sum_rsite = sum_rsite + rsite(i)
                end if
            end do
            write(*, *) n_step_tot, ' error:', rtot-sum_rsite
            rtot = sum_rsite
        end if

        ! 进行一次判断：当rtot小于1E-6时，认为无法发生事件，整个程序不再运行
        if (rtot < 1E-6) then
            write(*, *) "Error: rtot is too small! rtot=", rtot
            stop
        end if

        ! cal time
        call random_number(randum)
        dtime = -dlog(randum)/rtot
        ctime = ctime + dtime

        ! site pick
        call random_number(randum)
        rpoint = rtot * randum
        sum_rsite = 0.0
        do i = 1, nbulk
            if (site_type(i) == 2) then
                sum_rsite = sum_rsite + rsite(i)
                if (sum_rsite >= rpoint) exit
            end if
        end do
        ! 如果误差累积太大，则重新按新的rtot找ipick
        if (i > nbulk) then
            write(*, *) 'Warning: error accumulate:', rtot-sum_rsite
            rtot = sum_rsite
            rpoint = rtot * randum
            do i = 1, nbulk
                if(site_type(i) == 2) then
                    sum_rsite = sum_rsite + rsite(i)
                    if(sum_rsite >= rpoint) exit
                end if
            end do
        end if
        ipick = i

        ! event pick
        call random_number(randum)
        rpoint = rsite(ipick) * randum
        sum_revent = 0.0
        do k = 1, N_EVENT_TYPE
            if (enable_evt(k)) then
                sum_revent = sum_revent + revent(k, ipick)
                if (sum_revent >= rpoint) exit
            end if
        end do
        kpick = k
    
        ! j_cn & jpick
        if(.not. any(kpick == single_event)) then
            call random_number(randum)
            rpoint = revent(kpick, ipick) * randum
            sum_rneis = 0.0
            do j = 1, 12
                sum_rneis = sum_rneis + rneis(j, kpick, ipick)
                if (sum_rneis >= rpoint) exit
            end do
            j_cn = j
            jpick = nnsite(j_cn, ipick)
        end if
      
        n_step_tot = n_step_tot + 1

        SELECT CASE(kpick)
            CASE(1) ! atom jumping
                r_flag = 1
                cov_old = cov_type(ipick)
                site_type(jpick) = 2
                cov_type(jpick) = cov_type(ipick)
                site_type(ipick) = 3
                cov_type(ipick) = 5
                ! find all atoms of 2 neis of ipick and jpick and remove repeat
                natoms_ij = 0
                call Merge_snn(natoms_ij, nnnsite(:, ipick), nnnsite(:, jpick))
                ! update effcn and gcn
                do i = 1, 108
                    nn_i = natoms_ij(i)
                    if(nn_i > 0) then
                        effcn(nn_i) = 0
                        gcn(nn_i) = 0
                        if(site_type(nn_i) < 3) then
                            ni_effcn = 0
                            gcn_ni = 0.0
                            do j = 1, 12
                                nn_j = nnsite(j, nn_i)
                                if(site_type(nn_j) < 3) then
                                    ni_effcn = ni_effcn + 1
                                    do k = 1, 12
                                        nn_k = nnsite(k, nn_j)
                                        if(site_type(nn_k) < 3) gcn_ni = gcn_ni + 1.0
                                    end do
                                end if
                            end do
                            effcn(nn_i) = ni_effcn
                            gcn(nn_i) = gcn_ni / 12.0
                            ! whether there is adsorbate in the Site where cn >= maxCN_ADS
                            if(cov_type(nn_i) == 1 .and. ni_effcn >= maxCN_ADS) then
                                cov_type(nn_i) = 0
                                n_step_tot = n_step_tot + 1
                                n_step_each(3) = n_step_each(3) + 1
                            elseif(cov_type(nn_i) == 2 .and. ni_effcn >= maxCN_ADS) then
                                cov_type(nn_i) = 0
                                n_step_tot = n_step_tot + 1
                                n_step_each(5) = n_step_each(5) + 1
                            ! debug段
                            elseif(ni_effcn > 12) then
                                write(*, *) "Error: too many nn. ni=", nn_i, "ni_effcn", ni_effcn
                                stop
                            end if
                        end if
                    end if
                end do

                ! update site_type and cov_type
                do i = 1, 12
                    ! ipick, body -> surface
                    neis_i = nnsite(i, ipick)
                    if(site_type(neis_i) == 1 .and. effcn(neis_i) < 12) then
                        site_type(neis_i) = 2
                        cov_type(neis_i) = 0
                    end if
                    ! jpick, surface -> body
                    neis_j = nnsite(i, jpick)
                    if(site_type(neis_j) == 2 .and. effcn(neis_j) == 12) then
                        site_type(neis_j) = 1
                        cov_type(neis_j) = 5
                    end if
                end do
                ii_step = ii_step + 1
                n_step_each(1) = n_step_each(1) + 1

            CASE(2) ! ADS1 adsorption
                cov_type(ipick) = 1
                n_step_each(2) = n_step_each(2) + 1
                ADS_ADS1_TON(ipick) = ADS_ADS1_TON(ipick) + 1

            CASE(3) ! ADS1 desorption
                cov_type(ipick) = 0
                n_step_each(3) = n_step_each(3) + 1
                ! update the steps of this site covering by CO

            CASE(4) ! ADS2 adsorption
                cov_type(ipick) = 2
                n_step_each(4) = n_step_each(4) + 1
                ADS_ADS2_TON(ipick) = ADS_ADS2_TON(ipick) + 1
                
            CASE(5) ! ADS2 desorption
                cov_type(ipick) = 0
                n_step_each(5) = n_step_each(5) + 1

            CASE(6) ! ADS1 diffusion
                cov_type(ipick) = 0
                cov_type(jpick) = 1
                n_step_each(6) = n_step_each(6) + 1

            CASE(7) ! ADS2 diffusion
                cov_type(ipick) = 0
                cov_type(jpick) = 2
                n_step_each(7) = n_step_each(7) + 1

            CASE(8) ! reaction
                cov_type(ipick) = 0
                cov_type(jpick) = 0
                n_step_each(8) = n_step_each(8) + 1
                ! count TON
                R_TON(ipick) = R_TON(ipick) + 1
                R_TON(jpick) = R_TON(jpick) + 1
        END SELECT

        ! record
        i_nADS1 = 0
        i_nADS2 = 0
        j_nADS1 = 0
        j_nADS2 = 0
        if(kpick == 1) then  ! 原子跳事件：记录i、j位置近邻的吸附分子数量，记录到Ea文件里 
            do j = 1, 12
                nn_i = cov_type(nnsite(j, ipick))
                nn_j = cov_type(nnsite(j, jpick))
                if(nn_i == 1) i_nADS1 = i_nADS1 + 1
                if(nn_i == 2) i_nADS2 = i_nADS2 + 1
                if(nn_j == 1) j_nADS1 = j_nADS1 + 1
                if(nn_j == 2) j_nADS2 = j_nADS2 + 1
            end do
            write(14, "(f15.6, I5, 2I7, 2f11.4, 4I7, I13, e22.12)") Ea_r(1, j_cn, ipick), cov_old, &
            ijpick_cn(1, j_cn, ipick), ijpick_cn(2, j_cn, ipick), &
            ijpick_gcn(1, j_cn, ipick), ijpick_gcn(2, j_cn, ipick), &
            i_nADS1, i_nADS2, j_nADS1, j_nADS2, n_step_tot, ctime
        elseif(kpick == 8) then  ! 反应事件：记录i、j位置近邻的吸附分子数量，记录到ijk文件里 
            do j = 1, 12
                nn_i = cov_type(nnsite(j, ipick))
                nn_j = cov_type(nnsite(j, jpick))
                if(nn_i == 1) i_nADS1 = i_nADS1 + 1
                if(nn_i == 2) i_nADS2 = i_nADS2 + 1
                if(nn_j == 1) j_nADS1 = j_nADS1 + 1
                if(nn_j == 2) j_nADS2 = j_nADS2 + 1
            end do
            write(13, "(f15.6, 2I5, 2f11.4, 4I7, 2I13, e22.12)") &
            Ea_r(2, j_cn, ipick), effcn(ipick), effcn(jpick), gcn(ipick), gcn(jpick), &
            i_nADS1, i_nADS2, j_nADS1, j_nADS2, n_step_tot, n_step_each(1), ctime
        end if

        flag_record = 1
        if(mod(n_step_tot, record_int) == 0) then
            k = 0
            do i = 1, nbulk    ! 统计表面原子数量
                if(site_type(i) < 3 .and. cov_type(i) < 5) k = k + 1
            end do
            write(15, fmt_string1) ctime, k, n_step_tot, n_step_each
            write(16, fmt_string2) ctime, n_step_tot, rtot, revent_tot
            flag_record = 0
        end if

        if(r_flag == 1 .and. ii_step > 0 .and. mod(ii_step, struct_int) == 0) then   ! 记录结构信息
            j = 0
            k = 0
            write(11, '(I10, I10, e22.12)') ii_step, n_step_tot, ctime
            do i = 1, nbulk
                if(site_type(i) < 3) then
                    write(11, '(I10, I5, I5, F10.4)') atom_num(i), cov_type(i), effcn(i), gcn(i)
                    j = j + 1
                    if(cov_type(i) < 5) k = k + 1
                end if
            end do
            ! debug段
            if (j /= natoms) then
                write(*, *) "Warning: number of atoms changed!", n_step_tot, j, ipick, jpick
            end if
            if(flag_record == 1) then
                write(15, fmt_string1) ctime, k, n_step_tot, n_step_each
                write(16, fmt_string2) ctime, n_step_tot, rtot, revent_tot
            end if
            r_flag = 0
        end if


        ! initial and cal new relative r
        ! update rsite, revent, and rneis of all atoms in nnnsite of ipick/i+jpick
        if(any(k == single_event)) then   ! 如果是单位点事件，则把ipick及次近邻记录到一个数组
            update_list_i = 0
            update_list_i(1) = ipick
            update_list_i(2:56) = nnnsite(:, ipick)
        else                              ! 如果是双位点事件，则把ipick、jpick及它们的次近邻（去重后）记录到一个数组
            update_list_i = 0
            call Merge_snn(update_list_i, nnnsite(:, ipick), nnnsite(:, jpick))
        end if

        do i_nn = 1, 108
            i = update_list_i(i_nn)
            if(i /= 0) then
                rtot = rtot - rsite(i)
                rsite(i) = 0
                do k = 1, N_EVENT_TYPE
                    if(enable_evt(k)) then
                        ! initial
                        revent_tot(k) = revent_tot(k) - revent(k, i) 
                        revent(k, i) = 0
                        if(.not. any(k == single_event)) then   ! 双位点事件
                            do j = 1, 12
                                rneis(j, k, i) = 0
                            end do
                        end if
                        ! cal new r
                        if(site_type(i) == 2) then
                            if(any(k == single_event)) then    ! 单位点事件
                                r_site = 0
                                call rijk(r_site, i, k, 0, Ea_r)
                                revent(k, i) = revent(k, i) + r_site
                            else   ! 双位点事件
                                do j = 1, 12
                                    r_site = 0
                                    call rijk(r_site, i, k, j, Ea_r)
                                    rneis(j, k, i) = r_site
                                    revent(k, i) = revent(k, i) + r_site
                                end do
                            end if
                            rsite(i) = rsite(i) + revent(k, i)
                            revent_tot(k) = revent_tot(k) + revent(k, i)                
                        end if
                    end if
                end do
                rtot = rtot + rsite(i)
            end if
        end do
    end do   ! 结束总步数循环

    ! record final coordination information & step information
    k = 0
    write(11, '(I10, I10, e22.12)') ii_step, n_step_tot, ctime
    do i = 1, nbulk
        if(site_type(i) < 3) then
            write(11, '(I10, I5, I5, F10.4)') atom_num(i), cov_type(i), effcn(i), gcn(i)
            if(cov_type(i) < 5) k = k + 1
        end if
    end do
    write(15, fmt_string1) ctime, k, n_step_tot, n_step_each
    write(16, fmt_string2) ctime, n_step_tot, rtot, revent_tot

    close(11)
    ! close(12)
    close(13)
    close(14)
    close(15)
    close(16)

    ! write last_one (xyz & cov_type)
    open(4, file = "last_one.xyz")
    write(4, *) natoms
    write(4, *) "last_one.xyz", ii_step, n_step_tot, ctime
    do i = 1, nbulk
        if (site_type(i) < 3) then
            if (cov_type(i) == 0) then
                write(4, "(A5, 3f12.6, I5, 3I10)") 'Pt', xxx(i), yyy(i), zzz(i), cov_type(i), R_TON(i)
            elseif(cov_type(i) == 1) then
                write(4, "(A5, 3f12.6, I5, 3I10)") 'Au', xxx(i), yyy(i), zzz(i), cov_type(i), R_TON(i)
            elseif(cov_type(i) == 2) then
                write(4, "(A5, 3f12.6, I5, 3I10)") 'Ir', xxx(i), yyy(i), zzz(i), cov_type(i), R_TON(i)
            elseif(cov_type(i) == 5) then
                write(4, "(A5, 3f12.6, I5, 3I10)") 'Co', xxx(i), yyy(i), zzz(i), cov_type(i), R_TON(i)
            end if
        end if
    end do
    close(4)

    write(*, *) "====== END ======="
end program main



subroutine rijk(r_site, ipick, kpick, j_cn, Ea_r)
    use outer_data
    use struct_data
    
    implicit none

    integer,parameter :: r400 = selected_real_kind(r=400)
    real(kind=r400), intent(OUT) :: r_site
    integer(kind=4), intent(IN):: ipick, kpick, j_cn
    real(kind=8), intent(INOUT) :: Ea_r(2, 12, 300000)

    integer(kind=4) :: ni_ads, nj_ads
    integer(kind=4) :: jpick
    real(kind=8) :: rADS1_ads, rADS2_ads
    real(kind=8) :: gcn_ipick, gcn_jpick

    real(kind=8), parameter :: pi = 3.141592654
    real(kind=8), parameter :: kb = 8.6173324D-05
    real(kind=8), parameter :: h = 4.1356676D-15
    real(kind=8), parameter :: eV2J = 1.60217662D-19  ! 1 eV = 1.60217662D-19 J
    real(kind=8), parameter :: Na = 6.0221409D23
    real(kind=8), parameter :: Asite = (10D-10)**2
    real(kind=8), parameter :: p0 = 100000  ! p0 = 100kPa
    
    real(kind=8) :: s0ADS1, s0ADS2, mADS1, mADS2, SADS1_0, SADS2_0
    real(kind=8) :: Eads_i, Eads_j, dE, Ea
    real(kind=8) :: dS_ADS1, dS_ADS2, r_K_eq
    real(kind=4) :: pADS1, pADS2

    ! mADS1 = 28.01D-3/Na   ! kg/atom
    ! mADS2  = 15.999D-3/Na
    ! SADS1_0 = 85.142*(Temp**0.14709)/(Na*eV2J)
    ! SADS2_0 = 89.655*(Temp**0.14489)/(Na*eV2J)
    
    ! pADS1 = ppADS1 * p_tot
    ! pADS2 = ppADS2 * p_tot

    ! 记录i位点原本的信息
    ni_ads = cov_type(ipick)
    gcn_ipick = gcn(ipick)
    if(j_cn /= 0) then     ! 双位点事件，记录j位点原本的信息
        jpick = nnsite(j_cn, ipick)
        nj_ads = cov_type(jpick)
        gcn_jpick = gcn(jpick)
    end if
    
    SELECT CASE (kpick)
        CASE(1) ! atom jumping
            if(site_type(jpick) == 3) then
                call atom_jump(r_site, Ea_r, ipick, j_cn)
            end if
        
        ! CASE(2) ! ADS1 adsorption   !!!还没写！未经验证
        !     if(ni_ads == 0) then
        !         s0ADS1 = s0_ADS1_facet
        !         if(gcn_ipick >= 0 .and. gcn_ipick < 5.33) s0ADS1 = s0_ADS1_edge
        !         if(effcn(ipick) >= 10) s0ADS1 = 0.0
        !         r_site = (s0ADS1 * pADS1 * Asite)/sqrt(2. * pi * mADS1 * kb * eV2J * Temp)
        !     end if
    
        ! CASE(3) ! ADS1 desorption   !!!还没写！未经验证
        !     if(ni_ads == 1) then
        !         s0ADS1 = s0_ADS1_facet        
        !         ! CO ads
        !         if(gcn_ipick >= 0 .and. gcn_ipick < 5.33) s0ADS1 = s0_ADS1_edge
        !         if(effcn(ipick) >= 10) s0ADS1 = 0.0
        !         rADS1_ads = (s0ADS1 * pADS1 * Asite)/sqrt(2. * pi * mADS1 * kb * eV2J * Temp)
        !         ! CO equilibrium
        !         dS_ADS1 = 0 - (SADS1_0 - kb*log(pADS1/p0))
        !         Eads_i = 0.0
        !         call Eads_site(Eads_i, ipick, 1)
        !         r_K_eq = exp(-(Eads_i - Temp*dS_ADS1)/(kb*Temp))
        !         ! r of CO des
        !         r_site = rADS1_ads/(pADS1*r_K_eq)
        !     end if
            
        ! CASE(4) ! ADS2 adsorption   !!!还没写！未经验证
        !     if (ni_ads == 0 .and. nj_ads == 0) then
        !         s0ADS2 = s0_ADS2_facet
        !         if(gcn_ipick.ge. gcn_jpick) then
        !             if(gcn_ipick.ge.0 .and. gcn_ipick.lt.5.33) s0ADS2 = s0_ADS2_edge
        !             if(effcn(ipick).ge.10) s0ADS2 = 0.0
        !         else
        !             if(gcn_jpick.ge.0 .and. gcn_jpick.lt.5.33) s0ADS2 = s0_ADS2_edge
        !             if(effcn(jpick).ge.10) s0ADS2 = 0.0
        !         end if
        !         r_site = (s0ADS2 * pADS2 * Asite)/sqrt(2. * pi * mADS2 * kb * eV2J * Temp)
        !     end if

        ! CASE(5) ! ADS2 desorption   !!!还没写！未经验证
        !     if(ni_ads.eq.2 .and. nj_ads.eq.2) then
        !         ! O2 ads
        !         s0ADS2 = s0_ADS2_facet
        !         if(gcn_ipick.ge. gcn_jpick) then
        !             if(gcn_ipick.ge.0 .and. gcn_ipick.lt.5.33) s0ADS2 = s0_ADS2_edge
        !             if(effcn(ipick).ge.10) s0ADS2 = 0.0
        !         else
        !             if(gcn_jpick.ge.0 .and. gcn_jpick.lt.5.33) s0ADS2 = s0_ADS2_edge
        !             if(effcn(jpick).ge.10) s0ADS2 = 0.0
        !         end if
        !         rADS2_ads = (s0ADS2 * pADS2 * Asite)/sqrt(2. * pi * mADS2 * kb * eV2J * Temp)
        !         ! O2 equilibrium
        !         dS_ADS2 = 0 - (SADS2_0 - kb*log(pADS2/p0))
        !         Eads_i = 0.0
        !         call Eads_site(Eads_i, ipick, 1)
        !         Eads_j = 0.0
        !         call Eads_site(Eads_j, jpick, 1)
        !         r_K_eq = exp(-(Eads_i + Eads_j - Temp*dS_ADS2)/(kb*Temp))
        !         ! if(Eads_i+Eads_j .gt. 0) r_K_eq = exp(-(0 - Temp*dS_O2)/(kb*Temp))
        !         ! O2 des
        !         r_site = rADS2_ads/(pADS2*r_K_eq)
        !     end if
            
        CASE(6) ! ADS1 diffusion
            if(ni_ads == 1 .and. nj_ads == 0 .and. effcn(jpick) < maxCN_ADS) then
                Eads_i = 0.0
                call Eads_site(Eads_i, ipick, 1)
                cov_type(ipick) = 0
                cov_type(jpick) = 1
                Eads_j = 0.0
                call Eads_site(Eads_j, jpick, 1)
                dE = Eads_j - Eads_i
                if(dE <= 0) dE = 0.0
                Ea = dE + Ediff_ADS1
                r_site = (kb*Temp/h)*exp(-Ea/(kb*Temp))
                ! recovery
                cov_type(ipick) = 1
                cov_type(jpick) = 0
            end if

        CASE(7) ! ADS2 diffusion
            if(ni_ads == 2 .and. nj_ads == 0 .and. effcn(jpick) < maxCN_ADS) then
                Eads_i = 0.0
                call Eads_site(Eads_i, ipick, 1)
                cov_type(ipick) = 0
                cov_type(jpick) = 2
                Eads_j = 0.0
                call Eads_site(Eads_j, jpick, 1)
                dE = Eads_j - Eads_i
                if(dE <= 0) dE = 0.0
                Ea = dE + Ediff_ADS2
                r_site = (kb*Temp/h)*exp(-Ea/(kb*Temp))
                ! recovery
                cov_type(ipick) = 2
                cov_type(jpick) = 0
            end if

        CASE(8) ! reaction 
            if((ni_ads == 1 .and. nj_ads == 2) .or. (ni_ads == 2 .and. nj_ads == 1)) then
                Eads_i = 0.0
                call Eads_site(Eads_i, ipick, 1)
                Eads_j = 0.0
                call Eads_site(Eads_j, jpick, 1)
                Ea = BEP1_a*(Eads_i + Eads_j) + BEP1_b
                if(Ea <= 0) Ea = 0.0
                Ea_r(2, j_cn, ipick) = Ea
                r_site = (kb*Temp/h)*exp(-Ea/(kb*Temp))
            end if        
    END SELECT
end subroutine

subroutine Eads_site(Eads_s, site, num)
    use outer_data
    use struct_data
    implicit none

    integer(kind=4), intent(IN) :: site, num 
    real(kind=8), intent(OUT) :: Eads_s

    integer(kind=4) :: i, j
    integer(kind=4) :: nn_s, nn_j
    integer(kind=4) :: n_ADS1, n_ADS2, effcn_ni, effcn_nj

    n_ADS1 = 0
    n_ADS2 = 0

    if(num == 1) then
        effcn_ni = effcn(site)
        if(effcn_ni > 6) then    !统计周围排斥
            do i = 1, 12
                nn_s = nnsite(i, site)
                if(effcn(nn_s) > 6) then
                    if(cov_type(nn_s) == 1) n_ADS1 = n_ADS1 + 1
                    if(cov_type(nn_s) == 2) n_ADS2 = n_ADS2 + 1
                end if
            end do
        end if
    else if(num == 2) then
        effcn_ni = 0
        do i = 1, 12
            nn_s = nnsite(i, site)
            if(site_type(nn_s) < 3) effcn_ni = effcn_ni + 1
        end do
        if(effcn_ni > 6) then
            do i = 1, 12
                nn_s = nnsite(i, site)
                if(cov_type(nn_s) == 1 .or. cov_type(nn_s) == 2) then
                    effcn_nj = 0
                    do j = 1, 12
                        nn_j = nnsite(j, nn_s)
                        if(site_type(nn_j) < 3) effcn_nj = effcn_nj + 1
                    end do
                    if(effcn_nj > 6) then
                        if(cov_type(nn_s) == 1) n_ADS1 = n_ADS1 + 1
                        if(cov_type(nn_s) == 2) n_ADS2 = n_ADS2 + 1                      
                    end if
                end if
            end do
        end if
    end if
    ! 换用SISSO拟合的吸附能公式
    if(cov_type(site) == 1) Eads_s = -0.0245*(effcn_ni-gcn(site))*U_SHE*U_SHE &
        - 0.00215*gcn(site)*gcn(site)*U_SHE + 0.000771*effcn_ni*gcn(site)*gcn(site) - 1.103 &
        + (n_ADS1 * EADS1_ADS1 + n_ADS2 * EADS1_ADS2)
    ! H使用固定吸附能
    if(cov_type(site) == 2) Eads_s = -0.1 + (n_ADS1 * EADS1_ADS2 + n_ADS2 * EADS2_ADS2)
    ! if(cov_type(site) == 2) Eads_s = (0.246/(effcn_ni-gcn(site)) &
    !     - 0.00307*effcn_ni*gcn(site))*U_SHE - 0.279*effcn_ni/gcn(site) + 0.259 &
    !     + (n_ADS1 * EADS1_ADS2 + n_ADS2 * EADS2_ADS2)
    ! if(cov_type(site) == 1) Eads_s = (-0.243*gcn(site) + 0.168*effcn_ni + 0.168) * exp(U_SHE) &
    !     + 0.00119*effcn_ni*gcn(site)*gcn(site) - 1.17 &
    !     + (n_ADS1 * EADS1_ADS1 + n_ADS2 * EADS1_ADS2)
    ! if(cov_type(site) == 2) Eads_s = (-0.0112*gcn(site) + 0.0112) * exp(U_SHE) &
    !     - 0.198*effcn_ni/gcn(site) - 0.198*gcn(site) + 0.013*effcn_ni*gcn(site) + 0.729 &
    !     + (n_ADS1 * EADS1_ADS2 + n_ADS2 * EADS2_ADS2)
        
    ! if(cov_type(site) == 1) Eads_s = (Eads_ADS1_a * gcn(site) + Eads_ADS1_b) + (n_ADS1 * EADS1_ADS1 + n_ADS2 * EADS1_ADS2)
    ! if(cov_type(site) == 2) Eads_s = (Eads_ADS2_a * gcn(site) + Eads_ADS2_b) + (n_ADS1 * EADS1_ADS2 + n_ADS2 * EADS2_ADS2)
end subroutine

subroutine atom_jump(r_site, Ea_r, ipick, j_cn)
    use outer_data, only : Temp, U_SHE, maxCN_ADS, E_Mbond, E_Mcohe, corr_A1, corr_t1, corr_A2, corr_t2
    use struct_data
    use Array_oper, only : Merge_snn
    implicit none
    
    integer,parameter :: r400 = selected_real_kind(r=400)
    integer(kind=4), intent(IN) :: ipick, j_cn
    real(kind=r400), intent(OUT) :: r_site
    real(kind=8), intent(OUT) :: Ea_r(2, 12, 300000)
    
    integer(kind=4) :: i, j, k, jpick, ni_ads, cov_old
    integer(kind=4) :: natoms_ij(108)
    integer(kind=4) :: nn_i, neis_i, neis_j, effcn_jpick
    real(kind=8) :: Ea
    real(kind=8) :: U_ini, U_mid, U_fin, Eb_ini, Eb_fin, dE, Eb, Ef
    real(kind=8) :: Eads_i, Eads_ini, Eads_mid, Eads_fin
    real(kind=8) :: gcn_ni, gcn_old
    
    real(kind=8), parameter :: kb = 8.6173324D-05
    real(kind=8), parameter :: h = 4.1356676D-15

    
    ! record scn of ipick and jpick
    ni_ads = cov_type(ipick)
    jpick = nnsite(j_cn, ipick)
    natoms_ij = 0

    call Merge_snn(natoms_ij, nnnsite(:, ipick), nnnsite(:, jpick))

    Ea = 0.0
    ! initial state (adsorption energy, U, bonding energy)
    U_ini = 0.0
    Eads_ini = 0.0
    do i = 1, 108
        nn_i = natoms_ij(i)
        if(nn_i > 0) then
            if(site_type(nn_i) < 3) then
                ! U value
                gcn_ni = gcn(nn_i)
                U_ini = U_ini + E_Mcohe * (corr_A1 * exp(-gcn_ni/corr_t1) + corr_A2 * exp(-gcn_ni/corr_t2) + 1)
                ! Eads
                if(cov_type(nn_i) == 1 .or. cov_type(nn_i) == 2) then
                    Eads_i = 0.0
                    call Eads_site(Eads_i, nn_i, 1)
                    Eads_ini = Eads_ini + Eads_i
                end if
            end if
        end if
    end do

    Eb_ini = E_Mbond * effcn(ipick)
    ijpick_cn(1, j_cn, ipick) = effcn(ipick)
    ijpick_gcn(1, j_cn, ipick) = gcn(ipick)

    ! middle state (adsorption energy, U, bonding energy)
    site_type(ipick) = 3
    cov_type(ipick) = 5
    U_mid = 0.0
    Eads_mid = 0.0
    do i = 1, 108
        nn_i = natoms_ij(i)
        if(nn_i > 0) then
            if(site_type(nn_i) < 3) then
                gcn_old = gcn(nn_i)
                ! U value
                ! cal gcn of nn_i
                gcn_ni = 0
                do j = 1, 12
                    neis_i = nnsite(j, nn_i)
                    if(site_type(neis_i) < 3) then
                        do k = 1, 12
                            neis_j = nnsite(k, neis_i)
                            if(site_type(neis_j) < 3) gcn_ni = gcn_ni + 1
                        end do
                    end if
                end do
                gcn_ni = gcn_ni / 12.0
                gcn(nn_i) = gcn_ni
                U_mid = U_mid + E_Mcohe * (corr_A1 * exp(-gcn_ni/corr_t1) + corr_A2 * exp(-gcn_ni/corr_t2) + 1)
                ! Eads
                if(cov_type(nn_i) == 1 .or. cov_type(nn_i) == 2) then
                    Eads_i = 0.0
                    call Eads_site(Eads_i, nn_i, 2)
                    Eads_mid = Eads_mid + Eads_i
                end if
                gcn(nn_i) = gcn_old   ! recovery
            end if
        end if
    end do
    ! 用SISSO拟合的公式算吸附能
    ! if(ni_ads == 1) Eads_mid = Eads_mid + Eads_ADS1_b
    ! if(ni_ads == 2) Eads_mid = Eads_mid + Eads_ADS2_b
    if(ni_ads == 1) Eads_mid = Eads_mid - 1.103
    if(ni_ads == 2) Eads_mid = Eads_mid - 0.1

    
    ! final state (adsorption energy, missing U, bonding energy)
    site_type(jpick) = 2
    cov_type(jpick) = ni_ads
    U_fin = 0.0
    Eads_fin = 0.0
    do i = 1, 108
        nn_i = natoms_ij(i)
        if(nn_i > 0) then
            if(site_type(nn_i) < 3) then
                gcn_old = gcn(nn_i)
                ! U value
                ! cal gcn of nn_i
                gcn_ni = 0
                do j = 1, 12
                    neis_i = nnsite(j, nn_i)
                    if(site_type(neis_i) < 3) then
                        do k = 1, 12
                            neis_j = nnsite(k, neis_i)
                            if(site_type(neis_j) < 3) gcn_ni = gcn_ni + 1
                        end do
                    end if
                end do
                gcn_ni = gcn_ni / 12.0
                gcn(nn_i) = gcn_ni
                U_fin = U_fin + E_Mcohe * (corr_A1 * exp(-gcn_ni/corr_t1) + corr_A2 * exp(-gcn_ni/corr_t2) + 1)
                ! Eads
                if(cov_type(nn_i) == 1 .or. cov_type(nn_i) == 2) then
                    cov_old = cov_type(nn_i)
                    if(cov_type(nn_i) == 1 .and. effcn(nn_i) >= maxCN_ADS) cov_type(nn_i) = 0
                    Eads_i = 0.0
                    call Eads_site(Eads_i, nn_i, 2)
                    Eads_fin = Eads_fin + Eads_i
                    cov_type(nn_i) = cov_old
                end if
                gcn(nn_i) = gcn_old    ! recovery
            end if
        end if
        if(nn_i == jpick) ijpick_gcn(2, j_cn, ipick) = gcn_ni    ! recovery
    end do
    effcn_jpick = 0
    do i = 1, 12
        neis_j = nnsite(i, jpick)
        if(site_type(neis_j) < 3) effcn_jpick = effcn_jpick + 1
    end do
    Eb_fin = E_Mbond * effcn_jpick
    ijpick_cn(2, j_cn, ipick) = effcn_jpick
    
    ! change of initial & final state
    Eb = (U_ini - U_mid) + (Eads_ini - Eads_mid) + (Eb_ini - 0) ! breaking energy of initial state
    Ef = (U_fin - U_mid) + (Eads_fin - Eads_mid) + (Eb_fin - 0) ! formation energy of final state
    if(Eb + Ef == 0) then
        Ea = 0
        write (*, *) 'Warning!!!!!!!!!!Eb+Ef.eq.0'
    else
        dE = Ef - Eb
        Ea = dE - (Ef**4)/((Eb+Ef)*(Eb*Eb+Ef*Ef))
        if(Eb > 0 .and. Ef < 0) Ea = 0
        if(Eb < 0 .and. Ef > 0) Ea = dE
        ! 对前后CN一致的原子跳事件，为Ea计算增加四次项修正
        ! if(Eb < 0 .and. Ef < 0) Ea = dE - Ef*Ef*Ef*Ef/((Eb+Ef)*(Eb*Eb+Ef*Ef))
        ! if(Eb < 0 .and. Ef < 0 .and. effcn(ipick) == effcn_jpick) Ea = dE - (Ef**8)/((Eb+Ef)*(Eb*Eb+Ef*Ef)*(Eb**4+Ef**4))
        ! if(Eb < 0 .and. Ef < 0 .and. abs(Eb-Ef) < 0.001) Ea = dE - (Ef**8)/((Eb+Ef)*(Eb*Eb+Ef*Ef)*(Eb**4+Ef**4))
    end if

    ! return initial state
    site_type(ipick) = 2
    cov_type(ipick) = ni_ads
    site_type(jpick) = 3
    cov_type(jpick) = 5
    
    ! ! 进行一次判断：如果Ea过大，则认为该事件无法发生，程序
    ! if()
    Ea_r(1, j_cn, ipick) = Ea
    r_site = (kb*Temp/h)*exp(-Ea/(kb*Temp))
    
end subroutine
