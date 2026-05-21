program main
    implicit none

    character(len=2) :: elem
    character (len=100) :: line
    real(kind=4) :: xxx(300000), yyy(300000), zzz(300000), gcn, ctime
    integer(kind=4) :: n, nframes, nbulk, natoms, i, j, grid_num(300000), atom_num, cov_type, cn
    integer(kind=4) :: ii_step, n_step_tot, io



    open(1, file="Total_bulk.xyz", action="read", status="old")  
        read(1, *) nbulk
        read(1, *)
        do i = 1, nbulk
            grid_num(i) = i
            read (1, *) elem, xxx(i), yyy(i), zzz(i)
        end do
    close(1)

    open(2, file="ini.xyz", action="read", status="old")
        read(2, *) natoms
    close(2)


    ! 统计总行数，从而计算共有多少帧
    open(3, file="atom_inf_00.dat", action="read", status="old")
        n = 0
        do
            read(3, *, iostat=io) line
            if(io /= 0) exit
            n = n + 1
        end do
    close(3)

    if(mod(n, (natoms + 1)) /= 0) write(*, *) "Warning: the last frame is not complete!"
    nframes = n / (natoms + 1)
    write(*, *) "this traj file has ", nframes, "frames"



    open(3, file="atom_inf_00.dat", action="read", status="old")
    open(5, file="traj.xyz", action="write", status="new")

    do i = 1, nframes
        read(3, *) ii_step, n_step_tot, ctime
        write(5, *) natoms
        write(5, *) ii_step, n_step_tot, ctime
        do j = 1, natoms
            read(3, *) atom_num, cov_type, cn, gcn
            if(cov_type == 0) then
                write(5, *) 'Pt', xxx(atom_num), yyy(atom_num), zzz(atom_num), cn, gcn
            elseif(cov_type == 1) then
                write(5, *) 'Au', xxx(atom_num), yyy(atom_num), zzz(atom_num), cn, gcn
            elseif(cov_type == 2) then
                write(5, *) 'Ir', xxx(atom_num), yyy(atom_num), zzz(atom_num), cn, gcn
            elseif(cov_type == 5) then
                write(5, *) 'Co', xxx(atom_num), yyy(atom_num), zzz(atom_num), cn, gcn
            else
                write(*, *) "Error: cov_type wrong! frame: ", i, "line:", natoms
            end if
        end do
    end do

    close(3)
    close(5)

end program main