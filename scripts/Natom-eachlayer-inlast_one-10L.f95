program main
  implicit none

  real(kind=4), allocatable :: xxx(:),yyy(:),zzz(:), z_list(:)     ! 原子的xyz坐标
  real(kind=8) :: ctime        ! 读取的时间（原样读取写出）
  integer(kind=1), allocatable :: cov_type(:)
  integer(kind=4) :: nbulk, j, k, jump, n_step_tot, num_z_values, num_Mtotal
  integer(kind=4), allocatable :: z_count(:)
  character(kind=1, len=2) :: elem   ! 原子的元素种类
  character(kind=1, len=100) :: char_seq
  logical :: z_found
  

  ! 从last_one.xyz读取原子数
  open(10, file='last_one.xyz', status='old', action='read')
    read(10, *) nbulk
  close(10)
  allocate(xxx(nbulk))
  xxx = 0
  allocate(yyy(nbulk))
  yyy = 0
  allocate(zzz(nbulk))
  zzz = 0
  allocate(cov_type(nbulk))
  cov_type = 0
  allocate(z_list(nbulk))
  z_list = 0
  allocate(z_count(nbulk))
  z_count = 0
  
  ! 从last_one.xyz读取坐标
  open(11, file='n_atoms_highest_layer.dat', status='replace', action='write')
  open(10, file='last_one.xyz', status='old', action='read')

    read(10, *) nbulk
    read(10, *) char_seq, jump, n_step_tot, ctime
    write(*, *) jump
    num_z_values = 0
    do j = 1, nbulk
      read(10, *) elem, xxx(j), yyy(j), zzz(j), cov_type(j)
      z_found = .false.
      do k = 1, num_z_values
        if (abs(zzz(j) -z_list(k)) < 0.0001) then
          z_found = .true.
          z_count(k) = z_count(k) + 1
          exit
        end if
      end do
      if (.not. z_found) then  ! 增加这个zzz值到z_list里面
        num_z_values = num_z_values + 1
        z_list(num_z_values) = zzz(j)
        z_count(num_z_values) = 1
      end if
    end do

    ! 输出屏幕，用于测试
    write(*, *) "Z values and counts:"
    do k = 1, num_z_values
      print *, "Z =", z_list(k), "Count =", z_count(k)
    end do

    ! z_list按升序排列
    call sort_array(z_list, z_count, num_z_values)

    ! 输出到文件n_atoms_each_layer.dat
    write(11, '(I10, I10, E26.16, F10.2, I10, F10.2, I10, F10.2, I10, F10.2, I10, &
      F10.2, I10, F10.2, I10, F10.2, I10, F10.2, I10, F10.2, I10, F10.2, I10)') &
      n_step_tot, jump, ctime, z_list(num_z_values), z_count(num_z_values), &
      z_list(num_z_values-1), z_count(num_z_values-1), &
      z_list(num_z_values-2), z_count(num_z_values-2), &
      z_list(num_z_values-3), z_count(num_z_values-3), &
      z_list(num_z_values-4), z_count(num_z_values-4), &
      z_list(num_z_values-5), z_count(num_z_values-5), &
      z_list(num_z_values-6), z_count(num_z_values-6), &
      z_list(num_z_values-7), z_count(num_z_values-7), &
      z_list(num_z_values-8), z_count(num_z_values-8), &
      z_list(num_z_values-9), z_count(num_z_values-9)


  close(10)
  close(11)

  deallocate(xxx)
  deallocate(yyy)
  deallocate(zzz)
  deallocate(cov_type)
  deallocate(z_list)
  deallocate(z_count)

contains

  ! Subroutine to sort the array using insertion sort
  subroutine sort_array(arr, counts, size)
    real, intent(inout) :: arr(:)
    integer, intent(inout) :: counts(:)
    integer, intent(in) :: size
    real :: temp_z
    integer :: m, n, temp_count

    do m = 2, size
      temp_z = arr(m)
      temp_count = counts(m)
      n = m - 1
      do while (n > 0 .and. arr(n) > temp_z)
        arr(n + 1) = arr(n)
        counts(n + 1) = counts(n)
        n = n - 1
      end do
      arr(n + 1) = temp_z
      counts(n + 1) = temp_count
    end do
  end subroutine sort_array

end program
