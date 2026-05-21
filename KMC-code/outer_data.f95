module outer_data
    implicit none
    ! This module is used to read and save parameters form 'input'. All parameters is read-only.

    integer(kind=4), parameter :: N_EVENT_TYPE = 8    ! number of events
    ! 1: atom jumping;   2,3: ADS1 ads & des;   4,5: ADS2 ads & des;   6,7,: ADS1,ADS2 diff;   8: reaction(没写); 

    integer(kind=1), protected :: flag_restart, flag_reseed  
    integer(kind=8), protected :: nLoop, record_int, struct_int

    integer(kind=4), protected :: maxCN_ADS
    real(kind=8), protected :: s0_ADS1_edge, s0_ADS1_facet, s0_ADS2_edge, s0_ADS2_facet
    real(kind=4), protected :: Ediff_ADS1, Ediff_ADS2
    ! real(kind=8), protected :: Eads_ADS1_a, Eads_ADS1_b, Eads_ADS2_a, Eads_ADS2_b
    real(kind=8), protected :: E_Mbond, E_Mcohe, corr_A1, corr_t1, corr_A2, corr_t2
    real(kind=8), protected :: BEP1_a, BEP1_b
    real(kind=4), protected :: Temp, p_tot, ppADS1, ppADS2, U_SHE

    character(len=2) :: elem, specie
    logical :: enable_evt(N_EVENT_TYPE)
    real(kind=4), protected :: rx_slab, ry_slab, rz_slab
    real(kind=4), protected :: EADS1_ADS1, EADS1_ADS2, EADS2_ADS2


    contains 
    
    subroutine read_data
        open(1, file="input", action="read")
            read(1, *) ! operational control
            read(1, *) flag_restart, flag_reseed     ! 1st: 0:generate new struc, 1:read struc from existing last_one;     2nd: 0:generate new seed and record, 1:use old seed
            read(1, *) nLoop, record_int, struct_int ! steps
            read(1, *) specie
            read(1, *) enable_evt
            read(1, *) maxCN_ADS

            read(1, *) ! struc & energy parameter for metal
            read(1, *) elem
            if(specie == 'S') then 
                read(1, *) rx_slab, ry_slab, rz_slab
            elseif(specie == 'P') then
                read(1, *)   !!! 还没写！
            else
                write(*, *) "reading specie error! specie =", specie
                stop
            end if
            read(1, *) E_Mbond, E_Mcohe, corr_A1, corr_t1, corr_A2, corr_t2

            read(1, *) ! energy parameter for ADS    
            read(1, *) s0_ADS1_edge, s0_ADS1_facet ! sticking coffecient of ADS1
            read(1, *) s0_ADS2_edge, s0_ADS2_facet ! sticking coffecient of ADS2
            read(1, *) Ediff_ADS1, Ediff_ADS2      ! diffusion barrier
            ! read(1, *) Eads_ADS1_a, Eads_ADS1_b    ! Eads = a*gcn + b
            ! read(1, *) Eads_ADS2_a, Eads_ADS2_b
            read(1, *) EADS1_ADS1, EADS1_ADS2, EADS2_ADS2       ! Interaction between ADS1-ADS1; ADS1-ADS2; ADS2-ADS2  
            read(1, *) BEP1_a, BEP1_b              ! BEP coeffient for reaction

            read(1, *) ! environment 
            read(1, *) Temp, p_tot, U_SHE
            read(1, *) ppADS1, ppADS2 ! partial pressure of ADS1 & ADS2
        close(1)
    end subroutine read_data  


    subroutine init_random_seed(seed_flag)   ! random seed
        implicit none
      
        integer(kind=1), intent(IN) :: seed_flag
        integer(kind=4) :: i, n, clock
        integer, dimension(:), allocatable :: seed  ! ArrayList
      
        ! seed_flag = 0：generate new seed and record
        if(seed_flag == 0) then
            call random_seed(size = n)
            allocate(seed(n))
            call system_clock(count = clock)
            seed = clock + 37 * (/(i - 1, i=1, n)/)
            call random_seed(put = seed)
            open(1, file="seed.txt", status="replace", action="write")
                write(1, *) n
                write(1, *) clock
            close(1)

        ! seed_flag = 1：use old seed
        elseif (seed_flag == 1) then
            open(1, file="seed.txt", status="old", action="read")
                read (1, *) n
                read (1, *) clock
            close(1)
            allocate(seed(n))
            seed = clock + 37 * (/(i - 1, i=1, n)/)
            call random_seed(put = seed)
        end if
        deallocate(seed)
    end subroutine init_random_seed

end module outer_data