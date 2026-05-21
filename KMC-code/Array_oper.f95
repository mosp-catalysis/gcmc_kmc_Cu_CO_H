module Array_oper
    implicit none
    public :: Merge_snn

contains
    subroutine Merge_snn(merged_one, snn1, snn2)
        implicit none
        ! 合并+去重
        ! 利用了snn的单调递增特性
        integer(kind=4), intent(inout), dimension(:) :: merged_one
        integer(kind=4), intent(in), dimension(:) :: snn1, snn2
        
        integer(kind=4) :: s1, s2, q1
        q1 = 1
        s1 = 1
        s2 = 1
        do while(snn1(s1).ne.0 .or. snn2(s2).ne.0)
            if (snn1(s1) .eq. 0) then
                merged_one(q1) = snn2(s2)
                s2 = s2 + 1
            else if (snn1(s2) .eq. 0) then
                merged_one(q1) = snn1(s1)
                s1 = s1 + 1
            else if (snn1(s1) .lt. snn2(s2)) then
                merged_one(q1) = snn1(s1)
                s1 = s1 + 1
            else if (snn1(s1) .gt. snn2(s2)) then
                merged_one(q1) = snn2(s2)
                s2 = s2 + 1
            else if (snn1(s1) .eq. snn2(s2)) then
                merged_one(q1) = snn1(s1)
                s1 = s1 + 1
                s2 = s2 + 1
            end if
            q1 = q1 + 1
        end do
    end subroutine
end module Array_oper