module input_validation
    implicit none

    public :: first_missing_required_command

contains

    pure integer function first_missing_required_command(req_flags)
        logical, intent(in) :: req_flags(:)
        integer :: i

        first_missing_required_command = 0
        do i = 1, size(req_flags)
            if (.not. req_flags(i)) then
                first_missing_required_command = i
                return
            end if
        end do
    end function first_missing_required_command

end module input_validation
