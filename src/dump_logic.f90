module dump_logic

    implicit none

    public :: normalize_dump_every, should_write_dump

contains

    pure integer function normalize_dump_every(requested_every)
        integer, intent(in) :: requested_every

        normalize_dump_every = requested_every
        if (requested_every == 0) then
            normalize_dump_every = huge(1)
        end if
    end function normalize_dump_every

    pure logical function should_write_dump(timestep, dump_every, force_dump)
        integer, intent(in) :: timestep, dump_every
        logical, intent(in) :: force_dump

        should_write_dump = force_dump .or. (mod(timestep, dump_every) == 0)
    end function should_write_dump

end module dump_logic
