module dump_logic

    implicit none

    public :: normalize_dump_every, should_write_dump, strip_dump_extension

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

    pure function strip_dump_extension(filename) result(stripped_filename)
        character(len=*), intent(in) :: filename
        character(len=len(filename)) :: stripped_filename
        integer :: dot_index

        stripped_filename = filename
        dot_index = scan(filename, ".")
        if (dot_index > 0) then
            stripped_filename(dot_index:) = ""
        end if
    end function strip_dump_extension

end module dump_logic
