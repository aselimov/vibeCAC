module atom_types
    !This module contains the atom type variables
    use parameters
    use logger
    use str
    implicit none


    !Max number of atom types allowed 
    integer, parameter :: max_atom_types=10
    !number of different potential types, eam all count for the same
    integer, parameter :: max_pot_types=2
    !Atom type variables
    integer, save ::  natom_types
    real(kind=wp), dimension(max_atom_types), save  :: masses
    character(len=2), dimension(max_atom_types), save :: atom_names
    !Have to move some potential variables here
    integer, save, dimension(max_atom_types, max_atom_types) :: types_to_pot_type
    logical, dimension(max_pot_types) :: potential_types
    logical :: atom_types_set

    public
    contains

    subroutine init_atom_types 
        !Initialization subroutine
        atom_types_set = .false.
        types_to_pot_type=0
    end subroutine init_atom_types

    subroutine set_atom_types(num,names)
        !This subroutine sets the atom types for use with the potential functions
        integer, intent(in) :: num
        character(len=2), dimension(:) :: names
        integer :: i

        natom_types = num
        do i = 1, natom_types
            atom_names(i) = names(i)
        end do 
        atom_types_set = .true.

        return
    end subroutine set_atom_types

    subroutine parse_types(line)
        character(len=*), intent(in) :: line
        character(len=read_len) :: tmptxt
        integer :: i

        i = tok_count(line)
        natom_types = i-1
        read(line,*) tmptxt, (atom_names(i), i = 1, natom_types)

        atom_types_set = .true.
        return
    end subroutine parse_types

    subroutine parse_mass(line)
        character(len=*), intent(in) :: line
        character(len=read_len) :: tmptxt
        integer :: i
        real(kind=wp) :: m

        read(line,*) tmptxt, i, m

        masses(i) = m

        write(tmptxt, *) "Masses are now: ", masses(1:natom_types)
        call log_msg(tmptxt)
    
        return
    end subroutine parse_mass

    subroutine log_types
        !This command logs the atom_types to the log file
        integer :: i
        character(len = read_len) :: msg
        msg = ''
        do i = 1, natom_types
            write(msg, *) trim(adjustl(msg)), i, atom_names(i)
        end do

        call log_msg("Atom types are mapped as "//msg)
    end subroutine log_types
end module atom_types

