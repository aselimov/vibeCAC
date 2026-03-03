module read_dump
    use parameters
    use mpi
    implicit none
    real(kind=dp), allocatable :: read_atom_pos(:,:), read_atom_force(:,:), read_atom_energy(:)
    integer, allocatable :: read_atom_id(:)
    public 
    contains

    subroutine read_in_dump(readfile)
        character(len=*), intent(in) :: readfile
        character(len=1000) :: tmptxt
        integer :: i, j, natoms, stat

        open(unit=1, file=readfile, action='read', status='old', position='rewind')
        
        !Read header info
        read(1, *) tmptxt
        read(1, *) tmptxt
        read(1, *) tmptxt

        !Read atom number and setup the arrays
        read(1,*) natoms
        if (allocated(read_atom_id)) then 
            deallocate(read_atom_id,read_atom_pos,read_atom_force,read_atom_energy)
        end if
        allocate(read_atom_id(natoms),read_atom_pos(3,natoms),read_atom_force(3,natoms),read_atom_energy(natoms))


        !Read the rest of the useless info 
        read(1, *) tmptxt
        read(1, *) tmptxt
        read(1, *) tmptxt
        read(1, *) tmptxt
        read(1, *) tmptxt
        
        !Read in the atoms now
        do i= 1,natoms
            read(1, '(A)') tmptxt
            read(tmptxt, *, iostat=stat) read_atom_id(i), j, read_atom_pos(:, i), read_atom_force(:, i), read_atom_energy(i)
            if(stat > 0) then 
                print *, tmptxt
                call mpi_abort(mpi_comm_world, 1, ierr)
            end if
        end do
        close(1)
    end subroutine read_in_dump
end module read_dump
