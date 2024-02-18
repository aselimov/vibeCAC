module comms
    !This is code in charge of most interprocessor communications

    use mpi
    use parameters
    use elements
    use box
    use errors
    use math

    implicit none

    !Communicator parameters
    integer, save :: grank, pro_num, grid_comm, mpi_wp, & 
                     nreplica, ireplica, world, universe, uroots_comm, &
                     urank, upro_num

    character(len=100) :: rank_string

    !Processor arrays    
    integer ::  &
        num_pro(3), &
        grid_coords(3), &
        list_atom_num(2, 3), &
        list_atomap_num(2, 3)

    integer, allocatable, save :: &
        list_atom(:, :, :), list_atomap(:, :, :)
 
    real(kind = wp), save :: &
        pro_bd(6), &
        pro_bd_pre(6), &
        pro_bd_out(6), &
        pro_bd_in(6), &
        pro_bd_kj(6, 2, 3), &
        pro_length(3), &
        pro_length_out(3)

    logical, save :: &
        pro_gather(3), pro_newton(2), &
        pro_atomap(2), &
        list_atom_logic(2, 3), &
        list_atomap_logic(2, 3)


    !Scatter arrays for atoms
    integer, allocatable, save :: tag_atom_scatter(:), type_atom_scatter(:)
    real(kind=wp), allocatable, save :: r_atom_scatter(:,:), vel_atom_scatter(:,:)
    !Scatter arrays for elements
    integer, allocatable, save :: tag_ele_scatter(:), size_ele_scatter(:), etype_scatter(:), &
                                  basis_num_scatter(:), basis_type_scatter(:,:)
    real(kind=wp), allocatable, save :: r_scatter(:,:,:), vel_scatter(:,:,:)

    !arrays for sharing elements
    integer :: ele_shared_num
    integer, allocatable :: pro_shared_num(:), ele_id_shared(:)
    logical, allocatable :: who_has_ele(:)

    public 
    contains 
    subroutine processor_array
 
        implicit none
 
        allocate( &
            list_atom(seg_num, 2, 3), &
            stat = allostat &
        )
 
        return 
    end subroutine processor_array

    subroutine scatter_cg_array
        !Allocate scatter arrays for cg
        if(allocated(tag_ele_scatter)) then 
            deallocate(tag_ele_scatter, r_scatter, vel_scatter, size_ele_scatter, etype_scatter, &
                       basis_num_scatter, basis_type_scatter)
        end if
        allocate( &
            r_scatter(3, max_basisnum, node_num), &
            vel_scatter(3, max_basisnum, node_num), &
            tag_ele_scatter(ele_num), &
            size_ele_scatter(ele_num), &
            etype_scatter(ele_num), &
            basis_num_scatter(ele_num), &
            basis_type_scatter(max_basisnum, ele_num), &
            stat = allostat &
        )
   
    end subroutine scatter_cg_array

    subroutine dealloc_scatter_cg
        !Deallocate scatter arrays for cg
        deallocate( r_scatter, vel_scatter, tag_ele_scatter, size_ele_scatter, etype_scatter, basis_num_scatter, &
        basis_type_scatter, stat = allostat )
    end subroutine dealloc_scatter_cg

    subroutine scatter_at_array
        !Allocate scatter arrays for at
        if(allocated(tag_atom_scatter)) then 
            deallocate(tag_atom_scatter, type_atom_scatter, r_atom_scatter, vel_atom_scatter)
        end if
        allocate( &
            tag_atom_scatter(atom_num), type_atom_scatter(atom_num),  &
            r_atom_scatter(3,atom_num), vel_atom_scatter(3,atom_num), &
            stat = allostat &
        )
    end subroutine scatter_at_array

    subroutine dealloc_scatter_at
        !Deallocate scatter arrays for at
        deallocate( &
            tag_atom_scatter, type_atom_scatter,  &
            r_atom_scatter, vel_atom_scatter, &
            stat = allostat &
        )
    end subroutine dealloc_scatter_at

    subroutine precision_init
        !This function defines the mpi precison to be used

        select case(wp)
            case(-1, -2, -3)
                if(rank == root) then
                    print *, 'Error: The processor does not support the desired', &
                             ' real data type precision'
                    call mpi_abort(universe, 1, ierr)
                 end if

            case(selected_real_kind(15, 307))
                mpi_wp = mpi_double_precision

            case(selected_real_kind(33, 4931))
                mpi_wp = mpi_real16
 
            case default
               if(rank == root) then
                    print *, 'Error: Wp', wp, ' is not accepted'
                    call mpi_abort(universe, 1, ierr)
               end if
        end select

        return
    end subroutine precision_init

    subroutine comm_init
        !This is a simple subroutine which initializes the communicator information

        call mpi_init(ierr)
        universe=mpi_comm_world
        world=universe
        call mpi_comm_rank(mpi_comm_world, urank, ierr)
        write(rank_string, *) "urank ", urank !rank_string is just to make it easier to pass the rank to error subroutines
        call mpi_comm_size(mpi_comm_world, upro_num, ierr)
        call precision_init
        
        universe = mpi_comm_world
        world = universe
        uroots_comm = universe
        pro_num = upro_num
        rank = urank
        
    end subroutine
    
    subroutine multicomm_init
        ! initializes a multicomm communicator setup with nreplica sub-worlds

        ! check to make sure we can actually split the communicator
        if ( mod(upro_num, nreplica) /= 0 ) then
            print *, 'Error: Cannot evenly divide ', upro_num, ' processors into', &
                      nreplica, ' replica images.'
            call mpi_abort(universe, 1, ierr)
        endif 
        
        ! split comm world, make a universe communicator and save parameters
        ireplica = mod(urank, nreplica)
        call mpi_comm_split(universe, ireplica, urank, world, ierr)
        call mpi_comm_rank(world, rank, ierr)
        call mpi_comm_size(world, pro_num, ierr)

        ! create a communicator with all the world roots for easy bcasting
        call mpi_comm_split(universe, rank, urank, uroots_comm, ierr)

    end subroutine

    subroutine divide_domain
        !This subroutine divides the simulation domain among processors. It does rectilinear domain decomposition, 
        !but attempts to max the calculation points 


        integer :: i, k, ix, iy, iz, pro_num_try, i_this
        logical :: ireo
        integer, dimension(3) :: num_index
        integer, allocatable :: &
          num_pro_try(:, :), &
          num_pro_try_array(:, :)
        real(kind = wp), allocatable :: &
          area_try(:)

        if(pro_num == 1) then
            !If only one processor then the whole domain is assigned to the root processor
            num_pro(:) = 1
            grid_coords(:) = 0
            grank = 0
            do i = 1, 3
                pro_length(i) = box_length(i)
                pro_bd(2*i-1) = box_bd(2*i-1)
                pro_bd(2*i) = box_bd(2*i)
            end do
        else
            !Otherwise we have to split the domain

            pro_num_try = pro_num
            allocate(num_pro_try(pro_num_try, 3), stat = allostat) 
            if(allostat /= 0) call alloc_error(rank_string//" failed to allocate num_pro_try in divide domain",allostat)

            num_pro_try(:, :) = 0
            k = 0
            do ix = 1, pro_num
                if(mod(pro_num, ix) == 0) then
                    do iy = 1, pro_num / ix
                        if(mod(pro_num, ix * iy) == 0) then
                            iz = pro_num / (ix * iy) 
                            k = k + 1

                            if(k > pro_num_try) then
                                allocate( &
                                    num_pro_try_array(pro_num_try+20, 3), &
                                    stat = allostat &
                                )

                                if(allostat /= 0) call alloc_error(rank_string// &
                                                                   " failed to allocate num_pro_try_array in divide_domain", &
                                                                    allostat)

                                num_pro_try_array(:, :) = 0
                                num_pro_try_array(1:pro_num_try, :) = num_pro_try(:, :)
                                call move_alloc(num_pro_try_array, num_pro_try)
                                pro_num_try = pro_num_try + 20
                            end if

                            num_pro_try(k, :) = [ ix, iy, iz ]
                        end if
                    end do
                end if
            end do

            !compare all possible set of num_pro_try
            allocate(area_try(k),stat = allostat) 
            if(allostat /= 0) call alloc_error(rank_string // " failed to allocate area_try in divide_domain ", allostat)

            area_try(:) = 0.0_wp

            do i = 1, k

                num_index(:) = num_pro_try(i, :)

                if(period(1).eqv..true.) then
                    area_try(i) = area_try(i) + (num_index(1) + 1) &
                                * box_length(2) * box_length(3)
                else
                    area_try(i) = area_try(i) + (num_index(1) - 1) &
                                * box_length(2) * box_length(3)

                end if

                if(period(2).eqv..true.) then
                    area_try(i) = area_try(i) + (num_index(2) + 1) &
                                * box_length(3) * box_length(1)

                else
                    area_try(i) = area_try(i) + (num_index(2) - 1) &
                                * box_length(3) * box_length(1)

                end if

                if(period(3).eqv..true.) then
                    area_try(i) = area_try(i) + (num_index(3) + 1) &
                                * box_length(1) * box_length(2)

                else
                    area_try(i) = area_try(i) + (num_index(3) - 1) &
                                * box_length(1) * box_length(2)

                end if

            end do

            i_this = minloc(area_try, 1)

            num_pro(:) = num_pro_try(i_this, :)

            !Check to make sure everything has been setup correctly
            if(rank == root) then
                !Check to make sure the number of domains equals the total number of processors
                if(pro_num /= product(num_pro)) then
                    print *, 'Error: Wrong num_pro array, the product', &
                             product(num_pro), ' does not equal', pro_num
                    call mpi_abort(universe, 1, ierr)
                end if
            end if

            !create a cartesian grid
            ireo = .false.
            call mpi_cart_create(world, 3, num_pro, period, ireo, grid_comm, ierr)
            call mpi_comm_rank(grid_comm, grank, ierr)
            call mpi_cart_coords(grid_comm, grank, 3, grid_coords, ierr) 
        end if

         
        pro_length(:) = 0.0_wp
        pro_bd(:) = 0.0_wp

        do i = 1, 3
            !Define processor boundaries
            pro_length(i) = box_length(i) / num_pro(i)
            pro_bd(2*i-1) = box_bd(2*i-1) + grid_coords(i) * pro_length(i)
            pro_bd(2*i) = box_bd(2*i-1) + (grid_coords(i) + 1) * pro_length(i)

            if(grid_coords(i) == 0) then
                !If the processor has coordinate 0 it's the bottom one so set the boundary equal to the box boundary
                pro_bd(2*i-1) = box_bd(2*i-1)
            end if

            if(grid_coords(i) == num_pro(i) - 1) then
                !If the processor has coordinate num_pro(i) - 1 it's the top one so set the boundary equal to top boundary
                pro_bd(2*i) = box_bd(2*i)
            end if

            pro_length(i) = pro_bd(2*i) - pro_bd(2*i-1)
        end do

        return
    end subroutine divide_domain
    
    subroutine processor_bds
        !This subroutine defines processor boundaries
        integer :: i, j, k

        do i = 1, 3
            !Calculate the pro_length in each dimension and make sure it's larger than 2*rc_neigh
            pro_length(i) = box_length(i) / num_pro(i)
            if(pro_length(i) < 2.0_wp * rc_neigh) then
                print *, 'Error: pro_length along the', i, ' direction', &
                          pro_length(i), ' of rank', rank, ' is smaller than 2*rc_neigh:', &
                          2.0_wp * rc_neigh

                call mpi_abort(universe, 1, ierr)
            end if
        end do

        !define inner boundaries
        !Inner boundaries are used to check if atoms are within cutoff radius of processor
        !boundary requiring creation of ghost atoms.

        pro_bd_in(:) = pro_bd(:)
        do i = 1, 3
            pro_bd_in(2*i-1) = pro_bd(2*i-1) + rc_neigh
            pro_bd_in(2*i) = pro_bd(2*i) - rc_neigh
        end do

        !define outer boundaries
        pro_bd_out(:) = pro_bd(:)
        do i = 1, 3
            pro_bd_out(2*i-1) = pro_bd(2*i-1) - rc_neigh
            pro_bd_out(2*i) = pro_bd(2*i) + rc_neigh
            if(period(i).eqv..false.) then

                if(grid_coords(i) == 0) then
                    pro_bd_out(2*i-1) = box_bd(2*i-1)
                end if

                if(grid_coords(i) == num_pro(i) - 1) then
                    pro_bd_out(2*i) = box_bd(2*i)
                end if
            end if

            pro_length_out(i) = pro_bd_out(2*i) - pro_bd_out(2*i-1)

        end do

        !Now get processor kj boundaries. These are the inner and outer boundaries for every step in the process of sharing ghost
        !atoms/elemenmts
        do j = 1 ,3
            do k = 1, 2
                pro_bd_kj(:, k, j) = pro_bd(:)

                if(k == 1) then
                    pro_bd_kj(2*j-1, k, j) = pro_bd_in(2*j)

                else
                    pro_bd_kj(2*j, k, j) = pro_bd_in(2*j-1)
                end if
                if(j == 2) then
                    pro_bd_kj(1:2, k, j) = pro_bd_out(1:2)
            
                else if(j == 3) then
                    pro_bd_kj(1:4, k, j) = pro_bd_out(1:4)
                end if
            end do
        end do
        return
    end subroutine processor_bds

    subroutine update_proc_bd(opt1)
        !This subroutine updates processor boundaries if needed 
        logical, intent(in), optional :: opt1

        integer :: i, ip, ibasis
        real(kind = wp) :: max_posme(3), max_posall(3), min_posme(3), min_posall(3)
        logical :: resize_box, resize_all

        if(present(opt1)) then 
            resize_all=opt1
        else
            resize_all = .true.
        end if

        max_posme(:) = 0
        min_posme(:) = 0
        resize_box = .false.

        if(resize_all) then
            do i = 1, 3
                !Only need to resize if shrink wrapped in that boundary
                if(.not.period(i)) then 
                    !Get atomic max and min
                    if(atom_num_l > 0) then 
                        max_posme(i) = maxval(r_atom(i, 1:atom_num_l))
                        min_posme(i) = minval(r_atom(i, 1:atom_num_l))
                    end if

                    !Now get nodal max and min
                    do ip = 1, node_num_l
                        do ibasis = 1, basis_num(node_cg(ip))
                            max_posme(i) = max(max_posme(i), r(i, ibasis, ip)) 
                            min_posme(i) = min(min_posme(i), r(i, ibasis, ip))
                        end do
                    end do
                end if
            end do

            call mpi_allreduce(max_posme, max_posall, 3, mpi_wp, mpi_max, world, ierr)
            call mpi_allreduce(min_posme, min_posall, 3, mpi_wp, mpi_min, world, ierr)

            do i = 1, 3
                if (.not.period(i)) then 
                    if (max_posall(i) > box_bd(2*i)) then 
                        box_bd(2*i) = max_posall(i)+10D-4
                        resize_box = .true.
                    end if
                    if (min_posall(i) < box_bd(2*i-1)) then 
                        box_bd(2*i-1) = min_posall(i)-10D-4
                        resize_box = .true.
                    end if
                end if
                box_length(i) = box_bd(2*i) - box_bd(2*i-1)
            end do


        else
            !If not then we only check to resize the outer boundary
            do i = 1, 3
                !Only need to resize if shrink wrapped in that boundary
                if(.not.period(i)) then 
                    !Get atomic max and min
                    if(atom_num_l > 0) then 
                        max_posme(i) = maxval(r_atom(i, 1:atom_num_l))
                        min_posme(i) = minval(r_atom(i, 1:atom_num_l))
                    end if

                    !Now get nodal max and min
                    do ip = 1, node_num_l
                        do ibasis = 1, basis_num(node_cg(ip))
                            max_posme(i) = max(max_posme(i), r(i, ibasis, ip)) 
                            min_posme(i) = min(min_posme(i), r(i, ibasis, ip))
                        end do
                    end do
                end if
            end do
            
            call mpi_allreduce(max_posme, max_posall, 3, mpi_wp, mpi_max, world, ierr)
            call mpi_allreduce(min_posme, min_posall, 3, mpi_wp, mpi_min, world, ierr)

            !Now check to see if the max and min are greater than the box boundaries for shrink wrapped directions
            do i = 1, 3
                if(.not.period(i)) then 
                    if(min_posall(i) < box_bd(2*i-1)) then 
                        box_bd(2*i-1) = min_posall(i) - zerotol
                        !If the processor has coordinate 0 it's the bottom one so set the boundary equal to the box boundary
                        if(grid_coords(i) == 0) then
                            pro_bd(2*i-1) = box_bd(2*i-1)
                        end if
                    end if

                    if(max_posall(i) > box_bd(2*i)) then 
                        box_bd(2*i) = max_posall(i) + zerotol
                        if(grid_coords(i) == num_pro(i) - 1) then
                        !If the processor has coordinate num_pro(i) - 1 it's the top one so set the boundary equal to top boundary
                            pro_bd(2*i) = box_bd(2*i)
                        end if
                    end if
                end if
                box_length(i) = box_bd(2*i) - box_bd(2*i-1)
            end do
        end if

        !If we change the box size then resize processor boundaries
        if(resize_box) then 
            pro_length(:) = 0.0_wp
            pro_bd(:) = 0.0_wp

            do i = 1, 3
                !Define processor boundaries
                pro_length(i) = box_length(i) / num_pro(i)
                pro_bd(2*i-1) = box_bd(2*i-1) + grid_coords(i) * pro_length(i)
                pro_bd(2*i) = box_bd(2*i-1) + (grid_coords(i) + 1) * pro_length(i)

                if(grid_coords(i) == 0) then
                    !If the processor has coordinate 0 it's the bottom one so set the boundary equal to the box boundary
                    pro_bd(2*i-1) = box_bd(2*i-1)
                end if

                if(grid_coords(i) == num_pro(i) - 1) then
                    !If the processor has coordinate num_pro(i) - 1 it's the top one so set the boundary equal to top boundary
                    pro_bd(2*i) = box_bd(2*i)
                end if

                pro_length(i) = pro_bd(2*i) - pro_bd(2*i-1)
            end do
            call processor_bds
        end if
        return
    end subroutine update_proc_bd

    subroutine scatter_cg
        !Scatter coarse-grained elements to the different processors
        integer :: i, j, k, ie, je, le, jp, ip, inod, lp, iatomap, jatomap, latomap, &
                   seg_real, iatom, irank, jrank,  &
                   ibasis, etype_buff, basis_num_buff, size_ele_buff, &
                   ele_atomap_num, tag_ele_buff
              
        integer, dimension(3) :: info
        integer, dimension(max_basisnum) :: basis_type_buff
        integer, dimension(3,max_basisnum,ng_max_node) :: pb_node_buff
        integer, dimension(atomap_max) :: tag_atomap_buff, type_atomap_buff, &
                                              iatom_array, jatom_array

        !Add integer variables for setting up interpolation
        integer :: max_esize, unique_esizes(100), unique_enum, send_array(2)

        integer, dimension(pro_num) :: irank_array
        integer, dimension(3, max_basisnum, ng_max_node) :: pb_in

        real(kind = wp), dimension(3, max_basisnum) :: r_interp
        real(kind=wp), dimension(3, max_basisnum, ng_max_node) :: r_nodes, vel_nodes
        real(kind = wp), dimension(3,atomap_max) :: r_atomap_buff

        logical(kind = wp), dimension(atomap_max) :: need_ele
        logical(kind = wp), dimension(pro_num) :: logic_array

        !Variables for communication
        integer :: send_ni, send_nr 
        integer, dimension(4+max_basisnum+(1+3*max_basisnum)*ng_max_node+atomap_max) :: send_int
        real(kind=wp), dimension(2*3*max_basisnum*ng_max_node + 3*atomap_max) :: send_real

        integer, allocatable :: &
            who_ele(:), who_ele_all(:), who_element(:), who_element_all(:), who_ele_long(:), who_ele_long_all(:), &
            ele_id_buff(:)



        !First setup the atom interpolation
        max_esize = 0
        unique_enum = 0
        if (rank == root) then 
            eleloop:do i = 1, ele_num
                !Figure out max esize
                if(size_ele_scatter(i) > max_esize) max_esize = size_ele_scatter(i) 
                do j = 1, unique_enum
                    if (unique_esizes(j) == size_ele_scatter(i)) cycle eleloop
                end do
                unique_enum = unique_enum + 1
                unique_esizes(unique_enum) = size_ele_scatter(i)
            end do eleloop
            send_array(1) = unique_enum
            send_array(2) = max_esize
            call mpi_bcast(send_array, 2, mpi_integer, root, world, ierr)
            call mpi_bcast(unique_esizes(1:unique_enum), unique_enum, mpi_integer, root, &
                           world, ierr)
        else
            call mpi_bcast(send_array, 2, mpi_integer, root, world, ierr)
            unique_enum = send_array(1)
            max_esize = send_array(2)
            call mpi_bcast(unique_esizes(1:unique_enum), unique_enum, mpi_integer, root, &
                           world, ierr) 
        end if

        !Now setup the interpolation arrays
        call setup_interpolation(max_esize, unique_enum, unique_esizes(1:unique_enum))

        !Allocate some variables
        allocate(who_ele(ele_num), who_ele_long(pro_num*ele_num), stat=allostat) 
        if (allostat > 0) call alloc_error("Failure allocating who_ele in scatter_cg", allostat) 
        who_ele(:) = 0
        who_ele_long(:) = 0

        !If we only have one processor the code is simple
        if (pro_num == 1) then 

            le=0
            lp=0
            latomap=0

            do ie = 1, ele_num
            
                le = le + 1
                ele_glob_id(ie) = ie
                tag_ele(ie) = tag_ele_scatter(ie)
                size_ele(ie) = size_ele_scatter(ie)
                etype(ie) = etype_scatter(ie)
                basis_num(ie) = basis_num_scatter(ie)
                basis_type(:,ie) = basis_type_scatter(:,ie)
                who_ele(ie) = who_ele(ie) + 1
                who_ele_long(pro_num*(ie-1)+rank+1) = 1

                do inod = 1, ng_node(etype(ie))
 
                    lp = lp + 1
                    cg_node(inod, ie) = lp
                    node_cg(lp) = ie
                    
                    do ibasis = 1, basis_num(ie)

                        if(periodic) then
                            call cross_pb(r_scatter(:,ibasis,lp), info)
                            pb_node(:, ibasis, lp) = info(:)
                            pb_in(:, ibasis, inod) = info(:)
                        end if
                        r(:, ibasis, lp) = r_scatter(:,ibasis,lp)
                        r_nodes(:,ibasis,inod) = r(:,ibasis, lp)
                        if(need_vel) then
                            vel(:, ibasis, lp) = vel_scatter(:, ibasis, lp)
                        end if

                    end do

                end do

                do iatom = 1, (size_ele(ie)+1)**3
                    !Interpolate atoms
                    call interp_atom(iatom, size_ele(ie), etype(ie), pb_in, basis_num(ie), r_nodes, r_interp)

                    do ibasis = 1, basis_num(ie)
                        !Assign interpolated atoms to arrays
                        latomap = latomap + 1
                        tag_atomap(latomap) = latomap
                        type_atomap(latomap) = basis_type(ibasis,ie)
                        cg_atomap(basis_num(ie)*(iatom-1) + ibasis, ie) = latomap
                        r_atomap(:, latomap) = r_interp(:,ibasis)
                    end do
                end do

            end do

            if(le /= ele_num) then
                print *, 'Error: Le', le, ' should equal', ele_num
                call mpi_abort(universe, 1, ierr)

            else if(lp /= node_num) then
                print *, 'Error: Lp', lp, ' should equal', node_num
                call mpi_abort(universe, 1, ierr)

            else if(latomap /= atomap_num) then
                print *, 'Error: Latomap', latomap, ' should equal', atomap_num
                call mpi_abort(universe, 1, ierr)
            end if

        else
            !If we have more than one processor than we have to run more complex code
            ip = 0
            lp = 0
            le = 0
            iatomap = 0
            jatomap = 0
            latomap = 0
            
            !Initialize send counts for integer data and real data
            if (periodic) then 
                send_ni = 4+max_basisnum+(1+3*max_basisnum)*ng_max_node+atomap_max
            else 
                send_ni = 4+max_basisnum+ng_max_node+atomap_max
            end if

            if(need_vel) then 
                send_nr = 2*3*max_basisnum*ng_max_node + 3*atomap_max
            else
                send_nr = 3*max_basisnum*ng_max_node+3*atomap_max
            end if

            do ie = 1, ele_num

                tag_atomap_buff(:)  = 0
                type_atomap_buff(:) = 0
                basis_type_buff(:)  = 0
                pb_node_buff(:,:,:) = 0
                r_nodes(:,:,:)      = 0.0_wp
                vel_nodes(:,:,:)    = 0.0_wp
                r_atomap_buff(:,:)  = 0.0_wp

                if(rank == root) then

                    pb_in(:,:,:)            =   0
                    do inod = 1, ng_node(etype_scatter(ie))
                        ip = ip + 1

                        do ibasis = 1, basis_num_scatter(ie)
                            !Now get all of the correct nodal positions
                            if(periodic.eqv..true.) then
                                call cross_pb(r_scatter(:,ibasis,ip), info)
                                pb_in(:, ibasis, inod) = info
                            end if

                            !r_nodes will be used later interpolation
                            r_nodes(:,ibasis,inod) = r_scatter(:,ibasis,ip)
                            !vel nodes will be sent 
                            vel_nodes(:,ibasis,inod) = vel_scatter(:, ibasis, ip)

                        end do
                    end do

                    do iatom = 1, (size_ele_scatter(ie)+1)**3
                        !Now interpolate all of the atoms
                        !NOTE in this case iatom is the number of lattice points
                        call interp_atom(iatom, size_ele_scatter(ie), etype_scatter(ie), pb_in, basis_num_scatter(ie), &
                                         r_nodes, r_interp)

                            do ibasis = 1, basis_num_scatter(ie)
                                iatomap = iatomap + 1
                                tag_atomap_buff((iatom-1)*basis_num_scatter(ie) + ibasis) = iatomap
                                do i = 1,3
                                    r_atomap_buff(i,(iatom-1)*basis_num_scatter(ie) + ibasis) = r_interp(i,ibasis)
                                end do
                        end do
                    end do

                    if (need_vel) then 
                        !Now pack the element data 
                        call pack_ele_atomap(etype_scatter(ie),tag_ele_scatter(ie),size_ele_scatter(ie), basis_num_scatter(ie),&
                                             basis_type_scatter(:,ie), pb_in, r_nodes,  tag_atomap_buff, &
                                             r_atomap_buff, send_int, send_real, vel_nodes)                
                    else 
                        call pack_ele_atomap(etype_scatter(ie),tag_ele_scatter(ie),size_ele_scatter(ie), basis_num_scatter(ie),&
                                             basis_type_scatter(:,ie), pb_in, r_nodes, tag_atomap_buff, &
                                             r_atomap_buff, send_int, send_real)                
                    end if
                end if


                !Send integer data
                call mpi_bcast(send_int(1:send_ni), send_ni, mpi_integer, root, world, ierr)
                !Send real data
                call mpi_bcast(send_real(1:send_nr), send_nr, mpi_wp, root, world, ierr)


                !Now unpack the data into the buff arrays
                if (need_vel) then 
                    call unpack_ele_atomap(send_int, send_real, etype_buff, tag_ele_buff, size_ele_buff, basis_num_buff,&
                                           basis_type_buff, pb_node_buff, r_nodes, tag_atomap_buff, &
                                           type_atomap_buff, r_atomap_buff, vel_nodes) 
                else
                    call unpack_ele_atomap(send_int, send_real, etype_buff, tag_ele_buff, size_ele_buff, basis_num_buff,&
                                           basis_type_buff, pb_node_buff, r_nodes, tag_atomap_buff, &
                                           type_atomap_buff, r_atomap_buff) 
                end if

                need_ele(:) = .false.
                iatom_array(:) = 0

                !Loops through all interpolated atoms in the element and tag which ones belong to use
                ele_atomap_num = basis_num_buff*(size_ele_buff+1)**3
                do iatom = 1, ele_atomap_num
                    !NOTE: In this case iatom is each individual interpolated atom
                    if(in_block_bd(r_atomap_buff(:,iatom), pro_bd)) then
                        !If it's in the processor boundaries than mark it 
                        need_ele(iatom) = .true.
                        iatom_array(iatom) = 1
                    end if
              end do

              !Check to make sure all interpolated atoms have been grabbed
              seg_real = count(need_ele)

              call mpi_reduce(seg_real, i, 1, mpi_integer, &
                              mpi_sum, root, world, ierr)
              jatom_array(:) = 0
              call mpi_reduce(iatom_array(1:ele_atomap_num), jatom_array(1:ele_atomap_num), ele_atomap_num, mpi_integer, &
                              mpi_sum, root, world, ierr)

              if(rank == root) then
                    do iatom = 1, ele_atomap_num
                        if(jatom_array(iatom) /= 1) then
                            print *, 'Error: Jatom_array', jatom_array(iatom), &
                                     ' for iatom', iatom, ' should be 1 with position ', r_atomap_buff(:,iatom)
                            call mpi_abort(universe, 1, ierr)
                        end if
                    end do

                    if(i /= (ele_atomap_num)) then
                        print *, 'Error: Total number of atomaps', i, ' of element', &
                                  ie, ' should equal', ele_atomap_num
                        call mpi_abort(universe, 1, ierr)
                    end if
              end if

!             If the current element has atoms belonging to the processor than save information
!             regarding those atoms
              if(any(need_ele).eqv..true.) then
                    le = le + 1
                    
                    !Grow cg arrays if needed
                    if(le > ele_num_lr) call grow_cg_arrays(1)

                    lp = lp + ng_node(etype_buff)

                    if(lp > node_num_lr) call grow_cg_arrays(2, max(seg_num, seg_real))

                    latomap = latomap + seg_real

                    if(latomap > atomap_num_lr) call grow_cg_arrays(3, max(seg_num, seg_real))


                    !Save the element
                    ele_glob_id(le) = ie
                    tag_ele(le) = tag_ele_buff
                    size_ele(le) = size_ele_buff
                    etype(le) = etype_buff
                    basis_num(le) = basis_num_buff
                    basis_type(:,le) = basis_type_buff(:)
    
                    !Mark that we need it
                    who_ele(ie) = who_ele(ie) + 1
                    who_ele_long(pro_num*(ie-1)+rank+1) = 1
    
                    do inod = 1, ng_node(etype(le))
                        !Save the nodes of the elements
                        jp = lp - ng_node(etype(le)) + inod
                        node_cg(jp) = le
                        cg_node(inod, le) = jp

                        do ibasis = 1, basis_num(le)
                            do i = 1, 3

                                r(i, ibasis, jp) = r_nodes(i, ibasis, inod)

                                if(periodic) then
                                  pb_node(i, ibasis, jp) = pb_node_buff(i,ibasis,inod)
                                end if

                                if(need_vel) then
                                  vel(i, ibasis, jp) = vel_nodes(i, ibasis, inod)
                                end if

                          end do
                        end do

                    end do

                    !Runs through all atoms in the element, if the atom belongs to the current processor save
                    !the atom number to tag_atomap and the position to r_atomap, 
                
                    do iatom = 1, ele_atomap_num
                        if(need_ele(iatom).eqv..true.) then

                                !Save the atom information
                                jatomap = jatomap + 1
                                tag_atomap(jatomap) = tag_atomap_buff(iatom)
                                type_atomap(jatomap) = type_atomap_buff(iatom)
                                cg_atomap(iatom, le) = jatomap
                                r_atomap(:, jatomap) = r_atomap_buff(:,iatom)

                        end if
                    end do
                end if
            end do

!           debug

            if(jatomap /= latomap) then
                print *, 'Error: Total number of atomaps', jatomap, ' of rank', &
                          rank, ' should equal latomap', latomap
                call mpi_abort(universe, 1, ierr)
            end if

            !Total element count
            call mpi_reduce(le, i, 1, mpi_integer, mpi_sum, root, world, ierr)

            !Total node count
            call mpi_reduce(lp, j, 1, mpi_integer, mpi_sum, root, world, ierr) 

            !total atomap count
            call mpi_reduce(latomap, k, 1, mpi_integer, mpi_sum, root, world, ierr)

            if(rank == root) then

                if(i < ele_num) then
                    print *, 'Error: Total number of elements', i, ' should not be', &
                             ' smaller than ele_num', ele_num
                    call mpi_abort(universe, 1, ierr)

                else if(j < node_num) then
                    print *, 'Error: Total number of nodes', j, ' should not be', &
                             ' smaller than node_num', node_num
                    call mpi_abort(universe, 1, ierr)

                else if(k /= atomap_num) then

                    print *, 'Error: Total number of atomaps', k, ' should equal', &
                             ' atomap_num', atomap_num
                    call mpi_abort(universe, 1, ierr)
                end if

            end if
        end if
 
        ele_num_l = le
        node_num_l = lp
        atomap_num_l = latomap

        !debug who_ele

        if(sum(who_ele) /= ele_num_l) then
            print *, 'Error: sum(who_ele)', sum(who_ele), &
                     ' should equal', ele_num_l
            call mpi_abort(universe, 1, ierr)

        else if(maxval(who_ele) > 1) then
            print *, 'Error: Some elements are accounted for', &
                     ' more than once by rank', rank
            call mpi_abort(universe, 1, ierr)
        end if

        !debug who_ele_all

        allocate( who_ele_all(ele_num), stat = allostat)
        if(allostat /= 0) call alloc_error("Failure to allocate who_ele_all in scatter_cg",allostat)

        who_ele_all(:) = 0

        if(pro_num == 1) then
            who_ele_all(:) = who_ele(:)

        else
            call mpi_allreduce(who_ele, who_ele_all, ele_num, mpi_integer, mpi_sum, world, ierr)

        end if

        if(rank == root) then
            if(sum(who_ele_all) < ele_num) then
                print *, 'Error: sum(who_ele_all)', sum(who_ele_all), &
                         ' should not be smaller than', ele_num
                call mpi_abort(universe, 1, ierr)

            end if

        end if

        !who_ele_long_all
        allocate( who_ele_long_all(pro_num*ele_num), stat = allostat )
        if(allostat>0) call alloc_error('Failure to allocate who_ele_long_all', allostat)

        who_ele_long_all(:) = 0

        if(pro_num == 1) then
            who_ele_long_all(:) = who_ele_long(:)

        else
            call mpi_allreduce(who_ele_long, who_ele_long_all, pro_num*ele_num, &
                               mpi_integer, mpi_sum, world, ierr)
        end if

        !ele_shared_num and tag_ele_buff

        allocate( ele_id_buff(ele_num), stat = allostat )
        if(allostat /= 0) call alloc_error("Failure to allocate tag_ele_buff in scatter_cg", allostat)

        !Figure out which elements are shared
        ele_id_buff(:) = 0
        je = 0
        do ie = 1, ele_num
            if(who_ele_all(ie) > 1) then
                je = je + 1
                ele_id_buff(ie) = je
            end if
        end do
        ele_shared_num = je

        !ele_id_shared is an array which marks which elements are shared and among
        !how many other processors.

        allocate(ele_id_shared(ele_num_l), pro_shared_num(ele_num_l), stat = allostat) 
        if(allostat /= 0) call alloc_error("Failure to allocate ele_id_shared in scatter_cg", allostat)

        ele_id_shared(:) = 0
        pro_shared_num(:) = 0

        do ie = 1, ele_num_l
            le = ele_glob_id(ie)
            je = ele_id_buff(le)
            if(je /= 0) then
                ele_id_shared(ie) = je
            end if
        end do

        !who_has_ele distributes the elements. who_has_ele(ie) = true for the processor which is
        !in charge of calculating element number ie.
        allocate( who_has_ele(ele_num_l), who_element(ele_num), who_element_all(ele_num), stat = allostat) 
        if(allostat>0) call alloc_error("Failure allocating who_element arrays in scatter_cg", allostat)

        who_has_ele(:) = .false.
        who_element(:) = 0
        who_element_all(:) = 0
        irank_array(:) = [ (irank, irank = 1, pro_num) ]

        do ie = 1, ele_num_l

            le = ele_glob_id(ie)
            logic_array(:) = .false.

            do irank = 1, pro_num
                if(who_ele_long_all(pro_num*(le-1)+irank) == 1) then
                    logic_array(irank) = .true.
                end if
            end do

            pro_shared_num(ie) = count(logic_array)

            if((pro_shared_num(ie) > 1).and.(ele_id_shared(ie) == 0)) then
                print *, 'Error: When more than one processors share element', le, &
                         ' ele_id_shared(ie) should not be zero'
                call mpi_abort(universe, 1, ierr)
            end if

            jrank = maxloc(irank_array, 1, logic_array)
            if(rank == jrank - 1) then
                who_has_ele(ie) = .true.
                who_element(le) = 1
            end if

        end do

        !Get who_element_all
        if(pro_num == 1) then
             who_element_all(:) = who_element(:)
        else
          call mpi_reduce(who_element, who_element_all, ele_num, mpi_integer, mpi_sum, root, world, ierr)
        end if

        if(rank == root) then
            do ie = 1, ele_num

                if(who_element_all(ie) > 1) then
                    print *, 'Error: Element', ie, ' is included', &
                             ' in more than one,', who_element_all(ie), ' processors'
                    call mpi_abort(universe, 1, ierr)

                else if(who_element_all(ie) < 0) then
                    print *, 'Error: Element', ie, ' has a negative', &
                             ' who_element_all,', who_element_all(ie)
                    call mpi_abort(universe, 1, ierr)

                else if(who_element_all(ie) == 0) then
                    print *, 'Error: Element', ie, ' does not', &
                             ' belong to any processor'
                    call mpi_abort(universe, 1, ierr)

                end if
            end do
        end if

        if(pro_num == 1) then
            if(all(who_has_ele).eqv..false.) then
                print *, 'Error: All(who_has_ele) should be true', &
                         ' when there is only one processor'
                call mpi_abort(universe, 1, ierr)
            end if
        end if

        return

    end subroutine scatter_cg


    subroutine scatter_at
        !This subroutine scatters atomistic regions

        integer :: i, ia, ja, ka, la, seg_real, iatom, send_ni, send_nr, sent_nums, packet_num
        integer, dimension(3) :: info
        integer, dimension(seg_num) :: tag_buff, type_buff, iatom_array, jatom_array
        real(kind = wp), dimension(3,seg_num) :: r_buff, vel_buff
        integer, dimension(2*seg_num) :: send_int
        real(kind=wp), dimension(2*3*seg_num) :: send_real


        atom_num_lr = atom_num_l

        if(pro_num == 1) then
            !Code is only run on one processor (no MPI needed)
            type_atom(1:atom_num) = type_atom_scatter(:)

            do ia = 1, atom_num

                tag_atom(ia) = tag_atom_scatter(ia)

                if(periodic.eqv..true.) call cross_pb(r_atom_scatter(:,ia), info)
                r_atom(:, ia) = r_atom_scatter(:, ia)

            end do

            if(need_vel) then
                vel_atom(:, 1:atom_num) = vel_atom_scatter(:, :)
            end if

            atom_num_l = atom_num

        else
            !Code is run on multiple processors, MPI needed.

            ka = 0
            la = 0
            ia = 0
            tag_buff(:) = 0
            type_buff(:) = 0
            r_buff(:,:) = 0.0_wp
            vel_buff(:,:) = 0.0_wp
            atom_num_l = 0

            send_ni = 2*seg_num
            send_nr = 2*3*seg_num
            send_int=0
            send_real=0
            sent_nums = 0

            do while(sent_nums <  atom_num)

                if(rank == root) then

                    ia = ia +1
                    ka = ka + 1
                    tag_buff(ka) = tag_atom_scatter(ia)
                    type_buff(ka) = type_atom_scatter(ia)


                    r_buff(:,ka) = r_atom_scatter(:,ia)
                    if(periodic.eqv..true.) call cross_pb(r_buff(:,ka), info)
                    if(need_vel) vel_buff(:,ka) = vel_atom_scatter(:, ia)

                    if(ka == seg_num) then
                        seg_real = seg_num

                    else if(ia == atom_num) then
                        seg_real = mod(ia, seg_num)

                    else
                      cycle
                    end if

                    !Pack atom arrays
                    if (need_vel) then 
                        call pack_atoms(seg_num, tag_buff, type_buff, r_buff, send_int, send_real, vel_buff)
                    else
                        call pack_atoms(seg_num, tag_buff, type_buff, r_buff, send_int, send_real)
                    end if

                end if
                
                !Broadcast send arrays
                call mpi_bcast(send_int, send_ni, mpi_integer, root, world, ierr)
                call mpi_bcast(send_real, send_nr, mpi_wp, root, world, ierr)

                if(need_vel) then 
                    call unpack_atoms(seg_num, send_int, send_real, tag_buff, type_buff, r_buff, vel_buff)
                else
                    call unpack_atoms(seg_num, send_int, send_real, tag_buff, type_buff, r_buff)
                end if

                packet_num = 0
                iatom_array(:) = 0
                ja = 0
                do ka = 1, seg_num

                    sent_nums = sent_nums + 1

                    !First check to make sure the atom type is real
                    if(sent_nums>atom_num) then 
                        exit
                    else if(type_buff(ka) == 0) then 
                        print *, "Atom_type cannot be 0 in scatter_at for atom number ", sent_nums-seg_num+ka
                        call mpi_abort(1, universe, ierr)
                    end if

                    packet_num = packet_num + 1
                    !checks to see if the atom is within the processor boundaries

                    if(in_block_bd(r_buff(:,ka), pro_bd)) then

                        ja = ja + 1
                        la = la + 1
                        if(la > atom_num_lr) call grow_at_arrays

                        !Pull out values of properties belonging to atoms in the current processor
                        tag_atom(la) = tag_buff(ka)
                        type_atom(la) = type_buff(ka)
                        r_atom(:, la) = r_buff(:,ka)
                        if(need_vel) vel_atom(:, la) = vel_buff(:,ka)

                        !Iatom array just makes sure everything is only assigned once
                        iatom_array(ka) = 1

                        atom_num_l = atom_num_l + 1 
                    end if
                end do

                !Check to make sure everything was sent out correctly
                call mpi_reduce(ja, i, 1, mpi_integer, mpi_sum, root, world, ierr)
                call mpi_reduce(iatom_array(1:packet_num), jatom_array(1:packet_num), packet_num, mpi_integer, mpi_sum, root,&
                                world, ierr)

                if(rank == root) then 
                    do iatom = 1, packet_num
                        if (jatom_array(iatom) /= 1) then
                            print *, box_bd, r_buff(:,iatom)
                            print *, "Error: Jatom array should be 1 not ", jatom_array(iatom), " for iatom ", iatom
                            call mpi_abort(universe, 1, ierr)
                        end if
                    end do
                    if (i/=packet_num) then 
                        print *, packet_num, " atoms were sent but ", i, " atoms were claimed by the processors."
                        call mpi_abort(1, universe, ierr)
                    end if

                end if
                

                if(rank == root) ka = 0

            end do

        end if


        !Now make sure everything has been properly distributed 
        if(pro_num == 1) then
            i = atom_num_l
        else
            call mpi_reduce(atom_num_l, i, 1, mpi_integer, mpi_sum, root, world, ierr)

        end if

        if(rank == root) then
            if(i /= atom_num) then
                print *, 'Error: Total number of atoms', i, ' should equal', &
                         ' atom_num', atom_num
                call mpi_abort(universe, 1, ierr)
            end if

        end if

        return
    end subroutine scatter_at

    subroutine ghost_cg
        !Share ghost interpolated atoms from elements

        integer :: i, j, k, send_n, recv_n, iatomap, jatomap, latomap, &
                   seg_num_real, itag_t, itag_r, itag_n, itag_ty, &
                   ireq_t, ireq_r, ireq_n, ireq_ty, &
                   send_rank, recv_rank, send_ori, &
                   atomap_num_s, atomap_num_r, accept_count, accept_count_r

        logical :: send_l, recv_l

        integer, dimension(3) :: send_coords, recv_coords

        integer, dimension(mpi_status_size) :: mstatus

        real(kind = wp), dimension(3) :: r_in

        integer, allocatable :: &
          dir_atomap(:), dir_atomap_array(:), &
          tag_send_buff_array(:), tag_send_buff_ori(:), &
          tag_send_buff(:), tag_recv_buff(:), &
          tag_atomap_array(:), tag_send_array(:), tag_recv_array(:), &
          type_send_buff_array(:), type_send_buff_ori(:), &
          type_send_buff(:), type_recv_buff(:), &
          type_atomap_array(:), type_send_array(:), type_recv_array(:), &
          list_atomap_array(:, :, :), pos_send_buff_ori(:), &
          pos_send_buff(:), pos_send_array(:), accept_send(:), accept_recv(:), &
          accept_array(:)

        real(kind = wp), allocatable :: &
          r_send_buff_array(:), r_send_buff_ori(:), &
          r_send_buff(:), r_recv_buff(:), &
          r_atomap_array(:, :), r_send_array(:), r_recv_array(:)

        !Set seg_num_real which is basically array growth size
        seg_num_real = seg_num

        !Allocate original send_buff arrays
        allocate(tag_send_buff_ori(seg_num_real), &
                 type_send_buff_ori(seg_num_real),&
                 pos_send_buff_ori(seg_num_real),&
                 r_send_buff_ori(3*seg_num_real), &
                 stat = allostat)
        if(allostat /= 0) call alloc_error("Failure to allocate tag/r_send_buff_ori",allostat)

        !Initialize some variales
        send_n = 0
        tag_send_buff_ori(:) = 0
        pos_send_buff_ori(:) = 0
        type_send_buff_ori(:) = 0
        r_send_buff_ori(:) = 0.0_wp

        !Figure out which atomaps need to be sent
        do iatomap = 1, atomap_num_l
            r_in(:) = r_atomap(:, iatomap)

            !Check position of the atomap to see if it needs to be sent
            if(.not.in_block_bd(r_in, pro_bd_in)) then
                send_n = send_n + 1

                !Grow arrays if needed
                if(send_n > seg_num_real) then
                    allocate( tag_send_buff_array(seg_num_real+seg_num), &
                              pos_send_array(seg_num_real + seg_num), &  
                              type_send_buff_array(seg_num_real+seg_num), &
                              r_send_buff_array(3*(seg_num_real+seg_num)), &
                              stat = allostat )
                    if(allostat /= 0) call alloc_error("Failure to allocate tag/r_send_buff_array", allostat)

                    tag_send_buff_array(1:seg_num_real) = tag_send_buff_ori(:)
                    tag_send_buff_array(seg_num_real+1:) = 0
                    call move_alloc(tag_send_buff_array, tag_send_buff_ori)

                    pos_send_array(1:seg_num_real) = pos_send_buff_ori(:)
                    pos_send_array(seg_num_real+1:) = 0
                    call move_alloc(pos_send_array, pos_send_buff_ori)

                    type_send_buff_array(1:seg_num_real) = type_send_buff_ori(:)
                    type_send_buff_array(seg_num_real+1:) = 0
                    call move_alloc(type_send_buff_array, type_send_buff_ori)

                    r_send_buff_array(1:3*seg_num_real) = r_send_buff_ori(:)
                    r_send_buff_array(3*seg_num_real+1:) = 0.0_wp
                    call move_alloc(r_send_buff_array, r_send_buff_ori)
                    seg_num_real = seg_num_real + seg_num
                end if

                !Build send_buff_ori
                tag_send_buff_ori(send_n) = tag_atomap(iatomap)
                pos_send_buff_ori(send_n) = iatomap
                type_send_buff_ori(send_n) = type_atomap(iatomap)
                do i = 1, 3
                  r_send_buff_ori(3*(send_n-1)+i) = r_in(i)
                end do
            end if
        end do

        allocate(tag_send_buff(send_n), type_send_buff(send_n), r_send_buff(3*send_n), pos_send_buff(send_n), stat = allostat)
        if (allostat /= 0) call alloc_error("Failure to allocate tag/r_send_buff", allostat)

        !Initialize allocated buffs and recv_n
        recv_n = 0
        tag_send_buff(:) = 0
        type_send_buff(:) = 0
        r_send_buff(:) = 0.0_wp

        !Set recv_num
        if(pro_num == 1) then
          recv_n = send_n
        else
          call mpi_allreduce(send_n, recv_n, 1, mpi_integer, mpi_max, grid_comm, ierr)
        end if

        !Allocate arrays for receiving data
        allocate(tag_recv_buff(recv_n), type_recv_buff(recv_n), r_recv_buff(3*recv_n), stat = allostat)
        if(allostat /= 0) call alloc_error("Failure to allocate tag/r_recv_buff", allostat)

        !Also allocate arrays for communicating which atomaps were received
        allocate(accept_send(recv_n), accept_recv(recv_n), stat=allostat)
        if(allostat /= 0) call alloc_error("Failure to allocate accept arrays", allostat)

        !Allocate receive buff
        tag_recv_buff(:) = 0
        type_recv_buff(:) = 0
        r_recv_buff(:) = 0.0_wp
        accept_send(:) = 0
        accept_recv(:) = 0 

        !Initialize tags used for communications
        itag_r = 1
        itag_n = 2
        itag_t = 3
        itag_ty = 4
        ireq_r = 0
        ireq_n = 0
        ireq_t = 0
        ireq_ty = 0

        !Initialize some more variables
        send_rank = grank
        recv_rank = grank
        list_atomap_logic(:, :) = .true.

        !Allocate variables and initialize
        allocate(dir_atomap(atomap_num_lr), stat = allostat) 
        if(allostat /= 0) call alloc_error("Failure to allocate dir_atomap", allostat)

        !dir_atomap is the direction to send each atomap
        dir_atomap(:) = 4

        !Initialize the atomap counts where atomap_num_s is the send_count and atomap_num_r is the recv count
        jatomap = atomap_num_l
        atomap_num_s = send_n
        send_ori = send_n
        atomap_num_r = recv_n

        !For these loops, j is the box direction to communicate along and k is to direct either positive or negative comm.
        !The ghost approach here is described in https://doi.org/10.1016/j.commatsci.2017.11.051

        do j = 1, 3
            do k = 1, 2
                
                send_coords(:) = grid_coords(:)
                recv_coords(:) = grid_coords(:)
                send_l = .true.
                recv_l = .true.

                if(num_pro(j) == 1) then

                    !If we have one processor in the dimension then we don't do a normal send.
                    send_l = .false.
                    !Only send if periodic in that dimension
                    if(period(j).eqv..true.) then
                        recv_l = .true.
                        tag_recv_buff(1:send_ori) = tag_send_buff_ori(1:send_ori)
                        type_recv_buff(1:send_ori) = type_send_buff_ori(1:send_ori)
                        r_recv_buff(1:3*send_ori) = r_send_buff_ori(1:3*send_ori)
                        pos_send_buff(1:send_ori) = pos_send_buff_ori(1:send_ori)
                        latomap = send_ori

                        !Now loop over all the ghost atoms
                        do iatomap = atomap_num_l+1, jatomap
                            !If the ghost was received in the previous communication step add it to the send_Buff
                            if(dir_atomap(iatomap) < j) then
                                latomap = latomap + 1
                                !Resize arrays if needed
                                if(latomap > atomap_num_r) then

                                    allocate(tag_recv_array(atomap_num_r+seg_num), &
                                             type_recv_array(atomap_num_r+seg_num), &
                                             r_recv_array(3*(atomap_num_r+seg_num)), &
                                             stat = allostat)

                                    if(allostat /= 0) call alloc_error("Failure to allocate tag/r_recv_array", allostat)

                                    tag_recv_array(:) = 0
                                    tag_recv_array(1:atomap_num_r) = tag_recv_buff(:)
                                    call move_alloc(tag_recv_array, tag_recv_buff)

                                    type_recv_array(:) = 0
                                    type_recv_array(1:atomap_num_r) = type_recv_buff(:)
                                    call move_alloc(type_recv_array, type_recv_buff)

                                    r_recv_array(:) = 0.0_wp
                                    r_recv_array(1:3*atomap_num_r) = r_recv_buff(:)
                                    call move_alloc(r_recv_array, r_recv_buff)
                                    atomap_num_r = atomap_num_r + seg_num

                                end if

                                !Save the send atoms to the recieved array
                                tag_recv_buff(latomap) = tag_atomap(iatomap)
                                type_recv_buff(latomap) = type_atomap(iatomap)
                                do i = 1, 3
                                  r_recv_buff(3*(latomap-1)+i) = r_atomap(i, iatomap)
                                end do

                                !Add values to pos_send_buff
                                if(latomap > size(pos_send_buff)) then 
                                    allocate(pos_send_array(size(pos_send_buff)+ seg_num)) 
                                    pos_send_array(1:size(pos_send_buff)) = pos_send_buff
                                    pos_send_array(size(pos_send_buff)+1:) = 0
                                    call move_alloc(pos_send_array, pos_send_buff)
                                end if
                                pos_send_buff(latomap) = iatomap
                            end if
                        end do

                        recv_n = latomap
                        !Apply periodic boundaries 
                        do iatomap = 1, recv_n
                          r_recv_buff(3*(iatomap-1)+j) = r_recv_buff(3*(iatomap-1)+j) + (-1)**k * box_length(j)
                        end do

                    else
                        !Otherwise
                        recv_l = .false.
                        recv_n = 0
                    end if


                !If we have more than one processor in the dimension
                else

                    send_coords(j) = grid_coords(j) - (-1)**k
                    recv_coords(j) = grid_coords(j) + (-1)**k
                    !If periodic then we send to the next processor in every situation
                    if(period(j).eqv..true.) then
                        call mpi_cart_rank(grid_comm, send_coords, send_rank, ierr)
                        call mpi_cart_rank(grid_comm, recv_coords, recv_rank, ierr)
                    !Otherwise the first and last ones along each dimension have special rules as they don't send to one processro
                    !and don't receive from one processor at a specific step
                    else
                      if(k == 1) then
                          !Don't send if it's the last one and k = 1
                            if(grid_coords(j) < num_pro(j) - 1) then 
                                call mpi_cart_rank(grid_comm, send_coords, send_rank, ierr)
                            else
                                send_l = .false.
                                send_rank = -1
                            end if

                            !Don't receive if first one and k = 1
                            if(grid_coords(j) > 0) then
                                call mpi_cart_rank(grid_comm, recv_coords, recv_rank, ierr)
                            else
                                recv_l = .false.
                                recv_rank = -1
                            end if

                            else if(k == 2) then
                            !Don't send if first one and k = 2
                            if(grid_coords(j) > 0) then 
                                call mpi_cart_rank(grid_comm, send_coords, send_rank, ierr)
                            else
                                send_l = .false.
                                send_rank = -1
                            end if

                            !Don't receive if last one and k = 2
                            if(grid_coords(j) < num_pro(j) - 1) then
                                call mpi_cart_rank(grid_comm, recv_coords, recv_rank, ierr)
                            else
                                recv_l = .false.
                                recv_rank = -1
                            end if
                        end if 
                    end if
                    
                    !prepare send_buff
                    if(send_l.eqv..true.) then

                        tag_send_buff(1:send_ori) = tag_send_buff_ori(1:send_ori)
                        pos_send_buff(1:send_ori) = pos_send_buff_ori(1:send_ori)
                        type_send_buff(1:send_ori) = type_send_buff_ori(1:send_ori)
                        r_send_buff(1:3*send_ori) = r_send_buff_ori(1:3*send_ori)
                        latomap = send_ori

                        !Loop over all received ghost atoms
                        do iatomap = atomap_num_l+1, jatomap

                            !if we previously received the ghosts than we add them to the send_buff
                            if(dir_atomap(iatomap) < j) then

                                latomap = latomap + 1
                                !Grow arrays if eneded
                                if(latomap > atomap_num_s) then

                                  allocate( &
                                    tag_send_array(atomap_num_s+seg_num), &
                                    type_send_array(atomap_num_s+seg_num), &
                                    r_send_array(3*(atomap_num_s+seg_num)), &
                                    stat = allostat &
                                  )

                                  if(allostat /= 0) call alloc_error("Failure to allocate tag/r_send_array", allostat)

                                  tag_send_array(:) = 0
                                  tag_send_array(1:atomap_num_s) = tag_send_buff(:)
                                  call move_alloc(tag_send_array, tag_send_buff)

                                  type_send_array(:) = 0
                                  type_send_array(1:atomap_num_s) = type_send_buff(:)
                                  call move_alloc(type_send_array, type_send_buff)

                                  r_send_array(:) = 0.0_wp
                                  r_send_array(1:3*atomap_num_s) = r_send_buff(:)
                                  call move_alloc(r_send_array, r_send_buff)

                                  atomap_num_s = atomap_num_s + seg_num

                                end if

                                tag_send_buff(latomap) = tag_atomap(iatomap)
                                type_send_buff(latomap) = type_atomap(iatomap)
                                do i = 1, 3
                                  r_send_buff(3*(latomap-1)+i) = r_atomap(i, iatomap)
                                end do

                                if(latomap > size(pos_send_buff)) then 

                                    allocate(pos_send_array(size(pos_send_buff)+seg_num))
                                    pos_send_array(1:size(pos_send_buff)) = pos_send_buff
                                    pos_send_array(size(pos_send_buff)+1:) = 0
                                    call move_alloc(pos_send_array, pos_send_buff)
                                end if
                                pos_send_buff(latomap) = iatomap

                            end if
                        end do

                        send_n = latomap

!                       update r_send_buff for pbc
                        if(period(j).eqv..true.) then
                            do iatomap = 1, send_n
                                if(k == 1) then
                                    if(grid_coords(j) == num_pro(j) - 1) then
                                        r_send_buff(3*(iatomap-1)+j) = r_send_buff(3*(iatomap-1)+j) - box_length(j)
                                    end if
                                end if

                                if(k == 2) then

                                    if(grid_coords(j) == 0) then

                                        r_send_buff(3*(iatomap-1)+j) = r_send_buff(3*(iatomap-1)+j) + box_length(j)
                                    end if
                                end if
                            end do
                        end if
                    end if

                    !send/recv number
                    if(recv_l.eqv..true.) then
                      call mpi_irecv(recv_n, 1, mpi_integer, recv_rank, itag_n, grid_comm, ireq_n, ierr)
                    end if

                    if(send_l.eqv..true.) then
                      call mpi_send(send_n, 1, mpi_integer, send_rank, itag_n, grid_comm, ierr)
                    end if

                    if(recv_l.eqv..true.) then
                      call mpi_wait(ireq_n, mstatus, ierr)
                    end if

                    !resize recv_buff if needed
                    if(recv_l.eqv..true.) then
                        if(recv_n > atomap_num_r) then
                            deallocate( tag_recv_buff, type_recv_buff, r_recv_buff, stat = deallostat)
                            if(deallostat /= 0) call alloc_error("Failure to deallocate tag/r_recv_buff", allostat)
                            allocate( tag_recv_buff(recv_n), r_recv_buff(3*recv_n), type_recv_buff(recv_n), stat = allostat)
                            if(allostat /= 0) call alloc_error("Failure to allocate tag/r_recv_buff", allostat) 
                            !Initialize arrays
                            tag_recv_buff(:) = 0
                            type_recv_buff(:) = 0
                            r_recv_buff(:) = 0.0_wp
                            atomap_num_r = recv_n
                            
                            !Reallocate accept arrays
                            deallocate(accept_send)
                            allocate(accept_send(recv_n))
                            accept_send(:) = 0
                        end if
                    end if

                    !send/recv tag 
                    if(recv_l.eqv..true.) then
                        tag_recv_buff=0
                        call mpi_irecv(tag_recv_buff(1:recv_n), recv_n, mpi_integer, recv_rank, itag_t, grid_comm, ireq_t, ierr)
                    end if
                    if(send_l.eqv..true.) then
                        call mpi_send(tag_send_buff(1:send_n), send_n, mpi_integer, send_rank, itag_t, grid_comm, ierr)
                    end if
                    if(recv_l.eqv..true.) then
                      call mpi_wait(ireq_t, mstatus, ierr)
                    end if

                    !send/recv type
                    if(recv_l.eqv..true.) then
                        type_recv_buff=0
                        call mpi_irecv(type_recv_buff(1:recv_n), recv_n, mpi_integer, recv_rank, itag_ty, grid_comm, ireq_ty, ierr)
                    end if
                    if(send_l.eqv..true.) then
                        call mpi_send(type_send_buff(1:send_n), send_n, mpi_integer, send_rank, itag_ty, grid_comm, ierr)
                    end if
                    if(recv_l.eqv..true.) then
                        call mpi_wait(ireq_ty, mstatus, ierr)
                    end if

                    !Send and receive atomap positions
                    if(recv_l.eqv..true.) then
                        r_recv_buff=0.0_wp
                        call mpi_irecv(r_recv_buff(1:3*recv_n), 3*recv_n, mpi_wp, recv_rank, itag_r, grid_comm, ireq_r, ierr)
                    end if
                    if(send_l.eqv..true.) then
                        call mpi_send(r_send_buff(1:3*send_n), 3*send_n, mpi_wp, send_rank, itag_r, grid_comm, ierr)
                    end if
                    if(recv_l.eqv..true.) then
                        call mpi_wait(ireq_r, mstatus, ierr)
                    end if
                end if

                !check if each core needs the received atomaps
                if(recv_l.eqv..true.) then

                    accept_count = 0
                    accept_send(:) = 0
                    do iatomap = 1, recv_n

                        do i = 1, 3
                            r_in(i) = r_recv_buff(3*(iatomap-1)+i)
                        end do

                        !If it's in our processor boundaries then add it
                        if(in_block_bd(r_in, pro_bd_out)) then

                            jatomap = jatomap + 1

                            !Grow arrays if needed
                            if(jatomap > atomap_num_lr) then

                                allocate(tag_atomap_array(atomap_num_lr+seg_num), & 
                                         type_atomap_array(atomap_num_lr+seg_num), &
                                         r_atomap_array(3, atomap_num_lr+seg_num), &
                                         dir_atomap_array(atomap_num_lr+seg_num), &
                                         stat = allostat)

                                if(allostat /= 0) call alloc_error("Failure to allocate tag/r/dir_atomap_array", allostat)


                                tag_atomap_array(1:atomap_num_lr) = tag_atomap(:)
                                tag_atomap_array(atomap_num_lr+1:) = 0
                                call move_alloc(tag_atomap_array, tag_atomap)

                                type_atomap_array(1:atomap_num_lr) = type_atomap(:)
                                type_atomap_array(atomap_num_lr+1:) = 0
                                call move_alloc(type_atomap_array, type_atomap)

                                r_atomap_array(:, 1:atomap_num_lr) = r_atomap(:, :)
                                r_atomap_array(:, atomap_num_lr+1:) = 0.0_wp
                                call move_alloc(r_atomap_array, r_atomap)

                                dir_atomap_array(1:atomap_num_lr) = dir_atomap(:)
                                dir_atomap_array(atomap_num_lr+1:) = 0
                                call move_alloc(dir_atomap_array, dir_atomap)

                                atomap_num_lr = atomap_num_lr + seg_num
                            end if

                            tag_atomap(jatomap) = tag_recv_buff(iatomap)
                            type_atomap(jatomap) = type_recv_buff(iatomap)
                            r_atomap(:, jatomap) = r_in(:)
                            dir_atomap(jatomap) = j

                            !Now mark this atom as accepted
                            accept_count = accept_count + 1
                            if(accept_count > size(accept_send)) then 
                                allocate(accept_array(size(accept_send) + seg_num))
                                accept_array(1:size(accept_send)) = accept_send
                                accept_array(size(accept_send)+1:) = 0
                                call move_alloc(accept_array, accept_send)
                            end if

                            accept_send(accept_count) = iatomap
                        end if
                    end do
                end if
                list_atomap_logic(k, j) = send_l

                if((num_pro(j) == 1).and.(period(j))) then
                    !Resize if needed
                    if(accept_count > size(accept_recv)) then 
                        deallocate(accept_recv)
                        allocate(accept_recv(accept_count))
                    end if                    
                    accept_recv(1:accept_count) = accept_send(1:accept_count)
                    accept_count_r = accept_count
                else
                    !Now if we sent atoms this turn we have to receive the accepted ghost atoms to build our list, 
                    !First get the accepted counts
                    if(send_l) call mpi_irecv(accept_count_r, 1, mpi_integer, send_rank, 15, grid_comm, ireq_n, ierr)
                    if(recv_l) call mpi_send(accept_count, 1, mpi_integer, recv_rank, 15, grid_comm, ierr)
                    if(send_l) call mpi_wait(ireq_n, mstatus, ierr)

                    if(send_l) then 
                        !Resize if needed
                        if(accept_count_r > size(accept_recv)) then 
                            deallocate(accept_recv)
                            allocate(accept_recv(accept_count_r))
                        end if

                        !Now receive the accepted atomaps
                        call mpi_irecv(accept_recv, accept_count_r, mpi_integer, send_rank, 16, grid_comm, ireq_n, ierr)
                    end if

                    if(recv_l) call mpi_send(accept_send(1:accept_count), accept_count, mpi_integer, &
                                             recv_rank, 16, grid_comm, ierr)

                    if(send_l) call mpi_wait(ireq_n, mstatus, ierr)

                end if
                    
                if(send_l.or.((num_pro(j)==1).and.period(j))) then 
                    !Now build the list of sent atoms 
                    list_atomap_num(k,j) = accept_count_r

                    !initialize list_atomap which contains information on which atoms were sent on every dimension
                    if(.not.allocated(list_atomap)) then 
                        allocate(list_atomap(accept_count, 2, 3), stat = allostat)
                        if (allostat>0) call alloc_error("Failure to allocate list_atomap", allostat)
                    end if

                    
                    if (accept_count_r > size(list_atomap(:,k,j))) then 
                        allocate(list_atomap_array(accept_count_r, 2, 3))
                        list_atomap_array(1:size(list_atomap(:,k,j)),:,:) = list_atomap(1:size(list_atomap(:,k,j)),:,:)
                        list_atomap_array(size(list_atomap(:,k,j))+1,:,:) = 0
                        call move_alloc(list_atomap_array, list_atomap)
                    end if

                    do i = 1, accept_count_r
                        list_atomap(i, k, j) = pos_send_buff(accept_recv(i))
                    end do 
                else
                    list_atomap_num(k,j) = 0
                end if
            end do
        end do

        atomap_num_lg = jatomap

        !debug to make sure code is correct
        if(atomap_num_lr /= size(r_atomap, 2)) then
            print *, 'Error: Wrong atomap_num_lr', atomap_num_lr, ' which should equal size(r_atomap, 2)', size(r_atomap, 2)
            call mpi_abort(mpi_comm_world, 1, ierr)

        else if(atomap_num_lr /= size(tag_atomap)) then
            print *, 'Error: Wrong atomap_num_lr', atomap_num_lr, ' which should equal size(tag_atomap)', size(tag_atomap)
            call mpi_abort(mpi_comm_world, 1, ierr)

        else if(atomap_num_lr /= size(dir_atomap)) then
            print *, 'Error: Wrong atomap_num_lr', atomap_num_lr, ' which should equal size(dir_atomap)', size(dir_atomap)
            call mpi_abort(mpi_comm_world, 1, ierr)

        else if(atomap_num_lg > atomap_num_lr) then
            print *, 'Error: Wrong atomap_num_lg', atomap_num_lg, ' which is larger than', atomap_num_lr
            call mpi_abort(mpi_comm_world, 1, ierr)
        end if


!        list_atomap(:, :, :) = 0
!        list_atomap_num(:, :) = 0
!        seg_num_atomap = size(list_atomap, 1)
!
!        do j = 1, 3
!            do k = 1, 2
!
!                send_n = 0
!                pro_bd_temp(:) = pro_bd_kj(:, k, j)
!                do iatomap = 1, atomap_num_lg
!
!                    !If we send the atom then 
!                    r_in(:) = r_atomap(:, iatomap)
!                    if(in_block_bd(r_in, pro_bd_temp)) then
!
!                        send_n = send_n + 1
!
!                        if(send_n > seg_num_atomap) then
!                            allocate(list_atomap_array(seg_num_atomap+seg_num, 2, 3), stat = allostat)
!                            if(allostat /= 0) call alloc_error("Failure to allocate list_atomap_array", allostat)
!
!                            list_atomap_array(1:seg_num_atomap, :, :) = list_atomap(:, :, :)
!                            list_atomap_array(seg_num_atomap+1:, :, :) = 0
!                            call move_alloc(list_atomap_array, list_atomap)
!                            seg_num_atomap = seg_num_atomap + seg_num
!                        end if
!
!                        list_atomap(send_n, k, j) = iatomap
!                    end if
!                end do
!                list_atomap_num(k, j) = send_n
!            end do
!        end do

        return
    end subroutine ghost_cg
    
    subroutine ghost_at
        !This code send ghost atoms to processors that need them
        integer :: i, j, k, send_n, recv_n, ia, ja, la, &
                   seg_num_real, itag_t, itag_r, itag_n, itag_ty,&
                   ireq_t, ireq_r, ireq_n, ireq_ty,&
                   send_rank, recv_rank, send_ori, &
                   atom_num_s, atom_num_r, atom_num_ln, &
                   accept_count, accept_count_r

        logical :: send_l, recv_l

        integer, dimension(3) :: send_coords, recv_coords

        integer, dimension(mpi_status_size) :: mstatus

        real(kind = wp), dimension(3) :: r_in

        integer, allocatable :: dir_atom(:), dir_atom_array(:), &
                                tag_send_buff_array(:), tag_send_buff_ori(:), &
                                tag_send_buff(:), tag_recv_buff(:), &
                                tag_send_array(:), tag_recv_array(:), &
                                type_send_buff_array(:), type_send_buff_ori(:), &
                                type_send_buff(:), type_recv_buff(:), &
                                type_send_array(:), type_recv_array(:), &
                                list_atom_array(:, :, :), accept_array(:), &
                                pos_send_buff_ori(:), pos_send_buff(:), pos_send_array(:), &
                                accept_send(:), accept_recv(:)

        real(kind = wp), allocatable :: r_send_buff_array(:), r_send_buff_ori(:), &
                                        r_send_buff(:), r_recv_buff(:), &
                                        r_send_array(:), r_recv_array(:), &
                                        vel_array(:, :)

        !Initialize size variables and allocate variables
        seg_num_real = seg_num
        allocate(tag_send_buff_ori(seg_num_real), &
                 type_send_buff_ori(seg_num_real),&
                 pos_send_buff_ori(seg_num_real), &
                 r_send_buff_ori(3*seg_num_real), &
                 stat = allostat )
        if(allostat /= 0) call alloc_error("Failure to allocate tag/r_send_buff_ori", allostat)

        !Initialize other variables
        send_n = 0
        tag_send_buff_ori(:) = 0
        type_send_buff_ori(:) = 0
        pos_send_buff_ori(:) = 0 
        r_send_buff_ori(:) = 0.0_wp

        !Build atom send lists
        do ia = 1, atom_num_l

            r_in(:) = r_atom(:, ia)

            !Check to see if this needs to be sent
            if(.not.in_block_bd(r_in, pro_bd_in)) then

                !If we need to save it then check if we need to resize the arrays
                send_n = send_n + 1
                if(send_n > seg_num_real) then
                    allocate(tag_send_buff_array(seg_num_real+seg_num), &
                             type_send_buff_array(seg_num_real+seg_num), &
                             pos_send_array(seg_num_real + seg_num), &
                             r_send_buff_array(3*(seg_num_real+seg_num)), &
                             stat = allostat)
                    if(allostat /= 0) call alloc_error("Failure to allocate tag/r_send_buff_array", allostat)

                    tag_send_buff_array(1:seg_num_real) = tag_send_buff_ori(:)
                    tag_send_buff_array(seg_num_real+1:) = 0
                    call move_alloc(tag_send_buff_array, tag_send_buff_ori)

                    pos_send_array(1:seg_num_real) = pos_send_buff_ori(:)
                    pos_send_array(seg_num_real+1:) = 0 
                    call move_alloc(pos_send_array, pos_send_buff_ori)

                    type_send_buff_array(1:seg_num_real) = type_send_buff_ori(:)
                    type_send_buff_array(seg_num_real+1:) = 0
                    call move_alloc(type_send_buff_array, type_send_buff_ori)

                    r_send_buff_array(1:3*seg_num_real) = r_send_buff_ori(:)
                    r_send_buff_array(3*seg_num_real+1:) = 0.0_wp
                    call move_alloc(r_send_buff_array, r_send_buff_ori)

                    seg_num_real = seg_num_real + seg_num
                end if

                tag_send_buff_ori(send_n) = tag_atom(ia)
                type_send_buff_ori(send_n)= type_atom(ia)
                pos_send_buff_ori(send_n) = ia
                do i = 1, 3
                    r_send_buff_ori(3*(send_n-1)+i) = r_in(i)
                end do
            end if
        end do

        !Allocate send_buff variable which changes upon each send step
        allocate(tag_send_buff(send_n), &
                 type_send_buff(send_n),&
                 pos_send_buff(send_n), &
                 r_send_buff(3*send_n), &
                 stat = allostat)
        if(allostat /= 0) call alloc_error("Failure to allocate tag/r_send_buff", allostat)

        !Initialize new arrays
        tag_send_buff(:) = 0
        type_send_buff(:)= 0
        pos_send_buff(:) = 0 
        r_send_buff(:) = 0.0_wp

        if(pro_num == 1) then
            recv_n = send_n
        else
            call mpi_allreduce(send_n, recv_n, 1, mpi_integer, mpi_max, grid_comm, ierr)
        end if

        !Allocate arrays for receiving data
        allocate(tag_recv_buff(recv_n), type_recv_buff(recv_n), r_recv_buff(3*recv_n), stat = allostat ) 
        if(allostat /= 0) call alloc_error("Failure to allocate tag/r_recv_buff", allostat)

        !Also allocate arrays for communicating which atoms were received
        allocate(accept_send(recv_n), accept_recv(recv_n), stat=allostat)
        if(allostat /= 0) call alloc_error("Failure to allocate accept arrays", allostat)
        accept_send = 0
        accept_recv = 0


        !Initalize new variables and tags used for communications
        tag_recv_buff(:) = 0
        type_recv_buff(:)= 0
        r_recv_buff(:) = 0.0_wp
        itag_r = 1
        itag_n = 2
        itag_t = 3
        itag_ty= 4
        ireq_r = 0
        ireq_n = 0
        ireq_t = 0
        ireq_ty= 0
        send_rank = grank
        recv_rank = grank
        list_atom_logic(:, :) = .true.

        allocate(dir_atom(atom_num_lr), stat = allostat)
        if(allostat /= 0) call alloc_error("Failure to allocate dir_atom", allostat)

        dir_atom(:) = 4
        ja = atom_num_l
        atom_num_s = send_n
        send_ori = send_n
        atom_num_r = recv_n

        ! Check comments on ghost_cg for more information on the ghost scheme
        ! k == 1: positive
        ! k == 2: negative

        do j = 1, 3
            do k = 1, 2
                !Get send and receive coordinates
                send_coords(:) = grid_coords(:)
                recv_coords(:) = grid_coords(:)
                send_l = .true.
                recv_l = .true.

                !If only one processor in the dimension then we don't send ghosts unless periodic in that dim
                if(num_pro(j) == 1) then

                    send_l = .false.
                    if(period(j).eqv..true.) then

                        recv_l = .true.
                        tag_recv_buff(1:send_ori) = tag_send_buff_ori(1:send_ori)
                        type_recv_buff(1:send_ori)= type_send_buff_ori(1:send_ori)
                        r_recv_buff(1:3*send_ori) = r_send_buff_ori(1:3*send_ori)
                        pos_send_buff(1:send_ori) = pos_send_buff_ori(1:send_ori)
                        la = send_ori

                        do ia = atom_num_l+1, ja
                            !Check to add atoms which were previously sent 
                            if(dir_atom(ia) < j) then

                                la = la + 1

                                !Grow arrays if needed
                                if(la > atom_num_r) then
                                    allocate( tag_recv_array(atom_num_r+seg_num), &
                                              type_recv_array(atom_num_r+seg_num), &  
                                              r_recv_array(3*(atom_num_r+seg_num)), &
                                              stat = allostat )
                                    if(allostat /= 0) call alloc_error("Failure to allocate tag/r_recv_array", allostat)

                                    tag_recv_array(1:atom_num_r) = tag_recv_buff(:)
                                    tag_recv_array(atom_num_r+1:) = 0
                                    call move_alloc(tag_recv_array, tag_recv_buff)

                                    type_recv_array(1:atom_num_r) = type_recv_buff(:)
                                    type_recv_array(atom_num_r+1:) = 0
                                    call move_alloc(type_recv_array, type_recv_buff)

                                    r_recv_array(1:3*atom_num_r) = r_recv_buff(:)
                                    r_recv_array(3*atom_num_r+1:) = 0.0_wp
                                    call move_alloc(r_recv_array, r_recv_buff)

                                    atom_num_r = atom_num_r + seg_num
                                end if

                                tag_recv_buff(la) = tag_atom(ia)
                                type_recv_buff(la)= type_atom(ia)
                                do i = 1, 3
                                    r_recv_buff(3*(la-1)+i) = r_atom(i, ia)
                                end do

                                !Add values to pos_send_buff
                                if(la > size(pos_send_buff)) then 
                                    allocate(pos_send_array(size(pos_send_buff)+ seg_num)) 
                                    pos_send_array(1:size(pos_send_buff)) = pos_send_buff
                                    pos_send_array(size(pos_send_buff)+1:) = 0
                                    call move_alloc(pos_send_array, pos_send_buff)
                                end if
                                pos_send_buff(la) = ia
                            end if
                        end do

                        recv_n = la
                        do ia = 1, recv_n
                            r_recv_buff(3*(ia-1)+j) = r_recv_buff(3*(ia-1)+j) + (-1)**k * box_length(j)
                        end do
                    !If not periodic in that dimension then we don't send
                    else
                        recv_l = .false.
                        recv_n = 0
                    end if

                else
                    !Get Actual send and recv coords 
                    send_coords(j) = grid_coords(j) - (-1)**k
                    recv_coords(j) = grid_coords(j) + (-1)**k

                    if(period(j).eqv..true.) then
                      call mpi_cart_rank(grid_comm, send_coords, send_rank, ierr)
                      call mpi_cart_rank(grid_comm, recv_coords, recv_rank, ierr)
                    else
                        if(k == 1) then
                            if(grid_coords(j) < num_pro(j) - 1) then 
                                call mpi_cart_rank(grid_comm, send_coords, send_rank, ierr)
                            else
                                send_l = .false.
                                send_rank = -1
                            end if

                            if(grid_coords(j) > 0) then
                                call mpi_cart_rank(grid_comm, recv_coords, recv_rank, ierr)
                            else
                              recv_l = .false.
                              recv_rank = -1
                            end if

                        else if(k == 2) then

                            if(grid_coords(j) > 0) then 
                                call mpi_cart_rank(grid_comm, send_coords, send_rank, ierr)
                            else
                                send_l = .false.
                                send_rank = -1
                            end if

                            if(grid_coords(j) < num_pro(j) - 1) then
                                call mpi_cart_rank(grid_comm, recv_coords, recv_rank, ierr)
                            else
                                recv_l = .false.
                                recv_rank = -1
                            end if
                        end if 
                    end if

                    !prepare send_buff
                    if(send_l.eqv..true.) then

                        tag_send_buff(1:send_ori) = tag_send_buff_ori(1:send_ori)
                        pos_send_buff(1:send_ori) = pos_send_buff_ori(1:send_ori)
                        type_send_buff(1:send_ori)= type_send_buff_ori(1:send_ori)  
                        r_send_buff(1:3*send_ori) = r_send_buff_ori(1:3*send_ori)

                        la = send_ori
                        do ia = atom_num_l+1, ja
                            if(dir_atom(ia) < j) then

                                la = la + 1
                                if(la > atom_num_s) then

                                    allocate( tag_send_array(atom_num_s+seg_num), &
                                              type_send_array(atom_num_s+seg_num), &  
                                              r_send_array(3*(atom_num_s+seg_num)), &
                                              stat = allostat)

                                    if(allostat /= 0) call alloc_error("Failure to allocate tag/r_send_array", allostat)

                                    tag_send_array(1:atom_num_s) = tag_send_buff(:)
                                    tag_send_array(atom_num_s+1:) = 0
                                    call move_alloc(tag_send_array, tag_send_buff)

                                    type_send_array(1:atom_num_s) = type_send_buff(:)
                                    type_send_array(atom_num_s+1:) = 0
                                    call move_alloc(type_send_array, type_send_buff)

                                    r_send_array(1:3*atom_num_s) = r_send_buff(:)
                                    r_send_array(3*atom_num_s+1:) = 0.0_wp
                                    call move_alloc(r_send_array, r_send_buff)

                                    atom_num_s = atom_num_s + seg_num

                                end if

                                tag_send_buff(la) = tag_atom(ia)
                                type_send_buff(la) = type_atom(ia)
                                do i = 1, 3
                                    r_send_buff(3*(la-1)+i) = r_atom(i, ia)
                                end do

                                if(la > size(pos_send_buff)) then 
                                    allocate(pos_send_array(size(pos_send_buff)+seg_num))
                                    pos_send_array(1:size(pos_send_buff)) = pos_send_buff
                                    pos_send_array(size(pos_send_buff)+1:) = 0
                                    call move_alloc(pos_send_array, pos_send_buff)
                                end if
                                pos_send_buff(la) = ia
                            end if
                        end do

                        send_n = la
                        !update r_send_buff for pbc
                        if(period(j).eqv..true.) then
                            do ia = 1, send_n
                                if(k == 1) then
                                    if(grid_coords(j) == num_pro(j) - 1) then
                                          r_send_buff(3*(ia-1)+j) = r_send_buff(3*(ia-1)+j) - box_length(j)
                                    end if
                                end if

                                if(k == 2) then
                                    if(grid_coords(j) == 0) then
                                        r_send_buff(3*(ia-1)+j) = r_send_buff(3*(ia-1)+j) + box_length(j)
                                    end if
                                end if
                            end do
                        end if
                    end if

                    !send/recv number
                    if(recv_l.eqv..true.) then
                        call mpi_irecv(recv_n, 1, mpi_integer, recv_rank, itag_n, grid_comm, ireq_n, ierr)
                    end if

                    if(send_l.eqv..true.) then
                        call mpi_send(send_n, 1, mpi_integer, send_rank, itag_n, grid_comm, ierr)
                    end if
                    if(recv_l.eqv..true.) then
                        call mpi_wait(ireq_n, mstatus, ierr)
                    end if

                    !update recv_buff for size
                    if(recv_l.eqv..true.) then
                        if(recv_n > atom_num_r) then

                            deallocate(tag_recv_buff, type_recv_buff,r_recv_buff,  stat = deallostat)
                            if(deallostat /= 0) call alloc_error("Failure to deallocate tag/r_recv_buff", deallostat)

                            allocate(tag_recv_buff(recv_n),  type_recv_buff(recv_n), r_recv_buff(3*recv_n), stat = allostat)
                            if(allostat /= 0) call alloc_error("Failure to allocate tag/r_recv_buff", allostat)

                            tag_recv_buff(:) = 0
                            type_recv_buff(:)= 0    
                            r_recv_buff(:) = 0.0_wp
                            atom_num_r = recv_n

                            !Reallocate accept arrays
                            deallocate(accept_send)
                            allocate(accept_send(recv_n))
                        end if
                    end if

                    !send/recv tag 
                    if(recv_l.eqv..true.) then
                        call mpi_irecv(tag_recv_buff, recv_n, mpi_integer, recv_rank, itag_t, grid_comm, ireq_t, ierr)
                    end if
                    if(send_l.eqv..true.) then
                        call mpi_send(tag_send_buff, send_n, mpi_integer, send_rank, itag_t, grid_comm, ierr)
                    end if
                    if(recv_l.eqv..true.) then
                        call mpi_wait(ireq_t, mstatus, ierr)
                    end if
                    !send/recv type
                    if(recv_l.eqv..true.) then
                        call mpi_irecv(type_recv_buff, recv_n, mpi_integer, recv_rank, itag_ty, grid_comm, ireq_ty, ierr)
                    end if
                    if(send_l.eqv..true.) then
                        call mpi_send(type_send_buff, send_n, mpi_integer, send_rank, itag_ty, grid_comm, ierr)
                    end if
                    if(recv_l.eqv..true.) then
                        call mpi_wait(ireq_t, mstatus, ierr)
                    end if

                    !send/recv r
                    if(recv_l.eqv..true.) then
                        call mpi_irecv(r_recv_buff, 3*recv_n, mpi_wp, recv_rank, itag_r, grid_comm, ireq_r, ierr)
                    end if
                    if(send_l.eqv..true.) then
                        call mpi_send(r_send_buff, 3*send_n, mpi_wp, send_rank, itag_r, grid_comm, ierr)
                    end if
                    if(recv_l.eqv..true.) then
                        call mpi_wait(ireq_r, mstatus, ierr)
                    end if
                end if

                !check if each core needs the received atoms
                if(recv_l.eqv..true.) then

                    accept_count = 0 
                    accept_send = 0 
                    do ia = 1, recv_n

                        do i = 1, 3
                            r_in(i) = r_recv_buff(3*(ia-1)+i)
                        end do

                        !We need this atom as a ghost
                        if(in_block_bd(r_in, pro_bd_out)) then

                            ja = ja + 1
                            !Resize if needed
                            if(ja > atom_num_lr) then 

                                allocate(dir_atom_array(atom_num_lr+seg_num), stat = allostat) 
                                if(allostat /= 0) call alloc_error("Failure to allocate dir_atom_array", allostat)
                                dir_atom_array(1:atom_num_lr) = dir_atom(:)
                                dir_atom_array(atom_num_lr+1:) = 0
                                call move_alloc(dir_atom_array, dir_atom)

                                call grow_at_arrays
                              
                            end if

                            tag_atom(ja) = tag_recv_buff(ia)
                            type_atom(ja)= type_recv_buff(ia)
                            r_atom(:, ja) = r_in(:)
                            dir_atom(ja) = j

                            !Now mark this atom as accepted
                            accept_count = accept_count +1
                            if(accept_count > size(accept_send)) then 
                                allocate(accept_array(size(accept_send) + seg_num))
                                accept_array(1:size(accept_send)) = accept_send
                                accept_array(size(accept_send)+1:) = 0
                                call move_alloc(accept_array, accept_send)
                            end if

                            accept_send(accept_count) = ia

                        end if
                    end do
                end if
                list_atom_logic(k, j) = send_l
                if((num_pro(j) == 1).and.(period(j))) then
                    !Resize if needed
                    if(accept_count > size(accept_recv)) then 
                        deallocate(accept_recv)
                        allocate(accept_recv(accept_count))
                    end if                    
                    accept_recv(1:accept_count) = accept_send(1:accept_count)
                    accept_count_r = accept_count
                else
                    !Now if we sent atoms this turn we have to receive the accepted ghost atoms to build our list, 
                    !First get the accepted counts
                    if(send_l) call mpi_irecv(accept_count_r, 1, mpi_integer, send_rank, 15, grid_comm, ireq_n, ierr)
                    if(recv_l) call mpi_send(accept_count, 1, mpi_integer, recv_rank, 15, grid_comm, ierr)
                    if(send_l) call mpi_wait(ireq_n, mstatus, ierr)

                    if(send_l) then 
                        !Resize if needed
                        if(accept_count_r > size(accept_recv)) then 
                            deallocate(accept_recv)
                            allocate(accept_recv(accept_count_r))
                        end if

                        !Now receive the accepted atoms
                        call mpi_irecv(accept_recv, accept_count_r, mpi_integer, send_rank, 16, grid_comm, ireq_n, ierr)
                    end if

                    if(recv_l) call mpi_send(accept_send(1:accept_count), accept_count, mpi_integer, &
                                             recv_rank, 16, grid_comm, ierr)

                    if(send_l) call mpi_wait(ireq_n, mstatus, ierr)

                end if
                    
                if(send_l.or.((num_pro(j)==1).and.period(j))) then 
                    !Now build the list of sent atoms 
                    list_atom_num(k,j) = accept_count_r

                    !initialize list_atom which contains information on which atoms were sent on every dimension
                    if(.not.allocated(list_atom)) then 
                        allocate(list_atom(accept_count_r, 2, 3), stat = allostat)
                        if (allostat>0) call alloc_error("Failure to allocate list_atom", allostat)
                    end if

                    
                    if (accept_count_r > size(list_atom(:,k,j))) then 
                        allocate(list_atom_array(accept_count_r, 2, 3))
                        list_atom_array(1:size(list_atom(:,k,j)),:,:) = list_atom(1:size(list_atom(:,k,j)),:,:)
                        list_atom_array(size(list_atom(:,k,j))+1,:,:) = 0
                        call move_alloc(list_atom_array, list_atom)
                    end if

                    do i = 1, accept_count_r
                        list_atom(i, k, j) = pos_send_buff(accept_recv(i))
                    end do 
                else
                    list_atom_num(k,j) = 0
                end if
            end do
        end do
        atom_num_lg = ja


        !debug
        if(atom_num_lr /= size(r_atom, 2)) then
            print *, 'Error: Wrong atom_num_lr', atom_num_lr, &
                     ' which should equal size(r_atom, 2)', size(r_atom, 2)
            call mpi_abort(mpi_comm_world, 1, ierr)

        else if(atom_num_lr /= size(tag_atom)) then
            print *, 'Error: Wrong atom_num_lr', atom_num_lr, &
                     ' which should equal size(tag_atom)', size(tag_atom)
            call mpi_abort(mpi_comm_world, 1, ierr)

        else if(atom_num_lr /= size(dir_atom)) then
            print *, 'Error: Wrong atom_num_lr', atom_num_lr, &
                     ' which should equal size(dir_atom)', size(dir_atom)
            call mpi_abort(mpi_comm_world, 1, ierr)

        else if(atom_num_lg > atom_num_lr) then
            print *, 'Error: Wrong atom_num_lg', atom_num_lg, &
                     ' which is larger than', atom_num_lr
            call mpi_abort(mpi_comm_world, 1, ierr)
        end if

        !resize vel_atom
        if(need_vel) then
            atom_num_ln = size(vel_atom,2)
            allocate( vel_array(3, atom_num_lr), stat = allostat) 
            if(allostat /= 0) call alloc_error("Failure to allocate vel_array", allostat)
            vel_array(:, 1:atom_num_ln) = vel_atom(:, :)
            vel_array(:, atom_num_ln+1:) = 0.0_wp
            call move_alloc(vel_array, vel_atom)
        end if

        !list_atom
!        if(.not.allocated(list_atom)) allocate(list_atom(seg_num, 2, 3), stat = allostat)
!        list_atom(:, :, :) = 0
!        list_atom_num(:, :) = 0
!        seg_num_atom = size(list_atom, 1)
!
!        do j = 1, 3
!            do k = 1, 2
!
!                send_n = 0
!                pro_bd_temp(:) = pro_bd_kj(:, k, j)
!
!                do ia = 1, atom_num_lg
!
!                    r_in(:) = r_atom(:, ia)
!                    if(in_block_bd(r_in, pro_bd_temp)) then
!                        send_n = send_n + 1
!                        if(send_n > seg_num_atom) then
!
!                            allocate(list_atom_array(seg_num_atom+seg_num, 2, 3), stat = allostat)
!                            if(allostat /= 0) call alloc_error("Failure to allocate list_atom_array", allostat)
!                            list_atom_array(1:seg_num_atom, :, :) = list_atom(:, :, :)
!                            list_atom_array(seg_num_atom+1:, :, :) = 0
!                            call move_alloc(list_atom_array, list_atom)
!                            seg_num_atom = seg_num_atom + seg_num
!                        end if
!                        list_atom(send_n, k, j) = ia
!                    end if
!                end do
!                list_atom_num(k, j) = send_n
!            end do
!        end do
        
        return
    end subroutine ghost_at

    subroutine processor_atomap

        integer :: i, j, k, send_n, recv_n, iatomap, jatomap, latomap, send_n_max, &
                   itag_r, itag_t, itag_n, ireq_r, ireq_t, ireq_n, &
                   send_rank, recv_rank, atomap_num_r

        logical :: send_l, recv_l

        integer, dimension(3) :: send_coords, recv_coords

        integer, dimension(mpi_status_size) :: mstatus

        integer, allocatable :: &
          tag_send_buff(:), tag_recv_buff(:)

        real(kind = wp), allocatable :: &
          r_send_buff(:), r_recv_buff(:)

        send_n_max = maxval(list_atomap_num)

        allocate(tag_send_buff(send_n_max), tag_recv_buff(send_n_max), r_send_buff(3*send_n_max), &
          r_recv_buff(3*send_n_max), stat = allostat)
        if(allostat /= 0) call alloc_error("Failure to allocate send_buff arrays in processor_atomap", allostat)

        !Initialize variables
        tag_send_buff(:) = 0
        tag_recv_buff(:) = 0
        r_send_buff(:) = 0.0_wp
        r_recv_buff(:) = 0.0_wp
        itag_t = 1
        itag_r = 2
        itag_n = 3
        ireq_t = 0
        ireq_r = 0
        ireq_n = 0
        send_rank = grank
        recv_rank = grank
        atomap_num_r = send_n_max
        jatomap = atomap_num_l

        !Check the ghost code to see how this looping works
        ! k == 1: positive
        ! k == 2: negative
        do j = 1, 3
            do k = 1, 2

              send_coords(:) = grid_coords(:)
              recv_coords(:) = grid_coords(:)
              send_l = .true.
              recv_l = .true.

                if(num_pro(j) == 1) then
                    send_l = .false.

                    if(period(j).eqv..true.) then

                        recv_l = .true.
                        send_n = list_atomap_num(k, j)
                        recv_n = send_n

                        if(recv_n > atomap_num_r) then
                            deallocate(tag_recv_buff, r_recv_buff, stat = deallostat)
                            if(deallostat /= 0) call alloc_error("Failure to deallocate tag/r_recv_buff", allostat)
                            allocate(tag_recv_buff(recv_n), r_recv_buff(3*recv_n), stat = allostat)

                            tag_recv_buff(:) = 0
                            r_recv_buff(:) = 0.0_wp
                            atomap_num_r = recv_n

                        end if

                        do iatomap = 1, recv_n
                            latomap = list_atomap(iatomap, k, j)
                            tag_recv_buff(iatomap) = tag_atomap(latomap)
                            do i = 1, 3
                                r_recv_buff(3*(iatomap-1)+i) = r_atomap(i, latomap)
                            end do
                            r_recv_buff(3*(iatomap-1)+j) = r_recv_buff(3*(iatomap-1)+j) &
                                                       + (-1)**k * box_length(j)
                        end do

                    else
                        recv_l = .false.
                        recv_n = 0
                    end if
                else
                    send_coords(j) = grid_coords(j) - (-1)**k
                    recv_coords(j) = grid_coords(j) + (-1)**k

                    if(period(j).eqv..true.) then
                        call mpi_cart_rank(grid_comm, send_coords, send_rank, ierr)
                        call mpi_cart_rank(grid_comm, recv_coords, recv_rank, ierr)

                    else
                        if(k == 1) then
                            if(grid_coords(j) < num_pro(j) - 1) then 
                                call mpi_cart_rank(grid_comm, send_coords, send_rank, ierr)
                            else
                                send_l = .false.
                                send_rank = -1
                            end if

                            if(grid_coords(j) > 0) then
                                call mpi_cart_rank(grid_comm, recv_coords, recv_rank, ierr)
                            else
                                recv_l = .false.
                                recv_rank = -1
                            end if

                        else if(k == 2) then
                            if(grid_coords(j) > 0) then 
                                  call mpi_cart_rank(grid_comm, send_coords, send_rank, ierr)
                            else
                                  send_l = .false.
                                  send_rank = -1
                            end if
                            if(grid_coords(j) < num_pro(j) - 1) then
                                  call mpi_cart_rank(grid_comm, recv_coords, recv_rank, ierr)
                            else
                                  recv_l = .false.
                                  recv_rank = -1
                            end if
                        end if 
                    end if

                    !check send_l
                    if(send_l.neqv.list_atomap_logic(k, j)) then
                        print *, 'Error: send_l', send_l, ' should equal list_atomap_logic', &
                                  list_atomap_logic(k, j), ' for k', k, ' and j', j
                        call mpi_abort(mpi_comm_world, 1, ierr)
                    end if

                    !prepare send_buff
                    if(send_l.eqv..true.) then
                        send_n = list_atomap_num(k, j)
                        do iatomap = 1, send_n
                            latomap = list_atomap(iatomap, k, j)
                            tag_send_buff(iatomap) = tag_atomap(latomap)
                            do i = 1, 3
                                r_send_buff(3*(iatomap-1)+i) = r_atomap(i, latomap)
                            end do

                            !update send buff for pbc
                            if(period(j).eqv..true.) then
                                if(k == 1) then
                                    if(grid_coords(j) == num_pro(j) - 1) then
                                        r_send_buff(3*(iatomap-1)+j) = r_send_buff(3*(iatomap-1)+j) - box_length(j)
                                    end if
                                end if
                                if(k == 2) then
                                    if(grid_coords(j) == 0) then
                                        r_send_buff(3*(iatomap-1)+j) = r_send_buff(3*(iatomap-1)+j) + box_length(j)
                                    end if
                                end if
                            end if
                        end do
                    end if

                    !send/recv number
                    if(recv_l.eqv..true.) then
                        call mpi_irecv(recv_n, 1, mpi_integer, recv_rank, itag_n, grid_comm, ireq_n, ierr)
                    end if
                    if(send_l.eqv..true.) then
                        call mpi_send(send_n, 1, mpi_integer, send_rank, itag_n, grid_comm, ierr)
                    end if
                    if(recv_l.eqv..true.) then
                        call mpi_wait(ireq_n, mstatus, ierr)
                    end if

                    !update recv_buff
                    if(recv_l.eqv..true.) then

                        if(recv_n > atomap_num_r) then

                            deallocate(tag_recv_buff, r_recv_buff, stat = deallostat)
                            if(deallostat /= 0) call alloc_error("Failure to deallocate recvbuff in processor_atomap", allostat)

                            allocate(tag_recv_buff(recv_n), r_recv_buff(3*recv_n), stat = allostat)

                            if(allostat /= 0) call alloc_error("Failure to allocate recv_buff in processor_atomap", allostat)

                            tag_recv_buff(:) = 0
                            r_recv_buff(:) = 0.0_wp
                            atomap_num_r = recv_n
                        end if
                    end if

                    !send/recv tag and r
                    if(recv_l.eqv..true.) then
                        call mpi_irecv(tag_recv_buff, recv_n, mpi_integer, recv_rank, itag_t, grid_comm, ireq_t, ierr)
                    end if
                    if(send_l.eqv..true.) then
                        call mpi_send(tag_send_buff, send_n, mpi_integer, send_rank, itag_t, grid_comm, ierr)
                    end if
                    if(recv_l.eqv..true.) then
                        call mpi_wait(ireq_t, mstatus, ierr)
                    end if

                    if(recv_l.eqv..true.) then
                        call mpi_irecv(r_recv_buff, 3*recv_n, mpi_wp, recv_rank, itag_r, grid_comm, ireq_r, ierr)
                    end if
                    if(send_l.eqv..true.) then
                        call mpi_send(r_send_buff, 3*send_n, mpi_wp, send_rank, itag_r, grid_comm, ierr)
                    end if
                    if(recv_l.eqv..true.) then
                        call mpi_wait(ireq_r, mstatus, ierr)
                    end if
                end if

                !check if the tags match
                if(recv_l.eqv..true.) then
                    do iatomap = 1, recv_n

                        jatomap = jatomap + 1
                        if(tag_recv_buff(iatomap) /= tag_atomap(jatomap)) then
                            print *, 'Error: Rank', rank, ' received atomap tag', &
                                      tag_recv_buff(iatomap), ' of iatomap', iatomap, &
                                     ' from rank', recv_rank, ' which should equal tag_atomap', &
                                      tag_atomap(jatomap), ' of jatomap', jatomap, &
                                     ' when j is', j, ' and k is', k, ' in atomap'

                            call mpi_abort(mpi_comm_world, 1, ierr)
                        end if

                        do i = 1, 3
                            r_atomap(i, jatomap) = r_recv_buff(3*(iatomap-1)+i)
                        end do
                    end do
                end if
            end do
        end do

        !debug

        if(jatomap /= atomap_num_lg) then
            print *, 'Error: Wrong jatomap', jatomap, &
                     ' which should equal atomap_num_lg', atomap_num_lg
            call mpi_abort(mpi_comm_world, 1, ierr)
        end if
        return

    end subroutine processor_atomap

    subroutine processor_atomistic
 
        integer :: i, j, k, send_n, recv_n, ia, ja, la, send_n_max, &
                   itag_r, itag_t, itag_n, ireq_r, ireq_t, ireq_n, &
                   send_rank, recv_rank, atom_num_r
        logical :: send_l, recv_l
        integer, dimension(3) :: send_coords, recv_coords
        integer, dimension(mpi_status_size) :: mstatus
        integer, allocatable :: tag_send_buff(:), tag_recv_buff(:)
        real(kind = wp), allocatable :: r_send_buff(:), r_recv_buff(:)

        send_n_max = maxval(list_atom_num)
        allocate(tag_send_buff(send_n_max), tag_recv_buff(send_n_max), r_send_buff(3*send_n_max), &
                 r_recv_buff(3*send_n_max), stat = allostat) 
        if(allostat /= 0) call alloc_error("Failure to allocate tag/r_send/recv_buff", allostat)

        !Initialize variables
        tag_send_buff(:) = 0
        tag_recv_buff(:) = 0
        r_send_buff(:) = 0.0_wp
        r_recv_buff(:) = 0.0_wp
        itag_t = 1
        itag_r = 2
        itag_n = 3
        ireq_t = 0
        ireq_r = 0
        ireq_n = 0
        send_rank = grank
        recv_rank = grank
        atom_num_r = send_n_max
        ja = atom_num_l

        !Check ghost code for more comments 
        ! k == 1: positive
        ! k == 2: negative
        do j = 1, 3
            do k = 1, 2

                send_coords(:) = grid_coords(:)
                recv_coords(:) = grid_coords(:)
                send_l = .true.
                recv_l = .true.

                if(num_pro(j) == 1) then

                    send_l = .false.
                    if(period(j).eqv..true.) then

                        recv_l = .true.
                        send_n = list_atom_num(k, j)
                        recv_n = send_n

                        if(recv_n > atom_num_r) then

                            deallocate(tag_recv_buff, r_recv_buff, stat = deallostat) 
                            if(deallostat /= 0) call alloc_error("Failure to deallocate tag/r_recv_buff", deallostat)

                            allocate(tag_recv_buff(recv_n), r_recv_buff(3*recv_n), stat = allostat) 
                            if(allostat /= 0) call alloc_error("Failure to allocate tag/r_recv_buff", allostat)

                            tag_recv_buff(:) = 0
                            r_recv_buff(:) = 0.0_wp
                            atom_num_r = recv_n
                        end if

                        do ia = 1, recv_n
                            la = list_atom(ia, k, j)
                            tag_recv_buff(ia) = tag_atom(la)
                            do i = 1, 3
                                r_recv_buff(3*(ia-1)+i) = r_atom(i, la)
                            end do

                            r_recv_buff(3*(ia-1)+j) = r_recv_buff(3*(ia-1)+j) + (-1)**k * box_length(j)
                        end do
                    else
                      recv_l = .false.
                      recv_n = 0
                    end if

                else

                    send_coords(j) = grid_coords(j) - (-1)**k
                    recv_coords(j) = grid_coords(j) + (-1)**k
                    if(period(j).eqv..true.) then
                        call mpi_cart_rank(grid_comm, send_coords, send_rank, ierr)
                        call mpi_cart_rank(grid_comm, recv_coords, recv_rank, ierr)
                    else
                        if(k == 1) then
                            if(grid_coords(j) < num_pro(j) - 1) then 
                                call mpi_cart_rank(grid_comm, send_coords, send_rank, ierr)
                            else
                                send_l = .false.
                                send_rank = -1
                            end if
                            if(grid_coords(j) > 0) then
                                call mpi_cart_rank(grid_comm, recv_coords, recv_rank, ierr)
                            else
                                recv_l = .false.
                                recv_rank = -1
                            end if

                        else if(k == 2) then
                            if(grid_coords(j) > 0) then 
                                call mpi_cart_rank(grid_comm, send_coords, send_rank, ierr)
                            else
                                send_l = .false.
                                send_rank = -1
                            end if
                            if(grid_coords(j) < num_pro(j) - 1) then
                                call mpi_cart_rank(grid_comm, recv_coords, recv_rank, ierr)
                            else
                                recv_l = .false.
                                recv_rank = -1
                            end if
                        end if 
                    end if
                    !check send_l
                    if(send_l.neqv.list_atom_logic(k, j)) then
                        print *, 'Error: send_l', send_l, ' should equal list_atom_logic', &
                                  list_atom_logic(k, j), ' for k', k, ' and j', j
                        call mpi_abort(mpi_comm_world, 1, ierr)
                    end if

                    !prepare send_buff
                    if(send_l.eqv..true.) then

                        send_n = list_atom_num(k, j)

                        do ia = 1, send_n

                            la = list_atom(ia, k, j)
                            tag_send_buff(ia) = tag_atom(la)

                            do i = 1, 3
                                r_send_buff(3*(ia-1)+i) = r_atom(i, la)
                            end do

                            !update send buff for pbc
                            if(period(j).eqv..true.) then
                                if(k == 1) then
                                    if(grid_coords(j) == num_pro(j) - 1) then
                                        r_send_buff(3*(ia-1)+j) = r_send_buff(3*(ia-1)+j) - box_length(j) 
                                    end if
                                end if
                                if(k == 2) then
                                    if(grid_coords(j) == 0) then
                                        r_send_buff(3*(ia-1)+j) = r_send_buff(3*(ia-1)+j) + box_length(j)
                                    end if
                                end if
                            end if
                        end do
                    end if

                    !send/recv number
                    if(recv_l.eqv..true.) then
                        call mpi_irecv(recv_n, 1, mpi_integer, recv_rank, itag_n, grid_comm, ireq_n, ierr) 
                    end if
                    if(send_l.eqv..true.) then
                        call mpi_send(send_n, 1, mpi_integer, send_rank, itag_n, grid_comm, ierr) 
                    end if
                    if(recv_l.eqv..true.) then
                        call mpi_wait(ireq_n, mstatus, ierr)
                    end if

                    !update recv_buff
                    if(recv_l.eqv..true.) then

                        if(recv_n > atom_num_r) then

                            deallocate(tag_recv_buff, r_recv_buff, stat = deallostat)
                            if(deallostat /= 0) call alloc_error("Failure to deallocate recv_buff in processor_atomistic",allostat)

                            allocate(tag_recv_buff(recv_n), r_recv_buff(3*recv_n), stat = allostat)

                            tag_recv_buff(:) = 0
                            r_recv_buff(:) = 0.0_wp
                            atom_num_r = recv_n
                        end if
                    end if

                    !send/recv tag and r
                    if(recv_l.eqv..true.) then
                        call mpi_irecv(tag_recv_buff, recv_n, mpi_integer, recv_rank, itag_t, grid_comm, ireq_t, ierr)
                    end if
                    if(send_l.eqv..true.) then
                        call mpi_send(tag_send_buff, send_n, mpi_integer, send_rank, itag_t, grid_comm, ierr) 
                    end if
                    if(recv_l.eqv..true.) then
                        call mpi_wait(ireq_t, mstatus, ierr)
                    end if
                    itag_t = itag_t  + 1

                    if(recv_l.eqv..true.) then
                        call mpi_irecv(r_recv_buff, 3*recv_n, mpi_wp, recv_rank, itag_r, grid_comm, ireq_r, ierr) 
                    end if
                    if(send_l.eqv..true.) then
                        call mpi_send(r_send_buff, 3*send_n, mpi_wp, send_rank, itag_r, grid_comm, ierr) 
                    end if
                    if(recv_l.eqv..true.) then
                        call mpi_wait(ireq_r, mstatus, ierr)
                    end if

                    itag_r = itag_r + 1
                end if

                !check if the tags match
                if(recv_l.eqv..true.) then
                    do ia = 1, recv_n
                        ja = ja + 1
                        if(tag_recv_buff(ia) /= tag_atom(ja)) then
                            print *, 'Error: Rank', rank, ' received atom tag', &
                                      tag_recv_buff(ia), ' of ia', ia, &
                                     ' from rank', recv_rank, ' which should equal tag_atom', &
                                      tag_atom(ja), ' of ja', ja, &
                                     ' when j is', j, ' and k is', k, ' in atomistic'
                            call mpi_abort(mpi_comm_world, 1, ierr)
                        end if
                        do i = 1, 3
                            r_atom(i, ja) = r_recv_buff(3*(ia-1)+i)
                        end do
                    end do
                end if
            end do
        end do
        !debug
        if(ja /= atom_num_lg) then
            print *, 'Error: Wrong ja', ja, ' which should equal atom_num_lg', atom_num_lg
            call mpi_abort(mpi_comm_world, 1, ierr)
        end if

        return
    end subroutine processor_atomistic

    subroutine pack_ele_dump(send_etype, send_tag_ele, send_esize, mask, send_basis_num, send_basis_type, &
                        send_pb_node, send_r_nodes, send_energy, send_force, send_virial, send_vel, send_int, send_real)
        !This subroutine packs all data for one element into 2 arrays, a send_int array and a send_real array 
        !containing all information needed for elements

        integer, intent(in) :: send_etype, send_tag_ele, send_esize, send_basis_num, send_basis_type(max_basisnum),&
                               send_pb_node(3, max_basisnum, ng_max_node), mask
        real(kind=wp), intent(in) :: send_r_nodes(3,max_basisnum, ng_max_node), send_energy(max_basisnum, ng_max_node), &
                                     send_force(3, max_basisnum, ng_max_node), send_virial(3,3, max_basisnum, ng_max_node), &
                                     send_vel(3, max_basisnum, ng_max_node)
        integer, dimension(:), intent(out) :: send_int
        real(kind=wp), dimension(:), intent(out) :: send_real

        integer i, j, k, inod, ibasis
        real(kind=wp) :: r_out(3)
        
        !Initialize variables
        send_int(:) = 0  
        send_real(:) = 0.0_wp

        !Calculate send counts
        !First pack the send_int variable
        send_int(1) = send_etype
        send_int(2) = send_tag_ele
        send_int(3) = send_esize
        send_int(4) = send_basis_num
        send_int(5) = mask
        j = 5+send_basis_num
        send_int(6:j) = send_basis_type(1:send_basis_num)

        !Now pack the send_real variable
        j = 1
        do inod = 1, ng_node(send_etype)
            do ibasis = 1, send_basis_num
                send_real(j) = send_energy(ibasis,inod)
                j=j+1
                !Apply periodic boundaries if needed
                r_out = send_r_nodes(:, ibasis, inod)
                call restore_pb(r_out, send_pb_node(:,ibasis,inod))
                do i =1, 3
                    send_real(j) = r_out(i) 
                    j = j+1
                end do
                do i = 1, 3
                    send_real(j) = send_force(i, ibasis, inod) 
                    j = j+1
                end do
                do i = 1, 3
                    do k = 1,3
                        send_real(j) = send_virial(i,k,ibasis,inod)
                        j=j+1
                    end do
                end do
                
                if(need_vel) then
                    do i = 1, 3
                        send_real(j) = send_vel(i, ibasis, inod)
                        j=j+1
                    end do
                end if
            end do 
        end do
        
        return
    end subroutine pack_ele_dump

    subroutine unpack_ele_dump(recv_int, recv_real, recv_etype, recv_tag_ele, recv_esize, mask, recv_basis_num, &
                               recv_basis_type,  recv_r_nodes, recv_energy, recv_force, recv_virial, recv_vel)
        !This subroutine unpacks all data for one element into 2 arrays, a recv_int array and a recv_real array 
        !containing all information needed for elements
        integer, dimension(:), intent(in) :: recv_int
        real(kind=wp), dimension(:), intent(in) :: recv_real
        integer, intent(out) :: recv_etype, recv_tag_ele, recv_esize, mask, recv_basis_num, recv_basis_type(max_basisnum)
        real(kind=wp), intent(out) :: recv_r_nodes(3,max_basisnum, ng_max_node), recv_energy(max_basisnum, ng_max_node), &
                                     recv_force(3, max_basisnum, ng_max_node), recv_virial(3,3, max_basisnum, ng_max_node), &
                                     recv_vel(3, max_basisnum, ng_max_node)

        integer i, j, k, inod, ibasis

        !Calculate recv counts
        !First pack the recv_int variable
        recv_etype     = recv_int(1)  
        recv_tag_ele   = recv_int(2)  
        recv_esize     = recv_int(3)  
        recv_basis_num = recv_int(4)  
        mask = recv_int(5)
        j = 5+recv_basis_num
        recv_basis_type(1:recv_basis_num) = recv_int(6:j) 

        !Now pack the recv_real variable
        j = 1
        do inod = 1, ng_node(recv_etype)
            do ibasis = 1, recv_basis_num
            recv_energy(ibasis,inod) = recv_real(j) 
                j=j+1
                !Apply periodic boundaries if needed
                do i =1, 3
                    recv_r_nodes(i, ibasis, inod) = recv_real(j)  
                    j = j+1
                end do
                do i = 1, 3
                    recv_force(i, ibasis, inod) = recv_real(j) 
                    j = j+1
                end do
                do i = 1, 3
                    do k = 1,3
                         recv_virial(i,k,ibasis,inod) = recv_real(j)
                        j=j+1
                    end do
                end do
                if(need_vel) then
                    do i = 1, 3
                        recv_vel(i, ibasis, inod) = recv_real(j)
                        j=j+1
                    end do
                end if
            end do 
        end do
        
        return
    end subroutine unpack_ele_dump

    subroutine pack_atom_dump(send_tag, send_type, mask, send_r, send_energy, send_force, send_virial, send_vel, &
                              send_int, send_real)
        !Pack the atom information into the dump
        !This subroutine packs all data for one atom into 2 arrays, a send_int array and a send_real array 
        !containing all information needed for the atoms

        integer, intent(in) :: send_tag, send_type, mask
        real(kind=wp), intent(in) :: send_r(3), send_energy, send_force(3), send_virial(3,3), send_vel(3)

        integer, dimension(:), intent(out) :: send_int
        real(kind=wp), dimension(:), intent(out) :: send_real

        integer i, j, k
        
        !Initialize variables
        send_int(:) = 0  
        send_real(:) = 0.0_wp

        !First pack the send_int variable
        send_int(1) = send_tag
        send_int(2) = send_type
        send_int(3) = mask

        !Now pack the send_real variable
        j = 1
        send_real(j:j+2) = send_r(:)
        j =j+3
        send_real(j) = send_energy
        j = j+1
        send_real(j:j+2) = send_force(:)
        j = j+3
        do i = 1, 3
            do k = 1, 3
                send_real(j) = send_virial(i,k)
                j = j+1
            end do 
        end do 

        if(need_vel) send_real(j:j+2) = send_vel
        return
    end subroutine pack_atom_dump

    subroutine unpack_atom_dump(recv_int, recv_real, recv_tag, recv_type, mask, recv_r, recv_energy, recv_force, recv_virial, &
                                recv_vel)
        !Unpack the atom information into the dump
        !This subroutine unpacks all data for one atom from 2 arrays, a send_int array and a send_real array 
        integer, dimension(:), intent(in) :: recv_int
        real(kind=wp), dimension(:), intent(in) :: recv_real

        integer, intent(out) :: recv_tag, recv_type , mask
        real(kind=wp), intent(out) :: recv_r(3), recv_energy, recv_force(3), recv_virial(3,3), recv_vel(3)

        integer i, j, k
        !Initialize variables
        !First pack the recv_int variable
        recv_tag  = recv_int(1) 
        recv_type = recv_int(2) 
        mask = recv_int(3)

        !Now pack the recv_real variable
        j = 1
        recv_r(:) = recv_real(j:j+2) 
        j =j+3
        recv_energy = recv_real(j) 
        j = j+1
        recv_force = recv_real(j:j+2) 
        j = j+3
        do i = 1, 3
            do k = 1, 3
                recv_virial(i,k) = recv_real(j) 
                j = j+1
            end do 
        end do 

        if(need_vel) recv_vel = recv_real(j:j+2)
        return
    end subroutine unpack_atom_dump

    subroutine pack_ele_atomap(send_etype, send_tag_ele, send_esize, send_basis_num, send_basis_type, &
                        send_pb_node, send_r_nodes, send_tag_atomap, send_r_atomap, send_int, send_real, send_vel_nodes)

        !This subroutine packs the element and interpolated atom information into one send array
        !NOTE: For this code to interface properly with the unpack ele_atomap code, the atomaps should be arranged such that
        !All basis atoms belonging to one lattice point should be grouped together. So if we have 2 basis atoms
        !of different types, index 1 should be type 1, index 2 should be type 2, index 3 should be type 1, and 
        !index 4 should be type 4. That way we can get the types from the basis_type variable

        integer, intent(in) :: send_etype, send_tag_ele, send_esize, send_basis_num, send_basis_type(max_basisnum),&
                               send_pb_node(3, max_basisnum ,ng_max_node), send_tag_atomap(atomap_max)
        real(kind=wp), intent(in) :: send_r_nodes(3,max_basisnum, ng_max_node), send_r_atomap(3,atomap_max) 
        
        integer, dimension(:), intent(out) :: send_int
        real(kind=wp), dimension(:), intent(out) :: send_real
        real(kind=wp), intent(in), optional :: send_vel_nodes(3, max_basisnum, ng_max_node)

        integer :: ia, i, j, inod, ibasis, send_atomap_num

        logical :: pack_vel
        pack_vel=.false.
        if (present(send_vel_nodes)) pack_vel = .true.

        !Initialize variables
        send_int(:) = 0  
        send_real(:) = 0.0_wp
        send_atomap_num = send_basis_num*(send_esize+1)**3

        !Calculate send counts
        !First pack the send_int variable
        send_int(1) = send_etype
        send_int(2) = send_tag_ele
        send_int(3) = send_esize
        send_int(4) = send_basis_num
        j = 4+send_basis_num
        send_int(5:j) = send_basis_type(1:send_basis_num)
        i=1
        if(periodic) then 
            do inod = 1, ng_node(send_etype)
                do ibasis =1, send_basis_num
                    do i = 1,3
                        j=j+1
                        send_int(j) = send_pb_node(i,ibasis,inod)
                    end do
                end do
            end do
        end if
        !Now add the tag_atomap
        send_int(j+1:j+send_atomap_num) = send_tag_atomap(1:send_atomap_num)

        !Now pack the send_real variable
        j = 1
        do inod = 1,ng_node(send_etype)
            do ibasis = 1, send_basis_num
                do i =1, 3
                    send_real(j) = send_r_nodes(i, ibasis, inod) 
                    j = j+1
                    if(pack_vel) then 
                        send_real(j) = send_vel_nodes(i, ibasis,inod)
                        j = j+1
                    end if
                end do
            end do 
        end do

        !Now add r_atomap
        do ia = 1, send_atomap_num
            do i = 1, 3
                send_real(j) = send_r_atomap(i,ia)
                j=j+1
            end do 
        end do
        
        return
    end subroutine pack_ele_atomap

    subroutine unpack_ele_atomap(recv_int, recv_real, recv_etype, recv_tag_ele, recv_esize, recv_basis_num, &
                          recv_basis_type, recv_pb_node, recv_r_nodes, recv_tag_atomap, &
                          recv_type_atomap, recv_r_atomap, recv_vel_nodes)
        !This subroutine unpacks the arrays that are communicated 

        integer, dimension(4+max_basisnum+(1+3*max_basisnum)*ng_max_node + atomap_max), intent(in) :: recv_int
        real(kind=wp), dimension(2*3*max_basisnum*ng_max_node + 3*atomap_max), intent(in) :: recv_real

        integer, intent(out) :: recv_etype, recv_tag_ele, recv_esize, recv_basis_num, recv_basis_type(max_basisnum),&
                               recv_pb_node(3, max_basisnum, ng_max_node), &
                               recv_type_atomap(atomap_max), recv_tag_atomap(atomap_max)
        real(kind=wp), intent(out) :: recv_r_nodes(3,max_basisnum, ng_max_node), recv_r_atomap(3,atomap_max)  

        real(kind=wp), intent(out), optional :: recv_vel_nodes(3, max_basisnum, ng_max_node)
        

        integer :: ia, i, j, inod, ibasis, recv_atomap_num
        
        logical :: pack_vel
        pack_vel = .false.
        if (present(recv_vel_nodes)) pack_vel = .true.

        recv_basis_type(:)    = 0
        recv_pb_node(:,:,:)      = 0
        recv_r_nodes(:,:,:)   = 0.0_wp
        if(pack_vel) recv_vel_nodes(:,:,:) = 0.0_wp

        !First unpack the recv_int variable
        recv_etype = recv_int(1) 
        recv_tag_ele = recv_int(2) 
        recv_esize = recv_int(3) 
        recv_basis_num = recv_int(4) 
        j = 4+recv_basis_num
        recv_basis_type(1:recv_basis_num) = recv_int(5:j) 
        if(periodic) then 
            do inod = 1, ng_node(recv_etype)
                do ibasis =1, recv_basis_num
                    do i =1,3
                        j=j+1
                        recv_pb_node(i, ibasis, inod) = recv_int(j) 
                    end do 
                end do
            end do
        end if

        !Now unpack the interpolated atom tags
        recv_atomap_num = recv_basis_num*(recv_esize+1)**3
        recv_tag_atomap(1:recv_atomap_num) = recv_int(j+1:j+recv_atomap_num) 

        !Now unpack the recv_real variable
        j = 1
        do inod = 1,ng_node(recv_etype)
            do ibasis = 1, recv_basis_num
                do i =1, 3
                    recv_r_nodes(i, ibasis, inod) = recv_real(j) 
                    j = j+1
                    if(pack_vel) then 
                        recv_vel_nodes(i, ibasis,inod) = recv_real(j) 
                        j = j+1
                    end if
                end do
            end do 
        end do

        !Now unpack r_atomap and assign the atom types 
        do ia = 1, recv_atomap_num
            recv_type_atomap(ia) = recv_basis_type(ia-int((ia-1)/recv_basis_num)*recv_basis_num)
            do i = 1, 3 
                recv_r_atomap(i,ia) = recv_real(j) 
                j = j+1
            end do 
        end do       
        return
    end subroutine unpack_ele_atomap

    subroutine pack_atoms(num_pack, tag_buff, type_buff, r_buff, send_int, send_real, vel_buff)
        !This subroutine packs the atom information into send_int and send_real
        integer, intent(in) :: num_pack 
        integer, dimension(num_pack), intent(in) :: tag_buff, type_buff
        real(kind=wp), dimension(3,num_pack), intent(in) :: r_buff
        real(kind=wp), intent(in), optional :: vel_buff(3, num_pack)
        integer, dimension(2*num_pack), intent(out) :: send_int
        real(kind=wp), dimension(2*3*num_pack), intent(out) :: send_real

        integer :: i, j, k
        logical :: pack_vel

        pack_vel = .false.
        if (present(vel_buff)) pack_vel = .true.

        !Check to make sure vel_buff is provided if need_vel
        !First pack send_int
        send_int(1:num_pack) = tag_buff
        send_int(num_pack+1:) = type_buff

        !Now pack send_real
        k = 1 
        do j = 1, num_pack
            do i = 1,3
                send_real(k) = r_buff(i,j)
                k=k+1
                if(pack_vel) then
                    send_real(k) = vel_buff(i,j)
                    k=k+1
                end if
            end do
        end do

        return

    end subroutine pack_atoms

    subroutine unpack_atoms(num_pack, recv_int, recv_real, tag_buff, type_buff, r_buff, vel_buff)
        !This subroutine packs the atom information into recv_int and recv_real
        integer, intent(in) :: num_pack 
        integer, dimension(2*num_pack), intent(in) :: recv_int
        real(kind=wp), dimension(2*3*num_pack), intent(in) :: recv_real
        integer, dimension(num_pack), intent(out) :: tag_buff, type_buff
        real(kind=wp), dimension(3,num_pack), intent(out) :: r_buff 
        real(kind=wp), dimension(3,num_pack), intent(out), optional :: vel_buff

        integer :: i, j, k
        logical :: pack_vel

        !Figure out if we are packing the velocity
        pack_vel = .false.
        if (present(vel_buff)) pack_vel = .true.

        !First pack recv_int
        tag_buff=recv_int(1:num_pack)  
        type_buff=recv_int(num_pack+1:) 

        !Now pack recv_real
        k = 1 
        do j = 1, num_pack
            do i = 1,3
                r_buff(i,j)=recv_real(k) 
                k=k+1
                if(pack_vel) then 
                    vel_buff(i,j)=recv_real(k) 
                    k=k+1
                end if
            end do
        end do

        return

    end subroutine unpack_atoms

end module comms
