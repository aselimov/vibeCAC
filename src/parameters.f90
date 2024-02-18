module parameters 

    implicit none
 
    !Precision defining parameters
    integer, parameter :: dp = selected_real_kind(15, 307), &   ! double real
                          qp = selected_real_kind(33, 4931), &  ! quadrupole real
                          wp = dp
    integer, parameter :: read_len = 500 !Read 500 characters for every input command
    integer :: error

    !Integers for which hold error codes
    integer, save :: allostat, deallostat, ierr

    !Chunk size for arrays
    integer, parameter :: seg_num = 1024

    !Various parameters used for tolerance checks
    integer, save :: lim_large_int = huge(1)
    real(kind = wp), parameter :: lim_zero = epsilon(1.0_wp), lim_small = epsilon(1.0), lim_large = huge(1.0)
    real(kind=dp), parameter :: zerotol=1E-8
 
    !constant in physics
    real(kind = wp), parameter :: &
        boltzmann = 8.6173324e-5_wp, &
        avogadro = 6.02214129e23_wp, &
        electron_volt = 1.602176565e-19_wp, &
        amu = 1.6605402e-27_wp, & 
        pi = 3.14159265358979323846_wp, &
        !nktv2p is a conversion that lammps uses to get virial stress to atm
        nktv2p = 1602176.5_wp, &
        !ftm2v is another conversion factor from lammps used in fire
        ftm2v = 1.0_wp/1.0364269d-4

 
    !constant in motion equation and kinetic energy calculation, conversion from mv^2 to energy
    real(kind = wp), parameter :: const_motion = 0.00010364269_wp 

    !Parameters for simulator
    character(len = 20), save :: simulator
    logical,save :: need_virial = .false., first_run

    !Cutoff radiuses
    real(kind=wp), save :: rc_off, rc_min, rc_neigh, rc_sq

    !Comm variables
    integer, parameter :: root = 0
    integer :: rank
 
end module parameters
