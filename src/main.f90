program wannierise
    use param_parser
    use oracles
    use optimizer
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    integer Nk, Nb, Ne

    complex(dp), allocatable :: U(:,:,:)
    complex(dp), allocatable :: S(:,:,:,:)
    real(dp), allocatable :: w(:)
    integer, allocatable :: kplusb(:, :)
    complex(dp), allocatable :: grad_omega(:,:,:)
    character(:), allocatable :: dirname

    integer :: num_args, ix
    character(len=128), dimension(:), allocatable :: args

    ix = 1
    num_args = command_argument_count()
    allocate(args(num_args))
    call get_command_argument(ix,args(ix))

    dirname = trim(args(ix))

    call load_dimensions(Nk, Nb, Ne, dirname) 
    allocate(grad_omega(Ne, Ne, Nk))
    grad_omega = 0
    call load_matrices(S, w, kplusb, Nk, Nb, Ne, dirname) 
    ! This uses the SCDM gauge to begin with
    ! call load_gauge(U, Nk, Ne, dirname)

    ! This uses the identity gauge 
    call id_gauge(U, Nk, Ne)

    ! call omega_oracle(S, U, w, kplusb, Nk, Nb, Ne, omega, grad_omega, .false.)
    call gradient_descent(S, U, w, kplusb, Nk, Nb, Ne)

    call unload_matrices(S, w, kplusb)
    deallocate(U)
    deallocate(grad_omega)

end program wannierise


