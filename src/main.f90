program wannierise
    use param_parser
    use oracles
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    integer Nk, Nb, Ne

    complex(dp), allocatable :: U(:,:,:)
    complex(dp), allocatable :: S(:,:,:,:)
    real(dp), allocatable :: w(:)
    integer, allocatable :: kplusb(:, :)
    complex(dp), allocatable :: grad_omega(:,:,:)
    real(dp) :: omega
    character(len=:), allocatable :: dirname
    dirname = "../si_data"

    call load_dimensions(Nk, Nb, Ne, dirname) 

    ! allocate(U(Ne, Ne, Nk))
    allocate(grad_omega(Ne, Ne, Nk))
    grad_omega = 0

    call load_matrices(S, w, kplusb, Nk, Nb, Ne, dirname) 
    call load_gauge(U, Nk, Ne, dirname)
    ! call id_gauge(U, Nk, Ne)

    call omega_oracle(S, U, w, kplusb, Nk, Nb, Ne, omega, grad_omega)

    call unload_matrices(S, w, kplusb)
    deallocate(U)
    deallocate(grad_omega)

end program wannierise


