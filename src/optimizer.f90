module optimizer
    ! This module is for development purposes.
    ! Eventually the oracles should be exposed to 
    ! an optimizer in a higher level language.
    use oracles
    implicit none

contains
    subroutine gradient_descent(S, U, w, kplusb, Nk, Nb, Ne)
        complex(dp), intent(in) :: S(:,:,:,:)
        complex(dp), intent(inout) :: U(:,:,:)
        real(dp), intent(in) :: w(:)
        integer, intent(in) :: kplusb(:,:)
        integer, intent(in) :: Nk, Nb, Ne
        real(dp) :: omega(1)
        complex(dp), allocatable :: grad_omega(:,:,:)
        complex(dp), allocatable :: grad_work(:, :)
        complex(dp), allocatable :: Rmn(:,:,:,:)
        complex(dp), allocatable :: U_work(:)
        complex(dp), allocatable :: rho_hat(:,:)
        real(dp) :: grad_res, step_size
        complex(dp) alpha, beta, theta
        integer :: n_iter, ideg, size_u_work, Nj
        parameter ( alpha = 1, beta = 0, theta = 1)
        step_size = -0.05
        ideg = 4
        size_u_work = 4 * Ne * Ne + ideg + 1
        allocate(grad_omega(Ne, Ne, Nk))
        allocate(Rmn(Ne, Ne, Nk, Nb))
        allocate(grad_work(Ne, Ne))
        allocate(U_work(size_u_work))
        allocate(rho_hat(Ne, Nb))

        grad_res = 1e6

        n_iter = 0
        do while (grad_res > REAL(Nk, 8) * 1e-4)
            call f_oracle(S, Rmn, U, w, kplusb, rho_hat, Nk, Nb, Ne, Nj, omega)
            call grad_f_oracle(Rmn, w, rho_hat, Nk, Nb, Ne, Nj, grad_omega)
            call project(U, grad_omega, grad_work, Nk, Ne)
            grad_res = NORM2(abs(grad_omega))
            call retract(U, grad_omega, U_work, Nk, Ne, ideg, size_u_work, step_size)
            print *, n_iter, omega
            n_iter = n_iter + 1
        enddo

        deallocate(grad_omega)
        deallocate(Rmn)
        deallocate(grad_work)

    end subroutine
end module optimizer
