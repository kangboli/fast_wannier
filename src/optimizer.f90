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
        real(dp) :: grad_res, step_size
        complex(dp) alpha, beta, theta
        integer :: k
        parameter ( alpha = 1, beta = 0, theta = 1)
        step_size = -0.05
        allocate(grad_omega(Ne, Ne, Nk))

        grad_res = 1e6

        do while (grad_res > Nk * 1e-6)
            call omega_oracle(S, U, w, kplusb, Nk, Nb, Ne, omega, grad_omega)
            grad_omega = step_size * grad_omega
            call project(U, grad_omega, Nk, Ne)
            grad_res = NORM2(abs(grad_omega))
            call retract(U, grad_omega, Nk, Ne)
            print *, grad_res
            print *, omega
        enddo

    end subroutine
end module optimizer
