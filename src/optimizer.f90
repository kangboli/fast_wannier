module optimizer
    ! This module is for development purposes.
    ! Eventually the oracles should be exposed to 
    ! an optimizer in a higher level language.
    use oracles
    use expm
    implicit none

contains
    subroutine gradient_descent(S, U, w, kplusb, Nk, Nb, Ne)
        complex(dp), intent(in) :: S(:,:,:,:)
        complex(dp), intent(inout) :: U(:,:,:)
        real(dp), intent(in) :: w(:)
        integer, intent(in) :: kplusb(:,:)
        integer, intent(in) :: Nk, Nb, Ne
        real(dp) :: omega
        complex(dp), allocatable :: grad_omega(:,:,:)
        real(dp) :: grad_res, step_size
        step_size = -0.05
        allocate(grad_omega(Ne, Ne, Nk))

        grad_res = 1e6

        do while (grad_res > Nk * 1e-3)
            call omega_oracle(S, U, w, kplusb, Nk, Nb, Ne, omega, grad_omega)
            call project(grad_omega, Nk)
            call retract(U, grad_omega, Nk, Ne)

            grad_res = NORM2(abs(grad_omega))
            print *, omega
        enddo

    end subroutine
end module optimizer
