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
        integer :: r, k
        real(dp) :: grad_res, step_size
        complex(dp) :: alpha, beta, im
        integer :: lwsp, ideg, iexp, iflag, ns
        integer, allocatable :: ipiv(:)
        complex(dp), allocatable :: wsp(:)
        alpha = 1.0
        beta = 0.0
        step_size = -0.05
        im = COMPLEX(0, 1)
        ideg = 6
        lwsp = 4 * Ne * Ne + ideg + 1
        allocate(wsp(lwsp))
        allocate(ipiv(Ne))
        allocate(grad_omega(Ne, Ne, Nk))

        grad_res = 1e6

        do while (grad_res > Nk * 1e-3)
            call omega_oracle(S, U, w, kplusb, Nk, Nb, Ne, omega, grad_omega)

            do k = 1, Nk
                call ZGPADM(ideg, Ne, step_size, grad_omega(:, :, k), Ne, wsp, lwsp, ipiv, iexp, ns, iflag)
                grad_omega(:, :, k) = reshape(wsp(iexp:iexp+Ne*Ne-1), shape(grad_omega(:, :, k)))
                call ZGEMM('N', 'N', Ne, Ne, Ne, alpha, U(:, :, k), Ne, grad_omega(:, :, k), Ne, beta, U(:, :, k), Ne)
            enddo

            grad_res = NORM2(abs(grad_omega))
            print *, omega
        enddo

    end subroutine
end module optimizer
