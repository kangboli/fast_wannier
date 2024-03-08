module oracles
    use, intrinsic :: iso_fortran_env, only : dp => real64
    use expm
    implicit none


contains

    subroutine omega_oracle(S, U, w, kplusb, Nk, Nb, Ne, omega, grad_omega)
        integer, intent(in) :: Nk, Nb, Ne
        complex(dp), intent(in) :: S(Ne,Ne,Nk,Nb)
        complex(dp), intent(in) :: U(Ne,Ne,Nk)
        real(dp), intent(in) :: w(Nb)
        integer, intent(in) :: kplusb(Nk,Nb)
        real(dp), intent(out) :: omega(1)
        complex(dp), intent(inout) :: grad_omega(Ne,Ne,Nk)
        complex(dp) one, zero, theta
        complex(dp) scalar, dot
        parameter ( one = 1, zero = 0, theta = 1)
        complex(dp), allocatable :: M_work(:)
        complex(dp), allocatable :: Rmn(:,:,:,:)
        complex(dp), allocatable :: rho_hat(:,:)
        complex(dp), allocatable :: rho_hat_conj(:, :)
        integer k, b, n, p, q

        allocate(Rmn(Ne, Ne, Nk, Nb))
        allocate(M_work(Ne))
        allocate(rho_hat(Ne, Nb))
        ! allocate(rho_hat_conj(Ne, Nb))

        ! Compute the objective function and the gradient in one go.
        ! This is the main code to optimize and parallelize. 
        ! The two ZGEMM should go on GPUs.
        omega = 0
        grad_omega = 0
        do b = 1,Nb
            M_work = 0

            do k = 1,Nk
                call ZGEMM('N', 'N', Ne, Ne, Ne, one, S(:, :, k, b), Ne, U(:, :, kplusb(k, b)), Ne, zero, Rmn(:, :, k, b), Ne)
                ! call ZGEMM('C', 'N', Ne, Ne, Ne, one, U(:, :, k), Ne, Rmn(:, :, k, b), Ne, one, M_work, Ne)
                do n = 1, Ne
                    dot = DOT_PRODUCT(U(:, n, k), Rmn(:, n, k, b))
                    M_work(n) = M_work(n) + dot
                enddo

            enddo

            do n = 1, Ne
                rho_hat(n, b) = M_work(n) / Nk
                omega(1) = omega(1) + 2 * w(b) * (1 - abs(rho_hat(n, b)))
                ! rho_hat_conj(n, b) = conjg(rho_hat(n, b)) / abs(rho_hat(n, b))
                rho_hat(n, b) = conjg(rho_hat(n, b)) / abs(rho_hat(n, b))
            enddo

            do k = 1,Nk
                do q = 1, Ne
                    scalar = rho_hat(q, b) * w(b) 
                    call ZAXPY(Ne, scalar, Rmn(:, q, k, b), 1, grad_omega(:, q, k), 1)
                    ! grad_omega(:, q, k) = grad_omega(:, q, k) + Rmn(:, q, k, b) * rho_hat_conj(q, b) * w(b) 
                enddo
            enddo
        enddo

        grad_omega = (-2.0 / Nk) * grad_omega

        deallocate(Rmn)
        deallocate(M_work)
        deallocate(rho_hat)
        ! deallocate(rho_hat_conj)

    end subroutine 

    subroutine project(U, grad_omega, Nk, Ne)
        integer, intent(in) :: Nk, Ne
        complex(dp), intent(inout) :: U(Ne, Ne, Nk)
        complex(dp), intent(inout) :: grad_omega(Ne, Ne, Nk)
        complex(dp), allocatable :: grad_work(:, :)
        integer :: k
        complex(dp) alpha, beta, theta
        parameter ( alpha = 1, beta = 0, theta = 1)
        allocate(grad_work(Ne, Ne))

        do k = 1,Nk
            call ZGEMM('C', 'N', Ne, Ne, Ne, alpha, U(:, :, k), Ne, grad_omega(:, :, k), Ne, beta, grad_work, Ne)
            grad_omega(:, :, k) = grad_work - CONJG(TRANSPOSE(grad_work))
            ! grad_omega(:, :, k) = grad_omega(:, :, k) - CONJG(TRANSPOSE(grad_omega(:, :, k)))
        enddo

        deallocate(grad_work)
    end subroutine project

    subroutine retract(U, DeltaU, Nk, Ne)
        integer :: Nk, Ne
        complex(dp), intent(inout) :: U(Ne, Ne, Nk)
        complex(dp), intent(inout) :: DeltaU(Ne, Ne, Nk)
        integer :: size_u_work, ideg, iexp, iflag, ns
        integer, allocatable :: ipiv(:)
        complex(dp), allocatable :: U_work(:)
        real(dp) :: t
        integer :: r, k
        complex(dp) one, zero
        parameter ( one = 1, zero = 0)
        t = 1.0
        ideg = 4
        ! This is actually a lot of memory.
        size_u_work = 4 * Ne * Ne + ideg + 1

        allocate(U_work(size_u_work))
        allocate(ipiv(Ne))

        do k = 1, Nk
            ! U = U exp(delta U)
            ! Compute the matrix exponential 
            call ZGPADM(ideg, Ne, t, DeltaU(:, :, k), Ne, U_work, size_u_work, ipiv, iexp, ns, iflag)
            DeltaU(:, :, k) = reshape(U_work(iexp:iexp+Ne*Ne-1), shape(DeltaU(:, :, k)))

            ! Update the gauge.
            call ZGEMM('N', 'N', Ne, Ne, Ne, one, U(:, :, k), Ne, DeltaU(:, :, k), Ne, zero, U_work(iexp:iexp+Ne*Ne-1), Ne)
            U(:, :, k) = reshape(U_work(iexp:iexp+Ne*Ne-1), shape(U(:, :, k)))

            ! call ZGEMM('N', 'N', Ne, Ne, Ne, one, U(:, :, k), Ne, DeltaU(:, :, k), Ne, zero, U_work, Ne)
        enddo

        deallocate(U_work)
        deallocate(ipiv)

    end subroutine retract
end module oracles

