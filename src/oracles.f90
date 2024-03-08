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
        complex(dp) scalar
        parameter ( one = 1, zero = 0, theta = 1)
        complex(dp), allocatable :: M_work(:, :)
        complex(dp), allocatable :: Rmn(:,:,:,:)

        complex(dp), allocatable :: rho_hat(:,:)
        complex(dp), allocatable :: rho_hat_conj(:, :)
        integer k, b, n, p, q

        allocate(Rmn(Ne, Ne, Nk, Nb))
        allocate(M_work(Ne, Ne))
        allocate(rho_hat(Ne, Nb))
        allocate(rho_hat_conj(Ne, Nb))

        ! Compute the objective function and the gradient in one go.
        ! This is the main code to optimize and parallelize. 
        ! The two ZGEMM should go on GPUs.
        omega = 0
        grad_omega = 0
        do b = 1,Nb
            M_work = 0

            do k = 1,Nk
                call ZGEMM('N', 'N', Ne, Ne, Ne, one, S(:, :, k, b), Ne, U(:, :, kplusb(k, b)), Ne, zero, Rmn(:, :, k, b), Ne)
                call ZGEMM('C', 'N', Ne, Ne, Ne, one, U(:, :, k), Ne, Rmn(:, :, k, b), Ne, one, M_work, Ne)
            enddo

            do n = 1, Ne
                rho_hat(n, b) = M_work(n, n) / Nk
                rho_hat_conj(n, b) = conjg(rho_hat(n, b)) / abs(rho_hat(n, b))
                omega(1) = omega(1) + 2 * w(b) * (1 - abs(rho_hat(n, b)))
            enddo

            do k = 1,Nk
                do q = 1, Ne
                    scalar = rho_hat_conj(q, b) * w(b) 
                    call ZAXPY(Ne, scalar, Rmn(:, q, k, b), 1, grad_omega(:, q, k), 1)
                    ! grad_omega(:, q, k) = grad_omega(:, q, k) + Rmn(:, q, k, b) * rho_hat_conj(q, b) * w(b) 
                enddo
            enddo
        enddo

        grad_omega = (-2.0 / Nk) * grad_omega

        deallocate(Rmn)
        deallocate(M_work)
        deallocate(rho_hat)
        deallocate(rho_hat_conj)

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
        complex(dp), allocatable :: U_work(:, :)
        integer :: lwsp, ideg, iexp, iflag, ns
        integer, allocatable :: ipiv(:)
        complex(dp), allocatable :: wsp(:)
        real(dp) :: t
        integer :: r, k
        complex(dp) one, zero, theta
        parameter ( one = 1, zero = 0, theta = 1)
        t = 1.0
        ideg = 6
        lwsp = 4 * Ne * Ne + ideg + 1

        allocate(wsp(lwsp))
        allocate(ipiv(Ne))
        allocate(U_work(Ne, Ne))

        do k = 1, Nk
            call ZGPADM(ideg, Ne, t, DeltaU(:, :, k), Ne, wsp, lwsp, ipiv, iexp, ns, iflag)
            DeltaU(:, :, k) = reshape(wsp(iexp:iexp+Ne*Ne-1), shape(DeltaU(:, :, k)))
            call ZGEMM('N', 'N', Ne, Ne, Ne, one, U(:, :, k), Ne, DeltaU(:, :, k), Ne, zero, U_work, Ne)
            U(:, :, k) = U_work(:, :)
        enddo

        deallocate(wsp)
        deallocate(ipiv)
        deallocate(U_work)

    end subroutine retract
end module oracles

