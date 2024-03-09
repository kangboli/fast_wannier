module oracles
    use, intrinsic :: iso_fortran_env, only : dp => real64
    use expm
    implicit none


contains

    subroutine omega_oracle(S, Rmn, U, w, kplusb, Nk, Nb, Ne, omega, grad_omega, obj_only)
        integer, intent(in) :: Nk, Nb, Ne
        complex(dp), intent(in) :: S(Ne,Ne,Nk,Nb)
        complex(dp), intent(inout) :: Rmn(Ne,Ne,Nk,Nb)
        complex(dp), intent(in) :: U(Ne,Ne,Nk)
        real(dp), intent(in) :: w(Nb)
        integer, intent(in) :: kplusb(Nk,Nb)
        logical, intent(in) :: obj_only
        real(dp), intent(out) :: omega(1)
        complex(dp), intent(inout) :: grad_omega(Ne,Ne,Nk)
        complex(dp) one, zero, theta
        complex(dp) scalar, dot
        parameter ( one = 1, zero = 0, theta = 1)
        complex(dp), allocatable :: M_work(:)
        complex(dp), allocatable :: rho_hat(:,:)
        complex(dp), allocatable :: rho_hat_conj(:, :)
        integer k, b, n, p, q

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
                rho_hat(n, b) = conjg(rho_hat(n, b)) / abs(rho_hat(n, b))
            enddo

            if (obj_only) then
                deallocate(M_work)
                deallocate(rho_hat)
                return 
            end if

            do k = 1,Nk
                do q = 1, Ne
                    scalar = rho_hat(q, b) * w(b) 
                    call ZAXPY(Ne, scalar, Rmn(:, q, k, b), 1, grad_omega(:, q, k), 1)
                enddo
            enddo
        enddo

        grad_omega = (-2.0 / Nk) * grad_omega
        deallocate(M_work)
        deallocate(rho_hat)
    end subroutine 

    subroutine project(U, grad_omega, grad_work, Nk, Ne)
        integer, intent(in) :: Nk, Ne
        complex(dp), intent(inout) :: U(Ne, Ne, Nk)
        complex(dp), intent(inout) :: grad_omega(Ne, Ne, Nk)
        complex(dp), intent(inout) :: grad_work(Ne, Ne)
        integer :: k
        complex(dp) alpha, beta, theta
        parameter ( alpha = 1, beta = 0, theta = 1)

        do k = 1,Nk
            call ZGEMM('C', 'N', Ne, Ne, Ne, alpha, U(:, :, k), Ne, grad_omega(:, :, k), Ne, beta, grad_work, Ne)
            grad_omega(:, :, k) = grad_work - CONJG(TRANSPOSE(grad_work))
        enddo

    end subroutine project

    subroutine retract(U, DeltaU, U_work, Nk, Ne, ideg, size_u_work)
        integer :: Nk, Ne, ideg, size_u_work
        complex(dp), intent(inout) :: U(Ne, Ne, Nk)
        complex(dp), intent(in) :: DeltaU(Ne, Ne, Nk)
        complex(dp), intent(inout) :: U_work(size_u_work)
        complex(dp), allocatable :: U_work_2(:, :)

        integer :: iexp, iflag, ns
        integer, allocatable :: ipiv(:)
        real(dp) :: t
        integer :: r, k
        complex(dp) one, zero
        parameter ( one = 1, zero = 0)
        t = 1.0
        ! This is actually a lot of memory.
        allocate(ipiv(Ne))
        allocate(U_work_2(Ne, Ne))

        do k = 1, Nk
            ! U = U exp(delta U)
            ! Compute the matrix exponential 
            call ZGPADM(ideg, Ne, t, DeltaU(:, :, k), Ne, U_work, size_u_work, ipiv, iexp, ns, iflag)
            U_work_2 = reshape(U_work(iexp:iexp+Ne*Ne-1), shape(DeltaU(:, :, k)))

            ! Update the gauge.
            call ZGEMM('N', 'N', Ne, Ne, Ne, one, U(:, :, k), Ne, U_work_2, Ne, zero, U_work(iexp:iexp+Ne*Ne-1), Ne)
            U(:, :, k) = reshape(U_work(iexp:iexp+Ne*Ne-1), shape(U(:, :, k)))

            ! call ZGEMM('N', 'N', Ne, Ne, Ne, one, U(:, :, k), Ne, DeltaU(:, :, k), Ne, zero, U_work, Ne)
        enddo

        deallocate(ipiv)
        deallocate(U_work_2)

    end subroutine retract
end module oracles

