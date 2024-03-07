module oracles
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none


contains
    subroutine omega_oracle(S, U, w, kplusb, Nk, Nb, Ne, omega, grad_omega)
        complex(dp), intent(in) :: S(:,:,:,:)
        complex(dp), intent(in) :: U(:,:,:)
        real(dp), intent(in) :: w(:)
        integer, intent(in) :: kplusb(:,:)
        integer, intent(in) :: Nk, Nb, Ne
        real(dp), intent(out) :: omega
        complex(dp) alpha, beta, theta
        parameter ( alpha = 1, beta = 0, theta = 1)
        complex(dp), allocatable :: M_work(:, :)
        complex(dp), allocatable :: Mmn(:, :, :, :)
        complex(dp), allocatable :: M_work_2(:, :)
        complex(dp), allocatable :: S_work(:, :)
        complex(dp), allocatable :: U_k(:, :)
        complex(dp), allocatable :: U_kb(:, :)

        complex(dp), allocatable :: rho_hat(:, :)
        complex(dp), allocatable :: rho_hat_conj(:, :)
        complex(dp), intent(inout) :: grad_omega(:,:,:)
        integer k, b, n, p, q

        allocate(Mmn(Ne, Ne, Nk, Nb))
        allocate(M_work(Ne, Ne))
        allocate(M_work_2(Ne, Ne))
        allocate(rho_hat(Ne, Nb))
        allocate(rho_hat_conj(Ne, Nb))
        allocate(S_work(Ne, Ne))
        allocate(U_k(Ne, Ne))
        allocate(U_kb(Ne, Ne))

        ! Compute the objective function and the gradient in one go.
        ! This is the main code to optimize and parallelize. 
        ! The two ZGEMM should go on GPUs.
        omega = 0
        Mmn = 0
        do b = 1,Nb
            M_work = 0

            do k = 1,Nk
                M_work_2 = 0
                S_work(:, :) = S(:, :, k, b)
                U_k(:, :) = U(:, :, k)
                U_kb(:, :) = U(:, :, kplusb(k, b))
                call ZGEMM('N', 'N', Ne, Ne, Ne, alpha, S_work, Ne, U_kb, Ne, beta, M_work_2, Ne)
                ! call ZGEMM('C', 'N', Ne, Ne, Ne, alpha, U_k, Ne, M_work_2, Ne, theta, M_work, Ne)
                call ZGEMM('C', 'N', Ne, Ne, Ne, alpha, U_k, Ne, M_work_2, Ne, theta, Mmn(:, :, k, b), Ne)
                M_work(:, :) = M_work(:, :) + Mmn(:, :, k, b)
            enddo

            do n = 1, Ne
                rho_hat(n, b) = M_work(n, n) / Nk
                rho_hat_conj(n, b) = conjg(rho_hat(n, b)) / abs(rho_hat(n, b))
                omega = omega + 2 * w(b) * (1 - abs(rho_hat(n, b)))
            enddo

            do k = 1,Nk
                do q = 1, Ne
                    grad_omega(:, q, k) = grad_omega(:, q, k) + w(b) * Mmn(:, q, k, b) * rho_hat_conj(q, b)
                enddo
            enddo
        enddo

        do k = 1, Nk
            grad_omega(:, :, k) = grad_omega(:, :, k) - CONJG(TRANSPOSE(grad_omega(:, :, k)))
        enddo

        grad_omega = (-2.0 / Nk) * grad_omega

        ! print *, omega
        ! print *, grad_omega(:, :, 1)

    end subroutine 

end module oracles

            ! do n = 1,Ne
            !     do k = 1,Nk
            !         do p = 1,Ne
            !             do q = 1,Ne
            !                 rho_hat(n, b) = rho_hat(n, b) + CONJG(U(p,n,k)) * S(p,q,k,b) * U(q,n,kplusb(k, b)) / Nk
            !             enddo
            !         enddo
            !     enddo
            ! enddo


            
        ! do b = 1, Nb
        !     M_work = 0
        !     do k = 1,Nk
        !         M_work_2 = 0
        !         S_work(:, :) = S(:, :, k, b)
        !         U_k(:, :) = U(:, :, k)
        !         U_kb(:, :) = U(:, :, kplusb(k, b))
        !         call ZGEMM('N', 'N', Ne, Ne, Ne, alpha, S_work, Ne, U_kb, Ne, beta, M_work_2, Ne)
        !         call ZGEMM('C', 'N', Ne, Ne, Ne, alpha, U_k, Ne, M_work_2, Ne, theta, M_work, Ne)
        !         do q = 1, Ne
        !             grad_omega(:, q, k) = grad_omega(:, q, k) + w(b) * M_work(:, q) * rho_hat_conj(q, b)
        !         enddo
        !     enddo
        ! enddo
