module oracles
    use, intrinsic :: iso_fortran_env, only : dp => real64
    use expm
    implicit none

contains

    subroutine f_oracle(S, Rmn, U, w, kplusb, rho_hat, Nk, Nb, Ne, Nj, omega)
        ! """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        ! Description: This function computes the spread omega.
        ! It has the side effect of computing and saving Rmn & rho_hat, 
        ! which can be used to compute the gradient later to avoid recomputation.
        !
        ! S: the transition matrices.
        ! Rmn: buffer that saves S * U.
        ! w: the list of weights
        ! kplusb: the tabulated index mapping so that kplusb(i(k), i(b)) = i(k+b).
        ! rho_hat: the first modes of the density. 
        ! Nk, Nb, Ne, Nj: number of kpoints, b-vectors, electrons, and bands.
        !       Nj is for future use when implementing the disentanglement.
        ! """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        
        integer, intent(in) :: Nk, Nb, Ne, Nj
        complex(dp), intent(in) :: S(Ne,Ne,Nk,Nb)
        complex(dp), intent(inout) :: Rmn(Ne,Ne,Nk,Nb)
        complex(dp), intent(in) :: U(Ne,Ne,Nk)
        real(dp), intent(in) :: w(Nb)
        integer, intent(in) :: kplusb(Nk,Nb)
        complex(dp), intent(inout) :: rho_hat(Ne,Nb)
        real(dp), intent(inout) :: omega(1)
        complex(dp) one, zero, theta
        parameter ( one = 1, zero = 0, theta = 1)
        complex(dp), allocatable :: M_work(:)
        integer k, b, n

        allocate(M_work(Ne))
        omega = 0
        do b = 1,Nb
            M_work = 0
            do k = 1,Nk
                call ZGEMM('N', 'N', Ne, Ne, Ne, one, S(:, :, k, b), Ne, U(:, :, kplusb(k, b)), Ne, zero, Rmn(:, :, k, b), Ne)
                do n = 1, Ne
                    M_work(n) = M_work(n) + DOT_PRODUCT(U(:, n, k), Rmn(:, n, k, b))
                enddo
            enddo

            rho_hat(:, b) = M_work / Nk
            omega(1) = omega(1) + 2 * w(b) * (Ne - SUM(ABS(rho_hat(:, b))))
        enddo

        deallocate(M_work)
    end subroutine f_oracle

    subroutine grad_f_oracle(Rmn, w, rho_hat, Nk, Nb, Ne, Nj, grad_omega)
        ! """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        ! Description: Compute the gradient.
        ! Rmn: the intermediate result computed from f_oracle.
        ! w: the list of weights.
        ! rho_hat: the first mode of the density computed from f_oracle.
        ! grad_omega: the output gradient.
        ! Nk, Nb, Ne, Nj: the number of kpoints, b-vectors, electrons, and bands.
        ! """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        
        integer, intent(in) :: Nk, Nb, Ne, Nj
        complex(dp), intent(in) :: Rmn(Ne,Ne,Nk,Nb)
        real(dp), intent(in) :: w(Nb)
        complex(dp), intent(inout) :: rho_hat(Ne,Nb)
        complex(dp), intent(inout) :: grad_omega(Ne,Ne,Nk)
        integer k, b, n

        grad_omega = 0
        rho_hat = CONJG(rho_hat) / ABS(rho_hat)
        
        do b = 1,Nb
            rho_hat(:, b) = rho_hat(:, b) * w(b)
            do k = 1,Nk
                do n = 1, Ne
                    call ZAXPY(Ne, rho_hat(n, b), Rmn(:, n, k, b), 1, grad_omega(:, n, k), 1)
                enddo
            enddo
        enddo

        grad_omega = (-2.0 / REAL(Nk, 8)) * grad_omega
    end subroutine grad_f_oracle

    subroutine project(U, grad_omega, grad_work, Nk, Ne)
        integer, intent(in) :: Nk, Ne
        complex(dp), intent(inout) :: U(Ne, Ne, Nk)
        complex(dp), intent(inout) :: grad_omega(Ne, Ne, Nk)
        complex(dp), intent(inout) :: grad_work(Ne, Ne)
        integer :: k
        complex(dp) one, zero
        parameter ( one = 1, zero = 0)

        do k = 1,Nk
            call ZGEMM('C', 'N', Ne, Ne, Ne, one, U(:, :, k), Ne, grad_omega(:, :, k), Ne, zero, grad_work, Ne)
            grad_omega(:, :, k) = grad_work - CONJG(TRANSPOSE(grad_work))
        enddo

    end subroutine project

    subroutine retract(U, DeltaU, U_work, Nk, Ne, ideg, size_u_work, t)
        ! """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        ! Description: This is a fairly slow retraction. 
        ! This should be replaced with a faster retraction algorithm, so 
        ! there is no need to optimize the code for now.
        ! """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        
        integer :: Nk, Ne, ideg, size_u_work
        complex(dp), intent(inout) :: U(Ne, Ne, Nk)
        complex(dp), intent(in) :: DeltaU(Ne, Ne, Nk)
        complex(dp), intent(inout) :: U_work(size_u_work)
        real(dp), intent(in) :: t
        complex(dp), allocatable :: U_work_2(:, :)

        integer :: iexp, iflag, ns
        integer, allocatable :: ipiv(:)
        integer :: r, k
        complex(dp) one, zero
        parameter ( one = 1, zero = 0)
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

        enddo

        deallocate(ipiv)
        deallocate(U_work_2)

    end subroutine retract
end module oracles

