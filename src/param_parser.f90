module param_parser
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

contains
    subroutine load_dimensions(Nk, Nb, N, dirname)
        integer, intent(out) :: Nk, Nb, N
        integer :: io
        character(len=:), intent(in), allocatable :: dirname

        ! open(newunit=io, file="../si_data/dimensions.fdat", form="unformatted")
        open(newunit=io, file=dirname // "/dimensions.fdat", form="unformatted")
        read(io) Nk, Nb, N
        close(io)

    end subroutine load_dimensions

    subroutine load_matrices(S, w, kplusb, Nk, Nb, Ne, dirname)
        ! Fortran automatically deallocates arrays with intent(out).
        ! What a terrible idea.
        complex(dp), intent(out), allocatable :: S(:,:,:,:)
        real(dp), intent(out), allocatable :: w(:)
        integer, intent(out), allocatable :: kplusb(:,:)
        integer, intent(in) :: Nk, Nb, Ne
        character(len=:), intent(in), allocatable :: dirname
        integer :: io, i
        integer :: l1, l2, l3, l4

        allocate(S(Ne, Ne, Nk, Nb))
        allocate(w(Nb))
        allocate(kplusb(Nk, Nb))
        ! filepath = trim(dirname) // "/w_list.fdat"
        open(newunit=io, file=dirname // "/w_list.fdat", form="unformatted")
        ! open(newunit=io, file="../si_data/w_list.fdat", form="unformatted")
        read(io) (w(i), i=1,Nb)
        close(io)


        open(newunit=io, file=dirname // "/mmn.fdat", form="unformatted")
        ! open(newunit=io, file="../si_data/mmn.fdat", form="unformatted")
        ! read(io) ((((S(l1, l2, l3, l4), l4=1,Nk), l3=1,Nb), l2=1,Ne), l1=1,Ne)
        read(io) ((((S(l1, l2, l3, l4), l1=1,Ne), l2=1,Ne), l3=1,Nk), l4=1,Nb)
        close(io)

        open(newunit=io, file=dirname // "/kplusb.fdat", form="unformatted")
        ! open(newunit=io, file="../si_data/kplusb.fdat", form="unformatted")
        read(io) ((kplusb(l1, l2), l1=1,Nk), l2=1,Nb)
        close(io)

    end subroutine load_matrices

    subroutine load_gauge(U, Nk, Ne, dirname)
        character(len=:), intent(in), allocatable :: dirname
        complex(dp), intent(out), allocatable :: U(:,:,:)
        integer, intent(in) :: Nk, Ne
        integer :: l1, l2, l3
        integer :: io

        allocate(U(Ne, Ne, Nk))

        open(newunit=io, file=dirname // "/amn.fdat", form="unformatted")
        read(io) (((U(l1, l2, l3), l1=1,Ne), l2=1,Ne), l3=1,Nk)
        close(io)
    end subroutine load_gauge

    subroutine id_gauge(U, Nk, Ne)
        complex(dp), intent(out), allocatable :: U(:,:,:)
        integer, intent(in) :: Nk, Ne
        integer p, q, k

        allocate(U(Ne, Ne, Nk))

        do k = 1,Nk
            do p = 1, Ne
                do q = 1, Ne
                    if (p == q) then
                        U(p, q, k) = 1
                    else
                        U(p, q, k) = 0
                    endif
                enddo
            enddo
        enddo

    end subroutine id_gauge

    subroutine unload_matrices(S, w, kplusb)
        complex(dp), allocatable :: S(:,:,:,:)
        real(dp), allocatable :: w(:)
        integer, allocatable :: kplusb(:,:)

        deallocate(S)
        deallocate(w)
        deallocate(kplusb)
    end subroutine unload_matrices

end module param_parser
