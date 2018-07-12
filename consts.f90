module Consts
    implicit none
    private
    public :: Ne, Ip1, Ip2, im, Neq, Neq0, soft1, soft2 &
    & , CheckDelete

    complex, parameter :: im = dcmplx(0.0d0, 1.0d0)

    integer(8), parameter :: Ne = 5000000
    integer(4), parameter :: Neq = 8
    integer(4), parameter :: Neq0 = 4

    real(8), parameter :: Ip1 = - 0.9d0
    real(8), parameter :: Ip2 = - 2.0d0

    real(8), parameter :: soft1 = 0.75
    real(8), parameter :: soft2 = 0.01

contains
    subroutine CheckDelete(filename)
        implicit none
        character(*), intent(in) :: filename
        integer stat
        open(unit=1234, iostat=stat, file=filename, status='old')
        if (stat == 0) close(1234, status='delete')
    end subroutine CheckDelete
end module consts
