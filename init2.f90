subroutine Init2
    use Consts, only : Ne, Ip2, CheckDelete, soft1, soft2, Neq0 &
    & , CheckDelete
    use Laserparameters, only : Pi
    implicit none
    real(8), parameter :: R0 = 0.5d0! everage position for electron 2
    real(8), parameter :: tstart = 0.0d0, tend = 50.0d0

    real(8), allocatable :: x(:)[:], y(:)[:], px(:)[:], py(:)[:]

    real(8) :: Ek, p0
    real(8) :: rand1, rand2
    integer :: i, nim, rsize

    nim = num_images()
    if(this_image()==1) then
        call CheckDelete("data/gs2.bin")
    end if
    allocate(x(Ne)[*], y(Ne)[*], px(Ne)[*], py(Ne)[*])

    Ek = Ip2 + 2.0d0/sqrt(R0**2 + soft1)
    p0 = dsqrt(2*Ek)

    sync all
    call random_seed()

    do i=1, Ne
        call random_number(rand1)
        x(i) = R0 * dsin(2*Pi*rand1)
        y(i) = R0 * dcos(2*Pi*rand1)
        call random_number(rand2)
        px(i) = p0 * dsin(2*Pi*rand2)
        py(i) = p0 * dcos(2*Pi*rand2)

        call Ground(tstart, tend, x(i), y(i), px(i), py(i))
        if(this_image()==1 .and. mod(i, 1000) == 0) then
            print*,'i=', i
        end if
    end do

    sync all

    if (this_image()==1) then
        print*, 'R0=', R0
        print*, 'p0=', p0

        open(100, file = 'data/gs2.bin', access = 'stream')
        do i = 1, nim
           write(100) x(:)[i]
        end do
        do i = 1, nim
           write(100) y(:)[i]
        end do
        do i = 1, nim
           write(100) px(:)[i]
        end do
        do i = 1, nim
           write(100) py(:)[i]
        end do
        close(100)

        print*, 'written gs2.bin'
    end if

    deallocate(x, y, px, py)
end subroutine Init2
