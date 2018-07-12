subroutine asymptote !!!  2 Dimension Kepler Code
  use Consts, only : Ne, soft1, soft2, CheckDelete
  implicit none

  real(8), dimension(:, :), allocatable :: u
  real(8), dimension(:), allocatable :: act, w, ti
  real(8) :: r1, r2, energy1, energy2
  integer(8) :: i, count, nim, rsize

  nim = num_images()

  if (this_image() == 1) then
     allocate(u(Ne*nim, 9))
     print*, 'nim=', nim
     allocate(act(Ne*nim), ti(Ne*nim), w(Ne*nim))

     inquire(iolength = rsize) act(1)
     call CheckDelete("data/difs.bin")
!     call CheckDelete("data/digs.bin")

     ! open(100, file = 'data/gs.bin', access = 'direct', form = 'unformatted', recl = rsize*4*Ne*nim)
     ! read(100, rec = 1) z0, px0, ti, w
     ! close(100)

     open(100, file = 'data/fs.bin', access = 'stream')
     read(100) u(:, 1), u(:, 2), u(:, 3), u(:, 4), u(:, 5), u(:, 6), u(:, 7), u(:, 8)
     close(100)
     print*,'read fs'
     open(100, file = 'data/fsawt.bin', access = 'stream')
     read(100) u(:, 9), w, ti
     close(100)
     print*,'read fsawt'

     close(100)

     print*,'read finished'

     count = 0
     do i = 1, Ne*nim
        r1 = dsqrt( u(i, 1)**2 + u(i, 2)**2 + soft1 )
        r2 = dsqrt( u(i, 5)**2 + u(i, 6)**2 + soft1 )
        energy1 = (u(i, 3)**2 + u(i, 4)**2)/2.0D0 - 2.0d0/r1 + 0.5d0/dsqrt( (u(i,1)-u(i,5))**2 + (u(i,2)-u(i,6))**2 + soft2)
        energy2 = (u(i, 7)**2 + u(i, 8)**2)/2.0D0 - 2.0d0/r2 + 0.5d0/dsqrt( (u(i,1)-u(i,5))**2 + (u(i,2)-u(i,6))**2 + soft2)
        ! print*, 'r1,r2 = ', r1, r2
        ! print*, 'energy1,energy2 = ', energy1, energy2

        if(energy1 > 0 .and. energy2 >0) then
           count = count + 1

           call kepler(u(i, 1), u(i, 2), u(i, 3), u(i, 4), energy1, r1)
           call kepler(u(i, 5), u(i, 6), u(i, 7), u(i, 8), energy2, r2)
           u(count, :) = u(i, :)
           act(count) = u(i, 9)
           w(count) = w(i)
           ti(count) = ti(i)
        end if
!!!!!! the following transform is given by Liu Mingming
        !         am = x(i)*pz(i) - z(i)*px(i)

        !         rl(1) =   pz(i)*am - z_nuclear*x(i)/r
        !         rl(2) = - px(i)*am - z_nuclear*z(i)/r

        !         px(i) = p*(p*(-am*rl(2)) - rl(1))/(z_nuclear+p**2*am**2)
        !         pz(i) = p*(p*( am*rl(1)) - rl(2))/(z_nuclear+p**2*am**2)
     end do
     print*, 'count=', count
     print*, 'p(di)=', real(count)/real(Ne*nim)

     open( 100, file = 'data/difs.bin', access = 'stream')
     write(100) u(1:count, 1), u(1:count, 2), &
          u(1:count, 3), u(1:count, 4), u(1:count, 5),&
          u(1:count, 6), u(1:count, 7), u(1:count, 8)
     close(100)

     open( 100, file = 'data/difsawt.bin', access = 'stream')
     write(100) act(1:count), w(1:count), ti(1:count)
     close(100)


     deallocate(u)
     ! open( 100, file = 'data/digs.bin', access = 'direct', form = 'unformatted', recl = 4*rsize*count )
     ! write(100, rec = 1) z0(1:count), px0(1:count), ti(1:count), w(1:count)
     ! close(100)
     deallocate(act, ti, w)
  end if
end subroutine asymptote


subroutine kepler(x, z, px, pz, energy, r)
  implicit none
  real(8), intent(inout) :: x, z, px, pz, energy, r
  real(8) :: y, py, len, p
  real(8) :: lx, ly, lz, ax, ay, az
  real(8) :: am, rl(2)

  y = 0.0d0
  py = 0.0d0

  p = sqrt(2*energy)

  lx=y*pz-z*py
  ly=z*px-x*pz
  lz=x*py-y*px

  len=DSQRT(lx**2+ly**2+lz**2)

  ax=py*lz-pz*ly-x/r
  ay=pz*lx-px*lz-y/r
  az=px*ly-py*lx-z/r

  px=p*((p*(ly*az-lz*ay)-ax)/(1+(p**2)*(len**2)))
!  py=p*((p*(lz*ax-lx*az)-ay)/(1+(p**2)*(len**2)))
  pz=p*((p*(lx*ay-ly*ax)-az)/(1+(p**2)*(len**2)))

end subroutine kepler
