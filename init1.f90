subroutine Init1
  use Consts, only : Ne, Ip1, CheckDelete
  use Laserparameters, only : Ele0, Omega, PulseStart, PulseEnd
  implicit none

  real(8) , allocatable :: y(:)[:], px(:)[:], ti(:)[:], w(:)[:], t(:), et(:)

  integer :: nim, rsize
  real(8) :: vper, randp, randti
  real(8) :: ele, dt
  real(8) :: keldysh, Up
  integer i, count, nt

  nim = num_images()
  if (this_image()==1) call CheckDelete("data/gs1.bin")

  allocate(y(Ne)[*], px(Ne)[*], ti(Ne)[*], w(Ne)[*])

  vper = dsqrt( dabs(Ele0)/ dsqrt(2.0d0 * dabs(Ip1))) !!! the width of p

  call random_seed()
  do i=1, Ne
     keldysh = 15.0d0
     do while( keldysh > 10.0d0)
        call random_number(randti)
        ti(i) = PulseStart + randti*(PulseEnd - PulseStart)*8.0d0/8.0d0
        call GetE(ti(i), ele)
        Up = ( (ele/Omega)**2 )/4.D0
        keldysh = dsqrt( dabs(Ip1)/2.0d0/Up )
     end do

     call InitX(ele, y(i))
     call random_number(randp)

     px(i) = 5.0d0 * vper * (2.0d0 * randp - 1.0d0)

     call InitWeight(w(i), ele, dabs(px(i)))
  end do

  if (this_image()==1) then
     dt = 0.1d0
     nt = int((PulseEnd - PulseStart)/dt)
     allocate(t(nt), et(nt))
     do i=1,nt
        t(i) = i*dt
        call GetE(t(i), et(i))
     end do
     inquire(iolength = rsize) t(1)
     open(100, file = 'data/et.bin', action = 'write', access = 'direct', form = 'unformatted', recl = rsize*2*nt)
     write(100, rec = 1) t, et
     close(100)
     deallocate(t, et)

     print*, 'Ne=', Ne

     open(100, file = 'data/gs1.bin', access = 'stream')
     do i = 1, nim
        write(100) y(:)[i]
     end do
     do i = 1, nim
        write(100) px(:)[i]
     end do
     do i = 1, nim
        write(100) ti(:)[i]
     end do
     do i = 1, nim
        write(100) w(:)[i]
     end do
     close(100)

     print*,'written gs'
  end if

  deallocate(px, y, ti, w)

end subroutine Init1

subroutine InitWeight(w, ele, p)
  use Laserparameters, only : Pi
  use Consts, only : Ip1
  implicit none
  real(8), intent(in) :: p, ele
  real(8), intent(out) :: w
  real(8) :: w0, w1

  w0 = ( (((2.0d0*Ip1)**2)/dabs(ele))**(2.0d0/dsqrt(2.0d0 &
       *dabs(Ip1))-1.0d0) ) * dexp( -2.0d0 *( dsqrt(2.0d0*&
       dabs(Ip1))**3 )/3.0d0/dabs(ele) )
  w1 = dsqrt( 2.0d0*dabs(Ip1) )/dabs(ele) * dexp( -(p**2)*&
       dsqrt(2.0d0*dabs(Ip1)) / dabs(ele) )

  w = w0*w1

end subroutine InitWeight
