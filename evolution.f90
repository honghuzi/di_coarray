subroutine evolution
  use Consts, only : Neq, Ne, CheckDelete, Ip1, Ip2
  ! use statearrays
  use Laserparameters, only : PulseEnd, FreeProp

  implicit none

  ! real(8) , allocatable :: u(:,:)[:], ti(:)[:], w(:)[:], work(:)
  real(8), allocatable :: ti(:)[:], w(:)[:], act(:)[:], work(:), ti2(:)[:]
  real(8), allocatable :: x1(:)[:], y1(:)[:], px1(:)[:], py1(:)[:]
  real(8), allocatable :: x2(:)[:], y2(:)[:], px2(:)[:], py2(:)[:]
  ! real(8), allocatable :: stat(:)[:], ti2(:)[:]
  real(8) :: energy, r

  character(8) :: nowtime
  integer(4) i, j
  integer(4) nim
  integer(4) :: rsize

  nim = num_images()

  ! allocate(u(Ne, Neq)[*], ti(Ne)[*], w(Ne)[*])
  allocate(ti(Ne)[*], w(Ne)[*], act(Ne)[*])
  allocate(x1(Ne)[*], y1(Ne)[*], px1(Ne)[*], py1(Ne)[*])
  allocate(x2(Ne)[*], y2(Ne)[*], px2(Ne)[*], py2(Ne)[*])
  allocate(ti2(Ne)[*])

  allocate(work(Ne))

  ! the order, r1, r2, r3,  r4,  r5, r6, r7,  r8,  r9
  ! the order, x1, y1, px1, py1, x2, y2, px2, py2, act
  x1 = 0.0d0
  py1 = 0.0d0

  ti = 0.0d0
  w = 0.0d0
  act = 0.0d0
  sync all
  if (this_image()==1) then
     call CheckDelete('data/fs.bin')
     call CheckDelete('data/fsawt.bin')
     open(100, file = 'data/gs1.bin', access = 'stream')
     open(200, file = 'data/gs2.bin', access = 'stream')

     do i = 1, nim
        read(100) work
        y1(:)[i] = work
        read(200) work
        x2(:)[i] = work
     end do

     do i = 1, nim
        read(100) work
        px1(:)[i] = work
        read(200) work
        y2(:)[i] = work
     end do

     do i = 1, nim
        read(100) work
        ti(:)[i] = work
        read(200) work
        px2(:)[i] = work
     end do

     do i = 1, nim
        read(100) work
        w(:)[i] = work
        read(200) work
        py2(:)[i] = work
     end do

     close(100)
     close(200)
     print*,'all data loaded'
  end if
  ti2 = ti

  sync all
  do i=1, Ne
     call Propagate(ti(i), PulseEnd+FreeProp, x1(i), y1(i), px1(i), py1(i), &
          x2(i), y2(i), px2(i), py2(i), act(i), ti2(i))

     act(i) = act(i) + Ip1*ti(i) + Ip2*ti2(i)

     if(this_image()==1 .and. mod(i, 10000)==0) then
        call time(nowtime)
        print *,'we have calculated', real(i)/real(Ne)*100,'%'
        print *, 'time =', nowtime
     end if
  end do
!!!!!!!!!!!!!!!!
  sync all
!!!!!!!!!!!!! sync here to make sure that all local variables are calculated
!!!!!!!!!!!!! if we sync after the copy of local variables, part of the global
!!!!!!!!!!!!! variables might be lost.
  if (this_image()==1) then
     print*, 'calculated'
     open(100, file = 'data/fs.bin', access = 'stream')
     do i = 1, nim
        write(100) x1(:)[i]
     end do
     do i = 1, nim
        write(100) y1(:)[i]
     end do
     do i = 1, nim
        write(100) px1(:)[i]
     end do
     do i = 1, nim
        write(100) py1(:)[i]
     end do
     do i = 1, nim
        write(100) x2(:)[i]
     end do
     do i = 1, nim
        write(100) y2(:)[i]
     end do
     do i = 1, nim
        write(100) px2(:)[i]
     end do
     do i = 1, nim
        write(100) py2(:)[i]
     end do
     close(100)

     open(200, file = 'data/fsawt.bin', access = 'stream')
     do i = 1, nim
        write(200) act(:)[i]
     end do
     do i = 1, nim
        write(200) w(:)[i]
     end do
     do i = 1, nim
        write(200) ti(:)[i]
     end do
     close(200)
  end if

  deallocate(ti, w, x1, y1, px1, py1, x2, y2, px2, py2, act)
end subroutine evolution
