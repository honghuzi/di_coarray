subroutine Propagate(t1, t2, x1, y1, px1, py1, x2, y2, px2, py2, act, ti2)
  use Consts, only : Neq, soft1, soft2
  use Laserparameters, only : Pi
  implicit none
  real(8), intent(in) :: t1, t2
  real(8), intent(inout) :: x1, y1, px1, py1, x2, y2, px2, py2, act, ti2
  integer flag
  !==== parameters for subroutine dlsoda:

  external   dudt ! declared external in the calling program
  real(8), dimension(Neq) :: uwork
  real(8)  :: t
  real(8)  :: tout
  integer  :: itol
  real(8)  :: rtol
  real(8), dimension(Neq) :: atol
  integer :: itask, istate, iopt
  real(8) :: rwork(200) ! lrw = 200, at least 22 + Neq * max(16, Neq + 9)
  integer :: lrw
  integer :: iwork(40)  ! liw = 40,  at least  20 + Neq
  integer :: liw
  integer :: jt
  integer :: state2
  real(8) :: E2

  t       = t1     ! the initial value of the independent variable
  tout    = t2     ! first point where output is desired (.ne. t)
  !set tolerance parameter
  itol    = 2      ! 1 or 2 according as atol (below) is a scalar or array
  rtol    = 1.0d-11! relative tolerance parameter (scalar)
  atol    = 1.0d-11 ! absolute tolerance parameter (scalar or array)

  itask   = 1      ! 1 for normal computation of output values of uwork at t = tout
  istate  = 1      ! integer flag (input and output).  set istate = 1
  iopt    = 0      ! 0 to indicate no optional inputs used
  lrw     = 200   ! declared length of rwork (in user's dimension)
  liw     = 40     ! declared length of iwork (in user's dimension)
  jt      = 2      ! 2 means an internally generated (difference quotient) full jacobian (using Neq extra calls to f per df/dy value)

  state2 = 0
  E2 = 0.0
  uwork = [x1, y1, px1, py1, x2, y2, px2, py2]
  !------------------------------------------
  call dlsoda(dudt, Neq, uwork, t, tout, itol, rtol, atol, itask,&
       &istate, iopt, rwork, lrw, iwork, liw, dudt, jt)
  !------------------------------------------
  if (istate .ne. 2) then
     if(this_image()==1) then
        print*, 'the istate =', istate
     end if
     if((istate .eq. -1).or.(istate .eq. -4) ) then
        istate = 1
        continue
     else
        stop
     endif
  endif

  E2 = -2.0d0/sqrt(uwork(5)**2 + uwork(6)**2+soft1) + 0.5d0*(uwork(7)**2 + uwork(8)**2)&
       + 0.5d0/dsqrt((uwork(1)-uwork(5))**2+(uwork(2)-uwork(6))**2+soft2)
  if (state2 == 0 .and. E2 < 0.0) then
     act = act + 0.5d0*(uwork(3)**2 + uwork(4)**2) &
          - 2.0d0/dsqrt(uwork(1)**2+uwork(2)**2+soft1)
  elseif (state2 == 0) then
     state2 = 1
     ti2 = t
  else
     act = act + 0.5d0*(uwork(3)**2 + uwork(4)**2 + uwork(7)**2 + uwork(8)**2) &
          - 2.0d0/dsqrt(uwork(1)**2+uwork(2)**2+soft1) &
          - 2.0d0/dsqrt(uwork(5)**2+uwork(6)**2+soft1) &
          + 1.0d0/dsqrt((uwork(1)-uwork(5))**2+(uwork(2)-uwork(6))**2+soft2)
  end if

  ! ==== record the results:
  x1 = uwork(1)
  y1 = uwork(2)
  px1 = uwork(3)
  py1 = uwork(4)
  x2 = uwork(5)
  y2 = uwork(6)
  px2 = uwork(7)
  py2 = uwork(8)

end subroutine propagate

subroutine dudt(Neq, t, u, du)
  use Consts, only : soft1, soft2, Ip1
  real(8), intent(in) :: t
  real(8), dimension(Neq), intent(in) :: u
  real(8), dimension(Neq), intent(out) :: du
  real(8)  :: ele

  call GetE(t, ele)

  du(1) = u(3)
  du(2) = u(4)
  du(5) = u(7)
  du(6) = u(8)

  du(3) = -2.0d0*u(1)/(u(1)**2+u(2)**2+soft1)**1.5 + (u(1)-u(5))/((u(1)-u(5))**2+(u(2)-u(6))**2+soft2)**1.5
  du(4) = -2.0d0*u(2)/(u(1)**2+u(2)**2+soft1)**1.5 + (u(2)-u(6))/((u(1)-u(5))**2+(u(2)-u(6))**2+soft2)**1.5 - ele

  du(7) = -2.0d0*u(5)/(u(5)**2+u(6)**2+soft1)**1.5 - (u(1)-u(5))/((u(1)-u(5))**2+(u(2)-u(6))**2+soft2)**1.5
  du(8) = -2.0d0*u(6)/(u(5)**2+u(6)**2+soft1)**1.5 - (u(2)-u(6))/((u(1)-u(5))**2+(u(2)-u(6))**2+soft2)**1.5 - ele

end subroutine dudt
