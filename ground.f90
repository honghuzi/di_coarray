subroutine Ground(t1, t2, x, y, px, py)
    use Consts, only : Neq0
    implicit none
    real(8), intent(in) :: t1, t2
    real(8), intent(inout) :: x, y, px, py

    !==== parameters for subroutine dlsoda:

    external   dudtg ! declared external in the calling program
    real(8), dimension(Neq0) :: uwork
    real(8)  :: t
    real(8)  :: tout
    integer :: itol
    real(8)  :: rtol
    real(8), dimension(Neq0) :: atol
    integer :: itask, istate, iopt
    real(8)  :: rwork(200) ! lrw = 200, at least 22 + Neq0 * max(16, Neq0 + 9)
    integer :: lrw
    integer :: iwork(40)  ! liw = 40,  at least  20 + Neq0
    integer :: liw
    integer :: jt
    integer :: i

    uwork = [x, y, px, py]
    t       = t1     ! the initial value of the independent variable
    tout    = t2     ! first point where output is desired (.ne. t)
    !set tolerance parameter
    itol    = 2      ! 1 or 2 according as atol (below) is a scalar or array
    rtol    = 1.d-9  ! relative tolerance parameter (scalar)
    atol = 1.d-11  ! absolute tolerance parameter (scalar or array)

    itask   = 1      ! 1 for normal computation of output values of uwork at t = tout
    istate  = 1      ! integer flag (input and output).  set istate = 1
    iopt    = 0      ! 0 to indicate no optional inputs used
    lrw     = 200   ! declared length of rwork (in user's dimension)
    liw     = 40     ! declared length of iwork (in user's dimension)
    jt      = 2      ! 2 means an internally generated (difference quotient) full jacobian (using Neq0 extra calls to f per df/dy value)

    !------------------------------------------
    call dlsoda(dudtg, Neq0, uwork, t, tout, itol, rtol, atol, itask,&
        &istate, iopt, rwork, lrw, iwork, liw, dudtg, jt)
    !------------------------------------------
    if (istate .ne. 2) then
        print*, 'the istate =', istate
        if((istate .eq. -1).or.(istate .eq. -4) ) then
            istate =1
            continue
        else
            stop
        end if
    end if
    ! ==== record the results:
    x = uwork(1)
    y = uwork(2)
    px = uwork(3)
    py = uwork(4)

end subroutine Ground

subroutine dudtg(Neq0, t, u, du)
    use Consts, only : soft1
    real(8), intent(in) :: t
    real(8), dimension(Neq0), intent(in) :: u
    real(8), dimension(Neq0), intent(out) :: du
    real(8)  :: r, a

    r = dsqrt( u(1)**2 + u(2)**2 + soft1 )
    du(1) = u(3)
    du(2) = u(4)

    du(3) = -u(1)* 2.0d0 /r**3
    du(4) = -u(2)* 2.0d0 /r**3
end subroutine dudtg
