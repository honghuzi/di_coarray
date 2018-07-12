subroutine InitX(ele, y)
  implicit none
  real(8), intent(in) :: ele
  real(8), intent(out) :: y
  real(8) :: eta, ele2

  ele2 = ele

  if(ele2 < 0.0d0) then
     ele2 = -ele2
     call Eta_Turning_Point(ele2, eta)
     eta = -eta
  else
     call Eta_Turning_Point(ele2, eta)
  endif
  y = -eta/2.0d0

end subroutine InitX

! ######################################################################

subroutine Eta_Turning_Point(ele, eta)
  use Consts, only: Ip1, im
  implicit none
  real(8) :: ele, eta
  real(8) :: p, q
  real(8) :: delta
  real(8) :: s
  complex(8) :: c_s
  real(8) :: eta1, eta2, eta3
  complex(8) :: w, w_2

  w   = ( -1.0d0 + dsqrt(3.0d0)*im )/2.0d0
  w_2 = ( -1.0d0 - dsqrt(3.0d0)*im )/2.0d0
  p = -4.0d0/3.0d0*( (Ip1/ele)**2.0d0 ) + 2.0d0/ele
  q = 16.0d0/27.0d0*( (Ip1/ele)**3.0d0 ) - 4.0d0/3.0d0*Ip1/(ele**2.0d0) + 1.0d0/ele

  delta = (q/2.0d0)**2.0d0+(p/3.0d0)**3.0d0
  ! ==== only one solution:
  if( delta>=0.0d0 ) then
     s = dsqrt(delta)
     eta1 = (-q/2.0d0+s)**(1.0d0/3.0d0) + (-q/2.0d0-s)**(1.0d0/3.0d0) - 2.0d0/3.0d0*Ip1/ele
     eta2 = 0.0d0 ! over barrier ionization
     eta3 = 0.0d0 ! over barrier ionization
     eta  = 0.0d0

     ! ==== three solutions:

  elseif( delta<0.0d0) then
     c_s = im*dsqrt(-delta)
     eta1 = dble( ((-q/2.0d0+c_s)**(1.0d0/3.0d0)) + ((-q/2.0d0-c_s)**(1.0d0/3.0d0)) ) - 2.0d0/3.0d0*Ip1/ele
     eta2 = dble( w   * ((-q/2.0d0+c_s)**(1.0d0/3.0d0)) + w_2 * ((-q/2.0d0-c_s)**(1.0d0/3.0d0)) ) - 2.0d0/3.0d0*Ip1/ele
     eta3 = dble( w_2 * ((-q/2.0d0+c_s)**(1.0d0/3.0d0)) + w * ((-q/2.0d0-c_s)**(1.0d0/3.0d0)) ) - 2.0d0/3.0d0*Ip1/ele
     if(ele>0.0d0) eta = max(dabs(eta1), dabs(eta2), dabs(eta3))
     if(ele<0.0d0) eta =-max(dabs(eta1), dabs(eta2), dabs(eta3))

  endif

end subroutine Eta_Turning_Point

! ######################################################################

subroutine Xi_Turning_Point(ele, xi)
  use Consts, only: Ip1
  implicit none
  real(8) :: ele, xi
  real(8) :: p, q, r, theta
  p = -4.0d0/3.0d0*( (Ip1/ele)**2.0d0 ) - 2.0d0/ele
  q = -16.0d0/27.0d0*( (Ip1/ele)**3.0d0 ) - 4.0d0/3.0d0*Ip1/(ele**2.0d0) - 1.0d0/ele
  r = dsqrt( (-p/3.0d0)**3.0d0 )
  theta = 1.0d0/3.0d0*dacos(-q/2.0d0/r)
  xi = 2.0d0/3.0d0*Ip1/ele + 2.0d0*( r**(1.0d0/3.0d0) )*dcos(theta)

end subroutine Xi_Turning_Point

! ######################################################################

subroutine Z_Turning_Point(ele, z)
  use Consts, only: Ip1
  implicit none
  real(8) :: ele, z
  real(8) :: e1, e2
  real(8) :: delta
  if(ele>0.0d0) then
     delta = (Ip1/ele)**2.0d0 - 4.0d0/ele
     if(delta >=0.0d0) then
        e1 = ( Ip1/ele-dsqrt(delta) )/2.0d0 ! most left, choose this
        e2 = ( Ip1/ele+dsqrt(delta) )/2.0d0 !
        z  = e1
     elseif(delta < 0.0d0) then
        z  = 0.0d0
     endif
  elseif(ele<0.0d0) then
     delta = (Ip1/ele)**2.0d0 + 4.0d0/ele
     if(delta >=0.0d0) then
        e1 = ( Ip1/ele-dsqrt(delta) )/2.0d0 !
        e2 = ( Ip1/ele+dsqrt(delta) )/2.0d0 ! most right, choose this
        z  = e2
     elseif(delta < 0.0d0) then
        z  = 0.0d0
     endif
  endif

end subroutine Z_Turning_Point
