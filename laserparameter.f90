module Laserparameters
    implicit none
    private
    public :: Pi, Omega, Ele0, A0, LaserCycle, N_cycle, &
         PulseStart, PulseEnd, RampOnEnd, RampOffStart, &
         PulseCenter, Cep, FreeProp, Duration, RampOnTime,&
         RampOffTime

    real(8), parameter :: Pi = 3.141592653589
    real(8), parameter :: LightSpeed = 3.0d8
    real(8), parameter :: I_si = 1.0d15
    real(8), parameter :: Wavelength = 390.0d-9
    real(8), parameter :: Omega_si = 2 * Pi * LightSpeed / Wavelength
    real(8), parameter :: Ele_si = 2.742d3 * dsqrt(I_si)
! ==================================================
!  atomic units
    real(8), parameter :: Omega = Omega_si * 2.4189d-17
    real(8), parameter :: Ele0 = Ele_si / 5.1421d11
    real(8), parameter :: A0 = Ele0 / Omega

    real(8), parameter :: LaserCycle = 2 * Pi / Omega
    real(8), parameter :: N_cycle = 8.0d0
    real(8), parameter :: Duration = N_cycle * Lasercycle
    real(8), parameter :: PulseStart = 0.0d0
    real(8), parameter :: PulseEnd = Duration + PulseStart

! pulse shape
    real(8), parameter :: RampOnTime = LaserCycle * 2.0
    real(8), parameter :: RampOffTime = LaserCycle * 2.0

    real(8), parameter :: RampOnEnd = PulseStart + RampOnTime
    real(8), parameter :: RampOffStart = PulseEnd - RampOffTime

    real(8), parameter :: PulseCenter = 0.0d0
    real(8), parameter :: Cep = 0.0d0
    real(8), parameter :: FreeProp = 0.0d0 !3.0d0 * LaserCycle

end module

subroutine GetE(t, ele)
  use Laserparameters, only : PulseStart, PulseEnd, RampOnEnd, &
       RampOffStart, Omega, Cep, PulseCenter, Ele0, Duration, &
       Pi, RampOnTime, RampOffTime
    implicit none
    real(8), intent(in) :: t
    real(8), intent(out) :: ele
    real(8) :: envelope

    if(t<PulseStart .or. t>=PulseEnd) then
        envelope = 0.0d0
        else if(t>=PulseStart .and. t<RampOnEnd) then
            !envelope = (t - PulseStart)/(RampOnEnd - PulseStart)
           envelope = (dcos( Pi * (t - PulseCenter + 0.5d0 * (Duration &
                & -RampOffTime-RampOnTime)) / (2*RampOnTime ) ))**2

        elseif(t>=RampOnEnd .and. t<RampOffStart) then
            envelope = 1.0d0

        else
!            envelope = (PulseEnd - t)/(PulseEnd - RampOffStart)
           envelope = (dcos( Pi * (t - PulseCenter - 0.5d0 * (Duration &
                & -RampOffTime-RampOnTime)) / (2*RampOffTime ) ))**2

    end if
    ele = Ele0 * envelope * dcos(Omega*(t-PulseCenter) + Cep)
end subroutine GetE
