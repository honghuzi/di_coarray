Module Statearrays
  implicit none
  use Consts, only : Neq
  private
  public :: state2nd, ti2nd

  integer(8) :: state2nd(Neq*num_images())
  real(8) :: ti2nd(Neq*num_images())

end Module Statearrays
