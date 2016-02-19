module eos_types
  implicit none

  type :: eos_data
     double precision :: rho, e, p, t
     double precision :: drho_dp, de_dp
  end type eos_data
end module eos_types

