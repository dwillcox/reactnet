module eos
  use eos_types
  use ls_wrap_eos
  use polytrope_eos

  implicit none

  abstract interface
     subroutine eos_set_pars_iface(pfile_unit)
       integer, intent(in) :: pfile_unit
     end subroutine eos_set_pars_iface
  end interface
  abstract interface
     subroutine eos_p_iface(pres_target)
       double precision, intent(in) :: pres_target
     end subroutine eos_p_iface
  end interface
  abstract interface
     subroutine eos_rho_iface(dens_target)
       double precision, intent(in) :: dens_target
     end subroutine eos_rho_iface
  end interface
  abstract interface
     subroutine eos_compute_iface()
     end subroutine eos_compute_iface
  end interface
  

  procedure(eos_set_pars_iface), pointer :: eos_set_pars
  procedure(eos_p_iface), pointer ::  eos_p
  procedure(eos_rho_iface), pointer :: eos_rho
  procedure(eos_compute_iface), pointer :: eos_compute     
  type(eos_data), pointer :: eos_vars

  ! interface
  !    subroutine eos_set_pars(pfile_unit)
  !      integer, intent(in) :: pfile_unit
  !    end subroutine eos_set_pars
  !    subroutine eos_p(pres_target)
  !      double precision, intent(in) :: pres_target
  !    end subroutine eos_p
  !    subroutine eos_rho(rho_target)
  !      double precision, intent(in) :: rho_target
  !    end subroutine eos_rho
  !    subroutine eos_compute()
  !    end subroutine eos_compute
  ! end interface

  ! eos_set_pars => null()
  ! eos_p => null()
  ! eos_rho => null()
  ! eos_compute => null()
  ! eos_vars => null()
  
contains
  
  subroutine eos_init(pfile_unit)
    integer, intent(in) :: pfile_unit
    integer :: WHICH_EOS

    namelist /eosselect/ WHICH_EOS
    rewind(unit=pfile_unit)
    read(unit=pfile_unit, nml=eosselect)
    
    if (WHICH_EOS == 1) then
       ! Polytropic EOS
       eos_vars => poly_vars
       eos_set_pars => poly_init
       eos_p => poly_p
       eos_rho => poly_rho
       eos_compute => poly_compute
    else if (WHICH_EOS == 2) then
       ! Lattimer-Swesty EOS
       eos_vars => ls_vars
       eos_set_pars => ls_init
       eos_p => ls_p
       eos_rho => ls_rho
       eos_compute => ls_compute
    end if
    call eos_set_pars(pfile_unit)
  end subroutine eos_init
  
end module eos
