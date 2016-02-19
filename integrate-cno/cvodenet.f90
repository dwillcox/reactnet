module cvodenet

  use physical_constants
  use cvode_parameters
  use cno_nuc9

  implicit none

contains
  
  subroutine FCVFUN(T, Y, YDOT, IPAR, RPAR, IER) bind(C, name='fcvfun_')
    double precision, dimension(*), intent(in) :: Y, RPAR
    integer, dimension(*), intent(in) :: IPAR
    double precision, dimension(*), intent(out) :: YDOT
    double precision, intent(in) :: T
    integer, intent(out) :: IER

    double precision, dimension(net_cno%nrxn) :: rxn_rates
    integer :: i
    double precision :: ri
    
    IER = -1

    ! Calculate rates
    ! density is cv_pars%dens
    ! temperature is cv_pars%temp

    ri = 0.0d0
    do i = 1, net_cno%nrxn
       call net_cno%evaluate(cv_pars%temp, i, ri)
       rxn_rates(i) = ri
    end do
    
    YDOT(net_cno%ih1)  = -cv_pars%dens * ( &
         Y(net_cno%ic12) * Y(net_cno%ih1) * rxn_rates(net_cno%k_c12_pg_n13) + &
         Y(net_cno%ic13) * Y(net_cno%ih1) * rxn_rates(net_cno%k_c13_pg_n14) + &
         Y(net_cno%in14) * Y(net_cno%ih1) * rxn_rates(net_cno%k_n14_pg_o15) + &
         Y(net_cno%in15) * Y(net_cno%ih1) * rxn_rates(net_cno%k_n15_pa_c12) + &
         Y(net_cno%in13) * Y(net_cno%ih1) * rxn_rates(net_cno%k_n13_pg_o14))
    
    YDOT(net_cno%ihe4) = +cv_pars%dens * ( &
         Y(net_cno%in15) * Y(net_cno%ih1) * rxn_rates(net_cno%k_n15_pa_c12))
    
    YDOT(net_cno%ic12) = +cv_pars%dens * ( &
         -Y(net_cno%ic12) * Y(net_cno%ih1) * rxn_rates(net_cno%k_c12_pg_n13) + &
         Y(net_cno%in15) * Y(net_cno%ih1) * rxn_rates(net_cno%k_n15_pa_c12))

    YDOT(net_cno%ic13) = &
         -cv_pars%dens * Y(net_cno%ic13) * Y(net_cno%ih1) * rxn_rates(net_cno%k_c13_pg_n14) +&
         Y(net_cno%in13) * rxn_rates(net_cno%k_n13_b)

    YDOT(net_cno%in13) = &
         cv_pars%dens * Y(net_cno%ic12) * Y(net_cno%ih1) * rxn_rates(net_cno%k_c12_pg_n13) - &
         cv_pars%dens * Y(net_cno%in13) * Y(net_cno%ih1) * rxn_rates(net_cno%k_n13_pg_o14) - &
         Y(net_cno%in13) * rxn_rates(net_cno%k_n13_b)

    YDOT(net_cno%in14) = &
         cv_pars%dens * Y(net_cno%ic13) * Y(net_cno%ih1) * rxn_rates(net_cno%k_c13_pg_n14) - &
         cv_pars%dens * Y(net_cno%in14) * Y(net_cno%ih1) * rxn_rates(net_cno%k_n14_pg_o15) + &
         Y(net_cno%io14) * rxn_rates(net_cno%k_o14_b)

    YDOT(net_cno%in15) = &
         -cv_pars%dens * Y(net_cno%in15) * Y(net_cno%ih1) * rxn_rates(net_cno%k_n15_pa_c12) +&
         Y(net_cno%io15) * rxn_rates(net_cno%k_o15_b)

    YDOT(net_cno%io14) = &
         cv_pars%dens * Y(net_cno%in13) * Y(net_cno%ih1) * rxn_rates(net_cno%k_n13_pg_o14) - &
         Y(net_cno%io14) * rxn_rates(net_cno%k_o14_b)

    YDOT(net_cno%io15) = &
         cv_pars%dens * Y(net_cno%in14) * Y(net_cno%ih1) * rxn_rates(net_cno%k_n14_pg_o15) - &
         Y(net_cno%io15) * rxn_rates(net_cno%k_o15_b)

    YDOT(net_cno%ienuc) = 0.0d0
    do i = 1, net_cno%nnuc
       YDOT(net_cno%ienuc) = YDOT(net_cno%ienuc) + pc%N_Avogadro * net_cno%ebind(i) * YDOT(i)
    end do

    write(*,*) '______________________________'
    do i = 1, net_cno%nnuc
       write(*,*) 'YDOT(',i,'): ',YDOT(i)
    end do
    
    IER = 0 ! Successful
  end subroutine FCVFUN

  subroutine FCVDJAC(NEQ, T, Y, FY, DJAC, H_STEP, IPAR, RPAR, WK1, WK2, WK3, IER) bind(C, name='fcvdjac_')
    integer, intent(in) :: NEQ ! number of ODEs
    double precision, intent(in) :: T ! independent variable
    double precision, dimension(*), intent(in) :: Y, FY ! y and its derivative
    double precision, dimension(NEQ,*), intent(out) :: DJAC ! dense Jacobian
    double precision, intent(in) :: H_STEP ! current stepsize
    integer, dimension(*), intent(in) :: IPAR
    double precision, dimension(*), intent(in) :: RPAR
    double precision, dimension(NEQ), intent(in) :: WK1, WK2, WK3
    integer, intent(out) :: IER

    double precision, dimension(net_cno%nrxn) :: rxn_rates
    double precision :: ri
    integer :: i, j
    
    IER = -1
    
  
    do i = 1, NEQ
       do j = 1, NEQ
          DJAC(j, i) = 0.0d0
       end do
    end do

    if (T /= 0.0d0) then
       ri = 0.0d0
       do i = 1, net_cno%nrxn
          call net_cno%evaluate(cv_pars%temp, i, ri)
          rxn_rates(i) = ri
       end do

       ! DJAC(j, i) = d(YDOT(j))/dY(i)

       DJAC(1,1) = -cv_pars%dens * ( &
            Y(net_cno%ic12) * rxn_rates(net_cno%k_c12_pg_n13) + &
            Y(net_cno%ic13) * rxn_rates(net_cno%k_c13_pg_n14) + &
            Y(net_cno%in14) * rxn_rates(net_cno%k_n14_pg_o15) + &
            Y(net_cno%in15) * rxn_rates(net_cno%k_n15_pa_c12) + &
            Y(net_cno%in13) * rxn_rates(net_cno%k_n13_pg_o14))
       DJAC(1,3) = -cv_pars%dens * ( &
            Y(net_cno%ih1) * rxn_rates(net_cno%k_c12_pg_n13))
       DJAC(1,4) = -cv_pars%dens * ( &
            Y(net_cno%ic13) * Y(net_cno%ih1) * rxn_rates(net_cno%k_c13_pg_n14))
       DJAC(1,5) = -cv_pars%dens * ( &
            Y(net_cno%in13) * Y(net_cno%ih1) * rxn_rates(net_cno%k_n13_pg_o14))
       DJAC(1,6) = -cv_pars%dens * ( &
            Y(net_cno%in14) * Y(net_cno%ih1) * rxn_rates(net_cno%k_n14_pg_o15))
       DJAC(1,7) = -cv_pars%dens * ( &
            Y(net_cno%in15) * Y(net_cno%ih1) * rxn_rates(net_cno%k_n15_pa_c12))
       
       DJAC(2,1) = +cv_pars%dens * ( &
            Y(net_cno%in15) * rxn_rates(net_cno%k_n15_pa_c12))
       DJAC(2,7) = +cv_pars%dens * ( &
            Y(net_cno%ih1) * rxn_rates(net_cno%k_n15_pa_c12))

       DJAC(3,1) =        +cv_pars%dens * ( &
            -Y(net_cno%ic12) * rxn_rates(net_cno%k_c12_pg_n13) + &
            Y(net_cno%in15) * rxn_rates(net_cno%k_n15_pa_c12))
       DJAC(3,3) =        +cv_pars%dens * ( &
            -Y(net_cno%ih1) * rxn_rates(net_cno%k_c12_pg_n13))
       DJAC(3,7) =        +cv_pars%dens * ( &
            Y(net_cno%ih1) * rxn_rates(net_cno%k_n15_pa_c12))

       DJAC(4,1) =        -cv_pars%dens * Y(net_cno%ic13) * rxn_rates(net_cno%k_c13_pg_n14)
       DJAC(4,4) =        -cv_pars%dens * Y(net_cno%ih1) * rxn_rates(net_cno%k_c13_pg_n14)
       DJAC(4,5) =   rxn_rates(net_cno%k_n13_b)

       DJAC(5,1) =          cv_pars%dens * Y(net_cno%ic12) * rxn_rates(net_cno%k_c12_pg_n13) - &
            cv_pars%dens * Y(net_cno%in13) * rxn_rates(net_cno%k_n13_pg_o14)
       DJAC(5,3) =          cv_pars%dens * Y(net_cno%ih1) * rxn_rates(net_cno%k_c12_pg_n13)
       DJAC(5,5) =        -cv_pars%dens * Y(net_cno%ih1) * rxn_rates(net_cno%k_n13_pg_o14) &
            -rxn_rates(net_cno%k_n13_b)
       
       DJAC(6,1) =          cv_pars%dens * Y(net_cno%ic13) * rxn_rates(net_cno%k_c13_pg_n14)  &
            -cv_pars%dens * Y(net_cno%in14) * rxn_rates(net_cno%k_n14_pg_o15)
       DJAC(6,4) =          cv_pars%dens * Y(net_cno%ih1) * rxn_rates(net_cno%k_c13_pg_n14)
       DJAC(6,6) =           -cv_pars%dens * Y(net_cno%ih1) * rxn_rates(net_cno%k_n14_pg_o15)
       DJAC(6,8) =   +rxn_rates(net_cno%k_o14_b)

       DJAC(7,1) =          -cv_pars%dens * Y(net_cno%in15) * rxn_rates(net_cno%k_n15_pa_c12)
       DJAC(7,7) =         -cv_pars%dens * Y(net_cno%ih1) * rxn_rates(net_cno%k_n15_pa_c12)
       DJAC(7,9) =        rxn_rates(net_cno%k_o15_b)

       DJAC(8,1) =          cv_pars%dens * Y(net_cno%in13) * rxn_rates(net_cno%k_n13_pg_o14)
       DJAC(8,5) =        cv_pars%dens * Y(net_cno%ih1) * rxn_rates(net_cno%k_n13_pg_o14)
       DJAC(8,8) =             -rxn_rates(net_cno%k_o14_b)
       
       DJAC(9,1) =        cv_pars%dens * Y(net_cno%in14) * rxn_rates(net_cno%k_n14_pg_o15)
       DJAC(9,6) =        cv_pars%dens * Y(net_cno%ih1) * rxn_rates(net_cno%k_n14_pg_o15)
       DJAC(9,9) =             -rxn_rates(net_cno%k_o15_b)

    end if
       
    IER = 0 ! Success
  end subroutine FCVDJAC

  subroutine FCVROOTFN(T, Y, G, IPAR, RPAR, IER) bind(C, name='fcvrootfn_')
    double precision, intent(in) :: T
    double precision, dimension(*), intent(in) :: Y
    double precision, dimension(*), intent(inout) :: G
    integer, dimension(*), intent(in) :: IPAR
    double precision, dimension(*), intent(in) :: RPAR
    integer, intent(out) :: IER

    IER = -1
    G(1) = Y(net_cno%ih1)
    IER = 0
  end subroutine FCVROOTFN
  
end module cvodenet
