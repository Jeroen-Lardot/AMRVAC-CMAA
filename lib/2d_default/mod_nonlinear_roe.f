!> Module containing Roe solver for scalar nonlinear equation
module mod_nonlinear_roe
  use mod_physics_roe
  use mod_nonlinear_phys

  implicit none
  private

  public :: nonlinear_roe_init

contains

  subroutine nonlinear_roe_init()
    use mod_physics_roe

    nworkroe = 1

    phys_average         => nonlinear_average
    phys_get_eigenjump   => nonlinear_get_eigenjump
    phys_rtimes          => nonlinear_rtimes
  end subroutine nonlinear_roe_init

  subroutine nonlinear_average(wL, wR, x, ixmin1,ixmin2,ixmax1,ixmax2, idim,&
      wroe, workroe)
    use mod_global_parameters
    integer, intent(in)             :: ixmin1,ixmin2,ixmax1,ixmax2, idim
    double precision, intent(in)    :: wL(ixGlo1:ixGhi1,ixGlo2:ixGhi2, nw),&
        wR(ixGlo1:ixGhi1,ixGlo2:ixGhi2, nw)
    double precision, intent(inout) :: wroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2, nw)
    double precision, intent(inout) :: workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
        nworkroe)
    double precision, intent(in)    :: x(ixGlo1:ixGhi1,ixGlo2:ixGhi2, 1:2)

    wroe(ixmin1:ixmax1,ixmin2:ixmax2, rho_)=half*(wL(ixmin1:ixmax1,&
       ixmin2:ixmax2, rho_)+wR(ixmin1:ixmax1,ixmin2:ixmax2, rho_))
  end subroutine nonlinear_average

  subroutine nonlinear_get_eigenjump(wL, wR, wC, x, ixmin1,ixmin2,ixmax1,&
     ixmax2, il, idim, smalla, a, jump, workroe)

    ! Calculate the characteristic speed a and the jump in the
    ! characteristic variable in the idim direction within ixL.
    ! For a scalar equation the characteristic and conservative variables coincide
    ! The characteristic speed is just the velocity, but it should be averaged
    ! for the cell interfaces

    use mod_global_parameters

    integer, intent(in)                          :: ixmin1,ixmin2,ixmax1,&
       ixmax2, il, idim
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2, nw)       :: wL,&
        wR, wC
    double precision, dimension(ixGlo1:ixGhi1,&
       ixGlo2:ixGhi2)           :: smalla, a, jump, v
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
        nworkroe) :: workroe
    double precision, intent(in)                 :: x(ixGlo1:ixGhi1,&
       ixGlo2:ixGhi2, 1:2)
    integer                                      :: jxmin1,jxmin2,jxmax1,&
       jxmax2, ixCmin1,ixCmin2,ixCmax1,ixCmax2

    jxmin1=ixmin1+kr(idim,1);jxmin2=ixmin2+kr(idim,2)
    jxmax1=ixmax1+kr(idim,1);jxmax2=ixmax2+kr(idim,2);
    ixCmin1=ixmin1;ixCmin2=ixmin2; ixCmax1=jxmax1;ixCmax2=jxmax2;

    ! No entropy fix
    smalla(ixmin1:ixmax1,ixmin2:ixmax2)= -one
    ! The velocity is dependent of w in the nonlinear scalar equation,
    ! and thus depends on the location
    !> TODO: check this, for advection added argument to get velocity at cell edge!!!
    call nonlinear_get_v(wL, x, ixGlo1,ixGlo2,ixGhi1,ixGhi2, ixCmin1,ixCmin2,&
       ixCmax1,ixCmax2, idim, v)

    a(ixmin1:ixmax1,ixmin2:ixmax2)=(v(jxmin1:jxmax1,&
       jxmin2:jxmax2)+v(ixmin1:ixmax1,ixmin2:ixmax2))/2

    jump(ixmin1:ixmax1,ixmin2:ixmax2)=wR(ixmin1:ixmax1,ixmin2:ixmax2,&
        rho_)-wL(ixmin1:ixmax1,ixmin2:ixmax2, rho_)

  end subroutine nonlinear_get_eigenjump

  subroutine nonlinear_rtimes(q, w, ixmin1,ixmin2,ixmax1,ixmax2, iw, il, idim,&
      rq, workroe)

    ! Multiply q by R(il, iw), where R is the right eigenvalue matrix at wC.
    ! For a scalar equation the R matrix is unity

    use mod_global_parameters
    integer, intent(in)             :: ixmin1,ixmin2,ixmax1,ixmax2, iw, il,&
        idim
    double precision, intent(in)    :: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2, nw),&
        q(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
    double precision, intent(inout) :: rq(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
    double precision, intent(inout) :: workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
        nworkroe)

    rq(ixmin1:ixmax1,ixmin2:ixmax2)=q(ixmin1:ixmax1,ixmin2:ixmax2)

  end subroutine nonlinear_rtimes

end module mod_nonlinear_roe
