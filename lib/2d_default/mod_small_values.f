!> Module for handling problematic values in simulations, such as negative
!> pressures
module mod_small_values

  implicit none
  private

  !> How to handle small values
  character(len=20), public :: small_values_method = "error"

  !> Average over this many cells in each direction
  integer, public :: small_values_daverage = 1

  !> trace small values in the source file using traceback flag of compiler
  logical, public :: trace_small_values=.false.

  public :: small_values_error
  public :: small_values_average

contains

  subroutine small_values_error(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, w_flag, subname, smallw)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim)
    double precision, intent(in) :: smallw(1:nw)
    integer, intent(in)          :: w_flag(ixImin1:ixImax1,ixImin2:ixImax2)
    integer                      :: ix_bad(ndim), iw
    character(len=*), intent(in) :: subname

    ix_bad = maxloc(w_flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2)) + [ ixOmin1-1,&
       ixOmin2-1 ]

    if (.not.crash) then
      write(*,*) "Error: small value of ",&
          trim(prim_wnames(maxval(w_flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2)))),&
          " encountered when call ", subname
      write(*,*) "Iteration: ", it, " Time: ", global_time
      write(*,*) "Location: ", x(ix_bad(1),ix_bad(2), :)
      write(*,*) "Cell number: ", ix_bad(:)
      do iw = 1, nw
         write(*, '(A20,A,E14.7)') trim(cons_wnames(iw)), ": ", w(ix_bad(1),&
            ix_bad(2), iw)
      end do
      write(*,*) 'small value:',smallw(maxval(w_flag(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)))
      ! use erroneous arithmetic operation to crash the run
      if(trace_small_values) write(*,*) sqrt(smallw(1)-bigdouble)
      write(*,*) "Saving status at the previous time step"
      crash=.true.
    end if
  end subroutine small_values_error

  subroutine small_values_average(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, w, x, w_flag)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    integer, intent(in)             :: w_flag(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndim)
    integer                         :: iw, kxOmin1,kxOmin2,kxOmax1,kxOmax2,&
        ix1,ix2, i

    do ix2= ixOmin2,ixOmax2
    do ix1= ixOmin1,ixOmax1

    ! point with local failure identified by w_flag
    if (w_flag(ix1,ix2) /= 0) then
      ! verify in cube with border width small_values_daverage the presence of
      ! cells where all went ok
      do i = 1, max(small_values_daverage, 1)
        kxOmin1= max(ix1-i, ixOmin1);
        kxOmax1= min(ix1+i, ixOmax1);
        kxOmin2= max(ix2-i, ixOmin2);
        kxOmax2= min(ix2+i, ixOmax2);

        ! in case cells are fine within smaller cube than
        ! the userset small_values_daverage: use that smaller cube
        if (any(w_flag(kxOmin1:kxOmax1,kxOmin2:kxOmax2) == 0)) exit
      end do

      if (any(w_flag(kxOmin1:kxOmax1,kxOmin2:kxOmax2) == 0)) then
        ! within surrounding cube, cells without problem were found

        ! faulty cells are corrected by averaging here
        ! only average those which were ok and replace faulty cells
        do iw = 1, nw
          if (small_values_fix_iw(iw)) then
            w(ix1,ix2, iw) = sum(w(kxOmin1:kxOmax1,kxOmin2:kxOmax2, iw),&
                w_flag(kxOmin1:kxOmax1,kxOmin2:kxOmax2) == 0)/ &
               count(w_flag(kxOmin1:kxOmax1,kxOmin2:kxOmax2) == 0)
          end if
        end do
      else
        write(*,*) "no cells without error were found in cube of size",&
            small_values_daverage
        write(*,*) "at location:", x(ix1,ix2, 1:ndim)
        write(*,*) "at index:", ix1,ix2
        write(*,*) "w_flag(ix^D):", w_flag(ix1,ix2)
        write(*,*) "Saving status at the previous time step"
        crash=.true.
      end if
    end if
    enddo
    enddo

  end subroutine small_values_average

end module mod_small_values
