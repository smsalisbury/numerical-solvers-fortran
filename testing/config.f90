module config
implicit none
private

! Sets the working precision. Edit the wp
! equality to modify the precision.
!
! sp - single precision
! dp - double precision
! ep - extended precision
! qp - quadratic precision

integer,parameter::sp=selected_real_kind(6)
integer,parameter::dp=selected_real_kind(15)
integer,parameter::ep=selected_real_kind(18)
integer,parameter::qp=selected_real_kind(30)
integer,parameter::wp=sp


public wp,disp_wp

contains
  ! This subroutine displays the working precision of the current configuration.
  subroutine disp_wp
    write(*,*)wp
  end subroutine
end module config