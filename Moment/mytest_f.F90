module mytest_module
  use amrex_fort_module
  use amrex_constants_module
  use amrex_error_module
  use amrex_ebcellflag_module, only : is_single_valued_cell
  implicit none

contains

  subroutine set_intg_2d (lo, hi, Sx, Sx2, Sy, Sy2, &
       flag, flo, fhi, vol, vlo, vhi, ax, axlo, axhi, ay, aylo, ayhi, bcen, blo, bhi) &
       bind(c,name='set_intg_2d')
    integer, dimension(2) :: lo, hi, flo, fhi, vlo, vhi, axlo, axhi, aylo, ayhi, blo, bhi
    real(amrex_real), intent(inout) :: Sx  (  lo(1):  hi(1),  lo(2):  hi(2))
    real(amrex_real), intent(inout) :: Sx2 (  lo(1):  hi(1),  lo(2):  hi(2))
    real(amrex_real), intent(inout) :: Sy  (  lo(1):  hi(1),  lo(2):  hi(2))
    real(amrex_real), intent(inout) :: Sy2 (  lo(1):  hi(1),  lo(2):  hi(2))
    real(amrex_real), intent(in   ) :: vol ( vlo(1): vhi(1), vlo(2): vhi(2))
    real(amrex_real), intent(in   ) :: ax  (axlo(1):axhi(1),axlo(2):axhi(2))
    real(amrex_real), intent(in   ) :: ay  (aylo(1):ayhi(1),aylo(2):ayhi(2))
    real(amrex_real), intent(in   ) :: bcen( blo(1): bhi(1), blo(2): bhi(2),2)
    integer         , intent(in   ) :: flag( flo(1): fhi(1), flo(2): fhi(2))

    integer :: i, j
    real(amrex_real) :: axm, axp, aym, ayp, apnorm, apnorminv, anrmx, anrmy, bcx, bcy
    real(amrex_real) :: xmin, xmax, ymin, ymax
    real(amrex_real), parameter :: almostone = 1.d0 - 1.d2*epsilon(1._amrex_real)
    real(amrex_real), parameter :: twentyfourth = 1.d0/24.d0

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)

          if (is_single_valued_cell(flag(i,j)) .and. vol(i,j).le.almostone) then

             axm = ax(i,j)
             axp = ax(i+1,j)
             aym = ay(i,j)
             ayp = ay(i,j+1)

             apnorm = sqrt((axm-axp)**2 + (aym-ayp)**2)
             if (apnorm .eq. 0.d0) then
                print *, "xxxxx", i, j, axm, axp, aym, ayp, vol(i,j), bcen(i,j,:)
                call amrex_abort("amrex_mlndlap_set_connection: we are in trouble")
             end if

             apnorminv = 1.d0/apnorm
             anrmx = (axm-axp) * apnorminv  ! pointing to the wall
             anrmy = (aym-ayp) * apnorminv

             bcx = bcen(i,j,1)
             bcy = bcen(i,j,2)

             if (anrmx .eq. 0.d0) then
                Sx(i,j) = 0.d0
                Sx2(i,j) = twentyfourth*(axm+axp)
             else if (anrmy .eq. 0.d0) then
                Sx(i,j)  = 0.125d0     *(axp-axm) + anrmx*0.5d0*bcx**2
                Sx2(i,j) = twentyfourth*(axp+axm) + anrmx*third*bcx**3
             else
                if (anrmx .gt. 0.d0) then
                   xmin = -0.5d0 + min(aym,ayp)
                   xmax = -0.5d0 + max(aym,ayp)
                else
                   xmin = 0.5d0 - max(aym,ayp)
                   xmax = 0.5d0 - min(aym,ayp)
                end if
                Sx(i,j)  = 0.125d0     *(axp-axm) + (anrmx/abs(anrmy))*sixth  *(xmax**3-xmin**3)
                Sx2(i,j) = twentyfourth*(axp+axm) + (anrmx/abs(anrmy))*twelfth*(xmax**4-xmin**4)
             end if

             if (anrmy .eq. 0.d0) then
                Sy(i,j) = 0.d0
                Sy2(i,j) = twentyfourth*(aym+ayp)
             else if (anrmx .eq. 0.d0) then
                Sy(i,j)  = 0.125d0     *(ayp-aym) + anrmy*0.5d0*bcy**2
                Sy2(i,j) = twentyfourth*(ayp+aym) + anrmy*third*bcy**3
             else
                if (anrmy .gt. 0.d0) then
                   ymin = -0.5d0 + min(axm,axp)
                   ymax = -0.5d0 + max(axm,axp)
                else
                   ymin = 0.5d0 - max(axm,axp)
                   ymax = 0.5d0 - min(axm,axp)
                end if
                Sy(i,j)  = 0.125d0     *(ayp-aym) + (anrmy/abs(anrmx))*sixth  *(ymax**3-ymin**3)
                Sy2(i,j) = twentyfourth*(ayp+aym) + (anrmy/abs(anrmx))*twelfth*(ymax**4-ymin**4)
             end if

          end if
       end do
    end do

  end subroutine set_intg_2d

end module mytest_module
