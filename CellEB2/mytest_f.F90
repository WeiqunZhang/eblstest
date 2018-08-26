
module mytest
  use amrex_fort_module, only : amrex_real
  use amrex_constants_module, only : half, three, zero, seven
  use amrex_ebcellflag_module, only : is_regular_cell, is_covered_cell, is_single_valued_cell
  implicit none

contains

#if (AMREX_SPACEDIM == 2)
  subroutine mytest_set_phi_reg(lo, hi, phie, elo, ehi, rhs, rlo, rhi, dx) &
       bind(c,name='mytest_set_phi_reg')
    integer, dimension(2), intent(in) :: lo, hi, elo, ehi, rlo, rhi
    real(amrex_real), intent(in) :: dx(2)
    real(amrex_real), intent(inout) ::  phie(elo(1):ehi(1),elo(2):ehi(2))
    real(amrex_real), intent(inout) ::  rhs (rlo(1):rhi(1),rlo(2):rhi(2))

    integer :: i,j
    real(amrex_real) :: x, y, r2, theta

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          x = (i+half)*dx(1) - half
          y = (j+half)*dx(2) - half
          theta = atan2(x,y)
          r2 = x**2 + y**2
          phie(i,j) = r2**2 * cos(three*theta)
       end do
    end do
  end subroutine mytest_set_phi_reg

  subroutine mytest_set_phi_eb(lo, hi, phie, elo, ehi, phib, blo, bhi, &
       rhs, rlo, rhi, flag, flo, fhi, bcent, clo, chi, dx) &
       bind(c,name='mytest_set_phi_eb')
    integer, dimension(2), intent(in) :: lo, hi, elo, ehi, blo, bhi, &
         rlo, rhi, flo, fhi, clo, chi
    real(amrex_real), intent(in) :: dx(2)
    real(amrex_real), intent(inout) ::  phie(elo(1):ehi(1),elo(2):ehi(2))
    real(amrex_real), intent(inout) ::  phib(blo(1):bhi(1),blo(2):bhi(2))
    real(amrex_real), intent(inout) ::  rhs (rlo(1):rhi(1),rlo(2):rhi(2))
    real(amrex_real), intent(in   ) :: bcent(clo(1):chi(1),clo(2):chi(2),2)
    integer         , intent(in   ) ::  flag(flo(1):fhi(1),flo(2):fhi(2))

    integer :: i,j
    real(amrex_real) :: x, y, r2, theta

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (is_covered_cell(flag(i,j))) then
             phie(i,j) = zero
             phib(i,j) = zero
          else
             x = (i+half)*dx(1) - half
             y = (j+half)*dx(2) - half
             theta = atan2(x,y)
             r2 = x**2 + y**2
             phie(i,j) = r2**2 * cos(three*theta)
             rhs(i,j) = -seven * r2 * cos(three*theta)
             if (is_single_valued_cell(flag(i,j))) then
                x = (i+half+bcent(i,j,1))*dx(1) - half
                y = (j+half+bcent(i,j,2))*dx(2) - half
                theta = atan2(x,y)
                r2 = x**2 + y**2
                phib(i,j) = r2**2 * cos(three*theta)
             else
                phib(i,j) = zero
             end if
          end if
       end do
    end do
  end subroutine mytest_set_phi_eb
#else
#endif

end module mytest
