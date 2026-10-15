!-------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
! Some of the content of this file has been produced with the assistance of
! Met Office Github Copilot Enterprise
!-------------------------------------------------------------------------------
!> @brief Compute the true vertical wind from the fluxes of a W2 field.
!> @details Over sloping coordinate surfaces the W2v flux is not the vertical
!!          component of the wind. For lowest-order fields, the true vertical
!!          wind at each W2v point is calculated as
!!            w = F_v / dA_v * slope_factor + dz_dx * u_bar + dz_dy * v_bar
!!          where F_v / dA_v is the wind normal to the sloped face, and u_bar
!!          and v_bar are the horizontal winds (the W2h fluxes divided by their
!!          face areas) averaged to the W2v point from the cells either side.
module sci_split_vert_wind_kernel_mod

  use argument_mod,       only: arg_type,                                      &
                                GH_FIELD, GH_REAL,                             &
                                GH_READ, GH_WRITE,                             &
                                CELL_COLUMN
  use constants_mod,      only: r_double, r_single, i_def
  use fs_continuity_mod,  only: W2, W2h, W2v
  use kernel_mod,         only: kernel_type

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  type, public, extends(kernel_type) :: split_vert_wind_kernel_type
    private
    type(arg_type) :: meta_args(7) = (/                                        &
        arg_type(GH_FIELD, GH_REAL, GH_WRITE, W2v),                            &
        arg_type(GH_FIELD, GH_REAL, GH_READ,  W2h),                            &
        arg_type(GH_FIELD, GH_REAL, GH_READ,  W2v),                            &
        arg_type(GH_FIELD, GH_REAL, GH_READ,  W2),                             &
        arg_type(GH_FIELD, GH_REAL, GH_READ,  W2v),                            &
        arg_type(GH_FIELD, GH_REAL, GH_READ,  W2v),                            &
        arg_type(GH_FIELD, GH_REAL, GH_READ,  W2v)                             &
    /)
    integer :: operates_on = CELL_COLUMN
  end type split_vert_wind_kernel_type

  !-----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  public :: split_vert_wind_code

  ! Generic interface for real32 and real64 types
  interface split_vert_wind_code
    module procedure                                                           &
        split_vert_wind_code_r_single,                                         &
        split_vert_wind_code_r_double
  end interface

contains

!> @brief Compute the true vertical wind from the fluxes of a W2 field
!> @param[in]     nlayers      Number of layers
!> @param[in,out] w_vert       True vertical wind at W2v points
!> @param[in]     u_hori       Horizontal wind fluxes on W2h
!> @param[in]     w_flux       Vertical wind fluxes on W2v
!> @param[in]     dA_at_w2     Face areas at W2 points
!> @param[in]     dz_dx        Slope in the first horizontal direction
!> @param[in]     dz_dy        Slope in the second horizontal direction
!> @param[in]     slope_factor Ratio of sloped to horizontal face area
!> @param[in]     ndf_w2v      Num DoFs per cell for W2v
!> @param[in]     undf_w2v     Num DoFs in this partition for W2v
!> @param[in]     map_w2v      Map of base cell DoFs for W2v
!> @param[in]     ndf_w2h      Num DoFs per cell for W2h
!> @param[in]     undf_w2h     Num DoFs in this partition for W2h
!> @param[in]     map_w2h      Map of base cell DoFs for W2h
!> @param[in]     ndf_w2       Num DoFs per cell for W2
!> @param[in]     undf_w2      Num DoFs in this partition for W2
!> @param[in]     map_w2       Map of base cell DoFs for W2

! R_DOUBLE PRECISION
! ==================
subroutine split_vert_wind_code_r_double(                                      &
    nlayers,                                                                   &
    w_vert, u_hori, w_flux, dA_at_w2, dz_dx, dz_dy, slope_factor,              &
    ndf_w2v, undf_w2v, map_w2v,                                                &
    ndf_w2h, undf_w2h, map_w2h,                                                &
    ndf_w2, undf_w2, map_w2                                                    &
)

  use reference_element_mod, only: W, S, E, N, B

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ndf_w2v
  integer(kind=i_def), intent(in)    :: undf_w2v
  integer(kind=i_def), intent(in)    :: ndf_w2h
  integer(kind=i_def), intent(in)    :: undf_w2h
  integer(kind=i_def), intent(in)    :: ndf_w2
  integer(kind=i_def), intent(in)    :: undf_w2

  real(kind=r_double), intent(inout) :: w_vert(undf_w2v)
  real(kind=r_double), intent(in)    :: u_hori(undf_w2h)
  real(kind=r_double), intent(in)    :: w_flux(undf_w2v)
  real(kind=r_double), intent(in)    :: dA_at_w2(undf_w2)
  real(kind=r_double), intent(in)    :: dz_dx(undf_w2v)
  real(kind=r_double), intent(in)    :: dz_dy(undf_w2v)
  real(kind=r_double), intent(in)    :: slope_factor(undf_w2v)

  integer(kind=i_def), intent(in)    :: map_w2v(ndf_w2v)
  integer(kind=i_def), intent(in)    :: map_w2h(ndf_w2h)
  integer(kind=i_def), intent(in)    :: map_w2(ndf_w2)

  ! Internal variables
  integer(kind=i_def) :: k, k_below, k_above
  real(kind=r_double) :: u_cell(nlayers), v_cell(nlayers)
  real(kind=r_double) :: u_bar, v_bar

  ! Wind normal to the vertical faces, averaged to the centre of each cell.
  ! Lowest order only: the W2h DoFs are W, S, E, N and W2 has B, T after them.
  do k = 0, nlayers-1
    u_cell(k+1) = 0.5_r_double                                                 &
        * (u_hori(map_w2h(W)+k) / dA_at_w2(map_w2(W)+k)                        &
           + u_hori(map_w2h(E)+k) / dA_at_w2(map_w2(E)+k))
    v_cell(k+1) = 0.5_r_double                                                 &
        * (u_hori(map_w2h(S)+k) / dA_at_w2(map_w2(S)+k)                        &
           + u_hori(map_w2h(N)+k) / dA_at_w2(map_w2(N)+k))
  end do

  ! Level k has cell k below it and cell k+1 above it (1-based cell indices).
  ! At the bottom and top of the column the single cell is used twice.
  do k = 0, nlayers
    k_below = max(k, 1)
    k_above = min(k+1, nlayers)
    u_bar = 0.5_r_double * (u_cell(k_below) + u_cell(k_above))
    v_bar = 0.5_r_double * (v_cell(k_below) + v_cell(k_above))

    w_vert(map_w2v(1)+k) = w_flux(map_w2v(1)+k) / dA_at_w2(map_w2(B)+k)        &
        * slope_factor(map_w2v(1)+k)                                           &
        + dz_dx(map_w2v(1)+k) * u_bar + dz_dy(map_w2v(1)+k) * v_bar
  end do

end subroutine split_vert_wind_code_r_double

! R_SINGLE PRECISION
! ==================
subroutine split_vert_wind_code_r_single(                                      &
    nlayers,                                                                   &
    w_vert, u_hori, w_flux, dA_at_w2, dz_dx, dz_dy, slope_factor,              &
    ndf_w2v, undf_w2v, map_w2v,                                                &
    ndf_w2h, undf_w2h, map_w2h,                                                &
    ndf_w2, undf_w2, map_w2                                                    &
)

  use reference_element_mod, only: W, S, E, N, B

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ndf_w2v
  integer(kind=i_def), intent(in)    :: undf_w2v
  integer(kind=i_def), intent(in)    :: ndf_w2h
  integer(kind=i_def), intent(in)    :: undf_w2h
  integer(kind=i_def), intent(in)    :: ndf_w2
  integer(kind=i_def), intent(in)    :: undf_w2

  real(kind=r_single), intent(inout) :: w_vert(undf_w2v)
  real(kind=r_single), intent(in)    :: u_hori(undf_w2h)
  real(kind=r_single), intent(in)    :: w_flux(undf_w2v)
  real(kind=r_single), intent(in)    :: dA_at_w2(undf_w2)
  real(kind=r_single), intent(in)    :: dz_dx(undf_w2v)
  real(kind=r_single), intent(in)    :: dz_dy(undf_w2v)
  real(kind=r_single), intent(in)    :: slope_factor(undf_w2v)

  integer(kind=i_def), intent(in)    :: map_w2v(ndf_w2v)
  integer(kind=i_def), intent(in)    :: map_w2h(ndf_w2h)
  integer(kind=i_def), intent(in)    :: map_w2(ndf_w2)

  ! Internal variables
  integer(kind=i_def) :: k, k_below, k_above
  real(kind=r_single) :: u_cell(nlayers), v_cell(nlayers)
  real(kind=r_single) :: u_bar, v_bar

  ! Wind normal to the vertical faces, averaged to the centre of each cell.
  ! Lowest order only: the W2h DoFs are W, S, E, N and W2 has B, T after them.
  do k = 0, nlayers-1
    u_cell(k+1) = 0.5_r_single                                                 &
        * (u_hori(map_w2h(W)+k) / dA_at_w2(map_w2(W)+k)                        &
           + u_hori(map_w2h(E)+k) / dA_at_w2(map_w2(E)+k))
    v_cell(k+1) = 0.5_r_single                                                 &
        * (u_hori(map_w2h(S)+k) / dA_at_w2(map_w2(S)+k)                        &
           + u_hori(map_w2h(N)+k) / dA_at_w2(map_w2(N)+k))
  end do

  ! Level k has cell k below it and cell k+1 above it (1-based cell indices).
  ! At the bottom and top of the column the single cell is used twice.
  do k = 0, nlayers
    k_below = max(k, 1)
    k_above = min(k+1, nlayers)
    u_bar = 0.5_r_single * (u_cell(k_below) + u_cell(k_above))
    v_bar = 0.5_r_single * (v_cell(k_below) + v_cell(k_above))

    w_vert(map_w2v(1)+k) = w_flux(map_w2v(1)+k) / dA_at_w2(map_w2(B)+k)        &
        * slope_factor(map_w2v(1)+k)                                           &
        + dz_dx(map_w2v(1)+k) * u_bar + dz_dy(map_w2v(1)+k) * v_bar
  end do

end subroutine split_vert_wind_code_r_single

end module sci_split_vert_wind_kernel_mod
