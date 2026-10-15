!-------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
! Some of the content of this file has been produced with the assistance of
! Met Office Github Copilot Enterprise
!-------------------------------------------------------------------------------
!> @brief Compute the orographic slopes and slope factor at W2v locations.
!> @details For lowest-order W2v DoFs, computes the quantities needed to
!!          convert between the flux through a sloped horizontal face and
!!          the true vertical component of a vector:
!!            dz_dx = (rhat . J1) * |J2 x J3| / det(J)
!!            dz_dy = (rhat . J2) * |J3 x J1| / det(J)
!!            slope_factor = |J3| * |J1 x J2| / det(J)
!!          where Ji is the i-th column of the Jacobian (in Cartesian
!!          coordinates) and rhat = J3 / |J3| is the local vertical
!!          direction.
!!          On an orthogonal planar mesh, dz_dx and dz_dy are the slopes of
!!          the cell faces and slope_factor = sqrt(1 + dz_dx**2 + dz_dy**2).
module sci_calc_slopes_at_w2v_kernel_mod

  use argument_mod,       only: arg_type, func_type,                           &
                                GH_FIELD, GH_SCALAR,                           &
                                GH_REAL, GH_INTEGER,                           &
                                GH_READ, GH_WRITE,                             &
                                ANY_SPACE_1, ANY_DISCONTINUOUS_SPACE_3,        &
                                GH_BASIS, GH_DIFF_BASIS,                       &
                                CELL_COLUMN, GH_EVALUATOR
  use constants_mod,      only: r_def, i_def
  use fs_continuity_mod,  only: W2v
  use kernel_mod,         only: kernel_type

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  type, public, extends(kernel_type) :: calc_slopes_at_w2v_kernel_type
    private
    type(arg_type) :: meta_args(9) = (/                                        &
        arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, W2v),                       &
        arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, W2v),                       &
        arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, W2v),                       &
        arg_type(GH_FIELD*3, GH_REAL,    GH_READ,  ANY_SPACE_1),               &
        arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_3), &
        arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                             &
        arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                             &
        arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                             &
        arg_type(GH_SCALAR,  GH_REAL,    GH_READ)                              &
    /)
    type(func_type) :: meta_funcs(1) = (/                                      &
        func_type(ANY_SPACE_1, GH_BASIS, GH_DIFF_BASIS)                        &
    /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass :: calc_slopes_at_w2v_code
  end type calc_slopes_at_w2v_kernel_type

  !-----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  public :: calc_slopes_at_w2v_code

contains

!> @param[in]     nlayers        Number of layers
!> @param[in,out] dz_dx          Slope in the first horizontal direction
!> @param[in,out] dz_dy          Slope in the second horizontal direction
!> @param[in,out] slope_factor   Ratio of sloped to horizontal face area
!> @param[in]     chi1           Coordinates in the first direction
!> @param[in]     chi2           Coordinates in the second direction
!> @param[in]     chi3           Coordinates in the third direction
!> @param[in]     panel_id       Field giving the ID for mesh panels
!> @param[in]     coord_system   Finite-element coordinate system enumeration
!> @param[in]     geometry       Mesh geometry enumeration
!> @param[in]     topology       Mesh topology enumeration
!> @param[in]     scaled_radius  Scaled planetary radius
!> @param[in]     ndf_w2v        Num DoFs per cell for W2v
!> @param[in]     undf_w2v       Num DoFs in this partition for W2v
!> @param[in]     map_w2v        Map of base cell DoFs for W2v
!> @param[in]     ndf_chi        Num DoFs per cell for chi
!> @param[in]     undf_chi       Num DoFs in this partition for chi
!> @param[in]     map_chi        Map of base cell DoFs for chi
!> @param[in]     basis_chi      Coordinate basis at W2v nodal points
!> @param[in]     diff_basis_chi Coordinate diff basis at W2v nodal points
!> @param[in]     ndf_pid        Num DoFs per cell for panel_id
!> @param[in]     undf_pid       Num DoFs in this partition for panel_id
!> @param[in]     map_pid        Map of base cell DoFs for panel_id
subroutine calc_slopes_at_w2v_code(                                            &
    nlayers,                                                                   &
    dz_dx, dz_dy, slope_factor,                                                &
    chi1, chi2, chi3, panel_id,                                                &
    coord_system, geometry, topology, scaled_radius,                           &
    ndf_w2v, undf_w2v, map_w2v,                                                &
    ndf_chi, undf_chi, map_chi, basis_chi, diff_basis_chi,                     &
    ndf_pid, undf_pid, map_pid                                                 &
)

  use cross_product_mod,           only: cross_product
  use sci_coordinate_jacobian_mod, only: coordinate_jacobian

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ndf_w2v
  integer(kind=i_def), intent(in)    :: undf_w2v
  integer(kind=i_def), intent(in)    :: ndf_chi
  integer(kind=i_def), intent(in)    :: undf_chi
  integer(kind=i_def), intent(in)    :: ndf_pid
  integer(kind=i_def), intent(in)    :: undf_pid
  integer(kind=i_def), intent(in)    :: coord_system
  integer(kind=i_def), intent(in)    :: geometry
  integer(kind=i_def), intent(in)    :: topology
  real(kind=r_def),    intent(in)    :: scaled_radius

  real(kind=r_def),    intent(inout) :: dz_dx(undf_w2v)
  real(kind=r_def),    intent(inout) :: dz_dy(undf_w2v)
  real(kind=r_def),    intent(inout) :: slope_factor(undf_w2v)
  real(kind=r_def),    intent(in)    :: chi1(undf_chi)
  real(kind=r_def),    intent(in)    :: chi2(undf_chi)
  real(kind=r_def),    intent(in)    :: chi3(undf_chi)
  real(kind=r_def),    intent(in)    :: panel_id(undf_pid)

  integer(kind=i_def), intent(in)    :: map_w2v(ndf_w2v)
  integer(kind=i_def), intent(in)    :: map_chi(ndf_chi)
  integer(kind=i_def), intent(in)    :: map_pid(ndf_pid)
  real(kind=r_def),    intent(in)    :: basis_chi(1,ndf_chi,ndf_w2v)
  real(kind=r_def),    intent(in)    :: diff_basis_chi(3,ndf_chi,ndf_w2v)

  ! Internal variables
  integer(kind=i_def) :: df, k, ipanel
  real(kind=r_def)    :: chi1_e(ndf_chi), chi2_e(ndf_chi), chi3_e(ndf_chi)
  real(kind=r_def)    :: jacobian(3,3,1)
  real(kind=r_def)    :: dj(1)
  real(kind=r_def)    :: j3_norm
  real(kind=r_def)    :: rhat(3)

  ipanel = int(panel_id(map_pid(1)), i_def)

  ! Lowest order only: each cell has a bottom (df=1) and a top (df=2) W2v point.
  ! W2v is continuous in the vertical, so the top of one cell is the bottom of
  ! the cell above. Calculate every bottom point once, then the domain top.
  do k = 0, nlayers-1
    do df = 1, ndf_chi
      chi1_e(df) = chi1(map_chi(df) + k)
      chi2_e(df) = chi2(map_chi(df) + k)
      chi3_e(df) = chi3(map_chi(df) + k)
    end do

    ! Evaluate the Jacobian only at the first (bottom) W2v point
    call coordinate_jacobian(coord_system, geometry, topology, scaled_radius,  &
                             ndf_chi, 1_i_def, chi1_e, chi2_e, chi3_e,         &
                             ipanel, basis_chi(:,:,1:1),                       &
                             diff_basis_chi(:,:,1:1), jacobian, dj)

    ! rhat is the local vertical, taken from the third column of the Jacobian
    j3_norm = norm2(jacobian(:,3,1))
    rhat = jacobian(:,3,1) / j3_norm

    dz_dx(map_w2v(1)+k) = dot_product(rhat, jacobian(:,1,1))                   &
        * norm2(cross_product(jacobian(:,2,1), jacobian(:,3,1))) / dj(1)
    dz_dy(map_w2v(1)+k) = dot_product(rhat, jacobian(:,2,1))                   &
        * norm2(cross_product(jacobian(:,3,1), jacobian(:,1,1))) / dj(1)
    slope_factor(map_w2v(1)+k) = j3_norm                                       &
        * norm2(cross_product(jacobian(:,1,1), jacobian(:,2,1))) / dj(1)
  end do

  ! Top of the domain: the top (second) W2v point of the highest cell
  k = nlayers-1
  do df = 1, ndf_chi
    chi1_e(df) = chi1(map_chi(df) + k)
    chi2_e(df) = chi2(map_chi(df) + k)
    chi3_e(df) = chi3(map_chi(df) + k)
  end do

  call coordinate_jacobian(coord_system, geometry, topology, scaled_radius,    &
                           ndf_chi, 1_i_def, chi1_e, chi2_e, chi3_e,           &
                           ipanel, basis_chi(:,:,2:2),                         &
                           diff_basis_chi(:,:,2:2), jacobian, dj)

  j3_norm = norm2(jacobian(:,3,1))
  rhat = jacobian(:,3,1) / j3_norm

  dz_dx(map_w2v(2)+k) = dot_product(rhat, jacobian(:,1,1))                     &
      * norm2(cross_product(jacobian(:,2,1), jacobian(:,3,1))) / dj(1)
  dz_dy(map_w2v(2)+k) = dot_product(rhat, jacobian(:,2,1))                     &
      * norm2(cross_product(jacobian(:,3,1), jacobian(:,1,1))) / dj(1)
  slope_factor(map_w2v(2)+k) = j3_norm                                         &
      * norm2(cross_product(jacobian(:,1,1), jacobian(:,2,1))) / dj(1)

end subroutine calc_slopes_at_w2v_code

end module sci_calc_slopes_at_w2v_kernel_mod
