

module constants

  implicit none
  public
  save

  ! double precision

  integer :: dp
  parameter ( dp = selected_real_kind (15, 307) )

  ! common constants

  real(kind=dp) :: d0, d1, d2, d3
  parameter ( d0 = 0.0e0_dp, d1 = 1.0e0_dp )
  parameter ( d2 = 2.0e0_dp, d3 = 3.0e0_dp )

  complex(kind=dp) :: z0, z1, zi
  parameter ( z0 = (0.0e0_dp, 0.0e0_dp) )
  parameter ( z1 = (1.0e0_dp, 0.0e0_dp) )
  parameter ( zi = (0.0e0_dp, 1.0e0_dp) )

end module constants


