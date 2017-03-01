

module i2sint

  use constants, only : dp

  implicit none
  public

  ! sparse storage for two-electron integrals

  ! indices are stored as 2-byte integers
  ! the actual integral is stored as a double precision number

  type i2s
    integer*2 :: ii, ij, ik, il
    real(kind=dp) :: int
  end type i2s

end module i2sint


