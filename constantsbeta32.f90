module constants
    implicit none
    integer,parameter :: dp = selected_real_kind(p=15,r=200)
    !Constants are non-dimensional with non-dimensionalisation L = 4e5 m and T = 4e6 s
    real :: pi = 4*atan(1.)
    real :: tau = 8*atan(1.)
    ! real :: g = 3.92e8
    real :: Rd = 0.1
    real :: beta = 32
    ! real :: f = 400
end module constants
