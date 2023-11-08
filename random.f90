module random
    use constants
    implicit none

    contains
    function random_number1() result(r)
        real(dp) :: r

        call random_number(r)
    end function random_number1

    function random_numbers(num) result(r)
        real(dp), allocatable :: r(:)
        integer :: num

        allocate(r(num))
        call random_number(r)
    end function random_numbers

    function random_matrix(num,num2) result(r)
        real(dp), allocatable :: r(:,:)
        integer :: num, num2

        allocate(r(num,num2))
        call random_number(r)
    end function random_matrix
end module random