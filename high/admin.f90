module admin
    use constants
    implicit none

    contains
    function linspace(start_point, end_point, num_points) result(points_array)
        real :: start_point, end_point, increment, new_point
        integer :: num_points, counter
        real, dimension(:), allocatable :: points_array
        
        allocate(points_array(num_points))
        increment = (end_point - start_point)/(num_points-1.)
        new_point = start_point
        do counter = 1, num_points
            points_array(counter) = new_point
            new_point = new_point + increment
        end do

    end function linspace

    function MaduLawrence(k, N) result(ML)
        real(dp) :: k, ML, N
        ! k = t 
        ! mu = mu
        ! N = vu - 1
        ! N + 1 = vu
        ML = 0
        N = N + 1 
        if (k == 0) then
                ML = 0
        else if ((0 < k) .and. (k < 3)) then
                ML = exp(1-1/(1-(k-1)**2))
        else if ((3 .le. k) .and. (k .le. N-2)) then
                ML = 1 
        else if ((N-2 < k) .and. (k < N + 1)) then
                ML = exp(1-1/(1-(k-N+2)**2))
        else if (k==N+1) then
                ML = 0 
        end if 

    end function MaduLawrence
    
    function okubo_class(div, Okw) result(classif)
        real(dp) :: div, Okw, lambda1, lambda2
        integer :: classif

        classif = 0

        if(Okw > 0) then
            lambda1 = 0.5*(div + Okw**0.5)
            lambda2 = 0.5*(div - Okw**0.5)
            if((lambda1 > 0) .and. (lambda2 < 0)) then
                classif = 3
            else if(lambda1 < 0) then
                classif = 2
            else if(lambda2 > 0) then
                classif = 4
            end if
        else if(Okw < 0) then
            if(div < 0) then
                classif = 1
            else if (div > 0) then
                classif = 5
            else if (div == 0) then
                classif = 6
            end if
        end if

    end function okubo_class
end module admin
