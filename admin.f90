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

    function MaduLawrence(k, N, mu) result(ML)
        real(dp) :: k, ML, mu, a, b, N
        
        ML = 0
        if ((mu > 0) .and. (mu .le. 0.5)) then
            if (k == 0) then
                ML = 0
            else if (k < mu*N) then
                a = mu*N*(1/k + 1/(k-mu*N))
                ML = 1/(exp(a)+1)
            else if (k .le. (1-mu)*N) then
                ML = 1
            else if (k < N) then
                b = mu*N*(1/(N-k) + 1/((1-mu)*N-k))
                ML = 1/(exp(b)+1)
            else if (k == N) then
                ML = 0
            end if
        else
            ML = 1
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