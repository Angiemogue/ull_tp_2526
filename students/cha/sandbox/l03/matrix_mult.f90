program matrix_mult
    implicit none

    ! Matrix dimensions
    integer :: na, ma, nb, mb
    real :: A(3,4), B(3,3), C(3,4)
    integer :: i, j, k

    ! Define matrix sizes
    na = 3
    ma = 4
    nb = 3
    mb = 3

    ! Define matrices A and B
    A = reshape((/3.0,2.0,1.0, &
                  2.0,4.0,2.0, &
                  4.0,2.0,3.0, &
                  1.0,2.0,7.0/), (/na,ma/))

    B = reshape((/3.0,2.0,3.0, &
                  2.0,1.0,0.0, &
                  4.0,2.0,2.0/), (/nb,mb/))

    ! Check if B*A is possible (columns of B must equal rows of A)
    if (mb /= na) then
        print *, "ERROR: Cannot multiply B*A. Dimensions incompatible."
        print *, "B is ", nb, "x", mb, " and A is ", na, "x", ma
        stop
    end if

    ! Initialize result matrix
    C = 0.0

    ! Perform B*A multiplication
    do i = 1, nb
        do j = 1, ma
            do k = 1, mb
                C(i,j) = C(i,j) + B(i,k) * A(k,j)
            end do
        end do
    end do

    ! Print result
    print *, "Matrix B*A ="
    do i = 1, nb
        print '(4F6.1)', C(i,1:ma)
    end do

end program matrix_mult