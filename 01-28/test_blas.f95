
program test_blas
    implicit none

    !! Vector dot-product: result = x.y
    ! Declare variables
    real*8, dimension(4) :: x, y
    real*8 :: result
    integer :: n, incx, incy

    ! Declare the BLAS function
    real*8 :: ddot
    external :: ddot, dgemv

    !! Matrix-ector multiply: y <- alpha*A*x + beta*y
    real*8, dimension(4, 4) :: A
    real*8 :: alpha, beta
    integer :: m, lda

    !!===========================================================
    !! DDOT:
    ! Initialize vectors
    x = (/ 1.0, 2.0, 3.0, 4.0 /)
    y = (/ 5.0, 6.0, 7.0, 8.0 /)
    n = 4
    incx = 1
    incy = 1

    ! Call BLAS SDOT function
    result = ddot(n, x, incx, y, incy)

    ! Print the result
    print *, "Dot product: (x,y): ", result

    !! Now let's do a matrix multiplication
    !!===========================================================
    !! dgemv:
    ! Initialize the matrix and vectors
    A = reshape((/ 1.0, 2.0, 3.0, 4.0, &
                  5.0, 6.0, 7.0, 8.0, &
                  9.0, 10.0, 11.0, 12.0, &
                  13.0, 14.0, 15.0, 16.0 /), shape(A), order=(/ 2, 1 /))
    x = (/ 1.0, 1.0, 1.0, 1.0 /)
    y = (/ 0.0, 0.0, 0.0, 0.0 /)

    ! Parameters for SGEMV
    alpha = 1.0
    beta = 0.0
    m = 4            ! Number of rows in A
    n = 4            ! Number of columns in A
    lda = 4          ! Leading dimension of A
    incx = 1         ! Stride for x
    incy = 1         ! Stride for y

    ! Call BLAS _GEMV subroutine
    call dgemv('N', m, n, alpha, A, lda, x, incx, beta, y, incy)

    ! Print the result
    print *, "Matrix-vector product: ", y

end program test_blas

