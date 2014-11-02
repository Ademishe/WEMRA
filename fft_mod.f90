module fft_mod
implicit none
real*8, parameter     :: pi=3.141592653589793238460
contains

recursive subroutine fft(x)
    complex, dimension(:), intent(inout)    :: x
    complex                                 :: t
    integer                                 :: N
    integer                                 :: i
    complex, dimension(:), allocatable      :: even, odd
 
    N=size(x)
 
    if(N .le. 1) return
 
    allocate(odd((N+1)/2))
    allocate(even(N/2))

    odd =x(1:N:2)
    even=x(2:N:2)

    call fft(odd)
    call fft(even)

    do i=1,N/2
        t=exp( cmplx( 0.0, -2.0 * pi * real(i-1) / real(N) ) )*even(i)
        x(i) = odd(i) + t
        x(i+N/2) = odd(i) - t
    end do
 
    deallocate(odd)
    deallocate(even)
 
end subroutine fft
 
end module fft_mod
