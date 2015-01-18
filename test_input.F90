subroutine plus(A,N,B)
implicit none
integer, intent(in) :: N
double precision, dimension(1:N,1:N) :: A, B
write(*,*) A
B = A + A
end subroutine plus
