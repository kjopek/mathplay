      module dense
            implicit none

            contains

            subroutine cholesky(A, n)
                  real, dimension(n,n), intent(in out) :: A
                  integer, intent(in)                  :: n
                  integer                              :: i, j, k

                  do j = 1, n
                        do k = 1, j - 1
                              do i = j, n
                                    A(i,j) = A(i,j) - A(i,k) * A(j,k)
                              end do
                        end do

                        A(j,j) = sqrt(A(j,j))

                        do i = j + 1, n
                              A(i,j) = A(i,j) / A(j,j)
                        end do
                  end do
            end subroutine cholesky
      end module dense

      program test_cholesky
            use dense
            implicit none

            integer                       :: n, i, j
            real, pointer, dimension(:,:) :: A

            n = 3
            allocate(A(n,n))

            A(1,1) =   4.0
            A(2,1) =  12.0
            A(3,1) = -16.0
            A(1,2) =  12.0
            A(2,2) =  37.0
            A(3,2) = -43.0
            A(1,3) = -16.0
            A(2,3) = -43.0
            A(3,3) =  98.0

            call cholesky(A, n)

            do i = 1, n
                  write (*,*) A(i, 1:i)
            end do

            deallocate(A)
      end program test_cholesky



