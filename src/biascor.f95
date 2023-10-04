subroutine biascor( sp, a, nx, ny, nl, nd, res )

integer :: nx, ny, nl, nd, i, j, k, l
real( 8 ), dimension( nl, nx, ny, nd ) :: sp, res
real( 8 ), dimension( nd*nl, nd*nl )   :: A
real( 8 ), dimension( nd*nl )          :: tmp

res( :, :, :, : ) = 0

do i = 1, nx
    do j = 1, ny
        tmp(:) = 0
        do k = 1, nl
            do l = 1, nd
                tmp( k + (l-1)*nl ) = sp( k, i, j, l )
            end do
        end do
        tmp = MATMUL( a, tmp )
        do k = 1, nl
            do l = 1, nd
                res( k, i, j, l ) = tmp( k + (l-1)*nl )
            end do
        end do
    end do
end do

end subroutine biascor
