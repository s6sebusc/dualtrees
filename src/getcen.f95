subroutine getcen( pyramid, nx, ny, nl, res )

integer :: nx, ny, nl, i, j
real( 8 ), dimension( nl, nx, ny, 6 ) :: pyramid
real( 8 ), dimension( nx, ny, 3 ) :: res
real( 8 ), dimension( nl, 6 ) :: a, b, c, e
real( 8 ) ::  pi, ang


pi = 4.D0*DATAN(1.D0)

do i = 1,nl
    do j = 1,6
        ang = pi*60.D0/180.D0*DFLOAT(j-1)
        a( i,j ) = cos( ang )
        b( i,j ) = sin( ang )
        c( i,j ) = i
    end do
end do

do i = 1,nx
    do j = 1,ny
        e = pyramid( :,i,j,: )
        e = e / sum( e )
       
        res( i,j,1 ) = sum( a*e )
        res( i,j,2 ) = sum( b*e )
        res( i,j,3 ) = sum( c*e )
        
    end do
end do

end subroutine getcen

