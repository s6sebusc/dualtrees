subroutine getcen( pyramid, nx, ny, nl, res )

integer :: nx, ny, nl, i, j
real( 8 ), dimension( nx, ny, nl, 6 ) :: pyramid
real( 8 ), dimension( nx, ny, 3 ) :: res
real( 8 ), dimension( nl, 6 ) :: a, b, c, e
real( 8 ) :: x, y, z, rho, phi, pi, ang


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
        e = pyramid( i,j,:,: )
        e = e / sum( e )
        
        x = sum( a*e )
        y = sum( b*e )
        z = sum( c*e )
        
        rho = sqrt( x**2 + y**2 )
        phi = ( atan2( y,x )*180.D0/pi )/2.D0
        phi = phi + 15.D0
        if ( phi < 0. ) then
            phi =  phi + 180.D0
        end if
        
        res( i,j,1 ) = rho
        res( i,j,2 ) = phi
        res( i,j,3 ) = z
        
    end do
end do

end subroutine getcen

