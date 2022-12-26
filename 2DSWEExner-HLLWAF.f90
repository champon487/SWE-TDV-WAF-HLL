subroutine propagation(sl,sc,sr,u,v,h,g,snm,spec,dd,tsc,poro,hmin)
    implicit none
    double precision, intent(in)  :: u, v, h, g, snm, spec, dd, poro, tsc, hmin
    double precision, intent(out) :: sl, sc, sr

    double precision :: cc, cf, vv, tt, aa
    
!    cl = dsqrt(g*wh(i,j))
!    cr = dsqrt(g*wh(i+1,j))
!    ust = 0.5d0*(u(i,j)+u(i+1,j))+cl-cr
!    cst = 0.5*(cl+cr)+0.25*(u(i,j)-u(i+1,j))
!    sl(i,j) = min(u(i,j)-cl,ust-cst)
!    sr(i,j) = max(u(i+1,j)+cr,ust+cst)

    cc = dsqrt(g*h)
    vv = dsqrt(u**2.d0+v**2.d0)
    if( vv>1e-5 .and. h>hmin ) then
        cf = g*snm**2.d0/h**(1.d0/3.d0)
        tt = cf*vv**2.d0/(spec*g*dd)
    else
        cf = 0.d0
        tt = 0.d0
    end if

    if( tt>tsc ) then
!       aa = 4.d0*cf**1.5d0*dsqrt(dd/spec/g)*(1.d0-tsc/tt)**1.5d0*u/h*cc**2.d0*(3.d0*u**2.d0+v**2.d0)/(1.d0-poro)
       aa = 4.d0/(1.d0-poro)*(tt-tsc)**0.5d0*cc**2.d0*(1.5d0*cf*dsqrt(dd/spec/g)*u/h*(3.d0*u**2.d0+v**2.d0)/dsqrt(u**2.d0+v**2.d0)   &
            +(tt-tsc)*dsqrt(spec*g*dd**3.d0)*u*v**2.d0/h/(u**2.d0+v**2.d0)**1.5d0)
!       cf**1.5d0*dsqrt(dd/spec/g)*(1.d0-tsc/tt)**1.5d0*u/h*cc**2.d0*(3.d0*u**2.d0+v**2.d0)/(1.d0-poro)
    else
        aa = 0.d0
    end if

    if( u>=0.d0 ) then
        sr = u+cc
        sc = 0.5d0*(u-cc+dsqrt((u-cc)**2.d0+4.d0*aa/(u+cc)))
        sl = 0.5d0*(u-cc-dsqrt((u-cc)**2.d0+4.d0*aa/(u+cc)))
    else
        sr = 0.5d0*(u+cc+dsqrt((u+cc)**2.d0+4.d0*aa/(u-cc)))
        sc = 0.5d0*(u+cc-dsqrt((u+cc)**2.d0+4.d0*aa/(u-cc)))
        sl = u-cc
    end if

end subroutine

subroutine hllx1storder(f,ff,sl,sr,u,nx,ny)
    implicit none
    integer, intent(in) :: nx, ny
    double precision, dimension(0:nx,0:ny), intent(in)  :: ff, sl, sr, u
    double precision, dimension(0:nx,0:ny), intent(out) :: f

    integer :: i, j

!$omp do private(i,j) 
    do j=1,ny-1
        do i=0,nx-1
            if( sl(i,j)>0.d0 ) then
                f(i,j) = ff(i,j)
            else if( sr(i,j)<0.d0 ) then
                f(i,j) = ff(i+1,j)
            else
                f(i,j) = (ff(i,j)*sr(i,j)-ff(i+1,j)*sl(i,j)+sr(i,j)*sl(i,j)*(u(i+1,j)-u(i,j)))/(sr(i,j)-sl(i,j))
            end if
        end do
    end do

end subroutine

subroutine hllzx1storder(f,ff,sl,sc,sr,z,u,snm,dx,qx,h,nx,ny)
    implicit none
    integer, intent(in) :: nx, ny
    double precision, intent(in) :: snm, dx
    double precision, dimension(0:nx,0:ny), intent(in)  :: ff, sl, sc, sr, u, z, qx, h
    double precision, dimension(0:nx,0:ny), intent(out) :: f

    integer :: i, j
    double precision :: ust, hst, s0dx, zzz

!$omp do private(i,j,ust,hst,s0dx,zzz)
    do j=1,ny-1
        do i=0,nx-1
            if( sl(i,j)>0.d0 ) then
                f(i,j) = ff(i,j)
            else if( sr(i,j)<0.d0 ) then
                f(i,j) = ff(i+1,j)
            else
                ust = (u(i,j)+u(i+1,j))*0.5d0
                hst = (h(i,j)+h(i+1,j))*0.5d0
                s0dx = snm**2.*ust**2./hst**(4./3.)*dx
                zzz = (z(i+1,j)-z(i,j))+s0dx
                if( ust>0.d0 ) then
!                    f(i,j) = (ff(i,j)*sc(i,j)-ff(i+1,j)*sl(i,j)+sc(i,j)*sl(i,j)*zzz)/(sc(i,j)-sl(i,j))
                    f(i,j) = (ff(i,j)*sc(i,j)-ff(i+1,j)*sl(i,j)+sc(i,j)*sl(i,j)*(z(i+1,j)-z(i,j)))/(sc(i,j)-sl(i,j))
                else
!                    f(i,j) = (-ff(i,j)*sr(i,j)+ff(i+1,j)*sc(i,j)-sc(i,j)*sr(i,j)*zzz)/(sc(i,j)-sr(i,j))
                    f(i,j) = (-ff(i,j)*sr(i,j)+ff(i+1,j)*sc(i,j)-sc(i,j)*sr(i,j)*(z(i+1,j)-z(i,j)))/(sc(i,j)-sr(i,j))
                end if
            end if
        end do
    end do

end subroutine

subroutine wafx(ff,fsta,f,sl,sr,dx,dt,nx,ny)
    implicit none
    integer, intent(in) :: nx, ny
    double precision, intent(in) :: dx, dt
    double precision, dimension(0:nx,0:ny), intent(in) :: fsta, f, sl, sr
    double precision, dimension(0:nx,0:ny), intent(out) :: ff

    integer :: i, j

    do j=1,ny-1
        do i=0,nx-1
            if( sl(i,j)>0.d0 .or. sr(i,j)<0. ) then
                ff(i,j) = fsta(i,j)
            else
                ff(i,j) = 0.5d0*(f(i,j)*(1.d0+sl(i,j)*dt/dx)+fsta(i,j)*(sr(i,j)-sl(i,j))*dt/dx+f(i+1,j)*(1.d0-sr(i,j)*dt/dx) )
            end if
        end do
    end do

end subroutine


subroutine phical(phi,c,dql,dqc,dqr)
    implicit none
    double precision, intent(in)  :: dql,dqc,dqr,c
    double precision, intent(out) :: phi

    double precision :: r

    if (c>0) then
        if (dqc==0.d0) then
            r = 1.d0
        else
            r = dql/dqc
        end if
    else
        if (dqc==0.d0) then
            r = 1.
        else
            r = dqr/dqc
        end if
    end if
            
            
!    if (r<0.) then
!        phi = 1.
!    else
!        phi = 1.-(1.-dabs(c))*r*(1.+r)/(1.+r**2.)
!!        phi = 1.-(1.-np.abs(c))*2.*r/(1.+r)
!    end if

    ! stable
            
    if (r<0.) then
        phi = 1.d0
    else if ( r>0 .and. r<1. ) then
        phi = 1.d0-(1.d0-dabs(c))*r
    else
        phi = dabs(c)
    end if
            
!    phi = np.abs(c)
            
!    if (r<=0.d0) then
!        phi = 1.d0
!    else if (r>0. .and. r<0.5) then
!        phi = 1.d0-2.d0*(1.-dabs(c))*r
!    else if (r>0.5 .and. r<1) then
!        phi = dabs(c)
!    else if (r>1. .and. r<1.5) then
!        phi = 1.-(1.-dabs(c))*r
!    else
!        phi = 2.*dabs(c)-1.d0
!    end if

end subroutine

subroutine tvdwafx(ff,fsta,f,dq,sl,sr,dx,dt,nx,ny)
        implicit none
        integer, intent(in) :: nx, ny
        double precision, intent(in) :: dx, dt
        double precision, dimension(0:nx,0:ny), intent(in) :: fsta, f, sl, sr,dq
        double precision, dimension(0:nx,0:ny), intent(out) :: ff
    
        integer :: i, j
        double precision :: cr, cl, ull, urr, phi
    
!$omp single
        do j=1,ny-1
            i=0
!            if( sl(i,j)>0.d0 .or. sr(i,j)<0. ) then
                ff(i,j) = fsta(i,j)
!            else
!                ff(i,j) = 0.5d0*(f(i,j)*(1.d0+sl(i,j)*dt/dx)+fsta(i,j)*(sr(i,j)-sl(i,j))*dt/dx+f(i+1,j)*(1.d0-sr(i,j)*dt/dx) )
!            end if

            i=nx-1
!            if( sl(i,j)>0.d0 .or. sr(i,j)<0. ) then
                ff(i,j) = fsta(i,j)
!            else
!                ff(i,j) = 0.5d0*(f(i,j)*(1.d0+sl(i,j)*dt/dx)+fsta(i,j)*(sr(i,j)-sl(i,j))*dt/dx+f(i+1,j)*(1.d0-sr(i,j)*dt/dx) )
!            end if
        end do
!$omp end single

!$omp do private(i,j,cl,ull,cr,urr,phi)
        do j=1,ny-1
            do i=1,nx-2
                if( sl(i,j)>0.d0 .or. sr(i,j)<0. ) then
                    ff(i,j) = fsta(i,j)
                else
                    cl = sl(i,j)*dt/dx
                    call phical(phi,cl,dq(i-1,j),dq(i,j),dq(i+1,j))
                    ull = sign(1.d0,cl)*0.5d0*phi*(fsta(i,j)-f(i,j))
                    
                    cr = sr(i,j)*dt/dx
                    call phical(phi,cr,dq(i-1,j),dq(i,j),dq(i+1,j))
                    urr = sign(1.d0,cr)*0.5d0*phi*(f(i+1,j)-fsta(i,j))

                    ff(i,j) = 0.5d0*(f(i,j)+f(i+1,j))-ull-urr
                end if
            end do
        end do

end subroutine

subroutine hlly1storder(f,ff,sl,sr,u,nx,ny)
    implicit none
    integer, intent(in) :: nx, ny
    double precision, dimension(0:nx,0:ny), intent(in)  :: ff, sl, sr, u
    double precision, dimension(0:nx,0:ny), intent(out) :: f

    integer :: i, j

!$omp do private(i,j)
    do j=0,ny-1
        do i=1,nx-1
            if( sl(i,j)>0.d0 ) then
                f(i,j) = ff(i,j)
            else if( sr(i,j)<0.d0 ) then
                f(i,j) = ff(i,j+1)
            else
                f(i,j) = (ff(i,j)*sr(i,j)-ff(i,j+1)*sl(i,j)+sr(i,j)*sl(i,j)*(u(i,j+1)-u(i,j)))/(sr(i,j)-sl(i,j))
            end if
        end do
    end do

end subroutine

subroutine hllzy1storder(f,ff,sl,sc,sr,z,u,nx,ny)
    implicit none
    integer, intent(in) :: nx, ny
    double precision, dimension(0:nx,0:ny), intent(in)  :: ff, sl, sc, sr, u, z
    double precision, dimension(0:nx,0:ny), intent(out) :: f

    integer :: i, j
    double precision :: ust

!$omp do private(i,j,ust)
    do j=0,ny-1
        do i=1,nx-1
            if( sl(i,j)>0.d0 ) then
                f(i,j) = ff(i,j)
            else if( sr(i,j)<0.d0 ) then
                f(i,j) = ff(i,j+1)
            else
                ust = (u(i,j)+u(i,j+1))*0.5d0
                if( ust>0.d0 ) then
                    f(i,j) = (ff(i,j)*sc(i,j)-ff(i,j+1)*sl(i,j)+sc(i,j)*sl(i,j)*(z(i,j+1)-z(i,j)))/(sc(i,j)-sl(i,j))
                else
                    f(i,j) = (-ff(i,j)*sr(i,j)+ff(i,j+1)*sc(i,j)-sc(i,j)*sr(i,j)*(z(i,j+1)-z(i,j)))/(sc(i,j)-sr(i,j))
                end if
            end if
        end do
    end do

end subroutine

subroutine wafy(ff,fsta,f,sl,sr,dy,dt,nx,ny)
    implicit none
    integer, intent(in) :: nx, ny
    double precision, intent(in) :: dy, dt
    double precision, dimension(0:nx,0:ny), intent(in) :: fsta, f, sl, sr
    double precision, dimension(0:nx,0:ny), intent(out) :: ff

    integer :: i, j

    do j=0,ny-1
        do i=1,nx-1
            if( sl(i,j)>0.d0 .or. sr(i,j)<0. ) then
                ff(i,j) = fsta(i,j)
            else
                ff(i,j) = 0.5d0*(f(i,j)*(1.d0+sl(i,j)*dt/dy)+fsta(i,j)*(sr(i,j)-sl(i,j))*dt/dy+f(i,j+1)*(1.d0-sr(i,j)*dt/dy) )
            end if
        end do
    end do

end subroutine


subroutine tvdwafy(ff,fsta,f,dq,sl,sr,dy,dt,nx,ny)
    implicit none
    integer, intent(in) :: nx, ny
    double precision, intent(in) :: dy, dt
    double precision, dimension(0:nx,0:ny), intent(in) :: fsta, f, sl, sr, dq
    double precision, dimension(0:nx,0:ny), intent(out) :: ff

    integer :: i, j
    double precision :: cr, cl, ull, urr, phi

    do i=1,nx-1
        j=0
!        if( sl(i,j)>0.d0 .or. sr(i,j)<0. ) then
            ff(i,j) = fsta(i,j)
!        else
!            ff(i,j) = 0.5d0*(f(i,j)*(1.d0+sl(i,j)*dt/dy)+fsta(i,j)*(sr(i,j)-sl(i,j))*dt/dy+f(i,j+1)*(1.d0-sr(i,j)*dt/dy) )
!        end if

        j=ny-1
!        if( sl(i,j)>0.d0 .or. sr(i,j)<0. ) then
            ff(i,j) = fsta(i,j)
!        else
!            ff(i,j) = 0.5d0*(f(i,j)*(1.d0+sl(i,j)*dt/dy)+fsta(i,j)*(sr(i,j)-sl(i,j))*dt/dy+f(i,j+1)*(1.d0-sr(i,j)*dt/dy) )
!        end if
    end do

    do j=1,ny-2
        do i=1,nx-1
            if( sl(i,j)>0.d0 .or. sr(i,j)<0. ) then
                ff(i,j) = fsta(i,j)
            else
                cl = sl(i,j)*dt/dy
                call phical(phi,cl,dq(i,j-1),dq(i,j),dq(i,j+1))
                ull = sign(1.d0,cl)*0.5d0*phi*(fsta(i,j)-f(i,j))
                
                cr = sr(i,j)*dt/dy
                call phical(phi,cr,dq(i,j-1),dq(i,j),dq(i,j+1))
                urr = sign(1.d0,cr)*0.5d0*phi*(f(i,j+1)-fsta(i,j))

                ff(i,j) = 0.5d0*(f(i,j)+f(i,j+1))-ull-urr
            end if
        end do
    end do

end subroutine

subroutine tvdwafyz(ff,fsta,f,dq,sl,sc,sr,v,h,dy,dt,g,nx,ny)
implicit none
    integer, intent(in) :: nx, ny
    double precision, intent(in) :: dy, dt, g
    double precision, dimension(0:nx,0:ny), intent(in) :: fsta, f, sl, sc, sr, dq, v, h
    double precision, dimension(0:nx,0:ny), intent(out) :: ff

    integer :: i, j
    double precision :: cr, cl, ull, urr, phi, ust

!$omp single
    do i=1,nx-1
        j=0
        ust = 0.5d0*(v(i,j)+v(i,j+1))
        if( sl(i,j)>0.d0 .or. sr(i,j)<0. ) then
            ff(i,j) = fsta(i,j)
        else
            if (ust>0.d0) then
                ff(i,j) = (f(i,j)*(1.d0+sl(i,j)*dt/dy)+fsta(i,j)*dt/dy*(sc(i,j)-sl(i,j))+f(i,j+1)*(1.d0-sc(i,j)*dt/dy))*0.5d0
            else
                ff(i,j) = (f(i,j)*(1.d0+sc(i,j)*dt/dy)+fsta(i,j)*dt/dy*(sr(i,j)-sc(i,j))+f(i,j+1)*(1.d0-sr(i,j)*dt/dy))*0.5d0
            end if
        end if

        j=ny-1
        ust = 0.5d0*(v(i,j)+v(i,j+1))
        if( sl(i,j)>0.d0 .or. sr(i,j)<0. ) then
            ff(i,j) = fsta(i,j)
        else
            if (ust>0.d0) then
                ff(i,j) = (f(i,j)*(1.d0+sl(i,j)*dt/dy)+fsta(i,j)*dt/dy*(sc(i,j)-sl(i,j))+f(i,j+1)*(1.d0-sc(i,j)*dt/dy))*0.5d0
            else
                ff(i,j) = (f(i,j)*(1.d0+sc(i,j)*dt/dy)+fsta(i,j)*dt/dy*(sr(i,j)-sc(i,j))+f(i,j+1)*(1.d0-sr(i,j)*dt/dy))*0.5d0
            end if
        end if
    end do
!$omp end single

!$omp do private(i,j,ust,cl,cr,phi,ull,urr)
    do j=1,ny-2
        do i=1,nx-1
            if( sl(i,j)>0.d0 .or. sr(i,j)<0. ) then
                ff(i,j) = fsta(i,j)
            else
                ust = 0.5d0*(v(i,j)+v(i,j+1))
                if( ust>0.d0 ) then
                    cl = sl(i,j)*dt/dy
                    cr = sc(i,j)*dt/dy
                else
                    cl = sc(i,j)*dt/dy
                    cr = sr(i,j)*dt/dy
                end if

                call phical(phi,cl,dq(i,j-1),dq(i,j),dq(i,j+1))
                ull = sign(1.d0,cl)*0.5d0*phi*(fsta(i,j)-f(i,j))
                
                call phical(phi,cr,dq(i,j-1),dq(i,j),dq(i,j+1))
                urr = sign(1.d0,cr)*0.5d0*phi*(f(i,j+1)-fsta(i,j))

                ff(i,j) = 0.5d0*(f(i,j)+f(i,j+1))-ull-urr
            end if
        end do
    end do

end subroutine

subroutine boundary(h,qx,qy,z,dz,dis,wid,snm,ib,dy,nx,ny)
    implicit none
    integer, intent(in)  :: nx, ny
    double precision, intent(in) :: dis, wid, dy, snm, ib
    double precision, dimension(0:nx,0:ny), intent(inout) :: h, qx, qy, z, dz

    integer :: i, j
    
!$omp single
    do j=0,ny
        h( 0,j) = h(   1,j)
        h(nx,j) = h(nx-1,j)

!        qx( 0,j) = 0.d0
!        qx(nx,j) = 0.d0
!        qx( 0,j) = -qx(   1,j)
!        qx(nx,j) = -qx(nx-1,j)

        qx( 0,j) = dis/wid
!        if( j>0.4*ny .and. j<0.6*ny ) then
!            qx(0,j) = dis/(0.2*dble(ny)*dy)
!            h(0,j) = (snm*qx(0,j)/ib**0.5)**0.6d0
!        end if
        qx(nx,j) = qx(nx-1,j)

        qy( 0,j) = qy(   1,j)
        qy(nx,j) = qy(nx-1,j)

        dz( 0,j) = 0.d0
        dz(nx,j) = dz(nx-1,j)

        z( 0,j) = z( 0,j)+dz( 0,j)
        z(nx,j) = z(nx,j)+dz(nx,j)
    end do

    do i=0,nx
        h(i, 0) = h(i,   1)
        h(i,ny) = h(i,ny-1)

        qx(i, 0) = qx(i,   1)
        qx(i,ny) = qx(i,ny-1)

!        qy(i, 0) = 0.d0
!        qy(i,ny) = 0.d0
        qy(i, 0) = -qy(i,   1)
        qy(i,ny) = -qy(i,ny-1)

        dz(i, 0) = dz(i,   1)
        dz(i,ny) = dz(i,ny-1)

!        z(i, 0) = z(i, 0)+dz(i, 0)
!        z(i,ny) = z(i,ny)+dz(i,ny)
        
        z(i, 0) = z(i,   1)
        z(i,ny) = z(i,ny-1)
    end do
!$omp end single

end subroutine

subroutine out2paraview(x,y,z,h,u,v,z0,nx,ny,fpout,tt)
    implicit none
    integer, intent(in) :: nx, ny
    double precision, intent(in) :: tt
    double precision, dimension(0:nx,0:ny), intent(in)  :: x, y, z, h, u, v, z0
    character(20),intent(in) :: fpout

    integer ::  i, j

    open(100,file=fpout,status='unknown')
				
		write(100,'(a26)') '# vtk DataFile Version 3.0'
		write(100,'(a13,f10.3,a7)') 'Time= ', tt
		write(100,'(a5)') 'ASCII'
		write(100,'(a23)') 'DATASET STRUCTURED_GRID'
		write(100,'(a10,3i8)') 'DIMENSIONS', nx+1, ny+1, 1
		write(100,'(a6,i15,a7)') 'POINTS', (nx+1)*(ny+1), 'double'
			
		do j=0,ny
			do i=0,nx
				write(100,*) x(i,j), y(i,j), 0.d0
			end do
		end do
				
		write(100,'(a10,i15)') 'POINT_DATA', (nx+1)*(ny+1)
		write(100,'(a7,2x,a5,2x,a7,i4)') 'SCALARS', 'H', 'double', 1
		write(100,'(a20)') 'LOOKUP_TABLE default'
			
		do j=0,ny
			do i=0,nx
				write(100,*) z(i,j)+h(i,j)
			end do
		end do
        
		write(100,'(a7,2x,a5,2x,a7,i4)') 'SCALARS', 'hs', 'double', 1
		write(100,'(a20)') 'LOOKUP_TABLE default'
			
		do j=0,ny
			do i=0,nx
				write(100,*) h(i,j)
			end do
		end do

        write(100,'(a7,2x,a5,2x,a7,i4)') 'SCALARS', 'Z', 'double', 1
		write(100,'(a20)') 'LOOKUP_TABLE default'
			
		do j=0,ny
			do i=0,nx
				write(100,*) z(i,j)
			end do
		end do
        
        write(100,'(a7,2x,a5,2x,a7,i4)') 'SCALARS', 'DZ', 'double', 1
		write(100,'(a20)') 'LOOKUP_TABLE default'
			
		do j=0,ny
			do i=0,nx
				write(100,*) z(i,j)-z0(i,j)
			end do
		end do

        write(100,'(a23)') 'VECTORS velocity double'
        do j=0,ny
            do i=0,nx
                write(100,*) u(i,j),v(i,j),0.
            end do
        end do
				
	close(100) 

end subroutine

subroutine out2txt(x,y,z,h,u,v,z0,nx,ny,fpout,tt)
    implicit none
    integer, intent(in) :: nx, ny
    double precision, intent(in) :: tt
    double precision, dimension(0:nx,0:ny), intent(in)  :: x, y, z, h, u, v, z0
    character(20),intent(in) :: fpout

    integer ::  i, j

    open(101,file=fpout,status='unknown')
				
		write(101,*) tt, nx, ny

        do j=0,ny
            do i=0,nx
                write(101,*) x(i,j), y(i,j), z(i,j)-z0(i,j)
            end do
        end do
			
	close(101) 

end subroutine

program SWE_HLL

    implicit none 
    integer :: i, j, nx, ny, m, ii, iomp
    integer :: cal_t1,cal_t2,t_rate,t_max,t_diff
    double precision :: dx, dy, dt, rdx, rdy, ct, dt1, dt2, hmin
    double precision :: g, nu, dis, chlen, wid, snm, ib, spec, diam, poro, tuk, etime, bedtime, ib2, z00
    double precision :: h0, hmax, tt, optime, ss, cl, cr, ust, cst
    double precision :: pressure, roughness, advection, volume
    double precision :: tsc, tau, vv, dzdx, dzdy, dzdn, qbs, qbn, uc, vc, hc, frc, mu_s, sigl, sigr, pi, wdz
    double precision :: wsl1, wsc1, wsr1, wsl2, wsc2, wsr2
    double precision :: dql, dqc, dqr, ucc, cc, phi
    double precision, dimension(:,:), allocatable :: x, y, z, wz, z0, dz, h, wh, u, v, qx, wqx, qy, wqy, sr, sl, sc
    double precision, dimension(:,:), allocatable :: f1, f2, f3, ffx, ffy, qbsta
    double precision, dimension(:,:), allocatable :: qbx, qby, qbx1, qby1, qbx2, qby2, f4
    double precision, dimension(:,:), allocatable :: qsta, fsta, dh
    character(20) :: fpout, fptxt


    g  = 9.81d0
    nu = 1e-6
    pi = 3.14159d0

    dis     = 0.0025d0
    chlen   = 50.d0
    wid     = 0.45d0
    snm     = 0.0145d0
    ib      = 0.01d0
    spec    = 1.65d0
    diam    = 0.000755d0
    poro    = 0.4d0
    tsc     = 0.035d0
    mu_s    = 0.7d0
    hmin    = 0.0001d0

    ib2 = ib    !*5.

    nx = 2000
    ny = 45
    dx = chlen/dble(nx)
    dy = wid/dble(ny)
    dt = 0.02d0
    ct = 0.5d0
    
    rdx = 1.d0/dx
    rdy = 1.d0/dy
    
    tuk   = 60.
    bedtime = 60.
    etime = tuk*1000.

    iomp = 16

    allocate( x(0:nx,0:ny), y(0:nx,0:ny), z(0:nx,0:ny), wz(0:nx,0:ny), z0(0:nx,0:ny), dz(0:nx,0:ny), h(0:nx,0:ny), wh(0:nx,0:ny), qbsta(0:nx,0:ny) )
    allocate( u(0:nx,0:ny), v(0:nx,0:ny), qx(0:nx,0:ny), wqx(0:nx,0:ny), qy(0:nx,0:ny), wqy(0:nx,0:ny) )
    allocate( sr(0:nx,0:ny), sl(0:nx,0:ny), sc(0:nx,0:ny), ffx(0:nx,0:ny), ffy(0:nx,0:ny), f1(0:nx,0:ny), f2(0:nx,0:ny), f3(0:nx,0:ny) )
    allocate( qbx(0:nx,0:ny), qby(0:nx,0:ny), qbx1(0:nx,0:ny), qby1(0:nx,0:ny), qbx2(0:nx,0:ny), qby2(0:nx,0:ny), f4(0:nx,0:ny) )
    allocate( qsta(0:nx,0:ny), fsta(0:nx,0:ny), dh(0:nx,0:ny) )

    do j=0,ny
        x(0,j) = 0.d0
        do i=1,nx
            x(i,j) = x(i-1,j)+dx
        end do
    end do
    
    do i=0,nx
        y(i,0) = 0.d0
        do j=1,ny
            y(i,j) = y(i,j-1)+dy
        end do
    end do
       
!    z00  = chlen*ib
    z00 = chlen*ib*0.75 + chlen*ib2*0.25

    do j=0,ny
        z(0,j) = z00
        do i=1,nx
!!            z(i,j) = z(i-1,j)-ib*dx
            if ( i<0.3*nx .or. (0.4*nx<i .and. i<0.6*nx) .or. i>0.7*nx ) then
                z(i,j) = z(i-1,j)-ib*dx
                h(i,j) = (snm*dis/(wid*ib**0.5))**0.6
                qx(i,j) = dis/wid
                qy(i,j) = 0.d0
            else
                z(i,j) = z(i-1,j)-ib2*dx
                h(i,j) = (snm*dis/(wid*ib2**0.5))**0.6
                qx(i,j) = dis/wid
                qy(i,j) = 0.d0
            end if
        end do
    end do
    
    write(*,*) "Fr= ", h(0.5*nx,0.5*ny)**(2./3.)*ib**0.5/snm/dsqrt(g*h(0.5*nx,0.5*ny))

    do j=0,ny
        do i=0,nx
            z0(i,j) = z(i,j)
        end do
    end do

    do j=0,0.1*ny
        do i=0.1*nx,0.11*nx
            z(i,j) = z(i,j)+h(i,j)*0.2
        end do
    end do

!    do j=0,ny
!        z(0,j) = z00
!        do i=1,nx
!            if ( i<0.3*nx .or. (0.4*nx<i .and. i<0.6*nx) .or. i>0.7*nx ) then
!                z(i,j) = z(i-1,j)-ib*dx
!                h(i,j) = hmin
!                qx(i,j) = 0.d0
!                qy(i,j) = 0.d0
!            else
!                z(i,j) = z(i-1,j)-ib2*dx
!                h(i,j) = hmin
!                qx(i,j) = 0.d0
!                qy(i,j) = 0.d0
!            end if
!        end do
!    end do
    
!    do j=0.4*ny,0.6*ny
!        do i=0.48*nx,0.52*nx
!            z(i,j) = z(i,j)+h(i,j)*0.01
!        end do
!    end do

    dz = 0.d0

!    do j=0,ny
!        do i=0.1*nx,0.9*nx
!            z(i,j) = z(i,j)+0.005*cos(2.d0*pi*(x(i,j)-x(0.1*nx,j))/7.d0)*sin(2.d0*pi*(y(i,j)-wid*0.5)/(2.d0*wid))
!        end do
!    end do


!    do j=0,ny
!        do i=0,nx
!            h(i,j) = 2.d0
!            if(i>0.4*nx .and. i<0.6*nx) then
!                if(j>0.4*ny .and. j<0.6*ny) then
!                    h(i,j) = 10.d0
!                end if
!            end if
!!            if( x(i,j)<500.) then
!!                h(i,j) = 10.
!!            end if
!            qx(i,j) = 0.d0
!            qy(i,j) = 0.d0
!            
!            u(i,j) = qx(i,j)/h(i,j)
!            v(i,j) = qy(i,j)/h(i,j)
!        end do
!    end do

    call boundary(h,qx,qy,z,dz,dis,wid,snm,ib,dy,nx,ny)

    do j=0,ny
        do i=0,nx
            if( h(i,j)>hmin ) then
                u(i,j) = qx(i,j)/h(i,j)
                v(i,j) = qy(i,j)/h(i,j)
            else
                u(i,j) = 0.d0
                v(i,j) = 0.d0
            end if
        end do
    end do

    wh  = h
    wz  = z
    wqx = qx
    wqy = qy 
    
    tt = 0.
    optime = 0.
    
	m = 0
    fpout = 'vc000.vtk'
    fptxt = 'tt000.txt'

    call system_clock(cal_t1)	!h160105 計算開始時時刻

!$	call omp_set_num_threads(iomp)

!$omp parallel

    do
        
!$omp do private(i,j)
        do j=0,ny
            do i=0,nx
                h(i,j) = wh(i,j)
                z(i,j) = wz(i,j)
                qx(i,j) = wqx(i,j)
                qy(i,j) = wqy(i,j)
            end do
        end do
        
!$omp do private(i,j,vv,tau,qbs)
        do j=0,ny
            do i=0,nx
                if ( wh(i,j)>hmin ) then
                    u(i,j) = wqx(i,j)/wh(i,j)
                    v(i,j) = wqy(i,j)/wh(i,j)
                    
                    ffx(i,j) = wqx(i,j)**2.d0/wh(i,j)+0.5d0*g*wh(i,j)**2.d0     !+g*wh(i,j)*z(i,j)
                    ffy(i,j) = wqx(i,j)*wqy(i,j)/wh(i,j)

                    vv = dsqrt(u(i,j)**2.d0+v(i,j)**2.d0)
                    tau = snm**2.d0*vv**2.d0/(spec*diam*wh(i,j)**(1.d0/3.d0))
    !                dzdx = (-z(i-1,j)+z(i+1,j))*0.5d0*rdx
    !                dzdy = (-z(i,j-1)+z(i,j+1))*0.5d0*rdy
    !                dzdn = (-v(i,j)*dzdx+u(i,j)*dzdy)/vv

                    if( tau>tsc ) then
                        qbs = 4.d0*(tau-tsc)**1.5d0*dsqrt(spec*g*diam**3.d0)
    !                    qbn = -qbs*dsqrt(tsc/tau)/mu_s*dzdn
                        qbx1(i,j) = u(i,j)/vv*qbs
    !                    qbx2(i,j) = -v(i,j)/vv*qbn
                    else
                        qbx1(i,j) = 0.d0
    !                    qbx2(i,j) = 0.d0
                    end if
                else
                    u(i,j) = 0.d0
                    v(i,j) = 0.d0
                    ffx(i,j) = 0.d0
                    ffy(i,j) = 0.d0
                    qbx1(i,j) = 0.d0
                end if

            end do
        end do
        
!$omp do private(i,j,uc,vc,hc,vv,tau,dzdx,dzdy,dzdn,qbs,qbn)
        do j=1,ny-1
            do i=0,nx-1
                uc = (u(i,j)+u(i+1,j))*0.5d0
                vc = (v(i,j)+v(i+1,j))*0.5d0
                hc = (wh(i,j)+wh(i+1,j))*0.5d0

                vv = dsqrt(uc**2.d0+vc**2.d0)

                if ( hc>hmin .and. vv>1e-5 ) then
                    tau = snm**2.d0*vv**2.d0/(spec*diam*hc**(1.d0/3.d0))
                else
                    tau = 0.d0
                end if

                dzdx = (-z(i,j)+z(i+1,j))*rdx
                dzdy = (-(z(i,j-1)+z(i+1,j-1))+(z(i,j+1)+z(i+1,j+1)))*0.25d0*rdy
                dzdn = (-vc*dzdx+uc*dzdy)/vv

                if( tau>tsc ) then
                    qbs = 4.d0*(tau-tsc)**1.5d0*dsqrt(spec*g*diam**3.d0)
                    qbn = -qbs*dsqrt(tsc/tau)/mu_s*dzdn
!                    qbx2(i,j) = uc/vv*qbs-vc/vv*qbn
                    qbx2(i,j) = -vc/vv*qbn
                else
                    qbx2(i,j) = 0.d0
                end if

            end do
        end do

!$omp do private(i,j,wsl1,wsc1,wsr1,wsl2,wsc2,wsr2)
        do j=1,ny-1
            do i=0,nx-1
!                uc = (u(i,j)+u(i+1,j))*0.5d0
!                vc = (v(i,j)+v(i+1,j))*0.5d0
!                hc = (wh(i,j)+wh(i+1,j))*0.5d0
!                call propagation(sl(i,j),sc(i,j),sr(i,j),uc,vc,hc,g,snm,spec,diam,tsc,poro,hmin)
                call propagation(wsl1, wsc1, wsr1,u(i,j),v(i,j),wh(i,j),g,snm,spec,diam,tsc,poro,hmin)
                call propagation(wsl2, wsc2, wsr2,u(i+1,j),v(i+1,j),wh(i+1,j),g,snm,spec,diam,tsc,poro,hmin)
                sl(i,j) = wsl1  !(wsl1+wsl2)*0.5d0
                sc(i,j) = (wsc1+wsc2)*0.5d0
                sr(i,j) = wsr2  !(wsr1+wsr2)*0.5d0
            end do
        end do

!$omp single
        dt1 = 999.d0
        do j=1,ny-1
            do i=0,nx-1
                dt1 = min(dabs(dx/sl(i,j)), dabs(dx/sc(i,j)), dabs(dx/sr(i,j)),dt1)
            end do
        end do
!$omp end single
               
!        call hllx1storder(f1,wqx,sl,sr,wh+z,nx,ny)
!        call hllx1storder(f2,ffx,sl,sr,wqx,nx,ny)
        
!        call hllx1storder(qsta,wqx,sl,sr,wh+z,nx,ny)
        call hllx1storder(qsta,wqx,sl,sr,wh,nx,ny)
        call hllx1storder(fsta,ffx,sl,sr,wqx,nx,ny)

!!        call wafx(f1,qsta,wqx,sl,sr,dx,dt,nx,ny)
!!        call wafx(f2,fsta,ffx,sl,sr,dx,dt,nx,ny)

!$omp do private(i,j)
        do j=0,ny
            do i=0,nx-1
!                dh(i,j) = (wh(i+1,j)+z(i+1,j))-(wh(i,j)+z(i,j))
                dh(i,j) = (wh(i+1,j))-(wh(i,j))
            end do
        end do

        call tvdwafx(f1,qsta,wqx,dh,sl,sr,dx,dt,nx,ny)
        call tvdwafx(f2,fsta,ffx,dh,sl,sr,dx,dt,nx,ny)

!$omp do private(i,j,hc,ust,dqc,dql,dqr,cc,phi,ucc)
        do j=1,ny-1
            do i=0,nx-1
!                ust = 0.5d0*(u(i,j)+u(i+1,j))+dsqrt(g*wh(i,j))-dsqrt(g*wh(i+1,j))
!!                ust = 0.5d0*(u(i,j)+u(i+1,j))
!!                if(ust>0) then
!!                    f3(i,j) = f1(i,j)*v(i,j)
!!!                   f3(i,j) = ffy(i,j)
!!                else
!!                    f3(i,j) = f1(i,j)*v(i+1,j)
!!!                   f3(i,j) = ffy(i+1,j)
!!                end if
!                f3(i,j) = f1(i,j)*(0.5d0*(v(i,j)+v(i+1,j)+ust*dt/dx*(v(i,j)-v(i+1,j))))


                hc = (wh(i,j)+wh(i+1,j))*0.5d0
                if( hc>hmin ) then
                    ust = f1(i,j)/hc
!                    if( ust>0 ) then
!                        f3(i,j) = ust*wqy(i,j)
!                    else
!                        f3(i,j) = ust*wqy(i+1,j)
!                    end if

                    dqc = wh(i+1,j)-wh(i,j)
                    
                    if( i==0 ) then
                        dql = 0.d0
                    else
                        dql = wh(i,j)-wh(i-1,j)
                    end if

                    if( i==nx-1 ) then
                        dqr = 0.d0
                    else
                        dqr = wh(i+2,j)-wh(i+1,j)
                    end if

                    cc = ust*dt/dy
                    call phical(phi,cc,dql,dqc,dqr)
                    ucc = sign(1.d0,cc)*0.5d0*phi*(wqy(i+1,j)-wqy(i,j))

                    f3(i,j) = ust*( 0.5d0*(wqy(i,j)+wqy(i+1,j)) - ucc )

!                    f3(i,j) = ust*0.5d0*(wqy(i,j)+wqy(i+1,j)-ust*dt/dy*(wqy(i+1,j)-wqy(i,j)))
                else
                    f3(i,j) = 0.0d0
                end if

            end do
        end do
        
!!        call hllx1storder(f4,qbx1,sl,sc,wz,nx,ny)
!        call hllzx1storder(f4,qbx1,sl,sc,sr,wz,u,snm,dx,wqx,wh,nx,ny)
        call hllzx1storder(qbsta,qbx1,sl,sc,sr,wz,u,snm,dx,wqx,wh,nx,ny)

!        do j=1,ny-1
!            do i=0,nx-1
!                f4(i,j) = (qbx1(i,j)*(1.d0+sl(i,j)*dt/dx)+qbsta(i,j)*dt/dx*(sc(i,j)-sl(i,j))+qbx1(i+1,j)*(1.d0-sc(i,j)*dt/dx))*0.5d0
!            end do
!        end do

!$omp do private(i,j)
        do j=0,ny
            do i=0,nx-1
                dh(i,j) = z(i+1,j)-z(i,j)
            end do
        end do

        call tvdwafx(f4,qbsta,qbx1,dh,sl,sc,dx,dt,nx,ny)

!$omp do private(i,j)
        do j=1,ny-1
            do i=1,nx-1
                h(i,j) = wh(i,j)-(-f1(i-1,j)+f1(i,j))*rdx*dt
                if ( h(i,j)<=hmin ) h(i,j) = hmin
            end do
        end do
    
!$omp do private(i,j,pressure,roughness,sigl,sigr,advection)
        do j=1,ny-1
            do i=1,nx-1
                pressure = -g*wh(i,j)*(-z(i-1,j)+z(i+1,j))*rdx*0.5d0
                roughness = g*snm**2.d0*dsqrt(u(i,j)**2.d0+v(i,j)**2.d0)/wh(i,j)**(4.d0/3.d0)

                sigl = f2(i-1,j)-(sr(i-1,j))/(sr(i-1,j)-sl(i-1,j))*g*(wh(i,j)+wh(i-1,j))*0.5*(wz(i,j)-wz(i-1,j))
                sigr = f2(i,j)-(sl(i,j))/(sr(i,j)-sl(i,j))*g*(wh(i+1,j)+wh(i,j))*0.5*(wz(i+1,j)-wz(i,j))

                advection = (-sigl+sigr)*rdx
!                advection = (-f2(i-1,j)+f2(i,j))*rdx

!                qx(i,j) = (wqx(i,j)+(-advection+pressure)*dt)/(1.0+roughness*dt)
                qx(i,j) = (wqx(i,j)+(-advection)*dt)/(1.0+roughness*dt)
            end do
        end do

!$omp do private(i,j)
        do j=1,ny-1
            do i=1,nx-1
                qy(i,j) = wqy(i,j)-(-f3(i-1,j)+f3(i,j))*rdx*dt
            end do
        end do
        
        if( tt>bedtime ) then

!$omp do private(i,j,wdz)
            do j=1,ny-1
                do i=1,nx-1
!                    wdz = -((-f4(i-1,j)+f4(i,j))*rdx)*dt/(1.d0-poro)
                    wdz = -((-f4(i-1,j)+f4(i,j))*rdx+(-qbx2(i-1,j)+qbx2(i,j))*rdx)*dt/(1.d0-poro)
!                    wdz = -((-qbx2(i-1,j)+qbx2(i,j))*rdx)*dt/(1.d0-poro)
!                    z(i,j) = wz(i,j)-((-f4(i-1,j)+f4(i,j))*rdx+(-qbx2(i-1,j)+qbx2(i+1,j))*0.5*rdx)*dt/(1.d0-poro)
                    z(i,j) = wz(i,j)+wdz
                    dz(i,j) = wdz
!                    z(i,j) = wz(i,j)-((-qbx2(i-1,j)+qbx2(i,j))*rdx)*dt/(1.d0-poro)
!                    z(i,j) = wz(i,j)-((-f4(i-1,j)+f4(i,j))*rdx)*dt/(1.d0-poro)
                end do
            end do

        end if

        call boundary(h,qx,qy,z,dz,dis,wid,snm,ib,dy,nx,ny)
        
!$omp do private(i,j,vv,tau,qbs)
        do j=0,ny
            do i=0,nx
                if ( h(i,j)>hmin ) then
                    u(i,j) = qx(i,j)/h(i,j)
                    v(i,j) = qy(i,j)/h(i,j)
                    
                    ffx(i,j) = qx(i,j)*qy(i,j)/h(i,j)
                    ffy(i,j) = qy(i,j)**2.d0/h(i,j)+0.5d0*g*h(i,j)**2.d0    !+g*h(i,j)*z(i,j)
    
                    vv = dsqrt(u(i,j)**2.d0+v(i,j)**2.d0)
                    tau = snm**2.d0*vv**2.d0/(spec*diam*h(i,j)**(1.d0/3.d0))
!                    dzdx = (-z(i-1,j)+z(i+1,j))*0.5d0*rdx
!                    dzdy = (-z(i,j-1)+z(i,j+1))*0.5d0*rdy
!                    dzdn = (-v(i,j)*dzdx+u(i,j)*dzdy)/vv

                    if( tau>tsc ) then
                        qbs = 4.d0*(tau-tsc)**1.5d0*dsqrt(spec*g*diam**3.d0)
!                        qbn = -qbs*dsqrt(tsc/tau)/mu_s*dzdn
                        qby1(i,j) = v(i,j)/vv*qbs
    !                    qby2(i,j) = u(i,j)/vv*qbn
                    else
                        qby1(i,j) = 0.d0
    !                    qby2(i,j) = 0.d0
                    end if
                else
                    u(i,j) = 0.d0
                    v(i,j) = 0.d0
                    ffx(i,j) = 0.d0
                    ffy(i,j) = 0.d0
                    qby1(i,j) = 0.d0
                end if
            end do
        end do
             
!$omp do private(i,j,uc,vc,hc,vv,tau,dzdx,dzdy,dzdn,qbs,qbn)
        do j=0,ny-1
            do i=1,nx-1
                uc = (u(i,j)+u(i,j+1))*0.5d0
                vc = (v(i,j)+v(i,j+1))*0.5d0
                hc = (h(i,j)+h(i,j+1))*0.5d0
                
                vv = dsqrt(uc**2.d0+vc**2.d0)

                if ( hc>hmin .and. vv>1e-5 ) then
                    tau = snm**2.d0*vv**2.d0/(spec*diam*hc**(1.d0/3.d0))
                else
                    tau = 0.d0
                end if

                dzdx = (-(z(i-1,j)+z(i-1,j+1))+(z(i+1,j)+z(i+1,j+1)))*0.25d0*rdx
                dzdy = (-z(i,j)+z(i,j+1))*rdy
                dzdn = (-vc*dzdx+uc*dzdy)/vv

                if( tau>tsc ) then
                    qbs = 4.d0*(tau-tsc)**1.5d0*dsqrt(spec*g*diam**3.d0)
                    qbn = -qbs*dsqrt(tsc/tau)/mu_s*dzdn
!                    qby2(i,j) = vc/vv*qbs+uc/vv*qbn
                    qby2(i,j) = uc/vv*qbn
                else
                    qby2(i,j) = 0.d0
                end if
            end do
        end do

!$omp do private(i)
        do i=0,nx
            qby1(i, 0) = -qby1(i,   1)
            qby1(i,ny) = -qby1(i,ny-1)
!            qby1(i, 0) = 0.d0
!            qby1(i,ny) = 0.d0
        end do

!$omp do private(i,j,wsl1,wsc1,wsr1,wsl2,wsc2,wsr2)
        do j=0,ny-1
            do i=1,nx-1
!                uc = (u(i,j)+u(i,j+1))*0.5d0
!                vc = (v(i,j)+v(i,j+1))*0.5d0
!                hc = (h(i,j)+h(i,j+1))*0.5d0
!                call propagation(sl(i,j),sc(i,j),sr(i,j),vc,uc,hc,g,snm,spec,diam,tsc,poro,hmin)
                call propagation(wsl1, wsc1, wsr1,v(i,j),u(i,j),h(i,j),g,snm,spec,diam,tsc,poro,hmin)
                call propagation(wsl2, wsc2, wsr2,v(i,j+1),u(i,j+1),h(i,j+1),g,snm,spec,diam,tsc,poro,hmin)
                sl(i,j) = wsl1  !(wsl1+wsl2)*0.5d0
                sc(i,j) = (wsc1+wsc2)*0.5d0
                sr(i,j) = wsr2  !(wsr1+wsr2)*0.5d0
            end do
        end do

!$omp single
        dt2 = 999.d0
        do j=0,ny-1
            do i=1,nx-1
                dt2 = min(dabs(dy/sl(i,j)), dabs(dy/sc(i,j)), dabs(dy/sr(i,j)),dt2)
            end do
        end do
!$omp end single

!        call hlly1storder(f1,qy,sl,sr,h+z,nx,ny)
!        call hlly1storder(f3,ffy,sl,sr,qy,nx,ny)
        
!        call hlly1storder(qsta,qy,sl,sr,h+z,nx,ny)
        call hlly1storder(qsta,qy,sl,sr,h,nx,ny)
        call hlly1storder(fsta,ffy,sl,sr,qy,nx,ny)
        
!!        call wafy(f1,qsta,qy,sl,sr,dy,dt,nx,ny)
!!        call wafy(f3,fsta,ffy,sl,sr,dy,dt,nx,ny)

!$omp do private(i,j)
        do j=0,ny-1
            do i=0,nx
!                dh(i,j) = (h(i,j+1)+z(i,j+1))-(h(i,j)+z(i,j))
                dh(i,j) = (h(i,j+1))-(h(i,j))
            end do
        end do

        call tvdwafy(f1,qsta,qy,dh,sl,sr,dy,dt,nx,ny)
        call tvdwafy(f3,fsta,ffy,dh,sl,sr,dy,dt,nx,ny)

!$omp do private(i,j,hc,ust,dqc,dql,dqr,cc,phi,ucc)
        do j=0,ny-1
            do i=1,nx-1
!                ust = 0.5d0*(v(i,j)+v(i,j+1))+dsqrt(g*h(i,j))-dsqrt(g*h(i,j+1))
!!                ust = 0.5d0*(v(i,j)+v(i,j+1))
!!                if(ust>0) then
!!                    f2(i,j) = f1(i,j)*u(i,j)
!!!                    f2(i,j) = ffx(i,j)
!!                else
!!                    f2(i,j) = f1(i,j)*u(i,j+1)
!!!                    f2(i,j) = ffx(i,j+1)
!!                end if
!                f2(i,j) = f1(i,j)*(0.5d0*(u(i,j)+u(i,j+1)+ust*dt/dy*(u(i,j)-u(i,j+1))))

                                
                hc = (h(i,j)+h(i,j+1))*0.5d0
                if( hc>hmin ) then
                    ust = f1(i,j)/hc
!                    if( ust>0 ) then
!                        f2(i,j) = ust*qx(i,j)
!                    else
!                        f2(i,j) = ust*qx(i,j+1)
!                    end if

                    dqc = h(i,j+1)-h(i,j)
                    
                    if( j==0 ) then
                        dql = 0.d0
                    else
                        dql = h(i,j)-h(i,j-1)
                    end if

                    if( j==ny-1 ) then
                        dqr = 0.d0
                    else
                        dqr = h(i,j+2)-h(i,j+1)
                    end if

                    cc = ust*dt/dx
                    call phical(phi,cc,dql,dqc,dqr)
                    ucc = sign(1.d0,cc)*0.5d0*phi*(qx(i,j+1)-qx(i,j))

                    f2(i,j) = ust*( 0.5d0*(qx(i,j)+qx(i,j+1)) - ucc )

!                    f2(i,j) = ust*0.5d0*(qx(i,j)+qx(i,j+1)-ust*dt/dx*(qx(i,j+1)-qx(i,j)))
                else
                    f2(i,j) = 0.0d0
                end if
            end do
        end do

!!        call hlly1storder(f4,qby1,sl,sc,z,nx,ny)
!        call hllzy1storder(f4,qby1,sl,sc,sr,z,v,nx,ny)
        call hllzy1storder(qbsta,qby1,sl,sc,sr,z,v,nx,ny)

!        do j=0,ny-1
!            do i=1,nx-1
!                ust = 0.5d0*(v(i,j)+v(i,j+1))+dsqrt(g*h(i,j))-dsqrt(g*h(i,j+1))
!                if (ust>0.d0) then
!                    f4(i,j) = (qby1(i,j)*(1.d0+sl(i,j)*dt/dy)+qbsta(i,j)*dt/dy*(sc(i,j)-sl(i,j))+qby1(i,j+1)*(1.d0-sc(i,j)*dt/dy))*0.5d0
!                else
!                    f4(i,j) = (qby1(i,j)*(1.d0+sc(i,j)*dt/dy)+qbsta(i,j)*dt/dy*(sr(i,j)-sc(i,j))+qby1(i,j+1)*(1.d0-sr(i,j)*dt/dy))*0.5d0
!                end if
!            end do
!        end do

!$omp do private(i,j)
        do j=0,ny-1
            do i=0,nx
                dh(i,j) = z(i,j+1)-z(i,j)
            end do
        end do

        call tvdwafyz(f4,qbsta,qby1,dh,sl,sc,sr,v,h,dy,dt,g,nx,ny)
        
!$omp do private(i,j)
        do j=1,ny-1
            do i=1,nx-1
                wh(i,j) = h(i,j)-(-f1(i,j-1)+f1(i,j))*rdy*dt
                if ( wh(i,j)<=hmin ) wh(i,j) = hmin
            end do
        end do
           
!$omp do private(i,j,pressure,roughness,sigl,sigr,advection)
        do j=1,ny-1
            do i=1,nx-1
                pressure = -g*h(i,j)*(-z(i,j-1)+z(i,j+1))*rdy*0.5d0
                roughness = g*snm**2.d0*dsqrt(u(i,j)**2.d0+v(i,j)**2.d0)/h(i,j)**(4.d0/3.d0)

                sigl = f3(i,j-1)-(sr(i,j-1))/(sr(i,j-1)-sl(i,j-1))*g*(h(i,j)+h(i,j-1))*0.5*(z(i,j)-z(i,j-1))
                sigr = f3(i,j)-(sl(i,j))/(sr(i,j)-sl(i,j))*g*(h(i,j+1)+h(i,j))*0.5*(z(i,j+1)-z(i,j))

                advection = (-sigl+sigr)*rdy
!                advection = (-f3(i,j-1)+f3(i,j))*rdy
!                wqy(i,j) = (qy(i,j)+(-advection+pressure)*dt)/(1.0+roughness*dt)
                wqy(i,j) = (qy(i,j)+(-advection)*dt)/(1.0+roughness*dt)
            end do
        end do

!$omp do private(i,j)
        do j=1,ny-1
            do i=1,nx-1
                wqx(i,j) = qx(i,j)-(-f2(i,j-1)+f2(i,j))*rdy*dt
            end do
        end do
                
        if( tt>bedtime ) then

!$omp do private(i,j,wdz)
            do j=1,ny-1
                do i=1,nx-1
!                    wdz = -((-f4(i,j-1)+f4(i,j))*rdy)*dt/(1.d0-poro)
                    wdz = -((-f4(i,j-1)+f4(i,j))*rdy+(-qby2(i,j-1)+qby2(i,j))*rdy)*dt/(1.d0-poro)
!                    wdz = -((-qby2(i,j-1)+qby2(i,j))*rdy)*dt/(1.d0-poro)
!                    wz(i,j) = z(i,j)-((-f4(i,j-1)+f4(i,j))*rdy+(-qby2(i,j-1)+qby2(i,j+1))*0.5*rdy)*dt/(1.d0-poro)
                    wz(i,j) = z(i,j)+wdz
                    dz(i,j) = wdz
!                    wz(i,j) = z(i,j)-((-qby2(i,j-1)+qby2(i,j))*rdy)*dt/(1.d0-poro)
!                    wz(i,j) = z(i,j)-((-f4(i,j-1)+f4(i,j))*rdy)*dt/(1.d0-poro)
                end do
            end do

        end if

        call boundary(wh,wqx,wqy,wz,dz,dis,wid,snm,ib,dy,nx,ny)

!$omp single
        dt = min(dt1,dt2)*ct
            
        if (optime>tuk .or. m==0) then

            volume = 0.

            do j=1,ny-1
                do i=1,nx-1
                    volume = volume+h(i,j)
                end do
            end do

            write(*,'(f12.5,f10.6,e17.5)') tt, dt, volume

            write(fpout(3:5),'(i3)') m
            
            do ii=3,4
                if(fpout(ii:ii) == ' ') fpout(ii:ii) = '0'
            end do

            write(fptxt(3:5),'(i3)') m
            
            do ii=3,4
                if(fptxt(ii:ii) == ' ') fptxt(ii:ii) = '0'
            end do

            call out2paraview(x,y,z,h,u,v,z0,nx,ny,fpout,tt)
            call out2txt(x,y,z,h,u,v,z0,nx,ny,fptxt,tt)

            if( m>0 ) optime = optime-tuk
            m = m+1
        
        end if
        
        optime = optime+dt
        tt = tt+dt
!$omp end single

        if( tt>etime ) exit
        
    end do

!$omp end parallel

    call system_clock(cal_t2, t_rate, t_max)	!h160105 計算終了時時刻
	if ( cal_t2 < cal_t1 ) then
		t_diff = (t_max - cal_t1) + cal_t2 + 1
	else
		t_diff = cal_t2 - cal_t1
	endif
	write(*,*) "Calcuration time",real(t_diff/t_rate),"sec."
	write(*,*) "Calcuration time",real(t_diff/t_rate)/60.,"min."
	write(*,*) "Calcuration time",real(t_diff/t_rate)/3600.,"hour."

end program SWE_HLL